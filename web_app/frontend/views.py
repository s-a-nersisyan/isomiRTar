from . import frontend
from web_app import app
from web_app.api.db_functions import *


from flask import render_template, abort

import json
from scipy.stats import linregress, spearmanr
from statsmodels.stats.multitest import multipletests


@frontend.route("/", methods=["GET"])
def show_index():
    isomiRs = get_all_isomiRs()
    miRNAs = sorted(list({isomiR.split("|")[0] for isomiR in isomiRs}))
    genes = get_all_genes()
    cancers = get_all_cancers()
    return render_template(
        "index.html",
        miRNAs=miRNAs,
        isomiRs=isomiRs,
        genes=genes,
        cancers=cancers
    )


@frontend.route("/molecule/<molecule>", methods=["GET"])
def show_molecule(molecule):
    if not check_molecule(molecule):
        abort(400)
    
    # First, get molecule expression for boxplot
    df = get_molecule_expression_pan_cancer(molecule)
    order = df.groupby("cancer").quantile(0.75).sort_values("expression", ascending=False).index
    df = pd.concat([df.loc[df["cancer"] == cancer] for cancer in order])
    df["expression"] = 2**df["expression"] - 1
    cancers = df["cancer"].unique()
    at_least_one_highly_expressed = "*" in [cancer[-1] for cancer in cancers]
    
    expression = df.rename(columns={"cancer": "x", "expression": "y"}).to_dict(orient="list")
    expression["type"] = "box"
    expression = [expression]

    # Then, retrieve the results of isomiR target prediction
    if molecule.startswith("hsa-"):
        is_isomiR = True
        index_col = "target"
        targets_pan_cancer = get_molecule_targeting_pan_cancer(isomiR=molecule)
    else:
        is_isomiR = False
        index_col = "isomir"
        targets_pan_cancer = get_molecule_targeting_pan_cancer(target=molecule)
    
    # Sort miRDB/TargetScan predictions according to the scores
    targets_seq = targets_pan_cancer.set_index(index_col)[["mirdb_score", "targetscan_score"]]
    targets_seq = targets_seq.loc[~targets_seq.index.duplicated()]
    targets_seq["mirdb_score"] *= -1
    df1 = targets_seq.loc[(targets_seq["mirdb_score"].notna()) & (targets_seq["targetscan_score"].notna())]
    df2 = targets_seq.loc[(targets_seq["mirdb_score"].notna()) & (targets_seq["targetscan_score"].isna())]
    df3 = targets_seq.loc[(targets_seq["mirdb_score"].isna()) & (targets_seq["targetscan_score"].notna())]
    targets_seq = pd.concat([
        df1.sort_values(["mirdb_score", "targetscan_score"]),
        df2.sort_values("mirdb_score"),
        df3.sort_values("targetscan_score")
    ])
    targets_seq["mirdb_score"] *= -1
    targets_seq["mirdb_score"] = targets_seq["mirdb_score"].astype("Int64")
    targets_seq = targets_seq.astype("string").fillna("-")
        
    targets_conserved = []
    targets_pan_cancer = targets_pan_cancer.loc[
        (targets_pan_cancer["spearman_corr"] < -0.3) & (targets_pan_cancer["spearman_p_value"] < 0.05)
    ]
    for target, df in targets_pan_cancer.groupby(index_col):
        anti_corr_cancers = df["cancer"].to_list()
        if len(anti_corr_cancers) >= 1:
            if is_isomiR:
                targets_conserved.append((molecule, target, sorted(anti_corr_cancers)))
            else:
                targets_conserved.append((target, molecule, sorted(anti_corr_cancers)))

    targets_conserved = sorted(targets_conserved, key=lambda x: len(x[2]), reverse=True)

    return render_template(
        "molecule/main.html",
        page="molecule",
        molecule=molecule,
        expression=expression,
        cancers=cancers,
        at_least_one_highly_expressed=at_least_one_highly_expressed,
        is_isomiR=is_isomiR,
        targets_seq=targets_seq,
        targets_conserved=targets_conserved
    )


@frontend.route("/miRNA/<miRNA>", methods=["GET"])
def show_miRNA(miRNA):
    if not check_miRNA(miRNA):
        abort(400)

    dfs = get_miRNA_expression_pan_cancer(miRNA)
    expression = [] 
    for molecule, df in dfs:
        expression.append({
            "df": df,
            "name": molecule,
            "type": "box",
            "orientation": "h",
            "max_q75": df.groupby("cancer").quantile(0.75)["expression"].max()
        })
    
    # Get the most abundant isomiR
    i_max = sorted(list(range(len(expression))), key=lambda i: expression[i]["max_q75"])[-1]
    # Sort cancers
    order = expression[i_max]["df"].groupby("cancer").quantile(0.75).sort_values("expression", ascending=True).index
    for i in range(len(expression)):
        df = expression[i]["df"]
        df = pd.concat([
            df.loc[df["cancer"] == cancer] if cancer in list(df["cancer"]) else pd.DataFrame({"cancer": [cancer], "expression": [0.0]})
            for cancer in order
        ])
        del expression[i]["df"]
        expression[i] = {
            **expression[i],
            **df.rename(columns={"cancer": "y", "expression": "x"}).to_dict(orient="list"),
            "num_cancers": len(order)
        }

    # Get predicted targets for each isomiR and make pairwise intersections matrix
    isomiRs = [isomiR for isomiR, _ in dfs]
    targets = {}
    for isomiR in isomiRs:
        targets[isomiR] = set(get_isomiR_targets_unique(isomiR))
    
    intersection, union, jaccard = [], [], []
    for i, isomiR1 in enumerate(isomiRs):
        intersection.append([])
        union.append([])
        jaccard.append([])
        for j, isomiR2 in enumerate(isomiRs):
            intersection_ = len(targets[isomiR1] & targets[isomiR2])
            union_ = len(targets[isomiR1] | targets[isomiR2])
            jaccard_ = intersection_ / union_ if union_ != 0 else 0
            intersection[i].append(intersection_)
            union[i].append(union_)
            jaccard[i].append(jaccard_)

    return render_template(
        "miRNA/main.html",
        page="miRNA",
        miRNA=miRNA,
        expression=expression,
        isomiRs=isomiRs,
        intersection=intersection,
        union=union,
        jaccard=jaccard
    )


def df_to_network(interactions):
    edges = interactions.rename(columns={
        "isomir": "from",
        "target": "to"
    })
    
    nodes_from_size = edges[["from", "to"]].groupby("from").count().rename(columns={"to": "size"})
    nodes_to_size = edges[["from", "to"]].groupby("to").count().rename(columns={"from": "size"})
    nodes_size = pd.concat([nodes_from_size, nodes_to_size])

    nodes_from_median_tpm = (
        edges[["from", "isomir_median_tpm"]]
        .groupby("from").first()
        .rename(columns={"isomir_median_tpm": "median_tpm"})
    )
    nodes_to_median_tpm = (
        edges[["to", "target_median_tpm"]]
        .groupby("to").first()
        .rename(columns={"target_median_tpm": "median_tpm"})
    )
    nodes_median_tpm = pd.concat([nodes_from_median_tpm, nodes_to_median_tpm])
    nodes = nodes_size.join(nodes_median_tpm).reset_index().rename(columns={"index": "id"})
    nodes["median_tpm"] = (2**nodes["median_tpm"] - 1).round(1)

    nodes["label"] = nodes["id"]
    nodes["size"] = 24 + (nodes["size"] - 1) / (nodes["size"].max() - 1)*24*3
    
    nodes_dict = nodes.to_dict(orient="records")
    for i in range(len(nodes_dict)):
        if nodes_dict[i]["id"].startswith("hsa-"):
            nodes_dict[i]["font"] = {"size": 14*4}
            if nodes_dict[i]["id"].endswith("|0"):
                nodes_dict[i]["color"] = "#c7f464"
            else:
                nodes_dict[i]["color"] = "#ff6b6b"
            units = "RPM"
        else:
            nodes_dict[i]["color"] = "#95d4f3"
            units = "TPM"
        

        nodes_dict[i]["title"] = f"""{nodes_dict[i]["id"]}
        Median expression: {nodes_dict[i]["median_tpm"]} {units}"""
    
    edges_dict = edges.to_dict(orient="records")
    for i in range(len(edges_dict)):
        mirdb_score = int(edges_dict[i]["mirdb_score"]) if not pd.isna(edges_dict[i]["mirdb_score"]) else "-"
        targetscan_score = edges_dict[i]["targetscan_score"] if not pd.isna(edges_dict[i]["targetscan_score"]) else "-"
        edges_dict[i]["title"] = f"""From: {edges_dict[i]["from"]}
        To: {edges_dict[i]["to"]}
        miRDB score: {mirdb_score}
        TargetScan score: {targetscan_score}
        Spearman's correlation: {edges_dict[i]["spearman_corr"]}
        Spearman's p-value: {edges_dict[i]["spearman_p_value"]:.2e}
        Spearman's FDR: {edges_dict[i]["spearman_fdr"]:.2e}
        """
        for k in list(edges_dict[i].keys()):
            if k not in ["from", "to", "title"]:
                del edges_dict[i][k]
    
    return {"nodes": nodes_dict, "edges": edges_dict}


@frontend.route("/cancer/<cancer>", methods=["GET"])
def show_cancer(cancer):
    if not check_cancer(cancer):
        abort(400)

    interactions = get_significant_interactions(cancer)
    network = df_to_network(interactions)
    targets_summary = get_isomirs_targeting_summary_in_cancer(cancer)
    targets_summary = targets_summary.sort_values("isomir_median_tpm", ascending=False)
    targets_summary["isomir_median_tpm"] = (2**targets_summary["isomir_median_tpm"] - 1).round(1)

    return render_template(
        "cancer/main.html",
        page="cancer",
        cancer=cancer,
        interactions=interactions,
        network=network,
        targets_summary=targets_summary
    )


@frontend.route("/cancer_custom/<cancer>", methods=["GET"])
def show_cancer_custom(cancer):
    if cancer == "pan-cancer":
        network = get_pan_cancer_network()
        cancer = "COAD"  # Hot fix
        interactions = get_significant_interactions(cancer)
    else:
        interactions = pd.read_csv(f"web_app/api/cancer_custom/{cancer}.csv")
        interactions["isomir_median_tpm"] = 0
        interactions["target_median_tpm"] = 0
        interactions["mirdb_score"] = 0
        interactions["targetscan_score"] = 0
        interactions["spearman_corr"] = 0
        interactions["spearman_p_value"] = 0
        interactions["spearman_fdr"] = 0
        network = df_to_network(interactions)
    
    targets_summary = get_isomirs_targeting_summary_in_cancer(cancer)
    targets_summary = targets_summary.sort_values("isomir_median_tpm", ascending=False)
    targets_summary["isomir_median_tpm"] = (2**targets_summary["isomir_median_tpm"] - 1).round(1)
    
    return render_template(
        "cancer/main.html",
        page="cancer",
        cancer=cancer,
        interactions=interactions,
        network=network,
        targets_summary=targets_summary
    )


def html_p_value(p):
    if p == 0:
        return "0"

    float_str = '{0:.2e}'.format(p)
    base, exponent = float_str.split('e')
    
    if 2 >= int(exponent) >= -2:
        return '{:.3f}'.format(p)

    return r'{0} &#x2715; 10<sup>{1}</sup>'.format(base, int(exponent))


@frontend.route("/cancer_molecule/<cancer>/<molecule>", methods=["GET"])
def show_cancer_molecule(cancer, molecule):
    if not check_cancer(cancer) or not check_molecule(molecule):
        abort(400)
    
    # First, get expression data
    try:
        expression, highly_expressed = get_molecule_expression_in_cancer(molecule, cancer)
    except:
        abort(400)

    expression = [{
        "name": molecule,
        "x": expression.tolist(),
        "type": "violin"
    }]

    # Then, retrieve the results of isomiR target prediction
    if molecule.startswith("hsa-"):
        is_isomiR = True
        index_col = "target"
        targets = get_molecule_targeting_in_cancer(cancer, isomiR=molecule)
    else:
        is_isomiR = False
        index_col = "isomir"
        targets = get_molecule_targeting_in_cancer(cancer, target=molecule)
    
    targets["spearman_fdr_loc"] = multipletests(targets["spearman_p_value"], method="fdr_bh")[1]
    targets = targets.sort_values(["spearman_corr", "spearman_fdr_loc"]).set_index(index_col)
    targets["significant"] = (targets["spearman_corr"] < -0.3) & (targets["spearman_fdr_loc"] < 0.05)
    
    targets["spearman_p_value"] = [html_p_value(p) for p in targets["spearman_p_value"]]
    targets["spearman_fdr_glob"] = [html_p_value(p) for p in targets["spearman_fdr"]]
    targets["spearman_fdr_loc"] = [html_p_value(p) for p in targets["spearman_fdr_loc"]]
    
    targets["mirdb_score"] = targets["mirdb_score"].astype("Int64")
    targets["isomir_median_tpm"] = (2**targets["isomir_median_tpm"] - 1).round(1)
    targets["target_median_tpm"] = (2**targets["target_median_tpm"] - 1).round(1)
    targets = targets.astype("string").fillna("-")

    return render_template(
        "cancer_molecule/main.html",
        page="cancer_molecule",
        cancer=cancer,
        molecule=molecule,
        is_isomiR=is_isomiR,
        expression=expression,
        targets=targets
    )

    
@frontend.route("/cancer_isomir_target/<cancer>/<isomir>/<target>", methods=["GET"])
def show_cancer_isomir_target(cancer, isomir, target):
    if not check_cancer(cancer) or not check_molecule(isomir) or not check_molecule(target):
        abort(400)
    
    # Get expression arrays
    try:
        isomir_expr, highly_expressed = get_molecule_expression_in_cancer(isomir, cancer, units="tmm")
        target_expr, _ = get_molecule_expression_in_cancer(target, cancer, units="tmm")
    except:
        abort(400)
    
    # Get or calculate Spearman's correlation
    try:
        targeting = get_cancer_isomir_target(cancer, isomir, target)
        spearman_corr = targeting.spearman_corr
        spearman_p_value = targeting.spearman_p_value
    except:
        spearman_corr, spearman_p_value = spearmanr(isomir_expr, target_expr)
    
    spearman_corr = f"{spearman_corr:.2f}"
    spearman_p_value = html_p_value(spearman_p_value)
    
    # Calculate slope for regression line
    lr = linregress(isomir_expr, target_expr)
    x_min, x_max = isomir_expr.min(), isomir_expr.max()
    if x_min*lr.slope + lr.intercept < 0 and lr.slope != 0:
        x_min = -lr.intercept/lr.slope
    if x_max*lr.slope + lr.intercept < 0 and lr.slope != 0:
        x_max = -lr.intercept/lr.slope

    expression = [
        {
            "name": cancer,
            "x": isomir_expr.tolist(),
            "y": target_expr.tolist(),
            "type": "scatterplot",
            "mode": "markers"
        },
        {
            "x": [x_min, x_max],
            "y": [
                x_min*lr.slope + lr.intercept,
                x_max*lr.slope + lr.intercept
            ],
            "type": "scatterplot",
            "mode": "lines"
        }
    ]
    
    return render_template(
        "cancer_isomir_target/main.html",
        page="cancer_isomir_target",
        cancer=cancer,
        isomir=isomir,
        target=target,
        expression=expression,
        spearman_corr=spearman_corr,
        spearman_p_value=spearman_p_value
    )
