from . import frontend
from web_app import app
from web_app.api.db_functions import *


from flask import render_template
from flask import abort

import json
from scipy.stats import linregress, spearmanr


@frontend.route("/", methods=["GET"])
def show_index():
    '''
    TODO
    '''
    return render_template("index.html")


@frontend.route("/molecule/<molecule>", methods=["GET"])
def show_molecule(molecule):
    '''
    Show information about a specified molecule:
        - TODO: links to the neighboring molecules
    '''
    # TODO: check that molecule is in DB
    # First, get molecule expression for boxplot
    df = get_molecule_expression_pan_cancer(molecule)
    order = df.groupby("cancer").quantile(0.75).sort_values("expression", ascending=False).index
    df = pd.concat([df.loc[df["cancer"] == cancer] for cancer in order])
    df["expression"] = 2**df["expression"] - 1
    
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
    targets_pan_cancer = targets_pan_cancer.loc[targets_pan_cancer["spearman_corr"] < -0.3]
    for target, df in targets_pan_cancer.groupby(index_col):
        cancers = df["cancer"].to_list()
        if len(cancers) > 1:
            if is_isomiR:
                targets_conserved.append((molecule, target, sorted(cancers)))
            else:
                targets_conserved.append((target, molecule, sorted(cancers)))

    targets_conserved = sorted(targets_conserved, key=lambda x: len(x[2]), reverse=True)

    return render_template(
        "molecule/main.html",
        molecule=molecule,
        expression=expression,
        is_isomiR=is_isomiR,
        targets_seq=targets_seq,
        targets_conserved=targets_conserved
    )


@frontend.route("/miRNA/<miRNA>", methods=["GET"])
def show_miRNA(miRNA):
    '''
    Show information about a specified miRNA:
        - Pan-cancer expression distribution
        - Predicted miRDB/TargetScan targets
        - Pan-cancer anti-correlation patterns
        - TODO: cancer-universal targets
        - TODO: links to the neighboring miRNAs
    '''
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

    return render_template(
        "miRNA/main.html",
        expression=expression
    )


def df_to_network(interactions):
    interactions = interactions.rename(columns={
        "isomir": "from",
        "target": "to"
    })
    
    nodes_from_size = interactions[["from", "to"]].groupby("from").count().rename(columns={"to": "size"})
    nodes_to_size = interactions[["from", "to"]].groupby("to").count().rename(columns={"from": "size"})
    nodes_size = pd.concat([nodes_from_size, nodes_to_size])

    nodes_from_median_tpm = (
        interactions[["from", "isomir_median_tpm"]]
        .groupby("from").first()
        .rename(columns={"isomir_median_tpm": "median_tpm"})
    )
    nodes_to_median_tpm = (
        interactions[["to", "target_median_tpm"]]
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
    
    edges_dict = interactions.to_dict(orient="records")
    for i in range(len(edges_dict)):
        mirdb_score = int(edges_dict[i]["mirdb_score"]) if not pd.isna(edges_dict[i]["mirdb_score"]) else "-"
        targetscan_score = edges_dict[i]["targetscan_score"] if not pd.isna(edges_dict[i]["targetscan_score"]) else "-"
        edges_dict[i]["title"] = f"""From: {edges_dict[i]["from"]}
        To: {edges_dict[i]["to"]}
        miRDB score: {mirdb_score}
        TargetScan score: {targetscan_score}
        Spearman corr: {edges_dict[i]["spearman_corr"]}
        Spearman p-value: {edges_dict[i]["spearman_p_value"]:.2e}
        """
        for k in list(edges_dict[i].keys()):
            if k not in ["from", "to", "title"]:
                del edges_dict[i][k]
    
    return {"nodes": nodes_dict, "edges": edges_dict}


@frontend.route("/cancer/<cancer>", methods=["GET"])
def show_cancer(cancer):
    interactions = get_significant_interactions(cancer)
    network = df_to_network(interactions)
    targets_summary = get_isomirs_targeting_summary_in_cancer(cancer)
    targets_summary = targets_summary.sort_values("isomir_median_tpm", ascending=False)
    targets_summary["isomir_median_tpm"] = (2**targets_summary["isomir_median_tpm"] - 1).round(1)

    return render_template(
        "cancer/main.html",
        interactions=interactions,
        network=network,
        targets_summary=targets_summary
    )


@frontend.route("/cancer_custom/<cancer>", methods=["GET"])
def show_cancer_custom(cancer):
    return render_template("cancer_custom/main.html")


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
    '''
    Show information about a specified miRNA:
        - Pan-cancer expression distribution
        - Predicted miRDB/TargetScan targets
        - Pan-cancer anti-correlation patterns
        - TODO: cancer-universal targets
        - TODO: links to the neighboring miRNAs
    '''
    
    # First, get expression data
    expression, highly_expressed = get_molecule_expression_in_cancer(molecule, cancer)
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
    
    targets = targets.sort_values("spearman_corr").set_index(index_col)
    targets["significant"] = (targets["spearman_corr"] < -0.3) & (targets["spearman_p_value"] < 0.05)
    targets["spearman_p_value"] = [html_p_value(p) for p in targets["spearman_p_value"]]
    targets["mirdb_score"] = targets["mirdb_score"].astype("Int64")
    targets["isomir_median_tpm"] = (2**targets["isomir_median_tpm"] - 1).round(1)
    targets["target_median_tpm"] = (2**targets["target_median_tpm"] - 1).round(1)
    targets = targets.astype("string").fillna("-")

    return render_template(
        "cancer_molecule/main.html",
        is_isomiR=is_isomiR,
        expression=expression,
        targets=targets
    )

    
@frontend.route("/cancer_isomir_target/<cancer>/<isomir>/<target>", methods=["GET"])
def show_cancer_isomir_target(cancer, isomir, target):
    '''
    Show information about a specified miRNA:
        - Target Scores
        - Scatterplots
    '''
    
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
        cancer=cancer,
        isomir=isomir,
        target=target,
        expression=expression,
        spearman_corr=spearman_corr,
        spearman_p_value=spearman_p_value
    )
