from . import api
from web_app import db
from web_app.api.models import *

from sqlalchemy import and_, or_, not_, distinct
from flask import make_response, jsonify

import pandas as pd
import numpy as np


def check_molecule(molecule):
    if Molecules.query.filter(Molecules.molecule == molecule).count() == 1:
        return True
    else:
        return False

def check_miRNA(miRNA):
    if Molecules.query.filter(Molecules.molecule.startswith(miRNA)).count() > 0:
        return True
    else:
        return False


def check_cancer(cancer):
    if Cancers.query.filter(Cancers.cancer == cancer).count() == 1:
        return True
    else:
        return False


def get_all_isomiRs():
    results = Molecules.query.filter(Molecules.molecule.startswith("hsa")).all()
    return sorted([row.molecule for row in results], key=lambda x: (x.split("|")[0], int(x.split("|")[1])))


def get_all_genes():
    results = Molecules.query.filter(not_(Molecules.molecule.startswith("hsa"))).order_by(Molecules.molecule).all()
    return [row.molecule for row in results]


def get_all_cancers():
    return [row.cancer for row in Cancers.query.order_by(Cancers.cancer).all()]


def get_molecule_expression_pan_cancer(molecule, units="tpm", add_asterisk=True):
    '''
        Returns json {"cancer": [...], "expression": [...]}
    '''
    results = (
        Expression.query
        .filter(Expression.molecule == molecule)
        .order_by(Expression.cancer)
        .all()
    )
    dfs = []
    for res in results:
        df = pd.DataFrame(columns=["cancer", "expression"])
        if units == "tpm":
            df["expression"] = res.tpm
        else:
            df["expression"] = res.tmm

        df["cancer"] = res.cancer + ("*" if add_asterisk and res.highly_expressed else "")
        dfs.append(df)
    df = pd.concat(dfs)
    df["expression"] = df["expression"].astype("float64")

    return df


def get_miRNA_expression_pan_cancer(miRNA):
    '''
        Returns json {"cancer": [...], "expression": [...]}
    '''
    results = Molecules.query.filter(Molecules.molecule.startswith(miRNA)).all()
    isomiRs = sorted([i.molecule for i in results], key=lambda x: int(x.split("|")[1]))
    dfs = []
    for molecule in isomiRs:
        df = get_molecule_expression_pan_cancer(molecule, add_asterisk=False)
        df["expression"] = 2**df["expression"] - 1
        dfs.append((molecule, df))

    return dfs


def get_molecule_expression_in_cancer(molecule, cancer, units="tpm"):
    result = (
        Expression.query
        .filter(and_(
            Expression.molecule == molecule,
            Expression.cancer == cancer
        ))
    ).one()
    if units == "tpm":
        return np.array(result.tpm, dtype=np.float64), result.highly_expressed
    else:
        return np.array(result.tmm, dtype=np.float64), result.highly_expressed


def get_molecule_targeting_pan_cancer(isomiR=None, target=None, add_asterisk=True):
    if not isomiR and not target:
        raise Exception("isomiR or target should be passed")
    
    filt = (Targets_raw.isomir == isomiR) if isomiR else (Targets_raw.target == target)

    query = (
        Targets_raw.query
        .filter(filt)
        .join(Expression, and_(
            Targets_raw.isomir == Expression.molecule,
            Targets_raw.cancer == Expression.cancer
        ))
        .with_entities(
            *list(Targets_raw.__table__.c) + [Expression.highly_expressed]
        )
        .statement
    )
    df = pd.read_sql(query, db.engine)
    
    # Mark cancers, where isomiR is highly expressed
    df["cancer"] = [
        row["cancer"] + ("*" if add_asterisk and row["highly_expressed"] else "")
        for _, row in df.iterrows()
    ] 
    
    return df


def get_isomiR_targets_unique(isomiR):
    result  = (
        Targets_raw.query
        .filter(Targets_raw.isomir == isomiR)
        .with_entities(distinct(Targets_raw.target))
    ).all()

    return [p[0] for p in result]

def get_molecule_targeting_in_cancer(cancer, isomiR=None, target=None):
    if not isomiR and not target:
        raise Exception("isomiR or target should be passed")
    
    filt = (Targets_raw.isomir == isomiR) if isomiR else (Targets_raw.target == target)

    query = (
        Targets_raw.query
        .filter(and_(
            filt,
            Targets_raw.cancer == cancer
        ))
        .join(Expression, and_(
            Targets_raw.isomir == Expression.molecule,
            Targets_raw.cancer == Expression.cancer
        ))
        .with_entities(
            *list(Targets_raw.__table__.c) + [Expression.highly_expressed]
        )
        .statement
    )
    df = pd.read_sql(query, db.engine)
    
    return df


def get_significant_interactions(cancer):
    '''
    Return list of isomiR->target interactions
    with corr < -0.3, FDR < 0.05 and highly expressed isomiRs
    '''
    
    filt = and_(
        Targets_raw.cancer == cancer,
        Targets_raw.spearman_corr < -0.3,
        Targets_raw.spearman_fdr < 0.05
    )
    
    query = (
        Targets_raw.query
        .filter(filt)
        .join(Expression, and_(
            Targets_raw.isomir == Expression.molecule,
            Targets_raw.cancer == Expression.cancer,
            Expression.highly_expressed == True
        ))
        .with_entities(
            *list(Targets_raw.__table__.c)
        )
        .statement
    )
    df = pd.read_sql(query, db.engine)

    return df


def get_isomirs_targeting_summary_in_cancer(cancer):
    query = (
        Targets_summary.query
        .filter(Targets_summary.cancer == cancer)
        .statement
    )
    df = pd.read_sql(query, db.engine)

    return df


def get_cancer_isomir_target(cancer, isomir, target):
    result = (
        Targets_raw.query
        .filter(and_(
            Targets_raw.cancer == cancer,
            Targets_raw.isomir == isomir,
            Targets_raw.target == target,
        ))
    ).one()

    return result


def get_pan_cancer_network():
    # Find the number of cancers for each interaction
    query = """
        SELECT C.isomir, C.target, C.n_corr, D.n_he FROM 
            (SELECT A.isomir, A.target, count(A.isomir) as n_corr FROM
                (SELECT * FROM targets_raw WHERE spearman_corr < -0.3 AND spearman_p_value < 0.05) A
                INNER JOIN
                (SELECT molecule, cancer FROM expression WHERE highly_expressed = TRUE) B
                ON A.isomir = B.molecule AND A.cancer = B.cancer
            GROUP BY (A.isomir, A.target) ORDER BY n_corr DESC) C
            INNER JOIN
            (SELECT molecule, count(cancer) AS n_he FROM expression WHERE highly_expressed = TRUE GROUP BY molecule) D
            ON C.isomir = D.molecule
    """
    df = pd.read_sql(query, db.engine)
    df = df.loc[df["n_he"] >= 10]
    df = df.loc[df["n_corr"] >= df["n_he"] / 2]
    
    edges = df.rename(columns={
        "isomir": "from",
        "target": "to"
    })
    edges["frac"] = df["n_corr"] / df["n_he"]
    edges["width"] = 4 + (edges["n_corr"] - edges["n_corr"].min()) / (edges["n_corr"].max() - edges["n_corr"].min())*13

    nodes = pd.DataFrame({"id": edges["from"].unique().tolist() + edges["to"].unique().tolist()})
    nodes["label"] = nodes["id"]
    nodes["size"] = 24
    
    nodes_from_size = edges[["from", "to"]].groupby("from").count().rename(columns={"to": "size"})
    nodes_to_size = edges[["from", "to"]].groupby("to").count().rename(columns={"from": "size"})
    nodes_size = pd.concat([nodes_from_size, nodes_to_size])
    nodes = nodes_size.reset_index().rename(columns={"index": "id"})
    nodes["label"] = nodes["id"]
    nodes["size"] = 24 + (nodes["size"] - nodes["size"].min()) / (nodes["size"].max() - nodes["size"].min())*12
    
    nodes_dict = nodes.to_dict(orient="records")
    for i in range(len(nodes_dict)):
        nodes_dict[i]["font"] = {"size": 14*1.5}
        nodes_dict[i]["font"] = {
            "size": 25 + (nodes_dict[i]["size"] - nodes["size"].min()) / (nodes["size"].max() - nodes["size"].min())*5
        }
        if nodes_dict[i]["id"].startswith("hsa-"):
            if nodes_dict[i]["id"].endswith("|0"):
                nodes_dict[i]["color"] = "#c7f464"
            else:
                nodes_dict[i]["color"] = "#ff6b6b"
        else:
            nodes_dict[i]["color"] = "#95d4f3"
    
    edges_dict = edges.to_dict(orient="records")
    
    return {"nodes": nodes_dict, "edges": edges_dict}
