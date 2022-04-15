from . import api
from web_app import db
from web_app.api.models import *

from sqlalchemy import and_, or_, distinct
from flask import make_response, jsonify

import pandas as pd
import numpy as np


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
    with corr < -0.3, p < 0.05 and highly expressed isomiRs
    '''
    
    filt = and_(
        Targets_raw.cancer == cancer,
        Targets_raw.spearman_corr < -0.3,
        Targets_raw.spearman_p_value < 0.05
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
