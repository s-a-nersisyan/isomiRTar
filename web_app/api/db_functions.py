from . import api
from web_app import db
from web_app.api.models import *

from sqlalchemy import and_, or_
from flask import make_response, jsonify

import pandas as pd


def get_molecule_expression_pan_cancer(molecule, add_asterisk=True):
    '''
        Returns json {"cancer": [...], "expression": [...]}
    '''
    results = Expression.query.filter(Expression.molecule == molecule).all()
    dfs = []
    for res in results:
        df = pd.DataFrame(columns=["cancer", "expression"])
        df["expression"] = res.tpm
        df["cancer"] = res.cancer + ("*" if add_asterisk and res.highly_expressed else "")
        dfs.append(df)
    df = pd.concat(dfs)
    df["expression"] = df["expression"].astype("float64")

    return df


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
