from . import api
from web_app import db
from web_app.api.models import *

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


def get_isomiR_targets_pan_cancer(isomiR):
    query = Targets_raw.query.filter(Targets_raw.isomir == isomiR).statement
    df = pd.read_sql(query, db.engine)
    return df


def get_isomiR_highle_expressed_cancers(isomiR):
    query = (
        Expression.query
        .with_entities(Expression.cancer, Expression.highly_expressed)
        .filter(Expression.molecule == isomiR).statement
    )
    df = pd.read_sql(query, db.engine).set_index("cancer")
    return df
