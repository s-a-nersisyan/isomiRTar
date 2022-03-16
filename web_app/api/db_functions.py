from . import api
from web_app.api.models import Expression, Molecules

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
