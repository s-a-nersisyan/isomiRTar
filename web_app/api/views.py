from . import api
from web_app.api.models import Expression

from flask import make_response

import pandas as pd


@api.route("/molecule/<molecule>/expression", methods=["GET"])
def get_molecule_expression(molecule):
    '''
        Return csv table with pan-cancer molecule expression distribution:
        cancer,sample,expression
    '''
    # TODO: check that molecule is in the list of molecules
    # TODO: perhaps write some function to do this automatically?
    results = Expression.query.filter(Expression.molecule == molecule).all()
    dfs = []
    for res in results:
        df = pd.DataFrame(columns=["molecule", "cancer", "expression"])
        df["expression"] = res.tpm
        df["molecule"] = molecule
        df["cancer"] = res.cancer + ("*" if res.highly_expressed else "")
        dfs.append(df)
    df = pd.concat(dfs)
    df["expression"] = df["expression"].astype("float64")

    order = df[["cancer", "expression"]].groupby("cancer").median().sort_values("expression", ascending=False).index
    df = pd.concat([df.loc[df["cancer"] == cancer] for cancer in order])
    df["expression"] = 2**df["expression"] - 1
    
    return df.to_csv(index=None)

