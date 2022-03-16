from . import api
from web_app.api.db_functions import *
from web_app.api.models import Expression, Molecules

from flask import make_response, jsonify

import pandas as pd


@api.route("/molecule/<molecule>/expression", methods=["GET"])
def get_molecule_expression(molecule):
    '''
        Return csv table with pan-cancer molecule expression distribution:
        cancer,sample,expression
    '''
    # TODO: check that molecule is in the list of molecules
    df = get_molecule_expression_pan_cancer(molecule)
    order = df.groupby("cancer").quantile(0.75).sort_values("expression", ascending=False).index
    df = pd.concat([df.loc[df["cancer"] == cancer] for cancer in order])
    df["expression"] = 2**df["expression"] - 1
    
    resp = df.rename(columns={"cancer": "x", "expression": "y"}).to_dict(orient="list")
    resp["type"] = "box"
    resp = [resp]

    return jsonify(resp)


@api.route("/miRNA/<miRNA>/expression", methods=["GET"])
def get_miRNA_expression(miRNA):
    '''
    TODO
    '''
    # TODO: check that molecule is in the list of molecules
    results = Molecules.query.filter(Molecules.molecule.startswith(miRNA)).all()
    isomiRs = sorted([i.molecule for i in results], key=lambda x: int(x.split("|")[1]))
    resp = []
    for molecule in isomiRs:
        df = get_molecule_expression_pan_cancer(molecule, add_asterisk=False)
        df["expression"] = 2**df["expression"] - 1
        resp.append({
            "df": df,
            "name": molecule,
            "type": "box",
            "orientation": "h",
            "max_q75": df.groupby("cancer").quantile(0.75)["expression"].max()
        })
    
    # Get the most abundant isomiR
    i_max = sorted(list(range(len(resp))), key=lambda i: resp[i]["max_q75"])[-1]
    # Sort cancers
    order = resp[i_max]["df"].groupby("cancer").quantile(0.75).sort_values("expression", ascending=True).index
    for i in range(len(resp)):
        df = resp[i]["df"]
        df = pd.concat([
            df.loc[df["cancer"] == cancer] if cancer in list(df["cancer"]) else pd.DataFrame({"cancer": [cancer], "expression": [0.0]})
            for cancer in order
        ])
        del resp[i]["df"]
        resp[i] = {
            **resp[i],
            **df.rename(columns={"cancer": "y", "expression": "x"}).to_dict(orient="list"),
            "num_cancers": len(order)
        }
    
    return jsonify(resp)
