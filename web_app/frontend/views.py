from . import frontend
from web_app import app
from web_app.api.db_functions import *


from flask import render_template
from flask import request
from flask import Response
from flask import jsonify

#import pandas as pd
#import pickle as pkl
#from itertools import product
#import os
#import json


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
        - Pan-cancer expression distribution
        - Predicted miRDB/TargetScan targets
        - Pan-cancer anti-correlation patterns
        - TODO: cancer-universal targets
        - TODO: links to the neighboring molecules
    '''
    # TODO: check that molecule is in DB
    if molecule.startswith("hsa-"):
        is_isomiR = True
        targets_pan_cancer = get_isomiR_targets_pan_cancer(molecule)
        # Annotate cancers, where isomiR is highly expressed 
        isomiR_he = get_isomiR_highle_expressed_cancers(molecule)
        targets_pan_cancer = targets_pan_cancer.set_index("cancer").join(isomiR_he).reset_index()
        targets_pan_cancer["cancer"] = [
            row["cancer"] + ("*" if row["highly_expressed"] else "")
            for _, row in targets_pan_cancer.iterrows()
        ] 

        targets_seq = targets_pan_cancer.set_index("target")[["mirdb_score", "targetscan_score"]]
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
        for target, df in targets_pan_cancer.loc[targets_pan_cancer["spearman_corr"] < -0.3].groupby("target"):
            cancers = df["cancer"].to_list()
            if len(cancers) > 1:
                targets_conserved.append((target, sorted(cancers)))

        targets_conserved = sorted(targets_conserved, key=lambda x: len(x[1]), reverse=True)
    else:
        is_isomiR = False
        targets_seq = None
        targets_conserved = None
        pass

    return render_template(
        "molecule/main.html",
        molecule=molecule,
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

    return render_template("miRNA/main.html")
