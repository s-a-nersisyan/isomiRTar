from flask import session
import pandas as pd
import numpy as np

import sys
import os
from tqdm.auto import tqdm
from sqlalchemy.exc import NoResultFound

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append("/".join(script_dir.split("/")[:-1]))

from web_app import db
from web_app.api.models import Targets, Genes, Isomirnas, Comparasion_ts_1, Comparasion_mirdb_1


print("Reading comparasion file")
comparasion = pd.read_csv("comparasion.tsv", sep='\t', index_col=0)

print("Reading mirdb files")
mirdb_res_full = pd.read_csv("mirdb_res_full.tsv", sep='\t', index_col=0)
mirdb_res_80_upper = pd.read_csv("mirdb_res_80_upper.tsv", sep='\t', index_col=0)
mirdb_comparasion_1 = pd.read_csv("mirdb_comparasion_1.tsv", sep='\t', index_col=0)

print("Reading TargetScan files")
ts_res_full = pd.read_csv("ts_res_full.tsv", sep='\t', index_col=0)
ts_res_80_upper = pd.read_csv("ts_res_80_upper.tsv", sep='\t', index_col=0)
ts_comparasion_1 = pd.read_csv("ts_comparasion_1.tsv", sep='\t', index_col=0)

print("Creating DataBase shema")
db.drop_all()
db.create_all()

genes = sorted(set(mirdb_res_full['Gene Symbol'].dropna()) | set(ts_res_full['Gene Symbol'].dropna()))
isomirnas = comparasion.index.unique()

print(f"isomirnas len = {len(isomirnas)}, genes len = {len(genes)}")

print("Inserting genes table")
objects = [
    Genes(gene=gene) for gene in tqdm(genes)
]
db.session.bulk_save_objects(objects)
db.session.commit()

print("Inserting isomirnas, Comparasion_ts_1 and Comparasion_mirdb_1 tables")

mirdb_groupby = mirdb_res_full.groupby("isomiRNA")
ts_groupby = ts_res_full.groupby("isomiRNA")

def nan_to_int(num):
    return 0 if not num or np.isnan(num) else int(num)

def nan_to_float(num):
    return 0.0 if not num or np.isnan(num) else float(num)

objects = []
for isomirna, (ts_full, ts_short) in tqdm(comparasion.iterrows()):
    ts_tmp = ts_comparasion_1.loc[isomirna].copy()
    mirdb_tmp = mirdb_comparasion_1.loc[isomirna].copy()

    isom = Isomirnas(
        isomirna=isomirna,
        ts_full_mirdb_short_comparasion=ts_full,
        ts_short_mirdb_short_comparasion=ts_short
    )
    ts_comp = Comparasion_ts_1(
        isomirna = isom.isomirna,
        ts_count = nan_to_int(ts_tmp["isomiRNA|count"]),
        ts_count_1 = nan_to_int(ts_tmp["isomiRNA|+1|count"]),
        ts_intersection_1 = nan_to_float(ts_tmp["isomiRNA|+1|intersection"]),
        ts_union_1 = nan_to_float(ts_tmp["isomiRNA|+1|union"]),
    )
    mirdb_comp = Comparasion_mirdb_1(
        isomirna = isom.isomirna,
        mirdb_count = nan_to_int(mirdb_tmp["isomiRNA|count"]),
        mirdb_count_1 = nan_to_int(mirdb_tmp["isomiRNA|+1|count"]),
        mirdb_intersection_1 = nan_to_float(mirdb_tmp["isomiRNA|+1|intersection"]),
        mirdb_union_1 = nan_to_float(mirdb_tmp["isomiRNA|+1|union"]),
    )

    objects.extend([isom, ts_comp, mirdb_comp])
db.session.bulk_save_objects(objects)
db.session.commit()

print("Inserting targets of ts")
objects = []
for isomirna, (cwcs, gene, _) in tqdm(ts_res_full.iterrows()):
    try:
        objects.append(
            Targets(
                isomirna = db.session.query(Isomirnas.isomirna).filter_by(isomirna=isomirna).one()[0],
                gene = db.session.query(Genes.gene).filter_by(gene=gene).one()[0],
                target_score = 0,
                cwcs = cwcs,
            )
        )
    except NoResultFound:
        print(f"Can't find {isomirna} in isomirna table")
        continue
db.session.bulk_save_objects(objects)
db.session.commit()

print("Inserting targets of mirdb")
objects = []
for isomir, (target_score, gene, _) in tqdm(mirdb_res_full.iterrows()):
    target = Targets.query.filter(db.and_(
        Targets.isomirna == isomirna,
        Targets.gene == gene
        )).one_or_none()
    if target:
        target.target_score = target_score
    else:
        target = Targets(
            isomirna = db.session.query(Isomirnas.isomirna).filter_by(isomirna=isomirna).one()[0],
            gene = db.session.query(Genes.gene).filter_by(gene=gene).one()[0],
            target_score = target_score,
            cwcs = 0,
        )
    objects.append(target)

db.session.bulk_save_objects(objects)
db.session.commit()

print('Check select time')
import time
start = time.time()
db.session.query(Targets).filter(Targets.cwcs < -1).first()
print(time.time() - start)
