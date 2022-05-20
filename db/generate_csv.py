import os
import csv

import pandas as pd
import numpy as np

from tqdm import tqdm
from statsmodels.stats.multitest import multipletests


#################
# 1. Generate core tables: Molecules, Cancers, Expression
#################
'''
expr_path = "raw_text_tables/expression/TMM-RPM/"
he_path = "raw_text_tables/expression/highly_expressed_isomiRs/"
dfs = []
for fn in tqdm(os.listdir(expr_path)):
    if "tumor" not in fn:
        continue
    
    df = pd.read_csv(f"{expr_path}/{fn}", sep="\t", index_col=0)
    df.index.name = "molecule"
    df.columns = ["-".join(c.split("-")[1:-1]) for c in df.columns]
    
    df_isomiR = (2**df.loc[df.index.str.startswith("hsa-")] - 1).copy()
    df_gene = (2**df.loc[~df.index.str.startswith("hsa-")] - 1).copy()
    df_isomiR = np.log2(df_isomiR / df_isomiR.sum(axis=0) * 1e+6 + 1)
    df_gene = np.log2(df_gene / df_gene.sum(axis=0) * 1e+6 + 1)
    TPM = pd.concat([df_isomiR, df_gene])
    TPM = TPM.loc[df.index]

    df = df.round(1).astype("string")
    TPM = TPM.round(1).astype("string")
    
    df["tmm"] = "{" + df[df.columns].agg(",".join, axis=1) + "}"
    df["tpm"] = "{" + TPM[TPM.columns].agg(",".join, axis=1) + "}"
    df["cancer"] = fn.split("_")[0][5:]
    
    he = pd.read_csv(f"{he_path}/{fn}", sep="\t", index_col=0)
    he["highly_expressed"] = True
    df = df.join(he).fillna(False)
    
    df = df.reset_index()[["molecule", "cancer", "tmm", "tpm", "highly_expressed"]]
    dfs.append(df)

Expression = pd.concat(dfs)
Molecules = pd.DataFrame({"molecule": Expression["molecule"].unique()})
Cancers = pd.DataFrame({"cancer": Expression["cancer"].unique()})

for fn in ["Molecules", "Cancers", "Expression"]:
    eval(fn).to_csv(f"postgresql_data/csv/{fn}.csv", index=None, quoting=csv.QUOTE_ALL)
'''
#################
# 2. Generate target tables
#################

expr = pd.read_csv("postgresql_data/csv/Expression.csv")
expr["tpm"] = expr["tpm"].str[1:-1]
expr["median_tpm"] = [np.median([float(v) for v in t.split(",")]) for t in tqdm(expr["tpm"])]
expr.index = expr["molecule"] + "," + expr["cancer"]
expr = expr[["median_tpm"]]

miRDB_path = "raw_text_tables/miRDB"
TargetScan_path = "raw_text_tables/TargetScan"
dfs = []
for fn in tqdm(os.listdir(f"{miRDB_path}/corr_analysis")):
    if "tumor" not in fn or "raw" not in fn:
        continue
    fn2 = fn.replace("_raw", "")
    
    df1 = pd.read_csv(f"{miRDB_path}/corr_analysis/{fn}", sep="\t")
    df2 = pd.read_csv(f"{miRDB_path}/predicted_targets/{fn2}", sep="\t")
    df1.index = df1["isomiR"] + "," + df1["gene"]
    df2.index = df2["isomiR"] + "," + df2["Gene Symbol"]
    miRDB = df1[["corr", "p_value"]].join(df2[["Target Score"]])
    miRDB = miRDB.rename(columns={
        "corr": "spearman_corr_miRDB",
        "p_value": "spearman_p_value_miRDB",
        "Target Score": "mirdb_score"
    })
    
    df1 = pd.read_csv(f"{TargetScan_path}/corr_analysis/{fn}", sep="\t")
    df2 = pd.read_csv(f"{TargetScan_path}/predicted_targets/{fn2}", sep="\t")
    df1.index = df1["isomiR"] + "," + df1["gene"]
    df2.index = df2["isomiR"] + "," + df2["Gene Symbol"]
    TargetScan = df1[["corr", "p_value"]].join(df2[["CWCS"]])
    TargetScan = TargetScan.rename(columns={
        "corr": "spearman_corr_TargetScan",
        "p_value": "spearman_p_value_TargetScan",
        "CWCS": "targetscan_score"
    })
    
    df = miRDB.join(TargetScan, how="outer")
    miRDB_only = df["targetscan_score"].isna()
    df["spearman_corr"] = df["spearman_corr_TargetScan"]
    df.loc[miRDB_only, "spearman_corr"] = df.loc[miRDB_only, "spearman_corr_miRDB"]
    df["spearman_p_value"] = df["spearman_p_value_TargetScan"]
    df.loc[miRDB_only, "spearman_p_value"] = df.loc[miRDB_only, "spearman_p_value_miRDB"]
   
    df["isomir"] = df.index.str.split(",").str[0]
    df["target"] = df.index.str.split(",").str[1]
    df["cancer"] = fn.split("_")[0][5:]
    
    df.index = df["isomir"] + "," + df["cancer"]
    df = df.join(expr.rename(columns={"median_tpm": "isomir_median_tpm"}))
    df.index = df["target"] + "," + df["cancer"]
    df = df.join(expr.rename(columns={"median_tpm": "target_median_tpm"}))
    df = df.reset_index()
    del df["index"]
    
    df["spearman_fdr"] = multipletests(df["spearman_p_value"], method="fdr_bh")[1]

    df = df[["isomir", "target", "cancer", "mirdb_score", "targetscan_score", "spearman_corr", "spearman_p_value", "spearman_fdr", "isomir_median_tpm", "target_median_tpm"]]
    df["mirdb_score"] = df["mirdb_score"].astype("Int64").astype("string")
    df["targetscan_score"] = df["targetscan_score"].round(1).astype("string")
    df["spearman_corr"] = df["spearman_corr"].round(2).astype("string")
    df["isomir_median_tpm"] = df["isomir_median_tpm"].round(1).astype("string")
    df["target_median_tpm"] = df["target_median_tpm"].round(1).astype("string")
    df = df.fillna("NULL")
    dfs.append(df)

TargetsRaw = pd.concat(dfs)
TargetsRaw.to_csv(f"postgresql_data/csv/Targets_raw.csv", index=None)
