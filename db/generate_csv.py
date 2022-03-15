import os
import csv

import pandas as pd
import numpy as np

from tqdm import tqdm


#################
# 1. Generate core tables: Molecules, Cancers, Expression
#################

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
