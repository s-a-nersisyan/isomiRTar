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
    df = df.round(1).astype("string")
    df["expression"] = "{" + df[df.columns].agg(",".join, axis=1) + "}"
    df["cancer"] = fn.split("_")[0][5:]
    
    he = pd.read_csv(f"{he_path}/{fn}", sep="\t", index_col=0)
    he["highly_expressed"] = True
    df = df.join(he).fillna(False)
    
    df = df.reset_index()[["molecule", "cancer", "expression", "highly_expressed"]]
    dfs.append(df)

Expression = pd.concat(dfs)
Molecules = pd.DataFrame({"molecule": Expression["molecule"].unique()})
Cancers = pd.DataFrame({"cancer": Expression["cancer"].unique()})

for fn in ["Molecules", "Cancers", "Expression"]:
    eval(fn).to_csv(f"postgresql_data/csv/{fn}.csv", index=None, quoting=csv.QUOTE_ALL)
