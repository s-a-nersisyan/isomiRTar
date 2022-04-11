import sys
import os
import pandas as pd
from sqlalchemy import func, and_, or_

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append("/".join(script_dir.split("/")[:-1]))

from web_app import db
from web_app.api.models import *

db.create_all()

query1 = (
    Targets_raw.query
    .group_by(Targets_raw.isomir, Targets_raw.cancer)
    .with_entities(
        Targets_raw.isomir,
        Targets_raw.cancer,
        func.count(Targets_raw.target),
        func.max(Targets_raw.isomir_median_tpm)
    )
    .statement
)
df1 = pd.read_sql(query1, db.engine)
df1 = df1.rename(columns={"count_1": "num_targets_expressed", "max_1": "isomir_median_tpm"})

query2 = (
    Targets_raw.query
    .filter(and_(
        Targets_raw.spearman_corr < -0.3,
        Targets_raw.spearman_p_value < 0.05
    ))
    .group_by(Targets_raw.isomir, Targets_raw.cancer)
    .with_entities(
        Targets_raw.isomir,
        Targets_raw.cancer,
        func.count(Targets_raw.target),
    )
    .statement
)
df2 = pd.read_sql(query2, db.engine)
df2 = df2.rename(columns={"count_1": "num_targets_anticorrelated"})

query3 = (
    Expression.query
    .filter(Expression.molecule.startswith("hsa-"))
    .group_by(Expression.molecule, Expression.cancer, Expression.highly_expressed)
    .with_entities(
        Expression.molecule,
        Expression.cancer,
        Expression.highly_expressed,
    )
    .statement
)
df3 = pd.read_sql(query3, db.engine)

df1.index = df1.iloc[:, 0] + ";" + df1.iloc[:, 1]
df2.index = df2.iloc[:, 0] + ";" + df2.iloc[:, 1]
df3.index = df3.iloc[:, 0] + ";" + df3.iloc[:, 1]
df = df1.join(df2[["num_targets_anticorrelated"]], how="outer").join(df3[["highly_expressed"]], how="outer")
df["isomir"] = df.index.str.split(";").str[0]
df["cancer"] = df.index.str.split(";").str[1]
df = df.reset_index().drop(columns=["index"])
df = df[["isomir", "cancer", "num_targets_expressed", "num_targets_anticorrelated", "isomir_median_tpm", "highly_expressed"]]

df["num_targets_expressed"] = df["num_targets_expressed"].fillna(0).astype("int64")
df["num_targets_anticorrelated"] = df["num_targets_anticorrelated"].fillna(0).astype("int64")
df["isomir_median_tpm"] = df["isomir_median_tpm"].fillna(0).round(1).astype("string")
print(df)

df.to_sql("targets_summary", db.engine, index=False, if_exists="append", schema="public")
db.session.commit()
