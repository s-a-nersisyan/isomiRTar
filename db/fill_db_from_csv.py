import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append("/".join(script_dir.split("/")[:-1]))

from web_app import db

#db.drop_all()
db.create_all()

#for fn in ["Molecules", "Cancers", "Expression"]:
for fn in ["Targets_raw"]:
    print(fn)
    res = db.session.execute(f"""
    COPY {fn} FROM 'csv/{fn}.csv' WITH (FORMAT CSV, HEADER TRUE, NULL 'NULL');
    """)
    db.session.commit()
