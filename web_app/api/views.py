from . import api
from web_app.api.db_functions import *
from web_app.api.models import Expression, Molecules

from flask import make_response, jsonify

import pandas as pd


@api.route("/cancer_custom/<cancer>/network", methods=["GET"])
def get_cancer_custom_network(cancer):
    interactions = pd.read_csv(f"web_app/api/cancer_custom/{cancer}.csv")
    interactions = interactions.rename(columns={
        "isomir": "from",
        "target": "to"
    })
    
    nodes_from_size = interactions[["from", "to"]].groupby("from").count().rename(columns={"to": "size"})
    nodes_to_size = interactions[["from", "to"]].groupby("to").count().rename(columns={"from": "size"})
    nodes_size = pd.concat([nodes_from_size, nodes_to_size])
    nodes = nodes_size.reset_index().rename(columns={"index": "id"})

    nodes["label"] = nodes["id"]
    nodes["size"] = 24 + (nodes["size"] - 1) / (nodes["size"].max() - 1)*24*3
    nodes_dict = nodes.to_dict(orient="records")
    for i in range(len(nodes_dict)):
        if nodes_dict[i]["id"].startswith("hsa-"):
            nodes_dict[i]["font"] = {"size": 14*4}
            if nodes_dict[i]["id"].endswith("|0"):
                nodes_dict[i]["color"] = "#c7f464"
            else:
                nodes_dict[i]["color"] = "#ff6b6b"
        else:
            nodes_dict[i]["color"] = "#95d4f3"
    
    edges_dict = interactions.to_dict(orient="records")
    
    resp = {"nodes": nodes_dict, "edges": edges_dict}

    return jsonify(resp)
