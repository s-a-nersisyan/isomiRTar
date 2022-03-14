from . import api

from flask import make_response

import pandas as pd


@api.route("/isomiR/<isomiR>/expression", methods=["GET"])
def get_isomiR_expression(isomiR):
    '''
        Return csv table with pan-cancer isomiR expression distribution:
        cancer,sample,expression
    '''
    df = pd.read_csv("https://raw.githubusercontent.com/plotly/datasets/master/violin_data.csv")
    return df.to_csv(index=None)

