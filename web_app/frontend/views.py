from . import frontend
from web_app import app
#from web_app.api.models import HLAAllelesPeptides


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


@frontend.route("/isomiR/<isomiR>", methods=["GET"])
def show_isomiR(isomiR):
    '''
    Show information about a specified isomiR:
        - Pan-cancer expression distribution
        - Predicted miRDB/TargetScan targets
        - Pan-cancer anti-correlation patterns
        - TODO: cancer-universal targets
        - TODO: links to the neighboring isomiRs
    '''

    return render_template("isomiR.html")
