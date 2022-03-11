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
    # Load table with available COVID-19 lineages
    
    return render_template("index.html")
