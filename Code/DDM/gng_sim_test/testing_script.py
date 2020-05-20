"""
Created on Thu Apr 30 15:33:03 2020

@author: natehall
"""

## Read dependencies
import hddm
import os
#import pandas as pd
#from pandas import Series
#from matplotlib.backends.backend_agg import FigureCanvasAgg
#import matplotlib.pyplot as plt
import pymp
from patsy import dmatrix
#import kabuki
#from kabuki.analyze import gelman_rubin
#import pickle
import numpy as np
#import sys
import argparse
import os.path


# test loading of packages in hddm_local venv and execution of parallel looping in pymp.

x = ['v', 'v_st']#, 'v_sv', 'v_sv_st']
nchains = 5

##configure pymp for parallel processing
pymp.config.nested=True

with pymp.Parallel(len(x)) as p:
    for index in range(0, len(mod_dict_torun)):
        with pymp.Parallel(nchains) as ch:
            for ch_index in ch.range(0,nchains):
                print(x[index])
                print(nchains[ch_index])
