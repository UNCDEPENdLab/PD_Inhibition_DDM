#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 12:45:27 2020

@author: natehall
"""

######################################################################
## test reading args on ACI as I haven't done this before. 
######################################################################
##Read dependencies 
#import hddm
#import os
#import pandas as pd
#from pandas import Series
#from matplotlib.backends.backend_agg import FigureCanvasAgg
#import matplotlib.pyplot as plt 
#import pymp
#from patsy import dmatrix
#import kabuki
#from kabuki.analyze import gelman_rubin
#import pickle
#import numpy as np
#import sys
import argparse
#######################

## in case the user inputs nothing for the models they would like to have tested. These are the defaults
supported_flanker = ['v', 'vsv', 'v_block', 'v_blocksv', 'vst', 'vsvst', 'v_blockst', 'v_blocksvst']
supported_recent_probes = ['v', 'vsv', 'vst', 'vsvst']


# setup parser and parse inputs

parser = argparse.ArgumentParser(description='Run Hierarchical Bayesian Drift Diffusion Model on experimental task data: http://ski.clps.brown.edu/hddm_docs/index.html')
parser.add_argument('rawdir', nargs=1, type=str, help = "path to raw dataframe")
parser.add_argument('outputdir', nargs=1, type=str, help = "path to output directory")
parser.add_argument('task', nargs =1, help = "Experimental task to run HDDM on. The models to be run are coded into the function, so to get it to work for your needs you may need to mess with the HDDM  regressor calls within the function itself.", type = str)
#if par 
parser.add_argument('-models', '-m', nargs='+', help = "Different parameterizations of the DDM models to run. Currently supports " + str(supported_recent_probes) + " for recent probes and " + str(supported_flanker) + " for flanker. If this argument is left blank, will evaluate all of these. More coming soon.")
parser.add_argument('-code', '-c', help="coding scheme: stimulus(stim) or accuracy(acc)", type= str, default = "acc")
parser.add_argument('-nchains', '-nc', help="Number of unique MCMC chains to run. N.B. these are run in parallel so be careful the requested processors matches. Defaults to 5", type= int, default= 5)
parser.add_argument('-nsamp', '-ns', help = "total number of draws from MCMC posterior. Defaults to 2000", type = int, default = 2000)
parser.add_argument('-nburn', '-nb', help = "number of burn in samples. Defaults to 500", default = 500)

args = parser.parse_args()

# reassign parsed arguments to variables outside of args

if args.models==None:
    print(args.task)
#    print(str(args.task))
    if str(args.task)=="['flanker']":
        models = supported_flanker
        task = 'flanker'
    elif str(args.tas)=="['recent_probes']":
        models = supported_recent_probes
        task = 'recent_probes'
elif str(args.task)=="['flanker']":
    task = 'flanker'
elif str(args.task)=="['recent_probes']":
    task = 'recent_probes'
        
print 'Processing ' + task + 'according to these parameterizations: ' + str(models)      
rawdir = str(args.rawdir).strip('[]')
print 'Input: ' + rawdir
outputdir = str(args.outputdir).strip('[]')
print 'Outputdir: ' + outputdir

code = args.code
print 'Response coding: ' + code

nchains = args.nchains
print 'Nchains: ' + str(nchains)

nsamp = args.nsamp
print 'Nsamples: ' + str(nsamp)

nburn = args.nburn
print 'Nburn: ' + str(nburn)

