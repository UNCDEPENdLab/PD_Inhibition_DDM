#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 12:45:27 2020

@author: natehall
"""
import argparse,  hddm, csv, os
from pathlib import Path
import pandas as pd
import itertools

## in case the user inputs nothing for the models they would like to have tested. These are the defaults
supported_flanker = ['v', 'vsv', 'v_block', 'v_blocksv', 'vst', 'vsvst', 'v_blockst', 'v_blocksvst']
supported_recent_probes = ['v', 'vsv', 'vst', 'vsvst']


# setup parser and parse inputs

parser = argparse.ArgumentParser(description='Run Hierarchical Bayesian Drift Diffusion Model on experimental task data: http://ski.clps.brown.edu/hddm_docs/index.html')
#parser.add_argument('rawdf', nargs=1, type=argparse.FileType('r'), help = "path to raw dataframe")
parser.add_argument('rawdf', nargs=1, type=str, help = "path to raw dataframe")
parser.add_argument('outputdir', nargs=1, type=str, help = "path to output directory")
parser.add_argument('task', nargs =1, help = "Experimental task to run HDDM on. The models to be run are coded into the function, so to get it to work for your needs you may need to mess with the HDDM  regressor calls within the function itself.", type = str)
#if par 
parser.add_argument('-models', '-m', nargs='+', help = "Different parameterizations of the DDM models to run. Currently supports " + str(supported_recent_probes) + " for recent probes and " + str(supported_flanker) + " for flanker. If this argument is left blank, will evaluate all of these. More coming soon.")
parser.add_argument('-code', '-c', help="coding scheme: stimulus(stim) or accuracy(acc)", type= str, default = "acc")
parser.add_argument('-nchains', '-nc', help="Number of unique MCMC chains to run. N.B. these are run in parallel so be careful the requested processors matches. Defaults to 5", type= int, default= 5)
parser.add_argument('-nsamp', '-ns', help = "total number of draws from MCMC posterior. Defaults to 2000", type = int, default = 2000)
parser.add_argument('-nburn', '-nb', help = "number of burn in samples. Defaults to 500", default = 500)




args = parser.parse_args()
print args

##str(args.rawdir).strip('[]')
#data = os.open(str(args.rawdf).strip('[]'))
#print data
#
#
#import os
#data = os.open('/Users/natehall/github_repos/PD_Inhibition_DDM/Data/preprocessed/flanker_accCode.csv', os.O_RDONLY)
#print data

#df = pd.read_csv('/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Data/preprocessed/flanker_accCode.csv', sep=',')
#print(df)

#rawdf = '/Users/natehall/ics/Nate/PD_Inhibition_DDM/Data/preprocessed/flanker_accCode.csv'

rawdf = str(args.rawdf).strip('[]')
#print rawdf + '\n'
#
#rawdf.replace("'", "")
#print rawdf

rawdf = rawdf.replace("'", "")
print rawdf

data = hddm.load_csv(rawdf)
print data
#print '\n'
#print type(rawdf).__name__ + '\n'
#
#
#print Path(rawdf).is_file()
#
#print os.path.isfile(rawdf)
#
#rawdf_new = '/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Data/preprocessed/flanker_accCode.csv'
#print rawdf_new 
#
#print '\n'
#print type(rawdf).__name__ + '\n'
#
#
#print Path(rawdf_new).is_file()
#
#print os.path.isfile(rawdf_new)
#
#
#df = pd.read_csv(rawdf_new, sep=',')
#print(df)

df = pd.read_csv(rawdf, sep=',')
print(df)



#df = pd.read_csv(rawdf, sep=',')
#print(df)

#def csvParser(f):
#   with f:
#       csv.reader(f)
#
#data = csvParser(args.rawdf)
#print data

#data = pd.read_csv(args.rawdir)
#print data
#
#csv_reader = csv.reader(
#    io.open(args.csv_file_path, "r", encoding=args.encoding),
#    delimiter=",",
#    quotechar='"'
#)

#
#with open( args.rawdir ) as csvFile:
#    reader = rawdir.reader(csvFile, delimiter=',')
#
#data = csvParser(rawdir)

#hddm.load_csv(rawdir)
    
    
    
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

#######################

#
#if args.models==None:
#    print(args.task)
##    print(str(args.task))
#    if str(args.task)=="['flanker']":
#        models = supported_flanker
#        task = 'flanker'
#    elif str(args.tas)=="['recent_probes']":
#        models = supported_recent_probes
#        task = 'recent_probes'
#elif str(args.task)=="['flanker']":
#    task = 'flanker'
#elif str(args.task)=="['recent_probes']":
#    task = 'recent_probes'
#        
#print 'Processing ' + task + 'according to these parameterizations: ' + str(models)      
#rawdir = str(args.rawdir).strip('[]')
#print 'Input: ' + rawdir
#outputdir = str(args.outputdir).strip('[]')
#print 'Outputdir: ' + outputdir
#
#code = args.code
#print 'Response coding: ' + code
#
#nchains = args.nchains
#print 'Nchains: ' + str(nchains)
#
#nsamp = args.nsamp
#print 'Nsamples: ' + str(nsamp)
#
#nburn = args.nburn
#print 'Nburn: ' + str(nburn)
#
