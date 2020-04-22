#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 12:45:27 2020

@author: natehall
"""

######################################################################
## PD_Inhibition run models through HDDM based on inputs to a bash script (typically run as a PBS job)
######################################################################
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

#######################



##############################################
## setup parser and parse inputs
##############################################

## in case the user inputs nothing for the models they would like to have tested. These are the defaults
supported_flanker = ['v', 'vsv', 'v_block', 'v_blocksv', 'vst', 'vsvst', 'v_blockst', 'v_blocksvst']
supported_recent_probes = ['v', 'vsv', 'vst', 'vsvst']
 

parser = argparse.ArgumentParser(description='Run Hierarchical Bayesian Drift Diffusion Model on experimental task data: http://ski.clps.brown.edu/hddm_docs/index.html\n\nSample call run on the command line: <python run_hddm.py flanker_accCode.csv /gpfs/group/mnh5174/default/Nate/HDDM_outputs_PD_Inhibition/samp1000/flanker/clean_sample flanker -m v -nc 10 -nb 5000 -ns 20000>')
parser.add_argument('rawdf', nargs=1, type=str, help = "path to raw dataframe (must be csv that match HDDM naming conventions)")
parser.add_argument('outputdir', nargs=1, type=str, help = "path to output directory")
parser.add_argument('task', nargs =1, help = "Experimental task to run HDDM on. The models to be run are coded into the function, so to get it to work for your needs you may need to mess with the HDDM  regressor calls within the function itself.", type = str)
#if par 
parser.add_argument('-models', '-m', nargs='+', type =str, help = "N.B. CURRENTLY ONLY SUPPORTS ONE MODEL PARAMETERIZATION Different parameterizations of the DDM models to run. Currently supports " + str(supported_recent_probes) + " for recent probes and " + str(supported_flanker) + " for flanker. If this argument is left blank, will evaluate all of these. More coming soon.")
parser.add_argument('-code', '-c', help="coding scheme: stimulus(stim) or accuracy(acc)", type= str, default = "acc")
parser.add_argument('-nchains', '-nc', help="Number of unique MCMC chains to run. N.B. these are run in parallel so be careful the requested processors matches. Defaults to 1", type= int, default= 1)
parser.add_argument('-nsamp', '-ns', help = "total number of draws from MCMC posterior. Defaults to 2000", type = int, default = 2000)
parser.add_argument('-nburn', '-nb', help = "number of burn in samples. Defaults to 500", type = int, default = 500)
#parser.add_argument('-verbose', '-v', help = "Logical. whether or not to narrate steps of processing. Defaults to true")

print '\n################################################ \n################################################ \nParsing arguments: \n'

args = parser.parse_args()

#if args.verbose==None:
#    verb = True
#else: 
#    verb = False

print(args)

# reassign parsed arguments to variables outside of args

if args.models==None:
    print(args.task)
    if str(args.task)=="['flanker']":
        models = supported_flanker
    elif str(args.tas)=="['recent_probes']":
        models = supported_recent_probes
else:
    models = str(args.models).strip('[]')
#    models = args.models
    

models = models.replace("'","")
models = [models]


task = str(args.task).strip('[]')  
task = task.replace("'","") 
print 'Processing ' + task + ' according to these parameterizations: ' + str(models)      
rawdf = str(args.rawdf).strip('[]')
rawdf = rawdf.replace("'", "") # this took way too long to figure out.

print 'Input: ' + rawdf
outputdir = str(args.outputdir).strip('[]')
outputdir = outputdir.replace("'", "")
print 'Outputdir: ' + outputdir

code = args.code
print 'Response coding: ' + code

nchains = args.nchains
print 'Nchains: ' + str(nchains)

nsample = args.nsamp
print 'Nsamples: ' + str(nsample)

nburn = args.nburn
print 'Nburn: ' + str(nburn) + '\n################################################\n################################################\n'

#if args.verbose==None:
#    verb = True
#else: 
#    verb = False

########################
 


##############################################
## load data and flip incorrect RTs if accuracy coding
##############################################
print 'loading raw RT data from: ' + rawdf
data = hddm.load_csv(rawdf)

print 'data structure: '
print data.head(5)



if code == 'acc':
    print '\nflipping RT errors\n'
    data = hddm.utils.flip_errors(data)

## navigate to outputdir
os.chdir(outputdir) 
print 'Navigated to ' + os.getcwd() 

print 'Directory contents:\n' 
print os.listdir('.')

###define z link function to fix values of z to be 0-1 if needed.
def z_link_func(x, data=data):
    condition = (np.asarray(dmatrix('0 + C(s, [[1], [-1]])',
                               {'s': data.condition.ix[x.index]}))
    )
    return 1 / (1 + np.exp(-(x * condition)))
   
#######################



##############################################
## generate design matrices for HDDMRegressor 
##############################################
##### 

print '\nGenerating design matrices for HDDMRegressor'

# define regressions that can be estimated for all tasks
sv = {'model': "sv ~ 1 ", 'link_func': lambda x: x}
# z = {'model': "z ~ 1", 'link_func': z_link_func}
# sz = {'model': "sz ~ 1", 'link_func': lambda x: x}
st = {'model': "st ~ 1", 'link_func': lambda x: x}

if task == 'flanker':
    
    v = {'model': "v ~ 1  + C(stim, Treatment(0))", 'link_func': lambda x: x}
    v_block = {'model': "v ~ 1  + C(stim, Treatment(0)) * C(block, Treatment(0))", 'link_func': lambda x: x}
    
    ##### create dict to index that includes different combinations of the design matrices from above
    
    mod_dict = {'v':hddm.HDDMRegressor(data, v, group_only_regressors = False),
    'vsv':hddm.HDDMRegressor(data, [v,sv], include = 'sv', group_only_regressors = False),
    'v_block':hddm.HDDMRegressor(data, v_block, group_only_regressors = False),
    'v_blocksv':hddm.HDDMRegressor(data, [v_block,sv], include = 'sv',group_only_regressors = False),
    'vst':hddm.HDDMRegressor(data, [v,st], include = 'st', group_only_regressors = False),
    'vsvst':hddm.HDDMRegressor(data, [v,sv,st], include = ('sv','st'), group_only_regressors = False),
    'v_blockst':hddm.HDDMRegressor(data, [v_block, st], include = 'st', group_only_regressors = False),
    'v_blocksvst':hddm.HDDMRegressor(data, [v_block,sv,st], include = ('sv', 'st'),group_only_regressors = False)}
    
    # 'vz_reg':hddm.HDDMRegressor(data, [v,z], include = 'z')}
    #'vsvz_reg':hddm.HDDMRegressor(data, [v,sv,z], include = ('sv','z')),
    #'v_blockz_reg':hddm.HDDMRegressor(data, [v_block,z], include = 'z'),
    #'v_blocksvz_reg':hddm.HDDMRegressor(data, [v_block,sv,z], include = ('sv', 'z'))}
    # mod_dict = {'v_blocksvst_reg':hddm.HDDMRegressor(data, [v_block,sv,st], include = ('sv', 'st'),group_only_regressors = False)}
elif task == 'recent_probes':
    
    v = {'model': "v ~ 1  + C(Condition, Treatment('positive'))", 'link_func': lambda x: x}
    
    ##### create dict to index that includes different combinations of the design matrices from above
    
    mod_dict = {'v':hddm.HDDMRegressor(data, v, group_only_regressors = False),
    'vsv':hddm.HDDMRegressor(data, [v,sv], include = 'sv', group_only_regressors = False),
    'vst':hddm.HDDMRegressor(data, [v,st], include = 'st', group_only_regressors = False),
    'vsvst':hddm.HDDMRegressor(data, [v,sv,st], include = ('sv','st'), group_only_regressors = False)}
    # 'vz_reg':hddm.HDDMRegressor(data, [v,z], include = 'z')}
    #'vsvz_reg':hddm.HDDMRegressor(data, [v,sv,z], include = ('sv','z')),
    #'v_blockz_reg':hddm.HDDMRegressor(data, [v_block,z], include = 'z'),
    #'v_blocksvz_reg':hddm.HDDMRegressor(data, [v_block,sv,z], include = ('sv', 'z'))}

#######################

# allows for DDM to drop certain parameterizations if not requested.
mod_dict_torun = {key: mod_dict[key] for key in models}

print 'final models getting passed to HDDM: ' 
print mod_dict_torun 
##############################################
## parallel loop over models and number of chains. Sample and save 
##############################################
##### 

###configure pymp for parallel processing
pymp.config.nested=True

#N.B. 4/9/20 per MH's suggestion, it is advised that this script be called to run one model and parallelize over the number of chains requested.
## however, the function will still support running multiple models, though this will now be executed serially.

print '\n' + 'COMMENCE SAMPLING\n'




#with pymp.Parallel(len(mod_dict_torun)) as p:
for index in range(0, len(mod_dict_torun)):
#    for ch_index in range(0,nchains):
    with pymp.Parallel(nchains) as ch:     
        for ch_index in ch.range(0,nchains):
            model = mod_dict_torun.values()[index]
            print model
            dbname = mod_dict_torun.keys()[index] + '_chain'+str(ch_index)+ '_'+ code +'Code.db'
            print dbname
            modelname = mod_dict_torun.keys()[index] + '_chain'+str(ch_index)+ '_'+ code +'Code.model'
            print modelname
            
#            nsample = 10
#            nburn = 1
            
            model.sample(nsample, burn = nburn, dbname = dbname, db = 'pickle')
            model.save(modelname)


