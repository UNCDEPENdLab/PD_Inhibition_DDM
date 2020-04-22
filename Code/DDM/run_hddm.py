,st,sv#!/usr/bin/env python2
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
supported_flanker = ['v', 'vst', 'v_block', 'v_blockst', 'vst', 'vstst', 'v_blockst', 'v_blockstst']
supported_recent_probes = ['v', 'vst', 'vst', 'vstst']


parser = argparse.ArgumentParser(description='Run Hierarchical Bayesian Drift Diffusion Model on experimental task data: http://ski.clps.brown.edu/hddm_docs/index.html\n\nSample call run on the command line: <python run_hddm.py flanker_accCode.cst /gpfs/group/mnh5174/default/Nate/HDDM_outputs_PD_Inhibition/samp1000/flanker/clean_sample flanker -m v -nc 10 -nb 5000 -ns 20000>')
parser.add_argument('rawdf', nargs=1, type=str, help = "path to raw dataframe (must be cst that match HDDM naming conventions)")
parser.add_argument('outputdir', nargs=1, type=str, help = "path to output directory")
parser.add_argument('task', nargs =1, help = "Experimental task to run HDDM on. The models to be run are coded into the function, so to get it to work for your needs you may need to mess with the HDDM  regressor calls within the function itself.", type = str)
parser.add_argument('-models', '-m', nargs='+', type =str, help = "N.B. CURRENTLY ONLY SUPPORTS ONE MODEL PARAMETERIZATION Different parameterizations of the DDM models to run. Currently supports " + str(supported_recent_probes) + " for recent probes and most combinations of: stim, stimblock (congruent vs. incongruent by block), trial number, trial in a specific run, previous rt, sv, and st for flanker. If this argument is left blank, will evaluate simple variants of these models, though it is suggested to specify one model, as multiple models will be run in parallel. More coming soon.")
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
data = hddm.load_cst(rawdf)

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
st = {'model': "st ~ 1 ", 'link_func': lambda x: x}
# z = {'model': "z ~ 1", 'link_func': z_link_func}
# sz = {'model': "sz ~ 1", 'link_func': lambda x: x}
st = {'model': "st ~ 1", 'link_func': lambda x: x}

if task == 'flanker':

    #### Simplest model: stim type only IV
    v = {'model': "v ~ 1  + C(stim, Treatment('congruent'))", 'link_func': lambda x: x}
    v_stimblock = = {'model': "v ~ 1  + C(stimblock, Treatment('congruent'))", 'link_func': lambda x: x}

    #### single IV (not including stim)
    # stimulus x block ixn
    v_block = {'model': "v ~ 1  + C(stim, Treatment('congruent')) * C(block, Treatment('most_con'))", 'link_func': lambda x: x}
    # standard stimulus contrasts
    v_trial = {'model': "v ~ 1  + C(stim, Treatment('congruent')+trial_z )", 'link_func': lambda x: x}
    v_runtrial = {'model': "v ~ 1  + C(stim, Treatment('congruent')+run_trial_z )", 'link_func': lambda x: x}
    v_prev_rt = {'model': "v ~ 1  + C(stim, Treatment('congruent')+ v_prev_rt)", 'link_func': lambda x: x}
    # stimblock contrasts
    v_stimblock_trial = {'model': "v ~ 1  + C(stimblock, Treatment('congruent')+trial_z )", 'link_func': lambda x: x}
    v_stimblock_runtrial = {'model': "v ~ 1  + C(stimblock, Treatment('congruent')+run_trial_z )", 'link_func': lambda x: x}
    v_stimblock_prev_rt = {'model': "v ~ 1  + C(stimblock, Treatment('congruent')+ v_prev_rt)", 'link_func': lambda x: x}

    #### two IVs (not including stim)
    # stimulus x block ixn
    v_block_trial = {'model': "v ~ 1  + C(stim, Treatment('congruent')) * C(block, Treatment('most_con')) + trial_z", 'link_func': lambda x: x}
    v_block_runtrial = {'model': "v ~ 1  + C(stim, Treatment('congruent')) * C(block, Treatment('most_con')) + run_trial_z", 'link_func': lambda x: x}
    v_block_prev_rt = {'model': "v ~ 1  + C(stim, Treatment('congruent')) * C(block, Treatment('most_con')) + prev_rt", 'link_func': lambda x: x}
    # standard stimulus contrasts
    v_trial_runtrial = {'model': "v ~ 1  + C(stim, Treatment('congruent')+trial_z +run_trial_z)", 'link_func': lambda x: x}
    v_trial_prev_rt = {'model': "v ~ 1  + C(stim, Treatment('congruent')+trial_z +prev_rt)", 'link_func': lambda x: x}
    v_runtrial_prev_rt = {'model': "v ~ 1  + C(stim, Treatment('congruent')+trial_z +run_trial_z)", 'link_func': lambda x: x}
    # stimblock contrasts
    v_stimblock_trial_runtrial = {'model': "v ~ 1  + C(stimblock, Treatment('congruent')+trial_z +run_trial_z)", 'link_func': lambda x: x}
    v_stimblock_trial_prev_rt = {'model': "v ~ 1  + C(stimblock, Treatment('congruent')+trial_z +prev_rt)", 'link_func': lambda x: x}
    v_stimblock_runtrial_prev_rt = {'model': "v ~ 1  + C(stimblock, Treatment('congruent')+trial_z +run_trial_z)", 'link_func': lambda x: x}

    #### three IVs (not including stim)
    # stimulus x block ixn
    v_block_trial_runtrial = {'model': "v ~ 1  + C(stim, Treatment('congruent')) * C(block, Treatment('most_con')) + trial_z + run_trial_z", 'link_func': lambda x: x}
    v_block_trial_prev_rt = {'model': "v ~ 1  + C(stim, Treatment('congruent')) * C(block, Treatment('most_con')) + trial_z + prev_rt", 'link_func': lambda x: x}
    v_block_runtrial_prev_rt = {'model': "v ~ 1  + C(stim, Treatment('congruent')) * C(block, Treatment('most_con')) + run_trial_z + prev_rt", 'link_func': lambda x: x}
    # stimblock contrasts
    v_stimblock_trial_runtrial = {'model': "v ~ 1  + C(stimblock, Treatment('congruent'))  + trial_z + run_trial_z", 'link_func': lambda x: x}
    v_stimblock_trial_prev_rt = {'model': "v ~ 1  + C(stimblock, Treatment('congruent'))  + trial_z + prev_rt", 'link_func': lambda x: x}
    v_stimblock_runtrial_prev_rt = {'model': "v ~ 1  + C(stimblock, Treatment('congruent')) + run_trial_z + prev_rt", 'link_func': lambda x: x}

    #### all four IVs (not including stim)
    # stimulus x block ixn
    v_block_trial_runtrial_prev_rt = {'model': "v ~ 1  + C(stim, Treatment('congruent')) * C(block, Treatment('most_con')) + trial_z + run_trial_z + prev_rt", 'link_func': lambda x: x}
    # stimblock contrasts
    v_stimblock_trial_runtrial_prev_rt = {'model': "v ~ 1  + C(stimblock, Treatment('congruent'))  + trial_z + run_trial_z + prev_rt", 'link_func': lambda x: x}

    ###################################
    ##### create dict to index that includes different combinations of the design matrices from above
    ###################################
    mod_dict = {
    ##################### start with no st or st
    #### Simplest model: stim type only IV
    'v':hddm.HDDMRegressor(data, v, group_only_regressors = False),
    'v_stimblock':hddm.HDDMRegressor(data, v_stimblock, group_only_regressors = False),
    #### single IV (not including stim)
    # stimulus x block ixn
    'v_block':hddm.HDDMRegressor(data, v_block, group_only_regressors = False),
    # standard stimulus contrasts
    'v_trial':hddm.HDDMRegressor(data, v_trial, group_only_regressors = False),
    'v_runtrial':hddm.HDDMRegressor(data, v_runtrial, group_only_regressors = False),
    'v_prev_rt':hddm.HDDMRegressor(data, v_prev_rt, group_only_regressors = False),
    # stimblock contrasts
    'v_stimblock_trial':hddm.HDDMRegressor(data, v_stimblock_trial, group_only_regressors = False),
    'v_stimblock_runtrial':hddm.HDDMRegressor(data, v_stimblock_runtrial, group_only_regressors = False),
    'v_stimblock_prev_rt':hddm.HDDMRegressor(data, v_stimblock_prev_rt, group_only_regressors = False),
    #### two IVs (not including stim)
    # stimulus x block ixn
    'v_block_trial':hddm.HDDMRegressor(data, v_block_trial, group_only_regressors = False),
    'v_block_runtrial':hddm.HDDMRegressor(data, v_block_runtrial, group_only_regressors = False),
    'v_block_prev_rt':hddm.HDDMRegressor(data, v_block_prev_rt, group_only_regressors = False),
    # standard stimulus contrasts
    'v_trial_runtrial':hddm.HDDMRegressor(data, v_trial_runtrial, group_only_regressors = False),
    'v_trial_prev_rt':hddm.HDDMRegressor(data, v_trial_prev_rt, group_only_regressors = False),
    'v_runtrial_prev_rt':hddm.HDDMRegressor(data, v_runtrial_prev_rt, group_only_regressors = False),
    # stimblock contrasts
    'v_stimblock_trial_runtrial':hddm.HDDMRegressor(data, v_stimblock_trial_runtrial, group_only_regressors = False),
    'v_stimblock_trial_prev_rt':hddm.HDDMRegressor(data, v_stimblock_trial_prev_rt, group_only_regressors = False),
    'v_stimblock_runtrial_prev_rt':hddm.HDDMRegressor(data, v_stimblock_runtrial_prev_rt, group_only_regressors = False),
    #### three IVs (not including stim)
    # stimulus x block ixn
    'v_block_trial_runtrial':hddm.HDDMRegressor(data, v_block_trial_runtrial, group_only_regressors = False),
    'v_block_trial_prev_rt':hddm.HDDMRegressor(data, v_block_trial_prev_rt, group_only_regressors = False),
    'v_block_runtrial_prev_rt':hddm.HDDMRegressor(data, v_block_runtrial_prev_rt, group_only_regressors = False),
    # stimblock contrasts
    'v_stimblock_trial_runtrial':hddm.HDDMRegressor(data, v_stimblock_trial_runtrial, group_only_regressors = False),
    'v_stimblock_trial_prev_rt':hddm.HDDMRegressor(data, v_stimblock_trial_prev_rt, group_only_regressors = False),
    'v_stimblock_runtrial_prev_rt':hddm.HDDMRegressor(data, v_stimblock_runtrial_prev_rt, group_only_regressors = False),
    #### all four IVs (not including stim)
    # stimulus x block ixn
    'v_block_trial_runtrial_prev_rt':hddm.HDDMRegressor(data, v_block_trial_runtrial_prev_rt, group_only_regressors = False),
    # stimblock contrasts
    'v_stimblock_trial_runtrial_prev_rt':hddm.HDDMRegressor(data, v_stimblock_trial_runtrial_prev_rt, group_only_regressors = False),

    ##################### add st
    #### Simplest model: stim type only IV
    'v_st':hddm.HDDMRegressor(data, [v,st], include = 'st', group_only_regressors = False),
    'v_stimblock_st':hddm.HDDMRegressor(data, [v_stimblock,st], include = 'st', group_only_regressors = False),
    #### single IV (not including stim)
    # stimulus x block ixn
    'v_block_st':hddm.HDDMRegressor(data, [v_block,st], include = 'st', group_only_regressors = False),
    # standard stimulus contrasts
    'v_trial_st':hddm.HDDMRegressor(data, [v_trial,st], include = 'st', group_only_regressors = False),
    'v_runtrial_st':hddm.HDDMRegressor(data, [v_runtrial,st], include = 'st', group_only_regressors = False),
    'v_prev_rt_st':hddm.HDDMRegressor(data, [v_prev_rt,st], include = 'st', group_only_regressors = False),
    # stimblock contrasts
    'v_stimblock_trial_st':hddm.HDDMRegressor(data, [v_stimblock_trial,st], include = 'st', group_only_regressors = False),
    'v_stimblock_runtrial_st':hddm.HDDMRegressor(data, [v_stimblock_runtrial,st], include = 'st', group_only_regressors = False),
    'v_stimblock_prev_rt_st':hddm.HDDMRegressor(data, [v_stimblock_prev_rt,st], include = 'st', group_only_regressors = False),
    #### two IVs (not including stim)
    # stimulus x block ixn
    'v_block_trial_st':hddm.HDDMRegressor(data, [v_block_trial,st], include = 'st', group_only_regressors = False),
    'v_block_runtrial_st':hddm.HDDMRegressor(data, [v_block_runtrial,st], include = 'st', group_only_regressors = False),
    'v_block_prev_rt_st':hddm.HDDMRegressor(data, [v_block_prev_rt,st], include = 'st', group_only_regressors = False),
    # standard stimulus contrasts
    'v_trial_runtrial_st':hddm.HDDMRegressor(data, [v_trial_runtrial,st], include = 'st', group_only_regressors = False),
    'v_trial_prev_rt_st':hddm.HDDMRegressor(data, [v_trial_prev_rt,st], include = 'st', group_only_regressors = False),
    'v_runtrial_prev_rt_st':hddm.HDDMRegressor(data, [v_runtrial_prev_rt,st], include = 'st', group_only_regressors = False),
    # stimblock contrasts
    'v_stimblock_trial_runtrial_st':hddm.HDDMRegressor(data, [v_stimblock_trial_runtrial,st], include = 'st', group_only_regressors = False),
    'v_stimblock_trial_prev_rt_st':hddm.HDDMRegressor(data, [v_stimblock_trial_prev_rt,st], include = 'st', group_only_regressors = False),
    'v_stimblock_runtrial_prev_rt_st':hddm.HDDMRegressor(data, [v_stimblock_runtrial_prev_rt,st], include = 'st', group_only_regressors = False),
    #### three IVs (not including stim)
    # stimulus x block ixn
    'v_block_trial_runtrial_st':hddm.HDDMRegressor(data, [v_block_trial_runtrial,st], include = 'st', group_only_regressors = False),
    'v_block_trial_prev_rt_st':hddm.HDDMRegressor(data, [v_block_trial_prev_rt,st], include = 'st', group_only_regressors = False),
    'v_block_runtrial_prev_rt_st':hddm.HDDMRegressor(data, [v_block_runtrial_prev_rt,st], include = 'st', group_only_regressors = False),
    # stimblock contrasts
    'v_stimblock_trial_runtrial_st':hddm.HDDMRegressor(data, [v_stimblock_trial_runtrial,st], include = 'st', group_only_regressors = False),
    'v_stimblock_trial_prev_rt_st':hddm.HDDMRegressor(data, [v_stimblock_trial_prev_rt,st], include = 'st', group_only_regressors = False),
    'v_stimblock_runtrial_prev_rt_st':hddm.HDDMRegressor(data, [v_stimblock_runtrial_prev_rt,st], include = 'st', group_only_regressors = False),
    #### all four IVs (not including stim)
    # stimulus x block ixn
    'v_block_trial_runtrial_prev_rt_st':hddm.HDDMRegressor(data, [v_block_trial_runtrial_prev_rt,st], include = 'st', group_only_regressors = False),
    # stimblock contrasts
    'v_stimblock_trial_runtrial_prev_rt_st':hddm.HDDMRegressor(data, [v_stimblock_trial_runtrial_prev_rt,st], include = 'st', group_only_regressors = False),


    ##################### add st
    #### Simplest model: stim type only IV
    'v_st':hddm.HDDMRegressor(data, [v,st], include = 'st', group_only_regressors = False),
    'v_stimblock_st':hddm.HDDMRegressor(data, [v_stimblock,st], include = 'st', group_only_regressors = False),
    #### single IV (not including stim)
    # stimulus x block ixn
    'v_block_st':hddm.HDDMRegressor(data, [v_block,st], include = 'st', group_only_regressors = False),
    # standard stimulus contrasts
    'v_trial_st':hddm.HDDMRegressor(data, [v_trial,st], include = 'st', group_only_regressors = False),
    'v_runtrial_st':hddm.HDDMRegressor(data, [v_runtrial,st], include = 'st', group_only_regressors = False),
    'v_prev_rt_st':hddm.HDDMRegressor(data, [v_prev_rt,st], include = 'st', group_only_regressors = False),
    # stimblock contrasts
    'v_stimblock_trial_st':hddm.HDDMRegressor(data, [v_stimblock_trial,st], include = 'st', group_only_regressors = False),
    'v_stimblock_runtrial_st':hddm.HDDMRegressor(data, [v_stimblock_runtrial,st], include = 'st', group_only_regressors = False),
    'v_stimblock_prev_rt_st':hddm.HDDMRegressor(data, [v_stimblock_prev_rt,st], include = 'st', group_only_regressors = False),
    #### two IVs (not including stim)
    # stimulus x block ixn
    'v_block_trial_st':hddm.HDDMRegressor(data, [v_block_trial,st], include = 'st', group_only_regressors = False),
    'v_block_runtrial_st':hddm.HDDMRegressor(data, [v_block_runtrial,st], include = 'st', group_only_regressors = False),
    'v_block_prev_rt_st':hddm.HDDMRegressor(data, [v_block_prev_rt,st], include = 'st', group_only_regressors = False),
    # standard stimulus contrasts
    'v_trial_runtrial_st':hddm.HDDMRegressor(data, [v_trial_runtrial,st], include = 'st', group_only_regressors = False),
    'v_trial_prev_rt_st':hddm.HDDMRegressor(data, [v_trial_prev_rt,st], include = 'st', group_only_regressors = False),
    'v_runtrial_prev_rt_st':hddm.HDDMRegressor(data, [v_runtrial_prev_rt,st], include = 'st', group_only_regressors = False),
    # stimblock contrasts
    'v_stimblock_trial_runtrial_st':hddm.HDDMRegressor(data, [v_stimblock_trial_runtrial,st], include = 'st', group_only_regressors = False),
    'v_stimblock_trial_prev_rt_st':hddm.HDDMRegressor(data, [v_stimblock_trial_prev_rt,st], include = 'st', group_only_regressors = False),
    'v_stimblock_runtrial_prev_rt_st':hddm.HDDMRegressor(data, [v_stimblock_runtrial_prev_rt,st], include = 'st', group_only_regressors = False),
    #### three IVs (not including stim)
    # stimulus x block ixn
    'v_block_trial_runtrial_st':hddm.HDDMRegressor(data, [v_block_trial_runtrial,st], include = 'st', group_only_regressors = False),
    'v_block_trial_prev_rt_st':hddm.HDDMRegressor(data, [v_block_trial_prev_rt,st], include = 'st', group_only_regressors = False),
    'v_block_runtrial_prev_rt_st':hddm.HDDMRegressor(data, [v_block_runtrial_prev_rt,st], include = 'st', group_only_regressors = False),
    # stimblock contrasts
    'v_stimblock_trial_runtrial_st':hddm.HDDMRegressor(data, [v_stimblock_trial_runtrial,st], include = 'st', group_only_regressors = False),
    'v_stimblock_trial_prev_rt_st':hddm.HDDMRegressor(data, [v_stimblock_trial_prev_rt,st], include = 'st', group_only_regressors = False),
    'v_stimblock_runtrial_prev_rt_st':hddm.HDDMRegressor(data, [v_stimblock_runtrial_prev_rt,st], include = 'st', group_only_regressors = False),
    #### all four IVs (not including stim)
    # stimulus x block ixn
    'v_block_trial_runtrial_prev_rt_st':hddm.HDDMRegressor(data, [v_block_trial_runtrial_prev_rt,st], include = 'st', group_only_regressors = False),
    # stimblock contrasts
    'v_stimblock_trial_runtrial_prev_rt_st':hddm.HDDMRegressor(data, [v_stimblock_trial_runtrial_prev_rt,st], include = 'st', group_only_regressors = False),


    ##################### add st and sv
    #### Simplest model: stim type only IV
    'v_sv_st':hddm.HDDMRegressor(data, [v,st,sv], include = ('st','sv'), group_only_regressors = False),
    'v_stimblock_sv_st':hddm.HDDMRegressor(data, [v_stimblock,st,sv], include = ('st','sv'), group_only_regressors = False),
    #### single IV (not including stim)
    # stimulus x block ixn
    'v_block_sv_st':hddm.HDDMRegressor(data, [v_block,st,sv], include = ('st','sv'), group_only_regressors = False),
    # standard stimulus contrasts
    'v_trial_sv_st':hddm.HDDMRegressor(data, [v_trial,st,sv], include = ('st','sv'), group_only_regressors = False),
    'v_runtrial_sv_st':hddm.HDDMRegressor(data, [v_runtrial,st,sv], include = ('st','sv'), group_only_regressors = False),
    'v_prev_rt_st':hddm.HDDMRegressor(data, [v_prev_rt,st,sv], include = ('st','sv'), group_only_regressors = False),
    # stimblock contrasts
    'v_stimblock_trial_sv_st':hddm.HDDMRegressor(data, [v_stimblock_trial,st,sv], include = ('st','sv'), group_only_regressors = False),
    'v_stimblock_runtrial_sv_st':hddm.HDDMRegressor(data, [v_stimblock_runtrial,st,sv], include = ('st','sv'), group_only_regressors = False),
    'v_stimblock_prev_rt_sv_st':hddm.HDDMRegressor(data, [v_stimblock_prev_rt,st,sv], include = ('st','sv'), group_only_regressors = False),
    #### two IVs (not including stim)
    # stimulus x block ixn
    'v_block_trial_sv_st':hddm.HDDMRegressor(data, [v_block_trial,st,sv], include = ('st','sv'), group_only_regressors = False),
    'v_block_runtrial_sv_st':hddm.HDDMRegressor(data, [v_block_runtrial,st,sv], include = ('st','sv'), group_only_regressors = False),
    'v_block_prev_rt_sv_st':hddm.HDDMRegressor(data, [v_block_prev_rt,st,sv], include = ('st','sv'), group_only_regressors = False),
    # standard stimulus contrasts
    'v_trial_runtrial_sv_st':hddm.HDDMRegressor(data, [v_trial_runtrial,st,sv], include = ('st','sv'), group_only_regressors = False),
    'v_trial_prev_rt_sv_st':hddm.HDDMRegressor(data, [v_trial_prev_rt,st,sv], include = ('st','sv'), group_only_regressors = False),
    'v_runtrial_prev_rt_sv_st':hddm.HDDMRegressor(data, [v_runtrial_prev_rt,st,sv], include = ('st','sv'), group_only_regressors = False),
    # stimblock contrasts
    'v_stimblock_trial_runtrial_sv_st':hddm.HDDMRegressor(data, [v_stimblock_trial_runtrial,st,sv], include = ('st','sv'), group_only_regressors = False),
    'v_stimblock_trial_prev_rt_sv_st':hddm.HDDMRegressor(data, [v_stimblock_trial_prev_rt,st,sv], include = ('st','sv'), group_only_regressors = False),
    'v_stimblock_runtrial_prev_rt_sv_st':hddm.HDDMRegressor(data, [v_stimblock_runtrial_prev_rt,st,sv], include = ('st','sv'), group_only_regressors = False),
    #### three IVs (not including stim)
    # stimulus x block ixn
    'v_block_trial_runtrial_sv_st':hddm.HDDMRegressor(data, [v_block_trial_runtrial,st,sv], include = ('st','sv'), group_only_regressors = False),
    'v_block_trial_prev_rt_sv_st':hddm.HDDMRegressor(data, [v_block_trial_prev_rt,st,sv], include = ('st','sv'), group_only_regressors = False),
    'v_block_runtrial_prev_rt_sv_st':hddm.HDDMRegressor(data, [v_block_runtrial_prev_rt,st,sv], include = ('st','sv'), group_only_regressors = False),
    # stimblock contrasts
    'v_stimblock_trial_runtrial_sv_st':hddm.HDDMRegressor(data, [v_stimblock_trial_runtrial,st,sv], include = ('st','sv'), group_only_regressors = False),
    'v_stimblock_trial_prev_rt_sv_st':hddm.HDDMRegressor(data, [v_stimblock_trial_prev_rt,st,sv], include = ('st','sv'), group_only_regressors = False),
    'v_stimblock_runtrial_prev_rt_sv_st':hddm.HDDMRegressor(data, [v_stimblock_runtrial_prev_rt,st,sv], include = ('st','sv'), group_only_regressors = False),
    #### all four IVs (not including stim)
    # stimulus x block ixn
    'v_block_trial_runtrial_prev_rt_sv_st':hddm.HDDMRegressor(data, [v_block_trial_runtrial_prev_rt,st,sv], include = ('st','sv'), group_only_regressors = False),
    # stimblock contrasts
    'v_stimblock_trial_runtrial_prev_rt_sv_st':hddm.HDDMRegressor(data, [v_stimblock_trial_runtrial_prev_rt,st,sv], include = ('st','sv'), group_only_regressors = False)
    }




elif task == 'recent_probes':

    v = {'model': "v ~ 1  + C(Condition, Treatment('positive'))", 'link_func': lambda x: x}

    ##### create dict to index that includes different combinations of the design matrices from above

    mod_dict = {'v':hddm.HDDMRegressor(data, v, group_only_regressors = False),
    'vst':hddm.HDDMRegressor(data, [v,st], include = 'st', group_only_regressors = False),
    'vst':hddm.HDDMRegressor(data, [v,st], include = 'st', group_only_regressors = False),
    'vstst':hddm.HDDMRegressor(data, [v,st,st], include = ('st','st'), group_only_regressors = False)}


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
