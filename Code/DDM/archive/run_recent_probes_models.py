#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 11:56:44 2020

@author: natehall
"""
######################################################################
## PD_Inhibition run recent probes models through HDDM 
######################################################################
##Read dependencies - some borrowed from Sophie
import hddm
import os
import pandas as pd
from pandas import Series
from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.pyplot as plt 
import pymp
from patsy import dmatrix
import kabuki
from kabuki.analyze import gelman_rubin
import pickle
import numpy as np
import sys

#######################
## setup fitting parameters
ics = 1 #binary for working on local computer(0) or on ICS(1). Note that if working on local computer ICS must be mounted. 
code = 'acc' #'stim' #coding scheme: stimulus(stim) or accuracy(acc) 
nchains = 5
short_chains = 1 # binary for ICS fitting models with fewer samples in each chain. This option is meant to get rough posteriors more quickly, though longer chains may be required eventually

# vestigial from earlier days

##Create a string for SNAP (Schedule for Nonadaptive and Adaptive Personality) dimensions 
#SNAP_dim = ('DISINH', 'SUICPRON', 'NEGTEMP', 'MISTRUST', 'MANIP', 'AGG', 'SELFHARM', 'ECCPERC', 'DEPEN', 'POSTEMP', 'EXHIB', 'ENTITL', 'DETACH', 'IMPUL', 'PROPER', 'HDWK','DISINHP', 'LOSLFEST')

#nmodels=6 #number of distinct models to run below
#mods = 1 #estimate baseline models?
#SNAP_mods = 0 #estimate dimensional models

######################################################################
##Read in data and set up burn and sample quantities
if ics == 0:
    #local
#    basedir = '/Users/natehall/github_repos/PD_Inhibition_DDM'
    basedir = '/Users/natehall/ics/Nate/PD_Inhibition_DDM'   # mounted ACI directory
    os.chdir(basedir)
    if code == 'stim':
        data = hddm.load_csv('Data/preprocessed/rp_stimCode.csv')
    elif code == 'acc':
        data = hddm.load_csv('Data/preprocessed/rp_accCode.csv')
        data = hddm.utils.flip_errors(data)
    outputdir = basedir+'/Outputs/practice_hddm/recent_probes'
    nsample = 10
    nburn = 3
elif ics == 1:
    basedir = '/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM'
    os.chdir(basedir)
    if code == 'stim':
        data = hddm.load_csv('Data/preprocessed/rp_stimCode.csv')
    elif code == 'acc':
        data = hddm.load_csv('Data/preprocessed/rp_accCode.csv')
        data = hddm.utils.flip_errors(data)
    if short_chains:
        outputdir = basedir+'/Outputs/practice_hddm/recent_probes'
        nsample = 2000
        nburn = 100
    else:
        outputdir = basedir+'/Outputs/full_hddm/recent_probes'
        nsample = 20000
        nburn = 5000

os.chdir(outputdir)    


## examine RT distributions across subjects
fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in data.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)

plt.savefig(basedir + '/Figures/rp_rt_dists_signflipped.pdf')


##configure pymp for parallel processing
pymp.config.nested=True

##define z link function to fix values of z to be 0-1 if needed.
def z_link_func(x, data=data):
    condition = (np.asarray(dmatrix('0 + C(s, [[1], [-1]])',
                               {'s': data.condition.ix[x.index]}))
    )
    return 1 / (1 + np.exp(-(x * condition)))

##### generate design matrices for HDDMRegressor
v = {'model': "v ~ 1  + C(Condition, Treatment('positive'))", 'link_func': lambda x: x}
sv = {'model': "sv ~ 1 ", 'link_func': lambda x: x}
st = {'model': "st ~ 1", 'link_func': lambda x: x}


##### create dict to index that includes different combinations of the design matrices from above
mod_dict = {'v_reg':hddm.HDDMRegressor(data, v, group_only_regressors = False),
'vsv_reg':hddm.HDDMRegressor(data, [v,sv], include = 'sv', group_only_regressors = False),
'vst_reg':hddm.HDDMRegressor(data, [v,st], include = 'st', group_only_regressors = False),
'vsvst_reg':hddm.HDDMRegressor(data, [v,sv,st], include = ('sv','st'), group_only_regressors = False)}
# 'vz_reg':hddm.HDDMRegressor(data, [v,z], include = 'z')}
#'vsvz_reg':hddm.HDDMRegressor(data, [v,sv,z], include = ('sv','z')),
#'v_blockz_reg':hddm.HDDMRegressor(data, [v_block,z], include = 'z'),
#'v_blocksvz_reg':hddm.HDDMRegressor(data, [v_block,sv,z], include = ('sv', 'z'))}
                                        
##### parallel loop over models and number of chains for gelman rubin statistic. Sample and save 
with pymp.Parallel(len(mod_dict)) as p:
    with pymp.Parallel(nchains) as ch:
        for index in p.range(0, len(mod_dict)):
            for ch_index in ch.range(0,nchains):
                model = mod_dict.values()[index]
                model.sample(nsample, burn = nburn, dbname = mod_dict.keys()[index] + '_rp_chain'+str(ch_index)+'.db', db = 'pickle')
                model.save(mod_dict.keys()[index] + '_rp_chain'+str(ch_index)+'.model')


