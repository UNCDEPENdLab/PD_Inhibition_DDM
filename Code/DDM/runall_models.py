#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 11:36:50 2017

@author: nth7
"""
######################################################################
##PD_Inhibition run flanker models through HDDM 
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

##Create a string for SNAP (Schedule for Nonadaptive and Adaptive Personality) dimensions 
#SNAP_dim = ('DISINH', 'SUICPRON', 'NEGTEMP', 'MISTRUST', 'MANIP', 'AGG', 'SELFHARM', 'ECCPERC', 'DEPEN', 'POSTEMP', 'EXHIB', 'ENTITL', 'DETACH', 'IMPUL', 'PROPER', 'HDWK','DISINHP', 'LOSLFEST')

ics = 0 #binary for working on local computer(0) or on ICS(1). Note that if working on local computer ICS must be mounted. 
code = 'stim' #coding scheme: stimulus(stim) or accuracy(acc) 
nchains = 5

# vestigial from earlier days
#nmodels=6 #number of distinct models to run below
#mods = 1 #estimate baseline models?
#SNAP_mods = 0 #estimate dimensional models

######################################################################
##Read in data and set up burn and sample quantities
if ics == 0:
    #local
    os.chdir('/Users/nth7/PD_Inhibition_DDM')
    if code == 'stim':
        data = hddm.load_csv('Data/preprocessed/flank_use.csv')
    elif code == 'acc':
        data = hddm.load_csv('Flanker_reduced_acc_coded_nooutliers.csv')
    outputdir = ('practice_files') 
    nsample = 50
    nburn = 3
    os.chdir('/mnt/ics/PD_Inhibition_DDM_bashfiles_outputs/practice_files')
elif ics == 1:
    os.chdir('/gpfs/group/mnh5174/default/DEPENd_Box/Projects/PD_Inhibition_DDM')
    if code == 'stim':
        data = hddm.load_csv('Flanker_reduced_stim_coded_nooutliers.csv')
    elif code == 'acc':
        data = hddm.load_csv('Flanker_reduced_acc_coded_nooutliers.csv')
    #outputdir = ('practice_files') 
    nsample = 20000
    nburn = 5000
    os.chdir('/mnt/ics/PD_Inhibition_DDM_bashfiles_outputs/full_models')

##configure pymp for parallel processing
pymp.config.nested=True

##define z link function to fix values of z to be 0-1
def z_link_func(x, data=data):
    condition = (np.asarray(dmatrix('0 + C(s, [[1], [-1]])',
                               {'s': data.condition.ix[x.index]}))
    )
    return 1 / (1 + np.exp(-(x * condition)))


##### Estimate models: step 1 simple base models
v = {'model': "v ~ 1  + C(stim, Treatment(0))", 'link_func': lambda x: x}
#v = {'model': "v ~ 1  ", 'link_func': lambda x: x}
v_block = {'model': "v ~ 1  + C(stim, Treatment(0)) * C(CongruentBlock, Treatment(0))", 'link_func': lambda x: x}
sv = {'model': "sv ~ 1 ", 'link_func': lambda x: x}
z = {'model': "z ~ 1", 'link_func': z_link_func}
sz = {'model': "sz ~ 1", 'link_func': lambda x: x}
a = {'model': "a ~ 1", 'link_func': lambda x: x}
t = {'model': "t ~ 1", 'link_func': lambda x: x}
st = {'model': "st ~ 1", 'link_func': lambda x: x}


mod_dict = {'v_reg':hddm.HDDMRegressor(data, v),
'vsv_reg':hddm.HDDMRegressor(data, [v,sv], include = 'sv'),
'vz_reg':hddm.HDDMRegressor(data, [v,z], include = 'z')}#,
#'vsvz_reg':hddm.HDDMRegressor(data, [v,sv,z], include = ('sv','z')),
#'v_block_reg':hddm.HDDMRegressor(data, v_block),
#'v_blocksv_reg':hddm.HDDMRegressor(data, [v_block,sv], include = 'sv'),
#'v_blockz_reg':hddm.HDDMRegressor(data, [v_block,z], include = 'z'),
#'v_blocksvz_reg':hddm.HDDMRegressor(data, [v_block,sv,z], include = ('sv', 'z'))}
                                        
mod_dict
all_models = {}

with pymp.Parallel(len(mod_dict)) as m:
    with pymp.Parallel(nchains) as ch:
            for mod in m.range(0, len(mod_dict)):
                print mod
                models = []
                for chain in ch.range(0, nchains):
                    mod_dict[mod].sample(nsample, burn = nburn, dbname = 'flanker/v_reg.db', db = 'pickle')    
                    models.append(mod_dict[mod])
                all_models[mod] = models
      
models = []
with pymp.Parallel(len(mod_dict)) as m:
    for mod in mod_dict:
        models.append(mod)
        
        
        
for mod in mod_dict:
    print mod
    models = []
    for chain in range(nchains):
        mod_dict[mod].sample(nsample, burn = nburn, dbname = 'flanker/v_reg.db', db = 'pickle')    
        models.append(mod_dict[mod])
    all_models[mod] = models


all_models.save('flanker/all_models_simple_multchains')



with pymp.Parallel(len(SNAP_dim)) as p:
            with pymp.Parallel(nmodels) as pmodel:
                #for predict in p.iterate(SNAP_dim):    
                for index in p.range(0, len(SNAP_dim)):
                    for m in pmodel.range(0, nmodels):




stats = mod_dict[mod].gen_stats()











#sv_reg = hddm.HDDMRegressor(data, sv, include = 'sv')
#sv_reg.sample(nsample, burn = nburn, dbname = 'baseline/sv_reg.db', db = 'pickle')
#sv_reg.save('baseline/sv_reg.model')
#
#a_reg = hddm.HDDMRegressor(data, a)
#a_reg.sample(nsample, burn = nburn, dbname = 'baseline/a_reg.db', db = 'pickle')
#a_reg.save('baseline/a_reg.model')
#
#aC_reg = hddm.HDDMRegressor(data, aC)
#aC_reg.sample(nsample, burn = nburn, dbname = 'baseline/aC_reg.db', db = 'pickle')
#aC_reg.save('baseline/aC_reg.model')
#
#av_comb_reg = hddm.HDDMRegressor(data, av_comb)
#av_comb_reg.sample(nsample, burn = nburn, dbname = 'baseline/av_comb_reg.db', db = 'pickle')
#av_comb_reg.save('baseline/av_comb_reg.model')
#
#asv_comb_reg = hddm.HDDMRegressor(data, asv_comb, include = 'sv')
#asv_comb_reg.sample(nsample, burn = nburn, dbname = 'baseline/asv_comb_reg.db', db = 'pickle')
#asv_comb_reg.save('baseline/asv_comb_reg.model')
#   
#aCv_comb_reg = hddm.HDDMRegressor(data, aCv_comb)
#aCv_comb_reg.sample(nsample, burn = nburn, dbname = 'baseline/aCv_comb_reg.db', db = 'pickle')
#aCv_comb_reg.save('baseline/aCv_comb_reg.model')  
#
#aCsv_comb_reg = hddm.HDDMRegressor(data, aCsv_comb, include = 'sv')
#aCsv_comb_reg.sample(nsample, burn = nburn, dbname = 'baseline/aCsv_comb_reg.db', db = 'pickle')
#aCsv_comb_reg.save('baseline/aCsv_comb_reg.model')
#
#avz_comb_reg = hddm.HDDMRegressor(data, avz_comb, include = 'z')
#avz_comb_reg.sample(nsample, burn = nburn, dbname = 'baseline/avz_comb_reg.db', db = 'pickle')
#avz_comb_reg.save('baseline/avz_comb_reg.model')
#
#asvz_comb_reg = hddm.HDDMRegressor(data, asvz_comb, include = ('z', 'sv'))
#asvz_comb_reg.sample(nsample, burn = nburn, dbname = 'baseline/asvz_comb_reg.db', db = 'pickle')
#asvz_comb_reg.save('baseline/asvz_comb_reg.model')
#
#aCvz_comb_reg = hddm.HDDMRegressor(data, aCvz_comb, include = 'z')
#aCvz_comb_reg.sample(nsample, burn = nburn, dbname = 'baseline/aCvz_comb_reg.db', db = 'pickle')
#aCvz_comb_reg.save('baseline/aCvz_comb_reg.model')
#
#aCsvz_comb_reg = hddm.HDDMRegressor(data, aCsvz_comb, include = ('z', 'sv'))
#aCsvz_comb_reg.sample(nsample, burn = nburn, dbname = 'baseline/aCsvz_comb_reg.db', db = 'pickle')
#aCsvz_comb_reg.save('baseline/aCsvz_comb_reg.model')

# N.B. 3/3/2020 incorporating SNAP dims into these models has been deprecated. We will look at these differences on the post-estimation side.

##Estimate models on SNAP dimensions
#if SNAP_mods == 1:
#    if code == 'stim':
#        with pymp.Parallel(len(SNAP_dim)) as p:
#            with pymp.Parallel(nmodels) as pmodel:
#                #for predict in p.iterate(SNAP_dim):    
#                for index in p.range(0, len(SNAP_dim)):
#                    for m in pmodel.range(0, nmodels):
#                        predict = SNAP_dim[index]
#                        x_v = {'model': "v ~ 1 + " + predict + " + C(condition)", 'link_func': lambda x: x}
#                        x_v_ixn = {'model': "v ~ 1 + " + predict + " * C(condition)", 'link_func': lambda x: x}
#                        x_sv = {'model': "sv ~ 1 + " + predict + " + C(condition)", 'link_func': lambda x: x}
#                        x_sv_ixn = {'model': "sv ~ 1 + " + predict + " * C(condition)", 'link_func': lambda x: x}
#                        x_a = {'model': "a ~ 1 + " + predict, 'link_func': lambda x: x}
#                        x_z = {'model': "z ~ 1 + " + predict + " + C(condition)", 'link_func': z_link_func}        
#                        
#                        x_av_comb = [x_a, x_v]
#                        x_av_ixn_comb = [x_a, x_v_ixn]
#                        x_asv_comb = [x_a, x_sv]
#                        x_asv_ixn_comb = [x_a, x_sv_ixn]
#                        x_avz_comb = [x_a, x_v, x_z]
#                        x_avz_ixn_comb = [x_a, x_v_ixn, x_z]
#                        x_asvz_comb = [x_a, x_sv, x_z]
#                        x_asvz_ixn_comb = [x_a, x_sv_ixn, x_z]
#                        
#                        
#                        if m == 0:
#                            ##DRIFT RATE (V)
#                            x_stim_v_reg = hddm.HDDMRegressor(data, x_stim_v)
#                            x_stim_v_reg.sample(nsample, burn = nburn, dbname = predict + '_stim_v_reg.db', db = 'pickle')
#                            x_stim_v_reg.save(predict + '_stim_v_reg.model')
#                        elif m == 1:    
#                            ###DECISION THRESHOLD (A)
#                            x_stim_a_reg = hddm.HDDMRegressor(data, x_stim_a)
#                            x_stim_a_reg.sample(nsample, burn = nburn, dbname = predict + '_stim_a_reg.db', db = 'pickle')
#                            x_stim_a_reg.save(predict + '_stim_a_reg.model')
#                        elif m==2:
#                        #if m==2:
#                            #####COMBINED A and V model
#                            x_stim_comb_reg = hddm.HDDMRegressor(data, x_stim_comb)
#                            x_stim_comb_reg.sample(nsample, burn = nburn, dbname = predict + '_stim_comb_reg.db', db = 'pickle')
#                            x_stim_comb_reg.save(predict + '_stim_comb_reg.model')
#                        elif m == 3:
#                            ##DRIFT RATE (V)
#                            x_v_reg = hddm.HDDMRegressor(data, x_v)
#                            x_v_reg.sample(nsample, burn = nburn, dbname = predict + '_v_reg.db', db = 'pickle')
#                            x_v_reg.save(predict + '_v_reg.model')
#                        elif m == 4:    
#                            ###DECISION THRESHOLD (A)
#                            x_a_reg = hddm.HDDMRegressor(data, x_a)
#                            x_a_reg.sample(nsample, burn = nburn, dbname = predict + '_a_reg.db', db = 'pickle')
#                            x_a_reg.save(predict + '_a_reg.model')
#                        elif m==5:
#                            #####COMBINED A and V model
#                            x_comb_reg = hddm.HDDMRegressor(data, x_comb)
#                            x_comb_reg.sample(nsample, burn = nburn, dbname = predict + '_comb_reg.db', db = 'pickle')
#                            x_comb_reg.save(predict + '_comb_reg.model')