#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 11:36:50 2017

@author: nth7
"""
######################################################################
##PD_Inhibition run flanker models through HDDM 
######################################################################
##Read dependencies
import numpy as np
import hddm
import os
import pandas as pd
from pandas import Series
import matplotlib.pyplot as plt
import pymp
from patsy import dmatrix

##Create a string for SNAP (Schedule for Nonadaptive and Adaptive Personality) dimensions 
SNAP_dim = ('DISINH', 'SUICPRON', 'NEGTEMP', 'MISTRUST', 'MANIP', 'AGG', 'SELFHARM', 'ECCPERC', 'DEPEN', 'POSTEMP', 'EXHIB', 'ENTITL', 'DETACH', 'IMPUL', 'PROPER', 'HDWK','DISINHP', 'LOSLFEST')

ics = 0 #binary for working on local computer(0) or on ICS(1). Note that if working on local computer ICS must be mounted. 
code = 'stim' #coding scheme: stimulus(stim) or accuracy(acc) 
nmodels=6 #number of distinct models to run below
baseline_mods = 1 #estimate baseline models?
SNAP_mods = 0 #estimate dimensional models

######################################################################
##Read in data and set up burn and sample quantities
if ics == 0:
    #local
    os.chdir('/mnt/ics/DEPENd_Box/Projects/PD_Inhibition_DDM')
    if code == 'stim':
        data = hddm.load_csv('Flanker_reduced_stim_coded_nooutliers.csv')
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
    nsample = 30000
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


#####Estimate baseline models
if baseline_mods == 1:
    baseline_v = {'model': "v ~ 1  + C(condition)", 'link_func': lambda x: x}
    baseline_sv = {'model': "sv ~ 1  + C(condition)", 'link_func': lambda x: x}
    baseline_a = {'model': "a ~ 1", 'link_func': lambda x: x}
    baseline_aC = {'model': "a ~ 1  + C(condition)", 'link_func': lambda x: x}   #theoretical interpretation difficult, but test out baseline
    baseline_z = {'model': "z ~ 1 + C(condition)", 'link_func': z_link_func}

    
    baseline_av_comb = [baseline_a, baseline_v]
    baseline_asv_comb = [baseline_a, baseline_sv]
    baseline_aCv_comb = [baseline_aC, baseline_v]
    baseline_aCsv_comb = [baseline_aC, baseline_sv]
    baseline_avz_comb = [baseline_a, baseline_v, baseline_z]
    baseline_asvz_comb = [baseline_a, baseline_sv, baseline_z] 
    baseline_aCvz_comb = [baseline_aC, baseline_v, baseline_z]
    baseline_aCsvz_comb = [baseline_aC, baseline_sv, baseline_z]
    
    
    baseline_v_reg = hddm.HDDMRegressor(data, baseline_v)
    baseline_v_reg.sample(nsample, burn = nburn, dbname = 'baseline/baseline_v_reg.db', db = 'pickle')
    baseline_v_reg.save('baseline/baseline_v_reg.model')
    
    baseline_sv_reg = hddm.HDDMRegressor(data, baseline_sv, include = 'sv')
    baseline_sv_reg.sample(nsample, burn = nburn, dbname = 'baseline/baseline_sv_reg.db', db = 'pickle')
    baseline_sv_reg.save('baseline/baseline_sv_reg.model')
    
    baseline_a_reg = hddm.HDDMRegressor(data, baseline_a)
    baseline_a_reg.sample(nsample, burn = nburn, dbname = 'baseline/baseline_a_reg.db', db = 'pickle')
    baseline_a_reg.save('baseline/baseline_a_reg.model')
    
    baseline_aC_reg = hddm.HDDMRegressor(data, baseline_aC)
    baseline_aC_reg.sample(nsample, burn = nburn, dbname = 'baseline/baseline_aC_reg.db', db = 'pickle')
    baseline_aC_reg.save('baseline/baseline_aC_reg.model')
    
    baseline_av_comb_reg = hddm.HDDMRegressor(data, baseline_av_comb)
    baseline_av_comb_reg.sample(nsample, burn = nburn, dbname = 'baseline/baseline_av_comb_reg.db', db = 'pickle')
    baseline_av_comb_reg.save('baseline/baseline_av_comb_reg.model')
    
    baseline_asv_comb_reg = hddm.HDDMRegressor(data, baseline_asv_comb, include = 'sv')
    baseline_asv_comb_reg.sample(nsample, burn = nburn, dbname = 'baseline/baseline_asv_comb_reg.db', db = 'pickle')
    baseline_asv_comb_reg.save('baseline/baseline_asv_comb_reg.model')
   
    baseline_aCv_comb_reg = hddm.HDDMRegressor(data, baseline_aCv_comb)
    baseline_aCv_comb_reg.sample(nsample, burn = nburn, dbname = 'baseline/baseline_aCv_comb_reg.db', db = 'pickle')
    baseline_aCv_comb_reg.save('baseline/baseline_aCv_comb_reg.model')  
    
    baseline_aCsv_comb_reg = hddm.HDDMRegressor(data, baseline_aCsv_comb, include = 'sv')
    baseline_aCsv_comb_reg.sample(nsample, burn = nburn, dbname = 'baseline/baseline_aCsv_comb_reg.db', db = 'pickle')
    baseline_aCsv_comb_reg.save('baseline/baseline_aCsv_comb_reg.model')
    
    baseline_avz_comb_reg = hddm.HDDMRegressor(data, baseline_avz_comb, include = 'z')
    baseline_avz_comb_reg.sample(nsample, burn = nburn, dbname = 'baseline/baseline_avz_comb_reg.db', db = 'pickle')
    baseline_avz_comb_reg.save('baseline/baseline_avz_comb_reg.model')
    
    baseline_asvz_comb_reg = hddm.HDDMRegressor(data, baseline_asvz_comb, include = ('z', 'sv'))
    baseline_asvz_comb_reg.sample(nsample, burn = nburn, dbname = 'baseline/baseline_asvz_comb_reg.db', db = 'pickle')
    baseline_asvz_comb_reg.save('baseline/baseline_asvz_comb_reg.model')

    baseline_aCvz_comb_reg = hddm.HDDMRegressor(data, baseline_aCvz_comb, include = 'z')
    baseline_aCvz_comb_reg.sample(nsample, burn = nburn, dbname = 'baseline/baseline_aCvz_comb_reg.db', db = 'pickle')
    baseline_aCvz_comb_reg.save('baseline/baseline_aCvz_comb_reg.model')

    baseline_aCsvz_comb_reg = hddm.HDDMRegressor(data, baseline_aCsvz_comb, include = ('z', 'sv'))
    baseline_aCsvz_comb_reg.sample(nsample, burn = nburn, dbname = 'baseline/baseline_aCsvz_comb_reg.db', db = 'pickle')
    baseline_aCsvz_comb_reg.save('baseline/baseline_aCsvz_comb_reg.model')

##Estimate models on SNAP dimensions
if SNAP_mods == 1:
    if code == 'stim':
        with pymp.Parallel(len(SNAP_dim)) as p:
            with pymp.Parallel(nmodels) as pmodel:
                #for predict in p.iterate(SNAP_dim):    
                for index in p.range(0, len(SNAP_dim)):
                    for m in pmodel.range(0, nmodels):
                        predict = SNAP_dim[index]
                        x_v = {'model': "v ~ 1 + " + predict + " + C(condition)", 'link_func': lambda x: x}
                        x_v_ixn = {'model': "v ~ 1 + " + predict + " * C(condition)", 'link_func': lambda x: x}
                        x_sv = {'model': "sv ~ 1 + " + predict + " + C(condition)", 'link_func': lambda x: x}
                        x_sv_ixn = {'model': "sv ~ 1 + " + predict + " * C(condition)", 'link_func': lambda x: x}
                        x_a = {'model': "a ~ 1 + " + predict, 'link_func': lambda x: x}
                        x_z = {'model': "z ~ 1 + " + predict + " + C(condition)", 'link_func': z_link_func}        
                        
                        x_av_comb = [x_a, x_v]
                        x_av_ixn_comb = [x_a, x_v_ixn]
                        x_asv_comb = [x_a, x_sv]
                        x_asv_ixn_comb = [x_a, x_sv_ixn]
                        x_avz_comb = [x_a, x_v, x_z]
                        x_avz_ixn_comb = [x_a, x_v_ixn, x_z]
                        x_asvz_comb = [x_a, x_sv, x_z]
                        x_asvz_ixn_comb = [x_a, x_sv_ixn, x_z]
                        
                        
                        if m == 0:
                            ##DRIFT RATE (V)
                            x_stim_v_reg = hddm.HDDMRegressor(data, x_stim_v)
                            x_stim_v_reg.sample(nsample, burn = nburn, dbname = predict + '_stim_v_reg.db', db = 'pickle')
                            x_stim_v_reg.save(predict + '_stim_v_reg.model')
                        elif m == 1:    
                            ###DECISION THRESHOLD (A)
                            x_stim_a_reg = hddm.HDDMRegressor(data, x_stim_a)
                            x_stim_a_reg.sample(nsample, burn = nburn, dbname = predict + '_stim_a_reg.db', db = 'pickle')
                            x_stim_a_reg.save(predict + '_stim_a_reg.model')
                        elif m==2:
                        #if m==2:
                            #####COMBINED A and V model
                            x_stim_comb_reg = hddm.HDDMRegressor(data, x_stim_comb)
                            x_stim_comb_reg.sample(nsample, burn = nburn, dbname = predict + '_stim_comb_reg.db', db = 'pickle')
                            x_stim_comb_reg.save(predict + '_stim_comb_reg.model')
                        elif m == 3:
                            ##DRIFT RATE (V)
                            x_v_reg = hddm.HDDMRegressor(data, x_v)
                            x_v_reg.sample(nsample, burn = nburn, dbname = predict + '_v_reg.db', db = 'pickle')
                            x_v_reg.save(predict + '_v_reg.model')
                        elif m == 4:    
                            ###DECISION THRESHOLD (A)
                            x_a_reg = hddm.HDDMRegressor(data, x_a)
                            x_a_reg.sample(nsample, burn = nburn, dbname = predict + '_a_reg.db', db = 'pickle')
                            x_a_reg.save(predict + '_a_reg.model')
                        elif m==5:
                            #####COMBINED A and V model
                            x_comb_reg = hddm.HDDMRegressor(data, x_comb)
                            x_comb_reg.sample(nsample, burn = nburn, dbname = predict + '_comb_reg.db', db = 'pickle')
                            x_comb_reg.save(predict + '_comb_reg.model')