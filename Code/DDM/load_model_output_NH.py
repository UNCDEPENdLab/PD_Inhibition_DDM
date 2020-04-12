#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 09:35:09 2017

@author: nth7
"""

#### load dependencies
import numpy as np
import hddm
import os
import pandas as pd
from pandas import Series
import matplotlib.pyplot as plt
import pymp


####################################################################################################################################
##Write a loop that goes through ics folder and reads output into a dataframe

SNAP_dim = ('DISINH', 'SUICPRON', 'NEGTEMP', 'MISTRUST', 'MANIP', 'AGG', 'SELFHARM', 'ECCPERC', 'DEPEN', 'POSTEMP', 'EXHIB', 'ENTITL', 'DETACH', 'IMPUL', 'PROPER', 'HDWK','DISINHP', 'LOSLFEST')

##
os.chdir('/Users/nth7/ics/PD_Inhibition_DDM/model_outputs_nth_full_2/')

#AGG_comb = hddm.load('AGG_comb_reg.model')

data_model_outputs = pd.DataFrame()

DICs = pd.DataFrame()

models = ('stim_v_reg', 'stim_a_reg', 'stim_comb_reg', 'v_reg', 'a_reg', 'comb_reg')

DICs = pd.DataFrame()
for pred in SNAP_dim:
    for mod in models:
        thism = hddm.load(pred + '_' + mod + '.model')
        DICs = DICs.append(pd.DataFrame({"pred": pred, "model":mod, "DIC":[thism.dic]}))
DICs.to_csv('DICs_flanker_SNE.csv')       

b_a = hddm.load('baseline_a_reg.model')
b_a.dic
b_v = hddm.load('baseline_v_reg.model')
b_v.dic
b_comb = hddm.load('baseline_comb_reg.model')
b_comb.dic

####################################################################################################################################
###adapt code below for above double for loop
#Pull traces and DICs from fitted models
nmodels=18
for m in SNAP_dim:
    posteriors = pd.DataFrame()
    ###models that simply vary V by between subjs factor
    v_model = hddm.load(m + '_v_reg.model')
    #store trace of the posterior of v in a model that predicts v based on between subjs factor
    posteriors['v_model'] = v_model.nodes_db.node['v_' + m].trace()
    #store trace of the posterior of condition effect in a model that predicts v based on between subjs factor
    posteriors['v_C_model'] = v_model.nodes_db.node["v_C(con, Treatment('Congruent'))[T.Incongruent]"].trace()   
    ###models that simply vary A by between subjs factor
    a_model = hddm.load(m + '_a_reg.model')
    posteriors['a_model'] = a_model.nodes_db.node['a_' + m].trace()
    posteriors['a_C_model'] = a_model.nodes_db.node["a_C(con, Treatment('Congruent'))[T.Incongruent]"].trace()    
    ##combined A + V models by between subjs factor
    comb_model = hddm.load(m + '_comb_reg.model')
    posteriors['comb_model_a'] = comb_model.nodes_db.node['a_' + m].trace()
    posteriors['comb_model_a_C'] = comb_model.nodes_db.node["a_C(con, Treatment('Congruent'))[T.Incongruent]"].trace()
    posteriors['comb_model_v'] = comb_model.nodes_db.node['v_' + m].trace()
    posteriors['comb_model_v_C'] = comb_model.nodes_db.node["v_C(con, Treatment('Congruent'))[T.Incongruent]"].trace()
   
    ##models that simply vary V by between subjs factor, condition, and their ixn
    v_stim_model = hddm.load(m + '_stim_v_reg.model')
    posteriors['v_stim_model_C'] = v_stim_model.nodes_db.node["v_C(con, Treatment('Congruent'))[T.Incongruent]"].trace()
    posteriors['v_stim_model_v'] = v_stim_model.nodes_db.node['v_' + m].trace()
    posteriors['v_stim_model_ixn'] = v_stim_model.nodes_db.node['v_' + m + ":C(con, Treatment('Congruent'))[T.Incongruent]"].trace()
    ##models that simply vary A by between subjs factor, condition, and their ixn
    a_stim_model = hddm.load(m + '_stim_a_reg.model')
    posteriors['a_stim_model_C'] = a_stim_model.nodes_db.node["a_C(con, Treatment('Congruent'))[T.Incongruent]"].trace()
    posteriors['a_stim_model_a'] = a_stim_model.nodes_db.node['a_' + m].trace()
    posteriors['a_stim_model_ixn'] = a_stim_model.nodes_db.node['a_' + m + ":C(con, Treatment('Congruent'))[T.Incongruent]"].trace()
    ##combined A + V models by between subjs factor, condition, and their ixn
    ###6/5/17: combined model does not include stimulus interactions (x_comb = x_stim_comb = [x_a, x_v] in SNAP_dim_flanker_mnh_norelpath.py on ICS)
    comb_stim_model = hddm.load(m + '_stim_comb_reg.model')
    #posteriors['a_stim_model_C_' + m] = a_stim_model.nodes_db.node["a_C(con, Treatment('Congruent'))[T.Incongruent]"].trace()
    ###code based on SNAP dimension
    posteriors['dimension'] = np.repeat(m, len(posteriors))
    ##append to previous SNAP dimensions data frame
    data_model_outputs = data_model_outputs.append(posteriors)

data_model_outputs.to_csv('flanker_model_outputs_SNE.csv')

