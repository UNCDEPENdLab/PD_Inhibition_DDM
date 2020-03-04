#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 17:08:49 2017

@author: nth7
"""

import numpy as np
import hddm
import os
import pandas as pd
from pandas import Series
import matplotlib.pyplot as plt
import pymp

os.chdir('/Users/nth7/ics/PD_Inhibition_DDM')
data = hddm.load_csv('Flanker_reduced_nooutliers_mnh.csv')
#outputdir = ('model_outputs_flanker_nth_full')



nsample = 5000
nburn = 500

os.chdir('/Users/nth7/ics/PD_Inhibition_DDM/model_outputs_nth_full_2')

baseline_v = {'model': "v ~ 1  + C(con, Treatment('Congruent'))", 'link_func': lambda x: x}
baseline_a = {'model': "a ~ 1 + C(con, Treatment('Congruent'))", 'link_func': lambda x: x}
baseline_comb = [baseline_a, baseline_v]

baseline_v_reg = hddm.HDDMRegressor(data, baseline_v)
baseline_v_reg.sample(nsample, burn = nburn, dbname = 'baseline_v_reg.db', db = 'pickle')
baseline_v_reg.save('baseline_v_reg.model')

baseline_a_reg = hddm.HDDMRegressor(data, baseline_a)
baseline_a_reg.sample(nsample, burn = nburn, dbname = 'baseline_a_reg.db', db = 'pickle')
baseline_a_reg.save('baseline_a_reg.model')

baseline_comb_reg = hddm.HDDMRegressor(data, baseline_comb)
baseline_comb_reg.sample(nsample, burn = nburn, dbname = 'baseline_comb_reg.db', db = 'pickle')
baseline_comb_reg.save('baseline_comb_reg.model')


KJFLSDKLFSDAJFLAJSDJADSLKFJAL;SDKFJASLD;KJASDL;FKJASDL;KFJASD;LKFJASDFLKJSDFKLSJDF