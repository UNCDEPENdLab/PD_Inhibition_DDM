#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 12:46:07 2020

@author: natehall
"""

##Read dependencies - some borrowed from Sophie
import hddm
import os
from os import path
import pandas as pd
from pandas import Series
from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.pyplot as plt 
import pymp
from patsy import dmatrix
import kabuki
from kabuki.analyze import gelman_rubin
from kabuki.analyze import check_geweke
import pickle
import numpy as np
import sys
#import time

#time.perf_counter()
#
#
#time.perf_counter()

######### model checking and diagnostics'

ics = 1
task = 'recent_probes'
#task = 'flanker'
full_samples = 1 # denotes the number of samples. 0 is ~2000, and 1 will be full sampling around 20,000 draws.
nchains = 5


# which models need to be checked?

if task == 'flanker':
    models = ['v_reg','vsv_reg','v_block_reg','v_blocksv_reg','vst_reg','vsvst_reg','v_blockst_reg','v_blocksvst_reg']
    models = ['v_reg','vsv_reg','v_block_reg','v_blocksv_reg','vsvst_reg']
#    models = ['v_blockst_reg','v_blocksvst_reg']
elif task == 'recent_probes':
    models = ['v_reg','vsv_reg','vst_reg','vsvst_reg']

if ics:
    basedir = '/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM'
else:
    basedir = '/Users/natehall/ics/Nate/PD_Inhibition_DDM'

#os.chdir(basedir)

if full_samples:
    outputdir = basedir + '/Outputs/full_hddm/'+task
else:
    outputdir = basedir + '/Outputs/practice_hddm/'+task

os.chdir(outputdir)

# overwrite task naming so models will load properly for recent probes
if task == 'recent_probes':
    task = 'rp'


#### loop over models and chains to assess convergence.

# 3/22/20: check model convergence via gelman-rubin r-hat statistic, export traces of one chain, and save table of DIC values for model comparison.

dics = []
for mod in models:
    print "Assessing model fit for: ", mod

    mods = []
#    with pymp.Parallel(nchains) as ch:
    for chain in range(0,nchains):
        print chain
        this_model = hddm.load(mod+'_'+task+'_chain'+str(chain)+'.model') 
        mods.append(this_model)
#    sub_dict["models"] = models
    gel_rub =pd.DataFrame(gelman_rubin(mods).items(), columns = ['parameter', 'rhat']) 
    gel_rub.to_csv(outputdir+'/diagnostics/gr_' + mod + '.csv')
    dic = this_model.dic
    print dic
    dics.append(dic)
    traces = this_model.get_traces()
    traces.to_csv(outputdir+'/diagnostics/'+mod+'_traces.csv')
  
    
dics_export = {'model': models,
               'DIC': dics}

dics_exp = pd.DataFrame(dics_export, columns = ['model', 'DIC'])    
#dics = pd.DataFrame(dics, models)
dics_exp.to_csv(outputdir+'/diagnostics/dics_all.csv')









#
##    gr_descriptives = gel_rub.describe()
##    this_model.plot_posteriors()
##    x = this_model.print_stats()
##    this_model.plot_posterior_predictive()
#    
#    
##    uber_dict[mod] = sub_dict
##
##f = open("all_chains_gelman_rubins.pkl","wb")
##pickle.dump(uber_dict,f)
##f.close()    
#x = np.zeros(models)
#x[0] = dic
#x
#
#x = pd.DataFrame(np.zeros(len(models)), models)#, columns = 'model')
#x[0,1] = dic
#
#x = pd.DataFrame()
#
#x = []
#x.append(dic)
#x.append(123)
