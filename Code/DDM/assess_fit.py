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
import csv
#import time

#time.perf_counter()
#
#
#time.perf_counter()

######### model checking and diagnostics'

## this should follow similar naming conventions to the output of the R master script.
ics = 0
#task = 'recent_probes'
tasks = ['flanker', 'recent_probes']
full_sample = ['full_sample', 'clean_sample'] 
nsamp = ['samp2000', 'samp10000']
code = "acc"

if ics:
    basedir = '/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM'
else:
    basedir = '/Users/natehall/ics/Nate/PD_Inhibition_DDM'


# which models need to be checked?

#pull completed log and see who needs to be diagnosed
log = pd.read_csv(basedir + '/Code/DDM/pbs_outputs/PD_Inhibition_DDM_job_info_log_completed.csv')

log_finished = log[(log['outputs_located']== 'all')]


#### loop over models and chains to assess convergence.

# 3/22/20: check model convergence via gelman-rubin r-hat statistic, export traces of one chain, and save table of DIC values for model comparison.
# 5/6/20: expand to encompass log structure for those models that have finished running completely.

for samp in nsamp:
    print samp
    for task in tasks:
        print task
        for samp_size in full_sample:
            print samp_size
            log_trimmed = log_finished[(log_finished['NSAMP'] == samp) & (log_finished['TASK'] == task) & (log_finished['SAMPLE'] == samp_size)]
            models = list(log_trimmed['MODEL'])
            
            outputdir = basedir + '/../HDDM_outputs_PD_Inhibition/' + samp + '/' + task + '/' + samp_size + '/model_objects' 
            os.chdir(outputdir)
            print 'Navigated to ' + os.getcwd()
            
            print 'Directory contents:\n'
            print os.listdir('.')
            
            dics = []
            for mod in models:
                print "Assessing model fit for: ", mod
                
                nchains = int(log_trimmed['NCHAINS'][log_trimmed['MODEL'] == mod])
                mods = []
            #    with pymp.Parallel(nchains) as ch:
                for chain in range(0,nchains):
                    print chain
                    this_model = hddm.load(mod+'_chain'+str(chain) + '_'+ code + 'Code.model') 
                    mods.append(this_model)
            #    sub_dict["models"] = models
                gel_rub =pd.DataFrame(gelman_rubin(mods).items(), columns = ['parameter', 'rhat']) 
                gel_rub.to_csv(outputdir+'/../diagnostics/gr_' + mod + '.csv')
                dic = this_model.dic
                print dic
                dics.append(dic)
                traces = this_model.get_traces()
                traces.to_csv(outputdir+'/diagnostics/'+mod+'_traces.csv')
#              
#                
            dics_export = {'model': models,
                           'DIC': dics}
            
            dics_exp = pd.DataFrame(dics_export, columns = ['model', 'DIC'])    
            #dics = pd.DataFrame(dics, models)
            dics_exp.to_csv(outputdir+'/../diagnostics/dics_all.csv')









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
