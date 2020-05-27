#!/Users/natehall/ics/Nate/PD_Inhibition_DDM/Code/DDM/gng_sim_test/hddm_local/bin
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 12:46:07 2020

@author: natehall
"""

##Read dependencies - some borrowed from Sophie
import hddm
import os
#from os import path
import pandas as pd
#from pandas import Series
#from matplotlib.backends.backend_agg import FigureCanvasAgg
#import matplotlib.pyplot as plt 
import pymp
#from patsy import dmatrix
#import kabuki
from kabuki.analyze import gelman_rubin
#from kabuki.analyze import check_geweke
#import pickle
#import numpy as np
#import sys

#import time
from kimchi import convert #for re-pickling

#time.perf_counter()
#
#
#time.perf_counter()

######### model checking and diagnostics'


####################################
######### input parameters
## this should follow similar naming conventions to the output of the R master script.
ics = 1
tasks = ['flanker']#, 'recent_probes', 'go_nogo']
full_sample = ['full_sample']#, 'clean_sample'] 
#nsamp = ['samp2000']#, 'samp10000']
nsamp = ['samp10000']
code = "acc"
#date = "2020-05-18" #in case you want to just get models from a specific date 
use_log = 0
use_model_csv = 1
####################################
#########





if ics:
    basedir = '/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM'
else:
    basedir = '/Users/natehall/ics/Nate/PD_Inhibition_DDM'


# which models need to be checked?

#pull completed log and see who needs to be diagnosed


if use_log:
    log = pd.read_csv(basedir + '/Code/DDM/pbs_outputs/PD_Inhibition_DDM_job_info_log_completed.csv')
    log_finished = log[(log['outputs_located']== 'all')]
elif use_model_csv:
    print('reading csv')
    check_these = pd.read_csv(basedir + '/Data/cache/check_these_models.csv')
else:
    print('user must rework script so we know which diagnostics to grab.')

try:
    log_finished = log_finished[(log_finished['DAY_SUB'] == date)]
except NameError:
    print('no date defined')

 

  

#### loop over models and chains to assess convergence.

# 3/22/20: check model convergence via gelman-rubin r-hat statistic, export traces of one chain, and save table of DIC values for model comparison.
# 5/6/20: expand to encompass log structure for those models that have finished running completely.
# 5/25/20: since 

for samp in nsamp:
    print(samp)
    for task in tasks:
        print(task)
        for samp_size in full_sample:
            print(samp_size)
            
            if use_log:
                log_trimmed = log_finished[(log_finished['NSAMP'] == samp) & (log_finished['TASK'] == task) & (log_finished['SAMPLE'] == samp_size)]
                models = set(list(log_trimmed['MODEL'])) #set function should get unique values to not redo model grabbing
            elif use_model_csv:
                x = check_these[(check_these['task'] == task) & (check_these['samples'] == samp)]
                models = set(list(x['model']))
            
            outputdir = basedir + '/../HDDM_outputs_PD_Inhibition/' + samp + '/' + task + '/' + samp_size + '/model_objects' 
            os.chdir(outputdir)
            print('Navigated to ' + os.getcwd())
            
            print('Directory contents:\n')
            contents = os.listdir('.')
            print(contents)
            
            ## initialize dic df to append to.
            try:
                dics = pd.read_csv(outputdir+'/../diagnostics/dics_all.csv')[['model', 'DIC']]
            except:
                dics = []
               
     
            for mod in models:
                print "Assessing model fit for: ", mod
                
                if use_log:
                    nchains = max(log_trimmed['NCHAINS'][log_trimmed['MODEL'] == mod])
                elif use_model_csv:
#                    chain_list = [x for x in contents if x.startswith(mod + '_chain') & x.endswith('.model')] 
#                    res = [int(sub.split('_')[1]) for sub in chain_list]
#                    # using join() + isnumeric() + list comprehension + map() 
#                    # Extracting numbers from list of strings 
#                    res = list(map(lambda sub:int(''.join( 
#                          [ele for ele in sub if ele.isnumeric()])), chain_list)) 
                    nchains = 10
                
                mods = []
                
            #    with pymp.Parallel(nchains) as ch:
                
                for chain in range(0,nchains):
                    print chain
                    try:
                        this_model = hddm.load(mod+'_chain'+str(chain) + '_'+ code + 'Code.model') 
                        mods.append(this_model)
                        traces = this_model.get_traces()
                        traces.to_csv(outputdir+'/../diagnostics/'+mod+'_traces_'+ str(chain) +'.csv')
                    except Exception as e:
                        print(e)
                        if str(e) == 'unsupported pickle protocol: 3':
                            try:
                                convert(mod+'_chain'+str(chain) + '_'+ code + 'Code.model')
                                mods.append(this_model)
                                traces = this_model.get_traces()
                                traces.to_csv(outputdir+'/../diagnostics/'+mod+'_traces_'+ str(chain) +'.csv')
                            except Exception as ee:
                                print(ee)
                            
                        

                try:
                    gel_rub =pd.DataFrame(gelman_rubin(mods).items(), columns = ['parameter', 'rhat']) 
                    gel_rub.to_csv(outputdir+'/../diagnostics/gr_' + mod + '.csv')
                    dic = this_model.dic
                    print dic
                    dics = dics.append(pd.DataFrame([[mod,dic]], columns = ['model', 'DIC']))
                    
                    dics.to_csv(outputdir+'/../diagnostics/dics_all.csv')                
                except Exception as e:
                    print('Ran into issues with these settings: ' + samp + ' ' + task + ' ' + samp_size)
                
#              
#                
#            dics_export = {'model': models,
#                           'DIC': dics}
#            
#            dics_exp = pd.DataFrame(dics_export, columns = ['model', 'DIC'])    
            #dics = pd.DataFrame(dics, models)










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
