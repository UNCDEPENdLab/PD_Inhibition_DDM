# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 08:43:45 2017

@author: nth7
"""

###rerun between subjects effects using proper between subjects effects

##load modules
import numpy as np
import hddm
import os
import pandas as pd
from pandas import Series
import matplotlib.pyplot as plt


# Load data from csv file into a NumPy structured array
os.chdir('/Users/nth7/Box Sync/EF SNAP HDDM')
data = hddm.load_csv('Flanker_reduced.csv')
data.head(10)


####################################################################################################################################
##look at RT dists of each individual subject
data = hddm.utils.flip_errors(data) 
SSfig = plt.figure()
ax = SSfig.add_subplot(111, xlabel = 'RT', ylabel = 'count', title = 'RT distributions')
for i, subj_data in data.groupby('subj_index'):
    subj_data.rt.hist(bins = 20, histtype = 'step', ax=ax)
    
plt.savefig('hddm_flanker_subj_RTs.pdf')
####################################################################################################################################
##tailor an individual heirarchical DDM around our flanker dataset. 
simp = hddm.HDDM(data)
simp.find_starting_values
simp.sample(5000, burn = 50, thin = 3, dbname = 'traces.db', db = 'pickle')

simpStats = simp.gen_stats()
simpStats
simp.plot_posteriors()

#Gelman-Rubin R Statistic,values should be close to 1 indicating that multiple samplings of the different chains are essentially identical
models = list() 
for i in range(5):
    m = hddm.HDDM(data) 
    m.find_starting_values() 
    m.sample(5000, burn = 50, thin = 3, dbname = 'traces.db', db = 'pickle')
    models.append(m)
hddm.analyze.gelman_rubin(models)
models.save('gelman_rubin_models')
m_within_subj_simp.save('model_within_subj_simp')

####################################################################################################################################
##separate v and a by stimulus 
simp_a_stim = hddm.HDDM(data, depends_on={'a': 'stim'})
simp_a_stim.find_starting_values()
simp_a_stim.sample(10000, burn = 1000, thin = 3, dbname = 'traces.db', db = 'pickle')

a_C, a_I = simp_a_stim.nodes_db.node[['a(CONGRUENT)', 'a(INCONGRUENT)']]
hddm.analyze.plot_posterior_nodes([a_C, a_I])
plt.xlabel('threshold')
plt.ylabel('Posterior probability')
plt.title('Posterior of threshold values by stimulus grouping')
plt.savefig('posteriorDist_simp_a_stim.pdf')

simp_v_stim = hddm.HDDM(data, depends_on={'v': 'stim'})
simp_v_stim.find_starting_values()
simp_v_stim.sample(10000, burn = 1000, thin = 3, dbname = 'traces.db', db = 'pickle')

v_C, v_I = simp_v_stim.nodes_db.node[['v(CONGRUENT)', 'v(INCONGRUENT)']]
hddm.analyze.plot_posterior_nodes([v_C, v_I])
plt.xlabel('drift rate')
plt.ylabel('Posterior probability')
plt.title('Posterior of drift rates by stimulus grouping')
plt.savefig('posteriorDist_simp_v_stim.pdf')

####Doing significance testing directly on the posterior of a and v including proportions and information criteria
print ("P(CONGRUENT > INCONGRUENT)_a = ", (a_C.trace() > a_I.trace()).mean())
print ("P(CONGRUENT > INCONGRUENT)_v = ", (v_C.trace() > v_I.trace()).mean())
print ("Lumped model DIC: %f" % simp.dic)
print ("Lumped model BIC: %f" % simp.bic)
print ("Lumped model AIC: %f" % simp.aic)
print ("threshold by stimulus DIC: %f" % simp_a_stim.dic)
print ("threshold by stimulus BIC: %f" % simp_a_stim.bic)
print ("threshold by stimulus AIC: %f" % simp_a_stim.aic)
print ("drift rate by stimulus DIC: %f" % simp_v_stim.dic)
print ("drift rate by stimulus BIC: %f" % simp_v_stim.bic)
print ("drift rate by stimulus AIC: %f" % simp_v_stim.aic)
####################################################################################################################################
##modelling within-subjects effects (rather than sampling from separate gorup priors with depends_on)
from patsy import dmatrix
dmatrix("C(stim, Treatment('CONGRUENT'))", data.head(10))

within_v = hddm.HDDMRegressor(data, "v ~ C(stim, Treatment('CONGRUENT'))")
within_v.find_starting_values()
within_v.sample(10000, burn = 100, thin = 3, dbname = 'traces.db', db = 'pickle')

v_C, v_I = within_v.nodes_db.ix[["v_Intercept", "v_C(stim, Treatment('CONGRUENT'))[T.INCONGRUENT]"], 'node']
hddm.analyze.plot_posterior_nodes([v_C, v_I])
plt.xlabel('drift-rate')
plt.ylabel('Posterior probability')
plt.title('Group mean posteriors of within-subject drift-rate effects')
plt.savefig('posteriorDist_within_v.pdf')

within_a = hddm.HDDMRegressor(data, "a ~ C(stim, Treatment('CONGRUENT'))")
within_a.find_starting_values()
within_a.sample(10000, burn = 100, thin = 3, dbname = 'traces.db', db = 'pickle')

a_C, a_I = within_a.nodes_db.ix[["a_Intercept", "a_C(stim, Treatment('CONGRUENT'))[T.INCONGRUENT]"], 'node']
hddm.analyze.plot_posterior_nodes([a_C, a_I])
plt.xlabel('threshold')
plt.ylabel('Posterior probability')
plt.title('Group mean posteriors of within-subject threshold effects')
plt.savefig('posteriorDist_within_a.pdf')

####################################################################################################################################
#setting up HDDM between person regression models
####################################################################################################################################

#Z link function fixes z to values within 0 and 1
def z_link_func(x, data=data):
    stim = (np.asarray(dmatrix('0 + C(s, [[1], [-1]])',
                               {'s': data.stimulus.ix[x.index]}))
    )
    return 1 / (1 + np.exp(-(x * stim)))

###start by fitting a simple between subjects design predicting drift rate across stimuli by suicide proneness, add in conditions later
bw_sui_v = {'model': 'v ~ 1 + SUICPRON', 'link_func': lambda x: x}
bw_sui_v_reg = hddm.HDDMRegressor(data, bw_sui_v)
bw_sui_v_reg.sample(30000, burn = 3000, thin =5, dbname = 'traces.db', db = 'pickle')
bw_sui_v_reg.gen_stats()
bw_sui_v_reg.plot_posteriors()
bw_sui_v_reg.plot_posterior_predictive()

###okay, now never go back to this^^^
####################################################################################################################################
####SUICIDE PRONENESS
###now, with the same model as above, add variation within individuals at the level of stimulus condition (congruent vs incongruent)
bw_wi_sui_v = {'model': "v ~ 1 + SUICPRON + C(stim, Treatment('CONGRUENT'))", 'link_func': lambda x: x}
bw_wi_sui_v_reg = hddm.HDDMRegressor(data, bw_wi_sui_v)
bw_wi_sui_v_reg.sample(30000, burn = 3000, thin =5, dbname = 'bw_wi_sui_v_reg.db', db = 'pickle')
SUI_V = bw_wi_sui_v_reg.gen_stats()
SUI_V
bw_wi_sui_v_reg.plot_posteriors()
bw_wi_sui_v_reg.plot_posterior_predictive()

###start by fitting a simple between subjects design predicting drift rate across stimuli by suicide proneness, add in conditions later
bw_wi_sui_a = {'model': "a ~ 1 + SUICPRON + C(stim, Treatment('CONGRUENT'))", 'link_func': lambda x: x}
bw_wi_sui_a_reg = hddm.HDDMRegressor(data, bw_wi_sui_a)
bw_wi_sui_a_reg.sample(30000, burn = 3000, thin =5, dbname = 'bw_wi_sui_a_reg.db', db = 'pickle')
SUI_A = bw_wi_sui_a_reg.gen_stats()
SUI_A
bw_wi_sui_a_reg.plot_posteriors()
bw_wi_sui_a_reg.plot_posterior_predictive()

#####COMBINED A and V model
sui_comb = [bw_wi_sui_a, bw_wi_sui_v]
sui_comb_reg = hddm.HDDMRegressor(data, sui_comb)
sui_comb_reg.sample(40000, burn = 3000, thin =5, dbname = 'sui_comb_reg.db', db = 'pickle')
SUI_AV = sui_comb_reg.gen_stats()
SUI_AV
sui_comb_reg.plot_posteriors()
sui_comb_reg.plot_posterior_predictive()

####DISINHIBITION
###now, with the same model as above, add variation within individuals at the level of stimulus condition (congruent vs incongruent)
bw_wi_disin_v = {'model': "v ~ 1 + DISINH + C(stim, Treatment('CONGRUENT'))", 'link_func': lambda x: x}
bw_wi_disin_v_reg = hddm.HDDMRegressor(data, bw_wi_disin_v)
bw_wi_disin_v_reg.sample(30000, burn = 3000, thin =5, dbname = 'bw_wi_disin_v_reg.db', db = 'pickle')
DISIN_V = bw_wi_disin_v_reg.gen_stats()
DISIN_V
bw_wi_disin_v_reg.plot_posteriors()
bw_wi_disin_v_reg.plot_posterior_predictive()

###start by fitting a simple between subjects design predicting drift rate across stimuli by suicide proneness, add in conditions later
bw_wi_disin_a = {'model': "a ~ 1 + DISINH + C(stim, Treatment('CONGRUENT'))", 'link_func': lambda x: x}
bw_wi_disin_a_reg = hddm.HDDMRegressor(data, bw_wi_disin_a)
bw_wi_disin_a_reg.sample(30000, burn = 3000, thin =5, dbname = 'bw_wi_disin_a_reg.db', db = 'pickle')
DISIN_A = bw_wi_disin_a_reg.gen_stats()
DISIN_A
bw_wi_disin_a_reg.plot_posteriors()
bw_wi_disin_a_reg.plot_posterior_predictive()

#####COMBINED A and V model
disin_comb = [bw_wi_disin_a, bw_wi_disin_v]
disin_comb_reg = hddm.HDDMRegressor(data, sui_comb)
sui_comb_reg.sample(40000, burn = 3000, thin =5, dbname = 'sui_comb_reg.db', db = 'pickle')
SUI_AV = sui_comb_reg.gen_stats()
SUI_AV
sui_comb_reg.plot_posteriors()
sui_comb_reg.plot_posterior_predictive()



























####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
##Misc code to be implemented above
simp = {'model': "a", 'link_func': lambda x: x}

simp_reg = hddm.HDDMRegressor(data, asimp, group_only_regressors=['true'])
simp_reg.find_starting_values()
simp_reg.sample(30000, burn = 1000, thin = 5, dbname = 'traces.db', db = 'pickle')

simp_reg.print_stats() # check values of reg coefficients against the generated ones
simp_reg.plot_posteriors()
simp_reg.plot_posterior_predictive()
    
    
####################################################################################################################################
##simple within-subjects model of a and v (no between subjects factors, simply condition)
a_reg_within = {'model': "a ~ 1 + stim", 'link_func': lambda x: x}
v_reg_within = {'model': "v ~ 1 + stim", 'link_func': lambda x: x}

reg_comb_within = [a_reg_within, v_reg_within]

m_reg_within = hddm.HDDMRegressor(data, reg_comb_within, group_only_regressors = ['true'])
m_reg_within.find_starting_values()
m_reg_within.sample(30000, burn = 5000, thin =5, dbname = 'traces.db', db = 'pickle')
m_reg_within.print_stats()
m_reg_within.plot_posteriors()

    
    
    
#z_reg = {'model': 'z ~ 1 + SUICPRON', 'link_func': z_link_func}
#z_reg_within = {'model': 'z ~ 1 + SUICPRON + C(condition)', 'link_func': z_link_func}
a


a_reg = {'model': 'a ~ 1 + SUICPRON', 'link_func': lambda x: x}
a_reg_within = {'model': "a ~ 1 + SUICPRON + C(stim, Treatment('CRFR'))", 'link_func': lambda x: x}

v_reg = {'model': 'v ~ 1 + SUICPRON', 'link_func': lambda x: x}
v_reg_within = {'model': "v ~ 1 + SUICPRON + C(stim, Treatment('CRFR'))", 'link_func': lambda x: x}


reg_comb = [a_reg, v_reg]
#m_reg = hddm.HDDMRegressor(data_group, reg_comb, group_only_regressors=['true']) 

m_reg = hddm.HDDMRegressor(data, reg_comb, group_only_regressors=['true'])
m_reg.find_starting_values()
m_reg.sample(10000, burn = 1000, thin = 5, dbname = 'traces.db', db = 'pickle')

m_reg.print_stats() # check values of reg coefficients against the generated ones
m_reg.plot_posteriors()
m_reg.plot_posterior_predictive()





####################################################################################################################################
#Gelman-Rubin R Statistic,values should be close to 1 indicating that multiple samplings of the different chains are essentially identical
models = list() 
for i in range(5):
    m = hddm.HDDM(data) 
    m.find_starting_values() 
    m.sample(5000, burn=20) 
    models.append(m)
hddm.analyze.gelman_rubin(models)
model.plot_posteriors()
model.plot_posterior_predictive()