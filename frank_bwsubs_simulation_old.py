'''
Created on Feb 8, 2017

@author: michael frank (amended by MNH)
'''
import hddm
from numpy import mean, std
#from patsy import dmatrix  
import numpy as np         
from pandas import Series
import pandas as pd
import os as os
import matplotlib.pyplot as plt
#from matplotlib.pyplot import scatter

os.chdir('/Users/michael/Box_Sync/DEPENd/Projects/PD_Inhibition_DDM')

## Generate data

beta_a = 0.4  # between subjects scaling factor - adjust this and should be able to recover its true value
a_int = 1  # set intercept within range of empirical priors 
v_int = 0.3
x_range = range(11)  # set values of between subject measures, here from 0 to 10
trials_per_level = 200    
subjs_per_bin = 10

data_group = pd.DataFrame()  # empty df to append to  

for x in x_range:
    xx = (x - mean(x_range)) / std(x_range)  # z-score the x factor
    a = a_int + beta_a * xx  #  indiv subj param values that are centered on intercept but deviate from it up or down by z-scored x
    # v = v_int+ beta_a*xx  # can also do for drift, here using same beta coeff
    
    # parvec = {'v':.3, 'a':a, 't':.3, 'sv':0, 'z':.5, 'sz':0, 'st':0}
    # parvec2 = {'v':.3, 'a':a+.5, 't':.3, 'sv':0, 'z':.5, 'sz':0, 'st':0}
    parvec = {'v':.3, 'a':a, 't':.3}  # set a to value set by regression, here v is set to constant

    # note that for subjs_per_bin > 1, these are just the mean values of the parameters; indiv subjs within bin are sampled from distributions with the given means, but can still differ within bin around those means. 
    #not including sv, sz, st in the statement ensures those are actually 0.

    data_a, params_a = hddm.generate.gen_rand_data({'level1': parvec}, size=trials_per_level, subjs=subjs_per_bin)
    
    # can also do with two levels of within-subj conditions
    # data_a, params_a = hddm.generate.gen_rand_data({'level1': parvec,'level2': parvec2}, size=trials_per_level, subjs=subjs_per_bin)

    data_a['z_x'] = Series(xx * np.ones((len(data_a))), index=data_a.index)  # store the (z-scored) between subjects factor in the data frame, same value for every row for each subject in the bin
    data_a['x'] = Series(x * np.ones((len(data_a))), index=data_a.index)  # store the (z-scored) between subjects factor in the data frame, same value for every row for each subject in the bin
    data_a['a_population'] = Series(a * np.ones((len(data_a))), index=data_a.index)  # store the (z-scored) between subjects factor in the data frame, same value for every row for each subject in the bin
    data_a['subj_idx'] = Series(x * subjs_per_bin + data_a['subj_idx'] * np.ones((len(data_a))), index=data_a.index)  # store the correct subj_idx when iterating over multiple subjs per bin
     
    # concatenate data_a with group data
    data_a_df = pd.DataFrame(data=[data_a])
    data_group = data_group.append([data_a], ignore_index=True)

data_group.to_csv('data_group.csv')

## Recover

a_reg = {'model': 'a ~ 1 + z_x', 'link_func': lambda x: x}
# a_reg_within = {'model': 'a ~ 1+x + C(condition)', 'link_func': lambda x: x}
# for including and estimating within subject effects of  condition

v_reg = {'model': 'v ~ 1 + z_x', 'link_func': lambda x: x}
reg_comb = [a_reg, v_reg]
# m_reg = hddm.HDDMRegressor(data_group, reg_comb, group_only_regressors=['true']) 

m_reg = hddm.HDDMRegressor(data_group, a_reg, group_only_regressors=['true'])
m_reg.find_starting_values()
m_reg.sample(3000, burn=500, dbname='a_bwsubs_t200.db', db='pickle')
m_reg.save('a_bwsubs_model_t200')

m_reg.print_stats()  # check values of reg coefficients against the generated ones

m_reg = hddm.load('a_bwsubs_model')
data_group = pd.read_csv('data_group.csv')

#look at correlation of recovered parameter with original
subjdf = data_group.groupby('subj_idx').first().reset_index()

## check for residual correlation with x 
a_int_recovered =[]
pp=[]

from scipy import stats
for i in range(0,(1+max(x_range))*subjs_per_bin):
    a='a_Intercept_subj.'
    a+=str(i)
    a+='.0'
    xx=i//subjs_per_bin
    p1=m_reg.nodes_db.node[a] 
    a_int_recovered.append(p1.trace().mean()) 
    pp.append(xx)
    
a_x_recovered = m_reg.nodes_db.node['a_x'].trace().mean()

subjdf.apred = a_int_recovered + a_x_recovered * subjdf.z_x
stats.pearsonr(subjdf.apred,subjdf.a_population) # correlation between predicted a value and population a value
plt.scatter(subjdf.apred,subjdf.a_population) #predicted versus observed a

stats.pearsonr(a_int_recovered,subjdf.x) # should be zero correlation if entirely accounted for by x regressor - this correlation should instead be positive if x is removed from the model fit

fig = plt.figure()
fig.set_size_inches(5, 4)
plt.scatter(pp,a_int_recovered)
plt.show()
plt.savefig('residual correlation between ranef and ax.png', dpi=300, format='png')
print('Pearson correlation between a_x and x') 
print(stats.pearsonr(pp,a_int_recovered)) # should be zero correlation if entirely accounted for by x regressor - this correlation should instead be positive if x is removed from the model fit
