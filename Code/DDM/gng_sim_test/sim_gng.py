#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 15:33:03 2020

@author: natehall
"""

#### simulation and recovery test for go/no-go task using the HDDMRegressor framework.


import hddm
import numpy as np
import pandas as pd



## for reference, this is what our data looks like:
#mydat = hddm.load_csv( "~/github_repos/PD_Inhibition_DDM/Data/preprocessed/go_nogo_clean_sample_nafilt_accCode.csv")


## parameters for simulation

v_conditions = {'OneGo': 2, #'ThreeGo': 1.75, 'FiveGo': 1.5, 
                'SevenGo':1.25,
                'OneNoGo': 2, #'ThreeNoGo': 1.75, 'FiveNoGo': 1.5,
                'SevenNoGo':1.25} # drift rate will vary by condition
a_int = 1.5
t_int = .3
n_conditions = 12 # number of each condition
n_subs = 100



# In this scenario all conditions have the same number of trials, which is quite different than what we will be feeding HDDM, though it's probably worth starting here.
# In our data, the number of gos increases in conditions that are not "OneGo". Thus, there will be 3x12 = 36 go trials in the "ThreeGo" condition and 5x12 = 60 go trials in the "FiveGo" condition. 
# For now, let's not worry about this.
 
print('simulating data') 
data, params = hddm.generate.gen_rand_data(params={        
                                        'OneNoGo': {'v': v_conditions['OneNoGo'], 'a': a_int, 't': t_int},
#                                        'ThreeNoGo': {'v': v_conditions['ThreeNoGo'], 'a': a_int, 't': t_int}, # keep things simple initially
#                                        'FiveNoGo': {'v': v_conditions['FiveNoGo'], 'a': a_int, 't': t_int},
                                        'SevenNoGo': {'v': v_conditions['SevenNoGo'], 'a': a_int, 't': t_int},
                                        
                                        'OneGo': {'v': v_conditions['OneGo'], 'a': a_int, 't': t_int},
#                                        'ThreeGo': {'v': v_conditions['ThreeGo'], 'a': a_int, 't': t_int},
#                                        'FiveGo': {'v': v_conditions['FiveGo'], 'a': a_int, 't': t_int},
                                        'SevenGo': {'v': v_conditions['SevenGo'], 'a': a_int, 't': t_int}
                                        },  size = n_conditions, subjs = n_subs)


# search through conditions and relabel rts to correspond to no-gos (signs will then get flipped to designate correct vs incorrect no-gos)

data['rt'] = np.where((data['response'] == 1 & pd.Series(data['condition']).str.contains('No')),# or # correct no-go
#                (data['response'] == 0 & ~pd.Series(data['condition']).str.contains('No')), # incorrect go (omission error)
                999, data['rt'])

# I have no idea for the life of me how the above works, but it seems to be fine after checking so I won't dwell on it for now. moving on...

# flip errors
data = hddm.utils.flip_errors(data)

# standard HDDM: no regression
#v_standard = hddm.HDDM(data, depends_on={'v': 'condition'})

# Don't use. Produced the following error:
# RuntimeWarning: invalid value encountered in double_scalars
#  tmp2 = (x - v) * (fx - fw)
#Warning: Powell optimization failed. Falling back to simplex.

#v_standard.sample(100, burn=10, dbname='gng_standard.db', db='pickle') # for quicker debugging.
#v_standard.sample(2000, burn=400, dbname='gng_standard.db', db='pickle')
#v_standard.save('gng_standard_model')

#import pdb; pdb.set_trace()
# load fitted model
#v_standard = hddm.load('gng_standard_model')



# simple regressoion model

v_regression = {'model': "v ~ 1  + C(condition, Treatment('OneGo'))", 'link_func': lambda x: x}
v_reg = hddm.HDDMRegressor(data, v_regression, group_only_regressors=['true'])

#v_reg.find_starting_values()
v_reg.sample(10, burn=2, dbname='gng_regressor.db', db='pickle')

v_reg.save('gng_regressor')
import pdb; pdb.set_trace()


