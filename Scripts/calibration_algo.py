#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 11:23:35 2017

@author: login
"""
## Load required libraries 

pd.options.mode.chained_assignment = None
from copy import copy
import cPickle as pickle
import pandas as pd
import numpy as np
import os, sys
from calibration_functions import run_model, sample_params, best_val
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count


def run_calibration(idx, results_loc, save_file, settings, bool_key, obs_data, n_trials, mod_type, stds, last_score, converged_):
    # Calculate the number of threads possible
    thread_est = cpu_count()
    
    # Determine if it is initial run or a restart
    if idx == 1:
        print "Starting new optimization with settings for {}".format(mod_type)
        init_bool = True
    else:
        init_bool = False
    
    # determine which variables are left to optimize
    not_optimized = settings[settings.ix[:, bool_key]].index.tolist()
    # create a directory for each model run output
    lake_dirs = [os.path.join(results_loc, "lake_"+str(i))  for i in range(1,n_trials+1)]
    # create a filename for the output
    result_str = ["run_{}.mat".format(i) for i in range(1,n_trials+1)]
    # create paths for the output
    result_paths = [os.path.join(j, i) for i, j in zip(result_str, lake_dirs)]
    # fill a df with param choices 
    parameter_df = sample_params(settings, init_bool, n_trials, not_optimized, mod_type, stds)

    # run the model n_trials number of times if it hasn't converged
    if not converged_:
        score_vector = pd.Series(index=parameter_df.index, data=np.zeros((n_trials)))

        model_args = []
        for trial_n in range(1, n_trials+1):
            arg1 = parameter_df.ix[trial_n, :]
            arg2 = result_paths[trial_n-1]
            arg3 = obs_data.copy()
            arg4 = 'gene_objective'
            model_args.append((arg1, arg2, arg3, arg4))
        
        parallelism = True
        
        if not parallelism:
            score_list = []
            for d, p, o, s in model_args:
                print os.path.basename(p)
                score_list.append(run_model((d, p, o, s)))   
        else:
            score_list = []
            pool = ThreadPool(thread_est)
            results = pool.map(run_model, model_args)
            pool.close()
            pool.join()
            score_list=results
        
        score_vector[range(1, n_trials+1)] = score_list

    else:
        score_vector = pd.Series(index=parameter_df.index,
                                 data = np.ones((n_trials))*last_score)

parameter_df['score'] = score_vector
parameter_df.to_csv(parameter_df_f)
     
# Step 4: Perform an F-test to identify the most sensitive parameter & fix it
# Note: F-test is insensitive to range, but the model is most sensitive to the 
# range of the input values

from sklearn.feature_selection import f_regression
from sklearn.linear_model import LinearRegression

x_df = parameter_df.ix[:, to_optimize]
F_vals, p_vals = f_regression(x_df.values, score_vector.values)
p_series = pd.Series(index=x_df.columns, data=p_vals)

# check to see if no parameters make a difference
if p_series.isnull().sum() == len(p_series):
    p_series[0] = 0.05
    convergence_bool = True

winner = p_series[p_series == p_series.min()].index[0]

# Step 5: Take the best value for that parameter and those not requiring optimization
# Ensure this run beats the last run
new_params = {}

if last_score > parameter_df.score.max():
    previous_n = int(sys.argv[2])-1
    previous_pdf_f = os.path.join(results_loc, "trial_data_{}.csv".format(previous_n) )
    previous_pdf = pd.read_csv(previous_pdf_f, index_col=0)
    opt_value = best_val(previous_pdf, winner)
    print "Previous run had a better value, using that instead"
else:
    opt_value = best_val(parameter_df, winner)

new_params[winner] = opt_value

if not convergence_bool:
    print "{} is the most sensitive variable".format(winner)
else:
    print "Converged at maximum, no further optimization necessary"
    
to_optimize.remove(winner)
not_optimized.append(winner)

for n_o in not_optimized:
    new_params[n_o] = best_val(parameter_df, n_o)

# Step 6: Perform regression on the remaining optimizable params 
std_deviations = {}
for t_o in to_optimize:
    lr = LinearRegression()
    y = parameter_df.ix[:, 'score'].values
    x = parameter_df.ix[:, t_o].values
    lr.fit(y.reshape(-1, 1),x)
    best_guess = lr.predict(np.array([[y.max()]]))
    new_params[t_o] = best_guess
    std_deviations[t_o] = abs(best_guess - params[t_o])

this_score = parameter_df.score.max()
improvement = (this_score-last_score)
print "{:.2f}, a change of {:.1%}".format(this_score, improvement)

for o_key in params.keys():
    params[o_key] = new_params[o_key]

dump_f = os.path.join(results_loc, dump_bn)
dump_c = (std_deviations, to_optimize, this_score, params, convergence_bool)
with open(dump_f, "wb") as dump_h:
    pickle.dump( dump_c, dump_h )
