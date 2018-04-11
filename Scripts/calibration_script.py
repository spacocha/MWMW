#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 11:23:35 2017

@author: login
"""
## Load required libraries 


from copy import copy
import cPickle as pickle
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import os, sys

from calibration_functions import run_model, fill_param_dict, best_val
from calibration_functions import load_param_dict, load_optimization_var_list
from observed_data import load_chem_data, load_gene_data
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count

## Read repo location on disk into memory and command line arguments
home_dir = os.getcwd()
if len(sys.argv) > 1:
    save_file = sys.argv[1]
    dump_bn = sys.argv[3]
    n_samplings = int(sys.argv[4])
    results_loc = os.path.join(home_dir, sys.argv[5])
    parameter_df_f = os.path.join(results_loc, "trial_data_{}.csv".format(sys.argv[2]) )
    optimize_option = sys.argv[6]
else:
    save_file = 'init'
    dump_bn = 'save01.p'
    n_samplings = int('4')
    results_loc = os.path.join(home_dir, 'testing')
    parameter_df_f = os.path.join(results_loc, "trial_data_1.csv")
    optimize_option = 'new_model'

if not os.path.exists(results_loc):
    os.mkdir(results_loc)

## Calculate the number of threads possible and whether previous run completed
thread_est = cpu_count()

results_loc_contents = os.listdir(results_loc)
red_flags = [i for i in results_loc_contents if i.startswith("lake_")]

if len(red_flags) > 0:
    sys.exit("cleanup needed")
    
# Step 0: Check for a save file to determine if it is initial or a restart

if save_file == 'init':
    print "Starting new optimization with settings for {}".format(optimize_option)
    init_bool = True
    stds = None
    if optimize_option == 'new_model':
        params = load_param_dict('midpoint')
    elif optimize_option == 'old_model':
        params = load_param_dict('default-midpoint')

    last_score = -1.01
    convergence_bool = False
    to_optimize = load_optimization_var_list(optimize_option)
    for p in to_optimize:
        print "\tWill optimize {} starting at {}".format(p, params[p])
    # Step 2: Define the parameter space for each variable
else:
    init_bool = False
    pickle_f = os.path.join(results_loc, save_file)
    with open( pickle_f, "rb" ) as pickle_h:
        pickle_pack = pickle.load(pickle_h)
    stds, to_optimize, last_score, params, convergence_bool = pickle_pack
    
not_optimized = [k for k in params.keys() if k not in to_optimize]
assert len(params.keys()) == (len(to_optimize) + len(not_optimized))

# Step 1: Load the chem observations & practice data

chem_dir = os.path.join(home_dir, 'ChemData')
obs_df, gene_df = load_chem_data(), load_gene_data()
super_obs = obs_df.join(gene_df)

# Step 3: Exacute n number of samplings per peramater, preload a df with this

lake_dirs = [os.path.join(results_loc, "lake_"+str(i))  for i in range(1,n_samplings+1)]
result_str = ["run_{}.mat".format(i) for i in range(1,n_samplings+1)]
result_paths = [os.path.join(j, i) for i, j in zip(result_str, lake_dirs)]
parameter_df = fill_param_dict(params, init_bool, n_samplings, 0.0, to_optimize, stds)

# Step 3.5 Run the model n number of times
if not convergence_bool:
    score_vector = pd.Series(index=parameter_df.index,data = np.zeros((n_samplings)))

    model_args = []
    for trial_n in range(1, n_samplings+1):
        arg1 = parameter_df.ix[trial_n, :]
        arg2 = result_paths[trial_n-1]
        arg3 = copy(super_obs)
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
    
    score_vector[range(1, n_samplings+1)] = score_list

else:
    score_vector = pd.Series(index=parameter_df.index,
                             data = np.ones((n_samplings))*last_score)

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
