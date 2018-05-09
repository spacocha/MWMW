#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 11:23:35 2017

@author: login
"""
## Load required libraries 

from calibration_functions import cal_iteration, run_model
import sys, os
import cPickle as pickle
import pandas as pd
import numpy as np

# this is where all the input data & output data lives
results_loc = "../Data/calibration_results2"
in_data_loc = "../Data/calibration_data"

# this is a table with 2n+2 columns (n = # of models)
# odd columns hold initial conditions for parameters
# even columns specify whether a particular param should be optimized
# the last two columns defines the upper and lower limits
settings = pd.read_csv(in_data_loc+"/calibration_run_settings.tsv", index_col=0, sep="\t")

# load the chem observations & practice data
chem_df = pd.read_csv(in_data_loc+"/chem_data.tsv", sep="\t", index_col=[0,1], parse_dates=True)
gene_df = pd.read_csv(in_data_loc+"/gene_data.tsv", sep="\t", index_col=[0,1], parse_dates=True)
obs_data = chem_df.join(gene_df)

# these are the number of trials in each run 
n_samplings = [1600, 1600,  1600,  800,  800,  800,  800,  800,  400,  400,  400,  400,  100,  100,  100,  100, 100]
n_samplings = [10]*17
# make the output if necessary
if not os.path.exists(results_loc):
    os.mkdir(results_loc)

# these two will hold dataframes showing how the std of params change 
# & the actual settings/scores change for each model
dev_trace, param_trace = {}, {}

for model_run_col in [0, 2]:
    # pickup model type from column name
    mod_type = settings.columns[model_run_col]
    # create a new location for this data
    results_subdir = results_loc + "/" + mod_type
    if not os.path.exists(results_subdir):
        os.mkdir(results_subdir)

    # pickup adjacent colname 
    bool_key = settings.columns[model_run_col+1]

    # create a list of parameters to search
    to_optimize = settings[settings[bool_key]].index

    # the algorithm will run an iteration for each parameter
    iterations = settings[bool_key].sum() + 1
    assert iterations <= len(n_samplings)

    # these are the data saved between runs (in case of crash)
    save_files = ['save'+"{:02}.p".format(i) for i in range(1,iterations+1)]
    
    # preallocate space for deviations per iteration & parameter value choices per iteration
    dev_trace[mod_type] = pd.DataFrame(index=to_optimize, columns=range(iterations))
    param_trace[mod_type] = pd.DataFrame(index=list(to_optimize)+["score"], columns=range(iterations+1))
    
    # initial run conditions
    last_score = run_model((settings.ix[:, mod_type], 
                            results_subdir+"/baseline/baseline_run.mat", 
                            obs_data.copy(), 
                            'gene_objective'))
    print "Initial score for {} is {:.2f}".format(mod_type, last_score)

    # load initial run score & parameter values 
    param_trace[mod_type].ix[to_optimize, 0] = settings.ix[to_optimize, mod_type]
    param_trace[mod_type].ix['score', 0] = last_score
    
    stds, converged_ = None, False

    for idx in range(1,iterations+1):
        trial_result = cal_iteration(idx, results_subdir, save_files[idx-1], settings, bool_key, obs_data, 
                                     n_samplings[idx-1], mod_type, stds, last_score, converged_, param_trace)
        stds, settings, last_score, winner, converged_ = trial_result
        dev_trace[mod_type].ix[stds.keys(), idx-1] = pd.Series(stds)
        param_trace[mod_type].ix[to_optimize, idx] = settings.ix[to_optimize, mod_type]
        param_trace[mod_type].ix['score', idx] = last_score
        
    print "Completed calibration of {}".format(mod_type)

dump_f = os.path.join(results_loc, 'final_result.p')
dump_c = (param_trace, dev_trace)
with open(dump_f, "wb") as dump_h:
    pickle.dump( dump_c, dump_h )

