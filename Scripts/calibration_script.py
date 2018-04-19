#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 11:23:35 2017

@author: login
"""
## Load required libraries 

import sys, os
import pandas as pd

# this is where all the input data & output data lives
results_loc = "../Data/calibration_results"
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
n_samplings = [1600, 1500,  800,  800,  800,  800,  800,  400,  400,  400,  400,  400,  100,  100,  100,  100]

# make the output if necessary
if not os.path.exists(results_loc):
    os.mkdir(results_loc)

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
    iterations = settings[bool_key].sum()
    assert iterations <= len(n_samplings)
    # these are the data saved between runs (in case of crash)
    save_files = ['save'+"{:02}.p".format(i) for i in range(1,iterations+1)]
    # initial run conditions
    stds = None
    last_score = -1.01
    converged_ = False

    for idx in xrange(1,iterations+1):
        trail_result = run_calibration(idx, results_subdir, save_files[idx-1], settings, bool_key, obs_data, 
        	                           n_samplings[idx-1], mod_type, stds, last_score, converged_)
        stds, last_score, convergence_bool, settings = trail_result






