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
settings = pd.read_csv("../Data/calibration_data/calibration_run_settings.tsv", index_col=0, sep="\t")

# load the chem observations & practice data
chem_df = pd.read_csv(in_data_loc+"/chem_data.tsv", sep="\t", index_col=[0,1])
gene_df = pd.read_csv(in_data_loc+"/gene_data.tsv", sep="\t", index_col=[0,1])
super_obs = chem_df.join(gene_df)

# the algorithm will run an iteration for each parameter
iterations = settings.new_model_bool.sum()

# these are the data saved between runs (in case of crash)
save_files = ['save'+"{:02}.p".format(i) for i in range(1,iterations+1)]

# these are the data loaded between runs (necessary?)
load_files = ["init"] + save_files[:-1]

# these are the number of trials in each run 
n_samplings = [1600, 1500,  800,  800,  800,  800,  800,  400,  400,  400,  400,  400,  100,  100,  100,  100]


# these are types of optimizations we will run. 
# (new = with genes) and (old = w/o genes)
optimize_option = init_conds.columns.tolist()

# make the output if necessary
if not os.path.exists(results_loc):
    os.mkdir(results_loc)

