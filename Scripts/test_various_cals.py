
from calibration_functions import cal_iteration, run_model
import cPickle as pickle
import sys, os
import pandas as pd
import numpy as np

results_loc = "../Data/calibration_results"
in_data_loc = "../Data/calibration_data"
settings = pd.read_csv(in_data_loc+"/calibration_run_settings.tsv", index_col=0, sep="\t")
old_settings = pd.read_csv(in_data_loc+"/old_params.csv", index_col=0)


chem_df = pd.read_csv(in_data_loc+"/chem_data.tsv", sep="\t", index_col=[0,1], parse_dates=True)
gene_df = pd.read_csv(in_data_loc+"/gene_data.tsv", sep="\t", index_col=[0,1], parse_dates=True)
obs_data = chem_df.join(gene_df)

mod_type = settings.columns[0]
results_subdir = results_loc + "/" + mod_type

oldest_baseline = old_settings.ix[settings.index, 'Previous model parameters']
last_score = run_model((oldest_baseline, results_subdir+"/baseline/baseline_run.mat", obs_data.copy(), 'print_objective'))
older_baseline = old_settings.ix[settings.index, 'Without gene calibration']
last_score = run_model((older_baseline, results_subdir+"/baseline/baseline_run.mat", obs_data.copy(), 'print_objective'))
old_baseline = old_settings.ix[settings.index, 'With gene calibration']
last_score = run_model((old_baseline, results_subdir+"/baseline/baseline_run.mat", obs_data.copy(), 'print_objective'))