from calibration_functions import cal_iteration, run_model
import os
import pandas as pd

d1 = "../Data/calibration_results3/old_model"
d2 = "../Data/calibration_results3/new_model"

f1s = os.listdir(d1)
p1s = [os.path.join(d1, i) for i in f1s if i.endswith(".csv")]
csv1s = pd.concat([pd.read_csv(i, index_col=0) for i in p1s], ignore_index=1)
old_row = csv1s.ix[csv1s.score.argmax(), :]
old_row.name = "old_model"

f2s = os.listdir(d2)
p2s = [os.path.join(d2, i) for i in f2s if i.endswith(".csv")]
csv2s = pd.concat([pd.read_csv(i, index_col=0) for i in p2s], ignore_index=1)
new_row = csv2s.ix[csv2s.score.argmax(), :]
new_row.name = "new_model"

param_settings = pd.concat((old_row, new_row), 1)
new_settings = param_settings.ix[param_settings.index[:-1], :]

in_data_loc = "../Data/calibration_data"
settings = pd.read_csv(in_data_loc+"/calibration_run_settings.tsv", index_col=0, sep="\t")
chem_df = pd.read_csv(in_data_loc+"/chem_data.tsv", sep="\t", index_col=[0,1], parse_dates=True)
gene_df = pd.read_csv(in_data_loc+"/gene_data.tsv", sep="\t", index_col=[0,1], parse_dates=True)
obs_data = chem_df.join(gene_df)

old_settings = pd.read_csv(in_data_loc+"/old_params.csv", index_col=0)
all_settings = pd.concat((new_settings, old_settings), 1)
all_settings.to_csv("../Data/final_calibration/final_params.csv")

mod_type = 'old_model'
last_score = run_model((new_settings.ix[:, mod_type], "../Data/final_calibration/old_final/old_run.mat", obs_data.copy(), 'gene_objective'))

mod_type = 'new_model'
last_score = run_model((new_settings.ix[:, mod_type], "../Data/final_calibration/new_final/new_run.mat", obs_data.copy(), 'gene_objective'))

mod_type = 'new_model'
oldest_baseline = old_settings.ix[settings.index, 'Previous model parameters']
last_score = run_model((oldest_baseline, "../Data/final_calibration/old_no_final/oldno_run.mat", obs_data.copy(), 'gene_objective'))

older_baseline = old_settings.ix[settings.index, 'Without gene calibration']
last_score = run_model((older_baseline, "../Data/final_calibration/old_chem_final/chem_run.mat", obs_data.copy(), 'gene_objective'))

old_baseline = old_settings.ix[settings.index, 'With gene calibration']
last_score = run_model((old_baseline, "../Data/final_calibration/old_gene_final/gene_run.mat", obs_data.copy(), 'gene_objective'))