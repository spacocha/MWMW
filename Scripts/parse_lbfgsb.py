import numpy as np
import pandas as pd
from scipy.optimize import minimize
from calibration_functions import run_model
import os

new_copied_res = np.array([  1.70396396e+02,   1.00000000e-01,   5.24324324e+04,
                             1.65535135e+04,   1.31594823e-01,   1.20969339e-01,
                             0.00000000e+00,   1.30891892e+02,   8.77371116e-01,
                             1.59224783e-01,   3.41341341e+01,   2.12939671e+00,
                             1.31173794e+00,   4.18799759e+00,   2.84775087e-01,
                             3.25698820e+00])
old_copied_res = np.array([ 1.67691692e+02,   4.63463504e-02,   9.90090090e+04,
                            1.55070270e+04,   6.76676625e-02,   1.46146135e-01,
                            7.77477470e+01,   2.74774758e+00,   2.21622388e+00,
                            5.12912874e-01,   5.98822777e+00,   2.79079058e-01,
                            8.02802797e+00])

in_data_loc = "../Data/calibration_data"
settings = pd.read_csv(in_data_loc+"/calibration_run_settings.tsv", index_col=0, sep="\t")
chem_df = pd.read_csv(in_data_loc+"/chem_data.tsv", sep="\t", index_col=[0,1], parse_dates=True)
gene_df = pd.read_csv(in_data_loc+"/gene_data.tsv", sep="\t", index_col=[0,1], parse_dates=True)
obs_data = chem_df.join(gene_df)

mod_types = [0, 2]
mod_res = [new_copied_res, old_copied_res]

these_settings = settings.copy()
for x, mi in zip(mod_res, mod_types):
    mt = settings.columns[mi]
    mb = settings.columns[mi+1]
    to_optimize = list(settings.ix[settings.ix[:, mb], mt].index)
    new_vals = pd.Series({i:j for i, j in zip(to_optimize, x)})
    these_settings.loc[new_vals.index, mt] = new_vals
    these_settings = these_settings.ix[settings.index, settings.columns]
    last_score = run_model((these_settings.ix[:, mt], "../Data/final_calibration/final/final_run.mat", obs_data.copy(), 'gene_objective'))

new_res_scores = [("Fe-", 0.434782138697),
                  ("N+", 0.798278269577),
                  ("O", 0.27726432923),
                  ("S+", 0.747253228853),
                  ("Ammonia Oxidation (oxygen)", 0.587060673503),
                  ("C1 Oxidation (Sum)", -0.730431444158),
                  ("Denitrification (Sum)", -1.2492263305),
                  ("Iron Reduction", -0.22593009868),
                  ("Methanogenesis", 1.0),
                  ("Sulfate Reduction + Sulfur Oxidation (Sum)", 0.596767942531)]

old_res_scores = [("Fe-", 0.493773599971),
                  ("N+", 0.79108783725),
                  ("O", 0.060216731052),
                  ("S+", 0.84738517233),
                  ("Ammonia Oxidation (oxygen)", 0.472693430034),
                  ("C1 Oxidation (Sum)", -1.05188841441),
                  ("Denitrification (Sum)", -1.01351297295),
                  ("Iron Reduction", -0.697225912725),
                  ("Methanogenesis", 1.0),
                  ("Sulfate Reduction + Sulfur Oxidation (Sum)", -0.498723281337)]

ts2 = these_settings.drop(["new_model_bool", "old_model_bool"], 1)
avgs = ['chem average', "rate average", "model average"]
score_df = pd.DataFrame(index=list(zip(*new_res_scores)[0])+avgs, columns=ts2.columns)

for idx in range(len(old_res_scores)):
    idx_str = old_res_scores[idx][0]
    score_df.ix[idx_str, 'old_model'] = old_res_scores[idx][1]
    score_df.ix[idx_str, 'new_model'] = new_res_scores[idx][1]

all_cols = list(zip(*new_res_scores)[0])
chem_cols, proc_cols = all_cols[:4], all_cols[4:]
mod_cols = ["new_model", "old_model"]
score_df.loc["chem average", mod_cols] = score_df.ix[chem_cols, mod_cols].mean()
score_df.loc["rate average", mod_cols] = score_df.ix[proc_cols, mod_cols].mean()
score_df.loc["model average", mod_cols] = score_df.ix[all_cols, mod_cols].mean()
ts3 = ts2.append(score_df)
ts3.to_csv('../Data/final_calibration/scores_settings.csv')