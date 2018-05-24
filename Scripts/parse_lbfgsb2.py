import numpy as np
import pandas as pd
from scipy.optimize import minimize
from calibration_functions import run_model
import os

old_copied_res = np.array([  1.67691373e+02,   5.92937669e-02,   9.90090090e+04,
         1.55070270e+04,   6.61210182e-02,   1.47659072e-01,
         7.77478986e+01,   2.74743708e+00,   2.22209918e+00,
         5.08921069e-01,   5.98918069e+00,   2.83666731e-01,
         8.02864567e+00])

new_copied_res = np.array([  1.70395880e+02,   9.32218103e-02,   5.24324324e+04,
         1.65535135e+04,   9.65393666e-02,   6.39060123e-02,
         8.68327122e-03,   1.30895343e+02,   6.50887270e-01,
         1.68072590e-01,   3.41371397e+01,   2.00946863e+00,
         1.55067447e+00,   3.95156918e+00,   1.63120493e-01,
         3.00092414e+00])

in_data_loc = "../Data/calibration_data"
settings = pd.read_csv(in_data_loc+"/calibration_run_settings.tsv", index_col=0, sep="\t")
chem_df = pd.read_csv(in_data_loc+"/chem_data.tsv", sep="\t", index_col=[0,1], parse_dates=True)
gene_df = pd.read_csv(in_data_loc+"/gene_data.tsv", sep="\t", index_col=[0,1], parse_dates=True)
obs_data = chem_df.join(gene_df)

mod_types = [0, 2]
mod_res = [new_copied_res, old_copied_res]

these_settings = settings.copy()
score_list = []

#for x, mi in zip(mod_res, mod_types):
x = mod_res[0]
mi = mod_types[0]
mt = settings.columns[mi]
mb = settings.columns[mi+1]
to_optimize = list(settings.ix[settings.ix[:, mb], mt].index)
assert len(to_optimize) == len(x)
new_vals = pd.Series({i:j for i, j in zip(to_optimize, x)})
these_settings.loc[new_vals.index, mt] = new_vals
these_settings = these_settings.ix[settings.index, settings.columns]

old_settings = pd.read_csv(in_data_loc+"/old_params.csv", index_col=0)
score_list.append(run_model((old_settings.ix[:, "Previous model parameters"],
                             "../Data/final_calibration/finalb/finalb_run.mat", 
                             obs_data.copy(), 
                             'return_scores')))

score_list.append(run_model((these_settings.ix[:, mt], "../Data/final_calibration/finala/finala_run.mat", obs_data.copy(), 'return_scores')))



all_scores = pd.concat(score_list, 1)
all_scores.columns = ['pub vals', 'opt vals']
all_scores.to_csv("../Data/final_calibration/final_scores_REAL.tsv", sep="\t", index_label="Process")
params_used = pd.concat((old_settings.ix[:, "Previous model parameters"], these_settings.ix[:, mt]), 1)
params_used.to_csv("../Data/final_calibration/final_param_values_REAL.tsv", sep="\t", index_label="Params")