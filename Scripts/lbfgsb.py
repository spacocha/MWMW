import numpy as np
import pandas as pd
from scipy.optimize import minimize
from calibration_functions import run_model

in_data_loc = "../Data/calibration_data"
chem_df = pd.read_csv(in_data_loc+"/chem_data.tsv", sep="\t", index_col=[0,1], parse_dates=True)
gene_df = pd.read_csv(in_data_loc+"/gene_data.tsv", sep="\t", index_col=[0,1], parse_dates=True)
obs_data = chem_df.join(gene_df)

settings = pd.read_csv(in_data_loc+"/calibration_run_settings.tsv", index_col=0, sep="\t")
x0_ = settings.ix[settings.new_model_bool, "new_model"].tolist()
bounds_ = settings.ix[settings.new_model_bool, ["lower_limit", "upper_limit"]].apply(tuple, axis=1).tolist()
to_optimize = list(settings.ix[settings.new_model_bool, "new_model"].index)

def random_function(x):
    new_vals = pd.Series({i:j for i, j in zip(to_optimize, x)})
    these_settings = settings.copy()
    these_settings.loc[new_vals.index, 'new_model'] = new_vals
    newvs = these_settings.ix[new_vals.index, 'new_model'].values
    oldvs = settings.ix[new_vals.index, 'new_model'].values
    print "{:2%}".format(np.divide(abs(newvs - oldvs), oldvs).sum())
    last_score = run_model((these_settings.ix[:, 'new_model'], 
                            "../Data/final_calibration/old_final/old_run.mat", 
                            obs_data.copy(), 
                            'gene_objective'))
    return last_score*-1

z = minimize(random_function,
             x0_,
             method='L-BFGS-B',
             jac=False,
             bounds=bounds_,
             options={'disp': True,
                      'maxls': 20,
                      'iprint': 101,
                      'gtol': 1e-05,
                      'eps': 5e-05,
                      'maxiter': 200,
                      'ftol': 1e-05,
                      'maxcor': 10,
                      'maxfun': 1500})

print z

# 1e-5 (100 trials)
#array([  9.65000000e+01,   8.43558520e-02,   5.50000000e+04,
#         6.60000000e+03,   4.00000000e-01,   3.63691447e-01,
#         2.85504086e-01,   5.00000000e+01,   2.18617308e+00,
#         4.36009330e-01,   5.00000000e+01,   2.61605982e+00,
#         1.22410136e+00,   1.76502268e+00,   6.76877484e-01,
#         3.27699418e+00]) = 2.15
# 1e-6 (1500 trials)
#array([  9.65000000e+01,   3.24702142e-02,   5.50000000e+04,
#         6.60000000e+03,   3.97653239e-01,   2.34708464e-01,
#         3.88898357e-03,   5.00000000e+01,   1.00000000e+00,
#         5.00000000e-01,   5.00000000e+01,   1.50000000e+00,
#         1.00000000e+00,   4.58860532e-05,   1.65434361e-01,
#         4.54584061e-05] = 2.27
