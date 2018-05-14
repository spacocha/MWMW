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
                      'eps': 1e-05,
                      'maxiter': 100,
                      'ftol': 2.220446049250313e-09,
                      'maxcor': 10,
                      'maxfun': 100})

print z