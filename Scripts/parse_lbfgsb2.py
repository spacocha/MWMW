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
for x, mi in zip(mod_res, mod_types):
    mt = settings.columns[mi]
    mb = settings.columns[mi+1]
    to_optimize = list(settings.ix[settings.ix[:, mb], mt].index)
    assert len(to_optimize) == len(x)
    new_vals = pd.Series({i:j for i, j in zip(to_optimize, x)})
    these_settings.loc[new_vals.index, mt] = new_vals
    these_settings = these_settings.ix[settings.index, settings.columns]
    score_list.append(run_model((these_settings.ix[:, mt], "../Data/final_calibration/finala/finala_run.mat", obs_data.copy(), 'return_scores')))


new_res_scores = [("Fe-", 0.415706976913),
                  ("N+", 0.818041113592),
                  ("O", 0.275122702656),
                  ("S+", 0.773227738483),
                  ("Ammonia Oxidation (oxygen)", 0.647205634013),
                  ("C1 Oxidation (Sum)", -0.707954645044),
                  ("Denitrification (Sum)", -1.85644880163),
                  ("Iron Reduction", -0.200903221249),
                  ("Methanogenesis", 1.0),
                  ("Sulfate Reduction + Sulfur Oxidation (Sum)", 0.638621449069)]

old_res_scores = [("Fe-", 0.491221472257),
                  ("N+", 0.769363903053),
                  ("O", 0.0579604552629),
                  ("S+", 0.845973018746),
                  ("Ammonia Oxidation (oxygen)", 0.733623630935),
                  ("C1 Oxidation (Sum)", -1.02805034552),
                  ("Denitrification (Sum)", -0.573224365163),
                  ("Iron Reduction", -0.417504662629),
                  ("Methanogenesis", 1.0),
                  ("Sulfate Reduction + Sulfur Oxidation (Sum)", -0.485158926779)]

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
ts3.to_csv('../Data/final_calibration/scores_settings2.csv')
