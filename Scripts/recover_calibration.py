from calibration_functions import cal_iteration, run_model
import cPickle as pickle
import sys, os
import pandas as pd
import numpy as np

results_loc = "../Data/calibration_results"
in_data_loc = "../Data/calibration_data"
settings = pd.read_csv(in_data_loc+"/calibration_run_settings.tsv", index_col=0, sep="\t")
n_samplings = [1600, 1600,  1600,  800,  800,  800,  800,  800,  400,  400,  400,  400,  100,  100,  100,  100, 100]
chem_df = pd.read_csv(in_data_loc+"/chem_data.tsv", sep="\t", index_col=[0,1], parse_dates=True)
gene_df = pd.read_csv(in_data_loc+"/gene_data.tsv", sep="\t", index_col=[0,1], parse_dates=True)
obs_data = chem_df.join(gene_df)

dev_trace, param_trace = {}, {}

for model_run_col in [0, 2]:
    mod_type = settings.columns[model_run_col]
    results_subdir = results_loc + "/" + mod_type
    bool_key = settings.columns[model_run_col+1]
    to_optimize = settings[settings[bool_key]].index
    iterations = settings[bool_key].sum() + 1
    save_files = ['save'+"{:02}.p".format(i) for i in range(1,iterations+1)]
    
    dev_trace[mod_type] = pd.DataFrame(index=to_optimize, columns=range(iterations))
    param_trace[mod_type] = pd.DataFrame(index=list(to_optimize)+["score"], columns=range(iterations+1))
    
    last_score = run_model((settings.ix[:, mod_type], results_subdir+"/baseline/baseline_run.mat", obs_data.copy(), 'gene_objective'))
    print "Initial score for {} is {:.2f}".format(mod_type, last_score)
    
    # load initial run score & parameter values 
    param_trace[mod_type].ix[to_optimize, 0] = settings.ix[to_optimize, mod_type]
    param_trace[mod_type].ix['score', 0] = last_score
    
    for idx in range(1,iterations+1):
        dump_f = os.path.join(results_subdir, save_files[idx-1])
        with open(dump_f, "rb" ) as p_fh:
            trial_result = pickle.load(p_fh)
        stds, settings, last_score, winner, converged_ = trial_result
        final_cal_f = os.path.join(results_subdir, mod_type, "best_runs", str(idx)+"_best.mat")
        run_model((settings.ix[:, mod_type], final_cal_f, obs_data.copy(), 'print_objective'))
        dev_trace[mod_type].ix[stds.keys(), idx] = pd.Series(stds)
        param_trace[mod_type].ix[to_optimize, idx] = settings.ix[to_optimize, mod_type]
        param_trace[mod_type].ix['score', idx] = last_score
        
    print "Completed calibration of {}".format(mod_type)
    

dump_f = os.path.join(results_loc, 'final_result.p')
dump_c = (param_trace, dev_trace)
with open(dump_f, "wb") as dump_h:
    pickle.dump( dump_c, dump_h )
