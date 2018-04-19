#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 12:09:17 2017

@author: login
"""
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import os
from calibration_functions import copyDirectory, apply_conc_multiplier
from calibration_functions import score_results, importratesandconcs_mod
from calibration_functions import load_param_dict, fill_param_dict, min_max_scale_df
from calibration_functions import recreate_opt_sequence, combine_similar_procs
from observed_data import load_chem_data, load_gene_data
import subprocess as sp

# We can test against the published values or an optimized version thereupon

#test = 'original_baseline'
#test = 'optimized_baseline'
test = 'chem_only_baseline'
new_model_dir = 'gene_test_3'
labels = ["optimum", "baseline"]
params = load_param_dict('default')

# If we are using the published values, we can skip searching a second
# optimization folder and just pull out the parameter values we need
if test == 'optimized_baseline':
    old_model_dir = 'gene_test_2'   
    mod_dirs = [new_model_dir, old_model_dir]
elif test == 'chem_only_baseline':
    old_model_dir = 'n_13600'
    mod_dirs = [new_model_dir, old_model_dir]
    original_df = fill_param_dict(params, True, 1, 0, [], None)
    original_srs = original_df.ix[1, :]

# If we are scoring against an optimized version of the old model, this loop runs
# twice
for dir_fn, mod_type in zip(mod_dirs, labels):
    data_dir = os.path.join(os.getcwd(), dir_fn)
#    var_trace = recreate_opt_sequence(data_dir, mod_type)
    print "\nVarible Importance ({})".format(mod_type)
    print "==================================="
#    for _idx_, var in enumerate(var_trace):
#        print "\t{}. {}".format(_idx_+1, list(var)[0])
        
    df_paths = [os.path.join(data_dir, i) for i in os.listdir(data_dir) if i.endswith(".csv")]
    df_list = {i:pd.read_csv(j, index_col=0) for i, j in enumerate(sorted(df_paths))}
    
    n_total_trials = 0
    col_n = df_list[0].shape[1]
    
    for v in df_list.values():
        assert v.shape[1] == col_n    
        n_total_trials += v.shape[0]
    
    master_table = df_list[0].copy()
    
    
    for k in range(1,max(df_list.keys())+1):
        master_table = pd.concat([master_table, df_list[k]], ignore_index=True)
    
    print mod_type, master_table.shape
    
    from sklearn.linear_model import RidgeCV
    from sklearn.preprocessing import StandardScaler
    
    sklearn_ss = StandardScaler()
    
    master_mat = sklearn_ss.fit_transform(master_table.values)
    master_table_std = pd.DataFrame(index=master_table.index,
                                    columns=master_table.columns,
                                    data = master_mat)
    
    y_std = master_table_std.score.values
    x_std = master_table_std.drop(['score'], 1).values
    
    
    ridge_model = RidgeCV(alphas=(0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 50.),
                          normalize=False)
    ridge_model.fit(x_std, y_std)
    
    var_import = pd.Series(index=master_table.columns[:-1], data=ridge_model.coef_)
    
    
    print "\nVarible:Score Relationship ({})".format(mod_type)
    print "==================================="
    print var_import
    print "\nOptimal Parameter Set({})".format(mod_type)
    print "==================================="

    if mod_type == 'optimum':
        bool1 = master_table.ma_op_s_n_rate_const > 0
        bool2 = master_table.ix[:, 'N-'] > 0
        bool3 = master_table.ix[:, 'S-'] > 0
        bool4 = master_table.oxygen_source > 0
        bool5 = master_table.carbon_precip > 0
        bool6 = master_table.ma_op_ch4_s_rate_const > 0
        bool7 = master_table.ma_op_ch4_n_rate_const > 0
        bool8 = master_table.ma_op_o_s_rate_const > 0 
        bool9 = master_table.methane_source > 0
    
        bool_list = [bool1, bool2, bool3, bool4, bool5, bool6, bool7, bool8, bool9]
        super_bool1 = bool1 & bool2 & bool3 & bool4 & bool5 & bool6
        super_bool = super_bool1 & bool7 & bool8 & bool9
        obedience_check = (super_bool.sum()/float(len(super_bool)))
        if obedience_check == 1.:
            print "All trials obeyed rules of new model"
            best_rows = master_table
        else:
            print "{:.1%} of trials obeyed rules of new model".format(obedience_check)
            best_rows = master_table[super_bool]
        best_row_idx = best_rows[best_rows.score == best_rows.score.max()].index
        best_row = master_table.ix[best_row_idx[0], master_table.columns[:-1]]
        local_high = best_row
    elif mod_type == 'baseline':
        best_row_bool = master_table.score == master_table.score.max()
        best_row_idx = master_table[best_row_bool].index
        baseline_srs = master_table.ix[best_row_idx[0], master_table.columns[:-1]]
        local_high = baseline_srs
        
    for param in local_high.index:
        opt_val = local_high[param]
        old_val = params[param]
        opt_old_ratio = opt_val/float(old_val) - 1
        print "{}: {:g} changed {:.3%}".format(param, opt_val, opt_old_ratio)

best_lakes_dir = os.path.join(os.getcwd(), 'best_lakes_4')
os.mkdir(best_lakes_dir)

base_mat_f = os.path.join(best_lakes_dir, 'original_model', 'original_model.mat')
chem_only_f = os.path.join(best_lakes_dir, 'chem_only', 'chem_only.mat')
opt_mat_f = os.path.join(best_lakes_dir, 'chem_gene', 'chem_gene_2.mat')

df_pair = [best_row, baseline_srs, original_srs]
model_dirs = [opt_mat_f, chem_only_f, base_mat_f]
labels = ["optimum", "chem_only", "baseline"]
parameter_sets = pd.concat(df_pair, axis=1, verify_integrity=True)
output_labels = ["gene_chem_optimum", "chem_only_optimum", "previous_baseline"]
parameter_sets.columns = output_labels
parameter_sets.to_csv(os.path.join(best_lakes_dir, "parameter_sets.tsv"), sep="\t")

for subdf, out_f, label in zip(df_pair, model_dirs, labels):
    print "\nFit Metrics ({})".format(label)
    print "==================================="
    
    obs_df, gene_df = load_chem_data(), load_gene_data()
    obs_data = obs_df.join(gene_df)
    run_cmd = "/Applications/MATLAB_R2016b.app/bin/"

    # tell subprocess where and what to execute
    model_loc = os.path.dirname(out_f)
    source_dir = os.path.join(os.getcwd(), 'lake_model')
    copyDirectory(source_dir, model_loc)

    input_args_loc = os.path.join(model_loc, 'baseline.txt')
    # write out the parameter set into the right location
    subdf.T.to_csv(input_args_loc, header=False, index=False, float_format='%g')
    init_val_f = os.path.join(model_loc, "concs0.txt")
    apply_conc_multiplier(subdf, init_val_f)
    #matlab -nodisplay -nojvm -nosplash -nodesktop -r 'calibration_kaw; exit'
    run_cmd = run_cmd +"matlab -nodisplay -nojvm -nosplash -nodesktop " 
    run_cmd = run_cmd +"-r 'calibration_kaw; exit'"

    # what is the output file name ? 
    output_loc = os.path.join(model_loc, 'outFile.txt')
    
    with open(output_loc, 'w') as out_h:
        out_h.write(out_f)
    # run the model 
    p = sp.Popen(run_cmd, cwd=model_loc, shell=True, 
                 stderr=sp.PIPE, stdout=sp.PIPE)
    stdout, stderr = p.communicate()
    
    # pull results & return them to memory
    
    results_df = importratesandconcs_mod(out_f, 'full df')
    
#    r2_chem = score_results(obs_data, results_df, 'conc_objective')
    r2_genes = score_results(obs_data, results_df, 'gene_objective')
        
#    print "{}\n\tR2 (genes): {}\n\tR2 (chem): {}".format(label, r2_genes, r2_chem)


new_model_df = importratesandconcs_mod(opt_mat_f, 'full df')

# This section subsets the data for plotting the chemistry time series fit plot
obs_chem_df = load_chem_data()

new_chem_df = new_model_df.ix[:, list(obs_chem_df.columns)]

new_o_cols = [i+" (obs)" for i in obs_chem_df.columns]
obs_chem_df.columns = new_o_cols
obs_chem_t = obs_chem_df.T

new_m_cols = [i+" (new)" for i in new_chem_df.columns]
new_chem_df.columns = new_m_cols
new_chem_sub = new_chem_df[new_chem_df.index.isin(obs_chem_df.index)].T
                           
all_chem_data = pd.concat([new_chem_sub, obs_chem_t], 
                     verify_integrity=True)

plot_fn = 'rates_n_concentrations.pdf'
plot_path = os.path.join(best_lakes_dir, plot_fn)
pp = PdfPages(plot_path)
all_data_new_idx = all_chem_data.copy()
all_data_new_idx.columns = np.arange(85)
styles = ['c-', 'b-', 'r-', 'g-','c--', 'b--', 'r--', 'g--']
all_data_new_idx.T.plot(figsize=(12, 14), linewidth=3., style=styles)
plt.savefig(pp, format='pdf')
plt.clf()

# This section subsets the data for plotting the individual time point gene set fits
obs_gene_df = load_gene_data()
new_gene_df = combine_similar_procs(new_model_df, obs_gene_df)

new_gene_df = new_gene_df.ix[:, list(obs_gene_df.columns)]


new_o_cols = [i+" (obs)" for i in obs_gene_df.columns]
obs_gene_df.columns = new_o_cols
obs_gene_t = obs_gene_df.T

new_m_cols = [i+" (new)" for i in new_gene_df.columns]
new_gene_df.columns = new_m_cols
new_gene_sub = new_gene_df[new_gene_df.index.isin(obs_gene_df.index)].T
                           
all_gene_data = pd.concat([new_gene_sub, obs_gene_t], 
                           verify_integrity=True)

all_gene_data_idx = all_gene_data.copy()
depths_ = np.array([1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 20, 21, 22])
all_gene_data_idx.columns = depths_
title_list = list(load_gene_data().columns)

for title_ in title_list:
    idx_subset = [i for i in list(all_gene_data_idx.index) if title_ in i]
    assert len(idx_subset) == 2
    data_subset = all_gene_data_idx.ix[idx_subset, :]
    data_subset_noNan = data_subset.drop([1,3,5], axis=1).T
    data_subset_std = min_max_scale_df(data_subset_noNan).T
    styles = ['c-', 'b-', 'r-']
    data_subset_std.T.plot(figsize=(12, 14), linewidth=3., style=styles, 
                           title=title_)
    plt.savefig(pp, format='pdf')
    plt.clf()

results_dict = importratesandconcs_mod(opt_mat_f, None)
results_df = importratesandconcs_mod(opt_mat_f, 'full df')
    
# Here we can make a pdf of all the dataframes after stretching them 
del results_dict['Null']

plt.rcParams['figure.figsize']=(16,12)
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14

# start a loop here
for key, result_df in results_dict.items():
    top_lim = result_df.ix[:, 10:].max().max()
    result_df[result_df < 0.001] = 0
    result_df[result_df > top_lim] = top_lim
    plt.figure(1)
    plt.xticks(rotation='vertical')
    new_cols = [i.strftime('%m.%d') for i in result_df.columns]
    result_df.columns = new_cols
    ax = sns.heatmap(result_df, cmap="YlGnBu")
    plt.title(key)
    labels = ax.get_xticklabels()
    new_labels = []
    for idx, i in enumerate(labels):
        if (idx % 3) == 0:
            new_labels.append(i)
        else:
            new_labels.append("")
            
    ax.set_xticklabels(new_labels)
    plt.tight_layout()
    plt.savefig(pp, format='pdf')
    plt.clf()


pp.close()




