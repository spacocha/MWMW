#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 19:44:45 2017

@author: login
"""
import matplotlib.pyplot as plt
import subprocess32 as sp
import cPickle as pickle
import pandas as pd
import platform, shutil, sys, os
from scipy.interpolate import interp1d
from collections import OrderedDict
import numpy as np
from scipy.io import loadmat
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.feature_selection import f_regression
from sklearn.linear_model import LinearRegression

from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count

def cal_iteration(idx, res_loc, save_f, settings, bool_key, obs_data, n_trials, mod_type, stds, last_score, conv_, param_trace):
    # Calculate the number of threads possible
    thread_est = cpu_count()
    # Determine if it is initial run or a restart
    if idx == 1:
        print "Starting new opt with settings for {} using {} cores".format(mod_type, thread_est)
        init_bool = True
    else:
        init_bool = False

    # determine which variables are left to optimize
    not_optimized = settings[settings.ix[:, bool_key]].index.tolist()
    # create a directory for each model run output
    lake_dirs = [os.path.join(res_loc, "lake_"+str(i))  for i in range(1,n_trials+1)]
    # create a filename for the output
    result_str = ["run_{}.mat".format(i) for i in range(1,n_trials+1)]
    # create paths for the output
    result_paths = [os.path.join(j, i) for i, j in zip(result_str, lake_dirs)]
    # fill a df with param choices 
    parameter_df = sample_params(settings, init_bool, n_trials, not_optimized, mod_type, stds)
    
    print "Search space covereage:"
    for no_v in not_optimized:
        srch_spce = settings.ix[no_v, "upper_limit"] - settings.ix[no_v, "lower_limit"]
        samp_spce = parameter_df.ix[:, no_v].max() - parameter_df.ix[:, no_v].min()
        cov_pct = samp_spce / float(srch_spce)
        print "\t{}: {:.2%}".format(no_v, cov_pct)


    # run the model n_trials number of times if it hasn't converged
    if not conv_:
        score_vector = pd.Series(index=parameter_df.index, data=np.zeros((n_trials)))

        model_args = []
        for trial_n in range(1, n_trials+1):
            arg1 = parameter_df.ix[trial_n, :]
            arg2 = result_paths[trial_n-1]
            arg3 = obs_data.copy()
            arg4 = 'gene_objective'
            model_args.append((arg1, arg2, arg3, arg4))
        
        parallelism = True
        
        if not parallelism:
            score_list = []
            for d, p, o, s in model_args:
                print os.path.basename(p)
                score_list.append(run_model((d, p, o, s)))   
        else:
            pool = ThreadPool(thread_est)
            score_list = pool.map(run_model, model_args)
            pool.close()
            pool.join()
            
        score_vector[range(1, n_trials+1)] = score_list
    else:
        score_vector = pd.Series(index=parameter_df.index,
                                 data = np.ones((n_trials))*last_score)
    
    parameter_df_f = os.path.join(res_loc, "trial_data_{}.csv".format(idx))
    parameter_df['score'] = score_vector
    parameter_df.to_csv(parameter_df_f)
     
    # Step 4: Perform an F-test to identify the most sensitive parameter & fix it
    # Note: F-test is insensitive to range, but the model is most sensitive to the 
    # range of the input values
    x_df = parameter_df.ix[:, not_optimized]
    F_vals, p_vals = f_regression(x_df.values, score_vector.values)
    p_series = pd.Series(index=x_df.columns, data=p_vals)

    # check to see if no parameters make a difference
    print "{} variables were insensitive".format(p_series.isnull().sum())
    if p_series.isnull().sum() == len(p_series):
        p_series[0] = 0.05
        conv_ = True
    
    winner = p_series[p_series == p_series.min()].index[0]

    # Step 5: Take the best value for that parameter and those not requiring optimization
    # Ensure this run beats the last run
    param_hist = param_trace[mod_type]
    best_score = param_hist.ix["score", :].max()

    if best_score > parameter_df.score.max():
        best_run = param_hist.T[param_hist.T.score == best_score]
        opt_value = best_run.ix[:, winner].values[-1]
        this_score = best_score
        print "Previous run had a better value, using that instead"
    else:
        opt_value = best_val(parameter_df, winner)
        this_score = parameter_df.score.max()

    improvement = (this_score-last_score)
    print "{:.2f}, a change of {:.1%}".format(this_score, improvement)

    if not conv_:
        print "{} is the most sensitive variable".format(winner)
    else:
        print "Converged at maximum, no further optimization necessary"
    
    if idx != 1:
        print "Real Iteration"
        not_optimized.remove(winner)
        settings.ix[winner, mod_type] = opt_value
        settings.ix[winner, bool_key] = False
    else:
        print "All parameters will be regressed (burn in iteration)"
    
    # Step 6: Perform regression on the remaining optimizable params 
    std_deviations = {}
    for t_o in not_optimized:
        lr = LinearRegression()
        y = parameter_df.ix[:, 'score'].values.reshape(-1, 1)
        x = parameter_df.ix[:, t_o].values.reshape(-1, 1)
        lr.fit(y,x)
        best_guess = lr.predict(np.array([[y.max()]]))
        std_deviations[t_o] = abs(best_guess - settings.ix[t_o, mod_type])[0][0]*2
        settings.ix[t_o, mod_type] = best_guess[0][0]
    
    dump_f = os.path.join(res_loc, save_f)
    dump_c = (std_deviations, settings, this_score, winner, conv_)
    with open(dump_f, "wb") as dump_h:
        pickle.dump( dump_c, dump_h )
    return dump_c

def combine_similar_procs(mod_df, obs_df):
    mod_df2, obs_df2 = mod_df.copy(), obs_df.copy()
    
    mo_cols = ['C1 Oxidation (nitrate)', 
               'C1 Oxidation (oxygen)',
               'C1 Oxidation (sulfate)']
    
    mod_df2.loc[:, 'C1 Oxidation (Sum)'] = mod_df2.loc[:, mo_cols].sum(1)
    
    srso_cols = ['Sulfate Reduction + Sulfur Oxidation (Sum)',
                 'Sulfur Oxidation (nitrate)', 
                 'Sulfur Oxidation (oxygen)',
                 'Sulfate Reduction']
                 
    mod_df2.loc[:, srso_cols[0]] = mod_df2.loc[:, srso_cols[1:]].sum(1)
    
    dn_cols = ['Denitrification (Sum)', 
               'Denitrification',
               'Iron Oxidation (nitrate)', 
               'C1 Oxidation (nitrate)',
               'Sulfur Oxidation (nitrate)']
    
    mod_df2.loc[:, dn_cols[0]] = mod_df2.loc[:, dn_cols[1:]].sum(1)
    
    # drop all model data that is not scorable
    grp1 = set(obs_df2.columns)
    grp2 = set(mod_df2.columns)
    odd_cols = list(grp1.symmetric_difference(grp2))
    matched_model = mod_df2.drop(odd_cols, axis=1)
    return matched_model
    
def recreate_opt_sequence(path, opt_opt):
    if opt_opt == 'optimum':
        opt_type = 'new_model'
    else:
        opt_type = 'old_model'
        
    starting_point = set(load_optimization_var_list(opt_type))
    pickle_files = sorted([i for i in os.listdir(path) if i.endswith(".p")])
    pickle_paths = [os.path.join(path, i) for i in pickle_files]
    var_trace = []
    for pp in pickle_paths:
        with open( pp, "rb" ) as pickle_h:
            pickle_pack = pickle.load(pickle_h)
        to_optimize = pickle_pack[1]
        dropped = starting_point.symmetric_difference(set(to_optimize))
        assert len(dropped) == 1
        var_trace.append(dropped)
        starting_point.remove(list(dropped)[0])
    return var_trace

def importratesandconcs_mod(path_, type_str=None):
    """
    1. Create a date index for the new dataframes
    2. Create new dictionaries to hold the new dataframes
    3. Unload each DF one at a time 
    4. Interpolate each depth vector along new axis
    5. Load into new numpy array
    6. Assign date index & numpy array to new dataframe object
    7. Reload new dataframe into new dictionary, accessible by name string
    8. Return newly minted dictionaries
    """
    conc_idxs = OrderedDict()
    rate_idxs = OrderedDict()
    
    conc_idxs[0] = "O"
    conc_idxs[1] = "C"
    conc_idxs[2] = "N+"
    conc_idxs[3] = "N-"
    conc_idxs[4] = "S+"
    conc_idxs[5] = "S-"
    conc_idxs[6] = "Fe+"
    conc_idxs[7] = "Fe-"
    conc_idxs[8] = "CH4"
    conc_idxs[9] = "Null"
    rate_idxs[0] = "Iron Oxidation (oxygen)"
    rate_idxs[1] = "Ammonia Oxidation (oxygen)"
    rate_idxs[2] = "Sulfur Oxidation (oxygen)"
    rate_idxs[3] = "Iron Oxidation (nitrate)"
    rate_idxs[4] = "Sulfur Oxidation (nitrate)"
    rate_idxs[5] = "C1 Oxidation (oxygen)"
    rate_idxs[6] = "C1 Oxidation (nitrate)"
    rate_idxs[7] = "C1 Oxidation (sulfate)"
    rate_idxs[8] = "Aerobic Heterotrophy"
    rate_idxs[9] = "Denitrification"
    rate_idxs[10] = "Iron Reduction"
    rate_idxs[11] = "Sulfate Reduction"
    rate_idxs[12] = "Methanogenesis"
    
    if os.path.exists(path_):
        mat = loadmat(path_)
        concs_ = mat['concs_history']
        rates_ = mat['rates_history']
    else:
        concs_ = np.zeros((100, 17, 10))
        rates_ = np.zeros((100, 17, 13))
        
    if concs_.shape[0] != 100:
        concs_ = np.zeros((100, 17, 10))
        rates_ = np.zeros((100, 17, 13))
        
    conc_idxs_inv = {v: k for k, v in conc_idxs.iteritems()}
    rate_idxs_inv = {v: k for k, v in rate_idxs.iteritems()}
    
    mat_dict = {}
    inputs =  [concs_, rates_]
    translators = [conc_idxs_inv, rate_idxs_inv]
                
    for i_arr, t_dict in zip(inputs, translators):
        for name, idx_z in t_dict.items():
            mat_dict[name] = pd.DataFrame(data=i_arr[:, :, idx_z].T, 
                                          columns=range(0,100),
                                          index=range(6,23))
    n_days = 146
    start_date, end_date = '03/23/2013', '08/15/2013'
    dr = pd.date_range(start_date, end_date)
    assert len(dr) == n_days
    
    new_mat_dict = {}
    for a_spec in mat_dict.keys():
        this_df = mat_dict[a_spec]
        depths, n_slices = this_df.shape
        assert n_slices < n_days
        idx = np.arange(n_slices)
        new_interval = max(idx) / float(n_days)
        new_columns = np.arange(idx.min(), idx.max(), new_interval)
        new_df_data = np.zeros((depths, len(new_columns)))
        for depth in xrange(depths):
            a_vector = this_df.ix[depth+6, :].values
            f2 = interp1d(idx, a_vector, kind='cubic')
            new_df_data[depth, :] = f2(new_columns)
        new_df = pd.DataFrame(data=new_df_data.T, columns=np.arange(6,6+depths),
                              index=dr)
        if type_str and type_str == 'full df':
            new_mat_dict[a_spec] = new_df.T.unstack()
        else:
            new_mat_dict[a_spec] = new_df.T
    
    all_cols = sorted(new_mat_dict.keys())
    if type_str and type_str == 'full df':
        full_idx = new_mat_dict[all_cols[0]].index
        full_df = pd.DataFrame(index=full_idx, columns=all_cols)
        
        for name in all_cols:
            full_df.ix[:, name] = new_mat_dict[name] 
        
        total_proc = full_df.ix[:, rate_idxs_inv.keys()].sum(axis=1)
        full_df_norm = full_df.copy()
        for key in rate_idxs_inv.keys():
            full_df_norm[key] = full_df[key].divide(total_proc, axis=0)
        
        return full_df
    else:
        return new_mat_dict 
    
def copyDirectory(src, dest):
    try:
        shutil.copytree(src, dest)
    except shutil.Error as e:
        print "error: %s" % e
        raise OSError('doink')
    except OSError as e:
        print "error: %s" % e
        raise OSError('doink')

def turn_mat_into_single_df(mat):
    columns = ['nitrate', 'sulfate', 'iron(II)', 'oxygen' ]
    z_idxs = [2, 4, 7, 0 ]
    mat_dict = {}
    for idx, col in zip(z_idxs, columns):
        mat_dict[col] = pd.DataFrame(data=mat[:,:,idx].T, columns=range(0,100),
                                     index=range(6,23))

    n_days = 146
    start_date, end_date = '03/23/2013', '08/15/2013'
    dr = pd.date_range(start_date, end_date)
    assert len(dr) == n_days

    new_mat_dict = {}
    for a_spec in mat_dict.keys():
        this_df = mat_dict[a_spec]
        depths, n_slices = this_df.shape
        assert n_slices < n_days
        idx = np.arange(n_slices)
        new_interval = max(idx) / float(n_days)
        new_columns = np.arange(idx.min(), idx.max(), new_interval)
        new_df_data = np.zeros((depths, len(new_columns)))
        for depth in xrange(depths):
            a_vector = this_df.ix[depth+6, :].values
            f2 = interp1d(idx, a_vector, kind='cubic')
            new_df_data[depth, :] = f2(new_columns)
        new_df = pd.DataFrame(data=new_df_data.T, columns=np.arange(6,6+depths),
                              index=dr)
        new_mat_dict[a_spec] = new_df.T.unstack()
    
    all_cols = sorted(new_mat_dict.keys())
    full_idx = new_mat_dict[all_cols[0]].index

    full_df = pd.DataFrame(index=full_idx, columns=all_cols)
    
    for name in all_cols:
        full_df.ix[:, name] = new_mat_dict[name]

    return full_df
    
def standard_scale_df(df):
    sklearn_ss = MinMaxScaler()
    
    if type(df) == type(pd.DataFrame()):
        std_data = sklearn_ss.fit_transform(df.values)
        return pd.DataFrame(data=std_data, index=df.index, columns=df.columns)
        
    elif type(df) == type(pd.Series()):
        std_data = sklearn_ss.fit_transform(df.values.reshape(-1, 1))
        return pd.Series(data=std_data.flatten(), index=df.index)
    else:
        sys.exit("Unrecognized data type for scaling")
       
def min_max_scale_df(df):
    sklearn_ss = MinMaxScaler()    
    
    if type(df) == type(pd.DataFrame()):
        std_data = sklearn_ss.fit_transform(df.values)
        return pd.DataFrame(data=std_data, index=df.index, columns=df.columns)
        
    elif type(df) == type(pd.Series()):
        std_data = sklearn_ss.fit_transform(df.values.reshape(-1, 1))
        return pd.Series(data=std_data.flatten(), index=df.index)
    else:
        sys.exit("Unrecognized data type for scaling")
        
        
def score_results(obs_df_, data_df_, score_type):
    # temporal subsetting
    obs_df, data_df = obs_df_.copy(), data_df_.copy()
    bool1 = data_df.index.isin(obs_df.index)
    bool2 = obs_df.index.isin(data_df.index)
    assert (bool2.sum() > 0) and (bool1.sum() > 0)
    sub_data_df = data_df[bool1]
    sub_obs_df = obs_df[bool2]
    sub_data_df = combine_similar_procs(sub_data_df, sub_obs_df)
    assert sub_obs_df.shape == sub_data_df.shape

    # drop all columns that aren't related to the specific objective
    if score_type == 'gene_objective':
        pass
    elif score_type == 'print_objective':
        pass
    elif score_type == 'conc_objective':
        to_keep = set(['Fe-', 'N+', 'S+', 'O'])
        all_cols = set(sub_obs_df.columns)
        to_drop = to_keep.symmetric_difference(all_cols)
        sub_data_df.drop(to_drop, axis=1, inplace=True)
        sub_obs_df.drop(to_drop, axis=1, inplace=True)
        
    r2_array = np.zeros((len(sub_obs_df.columns),))
    
    for idx, col_ in enumerate(sub_obs_df.columns):
        obs_vec = sub_obs_df.ix[:, col_].dropna()
        mod_vec = sub_data_df.ix[:, col_].dropna()
        
        bool_11 = obs_vec.index.isin(mod_vec.index)
        bool_22 = mod_vec.index.isin(obs_vec.index)
        
        obs_vec_nn = obs_vec[bool_11]
        mod_vec_nn = mod_vec[bool_22]
        
        obs_vec_std = standard_scale_df(obs_vec_nn)
        data_vec_std = standard_scale_df(mod_vec_nn)

        obs_vals = obs_vec_std.values        
        model_vals = data_vec_std.values
        r2_array[idx] = abs(model_vals - obs_vals).mean()*-1
        #r2_array[idx] = r2_score(model_vals, obs_vals)
        #print "{}: {}".format(col_, r2_score(model_vals, obs_vals))
        
    #r2_array[r2_array < -1.] = -1.

    if score_type == 'print_objective':
        return None
    else:
        print r2_array.sum()
        return r2_array.sum()

def run_model(arg_tuple):
    subdf, out_f, obs_data, score_type = arg_tuple
    print os.path.basename(out_f).split(".")[0]
    if platform.system() == 'Linux':
        run_cmd = ""
    else:
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
        out_h.write(os.path.abspath(out_f))
    
    # run the model
    p = sp.Popen(run_cmd, cwd=model_loc, shell=True, stderr=sp.PIPE, stdout=sp.PIPE)
    try:
        stdout, stderr = p.communicate(timeout=int(600.*5))
        results_df = importratesandconcs_mod(out_f, 'full df')
    except sp.TimeoutExpired:
        p.kill()
        outs, errs = p.communicate()
        results_df = obs_data*0.0
        print "Timeout Exception Triggered, Process Killed"

    # pull results & return them to memory

    r2 = score_results(obs_data, results_df, score_type)
    shutil.rmtree(model_loc)
    return r2
    
def sample_params(settings, init_bool, n_trials, to_optimize, mod_type, stds):
    # create a matrix with trials in rows & parameter names in columns
    parameter_mat = np.zeros((n_trials, len(settings.index)))
    # use this in a dataframe 
    parameter_df = pd.DataFrame(index=np.arange(1,n_trials+1), columns=settings.index, data=parameter_mat)
    for key in settings.index:
        # start at the defaults
        this_mean = settings.ix[key, mod_type]
        # only randomize for paramters that need optimization
        if key in to_optimize:
            if init_bool:
                low_, high_ = settings.ix[key, ['lower_limit', 'upper_limit']].values
                this_sample = abs(np.linspace(low_, high_, n_trials))
            else:
                this_std = stds[key]
                this_sample = abs(np.random.normal(this_mean, abs(this_std)+np.finfo(float).eps, 
                                                  (n_trials,1)))
        else:
            this_sample = np.ones((n_trials,1))*this_mean
            
        np.random.shuffle(this_sample)
        parameter_df.ix[:, key] = this_sample
    return parameter_df

def best_val(param_df, variable):
    best_score = param_df.score.max()
    best_bool = param_df.score == best_score
    best_idx = param_df[best_bool].index
    return param_df.ix[best_idx, variable].values[0]

def apply_conc_multiplier(param_subdf, f_name):
    columns = ["S-", "C", "N-", "CH4"]
    col_no = [5, 1, 3, 8]
    multiplier = [param_subdf[i] for i in columns]
    conc0 = pd.read_csv(f_name, header=None)
    for c_num, mult in zip(col_no, multiplier):
        conc0.ix[:, c_num] *= mult
    conc0.to_csv(f_name, float_format="%g", header=False, index=False)
    return None

