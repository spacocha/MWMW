#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 19:44:45 2017

@author: login
"""
from scipy.interpolate import interp1d
from collections import OrderedDict
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
from scipy.io import loadmat
import shutil, sys, os
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import cPickle as pickle

def combine_similar_procs(mod_df, obs_df):
    mod_df2, obs_df2 = mod_df.copy(), obs_df.copy()
    
    mo_cols = ['Methane Oxidation (nitrate)', 
               'Methane Oxidation (oxygen)',
               'Methane Oxidation (sulfate)']
    
    mod_df2.loc[:, 'Methane Oxidation (Sum)'] = mod_df2.loc[:, mo_cols[0]] + \
                                                mod_df2.loc[:, mo_cols[1]] + \
                                                mod_df2.loc[:, mo_cols[2]] 
    
    srso_cols = ['Sulfate Reduction + Sulfur Oxidation (Sum)',
                 'Sulfur Oxidation (nitrate)', 'Sulfur Oxidation (oxygen)',
                 'Sulfate Reduction']
                 
    mod_df2.loc[:, srso_cols[0]] = mod_df2.loc[:, srso_cols[1]] + \
                                   mod_df2.loc[:, srso_cols[2]] + \
                                   mod_df2.loc[:, srso_cols[3]]
    
    dn_cols = ['Denitrification (Sum)', 'Denitrification',
               'Iron Oxidation (nitrate)', 'Methane Oxidation (nitrate)',
               'Sulfur Oxidation (nitrate)']
    
    mod_df2.loc[:, dn_cols[0]] = mod_df2.loc[:, dn_cols[1]] + \
                                 mod_df2.loc[:, dn_cols[2]] + \
                                 mod_df2.loc[:, dn_cols[3]] + \
                                 mod_df2.loc[:, dn_cols[4]]
    
    # drop all model data that is not scorable
    grp1 = set(obs_df2.columns)
    grp2 = set(mod_df2.columns)                 
    odd_cols = list(grp1.symmetric_difference(grp2))
    mod_df2.drop(odd_cols, axis=1, inplace=True)
    return mod_df2
    

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

def load_optimization_var_list(opt_option_str):
    """
    If 'old_model' is specified, only parameters available in the old model can be loaded
    for optimization, but if the 'new_model' is specified, then the additional precipitation
    and mass action rate constant can be optimized. 
    """
    if opt_option_str == 'old_model':
        to_optimize = [('nitrogen_ratio'), ('carbon_source_from_5e4'), ('diffusion_constant'), 
                       ('primary_ox_rate_const'), ('C'), ('N-'), ('CH4'), ('S-'), 
                       ('nitrogen_source'), ('oxygen_source'), ('fe_precipitation'), 
                       ('carbon_precip'), ('ma_op_fe_n_rate_const') ]

    elif opt_option_str == 'new_model':
        to_optimize = [('ma_op_s_n_rate_const'), ('nitrogen_ratio'), ('carbon_source_from_5e4'), 
                       ('diffusion_constant'), ('primary_ox_rate_const'), ('C'), 
                       ('ma_op_fe_n_rate_const'), ('N-'), ('CH4'), ('sminus_precipitation'), 
                       ('ma_op_ch4_n_rate_const'), ('S-'), ('nitrogen_source'), ('oxygen_source'), 
                       ('fe_precipitation'), ('carbon_precip')]
    else:
        sys.exit("illegal optimization option specified on command line")

    return to_optimize


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
    rate_idxs[5] = "Methane Oxidation (oxygen)"
    rate_idxs[6] = "Methane Oxidation (nitrate)"
    rate_idxs[7] = "Methane Oxidation (sulfate)"
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
        sys.exit("do")
    except OSError as e:
        print "error: %s" % e
        sys.exit("rae")

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
    sklearn_ss = StandardScaler()    
    
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
    sub_data_df = data_df[bool1]
    sub_obs_df = obs_df[bool2]
    sub_data_df = combine_similar_procs(sub_data_df, sub_obs_df)
    
    # drop all columns that aren't related to the specific objective
    if score_type == 'gene_objective':
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
        
        r2_array[idx] = r2_score(model_vals, obs_vals)
        print "{}: {}".format(col_, r2_array[idx])
        
    r2_array[r2_array < -1.] = -1.
    print r2_array.sum()
    return r2_array.sum()

import subprocess as sp
import platform

def run_model(arg_tuple):
    
    subdf, out_f, obs_data, score_type = arg_tuple
    print os.path.basename(out_f)
    
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
        out_h.write(out_f)
    # run the model 
    p = sp.Popen(run_cmd, cwd=model_loc, shell=True, 
                 stderr=sp.PIPE, stdout=sp.PIPE)
    stdout, stderr = p.communicate()

    # pull results & return them to memory
    
    results_df = importratesandconcs_mod(out_f, 'full df')
    
    r2 = score_results(obs_data, results_df, score_type)
    shutil.rmtree(model_loc)
    return r2

def param_lims():
    limits = {'carbon_source_from_5e4': (1e4, 1e5),
              'sminus_precipitation': (0, 0.4),
              'ma_op_s_n_rate_const': (0, 0.5),
              'ma_op_ch4_n_rate_const': (0, 100),
              'S-': (0,1.36),
              'nitrogen_source': (0, 193),
              'nitrogen_ratio': (0, 0.1),
              'oxygen_source': (6.6e2, 6.6e4),
              'methane_source': (0, 8e5),
              'fe_precipitation': (0, 0.4),
              'carbon_precip': (0, 0.4),
              'diffusion_constant': (5, 158),
              'ma_op_fe_n_rate_const': (0.6, 3), 
              'primary_ox_rate_const': (3.0e-5,3.0),
              'C':(0, 2.8), 'N-': (0, 8.24), 'CH4': (0, 10.)}
    return limits

def sample_params(params, init_bool, n_samplings, zero_pct, to_optimize, stds):
    
    parameter_mat = np.zeros ( (n_samplings, len(params.keys()))) 
    parameter_df = pd.DataFrame(index=np.arange(1,n_samplings+1),
                                columns=params.keys(), data=parameter_mat)
    for idx, key in enumerate(params.keys()):
        # start at the defaults
        this_mean = params[key]

        # only randomize for paramters that need optimization
        if key in to_optimize:
            if init_bool:
                low_, high_ = limits_[key]
                this_sample = abs(np.random.uniform(low_, high_, (n_samplings,1)))
            else:
                this_std = stds[key]
                this_sample = abs(np.random.normal(this_mean, abs(this_std)+np.finfo(float).eps, 
                                                   (n_samplings,1)))
        else:
            this_sample = np.ones((n_samplings,1))*this_mean
            
        if key != 't_max':
            this_sample = np.vstack((this_sample, np.zeros((n_zeros,1))))
        else:
            this_sample = np.vstack((this_sample, np.ones((n_zeros,1))*this_mean ))
        np.random.shuffle(this_sample)
        
        # create dataframe
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
    
def load_param_dict(type_str):
    """
    This function loads and ordered dict containing the parameters supplied to the lake
    model, to be written out into a csv in the specified order. The options include 'midpoint'
    which starts the optmization at the midpoint of the proscribed limits for all parameters, 
    'default', which loads the default parameters supplied by Sarah at the start of all of this,
    or 'midpoint-default' which starts the optimization at the midpoints of all the parameters, but
    turns off all the new precipitation & mass action rate constants that were added based on the 
    analysis of the metagenomic model
    """
    
    if type_str == 'midpoint':
        params['oxygen_bubble_rate'] = 0.0
        params['nitrogen_source'] = 96.5
        params['nitrogen_ratio'] = 0.05
        params['carbon_source_from_5e4'] = 5.5e4
        params['oxygen_source'] = 6.6e3
        params['methane_source'] = 2830.0
        params['t_max'] = 0.4
        params['sminus_precipitation'] = 0.2
        params['fe_precipitation'] = 0.2
        params['carbon_precip'] = 0.2
        params['diffusion_constant'] = 50
        params['ma_op_o_fe_rate_const'] = 10
        params['ma_op_o_n_rate_const'] = 5.0
        params['ma_op_o_s_rate_const'] = 0.16
        params['ma_op_fe_n_rate_const'] = 1.0
        params['ma_op_s_n_rate_const'] = 0.25
        params['ma_op_ch4_o_rate_const'] = 1e4
        params['ma_op_ch4_n_rate_const'] = 50.0
        params['ma_op_ch4_s_rate_const'] = 0.01
        params['primary_ox_rate_const'] = 1.5
        params['c_lim_o'] = 20.0
        params['c_lim_n'] = 5.0
        params['c_lim_fe'] = 0.1
        params['c_lim_s'] = 30
        params['C'] = 1.
        params['N-'] = 1.
        params['S-'] = 1.
        params['CH4'] = 1.
    elif type_str == 'default':
        params['oxygen_bubble_rate'] = 0.0
        params['nitrogen_source'] = 0.0
        params['nitrogen_ratio'] = 0.1
        params['carbon_source_from_5e4'] = 9.4e4
        params['oxygen_source'] = 6.6e3
        params['methane_source'] = 2830.0
        params['t_max'] = 0.4
        params['sminus_precipitation'] = 0.
        params['fe_precipitation'] = 0.3
        params['carbon_precip'] = 0.3
        params['diffusion_constant'] = 50
        params['ma_op_o_fe_rate_const'] = 10
        params['ma_op_o_n_rate_const'] = 5.0
        params['ma_op_o_s_rate_const'] = 0.16
        params['ma_op_fe_n_rate_const'] = 1.0
        params['ma_op_s_n_rate_const'] = 0.
        params['ma_op_ch4_o_rate_const'] = 1e4
        params['ma_op_ch4_n_rate_const'] = 0.
        params['ma_op_ch4_s_rate_const'] = 0.01
        params['primary_ox_rate_const'] = 1.
        params['c_lim_o'] = 20.0
        params['c_lim_n'] = 5.0
        params['c_lim_fe'] = 0.1
        params['c_lim_s'] = 30
        params['C'] = 1.
        params['N-'] = 1.
        params['S-'] = 1.
        params['CH4'] = 1.
    elif 'default-midpoint':
        params['oxygen_bubble_rate'] = 0.0
        params['nitrogen_source'] = 0.0
        params['nitrogen_ratio'] = 0.05
        params['carbon_source_from_5e4'] = 5.5e4
        params['oxygen_source'] = 6.6e3
        params['methane_source'] = 2830.0
        params['t_max'] = 0.4
        params['sminus_precipitation'] = 0.0
        params['fe_precipitation'] = 0.2
        params['carbon_precip'] = 0.2
        params['diffusion_constant'] = 50
        params['ma_op_o_fe_rate_const'] = 10
        params['ma_op_o_n_rate_const'] = 5.0
        params['ma_op_o_s_rate_const'] = 0.16
        params['ma_op_fe_n_rate_const'] = 1.0
        params['ma_op_s_n_rate_const'] = 0.0
        params['ma_op_ch4_o_rate_const'] = 1e4
        params['ma_op_ch4_n_rate_const'] = 0.0
        params['ma_op_ch4_s_rate_const'] = 0.01
        params['primary_ox_rate_const'] = 1.5
        params['c_lim_o'] = 20.0
        params['c_lim_n'] = 5.0
        params['c_lim_fe'] = 0.1
        params['c_lim_s'] = 30
        params['C'] = 1.
        params['N-'] = 1.
        params['S-'] = 1.
        params['CH4'] = 1.
    else:
        sys.exit("'midpoint' or 'default' are valid types for param dicts")
    return params

