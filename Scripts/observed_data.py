#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 19:27:43 2017

@author: login
"""
from sklearn.preprocessing import MinMaxScaler
import os
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np

def load_chem_data():
    """This function loads observed oxidized chemical
    specis"""
    chem_base = os.path.join(os.getcwd(), "ChemData")
    
    files = [i for i in os.listdir(chem_base) if "NEW" in i]
    file_ps = [os.path.join(chem_base, i) for i in files]
    mats = {}
    
    for f, p in zip(files, file_ps):
        if f.endswith("csv"):
            temp_df = pd.read_csv(p, index_col=0).T
            new_idx = pd.to_datetime(temp_df.index)
            temp_df.index = new_idx
            mats[f] = temp_df.stack()
        elif f.endswith("xlsx"):
            mats[f] = pd.read_excel(p, index_col=0).T.stack()
            
    del mats['oxygen_mg_NEW.xlsx']
    del mats['sulfate_mg_NEW.xlsx']
    del mats['nitrate_mg_NEW.xlsx']
    
    columns_ = ['Fe-', 'N+', 'O', 'S+']
    match_keys = ['iron_um_NEW.xlsx', 'nitrate_um_NEW.xlsx', 
                  'oxygen_um_NEW.xlsx', 'sulfate_um_NEW.csv']
                  
    n_rows = len(mats['oxygen_um_NEW.xlsx'].index)
    n_cols = len(columns_)
    
    obs_df = pd.DataFrame(index=mats['oxygen_um_NEW.xlsx'].index,
                          columns=columns_, data=np.zeros((n_rows, n_cols)) )
    
    for c, mk in zip(columns_, match_keys):
        obs_df.ix[:, c] = mats[mk]
    
    obs_df.index.names = ['Date', 'Depth']

    return obs_df
    
def load_gene_data():
    chem_base = os.path.join(os.getcwd(), "ChemData")
    gene_dists_f = os.path.join(chem_base, 'final_process_genes_2.xlsx')
    raw_dists = pd.read_excel(gene_dists_f, index_col=1)
    real_dists = raw_dists[raw_dists.Process != 'None']
    dists_only = real_dists.drop(['KEGG ID', 'Process'], axis=1)
    
    real_dists_mat = dists_only.T.values
    X = MinMaxScaler().fit_transform(real_dists_mat).T
    
    real_dists_std = pd.DataFrame(data=X, columns=dists_only.columns,
                                  index=dists_only.index)
    real_dists_std['Process'] = real_dists.Process
    process_df_sq = real_dists_std.groupby(by='Process').mean()
    process_df = process_df_sq.unstack()
    process_df = process_df.reset_index()
    process_df.columns = ['Depth', 'Process', 'Value']
    process_df['Date'] = pd.to_datetime(['08-12-2013']*process_df.shape[0])
    process_df = process_df.set_index(['Date', 'Depth'])
    process_df = process_df.pivot_table(values='Value', 
                           index=process_df.index, 
                           columns=['Process'])
    process_df.index = pd.MultiIndex.from_tuples(process_df.index)
    process_df.index.names = ['Date', 'Depth']
    return process_df
    
