from scipy.stats import ttest_rel
import pandas as pd
from scipy.spatial.distance import pdist, squareform
import numpy as np
import os, sys

def distance_extractor(df1, df2, dist_type, check_integrity=True):
    assert df1.shape[0] > df1.shape[1]
    assert df1.shape == df2.shape
    if check_integrity:
        assert (df2.index == df1.index).sum() == len(df1.index)
    assert (df2.columns == df1.columns).sum() == df2.shape[1]
    mat1, mat2 = df1.values, df2.values
    super_mat = np.vstack((mat1, mat2))
    full_dist =  squareform(pdist(super_mat, metric=dist_type))
    paired_dists = [full_dist[i, i+mat1.shape[0]] for i in range(mat1.shape[0])]
    return pd.Series(index=df1.index, data=paired_dists)

def measure_separation(dist_type, df_A, df_B, return_data=False):
    train_df, test_df = df_A.copy().reset_index(level=0), df_B.copy().reset_index(level=0)
    dist_choice = dist_type
    test_seqs = set(test_df.OTU)
    comprable_set = list(test_seqs & set(train_df.OTU))
    abund_cols = test_df.columns[1:12]
    test_c_df = test_df[test_df.OTU.isin(comprable_set)]
    train_c_df = train_df[train_df.OTU.isin(comprable_set)]
    test_set_df, train_set_df = test_c_df.set_index("OTU"), train_c_df.set_index("OTU")
    test_matched, train_matched = test_set_df.ix[comprable_set, abund_cols], train_set_df.ix[comprable_set, abund_cols]
    abund_dists = distance_extractor(train_matched, test_matched, dist_choice, True)

    rand_tr_idxs = np.random.choice(train_df.index, size=(len(comprable_set),), replace=False)
    rand_te_idxs = np.random.choice(test_df.index, size=(len(comprable_set),), replace=False)
    rand_tr_df = train_df.ix[rand_tr_idxs, :].set_index("OTU")
    rand_te_df = test_df.ix[rand_te_idxs, :].set_index("OTU")
    train_rand = rand_tr_df.ix[:, abund_cols]
    test_rand = rand_te_df.ix[:, abund_cols]

    rand_dists = distance_extractor(train_rand, test_rand, dist_choice, False)

    if return_data:
    	return (abund_dists, rand_dists)
    else:
        all_dists = np.hstack((rand_dists.values, abund_dists.values))
        mean_dist = -1.0*(rand_dists.values.mean() - abund_dists.values.mean())/ all_dists.mean()
        median_dist = -1.0*(np.median(rand_dists.values) - np.median(abund_dists.values)) / np.median(all_dists)
        _ , pval = ttest_rel(a=rand_dists.values, b=abund_dists.values)
        sumRatio = abund_dists.sum() / rand_dists.sum()
        return (mean_dist, median_dist, pval, sumRatio)
