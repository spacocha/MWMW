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

def taxa_stepper(bin_label, distdf, bintaxa, testtaxa, tol_):
    # slice row of distances and taxa classes assoc. w/ one bin 
    this_bin_matches = distdf.ix[bin_label, :]
    this_bin_taxa = bintaxa.ix[bin_label, :]
    # count the depth of phylogenetic classification (down to genus)
    score_denom = float(6 - list(this_bin_taxa.apply(len)).count(0))
    # iterate through each taxa level
    for t_label, t_name in zip(this_bin_taxa.index, this_bin_taxa.values):
        if t_name != "":
            # pull out each amplicon name with matching taxa label
            step_1_t = testtaxa[testtaxa.ix[:, t_label] == t_name].index
            # subtract 1/`score_denom`all of their the distance values
            this_bin_matches[step_1_t]-=(1./score_denom)
        else:
            break
    # find the minimum score (best hit)
    the_best_score = this_bin_matches.min()
    # extract all names assoc with that hit (ties possible)
    best_hitters = list(this_bin_matches[this_bin_matches < (the_best_score+tol_)].index)
    # remove the best hits and make a copy of the distances
    nex_best_vec = this_bin_matches.drop(best_hitters)
    # pull out the next best hit from the shortened vector
    next_best_score = nex_best_vec.min()
    # calculate the scaled difference to the next best score 
    diff = abs(the_best_score)-abs(next_best_score)
    weight = (diff)+abs(the_best_score)
    # return the bin label, the winners, the confidence space, the modified ditances, and the best score
    return (bin_label, best_hitters, weight, this_bin_matches, the_best_score, diff)

def taxa_iterator(label_list, dist_df_, bin_taxa, test_taxa_, tol_):
    matched_tags, dist_vectors = [], []
    while len(label_list) > 0:
        local_dists, local_tt = dist_df_.copy(), test_taxa_.copy()
        practice_hits = [taxa_stepper(ll, local_dists, bin_taxa, local_tt, tol_) for ll in label_list]
        bin_matched, tags_matched, weight, this_bin_matchd, best_score, diff = sorted(practice_hits, key=lambda k: k[2], reverse=True)[0]
        dist_vectors.append(this_bin_matchd)
        label_list.remove(bin_matched)
        local_dists.drop(tags_matched, axis=1, inplace=True)
        local_tt.drop(tags_matched, axis=0, inplace=True)
        am_ = dist_df_.ix[bin_matched, :].min()
        bm_r = dist_df_.ix[bin_matched, :]
        raw_hits = list(bm_r[bm_r < (am_+tol_)].index)
        mt_dict = {"Bin":bin_matched, "Tag":tags_matched, "Weight":weight, 
                   "Min": best_score, "RawMin": am_, "rTags": raw_hits, "Diff": diff}
        matched_tags.append(mt_dict)
        print "{} matched to {} tags, {} iterations remaining".format(bin_matched, tags_matched, len(label_list))
    return (matched_tags, dist_vectors)

def join_bins_and_tags(bin_df, tag_df, tol_, splitN, replacement):
    bin_abunds = bin_df.ix[:, bin_df.columns[:splitN]]
    tag_abunds = tag_df.ix[:, tag_df.columns[:splitN]]
    bin_taxa = bin_df.ix[:, bin_df.columns[splitN:]]
    tag_taxa = tag_df.ix[:, tag_df.columns[splitN:]]
    mat1, mat2 = bin_abunds.values, tag_abunds.values
    super_mat = np.vstack((mat1, mat2))
    full_dist =  squareform(pdist(super_mat, metric='euclidean'))
    bin_to_tags = full_dist[:bin_abunds.shape[0], bin_abunds.shape[0]:]
    dist_df = pd.DataFrame(index=bin_abunds.index, columns=tag_abunds.index, data=bin_to_tags)
    raw_dist_df = dist_df.copy()
    label_list = bin_df.index.tolist()
    if replacement:
        match_data = [taxa_stepper(ll, dist_df, bin_taxa, tag_taxa, tol_) for ll in label_list]
        bin_labels, hit_lists, weights, dist_vectors, best_scores, diffs = zip(*match_data)
        matched_tags = []
        for bm, tm, w, bs, ds in zip(bin_labels, hit_lists, weights, best_scores, diffs):
            am_ = raw_dist_df.ix[bm, :].min()
            bm_r = raw_dist_df.ix[bm, :]
            raw_hits = list(bm_r[bm_r < (am_+tol_)].index)
            matched_tags.append({"Bin":bm, "Tag":tm, "Weight":w, "Min": bs, "RawMin": am_, "rTags": raw_hits, "Diff": ds})
    else:
        matched_tags, dist_vectors = taxa_iterator(label_list, dist_df, bin_taxa, tag_taxa, tol_)
    mod_dist_df = pd.concat(dist_vectors, axis=1, keys=[s.name for s in dist_vectors]).T
    to_write = pd.DataFrame(matched_tags).set_index("Bin", verify_integrity=1).sort_values(["Weight"], ascending=False)
    to_write["refined_matches"] = to_write.Tag.apply(len)
    to_write["raw_matches"] = to_write.rTags.apply(len)
    return (mod_dist_df, to_write, raw_dist_df)

def make_and_write_reports(tag_df, bin_df, to_write_df, write_dir, fn_list, dist_mats):
    path_list = [os.path.join(write_dir, f) for f in fn_list]
    to_write_df.to_csv(path_list[0], sep="\t", index_label="Bin", index=1, header=1)
    dist_mats[0].to_csv(path_list[2], sep="\t", index=1, header=1)
    dist_mats[1].to_csv(path_list[3], sep="\t", index=1, header=1)
    full_rep_vecs = []
    for bin_no in to_write_df.index:
        this_bin_vec = bin_df.ix[bin_no, :].copy()
        this_bin_vec.loc['Group'] = bin_no
        full_rep_vecs.append(this_bin_vec)
        for match_no in to_write_df.ix[bin_no, 'Tag']:
            this_tag_vec = tag_df.ix[match_no, :].copy()
            this_tag_vec.loc['Group'] = bin_no
            full_rep_vecs.append(this_tag_vec)
    full_rep = pd.concat(full_rep_vecs, axis=1, keys=[s.name for s in full_rep_vecs]).T
    fr_cols = [full_rep.columns[-1]] + full_rep.columns[:-1].tolist()
    full_rep_ro = full_rep[fr_cols]
    full_rep_ro.to_csv(path_list[1], sep="\t", index=1, header=1)
    return
