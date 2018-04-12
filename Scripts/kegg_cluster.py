"""
# this option reads all the individual flat files in the specified dir into a single csv
python kegg_cluster.py -f ../Data/KEGG_Annotations/KO_Database/TO_Files

# this option was for debugging
python kegg_cluster.py -t ../Data/KEGG_Annotations/KO_Database/KO_by_TO_Table.tsv 

# this option reads in the output of [one of the other scripts in this repo] and reformats
# to match the output of the -f option above
python kegg_cluster.py -b ../Data/KEGG_Annotations/Aggregated_Annotations.tsv

# this option collects the output of the '-b' and '-f' and ....
python kegg_cluster.py -bt

"""

import os, sys, urllib2, time
import pandas as pd
import numpy as np
from itertools import chain
from scipy.spatial.distance import pdist, squareform, cdist
from scipy.stats import spearmanr
from sklearn.preprocessing import normalize

def unit_norm_row_col(df):
    mat = df.values.T
    print df.index[0], df.shape, df.columns[0]
    nmat, norms_ = normalize(mat, norm='l1', axis=0, return_norm=True)
    print norms_.shape, "normalizing along columns"
    nmat_2, norms_ = normalize(nmat, norm='l1', axis=1, return_norm=True)
    print norms_.shape, "normalizing along rows"
    norm_df = pd.DataFrame(index=df.index, columns=df.columns, data=nmat_2.T)
    return norm_df

args = sys.argv

if '-f' in args:
    to_tuples = [(os.path.join(args[2], i), i[:-4]) for i in sorted(os.listdir(args[1]))]
    def read_kegg_table(TO_tuple):
        TO_path, TO_number = TO_tuple
        with open(TO_path, 'r') as fh:
            content = [i.split(":")[-1] for i in fh.read().split("\n") if i != ""]
        kos, cnts = np.unique(content, return_counts=1)
        return pd.Series(data=cnts, index=kos, dtype=np.float64, name=TO_number)

    TO_srs_lst = [read_kegg_table(to_tple) for to_tple in to_tuples]
    full_rep = pd.concat(TO_srs_lst, axis=1, keys=[s.name for s in TO_srs_lst]).T
    full_rep.to_csv("../Data/KEGG_Annotations/KO_Database/KO_by_TO_Table.tsv", sep="\t", index=True, header=True)

if '-t' in args:
    full_rep = pd.read_csv(args[2], sep="\t", index_col=0).fillna(0)

# read in bin data
if '-b' in args:
    cols_ = ["Protein_Name", "Bin_Name", "KO_Annot_1", "KO_Annot_2", "KO_Annot_3", "K_number"]
    kegg_df = pd.read_csv(args[2], sep="\t", usecols=cols_, low_memory=False)
    bin_ko_lists = []
    for b in kegg_df.Bin_Name.unique():
        print b
        bko_df = kegg_df[kegg_df.Bin_Name == b]
        b_list = []
        for _ix_ in bko_df.index:
            a_rec = bko_df.ix[_ix_, bko_df.columns[-4:]].dropna().unique().tolist()
            b_list += list(set(chain.from_iterable([i.split(",") for i in a_rec])))
        kos, cnts = np.unique(b_list, return_counts=1)
        bin_ko_lists.append(pd.Series(data=cnts, index=kos, dtype=np.uint64, name=b))
    full_rep = pd.concat(bin_ko_lists, axis=1, keys=[s.name for s in bin_ko_lists]).T
    full_rep.to_csv("../Data/KEGG_Annotations/KO_by_Bin_Table.tsv", sep="\t", index=True, header=True)
    
# find nearest neighbors in correlation distance

if '-bt' in args:
    bin_ko_df_fn = "../Data/KEGG_Annotations/KO_by_Bin_Table.tsv"
    TO_ko_df_fn = "../Data/KEGG_Annotations/KO_Database/KO_by_TO_Table.tsv"
    full_rep_bin = pd.read_csv(bin_ko_df_fn, sep="\t", index_col=0).fillna(0)
    full_rep_TO = pd.read_csv(TO_ko_df_fn, sep="\t", index_col=0).fillna(0)
    bin_ko_labels, to_ko_labels = set(full_rep_bin.columns), set(full_rep_TO.columns)
    common_kos = bin_ko_labels & to_ko_labels
    bin_only_kos = bin_ko_labels - to_ko_labels
    sym_dis_kos = bin_ko_labels ^ to_ko_labels
    print "KOs found in all genomes {}".format(full_rep_TO.shape[1])
    print "KOs found in the bins: {}".format(full_rep_bin.shape[1])
    print "Found in both data sets: {}".format(len(common_kos))
    print "Found only in the bins: {}".format(len(bin_only_kos))
    TO_df_pre, Bin_df_pre = full_rep_TO.ix[:, common_kos], full_rep_bin.ix[:, common_kos]
    TO_df, Bin_df = unit_norm_row_col(TO_df_pre), unit_norm_row_col(Bin_df_pre)


    assert TO_df.isnull().sum().sum() == 0
    assert Bin_df.isnull().sum().sum() == 0
    
    def match_a_bin(binname, bindf, todf, dist_str, n_matches):
        if dist_str == 'spearman_corr':
            a_bin_annots = bindf.ix[binname, :].values
            sp_dists = [spearmanr( np.vstack((a_bin_annots, todf.ix[i,:].values)).T)[0] for i in todf.index]
            sp_d_arr = np.array(sp_dists)
            match_df = pd.Series(index=todf.index, data=sp_d_arr).sort_values(ascending=False)[:n_matches]
        else:
            a_bin_annots = bindf.ix[binname, :].values.reshape(1, bindf.shape[1])
            data_dist_mat = cdist(a_bin_annots, todf.values, dist_str).flatten()
            match_df = pd.Series(index=todf.index, data=data_dist_mat).sort_values()[:n_matches]
            
        match_tuples = zip(match_df.index, match_df.values)
        match_flat = list(chain.from_iterable(match_tuples))
        return match_flat

    n_matches_ = 5
    spr_matches, sce_matches = [], []
    for i_bin in Bin_df.index:
        print i_bin, "is being matched"
        sce_matches.append(match_a_bin(i_bin, Bin_df, TO_df, 'seuclidean', n_matches_))
        spr_matches.append(match_a_bin(i_bin, Bin_df, TO_df, 'spearman_corr', n_matches_))

    output_columns = [('Hit {}'.format(i), 'Score {}'.format(i)) for i in xrange(1,n_matches_+1)]
    out_cols = list(chain.from_iterable(output_columns))
    non_para_idx = [i+"_a" for i in Bin_df.index]
    para_idx = [i+"_b" for i in Bin_df.index]
    para_df = pd.DataFrame(index=para_idx, data=np.array(sce_matches), columns=out_cols)
    nonpara_df = pd.DataFrame(index=non_para_idx, data=np.array(spr_matches), columns=out_cols)

    hit_cols = [i for i in out_cols if "Hit" in i]
    score_cols = out_cols[0]
    print (nonpara_df.ix[:, score_cols] == para_df.ix[:, score_cols]).sum(), "aggrements between distance metrics"
    non_para_matches = set(nonpara_df.ix[:, hit_cols].values.flatten())
    non_para_matches.update(set(para_df.ix[:, hit_cols].values.flatten()))

    print len(non_para_matches), "TO names will be fetched"

    def fetch_name(this_to):
        pre_url = "http://rest.kegg.jp/find/genome/"
        data = urllib2.urlopen(pre_url+this_to).read()
        this_tos_name = data.strip().split(";")[-1].strip()
        print ": ".join([this_to, this_tos_name])
        time.sleep(2.)
        return ": ".join([this_to, this_tos_name])
    
    to_name_dict = {u:fetch_name(u) for u in list(non_para_matches)}
    fix_name = lambda x: to_name_dict[x]

    para_df['metric'] = ['standardized euclidean dist'] * para_df.shape[0]
    nonpara_df['metric'] = ['spearman correlation'] * nonpara_df.shape[0]

    for hc in hit_cols:
        nonpara_df.ix[:, hc] = nonpara_df.ix[:, hc].apply(fix_name)
        para_df.ix[:, hc] = para_df.ix[:, hc].apply(fix_name)

    assert (para_df.columns == nonpara_df.columns).sum() == len(nonpara_df.columns) == len(para_df.columns)

    output_df = para_df.append(nonpara_df).sort_index()
    output_df.to_csv("../Data/min_ko_distance_genomes_normed2x.tsv", sep="\t", index_label="Bin")
    




