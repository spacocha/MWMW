"""
python kegg_cluster.py -f ../Data/KEGG_Annotations/KO_Database/TO_Files
python kegg_cluster.py -t ../Data/KEGG_Annotations/KO_Database/KO_by_TO_Table.tsv 
python kegg_cluster.py -b ../Data/KEGG_Annotations/Aggregated_Annotations.tsv
"""

import os, sys
import pandas as pd
import numpy as np
from itertools import chain
from scipy.spatial.distance import pdist, squareform, cdist

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
    TO_df, Bin_df = full_rep_TO.ix[:, common_kos], full_rep_bin.ix[:, common_kos]
    assert TO_df.isnull().sum().sum() == 0
    assert Bin_df.isnull().sum().sum() == 0
    bin_26 = Bin_df.ix['bin_26', :].values.reshape(1, len(common_kos))
    dist_26_seuc = pd.Series(data=cdist(bin_26, TO_df.values, 'seuclidean').flatten(), index=TO_df.index)
    dist_26_euc = pd.Series(data=cdist(bin_26, TO_df.values, 'euclidean').flatten(), index=TO_df.index)
    dist_26_cos = pd.Series(data=cdist(bin_26, TO_df.values, 'cosine').flatten(), index=TO_df.index)
    dist_26_cor = pd.Series(data=cdist(bin_26, TO_df.values, 'correlation').flatten(), index=TO_df.index)
    # T00784 ecy, ECOLX, 409438; Escherichia coli O152:H28 SE11 (commensal strain)
    # T00944 ebr, 413997; Escherichia coli B REL606
    print dist_26_seuc.argmin()
    bin_52 = Bin_df.ix['bin_52', :].values.reshape(1, len(common_kos))
    dist_52s = pd.Series(data=cdist(bin_52, TO_df.values, 'seuclidean').flatten(), index=TO_df.index)
    # T00647 Microcystis aeruginosa NIES-843
    def match_a_bin(binname, bindf, todf, col_set):
        a_bin_annots = bindf.ix[binname, :].values.reshape(1, len(col_set))
        bin_dists = pd.Series(data=cdist(a_bin_annots, todf.values, 'seuclidean').flatten(), index=todf.index)
        print "{}\t{}\t{}\t".format(binname, bin_dists.argmin(), bin_dists.min())
        return (binname, bin_dists.argmin(), bin_dists.min())

    matches = [match_a_bin(i, Bin_df, TO_df, common_kos) for i in Bin_df.index]




    #dist_mat = squareform(pdist(full_df.values, 'correlation'))
    # save this as a pickle
    # 
    # plot distributions
    




