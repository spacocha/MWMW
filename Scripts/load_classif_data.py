"""
This is a module containing functions for loading in the 
abundance and taxanomic matrixes required to assign 16S to bins.

You can load in each normalized matrix with a specific function

# consider: $MARCCDIR/work/sprehei1/Keith_Files/kaw/xMetaPipeline/bin/taxFormatter.py
# consider: $MARCCDIR/work/sprehei1/Keith_Files/kaw/xMetaPipeline/data/centrifuge_*2.tsv
# consider: $MARCCDIR/work/sprehei1/Keith_Files/kaw/xMetaPipeline/data/hiseq_ko_taxonomy2.txt
# creates a list of root samples to create and pairs to merge to create them

"""

import time
import numpy as np
import pandas as pd
import sys
from sklearn.preprocessing import normalize

### Global variable definitions
otu_taxa_file = "../Data/16S_Info/unique.dbOTU.nonchimera.ns.fasta.rdp2.txt"

col_corresp = {"A4_3mMystic": "SB081213TAWMD03VV4TM",
               "A5_5mMystic": "SB081213TAWMD05VV4TM",
               "A6_7mMystic": "SB081213TAWMD07VV4TM",
               "A7_9mMystic": "SB081213TAWMD09VV4TM",
               "A8_11mMystic": "SB081213TAWMD11VV4TM",
               "B1_13mMystic": "SB081213TAWMD13VV4TM",
               "B2_15mMystic": "SB081213TAWMD15VV4TM",
               "B3_17mMystic": "SB081213TAWMD17VV4TM",
               "B5_20mMystic": "SB081213TAWMD20VV4TM",
               "B6_21mMystic": "SB081213TAWMD21VV4TM",
               "B7_22mMystic": "SB081213TAWMD22VV4TM"}

### Function Definitions

# remove undetected sequences
drop_zero_cols = lambda df: df.T[df.sum() != 0].T

def parseBiosample(df):
    """
    This function accepts a df with samples in rows and OTUs in columns.
    It extract metadata from "biosample" string and adds them as columns.
    """
    biosamples = list(df.index)
    dates_, primers_, kits_, replicates_, depths_ = [], [], [], [], []
    for bs in biosamples:
        if bs[:2] == "SB":
            clipped = bs[2:]
        else:
            sys.exit("Non-SB start to biosample")
        
        if 'TAWMD' in clipped:
            date, rest = clipped.split("TAWMD")
        elif "TAWWD" in clipped:
            date, rest = clipped.split("TAWWD")
        else:
            sys.exit("Bridge Sequence not Detected")
            
        if "VV4T" in rest:
            primer = "VV4"
            depth, rest2 = rest.split("VV4T")
        elif "VV4V5T" in rest:
            primer = "V4V5"
            depth, rest2 = rest.split("VV4V5T")
        else:
            sys.exit("Primer Sequence not Detected")
        
        if rest2[0] == "M":
            kit = "MolBio"
        elif rest2[0] == "Q":
            kit = "Qiagen"
        elif rest2[:4] == "filt" and rest2[4] == "M":
            kit = "MolBio"
        elif rest2[:2] == "NA":
            kit = "NA"
        else:
            print clipped
            print rest2[0]
            sys.exit("kit type not detected")
        
        if rest2[-2] == "R":
            replicate = rest2[-1]
        else:
            sys.exit("replicate signifier not detected")
        
        if depth == '015':
            depth = '01.5'
        
        dates_.append(date); primers_.append(primer); kits_.append(kit); replicates_.append(replicate); depths_.append(depth)
    
    df['date'], df['primers'], df['kit'], df['replicates'], df['depth'] = dates_, primers_, kits_, replicates_, depths_
    return df

def append_indv_col_taxa_to_df(taxa_df, fulldf):
    """
    this appends the taxa columns to the abundances columns
    """
    taxa_subdf = taxa_df.ix[list(fulldf.index), :]
    assert taxa_subdf.isnull().sum().sum() == 0
    fullerdf = pd.concat([fulldf.copy(), taxa_subdf], axis=1, verify_integrity=True)
    return fullerdf

def combine_replicate_pairs(replicate_trios, old_df):
    """
    This merges (if requested) l1 normalized replicates, and reports minimum correlation
    """
    roots, p1s, p2s = zip(*replicate_trios)
    new_df = pd.DataFrame(index=roots, columns=old_df.columns)
    corr_scores = []
    for root, p1, p2 in replicate_trios:
        if p1 in old_df.index and p2 in old_df.index:
            sub_old_df = old_df.ix[[p1, p2], :] / 2.0
            corr_scores.append(sub_old_df.T.corr().ix[p1, p2])
        elif p1 in old_df.index:
            sub_old_df = old_df.ix[[p1], :]
        elif p2 in old_df.index:
            sub_old_df = old_df.ix[[p2], :]
        else:
            raise ValueError("illegal column detected")
        new_df.ix[root, :] = sub_old_df.div(sub_old_df.sum(1), 0).sum()
        
    print "Minimum replicate correlation: {}".format(np.min(corr_scores).round(5))
    rescale_values = (1.0 / new_df[new_df > 0].min(1))
    rescaled_df = new_df.mul(rescale_values, 0)
    return rescaled_df

def drop_low_conf_and_IS(label, confi_label, cutoff, df):
    """
    Takes a DF containing fullrank RDP classifications and removes
    low confidence assignments and any gaps in the hierarchy noted
    by the phrase Incertae Sedis
    """
    df_copy = df.copy()
    df_bool = df.ix[:,confi_label] < float(cutoff)
    df_copy.ix[df_bool, label] = ""
    df_bool2 = df_copy.ix[:, label].str.contains("Incertae")
    df_copy.ix[df_bool2, label] = ""
    return df_copy

def clean_rdp_txt(otu_taxa_file, con_cutoff):
    """
    This takes the raw fixed rank RDP file and returns a cleaned
    up 6 column hierarchy with only high confidence classifications
    """
    with open(otu_taxa_file, "r") as fh:
        raw_taxa_otu = [ i for i in fh.read().replace("%",".").replace('"', "").split("\n") if i !=""][6:]
    taxa_cols = raw_taxa_otu[0].split(";")
    raw_taxa_otu.remove(raw_taxa_otu[0])
    taxa_splits = [i.split(";") for i in raw_taxa_otu]
    tax_arr = np.array(taxa_splits)
    taxa_df_raw = pd.DataFrame(index=tax_arr[:,0], columns = taxa_cols[2:], data=tax_arr[:, 2:])
    taxa_df_raw.sort_index(inplace=True)
    taxa_df_dt = taxa_df_raw.apply(pd.to_numeric, errors='ignore')
    taxa_df_gd = drop_low_conf_and_IS('Genus', 'G_con', con_cutoff, taxa_df_dt)
    taxa_df_fd = drop_low_conf_and_IS("Family",  "F_con", con_cutoff, taxa_df_gd)
    taxa_df_od = drop_low_conf_and_IS("Order",  "O_con", con_cutoff, taxa_df_fd)
    taxa_df_cd = drop_low_conf_and_IS("Class",  "C_con", con_cutoff, taxa_df_od)
    taxa_df_pd = drop_low_conf_and_IS("Phylum",  "P_con", con_cutoff, taxa_df_cd)
    taxa_df_kd = drop_low_conf_and_IS("Kingdom",  "K_con", con_cutoff, taxa_df_pd)
    # drop blank species column
    taxa_df = taxa_df_kd.drop([i for i in taxa_df_kd.columns if "_con" in i ], axis=1)
    return taxa_df

taxa_df = clean_rdp_txt(otu_taxa_file, 50.0)

def l1l2clr_norm(df, n_type, psct=None):
    """
    Accepts and returns df with [n_samples, n_features] after performing 
    either center-log transform or 'l1' or 'l2' normalization, as specified.
    """
    if n_type.startswith("l"):
        mat = df.values.astype(float)
        data_ = normalize(mat, norm=n_type, axis=1, copy=True)
    elif n_type == 'clr' and psct:
        mat = df.values.astype(float) + psct
        mat = (mat / mat.sum(axis=1, keepdims=True)).squeeze()
        lmat = np.log(mat.astype(float))
        gm = lmat.mean(axis=-1, keepdims=True)
        data_ = (lmat - gm).squeeze()
    else:
        raise ValueError("l1/l2/clr must be specified and psuedocount value as well if clr")
    return pd.DataFrame(index=df.index, columns=df.columns, data=data_)

def import_amplicon_matrix(norm_type, write_bool, test_label_set=None, psct_val=None):
    otu_matrix_file = "../Data/16S_Info/unique.dbOTU.nonchimera.mat.rdp.local.tsv"
    # load in test data
    otu_table = pd.read_csv(otu_matrix_file, sep="\t", index_col=0)
    # drop taxa series and bad columns first
    otu_time_table = otu_table.drop([otu_table.columns[-1]], axis=1).sort_index()
    do_not_use = ['SB100912TAWMD14VV4TMR1', 'SB061713TAWMD22VV4TMR1', 'SB011413TAWMD22VV4TMR1', 'SB011413TAWMDSBVV4TMR1', 'SB011413TAWMDEBVV4TMR1', 'SB011413TAWMD22VV4TMR2', 'SB011413TAWMDSBVV4TMR2']
    otu_time_table.drop(do_not_use, axis=1, inplace=True)
    # add metadata and drop undated columns
    samps_by_feats = parseBiosample(otu_time_table.T)
    undated = samps_by_feats[samps_by_feats.date == "NA"].index
    samps_by_feats.drop(undated, inplace=True)
    # parse dates
    samps_by_feats['date'] = pd.to_datetime(samps_by_feats['date'])
    # Drop samples of unknown provenance 
    squishy_depths = ['bottom', 'SB', 'control1', 'control3', 'control6', 'mid1', 'mid2', 'mid3', 'river', 'upper', 'River', 'Neg', 'NA', 'EB', 'FB', 'CR', 'Dock', 'Bridge']
    manual_date = np.datetime64("2013-08-12")
    putative_rep_trios = [(i, i+"R1", i+"R2") for i in sorted(col_corresp.values())]
    off_target = samps_by_feats[samps_by_feats.depth.isin(squishy_depths)].index
    only_depths = samps_by_feats.drop(off_target)
    # filter rows by date and drop metadata columns
    aug_abund_df = only_depths[only_depths.date == manual_date]
    row_set_abunds = aug_abund_df.ix[:, aug_abund_df.columns[:-5]]
    combo_train_df = combine_replicate_pairs(putative_rep_trios, row_set_abunds)
    combo_train = drop_zero_cols(combo_train_df)

    if test_label_set:
        train_label_set = set(list(combo_train.columns))
        print "Train labels are a subset of test labels: {}".format(test_label_set.issubset(train_label_set))
        lost_labels = test_label_set - train_label_set
        print "\t{} labels not found:".format(len(lost_labels))

    if norm_type != "raw":
        normed_otus = l1l2clr_norm(combo_train, norm_type, psct_val)
        print "Performed {} scaling on matched amplicons size {}".format(norm_type, normed_otus.shape)
    else:
        normed_otus = combo_train
        print "Performed no scaling on matched amplicons size {}".format(normed_otus.shape)

    train_df = append_indv_col_taxa_to_df(taxa_df, normed_otus.T)

    if write_bool:
        train_fname="otu_abund_tax_train.tsv"
        train_df.to_csv("../Data/16S_Info/"+train_fname, sep="\t", header=True, index=True)

    return train_df

def import_matched_amplicons(norm_type, filter_samples, write_bool, return_labels, psct_val=None):
    """
    """
    shotgun_abunds = "../Data/16S_Info/dbOTUbySample.tsv"
    smg_table = pd.read_csv(shotgun_abunds, sep="\t", index_col=0)
    if filter_samples:
        replicate_trios = [(i, i+"_1.sam", i+"_2.sam") for i in sorted(col_corresp.keys())]
    else:
        replicate_trios = [(i, i+"_1.sam", i+"_2.sam") for i in sorted(smg_table.bSample.unique())]

    just_tags = smg_table.ix[:, smg_table.columns[2:]]
    combo_df = combine_replicate_pairs(replicate_trios, just_tags)
    sample_select = drop_zero_cols(combo_df).T

    if filter_samples:
        new_cols = [col_corresp[i] for i in sample_select.columns]
        sample_select.columns = new_cols

    print "{} otus detected in shotgun libraries".format(sample_select.shape[0])
    test_label_set = set(list(sample_select.index))

    sample_select_t = sample_select.T
    if norm_type != "raw":
        normed_tags = l1l2clr_norm(sample_select_t, norm_type, psct_val)
        print "Performed {} scaling on matched amplicons size {}".format(norm_type, normed_tags.shape)
    else:
        normed_tags = sample_select_t
        print "Performed no scaling on matched amplicons size {}".format(normed_tags.shape)

    test_df = append_indv_col_taxa_to_df(taxa_df, normed_tags.T)
    test_df.index.name = "OTU"

    if write_bool:
        test_data_fname = "smg_abund_tax_test.tsv"
        test_df.to_csv("../Data/16S_Info/"+test_data_fname, sep="\t", header=True, index=True)
    
    if return_labels:
        return (test_label_set, test_df)
    else:
        return test_df

def import_bin_abundances(normalize_to_reads):
    """
    Imports bin abundances produced by metaWRAP. Normalizes them as recommended by 
    program's authors
    """
    bin_abund_file = "../Data/Bin_Abundances/bin_abundance.tab"
    bin_normalizer_f = "../Data/Bin_Abundances/sample_read_count.tab"
    abund_df = pd.read_csv(bin_abund_file, sep="\t")
    sorted_cols = sorted(list(abund_df.columns)[1:])
    fix_bin_names = lambda x: x.replace("b", "B").replace(".", " ")
    abund_df.ix[:, "Genomic bins"] = abund_df.ix[:, "Genomic bins"].apply(fix_bin_names)
    abund_df.ix[:, "Bin_num"] = abund_df.ix[:, "Genomic bins"].apply(lambda x: int(x.split()[-1]))
    sorted_df = abund_df.sort_values(["Bin_num"]).set_index("Genomic bins").ix[:, sorted_cols].T
    if normalize_to_reads:
        denom_df = pd.read_csv(bin_normalizer_f, sep="\t", index_col=0).sort_index()
        assert (sorted_df.index == denom_df.index).sum() == len(denom_df.index)
        normed_df = pd.DataFrame(index=sorted_df.index, 
                                 columns=sorted_df.columns, 
                                 data=sorted_df.values/denom_df.values)
        normed_scaled_df = normed_df * 1e6
        return normed_scaled_df
    else:
        return sorted_df

def check_taxa_hierarchies(bin_taxa_lr, taxa_df):
    bin_taxa_sets = [set(bin_taxa_lr.ix[:, i].tolist()) for i in bin_taxa_lr.columns]
    otu_taxa_sets = [set(taxa_df.ix[:, i].tolist()) for i in taxa_df.columns]
    taxa_is_subset = [(s.issubset(t)) for s, t in zip(bin_taxa_sets, otu_taxa_sets)]
    if sum(taxa_is_subset) == len(taxa_is_subset):
        print "All bin taxa levels are in otu taxa hierarchy!"
    else:
        taxa_diff = {n:(s - t) for n, s, t in zip(bin_taxa_lr.columns, bin_taxa_sets, otu_taxa_sets)}
        print taxa_diff
    return

def import_bin_data(norm_type, filter_samples, write_bool, check_taxa=True, psct_val=None):
    bin_taxa_file = "../Data/bin_taxonomy_edit.tsv"
    bin_taxa_raw = pd.read_csv(bin_taxa_file, sep="\t", index_col=0)
    good_cols = [i for i in bin_taxa_raw.columns if "%" not in i]
    good_cols.remove('Species')
    bin_taxa_lr = bin_taxa_raw.ix[:, good_cols].fillna("")
    
    if check_taxa:
        check_taxa_hierarchies(bin_taxa_lr, taxa_df)

    bin_table_ = import_bin_abundances(True).round(1)

    if filter_samples:
        col_sel_bins = bin_table_.T.ix[:, sorted(col_corresp.keys())]
        new_b_cols = [col_corresp[i] for i in col_sel_bins.columns]
        col_sel_bins.columns = new_b_cols
    else:
        col_sel_bins = bin_table_.T

    if norm_type != "raw":
        normed_bins = l1l2clr_norm(col_sel_bins.T, norm_type, psct_val)
        print "Performed {} scaling on bins size {}".format(norm_type, normed_bins.shape)
    else:
        normed_bins = col_sel_bins.T
        print "Performed no scaling on bins size {}".format(norm_type, normed_bins.shape)

    bin_df = append_indv_col_taxa_to_df(bin_taxa_lr, normed_bins.T)

    if write_bool:
        bin_data_fname = "bin_abund_tax_test.tsv"
        bin_df.to_csv("../Data/16S_Info/"+bin_data_fname, sep="\t", header=True, index=True)

    return bin_df