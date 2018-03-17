
"""
sA=../Data/16S_Info/dbOTUbySample.tsv
oMF=../Data/16S_Info/unique.dbOTU.nonchimera.mat.rdp.local.tsv
bAF=../Data/Bin_Abundances/bin_counts_normed.tsv
bTF=../Data/intermediate_files/bin_taxonomy_edit.tsv
oTF=../Data/16S_Info/unique.dbOTU.nonchimera.ns.fasta.rdp2.txt

python create_classification_data.py CLR $sA $oMF $bAF $bTF $oTF
python create_classification_data.py UNIT $sA $oMF $bAF $bTF $oTF
"""

import time
import numpy as np
import pandas as pd
import sys

drop_zero_cols = lambda df: df.T[df.sum() != 0].T

def parseBiosample(df):
    """
    This function accepts a df with samples in rows and OTUs in columns
    and parses the biosample key 
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
        
        dates_.append(date)
        primers_.append(primer)
        kits_.append(kit)
        replicates_.append(replicate)
        depths_.append(depth)
    
    df['date'] = dates_
    df['primers'] = primers_
    df['kit'] = kits_
    df['replicates'] = replicates_
    df['depth'] = depths_
    return df

def append_indv_col_taxa_to_df(taxa_df, fulldf):
    """
    takes a Series of seq's as index and RDP strings as values
    takes a df with a matching seq index
    adds 6 columns of the classification hierarchy extracted from strings
    (species column is dropped)
    """
    assert list(taxa_df.index) == list(fulldf.index)
    fullerdf = pd.concat([fulldf.copy(), taxa_df], axis=1, verify_integrity=True)
    return fullerdf

def combine_replicate_pairs(replicate_trios, old_df, comp_type):
    """
    Scales and sums replicates in old_df, returning a new df with all replicates merged
    Accepts a list of tuples containing index of new row and indexes of the two rows in old df to combine 
    """
    roots, p1s, p2s = zip(*replicate_trios)
    new_df = pd.DataFrame(index=roots, columns=old_df.columns)
    for root, p1, p2 in replicate_trios:
        if p1 in old_df.index and p2 in old_df.index:
            sub_old_df = old_df.ix[[p1, p2], :] / 2.0
            print root, "pearson corr:", sub_old_df.T.corr().ix[p1, p2]
        elif p1 in old_df.index:
            sub_old_df = old_df.ix[[p1], :]
        elif p2 in old_df.index:
            sub_old_df = old_df.ix[[p2], :]
        else:
            sys.exit("illegal column detected")
        if comp_type == "Raw":
            new_df.ix[root, :] = sub_old_df.div(sub_old_df.sum(1), 0).sum()
        elif comp_type == "Normalized":
            new_df.ix[root, :] = sub_old_df.sum()
        else:
            sys.exit("third argument must be `Raw` or `Normalized`")
    if comp_type == 'Raw':
        rescale_values = (1.0 / new_df[new_df > 0].min(1))
        rescaled_df = new_df.mul(rescale_values, 0)
        return rescaled_df
    else:
        return new_df

def unit_normalize(df):
    mat = df.values
    nmat = mat/mat.sum(1)[:, np.newaxis]
    return pd.DataFrame(index=df.index, columns=df.columns, data=nmat)

def clr(df, psct):
    """
    Accepts and returns a dataframe with all numeric values
    Performs center-log ratio transform by row
    """
    mat = df.values + psct
    mat = (mat / mat.sum(axis=1, keepdims=True)).squeeze()
    lmat = np.log(mat.astype(float))
    gm = lmat.mean(axis=-1, keepdims=True)
    data_ = (lmat - gm).squeeze()
    return pd.DataFrame(index=df.index, columns=df.columns, data=data_)

def time_k(old_stamp, deed_done):
    """
    Accepts and returns a time stamp
    prints the time passed since the old one and a string
    Returns a new time stamp
    """
    time_one_stamp = time.time()
    time_one = time_one_stamp - old_stamp
    print deed_done+(" (%05.3f s)" % (time_one))
    return time_one_stamp

def drop_low_conf_and_IS(label, confi_label, cutoff, df):
    df_copy = df.copy()
    df_bool = df.ix[:,confi_label] < float(cutoff)
    print df_bool.sum(), "{}'s will be erased for low cutoff".format(label)
    df_copy.ix[df_bool, label] = ""
    df_bool2 = df_copy.ix[:, label].str.contains("Incertae")
    print df_bool2.sum(), "{}'s will be erased for Incertae Sedis".format(label)
    df_copy.ix[df_bool2, label] = ""
    return df_copy

shotgun_abunds = "../Data/16S_Info/dbOTUbySample.tsv"
otu_matrix_file = "../Data/16S_Info/unique.dbOTU.nonchimera.mat.rdp.local.tsv"
bin_abund_file = "../Data/Bin_Abundances/bin_abundance.tab"
bin_normalizer_f = "../Data/Bin_Abundances/sample_read_count.tab"
bin_taxa_file = "../Data/intermediate_files/bin_taxonomy_edit.tsv"
otu_taxa_file = "../Data/16S_Info/unique.dbOTU.nonchimera.ns.fasta.rdp2.txt"
# correspondence between sample names across data sets
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
print "Selecting data assigned to {}".format(manual_date)

def import_bin_table():
    abund_df = pd.read_csv(bin_abund_file, sep="\t")
    sorted_cols = sorted(list(abund_df.columns)[1:])
    fix_bin_names = lambda x: x.replace("b", "B").replace(".", " ")
    abund_df.ix[:, "Genomic bins"] = abund_df.ix[:, "Genomic bins"].apply(fix_bin_names)
    abund_df.ix[:, "Bin_num"] = abund_df.ix[:, "Genomic bins"].apply(lambda x: int(x.split()[-1]))
    sorted_df = abund_df.sort_values(["Bin_num"]).set_index("Genomic bins").ix[:, sorted_cols].T
    denom_df = pd.read_csv(bin_normalizer_f, sep="\t", index_col=0).sort_index()
    assert (sorted_df.index == denom_df.index).sum() == len(denom_df.index)
    normed_df = pd.DataFrame(index=sorted_df.index, 
                             columns=sorted_df.columns, 
                             data=sorted_df.values/denom_df.values)
    normed_scaled_df = normed_df * 1e6
    return normed_scaled_df

smg_table = pd.read_csv(shotgun_abunds, sep="\t", index_col=0)

# split off taxa series and parse it into a dataframe
bin_taxa_raw = pd.read_csv(bin_taxa_file, sep="\t", index_col=0)

def clean_rdp_txt(otu_taxa_file):
    with open(otu_taxa_file, "r") as fh:
        raw_taxa_otu = [ i for i in fh.read().replace("%",".").replace('"', "").split("\n") if i !=""][6:]

    taxa_cols = raw_taxa_otu[0].split(";")
    raw_taxa_otu.remove(raw_taxa_otu[0])
    taxa_splits = [i.split(";") for i in raw_taxa_otu]
    tax_arr = np.array(taxa_splits)
    taxa_df_raw = pd.DataFrame(index=tax_arr[:,0], columns = taxa_cols[2:], data=tax_arr[:, 2:])
    taxa_df_raw.sort_index(inplace=True)
    taxa_df_dt = taxa_df_raw.apply(pd.to_numeric, errors='ignore')
    taxa_df_gd = drop_low_conf_and_IS('Genus', 'G_con', 50., taxa_df_dt)
    taxa_df_fd = drop_low_conf_and_IS("Family",  "F_con", 50., taxa_df_gd)
    taxa_df_od = drop_low_conf_and_IS("Order",  "O_con", 50., taxa_df_fd)
    taxa_df_cd = drop_low_conf_and_IS("Class",  "C_con", 50., taxa_df_od)
    taxa_df_pd = drop_low_conf_and_IS("Phylum",  "P_con", 50., taxa_df_cd)
    taxa_df_kd = drop_low_conf_and_IS("Kingdom",  "K_con", 50., taxa_df_pd)

    # drop blank species column
    taxa_df = taxa_df_kd.drop([i for i in taxa_df_kd.columns if "_con" in i ], axis=1)
    return taxa_df

good_cols = [i for i in bin_taxa_raw.columns if "%" not in i]
good_cols.remove('Species')
bin_taxa_lr = bin_taxa_raw.ix[:, good_cols].fillna("")

# make sets of both classifications at each hierarchy
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

def import_dist_matrices(norm_choice, write_bool, bin_or_train="train"):
    norm_type = norm_choice
    time_one = time.time()
    taxa_df = clean_rdp_txt(otu_taxa_file)
    check_taxa_hierarchies(bin_taxa_lr, taxa_df)
    # need to add this
    # consider: $MARCCDIR/work/sprehei1/Keith_Files/kaw/xMetaPipeline/bin/taxFormatter.py
    # consider: $MARCCDIR/work/sprehei1/Keith_Files/kaw/xMetaPipeline/data/centrifuge_*2.tsv
    # consider: $MARCCDIR/work/sprehei1/Keith_Files/kaw/xMetaPipeline/data/hiseq_ko_taxonomy2.txt
    # creates a list of root samples to create and pairs to merge to create them
    if bin_or_train == 'bin':
        replicate_trios = [(i, i+"_1.sam", i+"_2.sam") for i in sorted(smg_table.bSample.unique())]
        print replicate_trios[0]
        print replicate_trios[-1]
    else:
        replicate_trios = [(i, i+"_1.sam", i+"_2.sam") for i in sorted(col_corresp.keys())]
    # drops any metadata leaving only abundances
    just_tags = smg_table.ix[:, smg_table.columns[2:]]
    # sums scaled replicates and rescales
    combo_df = combine_replicate_pairs(replicate_trios, just_tags, "Raw")
    time_two = time_k(time_one, "Combined replicates {}")
    # drops non-matching samples and flips and renames matching ones
    sample_select = drop_zero_cols(combo_df).T
    if bin_or_train != 'bin':
        new_cols = [col_corresp[i] for i in sample_select.columns]
        sample_select.columns = new_cols
        
    # creates a list of all samples that survived culling
    test_label_set = set(list(sample_select[sample_select.sum(1) != 0].index))
    print "{} otus detected in shotgun libraries".format(len(test_label_set))
    # performs CLR on test after reflipping
    sample_select_t = sample_select.T
    if norm_type == "UNIT":
        clr_tags = unit_normalize(sample_select_t)
        time_three = time_k(time_two, "Performed Unit scaling on test size {}".format(clr_tags.shape))
    elif norm_type == "CLR":
        clr_tags = clr(sample_select_t, 0.001)
        time_three = time_k(time_two, "Performed CLR on test size {}".format(clr_tags.shape))
    elif norm_type == "RAW":
        clr_tags = sample_select_t
        time_three = time_k(time_two, "Skipped scaling on test size {}".format(clr_tags.shape))

    if norm_type == "CLR":
        off_target = samps_by_feats[samps_by_feats.depth.isin(squishy_depths)].index
        only_depths = samps_by_feats.drop(off_target)
        # filter rows by date and drop metadata columns
        aug_abund_df = only_depths[only_depths.date == manual_date]
        row_set_abunds = aug_abund_df.ix[:, aug_abund_df.columns[:-5]]
        combo_train_df = combine_replicate_pairs(putative_rep_trios, row_set_abunds, "Raw")
        time_eight = time_k(time_three, "Combined replicates")
        combo_train = drop_zero_cols(combo_train_df)
        train_label_set = set(list(combo_train.T[combo_train.sum(0) != 0].index))
        if not test_label_set.issubset(train_label_set):
            print "Train labels are a subset of test labels: {}".format(test_label_set.issubset(train_label_set))
            lost_labels = test_label_set - train_label_set
            print "These {} labels not found:\n{}".format(len(lost_labels), lost_labels)
        combo_trained = clr(combo_train, 0.001)
        t_twelve = time_k(time_eight, "Performed CLR on train of size {}".format(combo_train.shape))
    elif norm_type == "UNIT":
        aug_abund_df = samps_by_feats[samps_by_feats.date == manual_date]
        row_set_abunds = aug_abund_df.ix[:, aug_abund_df.columns[:-5]]
        norm_train_df = unit_normalize(row_set_abunds)
        off_target = norm_train_df[aug_abund_df.depth.isin(squishy_depths)].index
        only_depths = norm_train_df.drop(off_target)
        combo_train_df = combine_replicate_pairs(putative_rep_trios, only_depths, "Normalized")
        time_eight = time_k(time_three, "Combined replicates")
        combo_trained = drop_zero_cols(combo_train_df)
        train_label_set = set(list(combo_trained.T[combo_trained.sum(0) != 0].index))
        if not test_label_set.issubset(train_label_set):
            print "Train labels are a subset of test labels: {}".format(test_label_set.issubset(train_label_set))
            lost_labels = test_label_set - train_label_set
            print "These {} labels not found:\n{}".format(len(lost_labels), lost_labels)
        t_twelve = time_k(time_eight, "Performed Unit normalization on train of size {}".format(combo_trained.shape))
    elif norm_type == "RAW":
        aug_abund_df = samps_by_feats[samps_by_feats.date == manual_date]
        row_set_abunds = aug_abund_df.ix[:, aug_abund_df.columns[:-5]]
        norm_train_df = row_set_abunds
        off_target = norm_train_df[aug_abund_df.depth.isin(squishy_depths)].index
        only_depths = norm_train_df.drop(off_target)
        combo_train_df = combine_replicate_pairs(putative_rep_trios, only_depths, "Raw")
        time_eight = time_k(time_three, "Combined replicates")
        combo_trained = drop_zero_cols(combo_train_df)
        train_label_set = set(list(combo_trained.T[combo_trained.sum(0) != 0].index))
        if not test_label_set.issubset(train_label_set):
            print "Train labels are a subset of test labels: {}".format(test_label_set.issubset(train_label_set))
            lost_labels = test_label_set - train_label_set
            print "These {} labels not found:\n{}".format(len(lost_labels), lost_labels)
        t_twelve = time_k(time_eight, "Skipped normalization on train of size {}".format(combo_trained.shape))

    # create train tsv file after appending taxa series
    taxa_df_train = taxa_df.ix[list(combo_trained.T.index), :]
    train_df = append_indv_col_taxa_to_df(taxa_df_train, combo_trained.T)

    # create test tsv file after appending taxa columns
    test_data_fname = "smg_abund_tax_test.tsv"
    taxa_df_test = taxa_df.ix[list(clr_tags.T.index), :]
    test_df = append_indv_col_taxa_to_df(taxa_df_test, clr_tags.T)

    bin_table_ = import_bin_table()
    if bin_or_train == 'bin':
        col_sel_bins = bin_table_.T
    else:
        col_sel_bins = bin_table_.T.ix[:, sorted(col_corresp.keys())]
        new_b_cols = [col_corresp[i] for i in col_sel_bins.columns]
        col_sel_bins.columns = new_b_cols
    if norm_type == "CLR":
        clr_bins = clr(col_sel_bins.T, 0.001)
        t_thirt = time_k(t_twelve, "Performed CLR on bin data ({})".format(clr_bins.shape))
    elif norm_type == "UNIT":
        clr_bins = unit_normalize(col_sel_bins.T)
        t_thirt = time_k(t_twelve, "Performed unit scaling on bin data ({})".format(clr_bins.shape))
    elif norm_type == "RAW":
        clr_bins = col_sel_bins.T
        t_thirt = time_k(t_twelve, "Performed unit scaling on bin data ({})".format(clr_bins.shape))

    bin_data_fname = "bin_abund_tax_test.tsv"
    taxa_bin_df = bin_taxa_lr.ix[list(clr_bins.T.index), :]
    bin_df = append_indv_col_taxa_to_df(taxa_bin_df, clr_bins.T)
    if write_bool:
        train_df.to_csv("../Data/16S_Info/otu_abund_tax_train.tsv", sep="\t", header=True, index=True)
        test_df.to_csv("../Data/16S_Info/"+test_data_fname, sep="\t", header=True, index=True)
        bin_df.to_csv("../Data/16S_Info/"+bin_data_fname, sep="\t", header=True, index=True)

    return (train_df, test_df, bin_df)

