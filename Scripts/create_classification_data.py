import time
time_start = time.time()
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

def append_indv_col_taxa_to_df(taxa_srs, fulldf):
    """
    takes a Series of seq's as index and RDP strings as values
    takes a df with a matching seq index
    adds 6 columns of the classification hierarchy extracted from strings
    (species column is dropped)
    """
    assert list(taxa_srs.index) == list(fulldf.index)
    taxa_list = list(taxa_srs.values)
    taxa_splits = [[j.split("__")[-1] for j in i.split(";")] for i in taxa_list]
    edited_splits = []
    tax_arr = np.array(taxa_splits)
    taxa_df = pd.DataFrame(index=taxa_srs.index, columns = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"], data=tax_arr)
    fullerdf = pd.concat([fulldf.copy(), taxa_df], axis=1, verify_integrity=True)
    fullerdf.drop(['Species'], axis=1, inplace=True)
    return fullerdf

def combine_replicate_pairs(replicate_trios, old_df):
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
        new_df.ix[root, :] = sub_old_df.div(sub_old_df.sum(1), 0).sum()
    rescale_values = (1.0 / new_df[new_df > 0].min(1))
    rescaled_df = new_df.mul(rescale_values, 0)
    return rescaled_df

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

time_one = time_k(time_start, "Loaded packages")

# need to add this
if len(sys.argv) < 3:
    shotgun_abunds = "../Data/16S_Info/dbOTUbySample.tsv"
    otu_matrix_file = "../Data/16S_Info/unique.dbOTU.nonchimera.mat.rdp.local.tsv"
    bin_abund_file = "../Data/Bin_Abundances/bin_counts_normed.tsv"
    bin_taxa_file = "../Data/intermediate_files/bin_taxonomy_edit.tsv"
    otu_taxa_file = "../Data/16S_Info/unique.dbOTU.nonchimera.ns.fasta.rdp2.txt"
    # need to add new taxonomy
    # need to edit it to remove propogated classifications 
    # need to remove gaps in classifications 
    

else:
    bin_abunds = sys.argv[-3]
    shotgun_abunds = sys.argv[-2]
    otu_matrix_file = sys.argv[-1]

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

bin_table = pd.read_csv(bin_abund_file, sep="\t", index_col=0)
col_sel_bins = bin_table.ix[:, sorted(col_corresp.keys())]
new_b_cols = [col_corresp[i] for i in col_sel_bins.columns]
col_sel_bins.columns = new_b_cols
bin_taxa_raw = pd.read_csv(bin_taxa_file, sep="\t", index_col=0)
good_cols = [i for i in bin_taxa_raw.columns if "%" not in i]
good_cols.remove('Species')
bin_taxa_lr = bin_taxa_raw.ix[:, good_cols].fillna("")

# reads in SMG data
smg_table = pd.read_csv(shotgun_abunds, sep="\t", index_col=0)
time_two = time_k(time_one, "Read in test (shotgun) data")
# creates a list of root samples to create and pairs to merge to create them
replicate_trios = [(i, i+"_1.sam", i+"_2.sam") for i in sorted(col_corresp.keys())]
# drops any metadata leaving only abundances
just_tags = smg_table.ix[:, smg_table.columns[2:]]
# sums scaled replicates and rescales
combo_df = combine_replicate_pairs(replicate_trios, just_tags)
time_two_half = time_k(time_two, "Combined replicates {}")
# drops non-matching samples and flips and renames matching ones
sample_select = drop_zero_cols(combo_df).T
new_cols = [col_corresp[i] for i in sample_select.columns]
sample_select.columns = new_cols
# creates a list of all samples that survived culling
test_label_set = set(list(sample_select[sample_select.sum(1) != 0].index))
print "{} otus detected in shotgun libraries".format(len(test_label_set))
# performs CLR on test after reflipping
sample_select_t = sample_select.T
clr_tags = clr(sample_select_t, 0.001)
time_three = time_k(time_two_half, "Performed CLR on test size {}".format(clr_tags.shape))

# load in test data
otu_table = pd.read_csv(otu_matrix_file, sep="\t", index_col=0)
# sort index here
time_four = time_k(time_three, "Read in train (OTU) data")

# split off taxa series and parse it into a dataframe
with open(otu_taxa_file, "r") as fh:
    raw_taxa_otu = [ i for i in fh.read().replace("%",".").replace('"', "").split("\n") if i !=""][6:]
taxa_cols = raw_taxa_otu[0].split(";")
raw_taxa_otu.remove(raw_taxa_otu[0])
taxa_splits = [i.split(";") for i in raw_taxa_otu]
tax_arr = np.array(taxa_splits)
taxa_df_raw = pd.DataFrame(index=tax_arr[:,0], columns = taxa_cols[2:], data=tax_arr[:, 2:])
# sort index here
taxa_df_dt = taxa_df_raw.apply(pd.to_numeric, errors='ignore')

def drop_low_conf_and_IS(label, confi_label, cutoff, df):
    df_copy = df.copy()
    df_bool = df.ix[:,confi_label] < float(cutoff)
    print df_bool.sum(), "{}'s will be erased for low cutoff".format(label)
    df_copy.ix[df_bool, label] = ""
    df_bool2 = df_copy.ix[:, label].str.contains("Incertae")
    print df_bool2.sum(), "{}'s will be erased for Incertae Sedis".format(label)
    df_copy.ix[df_bool2, label] = ""
    return df_copy

taxa_df_gd = drop_low_conf_and_IS('Genus', 'G_con', 50., taxa_df_dt)
taxa_df_fd = drop_low_conf_and_IS("Family",  "F_con", 50., taxa_df_gd)
taxa_df_od = drop_low_conf_and_IS("Order",  "O_con", 50., taxa_df_fd)
taxa_df_cd = drop_low_conf_and_IS("Class",  "C_con", 50., taxa_df_od)
taxa_df_pd = drop_low_conf_and_IS("Phylum",  "P_con", 50., taxa_df_cd)
taxa_df_kd = drop_low_conf_and_IS("Kingdom",  "K_con", 50., taxa_df_pd)

# drop blank species column
taxa_df = taxa_df_kd.drop([i for i in taxa_df_kd.columns if "_con" in i ], axis=1)
# make sets of both classifications at each hierarchy
bin_taxa_sets = [set(bin_taxa_lr.ix[:, i].tolist()) for i in bin_taxa_lr.columns]
otu_taxa_sets = [set(taxa_df.ix[:, i].tolist()) for i in taxa_df.columns]
taxa_is_subset = [(s.issubset(t)) for s, t in zip(bin_taxa_sets, otu_taxa_sets)]
if sum(taxa_is_subset) == len(taxa_is_subset):
    print "All bin taxa levels are in otu taxa hierarchy!"
else:
    taxa_diff = {n:(s - t) for n, s, t in zip(bin_taxa_lr.columns, bin_taxa_sets, otu_taxa_sets)}
    print taxa_diff

# drop taxa series and bad columns first
otu_time_table = otu_table.drop([otu_table.columns[-1]], axis=1)
do_not_use = ['SB100912TAWMD14VV4TMR1', 'SB061713TAWMD22VV4TMR1', 'SB011413TAWMD22VV4TMR1', 'SB011413TAWMDSBVV4TMR1', 'SB011413TAWMDEBVV4TMR1', 'SB011413TAWMD22VV4TMR2', 'SB011413TAWMDSBVV4TMR2']
otu_time_table.drop(do_not_use, axis=1, inplace=True)

# add metadata and drop undated columns
samps_by_feats = parseBiosample(otu_time_table.T)
undated = samps_by_feats[samps_by_feats.date == "NA"].index
samps_by_feats.drop(undated, inplace=True)
# parse dates
samps_by_feats['date'] = pd.to_datetime(samps_by_feats['date'])
time_seven = time_k(time_four, "Parsed biosample & dates")

# Drop samples of unknown provenance 
squishy_depths = ['bottom', 'SB', 'control1', 'control3', 'control6', 'mid1', 'mid2', 'mid3', 'river', 'upper', 'River', 'Neg', 'NA', 'EB', 'FB', 'CR', 'Dock', 'Bridge']
off_target = samps_by_feats[samps_by_feats.depth.isin(squishy_depths)].index
only_depths = samps_by_feats.drop(off_target)
manual_date = np.datetime64("2013-08-12")
# filter rows by date and drop metadata columns
print "Selecting data assigned to {}".format(manual_date)
aug_abund_df = only_depths[only_depths.date == manual_date]
print "Dropping {}".format(list(aug_abund_df.columns[-5:]))
row_set_abunds = aug_abund_df.ix[:, aug_abund_df.columns[:-5]]
putative_rep_trios = [(i, i+"R1", i+"R2") for i in sorted(col_corresp.values())]
combo_train_df = combine_replicate_pairs(putative_rep_trios, row_set_abunds)
time_eight = time_k(time_seven, "Combined replicates")
combo_train = drop_zero_cols(combo_train_df)
print combo_train.shape
train_label_set = set(list(combo_train.T[combo_train.sum(0) != 0].index))
if not test_label_set.issubset(train_label_set):
    print "Train labels are a subset of test labels: {}".format(test_label_set.issubset(train_label_set))
    lost_labels = test_label_set - train_label_set
    print "These {} labels not found:\n{}".format(len(lost_labels), lost_labels)
clr_train = clr(combo_train, 0.001)
time_nine = time_k(time_eight, "Performed CLR on train of size {}".format(combo_train.shape))

# create train tsv file after appending taxa series
taxa_series_train = taxa_series[list(clr_train.T.index)]
train_df = append_indv_col_taxa_to_df(taxa_series_train, clr_train.T)
time_ten = time_k(time_nine, "Added taxa to train")
train_df.to_csv("../Data/16S_Info/otu_abund_tax_train.tsv", sep="\t", header=True, index=True)
time_eleven = time_k(time_ten, "Wrote out train data ({})".format(train_df.shape))

# create test tsv file after appending taxa columns
test_data_fname = "smg_abund_tax_test.tsv"
taxa_series_test = taxa_series[list(clr_tags.T.index)]
test_df = append_indv_col_taxa_to_df(taxa_series_test, clr_tags.T)
test_df.to_csv("../Data/16S_Info/"+test_data_fname, sep="\t", header=True, index=True)
_ = time_k(time_eleven, "Wrote out test data ({})".format(test_df.shape))





sys.exit("Done. Exiting")