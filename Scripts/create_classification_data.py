
"""
sA=../Data/16S_Info/dbOTUbySample.tsv
oMF=../Data/16S_Info/unique.dbOTU.nonchimera.mat.rdp.local.tsv
bAF=../Data/Bin_Abundances/bin_counts_normed.tsv
bTF=../Data/intermediate_files/bin_taxonomy_edit.tsv
oTF=../Data/16S_Info/unique.dbOTU.nonchimera.ns.fasta.rdp2.txt

python create_classification_data.py CLR $sA $oMF $bAF $bTF $oTF
python create_classification_data.py UNIT $sA $oMF $bAF $bTF $oTF
"""
# make sets of both classifications at each hierarchy


def import_dist_matrices(norm_choice, write_bool, bin_or_train="train"):


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
            print "\t{} labels not found:".format(len(lost_labels))
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
            print "\t{} labels not found:".format(len(lost_labels))
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
            print "\t{} labels not found:".format(len(lost_labels))
        t_twelve = time_k(time_eight, "Skipped normalization on train of size {}".format(combo_trained.shape))

    # create train tsv file after appending taxa series
    taxa_df_train = taxa_df.ix[list(combo_trained.T.index), :]

    
    
    train_df = append_indv_col_taxa_to_df(taxa_df, combo_trained.T)
    if write_bool:
        train_fname="otu_abund_tax_train.tsv"
        train_df.to_csv("../Data/16S_Info/"+train_fname, sep="\t", header=True, index=True)

    

