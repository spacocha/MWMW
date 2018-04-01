from load_classif_data import *
from analyze_classif_data import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
import os
from itertools import chain

rrna_preloc_bins = ["Bin 19", "Bin 31", "Bin 32", "Bin 55", "Bin 9"]
otu_match_seqs = ["seq39", "seq110", "seq359", "seq144", "seq99"]
pos_bins = ["Bin 18", "Bin 26", "Bin 52"]
bin_df_matched = import_bin_data('l1', 'rc', True, False)
otu_df_matched = import_amplicon_matrix('l1', 'rc', False, None, None)
filt_tol2, splitN2, splitN1 = 0.0251, 11, 17
otu_mdist_df, otu_matches, otu_rdist_df = join_bins_and_tags(bin_df_matched, otu_df_matched, filt_tol2, splitN2, True)
print sum_error(otu_match_seqs, rrna_preloc_bins, otu_mdist_df)
print otu_matches.ix[rrna_preloc_bins, :].T

summary_fn = "summary_of_filtered_matches.tsv"
rep_fn = "abundances_and_taxa_of_matches.tsv"
d_fn_raw = "raw_distances_of_abundances.tsv"
d_fn_mod = "mod_distances_of_tax_and_abund.tsv"
write_dir = "../Data/Metagenomic_OTUs/amplicon_lib_results"

filename_list = [summary_fn, rep_fn,  d_fn_raw, d_fn_mod]
dist_mat_list = [otu_rdist_df, otu_mdist_df]
make_and_write_reports(otu_df_matched, bin_df_matched, otu_matches, write_dir, filename_list, dist_mat_list)


filt_tol3 = 0.125
f_bin_df2 = import_bin_data('l1', 'rc', False, False)
ufilt_mgOTUs = load_mgOTU_data(filtered_data=False, norm_type='l1', norm_axes='rc', psct_val=None, check_taxa=False)
unfilt_match_seqs = ["seq539", "seq412", "seq401", "seq595", "seq28"]
uf_mdist_df, ufilt_matches, uf_rdist_df = join_bins_and_tags(f_bin_df2, ufilt_mgOTUs, filt_tol3, splitN1, True)
print sum_error(unfilt_match_seqs, rrna_preloc_bins, uf_mdist_df)
write_dir = "../Data/Metagenomic_OTUs/unfiltered_results"

filename_list = [summary_fn, rep_fn,  d_fn_raw, d_fn_mod]
dist_mat_list = [uf_rdist_df, uf_mdist_df]
make_and_write_reports(ufilt_mgOTUs, f_bin_df2, ufilt_matches, write_dir, filename_list, dist_mat_list)
print ufilt_matches.ix[pos_bins, :].T
