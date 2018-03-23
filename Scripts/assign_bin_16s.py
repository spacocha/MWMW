from load_classif_data import *
from analyze_classif_data import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
import os
from itertools import chain
f_bin_df1 = import_bin_data('l1', 'rc', False, False)
filt_mgOTUs = load_mgOTU_data(filtered_data=True, norm_type='l1', norm_axes='rc', psct_val=None, check_taxa=False)
filt_tol1, splitN1 = 0.019, 17
f_mdist_df, filt_matches, f_rdist_df = join_bins_and_tags(f_bin_df1, filt_mgOTUs, filt_tol1, splitN1, True)
test_bins = ["Bin 18", "Bin 26", "Bin 52", "Bin 19", "Bin 31", "Bin 32", "Bin 55", "Bin 9"]
filt_matches.ix[test_bins, :].T
bin_df_matched = import_bin_data('l1', 'rc', True, False)
otu_df_matched = import_amplicon_matrix('l1', 'rc', False, None, None)
filt_tol2, splitN2 = 0.05, 11
otu_mdist_df, otu_matches, otu_rdist_df = join_bins_and_tags(bin_df_matched, otu_df_matched, filt_tol2, splitN2, True)
pos_bins = ["Bin 18", "Bin 26", "Bin 52"]
otu_matches.ix[pos_bins, :].T
f_bin_df2 = import_bin_data('l1', 'rc', False, False)
ufilt_mgOTUs = load_mgOTU_data(filtered_data=False, norm_type='l1', norm_axes='rc', psct_val=None, check_taxa=False)
uf_mdist_df, ufilt_matches, uf_rdist_df = join_bins_and_tags(f_bin_df2, ufilt_mgOTUs, filt_tol2, splitN1, True)
ufilt_matches.ix[pos_bins, :].T

summary_fn = "summary_of_filtered_matches.tsv"
rep_fn = "abundances_and_taxa_of_matches.tsv"
d_fn_raw = "raw_distances_of_abundances.tsv"
d_fn_mod = "mod_distances_of_tax_and_abund.tsv"

write_dir = "../Data/Metagenomic_OTUs/filtered_results"
filename_list = [summary_fn, rep_fn,  d_fn_raw, d_fn_mod]
dist_mat_list = [f_raw_dist_df, mod_dist_df]
make_and_write_reports(filt_mgOTUs, bin_df, filt_matches, write_dir, filename_list, dist_mat_list)


"""
>bin_19_GGINHOFJ_01299_seq539_seq39
>bin_31_GEBKPEPI_03035_seq412_seq110
>bin_32_HMLACPEJ_00762_seq401_seq359
>bin_55_JHAHOKBP_00807_seq595_seq144
>bin_9_MCICBKJM_02750_seq28_seq99

X = [dist_df.ix[lab, :].values for lab in label_list]

fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(13,8))

for ax, cnt, col in zip(axes.ravel(), range(3), list(('red', 'green', 'blue'))):
    matched_dists = X[cnt]
    min_b = np.floor(np.min(matched_dists))
    max_b = np.ceil(np.max(matched_dists))
    bins = np.linspace(min_b, max_b, 80)
    ax.hist(matched_dists, color=col, bins=bins, alpha=0.5)
    ylims = ax.get_ylim()
    ax.set_ylim([0, 20])
    ax.set_xlabel('Distance')
    ax.set_title(label_list[cnt])
    # hide axis ticks
    ax.tick_params(axis="both", which="both", bottom="off", top="off", 
                   labelbottom="on", left="off", right="off", labelleft="on")
    # remove axis spines
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)

axes[2].set_ylabel('count')
fig.tight_layout()
fig.savefig("../Data/16S_Info/pos_control_distance_hists.png")
"""


