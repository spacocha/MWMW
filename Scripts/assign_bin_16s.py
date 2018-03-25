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
filt_map = '../Data/Metagenomic_OTUs/filtered_data/wd3_unique.map'
map_df = pd.read_csv(filt_map, sep="\t", header=None, index_col=1, usecols=[0,1])
filt_match_seqs = [map_df.ix[i][0] for i in otu_match_seqs]

#f_bin_df1 = import_bin_data('l1', 'rc', False, False)
#filt_mgOTUs = load_mgOTU_data(filtered_data=True, norm_type='l1', norm_axes='rc', psct_val=None, check_taxa=False)
filt_tol1, splitN1 = 0.064, 17
#f_mdist_df, filt_matches, f_rdist_df = join_bins_and_tags(f_bin_df1, filt_mgOTUs, filt_tol1, splitN1, True)
#pos_bins = ["Bin 18", "Bin 26", "Bin 52"]

def sum_error(matches, bins, dists):
    errors = [dists.ix[j, i] - dists.ix[j, :].min() for i, j in zip(matches, bins)]
    for e, b in zip(errors, bins):
        print b, round(float(e), 3)
    return sum(errors)

#sum_error(filt_match_seqs, rrna_preloc_bins, f_mdist_df)

# filt_matches.ix[pos_bins, :].T

bin_df_matched = import_bin_data('l1', 'rc', True, False)
otu_df_matched = import_amplicon_matrix('l1', 'rc', False, None, None)
filt_tol2, splitN2 = 0.0251, 11
otu_mdist_df, otu_matches, otu_rdist_df = join_bins_and_tags(bin_df_matched, otu_df_matched, filt_tol2, splitN2, True)
print sum_error(otu_match_seqs, rrna_preloc_bins, otu_mdist_df)
#otu_matches.ix[rrna_preloc_bins, :].T

summary_fn = "summary_of_filtered_matches.tsv"
rep_fn = "abundances_and_taxa_of_matches.tsv"
d_fn_raw = "raw_distances_of_abundances.tsv"
d_fn_mod = "mod_distances_of_tax_and_abund.tsv"
write_dir = "../Data/Metagenomic_OTUs/amplicon_lib_results"

filename_list = [summary_fn, rep_fn,  d_fn_raw, d_fn_mod]
dist_mat_list = [otu_rdist_df, otu_mdist_df]
make_and_write_reports(otu_df_matched, bin_df_matched, otu_matches, write_dir, filename_list, dist_mat_list)


#filt_tol3 = 0.176
#f_bin_df2 = import_bin_data('l1', 'rc', False, False)
#ufilt_mgOTUs = load_mgOTU_data(filtered_data=False, norm_type='l1', norm_axes='rc', psct_val=None, check_taxa=False)
#uf_mdist_df, ufilt_matches, uf_rdist_df = join_bins_and_tags(f_bin_df2, ufilt_mgOTUs, filt_tol3, splitN1, True)
#unfilt_match_seqs = ["seq539", "seq412", "seq401", "seq595", "seq28"]
#sum_error(unfilt_match_seqs, rrna_preloc_bins, uf_mdist_df)
#write_dir = "../Data/Metagenomic_OTUs/unfiltered_results"

#filename_list = [summary_fn, rep_fn,  d_fn_raw, d_fn_mod]
#dist_mat_list = [uf_rdist_df, uf_mdist_df]
#make_and_write_reports(ufilt_mgOTUs, f_bin_df2, ufilt_matches, write_dir, filename_list, dist_mat_list)
# bin 55 just doesn't match very closely with seq595. The abundance ratios are fundamentally different.
# Using only sample normalization doesn't help and using different distance measures doesn't help
# there are no other sequences similar to 595 which could be absorbing counts
# bin 9 is even further away from seq28. Bin9 comes up primarily in the 'A7_9mMystic' (44%) and seq28
# is entirely absent from that sample. While seq28 primarily comes up in 'A6_7mMystic' (40%) and Bin9 only
# is seen 6% of the time in that sample. 
# To catch the matches to the positive control sample, I turned the tolerance to 0.12, but I had to adjust it
# up to 0.175 for bin 55 and would have to move it all the way up to 0.33 to catch Bin 9. 
# ufilt_matches.ix[pos_bins, :].T



"""
>bin_19_GGINHOFJ_01299_seq539_
>bin_31_GEBKPEPI_03035_seq412_
>bin_32_HMLACPEJ_00762_seq401_
>bin_55_JHAHOKBP_00807_seq595_
>bin_9_MCICBKJM_02750_seq28_

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


