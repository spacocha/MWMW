from load_classif_data import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
import os

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

def taxa_iterator(label_list, dist_df_, bin_taxa, test_taxa_):
    matched_tags, dist_vectors = [], []
    while len(label_list) > 0:
        local_dists, local_tt = dist_df_.copy(), test_taxa_.copy()
        practice_hits = [taxa_stepper(ll, local_dists, bin_taxa, local_tt, 0.019) for ll in label_list]
        bin_matched, tags_matched, weight, this_bin_matchd, best_score, diff = sorted(practice_hits, key=lambda k: k[2], reverse=True)[0]
        dist_vectors.append(this_bin_matchd)
        label_list.remove(bin_matched)
        local_dists.drop(tags_matched, axis=1, inplace=True)
        local_tt.drop(tags_matched, axis=0, inplace=True)
        matched_tags.append({"Bin":bin_matched, "Tag":tags_matched, "Weight":weight, "Min": best_score, "Diff": diff})
        print "{} matched to {} tags, {} iterations remaining".format(bin_matched, tags_matched, len(label_list))
    return (matched_tags, dist_vectors)


def join_bins_and_tags(bin_df, tag_df):
    bin_abunds = bin_df.ix[:, bin_df.columns[:17]]
    tag_abunds = tag_df.ix[:, tag_df.columns[:17]]
    bin_taxa = bin_df.ix[:, bin_df.columns[17:]]
    tag_taxa = tag_df.ix[:, tag_df.columns[17:]]
    mat1, mat2 = bin_abunds.values, tag_abunds.values
    super_mat = np.vstack((mat1, mat2))
    full_dist =  squareform(pdist(super_mat, metric='euclidean'))
    bin_to_tags = full_dist[:bin_abunds.shape[0], bin_abunds.shape[0]:]
    dist_df = pd.DataFrame(index=bin_abunds.index, columns=tag_abunds.index, data=bin_to_tags)
    raw_dist_df = dist_df.copy()
    label_list = bin_df.index.tolist()
    matched_tags, dist_vectors = taxa_iterator(label_list, dist_df, bin_taxa, tag_taxa)
    to_write = pd.DataFrame(matched_tags).set_index("Bin", verify_integrity=1).sort_values(["Weight"], ascending=False)
    return (dist_vectors, to_write, raw_dist_df)


#tag_df = import_matched_amplicons('l1', False, False, False, psct_val=None)

# TODO: test four transforms on this data to improve resolution

#bin_df = import_bin_data('l1', 'r', False, False, check_taxa=True, psct_val=None)
#filt_mgOTUs = load_mgOTU_data(filtered_data=True, norm_type='l1', norm_axes='r', psct_val=None, check_taxa=False)

bin_df_rc = import_bin_data('l1', 'rc', False, False, check_taxa=True, psct_val=None)
filt_mgOTUs_rc = load_mgOTU_data(filtered_data=True, norm_type='l1', norm_axes='rc', psct_val=None, check_taxa=False)

#bin_df_cr = import_bin_data('l1', 'cr', False, False, check_taxa=True, psct_val=None)


filt_dist_vecs, filt_matches, f_raw_dist_df = join_bins_and_tags(bin_df_rc, filt_mgOTUs_rc)
filt_matches["n_matches"] = filt_matches.Tag.apply(len)
mod_dist_df = pd.concat(filt_dist_vecs, axis=1, keys=[s.name for s in filt_dist_vecs]).T
print mod_dist_df.shape, f_raw_dist_df.shape

#tax_cols = unfilt_mgOTUs.columns[17:]
#abund_cols = unfilt_mgOTUs.columns[:17]
#print uf_mod_dist_df.shape, uf_raw_dist_df.shape

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

summary_fn = "summary_of_filtered_matches.tsv"
rep_fn = "abundances_and_taxa_of_matches.tsv"
d_fn_raw = "raw_distances_of_abundances.tsv"
d_fn_mod = "mod_distances_of_tax_and_abund.tsv"

write_dir = "../Data/Metagenomic_OTUs/filtered_results"
filename_list = [summary_fn, rep_fn,  d_fn_raw, d_fn_mod]
dist_mat_list = [f_raw_dist_df, mod_dist_df]
make_and_write_reports(filt_mgOTUs, bin_df, filt_matches, write_dir, filename_list, dist_mat_list)

#rep_fn = "abundances_and_taxa_of_matches.tsv"
#write_dir = "../Data/Metagenomic_OTUs"
#path_name = os.path.join(write_dir, rep_fn)
#to_write_df = filt_matches.copy()
#tag_df = filt_mgOTUs.copy()
#full_rep_ro.to_csv(path_name, sep="\t", index=1, header=1)


#write_dir2 = "../Data/Metagenomic_OTUs/unfiltered_results"
#unfilt_mgOTUs = load_mgOTU_data(filtered_data=False, norm_type='l1', psct_val=None, check_taxa=True)
#bin_df = import_bin_data('l1', False, False, check_taxa=True, psct_val=None)
#ufilt_dist_vecs, ufilt_matches, uf_raw_dist_df = join_bins_and_tags(bin_df, unfilt_mgOTUs)
#uf_mod_dist_df = pd.concat(ufilt_dist_vecs, axis=1, keys=[s.name for s in ufilt_dist_vecs]).T
#make_and_write_reports(unfilt_mgOTUs, bin_df, ufilt_matches, write_dir2, summary_fn, rep_fn, d_fn_raw, d_fn_mod, uf_raw_dist_df, uf_mod_dist_df)

"""
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


