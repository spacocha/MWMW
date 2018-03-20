from load_classif_data import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import normalize

tag_df = import_matched_amplicons('l1', False, False, False, psct_val=None)
#tag_df2 = import_matched_amplicons('raw', False, False, False, psct_val=None)

# TODO: test four transforms on this data to improve resolution
bin_df = import_bin_data('l1', False, False, check_taxa=True, psct_val=None)

bin_abunds = bin_df.ix[:, bin_df.columns[:17]]
tag_abunds = tag_df.ix[:, tag_df.columns[:17]]
#tag_abunds2 = tag_df2.ix[:, tag_df2.columns[:17]]
bin_taxa = bin_df.ix[:, bin_df.columns[17:]]
tag_taxa = tag_df.ix[:, tag_df.columns[17:]]

mat1, mat2 = bin_abunds.values, tag_abunds.values
super_mat = np.vstack((mat1, mat2))
full_dist =  squareform(pdist(super_mat, metric='cosine'))

bin_to_tags = full_dist[:bin_abunds.shape[0], bin_abunds.shape[0]:]
dist_df = pd.DataFrame(index=bin_abunds.index, columns=tag_abunds.index, data=bin_to_tags)

label_list = bin_df.index.tolist()

def taxa_stepper(bin_label, distdf, bintaxa, testtaxa, sigfigs):
    this_bin_matches = distdf.ix[bin_label, :].copy().round(sigfigs)
    this_bin_taxa = bintaxa.ix[bin_label, :].copy()
    tt_cp = testtaxa.copy()
    score_denom = float(6 - list(this_bin_taxa.apply(len)).count(0))
    for t_label, t_name in zip(this_bin_taxa.index, this_bin_taxa.values):
        if t_name != "":
            step_1_t = tt_cp[tt_cp.ix[:, t_label] == t_name].index
            this_bin_matches[step_1_t]-=(1./score_denom)
        else:
            break
    the_best_score = this_bin_matches.min()
    best_hitters = list(this_bin_matches[this_bin_matches == the_best_score].index)
    this_bin_matches.drop(best_hitters, inplace=True)
    next_best_score = this_bin_matches.min()
    next_hitters = list(this_bin_matches[this_bin_matches == next_best_score].index)
    diff_score = abs(the_best_score-next_best_score)/abs(next_best_score)
    return (bin_label, best_hitters, diff_score)

def taxa_iterator(label_list, dist_df_, bin_taxa, test_taxa_):
    local_dists, local_tt = dist_df_.copy(), test_taxa_.copy()
    matched_tags = []
    while len(label_list) > 0:
        practice_hits = [taxa_stepper(ll, local_dists, bin_taxa, local_tt, 15) for ll in label_list]
        bin_matched, tags_matched, winning_score = sorted(practice_hits, key=lambda k: k[2], reverse=True)[0]
        label_list.remove(bin_matched)
        local_dists.drop(tags_matched, axis=1, inplace=True)
        local_tt.drop(tags_matched, axis=0, inplace=True)
        matched_tags.append({"Bin":bin_matched, "Tag":tags_matched, "Difference":winning_score})
        print "{} matched to {} tags, {} iterations remaining".format(bin_matched, tags_matched, len(label_list))
    return matched_tags

matched_tags = taxa_iterator(label_list, dist_df, bin_taxa, tag_taxa)
to_write = pd.DataFrame(matched_tags).set_index("Bin", verify_integrity=1).sort_values(["Difference"], ascending=False)
match_file_n = "../Data/16S_Info/bin_tag_matches.tsv"
to_write.to_csv(match_file_n, sep="\t", index_label="Bin", index=1, header=1)











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


