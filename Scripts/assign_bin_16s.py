from create_classification_data import import_dist_matrices
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform

train_df, test_df, bin_df = import_dist_matrices("UNIT", write_bool=False, bin_or_train='bin')

bin_abunds = bin_df.ix[:, bin_df.columns[:17]]
test_abunds = test_df.ix[:, bin_df.columns[:17]]

bin_taxa = bin_df.ix[:, bin_df.columns[17:]]
test_taxa = test_df.ix[:, bin_df.columns[17:]]

mat1, mat2 = bin_abunds.values, test_abunds.values
super_mat = np.vstack((mat1, mat2))
full_dist =  squareform(pdist(super_mat, metric='cosine'))

bin_to_tags = full_dist[:bin_abunds.shape[0], bin_abunds.shape[0]:]
dist_df = pd.DataFrame(index=bin_abunds.index, columns=test_abunds.index, data=bin_to_tags)

label_list = ['Bin 26', 'Bin 18', "Bin 52"]

def taxa_stepper(bin_label, dist_df, bin_taxa):
    this_bin_matches = dist_df.ix[bin_label, :].copy().round(12)
    this_bin_taxa = bin_taxa.ix[bin_label, :].copy()
    for t_label, t_name in zip(this_bin_taxa.index, this_bin_taxa.values):
        if t_name != "":
            step_1_t = test_taxa[test_taxa.ix[:, t_label] == t_name].index
            this_bin_matches[step_1_t]-=1
        else:
            break
    the_best_score = this_bin_matches.min()
    print bin_label
    best_hitters = list(this_bin_matches[this_bin_matches == the_best_score].index)
    print "\tTop", the_best_score, repr(best_hitters)
    this_bin_matches.drop(best_hitters, inplace=True)
    next_best_score = this_bin_matches.min()
    next_hitters = list(this_bin_matches[this_bin_matches == next_best_score].index)
    print "\tNext", next_best_score, repr(next_hitters)
    return (best_hitters, the_best_score, next_hitters, next_best_score)

practice_hits = [taxa_stepper(ll, dist_df, bin_taxa) for ll in label_list]

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


