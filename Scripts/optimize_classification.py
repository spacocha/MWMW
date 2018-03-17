from create_classification_data import import_dist_matrices
from analyze_classification_data import measure_separation
import numpy as np
import pandas as pd
from itertools import product
import matplotlib.pyplot as plt

dist_types = ['braycurtis', 'cityblock', 'correlation', 'cosine', 'euclidean', 'mahalanobis', 'seuclidean']
norm_choices = ["CLR", "UNIT", "RAW"]

dist_list, norm_list = zip(*[(i, j) for i, j in product(dist_types, norm_choices)])
trial_no = len(dist_list)
first_two = np.array([list(dist_list), list(norm_list)]).T
last_three = np.zeros((trial_no, 3))
data_ = np.hstack((first_two, last_three))
data_cols = ['dist_met', 'norm_choice', "mean_dist", "median_dist", "pval", "sumRatio"]
results = pd.DataFrame(index=np.arange(trial_no)+1, columns=data_cols, data=data_)
idxd_results = results.set_index(["dist_met", "norm_choice"])

for nc in norm_choices:
    ordered_mats = import_dist_matrices(nc, write_bool=False)
    for dt in dist_types:
        separations = measure_separation(dt, matrices_in_order=ordered_mats, return_data=False)
        idxd_results.loc[dt].loc[nc] = np.array(list(separations))

# braycurtis & UNIT
#raw_mats = import_dist_matrices("RAW", write_bool=False)
#raw_dists = measure_separation("euclidean", matrices_in_order=raw_mats, return_data=True)
# correlation & UNIT
#clr_mats = import_dist_matrices("CLR", write_bool=False)
#clr_dists = measure_separation("cityblock", matrices_in_order=clr_mats, return_data=True)
# cosine & UNIT
unit_mats = import_dist_matrices("UNIT", write_bool=False)
dist1 = measure_separation("braycurtis", matrices_in_order=unit_mats, return_data=True)
dist2 = measure_separation("correlation", matrices_in_order=unit_mats, return_data=True)
dist3 = measure_separation("cosine", matrices_in_order=unit_mats, return_data=True)

# mabye also  / euclidean & CLR / cosine & raw

X = [dist1, dist2, dist3]
label_list = ['Unit+BrayCurtis', 'Unit+Correlation', 'Unit+Cosine']
class_list = ["Matched", "Random"]

fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(12,8))

for ax, cnt in zip(axes.ravel(), range(3)):
    matched_dists, null_model = X[cnt]
    all_dists = np.hstack((matched_dists.values, null_model.values))
    min_b = np.floor(np.min(all_dists))
    max_b = np.ceil(np.max(all_dists))
    bins = np.linspace(min_b, max_b, 25)
    for lab, col in zip([0,1], ('red', 'green')):
        ax.hist(X[cnt][lab].values, color=col, label=class_list[lab], bins=bins, alpha=0.5)
    ylims = ax.get_ylim()
    # plot annotation
    leg = ax.legend(loc='upper right', fancybox=True, fontsize=8)
    leg.get_frame().set_alpha(0.5)
    ax.set_ylim([0, 150])
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
fig.savefig("../Data/16S_Info/best_distance_and_norming.png")
