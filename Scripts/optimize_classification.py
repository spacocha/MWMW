import pandas as pd
from load_classif_data import *
from analyze_classification_data import *
import numpy as np
from itertools import product
import matplotlib.pyplot as plt

dist_types = ['braycurtis', 'cityblock', 'correlation', 'cosine', 'euclidean', 'mahalanobis', 'seuclidean']
norm_choices = ['l1', 'l2', 'raw']+['clr'+str(i) for i in range(7)]
three_axes = [dist_types, norm_choices, norm_choices]
dist_list, test_tf, train_tf = zip(*[(i, j, k) for i, j, k in product(*three_axes)])
trial_no = len(dist_list)
first_cols = np.array([list(dist_list), list(test_tf), list(train_tf)]).T
last_cols = np.zeros((trial_no, 6))
data_ = np.hstack((first_cols, last_cols))
data_cols = ['dist_met', 'test_norm', 'trial_norm', "te_pct_val", "tr_pct_val",
             "mean_dist", "median_dist", "pval", "sumRatio"]
results = pd.DataFrame(index=np.arange(trial_no)+1, columns=data_cols, data=data_)

# extract pseudocount value from string
enter_pscts = lambda x:float(10.**(int(x[-1])*-1.))

for tn, tp in zip(data_cols[1:3], data_cols[3:5]):
    clr_pos = results.ix[:, tn].str.contains("clr")
    results.ix[clr_pos, tp] = results.ix[clr_pos, tn].apply(enter_pscts)
    results.ix[~clr_pos, tp] = None

# preload distance matrices
preload_mats = {"test":{}, "train":{}}
for nc_nc in norm_choices:
    if nc_nc.startswith("clr"):
        a_pct_val = enter_pscts(nc_nc)
        actual_nc = "clr"
    else:
        a_pct_val = None
        actual_nc = nc_nc
    preload_mats["test"][nc_nc] = import_matched_amplicons(actual_nc, True, False, False, psct_val=a_pct_val)
    preload_mats["train"][nc_nc] = import_amplicon_matrix(actual_nc, False, None, psct_val=a_pct_val)

for r_idx in results.index.tolist():
    dt, te_nc, tr_nc = results.ix[r_idx, :].tolist()[:3]
    a_test_df = preload_mats["test"][te_nc]
    a_train_df = preload_mats["train"][tr_nc]
    separations = measure_separation(dt, a_test_df, a_train_df, return_data=False)
    results.ix[r_idx, data_cols[-4:]] = np.array(list(separations))

results["neg_log_pval"] = results.pval.apply(lambda x: -1.0*np.log(x))
results['pval_rank'] = results.neg_log_pval.rank(ascending=0)
results['median_rank'] = results.median_dist.rank(ascending=1)
results['mean_rank'] = results.mean_dist.rank(ascending=1)

results['rank_sum'] = results.ix[:, ["pval_rank", "median_rank", "mean_rank"]].sum(1)
results.sort_values(['rank_sum']).head()

# braycurtis & UNIT
#raw_mats = import_dist_matrices("RAW", write_bool=False)
#raw_dists = measure_separation("euclidean", matrices_in_order=raw_mats, return_data=True)
# correlation & UNIT
#clr_mats = import_dist_matrices("CLR", write_bool=False)
#clr_dists = measure_separation("cityblock", matrices_in_order=clr_mats, return_data=True)
# cosine & UNIT
"""
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
"""