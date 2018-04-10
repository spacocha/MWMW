# ffn match Rhodoferax ferrireducens T118
# faa match Geobacter sulfurreducens PCA

import sys, os
import pandas as pd
import numpy as np

match_dir = sys.argv[-1]
# match_dir = "../Data/KEGG_Annotations/Iron_Reducers/match_dir"
faa_fs = [os.path.join(match_dir, i) for i in sorted(os.listdir(match_dir)) if i.endswith(".faa.hits")]
ffn_fs = [os.path.join(match_dir, i) for i in sorted(os.listdir(match_dir)) if i.endswith(".ffn.hits")]

# path_ = faa_fs[1]

columns_ = ["bin", "pid", "start", "end", "length", "annotation"]

def make_hits_unique(df_):
    indexer = df_.pid.unique()
    good_idxs = []
    for idx1 in indexer:
        # subdivide by protein
        subdf = df_[df_.pid == idx1]
        a_g_i = subdf[subdf.length == subdf.length.max()].index[0]
        good_idxs.append(a_g_i)
        # only take the longest
    return df_.ix[good_idxs, :]

def open_read_split(path_):
    bin_tag = os.path.basename(path_).split(".")[0]
    if path_.endswith(".faa.hits"):
        anno_tag = "fe_red_geobacter"
    elif path_.endswith(".ffn.hits"):
        anno_tag = "fe_red_rhodoferax"
    else:
        raise ValueError("Illegal file type {}".format(path_))
    
    with open(path_, "r") as fh:
        content = [i for i in fh.read().split("\n") if i !=""]
    
    if content == []:
        return None
    else:
        data_ = np.array([[bin_tag]+i.split("\t")[:4]+[anno_tag] for i in content])
        frame_ = pd.DataFrame(data=data_, columns=columns_)
        return make_hits_unique(frame_)

parsed_blasts = [open_read_split(i) for i in faa_fs] + [open_read_split(i) for i in ffn_fs]
good_data = [x for x in parsed_blasts if x is not None]
master_df = pd.concat(good_data, ignore_index=True).drop(["start", "end", "length"], axis=1)
master_df.to_csv("../Data/Aggregate_Iron_Annots.tsv"), sep="\t")
# output bin, protein, fe_red_geobacter/fe_red_rhodoferax

