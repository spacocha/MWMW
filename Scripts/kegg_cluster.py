"""
python kegg_cluster.py -f ../Data/KEGG_Annotations/KO_Database/TO_Files
python kegg_cluster.py -t ../Data/KEGG_Annotations/KO_Database/KO_by_TO_Table.tsv
"""

import os, sys
import pandas as pd
import numpy as np

args = sys.argv
if '-f' in args:
    to_tuples = [(os.path.join(args[2], i), i[:-4]) for i in sorted(os.listdir(args[1]))]
    def read_kegg_table(TO_tuple):
        TO_path, TO_number = TO_tuple
        with open(TO_path, 'r') as fh:
            content = [i.split(":")[-1] for i in fh.read().split("\n") if i != ""]
        kos, cnts = np.unique(content, return_counts=1)
        return pd.Series(data=cnts, index=kos, dtype=np.float64, name=TO_number)

    TO_srs_lst = [read_kegg_table(to_tple) for to_tple in to_tuples]
    full_rep = pd.concat(TO_srs_lst, axis=1, keys=[s.name for s in TO_srs_lst]).T
    full_rep.to_csv("../Data/KEGG_Annotations/KO_Database/KO_by_TO_Table.tsv", sep="\t", index=True, header=True)
elif '-t' in args:
	full_rep = pd.read_csv(args[2], sep="\t", index_col=0).fillna(0)

# read in bin data
# find nearest neighbors in correlation distance
# 