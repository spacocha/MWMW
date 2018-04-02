#./configure --enable-python-binding=/home-3/karoraw1@jhu.edu/scratch/miniconda2/lib/python2.7

import os
import pandas as pd
import numpy as np
fs = [i for i in os.listdir(os.getcwd()) if i.endswith(".txt")]
ref_bins = [i.split("_")[2][-5:] for i in fs]
query_bins = []
query_frac_aln, ref_frac_aln = [], []
for qr in fs:
    splits_ = qr.split("_")[:3]
    splits_[2] = splits_[2][:-5]
    query_bins.append(".".join(splits_))

for this_f in fs:
    with open(this_f) as fh:
        content = [i.split("\t") for i in fh.read().split("\n")[3:] if i != ""]

    content_cols = content[0][:-1]
    content_cols.extend(["Ref_Name", "Query_Name"])
    ref_covd, ref_total, q_covd, q_total = 0, 0, 0, 0
    data_ = np.array(content[1:])

    df = pd.DataFrame(data=data_, columns=content_cols)

    for this_contig in df.Ref_Name.unique():
        ref_covd += df[df.Ref_Name == this_contig].ix[:, "[LEN 1]"].astype(float).sum()
        ref_total +=  df[df.Ref_Name == this_contig].ix[:, "[LEN R]"].astype(float).unique()[0]

    for that_contig in df.Query_Name.unique():
        q_covd += df[df.Query_Name == this_contig].ix[:, "[LEN 2]"].astype(float).sum()/100.0
        q_total += df[df.Query_Name == this_contig].ix[:, "[LEN Q]"].astype(float).unique()[0]

    query_frac_aln.append(q_covd/q_total)
    ref_frac_aln.append(ref_covd/ref_total)

data_final = np.array([ref_bins, ref_frac_aln, query_bins, query_frac_aln]).T
final_cols = ["ref_bins", "ref_frac_alnd", "query_bins", "query_frac_aln"]
final_df = pd.DataFrame(data=data_final, columns=final_cols)
final_df.to_csv("reciprocal_algnment.tsv", sep="\t", index=False, header=True)


