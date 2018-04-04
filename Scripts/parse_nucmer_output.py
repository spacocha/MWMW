#./configure --enable-python-binding=/home-3/karoraw1@jhu.edu/scratch/miniconda2/lib/python2.7

import os
import pandas as pd
import numpy as np

coord_dir = "../Data/Thrash_Libs/coords"
my_bins_dir = "../Data/Thrash_Libs/our_bins"
their_bins_dir = "../Data/Thrash_Libs/their_bins"
sensitivity_file_1 = "../Data/Thrash_Libs/thrash_in_maxbin/bin_specificity.tsv"
sensitivity_file_2 = "../Data/Thrash_Libs/max_bin_in_thrash/bin_specificity.tsv"

fs = [i for i in os.listdir(coord_dir) if i.endswith(".txt")]
ref_bins = [i.split("_")[-3] for i in fs]
query_bins = [".".join(i.split("_")[:2]) for i in fs]
ps = [os.path.join(coord_dir, i) for i in fs]
ref_ps = [os.path.join(their_bins_dir, i+".assembled.fna") for i in ref_bins]
query_ps = [os.path.join(my_bins_dir, i+".fa") for i in query_bins]

assert sum([os.path.exists(i) for i in ref_ps]) == len(ref_ps)
assert sum([os.path.exists(i) for i in query_ps]) == len(query_ps)

def genome_size(path_):
    head_check = lambda x: ">" not in x
    with open(path_, "r") as fh:
        content = [i.strip() for i in fh.read().split("\n") if i != ""]
    return float(len("".join(filter(head_check, content))))

ref_g_size = [genome_size(i) for i in ref_ps]
query_g_size = [genome_size(i) for i in query_ps]

query_aln_len, ref_aln_len = [], []
query_frac_aln, ref_frac_aln = [], []
for idx, this_f in enumerate(ps):
    print this_f
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
        q_covd += df[df.Query_Name == that_contig].ix[:, "[LEN 2]"].astype(float).sum()
        q_total += df[df.Query_Name == that_contig].ix[:, "[LEN Q]"].astype(float).unique()[0]
    
    query_frac_aln.append(q_covd/q_total)
    ref_frac_aln.append(ref_covd/ref_total)
    query_aln_len.append(q_covd)
    ref_aln_len.append(ref_covd)

data_final = np.array([ref_bins, ref_frac_aln, ref_aln_len, ref_g_size, query_frac_aln, query_aln_len, query_g_size]).T
final_cols = ["ref_bins", "ref_frac_alnd", "ref_alnd_length", "ref_g_size", "query_frac_aln", "query_alnd_length", "query_g_size"]
final_df = pd.DataFrame(data=data_final, columns=final_cols, index=query_bins)
print final_df.head()
final_df = final_df.apply(pd.to_numeric, errors='ignore')
final_df["ref_covered"] = final_df.ref_alnd_length / final_df.ref_g_size
final_df["query_covered"] = final_df.query_alnd_length / final_df.query_g_size
print final_df.ix[:, ["ref_frac_alnd", "ref_covered"]].corr()
final_df_pre = final_df.drop(["ref_frac_alnd", "ref_alnd_length", "query_frac_aln"], axis=1)


sens_srs = pd.read_csv(sensitivity_file_1, sep="\t", header=None, index_col=1)
final_df2 = pd.concat([final_df_pre, sens_srs], axis=1, verify_integrity=True)
final_df2.columns = list(final_df2.columns)[:-2] + ['mash_match_1', 'specificity of match to our bin']
assert (final_df2.ref_bins == final_df2.mash_match_1).sum() == final_df2.shape[0]

sens_srs2 = pd.read_csv(sensitivity_file_2, sep="\t", header=None, index_col=0)
final_df3 = pd.concat([final_df2, sens_srs2], axis=1, verify_integrity=True)
final_df3.columns = list(final_df3.columns)[:-2] + ['mash_match_2', 'specificity of match to their bin']
assert (final_df3.ref_bins == final_df3.mash_match_2).sum() == final_df3.shape[0]

final_df4 = final_df3.drop(['mash_match_1', "mash_match_2"], axis=1)
final_df4.columns = ['their bin', 'their bin length', 'aligned length', 'our bin length',
                     "fraction of their bin aligned", "fraction of our bin aligned"] + list(final_df4.columns[-2:])
print final_df4.head()
final_df4.to_csv("../Data/Thrash_Libs/reciprocal_algnment.tsv", sep="\t", index=False, header=True, index_label="our bin")


print "ref >90% covered: {}".format((final_df.ref_covered > 0.9).sum())
print "ref >75% covered: {}".format((final_df.ref_covered > 0.75).sum())
print "ref <75%% covered: {}".format((final_df.ref_covered < 0.75).sum())
print "query >90% covered: {}".format((final_df.query_covered > 0.9).sum())
print "query >75% covered: {}".format((final_df.query_covered > 0.75).sum())
print "query <75%% covered: {}".format((final_df.query_covered < 0.75).sum())

