import os, sys
import pandas as pd

data_dir = sys.argv[-1]
type_name = sys.argv[-2]

fs = [i for i in os.listdir(data_dir) if i.endswith(".srt.tsv")]
sfs = sorted(fs)
bins_fs = [i.split(".srt")[0] for i in sfs]

if type_name == 'ours':
    bins = [i.split(".fa")[0] for i in sfs]
elif type_name == 'theirs':
    bins = [i.split(".")[0] for i in sfs]

sps = [os.path.join(data_dir, i) for i in sfs]

def parse_containment(fn):
    with open(fn) as fh:
        best = fh.read().split("\n")[0].split("\t")[-2]
    return best

best_match = [parse_containment(i) for i in sps]

sample_sheet_f = os.path.join(data_dir, "master_containment.tsv")
with open(sample_sheet_f, "w") as ofh:
    for i, j in zip(bins_fs, best_match):
        k="_".join(i.split(".")[:-1])
        l=j.split(".")[0]
        ofh.write("{}\t{}\t{}\t{}\n".format(i, j, k, l))


fix_name= lambda x: x.split(".")[0]
fix_hash= lambda x: float(x.split("/")[0]) / float(x.split("/")[1])

def measure_specificity(path_, bin_):
    df = pd.read_csv(path_, sep="\t", header=None)
    df['ref_match_name'] = df.ix[:, 4].apply(fix_name)
    df['hash_ratios'] = df.ix[:, 1].apply(fix_hash)
    sensitivity_ = df.hash_ratios.max() / df.hash_ratios.sum()
    name_ = df.ix[0, 'ref_match_name']
    return (bin_, name_, sensitivity_)

spec_vec = [measure_specificity(this_path, this_bin) for this_path, this_bin in zip(sps, bins)]
spec_sheet_f = os.path.join(data_dir, "bin_specificity.tsv")
with open(spec_sheet_f, "w") as ofh:
    for i, j, k in spec_vec:
        ofh.write("{}\t{}\t{}\n".format(i, j, k))

