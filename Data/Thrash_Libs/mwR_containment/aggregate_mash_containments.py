import os
fs = [i for i in os.listdir(os.getcwd()) if i.endswith(".srt.tsv")]
sfs = sorted(fs)
bins = [i.split(".srt")[0] for i in sfs]
def parse_containment(fn):
    with open(fn) as fh:
        best = fh.read().split("\n")[0].split("\t")[-2]
    return best
best_match = [parse_containment(i) for i in sfs]

with open("master_containment.tsv", "w") as ofh:
    for i, j in zip(bins, best_match):
        k="_".join(i.split(".")[:-1])
        l=j.split(".")[0]
        ofh.write("{}\t{}\t{}\t{}\n".format(i, j, k, l))

