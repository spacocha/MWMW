import pandas as pd

f = "bin_stats.analyze.tsv"
with open(f, "r") as fh:
    content = fh.read().split("\n")[:-1]

def prettify_checkm(test):
	name = "Bin "+test.split("\t")[0].split(".")[1]
	data = {i.split(":")[0].replace("'","").strip(): float(i.split(":")[1]) for i in test.split("\t")[-1][1:-1].split(",")}
	return (name, data)

index_, data_ = zip(*[prettify_checkm(l) for l in content])
pd.DataFrame(index=list(index_), data=list(data_)).to_csv("checkm_stats.tsv", sep="\t")