import os, sys
import pandas as pd

# read in sample names (everything before first period in file name)
# read in corresponding data vector
q_fs = [(i.split(".")[0], pd.read_csv(i,sep="\t", index_col=0)) for i in sorted(os.listdir(os.getcwd())) if i.endswith(".quant.counts")]

# make column names of each data vector the sample names
for i, j in q_fs:
    j.columns = [i]

# concatenate individaul columns into single frame
full_df = pd.concat(zip(*q_fs)[1], 1)

# write out to command line specified path
full_df.to_csv(sys.argv[-1]+"/gene_abundances_raw_wFe.tsv", sep="\t", index_label="Gene_Sequence")
