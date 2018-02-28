import pandas as pd

s_counts = pd.read_csv("../Data/Bin_Abundances/sample_read_count.tab", sep="\t", index_col=0).sort_index()
df2 = pd.read_csv("../Data/Bin_Abundances/bin_abundance.tab", sep="\t", index_col=0).T.sort_index().T
row_normed = df2.div(df2.sum())
norm_values = s_counts/s_counts.min()
sample_normed = row_normed.rmul(norm_values.ix[:, 'reads'], axis=1)
sample_normed.to_csv("../Data/Bin_Abundances/bin_counts_normed.tsv", sep="\t", index_label="Genomic bins")

