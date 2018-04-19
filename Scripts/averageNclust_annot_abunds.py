from sklearn.preprocessing import MinMaxScaler
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from itertools import product
import matplotlib.pyplot as plt
import time
import pandas as pd
import urllib2
import numpy as np
from kegg_cluster import unit_norm_row_col

# read in protein annotation lookup table
protein_fxn_df_raw = pd.read_csv("../Data/KEGG_Annotations/bin_protein_annotation.tsv", sep="\t")

# fix wierd name I made, leaving just the KO version
wierd_idx = protein_fxn_df_raw[protein_fxn_df_raw.Annotation == 'K04561--nitric_oxide_reductase_norB'].index
protein_fxn_df_raw.ix[wierd_idx, "Annotation"] = ["K04561"]*len(wierd_idx)

# ~ 56 proteins got annotated as multiple things, since all the annotation categories we retain cover
# distinct functions, I just dumped them all as spurious. It is 1.07% of the annotations
prot_id_vec = protein_fxn_df_raw.ix[:, "ProteinID"]
uniq_ids, uniq_cnts = np.unique(prot_id_vec.values, return_counts=1)
suspect_annotations = uniq_ids[uniq_cnts == 2]
susp_annot = protein_fxn_df_raw.ix[:, "ProteinID"].isin(suspect_annotations) == False
protein_fxn_df = protein_fxn_df_raw[susp_annot]

# extract IDs assocated w/ iron reduction
iron_annots = [i for i in protein_fxn_df.Annotation.unique() if 'fe_red_' in i]

# extract all IDs called by carbohydrate active enzyme hmm
CAZy_annots = [i for i in protein_fxn_df.Annotation.unique() if 'CAZy_' in i]

# extract selected KO numbers
ko_annots = [i for i in protein_fxn_df.Annotation.unique() if i.startswith("K")]

# extract mixed bag of Banfield HMMs
meta_hmms = list(set(protein_fxn_df.Annotation.unique()) - set(ko_annots+CAZy_annots+iron_annots))

# create a list of all the annotations we want to look into
kept_annots = meta_hmms+ko_annots+iron_annots
print len(kept_annots), "annotations will be analyzed"

# fetch informative labels
def download_KO_names(this_ko):
    pre_url = "http://rest.kegg.jp/find/ko/"
    data = urllib2.urlopen(pre_url+this_ko).read()
    this_kos_name = data.split("\t")[1].strip()
    print data
    time.sleep(0.5)
    return this_kos_name

# create a lookup dict for KO gene names
ko_names = {i:download_KO_names(i) for i in ko_annots}

for ka in kept_annots:
    if not ka.startswith("K"):
        ko_names[ka] = ""

# create column of strings to match to protein abundances
protein_fxn_df['gene_idx'] = protein_fxn_df.Bin + "_" + protein_fxn_df.ProteinID
# NOTE: each instance of "gene_idx" may be annotated multiple times

# drop annotations corresponding to proteins that were flagged as duplicates by Salmon
#protein_fxn_df.drop([8814, 11137, 12395], axis=0, inplace=True)
protein_fxn_df_filt = protein_fxn_df[protein_fxn_df.Annotation.isin(kept_annots)]
protein_fxn_df_filt.reset_index(inplace=True)

# load protein abundances (excluding control columns) & normalize by sample and then by row
prot_abunds = pd.read_csv("../Data/KEGG_Annotations/gene_abundances_raw_wFe.tsv", sep="\t", usecols=range(15))
norm = 1
if norm == True:
    normed_abunds = unit_norm_row_col(prot_abunds.ix[:, prot_abunds.columns[1:]].T).T
    normed_abunds['Gene_Sequence'] = prot_abunds.Gene_Sequence
else:
    normed_abunds = prot_abunds

# filter to retain only those in `kept_annots` & drop more duplicate proteins that went unmapped 
prot_abunds_filt_pre = normed_abunds[normed_abunds.Gene_Sequence.isin(protein_fxn_df_filt.gene_idx.unique())]
#prot_abunds_filt_pre.drop([13292, 11125, 5915, 11093, 10488, 10672, 10772, 10914, 6247], axis=0, inplace=True)
prot_abunds_filt = prot_abunds_filt_pre.reset_index()

# ensure 1-to-1 mapping between annotation lookup and values read back from salmon output
## create an empty dataframe 
new_data = np.empty((len(protein_fxn_df_filt.index), len(prot_abunds_filt.columns)-2))
new_data[:] = np.nan
numeric_cols = list(prot_abunds_filt.drop(["index", "Gene_Sequence"], 1).columns)
to_append = pd.DataFrame(index=protein_fxn_df_filt.index, data=new_data, columns=sorted(numeric_cols))

## fill dataframe by going protein by protein through all annotations submitted to salmon
for d_n, anno_n in product(numeric_cols, protein_fxn_df_filt.index):
    this_gi = protein_fxn_df_filt.ix[anno_n, "gene_idx"]
    this_slice = prot_abunds_filt[prot_abunds_filt.Gene_Sequence == this_gi]
    try:
        # make sure multiple/missing entries are caught
        assert this_slice.shape[0] == 1
    except AssertionError:
        print this_slice.ix[this_slice.index[1], "index"]
        raise ValueError("There are wierd duplicate gene abundances, unrelated to colliding annotations")
    to_append.ix[anno_n, d_n] = this_slice[d_n][this_slice.index[0]]

assert to_append.isnull().sum().sum() == 0

# slap abundance columns alongside protein-annotation lookups
df = pd.concat([protein_fxn_df_filt, to_append], axis=1)
df_final = df.drop(['index'], axis=1)

# create a dictionary of functions to apply to specific columns (numerics)
fxn_dict1 = {i:"mean" for i in df_final.columns[4:]}
fxn_dict2 = {i:"std" for i in df_final.columns[4:]}

# apply them to columns, and then put new columns back in depth order
anno_df = df_final.groupby('Annotation').agg(fxn_dict1)
anno_df = anno_df[sorted(list(anno_df.columns))]
resid_df = df_final.groupby('Annotation').agg(fxn_dict2)
resid_df = resid_df[sorted(list(anno_df.columns))]

# quanitify the variance among unique protein sequences within each annotation per depth
dispersion_srs = resid_df.sum(1)
# note: this is relatively uniform 

# count the number of unique proteins per annotation
anno_N = df_final.groupby('Annotation').count().Bin
anno_N.name = "Unique_Copies"

# add gene names to KO labels
replace_names = lambda x: " ".join([ko_names[x].split(";")[0], x]).strip()
anno_df2 = anno_df.reset_index()
anno_df2['Annotation'] = anno_df2.Annotation.apply(replace_names)
anno_df2.set_index("Annotation", inplace=True, verify_integrity=True)
print anno_df2.head()

# cluster using euclidean distance
to_cluster = anno_df2.values
row_clusters = linkage(to_cluster, method='ward', metric='euclidean')

plt.figure(1, figsize=(8,8), dpi=1000)
row_dendr = dendrogram(row_clusters, labels=list(anno_df2.index))
plt.tight_layout()
plt.gcf()
plt.savefig("../Data/KEGG_Annotations/Averaged_Select_Annots/annotation_clusters_euclid.png", dpi=1000)

plt.figure(2, figsize=(8,8), dpi=1000)
row_clusters2 = linkage(to_cluster, method='complete', metric='correlation')
row_dendr2 = dendrogram(row_clusters2, labels=list(anno_df2.index))
plt.tight_layout()
plt.gcf()
plt.savefig("../Data/KEGG_Annotations/Averaged_Select_Annots/annotation_clusters_corr.png", dpi=1000)

# reorder rows with respect to the clustering
row_dendr2 = dendrogram(row_clusters2, labels=list(anno_df2.index), no_plot=True)
reordered_index = [anno_df2.index[i] for i in row_dendr2['leaves']]
df_rowclust = anno_df2.ix[reordered_index, :]

# plot heatmap
plt.clf()
fig = plt.figure(3, figsize=(8,8), dpi=1000)
ax_is = fig.add_subplot(111)
a_cmap = sns.cubehelix_palette(dark=0, light=1, as_cmap=True)
cax = sns.heatmap(df_rowclust, ax=ax_is, cmap=a_cmap)
plt.tight_layout()
plt.gcf()
plt.savefig("../Data/KEGG_Annotations/Averaged_Select_Annots/raw_corr_heatmap.png", dpi=1000)

# Write out raw data for sarah & other plotting routines
anno_df3 = anno_df.copy()
anno_df3['dispersion'] = dispersion_srs
anno_df3['unique_copies'] = anno_N
anno_df4 = anno_df3.reset_index()
anno_df4['Annotation'] = anno_df4.Annotation.apply(replace_names)
anno_df4.set_index("Annotation", inplace=True, verify_integrity=True)
data_table_fn = "../Data/KEGG_Annotations/Averaged_Select_Annots/Annotation_Abunances_Dispersion_and_Counts.tsv"
anno_df4.to_csv(data_table_fn, sep="\t", index=True, header=True)

# define anchors for each process

from collections import OrderedDict
gene_proc_map = OrderedDict()
gene_proc_map['hao K10535'] = "Ammonia Oxidation (oxygen)"
gene_proc_map['nirK K00368'] = "Ammonia Oxidation (oxygen)"
gene_proc_map['nosZ K00376'] = "Denitrification (Sum)"
gene_proc_map['nirD K00363'] = "Denitrification (Sum)"
gene_proc_map['dsrA K11180'] = "Sulfate Reduction + Sulfur Oxidation (Sum)"
gene_proc_map['dsrB K11181'] = "Sulfate Reduction + Sulfur Oxidation (Sum)"
gene_proc_map["mcrA K00399"] = "Methanogenesis"
gene_proc_map["mxaG K16255"] = "C1 oxidation (Sum)"
gene_proc_map["mxaC K16257"] = "C1 oxidation (Sum)"
gene_proc_map["mxaL K16259"] = "C1 oxidation (Sum)"
gene_proc_map["fe_red_geobacter"] = "Iron Reduction"
gene_proc_map["fe_red_rhodoferax"] = "Iron Reduction"


anno_df_final = anno_df2.ix[gene_proc_map.keys(), :]
to_cluster_final = anno_df_final.values
clean_clusters = linkage(to_cluster_final, method='complete', metric='correlation')
final_dendro = dendrogram(clean_clusters, labels=list(anno_df_final.index), no_plot=True)
final_ro_index = [anno_df_final.index[i] for i in final_dendro['leaves']]
final_rowclust = anno_df_final.ix[final_ro_index, :]

# plot heatmap
plt.clf()
fig = plt.figure(4, figsize=(8,8), dpi=1000)
ax_is = fig.add_subplot(111)
a_cmap = sns.cubehelix_palette(dark=0, light=1, as_cmap=True)
cax = sns.heatmap(final_rowclust, ax=ax_is, cmap=a_cmap)
plt.tight_layout()
plt.gcf()
plt.savefig("../Data/KEGG_Annotations/Averaged_Select_Annots/final_proc_genes_heatmap.png", dpi=1000)

fpg_fn = "../Data/calibration_data/final_proc_genes_metaWRAP.xlsx"
to_write = anno_df_final.copy()
# methanogenesis annotations are super suspect. setting them to a single flat value 
# allowing calibration to dial them  to wherever they need to be 
to_write.ix["mcrA K00399", :] *= 0
to_write.ix["mcrA K00399", :] += float(1)/14
to_write['Process'] = pd.Series(gene_proc_map)
to_write = to_write.ix[:, [to_write.columns[-1]]+ list(to_write.columns[:-1])]
to_write.to_excel(fpg_fn, sheet_name='Relative Abundances', index_label="Annotation", verbose=True)


# exact formatting for calibration
dists_only = to_write.drop(['Process', 'A1_FB_ShallowInflow'], axis=1)
real_dists_mat = dists_only.T.values
X = MinMaxScaler().fit_transform(real_dists_mat).T
real_dists_std2 = pd.DataFrame(data=X, columns=dists_only.columns, index=dists_only.index)
real_dists_std2.columns = [int(i.split("_")[1].split("m")[0]) for i in real_dists_std2.columns]
real_dists_std2['Process'] = to_write.Process
process_df_sq = real_dists_std2.groupby(by='Process').mean()
process_df = process_df_sq.unstack()
process_df = process_df.reset_index()
process_df.columns = ['Depth', 'Process', 'Value']
process_df['Date'] = pd.to_datetime(['08-12-2013']*process_df.shape[0])
process_df = process_df.set_index(['Date', 'Depth'])
process_df = process_df.pivot_table(values='Value', index=process_df.index, columns=['Process'])
process_df.index = pd.MultiIndex.from_tuples(process_df.index)
process_df.index.names = ['Date', 'Depth']
process_df.to_csv("../Data/calibration_data/gene_data.tsv", sep="\t", index=True, header=True)
process_df.to_csv("gene_data.tsv", sep="\t", index=True, header=True)
