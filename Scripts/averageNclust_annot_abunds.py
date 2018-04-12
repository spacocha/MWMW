import matplotlib.pyplot as plt
import time
import pandas as pd
import urllib2
import numpy as np
from scipy.cluster.hierarchy import linkage
from kegg_cluster import unit_norm_row_col


# read in protein annotation lookup table
protein_fxn_df = pd.read_csv("../Data/KEGG_Annotations/bin_protein_annotation.tsv", sep="\t")

# fix wierd name I made leaving just KO
wierd_idx = protein_fxn_df[protein_fxn_df.Annotation == 'K04561--nitric_oxide_reductase_norB'].index
protein_fxn_df.ix[wierd_idx, "Annotation"] = ["K04561"]*len(wierd_idx)

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

# create column of strings to match to protein abundances
protein_fxn_df['gene_idx'] = protein_fxn_df.Bin + "_" + protein_fxn_df.ProteinID
# NOTE: each instance of "gene_idx" may be annotated multiple times

# drop annotations corresponding to proteins that were flagged as duplicates by Salmon
protein_fxn_df.drop([8814, 11137, 12395], axis=0, inplace=True)
protein_fxn_df_filt = protein_fxn_df[protein_fxn_df.Annotation.isin(kept_annots)]
protein_fxn_df_filt.reset_index(inplace=True)

# load protein abundances (excluding control columns) & normalize by sample and then by row
prot_abunds = pd.read_csv("../Data/KEGG_Annotations/gene_abundances_raw_wFe.tsv", sep="\t", usecols=range(15))
normed_abunds = unit_norm_row_col(prot_abunds.ix[:, prot_abunds.columns[1:]].T).T
normed_abunds['Gene_Sequence'] = prot_abunds.Gene_Sequence

# filter to retain only those in `kept_annots` & drop more duplicate proteins that went unmapped 
prot_abunds_filt_pre = normed_abunds[normed_abunds.Gene_Sequence.isin(protein_fxn_df_filt.gene_idx.unique())]
prot_abunds_filt_pre.drop([13292, 11125, 5915, 11093, 10488, 10672, 10772, 10914, 6247], axis=0, inplace=True)
prot_abunds_filt = prot_abunds_filt_pre.reset_index()

# append abundance to annotation df 
new_data = np.empty((len(protein_fxn_df_filt.index), len(prot_abunds_filt.columns)-2))
new_data[:] = np.nan
numeric_cols = list(prot_abunds_filt.drop(["index", "Gene_Sequence"], 1).columns)
to_append = pd.DataFrame(index=protein_fxn_df_filt.index, data=new_data, columns=sorted(numeric_cols))
for d_n in numeric_cols:
    for anno_n in protein_fxn_df_filt.index:
        this_gi = protein_fxn_df_filt.ix[anno_n, "gene_idx"]
        this_slice = prot_abunds_filt[prot_abunds_filt.Gene_Sequence == this_gi]
        try:
            assert this_slice.shape[0] == 1
        except AssertionError:
            if this_gi == 'unbinned_IJDLJFPM_332458':
                pass
            else:
                print this_slice
                assert this_slice.shape[0] == 1
        to_append.ix[anno_n, d_n] = this_slice[d_n][this_slice.index[0]]

# collapse to include only mean of process per depth, average dispersion within annotations, 

df = pd.concat([protein_fxn_df_filt, to_append], axis=1)
df_final = df.drop(['index'], axis=1)
fxn_dict1 = {i:"mean" for i in df_final.columns[4:]}
fxn_dict2 = {i:"std" for i in df_final.columns[4:]}
anno_df = df_final.groupby('Annotation').agg(fxn_dict1)
anno_df = anno_df[sorted(list(anno_df.columns))]
resid_df = df_final.groupby('Annotation').agg(fxn_dict2)
resid_df = resid_df[sorted(list(anno_df.columns))]
dispersion_srs = resid_df.sum(1)/float(resid_df.shape[1])
anno_n = df_final.groupby('Annotation').count().Bin
anno_n.name = "Unique_Copies"

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

replace_names = lambda x: " ".join([ko_names[x].split(";")[0], x]).strip()
anno_df2 = anno_df.reset_index()
anno_df2['Annotation'] = anno_df2.Annotation.apply(replace_names)
anno_df2.set_index("Annotation", inplace=True, verify_integrity=True)
print anno_df2.head()
to_cluster = anno_df2.values

row_clusters = linkage(to_cluster, method='ward', metric='euclidean')
from scipy.cluster.hierarchy import dendrogram

plt.figure(1, figsize=(8,8), dpi=1000)
row_dendr = dendrogram(row_clusters, labels=list(anno_df2.index))
plt.tight_layout()
plt.gcf()
plt.savefig("../Data/KEGG_Annotations/Averaged_Select_Annots/annotation_clusters_euclid.png", dpi=1000)

plt.figure(2, figsize=(8,8), dpi=1000)
row_clusters2 = linkage(to_cluster, method='complete', metric='correlation')
row_dendr = dendrogram(row_clusters2, labels=list(anno_df2.index))
plt.tight_layout()
plt.gcf()
plt.savefig("../Data/KEGG_Annotations/Averaged_Select_Annots/annotation_clusters_corr.png", dpi=1000)

# Write out raw data for sarah & other plotting routines
anno_df3 = anno_df.copy()
anno_df3['dispersion'] = dispersion_srs
anno_df3['unique_copies'] = anno_n
anno_df4 = anno_df3.reset_index()
anno_df4['Annotation'] = anno_df4.Annotation.apply(replace_names)
anno_df4.set_index("Annotation", inplace=True, verify_integrity=True)

anno_df4.to_csv("../Data/KEGG_Annotation/Annotation_Abunances_Dispersion_and_Counts.tsv", sep="\t", index=True, header=True)

# define anchors for each process

#'K10535': "Ammonia Oxidation (oxygen)" # hao
#'K00376': "Denitrification (Sum)"      # nosZ 
#'K11180': "Sulfate Reduction + Sulfur Oxidation (Sum)" # dsrA
#'K11181': "Sulfate Reduction + Sulfur Oxidation (Sum)" # dsrB
#"K00399": "Methanogenesis" # mcrA
#"fe_red_geobacter": "Iron Reduction"
#"fe_red_rhodoferax": "Iron Reduction"


