import os
import pandas as pd

# Step 0: Make a list of the models used and select classifications 

select_K_nums = ["K08738", "K14028", "K10535", "K10944", "K10945", "K10946", "K00370"
                 "K00362", "K00368", "K00376", "K02567", "K03385", "K04561", "K15864"
                 "K00399", "K00400", "K00401", "K00402", "K17222", "K17223", "K17224"
                 "K17225", "K17226", "K17227", "K00394", "K00958", "K11180", "K11181"]

metabolic_hmms = ["acetate_citrate_lyase_aclA", "acetate_citrate_lyase_aclB", 
                  "carbon_monoxide_dehydrogenase_coxM", "carbon_monoxide_dehydrogenase_coxS", 
                  "hydrazine_oxidase_hzoA", "hydrazine_synthase_hzsA", "Hydrogenase_Group_1", 
                  "Hydrogenase_Group_2a", "Hydrogenase_Group_2b", "Hydrogenase_Group_3a", 
                  "Hydrogenase_Group_3b", "Hydrogenase_Group_3c", "Hydrogenase_Group_3d", 
                  "Hydrogenase_Group_4", "methanol_dehydrogenase_pqq_xoxF_mxaF", 
                  "nitric_oxide_reductase_norB", "nitric_oxide_reductase_norC", 
                  "nitrite_oxidoreductase_nxrA", "nitrite_oxidoreductase_nxrB", 
                  "nitrite_reductase_nirS", "rubisco_form_I", "rubisco_form_II", 
                  "rubisco_form_III", "rubisco_form_II_III", "rubisco_form_IV", 
                  "sulfide_quinone_oxidoreductase_sqr", "sulfur_dioxygenase_sdo", 
                  "thiosulfate_reductase_phsA"] 

cazy_hmm = ["dbCAN-fam-HMMs.v6"]

# Step 1: Read in GFF for each bin and initialize a dict to hold the raw data

gff_dir = "/home-3/karoraw1@jhu.edu/scratch/metaWRAP_Out/Mystic_Bin_Fxns/bin_funct_annotations"
gff_files = [gff_f for gff_f in sorted(os.listdir(gff_dir)) if gff_f.endswith(".gff")]
gff_paths = [os.path.join(gff_dir, an_gff) for an_gff in gff_files]
bin_names = [gff.split(".")[0] for gff in gff_files]
assert sum([os.path.exists(i) for i in gff_paths]) == len(gff_paths)
 
# Step 2: Define a program to parse each line of a GFF for the following:
#         If it matches a Banfield HMM, just get the "locus_tag=IJDLJFPM_198476;"
#         If it matches the dbCAN HMM, also get the specific subtype e.g. : "protein motif:dbCAN-fam-HMMs.v6:GT41.hmm;"

def parse_gff_line(r):
    assert len(r) == 9
    info_string = r[-1]
    dbcan_check = "dbCAN-fam-HMMs.v6" in info_string
    mhmm_check = sum([i in info_string for i in metabolic_hmms]) > 0
    if mhmm_check or dbcan_check:
        prot_name = info_string.split("locus_tag=")[-1].split(";")[0]
    else:
        return None
    if dbcan_check:
        cazy_fam = "CAZy_"+ info_string.split("dbCAN-fam-HMMs.v6:")[-1].split(".hmm")[0]
        return (prot_name, cazy_fam)
    elif mhmm_check:
        mhmm = info_string.split("protein motif:")[-1].split(":")[0]
        return (prot_name, mhmm)

# Step 3: Open each GFF, parse each line, create a set of all annotation types, and simply collect the data
all_hit_types, all_hmm_hits = set(), list()
for b, gff_p in zip(bin_names, gff_paths):
    print "Reading {}'s annotations from:\n\t{}".format(b, gff_p)
    with open(gff_p, "r") as th:
        content = th.read().split("\n")
    recs = [c.split("\t") for c in content if c != ""]
    print "\t{} annotations detected".format(len(recs))
    parsed_annots = [parse_gff_line(r) for r in recs]
    spec_annots = filter(None, parsed_annots)
    if len(spec_annots) != 0:
        all_hit_types.update(zip(*spec_annots)[1])
    print "\tLogging {} environmental HMM hits".format(len(spec_annots))
    print "{} types of env genes detected so far".format(len(all_hit_types))
    for sa in spec_annots:
        all_hmm_hits.append([b, sa[0], sa[1]])

# Deliverable 1: a CSV with bin number, locus tag, and annotation
bulk_env_gene_df = pd.DataFrame(data=all_hmm_hits, 
                                index=range(len(all_hmm_hits)),
                                columns=["Bin", "ProteinID", "Annotation"])
bulk_fname = "All_HMM_Hits_raw.tsv"
bulk_env_gene_df.to_csv(bulk_fname, sep="\t", index=False)
print "Writing {} records total to {}".format(bulk_env_gene_df.shape[0], bulk_fname)

# Deliverable 2: a CSV with bin number and counts of each type of annotation 

from collections import Counter
print "Initializing counts to 0 for all bins"
all_hmm_counts = []
for b in bin_names:
    hmm_counts = Counter()
    for hmm_ in all_hit_types:
        hmm_counts[hmm_] = 0
    all_hmm_counts.append(hmm_counts)

for idx, b in enumerate(bin_names):
    sub_df = bulk_env_gene_df[bulk_env_gene_df.Bin == b]
    print "Adding {} hits from GFF".format(sub_df.shape[0])
    all_hmm_counts[idx].update( list( sub_df.Annotation ) ) 

annots_by_bin = pd.DataFrame(index=bin_names, columns=all_hit_types, data=all_hmm_counts)

annots_name = "HMM_counts_by_bin.tsv"
annots_by_bin.to_csv(annots_name, sep="\t")
print "Wrote HMM tallies by bin to {}".format(annots_name)
