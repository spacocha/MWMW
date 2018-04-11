"""
python select_gene_seqs.py ../Data/KEGG_Annotations

step 1: make a table of bin number, protein id, and select annotation type (k num or hmm) 
step 2: parse GFF to get contig number & locus information (add +20 bp on each side if possible)
step 3: parse *.fna to pull trimmed sequence
step 4: write each to single FASTA
step 5: index with salmon
step 6: map reads back 
"""

# Step 1:

import pandas as pd
import os, sys
## import
hmm_df = pd.read_csv(sys.argv[-1]+"/All_HMM_Hits_raw.tsv", sep="\t")
knum_df = pd.read_csv(sys.argv[-1]+"/Select_Annotations.tsv", sep="\t", usecols=[0, 1, 5])
iron_df = pd.read_csv(sys.argv[-1]+"/Aggregate_Iron_Annots.tsv", sep="\t", index_col=0)

## preprocess
def edit_pname(p_str):
    if 'unbinned' in p_str:
        return p_str.replace("unbinned_", "")
    elif p_str.startswith("bin_"):
        return "_".join(p_str.split("_")[2:])
    else:
        raise ValueError('A strangely formatted Protein Name was found')

## merge
knum_df["ProteinID"] = knum_df.Protein_Name.apply(edit_pname)
all_annots_df = knum_df.ix[:, ["Bin_Name", "ProteinID", "K_number_1"]].copy()
all_annots_df.columns = list(hmm_df.columns)
iron_df.columns = list(hmm_df.columns)
some_annots_df = all_annots_df.append(hmm_df, ignore_index=True)
all_annots_df = some_annots_df.append(iron_df, ignore_index=True)
## edit redundancy
dup_entries = all_annots_df[all_annots_df.ProteinID == 'IJDLJFPM_212255']
merged_annotation = "--".join(dup_entries.Annotation.tolist())
all_annots_df.ix[dup_entries.index[0], "Annotation"] = merged_annotation
all_annots_df.drop([dup_entries.index[1]], inplace=True)
print all_annots_df.head()
all_annots_df.to_csv("../Data/KEGG_Annotations/bin_protein_annotation.tsv", sep="\t", index=False)
sys.exit()
## make space for new data
all_annots_df["Start"] = pd.Series(data=[None]*all_annots_df.shape[0], index=all_annots_df.index)
all_annots_df["End"] = pd.Series(data=[None]*all_annots_df.shape[0], index=all_annots_df.index)
all_annots_df["Contig"] = pd.Series(data=[None]*all_annots_df.shape[0], index=all_annots_df.index)
all_annots_df["Sequence"] = pd.Series(data=[None]*all_annots_df.shape[0], index=all_annots_df.index)
all_annots_df["Insert_Len"] = pd.Series(data=[None]*all_annots_df.shape[0], index=all_annots_df.index)
seq_path = "/home-3/karoraw1@jhu.edu/scratch/metaWRAP_Out/Mystic_Bin_Fxns/prokka_out"
make_fasta_p = lambda x: os.path.join(seq_path, x, x+".fna")

## add sequence file locations
all_annots_df["FASTA"] =  all_annots_df.Bin.apply(make_fasta_p)
all_annots_df.reset_index(drop=True, inplace=True)

# Step 2:
gff_dir = "/home-3/karoraw1@jhu.edu/scratch/metaWRAP_Out/Mystic_Bin_Fxns/bin_funct_annotations"
gff_files = [gff_f for gff_f in sorted(os.listdir(gff_dir)) if gff_f.endswith(".gff")]
gff_paths = [os.path.join(gff_dir, an_gff) for an_gff in gff_files]
bin_names = [gff.split(".")[0] for gff in gff_files]
assert sum([os.path.exists(i) for i in gff_paths]) == len(gff_paths)

def parse_gff_line(r):
    assert len(r) == 9
    prot_name = r[-1].split("locus_tag=")[-1].split(";")[0]
    contig_name = r[0]
    start_n, end_n = int(r[3]), int(r[4])
    return (prot_name, {"cID": contig_name, "start": start_n, "end": end_n})

# b, gff_p = bin_names[0], gff_paths[0]
for b, gff_p in zip(bin_names, gff_paths):
    good_hits = all_annots_df[all_annots_df.Bin == b].index
    print "Reading {}'s annotations from:\n\t{}".format(b, gff_p)
    with open(gff_p, "r") as th:
        content = th.read().split("\n")
    recs = [c.split("\t") for c in content if c != ""]
    print "\t{} annotations detected".format(len(recs))
    parsed_annots = [parse_gff_line(r) for r in recs]
    anno_dic = {i[0]:i[1] for i in parsed_annots}
    assert len(anno_dic.keys()) == len(parsed_annots)
    for idx in good_hits:
        pID = all_annots_df.ix[idx, 'ProteinID']
        all_annots_df.ix[idx, "Start"] = anno_dic[pID]['start']
        all_annots_df.ix[idx, "End"] = anno_dic[pID]['end']
        all_annots_df.ix[idx, "Contig"] = anno_dic[pID]['cID']
        if idx == 369:
            print idx
            print all_annots_df.ix[idx, :]
bin_fnas = sorted(all_annots_df.FASTA.unique())
bin_names = [os.path.basename(i).split(".")[0] for i in bin_fnas]

def pull_sequences(f):
    sequences = {}
    with open(f, "r") as fh:
        pre_content = [pc for pc in fh.read().split("\n") if pc != ""]
    clean_lines = "\n".join([i[0] + i[1:].replace(">", "^") for i in pre_content])
    content = [l for l in clean_lines.split(">") if l != ""]
    for rec in content:
        rec_lines = rec.split("\n")
        sequences[rec_lines[0]] = "".join(rec_lines[1:])
    return sequences

def cut_gene_seq(s_, e_, seq_, pad_n):
    half_pad = pad_n/2
    est_len = ((e_-s_)+pad_n)
    if est_len < len(seq_):
        if s_ >= half_pad and ((e_+half_pad) < len(seq_)):
            return seq_[(s_-half_pad):(e_+half_pad)]
        elif s_ < half_pad and ((e_+(pad_n-s_)) < len(seq_)):
            return seq_[:(e_+(pad_n-s_))]
        elif (s_ >= (pad_n - (len(seq_)-e_)))  and ((e_+half_pad) >= len(seq_)-1):
            return seq_[(s_ - (pad_n - (len(seq_)-e_))):]
    else:
        return seq_

#b, p = bin_names[0], bin_fnas[0]
for b, p in zip(bin_names, bin_fnas):
    print "Pulling sequences for {}".format(b)
    seq_dict = pull_sequences(p)
    good_hits_ = all_annots_df[all_annots_df.Bin == b].index
    for idx in good_hits_:
        s_ = all_annots_df.ix[idx, "Start"]
        e_ = all_annots_df.ix[idx, "End"]
        cID = all_annots_df.ix[idx, 'Contig']
        seq_ = seq_dict[cID]
        est_len = e_-s_
        if est_len < 250:
            oligo = cut_gene_seq(s_, e_, seq_, 200)
        else:
            oligo = cut_gene_seq(s_, e_, seq_, 100)
        if not oligo:
            sys.exit("length check")
        ins_len = len(oligo)
        all_annots_df.ix[idx, "Sequence"] = oligo
        all_annots_df.ix[idx, "Insert_Len"] = ins_len
    print "\t Min seq length is {}".format(all_annots_df.ix[ good_hits_ , "Insert_Len"].min())

print "{} sequences were not found".format(all_annots_df.Sequence.isnull().sum())

all_annots_df.to_csv(sys.argv[-1]+"/Integrated_Annotations_with_Seqs.tsv", sep="\t", index=False)

with open(sys.argv[-1]+"/Annotated_Gene_Seqs_wFe.fa", "w") as ags_h:
    for idx in all_annots_df.index:
        header = ">"+"_".join(list(all_annots_df.ix[idx, ["Bin","ProteinID"]]))
        sequence = all_annots_df.ix[idx, "Sequence"]
        ags_h.write(header+"\n")
        ags_h.write(sequence+"\n")





#salmon index -p $threads -t $assembly -i ${out}/assembly_index

#salmon quant -i ${out}/assembly_index --libType IU -1 $reads_1 -2 $reads_2 -o ${out}/alignment_files/${sample}.quant
