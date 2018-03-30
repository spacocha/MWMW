import pandas as pd
import sys
from collections import Counter

fasta_file = sys.argv[-1]+"/AllProteinsRenamed.faa"
with open(fasta_file, "r") as ff_h:
    headers = [i[1:].split()[0] for i in ff_h.read().split("\n") if i != "" and i.startswith(">")]

annotation_file = sys.argv[-1]+"/mystic_kegg_annotations.txt"
with open(annotation_file, "r") as af_h:
    annots = [i.split("\t") for i in af_h.read().split("\n") if i != ""]

real_annots = [i for i in annots if len(i) > 1]
stacked_annots = {h:{1: None, 2: None, 3: None} for h in headers}
for ra in real_annots:
    this_dict = stacked_annots[ra[0]]
    if not this_dict[1]:
        stacked_annots[ra[0]][1] = ra[1]
    elif not this_dict[2]:
        stacked_annots[ra[0]][2] = ra[1]
    elif not this_dict[3]:
        stacked_annots[ra[0]][3] = ra[1]
    else:
        sys.exit("Need more rows")

taxonomy_file = sys.argv[-1]+"/mystic_kegg_taxonomy.tsv"
with open(taxonomy_file, "r") as th:
    recs = [i.split("\t") for i in th.read().split("\n") if i != ""]

rec_dicts = {}
for r in recs:
    assert len(r) == 7
    assert r[0][:5] == "user:"

    if r[1] != "":
        prot_name = r[0][5:]
        rec_dict = {}
        if "," in r[1]:
            assert len(r[1].split(",")) < 4
            rec_dict['K-number_1'] = r[1].split(",")[0]
            rec_dict['K-number_2'] = r[1].split(",")[1]
            if len(r[1].split(",")) == 3: 
                rec_dict['K-number_3'] = r[1].split(",")[2]
            else:
                rec_dict['K-number_3'] = ""
        else:
            rec_dict['K-number_1'] = r[1]
            rec_dict['K-number_2'] = ""
            rec_dict['K-number_3'] = ""

        rec_dict['Kingdom'] = r[2]
        rec_dict['Class'] = r[3]
        rec_dict['Genus'] = r[4]
        rec_dict['KEGG GenesID'] = r[5]
        rec_dict['GHOSTX score'] =  r[6]
        rec_dicts[prot_name] = rec_dict

master_dicts = []

for h in headers:
    master_dict = {}
    master_dict["Protein_Name"] = h

    if h[:3] == "bin":
        master_dict["Bin_Name"] = "_".join(h.split("_")[:2])
    elif h[:8] == "unbinned":
        master_dict["Bin_Name"] = 'unbinned'
    else:
        sys.exit("illegal name")
    
    master_dict["KO_Annot_1"] = stacked_annots[h][1]
    master_dict["KO_Annot_2"] = stacked_annots[h][2]
    master_dict["KO_Annot_3"] = stacked_annots[h][3]

    if rec_dicts.has_key(h):
        master_dict['K_number_1'] = rec_dicts[h]['K-number_1']
        master_dict['K_number_2'] = rec_dicts[h]['K-number_2']
        master_dict['K_number_3'] = rec_dicts[h]['K-number_3']
        master_dict['Kingdom'] = rec_dicts[h]['Kingdom']
        master_dict['Class'] = rec_dicts[h]['Class']
        master_dict['Genus'] = rec_dicts[h]['Genus']
        master_dict['KEGG_GenesID'] = rec_dicts[h]['KEGG GenesID']
        master_dict['GHOSTX_score'] = float(rec_dicts[h]['GHOSTX score'])

    master_dicts.append(master_dict)

master_df = pd.DataFrame(master_dicts)
b1 = master_df.K_number_1.isnull()
b2 = master_df.KO_Annot_1.isnull()
clean_master = master_df[~(b1 & b2)]
b3 = clean_master.K_number_1.isnull()
b4 = clean_master.KO_Annot_1.isnull()

print "All annotation classification scores (n={})".format(clean_master.shape[0])
print "\tMean: {}".format(clean_master.GHOSTX_score.mean())
print "\tSTD: {}".format(clean_master.GHOSTX_score.std())
print "\tMax: {}".format(clean_master.GHOSTX_score.max())
print "\tMin: {}".format(clean_master.GHOSTX_score.min())

just_K_nums = clean_master[~b3 & b4]

print "Taxonomy File (only) K-numbers (n={})".format(just_K_nums.shape[0])
print "\tMean: {}".format(just_K_nums.GHOSTX_score.mean())
print "\tSTD: {}".format(just_K_nums.GHOSTX_score.std())
print "\tMax: {}".format(just_K_nums.GHOSTX_score.max())
print "\tMin: {}".format(just_K_nums.GHOSTX_score.min())

just_KOs = clean_master[~b4 & b3]
print "Annotation File (only) K-numbers (n={})".format(just_KOs.shape[0])

ordered_cols = ["Protein_Name", "Bin_Name", "KO_Annot_1","KO_Annot_2", "KO_Annot_3",
                "K_number_1", "K_number_2", "K_number_3", "Kingdom", "Class", "Genus", 
                "GHOSTX_score", "KEGG_GenesID"]

print "Filtering out GHOSTX scores below 100, see code for citation link"
# https://www.hindawi.com/journals/bmri/2016/8124636/

clean_master = clean_master[clean_master.GHOSTX_score > 100]

clean_master = clean_master[ordered_cols]
clean_master.to_csv(sys.argv[-1]+"/Aggregated_Annotations.tsv", sep="\t", index=False)

select_K_nums = ["K00362", "K00363", "K00368", "K00370", "K00371", "K00376", "K00394", 
                 "K00399", "K00400", "K00401", "K00402", "K00958", "K02305", "K02567", 
                 "K02568", "K03385", "K04561", "K07218", "K08738", "K10535", "K10944", 
                 "K10945", "K10946", "K11180", "K11181", "K14028", "K15864", "K16255", 
                 "K16257", "K16258", "K16259", "K17222", "K17223", "K17224", "K17225", 
                 "K17226", "K17227"]

select_bools_1 = [clean_master.KO_Annot_1 == i for i in select_K_nums]
select_bools_2 = [clean_master.KO_Annot_2 == i for i in select_K_nums]
select_bools_3 = [clean_master.KO_Annot_3 == i for i in select_K_nums]
select_bools_4 = [clean_master.K_number_1 == i for i in select_K_nums]
select_bools_5 = [clean_master.K_number_2 == i for i in select_K_nums]
select_bools_6 = [clean_master.K_number_3 == i for i in select_K_nums]

select_bools = select_bools_1 + select_bools_2 + select_bools_3 + select_bools_4 + select_bools_5 + select_bools_6

select_master = clean_master.copy()

print "Select boolean filters length = {}".format(len(select_bools))

for idx, sb in enumerate(select_bools):
    select_master["Filter_{}".format(idx)] = sb

print "Filter columns added equal {} - {}".format(select_master.shape[1],
                                                  clean_master.shape[1])

select_master = select_master[select_master.select_dtypes([bool]).any(1)]   

select_master.drop([i for i in list(select_master.columns) if "Filter_" in i], axis=1, inplace=True)

print "Proteins annotated as select classes (n={})".format(select_master.shape[0])
print "\tMean: {}".format(select_master.GHOSTX_score.mean())
print "\tSTD: {}".format(select_master.GHOSTX_score.std())
print "\tMax: {}".format(select_master.GHOSTX_score.max())
print "\tMin: {}".format(select_master.GHOSTX_score.min())

select_master.to_csv(sys.argv[-1]+"/Select_Annotations.tsv", sep="\t", index=False)

bin_numbers = sorted(list(select_master.Bin_Name.unique()))

countdown, k_counts_by_bin = select_master.shape[0], []
for a_bin in bin_numbers:
    k_counter = Counter()
    for k_num in select_K_nums:
        k_counter[k_num] = 0
    k_counts_by_bin.append(k_counter)
        
for idx, this_bin in enumerate(bin_numbers):    
    bin_df = select_master[select_master.Bin_Name == this_bin]
    countdown = countdown - bin_df.shape[0]
    print "{} has {} recs, {} remaining".format(this_bin, bin_df.shape[0], countdown)    
    class_1 = list(bin_df.K_number_1.dropna())
    class_2 = list(bin_df.K_number_2.dropna())
    class_3 = list(bin_df.K_number_3.dropna())
    print "Precount: {}".format(sum(k_counts_by_bin[idx].values()))
    k_counts_by_bin[idx].update(class_1)
    k_counts_by_bin[idx].update(class_2)
    k_counts_by_bin[idx].update(class_3)
    print "Post-count: {}".format(sum(k_counts_by_bin[idx].values()))

bin_k_df = pd.DataFrame(index=bin_numbers,data=k_counts_by_bin, columns=select_K_nums)
bin_k_df.to_csv(sys.argv[-1]+"/Select_Ks_By_Bin.tsv", sep="\t", index_label="Bin")

