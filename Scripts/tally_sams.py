import os, sys
import pandas as pd
import numpy as np

target_dir = sys.argv[-1]
sam_fs = [os.path.join(target_dir, i) for i in os.listdir(target_dir) if i.endswith(".sam")]
# sam_fs = sorted(os.listdir(os.getcwd()))
def parse_sams(a_sam_f):
    print os.path.basename(a_sam_f)
    known_bollocks = set(["seq50099", 'seq178877', 'seq86424', "seq1097", "seq421", "seq120556"])
    column_names = ["Query", "Flag", "Reference", "Position_on_Ref", "MAPQ",
                    "CIGAR", "RNEXT", "PNEXT", "TLEN", "Seq", "Qual", "Alignment_Score",
                    "score_of_alt_align", "ambig_bases", "mismatches", "gap_opens", 
                    "gap_extensions", "edit distance", "mismatch_str", "pairing"]
    df = pd.read_csv(a_sam_f, sep="\t", comment="@", header=None, names=range(20))
    bad_rows = df[df.ix[:,12].str.contains("XN")].index
    df2=df.copy()
    df2.ix[bad_rows, 13:19] = df.ix[bad_rows, 12:18].values
    df2.ix[bad_rows, 12] = df.ix[bad_rows, 19].values
    df2.columns = column_names
    df2.drop(["Seq", "MAPQ"], axis=1, inplace=True)
    drop_checks = ['Flag']+list(df2.columns[4:])
    unique_types = [4, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1]
    for dc, ut in zip(drop_checks, unique_types):
        try:
            assert len(df2.ix[:, dc].unique()) == ut
        except AssertionError:
            if dc == "Flag":
                assert set(list(df2.ix[:, dc].unique())) <= set([16, 256, 272, 0])  
            elif dc == "score_of_alt_align":
            	assert len(df2.ix[:, dc].unique()) <= ut
            else:
                print dc, df2.ix[:, dc].unique()
                raise AssertionError("Flags are weird")
        if dc == "Flag":
            found_flags = set(sorted(list(df2.ix[:, "Flag"].unique())))
            try:
                assert found_flags <= set([0,16,256,272])
            except AssertionError:
                print found_flags, a_sam_f
                raise AssertionError("Flags are weird")
        df2.drop([dc], axis=1, inplace=True)
    try:
        assert len(np.unique(df2.Reference.values)) == df2.Reference.shape[0]
    except AssertionError:
        escape_flag = 0
        print len(np.unique(df2.Reference.values)), "==", df2.Reference.shape[0]
        rv_ns, rv_cs = np.unique(df2.Reference.values, return_counts=True)
        for rv_n, rv_c in zip(list(rv_ns), list(rv_cs)):
            if rv_c != 1:
                subdf = df2[df2.Reference == rv_n]
                rogue_seqs = list(subdf.Query.values)
                if len(set(rogue_seqs) & known_bollocks) > 0:
                    escape_flag = 1
                else:
                    print rogue_seqs
                    raise AssertionError("Multirec alignments are weird")
        if escape_flag == 0:
            raise AssertionError("Multirec alignments are weird")
        else:
            print "Warning: Nearly identical amplicons detected {}".format(rogue_seqs)
    amplicon_names, counts = np.unique(df2.Query.values, return_counts=True)
    sample_name = a_sam_f[:-4]
    orientation = sample_name[-1]
    amplicon_counts = {i:j for i, j in zip(amplicon_names, counts)}
    amplicon_counts['bSample'], amplicon_counts['cOrientation'],  = sample_name[:-2], orientation
    amplicon_counts["aIndex"] = a_sam_f
    return amplicon_counts

parsed_sams = [parse_sams(i) for i in sam_fs if i.endswith(".sam")]
final_df = pd.DataFrame(parsed_sams)
final_df = final_df.fillna(0)
sorted_cols = sorted(list(final_df.columns))
final_df.columns = sorted_cols
final_df.to_csv(sys.argv[-2], sep="\t", index=False)
