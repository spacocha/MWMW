import os, sys
import pandas as pd

bedf = pd.read_csv(sys.argv[-1], sep="\t", header=None)

extract_mate = lambda x: x[-1]
just_header = lambda y: y[:-2]

bedf['mateN'] = bedf.ix[:,3].apply(extract_mate)
bedf['readID'] = bedf.ix[:,3].apply(just_header)
bedf['midpoint'] = (pd.to_numeric(bedf.ix[:,2]) - pd.to_numeric(bedf.ix[:,1]) )/2

aligned_headers = set(bedf.readID.tolist())
print "{} unique headers".format(len(aligned_headers))


# go through each pair and if their midpoints are within 1kb
# drop one, if not keep both

def pair_check(a_header, df):
    sub_df = df[df.readID == a_header]
    if sub_df.shape[0] == 2:
        dists = sub_df.midpoint.tolist()
        inter_dist = abs(dists[0] - dists[1])
        if inter_dist >= 1000:
            return None
        else:
            return sub_df.index[-1]
    elif sub_df.shape[0] == 1:
        return None
    else:
        sys.exit("multiple alignments detected")



returned_idxs = [pair_check(i, bedf) for i in list(aligned_headers)]
bad_idxs = [x for x in returned_idxs if x is not None]

print "{}/{} bad indexes detected".format(len(bad_idxs), len(returned_idxs))

bedf.drop(bad_idxs, axis=0, inplace=True)
bedf.drop(['mateN','readID','midpoint'], axis=1, inplace=True)
new_name = ".".join([sys.argv[-1].split(".")[0], "pe", sys.argv[-1].split(".")[-1]])
bedf.to_csv(new_name, sep="\t", header=False, index=False)

print "Wrote to {}".format(new_name)

