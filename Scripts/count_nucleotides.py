import sys
from Bio import SeqIO

counter_ = 0
with open(sys.argv[-1], "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        counter_ += len(record.seq)


print "{} bp read".format(counter_)