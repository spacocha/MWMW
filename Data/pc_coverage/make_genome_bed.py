import numpy as np
import pandas as pd

genome_len=6568739
p1_len=106999
p2_len=93527
p3_len=66926
p4_len=44754
p5_len=38784


ordered_regions = [genome_len, p1_len, p2_len, p3_len, p4_len, p5_len]
label_prefix = ["NC_015703.1", "NC_015693.1", "NC_015704.1", 
                "NC_015694.1", "NC_015705.1", "NC_015695.1"]

start_intervals, end_intervals = [], []
row_labels = []
for idx, o_r in enumerate(ordered_regions):
    x = np.arange(1,(o_r/1000)+2)*1000
    y = x-999
    x[-1] = o_r
    row_labels = row_labels + [label_prefix[idx]]*len(x)
    print label_prefix[idx]
    print y[0], x[0]
    print y[1], x[1]
    print y[-2], x[-2]
    print y[-1], x[-1]
    end_intervals.append(x)
    start_intervals.append(y)

print ""
end_int_arr = np.concatenate(end_intervals)
start_int_arr = np.concatenate(start_intervals)
label_arr = np.array(row_labels)
data_ = np.vstack((label_arr, start_int_arr, end_int_arr))
print end_int_arr.shape, label_arr.shape, data_.shape
print data_[:3, :]
bedf = pd.DataFrame(data_.T)
print bedf.head()
bedf.to_csv("r_slithy_genome.bed", sep="\t", index=False, header=False)
