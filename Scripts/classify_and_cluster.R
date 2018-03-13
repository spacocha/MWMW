# classify_and_cluster.R train_df.tsv test_df.tsv bin_df.tsv
# Written by Keith Arora-Williams
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
   write("Three arguments must be provided train, test, and bin data", stderr())
   test_file = ""
   train_file = ""
   bin_file = ""
}  else {
   train_file = args[1]
   test_file = args[2]
   bin_file = args[3]
}

#df = read.table(args[1], header=TRUE)
#num_vars = which(sapply(df, class)=="numeric")
#df_out = df[ ,num_vars]
#write.table(df_out, file=args[2], row.names=FALSE)