# classify_and_cluster.R train_df.tsv test_df.tsv bin_df.tsv
# Written by Keith Arora-Williams
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
write("Using default train, test, and bin data", stderr())
#setwd("/Users/login/Documents/MysticLakeBins/MWMW/Scripts")
test_file = "../Data/16S_Info/smg_abund_tax_test.tsv"
train_file = "../Data/16S_Info/otu_abund_tax_train.tsv"
bin_file = "../Data/16S_Info/bin_abund_tax_test.tsv"
}  else {
   train_file = args[1]
   test_file = args[2]
   bin_file = args[3]
}

test_df = read.table(test_file, header=T, sep="\t")
train_df = read.table(train_file, header=T, sep="\t")
bin_df = read.table(bin_file, header=T, sep="\t")
colnames(test_df)[1] <- "OTU"

# https://cran.r-project.org/web/packages/PCAmixdata/vignettes/PCAmixdata.html

library(PCAmixdata)
bin_df = read.table(bin_file, header=T, sep="\t", row.names=1)
bin_df$Kingdom <- NULL
edit_df = bin_df[-c(26, 18, 52), ]

split <- splitmix(edit_df)
X1 <- split$X.quanti
X2 <- split$X.quali
res.pcamix <- PCAmix(X.quanti=X1, X.quali=X2, rename.level=TRUE, graph=FALSE)
par(mfrow=c(2,2))
plot(res.pcamix,choice="ind",coloring.ind=X2$Phylum,label=T, posleg="topright", main="Observations")
# this shows that Bins (43, 40) 4, 59, (29, 82,) (64 and 39) should be easier to identify
plot(res.pcamix,choice="cor",main="Numerical variables")
# this shows the primary axes of variation correspond to upper/lower water column & a control sample
bin_df = read.table(bin_file, header=T, sep="\t", row.names=1)
bin_df$Kingdom <- NULL
bin_df$B9_PosControlEColiMystic <- NULL
bin_df$H2_PosControlEBAlly <- NULL
bin_df$H3_PosControlEMAlly <- NULL
edit_df = bin_df[-c(26, 18, 52), ]
split <- splitmix(edit_df)
X1 <- split$X.quanti
X2 <- split$X.quali
res.pcamix2 <- PCAmix(X.quanti=X1, X.quali=X2, rename.level=TRUE, graph=FALSE)
head(res.pcamix2$eig)
plot(res.pcamix2,choice="ind",coloring.ind=X2$Phylum,label=T, posleg="topright", main="Observations")
# this shows that Bin groups above can be extended to:
# (43, 40)  (29, 82, 80) (4, 59, 78, 64 and 39) should be easier to identify
plot(res.pcamix2,choice="cor",main="Numerical variables")
plot(res.pcamix2,choice="sqload",coloring.var=T, leg=TRUE, posleg="topright", main="All variables")
set.seed(42)











