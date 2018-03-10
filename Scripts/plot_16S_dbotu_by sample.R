setwd("/Users/login/Documents/MysticLakeBins/MWMW/Data")
otu_df = read.table("dbOTUbySample.tsv", row.names=1, header=T)
just_data = otu_df[colnames(otu_df)[3:714]]
corrplot(cor(t(just_data)), method="circle", tl.cex=0.4, type="upper", order="hclust")
