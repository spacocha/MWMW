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

train_x = train_df[,colnames(train_df)[2:18]]
train_y = train_df[,colnames(train_df)[1]]
test_x = test_df[,colnames(train_df)[2:18]]
test_y = test_df[,colnames(train_df)[1]]
#tune_df = rbind(train_df, test_df)
#train_idx = list(1:9953)
#test_idx = list(9954:10525)
#test_idxs = rep(test_idx, 10)
#train_idxs = rep(train_idx, 10)

library(randomForest)
library(caret)
set.seed(42)

tunedRF = randomForest(train_x, y=train_y,  xtest=test_x, ytest=test_y, importance=TRUE, ntree=1000)
y_hat = predict(tunedRF, newdata=data_test)
tab = table(y_hat, data_test$quality)
error = 1-sum(diag(tab))/sum(tab)

#mydata <- data
#wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
#for (i in 2:15) wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
#plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")