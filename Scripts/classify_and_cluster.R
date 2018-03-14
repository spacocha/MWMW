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

library(caret)
library(randomForest)

control <- trainControl(method = "repeatedcv", number = 10, repeats = 2,
                            classProbs = FALSE, verboseIter = TRUE,
                            preProcOptions=list(na.remove=TRUE,verbose=TRUE))
tunegrid <- expand.grid(.mtry=c(1:10))

rf_gridsearch <- train(quality ~ . , data=data_train,
                           method="rf", metric="Accuracy",
                           tuneGrid=tunegrid, trControl=control)
print(rf_gridsearch)                                                  
plot(rf_gridsearch)

# Solve a randomForest with the tuned value for mtry
tunedRF = randomForest(quality ~ . , data=data, subset=train, mtry=2, importance=TRUE, ntree=1000)
y_hat = predict(tunedRF, newdata=data_test)
tab = table(y_hat, data_test$quality)
error = 1-sum(diag(tab))/sum(tab)

#mydata <- data
#wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
#for (i in 2:15) wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
#plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")