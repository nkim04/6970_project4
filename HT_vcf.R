library(rpart)
library(rpart.plot)
library(rattle)
library(caret)
library(pROC)
library(randomForest)

counts <- read.csv("counts.csv")
rownames(counts) <- counts[,1]
counts <- counts[,-1]
SNP_freq <- apply(counts[,-1], 2, function(x) {sum(x)/(2*nrow(counts))})
counts <- counts[, c("Population",names(SNP_freq[SNP_freq < 0.999]))]
### Decision Trees and Bagging
set.seed(4242) 

train.index <- sample(1:nrow(counts), round(0.70*nrow(counts),0))

B <- 300

n <- length(train.index)

L_tree <-  list()
for(j in 1:300){
  idx = sample(train.index, size=n, replace=TRUE)
  L_tree[[j]] = rpart(as.factor(Population)~., counts[idx,],minsplit=1 , cp=0, method = "class")
}

#Example Decision Trees
par(mfrow=c(1,2))
fancyRpartPlot(L_tree[[30]], main = "Bagging Decision Tree: tree 30", sub = "")
fancyRpartPlot(L_tree[[275]], main = "Bagging Decision Tree: tree 275", sub = "")

## Average results
#predictions on test data for all trees
predict_test_dat <- function(i){
  predict(L_tree[[i]],newdata=counts[-train.index,],type="prob")[,2]}

out <- lapply(1:300,predict_test_dat)

pred_mat <- matrix(unlist(out),ncol=length(out))

dim(pred_mat) ## 302 observations in the test set!

pred_test <- apply(pred_mat,1,mean)
table(counts[-train.index,]$Population, round(pred_test))
## Predictive Performance

bag_roc <- roc(counts$Population[-train.index],pred_test)
bag_roc

#Single tree for comparison
train.control <- trainControl(method = "repeatedcv",
                              number = 10, ## 10-fold CV
                              repeats = 3,## repeated three times
                              summaryFunction = twoClassSummary, 
                              classProbs = TRUE)


rpartFit <- train(Population ~ ., data = counts[train.index,], 
                               method = "rpart2", 
                               tuneLength = 10,
                               trControl = train.control,
                               metric = "ROC"
)
par(mfrow=c(1,1))
fancyRpartPlot(rpartFit$finalModel, main = "Single Decision Tree for European and East Asian Superpopulations", sub = "Using ALDH2, CREB1, OCA2 and SLC45A2 SNPs")  

testPrdsTree <- predict(rpartFit$finalModel, counts[-train.index,],type = "class")
table(counts[-train.index,]$Population, testPrdsTree) 

single_roc <- roc(counts$Population[-train.index],as.numeric(testPrdsTree))
single_roc
plot(single_roc, main = "AUC for Single Decision Tree")
single_roc$auc

#Compare AUC Plots
par(mfrow=c(1,1))
plot(bag_roc, main = "AUC for Bagging Decision Trees")
bag_roc$auc 
