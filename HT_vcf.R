library(rpart)
library(rpart.plot)
library(rattle)
library(caret)
library(pROC)

gene_counts <- read.csv("counts.csv")

set.seed(1245)
test_index <- sample.int(dim(gene_counts)[1], round(dim(gene_counts)[1]*.15),replace=FALSE)

train.control <- trainControl(method = "repeatedcv",
                              number = 10, ## 10-fold CV
                              repeats = 3,## repeated three times
                              summaryFunction = twoClassSummary, 
                              classProbs = TRUE)


rpartFit <- train(Population ~ ., data = gene_counts[-test_index,-1], 
                               method = "rpart2", 
                               tuneLength = 10,
                               trControl = train.control,
                               metric = "ROC"
)

rpartFit$results
fancyRpartPlot(rpartFit$finalModel, main = "Decision Tree for European and East Asian Superpopulations", sub = "Using ALDH2, CREB1, OCA2 and SLC45A2 SNPs")  

trainPrds <- predict(rpartFit$finalModel,gene_counts[-test_index,-1],type = "class")
table(gene_counts[-test_index,-1]$Population,trainPrds)

testPrdsTree <- predict(rpartFit$finalModel, gene_counts[test_index,-1],type = "class")
table(gene_counts[test_index,-1]$Population, testPrdsTree) 

hist(predict(rpartFit$finalModel, gene_counts[test_index,-1],type = "prob")[,2] )
sn.sp <- function(mat){
  sn <- mat[2,2]/sum(mat[2,])
  sp <- mat[1,1]/sum(mat[1,])
  return(unlist(list(sensitivity=sn, specificity=sp)))
}
sn.sp(table(gene_counts[test_index,]$Population, testPrdsTree))

y <- as.integer(as.factor(gene_counts$Population))-1
auc.test_tree <- roc(y[test_index],as.numeric(testPrdsTree))
auc.test_tree
plot(auc.test_tree)
