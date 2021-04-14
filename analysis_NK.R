# Preliminary analyses


library(tidyverse)
library(HardyWeinberg)
library(glmnet)
library(pROC)
library(caret)


##
### DATA CONSOLIDATION-----
##


ALDH2_counts <- read_csv("ALDH2_counts.csv") %>% 
  remove_rownames() %>% 
  column_to_rownames(var="X1")


OCA2_counts <- read_csv("OCA2_counts.csv") %>% 
  remove_rownames() %>% 
  column_to_rownames(var="X1")


CREB1_counts <- read_csv("CREB1_counts.csv") %>% 
  remove_rownames() %>% 
  column_to_rownames(var="X1")


SLC45A2_counts <- read_csv("SLC45A2_counts.csv") %>% 
  remove_rownames() %>% 
  column_to_rownames(var="X1")


sum(rownames(ALDH2_counts) %in% rownames(CREB1_counts))
sum(rownames(CREB1_counts) %in% rownames(OCA2_counts))
sum(rownames(OCA2_counts) %in% rownames(SLC45A2_counts))


counts <- cbind(ALDH2_counts, 
                CREB1_counts[2:ncol(CREB1_counts)],
                OCA2_counts[2:ncol(OCA2_counts)],
                SLC45A2_counts[2:ncol(SLC45A2_counts)])


rm(list=ls())
#write.csv(counts, file = "counts.csv")
counts <- read.csv("counts.csv")


##
### HARDY WEINBERG EQUILIBRIUM TEST-----
##


# Remove non-biallelic SNPs
multi_names <- colnames(counts)[grep("\\.2", colnames(counts))] %>% 
  substr(., 1, nchar(.)-2)


multi_index <- NULL
for (i in seq(length(multi_names))) {
  multi_index <- c(multi_index, grep(multi_names[i], colnames(counts)))
}
biallelic_counts <- t(counts[, -c(1, 2, multi_index)])
  

# Get genotype counts
AA <- apply(biallelic_counts, 1, function(x) sum(x == 0)) # Homozygous reference
AB <- apply(biallelic_counts, 1, function(x) sum(x == 1)) # Heterozygous
BB <- apply(biallelic_counts, 1, function(x) sum(x == 2)) # Homozygous alternate
genotypes <- cbind(AA, AB, BB)


# Plot genotype frequencies (blue line = HWE)
HWGenotypePlot(genotypes, plottype = 1, pch = 19, cex = 0.8)
HWGenotypePlot(genotypes, plottype = 2, pch = 19, cex = 0.8)


# Except for a few alleles, HWE can be roughly assumed


rm(list=setdiff(ls(), c("counts")))


###
#### LOGISTIC REGRESSION-----
###


# Integer encoding for Population
counts <- read.csv("counts.csv") %>% 
  remove_rownames() %>% 
  column_to_rownames(var="X")
counts[,1] <- ifelse(counts$Population =="EAS", 0, 1)
y <- counts[,1]
table(counts[,1])


# Get the model matrix for glmnet, remove intercept column
X <- model.matrix(lm(Population ~ ., data = counts))[,-1]


# Are the predictors standardized?
plot(as.ts(apply(X,2,sd)), xlab="covariate", ylab="sd")


# Create test, train, validation sets
set.seed(12345) 
test.index <- sample.int(nrow(X), ceiling(nrow(X)*.2), replace = FALSE)
train.x <- X[-test.index,]
train.y <- y[-test.index]
test.x <- X[test.index,]
test.y <- y[test.index]


set.seed(12345)
validation.index <- sample.int(nrow(train.x), ceiling(nrow(train.x)*.5), replace = FALSE)
validation.x <- train.x[validation.index,]
validation.y <- train.y[validation.index]
train.x <- train.x[-validation.index,]
train.y <- train.y[-validation.index]


rm(test.index, validation.index)


# Train logistic Lasso
set.seed(12345) 
train.model <- cv.glmnet(train.x, train.y, family = "binomial", type.measure = "auc", nfolds = 5, keep = TRUE)
setfolds <- train.model$foldid


plot(train.model, )


# Observe key aspects of the glmnet object
rbind(c("lambda.min", "lambda.1se"),
      round(cbind(train.model$lambda.min, train.model$lambda.1se),5),
      train.model$nzero[train.model$index+1],
      train.model$cvm[train.model$index+1])

kept <- coef(train.model) %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("Feature") %>% 
  .[-which(.[,2] == 0),]
kept


# Validate the best lambda value with lambda.min and lambda.1se
predict.min <- predict(train.model, newx = validation.x, type = "response", s = train.model$lambda.min)[,1]
predict.1se <- predict(train.model, newx = validation.x, type = "response", s = train.model$lambda.1se)[,1]


# Calculate and plot prediction accuracy (ROC)
auc.min <- roc(validation.y, predict.min); auc.min$auc
auc.1se <- roc(validation.y, predict.1se); auc.1se$auc


par(mfrow=c(1,2))
plot(auc.min, main = "Logistic regression with lambda.min")
plot(auc.1se, main = "Logistic regression with lambda.1se")
par(mfrow=c(1,1))


# Create confusion matrices for the validation set using 0.5 threshold
validation.y[validation.y == 1] <- "TRUE"
validation.y[validation.y == 0] <- "FALSE"
confusionMatrix(data = as.factor(predict.min > 0.5), reference = as.factor(validation.y))
confusionMatrix(data = as.factor(predict.1se > 0.5), reference = as.factor(validation.y))


# They both perform extremely well, so go forward with the more parsimonious lambda.1se model on the test set
predict.test <- predict(train.model, newx = test.x, type = "response", s = train.model$lambda.1se)[,1]
auc.test <- roc(test.y, predict.test); auc.test$auc


# Create confusion matrix for test set using 0.5 threshold
test.y[test.y == 1] <- "TRUE"
test.y[test.y == 0] <- "FALSE"
confusionMatrix(data = as.factor(predict.test > 0.5), reference = as.factor(test.y))


# Find optimal thresholds on the train set
snsp.train <- cbind(auc.1se$sensitivities, auc.1se$specificities)
head(snsp.train)
indx <- which.max(apply(snsp.train, 1, min))
snsp.train[indx,]

auc.1se$thresholds[indx]


# Observe on the validation set
confusionMatrix(data = as.factor(predict.min > 0.47), reference = as.factor(validation.y))
confusionMatrix(data = as.factor(predict.1se > 0.47), reference = as.factor(validation.y))


# Find optimal thresholds for test set
snsp.test <- cbind(auc.test$sensitivities, auc.test$specificities)
indx2 <- which.max(apply(snsp.test, 1, min))  ### min-max approach!
snsp.test[indx2,]

auc.test$thresholds[indx2]


# Calculate and plot prediction accuracy (ROC) for test set with optimized 
confusionMatrix(data = as.factor(predict.test > 0.47), reference = as.factor(test.y))


plot(auc.test, main = "Logistic regression on test data")
abline(h = snsp.test[indx2, 1], v = snsp.test[indx2, 2], col = 'blue', lty = 2)


