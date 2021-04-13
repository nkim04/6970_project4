# Preliminary analyses


library(tidyverse)
library(HardyWeinberg)
library(factoextra)
library(PCAtools)
library(pca3d)
library(cluster)



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
#### Logistic Regression-----
###


library(glmnet)


# Integer encoding for Population
counts[,1] <- ifelse(counts$Population =="EAS", 0, 1)
y <- counts[,1]
table(counts[,1])


# Get the model matrix for glmnet, remove intercept column
f0 <- lm(Population ~ ., data = counts)
X <- model.matrix(f0)[,-1]
class(X)


# Create test, train, validation sets
set.seed(12345) 
test.index <- sample.int(nrow(X), ceiling(nrow(X)*.2), replace = FALSE)
train.x <- X[-test.index,]
train.y <- y[-test.index]
test.x <- X[test.index,]
test.y <- y[test.index]


validation.index <- sample.int(nrow(train.x), ceiling(nrow(train.x)*.5), replace = FALSE)
validation.x <- train.x[validation.index,]
validation.y <- train.y[validation.index]
train.x <- train.x[-validation.index,]
train.y <- train.y[-validation.index]


rm(train.index, validation.index)


#Train logistic Lasso
train.model <- cv.glmnet(train.x, train.y, family = "binomial", type.measure = "auc", nfolds = 5, keep = TRUE)
setfolds <- train.model$foldid


