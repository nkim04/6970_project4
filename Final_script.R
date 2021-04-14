###Overall Analysis ###

#### Load Packages ####

library(vcfR)
library(adegenet)
library(spdep)
library(tidyverse)
library(HardyWeinberg)
library(factoextra)
library(PCAtools)
library(pca3d)
library(cluster)
library(rpart)
library(rpart.plot)
library(rattle)
library(caret)
library(pROC)
library(glmnet)

##
#### LOAD AND FILTER DATA ####
##

# ALDH2
ALDH2_EAS_file <- "12.111766933-111817532_ALDH2_EAS.vcf.gz"
ALDH2_EUR_file <- "12.111766933-111817532_ALDH2_EUR.vcf.gz"


# OCA2
OCA2_EAS_file <- "15.27719008-28099342_OCA2_EAS.vcf.gz"
OCA2_EUR_file <- "15.27719008-28099342_OCA2_EUR.vcf.gz"


# CREB1
CREB1_EAS_file <- "2.207529943-207605988_CREB1_EAS.vcf.gz"
CREB1_EUR_file <- "2.207529943-207605988_CREB1_EUR.vcf.gz"


# SLC45A2
SLC45A2_EAS_file <- "5.33944623-33984693_SLC45A2_EAS.vcf.gz"
SLC45A2_EUR_file <- "5.33944623-33984693_SLC45A2_EUR.vcf.gz"

genes = c("ALDH2", "OCA2", "CREB1", "SLC45A2")


for (gene in genes) {
  
  
  # Import files
  EAS_vcf <- read.vcfR(get(paste0(gene, "_EAS_file")), verbose = FALSE)
  EUR_vcf <- read.vcfR(get(paste0(gene, "_EUR_file")), verbose = FALSE)
  
  
  # Get matrices of allele counts
  EAS_counts <- vcfR2genind(EAS_vcf)@tab
  EUR_counts <- vcfR2genind(EUR_vcf)@tab
  
  
  # There are different numbers of alleles in the populations
  # Let's see which ones overlap, and put them in one df
  SNP_counts <- bind_rows(as.data.frame(EAS_counts), as.data.frame(EUR_counts)) %>% 
    .[, colSums(is.na(.)) == 0]
  
  
  # Get allele frequencies and filter out the rare (<0.001) ones from our original
  # Also filter out frequencies equal to 1 (same across populations)
  SNP_freq <- apply(SNP_counts, 2, function(x) {sum(x)/(2*nrow(SNP_counts))})
  SNP_counts <- SNP_counts[, names(SNP_freq[SNP_freq > 0.001 & SNP_freq != 1])] %>% 
    cbind(Population = c(rep("EAS", nrow(EAS_counts)), rep("EUR", nrow(EUR_counts))), .) # Add population labels as the first column
  
  
  # For bi-allelic SNPs, remove their alternate alleles ("1") - this is redundant information
  # For the multi-allelic SNP, still remove the alternate "1" and keep all othe alleles - they are still informative
  SNP_counts <- SNP_counts[, -grep(".1$", colnames(SNP_counts))]
  
  #return the dataframe for each gene
  assign(paste0(gene,"_counts"), SNP_counts)
  
}

#combine each of the individual genes into one overall dataframe 
sum(rownames(ALDH2_counts) %in% rownames(CREB1_counts))
sum(rownames(CREB1_counts) %in% rownames(OCA2_counts))
sum(rownames(OCA2_counts) %in% rownames(SLC45A2_counts))


counts <- cbind(ALDH2_counts, 
                CREB1_counts[2:ncol(CREB1_counts)],
                OCA2_counts[2:ncol(OCA2_counts)],
                SLC45A2_counts[2:ncol(SLC45A2_counts)])

#remove everything except final dataframe
#cleans up environment
rm(list=setdiff(ls(), "counts"))

##
### HARDY WEINBERG EQUILIBRIUM TEST ####
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

##
#### LOGISTIC REGRESSION ####
##

# Read and scale data
counts_raw <- read.csv("counts.csv") %>% 
  remove_rownames() %>% 
  column_to_rownames(var="X")
counts <- cbind("Population" = counts_raw[,1], scale(counts_raw[, 2:ncol(counts_raw)])) %>% 
  as.data.frame()
counts[,2:ncol(counts)] <- sapply(counts[,2:ncol(counts)],as.numeric)

table(round(apply(counts[2:ncol(counts)], 2, mean), 2))
table(round(apply(counts[2:ncol(counts)], 2, sd), 2))


# Integer encoding for Population
counts[,1] <- ifelse(counts[,1] =="EAS", 0, 1)
y <- counts[,1]
table(y)


# Get the model matrix for glmnet, remove intercept column
X <- model.matrix(lm(Population ~ ., data = as.data.frame(counts)))[,-1]


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


plot(train.model)


# Observe key aspects of the glmnet object
rbind(c("lambda.min", "lambda.1se"),
      round(cbind(train.model$lambda.min, train.model$lambda.1se),5),
      train.model$nzero[train.model$index+1],
      train.model$cvm[train.model$index+1])

selected_features <- coef(train.model) %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("Feature") %>% 
  .[-which(.[,2] == 0),]
selected_features <- cbind(selected_features[-1,], abs(selected_features[-1,2]))
selected_features[order(selected_features[, 3], decreasing = T),]


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


# Influential features
coef(train.model, s = train.model$lambda.1se)
## In vector form:
coef.min <- coef(cvx, s = cvx$lambda.1se)[,1]
coef.min[coef.min != 0] ## The selected ones only
names(coef.min[coef.min != 0])[-1]

##
#### PCA ####
##

#Perform PCA on counts matrix
pca_allele_counts<-prcomp(allele_counts_df[,3:ncol(allele_counts_df)],center=T)

#Plot eigenvalues in a scree plot to show variances explained by each PC. 
fviz_eig(pca_allele_counts, main="Screen Plot of first 10 PCs",addlabels=T,ggtheme=theme_minimal())
findElbowPoint(pca_allele_counts$sdev^2) #35

#PCA biplot of PC1 & PC2 colored by superpopulation
fviz_pca_biplot(pca_allele_counts,
                select.var=list(contrib=10),
                label="var",
                repel=TRUE,
                addEllipses=F,
                col.ind=allele_counts_df$Population,
                #col.ind=allele_counts_df$Population,
                legend.title="Population",
                labelsize=3,
                max.overlaps=3,
                title="PCA Biplot for Allele Counts, Colored by Superpopulation",
                subtitle="Depicting Scores and Top 10 Loadings for the first two PCs",
                alpha.var=0.3,
                pointsize=2,
                alpha.ind=0.8)+
  scale_shape_manual(values=c(9,19,20,21,22,23,24,8,2,3))

#3D plot colored by superpopulation
pca3d(pca_allele_counts,biplot=TRUE,biplot.vars=1,show.ellipses=F,group=allele_counts_df$Superpopulation.code)

#Combine PCA scores with metadata (sample names and population) to plot other PCs
pca_allele_counts_df=data.frame(pca_allele_counts$x)
pca_meta=cbind(X=allele_counts_df$X,Population=allele_counts_df$Population,pca_allele_counts_df)

# Might other PCs have more discrimatory power for superpopulations?
#Scatterplot of 1st and 3rd PCs colored by superpopulation
ggplot(pca_meta,aes(x=PC1,y=PC3,color=Population))+
  geom_point()


## Hierarchical Clustering on Original Data----------

#Convert counts+metadata df to matrix to permit duplicate rownames
allele_counts_matrix<-as.matrix(allele_counts_df)

#Find column index for superpopulation.code
#grep("Superpopulation.code", colnames(allele_counts_matrix)) #1549

#Use Population as row names
rownames(allele_counts_matrix)<-allele_counts_matrix[,2]

#Remove columns with metadata
#colnames(allele_counts_matrix)[1545:1552]
#allele_counts_matrix<-allele_counts_matrix[,-c(1,seq(1545,1552))]
#allele_counts_matrix<-allele_counts_matrix[,-1]

#note, can't scale- will introduce NAs in all 0 cols
#scale_allele_counts<-scale(allele_counts) 

#Create distance matrix for rows (samples); (exclude cols 'X' (sample name') and 'Population)
dist_samples<-get_dist(allele_counts_matrix[,3:ncol(allele_counts_matrix)],method="binary",stand=FALSE)

#Perform agglomerative clustering (using complete linkage and default euclidean distance)
agnes<-agnes(abs(dist_samples), method = "complete")

#Create object of class 'dendrogram' to facilitate manipulation
agnes_plot<- as.dendrogram(agnes)

#Assign colors to the labels of the dendrogram
library(dendextend)
colors_to_use <- ifelse(allele_counts_df$Population=="EUR", 1, 2)
colors_to_use <- colors_to_use[order.dendrogram(agnes_plot)]
labels_colors(agnes_plot) <- colors_to_use

#Plot whole dendrogram (colored by leaf label)
plot(agnes_plot, main="Dendrogram of Samples Colored by Population (on original data)")

#Plot dendrogram with split branches (colored by Population); easier to visualize
par(mfrow=c(5,1))
plot(cut(agnes_plot, h=0.11)$upper, 
     main="Upper tree of cut at h=0.11")
plot(cut(agnes_plot, h=0.11)$lower[[1]], 
     main="First branch of lower tree with cut at h=45")
plot(cut(agnes_plot, h=0.11)$lower[[2]], 
     main="Second branch of lower tree with cut at h=45")
plot(cut(agnes_plot, h=0.11)$lower[[3]], 
     main="Third branch of lower tree with cut at h=45")
plot(cut(agnes_plot, h=0.11)$lower[[4]], 
     main="Fourth branch of lower tree with cut at h=45")


## Hierarchical Clustering on top three PCs------------

#Convert counts+metadata df to matrix to permit duplicate rownames
pca_meta_matrix<-as.matrix(pca_meta)

#Use superpopulation.code as row names
rownames(pca_meta_matrix)<-pca_meta_matrix[,2]

#Remove column with superpopulation.code
#pca_meta_matrix<-pca_meta_matrix[,-1]

#Subset just first 3 PCs & metadata
pca_meta_matrix_3<-pca_meta_matrix[,1:5]

#Create distance matrix for rows (samples)
dist_samples_pca<-get_dist(pca_meta_matrix_3[,3:5],method="euclidean",stand=FALSE)

#Plot dendrogram (super slow)
#fviz_dist(dist_samples)

#Perform agglomerative clustering (using complete linkage and default euclidean distance)
agnes_pca<-agnes(abs(dist_samples_pca), method = "complete")

#Create object of class 'dendrogram' to facilitate manipulation
agnes_pca_plot<- as.dendrogram(agnes_pca)

#Assign colors to the labels of the dendrogram
library(dendextend)
colors_to_use <- ifelse(allele_counts_df$Population=="EUR", 1, 2)
colors_to_use <- colors_to_use[order.dendrogram(agnes_pca_plot)]
labels_colors(agnes_pca_plot) <- colors_to_use

#Plot whole dendrogram (colored by leaf label)
plot(agnes_pca_plot, main="Dendrogram of Samples Colored by Population (on top 3 PCs)")

#Plot dendrogram with split branches (colored by Population); easier to visualize
par(mfrow=c(4,1))
plot(cut(agnes_pca_plot, h=45)$upper, 
     main="Upper tree of cut at h=45")
plot(cut(agnes_pca_plot, h=45)$lower[[1]], 
     main="First branch of lower tree with cut at h=45")
plot(cut(agnes_pca_plot, h=45)$lower[[2]], 
     main="Second branch of lower tree with cut at h=45")
plot(cut(agnes_pca_plot, h=45)$lower[[3]], 
     main="Third branch of lower tree with cut at h=45")



#Hierarchical clustering comparison ------

#Comparison of dendrograms of original data vs top 3 PCs
par(mfrow=c(2,1))
plot(agnes_plot, main="Dendrogram of Samples Colored by Population (on original data)")
plot(agnes_pca_plot, main="Dendrogram of Samples Colored by Population (on top 3 PCs)")


#Non-hierarchical clustering on Original Data-------

#Note: use dataframe not matrix, because 1/0 columns needs to be integer/numberic not character

#Find optimal number of clusters 
fviz_nbclust(allele_counts_df[,3:ncol(allele_counts_df)], kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method")

#k-means for k=2:3
kmean2 <- kmeans(allele_counts_matrix[,3:ncol(allele_counts_df)],2,nstart=25)
kmean3 <- kmeans(allele_counts_matrix[,3:ncol(allele_counts_df)],3,nstart=25)

#Plot 2-means
fviz_cluster(kmean2,
             data=allele_counts_df[,3:ncol(allele_counts_df)],
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             stand=FALSE,
             geom = c("text","point"),
             labelsize=6,
             ellipse.type = "convex", 
             ellipse.alpha=0.1,
             show.clust.cent = TRUE,
             outlier.color = "black",
             main="PCA Plot of Samples Clustered by 2-medoids",
             ggtheme = theme_bw()
)

clusplot(dat, kmean2$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)

?fviz_cluster
#Create DF of cluster assignments and real superpop.code
kmean_df<-data.frame(Cluster.3kmean=kmean3$cluster,Population=allele_counts_df$Population)
kmean_df$Cluster.2kmean<-kmean2$cluster

#Histogram of samples per ea. cluster for k=2, colored by superpopulation.code
ggplot(data = kmean_df, aes(y = Cluster.2kmean)) +
  geom_bar(aes(fill = Population)) +
  ggtitle("Count of Samples per Cluster by Type (2-med") +
  labs(y="Cluster")+
  theme(plot.title = element_text(hjust = 0.5))


#Non-hierarchical clustering on PCA-------

#Making numeric values finite
pca_meta[,3:ncol(pca_meta)]<-signif(pca_meta[,3:ncol(pca_meta)],15)

#Subset just first 3 PCs; performed worse
#pca_meta_3<-pca_meta[,1:4]

#Find and visualize optimal number of clusters using elbow method
fviz_nbclust(pca_meta[,3:ncol(pca_meta)], kmeans, method = "wss",print.summary=TRUE) + geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method (clustering samples)")

#k-means for k=2:3
kmean2_pca <- kmeans(pca_meta[,3:ncol(pca_meta)],2,nstart=25)
kmean3_pca <- kmeans(pca_meta[,3:ncol(pca_meta)],3,nstart=25)

#Plot 2-means
fviz_cluster(kmean2_pca,
             data=pca_meta[,3:ncol(pca_meta)],
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = c("text","point"),
             labelsize=6,
             ellipse.type = "convex", 
             ellipse.alpha=0.1,
             show.clust.cent = TRUE,
             main="PCA Plot of Samples Clustered by 2-medoids",
             ggtheme = theme_bw()
)


#Plot 3-means
fviz_cluster(kmean3_pca,
             data=pca_meta[,3:ncol(pca_meta)],
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = c("text","point"),
             labelsize=6,
             ellipse.type = "convex", 
             ellipse.alpha=0.1,
             outlier.color = "black",
             show.clust.cent=T,
             main="PCA Plot of Samples Clustered by 3-medoids",
             ggtheme = theme_bw()
)


#Combine clustering assignments with k-means cluster assignments
kmean_pca_df<-data.frame(Cluster.2kmean=as.factor(kmean2_pca$cluster),Cluster.3kmean=as.factor(kmean3_pca$cluster),pca_meta)

#Histogram of samples per ea. cluster for k=2, colored by superpopulation.code
ggplot(data = kmean_pca_df, aes(y = Cluster.2kmean)) +
  geom_bar(aes(fill = Population)) +
  ggtitle("Count of Samples per Cluster by Type (2-med") +
  labs(y="Cluster")+
  theme(plot.title = element_text(hjust = 0.5))
)

#Scatterplot of PC1 & PC2 colored by population and shaped by cluster assignment in 2-kmeans
ggplot(kmean_pca_df,aes(x=PC1,y=PC2,color=Population,shape=Cluster.2kmean))+
  geom_point()


# Logistic regression (on all PCs) with alpha=1------------

library(glmnet)

#Encode response variable (Population) as numeric: 0 (EAS) or 1 (EUR)
y <- ifelse(pca_meta$Population =="EUR", 0, 1)
table(y)
pca_meta_coded<-pca_meta
pca_meta_coded$Population<-y

#Fit the linear model  (exclude the metadata column that's neither response nor predictor)
f0 <- lm(Population ~ ., data=pca_meta_coded[,-1])

#Create a model matrix and Remove intercept
X <- model.matrix(f0)[,-1]

set.seed(1545) 

#Create test indices (20%)
test.index <- sample.int(dim(X)[1],round(dim(X)[1]*.2), replace = FALSE)

#Train logistic Lasso on test set (performs 10-CV to tune)
cvx <- cv.glmnet(X[-test.index,],y[-test.index], nfolds = 5, family="binomial", alpha=1, type.measure = "auc")
plot(cvx)

## Predict values for test and train set (using observed values)
prds.train <- predict(cvx,newx = X[-test.index,], type = "response", s=cvx$lambda.min)[,1]
prds.test <- predict(cvx,newx = X[test.index,], type = "response", s=cvx$lambda.min)[,1]

# Calculate and plot prediction accuracy for train set and test set
library(pROC)
auc.train <- roc(y[-test.index],prds.train) 

par(mfrow=c(1,2))
plot(auc.train)

auc.test <- roc(y[test.index],prds.test)
auc.test
plot(auc.test)
par(mfrow=c(1,1))

#Create Confusion matrix for training set if threshold is set to 0.5% probability
#3 incorrect classifications 
conf.mat1<- table(y=y[-test.index],yhat=as.numeric(prds.train>.5))

## A function to compute Sensitivity and Specificity
sn.sp <- function(mat){
  sn <- mat[2,2]/sum(mat[2,])
  sp <- mat[1,1]/sum(mat[1,])
  return(unlist(list(sensitivity=sn, specificity=sp)))
}

#Find optimal thresholds for training set
snsp.train <- cbind(auc.train$sensitivities,auc.train$specificities)
head(snsp.train)
indx <- which.max(apply(snsp.train,1,min))  ### min-max approach!
snsp.train[indx,]

cutoff <- auc.train$thresholds[indx]
cutoff

#Find optimal thresholds for test set
snsp.test <- cbind(auc.test$sensitivities,auc.test$specificities)
indx2 <- which.max(apply(snsp.test,1,min))  ### min-max approach!
snsp.test[indx2,]

cutoff2 <- auc.test$thresholds[indx2]
cutoff2

#Plot ROC curves for training and test @ optimal thresholds
par(mfrow=c(1,2))
plot(auc.train)
abline(h=snsp.train[indx,1],v=snsp.train[indx,2], col='blue', lty=2)
plot(auc.test)
abline(h=snsp.test[indx2,1],v=snsp.test[indx2,2], col='blue', lty=2)
par(mfrow=c(1,1))

#Sensitivity and specificity at threshold for train set
sn.sp(table(y=y[-test.index],yhat=as.numeric(prds.train>cutoff)))

#Sensitivity and specificity at threshold for test set
sn.sp(table(y=y[ test.index],yhat=as.numeric(prds.test >cutoff2)))



#Classification Tree on PCA data-----
library(rpart)
library(rattle)

fancyRpartPlot(classtree_pca)

par(mfrow=c(1,2))
plot(classtree_pca,margin=0.1)
text(classtree_pca,use.n=T,cex=1.3)

### Decision Trees and Bagging
#Previously very rare SNPs were filtered out but very common SNPs (>0.999) are also not of much use during classification
#SNPs that have a high frequency means a lot of the samples, from both super populations, have that genotype
#Differences in SNPs are what allow us to classify, therefore common SNPs will be removed 
counts <- read.csv("counts.csv")
rownames(counts) <- counts[,1]
counts <- counts[,-1]
SNP_freq <- apply(counts[,-1], 2, function(x) {sum(x)/(2*nrow(counts))})
counts <- counts[, c("Population",names(SNP_freq[SNP_freq < 0.999]))]

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
