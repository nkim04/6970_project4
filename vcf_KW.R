library(vcfR)
library(adegenet)
library(spdep)
library(factoextra)
library(pca3d)
library(PCAtools)
library(cluster)
library(caret)
library(dendextend)
library(tidyverse)

#### For ALDH2 (alcohol) in Asians and Europeans

#This VCF into vcf object
#503 samples, 1316 variants, no NAs
vcf <- read.vcfR( "12.111766933-111817532.ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf", verbose = FALSE )

#Get the vcf data in tabular format (easier to navigate)
vcf_df <- cbind(getFIX(vcf), vcf@gt)

#To get a matrix of allele counts, turn into a genind object
my_genind <- vcfR2genind(vcf)

#Get a matrix of allele counts
all_allele_counts <- my_genind@tab

#Get allele frequencies, and count those with lower than 0.001 frequency
all_allele_freqs <- apply(all_allele_counts, 2, function(x) {sum(x)/(2*nrow(all_allele_counts))})
all_allele_freqs[all_allele_freqs < 0.001] #to get rare alleles and their frequencies

allele_freqs <- all_allele_freqs[all_allele_freqs > 0.001]
length(all_allele_freqs); length(allele_freqs)

#removing rare alleles
allele_counts <- all_allele_counts[, all_allele_freqs > 0.001]
dim(allele_counts)

#Read in metadata for asian & european superpops; combine rows
meta_EUR<-read.table("igsr_samples.tsv",header=T,sep="\t")
meta_EAS<-read.table("igsr_samples_EAS.tsv",header=T,sep="\t")
meta_EUR_EAS<-bind_rows(meta_EUR,meta_EAS)

#Convert matrix of allele counts to dataframe, place sample names as 1st column, and combine with metadata
allele_counts_df<-as.data.frame(allele_counts)
allele_counts_df<-allele_counts_df %>% 
  rownames_to_column(var="Sample.name") %>% 
  left_join(y=meta_EUR_EAS,by="Sample.name")

#View sample sizes for each superpop (504 EAS, 503 EUR)
allele_counts_df %>% 
  group_by(Superpopulation.code) %>% 
  summarize(n())

#Alternatively, using counts.csv
allele_counts_df <- read.csv("counts.csv")

## PCA ------

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

#Scatterplot of PC1 & PC2 colored by population and shaped by cluster assignment in k=2
ggplot(kmean_pca_df,aes(x=PC1,y=PC2,color=Population,shape=Cluster.2kmean))+
  geom_point()

#Function to compute sensitivity and specificity
sn.sp <- function(mat){
  sn <- mat[2,2]/sum(mat[2,])
  sp <- mat[1,1]/sum(mat[1,])
  return(unlist(list(sensitivity=sn, specificity=sp)))
}

#Confusion Matrix, sensitivity, and specificity for k=2
table(kmean_pca_df$Cluster.2kmean,kmean_pca_df$Population)
sn.sp(table(kmean_pca_df$Cluster.2kmean,kmean_pca_df$Population))


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


#Natalie's Folds----------------

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

#-------------------


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



#K-means using folds from logistic regression----

#Create list of 10 folds
folds_for_kmean2<-createFolds(pca_meta$Population,k=10)

#Update folds based on 'setfolds'
for(i in 1:10){folds_for_kmean2[[i]]<-which(train.model$foldid==i)}


#Create function to do k-means (k=2) on perturbed datasets
kmeans_folds_fn<-function(x){
  fold<-folds_for_kmean2[[x]]
  pca_meta_subset<-pca_meta[,3:ncol(pca_meta)]
  train<-pca_meta_subset[-fold,]
  #print(length((pca_meta_subset[-fold])))
  km<-kmeans(x=pca_meta_subset[-fold,],centers=2,nstart=25)
  print(table(as.factor(km$cluster),pca_meta[-fold,]$Population))
  sn.sp(table(as.factor(km$cluster),pca_meta[-fold,]$Population))
  }

#Run k=means for each -fold
kmeans_folds<-c()
for (i in 1:10) {kmeans_folds[[i]]<-kmeans_folds_fn(i)}

#ERROR:
#Error in kmeans(pca_meta_subset[-fold, ], 2, nstart = 25) : more cluster centers than distinct data points. 