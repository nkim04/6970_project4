library(vcfR)
library(adegenet)
library(spdep)
library(factoextra)
library(pca3d)
library(PCAtools)
library(cluster)
library(caret)
library(dendextend)
library(glmnet)
library(tidyverse)

#Read in data
allele_counts_df <- read.csv("counts.csv")


## PCA ------

#Perform PCA on counts matrix
pca_allele_counts<-prcomp(allele_counts_df[,3:ncol(allele_counts_df)],center=T)

#Plot eigenvalues in a scree plot to show variances explained by each PC. 
fviz_eig(pca_allele_counts, main="Screen Plot of first 10 PCs",addlabels=T,ggtheme=theme_minimal())
findElbowPoint(pca_allele_counts$sdev^2) #35

#PCA biplot of PC1 & PC2 colored by Population
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
                title="PCA Biplot for Allele Counts, Colored by Population",
                subtitle="Depicting Scores and Top 10 Loadings for the first two PCs",
                alpha.var=0.3,
                pointsize=2,
                alpha.ind=0.8)+
scale_shape_manual(values=c(9,19,20,21,22,23,24,8,2,3))

#3D plot colored by Population
pca3d(pca_allele_counts,biplot=TRUE,biplot.vars=1,show.ellipses=F,group=allele_counts_df$Population)

#Combine PCA scores with metadata (sample names and population) to plot other PCs
pca_allele_counts_df=data.frame(pca_allele_counts$x)
pca_meta=cbind(X=allele_counts_df$X,Population=allele_counts_df$Population,pca_allele_counts_df)

# Might other PCs have more discrimatory power for Populations?
#Scatterplot of 1st and 3rd PCs colored by Population
ggplot(pca_meta,aes(x=PC1,y=PC3,color=Population))+
  geom_point()



## Hierarchical Clustering on Original Data----------


#Convert counts+metadata df to matrix to permit duplicate rownames
allele_counts_matrix<-as.matrix(allele_counts_df)

#Use Population as row names
rownames(allele_counts_matrix)<-allele_counts_matrix[,2]

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

agnes_plot<-color_branches(agnes_plot, k=4,col = c("deepskyblue","orange","darkolivegreen3","purple"))

#Plot whole dendrogram (colored by leaf label)
plot(agnes_plot, main="Dendrogram of Samples Colored by Population (on original data)")

#Plot dendrogram with split branches (colored by Population); easier to visualize
par(mfrow=c(5,1))
plot(cut(agnes_plot, h=0.061)$upper, 
     main="Dendrogram of Samples Colored by Population (on original data); Upper tree of cut at h=0.061")
plot(cut(agnes_plot, h=0.061)$lower[[1]], 
     main="First branch of lower tree with cut at h=0.061")
plot(cut(agnes_plot, h=0.061)$lower[[2]], 
     main="Second branch of lower tree with cut at h=0.061")
plot(cut(agnes_plot, h=0.061)$lower[[3]], 
     main="Third branch of lower tree with cut at h=0.061")
plot(cut(agnes_plot, h=0.061)$lower[[4]], 
     main="Fourth branch of lower tree with cut at h=0.061")

#View Counts (if split into 2 clusters)
clusterCut <- dendextend::cutree(agnes, 2)
clusterCut <- dendextend::cutree(agnes,k=4)
table(Dend=clusterCut, True=allele_counts_df$Population)




## Hierarchical Clustering on top three PCs------------

#Convert counts+metadata df to matrix to permit duplicate rownames
pca_meta_matrix<-as.matrix(pca_meta)

#Use Population as row names
rownames(pca_meta_matrix)<-pca_meta_matrix[,2]

#Create distance matrix for rows (samples)
dist_samples_pca<-get_dist(pca_meta_matrix[,3:5],method="euclidean",stand=FALSE)

#Perform agglomerative clustering (using complete linkage and default euclidean distance)
agnes_pca<-agnes(abs(dist_samples_pca), method = "average")

#Create object of class 'dendrogram' to facilitate manipulation
agnes_pca_plot<- as.dendrogram(agnes_pca)

#Assign colors to the labels of the dendrogram
colors_to_use <- ifelse(allele_counts_df$Population=="EUR", 1, 2)
colors_to_use <- colors_to_use[order.dendrogram(agnes_pca_plot)]
labels_colors(agnes_pca_plot) <- colors_to_use

#Color branches
agnes_pca_plot<-color_branches(agnes_pca_plot, k=2,col=c("deepskyblue","darkolivegreen3"))

#Plot whole dendrogram (colored by leaf label)
plot(agnes_pca_plot, main="Dendrogram of Samples Colored by Population (on top 3 PCs)")

#Plot dendrogram with split branches (colored by Population); easier to visualize
par(mfrow=c(4,1))
plot(cut(agnes_pca_plot, h=45)$upper, 
     main="Dendrogram of Samples Colored by Population (on top 3 PCs); Upper tree of cut at h=45")
plot(cut(agnes_pca_plot, h=45)$lower[[1]], 
     main="First branch of lower tree with cut at h=45")
plot(cut(agnes_pca_plot, h=45)$lower[[2]], 
     main="Second branch of lower tree with cut at h=45")
plot(cut(agnes_pca_plot, h=45)$lower[[3]], 
     main="Third branch of lower tree with cut at h=45")

#View counts (if split into 3 clusters)
clusterCut_pca <- dendextend::cutree(agnes_pca,k=2)
table(Dend=clusterCut_pca, True=allele_counts_df$Population)


#Hierarchical clustering comparison ------

#Comparison of dendrograms of original data vs top 3 PCs
par(mfrow=c(2,1))
plot(agnes_plot, main="Dendrogram of Samples Colored by Population (on original data)")
plot(agnes_pca_plot, main="Dendrogram of Samples Colored by Population (on top 3 PCs)")





# Logistic regression Setup (From Natalie's script)------

#Encode response variable (Population) as numeric: 0 (EAS) or 1 (EUR)
y <- ifelse(pca_meta$Population =="EUR", 0, 1)
table(y)
pca_meta_coded<-pca_meta
pca_meta_coded$Population<-y

#Fit the linear model  (exclude the metadata column that's neither response nor predictor)
f0 <- lm(Population ~ ., data=pca_meta_coded[,-1])

#Create a model matrix and Remove intercept
X <- model.matrix(f0)[,-1]

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



#Non-hierarchical (k-means) clustering on top 3 PCs-------

#Making numeric values finite
pca_meta[,3:ncol(pca_meta)]<-signif(pca_meta[,3:ncol(pca_meta)],15)

#Find and visualize optimal number of clusters using elbow method
fviz_nbclust(pca_meta[,3:5], kmeans, method = "wss",print.summary=TRUE) + geom_vline(xintercept = 2, linetype = 2)+labs(subtitle = "Elbow method (clustering samples)")

fviz_nbclust(pca_meta[,3:5], kmeans, method = "gap_stat",print.summary=TRUE) + 
  labs(subtitle = "Gap statistic method (clustering samples)")

#k-means for k=2:3
kmean2_pca <- kmeans(pca_meta[,3:5],2,nstart=25)
kmean3_pca <- kmeans(pca_meta[,3:5],3,nstart=25)
kmean4_pca <- kmeans(pca_meta[,3:5],4,nstart=25)

#Combine clustering assignments with k-means cluster assignments
kmean_pca_df<-data.frame(Cluster.2kmean=as.factor(kmean2_pca$cluster),Cluster.3kmean=as.factor(kmean3_pca$cluster),Cluster.4kmean=as.factor(kmean4_pca$cluster),pca_meta)

#Histogram of samples per ea. cluster for k=2, colored by Population
ggplot(data = kmean_pca_df, aes(y = Cluster.2kmean)) +
  geom_bar(aes(fill = Population)) +
  ggtitle("Count of Samples per Cluster by Population (2-means)") +
  labs(y="Cluster")+
  theme(plot.title = element_text(hjust = 0.5))

#Histogram of samples per ea. cluster for k=3, colored by Population
ggplot(data = kmean_pca_df, aes(y = Cluster.3kmean)) +
  geom_bar(aes(fill = Population)) +
  ggtitle("Count of Samples per Cluster by Population (3-means)") +
  labs(y="Cluster")+
  theme(plot.title = element_text(hjust = 0.5))

#Histogram of samples per ea. cluster for k=4, colored by Population
ggplot(data = kmean_pca_df, aes(y = Cluster.4kmean)) +
  geom_bar(aes(fill = Population)) +
  ggtitle("Count of Samples per Cluster by Population (4-means)") +
  labs(y="Cluster")+
  theme(plot.title = element_text(hjust = 0.5))

#Scatterplot of PC1 & PC2 colored by population and shaped by cluster assignment in k=2:4
library(ggpubr)
ggscatter(kmean_pca_df,x="PC1",y="PC2",color="Cluster.2kmean",shape="Population",ellipse=TRUE,ellipse.type="convex",legend="right",repel=TRUE,size=1.2, xlab="Dim 1 (17.1%)", ylab="Dim 2 (10.5%)")+ stat_mean(aes(color = Cluster.2kmean), size = 4)

ggscatter(kmean_pca_df,x="PC1",y="PC2",color="Cluster.3kmean",shape="Population",ellipse=TRUE,ellipse.type="convex",legend="right",repel=TRUE,size=1.2, xlab="Dim 1 (17.1%)", ylab="Dim 2 (10.5%)")+ stat_mean(aes(color = Cluster.2kmean), size = 4)

ggscatter(kmean_pca_df,x="PC1",y="PC2",color="Cluster.4kmean",shape="Population",ellipse=TRUE,ellipse.type="convex",legend="right",repel=TRUE,size=1.2, xlab="Dim 1 (17.1%)", ylab="Dim 2 (10.5%)")+ stat_mean(aes(color = Cluster.2kmean), size = 4)

#Function to compute sensitivity and specificity
sn.sp <- function(mat){
  sn <- mat[2,2]/sum(mat[2,])
  sp <- mat[1,1]/sum(mat[1,])
  return(unlist(list(sensitivity=sn, specificity=sp)))
}

#Confusion Matrix, sensitivity, and specificity for k=2
table(kmean_pca_df$Population,kmean_pca_df$Cluster.2kmean)
sn.sp(table(kmean_pca_df$Population,kmean_pca_df$Cluster.2kmean))

#Confusion Matrix for k=3
table(kmean_pca_df$Population,kmean_pca_df$Cluster.3kmean)

#Confusion Matrix for k=4
table(kmean_pca_df$Population,kmean_pca_df$Cluster.4kmean)



#K-means on top 3 PCs using folds from logistic regression----

#Create list of 5 folds
folds_for_kmean2<-createFolds(pca_meta$Population,k=5)

#Update folds based on 'setfolds'
for(i in 1:5){folds_for_kmean2[[i]]<-which(train.model$foldid==i)}

#Create function to do k-means (k=2) on perturbed datasets
kmeans_folds_fn<-function(x){
  fold<-folds_for_kmean2[[x]]
  pca_meta_subset<-pca_meta[,3:ncol(pca_meta)]
  train<-pca_meta_subset[-fold,]
  #print(length((pca_meta_subset[-fold])))
  km<-kmeans(x=pca_meta_subset[-fold,],centers=2,nstart=25)
  print(table(pca_meta[-fold,]$Population,as.factor(km$cluster)))
  }

#Run k=means for each -fold
kmeans_folds<-c()
for (i in 1:5) {kmeans_folds[[i]]<-kmeans_folds_fn(i)}

