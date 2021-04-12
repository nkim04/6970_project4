library(vcfR)
library(adegenet)
library(spdep)
library(factoextra)
library(pca3d)
library(PCAtools)
library(cluster)
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



## PCA ------

#Perform PCA on counts matrix
pca_allele_counts<-prcomp(allele_counts,center=T)

#Plot eigenvalues in a scree plot to show variances explained by each PC. 
fviz_eig(pca_allele_counts, main="Screen Plot of first 10 PCs",addlabels=T,ggtheme=theme_minimal())
findElbowPoint(pca_allele_counts$sdev^2) #10

#PCA biplot of PC1 & PC2 colored by superpopulation
fviz_pca_biplot(pca_allele_counts,
                select.var=list(contrib=10),
                label="var",
                repel=TRUE,
                addEllipses=F,
                col.ind=allele_counts_df$Superpopulation.code,
                #col.ind=allele_counts_df$Population.code,
                legend.title="SuperPopulation",
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

#Combine PCA scores and labels (superpopulation code) to plot other PCs
pca_allele_counts_df=data.frame(pca_allele_counts$x)
pca_meta=cbind(Superpopulation.code=allele_counts_df$Superpopulation.code,pca_allele_counts_df)

# Might other PCs have more discrimatory power for superpopulations?
#Scatterplot of 2rd and 3th PCs colored by superpopulation
ggplot(pca_meta,aes(x=PC2,y=PC3,color=Superpopulation.code))+
  geom_point()



## Hierarchical Clustering on Original Data----------

#Convert counts+metadata df to matrix to permit duplicate rownames
pca_meta_matrix<-as.matrix(pca_meta)

#Find column index for superpopulation.code
grep("Superpopulation.code", colnames(allele_counts_matrix)) #1549

#Use superpopulation.code as row names
rownames(allele_counts_matrix)<-allele_counts_matrix[,1549]

#Remove columns with metadata
colnames(allele_counts_matrix)[1545:1552]
allele_counts_matrix<-allele_counts_matrix[,-c(1,seq(1545,1552))]

#note, can't scale- will introduce NAs in all 0 cols
#scale_allele_counts<-scale(allele_counts) 

#Create distance matrix for rows (samples)
dist_samples<-get_dist(allele_counts_matrix,method="binary",stand=FALSE)

#Visualize distance matrix (if you dare..... super slow)
#fviz_dist(dist_samples)

#Agglomerative clustering
plot(agnes(abs(dist_samples), method = "complete"), main="using agnes()",which.plots=2)

#Find optimal number of clusters 
fviz_nbclust(allele_counts_matrix, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method")

## Hierarchical Clustering on 1st 3 PCs------------

#Convert counts+metadata df to matrix to permit duplicate rownames
pca_meta_matrix<-as.matrix(pca_meta)

#Use superpopulation.code as row names
rownames(pca_meta_matrix)<-pca_meta_matrix[,1]

#Remove column with superpopulation.code
pca_meta_matrix<-pca_meta_matrix[,-1]

#Subset just first 3 PCs
pca_meta_matrix_3<-pca_meta_matrix[,1:3]

#Create distance matrix for rows (samples)
dist_samples_pca<-get_dist(pca_meta_matrix_3,method="euclidean",stand=FALSE)

#Plot dendrogram (super slow)
#fviz_dist(dist_samples)



#Non-hierarchical clustering on PCA-------

#Making numeric values finite
pca_meta[,-1]<-signif(pca_meta[,-1],10)

#Subset just first 3 PCs; performed worse
#pca_meta_3<-pca_meta[,1:4]

#Find and visualize optimal number of clusters using elbow method
fviz_nbclust(pca_meta[,-1], kmeans, method = "wss",print.summary=TRUE) + geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method (clustering samples)")

#k-means for k=2:3
kmean2 <- kmeans(pca_meta[,-1],2,nstart=25)
kmean3 <- kmeans(pca_meta[,-1],3,nstart=25)

#Plot 2-means
fviz_cluster(kmean2,
             data=pca_meta[,-1],
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
fviz_cluster(kmean3,
             data=pca_meta[,-1],
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

#Create DF of cluster assignments and real superpop.code
kmean.df<-data.frame(Cluster.3kmean=kmean3$cluster,Superpopulation.code=pca_meta$Superpopulation.code)
kmean.df$Cluster.2kmean<-kmean2$cluster

#Histogram of samples per ea. cluster for k=2, colored by superpopulation.code
ggplot(data = kmean.df, aes(y = Cluster.2kmean)) +
  geom_bar(aes(fill = Superpopulation.code)) +
  ggtitle("Count of Samples per Cluster by Type (2-med") +
  labs(y="Cluster")+
  theme(plot.title = element_text(hjust = 0.5))



# Logistic regression------------

library(glmnet)

#encode superpopulation.code as 0(EAS) 1 (EUR)
y <- ifelse(pca_meta$Superpopulation.code =="EUR", 0, 1)
table(y)
pca_meta_coded<-pca_meta
pca_meta_coded[,1]<-y


#model matrix
f0 <- lm(Superpopulation.code ~ ., data=pca_meta_coded)

#remove intercept
X <- model.matrix(f0)[,-1]
class(X)

set.seed(1545) 

#Create test indices (20%)
test.index <- sample.int(dim(X)[1],round(dim(X)[1]*.2), replace = FALSE)

#Train logistic Lasso
cvx <- cv.glmnet(X[-test.index,],y[-test.index], nfolds = 5, family="binomial", alpha=1, type.measure = "auc")
plot(cvx)

## Predict values for test and train sets
prds.train <- predict(cvx,newx = X[-test.index,], type = "response", s=cvx$lambda.min)[,1]
prds.test <- predict(cvx,newx = X[test.index,], type = "response", s=cvx$lambda.min)[,1]

# Calculate and plot prediction accuracy
library(pROC)
auc.train <- roc(y[-test.index],prds.train) 
auc.train

par(mfrow=c(1,2))
plot(auc.train)
auc.test <- roc(y[test.index],prds.test)
auc.test
plot(auc.test)
par(mfrow=c(1,1))

#Create Confusion matrix for training set
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

#Plot ROC curves for training and test
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


