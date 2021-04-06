library(vcfR)
library(ggfortify)
library(plotly)
library(cluster)
library(rpart)
library(rpart.plot)
library(rattle)
#library(adegenet)
#library(spdep)

data <- read.delim("igsr_samples.tsv")
vcf <- read.vcfR( "7.22725889-22732002.ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", verbose = FALSE )
class(vcf)
head(vcf)

head(getFIX(vcf),20)
vcf@gt[1:10, 1:10]

info <- cbind(getFIX(vcf),vcf@gt)
geneID <- vcfR2genind(vcf)
alleles <- geneID@tab
#remove columns for allele 1 since they are just opposites of their respective allele 0 columns
alleles <- alleles[, -grep(".1$", colnames(alleles))]
allele_freqs <- apply(alleles, 2, function(x) {sum(x)/(2*nrow(alleles))})

#removed SNPs where frequency = 1 because that means all individuals have the reference alleles
for(i in 1:length(allele_freqs)){
  names <- NULL
  if(allele_freqs[i] == 0){
    names <- c(names,names(allele_freqs[i]))
  }
  alleles <- alleles[ , !(colnames(alleles) %in% names)]
}

#PCA
New_table <- data[data$Sample.name%in% rownames(alleles),]
allele_pca <- prcomp(alleles)
autoplot(allele_pca, data = New_table, colour = "Population.code")

#sorted by sample names so population tag will line up with the right sample
alleles2 <- alleles[sort(rownames(alleles)),]
sorted_table <- New_table[order(New_table$Sample.name),]
rownames(alleles2) <- sorted_table$Population.code

#tried agglomerative and divisive clustering
par(mfrow=c(2,1))
plot(agnes(dist(alleles2), method = "single"), main="using agnes(brand-complete)",which.plots=2, cex = 0.5)
plot(diana(dist(alleles2)),which.plots=2, main="using diana(brand)", cex=0.5)

#tried decision tree
alleles3 <- as.data.frame(alleles[sort(rownames(alleles)),])
alleles3$Population.code <- sorted_table$Population.code
pop_tree=rpart(factor(alleles3$Population.code)~.,data=alleles3, minsplit=1)

par(mfrow=c(1,1))
plot(pop_tree,margin=0.1, )
text(pop_tree,use.n=T,cex=1.3)

prune.cp <- function(cptable){
  
  CP <- cptable[,1]
  cp <- sqrt(CP * c(Inf, CP[-length(CP)])) 	### inspect the code inside plotcp()!
  xerr <- cptable[,4]
  minrow <- which.min(xerr)
  
  xerr.1se <- sum(cptable[minrow,4:5])
  index <- min(which(xerr < xerr.1se))
  
  cp.1se <- cp[index]
  return(as.numeric(cp.1se) )}

prune.cp (pop_tree$cptable)
pop_tree2 <- prune(pop_tree, cp = prune.cp (pop_tree$cptable) )

par(mfrow=c(1,1))
plot(pop_tree2,margin=.05)
text(pop_tree2,use.n=T)
#would likely make the most sense to make a large dataframe of all the SNPs from all the genes
#possibly combination of genes will make clustering more distinct to population