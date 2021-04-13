library(rpart)
library(rpart.plot)
library(rattle)

gene_counts <- read_csv("counts.csv")

#tried decision tree
pop_tree=rpart(factor(gene_counts$Population)~.,data=gene_counts[,-1],minsplit =1)

par(mfrow=c(1,1))
plot(pop_tree,margin=0.1, )
text(pop_tree,use.n=T,cex=1.3)

summary(pop_tree) ## detailed steps
printcp(pop_tree) ## short summary

fancyRpartPlot(pop_tree, main = "Decision Tree for European and East Asian Superpopulations", sub = "Using ALDH2, CREB1, OCA2 and SLC45A2 SNPs")

prds <- predict(pop_tree,type='class')
table(gene_counts$Population, prds)


