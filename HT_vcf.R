library(rpart)
library(rpart.plot)
library(rattle)

gene_counts <- read_csv("counts.csv")

#tried decision tree
pop_tree=rpart(factor(gene_counts$Population)~.,data=gene_counts[,-1], minsplit=1)

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
