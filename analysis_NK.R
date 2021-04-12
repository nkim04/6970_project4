# Preliminary analyses


library(tidyverse)
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

write.csv(counts, file = "counts.csv")

scaled <- scale(counts[2:ncol(counts)])


labeled <- as.data.frame(scaled)
rownames(labeled) <- paste0(counts[,1], seq(nrow(labeled)))


###
#### PCA-----
###


counts_pca <- prcomp(scaled)


# Visualizations
# Scree plot
fviz_eig(counts_pca, main = "Screen Plot of first 10 PCs", addlabels = T) # Visually, there is an elbow between PCA3 and 4
# Contribution plots
fviz_contrib(counts_pca, choice = "var", axes = 1, top = 20) +
  ggtitle("Contribution of genes to \nPC1 (Top 20)") +
  theme(plot.title = element_text(size=10))
fviz_contrib(counts_pca, choice = "var", axes = 2, top = 20) +
  ggtitle("Contribution of genes to \nPC2 (Top 20)") +
  theme(plot.title = element_text(size=10))
# Biplots
fviz_pca_biplot(counts_pca,
                select.var = list(contrib=5),
                label = "var",
                repel = TRUE,
                col.ind = counts$Population,
                addEllipses = F,
                legend.title = "Super-population",
                labelsize = 3,
                title = "PCA Biplot of SNPs")

pca3d(counts_pca, 
      biplot = TRUE, 
      biplot.vars = 1, 
      show.ellipses = T,
      group = counts$Population)

