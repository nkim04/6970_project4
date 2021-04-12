library(tidyverse)
library(vcfR)
library(adegenet)
library(spdep)


###
#### PROJECT PREMISE-----
###


## OBJECTIVE: Use multiple clustering and classification methods to sort super populations
# Classification: Ridge, LASSO, elastic w/ threshold; logistic regression; bagging
# Clustering: Non-hierarchical; hierarchical


# Alcohol digestion
# ALDH2 12:111766933-111817532

# Hair colour
# OCA2 15:27719008-28099342

# Bipolar disorder
# CREB1 2:207529943-207605988

# Solute carrier family 45 member 2
# SLC45A2 5:33944623-33984693


###
#### DATA IMPORT + PREPROCESSING-----
###


getwd()


# Import vcf files
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
  
  
  # Save the resulting df to file to avoid having to re-do this workflow
  write.csv(SNP_counts, file = paste0(gene, "_counts.csv"))
  

}


rm(list=ls())


