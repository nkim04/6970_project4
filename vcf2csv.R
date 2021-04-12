library(tidyverse)
library(vcfR)
library(adegenet)
library(spdep)


###
#### PROJECT PREMISE-----
###


## OBJECTIVE: classification OR clustering into super populations
# Classification: Ridge, LASSO, elastic w/ threshold; logistic regression; bagging
# Clustering: Non-hierarchical; hierarchical

# Lactose intolerance
# LCT 2:135787850-135837195
# MCM6 2:135839626-135876443

# Hair colour
# OCA2 15:27719008-28099342

###
#### DATA IMPORT + PREPROCESSING-----
###


getwd()


# Import vcf files
# LCT
LCT_EAS_file <- "2.135787850-135837195_LCT_EAS.vcf.gz"
LCT_EUR_file <- "2.135787850-135837195_LCT_EUR.vcf.gz"


# MCM6
MCM6_EAS_file <- "2.135839626-135876443_MCM6_EAS.vcf.gz"
MCM6_EUR_file <- "2.135839626-135876443_MCM6_EUR.vcf.gz"


# OCA2
OCA2_EAS_file <- "15.27719008-28099342_OCA2_EAS.vcf.gz"
OCA2_EUR_file <- "15.27719008-28099342_OCA2_EUR.vcf.gz"


genes = c("LCT", "MCM6", "OCA2")


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





