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
# OCA2 15:135839626-135876443


###
#### DATA IMPORT + PREPROCESSING-----
###


getwd()


# Import vcf files
# LCT
LCT_EAS_vcf <- read.vcfR( "~/Downloads/2.135787850-135837195_LCT_EAS.vcf.gz", verbose = FALSE )
LCT_EUR_vcf <- read.vcfR( "~/Downloads/2.135787850-135837195_LCT_EUR.vcf.gz", verbose = FALSE )


# MCM6
MCM6_EAS_vcf <- read.vcfR( "~/Downloads/2.135839626-135876443_MCM6_EAS.vcf.gz", verbose = FALSE )
MCM6_EUR_vcf <- read.vcfR( "~/Downloads/2.135839626-135876443_MCM6_EUR.vcf.gz", verbose = FALSE )


# OCA2
OCA2_EAS_vcf <- read.vcfR( "~/Downloads/15.135839626-135876443_OCA2_EAS.vcf.gz", verbose = FALSE )
OCA2_EUR_vcf <- read.vcfR( "~/Downloads/15.135839626-135876443_OCA2_EUR.vcf.gz", verbose = FALSE )


# Optional: visually inspect the data (shown for one file)
#head(LCT_EAS_vcf) # View different sections
#head(getFIX(LCT_EAS_vcf)) # View fixed information
#LCT_EAS_vcf@gt[1:6, 1:4] # View genotype information
#EAS_LCT_df <- cbind(getFIX(LCT_EAS_vcf), LCT_EAS_vcf@gt) # Optional: tabular form



# Get matrices of allele counts
# LCT
LCT_EAS_counts <- vcfR2genind(LCT_EAS_vcf)@tab; dim(LCT_EAS_counts) #504 ind.; 1624 SNPs
LCT_EUR_counts <- vcfR2genind(LCT_EUR_vcf)@tab; dim(LCT_EUR_counts) #503 ind.; 1549 SNPs


# There are different numbers of alleles in the populations - let's see which ones overlap, and put them in one df
# LCT
length(which(colnames(LCT_EAS_counts) %in% colnames(LCT_EUR_counts))) # There are 1385 overlapping SNPs
length(which(rownames(LCT_EAS_counts) %in% rownames(LCT_EUR_counts))) # No individuals overlap
LCT_counts <- bind_rows(as.data.frame(LCT_EAS_counts), as.data.frame(LCT_EUR_counts)) %>% 
  .[, colSums(is.na(.)) == 0]; dim(LCT_counts)  # 1007 individuals, 1385 SNPs


# Get allele frequencies and filter out the rare (<0.001) ones from our original
# Also filter out frequencies equal to 1 (same across populations)
LCT_freq <- apply(LCT_counts, 2, function(x) {sum(x)/(2*nrow(LCT_counts))})
table(LCT_freq[LCT_freq < 0.001]) # Three are rare
table(LCT_freq[LCT_freq == 1]) #771 are the same across populations
LCT_counts <- LCT_counts[, names(LCT_freq[LCT_freq > 0.001 & LCT_freq != 1])] %>% 
  cbind(Population = c(rep("EAS", nrow(LCT_EAS_counts)), rep("EUR", nrow(LCT_EUR_counts))), .) # Add population labels as the first column


rm(LCT_EAS_vcf, LCT_EAS_counts, LCT_EUR_vcf, LCT_EUR_counts, LCT_freq)


# Find out how many alleles exist for each SNP
str_sub(colnames(LCT_counts[2:ncol(LCT_counts)]), start= -1) %>% table()
# There are 536 SNPs that report only the reference allele ("0")
# There are 133 SNPs that report the reference AND alternate allele ("1")
# Only one SNP is multi-allelic ("2")


# For bi-allelic SNPs, remove their alternate alleles ("1") - this is redundant information
# For the multi-allelic SNP, keep the "0" reference allele, and the "2" alternate allele - it is still informative
# In other words, remove all SNPs ending with ".1"
colnames(LCT_counts)[grep(".2$", colnames(LCT_counts))] # For reference, 2_135834838.2
LCT_counts <- LCT_counts[, -grep(".1$", colnames(LCT_counts))]


# Save the resulting df to file to avoid having to re-do this workflow
write.csv(LCT_counts, file = "LCT_counts.csv")






