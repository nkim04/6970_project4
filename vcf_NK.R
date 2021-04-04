library(vcfR)
library(adegenet)
library(spdep)

##
### DATA EXPLORATION-----
##



#This reads the file into a vcf object
vcf <- read.vcfR( "7.22725889-22732002.ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", verbose = FALSE )
class(vcf)

#Observe the data
#Head of vcf file
head(vcf)
#Head of fixed info
head(getFIX(vcf))
#Head of genotype info
vcf@gt[1:6, 1:4]

#Get the vcf data in tabular format (easier to navigate)
vcf_df <- cbind(getFIX(vcf), vcf@gt)

#To get a matrix of allele counts, turn into a genind object
my_genind <- vcfR2genind(vcf)

#Get a matrix of allele counts
#Individuals are rows, genotypes are columns
#Genotypes are given by [chromosome_number].[start_position].[allele]
#Note that if all the individuals are e.x. "0" at locus 7_22730643, only 7_22730643.0 will be shown 
all_allele_counts <- my_genind@tab
dim(allele_counts)

#You can also look at the number of alleles represented at each locus
my_genind@loc.n.all
range(my_genind@loc.n.all)



##
### DATA EXPLORATION-----
##



#Get allele frequencies, and count those with lower than 0.001 frequency
all_allele_freqs <- apply(all_allele_counts, 2, function(x) {sum(x)/(2*nrow(allele_counts))})
all_allele_freqs[all_allele_freqs < 0.001] #to get rare alleles and their frequencies

allele_freqs <- all_allele_freqs[all_allele_freqs > 0.001]
length(all_allele_freqs); length(allele_freqs) #There are 17 rare alleles and 220 after filtering.

allele_counts <- all_allele_counts[, all_allele_freqs > 0.001]
dim(allele_counts)

rm(all_allele_counts, all_allele_freqs)
