install.packages("vcfR")
library(vcfR)

#This reads the file into a vcf object
vcf <- read.vcfR( "~/Downloads/7.22725889-22732002.ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", verbose = FALSE )
class(vcf)

#I tried this out, but it doesn't seem right. Apparently adding annotations and sequences help interpret this plot, but I'm thinking it's a dead end.
chrom <- create.chromR(name='Supercontig', vcf=vcf)
plot(chrom)

#This resources was helpful to understand vcf format in R a bit better: https://cran.r-project.org/web/packages/vcfR/vignettes/vcf_data.html
#Head of vcf file
head(vcf)

#Head of fixed info
head(getFIX(vcf))

#Head of genotype info
vcf@gt[1:6, 1:4]
