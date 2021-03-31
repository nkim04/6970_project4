# BINF6970: Assignment 4
Due Monday, April 12, 11:59 pm 



# Analysis of Genotype Data Available from the 1000 Genomes Project



Task Description

1. Use genotype data that is available in the 1000 Genomes project.

2. Use various statistical analysis tools of your choice as covered by the end this course: linear regression model; classification methods such as logistic regression, classification trees, random forests; variable selection methods such as LASSO-type selection and using decision trees-based pruning; PCA and clustering analysis: hierarchical and non-hierarchical clustering. You are welcome to combine other methods not covered in this course with the ones you choose to apply form the list above.

3. Make up your own objective of the study by exploring the IGSR website. Examples of objective include:

i)	Cluster a sample of individuals according to their genetic information.
ii)	Perform gene-based association analysis with some traits of interest such as Skin color, ethnicities, etc.
iii)	Classify individuals to different ethnic group according to their genetic information.

Page limit: 5 pages, including title, names, abstract, main text, references. Figures and tables that you choose to use should appear in a natural flow in the main text.

Some Practical Hints:

a) Remove SNPs that are really rare across populations: <0.001.
b) You can back check those equilibriums covered in the course: HWE, LD structure, etc, in your data analysis.
c) Quality is more important than quantity.
d) Please limit the number of SNPs in your analysis; don’t do the whole genome/genome-wide analysis. You may focus on a few genes/regions.

The grading of your group project will be according to the following grading scheme, where each area is worth 5 marks that make up a total of 40 marks.

1. Background and objective
2. Details and description of data
3. Detail of data analysis methods
4. Summary of results
5. Conclusions and discussions
6. Citations (appropriateness) and bibliography (accuracy and format of the reference)
7. Appendix of R-code
8. Quality of the presentation (grammar, spelling, and coherence)



# Additional Guide



Getting 1000 Genomes Project data.

Important website: http://www.internationalgenome.org

IGSR: The international Genome Sample Resource (since January 2015), support the ongoing 1000 Genomes Project.

Three main aims:

1. Ensure maximal usefulness and relevance of the existing 1000 Genomes data resources
2. Extend the resource for the existing populations
3. Expand the resource to new populations

Sampling process:

1. The samples should be collected from adults.
2. In order to achieve 100 unrelated cell lines from a population, it is recommended to collect primary material from at least 130 unrelated individuals. Data from mother-father-adult child trios are valuable for assessing data quality and producing haplotypes, so collecting samples in trios and establishing cell lines from the offspring is strongly encouraged.

Data production

1. The samples (including from any adult children in trios) should be genotyped across the genome using a high density genotyping array with over 500,000 markers.
2. At least 90 unrelated samples should be sequenced. This sequencing can follow two different strategies. Low coverage whole genome sequencing (minimum 3x aligned non duplicated coverage per sample) and exome sequencing (at least 70% of exome target bases covered to 30x of higher), this strategy was used by the main 1000 Genomes Project, or high coverage whole genome sequencing (minimum 30x aligned non duplicated coverage per sample).


The 1000 Genomes Project (2008-2015), the IGSR was set up to support the 1000 Genomes Project.

The 1000 Genome adapted the HAPMAP project that originally has samples from  three populations: Caucasian (CEU), African (YRI), East Asian. Each population has a sample of size 90 consisting 30 tio families: Mom, dad, child. C samples were collected from Central European descendancy in Utah, US, African samples were collected from Yoruba ethnic group, mainly in western Nigeria, Africa, East Asian samples were collected from Chinese Han in Beijing, China (45) and Japaneses in Tokyo, Japan (45).

Current 1000 Genome project covers a much larger number of populations.

Six super-populations: African, American,East Asian, European, South Asian. Each of these superpopulation consists of multiple subpopulations that comprise a total of 26 populations.


Example, I am interested in a gene, MC1R, that involves in determining the skin color or hair color of human. To find information regarding MC1R gene use the NCBI webpage
https://www.ncbi.nlm.nih.gov

From there, search the MC1R under the gene category.

https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/

https://www.ncbi.nlm.nih.gov/gene/4157

You will find that the MC1R (melanocortin 1 receptor) gene in Human is in chromosome 16 in the location of 89984287 to 89987385 bp and so, has a length about 3098bp.

Now, we can use the data slice to extract the SNPs variants data in the MC1R gene using the following link
http://grch37.ensembl.org/Homo_sapiens/Tools/DataSlicer?db=core

For example, extract the SNPs in the MCIR genes for both the FIN (Finnish) and BEB (Bengali) subpopulations, in which, FIN is subpopulation within European and BEB is a subpopulation within South Asia.



-----------------------------------------------------



Other website that I had been to, information that I had collected.

Skin Color:
https://en.wikipedia.org/wiki/Human_skin_color

Genes involved in determining the skin tone among people may include: MC1R, KITLG, ASIP, SLC24A5, SLC24A2 (or MATP), TYR, OCA2.

Hair color:
https://ghr.nlm.nih.gov/primer/traits/haircolor

Genes involved in determining the hair color among people may include: ASIP, DTNBP1, GPR143, HPS3, KITLG, MC1R, MLPH, MYO5A, MYO7A, OCA2, SLC45A2, SLC24A5, TYRP1, TYR, ERCC6, GNAS, HERC2, IRF4, OBSCN, SLC24A4, TPCN2, and MITF



——————————
R function

library(MASS)

#open the help page of the function chisq.test.
>?chisq.test
> x=matrix(c(118, 78, 82, 122), ncol=2)
> x
     [,1] [,2]
[1,]  118   82
[2,]   78  122
> chisq.test(x)

	Pearson's Chi-squared test with Yates' continuity
	correction

data:  x
X-squared = 15.2161, df = 1, p-value = 9.588e-05

>chisq.test(x)$p.value
[1] 9.588316e-05




