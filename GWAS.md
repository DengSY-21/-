# Getting the data
We will be using publicly available GWAS data from the PennCATH study of genetic risk factors for coronary artery disease. <br>
<a href="https://www.mtholyoke.edu/courses/afoulkes/Data/GWAStutorial/GWASTutorial_Files.zip" title="Data">Data(zip)</a> <br>
<a href="http://www.ncbi.nlm.nih.gov/pubmed/21239051" title="Paper">Paper</a> <br>
Download and unzip the data; you can read the paper as well if you wish. In what follows, I will assume that the unzipped data files are in a folder called data; if you store them somewhere else, change the directory references. <br>
# File formats
The data are given in “PLINK” format, which is the most common format for chip-based GWAS data (as of this writing!). **PLINK** is an open-source whole genome association analysis toolset designed to perform a range of basic large-scale analyses in a computationally efficient manner. It is worth knowing how to use PLINK, although you can also do most of these things in R. <br>
## .fam
This contains information on the subjects:
```
(fam <- fread('data/GWAStutorial.fam'))
#          V1 V2 V3 V4 V5 V6
#    1: 10002  1  0  0  1  1
#    2: 10004  1  0  0  2  1
#    3: 10005  1  0  0  1  2
#    4: 10007  1  0  0  1  1
#    5: 10008  1  0  0  1  2
#   ---                     
# 1397: 11591  1  0  0  2  0
# 1398: 11592  1  0  0  1  2
# 1399: 11593  1  0  0  1  1
# 1400: 11594  1  0  0  1  1
# 1401: 11596  1  0  0  1  0
```
There are 1401 rows, one for each subject. The six colums are: <br>
1. Family ID
2. Individual ID
3. Paternal ID
4. Maternal ID
5. Sex (1=male; 2=female; other=unknown)
6. Phenotype
In this data set, columns 2-4 are unimportant. In general, they are used to specify pedigrees (e.g., subject 3 is the daughter of subjects 1 and 2). In this study, however, none of the subjects are related, so the only column that is important is the first, which records the subject’s unique ID. <br>
Phenotype is typically used to record case-control status or something like that, but it is also quite common to just record clinical/biological information in a separate spreadsheet, which is what was done here. <br>
```
(clinical <- fread('data/GWAStutorial_clinical.csv'))
#       FamID CAD sex age  tg hdl ldl
#    1: 10002   1   1  60  NA  NA  NA
#    2: 10004   1   2  50  55  23  75
#    3: 10005   1   1  55 105  37  69
#    4: 10007   1   1  52 314  54 108
#    5: 10008   1   1  58 161  40  94
#   ---                              
# 1397: 11591   0   2  59  34  44  89
# 1398: 11592   1   1  45  69 101  77
# 1399: 11593   1   1  59  77  27  41
# 1400: 11594   1   1  30  NA  NA  NA
# 1401: 11596   0   1  64 224  35  96
```
As you can see, we’ve got the **FamID** to match this spreadsheet up with the genetic data, the disease status (**CAD=1** means that the subject has coronary artery disease), and some covariates (age, triglycerides, HDL and LDL cholesterol levels). <br>
## .bim
The **.bim** file, by contrast, contains information on the genetic loci (SNPs):
```
(bim <- fread('data/GWAStutorial.bim'))
#         V1         V2 V3        V4 V5 V6
#      1:  1 rs10458597  0    564621 -9  C
#      2:  1 rs12565286  0    721290  G  C
#      3:  1 rs12082473  0    740857  T  C
#      4:  1  rs3094315  0    752566  C  T
#      5:  1  rs2286139  0    761732  C  T
#     ---                                 
# 861469: 10  rs2104725  0  84962606  A  G
# 861470:  7  rs5953712  0 119649910 -9  A
# 861471:  2  rs5951861  0 183267045 -9  A
# 861472:  1 rs12396794  0 159242567 -9  T
# 861473:  3  rs5970564  0 104183552  G  A
```
As you can see, we have 861473 rows here, one for each SNP measured in the study. The columns are:
1. chromosome (1-22, X, Y or 0 if unplaced)
2. rs# or snp identifier
3. Genetic distance (morgans)
4. Base-pair position (bp units)
5. Allele 1 (usually minor)
6. Allele 2 (usually major)
It is pretty common for column 3 to be ignored, as it is here. <br>
So, for example, the file tells us that genetic locus rs12565286 is located 721290 bases into chromosome 1, and that most people have a C there, but some have a G. <br>
## .bed
Finally, the **.bed** file, which has all the data. This is by far the largest of the three files, as it contains the entire 1401 by 861473 matrix of genotype calls for every subject and every locus. To keep things manageable, this file is encoded in a special binary format – i.e., you can’t just read it in through normal means. <br>
To access it, you’ll have to use specialized applications. I’ll discuss two, an R package called **snpStats** and a command-line interface (CLI) called PLINK.
# Software
## snpStats
This is a Bioconductor package. So, you’ll have to install it via **BiocManager** <br>
```
# install.packages('BiocManager')
# BiocManager::install('snpStats')
library(snpStats)
```
To read in data, there is the **read.plink()** function: <br>
```
obj <- read.plink('data/GWAStutorial')
```
The function assumes that all the files have the same base filename, and differ only in their extension. If this is not the case, then you need to specify the filenames for the **.bed**, **.bim**, and **.fam files** separately. <br>
From here, **snpStats** has a lot of functions. For example, here’s a plot (there are 1401 points, one for each subject) of whether the call rate (% of genotype calls that are non-missing) is related to the heterozygosity rate (% of loci that are called AB, as opposed to AA or BB): <br>
```
plot(row.summary(obj$genotypes)[c(1,3)])
```
Feel free to read the **snpStats** documentation and explore for yourself, but one standard thing that one is always interested in is to simply convert various SNPs to a regular numeric matrix so that you can analyze them using standard R tools. For example, let’s do a Fisher’s exact test to see whether CAD is associated with SNP 143:
```
x <- as(obj$genotypes[,143], 'numeric')
fisher.test(drop(x), clinical$CAD)
# 
#   Fisher's Exact Test for Count Data
# 
# data:  drop(x) and clinical$CAD
# p-value = 0.8043
# alternative hypothesis: two.sided
```
A GWAS is then basically just a big loop where we repeat this analysis for every single SNP (although there are of course statistical issues that come up in doing so). <br>
Side note: In general, code like the above is risky, as it assumes that the clinical spreadsheet is in the same order as the **.fam** and **.bed** files. This happens to be the case here:
```
all.equal(rownames(x), as.character(clinical$FamID))
# [1] TRUE
```
But you should get in the habit of explicitly checking for things like this by including lines like this in your code: <br>
```
stopifnot(all.equal(rownames(obj$genotypes), as.character(clinical$FamID)))
```
# Quality control
The first step in any GWAS is to examine the data for potential problems. You don’t want to carry out a GWAS, think you have an exciting result, then discover that it was all just an artifact of bad data. This is a fairly “clean” data set, so it’s not really ideal for showing these steps, but I’ll go through them anyway. There is also a sample data set in the **snpStats** package with some (at least one) bad samples; might be worth checking that one out as well.
```
data(for.exercise)
snps.10
# A SnpMatrix with  1000 rows and  28501 columns
# Row names:  jpt.869 ... ceu.464 
# Col names:  rs7909677 ... rs12218790
```
Most of these QC steps involve calculating summaries at the individual (“row”) level or the SNP (“column”) level: <br>
```
rs <- row.summary(obj$genotypes)
cs <- col.summary(obj$genotypes)
ggbox <- function (X, xlab = "ind", ylab = "values") {
    if (!is.data.frame(X)) X <- as.data.frame(X)
    ggplot2::ggplot(utils::stack(X), ggplot2::aes_string("ind", 
        "values")) + ggplot2::geom_boxplot() + ggplot2::xlab(xlab) + 
        ggplot2::ylab(ylab)
}
```
## Chromosome check
This isn’t exactly a QC step, but extremely helpful to do as a first step when getting any genetic data: what chromosomes are the SNPs on?
```
table(obj$map$chromosome)
# 
#     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17 
# 71038 73717 60565 55675 56178 54176 46391 48299 41110 47930 44213 42124 34262 28054 25900 27591 19939 
#    18    19    20    21    22 
# 26231 11482 22753 12463 11382
```
For the most part, the chromosomes are ordered by size, so chromosomes 1 and 2 are much bigger than (and have many more SNPs than) chromosomes 21 and 22. This particular data set only contains SNPs from the “autosomal” chromosomes (1-22); we do not have any data on X and Y chromosomes. Also, there are no “strange” chromosomes, such as MT (the mitochondrial chromosome), XY (the pseudo-autosomal region of chromosome X), and 0 (SNPs that cannot be mapped to any chromosomes, which can happen for a variety of reasons). Depending on the scientific goal, we may wish to subset our analysis down to include only the autosomal chromosomes, although this isn’t necessary here since it has already been done. <br>
## Missing data
Any SNP with a lot of missing data is probably questionable; these SNPs are often excluded from analysis (although we will talk about other approaches later). Likewise, any sample with lots of missing data suggests that there may be issues with the processing of that sample.
```
ggbox(rs$Call.rate, 'Individuals', 'Call rate')
ggbox(cs$Call.rate, 'SNPs', 'Call rate')
```
Individuals look good – SNPs, on the other hand, there are definitely some SNPs with lots of missing values. A common practice is to exclude SNPs with >5% or >10% missing data. We’ll actually do the subsetting a little later, right now we’re just exploring. <br>
## Minor allele frequency
Minor allele frequency is the percent of alleles that belong to less common category. For example:
```
(Tab <- table(as(obj$genotypes[,143], 'numeric')))
# 
#    0    1    2 
#   17  247 1121
(2*Tab[1] + Tab[2]) / (2*sum(Tab))
#        0 
# 0.101444
cs[143,]$MAF
# [1] 0.101444
```
Excluding SNPs on the basis of minor allele frequency is a bit controversial. It’s done, and it makes sense, but has nothing to do with quality control – there is no reason to think there are any errors in the data. The main justification is statistical: <br>
1. If MAF is low, power is low (i.e., don’t spend multiple testing corrections on tests that are unlikely to find anything anyway)
2. Some statistical methods perform badly with low MAF (e.g., the *chi2*-test)
An appropriate cutoff definitely depends on sample size – the larger the sample, the greater your ability to include rare SNPs. Let’s look at the distributon of MAFs:
```
hist(cs$MAF, breaks=seq(0, 0.5, 0.01), border='white', col='gray', las=1)
```
With a sample size of 1401, I would say a reasonable MAF would be something like 0.001 (0.1%).
```
# How many SNPs would this exclude?
table(cs$MAF < 0.001)
# 
#  FALSE   TRUE 
# 801054  60418
# Would we really learn anything from analyzing a SNP like this?
table(as(obj$genotypes[,62], 'numeric'))
# 
#    0    2 
#    1 1024
```
Finally, it is worth noting that no matter what the sample size, monomorphic SNPs (i.e., SNPs that show no genetic variation whatsoever in the sample) are usually problematic and should always be removed. Some code crashes when monomorphic SNPs are included; even if this weren’t the case, these SNPs cannot possibly be informative in a genome-wide association study.
```
table(cs$MAF == 0)                 # >26000 monomorphic SNPs
# 
#  FALSE   TRUE 
# 835319  26153
obj$map[head(which(cs$MAF==0)),]   # Note that "allele.1' is missing for these SNPs
#            chromosome   snp.name cM position allele.1 allele.2
# rs10458597          1 rs10458597 NA   564621     <NA>        C
# rs307378            1   rs307378 NA  1268847     <NA>        G
# rs12354350          1 rs12354350 NA  1411876     <NA>        C
# rs12563290          1 rs12563290 NA  2028737     <NA>        G
# rs3122920           1  rs3122920 NA  2449711     <NA>        G
# rs897620            1   rs897620 NA  2787707     <NA>        C
```
## Sex check
In general, since we have genetic data on the individuals in the sample, including the X chromosome, we can determine (or at least, estimate) their “genetic sex” and compare that to the sex that is recorded in their clinical information. A discrepancy is very troubling, as it may be the result of a sample being switched or mis-labeled (there are other explanations as well). <br>
In general, though, you have to use PLINK for this. The relevant command is called **--check-sex**. <br>
Note that there are no discrepancies between the sex recorded in **clinical.csv** and the one recorded in the **.fam** file:
```
table(obj$fam$sex, clinical$sex)
#    
#       1   2
#   1 937   0
#   2   0 464
```
## Hardy-Weinberg equilibrium
The <a href="https://en.wikipedia.org/wiki/Hardy-Weinberg_principle" title="Hardy-Weinberg principle">Hardy-Weinberg principle</a> states that under the assumption of random mating, the distribution of genotypes should follow a binomial distribution with probability π equal to the MAF. If this doesn’t happen, this is an indication that either: <br>
1. There was a genotyping error for this SNP, or <br>
2. Mating is not random <br>
In the real world, mating is of course not random, making it difficult to exclude SNPs on the basis of HWE. The usual recommendation is to exclude a SNP only if HWE is hugely violated (e.g., p<10−10 for a test of whether the data follow a binomial distribution).
```
ggbox(cs$z.HWE)  # Mostly near zero, but some huge outliers
# Warning: Removed 26154 rows containing non-finite values (stat_boxplot).
p_hwe <- 2*pnorm(-abs(cs$z.HWE))
table(p_hwe < 10^(-10))
# 
#  FALSE   TRUE 
# 811689  23630
# This seems utterly bizarre -- why would there be so many A/B's, but
# no A/A's or B/B's?  Something is definitely wrong:
table(as(obj$genotypes[,which.max(cs$z.HWE)], 'character'))
# 
#  A/B  B/B   NA 
# 1382    2   17
```
## Heterozygosity
A somewhat similar idea, but applied to individuals instead of SNPs (if an individual had a ton of A/B calls but no A/A or B/B calls, or vice versa, that would indicate something was wrong):
```
ggbox(rs$Heterozygosity)
```
## Pipeline
OK, now that we’ve surveyed all these QC concepts, let’s actually do some filtering. I might choose to do something like this:
```
keep <- cs$MAF > 0.001 &
  cs$Call.rate > 0.9 &
  abs(cs$z.HWE) < 6.5
table(keep)
# keep
#  FALSE   TRUE 
# 108798 752675
```
So, we’re getting rid of about 100,000 SNPs and keeping about 750,000. <br>
Now, let’s actually do the subsetting. **IMPORTANT**: The key thing to remember here is that you need to subset both the **map** object and the genotypes – if you don’t, these objects will no longer match and you will end up with a devastating row mismatch problem.
```
(obj$genotypes <- obj$genotypes[, keep])
# A SnpMatrix with  1401 rows and  752675 columns
# Row names:  10002 ... 11596 
# Col names:  rs12565286 ... rs28729663
obj$map <- obj$map[keep, ]
```
In principle, we might also be throwing some subjects out at this point, but in this particular example, none of the subjects looked questionable. Again, if throwing away subjects, you need to remember to also subset the fam object and the clinical data table. <br>
Let’s save this QC’d data set for future use in downstream analyses:
```
saveRDS(obj, 'data/gwas-qc.rds')
```
# Imputation
To begin, read in the qc (quality controlled) data from earlier step
```
# Load our packages  (same ones mentioned in the data module)
library(snpStats)
library(SNPRelate)
library(data.table)
library(magrittr) # for the pipe operator '%>%'

## Register cores for parallel processing - very helpful if you're on a laptop
library(doSNOW)
registerDoSNOW(makeCluster(4))

# Read in data object created in previous module 
obj <- readRDS('data/gwas-qc.rds')
obj$genotypes
# A SnpMatrix with  1401 rows and  752675 columns
# Row names:  10002 ... 11596 
# Col names:  rs12565286 ... rs28729663
# double check dimensions 
dim(obj$map)
# [1] 752675      6
cs <- col.summary(obj$genotypes) # hold onto this --- will need it later 
```
In a genetics research context, most observations (e.g patients or subjects) will have missing values for at least one SNP. A common method of dealing with missing SNP data is imputation. <br>
## Why impute?
There are two main reasons that one would use imputation: <br>
1. To replace missing SNP values with what these values are predicted to be, based upon a person’s (or subject’s) available SNP values near the loci of interest. For instance, suppose that in a given data set, patient A is missing SNP 123. This value for patient A could be imputed based on the patient’s other SNP values at loci near 123. <br>
2. To infer values for SNPs that were not measured at all for any patients (subjects). This would be the case if one was to merge data from studies that examined different loci. For instance, suppose I am merging data from studies A and B. Suppose further that study A measured SNP 123 for all patients, but study B did not. In the merged data, I may need to impute values for SNP 123 for all patients in the study B data. <br>

For the purposes of this tutorial, let us limit ourselves to scenario (1). <br>

Recall that in the QC step of our analysis, we excluded SNPs with ≥90% missingness. However, there may still be SNPs with some missingness. SNPs that are not missing are described as “called.” The **call rate** is the proportion of genotypes that are called (see **snpStats::col.summary()** documentation for details). Therefore, a call rate of 1 indicates that a SNP has no missing values. <br>
By examining the call rate information, we will first check how many SNPs in our qc’d data set have some missing data:
```
table(cs$Call.rate == 1)
# 
#  FALSE   TRUE 
# 605524 147151
```
