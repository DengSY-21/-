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
This tells us that 147151 SNPs have no missingness, but 605524 still have some missingness (albeit less than 10%.) As a first step, we will try to impute values for these SNPs using the **snp.imputation()** function from **snpStats**. snp.imputation() has numerous options that can be tweaked according to the needs of a specific problem. We will perform a basic imputation for now; see the R documentation for more details.

The package **snpStats** uses a two step imputation procedure. First, the function determines a set of “tag” SNPS. These tag SNPs are used to predict the missing SNP values and to generate prediction rules for the missing SNPs. Second, these prediction rules are applied to the supplied genotype matrix where missing SNP values are imputed.

**N.B** In the case where there is insufficient data or a lack of tagging SNPs, it is possible for the generated prediction rules to fail at yielding predictions. We will see this occur as we go through our example.
## Implementation/tools
To implement imputation, we will use a three step approach:
1. Determine tag SNPs
2. Use tag SNPs to generate prediction rules
3. Apply these prediction rules to our genotype matrix and “fill in the blanks”
### Determine the ‘tag’ SNPs
A SNP is called a ‘tag’ SNP if it is being used to represent (or mark) a specific haplotype. Typically, a tag SNP is in a region of the genome with high linkage disequilibrium. As you will recall, the areas of the genome where there appears to be non-random association of alleles in the population are the areas of our interest. <br>

We want to find these tag SNPs and use them to help us impute missing values. There are many algorithms that can be used to identify tag SNPs - a deep dive into this will take you into computational complexity theory. Check out <a href="https://en.wikipedia.org/wiki/Tag_SNP#Steps_for_tag_SNP_selection" title="the Wikipedia page">the Wikipedia page</a> if you want to take that plunge – for our purposes here, we will use the same function in the **snpStats** package to both identify tag SNPs and generate prediction rules. <br>
### Use tag SNPs to generate prediction rules
As mentioned above, I use the **snp.imputation()** function to both identify the tag SNPs and generate the prediction rules for imputing the missing values:
```
?snp.imputation # check out the help file -- there is a lot here 

# determine tagging SNPs. Note: this can take a few minutes
rules <- snpStats::snp.imputation(obj$genotypes, minA=0)
# NB: minA is a threshold of the amount of existing data needed to impute missing values. Higher minA is a more stringent threshold. Here, we are setting the most loose threshold possible - the default threshold value is 5.
```
### Fill in the blanks
Now that we have tagged important SNPs and created rules for imputation, we can actually implement the imputation with the **impute.snps()** function. This will “fill in the blanks” in our data set, decreasing the number of missing values.
```
rules_imputed_numeric <- impute.snps(rules, obj$genotypes, as.numeric = TRUE)
# returns numeric matrix (see help documentation)

# compare this column summary to the numeric format 
# NB: using `apply()` exhausts memory, but `sapply()` will work: 
call_rates <- sapply(X = 1:ncol(rules_imputed_numeric),
                    FUN = function(x){sum(!is.na(rules_imputed_numeric[,x]))/nrow(rules_imputed_numeric)})
```
We can look at the *R2* values to check the imputation quality. This vignette has additional information about accessing the *R2* values and evaluating imputation quality with snpstats: <br>
<a href="https://www.bioconductor.org/packages/release/bioc/vignettes/snpStats/inst/doc/imputation-vignette.pdf" title="Imputation Vignette">Imputation Vignette</a> <br>

Notice that even after going through the imputation process, there are still 239687 missing values in this data set. This is not unusual for genetics data. It is not uncommon for there to be SNPs that are both missing a notable amount of values and located “far” from surrounding SNPs. In such situations, it is not possible to impute values – we do not know enough to impute a value in these cases. However, we also know that we cannot have any missing values in a regression model (which is where we are headed in our analysis). So, for the missing values that remain after imputation, we can use this case-by-case approach: <br>
1. Does the SNP have $ > 50 %$ missingness? If so, exclude it from the analysis. We do not know enough to impute a value, and there is not enough information in this SNP for us to learn anything about our outcome(s) of interest.
2. Does the matrix of SNP data fit into my computer’s memory? If so, then I can do a simple mean imputation for SNPs with $ %$ missingness. That is, take the mean value of that SNP across all the genotypes, and use this mean value to “fill in the blanks.” I give an example of this in just a bit.
3. If the matrix of SNP data is too large for my computer’s memory, I can use some functions from the package **SNPRelate** to work with the SNP data without storing it as a matrix in my memory (e.g as an object in my global environment)

Let’s see how many SNPs in our example data have $ %$ missingness.
```
# how many of the SNPs with some missingness are missing for <= 50% of observations? 

sum(call_rates >= 0.5)
# [1] 752675
```
So we notice that all of our SNPs with remaining missing data are missing values for no more than half of the patients in the study. I will use the simple mean imputation to address this missingness - no more SNPs need to be eliminated from the analysis.
A couple of notes about this approach:
1. The cutoff of 50% is an arbitrary choice on my part. You could choose 60% or 75% of you wanted to… it may even be best to examine what happens to the results for your specific data set across several cutoff values.
2. Of course, the simple mean imputation only applies when you are talking about a continuous trait. For a categorical outcome, you would need another approach
## When the SNP matrix fits into memory
Again, many statistical methods we may want to apply to these data which cannot handle any missingness. As mentioned above, one simplistic yet reasonable thing to do for these values is to replace them with their HWE expected value (i.e. the average of that SNP).
### Mean imputation
```
# identify which SNPs have missingness
to_impute <- which(call_rates < 1)

# Now, I will try to perform the mean imputation on the numeric matrix 

#' A function for simple mean imputation of continuous SNP data
#' @param j A column of data from a SNP matrix
#' @return j A mean-imputed version of that data
impute_mean <- function(j){
  # identify missing values in a numeric vector
  miss_idx <- which(is.na(j)) 
  # replace missing values with that SNP mean
  j[miss_idx] <- mean(j, na.rm = TRUE) 
  
  return(j)
}
```
```
# Create the fully imputed matrix - I am about to "fill in" all the blanks with 
#   the function I just wrote
fully_imputed_numeric <- rules_imputed_numeric

# Apply function to the columns (SNPs) where there is missingness
fully_imputed_numeric[,to_impute] <- apply(X = rules_imputed_numeric[,to_impute],
                                           MARGIN = 2,
                                           FUN = impute_mean)

# now, for the sake of saving memory, remove the rules_imputed_numeric - won't need this again 
rm(rules_imputed_numeric)
```
A brief note for those running this on macOS: when I first tried to run this **apply(...)** statement on my MacBook Pro, I got the error **vector memory exhausted (limit reached?)** several times. To address this issue, here is what worked for me:
1. In the console, run usethis::edit_r_environ() to open the .Renviron file
2. Edit that file by changing the R_MAX_SIZE argument to be something large, like 100Gb (‘large’, of course, is relative to the computer in use)
3. Close the file, restart R, and try running the code again.
Now, I can hold onto the **fully_imputed_numeric** matrix and use this to examine the data for population structure (see the next module).

We can check missingness using base **R** functions in multiple chunks **R** can handle. This can take a while, but it will reassure us that we are ready to move on to something like principal component analysis (see next section).
```
# make sure there are no missing values remaining / count missing values
missing <- 0
chunks <- ceiling(nrow(fully_imputed_numeric) / 100) # I'm breaking this up using 100 based on on trial and error but this can be tweaked.
start <- 1
for (i in 1:chunks){
  stop <- min(i*100, nrow(fully_imputed_numeric))
  missing <- missing + sum(is.na(fully_imputed_numeric[start:stop,]))
  start <- stop + 1
}
missing # should be 0



# Check for Inf values 
inf <- 0
start <- 1
for (i in 1:chunks) {
  stop <- min(i*100, nrow(fully_imputed))
  inf <- inf + sum(!(is.finite(fully_imputed[start:stop,])))
  start <- stop + 1
}

inf # should be 0 
```
Let’s save this fully imputed data set for future use in downstream analyses:
```
# saveRDS(obj, 'data/gwas-imp.rds')
saveRDS(fully_imputed_numeric, 
        "data/fully_imputed_numeric.rds")
```
# SNP testing (marginal approach)
## Set up
In this module, we work through an example of a GWAS analysis using a marginal approach. Our approach is ‘marginal’ in the sense that we are testing the relationship between the outcome and the genetic data by going one SNP at a time. This is (by far) the most pervasive approach in the existing literature.

This module does not require that you have worked through the previous modules; however, I use language from the quality control, imputation, and population structure modules throughout my explanation of this marginal approach.

As always, we begin by loading the libraries with the tools we need for analysis.
```
library(data.table)
library(magrittr)
library(qqman)
library(snpStats)
library(dplyr)
```
Next, we load the data. Start by loading the clinical data, as this has the outcome (coronary artery disease (‘CAD’)) we need for our models.
```
clinical <- fread("data/GWAStutorial_clinical.csv")
# str(clinical) # if you need to remind yourself of what is here
```
Next, we need to load the genetic data. For this section, we will work with the quality controlled (QC’d) data from the **SNPRelate** package (see module 1). We will also need the “.bim” file from the original data for making plots.
```
# load QC'd data:
qc_dat <- readRDS('data/gwas-qc.rds')
# load the bim file
bim <- fread('data/GWAStutorial.bim')
```
If you completed the population structure module, you should load the principal components as well - we need them for our analysis. If you did not work through the population structure module, you can skip this step.
```
# load principal components 
PCs <- readRDS(file = "data/PCs_base.rds") %>%
  as.data.frame()

names(PCs) <- paste0("PC", 1:ncol(PCs))
```
## Implement tests
We begin our analysis with a logistic regression model which uses the predictors sex and age to study the binary outcome CAD. With the function **snp.rhs.tests()**, we will perform a large set of logistic regression tests and save the p-values in a vector called **assoc_p_vals**. In addition to adjusting for sex and age, we will adjust for the first 4 principal components, as per the observations we made about the population structure in our data.
```
assoc_test <- snpStats::snp.rhs.tests(
  formula = clinical$CAD ~ clinical$sex + clinical$age + 
    PCs$PC1 + PCs$PC2 + PCs$PC3 + PCs$PC4,
  family   = "binomial",
  link = "logit", # echoes the default settings
  data = qc_dat$fam,
  snp.data = qc_dat$genotypes
  )

assoc_p_vals <- p.value(assoc_test)
```
Alternatively, one could implement this testing using the imputed data. There is an optional **rules** argument in the **snp.rhs.tests()** function that allows the user to include the rules used for imputation. I will forgo doing this here – simply want to point out that this approach is possible.
## Visualize results
There are two methods for visualizing the results of a GWAS analysis: **qq plots** and **Manhattan plots**. Both of these tools are meant to highlight the genetic variants (SNPs) with the smallest corresponding p-values.

Due to their magnitude, p-values from GWAS studies are typically illustrated on the log-transformed scale.
### Using a qqplot
As a first look at our results, we will examine a qq-plot. This plot compares the observed p-values from our SNP tests with the p-values we would expect if there were no associations. To make our qq-plot, we will use the **qq()** function from the package **qqman**.
```
qq_plot <- qqman::qq(assoc_p_vals)
(qq_plot)
# NULL
```
If the observed p-values follow exactly the distribution we would expect under the global null hypothesis, then all points will be plotted on the red line.

We notice in this qq-plot that some of the p-values (plotted as points) diverge from the red line to the left. This indicates that some observed values are below what is expected. From a statistical point of view, this is pretty promising. Having some p-values that are smaller than expected indicates that there are likely to be significant results in the data.
### Using a Manhattan plot
A <a href="https://en.wikipedia.org/wiki/Manhattan_plot" title="Manhattan plot">Manhattan plot</a> makes the smallest p-values “pop” by plotting all p-values so that the smallest ones are the highest. We can create a Manhattan plot of our results using the **manhattan()** function (from the aforementioned **qqman** package)
```
# format data to have readable labels. 
manh_data <- data.frame(
    SNP = assoc_test@snp.names, # NB: for S4 objects in R, use the "@" to access items
    P = assoc_p_vals
  ) %>%
  left_join(bim, by = c("SNP" = "V2")) %>%
  rename(
    CHR = V1, 
    BP = V4
  ) # recall that 'bim' files have a standardized format, so the column order is 
# always the same 

manh_plot <- manhattan(manh_data, ylim = c(0, 20))

(unlist(manh_plot))
#   xpd 
# FALSE
```
Here, the x-axis of the plot is divided into “bins”, where each bin represents a chromosome. The y axis shows the negative, log-transformed p-values. Each of the p-values from our analysis is plotted in the bin corresponding to the chromosome of the gene in that particular test, and the height of the point directly correlates with the significance of the test.

The goal of Manhattan plots are to help identify SNPs (or a region of SNPs) that are associated with a phenotype of interest. The blue and red lines represent values for **−log10** transformed p-values at two specified thresholds of “significance.” When we do have SNPs with a p-value that exceeds these lines, we are often interested in the one that is the highest in a given region.

I could improve the Manhattan plot above by adding annotations which indicate the SNPs with the smallest p-values:
```
# NB:  5 × 10e−8 is a common threshold for significance in GWAS studies, 
#   whereas 5 x 10e-6 is a common threshold for "suggestive" results
# signif_threshold <- 5e-8 # this would be a more stringent alternative 
suggest_threshold <- 5e-6 
manh2 <- manhattan(manh_data, ylim = c(0, 10), annotatePval = suggest_threshold)

(unlist(manh2))
#  xpd 
# TRUE
```
Based on the Manhattan plots, one region of interest could be around rs9632884 on Chromosome 9 - this SNP is above the suggestive threshold (indicated by the blue line), and there are a lot of other notable results clustered near this SNP. If I wanted to highlight a region of interest, I could do this to highlight the specific genes in this region, or <a href="https://www.genome.gov/genetics-glossary/Locus" title="locus">locus</a>. This will also allow us to practice the final piece of functionality from **qqman**:
```
bp_center <- manh_data %>% # NB: 'bp' stands for 'base pair'
  filter(SNP == "rs9632884") %>%
  pull(BP)

bp_range <- c(-1, 1) * 100000 + bp_center

snps_highlight <- manh_data %>%
  filter(BP >= bp_range[1], BP <= bp_range[2], CHR == 9) %>%
  pull(SNP)

manh3 <- manhattan(
  manh_data, 
  ylim = c(0, 10), 
  annotatePval = .000005, 
  highlight = snps_highlight
)

(unlist(manh3))
#  xpd 
# TRUE
```
# SNP testing (joint approach)
The objective of this module is to test for significant SNPs using a joint approach. For this, we will use the **penalizedLMM** R package, which is available <a href="https://github.com/areisett/penalizedLMM" title="on GitHub">on GitHub</a>. Once again, we will begin by loading the necessary libraries:
```
library(data.table)
library(magrittr)
library(qqman)
library(snpStats)
library(dplyr)

# devtools::install_github("areisett/penalizedLMM")
library(penalizedLMM)
```
Next, we load the data. Start by loading the clinical data, as this has the outcome (coronary artery disease (‘CAD’)) we need for our models.
```
clinical <- fread("data/GWAStutorial_clinical.csv")
# str(clinical) # if you need to remind yourself of what is here
```
We also need to load the genetic data. For this section, we will work with the quality controlled (QC’d) data from the **SNPRelate** package (see module 1). We will also need the “.bim” file from the original data for making plots. Finally, we need our design matrix *X* with no missing values (this is the *X* we obtained by our imputation procedures).
```
# load QC'd data:
qc_dat <- readRDS('data/gwas-qc.rds')
# load the bim file
bim <- fread('data/GWAStutorial.bim')
# load design matrix of SNP data (you would need to complete the imputation module first)
X <- readRDS(file = "data/fully_imputed_numeric.rds")
```
If you completed the population structure module, you should load the principal components as well - we need them for our analysis. If you did not work through the population structure module, you can skip this step.
```
# load principal components 
PCs <- readRDS(file = "data/PCs_base.rds") %>%
  as.data.frame()

names(PCs) <- paste0("PC", 1:ncol(PCs))
```
## Constructing the model
```
# the CAD outcome is binary (0/1), so it may not be the best place to begin for the 
#   joint model example 
joint_model <- plmm(X = X,
                    y = clinical$CAD,
                    intercept = FALSE)
```
## Examining the results
## Comparing the results of the marginal and joint approaches
