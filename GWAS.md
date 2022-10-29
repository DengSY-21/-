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
1.Family ID
2.Individual ID
3.Paternal ID
4.Maternal ID
5.Sex (1=male; 2=female; other=unknown)
6.Phenotype
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
1.chromosome (1-22, X, Y or 0 if unplaced)
2.rs# or snp identifier
3.Genetic distance (morgans)
4.Base-pair position (bp units)
5.Allele 1 (usually minor)
6.Allele 2 (usually major)
It is pretty common for column 3 to be ignored, as it is here. <br>
So, for example, the file tells us that genetic locus rs12565286 is located 721290 bases into chromosome 1, and that most people have a C there, but some have a G. <br>
## .bed
Finally, the **.bed** file, which has all the data. This is by far the largest of the three files, as it contains the entire 1401 by 861473 matrix of genotype calls for every subject and every locus. To keep things manageable, this file is encoded in a special binary format – i.e., you can’t just read it in through normal means. <br>
To access it, you’ll have to use specialized applications. I’ll discuss two, an R package called **snpStats** and a command-line interface (CLI) called PLINK.
# Software
## snpStats
