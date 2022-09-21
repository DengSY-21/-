# GATK calling snp
## 前期准备
### 下载参考基因组（推荐hg19或hg38）
从NCBI/ensembl/UCSC下载人类基因组  
for i in $(seq 1 22) X Y M;  
do echo $i;  
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr${i}.fa.gz;  
done  
### bwa比对到参考基因组
