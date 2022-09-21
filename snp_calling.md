# GATK calling snp
## 前期准备
### 下载参考基因组（推荐hg19或hg38）
#### 1.从NCBI/ensembl/UCSC下载人类基因组  
```
for i in $(seq 1 22) X Y M;  
do echo $i;  
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr${i}.fa.gz;  
done
```
#### 2.建立参考基因组索引文件.fai
```
samtools faidx ref.fa
picard CreateSequenceDictionary R=ref.fa O=ref.dict
```
### 原始测序数据质控
```
fastqc -t 10 ./*_1.fastq.gz -o  ./_1.qc.out
fastqc -t 10 ./*_2.fastq.gz -o  ./_2.qc.out
```
根据质控结果选择直接比对或是进行过滤、去接头等
### bwa比对到参考基因组
#### 1.bwa建立参考序列索引
```
bwa-mem2 index ref.fa -p ref
```
完成之后，会看到类似如下几个以ref.fa为前缀的文件：ref.fa.0123, ref.fa.amb, ref.fa.ann, ref.fa.bwt.2bit.64, ref.fa.pac
#### 2.将测序结果比对到参考基因组
```
bwa mem ref.fa sample_1.fastq.gz sample_2.fastq.gz -R '@RG\tID:sample\tLB:sample\tSM:sample\tPL:ILLUMINA' \
	2>sample_map.log | samtools sort -@ 20 -O bam -o sample.sorted.bam 1>sample_sort.log 2>&1
```
