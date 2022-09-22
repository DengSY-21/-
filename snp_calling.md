# GATK calling germline variants
## 前期准备
### 1.下载参考基因组（推荐hg19或hg38）
#### a.从NCBI/ensembl/UCSC下载人类基因组  
```
for i in $(seq 1 22) X Y M;  
do echo $i;  
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr${i}.fa.gz;  
done
```
#### b.建立参考基因组索引文件.fai
```
samtools faidx ref.fa
picard CreateSequenceDictionary R=ref.fa O=ref.dict
```
### 2.原始测序数据质控
```
fastqc -t 10 ./*_1.fastq.gz -o  ./_1.qc.out
fastqc -t 10 ./*_2.fastq.gz -o  ./_2.qc.out
```
根据质控结果选择直接比对或是进行过滤、去接头等
### 3.bwa比对到参考基因组
#### a.bwa建立参考序列索引
```
bwa-mem2 index ref.fa -p ref
```
完成之后，会看到类似如下几个以ref.fa为前缀的文件：ref.fa.0123, ref.fa.amb, ref.fa.ann, ref.fa.bwt.2bit.64, ref.fa.pac
#### b.将测序结果比对到参考基因组
```
bwa mem ref.fa sample_1.fastq.gz sample_2.fastq.gz -R '@RG\tID:sample\tLB:sample\tSM:sample\tPL:ILLUMINA' | samtools sort -@ 20 -O bam -o sample.sorted.bam
samtools index sample.sorted.bam
```
#### c.比对结果bam文件优化之假阳性reads去重复
```
gatk MarkDuplicates -I sample.sorted.bam -O sample.dedup.bam -M sample.dedup_metrics.txt
samtools index sample.dedup.bam
```
### 4.Varints calling
现在已经有很多团队开发过call variants的软件，这里我们以GATK为例说明
```
#先对每个样本生成gvcf文件
gatk HaplotypeCaller -R ref.fa -I sample.dedup.bam -ERC GVCF --minimum-mapping-quality 30 -O sample.g.vcf.gz

#然后合并所有样本的gvcf文件
gatk CombineGVCFs -R ref.fa --variant sample1.g.vcf --variant sample2.g.vcf --variant sample3.g.vcf --variant sample4.g.vcf -O SRR4.g.vcf

#将gvcf文件转换为vcf文件
gatk GenotypeGVCFs -R ref.fa -V SRR4.g.vcf -stand-call-conf 5 -O SRR4.vcf
```
