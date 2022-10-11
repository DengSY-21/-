# GATK calling variants
## 前期准备
### 1.下载参考基因组
#### a.从NCBI/ensembl/UCSC下载人类基因组  
```
for i in $(seq 1 22) X Y M;  
do echo $i;  
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr${i}.fa.gz;  
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
完成之后，会看到类似如下几个以ref为前缀的文件：ref.0123, ref.amb, ref.ann, ref.bwt.2bit.64, ref.pac
#### b.将测序结果比对到参考基因组
```
bwa mem ref.fa sample_1.fastq.gz sample_2.fastq.gz -R '@RG\tID:sample\tLB:sample\tSM:sample\tPL:ILLUMINA' | samtools sort -@ 20 -O bam -o sample.sorted.bam
samtools index sample.sorted.bam
```
#### c.假阳性reads去重复
```
gatk MarkDuplicates -I sample.sorted.bam -O sample.dedup.bam -M sample.dedup_metrics.txt
samtools index sample.dedup.bam
```
#### d.校正碱基质量分数
校正的前提是测序错误碱基的可能性远远大于基因组变异的概率，或者说物种基因变异很小；把所有与参考基因组不一致的碱基视为测序错误导致。所以这里强调是碱基质量分数校准这一步适合于变异概率很少，并且有已知参考变异数据库的物种基因组。因此这一步基本上只适合人类的测序数据。
```
gatk BaseRecalibrator -I sample.dedup.bam -R ref.fa --known-sites dbsnp_146.hg38.vcf.gz Mills_and_1000G_gold_standard.indels.hg38.vcf.gz 
--known-sites 1000G_phase1.snps.high_confidence.hg38.vcf.gz -O sample.recal_data.table
gatk ApplyBQSR --bqsr-recal-file sample.recal_data.table -R ref.fa -I sample.dedup.bam -O sample.recal.bam
```
## Germline variants calling
### 1.GATK
GATK(全称 The Genome Analysis Toolkit)是Broad Institute开发的用于二代重测序数据分析的一款软件，是经常被使用的 variants calling 的软件之一。https://gatk.broadinstitute.org/hc/en-us
```
#先对每个样本生成gvcf文件
gatk HaplotypeCaller -R ref.fa -I sample.recal.bam -ERC GVCF --minimum-mapping-quality 30 -O sample.g.vcf.gz

#然后合并所有样本的gvcf文件
gatk CombineGVCFs -R ref.fa --variant sample1.g.vcf --variant sample2.g.vcf --variant sample3.g.vcf --variant sample4.g.vcf -O SRR4.g.vcf

#将gvcf文件转换为vcf文件
gatk GenotypeGVCFs -R ref.fa -V SRR4.g.vcf -stand-call-conf 5 -O SRR4.vcf
```
### 2.Strelka2
Kim, S., Scheffler, K., Halpern, A.L. et al. Strelka2: fast and accurate calling of germline and somatic variants. Nat Methods 15, 591–594 (2018). https://doi.org/10.1038/s41592-018-0051-x
由 illumina 公司开发，用于突变检测，可以检测 somatic 和 germline ，通常来说，该软件对于小片段的 indel 检测效果比 Mutect2 更好。Strelka2 introduces a novel mixture-model-based estimation of insertion/deletion error parameters from each sample, an efficient tiered haplotype-modeling strategy, and a normal sample contamination model to improve liquid tumor analysis.
```
configureStrelkaGermlineWorkflow.py --bam sample.recal.bam --referenceFasta ref.fa --runDir /germline
```
## Somatic variants calling
### 1.GATK
来自GATK官网的例子，使用Mutect2 call HCC1143肿瘤样本体细胞突变。The command calls somatic variants in the tumor sample and uses a matched normal, a panel of normals (PoN) and a population germline variant resource.
```
gatk Mutect2 \
    -R hg38/Homo_sapiens_assembly38.fasta \
    -I tumor.bam \
    -I normal.bam \
    -tumor HCC1143_tumor \
    -normal HCC1143_normal \
    -pon resources/chr17_pon.vcf.gz \
    --germline-resource resources/chr17_af-only-gnomad_grch38.vcf.gz \
    --af-of-alleles-not-in-resource 0.0000025 \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    -L chr17plus.interval_list \
    -O 1_somatic_m2.vcf.gz \
    -bamout 2_tumor_normal_m2.bam
 ```
### 2.Strelka2
```
# The candidate indel file ${MANTA_ANALYSIS_PATH}/results/variants/candidateSmallIndels.vcf.gz is a recommended best practice but not required. To generate these candidate indels the corresponding configuration for Manta is:
${STRELKA_INSTALL_PATH}/bin/configManta.py --normalBam HCC1187BL.bam --tumorBam HCC1187C.bam --referenceFasta hg19.fa --runDir ${MANTA_ANALYSIS_PATH}
# Somatic analysis
${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py --normalBam HCC1187BL.bam --tumorBam HCC1187C.bam --referenceFasta hg19.fa --indelCandidates ${MANTA_ANALYSIS_PATH}/results/variants/candidateSmallIndels.vcf.gz --runDir ${STRELKA_ANALYSIS_PATH}
```
