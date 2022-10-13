# 比较基因组学
## 同源基因的查找
### OrthoFinder
OrthoFinder是比较基因组学中的实用的，运行快速，准确的全面的工具。它的主要功能是，找到了正交群和直系同源物，推断出所有正交群的根基因树，并识别那些基因树中的所有基因重复事件。它还为所分析的物种推断出有根的物种树，并将基因重复事件从基因树比对到物种树的分支中。<br>
正式运行Orthofinder，相当简单的操作,-f 输入目录，里面包含你需要运行的蛋白质fasta文件， -t 所用到的CPU数目。
```
orthofinder -f Mycoplasma/ -t 2
```
## 多序列比对
### MUSCLE
MUSCLE(Multiple Protein Sequence Alignment)是在2004年公布的一款蛋白质水平多序列比对的开源软件，在速度和精度上都优于ClustalW。比对速度快。因此在进行多序列比对的时候，大多数情况下可以优先使用Muscle。例如Mega等软件里面也集成了muscle的多序列比对。Muscle同样可以用于DNA的多序列的比对。使用起来十分方便，大多数情况下用户只需要指定输入输出文件即可。它是适合于多序列比对，多个序列之间具有同源关系，并且具有同一方向，muscle需要对序列进行拉伸，例如适合同一个或多个看家基因、16s等放在一起比对。我们同样也可以将多样品的SNP结果连接起来进行多序列比对，做系统发育分析。所以对于我们使用muscle主要需要做的就是将输入文件格式化为满足muscle输入即可，主要就是fasta格式。<br>
-clw 输出CLUSTALW 格式的结果
```
muscle -in multi.fasta -out ex1.mfa
muscle -in multi.fasta -out ex1.clw -clw
```
## 调取保守区域，并收尾连接，形成supergene
### Gblocks
Gblocks它可以消除DNA或蛋白质序列中排列不一致的位置和不同的区域。这些位置可能不是同源的，或者可能已经被多个取代所饱和，在进行系统发育分析之前消除比较方便。Gblocks选择块的方式与通常手工操作的方式类似，但遵循可重复的一组条件。对于缺少大段连续的非保守位置、间隙位置缺乏或密度低、侧翼位置保守程度高的情况，选择的块体必须满足一定的要求，使最终的对齐更适合于系统发育分析。Gblocks输出几个文件以可视化所选的块。<br>
-t=p/d/c  序列类型p/d/c分别对应蛋白质、DNA、密码子 <br>
-b1   保守位置的最小序列数(50%的序列数+ 1)，应该大于序列数的一半，小于或等于序列总数的任意整数 <br>
-b2   最小数量的序列为侧翼位置，(序列数的85%)等于或大于保守位置序列的最小数目的任意整数 <br>
```
Gblocks proteins.fasta -t=p   # 使用蛋白质序列
```
## 进化树构建
### RAxML
RAxML (Random Axelerated Maximum Likelikhood) 能使用多线程或并行化使用最大似然法构建进化树。<br>

```
raxmlHPC-PTHREADS-SSE3 -f a -m GTRGAMMA -p 12345 -x 12345 -# 10 -s 20k.phy -n chr001.raxml -T 30
-f a
此参数用于选择 RAxML 运算的算法。可以设定的值非常之多。 a 表示执行快速 Bootstrap 分析并搜索最佳得分的 ML 树。
-x 12345
指定一个 int 数作为随机种子，以启用快速 Bootstrap 算法。
-p 12345
指定一个随机数作为 parsimony inferences 的种子。
-# 10
指定 bootstrap 的次数。
-m GTRGAMMA
指定核苷酸或氨基酸替代模型。 "PROT" 表示氨基酸替代模型，“GTR”表示碱基替代模型； GAMMA 表示使用 GAMMA 模型； X 表示使用最大似然法估计碱基频率。
-s 20k.phy
指定输入文件。phy 格式的多序列比对结果。软件包中包含一个程序来将 fasta 格式转换为 phy 格式。也可以通过Tassel或者Mega转换格式：vcf-phylip
-n chr001.raxml
输出文件的后缀为 .chr001.raxml 。
-T 30
指定多线程运行的 CPUs 。
```
## 分化时间分析 divergence time
### PAML
PAML是一个用最大似然法来对DNA和蛋白质序列进行系统发育分析的软件包，由伦敦大学的杨子恒教授开发并更新。PAML的功能很多，但目前主要用于计算密码子同义替换和非同义替换的比率omega值，从而预测氨基酸序列在进化过程中所受的选择压力。<br>
## 基因扩张收缩分析
