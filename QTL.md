# QTL分析
一般使用R/qtl <br>
Broman KW, Wu H, Sen S, Churchill GA (2003) R/qtl: QTL mapping in experimental crosses. Bioinformatics
19:889-890 <br>
## Example : Hypertension
```
install.packages("qtl")
library(qtl)
# 加载数据
data(hyper)
# getting summary information on the data
summary(hyper)
nind(hyper)
nphe(hyper)
nchr(hyper)
totmar(hyper)
nmar(hyper)
```
The function drop.nullmarkers may be used to remove markers that have no genotype data (such as the marker on chr 14). A call to totmar will show that there are now 173 markers (rather than 174, as there were initially).
```
hyper <- drop.nullmarkers(hyper)
totmar(hyper)
```
Estimate recombination fractions between all pairs of markers. This also calculates LOD scores for the test of H0: r = 1/2.
```
hyper <- est.rf(hyper)
```
Re-estimate the genetic map (keeping the order of markers fixed)
```
newmap <- est.map(hyper, error.prob=0.01)
```
We now turn to the identification of genotyping errors.
```
hyper <- calc.errorlod(hyper, error.prob=0.01)
```
QTL mapping <br>
The core of R/qtl is a set of functions which make use of the hidden Markov model (HMM) technology to calculate QTL genotype probabilities, to simulate from the joint genotype distribution and to calculate the most likely sequence of underlying genotypes (all conditional on the observed marker data). <br>
The function calc.genoprob calculates QTL genotype probabilities, conditional on the available marker data. These are needed for most of the QTL mapping functions. The argument step indicates the step size (in cM) at which the probabilities are calculated, and determines the step size at which later LOD scores are calculated.
```
hyper <- calc.genoprob(hyper, step=1, error.prob=0.01)
out.em <- scanone(hyper)
out.hk <- scanone(hyper, method="hk")
hyper <- sim.geno(hyper, step=2, n.draws=16, error.prob=0.01)
out.imp <- scanone(hyper, method="imp")
```
The function scantwo performs a two-dimensional genome scan with a two-QTL model. For every pair of positions, it calculates a LOD score for the full model (two QTL plus interaction) and a LOD score for the additive model (two QTL but no interaction). This be quite time consuming, and so you may wish to do the calculations on a coarser grid.
```
hyper <- calc.genoprob(hyper, step=5, error.prob=0.01)
out2.hk <- scantwo(hyper, method="hk")
# One can also use method="em" or method="imp", but they are even more time consuming.
```
