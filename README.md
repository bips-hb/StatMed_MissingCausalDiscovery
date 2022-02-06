# StatMed_MissingCausalDiscovery

This repository contains the files necessary to reproduce all simulation
experiments in Witte J, Foraita R, Didelez V: Multiple imputation and
test-wise deletion for causal discovery with cohort data. Statistics in
Medicine, 2022.

Author: Janine Witte (wittejanine@gmx.de)

The following files contain functions from the micd R-package:
- disCItwd.R
- disMItest.R
- gaussCItwd.R
- gaussMItest.R
- getSuff.R
- make_formulas_saturated.R
- mixCItwd.R
- mixCItest.R
- mixMItest.R

The following file contains an additional convenience function:
- streetparty.R

The following files contain functions from the MVPC repository
(www.github.com/TURuibo/MVPC):
- Tu_CITest.R
- Tu_MissingValuePC.R
- Tu_myfixes.R (modified for better error handling)

The following files contain the graphs and data-generating mechanisms from the
bnlearn Bayesian Network Repository (www.bnlearn.com/bnrepository):
- asia.rds
- ecoli70.rds
- healthcare.rds
- magic-niab.rds
- mehra.rds
- sachs.rds


## How to replicate the simulation experiments:

1) In order to replicate Illustration 1, run CurvePlotCluster.R and
CurvePlot.R. We ran CurvePlotCluster.R on a high-performing computing cluster
with 240 nodes; the runtime was about 13 hours. CurveNullCluster.R contains
code for replicating Illustration 1 but under the null hypothesis of
conditional independence between X and Y given Z (see Online Supplement).

2) In order to replicate Illustration 2, run TilePlotCluster.R and TilePlot.R.
We ran TilePlotCluster.R on a high-performing computing cluster with 240
nodes; the runtime was a bit more than 1 hour.

3) In order to replicate the main simulation study, first run all of
sim_ASIA.R, sim_ECOLI.R, sim_ECOLI_large.R, sim_HEALTHCARE.R, sim_MAGIC.R,
sim_MEHRA.R, sim_SACHS.R. The scripts produce the results in raw format
(estimated adjacency matrices). We ran the scripts on a high-performing
computing cluster with 240 nodes, where they took about 3 weeks to complete.
Next, run all of ana_ASIA.R, ana_ECOLI.R, ana_ECOLI_large.R, ana_HEALTHCARE.R,
ana_MAGIC.R, ana_MEHRA.R, ana_SACHS.R to summarise the results in tabular
format. This can be done on a local machine. plot_summary.R produces Figure 7
of the paper, and ResultsPlot.R produces Figures 8 and 9.


## How to replicate the real data analysis:

Unfortunately, we are not allowed to share the real data. However, the code
for replicating the main analysis is in RealDataAnalysis.R and the code for
replicating the bootstrap analyses is in RealDataBootstrap.R.


## Package versions

We used R versions (on the high-performing computing cluster) and  4.1.0 (on
the local machine).

We used the following versions of the R packages (if two numbers are given,
the first one is for the high-performance computing cluster, the second for
the local machine):
- abind 1.4-5
- backports 1.1.10 / 1.2.1
- bdsmatrix 1.3-4
- BiocGenerics 0.34.0 / 0.38.0
- bnlearn 4.6.1
- boot 1.3-25
- broom 0.7.1 / 0.7.8
- CALIBERrfimpute 1.0.1 / 1.0-5
- class 7.3-17
- clue 0.3.57 / 0.3-59
- cluster 2.1.0 / 2.1.2
- colorspace 2.0-1
- compiler 4.0.2 / 4.1.0
- corpcor 1.6.9
- crayon 1.3.4 / 1.4.1
- data.table 1.13.0
- DEoptimR 1.0-8 / 1.0-9
- DescTools 0.99.41
- dplyr 1.0.2 / 1.0.7
- e1071 1.7.6
- ellipsis 0.3.1 / 0.3.2
- Exact 2.1
- expm 0.999-6
- fansi 0.5.0
- fastICA 1.2-2
- generics 0.0.2 / 0.1.0
- ggplot2 3.3.5
- glue 1.4.2
- ggm 2.5
- gld 2.6.2
- graph 1.68.0 / 1.71.2
- grid 4.0.2 / 4.1.0
- gtable 0.3.0
- igraph 1.2.5 / 1.2.6
- KernSmooth 2.23-17
- ks 1.12.0
- lattice 0.20-41 / 0.20-44
- lifecycle 0.2.0 / 1.0.0
- lmom 2.8
- magrittr 1.5 / 2.0.1
- MASS 7.3.51.6 / 7.3-54
- Matrix 1.2-18
- mclust 5.4.7
- micd 1.2.0
- mice 3.11.0 / 3.13.0
- mipfp 3.2.1
- munsell 0.5.0
- mvtnorm 1.1-1 / 1.1-2
- naniar 0.6.1
- parallel 4.0.2 / 4.1.0
- pcalg 2.6.12 / 2.7-3
- pillar 1.4.6 / 1.6.1
- pkgconfig 2.0.3
- proxy 0.4-25
- purrr 0.3.4
- R6 2.4.1 / 2.5.0
- randomForest 4.6-14
- RBGL 1.64.0 / 1.68.0
- Rcpp 1.0.5 / 1.0.7
- RcppZiggurat 0.1.5 / 0.1.6
- Rfast 2.0.1 / 2.0.3
- Rgraphviz 2.36.0
- rlang 0.4.7 / 0.4.11
- Rmpi 0.6.9
- robustbase 0.93-6 / 0.93-7
- rootSolve 1.8.2.1
- rstudioapi 0.13
- scales 1.1.1
- sfsmisc 1.1-7 / 1.1-11
- snow 0.4.3
- stats4 4.0.2 / 4.1.0
- tibble 3.0.3 / 3.1.2
- tidyr 1.1.2 / 1.1.3
- tidyselect 1.1.0 / 1.1.1
- tools 4.1.0
- tpc 0.5
- utf8 1.2.1
- vctrs 0.3.4 / 0.3.8
- visdat 0.5.3
- weights 1.0.1
- withr 2.4.2
