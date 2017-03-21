# mMPA
An R package for implementing marker-assisted mini-pooling with algorithm (mMPA)

## Installation 

The latest version of the `mMPA` package is available at GitHub [taotliu/mMPA](http://github.com/taotliu/mMPA). It requires the `devtools` package to be installed in R. If you do not have `devtools` in your R program, use the code  `install.packages("devtools")` to install the `devtools` package first. Then run the following codes to install the `mMPA` package. 

```R
devtools::install_github("taotliu/mMPA")
library(mMPA)
```

## Example 

The following R code example demonstrates the use of the `mMPA` package. 

### Estimate the average number of assays required by mMPA 

Let us assume that blood samples of `n = 300` HIV+ individuals are collected for HIV viral load (VL) testing. We simulate the VL test results using a Gamma (shape = 2.8, scale = 150) distribution, and generate their corresponding risk scores by adding a uniform random noise to the percentile of VL. The resulting VL has a median of `392` and an interquartile range (IQR) from `224` to `565`; the resulting risk score and VL have a Spearman’s correlation of `0.69`. 

```R
 > n = 300
 > set.seed(100)
 > pvl = rgamma(n, shape = 2.8, scale = 150)
 
 > summary(pvl)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   53      224     392     424     565    1373
 > riskscore = (rank(pvl)/n) * 0.5 + runif(n) * 0.5
 > cor(pvl, riskscore, method = "spearman")
 [1] 0.69
```

We use mMPA to do a pooled VL testing with a pool size of `K = 5`. A total of 60 pools are formed. 

```R
 > # Pool size K is set to 5
 > K = 5
 > # so, the number of pools = 60
 > n.pool  = n/K; n.pool
 [1] 60
``` 
Of course, there are many ways to form pools. Using Monte Carlo simulation, we permute the data `perm_num = 100` time to mimic situations that the individuals came to the clinics in different orders. Thus, different choices of five blood samples are pooled. 

The `mMPA` package includes a function called `pooling_mc(v, s, K, perm_num, method, ...)`, which takes five main arguments as function inputs: Values of test results (`v`), corresponding risk scores (`s`), pool size (`K`), the number of Monte Carlo simulations (`perm_num`), and the method for pooling (which by default use `method = "mMPA"`). The function outputs the total number of VL assays needed for each of the 60 pools from each permutation. 

```R
 > foo = pooling_mc(pvl, riskscore, K, perm_num = 100)
```
 
The output `foo` is a 60x100 matrix, of which each column stores the numbers of VL tests needed by the 60 pools that are formed for each permutation. 

The average number of VL tests needed per pool is then calculated to be `3.35`. 

```R
 > mean(foo)
 [1] 3.35
```

The average number of VL tests needed per individual is then calculated as `0.67`.
```R
> mean(foo)/K
 [1] 0.67
``` 
So the Average number of VL Tests Required per 100 individuals (ATR) is estimated to be 67.  

### Comparison with other pooling algorithm 

If we use mini-pooling (MP) for VL testing, we need an average of `1.19` assays for each individual. 

```R
> foo_mp = pooling_mc(pvl, riskscore, perm_num = 100, method = "minipool")
> mean(foo_mp)
[1] 5.96

> mean(foo_mp)/K
[1] 1.19
> 
```

If we use mini-pooling with algorithm (MPA) (c.f. May et al, 2010), we need `0.79` assay per individual on average. 

```R
> foo_mpa = pooling_mc(pvl, riskscore, perm_num = 100, method = "mpa")
> mean(foo_mpa)
[1] 3.94
> mean(foo_mpa)/K
[1] 0.79
```

Graphically, the efficiency of the three pooling algorithms is illustrated by the following graph. 

```R
boxplot(cbind(MP=apply(foo_mp, 2, mean),
              MPA=apply(foo_mpa, 2, mean),
              mMPA=apply(foo, 2, mean))/K*100,
        ylab = "Average number of assays required per 100 individuals")
```
![](fig/pooling_comp.png)

## Contact

Tao Liu, PhD
tliu@stat.brown.edu
