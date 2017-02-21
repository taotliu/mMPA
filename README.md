# mMPA
An R package for implementing marker-assisted mini-pooling with algorithm (mMPA)

## Installation 

The latest version of the `mMPA` package is available at GitHub [taotliu/mMPA](http://github.com/taotliu/mMPA). It requires the `devtools` package to install it in R. If you do not have `devtools` in your R program, use the code  `install.packages("devtools")` to install the `devtools` package first. Then run the following codes to install the `mMPA` package. 

```R
devtools::install_github("taotliu/mMPA")
library(mMPA)
```

## Example 

The following R code example demonstrates the use of the `mMPA` package. 

### Estimate the average number of assays required by mMPA 

Let us assume that blood samples of `n = 300` HIV+ individuals are collected for HIV viral load (VL) testing. We simulate the results of VL tests using a Gamma (shape = 2.8, scale = 150) distribution, and generate their corresponding risk scores by adding a uniform random noise to the percentile of VL. The resulting VL has a median of `392` and an inter-quantile from `224` to `565`, and the resulting risk score and VL have a Spearman’s correlation of `0.69`. 

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

We use mMPA to do a pooled VL testing with a pool size of `K = 5`. 
A total of 60 pools are formed. 

```R
 > # Pool size K is set to 5
 > K = 5
 > # so, the number of pools = 60
 > n.pool  = n/K; n.pool
 [1] 60
``` 
Of course, there are many ways to form pools. We use Monte Carlo simulation and permute the data `perm_num = 100` time to mimic situations that the individuals came to HIV clinics in different orders. Thus, different choices of five blood samples are pooled. 

The `mMPA` package includes a function called `pooling_mc(v, s, K, perm_num, vf_cut, ...)`, which takes four main arguments as input: Values of test results (`v`), corresponding risk scores (`s`), pool size (`K`), the number of Monte Carlo simulations (`perm_num`), and the method for pooling (which by default use mMPA). The function outputs the total number of VL assays needed for each of the 60 pools from each permutation. 

```R
 > foo = pooling_mc(pvl, riskscore, K, perm_num = 100)
```
 
The output `foo` is a 60x100 matrix, where each column stores the numbers of VL tests needed by the 60 pools that are formed for each permutation. 

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
So the average number of VL tests needed per 100 individuals (ATR) is estimated to be 67.  

### Comparison with other pooling algorithm 

If we choose to use mini-pooling to carry out VL pooling, we need `1.19` assays per individual on average. 

```R
> foo_mp = pooling_mc(pvl, riskscore, perm_num = 100, method = "minipool")
> mean(foo_mp)
[1] 5.96

> mean(foo_mp)/K
[1] 1.19
> 
```

If we choose to use MPA (May et al, 2010) to carry out VL pooling, we need `0.79` assays per individual on average. 

```R
> foo_mpa = pooling_mc(pvl, riskscore, perm_num = 100, method = "mpa")
> mean(foo_mpa)
[1] 3.94
> mean(foo_mpa)/K
[1] 0.79
```

Graphically  

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
