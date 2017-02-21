# mMPA
R package for implementing marker-assisted mini-pooling with algorithm (mMPA)

## Installation 

The latest version of the `mMPA` package is available at GitHub [taotliu/mMPA](http://github.com/taotliu/mMPA). It requires the `devtools` package to install it in R. If you have not had `devtools` on your computer, using the following code  `install.packages("devtools")` to install the `devtools` package. Then run the following codes to install the `mMPA` package. 

```R
devtools::install_github("taotliu/mMPA")
library(mMPA)
```

## Example 

### Estimate the average number of assays required by mMPA 

The following R code example demonstrates the use of the `mMPA` package. 

Let us assume that blood samples of `n = 300` individuals are collected for HIV viral load (VL) testing. We simulate the results of VL tests using a Gamma (shape = 2.8, scale = 150) distribution, and generate the corresponding risk scores by adding a uniform random noise to the percentile of VL. The resulting VL has a median of `392` and an inter-quantile from `224` to `565`, and the resulting risk score has a Spearman’s correlation of `0.69` with the VL. 

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

We will do a pooled testing using a pool size of `K = 5`. So a total of 60 pools will be formed. 
```R
 > # Pool size K is set to 5
 > K = 5
 > # so, the number of pools = 60
 > n.pool  = n/K; n.pool
 [1] 60
``` 
Of course, there are many ways to form pools. We will permute the data `perm_num = 100` to mimic the situations that the individuals came to the clinic in different orders and different choices of five blood samples are pooled. 

The `mMPA` package includes a function called `mmpa(v, s, K, perm_num, vf_cut, ...)`, which takes four main arguments as input: Values of test results (`v`), corresponding risk scores (`s`), pool size (`K`), the number of Monte Carlo simulations (`perm_num`), and threshold for defining test positivity (`vf_cut`).  The function outputs the total number of VL assays needed for each of the 60 pools for each permutation. 

```R
 > foo = mmpa(pvl, riskscore, K, perm_num = 100, vf_cut = 1000)
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



