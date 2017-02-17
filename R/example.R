# Pool size K is set to 5
K=5; n = 300;
# so, the number of pools = 60
n.pool  = n/K; n.pool
## number of resamples
n.resample = 100
set.seed(100)
pvl = rgamma(n, shape = 2.8, scale = 150)
summary(pvl)
riskscore = (rank(pvl)/n) * 0.5 + runif(n) * 0.5
cor(pvl, riskscore); cor(pvl, riskscore, method = "spearman")

#foo = cmpa(pvl, riskscore, K, perm_num= n.resample, vf_cut = 1000)
#mean(foo); mean(foo)/K
