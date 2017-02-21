n = 300;
set.seed(100)
pvl = rgamma(n, shape = 2.8, scale = 150)
summary(pvl)

riskscore = (rank(pvl)/n) * 0.5 + runif(n) * 0.5
cor(pvl, riskscore, method = "spearman")

# Pool size K is set to 5
K=5;
# so, the number of pools = 60
n.pool  = n/K; n.pool

foo = pooling_mc(pvl, riskscore, perm_num = 100)
mean(foo)
mean(foo)/K


foo_mp = pooling_mc(pvl, riskscore, perm_num = 100, method = "minipool")
mean(foo_mp)
mean(foo_mp)/K

foo_mpa = pooling_mc(pvl, riskscore, perm_num = 100, method = "mpa")
mean(foo_mpa)
mean(foo_mpa)/K


boxplot(cbind(MP=apply(foo_mp, 2, mean),
              MPA=apply(foo_mpa, 2, mean),
              mMPA=apply(foo, 2, mean))/K*100,
        ylab = "Average number of assays required per 100 individuals")


