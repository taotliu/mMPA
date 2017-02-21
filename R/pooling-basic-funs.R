#' Monte Carlo Simulation for Estimating the Number of Assays Needed
#' when Using mMPA
#'
#' This function uses Monte Carlo to simulate different orders in
#' which the samples would be collected to form pools. Unlike the
#' function \code{mmpa0} which calculates the number of assays
#' needed for pools that are formed following the exact order
#' of the samples that are listed in the data, the function
#' \code{mmp} permutes the data many (\code{perm_num}) times
#' and thus can be used to estimate the average number of
#' assays required (ATR) per individual without depending on a
#' specific configuration of forming pools. To obtain a reliable
#' estimate of ATR, a large number of \code{perm_num} should be used.
#'
#' See \link{mmpa0} for detail descriptions
#'
#' @inheritParams mmpa0
#' @param perm_num The number of permutation to be used for the calculation;
#' default is \code{100}.
#' @param msg Message generated during calculation; default is \code{FALSE}.
#' @return
#' The outcome is a matrix of dimension \code{num_pool} by \code{perm_num}.
#' The row number is the number of pools (\code{num_pool}) for each permutation,
#' which is
#' determined by the sample size \code{N} and pool size \code{K}; \code{num_pool
#' = N\%/\%K}. The column number is the number of
#' permutations (\code{num_pool}).
#' @keywords Pooling.
#' @export
#' @seealso \link{mmpa0}
#' @examples
#' d = Simdata
#' V = d$VL # Viral Load
#' S = d$S # Risk Score
#' mmpa(V, S, K = 3, perm_num = 3)
#' foo; table(foo)


pooling_mc = function(v, # vector of true VL
                s, # vector of risk score in same length
                K = 5, # pool size
                vf_cut = 1000, # cutoff for individual viral failure
                lod = 0, # vector of true VL of those undetectable
                perm_num = 100,
                msg = F
){
  ##########################
  n = length(v)
  pool.n = (n%/%K)*K
  #### permutation function
  sample.foo = function(x){sample(x, pool.n, replace = F)}
  ####
  #permindex = sample(rep(1:n, perm_num), perm_num*pool.n, replace = T)
  permindex = apply(matrix(rep(1:n, perm_num), ncol = perm_num), 2, sample.foo)
  ###########################
  v0 = v[permindex]
  s0 = s[permindex]
  impafoo = mmpa0(v0, s0, K, vf_cut, lod, msg = msg)
  if(msg) cat("For pool size of", K, "the average number of assays needed per pool is",
              mean(impafoo), ", \ni.e. average number of assays per subject is",
              mean(impafoo)/K, ".")
  matrix(impafoo, ncol = perm_num)
}






#' Number of Assays Needed using Marker-Assisted Mini-Pooling with Algorithm
#'
#' Function \code{mmpa0} calculates the number of assays, when using mMPA, for
#' pools that are formed following the order of individual samples in the data.
#'
#' For a given sample (v_i, s_i), i = 1, ..., N, the first \code{K} samples are combined to
#' form a pool, the next \code{K} samples are combined to form the second
#' pool, and so on. If the number of samples for the last pool is less than
#' \code{K}, these remaining samples are not used to form a pool (i.e.
#' not included
#' in the calculation) . Therefore, a total of
#' \code{N\%/\%K} pools are formed. The function calculates the number of
#' assays needed for each of these pools.
#'
#' @inheritParams mpa0
#' @param s A vector of risk scores; \code{v} and \code{s} must have the same
#' length.
#'
#' @return
#' A vectorof length \code{N\%/\%K} for the numbers of assays needed for all pools
#' that are formed.
#' @keywords mMPA.
#' @seealso \link{mmpa}
#' @export
#' @examples
#' K=5; n = 50;
#' n.pool  = n/K; n.pool
#' > [1] 10
#' set.seed(100)
#' pvl = rgamma(n, shape = 2.8, scale = 150)
#' riskscore = (rank(pvl)/n) * 0.5 + runif(n) * 0.5
#' mmpa0(v = pvl, s = riskscore)
#' > A total of 10 pools are formed.
#' > The numbers of assays required by these pools are:
#' > [1] 3 3 4 4 2 3 3 4 3 3
#'
mmpa = function(
  v,
  s,
  K = 5,
  vf_cut = 1000,
  lod = 0,
  msg = T
){
  n = length(v)
  if (length(s) != n) {
    warning("v and s have different length!")
    n = min(n, length(s))
  }
  if (sum(c(is.na(v), is.na(s)))>0) {
    stop("v or s contains missing value!")
  }
  num_pool = n %/% K
  foo0 = n - num_pool*K
  if(msg) {
    cat("A total of", num_pool, "pools are formed.\n")
    cat("The numbers of assays required by these pools are: \n")
  }
  if(foo0 !=0) {
    if(msg)
      warning("The last ", foo0, " samples(s) are dropped to make the ", num_pool, " pools!")
    v = v[1:(num_pool*K)]
    s = s[1:(num_pool*K)]
  }
  ### create a matrix of size  K x num.pool.
  ### i.e. each column is a pool of size K
  #v_mat = matrix(v, nrow = K)
  s_mat = matrix(s, nrow = K)
  ### re-order VL based on the rank of S
  ### s.t. the 1st row has the lowest risk score
  ### and the Kth row the highest
  order_mat = apply(s_mat, 2, order)
  foo = c(order_mat) + rep(0:(num_pool-1), rep(K, num_pool)) * K
  v_mat = matrix(v[foo], nrow = K)
  #print(matrix(s_mat[foo], nrow=K))
  #print(v_mat)
  t_mat = apply(v_mat, 2, cumsum)
  #print(t_mat)
  t0_mat = 0
  if(lod > 0){
    v0_mat = v_mat
    v0_mat[v0_mat >= lod] = 0
    t0_mat = apply(v0_mat, 2, function(x) sum(x)-cumsum(x))
  }
  return(1 + apply((t_mat+t0_mat)[-1,] > vf_cut, 2, sum))
}


#' Number of Assays Needed using Mini-Pooling with Algorithm
#'
#' Function \code{mpa0} calculates the number of assays, when using MPA, for
#' pools that are formed following the order of individual samples in the data.
#'
#' For a given sample v_i, i = 1, ..., N, the first \code{K} samples are combined to
#' form a pool, the next \code{K} samples are combined to form the second
#' pool, and so on. If the number of samples for the last pool is less than
#' \code{K}, these remaining samples are not used to form a pool (i.e.
#' not included
#' in the calculation) . Therefore, a total of
#' \code{N\%/\%K} pools are formed. The function calculates the number of
#' assays needed for each of these pools.
#'
#' @param v A vector of non-negative numerical assay results.
#' @param K Pool size; default is \code{K = 5}.
#' @param vf_cut Cutoff value for defining positive cases;
#' default is \code{vf_cut = 1000}.
#' @param lod A vector of lower limits of detection or a scalar if the limits are the
#' same; default is \code{lod = 0}.
#' @param msg Message generated during calculation; default is \code{TRUE}.
#' @return
#' A vectorof length \code{N\%/\%K} for the numbers of assays needed for all pools
#' that are formed.
#' @keywords mMPA.
#' @seealso \link{mmpa}
#' @export
#' @examples
#' K=5; n = 50;
#' n.pool  = n/K; n.pool
#' > [1] 10
#' set.seed(100)
#' pvl = rgamma(n, shape = 2.8, scale = 150)
#' mpa0(v = pvl)
#' > A total of 10 pools are formed.
#' > The numbers of assays required by these pools are:
#' > [1] 3 3 4 4 2 5 4 4 4 4
#'
mpa <- function(
  v,
  K = 5,
  vf_cut = 1000,
  lod = 0,
  msg = T
){
  s = c(1:length(v))
  mmpa0(v, s, K, vf_cut, lod, msg)
}





#' Number of Assays Needed using Mini-Pooling
#'
#' Function \code{minipool(...)} calculates the number of assays required, when
#' using mini-pooling, for
#' pools that are formed following the order that individual samples appear in the data.
#'
#' Suppose that N samples are collected for pooled testing, the first
#' \code{K} samples are combined to
#' form a pool, the next \code{K} samples are combined to form the second
#' pool, and so on. If the number of samples for the last pool is less than
#' \code{K}, these remaining samples are not used to form a pool (i.e.
#' not included in the calculation). Therefore, a total of
#' \code{N\%/\%K} pools are formed. The function calculates the number of
#' assays needed for each of these pools. For mini-pooling, if a pool is
#' negative, no further tests are needed and all samples in the pool
#' are concluded as negative; so the total number of
#' test required is one. Otherwise if the pool is tested positive, all
#' individual samples in the pool are tested and the total number of assays
#' required is (K + 1).
#'
#' @inheritParams mpa
#'
#' @return
#' A vectorof length \code{N\%/\%K} for the numbers of assays needed for all pools
#' that are formed.
#' @keywords mMPA.
#' @seealso \link{mpa}, \link{mmpa}, \link{pooling_mc}
#' @export
#' @examples
#' K=5; n = 50;
#' n.pool  = n/K; n.pool
#' > [1] 10
#' set.seed(100)
#' pvl = rgamma(n, shape = 2.8, scale = 150)
#' riskscore = (rank(pvl)/n) * 0.5 + runif(n) * 0.5
#' mmpa0(v = pvl, s = riskscore)
#' > A total of 10 pools are formed.
#' > The numbers of assays required by these pools are:
#' > [1] 6 6 6 6 6 6 6 6 6 6
#'
minipool <- function(
  v,
  K = 5,
  vf_cut = 1000,
  lod = 0,
  msg = T
){
  s = c(1:length(v))
  foo = mmpa(v, s, K, vf_cut, lod, msg)
  foo[foo>1] = K+1
  foo
}
