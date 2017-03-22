#' Monte Carlo Simulation for Estimating the Number of Assays Required
#' when Using Pooled Testing
#'
#' This function uses Monte Carlo (MC) to simulate different orders in
#' which the samples would be collected to form pools. Unlike the
#' function \code{minipool}, \code{mpa}, and \code{mmpa} that calculate
#' the number of assays
#' needed for pools that are formed following the exact order
#' of the samples that are listed in the data, this function
#' \code{pooling_mc} permutes the data many (\code{perm_num}) times
#' so as to estimate the average number of
#' assays required (ATR) per individual. Using MC avoids the dependence
#' on any
#' specific ordering of forming pools.
#'
#' @inheritParams mmpa
#' @param s A vector of risk scores; \code{s} must be provided if \code{method = "mmpa"}
#' and have the same length as \code{v}. The risk score \code{s} needs to be positively
#' associated with \code{v}; othewise an error message will be generated.
#' @param method Method that is used for pooled testing; must be one of \code{minipool},
#' \code{mpa}, and \code{mmpa}. By default, \code{method = "mmpa"}.
#' @param perm_num The number of permutation to be used for the calculation;
#' default is \code{100}.
#' @param msg Message generated during calculation; default is \code{FALSE}.
#' @return
#' The outcome is a matrix of dimension \code{num_pool} by \code{perm_num}.
#' The row number is the number of pools (\code{num_pool}) from each permutation
#' of the data, which is
#' determined by the sample size \code{N} and pool size \code{K}; \code{num_pool
#' = N\%/\%K}. The column number is the number of
#' permutations (\code{num_pool}).
#' @keywords Pooling.
#' @references
#' Bilder CR, Tebbs JM, Chen P. Informative retesting. Journal of the American
#' Statistical Association. 2010;105(491):942-955.
#'
#' May S, Gamst A, Haubrich R, Benson C, Smith DM. Pooled nucleic acid testing to
#' identify antiretroviral treatment failure during HIV infection. Journal of Acquired
#' Immune Deficiency Syndromes. 2010;53(2):194-201.
#'
#' Dorfman R. The detection of defective members of large populations. The
#' Annals of Mathematical Statistics. 1943;14(4):436-440.
#'
#' Our manuscript; to be added.
#' @export
#' @seealso \link{minipool}, \link{mpa}, \link{mmpa}
#' @examples
#' ### sample size = 300
#' n = 300;
#' set.seed(100)
#' pvl = rgamma(n, shape = 2.8, scale = 150)
#' summary(pvl)
#' #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#' #   53      225      392     424    564    1373
#' riskscore = (rank(pvl)/n) * 0.5 + runif(n) * 0.5
#' cor(pvl, riskscore, method = "spearman")
#' #   [1] 0.69
#' ### Pool size K is set to 5
#' K=5;
#' ### so, the number of pools = 60
#' n.pool  = n/K; n.pool
#' #   [1] 60
#' foo = pooling_mc(pvl, riskscore, perm_num = 100)
#' ### Average number of assays needed per pool for each of the 100
#' ### permutations of the data
#' apply(foo, 2, mean)
#' #   [1] 3.43 3.33 3.35 3.47 3.37 3.33 3.37 3.27 3.43 3.28 3.32 3.35 3.35 3.37
#' #   [15] 3.38 3.37 3.30 3.43 3.28 3.38 3.42 3.35 3.35 3.48 3.30 3.47 3.40 3.35
#' #  [29] 3.25 3.30 3.38 3.43 3.25 3.45 3.35 3.33 3.42 3.38 3.40 3.33 3.32 3.38
#' #   [43] 3.33 3.37 3.37 3.33 3.35 3.38 3.38 3.30 3.30 3.33 3.37 3.32 3.30 3.40
#' #   [57] 3.37 3.42 3.30 3.37 3.38 3.32 3.45 3.38 3.37 3.50 3.33 3.40 3.28 3.37
#' #   [71] 3.23 3.33 3.23 3.42 3.32 3.32 3.45 3.35 3.32 3.32 3.33 3.33 3.30 3.38
#' #   [85] 3.37 3.33 3.33 3.20 3.37 3.33 3.30 3.40 3.40 3.32 3.33 3.37 3.40 3.38
#' #   [99] 3.30 3.33
#' ### Estimated average number of assays needed per pool
#' mean(foo)
#' #   3.35
#' ### Estimated average number of assays needed per individual
#' mean(foo)/K
#' #   [1] 0.67
pooling_mc = function(v, # vector of true VL
                s = NULL, # vector of risk score in same length
                K = 5, # pool size
                vf_cut = 1000, # cutoff for individual viral failure
                lod = 0, # vector of true VL of those undetectable
                method = "mmpa",
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
  if(method == "mmpa"){
    impafoo = mmpa(v0, s0, K, vf_cut, lod, msg = msg)
  } else if(method == "mpa"){
    impafoo = mpa(v0, K, vf_cut, lod, msg = msg)
  } else if(method == "minipool") {
    impafoo = minipool(v0, K, vf_cut, lod, msg = msg)
  } else {
    stop("The function argument method must be one of minipool, mpa, and mmpa.")
  }

  if(msg) cat("For pool size of", K, "the average number of assays needed per pool is",
              mean(impafoo), ", \ni.e. average number of assays per subject is",
              mean(impafoo)/K, ".")
  matrix(impafoo, ncol = perm_num)
}






#' Number of Assays Required using Marker-Assisted Mini-Pooling with Algorithm
#' (mMPA)
#'
#' Function \code{mmpa(...)} calculates the number of assays required, when
#' using mMPA, for
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
#' @references
#'
#' Bilder CR, Tebbs JM, Chen P. Informative retesting. Journal of the American
#' Statistical Association. 2010;105(491):942-955.
#'
#' May S, Gamst A, Haubrich R, Benson C, Smith DM. Pooled nucleic acid testing
#' to identify antiretroviral treatment failure during HIV infection. Journal
#' of Acquired Immune Deficiency Syndromes. 2010;53(2):194-201.
#'
#' Our manuscript (under review);  to be added.
#'
#' @inheritParams mpa
#' @param s A vector of risk scores; \code{s} must have the same
#' length as \code{v}. The risk score \code{s} needs to be positively
#' associated with \code{v}; othewise an error message will be generated.
#'
#' @return
#' A vectorof length \code{N\%/\%K} for the numbers of assays needed for all pools
#' that are formed.
#' @keywords mMPA.
#' @seealso \link{minipool}, \link{mpa}, \link{pooling_mc}
#' @export
#' @examples
#' K=5; n = 50;
#' n.pool  = n/K; n.pool
#' #   [1] 10
#' set.seed(100)
#' pvl = rgamma(n, shape = 2.8, scale = 150)
#' riskscore = (rank(pvl)/n) * 0.5 + runif(n) * 0.5
#' mmpa(v = pvl, s = riskscore)
#' #   A total of 10 pools are formed.
#' #   The numbers of assays required by these pools are:
#' #   [1] 3 3 4 4 2 3 3 4 3 3
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


#' Number of Assays Needed using Mini-Pooling with Algorithm (MPA)
#'
#' Function \code{mpa(...)} calculates the number of assays required, when using MPA, for
#' pools that are formed following the order of individual samples in the data.
#'
#' For a given sample v_i, i = 1, ..., N, the first \code{K} samples v_1,
#' ..., v_5 are combined to form a pool, the next \code{K} samples v_6, ...,
#' v_10 are combined to form the second
#' pool, and so on. If the number of samples for the last pool is less than
#' \code{K}, these remaining samples are not used to form a pool (i.e.
#' not included in the calculation) . Therefore, a total of
#' \code{N\%/\%K} pools are formed. The function calculates the number of
#' assays needed for each of these pools. See May et al (2010).
#'
#' @references May, S., Gamst, A., Haubrich, R., Benson, C., & Smith,
#' D. M. (2010). Pooled nucleic acid testing to identify antiretroviral
#' treatment failure during HIV infection. Journal of acquired immune
#' deficiency syndromes (1999), 53(2), 194.
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
#' @seealso \link{minipool}, \link{mmpa}, \link{pooling_mc}
#' @export
#' @examples
#' K=5; n = 50;
#' n.pool  = n/K; n.pool
#' #   [1] 10
#' set.seed(100)
#' pvl = rgamma(n, shape = 2.8, scale = 150)
#' mpa(v = pvl)
#' #  A total of 10 pools are formed.
#' #   The numbers of assays required by these pools are:
#' #   [1] 3 3 4 4 2 5 4 4 4 4
#'
mpa <- function(
  v,
  K = 5,
  vf_cut = 1000,
  lod = 0,
  msg = T
){
  s = c(1:length(v))
  mmpa(v, s, K, vf_cut, lod, msg)
}





#' Number of Assays Needed using Mini-Pooling
#'
#' Function \code{minipool(...)} calculates the number of assays required, when
#' using mini-pooling, for
#' pools that are formed following the order that individual samples appear in the data.
#'
#' Suppose that N samples are collected for pooled testing. The first
#' \code{K} samples are combined to
#' form a pool, the next \code{K} samples are combined to form the second
#' pool, and so on. If the number of samples for the last pool is less than
#' \code{K}, these remaining samples are not used to form a pool (i.e.
#' not included in the calculation). Therefore, a total of
#' \code{N\%/\%K} pools are formed. The function calculates the number of
#' assays needed for each of these pools. For mini-pooling, if a pool is
#' negative, no further tests are needed and all samples in the pool
#' are concluded as being negative; so the total number of
#' assays required is one. Otherwise if the pool is tested positive, all
#' individual samples in the pool are tested and the total number of assays
#' required is \code{(K + 1)}.
#'
#' @inheritParams mpa
#'
#' @return
#' A vectorof length \code{N\%/\%K} for the numbers of assays needed for all pools
#' that are formed.
#' @keywords mini-pooling.
#' @references
#' Dorfman R. The detection of defective members of large populations. The
#' Annals of Mathematical Statistics. 1943;14(4):436-440.
#' @seealso \link{mpa}, \link{mmpa}, \link{pooling_mc}
#' @export
#' @examples
#' K=5; n = 50;
#' n.pool  = n/K; n.pool
#' #  [1] 10
#' set.seed(100)
#' pvl = rgamma(n, shape = 2.8, scale = 150)
#' minipool(pvl)
#' #  A total of 10 pools are formed.
#' #  The numbers of assays required by these pools are:
#' #  [1] 6 6 6 6 6 6 6 6 6 6
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
