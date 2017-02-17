#' Number of Assays Need for Marker-Assisted Mini-Pooling with Algorithm
#'
#' This function uses Monte Carlo to compute the average number of assays
#' needed per pool using mMPA.
#'
#' alksdf detials
#'
#' @inheritParams mmpa0
#' @param v Vector of numerical assay results.
#' @param s Vector of risk score with the same length of viral load.
#' @param K Pool size, default is \code{K=5}.
#' @param vf_cut Cutoff for defining disease positive, default is \code{vf_cut = 1000}.
#' @param lod Vector of lower limited of detection, default is \code{lod = 0}.
#' @param perm_num The number of permutation to be used for the calculation, default is 100.
#' @param msg Message generated during calculation.
#' @return
#' The average number of assays needed per pool and per subject.
#' @keywords Pooling.
#' @export
#' @examples
#' d = Simdata
#' V = d$VL # Viral Load
#' S = d$S # Risk Score
#' mmpa(V, S, K = 3, perm_num = 3)
#' foo; table(foo)


mmpa = function(v, # vector of true VL
                s, # vector of risk score in same length
                K = 5, # pool size
                vf_cut = 500, # cutoff for individual viral failure
                lod = 0, # vector of true VL of those undetectable
                perm_num = 100,
                msg = T
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
  impafoo = mmpa0(v0, s0, K, vf_cut, lod)
  if(msg) cat("For pool size of", K, "the average number of assays needed per pool is",
              mean(impafoo), ", \ni.e. average number of assays per subject is",
              mean(impafoo)/K, ".")
  matrix(impafoo, ncol = perm_num)
}





#' Number of Assays Need for Marker-Assisted Mini-Pooling with Algorithm
#'
#' This function uses Monte Carlo to compute the average number of assays
#' needed per pool using mMPA.
#'
#' alksdf detials
#'
#' @param v Vector of numerical assay results.
#' @param s Vector of risk score with the same length of viral load.
#' @param K Pool size, default is \code{K=5}.
#' @param vf_cut Cutoff for defining disease positive, default is \code{vf_cut = 1000}.
#' @param lod Vector of lower limited of detection, default is \code{lod = 0}.
#' @return
#' The average number of assays needed per pool and per subject.
#' @keywords Pooling.
#' @export
#' @examples
#' d = Simdata
#' V = d$VL # Viral Load
#' S = d$S # Risk Score
#' mmpa(V, S, K = 3, perm_num = 3)
#' foo; table(foo)

mmpa0 = function(
  v, # vector of VL
  s, # vector of risk score in same length
  K = 5, # pool size
  vf_cut = 500, # cutoff for individual viral failure
  lod = 0 # vector of true VL of those undetectable
){
  n = length(v)
  if (length(s) != n) {
    warning("V and S have different length!")
    n = min(n, length(s))
  }
  num_pool = n %/% K
  foo0 = n - num_pool*K
  if(foo0 !=0) {
    #warning("Last (", foo0, ") observation(s) are dropped to make (", num_pool, ") pools!")
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

