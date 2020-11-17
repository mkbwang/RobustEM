#' Robust EM Algorithm
#'
#' Calculate the mean and covariance matrix of each cluster based on the data matrix.
#'
#' @param datamat
#'
#' The matrix of points to be clustered.
#'
#' @param cluster
#'
#' Number of clusters
#'
#'
#' @param lambda
#'
#' Regularization Parameter
#'
#' @return
#'A list of results, include:
#'
#' Updated mean vectors for each cluster
#'
#' Updated covariances for each cluster
#'
#' Probability of each point in each cluster(a cluster*n matrix)
#'
#' Number of clusters (cluster)
#'
#' Dimension (d)
#'
#' Number of points (n)
#'
#' @export
#'
#' @examples
#' sim_info <- simMultGauss(n = 120, d = 2, cluster = 6, out_perc = 0.03, out_mag = 4)
#' result <- robustEM(sim_info[["simdata"]], cluster = 6)
#'
robustEM <- function(datamat, cluster, lambda = 3){
  # first do hierarchical clustering
  initial_info <- initial_hier(datamat, cluster)

  # then do EM algorithm
  result <- EM_alg_GMM(datamat, cluster, lambda = lambda,  inits = initial_info)

  class(result) <- "robustEM"

  return(result)
}
