
#' Hierarchical Clustering
#'
#' Carry out hierarchical clustering on the raw data.
#' Use the outcome as the initial value for EM algorithm. For complete example please check
#' out the documentation for function robustEM.
#'
#' @param sampleMat
#' A matrix with number of rows equal to the number of points and
#' the number of columns equal to the dimension.
#'
#' @param cluster
#' Prespecified number of clusters
#'
#' @return
#' A list with:
#'
#' initial means (matrix with dimension c*d)
#'
#' initial covariance (a list with c matrices of dimension d*d)
#'
#' initial proportion of clusters tau (a vector with c elements)
#'
#'
#' @export
#' @importFrom mclust hc hclass
#' @importFrom dplyr mutate group_by summarise_all funs select %>%
#' @importFrom stats cov
#'
initial_hier = function(sampleMat, cluster) {
  hclass1 = hclass(hc(sampleMat), G=cluster)
  mat_withlabels <- as.data.frame(sampleMat) %>% mutate(class = hclass1)
  initial_means = mat_withlabels %>% group_by(class) %>%
    summarise_all(funs(mean)) %>% select(-class) %>% as.matrix()
  initial_cov = lapply(1:cluster, function(x) {sampleMat[which(hclass1 == x), ] %>% cov()})
  initial_tau = sapply(1:cluster, function(x) {sum(hclass1 == x) / length(hclass1)})
  return(list(initial_means, initial_cov, initial_tau))
}
