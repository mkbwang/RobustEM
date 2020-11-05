
#' Carry out hierarchical clustering on the raw data.
#' Use the outcome as the initial value for EM algorithm.
#'
#' @param sampleMat
#' A matrix with number of rows equal to the number of points and
#' the number of columns equal to the dimension.
#'
#' @param c
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
#' @examples
#' sim_info <- simMultGauss(n = 120, d = 2, c = 6, out_perc = 0.03, out_mag = 4)
#' initial_info <- initial_hier(sim_info[["simdata"]], c=6)
#'
initial_hier = function(sampleMat, c) {
  hclass1 = hclass(hc(sampleMat), G=c)
  mat_withlabels <- as.data.frame(sampleMat) %>% mutate(class = hclass1)
  initial_means = mat_withlabels %>% group_by(class) %>%
    summarise_all(funs(mean)) %>% select(-class) %>% as.matrix()
  initial_cov = lapply(1:c, function(x) {sampleMat[which(hclass1 == x), ] %>% cov()})
  initial_tau = sapply(1:c, function(x) {sum(hclass1 == x) / length(hclass1)})
  return(list(initial_means, initial_cov, initial_tau))
}
