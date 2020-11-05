
#' Simulate multivariate gaussian distribution with some outliers.
#'
#' @param n
#' Number of points
#'
#' @param d
#' Number of dimensions
#'
#' @param out_perc
#' The proportion of outliers
#'
#' @param out_mag
#' The magnitude of the outliers in terms of the covariance
#'
#' @param cov_scale
#' Constant parameter for setting the covariance of the non-outliers. Default to be 1.
#'
#' @return
#' A list with:
#'
#' mu - mean vector
#'
#' sigma - covariance matrix
#'
#' gauss - matrix of points
#'
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats runif
#'
#' @examples
#' simulation <- multivarGaussian(n = 200, d = 3, out_perc = 0.03, out_mag = 4)
#'
multivarGaussian = function(n, d, out_perc, out_mag, cov_scale = 1){
  mu = runif(d,1,20*(1/d))

  sigma = matrix(runif(d*d, -1, 1), d, d)
  sigma = (sigma + t(sigma))/2
  eigs = eigen(sigma)$values
  if (min(eigs) <= 0.5) {sigma = sigma - (min(eigs) - 0.5) * diag(d)}

  sigma = cov_scale * sigma * runif(1,1,10)
  sigma_out = cov_scale * sigma * out_mag
  if (out_perc == 0 ) {
    gauss = mvrnorm(n, mu = mu, Sigma = sigma) #, tol = 1)
  }
  else {
    n_inlier <- round(n*(1-out_perc))
    n_outlier <- n - n_inlier

    # Else, pull from an MVN with the specified covariance matrix
    gauss1 = mvrnorm(n_inlier, mu = mu, Sigma = sigma) #, tol = 1)

    # Then pull the outliers from an MVN with the scaled up covariance matrix
    gauss2 = mvrnorm(n_outlier, mu = mu, Sigma = sigma_out) #, tol= 1)

    # Combine for one dataset
    gauss = rbind(gauss1, gauss2)
  }

  mvGauss = list(mu = mu, sigma = sigma, gauss = gauss)
  return(mvGauss)
}


#' Simulate data points from several different gaussian distributions. The number
#' of points from different distributions are the same.
#'
#' @param n
#' Number of points in each Gaussian distribution
#'
#' @param d
#' Dimension of each point
#'
#' @param c
#' Number of clusters
#'
#' @param out_perc
#' Proportion of outliers in each cluster
#'
#'
#' @param out_mag
#' Magitude of covariance difference between outliers and non-outliers
#'
#' @param cov_scale
#' Covariance Scaling constant for the covariance of non-outliers
#'
#' @return
#' A list with:
#'
#' All the means of each cluster in a c*d matrix
#'
#' All the covariances of each cluster in a list of matrices(length c, each matrix d*d)
#'
#' All the simulated data points(n*c rows, d columns)
#'
#' @export
#'
#'
#' @examples
#' sim_info <- simMultGauss(n = 120, d = 2, c = 6, out_perc = 0.03, out_mag = 4)
#'
#'
simMultGauss = function(n, d, c, out_perc, out_mag, cov_scale = 1){
  samples_simMultGauss = replicate(c, multivarGaussian(n = n, d = d,
                                                       out_perc = out_perc, out_mag = out_mag, cov_scale))

  sampleMu = do.call(rbind, samples_simMultGauss[1,])
  sampleSigma = lapply(samples_simMultGauss[2,], function(y) as.matrix(y))
  simSamp = do.call(rbind, samples_simMultGauss[3,])

  return(list(mus = sampleMu, sigmas = sampleSigma, simdata = simSamp))
}
