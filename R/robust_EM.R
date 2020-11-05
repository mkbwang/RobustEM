

#' Robust EM Algorithm for Gaussian Mixture Models. Given the prespecified number of clusters,
#' the algorithm will calculate the mean and covariance of each cluster.
#'
#' @param sampleMat
#' All the points in a matrix
#'
#' @param c
#' Number of clusters
#'
#' @param lambda
#' Regularization Parameter
#'
#' @param inits
#' List of initial values, include:
#'
#' mean (a c*d matrix)
#'
#' Covariance (a list with c matrices of dimension d*d)
#'
#' proportion of clusters (a vector with c elements adding up to one)
#'
#'
#' @return
#' A list of results, include:
#'
#' Updated mean vectors for each cluster
#'
#' Updated covariances for each cluster
#'
#' Probability of each point in each cluster(a c*n matrix)
#'
#' Number of clusters (c)
#'
#' Dimension (d)
#'
#' Number of points (n)
#'
#'
#' @export
#'
#' @examples
#' sim_info <- simMultGauss(n = 120, d = 2, c = 6, out_perc = 0.03, out_mag = 4)
#' initial_info <- initial_hier(sim_info[["simdata"]], c=6)
#' result <- EM_alg_GMM(sim_info[['simdata']], c=6, inits = initial_info)
#'
EM_alg_GMM = function(sampleMat, c, lambda = 10, inits) {

  ## Observations
  x = sampleMat
  n = nrow(x) # number of observations
  d = ncol(x) # number of dimensions

  tau = matrix(inits[[3]], c, 1) # initial cluster proportion
  mu = inits[[1]] # initial mean
  sigma = inits[[2]] # initial covariance

  T_mat = matrix(0,c,n) # probability of each point in each cluster
  max_it = 200

  for (l in unique(c(Inf, lambda^(5:1)))) { #unique
    old_mu = mu

    for (it in 1:max_it) {
      old_T_mat = T_mat

      # E STEP: Construct the vector of the denom for T_mat for each i
      for (j in 1:c) {
        #print(sigma[[j]])
        if (det(sigma[[j]]) < 1e-7) {
          sigma[[j]] = diag(d) * 1e-1
        }
        L = chol(sigma[[j]], pivot = TRUE)
        y = solve(t(L), t(x - rep(mu[j,],each=n)))
        T_mat[j,] = log(tau[j]) - 0.5 * colSums(y^2) - log(sqrt(2*pi*abs(prod(diag(L))))) # last term needs to be cleared up
      }

      for (i in 1:n) {
        scale = max(T_mat[,i])
        normalizer <- sum(exp(T_mat[,i] - scale))
        T_mat[,i] = exp(T_mat[,i] - scale) / normalizer
      }

      # M STEP: Update tau, mu and sigma

      # Update tau
      for (j in 1:c) {
        e = matrix(0, n, d)
        tau[j,] = (1/n)*sum(T_mat[j,])

        x1 = x - rep(mu[j, ], each = n)
        sigma_inv = solve(as.matrix(sigma[[j]]))

        x1_sigma_inv = x1 %*% sigma_inv
        A = matrix(rowSums(x1_sigma_inv * x1), n, 1)

        for (i in 1:n) {
          if (A[i,] < l) { e[i,] = 0 }
          else { e[i,] = x[i,] - mu[j,] }
        }

        # Update mu
        mu[j,] = colSums(T_mat[j,]%*%(x-e))/sum(T_mat[j,])

        # Update sigma
        num = matrix(0,d,d)

        # Only use those points with error 0 to calc to the covariance martix
        indices = which(e[,1]==0)
        denom = sum(T_mat[j,indices])
        x_prime = x[indices,]-e[indices,]-matrix(mu[j,], ncol=d, nrow=length(indices), byrow=T)
        print(dim(x_prime))
        print(dim(diag(T_mat[j,indices])))
        num = t(x_prime)%*%diag(T_mat[j,indices])%*%x_prime
        sigma[[j]] = matrix(num/denom, d, d)

        if (sum(T_mat[j,]!=0) <= 3) {
          mu[j,] = matrix(0, 1, d)
          sigma[[j]] = diag(d) * 1e-2
        }


      }
      if (max(abs(old_T_mat - T_mat)) < 1e-10) break
    }
    if(max(abs(old_mu - mu)) < 1e-6) break

  }

  returnList = list(mu, sigma, T_mat, tau, c, d, n)
  names(returnList) = c("mu", "sigma", "T_mat", "tau", "c", "d", "n")
  return(returnList)

}
