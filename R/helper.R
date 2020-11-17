
#' @export
summary.robustEM <- function(object, ...){
  t_mat <- object$T_mat
  point_cluster <- apply(t_mat, 2, which.max)

  output <- list("Cluster Point Count" = table(point_cluster),
                 "Cluster Mean" = object$mu)

  return(output)
}

#' @export
#' @import ggplot2
#' @importFrom mixtools ellipse
#' @importFrom dplyr %>%
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics plot.new
#'
plot.robustEM <- function(x, ...){
  raw_data <- x$raw # raw data

  Axis1 <- NULL
  Axis2 <- NULL

  meanvec <- as.data.frame(x$mu) # mean vector
  colnames(meanvec) <- c("Axis1", "Axis2")

  covars <- x$sigma # covariance matrices

  cluster_result <- data.frame(raw_data)
  colnames(cluster_result) <- c("Axis1", "Axis2")

  plot.new()
  ellipses <- lapply(1:nrow(meanvec), function(j) ellipse(mu=x$mu[j,],
                                                       sigma=covars[[j]], npoints = 300, newplot=F))

  covar_plot_f_EM <- function(j){
    boundary <- data.frame(ellipses[[j]])
    colnames(boundary) <- c("Axis1", "Axis2")
    return(geom_path(data = boundary, aes(x=Axis1, y=Axis2), size=0.5, color = 'red'))
  }

  covar_EM_plots <- lapply(1:nrow(x$mu), covar_plot_f_EM)

  final_plot <- ggplot() +
    geom_point(data = cluster_result, aes(x=Axis1, y=Axis2), size = 0.8, alpha = 0.5) +
    geom_point(data = meanvec, aes(x=Axis1, y=Axis2), shape=23, size = 3, stroke = 2, color = 'red') +
    covar_EM_plots

  return(final_plot)
}
