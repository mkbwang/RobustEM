
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
  Cluster_number <- NULL

  meanvec <- as.data.frame(x$mu) # mean vector
  colnames(meanvec) <- c("Axis1", "Axis2")
  meanvec$Cluster_number <- seq(1, nrow(meanvec))
  meanvec$Cluster_number <- as.factor(meanvec$Cluster_number)

  covars <- x$sigma # covariance matrices
  t_mat <- x$T_mat # probability of each point in each cluster

  hard_assign <- apply(t_mat, 2, which.max) %>% as.vector()

  cluster_result <- data.frame(cbind(raw_data, hard_assign))
  colnames(cluster_result) <- c("Axis1", "Axis2", "Cluster_number")
  cluster_result$Cluster_number <- as.factor(cluster_result$Cluster_number)

  plot.new()
  ellipses <- lapply(1:nrow(meanvec), function(j) ellipse(mu=x$mu[j,],
                                                       sigma=covars[[j]], npoints = 300, newplot=F))

  boundaries <- do.call(rbind, ellipses) %>% as.data.frame()
  colnames(boundaries) <- c("Axis1", "Axis2")
  boundaries$Cluster_number <- rep(seq(1,nrow(meanvec)), each = 300)
  boundaries$Cluster_number <- as.factor(boundaries$Cluster_number)

  final_plot <- ggplot() +
    geom_point(data = cluster_result, aes(x=Axis1, y=Axis2, color = Cluster_number), size = 0.8, alpha = 0.5) +
    geom_point(data = meanvec, aes(x=Axis1, y=Axis2, color = Cluster_number), shape=23, size = 3, stroke = 2) +
    geom_path(data = boundaries, aes(x=Axis1, y=Axis2, color = Cluster_number), size = 0.5)+
    scale_colour_brewer(palette="Dark2")

  return(final_plot)
}
