
#' @export
summary.robustEM <- function(object, ...){
  t_mat <- object$T_mat
  point_cluster <- apply(t_mat, 2, which.max)

  output <- list("Cluster Point Count" = table(point_cluster),
                 "Cluster Mean" = object$mu)

  return(output)
}
