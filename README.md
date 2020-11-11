
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RobustEM

<!-- badges: start -->

<!-- badges: end -->

The goal of RobustEM is to carry out clustering on high dimensional
points using expectation-maximization algorithm which is robust against
outliers.

## Usage Example

``` r
library(RobustEM)
```

The example here simulate 720 points belonging to 6 clusters. Each
cluster has 120 points. All the points have two dimensions. All the
clusters have approximately 6% of points as outliers. The outliers are
generated here as having the same mean but a covariance matrix with much
larger elements. This customized function `simMultGauss` not only
generates the data points but also includes detailed cluster information
(mean and covariance of each
cluster).

``` r
sim_info <- simMultGauss(n = 120, d = 2, cluster = 6, out_perc = 0.03, out_mag = 4)
```

We can use function `robustEM` to cluster the points.

``` r
result <- robustEM(sim_info[["simdata"]], cluster = 6)
```

I have written a customized summary function to summarize the cluster
mean and the number of points in each cluster.

``` r
summary(result)
#> $`Cluster Point Count`
#> point_cluster
#>   1   2   3   4   5   6 
#> 180 101 121 138 102  78 
#> 
#> $`Cluster Mean`
#>             V1        V2
#> [1,]  6.568047 22.426874
#> [2,] 14.327687 22.157949
#> [3,]  5.391978 36.699415
#> [4,] 29.668266  8.368254
#> [5,] 35.256049 12.657038
#> [6,] 14.900916 26.573123
```
