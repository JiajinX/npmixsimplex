# Find smoothing factor λ by max pseudo-log-likelihood in
# KDE of Aitchison method on alr transformation
# data: data on simplex
Halr <- function(data) {
  adata <- unclass(compositions::alr(data)) # simplex -> alr space
  jac <- exp(rowSums(log(data)))
  cm <- cov(adata) # Covariance matrix on alr space
  r <- optimise(function(lamb, data, cm, jac) {
    dm <- npden:::outer.dmvn(data, data, Sigma = lamb * cm)
    dm <- t(t(dm) / jac)
    diag(dm) <- 0
    sum(log(rowSums(dm)))
  }, interval = c(1e-3, 1), data = adata, cm = cm, jac = jac,
  maximum = TRUE, tol = 1e-4)
  H <- r$maximum * cm
  attr(H, "smoothing_factor") <- r$maximum
  H
}

# Find smoothing factor λ by max pseudo-log-likelihood in KDE of Aitchison,
# using Dirchlet kernel on simplex space.
hkdeDir <- function(data) {
  r <- optimise(function(lamb, data) {
    a <- data / lamb + 1
    dm <- ddirichlet(data, a, outer = TRUE)
    diag(dm) <- 0
    dtest <- rowSums(dm)
    sum(log(dtest[dtest > 0]))
    # Avoid zero density test on boundary but training set not
  }, interval = c(1e-3, 1), data = data, maximum = TRUE, tol = 1e-4)
  r$maximum
}
