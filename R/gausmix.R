# Generate random data from Gaussian mixture
rmvnormix <- function(n, mix) {
  if (class(mix)[1] == "mvnmix") {
    npden::rmvnmix(n, mix)
  } else {
    ks::rmvnorm.mixt(n, mix$mu, mix$Sigma, props = mix$pr)
  }
}

# Density value of Gaussian mixture
dmvnormix <- function(x, mix, log = FALSE) {
  if (class(mix)[1] == "mvnmix") {
    d <- npden::dmvnmix(x, mix, log = log)
  } else {
    d <- ks::dmvnorm.mixt(x, mix$mu, mix$Sigma, props = mix$pr)
    if (log) d <- log(d)
  }
  d
}
