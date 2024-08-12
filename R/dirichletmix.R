## ======= Basic function for Dirichlet, Beta (mixture) distribution ======= ##

# Sort components of Dirichlet mixture
# by: index of column used to sort
sort.dirichletmix <- function(mix, by = 1, decreasing = FALSE) {
  if (length(mix$pr) <= 1) return(mix)
  o <- order(mix$alpha[, by], decreasing = decreasing)
  mix$alpha <- mix$alpha[o, ]
  mix$pr <- mix$pr[o]
  attr(mix, "mode") <- attr(mix, "mode")[o, ]
  mix
}

# Create a class called "dirichletmix" to set parameters of Dirichlet mixture
dirichletmix <- function(a = rep(2, 3), pr = 1, sort = TRUE, by = 1) {
  if (is.vector(a)) a <- t(a)
  m <- nrow(a)
  D <- ncol(a)
  a0 <- sum(a[1, ])
  pr <- rep(pr, len = m)
  pr <- pr / sum(pr)
  h <- 1 / (a0 - D) # bandwidth
  at <- a - 1 # To avoid some numeric issue in R
  if (any(at < -1e-8)) stop("Alpha should larger than 1")
  at[at < 0] <- 0
  a <- at / rowSums(at) * (a0 - D) + 1
  mix <- list("alpha" = a, "pr" = pr, "h" = h)
  attr(mix, "mode") <- (a - 1) * h
  class(mix) <- "dirichletmix"
  if (sort) mix <- sort(mix, by = by)
  mix
}

# Generate random data from Dirichlet mixture distribution
rdirichletmix <- function(n, mix) {
  m <- nrow(mix$alpha)
  D <- ncol(mix$alpha)
  idx <- sample.int(m, n, prob = mix$pr, replace = TRUE)
  data <- matrix(nrow = n, ncol = D)
  start <- 1 # index
  for (i in seq_len(m)) {
    cnt <- sum(idx == i)
    end <- start + cnt - 1
    data[start:end, ] <- gtools::rdirichlet(cnt, mix$alpha[i, ])
    start <- end + 1
  }
  data[sample(n), ]
}

# outer = FALSE: calculate vector log D(x_i; α_i) with length n
# outer = TRUE: calculate matrix log D(x_i; α_j) with n * m
ddirichlet <- function(x, a, log = FALSE, outer = FALSE) {
  if (is.vector(x)) x <- t(x)
  if (is.vector(a)) a <- t(a)
  n <- nrow(x)
  m <- nrow(a)
  logx <- log(x)
  logmultibeta <- lgamma(rowSums(a)) - rowSums(lgamma(a)) # Constant part
  if (!outer) { # Got a vector by adopting cycling rule on x, a
    a <- a[rep_len(1:m, n), ]
    logmultibeta <- rep_len(logmultibeta, n)
  } else { # Got a matrix like outer function
    a <- a[rep(1:m, each = n), ]
    logmultibeta <- rep(logmultibeta, each = n)
    logx <- logx[rep(1:n, m), ]
  }
  xa <- (a - 1) * logx
  xa[is.na(xa)] <- 0
  d <- logmultibeta + rowSums(xa)
  if (outer) dim(d) <- c(n, m)
  if (!log) d <- exp(d)
  d
}

# Calculate c(log sum(π_j * D(x_1; α_j)), ..., log sum(π_j * D(x_n; α_j)))
# Length of vector = the number of data
ddirichletmix <- function(x, mix, log = FALSE) {
  if (is.vector(x)) x <- t(x)
  j <- mix$pr == 0
  if (any(j)) {
    mix$alpha <- mix$alpha[!j, , drop = FALSE]
    mix$pr <- mix$pr[!j]
  }
  logd <- ddirichlet(x, mix$alpha, log = TRUE, outer = TRUE) +
    rep(log(mix$pr), rep(nrow(x), length(mix$pr)))
  # log D(x_i;α_j) + log π_j = log [π_j * D(x_i;α_j)]
  ma <- lsei::matMaxs(logd)
  ma[ma == -Inf] <- 0 # Otherwise -Inf-(-Inf) = NaN ???
  pid <- exp(logd - ma)
  pis <- rowSums(pid)
  r <- ma + log(pis)
  # c(log sum(π_j * D(x_1; α_j)), ..., log sum(π_j * D(x_n; α_j)))
  if (!log) r <- exp(r)
  r
}

# Log-likelihood of Dirichlet mixture
logLik.dirichletmix <- function(mix, x) {
  sum(ddirichletmix(x, mix, log = TRUE))
}

# Create a class called "betamix" to set parameters of beta mixture
betamix <- function(shape = c(1, 1), pr = 1) {
  if (is.vector(shape)) shape <- t(shape)
  if (ncol(shape) != 2) stop("wrong shape parameters")
  m <- nrow(shape)
  pr <- rep(pr, len = m)
  pr <- pr / sum(pr)
  sumshape <- sum(shape[1, ])
  merge <- tapply(pr, shape[, 1], sum) # Merge components with the same shape
  a <- as.numeric(names(merge))
  shape <- cbind(a, sumshape - a)
  colnames(shape) <- names(merge) <- NULL
  mix <- list("shape" = shape, "pr" = merge)
  attr(mix, "mode") <- (a - 1) / (sumshape - 2)
  mix
}

# Calculate density of beta mixture
dbetamix <- function(x, mix) {
  j <- mix$pr == 0
  if (any(j)) {
    mix$shape <- mix$shape[!j, , drop = FALSE]
    mix$pr <- mix$pr[!j]
  }
  n <- length(x)
  m <- nrow(mix$shape)
  shape <- mix$shape[rep(1:m, rep(n, m)), ]
  outerd <- dbeta(x[rep(1:n, m)], shape[, 1], shape[, 2], log = TRUE)
  dim(outerd) <- c(n, m)
  logd <- outerd + rep(log(mix$pr), rep(n, m))
  rowSums(exp(logd))
}
