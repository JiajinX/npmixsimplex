# Sort components of Dirichlet mixture
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
  h <- 1 / (a0 - D)
  mix <- list("alpha" = a, "pr" = pr, "h" = h)
  attr(mix, "mode") <- (a - 1) / (a0 - D)
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
  s <- 1
  for (i in 1:m) {
    cnt <- sum(idx == i)
    e <- s + cnt - 1
    data[s:e, ] <- gtools::rdirichlet(cnt, mix$alpha[i, ])
    s <- e + 1
  }
  data[sample(n), ]
}

# Calculate matrix log D(x_i; α_j) with n * m
outer.dirichlet <- function(x, a, log = FALSE) {
  if (is.vector(x)) x <- t(x)
  if (is.vector(a)) a <- t(a)
  n <- nrow(x)
  m <- nrow(a)
  d <- compositions::dDirichlet(x[rep(1:n, m), , drop = FALSE],
                                a[rep(1:m, rep(n, m)), , drop = FALSE],
                                log = log) # log f(x_i; α_j)
  dim(d) <- c(n, m) # row: x_i; col: α_j; dim = n * m
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
  logd <- outer.dirichlet(x, mix$alpha, log = TRUE) +
    rep(log(mix$pr), rep(nrow(x), length(mix$pr)))
  # log D(x_i;α_j) + log π_j = log [π_j * D(x_i;α_j)]
  ma <- lsei::matMaxs(logd)
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
