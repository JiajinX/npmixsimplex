library(gtools)
library(compositions)
library(viridis)
library(npden)

# Dirichlet mixture parameters
dirichletmix <- function(a = rep(2, 3), pr = 1, h) {
  if (is.vector(a)) a <- t(a)
  k <- nrow(a)
  pr <- rep(pr, len = k)
  pr <- pr / sum(pr)
  mix <- list("alpha" = a, "pr" = pr, "h" = h)
  class(mix) <- "dirichletmix"
  mix
}

# generate random data from Dirichlet mixture distribution
rdirichletmix <- function(n, mix, threshold2zero = NULL) {
  m <- nrow(mix$alpha)
  i <- sample.int(m, n, prob = mix$pr, replace = TRUE)
  data <- list()
  for (j in 1:m) { # number of components
    data[[j]] <- rdirichlet(sum(i == j), alpha = mix$alpha[j, ])
  }
  data <- do.call(rbind, data)
  data
}

# output: matrix log f(x_i; α_j) with n * m
outer.dirichlet <- function(x, a, log = FALSE) {
  if (is.vector(x)) x <- t(x)
  if (is.vector(a)) a <- t(a)
  n <- nrow(x)
  m <- nrow(a)
  d <- dDirichlet(x[rep(1:n, m), , drop = FALSE],
                  a[rep(1:m, rep(n, m)), , drop = FALSE],
                  log = log) # log f(x_i; α_j)
  dim(d) <- c(n, m) # row: x_i; col: α_j; dim = n * m
  d
}

# output: vector [log sum(π_j*f(x_1;α_j)), ..., log sum(π_j*f(x_n;α_j))]
ddirichletmix <- function(x, mix, log = FALSE) {
  if (is.vector(x)) x <- t(x)
  j <- mix$pr == 0
  if (any(j)) {
    mix$alpha <- mix$alpha[!j, , drop = FALSE]
    mix$pr <- mix$pr[!j]
  }
  logd <- outer.dirichlet(x, mix$alpha, log = TRUE) +
    rep(log(mix$pr), rep(nrow(x), length(mix$pr)))
  # log f(x_i;α_j) + log π_j = log [π_j * f(x_i;α_j)]
  ma <- matMaxs(logd)
  pid <- exp(logd - ma)
  pis <- rowSums(pid)
  r <- ma + log(pis)
  # [log sum(π_j*f(x_1;α_j)), ..., log sum(π_j*f(x_n;α_j))]
  if (!log) r <- exp(r)
  r
}

simplexpanel <- function() {
  sidel <- 2 / sqrt(3)
  liml <- 0.1
  labl <- 0.06
  plot.new()
  plot.window(xlim = c(-liml, sidel + liml), ylim = c(-liml, 1 + liml),
              xaxs = "i", yaxs = "i", asp = 1)

  text(sidel/2, 1 + labl, expression(x[1]))
  text(-labl/sidel, -labl/2, expression(x[2]))
  text(sidel + labl/sidel, -labl/2, expression(x[3]))
  segments(0, 0, sidel, 0)
  segments(0, 0, sidel/2, 1)
  segments(sidel/2, 1, sidel, 0)

  margin <- 0.04
  gap <- c(5, 10, 25, 50, 75, 90, 95) / 100
  temp <- gap / sqrt(3)
  segments(temp, gap, sidel - temp, gap, lty = "dotted", col = "grey")
  segments(temp, gap, 2 * temp, 0, lty = "dotted", col = "grey")
  segments(sidel - temp, gap, sidel - 2 * temp, 0, lty = "dotted", col = "grey")
  text(temp - margin, gap, gap * 100, cex = 0.6)
  text(sidel - temp + margin, gap, (1 - gap) * 100, cex = 0.6)
  text(2 * temp, -margin, (1 - gap) * 100, cex = 0.6)
}

points.simplex <- function(data, add = FALSE, col = 2, pch = 1, cex = 0.5) {
  if (!add) simplexpanel()
  data <- cbind((data[, 3] * 2 + data[, 1])/sqrt(3), data[, 1])
  points(data, col = col, pch = pch, cex = cex)
}

points.dirichletmix <- function(mix, add = FALSE,
                                col = 4, pch = 16, cex = 0.6) {
  data <- (mix$alpha - 1) * mix$h
  points.simplex(data, add = add, col = col, pch = pch, cex = cex)
}

contour.dirichletmix <- function(mix, data = NULL, tfback = TRUE, nlev = 15) {
  x1 <- seq(0, 2/sqrt(3), l = 100)
  x2 <- seq(0, 1, l = 100)
  x <- as.matrix(expand.grid(x1, x2))
  x <- cbind(x[, 2], 1 - (sqrt(3) * x[, 1] + x[, 2])/2)
  x <- cbind(x, 1 - rowSums(x))
  i <- rowSums(x > 0) == ncol(x)

  d <- numeric(length = nrow(x))
  d[i] <- ddirichletmix(x[i, ], mix)
  dm <- matrix(d, nrow = length(x1), ncol = length(x2))

  simplexpanel()
  contour(x1, x2, dm, levels = seq(0, max(d), l = nlev + 1),
          xlab = "x", ylab = "y", col = viridis(nlev), add = TRUE)

  if (!is.null(data)) {
    points.simplex(data, add = TRUE)
  }
}

# initial.mvnmix
initial.simple <- function(data, h, split = 3) {
  l <- seq(0, 1, length = split * 3 + 1)
  loc <- list(l[1:split * 3 - 1], l[1:(split - 1) * 3])
  p <- lapply(loc, function(x) {
    p <- as.matrix(expand.grid(x, x))
    p[rowSums(p) < 1, ]
  })
  p <- do.call(rbind, p)
  p <- cbind(p, 1 - rowSums(p))
  dimnames(p) <- NULL
  m <- nrow(p)

  D <- ncol(data)
  n <- nrow(data)
  dsq <- (data[rep(1:n, m), -D] - p[rep(1:m, each = n), -D])^2
  # square difference between each data point and support point
  d2mu <- rowSums(dsq)
  dim(d2mu) <- c(n, m)
  # each column means the sum square from one data point to one support point
  imu <- apply(d2mu, 1, which.min) # find the closest support point of each data
  nj <- tabulate(imu, nbins = m)
  pr <- nj / sum(nj)
  j <- pr != 0
  dirichletmix(p[j, , drop = FALSE] / h + 1, pr[j], h)
}

# output: log-likelyhood
logLik.dirichletmix <- function(mix, x) {
  sum(ddirichletmix(x, mix, log = TRUE))
}

iter.alpha <- function(x, mix, p, tol = 1e-4, maxit = 100) {
  n <- nrow(x)
  D <- ncol(x)
  a <- mix$alpha
  m <- nrow(a)
  a0 <- sum(a[1, ])
  ll.old <- logLik(mix, data)
  for (j in 1:maxit) {
    dxp <- (log(x)[rep(1:n, m), ] - digamma(a)[rep(1:m, each = n), ]) * c(p)
    # p_ij * [- ψ(α_jk) + log(x_ik)]
    # xk = log(data)[rep(1:n, m),] - digamma(mix$alpha)[rep(1:m,each = n),]
    # xk = array(t(xk), dim=c(n,D,m)) # problem
    # dxp = sweep(xk, c(1,3), p, "*")
    dim(dxp) <- c(n, m * D)
    g <- matrix(colSums(dxp), nrow = m)
    txp <- trigamma(a)[rep(1:m, each = n), ] * c(p) # p_ij * ψ'(α_jk)
    dim(txp) <- c(n, m * D)
    H <- matrix(-colSums(txp), nrow = m)
    for (i in 1:m) {
      h <- diag(H[i, ])
      P <- h %*% (a[i, ] - 1) - g[i, ]
      Q <- -h
      a[i, ] <- pnnqp(Q, P, sum = a0 - D)$x + 1
    }
    mix$alpha <- a
    ll <- logLik(mix, data)
    if (ll <= ll.old + tol) break
    ll.old <- ll
  }
  a
}

EM <- function(data, mix, tem = 5) {
  n <- nrow(data)
  for (iem in 1:tem) {
    logd <- outer.dirichlet(data, mix$alpha, log = TRUE) +
      rep(log(mix$pr), rep(n, length(mix$pr)))
    # log f(x_i;α_j) + log π_j = log [π_j * f(x_i;α_j)]
    ma <- matMaxs(logd)
    pi.d <- exp(logd - ma)
    pi.s <- rowSums(pi.d)
    p <- pi.d / pi.s # p_{ij}
    # -ma or not will not change p
    # which is the same as exp(logd)/rowSums(exp(logd)).
    mix$pr <- colMeans(p) # π_j'
    j <- mix$pr < 1e-10
    if (any(j)) { # delete support points with pr close to 0
      mix$pr <- mix$pr[!j]
      mix$pr <- mix$pr / sum(mix$pr)
      mix$alpha <- mix$alpha[!j, ]
      p <- p[, !j, drop = FALSE]
    }
    mix$alpha <- iter.alpha(data, mix, p) # update alpha
    print(logLik(mix, data))
    # points(mix)
  }
  mix
}

set.seed(1019)
h <- 0.1
mu0 <- rbind(c(0.2, 0.5, 0.3), c(0.7, 0.1, 0.2))
mix0 <- dirichletmix(a = mu0 / h + 1, h = h)
data <- rdirichletmix(50, mix0)
points.simplex(data)
contour(mix0, data)

h <- 0.1
mix <- initial.simple(data, h)
ll <- logLik(mix, data)
mix <- EM(data, mix) # no need combine support points in EM step
ll <- logLik(mix, data)
