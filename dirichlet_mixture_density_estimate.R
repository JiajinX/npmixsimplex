# Update dirichletmix function without h as input
# Function initial.simple with fix support points
# Add back tracking in iter.alpha function
# New function of iter.alpha by Newton-Raphson method
# Modal EM:
# gmultinomix: generate initial points for alpha
# MEM: find the mode of gradient function
# contour.gradient: show the density of gradient and
# tracks from initial points to modes

library(gtools)
library(compositions)
library(viridis)
library(npden)

# Dirichlet mixture parameters
dirichletmix <- function(a = rep(2, 3), pr = 1) {
  if (is.vector(a)) a <- t(a)
  m <- nrow(a)
  D <- ncol(a)
  a0 <- sum(a[1, ])
  pr <- rep(pr, len = m)
  pr <- pr / sum(pr)
  h <- 1 / (a0 - D)
  mix <- list("alpha" = a, "pr" = pr, "h" = h)
  attr(mix, "mode") = (a - 1) / (a0 - D)
  class(mix) <- "dirichletmix"
  mix
}

# generate random data from Dirichlet mixture distribution
rdirichletmix <- function(n, mix) {
  m <- nrow(mix$alpha)
  D <- ncol(mix$alpha)
  idx <- sample.int(m, n, prob = mix$pr, replace = TRUE)
  data <- matrix(nrow = n, ncol = D)
  s = 1
  for (i in 1:m) {
    cnt <- sum(idx == i)
    e <- s + cnt - 1
    data[s:e, ] <- rdirichlet(cnt, mix$alpha[i, ])
    s <- e + 1
  }
  data[sample(n), ]
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

multinomix <- function(p, size, pr = 1) {
  if (is.vector(p)) p <- t(p)
  m <- nrow(p)
  pr <- rep(pr, len = m)
  pr <- pr / sum(pr)
  mix <- list("p" = p, "pr" = pr, "size" = size)
  class(mix) <- "multinomix"
  mix
}

rmultinomix <- function(n, mix) {
  m <- nrow(mix$p)
  D <- ncol(mix$p)
  size <- ceiling(mix$size)
  tb <- table(sample.int(m, n, prob = mix$pr, replace = TRUE))
  cl <- c(0, cumsum(tb))
  point <- matrix(nrow = D, ncol = n)
  for (i in 1:length(tb)) {
    point[, (cl[i] + 1):cl[i + 1]] <-
      rmultinom(tb[i], size, data[as.integer(names(tb[i])), ])
  }
  t(point)
}

gmultinomix <- function(n, mix) {
  point <- rmultinomix(n, mix)
  noise <- matrix(runif(length(point), min = -0.5, max = 0.5), nrow = n)
  noise <- noise - rowMeans(noise)
  point <- point + noise
  point[point < 0] <- 0
  point[point > mix$size] <- mix$size
  point / rowSums(point) * mix$size
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
  l <- matrix(c(0, 0, sidel, 0,
                0, 0, sidel/2, 1,
                sidel/2, 1, sidel, 0),
              nrow = 4)
  segments(l[1, ], l[2, ], l[3, ], l[4, ])

  margin <- 0.04
  gap <- c(5, 10, 25, 50, 75, 90, 95) / 100
  temp <- gap / sqrt(3)
  segments(temp, gap, sidel - temp, gap, lty = 3, col = "grey")
  segments(temp, gap, 2 * temp, 0, lty = 3, col = "grey")
  segments(sidel - temp, gap, sidel - 2 * temp, 0, lty = 3, col = "grey")
  text(temp - margin, gap, gap * 100, cex = 0.6)
  text(sidel - temp + margin, gap, (1 - gap) * 100, cex = 0.6)
  text(2 * temp, -margin, (1 - gap) * 100, cex = 0.6)
}

points.simplex <- function(data, size = 1, add = FALSE,
                           col = 2, pch = 1, cex = 0.5) {
  if (!add) simplexpanel()
  data <- cbind((data[, 3] * 2 + data[, 1])/sqrt(3), data[, 1]) / size
  points(data, col = col, pch = pch, cex = cex)
}

points.dirichletmix <- function(mix, add = FALSE,
                                col = 4, pch = 16, cex = 0.6) {
  data <- attr(mix, "mode")
  points.simplex(data, add = add, col = col, pch = pch, cex = cex)
}

lines.simplex <- function(x, y, size = 1, arrow = FALSE, add = FALSE,
                          lty = 3, col = "gray40", length = 0.05) {
  if (!add) simplexpanel()
  x <- cbind((x[, 3] * 2 + x[, 1])/sqrt(3), x[, 1]) / size
  y <- cbind((y[, 3] * 2 + y[, 1])/sqrt(3), y[, 1]) / size
  if (arrow) {
    suppressWarnings(
      arrows(x[, 1], x[, 2], y[, 1], y[, 2],
             col = col, lty = lty, length = length))
  } else {
    segments(x[, 1], x[, 2], y[, 1], y[, 2], col = col, lty = lty)
  }
}

contour.dirichletmix <- function(mix, data = NULL, lv = NULL, nlev = 15) {
  x1 <- seq(0, 2/sqrt(3), l = 100)
  x2 <- seq(0, 1, l = 100)
  x <- as.matrix(expand.grid(x1, x2))
  x <- cbind(x[, 2], 1 - (sqrt(3) * x[, 1] + x[, 2])/2)
  x <- cbind(x, 1 - rowSums(x))
  i <- rowSums(x > 0) == ncol(x)

  d <- numeric(length = nrow(x))
  d[i] <- ddirichletmix(x[i, ], mix)
  dm <- matrix(d, nrow = length(x1), ncol = length(x2))
  if (is.null(lv)) lv <- seq(0, max(d), l = nlev + 1)

  simplexpanel()
  contour(x1, x2, dm, levels = lv, xlab = "x", ylab = "y",
          col = viridis(nlev), add = TRUE)

  if (!is.null(data)) {
    points.simplex(data, add = TRUE)
  }

  attr(d, "level") <- lv
  d
}

contour.gradient <- function(data, w, size, lv = NULL, nlev = 15) {
  x1 <- seq(0, 2/sqrt(3), l = 100)
  x2 <- seq(0, 1, l = 100)
  x <- as.matrix(expand.grid(x1, x2))
  x <- cbind(x[, 2], 1 - (sqrt(3) * x[, 1] + x[, 2])/2)
  x <- cbind(x, 1 - rowSums(x))
  i <- rowSums(x > 0) == ncol(x)

  d <- numeric(length = nrow(x))
  logd <- outer.dirichlet(data, x[i,] * size, log = TRUE) - log(w)
  d[i] <- colSums(exp(logd)) # - nrow(data)
  # Gradient value of alpha, no need to -n (nrow(data)) here
  dm <- matrix(d, nrow = length(x1), ncol = length(x2))
  if (is.null(lv)) lv <- seq(min(d), max(d), l = nlev + 1)

  simplexpanel()
  sidel <- 2 / sqrt(3)
  l1 <- 1 / size
  l2 <- l1 * sqrt(3)
  l3 <- l1 * 2
  l <- matrix(c(l2, l1, sidel - l2, l1,
                l2, l1, sidel/2, 1 - l3,
                sidel/2, 1 - l3, sidel - l2, l1),
              nrow = 4)
  segments(l[1, ], l[2, ], l[3, ], l[4, ], col = "#0000CD", lty = 2)
  contour(x1, x2, dm, levels = lv, xlab = "x", ylab = "y",
          col = viridis(nlev), add = TRUE)
}

print_verbose <- function(verbose, i, ll, mix) {
  if (verbose >= 1) cat("\n## Iteration ", i, ":\nLog-likelihood = ",
                       as.character(signif(ll[1], 16)), "\n", sep = "")
  if (verbose >= 2) {cat("Component proportions:\n"); print(mix$pr)}
  if (verbose >= 3) {cat("Component alpha:\n"); print(mix$alpha)}
  if (verbose >= 4) {cat("Component mode:\n"); print(attr(mix, "mode"))}
  if (verbose >= 5) {cat("Component bandwith:\n"); print(mix$h)}
}

print.dirichletmix <- function(mix) {
  i <- order(mix$alpha[, 1])
  mix$alpha <- mix$alpha[i, ]
  mix$pr <- mix$pr[i]
  attr(mix, "mode") <- attr(mix, "mode")[i, ]
  class(mix) <- NULL
  print(mix)
}

# initial.mvnmix
initial.simple <- function(data, np, h) {
  split <- ceiling(sqrt(np))
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
  drop <- m - np
  j <- order(pr)[(drop + 1):m] # drop points with small pr
  dirichletmix(p[j, , drop = FALSE] / h + 1, pr[j])
}

# output: log-likelyhood
logLik.dirichletmix <- function(mix, x) {
  sum(ddirichletmix(x, mix, log = TRUE))
}

calfun.alpha <- function(x, a, p, deriv = 0) {
  n <- nrow(x)
  D <- ncol(x)
  m <- nrow(a)
  xp <-
    switch(deriv + 1,
           rowSums((a %*% t(log(x)) - rowSums(lgamma(a))) * t(p)),
           (log(x)[rep(1:n, m), ] - digamma(a)[rep(1:m, each = n), ]) * c(p),
           # p_ij * [log(x_ik) - ψ(α_jk)]
           -trigamma(a)[rep(1:m, each = n), ] * c(p)
           # -p_ij * ψ'(α_jk)
    )
  if (deriv != 0) {
    dim(xp) <- c(n, m * D)
    xp <- matrix(colSums(xp), nrow = m)
  }
  xp
}

# Iteration alpha by using pnnqp
iter.alpha <- function(x, a, p, tol = 1e-4, maxit = 100) {
  D <- ncol(x)
  a.old <- a
  m <- nrow(a)
  a0 <- sum(a[1, ])
  fa <- fa.old <- calfun.alpha(x, a, p)
  # cat("\nalpha-function:"); print(as.character(signif(sum(fa), 16)))
  for (j in 1:maxit) {
    g <- calfun.alpha(x, a, p, deriv = 1)
    H <- calfun.alpha(x, a, p, deriv = 2)
    for (i in 1:m) {
      h <- diag(H[i, ])
      P <- h %*% (a[i, ] - 1) - g[i, ]
      Q <- -h
      a[i, ] <- a.up <- pnnqp(Q, P, sum = a0 - D)$x + 1
      dea <- sum(g[i, ] * (a[i, ] - a.old[i, ])) * 0.333
      alpha <- 1
      for (k in 1:50) { # back tracking
        fa[i] <- sum((a[i, ] %*% t(log(x)) - sum(lgamma(a[i, ]))) * p[, i])
        if (fa[i] >= fa.old[i] + dea * alpha) break
        alpha <- alpha * 0.5
        a[i, ] <- (1 - alpha) * a.old[i, ] + alpha * a.up
      }
    }
    # cat("\nalpha-function:"); print(as.character(signif(sum(fa), 16)))
    # if (sum(fa) < sum(fa.old)) cat("Decreasing-----------------------------------\n")
    if (sum(fa) <= sum(fa.old) + tol) break
    fa.old <- fa
    a.old <- a
  }
  a
}

iter.alpha <- function(x, a, p, tol = 1e-4, maxit = 100) {
  D <- ncol(x)
  m <- nrow(a)
  a0 <- sum(a[1, ])
  fa.old <- fa <- sum(calfun.alpha(x, a, p))
  # cat("\nalpha-function:"); print(as.character(signif(fa, 16)))
  for (j in 1:maxit) {
    g <- calfun.alpha(x, a, p, deriv = 1)
    H <- calfun.alpha(x, a, p, deriv = 2)
    for (i in 1:m) {
      h <- diag(H[i, ])
      P <- h %*% (a[i, ] - 1) - g[i, ]
      Q <- -h
      a[i, ] <- pnnqp(Q, P, sum = a0 - D)$x + 1
    }
    fa <- sum(calfun.alpha(x, a, p))
    # cat("\nalpha-function:"); print(as.character(signif(fa, 16)))
    # if (fa < fa.old) cat("Decreasing-----------------------------------\n")
    if (fa <= fa.old + tol) break
    fa.old <- fa
  }
  a
}

stepup <- function(x, update) {
  a <- x + update
  j <- which(rowSums(a < 1) > 0) # Find components with alpha' < 1
  if (!identical(j, integer(0))) {
    # If a_k - 1 < 1e-10; let a_k = 1
    sx <- x[j, , drop = FALSE]
    supdate <- update[j, , drop = FALSE]
    lf <- (sx - 1) / -supdate
    # For a_k = 1, update_k will be >= 0
    # when a_k = 1 & update_k = 0, lf = NAN
    # when a_k = 1 & update_k > 0, lf = 0
    lf <- apply(lf, 1, function(r) min(r[r > 0], na.rm = TRUE))
    a[j, ] <- sx + lf * supdate # Back tracking, pull back to boundary
  }
  a
}

# Iteration alpha by using Newton-Raphson method
iter.alpha <- function(x, a, p, maxit = 100, tol = 1e-4) {
  fa.old <- fa <- sum(calfun.alpha(x, a, p))
  # cat("\nalpha-function:"); print(as.character(signif(fa, 16)))
  for (i in 1:maxit) {
    isbound <- a == 1
    # TRUE for "1_set" (on the boundary), FALSE for ">1_set"
    g <- calfun.alpha(x, a, p, deriv = 1)
    H <- calfun.alpha(x, a, p, deriv = 2)
    iH <- iH_m <- 1 / H
    iHg <- iHg_m <- g * iH

    # For ">1_set", no "1_set" elements added
    iH_m[isbound] <- iHg_m[isbound] <- NA
    siH_m <- rowSums(iH_m, na.rm = TRUE)
    siHg_m <- rowSums(iHg_m, na.rm = TRUE)
    lamb <- siHg_m / siH_m # Lambda for all ">1_set" elements

    # Add each "1_set" element into ">1_set" to calculate gradient
    j <- which(rowSums(isbound) > 0) # Find rows have alpha on boundary
    if (!identical(j, integer(0))) {
      iH_1 <- iH[j, , drop = FALSE]
      iHg_1 <- iHg[j, , drop = FALSE]
      iH_1[!isbound[j, ]] <- iHg_1[!isbound[j, ]] <- NA
      siH_1 <- siH_m[j] + iH_1
      siHg_1 <- siHg_m[j] + iHg_1
      lamb_1 <- siHg_1 / siH_1
      trend <- g[j, ] - lamb_1
      trend[trend < 0] <- NA
      rmax <- apply(trend, 1, function(r) {
        rmax <- which.max(r)
        if (identical(rmax, integer(0))) rmax <- NA
        rmax
      })
      maxloc <- cbind(seq(length(j)), j, rmax)
      maxloc <- maxloc[!is.na(maxloc[,3]), , drop = FALSE]
      # Find one in "1_set" with largest gradient value (>0) and move to ">1_set"

      isbound[maxloc[, -1, drop = FALSE]] <- FALSE
      lamb[maxloc[, 2]] <- lamb_1[maxloc[, -2, drop = FALSE]]
    }
    iH[isbound] <- iHg[isbound] <- 0
    lambH <- lamb * iH
    update <- lambH - iHg
    a <- stepup(a, update)
    fa <- sum(calfun.alpha(x, a, p))
    # cat("\nalpha-function:"); print(as.character(signif(fa, 16)))
    # if (fa < fa.old) cat("Decreasing-----------------------------------\n")
    if (fa <= fa.old + tol) break
    fa.old <- fa
  } # The alpha-function seems to be monotonically increasing for each iteration
  a
}

EM <- function(data, mix, tem = 5, tol = 1e-4, verbose = 0) {
  n <- nrow(data)
  ll <- ll.old <- logLik(mix, data)
  for (iem in 1:tem) {
    if (verbose > 0) print_verbose(verbose, iem - 1, ll, mix)
    logd <- outer.dirichlet(data, mix$alpha, log = TRUE) +
      rep(log(mix$pr), rep(n, length(mix$pr)))
    # log f(x_i;α_j) + log π_j = log [π_j * f(x_i;α_j)]
    ma <- matMaxs(logd)
    pi.d <- exp(logd - ma)
    pi.s <- rowSums(pi.d)
    p <- pi.d / pi.s # p_{ij} log scale !!!!
    # -ma or not will not change p
    # which is the same as exp(logd) / rowSums(exp(logd)).
    mix$pr <- colMeans(p) # π_j'
    j <- mix$pr < 1e-10
    if (any(j)) { # delete support points with pr close to 0
      mix$pr <- mix$pr[!j]
      mix$pr <- mix$pr / sum(mix$pr)
      mix$alpha <- mix$alpha[!j, ]
      p <- p[, !j, drop = FALSE]
    }
    mix$alpha <- iter.alpha(data, mix$alpha, p)
    attr(mix, "mode") = (mix$alpha - 1) * mix$h
    ll <- logLik(mix, data)
    if (ll <= ll.old + tol) break
    ll.old <- ll
  }
  if (verbose > 0) print_verbose(verbose, iem, ll, mix)
  mix
}

MEM <- function(data, a, mix, plot = FALSE, tem = 5, tol = 1e-4) {
  n <- nrow(data)
  w <- pmax(ddirichletmix(data, mix), 1e-100)
  size <- sum(mix$alpha[1, ])
  if (plot) contour.gradient(data, w, size)
  g <- g.old <- -Inf
  a.old <- a
  for (i in 1:tem) {
    logd <- outer.dirichlet(data, a, log = TRUE) - log(w)
    g <- sum(colSums(exp(logd)) - n)
    if (g <= g.old + tol) break
    g.old <- g

    ma <- matMaxs(logd, dim = 2)
    pi.d <- exp(logd - rep(ma, each = n))
    pi.s <- colSums(pi.d)
    p <- pi.d / rep(pi.s, each = n)
    a <- iter.alpha(data, a, p)
    if (plot) {
      arrow <- i == 1
      lty <- ifelse(arrow, 1, 3)
      lines.simplex(a.old, a, size, arrow = arrow, add = TRUE, lty = lty)
    }
    a.old <- a
  }
  if (plot) points.simplex(a, size, add = TRUE, pch = 19, cex = 0.65)
  a
}

set.seed(1019)
h <- 0.1
n <- 50
mu0 <- rbind(c(0.2, 0.5, 0.3), c(0.7, 0.1, 0.2))
mix0 <- dirichletmix(a = mu0 / h + 1)
data <- rdirichletmix(n, mix0)
# points.simplex(data)
# contour(mix0, data)

h <- 0.1
mix <- initial.simple(data, np = 2, h)
mix <- EM(data, mix, verbose = 1) # no need combine support points in EM step
ll <- logLik(mix, data)

# Finding new support points
nt <- 20
w <- pmax(ddirichletmix(data, mix), 1e-100) # f(x_i;G) for robust
gmix <- multinomix(data, 1/mix$h, 1/w)
t <- rbind(mix$alpha, gmultinomix(n = nt, gmix) + 1) # nt = 20 random grid
t <- MEM(data, t, mix, tem = 100, plot = TRUE)
# Because the limitation of alpha >= 1, some results points are not in center

# test EM algorithm ------------------------------------------------------------

cs_mix <- function(no) {
  p <- switch(no,
              list(mu = rbind(c(0.2, 0.5, 0.3), c(0.7, 0.1, 0.2),
                              c(0.4, 0.4, 0.2), c(0.1, 0.7, 0.2)),
                   pr = c(0.1, 0.4, 0.2, 0.3), h = 0.1),
              list(mu = rbind(c(0.2, 0.5, 0.3), c(0.7, 0.1, 0.2),
                              c(0.4, 0.4, 0.2), c(0.9, 0.05, 0.05)),
                   pr = c(0.1, 0.4, 0.2, 0.3), h = 0.1),
              list(mu = rbind(c(0, 0.2, 0.8), c(1, 0, 0),
                              c(1/3, 1/3, 1/3), c(0.6, 0.3, 0.1)),
                   pr = c(0.1, 0.2, 0.3, 0.4), h = 0.1),
              list(mu = rbind(c(0.2, 0.5, 0.3), c(0.7, 0.1, 0.2),
                              c(0.4, 0.4, 0.2)),
                   pr = c(0.5, 0.2, 0.3), h = 0.2)
              )
  mix0 <- dirichletmix(a = p$mu / p$h + 1, pr = p$pr)
  mix0
}

compare_em <- function(mix0, mix) {
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(2, 2), mar = rep(0, 4))
  points(mix0)
  d0 <- contour(mix0)
  points(mix)
  contour(mix, lv = attr(d0, "level"))
  par(old_par)
}


no <- 3
n <- 100
mix0 <- cs_mix(no)
data <- rdirichletmix(n, mix0)

mix <- initial.simple(data, np = 4, h = 0.1)
mix <- EM(data, mix, verbose = 1)
logLik(mix, data)

nt <- 20
w <- pmax(ddirichletmix(data, mix), 1e-100) # f(x_i;G) for robust
gmix <- multinomix(data, 1/mix$h, 1/w)
t <- rbind(mix$alpha, gmultinomix(n = nt, gmix) + 1) # nt = 20 random grid
t <- MEM(data, t, mix, tem = 100, plot = TRUE)

