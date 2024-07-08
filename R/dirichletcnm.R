# Generate initial support points of Dirichlet mixture
initial.dirichletmix <- function(data, np, h, by) {
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

  n <- nrow(data)
  d2mu <- rowSums((data[rep(1:n, m), ] - p[rep(1:m, each = n), ])^2)
  # square difference between each data point and support point
  dim(d2mu) <- c(n, m)
  # each column means the sum square from one data point to one support point
  imu <- apply(d2mu, 1, which.min) # find the closest support point of each data
  nj <- tabulate(imu, nbins = m)
  pr <- nj / sum(nj)
  drop <- m - np
  j <- order(pr)[(drop + 1):m] # drop points with small pr
  dirichletmix(p[j, , drop = FALSE] / h + 1, pr[j], by = by)
}

# EM: objective function needs to be maximized by iteration alpha
# Also with the first and second derivative of objective function
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

# EM: updating alpha by using pnnqp method and back tracking with Armijo rule
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
      a[i, ] <- a.up <- lsei::pnnqp(Q, P, sum = a0 - D)$x + 1
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
    # if (sum(fa) < sum(fa.old)) cat("Decreasing----------------------------\n")
    if (sum(fa) <= sum(fa.old) + tol) break
    fa.old <- fa
    a.old <- a
  }
  a
}

# EM: updating alpha by using pnnqp method
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
      a[i, ] <- lsei::pnnqp(Q, P, sum = a0 - D)$x + 1
    }
    fa <- sum(calfun.alpha(x, a, p))
    # cat("\nalpha-function:"); print(as.character(signif(fa, 16)))
    # if (fa < fa.old) cat("Decreasing--------------------------------------\n")
    if (fa <= fa.old + tol) break
    fa.old <- fa
  }
  a
}

# Cheack if alpha points out of boundary after updating and do back tracking
stepup <- function(x, update) {
  a <- x + update
  j <- which(rowSums(a < 1) > 0) # Find components with alpha' < 1
  if (!identical(j, integer(0))) {
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

# EM: updating alpha by using Newton-Raphson method
iter.alpha <- function(x, a, p, maxit = 100, tol = 1e-4) {
  fa.old <- fa <- sum(calfun.alpha(x, a, p))
  # cat("\nalpha-function:"); print(as.character(signif(fa, 16)))
  for (i in 1:maxit) {
    a[abs(a - 1) < 1e-10] <- 1 # round ?????
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
      rmax <- apply(trend, 1, which.max)
      maxloc <- cbind(seq_along(j), j, rmax)
      maxloc <- maxloc[trend[maxloc[, -2, drop = FALSE]] > 0, , drop = FALSE]
      # Find one in "1_set" with largest gradient positive value
      # and move to ">1_set"

      isbound[maxloc[, -1, drop = FALSE]] <- FALSE
      lamb[maxloc[, 2]] <- lamb_1[maxloc[, -2, drop = FALSE]]
    }
    iH[isbound] <- iHg[isbound] <- 0
    lambH <- lamb * iH
    update <- lambH - iHg
    a <- stepup(a, update)
    fa <- sum(calfun.alpha(x, a, p))
    # cat("\nalpha-function:"); print(as.character(signif(fa, 16)))
    # if (fa < fa.old) cat("Decreasing--------------------------------------\n")
    if (fa <= fa.old + tol) break
    fa.old <- fa
  } # The alpha-function seems to be monotonically increasing for each iteration
  a
}

# EM algorithm of Dirichlet mixture
EM <- function(data, mix, tem = 5, tol = 1e-4, verbose = 0) {
  n <- nrow(data)
  ll <- ll.old <- logLik(mix, data)
  for (i in 1:tem) {
    if (verbose > 0) print_verbose(verbose, i - 1, ll, mix)
    logd <- outer.dirichlet(data, mix$alpha, log = TRUE) +
      rep(log(mix$pr), rep(n, length(mix$pr)))
    # log f(x_i;α_j) + log π_j = log [π_j * f(x_i;α_j)]
    ma <- lsei::matMaxs(logd)
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
    attr(mix, "mode") <- (mix$alpha - 1) * mix$h
    ll <- logLik(mix, data)
    if (ll <= ll.old + tol) break
    ll.old <- ll
  }
  if (verbose > 0) print_verbose(verbose, i, ll, mix)
  mix
}

# Generate initial alpha points of Modal EM
gmultinomix <- function(n, size, pr, mass) {
  if (is.vector(pr)) pr <- t(pr)
  m <- nrow(pr)
  D <- ncol(pr)
  mass <- rep(mass, len = m)
  mass <- mass / sum(mass)
  insz <- ceiling(size)
  tb <- table(sample.int(m, n, prob = mass, replace = TRUE))
  cl <- c(0, cumsum(tb))
  point <- matrix(nrow = D, ncol = n)
  for (i in seq_along(tb)) {
    point[, (cl[i] + 1):cl[i + 1]] <-
      rmultinom(tb[i], insz, pr[as.integer(names(tb[i])), ])
  }
  point <- t(point)
  noise <- matrix(runif(length(point), min = -0.5, max = 0.5), nrow = n)
  noise <- noise - rowMeans(noise)
  point <- point + noise
  point[point < 0] <- 0
  point[point > size] <- size
  point / rowSums(point) * size
}

# Modal EM algorithm of the gradient function of
# Dirichlet mixture's log-likelihood
MEM <- function(data, a, w, plot = FALSE, tem = 5, tol = 1e-4) {
  n <- nrow(data)
  size <- sum(a[1, ])
  if (plot) contour.gradient(data, w, size)
  sg <- sg.old <- -Inf
  a.old <- a
  for (i in 1:tem) {
    logd <- outer.dirichlet(data, a, log = TRUE) - w
    sg <- sum(colSums(exp(logd)))
    if (sg <= sg.old + tol) break
    sg.old <- sg

    ma <- lsei::matMaxs(logd, dim = 2)
    pi.d <- exp(logd - rep(ma, each = n))
    pi.s <- colSums(pi.d)
    p <- pi.d / rep(pi.s, each = n)
    a <- iter.alpha(data, a, p)
    if (plot) segments.simplex(a.old, a, size, add = TRUE)
    a.old <- a
  }
  if (plot) points.simplex(a, size, add = TRUE, pch = 19, cex = 0.65)
  logd <- outer.dirichlet(data, a, log = TRUE) - w
  g <- colSums(exp(logd)) - n
  list("alpha" = a, "gradient" = g)
}

# Select new support points with high gradient value
select_spt <- function(tg, a, by = 1) {
  o <- order(a[, by])
  a <- a[o, , drop = FALSE]
  mids <- (a[-1, by] + a[-nrow(a), by]) * 0.5
  # Middle points of each pair support points on one coordinate
  t <- tg$alpha
  o <- order(t[, by])
  t <- t[o, , drop = FALSE] # order t by increasing
  index <- lsei::indx(t[, by], mids) + 1
  # Split t into several groups by interval of mids pairs
  g <- tg$gradient[o]
  jj <- aggregate(g, by = list(group = index), which.max)
  # Find the index of largest gradient value of each group
  j <- match(jj$group, index) + jj$x - 1 # transform index
  np <- t[j, , drop = FALSE]
  gnp <- g[j]
  list("alpha" = np, "gradient" = gnp)
  # New support points with their gradient values
}

# Update mass of components of Dirichlet mixture
CNM <- function(np, mix, data, ll, w, by = 1) {
  if (missing(w)) w <- ddirichletmix(data, mix, log = TRUE) # log f(x_i; G)
  mix.raw <- dirichletmix(rbind(mix$alpha, np), c(mix$pr, double(length(np))),
                          by = by)
  # Add new support points and set their mass = 0
  logd <- outer.dirichlet(data, mix.raw$alpha, log = TRUE) # f(x_i; α_j)
  s <- exp(logd - w) # First-derivative of ll with respect to π
  g <- colSums(s)
  r <- lsei::pnnls(s, 2, sum = 1)
  # min||Sπ'-2||^2: a = S, b = 2 with π'1 = 1, π' >= 0
  sol <- r$x / sum(r$x) # ensure sum to 1
  mix.up <- mix.raw
  mix.up$pr <- sol # update π
  dea <- sum(g * (mix.up$pr - mix.raw$pr)) * .333
  mix <- mix.up
  ll.old <- ll
  alpha <- 1
  for (i in 1:50) {
    ll <- logLik(mix, data)
    if (ll >= ll.old + dea * alpha) break
    alpha <- alpha * 0.5
    mix$pr <- (1 - alpha) * mix.raw$pr + alpha * mix.up$pr
  } # Back tracking with Armijo rule
  list("mix" = mix, "ll" = ll)
}

# Drop and merge the components of Dirichlet mixture
collapse_spt <- function(mix, tol = 0.1, by = 1) {
  m <- length(mix$pr)
  if (m == 1) return(mix)
  mix <- sort(mix, by = by)
  j <- mix$pr != 0
  mix$alpha <- mix$alpha[j, , drop = FALSE]
  mix$pr <- mix$pr[j]
  if (m == 1 || tol <= 0) return(mix)

  tol2 <- tol * tol
  m <- nrow(mix$alpha)
  d2 <- rowSums((mix$alpha[rep(1:m, m), ] - mix$alpha[rep(1:m, each = m), ])^2)
  # distance of each pair of alpha (change the type of distance ?????)
  dim(d2) <- c(m, m)
  d2l <- d2[lower.tri(d2)] # Extract the values of the lower triangular matrix
  d2i <- cbind(d2l, t(combn(m, 2))) # Add index of pairs
  d <- d2i[order(d2i[, 1]), , drop = FALSE]

  a <- mix$alpha
  pr <- mix$pr
  index <- rep(TRUE, m)
  for (i in seq_along(d)) {
    if (d[i, 1] > tol2) break # If distance > tol2, no treatment
    ip <- d[i, 2:3] # indices of the pair
    if (all(index[ip])) {
      mp <- sum(pr[ip])
      a[ip[1], ] <- colSums(a[ip, , drop = FALSE] * pr[ip] / mp)
      # Weighted averaging of alpha; weighted.mean.array ?????
      pr[ip[1]] <- mp # Merge alpha and p to first index
      index[ip[2]] <- FALSE # Drop second index
    }
  }
  mix$alpha <- a[index, , drop = FALSE]
  mix$pr <- pr[index]
  attr(mix, "mode") <- (mix$alpha - 1) * mix$h
  mix
}

# Overall algorithm
cnmdm <- function(data, h, nump, verbose = 0,
                  maxit = 100, nt = 20, tol = 1e-4, ctol = 0.2) {
  if (!is.numeric(data)) data <- as.matrix(data)
  # if (missing(by))
  by <- which.min(colSums(apply(data, 2, duplicated)))
  # Initial the mixture
  mix <- initial.dirichletmix(data, nump, h, by)
  ll <- ll.old <- logLik(mix, data)
  convergence <- 1
  # tem <- tems <- 5
  tem <- 5
  for (i in 1:maxit) {
    if (verbose > 0) print_verbose(verbose, i - 1, ll, mix)
    # EM algorithm
    mix <- EM(data, mix, tem = tem) # no need combine support points in EM step
    # if (tem <= tems) ll <- logLik(mix, data) # ???
    ll <- logLik(mix, data)
    # Expand new support points
    w <- ddirichletmix(data, mix, log = TRUE) # log f(x_i;G) for robust
    t <- rbind(mix$alpha, gmultinomix(nt, 1 / mix$h, data, 1 / exp(w)) + 1)
    tg <- MEM(data, t, w) # Blue line is the limitation of alpha >= 1
    st <- select_spt(tg, mix$alpha, by)
    np <- st$alpha
    # points.simplex(np, size = sum(np[1, ]), add = TRUE,
    # col = "blue", pch = 19, cex = 0.65)
    # Update components' mass
    mixll <- CNM(np, mix, data, ll, w, by)
    mix <- mixll$mix
    ll <- mixll$ll
    # Collapse support points
    mix <- collapse_spt(mix, by = by, tol = ctol) # tol = ???
    if (ll <= ll.old + tol) {
      if (tem == 5) {
        tem <- 100
      } else { # if (tem == 100) { # ???
        convergence <- 0
        break
      }
    } else {
      tem <- 5
    }
    ll.old <- ll
  }
  if (verbose > 0) print_verbose(verbose, i, ll, mix)
  structure(list(num.iterations = i, convergence = convergence,
                 max.gradient = max(st$gradient), ll = ll, mix = mix, h = h),
            class = "cnmdm")
}
