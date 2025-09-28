#  Step function of nonparamatrix Dirichlet mixture density estimation

# Check if points out of simplex boundary after updating and do back tracking
stepup <- function(old, update, bound = 0) {
  if (is.vector(old)) old <- t(old)
  if (is.vector(update)) update <- t(update)
  new <- old + update
  j <- which(rowSums(new < bound) > 0) # Find components with new < boundary
  if (!identical(j, integer(0))) { # Pull new back to boundary
    sold <- old[j, , drop = FALSE]
    supdate <- update[j, , drop = FALSE]
    lf <- (sold - bound) / -supdate
    # Case of update alpha:
    # For a_k = 1, update_k will be >= 0:
    # when a_k = 1 & update_k = 0, lf = NAN
    # when a_k = 1 & update_k > 0, lf = 0
    lf[lf < 0 | is.na(lf)] <- Inf
    # lf < 0: a_k > 1 & update_k > 0
    # lf = NA: update_k = 0
    rmin <- apply(lf, 1, function(r) which.min(r))
    minloc <- cbind(seq_along(j), j, rmin)
    new[j, ] <- sold + lf[minloc[, -2, drop = FALSE]] * supdate
    # Back tracking, pull back to boundary
    new[minloc[, -1, drop = FALSE]] <- bound
  }
  new
}

# Initialize Dirichlet mixture by grouping data (simple clustering)
# h: bandwidth; by: index of column used to sort components
initial_dirichletmix <- function(data, h, by) {
  n <- nrow(data)
  D <- ncol(data)
  a0 <- 1 / h + D
  cst <- 2 * sqrt(a0 + 1) # Constant value in transformation
  tdata <- cst * asin(sqrt(data))
  # Variance stabilization transformation
  tmu <- cst * asin(sqrt(colMeans(data)))
  # tmu should be the first center of cluster groups
  d2 <- d2mu <- colSums((t(tdata) - tmu)^2) # Distance from tdata to tmu
  d <- sqrt(d2)
  index <- 1:n
  pts <- matrix(nrow = n, ncol = D)
  # r2 <- ifelse(D > 10, qchisq(0.95, df = D), 9)
  # exout <- sqrt(r2) * 4 / 3
  r2 <- 9 # Distance threshold
  exout <- 4
  # tfsimplexpanel(a0)
  # plot3d(t(tmu), col = 1, add = T, size = 6)
  for (k in seq_len(n + 1)) { # rep = number of data points + 1 (tmu)
    i <- d2 > r2
    # All points with distance < r2 will be in a cluster group with
    # pts[k, ] as the center and r2 as the radius of circle
    # plot3d(tdata[index, , drop = F][!i, , drop = F], add = T,
    #        col = col[k %% length(col) + 1], size = 4)
    index <- index[i] # Uncovered data points
    if (length(index) == 0) break # If no uncovered points then return
    imin <- index[which.min(d2mu[index])]
    update <- (tdata[imin, ] - tmu) * sqrt(1 + exout / d[imin])
    pts[k, ] <- tmu + update
    # Find the closet uncovered point to the tmu, expand this point
    # outward a certain distance as the new cluster center
    if (any(sum(pts[k, ] < 0))) {
      # If points outside the boundary, pull it back
      lf <- tmu / -update
      lf <- min(lf[lf >= 0])
      pts[k, ] <- tmu + lf * update # numeric issue for results close to 1 but < 1
    }
    d2 <- colSums((t(tdata)[, index, drop = FALSE] - pts[k, ])^2)
    # Distance for uncovered points to new cluster center
  }
  tmu <- if (k == 1) t(tmu) else rbind(tmu, pts[1:(k - 1), , drop = FALSE])
  # Set tmu and pts as center of clusters
  rownames(tmu) <- NULL
  # plot3d(tmu, col = 1, add = T, size = 6)
  d2mu <- rowSums((tdata[rep(1:n, k), ] - tmu[rep(1:k, each = n), ])^2)
  dim(d2mu) <- c(n, k)
  imu <- apply(d2mu, 1, which.min)
  nj <- tabulate(imu, nbins = k)
  pr <- nj / sum(nj) # Calculate pr by distance of data points and tmu
  j <- pr != 0
  mu <- sin(tmu[j, , drop = FALSE] / cst)^2 # Transform back
  mu <- mu / rowSums(mu) # mu should > 1 / a0
  update <- (mu * D - 1) * h # Update from mu to mode
  # Transform tmu back to the simplex, which are the means (a_i / a0) of the
  # Dirichlet components, then got the modes by means
  mode <- stepup(mu, update)
  alpha <- mode / h + 1
  mix <- dirichletmix(alpha, pr[j], by = by)
  attr(mix, "ll") <- logLik(mix, data)
  mix
}

# Objective function needs to be maximized by iteration alpha (EM -> iter_alpha)
# Also with the first and second derivative of objective function
# x: data; a: alpha; p: value get from EM
# hete: FALSE - a0 is fixed (used in MDE)
#       TRUE - no restrict of a0 (used in finite mixture)
calfun_alpha <- function(x, a, p, deriv = 0, hete = FALSE) {
  n <- nrow(x)
  D <- ncol(x)
  m <- nrow(a)
  logx <- attr(x, "log")
  if (is.null(logx)) logx <- pmax(log(x), log(1e-100))
  if (deriv != 0) dim(p) <- NULL
  xp <-
    switch(deriv + 1,
      rowSums((tcrossprod(a, logx) - rowSums(lgamma(a))) * t(p)),
      # Original objective function
      (logx[rep(1:n, m), ] - digamma(a)[rep(1:m, each = n), ]) * p,
      # p_ij * [log(x_ik) - ψ(α_jk)] - first derivative
      -trigamma(a)[rep(1:m, each = n), ] * p
      # -p_ij * ψ'(α_jk) - second derivative
    )
  if (deriv != 0) {
    dim(xp) <- c(n, m * D)
    xp <- colSums(xp)
    dim(xp) <- c(m, D)
  }
  if (hete) {
    a0 <- rowSums(a)
    dim(p) <- c(n, m)
    at <- switch(deriv + 1, # Add a0 term when a0 is not fixed
                 rowSums(t(p) * lgamma(a0)),
                 rowSums(t(p) * digamma(a0)),
                 rowSums(t(p) * trigamma(a0)))
    xp <- if (deriv == 2) {
      lapply(1:m, function(i) diag(xp[i, ]) + at[i])
    } else xp + at
  }
  xp
}

# Updating alpha by using Newton-Raphson method (EM)
# x: data; a: alpha; p: value get from EM
# hete: FALSE - a0 is fixed (used in MDE)
#       TRUE - no restrict of a0 (used in finite mixture)
iter_alpha <- function(x, a, p, hete = FALSE,
                       maxit = 100, maxbt = 50, tol = 1e-4) {
  m <- nrow(a)
  obj.old <- obj <- calfun_alpha(x, a, p, hete = hete)
  fa.old <- fa <- sum(obj)
  # cat("\nalpha-function:"); print(fa, digits = 16)
  for (i in seq_len(maxit)) {
    isbound <- a == 1
    # TRUE for "1_set" (on the boundary), FALSE for ">1_set"
    g <- calfun_alpha(x, a, p, deriv = 1, hete = hete)
    H <- calfun_alpha(x, a, p, deriv = 2, hete = hete)

    if (hete) { # If no restrict of constant a0
      a.old <- a
      a <- t(sapply(1:m, function(comp) {
        P <- H[[comp]] %*% (a[comp, ] - 1) - g[comp, ]
        Q <- -H[[comp]]
        lsei::pnnqp(Q, P)$x + 1
      }))
    } else { # a0 is a constant value
      # First and second derivative of objective function
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
        # Set the elements not = 1 to NA
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
      a <- stepup(a, update, bound = 1) # Update alpha
    }
    obj <- calfun_alpha(x, a, p, hete = hete)

    if (hete) { # Back Tracking part for heterogeneous finite mixture
      ibt <- which(obj < obj.old) # The True components need back tracking
      if (!identical(ibt, integer(0))) {
        dea <- 0.333 *
          rowSums(g[ibt, , drop = FALSE] *
                    (a[ibt, , drop = FALSE] - a.old[ibt, , drop = FALSE]))
        alpha <- 1
        a.up <- a
        for (k in seq_len(maxbt)) { # back tracking
          obj[ibt] <- calfun_alpha(x, a[ibt, , drop = FALSE],
                                   p[, ibt, drop = FALSE], hete = hete)
          idx <- which(obj[ibt] < obj.old[ibt] + dea * alpha)
          # TRUE for components still need back tracking
          if (identical(idx, integer(0))) break
          ibt <- ibt[idx]
          dea <- dea[idx]
          alpha <- alpha * 0.5
          a[ibt, ] <- (1 - alpha) * a.old[ibt, , drop = FALSE] +
            alpha * a.up[ibt, , drop = FALSE]
        }
      }
    }
    fa <- sum(obj)
    # cat("\nalpha-function:"); print(fa, digits = 16)
    if (fa <= fa.old + tol) break
    obj.old <- obj
    fa.old <- fa
  } # The alpha-function seems to be monotonically increasing for each iteration
  attr(a, "logmbeta") <- lgamma(rowSums(a)) - rowSums(lgamma(a))
  a
}

# EM algorithm for finite Dirichlet mixture
EM <- function(data, mix, tem = 5, tol = 1e-4, verbose = 0) {
  n <- nrow(data)
  cll <- tem > 5
  if (cll) {  # If tem <= 5, no need to calculate log-likelihood
    ll <- attr(mix, "ll")
    if (is.null(ll)) ll <- logLik(mix, data)
  }
  for (i in seq_len(tem)) {
    # if (verbose > 0 && cll) print_verbose(verbose, i - 1, ll, mix)
    logd <- ddirichlet(data, mix$alpha, log = TRUE, outer = TRUE) +
      rep(log(mix$pr), rep(n, length(mix$pr)))
    # log f(x_i;α_j) + log π_j = log [π_j * f(x_i;α_j)]
    ma <- lsei::matMaxs(logd)
    pi.d <- exp(logd - ma)
    pi.s <- rowSums(pi.d)
    p <- pi.d / pi.s
    # -ma or not will not change p, which is the same as
    # exp(logd) / rowSums(exp(logd))
    mix$pr <- colMeans(p) # π_j'
    j <- mix$pr < 1e-10
    if (any(j)) { # Delete support points with pr close to 0
      mix$pr <- mix$pr[!j]
      mix$pr <- mix$pr / sum(mix$pr)
      logmbeta <- attr(mix$alpha, "logmbeta")[!j]
      mix$alpha <- mix$alpha[!j, , drop = FALSE]
      attr(mix$alpha, "logmbeta") <- logmbeta
      p <- p[, !j, drop = FALSE]
    }
    mix$alpha <- iter_alpha(data, mix$alpha, p) # Update alpha
    if (cll) { # If tem <= 5, no need to check convergence
      ll.old <- ll
      ll <- logLik(mix, data)
      if (ll <= ll.old + tol) break
    }
  }
  # if (verbose > 0 && cll) print_verbose(verbose, i, ll, mix)
  if (!cll) ll <- logLik(mix, data)
  attr(mix, "mode") <- (mix$alpha - 1) * mix$h
  attr(mix, "ll") <- ll
  mix
}

# Generate initial alpha points of Modal EM by Multinomial mixture distribution,
# add the noise to make initial points more smooth and scatter
initial_MEM <- function(n, mix, data, mass) {
  D <- ncol(data)
  size <- 1 / mix$h
  pt <- rmultinomix(n, size, data, mass)
  noise <- runif(n * D, min = -0.5, max = 0.5)
  dim(noise) <- c(n, D)
  noise <- matrix(runif(length(pt), min = -0.5, max = 0.5), nrow = n)
  noise <- noise - rowMeans(noise)
  pt <- pt + noise
  # We want the generate x as the initial of alpha of Dirichlet mixture,
  # but x generate by Multinomial are integer. Thus add some noise of uniform
  # distribution to make x become real number
  pt[pt < 0] <- 0
  pt[pt > size] <- size # Pull back to boundary
  pt <- pt / rowSums(pt) * size + 1 # Ensure sum to size (a0)
  logmbeta <- lgamma(rowSums(pt)) - rowSums(lgamma(pt))
  ip <- rbind(mix$alpha, pt)
  attr(ip, "logmbeta") <- c(attr(mix$alpha, "logmbeta"), logmbeta)
  ip
}

# Modal EM algorithm of the gradient function of
# Dirichlet mixture's log-likelihood
# a: alpha; w: log f(x_i; G)
MEM <- function(data, a, w, plot = FALSE, tem = 5, tol = 1e-4, axl = 100) {
  n <- nrow(data)
  size <- sum(a[1, ]) # Calulate by h ??
  if (plot) {
    contour_gradient(data, w, size, axl = axl)
    points_simplex(a, size = size, add = TRUE,
                   pch = 20, col = "gray50", cex = 0.8)
  }
  # Plot contour of gradient of log-likelihood
  sg <- sg.old <- -Inf
  a.old <- a
  convergence <- FALSE
  for (i in seq_len(tem)) {
    logd <- ddirichlet(data, a, log = TRUE, outer = TRUE) - w
    # log f(x_i; a) - log f(x_i; G)
    g <- colSums(exp(logd)) # Sum_i f(x_i; a) / log f(x_i; G)
    sg <- sum(g) # Gradient value = sg - n (n is a constant, omit here)
    if (sg <= sg.old + tol) {
      convergence <- TRUE
      break
    }
    sg.old <- sg

    ma <- lsei::matMaxs(logd, dim = 2)
    pi.d <- exp(logd - rep(ma, each = n))
    pi.s <- colSums(pi.d)
    p <- pi.d / rep(pi.s, each = n)
    a <- iter_alpha(data, a, p) # Update alpha
    if (plot) segments_simplex(a.old, a, size, add = TRUE)
    # Plot the tack of alpha from initial to local maximum of gradient function
    a.old <- a
  }
  if (plot) points_simplex(a, size = size, col = "blue",
                           add = TRUE, pch = 19, cex = 0.8)
  # Final alpha after updating
  if (!convergence) {
    logd <- ddirichlet(data, a, log = TRUE, outer = TRUE) - w
    g <- colSums(exp(logd))
  }
  attr(a, "gradient") <- g - n # Gradient value
  a
}

# Select new support points with high gradient value
# t: points get from Modal EM; a: alpha
select_points <- function(t, a, by = 1) {
  a <- a[order(a[, by]), , drop = FALSE]
  mids <- (a[-1, by] + a[-nrow(a), by]) * 0.5
  # Middle points of each pair support points on one coordinate
  o <- order(t[, by])
  g <- attr(t, "gradient")[o]
  logmbeta <- attr(t, "logmbeta")[o]
  t <- t[o, , drop = FALSE] # order t by increasing
  index <- lsei::indx(t[, by], mids) + 1
  # Split t into several groups by interval of mids pairs
  jj <- aggregate(g, by = list(group = index), which.max)
  # Find the index of largest gradient value of each group,
  # these points should be new support points
  j <- match(jj$group, index) + jj$x - 1 # transform index
  np <- t[j, , drop = FALSE]
  gnp <- g[j]
  attr(np, "gradient") <- gnp
  attr(np, "logmbeta") <- logmbeta[j]
  np # New support points with their gradient values
}

# Update mass of components of Dirichlet mixture
# np: new support points; w: log f(x_i; G)
CNM <- function(np, mix, data, w, by = 1, maxbt = 50) {
  if (missing(w)) w <- ddirichletmix(data, mix, log = TRUE) # log f(x_i; G)
  ll <- attr(mix, "ll")
  if (is.null(ll)) ll <- logLik(mix, data)
  pt <- rbind(mix$alpha, np)
  attr(pt, "logmbeta") <- c(attr(mix$alpha, "logmbeta"), attr(np, "logmbeta"))
  mix.raw <-
    dirichletmix(pt, c(mix$pr, double(length(np))), by = by)
  # Add new support points and set their mass = 0
  logd <- ddirichlet(data, mix.raw$alpha, log = TRUE, outer = TRUE)
  # f(x_i; α_j)
  s <- exp(logd - w) # First-derivative of ll with respect to π
  g <- colSums(s)
  r <- lsei::pnnls(s, 2, sum = 1)
  # min||Sπ'-2||^2: a = S, b = 2 with π'1 = 1, π' >= 0
  sol <- r$x / sum(r$x) # Ensure sum to 1
  mix.up <- mix.raw
  mix.up$pr <- sol # Update π
  dea <- sum(g * (mix.up$pr - mix.raw$pr)) * 0.333 # 0.333 * first derivative
  mix <- mix.up
  alpha <- 1
  for (i in seq_len(maxbt)) { # Back tracking with Armijo rule
    llt <- logLik(mix, data)
    if (llt >= ll + dea * alpha) break
    alpha <- alpha * 0.5
    mix$pr <- (1 - alpha) * mix.raw$pr + alpha * mix.up$pr
  }
  attr(mix, "ll") <- llt
  mix
}

# Merge the similar components and drop tiny components (collapse_points)
# prec: threshold of merge and drop
unique.dirichletmix <- function(mix, prec = 0) {
  if (nrow(mix$alpha) == 1) return(mix)
  prec <- rep(prec, len = 2)
  # prec[1] for difference in each pair of support points, prec[2] for pr
  j <- mix$pr <= prec[2]
  a <- mix$alpha[!j, , drop = FALSE]
  pr <- mix$pr[!j] # Drop components with small pr
  if (prec[1] == 0) return(dirichletmix(a, pr))
  repeat {
    if (nrow(a) == 1) break
    d <- abs(a[-1, , drop = FALSE] - a[-nrow(a), , drop = FALSE])
    # Calculate difference of each pair of support points
    i <- lsei::matMaxs(d) <= prec[1] # Will use odd indexes of TRUE
    if (sum(i) == 0) break
    noti <- c(!i, TRUE) # not i
    even <- seq(2, by = 2, len = length(i))
    i[even] <- i[even] & noti[even - 1] & noti[even + 1]
    # Use evens if only isolated
    ind <- which(i) # Indexes of points need to be merged
    pr2 <- pr[ind] + pr[ind + 1]
    a[ind, ] <- (a[ind, , drop = FALSE] * pr[ind] +
                   a[ind + 1, , drop = FALSE] * pr[ind + 1]) / pr2
    # Collapse neighboring pairs to the first
    pr[ind] <- pr2
    a <- a[-(ind + 1), , drop = FALSE]
    pr <- pr[-(ind + 1)] # Remove the second of a pair
  }
  dirichletmix(a, pr)
}

# Drop and merge the components of Dirichlet mixture by compare log-likelihood
collapse_points <- function(mix, data, tol = 1e-3) {
  ll <- attr(mix, "ll")
  if (is.null(ll)) ll <- logLik(mix, data)
  j <- mix$pr == 0 # Drop components with mass = 0
  if (any(j)) {
    logmbeta <- attr(mix$alpha, "logmbeta")[!j]
    mix$alpha <- mix$alpha[!j, , drop = FALSE]
    mix$pr <- mix$pr[!j]
    attr(mix$alpha, "logmbeta") <- logmbeta
    attr(mix, "mode") <- attr(mix, "mode")[!j, , drop = FALSE]
    ll <- logLik(mix, data) # Update log-likelihood after drop components
  }
  m <- length(mix$pr)
  if (m == 1 || tol <= 0) {
    attr(mix, "ll") <- ll
    return(mix)
  }
  o <- do.call(order, split(mix$alpha, col(mix$alpha))) # Order by each column
  logmbeta <- attr(mix$alpha, "logmbeta")[o]
  mix$alpha <- mix$alpha[o, , drop = FALSE]
  mix$pr <- mix$pr[o]
  attr(mix$alpha, "logmbeta") <- logmbeta
  prec <- 1e-10 # Initial precision value
  repeat {
    if (nrow(mix$alpha) <= 1) break
    d <- abs(mix$alpha[-1, , drop = FALSE] -
               mix$alpha[-nrow(mix$alpha), , drop = FALSE])
    prec <- 10 * max(min(lsei::matMaxs(d)), prec)
    mixt <- unique(mix, prec = c(prec, 0))
    # Merge all pairs of components with difference less than
    # the current smallest difference * 10
    attr(mixt, "ll") <- logLik(mixt, data)
    if (attr(mixt, "ll") >= ll - tol) mix <- mixt
    else break
  }
  mix
}

# Overall algorithm under a fixed bandwidth h
cnmdm <- function(data, h,
                  verbose = 0, by = NULL, maxit = 200, nt = 20, tol = 1e-4) {
  if (!is.numeric(data)) data <- as.matrix(data)
  attr(data, "log") <- pmax(log(data), log(1e-100))
  if (is.null(by)) by <- which.min(colSums(apply(data, 2, duplicated)))
  # Initial the mixture
  mix <- initial_dirichletmix(data, h, by)
  convergence <- 1
  tem <- 5
  for (i in seq_len(maxit)) {
    ll.old <- attr(mix, "ll")
    if (verbose > 0) print_verbose(verbose, i - 1, mix)
    # EM algorithm
    mix <- EM(data, mix, tem = tem)
    # no need combine support points in EM step
    # Expand new support points
    w <- ddirichletmix(data, mix, log = TRUE) # log f(x_i;G)
    t <- initial_MEM(nt, mix, data, exp(-w))
    t <- MEM(data, t, w)
    np <- select_points(t, mix$alpha, by)
    # points_simplex(np, size = sum(np[1, ]), add = TRUE,
    # col = "blue", pch = 19, cex = 0.65)
    # Update mass of components
    mix <- CNM(np, mix, data, w, by)
    # Collapse support points
    mix <- collapse_points(mix, data)
    ll <- attr(mix, "ll")

    # If log-likelihood is converge, then execute the algorithm one last time
    # If log-likelihood is still converge after last time then break
    if (ll <= ll.old + tol) {
      if (tem == 100) {
        convergence <- 0
        break
      }
      if (tem == 5) tem <- 100
    } else {
      tem <- 5
    }
  }
  if (verbose > 0) print_verbose(verbose, i, mix)
  attr(mix, "ll") <- attr(mix$alpha, "logmbeta") <- NULL
  attr(mix, "mode") <- (mix$alpha - 1) * mix$h
  structure(list(num.iterations = i, convergence = convergence,
                 maxgrad = max(attr(np, "gradient")), ll = ll, mix = mix),
            class = "cnmdm")
}
