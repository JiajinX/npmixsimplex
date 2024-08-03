## ======================= Step function of algorithm ======================= ##

# Initialize Dirichlet mixture by grouping data (simple clustering)
# h: bandwidth; by: index of column used to sort components
initial_dirichletmix <- function(data, h, by) {
  n <- nrow(data)
  D <- ncol(data)
  a0 <- 1 / h + D
  tdata <- (asin(sqrt((data / h + 1) / a0)) - asin(sqrt(1 / a0)))
  # Variance stabilization transformation
  side <- (asin(sqrt((1 / h + 1) / a0)) - asin(sqrt(1 / a0)))
  tdata <- tdata / side
  # Normalize the tdata into a same size space (independent of h)
  tmu <- colMeans(tdata) # Or transform(colMeans(data)) ???
  # tmu should be the first center of cluster groups
  d2 <- d2mu <- colSums((t(tdata) - tmu)^2) # Distance from tdata to tmu
  d <- sqrt(d2)
  index <- 1:n
  pts <- matrix(nrow = n, ncol = D)
  for (k in seq_len(n)) {
    i <- d2 > h
    # All points with distance < h will be in a cluster group with
    # pts[k, ] as the center and h as the radius of circle
    index <- index[i] # Uncovered data points
    if (length(index) == 0) break
    imin <- index[which.min(d2mu[index])]
    update <- (tdata[imin, ] - tmu) * sqrt(1 + 4 / 9 * h  / d[imin]) # 4/9 ???
    pts[k, ] <- tmu + update
    # Find the closet uncovered point to the tmu, expand this point
    # outward a certain distance as the new cluster center
    if (any(sum(pts[k, ] < 0))) { # If points outside the boundary, pull it back
      lf <- tmu / -update
      lf <- min(lf[lf > 0])
      pts[k, ] <- tmu + lf * update
    }
    d2 <- colSums((t(tdata)[, index, drop = FALSE] - pts[k, ])^2)
  }
  tmu <- if (k == 1) tmu else rbind(tmu, pts[1:(k - 1), , drop = FALSE])
  # Set tmu as center of cluster
  rownames(tmu) <- NULL
  d2mu <- rowSums((data[rep(1:n, k), ] - tmu[rep(1:k, each = n), ])^2)
  dim(d2mu) <- c(n, k)
  imu <- apply(d2mu, 1, which.min)
  nj <- tabulate(imu, nbins = k)
  pr <- nj / sum(nj) # Calculate pr by distance of data points and tmu
  j <- pr != 0
  emu <- sin(tmu[j, , drop = FALSE] * side + asin(sqrt(1 / a0)))^2
  emu <- emu / rowSums(emu) # emu should > 1 / a0
  mode <- (emu * a0 - 1) * h
  # Transform tmu back to the simplex, which are the means (a_i / a0) of the
  # Dirichlet components, then got the modes by means
  i <- which(rowSums(mode < 0) > 0)
  if (!identical(i, integer(0))) { # If mode outside the boundary (mode < 0)
    smode <- mode[i, , drop = FALSE]
    smode[smode < 0] <- 0
    mode[i, ] <- smode / rowSums(smode)
  }
  alpha <- mode / h + 1
  dirichletmix(alpha, pr[j], by = by) # Drop tiny components ???
}

# Objective function needs to be maximized by iteration alpha (EM -> iter_alpha)
# Also with the first and second derivative of objective function
# x: data; a: alpha; p: value get from EM
calfun_alpha <- function(x, a, p, deriv = 0) {
  n <- nrow(x)
  D <- ncol(x)
  m <- nrow(a)
  logx <- attr(x, "log")
  if (is.null(logx)) logx <- log(x)
  xp <-
    switch(deriv + 1,
      rowSums((a %*% t(logx) - rowSums(lgamma(a))) * t(p)),
      # Original objective function
      (logx[rep(1:n, m), ] - digamma(a)[rep(1:m, each = n), ]) * c(p),
      # p_ij * [log(x_ik) - ψ(α_jk)] - first derivative
      -trigamma(a)[rep(1:m, each = n), ] * c(p)
      # -p_ij * ψ'(α_jk) - second derivative
    )
  if (deriv != 0) {
    dim(xp) <- c(n, m * D)
    xp <- matrix(colSums(xp), nrow = m)
  }
  xp
}

# Check if alpha points out of boundary after updating and do back tracking
# (EM -> iter_alpha)
stepup <- function(a, update) {
  a.new <- a + update
  # cat("\nalpha:\n"); print(a, digits = 16)
  # cat("\nupdate:\n"); print(update, digits = 16)
  j <- which(rowSums(a.new < 1) > 0) # Find components with alpha' < 1
  if (!identical(j, integer(0))) { # Pull alpha back to 1
    sa <- a[j, , drop = FALSE]
    supdate <- update[j, , drop = FALSE]
    lf <- (sa - 1) / -supdate
    # For a_k = 1, update_k will be >= 0:
    # when a_k = 1 & update_k = 0, lf = NAN
    # when a_k = 1 & update_k > 0, lf = 0
    lf[lf < 0 | is.na(lf)] <- Inf
    # lf < 0: a_k > 1 & update_k > 0
    # lf = NA: update_k = 0
    rmin <- apply(lf, 1, function(r) which.min(r))
    minloc <- cbind(seq_along(j), j, rmin)
    a.new[j, ] <- sa + lf[minloc[, -2, drop = FALSE]] * supdate
    # Back tracking, pull back to boundary
    a.new[minloc[, -1, drop = FALSE]] <- 1
  }
  # cat("\na.new:\n"); print(a.new, digits = 16)
  a.new
}

# Updating alpha by using Newton-Raphson method (EM)
# x: data; a: alpha; p: value get from EM
iter_alpha <- function(x, a, p, maxit = 100, tol = 1e-4) {
  attr(x, "log") <- log(x)
  fa.old <- fa <- sum(calfun_alpha(x, a, p))
  # cat("\nalpha-function:"); print(as.character(signif(fa, 16)))
  # cat("\nalpha:\n"); print(a, digits = 16)
  for (i in seq_len(maxit)) {
    isbound <- a == 1
    # TRUE for "1_set" (on the boundary), FALSE for ">1_set"
    g <- calfun_alpha(x, a, p, deriv = 1)
    H <- calfun_alpha(x, a, p, deriv = 2)
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
    a <- stepup(a, update) # Update alpha
    fa <- sum(calfun_alpha(x, a, p)) # ??? bug: fa -> -Inf
    # Decreasing from -15151.1323580999 to -Inf
    if (fa < fa.old) cat("Decreasing from", as.character(signif(fa.old, 16)),
                         "to", as.character(signif(fa, 16)), "\n")
    # cat("\nfa.old:"); print(as.character(signif(fa.old, 16)))
    # cat("\nfa:"); print(as.character(signif(fa, 16)))
    if (fa <= fa.old + tol) break
    fa.old <- fa
  } # The alpha-function seems to be monotonically increasing for each iteration
  a
}

# EM algorithm of Dirichlet mixture
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
      mix$alpha <- mix$alpha[!j, ]
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

# Generate initial alpha points of Modal EM by Multinomial mixture distribution
# n: number of points want to generate; size: sum of x (= a0);
# pr: parameters of multinomial; mass: mass of Multinomial mixture components
gmultinomix <- function(n, size, pr, mass) {
  if (is.vector(pr)) pr <- t(pr)
  m <- nrow(pr)
  D <- ncol(pr)
  mass <- rep(mass, len = m)
  mass <- mass / sum(mass)
  insz <- ceiling(size) # Must be an integer
  tb <- table(sample.int(m, n, prob = mass, replace = TRUE))
  cl <- c(0, cumsum(tb)) # Index
  point <- matrix(nrow = D, ncol = n)
  for (i in seq_along(tb)) {
    point[, (cl[i] + 1):cl[i + 1]] <-
      rmultinom(tb[i], insz, pr[as.integer(names(tb[i])), ])
    # Generate x from Multinomial mixture distribution
  }
  point <- t(point)
  noise <- matrix(runif(length(point), min = -0.5, max = 0.5), nrow = n)
  noise <- noise - rowMeans(noise)
  point <- point + noise
  # We want the generate x as the initial of alpha of Dirichlet mixture,
  # but x generate by Multinomial are integer. Thus add some noise of uniform
  # distribution to make x become real number
  point[point < 0] <- 0
  point[point > size] <- size # Pull back to boundary
  point / rowSums(point) * size # Ensure sum to size (a0)
}

# Modal EM algorithm of the gradient function of
# Dirichlet mixture's log-likelihood
# a: alpha; w: log f(x_i; G)
MEM <- function(data, a, w, plot = FALSE, tem = 5, tol = 1e-4) {
  n <- nrow(data)
  size <- sum(a[1, ])
  if (plot) contour_gradient(data, w, size)
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
    if (plot) segments.simplex(a.old, a, size, add = TRUE)
    # Plot the tack of alpha from initial to local maximum of gradient function
    a.old <- a
  }
  if (plot) points.simplex(a, size, add = TRUE, pch = 19, cex = 0.65)
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
  np # New support points with their gradient values
}

# Update mass of components of Dirichlet mixture
# np: new support points; w: log f(x_i; G)
CNM <- function(np, mix, data, w, by = 1, maxbt = 50) {
  if (missing(w)) w <- ddirichletmix(data, mix, log = TRUE) # log f(x_i; G)
  ll <- attr(mix, "ll")
  if (is.null(ll)) ll <- logLik(mix, data)
  mix.raw <- dirichletmix(rbind(mix$alpha, np), c(mix$pr, double(length(np))),
                          by = by)
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
    mix$alpha <- mix$alpha[!j, , drop = FALSE]
    mix$pr <- mix$pr[!j]
    ll <- logLik(mix, data) # Update log-likelihood after drop components
  }
  m <- length(mix$pr)
  if (m == 1 || tol <= 0) return(list("mix" = mix, "ll" = ll))
  o <- do.call(order, split(mix$alpha, col(mix$alpha))) # Order by each column
  mix$alpha <- mix$alpha[o, , drop = FALSE]
  mix$pr <- mix$pr[o]
  prec <- 1e-10 # Initial precision value
  repeat {
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

# Overall algorithm
cnmdm <- function(data, h, verbose = 0, maxit = 100, nt = 20, tol = 1e-4) {
  if (!is.numeric(data)) data <- as.matrix(data)
  # if (missing(by))
  by <- which.min(colSums(apply(data, 2, duplicated)))
  # Initial the mixture
  mix <- initial_dirichletmix(data, h, by)
  attr(mix, "ll") <- logLik(mix, data)
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
    t <- rbind(mix$alpha, gmultinomix(nt, 1 / mix$h, data, 1 / exp(w)) + 1)
    t <- MEM(data, t, w)
    np <- select_points(t, mix$alpha, by)
    # points.simplex(np, size = sum(np[1, ]), add = TRUE,
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
  structure(list(num.iterations = i, convergence = convergence,
                 max.gradient = max(attr(np, "gradient")), ll = ll, mix = mix),
            class = "cnmdm")
}
