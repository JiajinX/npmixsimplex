# Initial h can get from the Dirichlet distribution fit to the data
initial_a0 <- function(x, maxit = 20, tol = 1e-3) {
  logx <- attr(x, "log")
  if (is.null(logx)) logx <- pmax(log(x), log(1e-100))
  elog <- colMeans(logx)
  D <- length(elog)
  a <- rep(2, D)
  for (i in 1:maxit) {
    logC <- digamma(sum(a)) - digamma(a) # Constant part in log
    g <- logC + elog
    negH <- diag(trigamma(a)) - trigamma(sum(a)) # -H
    update <- solve(negH, g)
    at <- a + update
    if (any(at < 1)) {
      lf <- (a - 1) / -update
      lf <- lf[lf > 0] # Find all update < 0 and a < 1 after updating
      at <- a + min(lf) / 2 * update
    }
    if (abs(sum(at) - sum(a)) < tol) break
    a <- at
  }
  sum(at)
}

dmde <- function(data, bwf = 10:1 / 10, tol = 1e-4, ictol = 10,
                 criterion = c("AIC","AICc","BIC"),
                 exhaustive = FALSE, verbose = 0, by,
                 plot = FALSE) {
  criterion <- match.arg(criterion)
  bwf <- unique(sort(bwf, decreasing = TRUE))
  if (missing(by)) by <- which.min(colSums(apply(data, 2, duplicated)))
  k <- length(bwf)
  ic <- double(k)
  icmin <- Inf
  names(icmin) <- criterion
  n <- nrow(data) # sample size
  D <- ncol(data) # dimension
  d <- D - 1
  volf <- bwf^(1/d)
  ia0 <- initial_a0(data) # Initial a0
  info <- matrix(nrow = k, ncol = 4,
                 dimnames = list(1:k, c("bwf", "h", "ll", criterion)))
  for (i in seq_len(k)) {
    h <- 1 / ((ia0 + 1) / volf[i]^2 - 1 - D) # a0 + 1 \propto 1 / sd^2
    r <- cnmdm(data, h = h, by = by, tol = tol)
    p <- D * nrow(r$mix$alpha) - 1 # (Dim - 1) * n_Comp + n_Comp - 1
    ic[i] <- -2 * r$ll +
      switch(criterion,
             AICc = {2 * p + 2 * p * (p + 1) / max(n - p - 1, 0)},
             AIC = {2 * p},
             BIC = {log(n) * p}
             )
    info[i, ] <- c(bwf[i], h, r$ll, ic[i])
    if (k == 1 || ic[i] < icmin) {
      icmin <- ic[i]
      imin <- i
      mix <- r$mix
      ll <- r$ll
      maxgrad <- r$max.gradient
    }
    if (plot) pairs(r$mix, data)
    if (verbose >= 1) {
      cat("Iter.#", i,
          ":\tbwf = ", sprintf("%.3f", bwf[i]),
          "\th = ", sprintf("%.3f", h),
          "\tll = ", r$ll, "\t", criterion, " = ", ic[i], "\n", sep = "")
    }
    if (verbose >= 2) print(r$mix)
    if (ic[i] == Inf || !exhaustive && ic[i] >= icmin + ictol) break
  }

  fit <- list(info = info[1:i, ], bwf = bwf[imin],
              maxgrad = maxgrad, ll = ll, criterion = icmin, mix = mix)
  class(fit) = "dmde"
  fit
}

# Calculate CVISE and CVKLD for different bandwidth
CVfold <- function(train, test, bw, nsamp = 1e5,
                   tol = 1e-4, by = NULL, verbose) {
  bw <- unique(sort(bw, decreasing = TRUE))
  if (is.null(by)) by <- which.min(colSums(apply(train, 2, duplicated)))
  k <- length(bw)
  ise <- kl <- double(k)
  for (i in seq_len(k)) { # For each bandwidth
    if (verbose > 0) {
      if (i == k) cat(i, "\n") else cat(i, " ")
    }
    r <- cnmdm(train, h = bw[i], by = by, tol = tol)$mix
    # Dirichlet mixture estimator
    samp <- genemix(nsamp, r) # Generate 1e5 sample
    term1 <- predmix(samp, r)
    term2 <- predmix(test, r)
    # ise[i] <- mean(attr(term1, "euc")) - 2 * mean(attr(term2, "euc"))
    ise[i] <- mean(term1) - 2 * mean(term2)
    kl[i] <- -mean(log(term2[term2 > 0]))
    # Avoid zero density when test on the boundary but support points not
  }
  list(ise = ise, kl = kl)
}

# Model selection with CVISE and CVKLD
dmdeCV <- function(data, k = 5, bwf = 10:1 / 10, verbose = 0) {
  D <- ncol(data)
  d <- D - 1
  ia0 <- initial_a0(data) # Initial of a0
  bw <- 1 / ((ia0 + 1) / bwf^(2/d) - 1 - D)
  # Bandwidth candidates by initial a0 and bandwidth factor
  fold <- sample(rep(1:k, length.out = nrow(data))) # Split data into k-folds
  cvise <- matrix(nrow = k, ncol = length(bwf))
  colnames(cvise) <- bwf
  cvkl <- cvise # k * length(bw)
  for (i in seq_len(k)) {
    if (verbose > 0) cat("# K-fold:", i, "\n")
    train <- data[fold != i, ]
    test <- data[fold == i, ]
    ms <- CVfold(train, test, bw = bw, verbose = verbose)
    cvise[i, ] <- ms$ise
    cvkl[i, ] <- ms$kl
  }
  # Select the bandwidth with the best overall performance in CV
  h_ise <- bw[which.min(colMeans(cvise))]
  h_kl <- bw[which.min(colMeans(cvkl))]
  if (h_ise == h_kl) {
    mix_ise <- mix_kl <- cnmdm(data, h = h_ise)
  } else {
    mix_ise <- cnmdm(data, h = h_ise)
    mix_kl <- cnmdm(data, h = h_kl)
  }
  re <- list(mix_ise = mix_ise, mix_kl = mix_kl)
  names(bw) <- round(bwf, 2)
  attr(re, "bandwidth") <- bw
  attr(re, "isemtx") <- cvise
  attr(re, "klmtx") <- cvkl
  re
}
