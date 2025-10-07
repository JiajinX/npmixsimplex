# Create proposal density q as a smoothed version of the true density
# and generate large samples of q
IE_samp <- function(mix0, n = 1e5, fct = 1.2) {
  clst <- class(mix0)
  f <- mix0 # f is similar to mix0 with heavier tail
  if (clst == "kde") {
    diag(f$H) <- diag(f$H) * fct
    f <- ks::kde(f$x, H = f$H)
  } else if (clst == "mvnmix") {
    sig <- t(f$Sigma) %*% f$Sigma
    diag(sig) <- diag(sig) * fct
    f$Sigma <- chol(sig)
  } else if (clst == "mvnormix") {
    m <- length(mix0$pr)
    D <- ncol(mix0$Sigma)
    spSig <- split(mix0$Sigma, rep(1:m, each = D)) # Split sigma
    heavySig <- lapply(spSig, function(x) {
      dim(x) <- c(D, D)
      diag(x) <- diag(x) * fct # ??? + or * for heavier tail
      x
    })
    f$Sigma <- do.call(rbind, heavySig)
  } else { # Dirichlet mixture
    mode <- attr(mix0, "mode")
    h <- fct / (rowSums(mix0$alpha) - ncol(mix0$alpha)) # Larger bandwidth
    f$alpha <- mode / h + 1
  }
  samp <- genemix(n, f)
  d <- predmix(samp, mix0)
  attr(f, "samp") <- samp[d > 1e-10, ] # Delete points with too low density
  # attr(f, "samp") <- genemix(n, f) # Sampling from heavier tail mixture
  attr(f, "dsamp") <- predmix(attr(f, "samp"), f)
  attr(mix0, "dsamp") <- predmix(attr(f, "samp"), mix0)
  list(mix0 = mix0, mix0_ht = f) # Density of mix0 and heavier tail mixture
}

# Calculate ISE and IAE
IE <- function(mix, est, n = 1e5) {
  if (is.null(mix$mix0_ht)) mix <- IE_samp(mix, n = n)
  mix0 <- mix$mix0
  mix0_ht <- mix$mix0_ht
  data <- attr(mix0_ht, "samp")
  td <- attr(mix0, "dsamp")
  hd <- attr(mix0_ht, "dsamp") # Density of heavier tail
  ed <- predmix(data, est)
  euctd <- attr(td, "euc")
  euced <- attr(ed, "euc")
  euchd <- attr(hd, "euc")

  iae <- mean(abs(td - ed) / hd)
  iseR <- mean((euctd - euced)^2 / euchd)
  iseS <- mean((td - ed)^2 / hd)
  list(iae = iae, iseR = iseR, iseS = iseS)
}

# Generate large samples from the true density
KLD_samp <- function(mix0, n = 1e5) {
  attr(mix0, "samp") <- genemix(n, mix0)
  logd <- predmix(attr(mix0, "samp"), mix0, log = TRUE) # Should return log
  mlog <- mean(logd)
  attr(mix0, "mlog") <- mlog
  mix0
}

# Calculate KLD
KLD <- function(mix0, est, n = 1e5) {
  if (is.null(attr(mix0, "samp"))) mix0 <- KLD_samp(mix0, n = n)
  data <- attr(mix0, "samp")
  tmlog <- attr(mix0, "mlog")
  loged <- predmix(data, est, log = TRUE)
  emlog <- mean(loged)
  tmlog - emlog
}
