# Create multivariate normal distribution from Chacon'S parameters
# D: dimension of the mixture
chaconmix <- function(no, D = 2) {
  if (D < 2) stop("Dimension setting error.")
  mix <- switch(no,
                # 1.
                list(mu = rbind(c(0, 0)),
                     Sigma = ks::invvech(c(1 / 4, 0, 1)),
                     pr = 1),
                # 2.
                list(mu = rbind(c(0, 0)),
                     Sigma = ks::invvech(c(1, 9 / 10, 1)),
                     pr = 1),
                # 3.
                list(mu = matrix(rep(c(3, -3), 8) *
                                   rep((1 - (4 / 5)^(0:7)), each = 2),
                                 byrow = TRUE, nrow = 8),
                     Sigma = matrix(rep(ks::invvech(c(1, -9 / 10, 1)), 8) *
                                      rep((4 / 5)^(2 * (0:7)), each = 4),
                                    byrow = TRUE, nrow = 8 * 2),
                     pr = rep(1, 8) / 8),
                # 4.
                list(mu = rbind(c(-3 / 2, 3 / 2), c(3 / 2, -3 / 2), c(0, 0)),
                     Sigma = rbind(diag(9 / 16, 2), diag(9 / 16, 2),
                                   ks::invvech(c(4 / 5, -18 / 25, 4 / 5)) * 9 / 16),
                     pr = c(4, 4, 3) / 11),
                # 5.
                list(mu = rbind(c(0, 0), c(0, 0)),
                     Sigma = rbind(ks::invvech(c(1, 1, 4)),
                                   ks::invvech(c(4 / 9, -1 / 3, 4 / 9))),
                     pr = c(2, 1) / 3),
                # 6.
                list(mu = rbind(c(1, 0), c(-1, 0)),
                     Sigma = rbind(diag(4 / 9, 2), diag(4 / 9, 2)),
                     pr = c(1, 1) / 2),
                # 7.
                list(mu = rbind(c(-1, 1), c(1, -1)),
                     Sigma = rbind(ks::invvech(c(4 / 9, 4 / 15, 4 / 9)),
                                   ks::invvech(c(4 / 9, 4 / 15, 4 / 9))),
                     pr = c(1, 1) / 2),
                # 8.
                list(mu = rbind(c(1, -1), c(-1, 1)),
                     Sigma = rbind(ks::invvech(c(4 / 9, 14 / 45, 4 / 9)), diag(4 / 9, 2)),
                     pr = c(1, 1) / 2),
                # 9.
                list(mu = rbind(c(-1, 0), c(1, 2 / sqrt(3)), c(1, -2 / sqrt(3))),
                     Sigma = rbind(ks::invvech(c(9, 63 / 10, 49 / 4)),
                                   ks::invvech(c(9, 0, 49 / 4)),
                                   ks::invvech(c(9, 0, 49 / 4))) / 25,
                     pr = c(3, 3, 1) / 7),
                # 10.
                list(mu = matrix(c(rep(0, 4), c(-1, -1, 1, -1, -1, 1, 1, 1)),
                                 byrow = TRUE, nrow = 6),
                     Sigma = rbind(diag(2),
                                   matrix(rep(diag(1 / 16, 2), 5),
                                          byrow = TRUE, nrow = 5 * 2)),
                     pr = c(1 / 2, rep(1 / 10, 5))),
                # 11.
                list(mu = rbind(c(-3 / 2, 0), c(3 / 2, 0), cbind(-1:1 - 3 / 2, -1:1),
                                c(0, 0), cbind(-1:1 + 3 / 2, -1:1)),
                     Sigma = matrix(c(rep(ks::invvech(c(4 / 9, 4 / 15, 4 / 9)), 2),
                                      rep(ks::invvech(c(1 / 15, 1 / 25, 1 / 15)), 3) / 15,
                                      ks::invvech(c(1, 3 / 5, 1)) / 9,
                                      rep(ks::invvech(c(1 / 15, 1 / 25, 1 / 15)), 3) / 15),
                                    byrow = TRUE, nrow = 9 * 2),
                     pr = c(12 / 25, 12 / 25, rep(1 / 350, 3), 8 / 350, rep(1 / 350, 3))),
                # 12.
                list(mu = matrix(c(rep(0, 4), c(1, 1, -1, 1, -1, -1, 1, -1)),
                                 byrow = TRUE, nrow = 6),
                     Sigma = rbind(diag(2),
                                   ks::invvech(c(1, -9 / 10, 1)) / 16,
                                   ks::invvech(c(1, -9 / 10, 1)) / 4,
                                   diag(1 / 8, 2),
                                   ks::invvech(c(1, -9 / 10, 1)) / 8,
                                   diag(1 / 16, 2)),
                     pr = c(1 / 2, 3 / 40, 1 / 5, rep(3 / 40, 3)))
  )
  if (D > 2) {
    m <- length(mix$pr)
    newcol <- matrix(0, nrow = m, ncol = D - 2)
    mix$mu <- cbind(mix$mu, newcol) # Add 0 for mu of each component
    spSig <- split(mix$Sigma, rep(1:m, each = 2))
    exspSig <- lapply(spSig, function(x) {
      dim(x) <- c(2, 2)
      eigenvalue <- eigen(x)$values
      neweig <- seq(min(eigenvalue), max(eigenvalue), len = D)[-c(1, D)]
      # Add eigen value by sequence
      exmtrix <- matrix(0, nrow = D, ncol = D)
      exmtrix[1:2, 1:2] <- x # Original Sigma
      diag(exmtrix)[3:D] <- neweig
      exmtrix
    })
    mix$Sigma <- do.call(rbind, exspSig)
  }
  class(mix) <- "mvnormix"
  mix
}

Dir_mix <- function(no, D = 3) {
  if (D < 3) stop("Dimension setting error.")
  p <- switch(no,
              # Same bandwidth
              # 1.
              list(mu = rep(1, 3) / 3, pr = 1, h = 0.1), # Centered
              # 2.
              list(mu = c(0, 0.3, 0.7),  pr = 1, h = 0.25), # Boundary
              # 3.
              list(mu = c(0, 0, 1),  pr = 1, h = 0.2), # Cornered
              # 4.
              list(mu = rbind(c(0.5, 0.3, 0.2), c(0.2, 0.3, 0.5)),
                   pr = rep(1, 2) / 2, h = 0.05), # Bimodal
              # 5.
              list(mu = rbind(c(0.02, 0.8, 0.18), c(0.65, 0, 0.35), c(0, 0, 1)),
                   pr = c(1, 5, 12) / 18, h = 0.15), # Separated Trimodal
              # 6.
              {
                i <- 1:6
                l <- length(i)
                x1 <- i / l
                r <- 1 - x1
                mu <- cbind(x1, r * x1, r^2)
                colnames(mu) <- NULL
                list(mu = mu, pr = rep(1, l) / l, h = 0.2)
              },

              # Different bandwidth
              # 7.
              list(mu = rbind(c(0.8, 0, 0.2), c(0.6, 0.4, 0)),
                   pr = c(0.3, 0.7), h = c(0.1, 0.2)),
              # 8.
              list(mu = rbind(c(1/3, 1/3 + 0.08, 1/3 - 0.08),
                              c(1/3, 0.1, 0.9 - 1/3),
                              c(1/3 - 0.2, 1/3 + 0.18, 1/3 + 0.02),
                              c(0.44, 0.42, 0.14)),
                   pr = c(3, 7, 4, 1) / 15, h = c(0.1, 0.05, 0.035, 0.04)),
              # 9.
              list(mu = rbind(c(1, 0, 0), c(0.01, 0.98, 0.01), c(0, 0.01, 0.99),
                              c(0, 0.5, 0.5), c(0.3, 0, 0.7), c(0.6, 0.4, 0)),
                   pr = c(rep(0.4 / 3, 3), rep(0.6 / 3, 3)),
                   h = c(0.2, 0.1, 0.06, 0.06, 0.2, 0.1)),
              # 10.
              {
                i <- 1:4
                l <- length(i)
                x1 <- i / l * 0.6
                xr <- 1 - x1
                xmu <- cbind(x1, xr * 1 / 8, xr * 7 / 8)
                y1 <- i / l * 0.3
                yr <- 1 - y1
                ymu <- cbind(y1, yr * 9 / 10, yr * 1 / 10)
                colnames(xmu) <- colnames(ymu) <- NULL
                list(mu = rbind(xmu, ymu), pr = rep(1, 2 * l) / (2 * l),
                     h = c(0.5 * rev(x1), rev(y1)))
              },
              # 11.
              {
                i <- 1:6
                l <- length(i)
                x1 <- i / (l + 1.5)
                xr <- 1 - x1
                xmu <- cbind(x1, xr * x1, xr^2)
                ymu <- cbind(xmu[, 2] + (1 / 7.5)^2,
                             xmu[, 1] - (1 / 7.5)^2,
                             xmu[, 3])
                colnames(xmu) <- colnames(ymu) <- NULL
                list(mu = rbind(xmu, ymu),
                     pr = c(rep(0.45 / l, l), rep(0.55 / l, l)),
                     h = c(sqrt(x1) * 0.2, sqrt(rev(x1)) * 0.15))
              },
              # 12.
              {
                lx <- 3
                x3 <- seq(1/3, 1, len = lx)
                xmu <- cbind(1 - x3, x3, 0)
                i <- 1:4
                li <- length(i)
                y1 <- i / li
                ry <- 1 - y1
                ymu <- cbind(y1, ry * y1, ry^2)
                z <- ymu[1, ] / 2
                zmu <- c(z[1], 1 - z[1] - z[3], z[3])
                colnames(xmu) <- colnames(ymu) <- rownames(zmu) <- NULL
                list(mu = rbind(xmu, ymu, zmu),
                     pr = c(rep(0.25 / lx, lx), rep(0.6 / li, li), 0.15),
                     h = c(seq(13, by = -3, len = lx),
                           seq(10, by = -1.6, len = li), 6) / 100)
              }

  )
  if (D > 3) {
    if (is.vector(p$mu)) p$mu <- t(p$mu)
    p$mu <- p$mu * 3 / D
    newcol <- matrix((1 - 3 / D) / (D - 3), nrow = nrow(p$mu), ncol = (D - 3))
    p$mu <- cbind(p$mu, newcol)
  }
  dirichletmix(a = p$mu / p$h + 1, pr = p$pr)
}

plotestmix <- function(type = c("chacon", "dirichlet"), no = 1:12,
                       prop = FALSE) {
  type <- match.arg(type)
  n <- length(no)
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(ceiling(n / 4), 4))
  for (i in no) {
    mix0 <- if(type == "chacon") chaconmix(i) else Dir_mix(i)
    if (prop) data <- genemix(1000, mix0)
    ifelse(type == "chacon",
           contour_r2s(mix0, data, prop = prop),
           plot(mix0, data, prop = prop))
  }
  par(old_par)
}

# ms_model: function for simulation
# type: type of simulation models: chacon or dirichlet
# no: simulation No.: 1~12
# D: dimension; n: generated sample size
# samp: Monte Carlo sample size
# rep: the number of experimental repetitions
# Output: loss value (MIAE, MISE on Euclidean space, KLD), running time
ms_model <- function(type = c("dirichlet", "chacon"),
                     no, D = 3, n = 100, samp = 1e4, rep = 1) {
  type <- match.arg(type)
  cn <- c("no", "ilr-PI", "ilr-CV", "ilr-MDE",
          "Dir-KDE", "Dir-AIC", "Dir-CVISE", "Dir-CVKLD") # test estimators
  nest <- length(cn)
  miae <- matrix(nrow = rep, ncol = nest)
  colnames(miae) <- cn
  tm <- kld <- miseR <- miae
  tm <- tm[, -nest, drop = FALSE]
  colnames(tm)[nest - 1] <- "Dir-CV"
  mix0 <- switch(type,
                 chacon = chaconmix(no, D = D - 1),
                 dirichlet = Dir_mix(no, D = D))
  mix0_ie <- IE_samp(mix0, n = samp)
  mix0_kld <- KLD_samp(mix0, n = samp)
  for (i in 1:rep) {
    print(i)
    data <- genemix(n, mix0)
    idata <- attr(data, "ilr")

    t1 <- Sys.time()
    f_plug <- ks::kde(idata) # KDE-ilr: Plug_in
    t2 <- Sys.time()
    f_cv <- ks::kde(idata, H = ks::Hlscv(idata)) # KDE-ilr: CV
    t3 <- Sys.time()
    f_imde <- npden::mmde(idata, 10:1 / 10)$mix # MDE-ilr
    t4 <- Sys.time()

    f_kdir <- dirichletmix(a = data / hkdeDir(data) + 1) # KDE-Dir
    t5 <- Sys.time()
    f_dmde <- dmde(data, ictol = 40)$mix
    t6 <- Sys.time()
    recv <- dmdeCV(data)
    t7 <- Sys.time()
    f_cvise <- recv$mix_ise$mix
    f_cvkl <- recv$mix_kl$mix

    ftotal <- list(f_plug = f_plug, f_cv = f_cv, f_imde = f_imde,
                   f_kdir = f_kdir, f_dmde = f_dmde,
                   f_cvise = f_cvise, f_cvkl = f_cvkl)
    ielst <- lapply(ftotal, function(f) IE(mix0_ie, f))
    miae[i, ] <- c(no, sapply(ielst, `[[`, "iae")) # IAE,
    miseR[i, ] <- c(no, sapply(ielst, `[[`, "iseR")) # ISE on real space
    kld[i, ] <- c(no, sapply(ftotal, function(f) KLD(mix0_kld, f)))
    tm[i, ] <- c(no, diff(c(t1, t2, t3, t4, t5, t6, t7)))
  }
  list(miae = miae, miseR = miseR, kld = kld, time = tm)
}

