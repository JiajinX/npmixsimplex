# Create multivariate normal distribution from Chaco'S parameters
chaconmix <- function(no) {
  mix <- switch(no,
    list(mus = rbind(c(0, 0)),
         Sigmas = ks::invvech(c(1 / 4, 0, 1)),
         pr = 1),
    list(mus = rbind(c(0, 0)),
         Sigmas = ks::invvech(c(1, 9 / 10, 1)),
         pr = 1),
    list(mus = matrix(rep(c(3, -3), 8) *
                        rep((1 - (4 / 5)^(0:(8 - 1))), each = 2),
                      byrow = TRUE, nrow = 8),
         Sigmas = matrix(rep(ks::invvech(c(1, -9 / 10, 1)), 8) *
                           rep((4 / 5)^(2 * (0:(8 - 1))), each = 4),
                         byrow = TRUE, nrow = 8 * 2),
         pr = rep(1, 8) / 8),
    list(mus = rbind(c(-3 / 2, 3 / 2), c(3 / 2, -3 / 2), c(0, 0)),
         Sigmas = rbind(diag(9 / 16, 2), diag(9 / 16, 2),
                        ks::invvech(c(4 / 5, -18 / 25, 4 / 5)) * 9 / 16),
         pr = c(4, 4, 3) / 11),
    list(mus = rbind(c(0, 0), c(0, 0)),
         Sigmas = rbind(ks::invvech(c(1, 1, 4)),
                        ks::invvech(c(4 / 9, -1 / 3, 4 / 9))),
         pr = c(2, 1) / 3),
    list(mus = rbind(c(1, 0), c(-1, 0)),
         Sigmas = rbind(diag(4 / 9, 2), diag(4 / 9, 2)),
         pr = c(1, 1) / 2),
    list(mus = rbind(c(-1, 1), c(1, -1)),
         Sigmas = rbind(ks::invvech(c(4 / 9, 4 / 15, 4 / 9)),
                        ks::invvech(c(4 / 9, 4 / 15, 4 / 9))),
         pr = c(1, 1) / 2),
    list(mus = rbind(c(1, -1), c(-1, 1)),
         Sigmas = rbind(ks::invvech(c(4 / 9, 14 / 45, 4 / 9)), diag(4 / 9, 2)),
         pr = c(1, 1) / 2),
    list(mus = rbind(c(-1, 0), c(1, 2 / sqrt(3)), c(1, -2 / sqrt(3))),
         Sigmas = rbind(ks::invvech(c(9, 63 / 10, 49 / 4)),
                        ks::invvech(c(9, 0, 49 / 4)),
                        ks::invvech(c(9, 0, 49 / 4))) / 25,
         pr = c(3, 3, 1) / 7),
    list(mus = matrix(c(rep(0, 4), c(-1, -1, 1, -1, -1, 1, 1, 1)),
                      byrow = TRUE, nrow = 6),
         Sigmas = rbind(diag(2),
                        matrix(rep(diag(1 / 16, 2), 5),
                               byrow = TRUE, nrow = 5 * 2)),
         pr = c(1 / 2, rep(1 / 10, 5))),
    list(mus = rbind(c(-3 / 2, 0), c(3 / 2, 0), cbind(-1:1 - 3 / 2, -1:1),
                     c(0, 0), cbind(-1:1 + 3 / 2, -1:1)),
         Sigmas = matrix(c(rep(ks::invvech(c(4 / 9, 4 / 15, 4 / 9)), 2),
                           rep(ks::invvech(c(1 / 15, 1 / 25, 1 / 15)), 3) / 15,
                           ks::invvech(c(1, 3 / 5, 1)) / 9,
                           rep(ks::invvech(c(1 / 15, 1 / 25, 1 / 15)), 3) / 15),
                         byrow = TRUE, nrow = 9 * 2),
         pr = c(12 / 25, 12 / 25, rep(1 / 350, 3), 8 / 350, rep(1 / 350, 3))),
    list(mus = matrix(c(rep(0, 4), c(1, 1, -1, 1, -1, -1, 1, -1)),
                      byrow = TRUE, nrow = 6),
         Sigmas = rbind(diag(2),
                        ks::invvech(c(1, -9 / 10, 1)) / 16,
                        ks::invvech(c(1, -9 / 10, 1)) / 4,
                        diag(1 / 8, 2),
                        ks::invvech(c(1, -9 / 10, 1)) / 8,
                        diag(1 / 16, 2)),
         pr = c(1 / 2, 3 / 40, 1 / 5, rep(3 / 40, 3)))
  )
  class(mix) <- "mvnormix"
  mix
}

rmvnormix <- function(n, mix) {
  ks::rmvnorm.mixt(n, mix$mus, mix$Sigmas, props = mix$pr)
}

dmvnormix <- function(x, mix) {
  ks::dmvnorm.mixt(x, mix$mus, mix$Sigmas, props = mix$pr)
}

# Log-transform multivariate normal density on simplex
contour.mvnormix <- function(mix, data = NULL,
                             lv = NULL, plev = NULL, axl = 100) {
  x1 <- seq(0, 2 / sqrt(3), l = axl)
  x2 <- seq(0, 1, l = axl)
  x <- as.matrix(expand.grid(x1, x2))
  x <- cbind(x[, 2], 1 - (sqrt(3) * x[, 1] + x[, 2]) / 2)
  x <- cbind(x, 1 - rowSums(x))
  i <- rowSums(x > 0) == ncol(x)

  d <- numeric(length = nrow(x))
  sqrtD <- sqrt(ncol(x))
  ilrx <- unclass(compositions::ilr(x[i, ]))
  d[i] <- dmvnormix(ilrx, mix) / exp(rowSums(log(x[i, ]))) / sqrtD
  dm <- matrix(d, nrow = axl, ncol = axl)
  if (is.null(lv)) {
    if (is.null(plev)) {
      plev <- c(50, 65, 75, 80, 85, 90, 94, 97, 99, 100) / 100
    }
    lv <- quantile(d[i], probs = plev) # Proportion levels
  }
  simplexpanel()
  if (!is.null(data)) points.simplex(data, col = "gray70", add = TRUE)
  contour(x1, x2, dm, levels = lv, xlab = "x", ylab = "y",
          col = viridis::viridis(length(lv)), add = TRUE)
  lv
}

# Plot gaussian mixture distribution density on real space
contourR <- function(mix, data = NULL, x1 = c(-3, 3), x2 = c(-3, 3),
                     lv = NULL, plev = NULL, axl = 100) {
  x1 <- seq(x1[1], x1[2], l = axl)
  x2 <- seq(x2[1], x2[2], l = axl)
  x <- as.matrix(expand.grid(x1, x2))
  d <- dmvnormix(x, mix)
  dim(d) <- c(axl, axl)
  if (is.null(lv)) {
    if (is.null(plev)) {
      plev <- c(50, 65, 75, 80, 85, 90, 94, 97, 99, 100) / 100
    }
    lv <- quantile(d, probs = plev) # Proportion levels
  }
  contour(x1, x2, d, levels = lv, xlab = "x", ylab = "y",
          col = viridis::viridis(length(lv)))
  if (!is.null(data)) points(data, col = "gray70", cex = 0.5)
  points(mix$mus, col = 2, pch = 16, cex = 0.6)
  lv
}

# Plot Dirichlet mixture density contour on real space
contourR_dirichletmix <- function(mix, data = NULL, ilr = TRUE,
                                  lv = NULL, plev = NULL, axl = 100) {
  mode <- attr(mix, "mode")
  tdata <- if (ilr) compositions::ilr(data) else compositions::alr(data)
  tmode <- if (ilr) compositions::ilr(mode) else compositions::alr(mode)
  rg <- apply(tdata, 2, range)
  x1 <- seq(rg[1, 1], rg[2, 1], l = 100)
  x2 <- seq(rg[1, 2], rg[2, 2], l = 100)
  x <- as.matrix(expand.grid(x1, x2))
  tx <- if (ilr) compositions::ilrInv(x) else compositions::alrInv(x)
  d <- ddirichletmix(tx, mix) * exp(rowSums(log(tx))) # ???
  if (ilr) d <- d * sqrt(ncol(tx))
  dm <- matrix(d, nrow = axl, ncol = axl)
  if (is.null(lv)) {
    if (is.null(plev)) {
      plev <- c(50, 65, 75, 80, 85, 90, 94, 97, 99, 100) / 100
    }
    lv <- quantile(d, probs = plev) # Proportion levels
  }
  contour(x1, x2, dm, levels = lv, xlab = "x", ylab = "y",
          col = viridis::viridis(length(lv)))
  if (!is.null(data)) points(tdata, col = "gray70", cex = 0.5)
  points(tmode, col = 2, pch = 16, cex = 0.6)
  lv
}
