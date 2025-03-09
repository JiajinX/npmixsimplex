
library(manipulate)
library(sphunif)
library(DirStats)
library(rgl)
library(viridis)
library(rotasym)
library(testthat)
library(extraDistr)

## Density over positive orthant of S^{p-1}

d_orth <- function(x, mu, kappa = 1, g = function(x) abs(x)) {
  exp(
  dens <- colSums(g(t(x) * (kappa * drop(mu))))
  )
  d <- length(mu) - 1
  cons <- d / (2 * rotasym::w_p(p = d) * sum(abs(kappa * mu)))
  cons * dens
}
rotasym::w_p(p = 3) * mean(d_orth(r_unif_sph(n = 1e5, p = 3, M = 1)[, , 1],
                                  mu = c(sqrt(2), sqrt(2), 0) / 2,
                                  kappa = c(1, 2, 3)))

gr <- expand.grid(th = seq(0, 2 * pi, l = 201), ph = seq(0, pi, l = 101))
x <- to_sph(gr$th, gr$ph)

mu <- r_unif_sphere(n = 1, p = 3)
dens <- d_orth(x, mu = mu, kappa = 1)
cuts <- cut(dens, breaks = seq(0, 0.15, l = 10))
col <- viridis(length(levels(cuts)))[cuts]
plot3d(x, col = col, size = 10, xlab = "x", ylab = "y", zlab = "z")
       #xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), zlim = c(-1.1, 1.1))
e <- expand.grid(c(-1, 1), c(-1, 1), c(-1, 1))
for (j in 1:nrow(e)) {
  points3d(1.01 * mu * e[j, ], col = "red", size = 15)
}
planes3d(a = 1, b = 0, c = 0, d = 0, col = "lightblue", alpha = 0.75)
planes3d(a = 0, b = 1, c = 0, d = 0, col = "lightblue", alpha = 0.75)
planes3d(a = 0, b = 0, c = 1, d = 0, col = "lightblue", alpha = 0.75)

## Orthogonal basis over the positive orthant of S^{p-1}

p <- 3
k <- 5
ms <- as.matrix(do.call(expand.grid, c(rep(list(0:k), p - 1), list(0:1))))
even_ms <- apply(ms %% 2 == 0, 1, all)
ms_even <- ms[even_ms, , drop = FALSE]

f <- 0
for (i in seq_len(nrow(ms_even))) {
  # f <- f + ((-1)^(i + 1) / i^3) * g_i_k(x = x, m = ms_even[i, ])
  f <- f + runif(1) * g_i_k(x = x, m = ms_even[i, ])
}
col <- viridis(nrow(x))[rank(f)]
plot3d(x, col = col, size = 10, xlab = "x", ylab = "y", zlab = "z")
#xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), zlim = c(-1.1, 1.1))
planes3d(a = 1, b = 0, c = 0, d = 0, col = "lightblue", alpha = 0.75)
planes3d(a = 0, b = 1, c = 0, d = 0, col = "lightblue", alpha = 0.75)
planes3d(a = 0, b = 0, c = 1, d = 0, col = "lightblue", alpha = 0.75)

min(f)
all(f > 0)

test_that("m such that all m_i are even gives orthant symmetry (if and only if)", {

  for (p in 2:6) {

    x <- r_unif_sph(n = 1e2, p = p, M = 1)[, , 1]
    x_sh <- x * matrix(sample(c(-1, 1), size = length(x), replace = TRUE),
                       ncol = ncol(x), nrow = nrow(x))

    for (k in 1:10) {

      ms <- as.matrix(do.call(expand.grid, c(rep(list(0:k), p - 1), list(0:1))))
      ms <- ms[rowSums(ms) == k, ]
      even_ms <- apply(ms %% 2 == 0, 1, all)
      ms_even <- ms[even_ms, , drop = FALSE]
      ms_odd <- ms[!even_ms, , drop = FALSE]

      for (i in seq_len(nrow(ms_even))) {

        m <- ms_even[i, ]
        expect_true(max(abs(g_i_k(x = x, m = m) -
                              g_i_k(x = x_sh, m = m))) < 1e-10)

      }
      for (i in seq_len(nrow(ms_odd))) {

        m <- ms_odd[i, ]
        expect_false(max(abs(g_i_k(x = x, m = m) -
                               g_i_k(x = x_sh, m = m))) < 1e-10)

      }

    }

  }

})

## Density estimation

# Simulate some data on the simplex
n <- 1e4
p <- 3
alpha <- rep(2, p)
samp <- rdirichlet(n = n, alpha = alpha)
plot(samp, pch = 16, cex = 0.5)

# Estimator
f_hat_orth <- function(x, X, k, simplex = FALSE, show_m = FALSE) {

  # Dimension
  p <- ncol(X)
  stopifnot(p == ncol(x))

  # Check positive orthant sphere or simplex constraints
  stopifnot(all(X >= 0), all(X <= 1))
  if (simplex) {

    stopifnot(all(abs(rowSums(X) - 1) < 1e-10))
    stopifnot(all(abs(rowSums(x) - 1) < 1e-10))
    X_sph <- sqrt(X)
    x_sph <- sqrt(x)

  } else {

    stopifnot(all(abs(rowSums(X^2) - 1) < 1e-10))
    stopifnot(all(abs(rowSums(x^2) - 1) < 1e-10))
    X_sph <- X
    x_sph <- x

  }

  # Generate ms_even from order of the polynomial k (take only even k's)
  k <- ifelse(k %% 2, k - 1, k)
  ms <- as.matrix(do.call(expand.grid, c(rep(list(0:k), p - 1), list(0:1))))
  ms <- ms[rowSums(ms) <= k, , drop = FALSE]
  even_ms <- apply(ms %% 2 == 0, 1, all)
  ms_even <- ms[even_ms, , drop = FALSE]
  if (show_m) {
    message("Number of even m's: ", nrow(ms_even))
    print(ms_even)
  }

  # Estimate coefficients
  coefs_even <- numeric(nrow(ms_even))
  for (i in seq_len(nrow(ms_even))) {

    coefs_even[i] <- mean(g_i_k(x = X_sph, m = ms_even[i, ]))

  }
  coefs_even <- coefs_even / rotasym::w_p(p = p)

  # Construct density estimate
  f <- 0
  for (i in seq_len(nrow(ms_even))) {

    f <- f + coefs_even[i] * g_i_k(x = x_sph, m = ms_even[i, ])

  }

  # Return density on the simplex (instead of on the positive orthant)?
  if (simplex) {

    f <- 2 * f / apply(x_sph, 1, prod)
    # The 1 / 2^(p-1) cancels partially with the symmetrization factor 2^p
    # Recall the 2^(p-1) because one variable is redundant

  } else {

    f <- f * 2^p # Symmetrization factor

  }
  return(f)

}


test_that("f_hat_orth approximates well a Dirichlet density", {

  p <- 3
  n <- 100
  samp_simp <- rdirichlet(n = n, alpha = 1 + 1:p)
  samp_unif_simp <- rdirichlet(n = 1e4, alpha = rep(1, p))
  dens <- f_hat_orth(x = samp_unif_simp, X = samp_simp, k = 10,
                         simplex = TRUE)
  dens0 <- ddirichlet(x = samp_unif_simp, alpha = 1 + 1:p)
  plot(samp_unif_simp[, -3], col = viridis(10)[cut(dens, breaks = 10)])
  plot(samp_unif_simp[, -3], col = viridis(10)[cut(dens0, breaks = 10)])

})

test_that("f_hat_orth integrates one on the positive orthant of S^{p-1}", {

  for (p in 3:6) {

    skip_if(p == 2, message = "p = 2 TODO")
    n <- 100
    samp_sph <- sqrt(rdirichlet(n = n, alpha = rep(2, p)))
    int <- int_sph_MC(f = f_hat_orth, p = p, M = 1e4, X = samp_sph, k = 4,
                      simplex = FALSE) / 2^p
    expect_equal(int, 1, tolerance = 0.1)

  }

})

test_that("f_hat_orth integrates one on the simplex ∆^{p-1}", {

  for (p in 3:6) {

    skip_if(p == 2, message = "p = 2 TODO")
    n <- 100
    samp_simp <- rdirichlet(n = n, alpha = rep(2, p))
    samp_unif_simp <- rdirichlet(n = 1e4, alpha = rep(1, p))
    area <- 1 / ddirichlet(x = rep(1 / p, p), alpha = rep(1, p))
    int <- mean(f_hat_orth(x = samp_unif_simp, X = samp_simp, k = 4,
                           simplex = TRUE)) * area
    expect_equal(int, 1, tolerance = 0.1)

  }

})

p <- 3
n <- 1e4
k <- 4
alpha <- (1:p)^2
samp_sph <- rbind(sqrt(rdirichlet(n = n, alpha = alpha)),
                  sqrt(rdirichlet(n = n, alpha = rev(alpha))))
samp_simp <- samp_sph^2
f_hat_sph <- f_hat_orth(x = samp_sph, X = samp_sph, k = k, show_m = TRUE)
f_hat_simp <- f_hat_orth(x = samp_simp, X = samp_simp, k = k, show_m = TRUE,
                         simplex = TRUE)

nbreaks <- 10
probs <- seq(0, 1, l = nbreaks + 1)
f_dir_sph <- function(x, alpha) {
  2^(ncol(x) - 1) * ddirichlet(x = x^2, alpha = alpha) * apply(x, 1, prod)
}
f_dir_simp <- function(x, alpha) ddirichlet(x = x, alpha = alpha)

f_sph <- 0.5 * f_dir_sph(x = samp_sph, alpha = alpha) +
  0.5 * f_dir_sph(x = samp_sph, alpha = rev(alpha))
f_simp <- 0.5 * f_dir_simp(x = samp_simp, alpha = alpha) +
  0.5 * f_dir_simp(x = samp_simp, alpha = rev(alpha))

brks <- c(-Inf, seq(0, 1.1 * max(f_sph), l = nbreaks), Inf)
col_f_dir <- c("red", viridis(nbreaks), "black")[cut(f_sph, breaks = brks)]
col_f_hat <- c("red", viridis(nbreaks), "black")[cut(f_hat_sph, breaks = brks)]
plot(samp_sph[, 1:2], col = col_f_dir, pch = 16, cex = 0.5, xlim = c(0, 1), ylim = c(0, 1))
plot(samp_sph[, 1:2], col = col_f_hat, pch = 16, cex = 0.5, xlim = c(0, 1), ylim = c(0, 1))
table(f_hat_sph >= 0)
summary(f_sph)
summary(f_hat_sph)

brks <- c(-Inf, seq(0, 1.1 * max(f_simp), l = nbreaks), Inf)
col_f_dir <- c("red", viridis(nbreaks), "black")[cut(f_simp, breaks = brks)]
col_f_hat <- c("red", viridis(nbreaks), "black")[cut(f_hat_simp, breaks = brks)]
pairs(samp_simp[, 1:2], col = col_f_dir, pch = 16, cex = 0.5, xlim = c(0, 1), ylim = c(0, 1))
pairs(samp_simp[, 1:2], col = col_f_hat, pch = 16, cex = 0.5, xlim = c(0, 1), ylim = c(0, 1))
table(f_hat_simp >= 0)
summary(f_hat_simp)

## Visualize Dirichlet distribution

x <- seq(0, 1, length.out = 200)
y <- seq(0, 1, length.out = 200)
xyz <- as.matrix(expand.grid(x, y))
xyz <- cbind(xyz, 1 - rowSums(xyz))
nlev <- 15

manipulate({

  alpha <- c(alpha1, alpha2, alpha3)
  z <- matrix(ddirichlet(x = xyz, alpha = alpha), nrow = length(x),
              ncol = length(y))
  contour(x, y, z, levels = seq(0, max(z, na.rm = TRUE), l = nlev + 1),
          xlab = "x", ylab = "y", col = viridis(nlev),
          main = "Dirichlet density")
  abline(a = 1, b = -1)
  abline(h = 0)
  abline(v = 0)

},
alpha1 = slider(0, 10, step = 0.01, initial = 1),
alpha2 = slider(0, 10, step = 0.01, initial = 2),
alpha3 = slider(0, 10, step = 0.01, initial = 3)
)

manipulate({

  stopifnot(mu1 + mu2 <= 1)
  j <- rep(1, 3)
  mu <- c(mu1, mu2, 1 - mu1 - mu2)
  alpha <- j + mu / h
  z <- matrix(ddirichlet(x = xyz, alpha = alpha), nrow = length(x),
              ncol = length(y))
  contour(x, y, z, levels = seq(0, max(z, na.rm = TRUE), l = nlev + 1),
          xlab = "x", ylab = "y", col = viridis(nlev),
          main = "Dirichlet density")
  points(mu1, mu2, col = 2, pch = 16)
  abline(a = 1, b = -1)
  abline(h = 0)
  abline(v = 0)

},
h = slider(0, 1, step = 0.01, initial = 0.15),
mu1 = slider(0, 1, step = 0.01, initial = 0.72),
mu2 = slider(0, 1, step = 0.01, initial = 0.25)
)


manipulate({

  stopifnot(mu1 + mu2 <= 1)
  j <- rep(1, 3)
  mu <- c(mu1, mu2, 1 - mu1 - mu2)
  alpha <- j / (3 * h)

  x_minus_mu <- scale(xyz, center = FALSE, scale = mu)
  x_minus_mu <- x_minus_mu / rowSums(x_minus_mu)

  z <- matrix(ddirichlet(x = x_minus_mu, alpha = alpha), nrow = length(x),
              ncol = length(y))
  contour(x, y, z, levels = seq(0, max(z, na.rm = TRUE), l = nlev + 1),
          xlab = "x", ylab = "y", col = viridis(nlev),
          main = "Dirichlet density")
  points(mu1, mu2, col = 2, pch = 16)
  abline(a = 1, b = -1)
  abline(h = 0)
  abline(v = 0)

},
h = slider(0, 1, step = 0.01, initial = 0.15),
mu1 = slider(0, 1, step = 0.01, initial = 0.72),
mu2 = slider(0, 1, step = 0.01, initial = 0.25)
)
