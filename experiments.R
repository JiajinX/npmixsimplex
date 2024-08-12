invisible(lapply(dir("R", full.names = TRUE), source))

cs_mix <- function(no) {
  p <- switch(no,
    list(mu = rbind(c(0.2, 0.5, 0.3), c(0.7, 0.1, 0.2)),
         pr = c(0.5, 0.5), h = 0.1),
    list(mu = rbind(c(0.2, 0.5, 0.3), c(0.7, 0.1, 0.2),
                    c(0.4, 0.4, 0.2), c(0.1, 0.7, 0.2)),
         pr = c(0.1, 0.4, 0.2, 0.3), h = 0.1),
    list(mu = rbind(c(0, 0.2, 0.8), c(1, 0, 0),
                    c(1 / 3, 1 / 3, 1 / 3), c(0.6, 0.3, 0.1)),
         pr = c(0.1, 0.2, 0.3, 0.4), h = 0.1),
    list(mu = rbind(c(0.2, 0.3, 0.35, 0.15), c(0.45, 0.05, 0.2, 0.3)),
         pr = c(0.4, 0.6), h = 0.1),
    list(mu = rbind(c(0, 0, 0.7, 0.3), c(0.5, 0.45, 0.05, 0),
                    c(1, 0, 0, 0), c(0.1, 0.1, 0.25, 0.55)),
         pr = c(0.3, 0.5, 0.1, 0.1), h = 0.1),
    list(mu = rbind(c(0, 0.1, 0.5, 0.1, 0.3), c(0.5, 0.05, 0.05, 0, 0.4)),
         pr = c(0.8, 0.2), h = 0.1)
  )
  mix0 <- dirichletmix(a = p$mu / p$h + 1, pr = p$pr)
  mix0
}

# -----------------------------
# set.seed(1019)
no <- 1
n <- 500
h <- 0.1
mix0 <- cs_mix(no)
data <- rdirichletmix(n, mix0)
Rprof()
r <- cnmdm(data, h = h, verbose = 1)
Rprof(NULL)
summaryRprof(memory = "none")$by.total
mix <- r$mix
mix

lv <- pairs(mix0, data)
pairs(mix, data)
# -----------------------------
# Test Chacon mixture
no <- 12
n <- 500
h <- 0.1
mix0 <- chaconmix(no)
data <- rmvnormix(n, mix0)
data <- compositions::ilrInv(data)
data <- unclass(data)
r <- cnmdm(data, h = 0.01, verbose = 1)
mix <- r$mix

lv <- contour(mix0, data)
plot(mix, data = data, type = "both", lv = lv)
# -----------------------------
# Test contour plot
no <- 2
n <- 5e6
mix0 <- cs_mix(no)
data <- rdirichletmix(n, mix0)
lv <- plot(mix0)
contour_sample(data, split = 60)

mix0 <- chaconmix(no)
data <- rmvnormix(n, mix0)
data <- compositions::ilrInv(data)
data <- unclass(data)
contour(mix0)
contour_sample(data, split = 80)
# -----------------------------
# Test of time

# For d = 2, 4, 8, 16, 32, 64, 128
# Fixed n = 100 and h = 0.1
size <- 100
D <- 2^(1:7) + 1
H <- 0.1
loop <- 1
dimtime <- matrix(nrow = loop, ncol = length(D))
colnames(dimtime) <- D - 1
for (i in seq_along(D)) {
  for (j in seq_len(loop)) {
    d <- D[i]
    h <- size^(-1 / (4 + d))
    mix <- dirichletmix(rep(1 / d, d) / H + 1)
    data <- rdirichletmix(size, mix)
    start_time <- Sys.time()
    cnmdm(data, h = H, maxit = 1)
    end_time <- Sys.time()
    dimtime[j, i] <- end_time - start_time
  }
}
boxplot(dimtime)
lines(colMeans(dimtime), col = 2, type = "o", pch = 16)
boxplot(log2(dimtime))
lines(colMeans(log2(dimtime)), col = 2, type = "o", pch = 16)

# For n = 64, 128, 256, 512, 1042, 2048, 4096
# Fixed D = 2 and h = 0.1
size <- 2^(6:12)
D <- 2 + 1
H <- 0.1
loop <- 10
sizetime <- matrix(nrow = loop, ncol = length(size))
colnames(sizetime) <- size
for (i in seq_along(size)) {
  for (j in seq_len(loop)) {
    n <- size[i]
    h <- n^(-1 / (4 + D))
    mix <- dirichletmix(rep(1 / D, D) / H + 1)
    data <- rdirichletmix(n, mix)
    start_time <- Sys.time()
    cnmdm(data, h = H, maxit = 1)
    end_time <- Sys.time()
    sizetime[j, i] <- end_time - start_time
  }
}
boxplot(sizetime)
lines(colMeans(sizetime), col = 2, type = "o", pch = 16)
boxplot(log2(sizetime))
lines(colMeans(log2(sizetime)), col = 2, type = "o", pch = 16)
# -----------------------------
a0 <- 30
D <- 3
h <- 1 / (a0 - D)
a <- rbind(c(10, 10), c(1, 1))
mix <- dirichletmix(cbind(a, a0 - rowSums(a)))
data <- rdirichletmix(1000, mix)
mu <- colMeans(data)

ax <- seq(0, 1, l = 100)
x <- as.matrix(expand.grid(ax, ax))
x <- cbind(x, 1 - rowSums(x))
i <- rowSums(x >= 0) == 3
x <- x[i, ]
# points_simplex(x, col = 1)
points_simplex(data, col = 2)

y <- asin(sqrt(x)) * 2 * sqrt(a0 + 1)
tdata <- 2 * sqrt(a0 + 1) * asin(sqrt(data))
plot3d(y, col = 2)
# plot3d(tdata, col = 2, add = T)
bd <- (c(1, 0, 0) / h + 1) / a0
tbd <- asin(sqrt(bd)) * 2 * sqrt(a0 + 1)
plot3d(t(tbd), type = 's', radius = 3, add = T, alpha = 0.2)
# bound of mu to ensure mode is inside simplex

points_simplex(c(1, 0, 0), pch = 16)
points_simplex((c(1, 0, 0) * a0 - 1) / (a0 - D), add = T, col = "blue", pch = 16)
# mode when mu is (1, 0, 0)

tfsimplexpanel <- function(a0, simplex = FALSE) {
  cst <- 2 * sqrt(a0 + 1) # Constant value in transformation
  lim <- cst * asin(sqrt(1))
  rg = c(-0.01, lim)
  x <- y <- seq(0, lim, length.out = 100)
  z <- outer(x, y, function(x, y) {
    temp <- 1 - sin(x / cst)^2 - sin(y / cst)^2
    temp[temp < 0] <- 0
    cst * asin(sqrt(temp))
  })
  persp3d(x, y, z, col = "lightblue", alpha = 0.6,
          xlim = rg, ylim = rg, zlim = rg)
}

D <- 3:64
lp <- 10
r1 <- matrix(nrow = lp, ncol = length(D))
for (i in seq_along(D)) {
  d <- D[i]
  a0 <- d * 10
  h <- 1 / (a0 - d)
  mix0 <- dirichletmix(rep(10, d))
  for (j in seq_len(lp)) {
    data <- rdirichletmix(1000, mix0)
    cst <- 2 * sqrt(a0 + 1) # Constant value in transformation
    tdata <- cst * asin(sqrt(data)) # Variance stabilization transformation
    tmu <- cst * asin(sqrt(colMeans(data)))
    # tmu should be the first center of cluster groups
    d2 <- colSums((t(tdata) - tmu)^2) # Distance from tdata to tmu
    r1[j, i] <- quantile(d2, 0.95)
  }
}
m1 <- colMeans(r1) # +1.3 each time
sqrt(m1/D)
fit1 <- lm(m1 ~ D + I(D^2))
plot(D, m1, type = 'l')
lines(D, predict(fit1, newdata = data.frame(x = D)), col = "blue")

D <- 3:64
lp <- 10
r2 <- matrix(nrow = lp, ncol = length(D))
for (i in seq_along(D)) {
  d <- D[i]
  a0 <- d * 10
  h <- 1 / (a0 - d)
  mix0 <- dirichletmix(c(a0 - d + 1, rep(1, d - 1)))
  for (j in seq_len(lp)) {
    data <- rdirichletmix(1000, mix0)
    cst <- 2 * sqrt(a0 + 1) # Constant value in transformation
    tdata <- cst * asin(sqrt(data)) # Variance stabilization transformation
    tmu <- cst * asin(sqrt(colMeans(data)))
    # tmu should be the first center of cluster groups
    d2 <- colSums((t(tdata) - tmu)^2) # Distance from tdata to tmu
    r2[j, i] <- quantile(d2, 0.95)
  }
}
m2 <- colMeans(r2)
plot(D, m2, type = 'l')
fit2 <- lm(m2 ~ D + I(D^2))
lines(D, predict(fit2, newdata = data.frame(D)), col = "red")

D <- 3
h <- 0.01
a0 <- 1 / h + D
mix0 <- dirichletmix(rep(a0 / D, D))
# mix0 <- dirichletmix(c(a0 - D + 1, rep(1, D - 1)))
cst <- 2 * sqrt(a0 + 1) # Constant value in transformation
t <- 100
r <- rep(NA, t)
for (i in 1:t) {
  data <- rdirichletmix(1000, mix0)
  tdata <- cst * asin(sqrt(data)) # Variance stabilization transformation
  tmu <- cst * asin(sqrt(colMeans(data)))
  d2 <- colSums((t(tdata) - tmu)^2)
  r[i] <- quantile(d2, 0.95)
}
print(mean(r))

#                                               D = 3,   10,   30,   50
# h = 0.1, mode in center, with D increasing,   r = 9.3, 18.2, 40.5, 61.5
# h = 0.1, mode in corner, with D increasing,   r = 8.6, 17.4, 39.5, 60.4

#                                               D = 3,   10,   30,   50
# h = 0.01, mode in center, with D increasing,  r = 9,   18.6, 42.9, 65.3
# h = 0.01, mode in corner, with D increasing,  r = 7.9, 16.9, 39.4, 60.5

#                                                D = 3,   10,   30,   50
# h = 0.001, mode in center, with D increasing,  r = 8.9, 18.7, 43.8, 67.3
# h = 0.001, mode in corner, with D increasing,  r = 7.9, 16.9, 39.4, 60.6
