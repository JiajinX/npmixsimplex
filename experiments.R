invisible(lapply(dir("R", full.names = TRUE), source))

cs_mix <- function(no) {
  p <- switch(no,
    list(mu = rbind(c(0.2, 0.5, 0.3), c(0.7, 0.1, 0.2)),
         pr = c(0.5, 0.5), h = 0.1),
    list(mu = rbind(c(0.5, 0.05, 0.45), c(0.95, 0.01, 0.04)),
         pr = c(0.3, 0.7), h = 0.1),
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

# Set some of data value less than tol equal to 0
setbound_data <- function(data, tol = 0.01) {
  close0 <- data < tol
  n0 <- sample(sum(close0), 1) # randomly choose n0 data to be 0
  idx0 <- sample(which(close0), n0)
  data[idx0] <- 0
  data / rowSums(data)
}

# -----------------------------
# set.seed(1019)
no <- 2
n <- 500
h <- 0.1
mix0 <- cs_mix(no)
data <- rdirichletmix(n, mix0)
data <- setbound_data(data)
r <- dmde(data)
# Rprof()
# r <- cnmdm(data, h = h, verbose = 1)
# Rprof(NULL)
# summaryRprof(memory = "none")$by.total
mix <- r$mix
mix

lv <- pairs(mix0, data)
pairs(mix, data, lv = lv)

# -----------------------------
# Test Chacon mixture
no <- 2
n <- 1000
mix0 <- chaconmix(no)
data <- rmvnormix(n, mix0)
data <- compositions::ilrInv(data)
data <- unclass(data)
r <- dmde(data)
mix <- r$mix

lv <- contour(mix0, data)
plot(mix, data = data, lv = lv)

# -----------------------------
# Test of time

feature <-
  c("dimentions", "size", "loop", "h", "iter_number",
    "initial_points", "final_points", "ll", "gradient", "total", "iteration",
    "initial", "EM", "initial_MEM", "MEM", "select", "CNM", "drop")

# For d = 2, 4, 8, 16, 32, 64
# Fixed n = 100 and h = 0.1
size <- 100
D <- 2^(1:6) + 1
H <- 0.1
loop <- 10
dinfo <- data.frame(matrix(nrow = length(D) * loop, ncol = 18))
colnames(dinfo) <- feature
k <- 1
for (i in seq_along(D)) {
  for (j in seq_len(loop)) {
    print(k)
    d <- D[i]
    # h <- size^(-1 / (4 + d)) * 0.2 # C = 0.02
    mix <- dirichletmix(rep(1 / d, d) / H + 1)
    data <- rdirichletmix(size, mix)
    r <- cnmdm(data, h = H / 2)
    dinfo[k, ] <- c(d, size, j, H / 2, r)
    k <- k + 1
  }
}

# For n = 64, 128, 256, 512, 1042, 2048, 4096
# Fixed D = 2 and h = 0.1
size <- 2^(6:12)
D <- 2 + 1
H <- 0.1
loop <- 10
ninfo <- data.frame(matrix(nrow = length(size) * loop, ncol = 18))
colnames(ninfo) <- feature
k <- 1
for (i in seq_along(size)) {
  for (j in seq_len(loop)) {
    print(k)
    n <- size[i]
    h <- n^(-1 / (4 + D)) * 0.3
    mix <- dirichletmix(rep(1 / D, D) / H + 1)
    data <- rdirichletmix(n, mix)
    r <- cnmdm(data, h = h, maxit = 200, verbose = 0)
    ninfo[k, ] <- c(D, n, j, h, r)
    k <- k + 1
  }
}

# -----------------------------
#                                               D = 3,   10,   30,   50
# h = 0.1, mode in center, with D increasing,   r = 9.3, 18.2, 40.5, 61.5
# h = 0.1, mode in corner, with D increasing,   r = 8.6, 17.4, 39.5, 60.4

#                                               D = 3,   10,   30,   50
# h = 0.01, mode in center, with D increasing,  r = 9,   18.6, 42.9, 65.3
# h = 0.01, mode in corner, with D increasing,  r = 7.9, 16.9, 39.4, 60.5

#                                                D = 3,   10,   30,   50
# h = 0.001, mode in center, with D increasing,  r = 8.9, 18.7, 43.8, 67.3
# h = 0.001, mode in corner, with D increasing,  r = 7.9, 16.9, 39.4, 60.6

# D <- 3
# h <- 0.01
# a0 <- 1 / h + D
# mix0 <- dirichletmix(rep(a0 / D, D))
# cst <- 2 * sqrt(a0 + 1)
# data <- rdirichletmix(10000, mix0)
# tdata <- cst * asin(sqrt(data)) # Variance stabilization transformation
# tmu <- cst * asin(sqrt(colMeans(data)))
# ic = 1
# plot(density(tdata[, ic]), type = 'l')
# x <- seq(tmu[ic] - 4, tmu[ic] + 4, 0.1)
# lines(x, dnorm(x, tmu[ic]), col = 2)

# d = 3:100
# r2 = qchisq(0.95, df = 3:100); r = sqrt(r2 / d)
# dnorm(sqrt(r2 / d))^d / dnorm(0)^d
