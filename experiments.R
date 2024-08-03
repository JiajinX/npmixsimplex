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
                    c(1, 0 ,0 ,0), c(0.1, 0.1, 0.25, 0.55)),
         pr = c(0.3, 0.5, 0.1, 0.1), h = 0.1),
    list(mu = rbind(c(0, 0.1, 0.5, 0.1, 0.3), c(0.5, 0.05, 0.05, 0, 0.4)),
         pr = c(0.8, 0.2), h = 0.1)
  )
  mix0 <- dirichletmix(a = p$mu / p$h + 1, pr = p$pr)
  mix0
}

# -----------------------------
# set.seed(1019)
no <- 6
n <- 500
h <- 0.1
mix0 <- cs_mix(no)
data <- rdirichletmix(n, mix0)
# Rprof()
r <- cnmdm(data, h = h, verbose = 1)
# Rprof(NULL)
# summaryRprof(memory = "none")$by.total
mix <- r$mix
mix

lv <- pairs(mix0, data)
pairs(mix, data)
# -----------------------------
# Test Chacon mixture
no <- 2
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
