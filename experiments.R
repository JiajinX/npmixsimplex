invisible(lapply(dir("R", full.names = TRUE), source))

ISE <- function(tmix, emix, n = 10000) {
  i <- sample.int(n = 2, size = n, prob = c(0.5, 0.5), replace = TRUE)
  data1 <- rdirichletmix(sum(i == 1), tmix)
  data2 <- rdirichletmix(sum(i == 2), emix)
  data <- rbind(data1, data2)
  f1 <- ddirichletmix(data, tmix)
  f2 <- ddirichletmix(data, emix)
  ise <- 2 * mean((f1 - f2)^2 / (f1 + f2))
  ise
}

compr_mix <- function(mix0, mix) {
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(2, 2), mar = rep(0, 4))
  points(mix0)
  d0 <- contour(mix0)
  points(mix)
  contour(mix, lv = attr(d0, "level"))
  par(old_par)
  ISE(mix0, mix)
}

cs_mix <- function(no) {
  p <- switch(no,
    list(mu = rbind(c(0.2, 0.5, 0.3), c(0.7, 0.1, 0.2)),
         pr = c(0.5, 0.5), h = 0.1),
    list(mu = rbind(c(0.2, 0.5, 0.3), c(0.7, 0.1, 0.2),
                    c(0.4, 0.4, 0.2), c(0.1, 0.7, 0.2)),
         pr = c(0.1, 0.4, 0.2, 0.3), h = 0.1),
    list(mu = rbind(c(0.2, 0.5, 0.3), c(0.7, 0.1, 0.2),
                    c(0.4, 0.4, 0.2), c(0.9, 0.05, 0.05)),
         pr = c(0.1, 0.4, 0.2, 0.3), h = 0.1),
    list(mu = rbind(c(0, 0.2, 0.8), c(1, 0, 0),
                    c(1 / 3, 1 / 3, 1 / 3), c(0.6, 0.3, 0.1)),
         pr = c(0.1, 0.2, 0.3, 0.4), h = 0.1),
    list(mu = rbind(c(0.2, 0.5, 0.3), c(0.7, 0.1, 0.2), c(0.4, 0.4, 0.2)),
         pr = c(0.5, 0.2, 0.3), h = 0.2)
  )
  mix0 <- dirichletmix(a = p$mu / p$h + 1, pr = p$pr)
  mix0
}

# set.seed(1019)
no <- 1
n <- 500
mix0 <- cs_mix(no)
data <- rdirichletmix(n, mix0)
# points.simplex(data)
# contour(mix0, data)
r <- cnmdm(data, h = 0.1, nump = 2, verbose = 1)
mix <- r$mix
compr_mix(mix0, mix)
