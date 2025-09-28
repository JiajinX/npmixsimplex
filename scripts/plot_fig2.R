invisible(lapply(dir("R", full.names = TRUE), source))
pth <- "results/figures/3.3-"
axl <- 500

set.seed(1)
mix <- Dir_mix(4)
# contour_r2s(mix)
data <- genemix(500, mix)
# plot(mix, data = data, prop = TRUE)
h <- 0.05
by <- 1

# pdf(paste0(pth, "gradient.pdf"), width = 6, height = 5)
mix1 <- initial_dirichletmix(data, h, by = 1)
mix1 <- EM(data, mix1, tem = 5)
plot(mix1, data = data, prop = TRUE, axl = axl)
# no need combine support points in EM step
# Expand new support points
w <- ddirichletmix(data, mix1, log = TRUE) # log f(x_i;G)
t <- initial_MEM(10, mix1, data, exp(-w))
t <- MEM(data, t, w, plot = TRUE, tem = 100, axl = axl)
np <- select_points(t, mix1$alpha, by = by)
mix2 <- CNM(np, mix1, data, w, by)
mix2 <- collapse_points(mix2, data)
plot(mix2, data = data, prop = TRUE, mode = FALSE, axl = axl)
common <- c(1, 3, 4, 5, 7)
points_simplex(attr(mix2, "mode")[common, ], add = TRUE, pch = 16, cex = 0.8)
points_simplex(attr(mix2, "mode")[-common, ],
               add = TRUE, pch = 16, cex = 0.8, col = "blue")
# dev.off()
