invisible(lapply(dir("R", full.names = TRUE), source))
pth <- "results/figures/4-simulation/"
size <- 1000
axl <- 300
lv <-  seq(0.1, 0.85, by = 0.25)
cols <- paletteer::paletteer_c("ggthemes::Blue-Teal", 3)

# ilr-Gaussian mixtures for d = 2
# pdf(paste0(pth, "chacon.pdf"), width = 6, height = 5)
for (i in 1:12) {
  print(i)
  mix <- chaconmix(i)
  data <- genemix(size, mix)
  contour_r2s(mix, data, prop = TRUE, axl = axl, lv = lv)
}
# dev.off()

# ilr-Gaussian mixtures for d = 3
# rgl::par3d(windowRect = c(0, 0, 1500, 1500))
for (i in 1:12) {
  print(i)
  mix <- chaconmix(i, D = 3)
  data <- genemix(size, mix)
  contour3D(mix, data = data, prop = TRUE,
            n = 100, alo = 0.2, ahi = 0.2, col = cols)
  # rgl::rgl.snapshot(paste0(pth, "chacon", i, ".png"), fmt = "png")
}

# Dirichlet mixtures for d = 2
# pdf(paste0(pth, "dirmix.pdf"), width = 6, height = 5)
for (i in 1:12) {
  print(i)
  mix <- Dir_mix(i)
  data <- genemix(size, mix)
  plot(mix, data, prop = TRUE, axl = axl, lv = lv)
}
# dev.off()

# Dirichlet mixtures for d = 3
# rgl::par3d(windowRect = c(0, 0, 1500, 1500))
for (i in 1:12) {
  print(i)
  mix <- Dir_mix(i, D = 4)
  data <- genemix(1000, mix)
  contour3D(mix, data = data, prop = TRUE,
            n = 100, alo = 0.2, ahi = 0.2, col = cols)
  # rgl::rgl.snapshot(paste0(pth, "dirmix", i, ".png"), fmt = "png")
}

# Simulation time
estime <- readRDS("results/raw/simulation/estime")
logtime <- lapply(estime, log10)
cols <- c("#3366CC", "#0099C6", "#FF9900",
          "#109618", "#990099", "#DC3912")
cols <- transcols(cols, 0.5)
pdf(paste0(pth, "time.pdf"), width = 7, height = 6)
plot.new()
plot.window(xlim = c(1.7, 6.3), ylim = c(-1, 4.3))
axis(1, at = 2:6)
axis(2, at = seq(-1, 4, by = 1),
     labels = c("0.1", "1", "10", expression(10^2, 10^3, 10^4)), las = 2)
abline(h = -1:5, col = "grey90", lwd = 0.8)
abline(v = 2:7,  col = "grey90", lwd = 0.8)
box()
title(xlab = expression(italic(d)), ylab = "Time (seconds)", cex.lab = 1.4)
for (i in seq_along(logtime)) {
  m <- colMeans(logtime[[i]])
  points(2:6, m, col = cols[i], pch = 16, cex = 1.2, lwd = 1.6)
  lines(2:6, m, col = cols[i], lwd = 1.6)
}
legname <- names(logtime)
legend("topleft", bg = NA, bty = "n", inset = c(0.02, 0.02),
       legend = legname, col = cols, pch = 16,
       lwd = 1.6, pt.cex = 1.2, lty = 1)
dev.off()





