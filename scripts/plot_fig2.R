invisible(lapply(dir("R", full.names = TRUE), source))
pth <- "results/figures/3.3-"

mix <- Dir_mix(4)
data <- genemix(500, mix)
h <- 0.05

mix <- initial_dirichletmix(data, h, by = 1)
mix <- EM(data, mix, tem = 5)
# no need combine support points in EM step
# Expand new support points
w <- ddirichletmix(data, mix, log = TRUE) # log f(x_i;G)
t <- initial_MEM(20, mix, data, exp(-w))
# pdf(paste0(pth, "gradient.pdf"), width = 6, height = 5)
t <- MEM(data, t, w, plot = TRUE, tem = 100)
# dev.off()
