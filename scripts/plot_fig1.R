invisible(lapply(dir("R", full.names = TRUE), source))
pth <- "results/figures/2.1-"

axl <- 300
lv <- seq(0.1, 0.85, by = 0.25)
a <- matrix(rep(c(5, 13, 50), each = 3) / 3, byrow = TRUE, nrow = 3)

# pdf(paste0(pth, "dirichlet.pdf"), width = 6, height = 6)
for (i in seq_len(nrow(a))) {
  mix <- dirichletmix(a = a[i, ])
  data <- genemix(500, mix)
  plot(mix, data, prop = TRUE, axl = axl, lv = lv)
}
# dev.off()
