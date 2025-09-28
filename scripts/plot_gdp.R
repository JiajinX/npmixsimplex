invisible(lapply(dir("R", full.names = TRUE), source))
pth <- "results/figures/5.1-"
axl <- 500

gdpraw <- readRDS("data/processed/gdp2023")
gdp <- merge_gdp(gdpraw, D = 3)
gdp[, -1] <- gdp[, -1] / rowSums(gdp[, -1])
data <- as.matrix(gdp[, -1]) # Three-part gdp compositions
# Data excluding zeros
data <- data[rowSums(data == 0) == 0, ] # Excluding zeros
# Model fitting
set.seed(10)
f_dircv <- dmdeCV(data, bwf = 5:1 / 5, verbose = 1)
f_cvkld <- f_dircv$mix_kl$mix
# Density by grid search
x1 <- seq(0, 0.1, len = 100)
x2 <- seq(0, 1, len = 1000)
gd <- expand.grid(x1, x2)
gd <- cbind(gd, 1 - rowSums(gd))
gd <- gd[rowSums(gd >= 0) == 3, ] # Grid points on the simplex
dbound <- predmix(gd, f_cvkld) # Densities of the grid points
maxbp <- gd[which.max(dbound), ] # Grid point with max density

# Scatter plot of GDP compositions
# pdf(paste0(pth, "GDPdot.pdf"), width = 6, height = 5)
simplexpanel(lab = list("Ag", "Se", "In"), labcex = 1.3)
points_simplex(gdp[, -1], pch = 20,
               col = transcols("#8C8C8C", 0.5), add = TRUE, cex = 1)
# dev.off()

# Density contour plot for GDP compositions (without 0) by Dir-CVKLD
# pdf(paste0(pth, "GDPden.pdf"), width = 6, height = 5)
lv <- seq(0.1, 0.7, length = 4)
simplexpanel(lab = list("Ag", "Se", "In"), labcex = 1.3)
plot(f_cvkld, data = data, prop = TRUE,
     point = FALSE, mode = FALSE, drawlabels = TRUE,
     axl = axl, lv = lv, add = TRUE)
points_simplex(attr(f_cvkld, "mode"), pch = 16,
               col = 2, add = TRUE, cex = 0.7)
points_simplex(maxbp, pch = 19,
               col = transcols("#3A5FCD", 0.7), add = TRUE, cex = 0.8)
# dev.off()

# Classification into two groups
gdp_income <- cbind(gdp, Income = gdpraw$Income)
gdpl <- as.matrix(gdp_income[gdp_income$Income %in% c("L", "LM"), 2:4])
gdph <- as.matrix(gdp_income[gdp_income$Income %in% c("UM", "H"), 2:4])
nl <- nrow(gdpl)
nh <- nrow(gdph)
prior <- c(nl, nh) / (nl + nh)
# Fitting models of each class
mdl <- dmdeCV(gdpl, bwf = 10:1 / 10, verbose = 1)$mix_kl$mix
mdh <- dmdeCV(gdph, bwf = 10:1 / 10, verbose = 1)$mix_kl$mix
# pdf(paste0(pth, "GDPcls.pdf"), width = 6, height = 5)
plev = seq(0.9, 0.99, len = 3)
cls_income(gdp_income, pch = 20, cex = 0.8, alpha = 0.5)
plot(mdl, point = FALSE, mode = FALSE, drawlabels = TRUE, plev = plev,
     col = transcols("#756142", 0.8), axl = axl, add = TRUE, labcex = 0.8,
     scale = prior[1])
plot(mdh, point = FALSE, mode = FALSE, drawlabels = TRUE, plev = plev,
     col = transcols("#454975", 0.8), axl = axl, add = TRUE, labcex = 0.8,
     scale = prior[2])
contour_bondary(mdl, mdh, prior = prior, add = TRUE,
                ax1 = c(0.4, 2 / sqrt(3) - 0.1), ax2 = c(0, 0.15),
                lty = 2, lwd = 2)
# dev.off()




