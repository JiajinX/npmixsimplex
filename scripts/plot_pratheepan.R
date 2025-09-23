invisible(lapply(dir("R", full.names = TRUE), source))
pth <- "results/figures/5.3-pratheepan/"
mdpth <- "models/pratheepan/"
datapth <- "data/processed/"
repth <- "results/raw/"
fm_col <- "data/raw/Pratheepan_Dataset/FamilyPhoto/"

traskin <- readRDS(paste0(datapth, "pratheepan_traskin"))
traoth <- readRDS(paste0(datapth, "pratheepan_traoth"))
gds <- readRDS(paste0(repth, "pratheepan_gridsearch")) # Grid search results
op <- unlist(gds[which.max(gds$F2), 1:4]) # Find the optimal tuning parameters
skinmd <- readRDS(paste0(mdpth, "skinmd", op[1]))
othmd <- readRDS(paste0(mdpth, "othmd", op[2]))
gam <- op[3]
qt <- op[4]
dqt <- quantile(predmix(to_comp(traskin), skinmd$mix), probs = qt)
# pdf(paste0(pth, "pratheepan.pdf"), height = 6, width = 6)
simplexpanel(lab = list("R", "G", "B"))
points_simplex(to_comp(traoth), col = transcols("#CCCCCC", 0.2),
               pch = 20, add = TRUE)
points_simplex(to_comp(traskin), col = transcols("#E9C39B", 0.2),
               pch = 20, add = TRUE)
contour_bondary(skinmd$mix, othmd$mix, prior = c(gam, 1), threshold = dqt,
                level = 1e-3, axl = 500, add = TRUE)
legend(x = 0.85, y = 1, legend = c("non-skin", "skin"),
       col = c("gray80", "#E9C39B"), pch = 20, cex = 1, pt.cex = 1.2)
# dev.off()

fm <- magick::image_read(list.files(fm_col, full.names = TRUE))
fmname <- gsub(".jpg", "", list.files(fm_col))
fmtest <- readRDS(paste0(datapth, "pratheepan_test"))
testre <- readRDS(paste0(repth, "pratheepan_test_results"))
idx <- c(4, 12, 10, 3) # idx = 1,...,14
for (i in idx) {
  print(i)
  test <- fmtest[[i]]
  premsk <- getskin(test, skinmd$mix, othmd$mix, gamma = gam, dqt = dqt)
  showpx(test)
  showpx(test, premsk, type = "skin")
  mask <- showpx(test, premsk, type = "mask")
  # magick::image_write(mask,
  #                     path = paste0(pth, "mask", i, ".png"), format = "png")
}
