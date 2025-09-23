invisible(lapply(dir("R", full.names = TRUE), source))
library(magick)

# Path
fm_col <- "data/raw/Pratheepan_Dataset/FamilyPhoto/"
fm_msk <- "data/raw/Ground_Truth/GroundT_FamilyPhoto/"
# # Path + file name
fm <- image_read(list.files(fm_col, full.names = TRUE))
fmt <- image_read(list.files(fm_msk, full.names = TRUE))
# file name
fmname <- gsub(".jpg", "", list.files(fm_col))
# Split to 70% training and 30% test data
set.seed(422)
fm_trid <- sample(seq_along(fmname),  floor(length(fmname) * 0.7))
fm_tsid <- setdiff(seq_along(fmname), fm_trid) # Training and test index
# Get information of training and test images
fmtrain <- lapply(fm_trid, getimgpx)
fmtest <- lapply(fm_tsid, getimgpx)
if (FALSE) saveRDS(fmtest, "data/processed/pratheepan_test")
# Skin and non-skin pixels of training images
skinpx <- lapply(fmtrain, function(img) img$rgb[img$msk, ])
othpx <- lapply(fmtrain, function(img) img$rgb[!img$msk, ])
# Sampling from pixels
spsize <- 1e4
# Sampling index for each image in the skin and non-skin classes
traskin_idx <- sampleidx(skinpx, spsize)
traoth_idx <- sampleidx(othpx, spsize)
# Get the training pixels of skin class and non-skin classes
traskin <- samplepx(skinpx, traskin_idx)
traoth <- samplepx(othpx, traoth_idx)
if (FALSE) {
  saveRDS(traskin, "data/processed/pratheepan_traskin")
  saveRDS(traoth, "data/processed/pratheepan_traoth")
}

# Training models
bwf <- seq(0.4, 0.1, by = -0.1)
skinmdlst <- othmdlst <- vector("list", length(bwf))
names(skinmdlst) <- names(othmdlst) <- bwf
for (i in bwf) {
  i <- as.character(i)
  skinpth <- paste0("models/pratheepan/skinmd", i)
  othpth <- paste0("models/pratheepan/othmd", i)
  if (file.exists(skinpth)) {
    skinmdlst[[i]] <- readRDS(skinpth)
  } else {
    skinmdlst[[i]] <- cnmdm(to_comp(traskin), h = i)
    saveRDS(skinmdlst[[i]], skinpth)
  }
  if (file.exists(othpth)) {
    othmdlst[[i]] <- readRDS(othpth)
  } else {
    othmdlst[[i]] <- cnmdm(to_comp(traoth), h = i)
    saveRDS(othmdlst[[i]], othpth)
  }
}

# Select tuning parameters on validation set
pth <- "results/raw/pratheepan_gridsearch"
if (file.exists(pth)) {
  gds <- readRDS(pth)
} else {
  valskin <- samplepx(skinpx, traskin_idx, ex = TRUE) # Skin pixels
  valoth <- samplepx(othpx, traskin_idx, ex = TRUE) # Non-skin pixels
  val <- rbind(valskin, valoth)
  valmsk <- c(rep(TRUE, nrow(valskin)), rep(FALSE, nrow(valoth))) # Mask
  # Sampling from validation set
  valsize <- 1e5
  validx <- sample(seq_along(valmsk), valsize) # Sampling index
  valdata <- to_comp(val[validx, ])
  valmsk <- valmsk[validx]
  # Tuning parameter space
  skinbwf <- othbwf <- bwf
  vga <- seq(0.1, 2, by = 0.05) # Gamma
  vqt <- seq(0, 0.2, by = 0.05) # p
  tun <- expand.grid(skinbwf, othbwf, vga, vqt)
  colnames(tun) <- c("skinbwf", "othbwf", "vga", "vqt")
  cmlst <- list()
  case <- nrow(tun)
  for (idx in seq_len(case)) {
    print(idx)
    para <- tun[idx, ]
    skinmd <- skinmdlst[[as.character(para[1])]]
    othmd <- othmdlst[[as.character(para[2])]]
    gam <- para[[3]]
    dqt <- para[[4]]
    cmlst[[idx]] <- tuncm(gam, skinmd, othmd, dqt, valdata, valmsk)
    # Calculate confusion matrix for each combination of parameters
  }
  clsmetric <- t(sapply(cmlst, cminfo))
  colnames(clsmetric) <- c("TN", "FN", "FP", "TP",
                           "precesion", "recall", "F1", "F2")
  gds <- cbind(tun, clsmetric)
  saveRDS(gds, pth)
}

# Select a model with the optimal tuning parameters by maximize F2-score
op <- unlist(gds[which.max(gds$F2), 1:4])
skinmd <- readRDS(paste0("models/pratheepan/skinmd", op[1]))
othmd <- readRDS(paste0("models/pratheepan/othmd", op[2]))
gam <- op[3]
qt <- op[4]
dqt <- quantile(predmix(to_comp(traskin), skinmd$mix), probs = qt)

# Evaluate test images
pth <- "results/raw/pratheepan_test_results"
if (file.exists(pth)) {
  testre <- readRDS(pth)
} else {
  testre <- matrix(nrow = 0, ncol = 8)
  colnames(testre) <- c("TN", "FN", "FP", "TP",
                        "precision", "recall", "F1", "F2")
  for (test in fmtest) {
    print(test$name)
    premsk <- getskin(test, skinmd$mix, othmd$mix, gamma = gam, dqt = dqt)
    testre <- rbind(testre, testimg(test, premsk))
  }
  if (FALSE) saveRDS(testre, pth)
}

subtestre <- t(testre[c(4, 12, 10, 3), c("precision", "recall", "F2")])
colnames(subtestre) <- c("a", "b", "c", "d")
if (FALSE) write.csv(subtestre, "results/tables/pratheepan.csv")



