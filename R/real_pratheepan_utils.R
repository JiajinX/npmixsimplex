# Get RGB or gray values of all pixels in an image
getpx <- function(img) {
  info <- image_info(img)
  format <- info$format
  channels <- ifelse(format == "JPEG", "rgb", "gray")
  imgpx <- magick::image_data(img, channels = channels)
  if (channels == "rgb") {
    RGB <- cbind(as.integer(imgpx[1, , ]),
                 as.integer(imgpx[2, , ]),
                 as.integer(imgpx[3, , ]))
    RGB[RGB == 0] <- 1
    return(RGB)
  } else as.integer(imgpx[1, , ]) == 255
}

# Get the size, rgb, grey values of an image
getimgpx <- function(no) {
  info <- image_info(fm[no])
  list(name = fmname[no], size = c(info$width, info$height),
       rgb = getpx(fm[no]), msk = getpx(fmt[no]))
}

# Show the sub-pixels of an image given a mask
showskin <- function(img, iskin) {
  if (class(img) != "magick-image") img <- imgpx2img(img)
  if (is.logical(iskin)) iskin <- ifelse(iskin, 255, 0)
  img_rbga <- magick::image_data(img, channels = "rgba")
  img_rbga[4, , ] <- as.raw(iskin)
  plot(as.raster(aperm(img_rbga)))
}

# Show the mask
showmsk <- function(imgpx, iskin) {
  dim(iskin) <- imgpx$size
  mask <- magick::image_read(as.raster(t(iskin)))
  plot(mask)
  mask
}

showpx <- function(imgpx, iskin = NULL, type = c("org", "mask", "skin")) {
  type <- match.arg(type)
  if (is.null(iskin)) iskin <- imgpx$msk
  if (type == "org") {
    plot(imgpx2img(imgpx))
  } else if (type == "skin") {
    showskin(imgpx, iskin)
  } else showmsk(imgpx, iskin)
}

imgpx2img <- function(imgpx) {
  fm[which(fmname == imgpx$name)]
}

imgpx2msk <- function(imgpx) {
  fmt[which(fmname == imgpx$name)]
}

# Get the sample index from each image with proportion
sampleidx <- function(pxlist, n) {
  en <- sapply(pxlist, nrow)
  spen <- round(en / sum(en) * n)
  Map(sample, en, spen)
}

# Get the sample pixels given by the sample index
samplepx <- function(pxlist, idx, ex = FALSE) {
  do.call(rbind, Map(function(x, i) {
    if (ex) x[-i, , drop = FALSE] else x[i, , drop = FALSE]
  }, pxlist, idx))
}

to_comp <- function(data) data / rowSums(data)

# Get the confusion matrix given by the specified tuning parameters
tuncm <- function(gamma, fskin, foth, dqt, data, tmsk) {
  dskin <- predmix(data, fskin$mix)
  doth <- predmix(data, foth$mix)
  iskin1 <- gamma * dskin > doth
  iskin2 <- dskin > dqt
  table(tmsk, iskin1 & iskin2)
}

cminfo <- function(cm) {
  TN <- cm[1, 1]
  FN <- cm[2, 1]
  FP <- cm[1, 2]
  TP <- cm[2, 2]
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  F1 <- 2 * precision * recall / (precision + recall)
  F2 <- 5 * precision * recall / (4 * precision + recall)
  c(TN, FN, FP, TP, precision, recall, F1, F2)
}

getskin <- function(imgpx, skinmd, othmd, gamma, dqt) {
  dskin <- predmix(to_comp(imgpx$rgb), skinmd)
  doth <- predmix(to_comp(imgpx$rgb), othmd)
  mask_cls <- dskin * gamma > doth
  mask_qt <- dskin > dqt
  mask_cls & mask_qt
}

testimg <- function(imgpx, mask) {
  cm <- table(imgpx$msk, mask)
  cminfo(cm)
}
