# Transform digit data to compositions
transpx <- function(digit) {
  if (is.vector(digit)) digit <- t(digit)
  digit[, -1] <- (digit[, -1] + 1) / 2
  digit[, -1] <- digit[, -1] / rowSums(digit[, -1])
  digit
}

# Plot each digit row
plot_digit <- function(digit) {
  x <- digit[-1]
  n <- sqrt(length(x))
  dim(x) <- c(n, n)
  image(x[, n:1], col = gray((32:0)/32), axes = FALSE, asp = 1)
  box()
  print(digit[1])
}

# 16*16 -> 16 features by 4*4 window: Average grey value of each window
mergepx <- function(digit, pexnum = 16) {
  x <- digit[, -1]
  colidx <- matrix(1:256, nrow = 16, ncol = 16, byrow = TRUE)
  if (pexnum == 16) {
    wind <- matrix(rep(rep(1:4, each = 4), 16),
                   nrow = 16, ncol = 16, byrow = TRUE) +
      rep(c(0:3) * 4, each = 4)
  } else if (pexnum == 4) {
    wind <- matrix(rep(rep(1:2, each = 8), 16),
                   nrow = 16, ncol = 16, byrow = TRUE) + rep(c(0,2), each = 8)
  }
  mx <- matrix(nrow = nrow(digit), ncol = pexnum)
  for (i in 1:pexnum) {
    idx <- colidx[wind == i]
    mx[, i] <- rowSums(x[, idx])
  }
  cbind(digit[, 1], mx)
}

# Fitting models for each class
cls_digit <- function(train,
                      est = c("ilr-MDE", "Dir-KDE", "Dir-AIC", "Dir-CV")) {
  modelst <- list()
  data <- train[, -1]
  bwf <- 10:1 / 10
  if (est == "ilr-MDE") {
    data <- shift0(data, type = "mult")
    data <- attr(data, "ilr")
    bwf <- 100:90 / 100
  }
  for (digit in 0:9) {
    cat("Digit:", digit, "\n")
    idx <- train[, 1] == digit
    subdata <- data[idx, ]
    modelst[[paste0("digit", digit)]] <- genest(subdata, est, bwf)
  }
  modelst
}

# Predict labels
prelabel <- function(modelst, test, prior) {
  d <- lapply(modelst, function(model) predmix(test[, -1], model))
  md <- do.call(cbind, d)
  post <- prior * t(md)
  apply(post, 2, which.max) - 1
}

# Calculate accuracy, F1-macro, confusion matrix
evalcls <- function(pred, true) {
  acc <- mean(true == pred)
  cm <- table(Predicted = pred, Actual = true)
  f1micro <- MLmetrics::F1_Score_micro(true, pred)
  f1macro <- MLmetrics::F1_Score_macro(true, pred)
  list(acc = acc, f1micro = f1micro, f1macro = f1macro, cm = cm)
}

# Compare F1-macro of based on prediction labels of two different models
compF1 <- function(true, pred1, pred2, f1 = c("micro", "macro"), t = 1e4) {
  f1 <- ifelse(f1 == "micro",
               MLmetrics::F1_Score_micro, MLmetrics::F1_Score_macro)
  len <- length(true)
  diff <- rep(NA, t)
  for (i in seq_len(t)) {
    print(i)
    idx <- sample(len, len, replace = TRUE)
    diff[i] <- f1(true[idx], pred1[idx]) - f1(true[idx], pred2[idx])
  }
  quantile(diff, c(0.05, 0.95))
}

usps_margin <- function(models, digits, data, comb, lv = 0.2, axl = 100,
                        cols = paletteer::paletteer_d("ggsci::category10_d3")) {
  axe <- sapply(comb, function(i) paste(i, collapse = ","))
  simplexpanel(lab = list(bquote(px[.(axe[3])]),
                          bquote(px[.(axe[1])]),
                          bquote(px[.(axe[2])])), labcex = 1.3)
  for (i in digits) {
    subdata <- data[data[, 1] == i, -1]
    model <- models[[paste0("digit", i)]]
    plot(model, subdata, x1 = comb[[2]], x2 = comb[[3]],
         mode = FALSE, prop = TRUE, lv = lv, col = cols[i + 1], add = TRUE,
         axl = axl, drawlabels = FALSE)
  }
  legend(x = 0.9, y = 1, legend = digits, col = cols[digits + 1],
         lty = 1, cex = 1)
}
