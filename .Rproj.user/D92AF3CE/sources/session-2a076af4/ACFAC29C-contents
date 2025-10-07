invisible(lapply(dir("R", full.names = TRUE), source))
library(npden)

# Read data
train <- readRDS("data/raw/usps_train")
test <- readRDS("data/raw/usps_test")

# Transform to compositional data
ctrain <- transpx(train)
ctest <- transpx(test)
# Dimension reduce to d = 16
train16 <- mergepx(ctrain)
test16 <- mergepx(ctest)
# Fitting models
ests <- c("ilrMDE", "DirKDE", "DirAIC", "DirCV")
for (est in ests) { # Read models
  pth <- paste0("models/usps/", est)
  if (file.exists(pth)) {
    assign(est, readRDS(pth))
  } else { # If file not exists, fitting models
    assign(est, cls_digit(train16, est = est))
    saveRDS(get(est), pth)
  }
}
DirCVISE <- lapply(DirCV, function(x) x$mix_ise$mix)
DirCVKLD <- lapply(DirCV, function(x) x$mix_kl$mix)
est_model <- list(ilrMDE, DirKDE, DirAIC, DirCVISE, DirCVKLD)
estnames  <- c("ilrMDE", "DirKDE", "DirAIC", "DirCVISE", "DirCVKLD")
names(est_model) <- estnames

# Prediction
prior <- as.numeric(prop.table(table(train16[, 1])))
zrtest16 <- test16
zrtest16[, -1] <- shift0(test16[, -1], type = "mult")
# Test after multivariate zero replacement (for ilr-MDE)
clslab <- lapply(names(est_model), function(model) {
  tst <- if (model == "ilrMDE") zrtest16 else test16
  prelabel(est_model[[model]], tst, prior)
}) # Prediction labels
names(clslab) <- estnames
clsms <- lapply(clslab, function(i) evalcls(i, test16[, 1])) # Measurement
if (FALSE) write.csv(t(sapply(clsms, `[[`, "f1micro")),
                     "results/tables/usps.csv", row.names = FALSE)
# Significant test
compF1(test16[, 1], clslab$DirKDE, clslab$DirCVKLD, f1 = "micro")

