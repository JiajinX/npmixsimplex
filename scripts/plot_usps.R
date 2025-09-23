invisible(lapply(dir("R", full.names = TRUE), source))
pth <- "results/figures/5.2-"
axl <- 300

# Training data
train <- readRDS("data/raw/usps_train")
ctrain <- transpx(train)
train16 <- mergepx(ctrain)
# Model
DirCV <- readRDS("models/usps/DirCV")
DirCVKLD <- lapply(DirCV, function(x) x$mix_kl$mix)

# pdf(paste0(pth, "digit.pdf"), width = 8, height = 5)
usps_margin(DirCVKLD, digits = c(0,1,6,7,8,9), train16, axl = axl,
            comb = list(c(1,2,3,4,5,8),
                        c(6,7,10,11),
                        c(9,12,13,14,15,16)))
usps_margin(DirCVKLD, digits = c(2,5), train16, axl = axl,
            comb = list(c(5,6),
                        c(9,10),
                        c(1,2,3,4,7,8,11,13,14,15,16)))
usps_margin(DirCVKLD, digits = c(3, 4), train16, axl = axl,
            comb = list(c(1,2,5,9,13,14),
                        c(6,7,10,11),
                        c(3,4,8,12,15,16)))
# dev.off()
