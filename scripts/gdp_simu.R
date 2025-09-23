invisible(lapply(dir("R", full.names = TRUE), source))

D <- 3 # Dimension: 3 or 6
ref <- "ilr-PI" # Reference density: ilr-PI or Dir-CVKLD
gdp <- readRDS("data/processed/gdp2023")
gdp <- merge_gdp(gdp, D = D)
data <- as.matrix(gdp[, -1] / rowSums(gdp[, -1]))
mix0 <- est_refd(ref, data)

threshold <- 0.01
shift <- "no" # no or zero replacement method: achson, mult, lrEM
rep <- 5 # Repetition: 1000
re <- ms_gdp(mix0, n = nrow(gdp),
             threshold = threshold, shift = shift, rep = rep)
miae <- re$miae
miseR <- re$miseR
kld <- re$kld
# All the simulation results are saved in path: "results/raw/gdp"

# Summary all simulation results
pth <- "results/raw/gdp/"
simufiles <- list.files(pth, full.names = FALSE)
types <- c("ilr-PI", "Dir-CVKLD")
for (ref in types) {
  files <- grep(ref, simufiles, value = TRUE)
  refre <- lapply(files, function(file) {
    f <- readRDS(paste0(pth, file))
    fmean <- colMeans(f)
    fse <- apply(f, 2, sd) / sqrt(nrow(f))
    info <- strsplit(file, "_")[[1]]
    d <- as.numeric(gsub("\\D", "", info[1])) - 1
    loss <- toupper(info[4])
    if (info[3] == "no") {
      fmean <- c(rep(NA, 3), fmean)
      fse <- c(rep(NA, 3), fse)
      names(fmean)[1:3] <- names(fse)[1:3] <- c("ilr-PI", "ilr-CV", "ilr-MDE")
      info[3] <- "zero"
    }
    fmean <- data.frame(d = d, loss = loss, zerorep = info[3], t(fmean))
    fse <- data.frame(d = d, loss = loss, zerorep = info[3], t(fse))
    list(mean = fmean, se = fse)
  })
  mean <- do.call(rbind, lapply(refre, `[[`, "mean"))
  se <- do.call(rbind, lapply(refre, `[[`, "se"))
  colnames(mean) <- colnames(se) <- gsub("\\.", "-", colnames(mean))
  colnames(mean)[8] <- colnames(se)[8] <- "Dir-AIC"
  filename <- paste0("results/tables/gdp/", ref)
  if (FALSE) write.csv(mean, paste0(filename, "_mean.csv"), row.names = FALSE)
  if (FALSE) write.csv(se, paste0(filename, "_se.csv"), row.names = FALSE)
}


