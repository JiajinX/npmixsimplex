invisible(lapply(dir("R", full.names = TRUE), source))

no <- 1 # Simulation No.: 1~12
type <- "chacon" # Type of simulation models: chacon or dirichlet
n <- 100 # Set generated sample size: 100 or 500
D <- 3 # Dimension: 3 or 6
rep <- 3 # Repetition: 1000
simu <- ms_model(type = type, no = no, n = n, D = D, rep = rep)
miae <- simu$miae
miseR <- simu$miaeR
kld <- simu$kld
time <- simu$time
# All the simulation results are saved in path: "results/raw/simulation"

# Summary all simulation results
pth <- "results/raw/simulation/"
simufiles <- list.files(pth, full.names = FALSE)
types <- c("chacon", "dirichlet")
loss <- c("miae", "miseR", "kld")
tb <- expand.grid(types = types, loss = loss, stringsAsFactors = FALSE)
for (i in seq_len(nrow(tb))) {
  print(i)
  key <- tb[i, ]
  files <- grep(paste0("^", key[1], ".*", key[2], "$"), simufiles, value = TRUE)
  tbre <- lapply(files, function(file) {
    f <- readRDS(paste0(pth, file))
    fmean <- t(summaryre(f))
    fse <- t(summaryre(f, type = "se"))
    ds <- as.numeric(regmatches(file, gregexpr("\\d+", file))[[1]])
    fmean <- cbind(ds[1] - 1, ds[2], fmean)
    fse <- cbind(ds[1] - 1, ds[2], fse)
    list(mean = fmean, se = fse)
  })
  mean <- do.call(rbind, lapply(tbre, `[[`, "mean"))
  se <- do.call(rbind, lapply(tbre, `[[`, "se"))
  colnames(mean)[1:2] <- colnames(se)[1:2] <- c("d", "n")
  tbname <- paste0("results/tables/simulation/", key[1], "_", key[2])
  if (FALSE) write.csv(mean, paste0(tbname, "_mean.csv"))
  if (FALSE) write.csv(se, paste0(tbname, "_se.csv"))
}
