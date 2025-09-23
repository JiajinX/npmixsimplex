# log format -> wide format
clean_gdp <- function(gdp) {
  indicator <- unique(gdp$IndicatorName)[c(8:14)]
  gdp <- gdp[gdp$IndicatorName %in% indicator, ] # Select column
  gdp_wide <- tidyr::pivot_wider(gdp,
                                 names_from = IndicatorName,
                                 values_from = `2023`)
  colnames(gdp_wide) <- c("Area", "A2B", "C2E", "D", "F", "G2H", "I", "J2P")
  gdp_wide <- gdp_wide[rowSums(is.na(gdp_wide)) != length(indicator), ]
  gdp_wide[is.na(gdp_wide)] <- 0
  gdp_wide
}

# Merge gdp into d = 3, 6
merge_gdp <- function(gdp, D = c(3, 6)) {
  switch(as.character(D),
         "6" = gdp[, c("Area", "A2B", "C2E", "F", "G2H", "I", "J2P")],
         "3" = {
           # (Industry) A-B: Primary, C-F: Secondary, G-U: Tertiary
           gdp$C2F <- gdp$C2E + gdp$F
           gdp$G2P <- gdp$G2H + gdp$I + gdp$J2P
           gdp[, c("Area", "A2B", "C2F", "G2P")]
         },
         stop("D must be 3, 6."))
}

# Get reference density for GDP compositions
est_refd <- function(est = c("ilr-PI", "ilr-CV", "ilr-MDE",
                             "Dir-KDE", "Dir-AIC", "Dir-CVISE", "Dir-CVKLD"),
                     data) {
  if (est %in% c("ilr-PI", "ilr-CV", "ilr-MDE")) {
    # Gaussian-based estimator
    data <- data[rowSums(data == 0) == 0, ] # Excluding zeros
    idata <- ilrtra(data) # ilr transformation
    ref <- genest(idata, est)
  } else ref <- genest(data, est) # Dirichlet-based estimator
  ref
}

# Simulation results from the parametric bootstrap of GDP compositions
# mix0: reference density
# threshold: introduce zero
# shift: zero replacement method used for generated samples
ms_gdp <- function(mix0, n, threshold, shift, samp = 1e4, rep = 1) {
  cn <- c("ilr-PI", "ilr-CV", "ilr-MDE",
          "Dir-KDE", "Dir-AIC", "Dir-CVISE", "Dir-CVKLD")
  miae <- matrix(nrow = rep, ncol = length(cn))
  colnames(miae) <- cn
  if (shift == "no") miae <- miae[, 4:7, drop = FALSE]
  # Data with zeros can only allow Dir-based estimators
  kld <- miseR <- miae
  mix0_ie <- IE_samp(mix0, n = samp)
  mix0_kld <- KLD_samp(mix0, n = samp)
  for (i in 1:rep) {
    print(i)
    # Generate data
    data <- genemix(n, mix0)
    sdata <- art0(data, threshold = threshold) # Data contains zeros
    cat(sum(rowSums(sdata == 0)) / nrow(sdata), "\n")
    # Proportion of observations contain 0
    if (shift != "no") { # If zero replacement
      sdata <- shift0(sdata, type = shift) # Zero-replaced data
      isdata <- attr(sdata, "ilr")

      f_plug <- ks::kde(isdata) # KDE-ilr: Plug_in
      f_cv <- ks::kde(isdata, H = Hlscv(isdata)) # KDE-ilr: CV
      f_imde <- npden::mmde(isdata, 10:1 / 10)$mix # MDE-ilr
    }

    f_kdir <- dirichletmix(a = sdata / hkdeDir(sdata) + 1) # Dir-KDE
    f_dmde <- dmde(sdata, 10:1 / 10)$mix # Dir-AIC
    recv <- dmdeCV(sdata) # Dir-CV
    f_cvise <- recv$mix_ise$mix
    f_cvkl <- recv$mix_kl$mix

    ftotal <- if (shift != "no") {
      list(f_plug = f_plug, f_cv = f_cv, f_imde = f_imde,
           f_kdir = f_kdir, f_dmde = f_dmde,
           f_cvise = f_cvise, f_cvkl = f_cvkl)
    } else list(f_kdir = f_kdir, f_dmde = f_dmde,
                f_cvise = f_cvise, f_cvkl = f_cvkl)

    ielst <- lapply(ftotal, function(f) IE(mix0_ie, f))
    miae[i, ] <- sapply(ielst, `[[`, "iae") # IAE,
    miseR[i, ] <- sapply(ielst, `[[`, "iseR") # ISE on real space
    kld[i, ] <- sapply(ftotal, function(f) KLD(mix0_kld, f))
  }

  list(miae = miae, miseR = miseR, kld = kld)
}

# Plot gdp compositions colored by income group
cls_income <- function(gdp, pch = 1, cex = 1, alpha = 1) {
  groups <- unique(gdp$Income)
  gpnum <- length(groups)
  cols <- c("#D5695D", "#65A479", "#5D8CA8", "#D3BA68", "#333333")
  cols <- transcols(cols, alpha)
  mtx <- as.matrix(gdp[, c("A2B", "C2F", "G2P")])
  mtx <- mtx / rowSums(mtx)
  simplexpanel(lab = list("Ag", "Se", "In"))
  for (i in seq_along(groups)) {
    idx <- if (is.na(groups[i])) is.na(gdp$Income) else gdp$Income == groups[i]
    points_simplex(mtx[idx, ], col = cols[i], cex = cex, add = TRUE, pch = pch)
  }
  legend(x = 0.9, y = 1.1, title = "Income:", title.adj = 0.12,
         legend = c("Low", "Lower middle", "Upper middle", "High", "Undefined"),
         col = cols[c(1, 4, 2, 3, 5)], pch = pch, cex = 1)
}

