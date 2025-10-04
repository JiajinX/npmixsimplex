# ilr transformation of the compositional data
ilrtra <- function(x) {
  # if (any(x == 0)) warning("The ilr transformation can not handel 0.")
  y <- compositions::ilr(x)
  attributes(y) <- NULL
  dim(y) <- dim(x) - c(0, 1)
  y
}

# Generate sample from Gaussian or Dirichlet mixture
genemix <- function(n, mix) {
  cls <- class(mix)
  if (cls == "dirichletmix") { # For Dirichlet
    x <- rdirichletmix(n, mix)
    y <- ilrtra(x)
  } else { # For Gaussian mixture (also for kde as special case)
    y <- if (cls == "kde") ks::rkde(n, mix) else rmvnormix(n, mix)
    x <- unclass(compositions::ilrInv(y))
  }
  jacb <- exp(rowSums(log(x))) * sqrt(ncol(x))
  attr(x, "ilr") <- y
  attr(x, "jacb") <- jacb # no 0 values
  x # with ilr and jacobian
}

# Predict density given a Gaussian or Dirichlet mixture
predmix <- function(data, mix, log = FALSE) {
  if (is.vector(data)) data <- t(data)
  cls <- class(mix)
  ilrdata <- attr(data, "ilr")
  jacb <- attr(data, "jacb")
  if (is.null(ilrdata)) ilrdata <- ilrtra(data)
  logjacb <- if (is.null(jacb)) {
    rowSums(log(data)) + 1 / 2 * log(ncol(data))
  } else log(jacb) # log Jacobian

  if (cls == "dirichletmix") {
    logd <- ddirichletmix(data, mix, log = TRUE)
    logdeuc <- logd + logjacb
    if (log) {
      d <- logd
      attr(d, "euc") <- logdeuc
    } else {
      d <- exp(logd)
      attr(d, "euc") <- exp(logdeuc)
    }
  } else if (cls == "kde") { # ilr-KDE
    deuc <- predict(mix, x = ilrdata)
    deuc[deuc == 0] <- deuc[deuc == 0] + 1e-16 # Avoid log(density) == -Inf
    logdeuc <- log(deuc)
    logd <- logdeuc - logjacb
    d <- exp(logd)
    if (log) {
      d <- logd
      attr(d, "euc") <- logdeuc
    } else attr(d, "euc") <- deuc
  } else { # ilr-MDE
    logdeuc <- dmvnormix(ilrdata, mix, log = TRUE)
    logd <- logdeuc - logjacb
    d <- exp(logd)
    if (log) {
      d <- logd
      attr(d, "euc") <- logdeuc
    } else attr(d, "euc") <- exp(logdeuc)
  }
  d
}

genest <- function(data,
                   est = c("ilr-PI", "ilr-CV", "ilr-MDE",
                           "Dir-KDE", "Dir-AIC",
                           "Dir-CV", "Dir-CVISE", "Dir-CVKLD"),
                   bwf = 10:1 / 10) {
  switch(est,
         "ilr-PI" = ks::kde(data),
         "ilr-CV" = ks::kde(data, H = Hlscv(data)),
         "ilr-MDE" = npden::mmde(data, bwf = bwf)$mix,
         "Dir-KDE" = dirichletmix(a = data / hkdeDir(data) + 1),
         "Dir-AIC" = dmde(data, bwf = bwf)$mix,
         "Dir-CV" = dmdeCV(data),
         "Dir-CVISE" = dmdeCV(data)$mix_ise$mix,
         "Dir-CVKLD" = dmdeCV(data)$mix_kl$mix)
}

art0 <- function(x, threshold = 0.005) {
  x[x < threshold] <- 0
  x[x >= threshold & x < threshold * 3] <- threshold * 2
  x <- x / rowSums(x)
  attr(x, "jacb") <- attr(x, "ilr") <- NULL
  x
}

# Shift boundary data
shift_achson <- function(x, round_err = NULL) {
  x <- as.matrix(x)
  if (is.null(round_err)) round_err <- min(x[x > 0])
  if (length(round_err) != 1) stop("Wrong parameter.")
  D <- ncol(x)
  n0 <- rowSums(x == 0)
  uniqn0 <- unique(n0)
  uniqn0 <- uniqn0[uniqn0 != 0]
  for (c in uniqn0) {
    idx <- n0 == c
    temp <- x[idx, , drop = FALSE]
    temp[temp != 0] <- temp[temp != 0] - round_err * c * (c + 1) / D^2
    temp[temp == 0] <- round_err * (c + 1) * (D - c) / D^2
    x[idx, ] <- temp
  }
  x
}

shift0 <- function(data,
                   type = c("achson", "add", "mult", "lrEM", "lrDA"),
                   dl = NULL) {
  data <- as.matrix(data)
  data <- data / rowSums(data)
  if (is.null(dl)) {
    if (type == "achson") {
      dl <- min(data[data > 0])
    } else dl <- apply(data, 2, function(col) min(col[col > 0]))
  }
  shift <- switch(type,
    # Aitchison Method: round_error set by user
    "achson" = shift_achson(data, round_err = dl),
    "add" = { # Additive Replacement
      shift <- compositions::zeroreplace(data, d = dl)
      attr(shift, "Losts") <- NULL
      shift
    }, # 2/3 min
    # Multiplicative Replacement
    "mult" = zCompositions::multRepl(data, label = 0, dl = dl),
    "lrEM" = zCompositions::lrEM(data, label = 0, dl = dl),
    "lrDA" = zCompositions::lrDA(data, label = 0, dl = dl)
  )
  shift <- as.matrix(shift / rowSums(shift))
  attr(shift, "ilr") <- ilrtra(shift)
  attr(shift, "jacob") <- exp(rowSums(log(shift))) * sqrt(ncol(shift))
  shift
}

rMsr <- function(model, test, n = 1e4) {
  samp <- genemix(n, model)
  term1 <- predmix(samp, model) # g^2
  term2 <- predmix(test, model) # g*f ???? What if test contrain 0

  iseS <- mean(term1) - 2 * mean(term2) # Simplex ISE

  eucterm1 <- attr(term1, "euc")
  eucterm2 <- attr(term2, "euc")
  eucterm2[eucterm2 == 0] <- NA # Ignored boundary observations on simplex
  if (any(is.na(eucterm2))) {
    warning("Observations on the simplex boundary are ignored
             in the calculation of ISER.")
  }
  iseR <- mean(eucterm1) - 2 * mean(eucterm2, na.rm = TRUE)

  term2[term2 == 0] <- NA # Ignored observations with 0 densities
  if (any(is.na(term2))) {
    warning("Observations with zero density values are ignored
             in the calculation of KLD.")
  }
  kld <- -mean(log(term2), na.rm = TRUE)

  return(list(iseR = iseR, iseS = iseS, kld = kld))
}

CVreal <- function(data, mthd, k = 10) {
  nmthd <- length(mthd)
  shift <- nmthd == 7
  miseR <- matrix(nrow = k, ncol = length(mthd))
  colnames(miseR) <- mthd
  if (shift) { # For shift data, we use ilr for Gaussian-based methods
    kld <- miseR
    idata <- attr(data, "ilr")
    jacb <- attr(data, "jacb")
  } else miseS <- miseR
  fold <- sample(rep(1:k, length.out = nrow(data)))
  for (i in seq_len(k)) {
    cat("# Fold = ", i, "\n")
    train <- data[fold != i, ]
    test <- data[fold == i, ]
    if (shift) {
      itrain <- idata[fold != i, ]
      attr(test, "ilr") <- idata[fold == i, ] # ilr-test
      attr(test, "jacb") <- jacb[fold == i] # test jacobian
      f_plug <- kde(itrain)
      f_cv <- kde(itrain, H = Hlscv(unique(itrain)))
      f_imde <- mmde(itrain, 10:1 / 10)$mix
    }
    f_kdir <- dirichletmix(a = train / hkdeDir(train) + 1)
    f_dmde <- dmde(train, 10:1 / 10, window = 3,
                   ictol = 20, criterion = "AIC")$mix
    recv <- dmdeCV(train)
    f_cvise <- recv$mix_ise$mix
    f_cvkl <- recv$mix_kl$mix

    ftotal <- if (shift) {
      list(f_plug = f_plug, f_cv = f_cv, f_kdir = f_kdir,
           f_imde = f_imde, f_dmde = f_dmde,
           f_cvise = f_cvise, f_cvkl = f_cvkl)
    } else {
      list(f_kdir = f_kdir, f_dmde = f_dmde, f_cvise = f_cvise, f_cvkl = f_cvkl)
    }

    ms <- lapply(ftotal, function(f) rMsr(f, test, shift))
    miseR[i, ] <- sapply(ms, `[[`, "iseR")
    if (shift) {
      kld[i, ] <- sapply(ms, `[[`, "kld")
    } else miseS[i, ] <- sapply(ms, `[[`, "iseS")
  }
  if (shift) {
    return(list(miseR = colMeans(miseR), kld = colMeans(kld)))
  } else return(list(miseR = colMeans(miseR), miseS = colMeans(miseS)))
}

repCVreal <- function(data, shift = NULL, k = 10, rep = 1) {
  data <- as.matrix(data)
  data <- data / rowSums(data)
  mthd <- c("ilr-PI", "ilr-CV", "Dir-KDE",
            "ilr-MDE", "Dir-MDE", "Dir-CVISE", "Dir-CVKL")
  miseR <- matrix(nrow = rep, ncol = length(mthd))
  colnames(miseR) <- mthd
  org <- is.null(shift)
  if (!org) {
    data <- shift0(data, type = shift)
    attr(data, "ilr") <- ilrtra(data) # ilr
    attr(data, "jacb") <- exp(rowSums(log(data))) * sqrt(ncol(data)) # Jacobian
    kld <- miseR
  } else {
    mthd <- mthd[c(3, 5:7)]
    miseS <- miseR <- miseR[, c(3, 5:7), drop = FALSE]
  }
  for (i in seq_len(rep)) {
    cat("Rep = ", i, "\n")
    ms <- CVreal(data, mthd, k)
    miseR[i, ] <- ms$miseR
    if (!org) kld[i, ] <- ms$kld else miseS[i, ] <- ms$miseS
  }
  if (!org) {
    return(list(miseR = miseR, kld = kld))
  } else return(list(miseR = miseR, miseS = miseS))
}

CLSfold <- function(train, test, bw, tol = 1e-4, by = NULL, verbose) {
  bw <- unique(sort(bw, decreasing = TRUE))
  if (is.null(by)) by <- which.min(colSums(apply(train, 2, duplicated)))
  y <- attr(train, "y")
  uniy <- unique(y)
  trainlist <- lapply(uniy, function(i) train[y == i, , drop = FALSE])
  names(trainlist) <- uniy
  prior <- as.numeric(prop.table(table(y)))
  sapply(bw, function(h) {
    cat("bandwidth:", h, "\n")
    mixlist <- lapply(trainlist, function(itrain) {
      cat("+")
      cat(names(itrain))
      cnmdm(itrain, h = h, by = by, tol = tol)$mix
    }) # Fit models of each class data with same bandwidth h
    cat("\n")
    dm <- sapply(mixlist, function(mix) predmix(test, mix))
    post <- prior * t(dm)
    cls_est <- uniy[apply(post, 2, which.max)]
    mean(cls_est == attr(test, "y"))
  })
}

dmdeCLS <- function(data, y, k = 5, bwf = 10:1 / 10, verbose = 1) {
  D <- ncol(data)
  datalist <- lapply(unique(y), function(i) data[y == i, , drop = FALSE])
  names(datalist) <- unique(y)
  ia0vec <- sapply(datalist, function(i) initial_a0(i))
  ia0 <- max(ia0vec) # min?
  bw <- 1 / ((ia0 + 1) / bwf^2 - 1 - D)
  names(bw) <- round(bwf, 2)
  print(bw)
  fold <- sample(rep(1:k, length.out = nrow(data))) # Split data into k-folds
  acc <- matrix(nrow = k, ncol = length(bwf))
  colnames(acc) <- round(bw, 3)
  for (i in seq_len(k)) {
    if (verbose > 0) cat("# K-fold:", i, "\n")
    train <- data[fold != i, ]
    attr(train, "y") <- y[fold != i]
    test <- data[fold == i, ]
    attr(test, "y") <- y[fold == i]
    acc[i, ] <- CLSfold(train, test, bw = bw, verbose = verbose)
    print(acc)
  }
  h <- bw[which.max(colMeans(acc))]
  lapply(datalist, function(i) cnmdm(i, h = h)$mix)
}

classify <- function(method, train, target_train, test, target_test) {
  unitarget <- unique(c(target_train, target_test))
  if (method == "Dir-CV") { # Dir-CV will return 2 models
    dise <- dkl <- matrix(nrow = length(target_test), ncol = length(unitarget))
  } else d <- matrix(nrow = length(target_test), ncol = length(unitarget))
  if (method == "Dir-CLS") {
    print("Dir-CLS")
    models <- dmdeCLS(train, target_train)
    d <- sapply(models, function(mix) predmix(test, mix))
  } else {
    for (i in seq_along(unitarget)) {
      print(paste0(method, " --> ", unitarget[i]))
      x <- train[target_train == unitarget[i], ]
      model <- switch(method,
                      "ilr-PI" = kde(x),
                      "ilr-CV" = kde(x, H = Hlscv(unique(x))),
                      "Dir-KDE" = dirichletmix(a = x / hkdeDir(x) + 1),
                      "ilr-MDE" = mmde(x, 10:1 / 10)$mix,
                      "Dir-MDE" = dmde(x, 10:1 / 10, window = 3,
                                       ictol = 20, criterion = "AIC")$mix,
                      "Dir-CV" = dmdeCV(x))
      if (method == "Dir-CV") {
        dise[, i] <- predmix(test, model$mix_ise$mix)
        dkl[, i] <- predmix(test, model$mix_kl$mix)
      } else d[, i] <- predmix(test, model)
    }
  }
  prior <- as.numeric(prop.table(table(target_train)))
  if (method == "Dir-CV") {
    post_ise <- prior * t(dise)
    cls_ise <- unitarget[apply(post_ise, 2, which.max)]
    acc_ise <- mean(cls_ise == target_test)
    attr(acc_ise, "conf_mat") <- table(Predicted = cls_ise,
                                       Actual = target_test)
    post_kl <- prior * t(dkl)
    cls_kl <- unitarget[apply(post_kl, 2, which.max)]
    acc_kl <- mean(cls_kl == target_test)
    attr(acc_kl, "conf_mat") <- table(Predicted = cls_kl,
                                      Actual = target_test)
    return(c(acc_ise, acc_kl))
  } else {
    post <- prior * t(d)
    cls_est <- unitarget[apply(post, 2, which.max)]
    acc <- mean(cls_est == target_test)
    attr(acc, "conf_mat") <- table(Predicted = cls_est, Actual = target_test)
    return(acc)
  }
}

CVclassify <- function(data, target, mthd, k = 10) {
  isilr <- mthd %in% c("ilr-PI", "ilr-CV", "ilr-MDE")
  acc <- matrix(nrow = k, ncol = length(mthd))
  colnames(acc) <- mthd
  if (any(isilr)) { # For shift data, we use ilr for Gaussian-based methods
    idata <- attr(data, "ilr")
    jacb <- attr(data, "jacb")
  }
  fold <- sample(rep(1:k, length.out = nrow(data)))
  for (i in seq_len(k)) {
    cat("Fold = ", i, "\n")
    train <- data[fold != i, ]
    test <- data[fold == i, ]
    target_train <- target[fold != i]
    target_test <- target[fold == i]
    if (any(isilr)) {
      itrain <- idata[fold != i, ]
      attr(test, "ilr") <- idata[fold == i, ] # ilr-test
      attr(test, "jacb") <- jacb[fold == i] # test jacobian
      clsmthd <- mthd[isilr]
      for (j in clsmthd) {
        acc[i, j] <- classify(j, itrain, target_train, test, target_test)
      }
    }
    clsmthd <- mthd[!isilr]
    DirCV <- FALSE
    for (j in clsmthd) {
      if (j %in% c("Dir-CVISE", "Dir-CVKL")) {
        DirCV <- TRUE
      } else acc[i, j] <- classify(j, train, target_train, test, target_test)
    }
    if (DirCV) acc[i, c("Dir-CVISE", "Dir-CVKL")] <-
      classify("Dir-CV", train, target_train, test, target_test)
  }
  colMeans(acc)
}

repCVclassify <- function(data, target, mthd = NULL,
                          shift = NULL, k = 10, rep = 1) {
  data <- as.matrix(data)
  data <- data / rowSums(data)
  if (is.null(mthd)) mthd <- c("ilr-PI", "ilr-CV", "ilr-MDE",
                               "Dir-KDE", "Dir-MDE", "Dir-CVISE",
                               "Dir-CVKL", "Dir-CLS")
  acc <- matrix(nrow = rep, ncol = length(mthd))
  colnames(acc) <- mthd
  if (!is.null(shift)) data <- shift0(data, type = shift)
  for (i in seq_len(rep)) {
    cat("======= Rep", i, "======= \n")
    acc[i, ] <- CVclassify(data, target, mthd, k)
  }
  acc
}
