## ============== Print and plot function of Dirichlet mixture ============== ##

# Print information of Dirichlet mixture of each iterations
print_verbose <- function(verbose, i, mix) {
  if (verbose >= 1) cat("\n## Iteration ", i, ":\nLog-likelihood = ",
                        as.character(signif(attr(mix, "ll"), 16)), "\n",
                        sep = "")
  if (verbose >= 2) {
    cat("Component proportions:\n")
    print(mix$pr)
  }
  if (verbose >= 3) {
    cat("Component alpha:\n")
    print(mix$alpha)
  }
  if (verbose >= 4) {
    cat("Component mode:\n")
    print(attr(mix, "mode"))
  }
  if (verbose >= 5) {
    cat("Component bandwith:\n")
    print(mix$h)
  }
}

# Print Dirichlet mixture after sorting components
print.dirichletmix <- function(mix) {
  mix <- sort(mix)
  class(mix) <- NULL
  print(mix)
}

# Create a simplex panel
simplexpanel <- function(lab = list(expression(x[1]),
                                    expression(x[2]),
                                    expression(x[3])),
                         right = FALSE, title = NULL, faxis = 100,
                         margin = 0.05, boundlty = 1, labcex = 1.3) {
  par(mar = rep(1, 4))
  liml <- 0.1
  labl <- 0.06
  sidel <- ifelse(right, 1, 2 / sqrt(3))
  tcex <- 0.8
  labcex <- labcex
  plot.new()
  plot.window(xlim = c(-liml, sidel + liml), ylim = c(-liml, 1 + liml),
              xaxs = "i", yaxs = "i", asp = 1)
  title(main = title)
  if (right) {
    l <- matrix(c(0, 0, 0, 1,
                  0, 0, 1, 0,
                  0, 1, 1, 0),
                nrow = 4)
    gap <- seq(0, 1, by = 0.1)
    text(1 + labl, 0, lab[[1]], cex = labcex)
    text(0, 1 + labl, lab[[2]], cex = labcex)
    segments(gap[2:10], 0, y1 = 1 - gap[2:10], lty = 3, col = "grey")
    segments(0, gap[2:10], x1 = 1 - gap[2:10], lty = 3, col = "grey")
    text(-margin, gap, gap, cex = tcex)
    text(gap, -margin, gap, cex = tcex)
    # Draw right triangle, labs and axis value
  } else {
    l <- matrix(c(0, 0, sidel, 0,
                  0, 0, sidel / 2, 1,
                  sidel / 2, 1, sidel, 0),
                nrow = 4)
    gap <- c(10, 25, 50, 75, 90) / 100
    temp <- gap / sqrt(3)
    text(sidel / 2, 1 + labl, lab[[1]], cex = labcex)
    text(-labl / sidel, -labl / 2, lab[[2]], cex = labcex, adj = 1)
    text(sidel + labl / sidel, -labl / 2, lab[[3]], cex = labcex, adj = 0)
    segments(temp, gap, sidel - temp, gap, lty = 3, col = "grey")
    segments(temp, gap, 2 * temp, 0, lty = 3, col = "grey")
    segments(sidel - temp, gap, sidel - 2 * temp, 0, lty = 3, col = "grey")
    text(temp - margin, gap, round(gap * faxis, 1), cex = tcex)
    text(sidel - temp + margin, gap, round((1 - gap) * faxis, 1), cex = tcex)
    text(2 * temp, -0.05, round((1 - gap) * faxis, 1), cex = tcex)
    # Draw equilateral triangle, labs, axis value and auxiliary line
  }
  segments(l[1, ], l[2, ], l[3, ], l[4, ], col = "grey50", lty = boundlty)
}

# Create a panel used to draw marginal of one variable
# dmax: range of y is from 0 to dmax in marginal plot
marginpanel <- function(lab = expression(x), dmax, cex = cex) {
  par(mar = rep(c(1.2, 4), 2), mgp = c(0, 0.5, 0))
  plot.new()
  plot.window(xlim = c(-0.01, 1), ylim = c(0, dmax), xaxs = "i", yaxs = "i")
  axis(2, cex.axis = 0.8, las = 1, tck = -0.02)
  par(mgp = rep(0, 3))
  axis(1, cex.axis = 0.8, tck = 0)
  mtext(lab, side = 1, at = 1.05, line = -0.4, cex = cex, adj = 0)
}

# Set color transparency
transcols <- function(cols, alpha) {
  paste0(cols, as.hexmode(round(alpha * 255)))
}

# Draw a sequence of points on simplex
# size: sum of each observation, for compositions is default sum to 1
# right: TRUE for right triangle FALSE for equilateral triangle
points_simplex <- function(data, lab = NULL, size = 1, right = FALSE,
                           add = FALSE, col = 2, pch = 1, cex = 0.5) {
  if (!add) {
    if (is.null(lab)) {
      lab <- list(expression(x[1]), expression(x[2]), expression(x[3]))
    }
    simplexpanel(lab, right = right)
  }
  if (is.vector(data)) data <- t(data)
  if (!right) {
    data <- data / size
    data <- cbind((2 - 2 * data[, 2] - data[, 1]) / sqrt(3), data[, 1])
  } else {
    data <- data[, 1:2, drop = FALSE] / size
  }
  points(data, col = col, pch = pch, cex = cex)
}

# Draw contour lines, modes and support points on the simplex
# of Dirichlet mixture's density
# x1, x2: if x1/x2 is a value, it should be an index of column (one variable);
#         if x1/x2 is a vector, it should be a set of index (multiple variables),
#         function will sum these variables as one variable (marginal).
# x1 is the left bottom's variable, x2 is the right bottom's variable,
# and top will be the rest of variables (marginal).
# prop: if FALSE, contour lines based on grid density value;
#       if TRUE, contour lines based on sample points.
plot.dirichletmix <- function(mix, data = NULL,
                              point = FALSE, mode = TRUE, contour = TRUE,
                              x1 = 2, x2 = 3, right = FALSE, axl = 100,
                              prop = FALSE, lv = NULL,
                              plev = log(seq(exp(0.5), exp(0.975), len = 5)),
                              drawlabels = TRUE, boundlty = 1, labcex = 0.9,
                              col = NULL, add = FALSE, title = NULL,
                              scale = 1) {
  if (!identical(intersect(x1, x2), numeric(0)))
    stop("Duplicate variables in x1 and x2")
  D <- ncol(mix$alpha)
  if (length(c(x1, x2)) >= D)
    stop("The number of variables in x1 and x2 is greater than the provided")
  if ((point || prop) && is.null(data)) stop("No data provided")
  xD <- seq_len(D)
  x3 <- setdiff(xD, c(x1, x2))
  lab1 <- paste(x1, collapse = ",")
  lab2 <- paste(x2, collapse = ",")
  lab <- c(bquote(x[.(lab1)]), bquote(x[.(lab2)])) # Labs of simplex
  if (!right) {
    # Equilateral triangle need three labs, right triangle need two labs
    lab3 <- paste(x3, collapse = ",")
    lab <-  c(bquote(x[.(lab3)]), lab)
  }
  if (!add) simplexpanel(lab = lab, right = right,
                         title = title, boundlty = boundlty)
  group <- rep(3, D)
  group[x1] <- 1
  group[x2] <- 2 # Group all compositions by x1, x2
  if (point || prop) { # Transform data points into marginal simplex
    sdata <- t(rowsum(t(data), group))
    sdata <- sdata[, c(3, 1, 2)]
    pdata <- if (right) sdata[, -1] else sdata
  }
  # Draw data points
  if (point) points_simplex(pdata, right = right, add = TRUE, col = "gray70")
  if (mode || contour) { # Draw mode, contour
    a <- mix$alpha
    sa <- t(rowsum(t(a), group))
    if (!right) sa <- sa[, c(3, 1, 2)]
    smix <- dirichletmix(sa, pr = mix$pr) # Marginal distribution
    if (contour) {
      ax1 <- seq(0, 1, l = axl)
      ax2 <- if (right) seq(0, 1, l = axl) else seq(0, 2 / sqrt(3), l = axl)
      x <- as.matrix(expand.grid(ax2, ax1))
      if (!right) x <- cbind(x[, 2], 1 - (sqrt(3) * x[, 1] + x[, 2]) / 2)
      # Transform to right triangle space
      x <- cbind(x, 1 - rowSums(x))
      i <- rowSums(x >= 0) == 3
      # Generate grid points

      # d <- numeric(length = axl^2)
      d <- rep(NA, length = axl^2)
      d[i] <- ddirichletmix(x[i, ], smix) # Grid density
      if (!prop && scale != 1) d <- d * scale
      # Scaling density, such as the prior for classification
      if (prop) { # Contour lines based on sample points
        dp <- ddirichletmix(sdata, smix)
        d[i] <- colMeans(outer(dp, d[i], "<"))
        if (is.null(lv)) lv <- seq(10, 90, by = 20) / 100
        lb <- paste0((1 - lv) * 100, "%")
      } else if (is.null(lv)) { # Contour lines based on grid density value
        lv <- quantile(d[i], probs = plev)
        lb <- round(lv, 2)
      }
      dim(d) <- c(axl, axl)
      if (is.null(col)) col <- viridis::viridis(length(lv))
      contour(ax2, ax1, d, levels = lv, labels = lb, labcex = labcex,
              col = col, add = TRUE, drawlabels = drawlabels)
    }
    if (mode) { # Draw mode
      points_simplex(attr(smix, "mode"), right = right, add = TRUE,
                     pch = 19, cex = 0.8)
    }
  }
  lv
}

# Mapping density from simplex to the Euclidean space
contour_s2r <- function(mix, data = NULL,
                        point = FALSE, mode = TRUE, contour = TRUE,
                        x1 = c(-3, 3), x2 = c(-3, 3), axl = 100,
                        prop = FALSE, lv = NULL, drawlabels = TRUE,
                        plev = log(seq(exp(0.5), exp(0.975), len = 5)),
                        pcex = 0.5, pcol = "gray70", title = NULL) {
  if ((point || prop) && is.null(data)) stop("No data provided")
  plot.new()
  plot.window(xlim = x1, ylim = x2)
  axis(1)
  axis(2)
  title(main = title)
  box()
  if (point) points(data, col = pcol, cex = pcex)
  if (contour) {
    x1 <- seq(x1[1], x1[2], l = axl)
    x2 <- seq(x2[1], x2[2], l = axl)
    x <- as.matrix(expand.grid(x1, x2))
    sx <- unclass(compositions::ilrInv(x))
    jacb <- exp(rowSums(log(sx)) + 1 / 2 * log(ncol(sx)))
    d <- ddirichletmix(sx, mix) * jacb
    if (prop) {
      sdata <- unclass(compositions::ilrInv(data))
      jacb <- exp(rowSums(log(sdata)) + 1 / 2 * log(ncol(sdata)))
      dp <- ddirichletmix(sdata, mix) * jacb
      d <- colMeans(outer(dp, d, "<"))
      if (is.null(lv)) lv <- seq(10, 90, by = 20) / 100
      lb <- paste0((1 - lv) * 100, "%")
    } else if (is.null(lv)) {
      lv <- quantile(d, probs = plev) # Proportion levels
      lb <- round(lv, 2)
    }
    dim(d) <- c(axl, axl)
    contour(x1, x2, d, levels = lv, labels = lb, labcex = 0.9, add = TRUE,
            col = viridis::viridis(length(lv)), main = title,
            drawlabels = drawlabels)
  }
  if (mode) {
    if (length(mix$pr) > 15 && length(unique(mix$pr)) == 1) return(NULL)
    md <- attr(mix, "mode")
    md <- md[rowSums(md == 0) == 0, , drop = FALSE]
    rmd <- compositions::ilr(md)
    # md <- unclass(compositions::ilr(mix$alpha))
    attributes(rmd) <- NULL
    dim(rmd) <- dim(md) - c(0, 1)
    points(rmd, col = 2, pch = 19, cex = 0.8)
  }
  lv
}

# Just like pairs function of Dirichlet mixture with every two variables
# For high dimensional data, pairs plot marginal density of each group
# group: a list or vector to receive the group number of each variables, no
# group by default
# gprest: TRUE for group rest of variables as a new column, FALSE for don't
# group rest variables
pairs.dirichletmix <- function(mix, data = NULL, group = NULL, gprest = FALSE,
                               right = FALSE, axl = 100,
                               lv = NULL, plev = NULL) {
  ax1 <- seq(0, 1, l = axl)
  ax2 <- if (right) seq(0, 1, l = axl) else seq(0, 2 / sqrt(3), l = axl)
  x <- as.matrix(expand.grid(ax2, ax1))
  if (!right) x <- cbind(x[, 2], 1 - (sqrt(3) * x[, 1] + x[, 2]) / 2)
  x <- cbind(x, 1 - rowSums(x))
  i <- rowSums(x >= 0) == 3 # Generate grid points
  len <- axl^2

  a <- mix$alpha
  a0 <- sum(a[1, ])
  m <- length(mix$pr)
  if (!is.null(group)) {
    if (is.list(group)) {
      group <- lapply(group, unique)
      if (any(duplicated(unlist(group)))) {
        stop("Duplicate columns in different groups")
      } # Ensure each column belongs to only one group
      D <- ncol(a)
      vec <- rep(NA, D)
      for (grp in names(group)) {
        vec[group[[grp]]] <- as.integer(gsub("[^0-9]", "", grp))
      } # list -> vector
      isna <- is.na(vec)
      remain <- setdiff(seq_len(D), vec[!isna])
      vec[isna] <- if (gprest) remain[1] else remain[1:sum(isna)]
      # Fill the remaining ungrouped columns
      group <- vec
    }
    unigroup <- sort(unique(group))
    if (length(unigroup) < 3) stop("Number of groups should larger than 3")
    if (any(diff(unigroup) != 1)) stop("Error with group number")
    a <- t(rowsum(t(a), group)) # Group components
    tag <- sapply(unigroup,
                  function(x) paste(which(group == x), collapse = ","))
    if (!is.null(data)) data <- t(rowsum(t(data), group)) # Group data
  } else {
    tag <- seq_along(a[1, ])
  }
  D <- ncol(a) # Update D after grouping
  D2 <- D * D
  grid <- matrix(1:(D2), nrow = D, byrow = TRUE)

  if (is.null(lv)) {
    nulv <- TRUE
    lv <- rep(NA, length = D2)
    if (is.null(plev)) plev <- log(seq(exp(0.5), exp(0.975), len = 5))
    lv <- matrix(nrow = D2, ncol = length(plev))
    colnames(lv) <- paste0(plev * 100, "%")
  } else {
    nulv <- FALSE
  }
  nlev <- ncol(lv)

  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(D, D), oma = c(0, 1, 0.5, 1))
  for (k in seq_len(D2)) {
    print(k)
    idx <- which(grid == k, arr.ind = TRUE)
    ir <- idx[1]
    ic <- idx[2]
    if (ir == ic) { # Marginal density for each group
      smix <- betamix(cbind(a[, ir], a0 - a[, ir]), pr = mix$pr)
      d <- dbetamix(ax1, smix)
      smd <- attr(smix, "mode")
      if (nulv) lv[k] <- max(d)
      if (!is.null(data)) {
        hmax <- max(hist(data[, ir], breaks = 10, plot = FALSE)$density)
      } else hmax <- 0
      marginpanel(bquote(x[.(tag[ir])]), max(lv[k], hmax) + 0.1, cex = 1.5)
      if (!is.null(data)) {
        hist(data[, ir], breaks = 10, freq = FALSE, add = TRUE)
      }
      lines(ax1, d, col = "#0000CD", lwd = 1.4)
      segments(smd, rep(0, m), y1 = smix$pr * max(d), col = 2, lwd = 1.4)
      next
    } else if (ir < ic) { # Marginal contour for each two groups
      lab <- c(bquote(x[.(tag[ir])]), bquote(x[.(tag[ic])]))
      sa <- cbind(a[, ir], a[, ic], a0 - rowSums(a[, idx, drop = FALSE]))
      if (!right) {
        lab <- c(bquote(x[.(paste(tag[-idx], collapse = ","))]), lab)
        sa <- sa[, c(3, 1, 2)]
      }
      smix <- dirichletmix(sa, pr = mix$pr)
      d <- numeric(length = len)
      d[i] <- ddirichletmix(x[i, ], smix)
      dim(d) <- c(axl, axl)
      if (nulv) lv[k, ] <- quantile(d[i], probs = plev)
      simplexpanel(lab = lab, right = right)
      if (!is.null(data)) {
        sdata <- cbind(data[, ir], data[, ic])
        if (!right) sdata <- cbind(1 - rowSums(data[, idx, drop = FALSE]),
                                   sdata)
        points_simplex(sdata, right = right, add = TRUE, col = "gray70")
      }
      contour(ax2, ax1, d, levels = lv[k, ], labels = round(lv[k, ], 2),
              col = viridis::viridis(nlev), add = TRUE)
      points_simplex(attr(smix, "mode"), right = right, add = TRUE,
                     pch = 16)
    } else {
      plot.new()
    }
  }
  par(old_par)
  lv
}

# Draw line segments between pairs of points on simplex
segments_simplex <- function(x, y, size = 1, add = FALSE,
                             lty = 3, lwd = 1, col = "gray40") {
  if (is.vector(x)) x <- t(x)
  if (is.vector(y)) y <- t(y)
  if (!add) simplexpanel()
  x <- cbind((x[, 3] * 2 + x[, 1]) / sqrt(3), x[, 1]) / size
  y <- cbind((y[, 3] * 2 + y[, 1]) / sqrt(3), y[, 1]) / size
  segments(x[, 1], x[, 2], y[, 1], y[, 2], col = col, lty = lty, lwd = lwd)
}

# Draw contour lines on simplex of Dirichlet mixture's gradient function
# w: log f(x_i; G); size: sum of alpha, a0
contour_gradient <- function(data, w, size, lv = NULL, axl = 100,
                             plev = log(seq(exp(0.5), exp(0.975), len = 5))) {
  x1 <- seq(0, 2 / sqrt(3), l = axl)
  x2 <- seq(0, 1, l = axl)
  x <- as.matrix(expand.grid(x1, x2))
  x <- cbind(x[, 2], 1 - (sqrt(3) * x[, 1] + x[, 2]) / 2)
  x <- cbind(x, 1 - rowSums(x))
  i <- rowSums(x > 0) == ncol(x)

  n <- nrow(data)
  d <- rep(NA, nrow(x))
  logd <- ddirichlet(data, x[i, ] * size, log = TRUE, outer = TRUE) - w
  d[i] <- colSums(exp(logd)) - n # Gradient value of alpha
  if (is.null(lv)) lv <- quantile(d[i], probs = plev) # Proportion levels
  dim(d) <- c(length(x1), length(x2))

  lab <- list(expression(alpha[1]), expression(alpha[2]), expression(alpha[3]))
  simplexpanel(lab = lab, faxis = size, margin = 0.06)
  sidel <- 2 / sqrt(3)
  l1 <- 1 / size
  l2 <- l1 * sqrt(3)
  l3 <- l1 * 2
  l <- matrix(c(l2, l1, sidel - l2, l1,
                l2, l1, sidel / 2, 1 - l3,
                sidel / 2, 1 - l3, sidel - l2, l1),
              nrow = 4)
  segments(l[1, ], l[2, ], l[3, ], l[4, ], col = "#0000CD", lty = 2)
  contour(x1, x2, d, levels = lv, labels = round(lv, 2), labcex = 0.8,
          col = viridis::viridis(length(lv)), add = TRUE)
}

# Draw the classification boundary of two estimation on simplex
contour_bondary <- function(est1, est2, prior, threshold = NULL,
                            ax1 = c(0, 2 / sqrt(3)), ax2 = c(0, 1),
                            axl = 100, level = 0,
                            col = 2, lty = 1, lwd = 1, add = FALSE) {
  ax1 <- seq(ax1[1], ax1[2], length = axl)
  ax2 <- seq(ax2[1], ax2[2], length = axl)
  x <- as.matrix(expand.grid(ax1, ax2))
  x <- cbind(x[, 2], 1 - (sqrt(3) * x[, 1] + x[, 2]) / 2)
  x <- cbind(x, 1 - rowSums(x))
  i <- rowSums(x >= 0) == 3
  d <- rep(NA, length = axl^2)
  d1 <- predmix(x[i, ], est1) # Grid density
  d2 <- predmix(x[i, ], est2)
  if (!is.null(threshold)) d1[d1 < threshold] <- 0
  d[i] <- d1 * prior[1] - d2 * prior[2]
  dim(d) <- c(axl, axl)
  if (!add) simplexpanel()
  contour(ax1, ax2, d, levels = level, drawlabels = FALSE, add = TRUE,
          col = col, lty = lty, lwd = lwd)
}
