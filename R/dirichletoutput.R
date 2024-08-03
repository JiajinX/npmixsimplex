## ============== Plot and print function of Dirichlet mixture ============== ##

# Create a simplex panel
simplexpanel <- function(lab = list(expression(x[1]), expression(x[2]),
                                    expression(x[3])), right = FALSE) {
  par(mar = rep(0, 4))
  liml <- 0.1
  labl <- 0.06
  margin <- 0.04
  sidel <- ifelse(right, 1, 2 / sqrt(3))
  plot.new()
  plot.window(xlim = c(-liml, sidel + liml), ylim = c(-liml, 1 + liml),
              xaxs = "i", yaxs = "i", asp = 1)
  if (right) {
    l <- matrix(c(0, 0, 0, 1,
                  0, 0, 1, 0,
                  0, 1, 1, 0),
                nrow = 4)
    gap <- seq(0, 1, by = 0.1)
    text(1 + labl, 0, lab[[1]])
    text(0, 1 + labl, lab[[2]])
    segments(gap[2:10], 0, y1 = 1 - gap[2:10], lty = 3, col = "grey")
    segments(0, gap[2:10], x1 = 1 - gap[2:10], lty = 3, col = "grey")
    text(-margin, gap, gap, cex = 0.6)
    text(gap, -margin, gap, cex = 0.6)
    # Draw right triangle, labs and axis value
  } else {
    l <- matrix(c(0, 0, sidel, 0,
                  0, 0, sidel / 2, 1,
                  sidel / 2, 1, sidel, 0),
                nrow = 4)
    gap <- c(5, 10, 25, 50, 75, 90, 95) / 100
    temp <- gap / sqrt(3)
    text(sidel / 2, 1 + labl, lab[[1]])
    text(-labl / sidel, -labl / 2, lab[[2]])
    text(sidel + labl / sidel, -labl / 2, lab[[3]])
    segments(temp, gap, sidel - temp, gap, lty = 3, col = "grey")
    segments(temp, gap, 2 * temp, 0, lty = 3, col = "grey")
    segments(sidel - temp, gap, sidel - 2 * temp, 0, lty = 3, col = "grey")
    text(temp - margin, gap, gap * 100, cex = 0.6)
    text(sidel - temp + margin, gap, (1 - gap) * 100, cex = 0.6)
    text(2 * temp, -margin, (1 - gap) * 100, cex = 0.6)
    # Draw equilateral triangle, labs, axis value and auxiliary line
  }
  segments(l[1, ], l[2, ], l[3, ], l[4, ])
}

# Create a panel used to draw marginal of one variable
# dmax: range of y is from 0 to dmax in marginal plot
marginpanel <- function(lab = expression(x), dmax) {
  par(mar = rep(c(1.2, 4), 2), mgp = c(0, 0.5, 0))
  plot.new()
  plot.window(xlim = c(-0.01, 1), ylim = c(0, dmax), xaxs = "i", yaxs = "i")
  axis(2, cex.axis = 0.6, las = 1, tck = -0.02)
  par(mgp = rep(0, 3))
  axis(1, cex.axis = 0.6, tck = 0)
  mtext(lab, side = 1, at = 1.05, line = -0.4, cex = 0.7, adj = 0)
}

# Draw a sequence of points on simplex
# size: sum of each observation, for compositional data is defualt sum to 1
# right: TRUE for right triangle FALSE for equilateral triangle
points.simplex <- function(data, lab = NULL, size = 1, right = FALSE,
                           add = FALSE, col = 2, pch = 1, cex = 0.5) {
  if (!add) {
    if (is.null(lab)) {
      lab = list(expression(x[1]), expression(x[2]), expression(x[3]))
    }
    simplexpanel(lab, right = right)
  }
  if (is.vector(data)) data <- t(data)
  if (!right) {
    data <- cbind((2 - 2 * data[, 2] - data[, 1]) / sqrt(3), data[, 1]) / size
  } else data <- data[, 1:2, drop = FALSE] / size
  points(data, col = col, pch = pch, cex = cex)
}

# Draw contour lines and support points/mode
# on simplex of Dirichlet mixture's density
# x1, x2: can be value or vector, if x1/x2 is a value, it should be an index
# of column (one variable). x1 is the left bottom's variable, and x2 is the
# right bottom's variable. Top will be the rest of variables (marginal). If
# x1/x2 is a vector, it should be a set of index (multiple variables),
# function will sum these variables as one variable (marginal).
plot.dirichletmix <- function(mix, data = NULL,
                              type = c("both", "mode", "contour", "null"),
                              x1 = 2, x2 = 3, right = FALSE,
                              lv = NULL, plev = NULL, axl = 100) {
  if (!identical(intersect(x1,x2), numeric(0)))
    stop("Duplicate variables in x1 and x2")
  D <- ncol(mix$alpha)
  if (length(c(x1, x2)) >= D)
    stop("The number of variables in x1 and x2 is greater than the provided")
  type = match.arg(type)
  xD <- seq_len(D)
  x3 <- setdiff(xD, c(x1, x2))
  lab1 <- paste(x1, collapse = ",")
  lab2 <- paste(x2, collapse = ",")
  lab <- c(bquote(x[.(lab1)]), bquote(x[.(lab2)]))
  if (!right) {
    # Equilateral triangle need three labs, right triangle need two labs
    lab3 <- paste(x3, collapse = ",")
    lab <-  c(bquote(x[.(lab3)]), lab)
  }
  simplexpanel(lab = lab, right = right)
  group <- rep(3, D)
  group[x1] <- 1
  group[x2] <- 2 # Group all compositions by x1, x2
  if (!is.null(data)) { # Draw data
    sdata <- t(rowsum(t(data), group))
    sdata <- sdata[, c(3, 1, 2)]
    if (right) sdata <- sdata[, -1]
    points.simplex(sdata, right = right, add = TRUE, col = "gray70")
  }
  if (type != "null") { # Draw contour, mode
    ax1 <- seq(0, 1, l = axl)
    ax2 <- if (right) seq(0, 1, l = axl) else seq(0, 2 / sqrt(3), l = axl)
    x <- as.matrix(expand.grid(ax2, ax1))
    if (!right) x <- cbind(x[, 2], 1 - (sqrt(3) * x[, 1] + x[, 2]) / 2)
    # Transform to right triangle space
    x <- cbind(x, 1 - rowSums(x))
    i <- rowSums(x >= 0) == 3

    a <- mix$alpha
    a0 <- sum(a[1, ])
    sa <- t(rowsum(t(a), group))
    if (!right) sa <- sa[, c(3, 1, 2)]
    smix <- dirichletmix(sa, pr = mix$pr) # Marginal distribution
    d <- numeric(length = axl^2)
    d[i] <- ddirichletmix(x[i, ], smix)
    dm <- matrix(d, nrow = axl, ncol = axl)
    if (is.null(lv)) {
      if (is.null(plev)) {
        plev <- c(50, 65, 75, 80, 85, 90, 94, 97, 99, 100) / 100
      }
      lv <- quantile(d[i], probs = plev) # Proportion levels
    }
    if (type != "mode") {
      contour(ax2, ax1, dm, levels = lv, col = viridis::viridis(length(lv)),
              add = TRUE)
    }
    if (type != "contour") {
      points.simplex(attr(smix, "mode"), right = right,
                     add = TRUE, pch = 16, cex = 0.6)
    }
    lv
  }
}

# Just like pairs function of Dirichlet mixture with every two variables
# For high dimensional data, pairs plot marginal density of each group
pairs.dirichletmix <- function(mix, data = NULL, group = NULL, right = FALSE,
                               lv = NULL, plev = NULL, axl = 100) {
  ax1 <- seq(0, 1, l = axl)
  ax2 <- if (right) seq(0, 1, l = axl) else seq(0, 2 / sqrt(3), l = axl)
  x <- as.matrix(expand.grid(ax2, ax1))
  if (!right) x <- cbind(x[, 2], 1 - (sqrt(3) * x[, 1] + x[, 2]) / 2)
  x <- cbind(x, 1 - rowSums(x))
  i <- rowSums(x >= 0) == 3
  len <- axl^2

  a <- mix$alpha
  a0 <- sum(a[1, ])
  m <- length(mix$pr)
  if (!is.null(group)) {
    a <- t(rowsum(t(a), group))
    tag <- sapply(sort(unique(group)), function(x)
      paste(which(group == x), collapse = ","))
    if (!is.null(data)) data <- t(rowsum(t(data), group))
  } else tag <- seq_along(a[1,])
  D <- ncol(a)
  D2 <- D * D
  grid <- matrix(1:(D2), nrow = D, byrow = TRUE)

  if (is.null(lv)) {
    nulv <- TRUE
    lv <- rep(NA, length = D2)
    if (is.null(plev)) plev <- c(50, 65, 75, 80, 85, 90, 94, 97, 99, 100) / 100
    lv <- matrix(nrow = D2, ncol = length(plev))
    colnames(lv) <- paste0(plev * 100, "%")
  } else nulv <- FALSE
  nlev <- ncol(lv)

  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(D, D), oma = c(0, 1, 0.5, 1))
  for (k in seq_len(D2)) {
    idx <- which(grid == k, arr.ind = TRUE)
    ir <- idx[1]
    ic <- idx[2]
    if (ir == ic) { # Marginal density for each group
      smix <- betamix(cbind(a[, ir], a0 - a[, ir]), pr = mix$pr)
      d <- dbetamix(ax1, smix)
      smd <- attr(smix, "mode")
      if (nulv) lv[k] <- max(d)
      marginpanel(bquote(x[.(tag[ir])]), lv[k] + 0.1)
      if (!is.null(data)) {
        hist(data[, ir], breaks = 20, freq = FALSE, add = TRUE)
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
      dm <- matrix(d, nrow = axl, ncol = axl)
      if (nulv) lv[k, ] <- quantile(d[i], probs = plev)
      simplexpanel(lab = lab, right = right)
      if (!is.null(data)) {
        sdata <- cbind(data[, ir], data[, ic])
        if (!right) sdata <- cbind(1 - rowSums(data[, idx, drop = FALSE]),
                                   sdata)
        points.simplex(sdata, right = right, add = TRUE, col = "gray70")
      }
      contour(ax2, ax1, dm,
              levels = lv[k, ], col = viridis::viridis(nlev), add = TRUE)
      points.simplex(attr(smix, "mode"), right = right, add = TRUE,
                     pch = 16, cex = 0.6)
    } else plot.new()
  }
  par(old_par)
  lv
}

# Draw line segments between pairs of points on simplex
segments.simplex <- function(x, y, size = 1, add = FALSE,
                             lty = 3, col = "gray40") {
  if (!add) simplexpanel()
  x <- cbind((x[, 3] * 2 + x[, 1]) / sqrt(3), x[, 1]) / size
  y <- cbind((y[, 3] * 2 + y[, 1]) / sqrt(3), y[, 1]) / size
  segments(x[, 1], x[, 2], y[, 1], y[, 2], col = col, lty = lty)
}

# Draw contour lines on simplex of Dirichlet mixture's gradient function
# w: log f(x_i; G); size: sum of alpha, a0
contour_gradient <- function(data, w, size, lv = NULL, plev = NULL, axl = 100) {
  x1 <- seq(0, 2 / sqrt(3), l = axl)
  x2 <- seq(0, 1, l = axl)
  x <- as.matrix(expand.grid(x1, x2))
  x <- cbind(x[, 2], 1 - (sqrt(3) * x[, 1] + x[, 2]) / 2)
  x <- cbind(x, 1 - rowSums(x))
  i <- rowSums(x > 0) == ncol(x)

  n <- nrow(data)
  d <- numeric(length = nrow(x))
  logd <- ddirichlet(data, x[i, ] * size, log = TRUE, outer = TRUE) - w
  d[i] <- colSums(exp(logd)) - n # Gradient value of alpha
  dm <- matrix(d, nrow = length(x1), ncol = length(x2))
  if (is.null(lv)) {
    if (is.null(plev)) {
      plev <- c(50, 65, 75, 80, 85, 90, 94, 97, 99, 100) / 100
    }
    lv <- quantile(d[i], probs = plev) # Proportion levels
  }

  lab <- list(expression(alpha[1]), expression(alpha[2]), expression(alpha[3]))
  simplexpanel(lab = lab)
  sidel <- 2 / sqrt(3)
  l1 <- 1 / size
  l2 <- l1 * sqrt(3)
  l3 <- l1 * 2
  l <- matrix(c(l2, l1, sidel - l2, l1,
                l2, l1, sidel / 2, 1 - l3,
                sidel / 2, 1 - l3, sidel - l2, l1),
              nrow = 4)
  segments(l[1, ], l[2, ], l[3, ], l[4, ], col = "#0000CD", lty = 2)
  contour(x1, x2, dm, levels = lv, xlab = "x", ylab = "y",
          col = viridis::viridis(length(lv)), add = TRUE)
}

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
