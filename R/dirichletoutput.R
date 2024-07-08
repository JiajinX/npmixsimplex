## Plot and print function of Dirichlet mixture ##

# Create a simplex panel
simplexpanel <- function() {
  sidel <- 2 / sqrt(3)
  liml <- 0.1
  labl <- 0.06
  plot.new()
  plot.window(xlim = c(-liml, sidel + liml), ylim = c(-liml, 1 + liml),
              xaxs = "i", yaxs = "i", asp = 1)

  text(sidel / 2, 1 + labl, expression(x[1]))
  text(-labl / sidel, -labl / 2, expression(x[2]))
  text(sidel + labl / sidel, -labl / 2, expression(x[3]))
  l <- matrix(c(0, 0, sidel, 0,
                0, 0, sidel / 2, 1,
                sidel / 2, 1, sidel, 0),
              nrow = 4)
  segments(l[1, ], l[2, ], l[3, ], l[4, ])

  margin <- 0.04
  gap <- c(5, 10, 25, 50, 75, 90, 95) / 100
  temp <- gap / sqrt(3)
  segments(temp, gap, sidel - temp, gap, lty = 3, col = "grey")
  segments(temp, gap, 2 * temp, 0, lty = 3, col = "grey")
  segments(sidel - temp, gap, sidel - 2 * temp, 0, lty = 3, col = "grey")
  text(temp - margin, gap, gap * 100, cex = 0.6)
  text(sidel - temp + margin, gap, (1 - gap) * 100, cex = 0.6)
  text(2 * temp, -margin, (1 - gap) * 100, cex = 0.6)
}

# Draw a sequence of points on simplex
points.simplex <- function(data, size = 1, add = FALSE,
                           col = 2, pch = 1, cex = 0.5) {
  if (!add) simplexpanel()
  data <- cbind((data[, 3] * 2 + data[, 1]) / sqrt(3), data[, 1]) / size
  points(data, col = col, pch = pch, cex = cex)
}

# Draw support points/mode of Dirichlet mixture on simplex
points.dirichletmix <- function(mix, add = FALSE,
                                col = 4, pch = 16, cex = 0.6) {
  data <- attr(mix, "mode")
  points.simplex(data, add = add, col = col, pch = pch, cex = cex)
}

# Draw line segments between pairs of points on simplex
segments.simplex <- function(x, y, size = 1, add = FALSE,
                             lty = 3, col = "gray40") {
  if (!add) simplexpanel()
  x <- cbind((x[, 3] * 2 + x[, 1]) / sqrt(3), x[, 1]) / size
  y <- cbind((y[, 3] * 2 + y[, 1]) / sqrt(3), y[, 1]) / size
  segments(x[, 1], x[, 2], y[, 1], y[, 2], col = col, lty = lty)
}

# Draw contour lines on simplex of Dirichlet mixture's density
contour.dirichletmix <- function(mix, data = NULL, lv = NULL, nlev = 15) {
  x1 <- seq(0, 2 / sqrt(3), l = 100)
  x2 <- seq(0, 1, l = 100)
  x <- as.matrix(expand.grid(x1, x2))
  x <- cbind(x[, 2], 1 - (sqrt(3) * x[, 1] + x[, 2]) / 2)
  x <- cbind(x, 1 - rowSums(x))
  i <- rowSums(x > 0) == ncol(x)

  d <- numeric(length = nrow(x))
  d[i] <- ddirichletmix(x[i, ], mix)
  dm <- matrix(d, nrow = length(x1), ncol = length(x2))
  if (is.null(lv)) lv <- seq(0, max(d), l = nlev + 1)

  simplexpanel()
  contour(x1, x2, dm, levels = lv, xlab = "x", ylab = "y",
          col = viridis::viridis(nlev), add = TRUE)

  if (!is.null(data)) {
    points.simplex(data, add = TRUE)
  }

  attr(d, "level") <- lv
  d
}

# Draw contour lines on simplex of Dirichlet mixture's gradient function
contour.gradient <- function(data, w, size, lv = NULL, nlev = 15) {
  x1 <- seq(0, 2 / sqrt(3), l = 100)
  x2 <- seq(0, 1, l = 100)
  x <- as.matrix(expand.grid(x1, x2))
  x <- cbind(x[, 2], 1 - (sqrt(3) * x[, 1] + x[, 2]) / 2)
  x <- cbind(x, 1 - rowSums(x))
  i <- rowSums(x > 0) == ncol(x)

  n <- nrow(data)
  d <- numeric(length = nrow(x))
  logd <- outer.dirichlet(data, x[i, ] * size, log = TRUE) - w
  d[i] <- colSums(exp(logd)) - n # Gradient value of alpha
  dm <- matrix(d, nrow = length(x1), ncol = length(x2))
  if (is.null(lv)) lv <- seq(min(d), max(d), l = nlev + 1)

  simplexpanel()
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
          col = viridis::viridis(nlev), add = TRUE)
}

# Print information of Dirichlet mixture of each iterations
print_verbose <- function(verbose, i, ll, mix) {
  if (verbose >= 1) cat("\n## Iteration ", i, ":\nLog-likelihood = ",
                        as.character(signif(ll[1], 16)), "\n", sep = "")
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
