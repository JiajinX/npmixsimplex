# Density of KDE, MDE and mvnormix on simplex
contour_r2s <- function(model, data = NULL, jacobian = TRUE,
                        point = FALSE, mode = TRUE, contour = TRUE,
                        axl = 100, prop = FALSE, lv = NULL, drawlabels = TRUE,
                        plev = log(seq(exp(0.5), exp(0.975), len = 5)),
                        title = NULL, pcex = 0.5, pcol = "gray70",
                        add = FALSE, boundlty = 1) {
  if ((point || prop) && is.null(data)) stop("No data provided")
  if (!add) simplexpanel(title = title, boundlty = boundlty)
  if (point) points_simplex(data, col = pcol, cex = pcex, add = TRUE)
  if (contour || mode) {
    cls <- class(model)
    type <- cls[1] # MDE, KDE, mvnormix
    tf <- cls[2] # ilr or alr
    if (is.na(tf)) tf <- "ilr"

    if (contour) {
      x1 <- seq(0, 2 / sqrt(3), l = axl)
      x2 <- seq(0, 1, l = axl)
      x <- as.matrix(expand.grid(x1, x2))
      x <- cbind(x[, 2], 1 - (sqrt(3) * x[, 1] + x[, 2]) / 2)
      x <- cbind(x, 1 - rowSums(x))
      i <- rowSums(x > 0) == 3 # Grid in simplex
      jcb <- exp(rowSums(log(x[i, ])))
      if (tf == "ilr") {
        y <- compositions::ilr(x[i, ])
        jcb <- jcb * sqrt(3)
      } else y <- compositions::alr(x[i, ])
      y <- unclass(y)

      d <- numeric(length = nrow(x))
      d[i] <- if (type == "kde") predict(model, x = y) else dmvnormix(y, model)
      if (jacobian) d[i] <- d[i] / jcb
      if (prop) {
        jcb <- exp(rowSums(log(data)))
        if (tf == "ilr") {
          rdata <- compositions::ilr(data)
          jcb <- jcb * sqrt(3)
        } else rdata <- compositions::alr(data)
        rdata <- unclass(rdata)

        dp <- if (cls == "kde") {
          predict(model, x = rdata)
        } else dmvnormix(rdata, model)
        if (jacobian) dp <- dp / jcb
        d[i] <- colMeans(outer(dp, d[i], "<"))
        if (is.null(lv)) lv <- seq(10, 90, by = 20) / 100
        lb <- paste0((1 - lv) * 100, "%")
      } else if (is.null(lv)) {
        lv <- quantile(d[i], probs = plev) # Proportion levels
        lb <- round(lv, 2)
      }
      dim(d) <- c(axl, axl)
      contour(x1, x2, d, levels = lv, labels = lb, labcex = 0.9,
              col = viridis::viridis(length(lv)), add = TRUE,
              drawlabels = drawlabels)
    }

    if (type != "kde" && mode) {
      md <- model$mu
      md <- if (tf == "ilr") {
        compositions::ilrInv(md)
      } else compositions::alrInv(md)
      md <- unclass(md)
      points_simplex(md, add = TRUE, pch = 19, cex = 0.8)
    }
  }
  lv
}

# Plot Gaussian mixture distribution density on real space
contour_r <- function(model, data = NULL,
                      point = FALSE, mode = TRUE, contour = TRUE,
                      x1 = c(-3, 3), x2 = c(-3, 3),
                      prop = FALSE, lv = NULL, drawlabels = TRUE,
                      plev = log(seq(exp(0.5), exp(0.975), len = 5)),
                      axl = 100, title = NULL, pcex = 0.5, pcol = "gray70") {
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
    cls <- class(model)[1]
    d <- if (cls == "kde") predict(model, x = x) else dmvnormix(x, model)
    if (prop) {
      dp <- if (cls == "kde") {
        predict(model, x = data)
      } else dmvnormix(data, model)
      d <- colMeans(outer(dp, d, "<"))
      if (is.null(lv)) lv <- seq(10, 90, by = 20) / 100
      lb <- paste0((1 - lv) * 100, "%")
    } else if (is.null(lv)) {
      lv <- quantile(d, probs = plev) # Proportion levels
      lb <- round(lv, 2)
    }
    dim(d) <- c(axl, axl)
    contour(x1, x2, d, levels = lv, labels = lb, add = TRUE,
            drawlabels = drawlabels, col = viridis::viridis(length(lv)),
            main = title)
  }
  if (cls != "kde" && mode) {
    md <- model$mu
    points(md, col = 2, pch = 16, cex = 0.7)
  }
  lv
}
