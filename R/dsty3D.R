cubic <- function(x = c(-4, 4), y = c(-4, 4), z = c(-4, 4)) {
  p <- matrix(c(x[1], y[1], z[1],
                x[1], y[1], z[2],
                x[1], y[2], z[1],
                x[1], y[2], z[2],
                x[2], y[1], z[1],
                x[2], y[1], z[2],
                x[2], y[2], z[1],
                x[2], y[2], z[2]
              ), ncol = 3, byrow = TRUE)
  rgl::segments3d(c(p[1, ], p[2, ]))
  rgl::segments3d(c(p[1, ], p[3, ]))
  rgl::segments3d(c(p[1, ], p[5, ]))
  rgl::segments3d(c(p[2, ], p[4, ]))
  rgl::segments3d(c(p[2, ], p[6, ]))
  rgl::segments3d(c(p[3, ], p[4, ]))
  rgl::segments3d(c(p[3, ], p[7, ]))
  rgl::segments3d(c(p[4, ], p[8, ]))
  rgl::segments3d(c(p[5, ], p[7, ]))
  rgl::segments3d(c(p[5, ], p[6, ]))
  rgl::segments3d(c(p[6, ], p[8, ]))
  rgl::segments3d(c(p[7, ], p[8, ]))

  gap <- 0.5
  rgl::text3d(sum(x) / 2, y[1] - gap, z[1], text = "X")
  rgl::text3d(x[1] - gap, sum(y) / 2, z[1], text = "Y")
  rgl::text3d(x[1], y[1] - gap, sum(z) / 2, text = "Z")
}

# n: number of grids
contour_chacon3D <- function(no = 1, n = 50, range = 3,
                             plev = log(seq(exp(0.7), exp(0.975), len = 5)),
                             alo = 0.1, ahi = 0.5) {
  k <- length(plev)
  col <- viridis::viridis(k)
  al <- seq(alo, ahi, len = k)
  rg <- c(-range, range)
  x <- seq(rg[1], rg[2], len = n)
  grid <- expand.grid(x1 = x, x2 = x, x3 = x)
  mix <- chaconmix(no, D = 3)
  d <- dmvnormix(grid, mix)
  th <- quantile(d, probs = plev)
  chacon3D <- function(x, y, z, no) dmvnormix(cbind(x, y, z), mix)
  f <- function(x, y, z) chacon3D(x, y, z, no)
  misc3d::contour3d(f, th, x, x, x, color = col, alpha = al)
  cubic(x = rg, y = rg, z = rg)
}

# Barycentric coordinate -> Cartesian coordinates
# data: simplex data with 4-part components
b2c <- function(data) {
  if (is.vector(data)) data <- t(data)
  p <- matrix(c(sqrt(6) / 4, 3 / 4 * sqrt(2), 0,
                0, 0, 0,
                sqrt(6) / 2, 0, 0,
                sqrt(6) / 4, 1 / 4 * sqrt(2), 1),
              ncol = 3, byrow = TRUE)
  data %*% p
}

# Cartesian coordinates -> Barycentric coordinate
c2b <- function(data) {
  if (is.vector(data)) data <- t(data)
  p <- matrix(c(sqrt(6) / 4, 3 / 4 * sqrt(2), 0,
                0, 0, 0,
                sqrt(6) / 2, 0, 0,
                sqrt(6) / 4, 1 / 4 * sqrt(2), 1),
              ncol = 3, byrow = TRUE)
  inv_p <- solve(t(cbind(1, p)))
  tcrossprod(cbind(1, data), inv_p)
}

# Plot regular tetrahedron
tetrahedron <- function(add = FALSE) {
  p <- matrix(c(sqrt(6) / 4, 3 / 4 * sqrt(2), 0,
                0, 0, 0,
                sqrt(6) / 2, 0, 0,
                sqrt(6) / 4, 1 / 4 * sqrt(2), 1),
              ncol = 3, byrow = TRUE)
  if (add) open3d()
  rgl::segments3d(c(p[1, ], p[2, ]))
  rgl::segments3d(c(p[1, ], p[3, ]))
  rgl::segments3d(c(p[1, ], p[4, ]))
  rgl::segments3d(c(p[2, ], p[3, ]))
  rgl::segments3d(c(p[2, ], p[4, ]))
  rgl::segments3d(c(p[3, ], p[4, ]))

  # Add lab
  gap <- 0.05
  cex <- 1.5
  rgl::text3d(p[1, 1], p[1, 2] + gap, p[1, 3], text = expression(x[1]), cex = cex)
  rgl::text3d(p[2, 1] - gap, p[2, 2], p[2, 3], text = expression(x[2]), cex = cex)
  rgl::text3d(p[3, 1] + gap, p[3, 2], p[3, 3], text = expression(x[3]), cex = cex)
  rgl::text3d(p[4, 1], p[4, 2], p[4, 3] + gap, text = expression(x[4]), cex = cex)
}

# Points in regular tetrahedron
# data: simplex data with 4-part components
points_tetrahedron <- function(data, add = FALSE, col = 2) {
  coord <- b2c(data)
  if (!add) tetrahedron()
  rgl::points3d(coord[, 1], coord[, 2], coord[, 3], col = col)
}

contour3D <- function(mix, data = NULL, n = 50, prop = FALSE, lv = NULL,
                      col = NULL, alo = 0.1, ahi = 0.5, add = FALSE) {
  if (prop && is.null(data)) stop("No data provided")
  if (is.null(lv)) {
    lv <- if (prop) {
      seq(0.1, 0.9, len = 3)
    } else log(seq(exp(0.7), exp(0.975), len = 3))
  }
  k <- length(lv)
  if (is.null(col)) col <- viridis::viridis(k)
  al <- seq(alo, ahi, len = k)
  x <- seq(0, sqrt(6) / 2, len = n)
  y <- seq(0, 3 / 4 * sqrt(2), len = n)
  z <- seq(0, 1, len = n)
  if (!prop) {
    grid <- as.matrix(expand.grid(x1 = x, x2 = y, x3 = z))
    # Points by Cartesian coordinates
    coord <- c2b(grid) # Cartesian coordinates -> Barycentric coordinate
    coord <- coord / rowSums(coord)
    coord <- coord[rowSums(coord > 0) == 4, ]
    d <- predmix(coord, mix)
    lv <- quantile(d, probs = lv)
  } else dp <- predmix(data, mix)

  plot3D <- function(x, y, z) {
    coord <- c2b(cbind(x, y, z))
    coord <- coord / rowSums(coord)
    d <- numeric(length = nrow(coord))
    idx <- which(rowSums(coord > 0) == 4)
    d[idx] <- predmix(coord[idx, ], mix)
    if (prop) d[idx] <- colMeans(outer(dp, d[idx], "<"))
    d
  }
  f <- function(x, y, z) plot3D(x, y, z)
  misc3d::contour3d(f, lv, x, y, z, color = col, alpha = al, add = add)
  if (!add) tetrahedron()
  # view3d(theta = 0, phi = -3, zoom = 0.2)
  rgl::view3d(theta = 7, phi = 25, zoom = 0.3)
}


# contour3D(Dir_mix(10, D = 4), alo = 0.2, ahi = 0.2,
#           col = paletteer_c("ggthemes::Blue-Teal", 3),
#           lv = log(seq(exp(0.7), exp(0.975), len = 3)), n = 100)

# contour3D(chaconmix(6, D = 3), n = 80, alo = 0.1, ahi = 0.3,
#           lv = log(seq(exp(0.7), exp(0.975), len = 3)))

# par3d(windowRect = c(0, 0, 1500, 1500))
# for (i in 1:12) {
#   mix <- Dir_mix(i, D = 4)
#   # mix <- chaconmix(i, D = 3)
#   contour3D(mix, n = 100,
#             alo = 0.2, ahi = 0.2,
#             col = paletteer_c("ggthemes::Blue-Teal", 3),
#             lv = log(seq(exp(0.7), exp(0.975), len = 3)))
#   rgl.snapshot(paste0("dirichlet", i, ".png"), fmt = "png")
# }

# rgl.postscript("chacon1.pdf", fmt = "pdf", drawText = TRUE)

