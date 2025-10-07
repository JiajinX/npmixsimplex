summaryre <- function(x, type = c("mean", "se")) {
  type <- match.arg(type)
  df <- as.data.frame(x)
  lst <- split(df, df$no)
  do.call(rbind, lapply(lst, function(x){
    x <- x[, -1]
    if (type == "mean") {
      colMeans(x, na.rm = TRUE)
    } else apply(x, 2, sd, na.rm = TRUE) / sqrt(nrow(x))
  }))
}

letaxre <- function(re, prod = 1, rd = 2) {
  rd <- paste0("%.", rd, "f")
  mu <- sprintf(rd, colMeans(re) * prod)
  se <- sprintf(rd, apply(re, 2, sd) / sqrt(nrow(re)) * prod)
  txt <- paste0("$", mu, "_{", se, "}$")
  names(txt) <- colnames(re)
  txt
}

testre <- function(re, type = c("greater", "less")) {
  D <- ncol(re)
  pvalue <- rep(NA, D)
  names(pvalue) <- colnames(re)
  mu <- colMeans(re)
  idx <- ifelse(type == "greater", which.max(mu), which.min(mu))
  for (i in seq_len(D)) {
    tt <- t.test(re[, idx], re[, i], alternative = "greater", paired = TRUE)
    pvalue[i] <- tt$p.value
  }
  equal <- pvalue > 0.05
  attr(equal, "p-value") <- pvalue
  equal
}
