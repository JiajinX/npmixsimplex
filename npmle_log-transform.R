library(gtools)
library(scales)
library(viridis)
library(npden)
library(ks)
library(compositions)

boundarymix <- function(testno, h, pr) {
  mix = switch(testno,
               dirichletmix(u = rbind(c(0,0.5,0.5)), h = h, pr = 1), # center of boundary
               dirichletmix(u = rbind(c(0,1/4,3/4)), h = h, pr = 1), # boundary
               dirichletmix(u = rbind(c(1,0,0)), h = h, pr = 1), # corner
               dirichletmix(u = rbind(c(1,0,0), c(0,0,1)), h = h, pr = pr), # two corner
               dirichletmix(u = rbind(c(1,0,0), c(0,0.5,0.5)), h = h, pr = pr), # corner - opposite center boundary
               dirichletmix(u = rbind(c(1,0,0), c(0.5,0,0.5)), h = h, pr = pr), # corner - near center boundary
               dirichletmix(u = rbind(c(1,0,0), c(3/4,0,1/4)), h = h, pr = pr), # corner - near boundary
               dirichletmix(u = rbind(c(0,3/4,1/4), c(0.5,0,0.5)), h = h, pr = pr) # two boundary
  )
  mix
}

# Dirichlet mixture parameters
dirichletmix <- function(u, h, pr) {
  m = length(pr)
  u = matrix(u, nr = m)
  a = u / h + 1
  attr(a, 'mode') = u
  attr(a, 'h') = h
  mix = list('alpha' = a, 'pr' = pr)
  attr(mix, 'components') = m
  class(mix) = 'dirichletmix'
  mix
}

# generate random data from Dirichlet mixture distribution
# threshold2zero: if value is less than threshold, let it be 0 
rdirichletmix <- function(n, mix, threshold2zero = NULL) {
  m = attr(mix, 'components')
  i = sample.int(m, n, prob = mix$pr, replace = TRUE)
  data = list()
  for (j in 1:m) { # number of components
    data[[j]] = rdirichlet(sum(i==j), alpha = mix$alpha[j,])
  }
  data = do.call(rbind, data)
  
  if (!is.null(threshold2zero)) { # create points on boundary
    r = unique(which(data < threshold2zero, arr.ind = T)[,1])
    mx = data[r,]
    mx1 = mx2 = mx
    mx1[mx1 >= threshold2zero] = NA
    mx2[mx2 < threshold2zero] = NA
    add = rowSums(mx1, na.rm = T)/rowSums(is.na(mx1))
    mx3 = mx2 + add
    mx3[is.na(mx3)] = 0
    data[r,] = mx3
  }
  data
}

# density of Dirichlet mixture distribution
ddirichletmix <- function(x, mix) {
  D = ncol(x)
  if (ncol(x) + 1 == ncol(mix$alpha)) {
    D = D + 1
    x = cbind(x, 1-rowSums(x))
  }
  i = rowSums(x > 0) == D
  d = numeric(length = nrow(x))
  d[i] = apply(mix$alpha, 1, dDirichlet, x = x[i,]) %*% mix$pr
  d
}

# mixture parameters of chacon
chaconmix <- function(testno) {
  mix = switch(testno,
               list(mus = rbind(c(0,0)), 
                    Sigmas = invvech(c(1/4,0,1)),
                    pr = 1),
               list(mus = rbind(c(0,0)), 
                    Sigmas = invvech(c(1,9/10,1)),
                    pr = 1),
               list(mus = matrix(rep(c(3,-3), 8) * rep((1-(4/5)^(0:(8-1))), each = 2), 
                                 byrow = T, nrow = 8),
                    Sigmas = matrix(rep(invvech(c(1,-9/10,1)), 8) * rep((4/5)^(2*(0:(8-1))), each = 4),
                                    byrow = T, nrow = 8 * 2),
                    pr = rep(1, 8)/8),
               list(mus = rbind(c(-3/2,3/2), c(3/2,-3/2), c(0,0)),
                    Sigmas = rbind(diag(9/16,2), diag(9/16,2), invvech(c(4/5,-18/25,4/5)) * 9/16),
                    pr = c(4,4,3)/11),
               list(mus = rbind(c(0,0),c(0,0)),
                    Sigmas = rbind(invvech(c(1,1,4)), invvech(c(4/9,-1/3,4/9))),
                    pr = c(2,1)/3),
               list(mus = rbind(c(1,0), c(-1,0)),
                    Sigmas = rbind(diag(4/9,2), diag(4/9,2)),
                    pr = c(1,1)/2),
               list(mus = rbind(c(-1,1), c(1,-1)),
                    Sigmas = rbind(invvech(c(4/9,4/15,4/9)), invvech(c(4/9,4/15,4/9))),
                    pr = c(1,1)/2),
               list(mus = rbind(c(1,-1), c(-1,1)),
                    Sigmas = rbind(invvech(c(4/9,14/45,4/9)), diag(4/9,2)),
                    pr = c(1,1)/2),
               list(mus = rbind(c(-1,0), c(1,2/sqrt(3)), c(1,-2/sqrt(3))),
                    Sigmas = rbind(invvech(c(9,63/10,49/4)), invvech(c(9,0,49/4)),invvech(c(9,0,49/4)))/25,
                    pr = c(3,3,1)/7),
               list(mus = matrix(c(rep(0,4),c(-1,-1,1,-1,-1,1,1,1)),byrow = T, nrow = 6),
                    Sigmas = rbind(diag(2), matrix(rep(diag(1/16,2),5), byrow = T, nrow = 5*2)),
                    pr = c(1/2,rep(1/10,5))),
               list(mus = rbind(c(-3/2,0), c(3/2,0), cbind(-1:1-3/2, -1:1), c(0,0), cbind(-1:1+3/2, -1:1)),
                    Sigmas = matrix(c(rep(invvech(c(4/9,4/15,4/9)),2),
                                      rep(invvech(c(1/15,1/25,1/15)),3)/15,
                                      invvech(c(1,3/5,1))/9,
                                      rep(invvech(c(1/15,1/25,1/15)),3)/15),
                                    byrow = T, nrow = 9*2),
                    pr = c(12/25,12/25, rep(1/350,3), 8/350, rep(1/350,3))),
               list(mus = matrix(c(rep(0,4),c(1,1,-1,1,-1,-1,1,-1)),byrow = T, nrow = 6),
                    Sigmas = rbind(diag(2),
                                   invvech(c(1,-9/10,1))/16,
                                   invvech(c(1,-9/10,1))/4,
                                   diag(1/8,2),
                                   invvech(c(1,-9/10,1))/8,
                                   diag(1/16,2)),
                    pr = c(1/2,3/40,1/5,rep(3/40,3)))
  )
  attr(mix, 'components') = length(mix$pr)
  class(mix) = 'chaconmix'
  mix
}

rmvnormix <- function(n, mix) UseMethod("rmvnormix", mix)

# generate random data from chacon mixture
rmvnormix.default <- function(n, mix) {
  ilrx = rmvnorm.mixt(n, mix$mus, mix$Sigmas, props = mix$pr)
  x = ilrInv(ilrx) # multiple -1 or not ?
  alrx = alr(x)
  attr(x, 'ilrcoord') = unclass(ilrx)
  attr(x, 'alrcoord') = unclass(alrx)
  x
}

# generate random data from ilr transformed normal mixture
rmvnormix.ilrmix <- function(n, mix) {
  ilrx = rmvnmix(n, mix)
  x = ilrInv(ilrx) # multiple -1 or not ?
  attr(x, 'ilrcoord') = ilrx
  x
}

# generate random data from alr transformed normal mixture
rmvnormix.alrmix <- function(n, mix) {
  alrx = rmvnmix(n, mix)
  x = alrInv(alrx)
  attr(x, 'alrcoord') = alrx
  x
}

dmvnormix <- function(x, mix) UseMethod("dmvnormix", mix)

dmvnormix.default <- function(x, mix) dmvnorm.mixt(unclass(x), mix$mus, mix$Sigmas, props = mix$pr)

dmvnormix.mvnmix <- function(x, mix) dmvnmix(x, mix)

# calculate density of simpledx points
# tfback = F: normal model / tfback = T: logistic normal model 
dmixS <- function(x, mix, tf, tfback = T) { # calculate density value
  n = nrow(x) 
  D = ncol(x) + 1
  x = cbind(x, 1-rowSums(x))
  i = rowSums(x > 0) == D
  if (tf == 'alr') t = alr(x[i,]) else t = ilr(x[i,])
  dt = dmvnormix(t, mix)
  if (tfback == T) {
    dt = dt/exp(rowSums(log(x[i,])))
    if (tf == 'ilr') dt = dt/sqrt(D)
  }
  dx = numeric(length = n)
  dx[i] = dt
  dx
}

simplexpanel <- function() {
  sidel = 2/sqrt(3)
  liml = 0.1
  labl = 0.06
  plot.new()
  plot.window(xlim = c(-liml,sidel+liml), ylim = c(-liml,1+liml), 
              xaxs = 'i', yaxs = 'i', asp = 1)
  
  text(sidel/2, 1+labl, expression(x[1]))
  text(-labl/sidel, -labl/2, expression(x[2]))
  text(sidel+labl/sidel, -labl/2, expression(x[3]))
  segments(0, 0, sidel, 0)
  segments(0, 0, sidel/2, 1)
  segments(sidel/2, 1, sidel, 0)
  
  margin = 0.04
  gap = c(5,10,25,50,75,90,95)/100; temp = gap/sqrt(3)
  segments(temp, gap, sidel-temp, gap, lty = 'dotted', col = 'grey')
  segments(temp, gap, 2*temp, 0, lty = 'dotted', col = 'grey')
  segments(sidel-temp, gap, sidel-2*temp, 0, lty = 'dotted', col = 'grey')
  text(temp-margin, gap, gap*100, cex = 0.6)
  text(sidel-temp+margin, gap, (1-gap)*100, cex = 0.6)
  text(2*temp, -margin, (1-gap)*100, cex = 0.6)
}

# plot normal density on real space
contourR <- function(mix, data = NULL, x1 = c(-3,3), x2 = c(-3,3), nlev = 15) {
  x1 = seq(x1[1], x1[2], l = 100)
  x2 = seq(x2[1], x2[2], l = 100)
  x = as.matrix(expand.grid(x1, x2))
  d = dmvnormix(x, mix)
  contour(x1, x2, matrix(d, nrow = length(x1), ncol = length(x2)), 
          levels = seq(0, max(d), l = nlev + 1),
          xlab = "x", ylab = "y", col = viridis(nlev))
  if (!is.null(data)) points(data, col = 'red', cex = 0.5)
}

# (Logistic) normal density on simplex
contourS <- function(mix, data = NULL, tfback = T, nlev = 15) {
  x1 = seq(0, 2/sqrt(3), l = 100)
  x2 = seq(0, 1, l = 100)
  x = as.matrix(expand.grid(x1, x2))
  x = cbind(x[,2], 1 - (sqrt(3) * x[,1] + x[,2])/2)
  
  if (class(mix)[1] == 'dirichletmix') d = ddirichletmix(x, mix)
  else {
    tf = ifelse(class(mix)[1] == 'alrmix', 'alr', 'ilr')
    d = dmixS(x, mix, tf, tfback) # tf
  }
  dm <- matrix(d, nrow = length(x1),ncol = length(x2))
  
  simplexpanel()
  contour(x1, x2, dm, levels = seq(0, max(d), l = nlev + 1),
          xlab = "x", ylab = "y", col = viridis(nlev), add = T)
  # title(main = ifelse(tfback == T, 'Logistic Normal model', 'Normal model'))
  if (!is.null(data)) {
    data = cbind((data[,3]*2 + data[,1])/sqrt(3), data[,1])
    points(data, col = 'red', cex = 0.5)
  }
}

ISE <- function(tmix, emix) UseMethod("ISE", emix)

ISE.ilrmix <- function(tmix, emix, n = 100000) {
  i = sample.int(2, n, prob = c(0.5,0.5), replace = TRUE)
  data1 = rmvnormix(sum(i==1), tmix)
  data2 = rmvnormix(sum(i==2), emix)
  data = rbind(data1, data2)
  ilrdata = rbind(attr(data1, 'ilr'), attr(data2, 'ilr'))
  
  D = ncol(data)
  jacob = exp(rowSums(log(data))) * sqrt(D)
  
  f1 = dmvnormix(ilrdata, tmix)
  f2 = dmvnormix(ilrdata, emix)
  
  f1 = f1/jacob
  f2 = f2/jacob
  ISES = mean((f1 - f2)^2/(f1/2 + f2/2))
  ISES
}

ISE.alrmix <- function(tmix, emix, n = 100000) {
  i = sample.int(2, n, prob = c(0.5,0.5), replace = TRUE)
  data1 = rmvnormix(sum(i==1), tmix)
  data2 = rmvnormix(sum(i==2), emix)
  data = rbind(data1, data2)
  ilrdata = ilr(data)
  alrdata = alr(data)
  
  D = ncol(data)
  jacob = exp(rowSums(log(data)))
  
  f1 = dmvnormix(ilrdata, tmix)
  f2 = dmvnormix(alrdata, emix)
  
  f1 = f1/(jacob*sqrt(D))
  f2 = f2/jacob
  ISES = mean((f1 - f2)^2/(f1/2 + f2/2))
  ISES
}

ISED <- function(tmix, emix) UseMethod("ISED", emix)

ISED.ilrmix <- function(tmix, emix, n = 100000) {
  i = sample.int(2, n, prob = c(0.5,0.5), replace = TRUE)
  data1 = rdirichletmix(sum(i==1), tmix)
  data2 = rmvnormix(sum(i==2), emix)
  data = rbind(data1, data2)
  ilrdata = rbind(unclass(ilr(data1)), attr(data2, 'ilr'))
  
  D = ncol(data)
  jacob = exp(rowSums(log(data))) * sqrt(D)
  
  f1 = ddirichletmix(data, tmix)
  f2 = dmvnormix(ilrdata, emix)
  
  f2 = f2/jacob
  ISES = mean((f1 - f2)^2/(f1/2 + f2/2))
  ISES
}

ISED.alrmix <- function(tmix, emix, n = 100000) {
  i = sample.int(2, n, prob = c(0.5,0.5), replace = TRUE)
  data1 = rdirichletmix(sum(i==1), tmix)
  data2 = rmvnormix(sum(i==2), emix)
  data = rbind(data1, data2)
  alrdata = rbind(unclass(alr(data1)), attr(data2, 'alr'))
  
  D = ncol(data)
  jacob = exp(rowSums(log(data)))
  
  f1 = ddirichletmix(data, tmix)
  f2 = dmvnormix(alrdata, emix)
  
  f2 = f2/jacob
  ISES = mean((f1 - f2)^2/(f1/2 + f2/2))
  ISES
}

testchacon <- function(i, n, plot = 1) {
  mix = chaconmix(i)
  data = rmvnormix(n, mix)
  ilrdata = attr(data, 'ilrcoord')
  alrdata = attr(data, 'alrcoord')
  
  imodel = mmde(ilrdata); class(imodel$mix) = c('ilrmix', class(imodel$mix))
  amodel = mmde(alrdata); class(amodel$mix) = c('alrmix', class(amodel$mix))
  
  if (plot == 1 || plot == 2) {
    if (plot == 1) data = ilrdata = alrdata = NULL # no points in plot
    
    par(mfrow=c(2,2), mar=c(2.5,2.5,1,0.5))
    contourR(mix, ilrdata)
    contourR(imodel$mix, ilrdata)
    contourR(amodel$mix, alrdata)
    par(mfrow=c(1,1))
    
    par(mfrow=c(2,2), mar=rep(0,4))
    contourS(mix, data)
    contourS(imodel$mix, data)
    contourS(amodel$mix, data)
    par(mfrow=c(1,1))
    
    par(mfrow=c(2,2), mar=rep(0,4))
    contourS(mix, data, tfback = F)
    contourS(imodel$mix, data, tfback = F)
    contourS(amodel$mix, data, tfback = F)
    par(mfrow=c(1,1))
  }
  
  list('ilr' = ISE(mix, imodel$mix), 'alr' = ISE(mix, amodel$mix))
}

testbounday <- function(i, n, h, pr, plot = 1) {
  mix = boundarymix(i, h, pr)
  data = rdirichletmix(n, mix)
  # r = which(data == 0, arr.ind = T)[,1]
  ilrdata = unclass(ilr(data))
  alrdata = unclass(alr(data))
  
  # ilrdata1 = ilrdata[-r,]
  imodel = mmde(ilrdata); class(imodel$mix) = c('ilrmix', class(imodel$mix))
  amodel = mmde(alrdata); class(amodel$mix) = c('alrmix', class(amodel$mix))
  
  if (plot == 1 || plot == 2) {
    if (plot == 1) data = ilrdata = alrdata = NULL # no points in plot
    par(mfrow=c(2,2), mar=rep(0,4))
    contourS(mix, data)
    contourS(imodel$mix, data)
    contourS(amodel$mix, data)
    par(mfrow=c(1,1))
  }
  
  list('ilr' = ISED(mix, imodel$mix), 'alr' = ISED(mix, amodel$mix))
}

testchacon(9, 200, plot = 1)
testbounday(7, 300, h=c(0.3, 0.2), pr=c(0.4,0.6), plot = 2)




# n = 250
# B = 10
# iseresult = list()
# for (j in 1:12) {
#   r = matrix(nr = B, nc = 2)
#   colnames(r) = c('ilr','alr')
#   for (i in 1:B) {
#     print(c(j, i))
#     ISEvalue = testchacon(j, n, plot = 0)
#     r[i, 'ilr'] = ISEvalue$ilr
#     r[i, 'alr'] = ISEvalue$alr
#   }
#   iseresult[[j]] = r
# }

# 
# ISEresult = read.csv('output.csv')
# ISElist = split(ISEresult, ISEresult$mixture)
# par(mfrow=c(4,3))
# for (i in ISElist) {
#   boxplot(i[,-3])
# }
# par(mfrow=c(1,1))

# alrdata = transformData(iInv_data, 'alr')
# alrmus = alr(-1*ilrInv(mix$mus))
# clrvar = ilrvar2clr(mix$Sigmas)
# F = cbind(diag(2), rep(-1, 2)) # D = 3
# alrvar = F%*%clrvar%*%t(F)
# alrmix = list(mus = as.matrix(alrmus), Sigmas = alrvar, pr = mix$pr)