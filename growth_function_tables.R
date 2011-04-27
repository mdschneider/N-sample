#
# Compute scale factor versus Omega_m at matched linear growth
#
#  Author: Michael D. Schneider <mdschneider@me.com>
#
require(plyr)
source("cosmology.R")

p.default <- InitCosmologicalParameters()
# --- SCDM cosmology ----------
#p.default$Omega.m <- 1
#p.default$Omega.L <- 0
# -----------------------------
cat(sprintf("Default cosmological parameters:\n"))
print(p.default)

D.norm <- LinearGrowthFunction(1/21, p.default)
kF <- 2*pi/1000

kScaleFactorICs <- 1/21

#' Table of the linear growth function vs. Omega_m at fixed a
D1.omegam.table <- function(om.min=0.01, om.max=0.7, a=1) {
  om <-  seq(om.min, om.max, length.out=500)
  p <- InitCosmologicalParameters()
  D.om <- aaply(om, 1, function(x){p$Omega.m=x; LinearGrowthFunction(a, p) / LinearGrowthFunction(kScaleFactorICs, p)})
  return(data.frame(omega.m=om, D=D.om))
}

#' Table of the linear growth function vs. a
D1.a.table <- function(a.min=kScaleFactorICs, a.max=25) {
  a <- seq(a.min, a.max, length.out=500)
  p <- InitCosmologicalParameters()
  Da <- aaply(a, 1, function(x) {LinearGrowthFunction(x, p) / D.norm})
  return(data.frame(a=a, D=Da))
}

ComputeTables <- function() {
  Dm <- D1.omegam.table()
  Da <- D1.a.table()
  aspl <- splinefun(Da$D, Da$a)
  a.om <- aspl(Dm$D)
  return(data.frame(a=Da$a, a.om=a.om, D.a=Da$D, om=Dm$omega.m, D.om=Dm$D))
}

ComputeInverseTable <- function(Omega.min=0.01, Omega.max=0.5, a.min=0.1, a.max=100) {
  D.target <- LinearGrowthFunction(1, p.default) / LinearGrowthFunction(kScaleFactorICs, p.default)
  Omega.list <- seq(Omega.min, Omega.max, length.out=99)
  a.list <- c()
  p <- InitCosmologicalParameters()
  for (om in Omega.list) {
    a <- 10^seq(log10(a.min), log10(a.max), length.out=200)
    p$Omega.m <- om
    Da <- aaply(a, 1, function(x) {LinearGrowthFunction(x, p) / LinearGrowthFunction(kScaleFactorICs, p)})
    aspl <- splinefun(Da, a)
    a.new <- aspl(D.target)
    a.list <- c(a.list, a.new)
  }
  return(data.frame(a=a.list, Omega.m=Omega.list))
}

#' Create fig. 1 for the paper
PlotGrowthTables <- function(n.sigma=4.3, kthresh=8*kF) {
  tb <- ComputeTables()
  aspl <- splinefun(tb$om, tb$a.om)
  sigma <- sqrt(MassVarianceTopHatSpheres(2*pi/(256*kF), ln.kmin=log(kF), ln.kmax=log(kthresh)))
  omega.range <- p.default$Omega.m * c(1 - n.sigma*sigma, 1, 1 + n.sigma*sigma)
  a.range <- aspl(omega.range)
  print(a.range)
  print(sigma)
  print(omega.range)
  par(cex.lab=1.5)
  layout(matrix(c(1,2,0,3), nrow=2, ncol=2, byrow=TRUE))
  # D vs. a
  par(mar=c(0, 4.1, 4.1, 0))
  plot(tb$a, tb$D.a, type='l', xlim=c(0,6), ylim=c(0,25), lwd=3,
       xaxt='n',
       xlab="", ylab="")
  axis(side=1, at=seq(0,5), labels=seq(0,5), tick=TRUE)
  mtext("a", side=1, line=3, cex=1.3)
  mtext(expression(D(a, Omega[m]^{0}) / D(a==a[IC], Omega[m]^{0})), side=2, line=2, cex=1.3)
  # D vs. Omega_m
  par(mar=c(0, 0, 4.1, 4.1))
  plot(tb$om, tb$D.om, type='l', xlim=c(0,0.5), ylim=c(0,25), lwd=3,
       xaxt='n', yaxt='n',
       #xlab=expression(Omega[m]),
       xlab="", ylab="")
  abline(v=omega.range, lty=c(3,2,3), col=4, lwd=2)
  axis(side=4, at=seq(0, 25, by=5), labels=seq(0,25,by=5), tick=TRUE)
  axis(side=3, at=seq(0,5)/10, labels=seq(0,5)/10)
  mtext(expression(D(a==1, Omega[m]) / D(a==a[IC], Omega[m])), side=4, line=3, cex=1.3)
  # a vs. Omega_m
  par(mar=c(5.1, 0, 0, 4.1))
  plot(tb$om, tb$a.om, type='l', ylim=c(0,6), xlim=c(0,0.5), lwd=3,
       yaxt='n',
       ylab=expression(a),
       xlab=expression(Omega[m]))
  abline(v=omega.range, lty=c(3,2,3), col=4, lwd=2)
  abline(h=a.range, lty=c(3,2,3), col=3, lwd=2)
  axis(side=2, at=seq(0,5), labels=seq(0,5), tick=TRUE)
  mtext("a", side=2, line=3, cex=1.3)
}

#' Plot |a'-a| vs. Lbox with
#' ntresh = 8, delta x = 1000/256, and n.sigma = 4.3
#'
#' figure 2 for paper
PlotScaleFactorRangeVSLbox <- function(delta.x=1000/256, nthresh=8, n.sigma=4.3,
                                Omega.growth.min=0.156659) {
  tb <- ComputeTables()
  aspl <- splinefun(tb$om, tb$a.om)
  tbI <- ComputeInverseTable()
  ndx <- which(tbI$Omega.m >= 0.16)
  aspl.inv <- splinefun(tbI$Omega.m[ndx], tbI$a[ndx])

  layout(matrix(c(1,2), nrow=2, ncol=1))

  # ----- Loop over Lbox and compute maximum a' and |a'-a|
  lbox <- seq(1000, 3500, length.out=100)
  a.diff <- rep(0, length(lbox))
  a.inv.max <- rep(Inf, length(lbox))
  omega.min <- c()
  for (i in seq(length(lbox))) {
    kF <- 2*pi / lbox[i]
    kthresh <- nthresh * kF
    sigma <- sqrt(MassVarianceTopHatSpheres(delta.x, ln.kmin=log(kF), ln.kmax=log(kthresh)))
    omega.range <- p.default$Omega.m * c(1 - n.sigma*sigma, 1 + n.sigma*sigma)
    a.range <- aspl(omega.range)
    a.diff[i] <- max(abs(a.range-1))
    if (omega.range[1] > Omega.growth.min) {
      a.inv.max[i] <- aspl.inv(omega.range[1])
    }
    omega.min <- c(omega.min, omega.range[1])
  }

  # ----- Top panel: max|a'-a| vs. Lbox
  par(mar=c(0, 4.1, 4.1, 2.1))
  plot(lbox/1000, a.diff, type='l', ylim=c(0.1,3.5), log='y', lwd=3,
       xaxt='n',
       ylab="")
  mtext(expression(max(group("|",a^{"'"}-a,"|"))), side=2, line=2, cex=1.5)
  #abline(h=0.2, lty=2, lwd=2)
  text(x=2.6, y=2, labels="Mode addition", cex=1.5)

  # ----- Maximum a' for mode-subtraction growth matching
  par(mar=c(5.1, 4.1, 0, 2.1))
  plot(lbox/1000, a.inv.max, type='l', col=1, lty=1, lwd=3,
       ylim=c(0, 7),
       xlab="", ylab="")
  mtext(expression(L[box]~~(h^{-1}*Gpc)), side=1, line=3, cex=1.5)
  mtext("max(a')", side=2, line=2, cex=1.5)
  lspl <- splinefun(omega.min, lbox)
  text(x=2.6, y=6, labels="Mode subtraction", cex=1.5)
  cat(sprintf("minimum box size: %8.6f Gpc/h\n", lspl(Omega.growth.min)/1000))
  abline(v=lspl(Omega.growth.min)/1000, lty=2, lwd=2)
  abline(h=1, lty=1, col=1)
}

#' Make fig. 6 for paper
#'
PlotGrowthTablesInverse <- function(Omega.min=0.01, Omega.max=0.5, kthresh=8*kF,
                                    a.min=0.1, a.max=100,
                                    n.sigma=4.3,
                                    Omega.growth.min=0.156659) {
  D.target <- LinearGrowthFunction(1, p.default) / LinearGrowthFunction(kScaleFactorICs, p.default)
  cat(sprintf("D.target: %8.5g\n", D.target))

  layout(matrix(c(1,2), nrow=1, ncol=2))
  # a vs. D
  par(mar=c(5.1, 4.1, 4.1, 0), cex.lab=1.5)
  a <- 1
  plot(D.target, a, type='n', log='y',
       xaxt='n',
       xaxs='i', yaxs='i',
       xlim=c(0, 25), ylim=c(a.min, a.max),
       xlab=expression(D(a, Omega[m]) / D(a[IC], Omega[m])),
       ylab="a")
  
  Omega.list <- seq(Omega.min, Omega.max, length.out=99)
  a.list <- c()
  p <- p.default
  j <- 0
  for (om in Omega.list) {
    cat(sprintf("Omega_m: %5.4f\n", om))
    a <- 10^seq(log10(a.min), min(2, log10(8/om^2)), length.out=200)
    print(range(a))
    p$Omega.m <- om
    Da <- aaply(a, 1, function(x) {LinearGrowthFunction(x, p) / LinearGrowthFunction(kScaleFactorICs, p)})
    if (j %% 10 == 0) {
      lines(Da, a, lwd=2)
      print(om)
    }
    aspl <- splinefun(Da, a)
    a.new <- aspl(D.target)
    a.list <- c(a.list, a.new)
    j <- j + 1
  }
  abline(v=D.target, col=2, lwd=2)
  axis(side=1, at=seq(0,20, by=5), labels=seq(0,20, by=5))
  
  sigma <- sqrt(MassVarianceTopHatSpheres(2*pi/(256*kF), ln.kmin=log(kF), ln.kmax=log(kthresh), p=p.default))
  cat(sprintf("sigma: %8.6g\n", sigma))
  omega.range <- p.default$Omega.m * c(1 - n.sigma*sigma, 1, 1 + n.sigma*sigma)
  ndx <- which(Omega.list >= Omega.growth.min)
  aspl <- splinefun(Omega.list[ndx], a.list[ndx])
  a.range <- aspl(omega.range)
  abline(h=a.range, lty=c(3,2,3), col=3, lwd=2)
  print(omega.range)
  print(a.range)

  print(data.frame(a=a.list, omega=Omega.list))

  # a vs. Omega
  par(mar=c(5.1, 0, 4.1, 2.1))
  plot(Omega.list, a.list, type='l', yaxt='n', lwd=3, log='y',
       xaxs='i', yaxs='i',
       ylim=c(a.min, a.max),
       xlab=expression(Omega[m]))
  abline(v=omega.range, lty=c(3,2,3), col=4, lwd=2)
  abline(h=a.range, lty=c(3,2,3), col=3, lwd=2)
}
