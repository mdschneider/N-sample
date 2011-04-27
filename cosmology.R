require(splines)

## Constants
#' Hubble distance in Mpc/h
#' @export
kcOverH0 <- 2.99792458e3 # Mpc/h
#' Speed of light in Mpc/s
#' @export
kcMpcPerSec <- 9.71561189e-15 # Speed of light in Mpc/s
#' Hubble time in s/h
#' @export
ktH0 <- kcOverH0 / kcMpcPerSec # s/h
#' Critical density in (M_sun/h) (h/Mpc)^3
#' @export
krhoc <- 2.77519738e11   # (M_sun/h) (h/Mpc)^3


## Cosmological parameters
kDefaultCosmoParams <- list(h = 0.732,
                            Omega.m = 0.238,
                            Omega.L = 0.762,
                            Omega.b = 0.041,
                            ns = 0.958,
                            sigma.8 = 0.76)

#' Create new object of class "cosmoParams" with default cosmological parameter values
#'
#' @return Object of class "cosmoParams"
#' @export
InitCosmologicalParameters <- function() {
  p <- kDefaultCosmoParams
  class(p) <- c(class(p), "cosmoParams")
  return(p)
}

#' Pretty printing for object of class "cosmoParams"
#'
#' @return None
#' @export
print.cosmoParams <- function(p) {
  cat(sprintf("Hubble constant:      %5.3f km/s/Mpc\n", p$h*100))
  cat(sprintf("Matter density:       %5.3f\n", p$Omega.m))
  cat(sprintf("Dark energy density:  %5.3f\n", p$Omega.L))
  cat(sprintf("Baryon density:       %5.3f\n", p$Omega.b))
  cat(sprintf("Scalar index:         %5.3f\n", p$ns))
  cat(sprintf("sigma_8:              %5.3f\n", p$sigma.8))
}


# ------------------------------------------------------------
# Distance redshift relations
#
# units: h^(-1) Mpc
# ------------------------------------------------------------
#' Reduced Hubble parameter as a function of scale factor
#'
#' @param a Scale factor
#' @param p Object of class "cosmoParams"
#' @export
ReducedHubbleParameter <- function(a, p)
{
  return(sqrt( (p$Omega.m / a^3) +
              ((1 - p$Omega.m - p$Omega.L) / a^2) + p$Omega.L))
}

#' Redshift dependent matter density
OmegaM <- function(a, p) {
  return(p$Omega.m/a^3 / ReducedHubbleParameter(a, p)^2)
}

#' Comoving line of sight distance in Mpc/h
#'
#' Assumes Omega_k == 0.
#' @param z1 Initial redshift
#' @param z2 Final redshift
#' @param p Object of class "cosmoParams"
#' @export
DistanceLOS <- function(z1, z2, p, rel.tol=.Machine$double.eps^0.8, 
												abs.tol=.Machine$double.eps^0.8)
{
  I.chi <- integrate(f.chi, z1, z2, p=p,
                     rel.tol=rel.tol, abs.tol=abs.tol)
  return(I.chi$value * kcOverH0)
}

f.chi <- function(z, p)
{
  return( 1.0/ReducedHubbleParameter(1.0/(1.0+z), p) )
}

#' Luminosity distance from redshift 0 to z in Mpc/h
#'
#' Assumes Omega_k == 0.
#' @param z redshift
#' @param p Object of class "cosmoParams"
#' @author Michael D. Schneider <schneider42@llnl.gov>
#' @export
LuminosityDistance <- function(z,p) {
	return(DistanceLOS(0, z, p) * (1+z))
}

#' Angular diameter distance from redshift 0 to z
#'
#' Assumes Omega_k == 0.
#' @param z redshift
#' @param p Object of class "cosmoParams"
#' @author Michael D. Schneider <schneider42@llnl.gov>
#' @export
AngularDiameterDistance <- function(z, p) {
	return(DistanceLOS(0, z, p) / (1+z))
}

#' Lookback time from redshift z1 to redshift z2
#'
#' Assumes Omega_k == 0.
#' If z1=0, this is the usual lookback time.  
#' If z2=Inf this gives the age of the universe at z1.
#' @param z1 lower redshift
#' @param z2 upper redshift
#' @param pObject of class "cosmoParams"
#' @returnType numeric
#' @return The lookback time from z1 to z2
#' @author Michael D. Schneider <schneider42@llnl.gov>
#' @export
LookBackTime <- function(z1=0, z2=Inf, p, rel.tol=.Machine$double.eps^0.9, 
												 abs.tol=.Machine$double.eps^0.9) {
  I.t <- integrate(f.t, z1, z2, p=p, 
                   rel.tol=rel.tol, abs.tol=abs.tol)
  return(I.t$value * ktH0)
}

f.t <- function(z, p) {
  return(1 / ((1 + z) * ReducedHubbleParameter(1/(1+z), p))) 
}

ComovingVolume <- function(z.max=10, p) {
	I <- integrate(dV, 0, z.max, p=p, rel.tol=.Machine$double.eps^0.85, subdivisions=1000)
	return(I$value * 4 * pi)
}

dV <- function(z, p) {
	chi <- DistanceLOS(0, z, p=p)
	return(chi^2 * f.chi(z, p=p) * kcOverH0)
}

# ------------------------------------------------------------
# (unnormalized) linear growth function at scale factor a
# ------------------------------------------------------------

#' Linear growth function at scale factor a
#'
#' @param a Scale factor
#' @param p Object of class "cosmoParams"
LinearGrowthFunction <- function(a, p, rel.tol=1.e-10, abs.tol=1.e-10)
{
  I.growth <- integrate(f.growth, 0, a, p=p,
                        rel.tol=rel.tol, abs.tol=abs.tol)
  return(ReducedHubbleParameter(a,p) * I.growth$value * 2.5 * p$Omega.m)
}

f.growth <- function(a,p)
{
  return( (a / (p$Omega.m + (1 - p$Omega.m - p$Omega.L) * a +
              p$Omega.L * a^3))^1.5 )
}


# ------------------------------------------------------------
# Mass variance in top-hat spheres (for calculating e.g. sigma.8)
# ------------------------------------------------------------

## Fourier transform of top-hat window
TopHat <- function(k, R)
{
  y <- k*R
  return( (3/y^3) * (sin(y) - y*cos(y)) )
}

#' Mass variance in a top-hat sphere of radius R
#'
#' @param R Smoothing scale in Mpc/h
#' @param p Object of class "cosmoParams" (will use default params if neither p nor ps.spl are supplied)
#' @param ps.spl Spline in (ln k, ln Delta^2) as output from PowerSpectrumSpline (will generate from \link{PowerSpectrumDefault} if not supplied)
#' @param rel.tol Relative tolerance for integration of windowed power spectrum
#' @param abs.tol Absolute tolerance for integration of windowed power spectrum
#' @author Michael D. Schneider <mdschneider@me.com>
#' @returnType numeric
#' @return Mass variance at smoothing scale \code{R}
#' @export
MassVarianceTopHatSpheres <- function(R, p=NULL, ps.spl=NULL,
                                      ln.kmin=-4, ln.kmax=2,
                                      rel.tol=1.e-6, abs.tol=1.e-6)
{
  if (is.null(ps.spl)) {
    ps.spl <- PowerSpectrumDefault(p=p, return.spline=TRUE)
  }
  sigmasq <- c()
  for (i in 1:length(R)) {
    I.sigmasq <- integrate(f.sigma, ln.kmin, ln.kmax, R=R[i], ps.spl=ps.spl,
                           rel.tol=rel.tol, abs.tol=abs.tol)
    sigmasq <- c(sigmasq, I.sigmasq$value)
  }
  return(sigmasq)
}

f.sigma <- function(lnk, R, ps.spl)
{
  lnpk <- predict(ps.spl, lnk)$y
  return( exp(lnpk) * (TopHat(exp(lnk), R))^2 )
}

#' Create a spline of (ln k, ln Delta^2) given k, P(k)
#'
#' @param k Array of wavenumber values in h/Mpc
#' @param pk Array of power spectrum values at \code{k}
#' @author Michael D. Schneider <mdschneider@me.com>
#' @returnType npolySpine
#' @return Spline object output from \link{\code{interpSpline}}
#' @export
PowerSpectrumSpline <- function(k, pk)
{
  return(interpSpline(log(k), log(pk*k^3/(2.*pi^2))))
}

DeltasqMatterPower <- function(k, p, Transfer) {
  (k*kcOverH0)^(p$ns+3) * Transfer^2 * 1e-9
}

#' Default model for the matter power spectrum
#'
#' @param k Wavenumber in h/Mpc
#' @param p Object of class "cosmoParams" (will use default params if not supplied)
#' @param return.spline Logical.  Return the power spectrum evaluated at input wavenumbers
#' or a spline function?
#' @author Michael D. Schneider <mdschneider@me.com>
#' @returnType numeric or npolySpline
#' @return Matter power spectrum evaluated at wavenumber(s) \code{k}.
#' @export
PowerSpectrumDefault <- function(k=NULL, p=NULL, return.spline=FALSE) {
  if (is.null(p)) {
    p <- InitCosmologicalParameters()
  }
  # Set amplitude
  k.sigma8 <- exp(seq(-4,2, length.out=100))
  T.sigma8 <- TransferFunctionBBKS(k.sigma8, p)
  Deltasq.sigma8 <- DeltasqMatterPower(k.sigma8, p, T.sigma8)
  ps.spl <- PowerSpectrumSpline(k.sigma8, Deltasq.sigma8 * (2*pi^2)/k.sigma8^3)
  sigma8sq <- MassVarianceTopHatSpheres(R=8, ps.spl=ps.spl)
  if (return.spline) {
    return(PowerSpectrumSpline(k.sigma8, Deltasq.sigma8*p$sigma.8^2/sigma8sq*(2*pi^2)/k.sigma8^3))
  }
  else {
    if (is.null(k)) {
      stop("ERROR: need value for k in PowerSpectrumDefault")
    }
    # Return P(k) at input k
    Transfer <- TransferFunctionBBKS(k, p)
    Deltasq <- DeltasqMatterPower(k, p, Transfer) * p$sigma.8^2 / sigma8sq
    return(Deltasq * (2*pi^2) / k^3)
  }
}

# ------------------------------------------------------------
# Spherical collapse
# ------------------------------------------------------------
#' Linear overdensity for spherical collapse
#'
#' @param z Redshift
#' @param p Object of class "cosmoParams"
#' @author Michael D. Schneider <mdschneider@me.com>
#' @returnType numeric
#' @return delta_sc(z)
#' @export
deltaSC <- function(z, p) {
  a <- 1/(1+z)
  return(1.68647*(OmegaM(a, p))**0.0055)
}

#' Halo mass where linear spherical collapse overdensity equals mass variance in top-hat spheres
#'
#' @param z Redshift
#' @param p Cosmological parameters (with class "cosmoParams")
#' @author Michael D. Schneider <mdschneider@me.com>
#' @returnType numeric
#' @return Mstar(z)
#' @export
Mstar <- function(z, p=NULL, ps.spl=NULL) {
  if (is.null(p)) {
    p <- InitCosmologicalParameters()
  }
  dsc <- deltaSC(z, p)
  fsigma <- function(log10M) {
    R <- (3*10^log10M/(4*pi) / (krhoc * p$Omega.m * (1+z)^3))^(1/3)
    sigma <- sqrt(MassVarianceTopHatSpheres(R, p=p, ps.spl=ps.spl))
    return(dsc - sigma)
  }
  res <- uniroot(fsigma, lower=9, upper=15, tol=.Machine$double.eps^0.8, maxiter=10000)
  return(10^res$root)
}

# ------------------------------------------------------------
# Transfer functions
# ------------------------------------------------------------
#' BBKS model for the matter transfer function
#'
#' @param k Wavenumber in h/Mpc (scalar or array values supported)
#' @param p Cosmological parameters
#' @author Michael D. Schneider <mdschneider@me.com>
#' @returnType numeric
#' @return The BBKS matter transfer function evaluated at k
#' @export
TransferFunctionBBKS <- function(k, p) {
  q <- k / (p$Omega.m * p$h)
  A <- log(1+2.34*q)/(2.34*q)
  return(A*(1 + 3.89*q + (16.2*q)^2 + (5.47*q)^3 + (6.71*q)^4)^(-0.25))
}

