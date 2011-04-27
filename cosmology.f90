!> Routines relating to the cosmological model.
module cosmology
  use types_nsample
  implicit none

#ifdef SCDM
  real(dp), parameter :: Omega_m = 1.d0
  real(dp), parameter :: Omega_L = 0.d0
#else ! takahashi
  real(dp), parameter :: Omega_m = 0.238_dp
  real(dp), parameter :: Omega_L = 0.762_dp
#endif
  real(dp), parameter :: hubble = 0.732_dp

contains
	!-------------------------------------------------------------
	!> Normalized Hubble parameter as a function of scale factor.
	!! @param a Scale factor
	!<------------------------------------------------------------	
  pure function H_H0(a)
    real(dp), intent(in)   :: a
    real(dp)               :: H_H0
    H_H0 = sqrt(Omega_m/a**3 + Omega_L)
  end function H_H0

	!-------------------------------------------------------------
	!> Linear growth function as a function of scale factor.
	!! @param a Scale factor
	!<------------------------------------------------------------	
  function lingrowth(a)
    real(dp), intent(in)   :: a
    real(dp)               :: lingrowth
    lingrowth = rombint(f_D1,0._dp,a,1.d-12)*2.5*Omega_m*H_H0(a)
  end function lingrowth

	!-------------------------------------------------------------
	!> Logarithmic derivative of the linear growth function 
	!! with respect to the log scale factor.
	!! @param a Scale factor
	!<------------------------------------------------------------	
  function dlnDdlna(a)
    real(dp), intent(in) :: a
    real(dp)             :: dlnDdlna
    dlnDdlna = (-3.0+5.0*a/lingrowth(a))*(Omega_m/a**3)/(2.*H_H0(a)**2)
  end function dlnDdlna

	!-------------------------------------------------------------
	!> Logarithmic derivative of the linear growth function 
	!! with respect to the log matter density.
	!! @param a Scale factor
	!! @param
	!<------------------------------------------------------------	
  function dlnDdlnOmega(a)
    real(dp), intent(in) :: a
    real(dp)             :: dlnDdlnOmega
    real(dp)             :: I,G,term2,term3
    I = rombint(f_dD1dOm,0._dp,a,1.d-12)
    G = H_H0(a)
    term2 = Omega_m/(2.*(Omega_m+Omega_L*a**3))
    term3 = I*3.75*Omega_m**2*G/lingrowth(a)
    dlnDdlnOmega = 1._dp + term2 - term3
  end function dlnDdlnOmega

  pure function f_D1(a)
    real(dp), intent(in) :: a
    real(dp)             :: f_D1
    f_D1 = (a/(Omega_m + Omega_L*a**3))**1.5
  end function f_D1

  pure function f_dD1dOm(a)
    real(dp), intent(in) :: a
    real(dp)             :: f_dD1dOm
    real(dp)             :: G
    f_dD1dOm = a**1.5 / (Omega_m+Omega_L*a**3)**2.5
  end function f_dD1dOm


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function rombint(f,a,b,tol)
    !  Rombint returns the integral from a to b of using Romberg integration.
    !  The method converges provided that f(x) is continuous in (a,b).
    !  f must be real(dp) and must be declared external in the calling
    !  routine.  tol indicates the desired relative accuracy in the integral.
    !
    implicit none
    integer, parameter :: MAXITER=20
    integer, parameter :: MAXJ=5
    dimension g(MAXJ+1)
    real(dp) f
    external f
    real(dp) :: rombint
    real(dp), intent(in) :: a,b,tol
    integer :: nint, i, k, jmax, j
    real(dp) :: h, gmax, error, g, g0, g1, fourj
    !

    h=0.5d0*(b-a)
    gmax=h*(f(a)+f(b))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
10  i=i+1
    if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
         go to 40
    !  Calculate next trapezoidal rule approximation to integral.
    g0=0._dp
    do 20 k=1,nint
       g0=g0+f(a+(k+k-1)*h)
20     continue
       g0=0.5d0*g(1)+h*g0
       h=0.5d0*h
       nint=nint+nint
       jmax=min(i,MAXJ)
       fourj=1._dp
       do 30 j=1,jmax
          !  Use Richardson extrapolation.
          fourj=4._dp*fourj
          g1=g0+(g0-g(j))/(fourj-1._dp)
          g(j)=g0
          g0=g1
30        continue
          if (abs(g0).gt.tol) then
             error=1._dp-gmax/g0
          else
             error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
          go to 10
40        rombint=g0
          if (i.gt.MAXITER.and.abs(error).gt.tol)  then
             write(*,*) 'Warning: Rombint failed to converge; '
             write (*,*)'integral, error, tol:', rombint,error, tol
          end if

        end function rombint
end module cosmology
