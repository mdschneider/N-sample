MODULE types_nsample
  ! ----------------------------------------------------------------------
  ! This module defines types for the covariance resampling codes.
  ! It is modeled after healpix_types.F90 from the Healpix package.
  ! ----------------------------------------------------------------------
  INTEGER, PARAMETER, public :: i4b = SELECTED_INT_KIND(9)
#ifdef NO64BITS
  INTEGER, PARAMETER, public :: i8b = i4b
#else
  INTEGER, PARAMETER, public :: i8b = SELECTED_INT_KIND(16)
#endif
!!$  INTEGER, PARAMETER, public :: sp  = SELECTED_REAL_KIND(5,30)
  INTEGER, PARAMETER, public :: sp  = kind(1.0_4) !SELECTED_REAL_KIND(5,30)
!!$  INTEGER, PARAMETER, public :: dp  = SELECTED_REAL_KIND(12,200)
!!$  INTEGER, PARAMETER, public :: dp  = SELECTED_REAL_KIND(16,308)
  INTEGER, PARAMETER, public :: dp  = kind(1.0_8) 
  INTEGER, PARAMETER, public :: lgt = KIND(.TRUE.)
  INTEGER, PARAMETER, public :: csp = KIND((1.0_sp, 1.0_sp))
  INTEGER, PARAMETER, public :: cdp = KIND((1.0_dp, 1.0_dp))

  ! Numerical Constant (Double precision)
  REAL(kind=dp), PARAMETER, public :: PI    = 3.141592653589793238462643383279502884197_dp
  REAL(kind=dp), PARAMETER, public :: TWOPI = 6.283185307179586476925286766559005768394_dp
  REAL(kind=dp), PARAMETER, public :: SQRT2 = 1.41421356237309504880168872420969807856967_dp

  real(kind=DP), parameter, public :: RAD2DEG = 180.0_DP / PI
  real(kind=DP), parameter, public :: DEG2RAD = PI / 180.0_DP

  integer, parameter :: len_filename = 200

END MODULE types_nsample
