!> Handling of FFT plans and grid indexing
module fft_nsample
  use types_nsample
  use parameters
  implicit none

#ifndef MPI
#ifdef FFTW ! fftw3.f
  INTEGER FFTW_R2HC
  PARAMETER (FFTW_R2HC=0)
  INTEGER FFTW_HC2R
  PARAMETER (FFTW_HC2R=1)
  INTEGER FFTW_DHT
  PARAMETER (FFTW_DHT=2)
  INTEGER FFTW_REDFT00
  PARAMETER (FFTW_REDFT00=3)
  INTEGER FFTW_REDFT01
  PARAMETER (FFTW_REDFT01=4)
  INTEGER FFTW_REDFT10
  PARAMETER (FFTW_REDFT10=5)
  INTEGER FFTW_REDFT11
  PARAMETER (FFTW_REDFT11=6)
  INTEGER FFTW_RODFT00
  PARAMETER (FFTW_RODFT00=7)
  INTEGER FFTW_RODFT01
  PARAMETER (FFTW_RODFT01=8)
  INTEGER FFTW_RODFT10
  PARAMETER (FFTW_RODFT10=9)
  INTEGER FFTW_RODFT11
  PARAMETER (FFTW_RODFT11=10)
  INTEGER FFTW_FORWARD
  PARAMETER (FFTW_FORWARD=-1)
  INTEGER FFTW_BACKWARD
  PARAMETER (FFTW_BACKWARD=+1)
  INTEGER FFTW_MEASURE
  PARAMETER (FFTW_MEASURE=0)
  INTEGER FFTW_DESTROY_INPUT
  PARAMETER (FFTW_DESTROY_INPUT=1)
  INTEGER FFTW_UNALIGNED
  PARAMETER (FFTW_UNALIGNED=2)
  INTEGER FFTW_CONSERVE_MEMORY
  PARAMETER (FFTW_CONSERVE_MEMORY=4)
  INTEGER FFTW_EXHAUSTIVE
  PARAMETER (FFTW_EXHAUSTIVE=8)
  INTEGER FFTW_PRESERVE_INPUT
  PARAMETER (FFTW_PRESERVE_INPUT=16)
  INTEGER FFTW_PATIENT
  PARAMETER (FFTW_PATIENT=32)
  INTEGER FFTW_ESTIMATE
  PARAMETER (FFTW_ESTIMATE=64)
  INTEGER FFTW_ESTIMATE_PATIENT
  PARAMETER (FFTW_ESTIMATE_PATIENT=128)
  INTEGER FFTW_BELIEVE_PCOST
  PARAMETER (FFTW_BELIEVE_PCOST=256)
  INTEGER FFTW_NO_DFT_R2HC
  PARAMETER (FFTW_NO_DFT_R2HC=512)
  INTEGER FFTW_NO_NONTHREADED
  PARAMETER (FFTW_NO_NONTHREADED=1024)
  INTEGER FFTW_NO_BUFFERING
  PARAMETER (FFTW_NO_BUFFERING=2048)
  INTEGER FFTW_NO_INDIRECT_OP
  PARAMETER (FFTW_NO_INDIRECT_OP=4096)
  INTEGER FFTW_ALLOW_LARGE_GENERIC
  PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
  INTEGER FFTW_NO_RANK_SPLITS
  PARAMETER (FFTW_NO_RANK_SPLITS=16384)
  INTEGER FFTW_NO_VRANK_SPLITS
  PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
  INTEGER FFTW_NO_VRECURSE
  PARAMETER (FFTW_NO_VRECURSE=65536)
  INTEGER FFTW_NO_SIMD
  PARAMETER (FFTW_NO_SIMD=131072)
  INTEGER FFTW_NO_SLOW
  PARAMETER (FFTW_NO_SLOW=262144)
  INTEGER FFTW_NO_FIXED_RADIX_LARGE_N
  PARAMETER (FFTW_NO_FIXED_RADIX_LARGE_N=524288)
  INTEGER FFTW_ALLOW_PRUNING
  PARAMETER (FFTW_ALLOW_PRUNING=1048576)
  INTEGER FFTW_WISDOM_ONLY
  PARAMETER (FFTW_WISDOM_ONLY=2097152)
#endif FFTW
#else ! MPI
  public
  include "fftw_f77.i"

  ! ----- FFTW plans and slab decompisition dimensions
  integer(i8b)   :: planf, planr, planc2c
  integer(i4b)   :: local_nlast2_after_transform,local_last2_start_after_transform
  integer(i4b)   :: total_local_size
  !
  integer(i8b)   :: planfL, planrL, planc2cL
  integer(i4b)   :: local_nlast2_after_transformL,local_last2_start_after_transformL
  integer(i4b)   :: total_local_sizeL
#endif MPI
  integer(i4b)   :: local_nlast, local_last_start
  integer(i4b)   :: local_nlastL, local_last_startL

contains

	!-------------------------------------------------------------
	!> Return the FFT frequency in units of the fundamental frequency
  !! for an N-length FFT and array index i.
  !!
  !! 1 <= i <= N
  !! For 1 <= i <= N/2      return   0 <= fftfreq <= N/2 - 1
  !! For N/2 + 1 <= i <= N  return   -N/2 <= fftfreq <= -1
	!! @param i Grid index 
	!! @param N Number of grid points per dimension.
	!<------------------------------------------------------------	
  pure function fftfreq(i,N)
    integer(i4b), intent(in) :: i,N
    integer(i4b)             :: fftfreq
    integer(i4b)             :: im
    im = i-1
    if ( im .lt. N/2 ) then
       fftfreq = im
    else
       fftfreq = (im-N)
    end if
  end function fftfreq

#ifdef MPI
	!-------------------------------------------------------------
	!> Create fftw plans and slab decomposition for 
  !! forward and reverse transforms of 3D density 
  !! and/or velocity fields.
	!! 
	!! These are allocated early in the program and used as needed.
	!!
	!! @param rank MPI task number
	!! @param numtasks Number of MPI processes
	!! @param ng Number of grid points per dimension.
	!! @param ngL Number of grid points per dimension for grids 
	!! with only large-scale modes.
	!<------------------------------------------------------------	
  subroutine init_fftw(rank, numtasks, ng, ngL)
    integer(i4b), intent(in) :: rank, numtasks,ng
    integer(i4b), intent(in), optional :: ngL
    ! ----- Create plans
    call rfftw3d_f77_mpi_create_plan(planf, mpi_comm_world, ng, ng, ng,&
         fftw_real_to_complex, fftw_estimate)
    call rfftw3d_f77_mpi_create_plan(planr, mpi_comm_world, ng, ng, ng,&
         fftw_complex_to_real, fftw_estimate)
    ! c2c plan for drawing Gaussian fields
    if (.not. present(ngL)) then
       call fftw3d_f77_mpi_create_plan(planc2c, mpi_comm_world, ng, ng, ng,&
            fftw_backward, fftw_estimate)
    end if
    ! ----- Slab decomposition
    call rfftwnd_f77_mpi_local_sizes(planf, local_nlast, local_last_start,&
         local_nlast2_after_transform, local_last2_start_after_transform,&
         total_local_size)
    if (local_nlast /= ng/numtasks) then
       print *, rank,"  ERROR-deriv: fftw partition of f does not match my partition",local_nlast,ng/numtasks
    end if
    if (rank .eq. 0) then
       print *, "Initialized FFTW plans,  local_nlast: ", &
            & local_nlast,local_nlast2_after_transform
    end if
    ! ------------------------------------------------------------
    ! Initialize second set of plans for grid size ngL, if present
    ! ------------------------------------------------------------
    if (present(ngL)) then
       ! ----- Create plans
       call rfftw3d_f77_mpi_create_plan(planfL, mpi_comm_world, ngL, ngL, ngL,&
            fftw_real_to_complex, fftw_estimate)
       call rfftw3d_f77_mpi_create_plan(planrL, mpi_comm_world, ngL, ngL, ngL,&
            fftw_complex_to_real, fftw_estimate)
       ! c2c plan for drawing Gaussian fields
       call fftw3d_f77_mpi_create_plan(planc2cL, mpi_comm_world, ngL, ngL, ngL,&
            fftw_backward, fftw_estimate)
       ! ----- Slab decomposition
       call rfftwnd_f77_mpi_local_sizes(planfL, local_nlastL, local_last_startL,&
            local_nlast2_after_transformL, local_last2_start_after_transformL,&
            total_local_sizeL)
       if (local_nlastL /= ngL/numtasks) then
          print *, rank,"  ERROR-deriv: fftw partition of f does not match my partition",local_nlast,ng/numtasks
       end if
       if (rank .eq. 0) then
          print *, "Initialized FFTW plans,  local_nlast: ", &
               & local_nlastL,local_nlast2_after_transformL
       end if
    else
       planfL = planf
       planrL = planr
       planc2cL = planc2c
       local_nlastL = local_nlast
       local_last_startL = local_last_start
       local_nlast2_after_transformL = local_nlast2_after_transform
       local_last2_start_after_transformL = local_last2_start_after_transform
       total_local_sizeL = total_local_size
    end if
  end subroutine init_fftw

	!-------------------------------------------------------------
	!> Deallocate fftw plan memory
	!<------------------------------------------------------------	
  subroutine dealloc_fftw_plans()
    call rfftwnd_f77_mpi_destroy_plan(planf)
    call rfftwnd_f77_mpi_destroy_plan(planr)
    call fftwnd_f77_mpi_destroy_plan(planc2c)
  end subroutine dealloc_fftw_plans
#endif MPI


#ifndef FFTW
  ! ============================================================
  ! Wrappers for fftpack in the Sun performance library.
  !  These get used if code compiled without FFTW and without MPI flags.
  ! ============================================================
  subroutine fftpack_forward3d(iopt,N,x,y,ifac,trigs,work,scale)
    ! ------------------------------------------------------------
    ! r2c FFT:
    !  initialize if iopt = 0
    !  execute if    iopt = -1
    ! ------------------------------------------------------------
    use sunperf
    integer(i8b), intent(in)                   :: iopt
    integer(i4b), intent(in)                   :: N
    real(dp), dimension(N,N,N), intent(in)     :: x
    complex(cdp), dimension(N/2+1,N,N), intent(inout) :: y
    integer(i8b), dimension(:), intent(inout)  :: ifac
    real(dp), dimension(6*N), intent(inout)    :: trigs
    real(dp), dimension(:), intent(out)        :: work
    real(dp), intent(in)                       :: scale
    integer(i8b) :: N1,N2,N3,ldx1,ldx2,ldy1,ldy2,ierr

    N1 = int(N,i8b)
    N2 = int(N,i8b)
    N3 = int(N,i8b)
    ldx1 = N1
    ldx2 = N2
    ldy1 = N1/2+1
    ldy2 = N2
#ifdef NO64BITS
    call dfftz3(iopt,N1,N2,N3,scale,x,ldx1,ldx2,y,ldy1,ldy2, &
         & trigs,ifac,work,0,ierr)
#else
    call dfftz3_64(iopt,N1,N2,N3,scale,x,ldx1,ldx2,y,ldy1,ldy2, &
         & trigs,ifac,work,0,ierr)
#endif
    if (ierr .ne. 0) then
       print *, '##########'
       print *, "Error in dfftz3:",ierr
       select case(ierr)
       case(-1)
          print *, 'IOPT is not 0 or -1: ',iopt
       case(-2)
          print *, 'N1 < 0: ',N1
       case(-3)
          print *, 'N2 < 0: ',N2
       case(-4)
          print *, 'N3 < 0: ',N3
       case(-5)
          print *, '(LDX1 < N1) or (LDX not equal  2*LDY  when  X and Y are same array)'
       case(-6)
          print *, '(LDX2 < N2)'
       case(-7)
          print *, '(LDY1 < N1/2+1)'
       case(-8)
          print *, '(LDY2 < N2) or (LDY2 not equal  LDX2  when  X and Y are same array)'
       case(-9)
          print *, '(LWORK not equal 0) and (LWORK < (MAX(N,2*N2,2*N3) + 16*N3))'
       case(-10)
          print *, 'memory allocation failed'
       end select
       print *, '##########'
       stop
    end if
  end subroutine fftpack_forward3d

  subroutine fftpack_backward3d(iopt,N,x,y,ifac,trigs,work,scale)
    ! ------------------------------------------------------------
    ! c2r FFT:
    !  initialize if iopt = 0
    !  execute if    iopt = 1
    ! ------------------------------------------------------------
    use sunperf
    integer(i8b), intent(in)                   :: iopt
    integer(i4b), intent(in)                   :: N
    complex(cdp), dimension(N/2+1,N,N), intent(in) :: x
    real(dp), dimension(N+2,N,N), intent(inout)     :: y
    integer(i8b), dimension(:), intent(inout)  :: ifac
    real(dp), dimension(6*N), intent(inout)    :: trigs
    real(dp), dimension(:), intent(out)        :: work
    real(dp), intent(in)                       :: scale
    integer(i4b) :: N1,N2,N3,ldx1,ldx2,ldy1,ldy2,ierr

    N1 = int(N,i8b)
    N2 = int(N,i8b)
    N3 = int(N,i8b)
    ldx1 = N1/2+1
    ldx2 = N2
    ldy1 = 2*ldx1
    ldy2 = N2
#ifdef NO64BITS
    call zfftd3(iopt,N1,N2,N3,scale,x,ldx1,ldx2,y,ldy1,ldy2, &
         & trigs,ifac,work,0,ierr)
#else
    call zfftd3_64(iopt,N1,N2,N3,scale,x,ldx1,ldx2,y,ldy1,ldy2, &
         & trigs,ifac,work,0,ierr)
#endif
    if (ierr .ne. 0) then
       print *, '##########'
       print *, "Error in zfftd3:",ierr
       select case(ierr)
       case(-1)
          print *, 'IOPT is not 0 or 1'
       case(-2)
          print *, 'N1 < 0'
       case(-3)
          print *, 'N2 < 0'
       case(-4)
          print *, 'N3 < 0'
       case(-5)
          print *, '(LDX1 < N1/2+1)'
       case(-6)
          print *, '(LDX2 < N2)'
       case(-7)
          print *, 'LDY1 not equal 2*LDX1 when X and Y  are  same array'
       case(-8)
          print *, '(LDY1 < 2*LDX1) or (LDY1 is odd) when X and Y are not same array'
       case(-9)
          print *, '(LDY2 < N2) or (LDY2 not equal LDX2)  when  X and Y are same array'
       case(-10)
          print *, '(LWORK  not  equal  0)   and   ((LWORK   < MAX(N,2*N2,2*N3) + 16*N3)*NCPUS)'
       case(-11)
          print *, 'memory allocation failed'
       end select
       print *, '##########'
       stop
    end if
  end subroutine fftpack_backward3d

  subroutine fftpack_c2c3d(iopt,N,x,y,ifac,trigs,work,scale)
    ! ------------------------------------------------------------
    ! c2c FFT:
    !  initialize if   iopt = 0
    !  forward FFT if  iopt = -1
    !  backward FFT if iopt =  1
    ! ------------------------------------------------------------
    use sunperf
    integer(i8b), intent(in)                   :: iopt
    integer(i4b), intent(in)                   :: N
    complex(cdp), dimension(N,N,N), intent(in)    :: x
    complex(cdp), dimension(N,N,N), intent(inout) :: y
    integer(i8b), dimension(:), intent(inout)  :: ifac
    real(dp), dimension(6*N), intent(inout)    :: trigs
    real(dp), dimension(:), intent(out)        :: work
    real(dp), intent(in)                       :: scale
    integer(i8b) :: N1,N2,N3,ldx1,ldx2,ldy1,ldy2,ierr

    N1 = int(N,i8b)
    N2 = int(N,i8b)
    N3 = int(N,i8b)
    ldx1 = N1
    ldx2 = N2
    ldy1 = N1
    ldy2 = N2
#ifdef NO64BITS
    call zfftz3(iopt,N1,N2,N3,scale,x,ldx1,ldx2,y,ldy1,ldy2, &
         & trigs,ifac,work,0,ierr)
#else
    call zfftz3_64(iopt,N1,N2,N3,scale,x,ldx1,ldx2,y,ldy1,ldy2, &
         & trigs,ifac,work,0,ierr)
#endif
    if (ierr .ne. 0) then
       print *, '##########'
       print *, "Error in zfftz3:",ierr
       select case(ierr)
       case(-1)
          print *, 'IOPT is not 0, 1 or -1: ',iopt
       case(-2)
          print *, 'N1 < 0: ',N1
       case(-3)
          print *, 'N2 < 0: ',N2
       case(-4)
          print *, 'N3 < 0: ',N3
       case(-5)
          print *, '(LDX1 < N1)'
       case(-6)
          print *, '(LDX2 < N2)'
       case(-7)
          print *, '(LDY1 < N1 or (LDY1 not equal LDX1 when X and Y are same array)'
       case(-8)
          print *, '(LDY2 < N2) or (LDY2 not equal  LDX2  when  X and Y are same array)'
       case(-9)
          print *, '(LWORK not equal 0) and (LWORK < (MAX(N,2*N2,2*N3) + 16*N3)*2*NCPUS)'
       case(-10)
          print *, 'memory allocation failed'
       end select
       print *, '##########'
       stop
    end if
  end subroutine fftpack_c2c3d
#endif FFTW

!!$  ! ----- Definition of FFTW constants 
!!$  !       (copied from fftw_f77.i in the FFTW2.1.5 distribution)
!!$  !     This file contains PARAMETER statements for various constants
!!$  !     that can be passed to FFTW routines.  You should include
!!$  !     this file in any FORTRAN program that calls the fftw_f77
!!$  !     routines (either directly or with an #include statement
!!$  !     if you use the C preprocessor).
!!$
!!$  integer FFTW_FORWARD,FFTW_BACKWARD
!!$  parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)
!!$
!!$  integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
!!$  parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)
!!$
!!$  integer FFTW_ESTIMATE,FFTW_MEASURE
!!$  parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)
!!$
!!$  integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
!!$  parameter (FFTW_OUT_OF_PLACE=0)
!!$  parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)
!!$
!!$  integer FFTW_THREADSAFE
!!$  parameter (FFTW_THREADSAFE=128)
!!$
!!$  !     Constants for the MPI wrappers:
!!$  integer FFTW_TRANSPOSED_ORDER, FFTW_NORMAL_ORDER
!!$  integer FFTW_SCRAMBLED_INPUT, FFTW_SCRAMBLED_OUTPUT
!!$  parameter(FFTW_TRANSPOSED_ORDER=1, FFTW_NORMAL_ORDER=0)
!!$  parameter(FFTW_SCRAMBLED_INPUT=8192)
!!$  parameter(FFTW_SCRAMBLED_OUTPUT=16384)

  ! FFTW (only version 2 installed on cosma
!!$    call dfftw_init_threads
!!$    call dfftw_plan_with_nthreads(8)
!!$#ifdef 2DGRID
!!$    call dfftw_plan_dft_r2c_2d(plan,N,N,rho,rhoFT,FFTW_FORWARD,FFTW_ESTIMATE)
!!$#else
!!$    ! FFTW v 2.5.1
!!$    call fftw_f77_create_plan(plan,3,N,FFTW_FORWARD,FFTW_ESTIMATE)
!!$    ! FFTW v 3
!!$    call dfftw_plan_dft_r2c_3d(plan,N,N,N,rho,rhoFT,FFTW_FORWARD,FFTW_ESTIMATE)
!!$#endif
!!$    call dfftw_execute_dft_r2c(plan, rho, rhoFT)
!!$    ! ----- Erase the plan
!!$    call fftwnd_f77_destroy_plan(plan
!!$    call dfftw_destroy_plan(plan)

!!$    call dfftw_init_threads
!!$    call dfftw_plan_with_nthreads(8)
!!$#ifdef 2DGRID
!!$    call dfftw_plan_dft_c2r_2d(rplan,N,N,rhoFT,rho,FFTW_BACKWARD,FFTW_ESTIMATE)
!!$    call dfftw_plan_dft_r2c_2d(fplan,N,N,rho,resampFT,FFTW_FORWARD,FFTW_ESTIMATE)
!!$#else
!!$    call dfftw_plan_dft_c2r_3d(rplan,N,N,N,rhoFT,rho,FFTW_BACKWARD,FFTW_ESTIMATE)
!!$    call dfftw_plan_dft_r2c_3d(fplan,N,N,N,rho,resampFT,FFTW_FORWARD,FFTW_ESTIMATE)
!!$#endif
end module fft_nsample
