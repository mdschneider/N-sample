!-------------------------------------------------------------
!> Routines to resample the Fourier modes and apply the survey mask.
!<------------------------------------------------------------
module resample
  use types_nsample
  use io_nsample
  use parameters
  use fft_nsample, only: local_nlast, local_last_start, local_nlastL, local_last_startL
  use gridtools
  implicit none

  real(dp), private :: disp1, disp2, disp3, maxdisp1, maxdisp2, maxdisp3

contains

	!-------------------------------------------------------------
	!> Convert a gridded density field to perturbed scale factor 
  !! values based on matching the linear growth at perturbed 
  !! matter density.
  !!
  !! The scale factor is found by nearest neighbors interpolation 
  !! in the table of Omega_m(a), which is linearly spaced in Omega_m.
  !!
  !! The input density field is overwritten with the scale factor.
	!!
	!! @param delta Gridded density field, which is converted to gridded
	!! perturbed scale factor on output.
	!! @param ng Number of grid points per dimension in delta.
	!<------------------------------------------------------------
  subroutine convert_density_to_scale_factor(delta, ng)
    use cosmology
    use gridtools, only: print_gridarray
    integer, intent(in)                          :: ng
    real(dp), dimension(ng,ng,local_nlastL), intent(inout) :: delta
    real(dp)                                     :: ommin, omdelta, om
    integer :: i,j,k,ndx
    real :: t1,t2
    if (.not. atable_read) then
       print *, "ERROR: must call read_a_omega_table before convert_density_to_scale_factor"
       stop
    end if
    call cpu_time(t1)
    ommin = atable(1,1)
    omdelta = (atable(ntable,1) - ommin) / real(ntable-1,dp)
    delta = (((1.d0+delta)*Omega_m - ommin) / omdelta) + 1.d0
    call print_gridarray(delta, "omega index")
    do k = 1,local_nlastL
       do j = 1,ng
          do i = 1,ng
             om = delta(i,j,k)
             ndx = floor(om)
             if (ndx .lt. 1) then
                delta(i,j,k) = atable(1,2)
             else if (ndx .gt. ntable-1) then
                delta(i,j,k) = atable(ntable,2)
             else
                delta(i,j,k) = (om - ndx)*atable(ndx+1,2) + &
                     & (ndx+1-om)*atable(ndx,2)
!!$             delta(i,j,k) = atable(nint(delta(i,j,k)), 2)
             end if
          end do
       end do
    end do
    call cpu_time(t2)
    if (rank .eq. 0) print *, "Time for convert_density_to_scale_factor: ",t2-t1," seconds"
  end subroutine convert_density_to_scale_factor


	!-------------------------------------------------------------
	!> Resample modes with |k| < k_thresh and compute power spectra.
	!!
	!! This is the main routine to perform the mode resampling.
	!!
	!! @param ir_start Starting value for iteration of resamplings.
	!! @param nr Number of resamplings to perform.
	!! @param ng Number of grid points per dimension.
	!! @param ngL Number of grid points per dimension for grids with only
	!! large-scale Fourier modes.  This can be smaller than ng to 
	!! save memory.
	!! @param psnap Array with sorted particle positions from all the snapshots.
	!! @param outdir Directory for output power spectra.
	!! @param pk_theory Array with theoretical power spectrum for generating new 
	!! large-scale modes.
	!! @param iseed Seed for random number generator.
	!<------------------------------------------------------------	
  subroutine resample_modes(ir_start, nr, ng, ngL, psnap, outdir, pk_theory, iseed)
    use fft_nsample
    use gridtools, only: print_gridarray, psestimator, kurtosisEstimatorFromSmoothDensArray
    implicit none
    integer(i4b), intent(in)                                  :: ir_start, nr, ng, ngL
    real(dp), dimension(3, nsnap, npart/numtasks), intent(in) :: psnap
    character(len=*), intent(in)                              :: outdir
    real(dp), dimension(npk,2), intent(in)                    :: pk_theory
    integer(i4b), intent(inout)                               :: iseed
    integer(i4b)                             :: ir, npkout, rndx, stat
    real(dp), dimension(:,:,:), pointer      :: deltaL1, deltaL2
    real(dp), dimension(:), allocatable      :: pk
    character(len=10)                        :: filetag
    character(len=len_filename)              :: outfile
    real                                     :: t1,t2
    integer                                  :: seed(1)
    real(dp), dimension(:), allocatable      :: rnarray

    npkout = npkbins(ng)
    allocate(pk(0:npkout))
    call read_a_omega_table()

    seed(1) = iseed
    call random_seed(put=seed)

    ! Distributed kvecs: initialize global arrays
    call num_kvecs_per_pixel(ng)
    !
    do ir = ir_start,nr/2
       if (rank .eq. 0) then
          print *, "----------------------------------------"
          print *, "ir:",ir," / ",nr/2
       end if
       deallocate(deltaL1, stat=stat)
       nullify(deltaL1)
       deallocate(deltaL2, stat=stat)
       nullify(deltaL2)
       allocate(deltaL1(ngL+fftpad,ngL,local_nlastL),&
            & deltaL2(ngL+fftpad,ngL,local_nlastL))
#ifdef MPI
       call draw_modes(deltaL1, deltaL2, ngL, pk_theory)
#else ! MPI
       call draw_modes(deltaL1, deltaL2, ngL, pk_theory, iseed)
#endif MPI
       call print_gridarray(deltaL1(1:ngL,:,:), "deltaL1")
       call print_gridarray(deltaL2(1:ngL,:,:), "deltaL2")
       ! --- 1
       call gen_perturbed_density(deltaL1, psnap, ng, ngL)
       write(filetag,*) 2*ir-1
       call psfilename(outfile, filetag, outdir, ng)
       print *, rank,"  Writing power spectrum to file",trim(outfile)
       call power_spectrum_estimator(ng, pk, deltaL1, ir, nr, outfile)
       ! --- 2
       call gen_perturbed_density(deltaL2, psnap, ng, ngL)
       write(filetag,*) 2*ir
       call psfilename(outfile, filetag, outdir, ng)
       print *, rank,"  Writing power spectrum to file",trim(outfile)
       call power_spectrum_estimator(ng, pk, deltaL2, ir, nr, outfile)

       call print_gridarray(deltaL1(1:ng,:,:),"delta1")
       deallocate(deltaL1, stat=stat)
       if (stat /= 0) then 
          print *, rank,"  Unsuccessful deltaL1 deallocation:",stat, shape(deltaL1)
       else
          nullify(deltaL1)
       end if

       call print_gridarray(deltaL2(1:ng,:,:),"delta2")
       deallocate(deltaL2, stat=stat)
       if (stat /= 0) then
          print *, rank,"  Unsuccessful deltaL2 deallocation:",stat, shape(deltaL2)
       else
          nullify(deltaL2)
       end if
    end do

    deallocate(pk)
    if (allocated(numkpixels)) deallocate(numkpixels)
    if (allocated(kvecpixnum)) deallocate(kvecpixnum)
    if (allocated(nside)) deallocate(nside)
    if (allocated(rnarray)) deallocate(rnarray)
  end subroutine resample_modes

	!-------------------------------------------------------------
	!> From a large-scale mode realization, perturb particles and 
  !! grid to get the density with large-scale modes added.
  !! 
  !! @param delta On input this is the large-scale density (on ngL grid), 
  !!   on output it is the density of the perturbed density (on ng grid).
	!! @param psnap Array with sorted particle positions from all the snapshots.
	!! @param ng Number of grid points per dimension.
	!! @param ngL Number of grid points per dimension for grids with only
	!! large-scale Fourier modes.  This can be smaller than ng to 
	!! save memory.
	!<------------------------------------------------------------	
  subroutine gen_perturbed_density(delta, psnap, ng, ngL)
    use fft_nsample
    use gridtools, only: print_gridarray, part2grid
    use zeldovich
    real(dp), dimension(:,:,:), pointer :: delta
    !real(dp), dimension(:,:,:), intent(inout), allocatable       :: delta
    real(dp), dimension(3, nsnap, npart/numtasks), intent(in)    :: psnap
    integer(i4b), intent(in)                                     :: ng, ngL
    !
    real(dp), dimension(3, npart/numtasks)          :: coor
    integer(i4b), dimension(3)      :: ndx
    ! These arrays are held locally on each processor and use most of the memory
    !real(dp), dimension(3,ngL,ngL,ngL) :: disp
    real(dp), dimension(:,:,:,:), allocatable :: disp
    real(dp), dimension(ngL,ngL,ngL)   :: apert
    !
    real(dp)                        :: xscal, apertval
    real(dp)                        :: alow, ahigh, dlna, weight1, weight2
    integer(i4b)                    :: i, ilow, ihigh, stat
    integer(i4b), dimension(0:numtasks-1) :: apertndx, apertcounts
    real :: t1,t2,t3,t4
    integer(i4b), dimension(1) :: indx

    allocate(disp(3,ngL,ngL,ngL))

    if (size(delta(:,1,1)) /= ngL+fftpad .or. size(delta(1,:,1)) /= ngL &
         & .or. size(delta(1,1,:)) /= local_nlastL) then
       print *, rank,"  ERROR: bad dimensions for delta in gen_perturbed_density,",shape(delta)
    end if

    ! ----- Compute the Zeldovich displacements from deltaL
#ifdef ZELDOVICH
    call zeldovich_displacements(ngL, npart, delta, &
         & disp(:,:,:,(1+local_last_startL):(local_nlastL+local_last_startL)))
#endif ZELDOVICH
#ifdef MPI
    call mpi_allgather(mpi_in_place, 0, mpi_datatype_null, &
         & disp, 3*ngL**2*local_nlastL, mpi_double_precision, &
         & mpi_comm_world, ierr)
#endif MPI


    ! ----- Convert deltaL to perturbed scale factor values
    call convert_density_to_scale_factor(delta(1:ngL,:,:), ngL)
    call print_gridarray(delta(1:ngL,:,:), "pert. scale factor")
#ifdef MPI
    call mpi_allgather(delta(1:ngL,:,:), size(delta(1:ngL,:,:)), &
         & mpi_double_precision, &
         & apert, size(delta(1:ngL,:,:)), mpi_double_precision, & 
         & mpi_comm_world, ierr)
#else
    apert = delta(1:ngL,:,:)
#endif
    deallocate(delta, stat=stat)
    if (stat /= 0) print *, rank,"  ERROR deallocating delta in gen_perturbed_density"
    nullify(delta)

    ! ----- Compute perturbed particle positions
    xscal = boxsize/real(ngL)
    ipart: do i = 1,size(psnap(1,1,:))
#ifdef ZELDOVICH
       ! --- Add Zeldovich displacements to coordinates at a*
       coor(:,i) = psnap(:, isnap_a1+1, i)
       call sum_TSC(coor(:,i), disp, xscal, ngL)
#endif ZELDOVICH
       ! --- Find the perturbed scale factor, a', at this position
       apertval = interpgrid_TSC(coor(:,i), apert, xscal, ngL)
       ! --- Interpolate between snapshots to find coordinates at a'
       if (apertval .lt. snapshot_times(1)) then
          coor(:,i) = psnap(:,1,i)
!!$          print *, rank,"  *** a < amin:",i,apertval,snapshot_times(1)
       else if (apertval .gt. snapshot_times(nsnap)) then
          coor(:,i) = psnap(:,nsnap,i)
!!$          print *, rank,"  *** a > amax:",i,apertval,snapshot_times(nsnap)
       else
          indx = minloc(abs(apertval/snapshot_times - 1.d0)) ! index to nearest value in snapshot_times
          if (apertval .ge. snapshot_times(indx(1))) then
             ilow = indx(1)
          else
             ilow = indx(1) - 1
          end if
          ihigh = ilow + 1
          alow  = snapshot_times(ilow)
          ahigh = snapshot_times(ihigh)
          dlna = log(ahigh/alow)
          weight2 = log(apertval/alow) / dlna
          weight1 = log(ahigh/apertval) / dlna
          call interpolate_coordinates(coor(:,i), psnap(:,ilow,i), &
               & psnap(:,ihigh,i), weight1, weight2)
       end if
#ifdef ZELDOVICH
       ! --- Add Zeldovich displacement to interpolated particle position
       call sum_TSC(coor(:,i), disp, xscal, ngL)
#endif ZELDOVICH
    end do ipart


#ifdef DOTPLOT
    ! ----- Save particles for making dotplots
    call save_dotplot_particles(dotplotdir, "perturbed_positions.txt", &
         & 10.d0, coor, 1)
#endif DOTPLOT

    ! ----- Add particle weight to gridded density
    allocate(delta(ng+fftpad,ng,local_nlast))
    call part2grid(ng, npart/numtasks, boxsize, coor, delta)
    delta = delta*real(ng**3)/real(npart)

    deallocate(disp)
  end subroutine gen_perturbed_density

	!-------------------------------------------------------------
	!> Evaluate the model power spectrum at the given k index.
	!! @param kmag Wavenumber where the power spectrum should be 
	!! evaluated (h/Mpc).
	!! @param pk Array of log(k), log(Pk). Array dim Npk is set as 
	!! a parameter in Module parameters (sorry!).
	!<------------------------------------------------------------
  function PowerSpec(kmag, pk)
    real(dp), intent(in)                   :: kmag
    real(dp), dimension(Npk,2), intent(in) :: pk
    real(dp)                               :: PowerSpec

    real(dp)                               :: k0
    integer(i4b)                           :: k1, k2
#ifdef FITPACK
#ifdef PK_LOGK
    k0 = (log(kmag*kF) - pk(1,1))/(pk(2,1) - pk(1,1))
#else
    k0 = (kmag*kF-pk(1,1))/(pk(2,1) - pk(1,1))
#endif
    k1 = floor(k0)
    k2 = ceiling(k0)
    if (k1 .eq. k2) then
       PowerSpec = exp(pk(k1+1,2))/V
    else 
       PowerSpec = exp((k0-k1)*pk(k2+1,2) + (k2-k0)*pk(k1+1,2)) / V
    end if
#else ! FITPACK
    PowerSpec = exp( curv2(real(log(kmag*kF),sp), Npk, real(pk(:,1),sp), real(pk(:,2),sp), &
         & pkderivs, sigma_pkspl) ) / V
#endif FITPACK
  end function PowerSpec

  subroutine precompute_random_numbers(rnarray, rndx, ng, nr, iseed)
    ! ------------------------------------------------------------
    ! Precompute all the random numbers that will be needed for 
    ! the Gaussian mode draws on all processors
    ! ------------------------------------------------------------
    use fft_nsample
    use gridtools, only: draw_random_array
    real(dp), dimension(:), allocatable, intent(out) :: rnarray
    integer(i4b), intent(out)                        :: rndx
    integer(i4b), intent(in)                         :: ng, nr, iseed
    integer(i4b) :: count, i,j,k, maxcount
    real(dp) :: ki
    count = 0
    kloop: do k = 1,local_nlastL
       jloop: do j = 1,ng
          iloop: do i = 1,ng
             if (i .eq. 1 .and. j .eq. 1 .and. k .eq. 1) exit iloop
             ki = sqrt(real(fftfreq(i,ng)**2 + fftfreq(j,ng)**2 + &
                  & fftfreq(k+local_last_startL,ng)**2, dp))
             threshold: if (ki .le. real(nthresh,dp)) then
                count = count+1
             end if threshold
          end do iloop
       end do jloop
    end do kloop
#ifdef MPI
    call mpi_allreduce(count, maxcount, 1, mpi_integer, mpi_max, mpi_comm_world, ierr)
#else
    maxcount = count
#endif MPI
    allocate(rnarray(2*maxcount*nr))
    print *, rank,"  Size of random number array:",size(rnarray)
    call draw_random_array(rnarray, size(rnarray), iseed)
    rndx = 1
  end subroutine precompute_random_numbers


#ifdef MPI
	!-------------------------------------------------------------
	!> Generate a realization of the Fourier modes of a Gaussian
  !! density field by taking the real or imaginary part of a 
  !! complex Gaussian density field with twice the variance.
	!! (2 density fields are returned from the two complex components.)
	!! @param delta1 3D grid of Gaussian density.
	!! @param delta2 3D grid of Gaussian density.
	!! @param N Number of grid points per dimension
	!! @param pk Array with theoretical power spectrum for generating density fields.
	!! @param rndx Optional: index into pre-computed random number array, rnarray.
	!! @param rnarray Optional array of precomputed random numbers (not tested).
	!<------------------------------------------------------------	
  subroutine draw_modes(delta1, delta2, N, pk, rndx, rnarray)
    use ran_tools
    use fft_nsample

    real(dp), dimension(N+2,N,local_nlastL), intent(out) :: delta1, delta2
    integer(i4b), intent(in)                            :: N
    real(dp), dimension(Npk,2), intent(in)              :: pk
    integer(i4b), intent(inout), optional               :: rndx
    real(dp), dimension(:), intent(in), optional        :: rnarray
    !
    complex(cdp), dimension(N,N,local_nlastL)      :: delta, work
    integer(i4b)                                  :: i,j,k
    real(dp)                                      :: phiR,phiI,sigma,fac,kmag
    real(dp)                                      :: ki, r1,r2,r,phase
!!$    integer :: seed(1)
#ifdef CIC
    real(dp), dimension(N) :: mult
#endif
    real :: t1,t2

    call cpu_time(t1)

#ifdef VERBOSE
    print *, rank,"  Drawing new modes"
#endif

    delta = cmplx(0.d0, 0.d0)
    ! draw Fourier modes
    kloop: do k = 1,local_nlastL
       jloop: do j = 1,N
          iloop: do i = 1,N
             if (i .eq. 1 .and. j .eq. 1 .and. k .eq. 1) exit iloop
             ki = sqrt(real(fftfreq(i,N)**2 + fftfreq(j,N)**2 + &
                  & fftfreq(k+local_last_startL,N)**2, dp))
             threshold: if (ki .le. real(nthresh,dp)) then
                if (present(rnarray) .and. present(rndx)) then
                   r1 = rnarray(rndx)
                   r2 = rnarray(rndx+1)
                   rndx = rndx + 2
                else
                   call random_number(r2)
                   call random_number(r1)
                end if
                phase = TWOPI * r2
                r = sqrt(PowerSpec(ki, pk) * (-2.d0*log(r1)))
                delta(i,j,k) = r * cmplx(cos(phase), sin(phase))
             end if threshold
          end do iloop
       end do jloop
    end do kloop
!!$#ifdef CIC
!!$    ! ----- 3D CIC deconvolution
!!$    mult(1) = 1.d0
!!$    do i=2,N
!!$       kmag = fftfreq(i,N)*TWOPI/real(N,dp)/2.d0
!!$       mult(i) = (sin(kmag)/kmag)**2
!!$    end do
!!$    do k=1,local_nlastL
!!$       do j=1,N
!!$          do i=1,N
!!$             fac = mult(i)*mult(j)*mult(k+local_last_startL)
!!$             delta(i,j,k) = fac * delta(i,j,k)
!!$          end do
!!$       end do
!!$    end do
!!$#endif
    call fftwnd_f77_mpi(planc2cL, 1, delta, work, 0, fftw_normal_order)
    delta1 = 0.d0; delta2 = 0.d0
    delta1(1:N,:,:) = real(delta)
    delta2(1:N,:,:) = aimag(delta)
#ifdef VERBOSE
    print *, rank,"  Finished drawing modes"
#endif
    call cpu_time(t2)
    if (rank .eq. 0) print *, "Time for draw_modes: ",t2-t1," seconds"
  end subroutine draw_modes  
#else
	!-------------------------------------------------------------
	!> Generate a realization of the Fourier modes of a Gaussian
  !! density field by taking the real or imaginary part of a 
  !! complex Gaussian density field with twice the variance.
	!! (2 density fields are returned from the two complex components.)
	!! @param delta1 3D grid of Gaussian density.
	!! @param delta2 3D grid of Gaussian density.
	!! @param N Number of grid points per dimension
	!! @param pk Array with theoretical power spectrum for generating density fields.
	!! @param iseed Seed for random number generator.
	!<------------------------------------------------------------
  subroutine draw_modes(delta1, delta2, N, pk, iseed)
    use ran_tools
    use fft_nsample

    integer(i4b), intent(in)                      :: N
    real(dp), dimension(N, N, N), intent(out)     :: delta1, delta2
    real(dp), dimension(Npk,2), intent(in)        :: pk
    integer(i4b), intent(inout), optional         :: iseed
    complex(cdp), dimension(N,N,N)                :: phiFT,phi
    integer(i4b)                                  :: i,j,k,k1,k2
!!$    real(dp)                                      :: phiR,phiI,sigma
    real(dp)                                      :: r1,r2,phase,r
    real(dp)                                      :: ki,k0,delta

#ifdef FFTW
    integer(i8b)                                  :: planr
#else
    real(dp), dimension(2*(N+N+N))                :: trigsf ! trig. weights
    real(dp), dimension(:), allocatable           :: workf  ! workspace for fft
    integer(i8b), dimension(128*3)                :: ifacf  ! factors of N
#endif

    ! initialize ffts
#ifdef FFTW
    call dfftw_plan_dft_3d(planr, N, N, N, phiFT, phi, FFTW_ESTIMATE)
#else
    call fftpack_c2c3d(int(0,i8b),N,phiFT,phi,ifacf,trigsf,workf,1.0d0)
#endif

    phiFT = cmplx(0.d0, 0.d0)

    ! draw Fourier modes
    do k = 1,N
       do j = 1,N
          iloop: do i = 1,N
             if (i .eq. 1 .and. j .eq. 1 .and. k .eq. 1) exit iloop
             ki = sqrt(real(fftfreq(i,N)**2+fftfreq(j,N)**2+fftfreq(k,N)**2))
             threshold: if (ki .le. real(nthresh,dp)) then
!!$                sigma = sqrt(PowerSpec(real(nint(ki),dp), pk))
!!$                phiR = randgauss_boxmuller(iseed)*sigma
!!$                phiI = randgauss_boxmuller(iseed)*sigma
!!$                phiFT(i,j,k) = cmplx(phiR,phiI,kind=cdp)
                call random_number(r1)
                call random_number(r2)
                phase = TWOPI * r2
                r = sqrt(PowerSpec(ki, pk) * (-2.d0*log(r1)))
                phiFT(i,j,k) = r * cmplx(cos(phase), sin(phase))
             end if threshold
          end do iloop
       end do
    end do
    ! TODO: Add CIC deconvolution here
#ifdef FFTW
    call dfftw_execute(planr)
#else
    call fftpack_c2c3d(int(1,i8b),N,phiFT,phi,ifacf,trigsf,workf,1.0d0)
#endif
    delta1 = real(phi) !/ TWOPI**1.5
    delta2 = aimag(phi) !/ TWOPI**1.5
  end subroutine draw_modes
#endif MPI


  ! ------------------------------------------------------------
  ! Number of bins and bin indices for power spectrum estimator
  ! ------------------------------------------------------------
  function npkbins(ng)
    integer(i4b), intent(in) :: ng
    integer(i4b)             :: npkbins
    select case (pkbintype)
    case (1)
       npkbins = floor((ng*kF/2)/kbinwidth)
    case default
       npkbins = ng/2
    end select
  end function npkbins

  pure function pkbinndx(k)
    real(dp), intent(in) :: k
    integer(i4b)         :: pkbinndx
    select case (pkbintype)
    case (1)
       if (k .gt. 0.0) then
         pkbinndx = floor(k*kF/kbinwidth) + 1
       else
         pkbinndx = 0
       end if
    case default
       pkbinndx = nint(k)
    end select
  end function pkbinndx

  function k_frombinndx(i)
    integer(i4b), intent(in) :: i
    real(dp) :: k_frombinndx
    select case (pkbintype)
    case (1)
       k_frombinndx = (i-1)*kbinwidth + 0.5*kbinwidth
    case default
       k_frombinndx = i*TWOPI/boxsize
    end select
  end function k_frombinndx


#ifdef MPI
	!-------------------------------------------------------------
	!> Compute the power spectrum estimator from rhoFT by averaging 
  !! over k shells (or partial shells if k > k_thresh).
	!! @param ng Number of grid points per dimension.
	!! @param pkf Output power spectrum estimate.
	!! @param rho Gridded density field.
	!! @param ir Unused.
	!! @param Nr Unused.
	!! @param outfile Name of file for saving the power spectrum estimate.
	!! @param in_local_nlast Size of 3rd dimension of MPI distributed Fourier grid.
	!! @param in_local_last_start First index of distributed 3rd dimension of Fourier grid.
	!! @param in_planf Preallocated plan for FFT using FFTW v.2
	!<------------------------------------------------------------
  subroutine power_spectrum_estimator(ng, pkf, rho, ir, Nr, outfile, &
       & in_local_nlast, in_local_last_start, in_planf)
    use fft_nsample
    use gridtools, only: deconvolve_CIC
    integer(i4b), intent(in)                                :: ng,ir,Nr
!!$    real(dp), dimension(ng+2,ng,local_nlast), intent(inout) :: rho
    real(dp), dimension(:,:,:), intent(inout)               :: rho
    real(dp), dimension(0:), intent(inout)                  :: pkf
    character(len=*), intent(in), optional                  :: outfile
    integer(i4b), intent(in), optional :: in_local_nlast, in_local_last_start
    integer(i8b), intent(in), optional :: in_planf
    integer(i4b) :: tmp_local_nlast, tmp_local_last_start
    integer(i8b) :: tmp_planf
    integer(i4b), dimension(:), allocatable     :: nkf
    ! ----- Temporary power spectrum arrays
    real(dp), dimension(:), allocatable         :: pk
    integer(i4b), dimension(:), allocatable     :: nk    
    real(dp), dimension(ng)                     :: mult
    integer(i4b)                                :: npkout
    ! -----
    integer(i4b)                                :: nrside, unit_out
    integer(i4b)                                :: ii,jj,kk,ki,unit_out,i
    real(dp), dimension(:,:,:), allocatable     :: work
    real(dp)                                    :: fac,k
    real :: t1,t2

    if (present(in_local_nlast)) then
       tmp_local_nlast = in_local_nlast
    else
       tmp_local_nlast = local_nlast
    end if
    if (present(in_local_last_start)) then
       tmp_local_last_start = in_local_last_start
    else
       tmp_local_last_start = local_last_start
    end if
    if (present(in_planf)) then
       tmp_planf = in_planf
    else
       tmp_planf = planf
    end if

    npkout = npkbins(ng)
    allocate(nkf(0:npkout), pk(0:npkout), nk(0:npkout))

    ! ----- Do FT
    allocate(work(ng+fftpad, ng, tmp_local_nlast))
    call rfftwnd_f77_mpi(tmp_planf, 1, rho, work, 0, fftw_normal_order)

#ifdef CIC
    call deconvolve_CIC(ng, rho, tmp_local_nlast, tmp_local_last_start)
#endif CIC

    ! ----- Average over k shells to get 1D power spectrum
    pk = 0.0_dp
    nk = 0
    do kk=1,tmp_local_nlast
       do jj=1,ng
          do ii=1,ng !/2+1
                i = abs(fftfreq(ii,ng)) + 1
                k = sqrt(real(fftfreq(ii,ng)**2 + fftfreq(jj,ng)**2 + &
                     & fftfreq(kk+tmp_local_last_start,ng)**2))
                ki = pkbinndx(k)
                pkindex: if (ki .le. npkout) then ! max pk array index
                   pk(ki) = pk(ki) + rho(2*i-1,jj,kk)**2 + rho(2*i,jj,kk)**2
                   nk(ki) = nk(ki) + 1
                end if pkindex
          end do
       end do
    end do

    ! ----- Normalize by the number of modes in each shell
    pkf = 0.d0
    nkf = 0
    call mpi_allreduce(pk, pkf, size(pkf), mpi_double_precision, mpi_sum, &
         & mpi_comm_world, ierr)
    call mpi_allreduce(nk, nkf, size(nkf), mpi_integer, mpi_sum, &
         & mpi_comm_world, ierr)

    pkf = pkf*V/real(ng,dp)**6
    do ii = 0,npkout
       if (nkf(ii) .gt. 0) then
          pkf(ii) = pkf(ii)/nkf(ii)
       else
          pkf(ii) = 0.d0
       end if
    end do

    ! ----- Save to file if a filename argument was supplied
    if (rank .eq. 0) then
       if (present(outfile)) then
!!$          print *, rank,"  Writing power spectrum to ",trim(outfile)
          call get_free_unit(unit_out)
          open(unit_out,file=trim(outfile),status="replace",action="write")
          do i=1,npkout
             write(unit_out,'(2ES18.8,I10)') k_frombinndx(i), pkf(i), nkf(i)
          end do
          close(unit_out)
       end if
    end if
    deallocate(nkf, pk, nk)
  end subroutine power_spectrum_estimator
#else ! MPI
	!-------------------------------------------------------------
	!> Compute the power spectrum estimator from rhoFT by averaging 
  !! over k shells (or partial shells if k > k_thresh).
	!! @param ng Number of grid points per dimension.
	!! @param pkf Output power spectrum estimate.
	!! @param rho Gridded density field.
	!! @param ir Unused.
	!! @param Nr Unused.
	!! @param outfile Name of file for saving the power spectrum estimate.
	!<------------------------------------------------------------
  subroutine power_spectrum_estimator(ng, pk, rho, ir, Nr, outfile)
    use fft_nsample !, only: fftfreq
    integer(i4b), intent(in)                          :: ng,ir,Nr
    real(dp), dimension(ng,ng,ng), intent(in)         :: rho
    real(dp), dimension(0:), intent(inout)            :: pk
    character(len=*), intent(in), optional            :: outfile

    integer(i4b), dimension(:), allocatable           :: nk
    complex(cdp), dimension(ng/2+1,ng,ng)             :: rhoFT
!!$    real(dp), dimension(ng/2+1,ng,ng)                 :: power
    integer(i4b)                                      :: nrside
    integer(i4b)                                      :: ii,jj,kk,ki,i
    real(dp), dimension(ng)                           :: mult
    real(dp)                                          :: k
    integer(i4b)                                      :: unit_out, npkout

    real(dp), dimension(6*ng)                         :: trigs
    real(dp), dimension(:), allocatable               :: work
    integer(i8b), dimension(128*3)                    :: ifac
#ifdef FFTW
    integer(i8b) :: plan
#endif

    npkout = npkbins(ng)
    allocate(nk(0:npkout))

    ! ----- Forward FFT
#ifdef FFTW
    call dfftw_plan_dft_r2c_3d(plan,ng,ng,ng,rho,rhoFT,FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan,rho,rhoFT)
    call dfftw_destroy_plan(plan)
#else
    call fftpack_forward3d(int(0,i8b),ng,rho,rhoFT,ifac,trigs,work,1.0d0) ! initialize
    call fftpack_forward3d(int(-1,i8b),ng,rho,rhoFT,ifac,trigs,work,1.0d0) ! execute
#endif

#ifdef CIC
    ! ----- 3D CIC deconvolution
    mult(1) = 1.d0
    do ii=2,ng
       k = fftfreq(ii,ng)*TWOPI/real(ng,dp)/2.d0
       mult(ii) = 1.d0/(sin(k)/k)**2
    end do
    forall(ii=1:ng/2+1, jj=1:ng, kk=1:ng) rhoFT(ii,jj,kk) = &
         & mult(ii)*mult(jj)*mult(kk)*rhoFT(ii,jj,kk)
#endif

    ! ----- Average over k shells to get 1D power spectrum
    pk = 0.0_dp
    nk = 0
    do kk=1,ng
       do jj=1,ng
          do ii=1,ng !/2+1
             i = abs(fftfreq(ii,ng)) + 1
             k = sqrt(real(fftfreq(ii,ng)**2 + &
                  & fftfreq(jj,ng)**2 + fftfreq(kk,ng)**2))
             ki = pkbinndx(k)
             if (ki .le. npkout) then
                pk(ki) = pk(ki) + real(rhoFT(i,jj,kk) * conjg(rhoFT(i,jj,kk)),dp)
                nk(ki) = nk(ki) + 1
             end if
          end do
       end do
    end do
    print *, "Finished pk sum"

    ! ----- Normalize by the number of modes in each shell
    pk = pk*V/real(ng,dp)**6
    do ii=0,npkout
       if (nk(ii) .gt. 0) then
          pk(ii) = pk(ii) / nk(ii)
       else
          pk(ii) = 0.0
       end if
    end do
    print *, "Finished pk normalization"

    ! ----- Save to file if a filename argument was supplied
    if (present(outfile)) then
       print *, "Writing power spectrum to ",trim(outfile)
       call get_free_unit(unit_out)
       open(unit_out,file=trim(outfile),status="replace",action="write")
       do i=1,npkout
          write(unit_out,'(2ES18.8,I10)') k_frombinndx(i), pk(i), nk(i)
       end do
       close(unit_out)
    end if
    deallocate(nk)
  end subroutine power_spectrum_estimator
#endif

	!-------------------------------------------------------------
	!> Compute particle positions at the scale factor times indicated
  !! in apert. (Not used for mode addition.)
  !!
  !! Cycle through all snapshots, selecting appropriate particle
  !! positions from each one.
	!<------------------------------------------------------------	
  subroutine interp_part_pos(ng, apert, coor, disp, add_power)
    use parameters
    use zeldovich
    integer(i4b), intent(in)                    :: ng
    real(dp), dimension(ng,ng,ng), intent(in)   :: apert
    real(dp), dimension(3,npart), intent(inout) :: coor
    real(dp), dimension(3,ng,ng,ng), intent(in) :: disp
    logical, intent(in)                         :: add_power

    real(dp), dimension(:), allocatable :: asnap
    real(dp)                      :: a1, a2, weight1,weight2
    real(dp)                      :: xscal
    integer(i4b), dimension(3)    :: ndx1, ndx2
    integer(i4b) :: isnap_start, isnap_stop, unit_in, isnap, i, count,np_updated
    real(dp) :: amin, amax, dlna, lna1, lna2, asnap1, asnap2
    real(dp), dimension(:, :), allocatable :: coor1, coor2
    integer(i8b), dimension(:), allocatable :: ids1,ids2
    integer(i4b), dimension(:), allocatable :: ptr
    logical, dimension(npart) :: partcount
    real(dp), dimension(npart) :: a_particle

    real(dp) :: w1mean, w2mean
    w1mean = 0.d0
    w2mean = 0.d0

    xscal = boxsize / real(ng,dp)

    ! ----- Find the perturbed times associated with each particle, 
    !       with the particle positions evaluated at a = 1
    allocate(ids1(npart), ids2(npart))
    allocate(coor1(3,npart), coor2(3,npart))
    call read_snapshot(smallmodesdir, isnap_a1, coor1, ids1, asnap1)
    call save_dotplot_particles(dotplotdir, "smallmodes_a1.txt", 10.d0, coor1, 1)
#ifdef ZELDOVICH
    if (add_power) then
       call add_particle_displacements(ng, npart, asnap1, boxsize, disp, coor1)
       call save_dotplot_particles(dotplotdir, "smallmodes_a1_zeldovich.txt", 10.d0, coor1, 1)
    end if
#endif
    do i = 1,npart
       ! Use nearest-neighbor interpolation to find the perturbed scale
       ! factor value at the particle position.  The perturbed scale 
       ! factor should only vary on large scales so the interpolation 
       ! method should not matter.
       ndx1 = floor(coor1(:,i)/xscal) + 1
       a_particle(i) = apert(ndx1(1), ndx1(2), ndx1(3))       
    end do


    ! ----- Read times where snapshots were saved
    allocate(asnap(nsnap))
    call get_free_unit(unit_in)
    open(unit_in, file=snapshotlist, &
         & status="old", action="read")
    read(unit_in,*) asnap
    close(unit_in)    


    ! ----- Which snapshots to loop over?
    amin = minval(a_particle)
    amax = maxval(a_particle)
    isnap_start = 1
    do isnap = 1,nsnap
       if (asnap(isnap) .ge. amin) exit
       isnap_start = isnap
    end do
    do isnap = isnap_start, nsnap-1
       isnap_stop = isnap
       if (asnap(isnap) .gt. amax) exit
    end do
    print *, "min/max snapshots:",&
         isnap_start,isnap_stop,asnap(isnap_start),asnap(isnap_stop),amin,amax
!!$    deallocate(asnap)


    ! ----- Linearly interpolate particle positions between snapshots
    call read_snapshot(smallmodesdir, isnap_start-1, coor1, ids1, asnap1)
#ifdef ZELDOVICH
    call add_particle_displacements(ng, npart, asnap1, boxsize, disp, coor1)
#endif

    disp1 = 0.d0; disp2 = 0.d0; disp3 = 0.d0
    maxdisp1 = 0.d0; maxdisp2 = 0.d0; maxdisp3 = 0.d0
    np_updated = 0

    coor = 0.d0
    partcount = .false.
    snap: do isnap = isnap_start, isnap_stop-1
       print *, ""
       print *, "--------------------------------------------------"
       call read_snapshot(smallmodesdir, isnap, coor2, ids2, asnap2)
#ifdef ZELDOVICH
       call add_particle_displacements(ng, npart, asnap2, boxsize, disp, coor2)
#endif
       dlna = log(asnap2/asnap1) !log(asnap(isnap+1)/asnap(isnap))
       lna1 = log(asnap1)    !log(asnap(isnap))
       lna2 = log(asnap2)    !log(asnap(isnap+1))
       print *, "asnap 1,2: ",asnap1,asnap2,lna1,lna2,dlna
       particles: do i = 1,npart
          if (ids1(i) /= ids2(i)) print *, &
               & "ERROR: ids do not match:",i,ids1(i),ids2(i)
          a1 = a_particle(i)

!!$          if (a1 .lt. amin .and. isnap .eq. isnap_start) then
          if (a1 .lt. asnap(isnap_start) .and. isnap .eq. isnap_start) then
             coor(:,i) = coor1(:,i)
             w1mean = w1mean + 1.d0
             np_updated = np_updated + 1
             partcount(i) = .true.
             print *, "*** a1 < amin:",a1,amin
!!$          else if (a1 .ge. amax .and. isnap .eq. isnap_stop-1) then
          else if (a1 .ge. asnap(isnap_stop) .and. isnap .eq. isnap_stop-1) then
             coor(:,i) = coor2(:,i)
             w2mean = w2mean + 1.d0
             np_updated = np_updated + 1
             partcount(i) = .true.
             print *, "*** a1 > amax:",a1,amax
          else if (a1 .ge. asnap1 .and. a1 .lt. asnap2) then
             weight2 = (log(a1) - lna1) / dlna
             weight1 = (lna2 - log(a1)) / dlna
             if (weight1 .gt. 1.d0 .or. weight1 .lt. 0.d0) print *, "ERROR weight1: ",weight1
             if (weight2 .gt. 1.d0 .or. weight2 .lt. 0.d0) print *, "ERROR weight2: ",weight2
             w1mean = w1mean + weight1
             w2mean = w2mean + weight2
             call interpolate_coordinates(coor(:,i), coor1(:,i), coor2(:,i), &
                  & weight1, weight2)
             call check_particle_interpolation(coor(:,i),coor1(:,i),coor2(:,i))
             np_updated = np_updated + 1
             partcount(i) = .true.
          end if
       end do particles

       ! ----- Error checking and diagnostics
       print *, "Number of particles not interpolated: ", npart - count(partcount)
       print *, "  mean weights: ",w1mean/real(np_updated),w2mean/real(np_updated)
       print *, "  mean displacements: ",&
            & disp1/real(np_updated),disp2/real(np_updated),disp3/real(np_updated)
       print *, "  max  displacements: ",maxdisp1,maxdisp2,maxdisp3
       if (real(disp1,sp) .gt. real(disp3,sp) .or. &
            & real(disp2,sp) .gt. real(disp3,sp)) &
            & print *, "ERROR: bad displacements"
       if (real(maxdisp1,sp) .gt. real(maxdisp3,sp) .or. &
            & real(maxdisp2,sp) .gt. real(maxdisp3,sp)) &
            & print *, "ERROR: bad max displacements"

       ! ----- Reset quantities for next snap loop iteration
       w1mean = 0.d0; w2mean = 0.d0
       disp1 = 0.d0; disp2 = 0.d0; disp3 = 0.d0;
       maxdisp1 = 0.d0; maxdisp2 = 0.d0; maxdisp3 = 0.d0
       np_updated = 0
       coor1 = coor2
       ids1  = ids2
       asnap1 = asnap2
    end do snap

  end subroutine interp_part_pos

	!-------------------------------------------------------------
	!> Linearly interpolate between coor1 and coor2 with weights
  !! w1, w2 while taking account of the periodic box.
	!! (Not used for mode addition.)
	!<------------------------------------------------------------	
  subroutine interpolate_coordinates(coor, coor1, coor2, w1, w2)
    real(dp), dimension(3), intent(inout) :: coor
    real(dp), dimension(3), intent(in)    :: coor1,coor2
    real(dp), intent(in)                  :: w1,w2
    integer(i4b) :: i
    do i=1,3
       if (abs(coor1(i) - coor2(i)) .le. boxsize/2.d0) then
          coor(i) = w1*coor1(i) + w2*coor2(i)
       else if (coor1(i) .gt. boxsize/2.d0) then
          coor(i) = w1 * (coor1(i)-boxsize) + w2*coor2(i)
       else if (coor2(i) .gt. boxsize/2.d0) then
          coor(i) = w1*coor1(i) + w2 * (coor2(i) - boxsize)
       else
          print *, "ERROR: strange snapshot coordinates: ",coor1(i),coor2(i)
       end if
       ! make sure output coordinates are in the range [0,boxsize)
       if (coor(i) .lt. 0.d0 .and. coor(i) .ge. -boxsize) then
          coor(i) = coor(i) + boxsize
       else if (coor(i) .ge. boxsize .and. coor(i) .lt. 2*boxsize) then
          coor(i) = coor(i) - boxsize
       end if
    end do
  end subroutine interpolate_coordinates

	!-------------------------------------------------------------
	!> Error checking for interpolated particle coordinates.
	!! (Not used for mode addition.)
  !!
  !!  *disp* variables are defined globally within the module
	!<------------------------------------------------------------	
  subroutine check_particle_interpolation(coor,coor1,coor2)
    ! ------------------------------------------------------------
    ! 
    ! ------------------------------------------------------------
    real(dp), dimension(3), intent(in) :: coor,coor1,coor2
    real(dp), dimension(3) :: cmin,cmax
    real(dp) :: r1,r2,r3
    integer(i4b) :: i
    forall (i=1:3) cmin(i) = min(coor1(i), coor2(i))
    forall (i=1:3) cmax(i) = max(coor1(i), coor2(i))
!!$    do i=1,3
!!$       if (real(coor(i),sp) .gt. real(cmax(i),sp) .or. real(coor(i),sp) .lt. real(cmin(i),sp)) then
!!$          print *, "  ERROR: bad interpolation ",i,cmin(i),coor(i),cmax(i)
!!$       end if
!!$    end do
    r1 = particle_separation(coor,  coor1)
    r2 = particle_separation(coor,  coor2)
    r3 = particle_separation(coor1, coor2)
    if (r1 .gt. boxsize) r1 = r1 - boxsize
    if (r2 .gt. boxsize) r2 = r2 - boxsize
    if (r3 .gt. boxsize) r3 = r3 - boxsize
    disp1 = disp1 + r1
    disp2 = disp1 + r2
    disp3 = disp3 + r3
    maxdisp1 = max(maxdisp1, r1)
    maxdisp2 = max(maxdisp2, r2)
    maxdisp3 = max(maxdisp3, r3)
  end subroutine check_particle_interpolation


  function particle_separation(coor1, coor2)
    real(dp), dimension(3), intent(in) :: coor1, coor2
    real(dp) :: particle_separation, delta(3)
    integer(i4b) :: i
    forall (i=1:3) delta(i) = min(abs(coor1(i) - coor2(i)), &
         & abs(abs(coor1(i) - coor2(i)) - boxsize))
    particle_separation = sqrt(dot_product(delta, delta))
  end function particle_separation


  subroutine add_linear_time_evolution(coor, ng, da)
    ! ------------------------------------------------------------
    ! Add the linear shift in particle coordinates based on the 
    ! particle velocities
    ! ------------------------------------------------------------
    real(dp), dimension(6,npart), intent(inout) :: coor
    integer(i4b), intent(in)                    :: ng
    real(dp), dimension(ng,ng,ng), intent(in)   :: da

    real(dp)                                    :: f,xscal,g
    integer(i4b)                                :: i,ii,jj,kk,l,m,n
    integer(i4b), dimension(3,npart)            :: ndx
    real(dp), dimension(3,npart)                :: dist
    real(dp), dimension(npart)                  :: vol

    xscal = boxsize/real(ng)
    ndx=int(coor(1:3,:)/xscal)+1 ! Nearest grid point index to coordinate position
    dist=real(ndx,dp)-(coor(1:3,:)/xscal) ! 1 - distance from each coordinate component
    do ii=0,1
       do jj=0,1
          do kk=0,1
             vol=abs(real(ii,dp)-dist(1,:))*&
                  & abs(real(jj,dp)-dist(2,:))*&
                  & abs(real(kk,dp)-dist(3,:))
             do i = 1,npart
                l = modulo(ndx(1,i)+ii-1,ng)+1
                m = modulo(ndx(2,i)+jj-1,ng)+1
                n = modulo(ndx(3,i)+kk-1,ng)+1
                coor(1:3,i) = coor(1:3,i) + (da(l,m,n)*vol(i)*coor(4:6,i))
             end do
          end do
       end do
    end do
    where(coor .lt. 0.0) coor = coor + boxsize
    where(coor .ge. boxsize) coor = coor - boxsize
  end subroutine add_linear_time_evolution

end module resample
