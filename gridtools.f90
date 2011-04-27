! ------------------------------------------------------------
!> Routines for moving particles onto and off of grids
! ------------------------------------------------------------
module gridtools
  use types_nsample
  use parameters
  implicit none

  integer(i4b), dimension(8) :: npfile
  integer(i4b), dimension(2) :: npproc
  integer(i4b)               :: maxnp



contains

#ifdef MPI
	!-------------------------------------------------------------
	!> Print ranges of cubic grid array for debugging.
	!! @param a 3D array to print
	!! @param name Descriptive name for print output
	!<------------------------------------------------------------	
  subroutine print_gridarray(a,name)
    implicit none
    real(dp), dimension(:,:,:), intent(in) :: a
    character(len=*), intent(in) :: name
    real(dp) :: minv,maxv,meanv,varv
    call mpi_allreduce(minval(a), minv, 1, mpi_double_precision, &
         & mpi_min, mpi_comm_world, ierr)
    call mpi_allreduce(maxval(a), maxv, 1, mpi_double_precision, &
         & mpi_max, mpi_comm_world, ierr)
    call mpi_allreduce(sum(a), meanv, 1, mpi_double_precision, &
         & mpi_sum, mpi_comm_world, ierr)
    meanv = meanv/real(numtasks*size(a))
    call mpi_allreduce(sum((a-meanv)**2), varv, 1, mpi_double_precision, &
         & mpi_sum, mpi_comm_world, ierr)
    varv = varv/real(numtasks*size(a))
    if (rank .eq. 0) then
       print *, " "
       print *, "Total "//trim(name)//" ranges: ",minv,meanv,maxv,varv
    end if
  end subroutine print_gridarray
#else
	!-------------------------------------------------------------
	!> Print ranges of cubic grid array for debugging.
	!! @param a 3D array to print
	!! @param name Descriptive name for print output
	!<------------------------------------------------------------	
  subroutine print_gridarray(a,name)
    real(dp), dimension(:,:,:), intent(in) :: a
    character(len=*), intent(in)           :: name
    real(dp) :: meana
    print *, "Range of ",trim(name),":"
    meana = sum(a)/size(a)
    print '(4ES18.8)', minval(a), meana, maxval(a), sum((a-meana)**2)/size(a)
  end subroutine print_gridarray

  subroutine hist_grid_values(a, ng, binwidth, binmin, bincounts, nbins)
    ! ------------------------------------------------------------
    ! Sum the number of grid values in bins of width binwidth
    ! ------------------------------------------------------------
    integer(i4b), intent(in)                  :: ng, nbins
    real(dp), dimension(ng,ng,ng), intent(in) :: a
    real(dp), intent(in)                      :: binwidth, binmin
    integer(i8b), dimension(nbins), intent(inout) :: bincounts
    integer(i4b) :: i,j,k, ndx

    bincounts = 0_i8b
    do i=1,ng
       do j=1,ng
          do k=1,ng
             ndx = int((a(i,j,k) - binmin) / binwidth, i4b)
             if (ndx .le. nbins .and. ndx .gt. 0) bincounts(ndx) = bincounts(ndx) + 1_i8b
          end do
       end do
    end do
  end subroutine hist_grid_values
#endif MPI


#ifdef MPI
	!-------------------------------------------------------------
	!> Compute number density on a grid using CIC interpolation.
  !!
  !! The number of particles to read from coor on each processor 
  !! are stored in np_proc(numtasks).
  !!
	!! @param ng Number of grid points per dimension.
	!! @param np Unused.
	!! @param Box Simulation box size per dimension.
	!! @param coor Array with particle coordinates.
	!! @param rho Output array with gridded mass density field.
	!! @param vx Optional particle velocities in x-direction.
	!! @param vy Optional particle velocities in y-direction.
	!! @param vz Optional particle velocities in z-direction.		
	!<------------------------------------------------------------	
  subroutine part2grid(ng, np, Box, coor, rho, vx, vy, vz)
    integer(i4b), intent(in)                :: ng, np
    real(dp), intent(in)                    :: Box
    real(dp), dimension(:,:), intent(inout) :: coor
    real(dp), dimension(:,:,:), intent(out) :: rho
    real(dp), dimension(:,:,:), intent(out), optional :: vx,vy,vz

    real(dp)     :: xscal,vol,scal
    integer(i4b) :: ii,jj,kk,i,j,l,m,n,npp,ierr,tasknum,ng_per_proc, tag
    integer(i4b), dimension(3) :: ndx
    real(dp), dimension(3)     :: dist
    real(dp), dimension(:,:), allocatable :: tmp
    integer(i4b), dimension(2) :: reqs
    integer(i4b), dimension(mpi_status_size,2) :: stats

    ! --- FIXME: this is assumes a particular distribution of particles on processors
    integer(i4b), dimension(0:numtasks-1) :: np_proc
    
    np_proc = npart / numtasks
    ! --------------------

    xscal = Box/real(ng)
    rho = 0.d0

    if (present(vx) .and. present(vy) .and. present(vz)) then
       vx = 0.d0; vy = 0.d0; vz = 0.d0
       allocate(tmp(6,size(coor(1,:))))
    else
       allocate(tmp(3,size(coor(1,:))))
    end if
    if (size(tmp) /= size(coor)) then
       print *, "ERROR part2grid tmp/coor sizes do not match"
       call mpi_abort(ierr)
    end if

    tag = 1
    do j = 1,numtasks
       if (j .gt. 1) then ! cycle coor between processes
          call mpi_irecv(tmp, size(tmp), mpi_double_precision, &
               modulo(rank-1,numtasks), tag, mpi_comm_world, reqs(1), ierr)
          call mpi_isend(coor, size(coor), mpi_double_precision, &
               modulo(rank+1,numtasks), tag, mpi_comm_world, reqs(2), ierr)
          call mpi_waitall(2, reqs, stats, ierr)
          coor = tmp
       end if
       call mpi_barrier(mpi_comm_world, ierr)

       npp = np_proc(modulo(rank-(j-1),numtasks))
       ng_per_proc = ng/numtasks
       do ii=0,1
          do jj=0,1
             do kk=0,1
                do i=1,npp
                   ndx = int(coor(1:3,i)/xscal)+1
                   dist = real(ndx,dp) - (coor(1:3,i)/xscal)
                   vol = abs(real(ii)-dist(1))*abs(real(jj)-dist(2))*&
                        abs(real(kk)-dist(3))
                   l = modulo(ndx(1)+ii-1,ng)+1
                   m = modulo(ndx(2)+jj-1,ng) + 1
                   n = modulo(ndx(3)+kk-1,ng) + 1
                   tasknum = floor(real(n-1)/real(ng_per_proc))
                   n = n - (rank*ng_per_proc)
                   if (rank .eq. tasknum) then ! ignore particles with z coord outside range for this process
                      if (n .le. 0 .or. n > ng_per_proc) then
                         print *, rank,j,"  ERROR: n out of bounds: ",n,tasknum,ng_per_proc
                      end if
                      rho(l,m,n) = rho(l,m,n) + vol
                      if (present(vx) .and. present(vy) .and. present(vz)) then
                         vx(l,m,n) = vx(l,m,n) + vol*coor(4,i)
                         vy(l,m,n) = vy(l,m,n) + vol*coor(5,i)
                         vz(l,m,n) = vz(l,m,n) + vol*coor(6,i)
                      end if
                   end if
                end do  ! i = 1,npp
             end do     ! kk
          end do        ! jj
       end do           ! ii
    end do              ! j = 1,numtasks
    call mpi_barrier(mpi_comm_world, ierr)
    scal = real(ng**3)/real(sum(np_proc))
    rho = rho*scal
    if(present(vx) .and. present(vy) .and. present(vz)) then
       vx  = vx*scal
       vy  = vy*scal
       vz  = vz*scal
    end if
  end subroutine part2grid
#else !MPI
	!-------------------------------------------------------------
	!> Compute number density on a grid using CIC interpolation.
  !!
  !! The number of particles to read from coor on each processor 
  !! are stored in np_proc(numtasks).
  !!
	!! @param ng Number of grid points per dimension.
	!! @param np Unused.
	!! @param Box Simulation box size per dimension.
	!! @param coor Array with particle coordinates.
	!! @param rho Output array with gridded mass density field.
	!! @param vx Optional particle velocities in x-direction.
	!! @param vy Optional particle velocities in y-direction.
	!! @param vz Optional particle velocities in z-direction.		
	!<------------------------------------------------------------
  subroutine part2grid(ng,np,Box,coor,rho,vx,vy,vz)
    integer, intent(in)                        :: ng,np
    real(dp), dimension(ng,ng,ng), intent(inout) :: rho
    real(dp), dimension(ng,ng,ng), intent(out), optional :: vx,vy,vz
    real(dp), intent(in)                       :: Box
    real(dp), dimension(:,:), intent(inout)   :: coor

    real(dp)                                   :: xscal,scal
    integer                                    :: i,ii,jj,kk,l,m,n
    integer(i4b), dimension(3,size(coor(1,:))) :: ndx
    real(dp), dimension(3,size(coor(1,:)))     :: dist
    real(dp)                                   :: vol

    print *, ""
    print *, "---Assigning particles to grid"
    xscal = Box/real(ng)
    ndx=int(coor(1:3,:)/xscal, i4b) + 1 ! Nearest grid point index to coordinate position
    write(6,*) "Index ranges"
    write(6,*) minval(ndx(1,:)),minval(ndx(2,:)),minval(ndx(3,:))
    write(6,*) maxval(ndx(1,:)),maxval(ndx(2,:)),maxval(ndx(3,:))
    dist=ndx-(coor(1:3,:)/xscal) ! 1 - distance from each coordinate component
    rho=0.d0
    if (present(vx) .and. present(vy) .and. present(vz)) then
       vx = 0.d0; vy = 0.d0; vz = 0.d0
    end if
    do ii=0,1
       do jj=0,1
          do kk=0,1
             do i=1,size(coor(1,:)) !np
                vol=abs(real(ii,dp)-dist(1,i))*&
                     & abs(real(jj,dp)-dist(2,i))*&
                     & abs(real(kk,dp)-dist(3,i))
                l = modulo(ndx(1,i)+ii-1,ng)+1
                m = modulo(ndx(2,i)+jj-1,ng)+1
                n = modulo(ndx(3,i)+kk-1,ng)+1
                rho(l,m,n) = rho(l,m,n) + vol
                if(present(vx) .and. present(vy) .and. present(vz)) then
                   vx(l,m,n) = vx(l,m,n) + vol*coor(4,i)
                   vy(l,m,n) = vy(l,m,n) + vol*coor(5,i)
                   vz(l,m,n) = vz(l,m,n) + vol*coor(6,i)
                end if
             end do
          end do
       end do
    end do
    scal = real(ng**3)/real(np)
    rho = rho*scal
    if(present(vx) .and. present(vy) .and. present(vz)) then
       vx  = vx*scal
       vy  = vy*scal
       vz  = vz*scal
       where(rho>0.0) 
          vx=vx/rho
          vy=vy/rho
          vz=vz/rho
       end where
    end if
    print *, 'rho mean:',sum(rho)/ng**3
  end subroutine part2grid
#endif MPI

	!-------------------------------------------------------------
	!> Interpolate gridded values in disp to coor via CIC.
	!!
	!! Used to interplate Zel'dovich displacements.
	!!
	!! @param coor Coordinate point where interpolation is evaluated
	!! @param disp 4D array to be interpolated to coor
	!! @param xscal Scaling of coordinates to "box scale"
	!! @param ng Number of grid points per dimension.
	!<------------------------------------------------------------	
  subroutine sum_CIC(coor, disp, xscal, ng)
    real(dp), dimension(3), intent(inout)    :: coor
    real(dp), dimension(:,:,:,:), intent(in) :: disp
    real(dp), intent(in)                     :: xscal
    integer(i4b), intent(in)                 :: ng
    integer(i4b), dimension(3) :: ndx
    real(dp), dimension(3)     :: dist
    integer(i4b)               :: ii,jj,kk
    integer(i4b) :: l,m,n
    real(dp) :: vol
    ndx=floor(coor/xscal)+1 ! Nearest grid point index to coordinate position
    dist=real(ndx,dp)-(coor/xscal) ! 1 - distance from each coordinate component
    do ii=0,1
       do jj=0,1
          do kk=0,1
             vol=abs(real(ii,dp)-dist(1))*&
                  & abs(real(jj,dp)-dist(2))*&
                  & abs(real(kk,dp)-dist(3))
             l = modulo(ndx(1)+ii-1,ng)+1
             m = modulo(ndx(2)+jj-1,ng)+1
             n = modulo(ndx(3)+kk-1,ng)+1
             coor = coor + (disp(:,l,m,n)*vol)
          end do
       end do
    end do
    where(coor .lt. 0.0) coor = coor + boxsize
    where(coor .ge. boxsize) coor = coor - boxsize
  end subroutine sum_CIC

	!-------------------------------------------------------------
	!> Grid point weights for triangular shaped cloud interpolation
	!! @param ii Grid index
	!! @param dx Grid spacing
	!<------------------------------------------------------------	
  pure function tscweight(ii,dx)
    integer(i4b), intent(in) :: ii
    real(dp), intent(in)     :: dx
    real(dp) :: tscweight
       select case (ii)
       case (-1)
          tscweight = 0.5*(1.5 - (1.0 + dx))**2
       case (0)
          tscweight = 0.75 - dx**2
       case (1)
          tscweight = 0.5*(1.5 - (1.0 - dx))**2
       end select
  end function tscweight
 
	!-------------------------------------------------------------
	!> Interpolate gridded values in disp to coor via TSC.
	!!
	!! Used to interpolate Zel'dovich displacements.
	!!
	!! @param coor Coordinate point where interpolation is evaluated
	!! @param disp 4D array to be interpolated to coor
	!! @param xscal Scaling of coordinates to "box scale"
	!! @param ng Number of grid points per dimension.
	!<------------------------------------------------------------	
  subroutine sum_TSC(coor, disp, xscal, ng)
    real(dp), dimension(3), intent(inout)    :: coor
    real(dp), dimension(:,:,:,:), intent(in) :: disp
    real(dp), intent(in)                     :: xscal
    integer(i4b), intent(in)                 :: ng
    integer(i4b), dimension(3) :: ndx
    real(dp), dimension(3)     :: dist
    integer(i4b)               :: ii,jj,kk
    integer(i4b) :: l,m,n
    real(dp) :: wx,wy,wz
    ndx=floor(coor/xscal)+1 ! Nearest grid point index to coordinate position
    dist=(coor/xscal) - real(ndx-1,dp) ! distance from each coordinate component
    do ii=-1,1
       wx = tscweight(ii, dist(1))
       do jj=-1,1
          wy = tscweight(jj, dist(2))
          do kk=-1,1
             wz = tscweight(kk, dist(3))
             l = modulo(ndx(1)+ii-1+ng,ng)+1
             m = modulo(ndx(2)+jj-1+ng,ng)+1
             n = modulo(ndx(3)+kk-1+ng,ng)+1
             coor = coor + (disp(:,l,m,n)*wx*wy*wz)
          end do
       end do
    end do
    where(coor .lt. 0.0) coor = coor + boxsize
    where(coor .ge. boxsize) coor = coor - boxsize
  end subroutine sum_TSC

	!-------------------------------------------------------------
	!> Interpolate gridded values in rho to coor via TSC.
	!!
	!! Used to interpolate 3D grids (e.g. perturbed scale factor).
	!!
	!! @param coor Coordinate point where interpolation is evaluated
	!! @param disp 4D array to be interpolated to coor
	!! @param xscal Scaling of coordinates to "box scale"
	!! @param ng Number of grid points per dimension.
	!<------------------------------------------------------------
  function interpgrid_TSC(coor, rho, xscal, ng)
    real(dp), dimension(3), intent(in)       :: coor
    real(dp), dimension(:,:,:), intent(in)   :: rho
    real(dp), intent(in)                     :: xscal
    integer(i4b), intent(in)                 :: ng
    real(dp)                                 :: interpgrid_TSC
    integer(i4b), dimension(3) :: ndx
    real(dp), dimension(3)     :: dist
    integer(i4b)               :: ii,jj,kk
    integer(i4b) :: l,m,n
    real(dp) :: wx,wy,wz
    ndx=floor(coor/xscal)+1 ! Nearest grid point index to coordinate position
    dist=(coor/xscal) - real(ndx-1,dp) ! distance from each coordinate component
    interpgrid_TSC = 0.d0
    do ii=-1,1
       wx = tscweight(ii, dist(1))
       do jj=-1,1
          wy = tscweight(jj, dist(2))
          do kk=-1,1
             wz = tscweight(kk, dist(3))
             l = modulo(ndx(1)+ii-1+ng,ng)+1
             m = modulo(ndx(2)+jj-1+ng,ng)+1
             n = modulo(ndx(3)+kk-1+ng,ng)+1
             interpgrid_TSC = interpgrid_TSC + (rho(l,m,n)*wx*wy*wz)
          end do
       end do
    end do
  end function interpgrid_TSC

	!-------------------------------------------------------------
	!> Add i to ndx while imposing periodic boundary conditions.
	!<------------------------------------------------------------	
  pure function periodicGridIndex(ndx,i,ng) 
    integer, intent(in) :: ndx,i,ng
    integer             :: periodicGridIndex
    periodicGridIndex = ndx+i
    if (periodicGridIndex.ge.ng+1 .and. periodicGridIndex .lt. 2*ng) then 
       periodicGridIndex=periodicGridIndex-ng
    else if (periodicGridIndex.le.0 .and. periodicGridIndex .gt. -ng) then
       periodicGridIndex=periodicGridIndex+ng
    end if
  end function periodicGridIndex

	!-------------------------------------------------------------
	!> Compute the power spectrum of a gridded density field.
	!<------------------------------------------------------------	
  subroutine psestimator(ng, rho, ps, outfile)
    use fft_nsample
    integer(i4b), intent(in)                       :: ng
    real(dp), dimension(ng,ng,ng), intent(in)      :: rho
    real(dp), dimension(ng/2), intent(out)         :: ps
    character(len=*), optional                     :: outfile

    integer*8                              :: plan
    complex(cdp), dimension(ng/2+1,ng,ng)  :: rhoFT
    real(dp), dimension(0:ng-1)            :: pk,pkf
    integer(i4b), dimension(0:ng-1)        :: nk,nkf
    integer(i4b)                           :: ki,ii,jj,kk, i, unit_out,ios
    real(dp), dimension(ng)                :: mult
    real(dp)                               :: k

#ifndef FFTW
    real(dp), dimension(2*(ng+ng+ng))           :: trigs ! trig. weights
    real(dp), dimension(:), allocatable         :: work  ! workspace for fft
    integer(i8b), dimension(128*3)              :: ifac  ! factors of N
#endif

    print *, ""
    print *, "---Computing power spectrum estimator"

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
    print *, "  Deconvolving CIC window"
    mult(1) = 1.d0
    do ii=2,ng
       k = fftfreq(ii,ng)*TWOPI/real(ng,dp)/2.d0
       mult(ii) = 1.d0/(sin(k)/k)**2
    end do
    forall(ii=1:ng/2+1, jj=1:ng, kk=1:ng) rhoFT(ii,jj,kk) = &
         & mult(ii)*mult(jj)*mult(kk)*rhoFT(ii,jj,kk)
#endif

    ! ----- Bin to get 1D power spectrum
    pk = 0.0
    nk = 0
    print *, "  computing binned power spectrum...",ng
    do ii=1,ng/2+1
       do jj=1,ng
          do kk=1,ng
             ki = nint(sqrt(real(fftfreq(ii,ng)**2+fftfreq(jj,ng)**2+fftfreq(kk,ng)**2)))
             pk(ki) = pk(ki) + real(rhoFT(ii,jj,kk)*conjg(rhoFT(ii,jj,kk)),dp)
             nk(ki) = nk(ki) + 1
          end do
       end do
    end do
    print *, "  finished binned spectrum, normalizing..."
    ps = pk(1:ng/2)/nk(1:ng/2)
    ps = ps*V/real(ng,dp)**6

#ifdef FFTW
    call dfftw_cleanup()
#endif

    ! ----- Save to file if a filename argument was supplied
    if (present(outfile)) then
       print *, "Writing power spectrum to ",trim(outfile)
       call get_free_unit(unit_out)
       open(unit_out,file=trim(outfile),status="replace",action="write", &
            & iostat=ios, err=200)
       do i=1,ng/2
          write(unit_out,'(2ES18.8)',iostat=ios,err=210) i*2.*PI/boxsize,ps(i)
       end do
       close(unit_out)
    end if
    return
200 print *, "Error opening ps output file "//trim(outfile)
    print *, "ios = ",ios
    close(unit_out)
    stop
210 print *, "Error writing to file "//trim(outfile)
    print *, "ios = ",ios
    close(unit_out)
    stop
  end subroutine psestimator


#ifdef MPI
	!-------------------------------------------------------------
	!> Compute a density field, deltaL, with only those Fourier 
  !! modes with k <= nthresh*kF from delta.
	!<------------------------------------------------------------	
  subroutine select_largescale_modes(ng,delta)
    use fft_nsample
    integer(i4b), intent(in)                   :: ng
    real(dp), dimension(ng+fftpad,ng,local_nlastL), intent(inout) :: delta
    real(dp), dimension(ng+fftpad,ng,local_nlastL)                :: work
    integer(i4b)                               :: ii,jj,kk
    real(dp)                                   :: k

    ! ----- Forward FFT
    call rfftwnd_f77_mpi(planfL, 1, delta, work, 0, fftw_normal_order)
    do kk=1,local_nlastL
       do jj=1,ng
          do ii=1,ng/2+1
             k = sqrt(real(fftfreq(ii,ng)**2 + fftfreq(jj,ng)**2 + &
                  & fftfreq(kk+local_last_startL,ng)**2))
             if (k .gt. real(nthresh,dp)) then
                delta(2*ii-1,jj,kk) = 0.d0
                delta(2*ii,jj,kk) = 0.d0
             end if
          end do
       end do
    end do
    ! ----- Reverse FFT
    call rfftwnd_f77_mpi(planrL, 1, delta, work, 0, fftw_normal_order)
    delta = delta/ng**3
  end subroutine select_largescale_modes
#else ! MPI
	!-------------------------------------------------------------
	!> Compute a density field, deltaL, with only those Fourier 
  !! modes with k <= nthresh*kF from delta.
	!<------------------------------------------------------------
  subroutine select_largescale_modes(ng,delta,deltaL)
    use fft_nsample
    integer(i4b), intent(in)                   :: ng
    real(dp), dimension(ng,ng,ng), intent(in)  :: delta
    real(dp), dimension(ng,ng,ng), intent(out) :: deltaL
    complex(cdp), dimension(ng/2+1,ng,ng)      :: dFT
    integer(i4b)                               :: ii,jj,kk
    real(dp)                                   :: k
#ifdef FFTW
    integer(i8b)                            :: planf,planr
#else
    real(dp), dimension(ng+2,ng,ng)         :: tmp 
    real(dp), dimension(2*(ng+ng+ng))       :: trigsf,trigsb ! trig. weights
    real(dp), dimension(:), allocatable     :: workf,workb  ! workspace for fft
    integer(i8b), dimension(128*3)          :: ifacf,ifacb  ! factors of N
#endif FFTW

    ! ----- Forward FFT
#ifdef FFTW
    call dfftw_plan_dft_r2c_3d(planf,ng,ng,ng,delta,dFT,FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_3d(planr,ng,ng,ng,dFT,deltaL,FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(planf,delta,dFT)
#else
    ! initialize
    call fftpack_forward3d(int(0,i8b),ng,delta,dFT,ifacf,trigsf,workf,1.0d0)
    ! execute
    call fftpack_forward3d(int(-1,i8b),ng,delta,dFT,ifacf,trigsf,workf,1.0d0)
#endif FFTW
    do ii=1,ng/2+1
       do jj=1,ng
          do kk=1,ng
             k = sqrt(real(fftfreq(ii,ng)**2+fftfreq(jj,ng)**2+fftfreq(kk,ng)**2))
             if (k .gt. real(nthresh,dp)) dFT(ii,jj,kk) = cmplx(0._dp,0._dp,cdp)
          end do
       end do
    end do
    ! ----- Reverse FFT
#ifdef FFTW
    call dfftw_execute_dft_c2r(planr,dFT,deltaL)
    deltaL = deltaL/ng**3
#else
    ! initialize
    call fftpack_backward3d(int(0,i8b),ng,dFT,tmp,ifacb,trigsb,workb,1.0d0)
    ! execute
    call fftpack_backward3d(int(1,i8b),ng,dFT,tmp,ifacb,trigsb,workb,1.0d0/real(ng**3,dp))
    deltaL = tmp(1:ng,:,:)
#endif FFTW

#ifdef FFTW
    ! ----- Free plan memory
    call dfftw_destroy_plan(planf)
    call dfftw_destroy_plan(planr)
#endif FFTW
  end subroutine select_largescale_modes
#endif MPI

	!-------------------------------------------------------------
	!> Compute the derivative of f w.r.t. to coord (='x','y','z') and
  !! return the result in df.
	!<------------------------------------------------------------	
  subroutine deriv(df, f, coord, N)
    use fft_nsample
    integer(i4b), intent(in)                  :: N
    real(dp), dimension(N,N,N), intent(inout) :: df
    real(dp), dimension(N,N,N), intent(in)  :: f
    character(len=1), intent(in)            :: coord

    complex(cdp), dimension(N/2+1,N,N)      :: fFT
    real(dp), dimension(N/2+1,N,N)          :: k

#ifdef FFTW
    integer(i8b)                            :: planf,planr
#else
    real(dp), dimension(N+2,N,N)            :: dftmp 
    real(dp), dimension(2*(N+N+N))          :: trigsf,trigsb ! trig. weights
    real(dp), dimension(:), allocatable     :: workf,workb  ! workspace for fft
    integer(i8b), dimension(128*3)          :: ifacf,ifacb  ! factors of N
#endif

#ifdef FFTW
    call dfftw_plan_dft_r2c_3d(planf,N,N,N,f,fFT,FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_3d(planr,N,N,N,fFT,df,FFTW_ESTIMATE)
#endif

    ! ----- Choose the coordinate
    call fft_coord(k,coord,N)
    k = k*kF
    ! ----- Forward FFT
#ifdef FFTW
    call dfftw_execute(planf)
#else
    ! initialize
    call fftpack_forward3d(int(0,i8b),N,f,fFT,ifacf,trigsf,workf,1.0d0)
    ! execute
    call fftpack_forward3d(int(-1,i8b),N,f,fFT,ifacf,trigsf,workf,1.0d0)
#endif
    ! ----- Multiply by -i*kx, -i*ky, or -i*kz depending on coord value
    fFT = cmplx(aimag(fFT),-real(fFT),cdp) * k
    ! ----- Reverse FFT
#ifdef FFTW
    call dfftw_execute(planr)
    df = df/N**3
#else
    ! initialize
    call fftpack_backward3d(int(0,i8b),N,fFT,df,ifacb,trigsb,workb,1.0d0)
    ! execute
    call fftpack_backward3d(int(1,i8b),N,fFT,dftmp,ifacb,trigsb,workb,1.0d0/real(N**3,dp))
    df = dftmp(1:N,:,:)
#endif

#ifdef FFTW
    ! ----- Free plan memory
    call dfftw_destroy_plan(planf)
    call dfftw_destroy_plan(planr)
#endif
  end subroutine deriv

  subroutine fft_coord(k,coord,N)
    ! ------------------------------------------------------------
    ! Make a grid of wavenumber values
    ! ------------------------------------------------------------
    use fft_nsample, only: fftfreq
    integer(i4b), intent(in)                    :: N
    real(dp), dimension(N/2+1,N,N), intent(out) :: k
    character(len=1), intent(in)                :: coord
    integer(i4b)                                :: ii,jj,kk
    select case(coord)
    case('x')
       forall(ii=1:N/2+1, jj=1:N, kk=1:N) k(ii,jj,kk) = fftfreq(ii,N)
    case('y')
       forall(ii=1:N/2+1, jj=1:N, kk=1:N) k(ii,jj,kk) = fftfreq(jj,N)
    case('z')
       forall(ii=1:N/2+1, jj=1:N, kk=1:N) k(ii,jj,kk) = fftfreq(kk,N)
    case default
       write(6,*) "Error in deriv: invalid coord specified"
       write(6,*) coord
       stop
    end select
  end subroutine fft_coord


  subroutine deconvolve_CIC(ng, rho, in_local_nlast, in_local_last_start)
    ! ------------------------------------------------------------
    ! Deconvolve the 3D CIC window from the Fourier amplitudes in rho
    ! ------------------------------------------------------------
    use fft_nsample, only: local_nlast, local_last_start, fftfreq
    integer(i4b), intent(in) :: ng
!!$    real(dp), dimension(ng+2,ng,local_nlast), intent(inout) :: rho
    real(dp), dimension(:,:,:), intent(inout) :: rho
    integer(i4b), intent(in), optional :: in_local_nlast, in_local_last_start
    integer(i4b)                       :: tmp_local_nlast, tmp_local_last_start
    real(dp), dimension(ng)                     :: mult
    integer(i4b) :: ii,jj,kk
    real(dp) :: k, fac
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
    mult(1) = 1.d0
    do ii=2,ng
       k = fftfreq(ii,ng)*TWOPI/real(ng,dp)/2.d0
       mult(ii) = 1.d0/(sin(k)/k)**2
    end do
    do kk=1,tmp_local_nlast
       do jj=1,ng
          do ii=1,ng/2+1
             fac = mult(ii)*mult(jj)*mult(kk+tmp_local_last_start)
             rho(2*ii-1,jj,kk) = fac * rho(2*ii-1,jj,kk)
             rho(2*ii,jj,kk) = fac * rho(2*ii,jj,kk)
          end do
       end do
    end do
  end subroutine deconvolve_CIC


  subroutine draw_random_array(rnarray, n, iseed)
    ! ------------------------------------------------------------
    ! Draw an array of random numbers to be distributed among 
    ! MPI tasks (to make sure there is no repetition between
    ! the pseudo-random sequences on each processor).
    ! ------------------------------------------------------------
    integer(i4b), intent(in)            :: n, iseed
    real(dp), dimension(n), intent(out) :: rnarray
    real(dp), dimension(n*numtasks)     :: harvest
    integer, dimension(1) :: seed
    if (rank .eq. 0) then
       seed(1) = iseed
       call random_seed(put=seed)
       call random_number(harvest)
    end if
#ifdef MPI
    call mpi_scatter(harvest, n, mpi_double_precision, rnarray, n, mpi_double_precision, &
         & 0, mpi_comm_world, ierr)
#endif
  end subroutine draw_random_array


  subroutine smoothDensityField(densGrid, smoothGrid, ng, smoothLength)
    ! ------------------------------------------------------------
    ! Convolve the gridded density, densGrid, with a smoothing
    ! window and return an array  of smoothed density values in
    ! smoothGrid.  
    ! 
    ! Currently only smoothing in cubes of sidelength 
    ! gridsize * smoothLength is supported.
    ! ------------------------------------------------------------
    integer(i4b), intent(in)                  :: ng, smoothLength
    real(dp), dimension(ng,ng,ng), intent(in) :: densGrid
    real(dp), dimension(ng,ng,ng), intent(out)   :: smoothGrid
    integer(i4b) :: ishift, idim
    if (smoothLength .lt. 0 .or. smoothLength .gt. ng/2-1) then
       print *, "ERROR: bad value of smoothLength in smoothDensityField:",smoothLength,ng
       stop
    else
       if (smoothLength .eq. 0) then
          smoothGrid = densGrid
       else
          smoothGrid = 0.d0
          do ishift = -smoothLength, smoothLength, 1
             do idim = 1,3
                smoothGrid = smoothGrid + cshift(densGrid, shift=ishift, dim=idim)
             end do
          end do
       end if
       smoothGrid = smoothGrid / real((2*smoothLength+1)**3,dp)
    end if
  end subroutine smoothDensityField

end module gridtools
