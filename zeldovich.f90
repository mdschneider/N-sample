!
!> Compute Zeldovich displacements from a gridded density field
!! and add or subtract those displacements from a particle list
!! (using CIC interpolation from the grid points).
!
!  Created by Michael Schneider on 2010-01-05
!
module zeldovich
  use types_nsample
  use parameters

contains
#ifdef MPI
	!-------------------------------------------------------------
	!> Compute the Zeldovich displacements from the density field delta
	!<------------------------------------------------------------	
  subroutine zeldovich_displacements(ng,np,delta,disp)
    use fft_nsample
    use gridtools, only: print_gridarray
    integer(i4b), intent(in)                     :: ng,np
    real(dp), dimension(ng+fftpad,ng,local_nlastL), intent(in) :: delta
    real(dp), dimension(3,ng,ng,local_nlastL), intent(out):: disp

    real(dp), dimension(3)                  :: kv
    real(dp)                                :: ksq
    integer(i4b)                            :: ii,jj,kk,i1,i2
    real(dp), dimension(ng+fftpad,ng,local_nlastL) :: dx,dy,dz,work, delta_loc

#ifdef VERBOSE
    print *, rank,"  Computing Zeldovich displacements"
#endif
    delta_loc = delta

    ! ----- Compute gridded displacements via FFT
    call mpi_barrier(mpi_comm_world, ierr)
    call rfftwnd_f77_mpi(planfL, 1, delta_loc, work, 0, fftw_normal_order)
    do kk=1,local_nlastL
       do jj=1,ng
          do ii=1,ng/2+1
             kv(1)=real(fftfreq(ii,ng),dp)*kF
             kv(2)=real(fftfreq(jj,ng),dp)*kF
             kv(3)=real(fftfreq(kk+local_last_startL,ng),dp)*kF
             ksq = dot_product(kv,kv)
             i1 = 2*ii-1
             i2 = 2*ii
             if (ksq .eq. 0.0) then
                dx(i1:i2,jj,kk) = 0.d0
                dy(i1:i2,jj,kk) = 0.d0
                dz(i1:i2,jj,kk) = 0.d0
             else
                dx(i1:i2,jj,kk) = kv(1) * (/-delta_loc(i2,jj,kk), &
                     & delta_loc(i1,jj,kk)/) / ksq
                dy(i1:i2,jj,kk) = kv(2) * (/-delta_loc(i2,jj,kk), &
                     & delta_loc(i1,jj,kk)/) / ksq
                dz(i1:i2,jj,kk) = kv(3) * (/-delta_loc(i2,jj,kk), &
                     & delta_loc(i1,jj,kk)/) / ksq
             end if
          end do
       end do
    end do
    ! ----- Reverse FFT
    call rfftwnd_f77_mpi(planrL, 1, dx, work, 0, fftw_normal_order)
    dx = dx/real(ng,dp)**3
    call rfftwnd_f77_mpi(planrL, 1, dy, work, 0, fftw_normal_order)
    dy = dy/real(ng,dp)**3
    call rfftwnd_f77_mpi(planrL, 1, dz, work, 0, fftw_normal_order)
    dz = dz/real(ng,dp)**3
    call print_gridarray(dx(1:ng,:,:), "dx")
    call print_gridarray(dy(1:ng,:,:), "dy")
    call print_gridarray(dz(1:ng,:,:), "dz")
    ! ----- Assign output array
    disp(1,:,:,:) = dx(1:ng,:,:)
    disp(2,:,:,:) = dy(1:ng,:,:)
    disp(3,:,:,:) = dz(1:ng,:,:)
  end subroutine zeldovich_displacements
#else
	!-------------------------------------------------------------
	!> Compute the Zeldovich displacements from the density field delta
	!<------------------------------------------------------------
  subroutine zeldovich_displacements(ng,np,delta,disp)
    use fft_nsample
    integer(i4b), intent(in)                     :: ng,np
    real(dp), dimension(ng,ng,ng), intent(in)    :: delta
    real(dp), dimension(3,ng,ng,ng), intent(out) :: disp

    complex(cdp)                              :: ideltak
    real(dp), dimension(3)                    :: kv
    real(dp)                                  :: ksq
    integer(i4b)                              :: ii,jj,kk
    complex(cdp), dimension(3,ng/2+1,ng,ng)   :: dispFT
    complex(cdp), dimension(ng/2+1,ng,ng)     :: tmpFT,dFT
#ifdef FFTW
    integer(i8b)                            :: planf,planr
    real(dp), dimension(ng,ng,ng)           :: tmp
#else
    real(dp), dimension(ng+2,ng,ng)         :: tmp
    real(dp), dimension(2*(ng+ng+ng))       :: trigsf,trigsb ! trig. weights
    real(dp), dimension(:), allocatable     :: workf,workb  ! workspace for fft
    integer(i8b), dimension(128*3)          :: ifacf,ifacb  ! factors of N
#endif

    print *, ""
    print *, "Computing Zeldovich displacements"

#ifdef FFTW
    ! ----- FFTW plans 
    call dfftw_plan_dft_r2c_3d(planf,ng,ng,ng,delta,dFT,FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_3d(planr,ng,ng,ng,tmpFT,tmp,FFTW_ESTIMATE)
    ! ----- Compute gridded displacements via FFT
    call dfftw_execute_dft_r2c(planf,delta,dFT)
#else
    ! initialize
    call fftpack_forward3d(int(0,i8b),ng,delta,dFT,ifacf,trigsf,workf,1.0d0)
    ! execute
    call fftpack_forward3d(int(-1,i8b),ng,delta,dFT,ifacf,trigsf,workf,1.0d0)
#endif
    do kk=1,ng
       do jj=1,ng
          do ii=1,ng/2+1
             kv(1)=real(fftfreq(ii,ng),dp)*kF
             kv(2)=real(fftfreq(jj,ng),dp)*kF
             kv(3)=real(fftfreq(kk,ng),dp)*kF
             ksq = dot_product(kv,kv)
             if (ksq .eq. 0.0) then
                dispFT(:,ii,jj,kk) = cmplx(0.d0,0.d0,cdp)
             else
                ideltak = -cmplx(aimag(dFT(ii,jj,kk)),-real(dFT(ii,jj,kk)))/ksq
                dispFT(:,ii,jj,kk) = kv*ideltak
             end if
          end do
       end do
    end do
    ! ----- Reverse FFT
#ifndef FFTW
    call fftpack_backward3d(int(0,i8b),ng,tmpFT,tmp,ifacb,trigsb,workb,1.0d0)
#endif
    do ii=1,3
       tmpFT = dispFT(ii,:,:,:)
#ifdef FFTW
       call dfftw_execute_dft_c2r(planr,tmpFT,tmp)
       disp(ii,:,:,:) = tmp/ng**3
#else
       ! execute
       call fftpack_backward3d(int(1,i8b),ng,tmpFT,tmp,ifacb,trigsb,workb,1.0d0/real(ng**3,dp))
       disp(ii,:,:,:) = tmp(1:ng,:,:)
#endif
    end do
#ifdef FFTW
    ! ----- Free plan memory
    call dfftw_destroy_plan(planf)
    call dfftw_destroy_plan(planr)
#endif
  end subroutine zeldovich_displacements
#endif

	!-------------------------------------------------------------
	!> Interpolate displacements to particle positions and subtract
  !! from both positions and velocities in coor.
	!<------------------------------------------------------------	
  subroutine subtract_particle_displacements(ng,np,a,boxsize,disp,coor)
    use cosmology
    integer(i4b), intent(in)                    :: ng,np
    real(dp), intent(in)                        :: a,boxsize
    real(dp), dimension(3,ng,ng,ng), intent(in) :: disp
    real(dp), dimension(6,np), intent(inout)    :: coor

    real(dp)                                    :: f,g,xscal
    integer(i4b)                                :: i,ii,jj,kk,l,m,n
    integer(i4b), dimension(3,np)               :: ndx
    real(dp), dimension(3,np)                   :: dist
    real(dp), dimension(np)                     :: vol

    print *, ""
    print *, "---Subtracting displacements at particle positions"

    g = lingrowth(a) / lingrowth(1.d0)
    f = (dlnDdlna(a)/a - 1.d0)*H_H0(a)*g
    ! ----- Interpolate displacements to particle positions and subtract
    xscal = boxsize/real(ng)
    ndx=int(coor(1:3,:)/xscal)+1 ! Nearest grid point index to coordinate position
    dist=real(ndx,dp)-(coor(1:3,:)/xscal) ! 1 - distance from each coordinate component
    do ii=0,1
       do jj=0,1
          do kk=0,1
             vol=abs(real(ii,dp)-dist(1,:))*&
                  & abs(real(jj,dp)-dist(2,:))*&
                  & abs(real(kk,dp)-dist(3,:))
             do i = 1,np
                l = modulo(ndx(1,i)+ii-1,ng)+1
                m = modulo(ndx(2,i)+jj-1,ng)+1
                n = modulo(ndx(3,i)+kk-1,ng)+1
                coor(1:3,i) = coor(1:3,i) - (g*disp(:,l,m,n)*vol(i))
                coor(4:6,i) = coor(4:6,i) - (f*disp(:,l,m,n)*vol(i))
             end do
          end do
       end do
    end do
    where(coor(1:3,:) .lt. 0.0) coor(1:3,:) = coor(1:3,:) + boxsize
    where(coor(1:3,:) .ge. boxsize) coor(1:3,:) = coor(1:3,:) - boxsize
  end subroutine subtract_particle_displacements

	!-------------------------------------------------------------
	!> Interpolate displacements to particle positions and add to 
  !! particle positions within the specified range.
	!<------------------------------------------------------------	
  subroutine add_particle_displacements(ng,np,a,boxsize,disp,coor)
    use cosmology
    use gridtools, only: sum_CIC, sum_TSC
    integer(i4b), intent(in)                    :: ng,np
    real(dp), intent(in)                        :: a,boxsize
    real(dp), dimension(3,ng,ng,ng), intent(in) :: disp
    real(dp), dimension(:,:), intent(inout)     :: coor

    real(dp)                                    :: f,xscal,g
    integer(i4b)                                :: i,ii,jj,kk,l,m,n
    logical :: velocities

    velocities = .false.
    if (size(coor(:,1)) .eq. 6) velocities = .true.

    g = 1.d0 !lingrowth(a) / lingrowth(1.d0)
    f = (dlnDdlna(a)/a - 1.d0)*H_H0(a)*g
    xscal = boxsize/real(ng)

    print *, "add_particle_displacements, g:",g," f:",f," xscal:",xscal

    do i = 1,np
       call sum_TSC(coor(1:3,i), disp, xscal, ng)
    end do
!!$    if (velocities) then
!!$       do i = 1,np
!!$          call sum_TSC(coor(4:6,i), disp*f, xscal, ng)
!!$       end do
!!$    end if
  end subroutine add_particle_displacements

end module zeldovich
