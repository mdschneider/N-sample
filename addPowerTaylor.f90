!
!> Draw large-scale modes and add them to L-Gadget run without
!! large-scale power by adding the first term in the Taylor 
!! expansion of deltaS(x',a'), which depends on the growth 
!! rate of deltaS.  All calculations are performed using 
!! gridded densities (rather than particle positions).
!
!  Created by Michael Schneider on 2010-02-11
!
program addPowerTaylor
  use types_nsample
  use parameters
  use gridtools
  use zeldovich
  use resample
  use cosmology
  implicit none

  character(len=200)                       :: arg, densdir, outdir, outfile, infile
!!$  real(dp), dimension(ntable,2)            :: atable
  real(dp), dimension(:,:,:), allocatable  :: deltaL1, deltaL2, delta, tmp, vx, vy, vz, dv
  real(dp), dimension(:,:,:,:), allocatable :: disp
  integer(i4b)                             :: ng, iseed

  real(dp), dimension(npk,2)               :: pk
  real(dp), dimension(:), allocatable      :: ps

  real(dp) :: scale_factor
  logical  :: new_deltaL_draw


  real(dp), dimension(:,:), allocatable :: coor
  integer(i8b), dimension(:), allocatable :: pid
  real(dp) :: u
  integer :: i,j, unit_out

  ! Histogram 
  integer(i4b), parameter :: nbins = 2001 !1120
  real(dp), parameter :: binmin = -200.d0 !-235.d0
  integer(i8b), dimension(nbins) :: bincounts

  ! ----- Parse command line arguments
  if (command_argument_count() /= 3) then
    write(6,*) "Usage: addPowerTaylor ng outdir new_deltaL_draw"
    stop
  end if
  call get_command_argument(1,arg)
  read(arg,*) ng
  call get_command_argument(2,outdir)
  call get_command_argument(3,arg)
  read(arg,*) new_deltaL_draw
  print *, "ng:",ng
  print *, "outdir: ",trim(outdir)
  print *, "new_deltaL_draw: ",new_deltaL_draw

  call parameters_init()

  ! ----- Allocate output power spectrum array
  allocate(ps(ng/2))


  ! ------------------------------------------------------------
  ! Large-scale modes
  ! ------------------------------------------------------------
  ! ----- Read power spectrum from file
  call read_powerspectrum(pk)
  print *, "Finished reading theory power spectrum"

  ! ----- Get the large-scale modes on a grid
  allocate(deltaL1(ng,ng,ng))
  if (new_deltaL_draw) then
     ! ----- Draw large-scale modes
     iseed = -1
     allocate(deltaL2(ng,ng,ng))
     print *, "Drawing large-scale modes"
     call draw_modes(deltaL1, deltaL2, ng, pk, iseed)
     deallocate(deltaL2)
  else
     allocate(coor(3,sim%npart))
     allocate(delta(ng,ng,ng))
     allocate(pid(sim%npart))
     call read_snapshot(allmodesdir, 0, coor, pid, scale_factor)
     deallocate(pid)
     print *, "scale factor of allmodes input: ",scale_factor
     call part2grid(ng, sim%npart, sim%boxsize, coor, delta)
     delta = delta - 1.d0
     delta = delta * lingrowth(1.d0) / lingrowth(scale_factor)
     deallocate(coor)
     call select_largescale_modes(ng, delta, deltaL1)
     deallocate(delta)
  endif
  call print_gridarray(deltaL1, "deltaL")

!!$  ! ----- Save deltaL grid
!!$  outfile = trim(outdir)//"/deltaL"
!!$  print *, "Saving large-scale modes to file: ",trim(outfile)
!!$  call save_density(outfile, deltaL1, ng)

  ! ----- Compute power spectrum of large-scale modes
  call psfilename(outfile, "deltaL", outdir, ng)
  call psestimator(ng, deltaL1, ps, outfile)


  ! ----- Compute Zeldovich displacements at z = 0
  allocate(disp(3,ng,ng,ng))
  call zeldovich_displacements(ng, sim%npart, deltaL1, disp)
  call print_gridarray(disp(1,:,:,:), "x disp")
  call print_gridarray(disp(2,:,:,:), "y disp")
  call print_gridarray(disp(3,:,:,:), "z disp")


  ! ----- Add Zeldovich displacements to smallmodes simulation
  allocate(coor(6,sim%npart))
  allocate(pid(sim%npart))
  call read_snapshot(smallmodesdir, isnap_a1, coor, pid, scale_factor)
  deallocate(pid)
  print *, "scale factor of smallmodes input: ",scale_factor
  call add_particle_displacements(ng, sim%npart, scale_factor, sim%boxsize, disp, coor)
  deallocate(disp)


!!$  ! ----- Perturbed scale factor
!!$  call convert_density_to_scale_factor(deltaL1, ng)
!!$  deltaL1 = ((deltaL1/scale_factor) - 1.d0) / (H_H0(scale_factor)*100.d0)
!!$  call print_gridarray(deltaL1, "scale factor")  
!!$  call add_linear_time_evolution(coor, ng, tmp)


  ! ------------------------------------------------------------
  ! Small-scale growth rate
  ! ------------------------------------------------------------
  allocate(tmp(ng,ng,ng))
  allocate(delta(ng,ng,ng), vx(ng,ng,ng),vy(ng,ng,ng),vz(ng,ng,ng))
  call part2grid(ng, sim%npart, sim%boxsize, coor, delta, vx, vy, vz)
  call print_gridarray(delta, "input deltaS")
  call print_gridarray(vx, "vx")
  call print_gridarray(vy, "vy")
  call print_gridarray(vz, "vz")
  
  ! ----- Compute divergence of -(1+delta)*v in comoving coords via FFT
  print *, ""
  print *, "Computing derivatives"
  allocate(dv(ng,ng,ng))
  dv = 0.0_dp
  call deriv(tmp,vx*delta,'x',ng)
  dv = dv + tmp
  call deriv(tmp,vy*delta,'y',ng)
  dv = dv + tmp
  call deriv(tmp,vz*delta,'z',ng)
  dv = dv + tmp
  dv = dv/(H_H0(scale_factor)*100._dp) ! divide by H(a)
  call print_gridarray(dv,"dv")
  delta = delta - 1.d0

  dv = (deltaL1 * (dv * dlnDdlnOmega(scale_factor) / dlnDdlna(scale_factor)))
  call print_gridarray(dv, "full correction term")

  print *, "Computing histogram of full correction term"
  call hist_grid_values(dv / delta, ng, 0.2_dp, binmin, bincounts, nbins)
  write(90,*) bincounts


!!$  ! ----- Gradient of delta dotted into velocity
!!$  call deriv(tmp, delta, 'x', ng, kF)
!!$  dv = dv + tmp*vx/(H_H0(scale_factor)*100.d0)
!!$  call deriv(tmp, delta, 'y', ng, kF)
!!$  dv = dv + tmp*vy/(H_H0(scale_factor)*100.d0)
!!$  call deriv(tmp, delta, 'z', ng, kF)
!!$  dv = dv + tmp*vz/(H_H0(scale_factor)*100.d0)
!!$  deallocate(vx,vy,vz)
!!$  call print_gridarray(dv,"dv + gradient")
!!$
!!$  ! ----- Power spectrum of growth rate
!!$  call psfilename(outfile, "growthrate", outdir, ng)
!!$  call psestimator(ng, dv, ps, outfile)


  ! ------------------------------------------------------------
  ! Reconstruction
  ! ------------------------------------------------------------
  !tmp = delta + (deltaL1 * (1.d0 + dv * dlnDdlnOmega(scale_factor) / dlnDdlna(scale_factor)))
  tmp = delta + deltaL1 + dv
  call print_gridarray(tmp, "final delta")
  
  ! ----- Power spectrum of deltaS
  call psfilename(outfile, "deltaS", outdir, ng)
  call psestimator(ng, tmp, ps, outfile)

!!$  ! ----- Power spectrum of 1st order Taylor term
!!$  tmp = deltaL1*dv*dlnDdlnOmega(scale_factor)/dlnDdlna(scale_factor)
!!$  call psfilename(outfile, "TaylorTerm", outdir, ng)
!!$  call psestimator(ng, tmp, ps, outfile)

  stop

end program addPowerTaylor
