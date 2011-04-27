!
!> Draw large-scale modes and add them to L-Gadget run without
!! large-scale power by shifting particles using Zeldovich 
!! displacements at perturbed scale factor time and then gridding 
!! to get deltaS.
!!
!! This is a simple routine to add just a single large-scale mode realization.
!
!  Created by Michael Schneider on 2009-11-27
!
program addPowerSnapshots
  use types_nsample
  use parameters
  use gridtools
  use zeldovich
  use resample
  use cosmology, only: lingrowth
  use fft_nsample, only: local_nlast, local_last_start
  implicit none

  character(len=len_filename)              :: arg, densdir, outdir, outfile, infile
!!$  real(dp), dimension(ntable,2)            :: atable
  real(dp), dimension(:,:,:), allocatable  :: deltaL1, deltaL2, delta
  real(dp), dimension(:,:,:,:), allocatable :: disp
  integer(i4b)                             :: ng, ngL, iseed

  real(dp), dimension(npk,2)               :: pk
  real(dp), dimension(:), allocatable      :: ps

  real(dp) :: scale_factor
  logical  :: new_deltaL_draw, sorted_snapshots_algorithm, dotplot


  real(dp), dimension(:,:), allocatable :: coor
  integer(i8b), dimension(:), allocatable :: pid
  real(dp) :: u
  integer :: i,j, unit_out, npart_read
  real(dp), dimension(:,:,:), allocatable  :: psnap


  ! ----- Parse command line arguments
  if (command_argument_count() .ne. 4) then
     write(6,*) "Usage: addPowerSnapshots ng ngL outdir new_deltaL_draw"
     stop
  end if
  call get_command_argument(1,arg)
  read(arg,*) ng
  call get_command_argument(2,arg)
  read(arg,*) ngL
  call get_command_argument(3,outdir)
  call get_command_argument(4,arg)
  read(arg,*) new_deltaL_draw
  print *, "ng:",ng
  print *, "outdir: ",trim(outdir)
  print *, "new_deltaL_draw: ",new_deltaL_draw

  rank = 0
  numtasks = 1
  local_nlast = ng
  local_last_start = 0
  local_nlastL = ngL
  local_last_startL = 0

  call parameters_init()

  dotplot = .false.
  dotplotdir = trim(outdir)//"/dotplots"

  ! Use the algorithm that requires pre-sorted snapshot input?
  ! (Much faster!)
  sorted_snapshots_algorithm = .true.

  ! ----- Allocate output power spectrum array
  allocate(ps(ng/2))
!!$  allocate(ps(0:npkbins(ngL)))

  ! ----- Read theory power spectrum from file
  call read_powerspectrum(pk)
  print *, "Finished reading theory power spectrum"


!!$  ! ----- Read allmodes at a=1 just for saving dotplot particles
!!$  if (dotplot) then
!!$     allocate(coor(3,sim%npart))
!!$     call read_snapshot("/data/rw7/bskm46/work/covResampling/output/takahashi/old_allmodes", 15 ,coor, pid, scale_factor)
!!$     print *, "scale factor of allmodes input (at a=1):",scale_factor
!!$     call save_dotplot_particles(dotplotdir, "allmodes_a1.txt", 10.d0, coor, 1)
!!$     deallocate(coor)
!!$  end if


  ! ----- Get the large-scale modes on a grid
  allocate(deltaL1(ngL,ngL,ngL))
  if (new_deltaL_draw) then
     ! ----- Draw large-scale modes
     iseed = -1
     allocate(deltaL2(ngL,ngL,ngL))
     print *, "Drawing large-scale modes"
     call draw_modes(deltaL1, deltaL2, ngL, pk, iseed)
     deallocate(deltaL2)
  else
     allocate(coor(3,sim%npart))
     allocate(delta(ngL,ngL,ngL))
     ! Get deltaL from earliest snapshot and then evolve
     ! with the linear growth function.
     ! FIXME: change this to read the positions from the IC file
     allocate(pid(sim%npart))
     call read_snapshot(allmodesdir, 0, coor, pid, scale_factor, insim=sim)
     deallocate(pid)
     print *, "scale factor of allmodes input: ",scale_factor
     print *, "x coor range:",minval(coor(1,:)),maxval(coor(1,:))
     print *, "y coor range:",minval(coor(2,:)),maxval(coor(2,:))
     print *, "z coor range:",minval(coor(3,:)),maxval(coor(3,:))
     call part2grid(ngL, sim%npart, sim%boxsize, coor, delta)
     delta = delta - 1.d0
     delta = delta * lingrowth(1.d0) / lingrowth(scale_factor)
     deallocate(coor)
     call print_gridarray(delta, "input allmodes")
     call select_largescale_modes(ngL, delta, deltaL1)
     deallocate(delta)
  endif
  call print_gridarray(deltaL1, "deltaL")


!!$  ! ----- Save deltaL grid
!!$  outfile = trim(outdir)//"/deltaL_grid.txt"
!!$  print *, "Saving large-scale modes to file: ",trim(outfile)
!!$  call save_density(outfile, deltaL1, ng)
!!$

  ! ----- Compute power spectrum of large-scale modes
  call psfilename(outfile, "deltaL", outdir, ng)
  call psestimator(ngL, deltaL1, ps, outfile)
  call power_spectrum_estimator(ngL, ps, deltaL1, 1, 1, outfile)
  
  if (sorted_snapshots_algorithm) then
     ! ----- Read snapshot particle coordinates
!!$     call sort_snapshot_particles(32**3)
!!$     stop
     npart_read = sim%npart / numtasks
     print *, "npart_read:",npart_read
     allocate(psnap(3, nsnap, npart_read))
     call input_sorted_snapshot_particles(psnap, 1 + rank*npart_read, npart_read)
     print *, "psnap range:",minval(psnap),maxval(psnap)
     if (dotplot) then
        ! Save particles in smallmodes at a=1 for dotplotting
        allocate(coor(3,sim%npart))
        do i = 1,sim%npart
           coor(:,i) = psnap(:,isnap_a1+1,i)
        end do
        print *, "coor range:",minval(coor),maxval(coor)
        call save_dotplot_particles(dotplotdir, "smallmodes_a1.txt", &
             & 10.d0, coor, 1)
        deallocate(coor)
     end if
     call read_a_omega_table()
     call gen_perturbed_density(deltaL1, psnap, ng, ngL) ! FIXME: compile error on intel
  else
     ! ----- Compute Zeldovich displacements at z = 0
     allocate(disp(3,ng,ng,ng))
     call zeldovich_displacements(ng, sim%npart, deltaL1, disp)
     call print_gridarray(disp(1,:,:,:), "x disp")
     call print_gridarray(disp(2,:,:,:), "y disp")
     call print_gridarray(disp(3,:,:,:), "z disp")

     ! ----- Convert deltaL into perturbed scale factor grid
     print *, "Converting delta to perturbed scale factor"
     call read_a_omega_table()
     call convert_density_to_scale_factor(deltaL1, ng)
     call print_gridarray(deltaL1, "scale factor")  

     ! ----- Read snapshots at perturbed scale factor and 
     !       compute the density.
     print *, "Calculating particle positions at perturbed times..."
     allocate(coor(3,sim%npart))
     call interp_part_pos(ng, deltaL1, coor, disp, .true.)


     ! ##### TESTING ####################
!!$     ! ##### Test whether the Zeldovich move alone creates the correct 
!!$     !       large-scale power.
!!$     allocate(coor(3,sim%npart))
!!$     do i=1,3
!!$        do j=1,sim%npart
!!$           call random_number(u)
!!$           coor(i,j) = u * 1000.d0
!!$        end do
!!$     end do
!!$     ! Instead of random coordinates above, add Zeldovich displacements
!!$     ! to the deltaS coordinates at z = 0
!!$     allocate(pid(sim%npart))
!!$     call read_snapshot(smallmodesdir, isnap_a1, coor, pid, scale_factor)
!!$     deallocate(pid)
!!$     call add_particle_displacements(ng, sim%npart, 1.d0, sim%boxsize, disp, coor)
     ! ##### END TESTING ################

     ! ----- Compute gridded deltaS
!!$     call save_dotplot_particles(dotplotdir, "perturbed_positions.txt", &
!!$          & 10.d0, coor, 1)
     call part2grid(ng, sim%npart, sim%boxsize, coor, deltaL1)
     call print_gridarray(deltaL1, "deltaS")
  end if

!!$  ! ----- Set delta to zero except for 1 subcube
!!$  deltaL1(:,:,64:512) = 0.d0
!!$  deltaL1(:,64:512,:) = 0.d0
!!$  deltaL1(64:512,:,:) = 0.d0
!!$  deltaL1(64:512,64:512,:) = 0.d0
!!$  deltaL1(64:512,:,64:512) = 0.d0
!!$  deltaL1(:,64:512,64:512) = 0.d0
!!$  deltaL1(64:512,64:512,64:512) = 0.d0
!!$  deltaL1 = deltaL1 - 1.d0

  ! ----- Compute deltaS power spectrum
  call psfilename(outfile, "deltaS", outdir, ng)
  call psestimator(ng, deltaL1, ps, outfile)
  print *, "All done!"

  stop

end program addPowerSnapshots
