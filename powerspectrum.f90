!
!> Compute the power spectrum from an L-Gadget snapshot file.
!
!  Created by Michael Schneider on 2010-04-09
!
program powerspectrum
  use types_nsample
  use parameters
  use io_nsample
  use fft_nsample
  use gridtools, only: part2grid, print_gridarray
  use resample, only: npkbins, power_spectrum_estimator
  implicit none

  character(len=200)                       :: arg, outdir, outfile
  real(dp), dimension(:,:), allocatable    :: coor, incoor
  integer(i8b), dimension(npart)           :: pid
  real(dp)                                 :: scale_factor
  integer(i4b)                             :: ng, i,j, isnap
  real(dp), dimension(npk,2)               :: pk
  real(dp), dimension(:), allocatable      :: ps
  real(dp), dimension(:,:,:), allocatable  :: delta
  character(len=10) :: isnaplab
  logical :: icfile

#ifdef MPI
  ! ----- Initialize MPI
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,rank,ierr)
  call mpi_comm_size(mpi_comm_world,numtasks,ierr)
  print *, "Process ",rank," of ",numtasks," is alive"
#else
  rank = 0
  numtasks = 1
#endif MPI

  icfile = .false.


  ! ----- Parse command line arguments
  if (rank .eq. 0) then
     if (command_argument_count() /= 3) then
        write(6,*) "Usage: powerspectrum.x ng isnap outdir"
        stop
     end if
     call get_command_argument(1,arg)
     read(arg,*) ng
     call get_command_argument(2,arg)
     read(arg,*) isnap
     call get_command_argument(3,outdir)
  end if
#ifdef MPI  
  call mpi_bcast(ng, 1, mpi_integer, 0, mpi_comm_world, ierr)
  call mpi_bcast(isnap, 1, mpi_integer, 0, mpi_comm_world, ierr)
  call mpi_bcast(outdir, 200, mpi_character, 0, mpi_comm_world, ierr)
#endif

  ! ----- Allocate output power spectrum array
  allocate(ps(0:npkbins(ng)))


  ! ----- Initialize FFTW plans and MPI slab decomposition
#ifdef MPI
  call init_fftw(rank, numtasks, ng, ng)
  if (rank .eq. 0) print *, "Finished FFTW init"
#else !MPI
  local_nlast = ng
  local_last_start = 0
  local_nlastL = ng
  local_last_startL = 0
#endif MPI

  call parameters_init()

  ! ----- Read snapshot and distribute particles
  allocate(coor(3, npart/numtasks))
  if (rank .eq. 0) then
     allocate(incoor(3,npart))
     call read_snapshot(outdir, isnap, incoor, pid, scale_factor, icfile=icfile)
     print *, "scale factor of input snapshot: ",scale_factor
  end if
#ifdef MPI
  call mpi_bcast(scale_factor, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
  call mpi_scatter(incoor, 3*npart/numtasks, mpi_double_precision, &
       & coor, 3*npart/numtasks, mpi_double_precision, &
       & 0, mpi_comm_world, ierr)
  call mpi_barrier(mpi_comm_world, ierr)
#else
  coor = incoor
#endif MPI
  if (rank .eq. 0) deallocate(incoor)

  allocate(delta(ng+fftpad, ng, local_nlast))
  call part2grid(ng, npart, boxsize, coor(:,1:npart/numtasks), delta)
  deallocate(coor)
  call print_gridarray(delta, "snapshot density")


  ! ----- Compute power spectrum
  write(isnaplab,'(i3.3)') isnap
  if (pkbintype .eq. 1) then
     outfile = trim(outdir)//"/snapdir_"//trim(adjustl(isnaplab))// &
          & "/ps_takbins_"//trim(adjustl(isnaplab))//".txt"
  else
     if (icfile) then
        outfile = trim(outdir)// &
             & "/ps_"//trim(adjustl(isnaplab))//".txt"
     else
        outfile = trim(outdir)//"/snapdir_"//trim(adjustl(isnaplab))// &
             & "/ps_"//trim(adjustl(isnaplab))//".txt"
     end if
  end if
  call power_spectrum_estimator(ng, ps, delta, 1, 1, outfile)
  deallocate(delta, ps)  
  

  stop

end program powerspectrum
