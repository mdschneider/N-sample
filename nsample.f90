!
!> Main driver routine.
!!
!! Draw large-scale modes and add them to N-body simulation run 
!! without large-scale power.
program nsample
  use types_nsample
  use parameters
  use io_nsample
  use fft_nsample
  use resample

  integer(i4b)                             :: ng, ngL, nr, i, iseed
  character(len=len_filename)              :: outdir, arg
  
  real(dp), dimension(npk,2)               :: pk
  real(dp), dimension(:,:,:), allocatable  :: psnap


#ifdef MPI
  ! ----- Initialize MPI
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,rank,ierr)
  call mpi_comm_size(mpi_comm_world,numtasks,ierr)
!!$  print *, "Process ",rank," of ",numtasks," is alive"
#else
  rank = 0
  numtasks = 1
#endif MPI


  ! ----- Parse command line arguments
  if (rank .eq. 0) then
     if (command_argument_count() /= 6) then
        write(6,*) "Usage: nsample.x ng ngL nr ir_start iseed outdir"
        stop
     end if
     call get_command_argument(1,arg)
     read(arg,*) ng
     call get_command_argument(2,arg)
     read(arg,*) ngL
     call get_command_argument(3,arg)
     read(arg,*) nr
     call get_command_argument(4,arg)
     read(arg,*) ir_start
     call get_command_argument(5,arg)
     read(arg,*) iseed
     call get_command_argument(6,outdir)
  end if
#ifdef MPI
  call mpi_bcast(ng, 1, mpi_integer, 0, mpi_comm_world, ierr)
  call mpi_bcast(ngL, 1, mpi_integer, 0, mpi_comm_world, ierr)
  call mpi_bcast(nr, 1, mpi_integer, 0, mpi_comm_world, ierr)
  call mpi_bcast(ir_start, 1, mpi_integer, 0, mpi_comm_world, ierr)
  call mpi_bcast(iseed, 1, mpi_integer, 0, mpi_comm_world, ierr)
  !iseed = 2862933555777941757*iseed*rank + 1013904243
  iseed = modulo(53323793, (32674079 / (rank+1))) + iseed
  print *, rank,"  iseed: ",iseed
  call mpi_bcast(outdir, len_filename, mpi_character, 0, mpi_comm_world, ierr)
#endif MPI

  call parameters_init()


  ! ----- Read power spectrum from file
  if (rank .eq. 0) then
     call read_powerspectrum(pk)
     print *, "Finished reading theory power spectrum"
  end if
#ifdef MPI
  call mpi_bcast(pk, size(pk), mpi_double_precision, 0, &
       mpi_comm_world, ierr)


  ! ----- Initialize FFTW plans and MPI slab decomposition
  call init_fftw(rank, numtasks, ng, ngL)
  if (rank .eq. 0) print *, "Finished FFTW init"
#else !MPI
  local_nlast = ng
  local_last_start = 0
  local_nlastL = ngL
  local_last_startL = 0
#endif MPI


  ! ----- Read snapshot particle coordinates
!!$  call sort_snapshot_particles(npart)
!!$  stop
  npart_read = npart / numtasks
  allocate(psnap(3, nsnap, npart_read))
  call input_sorted_snapshot_particles(psnap, 1 + rank*npart_read, npart_read)


  ! ----- Draw new large-scale mode realizations and save power spectra
  call resample_modes(ir_start, nr, ng, ngL, psnap, outdir, pk, iseed)


  ! ----- Clean up
500 deallocate(psnap)
#ifdef MPI
  call dealloc_fftw_plans()
  call mpi_finalize(ierr)
#endif
  stop

end program nsample
