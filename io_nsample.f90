!
!> IO routines for density, mask, and Fourier modes.
!
module io_nsample
  use types_nsample
  use parameters
  implicit none

#ifdef FITPACK
  !--- Interfaces for spline routines from fitpack
  INTERFACE
     SUBROUTINE curv1(n,x,y,slp1,slpm,islpsw,yp,temp,sigma,ierr) ! fitpack.f
       integer n,islpsw,ierr
       real x(n),y(n),slp1,slpn,yp(n),temp(n),sigma
     END SUBROUTINE curv1
     FUNCTION curv2(t,n,x,y,yp,sigma) ! fitpack.f
       integer n
       real t,x(n),y(n),yp(n),sigma,curv2
     END FUNCTION curv2
  END INTERFACE
#endif FITPACK

contains

	!-------------------------------------------------------------
	!> Write the density field to file
	!<------------------------------------------------------------	
  subroutine save_density(outfile,rho,N)
    implicit none
    integer(i4b), intent(in)                   :: N
    character(len=*), intent(in)               :: outfile
    real(dp), dimension(:,:,:), intent(in)     :: rho
    integer(i4b)                               :: ios,unit_out
    character(len=10)                       :: nlab

    call get_free_unit(unit_out)
    open(unit=unit_out, file=outfile, form="unformatted", iostat=ios, &
         & status="replace", action="write", err=510)
!!    write(nlab,*) size(rho)
    write(unit=unit_out, iostat=ios, err=520) rho
    close(unit_out)
    return
    ! ----- IO error handling
510 print *, "Error opening output file "//trim(outfile)
    print *, "ios = ",ios
    close(unit_out)
    stop
520 print *, "Error writing to file "//trim(outfile)
    print *, "ios = ",ios
    close(unit_out)
    stop
  end subroutine save_density

  subroutine save_2Dslice(outfile,rho,N, istart, istop)
    ! ------------------------------------------------------------
    ! Save a 2D slice of rho to a text file
    ! ------------------------------------------------------------
    implicit none
    integer(i4b), intent(in)                   :: N, istart, istop
    character(len=*), intent(in)               :: outfile
    real(dp), dimension(N,N,N), intent(in)     :: rho
    integer(i4b)                               :: ios,unit_out,ii, kk
    character(len=10)                          :: nlab

    call get_free_unit(unit_out)
    open(unit=unit_out, file=outfile, iostat=ios, &
         & status="replace", action="write", err=610)
    write(nlab,*) N
    do kk = istart,istop
       do ii=1,N
          write(unit=unit_out, fmt='('//trim(adjustl(nlab))//'ES18.8)', &
               & iostat=ios, err=620) rho(:,ii,kk)
       end do
    end do
    close(unit_out)
    return
    ! ----- IO error handling
610 print *, "Error opening output file "//trim(outfile)
    print *, "ios = ",ios
    close(unit_out)
    stop
620 print *, "Error writing to file "//trim(outfile)
    print *, "ios = ",ios
    close(unit_out)
    stop
  end subroutine save_2Dslice

  subroutine read_powerspectrum(pk)
    ! ------------------------------------------------------------
    ! Read the power spectrum for the conditional mode distribution
    ! from the file pkfile
    ! ------------------------------------------------------------
    implicit none
    real(dp), dimension(Npk,2), intent(out) :: pk
    integer(i4b)                        :: ii, unit_in,ios,N
    real(dp) :: k,Delta,junk1,junk2, pkval
    
    print *, "Reading pk model from ",trim(pkfile)

    call get_free_unit(unit_in)
    open(unit=unit_in,file=pkfile, form="formatted", iostat=ios, &
         & status="old", action="read", err=610)
!!    read(unit_in,*) N,junk1
    do ii =1,size(pk(:,1))
       read(unit_in,*) k,pkval !Delta,junk1,junk2
#ifdef PK_LOGK
       pk(ii,1) = log(k)
#else
       pk(ii,1) = k
#endif
       pk(ii,2) = log(pkval) !Delta*2.d0*PI**2/k**3
    end do
    close(unit_in)
    return
610 print *, "Error opening power spectrum file "//trim(pkfile)
    print *, "ios = ",ios
    close(unit_in)
    return
  end subroutine read_powerspectrum

  subroutine init_powerspectrum_spline(pk)
    ! -------------------------------------------------------------
    ! Initialize a spline of the log of the power spectrum pk
    ! using fitpack.
    ! -------------------------------------------------------------
    real(dp), dimension(Npk,2), intent(in) :: pk
    integer(i4b) :: ierr
    real(sp), dimension(Npk) :: temp
    call curv1(Npk, real(pk(:,1),sp), real(pk(:,2),sp), 0.0, 0.0, 3, pkderivs, temp, sigma_pkspl, ierr)
    if (ierr /= 0) then
       print *, "ERROR: power spectrum spline initialization failed, ierr:",ierr
       stop
    end if
  end subroutine init_powerspectrum_spline

	!-------------------------------------------------------------
	!> Save the sample power spectra and mean sample power spectrum
  !! to text files.
	!<------------------------------------------------------------	
  subroutine save_powerspectra(N,Nr,meanpk,samppk)
    integer(i4b), intent(in)                  :: N,Nr
    real(dp), dimension(0:N/2), intent(in)    :: meanpk
    real(dp), dimension(Nr,0:N/2), intent(in) :: sampPk
    integer(i4b)                              :: unit_out,a,ios
    character(len=200)                        :: nent,outfile

    ! ----- save mean Pk to file
    outfile = "meanPk.txt"
    call get_free_unit(unit_out)
    open(unit_out,file=outfile, status="replace", action="write",&
         & iostat=ios, err=510)
    write(nent,*) N/2+1
    write(unit_out,'('//trim(adjustl(nent))//'ES18.8)', &
         & iostat=ios, err=520) meanpk
    close(unit_out)
    ! ----- save sample Pk's to file
    outfile = "samplePk.txt"
    call get_free_unit(unit_out)
    open(unit_out,file=outfile, status="replace", action="write", &
         & iostat=ios, err=510)
    do a=1,Nr
       write(unit_out,'('//trim(adjustl(nent))//'ES18.8)', &
            & iostat=ios, err=520) sampPk(a,:)
    end do
    close(unit_out)
    return
    ! ----- IO error handling
510 print *, "Error opening output file "//trim(outfile)
    print *, "ios = ",ios
    close(unit_out)
    return
520 print *, "Error writing to file "//trim(outfile)
    print *, "ios = ",ios
    close(unit_out)
    return
  end subroutine save_powerspectra


  subroutine save_dotplot_particles(outdir, filename, slabwidth, coor, axis)
    ! ------------------------------------------------------------
    ! Save particle positions in a slab of width slabwidth along 
    ! specified axis for later dotplotting.
    ! ------------------------------------------------------------
    character(len=*), intent(in) :: outdir, filename
    real(dp), intent(in) :: slabwidth
    real(dp), dimension(3,npart), intent(in) :: coor
    integer(i4b), intent(in) :: axis
    character(len=250) :: outfile
    integer(i4b) :: unit_out, i
    print *, "Saving dotplot particles for ",trim(filename)
    call get_free_unit(unit_out)
    outfile = trim(outdir)//"/"//trim(filename)
    open(unit_out, file=outfile, &
         & status="replace", action="write")
    do i = 1,npart
       if (coor(axis,i) .gt. 0.0 .and. coor(axis,i) .lt. slabwidth) then
          write(unit_out, fmt='(3ES18.8)') coor(:,i)
       end if
    end do
    close(unit_out)
  end subroutine save_dotplot_particles


	!-------------------------------------------------------------
	!> Read table of scale factor versus Omega_m at fixed linear
  !! growth function from file.
	!!
  !! Columns: Omega_m scale_factor
	!<------------------------------------------------------------	
  subroutine read_a_omega_table()
    use parameters
    integer :: unit_in, ios, i
    atable_read = .true.
    if (rank .eq. 0) then
       call get_free_unit(unit_in)
       open(unit=unit_in, file=atable_file, form="formatted", iostat=ios, &
            & status="old", action="read", err = 210)
       do i = 1,ntable
          ! Omega_m, a
          read(unit_in, *, err=220) atable(i,:)
       end do
       close(unit_in)
    end if
#ifdef MPI
    call mpi_bcast(atable, size(atable), mpi_double_precision, 0, &
         & mpi_comm_world, ierr)
#endif MPI
    return
210 print *, "Error opening input file "//trim(atable_file)
    print *, "ios = ",ios
    close(unit_in)
    stop
220 print *, "Error reading from file "//trim(atable_file)
    print *, "ios = ",ios
    close(unit_in)
    stop
  end subroutine read_a_omega_table

	!-------------------------------------------------------------
	!> Read table of scale factor versus Omega_m at matched inverse
  !! linear growth function from file.
	!!
  !! Columns: Omega_m scale_factor
	!<------------------------------------------------------------	
  subroutine read_a_omega_table_inverse()
    use parameters
    integer :: unit_in, ios, i
    atable_read = .true.
    if (rank .eq. 0) then
       call get_free_unit(unit_in)
       open(unit=unit_in, file=atable_inverse_file, form="formatted", iostat=ios, &
            & status="old", action="read", err = 230)
       do i = 1,ntable
          ! Omega_m, a
          read(unit_in, *, err=240) atable(i,:)
       end do
       close(unit_in)
    end if
#ifdef MPI
    call mpi_bcast(atable, size(atable), mpi_double_precision, 0, &
         & mpi_comm_world, ierr)
#endif MPI
    return
230 print *, "Error opening input file "//trim(atable_file)
    print *, "ios = ",ios
    close(unit_in)
    stop
240 print *, "Error reading from file "//trim(atable_file)
    print *, "ios = ",ios
    close(unit_in)
    stop
  end subroutine read_a_omega_table_inverse


	!-------------------------------------------------------------
	!> Read particle coordinates from snapshot files.
	!<------------------------------------------------------------	
  subroutine read_snapshot(snapdir, isnap, coor, idsout, asnap, icfile, insim)
    use read_gadgetmod
    character(len=*), intent(in)              :: snapdir
    integer, intent(in)                       :: isnap
    real(dp), dimension(:,:), intent(inout)   :: coor
    integer*8, dimension(:), intent(out)      :: idsout
    real(dp), intent(out)                     :: asnap
    logical, optional                         :: icfile
    type(nbody_simulation_params), optional   :: insim

    real, dimension(:), allocatable           :: x,y,z,vx,vy,vz
    integer*8, dimension(:), allocatable      :: ids
    integer(i4b), dimension(:), allocatable   :: pid
    type(nbody_simulation_params)             :: simloc
    type(headertype)                          :: header    
    character(len=250)                        :: infile
    character(len=3)                          :: slab
    integer                                   :: i,j, np
    logical :: velocities

    velocities = .false.
    if (size(coor(:,1)) .eq. 6) velocities = .true.

    if (present(insim)) then
       simloc = insim
    else
       simloc = sim
    end if
    print *, "read_snapshot: Assigned simloc",simloc%boxsize,simloc%npart

    write(slab,'(i3.3)') isnap
    if (present(icfile) .and. icfile) then
       infile = trim(snapdir)//"/ics"
    else
       infile = trim(snapdir)//"/snapdir_"//trim(adjustl(slab))//&
            &"/snapshot_"//trim(adjustl(slab))
    end if
    call read_header(infile, header, silent=.true.)
    asnap = 1.0/(1.0+header%redshift)
    print *, "snapshot scale factor:",asnap
    np = simloc%npart
    allocate(x(np), y(np), z(np), ids(np))
    if (size(coor(:,1)) .eq. 6) then
       velocities = .true.
       allocate(vx(np), vy(np), vz(np))
       call read_gadget(infile,2,np,x,y,z, vx,vy,vz, id8=ids, silent=.true.)
    else
       call read_gadget(infile,2,np,x,y,z, id8=ids, silent=.true.)
    end if
    print *, "Sorting particles "
    allocate(pid(simloc%npart))
!!    call indexxxi8(simloc%npart, ids, pid)
!!    do j = 1,simloc%npart
!!       idsout(j) = ids(pid(j))
!!       coor(1,j) = x(pid(j))
!!       coor(2,j) = y(pid(j))
!!       coor(3,j) = z(pid(j))
!!       if (velocities) then
!!          coor(4,j) = vx(pid(j))
!!          coor(5,j) = vy(pid(j))
!!          coor(6,j) = vz(pid(j))
!!       end if
!!    end do
    coor(1,:) = x
    coor(2,:) = y
    coor(3,:) = z
    if (velocities) then
       coor(4,:) = vx
       coor(5,:) = vy
       coor(6,:) = vz
    end if
    where(coor(1:3,:) .lt. 0.0) coor(1:3,:) = coor(1:3,:) + boxsize
    where(coor(1:3,:) .ge. boxsize) coor(1:3,:) = coor(1:3,:) - boxsize
    
!!    deallocate(x,y,z,ids,pid)
    deallocate(x,y,z, ids)
    if (velocities) deallocate(vx,vy,vz)
  end subroutine read_snapshot

	!-------------------------------------------------------------
	!> Create the filename for saving a power spectrum
	!<------------------------------------------------------------	
  subroutine psfilename(psfile, filetag, psfilepath, ng)
    character(len=len_filename), intent(out)  :: psfile
    character(len=*), intent(in)     :: filetag, psfilepath
    integer, intent(in)              :: ng
    character(len=20)                :: cng
    write(cng,*) ng
#ifdef CIC
    psfile = trim(psfilepath)//"/ps_"//trim(adjustl(filetag))//&
         &"_ng"//trim(adjustl(cng))//&
         &"_CIC.txt"
#else
    psfile = trim(psfilepath)//"/ps_"//trim(adjustl(filetag))//&
         &"_ng"//trim(adjustl(cng))//&
         &".txt"
#endif CIC
  end subroutine psfilename

	!-------------------------------------------------------------
	!> Read particle coordinates grouped so that the positions of 
  ! each particle for all saved snapshots are together.
	!<------------------------------------------------------------	
  subroutine input_sorted_snapshot_particles(coor, ipart_start, npart_read)
    use hdf5_wrapper
    integer(i4b), intent(in)    :: ipart_start,npart_read
    real(dp), dimension(3,nsnap,npart_read), intent(inout) :: coor
    integer(i8b), dimension(nsnap,npart_read) :: ids
    integer(i4b), dimension(3)       :: start, snapcount
    integer(i4b), dimension(2)       :: idstart, idcount
    integer(i4b)                :: ifile, ir, i,j
    character(len=250)          :: infile
    real(dp), dimension(nsnap)  :: scale_factors
    real :: t1,t2
    start = (/1_i4b, 1_i4b, ipart_start/)
    snapcount = (/3_i4b, nsnap, npart_read/)
    idstart = (/1_i4b, ipart_start/)
    idcount = (/nsnap, npart_read/)
    infile = trim(smallmodesdir)//"/sorted_snapshots.h5"
    do ir = 0, numtasks-1
       if (ir .eq. rank) then
          print *, rank,"  Reading sorted snapshot particles from file",trim(infile)
          call hdf5_open_file(ifile, infile, readonly=.true.)
          call cpu_time(t1)
          call hdf5_read_data(ifile, "coordinates", coor, start, snapcount)
          call hdf5_read_data(ifile, "ids", ids, idstart, idcount)
          call cpu_time(t2)
          print *, "  HDF5 read time: ",t2-t1," seconds"
          call hdf5_read_attribute(ifile, "coordinates/snapshot_times", snapshot_times)
          call hdf5_close_file(ifile)
#ifdef MPI
          call mpi_barrier(mpi_comm_world, ierr)
#endif MPI
       end if
    end do

    ! ----- Check that ids match for all snapshots
    do i = 1,npart_read
       do j = 1,nsnap-1
          if (ids(j,i) /= ids(j+1,i)) then
             print *, rank,"  ERROR: ids do not match ",j,j+1
          end if
       end do
    end do
#ifdef VERBOSE
    print *, rank,"  Finished reading snapshot particles"
#endif
  end subroutine input_sorted_snapshot_particles

	!-------------------------------------------------------------
	!> Sort particle coordinates from all snapshots into groups
  !! with the positions of each particle in all snapshots together.
  !!
  !! If npart_file < npart (in parameters.f90) and 
  !! npart / npart_file = integer, then the particles will be 
  !! replicated npart/npart_file times in each dimension to fill
  !! a larger box.
  !! If npart_file > npart, then abort.
  !!
  !! This routine has no MPI support.
	!! @param npart_file Number of particles per file
	!<------------------------------------------------------------	
  subroutine sort_snapshot_particles(npart_file)
    use hdf5_wrapper
    integer(i4b), dimension(3), parameter :: datasize=(/3,nsnap,npart/)
    integer(i4b), dimension(2), parameter :: idsize=(/nsnap,npart/)
    integer(i4b), parameter  :: nread = 2 ! Read this many snapshot files at once

    integer(i4b), intent(in) :: npart_file
    integer(i4b) :: j,j1,j2,j3, nrepl, nrepl_perdim
    real(dp) :: lbox_file

    integer(i4b), dimension(3)                :: start, snapcount
    integer(i4b), dimension(2)                :: idstart, idcount
    double precision, dimension(3,nread,npart)   :: coor
    double precision, dimension(3,npart_file)    :: coor_buf
    integer*8, dimension(npart_file)  :: ids
    integer*8, dimension(nread,npart) :: snapids
    real(dp)                       :: asnap
    real(dp), dimension(nsnap)     :: scale_factors
    integer(i4b)                   :: isnap, ifile, i
    character(len=250)             :: outfile
    real(sp) :: t1,t2

    type(nbody_simulation_params) :: insim

    print *, "Sorting snapshot particles..."
    print *, "npart_file:",npart_file,"  npart:",npart

    ! ----- Error checking
    if (npart_file .gt. npart) then
       print *, "ERROR: npart_file > npart",npart_file,npart
       stop
    end if
    if (npart_file .lt. npart) then 
       nrepl = npart / npart_file
       nrepl_perdim = int(real(nrepl)**(1.0/3.0))
       lbox_file = boxsize / nrepl_perdim
       insim%npart = npart_file
       insim%boxsize = lbox_file
       insim%V = insim%boxsize**3
       insim%kF = TWOPI / insim%boxsize
    else
       nrepl = 1
       insim = sim
    end if
    print *, "  nrepl:",nrepl,"  nrepl_perdim:",nrepl_perdim
    print *, "  lbox_file:",lbox_file


    print *, "Sorting snapshot particles from "//trim(smallmodesdir)
    ! ----- Initialize
    start = (/1, 1, 1/)
    snapcount = (/3, nread, npart/)
    idstart = (/1, 1/)
    idcount = (/nread, npart/)
    outfile = trim(smallmodesdir)//"/sorted_snapshots.h5"
    print *, "Saving sorting particle list to "//trim(outfile)
    call hdf5_create_file(ifile, outfile)
    ! ----- Loop over snapshots to read
    do isnap = 0,(nsnap-1)/nread
       ! --- Read snapshots
       call read_snapshot(smallmodesdir, nread*isnap, coor_buf, ids, asnap, &
            icfile=.false., insim = insim)
       print *, "Range of coor_buf:",minval(coor_buf),maxval(coor_buf)
       call assign_file_coords(coor, coor_buf, snapids, ids, &
            & nread, npart_file, 1, &
            & nrepl, nrepl_perdim, lbox_file)
       print *, "Range of assigned coordinates:",minval(coor),maxval(coor)
       scale_factors(nread*isnap+1) = asnap
       ! check if we've reached the end of the snapshots
       if (nread*isnap+2 .le. nsnap) then 
          call read_snapshot(smallmodesdir, nread*isnap+1, coor_buf, ids, asnap, &
               icfile=.false., insim = insim)
          print *, "Range of coor_buf:",minval(coor_buf),maxval(coor_buf)
          call assign_file_coords(coor, coor_buf, snapids, ids, &
               & nread, npart_file, 2, &
               & nrepl, nrepl_perdim, lbox_file)
          print *, "Range of assigned coordinates:",minval(coor),maxval(coor)
          scale_factors(nread*isnap+2) = asnap
       else
          snapcount(2) = 1
          idcount(1) = 1
       end if
       ! --- Update indices into total data array(s) in output file
       start(2)   = nread*isnap+1
       idstart(1) = nread*isnap+1
       ! --- Write to snapshots
       print *, "Writing to hdf5 file..."
       print *, start,snapcount
       print *, idstart,idcount
       call cpu_time(t1)
       if (nread*isnap+2 .le. nsnap) then
          call hdf5_write_data(ifile, "coordinates", coor, &
               & initial_size=datasize, start=start, count=snapcount, &
               & gzip = 2, overwrite=.false., extensible=.false.)
          call hdf5_write_data(ifile, "ids", snapids, &
               & initial_size=idsize, start=idstart, count=idcount, &
               & gzip = 2, overwrite=.false., extensible=.false.)
       else
          call hdf5_write_data(ifile, "coordinates", coor(:,1,:), &
               & initial_size=datasize, start=start, count=snapcount, &
               & gzip = 2, overwrite=.false., extensible=.false.)
          call hdf5_write_data(ifile, "ids", snapids(1,:), &
               & initial_size=idsize, start=idstart, count=idcount, &
               & gzip = 2, overwrite=.false., extensible=.false.)
       end if
       call cpu_time(t2)
       print *, "  HDF5 write time: ",t2-t1," seconds"
    end do
    ! ----- Save scale factor times at snapshots as an HDF5 attribute
    call hdf5_write_attribute(ifile, "coordinates/snapshot_times", scale_factors)
    call hdf5_close_file(ifile)
  end subroutine sort_snapshot_particles

!!  subroutine assign_file_coords(coor, coor_buf, snapids, ids, &
!!       & nread, npart_file, icoor, &
!!       & nrepl, nrepl_perdim, lbox_file)
!!    ! ------------------------------------------------------------
!!    ! Replicate coordinates and assign to output array
!!    ! ------------------------------------------------------------
!!    real(dp), dimension(3,nread,npart), intent(inout) :: coor
!!    real(dp), dimension(3,npart_file), intent(in)     :: coor_buf
!!    integer(i8b), dimension(nread,npart), intent(inout) :: snapids
!!    integer(i8b), dimension(npart_file), intent(in)     :: ids
!!    integer(i4b), intent(in) :: nread, npart_file, icoor, nrepl, nrepl_perdim
!!    real(dp), intent(in) :: lbox_file
!!    integer(i4b) :: i,j,j1,j2,j3, imin,imax
!!
!!    print *, "assign_file_coords - npart_file:",npart_file,"  nrepl_perdim:",nrepl_perdim,"  lbox_file:",lbox_file
!!    if (nrepl .gt. 1) then
!!       do j1 = 1,nrepl_perdim
!!          do j2 = 1,nrepl_perdim
!!             do j3 = 1,nrepl_perdim
!!                j = j3 + nrepl_perdim * (j2 - 1 + nrepl_perdim*(j1-1))
!!                imin = (j-1) * npart_file + 1
!!                imax = j * npart_file
!!                coor(1,icoor,imin:imax) = coor_buf(1,:) + (j1-1)*lbox_file
!!                coor(2,icoor,imin:imax) = coor_buf(2,:) + (j2-1)*lbox_file
!!                coor(3,icoor,imin:imax) = coor_buf(3,:) + (j3-1)*lbox_file
!!                snapids(icoor, imin:imax) = ids
!!             end do
!!          end do
!!       end do
!!    else
!!       forall (i=1:npart) coor(:,icoor,i) = coor_buf(:,i)
!!    end if
!!  end subroutine assign_file_coords
!!
!!
!!  subroutine saveKurtosisListFromResampledDensity(kurtosisList)
!!    ! ------------------------------------------------------------
!!    !
!!    ! ------------------------------------------------------------
!!    real(dp), dimension(:,0:), intent(in) :: kurtosisList
!!    integer(i4b) :: unit_out, ios, i
!!    character(len=200) :: nent, outfile, ranklab
!!    write(ranklab,*) rank
!!    outfile = "/data/rw7/bskm46/work/covResampling/output/nsamp/test/kurtosis/kurtosis_resamp_"//trim(adjustl(ranklab))//".txt"
!!    call get_free_unit(unit_out) 
!!    open(unit_out,file=outfile, status="replace", action="write", &
!!         & iostat=ios, err=840)
!!    write(nent,*) size(kurtosisList(:,1))
!!    do i=0,size(kurtosisList(1,:))-1
!!       write(unit_out,'('//trim(adjustl(nent))//'ES18.8)', &
!!            & iostat=ios, err=850) kurtosisList(:,i)
!!    end do
!!    close(unit_out)
!!    return
!!    ! ----- IO error handling
!!840 print *, "Error opening output file "//trim(outfile)
!!    print *, "ios = ",ios
!!    close(unit_out)
!!    return
!!850 print *, "Error writing to file "//trim(outfile)
!!    print *, "ios = ",ios
!!    close(unit_out)
!!    return
!!  end subroutine saveKurtosisListFromResampledDensity

end module io_nsample
