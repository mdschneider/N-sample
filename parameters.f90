!> Simulation parameters to be accessed globally
module parameters
  use types_nsample
  use parameters_nth8

#ifdef MPI
  include "mpif.h"
#endif MPI

  !> MPI task indices and error flag
  integer(i4b) :: rank, numtasks, ierr

#ifdef MPI
	!> Padding of 1st array dimensions to accomodate in-place c2r and r2c FFTs
  integer, parameter :: fftpad = 2 
#else
  integer, parameter :: fftpad = 0 
#endif MPI

  !> Type to hold parameters specific to an N-body simulation run
  type nbody_simulation_params
     real(dp)     :: boxsize, kF, V
     integer(i4b) :: npart
  end type nbody_simulation_params
  type(nbody_simulation_params) :: sim


  !> Binning scheme for computing power spectrum estimates
  !!       0: Grid bins
  !!       1: Linear binning using kbinwidth variable below
  !!          (Use this to match Takahashi et al.: kbinwidth = 0.01)
  integer(i4b), parameter :: pkbintype = 0
  !> Bin limits for output power spectrum
  !!       (Only used if pkbintype == 1)
  !! Width of linear k bins in units of grid spacing
  real(dp), parameter :: kbinwidth = 0.01 ! [def: kF] (h/Mpc) width of linearly spaced k bins

  
#ifdef SCDM
	!> Snapshot times for input simulation
  character(len=250), parameter :: snapshotlist = &
       & "/gpfs/data/bskm46/work/covResampling/input/gadget_scdm_snapshot_list_interpolation.txt"
  !> snapshot number where a=1 (for "smallmodes" simulation)
  integer(i4b), parameter :: isnap_a1 = 34
	!> snapshot number where a=1 (for "allmodes" simulation)
  integer(i4b), parameter :: allmodes_isnap_a1 = 34
!#else ! Takahashi settings - moved to parameters_nth*.f90
#endif SCDM


  ! ----- Theory power spectrum
  !character(len=100), parameter :: pkfile = "/gpfs/data/bskm46/work/covResampling/input/pk_takahashi.dat
  !integer(i4b), parameter       :: Npk = 2000
  character(len=100), parameter :: pkfile = "/gpfs/data/bskm46/work/covResampling/input/takahashi/linear_power_spectrum.txt"
  integer(i4b), parameter       :: Npk = 692
  real(sp), dimension(Npk) :: pkderivs
  real(sp), parameter :: sigma_pkspl = 1.0 ! tension for spline of power spectrum

#ifdef SCDM
	!> Directory where the "allmodes" snapshots are stored.
  character(len=200), parameter :: allmodesdir = &
       &"/gpfs/data/bskm46/work/covResampling/output/scdm/allmodes"
	!> Directory where the "smallmodes" snapshots are stored.
  character(len=200), parameter :: smallmodesdir = &
       &"/gpfs/data/bskm46/work/covResampling/output/scdm/smallmodes8"
!#else ! See parameters_nth*.f90 for other cases
#endif SCDM
  ! List of scale factors where snapshots were stored
  real(dp), dimension(nsnap) :: snapshot_times
  

  ! ----- Table of a versus omega_m at matched growth
#ifdef SCDM
	!> File with table of perturbed scale factor versus Omega_m at matched linear growth.
  character(len=200), parameter :: atable_file="/gpfs/data/bskm46/work/covResampling/input/omega_m_vs_a_table_matchedgrowthICs_SCDM.txt"
	!> File with atable versus Omega_m for at matched linear growth for mode subtraction.
  character(len=200), parameter :: atable_inverse_file="/gpfs/data/bskm46/work/covResampling/input/omega_m_vs_a_table_inverse_SCDM.txt"
#else
  character(len=200), parameter :: atable_file= &
       & "/gpfs/data/bskm46/work/covResampling/input/omega_m_vs_a_table_matchedgrowthICs_20091113.txt"
  character(len=200), parameter :: atable_inverse_file= &
       & "/gpfs/data/bskm46/work/covResampling/input/omega_m_vs_a_table_inverse_takahashi.txt"
#endif SCDM
  !> Number of entries in a(Omega_m) table
  integer(i4b), parameter       :: ntable = 500 
	!> Globally defined table for holding a(Omega_m)
  real(dp), dimension(ntable,2) :: atable
	!> Logical to test whether the a(Omega_m) table has been read from file.
  logical :: atable_read = .false.


  ! ----- Directory to save files for creating dotplots
  character(len=250) :: dotplotdir


  ! ----- Kurtosis computation from resampled density: how many grid locations to use?
  integer(i4b) :: ngsmooth = 4


contains

	!-------------------------------------------------------------
	!> Initialize default parameter values.
	!!
	!! This MUST be called near the start of the program!
	!<------------------------------------------------------------	
  subroutine parameters_init()
    sim%boxsize = 1000.0
    sim%npart   = 256**3
    sim%kF      = TWOPI / sim%boxsize
    sim%V       = boxsize**3
    print *, "Initilized N-body sim parameters",sim%boxsize,sim%npart
  end subroutine parameters_init

end module parameters
