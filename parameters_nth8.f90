!> Parameters specific to nth16 L-Gadget run
module parameters_nth8
  use types_nsample

  integer(i4b), parameter :: nthresh = 8
  integer(i4b), parameter :: nsnap = 43

#ifndef SCDM
  character(len=250), parameter :: snapshotlist = &
       & "/gpfs/data/bskm46/work/covResampling/input/gadget_takahashi_snapshot_list_interpolation_20100103_nth8.txt"

  ! snapshot number where a=1 for smallmodes run
  integer(i4b), parameter :: isnap_a1 = 24

  ! snapshot number where a=1 for allmodes run
  integer(i4b), parameter :: allmodes_isnap_a1 = 15

  character(len=200), parameter :: allmodesdir = &
       & "/gpfs/data/bskm46/work/covResampling/output/takahashi/old_allmodes"

  character(len=200), parameter :: smallmodesdir = &
       & "/gpfs/data/bskm46/work/covResampling/output/takahashi/smallmodes8_20100103"

!!$  character(len=200), parameter :: smallmodesdir = &
!!$       & "/gpfs/data/bskm46/work/covResampling/output/takahashi/lbox125/seed1"

!!$  character(len=200), parameter :: smallmodesdir = &
!!$       & "/gpfs/data/bskm46/gadget_runs/smallmodes8/samp1"
#endif SCDM
  

end module parameters_nth8
