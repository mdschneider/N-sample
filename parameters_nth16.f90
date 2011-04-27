!> Parameters specific to nth16 L-Gadget run
module parameters_nth16
  use types_nsample

  integer(i4b), parameter :: nthresh = 16
  integer(i4b), parameter :: nsnap = 64

  character(len=250), parameter :: snapshotlist = &
       & "/data/rw7/bskm46/work/covResampling/input/gadget_takahashi_snapshot_list_interpolation_20100322_nth16.txt"

  ! snapshot number where a=1 for smallmodes run
  integer(i4b), parameter :: isnap_a1 = 36

  ! snapshot number where a=1 for allmodes run
  integer(i4b), parameter :: allmodes_isnap_a1 = 15

  character(len=200), parameter :: allmodesdir = &
       & "/data/rw7/bskm46/work/covResampling/output/takahashi/old_allmodes"

  character(len=200), parameter :: smallmodesdir = &
       & "/data/rw7/bskm46/work/covResampling/output/takahashi/smallmodes16_20100322"

end module parameters_nth16
