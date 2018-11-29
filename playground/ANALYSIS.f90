module ANALYSIS
use PHYSICS
implicit none

!This module is inteded to be changed according to what the user
!wants to analyze from a grid that was made.
!Analysis is done in roughly the order they appear here

!Set the number of grids to be analyzed; will start at 001 and increment
!If this number is larger than the number of grids in the folder,
!then it will default to the number of grids in the folder
integer,parameter :: Ngrid_cap = 1
!The number of grids we will end up using (never more than Ngrid_cap)
integer :: Ngrid_total

!Set the number of children cells to be checked
integer,parameter :: Norder_cap = 1
!The number of children cells we will end up check (never more than Ngrid_cap)
integer :: Norder_total

!Set .true. to generate top-level heat maps for each complete grid
logical,parameter :: heatmap_flag = .true.

!Set .true. to generate the scattering angle plots of the trajectories
!that were generated with molecular dynamics for each grid
logical,parameter :: trueSA_flag = .false.
logical,parameter :: trueED_flag = .false.
integer,parameter :: SA_Nbins = 50
integer,parameter :: TRV_Nbins = 50

!Set .true. to generate trajectories with molecular dynamics
!and test them against the grid
!Note: if you ALREADY made these trajectories (stored in I0.6 files)
!then you can set this .false. to proceed with analysis and save computation
logical,parameter :: testtraj_flag = .true.

   !Set .true. if you are generating new trajectories and you want them
   !to be stored in new files (don't overwrite old files)
   !Only set .true. if, essentially, you are alright using old
   !data and want to save computation
   logical,parameter :: useolddata_flag = .false.

   !Set .true. if you want to use the initial conditions of a
   !previous test; specify which one in the variable that follows
   logical,parameter :: useoldinitialbonddata_flag = .false.
   integer,parameter :: initialbondname_length = 40
   character(initialbondname_length),parameter :: initialbondname = ""

   !Set how many trajectories will be generated for the test
   !If old data is being used, this number will be decreased internally
   integer,parameter :: Ntesttraj = 1000

   !Set .true. to generate the RMSD frequency plots for each
   !trajectory tested for each grid
   !(Not recommended for large Ntesttraj or Ngrid_total_cap)
   logical,parameter :: testtrajRMSD_flag = .false.

   !Set .true. to generate the checkTrajectory plots for each
   !trajectory tested for each grid
   !(Not recommended for large Ntesttraj or Ngrid_total_cap)
   logical,parameter :: testtrajDetailedRMSD_flag = .true.

   !Set .true. to generate a frequency plot of the percentage
   !of frames in each trajectory that were below some threshold RMSD
   !One plot will be produced per grid
   logical,parameter :: percentthreshold_flag = .true.
   integer,parameter :: RMSD_Nbins = 50

      !Set the threshold RMSD to be used for any rejection method
      !real(dp),parameter :: !threshold_rmsd! = !.200100d0
      real(dp),parameter :: threshold_rmsd = .00100d0
      real(dp),parameter :: default_rmsd = 1.000100d0

      !Set .true. to generate trajectories using md-calculated gradients
      !Otherwise, the program will use the above threshold as a rejection
      !method
      logical :: reject_flag = .true.

      !If reject_flag is false (and we are accepting frames) then
      !accept_first controls whether we use the first frame accepted
      !or use the frame that is closest in RMSD
      logical :: accept_first = .true.

      !If reject_flag is false (and we are accepting frames) then
      !accept_worst indicates to accept the frame with the
      !maximum RMSD (instead of the minimum RMSD)
      logical :: accept_worst = .false.

   !Set .true. to generate the scattering angle plots of
   !the trajectories for each grid
   logical,parameter :: testtrajSA_flag = .true.

   !Set .true. to generate the scattering angle plots of
   !the trajectories for each grid
   logical,parameter :: testheatmapSA_flag = .true.

   !Set .true. to generate the energy decomposition plots of
   !the trajectories for each grid
   logical,parameter :: testtrajTRV_flag = .true.

!Set .true. to ask the program to do a comparison of multiple
!trajectory sets (mostly for consistency checking)
logical,parameter :: comparison_flag = .false.
character(40),parameter :: comparison_SATRVname = ""
integer :: comparison_SATRVcolumn
real,parameter :: comparison_lowerlimit = 0.0d0
real,parameter :: comparison_upperlimit = 0.0d0

integer,parameter :: comparison_number = 1
character(11),parameter :: allprefixes = "placeholder"
integer,dimension(comparison_number),parameter :: alllengths = (/ 69 /)

end module ANALYSIS



