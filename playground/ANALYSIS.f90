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

!Set .true. to generate top-level heat maps for each grid
logical,parameter :: heatmap_flag = .true.

!Set .true. to generate the scattering angle plots of the trajectories
!that were generated with molecular dynamics for each grid
logical,parameter :: trueSA_flag = .false.

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

   !Set how many trajectories will be generated for the test
   !If old data is being used, this number will be decreased internally
   integer,parameter :: Ntesttraj = 1000

   !Set .true. to generate the RMSD frequency plots for each
   !trajectory tested for each grid
   !(Not recommended for large Ntesttraj or Ngrid_total_cap)
   logical,parameter :: testtrajRMSD_flag = .false.

   !Set .true. to generate a frequency plot of the percentage
   !of frames in each trajectory that were below some threshold RMSD
   !One plot will be produced per grid
   logical,parameter :: percentthreshold_flag = .true.

      !Set the threshold RMSD to be used for any rejection method
      !real(dp),parameter :: !threshold_rmsd! = !.200100d0
      real(dp),parameter :: threshold_rmsd = .00100d0

      !Set .true. to generate trajectories using md-calculated gradients
      !Otherwise, the program will use the above threshold as a rejection
      !method
      logical,parameter :: reject_flag = .true.

   !Set .true. to generate the scattering angle plots of
   !the trajectories for each grid
   logical,parameter :: testtrajSA_flag = .true.

!Set .true. to generate trajectories with gradient approximations
!produced from the grid (no MD at all)
logical,parameter :: approxtraj_flag = .false.

   !Set how many trajectories will be generated for the test
   integer,parameter :: Napproxtraj = 100

end module ANALYSIS



