module ANALYSIS
use PHYSICS
implicit none

!This module is inteded to be changed according to what the user
!wants to analyze from a grid that was made.
!Analysis is done in roughly the order they appear here

!The name of this experiment
integer,parameter :: expfolder_length = 7
character(expfolder_length),parameter :: expfolder = "exp001/"

!Set the number of grids to be analyzed; will start at 001 and increment
!If this number is larger than the number of grids in the folder,
!then it will default to the number of grids in the folder
integer,parameter :: Ngrid_cap = 1
!The number of grids we will end up using (never more than Ngrid_cap)
integer :: Ngrid_total
!$OMP THREADPRIVATE(Ngrid_total)

!The number of threads to use
integer, parameter :: Nthreads = 1

!Set the number of children cells to be checked
integer,parameter :: Norder_cap = 1

!Set .true. to generate top-level heat maps for each complete grid
logical,parameter :: heatmap_flag = .true.
character(14),parameter :: populationfile = "population.dat"
character(12),parameter :: covarmapfile = "covarmap.dat"

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
   integer,parameter :: initialbondfolder_length = 40
   character(initialbondfolder_length),parameter :: initialbondfolder = ""

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
   logical :: testtrajDetailedRMSD_flag = .false.


       !Various variables if we are tracking the interpolation
       !at various alpha ratios for this detailed trajectory
       integer,parameter :: Nalpha = 61
       real(dp) :: alpha_start = -7.0d0
       real(dp) :: alpha_end = 5.0d0
       logical :: logarithmic_alpha_flag = .true.
       character(9) :: alphafile = "alpha.dat"

   !This takes much more time but you can force the checkState subroutine
   !to also check the rmsd of frames in adjacent cells
   logical :: force_Neighbors = .true.

   logical :: force_Permutations = .false.

   logical :: memory_flag = .true.
   integer :: single_index_max

   !The point at which checkState stops checking neighbors
   !is determined here
   !The default (if none of these are set) is all zeros
   integer,parameter :: ssm_length = 2
   integer,dimension(ssm_length) :: ssm1 = (/ 24, 00 /)
   integer,dimension(ssm_length) :: ssm2 = (/ 00, 00 /)

   !Set .true. to generate a frequency plot of the percentage
   !of frames in each trajectory that were below some threshold RMSD
   logical,parameter :: percentthreshold_flag = .true.
   integer,parameter :: RMSD_Nbins = 50

   !One plot will be produced per grid if percentthrshold_key = 0
   !Otherwise, percentthreshold_key is translated into binary
   !The string of 1s and 0s dictate which number of grids will be
   !plotted and which won't
   ! ex. 27 = 11011 = 1, 2, 4, and 5 grids will be plotted
   integer,parameter :: percentthreshold_key = 0

      !Set the threshold RMSD to be used for any rejection method
      !real(dp),parameter :: !threshold_rmsd! = !.200100d0
      real(dp) :: outer_threshold_SI = .00100d0
      real(dp) :: inner_threshold_SI = 0.0d0
      real(dp),parameter :: R1_threshold_SI = 0.10d0
      real(dp),parameter :: R2_threshold_SI = 0.10d0
      integer :: Nsort = 1

      !The flag for diversity
      logical :: diversity_flag = .true.

      !Set .true. to generate trajectories using md-calculated gradients
      !Otherwise, the program will use the above threshold as a rejection
      !method
      logical :: reject_flag = .true.
      logical :: readtrajectory_flag = .false.

      !If reject_flag is false (and we are accepting frames) then
      !accept_first controls whether we use the first frame accepted
      !or use the frame that is closest in RMSD
      logical :: accept_first = .true.

      !If reject_flag is false (and we are accepting frames) then
      !accept_worst indicates to accept the frame with the
      !maximum RMSD (instead of the minimum RMSD)
      logical :: accept_worst = .false.

      !Set .true. if the real force calculations we do should be
      !added to the grid
      integer :: grid_addition = 0

      !Set .true. if you are continuing a previous analysis
      logical,parameter :: continue_analysis = .false.

      !$OMP THREADPRIVATE(reject_flag,accept_first,accept_worst)

         !Set .true. if interpolation should be used; that is to say
         !a weighted combination of acceptable frames are used to
         !calculate an approximate gradient
         logical :: interpolation_flag = .true.

         !Interpolation requires a scaling parameter for the weights
         !This is a positive, nonzero real number
         !(Now deprecated from introduction of 2nd order LS minimization)
         real(dp) :: interpolation_alpha1 = 2.0d0

         !For persistent data collection and analysis, a new file
         !is dedicated to interpolation results
         !Whether we record or not to this file is governed by the
         !gather_interpolation_flag
         logical :: gather_interpolation_flag = .true.
         character(17),parameter :: interpolationfile = "interpolation.dat"
         character(14),parameter :: interpolationfolder = "interpolation/"

            !Interpolation data can be checked at any point
            !whenever the interpolation counter reaches the
            !interpolation check, the data is checked
            !If the visual flag is true, then a visual is also
            !made when the data is checked
            integer :: interpolation_check = 50000
            integer :: interpolation_counter
            logical :: interpolation_check_visual = .false.

            real(dp) :: alpha_ratio = 1.0d0

            integer,parameter :: Ninterpolation_max = 40
            integer,parameter :: Ninterpolation_cap = 20

            integer,parameter :: Naccept_max = 20
            real(dp),allocatable :: trajRMSDbuffer(:,:)

            integer,parameter :: Nalpha_tries_max = 3
            real(dp),dimension(Nalpha_tries_max),parameter :: alpha_ratio_list = &
                    (/ 10.0d-2, 10.0d-1, 10.0d0 /)
!                   (/ 10.0d0, 10.0d-1, 10.0d0 /)

            !We can also make what I like to call
            !"consolidated" dropoffs that do an interpolation
            !for every N points up until N*M
            logical,parameter :: dropoff_flag = .false.
            character(11),parameter :: dropofffile = "dropoff.dat"
            integer,parameter :: dropoff_Npacket = 4
            integer,parameter :: dropoff_Mpacket = 5
            integer,parameter :: dropoff_NM = dropoff_Npacket*dropoff_Mpacket

            !For errorCheck11
            logical :: errorCheck11_flag = .false.
            character(8),parameter :: errorCheck11file = "tcut.dat"
            integer,parameter :: Ntcut = 6
            real(dp),parameter :: inner_start = 0.0000d0
            real(dp),parameter :: outer_start = 0.1000d0
            real(dp),parameter :: tcut_interval = 0.1000d0

            !For FCCheck
            logical :: FCCheck_flag = .false.

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
character(69),parameter :: comparison_file = ""
character(40),parameter :: comparison_SATRVname = ""
integer :: comparison_SATRVcolumn
real,parameter :: comparison_lowerlimit = 0.0d0
real,parameter :: comparison_upperlimit = 0.0d0

!Some variables used to calculate the KRP
!Useful to be global so that the convergence can
!communicate with the comparison
real :: lambda_penalty
real :: minsd_penalty

integer,parameter :: comparison_number = 1
character(expfolder_length),parameter :: allprefixes = expfolder
integer,dimension(comparison_number),parameter :: alllengths = (/ expfolder_length /)

character(11),parameter :: analysisfile = "placeholder"


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               getPrefixText
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!       PURPOSE
!               This subroutine uniquely describes an approximation method with a string
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!               prefix_text                     CHAR(12)                        The string that uniquely describes
!                                                                               the approximation method
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               reject_flag                     LOGICAL                         If true, then reject all approximations;
!                                                                               otherwise, use approximations
!               accept_first                    LOGICAL                         If true, then approximations use the first candidate;
!                                                                               otherwise, they use the best candidate
!               accept_worst                    LOGICAL                         If true, then approximations use the worst candidate;
!                                                                               otherwise, they use the best candidate
!
!               threshold_rmsd                  REAL(DP)                        The RMSD threshold which differentiates candidates
!                                                                               for approximation from others
!
!               reject_text                     CHAR(6)                         A string desribing the approximation method
!               Nthreshold_text                 CHAR(6)                         A string desribing the threshold of acceptance
!                                                                               for the approximation method

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getPrefixText(prefix_text)
use PARAMETERS
implicit none

!Grid Directory/File Formatting Strings
character(6) :: reject_text
character(6) :: Nthreshold_text
character(12),intent(out) :: prefix_text

write(Nthreshold_text,FMT=FMT6_pos_real0) outer_threshold_SI
if (reject_flag) then
        reject_text = "reject"
else
        if (accept_first) then
                 if (accept_worst) then
                         reject_text = "alphaW"
                 else
                         reject_text = "alphaA"
                 end if
        else
                 if (accept_worst) then
                         reject_text = "omegaW"
                 else
                         reject_text = "omegaA"
                 end if
        end if
end if

prefix_text = reject_text//Nthreshold_text

end subroutine getPrefixText


end module ANALYSIS



