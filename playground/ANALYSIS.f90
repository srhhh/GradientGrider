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
!$OMP THREADPRIVATE(Ngrid_total)

!The number of threads to use
integer, parameter :: Nthreads = 1

!Set the number of children cells to be checked
integer,parameter :: Norder_cap = 1

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
   logical :: testtrajDetailedRMSD_flag = .false.

   !This takes much more time but you can force the checkState subroutine
   !to also check the rmsd of frames in adjacent cells
   logical :: force_Neighbors = .true.

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

      !Set .true. if the real force calculations we do should be
      !added to the grid
      logical :: grid_addition = .true.

      !$OMP THREADPRIVATE(reject_flag,accept_first,accept_worst)

         !Set .true. if interpolation should be used; that is to say
         !a weighted combination of acceptable frames are used to
         !calculate an approximate gradient
         logical :: interpolation_flag = .false.

         !Interpolation requires a scaling parameter for the weights
         !This is a positive, nonzero real number
         real(dp) :: interpolation_alpha1 = 2.0d0

         !For persistent data collection and analysis, a new file
         !is dedicated to interpolation results
         !Whether we record or not to this file is governed by the
         !gather_interpolation_flag
         logical :: gather_interpolation_flag = .false.
         character(17),parameter :: interpolationfile = "interpolation.dat"

            !Interpolation data can be checked at any point
            !whenever the interpolation counter reaches the
            !interpolation check, the data is checked
            !If the visual flag is true, then a visual is also
            !made when the data is checked
            integer :: interpolation_check = 1000
            integer :: interpolation_counter
            logical :: interpolation_check_visual = .true.

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

!Some variables used to calculate the KRP
!Useful to be global so that the convergence can
!communicate with the comparison
real :: lambda_penalty
real :: minsd_penalty

integer,parameter :: comparison_number = 1
character(11),parameter :: allprefixes = "placeholder"
integer,dimension(comparison_number),parameter :: alllengths = (/ 69 /)



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

write(Nthreshold_text,FMT=FMT6_pos_real0) threshold_rmsd
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



