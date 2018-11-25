module PARAMETERS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      PATHS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Path to the trajectories
!character(28),parameter :: path1 = "/home/ruisun/proj/ruisun/B0/"
 character(25),parameter :: path_to_existing_frames = "/home/kazuumi/Desktop/B0/"
!Path to the f90s and bash scripts
!character(44),parameter :: path2 = "/home/kazuumi/lus/B0/branch1/GradientGrider/"
!character(37),parameter :: path2 = "/home/kazuumi/Desktop/GradientGrider/"
 character(48),parameter :: path_to_source = "/home/kazuumi/Desktop/GradientGrider/playground/"
!character(42),parameter :: path_to_gnuplot= "/lus/scratch/usr/ruisun/gnuplot-5.2.4/bin/"
 character(0),parameter :: path_to_gnuplot= ""



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      FILENAMES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!File that will keep the trajectory folder names
character(16),parameter :: trajectories = "trajectories.txt"
!File that writes the progress of the program
character(12),parameter :: progressfile = "progress.txt"
!File that writes the progress of a trajectory
character(14),parameter :: trajectoryfile = "trajectory.xyz"
!File that writes the rmsd retrieved from checkState every frame
character(14),parameter :: checkstatefile = "checkstate.dat"
!File that writes the progress of multiple trajectories
character(16),parameter :: trajectoriesfile = "trajectories.dat"
!File that has the parameters (self)
character(14),parameter :: parametersfile = "PARAMETERS.f90"
!File for any gnuplot scripting
character(11),parameter :: gnuplotfile = "gnuplotfile"
!File for storing trajectory data from across many grids
character(23),parameter :: cumulativefile = "cumulative_trajectories"
!File for storing the first and last frames of multiple trajectories
character(13),parameter :: timeslicefile = "timeslice.dat"
!File for storing information on the grid over multiple trajectories
character(15),parameter :: informaticsfile = "informatics.dat"
!File for storing scattering angle data over multiple trajectories
character(15),parameter :: SAfile = "SA_trajectories.dat"
!File for storing initial bonding data over multiple trajectories
character(22),parameter :: initialfile = "initbonds_trajectories.dat"
!File for storing TRV energy changes
character(16),parameter :: TRVfile = "TRV_trajectories.dat"
!File for both Eenrgy Decomposition and Scattering Angle over multiple trajectories
character(18),parameter :: SATRVfile = "SATRV_trajectories.dat"
character(24),parameter :: binnedSATRVfile = "binnedSATRV_trajectories.dat"

!Files where, when done, we save the counters with the numbers of frames
!inside each subcell, as well as other information
character(12),parameter :: counter0file = "counter0.txt"
character(12),parameter :: counter1file = "counter1.txt"
character(12),parameter :: counter2file = "counter2.txt"
character(12),parameter :: counter3file = "counter3.txt"
!File to write to for system calls
character(8),parameter :: temporaryfile1 = "tmp1.txt"
character(8),parameter :: temporaryfile2 = "tmp2.txt"
character(8),parameter :: temporaryfile3 = "tmp3.txt"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     FILE CHANNELS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Channel for the progressfile
integer,parameter :: progresschannel = 70
integer,parameter :: trajectorieschannel = 71
integer,parameter :: frameschannel = 72
integer,parameter :: filechannel1 = 73
integer,parameter :: filechannel2 = 74
integer,parameter :: filechannel3 = 75
integer,parameter :: filechannel4 = 76
integer,parameter :: gnuplotchannel = 77


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      FORMATTING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

logical,parameter :: unreadable_flag = .true.

character(14),parameter :: FMT1 = "(2(1x,F14.10))"          !For two variables
character(22) :: FMT2         !For six atoms, coords and gradient
character(22) :: FMT3         !For six atoms, coords or gradient
character(4),parameter :: FMT4 = "(I4)"                     !For subcell names
character(4),parameter :: FMT5 = "(I9)"                     !For number of substates
character(7),parameter :: FMT6 = "(F12.8)"                  !For RMSD
character(19),parameter :: FMT7 = "(30x,18(1x,F14.10))"     !For skipping the variables

character(4),parameter :: FMT5_variable = "(I5)"
character(4),parameter :: FMT8_counter = "(I8)"
character(6),parameter :: FMT6_pos_real0 = "(F6.5)"
character(6),parameter :: FMT6_pos_real1 = "(F6.4)"
character(6),parameter :: FMT6_neg_real1 = "(F6.3)"
character(6),parameter :: FMT6_pos_int = "(I0.6)"
character(6) :: Nbond_text
character(17) :: FMTinitial
character(6) :: Natom_text
character(19) :: FMTtimeslice
character(27),parameter :: FMTinformatics = "(2(F12.7),I6,2(I5),I8,F8.4)"
character(15),parameter :: FMTsa = "((F6.4),(F8.4))"
character(10),parameter :: FMTtrv = "(3(F11.6))"
character(32),parameter :: FMTdata = "((F6.4),(F8.4),2(F11.6),(F13.9))"
character(33),parameter :: FMTnow = "('Time: ',I2.2,':',I2.2,':',I2.2)"

integer,parameter :: trajectory_text_length = 5
integer,parameter :: scaling1_text_length = 3
integer,parameter :: scaling2_text_length = 3
integer,parameter :: overcrowd0_text_length = 5
integer,parameter :: Ngrid_text_length = 3
integer,parameter :: gridpath_length = 68
character(gridpath_length),parameter :: gridpath0 = ""
integer, parameter :: SAfiles_text_length = gridpath_length +&
                                                 4 + 12 + len(SAfile) + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!The number of atoms in the system
integer,parameter :: Natoms = 3
integer,parameter :: Ncoords = Natoms*3

!The number of variables in use
integer,parameter :: Nvar = 2
integer,parameter :: Nvar_next = Nvar + 1
integer,parameter :: Ncoordsvals = Nvar+Ncoords
integer,parameter :: Ncoordsvals_next = Ncoordsvals + 1
integer,parameter :: Ngradientcoordsvals = Nvar+Ncoords*2







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     GRID PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Norder_max controls how many children generation we will make
integer,parameter :: Norder_max = 2
logical,parameter :: order0_flag = (1 > Norder_max)
logical,parameter :: order1_flag = (2 > Norder_max)
logical,parameter :: order2_flag = (3 > Norder_max)

!The spacing is the spacing of the parent-level grid
!Because the variables may be unbound, we define the
!parent-level grid in terms of gridline spacing
real,parameter :: spacing1 = 0.01        !Angstroms
real,parameter :: spacing2 = 0.01        !Angstroms
real,parameter :: spacing3 = 0.1         !Not in use

!There are some outliers; making a maximum throws these away
!real,parameter :: max_var2 = 10.0	!Angstroms
real,parameter :: max_var1 = 10.0       !Angstroms
real,parameter :: max_var2 = 12.0	!Angstroms (cos(angle) + 1 \in [0,2])
!real,parameter :: max_var2 = 2.0	!Unitless (cos(angle) + 1 \in [0,2])
integer,parameter :: bounds1 = ceiling(max_var1/spacing1)+400 !if var > max_var every 500 steps
integer,parameter :: bounds2 = ceiling(max_var2/spacing2)+400 !if var > max_var every 500 steps
!integer,parameter :: bounds2 = ceiling(max_var2/spacing2)+4	!For the cos(theta) definition
!integer,parameter :: bounds2 = ceiling(max_var1/spacing1)+400 !This 4 is a cushioning because we only check
!Consequently, we know the maximum number of cells in the grid
integer, parameter :: counter0_max = bounds1*bounds2

!The threshold of "overcrowded" for a cell of order N
integer,parameter :: overcrowd0 = 50
integer,parameter :: overcrowd1 = 10000
integer,parameter :: overcrowd2 = 50
integer,parameter :: overcrowd3 = 50

!The scaling is the amount that is resolved for an overcrowded cell
! (10,10,10) = x1000 magnification
! For now, only two variables are used
integer,parameter :: scaling1_0 = 4
integer,parameter :: scaling2_0 = 4
integer,parameter :: resolution_0 = scaling1_0*scaling2_0
integer,dimension(1+Nvar),parameter :: SP0=(/scaling1_0,scaling2_0,resolution_0/)

integer,parameter :: scaling1_1 = 10
integer,parameter :: scaling2_1 = 10
integer,parameter :: resolution_1 = scaling1_1*scaling2_1
integer,dimension(1+Nvar),parameter :: SP1=(/scaling1_1,scaling2_1,resolution_1/)

integer,parameter :: scaling1_2 = 10
integer,parameter :: scaling2_2 = 10
integer,parameter :: resolution_2 = scaling1_2*scaling2_2
integer,dimension(1+Nvar),parameter :: SP2=(/scaling1_2,scaling2_2,resolution_2/)

!The formatting of the subcell of a particular order
character(6), parameter :: FMTorder0 = "(F9.2)"
character(6), parameter :: FMTorder1 = "(F9.4)"
character(6), parameter :: FMTorder2 = "(F9.5)"
character(6), parameter :: FMTorder3 = "(F9.6)"
integer,parameter :: FMTlength = 9


!Some useful constants to have calculated beforehand
!The variable multiplier is the length of the subcell of a particular order
real, parameter :: multiplier1_0 = spacing1
real, parameter :: multiplier2_0 = spacing2
real, parameter :: multiplier1_1 = spacing1/scaling1_0
real, parameter :: multiplier2_1 = spacing2/scaling2_0
real, parameter :: multiplier1_2 = spacing1/(scaling1_0*scaling1_1)
real, parameter :: multiplier2_2 = spacing2/(scaling2_0*scaling2_1)
real, parameter :: multiplier1_3 = spacing1/(scaling1_0*scaling1_1*scaling1_2)
real, parameter :: multiplier2_3 = spacing2/(scaling2_0*scaling2_1*scaling2_2)

real, parameter :: divisor1_0  = 1.0 / multiplier1_0
real, parameter :: divisor2_0  = 1.0 / multiplier2_0
real, parameter :: divisor1_1  = 1.0 / multiplier1_1
real, parameter :: divisor2_1  = 1.0 / multiplier2_1
real, parameter :: divisor1_2  = 1.0 / multiplier1_2
real, parameter :: divisor2_2  = 1.0 / multiplier2_2
real, parameter :: divisor1_3  = 1.0 / multiplier1_3
real, parameter :: divisor2_3  = 1.0 / multiplier2_3

real,dimension(2*Nvar),parameter :: MP0=(/multiplier1_0,multiplier2_0,&
					  divisor1_0,divisor2_0/)
real,dimension(2*Nvar),parameter :: MP1=(/multiplier1_1,multiplier2_1,&
					  divisor1_1,divisor2_1/)
real,dimension(2*Nvar),parameter :: MP2=(/multiplier1_2,multiplier2_2,&
					  divisor1_2,divisor2_2/)
real,dimension(2*Nvar),parameter :: MP3=(/multiplier1_3,multiplier2_3,&
					  divisor1_3,divisor2_3/)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    MEMORY OVERHEAD COST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!We need to estimate how many cells will be overcrowded in counter0
!Must not be greater than 99999*resolution
integer,parameter :: header1_max = 90000
integer,parameter :: counter1_max = header1_max*resolution_0
!We need to estimate how many subcells will be overcrowded in counter1
!Must not be greater than 99999*resolution
integer,parameter :: header2_max = 10
integer,parameter :: counter2_max = header2_max*resolution_1
!We need to estimate how many subcells will be overcrowded in counter2
!Must not be greater than 99999*resolution
integer,parameter :: header3_max = 10
integer,parameter :: counter3_max = header3_max*resolution_2

!The counter arrays hold information on the population of the subcell
!And the index to get to deeper subcells in the following counter
!If the value is XXXXYYYY, then XXXX is the key to the next
!counter and YYYY is the population of the current counter
integer,parameter :: population_max = 999
integer,parameter :: key_start = population_max + 1

!This takes much more time but you can force the checkState subroutine
!to also check the rmsd of frames in adjacent cells
!Not heavily tested in the newer updates so is deprecated (for now)
logical,parameter :: force_Neighbors = .false.

!If we want to add duplicate copies of a state (but permuted labels)
!Then set this to .true.
!Otherwise, all frames are relabelled to a specific format
!Note: if true, this will fill up the grid much faster so a larger
!value of overcrowd is recommended
logical,parameter :: force_Duplicates = .false.

!If force_Duplicate is set to false there is an additonal feature here
!for the frame to be added to have no relabelling
!Note: this makes it so that two indistinguishable frames
!      possibly do not have the same RMSD or even cell index
logical,parameter :: force_NoLabels = .false.



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   MULTI-GRID PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!The number of trajectories we add to a grid
integer,parameter :: Ntraj_max = 1000
!The maximum amount of time (seconds) we are willing to wait for a single trajectory to finish
real,parameter :: trajectory_CPU_time_max = 60.0
!The number of grids we will make
integer :: Ngrid_max = 1
!The number of trajectories to make before checking the grid-making progress
integer,parameter :: Ngrid_check_min = 1
integer,parameter :: Ngrid_check = max(Ntraj_max/10,Ngrid_check_min)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      COUNTERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer :: header1, header2, header3
integer, dimension(counter0_max) :: counter0
integer, dimension(counter1_max) :: counter1
integer, dimension(counter2_max) :: counter2
integer, dimension(counter3_max) :: counter3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  GRID FORMATTING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer :: Ngrid
character(gridpath_length+Ngrid_text_length+1) :: gridpath1
character(gridpath_length+Ngrid_text_length+1+5) :: gridpath2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    PLOT FORMATTING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

character(6) :: checkstateTrajectory
character(6) :: angle1descriptor,angle2descriptor,bond1descriptor,scatteringdescriptor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    TRAJECTORY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer :: Ntraj,Ntraj_allowed,Nfile,steps,heatmap_steps,Norder1
logical :: header_max_flag

!Set .true. to generate top-level heat map GIFS over time for each grid
logical,parameter :: heatmap_evolution_flag = .false.
integer,parameter :: heatmap_evolution_steps = 100





















end module PARAMETERS



