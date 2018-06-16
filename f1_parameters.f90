module f1_parameters
implicit none


!PATHS
!Path to the trajectories
!character(28),parameter :: path1 = "/home/ruisun/proj/ruisun/B0/"
character(25),parameter :: path1 = "/home/kazuumi/Desktop/B0/"
!Path to the f90s and bash scripts
!character(44),parameter :: path2 = "/home/kazuumi/lus/B0/branch1/GradientGrider/"
character(37),parameter :: path2 = "/home/kazuumi/Desktop/GradientGrider/"
!Path to the grid
!character(39),parameter :: path3 = "/home/kazuumi/lus/B0/branch1/grid_test/"
character(47),parameter :: path3 = "/home/kazuumi/Desktop/GradientGrider/grid_test/"
!Path for the following files (see below)
!character(29),parameter :: path4 = "/home/kazuumi/lus/B0/branch1/"
character(45),parameter :: path4 = "/home/kazuumi/Desktop/GradientGrider/f1_dump/"



!FILES
!File that will keep the trajectory folder names
character(19),parameter :: trajectories = "f1_trajectories.txt"
!File that writes the progress of the program
character(15),parameter :: progressfile = "f1_progress.txt"
!File that writes the progress of a trajectory
character(17),parameter :: trajectoryfile = "f1_trajectory.xyz"
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

!FILE CHANNELS
!Channel for the progressfile
integer,parameter :: progresschannel = 70
integer,parameter :: trajectorieschannel = 71
integer,parameter :: frameschannel = 72
integer,parameter :: filechannel1 = 73
integer,parameter :: filechannel2 = 74




!If we need to create a grid from the trajectories in path1, then set this true
logical :: start_from_scratch = .true.



!FORMATS
character(13),parameter :: FMT1 = "(3(1x,F11.6))"           !For three variables
character(14),parameter :: FMT2 = "(36(1x,F11.6))"          !For six atoms, full state
character(14),parameter :: FMT3 = "(18(1x,F11.6))"          !For six atoms, coords
character(4),parameter :: FMT4 = "(I4)"                     !For subcell names
character(4),parameter :: FMT5 = "(I9)"                     !For number of substates



!VARIABLES

!The number of atoms in the system
integer,parameter :: Natoms = 6

!The number of variables in use
integer,parameter :: Nvar = 3

!The spacing is the spacing of the parent-level grid
!Because the variables may be unbound, we define the
!parent-level grid in terms of gridline spacing
real,parameter :: spacing1 = 3.6        !Angstroms
real,parameter :: spacing2 = 3.6        !Angstroms
real,parameter :: spacing3 = 0.1        !Angstroms

!There are some outliers; making a maximum throws these away
real,parameter :: max_var1 = 40.0       !Angstroms
real,parameter :: max_var2 = 40.0       !Angstroms
!Consequently, we know the maximum number of cells in the grid
integer, parameter :: counter0_max = anint(max_var1/spacing1)&
                                     *anint(max_var2/spacing2)

!The threshold of "overcrowded" for a cell of order N
integer,parameter :: overcrowd0 = 25
integer,parameter :: overcrowd1 = 25
integer,parameter :: overcrowd2 = 25
integer,parameter :: overcrowd3 = 25

!The scaling is the amount that is resolved for an overcrowded cell
! (10,10,10) = x1000 magnification
! For now, only two variables are used
integer,parameter :: scaling1_0 = 3
integer,parameter :: scaling2_0 = 3
integer,parameter :: resolution_0 = scaling1_0*scaling2_0
integer,dimension(3),parameter :: SP0=(/scaling1_0,scaling2_0,resolution_0/)

integer,parameter :: scaling1_1 = 10
integer,parameter :: scaling2_1 = 10
integer,parameter :: resolution_1 = scaling1_1*scaling2_1
integer,dimension(3),parameter :: SP1=(/scaling1_1,scaling2_1,resolution_1/)

integer,parameter :: scaling1_2 = 10
integer,parameter :: scaling2_2 = 10
integer,parameter :: resolution_2 = scaling1_2*scaling2_2
integer,dimension(3),parameter :: SP2=(/scaling1_2,scaling2_2,resolution_2/)

!The formatting of the subcell of a particular order
character(6), parameter :: FMTorder0 = "(F9.2)"
character(6), parameter :: FMTorder1 = "(F9.3)"
character(6), parameter :: FMTorder2 = "(F9.4)"
character(6), parameter :: FMTorder3 = "(F9.5)"


!Some useful constants to have calculated beforehand
!The variable multiplier is the length of the subcell of a particular order
real, parameter :: multiplier1_0 = spacing1
real, parameter :: multiplier2_0 = spacing2
real, parameter :: multiplier1_1 = spacing1/scaling1_0
real, parameter :: multiplier2_1 = spacing1/scaling2_0
real, parameter :: multiplier1_2 = spacing1/(scaling1_0*scaling1_1)
real, parameter :: multiplier2_2 = spacing1/(scaling2_0*scaling2_1)
real, parameter :: multiplier1_3 = spacing1/(scaling1_0*scaling1_1*scaling1_2)
real, parameter :: multiplier2_3 = spacing1/(scaling2_0*scaling2_1*scaling2_2)

real, parameter :: divisor1_0  = 1.0 / multiplier1_0
real, parameter :: divisor2_0  = 1.0 / multiplier2_0
real, parameter :: divisor1_1  = 1.0 / multiplier1_1
real, parameter :: divisor2_1  = 1.0 / multiplier2_1
real, parameter :: divisor1_2  = 1.0 / multiplier1_2
real, parameter :: divisor2_2  = 1.0 / multiplier2_2
real, parameter :: divisor1_3  = 1.0 / multiplier1_3
real, parameter :: divisor2_3  = 1.0 / multiplier2_3





!We need to estimate how many cells will be overcrowded in counter0
!Must not be greater than 99999*resolution
integer,parameter :: counter1_max = 1000*resolution_0
!We need to estimate how many subcells will be overcrowded in counter1
!Must not be greater than 99999*resolution
integer,parameter :: counter2_max = 5000*resolution_1
!We need to estimate how many subcells will be overcrowded in counter2
!Must not be greater than 99999*resolution
integer,parameter :: counter3_max = 3000*resolution_2

!The counter arrays hold information on the population of the subcell
!And the index to get to deeper subcells in the following counter
!If the value is XXXXYYYY, then XXXX is the key to the next
!counter and YYYY is the population of the current counter
integer,parameter :: population_max = 999
integer,parameter :: key_start = population_max + 1

end module f1_parameters



