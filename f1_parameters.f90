module f1_parameters
implicit none


!PATHS
!Path to the trajectories
character(50),parameter :: path1 = "/home/ruisun/proj/ruisun/B0/"
!Path to the f90s and bash scripts
character(50),parameter :: path2 = "/home/kazuumi/lus/B0/branch1/GradientGrider/"
!Path to the grid
character(50),parameter :: path3 = "/home/kazuumi/lus/B0/branch1/grid2/"
!Path for the following files (see below)
character(50),parameter :: path4 = "/home/kazuumi/lus/B0/branch1/"



!FILES
!File that will keep the trajectory folder names
character(20),parameter :: trajectories = "f1_trajectories.txt" !Used to be file1
!File that writes the progress of the program
character(20),parameter :: progressfile = "f1_progress.txt" !Used to not exist
!Files where, when done, we save the counters with the numbers of frames
!inside each subcell, as well as other information
character(20),parameter :: counter0file = "counter0.txt"
character(20),parameter :: counter1file = "counter1.txt"
character(20),parameter :: counter2file = "counter2.txt"
character(20),parameter :: counter3file = "counter3.txt"
!File to write to for system calls
character(20),parameter :: temporaryfile1 = "tmp1.txt"
character(20),parameter :: temporaryfile2 = "tmp2.txt"
character(20),parameter :: temporaryfile3 = "tmp3.txt"




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

!The spacing is the dimension of the gridlines; currently var3 is not grided
real,parameter :: spacing1 = 1.0        !Angstroms
real,parameter :: spacing2 = 1.0        !Angstroms
real,parameter :: spacing3 = 0.1        !Angstroms

!The threshold of "overcrowded" for a cell of order N
integer,parameter :: overcrowd0 = 2000
integer,parameter :: overcrowd1 = 500
integer,parameter :: overcrowd2 = 100
integer,parameter :: overcrowd3 = 100

!The scaling is the amount that is resolved for an overcrowded grid
! (10,10,10) = x1000 magnification
! This might be a little too many
integer,parameter :: scaling1 = 10
integer,parameter :: scaling2 = 10
integer,parameter :: scaling3 = 10
integer,parameter :: resolution = scaling1*scaling2

!There are some outliers; making a maximum throws these away
real,parameter :: max_var1 = 80.0       !Angstroms
real,parameter :: max_var2 = 80.0       !Angstroms

!We need to estimate how many cells will be overcrowded in counter0
!Must not be greater than 999*resolution
integer,parameter :: counter1_max = 250*resolution
!We need to estimate how many subcells will be overcrowded in counter1
!Must not be greater than 999*resolution
integer,parameter :: counter2_max = 800*resolution
!We need to estimate how many subcells will be overcrowded in counter2
!Must not be greater than 999*resolution
integer,parameter :: counter3_max = 250*resolution

!When we get to this 'order' or deepness in the subcells, we have a good enough
!approximation; otherwise, we must calculate the gradient by hand.
!As example, a subcell r1_r2.dat inside no subfolder has order=0
integer,parameter :: acceptable_depth = 3

!The threshold for two numbers being equal
!Currently not in use
real,parameter :: dx = .000006

end module f1_parameters



