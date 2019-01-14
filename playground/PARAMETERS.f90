module PARAMETERS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!In a parallel system, most of these variables are shared among threads
!Except for those declared wtih !$OMP THREADPRIVATE

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

!If true, the grid has unformatted files; thus, keep this constant
!The processed files and check trajectory files are still formatted though
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

!These formats depend on the number of atoms and bonds in the system
!So they cannot be parameters and must be changed later
!These are always changed in the MAIN program
character(6) :: Nbond_text
character(19) :: FMTinitial
character(6) :: Natom_text
character(19) :: FMTtimeslice
!$OMP THREADPRIVATE(FMT2,FMT3,Nbond_text,FMTinitial,Natom_text,FMTtimeslice)

!Formatting for files that are processed AFTER grid making and checking
character(27),parameter :: FMTinformatics = "(2(F12.7),I6,2(I5),I8,F8.4)"
character(15),parameter :: FMTsa = "((F6.4),(F8.4))"
character(10),parameter :: FMTtrv = "(3(F11.6))"
character(32),parameter :: FMTdata = "((F6.4),(F8.4),2(F11.6),(F13.9))"
character(33),parameter :: FMTnow = "('Time: ',I2.2,':',I2.2,':',I2.2)"

!Because we want to not use trim and adjustl all the time, and also
!because we want nice-looking graphs, we want to store various numbers
!in strings (with formattings as specified above)
!But to adjust to different sizes of strings, we need to establish
!the length of certain strings
integer,parameter :: trajectory_text_length = 5
integer,parameter :: Ngrid_text_length = 3

!Gridpath is one variable that must be supplied by the user and is changed
!on the LOCAL version of a library; this is changed AUTOMATICALLY
!But within the library itself it should stay constant
integer,parameter :: gridpath_length = 68
character(gridpath_length),parameter :: gridpath0 = ""
integer, parameter :: SAfiles_text_length = gridpath_length +&
                                                 4 + 12 + len(SAfile) + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!The number of atoms in the system
!Natoms is another variable that must be supplied by the user and is changed
!on the LOCAL version of a library; this must be changed MANUALLY
integer,parameter :: Natoms = 3
integer,parameter :: Ncoords = Natoms*3

!The number of variables in use
!Right now, we only support a two-variable grid (easier to script and visualize)
integer,parameter :: Nvar = 2
integer,parameter :: Nvar_next = Nvar + 1
integer,parameter :: Ncoordsvals = Nvar+Ncoords
integer,parameter :: Ncoordsvals_next = Ncoordsvals + 1
integer,parameter :: Ngradientcoordsvals = Nvar+Ncoords*2







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     GRID PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Many of the variables below are supplied by the user and are changed
!on the LOCAL version of alibrary; this is done AUTOMATICALLY

!Norder_max controls how many children generation we will make
integer,parameter :: Norder_max = 2


!Some default values are good if we don't know at execution
!how many variables we plan to use

real,parameter,dimension(4) :: default_spacing = &
        (/ 0.01, 0.01, 0.01, 0.01 /)
real,parameter,dimension(4) :: default_maxvar = &
        (/ 10.0, 12.0, 10.0, 10.0 /)
integer,parameter,dimension(4,4) :: default_scaling = &
        reshape( &
                (/  4,  4,  4,  4, &
                    10, 10, 10, 10, &
                    10, 10, 10, 10, &
                    10, 10, 10, 10    /), &
                (/ 4, 4 /))
integer,parameter,dimension(4) :: default_overcrowd = &
        (/     50,     10000,      50,       50 /)
integer,dimension(4),parameter :: default_Norder_order = &
        (/ 1, 0, 2, 3 /)

!After making these default values, we decide how many we
!can use, then simply truncate the extra variables off

real,parameter,dimension(Nvar) :: var_spacing = &
        default_spacing(1:Nvar)
real,parameter,dimension(Nvar) :: var_maxvar = &
        default_maxvar(1:Nvar)
integer,parameter,dimension(Nvar,Norder_max+1) :: var_scaling = &
        default_scaling(1:Nvar,1:Norder_max+1)
integer,parameter,dimension(Norder_max+1) :: var_overcrowd = &
        default_overcrowd(1:Norder_max+1)
integer,dimension(Norder_max+1) :: Norder_order = &
        default_Norder_order(1:Norder_max+1)

!We also need to define some default formats for cells
!of different orders;
!The singleFMT strings are formats for each variable, which
!will be combined later to formats for each cell
!
!One crucial thing here is that each string is of the length
!specified in singleFMT_length

integer,parameter :: singleFMT_length = 6
character(singleFMT_length*4),parameter :: default_singleFMT = &
        "(F0.2)"// &
        "(F0.4)"// &
        "(F0.5)"// &
        "(F0.6)"

character(singleFMT_length*(Norder_max+1)),parameter :: var_singleFMT = &
        default_singleFMT(1:singleFMT_length*(Norder_max+1))

integer,parameter :: multipleFMT_length = (Nvar) * singleFMT_length + &
                                          5 * (Nvar - 1) + 7 + 2
character(multipleFMT_length*(Norder_max+1)) :: var_multipleFMT


!The "bounds", or the upper limit to the indexing of
!the variables, is given a little padding (400 in this case)

integer,parameter,dimension(Nvar) :: var_bounds = &
        (/ (ceiling(var_maxvar(increment) / var_spacing(increment)) + 400, increment=1,Nvar) /)

!The "resolution" is how many new cells are created when a cell
!of a certain order has divyUp called on it

integer,parameter,dimension(Norder_max+1) :: var_resolution = &
        (/ (product(var_scaling(:,increment),DIM=1), increment=1,Norder_max+1) /)

!And then some scaling factors useful to save calculation time later

real,dimension(Nvar,Norder_max+1) :: multiplier
real,dimension(Nvar,Norder_max+1) :: divisor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    MEMORY OVERHEAD COST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!This takes much more time but you can force the checkState subroutine
!to also check the rmsd of frames in adjacent cells
!Not heavily tested in the newer updates so is deprecated (for now)
logical,parameter :: force_Neighbors = .true.

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

!Ntraj_max, Ngrid_max, and Ngrid_check_min are three variables supplied by
!the user to change the grid-making procedure; this is updated AUTOMATICALLY

!The number of trajectories we add to a grid
integer,parameter :: Ntraj_max = 1000
!The maximum amount of time (seconds) we are willing to wait for a single trajectory to finish
real,parameter :: trajectory_CPU_time_max = 60.0
!The number of grids we will make
integer :: Ngrid_max = 1
!$OMP THREADPRIVATE(Ngrid_max)

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

integer,dimension(Norder_max+1) :: headers
integer,dimension(Norder_max+1) :: headers_old
integer,dimension(Norder_max+1) :: Norder_total

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      TRAVERSAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Traversal variables deal with counting how many cells are visited by a trajectory;
!The first column is counter0/1_index, the second column is thread (for parallel programming)

!This is not implemented for grid-making but is for grid-checking
!Right now, we only have output for traversal1 (number of subcells traversed)
!but we can easily implement traversal0 (number of parent cells traversed) output

!Traversal is slightly costly so if you're not interested in traversal,
!please flag this false MANUALLY

logical, parameter :: traversal_flag = .false.
integer, allocatable :: traversal0(:,:)
integer, allocatable :: traversal1(:,:)
!$OMP THREADPRIVATE(traversal0,traversal1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  GRID FORMATTING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Because Ngrid is a shared variable, but it is often manipulated within runTrajectory,
!we need each thread to have its own private copy of Ngrid
!(and consequently, gridpath1 and gridpath2)

integer :: Ngrid
character(gridpath_length+Ngrid_text_length+1) :: gridpath1
character(gridpath_length+Ngrid_text_length+1+5) :: gridpath2, gridpath3
!$OMP THREADPRIVATE(Ngrid,gridpath1,gridpath2,gridpath3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    PLOT FORMATTING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!These are some strings that are used often in analyzeScatteringAngle;
!because they are used often, I just define them here instead of redefining
!them each time inside of each subroutine

character(6) :: checkstateTrajectory
character(6) :: angle1descriptor,angle2descriptor,bond1descriptor,scatteringdescriptor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    TRAJECTORY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Every trajectory, we keep track of various metadata about the trajectory and file system in general like:
!          Ntraj: the "ID" of the trajectory
!  Ntraj_allowed: I've forgotten what this is
!          Nfile: the number of files in the filesystem
!          steps: the number of steps that have passed for a trajectory
!  heatmap_steps: the number of steps that have passed in the trajectory (for heatmap_evolution)
!        Norder1: the number of subcells (cells of order = 1) that were checked
!header_max_flag: a flag that signals that counter1 is full (need to check this later)

integer :: Ntraj,Ntraj_allowed,Nfile,steps,heatmap_steps,Norder1
logical :: header_max_flag
!$OMP THREADPRIVATE(Ntraj,Ntraj_allowed,Nfile,steps,Norder1)

!Set .true. to generate top-level heat map GIFS over time for each grid
logical,parameter :: heatmap_evolution_flag = .false.
integer,parameter :: heatmap_evolution_steps = 100





















end module PARAMETERS



