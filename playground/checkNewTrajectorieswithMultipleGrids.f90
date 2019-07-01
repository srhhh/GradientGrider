!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PROGRAM
!               checkNewTrajectorieswithMultipleGrids
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!		This program performs an analysis of a library of grids
!		specified in PARAMETERS, with the analysis guided by
!		logical flags in ANALYSIS
!
!		Often, new trajectories will have to be made with certain
!		rejection method and RMSD threshold; these sample trajectories
!		have each frame checked with the grids, and an approximate
!		gradient used if within threshold; these frames are not added
!		to the grids.
!
!		The graphs that output include: heat maps, scattering angle
!		plots, and percent rmsd threshold plots
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILECHANNELS                    ACTION
!
!               GNUPLOTCHANNEL                  OPEN, WRITE, CLOSE
!               FILECHANNEL1                    OPEN, WRITE, CLOSE
!		TRAJECTORIESCHANNEL		OPEN, WRITE, CLOSE
!		FILECHANNELS(1)          	none
!		FILECHANNELS(2:1+Ngrid_total)	OPEN, RUNTRAJECTORY, CLOSE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       CALLS                           MODULE
!
!               SYSTEM                          INTRINSIC
!               CPU_TIME                        INTRINSIC
!               SYSTEMCLOCK                     INTRINSIC
!
!               checkMultipleTrajectories       runTrajectory
!
!               getScatteringAngles1            analyzeScatteringAngleswithMultipleGrids
!               getScatteringAngles2            analyzeScatteringAngleswithMultipleGrids
!		analyzeHeatMaps1		analyzeHeatMapswithMultipleGrids
!		getRMSDThresholds1		analyzeRMSDThresholdwithMultipleGrids
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!		gridpath0//trajectories		DAT			A list of grids (001/, 002/, ...)
!		gridpath1//reject//		DAT			The min_rmsd retrieved for each step
!			threshold//TRAJ#				of a trajectory
!		gridpath0//Ngrid///reject//	DAT			A list of data from all trajectories with the
!		     threshold//trajectories				same suffix (Ngrid,reject,threshold)
!		gridpath1//HeatMap		JPG			The heat map of a grid
!		gridpath0//TrueSADist//		JPG			The scattering angle distribution of
!			TOTALTRAJ#					TOTALTRAJ# trajectories (multiple grids)
!		gridpath0//TestSADist//		JPG			The scattering angle distribution of all
!		     Ngrid//reject//threshold				trajectoreis with the same suffix
!		gridpath1//PercentRMSDDist//	JPG			The percent RMSD threshold distribution of
!		     Ngrid//reject//threshold				all trajectories with the same suffix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program checkNewTrajectorieswithMultipleGrids
use ANALYSIS
use PARAMETERS
use interactMultipleGrids
use runTrajectory
use analyzeHeatMapswithMultipleGrids
use analyzeRMSDThresholdwithMultipleGrids
use analyzeScatteringAngleswithMultipleGrids
implicit none

!Variables used to name the new directory and trajectories uniquely
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectory_text_length) :: Ntraj_text
character(6) :: Nthreshold_text
character(6) :: reject_text
character(12) :: short_prefix_text
character(12+Ngrid_text_length) :: prefix_text
character(150) :: old_filename, new_filename
character(1) :: answer

real(dp) :: scattering_angle
real(dp),dimension(3) :: TRVenergies1,TRVenergies2,dTRVenergies
real(dp),dimension(3,Natoms) :: coords_initial,velocities_initial
real(dp),dimension(3,Natoms) :: coords_final,velocities_final
integer :: seed,n,m,n_testtraj,initial_n_testtraj
real :: lowerlimit,upperlimit

!Variables
integer :: iostate
integer,allocatable :: filechannels(:)
integer :: OMP_GET_THREAD_NUM
logical :: return_flag = .false.
logical :: makeCheckTrajectory_flag = .false.
logical :: file_exists

integer :: traj_index
character(6) :: traj_index_text
integer :: inputchannel = 307

!Timing
real :: r1, r2, system_clock_rate
integer :: c1, c2, cr
integer,dimension(3) :: now
real :: grid_wall_time,checktrajectory_wall_time
character(10) :: grid_wall_time_text, checktrajectory_wall_time_text
integer :: trajectory_t0, trajectory_t1

!Incremental Integers
integer :: i, j, k


!Get the number of threads set up
call OMP_SET_NUM_THREADS(Nthreads)

!Get our timing system ready
call system_clock(c1,count_rate=cr)
system_clock_rate = 1.0/real(cr)

!$OMP PARALLEL DEFAULT(shared)&
!$OMP& PRIVATE(iostate)&
!$OMP& PRIVATE(variable_length_text,Ntraj_text,filechannels)&
!$OMP& FIRSTPRIVATE(Ngrid_text,reject_text)&
!$OMP& PRIVATE(c1,c2,cr,r1,r2,coords_initial,coords_final,velocities_initial,velocities_final)&
!$OMP& SHARED(return_flag,now,seed,system_clock_rate)&
!$OMP& SHARED(initial_n_testtraj,folder_text,prefix_text,short_prefix_text)&
!$OMP& SHARED(Nthreshold_text,old_filename,new_filename,answer)&
!$OMP& SHARED(var_multipleFMT,multiplier,divisor)&
!$OMP& PRIVATE(i,j,n)

!$OMP SINGLE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       RNG SETUP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!And we print the seed in case there's some bug that we need to reproduce
call system_clock(seed)
print *, ""
call itime(now)
write(6,FMT=FMTnow) now
print *, "   RNG seed: ", seed
print *, ""
seed = rand(seed)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       SCALING SETUP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!We also need to initialize a few scaling factors to
!speed up computation;
!These are multidimensional arrays that are difficult to
!create with a fortran array constructor
multiplier(:,1) = var_spacing(:)
do m = 1, Norder_max
do n = 1, Nvar
    multiplier(n,m+1) = var_spacing(n) / product(var_scaling(n,1:m),DIM=1)
end do
end do

do m = 1, Norder_max+1
do n = 1, Nvar
    divisor(n,m) = 1.0 / multiplier(n,m)
end do
end do

!Because the variables subcellsearchmax is in PARAMETERS,
!we must initialize it from whatever is given in ANALYSIS
call setSubcellSearchMax()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       FORMATTING SETUP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Cells of different orders (ex. parent vs child) have
!different formats which we obtain by combining together
!the formats of each variable according to which order it is
var_multipleFMT = ""
do m = 1, Norder_max+1
    var_multipleFMT = trim(adjustl(var_multipleFMT))//"("
    
    do n = 1, Nvar
        var_multipleFMT = trim(adjustl(var_multipleFMT))//&
        var_singleFMT(1+singleFMT_length*(m-1):&
                        singleFMT_length*m)
    
        if (n == Nvar) exit
    
        var_multipleFMT = trim(adjustl(var_multipleFMT))//&
                ',"_",'
    end do
    
    var_multipleFMT = trim(adjustl(var_multipleFMT))//',".dat")'
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       PRE-PROCESSING LIBRARY OVERVIEW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!All trajectory folders are formatted as I0.3 (3-digit integer)
!So search for these numbered folders and read them
call system("ls -p "//gridpath0//&
        " | grep '^[0123456789]*/' > "//&
        gridpath0//trajectories)

!Let's see how many grids are actually in the folder
Ngrid_max = 0
open(trajectorieschannel,file=gridpath0//&
        trajectories,action="read")
do
    read(trajectorieschannel,FMT="(A4)",&
            iostat=iostate) folder_text
    if (iostate /= 0) exit
    Ngrid_max = Ngrid_max + 1
end do
close(trajectorieschannel)

!Keep the cap to the number of grids in mind
Ngrid_total = min(Ngrid_cap, Ngrid_max)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       CHECKPOINT INDICATOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!If we want this program to be a "pick up where we
!left off last time" program we figure out how much
!progress has already been done for this experiment

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

inquire(file=gridpath5//initialfile,exist=file_exists)

!By default, we say we will remake all the trajectories
initial_n_testtraj = 1

!This is a "new experiment" if we are NOT continuing
!off from our last analysis OR if the experiment in
!question has NOT been started
if ((.not.(continue_analysis)).or.(.not.file_exists)) then

    !For a new experiment, remake all our data
    call system("rm -r "//gridpath5)
    call system("mkdir "//gridpath5)

!Otherwise, we need to figure out where to start
else
    print *, "     Deciding to continue previous analysis..."

    !First we check how many trajectory files we have
    !This system call puts them all onto one file
    !We need all grids 1 to Ngrid_total to have the
    !trajectories but we just check the last one
    !and we assume all previous grids have just as
    !many or more (usually a safe assumption)

    write(variable_length_text,FMT=FMT5_variable)&
            Ngrid_text_length
    write(Ngrid_text,FMT="(I0."//&
            trim(adjustl(variable_length_text))//")")&
            Ngrid_total
    write(variable_length_text,FMT=FMT5_variable)&
            trajectory_text_length
    call system("ls "//gridpath5//" | grep -E '^"//&
            Ngrid_text//"_[0123456789]{"//&
            trim(adjustl(variable_length_text))//&
            "}.dat' > "//gridpath5//trajectories)

    !Then we simply have to read the file to see
    !how many we have. This also assumes they were
    !numbered correctly (usually a safe assumption)
    open(trajectorieschannel,file=gridpath5//&
            trajectories,action="read")
    do
        read(trajectorieschannel,*,iostat=iostate)
        if (iostate /= 0) exit
        initial_n_testtraj = initial_n_testtraj + 1
    end do
    close(trajectorieschannel)

    !Second, we check how many lines are on the timeslice file
    Ntraj = 1
    open(trajectorieschannel,file=gridpath5//timeslicefile)
    do
        read(trajectorieschannel,*,iostat=iostate)
        if (iostate /= 0) exit
        Ntraj = Ntraj + 1
    end do
    close(trajectorieschannel)

    !If these two values are not congruent there has
    !been some sort of file corruption or deletion
    if (Ntraj /= initial_n_testtraj) then
        print *, ""
        print *, ""
        print *, ""
        print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print *, "          DATA INCONGRUENCE "
        print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        if (Ntraj > initial_n_testtraj) print *, "    trajectory files missing for "//&
                                                 expfolder//" ... .dat"
        if (Ntraj < initial_n_testtraj) print *, "    timeslice file corrupted for "//&
                                                 expfolder
        print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print *, ""
        
        if ((Ntraj < Ntesttraj).or.&
            (initial_n_testtraj < Ntesttraj)) then
        do
            print *, "    corrupted data interrupts trajectory indexing;"
            print *, "    overwrite corrupted data? otherwise the program exits    (y/n)"
            print *, ""
            read (*,*) answer
            if ((answer == "y").or.(answer == "n")) exit
        end do

        !The user may choose to quit
        if (answer == "n") then
            print *, ""
            call itime(now)
            write(6,FMT=FMTnow) now
            print *, "    Program exited disgracefully"
            print *, ""

            !Instead of immediately quitting, we set
            !a flag to signal a quit later
            return_flag = .true.

        !Or to just keep appending new
        !data onto the timeslice file
        else
            print *, "    starting off where the timeslice file left off!"
            print *, ""
            initial_n_testtraj = min(initial_n_testtraj,Ntraj)

        end if
        else
            print *, "    going along with analysis on uncorrupted data"
            print *, ""
        end if
        print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print *, ""
    end if

    !If all is well, then we might have
    !save some trajectories (and time!)
    if (.not.return_flag) then
    print *, "     Total Number of Trajectories Saved: ",&
            initial_n_testtraj - 1, " out of ", Ntesttraj
    print *, ""
    end if
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP END SINGLE

!Initialize some essential strings
!and format strings across all threads
!Unfortunately, I don't know a better
!way to do this
write(Nbond_text,FMT="(I0.6)") Nbonds
write(Natom_text,FMT="(I0.6)") Natoms
write(FMTinitial,FMT="(A19)")&
        "("//Nbond_text//"(6(F14.10)))"
write(FMTtimeslice,FMT="(A19)")&
        "("//Natom_text//"(12(F12.7)))"
write(FMT2,FMT="(A22)")&
        "("//Natom_text//"(6(1x,F14.10)))"
write(FMT3,FMT="(A22)")&
        "("//Natom_text//"(3(1x,F14.10)))"

!!$OMP BARRIER
!!$OMP SINGLE

print *, ""
call itime(now)
write(6,FMT=FMTnow) now
print *, "Analysis on directory ", gridpath0
print *, "Deciding on using ", Ngrid_total, " grids"
print *, "Deciding on using ", Nthreads, " threads"
print *, ""


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       PRE-CREATION LIBRARY ANALYSIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!The only pre-creation analysis
!we have now is for a heatmap

!This is for top-level heat map
!generation (for each of the grids)
if (heatmap_flag) then
    call itime(now)
    write(6,FMT=FMTnow) now
    print *, "   Making plot: ", "TopLevel_HeatMap"
    print *, ""

    call analyzeTopLevelHeatMaps("TopLevel_HeatMap")
end if

!!$OMP END SINGLE
!!$OMP BARRIER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       TRAJECTORY CREATION FLAG START
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!This is for checking trajectories against the grid
!Currently, this is the main use of this program
if (testtraj_flag) then

!New trajectories explore the phase space and
!discover new frames.
!If grid addition > 0, add these new frames
!to the specified grid
if (grid_addition > 0) then
    write(variable_length_text,FMT=FMT5_variable)&
            Ngrid_text_length
    write(Ngrid_text,FMT="(I0."//&
            trim(adjustl(variable_length_text))//")")&
            grid_addition
    gridpath2 = gridpath0//Ngrid_text//"/grid/"
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       INITIAL CONDITION DETERMINATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP SINGLE

!If we want to use the same initial conditions
!as an old set of trajectories open up whatever
!initial file that sample corresponds to
if ((useoldinitialbonddata_flag).and.&
    (.not.return_flag)) then
    open(trajectorieschannel,file=gridpath0//&
            initialbondfolder//intermediatefolder//initialfile)
    print *, ""
    print *, "     Deciding to use old initial conditions..."
    print *, ""

    !Check to see if we will run out of trajectories
    do n_testtraj = initial_n_testtraj, Ntesttraj
        read(trajectorieschannel,FMT=FMTinitial,iostat=iostate) &
                    ((INITIAL_BOND_DATA(j,i),j=1,6),i=1,Nbonds)

        !If so, we just exit out
        if (iostate /= 0) then
            print *, ""
            call itime(now)
            write(6,FMT=FMTnow) now
            print *, "Error! Not enough initial conditions in initialfile"
            print *, "Total number of initial conditions inside of it: ",&
                    n_testtraj-initial_n_testtraj
            print *, "Exiting analysis"
            print *, ""
       
            !Set the return flag true to
            !signal we want to quit
            return_flag = .true.
            exit
        end if
    end do

    close(trajectorieschannel)

    !And after we are done checking, open it up again
    !so we can read off initial conditions
    open(trajectorieschannel,file=gridpath0//&
            initialbondfolder//intermediatefolder//initialfile)
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!For bugtesting
if (makeCheckTrajectory_flag) then
    call InitialSampling3()

    testtrajDetailedRMSD_flag = .true.
    force_Neighbors = .true.

    subcellsearch_max1 = (/ 0, 0 /)
    subcellsearch_max2 = (/ 1, 1 /)

    call makeCheckTrajectoryGraphs
    call getRMSDDifferences1(gridpath0//prefix_text//"_PercentRMSDThreshold(1)")

    subcellsearch_max1 = (/ 0, 0 /)
    subcellsearch_max2 = (/ 2, 2 /)

    call makeCheckTrajectoryGraphs
    call getRMSDDifferences1(gridpath0//prefix_text//"_PercentRMSDThreshold(2)")

    subcellsearch_max1 = (/ 0, 0 /)
    subcellsearch_max2 = (/ 3, 3 /)

    call makeCheckTrajectoryGraphs
    call getRMSDDifferences1(gridpath0//prefix_text//"_PercentRMSDThreshold(3)")

    subcellsearch_max1 = (/ 0, 0 /)
    subcellsearch_max2 = (/ 4, 4 /)

    call makeCheckTrajectoryGraphs
    call getRMSDDifferences1(gridpath0//prefix_text//"_PercentRMSDThreshold(4)")

    subcellsearch_max1 = (/ 0, 0 /)
    subcellsearch_max2 = (/ 5, 5 /)

    call makeCheckTrajectoryGraphs
    call getRMSDDifferences1(gridpath0//prefix_text//"_PercentRMSDThreshold(5)")

    return_flag = .true.
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       TRAJECTORY CREATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP END SINGLE
!$OMP BARRIER

!If we are deciding to exit out, overflow
!initial_n_testtraj so that we don't make
!new trajectories
if (return_flag) then
    initial_n_testtraj = Ntesttraj + 1
else
    print *, ""
    call itime(now)
    write(6,FMT=FMTnow) now
    print *, "Trajectory Creation ... Start"
end if

!If we are checking traversal, allocate the necessary arrays
if (traversal_flag) allocate(traversal0(Ngrid_total,(var_bounds(1))**Nvar),&
                             traversal1(Ngrid_total,(var_bounds(2))**Nvar))

!If we are reading frames from a premade trajectory
!Open up the file to read the trajectories from
if (readtrajectory_flag) &
    open(trajectorieschannel,file=gridpath4//readtrajectoryfile)

!BEFORE going inside the parallel do loop, each thread should
!instantiate its own filechannels
allocate(filechannels(1+Ngrid_total))

!$OMP DO
do n_testtraj = initial_n_testtraj, Ntesttraj

    !If we're reusing old initial conditions,
    !that's easy; we just read off the file
    if (useoldinitialbonddata_flag) then
        !$OMP CRITICAL
        read(trajectorieschannel,&
                FMT=FMTinitial,iostat=iostate) &
                ((INITIAL_BOND_DATA(j,i),j=1,6),i=1,Nbonds)
        !$OMP END CRITICAL
    
    !Otherwise, we must create a random
    !initial trajectory from the ensemble
    !specified in PHYSICS and PARAMETERS
    else
        call InitialSampling3()
    end if

    !In most case, we will need to record the
    !initial conditions to a (potentially) new file
    !$OMP CRITICAL
    if ((.not.(useoldinitialbonddata_flag)).or.&
        (expfolder /= initialbondfolder)) then
        open(filechannel1,&
                file=gridpath5//initialfile,&
                position="append")
        write(filechannel1,FMTinitial)&
                ((INITIAL_BOND_DATA(j,i),j=1,6),i=1,Nbonds)
        close(filechannel1)
    end if
    !$OMP END CRITICAL

    !Each trajectory will have Ngrid_total
    !outputs; one for however many grids we use
    !The trajectory number will uniquely
    !identify one trajectory from another
    write(variable_length_text,FMT=FMT5_variable)&
            trajectory_text_length
    write(Ntraj_text,FMT="(I0."//&
            trim(adjustl(variable_length_text))//")")&
            n_testtraj
    
    !The grid number will uniquely identify one trajectory
    !Open all these files under each of their filechannels
    write(variable_length_text,FMT=FMT5_variable)&
            Ngrid_text_length
    filechannels(1) = 1000 + OMP_GET_THREAD_NUM()
    do Ngrid = 1, Ngrid_total
        write(Ngrid_text,FMT="(I0."//&
                trim(adjustl(variable_length_text))//")")&
                Ngrid

        !Our clever way to uniquely make a filechannel
        filechannels(1+Ngrid) = 1000 + 69 * Ngrid +&
                OMP_GET_THREAD_NUM()
        open(filechannels(1+Ngrid),file=gridpath5//&
                Ngrid_text//"_"//Ntraj_text//".dat")
    end do

    !Start timing to see how long each
    !trajectory takes to complete
    call system_clock(c1)
    call CPU_time(r1)
    print *, " Working on random new trajectory number: ",&
            Ntraj_text

    if (readtrajectory_flag) then
        read(trajectorieschannel,FMT=*) traj_index
        write(traj_index_text,FMT="(I0.6)") traj_index
        open(inputchannel,file=gridpath4//&
            readtrajectoryfolder//&
            traj_index_text,form="unformatted")
        call readTrajectory(&
                inputchannel,filechannels(1:1+Ngrid_total),&
                coords_initial,velocities_initial,&
                coords_final,velocities_final)

        close(inputchannel)
    else if (testtrajDetailedRMSD_flag) then
        call runTestTrajectory2(&
                filechannels(1:1+Ngrid_total),&
                coords_initial,velocities_initial,&
                coords_final,velocities_final)
    else

        !Then write the outputted RMSDS of
        !each trajectory onto those filechannels
        !Remark: checkMultipleGrids uses
        !filechannel1 to open files in the grid
        if (force_Permutations) then
            call runTrajectory_permute_cap(&
                    filechannels(1:1+Ngrid_total),&
                    coords_initial,velocities_initial,&
                    coords_final,velocities_final)
        else
            call runTrajectory_cap(&
                    filechannels(1:1+Ngrid_total),&
                    coords_initial,velocities_initial,&
                    coords_final,velocities_final)
        end if
!   call checkMultipleTrajectories(filechannels(1:1+Ngrid_total),&
!                              coords_initial,velocities_initial,&
!                              coords_final,velocities_final)
    end if

    call system_clock(c2)
    call CPU_time(r2)
    print *, " Trajectory "//Ntraj_text//&
            "     CPU Time: ", r2 - r1
    print *, " Trajectory "//Ntraj_text//&
            "    Wall Time: ", (c2 - c1) * system_clock_rate

    !Finally, close these files
    do Ngrid = 1, Ngrid_total
        close(filechannels(1+Ngrid))
    end do

    !$OMP CRITICAL
    !The only data that is recorded is the
    !first and last frame of the trajectory
    open(filechannel1,file=gridpath5//&
            timeslicefile,position="append")
    write(filechannel1,FMTtimeslice) &
                          ((coords_initial(i,j),i=1,3),j=1,Natoms),&
                          ((velocities_initial(i,j),i=1,3),j=1,Natoms),&
                          ((coords_final(i,j),i=1,3),j=1,Natoms),&
                          ((velocities_final(i,j),i=1,3),j=1,Natoms)
    close(filechannel1)

    !For traversal, we can see how many cells
    !the trajectory traversed merely by figuring
    !out how many nonzero values exist
    !(since 0s represent no traversal)
    if (traversal_flag) then
        open(filechannel1,file=gridpath5//&
                "traversal.dat",position="append")
        write(filechannel1,FMT=*) steps,&
                ((sum(traversal1(i,:))),i=1,Ngrid_total)
        close(filechannel1)
    end if

    !If we are doing a detailed test run then we
    !have a checkstate file that we can plot
    if (testtrajDetailedRMSD_flag) then

        !If there are eruptions (time rewinds) this gets
        !rid of them and puts the new checkstatefile
        !into truncatedcheckstatefile
        call processCheckstateFile()

        !This checks if the energy seems to be conserved
        call plotCheckstateEnergy("EnergyConservation"//Ntraj_text)

    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       DURING-CREATION TRAJECTORY ANALYSIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !If this flag is on, we are also checking the RMSD of a trajectory per timestep
    !This is very expensive so please don't turn this on
    !(and you have to manually turn it on in checkMultipleTrajectories anyway)
    if (testtrajRMSD_flag) then
        open(gnuplotchannel,file=gridpath5//gnuplotfile)
        write(gnuplotchannel,*) 'set term pngcairo size 600, 600'
        write(gnuplotchannel,*) 'set output "'//gridpath0//Ngrid_text//'/RMSD'//&
                                                reject_text//Ntraj_text//'.png"'
        write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
        write(gnuplotchannel,*) 'unset key'
        write(gnuplotchannel,*) 'bin_width = 0.001'
        write(gnuplotchannel,*) 'bin_number(x) = floor(x/bin_width)'
        write(gnuplotchannel,*) 'rounded(x) = bin_width * (bin_number(x) + 0.5)'
        write(gnuplotchannel,*) 'set xlabel "RMSD"'
        write(gnuplotchannel,*) 'set ylabel "Occurence"'
        write(gnuplotchannel,*) 'plot "'//gridpath4//Ngrid_text//"/"//&
                                Ntraj_text//'.dat'//&
                                '" u (rounded($1)):(1.0) smooth frequency with boxes'
        close(gnuplotchannel)
        call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)
    end if
    !$OMP END CRITICAL
    
end do
!$OMP END DO
print *, ""

if (readtrajectory_flag) close(trajectorieschannel)

deallocate(filechannels)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       TRAJECTORY CREATION FLAG END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Close the initial conditions file if
!we are reusing old ones
!$OMP SINGLE
if (useoldinitialbonddata_flag) close(trajectorieschannel)
!$OMP END SINGLE

!This closes that big, enclosing if
!statement on whether to make trajectories
end if

!$OMP END PARALLEL

!Now that we are out of the parallel section
!we can exit out if there was a problem
if (return_flag) then
    return
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       POST-CREATION INTERPOLATION ANALYSIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!If this is not a comparison analysis then
!most likely we just made some new trajectories
!that we want to check
if (.not.comparison_flag) then

        call processCheckstateFile()

if (testtrajDetailedRMSD_flag) then
    call getRMSDDifferences1(&
    gridpath4//"BinaryRMSDCheck")

    call plotEnergyConservationInformatics()
end if

if ((testtrajDetailedRMSD_flag).and.&
    (gather_interpolation_flag)) &
    call getAlphaErrorDistribution(&
    gridpath4//"AlphaErrorDistribution")

!But, note that we only have interpolations
!to check if we actually interpolated AND
!gathered data
if ((gather_interpolation_flag)) then

    !First process the raw data; this is
    !required for all the next steps and for
    !future comparisons
    call processInterpolationFile2()

    if (interpolation_flag) then

    call itime(now)
    write(6,FMT=FMTnow) now
    print *, "   Making plot: ", "TDDRED"
    print *, ""

    !For each interpolation analysis we need a
    !temporary comparison file that specifies
    !what we are analyzing

    !First is analysis of all of the data
    call system("rm "//gridpath0//comparison_file)
    call system('echo "0.0 0.0 100.0 100.0" >>'//gridpath0//comparison_file)
    call getRMSDinterpolation2("TDDRED_All")

    !Then, a region that is localised to what
    !I call "collision" or around the 3.0,
    !3.5 region for a var1,var2 system
    call system("rm "//gridpath0//comparison_file)
    call system('echo "2.990 3.780 0.55 0.55" >>'//gridpath0//comparison_file)
    call getRMSDinterpolation2("TDDRED_Collision")

    !And we also take a look at all data
    !that is outside of this region
    call system("rm "//gridpath0//comparison_file)
    call system('echo "2.990 3.780 -0.55 -0.55" >>'//gridpath0//comparison_file)
    call getRMSDinterpolation2("TDDRED_NonCollision")

    end if

    call itime(now)
    write(6,FMT=FMTnow) now
    print *, "   Making plot: ", "InterpolationOccurenceHeatmap"
    print *, ""

    !Finally, a heatmap of all these interpolation
    !occurences (with some restrictions)
    call plotInterpolationOccurenceHeatmap()
end if

!All other interpolation analysis is done as
!a comparison (if wanted, you can choose to
!compare only a single experiment)
else

call getRMSDErrorPlots("RMSDvsError")

!If this truly is an interpolation analysis then
!the word "interpolation" should be in the input
i = INDEX(comparison_SATRVname,"Interpolation")+13

!Currently, seven different kinds of analysis
!can be done
if ((i > 13).and.(&
    (trim(adjustl(comparison_SATRVname(i:)))=="TDD").or.&
    (trim(adjustl(comparison_SATRVname(i:)))=="ARD").or.&
    (trim(adjustl(comparison_SATRVname(i:)))=="IRD").or.&
    (trim(adjustl(comparison_SATRVname(i:)))=="RSV1D").or.&
    (trim(adjustl(comparison_SATRVname(i:)))=="RSV2D").or.&
    (trim(adjustl(comparison_SATRVname(i:)))=="AED").or.&
    (trim(adjustl(comparison_SATRVname(i:)))=="IED").or.&
    (trim(adjustl(comparison_SATRVname(i:)))=="RED")&
             )) then

    call itime(now)
    write(6,FMT=FMTnow) now
    print *, "   Making plot: ", trim(adjustl(comparison_SATRVname(i:)))
    print *, ""

    !This assumes the interpolation
    !processing was done already (at
    !the end of the experiment)
    call getInterpolationplot(trim(adjustl(comparison_SATRVname(i:))))

    print *, ""
    call itime(now)
    write(6,FMT=FMTnow) now
    print *, "Successfully exited analysis"
    print *, ""

    !Comparisons were designed not
    !to do any other analysis
    return
end if
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       POST-CREATION LIBRARY ANALYSIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Specify how many trajectories we will be analyzing
!We communicate to analyzeScatteringAngles and others
!through the shared variable Ntraj
Ntraj = Ntesttraj

!Initialize each maximum and minimum value
!These are also shared among modules
max_TranslationalEnergy = 0.0d0
max_absenergychange = 0.0
min_absenergychange = 1.0e9
max_relenergychange = 0.0
min_relenergychange = 1.0e9
max_rotenergychange = 0.0
min_rotenergychange = 1.0e9

!If we have a list of files we need to compare,
!then we need to comb through those
if (comparison_flag) then

    !For non-interpolation analyses, single
    !comparison cannot be done
    !(This may be deprecated)
    if (comparison_number < 1) then
        print *, ""
        call itime(now)
        write(6,FMT=FMTnow) now
        print *, "Number of trajectory sets being compared is inadequate"
        print *, "Terribly exited analysis"
        print *, ""

        return
    end if

    !The files with the desired data are
    !stored in allprefixes and the length
    !of each file is stored in alllengths

    !They must each individually be
    !processed to figure out what the
    !upper and lower bounds are
    call postProcess(allprefixes(1:alllengths(1))//&
            intermediatefolder)
    do i = 1, comparison_number-1
        call postProcess(allprefixes(1+sum(alllengths(1:i)):sum(alllengths(1:i+1)))//&
                intermediatefolder)
    end do

!Otherwise, we do not need to process
!all the files in the library
else
    !We only need to look at the files
    !if we are actually analyzing something
    if (percentthreshold_flag .or. testtrajSA_flag .or. testheatmapSA_flag) then
        call postProcess(expfolder//intermediatefolder)

    !Otherwise, we can just exit early
    else
        print *, ""
        call itime(now)
        write(6,FMT=FMTnow) now
        print *, "Successfully exited analysis"
        print *, ""
   
        return
    end if
end if

!We need to redo the postProcess (and
!convergence) of each grid as well in
!case the maximum and minimum of the
!trajectories just processed are different
if (comparison_flag .or. trueSA_flag .or.&
    testtrajSA_flag .or. testheatmapSA_flag) then

    call itime(now)
    write(6,FMT=FMTnow) now
    print *, "   Making plots: ", "Initial_SATRVDistribution"
    print *, "   Making plots: ", "Initial_HeatMap_SATRVDistribution"
    print *, ""

    old_filename = ""
    do Ngrid = 1, Ngrid_max
        !Prepare the name of the file; here we are
        !uniquely identifying the file by how many
        !trajectories we've accumulated and put in it
        Ntraj = Ngrid*Ntraj_max
        write(variable_length_text,FMT=FMT5_variable)&
                Ngrid_text_length
        write(Ngrid_text,FMT="(I0."//&
                trim(adjustl(variable_length_text))//")")&
                Ngrid
        write(variable_length_text,FMT=FMT5_variable)&
                trajectory_text_length
        write(Ntraj_text,FMT="(I0."//&
                trim(adjustl(variable_length_text))//")")&
                Ntraj

        !Post process the file; here we are telling
        !the subroutine that there are Ntraj_max
        !trajectories in the file
        Ntraj = Ntraj_max
        call postProcess(Ngrid_text//"/")

        !We create this unique file that has all
        !trajectories up until that grid so beforehand
        !we need to concatenate all the old
        !trajectories with the new file
        if (trim(adjustl(old_filename)) /= "") then
            new_filename = gridpath5//Ntraj_text//SATRVfile
            call system("cat "//trim(adjustl(old_filename))//" "//&
                        gridpath0//Ngrid_text//"/"//SATRVfile//" > "//new_filename)
            old_filename = new_filename

        !If this is the first grid, then we
        !just copy the new file to this unique file
        else
            old_filename = gridpath5//Ntraj_text//SATRVfile
            call system("cp "//gridpath0//Ngrid_text//"/"//SATRVfile//" "//trim(adjustl(old_filename)))
        end if

        !Then bin the data
        Ntraj = Ngrid*Ntraj_max
        call getScatteringAngles1(&
                Ntraj_text,"SATRVDistribution")
!               call getScatteringAngles1(&
!                       expfolder//intermediatefolder//Ntraj_text,&
!                       "SATRVDistribution")
    end do

    !Finally, we look at the scattering angle
    !and energy change distributions of the grid
    !with the maximums and minimums gathered from above
    !and see whether they converge;
    !these plots will be used for reference
    Ntraj = Ntesttraj
    write(variable_length_text,FMT=FMT5_variable)&
            Ngrid_text_length
    write(Ngrid_text,FMT="(I0."//&
            trim(adjustl(variable_length_text))//")")&
            Ngrid_total
    write(variable_length_text,FMT=FMT5_variable)&
            trajectory_text_length
    write(Ntraj_text,FMT="(I0."//&
            trim(adjustl(variable_length_text))//")")&
            Ntraj
    
    !And finally, make the convergence plots
    !with these lower and upper bounds and
    !binned data
    call getConvergenceImage(0.0,real(pi),&
            1, "ScatteringAngle")
    call getConvergenceImage(&
            min_absenergychange,max_absenergychange,&
            3, "AbsoluteEnergyChange")
    call getConvergenceImage(&
            min_relenergychange,max_relenergychange,&
            4, "RelativeEnergyChange")
    call getConvergenceImage(&
            min_rotenergychange,max_rotenergychange,&
            5, "RotationalEnergyChange")

end if

!In a comparison run, we are interested
!in a specific type of data like the
!scattering angle or whatsoever
if (comparison_flag) then

    !The name the user gives us tells us
    !  1) which column of the SATRV file to look in
    !  2) how to name our output file
    if (trim(adjustl(comparison_SATRVname)) == "ScatteringAngle") then
        comparison_SATRVcolumn = 1
        lowerlimit = 0.0
        upperlimit = real(pi)
    else if (trim(adjustl(comparison_SATRVname)) == "AbsoluteEnergyChange") then
        comparison_SATRVcolumn = 3
        lowerlimit = min_absenergychange
        upperlimit = max_absenergychange
    else if (trim(adjustl(comparison_SATRVname)) == "RelativeEnergyChange") then
        comparison_SATRVcolumn = 4
        lowerlimit = min_relenergychange
        upperlimit = max_relenergychange
    else if (trim(adjustl(comparison_SATRVname)) == "RotationalEnergyChange") then
        comparison_SATRVcolumn = 5
        lowerlimit = min_rotenergychange
        upperlimit = max_rotenergychange
    else
    end if

    !If need be, use the upper and lower
    !bounds dictated by the user
    if (comparison_upperlimit /= comparison_lowerlimit) then
        call getConvergenceImage(comparison_lowerlimit,comparison_upperlimit,&
                comparison_SATRVcolumn,trim(adjustl(comparison_SATRVname)))
    end if

    call itime(now)
    write(6,FMT=FMTnow) now
    print *, "   Making plot: ", "Comparison_"//&
            trim(adjustl(Ntraj_text))//trim(adjustl(comparison_SATRVname))
    print *, ""

    !We now need to bin each set of trajectories we
    !want to compare (done in getScatteringAngles1)
    call system("cp "//gridpath0//allprefixes(1:alllengths(1))//&
                intermediatefolder//SATRVfile//&
                " "//gridpath5//allprefixes(1:alllengths(1)-1)//SATRVfile)
    call getScatteringAngles1(allprefixes(1:alllengths(1)-1),&
            allprefixes(1:alllengths(1))//"SATRVdistribution")
    do i = 1, comparison_number-1
        call system("cp "//gridpath0//allprefixes(&
                1+sum(alllengths(1:i)):&
                sum(alllengths(1:i+1)))//&
                intermediatefolder//SATRVfile//&
                " "//gridpath5//allprefixes(&
                1+sum(alllengths(1:i)):&
                sum(alllengths(1:i+1))-1)//SATRVfile)
        call getScatteringAngles1(allprefixes(&
                1+sum(alllengths(1:i)):&
                sum(alllengths(1:i+1))-1),&
                                  allprefixes(&
                1+sum(alllengths(1:i)):&
                sum(alllengths(1:i+1))-1)//&
                "SATRVdistribution")
    end do

    !If the upper and lower limit are equal, that is
    !an internal sign that means that the user wants
    !to use the natural bounds of the distribution
    if (comparison_upperlimit /= comparison_lowerlimit) then
        call getComparedScatteringAngles(&
                 comparison_lowerlimit,comparison_upperlimit,&
                 "Comparison_"//trim(adjustl(Ntraj_text))//&
                 trim(adjustl(comparison_SATRVname)),&
                 comparison_SATRVcolumn,trim(adjustl(comparison_SATRVname)))

    !Otherwise, that means we must use the user-defined bounds
    else
        call getComparedScatteringAngles(lowerlimit,upperlimit,&
                 "Comparison_"//trim(adjustl(Ntraj_text))//&
                 trim(adjustl(comparison_SATRVname)),&
                 comparison_SATRVcolumn,trim(adjustl(comparison_SATRVname)))
    end if

    print *, ""
    call itime(now)
    write(6,FMT=FMTnow) now
    print *, "Successfully exited analysis"
    print *, ""

    return
end if



!In a run without comparison, usually we do
!some analysis on the trajectories afterwards
!The following calls do some data
!manipulation and make some figures

if (percentthreshold_flag) then
    call itime(now)
    write(6,FMT=FMTnow) now
    print *, "   Making plot: ", expfolder//"PercentRMSDThreshold"
    print *, ""

    !This subroutine takes an additional argument
    !that specifies (out of all the grids) which
    !ones to actually plot
    if (percentthreshold_key < 1) then
        call getRMSDThresholds1(expfolder//intermediatefolder,&
                "PercentRMSDThreshold")
    else
        call getRMSDThresholds1(expfolder//intermediatefolder,&
                "PercentRMSDThreshold",percentthreshold_key)
    end if
end if

if (testheatmapSA_flag) then
    call itime(now)
    write(6,FMT=FMTnow) now
    print *, "   Making plot: ", expfolder//"SATRVDistribution"
    print *, "   Making plot: ", expfolder//"HeatMap_SATRVDistribution"
    print *, ""

    call getScatteringAngles1("","SATRVDistribution")
end if

if (testtrajSA_flag) then
    call itime(now)
    write(6,FMT=FMTnow) now
    print *, "   Making plot: ", expfolder//"Final_SATRVDistribution"
    print *, ""
    
    call getScatteringAngles2(expfolder//intermediatefolder,&
            "Final_SATRVDistribution")
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print *, ""
call itime(now)
write(6,FMT=FMTnow) now
print *, "Successfully exited analysis"
print *, ""

return

end program checkNewTrajectorieswithMultipleGrids






subroutine makeCheckTrajectoryGraphs
use runTrajectory
use PARAMETERS
use FUNCTIONS
use VARIABLES
use ANALYSIS
use PHYSICS
use analyzeScatteringAngleswithMultipleGrids
use analyzeHeatMapswithMultipleGrids
implicit none

!Grid Directory/File Formatting Strings
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(6) :: reject_text
character(12) :: prefix_text
integer,allocatable :: filechannels(:)

!Trajectory Variables
integer :: iostate
integer :: order, neighbor_check, number_of_frames
real(dp) :: U, KE
real(dp) :: min_rmsd, min_rmsd_prime
real(dp),dimension(Nvar) :: vals

!Trajectory Bounds
integer :: trajectory_total_frames
integer :: trajectory_max_frames
integer :: trajectory_max_neighbor_check
real(dp) :: trajectory_min_rmsd
real(dp) :: trajectory_max_var1, trajectory_min_var1
real(dp) :: trajectory_max_var2, trajectory_min_var2
real(dp) :: trajectory_CPU_time,trajectory_wall_time

!Trajectory Output
real(dp),dimension(3,Natoms) :: coords_initial, velocities_initial
real(dp),dimension(3,Natoms) :: coords_final, velocities_final

!Timing related
integer :: OMP_GET_THREAD_NUM
real :: system_clock_rate
integer :: cr
integer :: trajectory_t0, trajectory_t1
real :: checktrajectory_wall_time
character(10) :: checktrajectory_wall_time_text
integer,dimension(3) :: now

!Initialize the clock
call system_clock(count_rate=cr)
system_clock_rate = 1.0/real(cr)

!The grid number will uniquely identify one trajectory
!Open all these files under filechannels
allocate(filechannels(1+Ngrid_total))
write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
filechannels(1) = 1000 + OMP_GET_THREAD_NUM()

do Ngrid = 1, Ngrid_total
        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
        filechannels(1+Ngrid) = 1000 + 69 * Ngrid + OMP_GET_THREAD_NUM()
        open(filechannels(1+Ngrid),file=gridpath4//Ngrid_text//"/dump.dat")
end do

reject_flag = (.not.(reject_flag))
call getPrefixText(prefix_text)
reject_text = prefix_text(1:6)

!Next, run the trajectory
call system_clock(trajectory_t0)
call checkMultipleTrajectories(filechannels(1:1+Ngrid_total),&
                               coords_initial,velocities_initial,&
                               coords_final,velocities_final)
call system_clock(trajectory_t1)
checktrajectory_wall_time = (trajectory_t1-trajectory_t0)*system_clock_rate

!To make a nice plot, we need to know various bounds of our data
!so the maximum and minimum are found by reading through each line
!of the files
trajectory_max_var1 = 0.0d0
trajectory_min_var1 = 1.0d9
trajectory_max_var2 = 0.0d0
trajectory_min_var2 = 1.0d9

trajectory_total_frames = 0
trajectory_max_frames = 0
trajectory_max_neighbor_check = 0
trajectory_min_rmsd = 1.0d9
open(filechannel1,file=gridpath0//intermediatefolder//checkstatefile)
do
        read(filechannel1,FMT=*,iostat=iostate) number_of_frames,order,neighbor_check,steps,&
                                                min_rmsd,min_rmsd_prime,vals(1),vals(2),U,KE
        if (iostate /= 0) exit
        trajectory_total_frames = trajectory_total_frames + 1

        trajectory_max_var1 = max(trajectory_max_var1,vals(1))
        trajectory_min_var1 = min(trajectory_min_var1,vals(1))
        trajectory_max_var2 = max(trajectory_max_var2,vals(2))
        trajectory_min_var2 = min(trajectory_min_var2,vals(2))

        trajectory_min_rmsd = min(trajectory_min_rmsd,min_rmsd)
        trajectory_max_frames = max(trajectory_max_frames,number_of_frames)
        trajectory_max_neighbor_check = max(trajectory_max_neighbor_check,&
                                            neighbor_check)
end do
close(filechannel1)

print *, ""
call itime(now)
write(6,FMT=FMTnow) now
print *, '   Making plot: "'//gridpath0//'checkTrajectory_'//reject_text//'.png"'
print *, ""

!Finally, plot the data obtained from this trajectory
open(gnuplotchannel,file=gridpath0//intermediatefolder//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath0//'checkTrajectory_'//reject_text//'.png"'
write(gnuplotchannel,*) 'set style line 1 lc rgb "red" pt 5'
write(gnuplotchannel,*) 'set style line 2 lc rgb "green" pt 7'
write(gnuplotchannel,*) 'set style line 3 lc rgb "blue" pt 13'
write(gnuplotchannel,*) 'set style line 4 lc rgb "orange" pt 9'
write(gnuplotchannel,*) 'set style line 5 lc rgb "yellow" pt 11'
write(gnuplotchannel,*) 'set style line 6 lc rgb "pink" pt 20'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set multiplot layout 6,1 margins 0.10,0.95,.1,.95 spacing 0,0 title '//&
        '"Test of Trajectory '//trim(adjustl(checkstateTrajectory))//' Against the Grid" font ",18" offset 0,3'
write(checktrajectory_wall_time_text,FMT="(F10.2)") checktrajectory_wall_time
write(gnuplotchannel,*) 'set label 3 "Total Wall Time: '//checktrajectory_wall_time_text//&
                                                 ' s" at screen 0.6, 0.940'

write(gnuplotchannel,*) 'max_var1 = ', trajectory_max_var1
write(gnuplotchannel,*) 'min_var1 = ', trajectory_min_var1
write(gnuplotchannel,*) 'max_var1 = ceil(max_var1)'
write(gnuplotchannel,*) 'min_var1 = floor(min_var1)'
write(gnuplotchannel,*) 'delta_var1 = (max_var1 - min_var1) / 4'
write(gnuplotchannel,*) 'max_var2 = ', trajectory_max_var2
write(gnuplotchannel,*) 'min_var2 = ', trajectory_min_var2
write(gnuplotchannel,*) 'max_var2 = ceil(max_var2)'
write(gnuplotchannel,*) 'min_var2 = floor(min_var2)'
write(gnuplotchannel,*) 'delta_var1 = (max_var1 - min_var1) / 4'
write(gnuplotchannel,*) 'delta_var2 = (max_var2 - min_var2) / 4'
write(gnuplotchannel,*) 'max_frames = ', trajectory_max_frames
write(gnuplotchannel,*) 'max_steps = ', trajectory_total_frames
write(gnuplotchannel,*) 'max_neighbor_check = ', max(trajectory_max_neighbor_check,4)
write(gnuplotchannel,*) 'min_rmsd = ', trajectory_min_rmsd

write(gnuplotchannel,*) 'steps_scaling = 1000'
write(gnuplotchannel,*) 'set xrange [0:max_steps/steps_scaling]'
write(gnuplotchannel,*) 'set xtics 0, 1, max_steps/steps_scaling'
write(gnuplotchannel,*) 'set format x ""'
write(gnuplotchannel,*) 'unset xlabel'

write(gnuplotchannel,*) 'set ylabel "Var1 (A)"'
write(gnuplotchannel,*) 'set yrange [min_var1-delta_var1*.25 : max_var1+delta_var1*.25]'
write(gnuplotchannel,*) 'set ytics min_var1, delta_var1, max_var2'
write(gnuplotchannel,*) 'plot "'//gridpath0//intermediatefolder//checkstatefile//&
                        '" u (($4)/steps_scaling):7 w lines'

write(gnuplotchannel,*) 'unset label 3'

write(gnuplotchannel,*) 'set ylabel "Var2 (A)"'
write(gnuplotchannel,*) 'set yrange [min_var2-delta_var2*.25 : max_var2+delta_var2*.25]'
write(gnuplotchannel,*) 'set ytics min_var2, delta_var2, max_var2'
write(gnuplotchannel,*) 'plot "'//gridpath0//intermediatefolder//checkstatefile//&
                        '" u (($4)/steps_scaling):8 w lines'

write(gnuplotchannel,*) 'set ylabel "Number of\nFrames Checked"'
write(gnuplotchannel,*) 'set yrange [-max_frames*.10:max_frames+max_frames*.10]'
write(gnuplotchannel,*) 'set ytics 0, floor(max_frames*.25), max_frames'
write(gnuplotchannel,*) 'plot "'//gridpath0//intermediatefolder//checkstatefile//&
                        '" u (($4)/steps_scaling):1 w points'

write(gnuplotchannel,*) 'set ylabel "Order of Cells\nChecked"'
write(gnuplotchannel,*) 'set yrange [-0.25:2.25]'
write(gnuplotchannel,*) 'set ytics 1'
write(gnuplotchannel,*) 'plot "'//gridpath0//intermediatefolder//checkstatefile//&
                        '" u (($4)/steps_scaling):2 w points'

write(gnuplotchannel,*) 'set ylabel "Number of Cells\nChecked"'
write(gnuplotchannel,*) 'set yrange [0:max_neighbor_check+max_neighbor_check*.10]'
write(gnuplotchannel,*) 'set ytics 1, floor(max_neighbor_check*.25), max_neighbor_check'
write(gnuplotchannel,*) 'plot "'//gridpath0//intermediatefolder//checkstatefile//&
                        '" u (($4)/steps_scaling):3 w points'

write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set format x'
write(gnuplotchannel,*) 'set xlabel "Timesteps (Thousands)"'

write(gnuplotchannel,*) 'set yrange [min_rmsd*0.9:',default_rmsd,']'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,*) 'set ylabel "Timestep\nRMSD (A)"'
write(gnuplotchannel,*) 'set ytics ("5e-1" .5, "5e-2" .05, "5e-3" .005,'//&
                                   '"5e-4" .0005, "5e-5" .00005)'
write(gnuplotchannel,*) 'plot "'//gridpath0//intermediatefolder//checkstatefile//&
                        '" u (($4)/steps_scaling):($6==$5?$5:1/0) w points lc rgb "blue",\'
write(gnuplotchannel,*) '     "'//gridpath0//intermediatefolder//checkstatefile//&
                        '" u (($4)/steps_scaling):($6>$5?$5:1/0) w points lc rgb "red",\'
write(gnuplotchannel,*) '     "'//gridpath0//intermediatefolder//checkstatefile//&
                        '" u (($4)/steps_scaling):($6<$5?$6:1/0) w points lc rgb "green"'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath0//intermediatefolder//gnuplotfile)

reject_flag = (.not.(reject_flag))
call getPrefixText(prefix_text)
reject_text = prefix_text(1:6)

!Next, run the trajectory
call system_clock(trajectory_t0)
call checkMultipleTrajectories(filechannels(1:1+Ngrid_total),&
                               coords_initial,velocities_initial,&
                               coords_final,velocities_final)
call system_clock(trajectory_t1)
checktrajectory_wall_time = (trajectory_t1-trajectory_t0)*system_clock_rate

!To make a nice plot, we need to know various bounds of our data
!so the maximum and minimum are found by reading through each line
!of the files
trajectory_max_var1 = 0.0d0
trajectory_min_var1 = 1.0d9
trajectory_max_var2 = 0.0d0
trajectory_min_var2 = 1.0d9

trajectory_total_frames = 0
trajectory_max_frames = 0
trajectory_max_neighbor_check = 0
trajectory_min_rmsd = 1.0d9
open(filechannel1,file=gridpath0//intermediatefolder//checkstatefile)
do
        read(filechannel1,FMT=*,iostat=iostate) number_of_frames,order,neighbor_check,steps,&
                                                min_rmsd,min_rmsd_prime,vals(1),vals(2),U,KE
        if (iostate /= 0) exit
        trajectory_total_frames = trajectory_total_frames + 1

        trajectory_max_var1 = max(trajectory_max_var1,vals(1))
        trajectory_min_var1 = min(trajectory_min_var1,vals(1))
        trajectory_max_var2 = max(trajectory_max_var2,vals(2))
        trajectory_min_var2 = min(trajectory_min_var2,vals(2))

        trajectory_min_rmsd = min(trajectory_min_rmsd,min_rmsd)
        trajectory_max_frames = max(trajectory_max_frames,number_of_frames)
        trajectory_max_neighbor_check = max(trajectory_max_neighbor_check,&
                                            neighbor_check)
end do
close(filechannel1)

print *, ""
call itime(now)
write(6,FMT=FMTnow) now
print *, '   Making plot: "'//gridpath0//'checkTrajectory_'//reject_text//'.png"'
print *, ""

!Finally, plot the data obtained from this trajectory
open(gnuplotchannel,file=gridpath0//intermediatefolder//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath0//'checkTrajectory_'//reject_text//'.png"'
write(gnuplotchannel,*) 'set style line 1 lc rgb "red" pt 5'
write(gnuplotchannel,*) 'set style line 2 lc rgb "green" pt 7'
write(gnuplotchannel,*) 'set style line 3 lc rgb "blue" pt 13'
write(gnuplotchannel,*) 'set style line 4 lc rgb "orange" pt 9'
write(gnuplotchannel,*) 'set style line 5 lc rgb "yellow" pt 11'
write(gnuplotchannel,*) 'set style line 6 lc rgb "pink" pt 20'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set multiplot layout 6,1 margins 0.10,0.95,.1,.95 spacing 0,0 title '//&
        '"Test of Trajectory '//trim(adjustl(checkstateTrajectory))//' Against the Grid" font ",18" offset 0,3'
write(checktrajectory_wall_time_text,FMT="(F10.2)") checktrajectory_wall_time
write(gnuplotchannel,*) 'set label 3 "Total Wall Time: '//checktrajectory_wall_time_text//&
                                                 ' s" at screen 0.6, 0.940'

write(gnuplotchannel,*) 'max_var1 = ', trajectory_max_var1
write(gnuplotchannel,*) 'min_var1 = ', trajectory_min_var1
write(gnuplotchannel,*) 'max_var1 = ceil(max_var1)'
write(gnuplotchannel,*) 'min_var1 = floor(min_var1)'
write(gnuplotchannel,*) 'delta_var1 = (max_var1 - min_var1) / 4'
write(gnuplotchannel,*) 'max_var2 = ', trajectory_max_var2
write(gnuplotchannel,*) 'min_var2 = ', trajectory_min_var2
write(gnuplotchannel,*) 'max_var2 = ceil(max_var2)'
write(gnuplotchannel,*) 'min_var2 = floor(min_var2)'
write(gnuplotchannel,*) 'delta_var1 = (max_var1 - min_var1) / 4'
write(gnuplotchannel,*) 'delta_var2 = (max_var2 - min_var2) / 4'
write(gnuplotchannel,*) 'max_frames = ', trajectory_max_frames
write(gnuplotchannel,*) 'max_steps = ', trajectory_total_frames
write(gnuplotchannel,*) 'max_neighbor_check = ', max(trajectory_max_neighbor_check,4)
write(gnuplotchannel,*) 'min_rmsd = ', trajectory_min_rmsd

write(gnuplotchannel,*) 'steps_scaling = 1000'
write(gnuplotchannel,*) 'set xrange [0:max_steps/steps_scaling]'
write(gnuplotchannel,*) 'set xtics 0, 1, max_steps/steps_scaling'
write(gnuplotchannel,*) 'set format x ""'
write(gnuplotchannel,*) 'unset xlabel'

write(gnuplotchannel,*) 'set ylabel "Var1 (A)"'
write(gnuplotchannel,*) 'set yrange [min_var1-delta_var1*.25 : max_var1+delta_var1*.25]'
write(gnuplotchannel,*) 'set ytics min_var1, delta_var1, max_var2'
write(gnuplotchannel,*) 'plot "'//gridpath0//intermediatefolder//checkstatefile//&
                        '" u (($4)/steps_scaling):7 w lines'

write(gnuplotchannel,*) 'unset label 3'

write(gnuplotchannel,*) 'set ylabel "Var2 (A)"'
write(gnuplotchannel,*) 'set yrange [min_var2-delta_var2*.25 : max_var2+delta_var2*.25]'
write(gnuplotchannel,*) 'set ytics min_var2, delta_var2, max_var2'
write(gnuplotchannel,*) 'plot "'//gridpath0//intermediatefolder//checkstatefile//&
                        '" u (($4)/steps_scaling):8 w lines'

write(gnuplotchannel,*) 'set ylabel "Number of\nFrames Checked"'
write(gnuplotchannel,*) 'set yrange [-max_frames*.10:max_frames+max_frames*.10]'
write(gnuplotchannel,*) 'set ytics 0, floor(max_frames*.25), max_frames'
write(gnuplotchannel,*) 'plot "'//gridpath0//intermediatefolder//checkstatefile//&
                        '" u (($4)/steps_scaling):1 w points'

write(gnuplotchannel,*) 'set ylabel "Order of Cells\nChecked"'
write(gnuplotchannel,*) 'set yrange [-0.25:2.25]'
write(gnuplotchannel,*) 'set ytics 1'
write(gnuplotchannel,*) 'plot "'//gridpath0//intermediatefolder//checkstatefile//&
                        '" u (($4)/steps_scaling):2 w points'

write(gnuplotchannel,*) 'set ylabel "Number of Cells\nChecked"'
write(gnuplotchannel,*) 'set yrange [0:max_neighbor_check+max_neighbor_check*.10]'
write(gnuplotchannel,*) 'set ytics 1, floor(max_neighbor_check*.25), max_neighbor_check'
write(gnuplotchannel,*) 'plot "'//gridpath0//intermediatefolder//checkstatefile//&
                        '" u (($4)/steps_scaling):3 w points'

write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set format x'
write(gnuplotchannel,*) 'set xlabel "Timesteps (Thousands)"'

write(gnuplotchannel,*) 'set yrange [min_rmsd*0.9:',default_rmsd,']'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,*) 'set ylabel "Timestep\nRMSD (A)"'
write(gnuplotchannel,*) 'set ytics ("5e-1" .5, "5e-2" .05, "5e-3" .005,'//&
                                   '"5e-4" .0005, "5e-5" .00005)'
write(gnuplotchannel,*) 'plot "'//gridpath0//intermediatefolder//checkstatefile//&
                        '" u (($4)/steps_scaling):($6==$5?$5:1/0) w points lc rgb "blue",\'
write(gnuplotchannel,*) '     "'//gridpath0//intermediatefolder//checkstatefile//&
                        '" u (($4)/steps_scaling):($6>$5?$5:1/0) w points lc rgb "red",\'
write(gnuplotchannel,*) '     "'//gridpath0//intermediatefolder//checkstatefile//&
                        '" u (($4)/steps_scaling):($6<$5?$6:1/0) w points lc rgb "green"'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath0//intermediatefolder//gnuplotfile)

deallocate(filechannels)
return

end subroutine makeCheckTrajectoryGraphs
