!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	PROGRAM
!		makeGridwithNewTrajectories
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	PURPOSE
!		The program creates a library (multiple grids) according to PARAMETERS
!		with MD simulations produced by runTrajectory, governed by PHYSICS
!		and with collective variables calculated according to VARIABLES;
!		each frame is added according to interactSingleGrid
!
!		The program may take a long time so there can be occasional
!		checks on a grid; these are governed by ANALYSIS
!		and the frequency of checks is governed by PARAMETERS;
!		the MD simulations used to check the grid are also governed
!		by runTrajectory which checks the grid with interactSingleGrid
!
!		When a grid is completed, its scattering angle plot is made
!		according to analyzeScatteringAngleswithMultipleGrids
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	FILECHANNELS			ACTION
!
!		GNUPLOTCHANNEL			OPEN, WRITE, CLOSE
!		PROGRESSCHANNEL			OPEN, WRITE, CLOSE
!		FILECHANNEL1			OPEN, WRITE, CLOSE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       CALLS                           MODULE
!
!		SYSTEM				INTRINSIC
!		CPU_TIME			INTRINSIC
!		SYSTEMCLOCK			INTRINSIC
!
!		addTrajectory			runTrajectory
!		checkTrajectory			runTrajectory
!
!		getScatteringAngles2		analyzeScatteringAngleswithMultipleGrids
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	OUTPUT				FILETYPE
!
!		ALOT				WOAH
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program makeGridwithNewTrajectories
use interactMultipleGrids
use runTrajectory
use PARAMETERS
use FUNCTIONS
use VARIABLES
use ANALYSIS
use PHYSICS
use analyzeScatteringAngleswithMultipleGrids
use analyzeRMSDThresholdwithMultipleGrids
use analyzeHeatMapswithMultipleGrids
implicit none

!Grid Directory/File Formatting Strings
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(trajectory_text_length) :: Ntraj_text
character(12) :: prefix_text

integer :: traj_index
character(6) :: traj_index_text
integer :: inputchannel = 307

character(15) :: vals_interpolation_text

!Trajectory Variables
real(dp) :: trajectory_CPU_time,trajectory_wall_time
real(dp) :: scattering_angle
real(dp),dimension(3) :: TRVenergies1,TRVenergies2,dTRVenergies
real(dp) :: totalEnergy
integer :: header1_old,header2_old,header3_old
integer :: max_header1_delta, max_header2_delta, max_header3_delta
integer,dimension(Norder_max) :: max_headers_delta
integer,allocatable :: filechannels(:)

!Trajectory Output
real(dp),dimension(3,Natoms) :: coords_initial, velocities_initial
real(dp),dimension(3,Natoms) :: coords_final, velocities_final

!For Consistency in the Trajectory Testing
real(dp),dimension(6,Nbonds) :: INITIAL_BOND_DATA_test

!Timing Variables
real :: r1,r2
integer :: seed,c1,c2,cr
real :: system_clock_rate
integer :: grid_t0, grid_t1
real :: grid_wall_time
integer,dimension(3) :: now

!Incremental Integers
integer :: n, m, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		LIBRARY (MULTI-GRID) INITIALIZATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Get a random seed and print it in case there's a problem you need to replicate
call system_clock(seed)
print *, ""
call itime(now)
write(6,FMT=FMTnow) now
print *, "     RNG Seed: ", seed
seed = rand(seed)

!Print statement
print *, ""
call itime(now)
write(6,FMT=FMTnow) now
print *, "Creation of directory ", gridpath0
print *, ""

write(Nbond_text,FMT="(I0.6)") Nbonds
write(Natom_text,FMT="(I0.6)") Natoms
write(FMTinitial,FMT="(A19)") "("//Nbond_text//"(6(F14.10)))"
write(FMTtimeslice,FMT="(A19)") "("//Natom_text//"(12(F12.7)))"
write(FMT2,FMT="(A22)") "("//Natom_text//"(6(1x,F14.10)))"
write(FMT3,FMT="(A22)") "("//Natom_text//"(3(1x,F14.10)))"  

!We need to initialize how we format the text
!This is big and unruly

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

!We also need to initialize a few scaling factors to
!speed up computation

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

!If there is no force approximation going on then
!set the subcell searching to zero
!if ((reject_flag).and.(.not.&
!        (gather_interpolation_flag))) ssm1 = 0

!Initialize the subcell search variables
call setSubcellSearchMax()



!Initialize the clock
call system_clock(count_rate=cr)
system_clock_rate = 1.0/real(cr)

call getPrefixText(prefix_text)

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

call system("mkdir "//gridpath5)

allocate(filechannels(1+Ngrid_max))

!The grid number will uniquely identify one trajectory
!Open all these files under filechannels
filechannels(1) = 1000 !+ OMP_GET_THREAD_NUM()

!Our clever way to uniquely make a filechannel
do m = 1, Ngrid_max
    filechannels(1+m) = 1000 + 69 * m
end do

!If we are reading frames from a premade trajectory
!Open up the file to read the trajectories from
if (readtrajectory_flag) &
    open(trajectorieschannel,file=gridpath4//readtrajectoryfile)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		GRID INITIALIZATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do Ngrid = 1, Ngrid_max

        !Time the grid creation
        call system_clock(grid_t0)

        !This formats the name of the grid directory (001/, 002/, ...)
        write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid

        !This names the paths formally
        gridpath1 = gridpath0//Ngrid_text//"/"
        gridpath2 = gridpath1//"grid/"

        !Gridpath1 is the directory 001/
        call system("mkdir "//gridpath1)

        !Gridpath2 is the directory housing the files with coordinates and gradients
        call system("mkdir "//gridpath2)

        !A special folder is used to house interpolation and other error files
!       call system("mkdir "//gridpath0//interpolationfolder)

        call itime(now)
        write(6,FMT=FMTnow) now
        print *, "           Making grid ", Ngrid_text 
        print *, ""

        !Inside the directory, we monitor the progress of the grid's creation
        open(progresschannel,file=gridpath1//progressfile)
        write(progresschannel,*) ""
        write(progresschannel,*) ""
        write(progresschannel,*) "Let's Start!"
        write(progresschannel,*) ""
        close(progresschannel)

        !We start off with zero trajectories
        Ntraj = 0

        !We start off with zero files
        Nfile = 0

        !We start off with no overcrowded files
        headers = 0
        headers_old = 0
        max_headers_delta = 0

        !For the checkstate, we specify we want to add
        !frames to the current grid
        grid_addition = Ngrid

        !If we are doing "grid-checking" by occasionally
        !doing a test trajectory, then we can reduce some
        !random variables by making them all start with
        !the same set of initial conditions.
        !
        !We store these in INITIAL_BOND_DATA_test:

        call InitialSampling3()
        INITIAL_BOND_DATA_test = INITIAL_BOND_DATA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		TRAJECTORY INITIALIZATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! TEST !!!!!!!!!!!
if (heatmap_evolution_flag) then
call addMultipleTrajectories()
exit
end if
!!!!!!!!!!!!!!!!!!!!



!!! ANOTHER TEST !!!
if (.false.) then
Ngrid_total = 1
call errorCheck2(filechannels)
exit
end if
!!!!!!!!!!!!!!!!!!!!

        do n = 1, Ntraj_max

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		TRAJECTORY ADDITION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !Get some random initial conditions for the trajectory
                !For now, we are only handling systems with H-H bonds

                call InitialSampling3()

                open(progresschannel,file=gridpath1//progressfile,position="append")
                write(progresschannel,*) ""
                write(progresschannel,*) ""
                write(progresschannel,*) "Starting trajectory ", n
                write(progresschannel,*) "  Initial Conditions: "
                do m = 1, Nbonds
                        write(progresschannel,*) "          BOND ", m, ":"
                        write(progresschannel,*) "                  Bond Distance: ", INITIAL_BOND_DATA(1,m)
                        write(progresschannel,*) "                   Bond Angle 1: ", INITIAL_BOND_DATA(4,m)
                        write(progresschannel,*) "                   Bond Angle 2: ", INITIAL_BOND_DATA(5,m)
                end do
                close(progresschannel)

                write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
                write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
                write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
                write(Ntraj_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") n

                !How many grids we can check, meaning all the grids
                !up until this grid
                Ngrid_total = Ngrid

                !If we are doing trajectory checking then open
                !a filechannel for each grid we are checking
                do m = 1, Ngrid_total
                    write(variable_length_text,FMT=FMT5_variable)&
                            Ngrid_text_length
                    write(Ngrid_text,FMT="(I0."//&
                            trim(adjustl(variable_length_text))//")")&
                            m 
                    open(filechannels(1+m),file=gridpath5//&
                            Ngrid_text//"_"//Ntraj_text//".dat")
                end do

                !We time how much time each trajectory takes, wall-time and CPU time
                call CPU_time(r1)
                call system_clock(c1)

                !Not much flexibility given to the case where we are
                !just reading in a trajectory
                if (readtrajectory_flag) then
                    read(trajectorieschannel,FMT=*) traj_index
                    write(traj_index_text,FMT="(I0.6)") traj_index
                    open(inputchannel,file=gridpath4//&
                        readtrajectoryfolder//&
                        traj_index_text,form="unformatted")
                    call readTrajectory(&
                            inputchannel,filechannels,&
                            coords_initial,velocities_initial,&
                            coords_final,velocities_final)
                    close(inputchannel)

                !Otherwise, we have some options
                !One option is intense grid checking
                else if (testtrajDetailedRMSD_flag) then

                    call runTestTrajectory2(&
                            filechannels(1:1+Ngrid_total),&
                            coords_initial,velocities_initial,&
                            coords_final,velocities_final)
                else

                    !If no grid checking then we keep
                    !cell searching to a minimum;
                    !we only check to divyUp
                    subcellsearch_max1 = 0

                    call runTrajectoryRewind1(&
                                filechannels(1:1+Ngrid_total),&
                                coords_initial,velocities_initial,&
                                coords_final,velocities_final)

!                   if (force_Permutations) then
!                       call runTrajectory_permute_cap(&
!                               filechannels(1:1+Ngrid_total),&
!                               coords_initial,velocities_initial,&
!                               coords_final,velocities_final)
!                   else
!                       call runTrajectory_cap(&
!                               filechannels(1:1+Ngrid_total),&
!                               coords_initial,velocities_initial,&
!                               coords_final,velocities_final)
!                   end if
                end if

!               !Otherwise, we have some options
!               else if (.not.(force_NoLabels)) then
!                   call runTrajectory1(filechannels,&
!                           coords_initial,velocities_initial,&
!                           coords_final,velocities_final)
!               else if (force_Duplicates) then
!                   call runTrajectory2(filechannels,&
!                           coords_initial,velocities_initial,&
!                           coords_final,velocities_final)
!               else if (force_Permutations) then
!                   call runTrajectory3(filechannels,&
!                           coords_initial,velocities_initial,&
!                           coords_final,velocities_final)
!               else
!                   print *, "error!"
!               end if

                call CPU_time(r2)
                call system_clock(c2)
                trajectory_CPU_time = r2 - r1
                trajectory_wall_time = (c2 -c1) * system_clock_rate

                !Close the filechannels once more
                do m = 1, Ngrid_total
                    close(filechannels(1+m))
                end do

                Ntraj = Ntraj + 1

                open(filechannel1,file=gridpath1//initialfile,&
                                  position="append")
                write(filechannel1,FMTinitial) ((INITIAL_BOND_DATA(l,m),l=1,6),m=1,Nbonds)
                close(filechannel1)


                open(progresschannel,file=gridpath1//progressfile,position="append")
                write(progresschannel,*) "Finished trajectory ", n
                write(progresschannel,*) "        Now we have ", Nfile, " files"
                write(progresschannel,*) "                      Wall Time: ",&
                                                trajectory_wall_time
                write(progresschannel,*) "                       CPU Time: ",&
                                                trajectory_CPU_time
                close(progresschannel)


                !There is an informatics files for data on the grid while creating
                open(filechannel1,file=gridpath1//informaticsfile,position="append")
!                write(filechannel1,FMTinformatics) trajectory_CPU_time/real(steps),trajectory_wall_time/real(steps), &
!                                                   Ntraj,header1-header1_old,header2-header2_old,Nfile,Norder1*100.0/steps
                write(filechannel1,FMTinformatics) trajectory_CPU_time/real(steps),trajectory_wall_time/real(steps), &
                                                   Ntraj,headers(1)-headers_old(1),headers(2)-headers_old(2),&
                                                   Nfile,Norder_total(2)*100.0/steps
                close(filechannel1)


                !This is a temporary file to store information on how fast
                !Children level cells are filling up
!                open(filechannel1,file=gridpath1//"maxframesofsubcells.dat",position="append")
!                write(filechannel1,FMT=*) Ntraj, min(maxval(counter1),overcrowd1), min(maxval(counter2),overcrowd2)
!                close(filechannel1)


                !There is a timeslice file for snapshots of the trajectory at the beginning and end
                open(filechannel1,file=gridpath1//timeslicefile,position="append")
                write(filechannel1,FMTtimeslice) &
                                                 ((coords_initial(l,m),l=1,3),m=1,Natoms),&
                                                 ((velocities_initial(l,m),l=1,3),m=1,Natoms),&
                                                 ((coords_final(l,m),l=1,3),m=1,Natoms),&
                                                 ((velocities_final(l,m),l=1,3),m=1,Natoms)
                close(filechannel1)

                do m = 1, Norder_max+1
                        max_headers_delta(m) = max(&
                        max_headers_delta(m),headers(m)-headers_old(m))
                        headers_old(m) = headers(m)
                end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                GRID CREATION MONITORING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !This big if-statement is if we want to monitor our grid creation
                !The check only happens every Ngrid_check
                !Right now it is spaced so that we get (at most) ten graphs over the period of its creation
                if ((modulo(n,Ngrid_check) == 0)) then
                        INITIAL_BOND_DATA = INITIAL_BOND_DATA_test

                        testtrajDetailedRMSD_flag = .true.
                        force_Neighbors = .true.

                        subcellsearch_max1 = (/ 0, 0 /)
                        subcellsearch_max2 = (/ 5, 5 /)

                        call makeCheckTrajectoryGraphs()
                        call getRMSDDifferences1(gridpath1//&
                                  prefix_text//"_"//Ntraj_text//"_PercentRMSDThreshold")

                        force_Neighbors = .false.
                        testtrajDetailedRMSD_flag = .false.
                end if

                !If it is taking too long then we stop
                if (trajectory_CPU_time > trajectory_CPU_time_max) then
                        Ntraj_allowed = min(Ntraj_allowed,Ntraj)
                        exit
                end if

        end do


	!Grid Timing
	call system_clock(grid_t1)
	grid_wall_time = (grid_t1-grid_t0)*system_clock_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		GRID FINAL ANALYSIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	open(progresschannel,file=gridpath1//progressfile,position="append")
	write(progresschannel,*) ""
	write(progresschannel,*) ""
	write(progresschannel,*) "Finished all trajectories for grid "//Ngrid_text//"/"
	write(progresschannel,*) "        Now we have ", Nfile, " files"
	write(progresschannel,*) "        Altogether the grid took this many seconds: ", grid_wall_time
	write(progresschannel,*) ""
	write(progresschannel,*) ""
	close(progresschannel)

!        call makeGridCreationGraph(Nfile,max_header1_delta,grid_wall_time)
        call makeGridCreationGraph(Nfile,max_headers_delta(1),grid_wall_time)

        do n = 3, 7
                write(vals_interpolation_text,&
                        FMT="(F7.3,'_',F7.3)") &
                        n*1.0+0.7,n*1.0+1.1
                if (interpolation_check_visual) call &
                        getRMSDinterpolation(&
                        (/n*1.0d0+0.7d0,n*1.0d0+1.1d0/),&
                        (/0.2d0,0.2d0/),&
                        vals_interpolation_text//&
                        "InterpolationCheck")
        end do

        if (interpolation_check_visual) call &
                getRMSDinterpolation(&
                (/0.0d0,0.0d0/),&
                (/1.0d3,1.0d3/),&
                "InterpolationCheck")

!        call makeMaxFramesOfSubcellsGraph()

        max_absenergychange = 0.0
        min_absenergychange = 1.0e9
        max_relenergychange = 0.0
        min_relenergychange = 1.0e9
        max_rotenergychange = 0.0
        min_rotenergychange = 1.0e9
        
        !Finally, do a post-creation timeslice-to-SA conversions here
        !We use the SA often so we do this at the beginning
        call postProcess(Ngrid_text//"/")

        !Also, make a scattering angle plot
        call getScatteringAngles2(Ngrid_text//"/","InitialScatteringAngleDistribution_"//Ngrid_text)

        !Also, make an initial bond distribution plot
        call getInitialimages(Ngrid_text//"/","InitialBondDistribution_"//Ngrid_text)
end do

if (readtrajectory_flag) close(trajectorieschannel)

deallocate(filechannels)

print *, ""
call itime(now)
write(6,FMT=FMTnow) now
print *, "Successfully exited grid creation"
print *, ""

end program makeGridwithNewTrajectories







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               makeCheckTrajectoryGraphs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!       PURPOSE
!               This subroutine uses whatever initial conditions are set and starts two trajectories with the grid
!               checking procedure: one rejecting all approximations and one using the user-defined approximation
!
!               This then plots some basic information about the trajectory
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               reject_flag                     LOGICAL                         If true, then reject all approximations;
!                                                                               otherwise, use approximations
!               prefix_text                     CHAR(12)                        A string desribing the approximation method and
!                                                                               the threshold of acceptance
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath1//checkTrajectory//    PNG                             A plot describing the evolution of a handful of
!                   prefix_text//_//#traj                                       variables over the course of a trajectory
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
character(trajectory_text_length) :: Ntraj_text
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
real :: system_clock_rate
integer :: cr
integer :: trajectory_t0, trajectory_t1
real :: checktrajectory_wall_time
character(10) :: checktrajectory_wall_time_text
integer,dimension(3) :: now

integer :: k

!Initialize the clock
call system_clock(count_rate=cr)
system_clock_rate = 1.0/real(cr)

call InitialSampling3()

write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
write(Ntraj_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ntraj

!The grid number will uniquely identify one trajectory
!Open all these files under filechannels
allocate(filechannels(1+Ngrid_total))
write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
filechannels(1) = 1000 !+ OMP_GET_THREAD_NUM()

do k = 1, Ngrid_total
        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") k
        filechannels(1+k) = 1000 + 69 * k !+ OMP_GET_THREAD_NUM()
        open(filechannels(1+k),file=gridpath0//Ngrid_text//"/dump.dat")
end do

reject_flag = (.not.(reject_flag))
call getPrefixText(prefix_text)

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
open(filechannel1,file=gridpath0//checkstatefile)
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
print *, '   Making plot: "'//gridpath0//'checkTrajectory_'//prefix_text//'.png"'
print *, ""

!Finally, plot the data obtained from this trajectory
open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath1//'checkTrajectory_'//Ntraj_text//"_"//prefix_text//'.png"'
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
!write(angle1descriptor,FMT=FMT6_neg_real1) initial_bond_angle1
!write(angle2descriptor,FMT=FMT6_pos_real1) initial_bond_angle2
!write(bond1descriptor,FMT=FMT6_pos_real1) initial_bond_distance
 write(checktrajectory_wall_time_text,FMT="(F10.2)") checktrajectory_wall_time
!write(gnuplotchannel,*) 'set label 1 "H2 Orientation: '//angle1descriptor//', '//angle2descriptor//&
!                                                ' radians" at screen 0.6, 0.955'
!write(gnuplotchannel,*) 'set label 2 "H2 Bond Length: '//bond1descriptor//&
!                                                ' A" at screen 0.6, 0.94'
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
write(gnuplotchannel,*) 'plot "'//gridpath0//checkstatefile//&
                        '" u (($4)/steps_scaling):7 w lines'

!write(gnuplotchannel,*) 'unset label 1'
!write(gnuplotchannel,*) 'unset label 2'
write(gnuplotchannel,*) 'unset label 3'

write(gnuplotchannel,*) 'set ylabel "Var2 (A)"'
write(gnuplotchannel,*) 'set yrange [min_var2-delta_var2*.25 : max_var2+delta_var2*.25]'
write(gnuplotchannel,*) 'set ytics min_var2, delta_var2, max_var2'
write(gnuplotchannel,*) 'plot "'//gridpath0//checkstatefile//&
                        '" u (($4)/steps_scaling):8 w lines'

!write(gnuplotchannel,*) 'set ylabel "Total Energy (eV)"'
!write(gnuplotchannel,*) 'set yrange [0:0.02]'
!write(gnuplotchannel,*) 'set ytics 0.005'
!write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//&
!                        '" u (($4)/steps_scaling):($9+$10) w lines'

write(gnuplotchannel,*) 'set ylabel "Number of\nFrames Checked"'
write(gnuplotchannel,*) 'set yrange [-max_frames*.10:max_frames+max_frames*.10]'
write(gnuplotchannel,*) 'set ytics 0, floor(max_frames*.25), max_frames'
write(gnuplotchannel,*) 'plot "'//gridpath0//checkstatefile//&
                        '" u (($4)/steps_scaling):1 w points'

write(gnuplotchannel,*) 'set ylabel "Order of Cells\nChecked"'
write(gnuplotchannel,*) 'set yrange [-0.25:2.25]'
write(gnuplotchannel,*) 'set ytics 1'
write(gnuplotchannel,*) 'plot "'//gridpath0//checkstatefile//&
                        '" u (($4)/steps_scaling):2 w points'

write(gnuplotchannel,*) 'set ylabel "Number of Cells\nChecked"'
write(gnuplotchannel,*) 'set yrange [0:max_neighbor_check+max_neighbor_check*.10]'
write(gnuplotchannel,*) 'set ytics 1, floor(max_neighbor_check*.25), max_neighbor_check'
write(gnuplotchannel,*) 'plot "'//gridpath0//checkstatefile//&
                        '" u (($4)/steps_scaling):3 w points'

write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set format x'
write(gnuplotchannel,*) 'set xlabel "Timesteps (Thousands)"'

write(gnuplotchannel,*) 'set yrange [min_rmsd*0.9:',default_rmsd,']'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,*) 'set ylabel "Timestep\nRMSD (A)"'
!write(gnuplotchannel,*) 'set ytics (".1" .1, ".05" .05, ".01" .01, ".001" .001, ".0001" .0001)'
write(gnuplotchannel,*) 'set ytics ("5e-1" .5, "5e-2" .05, "5e-3" .005,'//&
                                   '"5e-4" .0005, "5e-5" .00005)'
write(gnuplotchannel,*) 'plot "'//gridpath0//checkstatefile//&
                        '" u (($4)/steps_scaling):($6==$5?$5:1/0) w points lc rgb "blue",\'
write(gnuplotchannel,*) '     "'//gridpath0//checkstatefile//&
                        '" u (($4)/steps_scaling):($6>$5?$5:1/0) w points lc rgb "red",\'
write(gnuplotchannel,*) '     "'//gridpath0//checkstatefile//&
                        '" u (($4)/steps_scaling):($6<$5?$6:1/0) w points lc rgb "green"'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)

reject_flag = (.not.(reject_flag))
call getPrefixText(prefix_text)

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
open(filechannel1,file=gridpath0//checkstatefile)
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
print *, '   Making plot: "'//gridpath0//'checkTrajectory_'//prefix_text//'.png"'
print *, ""

!Finally, plot the data obtained from this trajectory
open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath1//'checkTrajectory_'//Ntraj_text//"_"//prefix_text//'.png"'
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
!write(angle1descriptor,FMT=FMT6_neg_real1) initial_bond_angle1
!write(angle2descriptor,FMT=FMT6_pos_real1) initial_bond_angle2
!write(bond1descriptor,FMT=FMT6_pos_real1) initial_bond_distance
 write(checktrajectory_wall_time_text,FMT="(F10.2)") checktrajectory_wall_time
!write(gnuplotchannel,*) 'set label 1 "H2 Orientation: '//angle1descriptor//', '//angle2descriptor//&
!                                                ' radians" at screen 0.6, 0.955'
!write(gnuplotchannel,*) 'set label 2 "H2 Bond Length: '//bond1descriptor//&
!                                                ' A" at screen 0.6, 0.94'
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
write(gnuplotchannel,*) 'plot "'//gridpath0//checkstatefile//&
                        '" u (($4)/steps_scaling):7 w lines'

!write(gnuplotchannel,*) 'unset label 1'
!write(gnuplotchannel,*) 'unset label 2'
write(gnuplotchannel,*) 'unset label 3'

write(gnuplotchannel,*) 'set ylabel "Var2 (A)"'
write(gnuplotchannel,*) 'set yrange [min_var2-delta_var2*.25 : max_var2+delta_var2*.25]'
write(gnuplotchannel,*) 'set ytics min_var2, delta_var2, max_var2'
write(gnuplotchannel,*) 'plot "'//gridpath0//checkstatefile//&
                        '" u (($4)/steps_scaling):8 w lines'

!write(gnuplotchannel,*) 'set ylabel "Total Energy (eV)"'
!write(gnuplotchannel,*) 'set yrange [0:0.02]'
!write(gnuplotchannel,*) 'set ytics 0.005'
!write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//&
!                        '" u (($4)/steps_scaling):($9+$10) w lines'

write(gnuplotchannel,*) 'set ylabel "Number of\nFrames Checked"'
write(gnuplotchannel,*) 'set yrange [-max_frames*.10:max_frames+max_frames*.10]'
write(gnuplotchannel,*) 'set ytics 0, floor(max_frames*.25), max_frames'
write(gnuplotchannel,*) 'plot "'//gridpath0//checkstatefile//&
                        '" u (($4)/steps_scaling):1 w points'

write(gnuplotchannel,*) 'set ylabel "Order of Cells\nChecked"'
write(gnuplotchannel,*) 'set yrange [-0.25:2.25]'
write(gnuplotchannel,*) 'set ytics 1'
write(gnuplotchannel,*) 'plot "'//gridpath0//checkstatefile//&
                        '" u (($4)/steps_scaling):2 w points'

write(gnuplotchannel,*) 'set ylabel "Number of Cells\nChecked"'
write(gnuplotchannel,*) 'set yrange [0:max_neighbor_check+max_neighbor_check*.10]'
write(gnuplotchannel,*) 'set ytics 1, floor(max_neighbor_check*.25), max_neighbor_check'
write(gnuplotchannel,*) 'plot "'//gridpath0//checkstatefile//&
                        '" u (($4)/steps_scaling):3 w points'

write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set format x'
write(gnuplotchannel,*) 'set xlabel "Timesteps (Thousands)"'

write(gnuplotchannel,*) 'set yrange [min_rmsd*0.9:',default_rmsd,']'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,*) 'set ylabel "Timestep\nRMSD (A)"'
!write(gnuplotchannel,*) 'set ytics (".1" .1, ".05" .05, ".01" .01, ".001" .001, ".0001" .0001)'
write(gnuplotchannel,*) 'set ytics ("5e-1" .5, "5e-2" .05, "5e-3" .005,'//&
                                   '"5e-4" .0005, "5e-5" .00005)'
write(gnuplotchannel,*) 'plot "'//gridpath0//checkstatefile//&
                        '" u (($4)/steps_scaling):($6==$5?$5:1/0) w points lc rgb "blue",\'
write(gnuplotchannel,*) '     "'//gridpath0//checkstatefile//&
                        '" u (($4)/steps_scaling):($6>$5?$5:1/0) w points lc rgb "red",\'
write(gnuplotchannel,*) '     "'//gridpath0//checkstatefile//&
                        '" u (($4)/steps_scaling):($6<$5?$6:1/0) w points lc rgb "green"'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)

deallocate(filechannels)
return

end subroutine makeCheckTrajectoryGraphs




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               makeGridCreationGraph
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!       PURPOSE
!               This subroutine plots the number of frames in the subcell that has the most frames over the
!               course of the grid's creation. This information is stored in a specific file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               Nfile                           INTEGER                         How many files are in the grid
!               max_header1_delta               INTEGER                         The maximum numer of times divyUp was
!                                                                               called in one trajectory
!               grid_wall_time                  REAL                            The amount of wall time that elapsed
!                                                                               over the creation of the grid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath1//                     DAT                             A file that holds information on various
!                   informaticsfile                                             file and counter variables of the grid
!                                                                               over the course of its creation
!               gridpath1//                     PNG                             A plot describing the evolution of some
!                   GridCreationGraph                                           file and counter information of the grid
!                                                                               over the course of its creation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine makeGridCreationGraph(Nfile_in, max_header1_delta_in, grid_wall_time_in)
use PARAMETERS
use FUNCTIONS
use VARIABLES
use ANALYSIS
use PHYSICS
implicit none

integer,intent(in) :: Nfile_in, max_header1_delta_in
real,intent(in) :: grid_wall_time_in
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(10) :: grid_wall_time_text
integer,dimension(3) :: now

print *, ""
call itime(now)
write(6,FMT=FMTnow) now
print *, '   Making plot: "'//gridpath1//'GridCreationGraph.png"'
print *, ""

write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid

!Now that we are done with a grid, we can see how much time the grid creation took
!There is also some other interesting data
open(gnuplotchannel,file=gridpath1//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath0//expfolder//&
                        'GridCreationGraph.png_'//Ngrid_text//'"'
write(gnuplotchannel,*) 'set style line 1 lc rgb "red" pt 5'
write(gnuplotchannel,*) 'set style line 2 lc rgb "green" pt 7'
write(gnuplotchannel,*) 'set style line 3 lc rgb "blue" pt 13'
write(gnuplotchannel,*) 'set style line 4 lc rgb "orange" pt 9'
write(gnuplotchannel,*) 'set style line 5 lc rgb "yellow" pt 11'
write(gnuplotchannel,*) 'set style line 6 lc rgb "pink" pt 20'
write(gnuplotchannel,*) 'unset xtics'
write(gnuplotchannel,*) 'tic_spacing = 0.20'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
write(gnuplotchannel,*) 'set multiplot layout 5,1 margins 0.15,0.95,.1,.95 spacing 0,0 title '//&
                                '"Logistics of Trajectory Addition for Grid '//Ngrid_text//'/" font ",18" offset 0,3'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'unset xlabel'

write(gnuplotchannel,*) 'set ylabel "Number of Files\n(Thousands)"'
write(gnuplotchannel,*) 'delta_Nfile = (',Nfile_in,' / 5.0)/1000.0'
write(gnuplotchannel,*) 'set yrange [-delta_Nfile*tic_spacing:delta_Nfile*(5.0+tic_spacing)]'
write(gnuplotchannel,*) 'set ytics 0, floor(delta_Nfile)'
write(gnuplotchannel,*) 'plot "'//gridpath1//informaticsfile//'" u 3:(($6)/1000.0) w lines'

write(gnuplotchannel,*) 'set ylabel "Number of Calls\nto DivyUp"'
write(gnuplotchannel,*) 'delta_header1 = (',max_header1_delta_in,' / 5.0)'
write(gnuplotchannel,*) 'if (delta_header1==0.0) { set yrange [-0.2:1.2];'
write(gnuplotchannel,*) 'set ytics 1'
write(gnuplotchannel,*) '} else {'
write(gnuplotchannel,*) 'set yrange [-delta_header1*tic_spacing:delta_header1*(5.0+tic_spacing)];'
write(gnuplotchannel,*) 'set ytics 0, floor(delta_header1)'
write(gnuplotchannel,*) '}'
write(gnuplotchannel,*) 'plot "'//gridpath1//informaticsfile//'" u 3:4 w lines, '//&
                           '"'//gridpath1//informaticsfile//'" u 3:5 w lines'

write(gnuplotchannel,*) 'set ylabel "Percentage of Frame\nSearches in Order 1"'
write(gnuplotchannel,*) 'delta_percentage = 20'
write(gnuplotchannel,*) 'set yrange [-delta_percentage*tic_spacing:100+delta_percentage*tic_spacing]'
write(gnuplotchannel,*) 'set ytics 0, delta_percentage, 100'
write(gnuplotchannel,*) 'plot "'//gridpath1//informaticsfile//'" u 3:7 w lines'

write(grid_wall_time_text,FMT="(F10.2)") grid_wall_time_in
write(gnuplotchannel,*) 'set label 1 "Total Wall Time (including grid checking): '//&
                                trim(adjustl(grid_wall_time_text))//' s" at graph 0.025,0.9'

write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set ytics autofreq'
write(gnuplotchannel,*) 'set format y "%.1e"'
write(gnuplotchannel,*) 'set ylabel "Wall Time\nPer Frame (ms)"'
write(gnuplotchannel,*) 'plot "'//gridpath1//informaticsfile//'" u 3:(($2)*1000.0) w lines'

write(gnuplotchannel,*) 'unset label 1'

write(gnuplotchannel,*) 'Ntraj_max = ', Ntraj_max
write(gnuplotchannel,*) 'set xrange [1:Ntraj_max]'
write(gnuplotchannel,*) 'set xtics nomirror 0,floor(Ntraj_max*.10), Ntraj_max'
write(gnuplotchannel,*) 'set xlabel "Trajectories"'

write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set ytics autofreq'
write(gnuplotchannel,*) 'set format y "%.1e"'
write(gnuplotchannel,*) 'set ylabel "CPU Time\nPer Frame (ms)"'
write(gnuplotchannel,*) 'plot "'//gridpath1//informaticsfile//'" u 3:(($1)*1000.0) w lines'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath1//gnuplotfile)

end subroutine makeGridCreationGraph


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               makeMaxFramesOfSubcellsGraph
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!       PURPOSE
!               This subroutine plots the number of frames in the subcell that has the most frames over the
!               course of the grid's creation. This information is stored in a specific file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath1//                     DAT                             A file that holds information on how many
!                   maxframesofsubcells                                         frames are in the most crowded cell of
!                                                                               each order
!               gridpath1//                     PNG                             A plot describing the number of frames in
!                   MaxFramesOfSubcells                                         the subcell that has the most frames over
!                                                                               the course of the grid creation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine makeMaxFramesOfSubcellsGraph()
use PARAMETERS
use FUNCTIONS
use VARIABLES
use ANALYSIS
use PHYSICS
implicit none

integer,dimension(3) :: now

print *, ""
call itime(now)
write(6,FMT=FMTnow) now
print *, '   Making plot: "'//gridpath1//'MaxFramesOfSubcells.png"'
print *, ""

!We can also glean some information on the concentration of frames in subcells by checking
!how many frames the most crowded subcell had over the course of the grid's creation
open(gnuplotchannel,file=gridpath1//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath1//'MaxFramesOfSubcells.png"'
write(gnuplotchannel,*) 'set style line 1 lc rgb "red" pt 5'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set xlabel "Trajectories"'
write(gnuplotchannel,*) 'set ylabel "Max Number of Frames in Subcells"'
write(gnuplotchannel,*) 'plot "'//gridpath1//'maxframesofsubcells.dat" u 1:2 t "Order 1" w lines,\'
write(gnuplotchannel,*) '     "'//gridpath1//'maxframesofsubcells.dat" u 1:3 t "Order 2" w lines'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath1//gnuplotfile)

end subroutine makeMaxFramesOfSubcellsGraph

