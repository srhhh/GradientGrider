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
character(12) :: prefix_text

!Trajectory Variables
real(dp) :: trajectory_CPU_time,trajectory_wall_time
real(dp) :: scattering_angle
real(dp),dimension(3) :: TRVenergies1,TRVenergies2,dTRVenergies
real(dp) :: totalEnergy
integer :: header1_old,header2_old,header3_old
integer :: max_header1_delta, max_header2_delta, max_header3_delta

!Trajectory Output
real(dp),dimension(3,Natoms) :: coords_initial, velocities_initial
real(dp),dimension(3,Natoms) :: coords_final, velocities_final

!Timing Variables
real :: r1,r2
integer :: seed,c1,c2,cr
real :: system_clock_rate
integer :: grid_t0, grid_t1
real :: grid_wall_time
character(10) :: grid_wall_time_text
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

!Initialize the clock
call system_clock(count_rate=cr)
system_clock_rate = 1.0/real(cr)

call getPrefixText(prefix_text)

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
        header1 = 1
	header1_old = 1
        max_header1_delta = 0
        header2 = 1
	header2_old = 1
        max_header2_delta = 0
        header3 = 1
	header3_old = 1
        max_header3_delta = 0

	!We start off with zero frames in all files
        do m = 1, counter0_max
                counter0(m) = 0
        end do
        do m = 1, counter1_max
                counter1(m) = 0
        end do
        do m = 1, counter2_max
                counter2(m) = 0
        end do
        do m = 1, counter3_max
                counter3(m) = 0
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		TRAJECTORY INITIALIZATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !$OMP PARALLEL
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,m,l,INITIAL_BOND_DATA,steps,&
        !                                       coords_initial,velocities_initial,coords_final,velocities_final,&
        !                                       trajectory_CPU_time,trajectory_wall_time,r1,r2,c1,c2)
        !$OMP DO





!!! TEST !!!
if (heatmap_evolution_flag) then
call addMultipleTrajectories()
exit
end if
!!!!!!!!!!!!

        do n = 1, Ntraj_max

                !Get some random initial conditions for the trajectory
                !For now, we are only handling systems with H-H bonds

                call InitialSampling3()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                GRID CREATION MONITORING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !This big if-statement is if we want to monitor our grid creation
                !The check only happens every Ngrid_check
                !Right now it is spaced so that we get (at most) ten graphs over the period of its creation
                if ((modulo(n,Ngrid_check) == 0)) then
                        call makeCheckTrajectoryGraphs()
                end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		TRAJECTORY ADDITION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!$OMP CRITICAL

                open(progresschannel,file=gridpath1//progressfile,position="append")
                write(progresschannel,*) ""
                write(progresschannel,*) ""
                write(progresschannel,*) "Starting trajectory ", n
                write(progresschannel,*) "  Initial Conditions: "
		do m = 1, Nbonds
	                write(progresschannel,*) "          BOND ", m, ":"
	                write(progresschannel,*) "  		Bond Distance: ", INITIAL_BOND_DATA(1,m)
	                write(progresschannel,*) "  		 Bond Angle 1: ", INITIAL_BOND_DATA(4,m)
	                write(progresschannel,*) "  		 Bond Angle 2: ", INITIAL_BOND_DATA(5,m)
		end do
                close(progresschannel)

		!$OMP END CRITICAL

		!We time how much time each trajectory takes, wall-time and CPU time
                call CPU_time(r1)
                call system_clock(c1)
                call addTrajectory(coords_initial,velocities_initial,coords_final,velocities_final)
                call CPU_time(r2)
                call system_clock(c2)
                trajectory_CPU_time = r2 - r1
                trajectory_wall_time = (c2 -c1) * system_clock_rate

		!If there have been a large number of subdivisions (so many that our array will go
		!out of bounds) then we stop; we also stop if it is taking too long
                if ((header1 == header1_max).or.&
                   (trajectory_CPU_time > trajectory_CPU_time_max)) then
			Ntraj_allowed = min(Ntraj_allowed,Ntraj)
			exit
		end if

                Ntraj = Ntraj + 1

		!$OMP CRITICAL

                open(filechannel1,file=gridpath0//Ngrid_text//"/Initial"//initialfile,&
                                  position="append")
                write(filechannel1,FMTinitial) ((INITIAL_BOND_DATA(l,m),l=1,6),m=1,Nbonds)
                close(filechannel1)

		!$OMP END CRITICAL

		!$OMP CRITICAL

                open(progresschannel,file=gridpath1//progressfile,position="append")
                write(progresschannel,*) "Finished trajectory ", n
                write(progresschannel,*) "        Now we have ", Nfile, " files"
                write(progresschannel,*) "                      Wall Time: ",&
                                                trajectory_wall_time
                write(progresschannel,*) "                       CPU Time: ",&
                                                trajectory_CPU_time
                close(progresschannel)

		!$OMP END CRITICAL

		!$OMP CRITICAL

		!There is an informatics files for data on the grid while creating
		open(filechannel1,file=gridpath1//"Initial"//informaticsfile,position="append")
		write(filechannel1,FMTinformatics) trajectory_CPU_time/real(steps),trajectory_wall_time/real(steps), &
		                                   Ntraj,header1-header1_old,header2-header2_old,Nfile,Norder1*100.0/steps
		close(filechannel1)

		!$OMP END CRITICAL

		!$OMP CRITICAL

		!This is a temporary file to store information on how fast
		!Children level cells are filling up
		open(filechannel1,file=gridpath1//"maxframesofsubcells.dat",position="append")
		write(filechannel1,FMT=*) Ntraj, min(maxval(counter1),overcrowd1), min(maxval(counter2),overcrowd2)
		close(filechannel1)

		!$OMP END CRITICAL

		!$OMP CRITICAL

		!There is a timeslice file for snapshots of the trajectory at the beginning and end
		open(filechannel1,file=gridpath1//"Initial"//timeslicefile,position="append")
		write(filechannel1,FMTtimeslice) &
                                                 ((coords_initial(l,m),l=1,3),m=1,Natoms),&
                                                 ((velocities_initial(l,m),l=1,3),m=1,Natoms),&
         				         ((coords_final(l,m),l=1,3),m=1,Natoms),&
                                                 ((velocities_final(l,m),l=1,3),m=1,Natoms)
                close(filechannel1)

		!$OMP END CRITICAL

                max_header1_delta = max(max_header1_delta,header1-header1_old)
                max_header2_delta = max(max_header2_delta,header2-header2_old)
                max_header3_delta = max(max_header3_delta,header3-header3_old)

		header1_old = header1
		header2_old = header2
		header3_old = header3
        end do

	!$OMP END PARALLEL
	!$OMP END DO NOWAIT


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

        call makeGridCreationGraph(Nfile,max_header1_delta,grid_wall_time)

        call makeMaxFramesOfSubcellsGraph()

        !Finally, we save all of the counters to their respective counter files in the folder
        open(filechannel1,file=trim(gridpath1)//counter0file)
        do n = 1, counter0_max
                write(filechannel1,FMT=FMT8_counter) counter0(n)
        end do
        close(filechannel1)
        
        open(filechannel1,file=trim(gridpath1)//counter1file)
        do n = 1, counter1_max
                write(filechannel1,FMT=FMT8_counter) counter1(n)
        end do
        close(filechannel1)
        
        open(filechannel1,file=trim(gridpath1)//counter2file)
        do n = 1, counter2_max
                write(filechannel1,FMT=FMT8_counter) counter2(n)
        end do
        close(filechannel1)
        
        open(filechannel1,file=trim(gridpath1)//counter3file)
        do n = 1, counter3_max
                write(filechannel1,FMT=FMT8_counter) counter3(n)
        end do
        close(filechannel1)
        
        max_absenergychange = 0.0
        min_absenergychange = 1.0e9
        max_relenergychange = 0.0
        min_relenergychange = 1.0e9
        max_rotenergychange = 0.0
        min_rotenergychange = 1.0e9
        
        !Finally, do a post-creation timeslice-to-SA conversions here
        !We use the SA often so we do this at the beginning
        call postProcess(Ngrid_text//"/Initial")

        !Also, make a scattering angle plot
        call getScatteringAngles2(Ngrid_text//"/Initial","InitialScatteringAngleDistribution_"//Ngrid_text)

        !Also, make an initial bond distribution plot
        call getInitialimages(Ngrid_text//"/Initial","InitialBondDistribution_"//Ngrid_text)
end do

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
character(12) :: prefix_text

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

!Initialize the clock
call system_clock(count_rate=cr)
system_clock_rate = 1.0/real(cr)

!Check first how the trajectory behaves with the OPPOSITE
!rejection method
reject_flag = (.not.(reject_flag))
call getPrefixText(prefix_text)

!Next, run the trajectory
!Remark: output of checkTrajectory is in the checkstatefile
!
call system_clock(trajectory_t0)
call checkTrajectory(coords_initial,velocities_initial,coords_final,velocities_final)
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
open(filechannel1,file=gridpath1//checkstatefile)
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

!To uniquely define this trajectory from other trajectories we check
!the trajectory number is concatenated to the image name
!
write(checkstateTrajectory,FMT=FMT6_pos_int) Ntraj

print *, ""
call itime(now)
write(6,FMT=FMTnow) now
print *, "   Making plot: "//gridpath1//'checkTrajectory_'//prefix_text//'_'//checkstateTrajectory//'.png"'
print *, ""

!Finally, plot the data obtained from this trajectory
open(gnuplotchannel,file=gridpath1//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath1//'checkTrajectory_'//prefix_text//&
                         '_'//checkstateTrajectory//'.png"'
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
write(gnuplotchannel,*) 'max_neighbor_check = ', trajectory_max_neighbor_check
write(gnuplotchannel,*) 'min_rmsd = ', trajectory_min_rmsd

write(gnuplotchannel,*) 'steps_scaling = 1000'
write(gnuplotchannel,*) 'set xrange [0:max_steps/steps_scaling]'
write(gnuplotchannel,*) 'set xtics 0, 1, max_steps/steps_scaling'
write(gnuplotchannel,*) 'set format x ""'
write(gnuplotchannel,*) 'unset xlabel'

write(gnuplotchannel,*) 'set ylabel "Var1 (A)"'
write(gnuplotchannel,*) 'set yrange [min_var1-delta_var1*.25 : max_var1+delta_var1*.25]'
write(gnuplotchannel,*) 'set ytics min_var1, delta_var1, max_var2'
write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//&
                        '" u (($4)/steps_scaling):7 w lines'

!write(gnuplotchannel,*) 'unset label 1'
!write(gnuplotchannel,*) 'unset label 2'
write(gnuplotchannel,*) 'unset label 3'

write(gnuplotchannel,*) 'set ylabel "Var2 (A)"'
write(gnuplotchannel,*) 'set yrange [min_var2-delta_var2*.25 : max_var2+delta_var2*.25]'
write(gnuplotchannel,*) 'set ytics min_var2, delta_var2, max_var2'
write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//&
                        '" u (($4)/steps_scaling):8 w lines'

!write(gnuplotchannel,*) 'set ylabel "Total Energy (eV)"'
!write(gnuplotchannel,*) 'set yrange [0:0.02]'
!write(gnuplotchannel,*) 'set ytics 0.005'
!write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//&
!                        '" u (($4)/steps_scaling):($9+$10) w lines'

write(gnuplotchannel,*) 'set ylabel "Number of\nFrames Checked"'
write(gnuplotchannel,*) 'set yrange [-max_frames*.10:max_frames+max_frames*.10]'
write(gnuplotchannel,*) 'set ytics 0, floor(max_frames*.25), max_frames'
write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//&
                        '" u (($4)/steps_scaling):1 w lines'

write(gnuplotchannel,*) 'set ylabel "Order of Cells\nChecked"'
write(gnuplotchannel,*) 'set yrange [-1.25:1.25]'
write(gnuplotchannel,*) 'set ytics 1'
write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//&
                        '" u (($4)/steps_scaling):2 w lines'

write(gnuplotchannel,*) 'set ylabel "Number of Cells\nChecked"'
write(gnuplotchannel,*) 'set yrange [0:max_neighbor_check+max_neighbor_check*.10]'
write(gnuplotchannel,*) 'set ytics 1, floor(max_neighbor_check*.25), max_neighbor_check'
write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//&
                        '" u (($4)/steps_scaling):3 w lines'

write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set format x'
write(gnuplotchannel,*) 'set xlabel "Timesteps (Thousands)"'

write(gnuplotchannel,*) 'set yrange [min_rmsd*0.9:',default_rmsd,']'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,*) 'set ylabel "Timestep\nRMSD (A)"'
!write(gnuplotchannel,*) 'set ytics (".1" .1, ".05" .05, ".01" .01, ".001" .001, ".0001" .0001)'
write(gnuplotchannel,*) 'set ytics ("5e-1" .5, "5e-2" .05, "5e-3" .005,'//&
                                   '"5e-4" .0005, "5e-5" .00005)'
write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//&
                        '" u (($4)/steps_scaling):5 w lines'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath1//gnuplotfile)



!Second, check the trajectory with the rejection method
!from before (as specified by the user)
reject_flag = (.not.(reject_flag))
call getPrefixText(prefix_text)

!Run the trajectory as usual
!Remark: output of checkTrajectory is in the checkstatefile
!
call system_clock(trajectory_t0)
call checkTrajectory(coords_initial,velocities_initial,coords_final,velocities_final)
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
open(filechannel1,file=gridpath1//checkstatefile)
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

!To uniquely define this trajectory from other trajectories we check
!the trajectory number is concatenated to the image name
!
write(checkstateTrajectory,FMT=FMT6_pos_int) Ntraj

print *, ""
call itime(now)
write(6,FMT=FMTnow) now
print *, "   Making plot: "//gridpath1//'checkTrajectory_'//prefix_text//'_'//checkstateTrajectory//'.png"'
print *, ""

!Finally, plot the data obtained from this trajectory
open(gnuplotchannel,file=gridpath1//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath1//'checkTrajectory_'//prefix_text//&
                         '_'//checkstateTrajectory//'.png"'
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
write(gnuplotchannel,*) 'max_neighbor_check = ', trajectory_max_neighbor_check
write(gnuplotchannel,*) 'min_rmsd = ', trajectory_min_rmsd

write(gnuplotchannel,*) 'steps_scaling = 1000'
write(gnuplotchannel,*) 'set xrange [0:max_steps/steps_scaling]'
write(gnuplotchannel,*) 'set xtics 0, 1, max_steps/steps_scaling'
write(gnuplotchannel,*) 'set format x ""'
write(gnuplotchannel,*) 'unset xlabel'

write(gnuplotchannel,*) 'set ylabel "Var1 (A)"'
write(gnuplotchannel,*) 'set yrange [min_var1-delta_var1*.25 : max_var1+delta_var1*.25]'
write(gnuplotchannel,*) 'set ytics min_var1, delta_var1, max_var2'
write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//&
                        '" u (($4)/steps_scaling):7 w lines'

!write(gnuplotchannel,*) 'unset label 1'
!write(gnuplotchannel,*) 'unset label 2'
write(gnuplotchannel,*) 'unset label 3'

write(gnuplotchannel,*) 'set ylabel "Var2 (A)"'
write(gnuplotchannel,*) 'set yrange [min_var2-delta_var2*.25 : max_var2+delta_var2*.25]'
write(gnuplotchannel,*) 'set ytics min_var2, delta_var2, max_var2'
write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//&
                        '" u (($4)/steps_scaling):8 w lines'

!write(gnuplotchannel,*) 'set ylabel "Total Energy (eV)"'
!write(gnuplotchannel,*) 'set yrange [0:0.02]'
!write(gnuplotchannel,*) 'set ytics 0.005'
!write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//&
!                        '" u (($4)/steps_scaling):($9+$10) w lines'

write(gnuplotchannel,*) 'set ylabel "Number of\nFrames Checked"'
write(gnuplotchannel,*) 'set yrange [-max_frames*.10:max_frames+max_frames*.10]'
write(gnuplotchannel,*) 'set ytics 0, floor(max_frames*.25), max_frames'
write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//&
                        '" u (($4)/steps_scaling):1 w lines'

write(gnuplotchannel,*) 'set ylabel "Order of Cells\nChecked"'
write(gnuplotchannel,*) 'set yrange [-1.25:1.25]'
write(gnuplotchannel,*) 'set ytics 1'
write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//&
                        '" u (($4)/steps_scaling):2 w lines'

write(gnuplotchannel,*) 'set ylabel "Number of Cells\nChecked"'
write(gnuplotchannel,*) 'set yrange [0:max_neighbor_check+max_neighbor_check*.10]'
write(gnuplotchannel,*) 'set ytics 1, floor(max_neighbor_check*.25), max_neighbor_check'
write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//&
                        '" u (($4)/steps_scaling):3 w lines'

write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set format x'
write(gnuplotchannel,*) 'set xlabel "Timesteps (Thousands)"'

write(gnuplotchannel,*) 'set yrange [min_rmsd*0.9:',default_rmsd,']'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,*) 'set ylabel "Timestep\nRMSD (A)"'
!write(gnuplotchannel,*) 'set ytics (".1" .1, ".05" .05, ".01" .01, ".001" .001, ".0001" .0001)'
write(gnuplotchannel,*) 'set ytics ("5e-1" .5, "5e-2" .05, "5e-3" .005,'//&
                                   '"5e-4" .0005, "5e-5" .00005)'
write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//&
                        '" u (($4)/steps_scaling):5 w lines'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath1//gnuplotfile)

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

subroutine makeGridCreationGraph(Nfile, max_header1_delta, grid_wall_time)
use PARAMETERS
use FUNCTIONS
use VARIABLES
use ANALYSIS
use PHYSICS
implicit none

integer,intent(in) :: Nfile, max_header1_delta
real,intent(in) :: grid_wall_time
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
integer,dimension(3) :: now

print *, ""
call itime(now)
write(6,FMT=FMTnow) now
print *, "   Making plot: "//gridpath1//'GridCreationGraph.png"'
print *, ""

write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid

!Now that we are done with a grid, we can see how much time the grid creation took
!There is also some other interesting data
open(gnuplotchannel,file=gridpath1//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath1//'GridCreationGraph.png"'
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
write(gnuplotchannel,*) 'delta_Nfile = (',Nfile,' / 5.0)/1000.0'
write(gnuplotchannel,*) 'set yrange [-delta_Nfile*tic_spacing:delta_Nfile*(5.0+tic_spacing)]'
write(gnuplotchannel,*) 'set ytics 0, floor(delta_Nfile)'
write(gnuplotchannel,*) 'plot "'//gridpath1//"Initial"//informaticsfile//'" u 3:(($6)/1000.0) w lines'

write(gnuplotchannel,*) 'set ylabel "Number of Calls\nto DivyUp"'
write(gnuplotchannel,*) 'delta_header1 = (',max_header1_delta,' / 5.0)'
write(gnuplotchannel,*) 'set yrange [-delta_header1*tic_spacing:delta_header1*(5.0+tic_spacing)]'
write(gnuplotchannel,*) 'set ytics 0, floor(delta_header1)'
write(gnuplotchannel,*) 'plot "'//gridpath1//"Initial"//informaticsfile//'" u 3:4 w lines, '//&
                           '"'//gridpath1//"Initial"//informaticsfile//'" u 3:5 w lines'

write(gnuplotchannel,*) 'set ylabel "Percentage of Frames\nAdded to Order 1"'
write(gnuplotchannel,*) 'delta_percentage = 20'
write(gnuplotchannel,*) 'set yrange [-delta_percentage*tic_spacing:100+delta_percentage*tic_spacing]'
write(gnuplotchannel,*) 'set ytics 0, delta_percentage, 100'
write(gnuplotchannel,*) 'plot "'//gridpath1//"Initial"//informaticsfile//'" u 3:7 w lines'

write(grid_wall_time_text,FMT="(F10.2)") grid_wall_time
write(gnuplotchannel,*) 'set label 1 "Total Wall Time (including grid checking): '//&
                                trim(adjustl(grid_wall_time_text))//' s" at graph 0.025,0.9'

write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set ytics autofreq'
write(gnuplotchannel,*) 'set format y "%.1e"'
write(gnuplotchannel,*) 'set ylabel "Wall Time\nPer Frame (ms)"'
write(gnuplotchannel,*) 'plot "'//gridpath1//"Initial"//informaticsfile//'" u 3:(($2)*1000.0) w lines'

write(gnuplotchannel,*) 'unset label 1'

write(gnuplotchannel,*) 'Ntraj_max = ', Ntraj_max
write(gnuplotchannel,*) 'set xrange [1:Ntraj_max]'
write(gnuplotchannel,*) 'set xtics nomirror 0,floor(Ntraj_max*.10), Ntraj_max'
write(gnuplotchannel,*) 'set xlabel "Trajectories"'

write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set ytics autofreq'
write(gnuplotchannel,*) 'set format y "%.1e"'
write(gnuplotchannel,*) 'set ylabel "CPU Time\nPer Frame (ms)"'
write(gnuplotchannel,*) 'plot "'//gridpath1//"Initial"//informaticsfile//'" u 3:(($1)*1000.0) w lines'
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
print *, "   Making plot: "//gridpath1//'MaxFramesOfSubcells.png"'
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

