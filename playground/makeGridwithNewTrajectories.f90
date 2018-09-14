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
character(6) :: reject_text
character(6) :: Nthreshold_text
character(12) :: prefix_text

!Trajectory Variables
real(dp) :: trajectory_CPU_time,trajectory_wall_time
real(dp) :: scattering_angle
real(dp),dimension(3) :: TRVenergies1,TRVenergies2,dTRVenergies
real(dp) :: totalEnergy
integer :: header1_old,header2_old,header3_old

!Trajectory Initialization Variables
real(dp) :: random_num1,random_num2,random_num3,random_r2,random_r3,i,j
real(dp) :: initial_bond_distance, initial_rotation_angle, initial_rotational_speed
real(dp) :: initial_bond_angle1, initial_bond_angle2
real(dp) :: initial_energy_H2,initial_vibrational_energy,initial_rotational_energy
real(dp) :: bond_period_elapsed
real(dp) :: probJ_max, J_factor3

!Trajectory Output
real(dp),dimension(3,Natoms) :: coords_initial, velocities_initial
real(dp),dimension(3,Natoms) :: coords_final, velocities_final

!Timing Variables
real :: r1,r2
integer :: seed,c1,c2,cr
real :: system_clock_rate
integer :: grid_t0, grid_t1, trajectory_t0, trajectory_t1
real :: grid_wall_time,checktrajectory_wall_time
character(10) :: grid_wall_time_text, checktrajectory_wall_time_text
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
write(FMTinitial,FMT="(A17)") "("//Nbond_text//"(6(F8.4)))"
write(FMTtimeslice,FMT="(A19)") "("//Natom_text//"(12(F12.7)))"
write(FMT2,FMT="(A22)") "("//Natom_text//"(6(1x,F14.10)))"
write(FMT3,FMT="(A22)") "("//Natom_text//"(3(1x,F14.10)))"  

!Initialize the clock
call system_clock(count_rate=cr)
system_clock_rate = 1.0/real(cr)

!We now do some formatting for the names of the files (for the occasional checks)
write(Nthreshold_text,FMT=FMT6_pos_real0) threshold_rmsd
if (reject_flag) then
	reject_text = "reject"
else
	reject_text = "accept"
end if
prefix_text = reject_text//Nthreshold_text

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
        header2 = 1
	header2_old = 1
        header3 = 1
	header3_old = 1

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
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,m,random_num1,random_num2,random_num3,random_r2,random_r3,&
	!				       initial_bond_angle1,initial_bond_angle2,initial_bond_distance,&
	!				       bond_period_elapsed,initial_vibrational_energy,J_factor3,probJ_max,&
	!				       initial_rotational_angle,initial_rotational_speed,&
	!				       INITIAL_BOND_DATA,steps,&
	!				       coords_initial,velocities_initial,coords_final,velocities_final,&
	!				       trajectory_CPU_time,trajectory_wall_time,r1,r2,c1,c2)
	!$OMP DO





!!! TEST !!!
!call addMultipleTrajectories()
!exit
!!!!!!!!!!!!




        do n = 1, Ntraj_max

		!Get some random initial conditions for the trajectory
		!For now, we are only handling systems with H-H bonds

		do m = 1, Nbonds

			!The orientation of the H2 should be some random
			!point in the unit sphere
			do
				!First get a random point in the unit cube
				!centered at zero
				random_num1 = rand() - 0.5d0
				random_num2 = rand() - 0.5d0
				random_num3 = rand() - 0.5d0
				random_r2 = random_num1**2 + random_num2**2
				random_r3 = random_r2 + random_num3**2

				!If the point lies outside of the cube, reject it
				if (random_r3 > 0.25d0) cycle
				random_r2 = sqrt(random_r2)

				!But if it lies in the sphere, use its direction (angles)
				initial_bond_angle1 = atan2(random_num1,random_num2)
				initial_bond_angle2 = atan2(random_r2,random_num3)
				exit
			end do

 			!The vibrational energy of the H2 should be some random value
			!that follows the boltzmann distribution at this temperature
			do
				!This picks a random value between zero and some very high upper limit
				random_num1 = rand() * upsilon_max
!				random_num2 = rand() * upsilon_factor2
				random_num2 = rand()

!				if (exp(-random_num1 * upsilon_factor1) < random_num2) cycle
				if (exp(-random_num1 * upsilon_factor1) * upsilon_factor2 < random_num2) cycle

				initial_vibrational_energy = (random_num1 + 0.5d0) * epsilon_factor
				exit
			end do
			initial_bond_distance = HOr0_hydrogen + sqrt(initial_vibrational_energy*2/HOke_hydrogen)
			J_factor3 = J_factor1 / (initial_bond_distance**2)
			probJ_max = sqrt(2*J_factor3) * exp(J_factor3*0.25d0 - 0.5d0)
!			J_factor3 = 85.3d0 / temperature

 			!The rotational energy of the H2 should be some random value
			!that follows the boltzmann distribution at this temperature
			do
				!This picks a random value between zero and some very high upper limit
				random_num1 = rand() * J_max
				random_num2 = rand() * probJ_max

				if ((2*random_num1 + 1.0d0) * J_factor3 * exp(-random_num1 * (random_num1 + 1.0d0) * &
                                    J_factor3) < random_num2) cycle

				initial_rotational_energy = (random_num1) * (random_num1 + 1.0d0) * J_factor2
				exit
			end do

			random_num1 = rand()
			initial_rotational_speed = sqrt(initial_rotational_energy/mass_hydrogen)
			initial_rotation_angle = random_num1*pi2 - pi
			bond_period_elapsed = rand()

			!All of this is stored for later use in the InitialSetup of runTrajectory
			INITIAL_BOND_DATA(:,m) = (/ initial_bond_distance,initial_rotational_speed,&
	                                        initial_rotation_angle,initial_bond_angle1,initial_bond_angle2,&
                                                bond_period_elapsed /)
		end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		GRID CREATION MONITORING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!This big if-statement is if we want to monitor our grid creation
		!The check only happens every Ngrid_check
		!Right now it is spaced so that we get (at most) ten graphs over the period of its creation
                if ((.true.).and.(modulo(n,Ngrid_check) == 0)) then

			!Remark: output of checkTrajectory is in the checkstatefile
			call system_clock(trajectory_t0)
!			call checkTrajectory(scattering_angle,TRVenergies1,TRVenergies1)
	                call checkTrajectory(coords_initial,velocities_initial,coords_final,velocities_final)
			call system_clock(trajectory_t1)
			checktrajectory_wall_time = (trajectory_t1-trajectory_t0)*system_clock_rate

			open(gnuplotchannel,file=gridpath1//gnuplotfile)
			write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
			write(checkstateTrajectory,FMT=FMT6_pos_int) Ntraj
			write(gnuplotchannel,*) 'set output "'//gridpath1//'checkTrajectory_'//prefix_text//&
			                         '_'//checkstateTrajectory//'.png"'
			write(gnuplotchannel,*) 'set style line 1 lc rgb "red" pt 5'
			write(gnuplotchannel,*) 'set style line 2 lc rgb "green" pt 7'
			write(gnuplotchannel,*) 'set style line 3 lc rgb "blue" pt 13'
			write(gnuplotchannel,*) 'set style line 4 lc rgb "orange" pt 9'
			write(gnuplotchannel,*) 'set style line 5 lc rgb "yellow" pt 11'
			write(gnuplotchannel,*) 'set style line 6 lc rgb "pink" pt 20'
			write(gnuplotchannel,*) 'set format x ""'
!			write(gnuplotchannel,*) 'unset xtics'
			write(gnuplotchannel,*) 'set tmargin 0'
			write(gnuplotchannel,*) 'set bmargin 0'
			write(gnuplotchannel,*) 'set lmargin 1'
			write(gnuplotchannel,*) 'set rmargin 1'
			write(gnuplotchannel,*) 'set multiplot layout 6,1 margins 0.15,0.95,.1,.95 spacing 0,0 title '//&
						'"Trajectory '//trim(adjustl(checkstateTrajectory))//'"'
			write(gnuplotchannel,*) 'unset key'
			write(gnuplotchannel,*) 'unset xlabel'
!			write(angle1descriptor,FMT=FMT6_neg_real1) initial_bond_angle1
!			write(angle2descriptor,FMT=FMT6_pos_real1) initial_bond_angle2
!			write(bond1descriptor,FMT=FMT6_pos_real1) initial_bond_distance
 			write(checktrajectory_wall_time_text,FMT="(F10.2)") checktrajectory_wall_time
!			write(gnuplotchannel,*) 'set label 1 "H2 Orientation: '//angle1descriptor//', '//angle2descriptor//&
!                                                ' radians" at screen 0.6, 0.955'
!			write(gnuplotchannel,*) 'set label 2 "H2 Bond Length: '//bond1descriptor//&
!                                                ' A" at screen 0.6, 0.94'
 			write(gnuplotchannel,*) 'set label 3 "Total Wall Time: '//checktrajectory_wall_time_text//&
                                                 ' s" at screen 0.6, 0.910'
			write(gnuplotchannel,*) 'set ylabel "Var1 (A)"'
			write(gnuplotchannel,*) 'set yrange [0:',max_var1,']'
			write(gnuplotchannel,*) 'set ytics 2'
			write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//'" u 4:7 w lines'
			write(gnuplotchannel,*) 'set ylabel "Var2 (A)"'
			write(gnuplotchannel,*) 'set yrange [0:',max_var2,']'
			write(gnuplotchannel,*) 'set ytics 2'
			write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//'" u 4:8 w lines'
!			write(gnuplotchannel,*) 'set ylabel "Total Energy (eV)"'
!			write(gnuplotchannel,*) 'set yrange [0:0.02]'
!			write(gnuplotchannel,*) 'set ytics 0.005'
!			write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//'" u 4:($9+$10) w lines'
			write(gnuplotchannel,*) 'set ylabel "Number of Frames Checked"'
			write(gnuplotchannel,*) 'set autoscale y'
			write(gnuplotchannel,*) 'set ytics autofreq'
			write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//'" u 4:1 w lines'
!			write(gnuplotchannel,*) 'unset label 1'
!			write(gnuplotchannel,*) 'unset label 2'
 			write(gnuplotchannel,*) 'unset label 3'
			write(gnuplotchannel,*) 'set ylabel "Order of Cell Checked"'
			write(gnuplotchannel,*) 'set yrange [-1.1:1.1]'
			write(gnuplotchannel,*) 'set ytics 1'
			write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//'" u 4:2 w lines'
			write(gnuplotchannel,*) 'set ylabel "Number of Cells Checked"'
			write(gnuplotchannel,*) 'set yrange [0:5]'
			write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//'" u 4:3 w lines'
!			write(gnuplotchannel,*) 'set xtics'
			write(gnuplotchannel,*) 'set ylabel "Timestep RMSD (A)"'
			write(gnuplotchannel,*) 'set xlabel "Timestep"'
			write(gnuplotchannel,*) 'set format x'
!			write(gnuplotchannel,*) 'set yrange [0:.2002]'
			write(gnuplotchannel,*) 'set autoscale y'
			write(gnuplotchannel,*) 'set logscale y'
			write(gnuplotchannel,*) 'set ytics (".1" .1, ".05" .05, ".01" .01, ".001" .001, ".0001" .0001)'
			write(gnuplotchannel,*) 'unset key'
			write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//'" u 4:5 w lines'
			close(gnuplotchannel)
			call system(path_to_gnuplot//"gnuplot < "//gridpath1//gnuplotfile)
			
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
		open(filechannel1,file=gridpath1//"maxframesorder1.dat",position="append")
		write(filechannel1,FMT=*) Ntraj, maxval(counter1), overcrowd1
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
	write(gnuplotchannel,*) 'set tmargin 0'
	write(gnuplotchannel,*) 'set bmargin 0'
	write(gnuplotchannel,*) 'set lmargin 1'
	write(gnuplotchannel,*) 'set rmargin 1'
	write(gnuplotchannel,*) 'set multiplot layout 5,1 margins 0.15,0.95,.1,.99 spacing 0,0 title "Creation of Grid '//Ngrid_text//'/"'
	write(gnuplotchannel,*) 'unset key'
	write(gnuplotchannel,*) 'unset xlabel'
	write(gnuplotchannel,*) 'set ylabel "Number of Files"'
	write(gnuplotchannel,*) 'plot "'//gridpath1//"Initial"//informaticsfile//'" u 3:6 w lines'
	write(gnuplotchannel,*) 'set ylabel "Number of Calls to DivyUp"'
	write(gnuplotchannel,*) 'plot "'//gridpath1//"Initial"//informaticsfile//'" u 3:4 w lines, '//&
	                           '"'//gridpath1//"Initial"//informaticsfile//'" u 3:5 w lines'
	write(gnuplotchannel,*) 'set ylabel "Percentage of Frames added to Order 1"'
	write(gnuplotchannel,*) 'set yrange [0:100]'
	write(gnuplotchannel,*) 'plot "'//gridpath1//"Initial"//informaticsfile//'" u 3:7 w lines'
	write(gnuplotchannel,*) 'set autoscale y'
	write(gnuplotchannel,*) 'set ylabel "Wall Time (sec)"'
	write(grid_wall_time_text,FMT="(F10.2)") grid_wall_time
	write(gnuplotchannel,*) 'set label 1 "Total Wall Time (including grid checking): '//&
                                trim(adjustl(grid_wall_time_text))//' s" at graph 0.025,0.9'
	write(gnuplotchannel,*) 'plot "'//gridpath1//"Initial"//informaticsfile//'" u 3:2 w lines'
	write(gnuplotchannel,*) 'set xtics'
	write(gnuplotchannel,*) 'set xlabel "Trajectories"'
	write(gnuplotchannel,*) 'set autoscale y'
	write(gnuplotchannel,*) 'unset label 1'
	write(gnuplotchannel,*) 'set ylabel "CPU Time (sec)"'
	write(gnuplotchannel,*) 'plot "'//gridpath1//"Initial"//informaticsfile//'" u 3:1 w lines'
	close(gnuplotchannel)
	call system(path_to_gnuplot//"gnuplot < "//gridpath1//gnuplotfile)
	
	
	open(gnuplotchannel,file=gridpath1//gnuplotfile)
	write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
	write(gnuplotchannel,*) 'set output "'//gridpath1//'MaxFramesOfOrder1.png"'
	write(gnuplotchannel,*) 'set style line 1 lc rgb "red" pt 5'
	write(gnuplotchannel,*) 'unset key'
	write(gnuplotchannel,*) 'set xlabel "Trajectories"'
	write(gnuplotchannel,*) 'set ylabel "Max Number of Frames in Order 1 Cells (Subcells)"'
	write(gnuplotchannel,*) 'plot "'//gridpath1//'maxframesorder1.dat" u 1:2 w lines'!,\'
!	write(gnuplotchannel,*) '     "'//gridpath1//'maxframesorder1.dat" u 1:3 w lines'
	close(gnuplotchannel)
	call system(path_to_gnuplot//"gnuplot < "//gridpath1//gnuplotfile)
	
	
	
	
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

