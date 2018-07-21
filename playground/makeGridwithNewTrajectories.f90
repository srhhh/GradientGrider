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
use analyzeScatteringAngleswithMultipleGrids
implicit none

!Grid Directory/File Formatting Strings
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(6) :: reject_text
character(6) :: Nthreshold_text

!Trajectory Variables
real(dp) :: trajectory_CPU_time,trajectory_wall_time
real(dp) :: speedH,speedH2,scattering_angle
real(dp),dimension(3) :: velocityH1,velocityH2
integer :: header1_old,header2_old,header3_old

!Trajectory Initialization Variables
real(dp) :: random_num1,random_num2,random_num3,random_r2,random_r3,i,j
real(dp) :: initial_bond_distance, initial_rotation_angle, initial_rotational_speed
real(dp) :: initial_bond_angle1, initial_bond_angle2
real(dp) :: initial_energy_H2,initial_vibrational_energy,initial_rotational_energy

!Timing Variables
real :: r1,r2
integer :: seed,c1,c2,cr
real :: system_clock_rate

!Incremental Integers
integer :: n, m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		LIBRARY (MULTI-GRID) INITIALIZATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Get a random seed and print it in case there's a problem you need to replicate
call system_clock(seed)
print *, ""
print *, "System clock seed: ", seed
seed = rand(seed)

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		GRID INITIALIZATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do Ngrid = 1, Ngrid_max

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

	!Whenever a folder is made, make a print statement
	print *, gridpath1
	print *, gridpath2

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

 			!The energy of the H2 should be some random value
			!that follows the boltzmann distribution at this temperature

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			WORK ON THIS LATER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			do
				!This picks a random value between zero and some very high upper limit
				random_num1 = rand()
				random_num2 = rand()
				initial_energy_H2 = (upsilon_max*random_num1 + 0.5d0)*upsilon_factor2

				!We calculate the probability of this value occuring in our distribution [0,1]
				!and we accept it if our second random number (also in range [0,1]) is below it
				if (random_num2 < temperature_scaling*exp(upsilon_max*random_num1*upsilon_factor1)) exit
			end do

			!The ratio of vib:rot energy of the H2
			!Right now, we set it to all vibrational
			random_num2 = 1.0d0
			initial_vibrational_energy = random_num2*initial_energy_H2
			initial_rotational_energy = initial_energy_H2 - initial_vibrational_energy

			!Right now, we make all vibrational energy be potential energy
			!by simply increasing the bond length to match the required energy
			initial_bond_distance = HOr0_hydrogen + sqrt(initial_vibrational_energy*2/HOke_hydrogen)
			random_num1 = rand()
			! RS: Please use the momentum of inertia and angular velocity to define rotation
			!	KF: must work on this in the future
			initial_rotational_speed = sqrt(initial_rotational_energy/mass_hydrogen)
			initial_rotation_angle = random_num1*pi2
	
			!All of this is stored for later use in the InitialSetup of runTrajectory
			INITIAL_BOND_DATA(Nbonds,:) = (/ initial_bond_distance,initial_rotational_speed,&
	                                        initial_rotation_angle,initial_bond_angle1,initial_bond_angle2 /)
		end do

                open(progresschannel,file=gridpath1//progressfile,position="append")
                write(progresschannel,*) ""
                write(progresschannel,*) ""
                write(progresschannel,*) "Starting trajectory ", Ntraj+1
                write(progresschannel,*) "  Initial Conditions: "
                write(progresschannel,*) "  		Bond Distance: ", initial_bond_distance
                write(progresschannel,*) "  		 Bond Angle 1: ", initial_bond_angle1
                write(progresschannel,*) "  		 Bond Angle 2: ", initial_bond_angle2
                close(progresschannel)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		GRID CREATION MONITORING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!This big if-statement is if we want to monitor our grid creation
		!The check only happens every Ngrid_check
		!Right now it is spaced so that we get (at most) ten graphs over the period of its creation
                if (modulo(n,Ngrid_check) == 0) then

			!Remark: ScatteringAngles2 checks the trajectoriesfile INSIDE the Ngrid/ subdirectory
			call getScatteringAngles2(Ngrid_text//"/Initial"//trajectoriesfile,8,9,10,&
                                                  "InitialScatteringAngleDistribution_"&
			                          //Ngrid_text//reject_text//Nthreshold_text)

			!Remark: output of checkTrajectory is in the checkstatefile
			call checkTrajectory(velocityH1,velocityH2)

			!We also record the scattering angle of the trajectory
			speedH = sqrt(velocityH1(1)**2 + velocityH1(2)**2 + velocityH1(3)**2)
			speedH2 = sqrt(velocityH2(1)**2 + velocityH2(2)**2 + velocityH2(3)**2)
			scattering_angle = acos(dot_product(velocityH1,velocityH2) / &
						           (speedH * speedH2))
			
			write(checkstateTrajectory,FMT=FMT6_pos_int) Ntraj
			open(gnuplotchannel,file=gridpath1//gnuplotfile)
			write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
			write(gnuplotchannel,*) 'set output "'//gridpath1//'checkTrajectory_'//reject_text//Nthreshold_text//&
			                         checkstateTrajectory//'.jpg"'
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
			write(gnuplotchannel,*) 'set multiplot layout 4,1 margins 0.15,0.95,.1,.95 spacing 0,0 title '//&
						'"Trajectory '//trim(adjustl(checkstateTrajectory))//'of '//gridpath0//'"'
			write(gnuplotchannel,*) 'unset key'
			write(gnuplotchannel,*) 'unset xlabel'
			write(angle1descriptor,FMT=FMT6_neg_real1) initial_bond_angle1
			write(angle2descriptor,FMT=FMT6_pos_real1) initial_bond_angle2
			write(bond1descriptor,FMT=FMT6_pos_real1) initial_bond_distance
			write(scatteringdescriptor,FMT=FMT6_pos_real1) scattering_angle
			write(gnuplotchannel,*) 'set label 1 "H2 Orientation: '//angle1descriptor//', '//angle2descriptor//&
                                                ' radians" at screen 0.7, 0.955'
			write(gnuplotchannel,*) 'set label 2 "H2 Bond Length: '//bond1descriptor//&
                                                ' A" at screen 0.7, 0.94'
			write(gnuplotchannel,*) 'set label 3 "Scattering Angle: '//scatteringdescriptor//&
                                                ' rand" at screen 0.7, 0.925'
			write(gnuplotchannel,*) 'set ylabel "Var1 (A)"'
			write(gnuplotchannel,*) 'set yrange [0:11]'
			write(gnuplotchannel,*) 'set ytics 2'
			write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//'" u 4:7 w lines'
			write(gnuplotchannel,*) 'unset label 1'
			write(gnuplotchannel,*) 'unset label 2'
			write(gnuplotchannel,*) 'unset label 3'
			write(gnuplotchannel,*) 'set ylabel "Var2 (A)"'
			write(gnuplotchannel,*) 'set yrange [0:11]'
			write(gnuplotchannel,*) 'set ytics 2'
			write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//'" u 4:8 w lines'
			write(gnuplotchannel,*) 'set ylabel "Total Energy (eV)"'
			write(gnuplotchannel,*) 'set yrange [0:0.02]'
			write(gnuplotchannel,*) 'set ytics 0.005'
			write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//'" u 4:($9+$10) w lines'
			write(gnuplotchannel,*) 'set xtics'
			write(gnuplotchannel,*) 'set ylabel "Timestep RMSD (A)"'
			write(gnuplotchannel,*) 'set xlabel "Timestep"'
			write(gnuplotchannel,*) 'set yrange [0:.2002]'
			write(gnuplotchannel,*) 'set ytics .02'
			write(gnuplotchannel,*) 'unset key'
			write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//'" u 4:5 w lines'
			close(gnuplotchannel)
			call system("gnuplot < "//gridpath1//gnuplotfile)
			
                end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		TRAJECTORY ADDITION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!We time how much time each trajectory takes, wall-time and CPU time
                call CPU_time(r1)
                call system_clock(c1)
                call addTrajectory(velocityH1,velocityH2)

                call CPU_time(r2)
                call system_clock(c2)
                trajectory_CPU_time = r2 - r1
                trajectory_wall_time = (c2 -c1) * system_clock_rate

		!If there have been a large number of subdivisions (so many that our array will go
		!out of bounds) then we stop; we also stop if it is taking too long
                if ((header1 == header1_max).or.&
                   (trajectory_CPU_time > trajectory_CPU_time_max)) exit

		!We also record the scattering angle of the trajectory
		speedH = sqrt(velocityH1(1)**2 + velocityH1(2)**2 + velocityH1(3)**2)
		speedH2 = sqrt(velocityH2(1)**2 + velocityH2(2)**2 + velocityH2(3)**2)
		scattering_angle = acos(dot_product(velocityH1,velocityH2) / &
					           (speedH * speedH2))

		!Otherwise, we consider this a successful trajectory addition
                Ntraj = Ntraj + 1
                open(progresschannel,file=gridpath1//progressfile,position="append")
                write(progresschannel,*) "Finished trajectory ", Ntraj
                write(progresschannel,*) "        Now we have ", Nfile, " files"
                write(progresschannel,*) "                      Wall Time: ",&
                                                trajectory_wall_time
                write(progresschannel,*) "                       CPU Time: ",&
                                                trajectory_CPU_time
                write(progresschannel,*) "               Scattering Angle: ",&
                                                scattering_angle
                close(progresschannel)

		!This is all recorded in the trajectoriesfile of the grid
                open(filechannel1,file=gridpath1//"Initial"//trajectoriesfile,position="append")
                write(filechannel1,*) Ntraj, header1-header1_old, header2-header2_old, Nfile,&
                                      trajectory_CPU_time/real(steps),trajectory_wall_time/real(steps),&
                                      Norder1*100.0/real(steps),&
				      scattering_angle,initial_bond_angle1,initial_bond_angle2
                close(filechannel1)
		header1_old = header1
		header2_old = header2
		header3_old = header3
        end do

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
	write(progresschannel,*) ""
	write(progresschannel,*) ""
	close(progresschannel)
	
	
	!Now that we are done with a grid, we can see how much time the grid creation took
	!There is also some other interesting data
	open(gnuplotchannel,file=gridpath1//gnuplotfile)
	write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
	write(gnuplotchannel,*) 'set output "'//gridpath1//'GridCreationGraph.jpg"'
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
	write(gnuplotchannel,*) 'set multiplot layout 5,1 margins 0.15,0.95,.1,.99 spacing 0,0 title "Trajectory '//Ngrid_text//'"'
	write(gnuplotchannel,*) 'unset key'
	write(gnuplotchannel,*) 'unset xlabel'
	write(gnuplotchannel,*) 'set ylabel "Number of Files"'
	write(gnuplotchannel,*) 'plot "'//gridpath1//"Initial"//trajectoriesfile//'" u 1:4 w lines'
	write(gnuplotchannel,*) 'set ylabel "Number of Calls to DivyUp"'
	write(gnuplotchannel,*) 'plot "'//gridpath1//"Initial"//trajectoriesfile//'" u 1:2 w lines, '//&
	                           '"'//gridpath1//"Initial"//trajectoriesfile//'" u 1:3 w lines'
	write(gnuplotchannel,*) 'set ylabel "Percentage of Frames added to Order 1"'
	write(gnuplotchannel,*) 'set yrange [0:100]'
	write(gnuplotchannel,*) 'plot "'//gridpath1//"Initial"//trajectoriesfile//'" u 1:7 w lines'
	write(gnuplotchannel,*) 'set autoscale y'
	write(gnuplotchannel,*) 'set ylabel "Wall Time (sec)"'
	write(gnuplotchannel,*) 'plot "'//gridpath1//"Initial"//trajectoriesfile//'" u 1:6 w lines'
	write(gnuplotchannel,*) 'set xtics'
	write(gnuplotchannel,*) 'set xlabel "Trajectories"'
	write(gnuplotchannel,*) 'set autoscale y'
	write(gnuplotchannel,*) 'set ylabel "CPU Time (sec)"'
	write(gnuplotchannel,*) 'plot "'//gridpath1//"Initial"//trajectoriesfile//'" u 1:5 w lines'
	close(gnuplotchannel)
	call system("gnuplot < "//gridpath1//gnuplotfile)
	
	!Also, make a scattering angle plot
	call getScatteringAngles2(Ngrid_text//"/"//trajectoriesfile,8,9,10,"InitialScatteringAngleDistribution_"&
	                          //Ngrid_text//reject_text//Nthreshold_text)
	
	
	
	
	
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
	
	
	
end do

end program makeGridwithNewTrajectories
