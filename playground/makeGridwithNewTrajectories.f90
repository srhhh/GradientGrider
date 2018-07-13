
program makeGridwithNewTrajectories

!I call all the modules here just to keep track of them for the makefile
use addNewTrajectorytoGrid
use checkNewTrajectorywithGrid
use checkGrid
use mapCellData
use addFrametoGrid
use PARAMETERS
use FUNCTIONS
use VARIABLES
use ANALYSIS
use PHYSICSFUNCTIONS
use analyzeScatteringAngleswithMultipleGrids
implicit none

!COLLISION PARAMETERS
real(dp) :: initial_bond_distance, initial_rotation_angle, initial_rotational_speed
real(dp) :: initial_bond_angle1, initial_bond_angle2
real(dp) :: initial_energy_H2,initial_vibrational_energy,initial_rotational_energy
real(dp) :: random_num1,random_num2,random_num3,random_r2,random_r3,i,j
integer :: seed,n,m

!Counter Parameters
integer :: header1 = 1
integer :: header2 = 1
integer :: header3 = 1
integer, dimension(counter0_max) :: counter0 = 0
integer, dimension(counter1_max) :: counter1 = 0
integer, dimension(counter2_max) :: counter2 = 0
integer, dimension(counter3_max) :: counter3 = 0

!Grid Naming Parameters
integer :: Ngrid
character(5) :: variable_length_text
character(resolution_text_length) :: resolution_text
character(overcrowd0_text_length) :: overcrowd0_text
character(trajectory_text_length) :: trajectory_text
character(Ngrid_text_length) :: Ngrid_text
character(gridpath_length) :: gridpath0
character(gridpath_length+len(Ngrid_text)+1) :: gridpath1
character(gridpath_length+len(Ngrid_text)+1+5) :: gridpath2
character(6) :: reject_text
character(6) :: Nthreshold_text

!Plot Naming Parameters
character(6) :: checkstateTrajectory
character(6) :: angle1descriptor,angle2descriptor,bond1descriptor,scatteringdescriptor

!Trajectory Parameters
real(dp) :: trajectory_CPU_time,trajectory_wall_time
integer :: Ntraj,Nfile,steps,Norder1
real(dp) :: speedH,speedH2,scattering_angle
real(dp),dimension(3) :: velocityH1,velocityH2

!Timing Parameters
real :: r1,r2
integer :: c1,c2,cr
real :: system_clock_rate

!Get a random seed and print it in case there's a problem
call system_clock(seed)
!seed = 96967942
print *, ""
print *, "System clock seed: ", seed
seed = rand(seed)

!Initialize the clock
call system_clock(count_rate=cr)
system_clock_rate = 1.0/real(cr)

write(Nthreshold_text,FMT="(F6.5)") threshold_rmsd
if (reject_flag) then
	reject_text = "reject"
else
	reject_text = "accept"
end if

!Make the multi-grid folder!
write(variable_length_text,FMT="(I5)") resolution_text_length
write(resolution_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") resolution_0
write(variable_length_text,FMT="(I5)") overcrowd0_text_length
write(overcrowd0_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") overcrowd0
write(variable_length_text,FMT="(I5)") trajectory_text_length
write(trajectory_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ntraj_max
gridpath0 = path5//resolution_text//"_"//overcrowd0_text//"_"//trajectory_text//"/"
call system("mkdir "//gridpath0)
call system("cp "//path5//parametersfile//" "//gridpath0//parametersfile)

!We start off with zero files
Nfile = 0

!For each grid, just spit out Ntraj_max trajectories into it
do Ngrid = 1, Ngrid_max

!Counter Parameters
header1 = 1
header2 = 1
header3 = 1
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

	!The grids will be named 001 with increments of 1
	write(variable_length_text,FMT="(I5)") Ngrid_text_length
        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
        gridpath1 = gridpath0//Ngrid_text//"/"
        gridpath2 = gridpath1//"grid/"

	!This will create a folder for the grid files and for its outputs
        call system("mkdir "//gridpath1)
        call system("mkdir "//gridpath2)

print *, gridpath1
print *, gridpath2

        open(progresschannel,file=gridpath1//progressfile)
        write(progresschannel,*) ""
        write(progresschannel,*) ""
        write(progresschannel,*) "Let's Start!"
        write(progresschannel,*) ""
        close(progresschannel)

	!We start off with zero trajectories
        Ntraj = 0
        do n = 1, Ntraj_max

!Get some random initial conditions for the trajectory
!The orientation of the H2
do
random_num1 = rand() - 0.5d0
random_num2 = rand() - 0.5d0
random_num3 = rand() - 0.5d0
random_r2 = random_num1**2 + random_num2**2
random_r3 = random_r2 + random_num3**2
if (random_r3 > 0.25d0) cycle
random_r2 = sqrt(random_r2)
initial_bond_angle1 = acos(random_num1 / random_r2)
initial_bond_angle2 = atan2(random_r2,random_num3)
exit
end do
!The energy of the H2
do
random_num1 = rand()
random_num2 = rand()
initial_energy_H2 = (upsilon_max*random_num1 + 0.5d0)*upsilon_factor2
if (random_num2 < temperature_scaling*exp(upsilon_max*random_num1*upsilon_factor1)) exit
end do
!The ratio of vib:rot energy of H2
random_num2 = 1.0d0
initial_vibrational_energy = random_num2*initial_energy_H2
initial_rotational_energy = initial_energy_H2 - initial_vibrational_energy
initial_bond_distance = HOr0_hydrogen + sqrt(initial_vibrational_energy*2/HOke_hydrogen)
random_num1 = rand()
initial_rotational_speed = sqrt(initial_rotational_energy/mass_hydrogen)
initial_rotation_angle = random_num1*pi2

                open(progresschannel,file=gridpath1//progressfile,position="append")
                write(progresschannel,*) ""
                write(progresschannel,*) ""
                write(progresschannel,*) "Starting trajectory ", Ntraj+1
                write(progresschannel,*) "  Initial Conditions: "
                write(progresschannel,*) "  		Bond Distance: ", initial_bond_distance
                write(progresschannel,*) "  		 Bond Angle 1: ", initial_bond_angle1
                write(progresschannel,*) "  		 Bond Angle 2: ", initial_bond_angle2
                close(progresschannel)


		!This big if-statement is if we want to monitor our grid creation
                if (modulo(n,Ntraj_max/10) == 0) then

			write(variable_length_text,FMT="(I5)") Ngrid_text_length
			write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
			call getScatteringAngles2(gridpath1,trajectoriesfile,8,9,10,"InitialScatteringAngleDistribution_"&
			                          //Ngrid_text//reject_text//Nthreshold_text)

			call checkTrajectory(initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
                           		     initial_bond_angle1,initial_bond_angle2,.false.,&
                           		     header1,header2,header3,counter0,counter1,counter2,counter3,gridpath2,&
					     velocityH1,velocityH2)

			!We also record the scattering angle of the trajectory
			speedH = sqrt(velocityH1(1)**2 + velocityH1(2)**2 + velocityH1(3)**2)
			speedH2 = sqrt(velocityH2(1)**2 + velocityH2(2)**2 + velocityH2(3)**2)
			scattering_angle = acos(dot_product(velocityH1,velocityH2) / &
						           (speedH * speedH2))

write(checkstateTrajectory,FMT="(I0.6)") Ntraj
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
			'"Trajectory '//checkstateTrajectory//'of '//gridpath0//'"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'unset xlabel'
write(angle1descriptor,FMT="(F6.3)") initial_bond_angle1
write(angle2descriptor,FMT="(F6.3)") initial_bond_angle2
write(bond1descriptor,FMT="(F6.4)") initial_bond_distance
write(scatteringdescriptor,FMT="(F6.3)") scattering_angle
write(gnuplotchannel,*) 'set label 1 "H2 Orientation: '//angle1descriptor//', '//angle2descriptor//' radians" at screen 0.7, 0.955'
write(gnuplotchannel,*) 'set label 2 "H2 Bond Length: '//bond1descriptor//' A" at screen 0.7, 0.94'
write(gnuplotchannel,*) 'set label 3 "Scattering Angle: '//scatteringdescriptor//' rand" at screen 0.7, 0.925'
write(gnuplotchannel,*) 'set ylabel "Var1 (A)"'
write(gnuplotchannel,*) 'set yrange [0:11]'
write(gnuplotchannel,*) 'set ytics 2'
write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//'" u 4:7 w lines'
write(gnuplotchannel,*) 'unset label 1'
write(gnuplotchannel,*) 'unset label 2'
write(gnuplotchannel,*) 'unset label 3'
!write(gnuplotchannel,*) 'set ylabel "Var1 Deviance (A)"'
!write(gnuplotchannel,*) 'set yrange [0:1.0]'
!write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//'" u 4:($11-$7) w lines'
write(gnuplotchannel,*) 'set ylabel "Var2 (A)"'
write(gnuplotchannel,*) 'set yrange [0:11]'
write(gnuplotchannel,*) 'set ytics 2'
write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//'" u 4:8 w lines'
!write(gnuplotchannel,*) 'set ylabel "Var2 Deviance (A)"'
!write(gnuplotchannel,*) 'set yrange [0:1.0]'
!write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//'" u 4:($12-$8) w lines'
!write(gnuplotchannel,*) 'set ylabel "Total RMSD (A)"'
!write(gnuplotchannel,*) 'set yrange [0:1.0e0]'
!write(gnuplotchannel,*) 'set ytics .2'
!write(gnuplotchannel,*) 'plot "'//gridpath1//checkstatefile//'" u 4:13 w lines'
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

		!We time how much time each trajectory takes, wall-time and CPU time
                call CPU_time(r1)
                call system_clock(c1)
                call addTrajectory(initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
		   initial_bond_angle1,initial_bond_angle2,&
		   header1,header2,header3,counter0,counter1,counter2,counter3,&
		   steps,Nfile,gridpath2,Norder1,velocityH1,velocityH2)

                call CPU_time(r2)
                call system_clock(c2)
                trajectory_CPU_time = r2 - r1
                trajectory_wall_time = (c2 -c1) * system_clock_rate

		!If there have been a large number of subdivisions (so many that our array will go
		!out of bounds) then we stop; we also stop if it is taking too long
                if ((header1 == header1_max).or.&
                   (trajectory_CPU_time > trajectory_CPU_time_max)) exit

		!Otherwise, we consider this a successful trajectory addition
                Ntraj = Ntraj + 1
                open(progresschannel,file=gridpath1//progressfile,position="append")
                write(progresschannel,*) "Finished trajectory ", Ntraj
                write(progresschannel,*) "        Now we have ", Nfile, " files"
                write(progresschannel,*) "                      Wall Time: ",&
                                                trajectory_wall_time
                write(progresschannel,*) "                       CPU Time: ",&
                                                trajectory_CPU_time
                close(progresschannel)

		!We also record the scattering angle of the trajectory
		speedH = sqrt(velocityH1(1)**2 + velocityH1(2)**2 + velocityH1(3)**2)
		speedH2 = sqrt(velocityH2(1)**2 + velocityH2(2)**2 + velocityH2(3)**2)
		scattering_angle = acos(dot_product(velocityH1,velocityH2) / &
					           (speedH * speedH2))

		!This is all recorded in the trajectoriesfile of the grid
                open(filechannel1,file=gridpath1//trajectoriesfile,position="append")
                write(filechannel1,*) Ntraj, header1, header2, Nfile,&
                                      trajectory_CPU_time,trajectory_wall_time,Norder1*1.0/real(steps),&
				      scattering_angle,initial_bond_angle1,initial_bond_angle2
                close(filechannel1)
        end do


open(progresschannel,file=gridpath1//progressfile,position="append")
write(progresschannel,*) ""
write(progresschannel,*) ""
write(progresschannel,*) "Finished all trajectories for grid "//Ngrid_text
write(progresschannel,*) "        Now we have ", Nfile, " files"
write(progresschannel,*) ""
write(progresschannel,*) ""
close(progresschannel)


!Here, we can see how much time the grid creation took
!There is also some other interesting data
open(gnuplotchannel,file=gridpath1//gnuplotfile)
write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath1//'CPUTime.jpg"'
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
write(gnuplotchannel,*) 'plot "'//gridpath1//trajectoriesfile//'" u 1:4 w lines'
write(gnuplotchannel,*) 'set ylabel "Number of Overcrowded Cells"'
write(gnuplotchannel,*) 'plot "'//gridpath1//trajectoriesfile//'" u 1:2 w lines, '//&
                           '"'//gridpath1//trajectoriesfile//'" u 1:3 w lines'
write(gnuplotchannel,*) 'set ylabel "Fraction of Frames added to Order 1"'
write(gnuplotchannel,*) 'set yrange [0:1.0]'
write(gnuplotchannel,*) 'plot "'//gridpath1//trajectoriesfile//'" u 1:7 w lines'
write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set ylabel "Wall Time (sec)"'
write(gnuplotchannel,*) 'plot "'//gridpath1//trajectoriesfile//'" u 1:6 w lines'
write(gnuplotchannel,*) 'set xtics'
write(gnuplotchannel,*) 'set xlabel "Trajectories"'
write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set ylabel "CPU Time (sec)"'
write(gnuplotchannel,*) 'plot "'//gridpath1//trajectoriesfile//'" u 1:5 w lines'
close(gnuplotchannel)

write(variable_length_text,FMT="(I5)") Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
call getScatteringAngles2(gridpath1,trajectoriesfile,8,9,10,"InitialScatteringAngleDistribution_"&
                          //Ngrid_text//reject_text//Nthreshold_text)






open(filechannel1,file=trim(gridpath1)//counter0file)
do n = 1, counter0_max
        write(filechannel1,FMT="(I8)") counter0(n)
end do
close(filechannel1)

open(filechannel1,file=trim(gridpath1)//counter1file)
do n = 1, counter1_max
        write(filechannel1,FMT="(I8)") counter1(n)
end do
close(filechannel1)

open(filechannel1,file=trim(gridpath1)//counter2file)
do n = 1, counter2_max
        write(filechannel1,FMT="(I8)") counter2(n)
end do
close(filechannel1)

open(filechannel1,file=trim(gridpath1)//counter3file)
do n = 1, counter3_max
        write(filechannel1,FMT="(I8)") counter3(n)
end do
close(filechannel1)

end do

end program makeGridwithNewTrajectories
