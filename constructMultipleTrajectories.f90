
program constructMultipleTrajectories

!I call all the modules here just to keep track of them for the makefile
use makeTrajectory
!use makeTrajectory3
!use checkCells4
!use ls_rmsd
!use mapCellData
use addCells5
use f2_parameters
use f1_functions
use f2_physics_parameters
use f2_variables
use decompose_velocities
implicit none

!COLLISION PARAMETERS
real :: initial_bond_distance, initial_rotation_angle, initial_rotational_speed
real :: initial_bond_angle1, initial_bond_angle2
real :: initial_energy_H2,initial_vibrational_energy,initial_rotational_energy
real :: random_num1,random_num2,i,j
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
integer,parameter :: Ngrid_max = 5
integer,parameter :: resolution_text_length = 4
integer,parameter :: overcrowd0_text_length = 5
integer,parameter :: overcrowd1_text_length = 5
integer,parameter :: Ngrid_text_length = 3
integer,parameter :: gridpath_length = resolution_text_length +&
                                       overcrowd0_text_length +&
                                       overcrowd1_text_length + len(path5) + 3
character(10) :: variable_length_text
character(resolution_text_length) :: resolution_text
character(overcrowd0_text_length) :: overcrowd0_text
character(overcrowd1_text_length) :: overcrowd1_text
character(Ngrid_text_length) :: Ngrid_text
character(gridpath_length) :: gridpath0
character(gridpath_length+len(Ngrid_text)+1) :: gridpath1
character(gridpath_length+len(Ngrid_text)+1+5) :: gridpath2

character(3) :: checkstateTrajectory
character(6) :: angle1descriptor,angle2descriptor,bond1descriptor

!Trajectory Parameters
real :: trajectory_CPU_time_max = 60.0
real :: trajectory_CPU_time,trajectory_wall_time
integer,parameter :: Ntraj_max = 500
integer :: Ntraj,Nfile,Norder1
real :: speedH,speedH2,scattering_angle
real,dimension(3) :: velocityH,velocityH2

!Timing Parameters
real :: r1,r2
integer :: c1,c2,cr
real :: system_clock_rate

!Initialize some stuff
seed = 669
 call system_clock(seed)
seed = rand(seed)
call system_clock(count_rate=cr)
system_clock_rate = 1.0/real(cr)
write(resolution_text,FMT="(I0.4)") resolution_0
write(overcrowd0_text,FMT="(I0.5)") overcrowd0
write(overcrowd1_text,FMT="(I0.5)") overcrowd1
gridpath0 = path5//resolution_text//"_"//overcrowd0_text//"_"//overcrowd1_text//"/"
call system("mkdir "//gridpath0)
call system("cp "//path2//"f2_parameters.f90 "//gridpath0//parametersfile)

!open(filechannel1,file=gridpath0//parametersfile,position="append")
!write(variable_length_text,FMT="(I10)") Ngrid_max
!write(filechannel1,*) "integer,parameter :: Ngrid_max = "//trim(adjustl(variable_length_text))
!write(variable_length_text,FMT="(I10)") Ngrid_text_length
!write(filechannel1,*) "integer,parameter :: Ngrid_text_length = "//trim(adjustl(variable_length_text))
!write(variable_length_text,FMT="(I10)") gridpath_length
!write(filechannel1,*) "integer,parameter :: gridpath_length = "//trim(adjustl(variable_length_text))
!write(filechannel1,*) "character(gridpath_length),parameter :: gridpath0 = "//gridpath0
!close(filechannel1)

Nfile = 0

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

        write(Ngrid_text,FMT="(I0.3)") Ngrid
        gridpath1 = gridpath0//Ngrid_text//"/"
        gridpath2 = gridpath1//"grid/"

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

        Ntraj = 0
        do n = 1, Ntraj_max

!Get some random initial conditions for the trajectory
random_num1 = rand()
random_num2 = rand()
initial_bond_angle1 = random_num1*pi2           !theta
initial_bond_angle2 = random_num2*pi2           !phi
do
random_num1 = rand()
random_num2 = rand()
initial_energy_H2 = (upsilon_max*random_num1 + 0.5)*upsilon_factor2
if (random_num2 < temperature_scaling*exp(upsilon_max*random_num1*upsilon_factor1)) exit
end do
random_num2 = 1.0
initial_vibrational_energy = random_num2*initial_energy_H2
initial_rotational_energy = initial_energy_H2 - initial_vibrational_energy
initial_bond_distance = HOr0_hydrogen + sqrt(initial_vibrational_energy*2/HOke_hydrogen)
random_num1 = rand()
initial_rotational_speed = sqrt(initial_rotational_energy/mass_hydrogen)
initial_rotation_angle = random_num1*2*pi

                if (.false.) then !(modulo(Ntraj,10) == 9) then

!                       call monitorTrajectory(initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
!                                  initial_bond_angle1,initial_bond_angle2,.false.,&
!                                  header1,header2,header3,counter0,counter1,counter2,counter3,gridpath2)

write(checkstateTrajectory,FMT="(I0.3)") Ntraj
open(filechannel1,file=path4//temporaryfile2)
write(filechannel1,*) 'set term jpeg size 1200,1200'
write(filechannel1,*) 'set output "'//gridpath1//'monitorTrajectory_'//checkstateTrajectory//'.jpg"'
write(filechannel1,*) 'set style line 1 lc rgb "red" pt 5'
write(filechannel1,*) 'set style line 2 lc rgb "green" pt 7'
write(filechannel1,*) 'set style line 3 lc rgb "blue" pt 13'
write(filechannel1,*) 'set style line 4 lc rgb "orange" pt 9'
write(filechannel1,*) 'set style line 5 lc rgb "yellow" pt 11'
write(filechannel1,*) 'set style line 6 lc rgb "pink" pt 20'
write(filechannel1,*) 'unset xtics'

if (.true.) then
write(filechannel1,*) 'set tmargin 0'
write(filechannel1,*) 'set bmargin 0'
write(filechannel1,*) 'set lmargin 1'
write(filechannel1,*) 'set rmargin 1'
write(filechannel1,*) 'set multiplot layout 6,1 margins 0.15,0.95,.1,.99 spacing 0,0 title '//&
			'"Trajectory '//checkstateTrajectory//'of '//gridpath0(len(path3)+1:len(gridpath0))//'"'
write(filechannel1,*) 'unset key'
write(filechannel1,*) 'unset xlabel'
write(angle1descriptor,FMT="(F6.4)") initial_bond_angle1
write(angle2descriptor,FMT="(F6.4)") initial_bond_angle2
write(bond1descriptor,FMT="(F6.4)") initial_bond_distance
write(filechannel1,*) 'set label 1 "H2 Orientation: '//angle1descriptor//', '//angle2descriptor//' radians" at screen 0.7, 0.975'
write(filechannel1,*) 'set label 2 "H2 Bond Length: '//bond1descriptor//' A" at screen 0.7, 0.955'
write(filechannel1,*) 'set ylabel "True Var1 (A)"'
write(filechannel1,*) 'set yrange [0:21]'
write(filechannel1,*) 'set ytics 5'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:7 w lines'
write(filechannel1,*) 'unset label 1'
write(filechannel1,*) 'unset label 2'
write(filechannel1,*) 'set ylabel "Var1 Deviance (A)"'
write(filechannel1,*) 'set yrange [0:1.0]'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:($11-$7) w lines'
write(filechannel1,*) 'set ylabel "True Var2 (A)"'
write(filechannel1,*) 'set yrange [0:21]'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:8 w lines'
write(filechannel1,*) 'set ylabel "Var2 Deviance (A)"'
write(filechannel1,*) 'set yrange [0:1.0]'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:($12-$8) w lines'
write(filechannel1,*) 'set ylabel "Total RMSD (A)"'
write(filechannel1,*) 'set yrange [0:1.0e0]'
write(filechannel1,*) 'set ytics .2'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:13 w lines'
write(filechannel1,*) 'set xtics'
write(filechannel1,*) 'set ylabel "Timestep RMSD (A)"'
write(filechannel1,*) 'set xlabel "Timestep"'
write(filechannel1,*) 'set yrange [0:1.0e-4]'
write(filechannel1,*) 'set ytics .00002'
write(filechannel1,*) 'unset key'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:14 w lines'
end if

close(filechannel1)
call system("gnuplot < "//path4//temporaryfile2)
call system("rm "//path4//temporaryfile2)

                end if

                call CPU_time(r1)
                call system_clock(c1)
                call addTrajectory(initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
		   initial_bond_angle1,initial_bond_angle2,&
		   header1,header2,header3,counter0,counter1,counter2,counter3,&
		   Nfile,gridpath2,Norder1,velocityH,velocityH2)

                call CPU_time(r2)
                call system_clock(c2)
                trajectory_CPU_time = r2 - r1
                trajectory_wall_time = (c2 -c1) * system_clock_rate

                if ((header1 == header1_max).or.&
                   (trajectory_CPU_time > trajectory_CPU_time_max)) exit

                Ntraj = Ntraj + 1
                open(progresschannel,file=gridpath1//progressfile,position="append")
                write(progresschannel,*) ""
                write(progresschannel,*) "Finished trajectory ", Ntraj
                write(progresschannel,*) "        Now we have ", Nfile, " files"
                write(progresschannel,*) "                      Wall Time: ",&
                                                trajectory_wall_time
                write(progresschannel,*) "                       CPU Time: ",&
                                                trajectory_CPU_time
                close(progresschannel)

		speedH = sqrt(velocityH(1)**2 + velocityH(2)**2 + velocityH(3)**2)
		speedH2 = sqrt(velocityH2(1)**2 + velocityH2(2)**2 + velocityH2(3)**2)
		scattering_angle = acos(dot_product(velocityH,velocityH2) / &
					           (speedH * speedH2))

                open(filechannel1,file=gridpath1//trajectoriesfile,position="append")
                write(filechannel1,*) Ntraj, header1, header2, Nfile,&
                                      trajectory_CPU_time,Norder1*1.0/real(Nsteps),&
				      scattering_angle
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



open(filechannel1,file=gridpath1//temporaryfile2)
write(filechannel1,*) 'set term jpeg size 1200,1200'
write(filechannel1,*) 'set output "'//gridpath0//'CPUTime_'//Ngrid_text//'.jpg"'
write(filechannel1,*) 'set style line 1 lc rgb "red" pt 5'
write(filechannel1,*) 'set style line 2 lc rgb "green" pt 7'
write(filechannel1,*) 'set style line 3 lc rgb "blue" pt 13'
write(filechannel1,*) 'set style line 4 lc rgb "orange" pt 9'
write(filechannel1,*) 'set style line 5 lc rgb "yellow" pt 11'
write(filechannel1,*) 'set style line 6 lc rgb "pink" pt 20'
write(filechannel1,*) 'unset xtics'
write(filechannel1,*) 'set tmargin 0'
write(filechannel1,*) 'set bmargin 0'
write(filechannel1,*) 'set lmargin 1'
write(filechannel1,*) 'set rmargin 1'
write(filechannel1,*) 'set multiplot layout 4,1 margins 0.15,0.95,.1,.99 spacing 0,0 title "Trajectory '//Ngrid_text//'"'
write(filechannel1,*) 'unset key'
write(filechannel1,*) 'unset xlabel'
write(filechannel1,*) 'set ylabel "Number of Files"'
write(filechannel1,*) 'plot "'//gridpath1//trajectoriesfile//'" u 1:4 w lines'
write(filechannel1,*) 'set ylabel "Number of Overcrowded Cells"'
write(filechannel1,*) 'plot "'//gridpath1//trajectoriesfile//'" u 1:2 w lines, '//&
                           '"'//gridpath1//trajectoriesfile//'" u 1:3 w lines'
write(filechannel1,*) 'set ylabel "Percent of Frames added to Order 1"'
write(filechannel1,*) 'set yrange [0:1.0]'
write(filechannel1,*) 'plot "'//gridpath1//trajectoriesfile//'" u 1:6 w lines'
write(filechannel1,*) 'set xtics'
write(filechannel1,*) 'set xlabel "timestep"'
write(filechannel1,*) 'set yrange autoscale'
write(filechannel1,*) 'set ylabel "CPU_time (sec)"'
write(filechannel1,*) 'plot "'//gridpath1//trajectoriesfile//'" u 1:5 w lines'
close(filechannel1)

call system("gnuplot < "//gridpath1//temporaryfile2)
call system("rm "//gridpath1//temporaryfile2)

open(filechannel1,file=gridpath1//temporaryfile2)
write(filechannel1,*) 'set term jpeg size 1200,1200'
write(filechannel1,*) 'set output "'//gridpath0//'ScatteringAngles_'//Ngrid_text//'.jpg"'
write(filechannel1,*) 'set style fill solid 1.0 noborder'
write(filechannel1,*) 'unset key'
write(filechannel1,*) 'bin_width = 0.001'
write(filechannel1,*) 'bin_number(x) = floor(x/bin_width)'
write(filechannel1,*) 'rounded(x) = bin_width * (bin_number(x) + 0.5)'
write(filechannel1,*) 'set xlabel "Scattering Angle"'
write(filechannel1,*) 'set ylabel "Occurence"'
write(filechannel1,*) 'plot "'//gridpath1//trajectoriesfile//'" u (rounded($7)):(7) smooth frequency with boxes'
close(filechannel1)

call system("gnuplot < "//gridpath1//temporaryfile2)
call system("rm "//gridpath1//temporaryfile2)























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

end program constructMultipleTrajectories
