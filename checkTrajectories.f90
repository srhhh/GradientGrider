program checkTrajectories

!I call all the modules here just to keep track of them for the makefile
use makeTrajectory2
use addCells5
use f2_parameters
use f1_functions
use f2_physics_parameters
use f2_variables
implicit none

!COLLISION PARAMETERS
real :: initial_bond_distance, initial_rotation_angle, initial_rotational_speed
real :: initial_bond_angle1, initial_bond_angle2
real :: initial_energy_H2,initial_vibrational_energy,initial_rotational_energy
integer :: seed,n,Ntraj,c1,c2,cr
character(50) :: checkstateDescriptor
character(3) :: checkstateTrajectory
character(6) :: angle1descriptor,angle2descriptor,bond1descriptor
real :: random_num1,random_num2,i,j,system_clock_rate

integer :: header1 = 1
integer :: header2 = 1
integer :: header3 = 1
integer, dimension(counter0_max) :: counter0 = 0
integer, dimension(counter1_max) :: counter1 = 0
integer, dimension(counter2_max) :: counter2 = 0
integer, dimension(counter3_max) :: counter3 = 0


open(progresschannel,file=path4//progressfile)
write(progresschannel,*) ""
write(progresschannel,*) ""
write(progresschannel,*) "Let's Start!"
write(progresschannel,*) ""
close(progresschannel)



!First, write the counters from their text files
open(filechannel1,file=trim(path4)//counter0file,status="old")
do n = 1, counter0_max
        read(filechannel1,FMT="(I8)") counter0(n)
end do
close(filechannel1)

open(filechanneannel1,*) 'set xtics'
write(filechannel1,*) 'set ylabel "Minimum RMSD (A)"'
write(filechannel1,*) 'set xlabel "Timestep"'
write(filechannel1,*) 'set yrange [0:2.0e-1]'
write(filechannel1,*) 'unset key'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:5 w lines'
close(filechannel1)
call system("gnuplot < "//path4//temporaryfile2)
call system("rm "//path4//temporaryfile2)

1,file=trim(path4)//counter1file,status="old")
do i = 1, counter1_max
        read(filechannel1,FMT="(I8)") counter1(n)
end do
close(filechannel1)

open(filechannel1,file=trim(path4)//counter2file,status="old")
do n = 1, counter2_max
        read(filechannel1,FMT="(I8)") counter2(n)
end do
close(filechannel1)

open(filechannel1,file=trim(path4)//counter3file,status="old")
do n = 1, counter3_max
        read(filechannel1,FMT="(I8)") counter3(n)
end do
close(filechannel1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		TRAJECTORY CHECKING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!The random number generator is setup
!This is a uniform distribution
seed = 69
seed = rand(seed)

!For timing purposes
call system_clock(c1,count_rate=cr)
system_clock_rate = real(cr)




do Ntraj = 1, 12

!Initialize a trajectory randomly
random_num1 = rand()
random_num2 = rand()
initial_bond_angle1 = random_num1*pi2			!theta
initial_bond_angle2 = random_num2*pi2			!phi

do
random_num1 = rand()
random_num2 = rand()

initial_energy_H2 = (upsilon_max*random_num1 + 0.5)*upsilon_factor2
if (random_num2 < temperature_scaling*exp(upsilon_max*random_num1*upsilon_factor1)) exit

end do

random_num2 = 1.0
initial_vibrational_energy = random_num2*initial_energy_H2
initial_rotational_energy = initial_energy_H2 - initial_vibrational_energy

initial_bond_distance = HOr0_hydrogen + sqrt(initial_vibrational_energy*2.0/HOke_hydrogen)

random_num1 = rand()
initial_rotational_speed = sqrt(initial_rotational_energy/mass_hydrogen)
initial_rotation_angle = random_num1*pi2

!Remark: module makeTrajectory2 is an augmented version that does not
!	 add a state to the grid or counters, so does not mess with the grid;
!	 it also produces a graph that tracks the closest frame's rmsd over time
call checkTrajectory(initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
		     initial_bond_angle1,initial_bond_angle2,&
		     header1,header2,header3,counter0,counter1,counter2,counter3)

write(checkstateDescriptor,FMT="(A50)") path3(len(path2)+1:len(path3)-1)
write(checkstateTrajectory,FMT="(I0.3)") Ntraj
open(filechannel1,file=path4//temporaryfile2)
write(filechannel1,*) 'set term jpeg size 1200,1200'
write(filechannel1,*) 'set output "'//path4//'checkstate_'//trim(adjustl(checkstateDescriptor))//'_'//checkstateTrajectory//'.jpg"'
write(filechannel1,*) 'set style line 1 lc rgb "red" pt 5'
write(filechannel1,*) 'set style line 2 lc rgb "green" pt 7'
write(filechannel1,*) 'set style line 3 lc rgb "blue" pt 9'
write(filechannel1,*) 'set style line 4 lc rgb "yellow" pt 14'
write(filechannel1,*) 'unset xtics'
write(filechannel1,*) 'set tmargin 0'
write(filechannel1,*) 'set bmargin 0'
write(filechannel1,*) 'set lmargin 1'
write(filechannel1,*) 'set rmargin 1'
write(filechannel1,*) 'set multiplot layout 6,1 margins 0.15,0.95,.1,.99 spacing 0,0'
write(filechannel1,*) 'unset key'
write(filechannel1,*) 'unset xlabel'
write(angle1descriptor,FMT="(F6.4)") initial_bond_angle1
write(angle2descriptor,FMT="(F6.4)") initial_bond_angle2
write(bond1descriptor,FMT="(F6.4)") initial_bond_distance
write(filechannel1,*) 'set label 1 "H2 Orientation: '//angle1descriptor//', '//angle2descriptor//' radians" at screen 0.7, 0.975'
write(filechannel1,*) 'set label 2 "H2 Bond Length: '//bond1descriptor//' A" at screen 0.7, 0.955'
write(filechannel1,*) 'set ylabel "Var1 (A)"'
write(filechannel1,*) 'set yrange [0:40]'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:6 w lines'
write(filechannel1,*) 'unset label 1'
write(filechannel1,*) 'unset label 2'
write(filechannel1,*) 'set ylabel "Var2 (A)"'
write(filechannel1,*) 'set yrange [0:40]'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:7 w lines'
write(filechannel1,*) 'set ylabel "Neighbor Check?"'
write(filechannel1,*) 'unset ytics'
write(filechannel1,*) 'set yrange [-0.2:1.2]'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:3 w lines, '//&
                           '"'//path4//checkstatefile//'" u 4:($2==0?$3:4/0) w p ls 3, '//&
                           '"'//path4//checkstatefile//'" u 4:(($3==1)&&($2<2)?$3:4/0) w p ls 1, '//&
                           '"'//path4//checkstatefile//'" u 4:(($3==1)&&($2>1)?$3:4/0) w p ls 2'
write(filechannel1,*) 'set ylabel "Granularity (order)"'
write(filechannel1,*) 'set ytics 1'
write(filechannel1,*) 'set yrange [-0.5:3.5]'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:2 w lines, '//&
                           '"'//path4//checkstatefile//'" u 4:($2==0?$2:4/0) w p ls 3, '//&
                           '"'//path4//checkstatefile//'" u 4:(($3==1)&&($2<2)?$2:4/0) w p ls 1, '//&
                           '"'//path4//checkstatefile//'" u 4:(($3==1)&&($2>1)?$2:4/0) w p ls 2'
write(filechannel1,*) 'set ylabel "Frames Checked"'
write(filechannel1,*) 'set ytics auto'
write(filechannel1,*) 'set yrange [0:300]'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:1 w lines, '//&
                           '"'//path4//checkstatefile//'" u 4:($2==0?$1:4/0) w p ls 3, '//&
                           '"'//path4//checkstatefile//'" u 4:(($3==1)&&($2<2)?$1:4/0) w p ls 1, '//&
                           '"'//path4//checkstatefile//'" u 4:(($3==1)&&($2>1)?$1:4/0) w p ls 2'
write(filechannel1,*) 'set xtics'
write(filechannel1,*) 'set ylabel "Minimum RMSD (A)"'
write(filechannel1,*) 'set xlabel "Timestep"'
write(filechannel1,*) 'set yrange [0:2.0e-1]'
write(filechannel1,*) 'unset key'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:5 w lines'
close(filechannel1)
call system("gnuplot < "//path4//temporaryfile2)
call system("rm "//path4//temporaryfile2)

end do


end program checkTrajectories
