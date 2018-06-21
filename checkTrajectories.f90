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
character(3) :: checkstateDescriptor
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
open(filechannel1,file=trim(path4)//"counter0.txt",status="old")
do n = 1, counter0_max
        read(filechannel1,FMT="(I8)") counter0(n)
end do
close(filechannel1)

open(filechannel1,file=trim(path4)//"counter1.txt",status="old")
do n = 1, counter1_max
        read(filechannel1,FMT="(I8)") counter1(n)
end do
close(filechannel1)

open(filechannel1,file=trim(path4)//"counter2.txt",status="old")
do n = 1, counter2_max
        read(filechannel1,FMT="(I8)") counter2(n)
end do
close(filechannel1)

open(filechannel1,file=trim(path4)//"counter3.txt",status="old")
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

i = sqrt(-2.0*log(rand()))
j = pi2*rand()
random_num1 = i*cos(j)
random_num2 = i*sin(j)
initial_energy_H2 = 0.5*mass_hydrogen*(velocity_scaling*random_num1)**2

random_num2 = rand()
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
call addTrajectory(initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
		   initial_bond_angle1,initial_bond_angle2,&
		   header1,header2,header3,counter0,counter1,counter2,counter3)

write(checkstateDescriptor,FMT="(I0.3)") Ntraj
open(filechannel1,file=path4//temporaryfile2)
write(filechannel1,*) 'set term jpeg'
write(filechannel1,*) 'set output "'//path4//'checkstate_newoutput_'//checkstateDescriptor//'.jpg"'
write(filechannel1,*) 'unset xtics'
write(filechannel1,*) 'set tmargin 0'
write(filechannel1,*) 'set bmargin 0'
write(filechannel1,*) 'set lmargin 1'
write(filechannel1,*) 'set rmargin 1'
write(filechannel1,*) 'set multiplot layout 3,1 margins 0.15,0.95,.1,.99 spacing 0,0'
write(filechannel1,*) 'unset key'
write(filechannel1,*) 'unset xlabel'
write(filechannel1,*) 'set ylabel "Var1 (A)"'
write(filechannel1,*) 'set yrange [0:40]'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 1:3 w lines'
write(filechannel1,*) 'set ylabel "Var2 (A)"'
write(filechannel1,*) 'unset key'
write(filechannel1,*) 'set yrange [0:40]'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 1:4 w lines'
write(filechannel1,*) 'set xtics'
write(filechannel1,*) 'set ylabel "Minimum RMSD (A)"'
write(filechannel1,*) 'set xlabel "Timestep"'
write(filechannel1,*) 'set yrange [0:2.0e-1]'
write(filechannel1,*) 'unset key'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 1:2 w lines'
close(filechannel1)
call system("gnuplot < "//path4//temporaryfile2)
call system("rm "//path4//temporaryfile2)

end do


end program checkTrajectories
