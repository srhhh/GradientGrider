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
real,parameter :: initial_translational_KE = (1.0)*eV/RU_energy
                                        !originally eV
real,parameter :: collision_distance = sqrt(cutoff_distance_squared)
real,parameter :: collision_skew = HOr0_hydrogen*0.0

real :: initial_bond_distance, initial_rotation_angle, initial_rotational_speed
real :: initial_bond_angle1, initial_bond_angle2
real :: initial_energy_H2,initial_vibrational_energy,initial_rotational_energy
integer :: seed,n,Ntraj,c1,c2,cr
real :: random_num1,random_num2,i,j,system_clock_rate

integer :: header1 = 1
integer :: header2 = 1
integer :: header3 = 1
integer, dimension(counter0_max) :: counter0 = 0
integer, dimension(counter1_max) :: counter1 = 0
integer, dimension(counter2_max) :: counter2 = 0
integer, dimension(counter3_max) :: counter3 = 0


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
seed = 699
seed = rand(seed)

!For timing purposes
call system_clock(c1,count_rate=cr)
system_clock_rate = real(cr)



!Initialize a trajectory randomly
random_num1 = rand()
random_num2 = rand()
initial_bond_angle1 = random_num1*2*pi			!theta
initial_bond_angle2 = acos(2*random_num2 - 1.0)		!phi

i = sqrt(-2*log(rand()))
j = 2*pi*rand()
random_num1 = i*cos(j)
initial_energy_H2 = 0.5*temperature_constant*random_num1**2

random_num2 = rand()
random_num2 = 1.0
initial_vibrational_energy = random_num2*initial_energy_H2
initial_rotational_energy = initial_energy_H2 - initial_vibrational_energy

initial_bond_distance = HOr0_hydrogen + sqrt(initial_vibrational_energy*2/HOke_hydrogen)

random_num1 = rand()
initial_rotational_speed = sqrt(initial_rotational_energy/mass_hydrogen)
initial_rotation_angle = random_num1*2*pi

!Remark: module makeTrajectory2 is an augmented version that does not
!	 add a state to the grid or counters, so does not mess with the grid;
!	 it also produces a graph that tracks the closest frame's rmsd over time
call addTrajectory(initial_translational_KE,collision_distance,collision_skew,&
		   initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
		   initial_bond_angle1,initial_bond_angle2,&
		   header1,header2,header3,counter0,counter1,counter2,counter3)






end program checkTrajectories
