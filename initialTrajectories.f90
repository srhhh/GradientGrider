program initialTrajectories

!I call all the modules here just to keep track of them for the makefile
use makeTrajectory
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

!The random number generator is setup
!This is a uniform distribution
seed = 699
seed = rand(seed)

!For timing purposes
call system_clock(c1,count_rate=cr)
system_clock_rate = real(cr)


!For 200 trajectories, for starters
do Ntraj = 1, 150




!The angle of the H2 bond relative to collision follows a uniform distribution
!Because of the uneven griding of a sphere by theta and phi, phi needs to have
!been picked from a normal distribution
random_num1 = rand()
random_num2 = rand()
initial_bond_angle1 = random_num1*2*pi			!theta
initial_bond_angle2 = acos(2*random_num2 - 1.0)		!phi

!A Box-Muller calculation is setup
!This is a normal distribution (will use generated number later)
i = sqrt(-2*log(rand()))
j = 2*pi*rand()
random_num1 = i*cos(j)


!The initial energy of H2 should follow a Boltzmann distribution
!First we get a Gaussian distribution by using Box-Muller and multiply the
!number N by sqrt(kT/m) to get a velocity
!Then simply calculate the KE as 0.5*m*v**2
initial_energy_H2 = 0.5*temperature_constant*random_num1**2

!This energy can be distributed as either vibrational or rotational
!Assume that both modes are equally likely (uniform distribution)
random_num2 = rand()

!for testing, set it to be all vibration or all rotation
random_num2 = 1.0
initial_vibrational_energy = random_num2*initial_energy_H2
initial_rotational_energy = initial_energy_H2 - initial_vibrational_energy

!Vibrational energy is reflected in the bond distance and
!movement of the hydrogens ANTIPARALLEL to each other and
!PARALLEL to the intermolecular axis
!For simplicity sake, assume at initial time, no KE, all PE
initial_bond_distance = HOr0_hydrogen + sqrt(initial_vibrational_energy*2/HOke_hydrogen)

!Rotational energy is reflected in the the movement of the hydrogens
!ANTIPARALLEL to each other and PERPENDICULAR to the intermolecular axis
!There is no rotational PE to my knowledge
random_num1 = rand()
initial_rotational_speed = sqrt(initial_rotational_energy/mass_hydrogen)
initial_rotation_angle = random_num1*2*pi

call addTrajectory(initial_translational_KE,collision_distance,collision_skew,&
		   initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
		   initial_bond_angle1,initial_bond_angle2,&
		   header1,header2,header3,counter0,counter1,counter2,counter3)




call system_clock(c2)
print *, ""
print *, "Trajectory ", Ntraj, " finished in ", (c2-c1)/system_clock_rate
c1 = c2

end do


!Finally, write the counters to a text file

open(filechannel1,file=trim(path4)//"counter0.txt",status="old")
do n = 1, counter0_max
        write(filechannel1,FMT="(I8)") counter0(n)
end do
close(filechannel1)

open(filechannel1,file=trim(path4)//"counter1.txt",status="old")
do n = 1, counter1_max
        write(filechannel1,FMT="(I8)") counter1(n)
end do
close(filechannel1)

open(filechannel1,file=trim(path4)//"counter2.txt",status="old")
do n = 1, counter2_max
        write(filechannel1,FMT="(I8)") counter2(n)
end do
close(filechannel1)

open(filechannel1,file=trim(path4)//"counter3.txt",status="old")
do n = 1, counter3_max
        write(filechannel1,FMT="(I8)") counter3(n)
end do
close(filechannel1)


end program initialTrajectories
