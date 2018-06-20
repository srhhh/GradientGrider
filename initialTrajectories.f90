! RS: This program does what

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
! RS: Move thhere stuff to f2_physics_parameters.f90 since they are parameters.

! relative translation energy eV between H and H2 in reduced unit
real,parameter :: initial_translational_KE = (1.0)*eV/RU_energy
! cut-off distance to halt the simulation after collision
! RS: please make this equal to the initial separation in reduced unit
real,parameter :: collision_distance = sqrt(cutoff_distance_squared)
! optimal H-H bond length in reduced unit?
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

! RS: I would use a less "psedo-" algorithm
! RS: the following code generates 3 random numbers at the same time
! RS: but I am sure you can make it generate one
!246 ! uniform distribution random number
!247       subroutine uniform_rand(x,n)
!248       implicit none
!249 
!250       integer, allocatable :: seed(:)
!251       integer :: un, istat, n, i
!252       real :: x(3,n)
!253 
!254       call random_seed(size = i)
!255       allocate(seed(i))
!256 
!257       open(newunit=un, file="/dev/urandom", access="stream", &
!258              form="unformatted", action="read", status="old", &
!259              iostat=istat)
!260       if (istat == 0) then
!261          read(un) seed
!262 !         write(6,*) seed
!263          close(un)
!264       endif
!265 
!266       CALL random_seed(put=seed)
!267       CALL RANDOM_NUMBER(x)
!268       return
!269       end subroutine

!The random number generator is setup
!This is a uniform distribution
seed = 699
seed = rand(seed)

!For timing purposes
! RS: Did you figure out the system_clock? What is it reporting?
call system_clock(c1,count_rate=cr)
system_clock_rate = real(cr)


!For 200 trajectories, for starters
! RS: 150?
do Ntraj = 1, 150

!The orientation of the H2 bond relative to center of mass vector* (CMV) follows a uniform distribution
! * CMV is the vector connectiong the center of mass of both molecules
!Because of the uneven griding of a sphere by theta and phi, phi needs to have
!been picked from a normal distribution
random_num1 = rand()
random_num2 = rand()
! RS: make 2*pi a constant 
! RS: again, be MORE CAREFUL with variable type
! RS: 2*pi -> 2.0*pi, you never want int*real unless you are certain about it!
initial_bond_angle1 = random_num1*2*pi			!theta
! RS: why don't just generate another random angle like you did for initial_bond_angle1?
initial_bond_angle2 = acos(2*random_num2 - 1.0)		!phi

!A Box-Muller calculation is setup to generate normal distribution (will use generated number later)
! RS: What are the types of i and j?? integer?? real?? not defined?? 
! RS: I abusoltely have no idea what the next three lines do.
! RS: It is definitely not Box-Muller. Box-Muller generates a pair of random numbers per time.
i = sqrt(-2*log(rand()))
j = 2*pi*rand()
random_num1 = i*cos(j)


! RS: the center of mass translational energy of H2 is assumed to be zero
! RS: the average vibrational energy of H2 is assumed to be RT/2
! RS: this energy is stored in the form of potential energy
! RS: the rotation energy of H2 is assume to be zero
!The initial energy of H2 should follow a Boltzmann distribution
!First we get a Gaussian distribution by using Box-Muller and multiply the
!number N by sqrt(kT/m) to get a velocity
!Then simply calculate the KE as 0.5*m*v**2
initial_energy_H2 = 0.5*temperature_constant*random_num1**2

! RS: this is not ture -- due to the partition function the energy is distributated
! RS: very unevenly
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
! RS: you mean away from each other?
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
