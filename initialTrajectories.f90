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

open(progresschannel,file=path4//progressfile)
write(progresschannel,*) ""
write(progresschannel,*) ""
write(progresschannel,*) "Let's Start!"
write(progresschannel,*) ""
close(progresschannel)

! RS: I would use a less "psedo-" algorithm
!		KF: note to self: make sure to use completely random
!                                 if I want an unbiased space
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
seed = 669
seed = rand(seed)

!For timing purposes
! RS: Did you figure out the system_clock? What is it reporting?
!		KF: note to self: save this comment
call system_clock(c1,count_rate=cr)
system_clock_rate = real(cr)


!For 100 trajectories, for starters

do Ntraj = 1, 100

!The orientation of the H2 bond relative to center of mass vector* (CMV) follows a uniform distribution
! * CMV is the vector connectiong the center of mass of both molecules
!Because of the uneven griding of a sphere by theta and phi, phi needs to have
!been picked from a normal distribution
random_num1 = rand()
random_num2 = rand()
initial_bond_angle1 = random_num1*pi2			!theta
!		KF: note: the second angle should follow an inverse trig function curve
!		    note-to-self: do uniform distribution for now
!                                 switch to real distribution later
initial_bond_angle2 = random_num2*pi2		!phi

!A Box-Muller calculation is setup to generate normal distribution (will use generated number later)
i = sqrt(-2*log(rand()))
j = 2*pi*rand()
random_num1 = i*cos(j)
random_num2 = i*sin(j)

! RS: the center of mass translational energy of H2 is assumed to be zero
! RS: the collision energy is the kinetic energy of H
! RS: the average vibrational energy of H2 is assumed to be xxx???
! RS: this energy is stored in the form of potential energy
! RS: the rotation energy of H2 is assume to be zero
!		KF: this I still need to perfect

!First we get a Gaussian distribution by using Box-Muller and multiply the
!number N by sqrt(kT/m) to get a velocity
!Then simply calculate the KE as 0.5*m*v**2
initial_energy_H2 = 0.5*temperature_constant*random_num1**2
! RS: How about 
! RS: initial_energy_H2 = temperature_constant*random_num1


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

call addTrajectory(initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
		   initial_bond_angle1,initial_bond_angle2,&
		   header1,header2,header3,counter0,counter1,counter2,counter3)




call system_clock(c2)
open(progresschannel,file=path4//progressfile,position="append")
write(progresschannel,*) ""
write(progresschannel,*) "Trajectory ", Ntraj, " finished in ", (c2-c1)/system_clock_rate
close(progresschannel)
c1 = c2

end do


!Finally, write the counters to a text file

open(filechannel1,file=trim(path4)//"counter0.txt")
do n = 1, counter0_max
        write(filechannel1,FMT="(I8)") counter0(n)
end do
close(filechannel1)

open(filechannel1,file=trim(path4)//"counter1.txt")
do n = 1, counter1_max
        write(filechannel1,FMT="(I8)") counter1(n)
end do
close(filechannel1)

open(filechannel1,file=trim(path4)//"counter2.txt")
do n = 1, counter2_max
        write(filechannel1,FMT="(I8)") counter2(n)
end do
close(filechannel1)

open(filechannel1,file=trim(path4)//"counter3.txt")
do n = 1, counter3_max
        write(filechannel1,FMT="(I8)") counter3(n)
end do
close(filechannel1)


end program initialTrajectories
