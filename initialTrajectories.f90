! RS: This program does what

program initialTrajectories

!I call all the modules here just to keep track of them for the makefile
use makeTrajectory
use makeTrajectory2
use addCells5
use f2_parameters
use f1_functions
use f2_physics_parameters
use f2_variables
implicit none

!COLLISION PARAMETERS
! RS: Move thhere stuff to f2_physics_parameters.f90 since they are parameters.

! RS: please make this equal to the initial separation in reduced unit
! optimal H-H bond length in reduced unit?
!		KF: I believe you call is 'cross-sectional distance'
!		    the distance from the intermolecular axis
!
!		KF: I changed the parametrization
!		    these are now in f2_physics_parameters

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

character(50) :: checkstateDescriptor
character(3) :: checkstateTrajectory
character(6) :: angle1descriptor,angle2descriptor,bond1descriptor


open(progresschannel,file=path4//progressfile)
write(progresschannel,*) ""
write(progresschannel,*) ""
write(progresschannel,*) "Let's Start!"
write(progresschannel,*) ""
close(progresschannel)

! RS: I would use a less "psedo-" algorithm
! RS: the following code generates 3 random numbers at the same time
! RS: but I am sure you can make it generate one
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
!		KF: I did not study it yet
!		KF: note to self: save this comment
call system_clock(c1,count_rate=cr)
system_clock_rate = real(cr)


do Ntraj = 1, 1000

!The orientation of the H2 bond relative to center of mass vector* (CMV) follows a uniform distribution
! * CMV is the vector connectiong the center of mass of both molecules
!Because of the uneven griding of a sphere by theta and phi, phi needs to have
!been picked from a normal distribution
random_num1 = rand()
random_num2 = rand()
! RS: make 2*pi a constant 
! RS: again, be MORE CAREFUL with variable type
! RS: 2*pi -> 2.0*pi, you never want int*real unless you are certain about it!
!		KF: very much resolved
initial_bond_angle1 = random_num1*pi2			!theta
! RS: why don't just generate another random angle like you did for initial_bond_angle1?
!		KF: note: the second angle should follow an inverse trig function curve
!		    note-to-self: do uniform distribution for now
!                                 switch to real distribution later
initial_bond_angle2 = random_num2*pi2		!phi

!A Box-Muller calculation is setup to generate normal distribution (will use generated number later)
! RS: What are the types of i and j?? integer?? real?? not defined?? 
! RS: I abusoltely have no idea what the next three lines do.
! RS: It is definitely not Box-Muller. Box-Muller generates a pair of random numbers per time.
!		KF: we talked about this.
!		KF: note-to-self: need Box-Muller PAIRS, not just one or the other
!i = sqrt(-2*log(rand()))
!j = 2*pi*rand()
!random_num1 = i*cos(j)
!random_num2 = i*sin(j)

! RS: the center of mass translational energy of H2 is assumed to be zero
! RS: the average vibrational energy of H2 is assumed to be RT/2
! RS: this energy is stored in the form of potential energy
! RS: the rotation energy of H2 is assume to be zero
!		KF: this I still need to perfect
!The initial energy of H2 should follow a Boltzmann distribution
!First we get a Gaussian distribution by using Box-Muller and multiply the
!number N by sqrt(kT/m) to get a velocity
!Then simply calculate the KE as 0.5*m*v**2
do
random_num1 = rand()
random_num2 = rand()

initial_energy_H2 = (upsilon_max*random_num1 + 0.5)*upsilon_factor2
if (random_num2 < temperature_scaling*exp(upsilon_max*random_num1*upsilon_factor1)) exit

end do

! RS: this is not ture -- due to the partition function the energy is distributated
! RS: very unevenly
!This energy can be distributed as either vibrational or rotational
!Assume that both modes are equally likely (uniform distribution)
!random_num2 = rand()

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


if (modulo(Ntraj,10)==0) then

call checkTrajectory(initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
		     initial_bond_angle1,initial_bond_angle2,.false.,&
		     header1,header2,header3,counter0,counter1,counter2,counter3)

write(checkstateDescriptor,FMT="(A50)") path3(len(path2)+1:len(path3)-1)
write(checkstateTrajectory,FMT="(I0.3)") Ntraj
open(filechannel1,file=path4//temporaryfile2)
write(filechannel1,*) 'set term jpeg size 1200,1200'
write(filechannel1,*) 'set output "'//path4//'overcrowdtest_'//&
			trim(adjustl(checkstateDescriptor))//'_'//checkstateTrajectory//'.jpg"'
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
write(filechannel1,*) 'set multiplot layout 6,1 margins 0.15,0.95,.1,.99 spacing 0,0 title "Trajectory ',Ntraj,'of '//&
			trim(adjustl(checkstateDescriptor))//'"'
write(filechannel1,*) 'unset key'
write(filechannel1,*) 'unset xlabel'
write(angle1descriptor,FMT="(F6.4)") initial_bond_angle1
write(angle2descriptor,FMT="(F6.4)") initial_bond_angle2
write(bond1descriptor,FMT="(F6.4)") initial_bond_distance
write(filechannel1,*) 'set label 1 "H2 Orientation: '//angle1descriptor//', '//angle2descriptor//' radians" at screen 0.7, 0.975'
write(filechannel1,*) 'set label 2 "H2 Bond Length: '//bond1descriptor//' A" at screen 0.7, 0.955'
write(filechannel1,*) 'set ylabel "Var1 (A)"'
write(filechannel1,*) 'set yrange [0:21]'
write(filechannel1,*) 'set ytics 5'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:7 w lines'
write(filechannel1,*) 'unset label 1'
write(filechannel1,*) 'unset label 2'
write(filechannel1,*) 'set ylabel "Var2 (A)"'
write(filechannel1,*) 'set yrange [0:21]'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:8 w lines'
write(filechannel1,*) 'set ylabel "Neighbor Check?"'
write(filechannel1,*) 'unset ytics'
write(filechannel1,*) 'set yrange [-0.2:1.2]'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:3 w lines, '//&
                           '"'//path4//checkstatefile//'" u 4:($2==0?$3:4/0) w p ls 3, '//&
                           '"'//path4//checkstatefile//'" u 4:(($3==1)&&($2<2)?$3:4/0) w p ls 1, '//&
!                          '"'//path4//checkstatefile//'" u 4:(($5>0.05)&&($3==0)?$3:4/0) w p ls 4, '//&
!                          '"'//path4//checkstatefile//'" u 4:(($5>0.05)&&($3==1)?$3:4/0) w p ls 5'
                           '"'//path4//checkstatefile//'" u 4:($1==1?$3:4/0) w p ls 4, '//&
                           '"'//path4//checkstatefile//'" u 4:($1==2?$3:4/0) w p ls 5'
write(filechannel1,*) 'set ylabel "Granularity (order)"'
write(filechannel1,*) 'set ytics 1'
write(filechannel1,*) 'set yrange [-0.5:3.5]'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:2 w lines, '//&
                           '"'//path4//checkstatefile//'" u 4:($2==0?$2:4/0) w p ls 3, '//&
                           '"'//path4//checkstatefile//'" u 4:(($3==1)&&($2<2)?$2:4/0) w p ls 1, '//&
                           '"'//path4//checkstatefile//'" u 4:(($3==1)&&($2>1)?$2:4/0) w p ls 2, '//&
!                          '"'//path4//checkstatefile//'" u 4:(($5>0.05)&&($3==0)?$2:4/0) w p ls 4, '//&
!                          '"'//path4//checkstatefile//'" u 4:(($5>0.05)&&($3==1)?$2:4/0) w p ls 5'
                           '"'//path4//checkstatefile//'" u 4:($1==1?$2:4/0) w p ls 4, '//&
                           '"'//path4//checkstatefile//'" u 4:($1==2?$2:4/0) w p ls 5'
write(filechannel1,*) 'set ylabel "Frames Checked"'
write(filechannel1,*) 'set ytics auto'
write(filechannel1,*) 'set yrange [0:75]'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:1 w lines, '//&
                           '"'//path4//checkstatefile//'" u 4:($2==0?$1:4/0) w p ls 3, '//&
                           '"'//path4//checkstatefile//'" u 4:(($3==1)&&($2<2)?$1:4/0) w p ls 1, '//&
                           '"'//path4//checkstatefile//'" u 4:(($3==1)&&($2>1)?$1:4/0) w p ls 2, '//&
!                          '"'//path4//checkstatefile//'" u 4:(($5>0.05)&&($3==0)?$1:4/0) w p ls 4, '//&
!                          '"'//path4//checkstatefile//'" u 4:(($5>0.05)&&($3==1)?$1:4/0) w p ls 5'
                           '"'//path4//checkstatefile//'" u 4:($1==1?$1:4/0) w p ls 4, '//&
                           '"'//path4//checkstatefile//'" u 4:($1==2?$1:4/0) w p ls 5'
write(filechannel1,*) 'set xtics'
write(filechannel1,*) 'set ylabel "Minimum RMSD (A)"'
write(filechannel1,*) 'set xlabel "Timestep"'
write(filechannel1,*) 'set yrange [0:2.0e-1]'
write(filechannel1,*) 'unset key'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 4:5 w lines, '//&
!                          '"'//path4//checkstatefile//'" u 4:(($3==1)&&($2<2)?$5:4/0) w p ls 1, '//&
!                          '"'//path4//checkstatefile//'" u 4:(($3==1)&&($2>1)?$5:4/0) w p ls 2, '//&
!                          '"'//path4//checkstatefile//'" u 4:(($5>0.05)&&($3==0)?$5:4/0) w p ls 4, '//&
!                          '"'//path4//checkstatefile//'" u 4:(($5>0.05)&&($3==1)?$5:4/0) w p ls 5'
                           '"'//path4//checkstatefile//'" u 4:($1==1?$5:4/0) w p ls 4, '//&
                           '"'//path4//checkstatefile//'" u 4:($1==2?$5:4/0) w p ls 5'
else
 write(filechannel1,*) 'set multiplot layout 2,1 margins 0.15,0.95,.1,.99 spacing 0,0 title "Trajectory ',Ntraj,'of '//&
 			trim(adjustl(checkstateDescriptor))//' RMSD checks vs RMSD"'
write(filechannel1,*) 'unset key'
write(angle1descriptor,FMT="(F6.4)") initial_bond_angle1
write(angle2descriptor,FMT="(F6.4)") initial_bond_angle2
write(bond1descriptor,FMT="(F6.4)") initial_bond_distance
write(filechannel1,*) 'set label 1 "H2 Orientation: '//angle1descriptor//', '//angle2descriptor//' radians" at screen 0.7, 0.975'
write(filechannel1,*) 'set label 2 "H2 Bond Length: '//bond1descriptor//' A" at screen 0.7, 0.955'
write(filechannel1,*) 'unset xtics'
write(filechannel1,*) 'set ylabel "Minimum RMSD Deviance of CheckState Method (CM) vs. Force Neighbor Check Method (FNCM) (A)"'
write(filechannel1,*) 'unset xlabel'
write(filechannel1,*) 'set autoscale y'
write(filechannel1,*) 'set xrange [0:250]'
write(filechannel1,*) 'unset key'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 1:($5-$6) w points'
!                          '"'//path4//checkstatefile//'" u 4:(($3==1)&&($2<2)?$5:4/0) w p ls 1, '//&
!                          '"'//path4//checkstatefile//'" u 4:(($3==1)&&($2>1)?$5:4/0) w p ls 2, '//&
!                          '"'//path4//checkstatefile//'" u 4:(($5>0.05)&&($3==0)?$5:4/0) w p ls 4, '//&
!                          '"'//path4//checkstatefile//'" u 4:(($5>0.05)&&($3==1)?$5:4/0) w p ls 5'
!                          '"'//path4//checkstatefile//'" u 4:($1==1?$5:4/0) w p ls 4, '//&
!                          '"'//path4//checkstatefile//'" u 4:($1==2?$5:4/0) w p ls 5'
write(filechannel1,*) 'set ylabel "Minimum RMSD of FNCM (A)"'
write(filechannel1,*) 'set xtics'
write(filechannel1,*) 'unset label 1'
write(filechannel1,*) 'unset label 2'
write(filechannel1,*) 'set xlabel "Frames Checked"'
write(filechannel1,*) 'set yrange [0:2.0e-1]'
write(filechannel1,*) 'plot "'//path4//checkstatefile//'" u 1:($5) w points'
!                          '"'//path4//checkstatefile//'" u 4:(($3==1)&&($2<2)?$5:4/0) w p ls 1, '//&
!                          '"'//path4//checkstatefile//'" u 4:(($3==1)&&($2>1)?$5:4/0) w p ls 2, '//&
!                          '"'//path4//checkstatefile//'" u 4:(($5>0.05)&&($3==0)?$5:4/0) w p ls 4, '//&
!                          '"'//path4//checkstatefile//'" u 4:(($5>0.05)&&($3==1)?$5:4/0) w p ls 5'
!                          '"'//path4//checkstatefile//'" u 4:($1==1?$5:4/0) w p ls 4, '//&
!                          '"'//path4//checkstatefile//'" u 4:($1==2?$5:4/0) w p ls 5'
end if


close(filechannel1)
call system("gnuplot < "//path4//temporaryfile2)
call system("rm "//path4//temporaryfile2)

end if

call system_clock(c1)
call addTrajectory(initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
		   initial_bond_angle1,initial_bond_angle2,&
		   header1,header2,header3,counter0,counter1,counter2,counter3)

call system_clock(c2)
open(progresschannel,file=path4//progressfile,position="append")
write(progresschannel,*) ""
write(progresschannel,*) "Trajectory ", Ntraj, " finished in ", (c2-c1)/system_clock_rate
close(progresschannel)

end do


!Finally, write the counters to a text file

open(filechannel1,file=trim(path4)//counter0file)
do n = 1, counter0_max
        write(filechannel1,FMT="(I8)") counter0(n)
end do
close(filechannel1)

open(filechannel1,file=trim(path4)//counter1file)
do n = 1, counter1_max
        write(filechannel1,FMT="(I8)") counter1(n)
end do
close(filechannel1)

open(filechannel1,file=trim(path4)//counter2file)
do n = 1, counter2_max
        write(filechannel1,FMT="(I8)") counter2(n)
end do
close(filechannel1)

open(filechannel1,file=trim(path4)//counter3file)
do n = 1, counter3_max
        write(filechannel1,FMT="(I8)") counter3(n)
end do
close(filechannel1)


end program initialTrajectories
