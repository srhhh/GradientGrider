program checkMultipleTrajectories

!I call all the modules here just to keep track of them for the makefile
use checkNewTrajectorywithGrid
use addFrametoGrid
use PARAMETERS
use FUNCTIONS
use PHYSICS
use VARIABLES
implicit none

!COLLISION PARAMETERS
real :: initial_bond_distance, initial_rotation_angle, initial_rotational_speed
real :: initial_bond_angle1, initial_bond_angle2
real :: initial_energy_H2,initial_vibrational_energy,initial_rotational_energy
character(50) :: checkstateDescriptor
character(3) :: checkstateTrajectory
character(6) :: angle1descriptor,angle2descriptor,bond1descriptor
real :: random_num1,random_num2,i,j

integer :: header1 = 1
integer :: header2 = 1
integer :: header3 = 1
integer, dimension(counter0_max) :: counter0 = 0
integer, dimension(counter1_max) :: counter1 = 0
integer, dimension(counter2_max) :: counter2 = 0
integer, dimension(counter3_max) :: counter3 = 0

!Grid Naming Parameters
integer :: Ngrid
character(resolution_text_length) :: resolution_text
character(overcrowd0_text_length) :: overcrowd0_text
character(overcrowd1_text_length) :: overcrowd1_text
character(Ngrid_text_length) :: Ngrid_text
character(gridpath_length) :: gridpath0
character(gridpath_length+len(Ngrid_text)+1) :: gridpath1
character(gridpath_length+len(Ngrid_text)+1+5) :: gridpath2

!Trajectory Parameters
real :: trajectory_CPU_time,trajectory_wall_time
integer :: seed,n,Ntraj,Nfile

!Timing Parameters
real :: r1,r2
integer :: c1,c2,cr
real :: system_clock_rate

!Initialize some stuff
seed = 669
seed = rand(seed)
call system_clock(count_rate=cr)
system_clock_rate = 1.0/real(cr)
write(resolution_text,FMT="(I0.4)") resolution_0
write(overcrowd0_text,FMT="(I0.5)") overcrowd0
write(overcrowd1_text,FMT="(I0.5)") overcrowd1
gridpath0 = path5//resolution_text//"_"//overcrowd0_text//"_"//overcrowd1_text//"/"
gridpath1 = gridpath0//"001/"


open(progresschannel,file=gridpath1//progressfile)
write(progresschannel,*) ""
write(progresschannel,*) ""
write(progresschannel,*) "Let's Start!"
write(progresschannel,*) ""
close(progresschannel)



!First, write the counters from their text files
open(filechannel1,file=trim(gridpath1)//counter0file,status="old")
do n = 1, counter0_max
        read(filechannel1,FMT="(I8)") counter0(n)
end do
close(filechannel1)

open(filechannel1,file=trim(gridpath1)//counter1file,status="old")
do n = 1, counter1_max
        read(filechannel1,FMT="(I8)") counter1(n)
end do
close(filechannel1)

open(filechannel1,file=trim(gridpath1)//counter2file,status="old")
do n = 1, counter2_max
        read(filechannel1,FMT="(I8)") counter2(n)
end do
close(filechannel1)

open(filechannel1,file=trim(gridpath1)//counter3file,status="old")
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
 call system_clock(seed)
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
		     initial_bond_angle1,initial_bond_angle2,.false.,&
		     header1,header2,header3,counter0,counter1,counter2,counter3,gridpath1)

write(checkstateDescriptor,FMT="(A50)") path3(len(path2)+1:len(path3)-1)
write(checkstateTrajectory,FMT="(I0.3)") Ntraj
open(filechannel1,file=gridpath1//temporaryfile2)
write(filechannel1,*) 'set term jpeg size 1200,1200'
write(filechannel1,*) 'set output "'//gridpath1//'checkstate_'//trim(adjustl(checkstateDescriptor))//&
                      '_'//checkstateTrajectory//'.jpg"'
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
write(filechannel1,*) 'set label 1 "H2 Orientation: '//angle1descriptor//', '//angle2descriptor//&
                      ' radians" at screen 0.7, 0.975'
write(filechannel1,*) 'set label 2 "H2 Bond Length: '//bond1descriptor//' A" at screen 0.7, 0.955'
write(filechannel1,*) 'set ylabel "Var1 (A)"'
write(filechannel1,*) 'set yrange [0:40]'
write(filechannel1,*) 'plot "'//gridpath1//checkstatefile//'" u 4:6 w lines'
write(filechannel1,*) 'unset label 1'
write(filechannel1,*) 'unset label 2'
write(filechannel1,*) 'set ylabel "Var2 (A)"'
write(filechannel1,*) 'set yrange [0:40]'
write(filechannel1,*) 'plot "'//gridpath1//checkstatefile//'" u 4:7 w lines'
if (.false.) then
write(filechannel1,*) 'set ylabel "Neighbor Check?"'
write(filechannel1,*) 'unset ytics'
write(filechannel1,*) 'set yrange [-0.2:1.2]'
write(filechannel1,*) 'plot "'//gridpath1//checkstatefile//'" u 4:3 w lines, '//&
                           '"'//gridpath1//checkstatefile//'" u 4:($2==0?$3:4/0) w p ls 3, '//&
                           '"'//gridpath1//checkstatefile//'" u 4:(($3==1)&&($2<2)?$3:4/0) w p ls 1, '//&
                           '"'//gridpath1//checkstatefile//'" u 4:(($3==1)&&($2>1)?$3:4/0) w p ls 2'
write(filechannel1,*) 'set ylabel "Granularity (order)"'
write(filechannel1,*) 'set ytics 1'
write(filechannel1,*) 'set yrange [-0.5:3.5]'
write(filechannel1,*) 'plot "'//gridpath1//checkstatefile//'" u 4:2 w lines, '//&
                           '"'//gridpath1//checkstatefile//'" u 4:($2==0?$2:4/0) w p ls 3, '//&
                           '"'//gridpath1//checkstatefile//'" u 4:(($3==1)&&($2<2)?$2:4/0) w p ls 1, '//&
                           '"'//gridpath1//checkstatefile//'" u 4:(($3==1)&&($2>1)?$2:4/0) w p ls 2'
end if
write(filechannel1,*) 'set ylabel "Frames Checked"'
write(filechannel1,*) 'set ytics auto'
write(filechannel1,*) 'set yrange [0:300]'
write(filechannel1,*) 'plot "'//gridpath1//checkstatefile//'" u 4:1 w lines, '//&
                           '"'//gridpath1//checkstatefile//'" u 4:($2==0?$1:4/0) w p ls 3, '//&
                           '"'//gridpath1//checkstatefile//'" u 4:(($3==1)&&($2<2)?$1:4/0) w p ls 1, '//&
                           '"'//gridpath1//checkstatefile//'" u 4:(($3==1)&&($2>1)?$1:4/0) w p ls 2'
write(filechannel1,*) 'set xtics'
write(filechannel1,*) 'set ylabel "Minimum RMSD (A)"'
write(filechannel1,*) 'set xlabel "Timestep"'
write(filechannel1,*) 'set yrange [0:2.0e-1]'
write(filechannel1,*) 'unset key'
write(filechannel1,*) 'plot "'//gridpath1//checkstatefile//'" u 4:5 w lines'
close(filechannel1)
call system("gnuplot < "//gridpath1//temporaryfile2)
call system("rm "//gridpath1//temporaryfile2)

end do


end program checkMultipleTrajectories
