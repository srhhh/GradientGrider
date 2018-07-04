program analyzeMultipleTrajectories
use f2_physics_parameters
use f2_parameters
use checkCells5
use ls_rmsd
use mapCellData
use makeTrajectory5
use makeTrajectory
implicit none

!Variables used to name the new directory uniquely
character(resolution_text_length) :: resolution_text
character(overcrowd0_text_length) :: overcrowd0_text
character(trajectory_text_length) :: trajectory_text
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(gridpath_length) :: gridpath0
character(trajectories_text_length*100) :: trajectories_text
character(6) :: Ntraj_text

!New Trajectory Parameters
integer,parameter :: Ntesttraj = 100
real :: initial_bond_distance, initial_rotation_angle, initial_rotational_speed
real :: initial_bond_angle1, initial_bond_angle2
real :: initial_energy_H2,initial_vibrational_energy,initial_rotational_energy
real :: random_num1,random_num2,i,j
integer :: seed,n,m,n_testtraj,initial_n_testtraj

!Variables
integer :: Ngrid,iostate,Ngrid_total
integer,allocatable :: filechannels(:)
real,allocatable :: percent_threshold_rmsd(:)
real :: min_rmsd
integer :: frames, step, total_threshold_rmsd
real,parameter :: threshold_rmsd = 1.0e-5

!We identify which grid we are working with by what parameters were used to build it
!In this case, the grid is uniquely determined by resolution0, overcrowd0, and overcrowd1
write(variable_length_text,FMT="(I5)") resolution_text_length
write(resolution_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") resolution_0
write(variable_length_text,FMT="(I5)") overcrowd0_text_length
write(overcrowd0_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") overcrowd0
write(variable_length_text,FMT="(I5)") trajectory_text_length
write(trajectory_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ntraj_max
gridpath0 = path5//resolution_text//"_"//overcrowd0_text//"_"//trajectory_text//"/"

!All trajectory folders are formatted as I0.3 (3-digit integer)
!So search for these numbered folders and read them
call system("ls -p "//gridpath0//" | grep '[0123456789]/' > "//gridpath0//trajectories)

!This is a slight 'hack'; for the 'true' scattering angle plots based on
!the trajectories FROM THE GRID, we only need to concatenate
!all the trajectories and read off a specific column
!To do this, we make a long string " file1 file2 file3 ... fileN"
!And only take whatever portion out we want (ex. " file1 file2")
open(trajectorieschannel,file=gridpath0//trajectories,action="read")
trajectories_text = ""
do Ngrid_total = 1, Ngrid_max
	read(trajectorieschannel,FMT="(A4)",iostat=iostate) folder_text
	if (iostate /= 0) exit
	trajectories_text = trim(trajectories_text)//" "//gridpath0//folder_text//&
			    trajectoriesfile
end do
close(trajectorieschannel)

!This subtraction by two is for two reasons:
!  1. We increment before we read, so we need to subtract out one increment here
!  2. The last folder may have a problem (it does in this case) so we subtract one out here
if (Ngrid_total < 3) return
Ngrid_total = Ngrid_total - 2

!We print the seed in case there's some bug that we need to reproduce
call system_clock(seed)
print *, "Working with system_clock seed: ", seed
print *, ""
seed = rand(seed)

!This part does not take long
print *, "Finished Initialization Part 1..."

!For now, we are putting a soft cap on how many trajectories we are using
!Each grid has Ntraj_max trajectories so we will use a maximum of 10 * Ntraj_max trajectories
Ngrid_total = min(10, Ngrid_total)
allocate(filechannels(Ngrid_total))

!First, we do the 'true' scattering angle plots
!This data was made during creation (or should have been!) so all we need to do
!is merge, read, and plot them
write(variable_length_text,FMT="(I5)") Ngrid_text_length
do Ngrid = 1, Ngrid_total

	!The folders are named starting from 001 by increments of 1
	write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid

	!The plots are named starting from Ntraj_max by increments of Ntraj_max (the number of trajectories)
	write(Ntraj_text,FMT="(I0.6)") Ngrid * Ntraj_max
	print *, " Working with this many trajectories: ", Ntraj_text

	!This system call concatenates all the data files from that previously made 'hack'
	!By doing that, we merge all the scattering angle data together
	call system("cat"//trajectories_text(1:Ngrid*trajectories_text_length)//" >> "//&
		    gridpath0//Ngrid_text//"/"//cumulativefile//Ntraj_text//".dat")

	!This is the gnuplot code to make the plots
	open(gnuplotchannel,file=gridpath0//gnuplotfile)
	write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
        write(gnuplotchannel,*) 'set output "'//gridpath0//'ScatteringAngle'//Ntraj_text//'"'
        write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
        write(gnuplotchannel,*) 'unset key'
        write(gnuplotchannel,*) 'bin_width = 0.001'
        write(gnuplotchannel,*) 'bin_number(x) = floor(x/bin_width)'
        write(gnuplotchannel,*) 'rounded(x) = bin_width * (bin_number(x) + 0.5)'
        write(gnuplotchannel,*) 'set xlabel "Scattering Angle"'
        write(gnuplotchannel,*) 'set ylabel "Occurence"'
        write(gnuplotchannel,*) 'plot "'//gridpath0//Ngrid_text//'/'//cumulativefile//Ntraj_text//&
				'.dat" u (rounded($7)):(7) smooth frequency with boxes'
	close(gnuplotchannel)

	!And then we just input it into gnuplot.exe
	call system("gnuplot < "//gridpath0//gnuplotfile)

end do

!This part does not take long
print *, "Finished Initialization Part 2..."

!I want this program to be a "pick up where we left off last time" program
!so I figure out how many new trajectories I already checked for RMSD
!All grids should have the same number of trajectories so I just look in grid 1
call system("ls "//gridpath0//"001/ | grep -E '^[0123456789]{6}.dat' > "//gridpath0//trajectories)
open(trajectorieschannel,file=gridpath0//trajectories,action="read")
initial_n_testtraj = 1
do Ngrid_total = 1, Ngrid_max
	read(trajectorieschannel,FMT="(A50)",iostat=iostate) trajectories_text
	if (iostate /= 0) exit
	initial_n_testtraj = initial_n_testtraj + 1
	print *, "     Already have trajectory: ", trim(adjustl(trajectories_text))
end do
close(trajectorieschannel)

!Now here we actually make and check these new trajectories
do n_testtraj = initial_n_testtraj, Ntesttraj

	!This is just the creation of the random initial trajectory
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

	!Each trajectory will have Ngrid_total outputs; one for however many grids we use
	!The trajectory number will uniquely identify one trajectory from another
	write(Ntraj_text,FMT="(I0.6)") n_testtraj
	print *, " Working on random new trajectory number: ", Ntraj_text

	!The grid number will uniquely identify one trajectory
	!Open all these files under filechannels
	do Ngrid = 1, Ngrid_total
		write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
		filechannels(Ngrid) = 1000 + 69 * Ngrid
		open(filechannels(Ngrid),file=gridpath0//Ngrid_text//"/"//Ntraj_text//".dat")
	end do

	!The write the outputted RMSDS of each trajectory onto those filechannels
	!Remark: checkCells5 uses filechannel1 to open files in the grid
	call checkMultipleTrajectories(initial_bond_distance,&
		    initial_rotational_speed,initial_rotation_angle,&
                    initial_bond_angle1,initial_bond_angle2,.false.,&
                    Ngrid_total,filechannels(1:Ngrid_max),gridpath0)

	!Finally, close them
	do Ngrid = 1, Ngrid_total
		close(filechannels(Ngrid))
	end do
end do

!This part has not yet finished yet
print *, "Finished Initialization Part 3..."
print *, ""

!Now, all we need to do is read the RMSDs obtained from each
!With some processing to get valuable data
do Ngrid = 1, Ngrid_total
write(Ngrid_text,FMT="(I0."//variable_length_text//")") Ngrid
print *, " Working on all grids up until grid number: ", Ngrid_text

!We will bin data by GRID, not by trajectory
!So we uniquely name each output .dat and graph by the grid number
open(filechannel1,file=gridpath0//"percent_rmsd"//Ngrid_text//".dat")
do n_testtraj = 1, Ntesttraj
	write(Ntraj_text,FMT="(I0.6)") n_testtraj

	!Read the trajectory (which has the rmsd) across all grids line-by-line
	iostate = 0
	frames = 0
	total_threshold_rmsd = 0
	open(filechannel2,file=gridpath0//Ngrid_text//"/"//Ntraj_text//".dat")

	do
		read(filechannel2,FMT=FMT6,iostat=iostate) min_rmsd
		if (iostate /= 0) exit
		frames = frames + 1

	!If the RMSD is below the threshhold we consider tally that
		if (min_rmsd < threshold_rmsd) total_threshold_rmsd = total_threshold_rmsd + 1
	end do
	close(filechannel2)

	!We want the percentage of frames that has an RMSD below the threshhold
	!So we keep track of the number of frames and divide by that
	percent_threshold_rmsd(n_testtraj) = total_threshold_rmsd * 1.0 / frames
	write(filechannel1,FMT="(I6,1x,F7.4,1x,I8)") n_testtraj, percent_threshold_rmsd(n_testtraj), frames
end do
close(filechannel1)

!Finally, plot the data
open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath0//'PercentRMSDThreshold'//Ngrid_text//'"'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'bin_width = 0.001'
write(gnuplotchannel,*) 'bin_number(x) = floor(x/bin_width)'
write(gnuplotchannel,*) 'rounded(x) = bin_width * (bin_number(x) + 0.5)'
write(gnuplotchannel,*) 'set xlabel "Scattering Angle"'
write(gnuplotchannel,*) 'set ylabel "Occurence"'
write(gnuplotchannel,*) 'plot "'//gridpath0//'percent_rmsd'//Ngrid_text//'.dat'//&
                        '" u (rounded($7)):(7) smooth frequency with boxes'
close(gnuplotchannel)

call system("gnuplot < "//gridpath0//gnuplotfile)

end do

end program analyzeMultipleTrajectories
