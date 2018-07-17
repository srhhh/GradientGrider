program checkNewTrajectorieswithMultipleGrids
use ANALYSIS
use PARAMETERS
use checkMultipleGrids
use mapCellData
use checkNewTrajectorywithMultipleGrids
use addNewTrajectorytoGrid
use analyzeHeatMapswithMultipleGrids
use analyzeRMSDThresholdwithMultipleGrids
use analyzeScatteringAngleswithMultipleGrids
implicit none

!Variables used to name the new directory uniquely
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectories_text_length*100) :: trajectories_text
character(6) :: Ntraj_text
character(6) :: Nthreshold_text
character(6) :: reject_text

!New Trajectory Parameters
real(dp) :: initial_bond_distance, initial_rotation_angle, initial_rotational_speed
real(dp) :: initial_bond_angle1, initial_bond_angle2
real(dp) :: initial_energy_H2,initial_vibrational_energy,initial_rotational_energy
real(dp) :: random_num1,random_num2,random_num3,random_r2,random_r3,i,j
real(dp) :: speedH, speedH2, scattering_angle
real(dp),dimension(3) :: velocityH,velocityH2
integer :: seed,n,m,n_testtraj,initial_n_testtraj

!Variables
integer :: Ngrid,iostate,Ngrid_total
integer,allocatable :: filechannels(:)
integer,dimension(counter0_max) :: counter0
real(dp),allocatable :: percent_threshold_rmsd(:)
real(dp) :: min_rmsd
integer :: frames, step, total_threshold_rmsd

!Timing
real :: r1, r2, system_clock_rate
integer :: c1, c2, cr

!All trajectory folders are formatted as I0.3 (3-digit integer)
!So search for these numbered folders and read them
call system("ls -p "//gridpath0//" | grep '[0123456789]/' > "//gridpath0//trajectories)

!Let's see how many grids are actually in the folder
open(trajectorieschannel,file=gridpath0//trajectories,action="read")
do Ngrid_total = 1, Ngrid_max
	read(trajectorieschannel,FMT="(A4)",iostat=iostate) folder_text
	if (iostate /= 0) exit
end do
close(trajectorieschannel)

!We increment before we read, so we need to subtract out one increment here
if (Ngrid_total < 2) return
Ngrid_total = Ngrid_total - 1

!Keep the cap to the number of grids in mind
Ngrid_total = min(Ngrid_cap, Ngrid_total)

print *, ""
print *, "Working on directory ", gridpath0
print *, "Deciding on using ", Ngrid_total, " grids"
print *, ""

!We may need a filechannel open for each grid
allocate(filechannels(Ngrid_total))

!And we print the seed in case there's some bug that we need to reproduce
call system_clock(seed)
print *, "Working with system_clock seed: ", seed
print *, ""
seed = rand(seed)

!This is for top-level heat map generation (from the grid)
if (heatmap_flag) call analyzeHeatMaps1(Ngrid_total)

!This is for scattering angle plots (from the grid)
if (trueSA_flag) call getScatteringAngles1(Ngrid_total,trajectoriesfile,8,"trueSA.jpg")





!Some intitialization stuff
if (reject_flag) then
	reject_text = "reject"
else
	reject_text = "accept"
end if

write(variable_length_text,FMT="(I5)") Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid_total
write(Ntraj_text,FMT="(F6.5)") threshold_rmsd
call system("rm "//gridpath0//Ngrid_text//reject_text//Ntraj_text//trajectoriesfile)


!This is for checking trajectories against the grid
!Currently, this is the main use of this program
if (testtraj_flag) then

!If we want this program to be a "pick up where we left off last time" program
!we figure out how many new trajectories I already checked for RMSD
!We need all grids 1 to Ngrid_total to have the trajectories so we just check the last one
initial_n_testtraj = 1
if (useolddata_flag) then
print *, "     Deciding to use old data..."
write(variable_length_text,FMT="(I5)") Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid_total
call system("ls "//gridpath0//Ngrid_text//"/ | grep -E '^[0123456789]{6}.dat' > "//gridpath0//trajectories)
open(trajectorieschannel,file=gridpath0//trajectories,action="read")
do
	read(trajectorieschannel,FMT="(A50)",iostat=iostate) trajectories_text
	if (iostate /= 0) exit
	initial_n_testtraj = initial_n_testtraj + 1
	print *, "     Already have trajectory: ", trim(adjustl(trajectories_text))
end do
print *, ""
close(trajectorieschannel)
end if

!Start timing because this part of the program takes a lot of time
!so we want to measure our progress
call system_clock(c1,count_rate=cr)
system_clock_rate = 1.0/real(cr)


write(Nthreshold_text,FMT="(F6.5)") threshold_rmsd
!Now here we actually make these new trajectories
do n_testtraj = initial_n_testtraj, Ntesttraj

	!This is just the creation of the random initial trajectory
	!The orientation of the H2
	do
	random_num1 = rand() - 0.5d0
	random_num2 = rand() - 0.5d0
	random_num3 = rand() - 0.5d0
	random_r2 = random_num1**2 + random_num2**2
	random_r3 = random_r2 + random_num3**2
	if (random_r3 > 0.25d0) cycle
	random_r2 = sqrt(random_r2)
	initial_bond_angle1 = acos(random_num1 / random_r2)
	initial_bond_angle2 = atan2(random_r2,random_num3)
	exit
	end do
	!The energy of the H2
        do
        random_num1 = rand()
        random_num2 = rand()
        initial_energy_H2 = (upsilon_max*random_num1 + 0.5d0)*upsilon_factor2
        if (random_num2 < temperature_scaling*exp(upsilon_max*random_num1*upsilon_factor1)) exit
        end do
        random_num2 = 1.0d0
        initial_vibrational_energy = random_num2*initial_energy_H2
        initial_rotational_energy = initial_energy_H2 - initial_vibrational_energy
        initial_bond_distance = HOr0_hydrogen + sqrt(initial_vibrational_energy*2/HOke_hydrogen)
        random_num1 = rand()
        initial_rotational_speed = sqrt(initial_rotational_energy/mass_hydrogen)
        initial_rotation_angle = random_num1*pi2

	!Each trajectory will have Ngrid_total outputs; one for however many grids we use
	!The trajectory number will uniquely identify one trajectory from another
	write(Ntraj_text,FMT="(I0.6)") n_testtraj
	
	call system_clock(c1)
	call CPU_time(r1)
	print *, " Working on random new trajectory number: ", Ntraj_text

	!The grid number will uniquely identify one trajectory
	!Open all these files under filechannels
	write(variable_length_text,FMT="(I5)") Ngrid_text_length
	do Ngrid = 1, Ngrid_total
		write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
		filechannels(Ngrid) = 1000 + 69 * Ngrid
		open(filechannels(Ngrid),file=gridpath0//Ngrid_text//"/"//&
                                              reject_text//Nthreshold_text//'_'//Ntraj_text//".dat")
	end do

	!Then write the outputted RMSDS of each trajectory onto those filechannels
	!Remark: checkMultipleGrids uses filechannel1 to open files in the grid
	call checkMultipleTrajectories(initial_bond_distance,&
		    initial_rotational_speed,initial_rotation_angle,&
                    initial_bond_angle1,initial_bond_angle2,.false.,threshold_rmsd,reject_flag,&
                    Ngrid_total,filechannels(1:Ngrid_max),gridpath0,velocityH,velocityH2)

	!Also let's see how long a single trajectory takes
	call system_clock(c2)
	call CPU_time(r2)
	print *, "        CPU Time: ", r2 - r1
	print *, "       Wall Time: ", (c2 - c1) * system_clock_rate

	!Finally, close them
	do Ngrid = 1, Ngrid_total
		close(filechannels(Ngrid))
	end do

        !We also record the scattering angle of the trajectory
        speedH = sqrt(velocityH(1)**2 + velocityH(2)**2 + velocityH(3)**2)
        speedH2 = sqrt(velocityH2(1)**2 + velocityH2(2)**2 + velocityH2(3)**2)
        scattering_angle = acos(dot_product(velocityH,velocityH2) / &
                                           (speedH * speedH2))

        !This is all recorded in the trajectoriesfile of the grid
        open(filechannel1,file=gridpath0//Ngrid_text//reject_text//Nthreshold_text//trajectoriesfile,position="append")
        write(filechannel1,*) r2-r1,(c2-c1)*system_clock_rate,scattering_angle, initial_bond_angle1, initial_bond_angle2
        close(filechannel1)


	if (testtrajRMSD_flag) then
		open(gnuplotchannel,file=gridpath0//gnuplotfile)
		write(gnuplotchannel,*) 'set term jpeg size 600, 600'
		write(gnuplotchannel,*) 'set output "'//gridpath0//Ngrid_text//'/RMSD'//&
                                                        reject_text//Ntraj_text//'.jpg"'
		write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
		write(gnuplotchannel,*) 'unset key'
		write(gnuplotchannel,*) 'bin_width = 0.001'
		write(gnuplotchannel,*) 'bin_number(x) = floor(x/bin_width)'
		write(gnuplotchannel,*) 'rounded(x) = bin_width * (bin_number(x) + 0.5)'
		write(gnuplotchannel,*) 'set xlabel "RMSD"'
		write(gnuplotchannel,*) 'set ylabel "Occurence"'
		write(gnuplotchannel,*) 'plot "'//gridpath0//Ngrid_text//"/"//reject_text//&
                                        Nthreshold_text//'_'//Ntraj_text//'.dat'//&
		                        '" u (rounded($1)):(1.0) smooth frequency with boxes'
		close(gnuplotchannel)
		call system("gnuplot < "//gridpath0//gnuplotfile)
	end if
	
end do
print *, ""

end if


write(variable_length_text,FMT="(I5)") Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid_total
write(Ntraj_text,FMT="(F6.5)") threshold_rmsd
if (percentthreshold_flag) call getRMSDThresholds1(Ngrid_total,1,&
                                                   "PercentRMSDThreshold_"//&
                                                   Ngrid_text//reject_text//Ntraj_text)
print *, "   Making plot: ", "PercentRMSDThreshold_"//Ngrid_text//reject_text//Ntraj_text
print *, ""

if (testtrajSA_flag) call getScatteringAngles2(Ngrid_text//reject_text//Ntraj_text//trajectoriesfile,&
                                               Ntesttraj,3,4,5,"TestScatteringAngleDistribution_"//&
                                               Ngrid_text//reject_text//Ntraj_text)
print *, "   Making plot: ", "TestScatteringAngleDistribution_"//Ngrid_text//reject_text//Ntraj_text
print *, ""

end program checkNewTrajectorieswithMultipleGrids
