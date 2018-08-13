!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PROGRAM
!               checkNewTrajectorieswithMultipleGrids
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!		This program performs an analysis of a library of grids
!		specified in PARAMETERS, with the analysis guided by
!		logical flags in ANALYSIS
!
!		Often, new trajectories will have to be made with certain
!		rejection method and RMSD threshold; these sample trajectories
!		have each frame checked with the grids, and an approximate
!		gradient used if within threshold; these frames are not added
!		to the grids.
!
!		The graphs that output include: heat maps, scattering angle
!		plots, and percent rmsd threshold plots
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILECHANNELS                    ACTION
!
!               GNUPLOTCHANNEL                  OPEN, WRITE, CLOSE
!               FILECHANNEL1                    OPEN, WRITE, CLOSE
!		TRAJECTORIESCHANNEL		OPEN, WRITE, CLOSE
!		FILECHANNELS(GRID#)		OPEN, RUNTRAJECTORY, CLOSE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       CALLS                           MODULE
!
!               SYSTEM                          INTRINSIC
!               CPU_TIME                        INTRINSIC
!               SYSTEMCLOCK                     INTRINSIC
!
!               checkMultipleTrajectories       runTrajectory
!
!               getScatteringAngles1            analyzeScatteringAngleswithMultipleGrids
!               getScatteringAngles2            analyzeScatteringAngleswithMultipleGrids
!		analyzeHeatMaps1		analyzeHeatMapswithMultipleGrids
!		getRMSDThresholds1		analyzeRMSDThresholdwithMultipleGrids
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!		gridpath0//trajectories		DAT			A list of grids (001/, 002/, ...)
!		gridpath1//reject//		DAT			The min_rmsd retrieved for each step
!			threshold//TRAJ#				of a trajectory
!		gridpath0//Ngrid///reject//	DAT			A list of data from all trajectories with the
!		     threshold//trajectories				same suffix (Ngrid,reject,threshold)
!		gridpath1//HeatMap		JPG			The heat map of a grid
!		gridpath0//TrueSADist//		JPG			The scattering angle distribution of
!			TOTALTRAJ#					TOTALTRAJ# trajectories (multiple grids)
!		gridpath0//TestSADist//		JPG			The scattering angle distribution of all
!		     Ngrid//reject//threshold				trajectoreis with the same suffix
!		gridpath1//PercentRMSDDist//	JPG			The percent RMSD threshold distribution of
!		     Ngrid//reject//threshold				all trajectories with the same suffix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program checkNewTrajectorieswithMultipleGrids
use ANALYSIS
use PARAMETERS
use interactMultipleGrids
use runTrajectory
use analyzeHeatMapswithMultipleGrids
use analyzeRMSDThresholdwithMultipleGrids
use analyzeScatteringAngleswithMultipleGrids
implicit none

!Variables used to name the new directory and trajectories uniquely
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(6) :: Ntraj_text
character(6) :: Nthreshold_text
character(6) :: reject_text
character(12) :: prefix_text

!New Trajectory Parameters
real(dp) :: initial_bond_distance, initial_rotation_angle, initial_rotational_speed
real(dp) :: initial_bond_angle1, initial_bond_angle2
real(dp) :: initial_energy_H2,initial_vibrational_energy,initial_rotational_energy
real(dp) :: random_num1,random_num2,random_num3,random_r2,random_r3
real(dp) :: scattering_angle
real(dp),dimension(3) :: TRVenergies1,TRVenergies2,dTRVenergies
real(dp),dimension(3,3) :: coords_initial,velocities_initial,coords_final,velocities_final
real(dp) :: totalEnergy
integer :: seed,n,m,n_testtraj,initial_n_testtraj

!Variables
integer :: iostate
integer,allocatable :: filechannels(:)

!Timing
real :: r1, r2, system_clock_rate
integer :: c1, c2, cr

!Incremental Integers
integer :: i, j, k


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       RNG SETUP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!And we print the seed in case there's some bug that we need to reproduce
call system_clock(seed)
print *, ""
print *, "   RNG seed: ", seed
print *, ""
seed = rand(seed)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       PRE-PROCESSING LIBRARY OVERVIEW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!All trajectory folders are formatted as I0.3 (3-digit integer)
!So search for these numbered folders and read them
call system("ls -p "//gridpath0//" | grep '[0123456789]/' > "//gridpath0//trajectories)

!Let's see how many grids are actually in the folder
Ngrid_max = 0
open(trajectorieschannel,file=gridpath0//trajectories,action="read")
do
	read(trajectorieschannel,FMT="(A4)",iostat=iostate) folder_text
	if (iostate /= 0) exit
	Ngrid_max = Ngrid_max + 1
end do
close(trajectorieschannel)

!Keep the cap to the number of grids in mind
Ngrid_total = min(Ngrid_cap, Ngrid_max)

print *, ""
print *, "Analysis on directory ", gridpath0
print *, "Deciding on using ", Ngrid_total, " grids"
print *, ""


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       PRE-CREATION LIBRARY ANALYSIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!This is for top-level heat map generation (for each of the grids)
if (heatmap_flag) then
	print *, "   Making plot: ", "TopLevel_HeatMap"
	print *, ""
	call analyzeHeatMaps2()
end if

!This is for scattering angle plots (from each of the grids)
if (trueSA_flag) then
	print *, "   Making plot: ", "trueScatteringAngleDistribution"
	print *, ""
	call getScatteringAngles1("Initial","trueScatteringAngleDistribution.jpg")
end if

!This is for the energy decomposition plots
if (trueED_flag) then
	print *, "   Making plot: ", "trueEnergyDecompositionDistribution"
	print *, ""
	call getTRVimages("Initial","trueEnergyDecompositionDistribution.jpg")
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       TRAJECTORY CREATION FLAG START
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!This is for checking trajectories against the grid
!Currently, this is the main use of this program
if (testtraj_flag) then

print *, ""
print *, "Trajectory Creation ... Start"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Some intitialization stuff
write(Nthreshold_text,FMT=FMT6_pos_real0) threshold_rmsd
if (reject_flag) then
	reject_text = "reject"
else
	reject_text = "accept"
end if
allocate(filechannels(Ngrid_total))
write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid_total

prefix_text = reject_text//Nthreshold_text

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       PRE-CREATION LIBRARY OVERVIEW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!If we want this program to be a "pick up where we left off last time" program
!we figure out how many new trajectories I already checked for RMSD

!By default, we say we will remake all the trajectories
initial_n_testtraj = 1

!And then we check if we want to reuse some old ones
if (useolddata_flag) then

        print *, "     Deciding to use old data..."

        !First we check how many trajectory (Ntraj_text format) files we have
        !This system call puts them all onto one file
	!We need all grids 1 to Ngrid_total to have the trajectories but we just check the last one
	!and we assume all previous grids have just as many or more (usually a safe assumption)
        write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid_total
        call system("ls "//gridpath0//Ngrid_text//"/ | grep -E '^"//&
                    reject_text//"\"//Nthreshold_text//"_[0123456789]{6}.dat' > "//gridpath0//trajectories)

        !Then we simply have to read the file to see how many we have
	!This also assumes they were numbered correctly (usually a safe assumption)
        open(trajectorieschannel,file=gridpath0//trajectories,action="read")
        do
                read(trajectorieschannel,*,iostat=iostate)
                if (iostate /= 0) exit
                initial_n_testtraj = initial_n_testtraj + 1
        end do
        close(trajectorieschannel)

        !Second, we check how many lines are on the trajectories file
        Ntraj = 1
        open(trajectorieschannel,file=gridpath0//Ngrid_text//prefix_text//trajectoriesfile)
        do
                read(trajectorieschannel,*,iostat=iostate)
                if (iostate /= 0) exit
                Ntraj = Ntraj + 1
        end do
        close(trajectorieschannel)

        !If these two values are not congruent there has been some sort of file corruption or deletion
        if (Ntraj /= initial_n_testtraj) then
                print *, ""
                print *, ""
                print *, ""
                print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                print *, "          DATA INCONGRUENCE "
                print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                if (Ntraj > initial_n_testtraj) print *, "    trajectory files missing for "//&
                                                         Ngrid_text//prefix_text//" ... .dat"
                if (Ntraj < initial_n_testtraj) print *, "    trajectories file corrupted for "//&
                                                         Ngrid_text//prefix_text
                print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                print *, ""
                print *, "    Program exited disgracefully"
                print *, ""
                return
        end if

        print *, "     Total Number of Trajectories Saved: ", initial_n_testtraj - 1, " out of ", Ntesttraj
        print *, ""

!If not...
else

        !We need to make a new one of these trajectory files altogether
        call system("rm "//gridpath0//Ngrid_text//prefix_text//trajectoriesfile)
end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       TRAJECTORY CREATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!Start timing because this part of the program takes a lot of time
!so we want to measure our progress
call system_clock(c1,count_rate=cr)
system_clock_rate = 1.0/real(cr)

!Now here we actually make these new trajectories
do n_testtraj = initial_n_testtraj, Ntesttraj

	!This is just the creation of the random initial trajectory
	do n = 1, Nbonds
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
                        !This picks a random value between zero and some very high upper limit
                        random_num1 = rand() * upsilon_max
                        random_num2 = rand() * upsilon_factor2
                        if (exp(-random_num1 * upsilon_factor1) < random_num2) cycle
                        initial_energy_H2 = (random_num1 + 0.5d0) * epsilon_factor
			exit
                end do
	        random_num2 = 1.0d0
	        initial_vibrational_energy = random_num2*initial_energy_H2
	        initial_rotational_energy = initial_energy_H2 - initial_vibrational_energy
	        initial_bond_distance = HOr0_hydrogen + sqrt(initial_vibrational_energy*2/HOke_hydrogen)
	        random_num1 = rand()
	        initial_rotational_speed = sqrt(initial_rotational_energy/mass_hydrogen)
	        initial_rotation_angle = random_num1*pi2

		totalEnergy = initial_translational_KE + initial_vibrational_energy + &
                              initial_rotational_energy

		INITIAL_BOND_DATA(n,:) = (/ initial_bond_distance,initial_rotational_speed,&
                           initial_rotation_angle,initial_bond_angle1,initial_bond_angle2 /)

	end do

	open(filechannel1,file=gridpath0//Ngrid_text//prefix_text//initialfile,&
                          position="append")
	write(filechannel1,FMTinitial) ((INITIAL_BOND_DATA(i,j),j=1,5),i=1,Nbonds)
        close(filechannel1)


	!Each trajectory will have Ngrid_total outputs; one for however many grids we use
	!The trajectory number will uniquely identify one trajectory from another
	write(Ntraj_text,FMT=FMT6_pos_int) n_testtraj
	
	call system_clock(c1)
	call CPU_time(r1)
	print *, " Working on random new trajectory number: ", Ntraj_text

	!The grid number will uniquely identify one trajectory
	!Open all these files under filechannels
	write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
	do Ngrid = 1, Ngrid_total
		write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
		filechannels(Ngrid) = 1000 + 69 * Ngrid
		open(filechannels(Ngrid),file=gridpath0//Ngrid_text//"/"//&
                                              prefix_text//'_'//Ntraj_text//".dat")
	end do

	!Then write the outputted RMSDS of each trajectory onto those filechannels
	!Remark: checkMultipleGrids uses filechannel1 to open files in the grid
	call checkMultipleTrajectories(filechannels(1:Ngrid_max),coords_initial,velocities_initial,&
                                                                 coords_final,velocities_final)

	!Also let's see how long a single trajectory takes
	call system_clock(c2)
	call CPU_time(r2)
	print *, "        CPU Time: ", r2 - r1
	print *, "       Wall Time: ", (c2 - c1) * system_clock_rate

	!Finally, close them
	do Ngrid = 1, Ngrid_total
		close(filechannels(Ngrid))
	end do

	!The only data that is recorded is the first and last frame of the trajectory
	open(filechannel1,file=gridpath0//Ngrid_text//prefix_text//timeslicefile,position="append")
	write(filechannel1,FMTtimeslice) &
                              ((coords_initial(i,j),i=1,3),j=1,Natoms),&
                              ((velocities_initial(i,j),i=1,3),j=1,Natoms),&
			      ((coords_final(i,j),i=1,3),j=1,Natoms),&
                              ((velocities_final(i,j),i=1,3),j=1,Natoms)
        close(filechannel1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       DURING-CREATION TRAJECTORY ANALYSIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	if (testtrajRMSD_flag) then
		open(gnuplotchannel,file=gridpath0//gnuplotfile)
		write(gnuplotchannel,*) 'set term pngcairo size 600, 600'
		write(gnuplotchannel,*) 'set output "'//gridpath0//Ngrid_text//'/RMSD'//&
                                                        reject_text//Ntraj_text//'.png"'
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
		call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)
	end if
	
end do
print *, ""


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       POST-CREATION LIBRARY ANALYSIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Finally, do a post-creation timeslice-to-SA conversions here
!We use the SA often so we do this at the beginning
write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid_total
call postProcess(Ngrid_text//prefix_text)

if (percentthreshold_flag) then
	print *, "   Making plot: ", "PercentRMSDThreshold_"//Ngrid_text//prefix_text
	print *, ""
	call getRMSDThresholds1(prefix_text,"PercentRMSDThreshold_"//Ngrid_text//prefix_text)
end if

Ntraj = Ntesttraj
if (testtrajSA_flag) then
	print *, "   Making plot: ", "TestScatteringAngleDistribution_"//Ngrid_text//prefix_text
	print *, ""
	call getScatteringAngles2(Ngrid_text//prefix_text,"TestScatteringAngleDistribution")
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       TRAJECTORY CREATION FLAG END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print *, ""
print *, "Successfully exited analysis"
print *, ""

end program checkNewTrajectorieswithMultipleGrids
