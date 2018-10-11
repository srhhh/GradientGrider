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
character(trajectory_text_length) :: Ntraj_text
character(6) :: Nthreshold_text
character(6) :: reject_text
character(12) :: prefix_text
character(150) :: old_filename, new_filename
character(1) :: answer

!New Trajectory Parameters
real(dp) :: initial_bond_distance, initial_rotation_angle, initial_rotational_speed
real(dp) :: initial_bond_angle1, initial_bond_angle2
real(dp) :: initial_energy_H2,initial_vibrational_energy,initial_rotational_energy
real(dp) :: bond_period_elapsed
real(dp) :: random_num1,random_num2,random_num3,random_r2,random_r3
real(dp) :: scattering_angle
real(dp),dimension(3) :: TRVenergies1,TRVenergies2,dTRVenergies
real(dp),dimension(3,Natoms) :: coords_initial,velocities_initial,coords_final,velocities_final
real(dp) :: probJ_max, J_factor3
integer :: seed,n,m,n_testtraj,initial_n_testtraj

!Variables
integer :: iostate
integer,allocatable :: filechannels(:)

!Timing
real :: r1, r2, system_clock_rate
integer :: c1, c2, cr
integer,dimension(3) :: now
real :: grid_wall_time,checktrajectory_wall_time
character(10) :: grid_wall_time_text, checktrajectory_wall_time_text
integer :: trajectory_t0, trajectory_t1

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
call itime(now)
write(6,FMT=FMTnow) now
print *, "   RNG seed: ", seed
print *, ""
seed = rand(seed)

write(Nbond_text,FMT="(I0.6)") Nbonds
write(Natom_text,FMT="(I0.6)") Natoms
write(FMTinitial,FMT="(A17)") "("//Nbond_text//"(6(F8.4)))"
write(FMTtimeslice,FMT="(A19)") "("//Natom_text//"(12(F12.7)))"
write(FMT2,FMT="(A22)") "("//Natom_text//"(6(1x,F14.10)))"
write(FMT3,FMT="(A22)") "("//Natom_text//"(3(1x,F14.10)))"

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
call itime(now)
write(6,FMT=FMTnow) now
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
	call itime(now)
	write(6,FMT=FMTnow) now
	print *, "   Making plot: ", "TopLevel_HeatMap"
	print *, ""
	call analyzeHeatMaps2()
end if

!This is for scattering angle plots (from each of the grids)
if (trueSA_flag) then
	call itime(now)
	write(6,FMT=FMTnow) now
	print *, "   Making plot: ", "InitialSATRVDistribution"
	print *, ""
	
	max_TranslationalEnergy = 0.0d0
        max_absenergychange = 0.0
        min_absenergychange = 1.0e9
        max_relenergychange = 0.0
        min_relenergychange = 1.0e9
	max_rotenergychange = 0.0
	min_rotenergychange = 1.0e9
	
	old_filename = ""
	do Ngrid = 1, Ngrid_total
		write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
		write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid

		!First, post process
	        Ntraj = Ntraj_max
        	call postProcess(Ngrid_text//"/Initial")

	        Ntraj = Ngrid*Ntraj_max
		write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
		write(Ntraj_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ntraj

	        if (trim(adjustl(old_filename)) /= "") then
	                new_filename = gridpath0//Ngrid_text//"/Initial"//Ntraj_text//SATRVfile
	                call system("cat "//trim(adjustl(old_filename))//" "//&
	                            gridpath0//Ngrid_text//"/Initial"//SATRVfile//" > "//new_filename)
	                old_filename = new_filename
	        else
	                old_filename = gridpath0//Ngrid_text//"/Initial"//Ntraj_text//SATRVfile
	                call system("cp "//gridpath0//Ngrid_text//"/Initial"//SATRVfile//" "//trim(adjustl(old_filename)))
	        end if

	        !Make an SA + TRV plot
	        call getScatteringAngles1(Ngrid_text//"/Initial"//Ntraj_text,"SATRVDistribution")
	end do
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
call itime(now)
write(6,FMT=FMTnow) now
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
        open(trajectorieschannel,file=gridpath0//Ngrid_text//prefix_text//timeslicefile)
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
                if (Ntraj < initial_n_testtraj) print *, "    timeslice file corrupted for "//&
                                                         Ngrid_text//prefix_text
                print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		print *, ""
		
		if ((Ntraj < Ntesttraj).or.(initial_n_testtraj < Ntesttraj)) then
		do
			print *, "    corrupted data interrupts trajectory indexing;"
			print *, "    overwrite corrupted data?    (y/n)"
			print *, ""
			read (*,*) answer
			if ((answer == "y").or.(answer == "n")) exit
		end do
		if (answer == "n") then
	                print *, ""
			call itime(now)
			write(6,FMT=FMTnow) now
	                print *, "    Program exited disgracefully"
	                print *, ""
	                return
		else
			initial_n_testtraj = Ntraj
		end if
		else
			print *, "    going along with analysis on uncorrupted data"
		end if
                print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		print *, ""
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
			initial_bond_angle1 = atan2(random_num1, random_r2)
			initial_bond_angle2 = atan2(random_r2,random_num3)
			exit
		end do

                do
                        random_num1 = rand() * upsilon_max
                        random_num2 = rand()
                        if (exp(-random_num1 * upsilon_factor1) * upsilon_factor2 < random_num2) cycle
                        initial_vibrational_energy = (random_num1 + 0.5d0) * epsilon_factor
                        exit
                end do
                initial_bond_distance = HOr0_hydrogen + sqrt(initial_vibrational_energy*2/HOke_hydrogen)
                J_factor3 = J_factor1 / (initial_bond_distance**2)
                probJ_max = sqrt(2*J_factor3) * exp(J_factor3*0.25d0 - 0.5d0)

                do
                        random_num1 = rand() * J_max
                        random_num2 = rand() * probJ_max
                        if ((2*random_num1 + 1.0d0) * J_factor3 * exp(-random_num1 * (random_num1 + 1.0d0) * &
                            J_factor3) < random_num2) cycle
                        initial_rotational_energy = (random_num1) * (random_num1 + 1.0d0) * J_factor2
                        exit
                end do

                random_num1 = rand()
                initial_rotational_speed = sqrt(initial_rotational_energy/mass_hydrogen)
                initial_rotation_angle = random_num1*pi2 - pi
                bond_period_elapsed = rand()

                !All of this is stored for later use in the InitialSetup of runTrajectory
                INITIAL_BOND_DATA(:,n) = (/ initial_bond_distance,initial_rotational_speed,&
                                        initial_rotation_angle,initial_bond_angle1,initial_bond_angle2,&
                                        bond_period_elapsed /)
	end do

	open(filechannel1,file=gridpath0//Ngrid_text//prefix_text//initialfile,&
                          position="append")
	write(filechannel1,FMTinitial) ((INITIAL_BOND_DATA(j,i),j=1,6),i=1,Nbonds)
        close(filechannel1)


	!Each trajectory will have Ngrid_total outputs; one for however many grids we use
	!The trajectory number will uniquely identify one trajectory from another
	write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
	write(Ntraj_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") n_testtraj
	
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

Ntraj = Ntesttraj
max_TranslationalEnergy = 0.0d0
max_absenergychange = 0.0
min_absenergychange = 1.0e9
max_relenergychange = 0.0
min_relenergychange = 1.0e9
max_rotenergychange = 0.0
min_rotenergychange = 1.0e9
call postProcess(Ngrid_text//prefix_text)

!We need to redo the postProcess (and convergence) of each grid
!in case the maximum and minimum of the trajectories just processed are different
if (.true.) then                     !(trueSA_flag) then
	call itime(now)
	write(6,FMT=FMTnow) now
	print *, "   Making plot: ", "InitialSATRVDistribution"
	print *, ""

	old_filename = ""
	do Ngrid = 1, Ngrid_total
		write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
		write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid

		!First, post process
	        Ntraj = Ntraj_max
        	call postProcess(Ngrid_text//"/Initial")

	        Ntraj = Ngrid*Ntraj_max
		write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
		write(Ntraj_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ntraj

	        if (trim(adjustl(old_filename)) /= "") then
	                new_filename = gridpath0//Ngrid_text//"/Initial"//Ntraj_text//SATRVfile
	                call system("cat "//trim(adjustl(old_filename))//" "//&
	                            gridpath0//Ngrid_text//"/Initial"//SATRVfile//" > "//new_filename)
	                old_filename = new_filename
	        else
	                old_filename = gridpath0//Ngrid_text//"/Initial"//Ntraj_text//SATRVfile
	                call system("cp "//gridpath0//Ngrid_text//"/Initial"//SATRVfile//" "//trim(adjustl(old_filename)))
	        end if

	        !Make an SA + TRV plot
	        call getScatteringAngles1(Ngrid_text//"/Initial"//Ntraj_text,"SATRVDistribution")
	end do
end if

!Finally, do a post-creation timeslice-to-SA conversions here
!We use the SA often so we do this at the beginning
write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid_total
Ntraj = Ntesttraj
write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
write(Ntraj_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ntraj

call getConvergenceImage(0.0, real(pi), 1, "ScatteringAngle")
call getConvergenceImage(min(min_relenergychange,min_absenergychange,min_rotenergychange),&
                         max(max_relenergychange,max_relenergychange,max_rotenergychange),4,"RelativeEnergyChange")
call getConvergenceImage(min(min_relenergychange,min_absenergychange,min_rotenergychange),&
                         max(max_relenergychange,max_relenergychange,max_rotenergychange),3,"AbsoluteEnergyChange")
call getConvergenceImage(min(min_relenergychange,min_absenergychange,min_rotenergychange),&
                         max(max_relenergychange,max_relenergychange,max_rotenergychange),5,"RotationalEnergyChange")

if (percentthreshold_flag) then
	call itime(now)
	write(6,FMT=FMTnow) now
	print *, "   Making plot: ", "PercentRMSDThreshold_"//Ngrid_text//prefix_text
	print *, ""
	call getRMSDThresholds1(prefix_text,"PercentRMSDThreshold_"//Ngrid_text//prefix_text)
end if

if (testtrajSA_flag) then
	call itime(now)
	write(6,FMT=FMTnow) now
	print *, "   Making plot: ", "TestScatteringAngleDistribution_"//Ngrid_text//prefix_text
	print *, ""
	
	call getScatteringAngles2(Ngrid_text//prefix_text,"TestScatteringAngleDistribution_"//&
                                  Ngrid_text//prefix_text)
end if

if (testheatmapSA_flag) then
	call itime(now)
	write(6,FMT=FMTnow) now
	print *, "   Making plot: ", "HeatMap_"//trim(adjustl(Ntraj_text))//"SATRVDistribution"
	print *, ""

	call getScatteringAngles1(Ngrid_text//prefix_text, "SATRVDistribution")
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
call itime(now)
write(6,FMT=FMTnow) now
print *, "Successfully exited analysis"
print *, ""

end program checkNewTrajectorieswithMultipleGrids
