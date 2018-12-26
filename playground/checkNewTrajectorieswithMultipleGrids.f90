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
!		FILECHANNELS(1)          	none
!		FILECHANNELS(2:1+Ngrid_total)	OPEN, RUNTRAJECTORY, CLOSE
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
character(12+Ngrid_text_length) :: prefix_text
character(150) :: old_filename, new_filename
character(1) :: answer

real(dp) :: scattering_angle
real(dp),dimension(3) :: TRVenergies1,TRVenergies2,dTRVenergies
real(dp),dimension(3,Natoms) :: coords_initial,velocities_initial,coords_final,velocities_final
integer :: seed,n,m,n_testtraj,initial_n_testtraj
real :: lowerlimit,upperlimit

!Variables
integer :: iostate
integer,allocatable :: filechannels(:)
integer :: OMP_GET_THREAD_NUM
logical :: return_flag = .false.

!Timing
real :: r1, r2, system_clock_rate
integer :: c1, c2, cr
integer,dimension(3) :: now
real :: grid_wall_time,checktrajectory_wall_time
character(10) :: grid_wall_time_text, checktrajectory_wall_time_text
integer :: trajectory_t0, trajectory_t1

!Incremental Integers
integer :: i, j, k


call OMP_SET_NUM_THREADS(Nthreads)

!Now here we actually make these new trajectories
!$OMP PARALLEL DEFAULT(none)&
!$OMP& PRIVATE(iostate)&
!$OMP& PRIVATE(variable_length_text,Ntraj_text,filechannels)&
!$OMP& FIRSTPRIVATE(Ngrid_text,reject_text)&
!$OMP& PRIVATE(c1,c2,cr,r1,r2,coords_initial,coords_final,velocities_initial,velocities_final)&
!$OMP& SHARED(return_flag,now,seed,system_clock_rate)&
!$OMP& SHARED(initial_n_testtraj,folder_text,prefix_text,Nthreshold_text,old_filename,new_filename,answer)&
!$OMP& PRIVATE(i,j,n)

!$OMP SINGLE

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

!$OMP END SINGLE

write(Nbond_text,FMT="(I0.6)") Nbonds
write(Natom_text,FMT="(I0.6)") Natoms
write(FMTinitial,FMT="(A19)") "("//Nbond_text//"(6(F14.10)))"
write(FMTtimeslice,FMT="(A19)") "("//Natom_text//"(12(F12.7)))"
write(FMT2,FMT="(A22)") "("//Natom_text//"(6(1x,F14.10)))"
write(FMT3,FMT="(A22)") "("//Natom_text//"(3(1x,F14.10)))"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       PRE-PROCESSING LIBRARY OVERVIEW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP CRITICAL

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

!$OMP END CRITICAL
!$OMP BARRIER

!$OMP SINGLE

print *, ""
call itime(now)
write(6,FMT=FMTnow) now
print *, "Analysis on directory ", gridpath0
print *, "Deciding on using ", Ngrid_total, " grids"
print *, "Deciding on using ", Nthreads, " threads"
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Some intitialization stuff
write(Nthreshold_text,FMT=FMT6_pos_real0) threshold_rmsd
if (reject_flag) then
	reject_text = "reject"
else
        if (accept_first) then
                 if (accept_worst) then
                         reject_text = "alphaW"
                 else
                         reject_text = "alphaA"
                 end if
        else
                 if (accept_worst) then
                         reject_text = "omegaW"
                 else
                         reject_text = "omegaA"
                 end if
        end if
end if


write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid_total

prefix_text = Ngrid_text//reject_text//Nthreshold_text

!$OMP END SINGLE
!$OMP BARRIER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       TRAJECTORY CREATION FLAG START
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!This is for checking trajectories against the grid
!Currently, this is the main use of this program
if (testtraj_flag) then

!$OMP SINGLE

print *, ""
call itime(now)
write(6,FMT=FMTnow) now
print *, "Trajectory Creation ... Start"

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
        write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
        call system("ls "//gridpath0//Ngrid_text//"/ | grep -E '^"//Ngrid_text//&
                    reject_text//"\"//Nthreshold_text//"_[0123456789]{"//trim(adjustl(variable_length_text))//&
                    "}.dat' > "//gridpath0//trajectories)

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
        open(trajectorieschannel,file=gridpath0//prefix_text//timeslicefile)
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
                                                         prefix_text//" ... .dat"
                if (Ntraj < initial_n_testtraj) print *, "    timeslice file corrupted for "//&
                                                         prefix_text
                print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		print *, ""
		
		if ((Ntraj < Ntesttraj).or.(initial_n_testtraj < Ntesttraj)) then
		do
			print *, "    corrupted data interrupts trajectory indexing;"
			print *, "    overwrite corrupted data? otherwise the program exits    (y/n)"
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

	                return_flag = .true.
		else
			print *, "    starting off where the timeslice file left off!"
			print *, ""
			initial_n_testtraj = min(initial_n_testtraj,Ntraj)

		end if
		else
			print *, "    going along with analysis on uncorrupted data"
			print *, ""
		end if
                print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		print *, ""
	end if

        if (.not.return_flag) then
        print *, "     Total Number of Trajectories Saved: ", initial_n_testtraj - 1, " out of ", Ntesttraj
        print *, ""
        end if

!If not...
else

        !We need to make a new one of these trajectory files altogether
        call system("rm "//gridpath0//prefix_text//trajectoriesfile)
        call system("rm "//gridpath0//prefix_text//timeslicefile)
        call system("rm "//gridpath0//prefix_text//"traversal.dat")
        if (prefix_text /= initialbondname) call system("rm "//gridpath0//prefix_text//initialfile)
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

!If we want to use the same initial conditions as an old set of trajectories
!open up whatever initial file that sample corresponds to
if ((useoldinitialbonddata_flag).and.(.not.return_flag)) then
        open(trajectorieschannel,file=gridpath0//initialbondname//initialfile)
        print *, ""
        print *, "     Deciding to use old initial conditions..."
        print *, ""

        !Check to see if we will run out of trajectories
        do n_testtraj = initial_n_testtraj, Ntesttraj
                read(trajectorieschannel,FMT=FMTinitial,iostat=iostate) &
                            ((INITIAL_BOND_DATA(j,i),j=1,6),i=1,Nbonds)

                !If so, we just exit out
                if (iostate /= 0) then
                        print *, ""
                        call itime(now)
                        write(6,FMT=FMTnow) now
                        print *, "Error! Not enough initial conditions in initialfile"
                        print *, "Total number of initial conditions inside of it: ", n_testtraj-initial_n_testtraj
                        print *, "Exiting analysis"
                        print *, ""
        
                        initial_n_testtraj = Ntesttraj+1
                        return_flag = .true.
                end if
        end do

        close(trajectorieschannel)
        open(trajectorieschannel,file=gridpath0//initialbondname//initialfile)
end if

!$OMP END SINGLE
!$OMP BARRIER

!If we are checking traversal, allocate the necessary arrays
if (traversal_flag) allocate(traversal0(Ngrid_total,counter0_max),&
                             traversal1(Ngrid_total,counter1_max))

!BEFORE going inside the parallel do loop, each thread should
!instantiate its own filechannels
allocate(filechannels(1+Ngrid_total))

!$OMP DO
do n_testtraj = initial_n_testtraj, Ntesttraj

        !If we're reusing old initial conditions, that's easy; we just read off the file
        if (useoldinitialbonddata_flag) then
                !$OMP CRITICAL
                read(trajectorieschannel,FMT=FMTinitial,iostat=iostate) &
                            ((INITIAL_BOND_DATA(j,i),j=1,6),i=1,Nbonds)
                !$OMP END CRITICAL
        
	!Otherwise, we must creat a random initial trajectory from the ensemble
        !specified in PHYSICS and PARAMETERS
        else
                call InitialSampling3()
        end if

        !In most case, we will need to record the initial conditions to a (potentially) new file
        !$OMP CRITICAL
        if ((.not.(useoldinitialbonddata_flag)).or.(prefix_text /= initialbondname)) then
        	open(filechannel1,file=gridpath0//prefix_text//initialfile,&
                                  position="append")
        	write(filechannel1,FMTinitial) ((INITIAL_BOND_DATA(j,i),j=1,6),i=1,Nbonds)
                close(filechannel1)
        end if
        !$OMP END CRITICAL

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
        filechannels(1) = 1000 + OMP_GET_THREAD_NUM()
	do Ngrid = 1, Ngrid_total
		write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
		filechannels(1+Ngrid) = 1000 + 69 * Ngrid + OMP_GET_THREAD_NUM()
		open(filechannels(1+Ngrid),file=gridpath0//Ngrid_text//"/"//&
                                              prefix_text//'_'//Ntraj_text//".dat")
	end do

	!Then write the outputted RMSDS of each trajectory onto those filechannels
	!Remark: checkMultipleGrids uses filechannel1 to open files in the grid
	call checkMultipleTrajectories(filechannels(1:1+Ngrid_total),&
                                                                     coords_initial,velocities_initial,&
                                                                     coords_final,velocities_final)

	!Also let's see how long a single trajectory takes
	call system_clock(c2)
	call CPU_time(r2)
	print *, " Trajectory "//Ntraj_text//"     CPU Time: ", r2 - r1
	print *, " Trajectory "//Ntraj_text//"    Wall Time: ", (c2 - c1) * system_clock_rate

	!Finally, close them
	do Ngrid = 1, Ngrid_total
		close(filechannels(1+Ngrid))
	end do

        !$OMP CRITICAL
	!The only data that is recorded is the first and last frame of the trajectory
	open(filechannel1,file=gridpath0//prefix_text//timeslicefile,position="append")
	write(filechannel1,FMTtimeslice) &
                              ((coords_initial(i,j),i=1,3),j=1,Natoms),&
                              ((velocities_initial(i,j),i=1,3),j=1,Natoms),&
			      ((coords_final(i,j),i=1,3),j=1,Natoms),&
                              ((velocities_final(i,j),i=1,3),j=1,Natoms)
        close(filechannel1)

        !For traversal, we can see how many cells the trajectory traversed
        !merely by figuring out how many nonzero values exist (since 0s represent no traversal)
        if (traversal_flag) then
                open(filechannel1,file=gridpath0//prefix_text//"traversal.dat",position="append")
                write(filechannel1,FMT=*) steps, ((sum(traversal1(i,:))),i=1,Ngrid_total)
                close(filechannel1)
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       DURING-CREATION TRAJECTORY ANALYSIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !If this flag is one, we are also checking the RMSD of a trajectory per timestpe
        !This is very expensive so please don't turn this on
        !(and you have to manually turn it on in checkMultipleTrajectories anyway)
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
		write(gnuplotchannel,*) 'plot "'//gridpath0//Ngrid_text//"/"//&
                                        prefix_text//'_'//Ntraj_text//'.dat'//&
		                        '" u (rounded($1)):(1.0) smooth frequency with boxes'
		close(gnuplotchannel)
		call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)
	end if
        !$OMP END CRITICAL
	
end do
!$OMP END DO
print *, ""

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       TRAJECTORY CREATION FLAG END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!If we did use old initial conditions, close that file up
!$OMP SINGLE
if (useoldinitialbonddata_flag) close(trajectorieschannel)
!$OMP END SINGLE

!This closes that big, enclosing if statement on whether to make trajectories
end if

!$OMP END PARALLEL

if (return_flag) then
        return
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       POST-CREATION LIBRARY ANALYSIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Specify how many trajectories we will be analyzing
!We communicate to analyzeScatteringAngles and others through the
!shared variable Ntraj
Ntraj = Ntesttraj

!Initialize each maximum and minimum value
!These are also shared among modules
max_TranslationalEnergy = 0.0d0
max_absenergychange = 0.0
min_absenergychange = 1.0e9
max_relenergychange = 0.0
min_relenergychange = 1.0e9
max_rotenergychange = 0.0
min_rotenergychange = 1.0e9

!If we have a list of files we need to compare, then we need to comb through those
if (comparison_flag) then

        if (comparison_number < 1) then
                print *, ""
                call itime(now)
                write(6,FMT=FMTnow) now
                print *, "Number of trajectory sets being compared is inadequate"
                print *, "Terribly exited analysis"
                print *, ""

                return
        end if

        !These files are stored in allprefixes and the length of each file
        !is stored in alllengths
        call postProcess(allprefixes(1:alllengths(1)))
        do i = 1, comparison_number-1
                call postProcess(allprefixes(1+sum(alllengths(1:i)):sum(alllengths(1:i+1))))
        end do

!Otherwise, we do not need to process all the files in the library
else
        !We only need to look at the files if we are actually analyzing something
        if (percentthreshold_flag .or. testtrajSA_flag .or. testheatmapSA_flag) then
             call postProcess(prefix_text)

        !Otherwise, we can just exit early
        else
                print *, ""
                call itime(now)
                write(6,FMT=FMTnow) now
                print *, "Successfully exited analysis"
                print *, ""
   
                return
        end if
end if

!We need to redo the postProcess (and convergence) of each grid as well
!in case the maximum and minimum of the trajectories just processed are different
if (comparison_flag .or. trueSA_flag .or. testtrajSA_flag .or. testheatmapSA_flag) then
	call itime(now)
	write(6,FMT=FMTnow) now
	print *, "   Making plots: ", "Initial_SATRVDistribution"
	print *, "   Making plots: ", "Initial_HeatMap_SATRVDistribution"
	print *, ""

	old_filename = ""
	do Ngrid = 1, Ngrid_max
                !Prepare the name of the file; here we are uniquely identifying the
                !file by how many trajectories we've accumulated and put in it
	        Ntraj = Ngrid*Ntraj_max
		write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
		write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
		write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
		write(Ntraj_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ntraj

		!Post process the file; here we are telling the subroutine that
                !there are Ntraj_max trajectories in the file
	        Ntraj = Ntraj_max
        	call postProcess(Ngrid_text//"/Initial")

                !We create this unique file that has all trajectories up until that grid
                !so beforehand we need to concatenate all the old trajectories with the new file
	        if (trim(adjustl(old_filename)) /= "") then
	                new_filename = gridpath0//Ngrid_text//"/Initial"//Ntraj_text//SATRVfile
	                call system("cat "//trim(adjustl(old_filename))//" "//&
	                            gridpath0//Ngrid_text//"/Initial"//SATRVfile//" > "//new_filename)
	                old_filename = new_filename

                !If this is the first grid, then we just copy the new file to this unique file
	        else
	                old_filename = gridpath0//Ngrid_text//"/Initial"//Ntraj_text//SATRVfile
	                call system("cp "//gridpath0//Ngrid_text//"/Initial"//SATRVfile//" "//trim(adjustl(old_filename)))
	        end if

	        !Then bin it
	        Ntraj = Ngrid*Ntraj_max
	        call getScatteringAngles1(Ngrid_text//"/Initial"//Ntraj_text,"SATRVDistribution")
	end do

        !Finally, we look at the scattering angle and energy change distributions of the grid
        !with the maximums and minimums gathered from above
        !and see whether they converge; these plots will be used for reference
        Ntraj = Ntesttraj
        write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid_total
        write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
        write(Ntraj_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ntraj
        
        call getConvergenceImage(0.0,real(pi), 1, "ScatteringAngle")
        call getConvergenceImage(min_absenergychange,max_absenergychange, 3, "AbsoluteEnergyChange")
        call getConvergenceImage(min_relenergychange,max_relenergychange, 4, "RelativeEnergyChange")
        call getConvergenceImage(min_rotenergychange,max_rotenergychange, 5, "RotationalEnergyChange")

end if

!In a comparison run, we don't look at any of the other flags
if (comparison_flag) then

        !The name the user gives us tells us
        !  1) which column of the SATRV file to look in
        !  2) how to name our output file
        !The user may also give user-defined bounds to the distribution
        if (trim(adjustl(comparison_SATRVname)) == "ScatteringAngle") then
                comparison_SATRVcolumn = 1
                lowerlimit = 0.0
                upperlimit = real(pi)
        else if (trim(adjustl(comparison_SATRVname)) == "AbsoluteEnergyChange") then
                comparison_SATRVcolumn = 3
                lowerlimit = min_absenergychange
                upperlimit = max_absenergychange
        else if (trim(adjustl(comparison_SATRVname)) == "RelativeEnergyChange") then
                comparison_SATRVcolumn = 4
                lowerlimit = min_relenergychange
                upperlimit = max_relenergychange
        else if (trim(adjustl(comparison_SATRVname)) == "RotationalEnergyChange") then
                comparison_SATRVcolumn = 5
                lowerlimit = min_rotenergychange
                upperlimit = max_rotenergychange
        else
        end if

        !If need be, use the upper and lower limit dictated by the user
        if (comparison_upperlimit /= comparison_lowerlimit) then
                call getConvergenceImage(comparison_lowerlimit,comparison_upperlimit,comparison_SATRVcolumn,&
                                         trim(adjustl(comparison_SATRVname)))
        end if

        call itime(now)
        write(6,FMT=FMTnow) now
	print *, "   Making plot: ", "Comparison_"//trim(adjustl(Ntraj_text))//trim(adjustl(comparison_SATRVname))
	print *, ""

        !We now need to bin each set of trajectories we want to compare
        call getScatteringAngles1(allprefixes(1:alllengths(1)),"")
        do i = 1, comparison_number-1
                call getScatteringAngles1(allprefixes(1+sum(alllengths(1:i)):sum(alllengths(1:i+1))),"")
        end do

        !If the upper and lower limit are equal, that is an internal sign that means
        !that the user wants to use the natural bounds of the distribution
        if (comparison_upperlimit /= comparison_lowerlimit) then
                call getComparedScatteringAngles(comparison_lowerlimit,comparison_upperlimit,&
                         "Comparison_"//trim(adjustl(Ntraj_text))//trim(adjustl(comparison_SATRVname)),&
                         comparison_SATRVcolumn,trim(adjustl(comparison_SATRVname)))

        !Otherwise, that means we must use the user-defined bounds
        else
                call getComparedScatteringAngles(lowerlimit,upperlimit,&
                         "Comparison_"//trim(adjustl(Ntraj_text))//trim(adjustl(comparison_SATRVname)),&
                         comparison_SATRVcolumn,trim(adjustl(comparison_SATRVname)))
        end if

        print *, ""
        call itime(now)
        write(6,FMT=FMTnow) now
        print *, "Successfully exited analysis"
        print *, ""

        return
end if



!In a run without comparison, usually we do some analysis on the trajectories afterwards
!The following calls do some data manipulation and make some figures

if (percentthreshold_flag) then
	call itime(now)
	write(6,FMT=FMTnow) now
	print *, "   Making plot: ", prefix_text//"_PercentRMSDThreshold"
	print *, ""

        if (percentthreshold_key < 1) then
	        call getRMSDThresholds1(prefix_text,prefix_text//"_PercentRMSDThreshold")
        else
	        call getRMSDThresholds1(prefix_text,prefix_text//"_PercentRMSDThreshold",&
                                        percentthreshold_key)
        end if
end if

if (testheatmapSA_flag) then
	call itime(now)
	write(6,FMT=FMTnow) now
	print *, "   Making plot: ", prefix_text//"_SATRVDistribution"
	print *, "   Making plot: ", prefix_text//"_HeatMap_SATRVDistribution"
	print *, ""

	call getScatteringAngles1(prefix_text, "SATRVDistribution")
end if

if (testtrajSA_flag) then
	call itime(now)
	write(6,FMT=FMTnow) now
	print *, "   Making plot: ", prefix_text//"_Final_SATRVDistribution"
	print *, ""
	
	call getScatteringAngles2(prefix_text,prefix_text//"_Final_SATRVDistribution")
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print *, ""
call itime(now)
write(6,FMT=FMTnow) now
print *, "Successfully exited analysis"
print *, ""

return

end program checkNewTrajectorieswithMultipleGrids
