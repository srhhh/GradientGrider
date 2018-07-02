program analyzeMultipleTrajectories
use f2_physics_parameters
use f2_parameters
use checkCells5
use ls_rmsd
use mapCellData
use makeTrajectory5
use makeTrajectory
implicit none

character(4) :: resolution_text
character(5) :: overcrowd0_text
character(5) :: overcrowd1_text
character(4) :: Ngrid_text
character(len(path5)+4+5+5+3) :: gridpath0 
character(23),parameter :: cumulativefile = "cumulative_trajectories"
integer, parameter :: trajectories_text_length = len(path5)+4+5+5+3 +&
					         4 + len(trajectoriesfile) + 1
character(trajectories_text_length*100) :: trajectories_text
character(6) :: Ntraj_text

integer,parameter :: Ntesttraj = 100
real,dimension(Ntesttraj) :: initial_bond_distance, initial_rotation_angle, initial_rotational_speed
real,dimension(Ntesttraj) :: initial_bond_angle1, initial_bond_angle2
real,dimension(Ntesttraj) :: initial_energy_H2,initial_vibrational_energy,initial_rotational_energy
real :: random_num1,random_num2,i,j
integer :: seed,n,m,n_testtraj

integer :: Ngrid,iostate,Ngrid_max
integer,allocatable :: grid_numbers(:)
integer :: Ntraj_max = 200
integer,allocatable :: filechannels(:)
real,allocatable :: min_rmsds(:),percent_threshold_rmsd(:)
real :: min_rmsd
integer :: frames, step, total_threshold_rmsd
real,parameter :: threshold_rmsd = 1.0e-5

write(resolution_text,FMT="(I0.4)") resolution_0
write(overcrowd0_text,FMT="(I0.5)") overcrowd0
write(overcrowd1_text,FMT="(I0.5)") overcrowd1
gridpath0 = path5//resolution_text//"_"//overcrowd0_text//"_"//overcrowd1_text//"/"

!All trajectory folders are formatted as a number
!So search for these numbered folders and read them
call system("ls -p "//gridpath0//" | grep '[0123456789]/' > "//gridpath0//trajectories)

open(trajectorieschannel,file=gridpath0//trajectories,action="read")
trajectories_text = ""
do Ngrid_max = 1, 100
	read(trajectorieschannel,FMT="(A4)",iostat=iostate) Ngrid_text
	if (iostate /= 0) exit
	trajectories_text = trim(trajectories_text)//" "//gridpath0//Ngrid_text//&
			    trajectoriesfile
end do
close(trajectorieschannel)

if (Ngrid_max < 3) return
Ngrid_max = Ngrid_max - 2

!Get Ntesttraj number of random initial conditions
call system_clock(seed)
print *, "Working with system_clock seed: ", seed
print *, ""
seed = rand(seed)
do n_testtraj = 1, Ntesttraj
        random_num1 = rand()
        random_num2 = rand()
        initial_bond_angle1(n_testtraj) = random_num1*pi2           !theta
        initial_bond_angle2(n_testtraj) = random_num2*pi2           !phi
        do
        random_num1 = rand()
        random_num2 = rand()
        initial_energy_H2(n_testtraj) = (upsilon_max*random_num1 + 0.5)*upsilon_factor2
        if (random_num2 < temperature_scaling*exp(upsilon_max*random_num1*upsilon_factor1)) exit
        end do
        random_num2 = 1.0
        initial_vibrational_energy(n_testtraj) = random_num2*initial_energy_H2(n_testtraj)
        initial_rotational_energy(n_testtraj) = initial_energy_H2(n_testtraj) - initial_vibrational_energy(n_testtraj)
        initial_bond_distance(n_testtraj) = HOr0_hydrogen + sqrt(initial_vibrational_energy(n_testtraj)*2/HOke_hydrogen)
        random_num1 = rand()
        initial_rotational_speed(n_testtraj) = sqrt(initial_rotational_energy(n_testtraj)/mass_hydrogen)
        initial_rotation_angle(n_testtraj) = random_num1*2*pi
end do


print *, "Finished Initialization Part 1..."

do Ngrid = 1, Ngrid_max
	write(Ngrid_text,FMT="(I0.3)") Ngrid
	Ngrid_text = trim(adjustl(Ngrid_text))//"/"
	write(Ntraj_text,FMT="(I0.6)") Ngrid * Ntraj_max

print *, " Working on trajectories: ", Ntraj_text

	call system("cat"//trajectories_text(1:Ngrid*trajectories_text_length)//" >> "//&
		    gridpath0//Ngrid_text//cumulativefile//Ntraj_text//".dat")

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
        write(gnuplotchannel,*) 'plot "'//gridpath0//Ngrid_text//cumulativefile//Ntraj_text//&
				'.dat" u (rounded($7)):(7) smooth frequency with boxes'
	close(gnuplotchannel)

	call system("gnuplot < "//gridpath0//gnuplotfile)

end do

print *, "Finished Initialization Part 2..."

Ngrid_max = min(10, Ngrid_max)
allocate(filechannels(Ngrid_max))
do n_testtraj = 1, Ntesttraj
	write(Ntraj_text,FMT="(I0.6)") n_testtraj

	do Ngrid = 1, Ngrid_max
		write(Ngrid_text,FMT="(I0.3)") Ngrid
		Ngrid_text = trim(adjustl(Ngrid_text))//"/"
		filechannels(Ngrid) = 1000 + 69 * Ngrid
		open(filechannels(Ngrid),file=gridpath0//Ngrid_text//Ntraj_text//".dat")
	end do

print *, " Working on trajectories: ", Ntraj_text

	!Remark: checkCells5 uses filechannel1 to open files in the gri
	open(filechannel2,file=gridpath0//Ntraj_text//".dat")
	call checkMultipleTrajectories(initial_bond_distance(n_testtraj),&
		    initial_rotational_speed(n_testtraj),initial_rotation_angle(n_testtraj),&
                    initial_bond_angle1(n_testtraj),initial_bond_angle2(n_testtraj),.false.,&
                    Ngrid_max,filechannels(1:Ngrid_max),gridpath0)
	close(filechannel2)

	do Ngrid = 1, Ngrid_max
		close(filechannels(Ngrid))
	end do
end do

print *, "Finished Initialization Part 3..."
print *, ""

do Ngrid = 1, Ngrid_max
write(Ngrid_text,FMT="(I0.3)") Ngrid
Ngrid_text = trim(adjustl(Ngrid_text))//"/"
print *, " Working on trajectories: ", Ngrid_text

open(filechannel1,file=gridpath0//"percent_rmsd"//Ngrid_text//".dat")
do n_testtraj = 1, Ntesttraj
	write(Ntraj_text,FMT="(I0.6)") n_testtraj
	!Read the trajectory (which has the rmsd) across all grids line-by-line
	iostate = 0
	frames = 0
	total_threshold_rmsd = 0
	open(filechannel2,file=gridpath0//Ngrid_text//Ntraj_text//".dat")
	do
		read(filechannel2,FMT="(F12.8)",iostat=iostate) min_rmsds(n)
		if (iostate /= 0) exit
		frames = frames + 1

		if (min_rmsd < threshold_rmsd) total_threshold_rmsd = total_threshold_rmsd + 1
	end do
	close(filechannel2)

	percent_threshold_rmsd(n_testtraj) = total_threshold_rmsd * 1.0 / frames
	write(filechannel1,FMT="(I6,1x,F7.4,1x,I8)") n_testtraj, percent_threshold_rmsd(n_testtraj), frames
end do
close(filechannel1)

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
