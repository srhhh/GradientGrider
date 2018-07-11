module analyzeScatteringAngleswithMultipleGrids
implicit none

contains


subroutine getScatteringAngles1(gridpath0,Ngrid_cap,DATfilename,scattering_angle_column,JPGfilename)
use PARAMETERS
implicit none
integer,parameter :: dp = kind(1.0d0)

!PATH TO MULTIPLE LIBRARIES
character(*), intent(in) :: gridpath0

!CAP TO NUMBER OF GRIDS ANALYZED
integer, intent(in) :: Ngrid_cap
integer :: Ngrid_total, Ngrid

!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
character(*), intent(in) :: DATfilename

!COLUMN OF DAT FILE WITH SCATTERING ANGLES
integer, intent(in) :: scattering_angle_column

!FORMAT OF JPG FILES TO BE MADE
character(*), intent(in) :: JPGfilename

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectories_text_length*100) :: trajectories_text
character(6) :: Ntraj_text

!I/O HANDLING
integer :: iostate


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
                            DATfilename
end do
close(trajectorieschannel)

!We increment before we read, so we need to subtract out one increment here
if (Ngrid_total < 2) return
Ngrid_total = Ngrid_total - 1

!For now, we are putting a soft cap on how many trajectories we are using
!Each grid has Ntraj_max trajectories so we will use a maximum of 4 * Ntraj_max trajectories
Ngrid_total = min(Ngrid_cap, Ngrid_total)
print *, "Working on directory ", gridpath0
print *, "Deciding on using ", Ngrid_total, " grids"
print *, ""


!This data was made during creation (or should have been!) so all we need to do
!is merge, read, and plot them
do Ngrid = 1, Ngrid_total

	write(variable_length_text,FMT="(I5)") Ngrid_text_length
        !The folders are named starting from 001 by increments of 1
        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
        print *, " Working on grid number: ", Ngrid_text

        !The plots are named starting from Ntraj_max by increments of Ntraj_max (the number of trajectories)
        write(Ntraj_text,FMT="(I0.6)") Ngrid * Ntraj_max

        !This system call concatenates all the data files from that previously made 'hack'
        !By doing that, we merge all the scattering angle data together
        call system("cat"//trajectories_text(1:Ngrid*trajectories_text_length)//" >> "//&
                    gridpath0//Ngrid_text//"/"//cumulativefile//Ntraj_text//".dat")

        !This is the gnuplot code to make the plots
        open(gnuplotchannel,file=gridpath0//gnuplotfile)
        write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
        write(gnuplotchannel,*) 'set output "'//gridpath0//JPGfilename//Ntraj_text//'"'
        write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
        write(gnuplotchannel,*) 'unset key'
        write(gnuplotchannel,*) 'bin_width = 0.001'
        write(gnuplotchannel,*) 'bin_number(x) = floor(x/bin_width)'
        write(gnuplotchannel,*) 'rounded(x) = bin_width * (bin_number(x) + 0.5)'
        write(gnuplotchannel,*) 'set xlabel "Scattering Angle"'
        write(gnuplotchannel,*) 'set ylabel "Occurence"'
	write(variable_length_text,FMT="(I5)") scattering_angle_column
        write(gnuplotchannel,*) 'plot "'//gridpath0//Ngrid_text//'/'//cumulativefile//Ntraj_text//&
                                '.dat" u (rounded($'//trim(adjustl(variable_length_text))//&
				')):(1.0) smooth frequency with boxes'
        close(gnuplotchannel)

        !And then we just input it into gnuplot.exe
        call system("gnuplot < "//gridpath0//gnuplotfile)

	!And remove any extra files
	call system("rm "//gridpath0//Ngrid_text//"/"//cumulativefile//Ntraj_text//".dat")
	call system("rm "//gridpath0//gnuplotfile)

end do

end subroutine getScatteringAngles1


subroutine getScatteringAngles2(gridpath0,DATfilename,scattering_angle_column,JPGfilename)
use PARAMETERS
use ANALYSIS
implicit none

!PATH TO MULTIPLE LIBRARIES
character(*), intent(in) :: gridpath0

!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
character(*), intent(in) :: DATfilename

!COLUMN OF DAT FILE WITH SCATTERING ANGLES
integer, intent(in) :: scattering_angle_column

!FORMAT OF JPG FILES TO BE MADE
character(*), intent(in) :: JPGfilename

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectories_text_length*100) :: trajectories_text
character(6) :: Ntraj_text


!This is the gnuplot code to make the plots
open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath0//JPGfilename//'"'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'bin_width = 0.001'
write(gnuplotchannel,*) 'bin_number(x) = floor(x/bin_width)'
write(gnuplotchannel,*) 'rounded(x) = bin_width * (bin_number(x) + 0.5)'
write(gnuplotchannel,*) 'set xlabel "Scattering Angle"'
write(gnuplotchannel,*) 'set ylabel "Occurence"'
write(variable_length_text,FMT="(I5)") scattering_angle_column
write(gnuplotchannel,*) 'plot "'//gridpath0//DATfilename//&
                        '" u (rounded($'//trim(adjustl(variable_length_text))//&
			')):(1.0) smooth frequency with boxes'
close(gnuplotchannel)

!And then we just input it into gnuplot.exe
call system("gnuplot < "//gridpath0//gnuplotfile)
call system("rm "//gridpath0//gnuplotfile)

end subroutine getScatteringAngles2


end module analyzeScatteringAngleswithMultipleGrids
