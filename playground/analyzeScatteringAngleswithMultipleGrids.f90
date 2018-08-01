!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       MODULE
!               analyzeScatteringAngleDistributions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!		This module simply plots scattering angles
!		from specified DAT files
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILECHANNELS                    ACTION
!
!               GNUPLOTCHANNEL                  OPEN, WRITE, CLOSE
!               TRAJECTORIESCHANNEL             OPEN, WRITE, CLOSE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       CALLS                           MODULE
!
!               SYSTEM                          INTRINSIC
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module analyzeScatteringAngleswithMultipleGrids
implicit none

contains


subroutine getScatteringAngles1(Ngrid_cap,Nsamples,DATfilename,scattering_angle_column,JPGfilename)
use PARAMETERS
implicit none
integer,parameter :: dp = kind(1.0d0)

!CAP TO NUMBER OF GRIDS ANALYZED
integer, intent(in) :: Ngrid_cap
integer :: Ngrid_total

!RANDOM TRAJECTORY SAMPLING
integer,intent(in) :: Nsamples
real,allocatable :: scatteringAngles
integer,allocatable :: originalIndexes(:)
integer,dimension(Ntesttraj,Nsamples) :: selectIntegers, selectIndexes
integer,dimension(Nsamples) :: selectIndex

!SCATTERING ANGLE BINNING
real :: scatteringAngle
real :: bin_width, bin_accept
real :: binMean, binSD
real,dimension(Nsamples) :: binSize

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
character(6) :: Ntraj_text,boxwidth_text

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


!This data was made during creation (or should have been!) so all we need to do
!is merge, read, and plot them
do Ngrid = 1, Ngrid_total

	write(variable_length_text,FMT="(I5)") Ngrid_text_length
        !The folders are named starting from 001 by increments of 1
        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
        print *, " Working on grid number: ", Ngrid_text

        !The plots are named starting from Ntraj_max by increments of Ntraj_max (the number of trajectories)
        write(Ntraj_text,FMT="(I6)") Ngrid * Ntraj_max

        !This system call concatenates all the data files from that previously made 'hack'
        !By doing that, we merge all the scattering angle data together
        call system("cat"//trajectories_text(1:Ngrid*(trajectories_text_length+7))//" >> "//&
                    gridpath0//Ngrid_text//"/"//cumulativefile//trim(adjustl(Ntraj_text))//".dat")

        !This is the gnuplot code to make the plots
        open(gnuplotchannel,file=gridpath0//gnuplotfile)
        write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
        write(gnuplotchannel,*) 'set output "'//gridpath0//trim(adjustl(Ntraj_text))//JPGfilename//'"'
        write(gnuplotchannel,*) 'set title "Scattering Angle Distribution of '//trim(adjustl(Ntraj_text))//&
                                ' trajectories of '//gridpath0//'"'
        write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
        write(gnuplotchannel,*) 'unset key'
	write(boxwidth_text,FMT="(F6.5)") 3.14159 * 10 / max(real(Ngrid*Ntraj_max),50.0)
	write(gnuplotchannel,*) 'set boxwidth '//boxwidth_text
	write(gnuplotchannel,*) 'bin_width = '//boxwidth_text
	write(gnuplotchannel,*) 'bin_number(x) = floor(x/bin_width)'
        write(gnuplotchannel,*) 'rounded(x) = bin_width * (bin_number(x) + 0.5)'
        write(gnuplotchannel,*) 'set xlabel "Scattering Angle"'
        write(gnuplotchannel,*) 'set ylabel "Occurence"'
        write(gnuplotchannel,*) 'set xrange [0:3.14159]'
        write(gnuplotchannel,*) 'set y autoscale'
	write(variable_length_text,FMT="(I5)") scattering_angle_column
        write(gnuplotchannel,*) 'plot "'//gridpath0//Ngrid_text//'/'//cumulativefile//&
                                trim(adjustl(Ntraj_text))//&
                                '.dat" u (rounded($'//trim(adjustl(variable_length_text))//&
				')):(1.0) smooth frequency with boxes'
        close(gnuplotchannel)

        !And then we just input it into gnuplot.exe
        call system("gnuplot < "//gridpath0//gnuplotfile)

end do

!We also want to know how much variance there is for a particular
!sample size (a sample size of Ntesttraj in this case)

!For this, we need to know how to bin the scattering angles of ALL the trajectories

!First, get all the trajectories in an internal array
allocate(scatteringAngles(Ngrid_total*Ntraj_max),originalIndexes(Ngrid_total*Ntraj_max))
open(filechannel1,file=gridpath0//Ngrid_test//"/"//cumulativefile//trim(adjustl(Ntraj_text))//".dat")
do i = 1, Ngrid_total*Ntraj_max
        read(filechannel1,FMT=*) (line_data(j),j=1,scattering_angle_column)
        originalIndexes(i) = (/ line_data(1) /)
        scatteringAngles(i) = (/ line_data(scattering_angle_column) /)
end do
close(filechannel1)

!Then we need to sort this big list
call qsort2(scatteringAngles,originalIndexes,Ngrid_total*Ntraj,1,1)

!Now we need some random samples of this big list
!Here, we just choose a random Ntesttraj number of integers within bounds of the array
do i = 1, Ntesttraj
        selectIndexes(i,:) = i
end do
do i = 1, Nsamples
        call chooseINT(Ntesttraj,1,Ngrid_total*Ntraj,selectIntegers(:,i))
        call qsort2(selectIntegers(:,i),selectIndexes(:,i),Ntesttraj,1,1)
end do

!Now we must do the laborious job of binning these
!Each bin has its own average and standard deviation based on how many samples we took
bin_width = pi2 / Nbins
bin_accept = bin_width
selectIndex = 1
open(filechannel1,file=gridpath0//"Adjusted"//cumulativefile)
do i = 1, Nbins
        binSize = 0.0
        do j = 1, Nsamples
                if (selectIndex(j) > Ntesttraj) cycle
                scatteringAngle = scatteringAngles(originalIndexes(&
                                  selectIntegers(j,selectIndexes(j,selectIndex(j)))))
                if (scatteringAngle > bin_accept) cycle
                selectIndex(j) = selectIndex(j) + 1
                binSize(j) = binSize(j) + 1.0
        end do
        binMean = sum(binSize) / Nsamples
        binSD = sqrt(sum(((binSize-binMean)**2))/(Nsamples - 1))

        write(filechannel1,*) bin_accept, binMean, binSD
        bin_accept = bin_accept + bin_width
end do
close(filechannel1)

!We have everything we need to draw the distribution with error bars
open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) "set terminal jpeg size 1200,1200"
write(variable_length_text,FMT="(I6)") Ntesttraj
write(gnuplotchannel,*) 'set output "'//gridpath0//'Adjustedto'//&
                        trim(adjustl(variable_length_text)//JPGfilename//'"'
write(gnuplotchannel,*) 'set title "Scattering Angle Distribution for Samples of '//&
                        trim(adjustl(variable_length_text))//' Trajectories"'
write(gnuplotchannel,*) 'set xlabel "Scattering Angle (rad)"'
write(gnuplotchannel,*) 'set xrange [0:6.28318]'
write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set ylabel "Frequency"'
write(gnuplotchannel,*) 'set ytics 1'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'plot "'//gridpath0//'Adjusted'//cumulativefile//'" u 1:2 w lines, \'
write(gnuplotchannel,*) '     "'//gridpath0//'Adjusted'//cumulativefile//'" u 1:2:3 w yerrorbars'
close(gnuplotchannel)


end subroutine getScatteringAngles1


subroutine getScatteringAngles2(DATfilename,scattering_angle_column,theta_column,phi_column,JPGfilename)
use PARAMETERS
use ANALYSIS
implicit none

!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
character(*), intent(in) :: DATfilename

!COLUMN OF DAT FILE WITH SCATTERING ANGLES
integer, intent(in) :: scattering_angle_column,theta_column,phi_column

!FORMAT OF JPG FILES TO BE MADE
character(*), intent(in) :: JPGfilename

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectories_text_length*100) :: trajectories_text
character(6) :: Ntraj_text
character(6) :: boxwidth_text


!This is the gnuplot code to make the plots
open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath0//JPGfilename//'"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set xrange [0:3.14159]'
write(gnuplotchannel,*) 'unset xlabel'
write(gnuplotchannel,*) 'unset xtics'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
write(boxwidth_text,FMT="(F6.5)") 3.14159 * 10 / max(real(Ntraj),50.0)
write(Ntraj_text,FMT="(I6)") Ntraj
write(gnuplotchannel,*) 'set multiplot layout 3,1 margins 0.15,0.95,.1,.9 spacing 0,0 title '//&
                        '"Angle Distribution of '//trim(adjustl(Ntraj_text))//' trajectories of '//gridpath0//'"'
write(Ntraj_text,FMT="(I6)") Ntesttraj * 10 / 10
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'set boxwidth '//boxwidth_text
write(gnuplotchannel,*) 'bin_width = '//boxwidth_text
write(gnuplotchannel,*) 'bin_number(x) = floor(x/bin_width)'
write(gnuplotchannel,*) 'rounded(x) = bin_width * (bin_number(x) + 0.5)'
write(gnuplotchannel,*) 'set ylabel "Scattering Angle Occurence"'
write(gnuplotchannel,*) 'set yrange [0:]'
write(variable_length_text,FMT="(I5)") scattering_angle_column
write(gnuplotchannel,*) 'plot "'//gridpath0//DATfilename//&
                        '" u (rounded($'//trim(adjustl(variable_length_text))//&
			')):(1.0) smooth frequency with boxes'
write(gnuplotchannel,*) 'set ylabel "Initial H2 Theta Occurence"'
write(gnuplotchannel,*) 'set yrange [0:]'
write(variable_length_text,FMT="(I5)") theta_column
write(gnuplotchannel,*) 'plot "'//gridpath0//DATfilename//&
                        '" u (rounded($'//trim(adjustl(variable_length_text))//&
			')):(1.0) smooth frequency with boxes'
write(gnuplotchannel,*) 'set xlabel "Angle (rad)"'
write(gnuplotchannel,*) 'set xtics'
write(gnuplotchannel,*) 'set ylabel "Initial H2 Phi Occurence"'
write(gnuplotchannel,*) 'set yrange [0:]'
write(variable_length_text,FMT="(I5)") phi_column
write(gnuplotchannel,*) 'plot "'//gridpath0//DATfilename//&
			'" u (rounded($'//trim(adjustl(variable_length_text))//&
			')):(1.0) smooth frequency with boxes'

close(gnuplotchannel)

!And then we just input it into gnuplot.exe
call system("gnuplot < "//gridpath0//gnuplotfile)
call system("rm "//gridpath0//gnuplotfile)

end subroutine getScatteringAngles2


end module analyzeScatteringAngleswithMultipleGrids
