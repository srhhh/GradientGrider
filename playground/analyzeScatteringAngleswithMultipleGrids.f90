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


subroutine getScatteringAngles1(DATfilename,scattering_angle_column,JPGfilename)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

!RANDOM TRAJECTORY SAMPLING
integer :: Nsamples

!SCATTERING ANGLE BINNING
integer :: Nsamples_max
integer :: scatteringAngle
real :: bin_width, binMean, binSD
real,allocatable :: binAverage(:,:),binRMSD(:)
integer,allocatable :: binTotal(:)
integer :: Nbins, binTally
real :: binThreshold = 0.1

!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
character(*), intent(in) :: DATfilename

!COLUMN OF DAT FILE WITH SCATTERING ANGLES
integer, intent(in) :: scattering_angle_column
real,dimension(10) :: line_data

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

!Incremental Integers
integer :: i, j


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
        write(gnuplotchannel,*) 'set autoscale y'
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
Nbins = Ntesttraj / 5
bin_width = pi / Nbins

!We will get a running average distribution by adding Ntesttraj trajectories at a time
!If the running average converges then we stop
Nsamples_max = (Ngrid_total * Ntraj_max) / Ntesttraj - 1
allocate(binAverage(Nbins,Nsamples_max+1),binTotal(Nbins,Nsamples_max+1),binRMSD(Nsamples_max))

binAverage(:,1) = 0.0
binTotal(:,:) = 0
Nsamples = 0
binTally = 0
open(filechannel1,file=gridpath0//Ngrid_text//"/"//cumulativefile//trim(Adjustl(Ntraj_text))".dat",action="read")
open(filechannel2,file=gridpath0//"RMSD"//cumulativefile//".dat",action="write")
do
        Nsamples = Nsamples + 1

        do i = 1, Ntesttraj
                read(filechannel1,FMT=*) line_data
                scatteringAngle = ceiling(line_data(scattering_angle_column) / bin_width)
                if (scatteringAngle > Nbins) scatteringAngle = Nbins
                binTotal(scatteringAngle,Nsamples+1) = binTotal(scatteringAngle,Nsamples+1) + 1
        end do

        binAverage(:,Nsamples+1) = (binAverage(:,Nsamples) * Nsamples + binTotal(:,Nsamples+1)) / (Nsamples)
        binRMSD(Nsamples) = sqrt(sum((binAverage(:,Nsamples+1)-binAverage(:,Nsamples))**2))/Nbins
        write(filechannel2,FMT=*) Nsamples*Ntesttraj, binRMSD(Nsamples), binThreshold

        if (binRMSD(Nsamples) < binThreshold) then
                binTally = binTally + 1
        else
                binTally = 0
        end if

        if (binTally == 5) exit
        if (Nsamples == Nsamples_max) then
                print *, "    No convergence for true scattering angle"
                print *, ""
                exit
        end if
end do
close(filechannel1)
close(filechannel2)

!Now we must do the laborious job of binning these
!Each bin has its own average and standard deviation based on how many samples we took
open(filechannel1,file=gridpath0//"Adjusted"//cumulativefile//".dat")
do i = 1, Nbins
        binMean = sum(binTotal(i,2:Nsamples+1)) * 1.0 / Nsamples
        binSD = sqrt(sum(((binTotal(i,2:Nsamples+1)-binMean)**2))/(Nsamples - 1))

        write(filechannel1,*) (i-0.5)*bin_width, binMean, binSD!/sqrt(real(Nsamples))
end do
close(filechannel1)

deallocate(binAverage,binTotal,binRMSD)

!We have everything we need to draw the distribution with error bars
open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) "set terminal jpeg size 1200,1200"
write(variable_length_text,FMT="(I5)") Ntesttraj*Nsamples
write(gnuplotchannel,*) 'set output "'//gridpath0//'Adjustedto'//&
                        trim(adjustl(variable_length_text))//JPGfilename//'"'
write(gnuplotchannel,*) 'set multiplot layout 2,1'
write(variable_length_text,FMT="(I5)") Ntesttraj
write(gnuplotchannel,*) 'set title "Convergence of the '//trim(adjustl(variable_length_text))//' Scattering Angle Distribution"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set xlabel "Number of Trajectories"'
write(gnuplotchannel,*) 'set ylabel "Distribution RMSD with Running Average"'
write(gnuplotchannel,*) 'set autoscale x'
write(variable_length_text,FMT="(I5)") Ntesttraj
write(gnuplotchannel,*) 'set xtics '//trim(adjustl(variable_length_text))
write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set label 1 "Convergence Threshold" at graph 0.02,0.08'
write(gnuplotchannel,*) 'plot "'//gridpath0//'RMSD'//cumulativefile//'.dat" u 1:2 w lines, \'
write(gnuplotchannel,*) '     "'//gridpath0//'RMSD'//cumulativefile//'.dat" u 1:3 w lines lc -1'
write(variable_length_text,FMT="(I5)") Ntesttraj
write(gnuplotchannel,*) 'set title "Final Scattering Angle Distribution for '//&
                        trim(adjustl(variable_length_text))//' Trajectories"'
write(gnuplotchannel,*) 'set xlabel "Scattering Angle (rad)"'
write(gnuplotchannel,*) 'set xrange [0:3.14159]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set ylabel "Frequency"'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'unset label 1'
write(gnuplotchannel,*) 'plot "'//gridpath0//'Adjusted'//cumulativefile//'.dat" u 1:2 w boxes, \'
write(gnuplotchannel,*) '     "'//gridpath0//'Adjusted'//cumulativefile//'.dat" u 1:2:3 w yerrorbars'
close(gnuplotchannel)

call system("gnuplot < "//gridpath0//gnuplotfile)

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
