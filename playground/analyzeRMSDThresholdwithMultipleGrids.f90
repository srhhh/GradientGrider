!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       MODULE
!               analyzeRMSDThresholdwithMultipleGrids
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This module simply plots the percentage of frames of a trajectory that are
!               below a certain threshold; these DAT files are prepared beforehand
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILECHANNELS                    ACTION
!
!		FILECHANNEL1			OPEN, WRITE, CLOSE
!		FILEHCANNEL2			OPEN, READ, CLOSE
!               GNUPLOTCHANNEL                  OPEN, WRITE, CLOSE
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
!		gridpath0//percentrmsd//	DAT			Temporary buffer to accumulate a minimum rmsd
!			Ngrid						across multiple grids
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module analyzeRMSDThresholdwithMultipleGrids
implicit none

contains


subroutine getRMSDThresholds1(prefix_filename,JPGfilename)
use PARAMETERS
use ANALYSIS
implicit none

!NUMBER OF TRAJECTORIES CHECKED
integer :: n_testtraj

!INTERNAL TALLY FOR RMSD BELOW A THRESHOLD
real, dimension(Ntesttraj) :: percent_threshold_rmsd
integer :: total_threshold_rmsd, frames

!ARRAY HOLDING DATA FROM FILE
real(dp) :: min_rmsd

!FORMAT OF JPG FILES TO BE MADE
character(*), intent(in) :: JPGfilename
character(*), intent(in) :: prefix_filename

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectory_text_length) :: Ntraj_text

!I/O HANDLING
integer :: iostate

!INTEGER INCREMENTALS
integer :: n



do Ngrid = 1, Ngrid_total
write(variable_length_text,"(I5)") Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid

!We will bin data by GRID, not by trajectory
!So we uniquely name each output .dat and graph by the grid number
open(filechannel1,file=gridpath0//"percent_rmsd"//Ngrid_text//".dat")
do n_testtraj = 1, Ntesttraj
        write(variable_length_text,"(I5)") trajectory_text_length
        write(Ntraj_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") n_testtraj

        !Read the trajectory (which has the rmsd) across all grids line-by-line
        frames = 0
        total_threshold_rmsd = 0
        open(filechannel2,file=gridpath0//Ngrid_text//"/"//prefix_filename//&
                               "_"//Ntraj_text//".dat")

        do
                read(filechannel2,FMT=*,iostat=iostate) min_rmsd
                if (iostate /= 0) exit
                frames = frames + 1

        !If the RMSD is below the threshhold we tally that
                if (min_rmsd < threshold_RMSD) total_threshold_rmsd = total_threshold_rmsd + 1
        end do
        close(filechannel2)

        !We want the percentage of frames that has an RMSD below the threshhold
        !So we keep track of the number of frames and divide by that
        percent_threshold_rmsd(n_testtraj) = total_threshold_rmsd * 100.0 / frames
        write(filechannel1,FMT="(I6,1x,F7.3,1x,I8)") n_testtraj, percent_threshold_rmsd(n_testtraj), frames
end do
close(filechannel1)

end do

!Finally, plot the data
open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath0//JPGfilename//'.jpg"'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
write(variable_length_text,"(I5)") Ngrid_text_length
write(Ngrid_text,FMT="(I"//trim(adjustl(variable_length_text))//")") Ngrid_total
write(gnuplotchannel,*) 'set multiplot layout '//trim(adjustl(Ngrid_text))//&
                        ',2 columnsfirst margins 0.1,0.95,.1,.9 spacing 0.1,0'&
                        //' title "Trajectory RMSD Distribution with '//prefix_filename//' method'
write(gnuplotchannel,*) 'set title "Percentages of Trajectories with RMSD Below Threshold"'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'unset xtics'
write(gnuplotchannel,*) 'unset xlabel'
write(gnuplotchannel,*) 'ymax = 100.0'
write(variable_length_text,"(I5)") RMSD_Nbins
write(gnuplotchannel,*) 'Nbins = '//trim(adjustl(variable_length_text))
write(gnuplotchannel,*) 'bin_width = ymax / Nbins'
write(gnuplotchannel,*) 'set boxwidth bin_width'
write(gnuplotchannel,*) 'min(x,y) = (x < y) ? x : y'
write(gnuplotchannel,*) 'bin_number(x) = min(floor(x/bin_width),Nbins-1)'
write(gnuplotchannel,*) 'rounded(x) = bin_width * (bin_number(x) + 0.5)'
write(gnuplotchannel,*) 'set xrange [0:ymax]'
write(gnuplotchannel,*) 'set ylabel "Occurence"'

do Ngrid = 1, Ngrid_total
if (Ngrid == Ngrid_total) then
	write(gnuplotchannel,*) 'set xtics'
	write(gnuplotchannel,*) 'set xlabel "Percentage of Frames with RMSD Below Threshold"'
end if
write(variable_length_text,"(I5)") Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
write(gnuplotchannel,*) 'plot "'//gridpath0//'percent_rmsd'//Ngrid_text//'.dat'//&
                        '" u (rounded($2)):(1.0) smooth frequency with boxes'
write(gnuplotchannel,*) 'unset title'
end do

write(gnuplotchannel,*) 'set title "Distribution of Trajectories, Percent RMSD vs. Trajectory Length"'
write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'unset xtics'
write(gnuplotchannel,*) 'unset xlabel'
write(gnuplotchannel,*) 'set xrange [0:ymax]'
write(gnuplotchannel,*) 'set ylabel "Length of Trajectory (Frames)"'

write(variable_length_text,"(I5)") Ngrid_text_length
do Ngrid = 1, Ngrid_total
if (Ngrid == Ngrid_total) then
	write(gnuplotchannel,*) 'set xtics'
	write(gnuplotchannel,*) 'set xlabel "Percentage of Frames with RMSD Below Threshold"'
end if
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
write(gnuplotchannel,*) 'plot "'//gridpath0//'percent_rmsd'//Ngrid_text//'.dat'//&
                        '" u 2:3 with points'
write(gnuplotchannel,*) 'unset title'
end do

close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)

write(variable_length_text,"(I5)") Ngrid_text_length
do Ngrid = 1, Ngrid_total
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
call system("rm "//gridpath0//"percent_rmsd"//Ngrid_text//".dat")
end do



end subroutine getRMSDThresholds1

end module analyzeRMSDThresholdwithMultipleGrids
