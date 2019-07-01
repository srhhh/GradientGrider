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
!               This module plots the percentage of frames of a trajectory that are
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
!               gridpath1//prefix_filename//    DAT                     Stores the RMSD of the frame accepted at
!                       _#traj                                          each timestep of some trajectory coresponding
!                                                                       to some prefix
!		gridpath0//percent_rmsd//	DAT			Temporary buffer to accumulate a minimum rmsd
!			Ngrid						across multiple grids
!               gridpath0//JPGfilename          JPG                     Generic format for image names
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module analyzeRMSDThresholdwithMultipleGrids
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               getRMSDThresholds1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine looks through the RMSD of the accepted frames over all trajectories
!               corresponding to prefix_filename, compares them to some threshold, and plots the
!               percentage of frames below the threshold as a distribution per number of grids used
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               prefix_filename                 CHARACTER(*)                    The prefix defining a set of trajectories
!               JPGfilename                     CHARACTER(*)                    The suffix we use for the output JPG
!               selection                       INTEGER                         An integer describing which grids to plot
!                                                                               by representing "plot" as a 1 and "no
!                                                                               plot" as a 0 in binary
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               Ntesttraj                       INTEGER                         The number of trajectories in this
!                                                                               set of trajectories
!               Ngrid_total                     INTEGER                         The number of grids used in this
!                                                                               set of trajectories
!               Ngrid_plotting                  INTEGER                         The number of grids to be plotted
!               grid_selection                  INTEGER,DIM(Ngrid_max)          An array with indexes of grids to plot
!               RMSD_Nbins                      INTEGER                         The number of bins for this distribution
!
!               threshold_rmsd                  REAL                            The RMSD threshold
!               frames                          INTEGER                         The number of frames read so far
!               total_threshold_rmsd            INTEGER                         The number of frames whose accepted frame
!                                                                               has an RMSD under the threshold so far
!               percent_threshold_rmsd          REAL                            The percentage of frames in a set of
!                                                                               trajectories whose accepted frame has an
!                                                                               RMSD under the threshold
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath0//prefix_filename//    DAT                     Stores the RMSD of the frame accepted at
!                       Ngrid_#traj                                     each timestep of some trajectory corresponding
!                                                                       to some prefix
!		gridpath5//percent_rmsd//	DAT			Temporary buffer to accumulate a minimum rmsd
!			Ngrid						across multiple grids
!               gridpath4//JPGfilename          PNG                     Generic format for image names
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine getRMSDThresholds1(prefix_filename,JPGfilename,selection)
use PARAMETERS
use ANALYSIS
implicit none

!NUMBER OF TRAJECTORIES CHECKED
integer :: n_testtraj

!NUMBER OF GRIDS TO BE PLOTTED
integer :: Ngrid_plotting

!INTERNAL TALLY FOR RMSD BELOW A THRESHOLD
real :: percent_threshold_rmsd
integer :: total_threshold_rmsd, frames

!ARRAY HOLDING DATA FROM FILE
real(dp) :: min_rmsd, min_rmsd_prime

!ARRAY HOLDING SELECTION OF GRIDS TO PLOT
integer,optional,intent(in) :: selection
integer,dimension(Ngrid_max) :: grid_selection
integer :: selection_current

!FORMAT OF JPG FILES TO BE MADE
character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: JPGfilename
character(*), intent(in) :: prefix_filename

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectory_text_length) :: Ntraj_text

integer :: i1,i2,i3,i4
real(dp) :: r1,r2,r3,r4
integer :: count1, count2, total1, total2

!I/O HANDLING
integer :: iostate

!INTEGER INCREMENTALS
integer :: n

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

!If we only want to plot a select number of grids,
!then that can be passed in by the user as the
!variable "selection"
if (present(selection)) then
    Ngrid_plotting = 0
    selection_current = selection
    do Ngrid = 1, Ngrid_max
        selection_current = selection_current&
                - 2 ** (Ngrid - 1)

        if (selection_current < 0) exit

        !If the binary representation of the
        !selection has a 1 in that digits place
        !then that is a signal by the user that
        !the grid should be plotted
        if (modulo(selection_current,2**Ngrid) == 0) then
            Ngrid_plotting = Ngrid_plotting + 1
            grid_selection(Ngrid_plotting) = Ngrid
        else
            selection_current = selection_current&
                    + 2 ** (Ngrid - 1)
        end if
    end do

!If not supplied, then we will just plot all
!grids up until Ngrid_total
else
    Ngrid_plotting = Ngrid_total
    do Ngrid = 1, Ngrid_plotting
        grid_selection(Ngrid) = Ngrid
    end do
end if


!Now we need to open up each grid's set of trajectories
do Ngrid = 1, Ngrid_plotting
write(variable_length_text,"(I5)")&
        Ngrid_text_length
write(Ngrid_text,FMT="(I0."//&
        trim(adjustl(variable_length_text))//")")&
        grid_selection(Ngrid)

!We will bin data by GRID, not by trajectory
!So we uniquely name each output .dat and graph by the
!grid number
open(filechannel1,file=gridpath5//&
        "percent_rmsd"//Ngrid_text//".dat")
do n_testtraj = 1, Ntesttraj
    write(variable_length_text,"(I5)")&
            trajectory_text_length
    write(Ntraj_text,FMT="(I0."//&
            trim(adjustl(variable_length_text))//")")&
            n_testtraj

    !Read the trajectory (which has the rmsd) across
    !all grids line-by-line
    frames = 0
    total_threshold_rmsd = 0
    open(filechannel2,file=gridpath0//&
            prefix_filename//Ngrid_text//&
            "_"//Ntraj_text//".dat")

    do
        !Tally how many frames we have in the trajectory
        read(filechannel2,FMT=*,iostat=iostate) min_rmsd
        if (iostate /= 0) exit
        frames = frames + 1

        !And if the RMSD is below the threshhold we
        !tally that separately
        if (.not.(accept_worst) .and. &
                 (min_rmsd < threshold_RMSD))&
                 total_threshold_rmsd = &
                 total_threshold_rmsd + 1
        if ((accept_worst) .and. &
                (min_rmsd > 0.0d0))&
                 total_threshold_rmsd = &
                 total_threshold_rmsd + 1
    end do
    close(filechannel2)

    !We want the percentage of frames that has an
    !RMSD below the threshhold so we keep track of
    !the number of frames and divide by that
    percent_threshold_rmsd = &
            total_threshold_rmsd * 100.0 / frames
    write(filechannel1,FMT="(I6,1x,F7.3,1x,I8)")&
            n_testtraj, percent_threshold_rmsd, frames
end do
close(filechannel1)

end do

!Finally, plot the data
open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath4//JPGfilename//'.jpg"'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
write(variable_length_text,"(I5)") Ngrid_text_length
write(Ngrid_text,FMT="(I"//trim(adjustl(variable_length_text))//")") Ngrid_plotting
write(gnuplotchannel,*) 'set multiplot layout '//trim(adjustl(Ngrid_text))//&
                        ',2 columnsfirst margins 0.1,0.95,.1,.9 spacing 0.1,0'&
                        //' title "Trajectory RMSD Distribution with '//expfolder//' method"'
write(gnuplotchannel,*) 'set title "Percentages of Trajectories with RMSD Below Threshold"'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'scaling = 1'
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
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'unset ylabel'

do Ngrid = 1, Ngrid_plotting
if (Ngrid == Ngrid_plotting) then
    write(gnuplotchannel,*) 'set label 1 "Occurence" at screen 0.01,0.45 rotate by 90'
    write(gnuplotchannel,*) 'set xtics'
    write(gnuplotchannel,*) 'set xlabel "Percentage of Frames with RMSD Below Threshold"'
end if

write(variable_length_text,"(I5)") Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") grid_selection(Ngrid)
write(gnuplotchannel,*) 'plot "'//gridpath5//'percent_rmsd'//Ngrid_text//'.dat'//&
                        '" u (rounded($2)):(1.0/scaling) smooth frequency with boxes'
write(gnuplotchannel,*) 'unset title'
end do

write(gnuplotchannel,*) 'set title "Distribution of Trajectories, Percent RMSD vs. Trajectory Length"'
write(gnuplotchannel,*) 'scaling = 1000'
write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'unset xtics'
write(gnuplotchannel,*) 'unset xlabel'
write(gnuplotchannel,*) 'set xrange [0:ymax]'
write(gnuplotchannel,*) 'unset ylabel'

write(variable_length_text,"(I5)") Ngrid_text_length
do Ngrid = 1, Ngrid_plotting
if (Ngrid == Ngrid_plotting) then
    write(gnuplotchannel,*) 'set label 2 "Length of Trajectory (Thousands of Frames)" at screen 0.51,0.40 rotate by 90'
    write(gnuplotchannel,*) 'set xtics'
    write(gnuplotchannel,*) 'set xlabel "Percentage of Frames with RMSD Below Threshold"'
end if

write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") grid_selection(Ngrid)
write(gnuplotchannel,*) 'plot "'//gridpath5//'percent_rmsd'//Ngrid_text//'.dat'//&
                        '" u 2:(($3)/scaling) with points'
write(gnuplotchannel,*) 'unset title'

end do

close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

write(variable_length_text,"(I5)") Ngrid_text_length
do Ngrid = 1, Ngrid_plotting
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") grid_selection(Ngrid)
call system("rm "//gridpath5//"percent_rmsd"//Ngrid_text//".dat")
end do


!frames = 0
!
!write(variable_length_text,"(I5)") Ngrid_text_length
!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") 1
!open(filechannel1,file=gridpath0//"percent_rmsd"//Ngrid_text//".dat")
!
!open(filechannel2,file=gridpath0//checkstatefile)
!!read(filechannel1) number_of_frames,order,neighbor_check,steps,&
!!                   min_rmsd,min_rmsd_prime,vals(1),vals(2),U,KE
!do 
!read(filechannel2,FMT=*,iostat=iostate) i1, i2, i3, i4,&
!                   min_rmsd,min_rmsd_prime,r1,r2,r3,r4
!if (iostate /= 0) exit
!write(filechannel1,FMT="(I6,1x,F9.6,1x,I8)") i4, max(min_rmsd,.00001), i1
!frames = frames + 1
!end do
!close(filechannel2)
!
!close(filechannel1)
!
!write(variable_length_text,"(I5)") Ngrid_text_length
!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") 2
!open(filechannel1,file=gridpath0//"percent_rmsd"//Ngrid_text//".dat")
!
!open(filechannel2,file=gridpath0//checkstatefile)
!!read(filechannel1) number_of_frames,order,neighbor_check,steps,&
!!                   min_rmsd,min_rmsd_prime,vals(1),vals(2),U,KE
!do 
!read(filechannel2,FMT=*,iostat=iostate) i1, i2, i3, i4,&
!                   min_rmsd,min_rmsd_prime,r1,r2,r3,r4
!if (iostate /= 0) exit
!write(filechannel1,FMT="(I6,1x,F9.6,1x,I8)") i4, max(min_rmsd_prime,.00001), i1
!end do
!close(filechannel2)
!
!close(filechannel1)
!
!total1 = 0
!total2 = 0
!
!write(variable_length_text,"(I5)") Ngrid_text_length
!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") 3
!open(filechannel1,file=gridpath0//"percent_rmsd"//Ngrid_text//".dat")
!
!open(filechannel2,file=gridpath0//checkstatefile)
!!read(filechannel1) number_of_frames,order,neighbor_check,steps,&
!!                   min_rmsd,min_rmsd_prime,vals(1),vals(2),U,KE
!do 
!read(filechannel2,FMT=*,iostat=iostate) i1, i2, i3, i4,&
!                   min_rmsd,min_rmsd_prime,r1,r2,r3,r4
!if (iostate /= 0) exit
!write(filechannel1,FMT="(I6,1x,F9.6,1x,F9.6,1x,I8)") i4, max(abs(min_rmsd - min_rmsd_prime),0.00001), &
!                                             min_rmsd - min_rmsd_prime, i1
!if ((min_rmsd-min_rmsd_prime) == 0.0d0) then
!else if ((min_rmsd-min_rmsd_prime) < 0.0d0) then
!        total1 = total1 + 1
!        total2 = total2 + 1
!else
!        total1 = total1 + 1
!end if
!end do
!close(filechannel2)
!
!close(filechannel1)
!
!print *, "total1 (number of nonzero RMSD differences): ", total1
!print *, "total2 (number of negative RMSD differences): ", total2
!
!
!write(variable_length_text,"(I5)") Ngrid_text_length
!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") 4
!open(filechannel1,file=gridpath0//"percent_rmsd"//Ngrid_text//".dat")
!
!open(filechannel2,file=gridpath0//checkstatefile)
!!read(filechannel1) number_of_frames,order,neighbor_check,steps,&
!!                   min_rmsd,min_rmsd_prime,vals(1),vals(2),U,KE
!do 
!read(filechannel2,FMT=*,iostat=iostate) i1, i2, i3, i4,&
!                   min_rmsd,min_rmsd_prime,r1,r2,r3,r4
!if (iostate /= 0) exit
!write(filechannel1,FMT="(I6,1x,F9.6,1x,I8)") i4, max(min_rmsd_prime/min(min_rmsd,1.0),.000001), i1
!end do
!close(filechannel2)
!
!close(filechannel1)
!
!
!
!
!Ngrid_plotting = 4
!
!!Finally, plot the data
!open(gnuplotchannel,file=gridpath0//gnuplotfile)
!write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
!write(gnuplotchannel,*) 'set output "'//gridpath0//JPGfilename//'.jpg"'
!write(gnuplotchannel,*) 'set tmargin 0'
!write(gnuplotchannel,*) 'set bmargin 0'
!write(gnuplotchannel,*) 'set lmargin 1'
!write(gnuplotchannel,*) 'set rmargin 1'
!write(variable_length_text,"(I5)") Ngrid_text_length
!write(Ngrid_text,FMT="(I"//trim(adjustl(variable_length_text))//")") Ngrid_plotting
!write(gnuplotchannel,*) 'set multiplot layout '//trim(adjustl(Ngrid_text))//&
!                        ',2 columnsfirst margins 0.1,0.95,.1,.9 spacing 0.1,0'&
!                        //' title "Trajectory RMSD Distribution with '//prefix_filename//' method"'
!write(gnuplotchannel,*) 'set title "RMSD Histogram Over the Trajectory"'
!write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
!write(gnuplotchannel,*) 'scaling = ', frames
!write(gnuplotchannel,*) 'unset key'
!write(gnuplotchannel,*) 'unset xtics'
!write(gnuplotchannel,*) 'unset xlabel'
!!write(gnuplotchannel,*) 'ymax = 100.0'
!write(gnuplotchannel,*) 'ymax =', 0.00001
!write(variable_length_text,"(I5)") RMSD_Nbins
!write(gnuplotchannel,*) 'Nbins = '//trim(adjustl(variable_length_text))
!!write(gnuplotchannel,*) 'bin_width = ymax / Nbins'
!write(gnuplotchannel,*) 'bin_width = (2.0 - log10(ymax)) / Nbins'
!write(gnuplotchannel,*) 'set boxwidth bin_width'
!write(gnuplotchannel,*) 'min(x,y) = (x < y) ? x : y'
!!write(gnuplotchannel,*) 'bin_number(x) = min(floor(x/bin_width),Nbins-1)'
!write(gnuplotchannel,*) 'bin_number(x) = min(floor(log10(x)/bin_width),Nbins-1)'
!write(gnuplotchannel,*) 'rounded(x) = bin_width * (bin_number(x) + 0.5)'
!!write(gnuplotchannel,*) 'set xrange [0:ymax]'
!write(gnuplotchannel,*) 'set xrange [log10(ymax):2.0]'
!write(gnuplotchannel,*) 'set yrange [0:0.2]'
!write(gnuplotchannel,*) 'unset ylabel'
!
!do Ngrid = 1, Ngrid_plotting
!if (Ngrid == Ngrid_plotting) then
!        write(gnuplotchannel,*) 'set label 1 "Occurence" at screen 0.01,0.45 rotate by 90'
!	write(gnuplotchannel,*) 'set xtics'
!	write(gnuplotchannel,*) 'set xlabel "log(RMSD)"'
!end if
!
!if (Ngrid == 3) then
!        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
!        write(gnuplotchannel,*) 'plot "'//gridpath0//'percent_rmsd'//Ngrid_text//'.dat'//&
!                                '" u (rounded($2)):(1.0/scaling) lc rgb "blue" '//&
!                                'smooth frequency with boxes,\'
!        write(gnuplotchannel,*) '     "'//gridpath0//'percent_rmsd'//Ngrid_text//'.dat'//&
!                                '" u ($3>0?(rounded($3)):1/0):(1.0/scaling) lc rgb "red" '//&
!                                'smooth frequency with boxes'
!        write(gnuplotchannel,*) 'unset title'
!else
!        write(variable_length_text,"(I5)") Ngrid_text_length
!        !write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") grid_selection(Ngrid)
!        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
!        write(gnuplotchannel,*) 'plot "'//gridpath0//'percent_rmsd'//Ngrid_text//'.dat'//&
!                                '" u (rounded($2)):(1.0/scaling) smooth frequency with boxes'
!        write(gnuplotchannel,*) 'unset title'
!end if
!end do
!
!write(gnuplotchannel,*) 'set title "RMSD Distribution vs. The Number of Frames Checked"'
!write(gnuplotchannel,*) 'scaling = 100'
!write(gnuplotchannel,*) 'set autoscale y'
!write(gnuplotchannel,*) 'unset xtics'
!write(gnuplotchannel,*) 'unset xlabel'
!!write(gnuplotchannel,*) 'set xrange [0:ymax]'
!write(gnuplotchannel,*) 'set xrange [log10(ymax):2.0]'
!write(gnuplotchannel,*) 'unset ylabel'
!
!write(variable_length_text,"(I5)") Ngrid_text_length
!do Ngrid = 1, Ngrid_plotting
!if (Ngrid == Ngrid_plotting) then
!        write(gnuplotchannel,*) 'set label 2 "Number of Frames Checked (Hundreds of Frames)" at screen 0.51,0.40 rotate by 90'
!	write(gnuplotchannel,*) 'set xtics'
!	write(gnuplotchannel,*) 'set xlabel "log(RMSD)"'
!end if
!!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") grid_selection(Ngrid)
!!write(gnuplotchannel,*) 'plot "'//gridpath0//'percent_rmsd'//Ngrid_text//'.dat'//&
!!                        '" u 2:(($3)/scaling) with points'
!
!if (Ngrid == 3) then
!        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
!        write(gnuplotchannel,*) 'plot "'//gridpath0//'percent_rmsd'//Ngrid_text//'.dat'//&
!                                '" u (log10($2)):(($4)/scaling) lc rgb "blue",\'
!        write(gnuplotchannel,*) '     "'//gridpath0//'percent_rmsd'//Ngrid_text//'.dat'//&
!                                '" u ($3>0?(log10($3)):1/0):(($4)/scaling) lc rgb "red"'
!        write(gnuplotchannel,*) 'unset title'
!else
!        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
!        write(gnuplotchannel,*) 'plot "'//gridpath0//'percent_rmsd'//Ngrid_text//'.dat'//&
!                                '" u (log10($2)):(($3)/scaling) with points'
!        write(gnuplotchannel,*) 'unset title'
!end if
!end do
!
!close(gnuplotchannel)
!
!call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)
!
!write(variable_length_text,"(I5)") Ngrid_text_length
!do Ngrid = 1, Ngrid_plotting
!!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") grid_selection(Ngrid)
!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
!!call system("rm "//gridpath0//"percent_rmsd"//Ngrid_text//".dat")
!end do


end subroutine getRMSDThresholds1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               getRMSDDifferences1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine looks at the RMSD acceptance rates between two trials, particularly
!               for different subcell neighbor checks (subcellsearch_max)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               PNGfilename                     CHARACTER(*)                    The suffix we use for the output PNG
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               Nbins                           INTEGER                         The number of bins for this distribution
!               threshold_rmsd                  REAL                            The RMSD threshold
!               frames                          INTEGER                         The number of frames read so far
!               number_of_frames_accepted       INTEGER                         The number of frames whose accepted frame has
!                                                                               an RMSD under the threshold so far
!               min_rmsd                        REAL                            The RMSD of the accepted frame (trial 1)
!               min_rmsd_prime                  REAL                            The RMSD of the accepted frame (trial 2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath5//checkstatefile       DAT                     Stores the RMSD of the accepted frames of
!                                                                       each trial (and other data)
!		gridpath5//percent_rmsd//	DAT			Temporary buffer to accumulate a minimum rmsd
!			Ngrid						across multiple grids
!               gridpath4//PNGfilename          PNG                     Generic format for image names
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getRMSDDifferences1(PNGfilename)
use PARAMETERS
use ANALYSIS
implicit none

!NUMBER OF FRAMES IN THE TRAJECTORY
integer :: frames

!RMSD
real(dp) :: min_rmsd, min_rmsd_prime
real(dp) :: min_min_rmsd
integer :: number_of_frames_accepted

!FORMAT OF PNG FILES TO BE MADE
character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: PNGfilename

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text

!OTHER VARIABLES IN THE CHECKSTATE FILE
integer :: i1,i2,i3,i4
real(dp) :: r1,r2,r3,r4

!I/O HANDLING
integer :: iostate

!HISTOGRAM VARIABLES
integer :: Nbins
real(dp) :: bin_width

!INTEGER INCREMENTALS
integer :: n

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

frames = 0
number_of_frames_accepted = 0
min_min_rmsd = default_rmsd

open(filechannel2,file=gridpath5//checkstatefile)
!read(filechannel1) number_of_frames,order,i&
!                   neighbor_check,steps,&
!                   min_rmsd,min_rmsd_prime,&
!                   vals(1),vals(2),U,KE
do 
    read(filechannel2,FMT=*,iostat=iostate)&
            i1, i2, i3, i4,&
            min_rmsd,min_rmsd_prime,&
            r1,r2,r3,r4
    if (iostate /= 0) exit

    min_min_rmsd = min(min_min_rmsd,&
            min_rmsd,min_rmsd_prime)
    if (min_rmsd_prime < threshold_rmsd) &
        number_of_frames_accepted =&
        number_of_frames_accepted + 1
    frames = frames + 1
end do
close(filechannel2)


Nbins = 50
bin_width = 2 * log10(default_rmsd/&
        min_min_rmsd) / Nbins

write(variable_length_text,"(I5)")&
        Ngrid_text_length
write(Ngrid_text,FMT="(I0."//&
        trim(adjustl(variable_length_text))//")")&
        4
open(filechannel1,file=gridpath5//&
        "percent_rmsd"//Ngrid_text//".dat")

open(filechannel2,file=gridpath5//checkstatefile)
!read(filechannel1) number_of_frames,order,&
!                   neighbor_check,steps,&
!                   min_rmsd,min_rmsd_prime,&
!                   vals(1),vals(2),U,KE
do 
    read(filechannel2,FMT=*,iostat=iostate)&
            i1, i2, i3, i4,&
            min_rmsd,min_rmsd_prime,&
            r1,r2,r3,r4
    if (iostate /= 0) exit

    write(filechannel1,FMT="(I6,1x,F9.6,1x,"//&
            "F9.6,1x,I8,1x,I8,1x,I8)") &
            i4, min_rmsd, min_rmsd_prime, i2, &
            nint(log10((min_rmsd_prime)/&
            (min_min_rmsd/default_rmsd)) / (0.5*bin_width)),&
            nint(log10((min_rmsd/min_rmsd_prime)/&
            (min_min_rmsd/default_rmsd)) / bin_width)
end do
close(filechannel2)
close(filechannel1)



open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//PNGfilename//'.png"'
write(gnuplotchannel,*) 'set title "RMSD Encountered in Binary Check"'
write(gnuplotchannel,*) 'set key left top'
write(gnuplotchannel,*) 'set logscale x'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,*) 'min_min_rmsd = ', min_min_rmsd
write(gnuplotchannel,*) 'set xrange [min_min_rmsd:',default_rmsd,']'
write(gnuplotchannel,*) 'set yrange [min_min_rmsd:',default_rmsd,']'
write(gnuplotchannel,*) 'set xlabel "RMSD (A) Encountered for the First Check"'
write(gnuplotchannel,*) 'set ylabel "RMSD (A) Encountered for the Second Check"'
write(gnuplotchannel,*) 'plot "'//gridpath5//"percent_rmsd"//Ngrid_text//'.dat" u '//&
                               '($4==0?$2:1/0):3 w p lc rgb "red" title "Both Order 0",\'
write(gnuplotchannel,*) '     "'//gridpath5//"percent_rmsd"//Ngrid_text//'.dat" u '//&
                               '($4==1?$2:1/0):3 w p lc rgb "blue" title "One Order 1",\'
write(gnuplotchannel,*) '     "'//gridpath5//"percent_rmsd"//Ngrid_text//'.dat" u '//&
                               '($4==2?$2:1/0):3 w p lc rgb "green" title "Both Order 1"'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//PNGfilename//'_1.png"'
write(gnuplotchannel,*) 'set title "Ratio of RMSD Encountered in Binary Check"'
write(gnuplotchannel,*) 'set key left top'
write(gnuplotchannel,*) 'scaling = ', frames
write(gnuplotchannel,*) 'xmin = ', min_min_rmsd / default_rmsd
write(gnuplotchannel,*) 'xmax = ', default_rmsd / min_min_rmsd
write(gnuplotchannel,*) 'set xlabel "RMSD (first check) / RMSD (second check)"'
write(gnuplotchannel,*) 'set ylabel "Occurence"'

write(gnuplotchannel,*) 'Nbins = ', Nbins
write(gnuplotchannel,*) 'bin_width = ', bin_width
write(gnuplotchannel,*) 'set boxwidth 1'
write(gnuplotchannel,*) 'min(x,y) = (x < y) ? x : y'
write(gnuplotchannel,*) 'set xrange [0:Nbins]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set style fill transparent solid 0.5'

write(gnuplotchannel,*) 'set xtics ('//&
                                           '"10^{-5}" (-5-log10(xmin))/bin_width, '//&
                                           '"10^{-4}" (-4-log10(xmin))/bin_width, '//&
                                           '"10^{-3}" (-3-log10(xmin))/bin_width, '//&
                                           '"10^{-2}" (-2-log10(xmin))/bin_width, '//&
                                           '"10^{-1}" (-1-log10(xmin))/bin_width, '//&
                                               '"10^0" (0-log10(xmin))/bin_width, '//&
                                               '"10^1" (1-log10(xmin))/bin_width, '//&
                                               '"10^2" (2-log10(xmin))/bin_width, '//&
                                               '"10^3" (3-log10(xmin))/bin_width, '//&
                                               '"10^4" (4-log10(xmin))/bin_width, '//&
                                               '"10^5" (5-log10(xmin))/bin_width)'

write(gnuplotchannel,*) 'plot "'//gridpath5//"percent_rmsd"//Ngrid_text//'.dat" u '//&
                               '6:($4==2?(1.0/scaling):0.0) '//&
                               'smooth frequency w boxes lc rgb "green" title "Both Order 1",\'
write(gnuplotchannel,*) '     "'//gridpath5//"percent_rmsd"//Ngrid_text//'.dat" u '//&
                               '6:($4==1?(1.0/scaling):0.0) '//&
                               'smooth frequency w boxes lc rgb "blue" title "One Order 1",\'
write(gnuplotchannel,*) '     "'//gridpath5//"percent_rmsd"//Ngrid_text//'.dat" u '//&
                               '6:($4==0?(1.0/scaling):0.0) '//&
                               'smooth frequency w boxes lc rgb "red" title "Both Order 0"'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//PNGfilename//'_2.png"'
write(gnuplotchannel,*) 'set title "Ratio of RMSD Encountered for Binary Check"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'scaling = ', frames
write(gnuplotchannel,*) 'xmin = ', min_min_rmsd
write(gnuplotchannel,*) 'xmax = ', default_rmsd
write(gnuplotchannel,*) 'set xlabel "RMSD (A) Encountered After the Second Check"'
write(gnuplotchannel,*) 'set ylabel "Occurence"'

write(variable_length_text,"(F5.2)") number_of_frames_accepted * 100.0 / frames
write(gnuplotchannel,*) 'set label 1 "Acceptance Rate: '//variable_length_text//'%" at graph 0.1,0.8'

write(gnuplotchannel,*) 'Nbins = ', Nbins
write(gnuplotchannel,*) 'bin_width = ', 0.5*bin_width
write(gnuplotchannel,*) 'set boxwidth 1'
write(gnuplotchannel,*) 'min(x,y) = (x < y) ? x : y'
write(gnuplotchannel,*) 'set xrange [-0.5:Nbins+.05]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set style fill transparent solid 0.5'

write(gnuplotchannel,*) 'set xtics ('//&
                                           '"10^{-5}" (-5-log10(xmin))/bin_width, '//&
                                           '"10^{-4}" (-4-log10(xmin))/bin_width, '//&
                                           '"10^{-3}" (-3-log10(xmin))/bin_width, '//&
                                           '"10^{-2}" (-2-log10(xmin))/bin_width, '//&
                                           '"10^{-1}" (-1-log10(xmin))/bin_width, '//&
                                               '"10^0" (0-log10(xmin))/bin_width, '//&
                                               '"10^1" (1-log10(xmin))/bin_width, '//&
                                               '"10^2" (2-log10(xmin))/bin_width, '//&
                                               '"10^3" (3-log10(xmin))/bin_width, '//&
                                               '"10^4" (4-log10(xmin))/bin_width, '//&
                                               '"10^5" (5-log10(xmin))/bin_width)'

write(gnuplotchannel,*) 'plot "'//gridpath5//"percent_rmsd"//Ngrid_text//'.dat" u '//&
                               '5:(1.0/scaling) '//&
                               'smooth frequency w boxes lc rgb "green"'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)



end subroutine getRMSDDifferences1


subroutine getAlphaErrorDistribution(PNGfilename)
use PARAMETERS
use ANALYSIS
implicit none

integer :: frames
integer,dimension(Nalpha) :: alpha_flagging
integer,dimension(Nalpha) :: alpha_binning
real(dp),dimension(Nalpha) :: alpha_array

!FORMAT OF PNG FILES TO BE MADE
character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: PNGfilename

!I/O HANDLING
integer :: iostate

!INTEGER INCREMENTALS
integer :: n

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

alpha_binning = 0
frames = 0

open(6666,file=gridpath5//alphafile)
do
    read(6666,FMT=*,iostat=iostate)&
            (alpha_flagging(n),n=1,Nalpha)

    if (iostate /= 0) exit

    frames = frames + 1
    alpha_binning = alpha_binning +&
            alpha_flagging

end do
close(6666)

do n = 1, Nalpha
    alpha_array(n) = alpha_start + (n-1-0.5) * &
        (alpha_end - alpha_start) / (Nalpha-1)
end do

if (logarithmic_alpha_flag) then
    do n = 1, Nalpha
        alpha_array(n) = &
             10.0d0 ** alpha_array(n)
    end do
end if

open(6666,file=gridpath5//"tmpA.dat")
do n = 1, Nalpha
    write(6666,FMT=*) alpha_array(n),&
            alpha_binning(n) * 1.0d2 / frames
end do
close(6666)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//PNGfilename//'.png"'
write(gnuplotchannel,*) 'set title "Blue Point Occurence by Alpha Ratio"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set yrange [0:', &
        min(100.0d0,maxval(alpha_binning)*110.0d0/frames),']'

if (logarithmic_alpha_flag) then
    write(gnuplotchannel,*) 'set logscale x'
    write(gnuplotchannel,FMT="(A,E16.6,':',E16.6,A)") &
            'set xrange [',&
            10.0d0 ** (alpha_start - (alpha_end-alpha_start)/Nalpha),&
            10.0d0 ** (alpha_end + (alpha_end-alpha_start)/Nalpha),']'
else
    write(gnuplotchannel,FMT="(A,E16.6,':',E16.6,A)") &
            'set xrange [',&
            (alpha_start - (alpha_end-alpha_start)/Nalpha),&
            (alpha_end + (alpha_end-alpha_start)/Nalpha),']'
end if


write(gnuplotchannel,*) 'set xlabel "Alpha Ratio"'
write(gnuplotchannel,*) 'set ylabel "Percentage of Frames with IE > AE"'
write(gnuplotchannel,*) 'set style histogram clustered gap 1'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'plot "'//gridpath5//'tmpA.dat" u 1:2 w boxes'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

end subroutine getAlphaErrorDistribution



subroutine getRMSDinterpolation(vals,delta_vals,PNGfilename)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!COLLECTIVE VARIABLES TO FILTER DATA
real(dp),dimension(Nvar),intent(in) :: vals, delta_vals
character(5) :: vals_interpolation_text

!NUMBER OF FRAMES IN THE DATA
real(dp),dimension(7) :: min_rmsd_vals, max_rmsd_vals
!real(dp),dimension(7) :: min_rsv, max_rsv, rsv

!RMSD
real(dp) :: rmsd_x, rmsd_y, rmsd_z, rmsd_fx
real(dp) :: min_rmsd_x, max_rmsd_x
real(dp) :: min_rmsd_y, max_rmsd_y
real(dp) :: min_rmsd_z, max_rmsd_z
real(dp) :: min_rmsd_fx, max_rmsd_fx

!FORMAT OF PNG FILES TO BE MADE
character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: PNGfilename

!FORMATTING OF PNG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text

!VARIABLES IN THE INTERPOLATION FILE
integer :: Ninterpolation
integer :: min_Ninterpolation, max_Ninterpolation
real :: vals1, vals2
!real(dp),dimension(2) :: temp_vals
real(dp),dimension(3) :: RMSDheatmap_coeff

!I/O HANDLING
integer :: iostate

!HISTOGRAM VARIABLES
integer :: Nbins

!INTEGER INCREMENTALS
integer :: n

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

!frames = 0
!tally1 = 0
!tally2 = 0
!
!min_rmsd_x = default_rmsd
!min_rmsd_y = default_rmsd
!min_rmsd_z = default_rmsd
!min_rmsd_fx = 1.0d9
!min_Ninterpolation = 1000
!max_rmsd_x = 0.0d0
!max_rmsd_y = 0.0d0
!max_rmsd_z = 0.0d0
!max_rmsd_fx = 0.0d0
!max_Ninterpolation = 0
!
!write(vals_interpolation_text,FMT="(F7.3,'_',F7.3)") vals(1),vals(2)
!open(filechannel1,file=gridpath0//vals_interpolation_text//interpolationfile)
!open(filechannel2,file=gridpath0//interpolationfile)
!do 
!        read(filechannel2,FMT=*,iostat=iostate) vals1, vals2, Ninterpolation, &
!                                       rmsd_y, rmsd_z, &
!                                       rmsd_x_prime, rmsd_fx_prime, &
!                                       rmsd_x, rmsd_fx
!        if (iostate /= 0) exit
!
!        if ((abs(vals1-vals(1)) > delta_vals(1)).or.&
!            (abs(vals2-vals(2)) > delta_vals(2))) cycle
!
!        min_rmsd_x = min(min_rmsd_x, rmsd_x)
!        max_rmsd_x = max(max_rmsd_x, rmsd_x)
!        min_rmsd_y = min(min_rmsd_y, rmsd_y)
!        max_rmsd_y = max(max_rmsd_y, rmsd_y)
!        min_rmsd_z = min(min_rmsd_z, rmsd_z)
!        max_rmsd_z = max(max_rmsd_z, rmsd_z)
!        min_rmsd_fx = min(min_rmsd_fx, rmsd_fx)
!        max_rmsd_fx = max(max_rmsd_fx, rmsd_fx)
!        min_Ninterpolation = min(min_Ninterpolation, Ninterpolation)
!        max_Ninterpolation = max(max_Ninterpolation, Ninterpolation)
!
!        frames = frames + 1
!!       if (Ninterpolation == 1) tally1 = tally1 + 1
!!       if (Ninterpolation == 2) tally2 = tally2 + 1
!
!        write(filechannel1,FMT=*) vals1, vals2, Ninterpolation,&
!               rmsd_y, rmsd_z, rmsd_x, rmsd_fx
!
!end do
!close(filechannel1)
!close(filechannel2)
!
!if (frames == 0) return
!
!Nbins = 100
!min_Nbin = 10
!max_Nbin = 0
!bin_width = log10((max_rmsd_fx/min_rmsd_x) / (min_rmsd_fx/max_rmsd_x)) / Nbins
!
!write(variable_length_text,"(I5)") Ngrid_text_length
!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") 4
!open(filechannel1,file=gridpath0//"percent_rmsd"//Ngrid_text//".dat")
!
!open(filechannel2,file=gridpath0//interpolationfile)
!do 
!        read(filechannel2,FMT=*,iostat=iostate) vals1, vals2, Ninterpolation, &
!                                       interpolation_alpha, threshold_rmsd_1, &
!                                       rmsd_x_prime, rmsd_fx_prime, &
!                                       rmsd_x, rmsd_fx
!        if (iostate /= 0) exit
!
!        if ((abs(vals1-vals(1)) > delta_vals(1)).or.&
!            (abs(vals2-vals(2)) > delta_vals(2))) cycle
!
!        Nbin = nint(log10((rmsd_fx/rmsd_x) / (min_rmsd_fx/max_rmsd_x)) / bin_width)
!
!        min_Nbin = min(min_Nbin,Nbin)
!        max_Nbin = max(max_Nbin,Nbin)
!
!        write(filechannel1,FMT=*) vals1, vals2, Nbin
!end do
!close(filechannel2)
!close(filechannel1)
!
!bin_width1 = log10(max_rmsd_z / min_rmsd_z) / Nbins
!bin_width2 = log10(max_rmsd_x / min_rmsd_x) / Nbins
!
!allocate(RMSDheatmap(Nbins,Nbins))
!RMSDheatmap = 1.0d-7
!
!write(variable_length_text,"(I5)") Ngrid_text_length
!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") 4
!open(filechannel2,file=gridpath0//interpolationfile)
!do 
!        read(filechannel2,FMT=*,iostat=iostate) vals1, vals2, Ninterpolation, &
!                                       rmsd_y, rmsd_z, &
!                                       rmsd_x_prime, rmsd_fx_prime, &
!                                       rmsd_x, rmsd_fx
!        if (iostate /= 0) exit
!
!        if ((abs(vals1-vals(1)) > delta_vals(1)).or.&
!            (abs(vals2-vals(2)) > delta_vals(2))) cycle
!
!        Nrmsd1 = nint(log10(rmsd_z/min_rmsd_z)/bin_width1)
!        if (Nrmsd1 < 1) Nrmsd1 = 1
!        Nrmsd2 = nint(log10(rmsd_x/min_rmsd_x)/bin_width2)
!        if (Nrmsd2 < 1) Nrmsd2 = 1
!
!        RMSDheatmap(Nrmsd1,Nrmsd2) = max(&
!                RMSDheatmap(Nrmsd1,Nrmsd2),rmsd_fx)
!end do
!close(filechannel2)
!
!Nheatmap = 0
!do Nrmsd1 = 1, Nbins
!        do Nrmsd2 = 1, Nbins
!                if (RMSDheatmap(Nrmsd1,Nrmsd2) > 1.0d-7) then
!                         Nheatmap = Nheatmap + 1
!                end if
!        end do
!end do
!
!allocate(A(Nheatmap,3),b(Nheatmap))
!Nheatmap = 0
!do Nrmsd1 = 1, Nbins
!        do Nrmsd2 = 1, Nbins
!                if (RMSDheatmap(Nrmsd1,Nrmsd2) > 1.0d-7) then
!                         Nheatmap = Nheatmap + 1
!                         A(Nheatmap,:) = (/ (Nrmsd1-0.5)*bin_width1,&
!                                 (Nrmsd2-0.5)*bin_width2,1.0d0/)
!                         b(Nheatmap) = RMSDheatmap(Nrmsd1,Nrmsd2)
!                end if
!        end do
!end do
!
!call LS(A,Nheatmap,3,b,RMSDheatmap_coeff)
!
!open(filechannel2,file=gridpath0//"heatmap_rmsd"//Ngrid_text//".dat")
!do Nrmsd1 = 1, Nbins
!        do Nrmsd2 = 1, Nbins
!                write(filechannel2,FMT=*) (Nrmsd1-0.5)*bin_width1,&
!                        (Nrmsd2-0.5)*bin_width2,RMSDheatmap(Nrmsd1,Nrmsd2),&
!                        max((Nrmsd1-0.5)*bin_width1*RMSDheatmap_coeff(1)+&
!                            (Nrmsd2-0.5)*bin_width2*RMSDheatmap_coeff(2)+&
!                                                    RMSDheatmap_coeff(3),1.0d-7)
!        end do
!        write(filechannel2,*) ""
!end do
!close(filechannel2)
!

Nbins = 100

call processInterpolationFile(vals,delta_vals,&
        RMSDheatmap_coeff,&
        min_rmsd_vals,max_rmsd_vals,&
        min_Ninterpolation,max_Ninterpolation)

write(vals_interpolation_text,FMT="(F7.3,'_',F7.3)") vals(1),vals(2)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,2400'
write(gnuplotchannel,*) 'set output "'//gridpath4//PNGfilename//'.png"'
write(gnuplotchannel,*) 'set title "RMSD Comparison of a Frame and Gradient with Interpolation"'
write(gnuplotchannel,*) 'set multiplot layout 3,1'
write(gnuplotchannel,*) 'set pm3d map'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set logscale x'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,FMT="(A,F7.3,',',F7.3,A)") &
        'set label 1 "Vals = (',vals(1),vals(2),')" at screen 0.1,0.900'
write(gnuplotchannel,FMT="(A,I6,A)") &
        'set label 2 "Ntraj = ', Ntraj, '" at screen 0.1,0.875'
!write(gnuplotchannel,*) 'min_x = ', min_rmsd_x
!write(gnuplotchannel,*) 'max_x = ', max_rmsd_x
!write(gnuplotchannel,*) 'min_y = ', min_rmsd_fx
!write(gnuplotchannel,*) 'max_y = ', max_rmsd_fx
write(gnuplotchannel,*) 'min_x = ', min_rmsd_vals(2)
write(gnuplotchannel,*) 'max_x = ', max_rmsd_vals(2)
write(gnuplotchannel,*) 'min_y = ', min_rmsd_vals(6)
write(gnuplotchannel,*) 'max_y = ', max_rmsd_vals(6)
write(gnuplotchannel,*) 'min_cx = ', min_Ninterpolation
write(gnuplotchannel,*) 'max_cx = ', max_Ninterpolation
write(gnuplotchannel,*) 'set xrange [min_x:max_x]'
write(gnuplotchannel,*) 'set yrange [min_y:max_y]'
write(gnuplotchannel,*) 'set cbrange [min_cx:max_cx+1]'
write(gnuplotchannel,*) 'set palette defined (min_cx "blue", max_cx "red")'
write(gnuplotchannel,*) 'set xlabel "Largest RMSD Between the Input Frame and the Interpolation Points"'
write(gnuplotchannel,*) 'set ylabel "RMSD Between the Output Gradient and the Interpolation"'
write(gnuplotchannel,*) 'set cblabel "Number of Frames Used to Interpolate"'
write(gnuplotchannel,*) 'set xtics ('//&
                                         '"1e-8" .00000001, '//&
                                         '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,*) 'set ytics ('//&
!                                        '"1e-8" .00000001, '//&
!                                        '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,*) 'plot "'//gridpath5//vals_interpolation_text//interpolationfile//'" u '//&
                               '(($5)):7:3 w p lw 6 palette'
write(gnuplotchannel,*) 'set logscale x'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,*) 'set xlabel "SQRT of Weighted Largest RMSD Squared Between the Input Frame and the Interpolation Points"'
write(gnuplotchannel,*) 'set ylabel "RMSD Between the Output Gradient and the Interpolation"'
!write(gnuplotchannel,*) 'min_x = ', min_rmsd_y
!write(gnuplotchannel,*) 'max_x = ', max_rmsd_y
write(gnuplotchannel,*) 'min_x = ', min_rmsd_vals(1)
write(gnuplotchannel,*) 'max_x = ', max_rmsd_vals(1)
write(gnuplotchannel,*) 'set xrange [0.5*sqrt(min_x):2*sqrt(max_x)]'
write(gnuplotchannel,*) 'plot "'//gridpath5//vals_interpolation_text//interpolationfile//'" u '//&
                               '(sqrt($4)):7:3 w p lw 6 palette'
write(gnuplotchannel,*) 'set logscale x'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,*) 'set xlabel "Weighted RMSD Between the Input Frame and the Interpolation Points"'
write(gnuplotchannel,*) 'set ylabel "RMSD Between the Output Gradient and the Interpolation"'
!write(gnuplotchannel,*) 'min_x = ', min_rmsd_z
!write(gnuplotchannel,*) 'max_x = ', max_rmsd_z
write(gnuplotchannel,*) 'min_x = ', min_rmsd_vals(5)
write(gnuplotchannel,*) 'max_x = ', max_rmsd_vals(5)
write(gnuplotchannel,*) 'set xrange [0.5*min_x:2*max_x]'
write(gnuplotchannel,*) 'plot "'//gridpath5//vals_interpolation_text//interpolationfile//'" u '//&
                               '(($6)):7:3 w p lw 6 palette'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 2400,1200'
write(gnuplotchannel,*) 'set output "'//gridpath4//"heatmap"//PNGfilename//'.png"'
write(gnuplotchannel,*) 'set title "RMSD Comparison of a Frame and Gradient with Interpolation"'
write(gnuplotchannel,*) 'set multiplot layout 1,2'
write(gnuplotchannel,*) 'set pm3d map'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,FMT="(A,F7.3,',',F7.3,A)") &
        'set label 1 "Vals = (',vals(1),vals(2),')" at screen 0.1,0.800 front'
write(gnuplotchannel,FMT="(A,I6,A)") &
        'set label 2 "Ntraj = ', Ntraj, '" at screen 0.1,0.775 front'
write(gnuplotchannel,*) 'xmin = ', min_rmsd_vals(1)
write(gnuplotchannel,*) 'xmax = ', max_rmsd_vals(1)
write(gnuplotchannel,*) 'ymin = ', min_rmsd_vals(5)
write(gnuplotchannel,*) 'ymax = ', max_rmsd_vals(5)
!write(gnuplotchannel,*) 'xmin = ', min_rmsd_z
!write(gnuplotchannel,*) 'xmax = ', max_rmsd_z
!write(gnuplotchannel,*) 'ymin = ', min_rmsd_x
!write(gnuplotchannel,*) 'ymax = ', max_rmsd_x
write(gnuplotchannel,*) 'min_cx = ', 1.0e-7
write(gnuplotchannel,*) 'max_cx = ', max(max_rmsd_vals(6),1.0e-7)
!write(gnuplotchannel,*) 'max_cx = ', max(max_rmsd_fx,1.0e-7)
write(gnuplotchannel,*) 'set cbrange [log10(.0000001):log10(.0001)]'
write(gnuplotchannel,*) 'set palette defined ('//&
                        'log10(.0000001) "white", '//&
                        'log10(.0000005) "yellow", '//&
                        'log10(.000001) "green", '//&
                        'log10(.000005) "cyan", '//&
                        'log10(.00001) "blue", '//&
                        'log10(.00005) "magenta", '//&
                        'log10(.0001) "red"'//&
                        ')'
write(gnuplotchannel,*) 'set xlabel "RMSD Variable 1 from Frame"'
write(gnuplotchannel,*) 'set ylabel "RMSD Variable 2 from Frame"'
write(gnuplotchannel,*) 'set cblabel "Maximum RMSD Between the Output Gradient and the Interpolation"'
write(gnuplotchannel,*) 'set xtics ('//&
                                           '"1e-7" log10(0.0000001/(xmin)), '//&
                                           '"1e-6" log10(0.000001/(xmin)), '//&
                                           '"1e-5" log10(0.00001/(xmin)), '//&
                                           '"1e-4" log10(0.0001/(xmin)), '//&
                                           '"1e-3" log10(0.001/(xmin)), '//&
                                           '"1e-2" log10(0.01/(xmin)), '//&
                                           '"1e-1" log10(0.1/(xmin)), '//&
                                               '"1e0" log10(1.0/(xmin)), '//&
                                               '"1e1" log10(10.0/(xmin)))'
write(gnuplotchannel,*) 'set ytics ('//&
                                           '"1e-7" log10(0.0000001/(ymin)), '//&
                                           '"1e-6" log10(0.000001/(ymin)), '//&
                                           '"1e-5" log10(0.00001/(ymin)), '//&
                                           '"1e-4" log10(0.0001/(ymin)), '//&
                                           '"1e-3" log10(0.001/(ymin)), '//&
                                           '"1e-2" log10(0.01/(ymin)), '//&
                                           '"1e-1" log10(0.1/(ymin)), '//&
                                               '"1e0" log10(1.0/(ymin)), '//&
                                               '"1e1" log10(10.0/(ymin)))'
write(gnuplotchannel,*) 'set cbtics ('//&
                                         '"1e-8" log10(.00000001), '//&
                                         '"5e-8" log10(.00000005), '//&
                                          '"1e-7" log10(.0000001), '//&
                                          '"5e-7" log10(.0000005), '//&
                                           '"1e-6" log10(.000001), '//&
                                           '"5e-6" log10(.000005), '//&
                                           '"1e-5"  log10(.00001), '//&
                                           '"5e-5"  log10(.00005), '//&
                                           '"1e-4"   log10(.0001), '//&
                                           '"5e-4"   log10(.0005), '//&
                                           '"1e-3"    log10(.001), '//&
                                           '"5e-3"    log10(.005), '//&
                                           '"1e-2"     log10(.01), '//&
                                           '"5e-2"     log10(.05), '//&
                                           '"1e-1"      log10(.1), '//&
                                           '"5e-1"      log10(.5), '//&
                                           ' "1.0"       log10(1), '//&
                                   ')'
write(gnuplotchannel,*) 'splot "'//gridpath5//'heatmap_rmsd.dat" u '//&
                                '1:2:(log10($3)) w image palette'
write(gnuplotchannel,FMT="(A)") 'unset label 1'
write(gnuplotchannel,FMT="(A,F9.5,A,F9.5,A,F9.5,A)") &
        'set label 1 "f(x,y) = ', RMSDheatmap_coeff(1),'x + ',&
                                  RMSDheatmap_coeff(2),'y + ',&
                                  RMSDheatmap_coeff(3),'" at screen 0.7,0.800 front'
write(gnuplotchannel,FMT="(A)") 'unset label 2'
write(gnuplotchannel,*) 'splot "'//gridpath5//'heatmap_rmsd.dat" u '//&
                                '1:2:(log10($4)) w image palette'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 2400,1200'
write(gnuplotchannel,*) 'set output "'//gridpath4//"heatmap_freq"//PNGfilename//'.png"'
write(gnuplotchannel,*) 'set multiplot layout 1,2'
write(gnuplotchannel,*) 'set title "RMSD Comparison of a Frame and Gradient with Interpolation"'
write(gnuplotchannel,*) 'set pm3d map'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,FMT="(A,F7.3,',',F7.3,A)") &
        'set label 1 "Vals = (',vals(1),vals(2),')" at screen 0.1,0.800 front'
write(gnuplotchannel,FMT="(A,I6,A)") &
        'set label 2 "Ntraj = ', Ntraj, '" at screen 0.1,0.775 front'
write(gnuplotchannel,*) 'xmin = ', min(min_rmsd_vals(4),min_rmsd_vals(6))
write(gnuplotchannel,*) 'xmax = ', max(max_rmsd_vals(4),max_rmsd_vals(6))
write(gnuplotchannel,*) 'ymin = xmin'
write(gnuplotchannel,*) 'ymax = xmax'
write(gnuplotchannel,*) 'set xrange [0:log10(xmax/xmin)]'
write(gnuplotchannel,*) 'set yrange [0:log10(ymax/ymin)]'
write(gnuplotchannel,*) 'set autoscale cb'
!write(gnuplotchannel,*) 'xmin = ', min_rmsd_z
!write(gnuplotchannel,*) 'xmax = ', max_rmsd_z
!write(gnuplotchannel,*) 'ymin = ', min_rmsd_x
!write(gnuplotchannel,*) 'ymax = ', max_rmsd_x
!write(gnuplotchannel,*) 'max_cx = ', max(max_rmsd_fx,1.0e-7)
write(gnuplotchannel,*) 'set xlabel "RMSD Between Approximate and Real Gradient (N=1)"'
write(gnuplotchannel,*) 'set ylabel "RMSD Between Approximate and Real Gradient (N>=1)"'
write(gnuplotchannel,*) 'set cblabel "Occurence"'
write(gnuplotchannel,*) 'set palette defined ('//&
                        '0 "white", '//&
                        '1 "yellow", '//&
                        '2 "green", '//&
                        '3 "cyan", '//&
                        '4 "blue", '//&
                        '5 "magenta", '//&
                        '6 "red"'//&
                        ')'
write(gnuplotchannel,*) 'set xtics ('//&
                                           '"1e-7" log10(0.0000001/(xmin)), '//&
                                           '"1e-6" log10(0.000001/(xmin)), '//&
                                           '"1e-5" log10(0.00001/(xmin)), '//&
                                           '"1e-4" log10(0.0001/(xmin)), '//&
                                           '"1e-3" log10(0.001/(xmin)), '//&
                                           '"1e-2" log10(0.01/(xmin)), '//&
                                           '"1e-1" log10(0.1/(xmin)), '//&
                                               '"1e0" log10(1.0/(xmin)), '//&
                                               '"1e1" log10(10.0/(xmin)))'
write(gnuplotchannel,*) 'set ytics ('//&
                                           '"1e-7" log10(0.0000001/(ymin)), '//&
                                           '"1e-6" log10(0.000001/(ymin)), '//&
                                           '"1e-5" log10(0.00001/(ymin)), '//&
                                           '"1e-4" log10(0.0001/(ymin)), '//&
                                           '"1e-3" log10(0.001/(ymin)), '//&
                                           '"1e-2" log10(0.01/(ymin)), '//&
                                           '"1e-1" log10(0.1/(ymin)), '//&
                                               '"1e0" log10(1.0/(ymin)), '//&
                                               '"1e1" log10(10.0/(ymin)))'
write(gnuplotchannel,*) 'splot "'//gridpath5//'heatmap_freq1_rmsd.dat" u '//&
                                '1:2:3 w image palette'
write(gnuplotchannel,*) 'unset pm3d'
write(gnuplotchannel,*) 'xmin = ', min_rmsd_vals(7)
write(gnuplotchannel,*) 'xmax = ', max_rmsd_vals(7)
write(gnuplotchannel,*) 'deltax = (log10(xmax/xmin)) / ', Nbins/5
write(gnuplotchannel,*) 'set boxwidth deltax'
write(gnuplotchannel,*) 'set style fill solid 1.0'
write(gnuplotchannel,*) 'set xrange [-deltax:log10(xmax/xmin)+deltax]'
write(gnuplotchannel,*) 'unset xtics'
write(gnuplotchannel,*) 'unset ytics'
write(gnuplotchannel,*) 'set ytics autofreq'
write(gnuplotchannel,*) 'set autoscale ymax'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set xtics ('//&
                                           '"1e-3" log10(0.001/(xmin)), '//&
                                           '"5e-3" log10(0.005/(xmin)), '//&
                                           '"1e-2" log10(0.01/(xmin)), '//&
                                           '"5e-2" log10(0.05/(xmin)), '//&
                                           '"1e-1" log10(0.1/(xmin)), '//&
                                           '"5e-1" log10(0.5/(xmin)), '//&
                                               '"1e0" log10(1.0/(xmin)), '//&
                                               '"5e0" log10(5.0/(xmin)), '//&
                                               '"1e1" log10(10.0/(xmin)), '//&
                                               '"5e1" log10(50.0/(xmin)), '//&
                                               '"1e2" log10(100.0/(xmin)), '//&
                                               '"5e2" log10(500.0/(xmin)))'
write(gnuplotchannel,*) 'set xlabel "Ratio of RMSD Between Approximate and Real Gradient (N=1) and (N>=1)"'
write(gnuplotchannel,*) 'set ylabel "Occurence"'
write(gnuplotchannel,*) 'plot "'//gridpath5//'heatmap_freq2_rmsd.dat" u '//&
                               '1:2 w boxes'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)




!open(gnuplotchannel,file=gridpath0//gnuplotfile)
!write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
!write(gnuplotchannel,*) 'set output "'//gridpath0//PNGfilename//'_1.png"'
!write(gnuplotchannel,*) 'set title "Ratio of RMSD Between the Frame and the Gradient"'
!write(gnuplotchannel,*) 'unset key'
!write(gnuplotchannel,*) 'scaling = ', frames
!write(gnuplotchannel,*) 'xmin = ', min_rmsd_fx / max_rmsd_x
!write(gnuplotchannel,*) 'xmax = ', max_rmsd_fx / min_rmsd_x
!write(gnuplotchannel,*) 'set xlabel "Ratio of RMSD"'
!write(gnuplotchannel,*) 'set ylabel "Occurence"'
!
!write(gnuplotchannel,*) 'Nbins = ', Nbins
!write(gnuplotchannel,*) 'bin_width = ', bin_width
!write(gnuplotchannel,*) 'set boxwidth 1'
!write(gnuplotchannel,*) 'set xrange [-0.5:Nbins+.05]'
!write(gnuplotchannel,*) 'set yrange [0:]'
!write(gnuplotchannel,*) 'set style fill transparent solid 0.5'
!
!write(gnuplotchannel,*) 'set xtics ('//&
!                                           '"1e-5" (-5-log10(xmin))/bin_width, '//&
!                                           '"1e-4" (-4-log10(xmin))/bin_width, '//&
!                                           '"1e-3" (-3-log10(xmin))/bin_width, '//&
!                                           '"1e-2" (-2-log10(xmin))/bin_width, '//&
!                                           '"1e-1" (-1-log10(xmin))/bin_width, '//&
!                                               '"1e0" (0-log10(xmin))/bin_width, '//&
!                                               '"1e1" (1-log10(xmin))/bin_width, '//&
!                                               '"1e2" (2-log10(xmin))/bin_width, '//&
!                                               '"1e3" (3-log10(xmin))/bin_width, '//&
!                                               '"1e4" (4-log10(xmin))/bin_width, '//&
!                                               '"1e5" (5-log10(xmin))/bin_width)'
!
!write(gnuplotchannel,*) 'plot "'//gridpath0//"percent_rmsd"//Ngrid_text//'.dat" u '//&
!                               '3:(1.0/scaling) '//&
!                               'smooth frequency w boxes lc rgb "green"'
!close(gnuplotchannel)
!
!call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)


!open(gnuplotchannel,file=gridpath0//gnuplotfile)
!write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
!write(gnuplotchannel,*) 'set output "'//gridpath0//PNGfilename//'_2.png"'
!write(gnuplotchannel,*) 'set title "RMSD Ratio Dependence on Variables"'
!write(gnuplotchannel,*) 'set pm3d map'
!write(gnuplotchannel,*) 'unset key'
!write(gnuplotchannel,*) 'min_x = ', min_rmsd_x
!write(gnuplotchannel,*) 'max_x = ', max_rmsd_x
!write(gnuplotchannel,*) 'min_y = ', min_rmsd_fx
!write(gnuplotchannel,*) 'max_y = ', max_rmsd_fx
!write(gnuplotchannel,*) 'min_cx = ', min_Nbin
!write(gnuplotchannel,*) 'max_cx = ', max_Nbin
!write(gnuplotchannel,*) 'set xrange [0:',var_maxvar(1),']'
!write(gnuplotchannel,*) 'set yrange [0:',var_maxvar(2),']'
!write(gnuplotchannel,*) 'set cbrange [min_cx:max_cx]'
!write(gnuplotchannel,*) 'set palette defined (min_cx "blue", max_cx "red")'
!write(gnuplotchannel,*) 'set xlabel "Var1 (A)"'
!write(gnuplotchannel,*) 'set ylabel "Var2 (A)"'
!write(gnuplotchannel,*) 'set cblabel "Ratio of RMSD"'
!
!write(gnuplotchannel,*) 'plot "'//gridpath0//"percent_rmsd"//Ngrid_text//'.dat" u '//&
!                               '1:2:3 w p lw 4 palette'
!close(gnuplotchannel)
!
!call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)

end subroutine getRMSDinterpolation

subroutine getRMSDErrorPlots(PNGfilename)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!NUMBER OF FRAMES IN THE DATA
real(dp),dimension(7) :: min_rmsd_vals, max_rmsd_vals

!FORMAT OF PNG FILES TO BE MADE
character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: PNGfilename
integer, dimension(comparison_number,3) :: counters

!I/O HANDLING
integer :: iostate

!INTEGER INCREMENTALS
integer :: n

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

call processMultipleInterpolationFiles(counters)

!    read(filechannel2,FMT=*,iostat=iostate) &
!            vals, Ninterpolation, &
!!           rsv1, rsv2, &
!!           rmsd_best, rmsd_interpolated, &
!!           error_best, error_interpolated
!            rsv(3), rsv(4), &
!            rsv(1), rsv(2), &
!            rsv(5), rsv(6)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 2400,1800'
write(gnuplotchannel,*) 'set output "'//gridpath4//PNGfilename//'_1.png"'
write(gnuplotchannel,*) 'set title "Error Convergence as Interpolations Get Closer" font ",32"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set logscale x'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,*) 'set xlabel "RMSD (A) Between Target and Interpolated Frame" font ",24"'
write(gnuplotchannel,*) 'set ylabel "Error (eV/A) Between Target and Interpolated Gradient" font ",24"'
write(gnuplotchannel,*) 'set xtics font ",16"'
write(gnuplotchannel,*) 'set ytics font ",16"'
write(gnuplotchannel,*) 'set xtics ('//&
                                         '"1e-8" .00000001, '//&
                                         '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,*) 'set ytics ('//&
!                                        '"1e-8" .00000001, '//&
!                                        '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,*) 'plot "'//gridpath5//"tmp"//interpolationfile//'" u '//&
                               '3:7 w p'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 2400,1800'
write(gnuplotchannel,*) 'set output "'//gridpath4//PNGfilename//'_2.png"'
write(gnuplotchannel,*) 'set title "Error Convergence as Candidates Get Closer" font ",32"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set logscale x'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,*) 'set xlabel "RMSD (A) Between Target and Candidate Frame" font ",24"'
write(gnuplotchannel,*) 'set ylabel "Error (eV/A) Between Target and Interpolated Gradient" font ",24"'
write(gnuplotchannel,*) 'set xtics font ",16"'
write(gnuplotchannel,*) 'set ytics font ",16"'
write(gnuplotchannel,*) 'set xtics ('//&
                                         '"1e-8" .00000001, '//&
                                         '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,*) 'set ytics ('//&
!                                        '"1e-8" .00000001, '//&
!                                        '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,*) 'plot "'//gridpath5//"tmp"//interpolationfile//'" u '//&
                               '2:6 w p'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 2400,1800'
write(gnuplotchannel,*) 'set output "'//gridpath4//PNGfilename//'_12.png"'
write(gnuplotchannel,FMT="(A)") 'set multiplot layout 1'//&
                        ',2 columnsfirst margins 0.1,0.95,.1,.9 spacing 0,0.1 title '//&
                        '"Error Convergence as Candidates/Interpolations Get Closer" font ",36" offset 0,9'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set logscale x'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,*) 'set xlabel "RMSD (A) Between Target and Candidate Frame"'//&
        ' font ",24" offset "0,3"'
write(gnuplotchannel,*) 'set ylabel "Error (eV/A) Between Target and Candidate/Interpolated Gradient"'//&
        ' font ",24" offset -3,0'
write(gnuplotchannel,*) 'set xtics font ",16"'
write(gnuplotchannel,*) 'set ytics font ",16"'
write(gnuplotchannel,*) 'set yrange [.001:3.0]'
write(gnuplotchannel,*) 'set xtics ('//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                   ')'
write(gnuplotchannel,*) 'set ytics ('//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,*) 'plot "'//gridpath5//"tmp"//interpolationfile//'" u '//&
                               '2:6 w p'

write(gnuplotchannel,*) 'unset ylabel'
write(gnuplotchannel,*) 'unset ytics'
write(gnuplotchannel,*) 'set xtics ('//&
                                         '"1e-8" .00000001, '//&
                                          '"1e-7" .0000001, '//&
                                           '"1e-6" .000001, '//&
                                           '"1e-5"  .00001, '//&
                                           '"1e-4"   .0001, '//&
                                           '"1e-3"    .001, '//&
                                           '"1e-2"     .01, '//&
                                           '"1e-1"      .1, '//&
                                   ')'
write(gnuplotchannel,*) 'set xlabel "RMSD (A) Between Target and Interpolated Frame"'//&
        ' font ",24" offset "0,3"'
write(gnuplotchannel,*) 'plot "'//gridpath5//"tmp"//interpolationfile//'" u '//&
                               '3:7 w p'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 2400,1800'
write(gnuplotchannel,*) 'set output "'//gridpath4//PNGfilename//'_CM2.png"'
write(gnuplotchannel,*) 'set title "Error Convergence as Candidates Get Closer" font ",32"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set logscale x'
write(gnuplotchannel,*) 'set logscale y'
write(gnuplotchannel,*) 'set xlabel "Coulomb Matrix Difference (1/A) Between Target and Candidate Frame" font ",24"'
write(gnuplotchannel,*) 'set ylabel "Error (eV/A) Between Target and Interpolated Gradient" font ",24"'
write(gnuplotchannel,*) 'set xtics font ",16"'
write(gnuplotchannel,*) 'set ytics font ",16"'
write(gnuplotchannel,*) 'set xtics ('//&
                                         '"1e-8" .00000001, '//&
                                         '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,*) 'set ytics ('//&
!                                        '"1e-8" .00000001, '//&
!                                        '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1e0"       1, '//&
                                   ')'
write(gnuplotchannel,*) 'plot "'//gridpath5//"tmp"//interpolationfile//'" u '//&
                               '9:6 w p'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)




end subroutine getRMSDErrorPlots



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               processInterpolationFile
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine looks at the processed interpolation file (truncated file)
!               and creates data files for plotting, excising anything outside of the
!               region of the phase space of interest
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               vals                            REAL,DIM(Nvar)                  The region of the phase space we are
!                                                                               interested in plotting
!               delta_vals                      REAL,DIM(Nvar)                  The size of the region of the phase
!                                                                               space of interest
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!               RMSDheatmap_coeff               REAL,DIM(3)                     The coefficients of the linear fit
!                                                                               produced by this data
!               min_rsv                         REAL,DIM(7)                     Minimums of the RMSD special variables
!               max_rsv                         REAL,DIM(7)                     Maximums of the RMSD special variables
!               min_Ninterpolation              INTEGER                         Minimum of the Ninterpolation
!               max_Ninterpolation              INTEGER                         Maximum of the Ninterpolation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               Nbins                           INTEGER                         The number of bins for the binned data
!               rsv_binwidth                    REAL,DIM(7)                     Width of the binning for each RSV
!               rsv_bin                         INTEGER,DIM(7)                  Temporary buffer for each RSV bin data
!               ERheatmap2D                     REAL(Nbins,Nins)                Heatmap that describes the occurence
!                                                                               of two RSV
!               ERheatmap1D                     REAL(Nbins/5)                   Heatmap that describes the occurence
!                                                                               of one RSV
!               IEheatmap2D                     REAL(Nbins,Nins)                Heatmap that describes the maximum RSV
!                                                                               (usually IE) observed given two RSV
!               Nheatmap                        INTEGER                         The number of unique nonzero entries
!                                                                               in IEheatmap2D
!               A                               REAL(Nheatmap,3)                The "independent variable matrix" for
!                                                                               IEheatmap2D used for linear fitting
!               b                               REAL(Nheatmap)                  The "dependent variable matrix" for
!                                                                               IEheatmap2D used for linear fitting
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!		gridpath5//truncated//          DAT			Holds processed interpolation data
!		       interpolationfile
!               gridpath4//JPGfilename          PNG                     Generic format for image names
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine processInterpolationFile(&
                vals,delta_vals,&
                RMSDheatmap_coeff,&
                min_rsv,max_rsv,&
                min_Ninterpolation,max_Ninterpolation)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!COLLECTIVE VARIABLES TO FILTER DATA
real(dp),dimension(Nvar),intent(in) :: vals, delta_vals

!NUMBER OF FRAMES IN THE DATA
real(dp),dimension(7),intent(out) :: min_rsv,max_rsv
real(dp),dimension(3),intent(out) :: RMSDheatmap_coeff
real(dp),dimension(7) :: rsv, rsv_binwidth
integer,dimension(7) :: rsv_bin

!FORMATTING OF PNG FILES
character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5

!VARIABLES IN THE INTERPOLATION FILE
integer :: Ninterpolation
integer,intent(out) :: min_Ninterpolation
integer,intent(out) :: max_Ninterpolation
real(dp),dimension(Nvar) :: current_vals

!I/O HANDLING
integer :: iostate

!HISTOGRAM VARIABLES
integer :: Nbins
integer :: i,j,k,l,Nheatmap
real(dp),allocatable :: IEheatmap2D(:,:)
integer,allocatable :: ERheatmap2D(:,:),ERheatmap1D(:)
real(dp),allocatable :: A(:,:), b(:)

!INTEGER INCREMENTALS
integer :: n

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

!The truncated interpolation file is more
!manageable
open(filechannel2,file=gridpath5//&
        "truncated"//interpolationfile)

!In particular, the first line describes the
!bounds of the data
read(filechannel2,FMT="(1x,3(I9,1x),2(I6,1x),"//&
        "2(7(G16.6,1x)))",iostat=iostate) &
        i,j,k,&
        min_Ninterpolation, max_Ninterpolation,&
        min_rsv, max_rsv

!The binning is done logarithmically and the
!number of bins is temporarily hardcoded
!The 2D heatmap looks better with higher
!resolution (compared to the 1D heatmap)
Nbins = 100
rsv_binwidth = log10(max_rsv/min_rsv)*1.0d0/Nbins
allocate(ERheatmap2D(Nbins,Nbins),&
         ERheatmap1D(Nbins/5),&
         IEheatmap2D(Nbins,Nbins))
ERheatmap2D = 0
ERheatmap1D = 0
IEheatmap2D = 1.0d-7

do 
    read(filechannel2,FMT=*,iostat=iostate) &
            current_vals, Ninterpolation, rsv
    if (iostate /= 0) exit

    !If the frame is not in the region of the phase
    !space of interest, then skip it
    if (any(abs(current_vals-vals) > &
            delta_vals)) cycle

    !Logarithms involve division by the minimum to
    !get a correct binning
    !The if statement is to make sure the edge cases
    !are counted in the distribution correctly
    rsv_bin = floor(log10(rsv/min_rsv)/rsv_binwidth) + 1
    do l = 1, 7
        if (rsv_bin(l) < 1) rsv_bin(l) = 1
        if (rsv_bin(l) > Nbins) rsv_bin(l) = Nbins
    end do

    !The variables chosen for the distributions are
    !not arbitrarily picked (but may be if wanted)
    ERheatmap2D(rsv_bin(5),rsv_bin(6)) =&
            ERheatmap2D(rsv_bin(5),rsv_bin(6)) + 1
    ERheatmap1D(rsv_bin(7)/5) =&
            ERheatmap1D(rsv_bin(7)/5) + 1

    !This last heatmap records the maximum of some
    !third value, not the total occurence (ideally
    !for maxmimum error checking)
    IEheatmap2D(rsv_bin(3),rsv_bin(4)) = max(&
            IEheatmap2D(rsv_bin(3),rsv_bin(4)),rsv(6))
end do
close(filechannel2)


!Then just bin the data normally

open(filechannel2,file=gridpath5//"ERheatmap2D.dat")
do i = 1, Nbins
    do j = 1, Nbins
        write(filechannel2,FMT=*)&
                (i-0.5)*rsv_binwidth(5),&
                (j-0.5)*rsv_binwidth(6),&
                ERheatmap2D(i,j)
    end do
    write(filechannel2,*) ""
end do
close(filechannel2)

open(filechannel2,file=gridpath5//"ERheatmap1D.dat")
do i = 1, Nbins/5
    write(filechannel2,FMT=*)&
            (i-0.5)*rsv_binwidth(7),&
            ERheatmap1D(i)
end do
close(filechannel2)

!For the third heatmap we are also interested in
!finding a fit for it, so we need to figure out
!how many nonzero entries are in it that we
!need to fit to

Nheatmap = 0
do i = 1, Nbins
    do j = 1, Nbins
        !We can see whether an entry is
        !nonzero by checking whether it is
        !greater than its initial value
        if (IEheatmap2D(i,j) > 1.0d-7) then
             Nheatmap = Nheatmap + 1
        end if
    end do
end do

!We store this data is matrices A and b
!A is the "independent data"
!b is the "dependent data"
!and we are trying to find a set of
!coefficients in matrix RMSDheatmap_coeff
!that, multiplied with A, produces a
!matrix close to b
allocate(A(Nheatmap,3),b(Nheatmap))
Nheatmap = 0
do i = 1, Nbins
    do j = 1, Nbins
        if (IEheatmap2D(i,j) > 1.0d-7) then
             Nheatmap = Nheatmap + 1

             !In this case, our fit is
             ! c2 * log10(rsv3) +
             !   c1 * log10(rsv4) +
             !                   c0 = log10(b)
             A(Nheatmap,:) = (/&
                     (i-0.5)*rsv_binwidth(3),&
                     (j-0.5)*rsv_binwidth(4),&
                     1.0d0/)
             b(Nheatmap) = log10(IEheatmap2D(i,j))
        end if
    end do
end do

!This fit is linear so can be obtained by doing
!a least squares regression on the above matrices
call LS(A,Nheatmap,3,b,RMSDheatmap_coeff)


!And then just plot the two sets of data (the
!interpolation and the fit) side by side

open(filechannel2,file=gridpath5//"IEheatmap2D.dat")
do i = 1, Nbins
    do j = 1, Nbins
        write(filechannel2,FMT=*)&
                (i-0.5)*rsv_binwidth(3),&
                (j-0.5)*rsv_binwidth(4),&
                IEheatmap2D(i,j),&
                10.0d0**(max((i-0.5)*rsv_binwidth(3)*RMSDheatmap_coeff(1)+&
                             (j-0.5)*rsv_binwidth(4)*RMSDheatmap_coeff(2)+&
                                            RMSDheatmap_coeff(3),-7.0d0))
    end do
    write(filechannel2,*) ""
end do
close(filechannel2)

end subroutine processInterpolationFile



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               processInterpolationFile2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine looks at the interpolation file, gets minimums and
!               maximums, generates secondary variables, and writes the data down in
!               a separate truncated file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               vals                            REAL,DIM(Nvar)                  The region of the phase space
!                                                                               (described by collective variables)
!               rsv                             REAL,DIM(7)                     The variables in order are:
!                                                                               1.RMSD of accept best frame
!                                                                               2.RMSD of interpolated frame
!                                                                               3.Special Variable 1
!                                                                               4.Special Variable 2
!                                                                               5.Error of accept best gradient
!                                                                               6.Error of interpolated gradient
!                                                                               7.Error of accept best / interpolated
!               min_rsv                         REAL,DIM(7)                     Minimums of the RMSD special variables
!               max_rsv                         REAL,DIM(7)                     Maximums of the RMSD special variables
!               Ninterpolation                  INTEGER                         Number of frames used in the
!                                                                               interpolation
!               min_Ninterpolation              INTEGER                         Minimum of the Ninterpolation
!               max_Ninterpolation              INTEGER                         Maximum of the Ninterpolation
!               frames                          INTEGER                         Number of frames in the interpolation
!                                                                               file
!               neginfinity_counter             INTEGER                         Number of frames where RSV5 = 0
!               posinfinity_counter             INTEGER                         Number of frames where RSV6 = 0
!               firstliner                      CHARACTER()                     A string with minimums and maximums
!                                                                               of the RSV and other data appended as
!                                                                               the first line of the truncated file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!		gridpath5//interpolationfile    DAT			Holds the raw interpolation data
!		gridpath5//truncated//          DAT			Holds processed interpolation data
!		       interpolationfile
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine processInterpolationFile2()
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!COLLECTIVE VARIABLES TO FILTER DATA
character(15) :: vals_interpolation_text
character(12) :: short_prefix_text

!VARIABLES IN THE INTERPOLATION FILE
real(dp),dimension(Nvar) :: vals
real(dp),dimension(9) :: rsv, min_rsv, max_rsv
integer :: min_Ninterpolation,max_Ninterpolation,Ninterpolation
integer :: Nswitch

integer :: frames,neginfinity_counter,posinfinity_counter
character(1+3*10+2*7+2*9*17) :: firstliner

character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5

!I/O HANDLING
integer :: iostate

!INTEGER INCREMENTALS
integer :: n

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

min_Ninterpolation = 1000
max_Ninterpolation = 0

min_rsv = 1.0d9
max_rsv = 0.0d0

frames = 0
neginfinity_counter = 0
posinfinity_counter = 0

!The data we read from the interpolation file, we do some
!minor processing and then record to the truncated interpolation
!file
open(filechannel1,file=gridpath5//"truncated"//interpolationfile)
open(filechannel2,file=gridpath5//interpolationfile)
do 
    !For each frame in the file, the data read is:
    !1.Its collective variables
    !2.Number of frames used for its interpolation
    !3.The RMSD special variables that results from it
    read(filechannel2,FMT=*,iostat=iostate) &
            vals, Ninterpolation, &
!           rsv1, rsv2, &
!           rmsd_best, rmsd_interpolated, &
!           min_CMdiff, interpolated_CMdiff, &
!           error_best, error_interpolated
            rsv(3), rsv(4), &
            rsv(1), rsv(2), &
            rsv(8), rsv(9), &
            rsv(5), rsv(6)
            
    if (iostate /= 0) exit

    rsv(5) = rsv(5) * RU_energy / eV
    rsv(6) = rsv(6) * RU_energy / eV

    !One of the processed variables that this
    !subroutine makes is the ratio between
    !RSV(5) and RSV(6) which means a few
    !conditionals need to be made to treat them

    !This may bias the data so these occurences
    !are recorded as negative and positive
    !infinity (since we will later take a
    !logarithm

    if ((rsv(5) == 0.0d0).and.(rsv(6)==0.0d0)) then
        rsv(7) = 1.0d0
        write(filechannel1,FMT=*) &
                vals, Ninterpolation, rsv

        min_Ninterpolation = min(min_Ninterpolation,&
                Ninterpolation)
        max_Ninterpolation = max(max_Ninterpolation,&
                Ninterpolation)

        do n = 1, 9
            max_rsv(n) = max(max_rsv(n),rsv(n))

            if (rsv(n) == 0.0d0) then
            else
                min_rsv(n) = min(min_rsv(n),rsv(n))
            end if
        end do
        
    else if (rsv(5) == 0.0d0) then
        neginfinity_counter = neginfinity_counter + 1

    else if (rsv(6) == 0.0d0) then
        posinfinity_counter = posinfinity_counter + 1

    else
        rsv(7) = rsv(5)/rsv(6)
        write(filechannel1,FMT=*) &
                vals, Ninterpolation, rsv

        min_Ninterpolation = min(min_Ninterpolation,&
                Ninterpolation)
        max_Ninterpolation = max(max_Ninterpolation,&
                Ninterpolation)

        do n = 1, 9
            max_rsv(n) = max(max_rsv(n),rsv(n))

            if (rsv(n) == 0.0d0) then
            else
                min_rsv(n) = min(min_rsv(n),rsv(n))
            end if
        end do
    
    end if

    !This subroutine also counts the total number of
    !frames in the file
    frames = frames + 1
end do
close(filechannel1)
close(filechannel2)

!All of this metadata on the bounds of the data,
!the positive/negative infinity bias, and the
!total number of frames is recorded in the first
!line of the file
write(firstliner,FMT="(A1,3(I9,1x),2(I6,1x),"//&
        "2(9(E16.6,1x)))") &
        "#",posinfinity_counter,&
        neginfinity_counter, frames, &
        min_Ninterpolation, max_Ninterpolation,&
        min_rsv, max_rsv

!Sed is used to deposit the first line onto the
!file
call system("sed -i '1i\"//firstliner//"' '"//&
        gridpath5//"truncated"//&
        interpolationfile//"'")

end subroutine processInterpolationFile2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               processMultipleInterpolationFiles
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine looks at multiple processed interpolation files (truncated
!               files), calculates lower and upper bounds, bins the data, and writes down the
!               data in separate files for later plotting
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!               counters                        INTEGER                         Number of frames and other counters
!                  (comparison_number,3)                                        in each interpolation file for
!                                                                               comparison
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               comparison_number               INTEGER                         The number of experiments to be
!                                                                               compared (not necessarily unique)
!               vals                            REAL,DIM(Nvar)                  The region of the phase space
!                                                                               (described by collective variables)
!               delta_vals                      REAL,DIM(Nvar)                  The size of the region of the phase
!                                                                               space of interest
!               rsv                             REAL,DIM(7)                     The variables in order are:
!                                                                               1.RMSD of accept best frame
!                                                                               2.RMSD of interpolated frame
!                                                                               3.Special Variable 1
!                                                                               4.Special Variable 2
!                                                                               5.Error of accept best gradient
!                                                                               6.Error of interpolated gradient
!                                                                               7.Error of accept best / interpolated
!               min_rsv                         REAL,DIM(7)                     Minimums of the RMSD special variables
!               max_rsv                         REAL,DIM(7)                     Maximums of the RMSD special variables
!               rsv_binwidth                    REAL,DIM(7)                     Width of the binning for each RSV
!               rsv_bin                         INTEGER,DIM(7)                  Temporary buffer for each RSV bin data
!               Ninterpolation                  INTEGER                         Number of frames used in the
!                                                                               interpolation
!               min_Ninterpolation              INTEGER                         Minimum of the Ninterpolation
!               max_Ninterpolation              INTEGER                         Maximum of the Ninterpolation
!               Ninterpolation_binwidth         REAL                            Width of the binning for
!                                                                               Ninterpolation
!               Ninterpolation_bin              INTEGER                         Temporary buffer for Ninterpolation
!                                                                               bin data
!               frames                          INTEGER                         Number of frames in the interpolation
!                                                                               file
!               frames_trials                   INTEGER                         Number of frames in each of the
!                  (comparison_number)                                          interpolation file for comparison
!               rsv_binning                     REAL,DIM(7,                     The RSV data binned (scaled so
!                                                 3*comparison_number,          to integrate to 1)
!                                                 Nbins)
!               Ninterpolation_binning          REAL,DIM(                       The Ninterpolation data binned
!                                                 3*comparison_number,          (scaled so to integrate to 1)
!                                                 Nbins)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!		gridpath0//comparison_file      DAT			Holds the parameters that describe
!                                                                       what regions of the phase space to
!                                                                       excise data from
!		gridpath0//allprefixes()//      DAT			Holds processed interpolation data
!		       intermediatefolder//
!		       interpolationfile
!		gridpath5//tmp//                DAT			Holds slightly more processed
!		       interpolationfile                                interpolation data
!		gridpath5//SATRVname//          DAT			Holds binned interpolation data
!		       _binning                                         for later binning
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine processMultipleInterpolationFiles(counters)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!COLLECTIVE VARIABLES TO FILTER DATA
real(dp),dimension(Nvar) :: vals, delta_vals, tmpvals
character(expfolder_length-1) :: otherexpfolder
character(70) :: longtext

!FORMATTING OF PNG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text

integer,dimension(comparison_number,3),intent(out) :: counters
integer,dimension(comparison_number) :: frames_trials

real(dp),dimension(9) :: rsv, min_rsv, max_rsv
real(dp),dimension(comparison_number,9) :: lower_rsv, upper_rsv
real(dp),dimension(comparison_number,9) :: arithmetic_mean_rsv,geometric_mean_rsv
integer :: min_Ninterpolation,max_Ninterpolation,Ninterpolation
integer,dimension(comparison_number) :: lower_Ninterpolation,upper_Ninterpolation
real(dp),dimension(comparison_number) :: arithmetic_mean_Ninterpolation
real(dp),dimension(comparison_number) :: geometric_mean_Ninterpolation

character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5

!I/O HANDLING
integer :: iostate

!HISTOGRAM VARIABLES
real(dp) :: Ninterpolation_binwidth
integer :: Ninterpolation_bin
real(dp), dimension(9) :: rsv_binwidth
integer, dimension(9) :: rsv_bin
integer :: frames, Nbins
real(dp),allocatable :: Ninterpolation_binning(:,:)
real(dp),allocatable :: rsv_binning(:,:,:)
integer :: starting_index
real(dp),allocatable :: zero_array(:)

real(dp),dimension(Nvar) :: vals_prev
integer,dimension(Nvar) :: deltavals
integer :: Nswitch
integer :: beyond_counter,ten_counter
integer,dimension(3) :: pos_counter,neg_counter

logical :: debug_flag = .false.

!INTEGER INCREMENTALS
integer :: n,m,l

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

counters = 0
frames_trials = 0
frames = 0

min_Ninterpolation = 1000000
max_Ninterpolation = 0
min_rsv = 1.0d19
max_rsv = 0.0d0

open(filechannel1,file=gridpath0//comparison_file)
open(filechannel3,file=gridpath5//"tmp"//interpolationfile)

!We open a file once per comparison number (even if
!it's a duplicate)
do n = 1, comparison_number

    !The comparison file has information regarding
    !what region of the phase space we want to plot
    read(filechannel1,FMT=*,iostat=iostate) vals, delta_vals
    if (iostate /= 0) then
        !If the amount of data on the comparison file does]
        !not match the comparison number than something
        !went wrong
        print *, ""
        print *, "Bad formatting in interpolation comparison arguments"
        print *, "   Must be 2*Nvar floats"
        print *, ""

        close(filechannel1)
        close(filechannel3)

        return
    end if

    !The string allprefixes has information regarding
    !what experiment we want to plot and the array
    !alllengths has information regarding how long
    !each string in it is
    if (n == 1) then
        starting_index = 0
    else
        starting_index = sum(alllengths(1:n-1))
    end if
    open(filechannel2,file=gridpath0//&
            allprefixes(starting_index+1:sum(alllengths(1:n)))//&
            intermediatefolder//"truncated"//interpolationfile)

    !The first line in the truncated interpolation file
    !has the bounds to the data
    read(filechannel2,FMT="(1x,3(I9,1x),2(I6,1x),"//&
            "2(9(G16.6,1x)))",iostat=iostate) &
            counters(n,1),counters(n,2),counters(n,3),&
            lower_Ninterpolation(n), upper_Ninterpolation(n),&
            lower_rsv(n,:),&
            upper_rsv(n,:)
!           (lower_rsv(n,i),i=1,7),&
!           (upper_rsv(n,i),i=1,7),&
!           (arithmetic_mean_rsv(n,i),i=1,7),&
!           (geometric_mean_rsv(n,i),i=1,7)

    pos_counter = 0
    neg_counter = 0
    beyond_counter = 0
    ten_counter = 0

    min_Ninterpolation = min(min_Ninterpolation,&
            lower_Ninterpolation(n))
    max_Ninterpolation = max(max_Ninterpolation,&
            upper_Ninterpolation(n))
    
    do m = 1, 9
        min_rsv(m) = min(min_rsv(m),lower_rsv(n,m))
        max_rsv(m) = max(max_rsv(m),upper_rsv(n,m))
    end do

    upper_Ninterpolation(n) = min_Ninterpolation
    upper_rsv(n,:) = min_rsv

    !We also excise any data that is not in the region
    !of the phase space of interest
    !All data is stored in the temporary interpolation
    !file; yes, there will be duplicates if the same
    !data is encountered across two trials
    do
        read(filechannel2,FMT=*,iostat=iostate)&
                tmpvals, Ninterpolation, rsv
        if (iostate /= 0) exit

        if (frames_trials(n) > 0) then
    
        deltavals = floor(abs(tmpvals-vals_prev)*divisor(:,1))
        deltavals = abs(floor(tmpvals*divisor(:,2)) - &
                        floor(vals_prev*divisor(:,2)))
        Nswitch = sum(deltavals)+1
    
        if ((Ninterpolation <= 9)) then
            ten_counter = ten_counter + 1
        else
!       if ((Nswitch < 4).and.(Ninterpolation > 9)) then
        if ((Nswitch < 4)) then
        if (rsv(7) > 1.0d0) then
            pos_counter(Nswitch) = pos_counter(Nswitch) + 1
        else
            neg_counter(Nswitch) = neg_counter(Nswitch) + 1
        end if
        else
            beyond_counter = beyond_counter + 1
        end if
        end if
    
        end if
    
        vals_prev = tmpvals

        tmpvals = abs(tmpvals-vals) - abs(delta_vals)

        if (all(delta_vals>0)) then
            if (any(tmpvals > 0)) cycle
        else
            if (all(tmpvals <= 0)) cycle
        end if

        upper_Ninterpolation(n) = max(upper_Ninterpolation(n),&
                Ninterpolation)
        do m = 1, 9
            upper_rsv(n,m) = max(upper_rsv(n,m),rsv(m))
        end do

        !Recording the number of frames in each
        !trial lets us know where each trial begins
        !and ends on this one big file
        frames_trials(n) = frames_trials(n) + 1
        write(filechannel3,FMT=*) Ninterpolation,rsv
    end do
    close(filechannel2)

    if (debug_flag) then
        print *, ""
        print *, allprefixes(starting_index+1:sum(alllengths(1:n)))
        do Nswitch = 1, 3
        print *, "n = ", Nswitch - 1, " neg / pos:", &
                neg_counter(Nswitch), "/", pos_counter(Nswitch), "%neg/total:", &
                neg_counter(Nswitch) *100.0 / (neg_counter(Nswitch) + pos_counter(Nswitch)), "%"
        end do
        print *, ""
        print *, "n > 2:", beyond_counter
        print *, ""
        print *, "Ninterpolation <= 9:", ten_counter
        print *, ""
    end if
end do

close(filechannel1)
close(filechannel3)

!We have to get the bounds to ALL the data
!across all experiments

!min_Ninterpolation = minval(lower_Ninterpolation)
!max_Ninterpolation = maxval(upper_Ninterpolation)
!
!do m = 1, 9
!    min_rsv(m) = minval(lower_rsv(:,m))
!    max_rsv(m) = maxval(upper_rsv(:,m))
!end do

!If any of the trials has no frames, there will
!probably be trouble plotting it
if (any(frames_trials == 0)) then
    print *, "MAJOR ERROR IN TRIALS SELECTED FOR "//&
             "INTERPOLATION ANALYSIS"
    frames_trials = 0
    return
end if

!The number of bins is dynamically chosen so that
!every bin (even at the beginning and end) has an
!integer number of Ninterpolation
!This means there is some fiddling with how the
!range of the data is factored
Nbins = max_Ninterpolation - min_Ninterpolation
do
    if (Nbins < 61) then
        exit
    else if (modulo(Nbins,2) == 0) then
        Nbins = Nbins / 2
    else if (modulo(Nbins,3) == 0) then
        Nbins = Nbins / 3
    else if (modulo(Nbins,5) == 0) then
        Nbins = Nbins / 5
    else
        exit
    end if
end do

!These fallbacks exist so that this dynamic binning
!doesn't get out of hand (but then each bin may not
!have an integer number of Ninterpolation)
if (Nbins == 0) Nbins = 50
if (Nbins > 100) Nbins = 100

!The binwidths are determined from the bounds of
!the data and the number of bins designated to it
Ninterpolation_binwidth = ceiling((max_Ninterpolation -&
        min_Ninterpolation) *1.0d0/ Nbins)
rsv_binwidth = log10(max_rsv/min_rsv)*1.0d0/Nbins

!There is a multiplication by three because three
!subsets of each distribution are differentiated
allocate(Ninterpolation_binning(3*comparison_number,Nbins),&
         rsv_binning(9,3*comparison_number,Nbins),&
         zero_array(3*comparison_number))
zero_array = 0.0d0
Ninterpolation_binning = 0.0d0
rsv_binning = 0.0d0

!do n = 1, comparison_number
!    upper_Ninterpolation(n) = min_Ninterpolation
!    upper_rsv(n,:) = min_rsv
!end do

arithmetic_mean_Ninterpolation = 0.0d0
geometric_mean_Ninterpolation = 0.0d0

arithmetic_mean_rsv = 0.0d0
geometric_mean_rsv = 0.0d0

!The data for ALL trials is stored in the temporary
!interpolation file
open(filechannel1,file=gridpath5//"tmp"//interpolationfile)
do n = 1, comparison_number

    !Each trial's data set is not separated so the
    !number of lines of data read must be counted
    do m = 1, frames_trials(n)
        read(filechannel1,FMT=*) Ninterpolation,rsv

!        upper_Ninterpolation(n) = max(upper_Ninterpolation(n),&
!                Ninterpolation)
        arithmetic_mean_Ninterpolation(n) = &
            arithmetic_mean_Ninterpolation(n) + Ninterpolation
        geometric_mean_Ninterpolation(n) = &
            geometric_mean_Ninterpolation(n) + log10(1.0d0*Ninterpolation)

        do l = 1, 9
!            upper_rsv(n,l) = max(upper_rsv(n,l),rsv(l))
            arithmetic_mean_rsv(n,l) = &
                arithmetic_mean_rsv(n,l) + rsv(l)

            if (rsv(l) == 0.0d0) then
                geometric_mean_rsv(n,l) = &
                    geometric_mean_rsv(n,l) + log10(min_rsv(n))
            else
                geometric_mean_rsv(n,l) = &
                    geometric_mean_rsv(n,l) + log10(rsv(l))
            end if
        end do
        

        !Read each line and bin them accordingly
    
        Ninterpolation_bin = floor((Ninterpolation&
                - min_Ninterpolation)&
                / Ninterpolation_binwidth) + 1
        if (Ninterpolation_bin < 1) Ninterpolation_bin = 1
        if (Ninterpolation_bin > Nbins) Ninterpolation_bin = Nbins

        rsv_bin = floor(log10(rsv/min_rsv)/rsv_binwidth) + 1
        do l = 1, 9
            if (rsv_bin(l) < 1) rsv_bin(l) = 1
            if (rsv_bin(l) > Nbins) rsv_bin(l) = Nbins
        end do

        !The first portion of the data is ALL of
        !the data; this will be in the back of
        !all other data
        if (.true.) then
            Ninterpolation_binning(3*(n-1)+1,Ninterpolation_bin) = &
                    Ninterpolation_binning(3*(n-1)+1,Ninterpolation_bin) + &
                    1.0d0/frames_trials(n)
            do l = 1, 9
                rsv_binning(l,3*(n-1)+1,rsv_bin(l)) =&
                        rsv_binning(l,3*(n-1)+1,rsv_bin(l)) +&
                        1.0d0/frames_trials(n)
            end do
        end if

        !The second portion of the data is that
        !which has this quality (can be
        !arbitrarily picked)
!       if (rsv(7) < 1.0d1) then
        if (rsv(6) <= upper_rsv(n,5)) then
            Ninterpolation_binning(3*(n-1)+2,Ninterpolation_bin) = &
                    Ninterpolation_binning(3*(n-1)+2,Ninterpolation_bin) + &
                    1.0d0/frames_trials(n)
            do l = 1, 9
                rsv_binning(l,3*(n-1)+2,rsv_bin(l)) =&
                        rsv_binning(l,3*(n-1)+2,rsv_bin(l)) +&
                        1.0d0/frames_trials(n)
            end do
        end if

        !The third portion of the data is that
        !which has this quality (can be
        !arbitrarily picked); this will be in
        !the front of all other data
!       if (rsv(7) <= 1.0d0) then
        if (rsv(6) <= 0.5d0 * upper_rsv(n,5)) then
            Ninterpolation_binning(3*n,Ninterpolation_bin) = &
                    Ninterpolation_binning(3*n,Ninterpolation_bin) + &
                    1.0d0/frames_trials(n)
            do l = 1, 9
                rsv_binning(l,3*n,rsv_bin(l)) =&
                        rsv_binning(l,3*n,rsv_bin(l)) +&
                        1.0d0/frames_trials(n)
            end do
        end if
    end do
end do
close(filechannel1)

do n = 1, comparison_number
    arithmetic_mean_rsv(n,:) = &
        arithmetic_mean_rsv(n,:) / frames_trials(n)
    geometric_mean_rsv(n,:) = &
        10.0d0 ** (geometric_mean_rsv(n,:)  / frames_trials(n))
    
    arithmetic_mean_Ninterpolation(n) = &
        arithmetic_mean_Ninterpolation(n) * 1.0d0 / &
        frames_trials(n)
    geometric_mean_Ninterpolation = &
        10.0d0 ** (geometric_mean_Ninterpolation / &
        frames_trials(n))
end do

!Finally, record the data for
!later plotting

open(filechannel1,file=gridpath5//"InterpolationTDD_binning.dat")

if (comparison_upperlimit > comparison_lowerlimit) then
    write(filechannel1,FMT="(A1,2(E16.6,1x))") "#",&
            comparison_lowerlimit, comparison_upperlimit
else
    write(filechannel1,FMT="(A1,2(E16.6,1x))") "#",&
            min_Ninterpolation - 0.5*Ninterpolation_binwidth, &
            max_Ninterpolation + 0.5*Ninterpolation_binwidth
end if

do n = 1, comparison_number
write(filechannel1,FMT="(A1,3(E16.6,1x))") &
        "#",1.0d0*upper_Ninterpolation(n), &
        arithmetic_mean_Ninterpolation(n),&
        geometric_mean_Ninterpolation(n)
end do

write(filechannel1,FMT=*) min_Ninterpolation + &
        Ninterpolation_binwidth*(-0.5), &
        zero_array
do n = 1, Nbins
    write(filechannel1,FMT=*) min_Ninterpolation + &
            Ninterpolation_binwidth*(n-0.5),&
            Ninterpolation_binning(:,n)
end do
write(filechannel1,FMT=*) min_Ninterpolation + &
        Ninterpolation_binwidth*(Nbins+0.5), &
        zero_array
close(filechannel1)

do l = 1, 9

!Accept Best RMSD Distribution
if (l == 1) open(filechannel1,file=gridpath5//&
        "InterpolationARD_binning.dat")
!Interpolation RMSD Distribution
if (l == 2) open(filechannel1,file=gridpath5//&
        "InterpolationIRD_binning.dat")
!RMSD Special Variable 1 Distribution
if (l == 3) open(filechannel1,file=gridpath5//&
        "InterpolationRSV1D_binning.dat")
!RMSD Special Variable 2 Distribution
if (l == 4) open(filechannel1,file=gridpath5//&
        "InterpolationRSV2D_binning.dat")
!Accept Best Error Distribution
if (l == 5) open(filechannel1,file=gridpath5//&
        "InterpolationAED_binning.dat")
!Interpolation Error Distribution
if (l == 6) open(filechannel1,file=gridpath5//&
        "InterpolationIED_binning.dat")
!Relative Error Distribution
if (l == 7) open(filechannel1,file=gridpath5//&
        "InterpolationRED_binning.dat")
!Accept Best Coulomb Matrix Difference Distribution
if (l == 8) open(filechannel1,file=gridpath5//&
        "InterpolationACMDD_binning.dat")
!Interpolation Coulomb Matrix Difference Distribution
if (l == 9) open(filechannel1,file=gridpath5//&
        "InterpolationICMDD_binning.dat")

if (comparison_upperlimit > comparison_lowerlimit) then
    write(filechannel1,FMT="(A1,2(E16.6,1x))") "#",&
            comparison_lowerlimit, comparison_upperlimit
else
    write(filechannel1,FMT="(A1,2(E16.6,1x))") "#",&
            log10(min_rsv(l)) - 0.5*rsv_binwidth(l), &
            log10(max_rsv(l)) + 0.5*rsv_binwidth(l)
end if

do n = 1, comparison_number
write(filechannel1,FMT="(A1,3(E16.6,1x))") &
        "#",log10(upper_rsv(n,l)), &
        log10(arithmetic_mean_rsv(n,l)),&
        log10(geometric_mean_rsv(n,l))
end do

write(filechannel1,FMT=*)&
        log10(min_rsv(l)) +&
        rsv_binwidth(l)*(-0.5d0),&
        zero_array
do n = 1, Nbins
    write(filechannel1,FMT=*)&
            log10(min_rsv(l)) +&
            rsv_binwidth(l)*(n-0.5d0),&
            rsv_binning(l,:,n)
end do
write(filechannel1,FMT=*)&
        log10(min_rsv(l)) +&
        rsv_binwidth(l)*(Nbins+0.5d0),&
        zero_array

close(filechannel1)
end do

deallocate(Ninterpolation_binning,rsv_binning,zero_array)

end subroutine processMultipleInterpolationFiles


subroutine plotInterpolationOccurenceHeatmap()
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!COLLECTIVE VARIABLES TO FILTER DATA
real,dimension(Nvar) :: tmpvals,var_minvar
integer,dimension(Nvar) :: max_val_bin, val_bins
integer :: max_bin

real(dp),dimension(9) :: rsv, min_rsv, max_rsv
integer :: min_Ninterpolation,max_Ninterpolation,Ninterpolation

character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5

!I/O HANDLING
integer :: iostate

!HISTOGRAM VARIABLES
real(dp),allocatable :: heatmap_binning(:,:)

!INTEGER INCREMENTALS
integer :: n,i,j

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

var_minvar = (/ 2.0, 2.0 /)

max_val_bin = int((var_maxvar - var_minvar)/var_spacing)

allocate(heatmap_binning(max_val_bin(1),&
                         max_val_bin(2)))
heatmap_binning = 0

open(filechannel1,file=gridpath5//&
        "truncated"//interpolationfile)
read(filechannel1,FMT="(1x,3(I9,1x),2(I6,1x),"//&
        "2(9(G16.6,1x)))",iostat=iostate) &
        i, j, n,&
        min_Ninterpolation, max_Ninterpolation,&
        min_rsv,max_rsv

do
    read(filechannel1,FMT=*,iostat=iostate)&
            tmpvals, Ninterpolation, rsv
    if (iostate /= 0) exit

    !Here we define what is NOT an occurence
    if (rsv(6) <= max_rsv(5)) cycle

    val_bins = floor((tmpvals-var_minvar)/var_spacing)

    do n = 1, Nvar
        if (val_bins(n) == 0) val_bins(n) = 1
        if (val_bins(n) > max_val_bin(n)) &
                val_bins(n) = max_val_bin(n)
    end do

    heatmap_binning(val_bins(1),val_bins(2)) = &
        heatmap_binning(val_bins(1),val_bins(2)) + 1

end do
close(filechannel1)

max_bin = maxval(heatmap_binning)

open(filechannel1,file=gridpath5//temporaryfile1)
do i = 1, max_val_bin(1)
    tmpvals(1) = var_minvar(1) + (i-1)*var_spacing(1)
    
    do j = 1, max_val_bin(2)
        tmpvals(2) = var_minvar(2) + (j-1)*var_spacing(2)
        
        write(filechannel1,FMT=*) tmpvals, &
                heatmap_binning(i,j)
    end do
    
    write(filechannel1,FMT=*) "" 
end do
close(filechannel1)

deallocate(heatmap_binning)



open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set terminal pngcairo size 1800,1800'
write(gnuplotchannel,*) 'set output "'//gridpath4//'InterpolationOccurenceHeatMap.png"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,FMT='(A,4(I0.5,A))') &
        'set palette defined (',&
        0, ' "green", ',&
        1, ' "blue", ',&
        max(max_bin/2,2), ' "yellow", ',&
        max(max_bin,3), ' "red")'
write(gnuplotchannel,FMT='(A,I0.5,A)') 'set cbrange [0:',max_bin,']'
write(gnuplotchannel,*) 'set cblabel "Number of Occurences" font ",18" offset 1,0'

write(gnuplotchannel,*) 'set title "Interpolation Occurence Heatmap of an H_2 - H_2 System" font ",32" offset 0,3'
write(gnuplotchannel,*) 'set xlabel "Var1 (A)" font ",28" offset 0,-2'
write(gnuplotchannel,*) 'set xtics 1 font ",24"'
write(gnuplotchannel,*) 'set xrange [', var_minvar(1), ':', var_maxvar(1),']'
write(gnuplotchannel,*) 'set ylabel "Var2 (A)" font ",28" offset -5,0'
write(gnuplotchannel,*) 'set ytics 1 font ",24"'
write(gnuplotchannel,*) 'set yrange [', var_minvar(2), ':',var_maxvar(2),']'
write(gnuplotchannel,*) 'set cbtics'
write(gnuplotchannel,*) 'set view map'
write(gnuplotchannel,*) 'set pm3d interpolate 1,1'
write(gnuplotchannel,*) 'splot "'//gridpath5//temporaryfile1//'" u 1:2:3 w pm3d'
close(gnuplotchannel)
call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

end subroutine plotInterpolationOccurenceHeatmap



subroutine getRMSDinterpolation2(PNGfilename)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!FORMAT OF PNG FILES TO BE MADE
character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: PNGfilename

integer, dimension(comparison_number,3) :: counters

integer :: iostate

!INTEGER INCREMENTALS
integer :: n

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

call processMultipleInterpolationFiles(counters)

if (any(counters(:,3) == 0)) then
    return
end if

open(gnuplotchannel,file=gridpath5//PNGfilename//"_"//gnuplotfile)
write(gnuplotchannel,FMT="(A,I6)") 'set term pngcairo size 2400,',1200*comparison_number
write(gnuplotchannel,*) 'set output ARG1."/../'//PNGfilename//'.png"'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
write(gnuplotchannel,FMT="(A,I0.2,A)") 'set multiplot layout ',comparison_number,&
                        ',2 columnsfirst margins 0.1,0.95,.1,.9 spacing 0.1,0 title '//&
                        '"Interpolation Error with Varying Threshold" font ",36" offset 0,3'
!write(gnuplotchannel,FMT="(A,F7.3,',',F7.3,A)") &
!        'set label 1 "Vals = (',vals(1),vals(2),')" at screen 0.3,0.950'
!write(gnuplotchannel,FMT="(A,I6,A)") &
!        'set label 2 "Ntraj = ', Ntraj_max, '" at screen 0.3,0.940'
!write(gnuplotchannel,FMT='(A,F9.4,A)') 'set label 3 "AlphaRatio = ',alpha_ratio, &
!        '" at screen 0.3,0.930'
write(gnuplotchannel,*) 'set ylabel "Frequency" font ",18"'
write(gnuplotchannel,*) 'set xtics nomirror'
write(gnuplotchannel,*) 'set grid xtics lw 2'
write(gnuplotchannel,*) 'unset xlabel'
write(gnuplotchannel,*) 'set format x ""'

if (comparison_number > 1) then
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '# Edit to title this plot (01)'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A,I0.2,A)") &
        'plot ARG1."/'//&
        'InterpolationTDD_binning.dat" u '//'1:',2,' w boxes t \'
write(gnuplotchannel,FMT="(A)") &
        '" "\'
write(gnuplotchannel,FMT="(A,F9.6,A)") &
        ' fs transparent solid ',1.0d0,' noborder'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'

do n = 2, comparison_number-1

write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A,I0.2,A)") '# Edit to title this plot (',n,')'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A,I0.2,A)") &
        'plot ARG1."/'//&
        'InterpolationTDD_binning.dat" u '//'1:',n+1,' w boxes t \'
write(gnuplotchannel,FMT="(A)") &
        '" "\'
write(gnuplotchannel,FMT="(A,F9.6,A)") &
        ' fs transparent solid ',1.0d0,' noborder'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'

end do
end if

write(gnuplotchannel,*) 'unset xtics'
write(gnuplotchannel,*) 'set xtics out nomirror'
write(gnuplotchannel,*) 'set format x'
write(gnuplotchannel,*) 'set xlabel "Number of Points Below Threshold (Ninterpolation)" font ",24"'

write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A,I0.2,A)") '# Edit to title this plot (',comparison_number,')'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A,I0.2,A)") &
        'plot ARG1."/'//&
        'InterpolationTDD_binning.dat" u '//'1:',comparison_number+1,' w boxes t \'
write(gnuplotchannel,FMT="(A)") &
        '" "\'
write(gnuplotchannel,FMT="(A,F9.6,A)") &
        ' fs transparent solid ',1.0d0,' noborder'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'

write(gnuplotchannel,*) 'unset xlabel'
write(gnuplotchannel,*) 'set xtics nomirror'
write(gnuplotchannel,*) 'set xtics ('//&
                                         '"" log10(.00000001), '//&
                                         '"" log10(.00000005), '//&
                                          '"" log10(.0000001), '//&
                                          '"" log10(.0000005), '//&
                                           '"" log10(.000001), '//&
                                           '"" log10(.000005), '//&
                                           '""  log10(.00001), '//&
                                           '""  log10(.00005), '//&
                                           '""   log10(.0001), '//&
                                           '""   log10(.0005), '//&
                                           '""    log10(.001), '//&
                                           '""    log10(.005), '//&
                                           '""     log10(.01), '//&
                                           '""     log10(.05), '//&
                                           '""      log10(.1), '//&
                                           '""      log10(.5), '//&
                                           ' ""       log10(1), '//&
                                           ' ""       log10(5), '//&
                                           ' ""      log10(10), '//&
                                           ' ""      log10(50), '//&
                                           ' ""     log10(100), '//&
                                           ' ""     log10(500), '//&
                                           ' ""    log10(1000), '//&
                                           ' ""    log10(5000), '//&
                                           ' ""   log10(10000), '//&
                                   ')'

if (comparison_number > 1) then
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        'plot ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',2,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',3,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',4,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder'

do n = 2, comparison_number-1

write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        'plot ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',3*(n-1)+2,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',3*(n-1)+3,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',3*(n-1)+4,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder'

end do
end if

write(gnuplotchannel,*) 'set xtics out nomirror'
write(gnuplotchannel,*) 'set xtics ('//&
                                         '"1e-8" log10(.00000001), '//&
                                         '"5e-8" log10(.00000005), '//&
                                          '"1e-7" log10(.0000001), '//&
                                          '"5e-7" log10(.0000005), '//&
                                           '"1e-6" log10(.000001), '//&
                                           '"5e-6" log10(.000005), '//&
                                           '"1e-5"  log10(.00001), '//&
                                           '"5e-5"  log10(.00005), '//&
                                           '"1e-4"   log10(.0001), '//&
                                           '"5e-4"   log10(.0005), '//&
                                           '"1e-3"    log10(.001), '//&
                                           '"5e-3"    log10(.005), '//&
                                           '"1e-2"     log10(.01), '//&
                                           '"5e-2"     log10(.05), '//&
                                           '"1e-1"      log10(.1), '//&
                                           '"5e-1"      log10(.5), '//&
                                           ' "1e0"       log10(1), '//&
                                           ' "5e0"       log10(5), '//&
                                           ' "1e1"      log10(10), '//&
                                           ' "5e1"      log10(50), '//&
                                           ' "1e2"     log10(100), '//&
                                           ' "5e2"     log10(500), '//&
                                           ' "1e3"    log10(1000), '//&
                                           ' "5e3"    log10(5000), '//&
                                           ' "1e4"   log10(10000), '//&
                                   ')'

write(gnuplotchannel,*) 'set xlabel "Relative Error of Accept Best to Interpolation" font ",24"'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        'plot ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',3*comparison_number-1,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',3*comparison_number,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/'//&
        'InterpolationRED_binning.dat" u '//'1:',3*comparison_number+1,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder'

close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot -c "//&
        gridpath5//PNGfilename//"_"//gnuplotfile//" "//&
        '"'//gridpath5(1:gridpath_length+expfolder_length+5-1)//'"')

end subroutine getRMSDinterpolation2




subroutine getInterpolationplot(PNGfilename)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!FORMAT OF PNG FILES TO BE MADE
character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: PNGfilename

integer, dimension(comparison_number,3) :: counters
real(dp),dimension(comparison_number) :: var_max,var_Amean,var_Gmean
real(dp) :: lowerlimit, upperlimit

integer :: iostate

!INTEGER INCREMENTALS
integer :: n

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

call processMultipleInterpolationFiles(counters)

if (all(counters == 0)) then
    return
end if

open(filechannel1,file=gridpath5//&
        trim(adjustl(comparison_SATRVname))//&
        "_binning.dat")
read(filechannel1,FMT="(1x,2(E16.6,1x))") &
    lowerlimit, upperlimit

do n = 1, comparison_number
    read(filechannel1,FMT="(1x,3(E16.6,1x))") &
        var_max(n), var_Amean(n), var_Gmean(n)
end do
close(filechannel1)

open(gnuplotchannel,file=gridpath5//gnuplotfile//"_"//&
        trim(adjustl(comparison_SATRVname)))
write(gnuplotchannel,FMT="(A,I6)") 'set term pngcairo size 1200,',400*comparison_number
write(gnuplotchannel,*) 'set encoding utf8'
write(gnuplotchannel,*) 'set output ARG1."/../'//PNGfilename//'.png"'
write(gnuplotchannel,*) 'PL=strlen(ARG1)'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
write(gnuplotchannel,FMT="(A,I0.2,A)") 'set multiplot layout ',comparison_number,&
                        ',1 columnsfirst margins 0.1,0.95,.1,.9 spacing 0.1,0 title '//&
                        '"" font ",32" offset 0,-3'
write(gnuplotchannel,*) 'set style histogram clustered gap 1'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'set ylabel "Occurence" font ",24"'
!write(gnuplotchannel,*) 'set ytics font ",16"'
write(gnuplotchannel,*) 'set format y ""'
write(gnuplotchannel,*) 'set xtics nomirror'
write(gnuplotchannel,*) 'set grid xtics lw 2'

write(gnuplotchannel,*) 'unset xlabel'
write(gnuplotchannel,*) 'set xtics nomirror font ",16"'

if (trim(adjustl(comparison_SATRVname)) /= 'InterpolationTDD') then
    write(gnuplotchannel,*) 'set xtics ('//&
                                         '"" log10(.00000001), '//&
                                         '"" log10(.00000005), '//&
                                          '"" log10(.0000001), '//&
                                          '"" log10(.0000005), '//&
                                           '"" log10(.000001), '//&
                                           '"" log10(.000005), '//&
                                           '""  log10(.00001), '//&
                                           '""  log10(.00005), '//&
                                           '""   log10(.0001), '//&
                                           '""   log10(.0005), '//&
                                           '""    log10(.001), '//&
                                           '""    log10(.005), '//&
                                           '""     log10(.01), '//&
                                           '""     log10(.05), '//&
                                           '""      log10(.1), '//&
                                           '""      log10(.5), '//&
                                           ' ""       log10(1), '//&
                                           ' ""       log10(5), '//&
                                           ' ""      log10(10), '//&
                                           ' ""      log10(50), '//&
                                           ' ""     log10(100), '//&
                                           ' ""     log10(500), '//&
                                           ' ""    log10(1000), '//&
                                           ' ""    log10(5000), '//&
                                           ' ""   log10(10000), '//&
                                   ')'
end if

write(gnuplotchannel,FMT="(A,E16.8,A,E16.8,A)")&
        'set xrange [',lowerlimit,':',upperlimit,']'

if (comparison_number > 1) then
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '# Edit to title this plot (01)'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") 'unset arrow'
write(gnuplotchannel,FMT="(A)") 'set format x ""'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_max(1),&
        ',graph 0 to ', var_max(1), ', graph 1 nohead front '//&
        'lw 2 lc rgb "black"'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_Amean(1),&
        ',graph 0 to ', var_Amean(1), ', graph 1 nohead front '//&
        'lw 2 lc rgb "blue"'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_Gmean(1),&
        ',graph 0 to ', var_Gmean(1), ', graph 1 nohead front '//&
        'lw 2 lc rgb "green"'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        'plot ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',2,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',3,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',4,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'

do n = 2, comparison_number-1

write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A,I0.2,A)") '# Edit to title this plot (',n,')'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") 'unset arrow'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_max(n),&
        ',graph 0 to ', var_max(n), ', graph 1 nohead front '//&
        'lw 2 lc rgb "black"'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_Amean(n),&
        ',graph 0 to ', var_Amean(n), ', graph 1 nohead front '//&
        'lw 2 lc rgb "blue"'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_Gmean(n),&
        ',graph 0 to ', var_Gmean(n), ', graph 1 nohead front '//&
        'lw 2 lc rgb "green"'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        'plot ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',3*(n-1)+2,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',3*(n-1)+3,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',3*(n-1)+4,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'

end do
end if

write(gnuplotchannel,*) 'set xtics out nomirror font ",16"'
if (trim(adjustl(comparison_SATRVname)) /= 'InterpolationTDD') then
    write(gnuplotchannel,*) 'set xtics ('//&
                                         '"1e-8" log10(.00000001), '//&
                                         '"5e-8" log10(.00000005), '//&
                                          '"1e-7" log10(.0000001), '//&
                                          '"5e-7" log10(.0000005), '//&
                                           '"1e-6" log10(.000001), '//&
                                           '"5e-6" log10(.000005), '//&
                                           '"1e-5"  log10(.00001), '//&
                                           '"5e-5"  log10(.00005), '//&
                                           '"1e-4"   log10(.0001), '//&
                                           '"5e-4"   log10(.0005), '//&
                                           '"1e-3"    log10(.001), '//&
                                           '"5e-3"    log10(.005), '//&
                                           '"1e-2"     log10(.01), '//&
                                           '"5e-2"     log10(.05), '//&
                                           '"1e-1"      log10(.1), '//&
                                           '"5e-1"      log10(.5), '//&
                                           ' "1e0"       log10(1), '//&
                                           ' "5e0"       log10(5), '//&
                                           ' "1e1"      log10(10), '//&
                                           ' "5e1"      log10(50), '//&
                                           ' "1e2"     log10(100), '//&
                                           ' "5e2"     log10(500), '//&
                                           ' "1e3"    log10(1000), '//&
                                           ' "5e3"    log10(5000), '//&
                                           ' "1e4"   log10(10000), '//&
                                   ')'
end if

write(gnuplotchannel,*) 'set xlabel "" font ",24"'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A,I0.2,A)") '# Edit to title this plot (',comparison_number,')'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") 'unset arrow'
write(gnuplotchannel,FMT="(A)") 'set format x "% g"'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_max(comparison_number),&
        ',graph 0 to ', var_max(comparison_number), ', graph 1 nohead front '//&
        'lw 2 lc rgb "black"'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_Amean(comparison_number),&
        ',graph 0 to ', var_Amean(comparison_number), ', graph 1 nohead front '//&
        'lw 2 lc rgb "blue"'
write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
        'set arrow from ', var_Gmean(comparison_number),&
        ',graph 0 to ', var_Gmean(comparison_number), ', graph 1 nohead front '//&
        'lw 2 lc rgb "green"'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        'plot ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',3*comparison_number-1,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',3*comparison_number,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder,\'
write(gnuplotchannel,FMT="(A,I2,A,F9.6,A)") &
        '     ARG1."/".'//&
        'ARG0[PL+14:]."_binning.dat" u '//'1:',3*comparison_number+1,' w boxes t "" '//&
        'fs transparent solid ',1.0d0,' noborder'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'
write(gnuplotchannel,FMT="(A)") '#'

close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot -c "//&
        gridpath5//gnuplotfile//"_"//trim(adjustl(comparison_SATRVname))//&
        ' "'//gridpath5(1:gridpath_length+expfolder_length+5-1)//'"')

end subroutine getInterpolationplot

subroutine processCheckstateFile()
use PARAMETERS
use FUNCTIONS
use ANALYSIS
implicit none

!VARIABLES IN THE INTERPOLATION FILE
real(dp),dimension(Nvar) :: vals
real(dp) :: min_rmsd, min_rmsd_prime
real(dp) :: U, KE
integer :: number_of_frames, order, neighbor_check

integer :: Naccept = 0
integer :: Nalpha_tries = 1
integer :: rewind_step = 1
logical :: Naccept_end_flag = .false.

integer,allocatable :: success_binning(:,:)
integer,allocatable :: failure_binning(:,:)

character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5

!I/O HANDLING
integer :: iostate

!INTEGER INCREMENTALS
integer :: n

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

allocate(success_binning(Naccept_max,Nalpha_tries_max),&
         failure_binning(Naccept_max,Nalpha_tries_max))
success_binning = 0
failure_binning = 0

call system("rm "//gridpath5//"truncated"//checkstatefile)

open(filechannel1,file=gridpath5//checkstatefile)
open(filechannel2,file=gridpath5//"truncated"//checkstatefile)
do
    read(filechannel1,iostat=iostate,FMT=*)&
            number_of_frames,order,&
            neighbor_check,steps,&
            min_rmsd,min_rmsd_prime,&
            vals(1),vals(2),U,KE
    if (iostate /= 0) exit

    if (Naccept_end_flag) then
        Naccept_end_flag = .false.

        if (steps == rewind_step) then
            failure_binning(Naccept,Nalpha_tries) =&
                failure_binning(Naccept,Nalpha_tries) + 1

            Nalpha_tries = Nalpha_tries + 1
            if (Nalpha_tries == Nalpha_tries_max + 1) then
                rewind_step = steps + 1
            end if

            Naccept = 0
        else if (steps == rewind_step+Naccept+1) then
            success_binning(Naccept,Nalpha_tries) =&
                success_binning(Naccept,Nalpha_tries) + 1

            Nalpha_tries = 1
            Naccept = 0

            close(filechannel2)
            call system("cat "//gridpath5//temporaryfile1//&
                    " >> "//gridpath5//"truncated"//checkstatefile)
            open(filechannel2,&
                    file=gridpath5//"truncated"//checkstatefile,&
                    position="append")

            rewind_step = steps
        else
            print *, "error at step", steps
        end if
    end if

    if (min_rmsd_prime < default_rmsd) then
        if (steps == rewind_step) then
            Naccept = 0

            if (Nalpha_tries == Nalpha_tries_max + 1) then
                Nalpha_tries = 1
            else
                open(filechannel3,file=gridpath5//temporaryfile1)
                write(filechannel3,FMT=*) number_of_frames,order,&
                        neighbor_check,steps,min_rmsd,min_rmsd_prime,&
                        vals(1),vals(2),U,KE
                close(filechannel3)
            end if
        else
            Naccept = Naccept + 1
            if (Naccept == Naccept_max) then
                Naccept_end_flag = .true.
            end if

            open(filechannel3,file=gridpath5//temporaryfile1,&
                 position="append")
            write(filechannel3,FMT=*) number_of_frames,order,&
                    neighbor_check,steps,min_rmsd,min_rmsd_prime,&
                    vals(1),vals(2),U,KE
            close(filechannel3)

        end if
    else
        if (steps == rewind_step) then
            rewind_step = steps + 1

            write(filechannel2,FMT=*) number_of_frames,order,&
                    neighbor_check,steps,min_rmsd,min_rmsd_prime,&
                    vals(1),vals(2),U,KE
        else
            Naccept = Naccept + 1
            Naccept_end_flag = .true.

            open(filechannel3,file=gridpath5//temporaryfile1,&
                 position="append")
            write(filechannel3,FMT=*) number_of_frames,order,&
                    neighbor_check,steps,min_rmsd,min_rmsd_prime,&
                    vals(1),vals(2),U,KE
            close(filechannel3)
        end if
    end if

end do
close(filechannel1)
close(filechannel2)

deallocate(success_binning,failure_binning)

end subroutine processCheckstateFile






end module analyzeRMSDThresholdwithMultipleGrids
