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


subroutine getScatteringAngles1(prefix_filename,JPGfilename)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

!SCATTERING ANGLE BINNING
integer :: Nsamples, Nsamples_max
real :: bin_width, binMean, binSD,binRMSD
real,allocatable :: binAverage(:)
integer,allocatable :: binTotal(:,:),sampleSize(:)
integer :: binTally
real :: binThreshold = 1.0

!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
character(*), intent(in) :: prefix_filename

!FORMAT OF JPG FILES TO BE MADE
character(*), intent(in) :: JPGfilename

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(5) :: variable_length_text1, variable_length_text2
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(6) :: Ntraj_text
character(150) :: old_filename

!SCATTERING ANGLES
real(dp),dimension(3,3) :: coords_initial,velocities_initial,coords_final,velocities_final
real(dp) :: ScatteringAngle_real, TranslationalEnergy_real
real(dp) :: TranslationalEnergy_max
integer :: scatteringAngle, TranslationalEnergy

!HEAT MAP VARIABLES
integer,allocatable :: angle_energy_bins(:,:)
real(dp) :: sizeEnergyBin,sizeAngleBin
integer :: energyBins = 100
integer :: angleBins = 100
integer :: occurence_max

!I/O HANDLING
integer :: iostate

!Incremental Integers
integer :: i, j, k


bin_width = pi / SA_Nbins
Nsamples_max = (Ngrid_total * Ntraj_max) / Ntesttraj
allocate(binTotal(SA_Nbins,Nsamples_max),binAverage(SA_Nbins),sampleSize(Nsamples_max))

!This data was made during creation (or should have been!) so all we need to do
!is merge, read, and plot them
TranslationalEnergy_max = 0.0d0
old_filename = ""
do Ngrid = 1, Ngrid_total

	write(variable_length_text,FMT="(I5)") Ngrid_text_length
        !The folders are named starting from 001 by increments of 1
        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid

        !The plots are named starting from Ntraj_max by increments of Ntraj_max (the number of trajectories)
        write(Ntraj_text,FMT="(I6)") Ngrid * Ntraj_max

	if (trim(adjustl(old_filename)) /= "") then
        	call system("cp "//trim(adjustl(old_filename))//" "//&
                            gridpath0//"/"//Ngrid_text//"/"//trim(adjustl(Ntraj_text))//SAfile)
        	old_filename = gridpath0//"/"//Ngrid_text//"/"//trim(adjustl(Ntraj_text))//SAfile
	else
        	old_filename = gridpath0//"/"//Ngrid_text//"/"//trim(adjustl(Ntraj_text))//SAfile
        	call system("rm "//trim(adjustl(old_filename)))
	end if


	open(filechannel1,file=gridpath0//Ngrid_text//"/"//prefix_filename//SAfile)
	open(filechannel2,file=trim(adjustl(old_filename)),position="append")
	do i = (Ngrid-1)*Ntraj_max+1, Ngrid*Ntraj_max
		read(filechannel1,FMTsa) ScatteringAngle_real, TranslationalEnergy_real

		TranslationalEnergy_max = max(TranslationalEnergy_max,TranslationalEnergy_real)

		scatteringAngle = ceiling(ScatteringAngle_real / bin_width)
                if (scatteringAngle > SA_Nbins) scatteringAngle = SA_Nbins

		write(filechannel2,FMT=*) scatteringAngle
	end do
	close(filechannel1)
	close(filechannel2)


        !This is the gnuplot code to make the plots
        open(gnuplotchannel,file=gridpath0//gnuplotfile)
        write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
        write(gnuplotchannel,*) 'set output "'//gridpath0//trim(adjustl(Ntraj_text))//JPGfilename//'"'
        write(gnuplotchannel,*) 'set title "Scattering Angle Distribution of '//trim(adjustl(Ntraj_text))//&
                                ' Trajectories"'
        write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
        write(gnuplotchannel,*) 'unset key'
	write(gnuplotchannel,*) 'pi = 3.14159265'
	write(variable_length_text,FMT="(F5.4)") bin_width
	write(gnuplotchannel,*) 'box_width = '//variable_length_text
	write(gnuplotchannel,*) 'set boxwidth box_width'
        write(gnuplotchannel,*) 'set xlabel "Scattering Angle"'
        write(gnuplotchannel,*) 'set ylabel "Occurence"'
        write(gnuplotchannel,*) 'set xrange [0:pi]'
        write(gnuplotchannel,*) 'set autoscale y'
        write(gnuplotchannel,*) 'plot "'//gridpath0//Ngrid_text//'/'//&
                                trim(adjustl(Ntraj_text))//SAfile//&
                                '" u (box_width*($1-0.5)):(1.0) smooth frequency with boxes'
        close(gnuplotchannel)

        !And then we just input it into gnuplot.exe
        call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)

end do


!Now we want a reference distribution
!This will be the distribution of ALL trajectories
binTotal = 0
sampleSize = 0
open(filechannel1,file=gridpath0//Ngrid_text//"/"//trim(adjustl(Ntraj_text))//SAfile,action="read")
do Nsamples = 1, Nsamples_max
        do i = 1, Ntesttraj
                read(filechannel1,FMT=*) scatteringAngle
                binTotal(scatteringAngle,Nsamples) = binTotal(scatteringAngle,Nsamples) + 1
        end do
        sampleSize(Nsamples) = Nsamples
end do
close(filechannel1)

do i = 1,SA_Nbins
        binAverage(i) = sum(binTotal(i,:)) * 1.0 / Nsamples_max
end do

!Now we want error bars
!For this, we will get a running average distribution by adding Ntesttraj trajectories at a time
!When the running average converges then we stop and see the variance
binTally = 0
open(filechannel1,file=gridpath0//"RMSD"//cumulativefile//".dat",action="write")
do Nsamples = 1, Nsamples_max
        binRMSD = 0.0
        do i = 1, SA_Nbins
                binRMSD = binRMSD + (sum(binTotal(i,1:Nsamples)) * 1.0 / Nsamples - binAverage(i))**2
        end do

        binRMSD = sqrt(binRMSD/SA_Nbins)
        write(filechannel1,FMT=*) Nsamples*Ntesttraj, binRMSD, binThreshold

        if (binRMSD < binThreshold) then
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

!Now we must do the laborious job of binning these
!Each bin has its own average and standard deviation based on how many samples we took
open(filechannel1,file=gridpath0//"Adjusted"//cumulativefile//".dat")
do i = 1, SA_Nbins
        binMean = binAverage(i)
        binSD = sqrt(sum((binTotal(i,1:Nsamples) - binMean)**2)/(Nsamples - 1))

        write(filechannel1,*) (i-0.5)*bin_width, binMean, binSD!/sqrt(real(Nsamples))
end do
close(filechannel1)

deallocate(binAverage,binTotal,sampleSize)

!We have everything we need to draw the distribution with error bars
open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) "set terminal jpeg size 1200,1200"
write(variable_length_text,FMT="(I5)") Ntesttraj
write(gnuplotchannel,*) 'set output "'//gridpath0//'Adjustedto'//&
                        trim(adjustl(variable_length_text))//JPGfilename//'"'
write(gnuplotchannel,*) 'set multiplot layout 2,1'
write(variable_length_text,FMT="(I5)") Ntesttraj
write(gnuplotchannel,*) 'set title "Convergence of the '//trim(adjustl(variable_length_text))//' Scattering Angle Distribution"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set xlabel "Number of Trajectories"'
write(gnuplotchannel,*) 'set ylabel "RMSD of the Distribution Over a Running Average"'
write(gnuplotchannel,*) 'set autoscale x'
write(variable_length_text,FMT="(I5)") Ntesttraj
write(gnuplotchannel,*) 'set xtics '//trim(adjustl(variable_length_text))
write(gnuplotchannel,*) 'set autoscale y'
write(variable_length_text1,FMT="(I5)") Ntesttraj
write(variable_length_text2,FMT="(F5.3)") binThreshold
write(gnuplotchannel,*) 'set label 1 "Convergence Threshold" at first '//&
                        variable_length_text1//','//variable_length_text2
write(gnuplotchannel,*) 'plot "'//gridpath0//'RMSD'//cumulativefile//'.dat" u 1:2 w lines, \'
write(gnuplotchannel,*) '     "'//gridpath0//'RMSD'//cumulativefile//'.dat" u 1:3 w lines lc -1'
write(variable_length_text,FMT="(I5)") Ntesttraj
write(gnuplotchannel,*) 'set title "Final Scattering Angle Distribution for '//&
                        trim(adjustl(variable_length_text))//' Trajectories"'
write(gnuplotchannel,*) 'pi = 3.14159265'
write(gnuplotchannel,*) 'set xrange [0:pi]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set xlabel "Scattering Angle (rad)"'
write(gnuplotchannel,*) 'set ylabel "Frequency"'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'unset label 1'
write(gnuplotchannel,*) 'plot "'//gridpath0//'Adjusted'//cumulativefile//'.dat" u 1:2 w boxes, \'
write(gnuplotchannel,*) '     "'//gridpath0//'Adjusted'//cumulativefile//'.dat" u 1:2:3 w yerrorbars'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)





allocate(angle_energy_bins(angleBins,energyBins))
angle_energy_bins = 0
sizeEnergyBin = TranslationalEnergy_max / (energyBins)
sizeAngleBin = pi2 / (angleBins)
do Ngrid = 1, Ngrid_total
        write(variable_length_text,FMT="(I5)") Ngrid_text_length
        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
        open(filechannel1,file=gridpath0//Ngrid_text//"/"//prefix_filename//SAfile)
        
        do j = 1, Ntraj_max
        
                read(filechannel1,FMTsa) ScatteringAngle_real, TranslationalEnergy_real
        
        	ScatteringAngle = ceiling(ScatteringAngle_real / sizeAngleBin)
        	TranslationalEnergy = ceiling((TranslationalEnergy_real)/ sizeEnergyBin)
        
        	if (ScatteringAngle > angleBins) ScatteringAngle = angleBins
        	if (ScatteringAngle == 0) ScatteringAngle = 1
        
        	if (TranslationalEnergy > energyBins) TranslationalEnergy = energyBins
        	if (TranslationalEnergy == 0) TranslationalEnergy = 1
        
        	angle_energy_bins(ScatteringAngle,TranslationalEnergy) = &
                           angle_energy_bins(ScatteringAngle,TranslationalEnergy) + 1
        end do
        close(filechannel1)
end do

occurence_max = maxval(angle_energy_bins)

open(filechannel1,file=gridpath0//temporaryfile1)
do i = 1, angleBins
	do j = 1, energyBins
		write(filechannel1,FMT=*) i*sizeAngleBin, j*sizeEnergyBin, angle_energy_bins(i,j)
	end do
	write(filechannel1,FMT=*) ""
end do
close(filechannel1)

!This is the gnuplot code to make the plots
open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath0//'HeatMap_'//trim(adjustl(Ntraj_text))//&
                        prefix_filename//JPGfilename//'"'
write(gnuplotchannel,*) 'set title "Scattering Angle Distribution" offset 0,2'
write(gnuplotchannel,*) 'set lmargin at screen 0.05'
write(gnuplotchannel,*) 'set rmargin at screen 0.85'
write(gnuplotchannel,*) 'set bmargin at screen 0.1'
write(gnuplotchannel,*) 'set tmargin at screen 0.9'
write(gnuplotchannel,*) 'set pm3d map'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set multiplot'
write(gnuplotchannel,*) 'set parametric'
write(gnuplotchannel,*) 'set isosamples 500'
write(gnuplotchannel,*) 'unset border'
write(gnuplotchannel,*) 'unset xtics'
write(gnuplotchannel,*) 'unset ytics'
write(gnuplotchannel,*) 'set angles radian'
write(gnuplotchannel,*) 'r = ', TranslationalEnergy_max
write(gnuplotchannel,*) 'rmax = r*1.1'
write(gnuplotchannel,*) 'pi = 3.14159'
!write(gnuplotchannel,*) 'set palette rgb 7, 5, 15'
write(variable_length_text,FMT="(I5)") occurence_max
write(gnuplotchannel,*) 'set cbrange [0:'//trim(adjustl(variable_length_text))//'/2]'
write(gnuplotchannel,*) 'set cblabel "Frequency"'
!write(gnuplotchannel,*) 'set palette positive nops_allcF maxcolors 0 gamma 1.5 color model XYZ '
!write(gnuplotchannel,*) 'set palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV'
!write(gnuplotchannel,*) 'set palette rgbformulae 3, 2, 2'
write(gnuplotchannel,*) 'set palette defined ( 0 0 1 0, 0.3333 0 0 1, 0.6667 1 0 0,\'
write(gnuplotchannel,*) '     1 1 0.6471 0 )'
write(gnuplotchannel,*) 'scaling_factor = 1.0'
write(gnuplotchannel,*) 'set xrange[-rmax:rmax]'
write(gnuplotchannel,*) 'set yrange[-rmax:rmax]'
write(gnuplotchannel,*) 'set colorbox user origin 0.9,0.1 size 0.03,0.8'
write(gnuplotchannel,*) 'splot "'//gridpath0//temporaryfile1//'" u ($2*cos($1)):($2*sin($1)):3'
write(gnuplotchannel,*) 'set style line 11 lc rgb "black" lw 2'
write(gnuplotchannel,*) 'set grid polar ls 11'
write(gnuplotchannel,*) 'set polar'
write(gnuplotchannel,*) 'set rrange[0:rmax]'
!write(gnuplotchannel,*) 'unset raxis'
!write(gnuplotchannel,*) 'set rlabel "Translational Energy Change (eV)"'
!write(gnuplotchannel,*) 'set rtics format "" scale 0'
write(gnuplotchannel,*) 'set rtics r / 4'
write(gnuplotchannel,*) 'unset parametric'
write(gnuplotchannel,*) 'set for [i=0:330:30] label at first (rmax*scaling_factor)*cos(i*pi/180),'//&
                        ' first (rmax*scaling_factor)*sin(i*pi/180)\'
write(gnuplotchannel,*) 'center sprintf(''%d'', i)'
write(gnuplotchannel,*) 'plot NaN w l'
write(gnuplotchannel,*) 'unset multiplot'
close(gnuplotchannel)

!And then we just input it into gnuplot.exe
call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)









allocate(angle_energy_bins(angleBins,energyBins))
angle_energy_bins = 0
sizeEnergyBin = TranslationalEnergy_max / (energyBins)
sizeAngleBin = pi2 / (angleBins)
do Ngrid = 1, Ngrid_total
        write(variable_length_text,FMT="(I5)") Ngrid_text_length
        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
        open(filechannel1,file=gridpath0//Ngrid_text//"/"//prefix_filename//SAfile)
        
        do j = 1, Ntraj_max
        
                read(filechannel1,FMTsa) ScatteringAngle_real, TranslationalEnergy_real
        
        	ScatteringAngle = ceiling(ScatteringAngle_real / sizeAngleBin)
        	TranslationalEnergy = ceiling((TranslationalEnergy_real)/ sizeEnergyBin)
        
        	if (ScatteringAngle > angleBins) ScatteringAngle = angleBins
        	if (ScatteringAngle == 0) ScatteringAngle = 1
        
        	if (TranslationalEnergy > energyBins) TranslationalEnergy = energyBins
        	if (TranslationalEnergy == 0) TranslationalEnergy = 1
        
        	angle_energy_bins(ScatteringAngle,TranslationalEnergy) = &
                           angle_energy_bins(ScatteringAngle,TranslationalEnergy) + 1
        end do
        close(filechannel1)
end do

occurence_max = maxval(angle_energy_bins)

open(filechannel1,file=gridpath0//temporaryfile1)
do i = 1, angleBins
	do j = 1, energyBins
		write(filechannel1,FMT=*) i*sizeAngleBin, j*sizeEnergyBin, angle_energy_bins(i,j)
	end do
	write(filechannel1,FMT=*) ""
end do
close(filechannel1)

!This is the gnuplot code to make the plots
open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath0//'HeatMap_'//trim(adjustl(Ntraj_text))//&
                        prefix_filename//JPGfilename//'"'
write(gnuplotchannel,*) 'set title "Scattering Angle Distribution" offset 0,2'
write(gnuplotchannel,*) 'set lmargin at screen 0.05'
write(gnuplotchannel,*) 'set rmargin at screen 0.85'
write(gnuplotchannel,*) 'set bmargin at screen 0.1'
write(gnuplotchannel,*) 'set tmargin at screen 0.9'
write(gnuplotchannel,*) 'set pm3d map'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set multiplot'
write(gnuplotchannel,*) 'set parametric'
write(gnuplotchannel,*) 'set isosamples 500'
write(gnuplotchannel,*) 'unset border'
write(gnuplotchannel,*) 'unset xtics'
write(gnuplotchannel,*) 'unset ytics'
write(gnuplotchannel,*) 'set angles radian'
write(gnuplotchannel,*) 'r = ', TranslationalEnergy_max
write(gnuplotchannel,*) 'rmax = r*1.1'
write(gnuplotchannel,*) 'pi = 3.14159'
!write(gnuplotchannel,*) 'set palette rgb 7, 5, 15'
write(variable_length_text,FMT="(I5)") occurence_max
write(gnuplotchannel,*) 'set cbrange [0:'//trim(adjustl(variable_length_text))//'/2]'
write(gnuplotchannel,*) 'set cblabel "Frequency"'
!write(gnuplotchannel,*) 'set palette positive nops_allcF maxcolors 0 gamma 1.5 color model XYZ '
!write(gnuplotchannel,*) 'set palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV'
!write(gnuplotchannel,*) 'set palette rgbformulae 3, 2, 2'
write(gnuplotchannel,*) 'set palette defined ( 0 0 1 0, 0.3333 0 0 1, 0.6667 1 0 0,\'
write(gnuplotchannel,*) '     1 1 0.6471 0 )'
write(gnuplotchannel,*) 'scaling_factor = 1.0'
write(gnuplotchannel,*) 'set xrange[-rmax:rmax]'
write(gnuplotchannel,*) 'set yrange[-rmax:rmax]'
write(gnuplotchannel,*) 'set colorbox user origin 0.9,0.1 size 0.03,0.8'
write(gnuplotchannel,*) 'splot "'//gridpath0//temporaryfile1//'" u ($2*cos($1)):($2*sin($1)):3'
write(gnuplotchannel,*) 'set style line 11 lc rgb "black" lw 2'
write(gnuplotchannel,*) 'set grid polar ls 11'
write(gnuplotchannel,*) 'set polar'
write(gnuplotchannel,*) 'set rrange[0:rmax]'
!write(gnuplotchannel,*) 'unset raxis'
!write(gnuplotchannel,*) 'set rlabel "Translational Energy Change (eV)"'
!write(gnuplotchannel,*) 'set rtics format "" scale 0'
write(gnuplotchannel,*) 'set rtics r / 4'
write(gnuplotchannel,*) 'unset parametric'
write(gnuplotchannel,*) 'set for [i=0:330:30] label at first (rmax*scaling_factor)*cos(i*pi/180),'//&
                        ' first (rmax*scaling_factor)*sin(i*pi/180)\'
write(gnuplotchannel,*) 'center sprintf(''%d'', i)'
write(gnuplotchannel,*) 'plot NaN w l'
write(gnuplotchannel,*) 'unset multiplot'
close(gnuplotchannel)

!And then we just input it into gnuplot.exe
call system("gnuplot < "//gridpath0//gnuplotfile)





end subroutine getScatteringAngles1


subroutine getScatteringAngles2(prefix_filename,JPGfilename)
use PARAMETERS
use ANALYSIS
implicit none

!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
character(*), intent(in) :: prefix_filename

!FORMAT OF JPG FILES TO BE MADE
character(*), intent(in) :: JPGfilename

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(6) :: Ntraj_text
character(6) :: boxwidth_text
logical :: grid_is_done

inquire(file=gridpath0//'Adjusted'//cumulativefile//'.dat',exist=grid_is_done)

!This is the gnuplot code to make the plots
open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath0//JPGfilename//'"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'pi = 3.14159265'
write(variable_length_text,FMT="(I0.5)") SA_Nbins
write(gnuplotchannel,*) 'box_width = pi / '//variable_length_text
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'bin_number(x) = floor(x/box_width)'
write(gnuplotchannel,*) 'rounded(x) = box_width * (bin_number(x) + 0.5)'
write(gnuplotchannel,*) 'set xrange [0:pi]'
write(gnuplotchannel,*) 'unset xlabel'
write(gnuplotchannel,*) 'unset xtics'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
write(Ntraj_text,FMT="(I6)") Ntraj
write(gnuplotchannel,*) 'set multiplot layout 3,1 margins 0.15,0.95,.1,.9 spacing 0,0 title '//&
                        '"Angle Distribution of '//trim(adjustl(Ntraj_text))//' Trajectories"'
write(gnuplotchannel,*) 'set style histogram clustered gap 1'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'set ylabel "Scattering Angle Occurence"'
write(gnuplotchannel,*) 'set yrange [0:]'

write(Ntraj_text,FMT="(I6)") Ntesttraj
if (grid_is_done) then
        write(gnuplotchannel,*) 'plot "'//gridpath0//trim(adjustl(Ntraj_text))//SAfile//&
                                '" u (box_width*($1-0.5)):(1.0) smooth frequency w boxes, \'
        write(gnuplotchannel,*) '     "'//gridpath0//'Adjusted'//cumulativefile//'.dat" u 1:2 w boxes'//&
                                       ' fs transparent solid 0.5 noborder, \'
        write(gnuplotchannel,*) '     "'//gridpath0//'Adjusted'//cumulativefile//'.dat" u 1:2:3 w yerrorbars'
else
        write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//SAfile//&
                                '" u (rounded($1)):(1.0) smooth frequency w boxes'
end if

write(gnuplotchannel,*) 'set ylabel "Initial H2 Theta Occurence"'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//initialfile//&
                        '" u (rounded(abs($4'//&
			'))):(1.0) smooth frequency with boxes'
write(gnuplotchannel,*) 'set xlabel "Angle (rad)"'
write(gnuplotchannel,*) 'set xtics'
write(gnuplotchannel,*) 'set ylabel "Initial H2 Phi Occurence"'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//initialfile//&
			'" u (rounded($5'//&
			')):(1.0) smooth frequency with boxes'

close(gnuplotchannel)

!And then we just input it into gnuplot.exe
call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)

end subroutine getScatteringAngles2
















subroutine getTRVimages(prefix_filename,JPGfilename)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
use PHYSICS
implicit none

!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
character(*), intent(in) :: prefix_filename
real,dimension(3) :: TRVenergies
integer,dimension(4) :: TRVenergies_rounded

!BOUNDS FOR THE ENERGY AND ANGULAR DATA
real(dp) :: translational_max, rotational_max, vibrational_max, rovibrational_max
real(dp) :: energy_max
real(dp) :: bin_width

!FORMAT OF JPG FILES TO BE MADE
character(*), intent(in) :: JPGfilename

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(5) :: variable_length_text1, variable_length_text2
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(6) :: Ntraj_text
character(150) :: old_filename

!I/O HANDLING
integer :: iostate

!Incremental Integers
integer :: i, j, k

rotational_max = 0.0d0
translational_max = 0.0d0
vibrational_max = 0.0d0
rovibrational_max = 0.0d0
old_filename = "" 

!This data was made during creation (or should have been!) so all we need to do
!is merge, read, and plot them
do Ngrid = 1, Ngrid_total

	write(variable_length_text,FMT="(I5)") Ngrid_text_length
        !The folders are named starting from 001 by increments of 1
        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid

        !The plots are named starting from Ntraj_max by increments of Ntraj_max (the number of trajectories)
        write(Ntraj_text,FMT="(I6)") Ngrid * Ntraj_max

	open(filechannel1,file=gridpath0//Ngrid_text//"/"//prefix_filename//TRVfile)
	do i = (Ngrid-1)*Ntraj_max+1, Ngrid*Ntraj_max
		read(filechannel1,FMTtrv) (TRVenergies(j),j=1,3)

		translational_max = max(translational_max,TRVenergies(1))
		rotational_max = max(rotational_max,TRVenergies(2))
		vibrational_max = max(vibrational_max,TRVenergies(3))
		rovibrational_max = max(rovibrational_max,TRVenergies(2)+TRVenergies(3))
	end do
	close(filechannel1)

	!Assume that the minimum value is zero
	!(Although that may not be the case for future initial conditions!)
	energy_max = max(rotational_max, translational_max, vibrational_max,rovibrational_max)
	bin_width = energy_max / (TRV_Nbins)


	if (trim(adjustl(old_filename)) /= "") then
        	call system("cp "//trim(adjustl(old_filename))//" "//&
                            gridpath0//"/"//Ngrid_text//"/"//trim(adjustl(Ntraj_text))//TRVfile)
        	old_filename = gridpath0//"/"//Ngrid_text//"/"//trim(adjustl(Ntraj_text))//TRVfile
	else
        	old_filename = gridpath0//"/"//Ngrid_text//"/"//trim(adjustl(Ntraj_text))//TRVfile
        	call system("rm "//trim(adjustl(old_filename)))
	end if

	open(filechannel1,file=gridpath0//Ngrid_text//"/"//prefix_filename//TRVfile)
	open(filechannel2,file=trim(adjustl(old_filename)),position="append")
	do i = 1, Ntraj_max
		read(filechannel1,FMTtrv) (TRVenergies(j),j=1,3)

		TRVenergies_rounded(1:3) = ceiling(TRVenergies / bin_width)
		TRVenergies_rounded(4) = ceiling((TRVenergies(2)+TRVenergies(3)) / bin_width)
		do j = 1, 4
			if (TRVenergies_rounded(j) > TRV_Nbins) TRVenergies_rounded(j) = TRV_Nbins
			if (TRVenergies_rounded(j) == 0) TRVenergies_rounded(j) = 1
		end do

		write(filechannel2,FMT=*) (TRVenergies_rounded(j),j=1,4)
	end do
	close(filechannel1)
	close(filechannel2)


        !This is the gnuplot code to make the plots
        open(gnuplotchannel,file=gridpath0//gnuplotfile)
        write(gnuplotchannel,*) 'set term jpeg size 1200,1200'
        write(gnuplotchannel,*) 'set output "'//gridpath0//trim(adjustl(Ntraj_text))//JPGfilename//'"'
        write(gnuplotchannel,*) 'set multiplot layout 4,1'
        write(gnuplotchannel,*) 'set title "Translational Energy Change Distribution of '//trim(adjustl(Ntraj_text))//&
                                ' Trajectories"'
        write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
        write(gnuplotchannel,*) 'unset key'
	write(variable_length_text,FMT="(I5)") TRV_Nbins
	write(gnuplotchannel,*) 'Nbins = '//variable_length_text
	write(variable_length_text,FMT="(F5.3)") energy_max
	write(gnuplotchannel,*) 'energy_max = '//variable_length_text
	write(gnuplotchannel,*) 'bin_width = energy_max / Nbins'
	write(gnuplotchannel,*) 'set boxwidth bin_width'
        write(gnuplotchannel,*) 'set xlabel "Energy Change (eV)"'
        write(gnuplotchannel,*) 'set ylabel "Occurence"'
        write(gnuplotchannel,*) 'set xrange [0:energy_max]'
        write(gnuplotchannel,*) 'set autoscale y'
        write(gnuplotchannel,*) 'plot "'//trim(adjustl(old_filename))//&
                                '" u (bin_width*(($1)-0.5)):(1.0) smooth frequency with boxes'
        write(gnuplotchannel,*) 'set title "Rovibrational Energy Change Distribution of '//trim(adjustl(Ntraj_text))//&
                                ' Trajectories"'
        write(gnuplotchannel,*) 'plot "'//trim(adjustl(old_filename))//&
                                '" u (bin_width*(($4)-0.5)):(1.0) smooth frequency with boxes'
        write(gnuplotchannel,*) 'set title "Rotational Energy Change Distribution of '//trim(adjustl(Ntraj_text))//&
                                ' Trajectories"'
        write(gnuplotchannel,*) 'plot "'//trim(adjustl(old_filename))//&
                                '" u (bin_width*(($2)-0.5)):(1.0) smooth frequency with boxes'
        write(gnuplotchannel,*) 'set title "Vibrational Energy Change Distribution of '//trim(adjustl(Ntraj_text))//&
                                ' Trajectories"'
        write(gnuplotchannel,*) 'plot "'//trim(adjustl(old_filename))//&
                                '" u (bin_width*(($3)-0.5)):(1.0) smooth frequency with boxes'
        close(gnuplotchannel)

        !And then we just input it into gnuplot.exe
        call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)

end do

end subroutine getTRVimages




subroutine getSAfiles(prefix_filename)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
use PHYSICS
implicit none

!FILENAMING
character(*),intent(in) :: prefix_filename

!SCATTERING ANGLE
real(dp) :: scatteringAngle

!TRAJECTORY DATA
real(dp),dimension(3,3) :: coords_initial,velocities_initial
real(dp),dimension(3,3) :: coords_final,velocities_final

!INCREMENTAL INTEGERS
integer :: i, j, k

	open(filechannel1,file=gridpath0//prefix_filename//timeslicefile)
	open(filechannel2,file=gridpath0//prefix_filename//SAfile)
	do k = 1, Ntraj_max
		read(filechannel1,FMTtimeslice) &
                        ((coords_initial(i,j),i=1,3),j=1,Natoms),&
                        ((velocities_initial(i,j),i=1,3),j=1,Natoms),&
                        ((coords_final(i,j),i=1,3),j=1,Natoms),&
                        ((velocities_final(i,j),i=1,3),j=1,Natoms)

		scatteringAngle = acos(dot_product(velocities_initial(:,1),velocities_final(:,1)) / &
                                  sqrt(sum(velocities_final(:,1)**2) * sum(velocities_final(:,1)**2)))
		write(filechannel2,FMTsa) scatteringAngle
	end do
	close(filechannel2)
	close(filechannel1)

end subroutine getSAfiles


subroutine getTRVfiles(prefix_filename)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
use PHYSICS
implicit none

!FILENAMING
character(*),intent(in) :: prefix_filename

!Translational, Rotational, Vibrational Energies
real(dp),dimension(3) :: TRVenergies1, TRVenergies2, dTRVenergies
real(dp),dimension(3) :: velocity_out
real(dp) :: bond_distance, rot_speed

!TRAJECTORY DATA
real(dp),dimension(Nbonds,5) :: initial_bonding_data
real(dp),dimension(3,3) :: coords_initial,velocities_initial
real(dp),dimension(3,3) :: coords_final,velocities_final

!BONDING DATA INDEXES
integer :: atom1, atom2

!INCREMENTAL INTEGERS
integer :: i, j, k

	open(filechannel1,file=gridpath0//prefix_filename//timeslicefile)
	open(filechannel2,file=gridpath0//prefix_filename//initialfile)
	open(filechannel3,file=gridpath0//prefix_filename//TRVfile)
	do k = 1, Ntraj_max
		read(filechannel1,FMTtimeslice) &
                        ((coords_initial(i,j),i=1,3),j=1,Natoms),&
                        ((velocities_initial(i,j),i=1,3),j=1,Natoms),&
                        ((coords_final(i,j),i=1,3),j=1,Natoms),&
                        ((velocities_final(i,j),i=1,3),j=1,Natoms)
		read(filechannel2,FMTinitial) ((initial_bonding_data(i,j),j=1,5),i=1,Nbonds)

		TRVenergies1 = 0.0d0
		TRVenergies2 = 0.0d0

		do i = 1, Nbonds
			bond_distance = initial_bonding_data(i,1)
			rot_speed = initial_bonding_data(i,2)
			TRVenergies1(2) = TRVenergies1(2) + mass_hydrogen*rot_speed
		        TRVenergies1(3) = TRVenergies1(3) + PotentialConstant2*bond_distance**2 + &
                                          PotentialConstant1*bond_distance + PotentialConstant0

			atom1 = BONDING_DATA(i,1)
			atom2 = BONDING_DATA(i,2)
			call decompose_two_velocities(coords_final(:,atom1:atom2),velocities_final(:,atom1:atom2),&
                                                      velocity_out,TRVenergies2)
		        TRVenergies2(3) = TRVenergies2(3) + &
                                          HOPotential(coords_final(:,atom1),coords_final(:,atom2))
		end do
		do i = 1, Natoms
			if (all(BONDING_DATA /= i)) then
                        TRVenergies1(1) = TRVenergies1(1) + KineticEnergy(velocities_initial(:,i))
                        TRVenergies2(1) = TRVenergies2(1) + KineticEnergy(velocities_final(:,i))
			end if
		end do

		write(filechannel3,FMTtrv) (dTRVenergies(i),i=1,3)
	end do
	close(filechannel2)
	close(filechannel1)

end subroutine getTRVfiles





subroutine postProcess(prefix_filename)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
use PHYSICS
implicit none

!FILENAMING
character(*),intent(in) :: prefix_filename

!Translational, Rotational, Vibrational Energies
real(dp),dimension(3) :: TRVenergies1, TRVenergies2, dTRVenergies
real(dp),dimension(3) :: velocity_out
real(dp) :: scatteringAngle, speed_out
real(dp) :: bond_distance, rot_speed

!TRAJECTORY DATA
real(dp),dimension(Nbonds,5) :: initial_bonding_data
real(dp),dimension(3,3) :: coords_initial,velocities_initial
real(dp),dimension(3,3) :: coords_final,velocities_final

!BONDING DATA INDEXES
integer :: atom1, atom2

!INCREMENTAL INTEGERS
integer :: i, j, k

	open(filechannel1,file=gridpath0//prefix_filename//timeslicefile)
	open(filechannel2,file=gridpath0//prefix_filename//initialfile)
	open(filechannel3,file=gridpath0//prefix_filename//TRVfile)
	open(filechannel4,file=gridpath0//prefix_filename//SAfile)
	do k = 1, Ntraj_max
		read(filechannel1,FMTtimeslice) &
                        ((coords_initial(i,j),i=1,3),j=1,Natoms),&
                        ((velocities_initial(i,j),i=1,3),j=1,Natoms),&
                        ((coords_final(i,j),i=1,3),j=1,Natoms),&
                        ((velocities_final(i,j),i=1,3),j=1,Natoms)
		read(filechannel2,FMTinitial) ((initial_bonding_data(i,j),j=1,5),i=1,Nbonds)

		TRVenergies1 = 0.0d0
		TRVenergies2 = 0.0d0

		do i = 1, Nbonds
			bond_distance = initial_bonding_data(i,1)
			rot_speed = initial_bonding_data(i,2)
			TRVenergies1(2) = TRVenergies1(2) + mass_hydrogen*rot_speed
		        TRVenergies1(3) = TRVenergies1(3) + PotentialConstant2*bond_distance**2 + &
                                          PotentialConstant1*bond_distance + PotentialConstant0

			atom1 = BONDING_DATA(i,1)
			atom2 = BONDING_DATA(i,2)
			call decompose_two_velocities(coords_final(:,atom1:atom2),velocities_final(:,atom1:atom2),&
                                                      velocity_out,TRVenergies2)
		        TRVenergies2(3) = TRVenergies2(3) + &
                                          HOPotential(coords_final(:,atom1),coords_final(:,atom2))
		end do
		do i = 1, Natoms
			if (all(BONDING_DATA /= i)) then
                        TRVenergies1(1) = TRVenergies1(1) + KineticEnergy(velocities_initial(:,i))
                        TRVenergies2(1) = TRVenergies2(1) + KineticEnergy(velocities_final(:,i))
			end if
		end do

		dTRVenergies = abs(TRVenergies1 - TRVenergies2)
		write(filechannel3,FMTtrv) (dTRVenergies(i),i=1,3)

		speed_out = sqrt(sum(velocities_final(:,1)**2))
		scatteringAngle = acos(dot_product(velocities_initial(:,1),velocities_final(:,1)) / &
                                  (sqrt(sum(velocities_initial(:,1)**2)) * speed_out))
		write(filechannel4,FMTsa) scatteringAngle, speed_out
	end do
	close(filechannel4)
	close(filechannel3)
	close(filechannel2)
	close(filechannel1)

end subroutine postProcess


end module analyzeScatteringAngleswithMultipleGrids
