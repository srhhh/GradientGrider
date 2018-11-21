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
use PHYSICS
use ANALYSIS
use FUNCTIONS
implicit none

!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
character(*), intent(in) :: prefix_filename

!FORMAT OF JPG FILES TO BE MADE
character(*), intent(in) :: JPGfilename

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(5) :: variable_length_text1, variable_length_text2
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectory_text_length) :: Ntraj_text
character(150) :: old_filename

!SCATTERING ANGLES
real(dp),dimension(3,3) :: coords_initial,velocities_initial,coords_final,velocities_final
real(dp) :: ScatteringAngle_real, TranslationalEnergy_real
integer :: ScatteringAngle, TranslationalEnergy
integer :: AbsEnergyChange, RelEnergyChange, RotEnergyChange

!HEAT MAP VARIABLES
integer,allocatable :: angle_energy_bins(:,:)
integer :: occurence_max
real :: bin_width
integer :: angle_ratio = angleBins / scatteringangleBins
integer :: angle_slice = 6!denominator for slicing the scattering angle plot

!I/O HANDLING
integer :: iostate

!Incremental Integers
integer :: i, j, k


!The plots are named starting from Ntraj_max by increments of Ntraj_max (the number of trajectories)
write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
write(Ntraj_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ntraj

allocate(angle_energy_bins(angleBins,energyBins))
angle_energy_bins = 0

sizeAbsEnergyBin = (max_absenergychange-min_absenergychange) / energychangeBins
sizeRelEnergyBin = (max_relenergychange-min_relenergychange) / energychangeBins
sizeRotEnergyBin = (max_rotenergychange-min_rotenergychange) / energychangeBins

sizeDeltaEnergyBin = (max(max_absenergychange,max_relenergychange,max_rotenergychange) - &
                      min(min_absenergychange,min_relenergychange,min_rotenergychange)) / energychangeBins
sizeEnergyBin = max_TranslationalEnergy / (energyBins)
sizeAngleBin = pi / (angleBins)

if (comparison_upperlimit == comparison_lowerlimit) then
else

if (trim(adjustl(comparison_SATRVname)) == "ScatteringAngle") then
        sizeAngleBin = (comparison_upperlimit - comparison_lowerlimit) / angleBins
else if (trim(adjustl(comparison_SATRVname)) == "RelativeEnergyChange") then
        sizeRelEnergyBin = (comparison_upperlimit - comparison_lowerlimit) / energychangeBins
else if (trim(adjustl(comparison_SATRVname)) == "AbsoluteEnergyChange") then
        sizeAbsEnergyBin = (comparison_upperlimit - comparison_lowerlimit) / energychangeBins
else if (trim(adjustl(comparison_SATRVname)) == "RotationalEnergyChange") then
        sizeRotEnergyBin = (comparison_upperlimit - comparison_lowerlimit) / energychangeBins
else
end if

end if

open(filechannel1,file=gridpath0//prefix_filename//SATRVfile)
open(filechannel2,file=gridpath0//prefix_filename//binnedSATRVfile)
do i = 1, Ntraj
        read(filechannel1,FMT=FMTdata,iostat=iostate) ScatteringAngle_real, TranslationalEnergy_real, &
						      abs_energychange, rel_energychange, rot_energychange

	ScatteringAngle = ceiling(ScatteringAngle_real / sizeAngleBin)
	TranslationalEnergy = ceiling((TranslationalEnergy_real)/ sizeEnergyBin)
	AbsEnergyChange = ceiling((abs_energychange)/sizeAbsEnergyBin)
	RelEnergyChange = ceiling((rel_energychange)/sizeRelEnergyBin)
	RotEnergyChange = ceiling((rot_energychange)/sizeRotEnergyBin)

	if (ScatteringAngle > angleBins) ScatteringAngle = angleBins
	if (ScatteringAngle == 0) ScatteringAngle = 1

	if (TranslationalEnergy > energyBins) TranslationalEnergy = energyBins
	if (TranslationalEnergy == 0) TranslationalEnergy = 1

	if (AbsEnergyChange > energychangeBins) AbsEnergyChange = energychangeBins
	if (AbsEnergyChange == 0) AbsEnergyChange = 1

	if (RelEnergyChange > energychangeBins) RelEnergyChange = energychangeBins
	if (RelEnergyChange == 0) RelEnergyChange = 1

	if (RotEnergyChange > energychangeBins) RotEnergyChange = energychangeBins
	if (RotEnergyChange == 0) RotEnergyChange = 1

	angle_energy_bins(ScatteringAngle,TranslationalEnergy) = &
                   angle_energy_bins(ScatteringAngle,TranslationalEnergy) + 1

	write(filechannel2,FMT=*) ceiling(scatteringAngle*1.0/angle_ratio), TranslationalEnergy, &
                                  AbsEnergyChange, RelEnergyChange, RotEnergyChange
end do
close(filechannel1)
close(filechannel2)


!This is the gnuplot code to make the plots
!This is only plotted if there is no comparison active

if (comparison_flag) return

open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath0//prefix_filename//JPGfilename//'.png"'
write(gnuplotchannel,*) 'set multiplot layout 4,1'
write(gnuplotchannel,*) 'set title "Scattering Angle Distribution of '//trim(adjustl(Ntraj_text))//&
                        ' Trajectories"'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'pi = 3.14159265'
write(variable_length_text,FMT="(F5.4)") sizeAngleBin*angle_ratio
write(gnuplotchannel,*) 'box_width = '//variable_length_text
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'set xlabel "Scattering Angle"'
write(gnuplotchannel,*) 'set ylabel "Occurence"'
write(gnuplotchannel,*) 'set xtics pi/2'
write(gnuplotchannel,*) "set format x '%.1P π'"
write(gnuplotchannel,*) 'set xrange [0:pi]'
write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//binnedSATRVfile//&
                        '" u (box_width*($1-0.5)):(1.0) smooth frequency with boxes'

write(gnuplotchannel,*) 'scaling = 1000'
write(gnuplotchannel,*) 'min_E = scaling * ', min_absenergychange
write(gnuplotchannel,*) 'max_E = scaling * ', max_absenergychange
write(gnuplotchannel,*) 'box_width = (max_E-min_E) /', energychangeBins
write(gnuplotchannel,*) 'set xrange [min_E:max_E]'
write(gnuplotchannel,*) 'set xtics min_E, box_width * 10, max_E'
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'set xlabel "Energy (eV)"'
write(gnuplotchannel,*) "set format x '%.3f'"
write(gnuplotchannel,*) 'set ylabel "Absolute Translational Energy Change"'
write(gnuplotchannel,*) 'set title "Absolute Translational Energy Change Distribution of '//trim(adjustl(Ntraj_text))//&
                        ' Trajectories"'
write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//binnedSATRVfile//&
                        '" u (box_width*($3-0.5)+min_E):(1.0) smooth frequency w boxes'
write(gnuplotchannel,*) 'set title "Relative Translational Energy Change Distribution of '//trim(adjustl(Ntraj_text))//&
                        ' Trajectories"'
write(gnuplotchannel,*) 'set ylabel "Relative Translational Energy Change"'
write(gnuplotchannel,*) 'min_E = scaling * ', min_relenergychange
write(gnuplotchannel,*) 'max_E = scaling * ', max_relenergychange
write(gnuplotchannel,*) 'box_width = (max_E-min_E) /', energychangeBins
write(gnuplotchannel,*) 'set xrange [min_E:max_E]'
write(gnuplotchannel,*) 'set xtics min_E, box_width * 10, max_E'
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//binnedSATRVfile//&
                        '" u (box_width*($4-0.5)+min_E):(1.0) smooth frequency w boxes'
write(gnuplotchannel,*) 'set title "Rotational Energy Change Distribution of '//trim(adjustl(Ntraj_text))//&
                        ' Trajectories"'
write(gnuplotchannel,*) 'set ylabel "Rotational Energy Change"'
write(gnuplotchannel,*) 'min_E = scaling * ', min_rotenergychange
write(gnuplotchannel,*) 'max_E = scaling * ', max_rotenergychange
write(gnuplotchannel,*) 'box_width = (max_E-min_E) /', energychangeBins
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'set xrange [min_E:max_E]'
write(gnuplotchannel,*) 'set xtics min_E, box_width * 10, max_E'
write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//binnedSATRVfile//&
                        '" u (box_width*($5-0.5)+min_E):(1.0) smooth frequency w boxes'
close(gnuplotchannel)

!And then we just input it into gnuplot.exe
call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)

occurence_max = maxval(angle_energy_bins)
open(filechannel1,file=gridpath0//temporaryfile1)
do i = 1, angleBins/angle_slice+1
	do j = 1, energyBins
		write(filechannel1,FMT=*) (i-1)*sizeAngleBin, (j-1)*sizeEnergyBin, angle_energy_bins(i,j)
	end do
	write(filechannel1,FMT=*) ""
end do
close(filechannel1)

!This is the gnuplot code to make the plots
open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath0//prefix_filename//'HeatMap_'//&
                        JPGfilename//'.png"'
write(gnuplotchannel,*) 'set title "Scattering Angle Distribution" offset 0,-20'
!write(gnuplotchannel,*) 'set lmargin at screen 0.05'
!write(gnuplotchannel,*) 'set rmargin at screen 0.85'
!write(gnuplotchannel,*) 'set bmargin at screen 0.1'
!write(gnuplotchannel,*) 'set tmargin at screen 0.7705'
write(gnuplotchannel,*) 'set lmargin at screen 0.05'
write(gnuplotchannel,*) 'set rmargin at screen 0.85'
write(gnuplotchannel,*) 'set bmargin at screen 0.1'
write(gnuplotchannel,*) 'set pm3d map'
write(gnuplotchannel,*) 'set pm3d corners2color c1'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set multiplot'
write(gnuplotchannel,*) 'unset border'
write(gnuplotchannel,*) 'unset xtics'
write(gnuplotchannel,*) 'unset ytics'
write(gnuplotchannel,*) 'set angles radian'
write(gnuplotchannel,*) 'r = ', max_TranslationalEnergy
write(gnuplotchannel,*) 'rmax = r*1.1'
write(gnuplotchannel,*) 'pi = 3.14159'
!write(gnuplotchannel,*) 'set palette rgb 7, 5, 15'
write(variable_length_text,FMT="(I5)") occurence_max
write(gnuplotchannel,*) 'set cbrange [0:'//trim(adjustl(variable_length_text))//']'
write(gnuplotchannel,*) 'set cblabel "Frequency"'
!write(gnuplotchannel,*) 'set palette positive nops_allcF maxcolors 0 gamma 1.5 color model XYZ '
!write(gnuplotchannel,*) 'set palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV'
!write(gnuplotchannel,*) 'set palette rgbformulae 3, 2, 2'
write(gnuplotchannel,*) 'set palette defined ( 0 0 1 0, 0.3333 0 0 1, 0.6667 1 0 0,\'
write(gnuplotchannel,*) '     1 1 0.6471 0 )'
write(gnuplotchannel,*) 'scaling_factor = 2.4'
write(gnuplotchannel,*) 'set size ratio -1'
write(gnuplotchannel,*) 'set xrange[0:scaling_factor*rmax]'
write(gnuplotchannel,*) 'set yrange[0:scaling_factor*rmax]'
write(gnuplotchannel,*) 'set colorbox user origin 0.9,0.1 size 0.03,0.5'
write(gnuplotchannel,*) 'set origin 0.0, 0.0'
write(gnuplotchannel,*) 'splot "'//gridpath0//temporaryfile1//'" u (scaling_factor*$2*cos($1)):'//&
                        '(scaling_factor*$2*sin($1)):3'
write(gnuplotchannel,*) 'set style line 11 lc rgb "black" lw 2'
!write(gnuplotchannel,*) 'set grid polar ls 11'
!write(gnuplotchannel,*) 'set polar'
!write(gnuplotchannel,*) 'set rrange[0:rmax]'
!write(gnuplotchannel,*) 'set yrange[0:rmax]'
!write(gnuplotchannel,*) 'set size ratio -1'
!write(gnuplotchannel,*) 'set bmargin at screen 0.2295'
!write(gnuplotchannel,*) 'set tmargin at screen 0.9'
!!write(gnuplotchannel,*) 'unset raxis'
!!write(gnuplotchannel,*) 'set rlabel "Translational Energy Change (eV)"'
!!write(gnuplotchannel,*) 'set rtics format "" scale 0'
!write(gnuplotchannel,*) 'unset title'
!write(gnuplotchannel,*) 'set rtics r / 4'
!write(gnuplotchannel,*) 'unset parametric'
!write(gnuplotchannel,*) 'set for [i=0:330:30] label at first (rmax*scaling_factor)*cos(i*pi/180),'//&
!                        ' first (rmax*scaling_factor)*sin(i*pi/180)\'
!write(gnuplotchannel,*) 'center sprintf(''%d'', i)'
!write(gnuplotchannel,*) 'plot NaN w l'
!write(gnuplotchannel,*) 'unset multiplot'
write(gnuplotchannel,*) 'set samples 1000'
write(gnuplotchannel,*) 'phi_max = pi/', angle_slice
write(gnuplotchannel,*) 'dphi = phi_max/3'
write(gnuplotchannel,*) 'dR = rmax/4'
write(gnuplotchannel,*) 'set xr [0:rmax]'
write(gnuplotchannel,*) 'set yr [0:rmax]'
write(gnuplotchannel,*) 'set xtics out nomirror'
write(gnuplotchannel,*) 'unset ytics'
write(gnuplotchannel,*) 'set style line 42 lc rgb "black" dt 3'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set origin 0.08, 0.09'
write(gnuplotchannel,*) 'set lmargin at screen 0.050'
write(gnuplotchannel,*) 'set rmargin at screen 0.77'
write(gnuplotchannel,*) 'set bmargin at screen 0.13'
write(gnuplotchannel,*) 'unset title'
write(gnuplotchannel,*) "set format x '%.3f'"
write(gnuplotchannel,*) 'set xtics 0, rmax/4, rmax'
write(gnuplotchannel,*) 'set xlabel "Translational Energy (eV)"'
write(gnuplotchannel,*) 'set for [i=0:phi_max/dphi+1] label at first (rmax*1.05)*cos(i*dphi),'//&
                        ' first (rmax*1.05)*sin(i*dphi) center sprintf(''%d^o'',i*dphi*1.02*180/pi)'
write(gnuplotchannel,*) 'plot for [i=1:ceil(rmax/dR)] "+" u (i*dR*cos($1*phi_max/rmax)):(i*dR*sin($1*phi_max/rmax))'//&
                        ' w l ls 42'! scale 0.8, 0.8'
write(gnuplotchannel,*) 'set origin 0.08, 0.09'
write(gnuplotchannel,*) 'plot for [i=0:phi_max/dphi+1] "+" u ($1*cos(i*dphi)):($1*sin(i*dphi))'//&
                        ' w l ls 42'! scale 0.8, 0.8'
!write(gnuplotchannel,*) 'plot \'
!write(gnuplotchannel,*) '    for [i=0:phi_max/dphi] (x>=0&&x<=rmax*cos(i*dphi))?tan(i*dphi)*x:1/0 w l ls 42'
!write(gnuplotchannel,*) 'set polar'
!write(gnuplotchannel,*) 'set trange [0:phi_max]'
!write(gnuplotchannel,*) 'unset raxis'
!write(gnuplotchannel,*) 'unset rtics'
!write(gnuplotchannel,*) 'plot \'
!write(gnuplotchannel,*) '    for [i=1:ceil(rmax/dR)] i*dR<=rmax?i*dR:1/0 w l ls 42'
!write(gnuplotchannel,*) 'unset raxis'
!write(gnuplotchannel,*) 'unset rtics'
close(gnuplotchannel)

!And then we just input it into gnuplot.exe
call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)


end subroutine getScatteringAngles1









subroutine getConvergenceImage(lowerlimit,upperlimit,SATRVcolumn,SATRVname)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

real,intent(in) :: lowerlimit, upperlimit
character(*),intent(in) :: SATRVname
integer,intent(in) :: SATRVcolumn
integer,dimension(5) :: SATRVdata

!SCATTERING ANGLE BINNING
integer :: Nsamples, Nsamples_max
real :: bin_width, binMean, binSD,binRMSD
real,allocatable :: binAverage(:), sampleKS(:), sampleRMSD(:)
integer,allocatable :: binCumulative(:,:), binTotal(:,:),sampleSize(:)
integer :: binTally
real :: binThreshold = 1.0

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(5) :: variable_length_text1, variable_length_text2
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectory_text_length) :: Ntraj_text

!INCREMENTAL INTEGER
integer :: i

bin_width = (upperlimit-lowerlimit) / scatteringangleBins
Nsamples_max = (Ngrid_max * Ntraj_max) / Ntesttraj
allocate(binCumulative(scatteringangleBins,Nsamples_max),&
         binTotal(scatteringangleBins,Nsamples_max),&
         binAverage(scatteringangleBins),&
         sampleSize(Nsamples_max),sampleKS(Nsamples_max),sampleRMSD(Nsamples_max))

!Now we want a reference distribution
!This will be the distribution of ALL trajectories
binTotal = 0
binCumulative = 0
sampleSize = 0
write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
write(Ntraj_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid_max * Ntraj_max
write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid_max
open(filechannel1,file=gridpath0//Ngrid_text//"/Initial"//Ntraj_text//binnedSATRVfile,action="read")
do Nsamples = 1, Nsamples_max
        do i = 1, Ntesttraj
                read(filechannel1,FMT=*) SATRVdata
                binTotal(SATRVdata(SATRVcolumn),Nsamples) = &
                         binTotal(SATRVdata(SATRVcolumn),Nsamples) + 1
        end do
        sampleSize(Nsamples) = Nsamples
end do
close(filechannel1)

do i = 1,scatteringangleBins
        binAverage(i) = sum(binTotal(i,:)) * 1.0 / Nsamples_max
end do

!Now we want error bars
!For this, we will get a running average distribution by adding Ntesttraj trajectories at a time
!When the running average converges then we stop and see the variance
binTally = 0
sampleKS = 0
sampleRMSD = 0
open(filechannel1,file=gridpath0//"RMSD"//SATRVname//cumulativefile//".dat",action="write")
do Nsamples = 1, Nsamples_max
        binRMSD = 0.0
        do i = 1, scatteringangleBins
                binRMSD = binRMSD + (sum(binTotal(i,1:Nsamples)) * 1.0 / Nsamples - binAverage(i))**2
                sampleRMSD(Nsamples) = sampleRMSD(Nsamples) + (binTotal(i,Nsamples) - binAverage(i))**2
                binCumulative(i,Nsamples) = sum(binTotal(1:i,Nsamples))
                sampleKS(Nsamples) = max(sampleKS(Nsamples),abs(binCumulative(i,Nsamples)- &
                                         sum(binAverage(1:i))))
        end do
        sampleRMSD(Nsamples) = sqrt(sampleRMSD(Nsamples) * 1.0 / (scatteringangleBins - 1.0))

        binRMSD = sqrt(binRMSD/scatteringangleBins)
        write(filechannel1,FMT=*) Nsamples*Ntesttraj, binRMSD, binThreshold,&
                                  sampleKS(Nsamples), sampleRMSD(Nsamples)

        if (binRMSD < binThreshold) then
                binTally = binTally + 1
        else
                binTally = 0
        end if

!        if (binTally == 5) exit
!        if (Nsamples == Nsamples_max) then
!                print *, "    No convergence for true scattering angle"
!                print *, ""
!                exit
!        end if
end do
close(filechannel1)

!Now we must do the laborious job of binning these
!Each bin has its own average and standard deviation based on how many samples we took
open(filechannel1,file=gridpath0//"Adjusted"//SATRVname//cumulativefile//".dat")
do i = 1, scatteringangleBins
        binMean = binAverage(i)
        binSD = sqrt(sum((binTotal(i,1:Nsamples_max)*1.0 - binMean)**2)/(Nsamples_max - 1))

        write(filechannel1,*) (i-0.5)*bin_width, binMean, binSD!/sqrt(real(Nsamples))
end do
close(filechannel1)

deallocate(binAverage,binTotal,sampleSize)

!We have everything we need to draw the distribution with error bars
open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) "set terminal pngcairo size 1200,1200"
write(variable_length_text,FMT="(I5)") Ntesttraj
write(gnuplotchannel,*) 'set output "'//gridpath0//'Convergence'//&
                        trim(adjustl(variable_length_text))//SATRVname//'.png"'
write(gnuplotchannel,*) 'set multiplot layout 2,2 columnsfirst'
write(variable_length_text,FMT="(I5)") Ntesttraj
write(gnuplotchannel,*) 'set title "Convergence of the '//trim(adjustl(variable_length_text))//&
                        ' '//SATRVname//' Distribution with '//trim(adjustl(Ngrid_text))//' Grids"'
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
write(gnuplotchannel,*) 'plot "'//gridpath0//'RMSD'//SATRVname//cumulativefile//'.dat" u 1:2 w lines, \'
write(gnuplotchannel,*) '     "'//gridpath0//'RMSD'//SATRVname//cumulativefile//'.dat" u 1:3 w lines lc -1'
write(gnuplotchannel,*) 'unset label 1'
write(variable_length_text,FMT="(I5)") Ntesttraj
write(gnuplotchannel,*) 'set title "Final '//SATRVname//' Distribution for '//&
                        trim(adjustl(variable_length_text))//' Trajectories"'
if (SATRVname == "ScatteringAngle") then
write(gnuplotchannel,*) 'scaling = 1'
write(gnuplotchannel,*) 'pi = 3.14159265'
write(gnuplotchannel,*) 'lowerlimit = ', lowerlimit
write(gnuplotchannel,*) 'upperlimit = ', upperlimit
write(gnuplotchannel,*) 'set xtics lowerlimit, 10*(upperlimit-lowerlimit)/',SA_Nbins,', upperlimit'
write(gnuplotchannel,*) "set format x '%.3P π'"
write(gnuplotchannel,*) 'set xrange [lowerlimit:upperlimit]'
write(gnuplotchannel,*) 'set xlabel "Scattering Angle (rad)"'
else
write(gnuplotchannel,*) 'scaling = 1000'
write(gnuplotchannel,*) 'E_min = scaling * ', lowerlimit
write(gnuplotchannel,*) 'E_max = scaling * ', upperlimit
write(gnuplotchannel,*) 'set xrange [E_min:E_max]'
write(gnuplotchannel,*) 'set xtics E_min, 10*(E_max-E_min)/',SA_Nbins,', E_max'
write(gnuplotchannel,*) "set format x '%.3f'"
write(gnuplotchannel,*) 'set xlabel "Energy (meV)"'
end if
write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set ytics autofreq'
write(gnuplotchannel,*) 'set boxwidth'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set ylabel "Frequency"'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'plot "'//gridpath0//'Adjusted'//SATRVname//cumulativefile//'.dat" u (scaling*($1)):2 w boxes, \'
write(gnuplotchannel,*) '     "'//gridpath0//'Adjusted'//SATRVname//cumulativefile//'.dat" u (scaling*($1)):2:3 w yerrorbars'
write(gnuplotchannel,*) 'set title "Distribution of the Kolmogorov-Smirnov Difference Among '//&
                        trim(adjustl(variable_length_text))//' '//SATRVname//' Samplings Across '//&
                        trim(adjustl(Ngrid_text))//' Grids"'
write(gnuplotchannel,*) 'set xlabel "Kolmogorov-Smirnov Difference"'
write(gnuplotchannel,*) 'set ylabel "Frequency"'
write(gnuplotchannel,*) 'minKS = ', max(0.0,minval(sampleKS)-2*(maxval(sampleKS)-minval(sampleKS))/Nsamples_max)
write(gnuplotchannel,*) 'maxKS = ', maxval(sampleKS) + 2*(maxval(sampleKS)-minval(sampleKS))/Nsamples_max
write(gnuplotchannel,*) 'set xrange [minKS:maxKS]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'box_width = (maxKS-minKS) / ',(2*Nsamples_max)
write(gnuplotchannel,*) 'set xtics minKS, (maxKS-minKS)/4, maxKS'
write(gnuplotchannel,*) "set format x '%.2f'"
write(gnuplotchannel,*) 'set ytics 1'
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'bin_number(x) = floor(x/box_width)'
write(gnuplotchannel,*) 'rounded(x) = box_width * (bin_number(x) + 0.5)'
write(gnuplotchannel,*) 'plot "'//gridpath0//'RMSD'//SATRVname//cumulativefile//'.dat"'//&
                        'u (rounded($4)):(1.0) smooth frequency with boxes'
write(gnuplotchannel,*) 'set title "Distribution of the Root Mean Square Difference Among '//&
                        trim(adjustl(variable_length_text))//' '//SATRVname//' Samplings Across '//&
                        trim(adjustl(Ngrid_text))//' Grids"'
write(gnuplotchannel,*) 'set xlabel "Root Mean Square Difference"'
write(gnuplotchannel,*) 'set ylabel "Frequency"'
write(gnuplotchannel,*) 'minRMSD = ', max(0.0,minval(sampleRMSD)-2*(maxval(sampleRMSD)-minval(sampleRMSD))/Nsamples_max)
write(gnuplotchannel,*) 'maxRMSD = ', maxval(sampleRMSD) + 2*(maxval(sampleRMSD)-minval(sampleRMSD))/Nsamples_max
write(gnuplotchannel,*) 'set xrange [minRMSD:maxRMSD]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'box_width = (maxRMSD-minRMSD) / ',(2*Nsamples_max)
write(gnuplotchannel,*) 'set xtics minRMSD, (maxRMSD-minRMSD)/4, maxRMSD'
write(gnuplotchannel,*) "set format x '%.2f'"
write(gnuplotchannel,*) 'set ytics 1'
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'bin_number(x) = floor(x/box_width)'
write(gnuplotchannel,*) 'rounded(x) = box_width * (bin_number(x) + 0.5)'
write(gnuplotchannel,*) 'plot "'//gridpath0//'RMSD'//SATRVname//cumulativefile//'.dat"'//&
                        'u (rounded($5)):(1.0) smooth frequency with boxes'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)

end subroutine getConvergenceImage




!
!
!subroutine getInitialimages(prefix_filename,JPGfilename)
!use PARAMETERS
!use ANALYSIS
!use FUNCTIONS
!use PHYSICS
!implicit none
!
!!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
!character(*),intent(in) :: prefix_filename
!character(*),intent(in) :: JPGfilename
!
!!FORMATTING OF JPG FILES
!character(5) :: variable_length_text
!character(5) :: variable_length_text1, variable_length_text2
!character(Ngrid_text_length) :: Ngrid_text
!character(Ngrid_text_length+1) :: folder_text
!character(6) :: Ntraj_text
!character(150) :: old_filename
!
!integer :: iostate
!real :: r0, rot
!real :: max_r0, max_rot
!real :: min_r0, min_rot
!real :: average_r0, average_rot
!real :: average_Evib, average_Erot
!integer :: total_bonds
!integer :: i
!
!average_r0 = 0.0
!average_rot = 0.0
!average_Evib = 0.0
!average_Erot = 0.0
!max_r0 = 0.0
!max_rot = 0.0
!min_r0 = 1.0e9
!min_rot = 1.0e9
!total_bonds = 0
!open(filechannel1,file=gridpath0//prefix_filename//initialfile)
!do
!	read(filechannel1,iostat=iostate,FMT=FMTinitial) INITIAL_BOND_DATA
!	if (iostate /= 0) exit
!	do i = 1, Nbonds
!		r0 = INITIAL_BOND_DATA(1,i)
!		rot = INITIAL_BOND_DATA(2,i)
!		max_r0 = max(max_r0,r0)
!		max_rot = max(max_rot,rot)
!		min_r0 = min(min_r0,r0)
!		min_rot = min(min_rot,rot)
!		average_r0 = average_r0 + r0
!		average_rot = average_rot + rot
!		average_Evib = average_Evib + (r0 - HOr0_hydrogen)**2
!		average_Erot = average_Erot + (rot**2)
!		total_bonds = total_bonds + 1
!	end do
!end do
!close(filechannel1)
!
!average_r0 = average_r0 / total_bonds
!average_rot = average_rot / total_bonds
!average_Evib = 0.5 * HOke_hydrogen * average_Evib / total_bonds
!average_Erot = average_Erot / (total_bonds * mass_hydrogen)
!
!open(gnuplotchannel,file=gridpath0//gnuplotfile)
!write(gnuplotchannel,*) 'set term pngcairo enhanced size 3600,1200'
!write(gnuplotchannel,*) 'set encoding utf8'
!write(gnuplotchannel,*) 'set output "'//gridpath0//JPGfilename//'.png"'
!write(gnuplotchannel,*) 'unset key'
!write(gnuplotchannel,*) 'pi = 3.14159265'
!write(gnuplotchannel,*) 'set style histogram clustered gap 1'
!write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
!write(gnuplotchannel,*) 'set label 1 "Average r0: ', average_r0, ' A" at screen 0.6,0.9'
!write(gnuplotchannel,*) 'set label 2 "Temperature: ', (RU_energy/kb)*average_Evib, ' K" at screen 0.6,0.85'
!write(gnuplotchannel,*) 'set label 3 "Average Rotational Speed: ', average_rot, ' A/fs" at screen 0.8,0.9'
!write(gnuplotchannel,*) 'set label 4 "Temperature: ', (RU_energy/kb)*average_Erot, ' K" at screen 0.8,0.85'
!write(Ntraj_text,FMT="(I6)") Ntraj
!write(gnuplotchannel,*) 'set multiplot layout ', Nbonds,',4 title '//&
!                        '"Initial Bond Distribution of '//trim(adjustl(Ntraj_text))//' Trajectories"'
!do i = 1, Nbonds
!write(gnuplotchannel,*) 'set ylabel "Initial H2 Theta Occurence"'
!write(gnuplotchannel,*) 'box_width = 2 * pi /', SA_Nbins
!write(gnuplotchannel,*) 'set boxwidth box_width'
!write(gnuplotchannel,*) 'bin_number(x) = floor(x/box_width)'
!write(gnuplotchannel,*) 'rounded(x) = box_width * (bin_number(x) + 0.5)'
!write(gnuplotchannel,*) 'set xrange [-pi:pi]'
!write(gnuplotchannel,*) 'set xlabel "Angle (rad)"'
!write(gnuplotchannel,*) 'set yrange [0:]'
!write(gnuplotchannel,*) 'set xtics pi/2'
!write(gnuplotchannel,*) "set format x '%.1P π'"
!write(Ntraj_text,FMT="(I6)") i*6 - 2
!write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//initialfile//&
!                        '" u (rounded($'//trim(adjustl(Ntraj_text))//&
!			')):(1.0) smooth frequency with boxes'
!write(gnuplotchannel,*) 'set ylabel "Initial H2 Phi Occurence"'
!write(gnuplotchannel,*) 'box_width = pi /', SA_Nbins
!write(gnuplotchannel,*) 'set boxwidth box_width'
!write(gnuplotchannel,*) 'set xrange [0:pi]'
!write(gnuplotchannel,*) 'set yrange [0:]'
!write(Ntraj_text,FMT="(I6)") i*6 - 1
!write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//initialfile//&
!			'" u (rounded($'//trim(adjustl(Ntraj_text))//&
!			')):(1.0) smooth frequency with boxes'
!write(gnuplotchannel,*) 'min_r0 = ', min_r0
!write(gnuplotchannel,*) 'max_r0 = ', max_r0
!write(gnuplotchannel,*) 'box_width = (max_r0-min_r0) /', SA_Nbins
!write(gnuplotchannel,*) 'set boxwidth box_width'
!write(gnuplotchannel,*) 'bin_number(x) = floor((x-min_r0)/box_width)'
!write(gnuplotchannel,*) 'rounded(x) = min_r0+box_width * (bin_number(x) + 0.5)'
!write(gnuplotchannel,*) 'set xlabel "Bond Distance (A)"'
!write(gnuplotchannel,*) 'set xtics min_r0, (max_r0-min_r0)/5, max_r0'
!write(gnuplotchannel,*) "set format x '%.4f'"
!write(gnuplotchannel,*) 'set ylabel "Initial H2 Bond Length Occurence"'
!write(gnuplotchannel,*) 'set xrange [min_r0:max_r0]'
!write(gnuplotchannel,*) 'set yrange [0:]'
!write(gnuplotchannel,*) 'set bmargin 3'
!write(Ntraj_text,FMT="(I6)") i*6 - 5
!write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//initialfile//&
!			'" u (rounded($'//trim(adjustl(Ntraj_text))//&
!			')):(1.0) smooth frequency with boxes'
!write(gnuplotchannel,*) 'min_rot = ', min_rot
!write(gnuplotchannel,*) 'max_rot = ', max_rot
!write(gnuplotchannel,*) 'box_width = (max_rot-min_rot) /', SA_Nbins
!write(gnuplotchannel,*) 'set boxwidth box_width'
!write(gnuplotchannel,*) 'bin_number(x) = floor((x-min_rot)/box_width)'
!write(gnuplotchannel,*) 'rounded(x) = min_rot+box_width * (bin_number(x) + 0.5)'
!write(gnuplotchannel,*) 'set xlabel "Rotational Speed (A/fs)"'
!write(gnuplotchannel,*) 'set xtics min_rot, (max_rot-min_rot)/5, max_rot'
!write(gnuplotchannel,*) 'set ylabel "Initial H2 Rotational Speed Occurence"'
!write(gnuplotchannel,*) 'set xrange [min_rot:max_rot]'
!write(gnuplotchannel,*) 'set yrange [0:]'
!write(Ntraj_text,FMT="(I6)") i*6 - 4
!write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//initialfile//&
!			'" u (rounded($'//trim(adjustl(Ntraj_text))//&
!			')):(1.0) smooth frequency with boxes'
!end do
!
!close(gnuplotchannel)
!
!!And then we just input it into gnuplot.exe
!call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)
!
!
!end subroutine getInitialimages
!




subroutine getScatteringAngles2(prefix_filename,JPGfilename)
use PARAMETERS
use PHYSICS
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
character(trajectory_text_length) :: Ntraj_text
character(6) :: boxwidth_text
logical :: grid_is_done

integer :: iostate
real :: speed_out, ScatteringAngle

inquire(file=gridpath0//'AdjustedScatteringAngle'//cumulativefile//'.dat',exist=grid_is_done)

!This is the gnuplot code to make the plots
open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo enhanced size 1200,1200'
write(gnuplotchannel,*) 'set encoding utf8'
write(gnuplotchannel,*) 'set output "'//gridpath0//JPGfilename//'.png"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'pi = 3.14159265'
write(gnuplotchannel,*) 'box_width = pi / ', SA_Nbins
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'bin_number(x) = floor(x/box_width)'
write(gnuplotchannel,*) 'rounded(x) = box_width * (bin_number(x) + 0.5)'
write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
write(Ntraj_text,FMT="(I0."//variable_length_text//")") Ntraj
write(gnuplotchannel,*) 'set multiplot layout 4,1 title '//&
                        '"Angle Distribution of '//trim(adjustl(Ntraj_text))//' Trajectories"'
write(gnuplotchannel,*) 'set style histogram clustered gap 1'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'set ylabel "Scattering Angle Occurence"'
write(gnuplotchannel,*) 'set xlabel "Angle (rad)"'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set xrange [0:pi]'
write(gnuplotchannel,*) 'set xtics pi/2'
write(gnuplotchannel,*) "set format x '%.1P π'"

if (grid_is_done) then
        write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//SATRVfile//&
                                '" u (rounded($1)):(1.0) smooth frequency w boxes, \'
        write(gnuplotchannel,*) '     "'//gridpath0//'AdjustedScatteringAngle'//&
				cumulativefile//'.dat" u 1:2 w boxes'//&
                                       ' fs transparent solid 0.5 noborder, \'
        write(gnuplotchannel,*) '     "'//gridpath0//'AdjustedScatteringAngle'//&
				cumulativefile//'.dat" u 1:2:3 w yerrorbars'
else
        write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//SATRVfile//&
                                '" u (rounded($1)):(1.0) smooth frequency w boxes'
end if
        write(gnuplotchannel,*) 'scaling = 1000'
        write(gnuplotchannel,*) 'min_E = scaling * ', min_absenergychange
        write(gnuplotchannel,*) 'max_E = scaling * ', max_absenergychange
        write(gnuplotchannel,*) 'box_width = (max_E-min_E) /', SA_Nbins
        write(gnuplotchannel,*) 'set xrange [min_E:max_E]'
        write(gnuplotchannel,*) 'set xtics min_E, box_width * 10, max_E'
        write(gnuplotchannel,*) 'set boxwidth box_width'
	write(gnuplotchannel,*) 'bin_number(x) = floor(scaling*x/box_width)'
        write(gnuplotchannel,*) 'rounded(x) = min_E + box_width * (bin_number(x) + 0.5)'
        write(gnuplotchannel,*) 'set xlabel "Energy (meV)"'
        write(gnuplotchannel,*) "set format x '%.3f'"
        write(gnuplotchannel,*) 'set ylabel "Absolute Translational Energy Change"'
if (grid_is_done) then
        write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//SATRVfile//&
                                '" u (rounded($3)):(1.0) smooth frequency w boxes, \'
        write(gnuplotchannel,*) '     "'//gridpath0//'AdjustedAbsoluteEnergyChange'//&
				cumulativefile//'.dat" u (scaling*($1)):2 w boxes'//&
                                       ' fs transparent solid 0.5 noborder, \'
        write(gnuplotchannel,*) '     "'//gridpath0//'AdjustedAbsoluteEnergyChange'//&
				cumulativefile//'.dat" u (scaling*($1)):2:3 w yerrorbars'
else
        write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//SATRVfile//&
                                '" u (rounded($3)):(1.0) smooth frequency w boxes'
end if
        write(gnuplotchannel,*) 'set ylabel "Relative Translational Energy Change"'
        write(gnuplotchannel,*) 'min_E = scaling * ', min_relenergychange
        write(gnuplotchannel,*) 'max_E = scaling * ', max_relenergychange
        write(gnuplotchannel,*) 'box_width = (max_E-min_E) /', SA_Nbins
        write(gnuplotchannel,*) 'set xrange [min_E:max_E]'
        write(gnuplotchannel,*) 'set xtics min_E, box_width * 10, max_E'
        write(gnuplotchannel,*) 'set boxwidth box_width'
	write(gnuplotchannel,*) 'bin_number(x) = floor(scaling*x/box_width)'
        write(gnuplotchannel,*) 'rounded(x) = min_E + box_width * (bin_number(x) + 0.5)'
if (grid_is_done) then
        write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//SATRVfile//&
                                '" u (rounded($4)):(1.0) smooth frequency w boxes, \'
        write(gnuplotchannel,*) '     "'//gridpath0//'AdjustedRelativeEnergyChange'//&
				cumulativefile//'.dat" u (scaling*($1)):2 w boxes'//&
                                       ' fs transparent solid 0.5 noborder, \'
        write(gnuplotchannel,*) '     "'//gridpath0//'AdjustedRelativeEnergyChange'//&
				cumulativefile//'.dat" u (scaling*($1)):2:3 w yerrorbars'
else
        write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//SATRVfile//&
                                '" u (rounded($4)):(1.0) smooth frequency w boxes'
end if
        write(gnuplotchannel,*) 'set ylabel "Rotational Translational Energy Change"'
        write(gnuplotchannel,*) 'min_E = scaling * ', min_rotenergychange
        write(gnuplotchannel,*) 'max_E = scaling * ', max_rotenergychange
        write(gnuplotchannel,*) 'box_width = (max_E-min_E) /', SA_Nbins
        write(gnuplotchannel,*) 'set xrange [min_E:max_E]'
        write(gnuplotchannel,*) 'set xtics min_E, box_width * 10, max_E'
        write(gnuplotchannel,*) "set format x '%.3f'"
        write(gnuplotchannel,*) 'set boxwidth box_width'
	write(gnuplotchannel,*) 'bin_number(x) = floor(scaling*x/box_width)'
        write(gnuplotchannel,*) 'rounded(x) = min_E + box_width * (bin_number(x) + 0.5)'
if (grid_is_done) then
        write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//SATRVfile//&
                                '" u (rounded($5)):(1.0) smooth frequency w boxes, \'
        write(gnuplotchannel,*) '     "'//gridpath0//'AdjustedRotationalEnergyChange'//&
				cumulativefile//'.dat" u (scaling*($1)):2 w boxes'//&
                                       ' fs transparent solid 0.5 noborder, \'
        write(gnuplotchannel,*) '     "'//gridpath0//'AdjustedRotationalEnergyChange'//&
				cumulativefile//'.dat" u (scaling*($1)):2:3 w yerrorbars'
else
        write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//SATRVfile//&
                                '" u (rounded($5)):(1.0) smooth frequency w boxes'
end if
close(gnuplotchannel)

!And then we just input it into gnuplot.exe
call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)

end subroutine getScatteringAngles2




subroutine getInitialimages(prefix_filename,JPGfilename)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
use PHYSICS
implicit none

!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
character(*),intent(in) :: prefix_filename
character(*),intent(in) :: JPGfilename

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(5) :: variable_length_text1, variable_length_text2
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(6) :: Ntraj_text
character(150) :: old_filename

integer :: iostate
real :: r0, rot
real :: max_r0, max_rot
real :: min_r0, min_rot
real :: average_r0, average_rot
real :: average_Evib, average_Erot
integer :: total_bonds
integer :: i

average_r0 = 0.0
average_rot = 0.0
average_Evib = 0.0
average_Erot = 0.0
max_r0 = 0.0
max_rot = 0.0
min_r0 = 1.0e9
min_rot = 1.0e9
total_bonds = 0
open(filechannel1,file=gridpath0//prefix_filename//initialfile)
do
	read(filechannel1,iostat=iostate,FMT=FMTinitial) INITIAL_BOND_DATA
	if (iostate /= 0) exit
	do i = 1, Nbonds
		r0 = INITIAL_BOND_DATA(1,i)
		rot = INITIAL_BOND_DATA(2,i)
		max_r0 = max(max_r0,r0)
		max_rot = max(max_rot,rot)
		min_r0 = min(min_r0,r0)
		min_rot = min(min_rot,rot)
		average_r0 = average_r0 + r0
		average_rot = average_rot + rot
		average_Evib = average_Evib + (r0 - HOr0_hydrogen)**2
		average_Erot = average_Erot + (rot**2)
		total_bonds = total_bonds + 1
	end do
end do
close(filechannel1)

average_r0 = average_r0 / total_bonds
average_rot = average_rot / total_bonds
average_Evib = 0.5 * HOke_hydrogen * average_Evib / total_bonds
average_Erot = average_Erot / (total_bonds * mass_hydrogen)

open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo enhanced size 3600,1200'
write(gnuplotchannel,*) 'set encoding utf8'
write(gnuplotchannel,*) 'set output "'//gridpath0//JPGfilename//'.png"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'pi = 3.14159265'
write(gnuplotchannel,*) 'set style histogram clustered gap 1'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'set label 1 "Average r0: ', average_r0, ' A" at screen 0.6,0.9'
write(gnuplotchannel,*) 'set label 2 "Temperature: ', (RU_energy/kb)*average_Evib, ' K" at screen 0.6,0.85'
write(gnuplotchannel,*) 'set label 3 "Average Rotational Speed: ', average_rot, ' A/fs" at screen 0.8,0.9'
write(gnuplotchannel,*) 'set label 4 "Temperature: ', (RU_energy/kb)*average_Erot, ' K" at screen 0.8,0.85'
write(Ntraj_text,FMT="(I6)") Ntraj
write(gnuplotchannel,*) 'set multiplot layout ', Nbonds,',4 title '//&
                        '"Initial Bond Distribution of '//trim(adjustl(Ntraj_text))//' Trajectories"'
do i = 1, Nbonds
write(gnuplotchannel,*) 'set ylabel "Initial H2 Theta Occurence"'
write(gnuplotchannel,*) 'box_width = 2 * pi /', SA_Nbins
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'bin_number(x) = floor(x/box_width)'
write(gnuplotchannel,*) 'rounded(x) = box_width * (bin_number(x) + 0.5)'
write(gnuplotchannel,*) 'set xrange [-pi:pi]'
write(gnuplotchannel,*) 'set xlabel "Angle (rad)"'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set xtics pi/2'
write(gnuplotchannel,*) "set format x '%.1P π'"
write(Ntraj_text,FMT="(I6)") i*6 - 2
write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//initialfile//&
                        '" u (rounded($'//trim(adjustl(Ntraj_text))//&
			')):(1.0) smooth frequency with boxes'
write(gnuplotchannel,*) 'set ylabel "Initial H2 Phi Occurence"'
write(gnuplotchannel,*) 'box_width = pi /', SA_Nbins
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'set xrange [0:pi]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(Ntraj_text,FMT="(I6)") i*6 - 1
write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//initialfile//&
			'" u (rounded($'//trim(adjustl(Ntraj_text))//&
			')):(1.0) smooth frequency with boxes'
write(gnuplotchannel,*) 'min_r0 = ', min_r0
write(gnuplotchannel,*) 'max_r0 = ', max_r0
write(gnuplotchannel,*) 'box_width = (max_r0-min_r0) /', SA_Nbins
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'bin_number(x) = floor((x-min_r0)/box_width)'
write(gnuplotchannel,*) 'rounded(x) = min_r0+box_width * (bin_number(x) + 0.5)'
write(gnuplotchannel,*) 'set xlabel "Bond Distance (A)"'
write(gnuplotchannel,*) 'set xtics min_r0, (max_r0-min_r0)/5, max_r0'
write(gnuplotchannel,*) "set format x '%.4f'"
write(gnuplotchannel,*) 'set ylabel "Initial H2 Bond Length Occurence"'
write(gnuplotchannel,*) 'set xrange [min_r0:max_r0]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set bmargin 3'
write(Ntraj_text,FMT="(I6)") i*6 - 5
write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//initialfile//&
			'" u (rounded($'//trim(adjustl(Ntraj_text))//&
			')):(1.0) smooth frequency with boxes'
write(gnuplotchannel,*) 'min_rot = ', min_rot
write(gnuplotchannel,*) 'max_rot = ', max_rot
write(gnuplotchannel,*) 'box_width = (max_rot-min_rot) /', SA_Nbins
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'bin_number(x) = floor((x-min_rot)/box_width)'
write(gnuplotchannel,*) 'rounded(x) = min_rot+box_width * (bin_number(x) + 0.5)'
write(gnuplotchannel,*) 'set xlabel "Rotational Speed (A/fs)"'
write(gnuplotchannel,*) 'set xtics min_rot, (max_rot-min_rot)/5, max_rot'
write(gnuplotchannel,*) 'set ylabel "Initial H2 Rotational Speed Occurence"'
write(gnuplotchannel,*) 'set xrange [min_rot:max_rot]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(Ntraj_text,FMT="(I6)") i*6 - 4
write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//initialfile//&
			'" u (rounded($'//trim(adjustl(Ntraj_text))//&
			')):(1.0) smooth frequency with boxes'
end do

close(gnuplotchannel)

!And then we just input it into gnuplot.exe
call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)


end subroutine getInitialimages



subroutine getComparedScatteringAngles(lowerlimit,upperlimit,imagename,SATRVcolumn,SATRVname)
use PARAMETERS
use PHYSICS
use ANALYSIS
implicit none

!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
integer,intent(in) :: SATRVcolumn
character(*), intent(in) :: SATRVname
integer,dimension(5) :: SATRVdata

!FORMAT OF JPG FILES TO BE MADE
character(*), intent(in) :: imagename

!Upper and Lower Limits for the plot
real,intent(in) :: upperlimit, lowerlimit

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectory_text_length) :: Ntraj_text
character(6) :: boxwidth_text
logical :: grid_is_done

integer :: iostate,i,j
real :: speed_out, ScatteringAngle

real :: referenceBins(scatteringangleBins)
real :: referenceMeans(scatteringangleBins)
real :: referenceSDs(scatteringangleBins)
integer :: binTotal(scatteringangleBins,comparison_number)
real :: comparisonRMSD, comparisonKS, comparisonCDF
real :: comparisonBin, comparisonMean, comparisonSD

open(filechannel1,file=gridpath0//"Adjusted"//SATRVname//cumulativefile//".dat")
do j = 1, scatteringangleBins
        read(filechannel1,FMT=*) referenceBins(j), referenceMeans(j), referenceSDs(j)
end do
close(filechannel1)

binTotal = 0
open(filechannel1,file=gridpath0//allprefixes(1:alllengths(1))//binnedSATRVfile)
do j = 1, Ntesttraj
        read(filechannel1,FMT=*) SATRVdata
        binTotal(SATRVdata(SATRVcolumn),1) = &
                 binTotal(SATRVdata(SATRVcolumn),1) + 1
end do
close(filechannel1)

do i = 1, comparison_number-1
        open(filechannel1,file=gridpath0//allprefixes(sum(alllengths(1:i))+1:sum(alllengths(1:i+1)))//&
                               binnedSATRVfile)
        do j = 1, Ntesttraj
                read(filechannel1,FMT=*) SATRVdata
                binTotal(SATRVdata(SATRVcolumn),i+1) = &
                         binTotal(SATRVdata(SATRVcolumn),i+1) + 1
        end do
        close(filechannel1)
end do

!This is the gnuplot code to make the plots
open(gnuplotchannel,file=gridpath0//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo enhanced size 1200,1200'
write(gnuplotchannel,*) 'set encoding utf8'
write(gnuplotchannel,*) 'set output "'//gridpath0//imagename//'.png"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
write(gnuplotchannel,*) 'unset xtics'
write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
write(Ntraj_text,FMT="(I0."//variable_length_text//")") Ntraj
write(gnuplotchannel,*) 'set multiplot layout ',comparison_number,&
                        ',1 columnsfirst margins 0.1,0.95,.1,.9 spacing 0.1,0 title '//&
                        '"'//SATRVname//' Distribution of '//trim(adjustl(Ntraj_text))//' Trajectories"'
write(gnuplotchannel,*) 'set style histogram clustered gap 1'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'set ylabel "Occurence"'
write(gnuplotchannel,*) 'set yrange [0:]'

if (SATRVname == "ScatteringAngle") then
        write(gnuplotchannel,*) 'scaling = 1'
        write(gnuplotchannel,*) 'pi = 3.14159265'
        write(gnuplotchannel,*) "lowerlimit = ", lowerlimit
        write(gnuplotchannel,*) "upperlimit = ", upperlimit
        write(gnuplotchannel,*) 'set xrange [lowerlimit:upperlimit]'
        write(gnuplotchannel,*) 'box_width = (upperlimit-lowerlimit) / ', SA_Nbins
        write(gnuplotchannel,*) 'set boxwidth box_width'
        write(gnuplotchannel,*) 'bin_number(x) = floor(x/box_width)'
        write(gnuplotchannel,*) 'rounded(x) = box_width * (bin_number(x) + 0.5)'
else
        write(gnuplotchannel,*) 'scaling = 1000'
        write(gnuplotchannel,*) 'lowerlimit = scaling * ', lowerlimit
        write(gnuplotchannel,*) 'upperlimit = scaling * ', upperlimit
        write(gnuplotchannel,*) 'set xrange [lowerlimit:upperlimit]'
        write(gnuplotchannel,*) 'box_width = (upperlimit - lowerlimit) / ', SA_Nbins
        write(gnuplotchannel,*) 'set boxwidth box_width'
        write(gnuplotchannel,*) 'bin_number(x) = floor(scaling * x/box_width)'
        write(gnuplotchannel,*) 'rounded(x) = box_width * (bin_number(x) + 0.5)'
end if

comparisonRMSD = 0
comparisonKS = 0
comparisonCDF = 0
do j = 1, scatteringangleBins
        comparisonCDF = comparisonCDF + 1.0*binTotal(j,1) - referenceMeans(j)
        comparisonKS = max(comparisonKS, abs(comparisonCDF))
        comparisonRMSD = comparisonRMSD + (1.0*bintotal(j,1) - referenceMeans(j))**2
end do

write(gnuplotchannel,*) 'set label 1 "'//allprefixes(1:alllengths(1))//'" at graph 0.85, 0.9'
write(variable_length_text,FMT="(F5.2)") comparisonKS
write(gnuplotchannel,*) 'set label 2" KSD: '//variable_length_text//'" at graph 0.85,0.825'
write(variable_length_text,FMT="(F5.2)") sqrt(comparisonRMSD * 1.0 / (scatteringangleBins - 1.0))
write(gnuplotchannel,*) 'set label 3" RMSD: '//variable_length_text//'" at graph 0.85,0.750'
write(variable_length_text,FMT=FMT5_variable) SATRVcolumn
write(gnuplotchannel,*) 'plot "'//gridpath0//allprefixes(1:alllengths(1))//&
                        SATRVfile//'" u (rounded($'//trim(adjustl(variable_length_text))//&
                        ')):(1.0) smooth frequency w boxes, \'
write(gnuplotchannel,*) '     "'//gridpath0//'Adjusted'//SATRVname//cumulativefile//&
                        '.dat" u (scaling*($1)):2 w boxes'//&
                        ' fs transparent solid 0.5 noborder, \'
write(gnuplotchannel,*) '     "'//gridpath0//'Adjusted'//SATRVname//cumulativefile//&
                        '.dat" u (scaling*($1)):2:3 w yerrorbars'
do i = 1, comparison_number-1
comparisonRMSD = 0
comparisonKS = 0
comparisonCDF = 0
do j = 1, scatteringangleBins
        comparisonCDF = comparisonCDF + 1.0*binTotal(j,i+1) - referenceMeans(j)
        comparisonKS = max(comparisonKS, abs(comparisonCDF))
        comparisonRMSD = comparisonRMSD + (1.0*binTotal(j,i+1) - referenceMeans(j))**2
end do

        write(gnuplotchannel,*) 'set label ', 3*i+1, '"'//allprefixes(sum(alllengths(1:i))+1:sum(alllengths(1:i+1)))//&
                                '" at graph 0.85, 0.9'
        write(variable_length_text,FMT="(F5.2)") comparisonKS
        write(gnuplotchannel,*) 'set label ',3*i+2,'" KSD: '//variable_length_text//'" at graph 0.85,0.825'
        write(variable_length_text,FMT="(F5.2)") sqrt(comparisonRMSD * 1.0 / (scatteringangleBins - 1.0))
        write(gnuplotchannel,*) 'set label ',3*i+3,'" RMSD: '//variable_length_text//'" at graph 0.85,0.750'
        write(gnuplotchannel,*) 'unset label ', 3*i-2
        write(gnuplotchannel,*) 'unset label ', 3*i-1
        write(gnuplotchannel,*) 'unset label ', 3*i
        if (i == comparison_number-1) then
                if (SATRVname == "ScatteringAngle") then
                        write(gnuplotchannel,*) 'set xlabel "Scattering Angle (rad)"'
                        write(gnuplotchannel,*) 'set xtics lowerlimit, 10*(upperlimit-lowerlimit)/',SA_Nbins,', upperlimit'
                        write(gnuplotchannel,*) "set format x '%.3P π'"
                else
                        write(gnuplotchannel,*) 'set xlabel "Energy Change (meV)"'
                        write(gnuplotchannel,*) 'set xtics lowerlimit, 10*box_width, upperlimit'
                        write(gnuplotchannel,*) "set format x '%.3f'"
                end if
        end if
        write(variable_length_text,FMT=FMT5_variable) SATRVcolumn
        write(gnuplotchannel,*) 'plot "'//gridpath0//allprefixes(sum(alllengths(1:i))+1:sum(alllengths(1:i+1)))//&
                                SATRVfile//'" u (rounded($'//trim(adjustl(variable_length_text))//&
                                ')):(1.0) smooth frequency w boxes, \'
        write(gnuplotchannel,*) '     "'//gridpath0//'Adjusted'//SATRVname//cumulativefile//&
                                '.dat" u (scaling*($1)):2 w boxes'//&
                                ' fs transparent solid 0.5 noborder, \'
        write(gnuplotchannel,*) '     "'//gridpath0//'Adjusted'//SATRVname//cumulativefile//&
                                '.dat" u (scaling*($1)):2:3 w yerrorbars'
end do
close(gnuplotchannel)

!And then we just input it into gnuplot.exe
call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)

end subroutine getComparedScatteringAngles





subroutine postProcess(prefix_filename)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
use PHYSICS
implicit none

!FILENAMING
character(*),intent(in) :: prefix_filename

!Translational, Rotational, Vibrational Energies
real(dp) :: scatteringAngle, speed_in, speed_out
real(dp) :: bond_distance, rot_speed
real(dp) :: relTranslationalEnergy1, relTranslationalEnergy2
real(dp) :: absTranslationalEnergy1, absTranslationalEnergy2
real(dp) :: RotationalEnergy1, RotationalEnergy2
real(dp) :: TotalEnergy1, TotalEnergy2

!TRAJECTORY DATA
real(dp),dimension(Nbonds,6) :: initial_bonding_data
real(dp),dimension(3,Natoms) :: coords_initial,velocities_initial
real(dp),dimension(3,Natoms) :: coords_final,velocities_final
real(dp),dimension(3) :: velocity_in, velocity_out
real(dp),dimension(3) :: velocity1, velocity2
real(dp),dimension(3) :: coords1, coords2
real(dp),dimension(3) :: velocity1_CM, velocity2_CM

!BONDING DATA INDEXES
integer :: atom1, atom2

!INCREMENTAL INTEGERS
integer :: i, j, k

	open(filechannel1,file=gridpath0//prefix_filename//timeslicefile)
	open(filechannel2,file=gridpath0//prefix_filename//initialfile)
	open(filechannel3,file=gridpath0//prefix_filename//SATRVfile)
	do k = 1, Ntraj
		read(filechannel1,FMTtimeslice) &
                        ((coords_initial(i,j),i=1,3),j=1,Natoms),&
                        ((velocities_initial(i,j),i=1,3),j=1,Natoms),&
                        ((coords_final(i,j),i=1,3),j=1,Natoms),&
                        ((velocities_final(i,j),i=1,3),j=1,Natoms)
		read(filechannel2,FMTinitial) ((initial_bonding_data(i,j),j=1,6),i=1,Nbonds)

		relTranslationalEnergy1 = 0.0d0
		relTranslationalEnergy2 = 0.0d0
		absTranslationalEnergy1 = 0.0d0
		absTranslationalEnergy2 = 0.0d0
		RotationalEnergy1 = 0.0d0
		RotationalEnergy2 = 0.0d0
		TotalEnergy1 = 0.0d0
		TotalEnergy2 = 0.0d0

		velocity1_CM = 0.0d0
		velocity2_CM = 0.0d0
		do i = 1, Natoms
			velocity1_CM = velocity1_CM + velocities_initial(:,i)
			velocity2_CM = velocity2_CM + velocities_final(:,i)
			TotalEnergy1 = TotalEnergy1 + KineticEnergy(velocities_initial(:,i))
			TotalEnergy2 = TotalEnergy2 + KineticEnergy(velocities_final(:,i))
		end do
		velocity1_CM = velocity1_CM / Natoms
		velocity2_CM = velocity2_CM / Natoms

		do i = 1, Nbonds
			atom1 = BONDING_DATA(i,1)
			atom2 = BONDING_DATA(i,2)

			velocity1 = (velocities_initial(:,atom1) + velocities_initial(:,atom2)) / 2
			velocity2 = (velocities_final(:,atom1) + velocities_final(:,atom2)) / 2

			absTranslationalEnergy1 = absTranslationalEnergy1 + 2*KineticEnergy(velocity1)
			absTranslationalEnergy2 = absTranslationalEnergy2 + 2*KineticEnergy(velocity2)
			relTranslationalEnergy1 = relTranslationalEnergy1 + &
                                                  2*KineticEnergy(velocity1-velocity1_CM)
			relTranslationalEnergy2 = relTranslationalEnergy2 + &
                                                  2*KineticEnergy(velocity2-velocity2_CM)

			coords1 = coords_initial(:,atom1)
			coords2 = coords_initial(:,atom2)
			TotalEnergy1 = TotalEnergy1 + HOPotential(coords1,coords2)
			RotationalEnergy1 = RotationalEnergy1 + RotationalEnergy(coords1,coords2,&
                                                                velocities_initial(:,atom1),velocities_initial(:,atom2))
                                           
			coords1 = coords_final(:,atom1)
			coords2 = coords_final(:,atom2)
			TotalEnergy2 = TotalEnergy2 + HOPotential(coords1,coords2)
			RotationalEnergy2 = RotationalEnergy2 + RotationalEnergy(coords1,coords2,&
                                                                velocities_final(:,atom1),velocities_final(:,atom2))
		end do

		velocity_in = 0.0d0
		velocity_out = 0.0d0
		do i = 1, Natoms
			if (all(BONDING_DATA /= i)) then
			velocity1 = velocities_initial(:,i)
			velocity2 = velocities_final(:,i)

			absTranslationalEnergy1 = absTranslationalEnergy1 + KineticEnergy(velocity1)
			absTranslationalEnergy2 = absTranslationalEnergy2 + KineticEnergy(velocity2)
			relTranslationalEnergy1 = relTranslationalEnergy1 + &
                                                  KineticEnergy(velocity1-velocity1_CM)
			relTranslationalEnergy2 = relTranslationalEnergy2 + &
                                                  KineticEnergy(velocity2-velocity2_CM)
			end if

			if (COLLISION_DATA(i) == 1) then
				velocity_in = velocity_in + velocities_initial(:,i)
			else
				velocity_out = velocity_out + velocities_final(:,i)
			end if
		end do
		velocity_in = velocity_in / (sum(COLLISION_DATA))
		velocity_out = velocity_out / (Natoms - sum(COLLISION_DATA))

		speed_out = sqrt(sum((velocity_out)**2))
		speed_in = sqrt(sum((velocity_in)**2))

		scatteringAngle = acos(dot_product(velocity_in,velocity_out) / &
                                  (speed_in * speed_out))

		abs_energychange = abs(absTranslationalEnergy2 - absTranslationalEnergy1)
		rel_energychange = abs(relTranslationalEnergy2 - relTranslationalEnergy1)
		rot_energychange = abs(RotationalEnergy2 - RotationalEnergy1)

		write(filechannel3,FMTdata) scatteringAngle, speed_out, &
					  abs_energychange, &
					  rel_energychange, &
					  rot_energychange

		max_TranslationalEnergy = max(speed_out,max_TranslationalEnergy)

		max_absenergychange = max(abs_energychange,max_absenergychange)
		min_absenergychange = min(abs_energychange,min_absenergychange)
		max_relenergychange = max(rel_energychange,max_relenergychange)
		min_relenergychange = min(rel_energychange,min_relenergychange)
		max_rotenergychange = max(rot_energychange,max_rotenergychange)
		min_rotenergychange = min(rot_energychange,min_rotenergychange)
	end do
	close(filechannel3)
	close(filechannel2)
	close(filechannel1)

end subroutine postProcess


end module analyzeScatteringAngleswithMultipleGrids
