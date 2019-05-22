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
!		This module processes outputs of trjaectories, gains information like the
!		scattering angle and energy change decomposition, and plots them
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILECHANNELS                    ACTION
!
!               GNUPLOTCHANNEL                  OPEN, WRITE, CLOSE
!               TRAJECTORIESCHANNEL             OPEN, WRITE, CLOSE
!               FILECHANNEL1                    OPEN, WRITE/READ, CLOSE
!               FILECHANNEL2                    OPEN, WRITE/READ, CLOSE
!               FILECHANNEL3                    OPEN, WRITE/READ, CLOSE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINES                     ARGUMENTS               KIND
!
!               getScatteringAngles1            prefix_filename         intent(in),character(*)
!                                               PNGfilename             intent(in),character(*)
!
!               getConvergenceImage             lowerlimit              intent(in),real(dp)
!                                               upperlimit              intent(in),real(dp)
!                                               SATRVcolumn             intent(in),integer
!                                               SATRVname               intent(in),character(*)
!
!               getScatteringAngles2            prefix_filename         intent(in),character(*)
!                                               PNGfilename             intent(in),character(*)
!
!               getInitialImages                prefix_filename         intent(in),character(*)
!                                               PNGfilename             intent(in),character(*)
!
!               getComparedScatteringAngles     lowerlimit              intent(in),real(dp)
!                                               upperlimit              intent(in),real(dp)
!                                               imagename               intent(in),character(*)
!                                               SATRVcolumn             intent(in),integer
!                                               SATRVname               intent(in),character(*)
!
!               postProcess                     prefix_filename         intent(in),character(*)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!               gridpath0//prefix_filename//    DAT                             Stores the initial conditions for the set of
!                 initialfile                                                   trajectories associated with this prefix
!               gridpath0//prefix_filename//    DAT                             Stores the initial and final frames of a
!                 timeslicefile                                                 set of trajectories corresponding to the
!                                                                               specified prefix
!               gridpath0//prefix_filename//    DAT                             Stores the scattering angle and energy change
!                 SATRVfile                                                     decomposition for a set of trajectories
!                                                                               corresponding to the specified prefix
!               gridpath0//prefix_filename//    DAT                             Same as the SATRV file but the values are
!                 binnedSATRVfile                                               stored as integers corresponding to some
!                                                                               bin in a specific binning
!               gridpath1//Initial//#traj//     DAT                             Stores the scattering angle and energy change
!                 binnedSATRVfile                                               decomposition for all trajectories in the library;
!                                                                               its is already binned
!               gridpath0//RMSD//SATRVname      DAT                             Stores the averages, devations, and other data
!                 cumulativefile                                                associated with a particular SATRV for all
!                                                                               trajectories in the library for some sampling
!               gridpath0//Adjusted//           DAT                             Similar to cumulative file but normalized and
!                 SATRVname                                                     dependent only on the number of sets and number
!                                                                               of trajectories per set
!               gridpath0//Convergence//        PNG                             The reference distribution of the library for
!                 SATRVname                                                     some sampling and some SATRV; also listed are
!                                                                               the RMSD, KRP, and running average
!               gridpath0//PNGfilename          PNG                             Generic format for image naming
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module analyzeScatteringAngleswithMultipleGrids
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               getScatteringAngles1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine does two things:
!                 1) It reads in frames from the file corresponding to the prefix provided,
!                    proceeding then to bin them and put them in a differently-name filed
!                    with the same prefix
!                 2) It creates a scattering angle heatmap of the set of trajectories
!                    corresponding to the prefix provided
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               prefix_filename                 CHARACTER(*)                    The prefix defining a set of trajectories
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
!               Ntraj                           INTEGER                         The number of trajectories in this
!                                                                               set of trajectories, or the "sampling"
!               energychangeBins                INTEGER                         The number of bins used for distributions
!                                                                               involving energy change
!               energyBins                      INTEGER                         The number of bins used for distributions
!                                                                               involving the final translational energy
!               angleBins                       INTEGER                         The number of bins used for distributions
!                                                                               involving scattering angle
!               angleratio                      INTEGER                         The ratio between the size of the scattering
!                                                                               angle distribution and heatmap;
!                                                                               The heatmap bins should be SMALLER
!               angle_energy_bins               INTEGER,                        The array storing binned scattering angle
!                                               DIM(angleBins,energyBins)       data
!
!               sizeAbsEnergyBin                REAL(DP)                        The size of a bin for the distribution
!                                                                               of absolute translational energy change
!               sizeRelEnergyBin                REAL(DP)                        The size of a bin for the distribution
!                                                                               of relative translational energy change
!               sizeRotEnergyBin                REAL(DP)                        The size of a bin for the distribution
!                                                                               of rotational energy change
!               sizeDeltaEnergyBin              REAL(DP)                        The size of a bin for the distribution
!                                                                               of all kinetic energy changes
!
!               sizeEnergyBin                   REAL(DP)                        The size of a bin for the distribution
!                                                                               of final translational energy
!               sizeAngleBin                    REAL(DP)                        The size of a bin for the distribution
!                                                                               of scattering angles
!               angleslice                      REAL(DP)                        The dividing factor that slices the
!                                                                               scattering angle heatmap into a more
!                                                                               visually pleasing figure
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath0//prefix_filename//    DAT                             Stores the scattering angle and energy change
!                 SATRVfile                                                     decomposition for a set of trajectories
!                                                                               corresponding to the specified prefix
!               gridpath0//prefix_filename//    DAT                             Same as the SATRV file but the values are
!                 binnedSATRVfile                                               stored as integers corresponding to some
!                                                                               bin in a specific binning
!               gridpath0//prefix_filename//    PNG                             Scattering angle and energy change decompositon
!                 PNGfilename                                                   distributions for a set of trajectories
!                                                                               corresponding to the specified prefix
!               gridpath0//prefix_filename//    PNG                             Scattering angle heatmap for a set of
!                 Heatmap_//PNGfilename                                         trajectories corresponding to the specified
!                                                                               prefix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine getScatteringAngles1(prefix_filename,PNGfilename)
use PARAMETERS
use PHYSICS
use ANALYSIS
use FUNCTIONS
implicit none

!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
character(*), intent(in) :: prefix_filename

!FORMAT OF PNG FILES TO BE MADE
character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: PNGfilename

!FORMATTING OF PNG FILES
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

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

!The plots are named starting with Ntraj (the number of trajectories)
write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
write(Ntraj_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ntraj

allocate(angle_energy_bins(angleBins,energyBins))
angle_energy_bins = 0

!The size of the bins, of course, depend on how many bins we have as
!well as the bounds of the distribution
sizeAbsEnergyBin = (max_absenergychange-min_absenergychange) / energychangeBins
sizeRelEnergyBin = (max_relenergychange-min_relenergychange) / energychangeBins
sizeRotEnergyBin = (max_rotenergychange-min_rotenergychange) / energychangeBins

sizeDeltaEnergyBin = (max(max_absenergychange,max_relenergychange,max_rotenergychange) - &
                      min(min_absenergychange,min_relenergychange,min_rotenergychange)) / energychangeBins
sizeEnergyBin = max_TranslationalEnergy / (energyBins)
sizeAngleBin = pi / (angleBins)

!If we are doing this for a comparison, we may have user-defined bounds
if (comparison_flag .and. (comparison_upperlimit /= comparison_lowerlimit)) then

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

!Now we can actually bin them
!Here we open up the SATRVfile (which has the observables) and the binnedSATRVfile
!(which has the observables after binning)
open(filechannel1,file=gridpath5//prefix_filename//SATRVfile)
open(filechannel2,file=gridpath5//prefix_filename//binnedSATRVfile)
do i = 1, Ntraj

        !First, we read the line
        read(filechannel1,FMT=FMTdata,iostat=iostate) ScatteringAngle_real, TranslationalEnergy_real, &
						      abs_energychange, rel_energychange, rot_energychange

        !Then we figure out what bin it should be in
	ScatteringAngle = ceiling(ScatteringAngle_real / sizeAngleBin)
	TranslationalEnergy = ceiling((TranslationalEnergy_real)/ sizeEnergyBin)
	AbsEnergyChange = ceiling((abs_energychange)/sizeAbsEnergyBin)
	RelEnergyChange = ceiling((rel_energychange)/sizeRelEnergyBin)
	RotEnergyChange = ceiling((rot_energychange)/sizeRotEnergyBin)

        !If the bin is out of bounds, put it back in bounds
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

        !For the scattering angle heatmap, we also store these values in another array
	angle_energy_bins(ScatteringAngle,TranslationalEnergy) = &
                   angle_energy_bins(ScatteringAngle,TranslationalEnergy) + 1

        !We have two separate bin sizes for the scattering angle heatmap and the
        !scattering angle distribution. The former have SMALLER bins
        !Thus, we can get the binning of the latter by dividing by some factor which
        !I call "angle_ratio": the ratio between the sizes of these two binnings
	write(filechannel2,FMT=*) ceiling((ScatteringAngle-0.5)/angle_ratio),&
                                  TranslationalEnergy,AbsEnergyChange,RelEnergyChange,RotEnergyChange
end do
close(filechannel1)
close(filechannel2)


!This is the gnuplot code to make the scattering angle and energy change (SATRV) plots
!This is only plotted if there is no comparison active

if (comparison_flag) then
        deallocate(angle_energy_bins)
        return
end if

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath4//PNGfilename//'.png"'
write(gnuplotchannel,*) 'set multiplot layout 4,1 margins 0.10,0.95,.1,.95 spacing 0,0.1'//&
                        'title "Scattering Angle and '//&
                        'Energy Change Distributions of '//trim(adjustl(Ntraj_text))//&
                        ' Trajectories" font ",18" offset 0,3'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'set ylabel "Occurence"'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'unset key'

write(gnuplotchannel,*) 'set xlabel "Scattering Angle (rad)"'
write(gnuplotchannel,*) 'pi = 3.14159265'
write(gnuplotchannel,*) 'box_width = ', sizeAngleBin*angle_ratio
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set xtics pi/2'
write(gnuplotchannel,*) "set format x '%.1P π'"
write(gnuplotchannel,*) 'set xrange [0:pi]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'plot "'//gridpath5//prefix_filename//binnedSATRVfile//&
                        '" u (box_width*($1-0.5)):(1.0) smooth frequency with boxes'

write(gnuplotchannel,*) 'scaling = 1000'

write(gnuplotchannel,*) 'set xlabel "Absolute Translation Energy Change (meV)"'
write(gnuplotchannel,*) 'min_E = scaling * ', min_absenergychange
write(gnuplotchannel,*) 'max_E = scaling * ', max_absenergychange
write(gnuplotchannel,*) 'box_width = (max_E-min_E) /', energychangeBins
write(gnuplotchannel,*) 'set xrange [min_E:max_E]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set xtics min_E, box_width * 10, max_E'
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) "set format x '%.3f'"
write(gnuplotchannel,*) 'plot "'//gridpath5//prefix_filename//binnedSATRVfile//&
                        '" u (box_width*($3-0.5)+min_E):(1.0) smooth frequency w boxes'

write(gnuplotchannel,*) 'set xlabel "Relative Translational Energy Change (meV)"'
write(gnuplotchannel,*) 'min_E = scaling * ', min_relenergychange
write(gnuplotchannel,*) 'max_E = scaling * ', max_relenergychange
write(gnuplotchannel,*) 'box_width = (max_E-min_E) /', energychangeBins
write(gnuplotchannel,*) 'set xrange [min_E:max_E]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set xtics min_E, box_width * 10, max_E'
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'plot "'//gridpath5//prefix_filename//binnedSATRVfile//&
                        '" u (box_width*($4-0.5)+min_E):(1.0) smooth frequency w boxes'

write(gnuplotchannel,*) 'set xlabel "Rotational Energy Change (meV)"'
write(gnuplotchannel,*) 'min_E = scaling * ', min_rotenergychange
write(gnuplotchannel,*) 'max_E = scaling * ', max_rotenergychange
write(gnuplotchannel,*) 'Nbins = ', energychangeBins
write(gnuplotchannel,*) 'box_width = (max_E-min_E) / Nbins'
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'set xrange [min_E:max_E]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set xtics min_E, box_width * 10, max_E'
write(gnuplotchannel,*) 'plot "'//gridpath5//prefix_filename//binnedSATRVfile//&
                        '" u (box_width*($5-0.5)+min_E):(1.0) smooth frequency w boxes'
close(gnuplotchannel)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

!For the scattering angle - final translational energy heatmap, all we need
!to do is write the values of angle_energy_bins onto a file and plot it
!For convenience, I also record the maximum frequency in the array
occurence_max = maxval(angle_energy_bins)
open(filechannel1,file=gridpath5//temporaryfile1)
do i = 1, angleBins/angle_slice+1
	do j = 1, energyBins
		write(filechannel1,FMT=*) (i-1)*sizeAngleBin, (j-1)*sizeEnergyBin, angle_energy_bins(i,j)
	end do
	write(filechannel1,FMT=*) ""
end do
close(filechannel1)
deallocate(angle_energy_bins)

!This is the gnuplot code to make the scattering angle - final translational energy heatmap
!This took A LOT of fine-tuning to make it look nice, so please don't change it
open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
write(gnuplotchannel,*) 'set output "'//gridpath4//'HeatMap_'//&
                        PNGfilename//'.png"'
write(gnuplotchannel,*) 'set title "Scattering Angle Distribution" font ",18" offset 0,-20'
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
write(gnuplotchannel,*) 'splot "'//gridpath5//temporaryfile1//'" u (scaling_factor*$2*cos($1)):'//&
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
call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)


end subroutine getScatteringAngles1







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               getConvergenceImage
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine takes some binning of the library's trajectories and finds
!               how it converges; it produces error bars along the way which are stored
!               in a separate file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               lowerlimit                      REAL(DP)                        The lower bound of the distribution
!               upperlimit                      REAL(DP)                        The upper bound of the distribution
!               SATRVcolumn                     INTEGER                         The column in the SATRV file this distribution
!                                                                               is recorded in
!               SATRVname                       CHARACTER(*)                    The name of the distribution (needs to be exact)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               Ntesttraj                       INTEGER                         The number of trajectories per sample
!                                                                               set of trajectories, or the "sampling"
!               scatteringangleBins             INTEGER                         The number of bins used for distributions
!                                                                               (all of them, in this case)
!               bin_width                       REAL(DP)                        The size of the binning
!               Nsamples_max                    INTEGER                         The number of samples the library has
!                                                                               in total
!
!               binTotal                        INTEGER,DIM(                    The size of each bin of the distribution
!                                               scatteringangleBins,            for every sample
!                                               Nsamples_max)
!               binCumulative                   INTEGER,DIM(                    The size of each bin of the cumulative
!                                               scatteringangleBins,            distribution for every sample; this is
!                                               Nsamples_max)                   used in the Kolmogorov-Smirnov difference
!               sampleSize                      INTEGER,DIM(                    The number of trajectories per sample;
!                                               Nsamples_max)                   (should be constant)
!
!               binAverage                      REAL,DIM(                       The average per bin over all samples;
!                                               scatteringangleBins)            this is used in the final distribution
!               binSD                           REAL,DIM(                       The standard deviation per bin over all samples;
!                                               scatteringangleBins)            this is used in the final distribution
!
!               sampleKS                        REAL,DIM(                       The Kolmogorov-Smirnov difference per
!                                               Nsamples_max)                   sample
!               sampleRMSD                      REAL,DIM(                       The root mean square difference per
!                                               Nsamples_max)                   sample
!               sampleKRP                       REAL,DIM(                       The KRP per sample
!                                               Nsamples_max)
!
!               lambda_penalty                  REAL                            The value of lambda used to
!                                                                               calculate the KRP
!               minsd_penalty                   REAL                            The minimum standard deviation used to
!                                                                               calculate the KRP
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath1//Initial//#traj//     DAT                             Stores the scattering angle and energy change
!                 binnedSATRVfile                                               decomposition for all trajectories in the library;
!                                                                               its is already binned
!               gridpath0//RMSD//SATRVname      DAT                             Stores the averages, devations, and other data
!                 cumulativefile                                                associated with a particular SATRV for all
!                                                                               trajectories in the library for some sampling
!               gridpath0//Adjusted//           DAT                             Similar to cumulative file but normalized and
!                 SATRVname                                                     dependent only on the number of sets and number
!                                                                               of trajectories per set
!               gridpath0//Convergence//        PNG                             The reference distribution of the library for
!                 SATRVname                                                     some sampling and some SATRV; also listed are
!                                                                               the RMSD, KRP, and running average
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
real :: bin_width, binMean, binRMSD
real,allocatable :: binAverage(:), binSD(:),sampleKS(:), sampleRMSD(:), sampleKRP(:)
integer,allocatable :: binCumulative(:,:), binTotal(:,:),sampleSize(:)
integer :: binTally
real :: binThreshold = 1.0

!FORMATTING OF PNG FILES
character(5) :: variable_length_text
character(5) :: variable_length_text1, variable_length_text2
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectory_text_length) :: Ntraj_text

character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5

!INCREMENTAL INTEGER
integer :: i

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

!Some initialization
bin_width = (upperlimit-lowerlimit) / scatteringangleBins
Nsamples_max = (Ngrid_max * Ntraj_max) / Ntesttraj
allocate(binCumulative(scatteringangleBins,Nsamples_max),&
         binTotal(scatteringangleBins,Nsamples_max),&
         binAverage(scatteringangleBins),binSD(scatteringangleBins),&
         sampleSize(Nsamples_max),sampleKS(Nsamples_max),&
         sampleKRP(Nsamples_max),sampleRMSD(Nsamples_max))

!The minimum standard deviation (used to calculate KRP) of a bin is
!set to be the theoretical standard deviation if a value for the
!observable occurs only once over the whole set. Thus:
!
! minBin = 1 / (Nsamples_max * Ntesttraj)
!
!  minSD = sqrt( SUM( (minBin - valueBin)^2 ) /
!                (Nsamples_max - 1)                   )
!
!        = sqrt( ((Nsamples_max - 1) * (minBin - 0)^2
!                        + (1) * (minBin - 1)^2 )) /
!                (Nsamples_max - 1)                   )
!        = sqrt( ((Nsamples_max - 1) * (1 / (Nsamples_max * Ntesttraj)^2)
!                        + ( 1 - 1/ (Nsamples_max * Ntesttraj))^2 ) /
!                (Nsamples_max - 1)                   )
!        = sqrt( ((Nsamples_max - 1) / (Nsamples_max * Ntesttraj)^2
!                        + ((Nsamples_max * Ntesttraj) - 1)^2 / (Nsamples_max * Ntesttraj)^2) /
!                (Nsamples_max - 1)                   )
!        = sqrt( ((Nsamples_max * Ntesttraj - 1)^2 + Nsamples_max - 1) /
!                (Nsamples_max - 1)                   ) / (Nsamples_max * Ntesttraj)
!
minsd_penalty = sqrt( ((Nsamples_max * Ntesttraj - 1)**2 + Nsamples_max - 1) * 1.0 / &
                      (Nsamples_max - 1) ) * 1.0 / (Nsamples_max * Ntesttraj)

!And the lambda used in the KRP is just set to be 2
!to mimic the RMSD
lambda_penalty = 2.0

!Now we want a reference distribution
!This will be the distribution of ALL trajectories
binTotal = 0
binCumulative = 0
sampleSize = 0
write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
write(Ntraj_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid_max * Ntraj_max
write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid_max
open(filechannel1,file=gridpath5//Ntraj_text//binnedSATRVfile,action="read")
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
        binSD(i) = sqrt(sum((binTotal(i,:)*1.0 - binAverage(i))**2)/(Nsamples_max - 1))
end do

!Now we want error bars
!For this, we will get a running average distribution by adding Ntesttraj trajectories at a time
!When the running average converges then we stop and see the variance
binTally = 0
sampleKS = 0
sampleRMSD = 0
sampleKRP = 0
open(filechannel1,file=gridpath5//"RMSD"//SATRVname//cumulativefile//".dat",action="write")
do Nsamples = 1, Nsamples_max
        binRMSD = 0.0
        do i = 1, scatteringangleBins
                binRMSD = binRMSD + (sum(binTotal(i,1:Nsamples)) * 1.0 / Nsamples - binAverage(i))**2
                sampleRMSD(Nsamples) = sampleRMSD(Nsamples) + (binTotal(i,Nsamples) - binAverage(i))**2
                binCumulative(i,Nsamples) = sum(binTotal(1:i,Nsamples))
                sampleKS(Nsamples) = max(sampleKS(Nsamples),abs(binCumulative(i,Nsamples)- &
                                         sum(binAverage(1:i))))
                sampleKRP(Nsamples) = sampleKRP(Nsamples) + (abs(binTotal(i,Nsamples) - binAverage(i)) /&
                                                            max(minsd_penalty,binSD(i)))**lambda_penalty
        end do
        sampleRMSD(Nsamples) = sqrt(sampleRMSD(Nsamples) * 1.0 / (scatteringangleBins - 1.0))
        sampleKRP(Nsamples) = (sampleKRP(Nsamples) * 1.0 / (scatteringangleBins - 1.0))**(1.0/lambda_penalty)

        binRMSD = sqrt(binRMSD/scatteringangleBins)
        write(filechannel1,FMT=*) Nsamples*Ntesttraj, binRMSD, binThreshold,&
                                  sampleKS(Nsamples), sampleRMSD(Nsamples), sampleKRP(Nsamples)

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
open(filechannel1,file=gridpath5//"Adjusted"//SATRVname//cumulativefile//".dat")
do i = 1, scatteringangleBins
        write(filechannel1,*) (i-0.5)*bin_width, binAverage(i), binSD(i)
end do
close(filechannel1)

!We have everything we need to draw the distribution with error bars
open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) "set terminal pngcairo size 1200,1200"
write(variable_length_text,FMT="(I5)") Ntesttraj
write(gnuplotchannel,*) 'set output "'//gridpath4//'Convergence'//&
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
write(gnuplotchannel,*) 'plot "'//gridpath5//'RMSD'//SATRVname//cumulativefile//'.dat" u 1:2 w lines, \'
write(gnuplotchannel,*) '     "'//gridpath5//'RMSD'//SATRVname//cumulativefile//'.dat" u 1:3 w lines lc -1'
write(gnuplotchannel,*) 'unset label 1'
write(variable_length_text,FMT="(I5)") Ntesttraj
write(gnuplotchannel,*) 'set title "Final '//SATRVname//' Distribution for '//&
                        trim(adjustl(variable_length_text))//' Trajectories"'
if (SATRVname == "ScatteringAngle") then
write(gnuplotchannel,*) 'scaling = 1'
write(gnuplotchannel,*) 'pi = 3.14159265'
write(gnuplotchannel,*) 'lowerlimit = ', lowerlimit
write(gnuplotchannel,*) 'upperlimit = ', upperlimit
write(gnuplotchannel,*) 'set xtics lowerlimit, 10*(upperlimit-lowerlimit)/',&
                        scatteringangleBins,', upperlimit'
write(gnuplotchannel,*) "set format x '%.3P π'"
write(gnuplotchannel,*) 'set xrange [lowerlimit:upperlimit]'
write(gnuplotchannel,*) 'set xlabel "Scattering Angle (rad)"'
else
write(gnuplotchannel,*) 'scaling = 1000'
write(gnuplotchannel,*) 'E_min = scaling * ', lowerlimit
write(gnuplotchannel,*) 'E_max = scaling * ', upperlimit
write(gnuplotchannel,*) 'set xrange [E_min:E_max]'
write(gnuplotchannel,*) 'set xtics E_min, 10*(E_max-E_min)/',scatteringangleBins,', E_max'
write(gnuplotchannel,*) "set format x '%.3f'"
write(gnuplotchannel,*) 'set xlabel "Energy (meV)"'
end if
write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set ytics autofreq'
write(gnuplotchannel,*) 'set boxwidth'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set ylabel "Frequency"'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'plot "'//gridpath5//'Adjusted'//SATRVname//cumulativefile//'.dat" u (scaling*($1)):2 w boxes, \'
write(gnuplotchannel,*) '     "'//gridpath5//'Adjusted'//SATRVname//cumulativefile//'.dat" u (scaling*($1)):2:3 w yerrorbars'
write(gnuplotchannel,*) 'set title "Distribution of the KRP Among\n'//&
                        trim(adjustl(variable_length_text))//' '//SATRVname//' Samplings Across '//&
                        trim(adjustl(Ngrid_text))//' Grids"'
write(gnuplotchannel,*) 'set xlabel "KRP"'
write(gnuplotchannel,*) 'set ylabel "Frequency"'
write(gnuplotchannel,*) 'minKRP = ', max(0.0,minval(sampleKRP)-2*(maxval(sampleKRP)-minval(sampleKRP))/Nsamples_max)
write(gnuplotchannel,*) 'maxKRP = ', maxval(sampleKRP) + 2*(maxval(sampleKRP)-minval(sampleKRP))/Nsamples_max
write(gnuplotchannel,*) 'set xrange [minKRP:maxKRP]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'box_width = (maxKRP-minKRP) / ',(2*Nsamples_max)
write(gnuplotchannel,*) 'set xtics minKRP, (maxKRP-minKRP)/4, maxKRP'
write(gnuplotchannel,*) "set format x '%.2f'"
write(gnuplotchannel,*) 'set ytics 1'
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'bin_number(x) = floor(x/box_width)'
write(gnuplotchannel,*) 'rounded(x) = box_width * (bin_number(x) + 0.5)'
write(gnuplotchannel,*) 'plot "'//gridpath5//'RMSD'//SATRVname//cumulativefile//'.dat"'//&
                        'u (rounded($6)):(1.0) smooth frequency with boxes'
write(gnuplotchannel,*) 'set title "Distribution of the Root Mean Square Difference Among\n'//&
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
write(gnuplotchannel,*) 'plot "'//gridpath5//'RMSD'//SATRVname//cumulativefile//'.dat"'//&
                        'u (rounded($5)):(1.0) smooth frequency with boxes'
close(gnuplotchannel)

deallocate(binAverage,binSD)
deallocate(binCumulative)
deallocate(binTotal)
deallocate(sampleKS,sampleRMSD,sampleSize,sampleKRP)

call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

end subroutine getConvergenceImage




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               getScatteringAngles2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine reads in frames from the file corresponding to the prefix provided,
!               and then proceeds to visualize them in a distribution with a reference distribution
!               if provided with one
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               prefix_filename                 CHARACTER(*)                    The prefix defining a set of trajectories
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
!               Ntraj                           INTEGER                         The number of trajectories in this
!                                                                               set of trajectories, or the "sampling"
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath0//prefix_filename//    DAT                             Stores the scattering angle and energy change
!                 SATRVfile                                                     decomposition for a set of trajectories
!                                                                               corresponding to the specified prefix
!               gridpath0//Adjusted//           DAT                             Stores the averages and deviations for the 
!                 SATRVname                                                     distribution of some SATRV for all the
!                                                                               trajectories in the library for some sampling
!               gridpath0//PNGfilename          PNG                             Distributions of scattering angle and energy
!                                                                               change decomposition for some set of
!                                                                               trajectories
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine getScatteringAngles2(prefix_filename,PNGfilename)
use PARAMETERS
use PHYSICS
use ANALYSIS
implicit none

!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
character(*), intent(in) :: prefix_filename

!FORMAT OF PNG FILES TO BE MADE
character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: PNGfilename

!FORMATTING OF PNG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectory_text_length) :: Ntraj_text
character(6) :: boxwidth_text
logical :: grid_is_done

integer :: iostate
real :: speed_out, ScatteringAngle

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

inquire(file=gridpath5//'AdjustedScatteringAngle'//cumulativefile//'.dat',exist=grid_is_done)

!This is the gnuplot code to make the plots
open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo enhanced size 1200,1200'
write(gnuplotchannel,*) 'set encoding utf8'
write(gnuplotchannel,*) 'set output "'//gridpath4//PNGfilename//'.png"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'pi = 3.14159265'
write(gnuplotchannel,*) 'Nbins = ', energychangeBins
write(gnuplotchannel,*) 'scaling = 1'
write(gnuplotchannel,*) 'box_width = pi / Nbins'
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'bin_number(x) = floor(x/box_width)'
write(gnuplotchannel,*) 'rounded(x) = box_width * (bin_number(x) + 0.5)'
write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
write(Ntraj_text,FMT="(I0."//variable_length_text//")") Ntraj
write(gnuplotchannel,*) 'set multiplot layout 4,1 margins 0.10,0.95,.1,.95 spacing 0,0.1'//&
                        'title "Scattering Angle and '//&
                        'Energy Change Distributions of '//trim(adjustl(Ntraj_text))//&
                        ' Trajectories" font ",18" offset 0,3'
write(gnuplotchannel,*) 'set style histogram clustered gap 1'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'set ylabel "Occurence"'
write(gnuplotchannel,*) 'set yrange [0:]'

write(gnuplotchannel,*) 'set xlabel "Scattering Angle (rad)"'
write(gnuplotchannel,*) 'set xrange [0:pi]'
write(gnuplotchannel,*) 'set xtics pi/2'
write(gnuplotchannel,*) "set format x '%.1P π'"

if (grid_is_done) then
        write(gnuplotchannel,*) 'plot "'//gridpath5//SATRVfile//&
                                '" u (scaling*$1>=pi?(pi-0.5*box_width):'//&
                                '(rounded($1))):(1.0) smooth frequency w boxes, \'
        write(gnuplotchannel,*) '     "'//gridpath5//'AdjustedScatteringAngle'//&
				cumulativefile//'.dat" u 1:2 w boxes'//&
                                       ' fs transparent solid 0.5 noborder, \'
        write(gnuplotchannel,*) '     "'//gridpath5//'AdjustedScatteringAngle'//&
				cumulativefile//'.dat" u 1:2:3 w yerrorbars'
else
        write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//SATRVfile//&
                                '" u (scaling*$1>=pi?(pi-0.5*box_width):'//&
                                '(rounded($1))):(1.0) smooth frequency w boxes'
end if
        write(gnuplotchannel,*) 'scaling = 1000'

        write(gnuplotchannel,*) 'set xlabel "Absolute Translational Energy Change (meV)"'
        write(gnuplotchannel,*) 'min_E = scaling * ', min_absenergychange
        write(gnuplotchannel,*) 'max_E = scaling * ', max_absenergychange
        write(gnuplotchannel,*) 'Nbins = ', energychangeBins
        write(gnuplotchannel,*) 'box_width = (max_E-min_E) / Nbins'
        write(gnuplotchannel,*) 'set xrange [min_E:max_E]'
        write(gnuplotchannel,*) 'set xtics min_E, box_width*Nbins/5, max_E'
        write(gnuplotchannel,*) 'set boxwidth box_width'
	write(gnuplotchannel,*) 'bin_number(x) = floor(scaling*x/box_width)'
        write(gnuplotchannel,*) 'rounded(x) = min_E + box_width * (bin_number(x) + 0.5)'
        write(gnuplotchannel,*) "set format x '%.3f'"
if (grid_is_done) then
        write(gnuplotchannel,*) 'plot "'//gridpath5//SATRVfile//&
                                '" u (scaling*$3>=max_E?(max_E-0.5*box_width):'//&
                                '(rounded($3))):(1.0) smooth frequency w boxes, \'
        write(gnuplotchannel,*) '     "'//gridpath5//'AdjustedAbsoluteEnergyChange'//&
				cumulativefile//'.dat" u (scaling*($1)):2 w boxes'//&
                                       ' fs transparent solid 0.5 noborder, \'
        write(gnuplotchannel,*) '     "'//gridpath5//'AdjustedAbsoluteEnergyChange'//&
				cumulativefile//'.dat" u (scaling*($1)):2:3 w yerrorbars'
else
        write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//SATRVfile//&
                                '" u (scaling*$3>=max_E?(max_E-0.5*box_width):'//&
                                '(rounded($3))):(1.0) smooth frequency w boxes'
end if
        write(gnuplotchannel,*) 'set xlabel "Relative Translational Energy Change (meV)"'
        write(gnuplotchannel,*) 'min_E = scaling * ', min_relenergychange
        write(gnuplotchannel,*) 'max_E = scaling * ', max_relenergychange
        write(gnuplotchannel,*) 'Nbins = ', energychangeBins
        write(gnuplotchannel,*) 'box_width = (max_E-min_E) / Nbins'
        write(gnuplotchannel,*) 'set xrange [min_E:max_E]'
        write(gnuplotchannel,*) 'set xtics min_E, box_width*Nbins/5, max_E'
        write(gnuplotchannel,*) 'set boxwidth box_width'
	write(gnuplotchannel,*) 'bin_number(x) = floor(scaling*x/box_width)'
        write(gnuplotchannel,*) 'rounded(x) = min_E + box_width * (bin_number(x) + 0.5)'
if (grid_is_done) then
        write(gnuplotchannel,*) 'plot "'//gridpath5//SATRVfile//&
                                '" u (scaling*$4>=max_E?(max_E-0.5*box_width):'//&
                                '(rounded($4))):(1.0) smooth frequency w boxes, \'
        write(gnuplotchannel,*) '     "'//gridpath5//'AdjustedRelativeEnergyChange'//&
				cumulativefile//'.dat" u (scaling*($1)):2 w boxes'//&
                                       ' fs transparent solid 0.5 noborder, \'
        write(gnuplotchannel,*) '     "'//gridpath5//'AdjustedRelativeEnergyChange'//&
				cumulativefile//'.dat" u (scaling*($1)):2:3 w yerrorbars'
else
        write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//SATRVfile//&
                                '" u (scaling*$4>=max_E?(max_E-0.5*box_width):'//&
                                '(rounded($4))):(1.0) smooth frequency w boxes'
end if
        write(gnuplotchannel,*) 'set xlabel "Rotational Energy Change (meV)"'
        write(gnuplotchannel,*) 'min_E = scaling * ', min_rotenergychange
        write(gnuplotchannel,*) 'max_E = scaling * ', max_rotenergychange
        write(gnuplotchannel,*) 'Nbins = ', energychangeBins
        write(gnuplotchannel,*) 'box_width = (max_E-min_E) / Nbins'
        write(gnuplotchannel,*) 'set xrange [min_E:max_E]'
        write(gnuplotchannel,*) 'set xtics min_E, box_width*Nbins/5, max_E'
        write(gnuplotchannel,*) "set format x '%.3f'"
        write(gnuplotchannel,*) 'set boxwidth box_width'
	write(gnuplotchannel,*) 'bin_number(x) = floor(scaling*x/box_width)'
        write(gnuplotchannel,*) 'rounded(x) = min_E + box_width * (bin_number(x) + 0.5)'
if (grid_is_done) then
        write(gnuplotchannel,*) 'plot "'//gridpath5//SATRVfile//&
                                '" u (scaling*$5>=max_E?(max_E-0.5*box_width):'//&
                                '(rounded($5))):(1.0) smooth frequency w boxes, \'
        write(gnuplotchannel,*) '     "'//gridpath5//'AdjustedRotationalEnergyChange'//&
				cumulativefile//'.dat" u (scaling*($1)):2 w boxes'//&
                                       ' fs transparent solid 0.5 noborder, \'
        write(gnuplotchannel,*) '     "'//gridpath5//'AdjustedRotationalEnergyChange'//&
				cumulativefile//'.dat" u (scaling*($1)):2:3 w yerrorbars'
else
        write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//SATRVfile//&
                                '" u (scaling*$5>=max_E?(max_E-0.5*box_width):'//&
                                '(rounded($5))):(1.0) smooth frequency w boxes'
end if
close(gnuplotchannel)

!And then we just input it into gnuplot.exe
call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

end subroutine getScatteringAngles2




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               getInitialImages
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine reads in the initial conditions of all the trajectories of
!               some set of trajectories and plots their distribution
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               prefix_filename                 CHARACTER(*)                    The prefix defining a set of trajectories
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
!               Ntraj                           INTEGER                         The number of trajectories in this
!                                                                               set of trajectories
!
!               max_r0                          REAL                            The maximum bond distance encountered
!                                                                               in this set of trajectories
!               min_r0                          REAL                            The minimum bond distance encountered
!                                                                               in this set of trajectories
!               average_r0                      REAL                            The average bond distance encountered
!                                                                               in this set of trajectories
!               max_rot                         REAL                            The maximum rotational speed encountered
!                                                                               in this set of trajectories
!               min_rot                         REAL                            The minimum rotational speed encountered
!                                                                               in this set of trajectories
!               average_rot                     REAL                            The average rotational speed encountered
!                                                                               in this set of trajectories
!
!               average_Evib                    REAL                            The average vibrational energy encountered
!                                                                               in this set of trajectories
!               average_Erot                    REAL                            The average rotational energy encountered
!                                                                               in this set of trajectories
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath0//prefix_filename//    DAT                             Stores the initial conditions for the set of
!                 initialfile                                                   trajectories associated with this prefix
!               gridpath0//PNGfilename          PNG                             Distributions of the initial conditions for the
!                                                                               set of trajectories associated with this prefix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine getInitialimages(prefix_filename,PNGfilename)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
use PHYSICS
implicit none

!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
character(*),intent(in) :: prefix_filename
character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5
character(*),intent(in) :: PNGfilename

!FORMATTING OF PNG FILES
character(5) :: variable_length_text
character(5) :: variable_length_text1, variable_length_text2
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(6) :: Ntraj_text
character(150) :: old_filename

integer,parameter :: initialBins = 50

integer :: iostate
real :: r0, rot
real :: max_r0, max_rot
real :: min_r0, min_rot
real :: average_r0, average_rot
real :: average_Evib, average_Erot
integer :: total_bonds
integer :: i

real :: N_upsilon, N_J
real :: expected_upsilon
real :: expected_evib
real :: expected_temperature

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

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
!	read(filechannel1,iostat=iostate,FMT=FMTinitial) INITIAL_BOND_DATA
	read(filechannel1,iostat=iostate,FMT=*) INITIAL_BOND_DATA
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

!The average bond length
average_r0 = average_r0 / total_bonds

!The average linear rotational speed
average_rot = average_rot / total_bonds

!The average vibrational energy
average_Evib = 0.5 * HOke_hydrogen * average_Evib / total_bonds

!The average rotational energy
average_Erot = average_Erot * mass_hydrogen / total_bonds

open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo enhanced size 3600,1200'
write(gnuplotchannel,*) 'set encoding utf8'
write(gnuplotchannel,*) 'set output "'//gridpath4//PNGfilename//'.png"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'pi = 3.14159265'
write(gnuplotchannel,*) 'set style histogram clustered gap 1'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'

!
!The average bond length is displayed on screen
!
write(Ntraj_text,FMT="(F6.4)") average_r0
write(gnuplotchannel,*) 'set label 1 "Average r0: '//Ntraj_text//' A" at screen 0.65,0.925'
write(gnuplotchannel,*) 'set label 1 front'

!
!The (vibrational) temperature is displayed on screen
!
!For the vibrational partition function, we can define a relationship between
!temperature and average vibrational energy
!McQuarrie, 1997 : section 18-4 : pg.741
!
! average_Evib = (kB) (theta_vib) ( (1/2) + (e^(theta_vib / T) - 1)^-1)                    (w/ ZPE)
!
! avearge_Evib = (kb) (theta_vib) / (e^(theta_vib / T) - 1)                                (w/o ZPE)
!
!Rearranging:
!
! T = theta_vib / log(1 + (2) (kB) (theta_vib) / ((2) (average_Evib) - (kB) (theta_vib)))  (w/ ZPE)
!
! T = theta_vib / log(1 + (kB) (theta_vib) / average_Evib)                                 (w/o ZPE)
!
!Thus, we can get the two equations below:
!
! Without ZPE:
!write(Ntraj_text,FMT="(F6.1)") theta_vib / log(1.0d0 + (2.0d0*(kb/RU_energy)*theta_vib) / &
!                                               (2.0d0*average_Evib ))
! With ZPE:
write(Ntraj_text,FMT="(F6.1)") theta_vib / log(1.0d0 + (2.0d0*(kb/RU_energy)*theta_vib) / &
                                               (2.0d0*average_Evib - (kb/RU_energy)*theta_vib))
write(gnuplotchannel,*) 'set label 2 "Temperature: '//Ntraj_text//' K" at screen 0.65,0.9'
write(gnuplotchannel,*) 'set label 2 front'

!
! BUG TESTING START
!print *, ""
!print *, "kb*theta_vib:", (kb/RU_energy)*theta_vib
!print *, "hv/2: ", 0.5d0*(hbar*pi2/(RU_energy*RU_time))*vib_frequency 
!print *, "evib: ", average_Evib
!print *, ""
!print *, "Ptotal:", upsilon_factor2*(1.0d0-exp(-upsilon_max*upsilon_factor1))/upsilon_factor1
!print *, "up_f2/up_f1:", upsilon_factor2/upsilon_factor1
!print *, "Pupsilon:", (-upsilon_max*exp(-upsilon_max*upsilon_factor1) + (1.0d0-exp(-upsilon_max*upsilon_factor1)) / &
!                          upsilon_factor1) * (upsilon_factor2 / upsilon_factor1)
!print *, "evib-hv/2: ", average_Evib - 0.50d0*(hbar*pi2/(RU_energy*RU_time))*vib_frequency 
!print *, "evib+hv/2: ", average_Evib + 0.50d0*(hbar*pi2/(RU_energy*RU_time))*vib_frequency 
!print *, "theta_vib: ", theta_vib
!print *, "vib_freq: ", vib_frequency
! BUG TESTING END
!

!
!Obtain the expectation value of the vibrational quantum
!number by integrating over the domain:
!
! <upsilon> = [ integral{0 -> upsilon_max}
!                 (upsilon) P(upsilon) d(upsilon) ] / N
!
!           = [ integral{0 -> upsilon_max}
!                 (upsilon) (upsilon_factor2) e^(-(upsilon) (upsilon_factor1)) d(upsilon) ] / N
!           = (upsilon_factor2 / N) [ integral{0 -> upsilon_max}
!                 (upsilon) e^(-(upsilon) (upsilon_factor1)) d(upsilon) ]
!           = (upsilon_factor2 / N) (
!                 [ (upsilon) (-upsilon_factor1)^-1 exp(-(upsilon) (upsilon_factor1)) ]{0 -> upsilon_max} -
!                 [ integral{0 -> upsilon_max} (-upsilon_factor)^-1 exp(-(upsilon) (upsilon_factor1)) d(upsilon) ]
!                                   )
!           = (upsilon_factor2 / (upsilon_factor1) (N)) (
!                 [ 0 - (upsilon_max) exp(-(upsilon_max) (upsilon_factor1)) ] +
!                 [ (-upsilon_factor1)^-1 exp(-(upsilon) (upsilon_factor1)) ]{0 -> upsilon_max}
!                                                       )
!           = (upsilon_factor2 / (upsilon_factor1) (N)) (
!                 -(upsilon_max) exp(-(upsilon_max) (upsilon_factor1)) -
!                 (upsilon_factor1)^-1 (exp(-(upsilon_max) (upsilon_factor1)) - 1)
!                                                       )
!           = (upsilon_factor2 / (upsilon_factor1) (N)) (
!                 -(upsilon_max) exp(-(upsilon_max) (upsilon_factor1)) +
!                 (upsilon_factor1)^-1 (1 - exp(-(upsilon_max) (upsilon_factor1)))
!                                                       )
!
!where N is the normalization constant of the integral:
!
!         N = integral{0 -> upsilon_max}
!                           P(upsilon) d(upsilon)
!
!           = [ integral{0 -> upsilon_max}
!                 (upsilon_factor2) e^(-(upsilon) (upsilon_factor1)) d(upsilon) ]
!           = (upsilon_factor2)
!             [ (-upsilon_factor1)^-1 exp(-(upsilon) (upsilon_factor1)) ]{0 -> upsilon_max}
!           = (upsilon_factor2 / upsilon_factor1)
!             [-exp(-(upsilon_max) (upsilon_factor1)) + 1]
!           = (1 - exp(-(upsilon_max) (upsilon_factor1))) (upsilon_factor2) / upsilon_factor1
!
N_upsilon = (1.0d0 - exp(-upsilon_max*upsilon_factor1)) * upsilon_factor2 / upsilon_factor1

expected_upsilon = (upsilon_factor2 / (upsilon_factor1 * N_upsilon)) * &
                   (     -upsilon_max*exp(-upsilon_max * upsilon_factor1) + &
                         (1.0d0 - exp(-upsilon_max * upsilon_factor1)) / upsilon_factor1    )

!
!Obtain the expectation value of the vibrational energy given
!the expecation value of the vibrational quantum number
!
! With ZPE:
!
! <Evib> = [ integral{0 -> upsilon_max}
!                 (h) (vib_frequency) (upsilon + 1/2) P(upsilon) d(upsilon) ] / N
!
!        = (h) (vib_frequency) [ integral{0 -> upsilon_max}
!                 (upsilon) P(upsilon) d(upsilon) ] / N
!           + (h) (vib_frequency) [ integral{0 -> upsilon_max}
!                       (1/2) P(upsilon) d(upsilon) ] / N
!        = (h) (vib_frequency) (<upsilon>) + (1/2) (h) (vib_frequency)
!        = (h) (vib_frequency) (<upsilon> + 1/2)
!
expected_evib = epsilon_factor1*(expected_upsilon+0.5d0)

! 
! Without ZPE:
!
! <Evib> = [ integral{0 -> upsilon_max}
!                 (h) (vib_frequency) (upsilon ) P(upsilon) d(upsilon) ] / N
!
!        = (h) (vib_frequency) [ integral{0 -> upsilon_max}
!                 (upsilon) P(upsilon) d(upsilon) ] / N
!        = (h) (vib_frequency) (<upsilon>)
!
!expected_evib = epsilon_factor1*expected_upsilon

!
!Obtain the expectation value of the temperature given
!the expecation value of the vibrational energy
!
!From before:
!
! With ZPE:
!
!   T = theta_vib / log( 1 + (2) (kB) (theta_vib) /
!                            ((2) (average_Evib) - (kB) (theta_vib)))
!
! <T> = theta_vib / log( 1 + (2) (kB) (theta_vib) /
!                            ((2) (<Evib>) - (kB) (theta_vib)))
!     = theta_vib / log( 1 + (2) (kB) (theta_vib) /
!                            ((2) (h) (vib_frequency) (<upsilon> + 1/2) - (kB) (theta_vib)))
!     = theta_vib / log( 1 + (2) (kB) (theta_vib) /
!                            ((2) (kB) (theta_vib) (<upsilon> + 1/2) - (kB) (theta_vib)))
!     = theta_vib / log( 1 + (2) (kB) (theta_vib) /
!                            (2) (kB) (theta_vib) (<upsilon>))
!     = theta_vib / log( 1 + 1 / <upsilon>)
!
expected_temperature = theta_vib / log(1.0d0 + 1.0d0 / expected_upsilon)

!
! Without ZPE:
!
!   T = theta_vib / log( 1 + (kB) (theta_vib) / average_Evib)
!
! <T> = theta_vib / log( 1 + (kB) (theta_vib) / <Evib>)
!     = theta_vib / log( 1 + (kB) (theta_vib) / (h) (vib_frequency) (<upsilon>))
!     = theta_vib / log( 1 + (kB) (theta_vib) / (kB) (theta_vib) (<upsilon>))
!     = theta_vib / log( 1 + 1 / <upsilon>)
!
!expected_temperature = theta_vib / log(1.0d0 + 2.0d0 / (2.0d0*expected_upsilon + 1.0d0))

!
! BUG TESTING START
!print *, ""
!print *, "expected upsilon: ", expected_upsilon
!print *, "expected evib: ", expected_evib
!print *, "expected temperature: ", expected_temperature
!print *, ""
! BUG TESTING END
!

!
!For later, also obtain the normalization constant of the corresponding
!integral but for the rotational quantum number (J)
!
!         N = integral{0 -> J_max}
!                           P(J) d(J)
!
!           = [ integral{0 -> J_max}
!                 ((2) (J) + 1) (J_factor1) e^(-(J_factor1) (J) (J + 1)) d(J) ]
!           = (J_factor1)
!             [ (-J_factor1)^-1 exp(-(J_factor1) (J) (J + 1)) ]{0 -> J_max}
!           = [-exp(-(J_factor1) (J_max) (J_max + 1)) + 1]
!
!           = 1 - exp(-(J_factor1) (J_max) (J_max + 1))
!
!Because J_factor1 depends on the bond length, use the average
!bond length
!
! J_factor1 = theta_rot / temperature
!           = hbar^2 / (moment_inertia) (kB) (temperature)
!           = hbar^2 / ((1/2) (mass_hydrogen) (bond_length)^2) (kB) (temperature)
!           = (2) (hbar^2) / (mass_hydrogen) (average_r0^2) (kB) (temperature)
!
!Thus:
!
! N = 1 - exp(-((2) (hbar^2) / (mass_hydrogen) (average_r0^2) (kB) (temperature))
!              (J_max) (J_max + 1))
!
N_J = (1.0d0 - exp(-(2.0d0 * ((hbar/(RU_energy * RU_time))**2) / &
                     (mass_hydrogen * (average_r0**2) * (kb / RU_energy) * temperature)) * &
                    (J_max) * (J_max + 1.0d0)))

!
!The average linear rotational speed is displayed on screen
!
write(Ntraj_text,FMT="(F6.3)") average_rot
write(gnuplotchannel,*) 'set label 3 "Average Rotational Speed: '//Ntraj_text//' A/fs" at screen 0.85,0.925'
write(gnuplotchannel,*) 'set label 3 front'

!
!The (rotational) temperature is displayed on screen
!
write(Ntraj_text,FMT="(F6.1)") (2.0d0/2.0d0)*(RU_energy/kb)*average_Erot
write(gnuplotchannel,*) 'set label 4 "Temperature: '//Ntraj_text//' K" at screen 0.85,0.9'
write(gnuplotchannel,*) 'set label 4 front'

!
!There may be multiple bonds in the system so multiple
!distributions may be plotted
!
write(Ntraj_text,FMT="(I6)") Ntraj
write(gnuplotchannel,*) 'set multiplot layout ',Nbonds,',4 '//&
                        'margins 0.025,0.975,.1,.95 spacing 0.05,0 '//&
                        'title "Initial Bond Distribution of '//trim(adjustl(Ntraj_text))//&
                        ' Trajectories" font ",24" offset 0,5'
do i = 1, Nbonds

!
! Bond Orientation (Theta) Distribution
!
! Plotted based on theta angles in the 4th column
!

write(gnuplotchannel,*) 'box_width = 2 * pi /', initialBins
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'bin_number(x) = floor(x/box_width)'
write(gnuplotchannel,*) 'rounded(x) = box_width * (bin_number(x) + 0.5)'

write(gnuplotchannel,*) 'set xlabel "Initial Bond Theta Angle (rad)" font ",18"'
write(gnuplotchannel,*) 'set xrange [-pi:pi]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) "set format x '%.1P π'"
write(gnuplotchannel,*) 'unset xtics'
if (i == Nbonds) write(gnuplotchannel,*) 'set xtics pi/2'

write(Ntraj_text,FMT="(I6)") i*6 - 2
write(gnuplotchannel,*) 'set ylabel "Occurence" font ",18"'
write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//initialfile//&
                        '" u (rounded($'//trim(adjustl(Ntraj_text))//&
			')):(1.0) smooth frequency with boxes'

!
! Bond Orientation (Phi) Distribution
!
! Plotted based on phi angles in the 5th column
!

write(gnuplotchannel,*) 'box_width = pi /', initialBins
write(gnuplotchannel,*) 'set boxwidth box_width'

write(gnuplotchannel,*) 'set xlabel "Initial Bond Phi Angle (rad)" font ",18"'
write(gnuplotchannel,*) 'set xrange [0:pi]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) "set format x '%.1P π'"
write(gnuplotchannel,*) 'unset xtics'
if (i == Nbonds) write(gnuplotchannel,*) 'set xtics pi/2'

write(Ntraj_text,FMT="(I6)") i*6 - 1
write(gnuplotchannel,*) 'set ylabel "Occurence" font ",18"'
write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//initialfile//&
			'" u (rounded($'//trim(adjustl(Ntraj_text))//&
			')):(1.0) smooth frequency with boxes'

!
! Vibrational Quantum Number Distribution
!
! Plotted based on bond distances in the 1st column
!

write(Ntraj_text,FMT="(I6)") i*6 - 5

!
!Obtain the vibrational quantum number of the bond given
!the bond distance when outstretched
!
! Evib = (1/2) (force_constant) (r0 - bond_distance)^2
!
! With ZPE:
!
! Evib = (epsilon_factor1) (upsilon + 1/2)
!
!Rearranging:
!
! upsilon = Evib / epsilon_factor1 - 1/2
!         = (1/2) (force_constant) (r0 - bond_distance)^2 / epsilon_factor1 - 1/2
!
write(gnuplotchannel,*) 'upsilon_factor1 = ', 0.5d0 * HOke_hydrogen / epsilon_factor1
write(gnuplotchannel,*) 'HOr0 = ', HOr0_hydrogen
write(gnuplotchannel,*) 'min_upsilon = upsilon_factor1*(HOr0 -',min_r0, ')**2 - 0.5'
write(gnuplotchannel,*) 'max_upsilon = upsilon_factor1*(HOr0 -',max_r0, ')**2 - 0.5'

!
! Without ZPE:
!
! Evib = (epsilon_factor1) (upsilon)
!
!Rearranging:
!
! upsilon = Evib / epsilon_factor1
!         = (1/2) (force_constant) (r0 - bond_distance)^2 / epsilon_factor1
!
!write(gnuplotchannel,*) 'min_upsilon = upsilon_factor1*(HOr0 -',min_r0, ')**2'
!write(gnuplotchannel,*) 'max_upsilon = upsilon_factor1*(HOr0 -',max_r0, ')**2'

write(gnuplotchannel,*) 'box_width = (max_upsilon-min_upsilon) /', initialBins
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'bin_number(x) = floor((x-min_upsilon)/box_width)'
write(gnuplotchannel,*) 'rounded(x) = min_upsilon+box_width * (bin_number(x) + 0.5)'
write(gnuplotchannel,*) 'set xlabel "Vibrational Quantum Number" font ",18"'
write(gnuplotchannel,*) 'set format x "%5.1e'
write(gnuplotchannel,*) 'unset xtics'
if (i == Nbonds) write(gnuplotchannel,*) 'set xtics min_upsilon, (max_upsilon-min_upsilon)/5, max_upsilon'
write(gnuplotchannel,*) 'set ylabel "Occurence" font ",18"'
write(gnuplotchannel,*) 'set xrange [min_upsilon:max_upsilon]'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'uf1 = ', upsilon_factor1
write(gnuplotchannel,*) 'uf2 = ', upsilon_factor2

!
!Plot the theoretical distribution for a specific number of trajectories
!and distribution box width
!
!The area under the curve should equal the area that all boxes occupy
!
!f(upsilon) = (1/N(upsilon)) * P(upsilon) * Area_Under_Curve(f)
!
!f(upsilon) = (1/N(upsilon)) * P(upsilon) * (Ntrajectories * box_width)
!
!           = (upsilon_factor2) e^(-(upsilon) (upsilon_factor1)) *
!                     (box_width) (Ntrajectories) / N(upsilon)
!
write(gnuplotchannel,*) 'f(x) = uf2 * exp(-x*uf1) * box_width *', Ntraj_max / N_upsilon
write(gnuplotchannel,*) 'set style line 1 linecolor rgb "red" linewidth 1.5'

!
! With ZPE:
!
write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//initialfile//&
			'" u (rounded((upsilon_factor1*(HOr0-$'//trim(adjustl(Ntraj_text))//&
			')**2)-0.5)):(1.0) smooth frequency with boxes,\'
!
! Without ZPE:
!
!write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//initialfile//&
!			'" u (rounded(upsilon_factor1*(HOr0-$'//trim(adjustl(Ntraj_text))//&
!			')**2)):(1.0) smooth frequency with boxes,\'
write(gnuplotchannel,*) '     f(x) with line linestyle 1'

!
! Rotational Quantum Number Distribution
!
! Plotted based on linear rotational speeds in the 2nd column
! as well as bond distances in the 1st column
!

write(variable_length_text,FMT="(I5)") i*6 - 5
write(Ntraj_text,FMT="(I6)") i*6 - 4

!
!Obtain the rotational quantum number of the bond given
!the bond distance when outstretched and the
!linera rotational speed
!
! Erot  = (2) (moment_inertia) (linear_rotational_speed)^2 / bond_length^2
!       = (2) ((1/2) (mass_hydrogen) (bond_length)^2) *
!             (linear_rotational_speed)^2 / bond_length^2
!       = (mass_hydrogen) (linear_rotational_speed)^2
!
! Erot = (epsilon_factor2) (J) (J + 1) / moment_inertia
!      = (epsilon_factor2) (J) (J + 1) / ((1/2) (mass_hydrogen) (bond_length)^2)
!
!Rearranging:
!
! 0 =J^2 + J - (1/2) (Erot) (mass_hydrogen) (bond_length)^2 / epsilon_factor2
!
! J = (-1/2) + (1/2) *
!              sqrt(1 + (2) (Erot) (mass_hydrogen)
!                           (bond_length)^2 / epsilon_factor2)
!   = (-1/2) + (1/2) *
!              sqrt(1 + (2) (mass_hydrogen)^2 (linear_rotational_speed)^2
!                           (bond_length)^2 / epsilon_factor2)
!   = (-1/2) + (1/2) *
!              sqrt(1 + ((2) (mass_hydrogen^2) / epsilon_factor2) *
!                       ((linear_rotational_speed) (bond_length))^2)
!
write(gnuplotchannel,*) 'J_factor1 = ', 2.0d0 * (mass_hydrogen**2) / epsilon_factor2
write(gnuplotchannel,*) 'min_J = -0.5 + 0.5*sqrt(1.0+J_factor1*', (min_rot*min_r0)**2, ')'
write(gnuplotchannel,*) 'max_J = -0.5 + 0.5*sqrt(1.0+J_factor1*', (max_rot*max_r0)**2, ')'
write(gnuplotchannel,*) 'box_width = (max_J-min_J) /', initialBins
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'bin_number(x) = floor((x-min_J)/box_width)'
write(gnuplotchannel,*) 'rounded(x) = min_J+box_width * (bin_number(x) + 0.5)'
write(gnuplotchannel,*) 'set xlabel "Rotational Quantum Number" font ",18"'
write(gnuplotchannel,*) "set format x '%.1f'"
write(gnuplotchannel,*) 'unset xtics'
if (i == Nbonds) write(gnuplotchannel,*) 'set xtics min_J, (max_J-min_J)/5, max_J'
write(gnuplotchannel,*) 'set ylabel "Occurence" font ",18"'
write(gnuplotchannel,*) 'set xrange [min_J:max_J]'
write(gnuplotchannel,*) 'set yrange [0:]'

!
!Plot the theoretical distribution for a specific number of trajectories
!and distribution box width
!
!The area under the curve should equal the area that all boxes occupy
!
!Assume for simplicity that all bonds are the length of the
!average bond length
!
!f(J) = (1/N(J)) * P(J) * Area_Under_Curve(f)
!
!f(J) = (1/N(J)) * P(J) * (Ntrajectories * box_width)
!
!     = ((2) (J) + 1) (J_factor1) e^(-(J_factor1) (J) (J + 1)) *
!                     (box_width) (Ntrajectories) / N(J)
!
write(gnuplotchannel,*) 'Jf1 = ', 2.0d0 * ((hbar/(RU_energy*RU_time))**2) / &
                                  (temperature * mass_hydrogen * (kb / RU_energy) * average_r0**2)
write(gnuplotchannel,*) 'f(x) = Jf1 * (2*x+1) * exp(-Jf1*x*(x+1)) * box_width * ', Ntraj_max / N_J

write(gnuplotchannel,*) 'plot "'//gridpath0//prefix_filename//initialfile//&
			'" u (rounded(-0.5+0.5*sqrt(1.0+J_factor1*(($'//trim(adjustl(Ntraj_text))//&
			')*($'//trim(adjustl(variable_length_text))//'))**2))):(1.0) '//&
                        'smooth frequency with boxes,\'
write(gnuplotchannel,*) '     f(x) with line linestyle 1'
end do

close(gnuplotchannel)

!And then we just input it into gnuplot.exe
call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)


end subroutine getInitialimages



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               getComparedScatteringAngles
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine reads all selected set of trajectories and compares
!               them with respect to some number of trajectories and some SATRV
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               lowerlimit                      REAL(DP)                        The lower bound of the distribution
!               upperlimit                      REAL(DP)                        The upper bound of the distribution
!               imagename                       CHARACTER(*)                    The name of the output image
!               SATRVcolumn                     INTEGER                         The column in the SATRV file this distribution
!                                                                               is recorded in
!               SATRVname                       CHARACTER(*)                    The name of the distribution (needs to be exact)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               Ntesttraj                       INTEGER                         The number of trajectories per sample
!                                                                               set of trajectories, or the "sampling"
!               Nbins                           INTEGER                         The number of bins used for distributions
!               comparison_number               INTEGER                         The number of sets of trajectories being
!                                                                               compared
!
!               binTotal                        INTEGER,DIM(                    The size of each bin of the distribution
!                                               Nbins,comparison_number         for each set of trajectories
!
!               referenceBins                   REAL,DIM(Nbins)                 The x-value corresponding to a bin
!               referenceMeans                  REAL,DIM(Nbins)                 The y-value corresponding to a bin
!               referenceSDs                    REAL,DIM(Nbins)                 The y-value uncertainty corresponding to a bin
!
!               comparisonCDF                   REAL                            The cumulative distribution function of some bin
!               comparisonKS                    REAL                            The Kolmogorov-Smirnov difference of some bin
!               comparisonRMSD                  REAL                            The RMSD of some distribution
!               comparisonKRP                   REAL                            The KRP of some distribution
!
!               lambda_penalty                  REAL                            The value of lambda used to
!                                                                               calculate the KRP
!               minsd_penalty                   REAL                            The minimum standard deviation used to
!                                                                               calculate the KRP
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath0//Adjusted//           DAT                             Stores the averages and deviations for the 
!                 SATRVname                                                     distribution of some SATRV for all the
!                                                                               trajectories in the library for some sampling
!               gridpath0//allprefixes(...)//   DAT                             Same as the SATRV file but the values are
!                 binnedSATRVfile                                               stored as integers corresponding to some
!                                                                               bin in a specific binning
!               gridpath0//imagename            PNG                             Compares distributions as selected in
!                                                                               allprefixes by some sampling and some SATRV
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine getComparedScatteringAngles(lowerlimit,upperlimit,imagename,SATRVcolumn,SATRVname)
use PARAMETERS
use PHYSICS
use ANALYSIS
implicit none

!FORMAT OF DAT FILES HOUSING SCATTERING ANGLES
integer,intent(in) :: SATRVcolumn
character(*), intent(in) :: SATRVname
integer,dimension(5) :: SATRVdata

!FORMAT OF PNG FILES TO BE MADE
character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5
character(*), intent(in) :: imagename

!Upper and Lower Limits for the plot
real,intent(in) :: upperlimit, lowerlimit

!FORMATTING OF PNG FILES
character(5) :: variable_length_text
character(8) :: difference_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(trajectory_text_length) :: Ntraj_text
character(6) :: boxwidth_text
logical :: grid_is_done

integer :: iostate,i,j
real :: speed_out, ScatteringAngle

real,allocatable :: referenceBins(:)
real,allocatable :: referenceMeans(:)
real,allocatable :: referenceSDs(:)
integer,allocatable :: binTotal(:,:)
real :: comparisonRMSD, comparisonKS, comparisonCDF, comparisonKRP
integer :: Nbins

gridpath4 = gridpath0//expfolder
gridpath5 = gridpath4//intermediatefolder

!Initialization
if (SATRVname == "ScatteringAngle") then
        Nbins = scatteringangleBins
else
        Nbins = energychangeBins
end if

allocate(referenceBins(Nbins),&
         referenceMeans(Nbins),&
         referenceSDs(Nbins))
allocate(binTotal(Nbins,comparison_number))

open(filechannel1,file=gridpath5//"Adjusted"//SATRVname//cumulativefile//".dat")
do j = 1, Nbins
        read(filechannel1,FMT=*) referenceBins(j), referenceMeans(j), referenceSDs(j)
end do
close(filechannel1)

binTotal = 0
open(filechannel1,file=gridpath5//allprefixes(1:alllengths(1)-1)//binnedSATRVfile)
do j = 1, Ntesttraj
        read(filechannel1,FMT=*) SATRVdata
        binTotal(SATRVdata(SATRVcolumn),1) = &
                 binTotal(SATRVdata(SATRVcolumn),1) + 1
end do
close(filechannel1)

do i = 1, comparison_number-1
        open(filechannel1,file=gridpath5//allprefixes(sum(alllengths(1:i))+1:sum(alllengths(1:i+1))-1)//&
                               binnedSATRVfile)
        do j = 1, Ntesttraj
                read(filechannel1,FMT=*) SATRVdata
                binTotal(SATRVdata(SATRVcolumn),i+1) = &
                         binTotal(SATRVdata(SATRVcolumn),i+1) + 1
        end do
        close(filechannel1)
end do

!This is the gnuplot code to make the plots
open(gnuplotchannel,file=gridpath5//gnuplotfile)
write(gnuplotchannel,*) 'set term pngcairo enhanced size 1200,1200'
write(gnuplotchannel,*) 'set encoding utf8'
write(gnuplotchannel,*) 'set output "'//gridpath4//imagename//'.png"'
write(gnuplotchannel,*) 'unset key'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
write(gnuplotchannel,*) 'unset xtics'
write(variable_length_text,FMT=FMT5_variable) trajectory_text_length
write(Ntraj_text,FMT="(I"//variable_length_text//")") Ntraj
write(gnuplotchannel,*) 'set multiplot layout ',comparison_number,&
                        ',1 columnsfirst margins 0.1,0.95,.1,.9 spacing 0.1,0 title '//&
                        '"'//SATRVname//' Distribution of '//trim(adjustl(Ntraj_text))//&
                        ' Trajectories" font ",32" offset 0,-3'
write(gnuplotchannel,*) 'set style histogram clustered gap 1'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
write(gnuplotchannel,*) 'set ylabel "Occurence" font ",24"'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set ytics font ",16"'

if (SATRVname == "ScatteringAngle") then
        write(gnuplotchannel,*) 'scaling = 1'
        write(gnuplotchannel,*) 'pi = 3.14159265'
else
        write(gnuplotchannel,*) 'scaling = 1000'
end if

write(gnuplotchannel,*) 'lowerlimit = scaling * ', lowerlimit
write(gnuplotchannel,*) 'upperlimit = scaling * ', upperlimit
write(gnuplotchannel,*) 'set xrange [lowerlimit:upperlimit]'
write(gnuplotchannel,*) 'box_width = (upperlimit - lowerlimit) / ', Nbins
write(gnuplotchannel,*) 'set boxwidth box_width'
write(gnuplotchannel,*) 'bin_number(x) = floor(scaling * x/box_width)'
write(gnuplotchannel,*) 'rounded(x) = box_width * (x - 0.5)'

comparisonRMSD = 0
comparisonKS = 0
comparisonCDF = 0
comparisonKRP = 0
do j = 1, Nbins
        comparisonCDF = comparisonCDF + 1.0*binTotal(j,1) - referenceMeans(j)
        comparisonKS = max(comparisonKS, abs(comparisonCDF))
        comparisonRMSD = comparisonRMSD + (1.0*bintotal(j,1) - referenceMeans(j))**2
        comparisonKRP = comparisonKRP + (abs(1.0*bintotal(j,1) - referenceMeans(j)) / &
                                        max(minsd_penalty,referenceSDs(j)))**lambda_penalty
end do

write(gnuplotchannel,*) 'set label 1 "'//allprefixes(1:alllengths(1)-1)//'" at graph 0.825, 0.9'
write(difference_text,FMT="(F8.4)") sqrt(comparisonRMSD * 1.0 / (Nbins - 1.0))
write(gnuplotchannel,*) 'set label 2" RMSD: '//difference_text//'" at graph 0.85,0.825'
write(difference_text,FMT="(F8.4)") (comparisonKRP * 1.0 / (Nbins - 1.0))**(1.0/lambda_penalty)
write(gnuplotchannel,*) 'set label 3" KRP: '//difference_text//'" at graph 0.85,0.750'
write(variable_length_text,FMT=FMT5_variable) SATRVcolumn
!write(gnuplotchannel,*) 'plot "'//gridpath0//allprefixes(1:alllengths(1))//&
!                        SATRVfile//'" u (rounded($'//trim(adjustl(variable_length_text))//&
!                        ')):(1.0) smooth frequency w boxes, \'
write(gnuplotchannel,*) 'plot "'//gridpath5//allprefixes(1:alllengths(1)-1)//&
                        binnedSATRVfile//'" u (rounded($'//trim(adjustl(variable_length_text))//&
                        ')):(1.0) smooth frequency w boxes, \'
write(gnuplotchannel,*) '     "'//gridpath5//'Adjusted'//SATRVname//cumulativefile//&
                        '.dat" u (scaling*($1)):2 w boxes'//&
                        ' fs transparent solid 0.5 noborder, \'
write(gnuplotchannel,*) '     "'//gridpath5//'Adjusted'//SATRVname//cumulativefile//&
                        '.dat" u (scaling*($1)):2:3 w yerrorbars'
do i = 1, comparison_number-1
        comparisonRMSD = 0
        comparisonKS = 0
        comparisonCDF = 0
        comparisonKRP = 0
        do j = 1, Nbins
                comparisonCDF = comparisonCDF + 1.0*binTotal(j,i+1) - referenceMeans(j)
                comparisonKS = max(comparisonKS, abs(comparisonCDF))
                comparisonRMSD = comparisonRMSD + (1.0*binTotal(j,i+1) - referenceMeans(j))**2
                comparisonKRP = comparisonKRP + (abs(1.0*bintotal(j,i+1) - referenceMeans(j)) / &
                                                max(minsd_penalty,referenceSDs(j)))**lambda_penalty
        end do

        write(gnuplotchannel,*) 'set label ', 3*i+1, '"'//allprefixes(sum(alllengths(1:i))+1:sum(alllengths(1:i+1))-1)//&
                                '" at graph 0.825, 0.9'
        write(difference_text,FMT="(F8.4)") sqrt(comparisonRMSD * 1.0 / (Nbins - 1.0))
        write(gnuplotchannel,*) 'set label ',3*i+2,'" RMSD: '//difference_text//'" at graph 0.85,0.825'
        write(difference_text,FMT="(F8.2)") (comparisonKRP * 1.0 / (Nbins - 1.0))**(1.0/lambda_penalty)
        write(gnuplotchannel,*) 'set label ',3*i+3,'" KRP: '//difference_text//'" at graph 0.85,0.750'
        write(gnuplotchannel,*) 'unset label ', 3*i-2
        write(gnuplotchannel,*) 'unset label ', 3*i-1
        write(gnuplotchannel,*) 'unset label ', 3*i
        if (i == comparison_number-1) then
                if (SATRVname == "ScatteringAngle") then
                        write(gnuplotchannel,*) 'set xlabel "Scattering Angle (rad)" font ",24"'
                        write(gnuplotchannel,*) 'set xtics lowerlimit, 10*(upperlimit-lowerlimit)/',&
                                                Nbins,', upperlimit font ",16"'
                        write(gnuplotchannel,*) "set format x '%.3P π'"
                else
                        write(gnuplotchannel,*) 'set xlabel "Energy Change (meV)" font ",24"'
                        write(gnuplotchannel,*) 'set xtics lowerlimit, 10*box_width, upperlimit font ",16"'
                        write(gnuplotchannel,*) "set format x '%.3f'"
                end if
        end if
        write(variable_length_text,FMT=FMT5_variable) SATRVcolumn
!        write(gnuplotchannel,*) 'plot "'//gridpath0//allprefixes(sum(alllengths(1:i))+1:sum(alllengths(1:i+1)))//&
!                                SATRVfile//'" u (rounded($'//trim(adjustl(variable_length_text))//&
!                                ')):(1.0) smooth frequency w boxes, \'
        write(gnuplotchannel,*) 'plot "'//gridpath5//allprefixes(sum(alllengths(1:i))+1:sum(alllengths(1:i+1))-1)//&
                                binnedSATRVfile//'" u (rounded($'//trim(adjustl(variable_length_text))//&
                                ')):(1.0) smooth frequency w boxes, \'
        write(gnuplotchannel,*) '     "'//gridpath5//'Adjusted'//SATRVname//cumulativefile//&
                                '.dat" u (scaling*($1)):2 w boxes'//&
                                ' fs transparent solid 0.5 noborder, \'
        write(gnuplotchannel,*) '     "'//gridpath5//'Adjusted'//SATRVname//cumulativefile//&
                                '.dat" u (scaling*($1)):2:3 w yerrorbars'
end do
close(gnuplotchannel)

deallocate(referenceBins,referenceMeans,referenceSDs)
deallocate(binTotal)

!And then we just input it into gnuplot.exe
call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

end subroutine getComparedScatteringAngles





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               postProcess
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine reads all initial and final frames for a set of trajectories specified
!               by the prefix_filename, then processes them into a few collective variables per
!               trajectory like scattering angle and energy change decompositions; these values are
!               then written onto a final file
!
!               While reading the files, it also searches for and records the minimums and maximums
!               found for each variable evaluated
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               prefix_filename                 CHARACTER(*)                    The prefix defining a set of trajectories
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               Ntraj                           INTEGER                         The number of trajectories to be read
!                                                                               for this set of trajectories
!               coords_initial                  REAL,DIM(3,Natoms)              The coordinates of the initial frame
!               velocities_initial              REAL,DIM(3,Natoms)              The velocities of the initial frame
!               coords_final                    REAL,DIM(3,Natoms)              The coordinates of the final frame
!               velocities_final                REAL,DIM(3,Natoms)              The velocities of the final frame
!
!               velocity_in                     REAL,DIM(3)                     The velocity of the incoming molecule
!               velocity_out                    REAL,DIM(3)                     The velocity of the outgoing molecule
!               speed_in                        REAL                            The speed of the incoming molecule
!               speed_out                       REAL                            The speed of the outgoing molecule
!               scatteringAngle                 REAL                            The angle between velocity_in and
!                                                                               velocity_out
!
!               velocity1_CM                    REAL,DIM(3)                     The velocity of the center of mass of
!                                                                               all atoms of the initial frame
!               velocity2_CM                    REAL,DIM(3)                     The velocity of the center of mass of
!                                                                               all atoms of the final frame
!               TotalEnergy1                    REAL                            Total kinetic energy of the initial frame
!               TotalEnergy2                    REAL                            Total kinetic energy of the final frame
!               RotationalEnergy1               REAL                            Rotational energy of the initial frame
!               RotationalEnergy2               REAL                            Rotational energy of the final frame
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath0//prefix_filename//    DAT                             Stores the initial and final frames of a
!                 timeslicefile                                                 set of trajectories corresponding to the
!                                                                               specified prefix
!               gridpath0//prefix_filename//    DAT                             Unused at the moment for this subroutine
!                 intialfile                                            
!               gridpath0//prefix_filename//    DAT                             Stores the scattering angle and energy change
!                 SATRVfile                                                     decomposition for a set of trajectories
!                                                                               corresponding to the specified prefix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine postProcess(prefix_filename)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
use PHYSICS
implicit none

!FILENAMING
character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5
character(*),intent(in) :: prefix_filename

!Translational, Rotational, Vibrational Energies
real(dp) :: scatteringAngle, speed_in, speed_out
real(dp) :: bond_distance, rot_speed
real(dp) :: relTranslationalEnergy1, relTranslationalEnergy2
real(dp) :: absTranslationalEnergy1, absTranslationalEnergy2
real(dp) :: RotationalEnergy1, RotationalEnergy2
real(dp) :: TotalEnergy1, TotalEnergy2

!TRAJECTORY DATA
!real(dp),dimension(Nbonds,6) :: initial_bonding_data
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
!		read(filechannel2,FMTinitial) ((initial_bonding_data(i,j),j=1,6),i=1,Nbonds)
!		read(filechannel2,FMT=*) ((initial_bonding_data(j,i),j=1,6),i=1,Nbonds)

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
