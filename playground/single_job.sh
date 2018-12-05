#!/bin/bash

echo ""
echo ""
echo "Now running the background job"
echo ""

###############################################################################################################################################
###############################################################################################################################################

#Enter the name of the original source code files

#The parameters file
oldPARAMETERS="PARAMETERS"

#The analysis file
oldANALYSIS="ANALYSIS"

#The physics file
oldPHYSICS="PHYSICS"

#The variables file
oldVARIABLES="VARIABLES"

#The makefile for grid creation
oldMAKEGRID="make_makeGridwithNewTrajectories2"

#The makefile for grid analysis
oldMAKEANALYSIS="make_checkNewTrajectorieswithMultipleGrids2"

###############################################################################################################################################
###############################################################################################################################################

#Enter the name of the source code files to be made in the new directory
#These must be different from the original source code file names

#The parameters file
newPARAMETERS=PARAMETERS_new

#The analysis file
newANALYSIS=ANALYSIS_new

#The makefile for grid creation
newMAKEGRID=make_makeGridwithNewTrajectories2_new

#The makefile for grid analysis
newMAKEANALYSIS=make_checkNewTrajectorieswithMultipleGrids2_new

###############################################################################################################################################
###############################################################################################################################################

#Change the below variables for whatever grid you want
#If you want to change other variables (ex. spacing1) you would
#have to make more sed statements below to change that

#How many children generations do you want?
Norder_max=1

#The ratio of child-level to parent-level cell spacings with respect to var1
scaling1_0=004
scaling1_1=004

#The ratio of child-level to parent-level cell spacings with respect to var2
scaling2_0=004
scaling2_1=004

#The number of frames a parent-level cell accepts before it is subdivided
overcrowd0=00050
overcrowd1=10000
overcrowd2=01010

#The number of trajectories simulated and added to a new grid
Ntraj_max=0700

#The number of grids to add to a new library
Ngrid_max=12

#Whether to add duplicate copies or use a labelling scheme
force_Duplicates=.false.

#Whether to use the labelling scheme or drop it altogether
force_NoLabels=.false.

#The default flags to be used for analyses
#Of course, you don't want all analyses to be the same so go down to each analysis and change
#what you want each individual one to do
heatmap_flag=.false.
trueSA_flag=.false.
trueED_flag=.false.
testtraj_flag=.true.
testtrajRMSD_flag=.false.
percentthreshold_flag=.false.
testtrajSA_flag=.false.
testtrajSAheatmap_flag=.false.
#threshold_rmsd=.200100d0
threshold_rmsd=.222000d0
threshold_rmsd1=.222000d0
threshold_rmsd2=.100000d0
threshold_rmsd3=.100000d0
threshold_rmsd4=.050000d0
threshold_rmsd5=.050000d0
reject_flag=.false.
accept_first=.false.
accept_worst=.false.
Ngrid_cap=1
Norder_cap=1
#Ngrid_cap=${Ngrid_max}
Ntrajectories=2

#These are flags relating to using old data
useolddata_flag=.false.
useoldinitialbonddata_flag=.false.
initialbondname="001omegaA.20000"

#If you have special set of parameters you want to compare, list them here
#These will be compared at each compilation

#These are the options for comparison:
#ScatteringAngle (rad)
#AbsoluteEnergyChange (eV)
#RelativeEnergyChange (eV)
#RotationalEnergyChange (eV)

#If the comparison lower and upper limits are the same, the program will
#use whatever the minimum and maximum is of the data (bad if outliers exist)
comparison_flag=none
comparison_lowerlimit="0.0d0"
comparison_upperlimit="0.200d0"

declare -a prefixes
#prefixes[0]="001accept.15000"
#prefixes[1]="002accept.15000"
#prefixes[2]="004accept.15000"
#prefixes[3]="008accept.15000"
#prefixes[4]="001reject.15000"
#prefixes[0]="004accept.05000"
#prefixes[1]="004accept.01000"
#prefixes[2]="004accept.00500"
#prefixes[0]="001omegaA.20000"
#prefixes[1]="001omegaA.15000"
#prefixes[2]="001omegaA.10000"
#prefixes[3]="001omegaA.05000"
prefixes[0]="001reject.05000"

###############################################################################################################################################
###############################################################################################################################################

#Enter what path the original source code is on and what you would like
#the library housing all the grids (a folder) to be called

#The name of the new library (folder)
#newGRID=HH_${scaling1_0}_${scaling2_0}_${overcrowd0}_${Ntraj_max}_1
newGRID="HH2_Oct16_label"

#If you want to make a new grid, set this to 1; otherwise, set it to zero
newGRID_flag=0

#How often you want to check the progress of the new grid's creation
#(has an intrinsic minimum of Ntraj_max/10)
#Just set this to a very large number if no progress checks are wanted
newGRID_check_min=100

#The number of post-grid analyses you would like done
#These are separate from the comparison and the post-grid-making analysis
Nanalyses=1

#The path that has the original source code
currentPATH=$(pwd)

#The file that keeps the times of the code
bashout="timing.dat"

#DO NOT TOUCH THESE
gridPATH=$currentPATH/$newGRID
newSOURCE=SOURCE
newPATH=$(pwd)/$newGRID/$newSOURCE

###############################################################################################################################################
###############################################################################################################################################
#		GRID CREATION AND SUMMARY ANALYSIS
###############################################################################################################################################
###############################################################################################################################################

#If you want to create a new grid, do that first
if [ $newGRID_flag -eq 1 ]
then

#If there is another folder of the same name delete that folder first
rm -r $currentPATH/$newGRID

#Make the directories, copy all the original source code, etc.
mkdir $currentPATH/$newGRID
mkdir $newPATH/
cp $currentPATH/*.f90 $newPATH/
cp $currentPATH/make_$(echo "*") $newPATH/

#Make changes to the parameters file as specified in the variables above
#Unless you want to change MORE variables, don't touch this
sed "s|Ntraj_max = [0-9]*|Ntraj_max = $Ntraj_max| 
     s|overcrowd0 = [0-9]*|overcrowd0 = $overcrowd0|
     s|overcrowd1 = [0-9]*|overcrowd1 = $overcrowd1|
     s|overcrowd2 = [0-9]*|overcrowd2 = $overcrowd2|
     s|scaling1_0 = [0-9]*|scaling1_0 = $scaling1_0|
     s|scaling2_0 = [0-9]*|scaling2_0 = $scaling2_0|
     s|scaling1_1 = [0-9]*|scaling1_1 = $scaling1_1|
     s|scaling2_1 = [0-9]*|scaling2_1 = $scaling2_1|
     s|Norder_max = [0-9]*|Norder_max = $Norder_max|
     s|Ngrid_max = [0-9]*|Ngrid_max = $Ngrid_max|
     s|trajectory_text_length = [0-9]*|trajectory_text_length = ${#Ntraj_max}|
     s|scaling1_text_length = [0-9]*|scaling1_text_length = ${#scaling1_0}|
     s|scaling2_text_length = [0-9]*|scaling2_text_length = ${#scaling2_0}|
     s|overcrowd0_text_length = [0-9]*|overcrowd0_text_length = ${#overcrowd0}|
     s|gridpath_length = .*|gridpath_length = $((${#gridPATH}+1))|
     s|Ngrid_check_min = .*|Ngrid_check_min = $newGRID_check_min|
     s|force_Duplicates = .*|force_Duplicates = $force_Duplicates|
     s|force_NoLabels = .*|force_NoLabels = $force_NoLabels|
     s|gridpath0 = .*|gridpath0 = \\&\\n\"$gridPATH/\"|
     s|$oldPARAMETERS\\.f90|$newPARAMETERS.f90|" <$currentPATH/$oldPARAMETERS.f90 >$newPATH/$newPARAMETERS.f90

#Make changes to the analysis file as specified in the variables above
#This first post-grid-making analysis does nothing but look at the data
#in the grid that was just made (heatmap, scattering angle, etc)
sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_cap/
     s|Norder_cap = [0-9]*|Norder_cap = $Norder_cap|
     s/heatmap_flag = .*/heatmap_flag = .true./
     s/trueSA_flag = .*/trueSA_flag = .true./
     s/trueED_flag = .*/trueED_flag = .true./
     s/testtraj_flag = .*/testtraj_flag = .false./
     s/useolddata_flag = .*/useolddata_flag = $useolddata_flag/
     s/useoldinitialbonddata_flag = .*/useoldinitialbonddata_flag = $useoldinitialbonddata_flag/
     s|initialbondname_length = .*|initialbondname_length = $((${#initialbondname}))|
     s/initialbondname = .*/initialbondname = \"$initialbondname\"/
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/comparison_flag = .*/comparison_flag = .false./
     s/percentthreshold_flag = .*/percentthreshold_flag = .false./
     s/threshold_rmsd = .*/threshold_rmsd = $threshold_rmsd/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s/accept_first = .*/accept_first = $accept_first/
     s/accept_worst = .*/accept_worst = $accept_worst/
     s/testheatmapSA_flag = .*/testheatmapSA_flag = .false./
     s/testtrajSA_flag = .*/testtrajSA_flag = .false./" <$currentPATH/$oldANALYSIS.f90 >$newPATH/$newANALYSIS.f90

#DO NOT TOUCH THIS
sed "s/$oldPARAMETERS\\.o/$newPARAMETERS.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEGRID >$newPATH/$newMAKEGRID

#DO NOT TOUCH THIS
sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS



#Now, all we need to do is go into the new folder and make the output files
cd $newPATH

make -f $newPATH/$newMAKEGRID

echo $(date) >> $newPATH/$bashout
/usr/bin/time -a -o $newPATH/$bashout -f "GRIDMAKING %E  %U  %S  %P  %O" ./a.out

make -f $newPATH/$newMAKEANALYSIS

echo $(date) >> $newPATH/$bashout
/usr/bin/time -a -o $newPATH/$bashout -f "ANALYSIS00 %E  %U  %S  %P  %O" ./a.out

fi

###############################################################################################################################################
###############################################################################################################################################


###############################################################################################################################################
###############################################################################################################################################
#		COMPARISON (USUALLY DONE AFTER ANALYSES)
###############################################################################################################################################
###############################################################################################################################################

#This prevents a new parameters file from being copied into the library
#While at the same time, updating all of the .f90 files
#Remember, the parameters file should never change after creation!
shopt -s extglob
cp $currentPATH/!($oldPARAMETERS|$newPARAMETERS|$oldPHYSICS|$oldVARIABLES)+(.f90) $newPATH/
cp $currentPATH/make_$(echo "*") $newPATH/
shopt -s extglob

#Whatever trajectories we want to check, we need to accumulate and pass
#that to the fortran code somehow
allprefixes=""
alllengths=()
for prefix in "${prefixes[@]}"
do
      allprefixes=$allprefixes$prefix
      alllengths+=("${#prefix}")
done
alllengths_statement=$(IFS=, ; echo "${alllengths[*]}")

if [ $comparison_flag == "ScatteringAngle" ] || [ $comparison_flag == "AbsoluteEnergyChange" ] || [ $comparison_flag == "RelativeEnergyChange" ] || [ $comparison_flag == "RotationalEnergyChange" ]
then

#Change the comparison analysis as specified in the variables above
sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_cap/
     s/heatmap_flag = .*/heatmap_flag = .false./
     s/trueSA_flag = .*/trueSA_flag = .false./
     s/trueED_flag = .*/trueED_flag = .false./
     s/testtraj_flag = .*/testtraj_flag = .false./
     s/useolddata_flag = .*/useolddata_flag = .false./
     s/useoldinitialbonddata_flag = .*/useoldinitialbonddata_flag = $useoldinitialbonddata_flag/
     s|initialbondname_length = .*|initialbondname_length = $((${#initialbondname}))|
     s/initialbondname = .*/initialbondname = \"$initialbondname\"/
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = .false./
     s/percentthreshold_flag = .*/percentthreshold_flag = .false./
     s/threshold_rmsd = .*/threshold_rmsd = $threshold_rmsd1/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s/accept_first = .*/accept_first = $accept_first/
     s/accept_worst = .*/accept_worst = $accept_worst/
     s/comparison_flag = .*/comparison_flag = .true./
     s/comparison_lowerlimit = .*/comparison_lowerlimit = $comparison_lowerlimit/
     s/comparison_upperlimit = .*/comparison_upperlimit = $comparison_upperlimit/
     s/comparison_SATRVname = .*/comparison_SATRVname = \"$comparison_flag\"/
     s/comparison_number = .*/comparison_number = ${#prefixes[@]}/
     s/character(11),parameter :: allprefixes = .*/character(${#allprefixes}),parameter :: allprefixes = \"$allprefixes\"/
     s|alllengths = .*|alllengths = (/$alllengths_statement/)|
     s/testheatmapSA_flag = .*/testheatmapSA_flag = .false./
     s/testtrajSA_flag = .*/testtrajSA_flag = .false./" <$currentPATH/$oldANALYSIS.f90 >$newPATH/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

cd $newPATH

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS

echo $(date) >> $newPATH/$bashout
/usr/bin/time -a -o $newPATH/$bashout -f "COMPARISON %E  %U  %S  %P  %O" ./a.out

fi

###############################################################################################################################################
###############################################################################################################################################

if [ $Nanalyses -lt 1 ]
then
exit
fi

###############################################################################################################################################
###############################################################################################################################################
#		FIRST ANALYSIS
###############################################################################################################################################
###############################################################################################################################################

#This prevents a new parameters file from being copied into the library
#While at the same time, updating all of the .f90 files
#Remember, the parameters file should never change after creation!
shopt -s extglob
cp $currentPATH/!($oldPARAMETERS|$newPARAMETERS|$oldPHYSICS|$oldVARIABLES)+(.f90) $newPATH/
cp $currentPATH/make_$(echo "*") $newPATH/
shopt -s extglob

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_cap/
     s/heatmap_flag = .*/heatmap_flag = $heatmap_flag/
     s/trueSA_flag = .*/trueSA_flag = $trueSA_flag/
     s/trueED_flag = .*/trueED_flag = $trueED_flag/
     s/testtraj_flag = .*/testtraj_flag = $testtraj_flag/
     s/useolddata_flag = .*/useolddata_flag = $useolddata_flag/
     s/useoldinitialbonddata_flag = .*/useoldinitialbonddata_flag = $useoldinitialbonddata_flag/
     s|initialbondname_length = .*|initialbondname_length = $((${#initialbondname}))|
     s/initialbondname = .*/initialbondname = \"$initialbondname\"/
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/percentthreshold_flag = .*/percentthreshold_flag = $percentthreshold_flag/
     s/threshold_rmsd = .*/threshold_rmsd = $threshold_rmsd1/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s/accept_first = .*/accept_first = $accept_first/
     s/accept_worst = .*/accept_worst = $accept_worst/
     s/comparison_flag = .*/comparison_flag = .false./
     s/testheatmapSA_flag = .*/testheatmapSA_flag = $testtrajSAheatmap_flag/
     s/testtrajSA_flag = .*/testtrajSA_flag = $testtrajSA_flag/" <$currentPATH/$oldANALYSIS.f90 >$newPATH/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

cd $newPATH

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS

echo $(date) >> $newPATH/$bashout
/usr/bin/time -a -o $newPATH/$bashout -f "ANALYSIS01 %E  %U  %S  %P  %O" ./a.out

###############################################################################################################################################
###############################################################################################################################################

if [ $Nanalyses -lt 2 ]
then
exit
fi

###############################################################################################################################################
###############################################################################################################################################
#		SECOND ANALYSIS
###############################################################################################################################################
###############################################################################################################################################

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_cap/
     s/heatmap_flag = .*/heatmap_flag = $heatmap_flag/
     s/trueSA_flag = .*/trueSA_flag = $trueSA_flag/
     s/trueED_flag = .*/trueED_flag = $trueED_flag/
     s/testtraj_flag = .*/testtraj_flag = $testtraj_flag/
     s/useolddata_flag = .*/useolddata_flag = $useolddata_flag/
     s/useoldinitialbonddata_flag = .*/useoldinitialbonddata_flag = $useoldinitialbonddata_flag/
     s|initialbondname_length = .*|initialbondname_length = $((${#initialbondname}))|
     s/initialbondname = .*/initialbondname = \"$initialbondname\"/
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/percentthreshold_flag = .*/percentthreshold_flag = $percentthreshold_flag/
     s/threshold_rmsd = .*/threshold_rmsd = $threshold_rmsd2/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s/accept_first = .*/accept_first = $accept_first/
     s/accept_worst = .*/accept_worst = $accept_worst/
     s/comparison_flag = .*/comparison_flag = .false./
     s/testheatmapSA_flag = .*/testheatmapSA_flag = $testtrajSAheatmap_flag/
     s/testtrajSA_flag = .*/testtrajSA_flag = $testtrajSA_flag/" <$currentPATH/$oldANALYSIS.f90 >$newPATH/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS

echo $(date) >> $newPATH/$bashout
/usr/bin/time -a -o $newPATH/$bashout -f "ANALYSIS02 %E  %U  %S  %P  %O" ./a.out

###############################################################################################################################################
###############################################################################################################################################

if [ $Nanalyses -lt 3 ]
then
exit
fi

###############################################################################################################################################
###############################################################################################################################################
#		THIRD ANALYSIS
###############################################################################################################################################
###############################################################################################################################################

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_cap/
     s/heatmap_flag = .*/heatmap_flag = $heatmap_flag/
     s/trueSA_flag = .*/trueSA_flag = $trueSA_flag/
     s/trueED_flag = .*/trueED_flag = $trueED_flag/
     s/testtraj_flag = .*/testtraj_flag = $testtraj_flag/
     s/useolddata_flag = .*/useolddata_flag = $useolddata_flag/
     s/useoldinitialbonddata_flag = .*/useoldinitialbonddata_flag = $useoldinitialbonddata_flag/
     s|initialbondname_length = .*|initialbondname_length = $((${#initialbondname}))|
     s/initialbondname = .*/initialbondname = \"$initialbondname\"/
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/percentthreshold_flag = .*/percentthreshold_flag = $percentthreshold_flag/
     s/threshold_rmsd = .*/threshold_rmsd = $threshold_rmsd3/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s/accept_first = .*/accept_first = $accept_first/
     s/accept_worst = .*/accept_worst = $accept_worst/
     s/comparison_flag = .*/comparison_flag = .false./
     s/testheatmapSA_flag = .*/testheatmapSA_flag = $testtrajSAheatmap_flag/
     s/testtrajSA_flag = .*/testtrajSA_flag = $testtrajSA_flag/" <$currentPATH/$oldANALYSIS.f90 >$newPATH/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS

echo $(date) >> $newPATH/$bashout
/usr/bin/time -a -o $newPATH/$bashout -f "ANALYSIS03 %E  %U  %S  %P  %O" ./a.out

###############################################################################################################################################
###############################################################################################################################################

if [ $Nanalyses -lt 4 ]
then
exit
fi

###############################################################################################################################################
###############################################################################################################################################
#		FOURTH ANALYSIS
###############################################################################################################################################
###############################################################################################################################################

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_cap/
     s/heatmap_flag = .*/heatmap_flag = $heatmap_flag/
     s/trueSA_flag = .*/trueSA_flag = $trueSA_flag/
     s/trueED_flag = .*/trueED_flag = $trueED_flag/
     s/testtraj_flag = .*/testtraj_flag = $testtraj_flag/
     s/useolddata_flag = .*/useolddata_flag = $useolddata_flag/
     s/useoldinitialbonddata_flag = .*/useoldinitialbonddata_flag = $useoldinitialbonddata_flag/
     s|initialbondname_length = .*|initialbondname_length = $((${#initialbondname}))|
     s/initialbondname = .*/initialbondname = \"$initialbondname\"/
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/percentthreshold_flag = .*/percentthreshold_flag = $percentthreshold_flag/
     s/threshold_rmsd = .*/threshold_rmsd = $threshold_rmsd4/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s/accept_first = .*/accept_first = $accept_first/
     s/accept_worst = .*/accept_worst = $accept_worst/
     s/comparison_flag = .*/comparison_flag = .false./
     s/testheatmapSA_flag = .*/testheatmapSA_flag = $testtrajSAheatmap_flag/
     s/testtrajSA_flag = .*/testtrajSA_flag = $testtrajSA_flag/" <$currentPATH/$oldANALYSIS.f90 >$newPATH/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS

echo $(date) >> $newPATH/$bashout
/usr/bin/time -a -o $newPATH/$bashout -f "ANALYSIS04 %E  %U  %S  %P  %O" ./a.out

###############################################################################################################################################
###############################################################################################################################################

if [ $Nanalyses -lt 5 ]
then
exit
fi

###############################################################################################################################################
###############################################################################################################################################
#		FIFTH ANALYSIS
###############################################################################################################################################
###############################################################################################################################################

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_cap/
     s/heatmap_flag = .*/heatmap_flag = $heatmap_flag/
     s/trueSA_flag = .*/trueSA_flag = $trueSA_flag/
     s/trueED_flag = .*/trueED_flag = $trueED_flag/
     s/testtraj_flag = .*/testtraj_flag = $testtraj_flag/
     s/useolddata_flag = .*/useolddata_flag = $useolddata_flag/
     s/useoldinitialbonddata_flag = .*/useoldinitialbonddata_flag = $useoldinitialbonddata_flag/
     s|initialbondname_length = .*|initialbondname_length = $((${#initialbondname}))|
     s/initialbondname = .*/initialbondname = \"$initialbondname\"/
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/percentthreshold_flag = .*/percentthreshold_flag = $percentthreshold_flag/
     s/threshold_rmsd = .*/threshold_rmsd = $threshold_rmsd5/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s/accept_first = .*/accept_first = $accept_first/
     s/accept_worst = .*/accept_worst = $accept_worst/
     s/comparison_flag = .*/comparison_flag = .false./
     s/testheatmapSA_flag = .*/testheatmapSA_flag = $testtrajSAheatmap_flag/
     s/testtrajSA_flag = .*/testtrajSA_flag = $testtrajSA_flag/" <$currentPATH/$oldANALYSIS.f90 >$newPATH/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS

echo $(date) >> $newPATH/$bashout
/usr/bin/time -a -o $newPATH/$bashout -f "ANALYSIS05 %E  %U  %S  %P  %O" ./a.out
