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

#The ratio of child-level to parent-level cell spacings with respect to var1
scaling1_0=004

#The ratio of child-level to parent-level cell spacings with respect to var2
scaling2_0=004

#The number of frames a parent-level cell accepts before it is subdivided
overcrowd0=00050

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
useolddata_flag=.true.
testtrajRMSD_flag=.false.
percentthreshold_flag=.true.
#threshold_rmsd=.200100d0
threshold_rmsd=.004000d0
threshold_rmsd1=.004000d0
threshold_rmsd2=.001000d0
threshold_rmsd3=.005000d0
threshold_rmsd4=.010000d0
threshold_rmsd5=.050000d0
reject_flag=.false.
accept_first=.false.
testtrajSA_flag=.true.
Ngrid_cap=4
#Ngrid_cap=${Ngrid_max}
Ntrajectories=350

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
newGRID_check_min=30

#The number of post-grid analyses you would like done
Nanalyses=1

#The path that has the original source code
currentPATH=$(pwd)

#DO NOT TOUCH THESE
gridPATH=$currentPATH/$newGRID
newSOURCE=SOURCE
newPATH=$(pwd)/$newGRID/$newSOURCE

###############################################################################################################################################
###############################################################################################################################################
#		GRID CREATION AND SUMMARY ANALYSIS
###############################################################################################################################################
###############################################################################################################################################

#Set this true if you want to create a new grid
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
     s|scaling1_0 = [0-9]*|scaling1_0 = $scaling1_0|
     s|scaling2_0 = [0-9]*|scaling2_0 = $scaling2_0|
     s|Ngrid_max = [0-9]*|Ngrid_max = $Ngrid_max|
     s|trajectory_text_length = [0-9]*|trajectory_text_length = ${#Ntraj_max}|
     s|scaling1_text_length = [0-9]*|scaling1_text_length = ${#scaling1_0}|
     s|scaling2_text_length = [0-9]*|scaling2_text_length = ${#scaling2_0}|
     s|overcrowd0_text_length = [0-9]*|overcrowd0_text_length = ${#overcrowd0}|
     s|gridpath_length = .*|gridpath_length = $((${#gridPATH}+1))|
     s|Ngrid_check_min = .*|Ngrid_check_min = $newGRID_check_min|
     s|force_Duplicates = .*|force_Duplicates = $force_Duplicates|
     s|force_NoLabels = .*|force_NoLabels = $force_NoLabels|
     s|gridpath0 = .*|gridpath0 = \"$gridPATH/\"|
     s|$oldPARAMETERS\\.f90|$newPARAMETERS.f90|" <$currentPATH/$oldPARAMETERS.f90 >$newPATH/$newPARAMETERS.f90

#Make changes to the analysis file as specified in the variables above
#This first analysis simulates new trajectories and checks them with the grid
#All of these variables can be changed;
#Ntesttraj is how many trajectorie to simulate (important!)
sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_cap/
     s/heatmap_flag = .*/heatmap_flag = .true./
     s/trueSA_flag = .*/trueSA_flag = .true./
     s/trueED_flag = .*/trueED_flag = .true./
     s/testtraj_flag = .*/testtraj_flag = .false./
     s/useolddata_flag = .*/useolddata_flag = $useolddata_flag/
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/percentthreshold_flag = .*/percentthreshold_flag = $percentthreshold_flag/
     s/threshold_rmsd = .*/threshold_rmsd = $threshold_rmsd/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s/accept_first = .*/accept_first = $accept_first/
     s/testtrajSA_flag = .*/testtrajSA_flag = $testtrajSA_flag/" <$currentPATH/$oldANALYSIS.f90 >$newPATH/$newANALYSIS.f90

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
./a.out

make -f $newPATH/$newMAKEANALYSIS
./a.out

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
cp $currentPATH/!($oldPARAMETERS|$newPARAMETERS)+(.f90) $newPATH/
cp $currentPATH/make_$(echo "*") $newPATH/
shopt -s extglob

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_cap/
     s/heatmap_flag = .*/heatmap_flag = $heatmap_flag/
     s/trueSA_flag = .*/trueSA_flag = $trueSA_flag/
     s/trueED_flag = .*/trueED_flag = $trueED_flag/
     s/testtraj_flag = .*/testtraj_flag = $testtraj_flag/
     s/useolddata_flag = .*/useolddata_flag = $useolddata_flag/
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/percentthreshold_flag = .*/percentthreshold_flag = $percentthreshold_flag/
     s/threshold_rmsd = .*/threshold_rmsd = $threshold_rmsd1/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s/accept_first = .*/accept_first = $accept_first/
     s/testtrajSA_flag = .*/testtrajSA_flag = $testtrajSA_flag/" <$currentPATH/$oldANALYSIS.f90 >$newPATH/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

cd $newPATH

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS
./a.out

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
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/percentthreshold_flag = .*/percentthreshold_flag = $percentthreshold_flag/
     s/threshold_rmsd = .*/threshold_rmsd = $threshold_rmsd2/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s/accept_first = .*/accept_first = $accept_first/
     s/testtrajSA_flag = .*/testtrajSA_flag = $testtrajSA_flag/" <$currentPATH/$oldANALYSIS.f90 >$newPATH/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS
./a.out

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
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/percentthreshold_flag = .*/percentthreshold_flag = $percentthreshold_flag/
     s/threshold_rmsd = .*/threshold_rmsd = $threshold_rmsd3/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s/accept_first = .*/accept_first = $accept_first/
     s/testtrajSA_flag = .*/testtrajSA_flag = $testtrajSA_flag/" <$currentPATH/$oldANALYSIS.f90 >$newPATH/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS
./a.out

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
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/percentthreshold_flag = .*/percentthreshold_flag = $percentthreshold_flag/
     s/threshold_rmsd = .*/threshold_rmsd = $threshold_rmsd4/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s/accept_first = .*/accept_first = $accept_first/
     s/testtrajSA_flag = .*/testtrajSA_flag = $testtrajSA_flag/" <$currentPATH/$oldANALYSIS.f90 >$newPATH/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS
./a.out

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
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntrajectories/
     s/testtrajRMSD_flag = .*/testtrajRMSD_flag = $testtrajRMSD_flag/
     s/percentthreshold_flag = .*/percentthreshold_flag = $percentthreshold_flag/
     s/threshold_rmsd = .*/threshold_rmsd = $threshold_rmsd5/
     s/reject_flag = .*/reject_flag = $reject_flag/
     s/accept_first = .*/accept_first = $accept_first/
     s/testtrajSA_flag = .*/testtrajSA_flag = $testtrajSA_flag/" <$currentPATH/$oldANALYSIS.f90 >$newPATH/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

make clean -f $newPATH/$newMAKEANALYSIS
make -f $newPATH/$newMAKEANALYSIS
./a.out
