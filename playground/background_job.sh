#!/bin/bash

oldPARAMETERS=PARAMETERS
newPARAMETERS=PARAMETERS_new

oldMAKEGRID=make_makeGridwithNewTrajectories
newMAKEGRID=make_makeGridwithNewTrajectories_new

oldANALYSIS=ANALYSIS
newANALYSIS=ANALYSIS_new

oldMAKEANALYSIS=make_checkNewTrajectorieswithMultipleGrids
newMAKEANALYSIS=make_checkNewTrajectorieswithMultipleGrids_new

##newGRID=0016_00050_0PATH
scaling1_0=004
scaling2_0=004
overcrowd0=00050
Ntraj_max=00200
Ngrid_max=1

newGRID=HH2${scaling1_0}_${scaling2_0}_${overcrowd0}_${Ntraj_max}
currentPATH=$(pwd)
gridPATH=$currentPATH/$newGRID
newSOURCE=SOURCE
newPATH=$(pwd)/$newGRID/$newSOURCE

if [ "0" -eq "0" ]
then

rm -r $currentPATH/$newGRID
mkdir $currentPATH/$newGRID
mkdir $newPATH/
cp $currentPATH/*.f90 $newPATH/
cp $currentPATH/make_$(echo "*") $newPATH/

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
     s|gridpath0 = .*|gridpath0 = \"$gridPATH/\"|
     s|$oldPARAMETERS\\.f90|$newPARAMETERS.f90|" <$currentPATH/$oldPARAMETERS.f90 >$newPATH/$newPARAMETERS.f90

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_max/
     s/heatmap_flag = \\.false\\./heatmap_flag = .true./
     s/trueSA_flag = \\.true\\./trueSA_flag = .false./
     s/testtraj_flag = \\.false\\./testtraj_flag = .true./
     s/useolddata_flag = \\.true\\./useolddata_flag = .false./
     s/Ntesttraj = [0-9]*/Ntesttraj = 100/
     s/testtrajRMSD_flag = \\.true\\./testtrajRMSD_flag = .false./
     s/percentthreshold_flag = \\.false\\./percentthreshold_flag = .true./
     s/threshold_rmsd = .*/threshold_rmsd = 0.0000d0/
     s/reject_flag = \\.false\\./reject_flag = .true./
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$currentPATH/$oldANALYSIS.f90 >$newPATH/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEGRID >$newPATH/$newMAKEGRID

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

cd $newGRID/

make -f $newPATH/$newMAKEGRID
make clean -f $newPATH/$newMAKEGRID
./a.out

make -f $newPATH/$newMAKEANALYSIS
make clean -f $newPATH/$newMAKEANALYSIS
./a.out

fi
exit

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_max/
     s/heatmap_flag = \\.false\\./heatmap_flag = .true./
     s/trueSA_flag = \\.false\\./trueSA_flag = .true./
     s/testtraj_flag = \\.false\\./testtraj_flag = .true./
     s/useolddata_flag = \\.true\\./useolddata_flag = .false./
     s/Ntesttraj = [0-9]*/Ntesttraj = 100/
     s/testtrajRMSD_flag = \\.true\\./testtrajRMSD_flag = .false./
     s/percentthreshold_flag = \\.false\\./percentthreshold_flag = .true./
     s/threshold_rmsd = .*/threshold_rmsd = 0.00000d0/
     s/reject_flag = \\.false\\./reject_flag = .true./
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$currentPATH/$oldANALYSIS.f90 >$newPATH/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS >$newPATH/$newMAKEANALYSIS

make -f $newPATH/$newMAKEANALYSIS
make clean -f $newPATH/$newMAKEANALYSIS
./a.out


#################################################################################################################################
#################################################################################################################################

exit

#################################################################################################################################
#################################################################################################################################


sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_max/
     s/heatmap_flag = \\.false\\./heatmap_flag = .true./
     s/trueSA_flag = \\.true\\./trueSA_flag = .false./
     s/testtraj_flag = \\.false\\./testtraj_flag = .true./
     s/useolddata_flag = \\.true\\./useolddata_flag = .false./
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntraj_max/
     s/testtrajRMSD_flag = \\.true\\./testtrajRMSD_flag = .false./
     s/percentthreshold_flag = \\.false\\./percentthreshold_flag = .true./
     s/threshold_rmsd = .*/threshold_rmsd = 0.00100d0/
     s/reject_flag = \\.true\\./reject_flag = .false./
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$currentPATH/$oldANALYSIS.f90 $newPATH/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS $newPATH/$newMAKEANALYSIS

make -f $newPATH/$newMAKEANALYSIS
make clean -f $newPATH/$newMAKEANALYSIS
./a.out

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_max/
     s/heatmap_flag = \\.false\\./heatmap_flag = .true./
     s/trueSA_flag = \\.true\\./trueSA_flag = .false./
     s/testtraj_flag = \\.false\\./testtraj_flag = .true./
     s/useolddata_flag = \\.true\\./useolddata_flag = .false./
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntraj_max/
     s/testtrajRMSD_flag = \\.true\\./testtrajRMSD_flag = .false./
     s/percentthreshold_flag = \\.false\\./percentthreshold_flag = .true./
     s/threshold_rmsd = .*/threshold_rmsd = 0.01000d0/
     s/reject_flag = \\.true\\./reject_flag = .false./
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$currentPATH/$oldANALYSIS.f90 $newPATH/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS $newPATH/$newMAKEANALYSIS

make -f $newPATH/$newMAKEANALYSIS
make clean -f $newPATH/$newMAKEANALYSIS
./a.out

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_max/
     s/heatmap_flag = \\.false\\./heatmap_flag = .true./
     s/trueSA_flag = \\.true\\./trueSA_flag = .false./
     s/testtraj_flag = \\.false\\./testtraj_flag = .true./
     s/useolddata_flag = \\.true\\./useolddata_flag = .false./
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntraj_max/
     s/testtrajRMSD_flag = \\.true\\./testtrajRMSD_flag = .false./
     s/percentthreshold_flag = \\.false\\./percentthreshold_flag = .true./
     s/threshold_rmsd = .*/threshold_rmsd = 0.10000d0/
     s/reject_flag = \\.true\\./reject_flag = .false./
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$currentPATH/$oldANALYSIS.f90 $newPATH/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS $newPATH/$newMAKEANALYSIS

make -f $newPATH/$newMAKEANALYSIS
make clean -f $newPATH/$newMAKEANALYSIS
./a.out

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = $Ngrid_max/
     s/heatmap_flag = \\.false\\./heatmap_flag = .true./
     s/trueSA_flag = \\.true\\./trueSA_flag = .false./
     s/testtraj_flag = \\.false\\./testtraj_flag = .true./
     s/useolddata_flag = \\.true\\./useolddata_flag = .false./
     s/Ntesttraj = [0-9]*/Ntesttraj = $Ntraj_max/
     s/testtrajRMSD_flag = \\.true\\./testtrajRMSD_flag = .false./
     s/percentthreshold_flag = \\.false\\./percentthreshold_flag = .true./
     s/threshold_rmsd = .*/threshold_rmsd = 0.20000d0/
     s/reject_flag = \\.true\\./reject_flag = .false./
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$currentPATH/$oldANALYSIS.f90 $newPATH/$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s|SOURCE = .*|SOURCE = $newPATH/|
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$currentPATH/$oldMAKEANALYSIS $newPATH/$newMAKEANALYSIS

make -f $newPATH/$newMAKEANALYSIS
make clean -f $newPATH/$newMAKEANALYSIS
./a.out
