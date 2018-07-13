#!/bin/bash

oldPARAMETERS=PARAMETERS
newPARAMETERS=PARAMETERS_new

oldMAKEGRID=make_makeGridwithNewTrajectories
newMAKEGRID=make_makeGridwithNewTrajectories_new

oldANALYSIS=ANALYSIS
newANALYSIS=ANALYSIS_new

oldMAKEANALYSIS=make_checkNewTrajectorieswithMultipleGrids
newMAKEANALYSIS=make_checkNewTrajectorieswithMultipleGrids_new

newGRID=0016_00050_0700

rm -r $newGrid

sed "s/Ntraj_max = [0-9]*/Ntraj_max = 700/
     s/overcrowd0 = [0-9]*/overcrowd0 = 50/
     s/scaling1_0 = [0-9]*/scaling1_0 = 4/
     s/scaling2_0 = [0-9]*/scaling2_0 = 4/
     s/Ngrid_max = [0-9]*/Ngrid_max = 3/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS.f90/" <$oldPARAMETERS.f90 >$newPARAMETERS.f90

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = 3/
     s/heatmap_flag = \\.false\\./heatmap_flag = .true./
     s/trueSA_flag = \\.true\\./trueSA_flag = .false./
     s/testtraj_flag = \\.false\\./testtraj_flag = .true./
     s/useolddata_flag = \\.true\\./useolddata_flag = .false./
     s/Ntesttraj = [0-9]*/Ntesttraj = 700/
     s/testtrajRMSD_flag = \\.true\\./testtrajRMSD_flag = .false./
     s/percentthreshold_flag = \\.false\\./percentthreshold_flag = .true./
     s/threshold_rmsd = .*/threshold_rmsd = 0.0000d0/
     s/reject_flag = \\.false\\./reject_flag = .true./
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$oldANALYSIS.f90 >$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS.f90/
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$oldMAKEGRID >$newMAKEGRID

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s/[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]/$newGrid/
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$oldMAKEANALYSIS >$newMAKEANALYSIS

make -f $newMAKEGRID
make clean -f $newMAKEGRID
./a.out

make -f $newMAKEANALYSIS
make clean -f $newMAKEANALYSIS
./a.out

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = 3/
     s/heatmap_flag = \\.false\\./heatmap_flag = .true./
     s/trueSA_flag = \\.true\\./trueSA_flag = .false./
     s/testtraj_flag = \\.false\\./testtraj_flag = .true./
     s/useolddata_flag = \\.true\\./useolddata_flag = .false./
     s/Ntesttraj = [0-9]*/Ntesttraj = 700/
     s/testtrajRMSD_flag = \\.true\\./testtrajRMSD_flag = .false./
     s/percentthreshold_flag = \\.false\\./percentthreshold_flag = .true./
     s/threshold_rmsd = .*/threshold_rmsd = 0.00010d0/
     s/reject_flag = \\.true\\./reject_flag = .false./
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$oldANALYSIS.f90 >$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s/[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]/$newGrid/
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$oldMAKEANALYSIS >$newMAKEANALYSIS

make -f $newMAKEANALYSIS
make clean -f $newMAKEANALYSIS
./a.out

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = 3/
     s/heatmap_flag = \\.false\\./heatmap_flag = .true./
     s/trueSA_flag = \\.true\\./trueSA_flag = .false./
     s/testtraj_flag = \\.false\\./testtraj_flag = .true./
     s/useolddata_flag = \\.true\\./useolddata_flag = .false./
     s/Ntesttraj = [0-9]*/Ntesttraj = 700/
     s/testtrajRMSD_flag = \\.true\\./testtrajRMSD_flag = .false./
     s/percentthreshold_flag = \\.false\\./percentthreshold_flag = .true./
     s/threshold_rmsd = .*/threshold_rmsd = 0.00100d0/
     s/reject_flag = \\.true\\./reject_flag = .false./
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$oldANALYSIS.f90 >$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s/[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]/$newGrid/
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$oldMAKEANALYSIS >$newMAKEANALYSIS

make -f $newMAKEANALYSIS
make clean -f $newMAKEANALYSIS
./a.out

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = 3/
     s/heatmap_flag = \\.false\\./heatmap_flag = .true./
     s/trueSA_flag = \\.true\\./trueSA_flag = .false./
     s/testtraj_flag = \\.false\\./testtraj_flag = .true./
     s/useolddata_flag = \\.true\\./useolddata_flag = .false./
     s/Ntesttraj = [0-9]*/Ntesttraj = 700/
     s/testtrajRMSD_flag = \\.true\\./testtrajRMSD_flag = .false./
     s/percentthreshold_flag = \\.false\\./percentthreshold_flag = .true./
     s/threshold_rmsd = .*/threshold_rmsd = 0.01000d0/
     s/reject_flag = \\.true\\./reject_flag = .false./
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$oldANALYSIS.f90 >$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s/[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]/$newGrid/
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$oldMAKEANALYSIS >$newMAKEANALYSIS

make -f $newMAKEANALYSIS
make clean -f $newMAKEANALYSIS
./a.out

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = 3/
     s/heatmap_flag = \\.false\\./heatmap_flag = .true./
     s/trueSA_flag = \\.true\\./trueSA_flag = .false./
     s/testtraj_flag = \\.false\\./testtraj_flag = .true./
     s/useolddata_flag = \\.true\\./useolddata_flag = .false./
     s/Ntesttraj = [0-9]*/Ntesttraj = 700/
     s/testtrajRMSD_flag = \\.true\\./testtrajRMSD_flag = .false./
     s/percentthreshold_flag = \\.false\\./percentthreshold_flag = .true./
     s/threshold_rmsd = .*/threshold_rmsd = 0.10000d0/
     s/reject_flag = \\.true\\./reject_flag = .false./
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$oldANALYSIS.f90 >$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s/[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]/$newGrid/
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$oldMAKEANALYSIS >$newMAKEANALYSIS

make -f $newMAKEANALYSIS
make clean -f $newMAKEANALYSIS
./a.out

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = 3/
     s/heatmap_flag = \\.false\\./heatmap_flag = .true./
     s/trueSA_flag = \\.true\\./trueSA_flag = .false./
     s/testtraj_flag = \\.false\\./testtraj_flag = .true./
     s/useolddata_flag = \\.true\\./useolddata_flag = .false./
     s/Ntesttraj = [0-9]*/Ntesttraj = 700/
     s/testtrajRMSD_flag = \\.true\\./testtrajRMSD_flag = .false./
     s/percentthreshold_flag = \\.false\\./percentthreshold_flag = .true./
     s/threshold_rmsd = .*/threshold_rmsd = 0.20000d0/
     s/reject_flag = \\.true\\./reject_flag = .false./
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$oldANALYSIS.f90 >$newANALYSIS.f90

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s/[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]/$newGrid/
     s/$oldANALYSIS\\.o/$newANALYSIS.o/
     s/$oldANALYSIS\\.f90/$newANALYSIS.f90/" <$oldMAKEANALYSIS >$newMAKEANALYSIS

make -f $newMAKEANALYSIS
make clean -f $newMAKEANALYSIS
./a.out
