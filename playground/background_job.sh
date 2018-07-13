#!/bin/bash

oldPARAMETERS=PARAMETERS
newPARAMETERS=PARAMETERS_new

oldMAKEGRID=make_makeGridwithNewTrajectories
newMAKEGRID=make_makeGridwithNewTrajectories_new

oldANALYSIS=ANALYSIS
ANALYSIS1=ANALYSIS_1
ANALYSIS2=ANALYSIS_2
ANALYSIS3=ANALYSIS_3
ANALYSIS4=ANALYSIS_4
ANALYSIS5=ANALYSIS_5
ANALYSIS6=ANALYSIS_6

oldMAKEANALYSIS=make_checkNewTrajectorieswithMultipleGrids
MAKEANALYSIS1=make_checkNewTrajectorieswithMultipleGrids1
MAKEANALYSIS2=make_checkNewTrajectorieswithMultipleGrids2
MAKEANALYSIS3=make_checkNewTrajectorieswithMultipleGrids3
MAKEANALYSIS4=make_checkNewTrajectorieswithMultipleGrids4
MAKEANALYSIS5=make_checkNewTrajectorieswithMultipleGrids5
MAKEANALYSIS6=make_checkNewTrajectorieswithMultipleGrids6

sed "s/Ntraj_max = [0-9]*/Ntraj_max = 10/
     s/overcrowd0 = [0-9]*/overcrowd0 = 50/
     s/scaling1_0 = [0-9]*/scaling1_0 = 4/
     s/scaling2_0 = [0-9]*/scaling2_0 = 4/
     s/Ngrid_max = [0-9]*/Ngrid_max = 1/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS.f90/" <$oldPARAMETERS.f90 >$newPARAMETERS.f90

sed "s/Ngrid_cap = [0-9]*/Ngrid_cap = 1/
     s/heatmap_flag = \\.false\\./heatmap_flag = .true./
     s/trueSA_flag = \\.true\\./trueSA_flag = .false./
     s/testtraj_flag = \\.false\\./testtraj_flag = .true./
     s/useolddata_flag = \\.true\\./useolddata_flag = .false./
     s/Ntesttraj = [0-9]*/Ntesttraj = 700/
     s/testtrajRMSD_flag = \\.true\\./testtrajRMSD_flag = .false./
     s/percentthreshold_flag = \\.false\\./percentthreshold_flag = .true./
     s/threshold_rmsd = .*/threshold_rmsd = 0.0000d0/
     s/reject_flag = \\.false\\./reject_flag = .true./
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$oldANALYSIS >$ANALYSIS1

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS.f90/
     s/$oldANALYSIS\\.o/$ANALYSIS1.o/
     s/$oldANALYSIS\\.f90/$ANALYSIS1.f90/" <$oldMAKEGRID >$newMAKEGRID

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s/[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]/0016_00050_0010/
     s/$oldANALYSIS\\.o/$ANALYSIS1.o/
     s/$oldANALYSIS\\.f90/$ANALYSIS1.f90/" <$oldMAKEANALYSIS >$MAKEANALYSIS1

make -f make_makeGridwithNewTrajectories
make clean -f make_makeGridwithNewTrajectories
./a.out

make -f make_checkNewTrajectorieswithMultipleGrids
make clean -f make_checkNewTrajectorieswithMultipleGrids
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
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$oldANALYSIS >$ANALYSIS2

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s/[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]/0016_00050_0010/
     s/$oldANALYSIS\\.o/$ANALYSIS2.o/
     s/$oldANALYSIS\\.f90/$ANALYSIS2.f90/" <$oldMAKEANALYSIS >$MAKEANALYSIS2

make -f make_checkNewTrajectorieswithMultipleGrids
make clean -f make_checkNewTrajectorieswithMultipleGrids
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
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$oldANALYSIS >$ANALYSIS3

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s/[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]/0016_00050_0010/
     s/$oldANALYSIS\\.o/$ANALYSIS2.o/
     s/$oldANALYSIS\\.f90/$ANALYSIS2.f90/" <$oldMAKEANALYSIS >$MAKEANALYSIS3

make -f make_checkNewTrajectorieswithMultipleGrids
make clean -f make_checkNewTrajectorieswithMultipleGrids
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
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$oldANALYSIS >$ANALYSIS4

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s/[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]/0016_00050_0010/
     s/$oldANALYSIS\\.o/$ANALYSIS2.o/
     s/$oldANALYSIS\\.f90/$ANALYSIS2.f90/" <$oldMAKEANALYSIS >$MAKEANALYSIS4

make -f make_checkNewTrajectorieswithMultipleGrids
make clean -f make_checkNewTrajectorieswithMultipleGrids
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
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$oldANALYSIS >$ANALYSIS5

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s/[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]/0016_00050_0010/
     s/$oldANALYSIS\\.o/$ANALYSIS2.o/
     s/$oldANALYSIS\\.f90/$ANALYSIS2.f90/" <$oldMAKEANALYSIS >$MAKEANALYSIS5

make -f make_checkNewTrajectorieswithMultipleGrids
make clean -f make_checkNewTrajectorieswithMultipleGrids
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
     s/testtrajSA_flag = \\.false\\./testtrajSA_flag = .true./" <$oldANALYSIS >$ANALYSIS6

sed "s/$oldPARAMETERS\\.o/$newPARAMETERS\\.o/
     s/$oldPARAMETERS\\.f90/$newPARAMETERS\\.f90/
     s/[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]/0016_00050_0010/
     s/$oldANALYSIS\\.o/$ANALYSIS2.o/
     s/$oldANALYSIS\\.f90/$ANALYSIS2.f90/" <$oldMAKEANALYSIS >$MAKEANALYSIS6

make -f make_checkNewTrajectorieswithMultipleGrids
make clean -f make_checkNewTrajectorieswithMultipleGrids
./a.out
