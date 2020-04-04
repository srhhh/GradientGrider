a.out : DOUBLE.o PARAMETERS.o FUNCTIONS.o PHYSICS.o VARIABLES.o ANALYSIS.o ls_rmsd_original.o SIMILARITY.o interactMultipleGrids.o mapCellData.o analyzeScatteringAngleswithMultipleGrids.o analyzeRMSDThresholdwithMultipleGrids.o analyzeHeatMapswithMultipleGrids.o runTrajectory.o makeGridwithNewTrajectories.o
	$(FC) -o a.out -O3 -fcheck="bound" DOUBLE.o PARAMETERS.o FUNCTIONS.o PHYSICS.o VARIABLES.o ANALYSIS.o ls_rmsd_original.o SIMILARITY.o interactMultipleGrids.o mapCellData.o analyzeScatteringAngleswithMultipleGrids.o analyzeRMSDThresholdwithMultipleGrids.o analyzeHeatMapswithMultipleGrids.o runTrajectory.o makeGridwithNewTrajectories.o

FC = gfortran
STEM = $(shell pwd)
SOURCE = $(STEM)/src/
EXPSOURCE = $(STEM)/startup/

DOUBLE.o : $(SOURCE)DOUBLE.f90
	$(FC) -c $(SOURCE)DOUBLE.f90
PARAMETERS.o : $(SOURCE)PARAMETERS.f90
	$(FC) -c $(SOURCE)PARAMETERS.f90
FUNCTIONS.o : $(SOURCE)FUNCTIONS.f90
	$(FC) -c $(SOURCE)FUNCTIONS.f90
VARIABLES.o : $(SOURCE)VARIABLES.f90
	$(FC) -c $(SOURCE)VARIABLES.f90
PHYSICS.o : $(SOURCE)PHYSICS.f90
	$(FC) -c $(SOURCE)PHYSICS.f90
ANALYSIS.o : $(EXPSOURCE)ANALYSIS.f90
	$(FC) -c $(EXPSOURCE)ANALYSIS.f90
SIMILARITY.o : $(SOURCE)SIMILARITY.f90
	$(FC) -c $(SOURCE)SIMILARITY.f90
ls_rmsd_original.o : $(SOURCE)ls_rmsd_original.f90
	$(FC) -c $(SOURCE)ls_rmsd_original.f90
#interactSingleGrid.o : $(SOURCE)interactSingleGrid.f90
#	$(FC) -c $(SOURCE)interactSingleGrid.f90
interactMultipleGrids.o : $(SOURCE)interactMultipleGrids.f90
	$(FC) -c $(SOURCE)interactMultipleGrids.f90
mapCellData.o : $(SOURCE)mapCellData.f90
	$(FC) -c $(SOURCE)mapCellData.f90
analyzeScatteringAngleswithMultipleGrids.o : $(SOURCE)analyzeScatteringAngleswithMultipleGrids.f90
	$(FC) -c $(SOURCE)analyzeScatteringAngleswithMultipleGrids.f90
analyzeRMSDThresholdwithMultipleGrids.o : $(SOURCE)analyzeRMSDThresholdwithMultipleGrids.f90
	$(FC) -fopenmp -c $(SOURCE)analyzeRMSDThresholdwithMultipleGrids.f90
analyzeHeatMapswithMultipleGrids.o : $(SOURCE)analyzeHeatMapswithMultipleGrids.f90
	$(FC) -c $(SOURCE)analyzeHeatMapswithMultipleGrids.f90
runTrajectory.o : $(SOURCE)runTrajectory.f90
	$(FC) -c $(SOURCE)runTrajectory.f90
makeGridwithNewTrajectories.o : $(SOURCE)makeGridwithNewTrajectories.f90
	$(FC) -c $(SOURCE)makeGridwithNewTrajectories.f90
clean :
	rm *.o *.mod
