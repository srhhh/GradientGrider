# GradientGrider
getCells.f90 searches through the trajectories (in your folder B0) and puts them into the grid folder
this makes use of f1_parameters, f1_variables
then getCells.f90 searches every subcell and subdivides it if it has too many states; this creates a new directory named after the subcell
this makes use of f1_paramters, f1_functions, and addCells

in the works:
getHeatMapData gets the data for a heatmap; very simple
ls_rmsd calculates the RMSD of two states (from Bill Research Group)
checkCells sees which cell a given state is in and procures the closest (in RMSD) state to it

everything else is not very intersting
