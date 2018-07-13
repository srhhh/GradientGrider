#!/bin/bash

make -f make_makeGridwithNewTrajectories
make clean -f make_makeGridwithNewTrajectories
./a.out

##make -f make_checkNewTrajectorieswithMultipleGrids
##make clean -f make_checkNewTrajectorieswithMultipleGrids
##./a.out

