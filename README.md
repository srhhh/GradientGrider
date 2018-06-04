To make the grid, enter in "make -f make_getCells" and start up the executable.
This has one main program: getCells.f90.
And four modules: f1_parameters.f90, f1_functions.f90, f1_variables.f90, addCells4.f90.

The main program goes roughly like follows:
Go into the directory containing the trajectories and pick a trajectory, one by one.
Use text manipulation to make it more fortran-friendly.
Read the formatted data frames; these are six-liners that describe an xyz coordinate and gradient.
Calculate the variables desired from these (right now, r1 and r2).
Bin a particular frame based on these variables and the desired spacing of the grid.
Add the state thoughout every cell and subcell of the grid it belongs to with "addState".
Repeat this for every frame of the trajectory, and every trajectory in the larger directory

"addState" is a subroutine of addCells4.f90
Taking as input the values of the variables and the xyz,gradient coordinates, bin the frame.
To 'bin' a frame, it simply appends onto a file the frame data.
The file is named to represent the bin; ex. 14.356, 3.759 -> 14.25_3.75.dat for first-order bins and 14.356_3.759 -> 14.350_3.750.dat for a second-order bin and so on.
To keep track of the number of frames in a bin, some array called counterN (where N is the number of places to the right of the decimal digit) is incremented by 1 every time a frame is added.
When a file has a specific number of states---overcrowd---"divyUp" is called.
If the file has more states than overcrowd, then a frame added to it also must be added to one of the file's subcells (ex. 14.3_3.7.dat is a subcell of 14._3..dat).
A subcell file with N digits to the right of the decimal place has its frame count in counterN, indexed by some indexN. IndexN can be calculated from the number in the N-th digit to the right of the decimal place (ex. 14.356, order=2 --> 5) and the first four digits of counterN-1(indexN-1).
This is possible because the number of frames in a cell is expected to be less than 10000, so occupies the last four digits of some integer string; the first four digits can be thus used to assign a unique indexN.
Not all possible subcells are accounted for in counterN; this means 'slots' are made on a first-come-first-serve basis.
So by using only the first four digits, 9999 unique indexes are allowed which is plenty for a resolution of 100; a new unique index is made for a frame that overcrowds and produces children and incremented (only once!) to counterN(indexN).
With such a procedure, a frame can be added to as many subcells as applicable.
As a temporary measure, indexer and counter are capped at order=3.

"divyUp" is a subroutine of addCells4.f90
Taking as input the truncated values of the variables, the subdirectory depth, and the current counterN/ it creates files with more digits than that of the original subcell.
It does this by first collecting all of the variable values and coordinates into an array.
Then it sorts the values, then it calculates which indexes belong to which subbin/subcell.
For each non-empty cell, it creates a new file r1_r2.dat where r1 and r2 now have N+1 digits to the right of the decimal.
The population of the new files are recorded in counterN+1 with IndexN+1.
IndexN+1 is calculated uniquely by multiplying indexN by the resolution and adding the unique integer between 0 and resolution for the subcell relative to other subcells (ex. for a 2x2 subcell, 1,1->1, 1,2->2, 2,1->3, 2,2->4).
This will insure the indexes of children of IndexN = a will not overlap with, let's say, indexes of children of IndexN = a + 1; each parent assigned an index will have distinct children from another parent.
If a rough estimate, K, is known of how many subcells in a particular order N will be overcrowded, then a tentative dimension for counterN would be something like K*resolution.

.

.

.

addCells.f90 is the original code (making use of .p files to keep track of populations)

.

addCells2.f90 is the secondary code, making use of arrays saved in subroutine addCells2 but separating the arrray holding indexes to go to a child subcell and the array holding the number of states of a parent cell.

.

addCells3.f90 is the tertiary code, making use of array counters that held both the key and the index to the child but with the archaic model of do-loop and a completely base-10 grid spacing model.

.

.
