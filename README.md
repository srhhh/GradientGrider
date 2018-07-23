
GRADIENT GRIDER

An MD simulation program that stores reaction data over time for gradient approximations.

After downloading the repository, open up the bash script and make changes acccording to what reaction dynamics is under study. Do not alter the original source code (.f90s and makefiles).

After the changes are saved, simply execute the bash script.

------------------------------------------------------------------------------------
A.		DEFINITIONS
------------------------------------------------------------------------------------

A LIBRARY is a collection of grids with the same parameters; each library is associated with a directory. For example, "HH2_004_004_00050_01000/" is the name of a directory of a library in my local repository. The name of the directory associated with a library is initialized in the bash script. A library must be built up before it is analyzed. All analysis output (jpgs,dat,etc.) is kept in the library. There is also a local version of the source code in a library, with changes made to reflect different paths and parameters.

A GRID is a collection of frames of the phase space of a reaction; the reaction, its dynamics, and the way the reaction was explored are detailed in parameters kept by the grid's library. For now, all grids are numbered. For example, "001/" is the name of a directory associated with a grid. All frames are stored in cells. Frames in the phase space are explored through MD trajectory simulation; some data about these trajectories are also stored in the grid.

A TRAJECTORY is a collection of data about one particular reaction, obtained by simulating the molecular dynamics of the reactants from some random starting configuration to some cutoff. Each trajectory can be described by some set of initial conditions and some method of accelerating a frame. A file associated with a trajectory will have the prefix part of the filename describe what the acceleration method is, and will have the suffix part of the filename describe what kind of data is listed. For example, "reject_.00100_000592.dat" is the 592nd trajectory of some set of trajectories which used the "reject_.00100" method; it is a dat file so it has, at least, information on the RMSD of each frame of the trajectory compared to some other frame in the grid.

A CELL is a collection of frames of the phase space, all related by some set of variables as detailed in the grid's library parameters. A file associated with a cell will have a filename fully described by what variables the frames are related by. For example: "6.9925_7.0075.dat" would be the filename of a two-variable cell, where the first variable is approximated by 6.9925 and the second variable is approximated by 7.0075. These files are kept in a subdirectory of a grid, for example, they may be kept in "001/grid/", the grid subdirectory of "001/". One could say all frames in the cell were 'binned' into that cell. A cell also has an order and size, as specified in the grid's library parameters. Two cells in a grid may be subcells, neighbors, or unrelated depending on these parameters.

The ORDER of a cell is its level of granularity or coarseness. A cell that is more granular is smaller, and a cell that is more coarse is larger. For example, let cell A be "6.99_7.00.dat" and let cell B be "6.9925_7.0075.dat". Because the first and second variables in cell B have more significant digits than those of cell A, cell B is more granular than cell A; that is to say, cell A is more coarse than cell B. We can also expect cell A to be larger than cell B.

The SIZE of a cell is its range for acceptance of frames. For example, if cell "6.99_7.04.dat" has size .01 by .01, then it will accept all frames which have a first variable in the range [6.99,7.00] and a second variable in the range [7.04,7.05]. All cells with the same size are of the same order.

A cell is a SUBCELL of another cell if all frames in the subcell would be accepted in the other cell. For example, let cell A be "6.99_7.00.dat" with size .01 by .01 and let cell B be "6.9925_7.0075.dat" with size .0025 by .0025. Clearly, if a frame has a first variable in range [6.9925,6.9950] and second variable in range [7.0075,7.0100], then it would also be in range [6.99,7.00] and [7.00,7.01], respectively. Thus, cell B is a subcell of cell A.

A cell is a NEIGHBOR of another cell if both cells are of the same order and their ranges overlap. For example, let cell A be "6.99_7.00.dat" and let cell B be "7.00_7.01.dat", both with size .01 by .01. Because they have the same size, they must be of the same order. Cell A has range [6.99,7.00] x [7.00,7.01] while cell B has range [7.00,7.01] x [7.01,7.02]. Looking at this on a Cartesian grid, it is clear to see these two ranges overlap on a single point: 7.00 x 7.01. Thus, these two cells are neighbors.

A FRAME is a point in the phase space of a reaction. Thus, a frame is associated with a set of coordinates and velocities. The gradient of a frame is dependent only on the coordinates so a frame is also associated with a gradient. Also, using the coordinates of a frame, a set of variables can be calculated, as specified in the grid's library parameters. Thus, a frame is also associated with a set of variables. If a frame is in a cell, then the variables, coordinates, and gradient associated with the frame are recorded in the file associated with that cell. Take note that no information on the velocities of the frame are recorded. Finally, if the cell is overcrowded, then there may be no information recorded in the cell, and instead the frame is recorded in the file associated with a subcell of that cell.

A cell is OVERCROWDED if there are more than a certain number of frames in the cell. The specific number is in the grid's library parameters and depends on the order of the cell.

The set of VARIABLES associated with a frame are calculated using the coordinates of the frame. These may be thought of as the collective variables of the reaction. The methods used to calculate them are specified in the grid's library parameters.




------------------------------------------------------------------------------------
B.		BUILDING UP A LIBRARY
------------------------------------------------------------------------------------

A library always starts with just source code. Although everything about every grid may be initialized---every cell, their sizes, their order---there are no frames in the beginning. Frames may either come from an I/O stream (future implementation) or from an MD simulation. And after a frame is created, it is input into a grid.

Because putting a frame into a grid requires the initialization of a grid (wasted computation time), grid's are filled up one-by-one. Each grid can be thought of as its own program. So to initialize a grid, first, the fortran code needs to be compiled, and second, the internal arrays of the program need to be allocated. There is an additional trouble if the grid already exists, meaning that some of the internal counter do not start at zero, and instead their data needs to be read from some text file written beforehand. Finally, when the program ends, these internal arrays need to be written to a text file (in case for future use) and need to be deallocated. All these points contribute to wasted computation time.

Each grid is its own directory and is autonomous except in two regards. One, the grid uses source code from its parent library. In particular, it uses its parent library's parameters. These govern things that should essentially be static within a grid, like the physics and collision parameters of a reaction, to things that should be static among a set of grids, like the grid spacing and variables. This would thus throw a fatal error if attempting to move a grid from one library to another library with different parameters. However, it would work without qualms if renamed or moved to another library with the same parameters. Two, the grids are named by increments of one. Right now they are name 001/, 002/, 003/, ... If they are renamed to something not numerical, the library would not be able to identify the grid.

Libraries have no information on how many or what type of grids they hold. They make two assumption. One, the grids are named numerically with leading zeroes: 001/, 002/, 003/, so the library would not be able to identify a grid named otherwise. Two, the grids were all made using parameters as specified in the library's source code. This is essential because not only may the file and data formatting be different between two libraries but so too may the internal array sizes and indexing, which would result in a fatal flaw or, worse, frames being added to the wrong cells.



------------------------------------------------------------------------------------
B.		BUILDING UP A GRID
------------------------------------------------------------------------------------

When the grid-creation program starts running, the number of files it makes and variables, frames, and gradients on those files it makes may seem daunting. I will describe a simple example of a trajectory, what kinds of variables might be used, an example of a simple grid spacing of these variables, and an example of a simple subdivision of a cell that gets too full.



Let the reaction be of a hydrogen atom colliding into an H2 molecule. Let variable one be the distance of the incident hydrogen from one of the H2 hydrogen atoms, and let variable two be the distance of the incident hydrogen from the other H2 hydrogen atoms. Suppose that we are in grid 001/ and we have just started a trajectory. Let the first frame of the trajectory look like so:
```
     (1)               (2) 
      H - - - - - - - - H
                         \
                          H (3)
```
Here let variable one (var1) be the distance between H(1) and H(2) or r12, and variable two (var2) be the distance between H(1) and H(2) or r13. Suppose r12 = 2.916677 A and r13 = 3.135552 A.

Suppose the spacing of cells of order 1 be 1.0 A. If we run a heat map of order 1 cells, we would get something like this:

```
   (A)  
     #########################
   5 # 0 # 0 # 0 # 0 # 0 # 0 #
     #########################
   4 # 0 # 0 # 0 # 0 # 0 # 0 #
     #########################
   3 # 0 # 0 # 0 # 0 # 0 # 0 #
     #########################
   2 # 0 # 0 # 0 # 0 # 0 # 0 #
     #########################
   1 # 0 # 0 # 0 # 0 # 0 # 0 #
     #########################
   0 # 0 # 0 # 0 # 0 # 0 # 0 #
     #########################
       0   1   2   3   4   5   (A)
```
Each square surrounded by pound symbols is a cell. In the beginning, all cells are empty so the heat map should be all zeroes; this is because no frames have been input. A file is created for a cell only if a frame is inside the cell. So at this point, no files exist for order 1 cells.

Given some parameters in the physics parametes file, the gradient of this frame can be calcualted. Suppose that the first order is order 1 and the formatting of order 1 cells is F9.1 (a fortran format). Then the variables, the frame, and its gradient will be deposited into a file called 2.0_3.0.dat; this file is 001/grid/. This is because we only add a frame to one cell and we add it to the lowest order cell that is not overcrowded. In addition, the number of frames in 2.0_3.0.dat are recorded in an internal array called counter1.

If we do a heat map of order 1 cells now, then we would get something like this:

```
   (A)  
     #########################
   5 # 0 # 0 # 0 # 0 # 0 # 0 #
     #########################
   4 # 0 # 0 # 0 # 0 # 0 # 0 #
     #########################
   3 # 0 # 0 # 1 # 0 # 0 # 0 #
     #########################
   2 # 0 # 0 # 0 # 0 # 0 # 0 #
     #########################
   1 # 0 # 0 # 0 # 0 # 0 # 0 #
     #########################
   0 # 0 # 0 # 0 # 0 # 0 # 0 #
     #########################
       0   1   2   3   4   5   (A)
```

This would have created one new file. However, suppose that after the frame was added, the number of frames in the cell was different. Suppose the heat map looked like this instead:

```
   (A)  
     #########################
   5 # 0 # 0 # 0 # 0 # 0 # 0 #
     #########################
   4 # 0 # 0 # 0 # 0 # 0 # 0 #
     #########################
   3 # 0 # 0 #10 # 0 # 0 # 0 #
     #########################
   2 # 0 # 0 # 0 # 0 # 0 # 0 #
     #########################
   1 # 0 # 0 # 0 # 0 # 0 # 0 #
     #########################
   0 # 0 # 0 # 0 # 0 # 0 # 0 #
     #########################
       0   1   2   3   4   5   (A)
```

where 2.0_3.0.dat now clearly has 10 frames in it. Suppose that the overcrowding limit of the library for order 1 cells is 50 frames. This would then prompt the execution of divyUp on 2.0_3.0. Suppose for this library, order 1 cells scale by 4 with respect to var1 and 4 with respect to var2. This would produce subcells of order 2 with size .25 by .25. Before divyUp was called, this is what a heat map of order 2 subcells of 2.0_3.0 would have looked like:

```
  (A)  
      #########################
   3. #     #     #     #     #
   75 #  0  #  0  #  0  #  0  #
      #########################
   3. #     #     #     #     #
   50 #  0  #  0  #  0  #  0  #
      #########################
   3. #     #     #     #     #
   25 #  0  #  0  #  0  #  0  #
      #########################
   3. #     #     #     #     #
   00 #  0  #  0  #  0  #  0  #
      #########################
       2.00  2.25  2.50  2.75  (A)
```
Again, all cells are empty so no files exist for these. Because these are order 2 cells, they had lower priority than order 1 cells so they had no frames added to them. Back to our original frame. We had var1 = 2.916677 A and var2 = 3.135552 A. With truncation the program divyUp would put them in cell 2.75_3.00. divyUp will also place all other 9 frames in whichever order 2 subcell of 2.0_3.0 they belong to.

Although the file system is great, it is strenuous to have to inquire for files, open the file, count the number of frames, see if it is full, and if it is, repeat for a cell of higher order. Instead, we flatten these multidimensional grids into one-dimensional arrays that have cell indexes point to their subcells if they are too full. For instance, let's say after some time the grid looks like this:
```
  (A)
     #########################
   5 #31 #32 #33 #34 #35 #36 #
     #########################
   4 #25 #26 #27 #28 #29 #30 #
     #########################
   3 #19 #20 #21 #22 #23 #24 #
     #########################
   2 #13 #14 #15 #16 #17 #18 #
     #########################
   1 # 7 # 8 # 9 #10 #11 #12 #
     #########################
   0 # 1 # 2 # 3 # 4 # 5 # 6 #
     #########################
       0   1   2   3   4   5   (A)
```
Here, this is not a heat map. A cell on the grid now correponds to an internal index. The number of frames of some index is recorded in counter1. Because the limit for overcrowding of order 1 cells is 10 frames, we may get something like:
```
            #######################################################
   index1     #    19    #    20    #    21    #    22    # 
            #######################################################
   counter1   # 00000005 # 00000002 # 00001010 # 00000007 # 
            #######################################################
```
The above indicates that indexes 19, 20, and 22 have 5, 2, and 7 frames respectively. Because the limit of overcrowding is 10 frames, the frame count is suspended at 10 frames and the entire count is incremented by some multiple of some factor. Suppose for this library, the factor is 1000 and the multiple is 1. Thus, index 21 (cell 2.0_3.0) which has value 1010, is one multiple of 1 x 1000, so it was the first cell of order 1 to be overcrowded. Another example illustrates this. In the next example, more frames have been added and almost all cells are overcrowded:
```
            #######################################################
   index1     #    19    #    20    #    21    #    22    # 
            #######################################################
   counter1   # 00016010 # 00004010 # 00001010 # 00030010 # 
            #######################################################
```
Indexes 19, 20, 21, and 22 were all overcrowded so their frame count suspended at 10 frames. When divyUp was called on them, the value of their counters were incremented by some multiple of 1000. Because index 21 was the first overcrowded, it has 1 x 1000. Therefore, from the information above we can glean that index 20 was the fourth overcrowded (4 x 1000), index 19 was 16th (16 x 1000), and index 22 was 30th (30 x 1000). 

But now we must describe all of the order 2 subcells that were created from divyUp. How we recover index2 is more complicated.
```
          1st Overcrowded               30th Overcrwoded
  (A)  
      #########################     #########################
   3. #     #     #     #     #     #     #     #     #     #
   75 # 13  # 14  # 15  # 16  #     # 13  # 14  # 15  # 16  #
      #########################     #########################
   3. #     #     #     #     #     #     #     #     #     #
   50 #  9  # 10  # 11  # 12  #     #  9  # 10  # 11  # 12  #
      #########################     #########################
   3. #     #     #     #     #     #     #     #     #     #
   25 #  5  #  6  #  7  #  8  #     #  5  #  6  #  7  #  8  #
      #########################     #########################
   3. #     #     #     #     #     #     #     #     #     #
   00 #  1  #  2  #  3  #  4  #     #  1  #  2  #  3  #  4  #
      #########################     #########################
       2.00  2.25  2.50  2.75  (A)   3.00  3.25  3.50  3.75  (A)
```
These are again not heatmaps but the internal indexes representing each cell. It is important to realize these are relative indexes. Just imagine all the other order 1 cells neighboring 2.0_3.0. For example, index 22 corresponds to 3.0_3.0 which, according to counter1 is the 30th overcrowded. We need to be able to distinguish the subcell index 4 of 2.0_3.0 from the subcell index 4 of 1.0_3.0, and so on. We do this by spacing out the indexes by some factor. In this case, the factor is 16, the product of the scaling factors of var1 and var2.
```
             Index1 = Relative Index				(order 1)

      Overcrowded # = Counter1 (index 1) / 1000

             Index2 = (Overcrowded # - 1)* 16 + Relative Index	(order 2)
```
```
            ##########################################################
   filename   # 2.50_3.00 # 2.75_3.00 # 2.00_3.25 # 2.25_3.25 # 
            ##########################################################
   index2     #      3    #      4    #      5    #      6    # 
            ##########################################################
   counter2   #  00000002 #  00000001 #  00000000 #  00000003 # 
            ##########################################################

            ##########################################################
   filename   # 3.50_3.00 # 3.75_3.00 # 3.00_3.25 # 3.25_3.25 # 
            ##########################################################
   index2     #    467    #    468    #    469    #    470    # 
            ##########################################################
   counter2   #  00000000 #  00000000 #  00000001 #  00000006 # 
            ##########################################################
```
As you could imagine, the same trick can be used for the cells of order 2, 3, 4, and so on; cells of order 3 would be kept track of in counter3 with index3. The number of counters to use is kept at 3 in the program; some copy-pasting would be required to manually increase the number of counters. However, it is usually sufficient just to increase the scaling factor and overcrowding limit. A scaling factor of 4 by 4 with an initial spacing of .01 by .01 and order 1 overcrowding limit of 50 is usually sufficient to hold 700 trajectories worth of internal subcell pointers.


------------------------------------------------------------------------------------
D.		ANALYZING A LIBRARY
------------------------------------------------------------------------------------

Analyzing a library is split up into three parts: initializing the analysis, creating data to analyze, and plotting the data. Initializing is straightforward. In the bashscript the same variable used to name the library to-be-made is used to identify the library to-be-analyzed. Because there is a copy of the source code with the grid parameters in the library, nothing needs to be changed on the spot except for whatever outputs are wanted for the analysis. These outputs are listed as flags that can be changed in the bashscript sed commands. Creating the data to analyze takes the bulk of the time; this is why there are a significant number of trajectory DAT files that are written to; these files can just be read later instead of having to simulate random ones over again. Finally, plotting the data is also straightforward. This is done with gnuplot and the gnuplot files are made on the fly in fortran.

Most times, it is only useful to look at data per trajectory. Some trajectory-independent pieces of data like the heat maps of various cell orders are good descriptors of the grid, but not particularly of the reaction. When creating trajectory data, it is important to sample a system realistically. This means randomizing the trajectory initial conditions so that they sample the system at the correct temperature and such. 

Because the frames used to build the grids were from trajectories (most likely), those trajectories could also be analyzed. These are called the "initial trajectories" and did not check frames alongside the grid because the grid was not done by that point. They can be considered a control group. For example, because they used no approximations, their scattering angle distribution could be considered the "true" scattering angle distribution of the simulated reaction.

All other trajectories used in analysis are made with every frame being checked alongside the grids. For example, if the analysis specifies that every frame is checked with four grids, the frame is first checked with grid 001/, then 002/, 003/, and 004/. Each returns a frame that is closeby to the original frame and the gradient is used from the frame with the lowest RMSD for an approximation. Depending on the rejection method and the RMSD threshold, the gradient approximation is either used for that frame to accelerate in time or the real gradient is calculated and used. Because the number of grids, rejection method, and RMSD threshold all affect the choice of gradient, these three are used as prefixes for trajectory filenames. For example, "003reject_.00100_00004.dat" is the 4th trajectory of a set of trajectories that use three grids, always reject the gradient approximation, but still check whether there is a frame that has less than .001 RMSD.

Say that the analysis requires 100 trajectories with so and so gradient approximations methods. First, the prefix that corresponds to this gradient approximation method is calculated (the "003reject_.00100" in the example above), and then the program looks for the 100th trajectory with such prefix ("003reject_.00100_00100.dat"). Trajectories are made sequentially so if the 100th trajectory exists, then so to does the first, second, and so on until the 99th.


------------------------------------------------------------------------------------
E.		CHECKING A FRAME ALONGSIDE A GRID
------------------------------------------------------------------------------------

A frame is checked alongside a grid much the same way that a frame is added to a grid. However, instead of adding the frame to the target cell, the target cell is read. All the frames in the cell are compared to the original frame. In the grid creation, the target cell is checked whether it is overcrowded by checking the internal counter and following pointers. Something similar can be done for grid checking. However, if multiple grids are being checked at once---in which it would be laborious to have to open and close every counter text file---this internal counter system is not used. Instead, a subcell of some order is inquired for in the grid; if it does not exist, the cell that is one order higher is then inquired for. So while the internal counter system, by nature, checks from top-down, the external inquiry system checks from bottom-up.

There is not yet a default value for the RMSD if there is no frame found in the grid closeby to the original frame; it is hard-coded to be .200100 A. So RMSD thresholds at this point must be below .2001 A.
