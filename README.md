# Introduction

Gradient Grider is a data storage, retrieval, and interpolation program written in Fortran. Data processing is tailored for user-defined collective variables and XYZ coordinates: ideal for machine learning energy gradients and use in molecular dynamics simulations. Hyperparameters can be adjusted by changing the default variable values in the PARAMETERS file.

As Gradient Grider needs data to train and validate, supplementary modules are included to run molecular dynamics simulations using classical mechanics. Trajectory simulation parameters can be adjusted by changing respective values in the PHYSICS file.

As Gradient Grider needs to output and visualize data, supplementary modules are included to post-process the data and analyze observables from classical trajectorys, including scattering angle and energy partitioning distributions. Gnuplot is employed for graph creation. Analysis output parameters can be adjusted by changing respective values in the ANALYSIS file.

For the prototype, the PHYSICS file has default values ready for bimolecular collisions of three case studies: H - H2, H2 - H 2 , and HBr + - CO 2 . An accompanying bash script enables users to change the most important hyperparameters and analysis parameters easily, storing results from experiments using different parameters in different folders. Execution of the program is reduced to a single command from the terminal.

# Data Storage

Collective variables must be defined

Let the reaction be of a hydrogen atom colliding into an H2 molecule. Let variable one be the distance of the incident hydrogen from one of the H2 hydrogen atoms, and let variable two be the distance of the incident hydrogen from the other H2 hydrogen atoms. Suppose this is what one frame of the trajectory looks like:

     (1)               (2) 
      H - - - - - - - - H
                         \
                          H (3)
Here variable one (var1) is the distance between H(1) and H(2) or r12, and variable two (var2) is the distance between H(1) and H(2) or r13. Suppose r12 = 2.916677 A and r13 = 3.135552 A.

Our grid is only interested in close-range interactions so var1 and var2 are capped at 6 A. We let the spacing of spots in the grid be 1.0 A. So our grid, without subdivisions, looks like this:

   (A)  
     #########################
   5 #   #   #   #   #   #   #
     #########################
   4 #   #   #   #   #   #   #
     #########################
   3 #   #   # x #   #   #   #
     #########################
   2 #   #   #   #   #   #   #
     #########################
   1 #   #   #   #   #   #   #
     #########################
   0 #   #   #   #   #   #   #
     #########################
       0   1   2   3   4   5   (A)
In this case, it is easy to see where our frame should be; this is indicated by the spot with an x. Suppose that the format for files with no subdivision is F9.1 (a fortran format). Then the variables, the frame, and its gradient will be deposited into a file called 2.0_3.0.dat. In addition, the number of frames in 2.0_3.0.dat are recorded in an internal array called counter0.

Suppose that at the instance that frame was added, the number of frames exceeded the capacity allotted to the spot. This is natural and it prompts the execution of divyUp on 2.0_3.0. Suppose for this grid, we say we want var1 and var2 to both scale by 4 on their first subdivision. Then we will get a subgrid like this:

```
  (A)  
      #########################
   3. #     #     #     #     #
   75 #     #     #     #     #
      #########################
   3. #     #     #     #     #
   50 #     #     #     #     #
      #########################
   3. #     #     #     #     #
   25 #     #     #     #     #
      #########################
   3. #     #     #     #     #
   00 #     #     #     #  x  #
      #########################
       2.00  2.25  2.50  2.75  (A)
```

We had var1 = 2.916677 A and var2 = 3.135552 A. And when we deposit them in a spot, we always truncate. Thus, we would have them in spot 2.75_3.00. However, when divyUp is called, all frames that were in 2.0_3.0 are deposited in the more granular spots, so not just 2.75_3.00 will get a frame

# Data Retrieval

Talk about similarity

# Interpolation

# Example Execution

```
make
```
