# Introduction

Gradient Grider is a data storage, retrieval, and interpolation program written in Fortran. Data processing is tailored for user-defined collective variables and XYZ coordinates: ideal for machine learning energy gradients and use in molecular dynamics simulations. Hyperparameters can be adjusted by changing the default variable values in the PARAMETERS file.

As Gradient Grider needs data to train and validate, supplementary modules are included to run molecular dynamics simulations using classical mechanics. Trajectory simulation parameters can be adjusted by changing respective values in the PHYSICS file.

As Gradient Grider needs to output and visualize data, supplementary modules are included to post-process the data and analyze observables from classical trajectorys, including scattering angle and energy partitioning distributions. Gnuplot is employed for graph creation. Analysis output parameters can be adjusted by changing respective values in the ANALYSIS file.

For the prototype, the PHYSICS file has default values ready for bimolecular collisions of three case studies: H - H2, H2 - H 2 , and HBr + - CO 2 . An accompanying bash script enables users to change the most important hyperparameters and analysis parameters easily, storing results from experiments using different parameters in different folders. Execution of the program is reduced to a single command from the terminal.

# Data Storage

Gradient Grider takes multidimensional data gathered from a program (namely, a molecular dynamics simulation) and stores it to a file system. Each data point consists of pairs of (x, f(x)), with the ultimate goal of learning and predicting the function f. Because the file system persists, it may be moved to a hard drive, reused, and repurposed for program execution on other computers and systems.

This means a lot of data must be stored permanantly.

This means an efficient storage (and consequently, retrival) system must be employed.

The user may define any number of collective variables to help store data. For simplicity, and without loss of generality, let there be two collective variables (as it is by default).

Let the system of interest be of a hydrogen atom colliding into an H2 molecule and the function of interest be that which takes the coordinates and returns the energy gradient. Let variable one be the distance of the incident hydrogen from one of the H2 hydrogen atoms, and let variable two be the distance of the incident hydrogen from the other H2 hydrogen atoms. Suppose this is what one frame of a trajectory looks like:

     (1)               (2) 
      H - - - - - - - - H
                         \
                          H (3)

Here variable one (var1) is the distance between H(1) and H(2) or r12, and variable two (var2) is the distance between H(1) and H(2) or r13. Suppose r12 = 2.916677 A and r13 = 3.135552 A.

For this system, only close-range interactions are of interest so var1 and var2 are capped at 6 A. And let the spacing of data in the grid be 1.0 A. So the grid, without subdivisions, looks like this:

```
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
```

In this case, it is easy to see where the frame should be; this is indicated by the spot with an x. Suppose that the format for files with no subdivision is F9.1 (a fortran format). Then the variables, the frame, and its gradient will be deposited into a file called 2.0_3.0.dat.

Suppose that at the instance that frame was added, the number of frames exceeded the capacity allotted to the spot. This is natural and it prompts the execution of divyUp on 2.0_3.0. Suppose for this grid, it is specified in the parameters to have var1 and var2 to both scale by 4 on their first subdivision. Then this produces a subgrid like this:

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

Originally, var1 = 2.916677 A and var2 = 3.135552 A. When these are deposited in a file in a grid, the file name results from truncation. Thus, this specific frame would be deposited in the file 2.75_3.00.dat. However, when a divyUp is called, all frames that were in 2.0_3.0 are deposited in the more granular spots, so not just 2.75_3.00 will get a frame

# Data Retrieval

Here, the accumulation of data in the grid is to enlargening the training set for a machine learning algorithm. Here, the ultimate goal is to predict values of f(x) given x. For molecular dynamics and the previous example, x is a frame of a trajectory. To predict f(x), similar frames must be acquired by checking the grid.

The longer a simulation progress, the larger the file system will grow and the more difficult it will be to check every data point in the grid. Collective variables were defined prior to the creation of the grid to help store the data and in a way so that no individual file within the larger grid is too large. This enables the search for frames within the grid that are similar to some requested frame--at least in terms of the collective variables.

After a search through the grid for frames that are similar enough--within a threshold--to the requested frame, additional pruning is done through similarity metrics. Because the energy of a frame is invariant to translation and rotation, and consequently the energy gradient is invariant up to a rotation, a metric which evaluates how similar frames are under these invariants is preferred. The root mean square distance (RMSD) and difference between Coloumb matrices (CMD) are two default similarity metrics used.

# Interpolation

The infrastructure setup from the previous two sections enables the machine learning algorithm to use frames similar to the requested frame to train. This is ideal for simple regression models. For more complicated models, the similarity thresholds can be relaxed and a larger, more diverse set of frames may be used to train.

A type of regression model is used here, where an interpolation of similar frames is used to approximate the requested frame, but with a regularization so that there is less overfitting. The amount of regularization, the similarity threshold, and other hyperparameters can be adjusted on-the-fly during the simulation. If there are not enough similar frames for a valid interpolation, an interpolation on the requested frame can be rejected.

When there is enough data in the grid, a subset of the training set can be used as a validation set. Error is measured with an L2-norm on the energy gradients after applying the appropriate rotations output from an RMSD call. In this way, error can be measured on-the-fly during a simulation. When a desired level of accuracy is required, a threshold on the error may be provided as well as a required level of sampling from the validation set. This enables Gradient Grider to sample preiodically the level of error encountered in various areas of the grid, pre-emptively rejecting frame requests in areas that have too much error uncertainty.

Whenever a frame submits a request for an energy gradient from Gradient Grider, either it is rejected or allowed to proceed. Following a rejection, the molecular dynamics simulation must provide its own energy gradient evaluation to allow the trajectory to proceed.

# Example Execution

Clone the master branch from github into your desired workspace.

```
make
```
