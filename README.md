Gradient Grider now has 3 different programs under its wings.

make_getCells is the makefile of getCells.f90. This program was designed specifically to extract coordinate and gradients from a folder 'B0' of trajectories involving the reaction of CH3I and F. Although the framework is for trajectories of B0, the reading in of frames is dictated by a sed command that can be altered for different trajectory formats. In any case, once a frame is read, variables are calculated that describe the frame by f1_parameters.f90. With the frame and its variable, this data is deposited into a grid; all grid additions and subdivision are supervised by addCells.f90.

make_test_checkCells is the makefile of test_checkCells.f90. This program checks for a frame with a low RMSD compared to some input frame; the program is designed specifically for grids formatted by addCells.f90. The thinking is that if the gradient of some frame is unknown but a similar frame--with low RMSD compared to it--has a known frame, that frame's gradient can be used instead.

make_makeTrajectory is the makefile of initialTrajectories.f90. This program randomly initializes an md program involving the collision of H2 and H, runs it, and deposits all frames into a grid with addCells.f90. The module makeTrajectory.f90 actually runs the md where as f2_physics_parameters.f90 determines the md parameters and f2_variables.f90 determines what variables are used to describe the frame.

Gradient Grider has 3 main modules/

addCells.f90 includes the subroutines addState and divyUp.

addState takes a frame, its gradient, and some variables that describe the frame, and with this data, finds a 'spot' (essentially, a file) for it in the grid, whereafter it writes the frame and gradients onto the spot. Each spot in the grid represents some value the variables used to describe the frame can take. How these spots are represented, how they are laid out, and the how full we let them be are determined by f1_parametesr.f90. When a spot is too full, the spot is subdivided into more spots, each which describe the value of the variable in more granularity; this is done through a call to divyUp.

divyUp takes a spot in the grid and subdivides it into more spots, each which are more granular than the original spot. It then deposits into those spots whichever frames belong to them from the orignal spot.

f1/f2_variables.f90 includes subroutines getDistanceSquared, getVar1, getVar2, and so on.

getVarN takes a frame and calculates some variable of that frame. Often time, a distance between two atoms is used as a variable, and describes fairly well the reaction coordinate.

f1/f2_parameters.f90 include parameters controlling grid depositing and subdivision. Here, files and filepaths are written. Formats for how to describe a variable (like how many decimal places for a particular subdivision) are described. Numbers describing how many frames is too many for a spot, how many spots to subdiide a full spot into, how granular the spacing between spots should be, and so on are described.







When one of the grid-creation programs starts running, the number of files it makes and varaibles, frames, and gradients on those files it makes may seem daunting. I will describe a simple example of a trajectory, what kinds of variables might be used, an example of a simple griding of these variables, and an example of a simple subdivision of a spot that gets too full.

Let the reaction be of a hydrogen atom colliding into an H2 molecule. Let variable one be the distance of the incident hydrogen from one of the H2 hydrogen atoms, and let variable two be the distance of the incident hydrogen from the other H2 hydrogen atoms. Suppose this is what one frame of the trajectory looks like:

     (1)               (2) 
      H - - - - - - - - H
                         \
                          H (3)

Here variable one (var1) is the distance between H(1) and H(2) or r12, and variable two (var2) is the distance between H(1) and H(2) or r13. Suppose r12 = 2.916677 A and r13 = 3.135552 A.

Our grid is only interested in close-range interactions so var1 and var2 are capped at 5 A. We let the spacing of spots in the grid be 1.0 A. So our grid, without subdivisions, looks like this:


# (A)  
#   #########################
# 5 #   #   #   #   #   #   #
#   #########################
# 4 #   #   #   #   #   #   #
#   #########################
# 3 #   #   # x #   #   #   #
#   #########################
# 2 #   #   #   #   #   #   #
#   #########################
# 1 #   #   #   #   #   #   #
#   #########################
# 0 #   #   #   #   #   #   #
#   #########################
#     0   1   2   3   4   5   (A)

In this case, it is easy to see where our frame should be; this is indicated by the spot with an x. Suppose that the format for files with no subdivision is F9.1 (a fortran format). Then the variables, the frame, and its gradient will be deposited into a file called 2.0_3.0.dat. In addition, the number of frames in 2.0_3.0.dat are recorded in an internal array called counter0.

Suppose that at the instance that frame was added, the number of frames exceeded the capacity allotted to the spot. This is natural and it prompts the execution of divyUp on 2.0_3.0. Suppose for this grid, we say we want var1 and var2 to both scale by 4 on their first subdivision. Then we will get a subgrid like this:

#(A)  
#    #########################
# 3. #     #     #     #     #
# 75 #     #     #     #     #
#    #########################
# 3. #     #     #     #     #
# 50 #     #     #     #     #
#    #########################
# 3. #     #     #     #     #
# 25 #     #     #     #     #
#    #########################
# 3. #     #     #     #     #
# 00 #     #     #     #  x  #
#    #########################
#     2.00  2.25  2.50  2.75  (A)

We had var1 = 2.916677 A and var2 = 3.135552 A. And when we deposit them in a spot, we always truncate. Thus, we would have them in spot 2.75_3.00. However, when divyUp is called, all frames that were in 2.0_3.0 are deposited in the more granular spots, so not just 2.75_3.00 will get a frame

Although the file system is great, it is strenuous to have to inquire for files, open the file, count the number of frames, see if it is full, and if it is, repeat for a file of a subdivision. Instead, we flatten these multidimensional grids into one-dimensional arrays that point to their subdivisions if they are too full. For instance:

#(A)
#   #########################
# 5 #31 #32 #33 #34 #35 #36 #
#   #########################
# 4 #25 #26 #27 #28 #29 #30 #
#   #########################
# 3 #19 #20 #21 #22 #23 #24 #
#   #########################
# 2 #13 #14 #15 #16 #17 #18 #
#   #########################
# 1 # 7 # 8 # 9 #10 #11 #12 #
#   #########################
# 0 # 1 # 2 # 3 # 4 # 5 # 6 #
#   #########################
#     0   1   2   3   4   5   (A)


A spot on the grid would correpond to an index. The number of frames of some index is recorded in counter0. Let's say that a file is too full when it has 10 frames, then after the subdivision, we would get:

#          #######################################################
# index0    ... #    19    #    20    #    21    #    22    # ...
#          #######################################################
# counter0  ... # 00000005 # 00000002 # 00001010 # 00000007 # ...
#          #######################################################

The above indicates that indexes 19, 20, and 22 gave 5, 2, and 7 frames respectively. Because index 10 was subdivided, the frame count is suspended at 10 frames and then entire count is incremented by some multiple of some factor. In this case the factor is 1000 and the multiple is 1 because index 21 is the first spot to have been subdivided. Another example illustrates this. In the next example, more frames have been added and almost all spots on the grid have been subdivided:

#          #######################################################
# index0    ... #    19    #    20    #    21    #    22    # ...
#          #######################################################
# counter0  ... # 00016010 # 00004010 # 00001010 # 00030010 # ...
#          #######################################################

Indexes 19, 20, 21, and 22 were all subdivided so their frame count suspended at 10 frames. When they were subdivided they were incremented by some multiple of 1000. Because index 21 was subdivided first, it has 1000. Therefore, from the information above we can glean that index 20 was subdivided fourth, index 19 was subdivided 16th, and index 22 was subdivided 30th. 

Now that we have described the condition of the top-level grid completely, we need to describe the subdivided spots. We essentially first describe a spot by how many subdivisions it has undergone. Because 2.75_3.00 is the result of one subdivision, it is in counter1, whereas 2.0_3.0 has no subdivision so it is in counter0. How we recover index1 is more complicated.

#         SUBDIVISION 1                  SUBDIVISION 30
#(A)  
#    #########################     #########################
# 3. #     #     #     #     #     #     #     #     #     #
# 75 # 13  # 14  # 15  # 16  #     # 13  # 14  # 15  # 16  #
#    #########################     #########################
# 3. #     #     #     #     #     #     #     #     #     #
# 50 #  9  # 10  # 11  # 12  #     #  9  # 10  # 11  # 12  #
#    #########################     #########################
# 3. #     #     #     #     #     #     #     #     #     #
# 25 #  5  #  6  #  7  #  8  #     #  5  #  6  #  7  #  8  #
#    #########################     #########################
# 3. #     #     #     #     #     #     #     #     #     #
# 00 #  1  #  2  #  3  #  4  #     #  1  #  2  #  3  #  4  #
#    #########################     #########################
#     2.00  2.25  2.50  2.75  (A)   3.00  3.25  3.50  3.75  (A)

Although we index 2.75_3.00 just as we did 2.0_3.0, here the index 4 is a relative index. Just imagine all the other first-order subdivisions above, below, to the left, and to the right of 2.0_3.0. For example, index 22 corresponds to 3.0_3.0 which, according to counter0 is the 30th subdivision. We need to be able to distinguish index 4 of subdivision 1 from index 4 of subdivison 30, and so on. We do this by spacing out the indexes by some factor. In this case, the factor is 16, the product of the scaling factors of var1 and var2.

             Index0 = Relative Index (no subdivision)

      Subdivision # = Counter0 (index 0) / 1000

             Index1 = (Subdivision # - 1)* 16 + Relative Index (1 subdivision)

#          ##########################################################
# filename  ... # 2.50_3.00 # 2.75_3.00 # 2.00_3.25 # 2.25_3.25 # ...
#          ##########################################################
# index1    ... #      3    #      4    #      5    #      6    # ...
#          ##########################################################
# counter1  ... #  00000002 #  00000001 #  00000000 #  00000003 # ...
#          ##########################################################

#          ##########################################################
# filename  ... # 3.50_3.00 # 3.75_3.00 # 3.00_3.25 # 3.25_3.25 # ...
#          ##########################################################
# index1    ... #    467    #    468    #    469    #    470    # ...
#          ##########################################################
# counter1  ... #  00000000 #  00000000 #  00000001 #  00000006 # ...
#          ##########################################################

As you could imagine, the same trick can be used to subdivide further; files with two subdivision would be kept track of in counter2 with index2. The number of counters to use is kept at 3 in the program; some copy-pasting would be required to manually increase the number of counters. However, it is usually sufficient just to increase the scaling factor and file cap. A scaling factor of 10 by 10 with a file cap (called overcrowd in the parameters) of 200 is usually sufficient for variables like distance.

checkCells takes advantage of this counter system instead of inquiring for files individually. Because it needs the counters, these counters are written down onto text files at the end of getCells.f90 and initialTrajectories.f90. Although the creation of the grid is important, an accurate counter system is necessary to decrease computational load.
