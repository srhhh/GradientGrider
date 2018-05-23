
module f1_parameters
implicit none
!integer :: Nvar,Nstates,Natoms,overcrowd
!real :: spacing1,spacing2,spacing3,scaling1,scaling2,scaling3,max_var1,max_var2,max_var3
!character(50) :: path1, path2
!character(20) :: file1



!Path to the files
character(50),parameter :: path1 = "/home/ruisun/proj/ruisun/B0/"
character(50),parameter :: path2 = "/home/kazuumi/lus/B0/"
character(50),parameter :: path3 = "/home/kazuumi/lus/B0/grid/"

!File that will keep the trajectory folder names
character(20),parameter :: file1 = "f1_trajectories.txt"

!Some parameters, Nlength is the initial guess for how many states; this can
!grow later
integer,parameter :: Natoms = 6
!integer,parameter :: Nlength = 10000
integer,parameter :: Nvar = 3

!The spacing is the spacing for the gridlines; currently var3 is not grided
real,parameter :: spacing1 = 1.0
real,parameter :: spacing2 = 1.0
real,parameter :: spacing3 = 0.1

!The number of states that we believe "overcrowds" the square
integer,parameter :: overcrowd = 10000

!The scaling is the amount that is resolved for a variable when there is
!overcrowding (10,10,10) = x1000 magnification
integer,parameter :: scaling1 = 10
integer,parameter :: scaling2 = 10
integer,parameter :: scaling3 = 10

!There are some outliers; making a maximum throws these away
real,parameter :: max_var1 = 80.0
real,parameter :: max_var2 = 80.0



end module f1_parameters



