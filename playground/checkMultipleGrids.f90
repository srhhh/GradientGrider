module checkMultipleGrids
implicit none


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              FUNCTION CHECKSTATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!INPUT: real,dim(3*Natoms)	coords                  "the to-be-checked frame"
!	logical 		force_Neighbors		"do we force it to check adjacent cells?"
!
!	character(*) 		path_to_directory	"the path to the directory containing the grids"
!	integer 		Ngrid_total		"the number of grids to check"
!	integer,dim(Ngrid_total) 	filechannels		"the files we write to for data on each grid"
!
!OUTPUT real, dim(6*Natoms) 	approx_gradient           "closest frame+gradient"
!       dp 			min_rmsd                "closest frame rmsd"
!
!	integer number_of frames			"the number of frames checked" (deprecated) 
!		order					"the order of the subcell checked" (deprecated)
!		neighbor_check				"did we check adjacent cells?" (deprecated)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	This checkstate does not rely on counters
!	Instead it inquires for some subcell in each grid
!	Of course, this assumes each grid uses the same parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine checkState(vals,coords,approx_gradient,min_rmsd,force_Neighbors,&
                      path_to_directory,Ngrid_total,filechannels,&
                      number_of_frames,order,neighbor_check)
use ANALYSIS
use VARIABLES
use PARAMETERS
use mapCellData
implicit none
integer :: i,j,k
integer,intent(out),optional :: order,number_of_frames,neighbor_check
integer :: var1_index0,var2_index0,var1_index1,var2_index1
integer :: index1_1,index1_2,index2_1,index2_2
integer :: min_position,population
logical :: stop_flag,flag1,flag2,flag3,flag4
logical,intent(in) :: force_Neighbors
real :: var1,var2,var1_cell,var2_cell,var1_round,var2_round
real :: var1_round0,var2_round0,var1_round1,var2_round1,var1_round2,var2_round2,var1_round3,var2_round3
real(dp), dimension(Nvar), intent(in) :: vals
real(dp), dimension(3,Natoms), intent(in) :: coords
real(dp), dimension(3,Natoms), intent(out) :: approx_gradient
real(dp), dimension(3,Natoms) :: candidate_gradient
real(dp), allocatable :: candidate_gradients(:,:,:), U(:,:,:)
real(dp), intent(out) :: min_rmsd
real(dp)  :: candidate_rmsd
real(dp), allocatable :: rmsds(:)
logical :: subcell_existence
character(100) :: subcell
character(FMTlength) ::  var1_filename0, var2_filename0, var1_filename1, var2_filename1
character(*),intent(in) :: path_to_directory
character(len(path_to_directory)+9) :: path_to_grid
integer :: Ngrid
integer, intent(in) :: Ngrid_total
character(3) :: Ngrid_text
integer, dimension(Ngrid_total),intent(in) :: filechannels
integer :: subcell0search_max,subcell1search_max

!New feature
!In order to decrease the number of frames checked
!If we need to look at adjacent cells for a frame
!(for instance, if the target cell is empty)
!Then we stop looking after we are subcellsearch_max cells away
subcell0search_max = 1
subcell1search_max = 1

!These variables always default to zero
number_of_frames = 0
neighbor_check = 0
number_of_frames = 0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 SUBCELL TARGETING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!var_cell keeps track of the child subcell portion of the decimal
!and throws away the parent subcell portion of the decimal
var1_cell = vals(1)
var2_cell = vals(2)

!var_index keeps track of the index of the child subcell
!inside of its parent subcell
var1_index0 = int(var1_cell * divisor1_0)
var2_index0 = int(var2_cell * divisor2_0)

!var_round is the variable rounded to the subcell length multiple
!This may have more or even less digits than required
var1_round0 = multiplier1_0 * var1_index0
var2_round0 = multiplier2_0 * var2_index0

!Repeat for the child subcell portion
var1_cell = var1_cell - var1_round0
var2_cell = var2_cell - var2_round0

var1_index1 = int(var1_cell * divisor1_1)
var2_index1 = int(var2_cell * divisor2_1)

var1_round1 = multiplier1_1 * var1_index1
var2_round1 = multiplier2_1 * var2_index1

var1_round = var1_round0 + var1_round1
var2_round = var2_round0 + var2_round1

!The child-level subcells will use var_round
write(var1_filename1,FMT=FMTorder1) var1_round
write(var2_filename1,FMT=FMTorder1) var2_round

!The parent-level subcells will use var_round0
write(var1_filename0,FMT=FMTorder0) var1_round0
write(var2_filename0,FMT=FMTorder0) var2_round0

!We set min_rmsd to be 0.200100 by default
min_rmsd = 0.200100d0

!Now, we start iterating over the grids
do Ngrid = 1, Ngrid_total

!We need a path to the grid
write(Ngrid_text,FMT="(I0.3)") Ngrid
path_to_grid = path_to_directory//Ngrid_text//"/grid/"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 ORDER 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Construct the subcell filename
subcell = path_to_grid//trim(adjustl(var1_filename1))//"_"//trim(adjustl(var2_filename1))

!We don't know how many frames may be in the child subcell
!But we do know there will be no more than overcrowd1
allocate(rmsds(overcrowd1),candidate_gradients(3,Natoms,overcrowd1),U(3,3,overcrowd1))

!See whether this subcell exists
inquire(file=trim(subcell)//".dat",exist=subcell_existence)

!For bug-testing
if (.false.) then
print *, ""
print *, ""
print *, "inquiring for: ", subcell
print *, "exists? ", subcell_existence
print *, ""
end if

!If the subcell exists, we will look in here for a closeby frame
if (subcell_existence) then
        !getRMSD reads off the coordinates and calculated the RMSD with
        !ls_rmsd module
	call getRMSD_dp(subcell,overcrowd1,coords,population,rmsds,candidate_gradients,U)
        
        !Using minloc locates the position in the array
        !with the lowest value rmsd
        min_position = minloc(rmsds(1:population),1)
        candidate_rmsd = rmsds(min_position)
        if (candidate_rmsd < min_rmsd) then
        	min_rmsd = candidate_rmsd
		approx_gradient = matmul(U(:,:,min_position),candidate_gradients(:,:,min_position))
        end if

	!However many frames are in the subcell we increment to number_of_frames
	number_of_frames = number_of_frames + population
	if (.not.(force_Neighbors)) then
		deallocate(rmsds,candidate_gradients,U)

		!If min_rmsd is lower than min_rmsd_prime, then we have a new candidate
		!For lowest RMSD
		write(filechannels(Ngrid),FMT=FMT6) min_rmsd
		cycle
	end if 
end if

!If there are no frames in this cell, we can look at one of its
!neighbors (better than having to look at the parent)
if ((.not.(subcell_existence)).or.(force_Neighbors)) then
        neighbor_check = 1

        !Once we go far enough outside, we can stop
        stop_flag = .false.

        !Integer i keeps track of how far away from the original
        !subcell we are; we look at cells on the 'diamond' surrounding
        !the original subcell
        do i = 1, subcell1search_max

!For bug-testing
if (.false.) then
print *, "Looking for neighbors: ", i, " cells away"
print *, "Scalings: ", scaling1_0, scaling2_0
print *, "Multipliers: ", multiplier1_1, multiplier2_1
print *, "var: ", var1, var2
print *, "var_round0: ", var1_round0, var2_round0
print *, "var_round1: ", var1_round1, var2_round1
end if

        !And integer j keeps track of where on the circumfrence
        !Of the diamond surrounding the subcell we are at
        do j = 0, i
        
                !Yes, there are redundancies for j = 0 and j = i
                !This is the first thing we can improve upon
                !(maybe with an if statement?)
                index1_1 = var1_index1 + j
                index1_2 = var1_index1 - j

                index2_1 = var2_index1 + i - j
                index2_2 = var2_index1 - i + j

                !This does the nitty-gritty of checking to see
                !If any of the indexes are out of bounds and 
                !obtaining their rmsds and coordinates
                call getNeighbors(scaling1_0,scaling2_0,multiplier1_1,multiplier2_1,FMTorder1,&
                                  index1_1,index1_2,index2_1,index2_2,&
                                  var1_round0,var2_round0,path_to_grid,overcrowd1,&
                                  coords,candidate_rmsd,candidate_gradient,number_of_frames)

                !Even if all neighbors are empty, it still returns a
                !minimum rmsd (default is 0.200100)
                !ANY frame is better than none so it stops after it
                !finds one
                if (candidate_rmsd < min_rmsd) then
                        min_rmsd = candidate_rmsd
                        approx_gradient = candidate_gradient
                        stop_flag = .true.
                end if
        end do

        !We don't stop looking immediately (we may find a better fit
        !somewhere else on the circumfrence) but we know we don't need
        !to look any farther because farther subcells will normally have
        !a larger RMSD
        if (stop_flag) then
		deallocate(rmsds,candidate_gradients,U)
		write(filechannels(Ngrid),FMT=FMT6) min_rmsd
		exit
	end if
        end do

	if (stop_flag) cycle
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 ORDER 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!We can use a smaller array size for parent-level cells
deallocate(rmsds,candidate_gradients,U)

allocate(rmsds(overcrowd0),candidate_gradients(3,Natoms,overcrowd0),U(3,3,overcrowd0))

!Rinse and repeat
subcell = path_to_grid//trim(adjustl(var1_filename0))//"_"//trim(adjustl(var2_filename0))

inquire(file=trim(subcell)//".dat",exist=subcell_existence)

!For bug-testing
if (.false.) then
print *, ""
print *, ""
print *, "inquiring for: ", trim(subcell)//".dat" 
print *, "exists? ", subcell_existence
print *, ""
end if

if (subcell_existence) then
	call getRMSD_dp(subcell,overcrowd0,coords,population,rmsds,candidate_gradients,U)

        min_position = minloc(rmsds(1:population),1)
        candidate_rmsd = rmsds(min_position)
	if (candidate_rmsd < min_rmsd) then
		min_rmsd = candidate_rmsd
		approx_gradient = matmul(U(:,:,min_position),candidate_gradients(:,:,min_position))
	end if

	number_of_frames = number_of_frames + population
end if

if ((.not.(subcell_existence)).or.(force_Neighbors)) then
        neighbor_check = 1

        stop_flag = .false.

        do i = 1, subcell0search_max

        do j = 0, i
        
                index1_1 = var1_index0 + j
                index1_2 = var1_index0 - j

                index2_1 = var2_index0 + i - j
                index2_2 = var2_index0 - i + j

                call getNeighbors(bounds1,bounds2,multiplier1_0,multiplier2_0,FMTorder1,&
                                  index1_1,index1_2,index2_1,index2_2,&
                                  0.0,0.0,path_to_grid,overcrowd0,&
                                  coords,candidate_rmsd,candidate_gradient,number_of_frames)

                if (candidate_rmsd < min_rmsd) then
                        min_rmsd = candidate_rmsd
                        approx_gradient = candidate_gradient
                        stop_flag = .true.
                end if
        end do

        if (stop_flag) then
		exit
	end if
        end do

end if

!To end, just deallocate and write the lowest RMSD frame
deallocate(rmsds,candidate_gradients,U)
write(filechannels(Ngrid),FMT=FMT6) min_rmsd

end do

end subroutine checkState

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          FUNCTION GETNEIGHBORS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!INPUT: integer		scaling1,scaling2                       "the number of subcells per parent cell"
!       real 		multiplier1,multiplier2                 "the spacing of the subcell"
!
!       integer 	index1_1,index1_2,index2_1,index2_2     "index of subcells to check"
!       real 		var1_round0, var2_round0		"the value of the parent subcell"
!
!       character(*) 		path_to_grid			"the path to the grid"
!       integer 		overcrowd			"the maximum number of frames we'll encounter in a subcell"
!       dp, dim(3,Natoms) 	coords_static                 	"comparison frame"
!
!OUTPUT dp			min_rmsd             		"min rmsd of subcells"
!       real, dim(6*Natoms) 	approx_gradient               	"coords of min rmsd frame"
!	integer 		number_of_frames		"keeps count of the number of frames"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getNeighbors(scaling1,scaling2,multiplier1,multiplier2,FMTorder,&
                        index1_1,index1_2,index2_1,index2_2,&
                        var1_round0,var2_round0,path_to_grid,overcrowd,&
                        coords,min_rmsd,approx_gradient,number_of_frames)
use ANALYSIS
use PARAMETERS
implicit none
logical :: flag1,flag2,flag3,flag4,subcell_existence
integer,intent(in) :: scaling1, scaling2
integer,intent(in) :: index1_1,index1_2,index2_1,index2_2
real,intent(in) :: multiplier1,multiplier2
real,intent(in) :: var1_round0,var2_round0
real :: var1_round,var2_round
!integer,intent(in) :: decimals1, decimals2
integer,intent(in) :: overcrowd
integer,intent(inout) :: number_of_frames
character(6),intent(in) :: FMTorder
character(FMTlength) :: var1_filename,var2_filename
character(100) :: subcell
character(*),intent(in) :: path_to_grid
real(dp),dimension(3,3,overcrowd) :: U
real(dp), intent(out) :: min_rmsd
real(dp), dimension(overcrowd) :: rmsds
real(dp), dimension(3,Natoms,overcrowd) :: candidate_gradients
integer :: min_position,population
real(dp), dimension(3,Natoms), intent(out) :: approx_gradient
real(dp),intent(in), dimension(3,Natoms) :: coords

!Default to 0.200100 rmsd
min_rmsd = 0.200100d0

!We need to check if any of the indexes (neighbors!) are out of bounds
flag1 = index1_1 < scaling1
flag2 = index1_2 .ge. 0
flag3 = index2_1 < scaling2
flag4 = index2_2 .ge. 0

!We check indexes in pairs (var1, var2)
if ((flag1) .and. (flag3)) then

        !Make the name of the subcell
	var1_round = var1_round0 + index1_1 * multiplier1
	var2_round = var2_round0 + index2_1 * multiplier2
        write(var1_filename,FMT=FMTorder) var1_round
        write(var2_filename,FMT=FMTorder) var2_round
        subcell = path_to_grid//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

	!Check if it exists
	inquire(file=trim(subcell)//".dat",exist=subcell_existence)

!For bug-testing
if (.false.) print *, "inquiring for: ", subcell

	!If it exists, get the rmsd of each frame in it
	if (subcell_existence) then

                !Calculate the RMSDs
		call getRMSD_dp(subcell,overcrowd,coords,population,rmsds,candidate_gradients,U)

                !Now, we want the state that is closest in terms of rmsd
                min_position = minloc(rmsds(1:population),1)

		!If it is lower than our original rmsd, we take it
		if (rmsds(min_position) < min_rmsd) then
                min_rmsd = rmsds(min_position)
		approx_gradient = matmul(U(:,:,min_position),candidate_gradients(:,:,min_position))
		end if
	end if

	!Increment however many frames we have
	number_of_frames = number_of_frames + population
end if

!Rinse and repeat
if ((flag1) .and. (flag4)) then

	var1_round = var1_round0 + index1_1 * multiplier1
	var2_round = var2_round0 + index2_2 * multiplier2
        write(var1_filename,FMT=FMTorder) var1_round
        write(var2_filename,FMT=FMTorder) var2_round
        subcell = path_to_grid//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

	inquire(file=trim(subcell)//".dat",exist=subcell_existence)
if (.false.) print *, "inquiring for: ", subcell
	if (subcell_existence) then

		call getRMSD_dp(subcell,overcrowd,coords,population,rmsds,candidate_gradients,U)

                min_position = minloc(rmsds(1:population),1)
		if (rmsds(min_position) < min_rmsd) then
                min_rmsd = rmsds(min_position)
		approx_gradient = matmul(U(:,:,min_position),candidate_gradients(:,:,min_position))
		end if
	end if

	number_of_frames = number_of_frames + population
end if

!Rinse and repeat
if ((flag2) .and. (flag3)) then

	var1_round = var1_round0 + index1_2 * multiplier1
	var2_round = var2_round0 + index2_1 * multiplier2
        write(var1_filename,FMT=FMTorder) var1_round
        write(var2_filename,FMT=FMTorder) var2_round
        subcell = path_to_grid//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

	inquire(file=trim(subcell)//".dat",exist=subcell_existence)
if (.false.) print *, "inquiring for: ", subcell
	if (subcell_existence) then

		call getRMSD_dp(subcell,overcrowd,coords,population,rmsds,candidate_gradients,U)

                min_position = minloc(rmsds(1:population),1)
		if (rmsds(min_position) < min_rmsd) then
                min_rmsd = rmsds(min_position)
		approx_gradient = matmul(U(:,:,min_position),candidate_gradients(:,:,min_position))
		end if
	end if

	number_of_frames = number_of_frames + population
end if

!Rinse and repeat
if ((flag2) .and. (flag4)) then

	var1_round = var1_round0 + index1_2 * multiplier1
	var2_round = var2_round0 + index2_2 * multiplier2
        write(var1_filename,FMT=FMTorder) var1_round
        write(var2_filename,FMT=FMTorder) var2_round
        subcell = path_to_grid//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

	inquire(file=trim(subcell)//".dat",exist=subcell_existence)
if (.false.) print *, "inquiring for: ", subcell
	if (subcell_existence) then

		call getRMSD_dp(subcell,overcrowd,coords,population,rmsds,candidate_gradients,U)

                min_position = minloc(rmsds(1:population),1)
		if (rmsds(min_position) < min_rmsd) then
                min_rmsd = rmsds(min_position)
		approx_gradient = matmul(U(:,:,min_position),candidate_gradients(:,:,min_position))
		end if
	end if

	number_of_frames = number_of_frames + population
end if



end subroutine getNeighbors



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               FUNCTION GETRMSD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  actually pretty self-explanatory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getRMSD(subcell,overcrowd,coords,population,rmsds,candidate_gradients,U)

use ls_rmsd
use PARAMETERS
implicit none
integer,intent(in) :: overcrowd
integer,intent(out) :: population
real,intent(out), dimension(3,Natoms,overcrowd) :: candidate_gradients
real,intent(out), dimension(overcrowd) :: rmsds
integer :: i,j,iostate
character(*),intent(in) :: subcell
real,intent(in), dimension(3,Natoms) :: coords
real, dimension(3,Natoms) :: coords2
real,dimension(3,3,overcrowd),intent(out) :: U

!Read off the states in this subcell
open(filechannel1,file=trim(subcell)//".dat")
population = 1
do
	!Once we've had enough, exit
        read(filechannel1,FMT=FMT7,advance="no",iostat=iostate) ((coords2(i,j),i=1,3),j=1,Natoms)
	if (iostate /= 0) exit
        read(filechannel1,FMT=FMT3) ((candidate_gradients(i,j,population),i=1,3),j=1,Natoms)

	!Get the RMSD
        call rmsd(Natoms,coords2,coords,rmsds(population),U(:,:,population))

	!Increment the position
	population = population + 1
end do

!The population is one less than the current position
population = population - 1
close(filechannel1)


end subroutine getRMSD



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               FUNCTION GETRMSD_DP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  actually pretty self-explanatory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getRMSD_dp(subcell,overcrowd,coords,population,rmsds,candidate_gradients,U)

use ls_rmsd_original
use PARAMETERS
implicit none
integer,intent(in) :: overcrowd
integer,intent(out) :: population
real(dp),intent(out), dimension(3,Natoms,overcrowd) :: candidate_gradients
real(dp),intent(out), dimension(overcrowd) :: rmsds
integer :: i,j,iostate
character(*),intent(in) :: subcell
real(dp),intent(in), dimension(3,Natoms) :: coords
real(dp) :: min_rmsd_dp
real(dp), dimension(3) :: x_center, y_center
real(dp), allocatable :: g(:,:)
real(dp),dimension(3,3,overcrowd),intent(out) :: U
real(dp), dimension(3,Natoms) :: coords2

!Read off the states in this subcell
open(filechannel1,file=trim(subcell)//".dat")
population = 1
do
	!Once we've had enough, exit
        read(filechannel1,FMT=FMT7,advance="no",iostat=iostate) ((coords2(i,j),i=1,3),j=1,Natoms)
	if (iostate /= 0) exit
        read(filechannel1,FMT=FMT3) ((candidate_gradients(i,j,population),i=1,3),j=1,Natoms)

        !Get the RMSD
        call rmsd_dp(Natoms,coords2,coords,1,U(:,:,population),x_center,y_center,&
                  rmsds(population))!,.false.,g)

        !Increment the position
        population = population + 1
end do

!The population is one less than the current position
population = population - 1
close(filechannel1)

end subroutine getRMSD_dp




end module checkMultipleGrids