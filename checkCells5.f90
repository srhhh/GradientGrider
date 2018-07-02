
module checkCells5
implicit none


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              FUNCTION CHECKSTATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!INPUT: real,dim(3*Natoms) coords                       "the to-be-checked frame"
!       integer, dim(...) counter'X'                    "input counters so as to not re-read everytime"
!OUTPUT real, dim(6*Natoms) closestCoords               "closest frame+gradient"
!       dp min_rmsd                                     "closest frame rmsd"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine checkState(coords,closestCoords,min_rmsd_prime,force_Neighbors,&
                      path_to_directory,Ngrid_max,filechannels,&
                      number_of_frames,order,neighbor_check)
use f2_variables
use f2_parameters
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
real, dimension(Ncoords), intent(in) :: coords
real, dimension(6*Natoms), intent(out) :: closestCoords
real, dimension(6*Natoms) :: candidateCoords
double precision, dimension(3,Natoms) :: rmsd_coords1,rmsd_coords2
double precision, dimension(3) :: x_center,y_center
double precision, intent(out) :: min_rmsd_prime
double precision  :: candidate_rmsd
double precision, allocatable :: neighbor_rmsds(:)
real, allocatable :: neighbor_coords(:,:)
double precision, allocatable :: U(:,:), g(:,:)
logical :: subcell_existence
character(100) :: subcell
character(FMTlength) ::  var1_filename0, var2_filename0, var1_filename1, var2_filename1
character(*),intent(in) :: path_to_directory
character(len(path_to_directory)+9) :: path_to_grid
integer :: Ngrid
integer, intent(in) :: Ngrid_max
character(3) :: Ngrid_text
double precision :: min_rmsd
integer, dimension(Ngrid_max),intent(in) :: filechannels
integer :: subcell0search_max,subcell1search_max

subcell0search_max = 1
subcell1search_max = 1

number_of_frames = 0
neighbor_check = 0

call getVar3(coords,Natoms,var1)
call getVar4(coords,Natoms,var2)
rmsd_coords1 = reshape(coords,(/3, Natoms/))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 ORDER 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!var_cell keeps track of the child subcell portion of the decimal
!and throws away the parent subcell portion of the decimal
var1_cell = var1
var2_cell = var2

!var_index keeps track of the index of the child subcell
!inside of its parent subcell
var1_index0 = int(var1_cell * divisor1_0)
var2_index0 = int(var2_cell * divisor2_0)

!var_round is the variable rounded to the subcell length multiple
!This may have more or even less digits than required
var1_round0 = multiplier1_0 * var1_index0
var2_round0 = multiplier2_0 * var2_index0

var1_cell = var1_cell - var1_round0
var2_cell = var2_cell - var2_round0

var1_index1 = int(var1_cell * divisor1_1)
var2_index1 = int(var2_cell * divisor2_1)

var1_round1 = multiplier1_1 * var1_index1
var2_round1 = multiplier2_1 * var2_index1

var1_round = var1_round0 + var1_round1
var2_round = var2_round0 + var2_round1

write(var1_filename1,FMT=FMTorder1) var1_round
write(var2_filename1,FMT=FMTorder1) var2_round

write(var1_filename0,FMT=FMTorder0) var1_round0
write(var2_filename0,FMT=FMTorder0) var2_round0

!Do this for every grid
min_rmsd_prime = 100.0
do Ngrid = 1, Ngrid_max

min_rmsd = 100.0
write(Ngrid_text,FMT="(I0.3)") Ngrid
path_to_grid = path_to_directory//Ngrid_text//"/grid/"
subcell = path_to_grid//trim(adjustl(var1_filename1))//"_"//trim(adjustl(var2_filename1))

allocate(neighbor_rmsds(overcrowd1))
allocate(neighbor_coords(overcrowd1,Nvar+6*Natoms))

inquire(file=trim(subcell)//".dat",exist=subcell_existence)

if (.false.) then
print *, ""
print *, ""
print *, "inquiring for: ", subcell
print *, "exists? ", subcell_existence
print *, ""
end if

if (subcell_existence) then
        !getRMSD reads off the coordinates and calculated the RMSD with
        !ls_rmsd module
        call getRMSD(subcell,overcrowd1,rmsd_coords1,population,neighbor_rmsds,neighbor_coords)
        
        !Using minloc locates the position in the array
        !with the lowest value rmsd
        min_position = minloc(neighbor_rmsds(1:population),1)
        candidate_rmsd = neighbor_rmsds(min_position)
        if (candidate_rmsd < min_rmsd) then
        	min_rmsd = candidate_rmsd
        	!The frame with the closest coordinates has this position
        	closestCoords = neighbor_coords(min_position,Nvar+1:Nvar+6*Natoms)
        end if

	number_of_frames = population
	if (.not.(force_Neighbors)) then
		deallocate(neighbor_rmsds)
		deallocate(neighbor_coords)
		if (min_rmsd < min_rmsd_prime) min_rmsd_prime = min_rmsd
		write(filechannels(Ngrid),FMT="(F12.8)") min_rmsd_prime
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
                                  rmsd_coords1,candidate_rmsd,candidateCoords,number_of_frames)

                !Even if all neighbors are empty, it still returns a
                !minimum rmsd (default is 100.0)
                !ANY frame is better than none so it stops after it
                !finds one
                if (candidate_rmsd < min_rmsd) then
                        min_rmsd = candidate_rmsd
                        closestCoords = candidateCoords
                        stop_flag = .true.
                end if
        end do

        !We don't stop looking immediately (we may find a better fit
        !somewhere else on the circumfrence) but we know we don't need
        !to look any farther because farther subcells will normally have
        !a larger RMSD
        if (stop_flag) then
		deallocate(neighbor_rmsds)
		deallocate(neighbor_coords)
		if (min_rmsd < min_rmsd_prime) min_rmsd_prime = min_rmsd
		write(filechannels(Ngrid),FMT="(F12.8)") min_rmsd_prime
		exit
	end if
        end do

	if (stop_flag) cycle
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 ORDER 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

deallocate(neighbor_rmsds)
deallocate(neighbor_coords)

allocate(neighbor_rmsds(overcrowd0))
allocate(neighbor_coords(overcrowd0,Nvar+6*Natoms))

subcell = path_to_grid//trim(adjustl(var1_filename0))//"_"//trim(adjustl(var2_filename0))

inquire(file=trim(subcell)//".dat",exist=subcell_existence)

if (.false.) then
print *, ""
print *, ""
print *, "inquiring for: ", trim(subcell)//".dat" 
print *, "exists? ", subcell_existence
print *, ""
end if

if (subcell_existence) then
        !Get the RMSDs of the frames inside of it
        call getRMSD(subcell,overcrowd0,rmsd_coords1,population,&
                        neighbor_rmsds,neighbor_coords)

        !Obtain the rmsd and coordinates of the closest frame
        min_position = minloc(neighbor_rmsds(1:population),1)
        candidate_rmsd = neighbor_rmsds(min_position)
	if (candidate_rmsd < min_rmsd) then
		min_rmsd = candidate_rmsd
        	closestCoords = neighbor_coords(min_position,Nvar+1:Nvar+6*Natoms)
	end if

	number_of_frames = population
!	if (.not.(force_Neighbors))
end if

!If there are no frames in this cell, we can look at one of its
!neighbors (no higher level subcell)
if ((.not.(subcell_existence)).or.(force_Neighbors)) then
        neighbor_check = 1

        !Once we go far enough outside, we can stop
        stop_flag = .false.

        !Integer i keeps track of how far away from the original
        !subcell we are; we look at cells on the 'diamond' surrounding
        !the original subcell
        do i = 1, subcell0search_max

        !And integer j keeps track of where on the circumfrence
        !Of the diamond surrounding the subcell we are at
        do j = 0, i
        
                !Yes, there are redundancies for j = 0 and j = i
                !This is the first thing we can improve upon
                !(maybe with an if statement?)
                index1_1 = var1_index0 + j
                index1_2 = var1_index0 - j

                index2_1 = var2_index0 + i - j
                index2_2 = var2_index0 - i + j

                !This does the nitty-gritty of checking to see
                !If any of the indexes are out of bounds and 
                !obtaining their rmsds and coordinates
                call getNeighbors(bounds1,bounds2,multiplier1_0,multiplier2_0,FMTorder1,&
                                  index1_1,index1_2,index2_1,index2_2,&
                                  0.0,0.0,path_to_grid,overcrowd0,&
                                  rmsd_coords1,candidate_rmsd,candidateCoords,number_of_frames)

                !Even if all neighbors are empty, it still returns a
                !minimum rmsd (default is 100.0)
                !ANY frame is better than none so it stops after it
                !finds one
                if (candidate_rmsd < min_rmsd) then
                        min_rmsd = candidate_rmsd
                        closestCoords = candidateCoords
                        stop_flag = .true.
                end if
        end do

        !We don't stop looking immediately (we may find a better fit
        !somewhere else on the circumfrence) but we know we don't need
        !to look any farther because farther subcells will normally have
        !a larger RMSD
        if (stop_flag) exit
        end do

end if

deallocate(neighbor_rmsds)
deallocate(neighbor_coords)
if (min_rmsd < min_rmsd_prime) min_rmsd_prime = min_rmsd
write(filechannels(Ngrid),FMT="(F12.8)") min_rmsd_prime

!print *, "   no nearest frame for: ", trim(subcell)
end do

end subroutine checkState

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          FUNCTION GETNEIGHBORS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!INPUT: integer scaling1,scaling2                       "subcell spacing"
!       integer index_start                             "index of parent subcell"
!       integer index1_1,index1_2,index2_1,index2_2     "index of subcells"
!       strings integer1,integer2                       "integer portion"
!       integer decimals1,decimals2                     "decimal portion"
!       integer order                                   "order of subcell"
!       integer counterN_max                            "dimension of counter"
!       integer, dim(counterN_max) counterN             "counter"
!       dp, dim(3,Natoms) coords_static                 "comparison frame"
!OUTPUT dp min_rmsd                                     "min rmsd of subcells"
!       real, dim(6*Natoms) closestCoords               "coords of min rmsd frame"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getNeighbors(scaling1,scaling2,multiplier1,multiplier2,FMTorder,&
                        index1_1,index1_2,index2_1,index2_2,&
                        var1_round0,var2_round0,path_to_grid,overcrowd,&
                        coords_static,min_rmsd,closestCoords,number_of_frames)

use f2_parameters
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
double precision, intent(out) :: min_rmsd
double precision, dimension(overcrowd) :: rmsds
real, dimension(overcrowd,Nvar+6*Natoms) :: coords
integer :: min_position,population
real, dimension(6*Natoms), intent(out) :: closestCoords
double precision,intent(in), dimension(3,Natoms) :: coords_static

!Default to 100.0 rmsd
min_rmsd = 100.0

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

	inquire(file=trim(subcell)//".dat",exist=subcell_existence)
if (.false.) print *, "inquiring for: ", subcell
	if (subcell_existence) then

                !Calculate the RMSDs
                call getRMSD(subcell,overcrowd,coords_static,population,rmsds,coords)

                !Now, we want the state that is closest in terms of rmsd
                min_position = minloc(rmsds(1:population),1)
		if (rmsds(min_position) < min_rmsd) then
                min_rmsd = rmsds(min_position)
                closestCoords = coords(min_position,Nvar+1:Nvar+6*Natoms)
		end if
	end if

	number_of_frames = number_of_frames + population
end if

!Rinse and repeat
if ((flag1) .and. (flag4)) then

        !Make the name of the subcell
	var1_round = var1_round0 + index1_1 * multiplier1
	var2_round = var2_round0 + index2_2 * multiplier2
        write(var1_filename,FMT=FMTorder) var1_round
        write(var2_filename,FMT=FMTorder) var2_round
        subcell = path_to_grid//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

	inquire(file=trim(subcell)//".dat",exist=subcell_existence)
if (.false.) print *, "inquiring for: ", subcell
	if (subcell_existence) then

                !Calculate the RMSDs
                call getRMSD(subcell,overcrowd,coords_static,population,rmsds,coords)

                !Now, we want the state that is closest in terms of rmsd
                min_position = minloc(rmsds(1:population),1)
		if (rmsds(min_position) < min_rmsd) then
                min_rmsd = rmsds(min_position)
                closestCoords = coords(min_position,Nvar+1:Nvar+6*Natoms)
		end if
	end if

	number_of_frames = number_of_frames + population
end if

!Rinse and repeat
if ((flag2) .and. (flag3)) then

        !Make the name of the subcell
	var1_round = var1_round0 + index1_2 * multiplier1
	var2_round = var2_round0 + index2_1 * multiplier2
        write(var1_filename,FMT=FMTorder) var1_round
        write(var2_filename,FMT=FMTorder) var2_round
        subcell = path_to_grid//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

	inquire(file=trim(subcell)//".dat",exist=subcell_existence)
if (.false.) print *, "inquiring for: ", subcell
	if (subcell_existence) then

                !Calculate the RMSDs
                call getRMSD(subcell,overcrowd,coords_static,population,rmsds,coords)

                !Now, we want the state that is closest in terms of rmsd
                min_position = minloc(rmsds(1:population),1)
		if (rmsds(min_position) < min_rmsd) then
                min_rmsd = rmsds(min_position)
                closestCoords = coords(min_position,Nvar+1:Nvar+6*Natoms)
		end if
	end if

	number_of_frames = number_of_frames + population
end if

!Rinse and repeat
if ((flag2) .and. (flag4)) then

        !Make the name of the subcell
	var1_round = var1_round0 + index1_2 * multiplier1
	var2_round = var2_round0 + index2_2 * multiplier2
        write(var1_filename,FMT=FMTorder) var1_round
        write(var2_filename,FMT=FMTorder) var2_round
        subcell = path_to_grid//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

	inquire(file=trim(subcell)//".dat",exist=subcell_existence)
if (.false.) print *, "inquiring for: ", subcell
	if (subcell_existence) then

                !Calculate the RMSDs
                call getRMSD(subcell,overcrowd,coords_static,population,rmsds,coords)

                !Now, we want the state that is closest in terms of rmsd
                min_position = minloc(rmsds(1:population),1)
		if (rmsds(min_position) < min_rmsd) then
                min_rmsd = rmsds(min_position)
                closestCoords = coords(min_position,Nvar+1:Nvar+6*Natoms)
		end if
	end if

	number_of_frames = number_of_frames + population
end if




end subroutine getNeighbors



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               FUNCTION GETRMSD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  actually pretty self-explanatory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getRMSD(subcell,overcrowd,coords_static,population,rmsds,coords)

use ls_rmsd
use f2_parameters
implicit none
integer,intent(in) :: overcrowd
integer,intent(out) :: population
real,intent(out), dimension(overcrowd,Nvar+6*Natoms) :: coords
double precision,intent(out), dimension(overcrowd) :: rmsds
integer :: i,j,iostate
character(*),intent(in) :: subcell
double precision,intent(in), dimension(3,Natoms) :: coords_static
double precision, dimension(3,Natoms) :: rmsd_coords2
double precision, allocatable :: U(:,:), g(:,:)
double precision, dimension(3) :: x_center,y_center

!Read off the states in this subcell
open(filechannel1,file=trim(subcell)//".dat")
population = 1
do
        read(filechannel1,FMT=FMT1,advance="no",iostat=iostate) (coords(population,j),j=1,Nvar)
	if (iostate /= 0) exit
        read(filechannel1,FMT=FMT2) (coords(population,j),j=Nvar+1,Nvar+6*Natoms)

        !Need to make the coordinates readable for the rmsd
        rmsd_coords2 = reshape(coords(population,Nvar+1:Nvar+3*Natoms),&
                                (/3, Natoms/))
        call rmsd(Natoms,coords_static,rmsd_coords2,0,U,&
                  x_center,y_center,rmsds(population),.false.,g)
	population = population + 1
end do
population = population-1
close(filechannel1)


end subroutine getRMSD




end module checkCells5
