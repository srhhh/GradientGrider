
module checkCells4
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

subroutine checkState(coords,closestCoords,min_rmsd,&
                      counter0,counter1,counter2,counter3)
use f2_variables
use f2_parameters
use mapCellData
implicit none
integer :: i,j,k
integer :: var1_index,var2_index,indexer,index1_1,index1_2,index2_1,index2_2
integer :: key0,key1,key2,key3
integer :: min_position
integer, dimension(counter0_max) :: counter0
integer, dimension(counter1_max) :: counter1
integer, dimension(counter2_max) :: counter2
integer, dimension(counter3_max) :: counter3
logical :: stop_flag,flag1,flag2,flag3,flag4
real :: var1,var2,var1_cell,var2_cell,var1_round,var2_round
real :: var1_round0,var2_round0,var1_round1,var2_round1,var1_round2,var2_round2,var1_round3,var2_round3
real, dimension(Ncoords), intent(in) :: coords
real, dimension(6*Natoms), intent(out) :: closestCoords
real, dimension(6*Natoms) :: candidateCoords
double precision, dimension(3,Natoms) :: rmsd_coords1,rmsd_coords2
double precision, dimension(3) :: x_center,y_center
double precision, intent(out) :: min_rmsd
double precision  :: candidate_rmsd
double precision, allocatable :: neighbor_rmsds(:)
real, allocatable :: neighbor_coords(:,:)
double precision, allocatable :: U(:,:), g(:,:)
character(50) :: subcell
character(9) ::  var1_filename, var2_filename

! Need a ridiculously large number
min_rmsd = 100.0

! Get the variables corresponding to frame
! RS: What do you mean by 'each'? Isn't the input only one frame?
! RS: calculating var again so you don't have to pass them?
!                       KF: the word 'each' is misleading, my bad
call getVar3(coords,Natoms,var1)
call getVar4(coords,Natoms,var2)

!The coordinates, as they are formatted in getCells and addCells, are the wrong
!shape for ls_rmsd. Thus, we reshape them first

! RS: ha! Should have read them in this format at the first place
!                       KF: ohhh maaaannnnnnn
rmsd_coords1 = reshape(coords,(/3, Natoms/))

! RS: I have some thoughts on the following -- Let's talk tomorrow


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 ORDER 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!var_cell keeps track of the child subcell portion of the decimal
!and throws away the parent subcell portion of the decimal
var1_cell = var1
var2_cell = var2

!var_index keeps track of the index of the child subcell
!inside of its parent subcell
var1_index = int(var1_cell * divisor1_0)
var2_index = int(var2_cell * divisor2_0)

!var_round is the variable rounded to the subcell length multiple
!This may have more or even less digits than required
var1_round0 = multiplier1_0 * var1_index
var2_round0 = multiplier2_0 * var2_index

!indexer accesses the counter (so as to keep track of population) and is uniquely determined by var_index
indexer = bounds1 * var2_index + var1_index + 1

!The population will be the digits on the right of the key
!The 'key' or index to the next counter are the remaining digits
key0 = counter0(indexer)

!Write to the progress file for bug-testing
!open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
!write(progresschannel,*) "Investigating key ", counter0(indexer), " with index ", indexer, " and order 0"
!write(progresschannel,*) ""
!close(progresschannel)

!If the key is zero, then that means divyUp was not called on it
!So there are no children subcells, so this subcell must be examined
if (key0 < overcrowd0) then

                !If the population of the parent cell is empty, we should just
                !give up really
                if (key0 == 0) return

                !Call mapCell to see the heat map of all parent cells
!               call mapCell(0,counter0,counter0_max,bounds1,bounds2)

                !Make the name of the subcell
		var1_round = var1_round0
		var2_round = var2_round0
                write(var1_filename,FMT=FMTorder0) var1_round
                write(var2_filename,FMT=FMTorder0) var2_round
                subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

                ! getRMSD reads off the coordinates and calculated the RMSD with
                !ls_rmsd module
                allocate(neighbor_rmsds(key0),&
                        neighbor_coords(key0,Nvar+6*Natoms))
                call getRMSD(subcell,key0,rmsd_coords1,&
                                neighbor_rmsds,neighbor_coords)

                !Using minloc locates the position in the array
                !with the lowest value rmsd
                min_position = minloc(neighbor_rmsds,1)
                min_rmsd = neighbor_rmsds(min_position)

                !The frame with the closest coordinates has this position
                closestCoords = neighbor_coords(min_position,Nvar+1:Nvar+6*Natoms)

                deallocate(neighbor_rmsds,neighbor_coords)
                return
end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 ORDER 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

var1_cell = var1_cell - var1_round0
var2_cell = var2_cell - var2_round0

var1_index = int(var1_cell * divisor1_1)
var2_index = int(var2_cell * divisor2_1)

var1_round1 = multiplier1_1 * var1_index
var2_round1 = multiplier2_1 * var2_index

! We need int(key/key_start) = the unique header1 value to access counter;
! And by multiplying by resolution (scaling1*scaling2), we assure that each
! subcell of the parent subcell gets its own unique index
indexer = resolution_0*int(key0/key_start-1) + scaling1_0*var2_index + var1_index + 1
key1 = counter1(indexer)

!open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
!write(progresschannel,*) "Investigating key ", counter1(indexer), " with index ", indexer, " and order 1"
!write(progresschannel,*) ""
!close(progresschannel)

if (key1 < overcrowd1) then

        !Call mapCell just to see the heatmap of this subcell
!       call mapCell(int(key0/key_start)-1,counter1,counter1_max,scaling1_0,scaling2_0)

        !If there are frames in this cell, then we retrieve those frames
        if (key1 > 0) then

                !Make the name of the subcell
		var1_round = var1_round0 + var1_round1
		var2_round = var2_round0 + var2_round1
                write(var1_filename,FMT=FMTorder1) var1_round
                write(var2_filename,FMT=FMTorder1) var2_round
                subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

                !Get the RMSDs of the frames inside of it
                allocate(neighbor_rmsds(key1),&
                        neighbor_coords(key1,Nvar+6*Natoms))
                call getRMSD(subcell,key1,rmsd_coords1,&
                                neighbor_rmsds,neighbor_coords)

                !Obtain the rmsd and coordinates of the closest frame
                min_position = minloc(neighbor_rmsds,1)
                min_rmsd = neighbor_rmsds(min_position)
                closestCoords = neighbor_coords(min_position,Nvar+1:Nvar+6*Natoms)

                deallocate(neighbor_rmsds,neighbor_coords)

        !If there are no frames in this cell, we can look at one of its
        !neighbors (better than having to look at the parent)
        else

		!We need to know the filename of the parent to get the filename of a child                
		var1_round = var1_round0
		var2_round = var2_round0

                !Once we go far enough outside, we can stop
                stop_flag = .false.

                !Integer i keeps track of how far away from the original
                !subcell we are; we look at cells on the 'diamond' surrounding
                !the original subcell
                do i = 1, scaling1_0

                !And integer j keeps track of where on the circumfrence
                !Of the diamond surrounding the subcell we are at
                do j = 0, i
                
                        !Yes, there are redundancies for j = 0 and j = i
                        !This is the first thing we can improve upon
                        !(maybe with an if statement?)
                        index1_1 = var1_index + j
                        index1_2 = var1_index - j

                        index2_1 = var2_index + i - j
                        index2_2 = var2_index - i + j

                        !This does the nitty-gritty of checking to see
                        !If any of the indexes are out of bounds and 
                        !obtaining their rmsds and coordinates
                        call getNeighbors(scaling1_0,scaling2_0,multiplier1_1,multiplier2_1,FMTorder1,&
                                          resolution_0*(int(key0/key_start)-1),index1_1,index1_2,index2_1,index2_2,&
                                          var1_round,var2_round,counter1,counter1_max,&
                                          rmsd_coords1,candidate_rmsd,candidateCoords)

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

        return
end if





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 ORDER 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

var1_cell = var1_cell - var1_round1
var2_cell = var2_cell - var2_round1

var1_index = int(var1_cell * divisor1_2)
var2_index = int(var2_cell * divisor2_2)

var1_round2 = multiplier1_2 * var1_index
var2_round2 = multiplier2_2 * var2_index

indexer = resolution_1*int(key1/key_start-1) + scaling1_1*var2_index + var1_index + 1
key2 = counter2(indexer)

!open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
!write(progresschannel,*) "Investigating key ", counter2(indexer), " with index ", indexer, " and order 2"
!write(progresschannel,*) ""
!close(progresschannel)

if (key2 < overcrowd2) then

!       call mapCell(key1/key_start-1,counter2,counter2_max,scaling1_1,scaling2_1)

!Force this to look at its neighbors (for testing)
!Remark: with normals parameters, the deepest subcell for a particular frame is
!usually of order 2 so this 

!       if (.false.) then !(key2 > 0) then
	if (key2 > 0) then

                !Make the name of the subcell
		var1_round = var1_round0 + var1_round1 + var1_round2
		var2_round = var2_round0 + var2_round1 + var2_round2
                write(var1_filename,FMT=FMTorder2) var1_round
                write(var2_filename,FMT=FMTorder2) var2_round
                subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

                allocate(neighbor_rmsds(key2),&
                        neighbor_coords(key2,Nvar+6*Natoms))
                call getRMSD(subcell,key2,rmsd_coords1,&
                                neighbor_rmsds,neighbor_coords)

                min_position = minloc(neighbor_rmsds,1)
                min_rmsd = neighbor_rmsds(min_position)
                closestCoords = neighbor_coords(min_position,Nvar+1:Nvar+6*Natoms)

                deallocate(neighbor_rmsds,neighbor_coords)
        else

!For bug-testing
!open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
!write(progresschannel,*) "Fetching neighbors..."
!close(progresschannel)

		var1_round = var1_round0 + var1_round1
		var2_round = var2_round0 + var2_round1

                stop_flag = .false.
                do i = 1, scaling1_1
                do j = 0, i                

                        index1_1 = var1_index + j
                        index1_2 = var1_index - j

                        index2_1 = var2_index + i - j
                        index2_2 = var2_index - i + j

                        call getNeighbors(scaling1_1,scaling2_1,multiplier1_2,multiplier2_2,FMTorder2,&
                                          resolution_1*(int(key1/key_start)-1),index1_1,index1_2,index2_1,index2_2,&
                                          var1_round,var2_round,counter2,counter2_max,&
                                          rmsd_coords1,candidate_rmsd,candidateCoords)

                        if (candidate_rmsd < min_rmsd) then
                                min_rmsd = candidate_rmsd
                                closestCoords = candidateCoords
                                stop_flag = .true.
                        end if

                end do
                if (stop_flag) exit
                end do

        end if

        return
end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 ORDER 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

var1_cell = var1_cell - var1_round2
var2_cell = var2_cell - var2_round2

var1_index = int(var1_cell * divisor1_3)
var2_index = int(var2_cell * divisor2_3)

var1_round3 = multiplier1_3 * var1_index
var2_round3 = multiplier2_3 * var2_index

indexer = resolution_2*int(key2/key_start-1) + scaling1_2*var2_index + var1_index + 1
key3 = counter3(indexer)

!open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
!write(progresschannel,*) "Investigating key ", counter3(indexer), " with index ", indexer, " and order 3"
!write(progresschannel,*) ""
!close(progresschannel)

if (key3 < overcrowd3) then

        if (key3 > 0) then

                !Make the name of the subcell
		var1_round = var1_round0 + var1_round1 + var1_round2 + var1_round3
		var2_round = var2_round0 + var2_round1 + var2_round2 + var2_round3
                write(var1_filename,FMT=FMTorder3) var1_round
                write(var2_filename,FMT=FMTorder3) var2_round
                subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

                allocate(neighbor_rmsds(key3),&
                        neighbor_coords(key3,Nvar+6*Natoms))
                call getRMSD(subcell,key3,rmsd_coords1,&
                                neighbor_rmsds,neighbor_coords)

                min_position = minloc(neighbor_rmsds,1)
                min_rmsd = neighbor_rmsds(min_position)
                closestCoords = neighbor_coords(min_position,Nvar+1:Nvar+6*Natoms)

                deallocate(neighbor_rmsds,neighbor_coords)
        else

		var1_round = var1_round0 + var1_round1 + var1_round2
		var2_round = var2_round0 + var2_round1 + var2_round2

                stop_flag = .false.
                do i = 1, scaling1_2
                do j = 0, i
                
                        index1_1 = var1_index + j
                        index1_2 = var1_index - j

                        index2_1 = var2_index + i - j
                        index2_2 = var2_index - i + j

                        call getNeighbors(scaling1_2,scaling2_2,multiplier1_3,multiplier2_3,FMTorder3,&
                                          resolution_2*(int(key2/key_start)-1),index1_1,index1_2,index2_1,index2_2,&
                                          var1_round,var2_round,counter3,counter3_max,&
                                          rmsd_coords1,candidate_rmsd,candidateCoords)

                        if (candidate_rmsd < min_rmsd) then
                                min_rmsd = candidate_rmsd
                                closestCoords = candidateCoords
                                stop_flag = .true.
                        end if
                end do
                if (stop_flag) exit
                end do

        end if

        return
end if

open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
write(progresschannel,*) ""
write(progresschannel,*) ""
write(progresschannel,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11"
write(progresschannel,*) "THE GRID IS NOT GRANULAR ENOUGH; A THIRD LEVEL SUBCELL IS OVERCROWDED!"
write(progresschannel,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11"
write(progresschannel,*) ""
write(progresschannel,*) ""
close(progresschannel)


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
                        index_start,index1_1,index1_2,index2_1,index2_2,&
                        var1_round0,var2_round0,counterN,counterN_max,&
                        coords_static,min_rmsd,closestCoords)

use f2_parameters
implicit none
logical :: flag1,flag2,flag3,flag4
integer,intent(in) :: scaling1, scaling2
integer,intent(in) :: index1_1,index1_2,index2_1,index2_2
real,intent(in) :: multiplier1,multiplier2
real,intent(in) :: var1_round0,var2_round0
real :: var1_round,var2_round
!integer,intent(in) :: decimals1, decimals2
integer,intent(in) :: index_start
integer :: indexer,population,key
integer,intent(in) :: counterN_max
integer,intent(in), dimension(counterN_max) :: counterN
character(6),intent(in) :: FMTorder
character(9) :: var1_filename,var2_filename
character(50) :: subcell
double precision, intent(out) :: min_rmsd
double precision, allocatable :: rmsds(:)
real, allocatable :: coords(:,:)
integer :: min_position
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

        !Calculate the index in counterN, just as in addCells.f90
        indexer = scaling1*(index2_1) + (index1_1) + 1
        key = counterN(index_start+indexer)
	population = modulo(key,key_start)

        !Only read off coordinates if there are any
        if (population > 0) then

                !Make the name of the subcell
		var1_round = var1_round0 + index1_1 * multiplier1
		var2_round = var2_round0 + index2_1 * multiplier2
                write(var1_filename,FMT=FMTorder) var1_round
                write(var2_filename,FMT=FMTorder) var2_round
                subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

                !Calculate the RMSDs
                allocate(rmsds(population),coords(population,Nvar+6*Natoms))
                call getRMSD(subcell,population,coords_static,rmsds,coords)

                !Now, we want the state that is closest in terms of rmsd
                min_position = minloc(rmsds,1)
                min_rmsd = rmsds(min_position)
                closestCoords = coords(min_position,Nvar+1:Nvar+6*Natoms)

                deallocate(rmsds,coords)
        end if
end if

!Rinse and repeat
if ((flag2) .and. (flag3)) then

        indexer = scaling1*(index2_1) + (index1_2) + 1
        key = counterN(index_start+indexer)
	population = modulo(key,key_start)

        if (population > 0) then

                !Make the name of the subcell
		var1_round = var1_round0 + index1_2 * multiplier1
		var2_round = var2_round0 + index2_1 * multiplier2
                write(var1_filename,FMT=FMTorder) var1_round
                write(var2_filename,FMT=FMTorder) var2_round
                subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

                allocate(rmsds(population),coords(population,Nvar+6*Natoms))
                call getRMSD(subcell,population,coords_static,rmsds,coords)

                min_position = minloc(rmsds,1)
                min_rmsd = rmsds(min_position)
                closestCoords = coords(min_position,Nvar+1:Nvar+6*Natoms)

                deallocate(rmsds,coords)
        end if
end if

if ((flag1) .and. (flag4)) then

        indexer = scaling1*(index2_2) + (index1_1) + 1
        key = counterN(index_start+indexer)
	population = modulo(key,key_start)

        if (population > 0) then

                !Make the name of the subcell
		var1_round = var1_round0 + index1_1 * multiplier1
		var2_round = var2_round0 + index2_2 * multiplier2
                write(var1_filename,FMT=FMTorder) var1_round
                write(var2_filename,FMT=FMTorder) var2_round
                subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))


                allocate(rmsds(population),coords(population,Nvar+6*Natoms))
                call getRMSD(subcell,population,coords_static,rmsds,coords)

                min_position = minloc(rmsds,1)
                min_rmsd = rmsds(min_position)
                closestCoords = coords(min_position,Nvar+1:Nvar+6*Natoms)

                deallocate(rmsds,coords)
        end if
end if

if ((flag2) .and. (flag4)) then

        indexer = scaling1*(index2_2) + (index1_2) + 1
        key = counterN(index_start+indexer)
	population = modulo(key,key_start)

        if (population > 0) then

                !Make the name of the subcell
		var1_round = var1_round0 + index1_2 * multiplier1
		var2_round = var2_round0 + index2_2 * multiplier2
                write(var1_filename,FMT=FMTorder) var1_round
                write(var2_filename,FMT=FMTorder) var2_round
                subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

                allocate(rmsds(population),coords(population,Nvar+6*Natoms))
                call getRMSD(subcell,population,coords_static,rmsds,coords)

                min_position = minloc(rmsds,1)
                min_rmsd = rmsds(min_position)
                closestCoords = coords(min_position,Nvar+1:Nvar+6*Natoms)

                deallocate(rmsds,coords)
        end if
end if

end subroutine getNeighbors



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               FUNCTION GETRMSD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  actually pretty self-explanatory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getRMSD(filename,population,coords_static,rmsds,coords)

use ls_rmsd
use f2_parameters
implicit none
integer,intent(in) :: population
real,intent(out), dimension(population,Nvar+6*Natoms) :: coords
double precision,intent(out), dimension(population) :: rmsds
integer :: i,j
character(50),intent(in) :: filename
double precision,intent(in), dimension(3,Natoms) :: coords_static
double precision, dimension(3,Natoms) :: rmsd_coords2
double precision, allocatable :: U(:,:), g(:,:)
double precision, dimension(3) :: x_center,y_center

!open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
!write(progresschannel,*) "Reading off subcell ", trim(filename), " with ", population, " frames"
!write(progresschannel,*) ""
!close(progresschannel)

!Read off the states in this subcell
open(filechannel1,file=trim(path3)//trim(filename)//".dat")
do i = 1, population
        read(filechannel1,FMT=FMT1,advance="no") (coords(i,j),j=1,Nvar)
        read(filechannel1,FMT=FMT2) (coords(i,j),j=Nvar+1,Nvar+6*Natoms)

        !Need to make the coordinates readable for the rmsd
        rmsd_coords2 = reshape(coords(i,Nvar+1:Nvar+3*Natoms),&
                                (/3, Natoms/))
        call rmsd(Natoms,coords_static,rmsd_coords2,0,U,&
                  x_center,y_center,rmsds(i),.false.,g)
        end do
close(filechannel1)


end subroutine getRMSD




end module checkCells4
