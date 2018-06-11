
module checkCells4
implicit none


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              FUNCTION CHECKSTATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!INPUT: real,dim(3*Natoms) coords                       "the to-be-checked frame"
!       integer, dim(...) counter'X'                    "so as to not re-read everytime"
!OUTPUT real, dim(6*Natoms) closestCoords               "closest frame+gradient"
!       dp min_rmsd                                     "closest frame rmsd"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine checkState(coords,closestCoords,min_rmsd,&
                      counter0,counter1,counter2,counter3)
use f1_variables
use f1_parameters
use mapCellData
implicit none
integer :: i,j,k,key0,key1,key2,key3,var1_new,var2_new
integer :: index_start,indexer,index1_1,index1_2,index2_1,index2_2
integer :: order,population,decimals1,decimals2
integer :: min_position
integer, dimension(nint(max_var1*max_var2)) :: counter0
integer, dimension(counter1_max) :: counter1
integer, dimension(counter2_max) :: counter2
integer, dimension(counter3_max) :: counter3
logical :: stop_flag,flag1,flag2,flag3,flag4
real :: var1,var2,var3, var1_old, var2_old
real, dimension(3*Natoms), intent(in) :: coords
real, dimension(6*Natoms), intent(out) :: closestCoords
real, dimension(6*Natoms) :: candidateCoords
double precision, dimension(3,Natoms) :: rmsd_coords1,rmsd_coords2
double precision, dimension(3) :: x_center,y_center
double precision, intent(out) :: min_rmsd
double precision  :: candidate_rmsd
double precision, allocatable :: neighbor_rmsds(:)
real, allocatable :: neighbor_coords(:,:)
double precision, allocatable :: U(:,:), g(:,:)
character(50) :: descriptor1, descriptor2, subcell
character(9) ::  descriptor3, descriptor4

! Need a ridiculously large number
min_rmsd = 100.0

! Get the variables corresponding to frame
! RS: What do you mean by 'each'? Isn't the input only one frame?
! RS: calculating var again so you don't have to pass them?
!                       KF: the word 'each' is misleading, my bad
call getVar1(coords,Natoms,var1)
call getVar2(coords,Natoms,var2)
call getVar3(coords,Natoms,var3)

!The coordinates, as they are formatted in getCells and addCells, are the wrong
!shape for ls_rmsd. Thus, we reshape them first

! RS: ha! Should have read them in this format at the first place
!                       KF: ohhh maaaannnnnnn
rmsd_coords1 = reshape(coords,(/3, Natoms/))

! RS: I have some thoughts on the following -- Let's talk tomorrow


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 ORDER 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

order = 0

!The first subcell will be the digits to the left of the decimal
!The 'integer part'
var1_old = anint(var1 - 0.5)
var2_old = anint(var2 - 0.5)
write(descriptor3,FMT="(F9.0)") var1_old
write(descriptor4,FMT="(F9.0)") var2_old
subcell = trim(adjustl(descriptor3))//"_"//trim(adjustl(descriptor4))

!To retrieve the index of this cell, we follow the same
!indexing scheme as in addCells.f90
var1_new = nint(var1_old)
var2_new = nint(var2_old)
indexer = nint(max_var1)*(var2_new-1) + var1_new

!The population will be the digits on the right of the key
!The 'key' or index to the next counter are the remaining digits
population = modulo(counter0(indexer),key_start)
key0 = counter0(indexer)/key_start

!Write to the progress file for bug-testing
open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
write(progresschannel,*) "Investigating key ", counter0(indexer), " with index ", indexer, " and order ", order
write(progresschannel,*) ""
close(progresschannel)

!If the key is zero, then that means divyUp was not called on it
!So there are no children subcells, so this subcell must be examined
if (key0 == 0) then

                !If the population of the parent cell is empty, we should just
                !give up really
                if (population == 0) return

                !Call mapCell to see the heat map of all parent cells
                call mapCell(0,counter0,counter0_max,nint(max_var1),nint(max_var2))

                ! getRMSD reads off the coordinates and calculated the RMSD with
                !ls_rmsd module
                allocate(neighbor_rmsds(population),&
                        neighbor_coords(population,Nvar+6*Natoms))
                call getRMSD(subcell,population,rmsd_coords1,&
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

order = order + 1

!The name of the subcell is a multiple of .25*.1**order
var1_old = anint(var1*scaling1_0-0.5)/scaling1_0
var2_old = anint(var2*scaling2_0-0.5)/scaling2_0

!However many multiples of .25 there are is the index of the variable
var1_new = modulo(nint(var1_old*scaling1_0),scaling1_0)
var2_new = modulo(nint(var2_old*scaling2_0),scaling2_0)

!The index of this subcell in counter1 comes from the key calculated earlier
index_start = resolution_0*(key0)
indexer = index_start + scaling1_0*var2_new + var1_new

!Population and key calculation as before
population = modulo(counter1(indexer),key_start)
key1 = counter1(indexer)/key_start

!Write to the progress file
open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
write(progresschannel,*) "Investigating key ", counter1(indexer), " with index ", indexer, " and order ", order
write(progresschannel,*) ""
close(progresschannel)

!If the key is zero, then that means divyUp was not called on it
!So there are no children subcells, so this subcell must be examined
if (key1 == 0) then

        !Call mapCell just to see the heatmap of this subcell
        call mapCell(index_start,counter1,counter1_max,scaling1_0,scaling2_0)

        !If there are frames in this cell, then we retrieve those frames
        if (population > 0) then

                !Make the name of the subcell
                write(descriptor3,FMT="(F9.2)") var1_old
                write(descriptor4,FMT="(F9.2)") var2_old
                subcell = trim(adjustl(descriptor3))//"_"//trim(adjustl(descriptor4))

                !Get the RMSDs of the frames inside of it
                allocate(neighbor_rmsds(population),&
                        neighbor_coords(population,Nvar+6*Natoms))
                call getRMSD(subcell,population,rmsd_coords1,&
                                neighbor_rmsds,neighbor_coords)

                !Obtain the rmsd and coordinates of the closest frame
                min_position = minloc(neighbor_rmsds,1)
                min_rmsd = neighbor_rmsds(min_position)
                closestCoords = neighbor_coords(min_position,Nvar+1:Nvar+6*Natoms)

                deallocate(neighbor_rmsds,neighbor_coords)

        !If there are no frames in this cell, we can look at one of its
        !neighbors (better than having to look at the parent)
        else
                
                !In this case, we need to make the names of the subcells
                !Start by retrieving the integer part
                write(descriptor3,FMT="(F9.0)") var1_old-0.5
                write(descriptor4,FMT="(F9.0)") var2_old-0.5
                descriptor3 = adjustl(descriptor3)
                descriptor4 = adjustl(descriptor4)

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
                        index1_1 = var1_new + j
                        index1_2 = var1_new - j

                        index2_1 = var2_new + i - j
                        index2_2 = var2_new - i + j

                        !This does the nitty-gritty of checking to see
                        !If any of the indexes are out of bounds and 
                        !obtaining their rmsds and coordinates
                        call getNeighbors(scaling1_0,scaling2_0,&
                                          index_start, index1_1,index1_2,index2_1,index2_2,&
                                          descriptor3,descriptor4,&
                                          0,0,order,counter1,counter1_max,&
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

!Code is basically the same but changing variables slightly to accomodate
!for the different scalings, resolutions, and counters

order = order + 1

var1_new = modulo(nint(var1*scaling1_0*scaling1_1-0.5),scaling1_1)
var2_new = modulo(nint(var2*scaling2_0*scaling2_1-0.5),scaling2_1)

index_start = resolution_1*(key1)
indexer = index_start + scaling1_1*var2_new + var1_new

population = modulo(counter2(indexer),key_start)
key2 = counter2(indexer)/key_start

open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
write(progresschannel,*) "Investigating key ", counter2(indexer), " with index ", indexer, " and order ", order
write(progresschannel,*) ""
close(progresschannel)

if (key2 == 0) then

        call mapCell(index_start,counter2,counter2_max,scaling1_1,scaling2_1)

!Force this to look at its neighbors (for testing)
!Remark: with normals parameters, the deepest subcell for a particular frame is
!usually of order 2 so this 
population = 0 

        if (population > 0) then

                var1_old = anint(var1*scaling1_0*scaling1_1-0.5)/(scaling1_0*scaling1_1)
                var2_old = anint(var2*scaling2_0*scaling2_1-0.5)/(scaling2_0*scaling2_1)

                write(descriptor3,FMT="(F9.3)") var1_old
                write(descriptor4,FMT="(F9.3)") var2_old
                subcell = trim(adjustl(descriptor3))//"_"//trim(adjustl(descriptor4))


                allocate(neighbor_rmsds(population),&
                        neighbor_coords(population,Nvar+6*Natoms))
                call getRMSD(subcell,population,rmsd_coords1,&
                                neighbor_rmsds,neighbor_coords)

                min_position = minloc(neighbor_rmsds,1)
                min_rmsd = neighbor_rmsds(min_position)
                closestCoords = neighbor_coords(min_position,Nvar+1:Nvar+6*Natoms)

                deallocate(neighbor_rmsds,neighbor_coords)
        else

!For bug-testing
open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
write(progresschannel,*) "Fetching neighbors..."
close(progresschannel)

                var1_old = anint(var1*scaling1_0-0.5)/(scaling1_0)
                var2_old = anint(var2*scaling2_0-0.5)/(scaling2_0)

                write(descriptor3,FMT="(F9.0)") var1_old-0.5
                write(descriptor4,FMT="(F9.0)") var2_old-0.5
                descriptor3 = adjustl(descriptor3)
                descriptor4 = adjustl(descriptor4)

                decimals1 = modulo(nint(var1_old*1000),1000)
                decimals2 = modulo(nint(var2_old*1000),1000)

!For bug-testing
open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
write(progresschannel,*) "   inside of cell ", trim(adjustl(descriptor3)), decimals1/10,&
                                  trim(adjustl(descriptor4)), decimals2/10
write(progresschannel,*) ""
close(progresschannel)

                stop_flag = .false.
                do i = 1, scaling1_1
                do j = 0, i
                
                        index1_1 = var1_new + j
                        index1_2 = var1_new - j

                        index2_1 = var2_new + i - j
                        index2_2 = var2_new - i + j

                        call getNeighbors(scaling1_1,scaling2_1,&
                                          index_start,index1_1,index1_2,index2_1,index2_2,&
                                          descriptor3,descriptor4,&
                                          decimals1,decimals2,order,counter2,counter2_max,&
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

!Code is basically the same but changing variables slightly to accomodate
!for the different scalings, resolutions, and counters

order = order + 1

var1_new = modulo(nint(var1*scaling1_0*scaling1_1**2-0.5),scaling1_1)
var2_new = modulo(nint(var2*scaling2_0*scaling2_1**2-0.5),scaling2_1)

index_start = resolution_1*(key2)
indexer = index_start + scaling1_1*var2_new + var1_new

population = modulo(counter3(indexer),key_start)
key3 = counter3(indexer)/key_start

open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
write(progresschannel,*) "Investigating key ", counter3(indexer), " with index ", indexer, " and order ", order
write(progresschannel,*) ""
close(progresschannel)

if (key3 == 0) then

        if (population > 0) then

                var1_old = anint(var1*scaling1_0*scaling1_1**2-0.5)/(scaling1_0*scaling1_1**2)
                var2_old = anint(var2*scaling2_0*scaling2_1**2-0.5)/(scaling2_0*scaling2_1**2)

                write(descriptor3,FMT="(F9.3)") var1_old
                write(descriptor4,FMT="(F9.3)") var2_old
                subcell = trim(adjustl(descriptor3))//"_"//trim(adjustl(descriptor4))


                allocate(neighbor_rmsds(population),&
                        neighbor_coords(population,Nvar+6*Natoms))
                call getRMSD(subcell,population,rmsd_coords1,&
                                neighbor_rmsds,neighbor_coords)

                min_position = minloc(neighbor_rmsds,1)
                min_rmsd = neighbor_rmsds(min_position)
                closestCoords = neighbor_coords(min_position,Nvar+1:Nvar+6*Natoms)

                deallocate(neighbor_rmsds,neighbor_coords)
        else
 
                var1_old = anint(var1*scaling1_0*scaling1_1-0.5)/(scaling1_0*scaling1_1)
                var2_old = anint(var2*scaling2_0*scaling2_1-0.5)/(scaling2_0*scaling2_1)
               
                write(descriptor3,FMT="(F9.0)") var1_old-0.5
                write(descriptor4,FMT="(F9.0)") var2_old-0.5
                descriptor3 = adjustl(descriptor3)
                descriptor4 = adjustl(descriptor4)

                decimals1 = modulo(nint(var1_old*10000),10000)
                decimals2 = modulo(nint(var2_old*10000),10000)

                stop_flag = .false.
                do i = 1, scaling1_1
                do j = 0, i
                
                        index1_1 = var1_new + j
                        index1_2 = var1_new - j

                        index2_1 = var2_new + i - j
                        index2_2 = var2_new - i + j

                        call getNeighbors(scaling1_1,scaling2_1,&
                                          index_start,index1_1,index1_2,index2_1,index2_2,&
                                          descriptor3,descriptor4,&
                                          decimals1,decimals2,order,counter2,counter2_max,&
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

print *, "this is a massive error"

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
subroutine getNeighbors(scaling1,scaling2,&
                        index_start,index1_1,index1_2,index2_1,index2_2,&
                        integer1,integer2,decimals1,decimals2,&
                        order,counterN,counterN_max,&
                        coords_static,min_rmsd,closestCoords)

use f1_parameters
implicit none
logical :: flag1,flag2,flag3,flag4
integer,intent(in) :: order,scaling1, scaling2
integer,intent(in) :: index1_1,index1_2,index2_1,index2_2
integer,intent(in) :: decimals1, decimals2
integer,intent(in) :: index_start
integer :: indexer,population
integer,intent(in) :: counterN_max
integer,intent(in), dimension(counterN_max) :: counterN
character(9),intent(in) :: integer1,integer2
character(9) :: descriptor1,descriptor2
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
        indexer = scaling1*(index2_1)+(index1_1)
        population = modulo(counterN(index_start+indexer),key_start)

        !Only read of coordinates if there are any
        if (population > 0) then

                !The integer portion of the filename is the same
                !but the decimal part is incremented
                !Because the decimal portion is from the PARENT
                !all we need to do is add by .25*.1**order
                write(descriptor1,FMT="(I9)") decimals1 + 25*index1_1
                write(descriptor2,FMT="(I9)") decimals2 + 25*index2_1
                subcell = trim(integer1)//trim(adjustl(descriptor1))//"_"//&
                              trim(integer2)//trim(adjustl(descriptor2))

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

        indexer = scaling1*(index2_1)+(index1_2)
        population = modulo(counterN(index_start+indexer),key_start)

        if (population > 0) then
                write(descriptor1,FMT="(I9)") decimals1 + 25*index1_2
                write(descriptor2,FMT="(I9)") decimals2 + 25*index2_1
                subcell = trim(integer1)//trim(adjustl(descriptor1))//"_"//&
                              trim(integer2)//trim(adjustl(descriptor2))

                allocate(rmsds(population),coords(population,Nvar+6*Natoms))
                call getRMSD(subcell,population,coords_static,rmsds,coords)

                min_position = minloc(rmsds,1)
                min_rmsd = rmsds(min_position)
                closestCoords = coords(min_position,Nvar+1:Nvar+6*Natoms)

                deallocate(rmsds,coords)
        end if
end if

if ((flag1) .and. (flag4)) then

        indexer = scaling1*(index2_2)+(index1_1)
        population = modulo(counterN(index_start+indexer),key_start)

        if (population > 0) then
                write(descriptor1,FMT="(I9)") decimals1 + 25*index1_1
                write(descriptor2,FMT="(I9)") decimals2 + 25*index2_2
                subcell = trim(integer1)//trim(adjustl(descriptor1))//"_"//&
                              trim(integer2)//trim(adjustl(descriptor2))

                allocate(rmsds(population),coords(population,Nvar+6*Natoms))
                call getRMSD(subcell,population,coords_static,rmsds,coords)

                min_position = minloc(rmsds,1)
                min_rmsd = rmsds(min_position)
                closestCoords = coords(min_position,Nvar+1:Nvar+6*Natoms)

                deallocate(rmsds,coords)
        end if
end if

if ((flag2) .and. (flag4)) then

        indexer = scaling1*(index2_2)+(index1_2)
        population = modulo(counterN(index_start+indexer),key_start)

        if (population > 0) then
                write(descriptor1,FMT="(I9)") decimals1 + 25*index1_2
                write(descriptor2,FMT="(I9)") decimals2 + 25*index2_2
                subcell = trim(integer1)//trim(adjustl(descriptor1))//"_"//&
                              trim(integer2)//trim(adjustl(descriptor2))

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
use f1_parameters
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

open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
write(progresschannel,*) "Reading off subcell ", trim(filename)
write(progresschannel,*) ""
close(progresschannel)

end subroutine getRMSD




end module checkCells4
