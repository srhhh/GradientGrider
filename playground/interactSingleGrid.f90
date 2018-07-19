module interactSingleGrid
use PARAMETERS
use DOUBLE
implicit none


contains


subroutine divyUp(subcellParent,FMTorder,var1_round,var2_round,SP,MP,&
                  indexN,counterN,lengthN,frames)
use PARAMETERS
use FUNCTIONS
implicit none

!SUBCELL FILENAMES
character(*), intent(in) :: subcellParent
character(50) :: subcellChild

!SUBCELL FILENAME FORMAT
character(6), intent(in) :: FMTorder
character(10) :: var1_filename, var2_filename

!VARIABLES OF SUBCELL TO BE DIVIDED
real,intent(in) :: var1_round, var2_round

!SUBCELL SCALING PARAMETERS
integer,dimension(1+Nvar),intent(in) :: SP
integer :: scaling1, scaling2, resolution
real,dimension(2*Nvar),intent(in) :: MP
real :: multiplier1, multiplier2, divisor1, divisor2

!NUMBER OF FRAMES IN SUBCELL
integer, intent(in) :: frames

!COUNTER HOUSING SUBCELLS TO BE MADE
integer, intent(in) :: lengthN
integer,dimension(lengthN),intent(out) :: counterN

!LOCATION OF SUBCELLS TO BE MADE IN COUNTER
integer, intent(in) :: indexN
!INDEX TRACKER FOR SUBCELLS TO BE MADE
integer :: index1, index2, recent_index1, recent_index2
integer :: indexer, recent_vals_index, sorting_index

!ARRAY HOUSING FRAMES
real(dp),dimension(frames,Nvar+6*Natoms) :: coords

!ARRAY HOUSING VARIABLES OF FRAMES
integer,dimension(frames,Nvar) :: vals

!ARRAY HOUSING ORIGINAL POSITIONS OF FRAMES IN SUBCELL
integer,dimension(frames,1) :: originalIndexes

!DO LOOP INCREMENTALS
integer :: i, j, k, l


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          INITIALIZATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

scaling1 = SP(1)
scaling2 = SP(2)
resolution = SP(3)

multiplier1 = MP(1)
multiplier2 = MP(2)
divisor1 = MP(3)
divisor2 = MP(4)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          FRAME RETRIEVAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(filechannel1,file=gridpath2//trim(subcellParent)//".dat")
do i=1, frames
        read(filechannel1,FMT=FMT1,advance="no") (coords(i,l),l=1,Nvar)
        read(filechannel1,FMT=FMT2)(coords(i,Nvar+l),l=1,6*Natoms)
        vals(i,1) = int((coords(i,1) - var1_round) * divisor1)
        if (vals(i,1) < 0) vals(i,1) = 0
        if (vals(i,1) .ge. scaling1_0 ) vals(i,1) = scaling1_0 - 1
        vals(i,2) = int((coords(i,2) - var2_round) * divisor2)
        if (vals(i,2) < 0) vals(i,2) = 0
        if (vals(i,2) .ge. scaling2_0 ) vals(i,2) = scaling2_0 - 1
        originalIndexes(i,1) = i
end do
close(filechannel1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         FRAME SORTING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call qsort2(vals,originalIndexes,frames,Nvar,1,frames,1)

recent_vals_index = 1
sorting_index = 0
do i = 1, frames
        if (vals(i,1) > sorting_index) then
                if (i > recent_vals_index) then
                        call qsort2(vals,originalIndexes,frames,Nvar,recent_vals_index,i-1,2)
                end if
                sorting_index = vals(i,1)
                recent_vals_index = i
        end if
end do
call qsort2(vals,originalIndexes,frames,Nvar,recent_vals_index,frames,2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         FRAME WRITING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sorting_index = 1
recent_vals_index = 1
recent_index1 = vals(1,1)
recent_index2 = vals(1,2)

do sorting_index = 2, frames
index1 = vals(sorting_index,1)
index2 = vals(sorting_index,2)

if ((index1/=recent_index1) .or. (index2/=recent_index2)) then
        write(var1_filename,FMT=FMTorder) var1_round + recent_index1*multiplier1
        write(var2_filename,FMT=FMTorder) var2_round + recent_index2*multiplier2
        subcellChild = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

        open(filechannel1,file=gridpath2//trim(subcellChild)//".dat",status="new")
        do i = recent_vals_index, sorting_index-1
                 write(filechannel1,FMT=FMT1,advance="no") (coords(originalIndexes(i,1),l),l=1,Nvar)
                 write(filechannel1,FMT=FMT2)(coords(originalIndexes(i,1),l),l=1+Nvar,Nvar+6*Natoms)
        end do
        close(filechannel1)

        Nfile = Nfile + 1
        indexer = resolution*indexN + scaling1*recent_index2 + recent_index1 + 1
        counterN(indexer) = sorting_index-recent_vals_index
        recent_vals_index = sorting_index

        recent_index1 = index1
        recent_index2 = index2
end if
end do

write(var1_filename,FMT=FMTorder) var1_round + recent_index1*multiplier1
write(var2_filename,FMT=FMTorder) var2_round + recent_index2*multiplier2
subcellChild = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

open(filechannel1,file=gridpath2//trim(subcellChild)//".dat",status="new")
do i = recent_vals_index, frames
        write(filechannel1,FMT=FMT1,advance="no") (coords(originalIndexes(i,1),l),l=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(originalIndexes(i,1),l),l=1+Nvar,Nvar+6*Natoms)
end do
close(filechannel1)

indexer = resolution*indexN + scaling1*recent_index2 + recent_index1 + 1
counterN(indexer) = frames-recent_vals_index+1
Nfile = Nfile + 1


end subroutine divyUp



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ADDSTATE FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      INPUT:      real :: vals(Nvar)            "variablies used to define cells/subcells"
!                  real :: coords(6*Natoms)      "coordinates and gradients of the system"
!             character :: path_to_grid          "the path to the grid"
!      IN/OUTPUT:  integer :: header1            "the number of orderN-1 cells that are overcrowded"
!                          :: header2            
!                          :: header3              
!                  integer :: counter0(:)        "the counting/key system for orderN cells"
!                             counter1(:)
!                             counter2(:)
!                             counter3(:)
!                  integer :: Nfile              "the running total number of files made"
!      OUTPUT:     logical :: header_max_flag    "did we run out of room in a counter?"
!                  integer :: order              "what subcell level did we just add to?"
!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      addState adds a state/frame (its coordinates and gradients stored in coords)
!      to all subcells (stored in a file) corresponding to var1, var2, var3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       FILECHANNELS USED:
!                       filechannel1 - file reading and writing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine addState(vals,coords,gradient)
use PARAMETERS
implicit none

!VARIABLES AFTER SUBTRACTING OFF ROUNDED VALUES
real :: var1_cell, var2_cell

!VARIABLES ROUNDED DOWN ACCORDING TO SUBCELL ORDER
real :: var1_round0, var1_round1, var1_round2, var1_round3, var1_round
real :: var2_round0, var2_round1, var2_round2, var2_round3, var2_round

!VARIABLES INDEX IN COUNTER
integer :: indexer, key, var1_index,var2_index

!VARIABLES OF THE FRAME
real(dp), dimension(Nvar), intent(in) :: vals

!COORDINATES AND GRADIENT OF THE GRAME
real(dp), dimension(3,Natoms), intent(in) :: coords,gradient

!SUBCELL FILENAMING
character(9) :: descriptor0,var1_filename,var2_filename
character(50) :: subcell

!DO LOOP INTEGER INCREMENTALS
integer :: i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               ORDER 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!var_cell keeps track of the child subcell portion of the decimal
!and throws away the parent subcell portion of the decimal
var1_cell = vals(1)
var2_cell = vals(2)

!var_index keeps track of the index of the child subcell
!inside of its parent subcell
var1_index = int(var1_cell * divisor1_0)
var2_index = int(var2_cell * divisor2_0)

!var_round is the decimal representation of var_index
var1_round0 = multiplier1_0 * var1_index
var2_round0 = multiplier2_0 * var2_index

!indexer accesses the counter (so as to keep track of population) and is uniquely determined by var_index
indexer = bounds1 * var2_index + var1_index + 1

!Increment by one to signify adding a frame;
!key keeps track of population AND the index of the next counter
key = counter0(indexer) + 1

!If its not overcrowded, we need to add frames
if (key < overcrowd0) then

        !keep track of the population
        counter0(indexer) = key

        !The filename is the sum of the decimal portions of each order
        var1_round = var1_round0
        var2_round = var2_round0
        write(var1_filename,FMT=FMTorder0) var1_round
        write(var2_filename,FMT=FMTorder0) var2_round
        subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

        !If this is the first time in the cell, this file has to be made
        if (key == 1) then
                descriptor0 = "new"
                Nfile = Nfile + 1
        else
                descriptor0 = "old"
        end if

        !Add the frame
        open(filechannel1,file=gridpath2//trim(subcell)//".dat",position="append",status=trim(descriptor0))
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT3,advance="no") ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
        close(filechannel1)

        !Return because we have not yet overcrowded the cell
        return

!If the cell is overcrowded (by that exact number) then we need to subdivide it
else if (key == overcrowd0) then

        !Adding by a large number (key_start) assures the portion of the integer
        !that holds the index of the next counter (header1 value)
        !is not incremented by an increment of key
        key = key + key_start*header1
        counter0(indexer) = key

        !The filename is the sum of the decimal portions of each order
        var1_round = var1_round0
        var2_round = var2_round0
        write(var1_filename,FMT=FMTorder0) var1_round
        write(var2_filename,FMT=FMTorder0) var2_round
        subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

        !For more information on how this works, see the subroutine above
        call divyUp(subcell,FMTorder1,var1_round,var2_round,SP0,MP1,&
                    header1-1,counter1,counter1_max,overcrowd0-1)

        !We still need to add the frame
        !Because we are adding the frame AFTER subdividing, we need to continue on to the
        !next order; that is why there is no return at the end of this conditional
        open(filechannel1,file=gridpath2//trim(subcell)//".dat",position="append",status="old")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT3,advance="no") ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
        close(filechannel1)

        !Incrementing header insures that the index this subcell is granted in the next counter is unique
        header1 = header1 + 1

!For testing purposes
if (header1 == header1_max) then
!print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!print *, "      HEADER 1 WILL BE OVERCROWDED SOON"
!print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
header_max_flag = .true.
end if

end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               ORDER 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Norder1 = Norder1 + 1

!Repeat exactly as before but with different scalings and such

var1_cell = var1_cell - var1_round0
var2_cell = var2_cell - var2_round0

var1_index = int(var1_cell * divisor1_1)
if (var1_index < 0) var1_index = 0
if (var1_index .ge. scaling1_0) var1_index = scaling1_0 - 1
var2_index = int(var2_cell * divisor2_1)
if (var2_index < 0) var2_index = 0
if (var2_index .ge. scaling2_0) var2_index = scaling2_0 - 1

var1_round1 = multiplier1_1 * var1_index
var2_round1 = multiplier2_1 * var2_index

! We need int(key/key_start) = the unique header1 value to access counter;
! And by multiplying by resolution (scaling1*scaling2), we assure that each
! subcell of the parent subcell gets its own unique index
indexer = resolution_0*int(key/key_start-1) + scaling1_0*var2_index + var1_index + 1
key = counter1(indexer) + 1

if (key < overcrowd1) then

        counter1(indexer) = key

        var1_round = var1_round0 + var1_round1
        var2_round = var2_round0 + var2_round1

        write(var1_filename,FMT=FMTorder1) var1_round
        write(var2_filename,FMT=FMTorder1) var2_round
        subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

        if (key == 1) then
                descriptor0 = "new"
                Nfile = Nfile + 1
        else
                descriptor0 = "old"
        end if

        open(filechannel1,file=gridpath2//trim(subcell)//".dat",position="append",status=trim(descriptor0))
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT3,advance="no") ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
        close(filechannel1)

        return

else if (key == overcrowd1) then

        key = key + key_start*header2
        counter1(indexer) = key

        var1_round = var1_round0 + var1_round1
        var2_round = var2_round0 + var2_round1

        write(var1_filename,FMT=FMTorder1) var1_round
        write(var2_filename,FMT=FMTorder1) var2_round
        subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

        call divyUp(subcell,FMTorder2,var1_round,var2_round,SP1,MP2,&
                    header2-1,counter2,counter2_max,overcrowd1-1)

        open(filechannel1,file=gridpath2//trim(subcell)//".dat",position="append",status="old")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT3,advance="no") ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
        close(filechannel1)

        header2 = header2 + 1

!For testing purposes
if (header2 == header2_max) then
!print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!print *, "      HEADER 2 WILL BE OVERCROWDED SOON"
!print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
header_max_flag = .true.
end if

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               ORDER 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

var1_cell = var1_cell - var1_round1
var2_cell = var2_cell - var2_round1

var1_index = int(var1_cell * divisor1_2)
var2_index = int(var2_cell * divisor2_2)

var1_round2 = multiplier1_2 * var1_index
var2_round2 = multiplier2_2 * var2_index

indexer = resolution_1*int(key/key_start-1) + scaling1_1*var2_index + var1_index + 1
key = counter2(indexer) + 1

if (key < overcrowd2) then

        counter2(indexer) = key

        var1_round = var1_round0 + var1_round1 + var1_round2
        var2_round = var2_round0 + var2_round1 + var2_round2

        write(var1_filename,FMT=FMTorder2) var1_round
        write(var2_filename,FMT=FMTorder2) var2_round
        subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

        if (key == 1) then
                descriptor0 = "new"
                Nfile = Nfile + 1
        else
                descriptor0 = "old"
        end if

        open(filechannel1,file=gridpath2//trim(subcell)//".dat",position="append",status=trim(descriptor0))
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT3,advance="no") ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
        close(filechannel1)

        return

else if (key == overcrowd2) then

        key = key + key_start*header3
        counter2(indexer) = key

        var1_round = var1_round0 + var1_round1 + var1_round2
        var2_round = var2_round0 + var2_round1 + var2_round2

        write(var1_filename,FMT=FMTorder2) var1_round
        write(var2_filename,FMT=FMTorder2) var2_round
        subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

        call divyUp(subcell,FMTorder3,var1_round,var2_round,SP2,MP3,&
                    header3-1,counter3,counter3_max,overcrowd2-1)

!Just for testing purposes
 open(progresschannel,file=trim(gridpath2)//trim(progressfile),position="append")
 write(progresschannel,*) "divyUp on subcell: ", trim(subcell)
 close(progresschannel)

        open(filechannel1,file=gridpath2//trim(subcell)//".dat",position="append",status="old")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT3,advance="no") ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
        close(filechannel1)

        header3 = header3 + 1

!For testing purposes
if (header3 == header3_max) then
print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
print *, "      HEADER 3 WILL BE OVERCROWDED SOON"
print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
end if
header_max_flag = .true.
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               ORDER 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

var1_cell = var1_cell - var1_round2
var2_cell = var2_cell - var2_round2

var1_index = int(var1_cell * divisor1_3)
var2_index = int(var2_cell * divisor2_3)

indexer = resolution_2*int(key/key_start-1) + scaling1_2*var2_index + var1_index + 1
key = counter3(indexer) + 1

if (key < overcrowd3) then

        counter3(indexer) = key

        var1_round3 = multiplier1_3 * var1_index
        var2_round3 = multiplier2_3 * var2_index

        var1_round = var1_round0 + var1_round1 + var1_round2 + var1_round3
        var2_round = var2_round0 + var2_round1 + var2_round2 + var2_round3

        write(var1_filename,FMT=FMTorder3) var1_round
        write(var2_filename,FMT=FMTorder3) var2_round
        subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

        if (key == 1) then
                descriptor0 = "new"
                Nfile = Nfile + 1
        else
                descriptor0 = "old"
        end if

        open(filechannel1,file=gridpath2//trim(subcell)//".dat",position="append",status=trim(descriptor0))
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT3,advance="no") ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
        close(filechannel1)

else if (key == overcrowd3) then

        counter3(indexer) = key

        var1_round = var1_round0 + var1_round1 + var1_round2 + var1_round3
        var2_round = var2_round0 + var2_round1 + var2_round2 + var2_round3

        write(var1_filename,FMT=FMTorder3) var1_round
        write(var2_filename,FMT=FMTorder3) var2_round
        subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

        open(filechannel1,file=gridpath2//trim(subcell)//".dat",position="append")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT3,advance="no") ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
        close(filechannel1)

!This means the current griding is probably not granular enough
print *, "A third-level subcell got overcrowded!"

else

        var1_round = var1_round0 + var1_round1 + var1_round2 + var1_round3
        var2_round = var2_round0 + var2_round1 + var2_round2 + var2_round3

        write(var1_filename,FMT=FMTorder3) var1_round
        write(var2_filename,FMT=FMTorder3) var2_round
        subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

        open(filechannel1,file=gridpath2//trim(subcell)//".dat",position="append")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT3,advance="no") ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
        close(filechannel1)

end if


return


end subroutine addState



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!n
!              FUNCTION CHECKSTATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!INPUT: real,dim(3*Natoms) coords                       "the to-be-checked frame"
!       integer, dim(...) counter'X'                    "input counters so as to not re-read everytime"
!OUTPUT real, dim(6*Natoms) approx_gradient               "closest frame+gradient"
!       dp min_rmsd                                     "closest frame rmsd"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine checkState(coords,gradient,min_rmsd,&
                      number_of_frames,order,neighbor_check)
use VARIABLES
use PARAMETERS
implicit none
integer :: i,j,k
integer,intent(out),optional :: order,number_of_frames,neighbor_check
integer :: var1_index,var2_index,indexer,index1_1,index1_2,index2_1,index2_2
integer :: key0,key1,key2,key3
logical :: stop_flag,flag1,flag2,flag3,flag4
real(dp) :: var1,var2
real :: var1_cell,var2_cell,var1_round,var2_round
real :: var1_round0,var2_round0,var1_round1,var2_round1,var1_round2,var2_round2,var1_round3,var2_round3
real(dp), dimension(Ncoords), intent(in) :: coords
real(dp), dimension(3,Natoms), intent(out) :: gradient
real(dp), dimension(3,3) :: U
real(dp), intent(inout) :: min_rmsd
real(dp)  :: old_min_rmsd
character(100) :: subcell
character(9) ::  var1_filename, var2_filename

old_min_rmsd = min_rmsd
neighbor_check = 0
number_of_frames = 0
stop_flag = .false.

! Get the variables corresponding to frame
call getVar3(coords,Natoms,var1)
call getVar4(coords,Natoms,var2)

!The coordinates, as they are formatted in getCells and addCells, are the wrong
!shape for ls_rmsd. Thus, we reshape them first

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

!If the key is zero, then that means divyUp was not called on it
!So there are no children subcells, so this subcell must be examined
if (key0 < overcrowd0) then

                !If the population of the parent cell is empty, we should just
                !give up really
                if (key0 == 0)  return

                !Make the name of the subcell
		var1_round = var1_round0
		var2_round = var2_round0
                write(var1_filename,FMT=FMTorder0) var1_round
                write(var2_filename,FMT=FMTorder0) var2_round
                subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

                call getRMSD(subcell,key0,coords,min_rmsd,gradient,U)

		number_of_frames = key0
		order = 0
		neighbor_check = 0
end if

if (min_rmsd < old_min_rmsd) then
	gradient = matmul(U,gradient)
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

if (key1 < overcrowd1) then

        !If there are frames in this cell, then we retrieve those frames
        if (key1 > 0) then
		neighbor_check = 0
		number_of_frames = key1

                !Make the name of the subcell
		var1_round = var1_round0 + var1_round1
		var2_round = var2_round0 + var2_round1
                write(var1_filename,FMT=FMTorder1) var1_round
                write(var2_filename,FMT=FMTorder1) var2_round
                subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

                call getRMSD(subcell,key1,coords,min_rmsd,gradient,U)
	end if

        !If there are no frames in this cell, we can look at one of its
        !neighbors (better than having to look at the parent)
        if ((key1 == 0) .or. (force_Neighbors)) then
		old_min_rmsd = min_rmsd
		neighbor_check = 1

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
                                          coords,min_rmsd,gradient,U,number_of_frames)

                        !Even if all neighbors are empty, it still returns a
                        !minimum rmsd (default is 100.0)
                        !ANY frame is better than none so it stops after it
                        !finds one
                        if (min_rmsd < old_min_rmsd) then
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

	order = 1
end if

if ((min_rmsd < old_min_rmsd).or.(stop_flag)) then
	gradient = matmul(U,gradient)
	return
end if


!Right now we are stopping after we look at order 0 (parent) and order 1 (child) cells
return

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

if (key2 < overcrowd2) then

	if (key2 > 0) then
		neighbor_check = 0
		number_of_frames = key2

                !Make the name of the subcell
		var1_round = var1_round0 + var1_round1 + var1_round2
		var2_round = var2_round0 + var2_round1 + var2_round2
                write(var1_filename,FMT=FMTorder2) var1_round
                write(var2_filename,FMT=FMTorder2) var2_round
                subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

                call getRMSD(subcell,key2,coords,min_rmsd,gradient,U)
	end if

        if ((key2 == 0) .or. (force_Neighbors)) then
		old_min_rmsd = min_rmsd
		neighbor_check = 1

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
                                          coords,min_rmsd,gradient,U,number_of_frames)

                        if (min_rmsd < old_min_rmsd) then
                                stop_flag = .true.
                        end if

                end do
                if (stop_flag) exit
                end do

        end if

	order = 2
end if

if ((min_rmsd < old_min_rmsd).or.(stop_flag)) then
	gradient = matmul(U,gradient)
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

if (key3 < overcrowd3) then

        if (key3 > 0) then
		neighbor_check = 0
		number_of_frames = key3

                !Make the name of the subcell
		var1_round = var1_round0 + var1_round1 + var1_round2 + var1_round3
		var2_round = var2_round0 + var2_round1 + var2_round2 + var2_round3
                write(var1_filename,FMT=FMTorder3) var1_round
                write(var2_filename,FMT=FMTorder3) var2_round
                subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

                call getRMSD(subcell,key3,coords,min_rmsd,gradient,U)
	end if

        if ((key3 == 0) .or. (force_Neighbors)) then
		old_min_rmsd = min_rmsd
		neighbor_check = 1

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
                                          coords,min_rmsd,gradient,U,number_of_frames)

                        if (min_rmsd < old_min_rmsd) then
                                stop_flag = .true.
                        end if
                end do
                if (stop_flag) exit
                end do

        end if

	order = 3
end if

if ((min_rmsd < old_min_rmsd).or.(stop_flag)) then
	gradient = matmul(U,gradient)
	return
end if

print *, "there is a problem"

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
!       real, dim(6*Natoms) approx_gradient               "coords of min rmsd frame"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getNeighbors(scaling1,scaling2,multiplier1,multiplier2,FMTorder,&
                        index_start,index1_1,index1_2,index2_1,index2_2,&
                        var1_round0,var2_round0,counterN,counterN_max,&
                        coords,min_rmsd,gradient,U,number_of_frames)

use PARAMETERS
implicit none
logical :: flag1,flag2,flag3,flag4
integer,intent(in) :: scaling1, scaling2
integer,intent(in) :: index1_1,index1_2,index2_1,index2_2
real,intent(in) :: multiplier1,multiplier2
real,intent(in) :: var1_round0,var2_round0
real :: var1_round,var2_round
integer,intent(in) :: index_start
integer :: indexer,population,key
integer,intent(in) :: counterN_max
integer,intent(in), dimension(counterN_max) :: counterN
integer,intent(inout) :: number_of_frames
character(6),intent(in) :: FMTorder
character(9) :: var1_filename,var2_filename
character(100) :: subcell
real(dp), intent(inout) :: min_rmsd
real(dp),dimension(3,Natoms),intent(inout) :: gradient
real(dp),dimension(3,3),intent(inout) :: U
real(dp),intent(in), dimension(3,Natoms) :: coords

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
		number_of_frames = number_of_frames + population

                !Make the name of the subcell
		var1_round = var1_round0 + index1_1 * multiplier1
		var2_round = var2_round0 + index2_1 * multiplier2
                write(var1_filename,FMT=FMTorder) var1_round
                write(var2_filename,FMT=FMTorder) var2_round
                subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

                call getRMSD(subcell,population,coords,min_rmsd,gradient,U)
        end if
end if

!Rinse and repeat
if ((flag2) .and. (flag3)) then

        indexer = scaling1*(index2_1) + (index1_2) + 1
        key = counterN(index_start+indexer)
	population = modulo(key,key_start)

        if (population > 0) then
		number_of_frames = number_of_frames + population

                !Make the name of the subcell
		var1_round = var1_round0 + index1_2 * multiplier1
		var2_round = var2_round0 + index2_1 * multiplier2
                write(var1_filename,FMT=FMTorder) var1_round
                write(var2_filename,FMT=FMTorder) var2_round
                subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

                call getRMSD(subcell,population,coords,min_rmsd,gradient,U)
        end if
end if

if ((flag1) .and. (flag4)) then

        indexer = scaling1*(index2_2) + (index1_1) + 1
        key = counterN(index_start+indexer)
	population = modulo(key,key_start)

        if (population > 0) then
		number_of_frames = number_of_frames + population

                !Make the name of the subcell
		var1_round = var1_round0 + index1_1 * multiplier1
		var2_round = var2_round0 + index2_2 * multiplier2
                write(var1_filename,FMT=FMTorder) var1_round
                write(var2_filename,FMT=FMTorder) var2_round
                subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

                call getRMSD(subcell,population,coords,min_rmsd,gradient,U)
        end if
end if

if ((flag2) .and. (flag4)) then

        indexer = scaling1*(index2_2) + (index1_2) + 1
        key = counterN(index_start+indexer)
	population = modulo(key,key_start)

        if (population > 0) then
		number_of_frames = number_of_frames + population

                !Make the name of the subcell
		var1_round = var1_round0 + index1_2 * multiplier1
		var2_round = var2_round0 + index2_2 * multiplier2
                write(var1_filename,FMT=FMTorder) var1_round
                write(var2_filename,FMT=FMTorder) var2_round
                subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

                call getRMSD(subcell,population,coords,min_rmsd,gradient,U)
        end if
end if

end subroutine getNeighbors



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               FUNCTION GETRMSD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  actually pretty self-explanatory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getRMSD(filename,population,coords,min_rmsd,gradient,U)

use ls_rmsd_original
use PARAMETERS
implicit none
integer,intent(in) :: population
real(dp),intent(inout), dimension(3,Natoms) :: gradient
real(dp),intent(inout), dimension(3,3) :: U
real(dp),intent(inout) :: min_rmsd
real(dp),dimension(3) :: x_center,y_center
integer :: i,j,k
character(*),intent(in) :: filename
real(dp),dimension(3,Natoms) :: candidate_gradient
real(dp),dimension(3,3) :: candidate_U
real(dp) :: candidate_min_rmsd
real(dp),intent(in), dimension(3,Natoms) :: coords
real(dp), dimension(3,Natoms) :: coords2

!Read off the states in this subcell
open(filechannel1,file=trim(filename)//".dat")
do i = 1, population
	read(filechannel1,FMT=FMT7,advance="no") ((coords2(j,k),j=1,3),k=1,Natoms)

        call rmsd_dp(Natoms,coords2,coords,1,candidate_U(:,:),&
                     x_center,y_center,candidate_min_rmsd)
	
	if (candidate_min_rmsd < min_rmsd) then
		min_rmsd = candidate_min_rmsd
		U = candidate_U
        	read(filechannel1,FMT=FMT3) ((gradient(j,k),j=1,3),k=1,Natoms)
	else
        	read(filechannel1,FMT=FMT3) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
	end if

end do
close(filechannel1)

end subroutine getRMSD




end module interactSingleGrid
