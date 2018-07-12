module addFrametoGrid
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	MODULE ADDFRAMETOGRID
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	FILECHANNELS USED:
!			filechannel1 - file reading and writing
!				(temporary)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      DIVYUP FUNCTION (QuickSort Alogrithm by Tony Hoare)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      INPUT:   character :: subcellParent			"filename of the 'to-be-divided' (parent, order N-1) cell/subcell"
!               character :: FMTorder				"format (decimal places) of the 'to-be-formed (children, order N) subcell"
!	  	real ::  var1_round,var2_round	    		"the number form of the subcell"
!		integer, dim(3) :: SP				"list of parameters for scaling and resolution"
!			SP(1) = scaling of variable1
!			SP(2) = scaling of variable2
!			SP(3) = product of SP(1),SP(2) 
!		real :: multiplier1,multiplier2	           	"the length of children subcells"
!               integer :: indexN               "number of overcrowded (excluding this 'to-be-divided' cell) cells in the parent level"
!                                               " ex. 'indexN * resolution' is the index of the first children subcell in the counter array
!               integer :: lengthN              "dimension of array counterN"
!               integer :: frames               "the number of frames in the parenet cell"
!
!      OUTPUT:  integer :: counterN             "an array of integers, which are the number of frames in the children subcells"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      This function grids an overcrowded cell/subcell into higher order subcells.
!      This function triggers only when the size of a cell/subcell reaches the overcrowd limit for the first time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	FILECHANNELS USED:
!			filechannel1 - file reading and writing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine divyUp(subcellParent,FMTorder,var1_round,var2_round,SP,MP,&
		  indexN,counterN,lengthN,frames,Nfile,path_to_grid)
use PARAMETERS
use FUNCTIONS
implicit none
integer,parameter :: dp = kind(1.0d0)

!SUBCELL FILENAMES
character(*), intent(in) :: subcellParent
character(50) :: subcellChild

!PATH TO GRID HOUSING SUBCELL FILES
character(*), intent(in) :: path_to_grid

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

!CUMULATIVE NUMBER OF FILES MADE
integer, intent(inout) :: Nfile

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

open(filechannel1,file=path_to_grid//trim(subcellParent)//".dat")
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

        open(filechannel1,file=path_to_grid//trim(subcellChild)//".dat",status="new")
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

open(filechannel1,file=path_to_grid//trim(subcellChild)//".dat",status="new")
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
!	      character :: path_to_grid		 "the path to the grid"
!      IN/OUTPUT:  integer :: header1            "the number of orderN-1 cells that are overcrowded"
!                          :: header2		 
!                          :: header3              
!                  integer :: counter0(:)        "the counting/key system for orderN cells"
!                             counter1(:)
!                             counter2(:)
!                             counter3(:)
!		   integer :: Nfile		 "the running total number of files made"
!      OUTPUT:	   logical :: header_max_flag	 "did we run out of room in a counter?"
!		   integer :: order		 "what subcell level did we just add to?"
!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      addState adds a state/frame (its coordinates and gradients stored in coords)
!      to all subcells (stored in a file) corresponding to var1, var2, var3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	FILECHANNELS USED:
!			filechannel1 - file reading and writing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine addState(vals,coords,gradient,&
                header1,header2,header3,&
                counter0,counter1,counter2,counter3,&
                Nfile,header_max_flag,path_to_grid,order)
use PARAMETERS
implicit none
integer,parameter :: dp = kind(1.0d0)

!CUMULATIVE NUMBER OF FILES MADE
integer, intent(inout) :: Nfile

!PATH TO GRID HOUSING SUBCELLS
character(*),intent(in) :: path_to_grid

!CUMULATIVE NUMBER OF ORDER-(N-1) SUBCELLS OVERCROWDED
integer,intent(out) :: header1
integer,intent(out) :: header2
integer,intent(out) :: header3

!INTERNAL COUNTER HOUSING SUBCELL INFORMATION
integer,dimension(counter0_max),intent(out) :: counter0
integer,dimension(counter1_max),intent(out) :: counter1
integer,dimension(counter2_max),intent(out) :: counter2
integer,dimension(counter3_max),intent(out) :: counter3

!VARIABLES AFTER SUBTRACTING OFF ROUNDED VALUES
real :: var1_cell, var2_cell

!VARIABLES ROUNDED DOWN ACCORDING TO SUBCELL ORDER
real :: var1_round0, var1_round1, var1_round2, var1_round3, var1_round
real :: var2_round0, var2_round1, var2_round2, var2_round3, var2_round

!VARIABLES INDEX IN COUNTER
integer :: indexer, key, var1_index,var2_index

!ORDER OF SUBCELL THAT FRAME WAS ADDED TO
integer,intent(out) :: order

!VARIABLES OF THE FRAME
real(dp), dimension(Nvar), intent(in) :: vals

!COORDINATES AND GRADIENT OF THE GRAME
real(dp), dimension(3,Natoms), intent(in) :: coords,gradient

!SUBCELL FILENAMING
character(9) :: descriptor0,var1_filename,var2_filename
character(50) :: subcell

!FLAG INDICATING IF A COUNTER IS FULL
logical,intent(out) :: header_max_flag

!DO LOOP INTEGER INCREMENTALS
integer :: i,j

header_max_flag = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               ORDER 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

order = 0

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
        open(filechannel1,file=path_to_grid//trim(subcell)//".dat",position="append",status=trim(descriptor0))
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
                    header1-1,counter1,counter1_max,overcrowd0-1,Nfile,path_to_grid)
 
	!We still need to add the frame
	!Because we are adding the frame AFTER subdividing, we need to continue on to the
	!next order; that is why there is no return at the end of this conditional
        open(filechannel1,file=path_to_grid//trim(subcell)//".dat",position="append",status="old")
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

!Repeat exactly as before but with different scalings and such

order = 1

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

        open(filechannel1,file=path_to_grid//trim(subcell)//".dat",position="append",status=trim(descriptor0))
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
                    header2-1,counter2,counter2_max,overcrowd1-1,Nfile,path_to_grid)
 
        open(filechannel1,file=path_to_grid//trim(subcell)//".dat",position="append",status="old")
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

order = 2

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

        open(filechannel1,file=path_to_grid//trim(subcell)//".dat",position="append",status=trim(descriptor0))
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
                    header3-1,counter3,counter3_max,overcrowd2-1,Nfile,path_to_grid)
 
!Just for testing purposes
 open(progresschannel,file=trim(path_to_grid)//trim(progressfile),position="append")
 write(progresschannel,*) "divyUp on subcell: ", trim(subcell)
 close(progresschannel)

        open(filechannel1,file=path_to_grid//trim(subcell)//".dat",position="append",status="old")
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

order = 3

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

        open(filechannel1,file=path_to_grid//trim(subcell)//".dat",position="append",status=trim(descriptor0))
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

        open(filechannel1,file=path_to_grid//trim(subcell)//".dat",position="append")
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

        open(filechannel1,file=path_to_grid//trim(subcell)//".dat",position="append")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT3,advance="no") ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
        close(filechannel1)

end if


return


end subroutine addState



end module addFrametoGrid
