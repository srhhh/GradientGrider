module addCells5
implicit none


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      DIVYUP FUNCTION (QuickSort Alogrithm by Tony Hoare)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      INPUT:   character :: subcellParent			"filename of the 'to-be-divided' (parent, order N-1) cell/subcell"
!               character :: FMTorder				"format (decimal places) of the 'to-be-formed (children, order N) subcell"
!	  	real ::  var1_round,var2_round	    		"the number form of the subcell"
!		integer, dim(3) :: SP				"list of parameters for scaling and resolution"
!			SP(1) = scaling of variable1
!			SP(2) = scaling of variable2
!			SP(3) = product of SP(1),SP(2) 
! RS: we should make a list explaining what those scaling parameters are
!			KF: alright
!		real :: multiplier1,multiplier2	           	"the length of children subcells"
! RS: can the above be part of the SP?
!			KF: I don't know; its a separate type so I don't think so
!               integer :: indexN               "number of overcrowded (excluding this 'to-be-divided' cell) cells in the parent level"
!                                               " ex. 'indexN * resolution' is the index of the first children subcell in the counter array
!               integer :: lengthN              "dimension of array counterN"
!               integer :: frames               "the number of frames in the parenet cell"
!
!      OUTPUT:  integer :: counterN             "an array of integers, which are the number of frames in the children subcells"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      This function grids an overcrowded cell/subcell into higher order subcells.
!      This function triggers only when the size of a cell/subcell reaches the overcrowd limit for the first time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RS: I like the new nomenclature! :D
!		KF: huzzahhhh

subroutine divyUp(subcellParent,FMTorder,var1_round,var2_round,SP,&
		  multiplier1,multiplier2,indexN,counterN,lengthN,frames)
use f2_parameters
use f1_functions
implicit none
integer, intent(in) :: indexN,lengthN,frames
integer, dimension(3),intent(in) :: SP
real,intent(in) :: multiplier1,multiplier2
integer, dimension(lengthN), intent(out) :: counterN
integer :: i,j,k,l,population,last_marker,Nmarkers,Nsortings,frames_plus_markers
integer :: indexer,counter_index
integer :: scaling1,scaling2,resolution
real,intent(in) :: var1_round,var2_round
real,allocatable :: coords(:,:),vals(:,:)
integer, allocatable :: originalIndexes(:,:)
character(50) :: subcellChild
character(50),intent(in) :: subcellParent
character(6) :: FMTorder
character(10) :: var1_filename,var2_filename

scaling1 = SP(1)
scaling2 = SP(2)
resolution = SP(3)
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          FRAME RETRIEVAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Open up the file, read the variables, coordinates, and gradients
!vals    ---  holds the criteria of sorting (three for now)
!coords  ---  holds the xyz coordinates and gradients (6*Natoms)
!        ---  coords will be sorted according to vals
!indexer ---  holds the indexes of coords
!        ---  indexer will be sorted according to vals
frames_plus_markers = frames+resolution-1
open(filechannel1,file=path3//trim(subcellParent)//".dat")
allocate(coords(frames,6*Natoms))
allocate(vals(frames_plus_markers,Nvar))
allocate(originalIndexes(frames_plus_markers,1))
do j=1, frames
        read(filechannel1,FMT=FMT1,advance="no",iostat=k) (vals(j,i),i=1,Nvar)
        read(filechannel1,FMT=FMT2) (coords(j,i),i=1,6*Natoms)
        originalIndexes(j,1) = j
end do
close(filechannel1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         FRAME SORTING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Add "markers" to know when a subcell begins/ends
!ex:    Let     scaling1 = scaling2 = 3
!               resolution = 9
!               var1 = int(N/3)
!               var2 = N%3
!		var1_round = 0.0
!		var2_round = 0.0
!		gap1 = 1.0
!		gap2 = 1.0
!
!               val1 = 0.0 + 1.0*var1
!               val2 = 0.0 + 1.0*var2
!
!       index N     (var1,var2)       (val1,val2)
!          1            0,1             0.0,1.0
!          2            0,2             0.0,2.0
!          3            1,0             1.0,0.0
!         ...           ...               ...
!          7            2,1             2.0,1.0
!          8            2,2             2.0,2.0

do i = 1, resolution-1
        counter_index = frames + i
        originalIndexes(counter_index,1) = counter_index
! 	If you want an integer, please make sure to add INT or ANIT, it would not increase the cost and it makes it safe
!			KF: resolved	
        vals(counter_index,:) = (/ var1_round+int(i/scaling2)*multiplier1,&
                                   var2_round+modulo(i,scaling2)*multiplier2 /)
end do

!Say the array vals looks something like this
!
!               ORIGINAL (16)
!                       . . ..   ..   ..   .  .   ..  ..  ..
!              -------@-------------------------------------@------> (var1)
!                    0.0                                   3.0
!              
!Then after the above we have something like this:
!
!               ORIGINAL + MARKERS (24)
!                       . . ..   ..   ..   .  .   ..  ..  ..xxXxxXxx
!              -------@-------------------------------------@------>  (var1)
!                    0.0                                   3.0
!
!Where x indicates a pair where var2 /= 0
!  and X indicates a pair where var2 == 0
!  and . indicates a frame
!

call qsort2(vals,originalIndexes,frames_plus_markers,Nvar,1,frames_plus_markers,1)

!Before quicksorting, we had something like this:
!
!               ORIGINAL + MARKERS
!                       . . ..   ..   ..   .  .   ..  ..  ..xxXxxXxx
!              -------@-------------------------------------@------> (var1)
!                    0.0                                   3.0
!
!Where x indicates a pair where var2 /= 0
!  and X indicates a pair where var2 == 0
!  and . indicates a frame
!
!After quicksorting, we have something like this:
!
!               VAR1 QSORT [1,24]
!                     xx. . ..Xxx..   ..   .  .Xxx..  ..  ..
!              -------@-------------------------------------@------>
!                    0.0                                   3.0
!

last_marker = 1
Nmarkers = 0
Nsortings = 1
do counter_index = 1, frames_plus_markers
        if (originalIndexes(counter_index,1) > frames) then
        Nmarkers = Nmarkers + 1
        if (Nmarkers == scaling2) then
                call qsort2(vals,originalIndexes,frames_plus_markers,Nvar,&
                           last_marker,counter_index-1,2)
                last_marker = counter_index
                Nmarkers = 0
                Nsortings = Nsortings+1
                if (Nsortings == scaling1) exit
        end if
        end if

end do
call qsort2(vals,originalIndexes,frames_plus_markers,Nvar,last_marker,frames_plus_markers,2)

!Before this second-variable quicksorting, we had something like this:
!
!               VAR1 QSORT [1,24]
!                     xx. . ..Xxx..   ..   .  .Xxx..  ..  ..
!              -------@-------------------------------------@------>
!                    0.0                                   3.0
!
!Where x indicates a pair where var2 /= 0
!  and X indicates a pair where var2 == 0
!  and . indicates a frame
!
!In the second-variable quicksort, we check whether a position corresponds to a X
!Every third x is an X, which indicates we have reached the end of a column.
!Then we need to sort all the . and x of that column
!In this case, the X are at 7 and 16, so we need to sort:
!       [1,7), [7,16), and [16,24]
!Columns [1,7) and [7,16) are done in the do-loop whereas the last one is
!done outside of the do-loop.
!
!After this second-variable quicksorting, we have something like this:
!
!               VAR2 QSORT [1,7), [7,16), [16,24]
!                       .x.x..X  .. x ..   .x .X  ..x ..x ..
!              -------@-------------------------------------@------>
!                    0.0                                   3.0
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          FRAME GRIDING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Now that vals and indexer are fully sorted, we need to identify subcells
! RS: What are p1 and p2, indexes of two adjacent markers?

!Any pair of markers (p1,p2) such that (p2-p1 > 1) is a subcell with a frame in
!it; these pairs are ordered in such a way so that the subcell (i,j) that
!corresponds with pair n is:
!       i = n / scaling2
!       j = n % scaling2
!And, of course, i and j start at zero.

last_marker = 1
Nmarkers = 0
do counter_index = 1, frames_plus_markers
        if (originalIndexes(counter_index,1) > frames) then
        population = counter_index - last_marker
        if (population > 0) then
                i = int(Nmarkers/scaling2)
                j = modulo(Nmarkers,scaling2)

                !var_round are only the decimal part of the PARENT cell
                !thus, we have to add on multiplier*var_index to get the decimal
                !portion of the CHILD subcell
                write(var1_filename,FMT=FMTorder) var1_round + i*multiplier1
                write(var2_filename,FMT=FMTorder) var2_round + j*multiplier2
                subcellChild = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

                !We write all of the frames onto the higher order subcell;
                !last_marker represents the first frame of the subcell
                !whereas counter_index is the position of the marker that ends
                !the subcell; thus, one is subctracted from it
                open(filechannel1,file=path3//trim(subcellChild)//".dat",status="new")
                do k = last_marker, counter_index-1
                        write(filechannel1,FMT=FMT1,advance="no") (vals(k,l),l=1,Nvar)
                        write(filechannel1,FMT=FMT2)(coords(originalIndexes(k,1),l),l=1,6*Natoms)                       
                end do 
                close(filechannel1)

                !And of course, we need to add this onto counterN
                !The indexing is exactly as in addState
                indexer = resolution*indexN + scaling1*j + i + 1
                counterN(indexer) = population

        end if

        !After we finish a pair of markers (p1,p2) we save p2 so as to make
        !the next pair (p2,p3); we add one to represent the fact that the first
        !frame of the subcell would be the position AFTER the marker
        last_marker = counter_index + 1
        Nmarkers = Nmarkers + 1

        end if

end do

!We do a separate write statement for the last subcell because there is no
!marker to represent the end; this will automatically have all of the frames
!that have not yet been written by the last_marker
!
!Exactly the same format as before; i,j are both at their last values
population = frames_plus_markers - last_marker
if (population > 0) then
        i = scaling1 - 1
        j = scaling2 - 1

        write(var1_filename,FMT=FMTorder) var1_round + i*multiplier1
        write(var2_filename,FMT=FMTorder) var2_round + j*multiplier2
        subcellChild = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))


        open(filechannel1,file=path3//trim(subcellChild)//".dat",status="new")
        do k = last_marker, frames+resolution-1
                write(filechannel1,FMT=FMT1,advance="no") (vals(k,l),l=1,Nvar)
                write(filechannel1,FMT=FMT2)(coords(originalIndexes(k,1),l),l=1,6*Natoms)                       
        end do 
        close(filechannel1)

        indexer = resolution*indexN + scaling1*j + i + 1
        counterN(indexer) = population    
end if

end subroutine divyUp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ADDSTATE FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      INPUT:      real :: vals(:)               "variablies used to define cells/subcells"
!                  real :: coords(:)             "coordinates and gradients of the system"
!      OUTPUT:     integer :: header1            "the number of orderN-1 cells that are overcrowded"
!                          :: header2		 
!                          :: header3              
!      IN/OUTPUT   integer :: counter0(:)        "the counting/key system for orderN cells"
!                             counter1(:)
!                             counter2(:)
!                             counter3(:)
!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      addState adds a state/frame (its coordinates and gradients stored in coords)
!      to all subcells (stored in a file) corresponding to var1, var2, var3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine addState(vals,coords,&
                header1,header2,header3,&
                counter0,counter1,counter2,counter3)
use f2_parameters
implicit none
integer :: indexer,i,j,key
integer,intent(out) :: header1
integer,intent(out) :: header2
integer,intent(out) :: header3
integer,dimension(counter0_max),intent(out) :: counter0
integer,dimension(counter1_max),intent(out) :: counter1
integer,dimension(counter2_max),intent(out) :: counter2
integer,dimension(counter3_max),intent(out) :: counter3
real :: var1_cell, var2_cell
real :: var1_round0, var1_round1, var1_round2, var1_round3, var1_round
real :: var2_round0, var2_round1, var2_round2, var2_round3, var2_round
integer :: var1_index,var2_index,order
real, dimension(Nvar), intent(in) :: vals
real, dimension(6*Natoms), intent(in) :: coords
character(9) :: descriptor0,var1_filename,var2_filename
character(50) :: subcell
integer :: c1,c2,cr,cm
real :: system_clock_rate

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
! RS: Again, this is only true if the above statement is true
var1_round0 = multiplier1_0 * var1_index
var2_round0 = multiplier2_0 * var2_index

!indexer accesses the counter (so as to keep track of population) and is uniquely determined by var_index
! RS: anint(max_var1/spacing1) as a constant stored in parameter
! RS: e.g. max_var=12 and spacing = 3.6, instead of four bins, you would like to have 3 bins?
! RS: the minimum value of indexer is zero -- i would suggest the index of an array starts from 1
!		KF: addressed with the new parameter 'bounds' which uses ceiling instead
!		KF: also, good point on indexer starting at zero/one
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
        else
                descriptor0 = "old"
        end if

	!Add the frame
        open(filechannel1,file=path3//trim(subcell)//".dat",position="append",status=trim(descriptor0))
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)

	!Return because we have not yet overcrowded the cell
	return

!If the cell is overcrowded (by that exact number) then we need to subdivide it
else if (key == overcrowd0) then

!For testing purposes
if (header1 == 900) then
print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
print *, "      HEADER 1 WILL BE OVERCROWDED SOON"
print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
end if

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
        call divyUp(subcell,FMTorder1,var1_round,var2_round,SP0,multiplier1_1,multiplier2_1,&
                    header1-1,counter1,counter1_max,overcrowd0-1)
 
	!We still need to add the frame
	!Because we are adding the frame AFTER subdividing, we need to continue on to the
	!next order; that is why there is no return at the end of this conditional
        open(filechannel1,file=path3//trim(subcell)//".dat",position="append",status="old")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)
      
        !Incrementing header insures that the index this subcell is granted in the next counter is unique
        header1 = header1 + 1

end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               ORDER 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Repeat exactly as before but with different scalings and such

order = 1

var1_cell = var1_cell - var1_round0
var2_cell = var2_cell - var2_round0

var1_index = int(var1_cell * divisor1_1)
var2_index = int(var2_cell * divisor2_1)

var1_round1 = multiplier1_1 * var1_index
var2_round1 = multiplier2_1 * var2_index

! We need int(key/key_start) = the unique header1 value to access counter;
! And by multiplying by resolution (scaling1*scaling2), we assure that each
! subcell of the parent subcell gets its own unique index
! RS: I would add int to make sure it comes out an integer
! RS: different compilers can mess things up
! RS: Again, will this not make the beginning of the array blank?
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
        else
                descriptor0 = "old"
        end if

        open(filechannel1,file=path3//trim(subcell)//".dat",position="append",status=trim(descriptor0))
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)

	return

else if (key == overcrowd1) then

!For testing purposes
if (header2 == 4900) then
print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
print *, "      HEADER 2 WILL BE OVERCROWDED SOON"
print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
end if


	key = key + key_start*header2
	counter1(indexer) = key

	var1_round = var1_round0 + var1_round1
	var2_round = var2_round0 + var2_round1

	write(var1_filename,FMT=FMTorder1) var1_round
	write(var2_filename,FMT=FMTorder1) var2_round
	subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

        call divyUp(subcell,FMTorder2,var1_round,var2_round,SP1,multiplier1_2,multiplier2_2,&
                    header2-1,counter2,counter2_max,overcrowd1-1)
 
        open(filechannel1,file=path3//trim(subcell)//".dat",position="append",status="old")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)
      
        header2 = header2 + 1

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
        else
                descriptor0 = "old"
        end if

        open(filechannel1,file=path3//trim(subcell)//".dat",position="append",status=trim(descriptor0))
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)

	return

else if (key == overcrowd2) then

!For testing purposes
if (header3 == 4900) then
print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
print *, "      HEADER 3 WILL BE OVERCROWDED SOON"
print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
end if


	key = key + key_start*header3
	counter2(indexer) = key

	var1_round = var1_round0 + var1_round1 + var1_round2
	var2_round = var2_round0 + var2_round1 + var2_round2

	write(var1_filename,FMT=FMTorder2) var1_round
	write(var2_filename,FMT=FMTorder2) var2_round
	subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

        call divyUp(subcell,FMTorder3,var1_round,var2_round,SP2,multiplier1_3,multiplier2_3,&
                    header3-1,counter3,counter3_max,overcrowd2-1)
 
!Just for testing purposes
 open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
 write(progresschannel,*) "divyUp on subcell: ", trim(subcell)
 close(progresschannel)

        open(filechannel1,file=path3//trim(subcell)//".dat",position="append",status="old")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)
      
        header3 = header3 + 1

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
        else
                descriptor0 = "old"
        end if

        open(filechannel1,file=path3//trim(subcell)//".dat",position="append",status=trim(descriptor0))
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)

else if (key == overcrowd3) then

	counter3(indexer) = key

	var1_round = var1_round0 + var1_round1 + var1_round2 + var1_round3
	var2_round = var2_round0 + var2_round1 + var2_round2 + var2_round3

	write(var1_filename,FMT=FMTorder3) var1_round
	write(var2_filename,FMT=FMTorder3) var2_round
	subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

        open(filechannel1,file=path3//trim(subcell)//".dat",position="append")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)

!This means the current griding is probably not granular enough
print *, "A third-level subcell got overcrowded!"

else

	var1_round = var1_round0 + var1_round1 + var1_round2 + var1_round3
	var2_round = var2_round0 + var2_round1 + var2_round2 + var2_round3

	write(var1_filename,FMT=FMTorder3) var1_round
	write(var2_filename,FMT=FMTorder3) var2_round
	subcell = trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

        open(filechannel1,file=path3//trim(subcell)//".dat",position="append")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)

end if


return


end subroutine addState



end module addCells5
