
module addCells4
implicit none


contains


!The following grids the cell into smaller cells
!This ASSUMES the cell is NOT already divided up

recursive subroutine divyUp(var1,var2,var3,order,&
                            scaling1,scaling2,scaling3,resolution,&
                            indexN,counterN,lengthN,overcrowdN)
use f1_parameters
use f1_functions
implicit none
integer, intent(in) :: order,indexN,lengthN,overcrowdN
integer, dimension(lengthN), intent(out) :: counterN
integer :: i,j,k,l,index_order,index1_1,index1_2,index2_1,index2_2
integer :: var1_NINT,var2_NINT
integer :: scaling1,scaling2,scaling3,resolution
logical :: state1,flag1
real,intent(in) :: var1,var2,var3
real :: gap1, gap2, gap3, var1_new, var2_new
real,allocatable :: coords(:,:),vals(:,:)
integer, dimension(scaling1) :: grid1
integer, dimension(scaling2) :: grid2
integer, allocatable :: indexer(:,:)
character(50) :: currentsubcell
character(50) :: subcell,line_data
character(1) :: descriptor1
character(10) :: descriptor2,descriptor3,descriptor4,descriptor5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          FRAME RETRIEVAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!First, get the filename of the cell to be divided

!If we are in a parent cell then we don't need any decimal places 
!Floor to the nearby integer
if (order == 0) then
        var1_new = anint(var1-0.5)
        var2_new = anint(var2-0.5)

        !We will make use of order (or zero) decimal places
        write(descriptor1,FMT="(I1)") order

!If we are in a subcell then we need decimal places
!Floor to the nearby multiple of .25*0.1**order
else
        var1_new = anint((scaling1_0*scaling1_1**(order-1))*var1-0.5)/(scaling1_0*scaling1_1**(order-1))
        var2_new = anint((scaling2_0*scaling2_1**(order-1))*var2-0.5)/(scaling2_0*scaling2_1**(order-1))

        !We will make use of order+1 decimal places
        ! ex. order=1 --> .00, .25, .50, .75, etc.
        !     order=2 --> .025, .125, .350, .500 etc.
        write(descriptor1,FMT="(I1)") order+1
end if

!The number of decimal places for the cell depends on the order (as set above)
write(descriptor2,FMT="(F9."//descriptor1//")") var1_new
write(descriptor3,FMT="(F9."//descriptor1//")") var2_new

!We don't know how many decimal places will be used (may be double digit or
!single digit) so we need to "adjustl" (remove trailing zeroes) and trim (remove
!trailing whitespace); this format is also followed in addState
subcell = trim(adjustl(descriptor2))//"_"//trim(adjustl(descriptor3))

!Open up the file,read the variables, coordinates, gradients
!vals    ---  holds the criteria of sorting (three for now)
!coords  ---  holds the xyz coordinates and gradients (6*Natoms)
!        ---  coords will be sorted according to vals
!indexer ---  holds the indexes of coords
!        ---  indexer will be sorted according to vals
open(72,file=trim(path3)//trim(subcell)//".dat")
allocate(vals(overcrowdN,Nvar))
allocate(coords(overcrowdN,6*Natoms))
allocate(indexer(overcrowdN,1))
do j=1, overcrowdN
        read(72,FMT=FMT1,advance="no",iostat=k) (vals(j,i),i=1,Nvar)
        read(72,FMT=FMT2) (coords(j,i),i=1,6*Natoms)
        indexer(j,1) = j
end do
close(72)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          FRAME SORTING/GRIDING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!We need to know what the spacing is of the subcells
!This is inversely proportional to the scaling
gap1 = spacing1/(scaling1_0*scaling1_1**(order))
gap2 = spacing2/(scaling2_0*scaling2_1**(order))

!The filenames of the subcells will be split into two parts:
!  the prefix is of form "xx." (the integer part)
!  the suffix is of form ".yyyyy" (the decimal part)
!Clearly, the integer part does not depend on the subcell
!but on the parent cell, so it is initialized here

! RS: You var1_new is already rounded 
! RS: "var1_new = anint(var1-0.5)"
! RS: Why are you doing it again?            KF: resolved (code error)

! RS: Wait... Did you fix it? shouldn't you just delete the following two lines
! RS: and move the "adjust1" to line 65? Am I missing something?
write(descriptor2,FMT="(F9.0)") var1_new - 0.5
write(descriptor3,FMT="(F9.0)") var2_new - 0.5

!And remove any leading zeroes
descriptor2 = adjustl(descriptor2)
descriptor3 = adjustl(descriptor3)

!We will need (order + 2) decimal places to represent the first level subcell from parent cell
!   ex. 4. --> 4.00, 4.25, 4.50, 4.75 (decimal places 0 -> 2, order 0 -> 1)
!then decimal places encrease by 1 per level (order) of subcell
!   ex. 4.25 --> 4.250, 4.275, 4.300, 4.325, 4.350, 4.375, 4.400, 4.425, 4.450, 4.475
!      (decimal places 2 -> 3, order 1 -> 2)             

write(descriptor1,FMT="(I1)") order+2

!Essentially this stores the digits of the number---aftering flooring--
!to the right of the decimal place
! ex. 4.25 --> 250   or   4. --> 00   or 3.2275 --> 22750
var1_NINT = nint((var1_new-floor(var1_new))*(10**(order+2)))
var2_NINT = nint((var2_new-floor(var2_new))*(10**(order+2)))

!Sort the indexed frames by the first variable (into columns); then grid it
!This sorts both vals and indexer; indexer can then be used to access coords

call qsort2(vals,indexer,overcrowdN,Nvar,1,overcrowdN,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   QSORT2 is a subroutine of F1_FUNCTIONS.F90
!   Takes vals (1st argument) with dimensions overcrowdN, Nvar (3rd, 4th argument)
!   And indexer (2nd argument) with dimensions overcrowdN, 1 (3rd argument)
!   And sorts them according to the 1-th column of vals (7th argument)
!   But only from indexes 1 to overcrowdN (5th, 6th argument)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call grider(grid1,vals,gap1,var1_new,scaling1,overcrowdN,Nvar,1,overcrowdN,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   GRIDER is a subroutine of F1_FUNCTIONS.F90
!   Outputs grid1 (1st argument) with dimension scaling1 (5th argument)
!   from vals (2nd argument) with dimensions overcrowdN, Nvar (6ht, 7th arguments)
!   The value of the first gridline (i=1) is at var1_new (4th argument)
!   And subsequent gridlines are multiples of gap1 (3rd argument)
!   This stops after scaling1 (5th argument) gridlines have been made
!   This only does gridings for indexes 1 to overcrowdN (8th, 9th arguments) of vals
!
!   The value of position i in grid1 corresponds to the index minus one of vals
!   with the largest value in the 1-th column (10th argument) but less than
!   the value of the gridline i. Let's say we get this output:
!   ex.   grid1 = [ 1, 10, 12, 20, 100, 200, 201, 201, 201, 501]
!                     for 500 frames in subcell 4.25_1.00 (grid spacing = 0.25)
!         * 4.250 would be the first gridline (i=1)
!         * 4.475 would be the last gridline (i=10)
!         * grid1(1) = 1 indicates that 1 - 1 = 0 is the index of the largest
!                      value lower than 4.250; this is out-of-bounds because
!                      no value of this cell should be less than 4.25
!         * grid1(2) = 10 indicates that 10 - 1 = 9 is the index of the largest
!                      value lower than 4.275; thus, 9 frames are in the first
!                      'cell' or 'column'
!         * grid1(9) = grid1(8) = grid1(7) = 201 indicates that index 200
!                      of vals is the largest value of all gridlines
!                      i=7,8,9 (or 4.400,4.425,4.450) LOWER than 4.400.
!                      This means cells/column 8 and 9 (4.400,4.425) are empty
!         * grid1(10) = 501 indicates that 501 - 1 = 500 is the index of the
!                      largest value lower than 4.475. Because this is the last
!                      frame, and it has a lower value than 4.475, cell/column
!                      10 is thus empty.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Great, this is very clear. Thank you

index1_1 = grid1(1)
do i = 1, scaling1
        !For grid1, we can think of indexes as "columns"
        !Column i is mapped to pair [p1,p2); p1=grid1(i), p2=grid1(i+1)
        !where p1 and p2 are indexes to vals, and
        !       vals(p1) is the FIRST frame in column i
        !       vals(p2) is the LAST frame in column i
        !       and everything between them is also in column i
        !For the last column, p2= overcrowdN+1 to indicate that
        !the remaining frames must be in it
        !if they were not already placed in a column
        if (i == scaling1) then
                index1_2 = overcrowdN+1
        else
                index1_2 = grid1(i+1)
        end if

        !If the column is empty, don't even bother
        if (index1_2 == index1_1) then
                index1_1 = index1_2
                cycle
        end if

        !Sort the indexed states in this grid element (into cells); then grid it
        call qsort2(vals,indexer,overcrowdN,Nvar,index1_1,index1_2-1,2)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   The key difference here is that we only want to sort values of the
        !   current column.
        !   Values are in the same column if their indexes are in range [p1,p2)
        !   p2 is not included in the cell so there is a minus one subtraction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call grider(grid2,vals,gap2,var2_new,scaling2,overcrowdN,Nvar,&
                        index1_1,index1_2-1,2)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   The key difference here is that we only want to grid based on values
        !   of the current column
        !   Similar index-choosing as in QSORT2 above
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        index2_1 = grid2(1)
        do j = 1, scaling2
                !For grid2, we can think of indexes as "cells" of the column
                !Cell j is mapped to pair [p1,p2); p1=grid1(j),p2=grid2(j+1)
                !where p1 and p2 are indexes to vals, and
                !       vals(p1) is the FIRST frame in cell j, column i
                !       vals(p2) is the LAST frame in cell j, column i
                !       and everything between them is also in cell j
                !For the last cell, p2 = index1_2 to indicate that
                !the remaining frames must be in it
                !if they were not already placed in a cell
                if (j == scaling2) then
                        index2_2 = index1_2
                else
                        index2_2 = grid2(j+1)
                end if                
        
                !If there's nothing in the cell, don't even bother
                if (index2_1 == index2_2) then
                        index2_1 = index2_2
                        cycle
                end if

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !             FRAME WRITING
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !Write these frames into their respective cell in the grid
                !The array index needs to start at 1 so a subtraction is
                !involved
                !As consequence of having a resolution of 16 (scaling of 4)
                !we require two decimal places so different subcells will have
                !filenames differing by a multiple of 25
                !Because subcells are children of the parent cell, their
                !numbering starts from var_INT of their parent
                write(descriptor4,FMT="(I"//trim(descriptor1)//"."//trim(descriptor1)//")") &
                        var1_NINT + (i-1)*(100)/scaling1_0
                write(descriptor5,FMT="(I"//trim(descriptor1)//"."//trim(descriptor1)//")") &
                        var2_NINT + (j-1)*(100)/scaling2_0
                descriptor4 = trim(descriptor2)//trim(adjustl(descriptor4))
                descriptor5 = trim(descriptor3)//trim(adjustl(descriptor5))
                subcell = trim(descriptor4)//"_"//trim(descriptor5)

                !We write all of the frames onto the higher order subcell
                open(72,file=trim(path3)//trim(subcell)//".dat",status="new")
                do k = index2_1, index2_2-1
                        write(72,FMT=FMT1,advance="no") (vals(k,l),l=1,Nvar)
                        write(72,FMT=FMT2)(coords(indexer(k,1),l),l=1,6*Natoms)
                end do
                close(72)

                !And we also want to keep track of how many frames are in this
                !new subcell; indexing is exactly as in addState
                !A minus one is involved because array indexes start at 1
                index_order = resolution*indexN + scaling1*(j-1) + i-1
                counterN(index_order) = index2_2 - index2_1

                !The next cell pair [p2,p3) starts at the end of [p1,p2)
                index2_1 = index2_2
        end do

        !The next column pair [p2,p3) starts at the end of [p1,p2)
        index1_1 = index1_2
end do

end subroutine divyUp





!addState adds a state/frame (its coordinates and gradients stored in coords)
!to all subcells corresponding to var1, var2, var3
!where var1, var2, var3 are stored in vals

subroutine addState(vals,coords,&
                header1,header2,header3,&
                counter0,counter1,counter2,counter3)
use f1_parameters
implicit none
integer :: indexer,i,j,population,key,overcrowd
integer,intent(out) :: header1
integer,intent(out) :: header2
integer,intent(out) :: header3
integer,dimension(counter0_max),intent(out) :: counter0
integer,dimension(counter1_max),intent(out) :: counter1
integer,dimension(counter2_max),intent(out) :: counter2
integer,dimension(counter3_max),intent(out) :: counter3
logical :: flag1
real :: var1, var2
integer :: var1_new,var2_new,order
real, dimension(Nvar), intent(in) :: vals
real, dimension(6*Natoms), intent(in) :: coords
character(50) :: descriptor0, descriptor1, descriptor2
character(9) :: descriptor3, descriptor4
character(50) :: subcell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               ORDER 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!First we get the integer out of the way, this is the 
!parent-level griding. Because this level of granularity is
!pointless, eventually we stop writing to the parent-level
!cells altogether.
!But it is still neceesary to get it out of the way
write(descriptor3,FMT="(F9.0)") vals(1)-0.5
write(descriptor4,FMT="(F9.0)") vals(2)-0.5

!We name the file after these digits (decimals included!)
!The parent-level is of order 0
subcell = trim(adjustl(descriptor3))//"_"//trim(adjustl(descriptor4))
order = 0

!Simply read off the integer
var1_new = nint(vals(1)-0.5)
var2_new = nint(vals(2)-0.5)

!Subtract one to not go out of bounds (potentially)
!And multiply by the 'grid length' to uniquely
!order each index pair (r1,r2)
indexer = anint(max_var1)*(var2_new-1) + var1_new

!The key is incremented to signify a frame is added
!(even though it may ultimately not be added)
key = counter0(indexer) + 1
population = modulo(key,key_start)

!Constantly having to write to the file (which gets big!) is time-consuming
!So after it is filled, it is no longer written to
!Additonally, counter0 could get 'overfilled' so it is not updated
!after it is overcrowded
if (population < overcrowd0) then
        counter0(indexer) = key

        !If this is the first time in the cell, this file has to be made
        if (population == 1) then
                descriptor0 = "new"
        else
                descriptor0 = "old"
        end if

        open(72,file=trim(path3)//trim(subcell)//".dat",position="append",status=trim(descriptor0))
        write(72,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(72,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(72)

        !Consequently, this is the end of the line. There is no deeper
        !subcell so the next frame can be added in now
        return

!Of course, if it is overcrowded, then we need to divyUp the cell
!This makes the smaller subcells out of the cell and additionally signals
!that next time, a frame can be added to a deeper subcell
else if (population == overcrowd0) then

        !To access a subcell of higher order ('deeper subcell')
        !We grant the cell a position in counter1 for each of its
        !potential cells (100 in this particular case)
        key = key + key_start*header1
        counter0(indexer) = key

        !For more information on how this works, see the subroutine above
        call divyUp(vals(1),vals(2),0.0,order,&
                    scaling1_0,scaling2_0,0,resolution_0,&
                    header1,counter1,counter1_max,overcrowd0-1)
 
        open(72,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(72,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(72,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(72)
      
        !Incrementing header insures that the position granted is unique
        header1 = header1 + 1

end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               ORDER 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Now, we go deeper; the 'depth' is represented by the order
order = order + 1

!This recovers the variable floored to the nearby multiple of 0.25
!    e.g.  1.2445 -> 1.00     or    1.4425 -> 1.25
var1 = anint(vals(1)*scaling1_0-0.5)/scaling1_0
var2 = anint(vals(2)*scaling2_0-0.5)/scaling2_0

!This recovers what integer var1, var2 correspond to
!    e.g.  1.00 -> 0     or    1.75 -> 3
var1_new = modulo(nint(var1*scaling1_0),scaling1_0)
var2_new = modulo(nint(var2*scaling2_0),scaling2_0)

!Here is the unqiue indexing method of this subroutine; the index to counter1
!comes from the key acquired through counter0; it is represented by the digits
!larger than key_start
indexer = resolution_0*(key/key_start) + scaling2_0*var2_new+var1_new

!And then make the name of the new subcell
write(descriptor3,FMT="(F9.2)") var1
write(descriptor4,FMT="(F9.2)") var2
subcell = trim(adjustl(descriptor3))//"_"//trim(adjustl(descriptor4))

!And just like in the previous step, we increment by one
!In this case, though, each frame is relevant so they are always stored
key = counter1(indexer) + 1
population = modulo(key,key_start)
counter1(indexer) = key

if (population < overcrowd1) then
        if (population == 1) then
                descriptor0 = "new"
        else
                descriptor0 = "old"
        end if

        open(72,file=trim(path3)//trim(subcell)//".dat",position="append",status=trim(descriptor0))
        write(72,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(72,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(72)

        return

!And if overcrowded, then call divyUp as expected
else if (population == overcrowd1) then 
        key = key + key_start*header2
        counter1(indexer) = key

        call divyUp(vals(1),vals(2),0.0,order,&
                    scaling1_1,scaling2_1,0,resolution_1,&
                    header2,counter2,counter2_max,overcrowd1-1)
        header2 = header2 + 1

        open(72,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(72,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(72,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(72)

else
        open(72,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(72,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(72,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(72)

end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               ORDER 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!And we go deeper once more
!Most of the steps are the same
!Remark: scaling is changed between order=1 and order=2 so that is why
!There is more than one scaling element;
!But all variables still act the same
order = order + 1
var1 = anint((scaling1_0*scaling1_1)*vals(1)-0.5)/(scaling1_0*scaling1_1)
var2 = anint((scaling2_0*scaling2_1)*vals(2)-0.5)/(scaling2_0*scaling2_1)

var1_new = modulo(nint(scaling1_0*scaling1_1*var1),scaling1_1)
var2_new = modulo(nint(scaling2_0*scaling2_1*var2),scaling2_1)

indexer = resolution_1*(key/key_start) + scaling2_1*var2_new+var1_new

write(descriptor3,FMT="(F9.3)") var1
write(descriptor4,FMT="(F9.3)") var2
subcell = trim(adjustl(descriptor3))//"_"//trim(adjustl(descriptor4))

key = counter2(indexer) + 1
population = modulo(key,key_start)
counter2(indexer) = key

if (population < overcrowd2) then
        if (population == 1) then
                descriptor0 = "new"
        else
                descriptor0 = "old"
        end if
 
        open(72,file=trim(path3)//trim(subcell)//".dat",position="append",status=trim(descriptor0))
        write(72,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(72,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(72)

        return

else if (population == overcrowd2) then 
        key = key + key_start*header3
        counter2(indexer) = key

        call divyUp(vals(1),vals(2),0.0,order,&
                    scaling1_1,scaling2_1,0,resolution_1,&
                    header3,counter3,counter3_max,overcrowd2-1)
        header3 = header3 + 1

        open(72,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(72,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(72,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(72)

else
        open(72,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(72,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(72,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(72)

end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               ORDER 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!And we go deeper once more
!Same as last time
order = order + 1
var1 = anint((scaling1_0*scaling1_1**2)*vals(1)-0.5)/(scaling1_0*scaling1_1**2)
var2 = anint((scaling2_0*scaling2_1**2)*vals(2)-0.5)/(scaling2_0*scaling2_1**2)

var1_new = modulo(nint((scaling1_0*scaling1_1**2)*var1),scaling1_1)
var2_new = modulo(nint((scaling2_0*scaling2_1**2)*var2),scaling2_1)

indexer = resolution_1*(key/key_start) + scaling2_1*var2_new+var1_new

write(descriptor3,FMT="(F9.4)") var1
write(descriptor4,FMT="(F9.4)") var2
subcell = trim(adjustl(descriptor3))//"_"//trim(adjustl(descriptor4))

key = counter3(indexer) + 1
population = modulo(key,key_start)
counter3(indexer) = key

if (population < overcrowd3) then
        if (population == 1) then
                descriptor0 = "new"
        else
                descriptor0 = "old"
        end if
 
        open(72,file=trim(path3)//trim(subcell)//".dat",position="append",status=trim(descriptor0))
        write(72,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(72,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(72)

        return

!Because this is the end of the line, a notice is put up for no further divyUps
else if (population == overcrowd3) then 
        open(72,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(72,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(72,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(72)

        open(80,file=trim(path4)//trim(progressfile),position="append")
        write(80,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(80,*) " FINAL LEVEL SUBCELL OVERCROWDED"
        write(80,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!"
        close(80)
else
        open(72,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(72,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(72,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(72)

end if

return






end subroutine addState



end module addCells4
