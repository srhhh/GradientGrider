
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
logical :: state1
real,intent(in) :: var1,var2,var3
real :: gap1, gap2, gap3, var1_new, var2_new
real,allocatable :: coords(:,:),vals(:,:)
integer, dimension(scaling1) :: grid1
integer, dimension(scaling2) :: grid2
integer, allocatable :: indexer(:,:)
character(50) :: currentsubcell
character(50) :: subcell,line_data
character(10) :: descriptor1,descriptor2,descriptor3,descriptor4,descriptor5


!If we are in a parent cell then we don't need any decimal places 
!Floor to the nearby integer
if (order == 0) then
        var1_new = anint(var1-0.5)
        var2_new = anint(var2-0.5)
        write(descriptor1,FMT="(I1)") order

!If we are in a subcell then we need decimal places
!Floor to the nearby multiple of .25*0.1**order
else
        var1_new = anint((scaling1_0*scaling1_1**(order-1))*var1-0.5)/(scaling1_0*scaling1_1**(order-1))
        var2_new = anint((scaling2_0*scaling2_1**(order-1))*var2-0.5)/(scaling2_0*scaling2_1**(order-1))
        write(descriptor1,FMT="(I1)") order+1
end if

!The smaller the subcell, the smaller the gap between gridline
gap1 = 1.0/(scaling1_0*scaling1_1**(order))
gap2 = 1.0/(scaling2_0*scaling2_1**(order))

!We need to know the name of the file so we can read its
!variables and coordinates into a local array
write(descriptor2,FMT="(F9."//trim(adjustl(descriptor1))//")") var1_new
write(descriptor3,FMT="(F9."//trim(adjustl(descriptor1))//")") var2_new
descriptor2 = adjustl(descriptor2)
descriptor3 = adjustl(descriptor3)
subcell = trim(descriptor2)//"_"//trim(descriptor3)

!Open up the file,read the variables, coordinates, gradients
!vals    ---  holds the to-be-sorted variables (three for now)
!coords  ---  holds the xyz coordinates and gradients (6*Natoms)
!indexer ---  holds the to-be-sorted indexes to access coords
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

!The heading digits are used to later reconstruct a subcell file name
write(descriptor2,FMT="(F9.0)") var1_new-0.5
write(descriptor3,FMT="(F9.0)") var2_new-0.5

!We will needer order + 2 decimal places to represent a smaller subcell
!write(descriptor1,FMT="(I1)") order+2
write(descriptor1,FMT="(I1)") order+2

!Essentially this stores the digits of the number---aftering flooring--
!to the right of the decimal place
!In this PARTICULAR case, every increment along the variable axis
!requires two decimal places for the first scaling
!   e.g. 1. -->  1.00, 1.25, 1.50, 1.75
!And one additional decimal place for subsequent scalings
!   e.g. 1.00 --> 1.025, 1.050, 1.075, ... , 1.225
var1_NINT = nint((var1_new-floor(var1_new))*(10**(order+2)))
var2_NINT = nint((var2_new-floor(var2_new))*(10**(order+2)))
descriptor2 = adjustl(descriptor2)
descriptor3 = adjustl(descriptor3)

!Sort the indexed frames by the first variable (into columns); then grid it
!This sorts both vals and indexer; indexer can then be used to access coords
call qsort(vals,indexer,overcrowdN,Nvar,1,overcrowdN,1)
call grider(grid1,vals,gap1,var1_new,scaling1,overcrowdN,Nvar,1,overcrowdN,1,scaling1,1)

index1_1 = grid1(1)
do i = 1, scaling1
        !We always accept these grid elements (cells) in pairs of indexed frames [p1,p2)
        !where indexed frame p2 is not in the cell i
        !For the last element (cell) the last indexed frame must be in it (Nstate)
        !For the first variable, we can think of these are "columns"
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
        call qsort(vals,indexer,overcrowdN,Nvar,index1_1,index1_2-1,2)
        call grider(grid2,vals,gap2,var2_new,scaling2,overcrowdN,Nvar,&
                        index1_1,index1_2-1,1,scaling2,2)
        index2_1 = grid2(1)
        do j = 1, scaling2
                !For the second variable, we can think of these as "cells"
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

                !Write these frames into their respective cell in the grid
                !The array index needs to start at 1 so a subtraction is
                !involved
                !In this PARTICULAR case, scaling_0 = 4 so 100/scaling_0 = 25
                !Each 'increment' along a variable-axis should thus increment
                !The prospective filename by 25; other parameter choices
                !will break this; see comment at top of do loops
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
                index_order = resolution*indexN + scaling1*(j-1) + i-1
                counterN(index_order) = index2_2 - index2_1

                !This is in the rare case that all frames are in one subcell
                !We can also exit the loop if this is the last frame
                if (index2_2 == overcrowdN+1) then
                if (index2_1 == 1) counterN(index_order) = counterN(index_order)-1
                exit
                end if

                !The next cell pair [p2,p3) starts at the end of [p1,p2)
                index2_1 = index2_2
        end do

        !The next cell pair [p2,p3) starts at the end of [p1,p2)
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
!! RS: Add comment:
!! RS: key_start in modulo(key,key_start) is defined in f1_parameters.f90
!! RS: key_start is defined as the max number of frames that the counter array will keep track
!! RS:   beyond which the digits are used as pointer to the children array
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
        !! RS: Add comment:
        !! RS: header1 = 1 from getCell.f90
        key = key + key_start*header1
        counter0(indexer) = key

        open(72,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(72,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(72,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(72)

        !For more information on how this works, see the subroutine above
        call divyUp(vals(1),vals(2),0.0,order,&
                    scaling1_0,scaling2_0,0,resolution_0,&
                    header1,counter1,counter1_max,overcrowd0)
       
        !Incrementing header insures that the position granted is unique
        header1 = header1 + 1

        !I do not call 'return' here because there is a small chance the 
        !subcell needs further subdividing (all frames in a cell were placed
        !in a single subcell); for divyUp to be called again, another
        !loop needs to go around
!! RS: Don't you need a "return" after this to quit from going deeper and add to order??
!! RS: This "else" is useless 
else
end if





!Now, we go deeper; the 'depth' is represented by the order
order = order + 1

!! RS: Add comment:
!! RS: scaling1_0 = scaling2_0 = 4
!! RS: scaling1_0*scaling2_0 are the number of subcells that is divided from their parent cell 

!This recovers the variable floored to the nearby multiple of 0.25
!    e.g.  1.2445 -> 1.00     or    1.4425 -> 1.25
var1 = anint(vals(1)*scaling1_0-0.5)/scaling1_0
var2 = anint(vals(2)*scaling2_0-0.5)/scaling2_0

!This recovers what integer var1, var2 correspond to
!    e.g.  1.00 -> 0     or    1.75 -> 3
var1_new = modulo(nint(var1*scaling1_0),scaling1_0)
var2_new = modulo(nint(var2*scaling2_0),scaling2_0)

!! RS: I believe your algorithm works, but how about:
!! RS: INT(MODULO(vals(1),spacing_0)/spacing_1)
!! RS: If A and P are of type REAL:
!! RS:     MODULO(A,P) has the value of A - FLOOR (A / P) * P.

!Here is the unqiue indexing method of this subroutine; the index to counter1
!comes from the key acquired through counter0; it is represented by the digits
!larger than key_start
!! RS: Add comment
!! RS: resolution_0 = scaling1_0*scaling2_0

!! RS: I am not sure if this indexer line will work. Please talk to me tomorrow
indexer = resolution_0*(key/key_start) + scaling2_0*var2_new+var1_new
!! RS: replace (key/key_start) by INT(key/key_start)
!! RS: shouldn't "scaling2_0*var2_new+var1_new" be "scaling2_0*(var2_new-1)+var1_new"?
!! RS: like you did for the parent cell?

!! RS: the following comment is not accurate
!! RS: it is not necessarily "making" -- I think "locating" could be a more accurate word
!And then simply make the name of the new subcell
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

        open(72,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(72,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(72,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(72)

        call divyUp(vals(1),vals(2),0.0,order,&
                    scaling1_1,scaling2_1,0,resolution_1,&
                    header2,counter2,counter2_max,overcrowd1)
        header2 = header2 + 1
!! RS: I think a "return" might be necessary
else
!! RS: what does it make you want to keep writing to the overcrowded file (but not at the parent level)?
        open(72,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(72,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(72,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(72)

end if



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

!! RS: Again, I think this can be replaced by an easier algorithm
!! RS: Please refer to line 324 to 327

indexer = resolution_1*(key/key_start) + scaling2_1*var2_new+var1_new
!! RS: "var2_new" should be "var2_new-1"

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

        open(72,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(72,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(72,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(72)

open(80,file=trim(path4)//trim(progressfile),position="append")
write(80,*) "   Subdividing: ", trim(subcell)
close(80)
!! RS: Why is this printing statement only showing up when order=2?
        call divyUp(vals(1),vals(2),0.0,order,&
                    scaling1_1,scaling2_1,0,resolution_1,&
                    header3,counter3,counter3_max,overcrowd2)
        header3 = header3 + 1
!! RS: return?
else
        open(72,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(72,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(72,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(72)

end if


!! Similar comments are previous level

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
