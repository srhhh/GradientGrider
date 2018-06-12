
module addCells4
implicit none


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      DIVYUP FUNCTION (QuickSort Alogrithm by Tony Hoare)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      INPUT:   real :: var1,var2,var3          "variablies used to define the 'to-be-divided' cells/subcells"
!               integer ::  order               "the order (N-1) of the 'to-be-divided' cell/subcell"
!               real :: scaling1                "the overcrowded cell/subcell will be divided into 'scaling1' parts along
!                       scaling2                    var1, 'scaling2' parts along 'var2', and 'scaling3' parts along 'var3'.
!                       scaling3                    the 'to-be-formed' cell/subcells have order N"
!               real :: resolution              "the number of subcells that the overcrowded cell/subcell will be divided into
!                                                resolution should be 'scaling1*scaling2*scaling3' but it is 
!                                                   'scaling1*scaling2' for now
!               integer :: indexN               "the index of the first 'to-be-formed' cell/subcells (order N) in 
!                                                   the counter array (see counterN)"
!               integer :: lengthN              "demision of array counterN (order N)"
!               integer :: overcrowdN           "the criteria of a cell be 'overcrowded' at order (N-1)"
!
!      OUTPUT:  integer :: counterN             "an array of integers, which are the number of frames in the 'to-be-formed' 
!                                                   cell/subcells (order N)"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      This function grids an overcrowded cell/subcell into higher order subcells.
!      This function triggers only when the size of a cell/subcell reaches the overcrowd limit for the first time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! RS: The variable "resolution" seems redundant
!                               KF: yes, although I don't want to have to
!                               calculate it again each time
! RS: OvercrowdN really should be overcrowd(N-1)
!                               KF: you got me!!!
!                               ... labelling it 'overcrowdN' is a little
!                               misleading, it is merely the size of the array
! RS: The trailing N in the varaible seems redundant and can be confusing
!                               KF: I use it to stress that the counter is
!                               dependent on the order, N, although it is
!                               certainly ambiguous here

recursive subroutine divyUp(var1,var2,var3,order,&
                            scaling1,scaling2,scaling3,resolution,&
                            indexN,counterN,lengthN,overcrowdN)
use f1_parameters
use f1_functions
implicit none
integer, intent(in) :: order,indexN,lengthN,overcrowdN
integer, dimension(lengthN), intent(out) :: counterN
integer :: i,j,k,l,population,last_marker,Nmarkers,index_order,counter_index
integer :: var1_NINT,var2_NINT
integer :: scaling1,scaling2,scaling3,resolution
real,intent(in) :: var1,var2,var3
real :: gap1, gap2, gap3, var1_new, var2_new
real,allocatable :: coords(:,:),vals(:,:)
integer, dimension(scaling1-1) :: grid1
integer, dimension(scaling2-1) :: grid2
integer, allocatable :: indexer(:,:)
character(50) :: subcell
character(1) :: descriptor1
character(10) :: descriptor2,descriptor3,descriptor4,descriptor5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          FRAME RETRIEVAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!First, get the filename of the cell to be divided

!If we are in a parent cell (order = 0) then we don't need any decimal places 
!Floor to the nearby integer
if (order == 0) then
        var1_new = anint(var1-0.5)
        var2_new = anint(var2-0.5)

        !We will make use of order (or zero) decimal places
        write(descriptor1,FMT="(I1)") order

!If we are in a subcell then we need decimal places
!Floor to the nearby multiple of .25*0.1**(order-1)
else
        var1_new = anint((scaling1_0*scaling1_1**(order-1))*var1-0.5)/(scaling1_0*scaling1_1**(order-1))
        var2_new = anint((scaling2_0*scaling2_1**(order-1))*var2-0.5)/(scaling2_0*scaling2_1**(order-1))

        !We will make use of order+1 decimal places
        ! ex. order=1 --> .00, .25, .50, .75, etc.
        !     order=2 --> .025, .050, .075, 0.100, .125, etc.
        write(descriptor1,FMT="(I1)") order+1
end if

!The number of decimal places for the cell depends on the order (as set above)
write(descriptor2,FMT="(F9."//descriptor1//")") var1_new
write(descriptor3,FMT="(F9."//descriptor1//")") var2_new

!We don't know how many decimal places will be used (may be double digit or
!single digit) so we need to "adjustl" (remove leading zeroes) and trim (remove
!trailing whitespace); this format is also followed in addState

subcell = trim(adjustl(descriptor2))//"_"//trim(adjustl(descriptor3))

!Open up the file, read the variables, coordinates, and gradients
!vals    ---  holds the criteria of sorting (three for now)
!coords  ---  holds the xyz coordinates and gradients (6*Natoms)
!        ---  coords will be sorted according to vals
!indexer ---  holds the indexes of coords
!        ---  indexer will be sorted according to vals
open(filechannel1,file=trim(path3)//trim(subcell)//".dat")
allocate(coords(overcrowdN,6*Natoms))
allocate(vals(overcrowdN+resolution-1,Nvar))
allocate(indexer(overcrowdN+resolution-1,1))
do j=1, overcrowdN
        read(filechannel1,FMT=FMT1,advance="no",iostat=k) (vals(j,i),i=1,Nvar)
        read(filechannel1,FMT=FMT2) (coords(j,i),i=1,6*Natoms)
        indexer(j,1) = j
end do
close(filechannel1)

! RS: Can you just use 'filechannel1' like this without defining it??
! RS: Could we avoid using filechannel'1' before filechannel is used? what does '1' refer to?
! RS: I noticed you used filechannel1 three times later. In this case, what differences do
! RS:    filechannel1 and '70' make?
!                               KF: it is in f1_parameters (73, I believe?)
! RS: Let's talk about naming tomorrow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          FRAME SORTING/GRIDING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!We need to know what the spacing is of the subcells
!This is inversely proportional to the scaling
! spacing1 and spacing 2 are the spacing of the parent-level grid (order = 0)
gap1 = spacing1/(scaling1_0*scaling1_1**(order))
gap2 = spacing2/(scaling2_0*scaling2_1**(order))

!Add markers to know when a subcell begins/ends
do i = 1, resolution-1
        counter_index = overcrowdN + i
        indexer(counter_index,1) = counter_index
        vals(counter_index,:) = (/ var1_new+(i/scaling2)*gap1,&
                                   var2_new+modulo(i,scaling2)*gap2,&
                                   0.0                                  /)
end do


!The filenames of the subcells will be split into two parts:
!  the prefix is of form "xx." (the integer part)
!  the suffix is of form ".yyyyy" (the decimal part)
!Clearly, the integer part does not depend on the subcell
!but on the parent cell, so it is initialized here

! RS: I understand this works for order > 0, what about order = 0?
! RS: will the value round like ANINT or just trimed like AINT?
!       KF: case 1: order == 0
!               var1_new = xx.00
!           case 2: order > 0
!               var1_new = xx.yyyy
!           in both cases, the decimal portion gets truncated
!               descriptor2 = 000000xx.
!                               I believe this works
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

!Sort the indexed frames by the first variable (into columns);
!Because we also added markers that have values corresponding to the gridlines,
!these markers will be sorted according to whether frames are before or after
!them. One pair of markers (p1,p2) would indicate that frames [p1,p2-1] are in
!the same column.
! RS: The comments do not seem to carry the same information as what we talked
! RS: By "one pair of markers", do you mean one 2D marker or two 1D markers?
! RS: What do you mean by "column"?
! RS: going to read qsort2 now, hopefully I can figure it out myself
!               KF: the cell is traversed one-dimensionally (counter_index)
!                   When a marker p2 is reached, we know we have reached the end
!                   of a subcell (and consequently the beginning of the next
!                   subcell); any frames in the range (p1,p2) are in the subcell
!                   corresponding to marker p2
!Qsort2 sorts both vals and indexer; indexer can then be used to access coords

! RS: I am very confused by this -- don't you add the markers in before you sort var1?
!               KF: I add the markers in before this point
!                       yeah... maybe we should talk about this
call qsort2(vals,indexer,overcrowdN+resolution-1,Nvar,1,overcrowdN+resolution-1,1)

!A marker can be identified by checking whether indexer(p) > overcrowd
!If so, then it is a marker. A pair of markers (p1,p2) would indicate a column
!but these need to be sorted further by var2 into cells. These are sorted the
!same way.
last_marker = 1
Nmarkers = 0
do counter_index = 1, overcrowdN+resolution-1
        if (indexer(counter_index,1) > overcrowdN) then
        Nmarkers = Nmarkers + 1
        if (Nmarkers == scaling2) then
                call qsort2(vals,indexer,overcrowdN+resolution-1,Nvar,&
                           last_marker,counter_index-1,2)
                last_marker = counter_index
                Nmarkers = 0
        end if
        end if

end do

!Now that vals and indexer are fully sorted, we need to identify subcells
!Any pair of markers (p1,p2) such that (p2-p1 > 1) is a subcell with a frame in
!it; these pairs are ordered in such a way so that the subcell (i,j) that
!corresponds with pair n is:
!       i = n / scaling2
!       j = n % scaling2
!And, of course, i and j start at zero.
last_marker = 1
Nmarkers = 0
do counter_index = 1, overcrowdN+resolution-1
        if (indexer(counter_index,1) > overcrowdN) then
        population = counter_index - last_marker
        if (population > 0) then
                i = (Nmarkers/scaling2)
                j = modulo(Nmarkers,scaling2)

                !var1_INT and var2_INT are only decimal part of the PARENT cell
                !thus, we have to add on .25*.1**order to get the decimal
                !portion of the CHILD subcell
                write(descriptor4,FMT="(I"//trim(descriptor1)//"."//trim(descriptor1)//")") &
                        var1_NINT + i*(100)/scaling1_0
                write(descriptor5,FMT="(I"//trim(descriptor1)//"."//trim(descriptor1)//")") &
                        var2_NINT + j*(100)/scaling2_0
                descriptor4 = trim(descriptor2)//trim(adjustl(descriptor4))
                descriptor5 = trim(descriptor3)//trim(adjustl(descriptor5))
                subcell = trim(descriptor4)//"_"//trim(descriptor5)

                !We write all of the frames onto the higher order subcell;
                !last_marker represents the first frame of the subcell
                !whereas counter_index is the position of the marker that ends
                !the subcell; thus, one is subctracted from it
                open(filechannel1,file=trim(path3)//trim(subcell)//".dat",status="new")
                do k = last_marker, counter_index-1
                        write(filechannel1,FMT=FMT1,advance="no") (vals(k,l),l=1,Nvar)
                        write(filechannel1,FMT=FMT2)(coords(indexer(k,1),l),l=1,6*Natoms)                       
                end do 
                close(filechannel1)

                !And of course, we need to add this onto counterN
                !The indexing is exactly as in addState
                index_order = resolution*indexN + scaling1*j + i
                counterN(index_order) = population
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
population = overcrowdN + resolution - 1 - last_marker
if (population > 0) then
        i = scaling1 - 1
        j = scaling2 - 1

        write(descriptor4,FMT="(I"//trim(descriptor1)//"."//trim(descriptor1)//")") &
                                var1_NINT + i*(100)/scaling1_0
        write(descriptor5,FMT="(I"//trim(descriptor1)//"."//trim(descriptor1)//")") &
                                var2_NINT + j*(100)/scaling2_0
        descriptor4 = trim(descriptor2)//trim(adjustl(descriptor4))
        descriptor5 = trim(descriptor3)//trim(adjustl(descriptor5))
        subcell = trim(descriptor4)//"_"//trim(descriptor5)

        open(filechannel1,file=trim(path3)//trim(subcell)//".dat",status="new")
        do k = last_marker, overcrowdN+resolution-1
                write(filechannel1,FMT=FMT1,advance="no") (vals(k,l),l=1,Nvar)
                write(filechannel1,FMT=FMT2)(coords(indexer(k,1),l),l=1,6*Natoms)                       
        end do 
        close(filechannel1)

        index_order = resolution*indexN + scaling1*j + i
        counterN(index_order) = population    
end if

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

        open(filechannel1,file=trim(path3)//trim(subcell)//".dat",position="append",status=trim(descriptor0))
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)

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
 
        open(filechannel1,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)
      
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

        open(filechannel1,file=trim(path3)//trim(subcell)//".dat",position="append",status=trim(descriptor0))
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)

        return

!And if overcrowded, then call divyUp as expected
else if (population == overcrowd1) then 
        key = key + key_start*header2
        counter1(indexer) = key

        call divyUp(vals(1),vals(2),0.0,order,&
                    scaling1_1,scaling2_1,0,resolution_1,&
                    header2,counter2,counter2_max,overcrowd1-1)
        header2 = header2 + 1

        open(filechannel1,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)

else
        open(filechannel1,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)

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
 
        open(filechannel1,file=trim(path3)//trim(subcell)//".dat",position="append",status=trim(descriptor0))
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)

        return

else if (population == overcrowd2) then 
        key = key + key_start*header3
        counter2(indexer) = key

        call divyUp(vals(1),vals(2),0.0,order,&
                    scaling1_1,scaling2_1,0,resolution_1,&
                    header3,counter3,counter3_max,overcrowd2-1)
        header3 = header3 + 1

        open(filechannel1,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)

else
        open(filechannel1,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)

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
 
        open(filechannel1,file=trim(path3)//trim(subcell)//".dat",position="append",status=trim(descriptor0))
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)

        return

!Because this is the end of the line, a notice is put up for no further divyUps
else if (population == overcrowd3) then 
        open(filechannel1,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)

        open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
        write(progresschannel,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(progresschannel,*) " FINAL LEVEL SUBCELL OVERCROWDED"
        write(progresschannel,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!"
        close(progresschannel)
else
        open(filechannel1,file=trim(path3)//trim(subcell)//".dat",position="append",status="old")
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT=FMT2)(coords(j),j=1,6*Natoms)
        close(filechannel1)

end if

return






end subroutine addState



end module addCells4
