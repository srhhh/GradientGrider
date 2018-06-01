
module addCells3
implicit none


contains


!The following grids the cell into smaller cells
!This ASSUMES the cell is NOT already divided up

recursive subroutine divyUp(var1,var2,var3,order,indexN,&
                                counterN,lengthN,overcrowdN)
use f1_parameters
use f1_functions
implicit none
integer, intent(in) :: order,indexN,lengthN,overcrowdN
integer, dimension(lengthN), intent(out) :: counterN
integer :: i,j,k,l,index_order,index1_1,index1_2,index2_1,index2_2
logical :: state1
real,intent(in) :: var1,var2,var3
real :: gap1, gap2, gap3, var1_new, var2_new
real,allocatable :: coords(:,:),vals(:,:)
integer, dimension(scaling1) :: grid1
integer, dimension(scaling2) :: grid2
integer, allocatable :: indexer(:,:)
character(50) :: currentsubcell
character(50) :: subcell,line_data
character(10) :: descriptor1,descriptor2,descriptor3,&
                descriptor4,descriptor5,descriptor6


!The size of the 'gaps' represents the length of the smaller subcells
!inside the larger (current!) subcell
gap1 = spacing1*(1.0/scaling1)**(order)
gap2 = spacing2*(1.0/scaling2)**(order)
gap3 = 1.0*(1.0/scaling3)**(order)

!This part, in a sense, extracts the current digit from vals(i)
!and stores in in vari_new;
!This part is not affected by reparametrizing scaling, I believe
var1_new = anint((scaling1**order)*var1-0.5)*(1.0/scaling1)**order
var2_new = anint((scaling2**order)*var2-0.5)*(1.0/scaling2)**order

!We also need to know the name of the file so we can read its
!variables and coordinates into a local array
write(descriptor1,FMT="(I1)") order
write(descriptor2,FMT="(F9."//trim(adjustl(descriptor1))//")") var1_new
write(descriptor3,FMT="(F9."//trim(adjustl(descriptor1))//")") var2_new
descriptor2 = adjustl(descriptor2)
descriptor3 = adjustl(descriptor3)
subcell = trim(descriptor2)//"_"//trim(descriptor3)

!Typically second-level subdivisions (where 6.42 --> 6.42x) don't happen
!often so I note this on the progressfile;
!third-level subdivisions have not been supported (inside addState)
if (order == 2) then
open(80,file=trim(path4)//trim(progressfile),position="append")
write(80,*) "   Subdividing: ", trim(subcell)
close(80)
end if

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

!Because we are subdividing, the gaps are now shorter
gap1 = gap1/scaling1
gap2 = gap2/scaling2

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
                write(descriptor4,FMT="(I1)") i-1
                write(descriptor5,FMT="(I1)") j-1
                descriptor4 = trim(descriptor2)//trim(adjustl(descriptor4))
                descriptor5 = trim(descriptor3)//trim(adjustl(descriptor5))
                subcell = trim(descriptor4)//"_"//trim(descriptor5)

                !We write all of the frames onto the next level subcell
                open(72,file=trim(path3)//trim(subcell)//".dat",status="new")
                do k = index2_1, index2_2-1
                        write(72,FMT=FMT1,advance="no") (vals(k,l),l=1,Nvar)
                        write(72,FMT=FMT2)(coords(indexer(k,1),l),l=1,6*Natoms)
                end do
                close(72)

                !And we also want to keep track of how many frames are in this
                !new subcell.
                index_order = resolution*indexN + scaling1*(j-1) + i-1
                counterN(index_order) = index2_2 - index2_1

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
integer,dimension(ceiling(max_var1*max_var2)),intent(out) :: counter0
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
character(50) :: subcell_old, subcell_new

!At first, we only want to look at the first order spacing
!So we truncate values by using ANINT (implicit) after
!subtracting by 0.5
!This hinges on the fact that the first-level spacing is 1 A
order = 0
write(descriptor1,FMT=FMT4) order
write(descriptor2,FMT=FMT4) order+1
write(descriptor3,FMT="(F9."//trim(descriptor1)//")") vals(1)-0.5
write(descriptor4,FMT="(F9."//trim(descriptor1)//")") vals(2)-0.5

!We name the file after these digits; assume population is zero
subcell_new = trim(adjustl(descriptor3))//"_"//trim(adjustl(descriptor4))
population = 0
overcrowd = overcrowd0

do

        !Check whether there is some subcell r1_r2.dat in path3/
        inquire(file=trim(path3)//trim(subcell_new)//".dat",exist=flag1)

        !If so, then we must append on to an existing file
        if (flag1) then
                descriptor0 = "old"

                !If this is the first time in the loop, then we need to
                !compute the index of r1_r2 in counter0.
                !This is simply 80*d2 + d1 (unless we change the parameters)
                !where d2 and d1 are the numbers .order. digits to the
                !right of the decimal palce
                if (order == 0) then
                        read(descriptor3,FMT="(I8)") var1_new
                        read(descriptor4,FMT="(I8)") var2_new
                        indexer = ceiling(max_var1)*(var2_new-1) + var1_new

                        !The number of frames is stored in counterN
                        !We increment the value to signify we are adding a frame
                        !And we measure the population as the last 5 digits
                        key = counter0(indexer) + 1
                        population = modulo(key,100000)
                        counter0(indexer) = key

                !Similar to the above case where order == 0.
                !However, we can access the index of r1_r2 in counterN by
                !simply taking the first three digits of the previous key
                !and adding on 10*d2 + d1 (unless we change the parameters)
                !That is to say, the key provided by counterN-1 of r1_r2 holds
                !information to accessing counterN of r1_r2
                else if (order == 1) then
                        read(descriptor3,FMT="(8x,(I1))") var1_new 
                        read(descriptor4,FMT="(8x,(I1))") var2_new
                        indexer = resolution*(key/100000) + scaling1*var2_new+var1_new
                        key = counter1(indexer) + 1
                        population = modulo(key,100000)
                        overcrowd = overcrowd1
                        counter1(indexer) = key

                !Exactly the same as above
                else if (order == 2) then
                        read(descriptor3,FMT="(8x,(I1))") var1_new 
                        read(descriptor4,FMT="(8x,(I1))") var2_new 
                        indexer = resolution*(key/100000) + scaling1*var2_new+var1_new
                        key = counter2(indexer) + 1
                        population = modulo(key,100000)
                        overcrowd = overcrowd2
                        counter2(indexer) = key

                !Exactly the same as above
                else if (order == 3) then
                        read(descriptor3,FMT="(8x,(I1))") var1_new 
                        read(descriptor4,FMT="(8x,(I1))") var2_new 
                        indexer = resolution*(key/100000) + scaling1*var2_new+var1_new
                        key = counter3(indexer) + 1
                        population = modulo(key,100000)
                        overcrowd = overcrowd3
                        counter3(indexer) = key

                !For now, we stop; we can go as deep as we want so long
                !as we have enough space for counterN arrays
                else
                        exit
                end if

        !If the file doesn't exist then we make it
        !This might be a problem if we don't have an overcrowded parent cell
        !but we prune out that possibility with a later conditional statement
        else
                descriptor0 = "new"

                !Very similar format as to the above statements
                !The only difference is that we define the new key here
                !And because we only add one frame, and because it cannot
                !be overcrowded, the key is always equal to 1
                if (order == 0) then
                        read(descriptor3,FMT="(I8)") var1_new
                        read(descriptor4,FMT="(I8)") var2_new
                        indexer = ceiling(max_var1)*(var2_new-1) + var1_new
                        key = 1
                        population = 1
                        counter0(indexer) = population

                else if (order == 1) then
                        read(descriptor3,FMT="(8x,(I1))") var1_new 
                        read(descriptor4,FMT="(8x,(I1))") var2_new 
                        indexer = resolution*(key/100000) + scaling1*var2_new+var1_new
                        key = 1
                        population = 1
                        counter1(indexer) = key
 
                else if (order == 2) then
                        read(descriptor3,FMT="(8x,(I1))") var1_new 
                        read(descriptor4,FMT="(8x,(I1))") var2_new 
                        indexer = resolution*(key/100000) + scaling1*var2_new+var1_new
                        key = 1
                        population = 1
                        counter2(indexer) = key

                else if (order == 3) then
                        read(descriptor3,FMT="(8x,(I1))") var1_new 
                        read(descriptor4,FMT="(8x,(I1))") var2_new 
                        indexer = resolution*(key/100000) + scaling1*var2_new+var1_new
                        key = 1
                        population = 1
                        counter3(indexer) = key

                else
                        exit
                end if
                
        end if

        !Add this state to its designated subcell in path3/subcell_new
        !subcell_new is of format r1_r2 where r = xx.yyyy where the number of
        !'y' elements is equal to .order.; the decimal is always included
        !the number of 'x' elements is variable

        open(72,file=trim(path3)//trim(subcell_new)//".dat",position="append",status=trim(descriptor0))
        write(72,FMT="(3(1x,F11.6))",advance="no") (vals(j),j=1,Nvar)
        write(72,FMT="(36(1x,F11.6))")(coords(j),j=1,6*Natoms)
        close(72)

        !By catching this scenario, we avoid the possibility of some frame
        !(r1,r2) being added to a 'deep' subcell when its parent is not
        !yet overcrowded.
        if (population < overcrowd) exit

        !If there is just the perfect number of frame, we need to divyUp
        !The format is the same for every order:
        !Array assignments are first-come-first-serve, with indexes assigned
        !by headerN+1 which is incremented whenever assigned.
        !Each frame may have .resolution. number of subcells so that many
        !indexes are put aside for them in counterN+1
        if (population == overcrowd) then 
                if (order == 0) then
                        key = key + 100000*header1
                        counter0(indexer) = key
                        call divyUp(vals(1),vals(2),0.0,order,header1,&
                                        counter1,250*resolution,overcrowd)
                        header1 = header1 + 1

!If this program fails midway then we want to know approximately when
if (modulo(header1,50) == 0) then
open(80,file=trim(path4)//trim(progressfile),position="append")
write(80,*) "         We will have ", header1, " overcrowded cells of order 0"
close(80)
end if

                else if (order == 1) then 
                        key = key + 100000*header2
                        counter1(indexer) = key
                        call divyUp(vals(1),vals(2),0.0,order,header2,&
                                        counter2,500*resolution,overcrowd)
                        header2 = header2 + 1

!If this program fails midway then we want to know approximately when
if (modulo(header2,100) == 0) then
open(80,file=trim(path4)//trim(progressfile),position="append")
write(80,*) "         We will have ", header2, " overcrowded subcells of order 1"
close(80)
end if

                else if (order == 2) then
                        key = key + 100000*header3
                        counter2(indexer) = key
                        call divyUp(vals(1),vals(2),0.0,order,header3,&
                                        counter3,250*resolution,overcrowd)
                        header3 = header3 + 1

!If this program fails midway then we want to know approximately when
if (modulo(header3,50) == 0) then
open(80,file=trim(path4)//trim(progressfile),position="append")
write(80,*) "         We will have ", header3, " overcrowded subcells of order 2"
close(80)
end if

                 else
 
        !There are two main ways our program will fail:
        !       1) a stupid error
        !       2) not allocating enough space for counterN (most likely)
        !If we start seeing this message pop-up or an array-out-of-bounds
        !segmentation fault, then we are probably either need more subcells
        !(more counterN) or a large value for overcrowding
open(80,file=trim(path4)//trim(progressfile),position="append")
write(80,*) "   Subcell: ", trim(subcell_new), " has gotten overcrowded; no more"//&
            &" resolution"
close(80)

                end if
                
                !There is a rare scenario when all frames in a cell will land
                !in a single subcell made from divyUp. In that case, we need to
                !loop again; however, I don't forsee this happening so I simply
                !exit and save maybe an extra millisecond.
                exit
        end if

        !Now that all that's done, we can go deeper
        order = order + 1

        !Essentially, this is like the truncation at the beginning
        !but for .order. number of digits to the right of the decimal
        !This part should be fine if we change scaling (from 10 to whatever)
        !but--again--the beginning part will not be okay
        var1 = anint((scaling1**order)*vals(1)-0.5)*(1.0/scaling1)**order
        var2 = anint((scaling2**order)*vals(2)-0.5)*(1.0/scaling2)**order

        !And then simply make the name of the new subcell
        write(descriptor1,FMT=FMT4) order
        write(descriptor3,FMT="(F9."//trim(adjustl(descriptor1))//")") var1
        write(descriptor4,FMT="(F9."//trim(adjustl(descriptor1))//")") var2
        subcell_new = trim(adjustl(descriptor3))//"_"//trim(adjustl(descriptor4))

end do



end subroutine addState



end module addCells3
