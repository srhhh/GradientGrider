
module addCells2
implicit none


contains


!The following grids the cell into smaller cells
!This ASSUMES the cell is NOT already divided up

recursive subroutine divyUp(var1,var2,var3,order,indexN,&
                                counterN,lengthN)
use f1_parameters
use f1_functions
implicit none
integer, intent(in) :: order,indexN,lengthN
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

var1_new = anint((scaling1**order)*var1-0.5)*(1.0/scaling1)**order
var2_new = anint((scaling2**order)*var2-0.5)*(1.0/scaling2)**order

write(descriptor1,FMT="(I1)") order
write(descriptor2,FMT="(F9."//trim(adjustl(descriptor1))//")") var1_new
write(descriptor3,FMT="(F9."//trim(adjustl(descriptor1))//")") var2_new
descriptor2 = adjustl(descriptor2)
descriptor3 = adjustl(descriptor3)
subcell = trim(descriptor2)//"_"//trim(descriptor3)

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
allocate(vals(overcrowd,Nvar))
allocate(coords(overcrowd,6*Natoms))
allocate(indexer(overcrowd,1))
do j=1, overcrowd
        read(72,FMT=FMT1,advance="no",iostat=k) (vals(j,i),i=1,Nvar)
        read(72,FMT=FMT2) (coords(j,i),i=1,6*Natoms)
        indexer(j,1) = j
end do
close(72)

!Because we are subdividing, the gaps are now shorter
gap1 = gap1/scaling1
gap2 = gap2/scaling2

!Sort the indexed states by the first variable (into columns); then grid it
!This sorts both vals and indexer; indexer can then be used to access coords
call qsort(vals,indexer,overcrowd,Nvar,1,overcrowd,1)
call grider(grid1,vals,gap1,var1_new,scaling1,overcrowd,Nvar,1,overcrowd,1,scaling1,1)

index1_1 = grid1(1)
do i = 1, scaling1
        !We always accept these grid elements (cells) in pairs of indexed states [p1,p2)
        !where indexed state p2 is not in the cell i
        !For the last element (cell) the last indexed state must be in it (Nstate)
        !For the first variable, we can think of these are "columns"
        if (i == scaling1) then
                index1_2 = overcrowd+1
        else
                index1_2 = grid1(i+1)
        end if

        !If the column is empty, don't even bother
        if (index1_2 == index1_1) then
                index1_1 = index1_2
                cycle
        end if

        !Sort the indexed states in this grid element (into cells); then grid it
        call qsort(vals,indexer,overcrowd,Nvar,index1_1,index1_2-1,2)
        call grider(grid2,vals,gap2,var2_new,scaling2,overcrowd,Nvar,&
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

                !Write these states into their respective cell in the grid
                !The array index needs to start at 1 so a subtraction is
                !involved
                write(descriptor4,FMT="(I1)") i-1
                write(descriptor5,FMT="(I1)") j-1
                descriptor4 = trim(descriptor2)//trim(adjustl(descriptor4))
                descriptor5 = trim(descriptor3)//trim(adjustl(descriptor5))
                subcell = trim(descriptor4)//"_"//trim(descriptor5)

                open(72,file=trim(path3)//trim(subcell)//".dat",status="new")
                do k = index2_1, index2_2-1
                        write(72,FMT=FMT1,advance="no") (vals(k,l),l=1,Nvar)
                        write(72,FMT=FMT2)(coords(indexer(k,1),l),l=1,6*Natoms)
                end do
                close(72)

                !And we also want to keep track of how many states are in this
                !new subcell.
                index_order = resolution*indexN + scaling1*(j-1) + i-1
                counterN(index_order) = index2_2 - index2_1

!               !If this single cell has all the states, call divyUp again
!               if (index2_2-index2_1 == overcrowd) then
!                       read(descriptor4,FMT="(F10."//trim(adjustl(descriptor6))//")") var1_new
!                       read(descriptor5,FMT="(F10."//trim(adjustl(descriptor6))//")") var2_new
!                       call divyUp(var1_new,var2_new,var3,&
!                                       order+1,counterN,indexerN,index_order)
!                       !This should technically not work since we're not going
!                       !deeper in terms of counterN, indexerN
!                       return
!               end if
               
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

subroutine addState(vals,coords)
use f1_parameters
implicit none
integer :: index_new,index_old,i,j,population
integer,save :: header1 = 1
integer,save :: header2 = 1
integer,save :: header3 = 1
integer :: overcrowd_p = overcrowd + 1
integer,dimension(ceiling(max_var1*max_var2)),save :: counter0 = 0
integer,dimension(ceiling(max_var1*max_var2)),save :: indexer1 = 0
integer,dimension(resolution*250),save :: counter1 = 0
integer,dimension(resolution*250),save :: indexer2 = 0
integer,dimension(resolution*500),save :: counter2 = 0
integer,dimension(resolution*500),save :: indexer3 = 0
integer,dimension(resolution*250),save :: counter3 = 0
logical :: flag1
real :: var1, var2, var1_old, var2_old, var3_old
integer :: var1_new,var2_new,var3_new,order
real, dimension(Nvar), intent(in) :: vals
real, dimension(6*Natoms), intent(in) :: coords
character(50) :: descriptor0, descriptor1, descriptor2
character(9) :: descriptor3, descriptor4
character(50) :: subcell_old, subcell_new

!At first, we only want to look at the first order spacing
!This hinges on the fact that the first-level spacing is 1 A
order = 0
write(descriptor1,FMT=FMT4) order
write(descriptor2,FMT=FMT4) order+1
write(descriptor3,FMT="(F9."//trim(descriptor1)//")") vals(1)-0.5
write(descriptor4,FMT="(F9."//trim(descriptor1)//")") vals(2)-0.5

subcell_new = trim(adjustl(descriptor3))//"_"//trim(adjustl(descriptor4))
population = 0

do

        !Check whether there is some subcell1 (var1,var2) in path3/
        !and reads the var1_var2.p file, increments it
        !to keep tabs on the number of states in var1_var2.dat
        inquire(file=trim(path3)//trim(subcell_new)//".dat",exist=flag1)
        if (flag1) then
                descriptor0 = "old"
                if (order == 0) then
                        read(descriptor3,FMT="(I8)") var1_new
                        read(descriptor4,FMT="(I8)") var2_new
                        index_new = ceiling(max_var1)*(var2_new-1) + var1_new
                        population = counter0(index_new) + 1
                        counter0(index_new) = population
                        index_old = index_new

                else if (order == 1) then
                        read(descriptor3,FMT="(8x,(I1))") var1_new 
                        read(descriptor4,FMT="(8x,(I1))") var2_new 
                        index_new = resolution*indexer1(index_old) + scaling1*(var2_new) + var1_new 
                        population = counter1(index_new) + 1
                        counter1(index_new) = population
                        index_old = index_new

                else if (order == 2) then
                        read(descriptor3,FMT="(8x,(I1))") var1_new 
                        read(descriptor4,FMT="(8x,(I1))") var2_new 
                        index_new = resolution*indexer2(index_old) + scaling1*(var2_new) + var1_new 
                        population = counter2(index_new) + 1
                        counter2(index_new) = population
                        index_old = index_new

                else if (order == 3) then
                        read(descriptor3,FMT="(8x,(I1))") var1_new 
                        read(descriptor4,FMT="(8x,(I1))") var2_new 
                        index_new = resolution*indexer3(index_old) + scaling1*(var2_new) + var1_new 
                        population = counter3(index_new) + 1
                        counter3(index_new) = population
                        index_old = index_new

                else
                        exit
                end if

        else
                descriptor0 = "new"

                if (order == 0) then
                        read(descriptor3,FMT="(I8)") var1_new
                        read(descriptor4,FMT="(I8)") var2_new
                        index_new = ceiling(max_var1)*(var2_new-1) + var1_new
                        population = 1
                        counter0(index_new) = population

                else if (population < overcrowd_p) then
                        !This means that the parent of this subcell has not yet
                        !been divyUp-ed. This should be the end of mosr calls
                        exit

                else if (order == 1) then
                        read(descriptor3,FMT="(8x,(I1))") var1_new 
                        read(descriptor4,FMT="(8x,(I1))") var2_new 
                        index_new = resolution*indexer1(index_old) + scaling1*(var2_new) + var1_new
                        population = 1
                        counter1(index_new) = population
 
                else if (order == 2) then
                        read(descriptor3,FMT="(8x,(I1))") var1_new 
                        read(descriptor4,FMT="(8x,(I1))") var2_new 
                        index_new = resolution*indexer2(index_old) + scaling1*(var2_new) + var1_new
                        population = 1
                        counter2(index_new) = population

                else if (order == 3) then
                        read(descriptor3,FMT="(8x,(I1))") var1_new 
                        read(descriptor4,FMT="(8x,(I1))") var2_new 
                        index_new = resolution*indexer3(index_old) + scaling1*(var2_new) + var1_new
                        population = 1
                        counter3(index_new) = population

                else
                        exit
                end if
                
        end if

        !Add this state to its designated subcell in path3/subcell1
        !subcell1 contains the parent subdirectories of subcell1 as well

        open(72,file=trim(path3)//trim(subcell_new)//".dat",position="append",status=trim(descriptor0))
        write(72,FMT="(3(1x,F11.6))",advance="no") (vals(j),j=1,Nvar)
        write(72,FMT="(36(1x,F11.6))")(coords(j),j=1,6*Natoms)
        close(72)

        !If there is just the perfect number of states, we need to divy it up
        if (population == overcrowd) then 
                if (order == 0) then
                        indexer1(index_old) = header1
                        call divyUp(vals(1),vals(2),0.0,order,header1,counter1,25*resolution)
                        header1 = header1 + 1

                else if (order == 1) then
                        indexer2(index_old) = header2
                        call divyUp(vals(1),vals(2),0.0,order,header2,counter2,500*resolution)
                        header2 = header2 + 1

                else if (order == 2) then
                        indexer3(index_old) = header3
                        call divyUp(vals(1),vals(2),0.0,order,header3,counter3,100*resolution)
                        header3 = header3 + 1

                 else
 
open(80,file=trim(path4)//trim(progressfile),position="append")
write(80,*) "   Subcell: ", trim(subcell_new), " has gotten overcrowded; no more"//&
            &" resolution"
close(80)

!               exit
                end if
                exit
        end if

        order = order + 1
        var1 = anint((scaling1**order)*vals(1)-0.5)*(1.0/scaling1)**order
        var2 = anint((scaling2**order)*vals(2)-0.5)*(1.0/scaling2)**order

        write(descriptor1,FMT=FMT4) order
        write(descriptor2,FMT=FMT4) order+1
        write(descriptor3,FMT="(F9."//trim(adjustl(descriptor1))//")") var1
        write(descriptor4,FMT="(F9."//trim(adjustl(descriptor1))//")") var2
        subcell_new = trim(adjustl(descriptor3))//"_"//trim(adjustl(descriptor4))

end do



end subroutine addState



end module addCells2
