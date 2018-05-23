
module checkCells
use ls_rmsd
use f1_parameters
use f1_functions
use f1_variables
implicit none


contains


!checkCells takes coordinates of a state, computes its related
!variables, then checks its cell in the database. This does
!not add it to the database yet; thus, it can be used in
!in a pseudo-neighborlist like Verlet

!It returns the path to the deepest subcell containing the state (file2)
!If order == 0, then file2 does not exist
!And it returns the number of states (population) in the deepest subcell

subroutine checkCells(coords,order,file2,population)
use f1_parameters
use f1_functions
use ls_rmsd
implicit none
integer :: state1,i
integer, intent(out) :: order, population
logical :: flag1
real :: gap1,gap2
integer :: var1_new,var2_new,var3_new,var1_old,var2_old,var3_old,order_new,line_num
real, dimension(3*Natoms), intent(in) :: coords
real, allocatable :: rmsds(:,:)
integer, allocatable :: indexer(:,:)
character(50) :: file2
character(50) :: line_data, descriptor1, descriptor2, file3,file4

!Used to temporarily store data
file4 = "tmp.txt"

!This will be the variable without truncation
call getVar1(coords,Natoms,var1_old)
call getVar2(coords,Natoms,var2_old)
call getVar3(coords,Natoms,var3_old)

!This is the spacing of the current subcell (gets shorter)
gap1 = spacing1
gap2 = spacing2

!This is the variable with truncation (to get the subcell)
var1_new = floor(var1/gap1)
var2_new = floor(var2/gap2)
var3_new = 1

!The order represents how 'deep' we are in directories
!Higher order means a shorter subcell length
!We return the population and order of the deepest subcell
order = 0
file2 = ""

!file3 stores the name of the current subcell of the state
!It is in the format var1_var2
write(descriptor1,FMT="(A50)") var1_new
write(descriptor2,FMT="(A50)") var2_new
file3 = trim(descriptor1)//"_"//trim(descriptor2)

do
        !Check whether var1_var2.dat is in path3/file3
        !If not, then we have found the deepest subcell
        !We return this
        inquire(file=trim(path3)//trim(file3)////file2//".dat",exist=flag1)
        if (flag1) then
                file2 = trim(file3)//".dat"
        else
                exit
        end if
 
        !Order represents how "deep" we are in subdirectories
        order = order+1

        !We check for a subcell in a subdirectory var1_var2
        call system("find "//trim(path3)//trim(file3)//"/ > "//trim(file4))
        inquire(file=trim(file4),size=i)
        
        !If there isn't a subdirectory for var1_var2, we have
        !also found the deepest subcell
        if (i == 0) then
                exit
        end if

        !In the other case where there IS a subdirectory then
        !we need to figure out what subcell the state is in
        var1_old = var1 - var1_new*gap1
        var2_old = var2 - var2_new*gap2
        gap1 = gap1/(scaling1)
        gap2 = gap2/(scaling2)
        var1_new = floor(var1_old/gap1)
        var2_new = floor(var2_old/gap2)

        !It will be named var1_var2
        write(descriptor1,FMT="(A50)") var1_new
        write(descriptor2,FMT="(A50)") var2_new
        file3 = file3//"/"//trim(descriptor1)//"_"//trim(descriptor2)
 
end do







!Let's see how large this subcell is
population = 0
if (order == 0) return


open(72,file=trim(path3)//trim(file2))
do
        read(72,FMT="(A50)",iostat=state1) line_data
        if (state1 /= 0) exit
        population = population + 1
end do
close(72)





!Now let's get the RMSD of each state
allocate(indexer(population,1),rmsds(population,1))
open(72,file=trim(path3)//trim(file2))
do i = 1, population
        read(72,FMT="",iostat=state1) line_data
        indexer(i) = i
        call rmsd(Natoms,)
end do
close(72)
call qsort(rmsds,indexer,population,1,1,population,1)






end subroutine checkCells



end module checkCells
