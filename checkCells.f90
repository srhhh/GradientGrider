
module checkCells
implicit none


contains


subroutine checkState(coords,closestCoords,min_rmsd,flag2_p,gradient_p)
use ls_rmsd
use addCells
use f1_functions
use f1_variables
use f1_parameters
implicit none
integer :: state1,i,j,k
integer, dimension(1) :: j_tuple
logical :: flag1,flag2
logical, optional, intent(in) :: flag2_p
real, dimension(3*Natoms), optional, intent(in) :: gradient_p
real, dimension(3*Natoms) :: gradient
real :: gap1,gap2, var1,var2,var3, var1_old, var2_old, var3_old
integer :: var1_new,var2_new,var3_new,order,line_num
real, dimension(Nvar) :: vals
real, dimension(3*Natoms), intent(in) :: coords
real, dimension(6*Natoms), intent(out) :: closestCoords
double precision, dimension(3,Natoms) :: rmsd_coords1,rmsd_coords2
double precision, dimension(3) :: x_center,y_center
double precision, intent(out) :: min_rmsd
double precision, allocatable :: neighbor_rmsds(:,:)
real, allocatable :: neighbor_coords(:,:)
double precision, allocatable :: U(:,:), g(:,:)
character(50) :: line_data, descriptor1, descriptor2, subcell0, subcell1

!We set this false if we want to add this as a new state (then gradient is
!required)
if (present(flag2_p).and.present(gradient_p)) then
        flag2 = flag2_p
        gradient = gradient_p
else
        flag2 = .false.
        gradient = gradient_p
end if

!First, get the values of each respective variable 
call getVar1(coords,Natoms,var1)
call getVar2(coords,Natoms,var2)
call getVar3(coords,Natoms,var3)

!We also need to reshape the coordinates for the rmsd module
rmsd_coords1 = reshape(coords,(/3, Natoms/))

!var_old changes, essentially only storing the trailing digits
var1_old = var1
var2_old = var2

!As we go deeper in the subcells, the length of the subcell gets smaller
!And var_new represents the index of the subcell in the grid
gap1 = spacing1
gap2 = spacing2
var1_new = floor(var1/gap1)
var2_new = floor(var2/gap2)
var3_new = 1

!The order keeps track of how deep we go into subcells
!It is used in the subroutine call for divyUp
order = 0

!The file is named by the indexes of the rounded variables
write(descriptor1,FMT=FMT4) var1_new
write(descriptor2,FMT=FMT4) var2_new
subcell0 = ""
subcell1 = trim(adjustl(descriptor1))//"_"//trim(adjustl(descriptor2))

do
        !Check whether there is some subcell1 (var1,var2) in path3/subcell0
        !and reads the var1_var2.p file, increments it
        !to keep tabs on the number of states in var1_var2.dat
        inquire(file=trim(path3)//trim(subcell1)//".dat",exist=flag1)
        if (flag1) then
                open(72,file=trim(path3)//trim(subcell1)//".p",status="old")
                read(72,FMT="(I9)") i
                if (flag2) then
                        rewind(72)
                        i = i + 1
                        write(72,FMT="(I9)") i
                        descriptor1 = "old"
                end if
                close(72)
        !If not, we specify two new files: var1_var2.dat, var1_var2.p
        !And we keep tabs on its number of states with var1_var2.p
        else
                if (flag2) then
                        open(72,file=trim(path3)//trim(subcell1)//".p",status="new")
                        i = 1
                        write(72,FMT=FMT5) i
                        close(72)
                        descriptor1 = "new"
                end if
        end if

        !Add this state to its designated subcell in path3/subcell1
        !subcell1 contains the parent subdirectories of subcell1 as well
        if (flag2) then
                open(72,file=trim(path3)//trim(subcell1)//".dat",position="append",status=trim(descriptor1))
                write(72,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
                write(72,FMT=FMT3,advance="no")(coords(j),j=1,3*Natoms)
                write(72,FMT=FMT3)(gradient(j),j=1,3*Natoms)
                close(72)
        end if

        !If there are not many states, this is the 'deepest' subcell;
        !Now we can get the rmsd of the neighbors
        if (i < overcrowd) then
                allocate(neighbor_rmsds(i,1),neighbor_coords(i,Nvar+6*Natoms))

                !Read off the states in this subcell
                open(72,file=trim(path3)//trim(subcell1)//".dat")
                do j = 1, i
                        read(72,FMT=FMT1,advance="no") &
                                (neighbor_coords(j,k),k=1,Nvar)
                        read(72,FMT=FMT2) &
                                (neighbor_coords(j,k),k=Nvar+1,Nvar+6*Natoms)

                        !Need to make the coordinates readable for the rmsd
                        rmsd_coords2 = reshape&
                                (neighbor_Coords(j,Nvar+1:Nvar+3*Natoms),&
                                (/3, Natoms/))
                        call rmsd(Natoms,rmsd_coords1,rmsd_coords2,0,U,&
                                x_center,y_center,neighbor_rmsds(j,1),.false.,g)
                end do
                close(72)

                !Now, we want the state that is closest in terms of rmsd
                j_tuple = minloc(neighbor_rmsds(:,1))
                min_rmsd = neighbor_rmsds(j_tuple(1),1)
                closestCoords = neighbor_coords(j_tuple(1),Nvar+1:Nvar+6*Natoms)

                !That is all that is required by this subroutine
                exit
        end if

        !We want a way to keep track of which 'digit' we are on
        var1_old = var1_old - var1_new*gap1
        var2_old = var2_old - var2_new*gap2

        !If there is just the perfect number of states, we need to divy it up
        if ((i == overcrowd).and.(.not.(flag2))) then 
                call divyUp(var1-var1_old,var2-var2_old,0.0,&
                                order,subcell0,i)
                exit
        end if

        !If there are many states, then we can go deeper (into a subdirectory)
        !So now we need to figure out what subcell the state is in

        !We lop off the heading digit used in the previous subcell
        !because we only make use of the trailing one
        !Consequently, we resolve and the subcell gap size shrinks
        gap1 = gap1/(scaling1)
        gap2 = gap2/(scaling2)
        var1_new = floor(var1_old/gap1)
        var2_new = floor(var2_old/gap2)

        !This new subcell will be named var1_var2
        !And now the path to get there is lengthened
        write(descriptor1,FMT=FMT4) var1_new
        write(descriptor2,FMT=FMT4) var2_new
        subcell0 = trim(subcell1)//"/"
        subcell1 = trim(subcell0)//trim(adjustl(descriptor1))//"_"//&
                        trim(adjustl(descriptor2))
 
        !Order represents how "deep" we are in subdirectories
        order = order+1
 
end do



end subroutine checkState



end module checkCells
