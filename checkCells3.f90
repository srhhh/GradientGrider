
module checkCells3
implicit none


contains


subroutine checkState(coords,closestCoords,min_rmsd)
use ls_rmsd
use addCells3
use f1_functions
use f1_variables
use f1_parameters
implicit none
integer :: state1, i,j,k,frames_max,order,Nstates
integer, dimension(1) :: j_tuple
logical :: flag1
real :: var1,var2,var3, var1_old, var2_old
real, dimension(3*Natoms), intent(in) :: coords
real, dimension(6*Natoms), intent(out) :: closestCoords
double precision, dimension(3,Natoms) :: rmsd_coords1,rmsd_coords2
double precision, dimension(3) :: x_center,y_center
double precision, intent(out) :: min_rmsd
double precision, allocatable :: neighbor_rmsds(:,:)
real, allocatable :: neighbor_coords(:,:)
double precision, allocatable :: U(:,:), g(:,:)
character(50) :: subcell_old, descriptor1, descriptor2, subcell_new
character(9) ::  descriptor3, descriptor4

!First, get the values of each respective variable 
call getVar1(coords,Natoms,var1)
call getVar2(coords,Natoms,var2)
call getVar3(coords,Natoms,var3)

!We also need to reshape the coordinates for the rmsd module
rmsd_coords1 = reshape(coords,(/3, Natoms/))

!The first subcell will be the digits to the left of the decimal
!Will not work well if the scaling/spacing changes
var1_old = anint(var1 - 0.5)
var2_old = anint(var2 - 0.5)

!The order keeps track of how deep we go into subcells
!It is used in the subroutine call for divyUp
order = 0
do

        !First, get the name of the subcelll
        write(descriptor1,FMT=FMT4) order
        write(descriptor2,FMT=FMT4) order+1
        write(descriptor3,FMT="(F9."//trim(descriptor1)//")") var1_old
        write(descriptor4,FMT="(F9."//trim(descriptor1)//")") var2_old
        subcell_new = trim(adjustl(descriptor3))//"_"//trim(adjustl(descriptor4))

        !Inquire whether it exists
        inquire(file=trim(path3)//trim(subcell_new)//".dat",exist=flag1)

        !If it does, we may be able to go deeper
        if (flag1) then
                subcell_old = subcell_new

        !If it doesn't, then we must use its parent subcell
        !AKA the subcell of the previous cycle
        else
                !Depending on the order (how deep we are) we can expect
                !different caps on how many states the subcell has
                if (order==0) then
                        print *, "no top-level cell with this state"
                        min_rmsd = 100
                        return
                else if (order==1) then
                        frames_max = overcrowd0
                else if (order==2) then
                        frames_max = overcrowd1
                else if (order==3) then
                        frames_max = overcrowd2
                else
                        frames_max = 1000  !Not overcrowd3 because currently we
                                           !do not subdivide past N = 3
                end if

open(70,file=trim(path4)//trim(progressfile),position="append")
write(70,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
write(70,*) ""
write(70,*) "For the closest frame, looking in"
write(70,*) trim(subcell_old)
write(70,*) ""
close(70)

                !Set the array to zero (because we won't traverse all of it)
                allocate(neighbor_rmsds(frames_max,1),&
                        neighbor_coords(frames_max,Nvar+6*Natoms))
                do i = 1, frames_max
                        neighbor_rmsds(i,1) = 100.0
                end do

                !Read off the states in this subcell
                open(72,file=trim(path3)//trim(subcell_old)//".dat")
                Nstates = 1
                do
                        read(72,FMT=FMT1,advance="no",iostat=state1) &
                                (neighbor_coords(Nstates,k),k=1,Nvar)
                        if (state1 /= 0) then
                                Nstates = Nstates - 1
                                exit
                        end if
                        read(72,FMT=FMT2) &
                                (neighbor_coords(Nstates,k),k=Nvar+1,Nvar+6*Natoms)

                        !Need to make the coordinates readable for the rmsd
                        rmsd_coords2 = reshape&
                                (neighbor_Coords(Nstates,Nvar+1:Nvar+3*Natoms),&
                                (/3, Natoms/))
                        call rmsd(Natoms,rmsd_coords1,rmsd_coords2,0,U,&
                                x_center,y_center,neighbor_rmsds(Nstates,1),.false.,g)
                        Nstates = Nstates + 1
                end do
                close(72)

                !Now, we want the state that is closest in terms of rmsd
                j_tuple = minloc(neighbor_rmsds(:,1))
                min_rmsd = neighbor_rmsds(j_tuple(1),1)
                closestCoords = neighbor_coords(j_tuple(1),Nvar+1:Nvar+6*Natoms)

open(70,file=trim(path4)//trim(progressfile),position="append")
write(70,*) "RMSD:"
write(70,*) ""
do i = 1, Nstates
write(70,*) neighbor_rmsds(i,1)
end do
write(70,*) ""
write(70,*) "closest one is index: ", j_tuple(1)
write(70,*) ""
write(70,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
write(70,*) ""
close(70)

                !That is all that is required by this subroutine
                exit
        end if

        !We can now check if there is a deeper subcell
        order = order+1
        var1_old = anint((scaling1**order)*var1 - 0.5)*(1.0/scaling1)**order
        var2_old = anint((scaling2**order)*var2 - 0.5)*(1.0/scaling2)**order
 
end do



end subroutine checkState



end module checkCells3
