module f1_variables
use f1_parameters
implicit none

!This module defines the variables we're using to 'grid' the system
!and organize the data

contains

!Variable one is the squared distance between atom 1 (fluorine) and 2 (carbon)
!you can make this more general by change v1 to two args, which are the index of the atoms
subroutine getVar1(coords,Natoms,var1)
implicit none
integer, intent(in) :: Natoms
real, dimension(3*Natoms), intent(in) :: coords
real, intent(out) :: var1



!Later, consider using the squared value instead (less computationally heavy)
call distance_squared(coords(1:3),coords(4:6),var1)
var1 =sqrt(var1)

end subroutine getVar1

!Variable two is the squared distance between atom 2 (carbon) and 6 (iodine)
subroutine getVar2(coords,Natoms,var2)
implicit none
integer, intent(in) :: Natoms
real, dimension(3*Natoms), intent(in) :: coords
real, intent(out) :: var2



!Later, consider using the squared value instead (less computationally heavy)
call distance_squared(coords(4:6),coords(16:18),var2)
var2 = sqrt(var2)

end subroutine getVar2

!Variable three is the cosine of the angle between atoms 1 - 2 - 6
!currently not used
!again, make this more general by use 3 args
subroutine getVar3(coords,Natoms,var3)
implicit none
integer, intent(in) :: Natoms
real, dimension(3*Natoms), intent(in) :: coords
real, intent(out) :: var3
real :: x1,y1,z1,x2,y2,z2,r1,r2

x1 = coords(16)-coords(4)
y1 = coords(17)-coords(5)
z1 = coords(18)-coords(6)
x2 = coords(1)-coords(4)
y2 = coords(2)-coords(5)
z2 = coords(3)-coords(6)
var3 = (x1*x2+y1*y2+z1*z2)/sqrt((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))
end subroutine getVar3

!Gets the squared distance between two sets of xyz xoordinates
subroutine distance_squared(coord1,coord2,d)
implicit none
real, dimension(3), intent(in) :: coord1, coord2
real, intent(out) :: d

d = (coord1(1)-coord2(1))**2 + (coord1(2)-coord2(2))**2 + &
        (coord1(3)-coord2(3))**2

end subroutine distance_squared

end module f1_variables
