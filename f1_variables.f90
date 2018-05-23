module f1_variables
use f1_parameters
implicit none


!This module defines the variables we're using to 'grid' the system
!and organize the data


contains

!Variable one is the squared distance between atom 1 (fluorine) and 2 (carbon)
subroutine getVar1(coords,Natoms,var1)
implicit none
integer, intent(in) :: Natoms
real, dimension(3*Natoms), intent(in) :: coords
real, intent(out) :: var1

var1 =sqrt( (coords(4)-coords(1))**2 + (coords(5)-coords(2))**2 + (coords(6)-coords(3))**2)
end subroutine getVar1

!Variable two is the squared distance between atom 2 (carbon) and 6 (iodine)
subroutine getVar2(coords,Natoms,var2)
implicit none
integer, intent(in) :: Natoms
real, dimension(3*Natoms), intent(in) :: coords
real, intent(out) :: var2

var2 = sqrt((coords(16)-coords(4))**2 + (coords(17)-coords(5))**2 + (coords(18)-coords(6))**2)
end subroutine getVar2

!Variable three is the cosine of the angle between atoms 1 - 2 - 6
subroutine getVar3(coords,Natoms,var3)
implicit none
integer, intent(in) :: Natoms
real, dimension(3*Natoms), intent(in) :: coords
real, intent(out) :: var3
real :: x1,y1,z1,x2,y2,z2

x1 = coords(16)-coords(4)
y1 = coords(17)-coords(5)
z1 = coords(18)-coords(6)
x2 = coords(1)-coords(4)
y2 = coords(2)-coords(5)
z2 = coords(3)-coords(6)
var3 = (x1*x2+y1*y2+z1*z2)/sqrt((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))
end subroutine getVar3






end module f1_variables


