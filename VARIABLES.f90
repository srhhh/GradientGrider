module VARIABLES
use PHYSICS
implicit none

!This module defines the variables we're using to 'grid' the system
!and organize the data

contains

!Variable one is the squared distance between the hydrogen and the midpoint of the bond
subroutine getVar1(coordsH,coords_bond,var1)
implicit none
real(dp), dimension(3),intent(in) :: coordsH,coords_bond
real(dp), intent(out) :: var1

!Later, consider using the squared value instead (less computationally heavy)
call getDistanceSquared(coordsH,coords_bond,var1)
var1 =sqrt(var1)

end subroutine getVar1

!Variable two is the cosine of the angle between the H - bond_midpoint vector and
!						 the bond_mindpoint - bond_endpoint vector
!but shifted over by one, so that all values are positive
subroutine getVar2(coords,Natoms,coords_bond,H_bond_distance,var2)
implicit none
integer, intent(in) :: Natoms
real(dp), intent(in) :: H_bond_distance
real(dp), dimension(3,Natoms), intent(in) :: coords
real(dp), dimension(3),intent(in) :: coords_bond
real(dp), intent(out) :: var2

call getDistanceSquared(coords_bond,coords(:,3),var2)
var2 = dot_product(coords(:,3)-coords_bond, coords_bond) / &
                 (H_bond_distance*sqrt(var2)) + 1


end subroutine getVar2

!Variables three and four are the distance between atoms 1 and 2 and 1 and 3, respectively
subroutine getVar3(coords,Natoms,var3)
implicit none
integer, intent(in) :: Natoms
real(dp), dimension(3,Natoms), intent(in) :: coords
real(dp), intent(out) :: var3

call getDistanceSquared(coords(:,1),coords(:,2),var3)
var3 = sqrt(var3)
end subroutine getVar3

subroutine getVar4(coords,Natoms,var4)
implicit none
integer, intent(in) :: Natoms
real(dp), dimension(3,Natoms), intent(in) :: coords
real(dp), intent(out) :: var4

call getDistanceSquared(coords(:,1),coords(:,3),var4)
var4 = sqrt(var4)
end subroutine getVar4




!Gets the squared distance between two sets of xyz xoordinates
subroutine getDistanceSquared(coord1,coord2,d)
implicit none
real(dp), dimension(3), intent(in) :: coord1, coord2
real(dp), intent(out) :: d

d = (coord1(1)-coord2(1))**2 + (coord1(2)-coord2(2))**2 + &
        (coord1(3)-coord2(3))**2

end subroutine getDistanceSquared

end module VARIABLES
