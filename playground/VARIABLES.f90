module VARIABLES
use DOUBLE
implicit none

!This module defines the variables we're using to 'grid' the system
!and organize the data

contains

!Variable one is the distance between the midpoints of two bonds (ONLY FOR H2 - H2)
subroutine getVar1(coords,Natoms,var1)
implicit none
integer,intent(in) :: Natoms
real(dp), dimension(3,Natoms),intent(in) :: coords
real(dp), intent(out) :: var1
real(dp) :: length1, length2, length3, length4

!call getDistanceSquared(coords(:,1)+coords(:,2),coords(:,3)+coords(:,4),var1)
!var1 = 0.5 * sqrt(var1)

call getDistanceSquared(coords(:,1),coords(:,3),length1)
call getDistanceSquared(coords(:,2),coords(:,4),length2)
call getDistanceSquared(coords(:,1),coords(:,4),length3)
call getDistanceSquared(coords(:,2),coords(:,3),length4)
var1 = sqrt(min(length1,length2,length3,length4))

end subroutine getVar1

!Variable two is the cosine of the angle between the first H2 bond vector and
!						 the second H2 bond vector
!but shifted over by one, so that all values are positive
subroutine getVar2(coords,Natoms,var2)
implicit none
integer, intent(in) :: Natoms
real(dp), dimension(3,Natoms), intent(in) :: coords
!real(dp), dimension(3) :: bond_vector1, bond_vector2, bond_CM_vector
!real(dp) :: bond1_length, bond2_length, bond_CM_length
!real(dp), dimension(4) :: bond_lengths       !NOTE: This is not generic yet
!integer :: bond_pair
real(dp) :: length1, length2, length3, length4
real(dp), intent(out) :: var2

call getDistanceSquared(coords(:,1),coords(:,3),length1)
call getDistanceSquared(coords(:,2),coords(:,4),length2)
call getDistanceSquared(coords(:,1),coords(:,4),length3)
call getDistanceSquared(coords(:,2),coords(:,3),length4)
var2 = sqrt(max(length1,length2,length3,length4))

!call getDistanceSquared(coords(:,1),coords(:,3),bond_lengths(1))
!call getDistanceSquared(coords(:,2),coords(:,4),bond_lengths(2))
!call getDistanceSquared(coords(:,1),coords(:,4),bond_lengths(3))
!call getDistanceSquared(coords(:,2),coords(:,3),bond_lengths(4))
!bond_pair = minloc(bond_lengths,dim=1)
!
!if (bond_pair == 1) then
!	bond_vector1 = coords(:,2) - coords(:,1)
!	bond_vector2 = coords(:,4) - coords(:,3)
!	bond1_length = sum(bond_vector1**2)
!	bond2_length = sum(bond_vector2**2)
!
!else if (bond_pair == 2) then
!	bond_vector1 = coords(:,1) - coords(:,2)
!	bond_vector2 = coords(:,3) - coords(:,4)
!	bond1_length = sum(bond_vector1**2)
!	bond2_length = sum(bond_vector2**2)
!
!else if (bond_pair == 3) then
!	bond_vector1 = coords(:,2) - coords(:,1)
!	bond_vector2 = coords(:,3) - coords(:,4)
!	bond1_length = sum(bond_vector1**2)
!	bond2_length = sum(bond_vector2**2)
!
!else
!	bond_vector1 = coords(:,1) - coords(:,2)
!	bond_vector2 = coords(:,4) - coords(:,3)
!	bond1_length = sum(bond_vector1**2)
!	bond2_length = sum(bond_vector2**2)
!
!end if
!
!var2 = 1.0d0 + dot_product(bond_vector1,bond_vector2)/sqrt(bond1_length*bond2_length)

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

!Variables five and six are the distance between atoms 1 and 3 and 2 and 4, respectively
subroutine getVar5(coords,Natoms,var5)
implicit none
integer, intent(in) :: Natoms
real(dp), dimension(3,Natoms), intent(in) :: coords
real(dp), intent(out) :: var5

call getDistanceSquared(coords(:,1),coords(:,3),var5)
var5 = sqrt(var5)
end subroutine getVar5

subroutine getVar6(coords,Natoms,var6)
implicit none
integer, intent(in) :: Natoms
real(dp), dimension(3,Natoms), intent(in) :: coords
real(dp), intent(out) :: var6

call getDistanceSquared(coords(:,2),coords(:,4),var6)
var6 = sqrt(var6)
end subroutine getVar6



!Gets the squared distance between two sets of xyz xoordinates
subroutine getDistanceSquared(coord1,coord2,d)
implicit none
real(dp), dimension(3), intent(in) :: coord1, coord2
real(dp), intent(out) :: d

d = (coord1(1)-coord2(1))**2 + (coord1(2)-coord2(2))**2 + &
        (coord1(3)-coord2(3))**2

end subroutine getDistanceSquared

end module VARIABLES
