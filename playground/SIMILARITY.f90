!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       MODULE
!               SIMILARITY
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This module contains variables dictating the physics of the simulations
!               as well as what ensemble we are sampling initial conditions from and
!               some subroutines used to smooth this along
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILECHANNELS                    ACTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINES                     ARGUMENTS               KIND
!
!               InitialSetup3                   coords                  intent(out),real(dp),dim(3,Natoms)
!                                               velocities              intent(out),real(dp),dim(3,Natoms)
!
!               BondedForce                     coords1                 intent(in),real(dp),dim(3)
!                                               coords2                 intent(in),real(dp),dim(3)
!                                               gradient1               intent(out),real(dp),dim(3,Natoms)
!                                               gradient1               intent(out),real(dp),dim(3,Natoms)
!                                               r                       intent(in),real(dp),optional
!
!               NonBondedForce                  coords1                 intent(in),real(dp),dim(3)
!                                               coords2                 intent(in),real(dp),dim(3)
!                                               gradient1               intent(out),real(dp),dim(3,Natoms)
!                                               gradient1               intent(out),real(dp),dim(3,Natoms)
!                                               r                       intent(in),real(dp),optional
!
!               cross                           A                       intent(in),real(dp),dim(3)
!                                               B                       intent(in),real(dp),dim(3)
!                                               AcrossB                 intent(out),real(dp),dim(3)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       CALLS                           MODULE
!
!               cross                           PHYSICS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



module SIMILARITY
use DOUBLE
use PARAMETERS

integer,parameter :: NSIs = 3
integer,parameter :: Ntransform = 1

real(dp),dimension(3,Natoms) :: coords_target
real(dp),dimension(Natoms,Natoms) :: CM_target

real(dp),dimension(3,Natoms,NSIs) :: coords_transformed
real(dp),dimension(3,3,NSIs) :: U_transformed
real(dp),dimension(NSIs) :: SIs

real(dp),dimension(NSIs) :: default_SIs
real(dp),allocatable :: SIbuffer1(:,:)

real(dp),dimension(NSIs) :: min_SIs
real(dp),dimension(NSIs) :: largest_SIs
real(dp),dimension(NSIs) :: largest_weighted_SIs
real(dp),dimension(NSIs) :: largest_weighted_SIs2
real(dp),dimension(NSIs) :: interpolated_SIs


contains


subroutine getSIs(coords_in,coords_out,U_out,SIs_out)
    use PARAMETERS
    implicit none
    real(dp),dimension(3,Natoms),intent(in) :: coords_in
    real(dp),dimension(3,Natoms),intent(out) :: coords_out
    real(dp),dimension(3,3),intent(out) :: U_out
    real(dp),dimension(NSIs),intent(out) :: SIs_out

    call getRMSD(coords_in,1)
    call getCMD(charges,coords_in,2)
    call getCMD(masses,coords_in,3)

    coords_out = coords_transformed(:,:,Ntransform)
    U_out = U_transformed(:,:,Ntransform)
    SIs_out = SIs

    return

end subroutine getSIs

subroutine setTarget(coords_in)
    use PARAMETERS
    use ANALYSIS
    implicit none
    real(dp),dimension(3,Natoms),intent(in) :: coords_in
    integer :: i,j

    default_SIs(1) = 01.000100d0
    default_SIs(2) = 20.000100d0
    default_SIs(3) = 40.000100d0

    if (accept_worst) then
        do i = 1, NSIs
            SIbuffer1(i,:) = 0.0d0
        end do
    else
        do i = 1, NSIs
            SIbuffer1(i,:) = default_SIs(i)
        end do
    end if

    coords_target = coords_in
    call getCoulombMatrix(coords_in,CM_target)

    return

end subroutine setTarget


subroutine getCoulombMatrix(coords_in,CM)
    use PARAMETERS
    implicit none
    real(dp),dimension(3,Natoms),intent(in) :: coords_in
    real(dp),dimension(Natoms,Natoms) :: CM

    integer :: i, j

    CM = 0.0d0
    do i = 1, Natoms
    do j = 1, Natoms
        if (i < j) then
            cycle
        else if (i == j) then
            CM(i,j) = 1.0d0
        else
            CM(i,j) = 1.0d0 / &
                sqrt(sum((coords_in(:,i)-coords_in(:,j))**2))
        end if
    end do
    end do

    return

end subroutine getCoulombMatrix

subroutine getCMD(coefficients,coords_in,NSI)
    use PARAMETERS
    implicit none
    integer,intent(in) :: NSI
    real,dimension(Natoms),intent(in) :: coefficients
    real(dp),dimension(3,Natoms),intent(in) :: coords_in
    real(dp),dimension(Natoms,Natoms) :: CM
    real(dp) :: CMdiff
    integer :: i,j

    call getCoulombMatrix(coords_in,CM)

    CMdiff = 0.0d0
    do i = 1, Natoms
    do j = 1, Natoms
        if (i > j) then
            CMdiff = CMdiff + &
                     coefficients(i) * coefficients(j) * &
                     (CM(i,j) - CM_target(i,j))**2
        end if

    end do
    end do

    SIs(NSI) = sqrt(CMdiff/Natoms)
    coords_transformed(:,:,NSI) = coords_in
    U_transformed(:,:,NSI) = reshape((/ 1, 0, 0, &
                                        0, 1, 0, &
                                        0, 0, 1 /),(/ 3, 3 /))

    return

end subroutine getCMD

subroutine getRMSD(coords_in,NSI)
    use ls_rmsd_original
    use PARAMETERS
    implicit none
    integer,intent(in) :: NSI
    real(dp),dimension(3,Natoms),intent(in) :: coords_in
    real(dp),dimension(3) :: x_center,y_center
    integer :: i

    call rmsd_dp(Natoms,coords_in,coords_target,&
                 1,U_transformed(:,:,NSI),&
                 x_center,y_center,SIs(NSI))

    do i = 1, 3
        coords_transformed(i,:,NSI) =&
            coords_in(i,:) - x_center(i)
    end do

    coords_transformed(:,:,NSI) = matmul(&
        U_transformed(:,:,NSI),&
        coords_transformed(:,:,NSI))

    do i = 1, 3
        coords_transformed(i,:,NSI) =&
            coords_transformed(i,:,NSI) + y_center(i)
    end do

    return

end subroutine getRMSD

end module SIMILARITY
