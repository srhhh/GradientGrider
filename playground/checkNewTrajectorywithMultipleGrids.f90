module checkNewTrajectorywithMultipleGrids
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	SUBROTUINE checkMultipleTrajectories
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!INPUT: real ::	initial_bond_distance		"the length of the H2 bond initially"
!		initial_rotational_speed	"the rotational speed of the H2"
!		initial_rotation_angle		"the direction the H2 is spinning (relative)"
!
!	real :: initial_bond_angle1		"the angle the H2 bond makes with the x-axis"
!		initial_bond_angle2		"the angle the H2 bond makes with the y-z plane"
!
!	logical            :: force_Neighbors		"the choice to check adjacent cells"
!	integer            :: Ngrid_total			"the number of grids to check"
!	int,dim(Ngrid_total) :: filechannels		"the files that we write RMSD calls to"
!	character(*)       :: path_to_directory		"the path to the directory with the grids"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Simulates a H --> H2 collision with the above parameters
!	Uses physical constants and parameters as supplied by f2_physics_parameters.f90
!	Checks a folder full of grids of files formatted by f2_parameters.f90
!	Checks each individual frame according to subroutines in checkCells5.f90
!	Makes use of a cross product which is supplied by f1_functions.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine checkMultipleTrajectories(initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
				     initial_bond_angle1,initial_bond_angle2,force_Neighbors,threshold_RMSD,&
				     reject_flag,Ngrid_total,filechannels,path_to_directory,velocityH,velocityH2)
        use PHYSICS
	use PHYSICSFUNCTIONS
	use PARAMETERS
	use VARIABLES
	use checkMultipleGrids
        implicit none

	!Coordinates, Velocities, and Variables
	real(dp),dimension(3,Natoms) :: coords,velocities
	real(dp),dimension(3,Natoms) :: gradient,approx_gradient
	real(dp),dimension(Nvar) :: vals
	real(dp),dimension(3),intent(out) :: velocityH, velocityH2
	real(dp),dimension(3) :: velocity_rotation,velocity_vibration,velocity_translation

	!Grid Parameters
	logical,intent(in) :: force_Neighbors
	real(dp),intent(in) :: threshold_RMSD
	logical,intent(in) :: reject_flag
        character(*),intent(in) :: path_to_directory
	integer,intent(in) :: Ngrid_total
	integer,dimension(Ngrid_total),intent(in) :: filechannels
	character(4) :: Ngrid_text

	!Collision Parameters
        real(dp),intent(in) :: initial_bond_distance,initial_rotational_speed,initial_rotation_angle
	real(dp),intent(in) :: initial_bond_angle1,initial_bond_angle2

	!Calculation Savers
        real(dp), parameter :: PotentialConstant0 =&
                0.5d0*HOke_hydrogen*HOr0_hydrogen**2
        real(dp), parameter :: PotentialConstant1 =&
                -HOke_hydrogen*HOr0_hydrogen
        real(dp), parameter :: PotentialConstant2 =&
                0.5d0*HOke_hydrogen
        real(dp), parameter :: AccelerationConstant0 =&
                HOke_hydrogen*(dt/mass_hydrogen)
        real(dp), parameter :: AccelerationConstant1 =&
                -HOke_hydrogen*(dt/mass_hydrogen)*HOr0_hydrogen
        real(dp), parameter :: AccelerationConstant2 =&
                2*MorseDe_hydrogen*Morsealpha_hydrogen*dt/mass_hydrogen

	!Various other variables
        real(dp) :: U, KE
	real(dp) :: min_rmsd,min_rmsd_prime
        integer :: TimeA,TimeB,step
	integer :: number_of_frames,order,neighbor_check,Ngrid

        !Initialize the scene
	call InitialSetup3(coords,velocities,&
                   	   initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
			   initial_bond_angle1,initial_bond_angle2)

        velocityH = velocities(:,1) - (velocities(:,1)+velocities(:,2)+velocities(:,3))/3

        open(filechannel1,file=path_to_directory//trajectoryfile)
        write(filechannel1,'(I1)') 3
        write(filechannel1,*) ""
        write(filechannel1,'(A1,3F10.6)') 'H',&
              coords(1,1), coords(2,1), coords(3,1)
        write(filechannel1,'(A1,3F10.6)') 'H',&
              coords(1,2), coords(2,2), coords(3,2)
        write(filechannel1,'(A1,3F10.6)') 'H',&
              coords(1,3), coords(2,3), coords(3,3)
        close(filechannel1)


	!Calculate the variables for use in acceleration
	call getVar3(coords,Natoms,vals(1))
	call getVar4(coords,Natoms,vals(2))

        !Accelerate the velcocities for a half step (verlet)
        call Acceleration(vals,coords,gradient,&
             AccelerationConstant0,AccelerationConstant1,AccelerationConstant2)

	velocities = velocities + 0.5d0 * gradient

        !Keep track of the time
        TimeA = time()
        do step = 1, Nsteps

if (modulo(step,50) == 1) then
                !Just to see progress, print something out every 500 steps
		!If we are too far out, stop the simulation

        open(filechannel1,file=path_to_directory//trajectoryfile,position="append")
        write(filechannel1,'(I1)') 3
        write(filechannel1,*) ""
        write(filechannel1,'(A1,3F10.6)') 'H',&
              coords(1,1), coords(2,1), coords(3,1)
        write(filechannel1,'(A1,3F10.6)') 'H',&
              coords(1,2), coords(2,2), coords(3,2)
        write(filechannel1,'(A1,3F10.6)') 'H',&
              coords(1,3), coords(2,3), coords(3,3)
        close(filechannel1)
                if (modulo(step,500) == 1) then      


 			if ((vals(1)>max_var1) .or.&
 			    (vals(2)>max_var2)) then
 				exit
 			end if
                endif
end if

		!Upate the coordinates with the velocities
		coords = coords + dt * velocities

		!Get the variables
		call getVar3(coords,Natoms,vals(1))
		call getVar4(coords,Natoms,vals(2))
 
		!Check for similar frames
		call checkState(vals,coords,approx_gradient,min_rmsd,&
				force_Neighbors,path_to_directory,Ngrid_total,filechannels,&
				number_of_frames,order,neighbor_check)

                !Update the gradient
		if ((min_rmsd .le. threshold_RMSD).or.(reject_flag)) then
                	call Acceleration(vals,coords,gradient,&
				AccelerationConstant0,AccelerationConstant1,AccelerationConstant2)
		else
			gradient = approx_gradient
		end if

		!Update the velocities
		velocities = velocities + gradient

        end do

        velocityH2 = velocities(:,1) - (velocities(:,1)+velocities(:,2)+velocities(:,3))/3

end subroutine checkMultipleTrajectories


subroutine Acceleration(vals,coords,gradient,&
		AccelerationConstant0,AccelerationConstant1,AccelerationConstant2)

        use PHYSICS
	use PARAMETERS
	use VARIABLES
        implicit none
	real(dp), dimension(3,Natoms), intent(in) :: coords
	real(dp), dimension(3,Natoms),intent(out) :: gradient
	real(dp), dimension(Nvar), intent(in) :: vals
        real(dp), dimension(3) :: coords_distance,velocity_change,coords_atom1_bond
	real(dp),intent(in) :: AccelerationConstant0,AccelerationConstant1,AccelerationConstant2
        real(dp) :: distance12,distance13,distance23,distance_squared
        real(dp) :: stretch_factor, distance_constant1, distance_constant2

	!Use the precalculated variables to our advantage
	distance12 = vals(1)
	distance13 = vals(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Morse Potential Derivative of Atoms 1 and 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!Save some calculation with this factor
	stretch_factor = exp(Morsealpha_hydrogen*(distance12-Morser0_hydrogen))

	!Calculate the velocity change
        coords_distance = coords(:,2) - coords(:,1)
        velocity_change = coords_distance * AccelerationConstant2 * ( &
                        stretch_factor - 1.0d0) * stretch_factor / distance12
        gradient(:,1) = velocity_change
        gradient(:,2) = -velocity_change

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Morse Potential Derivative of Atoms 1 and 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!Save some calculation with this factor
	stretch_factor = exp(Morsealpha_hydrogen*(distance13-Morser0_hydrogen))

	!Calculate the velocity change
        coords_distance = coords(:,3) - coords(:,1)
        velocity_change = coords_distance * AccelerationConstant2 * ( &
                        stretch_factor - 1.0d0) * stretch_factor / distance13
        gradient(:,1) = gradient(:,1) + velocity_change
        gradient(:,3) = -velocity_change

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Harmonic Oscillator of Atoms 2 and 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!Calculate the distance
        coords_distance = coords(:,3) - coords(:,2)
        distance_squared = coords_distance(1)**2 + coords_distance(2)**2 +&
                           coords_distance(3)**2
        distance23 = sqrt(distance_squared)

	!Calculate the velocity change
        velocity_change = coords_distance * ( &
                        AccelerationConstant0 + &
                        AccelerationConstant1/distance23 )
        gradient(:,2) = gradient(:,2) + velocity_change
        gradient(:,3) = gradient(:,3) - velocity_change

	return

end subroutine Acceleration



subroutine InitialSetup3(coords,velocities, &
                   	 initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
			 initial_bond_angle1,initial_bond_angle2)
        use PHYSICS
	use PARAMETERS
	use FUNCTIONS
        implicit none
	real(dp),dimension(3,Natoms),intent(out) :: coords,velocities
        real(dp), dimension(3) :: bond_vector,rotation_vector,rotation0_vector,rotation90_vector
        real(dp),intent(in) :: initial_bond_distance,initial_rotational_speed,initial_rotation_angle
	real(dp),intent(in) :: initial_bond_angle1,initial_bond_angle2

	!Figure out which direction the bond is going
        bond_vector = (/ sin(initial_bond_angle1)*cos(initial_bond_angle2),&
		         sin(initial_bond_angle1)*sin(initial_bond_angle2),&
		         cos(initial_bond_angle1) /)

	!And two orthogonal vectors that are both perpendicular to the bond
	rotation0_vector =  (/ bond_vector(2), -bond_vector(1), 0.0d0 /) / &
                               sqrt(bond_vector(2)**2 + bond_vector(1)**2)
	call cross(bond_vector,rotation0_vector,rotation90_vector)

	!Calculate the coordinates of the atoms via the distance, skew, and bond vector
	bond_vector = initial_bond_distance * bond_vector
        coords(:,1) = (/ 0.0d0, 0.0d0, 0.0d0 /)
        coords(:,2) = (/ collision_distance, collision_skew, 0.0d0 /) + bond_vector/2
        coords(:,3) = (/ collision_distance, collision_skew, 0.0d0 /) - bond_vector/2

	!Calculate the direction the bond is spinning given the two
	!orthogonal vectors
	rotation_vector = ( sin(initial_rotation_angle) * rotation0_vector + &
                            cos(initial_rotation_angle) * rotation90_vector ) * &
                          initial_rotational_speed

	!With the translational, vibrational, and rotational vectors and speeds
	!Calculate the velocities of the atoms
        velocities(:,1) = (/ sqrt(2*initial_translational_KE/mass_hydrogen), 0.0d0, 0.0d0 /)
        velocities(:,2) = rotation_vector
        velocities(:,3) = -rotation_vector

end subroutine InitialSetup3


real(dp) function MorsePotential(coords1,coords2)
        use PHYSICS
        implicit none
        integer :: i
        real(dp), dimension(3), intent(in) :: coords1,coords2
        real(dp), dimension(3) :: coords_distance
        real(dp) :: distance_squared,distance,stretch_factor

	!Calculate the distance
        coords_distance = coords2 - coords1
        distance_squared = coords_distance(1)**2 + coords_distance(2)**2 +&
                           coords_distance(3)**2

	distance = sqrt(distance_squared)
	stretch_factor = exp(Morsealpha_hydrogen*(distance-Morser0_hydrogen))
	MorsePotential = MorseDe_hydrogen*stretch_factor*(stretch_factor-2.0d0)

end function MorsePotential



real(dp) function HOPotential(coords1,coords2,&
		PotentialConstant0,PotentialConstant1,PotentialConstant2)
        use PHYSICS
        implicit none
        integer :: i
        real(dp), dimension(3), intent(in) :: coords1,coords2
        real(dp), dimension(3) :: coords_distance
	real(dp),intent(in) :: PotentialConstant0,PotentialConstant1,PotentialConstant2
        real(dp) :: distance_squared,distance

	!Calculate the distance
        coords_distance = coords2 - coords1
        distance_squared = coords_distance(1)**2 + coords_distance(2)**2 +&
                           coords_distance(3)**2

        distance = sqrt(distance_squared)
        HOpotential = PotentialConstant2*distance_squared + &
                      PotentialConstant1*distance + &
                      PotentialConstant0

end function HOPotential



real(dp) function KineticEnergy(velocity)
        use PHYSICS
        implicit none
        real(dp), dimension(3), intent(in) :: velocity
        real(dp) :: speed_squared

	!Compute the speed, then compute the kinetic energy
        speed_squared = velocity(1)**2 + velocity(2)**2 + velocity(3)**2
        KineticEnergy = 0.5d0*mass_hydrogen*speed_squared

end function KineticEnergy



end module checkNewTrajectorywithMultipleGrids
