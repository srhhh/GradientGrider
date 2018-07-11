module checkNewTrajectorywithGrid
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	SUBROTUINE addTrajectory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!INPUT:	real :: initial_translational_KE	"how fast the H is moving"
!		collision_distance		"how far away the H is"
!
!	real :: collision_skew			"the cross-sectional distance"
!
!	real ::	initial_bond_distance		"the length of the H2 bond initially"
!		initial_rotational_speed	"the rotational speed of the H2"
!		initial_rotation_angle		"the direction the H2 is spinning (relative)"
!
!	real :: initial_bond_angle1		"the angle the H2 bond makes with the x-axis"
!		initial_bond_angle2		"the angle the H2 bond makes with the y-z plane"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Simulates a H --> H2 collision with the above parameters
!	Uses physical constants and parameters as supplied by f2_physics_parameters.f90
!	Adds this data to a grid of files according to addState.f90
!	Makes use of a cross product which is supplied by f1_functions.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine checkTrajectory(initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
			   initial_bond_angle1,initial_bond_angle2,force_Neighbors,&
			   header1,header2,header3,counter0,counter1,counter2,counter3,path_to_grid)
        use PHYSICS
	use PARAMETERS
	use VARIABLES
	use addFrametoGrid
	use checkGrid
        implicit none

	!Coordinates, Velocities, and Variables
        real, dimension(3,Natoms) :: coords,velocities
	real, dimension(3,Natoms) :: gradient, approx_gradient
	real, dimension(Nvar) :: vals

	!Grid Parameters
	integer,intent(inout) :: header1,header2,header3
	integer,dimension(counter0_max),intent(out) :: counter0
	integer,dimension(counter1_max),intent(out) :: counter1
	integer,dimension(counter2_max),intent(out) :: counter2
	integer,dimension(counter3_max),intent(out) :: counter3
	logical,intent(in) :: force_Neighbors
        character(*),intent(in) :: path_to_grid

	!Collision Parameters
        real,intent(in) :: initial_bond_distance,initial_rotational_speed,initial_rotation_angle
	real,intent(in) :: initial_bond_angle1,initial_bond_angle2

	!Calculation Savers
        real, parameter :: PotentialConstant0 =&
                0.5*HOke_hydrogen*HOr0_hydrogen**2
        real, parameter :: PotentialConstant1 =&
                -HOke_hydrogen*HOr0_hydrogen
        real, parameter :: PotentialConstant2 =&
                0.5*HOke_hydrogen
        real, parameter :: AccelerationConstant0 =&
                HOke_hydrogen*(dt/mass_hydrogen)
        real, parameter :: AccelerationConstant1 =&
                -HOke_hydrogen*(dt/mass_hydrogen)*HOr0_hydrogen
        real, parameter :: AccelerationConstant2 =&
                2*MorseDe_hydrogen*Morsealpha_hydrogen*dt/mass_hydrogen

	!Various other variables
        real :: U, KE
	real :: min_rmsd,min_rmsd_prime
        integer :: TimeA,TimeB,step
	integer :: number_of_frames,order,neighbor_check

        !Initialize the scene
        call InitialSetup3(velocities,coords,&
                   	   initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
			   initial_bond_angle1,initial_bond_angle2)


        !Accelerate the velcocities for a half step (verlet)

	call getVar3(coords,Natoms,vals(1))
	call getVar4(coords,Natoms,vals(2))

        call Acceleration(vals,coords,gradient,&
             AccelerationConstant0,AccelerationConstant1,AccelerationConstant2)

        open(filechannel1,file=path_to_grid//trajectoryfile)
        write(filechannel1,'(I1)') 3
        write(filechannel1,*) ""
        write(filechannel1,'(A1,3F10.6)') 'H',&
              coords(1,1), coords(2,1), coords(3,1)
        write(filechannel1,'(A1,3F10.6)') 'H',&
              coords(1,2), coords(2,2), coords(3,2)
        write(filechannel1,'(A1,3F10.6)') 'H',&
              coords(1,3), coords(2,3), coords(3,3)
 	close(filechannel1)
  
	velocities = velocities + 0.5 * gradient

        !Keep track of the time
        TimeA = time()
open(filechannel2,file=path4//checkstatefile)
        do step = 1, Nsteps

		!Every 50 frames, print to an xyz file for visualization
                if (modulo(step,10) == 0) then
                        open(filechannel1,file=path4//trajectoryfile,position="append")
                        write(filechannel1,'(I1)') 3
                        write(filechannel1,*) ""
                        write(filechannel1,'(A1,3F10.6)') 'H',&
                                coords(1,1), coords(2,1), coords(3,1)
                        write(filechannel1,'(A1,3F10.6)') 'H',&
                                coords(1,2), coords(2,2), coords(3,2)
                        write(filechannel1,'(A1,3F10.6)') 'H',&
                                coords(1,3), coords(2,3), coords(3,3)
 			close(filechannel1)
                end if

                !Just to see progress, print something out every 500 steps
                if (modulo(step,500) == 1) then      

 			if ((vals(1)>max_var1) .or.&
 			    (vals(2)>max_var2)) then
 				exit
 			end if

!               	!Calculate the energies (done one-bby-one because there are few)
!               	U = MorsePotential(coords1,coords2)
!               	U = U + MorsePotential(coords1,coords3)
!               	U = U + HOPotential(coords2,coords3,&
!				PotentialConstant0,PotentialConstant1,PotentialConstant2)
!               	KE = KineticEnergy(velocity1)
!               	KE = KE + KineticEnergy(velocity2)
!               	KE = KE + KineticEnergy(velocity3)

!                       TimeB = time()
!                       print *, ""
!                       print *, "TIME STEP", step
!                       print *, "Time Passed: ", TimeB-TimeA
!                       print *, "KE: ", KE
!                       print *, "PE: ", U
!                       print *, "Total Energy: ", KE + U
!                       TimeA = TimeB
                endif

		!Upate the coordinates with the velocities
		coords = coords + dt * velocities

		!Get the variables
		call getVar3(coords,Natoms,vals(1))
		call getVar4(coords,Natoms,vals(2))
 
		!Check for similar frames
		min_rmsd = 100.0
		call checkState(coords,approx_gradient,min_rmsd_prime,.false.,&
				counter0,counter1,counter2,counter3,path_to_grid,&
				number_of_frames,order,neighbor_check)
		min_rmsd = min_rmsd_prime

		if (force_Neighbors) then
		call checkState(coords,approx_gradient,min_rmsd,.true.,&
				counter0,counter1,counter2,counter3,path_to_grid)
		end if

                !Update the gradient
                call Acceleration(vals,coords,gradient,&
			AccelerationConstant0,AccelerationConstant1,AccelerationConstant2)

write(filechannel2,FMT="(I8,1x,I1,1x,I1,1x)",advance="no") number_of_frames, order, neighbor_check
write(filechannel2,FMT="(I8,1x,E20.4,1x,E20.4,1x,F9.6,1x,F9.6)") step, min_rmsd, min_rmsd_prime, vals(1), vals(2)

		!Update the velocities
		velocities = velocities + gradient

        end do
close(filechannel2)

end subroutine checkTrajectory


subroutine Acceleration(vals,coords,gradient,&
		AccelerationConstant0,AccelerationConstant1,AccelerationConstant2)

        use PHYSICS
	use PARAMETERS
	use VARIABLES
        implicit none
	real, dimension(3,Natoms), intent(in) :: coords
	real, dimension(3,Natoms), intent(out) :: gradient
	real, dimension(3,Natoms) :: approx_gradient
	real, dimension(Nvar), intent(in) :: vals
        real, dimension(3) :: coords_distance,velocity_change,coords_atom1_bond
	real,intent(in) :: AccelerationConstant0,AccelerationConstant1,AccelerationConstant2
        real :: distance12,distance13,distance23,distance_squared
        real :: stretch_factor, distance_constant1, distance_constant2

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
                        stretch_factor - 1.0 ) * stretch_factor / distance12
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
                        stretch_factor - 1.0 ) * stretch_factor / distance13
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



subroutine InitialSetup3(velocities,coords, &
                   	 initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
			 initial_bond_angle1,initial_bond_angle2)
        use PHYSICS
	use PARAMETERS
	use FUNCTIONS
        implicit none
        real, dimension(3,Natoms), intent(out) :: velocities,coords
        real, dimension(3) :: bond_vector,rotation_vector,rotation0_vector,rotation90_vector
        real,intent(in) :: initial_bond_distance,initial_rotational_speed,initial_rotation_angle
	real,intent(in) :: initial_bond_angle1,initial_bond_angle2

	!Figure out which direction the bond is going
        bond_vector = (/ sin(initial_bond_angle1)*cos(initial_bond_angle2),&
		         sin(initial_bond_angle1)*sin(initial_bond_angle2),&
		         cos(initial_bond_angle1) /)

	!And two orthogonal vectors that are both perpendicular to the bond
	rotation0_vector =  (/ bond_vector(2), -bond_vector(1), 0.0 /) / &
                               sqrt(bond_vector(2)**2 + bond_vector(1)**2)
	call cross(bond_vector,rotation0_vector,rotation90_vector)

	!Calculate the coordinates of the atoms via the distance, skew, and bond vector
	bond_vector = initial_bond_distance * bond_vector
        coords(:,1) = (/ 0.0, 0.0, 0.0 /)
        coords(:,2) = (/ collision_distance, collision_skew, 0.0 /) + bond_vector/2
        coords(:,2) = (/ collision_distance, collision_skew, 0.0 /) - bond_vector/2

	!Calculate the direction the bond is spinning given the two
	!orthogonal vectors
	rotation_vector = ( sin(initial_rotation_angle) * rotation0_vector + &
                            cos(initial_rotation_angle) * rotation90_vector ) * &
                          initial_rotational_speed

	!With the translational, vibrational, and rotational vectors and speeds
	!Calculate the velocities of the atoms
        velocities(:,1) = (/ sqrt(2*initial_translational_KE/mass_hydrogen), 0.0, 0.0 /)
        velocities(:,2) = rotation_vector
        velocities(:,3) = -rotation_vector

open(progresschannel,file=path4//progressfile,position="append")
write(progresschannel,*) ""
write(progresschannel,*) "The parameters (in reduced units) are: "
write(progresschannel,*) ""
write(progresschannel,*) "Collision H - H2 Distance: ", collision_distance
write(progresschannel,*) "Initial H Velocity: ", velocities(1,1)
write(progresschannel,*) ""
write(progresschannel,*) "Initial H2 Bond distance: ", initial_bond_distance
write(progresschannel,*) "Initial H2 Rotation speed: ", initial_rotational_speed
write(progresschannel,*) "Initial H2 Rotation angle: ", initial_rotation_angle
close(progresschannel)

end subroutine InitialSetup3


real function MorsePotential(coords1,coords2)
        use PHYSICS
        implicit none
        integer :: i
        real, dimension(3), intent(in) :: coords1,coords2
        real, dimension(3) :: coords_distance
        real :: distance_squared,distance,stretch_factor

	!Calculate the distance
        coords_distance = coords2 - coords1
        distance_squared = coords_distance(1)**2 + coords_distance(2)**2 +&
                           coords_distance(3)**2

	!If its zero, something went wrong
        if (distance_squared == 0.0) then
                print *, "OVERLAP ERROR"
                print *, "Position of 1st atom:", coords1
                print *, "Position of 2nd atom: ", coords2
                call exit()
	end if

	distance = sqrt(distance_squared)
	stretch_factor = exp(Morsealpha_hydrogen*(distance-Morser0_hydrogen))
	MorsePotential = MorseDe_hydrogen*stretch_factor*(stretch_factor-2.0)

end function MorsePotential



real function HOPotential(coords1,coords2,&
		PotentialConstant0,PotentialConstant1,PotentialConstant2)
        use PHYSICS
        implicit none
        integer :: i
        real, dimension(3), intent(in) :: coords1,coords2
        real, dimension(3) :: coords_distance
	real,intent(in) :: PotentialConstant0,PotentialConstant1,PotentialConstant2
        real :: distance_squared,distance

	!Calculate the distance
        coords_distance = coords2 - coords1
        distance_squared = coords_distance(1)**2 + coords_distance(2)**2 +&
                           coords_distance(3)**2

	!If its zero, something went wrong
        if (distance_squared == 0.0) then
                print *, "OVERLAP ERROR"
                print *, "Position of 1st atom:", coords1
                print *, "Position of 2nd atom: ", coords2
                call exit()
	end if

        distance = sqrt(distance_squared)
        HOpotential = PotentialConstant2*distance_squared + &
                      PotentialConstant1*distance + &
                      PotentialConstant0

end function HOPotential



real function KineticEnergy(velocity)
        use PHYSICS
        implicit none
        real, dimension(3), intent(in) :: velocity
        real :: speed_squared

	!Compute the speed, then compute the kinetic energy
        speed_squared = velocity(1)**2 + velocity(2)**2 + velocity(3)**2
        KineticEnergy = 0.5*mass_hydrogen*speed_squared

end function KineticEnergy



end module checkNewTrajectorywithGrid
