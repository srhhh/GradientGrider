module makeTrajectory2
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

subroutine addTrajectory(initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
			 initial_bond_angle1,initial_bond_angle2,&
			 header1,header2,header3,counter0,counter1,counter2,counter3)
        use f2_physics_parameters
	use f2_parameters
	use f2_variables
	use addCells5
	use checkCells4
        implicit none

	!Coordinates, Velocities, and Variables
        real, dimension(Ncoords*2) :: coords_gradient,closestCoords
        real, dimension(3) :: velocity1,velocity2,velocity3
	real, dimension(Nvar) :: vals

	!Grid Parameters
	integer,intent(inout) :: header1,header2,header3
	integer,dimension(counter0_max),intent(out) :: counter0
	integer,dimension(counter1_max),intent(out) :: counter1
	integer,dimension(counter2_max),intent(out) :: counter2
	integer,dimension(counter3_max),intent(out) :: counter3

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
	double precision :: min_rmsd
        integer :: TimeA,TimeB,step

        !Initialize the scene
        call InitialSetup3(velocity1,velocity2,velocity3,&
                           coords_gradient(1:3),coords_gradient(4:6),coords_gradient(7:9),&
                   	   initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
			   initial_bond_angle1,initial_bond_angle2)


        !Accelerate the velcocities for a half step (verlet)

	call getVar3(coords_gradient(1:Ncoords),Natoms,vals(1))
	call getVar4(coords_gradient(1:Ncoords),Natoms,vals(2))

        call Acceleration(vals,coords_gradient,&
             AccelerationConstant0,AccelerationConstant1,AccelerationConstant2)

        open(filechannel1,file=path4//trajectoryfile)
        write(filechannel1,'(I1)') 3
        write(filechannel1,*) ""
        write(filechannel1,'(A1,3F10.6)') 'H',&
              coords_gradient(1), coords_gradient(2), coords_gradient(3)
        write(filechannel1,'(A1,3F10.6)') 'H',&
              coords_gradient(4), coords_gradient(5), coords_gradient(6)
        write(filechannel1,'(A1,3F10.6)') 'H',&
              coords_gradient(7), coords_gradient(8), coords_gradient(9)
 	close(filechannel1)
  
	velocity1 = velocity1 + 0.5*coords_gradient(Ncoords+1:Ncoords+3)
	velocity2 = velocity2 + 0.5*coords_gradient(Ncoords+4:Ncoords+6)
	velocity3 = velocity3 + 0.5*coords_gradient(Ncoords+7:Ncoords+9)

        !Keep track of the time
        TimeA = time()
open(filechannel2,file=path4//checkstatefile)
        do step = 1, Nsteps

		!Every 50 frames, print to an xyz file for visualization
                if (modulo(step,50) == 0) then
                        open(filechannel1,file=path4//trajectoryfile,position="append")
                        write(filechannel1,'(I1)') 3
                        write(filechannel1,*) ""
                        write(filechannel1,'(A1,3F10.6)') 'H',&
                                coords_gradient(1), coords_gradient(2), coords_gradient(3)
                        write(filechannel1,'(A1,3F10.6)') 'H',&
                                coords_gradient(4), coords_gradient(5), coords_gradient(6)
                        write(filechannel1,'(A1,3F10.6)') 'H',&
                                coords_gradient(7), coords_gradient(8), coords_gradient(9)
 			close(filechannel1)
                end if

		!Get the variables
		call getVar3(coords_gradient(1:Ncoords),Natoms,vals(1))
		call getVar4(coords_gradient(1:Ncoords),Natoms,vals(2))
 
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
		coords_gradient(1:3) = coords_gradient(1:3) + dt*velocity1
		coords_gradient(4:6) = coords_gradient(4:6) + dt*velocity2
		coords_gradient(7:9) = coords_gradient(7:9) + dt*velocity3

		!Check for similar frames
		call checkState(coords_gradient(1:Ncoords),closestCoords,min_rmsd,&
				counter0,counter1,counter2,counter3)

                !Update the gradient
                call Acceleration(vals,coords_gradient,&
			AccelerationConstant0,AccelerationConstant1,AccelerationConstant2)

write(filechannel2,FMT="(I8,1x,E20.4,1x,F9.6,1x,F9.6)") step, min_rmsd, vals(1), vals(2)

		!Update the velocities
		velocity1 = velocity1 + coords_gradient(Ncoords+1:Ncoords+3)
		velocity2 = velocity2 + coords_gradient(Ncoords+4:Ncoords+6)
		velocity3 = velocity3 + coords_gradient(Ncoords+7:Ncoords+9)

        end do
close(filechannel2)

end subroutine addTrajectory


subroutine Acceleration(vals,coords_gradient,&
		AccelerationConstant0,AccelerationConstant1,AccelerationConstant2)

        use f2_physics_parameters
	use f2_parameters
	use f2_variables
        implicit none
	real, dimension(6*Natoms), intent(inout) :: coords_gradient
	real, dimension(6*Natoms) :: closestCoords
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
        coords_distance = coords_gradient(4:6) - coords_gradient(1:3)
        velocity_change = coords_distance * AccelerationConstant2 * ( &
                        stretch_factor - 1.0 ) * stretch_factor / distance12
        coords_gradient(Ncoords+1:Ncoords+3) = velocity_change
        coords_gradient(Ncoords+4:Ncoords+6) = -velocity_change

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Morse Potential Derivative of Atoms 1 and 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!Save some calculation with this factor
	stretch_factor = exp(Morsealpha_hydrogen*(distance13-Morser0_hydrogen))

	!Calculate the velocity change
        coords_distance = coords_gradient(7:9) - coords_gradient(1:3)
        velocity_change = coords_distance * AccelerationConstant2 * ( &
                        stretch_factor - 1.0 ) * stretch_factor / distance13
        coords_gradient(Ncoords+1:Ncoords+3) = coords_gradient(Ncoords+1:Ncoords+3) + velocity_change
        coords_gradient(Ncoords+7:Ncoords+9) = -velocity_change

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Harmonic Oscillator of Atoms 2 and 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!Calculate the distance
        coords_distance = coords_gradient(7:9) - coords_gradient(4:6)
        distance_squared = coords_distance(1)**2 + coords_distance(2)**2 +&
                           coords_distance(3)**2
        distance23 = sqrt(distance_squared)

	!Calculate the velocity change
        velocity_change = coords_distance * ( &
                        AccelerationConstant0 + &
                        AccelerationConstant1/distance23 )
        coords_gradient(Ncoords+4:Ncoords+6) = coords_gradient(Ncoords+4:Ncoords+6) + velocity_change
        coords_gradient(Ncoords+7:Ncoords+9) = coords_gradient(Ncoords+7:Ncoords+9) - velocity_change

	return

end subroutine Acceleration



subroutine InitialSetup3(velocity1,velocity2,velocity3,&
                         coords1,coords2,coords3, &
                   	 initial_bond_distance,initial_rotational_speed,initial_rotation_angle,&
			 initial_bond_angle1,initial_bond_angle2)
        use f2_physics_parameters
	use f2_parameters
	use f1_functions
        implicit none
        real, dimension(3), intent(out) :: velocity1, velocity2, velocity3
        real, dimension(3), intent(out) :: coords1, coords2, coords3
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
        coords1 = (/ 0.0, 0.0, 0.0 /)
        coords2 = (/ collision_distance, collision_skew, 0.0 /) + bond_vector/2
        coords3 = (/ collision_distance, collision_skew, 0.0 /) - bond_vector/2

	!Calculate the direction the bond is spinning given the two
	!orthogonal vectors
	rotation_vector = ( sin(initial_rotation_angle) * rotation0_vector + &
                            cos(initial_rotation_angle) * rotation90_vector ) * &
                          initial_rotational_speed

	!With the translational, vibrational, and rotational vectors and speeds
	!Calculate the velocities of the atoms
        velocity1 = (/ sqrt(2*initial_translational_KE/mass_hydrogen), 0.0, 0.0 /)
        velocity2 = rotation_vector
        velocity3 = -rotation_vector

open(progresschannel,file=path4//progressfile,position="append")
write(progresschannel,*) ""
write(progresschannel,*) "The parameters (in reduced units) are: "
write(progresschannel,*) ""
write(progresschannel,*) "Collision H - H2 Distance: ", collision_distance
write(progresschannel,*) "Initial H Velocity: ", velocity1(1)
write(progresschannel,*) ""
write(progresschannel,*) "Initial H2 Bond distance: ", initial_bond_distance
write(progresschannel,*) "Initial H2 Rotation speed: ", initial_rotational_speed
write(progresschannel,*) "Initial H2 Rotation angle: ", initial_rotation_angle
close(progresschannel)

end subroutine InitialSetup3


real function MorsePotential(coords1,coords2)
        use f2_physics_parameters
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
        use f2_physics_parameters
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
        use f2_physics_parameters
        implicit none
        real, dimension(3), intent(in) :: velocity
        real :: speed_squared

	!Compute the speed, then compute the kinetic energy
        speed_squared = velocity(1)**2 + velocity(2)**2 + velocity(3)**2
        KineticEnergy = 0.5*mass_hydrogen*speed_squared

end function KineticEnergy



end module makeTrajectory2
