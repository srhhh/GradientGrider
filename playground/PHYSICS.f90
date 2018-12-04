module PHYSICS
use DOUBLE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 PHYSICAL CONSTANTS (all in SI)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(dp),parameter :: Na = 6.02214d23
real(dp),parameter :: kb = 1.38065d-23
real(dp),parameter :: Da = 1.66054d-27
real(dp),parameter :: c  = 2.99792e8
real(dp),parameter :: pi = 3.14159d0
real(dp),parameter :: pi2 = pi*2.0d0
real(dp),parameter :: eV = 1.60218d-19 
real(dp),parameter :: hbar = 1.05457d-34

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     REDUCED UNITS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(dp),parameter :: RU_length = 1.0d-10
real(dp),parameter :: RU_mass = Da
real(dp),parameter :: RU_time = 1.0d-15
real(dp),parameter :: RU_force = RU_mass*RU_length*RU_time**(-2)
real(dp),parameter :: RU_energy = RU_mass*(RU_length/RU_time)**2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    TIME PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(dp),parameter :: dt = (1.0d-17)/RU_time
integer,parameter :: Nsteps = 30000

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   ATOMIC PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(dp),parameter :: HOr0_hydrogen = (1.0d-10)*0.7412d0/RU_length
                                        ! originally in A
real(dp),parameter :: HOke_hydrogen =  (1.0d10*1.0d-8)*5.756d0/(RU_force/RU_length)
                                        ! originally in mdyn/A= 1.0e-8 N/A
real(dp),parameter :: Morser0_hydrogen = (1.0d-10)*3.25d0/RU_length
                                        ! originally in A
real(dp),parameter :: MorseDe_hydrogen = (1.0d0/(Na*2.39006d-4))*0.02347d0/RU_energy
                                        ! originally in kcal/mol
real(dp),parameter :: Morsealpha_hydrogen = -10.6d0
                                        ! unitless
real(dp),parameter :: mass_hydrogen = (.001d0/Na)*(1.00794d0)/RU_mass
                                        !originally g/mol

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 TEMPERATURE PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!real(dp), parameter :: temperature = 200.0d0
!					!Kelvin
!real(dp), parameter :: upsilon_max = 5.0d0
!real(dp), parameter :: frequency_factor  = sqrt(HOke_hydrogen/(0.5d0*mass_hydrogen))/pi2
!real(dp), parameter :: upsilon_factor1 = -pi2*(hbar/(RU_energy*RU_time))*frequency_factor/(kb*temperature/RU_energy)
!real(dp), parameter :: upsilon_factor2 = pi2*(hbar/(RU_energy*RU_time))*frequency_factor
!real(dp), parameter :: temperature_factor = exp(0.5d0*upsilon_factor1)
!real(dp), parameter :: temperature_scaling = (1.0d0-temperature_factor)*upsilon_max/temperature_factor

real(dp), parameter :: temperature = 300.0d0
					!Kelvin
real(dp), parameter :: upsilon_max = 5.0d0
					!Vibrational Quantum Number Cutoff
real(dp), parameter :: vib_frequency = sqrt(HOke_hydrogen/(0.5d0*mass_hydrogen))/pi2
integer, parameter  :: vib_period = floor((1.0 / vib_frequency) / dt)
real(dp), parameter :: theta_vib = (hbar/(RU_energy*RU_time)) * pi2 * HOke_hydrogen / (kb / RU_energy)
					!Vibrational Constant
real(dp), parameter :: upsilon_factor1 = theta_vib / temperature
real(dp), parameter :: upsilon_factor2 = 1.0d0 - exp(-upsilon_factor1)
real(dp), parameter :: epsilon_factor = (hbar/(RU_energy*RU_time)) * pi2 * vib_frequency
					!Some scaling factors
real(dp), parameter :: J_max = 10.0d0
					!Rotational Quantum Number Cutoff
real(dp), parameter :: J_factor1 = (hbar**2 / (kb * temperature * mass_hydrogen)) / ((RU_energy * RU_time**2))
real(dp), parameter :: J_factor2 = (hbar**2 / (mass_hydrogen)) / (((RU_energy * RU_time)**2))
					!Some scaling factors


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  COLLISION PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(dp),parameter :: initial_translational_KE = (1.0d0)*eV/RU_energy
!real(dp),parameter :: collision_distance = Morser0_hydrogen*5.0
 real(dp),parameter :: collision_distance = 8.5d0
real(dp),parameter :: collision_skew = HOr0_hydrogen*0.00d0

!For the COLLISION_DATA:
!Each atom index corresponds to a zero or a one
!A one indicates it is 'incoming'
!A zero indicates it is 'to-be-collided'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			BONDAGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!integer,parameter :: Nbonds = 1
!integer,dimension(Nbonds,2),parameter :: BONDING_DATA = reshape((/ 2, &
!                                                                   3 /),&
!                                                                                  (/ Nbonds, 2 /))
!integer,parameter :: Nvar_eff = 2
!integer,dimension(Nvar_eff,3),parameter :: BONDING_VALUE_DATA = reshape((/ 1, 1,    &
!                                                                           2, 3,    &
!                                                                           1, 2/),  &
!                                                                                  (/ Nvar_eff, 3 /))
!integer,dimension(3),parameter :: COLLISION_DATA = (/ 1, 0, 0 /)

!integer,parameter :: Nbonds = 2
!integer,dimension(Nbonds,2),parameter :: BONDING_DATA = reshape((/ 1, 3,   &
!                                                                   2, 4 /),        (/ Nbonds, 2 /))
!
!integer,parameter :: Nvar_eff = 2
!integer,dimension(Nvar_eff,3),parameter :: BONDING_VALUE_DATA = reshape((/ 1, 2,    &
!                                                                           3, 4,    &
!                                                                           1, 2 /), &
!                                                                                   (/ Nvar_eff, 3 /)) 
!integer,dimension(4),parameter :: COLLISION_DATA = (/ 1, 1, 0, 0 /)

!!H2 - H2 COLLISION STUFF
!integer,parameter :: Nbonds = 2
!integer,dimension(Nbonds,2),parameter :: BONDING_DATA = reshape((/ 1, 3,   &
!                                                                   2, 4 /),        (/ Nbonds, 2 /))
!
!integer,parameter :: Nvar_eff = 0
!integer,dimension(Nvar_eff,3) :: BONDING_VALUE_DATA
!integer,dimension(4),parameter :: COLLISION_DATA = (/ 1, 1, 0, 0 /)
!
!integer,parameter :: Nindistinguishables = 8
!integer,dimension(8,4),parameter :: INDISTINGUISHABLES = reshape((/ 1, 2, 1, 2, 3, 4, 3, 4, &
!							            2, 1, 2, 1, 4, 3, 4, 3, &
!							            3, 3, 4, 4, 1, 1, 2, 2, &
!							            4, 4, 3, 3, 2, 2, 1, 1 /), 	(/ 8, 4 /))
!integer,dimension(4) :: BOND_LABELLING_DATA = (/ 1, 2, 3, 4 /)
!integer,dimension(2,4),parameter :: BOND_LABELS_TENTATIVE = reshape((/ 1, 2, 1, 2, &
!								        3, 4, 4, 3 /),    (/2, 4 /))
!
!!H - H2 COLLISION STUFF
integer,parameter :: Nbonds = 1
integer,dimension(Nbonds,2),parameter :: BONDING_DATA = reshape((/ 2, &
                                                                   3 /),&
                                                                                  (/ Nbonds, 2 /))
integer,parameter :: Nvar_eff = 0
integer,dimension(Nvar_eff,3) :: BONDING_VALUE_DATA
integer,dimension(3),parameter :: COLLISION_DATA = (/ 1, 0, 0 /)

integer,parameter :: Nindistinguishables = 2
integer,dimension(2,3),parameter :: INDISTINGUISHABLES = reshape((/ 1, 1, &
								    2, 3, &
								    3, 2 /),    (/ 2, 3 /))
integer,dimension(3) :: BOND_LABELLING_DATA = (/ 1, 2, 3 /)
!
!integer,parameter :: Nindistinguishables = 1
!integer,dimension(1,3),parameter :: INDISTINGUISHABLES = reshape((/ 1, &
!								    2, &
!								    3 /),    (/ 1, 3 /))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                RANDOM COLLISION SETTINGS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(dp),dimension(6,Nbonds) :: INITIAL_BOND_DATA

real(dp) :: vibrational_max
real(dp) :: translational_max
real(dp) :: rotational_max
real(dp) :: rovibrational_max
real(dp) :: TranslationalEnergy_max

real :: max_absenergychange, min_absenergychange
real :: max_relenergychange, min_relenergychange
real :: max_rotenergychange, min_rotenergychange
real :: abs_energychange, rel_energychange, rot_energychange
real :: max_TranslationalEnergy

real(dp) :: sizeEnergyBin,sizeAngleBin
real(dp) :: sizeAbsEnergyBin,sizeRelEnergyBin,sizeRotEnergyBin
real(dp) :: sizeDeltaEnergyBin
integer,parameter :: energyBins = 100                  !This is for the SA heat map
integer,parameter :: angleBins = 200                   !Same
integer,parameter :: scatteringangleBins = 50          !This is for the regular histogram
                                                       !  (Note: make this a divisor of angleBins)
integer,parameter :: energychangeBins = 50             !Same



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



contains



subroutine decompose_two_velocities(coords,velocities,velocity_translation,TRVenergies)

implicit none
real(dp),dimension(3,2),intent(in) :: coords,velocities
real(dp),dimension(3),intent(out) :: TRVenergies
real(dp),dimension(3),intent(out) :: velocity_translation
real(dp),dimension(3) :: velocity_vibration,velocity_rotation
real(dp),dimension(3) :: parallel_vector,orthogonal_vector
real(dp),dimension(3) :: velocity_translation1,velocity_translation2,velocity_translation3
real(dp),dimension(3) :: velocity_rotation1,velocity_rotation2
real(dp),dimension(3,2) :: velocities_parallel
real(dp),dimension(3,2) :: velocities_orthogonal
real(dp),dimension(3,2) :: velocities_ortho2
real(dp),dimension(3,2) :: velocities_ortho1

parallel_vector = coords(:,1) - coords(:,2)
parallel_vector = parallel_vector / sqrt(sum(parallel_vector**2))

velocities_parallel(:,1) = dot_product(velocities(:,1),&
                                      parallel_vector) * parallel_vector
velocities_parallel(:,2) = dot_product(velocities(:,2),&
                                       parallel_vector) * parallel_vector

velocities_orthogonal = velocities - velocities_parallel

orthogonal_vector = velocities_orthogonal(:,1) + velocities_orthogonal(:,2)
if (all(orthogonal_vector /= 0.0d0)) orthogonal_vector = orthogonal_vector / sqrt(sum(orthogonal_vector**2))

velocities_ortho1(:,1) = dot_product(velocities_orthogonal(:,1),&
                                     orthogonal_vector) * orthogonal_vector
velocities_ortho1(:,2) = dot_product(velocities_orthogonal(:,2),&
                                     orthogonal_vector) * orthogonal_vector

velocities_ortho2 = velocities_orthogonal - velocities_ortho1

velocity_translation1 = 0.5 * (velocities_parallel(:,1) + velocities_parallel(:,2))
velocity_vibration = velocities_parallel(:,1) - velocity_translation1

velocity_translation2 = 0.5 * (velocities_ortho1(:,1) + velocities_ortho1(:,2))
velocity_rotation1 = velocities_ortho1(:,1) - velocity_translation2
velocity_translation3 = 0.5 * (velocities_ortho2(:,1) + velocities_ortho2(:,2))
velocity_rotation2 = velocities_ortho2(:,1) - velocity_translation3

velocity_translation = velocity_translation1 + &
                       velocity_translation2 + &
                       velocity_translation3
velocity_rotation = velocity_rotation1 + &
                    velocity_rotation2

!Note: these are INDIVIDUAL VELOCITIES
!So first we multiply by two and take a square root to get a speed
!Then we divide by two and take a square to get a kinetic energy
!These two calculations cancel each other out
TRVenergies(1) = sum(velocity_translation**2)
TRVenergies(2) = sum(velocity_rotation**2)
TRVenergies(3) = sum(velocity_vibration**2)

TRVenergies = mass_hydrogen*TRVenergies
!Note: now they are TOTAL KINETIC ENERGIES

end subroutine decompose_two_velocities

subroutine decompose_three_velocities(velocities,rel_velocities,velocity_center_of_mass)
implicit none
real(dp), dimension(9),intent(in) :: velocities
real(dp), dimension(9),intent(out) :: rel_velocities
real(dp), dimension(3),intent(out) :: velocity_center_of_mass

velocity_center_of_mass = (velocities(1:3) + velocities(4:6) +&
                           velocities(7:9)) / 3
rel_velocities(1:3) = velocities(1:3) - velocity_center_of_mass
rel_velocities(4:6) = velocities(4:6) - velocity_center_of_mass
rel_velocities(7:9) = velocities(7:9) - velocity_center_of_mass

end subroutine decompose_three_velocities


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     CROSS PRODUCT FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   INPUT:      real, dim (3) A                 "vectorA"
!               real, dim (3) B                 "vectorB"
!   OUTPUT:     real, dim (3) AcrossB           "the cross product"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cross(A,B,AcrossB)
implicit none
real(dp),dimension(3),intent(in) :: A, B
real(dp),dimension(3),intent(out) :: AcrossB

AcrossB = (/ A(2)*B(3) - A(3)*B(2), &
             A(3)*B(1) - A(1)*B(3), &
             A(1)*B(2) - A(2)*B(1)  /)

end subroutine cross


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     ENERGY, FORCE, ACCELERATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(dp) function KineticEnergy(velocity)
        implicit none
        real(dp), dimension(3), intent(in) :: velocity
        real(dp) :: speed_squared

        !Compute the speed, then compute the kinetic energy
        speed_squared = velocity(1)**2 + velocity(2)**2 + velocity(3)**2
        KineticEnergy = 0.5d0*mass_hydrogen*speed_squared

end function KineticEnergy

real(dp) function RotationalEnergy(coords1,coords2,velocity1,velocity2)
	implicit none
	real(dp),dimension(3),intent(in) :: coords1,coords2,velocity1,velocity2
	real(dp),dimension(3) :: CM_vector, bond_vector, angular_momentum1, angular_momentum2

	CM_vector = 0.5 * (coords1 + coords2)
	bond_vector = coords1 - CM_vector

	call cross(bond_vector,velocity1,angular_momentum1)
	call cross(-bond_vector,velocity2,angular_momentum2)

	RotationalEnergy = 0.25 * mass_hydrogen * sum((angular_momentum1 + &
                           angular_momentum2)**2) / (sum(bond_vector**2))

end function RotationalEnergy

subroutine VelocityBondAngle(coords1,coords2,velocity1,velocity2,theta1,theta2)
	implicit none
	real(dp),dimension(3),intent(in) :: coords1, coords2, velocity1, velocity2
	real(dp),intent(out) :: theta1, theta2
	real(dp),dimension(3) :: CM_vector, bond_vector
	real(dp) :: bond_length

	CM_vector = 0.5 * (coords1 + coords2)
	bond_vector = coords1 - CM_vector
	bond_length = sum(bond_vector**2)

	theta1 = acos(dot_product(bond_vector,velocity1) / sqrt(bond_length * sum(velocity1**2)))
	theta2 = acos(dot_product(-bond_vector,velocity2) / sqrt(bond_length * sum(velocity2**2)))

end subroutine VelocityBondAngle


real(dp) function MorsePotential(coords1,coords2)
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


real(dp) function HOPotential(coords1,coords2)
        implicit none
        integer :: i
        real(dp), dimension(3), intent(in) :: coords1,coords2
        real(dp), dimension(3) :: coords_distance
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



subroutine InitialSetup3(coords,velocities)
	use PARAMETERS
        implicit none
        real(dp), dimension(3,Natoms), intent(out) :: velocities,coords
        real(dp), dimension(3) :: bond_vector,rotation_vector,rotation0_vector,rotation90_vector
        real(dp) :: initial_bond_distance,initial_rotational_speed,initial_rotation_angle
        real(dp) :: initial_bond_angle1,initial_bond_angle2
	real(dp) :: incoming_speed
        integer :: i, atom1, atom2

	coords = 0.0d0
	velocities = 0.0d0

        do i = 1, Nbonds
                initial_bond_distance = INITIAL_BOND_DATA(1,i)
                initial_rotational_speed = INITIAL_BOND_DATA(2,i)
                initial_rotation_angle = INITIAL_BOND_DATA(3,i)
                initial_bond_angle1 = INITIAL_BOND_DATA(4,i)
                initial_bond_angle2 = INITIAL_BOND_DATA(5,i)

                atom1 = BONDING_DATA(i,1)
                atom2 = BONDING_DATA(i,2)

                !Figure out which direction the bond is going
                bond_vector = (/ sin(initial_bond_angle1)*cos(initial_bond_angle2),&
                                 sin(initial_bond_angle1)*sin(initial_bond_angle2),&
                                 cos(initial_bond_angle1) /)

                !And two orthogonal vectors that are both perpendicular to the bond
                rotation0_vector =  (/ bond_vector(2), -bond_vector(1), 0.0d0 /) / &
                                       sqrt(bond_vector(2)**2 + bond_vector(1)**2)
                call cross(bond_vector,rotation0_vector,rotation90_vector)

                !Calculate the coordinates of the atoms via the distance, skew, and bond vector
                bond_vector = initial_bond_distance * bond_vector / 2
                coords(:,atom1) = bond_vector
                coords(:,atom2) = -bond_vector

                !Calculate the direction the bond is spinning given the two
                !orthogonal vectors
                rotation_vector = ( sin(initial_rotation_angle) * rotation0_vector + &
                                    cos(initial_rotation_angle) * rotation90_vector ) * &
                                  initial_rotational_speed

                !With the translational, vibrational, and rotational vectors and speeds
                !Calculate the velocities of the atoms
                velocities(:,atom1) = rotation_vector
                velocities(:,atom2) = -rotation_vector
        end do

	incoming_speed = sqrt(2*initial_translational_KE/(mass_hydrogen*sum(COLLISION_DATA)))

	do i = 1, Natoms
		if (COLLISION_DATA(i) == 1) then
			velocities(:,i) = velocities(:,i) + &
                        (/ incoming_speed, 0.0d0, 0.0d0 /)
		else
			coords(:,i) = coords(:,i) + &
			(/ collision_distance, collision_skew, 0.0d0 /)
		end if
	end do

end subroutine InitialSetup3


subroutine InitialSetup4(coords,velocities)
	use PARAMETERS
        implicit none
        real(dp), dimension(3,Natoms), intent(out) :: velocities,coords
        real(dp), dimension(3) :: bond_vector,rotation_vector,rotation0_vector,rotation90_vector
        real(dp) :: random_num1,random_num2,random_num3,random_r2,random_r3
        real(dp) :: initial_bond_distance, initial_rotation_angle, initial_rotational_speed
        real(dp) :: initial_bond_angle1, initial_bond_angle2
        real(dp) :: initial_energy_H2,initial_vibrational_energy,initial_rotational_energy
        real(dp) :: bond_period_elapsed
        real(dp) :: probJ_max, J_factor3
	real(dp) :: incoming_speed
        integer :: i, atom1, atom2


	coords = 0.0d0
	velocities = 0.0d0

        do i = 1, Nbonds

                !The orientation of the H2 should be some random
                !point in the unit sphere
                do
                        !First get a random point in the unit cube
                        !centered at zero
                        random_num1 = rand() - 0.5d0
                        random_num2 = rand() - 0.5d0
                        random_num3 = rand() - 0.5d0
                        random_r2 = random_num1**2 + random_num2**2
                        random_r3 = random_r2 + random_num3**2

                        !If the point lies outside of the cube, reject it
                        if (random_r3 > 0.25d0) cycle
                        random_r2 = sqrt(random_r2)

                        !But if it lies in the sphere, use its direction (angles)
                        initial_bond_angle1 = atan2(random_num1,random_num2)
                        initial_bond_angle2 = atan2(random_r2,random_num3)
                        exit
                end do

                !The vibrational energy of the H2 should be some random value
                !that follows the boltzmann distribution at this temperature
                do
                        !This picks a random value between zero and some very high upper limit
                        random_num1 = rand() * upsilon_max
!                               random_num2 = rand() * upsilon_factor2
                        random_num2 = rand()

!                               if (exp(-random_num1 * upsilon_factor1) < random_num2) cycle
                        if (exp(-random_num1 * upsilon_factor1) * upsilon_factor2 < random_num2) cycle

                        initial_vibrational_energy = (random_num1 + 0.5d0) * epsilon_factor
                        exit
                end do
                initial_bond_distance = HOr0_hydrogen + sqrt(initial_vibrational_energy*2/HOke_hydrogen)
                J_factor3 = J_factor1 / (initial_bond_distance**2)
                probJ_max = sqrt(2*J_factor3) * exp(J_factor3*0.25d0 - 0.5d0)
!                J_factor3 = 85.3d0 / temperature

                !The rotational energy of the H2 should be some random value
                !that follows the boltzmann distribution at this temperature
                do
                        !This picks a random value between zero and some very high upper limit
                        random_num1 = rand() * J_max
                        random_num2 = rand() * probJ_max

                        if ((2*random_num1 + 1.0d0) * J_factor3 * exp(-random_num1 * (random_num1 + 1.0d0) * &
                            J_factor3) < random_num2) cycle

                        initial_rotational_energy = (random_num1) * (random_num1 + 1.0d0) * J_factor2
                        exit
                end do

                random_num1 = rand()
                initial_rotational_speed = sqrt(initial_rotational_energy/mass_hydrogen)
                initial_rotation_angle = random_num1*pi2 - pi
                bond_period_elapsed = rand()

                !Now we actually change the physical location of coordinates to reflect this bonding
                atom1 = BONDING_DATA(i,1)
                atom2 = BONDING_DATA(i,2)

                !Figure out which direction the bond is going
                bond_vector = (/ sin(initial_bond_angle1)*cos(initial_bond_angle2),&
                                 sin(initial_bond_angle1)*sin(initial_bond_angle2),&
                                 cos(initial_bond_angle1) /)

                !And two orthogonal vectors that are both perpendicular to the bond
                rotation0_vector =  (/ bond_vector(2), -bond_vector(1), 0.0d0 /) / &
                                       sqrt(bond_vector(2)**2 + bond_vector(1)**2)
                call cross(bond_vector,rotation0_vector,rotation90_vector)

                !Calculate the coordinates of the atoms via the distance, skew, and bond vector
                bond_vector = initial_bond_distance * bond_vector / 2
                coords(:,atom1) = bond_vector
                coords(:,atom2) = -bond_vector

                !Calculate the direction the bond is spinning given the two
                !orthogonal vectors
                rotation_vector = ( sin(initial_rotation_angle) * rotation0_vector + &
                                    cos(initial_rotation_angle) * rotation90_vector ) * &
                                  initial_rotational_speed

                !With the translational, vibrational, and rotational vectors and speeds
                !Calculate the velocities of the atoms
                velocities(:,atom1) = rotation_vector
                velocities(:,atom2) = -rotation_vector
        end do

	incoming_speed = sqrt(2*initial_translational_KE/(mass_hydrogen*sum(COLLISION_DATA)))

	do i = 1, Natoms
		if (COLLISION_DATA(i) == 1) then
			velocities(:,i) = velocities(:,i) + &
                        (/ incoming_speed, 0.0d0, 0.0d0 /)
		else
			coords(:,i) = coords(:,i) + &
			(/ collision_distance, collision_skew, 0.0d0 /)
		end if
	end do

end subroutine InitialSetup4





subroutine NonBondedForce(coords1,coords2,gradient1,gradient2,r)
        implicit none
        real(dp), dimension(3), intent(in) :: coords1,coords2
        real(dp), dimension(3), intent(inout) :: gradient1,gradient2
        real(dp),optional,intent(in) :: r
        real(dp), dimension(3) :: coords_distance,velocity_change
        real(dp) :: distance12,distance_squared,stretch_factor

        coords_distance = coords2 - coords1
        if (present(r)) then
                distance12 = r
        else
                distance_squared = coords_distance(1)**2 &
                                 + coords_distance(2)**2 &
                                 + coords_distance(3)**2
                distance12 = sqrt(distance_squared)
        end if

        stretch_factor = exp(Morsealpha_hydrogen*(distance12-Morser0_hydrogen))
        velocity_change = coords_distance * AccelerationConstant2 * ( &
                        stretch_factor - 1.0d0 ) * stretch_factor / distance12

        gradient1 = gradient1 + velocity_change
        gradient2 = gradient2 - velocity_change

        return

end subroutine NonBondedForce

subroutine BondedForce(coords1,coords2,gradient1,gradient2,r)
        implicit none
        real(dp), dimension(3), intent(in) :: coords1,coords2
        real(dp), dimension(3), intent(inout) :: gradient1,gradient2
        real(dp),optional,intent(in) :: r
        real(dp), dimension(3) :: coords_distance,velocity_change
        real(dp) :: distance12,distance_squared

        coords_distance = coords2 - coords1
        if (present(r)) then
                distance12 = r
        else
                distance_squared = coords_distance(1)**2 &
                                 + coords_distance(2)**2 &
                                 + coords_distance(3)**2
                distance12 = sqrt(distance_squared)
        end if

        velocity_change = coords_distance * ( &
                        AccelerationConstant0 + &
                        AccelerationConstant1/distance12 )

        gradient1 = gradient1 + velocity_change
        gradient2 = gradient2 - velocity_change

        return

end subroutine BondedForce






end module PHYSICS



