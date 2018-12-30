!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       MODULE
!               PHYSICS
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

real(dp),parameter :: HOr0_hydrogen = (1.0d-10)*0.7412d0/RU_length                         !RU_length
                                        ! originally in A
real(dp),parameter :: HOke_hydrogen =  (1.0d10*1.0d-8)*5.756d0/(RU_force/RU_length)        !RU_force/RU_length
                                        ! originally in mdyn/A= 1.0e-8 N/A
real(dp),parameter :: Morser0_hydrogen = (1.0d-10)*3.25d0/RU_length                        !RU_length
                                        ! originally in A
real(dp),parameter :: MorseDe_hydrogen = (1.0d0/(Na*2.39006d-4))*0.02347d0/RU_energy       !RU_energy
                                        ! originally in kcal/mol
real(dp),parameter :: Morsealpha_hydrogen = -10.6d0                                        ! -
                                        ! unitless
real(dp),parameter :: mass_hydrogen = (.001d0/Na)*(1.00794d0)/RU_mass                      !RU_mass
                                        !originally g/mol

        !Calculation Savers
        real(dp), parameter :: PotentialConstant0 =&                                       !RU_force*RU_length
                0.5d0*HOke_hydrogen*HOr0_hydrogen**2
        real(dp), parameter :: PotentialConstant1 =&                                       !RU_force
                -HOke_hydrogen*HOr0_hydrogen
        real(dp), parameter :: PotentialConstant2 =&                                       !RU_force/RU_length
                0.5d0*HOke_hydrogen
        real(dp), parameter :: AccelerationConstant0 =&                                    !RU_length/RU_time
                HOke_hydrogen*(dt/mass_hydrogen)
        real(dp), parameter :: AccelerationConstant1 =&                                    !RU_length^2/RU_time
                -HOke_hydrogen*(dt/mass_hydrogen)*HOr0_hydrogen
        real(dp), parameter :: AccelerationConstant2 =&                                    !RU_length^2/RU_time
                2*MorseDe_hydrogen*Morsealpha_hydrogen*dt/mass_hydrogen


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 TEMPERATURE PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!The temperature used to sample vibration and rotation (not the total temperature)
real(dp), parameter :: temperature = 300.0d0                                       !Kelvin

!The reduced mass of diatomic hydrogen
!
! reduced_mass = (mass1) (mass2) / (mass1 + mass2)
!
!        mass1 = mass2   for diatomics
!
! reduced_mass = (mass1) (mass1) / (mass1 + mass1)
!              = (1/2) (mass1)
!
real(dp), parameter :: reduced_mass = 0.5d0 * mass_hydrogen                        !RU_mass






!VIBRATION

!The vibrational frequency of diatomic hydrogen (based off of the force constant)
!McQuarrie, 1997 : section 18-4 : pg.740
!
! vib_frequency = sqrt (force constant / reduced mass) / 2 pi
!
real(dp), parameter :: vib_frequency = sqrt(HOke_hydrogen/reduced_mass)/pi2        !1/RU_time

!The vibrational period of the diatomic hydrogen, rounded down to the nearest timestep
!
! vib_period = 1 / vib_frequency                       (real number)
!
! vib period = floor(vib_period / dt)                  (integer)
!            = floor((1.0 / vib_frequency) / dt)
!
integer, parameter  :: vib_period = floor((1.0d0 / vib_frequency) / dt)            !RU_time

!The vibrational temperature of diatomic hydrogen
!McQuarrie, 1997 : section 18-4 : pg 740
!
! theta_vib = (h) (vib_frequency) / kb
!           = (hbar) (2) (pi) (vib_frequency) / kb
!
real(dp), parameter :: theta_vib = &
        (hbar/(RU_energy*RU_time)) * pi2 * vib_frequency / (kb / RU_energy)        !Kelvin

!A few variables are calculated ahead of time to save probability calculations:
!
! upsilon_factor1 = theta_vib / T
! upsilon_factor2 = 1 - e^(-theta_vib / T)
!
!The probability of a certain vibrational number (upsilon) existing is given
!in McQuarrie, 1997 : section 18-4 : pg 742
!
! P(upsilon) = (1 - e^(-theta_vib/T)) e^(-(upsilon) (theta_vib) / T)
!
! P(upsilon) = (1 - e^(-upsilon_factor1)) e^(-(upsilon) (upsilon_factor1))
!            = (upsilon_factor2) e^(-(upsilon) (upsilon_factor1))
!
real(dp), parameter :: upsilon_factor1 = 7.0d7 * theta_vib / temperature                   ! -
!real(dp), parameter :: upsilon_factor2 = 1.0d0 - exp(-upsilon_factor1)             ! -
real(dp), parameter :: upsilon_factor2 = 1.0d0 - 0.0d0!exp(-upsilon_factor1)             ! -

!The Vibrational Quantum Number Cutoff
!
!In this case, we want to make sure the smallest nonzero probability
!(assumed to be 10^-16 for double precision numbers)
!is a little more than the probability of upsilon_max
!
! 10^-16 = (upsilon_factor2) e^(-(upsilon_max) (upsilon_factor1))
!
! upsilon_max > log(upsilon_factor2 / 10^-16) / upsilon_factor1
!             > log((upsilon_factor2) (10^16)) / upsilon_factor1
!             > (log(upsilon_factor2) + 16 log(10)) / upsilon_factor1
!
! upsilon_max = 1.1 * (log(upsilon_factor2) + 16 log(10)) / upsilon_factor1
!
real(dp), parameter :: upsilon_max = &                                             ! -
        1.1d0 * (log(upsilon_factor2) + 16.0d0*log(10.0d0)) / upsilon_factor1

!A variable is calculated ahead of time to save energy calculations:
!
! epsilon_factor1 = (h) (vib_frequency)
!                 = (hbar) (2) (pi) (vib_frequency)
!
!The vibrational energy is computed using the vibrational number (upsilon)
!McQuarrie, 1997 : section 18-4 : pg.740
!
! Evib(upsilon) = (upsilon + 1/2) (h) (vib_frequency)
!
!We use a slightly different classical model that does not have the
!zero point energy (ZPE) associated:
!
! Evib(upsilon) = (upsilon) (h) (vib_frequency)
!
! Evib(upsilon) = (upsilon) (epsilon_factor1)
!
real(dp), parameter :: epsilon_factor1 = &                                          !RU_energy
        (hbar/(RU_energy*RU_time)) * pi2 * vib_frequency

!ROTATION

!The Rotational Quantum Number Cutoff
real(dp), parameter :: J_max = 10.0d0

!A variable is calculated ahead of time to save energy calculations:
!
! epsilon_factor2 = (hbar)^2 / 2
!
!The rotational energy is computed using the rotational number (J)
!McQuarrie, 1997 : section 18-5 : pg.743
!
! Erot(J) = (hbar^2) (J) (J + 1) / (2) (moment_inertia)
!
! Erot(J) = (hbar^2) (J) (J + 1) / (2) (moment_inertia)
!         = ((hbar^2) / 2) (J) (J + 1) / moment_inertia
!         = (epsilon_factor2) (J) (J + 1) / moment_inertia
!
real(dp), parameter :: epsilon_factor2 = &                                         !RU_energy
        (hbar/(RU_energy*RU_time))**2




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  COLLISION PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Self-explanatory

real(dp),parameter :: initial_translational_KE = (1.0d0)*eV/RU_energy
!real(dp),parameter :: collision_distance = Morser0_hydrogen*5.0
 real(dp),parameter :: collision_distance = 8.5d0
real(dp),parameter :: collision_skew = HOr0_hydrogen*0.00d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			BONDAGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!For the COLLISION_DATA:
!Each atom index corresponds to a zero or a one
!A one indicates it is 'incoming'
!A zero indicates it is 'to-be-collided'

!For the BONDING_DATA:
!Each bond corresponds to a pair in a column
!The number corresponds to the atom index of the atom

!For the BONDING_VALUE_DATA:
!If some variable in the acceleration subroutine (particularly
!the distance between two atoms) was calculated beforehand by,
!the VARIABLES module, increment Nvar_eff and signify between
!which atoms this was by 1s and 0s

!For the INDISTINGUISHABLES:
!Each separate (but indistinguishable) labeling scheme for a frame is listed
!in a column by the order of atom indexes that would result in a similar frame

!For the BOND_LABELS_TENTATIVE:
!Each pair of nonbonded atoms corresponds to a pair in a column
!The number corresponds to the atom index of the atom


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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                RANDOM COLLISION SETTINGS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Each trajectory has a set of initial conditions; this is stored in INITIAL_BOND_DATA
!It is only accessed by the InitialSetup subroutine
!Each thread in a parallel program will have separate initial conditions
!so we insure that each has its own private copy of it

real(dp),dimension(6,Nbonds) :: INITIAL_BOND_DATA
!$OMP THREADPRIVATE(INITIAL_BOND_DATA)

!In the check-grid program we often want to analyze distributions, and to analyze those,
!we need to know the bounds of our variables
!Here we declare variables for the minimum and maximum of each variable
!This can be over-rided by user input

real :: max_absenergychange, min_absenergychange
real :: max_relenergychange, min_relenergychange
real :: max_rotenergychange, min_rotenergychange
real :: abs_energychange, rel_energychange, rot_energychange
real :: max_TranslationalEnergy

!In the check-grid program, we also need to figure out the number and sizes of
!the bins in our distribution

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



!Legacy subroutine
!((NOT USED))
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

TRVenergies(1) = sum(velocity_translation**2)
TRVenergies(2) = sum(velocity_rotation**2)
TRVenergies(3) = sum(velocity_vibration**2)

TRVenergies = mass_hydrogen*TRVenergies

end subroutine decompose_two_velocities


!Legacy subroutine
!((NOT USED))
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
!     CROSS PRODUCT SUBROUTINE
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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               InitialSetup3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!       PURPOSE
!               This subroutine initialize the coordinates and velocities of a frame
!               according to INITIAL_BOND_DATA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!               coords                          REAL(DP),DIM(3,Natoms)          The coordinates of the frame
!               velocities                      REAL(DP),DIM(3,Natoms)          The velocities of the frame
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               initial_bond_distance           REAL(DP)                        The distance between two bonded atoms
!               initial_rotational_speed        REAL(DP)                        The speed at which two atoms are rotating
!               initial_rotation_angle          REAL(DP)                        The angle the rotation vector makes wtih a fixed axis
!                                                                               with respect to the plane perpendicular to the bond
!               initial_bond_angle1             REAL(DP)                        The angle the bond makes with one axis
!               initial_bond_angle2             REAL(DP)                        The angle the bond makes with another axis
!               incoming_speed                  REAL(DP)                        The speed each incoming atom has (for collision)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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



subroutine InitialSampling3()
        use PARAMETERS
        implicit none

        !Trajectory Initialization Variables
        real(dp) :: random_num1,random_num2,random_num3,random_r2,random_r3
        real(dp) :: initial_bond_distance, initial_rotation_angle, initial_rotational_speed
        real(dp) :: initial_bond_angle1, initial_bond_angle2
        real(dp) :: initial_energy_H2,initial_vibrational_energy,initial_rotational_energy
        real(dp) :: bond_period_elapsed
        real(dp) :: Jmp, prob_Jmp, J_factor1, theta_rot, moment_inertia
        real(dp) :: upsilon, J
        
        !Incremental Integers
        integer :: m

        !Get some random initial conditions for the trajectory
        !For now, we are only handling systems with H-H bonds

        do m = 1, Nbonds

                !The orientation of the H2 should be some random
                !point in the unit sphere
                do
                        !First get a random point in the unit cube
                        !centered at zero
                        random_num1 = rand() - 0.5d0
                        random_num2 = rand() - 0.5d0
                        random_num3 = rand() - 0.5d0

                        !Calculate the length of the point (squared) and the
                        !length of its projection onto the xy plane (squared)
                        random_r2 = random_num1**2 + random_num2**2
                        random_r3 = random_r2 + random_num3**2

                        !If the point lies outside of the sphere, reject it
                        if (random_r3 > 0.25d0) cycle
                        random_r2 = sqrt(random_r2)

                        !But if it lies in the sphere, use its direction (angles)
                        initial_bond_angle1 = atan2(random_num1,random_num2)
                        initial_bond_angle2 = atan2(random_r2,random_num3)
                        exit
                end do

                !The vibrational energy of the H2 should be some random value
                !picked from a distribution at this temperature
                do
                        !First, pick a random vibrational number (upsilon)
                        !between zero and some upper limit (upsilon_max)
                        upsilon = rand() * upsilon_max

                        !And some other number between 0 and 1 for
                        !metropolis rejection
                        random_num1 = rand()

                        !If the probability that this upsilon is less than that
                        !random rejection rate, then it is rejected
                        if (exp(-upsilon * upsilon_factor1) * upsilon_factor2 < random_num1) cycle
                        exit
                end do

                !The vibrational energy is computed using the vibrational number (upsilon)
                !McQuarrie, 1997 : section 18-4 : pg.740
                !
                ! Evib(upsilon) = (upsilon + 1/2) (h) (vib_frequency)
                !
                !We use a slightly different classical model that does not have the
                !zero point energy (ZPE) associated:
                !
                ! Evib(upsilon) = (upsilon) (h) (vib_frequency)
                !
                ! Evib(upsilon) = (upsilon) (epsilon_factor)

                initial_vibrational_energy = upsilon * epsilon_factor1                 !RU_energy
!                initial_vibrational_energy = (upsilon + 0.50d0) * epsilon_factor1      !RU_energy

                !The bond distance is derived from the vibrational energy assuming all
                !of the vibrational energy is potential (so it is fully stretched out)
                !McQuarries, 1997 : section 5-1 : pg.160
                !
                !          Evib = (1/2) (force_constant) (bond_distance - r0)^2
                !
                ! bond_distance = r0 + sqrt((2) (Evib) / force_constant)
                !
                initial_bond_distance = &                                              !RU_length
                        HOr0_hydrogen + sqrt(initial_vibrational_energy*2/HOke_hydrogen)

                !The moment of inertia for diatomic hydrogen
                !
                ! moment_inertia = (mass1) (radius1)^2 + (mass2) (radius2)^2
                !
                !        mass1 = mass2  ;  radius1 = radius2 = 0.5 bond_distance   for diatomics
                !
                ! moment_inertia = (2) (mass1)((0.5) (bond_distance))^2
                !                = (1/2) (mass1) (bond_distance)^2
                !

                moment_inertia = &                                                      !RU_mass*RU_length^2
                        0.5d0 * mass_hydrogen * (initial_bond_distance)**2

                !The rotational temperature of diatomic hydrogen
                !McQuarrie, 1997 : section 18-5 : pg 744
                !
                ! theta_rot = hbar^2 / (2) (moment_inertia) (kB)
                !
                !For homonuclear diatomics, theta_rot is always multiplied by 2 so
                !we use an "effective" theta_rot
                !McQuarrie, 1997 : section 18-5 : pg 747
                !
                ! theta_rot = hbar^2 / (moment_inertia) (kB)
                !
                theta_rot = &                                                           !Kelvin
                        ((hbar/(RU_energy*RU_time))**2) / (moment_inertia * (kb / RU_energy))

                !A few variables are calculated ahead of time to save probability calculations:
                !
                ! J_factor1 = theta_rot / temperature
                !
                !The probability of a certain rotational number (J) existing is given
                !in McQuarrie, 1997 : section 18-5 : pg 745
                !
                ! P(J) = ((2)(J) + 1) (theta_rot / T) e^(-(theta_rot) (J) (J + 1) / T)
                !
                ! P(J) = ((2) (J) + 1) (J_factor1) e^(-(J_factor1) (J) (J + 1))
                !
                !
                J_factor1 = theta_rot / temperature                                     ! -

                !The rotational quantum number (J) that is most likely to occur
                !McQuarrie, 1997 : section 18-5 : pg.746
                !
                ! Jmp = sqrt(T / (2) (theta_rot)) - 0.5
                !
                ! Jmp = sqrt((0.5) (T / theta_rot)) - 0.5
                !     = sqrt( 0.5 / J_factor1) - 0.5
                !
                Jmp = &                                                                 ! -
                        sqrt(0.50d0 / J_factor1) - 0.5d0

                !The probability that the most probable rotational quantum number (Jmp) occurs
                !
                ! P(Jmp) = ((2) (Jmp) + 1) (J_factor1) e^(-(J_factor1) (Jmp) (Jmp + 1))
                !
                ! P(Jmp) = ((2) (sqrt(0.5 / J_factor1) - 0.5) + 1) (J_factor1)
                !          e^(-(J_factor1) (sqrt(0.5 / J_factor1) - 0.5)
                !                          (sqrt(0.5 / J_factor1) - 0.5 + 1)
                !        = (2) sqrt(0.5 / J_factor1) (J_factor1)
                !          e^(-(J_factor1) (sqrt(0.5 / J_factor1) - 0.5)
                !                          (sqrt(0.5 / J_factor1) + 0.5)
                !        = (2) sqrt(0.5 / J_factor1) (J_factor1)
                !          e^(-(J_factor1) ((0.5 / J_factor1) - 0.25)
                !        = (2) sqrt((0.5) (J_factor1))
                !          e^(-(0.5 - (0.25) (J_factor1))
                !        = (2) sqrt((0.5) (J_factor1)) e^((0.25) (J_factor1) - 0.5)
                !
                prob_Jmp = &                                                            ! -
                        2 * sqrt(J_factor1 / 2) * exp(J_factor1/4 - 0.5d0)

                !The rotational energy of the H2 should be some random value
                !that follows the boltzmann distribution at this temperature
                do
                        !First pick a random rotational quantum number (J)
                        !between zero and some upper limit (J_max)
                        J = rand() * J_max

                        !And some random number between 0 and the highest
                        !probability (prob_Jmp) for metropolis rejection
                        random_num1 = rand() * prob_Jmp

                        !If the probability that this J is less than that
                        !random rejection rate, then it is rejected
                        if ((2*J + 1.0d0) * (J_factor1) * &
                            exp(-(J_factor1) * (J) * (J + 1.0d0)) < random_num1) cycle
                        exit
                end do

                !The rotational energy is computed using the rotational number (J)
                !McQuarrie, 1997 : section 18-5 : pg.743
                !
                ! Erot(J) = (hbar^2) (J) (J + 1) / (2) (moment_inertia)
                !
                ! Erot(J) = (hbar^2) (J) (J + 1) / (2) (moment_inertia)
                !         = ((hbar^2) / 2) (J) (J + 1) / moment_inertia
                !         = (epsilon_factor2) (J) (J + 1) / moment_inertia
                !
                initial_rotational_energy = &                                          !RU_energy
                        epsilon_factor2 * J * (J + 1.0d0) / moment_inertia

                !The linear rotational speed is computed from the rotational energy
                !McQuarrie, 1997 : section 5-8 : pg.174
                !
                !   Erot = (1/2) (moment_inertia) (angular_rotational_speed)^2
                !
                !   Erot = (1/2) (moment_inertia) (angular_rotational_speed)^2
                !        = (1/2) (moment_inertia) (linear_rotational_speed  / radius)^2
                !        = (1/2) (moment_inertia) (linear_rotational_speed / (bond_length / 2))^2
                !        = (2) (moment_inertia) (linear_rotational_speed)^2 / bond_length^2
                !
                ! linear_rotational_speed = sqrt((1/2) (E_rot) (initial_bond_distance)^2 / moment_inertia)
                ! linear_rotational_speed = sqrt(E_rot / mass_hydrogen)
                !
                initial_rotational_speed = &                                          !RU_length/RU_time
                        sqrt(initial_rotational_energy/mass_hydrogen)

                !The bond is rotating, meaning some part of the two atoms' velocity
                !vectors are perpendicular to the bond. Which direction it is going
                !should be random
                !
                !Two degrees of freedom can be decribed by a single angle
                !
                random_num1 = rand()
                initial_rotation_angle = random_num1*pi2 - pi                         !Radians

                !The bond in actuality may not be completely stretched out or
                !contracted so a random point in the period is picked to
                !start the bond at
                !
                !This is normalized to a value between 0 and 1
                !
                bond_period_elapsed = rand()                                          ! -

                !All of this is stored for later use in the InitialSetup of runTrajectory
                INITIAL_BOND_DATA(:,m) = (/ initial_bond_distance,initial_rotational_speed,&
                                        initial_rotation_angle,initial_bond_angle1,initial_bond_angle2,&
                                        bond_period_elapsed /)
        end do

end subroutine InitialSampling3






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               NonBondedForce
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!       PURPOSE
!               This subroutine gets the force between two nonbonded atoms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               coord1                          REAL(DP),DIM(3)                 The coordinates of the first atom
!               coord2                          REAL(DP),DIM(3)                 The coordinates of the second atom
!               r                               REAL(DP),OPTIONAL               The distance between the two atoms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!               gradient1                       REAL(DP),DIM(3)                 The force on the first atom
!               gradient2                       REAL(DP),DIM(3)                 The force on the second atom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               distance12                      REAL(DP)                        The distance between the two atoms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               BondedForce
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!       PURPOSE
!               This subroutine gets the force between two bonded atoms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               coord1                          REAL(DP),DIM(3)                 The coordinates of the first atom
!               coord2                          REAL(DP),DIM(3)                 The coordinates of the second atom
!               r                               REAL(DP),OPTIONAL               The distance between the two atoms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!               gradient1                       REAL(DP),DIM(3)                 The force on the first atom
!               gradient2                       REAL(DP),DIM(3)                 The force on the second atom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               distance12                      REAL(DP)                        The distance between the two atoms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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



