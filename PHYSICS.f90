module PHYSICS
implicit none

integer,parameter :: dp = kind(1.0d0)

!PHYSICAL CONSTANTS (all in SI)
real(dp),parameter :: Na = 6.02214d23
real(dp),parameter :: kb = 1.38065d-23
real(dp),parameter :: Da = 1.66054d-27
real(dp),parameter :: c  = 2.99792e8
real(dp),parameter :: pi = 3.14159d0
real(dp),parameter :: pi2 = pi*2.0d0
real(dp),parameter :: eV = 1.60218d-19 
real(dp),parameter :: hbar = 1.05457d-34

!REDUCED UNITS
real(dp),parameter :: RU_length = 1.0d-10
real(dp),parameter :: RU_mass = Da
real(dp),parameter :: RU_time = 1.0d-15
real(dp),parameter :: RU_force = RU_mass*RU_length*RU_time**(-2)
real(dp),parameter :: RU_energy = RU_mass*(RU_length/RU_time)**2

!TIME PARAMETERS
real(dp),parameter :: dt = (1.0d-17)/RU_time
integer,parameter :: Nsteps = 30000

!ATOMIC PARAMETERS
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

!TEMPERATURE PARAMETERS
real(dp), parameter :: temperature = 200.0d0
					!Kelvin
real(dp), parameter :: upsilon_max = 5.0d0
real(dp), parameter :: frequency_factor  = sqrt(HOke_hydrogen/(0.5d0*mass_hydrogen))/pi2
real(dp), parameter :: upsilon_factor1 = -pi2*(hbar/(RU_energy*RU_time))*frequency_factor/(kb*temperature/RU_energy)
real(dp), parameter :: upsilon_factor2 = pi2*(hbar/(RU_energy*RU_time))*frequency_factor
real(dp), parameter :: temperature_factor = exp(0.5d0*upsilon_factor1)
real(dp), parameter :: temperature_scaling = (1.0d0-temperature_factor)*upsilon_max/temperature_factor

!COLLISION PARAMETERS
real(dp),parameter :: initial_translational_KE = (1.0d0)*eV/RU_energy
!real(dp),parameter :: collision_distance = Morser0_hydrogen*5.0
 real(dp),parameter :: collision_distance = 9.5d0
real(dp),parameter :: collision_skew = HOr0_hydrogen*0.0d0

contains
subroutine decompose_two_velocities(coords,velocities,&
           velocity_translation,velocity_vibration,velocity_rotation)

implicit none
real(dp),dimension(3,2),intent(in) :: coords,velocities
real(dp),dimension(3),intent(out) :: velocity_translation,velocity_vibration,velocity_rotation
real(dp),dimension(3) :: parallel_vector,orthogonal_vector
real(dp),dimension(3) :: velocity_translation1,velocity_translation2,velocity_translation3
real(dp),dimension(3) :: velocity_rotation1,velocity_rotation2
real(dp),dimension(3,2) :: velocities_parallel
real(dp),dimension(3,2) :: velocities_orthogonal
real(dp),dimension(3,2) :: velocities_ortho2
real(dp),dimension(3,2) :: velocities_ortho1

parallel_vector = coords(:,1) - coords(:,2)
parallel_vector = parallel_vector / sqrt(parallel_vector(1)**2 +&
                                         parallel_vector(2)**2 +&
                                         parallel_vector(3)**2 )

velocities_parallel(:,1) = dot_product(velocities(:,1),&
                                      parallel_vector) * parallel_vector
velocities_parallel(:,2) = dot_product(velocities(:,2),&
                                       parallel_vector) * parallel_vector

velocities_orthogonal = velocities - velocities_parallel

orthogonal_vector = velocities_orthogonal(:,1) + velocities_orthogonal(:,2)
orthogonal_vector = orthogonal_vector / sqrt(orthogonal_vector(1)**2 +&
                                             orthogonal_vector(2)**2 +&
                                             orthogonal_vector(3)**2 )

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


end module PHYSICS



