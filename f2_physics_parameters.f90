module f2_physics_parameters
implicit none

!PHYSICAL CONSTANTS
real,parameter :: Na = 6.02214e23
real,parameter :: kb = 1.38065e-23
real,parameter :: Da = 1.66054e-27
real,parameter :: c  = 2.99792e8
real,parameter :: pi = 3.14159
real,parameter :: eV = 1.60218e-19 

!REDUCED UNITS
real,parameter :: RU_length = 1.0e-10
real,parameter :: RU_mass = Da
real,parameter :: RU_time = 1.0e-15
real,parameter :: RU_force = RU_mass*RU_length*RU_time**(-2)
real,parameter :: RU_energy = RU_mass*(RU_length/RU_time)**2

!TIME PARAMETERS
real,parameter :: dt = (1.0e-17)/RU_time
integer,parameter :: Nsteps = 20000

!ATOMIC PARAMETERS
real,parameter :: HOr0_hydrogen = (1.0e-10)*0.7412/RU_length
                                        !originally A
real,parameter :: HOke_hydrogen =  (1.0e10*1.0e-8)*5.756/(RU_force/RU_length)
                                        !originally mdyn/A= 1.0e-8 N/A
real,parameter :: Morser0_hydrogen = (1.0e-10)*3.25/RU_length
                                        !originally A
real,parameter :: MorseDe_hydrogen = (1.0/(Na*2.39006e-4))*0.02347/RU_energy
                                        !originally kcal/mol
real,parameter :: Morsealpha_hydrogen = -10.6
                                        !unitless
real,parameter :: mass_hydrogen = (.001/Na)*(1.00794)/RU_mass
                                        !originally g/mol
!CUTOFF PARAMETERS
real, parameter :: cutoff_distance_squared = 9*Morser0_hydrogen**2

!TEMPERATURE PARAMETERS
real, parameter :: temperature = 200.0
					!Kelvin
real, parameter :: temperature_constant = kb*temperature/RU_energy
real, parameter :: temperature_scaling = sqrt((2*pi*temperature_constant/mass_hydrogen)**3)
real, parameter :: average_KE = 1.5*temperature_constant

end module f2_physics_parameters


