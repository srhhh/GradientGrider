module f2_physics_parameters
implicit none

!PHYSICAL CONSTANTS (all in SI)
real,parameter :: Na = 6.02214e23
real,parameter :: kb = 1.38065e-23
real,parameter :: Da = 1.66054e-27
real,parameter :: c  = 2.99792e8
real,parameter :: pi = 3.14159
real,parameter :: pi2 = pi*2.0
real,parameter :: eV = 1.60218e-19 
real,parameter :: hbar = 1.05457e-34

!REDUCED UNITS
real,parameter :: RU_length = 1.0e-10
real,parameter :: RU_mass = Da
real,parameter :: RU_time = 1.0e-15
real,parameter :: RU_force = RU_mass*RU_length*RU_time**(-2)
real,parameter :: RU_energy = RU_mass*(RU_length/RU_time)**2

!TIME PARAMETERS
real,parameter :: dt = (1.0e-17)/RU_time
integer,parameter :: Nsteps = 30000

!ATOMIC PARAMETERS
! RS: why naming this "HOr0_hydrogen"? Do you mean "HHr0"?
!		KF: HO stands for harmonic oscillator
real,parameter :: HOr0_hydrogen = (1.0e-10)*0.7412/RU_length
                                        ! originally in A
! RS: why naming this "HOke_hydrogen"? Do you mean "HHke"?
real,parameter :: HOke_hydrogen =  (1.0e10*1.0e-8)*5.756/(RU_force/RU_length)
                                        ! originally in mdyn/A= 1.0e-8 N/A
real,parameter :: Morser0_hydrogen = (1.0e-10)*3.25/RU_length
                                        ! originally in A
real,parameter :: MorseDe_hydrogen = (1.0/(Na*2.39006e-4))*0.02347/RU_energy
                                        ! originally in kcal/mol
real,parameter :: Morsealpha_hydrogen = -10.6
                                        ! unitless
real,parameter :: mass_hydrogen = (.001/Na)*(1.00794)/RU_mass
                                        !originally g/mol
!CUTOFF PARAMETERS
! RS: This is an interesting way of defining it -- why don't just use things like 9*Morser0_hydrogen?
! RS: and don't do the sqrt operation in the main program?
!			KF: I no longer use a cutoff for Morse potential so I deleted this

!TEMPERATURE PARAMETERS
real, parameter :: temperature = 200.0
					!Kelvin
real, parameter :: upsilon_max = 5.0
real, parameter :: frequency_factor  = sqrt(HOke_hydrogen/(0.5*mass_hydrogen))/pi2
real, parameter :: upsilon_factor1 = -pi2*(hbar/(RU_energy*RU_time))*frequency_factor/(kb*temperature/RU_energy)
real, parameter :: upsilon_factor2 = pi2*(hbar/(RU_energy*RU_time))*frequency_factor
real, parameter :: temperature_factor = exp(0.5*upsilon_factor1)
!real, parameter :: temperature_scaling = (1.0-temperature_factor)*Emax/
real, parameter :: temperature_scaling = (1.0-temperature_factor)*upsilon_max/temperature_factor

!COLLISION PARAMETERS
real,parameter :: initial_translational_KE = (1.0)*eV/RU_energy
!real,parameter :: collision_distance = Morser0_hydrogen*5.0
 real,parameter :: collision_distance = 9.5
real,parameter :: collision_skew = HOr0_hydrogen*0.0


end module f2_physics_parameters



