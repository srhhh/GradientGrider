!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       MODULE
!               runTrajectory
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!		This module simulates an MD collision of two molecules
!		as described in PHYSICS
!
!		A frame advances in time to another frame with the
!		verlet method
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILECHANNELS                    ACTION
!
!               FILECHANNEL1                    OPEN, WRITE, CLOSE
!               FILECHANNEL2                    OPEN, CALL CHECKSTATE, WRITE, CLOSE
!               FILECHANNELS                    none
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	SUBROUTINES			ARGUMENTS		KIND
!
!		addTrajectory			coords_initial		intent(out),real(dp),dim(3,Natoms)
!						velocities_initial	intent(out),real(dp),dim(3,Natoms)
!						coords_final		intent(out),real(dp),dim(3,Natoms)
!						velocities_final	intent(out),real(dp),dim(3,Natoms)
!
!		checkaddTrajectory		coords_initial		intent(out),real(dp),dim(3,Natoms)
!						velocities_initial	intent(out),real(dp),dim(3,Natoms)
!						coords_final		intent(out),real(dp),dim(3,Natoms)
!						velocities_final	intent(out),real(dp),dim(3,Natoms)
!
!		checkTrajectory			coords_initial		intent(out),real(dp),dim(3,Natoms)
!						velocities_initial	intent(out),real(dp),dim(3,Natoms)
!						coords_final		intent(out),real(dp),dim(3,Natoms)
!						velocities_final	intent(out),real(dp),dim(3,Natoms)
!
!		checkMultipleTrajectories	filechannels		intent(in),dim(Ngrid_total+1),integer
!						coords_initial		intent(out),real(dp),dim(3,Natoms)
!						velocities_initial	intent(out),real(dp),dim(3,Natoms)
!						coords_final		intent(out),real(dp),dim(3,Natoms)
!						velocities_final	intent(out),real(dp),dim(3,Natoms)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       CALLS                           MODULE
!
!		getVarMaxMin			VARIABLES
!
!		InitialSetup3			PHYSICS
!
!		Acceleration			runTrajectory
!
!		addState			interactSingleGrid
!		checkState			interactSingleGrid
!		checkState			interactMultipleGrids
!
!		BondedForce			PHYSICS
!		NonBondedForce			PHYSICS
!
!		MorsePotential			PHYSICS
!		HOPotential			PHYSICS
!		KineticEnergy			PHYSICS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES	                          FILETYPE	SUBROUTINE				FMT
!
!		gridpath0//trajectory		XYZ		addTrajectory			'(A1,3F10.6)'
!		gridpath0//trajectory		XYZ		checkMultipleTrajectories	'(A1,3F10.6)'
!		gridpath1//trajectory		XYZ		checkTrajectory			'(A1,3F10.6)'
!		gridpath1//checkstatefile	DAT		checkTrajectory			*
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



module runTrajectory
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	SUBROUTINE
!		addTrajectory
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	PURPOSE
!		This subroutine simulates an MD collision by calculating the gradient
!		at every step with classical mechanics and parameters in PHYSICS
!
!		Every frame is added to the current grid through interactSingleGrid
!
!		The initial and final velocity of the incoming H/H2 is output to
!		gather scattering angle data
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	INPUT				KIND				DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	OUTPUT				KIND				DESCRIPTION
!
!               coords_initial                  REAL(DP),DIM(3,Natoms)          The coordinates of the initial frame
!               velocities_initial              REAL(DP),DIM(3,Natoms)          The velocities of the initial frame
!               coords_final                    REAL(DP),DIM(3,Natoms)          The coordinates of the final frame
!               velocities_final                REAL(DP),DIM(3,Natoms)          The velocities of the final frame
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	IMPORTANT VARIABLES		KIND				DESCRIPTION
!
!		vals				REAL(DP),DIM(Nvar)		The variables associated with a frame
!
!		coords				REAL(DP),DIM(3,Natoms)		The coordinates defining a frame
!		velocities			REAL(DP),DIM(3,Natoms)		The velocities of a frame
!		gradient			REAL(DP),DIM(3,Natoms)		The gradient associated with a frame
!
!		Norder1				INTEGER				The number of times a frame enters an order1 cell
!		steps				INTEGER				The number of timesteps advanced
!		header_max_flag			LOGICAL				The flag indicating whether counter1 is full or not
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES	                          FILETYPE			DESCRIPTION
!
!		gridpath0//trajectory		XYZ				Coordinates to view the trajectory over time;
!										currently just for bug-testing (set to false)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!subroutine addTrajectory(coords_initial,velocities_initial,coords_final,velocities_final)
!	use VARIABLES
!	use PHYSICS
!	use PARAMETERS
!	use interactSingleGrid
!        implicit none
!
!	!Coordinates, Velocities, and Variables
!	real(dp), dimension(3,Natoms) :: coords,gradient,velocities
!	real(dp), dimension(Nvar) :: vals
!	real(dp),dimension(3,Natoms),intent(out) :: coords_initial,velocities_initial
!	real(dp),dimension(3,Natoms),intent(out) :: coords_final,velocities_final
!	integer :: bond_index1, bond_index2
!
!	!Incremental Integer
!	integer :: n
!
!        !Initialize the scene
!        call InitialSetup3(coords,velocities)
!	Norder1 = 0
!
!	coords_initial = coords
!	velocities_initial = velocities
!
!	!Always calculate the variables before accelerating
!	!because we can reuse these calculations
!	call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
!
!        !Accelerate the velcocities for a half step (verlet)
!        call Acceleration(vals,coords,gradient)
!	call addState(vals,coords,gradient)
!
!	!Update the velocities
!	velocities = velocities + 0.5d0 * gradient
!
!	!To randomize the periods of the bond, I let the scene go on
!	!for a small period of time (need to standardize this later)
!	do n = 1, Nbonds
!		do steps = 1, int(INITIAL_BOND_DATA(6,n)*vib_period)
!			coords = coords + dt * velocities
!			call Acceleration(vals,coords,gradient)
!			velocities = velocities + gradient
!		end do
!
!		!And then reset the bond
!		coords(:,BONDING_DATA(n,1)) = coords_initial(:,BONDING_DATA(n,1))
!		coords(:,BONDING_DATA(n,2)) = coords_initial(:,BONDING_DATA(n,2))
!		velocities(:,BONDING_DATA(n,1)) = velocities_initial(:,BONDING_DATA(n,1))
!		velocities(:,BONDING_DATA(n,2)) = velocities_initial(:,BONDING_DATA(n,2))
!	end do
!
!	!Now we go into the mainloop
!	!We have a hard cap of Nsteps timesteps
!        do steps = 1, Nsteps
!
!		!Just for bug-testing
!                if (.false.) then !(modulo(steps,50) == 0) then
!                        open(filechannel1,file=gridpath0//trajectoryfile,position="append")
!                        write(filechannel1,'(I6)') Natoms
!                        write(filechannel1,*) ""
!			do n = 1, Natoms
!                        write(filechannel1,'(A1,3F10.6)') 'H',&
!                                coords(1,n), coords(2,n), coords(3,n)
!			end do
!			close(filechannel1)
!                end if
! 
!                !Check every 500 steps if we are out-of-bounds
!                if (modulo(steps,500) == 1) then      
! 			if ((vals(1)>max_var1) .or. (vals(2)>max_var2)) then
! 				exit
! 			end if
!                endif
!
!                !Update the coordinates with the velocities
!		coords = coords + dt * velocities
!
!		!Always calculate the variables before adding a frame or accelerating
!		call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
!
!                !Accelerate and update gradients
!                call Acceleration(vals,coords,gradient)
!
!		!Add the frame to the grid
!        	call addState(vals,coords,gradient)
!
!		!If there are too many subdivisions and counter0 will get out-of-bounds
!		!we would have to call this to exit
!                if (header_max_flag) exit
!
!		!Update the velocities
!		velocities = velocities + gradient
!        end do
!
!	coords_final = coords
!	velocities_final = velocities
!
!end subroutine addTrajectory
!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	SUBROUTINE
!		checkaddTrajectory
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	PURPOSE
!		This subroutine simulates an MD collision by calculating the gradient
!		at every step with classical mechanics and parameters in PHYSICS
!
!		Every frame is checked with the current grid through interactSingleGrid
!
!		The initial and final velocity of the incoming H/H2 is output to
!		gather scattering angle data
!
!		Data is output every frame to gather information on how the trajectory
!		interacts with the grid and whether it is accurate, realistic, etc.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	INPUT				KIND				DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	OUTPUT				KIND				DESCRIPTION
!
!               coords_initial                  REAL(DP),DIM(3,Natoms)          The coordinates of the initial frame
!               velocities_initial              REAL(DP),DIM(3,Natoms)          The velocities of the initial frame
!               coords_final                    REAL(DP),DIM(3,Natoms)          The coordinates of the final frame
!               velocities_final                REAL(DP),DIM(3,Natoms)          The velocities of the final frame
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	IMPORTANT VARIABLES		KIND				DESCRIPTION
!
!		vals				REAL(DP),DIM(Nvar)		The variables associated with a frame
!
!		coords				REAL(DP),DIM(3,Natoms)		The coordinates defining a frame
!		velocities			REAL(DP),DIM(3,Natoms)		The velocities of a frame
!		gradient			REAL(DP),DIM(3,Natoms)		The gradient associated with a frame
!		approx_gradient			REAL(DP),DIM(3,Natoms)		The gradient recovered for a frame from the grid
!
!		number_of_frames		INTEGER				The total amount of frames checked for a frame
!		order				INTEGER				The order of the cell checked for a frame
!		neighbor_check			INTEGER				Signifies whether neighbor cells were check
!										or not for a frame(0 or 1)
!
!		min_rmsd			REAL(DP)			The minimum rmsd recovered for a frame from the grid
!		min_rmsd_prime			REAL(DP)			The minimum rmsd recovered for a frame from the grid
!										after a forced neighbor check (may not apply)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES	                          FILETYPE			DESCRIPTION
!
!		gridpath1//trajectory		XYZ				Coordinates to view the trajectory over time;
!										currently it is written every 10 timesteps
!		gridpath1//checkstatefile	DAT				The DAT file containing information about the
!										trajectory over time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine checkaddTrajectory(filechannels,&
                              coords_initial,velocities_initial,&
                              coords_final,velocities_final)
        use PARAMETERS
        use PHYSICS
        use VARIABLES
        use ANALYSIS
        use interactMultipleGrids
        implicit none

        !Coordinates, Velocities, and Variables
        real(dp), dimension(3,Natoms) :: coords,coords_labelled,velocities
        real(dp), dimension(3,Natoms) :: gradient, approx_gradient
        real(dp), dimension(Nvar) :: vals
	real(dp), dimension(3,Natoms), intent(out) :: coords_initial, velocities_initial 
	real(dp), dimension(3,Natoms), intent(out) :: coords_final, velocities_final 
        integer,dimension(1+Ngrid_max) :: filechannels
        integer :: bond_index1, bond_index2

        !Various other variables
        real(dp) :: U, KE
        real(dp) :: min_rmsd,min_rmsd_prime
        integer :: number_of_frames,order,neighbor_check

        !Incremental Integer
        integer :: n

        !Initialize the scene
        call InitialSetup3(coords,velocities)
        Norder1 = 0
        Norder_total = 0

        coords_initial = coords
        velocities_initial = velocities

        !Always calculate the variables before accelerating
        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

        !Accelerate the velcocities for a half step (verlet)
        call Acceleration(vals,coords,gradient)

        !Update the velocities
        velocities = velocities + 0.5d0 * gradient

        !To randomize the periods of the bond, I let the scene go on
        !for a small period of time (need to standardize this later)
        do n = 1, Nbonds
                do steps = 1, int(INITIAL_BOND_DATA(6,n)*vib_period)
                        coords = coords + dt * velocities
                        call Acceleration(vals,coords,gradient)
                        velocities = velocities + gradient
                end do

                !And then reset the bond
                coords(:,BONDING_DATA(n,1)) = coords_initial(:,BONDING_DATA(n,1))
                coords(:,BONDING_DATA(n,2)) = coords_initial(:,BONDING_DATA(n,2))
                velocities(:,BONDING_DATA(n,1)) = velocities_initial(:,BONDING_DATA(n,1))
                velocities(:,BONDING_DATA(n,2)) = velocities_initial(:,BONDING_DATA(n,2))
        end do

        !We keep this file open for the whole trajectory (instead of
        !continually opening and closing) to keep data of each frame
if (testtrajDetailedRMSD_flag) open(filechannel2,file=gridpath1//checkstatefile)

        buffer1_size = 1 + var_overcrowd(1)
        allocate(valsbuffer1(Nvar,buffer1_size),&
                 coordsbuffer1(3,Natoms,buffer1_size),&
                 gradientbuffer1(3,Natoms,buffer1_size),&
                 Ubuffer1(3,3,buffer1_size),&
                 RMSDbuffer1(buffer1_size),&
                 approximation_index(buffer1_size))

        do steps = 1, Nsteps

                !Check every 500 steps to see if we are out-of-bounds
                if (modulo(steps,500) == 1) then
                        if (any(vals > var_maxvar)) exit
                endif

                !Upate the coordinates with the velocities
                coords = coords + dt * velocities

                !Always calculate the variables before checking a frame or accelerating
                call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

                !If the library forced no special labeling or duplicate additions, then
                !that means we don't do any special label switching when using the
                !approximated gradient and inputing the frame
                if ((force_Duplicates) .or. (force_NoLabels)) then
                        min_rmsd = default_rmsd
                        subcellsearch_max = (/ 0, 0 /)
                        call checkState_new(vals,coords,approx_gradient,min_rmsd,&
                                 filechannels,number_of_frames,order,neighbor_check)

                        if ((min_rmsd .ge. threshold_RMSD).or.(reject_flag)&
                            .or.(accept_worst)) then
                                call Acceleration(vals,coords,gradient)
                                if (grid_addition) call addState_new(vals,coords,gradient)
                        else
                                gradient = approx_gradient
                        end if

                !Otherwise, then we need to do the label switching
                else
                        !We need to input the frame using the labeling scheme
                        do n = 1, Natoms
                                coords_labelled(:,n) = coords(:,BOND_LABELLING_DATA(n))
                        end do

                       !Check for a frame in the grid
                       !Set the default value beforehand though
                       min_rmsd = default_rmsd
                       subcellsearch_max = subcellsearch_max1
                       call checkState_new(vals,coords_labelled,approx_gradient,min_rmsd,&
                                filechannels,number_of_frames,order,neighbor_check)

if (testtrajDetailedRMSD_flag) then

        min_rmsd_prime = default_rmsd

        subcellsearch_max = subcellsearch_max2
        call checkState_new(vals,coords_labelled,approx_gradient,min_rmsd_prime,&
                filechannels,number_of_frames,order,neighbor_check)

        !Calculate the potential and kinetic energy
        !(yes, making this nicer is next on the to-do list)
        !Right now, it is not generic (we only take into account 3 atoms)
        U = MorsePotential(coords(:,1),coords(:,2))
        U = U + MorsePotential(coords(:,1),coords(:,3))
        U = U + HOPotential(coords(:,2),coords(:,3))
        KE = KineticEnergy(velocities(:,1))
        KE = KE + KineticEnergy(velocities(:,2))
        KE = KE + KineticEnergy(velocities(:,3))

        !Finally write to the data file all the important data values
        write(filechannel2,FMT=*) number_of_frames,order,neighbor_check,steps,&
                                  min_rmsd,min_rmsd_prime,vals(1),vals(2),U,KE

        min_rmsd = min_rmsd_prime
end if

                       !Update the gradient with either the approximation or by acclerating
                       !This is dependent on the threshold and the rejection method
                       if ((min_rmsd .ge. threshold_RMSD).or.(reject_flag)&
                            .or.(accept_worst)) then
                                call Acceleration(vals,coords,gradient)
                                if (grid_addition) call addState_new(vals,coords,gradient)
                        else
                                !And we need to use the approximation after "unlabeling"
                                do n = 1, Natoms
                                         gradient(:,BOND_LABELLING_DATA(n)) = approx_gradient(:,n)
                                end do
                       end if
                end if

                !Update the velocities
                velocities = velocities + gradient
        end do

if (testtrajDetailedRMSD_flag) close(filechannel2)

        deallocate(valsbuffer1,coordsbuffer1,gradientbuffer1,Ubuffer1,RMSDbuffer1,&
                   approximation_index)

        coords_final = coords
        velocities_final = velocities

end subroutine checkaddTrajectory




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	SUBROUTINE
!		checkTrajectory
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	PURPOSE
!		This subroutine simulates an MD collision by calculating the gradient
!		at every step with classical mechanics and parameters in PHYSICS
!
!		Every frame is checked with the current grid through interactSingleGrid
!
!		The initial and final velocity of the incoming H/H2 is output to
!		gather scattering angle data
!
!		Data is output every frame to gather information on how the trajectory
!		interacts with the grid and whether it is accurate, realistic, etc.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	INPUT				KIND				DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	OUTPUT				KIND				DESCRIPTION
!
!               coords_initial                  REAL(DP),DIM(3,Natoms)          The coordinates of the initial frame
!               velocities_initial              REAL(DP),DIM(3,Natoms)          The velocities of the initial frame
!               coords_final                    REAL(DP),DIM(3,Natoms)          The coordinates of the final frame
!               velocities_final                REAL(DP),DIM(3,Natoms)          The velocities of the final frame
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	IMPORTANT VARIABLES		KIND				DESCRIPTION
!
!		vals				REAL(DP),DIM(Nvar)		The variables associated with a frame
!
!		coords				REAL(DP),DIM(3,Natoms)		The coordinates defining a frame
!		velocities			REAL(DP),DIM(3,Natoms)		The velocities of a frame
!		gradient			REAL(DP),DIM(3,Natoms)		The gradient associated with a frame
!		approx_gradient			REAL(DP),DIM(3,Natoms)		The gradient recovered for a frame from the grid
!
!		number_of_frames		INTEGER				The total amount of frames checked for a frame
!		order				INTEGER				The order of the cell checked for a frame
!		neighbor_check			INTEGER				Signifies whether neighbor cells were check
!										or not for a frame(0 or 1)
!
!		min_rmsd			REAL(DP)			The minimum rmsd recovered for a frame from the grid
!		min_rmsd_prime			REAL(DP)			The minimum rmsd recovered for a frame from the grid
!										after a forced neighbor check (may not apply)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES	                          FILETYPE			DESCRIPTION
!
!		gridpath1//trajectory		XYZ				Coordinates to view the trajectory over time;
!										currently it is written every 10 timesteps
!		gridpath1//checkstatefile	DAT				The DAT file containing information about the
!										trajectory over time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!subroutine checkTrajectory(coords_initial,velocities_initial,coords_final,velocities_final)
!        use PARAMETERS
!        use PHYSICS
!        use VARIABLES
!        use ANALYSIS
!        use interactSingleGrid
!        implicit none
!
!        !Coordinates, Velocities, and Variables
!        real(dp), dimension(3,Natoms) :: coords,coords_labelled,velocities
!        real(dp), dimension(3,Natoms) :: gradient, approx_gradient
!        real(dp), dimension(Nvar) :: vals
!	real(dp), dimension(3,Natoms), intent(out) :: coords_initial, velocities_initial 
!	real(dp), dimension(3,Natoms), intent(out) :: coords_final, velocities_final 
!	integer :: bond_index1, bond_index2
!
!        !Various other variables
!        real(dp) :: U, KE
!        real(dp) :: min_rmsd,min_rmsd_prime
!        integer :: number_of_frames,order,neighbor_check
!
!	!Incremental Integer
!	integer :: n
!
!        !Initialize the scene
!        call InitialSetup3(coords,velocities)
!
!	coords_initial = coords
!	velocities_initial = velocities
!
!	!Always calculate the variables before accelerating
!	call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
!
!        !Accelerate the velcocities for a half step (verlet)
!        call Acceleration(vals,coords,gradient)
!
!	!Update the velocities
!        velocities = velocities + 0.5d0 * gradient
!
!	!To randomize the periods of the bond, I let the scene go on
!	!for a small period of time (need to standardize this later)
!	do n = 1, Nbonds
!		do steps = 1, int(INITIAL_BOND_DATA(6,n)*vib_period)
!			coords = coords + dt * velocities
!			call Acceleration(vals,coords,gradient)
!			velocities = velocities + gradient
!		end do
!
!		!And then reset the bond
!		coords(:,BONDING_DATA(n,1)) = coords_initial(:,BONDING_DATA(n,1))
!		coords(:,BONDING_DATA(n,2)) = coords_initial(:,BONDING_DATA(n,2))
!		velocities(:,BONDING_DATA(n,1)) = velocities_initial(:,BONDING_DATA(n,1))
!		velocities(:,BONDING_DATA(n,2)) = velocities_initial(:,BONDING_DATA(n,2))
!	end do
!
!	!Get the trajectory file open for trajectory visualization
!        open(filechannel1,file=gridpath1//trajectoryfile)
!        write(filechannel1,'(I6)') Natoms
!        write(filechannel1,*) ""
!	do n = 1, Natoms
!        	write(filechannel1,'(A1,3F10.6)') 'H',&
!              		coords(1,n), coords(2,n), coords(3,n)
!        end do
!        close(filechannel1)
!
!	!We keep this file open for the whole trajectory (instead of
!	!continually opening and closing) to keep data of each frame
!	open(filechannel2,file=gridpath1//checkstatefile)
!        do steps = 1, Nsteps
!
!                !Every 10 frames, print to an xyz file for visualization
!                 if (modulo(steps,10) == 0) then
!                        open(filechannel1,file=gridpath1//trajectoryfile,position="append")
!                        write(filechannel1,'(I6)') Natoms
!                        write(filechannel1,*) ""
!			do n = 1, Natoms
!                        	write(filechannel1,'(A1,3F10.6)') 'H',&
!                              		coords(1,n), coords(2,n), coords(3,n)
!                        end do
!                        close(filechannel1)
!                 end if
!
!                !Check every 500 steps to see if we are out-of-bounds
!                if (modulo(steps,500) == 1) then
!                        if ((vals(1)>max_var1) .or.&
!                            (vals(2)>max_var2)) then
!                                exit
!                        end if
!                endif
!
!                !Upate the coordinates with the velocities
!                coords = coords + dt * velocities
!
!                !Always calculate the variables before checking a frame or accelerating
!		call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
!
!                !If the library forced no special labeling or duplicate additions, then
!                !that means we don't do any special label switching when using the
!                !approximated gradient and inputing the frame
!                if ((force_Duplicates) .or. (force_NoLabels)) then
!                        if (accept_worst) then
!                                min_rmsd = 0.0d0
!                        else
!                                min_rmsd = default_rmsd
!                        end if
!                        call checkState(vals,coords,approx_gradient,min_rmsd,&
!                                 number_of_frames,order,neighbor_check)
!
!                        if ((accept_worst) .and. (min_rmsd == 0.0d0)) then
!                                call Acceleration(vals,coords,gradient)
!                        else if ((min_rmsd .ge. threshold_RMSD).or.(reject_flag)) then
!                                call Acceleration(vals,coords,gradient)
!                        else
!				gradient = approx_gradient
!                        end if
!
!                !Otherwise, then we need to do the label switching
!                else
!                        !We need to input the frame using the labeling scheme
!                        do n = 1, Natoms
!                                coords_labelled(:,n) = coords(:,BOND_LABELLING_DATA(n))
!                        end do
!
!                       !Check for a frame in the grid
!                       !Set the default value beforehand though
!                       if (accept_worst) then
!                               min_rmsd = 0.0d0
!                       else
!                               min_rmsd = default_rmsd
!                       end if
!                       call checkState(vals,coords_labelled,approx_gradient,min_rmsd,&
!                                number_of_frames,order,neighbor_check)
!
!                       !Update the gradient with either the approximation or by acclerating
!                       !This is dependent on the threshold and the rejection method
!                       if ((accept_worst) .and. (min_rmsd == 0.0d0)) then
!                                call Acceleration(vals,coords,gradient)
!                       else if ((min_rmsd .ge. threshold_RMSD).or.(reject_flag)) then
!                                call Acceleration(vals,coords,gradient)
!                       else
!                                !And we need to use the approximation after "unlabeling"
!                                do n = 1, Natoms
!                                        gradient(:,BOND_LABELLING_DATA(n)) = approx_gradient(:,n)
!                                end do
!                       end if
!                end if
!
!		!Calculate the potential and kinetic energy
!		!(yes, making this nicer is next on the to-do list)
!		!Right now, it is not generic (we only take into account 3 atoms)
!	        U = MorsePotential(coords(:,1),coords(:,2))
!	        U = U + MorsePotential(coords(:,1),coords(:,3))
!	        U = U + HOPotential(coords(:,2),coords(:,3))
!	        KE = KineticEnergy(velocities(:,1))
!	        KE = KE + KineticEnergy(velocities(:,2))
!	        KE = KE + KineticEnergy(velocities(:,3))
!
!		!Finally write to the data file all the important data values
!		write(filechannel2,FMT=*) number_of_frames,order,neighbor_check,steps,&
!                                          min_rmsd,min_rmsd_prime,vals(1),vals(2),U,KE
!
!                !Update the velocities
!                velocities = velocities + gradient
!
!        end do
!	close(filechannel2)
!
!	coords_final = coords
!	velocities_final = velocities
!
!end subroutine checkTrajectory
!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	SUBROUTINE
!		checkMultipleTrajectories
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	PURPOSE
!		This subroutine simulates an MD collision by calculating the gradient
!		at every step with classical mechanics and parameters in PHYSICS
!
!		Every frame is checked with multiple grids through interactMultipleGrids
!
!		The initial and final velocity of the incoming H/H2 is output to
!		gather scattering angle data
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	INPUT				KIND				DESCRIPTION
!
!		filechannels			INTEGER,DIM(Ngrid_total+1)	These filechannels were opened beforehand to the
!										trajectory files that we output rmsd data to;
!										all that's left to do is to write to them every frame.
!                                                                               However, filechannel(1) is NOT opened beforehand
!                                                                               and is available for any random purpose
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	OUTPUT				KIND				DESCRIPTION
!
!               coords_initial                  REAL(DP),DIM(3,Natoms)          The coordinates of the initial frame
!               velocities_initial              REAL(DP),DIM(3,Natoms)          The velocities of the initial frame
!               coords_final                    REAL(DP),DIM(3,Natoms)          The coordinates of the final frame
!               velocities_final                REAL(DP),DIM(3,Natoms)          The velocities of the final frame
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	IMPORTANT VARIABLES		KIND				DESCRIPTION
!
!		vals				REAL(DP),DIM(Nvar)		The variables associated with a frame
!
!		coords				REAL(DP),DIM(3,Natoms)		The coordinates defining a frame
!		velocities			REAL(DP),DIM(3,Natoms)		The velocities of a frame
!		gradient			REAL(DP),DIM(3,Natoms)		The gradient associated with a frame
!		approx_gradient			REAL(DP),DIM(3,Natoms)		The gradient recovered for a frame from the grid
!
!		number_of_frames		INTEGER				The total amount of frames checked over all cells
!		order				INTEGER				The order of the cell checked for a frame
!		neighbor_check			INTEGER				Signifies whether neighbor cells were check or not (0 or 1)
!
!		min_rmsd			REAL(DP)			The minimum rmsd recovered for a frame from the grid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES	                          FILETYPE			DESCRIPTION
!
!		gridpath0//trajectory		XYZ				Coordinates to view the trajectory over time;
!										currently just for bug-testing (set to false)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine checkMultipleTrajectories(filechannels,&
                                                  coords_initial,velocities_initial,&
                                                  coords_final,velocities_final)
        use PHYSICS
        use PARAMETERS
        use VARIABLES
        use ANALYSIS
        use interactMultipleGrids
        implicit none

        !Coordinates, Velocities, and Variables
        real(dp),dimension(3,Natoms) :: coords,coords_labelled,velocities
        real(dp),dimension(3,Natoms) :: gradient,approx_gradient
        real(dp),dimension(3,Natoms),intent(out) :: coords_initial, velocities_initial
        real(dp),dimension(3,Natoms),intent(out) :: coords_final, velocities_final
        real(dp),dimension(Nvar) :: vals

        !Grid Parameters
        integer,dimension(1+Ngrid_total),intent(in) :: filechannels
        integer :: OMP_GET_THREAD_NUM

        !Various other variables
        real(dp) :: U, KE
        real(dp) :: min_rmsd, min_rmsd_prime
        integer :: number_of_frames,order,neighbor_check
        integer :: order0, order1

        !Incremental Integer
        integer :: n

        !Initialize the scene
        call InitialSetup3(coords,velocities)

        !For traversal
        if (traversal_flag) then
                traversal0 = 0
                traversal1 = 0
        end if

        if (testtrajDetailedRMSD_flag) open(filechannel2,file=gridpath0//checkstatefile)
        buffer1_size = 1 + var_overcrowd(1)
        allocate(valsbuffer1(Nvar,buffer1_size),&
                 coordsbuffer1(3,Natoms,buffer1_size),&
                 gradientbuffer1(3,Natoms,buffer1_size),&
                 Ubuffer1(3,3,buffer1_size),&
                 RMSDbuffer1(buffer1_size),&
                 approximation_index(buffer1_size))

        coords_initial = coords
        velocities_initial = velocities

        !Always calculate the variables before accelearting
        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

        !Accelerate the velcocities for a half step (verlet)
        call Acceleration(vals,coords,gradient)

        !Update the velocities
        velocities = velocities + 0.5d0 * gradient

        !To randomize the periods of the bond, I let the scene go on
        !for a small period of time (need to standardize this later)
        do n = 1, Nbonds
                do steps = 1, int(INITIAL_BOND_DATA(6,n)*vib_period)
                        coords = coords + dt * velocities
                        call Acceleration(vals,coords,gradient)
                        velocities = velocities + gradient
                end do

                !And then reset the bond
                coords(:,BONDING_DATA(n,1)) = coords_initial(:,BONDING_DATA(n,1))
                coords(:,BONDING_DATA(n,2)) = coords_initial(:,BONDING_DATA(n,2))
                velocities(:,BONDING_DATA(n,1)) = velocities_initial(:,BONDING_DATA(n,1))
                velocities(:,BONDING_DATA(n,2)) = velocities_initial(:,BONDING_DATA(n,2))
        end do

        !Start the main loop
        do steps = 1, Nsteps

                !Just for bug-testing
                if (.false.) then !(modulo(steps,50) == 1) then
                        open(filechannel1,file=gridpath0//trajectoryfile,position="append")
                        write(filechannel1,'(I6)') Natoms
                        write(filechannel1,*) ""
                        do n = 1, Natoms
                                write(filechannel1,'(A1,3F10.6)') 'H',&
                                              coords(1,n), coords(2,n), coords(3,n)
                        end do
                        close(filechannel1)
                end if

                !Check every 500 steps if we are out-of-bounds
                if (modulo(steps,500) == 1) then
                        if (any(vals > var_maxvar)) exit
                endif

                !Upate the coordinates with the velocities
                coords = coords + dt * velocities

                !Always calculate the variables before checking a frame or accelerating
                call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

                !If the library forced no special labeling or duplicate additions, then
                !that means we don't do any special label switching when using the
                !approximated gradient and inputing the frame
                if ((force_Duplicates) .or. (force_NoLabels)) then
                        if (accept_worst) then
                                min_rmsd = 0.0d0
                        else
                                min_rmsd = default_rmsd
                        end if
!                        call checkState(vals,coords,approx_gradient,min_rmsd,&
!                                 filechannels,number_of_frames,order,neighbor_check)
                        call checkState_new(vals,coords,approx_gradient,min_rmsd,&
                                 filechannels,number_of_frames,order,neighbor_check)

                        if ((accept_worst) .and. (min_rmsd == 0.0d0)) then
                                call Acceleration(vals,coords,gradient)
                        else if ((min_rmsd .ge. threshold_RMSD).or.(reject_flag)) then
                                call Acceleration(vals,coords,gradient)
                        else
                                gradient = approx_gradient
                        end if

                !Otherwise, then we need to do the label switching
                else
                        !We need to input the frame using the labeling scheme
                        do n = 1, Natoms
                                coords_labelled(:,n) = coords(:,BOND_LABELLING_DATA(n))
                        end do

                       !Check for a frame in the grid
                       !Set the default value beforehand though
                       if (accept_worst) then
                               min_rmsd = 0.0d0
                       else
                               min_rmsd = default_rmsd
                       end if
                       subcellsearch_max = subcellsearch_max1
                       call checkState_new(vals,coords_labelled,approx_gradient,min_rmsd,&
                                filechannels,number_of_frames,order0,neighbor_check)

if (testtrajDetailedRMSD_flag) then

        min_rmsd_prime = default_rmsd

        subcellsearch_max = subcellsearch_max2
        call checkState_new(vals,coords_labelled,approx_gradient,min_rmsd_prime,&
                filechannels,number_of_frames,order1,neighbor_check)

        !Calculate the potential and kinetic energy
        !(yes, making this nicer is next on the to-do list)
        !Right now, it is not generic (we only take into account 3 atoms)
        U = MorsePotential(coords(:,1),coords(:,2))
        U = U + MorsePotential(coords(:,1),coords(:,3))
        U = U + HOPotential(coords(:,2),coords(:,3))
        KE = KineticEnergy(velocities(:,1))
        KE = KE + KineticEnergy(velocities(:,2))
        KE = KE + KineticEnergy(velocities(:,3))

        !Finally write to the data file all the important data values
        write(filechannel2,FMT=*) number_of_frames,order0+order1,neighbor_check,steps,&
                                  min_rmsd,min_rmsd_prime,vals(1),vals(2),U,KE

        min_rmsd = min_rmsd_prime
end if

                       !Update the gradient with either the approximation or by acclerating
                       !This is dependent on the threshold and the rejection method
                       if ((accept_worst) .and. (min_rmsd == 0.0d0)) then
                               call Acceleration(vals,coords,gradient)
                       else if ((min_rmsd .ge. threshold_RMSD).or.(reject_flag)) then
                               call Acceleration(vals,coords,gradient)
                       else
                               !And we need to use the approximation after "unlabeling"
                               do n = 1, Natoms
                                        gradient(:,BOND_LABELLING_DATA(n)) = approx_gradient(:,n)
                               end do
                       end if
                end if

                !Update the velocities
                velocities = velocities + gradient

        end do

        if (testtrajDetailedRMSD_flag) close(filechannel2)
        deallocate(valsbuffer1,coordsbuffer1,gradientbuffer1,Ubuffer1,RMSDbuffer1,&
                   approximation_index)

        coords_final = coords
        velocities_final = velocities

end subroutine checkMultipleTrajectories





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	SUBROUTINE
!		Acceleration
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	PURPOSE
!		This subroutine calculates the gradient of a particular frame and
!		is governed by bond information and parameters in PHYSICS
!
!		Variables are supplied in case they have values useful for calcualtion
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	INPUT				KIND				DESCRIPTION
!
!		vals				REAL(DP),DIM(Nvar)		The variables associated with a frame
!		coords				REAL(DP),DIM(3,Natoms)		The coordinates defining a frame
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	OUTPUT				KIND				DESCRIPTION
!
!		gradient			REAL(DP),DIM(3,Natoms)		The gradient associated with a frame
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	IMPORTANT VARIABLES		KIND				DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES	                          FILETYPE			DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine Acceleration_generic(vals,coords,gradient)
	use PARAMETERS
        use PHYSICS
        implicit none

	!Frame information
        real(dp), dimension(Nvar), intent(in) :: vals
        real(dp), dimension(3,Natoms), intent(in) :: coords
        real(dp), dimension(3,Natoms), intent(out) :: gradient

	!Indexing of the atoms of the frame
        integer :: start_index1, start_index2
        integer :: index1, index2, bond_index1, bond_index2

	!Incremental integer
	integer :: i

	!Set the gradient to zero (we add on per pair of interactions)
        gradient = 0.0d0

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                 NON-BONDED INTERACTIONS
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!BONDING_DATA is arranged by ascending order of the first index
	!and the second index is always larger than the first index
	!So the first pair would naturally be (1,2)
        start_index1 = 1
        start_index2 = 2

	!We do a loop through every pair of atoms BETWEEN the bonds
        do i = 1, Nbonds

		!Figure out where the next pair of atoms forms a bond
		!And iterate over pairs of atoms until you reach that bond
                bond_index1 = BONDING_DATA(i,1)
                do index1 = start_index1, bond_index1

			!If we are almost there, then we only do pairs up until the bond index
                        if (index1 == bond_index1) then
                                bond_index2 = BONDING_DATA(i,2) - 1

			!Otherwise, we must consider all pairs of atoms
                        else
                                bond_index2 = Natoms
                        end if

                        do index2 = start_index2, bond_index2
                                !Remark: the optional 5th argument is the distance between
                                !        coord1 and coord2 (so it doesn't have to recalculate)
				!Right now, this is also not generic
                                call NonBondedForce(coords(:,index1),coords(:,index2),&
                                                gradient(:,index1),gradient(:,index2),vals(index2-1))
                        end do

			!Start off the next iteration from the index1 + 1
			!because index2 is always greater than index1
                        start_index2 = index1 + 2
                end do

		!Start off the next iteration over pairs of atoms
		!while skipping the pair of atoms that are bonded
                start_index1 = bond_index1
                start_index2 = bond_index2 + 2
        end do

	!Repeat the inner loop of the above; this loop takes into account the
	!pairs of atoms that have indexes larger than the last pair of atoms that are bonded
	do index1 = start_index1,Natoms
		do index2 = start_index2,Natoms
                        call NonBondedForce(coords(:,index1),coords(:,index2),&
                                        gradient(:,index1),gradient(:,index2),vals(i))
		end do
		start_index2 = index1 + 2
	end do


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                 BONDED INTERACTIONS
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!Finally, just calculate the forces between the bonded atoms
        do i = 1, Nbonds
                index1 = BONDING_DATA(i,1)
                index2 = BONDING_DATA(i,2)
                call BondedForce(coords(:,index1),coords(:,index2),&
                             gradient(:,index1),gradient(:,index2))
        end do

	return

end subroutine Acceleration_generic


subroutine Acceleration(vals,coords,gradient)
	use PARAMETERS
        use PHYSICS
        implicit none

	!Frame information
        real(dp), dimension(Nvar), intent(in) :: vals
        real(dp), dimension(3,Natoms), intent(in) :: coords
        real(dp), dimension(3,Natoms), intent(out) :: gradient

	!Indexing of the atoms of the frame
        integer :: index1, index2
        integer :: bond_index, bond_index1, bond_index2
        integer :: value_index, value_index1, value_index2
	logical :: bond_flag, value_flag

	!Set the gradient to zero (we add on per pair of interactions)
        gradient = 0.0d0

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                 NON-BONDED INTERACTIONS
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	bond_index = 1
	bond_flag = .true.
	bond_index1 = BONDING_DATA(1,1)
	bond_index2 = BONDING_DATA(1,2)

	value_index = 1
	if (Nvar_eff > 0) then
		value_flag = .true.
		value_index1 = BONDING_VALUE_DATA(1,1)
		value_index2 = BONDING_VALUE_DATA(1,2)
	else
		value_flag = .false.
	end if

	do index1 = 1, Natoms
		do index2 = index1+1, Natoms

			if ((bond_flag).and.(index1 == bond_index1).and.(index2 == bond_index2)) then
                                !Remark: the optional 5th argument is the distance between
                                !        coord1 and coord2 (so it doesn't have to recalculate)
				!Right now, this is also not generic
		                call BondedForce(coords(:,index1),coords(:,index2),&
		                             gradient(:,index1),gradient(:,index2))

				bond_index = bond_index + 1
				if (bond_index > Nbonds) then
					bond_flag = .false.
					cycle
				end if
				bond_index1 = BONDING_DATA(bond_index,1)
				bond_index2 = BONDING_DATA(bond_index,2)

			else
				if ((value_flag).and.(index1==value_index1).and.(index2==value_index2)) then
                                !Remark: the optional 5th argument is the distance between
                                !        coord1 and coord2 (so it doesn't have to recalculate)
				!Right now, this is also not generic
                                call NonBondedForce(coords(:,index1),coords(:,index2),&
                                                gradient(:,index1),gradient(:,index2),&
                                                vals(BONDING_VALUE_DATA(value_index,3)))

				value_index = value_index + 1
				if (value_index > Nvar_eff) then
					value_flag = .false.
					cycle
				end if
				value_index1 = BONDING_VALUE_DATA(value_index,1)
				value_index2 = BONDING_VALUE_DATA(value_index,2)

				else
                                call NonBondedForce(coords(:,index1),coords(:,index2),&
                                                gradient(:,index1),gradient(:,index2))
				end if
			end if
		end do
	end do
	
	return

end subroutine Acceleration



!
!subroutine addMultipleTrajectories()
!	use VARIABLES
!	use PHYSICS
!	use PARAMETERS
!	use interactSingleGrid
!	use analyzeHeatMapswithMultipleGrids
!        implicit none
!
!	!Coordinates, Velocities, and Variables
!	real(dp), dimension(3,Natoms,Ntraj_max) :: coords,gradient,velocities
!	real(dp), dimension(3,Natoms) :: coords_initial,velocities_initial
!	real(dp), dimension(Nvar,Ntraj_max) :: vals
!	logical,dimension(Ntraj_max) :: TRAJECTORIES_FLAG
!	integer :: current_cell
!	character(5) :: variable_length_text
!
!	!Incremental Integer
!	integer :: n,m
!
!        !Initialize the scene
!	Norder1 = 0
!	TRAJECTORIES_FLAG = .true.
!
!	if (heatmap_evolution_flag) then
!		call system("mkdir "//gridpath1//"png")
!		heatmap_steps = 0
!	end if
!
!	!aaa$OMP PARALLEL DO private(n)
!	do n = 1, Ntraj_max
!	        call InitialSetup3(coords(:,:,n),velocities(:,:,n))
!		call getVarsMaxMin(coords(:,:,n),Natoms,vals(:,n),Nvar,BOND_LABELLING_DATA)
!	        call Acceleration(vals(:,n),coords(:,:,n),gradient(:,:,n))
!	end do
!	!aaa$OMP END PARALLEL DO
!
!        !Accelerate the velcocities for a half step (verlet)
!	velocities = velocities + 0.5d0 * gradient
!
!	!To randomize the periods of the bond, I let the scene go on
!	!for a small period of time (need to standardize this later)
!	!aaa$OMP PARALLEL DO private(n,m,coords_initial,velocities_initial,steps)
!	do m = 1, Ntraj_max
!	coords_initial = coords(:,:,m)
!	velocities_initial = velocities(:,:,m)
!
!	do n = 1, Nbonds
!		do steps = 1, int(rand()*vib_period)
!			coords(:,:,m) = coords(:,:,m) + dt * velocities(:,:,m)
!			call Acceleration(vals(:,m),coords(:,:,m),gradient(:,:,m))
!			velocities(:,:,m) = velocities(:,:,m) + gradient(:,:,m)
!		end do
!
!		!And then reset the bond
!		coords(:,BONDING_DATA(n,1),m) = coords_initial(:,BONDING_DATA(n,1))
!		coords(:,BONDING_DATA(n,2),m) = coords_initial(:,BONDING_DATA(n,2))
!		velocities(:,BONDING_DATA(n,1),m) = velocities_initial(:,BONDING_DATA(n,1))
!		velocities(:,BONDING_DATA(n,2),m) = velocities_initial(:,BONDING_DATA(n,2))
!	end do
!	end do
!	!aaa$OMP END PARALLEL DO
!
!	!Now we go into the mainloop
!	!We have a hard cap of Nsteps timesteps
!	!aaa$OMP PARALLEL DO private(n,m,steps)
!        do steps = 1, Nsteps
!
!		!Just for bug-testing
!                if (.false.) then !(modulo(steps,50) == 0) then
!                        open(filechannel1,file=gridpath0//trajectoryfile,position="append")
!                        write(filechannel1,'(I6)') Natoms
!                        write(filechannel1,*) ""
!			do n = 1, Natoms
!                        write(filechannel1,'(A1,3F10.6)') 'H',&
!                                coords(1,n,1), coords(2,n,1), coords(3,n,1)
!			end do
!			close(filechannel1)
!                end if
!
!                if (heatmap_evolution_flag .and. (modulo(steps,heatmap_evolution_steps) == 0)) then
!			heatmap_steps = heatmap_steps + 1
!			write(variable_length_text,'(I0.5)') heatmap_steps
!			call analyzeHeatMaps2(counter0,counter1,"png/"//variable_length_text//".png")
!                end if
! 
!	do n = 1, Ntraj_max
!			if (.not.TRAJECTORIES_FLAG(n)) cycle
!
!	                !Check every 500 steps if we are out-of-bounds
!	                if (modulo(steps,500) == 1) then      
!	 			if ((vals(1,n)>max_var1) .or. (vals(2,n)>max_var2)) then
!	 				TRAJECTORIES_FLAG(n) = .false.
!	 			end if
!
!				if (.not.(any(TRAJECTORIES_FLAG))) header_max_flag = .true.
!	                endif
!	end do
!
!		!aaa$OMP PARALLEL
!		do n = 1, Ntraj_max
!			if (.not.TRAJECTORIES_FLAG(n)) cycle
!	
!	                !Update the coordinates with the velocities
!			coords(:,:,n) = coords(:,:,n) + dt * velocities(:,:,n)
!	
!			!Update the variables
!			call getVarsMaxMin(coords(:,:,n),Natoms,vals(:,n),Nvar,BOND_LABELLING_DATA)
!	
!	                !Accelerate and update gradients
!	                call Acceleration(vals(:,n),coords(:,:,n),gradient(:,:,n))
!	
!			!Add the frame to the grid
!	        	call addState(vals(:,n),coords(:,:,n),gradient(:,:,n))
!	
!			!If there are too many subdivisions and counter0 will get out-of-bounds
!			!we would have to call this to exit
!	                if (header_max_flag) exit
!	
!			!Update the velocities(:,:,n)
!			velocities(:,:,n) = velocities(:,:,n) + gradient(:,:,n)
!		end do
!		!aaa$OMP END PARALLEL
!
!	        if (header_max_flag) then
!			heatmap_steps = heatmap_steps + 1
!			write(variable_length_text,'(I0.5)') heatmap_steps
!			call analyzeHeatMaps2(counter0,counter1,"png/"//variable_length_text//".png")
!			exit
!		end if
!        end do
!	!aaa$OMP END PARALLEL DO
!
!	if (heatmap_evolution_flag) then
!		call system("convert -delay 20 -loop 0 "//gridpath1//"png/*.png "//gridpath1//"heatmap_evolution.gif")
!	end if
!
!end subroutine addMultipleTrajectories
!
!
!
!subroutine addMultipleTrajectories2()
!	use VARIABLES
!	use PHYSICS
!	use PARAMETERS
!	use interactSingleGrid
!	use analyzeHeatMapswithMultipleGrids
!        implicit none
!
!	!Coordinates, Velocities, and Variables
!	real(dp), dimension(3,Natoms,Ntraj_max) :: coords,gradient,velocities
!	real(dp), dimension(3,Natoms) :: coords_initial,velocities_initial
!	real(dp), dimension(Nvar,Ntraj_max) :: vals
!	integer,dimension(counter0_max) :: TRAJECTORIES0
!	integer,dimension(resolution_0) :: TRAJECTORIES1
!	integer,dimension(Ntraj_max) :: NEIGHBOR_LIST
!	logical,dimension(Ntraj_max) :: TRAJECTORIES_FLAG
!	integer :: current_cell
!	integer :: var1_index, var2_index,indexer
!	real :: var1_round0, var2_round0
!	character(5) :: variable_length_text
!
!	!Incremental Integer
!	integer :: n,m
!
!        !Initialize the scene
!	Norder1 = 0
!	TRAJECTORIES_FLAG = .true.
!
!	if (heatmap_evolution_flag) then
!		call system("mkdir "//gridpath1//"png")
!		heatmap_steps = 0
!	end if
!
!        !Initialize all the trajectories
!	!aaa$OMP PARALLEL DO private(n)
!	do n = 1, Ntraj_max
!	        call InitialSetup3(coords(:,:,n),velocities(:,:,n))
!		call getVarsMaxMin(coords(:,:,n),Natoms,vals(:,n),Nvar,BOND_LABELLING_DATA)
!	        call Acceleration(vals(:,n),coords(:,:,n),gradient(:,:,n))
!	end do
!	!aaa$OMP END PARALLEL DO
!
!        !Accelerate the velcocities for a half step (verlet)
!	velocities = velocities + 0.5d0 * gradient
!
!	!To randomize the periods of the bond, I let the scene go on
!	!for a small period of time (need to standardize this later)
!	!aaa$OMP PARALLEL DO private(n,m,steps)
!	do m = 1, Ntraj_max
!	coords_initial = coords(:,:,m)
!	velocities_initial = velocities(:,:,m)
!
!	do n = 1, Nbonds
!		do steps = 1, int(rand()*vib_period)
!			coords(:,:,m) = coords(:,:,m) + dt * velocities(:,:,m)
!			call Acceleration(vals(:,m),coords(:,:,m),gradient(:,:,m))
!			velocities(:,:,m) = velocities(:,:,m) + gradient(:,:,m)
!		end do
!
!		!And then reset the bond
!		coords(:,BONDING_DATA(n,1),m) = coords_initial(:,BONDING_DATA(n,1))
!		coords(:,BONDING_DATA(n,2),m) = coords_initial(:,BONDING_DATA(n,2))
!		velocities(:,BONDING_DATA(n,1),m) = velocities_initial(:,BONDING_DATA(n,1))
!		velocities(:,BONDING_DATA(n,2),m) = velocities_initial(:,BONDING_DATA(n,2))
!	end do
!	end do
!	!aaa$OMP END PARALLEL DO
!
!	TRAJECTORIES0 = 0
!	current_cell = 0
!	NEIGHBOR_LIST = 0
!	do n = 1, Ntraj_max
!		call getVarsMaxMin(coords(:,:,n),Natoms,vals(:,n),Nvar,BOND_LABELLING_DATA)
!
!		var1_index = int(vals(1,n) * divisor1_0)
!		var2_index = int(vals(2,n) * divisor2_0)
!		
!		var1_round0 = multiplier1_0 * var1_index
!		var2_round0 = multiplier2_0 * var2_index
!		
!		indexer = bounds1 * var2_index + var1_index + 1
!
!		if (TRAJECTORIES0(indexer) == 0) then
!			current_cell = current_cell + 1
!			NEIGHBOR_LIST(n) = -indexer
!			TRAJECTORIES0(indexer) = n
!		else
!			NEIGHBOR_LIST(n) = TRAJECTORIES0(indexer)
!			TRAJECTORIES0(indexer) = n
!		end if
!	end do
!
!	!Now we go into the mainloop
!        do
!
!		!Just for bug-testing
!                if (.false.) then !(modulo(steps,50) == 0) then
!                        open(filechannel1,file=gridpath0//trajectoryfile,position="append")
!                        write(filechannel1,'(I6)') Natoms
!                        write(filechannel1,*) ""
!			do n = 1, Natoms
!                        write(filechannel1,'(A1,3F10.6)') 'H',&
!                                coords(1,n,1), coords(2,n,1), coords(3,n,1)
!			end do
!			close(filechannel1)
!                end if
!
!                if (heatmap_evolution_flag .and. (modulo(steps,heatmap_evolution_steps) == 0)) then
!			heatmap_steps = heatmap_steps + 1
!			write(variable_length_text,'(I0.5)') heatmap_steps
!			call analyzeHeatMaps2(counter0,counter1,"png/"//variable_length_text//".png")
!                end if
! 
!                !Check every 500 steps if we are out-of-bounds
!                if (modulo(steps,500) == 1) then      
!		do n = 1, Ntraj_max
!			if (.not.TRAJECTORIES_FLAG(n)) cycle
!
!			if ((vals(1,n)>max_var1) .or. (vals(2,n)>max_var2)) then
!				TRAJECTORIES_FLAG(n) = .false.
!				NEIGHBOR_LIST(n) = 0
!			end if
!
!			if (.not.(any(TRAJECTORIES_FLAG))) header_max_flag = .true.
!		end do
!                endif
!
!	        !aaa$OMP PARALLEL DO
!		do
!			if (.not.TRAJECTORIES_FLAG(n)) cycle
!	
!	                !Update the coordinates with the velocities
!			coords(:,:,n) = coords(:,:,n) + dt * velocities(:,:,n)
!	
!			!Update the variables
!			call getVarsMaxMin(coords(:,:,n),Natoms,vals(:,n),Nvar,BOND_LABELLING_DATA)
!	
!	                !Accelerate and update gradients
!	                call Acceleration(vals(:,n),coords(:,:,n),gradient(:,:,n))
!	
!			!Add the frame to the grid
!	        	call addState(vals(:,n),coords(:,:,n),gradient(:,:,n))
!	
!			!If there are too many subdivisions and counter0 will get out-of-bounds
!			!we would have to call this to exit
!	                if (header_max_flag) exit
!	
!			!Update the velocities(:,:,n)
!			velocities(:,:,n) = velocities(:,:,n) + gradient(:,:,n)
!		end do
!	        !aaa$OMP END PARALLEL DO
!
!	        if (header_max_flag) then
!			heatmap_steps = heatmap_steps + 1
!			write(variable_length_text,'(I0.5)') heatmap_steps
!			call analyzeHeatMaps2(counter0,counter1,"png/"//variable_length_text//".png")
!			exit
!		end if
!        end do
!
!	if (heatmap_evolution_flag) then
!		call system("convert -delay 20 -loop 0 "//gridpath1//"png/*.png "//gridpath1//"heatmap_evolution.gif")
!	end if
!
!end subroutine addMultipleTrajectories2
!
!





end module runTrajectory
