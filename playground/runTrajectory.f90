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
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	SUBROUTINES			ARGUMENTS		KIND
!
!		addTrajectory			velocityH1		intent(out),real(dp)
!						velocityH2		intent(out),real(dp)
!
!		checkTrajectory			velocityH1		intent(out),real(dp)
!						velocityH2		intent(out),real(dp)
!
!		checkMultipleTrajectories	filechannels		intent(in),dim(Ngrid_total),integer
!
!						velocityH1		intent(out),real(dp)
!						velocityH2		intent(out),real(dp)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       CALLS                           MODULE
!
!		getVar5				VARIABLES
!		getVar6				VARIABLES
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
!		ScatteringAngle			REAL(DP)			The scatteringangle
!		TRVenergies1			REAL(DP),DIM(3)			The translation, rotational, and vibrational
!										speeds of the incoming H/H2 molecule
!										at the start of the trajectory
!		TRVenergies2			REAL(DP),DIM(3)			The translation, rotational, and vibrational
!										speeds of the incoming H/H2 molecule
!										at the end of the trajectory (post-collision)
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

subroutine addTrajectory(coords_initial,velocities_initial,coords_final,velocities_final)
	use VARIABLES
	use PHYSICS
	use PARAMETERS
	use interactSingleGrid
        implicit none

	!Coordinates, Velocities, and Variables
	real(dp), dimension(3,Natoms) :: coords,gradient,velocities
	real(dp), dimension(Nvar) :: vals
	real(dp),dimension(3,Natoms),intent(out) :: coords_initial,velocities_initial
	real(dp),dimension(3,Natoms),intent(out) :: coords_final,velocities_final
	integer :: bond_index1, bond_index2

	!Incremental Integer
	integer :: n

        !Initialize the scene
        call InitialSetup3(coords,velocities)
	Norder1 = 0

	coords_initial = coords
	velocities_initial = velocities

	!Always calculate the variables before accelerating
	!because we can reuse these calculations
	call getVar5(coords,Natoms,vals(1))
	call getVar6(coords,Natoms,vals(2))

        !Accelerate the velcocities for a half step (verlet)
        call Acceleration(vals,coords,gradient)
	call addState(vals,coords,gradient)

	!Update the velocities
	velocities = velocities + 0.5d0 * gradient

	!To randomize the periods of the bond, I let the scene go on
	!for a small period of time (need to standardize this later)
	do n = 1, Nbonds
		do steps = 1, int(INITIAL_BOND_DATA(n,6)*vib_period)
			coords = coords + dt * velocities
			call Acceleration(vals,coords,gradient)
			velocities = velocities + 0.5d0 * gradient
		end do

		!And then reset the bond
		coords(:,BONDING_DATA(n,1)) = coords_initial(:,BONDING_DATA(n,1))
		coords(:,BONDING_DATA(n,2)) = coords_initial(:,BONDING_DATA(n,2))
	end do

	!Now we go into the mainloop
	!We have a hard cap of Nsteps timesteps
        do steps = 1, Nsteps

		!Just for bug-testing
                if (.false.) then !(modulo(steps,50) == 0) then
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
 			if ((vals(1)>max_var1) .or. (vals(2)>max_var2)) then
 				exit
 			end if
                endif

                !Update the coordinates with the velocities
		coords = coords + dt * velocities

		!Always calculate the variables before adding a frame or accelerating
		call getVar5(coords,Natoms,vals(1))
		call getVar6(coords,Natoms,vals(2))

                !Accelerate and update gradients
                call Acceleration(vals,coords,gradient)

		!Add the frame to the grid
        	call addState(vals,coords,gradient)

		!If there are too many subdivisions and counter0 will get out-of-bounds
		!we would have to call this to exit
                if (header_max_flag) exit

		!Update the velocities
		velocities = velocities + gradient
        end do

	coords_final = coords
	velocities_final = velocities

end subroutine addTrajectory



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
!		ScatteringAngle			REAL(DP)			The scatteringangle
!		TRVenergies1			REAL(DP),DIM(3)			The translation, rotational, and vibrational
!										speeds of the incoming H/H2 molecule
!										at the start of the trajectory
!		TRVenergies2			REAL(DP),DIM(3)			The translation, rotational, and vibrational
!										speeds of the incoming H/H2 molecule
!										at the end of the trajectory (post-collision)
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


subroutine checkTrajectory(coords_initial,velocities_initial,coords_final,velocities_final)
        use PARAMETERS
        use PHYSICS
        use VARIABLES
        use ANALYSIS
        use interactSingleGrid
        implicit none

        !Coordinates, Velocities, and Variables
        real(dp), dimension(3,Natoms) :: coords,velocities
        real(dp), dimension(3,Natoms) :: gradient, approx_gradient
        real(dp), dimension(Nvar) :: vals
	real(dp), dimension(3,Natoms), intent(out) :: coords_initial, velocities_initial 
	real(dp), dimension(3,Natoms), intent(out) :: coords_final, velocities_final 
	integer :: bond_index1, bond_index2

        !Various other variables
        real(dp) :: U, KE
        real(dp) :: min_rmsd,min_rmsd_prime
        integer :: number_of_frames,order,neighbor_check

	!Incremental Integer
	integer :: n

        !Initialize the scene
        call InitialSetup3(coords,velocities)

	coords_initial = coords
	velocities_initial = velocities

	!Always calculate the variables before accelerating
        call getVar5(coords,Natoms,vals(1))
        call getVar6(coords,Natoms,vals(2))

        !Accelerate the velcocities for a half step (verlet)
        call Acceleration(vals,coords,gradient)

	!Update the velocities
        velocities = velocities + 0.5d0 * gradient

	!To randomize the periods of the bond, I let the scene go on
	!for a small period of time (need to standardize this later)
	do n = 1, Nbonds
		do steps = 1, int(INITIAL_BOND_DATA(n,6)*vib_period)
			coords = coords + dt * velocities
			call Acceleration(vals,coords,gradient)
			velocities = velocities + 0.5d0 * gradient
		end do

		!And then reset the bond
		coords(:,BONDING_DATA(n,1)) = coords_initial(:,BONDING_DATA(n,1))
		coords(:,BONDING_DATA(n,2)) = coords_initial(:,BONDING_DATA(n,2))
	end do

	!Get the trajectory file open for trajectory visualization
        open(filechannel1,file=gridpath1//trajectoryfile)
        write(filechannel1,'(I6)') Natoms
        write(filechannel1,*) ""
	do n = 1, Natoms
        	write(filechannel1,'(A1,3F10.6)') 'H',&
              		coords(1,n), coords(2,n), coords(3,n)
        end do
        close(filechannel1)

	!We keep this file open for the whole trajectory (instead of
	!continually opening and closing) to keep data of each frame
	open(filechannel2,file=gridpath1//checkstatefile)
        do steps = 1, Nsteps

                !Every 10 frames, print to an xyz file for visualization
                 if (modulo(steps,10) == 0) then
                        open(filechannel1,file=gridpath1//trajectoryfile,position="append")
                        write(filechannel1,'(I6)') Natoms
                        write(filechannel1,*) ""
			do n = 1, Natoms
                        	write(filechannel1,'(A1,3F10.6)') 'H',&
                              		coords(1,n), coords(2,n), coords(3,n)
                        end do
                        close(filechannel1)
                 end if

                !Check every 500 steps to see if we are out-of-bounds
                if (modulo(steps,500) == 1) then
                        if ((vals(1)>max_var1) .or.&
                            (vals(2)>max_var2)) then
                                exit
                        end if
                endif

                !Upate the coordinates with the velocities
                coords = coords + dt * velocities

                !Always calculate the variables before checking a frame or accelerating
                call getVar5(coords,Natoms,vals(1))
                call getVar6(coords,Natoms,vals(2))

                !Check for similar frames in the grid
		!Always set a default value; in this case, set min_rmsd_prime a default value
                min_rmsd_prime = default_rmsd
                call checkState(vals,coords,approx_gradient,min_rmsd_prime,&
                                number_of_frames,order,neighbor_check)
                min_rmsd = min_rmsd_prime

                !Update the gradient either with the approximation or by accelerating
		!This depends on the threshold and the rejection method
                if ((min_rmsd .ge. threshold_rmsd) .or. (reject_flag)) then
                        call Acceleration(vals,coords,gradient)
                else
                        gradient = approx_gradient
                end if

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

                !Update the velocities
                velocities = velocities + gradient

        end do
	close(filechannel2)

	coords_final = coords
	velocities_final = velocities

end subroutine checkTrajectory




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
!		filechannels			INTEGER,DIM(Ngrid_total)	These filechannels were opened beforehand to the
!										trajectory files that we output rmsd data to;
!										all that's left to do is to write to them every frame
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	OUTPUT				KIND				DESCRIPTION
!
!		ScatteringAngle			REAL(DP)			The scatteringangle
!		TRVenergies1			REAL(DP),DIM(3)			The translation, rotational, and vibrational
!										speeds of the incoming H/H2 molecule
!										at the start of the trajectory
!		TRVenergies2			REAL(DP),DIM(3)			The translation, rotational, and vibrational
!										speeds of the incoming H/H2 molecule
!										at the end of the trajectory (post-collision)
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
!		number_of_frames		INTEGER				The total amount of frames checked over all frames
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

subroutine checkMultipleTrajectories(filechannels,coords_initial,velocities_initial,&
                                                  coords_final,velocities_final)
        use PHYSICS
        use PARAMETERS
        use VARIABLES
	use ANALYSIS
        use interactMultipleGrids
        implicit none

        !Coordinates, Velocities, and Variables
        real(dp),dimension(3,Natoms) :: coords,velocities
        real(dp),dimension(3,Natoms) :: gradient,approx_gradient
        real(dp),dimension(3,Natoms),intent(out) :: coords_initial, velocities_initial
        real(dp),dimension(3,Natoms),intent(out) :: coords_final, velocities_final
	real(dp),dimension(Nvar) :: vals
	
	integer :: bond_index1, bond_index2

        !Grid Parameters
        integer,dimension(Ngrid_total),intent(in) :: filechannels

        !Various other variables
        real(dp) :: U, KE
        real(dp) :: min_rmsd,min_rmsd_prime
        integer :: number_of_frames,order,neighbor_check

	!Incremental Integer
	integer :: n

        !Initialize the scene
        call InitialSetup3(coords,velocities)

	coords_initial = coords
	velocities_initial = velocities

        !Always calculate the variables before accelearting
        call getVar5(coords,Natoms,vals(1))
        call getVar6(coords,Natoms,vals(2))

        !Accelerate the velcocities for a half step (verlet)
        call Acceleration(vals,coords,gradient)

	!Update the velocities
        velocities = velocities + 0.5d0 * gradient

	!To randomize the periods of the bond, I let the scene go on
	!for a small period of time (need to standardize this later)
	do n = 1, Nbonds
		do steps = 1, int(INITIAL_BOND_DATA(n,6)*vib_period)
			coords = coords + dt * velocities
			call Acceleration(vals,coords,gradient)
			velocities = velocities + 0.5d0 * gradient
		end do

		!And then reset the bond
		coords(:,BONDING_DATA(n,1)) = coords_initial(:,BONDING_DATA(n,1))
		coords(:,BONDING_DATA(n,2)) = coords_initial(:,BONDING_DATA(n,2))
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
                        if ((vals(1)>max_var1) .or.&
                            (vals(2)>max_var2)) then
                                exit
                        end if
                endif

                !Upate the coordinates with the velocities
                coords = coords + dt * velocities

                !Always calculate the variables before checking a frame or accelerating
                call getVar5(coords,Natoms,vals(1))
                call getVar6(coords,Natoms,vals(2))

		!Check for a frame in the grid
		!Set the default value beforehand though
                min_rmsd = default_rmsd
                call checkState(vals,coords,approx_gradient,min_rmsd,&
                         filechannels,number_of_frames,order,neighbor_check)

                !Update the gradient with either the approximation or by acclerating
		!This is dependent on the threshold and the rejection method
                if ((min_rmsd .ge. threshold_RMSD).or.(reject_flag)) then
                        call Acceleration(vals,coords,gradient)
                else
                        gradient = approx_gradient
                end if

                !Update the velocities
                velocities = velocities + gradient

        end do

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
	value_flag = .true.
	value_index1 = BONDING_VALUE_DATA(1,1)
	value_index2 = BONDING_VALUE_DATA(1,2)

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



end module runTrajectory
