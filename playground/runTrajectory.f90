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
use PARAMETERS
use ANALYSIS
implicit none

character(gridpath_length+expfolder_length) :: gridpath4
character(gridpath_length+expfolder_length+5) :: gridpath5

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
        use analyzeRMSDThresholdwithMultipleGrids
        implicit none

        !Coordinates, Velocities, and Variables
        real(dp), dimension(3,Natoms) :: coords,coords_labelled,velocities
        real(dp), dimension(3,Natoms) :: gradient,gradient_labelled
        real(dp), dimension(3,Natoms) :: approx_gradient,approx_gradient_prime
        real(dp), dimension(Nvar) :: vals
	real(dp), dimension(3,Natoms), intent(out) :: coords_initial, velocities_initial 
	real(dp), dimension(3,Natoms), intent(out) :: coords_final, velocities_final 
        integer,dimension(1+Ngrid_max) :: filechannels
        integer :: bond_index1, bond_index2

        !Various other variables
        real(dp) :: U, KE
        real(dp) :: min_rmsd,min_rmsd_prime
        integer :: number_of_frames,order,neighbor_check
        character(9) :: vals_interpolation_text

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
        if (testtrajDetailedRMSD_flag) open(filechannel2,file=&
                gridpath5//checkstatefile)
        if (gather_interpolation_flag) open(filechannel3,file=&
                gridpath5//interpolationfile,position="append")

        buffer1_size = 1 + var_overcrowd(1)
        allocate(valsbuffer1(Nvar,buffer1_size),&
                 coordsbuffer1(3,Natoms,buffer1_size),&
                 gradientbuffer1(3,Natoms,buffer1_size),&
                 Ubuffer1(3,3,buffer1_size),&
                 RMSDbuffer1(buffer1_size),&
                 approximation_index(buffer1_size))

!       if (interpolation_flag) then
                allocate(acceptable_frame_mask(buffer1_size),&
                         inputCLS(Ncoords+buffer1_size,buffer1_size))
                interpolation_counter = 0
!       end if

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

                        if (accept_worst) then
                                min_rmsd = 0.0d0
                        else
                                min_rmsd = default_rmsd
                        end if

                        subcellsearch_max = (/ 0, 0 /)
                        call checkState_new(vals,coords,approx_gradient,min_rmsd,&
                                 filechannels,number_of_frames,order,neighbor_check)

                        if ((min_rmsd .ge. threshold_RMSD).or.(reject_flag)&
                            .or.((accept_worst).and.(min_rmsd == 0.0d0))) then
                                call Acceleration(vals,coords,gradient)
                                if (grid_addition > 0) &
                                        call addState_new(vals,coords,gradient)
                        else
                                gradient = approx_gradient
                        end if

                !Otherwise, then we need to do the label switching
                else
                        !We need to input the frame using the labeling scheme
                        do n = 1, Natoms
                                coords_labelled(:,n) = coords(:,BOND_LABELLING_DATA(n))
!                               gradient_labelled(:,n) = gradient(:,BOND_LABELLING_DATA(n))
                        end do

                       !Check for a frame in the grid
                       !Set the default value beforehand though
                       if (accept_worst) then
                               min_rmsd = 0.0d0
                       else
                               min_rmsd = default_rmsd
                       end if

                       if (testtrajDetailedRMSD_flag) then
                               subcellsearch_max = subcellsearch_max1
                               if (gather_interpolation_flag) then 
                                       interpolation_flag = .false.
                               end if
                       end if
                       call checkState_new(vals,coords_labelled,approx_gradient,min_rmsd,&
                                filechannels,number_of_frames,order,neighbor_check)

if (testtrajDetailedRMSD_flag) then

        if (accept_worst) then
                min_rmsd_prime = 0.0d0
        else
                min_rmsd_prime = default_rmsd
        end if

        subcellsearch_max = subcellsearch_max2
        if (gather_interpolation_flag) interpolation_flag = .true.
        call checkState_new(vals,coords_labelled,approx_gradient_prime,min_rmsd_prime,&
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
else
        approx_gradient_prime = approx_gradient
end if

                       !Update the gradient with either the approximation or by acclerating
                       !This is dependent on the threshold and the rejection method
                       if ((min_rmsd .ge. threshold_RMSD).or.(reject_flag)&
                            .or.((accept_worst).and.(min_rmsd == 0.0d0))) then
                               call Acceleration(vals,coords,gradient)
                               if (grid_addition > 0) then
                                       do n = 1, Natoms
                                                gradient_labelled(:,n) = &
                                                        gradient(:,BOND_LABELLING_DATA(n))
                                                coords_labelled(:,n) = &
                                                        coords(:,BOND_LABELLING_DATA(n))
                                       end do
                                       call addState_new(vals,&
                                               coords_labelled,&
                                               gradient_labelled)
                               end if
                       else
                               !And we need to use the approximation after "unlabeling"
                               do n = 1, Natoms
                                        gradient(:,BOND_LABELLING_DATA(n)) = &
                                                approx_gradient(:,n)
                               end do
                               approx_gradient = gradient

                               do n = 1, Natoms
                                        gradient(:,BOND_LABELLING_DATA(n)) = &
                                                approx_gradient_prime(:,n)
                               end do
                               approx_gradient_prime = gradient

                               if ((interpolation_flag).and.(gather_interpolation_flag)) then

                                        call Acceleration(vals,coords,gradient)

                                        if (Ninterpolation > 0) then
                                        write(filechannel3,FMT=*) vals(1), vals(2), Ninterpolation,&
                                                    largest_weighted_rmsd2, largest_weighted_rmsd,&
                                                    min_rmsd,&
                                                    sqrt(sum((gradient-approx_gradient)**2)/Natoms),&
                                                    min_rmsd_prime,&
                                                    sqrt(sum((gradient-approx_gradient_prime)**2)/Natoms)

                                        interpolation_counter = interpolation_counter + 1

!                                       if ((sqrt(sum((gradient-approx_gradient_prime)**2)/Natoms)&
!                                               > .0001).and.(min_rmsd_prime < .005)) then
                                        if ((sqrt(sum((gradient-approx_gradient_prime)**2)/Natoms)&
                                                > .000001)) then
                                               write(6,FMT=*) ""
                                               write(6,FMT="(A,I7)") "Step: ", steps
                                               write(6,FMT="(A,F9.6,1x,F9.6)") &
                                                       "(mu1,E): ", min_rmsd_prime,&
                                               sqrt(sum((gradient-approx_gradient_prime)**2)/Natoms)
                                               write(6,FMT="(A,F9.6)") "rmsd:", largest_rmsd
                                               write(6,FMT="(A,F7.3,1x,F7.3)") "vals: ", vals
                                               write(6,FMT="(A,I5)") "Ninter: ", Ninterpolation
                                        end if

                                        end if

!                                       if ((interpolation_counter > interpolation_check)&
!                                           .and.(min_rmsd < 5.0e-5)) then
!                                       if (interpolation_counter == interpolation_check) then
                                        if (.false.) then
!                                       if (min_rmsd < 5.0e-5) then

                                                close(filechannel3)

!                                               write(vals_interpolation_text,&
!                                                     FMT="(F4.1,'_',F4.1)") vals
                                                if (interpolation_check_visual) call &
                                                        getRMSDinterpolation(0.1d0*anint(vals*10),&
!                                                       (/0.1d0,0.1d0/),&
!                                                       vals_interpolation_text//"InterpolationCheck")
                                                        (/0.1d3,0.1d3/),&
                                                        "InterpolationCheck")

                                                open(filechannel3,file=gridpath5//interpolationfile,&
                                                        position="append")

                                                interpolation_counter = 0

                                        end if
                               end if
                       end if
                end if

                !Update the velocities
                velocities = velocities + gradient
        end do

        if (testtrajDetailedRMSD_flag) close(filechannel2)
        if (gather_interpolation_flag) close(filechannel3)

        deallocate(valsbuffer1,coordsbuffer1,gradientbuffer1,Ubuffer1,RMSDbuffer1,&
                   approximation_index)

!       if (interpolation_flag) deallocate(acceptable_frame_mask,inputCLS)

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
        use analyzeRMSDThresholdwithMultipleGrids
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

        if (testtrajDetailedRMSD_flag) open(filechannel2,file=&
                gridpath5//checkstatefile)
        if (gather_interpolation_flag) open(filechannel3,file=&
                gridpath5//interpolationfile,position="append")

        buffer1_size = 1 + var_overcrowd(1)
        allocate(valsbuffer1(Nvar,buffer1_size),&
                 coordsbuffer1(3,Natoms,buffer1_size),&
                 gradientbuffer1(3,Natoms,buffer1_size),&
                 Ubuffer1(3,3,buffer1_size),&
                 RMSDbuffer1(buffer1_size),&
                 approximation_index(buffer1_size))

!       if (interpolation_flag) then
                allocate(acceptable_frame_mask(buffer1_size),&
                         inputCLS(Ncoords+buffer1_size,buffer1_size))
                interpolation_counter = 0
!       end if

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

                       if (testtrajDetailedRMSD_flag) then
                               subcellsearch_max = subcellsearch_max1
                               if (gather_interpolation_flag) then 
                                       interpolation_flag = .false.
                               end if
                       end if
                       call checkState_new(vals,coords_labelled,approx_gradient,min_rmsd,&
                                filechannels,number_of_frames,order0,neighbor_check)

if (testtrajDetailedRMSD_flag) then

        if (accept_worst) then
                min_rmsd_prime = 0.0d0
        else
                min_rmsd_prime = default_rmsd
        end if

        subcellsearch_max = subcellsearch_max2
        if (gather_interpolation_flag) interpolation_flag = .true.

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

                               if ((interpolation_flag).and.(gather_interpolation_flag)) then

                                        approx_gradient = gradient

                                        call Acceleration(vals,coords,gradient)

!                                       write(filechannel3,FMT=*) vals(1), vals(2), Ninterpolation,&
!                                                   largest_rmsd, threshold_rmsd, min_rmsd, &
!                                                   sqrt(sum((gradient-approx_gradient)**2)/Natoms)
                                        write(filechannel3,FMT=*) vals(1),vals(2),Ninterpolation,&
                                                candidate_rmsd,min_rmsd,&
                                                sqrt(sum((gradient-candidate_gradient)**2)/Natoms),&
                                                sqrt(sum((gradient-approx_gradient)**2)/Natoms)

                                        interpolation_counter = interpolation_counter + 1

                                        if (.false.) then
!                                       if (interpolation_counter == interpolation_check) then

                                                close(filechannel3)

                                                if (interpolation_check_visual) call &
                                                        getRMSDinterpolation((/5.0d0,5.5d0/),&
                                                        (/0.1d0,0.1d0/),"InterpolationCheck")

                                                open(filechannel3,file=gridpath5//interpolationfile,&
                                                        position="append")

                                                interpolation_counter = 0

                                        end if
                               end if
                       end if
                end if

                !Update the velocities
                velocities = velocities + gradient

        end do

        if (testtrajDetailedRMSD_flag) close(filechannel2)
        if (gather_interpolation_flag) close(filechannel3)

        deallocate(valsbuffer1,coordsbuffer1,gradientbuffer1,Ubuffer1,RMSDbuffer1,&
                   approximation_index)

!       if (interpolation_flag) &
                deallocate(acceptable_frame_mask,inputCLS)

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


subroutine errorCheck1(filechannels)
        use PARAMETERS
        use PHYSICS
        use VARIABLES
        use ANALYSIS
        use ls_rmsd_original
        use interactMultipleGrids
        use analyzeRMSDThresholdwithMultipleGrids
        implicit none

        !Coordinates, Velocities, and Variables
        real(dp), dimension(3,Natoms) :: coords,coords_labelled,velocities, delta_coords
        real(dp), dimension(3,Natoms) :: gradient,gradient_var,gradient_labelled,gradient_var_labelled
        real(dp), dimension(3,Natoms) :: approx_gradient,approx_gradient_prime
        real(dp), dimension(Nvar) :: vals
	real(dp), dimension(3,Natoms) :: coords_initial, velocities_initial 
	real(dp), dimension(3,Natoms) :: coords_final, velocities_final 
        integer :: bond_index1, bond_index2

        real(dp) :: handicap_rmsd
        real(dp), allocatable :: rmsd_x(:,:),rmsd_x_interpolated(:,:)
        real(dp), allocatable :: rmsd_fx(:,:),rmsd_fx_interpolated(:,:)
        real(dp),dimension(4) :: selected_means,selected_SDs
        real(dp) :: delta_length
        integer :: Ntest,Nsamples,Nsample
        integer :: Nanomaly
        character(3) :: Nanomaly_text
        real(dp), allocatable :: rmsd_weights(:,:,:), rmsd_fx_weights(:,:,:)
        real(dp), allocatable :: mean_weights(:), mean_rmsd_fx(:)

        !Various other variables
        real(dp) :: min_rmsd,min_rmsd_prime
        integer :: number_of_frames,order,neighbor_check
        character(9) :: vals_interpolation_text

        integer,dimension(1+Ngrid_total),intent(in) :: filechannels
        real(dp), dimension(3) :: x_center, y_center
        real(dp), allocatable :: g(:,:)
        real(dp),dimension(3,3) :: U
        real(dp),dimension(3,3) :: candidate_U

        !Incremental Integer
        integer :: i,n

        !Initialize the scene
        call InitialSetup3(coords,velocities)

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
if (gather_interpolation_flag) open(filechannel3,file=gridpath5//interpolationfile,position="append")

        buffer1_size = 1 + var_overcrowd(1)
        allocate(valsbuffer1(Nvar,buffer1_size),&
                 coordsbuffer1(3,Natoms,buffer1_size),&
                 gradientbuffer1(3,Natoms,buffer1_size),&
                 Ubuffer1(3,3,buffer1_size),&
                 RMSDbuffer1(buffer1_size),&
                 approximation_index(buffer1_size))

        if (interpolation_flag) then
                allocate(acceptable_frame_mask(buffer1_size),&
                         inputCLS(Ncoords+buffer1_size,buffer1_size))
                interpolation_counter = 0
        end if

        Nanomaly = 0
        Ntest = 10
        Nsamples = 10
        allocate(rmsd_x_interpolated(Ntest,Nsamples),rmsd_fx(Ntest,Nsamples),&
                 rmsd_fx_interpolated(Ntest,Nsamples),rmsd_x(Ntest,Nsamples),&
                 rmsd_weights(Ntest,Ntest,Nsamples),rmsd_fx_weights(Ntest,Ntest,Nsamples))
        allocate(temp_frame_weights(Ntest),temp_rmsd_weights(Ntest))
        rmsd_weights = 0.0d0

        do

        do steps = 1, 10000
                !Upate the coordinates with the velocities
                coords = coords + dt * velocities

                !Always calculate the variables before checking a frame or accelerating
                call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

                call Acceleration(vals,coords,gradient)

                !Update the velocities
                velocities = velocities + gradient
        end do

        if (any(vals > var_maxvar)) exit

        coords_final = coords
        velocities_final = velocities

        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
        call Acceleration(vals,coords,gradient)

        do n = 1, Natoms
                coords_labelled(:,n) = coords(:,BOND_LABELLING_DATA(n))
                gradient_labelled(:,n) = gradient(:,BOND_LABELLING_DATA(n))
        end do
        print *, "vals:", vals
        print *, ""

        open(filechannel2,file=gridpath5//"new"//errorcheckfile)

        do Nsample = 1, Nsamples

        print *, "Starting sampling:", Nsample
        call system("rm "//gridpath2//"*.dat")
        handicap_rmsd = threshold_rmsd* 0.1d0

        !Let's find one good point
        do
                do n = 1, Natoms
                do i = 1, 3
                        delta_coords(i,n) = rand() - 0.5d0
                end do
                end do

                delta_length = sqrt(sum(delta_coords**2))

                if (delta_length == 0.0d0) cycle

                coords = coords_final + delta_coords * (handicap_rmsd + &
                        (threshold_rmsd - handicap_rmsd)*rand()) / delta_length

                call rmsd_dp(Natoms,coords_final,coords,1,candidate_U,x_center,y_center,min_rmsd)

                if (min_rmsd > handicap_rmsd) cycle
                
                handicap_rmsd = min_rmsd

                call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
                call Acceleration(vals,coords,gradient_var)

                gridpath3 = gridpath2
                local_frame_count = 0
                call addState_new(vals,coords,gradient_var)

                gradient_var = matmul(candidate_U,gradient_var)

!               write(filechannel2,FMT="(I3,3(1x,F13.9))") 1,min_rmsd, &
!                       sqrt(sum((gradient_labelled - gradient_var_labelled)**2)/Natoms),&
!                       sqrt(sum((gradient_labelled - gradient_var_labelled)**2)/Natoms)
                rmsd_x(1,Nsample) = min_rmsd
                rmsd_x_interpolated(1,Nsample) = min_rmsd
                rmsd_fx_interpolated(1,Nsample) = &
                        sqrt(sum((gradient - gradient_var)**2)/Natoms)
                rmsd_fx(1,Nsample) = rmsd_fx_interpolated(1,Nsample)
                rmsd_weights(1,1,Nsample) = 1.0d0

                write(filechannel2,FMT="(F15.11,1x,F15.11)") min_rmsd,rmsd_fx(1,Nsample)
                exit
        end do

        !Now lets add lots of bad points (but still in threshold)
        subcellsearch_max = (/ 9, 9 /)
        interpolation_flag = .true.

        do steps = 2, Ntest
                do
                        do n = 1, Natoms
                        do i = 1, 3
                                delta_coords(i,n) = rand() - 0.5d0
                        end do
                        end do
        
                        delta_length = sqrt(sum(delta_coords**2))
        
                        if (delta_length >= 1.0d0) cycle
                        if (delta_length == 0.0d0) cycle

                        coords = coords_final + delta_coords * &
                                (min_rmsd + (threshold_rmsd - min_rmsd)*rand())

                        call rmsd_dp(Natoms,coords_final,coords,1,&
                                     candidate_U,x_center,y_center,min_rmsd)

                        if (min_rmsd >= threshold_rmsd) cycle

                        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
                        call Acceleration(vals,coords,gradient_var)

                        gradient_var = matmul(candidate_U,gradient_var)

                        write(filechannel2,FMT="(F15.11,1x,F15.11)") min_rmsd,&
                                sqrt(sum((gradient - gradient_var)**2)/Natoms)

                        if (min_rmsd <= handicap_rmsd) cycle
                        exit
                end do

                rmsd_x(steps,Nsample) = min_rmsd

                interpolation_flag = .true.
                call checkState_new(vals,coords,approx_gradient,min_rmsd,&
                                filechannels,number_of_frames,order,neighbor_check)
                call addState_new(vals,coords,gradient_var)

                interpolation_flag = .true.
                min_rmsd = default_rmsd
                temp_frame_weights = 0.0d0
                temp_rmsd_weights = 0.0d0
                call checkState_new(vals,coords_final,approx_gradient,min_rmsd,&
                                filechannels,number_of_frames,order,neighbor_check)

!               write(filechannel2,FMT="(I3,3(1x,F13.9))") Ninterpolation, min_rmsd, &
!                       sqrt(sum((gradient_labelled - approx_gradient)**2)/Natoms),&
!                       sqrt(sum((gradient_labelled - gradient_var_labelled)**2)/Natoms)

                rmsd_x_interpolated(steps,Nsample) = min_rmsd
                rmsd_fx_interpolated(steps,Nsample) = &
                        sqrt(sum((gradient - approx_gradient)**2)/Natoms)
                rmsd_fx(steps,Nsample) = sqrt(sum((gradient - gradient_var)**2)/Natoms)
                rmsd_weights(1:steps,steps,Nsample) = temp_frame_weights
                rmsd_fx_weights(1:steps,steps,Nsample) = temp_rmsd_weights

                if (maxval(rmsd_fx(1:steps,Nsample),dim=1) < rmsd_fx_interpolated(steps,Nsample)) then
                        Nanomaly = Nanomaly + 1
                        write(Nanomaly_text,FMT="(I3)") Nanomaly
                        open(filechannel3,file=gridpath5//"anomaly"//&
                                trim(adjustl(Nanomaly_text))//".xyz")
                        write(filechannel3,FMT="(I4)") Natoms*2
                        write(filechannel3,FMT="(A)") "Fixed Frame then Offset Frame"
                        do n = 1, Natoms
                                write(filechannel3,FMT="(A,3(1x,F9.5))") "H", coords_final(:,n)
                        end do
                        do n = 1, Natoms
                                write(filechannel3,FMT="(A,3(1x,F9.5))") "H", coords(:,n)
                        end do
                        close(filechannel3)
                end if
        end do

        end do

        close(filechannel2)

        open(filechannel2,file=gridpath5//errorcheckfile)
        do steps = 1, Ntest
                selected_means(1) = sum(rmsd_x_interpolated(steps,:))/Nsamples
                selected_means(2) = sum(rmsd_x(steps,:))/Nsamples
                selected_means(3) = sum(rmsd_fx_interpolated(steps,:))/Nsamples
                selected_means(4) = sum(rmsd_fx(steps,:))/Nsamples

                selected_SDs(1) = sqrt(sum((rmsd_x_interpolated(steps,:)-&
                                            selected_means(1))**2)/Nsamples)
                selected_SDs(2) = sqrt(sum((rmsd_x(steps,:)-&
                                            selected_means(2))**2)/Nsamples)
                selected_SDs(3) = sqrt(sum((rmsd_fx_interpolated(steps,:)-&
                                            selected_means(3))**2)/Nsamples)
                selected_SDs(4) = sqrt(sum((rmsd_fx(steps,:)-&
                                            selected_means(4))**2)/Nsamples)

                write(filechannel2,FMT="(I3,8(1x,F15.11))") steps, &
                        selected_means(1), selected_SDs(1),&
                        selected_means(2), selected_SDs(2),&
                        selected_means(3), selected_SDs(3),&
                        selected_means(4), selected_SDs(4)

        end do
        close(filechannel2)

        open(filechannel2,file=gridpath5//"weight"//errorcheckfile)
        do steps = 1, Ntest
                selected_means(1) = sum(rmsd_x_interpolated(steps,:))/Nsamples
                selected_means(2) = sum(rmsd_x(steps,:))/Nsamples
                selected_means(3) = sum(rmsd_fx_interpolated(steps,:))/Nsamples
                selected_means(4) = sum(rmsd_fx(steps,:))/Nsamples

                do n = 1, Ntest
                        mean_weights(n) = sum(rmsd_weights(n,steps,:))/Nsamples
                        mean_rmsd_fx(n) = sum(rmsd_fx_weights(n,steps,:))/Nsamples
                end do

                write(filechannel2,FMT=*) steps, &
                        selected_means(3), mean_weights, mean_rmsd_fx

        end do
        close(filechannel2)

        do Nsample = 1, Nsamples
                rmsd_fx_interpolated(:,Nsample) = rmsd_fx_interpolated(:,Nsample) / rmsd_fx(:,Nsample)
                rmsd_fx(:,Nsample) = rmsd_fx(1,Nsample) / rmsd_fx(:,Nsample)
        end do

        open(filechannel2,file=gridpath5//"relative"//errorcheckfile)
        do steps = 1, Ntest
                selected_means(1) = sum(rmsd_x_interpolated(steps,:))/Nsamples
                selected_means(2) = sum(rmsd_x(steps,:))/Nsamples
                selected_means(3) = sum(rmsd_fx_interpolated(steps,:))/Nsamples
                selected_means(4) = sum(rmsd_fx(steps,:))/Nsamples

                selected_SDs(1) = sqrt(sum((rmsd_x_interpolated(steps,:)-&
                                            selected_means(1))**2)/Nsamples)
                selected_SDs(2) = sqrt(sum((rmsd_x(steps,:)-&
                                            selected_means(2))**2)/Nsamples)
                selected_SDs(3) = sqrt(sum((rmsd_fx_interpolated(steps,:)-&
                                            selected_means(3))**2)/Nsamples)
                selected_SDs(4) = sqrt(sum((rmsd_fx(steps,:)-&
                                            selected_means(4))**2)/Nsamples)

                write(filechannel2,FMT="(I3,8(1x,F15.11))") steps, &
                        selected_means(1), selected_SDs(1),&
                        selected_means(2), selected_SDs(2),&
                        selected_means(3), selected_SDs(3),&
                        selected_means(4), selected_SDs(4)

        end do
        close(filechannel2)

        if (gather_interpolation_flag) close(filechannel3)

!       call plotErrorCheck1(vals,Nsamples)

        coords = coords_final
        velocities = velocities_final

        end do

        deallocate(valsbuffer1,coordsbuffer1,gradientbuffer1,Ubuffer1,RMSDbuffer1,&
                   approximation_index)
        if (interpolation_flag) deallocate(acceptable_frame_mask,inputCLS)

        deallocate(rmsd_x_interpolated,rmsd_fx,&
                   rmsd_fx_interpolated,rmsd_x,&
                   rmsd_weights,rmsd_fx_weights)
        deallocate(temp_frame_weights,temp_rmsd_weights)

end subroutine errorCheck1

subroutine errorCheck2(filechannels)
        use FUNCTIONS
        use PARAMETERS
        use PHYSICS
        use VARIABLES
        use ANALYSIS
        use ls_rmsd_original
        implicit none

        !Coordinates, Velocities, and Variables
        real(dp), dimension(3,Natoms) :: coords,coords_labelled,velocities, delta_coords
        real(dp), dimension(3,Natoms) :: gradient,gradient_var,gradient_labelled,gradient_var_labelled
        real(dp), dimension(3,Natoms) :: approx_gradient,approx_gradient_prime
        real(dp), dimension(Nvar) :: vals
	real(dp), dimension(3,Natoms) :: coords_initial, velocities_initial 
	real(dp), dimension(3,Natoms) :: coords_final, velocities_final 
        integer :: bond_index1, bond_index2

        real(dp) :: handicap_rmsd
        real(dp), allocatable :: rmsd_x(:,:),rmsd_x_interpolated(:,:)
        real(dp), allocatable :: rmsd_fx(:,:),rmsd_fx_interpolated(:,:)
        real(dp), allocatable :: rmsd_weights(:,:,:)
        real(dp),dimension(4) :: selected_means,selected_SDs
        real(dp) :: delta_length
        integer :: Ntest,Nsamples,Nsample
        integer :: Nanomaly,Nanomaly_index
        character(3) :: Nanomaly_text
        real(dp), allocatable :: mean_weights(:),mean_rmsd_fxs(:)
        real(dp), allocatable :: inputCLS(:,:),outputCLS(:)
        real(dp), allocatable :: gradient_steps(:,:,:)
        real(dp), allocatable :: frame_weights(:)
        real(dp), allocatable :: restraints(:,:),restraint_values(:)
        real(dp), allocatable :: minimized_differences2(:,:)
        real(dp) :: error1,error2
        real(dp), allocatable :: rmsd_x2_interpolated(:,:)
        real(dp) :: LSn,LSx,LSy,LSx2,LSxy,LSdet
        real(dp), allocatable :: LSa1(:),LSa2(:),LSerror(:),convergence(:)
        integer, allocatable :: dropoff(:)
        real(dp) :: dropoff_cutoff,dropoff_mean,dropoff_SD
        real(dp) :: convergence_mean,convergence_SD

        !Various other variables
        real(dp) :: min_rmsd,min_rmsd_prime
        integer :: number_of_frames,order,neighbor_check
        character(9) :: vals_interpolation_text

        integer,dimension(1+Ngrid_total),intent(in) :: filechannels
        real(dp), dimension(3) :: x_center, y_center
        real(dp), allocatable :: g(:,:)
        real(dp),dimension(3,3) :: U
        real(dp),dimension(3,3) :: candidate_U

        !Incremental Integer
        integer :: i,n

        !Initialize the scene
        call InitialSetup3(coords,velocities)

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

        Nanomaly = 0
        Ntest = 20
        Nsamples = 50
        allocate(rmsd_x_interpolated(Ntest,Nsamples),rmsd_fx(Ntest,Nsamples),&
                 rmsd_fx_interpolated(Ntest,Nsamples),rmsd_x(Ntest,Nsamples),&
                 rmsd_weights(Ntest,Ntest,Nsamples),rmsd_x2_interpolated(Ntest,Nsamples))
        allocate(mean_weights(Ntest),mean_rmsd_fxs(Ntest))
        allocate(inputCLS(Ncoords+Ntest,Ntest),gradient_steps(3,Natoms,Ntest))
        rmsd_weights = 0.0d0

        open(filechannel3,file=gridpath5//"convergence"//errorcheckfile)

        do

        Nanomaly_index = 0
        call system("rm "//gridpath5//"weight1"//errorcheckfile)

        do steps = 1, 1000
                !Upate the coordinates with the velocities
                coords = coords + dt * velocities

                !Always calculate the variables before checking a frame or accelerating
                call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

                call Acceleration(vals,coords,gradient)

                !Update the velocities
                velocities = velocities + gradient
        end do

        if (any(vals > var_maxvar)) exit

        coords_final = coords
        velocities_final = velocities

        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
        call Acceleration(vals,coords,gradient)

        do n = 1, Natoms
                coords_labelled(:,n) = coords(:,BOND_LABELLING_DATA(n))
                gradient_labelled(:,n) = gradient(:,BOND_LABELLING_DATA(n))
        end do
        print *, "vals:", vals

        open(filechannel2,file=gridpath5//"new"//errorcheckfile)

        do Nsample = 1, Nsamples

        handicap_rmsd = threshold_rmsd* 0.1d0

        !Let's find one good point
        do
                do n = 1, Natoms
                do i = 1, 3
                        delta_coords(i,n) = rand() - 0.5d0
                end do
                end do

                delta_length = sqrt(sum(delta_coords**2))

                if (delta_length == 0.0d0) cycle

                coords = coords_final + delta_coords * (handicap_rmsd + &
                        (threshold_rmsd - handicap_rmsd)*rand()) / delta_length

                call rmsd_dp(Natoms,coords,coords_final,1,candidate_U,x_center,y_center,min_rmsd)

                if (min_rmsd > handicap_rmsd) cycle
                
                handicap_rmsd = min_rmsd

                call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
                call Acceleration(vals,coords,gradient_var)

                gradient_var = matmul(candidate_U,gradient_var)

                do i = 1, 3
                        coords(i,:) = &
                        coords(i,:) - x_center(i)
                end do

                coords = matmul(candidate_U,coords)

                do i = 1, 3
                        coords(i,:) = &
                        coords(i,:) + y_center(i)
                end do

                inputCLS(1:Ncoords,1) = reshape(coords - coords_final,(/Ncoords/))
                gradient_steps(:,:,1) = gradient_var

                inputCLS(Ncoords+1,:) = 0.0d0
                inputCLS(Ncoords+1,1) = alpha_ratio*maxval(abs(inputCLS(1:Ncoords,1)))**2

                rmsd_x(1,Nsample) = min_rmsd
                rmsd_x_interpolated(1,Nsample) = min_rmsd
                rmsd_fx_interpolated(1,Nsample) = &
                        sqrt(sum((gradient - gradient_var)**2)/Natoms)
                rmsd_fx(1,Nsample) = rmsd_fx_interpolated(1,Nsample)

                rmsd_x2_interpolated(1,Nsample) = inputCLS(Ncoords+1,1)/alpha_ratio

                rmsd_weights(1,1,Nsample) = 1.0d0

                write(filechannel2,FMT="(F15.11,1x,F15.11)") min_rmsd,rmsd_fx(1,Nsample)
                exit
        end do

        !Now lets add lots of bad points (but still in threshold)
        do steps = 2, Ntest
                do
                        do n = 1, Natoms
                        do i = 1, 3
                                delta_coords(i,n) = rand() - 0.5d0
                        end do
                        end do
        
                        delta_length = sqrt(sum(delta_coords**2))
        
                        if (delta_length >= 1.0d0) cycle
                        if (delta_length == 0.0d0) cycle

                        coords = coords_final + delta_coords * &
                                (min_rmsd + (threshold_rmsd - min_rmsd)*rand())

                        call rmsd_dp(Natoms,coords,coords_final,1,&
                                     candidate_U,x_center,y_center,min_rmsd)

                        if (min_rmsd >= threshold_rmsd) cycle

                        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
                        call Acceleration(vals,coords,gradient_var)

                        gradient_var = matmul(candidate_U,gradient_var)

                        write(filechannel2,FMT="(F15.11,1x,F15.11)") min_rmsd,&
                                sqrt(sum((gradient - gradient_var)**2)/Natoms)

                        if (min_rmsd <= handicap_rmsd) cycle
                        exit
                end do

                do i = 1, 3
                        coords(i,:) = &
                        coords(i,:) - x_center(i)
                end do

                coords = matmul(candidate_U,coords)

                do i = 1, 3
                        coords(i,:) = &
                        coords(i,:) + y_center(i)
                end do

                inputCLS(1:Ncoords,steps) = reshape(coords - coords_final,(/Ncoords/))
                gradient_steps(:,:,steps) = gradient_var

                allocate(frame_weights(steps),&
                         outputCLS(Ncoords+steps),&
                         restraints(1,steps),&
                         restraint_values(1),&
                         minimized_differences2(steps,1))
                
                inputCLS(Ncoords+steps,:) = 0.0d0
                inputCLS(Ncoords+steps,steps) = alpha_ratio * maxval(abs(inputCLS(1:Ncoords,steps)))**2
                
                restraints = 1.0d0
                restraint_values = 1.0d0
                outputCLS(1:Ncoords) = 0.0d0
                outputCLS(Ncoords+1:Ncoords+steps) = 0.0d0
                
                call CLS2(inputCLS(1:Ncoords+steps,&
                          1:steps),Ncoords+steps,steps,&
                          restraints,1,restraint_values,&
                          outputCLS,frame_weights(1:steps))

                approx_gradient = 0.0d0
                coords = 0.0d0
                error2 = 0.0d0
                do n = 1, steps
                        approx_gradient = approx_gradient + frame_weights(n) *&
                                gradient_steps(:,:,n)
                        coords = coords + reshape((frame_weights(n)*inputCLS(1:Ncoords,n)),&
                                (/3,Natoms/))
                        error2 = error2 + (frame_weights(n)*inputCLS(Ncoords+n,n))**2
                end do
                error1 = sqrt(sum(coords**2)/Natoms)
                error2 = sqrt(error2)

                rmsd_x(steps,Nsample) = min_rmsd
                rmsd_x_interpolated(steps,Nsample) = error1
                rmsd_fx_interpolated(steps,Nsample) = &
                        sqrt(sum((gradient - approx_gradient)**2)/Natoms)
                rmsd_fx(steps,Nsample) = sqrt(sum((gradient - gradient_var)**2)/Natoms)

                rmsd_x2_interpolated(steps,Nsample) = error2

                rmsd_weights(1:steps,steps,Nsample) = frame_weights(1:steps)

                deallocate(frame_weights,outputCLS,restraints,restraint_values,minimized_differences2)

                if (maxval(rmsd_fx(1:steps,Nsample),dim=1) < rmsd_fx_interpolated(steps,Nsample)) then
                        Nanomaly_index = Nsample
                        Nanomaly = Nanomaly + 1
                        write(Nanomaly_text,FMT="(I3)") Nanomaly
                        open(filechannel4,file=gridpath5//"anomaly"//&
                                trim(adjustl(Nanomaly_text))//".xyz")
                        write(filechannel4,FMT="(I4)") Natoms*2
                        write(filechannel4,FMT="(A)") "Fixed Frame then Offset Frame"
                        do n = 1, Natoms
                                write(filechannel4,FMT="(A,3(1x,F9.5))") "H", coords_final(:,n)
                        end do
                        do n = 1, Natoms
                                write(filechannel4,FMT="(A,3(1x,F9.5))") "H", coords(:,n)
                        end do
                        close(filechannel4)
                end if
        end do

        end do

        close(filechannel2)

        open(filechannel2,file=gridpath5//errorcheckfile)
        do steps = 1, Ntest
                selected_means(1) = sum(rmsd_x_interpolated(steps,:))/Nsamples
                selected_means(2) = sum(rmsd_x(steps,:))/Nsamples
                selected_means(3) = sum(rmsd_fx_interpolated(steps,:))/Nsamples
                selected_means(4) = sum(rmsd_fx(steps,:))/Nsamples

                selected_SDs(1) = sqrt(sum((rmsd_x_interpolated(steps,:)-&
                                            selected_means(1))**2)/Nsamples)
                selected_SDs(2) = sqrt(sum((rmsd_x(steps,:)-&
                                            selected_means(2))**2)/Nsamples)
                selected_SDs(3) = sqrt(sum((rmsd_fx_interpolated(steps,:)-&
                                            selected_means(3))**2)/Nsamples)
                selected_SDs(4) = sqrt(sum((rmsd_fx(steps,:)-&
                                            selected_means(4))**2)/Nsamples)

                write(filechannel2,FMT="(I3,8(1x,F15.11))") steps, &
                        selected_means(1), selected_SDs(1),&
                        selected_means(2), selected_SDs(2),&
                        selected_means(3), selected_SDs(3),&
                        selected_means(4), selected_SDs(4)

        end do
        close(filechannel2)

        open(filechannel2,file=gridpath5//"linear"//errorcheckfile)
        do steps = 1, Ntest
                do Nsample = 1, Nsamples
                write(filechannel2,FMT="(I3,3(1x,F15.11))") steps, &
                        rmsd_x_interpolated(steps,Nsample),&
                        rmsd_x2_interpolated(steps,Nsample),&
                        rmsd_fx_interpolated(steps,Nsample)
                end do
        end do
        close(filechannel2)

        allocate(LSa1(Ntest),LSa2(Ntest),LSerror(Ntest))
        allocate(convergence(Nsample),dropoff(Nsample))

        open(filechannel2,file=gridpath5//"dropoff"//errorcheckfile)
        do Nsample = 1, Nsamples

        convergence(Nsample) = 0.0d0
        do steps = 1, Ntest/2
                LSn = Ntest - steps + 1
                LSx = 0.0d0
                LSx2 = 0.0d0
                LSy = 0.0d0
                LSxy = 0.0d0
                do i = steps, Ntest
                        LSx  =  LSx + exp(-dble(i))
                        LSy  =  LSy + log10(rmsd_fx_interpolated(i,Nsample))
                        LSxy = LSxy + log10(rmsd_fx_interpolated(i,Nsample))*exp(-dble(i))
                        LSx2 = LSx2 + exp(-dble(i))**2
                end do
                LSdet = 1.0d0 / (LSx2*LSn - LSx**2)
                LSa1(steps) = (LSx2*LSy - LSx*LSxy)*LSdet
                LSa2(steps) = (LSn*LSxy - LSx*LSy)*LSdet
                LSerror(steps) = 0.0d0
                do i = steps, Ntest
                       LSerror(steps) = LSerror(steps) + &
                               (LSa1(steps) + LSa2(steps)*exp(-dble(i)) - &
                               log10(rmsd_fx_interpolated(i,Nsample)))**2
                end do
                LSerror(steps) = LSerror(steps) / LSn
                if (LSerror(steps) == 0) then
                        convergence(Nsample) = LSa1(steps)
                        exit
                end if

!               print *, "   step:", steps
!               print *, "log(fx):", log10(rmsd_fx_interpolated(steps,Nsample))
!               print *, "   LSa1:", LSa1(steps)
!               print *, "   LSa2:", LSa2(steps)
!               print *, "LSerror:", LSerror(steps)
                convergence(Nsample) = convergence(Nsample) + &
                        LSa1(steps) / LSerror(steps)
        end do

        if (minval(LSerror(1:Ntest/2)) == 0) then
        else
                convergence(Nsample) = convergence(Nsample) /&
                        sum(LSerror(1:Ntest/2)**(-1))
        end if
        dropoff_cutoff = sqrt(rmsd_fx_interpolated(1,Nsample)*&
                                10.0d0**(convergence(Nsample)))
        dropoff(Nsample) = minloc(rmsd_fx_interpolated(1:Ntest,Nsample),1,&
                rmsd_fx_interpolated(1:Ntest,Nsample)>dropoff_cutoff) + 1

        write(filechannel2,FMT="(I3,1x,F15.11,1x,I3)") Nsample, &
                convergence(Nsample), dropoff(Nsample)
        end do
        close(filechannel2)

        dropoff_mean = sum(dropoff)*1.0d0/Nsamples
        convergence_mean = sum(convergence)/Nsamples

        dropoff_SD = sqrt(sum((dropoff-dropoff_mean)**2)/Nsamples)
        convergence_SD = sqrt(sum((convergence-convergence_mean)**2)/Nsamples)

        deallocate(LSa1,LSa2,LSerror,convergence,dropoff)

        open(filechannel2,file=gridpath5//"heatmapline"//errorcheckfile)
        do steps = 1, Ntest
                write(filechannel2,FMT="(I3,3(1x,F15.11))") steps, &
                        rmsd_x_interpolated(steps,min(Nsamples/2,69)),&
                        rmsd_x2_interpolated(steps,min(Nsamples/2,69)),&
                        rmsd_fx_interpolated(steps,min(Nsamples/2,69))
        end do
        close(filechannel2)

        open(filechannel2,file=gridpath5//"ratio"//errorcheckfile)
        do steps = 1, Ntest
                do Nsample = 1, Nsamples
                write(filechannel2,FMT="(I3,3(1x,F15.11))") steps, &
                        rmsd_fx_interpolated(steps,Nsample)/&
                        rmsd_fx(1,Nsample),&
                        rmsd_fx_interpolated(steps,Nsample)/&
                        rmsd_fx(steps,Nsample),&
                        rmsd_fx_interpolated(steps,Nsample)/&
                        maxval(rmsd_fx(1:steps,Nsample),dim=1)
                end do
        end do
        close(filechannel2)

        if (Nanomaly_index > 0) then
        open(filechannel2,file=gridpath5//"weight1"//errorcheckfile)
        do steps = 1, Ntest
                selected_means(1) = sum(rmsd_x_interpolated(steps,:))/Nsamples
                selected_means(2) = sum(rmsd_x(steps,:))/Nsamples
                selected_means(3) = sum(rmsd_fx_interpolated(steps,:))/Nsamples
                selected_means(4) = sum(rmsd_fx(steps,:))/Nsamples

                do n = 1, Ntest
!                       mean_weights(n) = sum(abs(rmsd_weights(n,steps,:)))/Nsamples
!                       mean_rmsd_fxs(n) = sum(rmsd_fx(n,:))/Nsamples
                        mean_weights(n) = abs(rmsd_weights(n,steps,Nanomaly_index))
                        mean_rmsd_fxs(n) = rmsd_fx(n,Nanomaly_index)
                end do

!               write(filechannel2,FMT=*) steps, &
!                       selected_means(3), mean_weights, mean_rmsd_fxs
                write(filechannel2,FMT=*) steps,&
                        rmsd_fx_interpolated(steps,Nanomaly_index),&
                        mean_weights, mean_rmsd_fxs

        end do
        close(filechannel2)
        end if

!       open(filechannel2,file=gridpath0//"weight2"//errorcheckfile)
!       do steps = 1, Ntest
!               selected_means(1) = sum(rmsd_x_interpolated(steps,:))/Nsamples
!               selected_means(2) = sum(rmsd_x(steps,:))/Nsamples
!               selected_means(3) = sum(rmsd_fx_interpolated(steps,:))/Nsamples
!               selected_means(4) = sum(rmsd_fx(steps,:))/Nsamples

!               do n = 1, Ntest
!                       mean_weights(n) = sum(rmsd_weights(steps,n,:))/Nsamples
!               end do

!               write(filechannel2,FMT=*) steps, &
!                       selected_means(3), mean_weights, rmsd_fx(:,Nsample)

!       end do
!       close(filechannel2)

        do Nsample = 1, Nsamples
                rmsd_fx_interpolated(:,Nsample) = rmsd_fx_interpolated(:,Nsample) / rmsd_fx(:,Nsample)
                rmsd_fx(:,Nsample) = rmsd_fx(1,Nsample) / rmsd_fx(:,Nsample)
        end do

        open(filechannel2,file=gridpath5//"relative"//errorcheckfile)
        do steps = 1, Ntest
                selected_means(1) = sum(rmsd_x_interpolated(steps,:))/Nsamples
                selected_means(2) = sum(rmsd_x(steps,:))/Nsamples
                selected_means(3) = sum(rmsd_fx_interpolated(steps,:))/Nsamples
                selected_means(4) = sum(rmsd_fx(steps,:))/Nsamples

                selected_SDs(1) = sqrt(sum((rmsd_x_interpolated(steps,:)-&
                                            selected_means(1))**2)/Nsamples)
                selected_SDs(2) = sqrt(sum((rmsd_x(steps,:)-&
                                            selected_means(2))**2)/Nsamples)
                selected_SDs(3) = sqrt(sum((rmsd_fx_interpolated(steps,:)-&
                                            selected_means(3))**2)/Nsamples)
                selected_SDs(4) = sqrt(sum((rmsd_fx(steps,:)-&
                                            selected_means(4))**2)/Nsamples)

                write(filechannel2,FMT="(I3,8(1x,F15.11))") steps, &
                        selected_means(1), selected_SDs(1),&
                        selected_means(2), selected_SDs(2),&
                        selected_means(3), selected_SDs(3),&
                        selected_means(4), selected_SDs(4)

        end do
        close(filechannel2)

        write(filechannel3,FMT=*) vals,dropoff_mean,dropoff_SD,&
                convergence_mean,convergence_SD

        call plotErrorCheck1(vals,Nsamples,&
                dropoff_mean,dropoff_SD,convergence_mean,convergence_SD)

        coords = coords_final
        velocities = velocities_final

        end do

        close(filechannel3)

        deallocate(rmsd_x_interpolated,rmsd_fx,&
                   rmsd_fx_interpolated,rmsd_x,&
                   rmsd_weights)
        deallocate(inputCLS,gradient_steps)

        call plotFinalErrorCheck1(Nsamples)

end subroutine errorCheck2

subroutine plotErrorCheck1(vals,Nsamples,&
        dropoff_mean,dropoff_SD,convergence_mean,convergence_SD)
        use FUNCTIONS
        use PARAMETERS
        use PHYSICS
        use VARIABLES
        use ANALYSIS
        implicit none

        real(dp), dimension(Nvar), intent(in) :: vals
        integer, intent(in) ::Nsamples
        real(dp), intent(in) :: dropoff_mean,dropoff_SD
        real(dp), intent(in) :: convergence_mean,convergence_SD
        real(dp) :: first_rmsd, first_rmsd_sd, max_rmsd, min_rmsd, max_weight
        real(dp) :: rmsd_x,rmsd_x_interpolated,rmsd_fx,rmsd_fx_interpolated
        real(dp) :: dum1,dum2,dum3,dum4
        real(dp) :: interpolated_rmsd
        real(dp) :: rmsd1,rmsd2,rmsd3
        real(dp) :: drmsd1,drmsd2,drmsd3
        real(dp) :: max_rmsd1,max_rmsd2,max_rmsd3
        real(dp) :: min_rmsd1,min_rmsd2,min_rmsd3
        real(dp), allocatable :: frame_weights(:), frame_rmsds(:)
        integer, allocatable :: rmsd1_bins(:),rmsd2_bins(:),rmsd3_bins(:)
        real(dp),allocatable :: RMSDheatmap(:,:), A(:,:), b(:)
        real(dp),allocatable :: RMSDheatmap_coeff(:,:), DE(:,:),coeff_root(:)
        real(dp),allocatable :: independent_variables(:,:),dependent_variables(:)
        real(dp) :: coeff_dist, E1,E2, dum2_previous, levenberg1
        integer :: Nheatmap,Ngif
        integer :: Nbins
        logical :: file_exists
        integer :: Niteration,Niterations
        real(dp) :: coeff_threshold,coeff_lambda
        integer :: n, m, iostate, Ninterpolation
        character(3) :: ntext

        max_rmsd = 0.0d0
        min_rmsd = default_rmsd
        Ninterpolation = 0

        open(filechannel2,file=gridpath5//errorcheckfile)
        do
                read(filechannel2,iostat=iostate,FMT=*) n, rmsd_x,dum1,&
                        rmsd_x_interpolated,dum2,&
                        rmsd_fx_interpolated,dum3,rmsd_fx,dum4
                if (iostate /= 0) exit

                if ((n == 1).and.(max_rmsd == 0.0d0)) then
                        first_rmsd = rmsd_fx
                        first_rmsd_sd = dum4
                end if

                max_rmsd = max(max_rmsd,rmsd_fx_interpolated+dum3,rmsd_fx+dum4)
                min_rmsd = min(min_rmsd,rmsd_fx_interpolated,rmsd_fx)
                Ninterpolation = n
        end do
        close(filechannel2)

        open(gnuplotchannel,file=gridpath5//gnuplotfile)
        write(gnuplotchannel,*) "set term pngcairo size 1200,1200"
        write(gnuplotchannel,FMT="(A,2F7.4,A)") &
                'set output "'//gridpath4, vals, '.png"'
        write(gnuplotchannel,*) 'set title "Error Convergence As More Points Interpolate"'
        write(gnuplotchannel,*) 'set xlabel "Ninterpolation"'
        write(gnuplotchannel,*) 'set ylabel "RMSD Between Interpolated and Target Gradient"'
        write(gnuplotchannel,FMT="(A,I5,A)") 'set label 1 "N = ', Nsamples, '" at screen 0.1,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,",",F7.4,A)') 'set label 2 "(var1,var2) = ',vals, '" at screen 0.2,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 3 "Threshhold = ',threshold_rmsd, ' A" at screen 0.5,0.925'
        write(gnuplotchannel,FMT='(A,F7.2,A,F7.4,A)') &
                'set label 4 "log(Convergence) = ',convergence_mean,&
                             '     Dropoff = ',dropoff_mean,'" at screen 0.1,0.900'
        write(gnuplotchannel,FMT='(A,F9.4,A)') 'set label 5 "AlphaRatio = ',alpha_ratio, '" at screen 0.50,0.900'
        write(gnuplotchannel,*) 'y0 = ', first_rmsd
        write(gnuplotchannel,*) 'y0_err = ', first_rmsd_sd
        write(gnuplotchannel,*) 'ymax = ', max_rmsd
        write(gnuplotchannel,*) 'ymin = ', min_rmsd
        write(gnuplotchannel,*) 'ydelta = log10(ymax/ymin)/8'
        write(gnuplotchannel,*) 'set yrange [ymin*10**(-ydelta):ymax*10**(ydelta)]'
        write(gnuplotchannel,*) 'xmax = ', Ninterpolation
        write(gnuplotchannel,*) 'set xrange [1:xmax]'
        write(gnuplotchannel,*) 'set logscale y'
        write(gnuplotchannel,*) 'set ytics ('//&
                                     '"1e-11" .00000000001, '//&
                                     '"5e-11" .00000000005, '//&
                                      '"1e-10" .0000000001, '//&
                                      '"5e-10" .0000000005, '//&
                                        '"1e-9" .000000001, '//&
                                        '"5e-9" .000000005, '//&
                                         '"1e-8" .00000001, '//&
                                         '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1.0"       1, '//&
                                   ')'
        write(gnuplotchannel,*) 'plot "'//gridpath5//errorcheckfile//'" u 1:6:7 w yerr lc "green" notitle,\'
        write(gnuplotchannel,*) '     "'//gridpath5//errorcheckfile//'" u 1:6 w lines lc "green" t "Interpolation",\'
        write(gnuplotchannel,*) '     "'//gridpath5//errorcheckfile//'" u 1:8:9 w yerr lc "red" notitle,\'
        write(gnuplotchannel,*) '     "'//gridpath5//errorcheckfile//'" u 1:8 w lines lc "red" t "Accept Current",\'
!       write(gnuplotchannel,*) '     "'//gridpath5//errorcheckfile//'" u 1:(y0):(y0_err) w yerr lc "blue" notitle,\'
        write(gnuplotchannel,*) '     "'//gridpath5//errorcheckfile//'" u 1:(y0+y0_err) w lines lc "blue" notitle,\'
        write(gnuplotchannel,*) '     "'//gridpath5//errorcheckfile//'" u 1:(y0-y0_err) w lines lc "blue" t "Accept Best"'
        close(gnuplotchannel)

        call system("gnuplot < "//gridpath5//gnuplotfile)

        max_rmsd = 0.0d0
        min_rmsd = 1.0d9
        Ninterpolation = 0

        open(filechannel2,file=gridpath5//"relative"//errorcheckfile)
        do
                read(filechannel2,iostat=iostate,FMT=*) n, rmsd_x,dum1,&
                        rmsd_x_interpolated,dum2,&
                        rmsd_fx_interpolated,dum3,rmsd_fx,dum4
                if (iostate /= 0) exit

                max_rmsd = max(max_rmsd,rmsd_fx_interpolated+dum3,rmsd_fx+dum4)
                min_rmsd = min(min_rmsd,rmsd_fx_interpolated,rmsd_fx)
                Ninterpolation = n
        end do
        close(filechannel2)

        open(gnuplotchannel,file=gridpath5//gnuplotfile)
        write(gnuplotchannel,*) "set term pngcairo size 1200,1200"
        write(gnuplotchannel,FMT="(A,2F7.4,A)") &
                'set output "'//gridpath5//"relative_", vals, '.png"'
        write(gnuplotchannel,*) 'set title "Error Convergence As More Points Interpolate"'
        write(gnuplotchannel,*) 'set xlabel "Ninterpolation"'
        write(gnuplotchannel,*) 'set ylabel "Relative RMSD Between Interpolated and Target Gradient"'
        write(gnuplotchannel,FMT="(A,I5,A)") 'set label 1 "N = ', Nsamples, '" at screen 0.1,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,",",F7.4,A)') 'set label 2 "(var1,var2) = ',vals, '" at screen 0.2,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 3 "Threshhold = ',threshold_rmsd, ' A" at screen 0.5,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 4 "AlphaRatio = ',alpha_ratio, '" at screen 0.50,0.895'
        write(gnuplotchannel,*) 'ymax = ', max_rmsd
        write(gnuplotchannel,*) 'ymin = ', min_rmsd
        write(gnuplotchannel,*) 'ydelta = log10(ymax/ymin)/8'
        write(gnuplotchannel,*) 'set yrange [ymin*10**(-ydelta):ymax*10**(ydelta)]'
        write(gnuplotchannel,*) 'xmax = ', Ninterpolation
        write(gnuplotchannel,*) 'set xrange [1:xmax]'
        write(gnuplotchannel,*) 'set logscale y'
        write(gnuplotchannel,*) 'set ytics ('//&
                                        '"1e-9" .000000001, '//&
                                        '"5e-9" .000000005, '//&
                                         '"1e-8" .00000001, '//&
                                         '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1.0"       1, '//&
                                           ' "5.0"       5, '//&
                                           ' "1e1"      10, '//&
                                           ' "5e1"      50, '//&
                                           ' "1e2"     100, '//&
                                           ' "5e2"     500, '//&
                                           ' "1e3"    1000, '//&
                                           ' "5e3"    5000, '//&
                                           ' "1e4"   10000, '//&
                                           ' "5e4"   50000, '//&
                                   ')'
        write(gnuplotchannel,*) 'plot "'//gridpath5//"relative"//errorcheckfile//&
                '" u 1:6:7 w yerr lc "green" notitle,\'
        write(gnuplotchannel,*) '     "'//gridpath5//"relative"//errorcheckfile//&
                '" u 1:6 w lines lc "green" t "Interpolation",\'
        write(gnuplotchannel,*) '     "'//gridpath5//"relative"//errorcheckfile//&
                '" u 1:(1.0) w lines lc "red" t "Accept Current",\'
        close(gnuplotchannel)

        call system("gnuplot < "//gridpath5//gnuplotfile)

        inquire(file=gridpath5//"weight1"//errorcheckfile,exist=file_exists)
        if (file_exists) then

        max_rmsd = 0.0d0
        min_rmsd = 1.0d9
        max_weight = 0.0d0

        allocate(frame_weights(Ninterpolation),frame_rmsds(Ninterpolation))

        open(filechannel2,file=gridpath5//"weight1"//errorcheckfile)
        do
                read(filechannel2,iostat=iostate,FMT=*) n, &
                        interpolated_rmsd,frame_weights,frame_rmsds
                if (iostate /= 0) exit

                max_weight = max(max_weight,maxval(frame_weights,dim=1))
                max_rmsd = max(max_rmsd,interpolated_rmsd)
                min_rmsd = min(min_rmsd,interpolated_rmsd)
        end do
        close(filechannel2)

!       max_rmsd = max(max_rmsd, maxval(frame_rmsds,dim=1))
!       min_rmsd = min(min_rmsd, minval(frame_rmsds,dim=1))

        open(gnuplotchannel,file=gridpath5//gnuplotfile)
        write(gnuplotchannel,*) "set term pngcairo size 1200,1200"
        write(gnuplotchannel,FMT="(A,2F7.4,A)") &
                'set output "'//gridpath5//"weight1_", vals, '.png"'
        write(gnuplotchannel,*) 'set title "Error Convergence As More Points Interpolate"'
        write(gnuplotchannel,*) 'set multiplot'
        write(gnuplotchannel,*) 'set size 1.0, 1.0'
        write(gnuplotchannel,*) 'set origin 0.0, 0.0'
        write(gnuplotchannel,*) 'set pm3d map'
        write(gnuplotchannel,*) 'unset key'
        write(gnuplotchannel,*) 'unset colorbox'
        write(gnuplotchannel,*) 'set xlabel "Ninterpolation"'
        write(gnuplotchannel,*) 'set ylabel "Absolute Value of Weight"'
        write(gnuplotchannel,*) 'set cblabel "RMSD Between Gradient of Frame and Target"'
        write(gnuplotchannel,FMT='(A,F7.4,",",F7.4,A)') 'set label 1 "(var1,var2) = ',vals, '" at screen 0.2,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 2 "Threshhold = ',threshold_rmsd, ' A" at screen 0.5,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 3 "AlphaRatio = ',alpha_ratio, '" at screen 0.50,0.900'
        write(gnuplotchannel,*) 'cbmax = ', max_rmsd
        write(gnuplotchannel,*) 'cbmin = ', min_rmsd
        write(gnuplotchannel,*) 'wmax = ', max_weight
        write(gnuplotchannel,*) 'cbdelta = log10(cbmax/cbmin)/8'
        write(gnuplotchannel,*) 'set cbrange [cbmin*10**(-cbdelta):cbmax*10**(cbdelta)]'
        write(gnuplotchannel,*) 'set palette defined (cbmin "blue", cbmax "red")'
        write(gnuplotchannel,*) 'Ninterpolation = ', Ninterpolation
        write(gnuplotchannel,*) 'set yrange [0:1.1*wmax]'
        write(gnuplotchannel,*) 'set xrange [1:Ninterpolation]'
        write(gnuplotchannel,*) 'set grid xtics front lw 2'
        write(gnuplotchannel,*) 'set grid'
        write(gnuplotchannel,*) 'set logscale cb'
!       write(gnuplotchannel,*) 'set cbtics ('//&
!                                       '"1e-9" .000000001, '//&
!                                       '"5e-9" .000000005, '//&
!                                        '"1e-8" .00000001, '//&
!                                        '"5e-8" .00000005, '//&
!                                         '"1e-7" .0000001, '//&
!                                         '"5e-7" .0000005, '//&
!                                          '"1e-6" .000001, '//&
!                                          '"5e-6" .000005, '//&
!                                          '"1e-5"  .00001, '//&
!                                          '"5e-5"  .00005, '//&
!                                          '"1e-4"   .0001, '//&
!                                          '"5e-4"   .0005, '//&
!                                          '"1e-3"    .001, '//&
!                                          '"5e-3"    .005, '//&
!                                          '"1e-2"     .01, '//&
!                                          '"5e-2"     .05, '//&
!                                          '"1e-1"      .1, '//&
!                                          '"5e-1"      .5, '//&
!                                          ' "1.0"       1, '//&
!                                          ' "5.0"       5, '//&
!                                          ' "1e1"      10, '//&
!                                          ' "5e1"      50, '//&
!                                          ' "1e2"     100, '//&
!                                          ' "5e2"     500, '//&
!                                          ' "1e3"    1000, '//&
!                                          ' "5e3"    5000, '//&
!                                          ' "1e4"   10000, '//&
!                                          ' "5e4"   50000, '//&
!                                  ')'
        write(gnuplotchannel,*) 'set palette defined ('//&
                                '0 "green", '//&
                                '0.800 "cyan", '//&
                                '0.850 "blue", '//&
                                '0.900 "violet", '//&
                                '0.950 "red", '//&
                                '1 "orange"'//&
                                ')'
        write(gnuplotchannel,*) 'plot "'//gridpath5//"weight1"//errorcheckfile//&
                '" u 1:3:Ninterpolation+3 w lines lw 2 palette,\'
        do n = 2, Ninterpolation-1
        write(gnuplotchannel,FMT="(A,I3,A,I3,A)") &
                '     "'//gridpath5//"weight1"//errorcheckfile//&
                '" u 1:',n+2,':',Ninterpolation+n+2,' w lines lw 2 palette,\'
        end do
        write(gnuplotchannel,*) '     "'//gridpath5//"weight1"//errorcheckfile//&
                '" u 1:Ninterpolation+2:Ninterpolation*2+2 w lines lw 2 palette'

        write(gnuplotchannel,*) 'set size 0.4, 0.4'
        write(gnuplotchannel,*) 'set origin 0.45, 0.5'
        write(gnuplotchannel,*) "set lmargin at screen 0.55"
        write(gnuplotchannel,*) "set rmargin at screen 0.95"
        write(gnuplotchannel,*) "set bmargin at screen 0.5"
        write(gnuplotchannel,*) "set tmargin at screen 0.9"

!       write(gnuplotchannel,*) 'set size 0.4, 0.4'
!       write(gnuplotchannel,*) 'set origin 0.45, 0.5'
!       write(gnuplotchannel,*) 'set isosample 100,100'
!       write(gnuplotchannel,*) 'unset title'
!       write(gnuplotchannel,*) 'unset xlabel'
!       write(gnuplotchannel,*) 'unset ylabel'
!       write(gnuplotchannel,*) 'set logscale y'
!       write(gnuplotchannel,*) 'unset colorbox'
!       write(gnuplotchannel,*) 'set yrange [cbmin:cbmax]'
!       write(gnuplotchannel,*) 'splot y t ""'

        write(gnuplotchannel,*) 'unset title'
        write(gnuplotchannel,*) 'unset xlabel'
        write(gnuplotchannel,*) 'unset ylabel'
        write(gnuplotchannel,*) 'unset xtics'
        write(gnuplotchannel,*) 'set logscale y'
        write(gnuplotchannel,*) 'set ylabel "Error of Interpolation" offset -0.5,0'
        write(gnuplotchannel,*) 'set ytics ('//&
                                      '"1e-10" .0000000001, '//&
                                      '"5e-10" .0000000005, '//&
                                        '"1e-9" .000000001, '//&
                                        '"5e-9" .000000005, '//&
                                         '"1e-8" .00000001, '//&
                                         '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1.0"       1, '//&
                                   ')'
        write(gnuplotchannel,*) 'set yrange [cbmax*10**(cbdelta):cbmin*10**(-cbdelta)]'
        write(gnuplotchannel,*) 'set isosample 100,100'
        write(gnuplotchannel,*) 'set pm3d map'
        write(gnuplotchannel,*) 'splot y'
        write(gnuplotchannel,*) 'set yrange [cbmin*10**(-cbdelta):cbmax*10**(cbdelta)]'
        write(gnuplotchannel,*) 'set xtics'
        write(gnuplotchannel,*) 'unset ytics'
        write(gnuplotchannel,*) 'unset ylabel'
        write(gnuplotchannel,*) 'plot "'//gridpath5//"weight1"//errorcheckfile//&
                '" u 1:2:2 w l lc "black"'
        close(gnuplotchannel)

        call system("gnuplot < "//gridpath5//gnuplotfile)

        end if

        open(gnuplotchannel,file=gridpath5//gnuplotfile)
        write(gnuplotchannel,*) "set term pngcairo size 1200,1200"
        write(gnuplotchannel,FMT="(A,2F7.4,A)") &
                'set output "'//gridpath5//"new", vals, '.png"'
        write(gnuplotchannel,*) 'set title "Error Convergence As Frames Get Closer"'
        write(gnuplotchannel,*) 'unset key'
        write(gnuplotchannel,*) 'set xlabel "RMSD Between Nearby and Target Frame"'
        write(gnuplotchannel,*) 'set ylabel "Error Between Nearby and Target Gradient"'
        write(gnuplotchannel,FMT="(A,I5,A)") 'set label 1 "N = ', Nsamples, '" at screen 0.1,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,",",F7.4,A)') 'set label 2 "(var1,var2) = ',vals, '" at screen 0.2,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 3 "Threshhold = ',threshold_rmsd, ' A" at screen 0.5,0.925'
        write(gnuplotchannel,*) 'set logscale xy'
        write(gnuplotchannel,*) 'set xtics ('//&
                                     '"1e-11" .00000000001, '//&
                                     '"5e-11" .00000000005, '//&
                                      '"1e-10" .0000000001, '//&
                                      '"5e-10" .0000000005, '//&
                                        '"1e-9" .000000001, '//&
                                        '"5e-9" .000000005, '//&
                                         '"1e-8" .00000001, '//&
                                         '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                           '"5e-2"     .05, '//&
                                           '"1e-1"      .1, '//&
                                           '"5e-1"      .5, '//&
                                           ' "1.0"       1, '//&
                                   ')'
        write(gnuplotchannel,*) 'set ytics ('//&
                                     '"1e-11" .00000000001, '//&
                                     '"5e-11" .00000000005, '//&
                                      '"1e-10" .0000000001, '//&
                                      '"5e-10" .0000000005, '//&
                                        '"1e-9" .000000001, '//&
                                        '"5e-9" .000000005, '//&
                                         '"1e-8" .00000001, '//&
                                         '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                   ')'
        write(gnuplotchannel,*) 'plot "'//gridpath5//"new"//errorcheckfile//'" u 1:2 w points'
        close(gnuplotchannel)

        call system("gnuplot < "//gridpath5//gnuplotfile)

        max_rmsd1 = 0.0d1
        min_rmsd1 = 1.0d9
        max_rmsd2 = 0.0d1
        min_rmsd2 = 1.0d9
        max_rmsd3 = 0.0d1
        min_rmsd3 = 1.0d9

        open(filechannel2,file=gridpath5//"linear"//errorcheckfile)
        do
                read(filechannel2,iostat=iostate,FMT=*) n, &
                        rmsd1,rmsd2,rmsd3
                if (iostate /= 0) exit

                max_rmsd1 = max(max_rmsd1,rmsd1)
                min_rmsd1 = min(min_rmsd1,rmsd1)
                max_rmsd2 = max(max_rmsd2,rmsd2)
                min_rmsd2 = min(min_rmsd2,rmsd2)
                max_rmsd3 = max(max_rmsd3,rmsd3)
                min_rmsd3 = min(min_rmsd3,rmsd3)
                Ninterpolation = n
        end do
        close(filechannel2)

        if (min_rmsd1 == 0.0d0) min_rmsd1 = 1.0d-11
        if (min_rmsd2 == 0.0d0) min_rmsd2 = 1.0d-11
        if (min_rmsd3 == 0.0d0) min_rmsd3 = 1.0d-11

        open(gnuplotchannel,file=gridpath5//gnuplotfile)
        write(gnuplotchannel,*) 'set term pngcairo size 1200,2400'
        write(gnuplotchannel,FMT="(A,2F7.4,A)") &
                'set output "'//gridpath5//"linear", vals, '.png"'
        write(gnuplotchannel,*) 'set title "Error Comparison of a Frame and Gradient with its Interpolation"'
        write(gnuplotchannel,*) 'set multiplot layout 2,1'
        write(gnuplotchannel,*) 'set pm3d map'
        write(gnuplotchannel,*) 'unset key'
        write(gnuplotchannel,FMT="(A,I5,A)") 'set label 1 "N = ', Nsamples, '" at screen 0.1,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,",",F7.4,A)') 'set label 2 "(var1,var2) = ',vals, '" at screen 0.2,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 3 "Threshhold = ',threshold_rmsd, ' A" at screen 0.5,0.925'
        write(gnuplotchannel,*) 'min_x = ', min_rmsd1
        write(gnuplotchannel,*) 'max_x = ', max_rmsd1
        write(gnuplotchannel,*) 'min_y = ', min_rmsd3
        write(gnuplotchannel,*) 'max_y = ', max_rmsd3
        write(gnuplotchannel,*) 'min_cx = ', 1
        write(gnuplotchannel,*) 'max_cx = ', Ninterpolation
        write(gnuplotchannel,*) 'set logscale x'
        write(gnuplotchannel,*) 'set logscale y'
        write(gnuplotchannel,*) 'set xrange [0.9*min_x:1.1*max_x]'
        write(gnuplotchannel,*) 'set yrange [0.9*min_y:1.1*max_y]'
        write(gnuplotchannel,*) 'set cbrange [min_cx:max_cx]'
        write(gnuplotchannel,*) 'set palette defined (min_cx "blue", max_cx "red")'
        write(gnuplotchannel,*) 'set xlabel "RMSD Variable 1"'
        write(gnuplotchannel,*) 'set ylabel "Error Between the Output Gradient and the Interpolation"'
        write(gnuplotchannel,*) 'set cblabel "Number of Frames Used to Interpolate"'

        write(gnuplotchannel,*) 'plot "'//gridpath5//"linear"//errorcheckfile//'" u '//&
                               '($2>0?$2:1/0):4:1 w p lw 6 palette'
        write(gnuplotchannel,*) 'set logscale x'
        write(gnuplotchannel,*) 'set logscale y'
        write(gnuplotchannel,*) 'set xlabel "RMSD Variable 2"'
        write(gnuplotchannel,*) 'set ylabel "Error Between the Output Gradient and the Interpolation"'
        write(gnuplotchannel,*) 'min_x = ', min_rmsd2
        write(gnuplotchannel,*) 'max_x = ', max_rmsd2
        write(gnuplotchannel,*) 'set xrange [0.9*min_x:1.1*max_x]'
        write(gnuplotchannel,*) 'plot "'//gridpath5//"linear"//errorcheckfile//'" u '//&
                               '($3>0?$3:1/0):4:1 w p lw 6 palette'
        close(gnuplotchannel)
        
        call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

        Nbins = 50
        allocate(RMSDheatmap(Nbins,Nbins))
        RMSDheatmap = 1.0d-12

        drmsd1 = log10(max_rmsd1/min_rmsd1)/Nbins
        drmsd2 = log10(max_rmsd2/min_rmsd2)/Nbins
        drmsd3 = log10(max_rmsd3/min_rmsd3)/Nbins
!       drmsd1 = (max_rmsd1-min_rmsd1)/Nbins
!       drmsd2 = (max_rmsd2-min_rmsd2)/Nbins
!       drmsd3 = (max_rmsd3-min_rmsd3)/Nbins
        
        open(filechannel2,file=gridpath5//"linear"//errorcheckfile)
        do
                read(filechannel2,FMT=*,iostat=iostate) Ninterpolation,&
                        rmsd1,rmsd2,rmsd3
                if (iostate /= 0) exit
        
                n = nint(log10(rmsd1/min_rmsd1)/drmsd1)
!               n = nint((rmsd1-min_rmsd1)/drmsd1)
                if (n < 1) n = 1
                m = nint(log10(rmsd2/min_rmsd2)/drmsd2)
!               m = nint((rmsd2-min_rmsd2)/drmsd2)
                if (m < 1) m = 1
        
                RMSDheatmap(min(n,Nbins),min(m,Nbins)) = max(&
                        RMSDheatmap(min(n,Nbins),min(m,Nbins)),rmsd3)
        end do
        close(filechannel2)

!       max_rmsd1 = log10(max_rmsd1)
!       max_rmsd2 = log10(max_rmsd2)
!       max_rmsd3 = log10(max_rmsd3)
!       min_rmsd1 = log10(min_rmsd1)
!       min_rmsd2 = log10(min_rmsd2)
!       min_rmsd3 = log10(min_rmsd3)

!       open(filechannel2,file=gridpath5//"heatmap_linear"//errorcheckfile)
!       do n = 1, Nbins
!       do m = 1, Nbins
!               write(filechannel2,FMT=*) (n-0.5)*drmsd1 + min_rmsd1,&
!                       (n-0.5)*drmsd2 + min_rmsd2, RMSDheatmap(n,m)
!       end do
!       end do
!       close(filechannel2)

        Nheatmap = 0
        do n = 1, Nbins
                do m = 1, Nbins
                        if (RMSDheatmap(n,m) > 1.0d-12) then
                                 Nheatmap = Nheatmap + 1
                        end if
                end do
        end do
        
!       allocate(A(Nheatmap,3),b(Nheatmap))
!       Nheatmap = 0
!       do n = 1, Nbins
!               do m = 1, Nbins
!                       if (RMSDheatmap(n,m) > 1.0d-12) then
!                                Nheatmap = Nheatmap + 1
!                                A(Nheatmap,:) = (/ (n-0.5)*drmsd1,&
!                                        (m-0.5)*drmsd2,1.0d0/)
!                                b(Nheatmap) = log10(RMSDheatmap(n,m))
!                                b(Nheatmap) = RMSDheatmap(n,m)
!                       end if
!               end do
!       end do
!       
!       call LS(A,Nheatmap,3,b,RMSDheatmap_coeff)

        allocate(independent_variables(Nheatmap,2),dependent_variables(Nheatmap))
        Nheatmap = 0
        do n = 1, Nbins
                do m = 1, Nbins
                        if (RMSDheatmap(n,m) > 1.0d-12) then
                                Nheatmap = Nheatmap + 1
!                               independent_variables(Nheatmap,:) = &
!                                       (/(n-0.5)*drmsd1,(m-0.5)*drmsd2/)
                                independent_variables(Nheatmap,:) = &
                                        (/min_rmsd1*10.0d0**((n-0.5)*drmsd1),&
                                          min_rmsd2*10.0d0**((m-0.5)*drmsd2)/)
                                dependent_variables(Nheatmap) = RMSDheatmap(n,m)
                        end if
                end do
        end do

        Niterations = 100
        coeff_threshold = 1.0d-12
        coeff_lambda = 1.0d0
        allocate(RMSDheatmap_coeff(Niterations+1,3),DE(Niterations+1,3),coeff_root(Nheatmap),&
                 A(Nheatmap+3,3),b(Nheatmap+3))

        RMSDheatmap_coeff(1,:) = (/1.0d-2,1.0d-2,0.0d0/)
        dum2 = 0.0d0
        do n = 1, Nheatmap
                dum1 = (RMSDheatmap_coeff(1,1)*independent_variables(n,1) +&
                        RMSDheatmap_coeff(1,2)*independent_variables(n,2) +&
                        RMSDheatmap_coeff(1,3))
                dum2 = dum2 + (log10(dependent_variables(n) / dum1))**2
        end do

        dum2_previous = dum2
        print *, "(k,a1,a2,a3):", 1, RMSDheatmap_coeff(1,:)
        print *, "       (k,E):", 1, dum2

        do Niteration = 1, Niterations
                DE(Niteration,1:3) = 0.0d0
                do n = 1, Nheatmap
                        dum1 = (RMSDheatmap_coeff(Niteration,1)*independent_variables(n,1) +&
                                RMSDheatmap_coeff(Niteration,2)*independent_variables(n,2) +&
                                RMSDheatmap_coeff(Niteration,3))**(-1)
                        DE(Niteration,1:3) = DE(Niteration,1:3) - &
                                2*(log10(dependent_variables(n) * dum1) * dum1) *&
                                (/ independent_variables(n,1),&
                                   independent_variables(n,2), 1.0d0 /)
                end do
                print *, "(k,DE1,DE2,DE3):", Niteration, DE(Niteration,:)
                do n = 1, Nheatmap
                        dum1 = (RMSDheatmap_coeff(Niteration,1)*independent_variables(n,1) +&
                                RMSDheatmap_coeff(Niteration,2)*independent_variables(n,2) +&
                                RMSDheatmap_coeff(Niteration,3))
                        coeff_root(n) = - dum1 / ( &
                                DE(Niteration,1)*independent_variables(n,1) +&
                                DE(Niteration,2)*independent_variables(n,2) +&
                                DE(Niteration,3) )
!                       print *, "root:", coeff_root(n)
                end do
                coeff_dist = maxval(coeff_root,dim=1)
                levenberg1 = 1.0d3
!               dum4 = coeff_dist * 1.1d0
                dum4 = 0.0d0
                print *, "lambda_min:", coeff_dist
                do
                        E1 = 0.0d0
                        E2 = 0.0d0
                        do n = 1, Nheatmap
                                dum1 = (RMSDheatmap_coeff(Niteration,1)*independent_variables(n,1) +&
                                        RMSDheatmap_coeff(Niteration,2)*independent_variables(n,2) +&
                                        RMSDheatmap_coeff(Niteration,3))
                                dum2 = (DE(Niteration,1)*independent_variables(n,1) +&
                                        DE(Niteration,2)*independent_variables(n,2) +&
                                        DE(Niteration,3))
                                dum3 = (dum1 + dum4*dum2)**(-1)
                                E1 = E1 - 2 * (log10(dependent_variables(n)*dum3)*&
                                        dum2*dum3)
                                E2 = E2 + 2*(1.0d0 + log10(dependent_variables(n)*dum3))*&
                                        (dum2*dum3)**2
                        end do

                        do
                                dum1 = E1 / (E2 * (1.0d0 + levenberg1))
!                               print *, "E2:", E2
!                               print *, "E1/E2:", dum1

                                if ((sign(1.0d0,dum4-dum1-coeff_dist) /= sign(1.0d0,dum4-coeff_dist))) then
                                        levenberg1 = levenberg1*10
                                        cycle
                                end if
!                               print *, "lambda:", dum4-dum1
        
                                RMSDheatmap_coeff(Niteration+1,:) = RMSDheatmap_coeff(Niteration,:) +&
                                                                    (dum4-dum1) * DE(Niteration,:)
        
                                dum2 = 0.0d0
                                do n = 1, Nheatmap
                                        dum1 = (RMSDheatmap_coeff(Niteration+1,1)*independent_variables(n,1) +&
                                                RMSDheatmap_coeff(Niteration+1,2)*independent_variables(n,2) +&
                                                RMSDheatmap_coeff(Niteration+1,3))
                                        dum2 = dum2 + (log10(dependent_variables(n) / dum1))**2
                                end do
!                               print *, "(k,a1,a2,a3):", Niteration, RMSDheatmap_coeff(Niteration+1,:)
!                               print *, "       (k,E):", Niteration, dum2

!                               call sleep(1)

                                if (dum2 > dum2_previous) then
                                        levenberg1 = levenberg1*10
                                        cycle
                                end if

!                               print *, "Successful levenberg"

!                               dum1 = E1 / (E2 + levenberg1)
                                dum1 = E1 / (E2 * (1.0d0 + levenberg1))
                                dum4 = dum4 - dum1
                                levenberg1 = levenberg1*0.1
                                dum2_previous = dum2
                                exit
                        end do

                        if ((abs(dum1/dum4) < 1.0d-12)) exit

!                       call sleep(1)
                end do
!               RMSDheatmap_coeff(Niteration+1,:) = RMSDheatmap_coeff(Niteration,:) +&
!                                                   dum4 * DE(Niteration,:)
!               do n = 1, Nheatmap
!                       dum1 = (RMSDheatmap_coeff(Niteration,1)*independent_variables(n,1) +&
!                               RMSDheatmap_coeff(Niteration,2)*independent_variables(n,2) +&
!                               RMSDheatmap_coeff(Niteration,3))**(-1)
!                       A(n,1:2) = independent_variables(n,:) * dum1
!                       A(n,3) = dum1
!                       A(n,:) = A(n,:)
!                       b(n) =  (log10(dependent_variables(n) * dum1) - &
!                                (A(n,1)*RMSDheatmap_coeff(Niteration,1)+&
!                                 A(n,2)*RMSDheatmap_coeff(Niteration,2)+&
!                                 A(n,3)*RMSDheatmap_coeff(Niteration,3)))
!               end do
!               do n = 1, 3
!                       A(Nheatmap+n,n) = coeff_lambda
!                       b(Nheatmap+n) = -RMSDheatmap_coeff(Niteration,n)
!               end do

!               call LS(A,Nheatmap,3,b,RMSDheatmap_coeff(Niteration+1,:))
!               dum2 = 0.0d0
!               do n = 1, Nheatmap
!                       dum1 = (RMSDheatmap_coeff(Niteration+1,1)*independent_variables(n,1) +&
!                               RMSDheatmap_coeff(Niteration+1,2)*independent_variables(n,2) +&
!                               RMSDheatmap_coeff(Niteration+1,3))
!                       dum2 = dum2 + (log10(dependent_variables(n) / dum1))**2
!               end do

                print *, "(k,a1,a2,a3):", Niteration+1, RMSDheatmap_coeff(Niteration+1,:)
                print *, "       (k,E):", Niteration+1, dum2

                exit
                if (maxval(abs(RMSDheatmap_coeff(Niteration+1,:)-&
                               RMSDheatmap_coeff(Niteration,:)),dim=1) < &
                               coeff_threshold) exit

!               call sleep(1)
!               call sleep(5)
        end do

        DE(Niteration+1,1:3) = 0.0d0
        do n = 1, Nheatmap
                dum1 = (RMSDheatmap_coeff(Niteration+1,1)*independent_variables(n,1) +&
                        RMSDheatmap_coeff(Niteration+1,2)*independent_variables(n,2) +&
                        RMSDheatmap_coeff(Niteration+1,3))**(-1)
                DE(Niteration+1,1:3) = DE(Niteration+1,1:3) - &
                        2*(log10(dependent_variables(n) * dum1) * dum1) *&
                        (/ independent_variables(n,1),&
                           independent_variables(n,2), 1.0d0 /)
        end do
        print *, "(k,DE1,DE2,DE3):", Niteration+1, DE(Niteration+1,:)

        print *, "Converged!"

        open(filechannel2,file=gridpath5//"heatmap_linear"//errorcheckfile)
        do n = 1, Nbins
                do m = 1, Nbins
                        dum1 = (n-0.5)*drmsd1
                        dum2 = (m-0.5)*drmsd2
                        write(filechannel2,FMT=*) dum1,&
!                               dum2,RMSDheatmap(n,m),&
                                dum2,log10(RMSDheatmap(n,m)),&
!                               10.0d0**(max(dum1*RMSDheatmap_coeff(Niteration+1,1)+&
!                                            dum2*RMSDheatmap_coeff(Niteration+1,2)+&
!                                                 RMSDheatmap_coeff(Niteration+1,3),-12.0d0))
                                max(log10(RMSDheatmap_coeff(Niteration+1,1)*min_rmsd1*10.0d0**(dum1)+&
                                          RMSDheatmap_coeff(Niteration+1,2)*min_rmsd2*10.0d0**(dum2)+&
                                          RMSDheatmap_coeff(Niteration+1,3)),-12.0d0)
                end do
                write(filechannel2,*) ""
        end do
        close(filechannel2)

        deallocate(A,b,RMSDheatmap)

        open(gnuplotchannel,file=gridpath5//gnuplotfile)
        write(gnuplotchannel,*) 'set term pngcairo size 2400,1200'
        write(gnuplotchannel,FMT="(A,2F7.4,A)") &
                'set output "'//gridpath5//"heatmap_linear", vals, '.png"'
        write(gnuplotchannel,*) 'set multiplot layout 1,2'
        write(gnuplotchannel,*) 'set title "Interpolation Error Comparison by RMSD Variables 1 and 2"'
        write(gnuplotchannel,*) 'set pm3d map'
        write(gnuplotchannel,*) 'unset key'
        write(gnuplotchannel,FMT="(A,I5,A)") 'set label 1 "N = ', Nsamples, '" at screen 0.05,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,",",F7.4,A)') 'set label 2 "(var1,var2) = ',vals, '" at screen 0.1,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 3 "Threshhold = ',threshold_rmsd, ' A" at screen 0.25,0.925'
        write(gnuplotchannel,*) 'xmin = ', log10(min_rmsd1)
        write(gnuplotchannel,*) 'xmax = ', log10(max_rmsd1)
        write(gnuplotchannel,*) 'ymin = ', log10(min_rmsd2)
        write(gnuplotchannel,*) 'ymax = ', log10(max_rmsd2)
        write(gnuplotchannel,*) 'set xrange [0:',drmsd1*Nbins,']'
        write(gnuplotchannel,*) 'set yrange [0:',drmsd2*Nbins,']'
        !write(gnuplotchannel,*) 'xmin = ', min_rmsd_z
        !write(gnuplotchannel,*) 'xmax = ', max_rmsd_z
        !write(gnuplotchannel,*) 'ymin = ', min_rmsd_x
        !write(gnuplotchannel,*) 'ymax = ', max_rmsd_x
        write(gnuplotchannel,*) 'min_cx = ', log10(min_rmsd3)
        write(gnuplotchannel,*) 'max_cx = ', log10(max_rmsd3)
        !write(gnuplotchannel,*) 'max_cx = ', max(max_rmsd_fx,1.0e-7)
        write(gnuplotchannel,*) 'set cbrange [min_cx:max_cx]'
        write(gnuplotchannel,*) 'set palette defined ('//&
                                'log10(.0000001) "white", '//&
                                'log10(.0000001) "yellow", '//&
                                'log10(.0000005) "yellow", '//&
                                'log10(.000001) "green", '//&
                                'log10(.000005) "cyan", '//&
                                'log10(.00001) "blue", '//&
                                'log10(.00005) "magenta", '//&
                                'log10(.0001) "red"'//&
                                ')'
!       write(gnuplotchannel,*) 'set palette defined ('//&
!                               '.0000001 "white", '//&
!                               '.0000005 "yellow", '//&
!                               '.000001 "green", '//&
!                               '.000005 "cyan", '//&
!                               '.00001 "blue", '//&
!                               '.00005 "magenta", '//&
!                               '.0001 "red"'//&
!                               ')'
        write(gnuplotchannel,*) 'set xlabel "RMSD Variable 1"'
        write(gnuplotchannel,*) 'set ylabel "RMSD Variable 2"'
        write(gnuplotchannel,*) 'set cblabel "Maximum Error Encountered Between Interpolated and Target Gradient"'
        write(gnuplotchannel,*) 'set xtics ('//&
                                                '"1e-9" log10(.000000001)-xmin, '//&
!                                               '"5e-9" log10(.000000005)-xmin, '//&
                                                 '"1e-8" log10(.00000001)-xmin, '//&
!                                                '"5e-8" log10(.00000005)-xmin, '//&
                                                  '"1e-7" log10(.0000001)-xmin, '//&
!                                                 '"5e-7" log10(.0000005)-xmin, '//&
                                                   '"1e-6" log10(.000001)-xmin, '//&
!                                                  '"5e-6" log10(.000005)-xmin, '//&
                                                   '"1e-5"  log10(.00001)-xmin, '//&
!                                                  '"5e-5"  log10(.00005)-xmin, '//&
                                                   '"1e-4"   log10(.0001)-xmin, '//&
!                                                  '"5e-4"   log10(.0005)-xmin, '//&
                                                   '"1e-3"    log10(.001)-xmin, '//&
!                                                  '"5e-3"    log10(.005)-xmin, '//&
                                                   '"1e-2"     log10(.01)-xmin, '//&
!                                                  '"5e-2"     log10(.05)-xmin, '//&
                                                   '"1e-1"      log10(.1)-xmin, '//&
!                                                  '"5e-1"      log10(.5)-xmin, '//&
                                                   ' "1.0"       log10(1)-xmin, '//&
                                           ')'
        write(gnuplotchannel,*) 'set ytics ('//&
                                                '"1e-9" log10(.000000001)-ymin, '//&
!                                               '"5e-9" log10(.000000005)-xmin, '//&
                                                 '"1e-8" log10(.00000001)-ymin, '//&
!                                                '"5e-8" log10(.00000005)-ymin, '//&
                                                  '"1e-7" log10(.0000001)-ymin, '//&
!                                                 '"5e-7" log10(.0000005)-ymin, '//&
                                                   '"1e-6" log10(.000001)-ymin, '//&
!                                                  '"5e-6" log10(.000005)-ymin, '//&
                                                   '"1e-5"  log10(.00001)-ymin, '//&
!                                                  '"5e-5"  log10(.00005)-ymin, '//&
                                                   '"1e-4"   log10(.0001)-ymin, '//&
!                                                  '"5e-4"   log10(.0005)-ymin, '//&
                                                   '"1e-3"    log10(.001)-ymin, '//&
!                                                  '"5e-3"    log10(.005)-ymin, '//&
                                                   '"1e-2"     log10(.01)-ymin, '//&
!                                                  '"5e-2"     log10(.05)-ymin, '//&
                                                   '"1e-1"      log10(.1)-ymin, '//&
!                                                  '"5e-1"      log10(.5)-ymin, '//&
                                                   ' "1.0"       log10(1)-ymin, '//&
                                           ')'
        write(gnuplotchannel,*) 'set cbtics ('//&
                                              '"1e-10" log10(.0000000001), '//&
!                                             '"5e-10" log10(.0000000005), '//&
                                                '"1e-9" log10(.000000001), '//&
!                                               '"5e-9" log10(.000000005), '//&
                                                 '"1e-8" log10(.00000001), '//&
!                                                '"5e-8" log10(.00000005), '//&
                                                  '"1e-7" log10(.0000001), '//&
!                                                 '"5e-7" log10(.0000005), '//&
                                                   '"1e-6" log10(.000001), '//&
!                                                  '"5e-6" log10(.000005), '//&
                                                   '"1e-5"  log10(.00001), '//&
!                                                  '"5e-5"  log10(.00005), '//&
                                                   '"1e-4"   log10(.0001), '//&
!                                                  '"5e-4"   log10(.0005), '//&
                                                   '"1e-3"    log10(.001), '//&
!                                                  '"5e-3"    log10(.005), '//&
                                                   '"1e-2"     log10(.01), '//&
!                                                  '"5e-2"     log10(.05), '//&
                                                   '"1e-1"      log10(.1), '//&
!                                                  '"5e-1"      log10(.5), '//&
                                                   ' "1.0"       log10(1), '//&
                                           ')'
        write(gnuplotchannel,*) 'splot "'//gridpath5//'heatmap_linear'//errorcheckfile//&
                '" u 1:2:3 w image palette'
        write(gnuplotchannel,*) 'set title "Linear Best Fit of Interpolation Error"'
        write(gnuplotchannel,*) 'unset colorbox'
        write(gnuplotchannel,FMT="(A,D12.5,A,D12.5,A,D12.5,A)") &
                'set label 4 "f(x,y) = ', RMSDheatmap_coeff(Niteration+1,1),'x + ',&
                                          RMSDheatmap_coeff(Niteration+1,2),'y + ',&
                                          RMSDheatmap_coeff(Niteration+1,3),'" at screen 0.75,0.850 front'
        write(gnuplotchannel,*) 'set xlabel "RMSD Variable 1 (x)"'
        write(gnuplotchannel,*) 'set ylabel "RMSD Variable 2 (y)"'
        write(gnuplotchannel,*) 'splot "'//gridpath5//'heatmap_linear'//errorcheckfile//&
                '" u 1:2:4 w image palette'
        close(gnuplotchannel)

        call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)
!
!        call system("rm -r "//gridpath5//"png/")
!        call system("mkdir "//gridpath5//"png")
!        Ngif = Ninterpolation
!        do n = 1, Ngif
!        write(ntext,FMT="(I0.3)") n
!
!        open(gnuplotchannel,file=gridpath5//gnuplotfile)
!        write(gnuplotchannel,*) 'set term pngcairo size 2400,1200'
!        write(gnuplotchannel,FMT="(A)") &
!                'set output "'//gridpath5//"png/"//ntext//'.png"'
!        write(gnuplotchannel,*) 'set multiplot'
!        write(gnuplotchannel,*) 'set size 0.5, 1.0'
!        write(gnuplotchannel,*) 'set origin 0.0, 0.0'
!        write(gnuplotchannel,*) 'set title "Interpolation Error Comparison by RMSD Variables 1 and 2"'
!        write(gnuplotchannel,*) 'set pm3d map'
!        write(gnuplotchannel,*) 'unset key'
!        write(gnuplotchannel,FMT="(A,I5,A)") 'set label 1 "N = ', Nsamples, '" at screen 0.055,0.925'
!        write(gnuplotchannel,FMT='(A,F7.4,",",F7.4,A)') 'set label 2 "(var1,var2) = ',vals, '" at screen 0.1,0.925'
!        write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 3 "Threshhold = ',threshold_rmsd, ' A" at screen 0.25,0.925'
!        write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 4 "AlphaRatio = ',alpha_ratio, '" at screen 0.25,0.900'
!        write(gnuplotchannel,*) 'xmin = ', log10(min_rmsd1)
!        write(gnuplotchannel,*) 'xmax = ', log10(max_rmsd1)
!        write(gnuplotchannel,*) 'ymin = ', log10(min_rmsd2)
!        write(gnuplotchannel,*) 'ymax = ', log10(max_rmsd2)
!        write(gnuplotchannel,*) 'set xrange [0:',drmsd1*Nbins,']'
!        write(gnuplotchannel,*) 'set yrange [0:',drmsd2*Nbins,']'
!        !write(gnuplotchannel,*) 'xmin = ', min_rmsd_z
!        !write(gnuplotchannel,*) 'xmax = ', max_rmsd_z
!        !write(gnuplotchannel,*) 'ymin = ', min_rmsd_x
!        !write(gnuplotchannel,*) 'ymax = ', max_rmsd_x
!        write(gnuplotchannel,*) 'min_cx = ', log10(min_rmsd3)
!        write(gnuplotchannel,*) 'max_cx = ', log10(max_rmsd3)
!        !write(gnuplotchannel,*) 'max_cx = ', max(max_rmsd_fx,1.0e-7)
!        write(gnuplotchannel,*) 'set cbrange [min_cx:max_cx]'
!        write(gnuplotchannel,*) 'set palette defined ('//&
!                                'log10(.0000001) "white", '//&
!                                'log10(.0000001) "yellow", '//&
!                                'log10(.0000005) "yellow", '//&
!                                'log10(.000001) "green", '//&
!                                'log10(.000005) "cyan", '//&
!                                'log10(.00001) "blue", '//&
!                                'log10(.00005) "magenta", '//&
!                                'log10(.0001) "red"'//&
!                                ')'
!!       write(gnuplotchannel,*) 'set palette defined ('//&
!!                               '.0000001 "white", '//&
!!                               '.0000005 "yellow", '//&
!!                               '.000001 "green", '//&
!!                               '.000005 "cyan", '//&
!!                               '.00001 "blue", '//&
!!                               '.00005 "magenta", '//&
!!                               '.0001 "red"'//&
!!                               ')'
!        write(gnuplotchannel,*) 'set xlabel "RMSD Variable 1"'
!        write(gnuplotchannel,*) 'set ylabel "RMSD Variable 2"'
!        write(gnuplotchannel,*) 'set cblabel "Maximum Error Encountered Between Interpolated and Target Gradient"'
!        write(gnuplotchannel,*) 'set xtics ('//&
!                                                '"1e-9" log10(.000000001)-xmin, '//&
!!                                               '"5e-9" log10(.000000005)-xmin, '//&
!                                                 '"1e-8" log10(.00000001)-xmin, '//&
!!                                                '"5e-8" log10(.00000005)-xmin, '//&
!                                                  '"1e-7" log10(.0000001)-xmin, '//&
!!                                                 '"5e-7" log10(.0000005)-xmin, '//&
!                                                   '"1e-6" log10(.000001)-xmin, '//&
!!                                                  '"5e-6" log10(.000005)-xmin, '//&
!                                                   '"1e-5"  log10(.00001)-xmin, '//&
!!                                                  '"5e-5"  log10(.00005)-xmin, '//&
!                                                   '"1e-4"   log10(.0001)-xmin, '//&
!!                                                  '"5e-4"   log10(.0005)-xmin, '//&
!                                                   '"1e-3"    log10(.001)-xmin, '//&
!!                                                  '"5e-3"    log10(.005)-xmin, '//&
!                                                   '"1e-2"     log10(.01)-xmin, '//&
!!                                                  '"5e-2"     log10(.05)-xmin, '//&
!                                                   '"1e-1"      log10(.1)-xmin, '//&
!!                                                  '"5e-1"      log10(.5)-xmin, '//&
!                                                   ' "1.0"       log10(1)-xmin, '//&
!                                           ')'
!        write(gnuplotchannel,*) 'set ytics ('//&
!                                                '"1e-9" log10(.000000001)-ymin, '//&
!!                                               '"5e-9" log10(.000000005)-xmin, '//&
!                                                 '"1e-8" log10(.00000001)-ymin, '//&
!!                                                '"5e-8" log10(.00000005)-ymin, '//&
!                                                  '"1e-7" log10(.0000001)-ymin, '//&
!!                                                 '"5e-7" log10(.0000005)-ymin, '//&
!                                                   '"1e-6" log10(.000001)-ymin, '//&
!!                                                  '"5e-6" log10(.000005)-ymin, '//&
!                                                   '"1e-5"  log10(.00001)-ymin, '//&
!!                                                  '"5e-5"  log10(.00005)-ymin, '//&
!                                                   '"1e-4"   log10(.0001)-ymin, '//&
!!                                                  '"5e-4"   log10(.0005)-ymin, '//&
!                                                   '"1e-3"    log10(.001)-ymin, '//&
!!                                                  '"5e-3"    log10(.005)-ymin, '//&
!                                                   '"1e-2"     log10(.01)-ymin, '//&
!!                                                  '"5e-2"     log10(.05)-ymin, '//&
!                                                   '"1e-1"      log10(.1)-ymin, '//&
!!                                                  '"5e-1"      log10(.5)-ymin, '//&
!                                                   ' "1.0"       log10(1)-ymin, '//&
!                                           ')'
!        write(gnuplotchannel,*) 'set cbtics ('//&
!                                              '"1e-10" log10(.0000000001), '//&
!!                                             '"5e-10" log10(.0000000005), '//&
!                                                '"1e-9" log10(.000000001), '//&
!!                                               '"5e-9" log10(.000000005), '//&
!                                                 '"1e-8" log10(.00000001), '//&
!!                                                '"5e-8" log10(.00000005), '//&
!                                                  '"1e-7" log10(.0000001), '//&
!!                                                 '"5e-7" log10(.0000005), '//&
!                                                   '"1e-6" log10(.000001), '//&
!!                                                  '"5e-6" log10(.000005), '//&
!                                                   '"1e-5"  log10(.00001), '//&
!!                                                  '"5e-5"  log10(.00005), '//&
!                                                   '"1e-4"   log10(.0001), '//&
!!                                                  '"5e-4"   log10(.0005), '//&
!                                                   '"1e-3"    log10(.001), '//&
!!                                                  '"5e-3"    log10(.005), '//&
!                                                   '"1e-2"     log10(.01), '//&
!!                                                  '"5e-2"     log10(.05), '//&
!                                                   '"1e-1"      log10(.1), '//&
!!                                                  '"5e-1"      log10(.5), '//&
!                                                   ' "1.0"       log10(1), '//&
!                                           ')'
!
!        write(gnuplotchannel,*) "set lmargin at screen 0.05"
!        write(gnuplotchannel,*) "set rmargin at screen 0.45"
!        write(gnuplotchannel,*) "set bmargin at screen 0.10"
!        write(gnuplotchannel,*) "set tmargin at screen 0.95"
!        write(gnuplotchannel,*) 'splot "'//gridpath5//'heatmap_linear'//errorcheckfile//&
!                '" u 1:2:3 w image palette'
!
!        write(gnuplotchannel,*) "unset title"
!        write(gnuplotchannel,*) "unset xlabel"
!        write(gnuplotchannel,*) "unset ylabel"
!        write(gnuplotchannel,*) 'unset xtics'
!        write(gnuplotchannel,*) 'unset ytics'
!        write(gnuplotchannel,*) 'plot "<(sed -n ''1,'//trim(adjustl(ntext))//&
!                'p'' '//gridpath5//"heatmapline"//errorcheckfile//&
!                ')" u (log10($2)-xmin):(log10($3)-ymin) w l lw 3 lc "black"'
!        write(gnuplotchannel,FMT="(A,F15.11,A)") 'f(x) = (ymax-ymin)/2 + x/(', alpha_ratio, '*(10**(xmin/ymin)))'
!        write(gnuplotchannel,*) 'plot f(x) w l lw 2 lc "black"'
!
!        write(gnuplotchannel,*) 'set xrange [1:',Ninterpolation,']'
!        write(gnuplotchannel,*) 'set logscale y'
!        write(gnuplotchannel,*) "unset lmargin"
!        write(gnuplotchannel,*) "unset rmargin"
!        write(gnuplotchannel,*) "unset bmargin"
!        write(gnuplotchannel,*) "unset tmargin"
!        write(gnuplotchannel,*) 'set xtics'
!        write(gnuplotchannel,*) 'set ytics'
!        write(gnuplotchannel,*) 'set xlabel "Ninterpolation"'
!
!        write(gnuplotchannel,*) 'set size 0.5, 0.5'
!        write(gnuplotchannel,*) 'set origin 0.5, 0.0'
!        write(gnuplotchannel,*) 'set title "Title 2"'
!        write(gnuplotchannel,*) 'set yrange [',min_rmsd1,':',max_rmsd1,']'
!        write(gnuplotchannel,*) 'set ylabel "RMSD Variable 1"'
!        write(gnuplotchannel,*) 'plot "<(sed -n ''1,'//trim(adjustl(ntext))//&
!                'p'' '//gridpath5//"heatmapline"//errorcheckfile//&
!                ')" u 1:2 w l lc "black" lw 3'
!
!        write(gnuplotchannel,*) 'set size 0.5, 0.5'
!        write(gnuplotchannel,*) 'set origin 0.5, 0.5'
!        write(gnuplotchannel,*) 'set title "Title 3"'
!        write(gnuplotchannel,*) 'set yrange [',min_rmsd2,':',max_rmsd2,']'
!        write(gnuplotchannel,*) 'set ylabel "RMSD Variable 2"'
!        write(gnuplotchannel,*) 'plot "<(sed -n ''1,'//trim(adjustl(ntext))//&
!                'p'' '//gridpath5//"heatmapline"//errorcheckfile//&
!                ')" u 1:3 w l lc "black" lw 3'
!        close(gnuplotchannel)
!
!        call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)
!        end do
!
!        print *, "Started making the gif!"
!        call system("convert -background white -alpha remove -layers OptimizePlus "//&
!                "-delay 20 -loop 0 "//gridpath0//"png/*.png "//&
!                gridpath5//"heatmap_trajectory.gif")
!        print *, "Finished making the gif!"
!        call sleep(5)
!
        max_rmsd1 = 0.0d1
        min_rmsd1 = 1.0d9
        max_rmsd2 = 0.0d1
        min_rmsd2 = 1.0d9
        max_rmsd3 = 0.0d1
        min_rmsd3 = 1.0d9
        
        Nbins = 50
        allocate(rmsd1_bins(Nbins),rmsd2_bins(Nbins),rmsd3_bins(Nbins))

        open(filechannel2,file=gridpath5//"ratio"//errorcheckfile)
        do
                read(filechannel2,iostat=iostate,FMT=*) n, &
                        rmsd1,rmsd2,rmsd3
                if (iostate /= 0) exit

                max_rmsd1 = max(max_rmsd,rmsd1)
                min_rmsd1 = min(min_rmsd,rmsd1)
                max_rmsd2 = max(max_rmsd,rmsd2)
                min_rmsd2 = min(min_rmsd,rmsd2)
                max_rmsd3 = max(max_rmsd,rmsd3)
                min_rmsd3 = min(min_rmsd,rmsd3)
                Ninterpolation = n
        end do
        close(filechannel2)

        if (min_rmsd1 == 0.0d0) min_rmsd1 = 1.0d-11
        if (min_rmsd2 == 0.0d0) min_rmsd2 = 1.0d-11
        if (min_rmsd3 == 0.0d0) min_rmsd3 = 1.0d-11

        drmsd1 = log10(max_rmsd1/min_rmsd1)/Nbins
        drmsd2 = log10(max_rmsd2/min_rmsd2)/Nbins
        drmsd3 = log10(max_rmsd3/min_rmsd3)/Nbins

        rmsd1_bins = 0
        rmsd2_bins = 0
        rmsd3_bins = 0
        open(filechannel2,file=gridpath5//"ratio"//errorcheckfile)
        do
                read(filechannel2,iostat=iostate,FMT=*) n, &
                        rmsd1,rmsd2,rmsd3
                if (iostate /= 0) exit

                n = max(0,nint(log10(rmsd1/min_rmsd1)/drmsd1 - 0.5))
                rmsd1_bins(min(n+1,Nbins)) = rmsd1_bins(min(n+1,Nbins)) + 1
                n = max(0,nint(log10(rmsd2/min_rmsd2)/drmsd2 - 0.5))
                rmsd2_bins(min(n+1,Nbins)) = rmsd2_bins(min(n+1,Nbins)) + 1
                n = max(0,nint(log10(rmsd3/min_rmsd3)/drmsd3 - 0.5))
                rmsd3_bins(min(n+1,Nbins)) = rmsd3_bins(min(n+1,Nbins)) + 1
        end do
        close(filechannel2)

        max_rmsd1 = log10(max_rmsd1)
        max_rmsd2 = log10(max_rmsd2)
        max_rmsd3 = log10(max_rmsd3)
        min_rmsd1 = log10(min_rmsd1)
        min_rmsd2 = log10(min_rmsd2)
        min_rmsd3 = log10(min_rmsd3)

        open(filechannel2,file=gridpath5//"binned_ratio"//errorcheckfile)
        do n = 1, Nbins
                write(filechannel2,FMT=*) &
                        min_rmsd1 + ((n-0.5)*drmsd1), rmsd1_bins(n),&
                        min_rmsd2 + ((n-0.5)*drmsd2), rmsd2_bins(n),&
                        min_rmsd3 + ((n-0.5)*drmsd3), rmsd3_bins(n)
        end do
        close(filechannel2)
        deallocate(rmsd1_bins,rmsd2_bins,rmsd3_bins)

        open(gnuplotchannel,file=gridpath5//gnuplotfile)
        write(gnuplotchannel,*) 'set term pngcairo size 1200,2400'
        write(gnuplotchannel,FMT="(A,2F7.4,A)") &
                'set output "'//gridpath5//"ratio", vals, '.png"'
        write(gnuplotchannel,*) 'set title "Error Comparison of a Frame and Gradient with its Interpolation"'
        write(gnuplotchannel,*) 'set multiplot layout 3,1'
        write(gnuplotchannel,*) 'unset key'
        write(gnuplotchannel,FMT="(A,I5,A)") 'set label 1 "N = ', Nsamples, '" at screen 0.1,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,",",F7.4,A)') 'set label 2 "(var1,var2) = ',vals, '" at screen 0.2,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 3 "Threshhold = ',threshold_rmsd, ' A" at screen 0.5,0.925'
        write(gnuplotchannel,*) 'set style fill transparent solid 0.5'
        write(gnuplotchannel,*) 'set ylabel "Frequency"'
        write(gnuplotchannel,*) 'set yrange [0:]'
        write(gnuplotchannel,*) 'max(x, y) = (x > y ? x : y)'

        write(gnuplotchannel,*) 'min_x = ', min_rmsd1
        write(gnuplotchannel,*) 'max_x = ', max_rmsd1
        write(gnuplotchannel,*) 'set boxwidth'
        write(gnuplotchannel,*) 'set xtics ('//&
                                           '"1e-7" log10(0.0000001), '//&
                                           '"1e-6" log10(0.000001), '//&
                                           '"1e-5" log10(0.00001), '//&
                                           '"1e-4" log10(0.0001), '//&
                                           '"1e-3" log10(0.001), '//&
                                           '"1e-2" log10(0.01), '//&
                                           '"1e-1" log10(0.1), '//&
                                               '"1e0" log10(1.0), '//&
                                               '"1e1" log10(10.0))'
        write(gnuplotchannel,*) 'set xrange [1.1*min_x:max(1.1*max_x,0.9*max_x)]'
        write(gnuplotchannel,*) 'set xlabel "Error Relative to Accept Best"'
        write(gnuplotchannel,*) 'plot "'//gridpath5//"binned_ratio"//errorcheckfile//'" u '//&
                               '1:2 w boxes lc rgb "green"'

        write(gnuplotchannel,*) 'min_x = ', min_rmsd2
        write(gnuplotchannel,*) 'max_x = ', max_rmsd2
        write(gnuplotchannel,*) 'set boxwidth'
        write(gnuplotchannel,*) 'set xtics ('//&
                                           '"1e-7" log10(0.0000001), '//&
                                           '"1e-6" log10(0.000001), '//&
                                           '"1e-5" log10(0.00001), '//&
                                           '"1e-4" log10(0.0001), '//&
                                           '"1e-3" log10(0.001), '//&
                                           '"1e-2" log10(0.01), '//&
                                           '"1e-1" log10(0.1), '//&
                                               '"1e0" log10(1.0), '//&
                                               '"1e1" log10(10.0))'
        write(gnuplotchannel,*) 'set xrange [1.1*min_x:max(1.1*max_x,0.9*max_x)]'
        write(gnuplotchannel,*) 'set xlabel "Error Relative to Current Error"'
        write(gnuplotchannel,*) 'plot "'//gridpath5//"binned_ratio"//errorcheckfile//'" u '//&
                               '3:4 w boxes lc rgb "blue"'

        write(gnuplotchannel,*) 'min_x = ', min_rmsd3
        write(gnuplotchannel,*) 'max_x = ', max_rmsd3
        write(gnuplotchannel,*) 'set boxwidth'
!       write(gnuplotchannel,*) 'set boxwidth max_x-min_x /', Nbins
        write(gnuplotchannel,*) 'set xtics ('//&
                                           '"1e-7" log10(0.0000001), '//&
                                           '"1e-6" log10(0.000001), '//&
                                           '"1e-5" log10(0.00001), '//&
                                           '"1e-4" log10(0.0001), '//&
                                           '"1e-3" log10(0.001), '//&
                                           '"1e-2" log10(0.01), '//&
                                           '"1e-1" log10(0.1), '//&
                                               '"1e0" log10(1.0), '//&
                                               '"1e1" log10(10.0))'
        write(gnuplotchannel,*) 'set xrange [1.1*min_x:max(1.1*max_x,0.9*max_x)]'
        write(gnuplotchannel,*) 'set xlabel "Error Relative to Maximum Error"'
        write(gnuplotchannel,*) 'plot "'//gridpath5//"binned_ratio"//errorcheckfile//'" u '//&
                               '5:6 w boxes lc rgb "red"'
        close(gnuplotchannel)
        
        call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

!       call sleep(5)

end subroutine plotErrorCheck1

subroutine plotFinalErrorCheck1(Nsamples)
        use FUNCTIONS
        use PARAMETERS
        use PHYSICS
        use VARIABLES
        use ANALYSIS
        implicit none

        integer, intent(in) :: Nsamples
        real(dp),dimension(Nvar) :: vals
        real(dp) :: dropoff_mean,dropoff_SD
        real(dp) :: convergence_mean,convergence_SD
        real(dp),dimension(Nvar) :: min_var
        real(dp) :: min_convergence,max_convergence
        integer :: n, m, iostate

        min_var = var_maxvar
        min_convergence = 1.0d9
        max_convergence = -1.0d9
        open(filechannel2,file=gridpath5//"convergence"//errorcheckfile)
        do
                read(filechannel2,iostat=iostate,FMT=*) vals,&
                        dropoff_mean,dropoff_SD,&
                        convergence_mean,convergence_SD
                if (iostate /= 0) exit

                do n = 1, Nvar
                        min_var(n) = min(min_var(n),vals(n))
                end do
                min_convergence = min(min_convergence,convergence_mean)
                max_convergence = max(max_convergence,convergence_mean)
        end do
        close(filechannel2)

        open(gnuplotchannel,file=gridpath5//gnuplotfile)
        write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
        write(gnuplotchannel,FMT="(A)") &
                'set output "'//gridpath5//'LocalErrorConvergence.png"'
        write(gnuplotchannel,*) 'set title "Local Error Convergence over a Trajectory"'
        write(gnuplotchannel,*) 'unset key'
        write(gnuplotchannel,FMT="(A,I5,A)") 'set label 1 "N = ', Nsamples, '" at screen 0.05,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 2 "Threshhold = ',threshold_rmsd, &
                ' A" at screen 0.15,0.925'
        write(gnuplotchannel,FMT='(A,F9.4,A)') 'set label 3 "AlphaRatio = ',alpha_ratio, &
                '" at screen 0.35,0.925'
        write(gnuplotchannel,*) 'xmin = ', floor(min_var(1))
        write(gnuplotchannel,*) 'xmax = ', var_maxvar(1)
        write(gnuplotchannel,*) 'ymin = ', floor(min_var(2))
        write(gnuplotchannel,*) 'ymax = ', var_maxvar(2)
        write(gnuplotchannel,*) 'set xrange [xmin:xmax]'
        write(gnuplotchannel,*) 'set yrange [ymin:ymax]'
        write(gnuplotchannel,*) 'min_cx = ', min_convergence
        write(gnuplotchannel,*) 'max_cx = ', max_convergence
        write(gnuplotchannel,*) 'set cbrange [min_cx:max_cx]'
        write(gnuplotchannel,*) 'set palette defined ('//&
                                'log10(.0000001) "white", '//&
                                'log10(.0000001) "yellow", '//&
                                'log10(.0000005) "yellow", '//&
                                'log10(.000001) "green", '//&
                                'log10(.000005) "cyan", '//&
                                'log10(.00001) "blue", '//&
                                'log10(.00005) "magenta", '//&
                                'log10(.0001) "red"'//&
                                ')'
        write(gnuplotchannel,*) 'set xlabel "Variable 1 (A)"'
        write(gnuplotchannel,*) 'set ylabel "Variable 2 (A)"'
        write(gnuplotchannel,*) 'set cblabel "Local Error Convergence"'
        write(gnuplotchannel,*) 'set cbtics ('//&
                                             '"1e-11" log10(.00000000001), '//&
!                                            '"5e-11" log10(.00000000005), '//&
                                              '"1e-10" log10(.0000000001), '//&
!                                             '"5e-10" log10(.0000000005), '//&
                                                '"1e-9" log10(.000000001), '//&
!                                               '"5e-9" log10(.000000005), '//&
                                                 '"1e-8" log10(.00000001), '//&
!                                                '"5e-8" log10(.00000005), '//&
                                                  '"1e-7" log10(.0000001), '//&
!                                                 '"5e-7" log10(.0000005), '//&
                                                   '"1e-6" log10(.000001), '//&
!                                                  '"5e-6" log10(.000005), '//&
                                                   '"1e-5"  log10(.00001), '//&
!                                                  '"5e-5"  log10(.00005), '//&
                                                   '"1e-4"   log10(.0001), '//&
!                                                  '"5e-4"   log10(.0005), '//&
                                                   '"1e-3"    log10(.001), '//&
!                                                  '"5e-3"    log10(.005), '//&
                                                   '"1e-2"     log10(.01), '//&
!                                                  '"5e-2"     log10(.05), '//&
                                                   '"1e-1"      log10(.1), '//&
!                                                  '"5e-1"      log10(.5), '//&
                                                   ' "1.0"       log10(1), '//&
                                           ')'
        write(gnuplotchannel,*) 'plot "'//gridpath5//'convergence'//errorcheckfile//&
                '" u 1:2:5 linecolor palette pt 7 ps 1'
        close(gnuplotchannel)

        call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

end subroutine plotFinalErrorCheck1



end module runTrajectory
