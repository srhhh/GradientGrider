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


subroutine runTrajectory1(filechannels,&
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

    !Various other variables
    real(dp) :: error1,error2
    real(dp) :: min_rmsd
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
    if (gather_interpolation_flag) open(filechannel3,file=&
            gridpath5//interpolationfile,position="append")

    !Allocate all buffers; initialize the buffer size to be
    !the maximum number of frames expected to be seen
    buffer1_size = 1 + var_overcrowd(1)
    allocate(valsbuffer1(Nvar,buffer1_size),&
             coordsbuffer1(3,Natoms,buffer1_size),&
             gradientbuffer1(3,Natoms,buffer1_size),&
             Ubuffer1(3,3,buffer1_size),&
             RMSDbuffer1(buffer1_size),&
             approximation_index(buffer1_size))

    allocate(Ntrajbuffer1(buffer1_size))

    !These buffers need not be allocated if interpolation
    !is not occuring
    allocate(acceptable_frame_mask(buffer1_size),&
             inputCLS(Ncoords+buffer1_size,buffer1_size))

    do steps = 1, Nsteps

        !Check every 500 steps to see if we are out-of-bounds
        if (modulo(steps,500) == 1) then
            if (any(vals > var_maxvar)) exit
        endif

        !Upate the coordinates with the velocities
        coords = coords + dt * velocities

        !Always calculate the variables before checking a frame or accelerating
        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

        !We need to input the frame using the labeling scheme
        do n = 1, Natoms
            coords_labelled(:,n) = coords(:,BOND_LABELLING_DATA(n))
        end do

        !Check for a frame in the grid
        !Set the default value beforehand
        if (accept_worst) then
            min_rmsd = 0.0d0
        else
            min_rmsd = default_rmsd
        end if
 
        !Check the library to get an
        !approximated gradient
        call checkState_new(vals,coords_labelled,approx_gradient,min_rmsd,&
                 filechannels,number_of_frames,order,neighbor_check)

        !Accept worst has not found an acceptable frame
        !if min_rmsd is zero
        !(deprecated)
        if ((accept_worst) .and. (min_rmsd == 0.0d0)) then
            call Acceleration(vals,coords,gradient)
        
        !If we are gather interpolation data...
        else if (gather_interpolation_flag) then
 
            !Save the approx_gradient with the
            !correct labelling
            do n = 1, Natoms
                 gradient(:,BOND_LABELLING_DATA(n)) = &
                     approx_gradient(:,n)
            end do
            approx_gradient = gradient
 
            !Save the candidate_gradient with the
            !correct labelling
            do n = 1, Natoms
                 gradient(:,BOND_LABELLING_DATA(n)) = &
                     candidate_gradient(:,n)
            end do
            candidate_gradient = gradient
 
            !Calculate the true gradient
            call Acceleration(vals,coords,gradient)
 
            !If a frame was found, record data on it
            if (Ninterpolation > 0) then
                error1 = sqrt(sum((gradient - &
                        candidate_gradient)**2)/Natoms)
                error2 = sqrt(sum((gradient - &
                        approx_gradient)**2)/Natoms)
 
                write(filechannel3,FMT=*) vals(1),vals(2),&
                        Ninterpolation,largest_weighted_rmsd2,&
                        largest_weighted_rmsd,candidate_rmsd,&
                        min_rmsd,error1,error2
            end if
 
            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
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
 
            !If the approximated frame is good enough
            !and we are not rejecting it, then use it
            if ((min_rmsd .ge. threshold_RMSD) &
                    .or. (reject_flag)) then
            else
                    gradient = approx_gradient
            end if
 
        !If no interpolation data is being gathered, then the
        !true gradient is not needed
 
        !If the frame is not good enough or we are rejecting
        !frames, calculate the true gradient
        else if ((min_rmsd .ge. threshold_RMSD) &
                .or. (reject_flag)) then
            call Acceleration(vals,coords,gradient)
 
            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
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
 
        !If the approximated frame is good enough then
        !we can use it (after relabelling)!
        else
            do n = 1, Natoms
                 gradient(:,BOND_LABELLING_DATA(n)) = &
                     approx_gradient(:,n)
            end do
        end if

        !Update the velocities
        velocities = velocities + gradient
    end do

    if (gather_interpolation_flag) close(filechannel3)

    !Deallocate the buffers
    deallocate(valsbuffer1,&
            coordsbuffer1,gradientbuffer1,&
            Ubuffer1,RMSDbuffer1,&
            approximation_index)

    deallocate(Ntrajbuffer1)

    deallocate(acceptable_frame_mask,inputCLS)

    !Output the final coordinates and velocities
    coords_final = coords
    velocities_final = velocities

end subroutine runTrajectory1

subroutine runTrajectory2(filechannels,&
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

    !Various other variables
    real(dp) :: error1,error2
    real(dp) :: min_rmsd
    integer :: number_of_frames,order,neighbor_check

    !Incremental Integer
    integer :: i,n

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
    if (gather_interpolation_flag) open(filechannel3,file=&
            gridpath5//interpolationfile,position="append")

    !Allocate all buffers; initialize the buffer size to be
    !the maximum number of frames expected to be seen
    buffer1_size = 2 + Ninterpolation_max
    allocate(valsbuffer1(Nvar,buffer1_size),&
             coordsbuffer1(3,Natoms,buffer1_size),&
             gradientbuffer1(3,Natoms,buffer1_size),&
             Ubuffer1(3,3,buffer1_size),&
             RMSDbuffer1(buffer1_size),&
             approximation_index(buffer1_size))

    allocate(Ntrajbuffer1(buffer1_size))

    !These buffers need not be allocated if interpolation
    !is not occuring
    allocate(acceptable_frame_mask(buffer1_size),&
             inputCLS(Ncoords+buffer1_size,buffer1_size))

    do steps = 1, Nsteps

        !Check every 500 steps to see if we are out-of-bounds
        if (modulo(steps,500) == 1) then
            if (any(vals > var_maxvar)) exit
        endif

        !Upate the coordinates with the velocities
        coords = coords + dt * velocities

        !Always calculate the variables before checking a frame or accelerating
        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

        !Check for a frame in the grid
        !Set the default value beforehand
        if (accept_worst) then
            min_rmsd = 0.0d0
        else
            min_rmsd = default_rmsd
        end if

        RMSDbuffer1 = min_rmsd
 
        !Check the library to get an
        !approximated gradient
        call checkState_new_cap(vals,coords,approx_gradient,min_rmsd,&
                 filechannels,number_of_frames,order,neighbor_check)

        !Accept worst has not found an acceptable frame
        !if min_rmsd is zero
        !(deprecated)
        if ((accept_worst) .and. (min_rmsd == 0.0d0)) then
            call Acceleration(vals,coords,gradient)
        
        !If we are gather interpolation data...
        else if (gather_interpolation_flag) then
 
            !Calculate the true gradient
            call Acceleration(vals,coords,gradient)
 
            !If a frame was found, record data on it
            if (Ninterpolation > 0) then
                error1 = sqrt(sum((gradient - &
                        candidate_gradient)**2)/Natoms)
                error2 = sqrt(sum((gradient - &
                        approx_gradient)**2)/Natoms)
 
                write(filechannel3,FMT=*) vals(1),vals(2),&
                        Ninterpolation,largest_weighted_rmsd2,&
                        largest_weighted_rmsd,candidate_rmsd,&
                        min_rmsd,error1,error2
            end if
 
            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
            if (grid_addition > 0) then
                do i = 1, Nindistinguishables
                    BOND_LABELLING_DATA = &
                        INDISTINGUISHABLES(i,:)
                    do n = 1, Natoms
                         gradient_labelled(:,n) = &
                             gradient(:,BOND_LABELLING_DATA(n))
                         coords_labelled(:,n) = &
                             coords(:,BOND_LABELLING_DATA(n))
                    end do
                    call addState_new(vals,&
                            coords_labelled,&
                            gradient_labelled)
                end do
            end if
 
            !If the approximated frame is good enough
            !and we are not rejecting it, then use it
            if ((min_rmsd .ge. threshold_RMSD) &
                    .or. (reject_flag)) then
            else
                    gradient = approx_gradient
            end if
 
        !If no interpolation data is being gathered, then the
        !true gradient is not needed
 
        !If the frame is not good enough or we are rejecting
        !frames, calculate the true gradient
        else if ((min_rmsd .ge. threshold_RMSD) &
                .or. (reject_flag)) then
            call Acceleration(vals,coords,gradient)
 
            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
            if (grid_addition > 0) then
                do i = 1, Nindistinguishables
                    BOND_LABELLING_DATA = &
                        INDISTINGUISHABLES(i,:)
                    do n = 1, Natoms
                         gradient_labelled(:,n) = &
                             gradient(:,BOND_LABELLING_DATA(n))
                         coords_labelled(:,n) = &
                             coords(:,BOND_LABELLING_DATA(n))
                    end do
                    call addState_new(vals,&
                            coords_labelled,&
                            gradient_labelled)
                end do
            end if
 
        !If the approximated frame is good enough then
        !we can use it (after relabelling)!
        else
            gradient = approx_gradient
        end if

        !Update the velocities
        velocities = velocities + gradient
    end do

    if (gather_interpolation_flag) close(filechannel3)

    !Deallocate the buffers
    deallocate(valsbuffer1,&
            coordsbuffer1,gradientbuffer1,&
            Ubuffer1,RMSDbuffer1,&
            approximation_index)

    deallocate(Ntrajbuffer1)

    deallocate(acceptable_frame_mask,inputCLS)

    !Output the final coordinates and velocities
    coords_final = coords
    velocities_final = velocities

end subroutine runTrajectory2

subroutine runTrajectory3(filechannels,&
                          coords_initial,velocities_initial,&
                          coords_final,velocities_final)
    use PARAMETERS
    use ls_rmsd_original
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

    !Various other variables
    real(dp) :: error1,error2
    real(dp),dimension(3) :: x_center, y_center
    real(dp) :: min_rmsd
    integer :: number_of_frames,order,neighbor_check

    real(dp) :: RMSD_prior
    real(dp) :: max_RMSDdelta, max_deltaRMSD
    real(dp),dimension(3,Natoms) :: coords_prior

    real(dp),allocatable :: libcoords(:,:,:), libgradients(:,:,:)
    integer,allocatable :: libNtraj(:)

    !Incremental Integer
    integer :: i,n,m

    max_RMSDdelta = 0.0d0
    max_deltaRMSD = 0.0d0

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
    if (gather_interpolation_flag) open(filechannel3,file=&
            gridpath5//interpolationfile,position="append")

    !Allocate all buffers; initialize the buffer size to be
    !the maximum number of frames expected to be seen
    buffer1_size = 2 + Ninterpolation_max
    allocate(valsbuffer1(Nvar,buffer1_size),&
             coordsbuffer1(3,Natoms,buffer1_size),&
             gradientbuffer1(3,Natoms,buffer1_size),&
             Ubuffer1(3,3,buffer1_size),&
             RMSDbuffer1(buffer1_size),&
             inputCLS(Ncoords+buffer1_size,buffer1_size))

    RMSDbuffer1 = default_rmsd
     
!   allocate(acceptable_frame_mask(buffer1_size),&
!            approximation_index(buffer1_size))

!   allocate(Ntrajbuffer1(buffer1_size))

    do steps = 1, Nsteps

        !Check every 500 steps to see if we are out-of-bounds
        if (modulo(steps,500) == 1) then
            if (any(vals > var_maxvar)) exit
        endif

        !Upate the coordinates with the velocities
        coords = coords + dt * velocities

        !Always calculate the variables before checking a frame or accelerating
        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

        !Check for a frame in the grid
        !Set the default value beforehand
        if (accept_worst) then
            min_rmsd = 0.0d0
        else
            min_rmsd = default_rmsd
        end if

        RMSDbuffer1 = min_rmsd
 
        !Check the library to get an
        !approximated gradient
        call checkState_new_permute_cap(vals,coords,approx_gradient,min_rmsd,&
                 filechannels,number_of_frames,order,neighbor_check)

!       error2 = 1.0d9
!       do i = 1, Nindistinguishables

!           BOND_LABELLING_DATA = INDISTINGUISHABLES(i,:)
    
!           !We need to input the frame using the labeling scheme
!           do n = 1, Natoms
!               coords_labelled(:,n) = coords(:,BOND_LABELLING_DATA(n))
!           end do
    
!           call rmsd_dp(Natoms,coords_labelled,&
!                   coords_initial,1,Ubuffer1(:,:,1),&
!                   x_center,y_center,error1)

!           if (error1 < error2) error2 = error1
!       end do

!       print *, ""
!       print *, steps
!       print *, vals
!       print *, error2

!       if (steps > 1) then
!           call rmsd_dp(Natoms,coords,&
!                   coords_initial,1,Ubuffer1(:,:,1),&
!                   x_center,y_center,error1)

!           if (abs(error1 - RMSD_prior) > max_deltaRMSD) &
!                   max_deltaRMSD = abs(error1 - RMSD_prior)
!           RMSD_prior = error1

!           call rmsd_dp(Natoms,coords,&
!                   coords_prior,1,Ubuffer1(:,:,1),&
!                   x_center,y_center,error1)

!           if (error1 > max_RMSDdelta) &
!                   max_RMSDdelta = error1
!           coords_prior = coords
!       else
!           call rmsd_dp(Natoms,coords,&
!                   coords_initial,1,Ubuffer1(:,:,1),&
!                   x_center,y_center,error1)

!           RMSD_prior = error1
!           coords_prior = coords
!       end if

        !Accept worst has not found an acceptable frame
        !if min_rmsd is zero
        !(deprecated)
        if ((accept_worst) .and. (min_rmsd == 0.0d0)) then
            call Acceleration(vals,coords,gradient)
        
        !If we are gather interpolation data...
        else if (gather_interpolation_flag) then
 
            !Calculate the true gradient
            call Acceleration(vals,coords,gradient)

            !If a frame was found, record data on it
            if (Ninterpolation > 0) then
                error1 = sqrt(sum((gradient - &
                        candidate_gradient)**2)/Natoms)
                error2 = sqrt(sum((gradient - &
                        approx_gradient)**2)/Natoms)
 
                write(filechannel3,FMT=*) vals(1),vals(2),&
                        Ninterpolation,largest_weighted_rmsd2,&
                        largest_weighted_rmsd,candidate_rmsd,&
                        min_rmsd,error1,error2

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!              if ((Ninterpolation > 10).and.(error1 < error2)) then
!              if ((Ninterpolation > 10)) then
!                  allocate(libcoords(Ninterpolation,3,Natoms),&
!                           libgradients(Ninterpolation,3,Natoms),&
!                           libNtraj(Ninterpolation))
!                  m = 0
!                  do n = 1, Ninterpolation
!                      do
!                      m = m + 1
!                      if (acceptable_frame_mask(m)) then
!                          libcoords(n,:,:) =&
!                                  coordsbuffer1(:,:,m)
!                          libgradients(n,:,:) =&
!                                  gradientbuffer1(:,:,m)
!                         !libNtraj(n) =&
!                         !        Ntrajbuffer1(m)
!                          exit
!                      else
!                      end if
!                      end do
!                  end do

!                  libNtraj = 0

!                  call errorCheck5(filechannels,&
!                          coords,gradient,&
!                          Ninterpolation,&
!                          libcoords,libgradients,&
!                          libNtraj)

!                  deallocate(libcoords,libgradients,libNtraj)
!              end if

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            end if

            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
            if (grid_addition > 0) then
                call addState_new(vals,&
                        coords,&
                        gradient)
            end if
 
            !If the approximated frame is good enough
            !and we are not rejecting it, then use it
            if ((min_rmsd .ge. threshold_RMSD) &
                    .or. (reject_flag)) then
            else
                    gradient = approx_gradient
            end if
 
        !If no interpolation data is being gathered, then the
        !true gradient is not needed
 
        !If the frame is not good enough or we are rejecting
        !frames, calculate the true gradient
        else if ((min_rmsd .ge. threshold_RMSD) &
                .or. (reject_flag)) then
            call Acceleration(vals,coords,gradient)
 
            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
            if (grid_addition > 0) then
                call addState_new(vals,&
                        coords,&
                        gradient)
            end if
 
        !If the approximated frame is good enough then
        !we can use it (after relabelling)!
        else
            gradient = approx_gradient
        end if

        !Update the velocities
        velocities = velocities + gradient
    end do

    if (gather_interpolation_flag) close(filechannel3)

    !Deallocate the buffers
    deallocate(valsbuffer1,&
            coordsbuffer1,gradientbuffer1,&
            Ubuffer1,RMSDbuffer1,&
            inputCLS)

!   deallocate(Ntrajbuffer1)

!   deallocate(acceptable_frame_mask,approximation_index)

    !Output the final coordinates and velocities
    coords_final = coords
    velocities_final = velocities

!   print *, "max_RMSDdelta"
!   print *, max_RMSDdelta
!   print *, "max_deltaRMSD"
!   print *, max_deltaRMSD

end subroutine runTrajectory3

subroutine runTestTrajectory(filechannels,&
                             coords_initial,velocities_initial,&
                             coords_final,velocities_final)
    use PARAMETERS
    use ls_rmsd_original
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

    !Various other variables
    real(dp) :: error1,error2
    real(dp),dimension(3) :: x_center, y_center
    real(dp),dimension(3,3) ::  U_prior
    real(dp) :: min_rmsd,min_rmsd_prime
    integer :: number_of_frames,order,neighbor_check
    real(dp) :: U, KE

    real(dp) :: RMSD_prior, error_prior
    real(dp) :: max_RMSDdelta, max_deltaRMSD
    real(dp),dimension(3,Natoms) :: coords_prior, gradient_prior
    real(dp) :: E, E_prior

    real(dp),allocatable :: libcoords(:,:,:), libgradients(:,:,:)
    integer,allocatable :: libNtraj(:)

    integer,dimension(Nalpha) :: alpha_flagging

    !Incremental Integer
    integer :: i,j,n,m

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
    if (gather_interpolation_flag) then
        open(filechannel2,file=&
            gridpath5//checkstatefile,position="append")
        open(filechannel3,file=&
            gridpath5//interpolationfile,position="append")
!       open(filechannel2,file=&
!           gridpath5//checkstatefile)
!       open(filechannel3,file=&
!           gridpath5//interpolationfile)
        open(6666,file=&
            gridpath5//alphafile,position="append")
    else
        open(filechannel2,file=gridpath5//checkstatefile)
    end if

    !Allocate all buffers; initialize the buffer size to be
    !the maximum number of frames expected to be seen
    buffer1_size = 2 + Ninterpolation_max
    allocate(valsbuffer1(Nvar,buffer1_size),&
             coordsbuffer1(3,Natoms,buffer1_size),&
             gradientbuffer1(3,Natoms,buffer1_size),&
             Ubuffer1(3,3,buffer1_size),&
             RMSDbuffer1(buffer1_size),&
             inputCLS(Ncoords+buffer1_size,buffer1_size))

    RMSDbuffer1 = default_rmsd
     
    do steps = 1, Nsteps

        !Check every 500 steps to see if we are out-of-bounds
        if (modulo(steps,500) == 1) then
            if (any(vals > var_maxvar)) exit
        endif

        !Upate the coordinates with the velocities
        coords = coords + dt * velocities

        !Always calculate the variables before checking a frame or accelerating
        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

        !Check for a frame in the grid
        !Set the default value beforehand
        if (accept_worst) then
            min_rmsd = 0.0d0
        else
            min_rmsd = default_rmsd
        end if
        
        order = 0

        RMSDbuffer1 = min_rmsd

        subcellsearch_max = subcellsearch_max1
 
        !Check the library to get an
        !approximated gradient
        call checkState_new_cap(vals,coords,approx_gradient,min_rmsd,&
                 filechannels,number_of_frames,m,neighbor_check)
        order = order + m

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      if ((Ninterpolation > 10)) then
!          allocate(libcoords(Ninterpolation,3,Natoms),&
!                   libgradients(Ninterpolation,3,Natoms))
!          do n = 1, Ninterpolation
!              libcoords(n,:,:) =&
!                      coordsbuffer1(:,:,n)
!              libgradients(n,:,:) =&
!                      gradientbuffer1(:,:,n)
!          end do

!          call errorCheck8(filechannels,&
!                  coords,gradient,&
!                  Ninterpolation,&
!                  libcoords,libgradients)

!          deallocate(libcoords,libgradients)
!      end if

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        min_rmsd = candidate_rmsd

        !Set a separate default value
        if (accept_worst) then
            min_rmsd_prime = 0.0d0
        else
            min_rmsd_prime = default_rmsd
        end if

        RMSDbuffer1 = min_rmsd_prime
    
        !The second search uses interpolation as it
        !usually is a wider search
        subcellsearch_max = subcellsearch_max2

        call checkState_new_permute_cap(vals,coords,&
                approx_gradient_prime,min_rmsd_prime,&
                filechannels,number_of_frames,m,neighbor_check)
        order = order + m

!       !For the B -> B' test
!       if (min_rmsd > threshold_rmsd) Ninterpolation = 0

        call getEnergies(Natoms,coords,velocities,U,KE)
        E = U + KE
    
        !Finally write to the data file all the important data values
        write(filechannel2,FMT=*) number_of_frames,order,neighbor_check,steps,&
                                  min_rmsd,min_rmsd_prime,vals(1),vals(2),U,KE

        !The better frame is usually found with second search
        !so this will be used for further data handling
        min_rmsd = min_rmsd_prime
        approx_gradient = approx_gradient_prime




        !Accept worst has not found an acceptable frame
        !if min_rmsd is zero
        !(deprecated)
        if ((accept_worst) .and. (min_rmsd == 0.0d0)) then
            call Acceleration(vals,coords,gradient)
        
        !If we are gather interpolation data...
        else if (gather_interpolation_flag) then
 
            !Calculate the true gradient
            call Acceleration(vals,coords,gradient)

!           !!!!!!!!!!!!!
!           ! TEST START
!           !!!!!!!!!!!!!

!           if (steps > 1) then
!               RMSD_prior = sqrt(sum((coords-coords_prior)**2)/Natoms)
!               error_prior = sqrt(sum((gradient-gradient_prior)**2)/Natoms)

!               print *, ""
!               print *, steps
!               print *, RMSD_prior
!               print *, error_prior
!           end if
!           coords_prior = coords
!           gradient_prior = gradient

!           !!!!!!!!!!!!!
!           ! TEST END
!           !!!!!!!!!!!!!

            !If a frame was found, record data on it
            if (Ninterpolation > 0) then
                error1 = sqrt(sum((gradient - &
                        candidate_gradient)**2)/Natoms)
                error2 = sqrt(sum((gradient - &
                        approx_gradient)**2)/Natoms)
 
                if ((steps > 1).and.(abs(E-E_prior) > E_threshold)) then
                    allocate(libcoords(Ninterpolation,3,Natoms),&
                             libgradients(Ninterpolation,3,Natoms))

                    do n = 1, Ninterpolation
                        libcoords(n,:,:) =&
                                coordsbuffer1(:,:,n)
                        libgradients(n,:,:) =&
                                gradientbuffer1(:,:,n)
                    end do

                    call errorCheck8(filechannels,&
                            coords,gradient,&
                            Ninterpolation,&
                            libcoords,libgradients,&
                            alpha_flagging)

                    deallocate(libcoords,libgradients)

!                   write(6666,FMT=*) alpha_flagging
                end if
 
                write(filechannel3,FMT=*) vals(1),vals(2),&
                        Ninterpolation,largest_weighted_rmsd2,&
                        largest_weighted_rmsd,candidate_rmsd,&
                        min_rmsd,error1,error2
            end if

            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
            if (grid_addition > 0) then
                call addState_new(vals,&
                        coords,&
                        gradient)
            end if
 
            !If the approximated frame is good enough
            !and we are not rejecting it, then use it
            if ((min_rmsd .ge. threshold_RMSD) &
                    .or. (reject_flag)) then
            else
                    gradient = approx_gradient
            end if
 
        !If no interpolation data is being gathered, then the
        !true gradient is not needed
 
        !If the frame is not good enough or we are rejecting
        !frames, calculate the true gradient
        else if ((min_rmsd .ge. threshold_RMSD) &
                .or. (reject_flag)) then
            call Acceleration(vals,coords,gradient)
 
            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
            if (grid_addition > 0) then
                call addState_new(vals,&
                        coords,&
                        gradient)
            end if
 
        !If the approximated frame is good enough then
        !we can use it (after relabelling)!
        else
            gradient = approx_gradient
        end if

!       !!!!!!!!!!!!!
!       ! TEST START
!       !!!!!!!!!!!!!

!       if ((reject_flag).and.(modulo(steps,2)==0)) then
!           call rmsd_dp(Natoms,coords,coords_prior,&
!                   1,U_prior,x_center,y_center,RMSD_prior)
!           gradient = matmul(U_prior,gradient_prior)
!!          gradient = gradient_prior
!       else
!           coords_prior = coords
!           gradient_prior = gradient
!       end if

!       !!!!!!!!!!!!!
!       ! TEST END
!       !!!!!!!!!!!!!

        !Update the velocities
        velocities = velocities + gradient
    end do

    close(filechannel2)

    if (gather_interpolation_flag) then
        close(filechannel3)
        close(6666)
    end if

    !Deallocate the buffers
    deallocate(valsbuffer1,&
            coordsbuffer1,gradientbuffer1,&
            Ubuffer1,RMSDbuffer1,&
            inputCLS)

    !Output the final coordinates and velocities
    coords_final = coords
    velocities_final = velocities

end subroutine runTestTrajectory

subroutine runTestTrajectory2(filechannels,&
                              coords_initial,velocities_initial,&
                              coords_final,velocities_final)
    use PARAMETERS
    use ls_rmsd_original
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

    !Various other variables
    real(dp) :: error1,error2
    real(dp),dimension(3) :: x_center, y_center
    real(dp),dimension(3,3) ::  U_prior
    real(dp) :: min_rmsd,min_rmsd_prime
    integer :: number_of_frames,order,neighbor_check
    real(dp) :: U, KE

    real(dp) :: RMSD_prior, error_prior
    real(dp) :: max_RMSDdelta, max_deltaRMSD
    real(dp),dimension(3,Natoms) :: coords_prior, gradient_prior
    real(dp),dimension(3,Natoms) :: velocities_prior
    real(dp) :: DE, E, E_prior, E_baseline
    integer :: Nalpha_tries = 1
    logical :: temp_reject_flag = .false.

    real(dp),allocatable :: libcoords(:,:,:), libgradients(:,:,:)
    integer,allocatable :: libNtraj(:)

    integer,dimension(Nalpha) :: alpha_flagging

    !Incremental Integer
    integer :: i,j,n,m

    !Initialize the scene
    call InitialSetup3(coords,velocities)
    Norder1 = 0
    Norder_total = 0
    allocate(trajRMSDbuffer(Ngrid_max,Naccept_max+1))

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
    if (gather_interpolation_flag) then
        open(filechannel2,file=&
            gridpath5//checkstatefile,position="append")
        open(filechannel3,file=&
            gridpath5//interpolationfile,position="append")
        open(6666,file=&
            gridpath5//alphafile,position="append")
    else
        open(filechannel2,file=gridpath5//checkstatefile)
    end if

    !Allocate all buffers; initialize the buffer size to be
    !the maximum number of frames expected to be seen
    buffer1_size = 2 + Ninterpolation_max
    allocate(valsbuffer1(Nvar,buffer1_size),&
             coordsbuffer1(3,Natoms,buffer1_size),&
             gradientbuffer1(3,Natoms,buffer1_size),&
             Ubuffer1(3,3,buffer1_size),&
             RMSDbuffer1(buffer1_size),&
             inputCLS(Ncoords+buffer1_size,buffer1_size))

    allocate(CMdiffbuffer1(buffer1_size))
    allocate(Ntrajbuffer1(buffer1_size))

    RMSDbuffer1 = default_rmsd

    steps = 1
    E_baseline = 0.0d0
    do
        coords = coords + dt * velocities
        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
        call Acceleration(vals,coords,gradient)

        call getEnergies(Natoms,coords,velocities,U,KE)
        E_baseline = E_baseline + U + KE

        !Finally write to the data file all the important data values
        write(filechannel2,FMT=*) 0,0,0,steps,&
                                  default_rmsd,default_rmsd,&
                                  vals(1),vals(2),U,KE

        do i = 1, Ngrid_total
            write(filechannels(1+i),FMT=FMT6)&
                default_rmsd
        end do

        velocities = velocities + gradient

        steps = steps + 1
        if (steps == Nsteps_baseline) then
            coords = coords + dt * velocities
            velocities = velocities + gradient

            call getEnergies(Natoms,coords,velocities,U,KE)
            E = U + KE
            E_baseline = E_baseline + E

            coords_prior = coords
            velocities_prior = velocities
            E_prior = E

            exit
        end if
    end do
    E_baseline = E_baseline / Nsteps_baseline
    Naccept = 0
    Nalpha_tries = 1
    alpha_ratio = alpha_ratio_list(1)
     
    do

        !Always calculate the variables before checking a frame or accelerating
        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

!       if (Nalpha_tries > 1) then
!               print *, "step", steps
!               print *, "Naccept", Naccept
!               print *, vals
!       end if

        !Check for a frame in the grid
        !Set the default value beforehand
        if (accept_worst) then
            min_rmsd = 0.0d0
        else
            min_rmsd = default_rmsd
        end if
        
        order = 0

        RMSDbuffer1 = min_rmsd

        subcellsearch_max = subcellsearch_max1
 
        !Check the library to get an
        !approximated gradient
        call checkState_new_cap(vals,coords,approx_gradient,min_rmsd,&
                 filechannels,number_of_frames,m,neighbor_check)
        order = order + m

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      if ((Ninterpolation > 10)) then
!          allocate(libcoords(Ninterpolation,3,Natoms),&
!                   libgradients(Ninterpolation,3,Natoms))
!          do n = 1, Ninterpolation
!              libcoords(n,:,:) =&
!                      coordsbuffer1(:,:,n)
!              libgradients(n,:,:) =&
!                      gradientbuffer1(:,:,n)
!          end do

!          call errorCheck8(filechannels,&
!                  coords,gradient,&
!                  Ninterpolation,&
!                  libcoords,libgradients)

!          deallocate(libcoords,libgradients)
!      end if

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        min_rmsd = candidate_rmsd

        !Set a separate default value
        if (accept_worst) then
            min_rmsd_prime = 0.0d0
        else
            min_rmsd_prime = default_rmsd
        end if

        RMSDbuffer1 = min_rmsd_prime
    
        !The second search uses interpolation as it
        !usually is a wider search
        subcellsearch_max = subcellsearch_max2

        call checkState_new_permute_cap(vals,coords,&
                approx_gradient_prime,min_rmsd_prime,&
                filechannels,number_of_frames,m,neighbor_check)
        order = order + m

!       !For the B -> B' test
!       if (min_rmsd > threshold_rmsd) Ninterpolation = 0
    
            !If the approximated frame is good enough
            !and we are not rejecting it, then use it
            if (((min_rmsd_prime .ge. threshold_RMSD) &
                    .or. (reject_flag)).and.&
                 (Naccept == 0)) then

                call Acceleration(vals,coords,gradient)

                do i = 1, Ngrid_total
                    write(filechannels(1+i),FMT=FMT6)&
                        trajRMSDbuffer(i,1)
                end do

!               coords_prior = coords
!               velocities_prior = velocities
!               E_prior = E
!   
!               !If we are adding to a grid, then relabel
!               !the coordinates and gradients, and do so
!               if (grid_addition > 0) then
!                   call addState_new(vals,&
!                           coords,&
!                           gradient)
!               end if

!               if ((any(vals > var_maxvar)).or.&
!                   (steps > Nsteps)) exit
 
            !If the approximated frame is good enough
            !and we are not rejecting it, then use it
            else if ((min_rmsd_prime .ge. threshold_RMSD) &
                    .or. (reject_flag).or.&
                     (Naccept == Naccept_max)) then

                call getEnergies(Natoms,coords,velocities,U,KE)
                E = U + KE

                    DE = abs(E-E_prior)/Naccept
                    if ((abs(E-E_baseline) > E_threshold)&
                        .or.(DE > DE_threshold)) then
        
                        coords = coords_prior
                        velocities = velocities_prior
                        steps = steps - Naccept
    
                        Nalpha_tries = Nalpha_tries + 1
                        if (Nalpha_tries > Nalpha_tries_max) then
                            alpha_ratio = alpha_ratio_list(1)
                        else
                            alpha_ratio = alpha_ratio_list(Nalpha_tries)
                        end if

                    else
                        do i = 1, Ngrid_total
                        do j = 1, Naccept
                            write(filechannels(1+i),FMT=FMT6)&
                                trajRMSDbuffer(i,j)
                        end do
                        end do
 
                        Nalpha_tries = 1
                        alpha_ratio = alpha_ratio_list(Nalpha_tries)

                        coords_prior = coords
                        velocities_prior = velocities
                        E_prior = E
            
                        !If we are adding to a grid, then relabel
                        !the coordinates and gradients, and do so
                        if (grid_addition > 0) then
                            call addState_new(vals,&
                                    coords,&
                                    gradient)
                        end if
        
                        if ((any(vals > var_maxvar)).or.&
                            (steps > Nsteps)) exit
 
                    end if

                Naccept = 0
                cycle

            else if (Nalpha_tries > Nalpha_tries_max) then

                call Acceleration(vals,coords,gradient)

                do i = 1, Ngrid_total
                    write(filechannels(1+i),FMT=FMT6)&
                        default_rmsd
                end do

                Nalpha_tries = 1

!               coords_prior = coords
!               velocities_prior = velocities
!               E_prior = E
!   
!               !If we are adding to a grid, then relabel
!               !the coordinates and gradients, and do so
!               if (grid_addition > 0) then
!                   call addState_new(vals,&
!                           coords,&
!                           gradient)
!               end if

!               if ((any(vals > var_maxvar)).or.&
!                   (steps > Nsteps)) exit
 
            else
                gradient = approx_gradient
                Naccept = Naccept + 1
            end if
 
!       !!!!!!!!!!!!!
!       ! TEST START
!       !!!!!!!!!!!!!

!       if ((reject_flag).and.(modulo(steps,2)==0)) then
!           call rmsd_dp(Natoms,coords,coords_prior,&
!                   1,U_prior,x_center,y_center,RMSD_prior)
!           gradient = matmul(U_prior,gradient_prior)
!!          gradient = gradient_prior
!       else
!           coords_prior = coords
!           gradient_prior = gradient
!       end if

!       !!!!!!!!!!!!!
!       ! TEST END
!       !!!!!!!!!!!!!

        call getEnergies(Natoms,coords,velocities,U,KE)
        E = U + KE

        !Finally write to the data file all the important data values
        write(filechannel2,FMT=*) number_of_frames,order,neighbor_check,steps,&
                                  min_rmsd,min_rmsd_prime,vals(1),vals(2),U,KE

        !The better frame is usually found with second search
        !so this will be used for further data handling
        min_rmsd = min_rmsd_prime
        approx_gradient = approx_gradient_prime

        if (gather_interpolation_flag) then

            call Acceleration(vals,coords,gradient)

!           !!!!!!!!!!!!!
!           ! TEST START
!           !!!!!!!!!!!!!

!           if (steps > 1) then
!               RMSD_prior = sqrt(sum((coords-coords_prior)**2)/Natoms)
!               error_prior = sqrt(sum((gradient-gradient_prior)**2)/Natoms)

!               print *, ""
!               print *, steps
!               print *, RMSD_prior
!               print *, error_prior
!           end if
!           coords_prior = coords
!           gradient_prior = gradient

!           !!!!!!!!!!!!!
!           ! TEST END
!           !!!!!!!!!!!!!

            !If a frame was found, record data on it
            if (Ninterpolation > 0) then
                error1 = sqrt(sum((gradient - &
                        candidate_gradient)**2)/Natoms)
                error2 = sqrt(sum((gradient - &
                        approx_gradient)**2)/Natoms)
 
!               if ((steps > 1).and.(abs(E-E_prior) > E_threshold)) then
!                   allocate(libcoords(Ninterpolation,3,Natoms),&
!                            libgradients(Ninterpolation,3,Natoms))

!                   do n = 1, Ninterpolation
!                       libcoords(n,:,:) =&
!                               coordsbuffer1(:,:,n)
!                       libgradients(n,:,:) =&
!                               gradientbuffer1(:,:,n)
!                   end do

!                   call errorCheck8(filechannels,&
!                           coords,gradient,&
!                           Ninterpolation,&
!                           libcoords,libgradients,&
!                           alpha_flagging)

!                   deallocate(libcoords,libgradients)

!                   write(6666,FMT=*) alpha_flagging
!               end if
 
                write(filechannel3,FMT=*) vals(1),vals(2),&
                        Ninterpolation,largest_weighted_rmsd2,&
                        largest_weighted_rmsd,candidate_rmsd,&
                        min_rmsd,error1,error2
            end if
        end if

        !Update the velocities
        velocities = velocities + gradient

        !Upate the coordinates with the velocities
        coords = coords + dt * velocities

        !Update the steps
        steps = steps + 1

        if (Naccept == 0) then
            coords_prior = coords
            velocities_prior = velocities
            E_prior = E

            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
            if (grid_addition > 0) then
                call addState_new(vals,&
                        coords,&
                        gradient)
            end if

            if ((any(vals > var_maxvar)).or.&
                (steps > Nsteps)) exit
        end if

    end do

    close(filechannel2)

    if (gather_interpolation_flag) then
        close(filechannel3)
        close(6666)
    end if

    !Deallocate the buffers
    deallocate(valsbuffer1,&
            coordsbuffer1,gradientbuffer1,&
            Ubuffer1,RMSDbuffer1,&
            inputCLS)

    deallocate(trajRMSDbuffer)

    deallocate(CMdiffbuffer1)
    deallocate(Ntrajbuffer1)

    !Output the final coordinates and velocities
    coords_final = coords
    velocities_final = velocities

end subroutine runTestTrajectory2





subroutine runTrajectoryRewind1(filechannels,&
           coords_initial,velocities_initial,&
           coords_final,velocities_final)
    use PARAMETERS
    use ls_rmsd_original
    use PHYSICS
    use VARIABLES
    use ANALYSIS
    use interactMultipleGrids
    use analyzeRMSDThresholdwithMultipleGrids
    implicit none

    !Coordinates, Velocities, and Variables
    real(dp), dimension(3,Natoms) :: coords,coords_labelled,velocities
    real(dp), dimension(3,Natoms) :: gradient,gradient_labelled
    real(dp), dimension(3,Natoms) :: approx_gradient
    real(dp), dimension(Nvar) :: vals
    real(dp), dimension(3,Natoms), intent(out) :: coords_initial, velocities_initial 
    real(dp), dimension(3,Natoms), intent(out) :: coords_final, velocities_final 
    integer,dimension(1+Ngrid_max) :: filechannels

    !Various other variables
    real(dp) :: error1,error2
    real(dp),dimension(3) :: x_center, y_center
    real(dp),dimension(3,3) ::  U_prior
    real(dp) :: min_rmsd
    integer :: number_of_frames,order,neighbor_check
    real(dp) :: U, KE

    real(dp) :: RMSD_prior, error_prior
    real(dp) :: max_RMSDdelta, max_deltaRMSD
    real(dp),dimension(3,Natoms) :: coords_prior, gradient_prior
    real(dp),dimension(3,Natoms) :: velocities_prior
    real(dp) :: DE, E, E_prior, E_baseline
    integer :: Nalpha_tries = 1

    real(dp),allocatable :: libcoords(:,:,:), libgradients(:,:,:)
    integer,allocatable :: libNtraj(:)

    integer,dimension(Nalpha) :: alpha_flagging

    !Incremental Integer
    integer :: i,j,n,m

    open(filechannel2,file=gridpath5//checkstatefile)




    !Initialize the scene
    call InitialSetup3(coords,velocities)
    Norder1 = 0
    Norder_total = 0
    allocate(trajRMSDbuffer(Ngrid_max,Naccept_max+1))

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
    if (gather_interpolation_flag) &
        open(filechannel3,file=&
            gridpath5//interpolationfile,position="append")

    !Allocate all buffers; initialize the buffer size to be
    !the maximum number of frames expected to be seen
    buffer1_size = 2 + Ninterpolation_max
    buffer2_size = 30
!   allocate(valsbuffer1(Nvar,buffer1_size),&
!            coordsbuffer1(3,Natoms,buffer1_size),&
!            gradientbuffer1(3,Natoms,buffer1_size),&
!            Ubuffer1(3,3,buffer1_size),&
!            RMSDbuffer1(buffer1_size),&
!            inputCLS(Ncoords+buffer1_size,buffer1_size))

!   allocate(CMdiffbuffer1(buffer1_size))
!   allocate(Ntrajbuffer1(buffer1_size))

    call setAllocations()

    RMSDbuffer1 = default_rmsd
    CMdiffbuffer1 = default_CMdiff

    steps = 1
    E_baseline = 0.0d0
    do
        coords = coords + dt * velocities
        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
        call Acceleration(vals,coords,gradient)

        call getEnergies(Natoms,coords,velocities,U,KE)
        E_baseline = E_baseline + U + KE

        !Finally write to the data file all the important data values
        write(filechannel2,FMT=*) 0,0,0,steps,&
                                  default_rmsd,default_rmsd,&
                                  vals(1),vals(2),U,KE

        do i = 1, Ngrid_total
            write(filechannels(1+i),FMT=FMT6)&
                default_rmsd
        end do

        velocities = velocities + gradient

        steps = steps + 1
        if (steps == Nsteps_baseline) then
            call getVarsMaxMin(coords,Natoms,vals,&
                    Nvar,BOND_LABELLING_DATA)

            coords = coords + dt * velocities
            velocities = velocities + gradient

            call getEnergies(Natoms,coords,velocities,U,KE)
            E = U + KE
            E_baseline = E_baseline + E

            coords_prior = coords
            velocities_prior = velocities
            E_prior = E

            exit
        end if
    end do
    E_baseline = E_baseline / Nsteps_baseline
    Naccept = 0
    Nalpha_tries = 1
    alpha_ratio = alpha_ratio_list(1)

    populationbuffer2 = -1
    do i = 1, Nvar
        previous_var_index(i) = int(vals(i)*divisor(i,Norder_order(1)+1))
    end do
     
    do
        !Always calculate the variables before checking a frame or accelerating
        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

        !Check for a frame in the grid
        !Set the default value beforehand
        if (accept_worst) then
            min_rmsd = 0.0d0
        else
            min_rmsd = default_rmsd
        end if
        
        RMSDbuffer1 = min_rmsd
        CMdiffbuffer1 = default_CMdiff

        subcellsearch_max = subcellsearch_max1
 
        call checkState_PCM(vals,coords,&
                approx_gradient,min_rmsd,&
                filechannels,number_of_frames,order,neighbor_check)
!       print *, ""
!       print *, "step", steps
!       if (Ninterpolation > 0) &
!       print *, RMSDbuffer1(1:Ninterpolation)

        !If the approximated frame is good enough
        !and we are not rejecting it, then use it
        if (((min_rmsd .ge. threshold_RMSD) &
                .or. (reject_flag)).and.&
             (Naccept == 0)) then

            call Acceleration(vals,coords,gradient)

            do i = 1, Ngrid_total
                write(filechannels(1+i),FMT=FMT6)&
                    trajRMSDbuffer(i,1)
            end do
 
        !If the approximated frame is good enough
        !and we are not rejecting it, then use it
        else if ((min_rmsd .ge. threshold_RMSD) &
                .or. (reject_flag).or.&
                 (Naccept == Naccept_max)) then

            call getEnergies(Natoms,coords,velocities,U,KE)
            E = U + KE

            DE = abs(E-E_prior)/Naccept
            if ((abs(E-E_baseline) > E_threshold)&
                .or.(DE > DE_threshold)) then

                coords = coords_prior
                velocities = velocities_prior
                steps = steps - Naccept

                Nalpha_tries = Nalpha_tries + 1
                if (Nalpha_tries > Nalpha_tries_max) then
                    alpha_ratio = alpha_ratio_list(1)
                else
                    alpha_ratio = alpha_ratio_list(Nalpha_tries)
                end if

            else
                do i = 1, Ngrid_total
                do j = 1, Naccept
                    write(filechannels(1+i),FMT=FMT6)&
                        trajRMSDbuffer(i,j)
                end do
                end do
 
                Nalpha_tries = 1
                alpha_ratio = alpha_ratio_list(Nalpha_tries)

                coords_prior = coords
                velocities_prior = velocities
                E_prior = E
    
                !If we are adding to a grid, then relabel
                !the coordinates and gradients, and do so
                if (grid_addition > 0) then
                    call addState_new(vals,&
                            coords,&
                            gradient)
                end if

                if ((any(vals > var_maxvar)).or.&
                    (steps > Nsteps)) exit
 
            end if

            Naccept = 0
            cycle

        else if (Nalpha_tries > Nalpha_tries_max) then
            call Acceleration(vals,coords,gradient)

            do i = 1, Ngrid_total
                write(filechannels(1+i),FMT=FMT6)&
                    default_rmsd
            end do

            Nalpha_tries = 1
        else
            gradient = approx_gradient
            Naccept = Naccept + 1
        end if

        if (gather_interpolation_flag) then

            call Acceleration(vals,coords,gradient)

            !If a frame was found, record data on it
            if (Ninterpolation > 0) then
                error1 = sqrt(sum((gradient - &
                        candidate_gradient)**2)/Natoms)
                error2 = sqrt(sum((gradient - &
                        approx_gradient)**2)/Natoms)
 
                write(filechannel3,FMT=*) vals(1),vals(2),&
                        Ninterpolation,largest_weighted_rmsd2,&
                        largest_weighted_rmsd,candidate_rmsd,&
                        min_rmsd,error1,error2
            end if
        end if

        !Update the velocities
        velocities = velocities + gradient

        !Upate the coordinates with the velocities
        coords = coords + dt * velocities

        !Update the steps
        steps = steps + 1

        !!!!!!!!!!!!!!!!!!!!
        ! TEST START
        !!!!!!!!!!!!!!!!!!!!
        call getEnergies(Natoms,coords,velocities,U,KE)
        E = U + KE

        !Finally write to the data file all the important data values
        write(filechannel2,FMT=*) number_of_frames,order,neighbor_check,steps,&
                                  min_rmsd,candidate_rmsd,vals(1),vals(2),U,KE
        !!!!!!!!!!!!!!!!!!!!
        ! TEST END
        !!!!!!!!!!!!!!!!!!!!


        if (Naccept == 0) then
            call getEnergies(Natoms,coords,velocities,U,KE)
            E = U + KE

            coords_prior = coords
            velocities_prior = velocities
            E_prior = E

            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
            if (grid_addition > 0) then
                call addState_new(vals,&
                        coords,&
                        gradient)
            end if

            if ((any(vals > var_maxvar)).or.&
                (steps > Nsteps)) exit
        end if

    end do

    if (gather_interpolation_flag) close(filechannel3)

    !Deallocate the buffers
!   deallocate(valsbuffer1,&
!           coordsbuffer1,gradientbuffer1,&
!           Ubuffer1,RMSDbuffer1,&
!           inputCLS)

!   deallocate(trajRMSDbuffer)

!   deallocate(CMdiffbuffer1)
!   deallocate(Ntrajbuffer1)

    call unsetAllocations()

    deallocate(trajRMSDbuffer)

    !Output the final coordinates and velocities
    coords_final = coords
    velocities_final = velocities



    close(filechannel2)

end subroutine runTrajectoryRewind1









subroutine runTrajectory_cap(filechannels,&
                             coords_initial,velocities_initial,&
                             coords_final,velocities_final)
    use PARAMETERS
    use ls_rmsd_original
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

    !Various other variables
    real(dp) :: error1,error2
    real(dp),dimension(3) :: x_center, y_center
    real(dp),dimension(3,3) ::  U_prior
    real(dp) :: min_rmsd,min_rmsd_prime
    integer :: number_of_frames,order,neighbor_check
    real(dp) :: U, KE

    real(dp) :: RMSD_prior, error_prior
    real(dp) :: max_RMSDdelta, max_deltaRMSD
    real(dp),dimension(3,Natoms) :: coords_prior, gradient_prior
    real(dp),dimension(3,Natoms) :: velocities_prior
    real(dp) :: DE, E, E_prior, E_baseline
    integer :: Nalpha_tries
    logical :: temp_reject_flag = .true.

    real(dp),allocatable :: libcoords(:,:,:), libgradients(:,:,:)
    integer,allocatable :: libNtraj(:)

    integer,dimension(Nalpha) :: alpha_flagging

    !Incremental Integer
    integer :: i,j,n,m

    !Initialize the scene
    call InitialSetup3(coords,velocities)
    Norder1 = 0
    Norder_total = 0
    subcellsearch_max = subcellsearch_max1
    allocate(trajRMSDbuffer(Ngrid_max,Naccept_max+1))

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
    if (gather_interpolation_flag) then
        open(filechannel3,file=&
            gridpath5//interpolationfile,position="append")
    end if

    !Allocate all buffers; initialize the buffer size to be
    !the maximum number of frames expected to be seen
    buffer1_size = 2 + Ninterpolation_max
    allocate(valsbuffer1(Nvar,buffer1_size),&
             coordsbuffer1(3,Natoms,buffer1_size),&
             gradientbuffer1(3,Natoms,buffer1_size),&
             Ubuffer1(3,3,buffer1_size),&
             RMSDbuffer1(buffer1_size),&
             inputCLS(Ncoords+buffer1_size,buffer1_size))

    RMSDbuffer1 = default_rmsd

    steps = 1
    E_baseline = 0.0d0
    do
        coords = coords + dt * velocities
        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
        call Acceleration(vals,coords,gradient)

        call getEnergies(Natoms,coords,velocities,U,KE)
        E_baseline = E_baseline + U + KE

        do i = 1, Ngrid_total
            write(filechannels(1+i),FMT=FMT6)&
                default_rmsd
        end do

        velocities = velocities + gradient

        steps = steps + 1
        if (steps > Nsteps_baseline) then
            coords_prior = coords
            velocities_prior = velocities - gradient
            gradient_prior = gradient

            E_prior = E

            exit
        end if
    end do
    E_baseline = E_baseline / Nsteps_baseline
    Naccept = 0
     
    do
        !Check every 500 steps to see if we are out-of-bounds
        if (modulo(steps,500) == 1) then
            if (any(vals > var_maxvar)) exit
        endif

        !Upate the coordinates with the velocities
        coords = coords + dt * velocities

        !Always calculate the variables before checking a frame or accelerating
        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

        !Check for a frame in the grid
        !Set the default value beforehand
        if (accept_worst) then
            min_rmsd = 0.0d0
        else
            min_rmsd = default_rmsd
        end if
        
        RMSDbuffer1 = min_rmsd

        call checkState_new_cap(vals,coords,&
                approx_gradient,min_rmsd,&
                filechannels,number_of_frames,&
                order,neighbor_check)

        call getEnergies(Natoms,coords,velocities,U,KE)
        E = U + KE

        if (temp_reject_flag) then
            call Acceleration(vals,coords,gradient)
 
            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
            if (grid_addition > 0) then
                call addState_new(vals,&
                        coords,&
                        gradient)
            end if

            temp_reject_flag = .false.

            if (Naccept > 0) then
                DE = abs(E-E_prior)/Naccept
                if ((abs(E-E_baseline) > E_threshold)&
                    .or.(DE > DE_threshold)) then
    
                    coords = coords_prior
                    velocities = velocities_prior
                    gradient = gradient_prior
    
                    steps = steps - Naccept - 1

                    Nalpha_tries = Nalpha_tries + 1

                    if (Nalpha_tries > Nalpha_tries_max) then
                        Nalpha_tries = 1
                        temp_reject_flag = .true.
                    end if
                else
                    do i = 1, Ngrid_total
                    do j = 1, Naccept+1
                        write(filechannels(1+i),FMT=FMT6)&
                            trajRMSDbuffer(i,j)
                    end do
                    end do
                end if

                Naccept = 0
            else
                Nalpha_tries = 1

                do i = 1, Ngrid_total
                    write(filechannels(1+i),FMT=FMT6)&
                        default_rmsd
                end do
            end if

            alpha_ratio = alpha_ratio_list(Nalpha_tries)
        
        !If we are gather interpolation data...
        else if (gather_interpolation_flag) then
 
            !Calculate the true gradient
            call Acceleration(vals,coords,gradient)

            !If a frame was found, record data on it
            if (Ninterpolation > 0) then
                error1 = sqrt(sum((gradient - &
                        candidate_gradient)**2)/Natoms)
                error2 = sqrt(sum((gradient - &
                        approx_gradient)**2)/Natoms)
 
                write(filechannel3,FMT=*) vals(1),vals(2),&
                        Ninterpolation,largest_weighted_rmsd2,&
                        largest_weighted_rmsd,candidate_rmsd,&
                        min_rmsd,error1,error2
            end if

            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
            if (grid_addition > 0) then
                call addState_new(vals,&
                        coords,gradient)
            end if
 
            !If the approximated frame is good enough
            !and we are not rejecting it, then use it
            if ((min_rmsd .ge. threshold_RMSD) &
                    .or. (reject_flag)) then

                if (Naccept > 0) then
                    DE = abs(E-E_prior)/Naccept
                    if ((abs(E-E_baseline) > E_threshold)&
                        .or.(DE > DE_threshold)) then
        
                        coords = coords_prior
                        velocities = velocities_prior
                        gradient = gradient_prior
        
                        steps = steps - Naccept - 1
    
                        Nalpha_tries = Nalpha_tries + 1
    
                        if (Nalpha_tries > Nalpha_tries_max) then
                            Nalpha_tries = 1
                            temp_reject_flag = .true.
                        end if
                    else
                        do i = 1, Ngrid_total
                        do j = 1, Naccept+1
                            write(filechannels(1+i),FMT=FMT6)&
                                trajRMSDbuffer(i,j)
                        end do
                        end do
                    end if

                    Naccept = 0
                else
                    Nalpha_tries = 1

                    do i = 1, Ngrid_total
                        write(filechannels(1+i),FMT=FMT6)&
                            trajRMSDbuffer(i,1)
                    end do
                end if

                alpha_ratio = alpha_ratio_list(Nalpha_tries)
            else
                gradient = approx_gradient

                Naccept = Naccept + 1
                if (Naccept == Naccept_max) temp_reject_flag = .true.
            end if
 
        !If no interpolation data is being gathered, then the
        !true gradient is not needed
 
        !If the frame is not good enough or we are rejecting
        !frames, calculate the true gradient
        else if ((min_rmsd .ge. threshold_RMSD) &
                .or. (reject_flag)) then
            call Acceleration(vals,coords,gradient)
 
            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
            if (grid_addition > 0) then
                call addState_new(vals,&
                        coords,&
                        gradient)
            end if

            if (Naccept > 0) then
                DE = abs(E-E_prior)/Naccept
                if ((abs(E-E_baseline) > E_threshold)&
                    .or.(DE > DE_threshold)) then
    
                    coords = coords_prior
                    velocities = velocities_prior
                    gradient = gradient_prior
    
                    steps = steps - Naccept - 1

                    Nalpha_tries = Nalpha_tries + 1

                    if (Nalpha_tries > Nalpha_tries_max) then
                        Nalpha_tries = 1
                        temp_reject_flag = .true.
                    end if
                else
                    do i = 1, Ngrid_total
                    do j = 1, Naccept+1
                        write(filechannels(1+i),FMT=FMT6)&
                            trajRMSDbuffer(i,j)
                    end do
                    end do
                end if

                Naccept = 0
            else
                Nalpha_tries = 1

                do i = 1, Ngrid_total
                    write(filechannels(1+i),FMT=FMT6)&
                        trajRMSDbuffer(i,1)
                end do
            end if

            alpha_ratio = alpha_ratio_list(Nalpha_tries)

        !If the approximated frame is good enough then
        !we can use it (after relabelling)!
        else
            gradient = approx_gradient

            Naccept = Naccept + 1
            if (Naccept == Naccept_max) temp_reject_flag = .true.
        end if

        if (Naccept == 0) then
            coords_prior = coords
            velocities_prior = velocities
            gradient_prior = gradient

            E_prior = E
        end if

        !Update the velocities
        velocities = velocities + gradient

        steps = steps + 1
        if (steps > Nsteps) exit
    end do

    if (gather_interpolation_flag) then
        close(filechannel3)
    end if

    !Deallocate the buffers
    deallocate(valsbuffer1,&
            coordsbuffer1,gradientbuffer1,&
            Ubuffer1,RMSDbuffer1,&
            inputCLS)

    deallocate(trajRMSDbuffer)

    !Output the final coordinates and velocities
    coords_final = coords
    velocities_final = velocities

end subroutine runTrajectory_cap

subroutine runTrajectory_permute_cap(filechannels,&
                                     coords_initial,velocities_initial,&
                                     coords_final,velocities_final)
    use PARAMETERS
    use ls_rmsd_original
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

    !Various other variables
    real(dp) :: error1,error2
    real(dp),dimension(3) :: x_center, y_center
    real(dp),dimension(3,3) ::  U_prior
    real(dp) :: min_rmsd,min_rmsd_prime
    integer :: number_of_frames,order,neighbor_check
    real(dp) :: U, KE

    real(dp) :: RMSD_prior, error_prior
    real(dp) :: max_RMSDdelta, max_deltaRMSD
    real(dp),dimension(3,Natoms) :: coords_prior, gradient_prior
    real(dp),dimension(3,Natoms) :: velocities_prior
    real(dp) :: DE, E, E_prior, E_baseline
    integer :: Nalpha_tries
    logical :: temp_reject_flag = .true.

    real(dp),allocatable :: libcoords(:,:,:), libgradients(:,:,:)
    integer,allocatable :: libNtraj(:)

    integer,dimension(Nalpha) :: alpha_flagging

    !Incremental Integer
    integer :: i,j,n,m

    !Initialize the scene
    call InitialSetup3(coords,velocities)
    Norder1 = 0
    Norder_total = 0
    subcellsearch_max = subcellsearch_max1
    allocate(trajRMSDbuffer(Ngrid_max,Naccept_max+1))

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
    if (gather_interpolation_flag) then
        open(filechannel3,file=&
            gridpath5//interpolationfile,position="append")
    end if

    !Allocate all buffers; initialize the buffer size to be
    !the maximum number of frames expected to be seen
    buffer1_size = 2 + Ninterpolation_max
    allocate(valsbuffer1(Nvar,buffer1_size),&
             coordsbuffer1(3,Natoms,buffer1_size),&
             gradientbuffer1(3,Natoms,buffer1_size),&
             Ubuffer1(3,3,buffer1_size),&
             RMSDbuffer1(buffer1_size),&
             inputCLS(Ncoords+buffer1_size,buffer1_size))

    RMSDbuffer1 = default_rmsd

    steps = 1
    E_baseline = 0.0d0
    do
        coords = coords + dt * velocities
        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
        call Acceleration(vals,coords,gradient)

        call getEnergies(Natoms,coords,velocities,U,KE)
        E_baseline = E_baseline + U + KE

        do i = 1, Ngrid_total
            write(filechannels(1+i),FMT=FMT6)&
                default_rmsd
        end do

        velocities = velocities + gradient

        steps = steps + 1
        if (steps > Nsteps_baseline) then
            coords_prior = coords
            velocities_prior = velocities - gradient
            gradient_prior = gradient

            E_prior = E

            exit
        end if
    end do
    E_baseline = E_baseline / Nsteps_baseline
    Naccept = 0
     
    do
        !Check every 500 steps to see if we are out-of-bounds
        if (modulo(steps,500) == 1) then
            if (any(vals > var_maxvar)) exit
        endif

        !Upate the coordinates with the velocities
        coords = coords + dt * velocities

        !Always calculate the variables before checking a frame or accelerating
        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

        !Check for a frame in the grid
        !Set the default value beforehand
        if (accept_worst) then
            min_rmsd = 0.0d0
        else
            min_rmsd = default_rmsd
        end if
        
        RMSDbuffer1 = min_rmsd

        call checkState_new_permute_cap(vals,coords,&
                approx_gradient,min_rmsd,&
                filechannels,number_of_frames,&
                order,neighbor_check)

        call getEnergies(Natoms,coords,velocities,U,KE)
        E = U + KE

        if (temp_reject_flag) then
            call Acceleration(vals,coords,gradient)
 
            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
            if (grid_addition > 0) then
                call addState_new(vals,&
                        coords,&
                        gradient)
            end if

            temp_reject_flag = .false.

            if (Naccept > 0) then
                DE = abs(E-E_prior)/Naccept
                if ((abs(E-E_baseline) > E_threshold)&
                    .or.(DE > DE_threshold)) then
    
                    coords = coords_prior
                    velocities = velocities_prior
                    gradient = gradient_prior
    
                    steps = steps - Naccept - 1

                    Nalpha_tries = Nalpha_tries + 1

                    if (Nalpha_tries > Nalpha_tries_max) then
                        Nalpha_tries = 1
                        temp_reject_flag = .true.
                    end if
                else
                    do i = 1, Ngrid_total
                    do j = 1, Naccept+1
                        write(filechannels(1+i),FMT=FMT6)&
                            trajRMSDbuffer(i,j)
                    end do
                    end do
                end if

                Naccept = 0
            else
                Nalpha_tries = 1

                do i = 1, Ngrid_total
                    write(filechannels(1+i),FMT=FMT6)&
                        default_rmsd
                end do
            end if

            alpha_ratio = alpha_ratio_list(Nalpha_tries)
        
        !If we are gather interpolation data...
        else if (gather_interpolation_flag) then
 
            !Calculate the true gradient
            call Acceleration(vals,coords,gradient)

            !If a frame was found, record data on it
            if (Ninterpolation > 0) then
                error1 = sqrt(sum((gradient - &
                        candidate_gradient)**2)/Natoms)
                error2 = sqrt(sum((gradient - &
                        approx_gradient)**2)/Natoms)
 
                write(filechannel3,FMT=*) vals(1),vals(2),&
                        Ninterpolation,largest_weighted_rmsd2,&
                        largest_weighted_rmsd,candidate_rmsd,&
                        min_rmsd,error1,error2
            end if

            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
            if (grid_addition > 0) then
                call addState_new(vals,&
                        coords,gradient)
            end if
 
            !If the approximated frame is good enough
            !and we are not rejecting it, then use it
            if ((min_rmsd .ge. threshold_RMSD) &
                    .or. (reject_flag)) then

                if (Naccept > 0) then
                    DE = abs(E-E_prior)/Naccept
                    if ((abs(E-E_baseline) > E_threshold)&
                        .or.(DE > DE_threshold)) then
        
                        coords = coords_prior
                        velocities = velocities_prior
                        gradient = gradient_prior
        
                        steps = steps - Naccept - 1
    
                        Nalpha_tries = Nalpha_tries + 1
    
                        if (Nalpha_tries > Nalpha_tries_max) then
                            Nalpha_tries = 1
                            temp_reject_flag = .true.
                        end if
                    else
                        do i = 1, Ngrid_total
                        do j = 1, Naccept+1
                            write(filechannels(1+i),FMT=FMT6)&
                                trajRMSDbuffer(i,j)
                        end do
                        end do
                    end if

                    Naccept = 0
                else
                    Nalpha_tries = 1

                    do i = 1, Ngrid_total
                        write(filechannels(1+i),FMT=FMT6)&
                            trajRMSDbuffer(i,1)
                    end do
                end if

                alpha_ratio = alpha_ratio_list(Nalpha_tries)
            else
                gradient = approx_gradient

                Naccept = Naccept + 1
                if (Naccept == Naccept_max) temp_reject_flag = .true.
            end if
 
        !If no interpolation data is being gathered, then the
        !true gradient is not needed
 
        !If the frame is not good enough or we are rejecting
        !frames, calculate the true gradient
        else if ((min_rmsd .ge. threshold_RMSD) &
                .or. (reject_flag)) then
            call Acceleration(vals,coords,gradient)
 
            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
            if (grid_addition > 0) then
                call addState_new(vals,&
                        coords,&
                        gradient)
            end if

            if (Naccept > 0) then
                DE = abs(E-E_prior)/Naccept
                if ((abs(E-E_baseline) > E_threshold)&
                    .or.(DE > DE_threshold)) then
    
                    coords = coords_prior
                    velocities = velocities_prior
                    gradient = gradient_prior
    
                    steps = steps - Naccept - 1

                    Nalpha_tries = Nalpha_tries + 1

                    if (Nalpha_tries > Nalpha_tries_max) then
                        Nalpha_tries = 1
                        temp_reject_flag = .true.
                    end if
                else
                    do i = 1, Ngrid_total
                    do j = 1, Naccept+1
                        write(filechannels(1+i),FMT=FMT6)&
                            trajRMSDbuffer(i,j)
                    end do
                    end do
                end if

                Naccept = 0
            else
                Nalpha_tries = 1

                do i = 1, Ngrid_total
                    write(filechannels(1+i),FMT=FMT6)&
                        trajRMSDbuffer(i,1)
                end do
            end if

            alpha_ratio = alpha_ratio_list(Nalpha_tries)

        !If the approximated frame is good enough then
        !we can use it (after relabelling)!
        else
            gradient = approx_gradient

            Naccept = Naccept + 1
            if (Naccept == Naccept_max) temp_reject_flag = .true.
        end if

        if (Naccept == 0) then
            coords_prior = coords
            velocities_prior = velocities
            gradient_prior = gradient

            E_prior = E
        end if

        !Update the velocities
        velocities = velocities + gradient

        steps = steps + 1
        if (steps > Nsteps) exit
    end do

    if (gather_interpolation_flag) then
        close(filechannel3)
    end if

    !Deallocate the buffers
    deallocate(valsbuffer1,&
            coordsbuffer1,gradientbuffer1,&
            Ubuffer1,RMSDbuffer1,&
            inputCLS)

    deallocate(trajRMSDbuffer)

    !Output the final coordinates and velocities
    coords_final = coords
    velocities_final = velocities

end subroutine runTrajectory_permute_cap







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
!		Every frame is checked with the current grid through interactMultipleGrids
!
!		The initial and final velocity of the incoming H/H2 is output to
!		gather scattering angle data
!
!               Every interpolated frame is output to gather interpolation data
!               according to ANALYSIS
!
!               Every frame whose gradient is calculated may be added to the library
!               according to PARAMETERS and grid_addition in ANALYSIS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	INPUT				KIND				DESCRIPTION
!
!		filechannels			INTEGER,DIM(Ngrid_max+1)	These filechannels were opened beforehand to the
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
!               buffer1_size                    INTEGER                         The size of the buffer for coordinates and
!                                                                               gradients to be interpolated; the buffer
!                                                                               size is dynamically increased
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

    !Allocate all buffers; initialize the buffer size to be
    !the maximum number of frames expected to be seen
    buffer1_size = 1 + var_overcrowd(1)
    allocate(valsbuffer1(Nvar,buffer1_size),&
             coordsbuffer1(3,Natoms,buffer1_size),&
             gradientbuffer1(3,Natoms,buffer1_size),&
             Ubuffer1(3,3,buffer1_size),&
             RMSDbuffer1(buffer1_size),&
             approximation_index(buffer1_size))

    allocate(Ntrajbuffer1(buffer1_size))

    !These buffers need not be allocated if interpolation
    !is not occuring
    allocate(acceptable_frame_mask(buffer1_size),&
             inputCLS(Ncoords+buffer1_size,buffer1_size))

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
        !(deprecated)
        if ((force_Duplicates) .or. (force_NoLabels)) then

            if (accept_worst) then
                min_rmsd = 0.0d0
            else
                min_rmsd = default_rmsd
            end if

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

        !Otherwise, we need to do the label switching
        else

            !We need to input the frame using the labeling scheme
            do n = 1, Natoms
                coords_labelled(:,n) = coords(:,BOND_LABELLING_DATA(n))
            end do

            !Check for a frame in the grid
            !Set the default value beforehand
            if (accept_worst) then
                min_rmsd = 0.0d0
            else
                min_rmsd = default_rmsd
            end if
 
            !If we want detailed RMSD output then we
            !search the library twice; the first time
            !with no interpolation
            if (testtrajDetailedRMSD_flag) then
                subcellsearch_max = subcellsearch_max1
                if (gather_interpolation_flag) then 
                    interpolation_flag = .false.
                end if
            end if

            !Check the library to get an
            !approximated gradient
            call checkState_new(vals,coords_labelled,approx_gradient,min_rmsd,&
                     filechannels,number_of_frames,order,neighbor_check)

!Here is the second search; not done otherwise
if (testtrajDetailedRMSD_flag) then

    !Set a separate default value
    if (accept_worst) then
        min_rmsd_prime = 0.0d0
    else
        min_rmsd_prime = default_rmsd
    end if

    !The second search uses interpolation as it
    !usually is a wider search
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

    !The better frame is usually found with second search
    !so this will be used for further data handling
    min_rmsd = min_rmsd_prime
else
    !A deault value is set to approx_gradient_prime for
    !further data handling
    approx_gradient_prime = approx_gradient
end if
           !Accept worst has not found an acceptable frame
           !if min_rmsd is zero
           !(deprecated)
           if ((accept_worst) .and. (min_rmsd == 0.0d0)) then
               call Acceleration(vals,coords,gradient)
           
           !If we are gather interpolation data...
           else if (gather_interpolation_flag) then

               !Save the approx_gradient with the
               !correct labelling
               do n = 1, Natoms
                    gradient(:,BOND_LABELLING_DATA(n)) = &
                            approx_gradient(:,n)
               end do
               approx_gradient = gradient

               !Save the approx_gradient_prime
               !with the correct labelling
               do n = 1, Natoms
                    gradient(:,BOND_LABELLING_DATA(n)) = &
                            approx_gradient_prime(:,n)
               end do
               approx_gradient_prime = gradient

               !Always calculate the true gradient
               call Acceleration(vals,coords,gradient)

               !If a frame was found, record data on it
               if (Ninterpolation > 0) then
               write(filechannel3,FMT=*) vals(1), vals(2), Ninterpolation,&
                           largest_weighted_rmsd2, largest_weighted_rmsd,&
                           min_rmsd,&
                           sqrt(sum((gradient-approx_gradient)**2)/Natoms),&
                           min_rmsd_prime,&
                           sqrt(sum((gradient-approx_gradient_prime)**2)/Natoms)
               end if

               !If we are adding to a grid, then relabel
               !the coordinates and gradients, and do so
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

               !If the approximated frame is good enough
               !and we are not rejecting it, then use it
               if ((min_rmsd .ge. threshold_RMSD).or.(reject_flag)) then
               else
                       gradient = approx_gradient
               end if

           !If no interpolation data is being gathered, then the
           !true gradient is not needed

           !If the frame is not good enough or we are rejecting
           !frames, calculate the true gradient
           else if ((min_rmsd .ge. threshold_RMSD).or.(reject_flag)) then
               call Acceleration(vals,coords,gradient)

               !If we are adding to a grid, then relabel
               !the coordinates and gradients, and do so
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

           !If the approximated frame is good enough then
           !we can use it (after relabelling)!
           else
               do n = 1, Natoms
                    gradient(:,BOND_LABELLING_DATA(n)) = &
                            approx_gradient(:,n)
               end do
           end if
        end if

        !Update the velocities
        velocities = velocities + gradient
    end do

    if (testtrajDetailedRMSD_flag) close(filechannel2)
    if (gather_interpolation_flag) close(filechannel3)

    !Deallocate the buffers
    deallocate(valsbuffer1,&
            coordsbuffer1,gradientbuffer1,&
            Ubuffer1,RMSDbuffer1,&
            approximation_index)

    deallocate(Ntrajbuffer1)

    deallocate(acceptable_frame_mask,inputCLS)

    !Output the final coordinates and velocities
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
!		Every frame is checked with the current grid through interactMultipleGrids
!
!		The initial and final velocity of the incoming H/H2 is output to
!		gather scattering angle data
!
!               Every interpolated frame is output to gather interpolation data
!               according to ANALYSIS
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
    use ls_rmsd_original
    use analyzeRMSDThresholdwithMultipleGrids
    implicit none

    !Coordinates, Velocities, and Variables
    real(dp),dimension(3,Natoms) :: coords,coords_labelled,velocities
    real(dp),dimension(3,Natoms) :: gradient,gradient_labelled,approx_gradient
    real(dp),dimension(3,Natoms),intent(out) :: coords_initial, velocities_initial
    real(dp),dimension(3,Natoms),intent(out) :: coords_final, velocities_final
    real(dp),dimension(Nvar) :: vals
    
    real(dp),dimension(3,Natoms) :: coords_mean
    real(dp),dimension(3) :: x_center,y_center
    real(dp) :: mu,sigma

    !Grid Parameters
    integer,dimension(1+Ngrid_total),intent(in) :: filechannels
    integer :: OMP_GET_THREAD_NUM

    !Various other variables
    real(dp) :: U, KE
    real(dp) :: min_rmsd, min_rmsd_prime
    real(dp) :: error1, error2
    integer :: number_of_frames,order,neighbor_check
    integer :: order0, order1

    real(dp),allocatable :: libcoords(:,:,:), libgradients(:,:,:)
    real(dp),dimension(3,Natoms) :: coords2,gradient2
    integer,allocatable :: libNtraj(:)

    !Incremental Integer
    integer :: n,m

    !Initialize the scene
    call InitialSetup3(coords,velocities)

    !For traversal
    if (traversal_flag) then
        traversal0 = 0
        traversal1 = 0
    end if

    !Open these files and keep them open if we want this
    !type of data
    if (testtrajDetailedRMSD_flag) open(filechannel2,file=&
            gridpath5//checkstatefile)
    if (gather_interpolation_flag) open(filechannel3,file=&
            gridpath5//interpolationfile,position="append")

    !Allocate all buffers; initialize the buffer size to be
    !the maximum number of frames expected to be seen
    buffer1_size = 1 + var_overcrowd(1)
    allocate(valsbuffer1(Nvar,buffer1_size),&
             coordsbuffer1(3,Natoms,buffer1_size),&
             gradientbuffer1(3,Natoms,buffer1_size),&
             Ubuffer1(3,3,buffer1_size),&
             RMSDbuffer1(buffer1_size),&
             approximation_index(buffer1_size))

    allocate(Ntrajbuffer1(buffer1_size))

    allocate(acceptable_frame_mask(buffer1_size),&
             inputCLS(Ncoords+buffer1_size,buffer1_size))

    !Output the initial coordinates and velocities
    coords_initial = coords
    velocities_initial = velocities

    !Always calculate the variables before accelearting
    call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

    do n = 1, Natoms
        coords2(:,n) = coords(:,BOND_LABELLING_DATA(n))
    end do
    
!   print *, ""
!   print *, 0
!   print *, vals
!   print *, BOND_LABELLING_DATA
!   print *, 0.0d0

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
!       if (modulo(steps,50) == 1) then
        if (.false.) then
            open(filechannel1,file=gridpath5//trajectoryfile,position="append")
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
        !(deprecated)
        if ((force_Duplicates) .or. (force_NoLabels)) then
            if (accept_worst) then
                min_rmsd = 0.0d0
            else
                min_rmsd = default_rmsd
            end if
            call checkState_new(vals,coords,approx_gradient,min_rmsd,&
                     filechannels,number_of_frames,order,neighbor_check)

            if ((accept_worst) .and. (min_rmsd == 0.0d0)) then
                call Acceleration(vals,coords,gradient)
            else if ((min_rmsd .ge. threshold_RMSD).or.(reject_flag)) then
                call Acceleration(vals,coords,gradient)
            else
                gradient = approx_gradient
            end if

        !Otherwise, then we do the label switching
        else
            !We need to input the frame using the labeling scheme
            do n = 1, Natoms
                coords_labelled(:,n) = coords(:,BOND_LABELLING_DATA(n))
            end do

            call rmsd_dp(Natoms,coords_labelled,coords2,&
                 1,Ubuffer1(:,:,1),x_center,y_center,&
                 min_rmsd)

            print *, ""
            print *, steps
            print *, vals
            print *, BOND_LABELLING_DATA
            print *, min_rmsd

            !Check for a frame in the grid
            !Set the default value beforehand
            if (accept_worst) then
                min_rmsd = 0.0d0
            else
                min_rmsd = default_rmsd
            end if

            !If we want detailed RMSD output then we
            !search the library twice; the first time
            !with no interpolation
            if (testtrajDetailedRMSD_flag) then
                subcellsearch_max = subcellsearch_max1
                if (gather_interpolation_flag) then 
                    interpolation_flag = .false.
                end if
            end if

            !Check the library to get an
            !approximated gradient
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

           !Accept worst has not found an acceptable frame
           !if min_rmsd is zero
           !(deprecated)
           if ((accept_worst) .and. (min_rmsd == 0.0d0)) then
               call Acceleration(vals,coords,gradient)

           !If we are gather interpolation data...
           else if (gather_interpolation_flag) then

               !Save the approx_gradient with the
               !correct labelling
               do n = 1, Natoms
                    gradient(:,BOND_LABELLING_DATA(n)) = &
                            approx_gradient(:,n)
               end do
               approx_gradient = gradient

               !Save the candidate_gradient with the
               !correct labelling
               do n = 1, Natoms
                    gradient(:,BOND_LABELLING_DATA(n)) = &
                            candidate_gradient(:,n)
               end do
               candidate_gradient = gradient

               !Calculate the true gradient
               call Acceleration(vals,coords,gradient)

               !If a frame was found, record data on it
               if (Ninterpolation > 0) then

               error1 = sqrt(sum((gradient-candidate_gradient)**2)/Natoms)
               error2 = sqrt(sum((gradient-approx_gradient)**2)/Natoms)

               write(filechannel3,FMT=*) vals(1), vals(2), Ninterpolation,&
                           largest_weighted_rmsd2, largest_weighted_rmsd,&
                           candidate_rmsd,min_rmsd,error1,error2

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               if ((Ninterpolation > 10).and.(error1 < error2)) then
                   allocate(libcoords(Ninterpolation,3,Natoms),&
                            libgradients(Ninterpolation,3,Natoms),&
                            libNtraj(Ninterpolation))
                   m = 0
!                  coords_mean = 0.0d0
                   do n = 1, Ninterpolation
                       do
                       m = m + 1
                       if (acceptable_frame_mask(m)) then
                           libcoords(n,:,:) =&
                                   coordsbuffer1(:,:,m)
                           libgradients(n,:,:) =&
                                   gradientbuffer1(:,:,m)
                           libNtraj(n) =&
                                   Ntrajbuffer1(m)

!                          coords_mean = coords_mean + &
!                                  libcoords(n,:,:)
                           exit
                       else
                       end if
                       end do
                   end do

!                  coords_mean = coords_mean / &
!                          Ninterpolation
!                  mu = sqrt(sum(coords_mean**2)/Natoms)
!                  
!                  sigma = 0.0d0
!                  do n = 1, Ninterpolation
!                      sigma = sigma + &
!                          sum((coords_mean-libcoords(n,:,:))**2)
!                  end do
!                  sigma = sqrt(sigma/Ninterpolation)

!                  do n = 1, Natoms
!                       gradient_labelled(:,n) = &
!                               gradient(:,BOND_LABELLING_DATA(n))
!                  end do

!   open(6666,file=gridpath5//"tmp_A.dat",&
!   position="append")
!   write(6666,FMT=*) vals(1), vals(2), mu, sigma, &
!                     coords_labelled, coords_mean
!   write(6666,FMT=*) mu, sigma, &
!                     coords_labelled, coords_mean
!   close(6666)
               
                   !errorCheck3 requires LABELLED data
!                  call errorCheck3(filechannels,&
!                          coords_labelled,&
!                          gradient_labelled,&
!                          Ninterpolation,&
!                          libcoords,libgradients)
!
!                   coords2 = coords + dt * &
!                           (velocities + gradient)
!                   call Acceleration(vals,coords2,gradient2)
!
!                   !errorCheck4 requires UNLABELLED data
!                   call errorCheck4(filechannels,&
!                           coords,gradient,&
!                           coords2,gradient2,&
!                           Ninterpolation/2,&
!                           Ninterpolation,&
!                           libcoords,libgradients)
!
                   !errorCheck5 requires UNLABELLED data
                   call errorCheck5(filechannels,&
                           coords,gradient,&
                           Ninterpolation,&
                           libcoords,libgradients,&
                           libNtraj)

                   deallocate(libcoords,libgradients,libNtraj)
               end if

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               end if

               !If the approximated frame is good enough
               !and we are not rejecting it, then use it
               if ((min_rmsd .ge. threshold_RMSD).or.(reject_flag)) then
               else
                   gradient = approx_gradient
               end if

           !If no interpolation data is being gathered, then the
           !true gradient is not needed

           !If the frame is not good enough or we are rejecting
           !frames, calculate the true gradient
           else if ((min_rmsd .ge. threshold_RMSD).or.(reject_flag)) then
               call Acceleration(vals,coords,gradient)

           !If the approximated frame is good enough then
           !we can use it (after relabelling)!
           else
               do n = 1, Natoms
                    gradient(:,BOND_LABELLING_DATA(n)) = &
                            approx_gradient(:,n)
               end do
           end if
        end if

        !Update the velocities
        velocities = velocities + gradient

    end do

    if (testtrajDetailedRMSD_flag) close(filechannel2)
    if (gather_interpolation_flag) close(filechannel3)

    !Deallocate the buffers
    deallocate(valsbuffer1,&
            coordsbuffer1,gradientbuffer1,&
            Ubuffer1,RMSDbuffer1,&
            approximation_index)

    deallocate(Ntrajbuffer1)

    deallocate(acceptable_frame_mask,inputCLS)

    !Output the final coordinates and velocities
    coords_final = coords
    velocities_final = velocities

end subroutine checkMultipleTrajectories







subroutine readTrajectory(filechannel_input,filechannels,&
                          coords_initial,velocities_initial,&
                          coords_final,velocities_final)
    use PHYSICS
    use PARAMETERS
    use VARIABLES
    use ANALYSIS
    use interactMultipleGrids
    use ls_rmsd_original
    use analyzeRMSDThresholdwithMultipleGrids
    implicit none

    !Coordinates, Velocities, and Variables
    real(dp),dimension(3,Natoms) :: coords,coords_labelled,velocities
    real(dp),dimension(3,Natoms) :: gradient,gradient_labelled,approx_gradient
    real(dp),dimension(3,Natoms),intent(out) :: coords_initial, velocities_initial
    real(dp),dimension(3,Natoms),intent(out) :: coords_final, velocities_final
    real(dp),dimension(Nvar) :: vals, vals_previous
    
    real(dp),dimension(3,Natoms) :: coords_mean
    real(dp),dimension(3) :: x_center,y_center
    real(dp) :: mu,sigma

    !Grid Parameters
    integer,dimension(1+Ngrid_total),intent(in) :: filechannels
    integer,intent(in) :: filechannel_input
    integer :: OMP_GET_THREAD_NUM

    !Various other variables
    real(dp) :: U, KE
    real(dp) :: min_rmsd, min_rmsd_prime
    real(dp) :: error1, error2
    integer :: number_of_frames,order,neighbor_check
    integer :: order0, order1
    integer :: iostate

    real(dp) :: RMSD_prior, error_prior
    real(dp) :: max_RMSDdelta, max_deltaRMSD
    real(dp),dimension(3,Natoms) :: coords_prior, gradient_prior
    real(dp),dimension(3,Natoms) :: velocities_prior
    real(dp) :: DE, E, E_prior, E_baseline
    integer :: Nalpha_tries = 1
    logical :: temp_reject_flag = .false.

    real(dp),allocatable :: libcoords(:,:,:), libgradients(:,:,:)
    real(dp),dimension(3,Natoms) :: coords2,gradient2
    integer,allocatable :: libNtraj(:)
    integer,dimension(Nalpha) :: alpha_flagging

    !Incremental Integer
    integer :: i,j,n,m

    allocate(trajRMSDbuffer(Ngrid_max,Naccept_max+1))

    !For traversal
    if (traversal_flag) then
        traversal0 = 0
        traversal1 = 0
    end if

    !Open these files and keep them open if we want this
    !type of data
    if (testtrajDetailedRMSD_flag) open(filechannel2,file=&
            gridpath5//checkstatefile)
    if (gather_interpolation_flag) open(filechannel3,file=&
            gridpath5//interpolationfile,position="append")

    !Allocate all buffers; initialize the buffer size to be
    !the maximum number of frames expected to be seen
    buffer1_size = 2 + Ninterpolation_max
    buffer2_size = 30

    call setAllocations()
!   allocate(valsbuffer1(Nvar,buffer1_size),&
!            coordsbuffer1(3,Natoms,buffer1_size),&
!            gradientbuffer1(3,Natoms,buffer1_size),&
!            Ubuffer1(3,3,buffer1_size),&
!            RMSDbuffer1(buffer1_size),&
!            approximation_index(buffer1_size))

!   allocate(CMdiffbuffer1(buffer1_size))
!   allocate(Ntrajbuffer1(buffer1_size))

!   allocate(acceptable_frame_mask(buffer1_size),&
!            inputCLS(Ncoords+buffer1_size,buffer1_size))

    read(filechannel_input,iostat=iostate) coords
    if (iostate /= 0) then
        print *, "uh oh: bad input"
        return
    end if
    read(filechannel_input) KE, U, E
    read(filechannel_input) gradient

    !Output the initial coordinates and velocities
    coords_initial = coords
    velocities_initial = velocities

    !Always calculate the variables before accelearting
    call getVarsHBrCO2(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

    populationbuffer2 = -1
    do i = 1, Nvar
        previous_var_index(i) = int(vals(i)*divisor(i,Norder_order(1)+1))
    end do 

    !Start the main loop
    steps = 2
    Naccept = 0
    do
    
        read(filechannel_input,iostat=iostate) coords
        if (iostate /= 0) exit
        read(filechannel_input) KE, U, E
        read(filechannel_input) gradient

        !Just for bug-testing
!       if (modulo(steps,50) == 1) then
!       if (.false.) then
        if (.true.) then
            open(filechannel1,file=gridpath5//trajectoryfile,position="append")
!           write(filechannel1,'(I6)') Natoms
!           write(filechannel1,*) ""
!           do n = 1, Natoms
!               write(filechannel1,'(A1,3F10.6)') 'H',&
!                             coords(1,n), coords(2,n), coords(3,n)
!           end do
            write(filechannel1,'(2F10.6)') vals
            close(filechannel1)
        end if

        !Always calculate the variables before accelearting
        call getVarsHBrCO2(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

        !Check for a frame in the grid
        !Set the default value beforehand
        if (accept_worst) then
            min_rmsd = 0.0d0
        else
            min_rmsd = default_rmsd
        end if
        
        RMSDbuffer1 = min_rmsd
        CMdiffbuffer1 = default_CMdiff

!       call checkState_new_permute_cap(&
!               vals,coords,&
!               approx_gradient,min_rmsd,&
!               filechannels,number_of_frames,&
!               order,neighbor_check)

        call checkState_PCM(&
                vals,coords,&
                approx_gradient,min_rmsd,&
                filechannels,number_of_frames,&
                order,neighbor_check)

        !If we are gather interpolation data...
        if (gather_interpolation_flag) then

            !If a frame was found, record data on it
!           if (Ninterpolation > 0) then
            do j = Ninterpolation, 1, -1
!           do j = min(Ninterpolation,1), 1, -1

                candidate_gradient = matmul(Ubuffer1(:,:,j),&
                        gradientbuffer1(:,:,j))
                candidate_rmsd = RMSDbuffer1(j)

                error1 = sqrt(sum((gradient - &
                        candidate_gradient)**2)/Natoms)
                error2 = sqrt(sum((gradient - &
                        approx_gradient)**2)/Natoms)

!               error1 = 0.0d0
!               error2 = 0.0d0
!               do i = 1, Natoms
!                   error1 = error1 + sum(((gradient(:,i)-&
!                           candidate_gradient(:,i))/&
!                           (masses(i)*(0.001/Na)/RU_mass))**2)
!                   error2 = error2 + sum(((gradient(:,i)-&
!                           approx_gradient(:,i))/&
!                           (masses(i)*(0.001/Na)/RU_mass))**2)
!               end do
!               error1 = sqrt(error1)*20*dt*(hartree/bohr)/(RU_energy/RU_length)
!               error2 = sqrt(error2)*20*dt*(hartree/bohr)/(RU_energy/RU_length)
 
                write(filechannel3,FMT=*) vals(1),vals(2),&
                        Ninterpolation,largest_weighted_rmsd2,&
                        largest_weighted_rmsd,candidate_rmsd,&
                        min_rmsd,CMdiffbuffer1(j),&
                        interpolated_CMdiff,error1,error2
!           end if
            end do

            !!!!!!!!!!!!!!!!!!!!!!
            ! alpha testing start
            !!!!!!!!!!!!!!!!!!!!!!
!           if (Ninterpolation > 0 ) then

!               ! Just to be reasonable:
!               Ninterpolation = min(Ninterpolation,20)

!               allocate(libcoords(Ninterpolation,3,Natoms),&
!                        libgradients(Ninterpolation,3,Natoms))

!               do n = 1, Ninterpolation
!                   libcoords(n,:,:) =&
!                           coordsbuffer1(:,:,n)
!                   libgradients(n,:,:) =&
!                           gradientbuffer1(:,:,n)
!               end do

!               call errorCheck8(filechannels,&
!                       coords,gradient,&
!                       Ninterpolation,&
!                       libcoords,libgradients,&
!                       alpha_flagging)

!               deallocate(libcoords,libgradients)
!           end if
            !!!!!!!!!!!!!!!!!!!!!!
            ! alpha testing end
            !!!!!!!!!!!!!!!!!!!!!!

            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
            if (grid_addition > 0) then
                call addState_new(vals,&
                        coords,gradient)
            end if
 
            !If the approximated frame is good enough
            !and we are not rejecting it, then use it
            if ((min_rmsd .ge. threshold_RMSD) &
                    .or. (reject_flag)) then

                if (Naccept > 0) then
                    DE = abs(E-E_prior)/Naccept
                    if ((abs(E-E_baseline) > E_threshold)&
                        .or.(DE > DE_threshold)) then
        
                        coords = coords_prior
                        velocities = velocities_prior
                        gradient = gradient_prior
        
                        steps = steps - Naccept - 1
    
                        Nalpha_tries = Nalpha_tries + 1
    
                        if (Nalpha_tries > Nalpha_tries_max) then
                            Nalpha_tries = 1
                            temp_reject_flag = .true.
                        end if
                    else
                        if (testtraj_flag) then
                        do i = 1, Ngrid_total
                        do j = 1, Naccept+1
                            write(filechannels(1+i),FMT=FMT6)&
                                trajRMSDbuffer(i,j)
                        end do
                        end do
                        end if
                    end if

                    Naccept = 0
                else
                    Nalpha_tries = 1

                    if (testtraj_flag) then
                    do i = 1, Ngrid_total
                        write(filechannels(1+i),FMT=FMT6)&
                            trajRMSDbuffer(i,1)
                    end do
                    end if
                end if

                alpha_ratio = alpha_ratio_list(Nalpha_tries)
            end if
 
        !If no interpolation data is being gathered, then the
        !true gradient is not needed
 
        !If the frame is not good enough or we are rejecting
        !frames, calculate the true gradient
        else
 
            !If we are adding to a grid, then relabel
            !the coordinates and gradients, and do so
            if (grid_addition > 0) then
                call addState_new(vals,&
                        coords,&
                        gradient)
            end if

            if (Naccept > 0) then
                DE = abs(E-E_prior)/Naccept
                if ((abs(E-E_baseline) > E_threshold)&
                    .or.(DE > DE_threshold)) then
    
                    coords = coords_prior
                    velocities = velocities_prior
                    gradient = gradient_prior
    
                    steps = steps - Naccept - 1

                    Nalpha_tries = Nalpha_tries + 1

                    if (Nalpha_tries > Nalpha_tries_max) then
                        Nalpha_tries = 1
                        temp_reject_flag = .true.
                    end if
                else
                    if (testtraj_flag) then
                    do i = 1, Ngrid_total
                    do j = 1, Naccept+1
                        write(filechannels(1+i),FMT=FMT6)&
                            trajRMSDbuffer(i,j)
                    end do
                    end do
                    end if
                end if

                Naccept = 0
            else
                Nalpha_tries = 1

                if (testtraj_flag) then
                do i = 1, Ngrid_total
                    write(filechannels(1+i),FMT=FMT6)&
                        trajRMSDbuffer(i,1)
                end do
                end if
            end if

            alpha_ratio = alpha_ratio_list(Nalpha_tries)
        end if

        steps = steps + 1
        if (steps > Nsteps) exit
    end do
    steps = steps - 1

    if (testtrajDetailedRMSD_flag) close(filechannel2)
    if (gather_interpolation_flag) close(filechannel3)

    !Deallocate the buffers
!   deallocate(valsbuffer1,&
!           coordsbuffer1,gradientbuffer1,&
!           Ubuffer1,RMSDbuffer1,&
!           approximation_index)

!   deallocate(Ntrajbuffer1)
!   deallocate(CMdiffbuffer1)

!   deallocate(acceptable_frame_mask,inputCLS)

    call unsetAllocations()

    !Output the final coordinates and velocities
    coords_final = coords
    velocities_final = velocities

    deallocate(trajRMSDbuffer)

end subroutine readTrajectory





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

subroutine setSubcellSearchMax()
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

integer :: i,j
integer :: single_index
integer,dimension(Nvar) :: var_index

do i = 1, min(Norder_max+1,ssm_length)
    subcellsearch_max(i) = ssm1(i)
    subcellsearch_max1(i) = ssm1(i)
    subcellsearch_max2(i) = ssm2(i)
end do

if (memory_flag) then
    
    !Initialize the memory buffer assuming
    !only subcellsearch_max1 and the first
    !order search are used
    j = subcellsearch_max(&
            Norder_order(1)+1)
    single_index_max = 0
    
    do i = 1, Nvar
    
        var_index = 0
    
        var_index(i) = j
        call getFlattened(Nvar,var_index,&
                single_index)
        single_index_max = max(single_index,&
                single_index_max)
    
        var_index(i) = -j
        call getFlattened(Nvar,var_index,&
                single_index)
        single_index_max = max(single_index,&
                single_index_max)
    
    end do

end if

return

end subroutine setSubcellSearchMax



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
    real(dp), dimension(3,Natoms) :: gradient,gradient_var
    real(dp), dimension(3,Natoms) :: gradient_labelled,gradient_var_labelled
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
    if (gather_interpolation_flag)&
            open(filechannel3,file=gridpath5//&
            interpolationfile,position="append")

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
    allocate(rmsd_x_interpolated(Ntest,Nsamples),rmsd_x(Ntest,Nsamples),&
             rmsd_fx_interpolated(Ntest,Nsamples),rmsd_fx(Ntest,Nsamples),&
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
    if (interpolation_flag)&
            deallocate(acceptable_frame_mask,inputCLS)

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

    integer :: Ntrials

    integer,dimension(1+Ngrid_total),intent(in) :: filechannels
    real(dp), dimension(3) :: x_center, y_center
    real(dp), allocatable :: g(:,:)
    real(dp),dimension(3,3) :: U
    real(dp),dimension(3,3) :: candidate_U

    !Incremental Integer
    integer :: i,n

    print *, ""
    print *, "Started Error Check 2!"
    print *, ""

    !Initialize the scene
    call InitialSetup3(coords,velocities)

    coords_initial = coords
    velocities_initial = velocities

print *, "coords:"
print *, coords

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
    Ntest = 50
    Nsamples = 10
    Ntrials = 0
    allocate(rmsd_x_interpolated(Ntest,Nsamples),rmsd_fx(Ntest,Nsamples),&
             rmsd_fx_interpolated(Ntest,Nsamples),rmsd_x(Ntest,Nsamples),&
             rmsd_weights(Ntest,Ntest,Nsamples),rmsd_x2_interpolated(Ntest,Nsamples))
    allocate(mean_weights(Ntest),mean_rmsd_fxs(Ntest))
    allocate(inputCLS(Ncoords+Ntest,Ntest),gradient_steps(3,Natoms,Ntest))
    rmsd_weights = 0.0d0

    open(filechannel3,file=gridpath5//"convergence"//errorcheckfile)




    !We can choose to start recording later in our
    !trajectory
    if (.true.) then
    do steps = 1, 3000
        !Upate the coordinates with the velocities
        coords = coords + dt * velocities

        !Always calculate the variables before checking a frame or accelerating
        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

        call Acceleration(vals,coords,gradient)

        !Update the velocities
        velocities = velocities + gradient
    end do
    end if

    !Start of the "trajectory"
    do

    Nanomaly_index = 0
!   call system("rm "//gridpath5//"weight1"//errorcheckfile)

    !We can choose not to move forward in the trajectory
    !if we are focused on only one point
    if (.false.) then
    do steps = 1, 1000
        !Upate the coordinates with the velocities
        coords = coords + dt * velocities

        !Always calculate the variables before checking a frame or accelerating
        call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

        call Acceleration(vals,coords,gradient)

        !Update the velocities
        velocities = velocities + gradient
    end do

    !In this case, we need a counter
    else
        Ntrials = Ntrials + 1
    end if

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

    handicap_rmsd = threshold_rmsd* 0.90d0

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
!       inputCLS(Ncoords+1,1) = alpha_ratio*maxval(abs(inputCLS(1:Ncoords,1)))**2
        inputCLS(Ncoords+1,1) = alpha_ratio*sum(inputCLS(1:Ncoords,1)**2)/Natoms

        rmsd_x(1,Nsample) = min_rmsd
        rmsd_x_interpolated(1,Nsample) = min_rmsd
        rmsd_fx_interpolated(1,Nsample) = &
                sqrt(sum((gradient - gradient_var)**2)/Natoms)
        rmsd_fx(1,Nsample) = rmsd_fx_interpolated(1,Nsample)

        rmsd_x2_interpolated(1,Nsample) = inputCLS(Ncoords+1,1)/alpha_ratio

        rmsd_weights(1,1,Nsample) = 1.0d0

        write(filechannel2,FMT="(F15.11,1x,F15.11)") min_rmsd,&
            rmsd_fx(1,Nsample)*(mass_hydrogen / dt) * RU_force / (hartree/bohr)
        exit
    end do

    !Now lets add lots of bad points (but still in threshold)
    do steps = 2, 10000
        do
            do n = 1, Natoms
            do i = 1, 3
                delta_coords(i,n) = rand() - 0.5d0
            end do
            end do
    
            delta_length = sqrt(sum(delta_coords**2))
    
            if (delta_length == 0.0d0) cycle

            coords = coords_final + delta_coords * &
                    (5*(threshold_rmsd)*rand()) / delta_length

            call rmsd_dp(Natoms,coords,coords_final,1,&
                         candidate_U,x_center,y_center,min_rmsd)

            if (min_rmsd >= threshold_rmsd) cycle
            if (min_rmsd <= 1.0d-6) cycle

            call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
            call Acceleration(vals,coords,gradient_var)

            gradient_var = matmul(candidate_U,gradient_var)

            write(filechannel2,FMT="(F15.11,1x,F15.11)") min_rmsd,&
                    sqrt(sum((gradient - gradient_var)**2)/Natoms)*&
                    (mass_hydrogen / dt) * RU_force / (hartree/bohr)

            exit
        end do
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
    
!           if (delta_length >= 1.0d0) cycle
            if (delta_length == 0.0d0) cycle

            coords = coords_final + delta_coords * &
                    (handicap_rmsd + 10*(threshold_rmsd - handicap_rmsd)*rand()) / delta_length

            call rmsd_dp(Natoms,coords,coords_final,1,&
                         candidate_U,x_center,y_center,min_rmsd)

            if (min_rmsd >= threshold_rmsd) cycle

            call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
            call Acceleration(vals,coords,gradient_var)

            gradient_var = matmul(candidate_U,gradient_var)

            if (min_rmsd <= handicap_rmsd) cycle

            write(filechannel2,FMT="(F15.11,1x,F15.11)") min_rmsd,&
                    sqrt(sum((gradient - gradient_var)**2)/Natoms)*&
                    (mass_hydrogen / dt) * RU_force / (hartree/bohr)

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
!       inputCLS(Ncoords+steps,steps) = alpha_ratio * maxval(abs(inputCLS(1:Ncoords,steps)))**2
        inputCLS(Ncoords+steps,steps) = alpha_ratio * sum(inputCLS(1:Ncoords,steps)**2)/Natoms
        
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

        rmsd_x2_interpolated(steps,Nsample) = error2 / alpha_ratio

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
                    (mass_hydrogen / dt) * RU_force / (hartree/bohr)*&
                selected_means(1),&
                    (mass_hydrogen / dt) * RU_force / (hartree/bohr)*&
                selected_SDs(1),&
                    (mass_hydrogen / dt) * RU_force / (hartree/bohr)*&
                selected_means(2),&
                    (mass_hydrogen / dt) * RU_force / (hartree/bohr)*&
                selected_SDs(2),&
                    (mass_hydrogen / dt) * RU_force / (hartree/bohr)*&
                selected_means(3),&
                    (mass_hydrogen / dt) * RU_force / (hartree/bohr)*&
                selected_SDs(3),&
                    (mass_hydrogen / dt) * RU_force / (hartree/bohr)*&
                selected_means(4),&
                    (mass_hydrogen / dt) * RU_force / (hartree/bohr)*&
                selected_SDs(4)

    end do
    close(filechannel2)

    open(filechannel2,file=gridpath5//"linear"//errorcheckfile)
    do steps = 1, Ntest
        do Nsample = 1, Nsamples
        write(filechannel2,FMT="(I3,3(1x,E16.8))") steps, &
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
        write(filechannel2,FMT="(I3,3(1x,E16.8))") steps, &
                rmsd_x_interpolated(steps,min(1+Nsamples/2,69)),&
                rmsd_x2_interpolated(steps,min(1+Nsamples/2,69)),&
                rmsd_fx_interpolated(steps,min(1+Nsamples/2,69))
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
            mean_weights(n) = abs(rmsd_weights(n,steps,Nanomaly_index))
            mean_rmsd_fxs(n) = rmsd_fx(n,Nanomaly_index)
        end do

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

!       do Nsample = 1, Nsamples
!               rmsd_fx_interpolated(:,Nsample) = rmsd_fx_interpolated(:,Nsample) / rmsd_fx(:,Nsample)
!               rmsd_fx(:,Nsample) = rmsd_fx(1,Nsample) / rmsd_fx(:,Nsample)
!       end do

    do Nsample = 1, Nsamples
        rmsd_fx_interpolated(:,Nsample) = rmsd_fx_interpolated(:,Nsample) / rmsd_fx(1,Nsample)
        rmsd_fx(:,Nsample) = rmsd_fx(:,Nsample) / rmsd_fx(1,Nsample)
    end do

    open(filechannel2,file=gridpath5//"relative"//errorcheckfile)
    do steps = 1, Ntest
!       selected_means(1) = sum(rmsd_x_interpolated(steps,:))/Nsamples
!       selected_means(2) = sum(rmsd_x(steps,:))/Nsamples
!       selected_means(3) = sum(rmsd_fx_interpolated(steps,:))/Nsamples
!       selected_means(4) = sum(rmsd_fx(steps,:))/Nsamples

        selected_means(1) = (maxval(rmsd_x_interpolated(steps,:)) +&
                             minval(rmsd_x_interpolated(steps,:)))/2
        selected_means(2) = (maxval(rmsd_x(steps,:)) +&
                             minval(rmsd_x(steps,:)))/2
        selected_means(3) = (maxval(rmsd_fx_interpolated(steps,:)) +&
                             minval(rmsd_fx_interpolated(steps,:)))/2
        selected_means(4) = (maxval(rmsd_fx(steps,:)) +&
                             minval(rmsd_fx(steps,:)))/2

!       selected_SDs(1) = sqrt(sum((rmsd_x_interpolated(steps,:)-&
!                                   selected_means(1))**2)/Nsamples)
!       selected_SDs(2) = sqrt(sum((rmsd_x(steps,:)-&
!                                   selected_means(2))**2)/Nsamples)
!       selected_SDs(3) = sqrt(sum((rmsd_fx_interpolated(steps,:)-&
!                                   selected_means(3))**2)/Nsamples)
!       selected_SDs(4) = sqrt(sum((rmsd_fx(steps,:)-&
!                                   selected_means(4))**2)/Nsamples)

        selected_SDs(1) = (maxval(rmsd_x_interpolated(steps,:))-&
                           minval(rmsd_x_interpolated(steps,:)))/2
        selected_SDs(2) = (maxval(rmsd_x(steps,:))-&
                           minval(rmsd_x(steps,:)))/2
        selected_SDs(3) = (maxval(rmsd_fx_interpolated(steps,:))-&
                           minval(rmsd_fx_interpolated(steps,:)))/2
        selected_SDs(4) = (maxval(rmsd_fx(steps,:))-&
                           minval(rmsd_fx(steps,:)))/2

        write(filechannel2,FMT="(I3,8(1x,F15.11))") steps, &
                selected_means(1), selected_SDs(1),&
                selected_means(2), selected_SDs(2),&
                selected_means(3), selected_SDs(3),&
                selected_means(4), selected_SDs(4)

    end do
    close(filechannel2)

    coords = coords_final
    velocities = velocities_final

    call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

    write(filechannel3,FMT=*) vals,dropoff_mean,dropoff_SD,&
            convergence_mean,convergence_SD

    call plotErrorCheck1(vals,Ntrials,Nsamples,&
            dropoff_mean,dropoff_SD,convergence_mean,convergence_SD)



    !We can end the trajectory after a certain number of
    !trials if the trajectory is going too slow (or is
    !not moving at all!)
    if (Ntrials == 100) exit

    end do

    close(filechannel3)

    deallocate(rmsd_x_interpolated,rmsd_fx,&
               rmsd_fx_interpolated,rmsd_x,&
               rmsd_weights)
    deallocate(inputCLS,gradient_steps)

!   call plotFinalErrorCheck1(Nsamples)

    print *, ""
    print *, "Finished Error Check 2!"
    print *, ""

end subroutine errorCheck2

subroutine errorCheck3(filechannels,coords_final,gradient_final,&
                       Ninterpolation,libcoords,libgradients)
    use FUNCTIONS
    use PARAMETERS
    use PHYSICS
    use VARIABLES
    use ANALYSIS
    use ls_rmsd_original
    implicit none

    !Coordinates, Velocities, and Variables
    real(dp), dimension(3,Natoms),intent(in) :: coords_final,gradient_final
    integer,intent(in) :: Ninterpolation
    real(dp), dimension(Ninterpolation,3,Natoms),intent(in) :: libcoords,libgradients
    real(dp), dimension(3,Natoms) :: coords_labelled,delta_coords
    real(dp), dimension(3,Natoms) :: gradient_var,gradient_labelled,gradient_var_labelled
    real(dp), dimension(3,Natoms) :: approx_gradient,approx_gradient_prime
    real(dp), dimension(Nvar) :: vals
    real(dp), dimension(3,Natoms) :: coords,gradient
    integer :: bond_index1, bond_index2

!   real(dp), allocatable :: rmsd_x(:,:),rmsd_x_interpolated(:,:)
!   real(dp), allocatable :: rmsd_fx(:,:),rmsd_fx_interpolated(:,:)
!   real(dp), allocatable :: rmsd_weights(:,:,:)
!   real(dp),dimension(4) :: selected_means,selected_SDs
    real(dp) :: delta_length
    real(dp),dimension(Ninterpolation,3,Natoms) :: randcoords
    integer :: Ntest,Nsamples,Nsample
    integer :: Nanomaly,Nanomaly_index
    character(3) :: Nanomaly_text
!   real(dp), allocatable :: mean_weights(:),mean_rmsd_fxs(:)
    real(dp), allocatable :: outputCLS(:)
    real(dp), dimension(3,Natoms,Ninterpolation) :: gradient_steps
    real(dp), dimension(Ncoords+Ninterpolation,Ninterpolation) :: inputCLS2
    real(dp), allocatable :: frame_weights(:)
    real(dp), allocatable :: restraints(:,:),restraint_values(:)
    real(dp), allocatable :: minimized_differences2(:,:)
    real(dp) :: error1,error2
!   real(dp), allocatable :: rmsd_x2_interpolated(:,:)
!   real(dp) :: LSn,LSx,LSy,LSx2,LSxy,LSdet
!   real(dp), allocatable :: LSa1(:),LSa2(:),LSerror(:),convergence(:)
!   integer, allocatable :: dropoff(:)
    real(dp) :: dropoff_cutoff,dropoff_mean,dropoff_SD
    real(dp) :: convergence_mean,convergence_SD

    !Various other variables
    real(dp) :: min_rmsd,min_rmsd_prime
    integer :: min_rmsd_index
    integer :: number_of_frames,order,neighbor_check
    character(9) :: vals_interpolation_text

    integer :: Ntrials

    integer,dimension(1+Ngrid_total),intent(in) :: filechannels
    real(dp), dimension(3) :: x_center, y_center
    real(dp), allocatable :: g(:,:)
    real(dp),dimension(3,3) :: U
    real(dp),dimension(3,3) :: candidate_U

    real(dp),dimension(3,Natoms) :: coords_mean
    real(dp) :: mu,sigma

    !Incremental Integer
    integer :: i,n

!   print *, ""
!   print *, "Started Error Check 3!"
!   print *, ""

!   !Initialize the scene
!   call InitialSetup3(coords,velocities)

!   coords_initial = coords
!   velocities_initial = velocities

!   !Always calculate the variables before accelerating
!   call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

!   !Accelerate the velcocities for a half step (verlet)
!   call Acceleration(vals,coords,gradient)

!   !Update the velocities
!   velocities = velocities + 0.5d0 * gradient

!   !To randomize the periods of the bond, I let the scene go on
!   !for a small period of time (need to standardize this later)
!   do n = 1, Nbonds
!       do steps = 1, int(INITIAL_BOND_DATA(6,n)*vib_period)
!           coords = coords + dt * velocities
!           call Acceleration(vals,coords,gradient)
!           velocities = velocities + gradient
!       end do

!       !And then reset the bond
!       coords(:,BONDING_DATA(n,1)) = coords_initial(:,BONDING_DATA(n,1))
!       coords(:,BONDING_DATA(n,2)) = coords_initial(:,BONDING_DATA(n,2))
!       velocities(:,BONDING_DATA(n,1)) = velocities_initial(:,BONDING_DATA(n,1))
!       velocities(:,BONDING_DATA(n,2)) = velocities_initial(:,BONDING_DATA(n,2))
!   end do

    Nanomaly = 0
!   Ntest = 30
    Ntest = Ninterpolation
    Nsamples = 1
    Ntrials = 0
!   allocate(rmsd_x_interpolated(Ntest,Nsamples),rmsd_fx(Ntest,Nsamples),&
!            rmsd_fx_interpolated(Ntest,Nsamples),rmsd_x(Ntest,Nsamples),&
!            rmsd_weights(Ntest,Ntest,Nsamples),rmsd_x2_interpolated(Ntest,Nsamples))
!   allocate(mean_weights(Ntest),mean_rmsd_fxs(Ntest))
!   allocate(gradient_steps(3,Natoms,Ntest))
!   allocate(inputCLS2(Ncoords+Ntest,Ntest))
!   rmsd_weights = 0.0d0

!   open(filechannel3,file=gridpath5//"convergence"//errorcheckfile)




!   !We can choose to start recording later in our
!   !trajectory
!   if (.true.) then
!   do steps = 1, 3000
!       !Upate the coordinates with the velocities
!       coords = coords + dt * velocities

!       !Always calculate the variables before checking a frame or accelerating
!       call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

!       call Acceleration(vals,coords,gradient)

!       !Update the velocities
!       velocities = velocities + gradient
!   end do
!   end if

!   !Start of the "trajectory"
!   do

    Nanomaly_index = 0
!   call system("rm "//gridpath5//"weight1"//errorcheckfile)

!   !We can choose not to move forward in the trajectory
!   !if we are focused on only one point
!   if (.false.) then
!   do steps = 1, 1000
!       !Upate the coordinates with the velocities
!       coords = coords + dt * velocities

!       !Always calculate the variables before checking a frame or accelerating
!       call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

!       call Acceleration(vals,coords,gradient)

!       !Update the velocities
!       velocities = velocities + gradient
!   end do

!   !In this case, we need a counter
!   else
!       Ntrials = Ntrials + 1
!   end if

!   if (any(vals > var_maxvar)) exit

!   coords_final = coords

    call getVarsMaxMin(coords_final,Natoms,vals,Nvar,BOND_LABELLING_DATA)
!   call Acceleration(vals,coords,gradient)

!   do n = 1, Natoms
!       coords_labelled(:,n) = coords(:,BOND_LABELLING_DATA(n))
!       gradient_labelled(:,n) = gradient(:,BOND_LABELLING_DATA(n))
!   end do
!   print *, "vals:", vals

!   open(filechannel2,file=gridpath5//"new"//errorcheckfile)

!   do Nsample = 1, Nsamples

!   handicap_rmsd = threshold_rmsd* 0.90d0

!   !Let's find one good point
!   do
!       do n = 1, Natoms
!       do i = 1, 3
!           delta_coords(i,n) = rand() - 0.5d0
!       end do
!       end do

!       delta_length = sqrt(sum(delta_coords**2))

!       if (delta_length == 0.0d0) cycle

!       coords = coords_final + delta_coords * (handicap_rmsd + &
!               (threshold_rmsd - handicap_rmsd)*rand()) / delta_length

!       call rmsd_dp(Natoms,coords,coords_final,1,candidate_U,x_center,y_center,min_rmsd)

!       if (min_rmsd > handicap_rmsd) cycle
!       
!       handicap_rmsd = min_rmsd

!       call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
!       call Acceleration(vals,coords,gradient_var)

!       gradient_var = matmul(candidate_U,gradient_var)

!       do i = 1, 3
!           coords(i,:) = &
!           coords(i,:) - x_center(i)
!       end do

!       coords = matmul(candidate_U,coords)

!       do i = 1, 3
!           coords(i,:) = &
!           coords(i,:) + y_center(i)
!       end do

!       inputCLS(1:Ncoords,1) = reshape(coords - coords_final,(/Ncoords/))
!       gradient_steps(:,:,1) = gradient_var

!       inputCLS(Ncoords+1,:) = 0.0d0
!       inputCLS(Ncoords+1,1) = alpha_ratio*maxval(abs(inputCLS(1:Ncoords,1)))**2

!       rmsd_x(1,Nsample) = min_rmsd
!       rmsd_x_interpolated(1,Nsample) = min_rmsd
!       rmsd_fx_interpolated(1,Nsample) = &
!               sqrt(sum((gradient - gradient_var)**2)/Natoms)
!       rmsd_fx(1,Nsample) = rmsd_fx_interpolated(1,Nsample)

!       rmsd_x2_interpolated(1,Nsample) = inputCLS(Ncoords+1,1)/alpha_ratio

!       rmsd_weights(1,1,Nsample) = 1.0d0

!       write(filechannel2,FMT="(F15.11,1x,F15.11)") min_rmsd,rmsd_fx(1,Nsample)
!       exit
!   end do

    !Now lets add lots of bad points (but still in threshold)
!   do steps = 2, Ntest
    min_rmsd_prime = 1.0d9
    coords_mean = 0.0d0
    do steps = 1, Ninterpolation
        do
            do n = 1, Natoms
            do i = 1, 3
                delta_coords(i,n) = rand() - 0.5d0
            end do
            end do
    
            delta_length = sqrt(sum(delta_coords**2))
    
            if (delta_length >= 1.0d0) cycle
            if (delta_length == 0.0d0) cycle

!           coords = coords_final + delta_coords * &
!                   (min_rmsd + (threshold_rmsd - min_rmsd)*rand())
            coords = coords_final + delta_coords * &
                    (0.0d0 + (threshold_rmsd - 0.0d0)*rand())

            call rmsd_dp(Natoms,coords,coords_final,1,&
                         candidate_U,x_center,y_center,min_rmsd)

            if (min_rmsd >= threshold_rmsd) cycle

!           call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)
!           call Acceleration(vals,coords,gradient_var)

!           gradient_var = matmul(candidate_U,gradient_var)

!           write(filechannel2,FMT="(F15.11,1x,F15.11)") min_rmsd,&
!                   sqrt(sum((gradient - gradient_var)**2)/Natoms)

!           if (min_rmsd <= handicap_rmsd) cycle
            exit
        end do

        randcoords(steps,:,:) = coords
        coords_mean = coords_mean + &
                coords
    end do

    coords_mean = coords_mean / &
        Ninterpolation
    mu = sqrt(sum((coords_final-coords_mean)**2)/Natoms)
    
    sigma = 0.0d0
    do n = 1, Ninterpolation
        sigma = sigma + &
            sum((coords_mean-randcoords(n,:,:))**2)
    end do
    sigma = sqrt(sigma/Ninterpolation)
    
    open(6666,file=gridpath5//"tmp_A.dat",&
    position="append")
!   write(6666,FMT=*) vals(1), vals(2), mu, sigma, &
!                     coords_labelled, coords_mean
    write(6666,FMT=*) mu, sigma, &
                      coords_final, coords_mean
    close(6666)


    min_rmsd_prime = 1.0d9
    coords_mean = 0.0d0
    do steps = 1, Ninterpolation
        coords(:,:) = libcoords(steps,:,:)
        coords_mean = coords_mean + &
                coords

        call rmsd_dp(Natoms,coords,coords_final,1,&
                     candidate_U,x_center,y_center,min_rmsd)

        if (min_rmsd < min_rmsd_prime) then
            min_rmsd_prime = min_rmsd
            min_rmsd_index = steps
        end if

        gradient_var = matmul(candidate_U,libgradients(steps,:,:))

        do i = 1, 3
            coords(i,:) = &
            coords(i,:) - x_center(i)
        end do

        coords = matmul(candidate_U,coords)

        do i = 1, 3
            coords(i,:) = &
            coords(i,:) + y_center(i)
        end do

        inputCLS2(1:Ncoords,steps) = reshape(coords - coords_final,(/Ncoords/))
        gradient_steps(:,:,steps) = gradient_var

        inputCLS2(Ncoords+steps,:) = 0.0d0
!       inputCLS(Ncoords+Ninterpolation,Ninterpolation) = alpha_ratio * &
!               maxval(abs(inputCLS(1:Ncoords,Ninterpolation)))**2
        inputCLS2(Ncoords+steps,steps) = alpha_ratio * &
                sum(inputCLS2(1:Ncoords,steps)**2)/Natoms

    end do

    allocate(frame_weights(Ninterpolation),&
             outputCLS(Ncoords+Ninterpolation),&
             restraints(1,Ninterpolation),&
             restraint_values(1),&
             minimized_differences2(Ninterpolation,1))
     
    restraints = 1.0d0
    restraint_values = 1.0d0
    outputCLS(1:Ncoords) = 0.0d0
    outputCLS(Ncoords+1:Ncoords+Ninterpolation) = 0.0d0
    
    call CLS2(inputCLS2(1:Ncoords+Ninterpolation,&
              1:Ninterpolation),Ncoords+Ninterpolation,Ninterpolation,&
              restraints,1,restraint_values,&
              outputCLS,frame_weights(1:Ninterpolation))

    Nsample = 1
    approx_gradient = 0.0d0
    coords = 0.0d0
    error2 = 0.0d0
    do n = 1, Ninterpolation
        approx_gradient = approx_gradient + frame_weights(n) *&
                gradient_steps(:,:,n)
        coords = coords + reshape((frame_weights(n)*inputCLS2(1:Ncoords,n)),&
                (/3,Natoms/))
        error2 = error2 + (frame_weights(n)*inputCLS2(Ncoords+n,n))**2
    end do
    error1 = sqrt(sum(coords**2)/Natoms)
    error2 = sqrt(error2)

    gradient_var = gradient_steps(:,:,min_rmsd_index)

!   rmsd_x(Ninterpolation,Nsample) = min_rmsd
!   rmsd_x_interpolated(Ninterpolation,Nsample) = error1
!   rmsd_fx_interpolated(Ninterpolation,Nsample) = &
!           sqrt(sum((gradient - approx_gradient)**2)/Natoms)
!   rmsd_fx(Ninterpolation,Nsample) = sqrt(sum((gradient - gradient_var)**2)/Natoms)

!   rmsd_x2_interpolated(Ninterpolation,Nsample) = error2 / alpha_ratio

!   rmsd_weights(1:Ninterpolation,Ninterpolation,Nsample) = frame_weights(1:Ninterpolation)

!   rmsd_x = min_rmsd
!   rmsd_x_interpolated = error1
!   rmsd_fx_interpolated = &
!           sqrt(sum((gradient-approx_gradient)**2)/Natoms)
!   rmsd_fx = &
!           sqrt(sum((gradient-gradient_var)**2)/Natoms)
    error1 = sqrt(sum((gradient_final-approx_gradient)**2)/Natoms)
    error2 = sqrt(sum((gradient_final-gradient_var)**2)/Natoms)
    if (error1 < error2) then
    print *, "!!!!!!!!!!!!!!!!!!!!!"
    print *, "  Small Error!"
    print *, "!!!!!!!!!!!!!!!!!!!!!"
    print *, "IE:", error1
    print *, "AE:", error2
    print *, "!!!!!!!!!!!!!!!!!!!!!"
    print *, "!!!!!!!!!!!!!!!!!!!!!"
    end if

    coords = coords_final

    call getVarsMaxMin(coords,Natoms,vals,Nvar,BOND_LABELLING_DATA)

!   write(6666,FMT=*) error1, error2

    coords_mean = coords_mean / &
        Ninterpolation
    mu = sqrt(sum((coords_final-coords_mean)**2)/Natoms)
    
    sigma = 0.0d0
    do n = 1, Ninterpolation
        sigma = sigma + &
            sum((coords_mean-libcoords(n,:,:))**2)
    end do
    sigma = sqrt(sigma/Ninterpolation)
    
    open(6666,file=gridpath5//"tmp_B.dat",&
    position="append")
!   write(6666,FMT=*) vals(1), vals(2), mu, sigma, &
!                     coords_labelled, coords_mean
    write(6666,FMT=*) mu, sigma, &
                      coords_final, coords_mean
    close(6666)

    deallocate(frame_weights,outputCLS,restraints,restraint_values,minimized_differences2)

!   if (maxval(rmsd_fx(1:Ninterpolation,Nsample),dim=1) < rmsd_fx_interpolated(Ninterpolation,Nsample)) then
!       Nanomaly_index = Nsample
!       Nanomaly = Nanomaly + 1
!       write(Nanomaly_text,FMT="(I3)") Nanomaly
!       open(filechannel4,file=gridpath5//"anomaly"//&
!               trim(adjustl(Nanomaly_text))//".xyz")
!       write(filechannel4,FMT="(I4)") Natoms*2
!       write(filechannel4,FMT="(A)") "Fixed Frame then Offset Frame"
!       do n = 1, Natoms
!           write(filechannel4,FMT="(A,3(1x,F9.5))") "H", coords_final(:,n)
!       end do
!       do n = 1, Natoms
!           write(filechannel4,FMT="(A,3(1x,F9.5))") "H", coords(:,n)
!       end do
!       close(filechannel4)
!   end if

!   close(filechannel2)

!   open(filechannel2,file=gridpath5//errorcheckfile)
!   do steps = 1, Ninterpolation
!       selected_means(1) = sum(rmsd_x_interpolated(steps,:))/Nsamples
!       selected_means(2) = sum(rmsd_x(steps,:))/Nsamples
!       selected_means(3) = sum(rmsd_fx_interpolated(steps,:))/Nsamples
!       selected_means(4) = sum(rmsd_fx(steps,:))/Nsamples

!       selected_SDs(1) = sqrt(sum((rmsd_x_interpolated(steps,:)-&
!                                   selected_means(1))**2)/Nsamples)
!       selected_SDs(2) = sqrt(sum((rmsd_x(steps,:)-&
!                                   selected_means(2))**2)/Nsamples)
!       selected_SDs(3) = sqrt(sum((rmsd_fx_interpolated(steps,:)-&
!                                   selected_means(3))**2)/Nsamples)
!       selected_SDs(4) = sqrt(sum((rmsd_fx(steps,:)-&
!                                   selected_means(4))**2)/Nsamples)

!       write(filechannel2,FMT="(I3,8(1x,F15.11))") steps, &
!               selected_means(1), selected_SDs(1),&
!               selected_means(2), selected_SDs(2),&
!               selected_means(3), selected_SDs(3),&
!               selected_means(4), selected_SDs(4)

!   end do
!   close(filechannel2)

!   open(filechannel2,file=gridpath5//"linear"//errorcheckfile)
!   do steps = 1, Ninterpolation
!       do Nsample = 1, Nsamples
!       write(filechannel2,FMT="(I3,3(1x,E16.8))") steps, &
!               rmsd_x_interpolated(steps,Nsample),&
!               rmsd_x2_interpolated(steps,Nsample),&
!               rmsd_fx_interpolated(steps,Nsample)
!       end do
!   end do
!   close(filechannel2)

!   allocate(LSa1(Ninterpolation),LSa2(Ninterpolation),LSerror(Ninterpolation))
!   allocate(convergence(Nsample),dropoff(Nsample))

!   open(filechannel2,file=gridpath5//"dropoff"//errorcheckfile)
!   do Nsample = 1, Nsamples

!   convergence(Nsample) = 0.0d0
!   do steps = 1, Ninterpolation/2
!       LSn = Ninterpolation - steps + 1
!       LSx = 0.0d0
!       LSx2 = 0.0d0
!       LSy = 0.0d0
!       LSxy = 0.0d0
!       do i = steps, Ninterpolation
!           LSx  =  LSx + exp(-dble(i))
!           LSy  =  LSy + log10(rmsd_fx_interpolated(i,Nsample))
!           LSxy = LSxy + log10(rmsd_fx_interpolated(i,Nsample))*exp(-dble(i))
!           LSx2 = LSx2 + exp(-dble(i))**2
!       end do
!       LSdet = 1.0d0 / (LSx2*LSn - LSx**2)
!       LSa1(steps) = (LSx2*LSy - LSx*LSxy)*LSdet
!       LSa2(steps) = (LSn*LSxy - LSx*LSy)*LSdet
!       LSerror(steps) = 0.0d0
!       do i = steps, Ninterpolation
!          LSerror(steps) = LSerror(steps) + &
!                  (LSa1(steps) + LSa2(steps)*exp(-dble(i)) - &
!                  log10(rmsd_fx_interpolated(i,Nsample)))**2
!       end do
!       LSerror(steps) = LSerror(steps) / LSn
!       if (LSerror(steps) == 0) then
!           convergence(Nsample) = LSa1(steps)
!           exit
!       end if

!       convergence(Nsample) = convergence(Nsample) + &
!               LSa1(steps) / LSerror(steps)
!   end do

!   if (minval(LSerror(1:Ninterpolation/2)) == 0) then
!   else
!       convergence(Nsample) = convergence(Nsample) /&
!               sum(LSerror(1:Ninterpolation/2)**(-1))
!   end if
!   dropoff_cutoff = sqrt(rmsd_fx_interpolated(1,Nsample)*&
!                           10.0d0**(convergence(Nsample)))
!   dropoff(Nsample) = minloc(rmsd_fx_interpolated(1:Ninterpolation,Nsample),1,&
!           rmsd_fx_interpolated(1:Ninterpolation,Nsample)>dropoff_cutoff) + 1

!   write(filechannel2,FMT="(I3,1x,F15.11,1x,I3)") Nsample, &
!           convergence(Nsample), dropoff(Nsample)
!   end do
!   close(filechannel2)

!   dropoff_mean = sum(dropoff)*1.0d0/Nsamples
!   convergence_mean = sum(convergence)/Nsamples

!   dropoff_SD = sqrt(sum((dropoff-dropoff_mean)**2)/Nsamples)
!   convergence_SD = sqrt(sum((convergence-convergence_mean)**2)/Nsamples)

!   deallocate(LSa1,LSa2,LSerror,convergence,dropoff)

!   open(filechannel2,file=gridpath5//"heatmapline"//errorcheckfile)
!   do steps = 1, Ninterpolation
!       write(filechannel2,FMT="(I3,3(1x,E16.8))") steps, &
!               rmsd_x_interpolated(steps,min(1+Nsamples/2,69)),&
!               rmsd_x2_interpolated(steps,min(1+Nsamples/2,69)),&
!               rmsd_fx_interpolated(steps,min(1+Nsamples/2,69))
!   end do
!   close(filechannel2)

!   open(filechannel2,file=gridpath5//"ratio"//errorcheckfile)
!   do steps = 1, Ninterpolation
!       do Nsample = 1, Nsamples
!       write(filechannel2,FMT="(I3,3(1x,F15.11))") steps, &
!               rmsd_fx_interpolated(steps,Nsample)/&
!               rmsd_fx(1,Nsample),&
!               rmsd_fx_interpolated(steps,Nsample)/&
!               rmsd_fx(steps,Nsample),&
!               rmsd_fx_interpolated(steps,Nsample)/&
!               maxval(rmsd_fx(1:steps,Nsample),dim=1)
!       end do
!   end do
!   close(filechannel2)

!   if (Nanomaly_index > 0) then
!   open(filechannel2,file=gridpath5//"weight1"//errorcheckfile)
!   do steps = 1, Ninterpolation
!       selected_means(1) = sum(rmsd_x_interpolated(steps,:))/Nsamples
!       selected_means(2) = sum(rmsd_x(steps,:))/Nsamples
!       selected_means(3) = sum(rmsd_fx_interpolated(steps,:))/Nsamples
!       selected_means(4) = sum(rmsd_fx(steps,:))/Nsamples

!       do n = 1, Ninterpolation
!           mean_weights(n) = abs(rmsd_weights(n,steps,Nanomaly_index))
!           mean_rmsd_fxs(n) = rmsd_fx(n,Nanomaly_index)
!       end do

!       write(filechannel2,FMT=*) steps,&
!               rmsd_fx_interpolated(steps,Nanomaly_index),&
!               mean_weights, mean_rmsd_fxs

!   end do
!   close(filechannel2)
!   end if

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

!       do Nsample = 1, Nsamples
!               rmsd_fx_interpolated(:,Nsample) = rmsd_fx_interpolated(:,Nsample) / rmsd_fx(:,Nsample)
!               rmsd_fx(:,Nsample) = rmsd_fx(1,Nsample) / rmsd_fx(:,Nsample)
!       end do

!   do Nsample = 1, Nsamples
!       rmsd_fx_interpolated(:,Nsample) = rmsd_fx_interpolated(:,Nsample) / rmsd_fx(1,Nsample)
!       rmsd_fx(:,Nsample) = rmsd_fx(:,Nsample) / rmsd_fx(1,Nsample)
!   end do

!   open(filechannel2,file=gridpath5//"relative"//errorcheckfile)
!   do steps = 1, Ninterpolation
!       selected_means(1) = sum(rmsd_x_interpolated(steps,:))/Nsamples
!       selected_means(2) = sum(rmsd_x(steps,:))/Nsamples
!       selected_means(3) = sum(rmsd_fx_interpolated(steps,:))/Nsamples
!       selected_means(4) = sum(rmsd_fx(steps,:))/Nsamples

!       selected_means(1) = (maxval(rmsd_x_interpolated(steps,:)) +&
!                            minval(rmsd_x_interpolated(steps,:)))/2
!       selected_means(2) = (maxval(rmsd_x(steps,:)) +&
!                            minval(rmsd_x(steps,:)))/2
!       selected_means(3) = (maxval(rmsd_fx_interpolated(steps,:)) +&
!                            minval(rmsd_fx_interpolated(steps,:)))/2
!       selected_means(4) = (maxval(rmsd_fx(steps,:)) +&
!                            minval(rmsd_fx(steps,:)))/2

!       selected_SDs(1) = sqrt(sum((rmsd_x_interpolated(steps,:)-&
!                                   selected_means(1))**2)/Nsamples)
!       selected_SDs(2) = sqrt(sum((rmsd_x(steps,:)-&
!                                   selected_means(2))**2)/Nsamples)
!       selected_SDs(3) = sqrt(sum((rmsd_fx_interpolated(steps,:)-&
!                                   selected_means(3))**2)/Nsamples)
!       selected_SDs(4) = sqrt(sum((rmsd_fx(steps,:)-&
!                                   selected_means(4))**2)/Nsamples)

!       selected_SDs(1) = (maxval(rmsd_x_interpolated(steps,:))-&
!                          minval(rmsd_x_interpolated(steps,:)))/2
!       selected_SDs(2) = (maxval(rmsd_x(steps,:))-&
!                          minval(rmsd_x(steps,:)))/2
!       selected_SDs(3) = (maxval(rmsd_fx_interpolated(steps,:))-&
!                          minval(rmsd_fx_interpolated(steps,:)))/2
!       selected_SDs(4) = (maxval(rmsd_fx(steps,:))-&
!                          minval(rmsd_fx(steps,:)))/2

!       write(filechannel2,FMT="(I3,8(1x,F15.11))") steps, &
!               selected_means(1), selected_SDs(1),&
!               selected_means(2), selected_SDs(2),&
!               selected_means(3), selected_SDs(3),&
!               selected_means(4), selected_SDs(4)

!   end do
!   close(filechannel2)

!   write(filechannel3,FMT=*) vals,dropoff_mean,dropoff_SD,&
!           convergence_mean,convergence_SD

!   call plotErrorCheck1(vals,Ntrials,Nsamples,&
!           dropoff_mean,dropoff_SD,convergence_mean,convergence_SD)



!   !We can end the trajectory after a certain number of
!   !trials if the trajectory is going too slow (or is
!   !not moving at all!)
!   if (Ntrials == 100) exit

!   end do

!   close(filechannel3)

!   deallocate(rmsd_x_interpolated,rmsd_fx,&
!              rmsd_fx_interpolated,rmsd_x,&
!              rmsd_weights)
!   deallocate(inputCLS2,gradient_steps)

!   call plotFinalErrorCheck1(Nsamples)

!   print *, ""
!   print *, "Finished Error Check 3!"
!   print *, ""

end subroutine errorCheck3

subroutine errorCheck4(filechannels,coords1,gradient1,&
                       coords2,gradient2,&
                       Ninterpolation1,Ninterpolation2,&
                       libcoords,libgradients)
    use FUNCTIONS
    use PARAMETERS
    use PHYSICS
    use VARIABLES
    use ANALYSIS
    use ls_rmsd_original
    implicit none

    !Coordinates, Velocities, and Variables
    real(dp), dimension(3,Natoms),intent(in) :: coords1,gradient1
    real(dp), dimension(3,Natoms),intent(in) :: coords2,gradient2
    integer,intent(in) :: Ninterpolation1,Ninterpolation2
    real(dp), dimension(Ninterpolation2,3,Natoms),intent(in) :: libcoords,libgradients
    real(dp), dimension(3,Natoms) :: coords_labelled,delta_coords
    real(dp), dimension(3,Natoms) :: gradient_var,gradient_labelled,gradient_var_labelled
    real(dp), dimension(3,Natoms) :: approx_gradient,approx_gradient_prime
    real(dp), dimension(Nvar) :: vals
    real(dp), dimension(3,Natoms) :: coords,gradient
    integer :: bond_index1, bond_index2

!   real(dp), allocatable :: rmsd_x(:,:),rmsd_x_interpolated(:,:)
!   real(dp), allocatable :: rmsd_fx(:,:),rmsd_fx_interpolated(:,:)
!   real(dp), allocatable :: rmsd_weights(:,:,:)
!   real(dp),dimension(4) :: selected_means,selected_SDs
    real(dp) :: delta_length
    real(dp),dimension(Ninterpolation2,3,Natoms) :: randcoords
    integer :: Ntest,Nsamples,Nsample
    integer :: Nanomaly,Nanomaly_index
    character(3) :: Nanomaly_text
!   real(dp), allocatable :: mean_weights(:),mean_rmsd_fxs(:)
    real(dp), allocatable :: outputCLS(:)
    real(dp), dimension(3,Natoms,Ninterpolation2) :: gradient_steps
    real(dp), dimension(Ncoords+Ninterpolation2,Ninterpolation2) :: inputCLS2
    real(dp), allocatable :: frame_weights(:)
    real(dp), allocatable :: restraints(:,:),restraint_values(:)
    real(dp), allocatable :: minimized_differences2(:,:)
    real(dp) :: error1,error2,error3
!   real(dp), allocatable :: rmsd_x2_interpolated(:,:)
!   real(dp) :: LSn,LSx,LSy,LSx2,LSxy,LSdet
!   real(dp), allocatable :: LSa1(:),LSa2(:),LSerror(:),convergence(:)
!   integer, allocatable :: dropoff(:)
    real(dp) :: dropoff_cutoff,dropoff_mean,dropoff_SD
    real(dp) :: convergence_mean,convergence_SD

    !Various other variables
    real(dp) :: min_rmsd,min_rmsd_prime
    integer :: min_rmsd_index
    integer :: number_of_frames,order,neighbor_check
    character(9) :: vals_interpolation_text

    integer :: Ntrials

    integer,dimension(1+Ngrid_total),intent(in) :: filechannels
    real(dp), dimension(3) :: x_center, y_center
    real(dp), allocatable :: g(:,:)
    real(dp),dimension(3,3) :: U
    real(dp),dimension(3,3) :: candidate_U

    !Incremental Integer
    integer :: i,n

!   print *, ""
!   print *, "Started Error Check 4!"
!   print *, ""

    Nanomaly = 0
    Ntest = Ninterpolation1
    Nsamples = 1
    Ntrials = 0

    Nanomaly_index = 0
!   call system("rm "//gridpath5//"weight1"//errorcheckfile)

    allocate(frame_weights(Ninterpolation2),&
             outputCLS(Ncoords+Ninterpolation2),&
             restraints(1,Ninterpolation2),&
             restraint_values(1),&
             minimized_differences2(Ninterpolation2,1))

    restraints = 1.0d0
    restraint_values = 1.0d0
    outputCLS(1:Ncoords) = 0.0d0
    outputCLS(Ncoords+1:Ncoords+Ninterpolation2) = 0.0d0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      Dropoff Part 1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Now lets add lots of bad points (but still in threshold)
!   do steps = 2, Ntest
    open(6666,file=gridpath5//"tmp_A.dat")
    min_rmsd_prime = 1.0d9
    do steps = 1, Ninterpolation1
        do
            do n = 1, Natoms
            do i = 1, 3
                delta_coords(i,n) = rand() - 0.5d0
            end do
            end do
    
            delta_length = sqrt(sum(delta_coords**2))
    
            if (delta_length >= 1.0d0) cycle
            if (delta_length == 0.0d0) cycle

            coords = coords1 + delta_coords * &
                    (0.0d0 + (threshold_rmsd - 0.0d0)*rand())

            call rmsd_dp(Natoms,coords,coords1,1,&
                         candidate_U,x_center,y_center,min_rmsd)

            if (min_rmsd >= threshold_rmsd) cycle
            randcoords(steps,:,:) = coords
    
            if (min_rmsd < min_rmsd_prime) then
                min_rmsd_prime = min_rmsd
                min_rmsd_index = steps
            end if
    
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
    
            inputCLS2(1:Ncoords,steps) = reshape(coords - coords1,(/Ncoords/))
            gradient_steps(:,:,steps) = gradient_var
    
            inputCLS2(Ncoords+steps,:) = 0.0d0
            inputCLS2(Ncoords+steps,steps) = alpha_ratio * &
                    sum(inputCLS2(1:Ncoords,steps)**2)/Natoms
            exit
        end do
    
        call CLS2(inputCLS2(1:Ncoords+steps,&
                  1:steps),Ncoords+steps,steps,&
                  restraints,1,restraint_values,&
                  outputCLS(1:Ncoords+steps),frame_weights(1:steps))
    
        approx_gradient = 0.0d0
        do n = 1, steps
            approx_gradient = approx_gradient + frame_weights(n) *&
                    gradient_steps(:,:,n)
        end do
    
        gradient_var = gradient_steps(:,:,min_rmsd_index)
    
        !Interpolated Error
        error1 = sqrt(sum((gradient1-approx_gradient)**2)/Natoms)
        !Accept Best Error
        error2 = sqrt(sum((gradient1-gradient_var)**2)/Natoms)
        !Accept Current Error
        error3 = sqrt(sum((gradient1-gradient_steps(:,:,steps))**2)/Natoms)

        write(6666,FMT=*) steps, error1, error2, error3

    end do
    close(6666)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      Dropoff Part 2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Now lets add lots of bad points (but still in threshold)
!   do steps = 2, Ntest
    open(6666,file=gridpath5//"tmp_B.dat")
    min_rmsd_prime = 1.0d9
    do steps = 1, Ninterpolation2
        do
            if (steps > Ninterpolation1) then
                do n = 1, Natoms
                do i = 1, 3
                    delta_coords(i,n) = rand() - 0.5d0
                end do
                end do
        
                delta_length = sqrt(sum(delta_coords**2))
        
                if (delta_length >= 1.0d0) cycle
                if (delta_length == 0.0d0) cycle
    
                coords = coords2 + delta_coords * &
                        (0.0d0 + (threshold_rmsd - 0.0d0)*rand())
    
                call rmsd_dp(Natoms,coords,coords2,1,&
                             candidate_U,x_center,y_center,min_rmsd)
    
                if (min_rmsd >= threshold_rmsd) cycle
            else
                coords = randcoords(steps,:,:)

                call rmsd_dp(Natoms,coords,coords2,1,&
                             candidate_U,x_center,y_center,min_rmsd)
            end if
    
            if (min_rmsd < min_rmsd_prime) then
                min_rmsd_prime = min_rmsd
                min_rmsd_index = steps
            end if
    
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
    
            inputCLS2(1:Ncoords,steps) = reshape(coords - coords2,(/Ncoords/))
            gradient_steps(:,:,steps) = gradient_var
    
            inputCLS2(Ncoords+steps,:) = 0.0d0
            inputCLS2(Ncoords+steps,steps) = alpha_ratio * &
                    sum(inputCLS2(1:Ncoords,steps)**2)/Natoms
            exit
        end do
    
        if (steps >= Ninterpolation1) then
            call CLS2(inputCLS2(1:Ncoords+steps,&
                      1:steps),Ncoords+steps,steps,&
                      restraints,1,restraint_values,&
                      outputCLS(1:Ncoords+steps),frame_weights(1:steps))
        
            approx_gradient = 0.0d0
            do n = 1, steps
                approx_gradient = approx_gradient + frame_weights(n) *&
                        gradient_steps(:,:,n)
            end do
        
            gradient_var = gradient_steps(:,:,min_rmsd_index)
        
            !Interpolated Error
            error1 = sqrt(sum((gradient2-approx_gradient)**2)/Natoms)
            !Accept Best Error
            error2 = sqrt(sum((gradient2-gradient_var)**2)/Natoms)
            !Accept Current Error
            error3 = sqrt(sum((gradient2-gradient_steps(:,:,steps))**2)/Natoms)
    
            write(6666,FMT=*) steps, error1, error2, error3
        end if

    end do
    close(6666)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      First Error
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call getVarsMaxMin(coords1,Natoms,vals,Nvar,BOND_LABELLING_DATA)

    do steps = 1, Ninterpolation1

        do i = 1, Natoms
            coords(:,i) = libcoords(steps,:,BOND_LABELLING_DATA(i))
        end do

        call rmsd_dp(Natoms,coords,coords1,1,&
                     candidate_U,x_center,y_center,min_rmsd)

        gradient_var = matmul(candidate_U,libgradients(steps,:,:))

        do i = 1, 3
            coords(i,:) = &
            coords(i,:) - x_center(i)
        end do

        coords = matmul(candidate_U,coords)

        do i = 1, 3
            coords(i,:) = &
            coords(i,:) + y_center(i)
        end do

        inputCLS2(1:Ncoords,steps) = reshape(coords - coords1,(/Ncoords/))

        do i = 1, Natoms
            gradient_steps(:,i,steps) = gradient_var(:,BOND_LABELLING_DATA(i))
        end do

        inputCLS2(Ncoords+steps,:) = 0.0d0
!       inputCLS(Ncoords+Ninterpolation,Ninterpolation) = alpha_ratio * &
!               maxval(abs(inputCLS(1:Ncoords,Ninterpolation)))**2
        inputCLS2(Ncoords+steps,steps) = alpha_ratio * &
                sum(inputCLS2(1:Ncoords,steps)**2)/Natoms

    end do
     
    restraints = 1.0d0
    restraint_values = 1.0d0
    outputCLS(1:Ncoords) = 0.0d0
    outputCLS(Ncoords+1:Ncoords+Ninterpolation1) = 0.0d0
    
    call CLS2(inputCLS2(1:Ncoords+Ninterpolation1,&
              1:Ninterpolation1),Ncoords+Ninterpolation1,Ninterpolation1,&
              restraints,1,restraint_values,&
              outputCLS(1:Ncoords+Ninterpolation1),&
              frame_weights(1:Ninterpolation1))

    approx_gradient = 0.0d0
    do n = 1, Ninterpolation1
        approx_gradient = approx_gradient + frame_weights(n) *&
                gradient_steps(:,:,n)
    end do

    !First Error
    error1 = sqrt(sum((gradient1-approx_gradient)**2)/Natoms)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      Second Error
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call getVarsMaxMin(coords1,Natoms,vals,Nvar,BOND_LABELLING_DATA)

    do steps = 1, Ninterpolation1

        do i = 1, Natoms
            coords(:,i) = libcoords(steps,:,BOND_LABELLING_DATA(i))
        end do

        call rmsd_dp(Natoms,coords,coords2,1,&
                     candidate_U,x_center,y_center,min_rmsd)

        gradient_var = matmul(candidate_U,libgradients(steps,:,:))

        do i = 1, 3
            coords(i,:) = &
            coords(i,:) - x_center(i)
        end do

        coords = matmul(candidate_U,coords)

        do i = 1, 3
            coords(i,:) = &
            coords(i,:) + y_center(i)
        end do

        inputCLS2(1:Ncoords,steps) = reshape(coords - coords2,(/Ncoords/))

        do i = 1, Natoms
            gradient_steps(:,i,steps) = gradient_var(:,BOND_LABELLING_DATA(i))
        end do

        inputCLS2(Ncoords+steps,:) = 0.0d0
!       inputCLS(Ncoords+Ninterpolation,Ninterpolation) = alpha_ratio * &
!               maxval(abs(inputCLS(1:Ncoords,Ninterpolation)))**2
        inputCLS2(Ncoords+steps,steps) = alpha_ratio * &
                sum(inputCLS2(1:Ncoords,steps)**2)/Natoms

    end do
     
    restraints = 1.0d0
    restraint_values = 1.0d0
    outputCLS(1:Ncoords) = 0.0d0
    outputCLS(Ncoords+1:Ncoords+Ninterpolation1) = 0.0d0
    
    call CLS2(inputCLS2(1:Ncoords+Ninterpolation1,&
              1:Ninterpolation1),Ncoords+Ninterpolation1,Ninterpolation1,&
              restraints,1,restraint_values,&
              outputCLS(1:Ncoords+Ninterpolation1),&
              frame_weights(1:Ninterpolation1))

    approx_gradient = 0.0d0
    do n = 1, Ninterpolation1
        approx_gradient = approx_gradient + frame_weights(n) *&
                gradient_steps(:,:,n)
    end do

    !Second Error
    error2 = sqrt(sum((gradient2-approx_gradient)**2)/Natoms)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      Third Error
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call getVarsMaxMin(coords2,Natoms,vals,Nvar,BOND_LABELLING_DATA)

    do steps = 1, Ninterpolation2

        do i = 1, Natoms
            coords(:,i) = libcoords(steps,:,BOND_LABELLING_DATA(i))
        end do

        call rmsd_dp(Natoms,coords,coords2,1,&
                     candidate_U,x_center,y_center,min_rmsd)

        gradient_var = matmul(candidate_U,libgradients(steps,:,:))

        do i = 1, 3
            coords(i,:) = &
            coords(i,:) - x_center(i)
        end do

        coords = matmul(candidate_U,coords)

        do i = 1, 3
            coords(i,:) = &
            coords(i,:) + y_center(i)
        end do

        inputCLS2(1:Ncoords,steps) = reshape(coords - coords2,(/Ncoords/))

        do i = 1, Natoms
            gradient_steps(:,i,steps) = gradient_var(:,BOND_LABELLING_DATA(i))
        end do

        inputCLS2(Ncoords+steps,:) = 0.0d0
!       inputCLS(Ncoords+Ninterpolation,Ninterpolation) = alpha_ratio * &
!               maxval(abs(inputCLS(1:Ncoords,Ninterpolation)))**2
        inputCLS2(Ncoords+steps,steps) = alpha_ratio * &
                sum(inputCLS2(1:Ncoords,steps)**2)/Natoms

    end do
     
    restraints = 1.0d0
    restraint_values = 1.0d0
    outputCLS(1:Ncoords) = 0.0d0
    outputCLS(Ncoords+1:Ncoords+Ninterpolation2) = 0.0d0
    
    call CLS2(inputCLS2(1:Ncoords+Ninterpolation2,&
              1:Ninterpolation2),Ncoords+Ninterpolation2,Ninterpolation2,&
              restraints,1,restraint_values,&
              outputCLS(1:Ncoords+Ninterpolation2),&
              frame_weights(1:Ninterpolation2))

    approx_gradient = 0.0d0
    do n = 1, Ninterpolation2
        approx_gradient = approx_gradient + frame_weights(n) *&
                gradient_steps(:,:,n)
    end do

    !Third Error
    error3 = sqrt(sum((gradient2-approx_gradient)**2)/Natoms)

    deallocate(frame_weights,outputCLS,restraints,restraint_values,minimized_differences2)








    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      Plotting
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call getVarsMaxMin(coords2,Natoms,vals,Nvar,BOND_LABELLING_DATA)

    open(gnuplotchannel,file=gridpath5//gnuplotfile)
    write(gnuplotchannel,*) "set term pngcairo size 1200,1200"
    write(gnuplotchannel,FMT="(A,2F7.4,A,I0.5,A)") &
            'set output "'//gridpath4, vals,'.png"'
    write(gnuplotchannel,*) 'set title "Error Convergence As More Points Interpolate"'
    write(gnuplotchannel,*) 'set xlabel "Ninterpolation"'
    write(gnuplotchannel,*) 'set ylabel "RMSD Between Interpolated and Target Gradient"'
    write(gnuplotchannel,FMT="(A,I5,A)") 'set label 1 "N = ',1, '" at screen 0.1,0.925'
    write(gnuplotchannel,FMT='(A,F7.4,",",F7.4,A)') 'set label 2 "(var1,var2) = ',vals, '" at screen 0.2,0.925'
    write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 3 "Threshhold = ',threshold_rmsd, ' A" at screen 0.5,0.925'
    write(gnuplotchannel,FMT='(A,E16.8,A)') 'set label 5 "AlphaRatio = ',alpha_ratio, '" at screen 0.50,0.900'
    write(gnuplotchannel,*) 'xmax = ', Ninterpolation2
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
    write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_A.dat" u 1:2 w lines lw 3 lc "green" t "Interpolation",\'
    write(gnuplotchannel,*) '     "'//gridpath5//'tmp_A.dat" u 1:4 w lines lw 3 lc "red" t "Accept Current",\'
    write(gnuplotchannel,*) '     "'//gridpath5//'tmp_A.dat" u 1:3 w lines lw 3 lc "blue" t "Accept Best",\'
    write(gnuplotchannel,*) '     "'//gridpath5//'tmp_B.dat" u 1:2 w lines lw 5 lc "green" t "",\'
    write(gnuplotchannel,*) '     "'//gridpath5//'tmp_B.dat" u 1:4 w lines lw 5 lc "red" t "",\'
    write(gnuplotchannel,*) '     "'//gridpath5//'tmp_B.dat" u 1:3 w lines lw 5 lc "blue" t "",\'
    write(gnuplotchannel,*) '     "<echo ''',Ninterpolation1, ' ',&
                                             error1,'''" w points ls 1 lw 8 lc "black" t "",\'
    write(gnuplotchannel,*) '     "<echo ''',Ninterpolation1, ' ',&
                                             error2,'''" w points ls 2 lw 6 lc "green" t "",\'
    write(gnuplotchannel,*) '     "<echo ''',Ninterpolation2, ' ',&
                                             error3,'''" w points ls 2 lw 6 lc "green" t ""'
    close(gnuplotchannel)

    call system("gnuplot < "//gridpath5//gnuplotfile)



!   print *, ""
!   print *, "Finished Error Check 4!"
!   print *, ""

end subroutine errorCheck4

subroutine errorCheck5(filechannels,coords1,gradient1,&
                       Ninterpolation,&
                       libcoords,libgradients,libNtraj)
    use FUNCTIONS
    use PARAMETERS
    use PHYSICS
    use VARIABLES
    use ANALYSIS
    use ls_rmsd_original
    implicit none

    !Coordinates, Velocities, and Variables
    real(dp), dimension(3,Natoms),intent(in) :: coords1,gradient1
    real(dp), dimension(3,Natoms) :: coords2,gradient2
    integer,intent(in) :: Ninterpolation
    real(dp), dimension(Ninterpolation,3,Natoms),intent(in) :: libcoords,libgradients
    integer,dimension(Ninterpolation) :: libNtraj
    real(dp), dimension(3,Natoms) :: coords_labelled,delta_coords
    real(dp), dimension(3,Natoms) :: gradient_var,gradient_labelled,gradient_var_labelled
    real(dp), dimension(3,Natoms) :: gradient_var_min, gradient_var_max
    real(dp), dimension(3,Natoms) :: approx_gradient,approx_gradient_prime
    real(dp), dimension(Nvar) :: vals
    real(dp), dimension(3,Natoms) :: coords,gradient
    integer :: bond_index1, bond_index2

!   real(dp), allocatable :: rmsd_x(:,:),rmsd_x_interpolated(:,:)
!   real(dp), allocatable :: rmsd_fx(:,:),rmsd_fx_interpolated(:,:)
!   real(dp), allocatable :: rmsd_weights(:,:,:)
!   real(dp),dimension(4) :: selected_means,selected_SDs
    real(dp) :: delta_length
    real(dp),dimension(Ninterpolation,3,Natoms) :: randcoords
    integer :: Ntest,Nsamples,Nsample
    integer :: Nanomaly,Nanomaly_index
    character(3) :: Nanomaly_text
!   real(dp), allocatable :: mean_weights(:),mean_rmsd_fxs(:)
    real(dp), allocatable :: outputCLS(:)
    real(dp), dimension(3,Natoms,Ninterpolation) :: gradient_steps
    real(dp), dimension(Ncoords+Ninterpolation,Ninterpolation) :: inputCLS2
    real(dp), allocatable :: frame_weights(:)
    real(dp), allocatable :: restraints(:,:),restraint_values(:)
    real(dp), allocatable :: minimized_differences2(:,:)
    real(dp) :: error1,error2,error3,error4,rmsd1
!   real(dp), allocatable :: rmsd_x2_interpolated(:,:)
!   real(dp) :: LSn,LSx,LSy,LSx2,LSxy,LSdet
!   real(dp), allocatable :: LSa1(:),LSa2(:),LSerror(:),convergence(:)
!   integer, allocatable :: dropoff(:)
    real(dp) :: dropoff_cutoff,dropoff_mean,dropoff_SD
    real(dp) :: convergence_mean,convergence_SD

    !Various other variables
    real(dp) :: min_rmsd,min_rmsd_prime,max_rmsd_prime
    integer :: min_rmsd_index,max_rmsd_index
    integer :: number_of_frames,order,neighbor_check
    character(9) :: vals_interpolation_text

    integer :: Ntrials

    integer :: Nbins
    integer :: min_rmsd_bin
    real(dp) :: min_rmsd_binwidth
    real(dp),allocatable :: min_rmsd_binning(:)
    real(dp),dimension(Ncoords+Ninterpolation,1) :: cost_final
    real(dp) :: alpha_test,total_cost

    integer :: Nunique, unique_counter
!   real(dp),dimension(3,Natoms) :: meancoords
    real(dp),dimension(Ncoords) :: meancoords
    real(dp),dimension(Ninterpolation) :: mu, sigma
    integer,dimension(Ninterpolation) :: uniqueNtraj
    logical :: unique_flag

    integer,dimension(1+Ngrid_total),intent(in) :: filechannels
    real(dp), dimension(3) :: x_center, y_center
    real(dp), allocatable :: g(:,:)
    real(dp),dimension(3,3) :: U
    real(dp),dimension(3,3) :: candidate_U

    !Incremental Integer
    integer :: i,n
    integer :: step

!   print *, ""
!   print *, "Started Error Check 5!"
!   print *, ""

    Nanomaly = 0
    Ntest = Ninterpolation
    Nsamples = 1
    Ntrials = 0

    Nanomaly_index = 0
!   call system("rm "//gridpath5//"weight1"//errorcheckfile)

    allocate(frame_weights(Ninterpolation),&
             outputCLS(Ncoords+Ninterpolation),&
             restraints(1,Ninterpolation),&
             restraint_values(1),&
             minimized_differences2(Ninterpolation,1))

    restraints = 1.0d0
    restraint_values = 1.0d0
    outputCLS(1:Ncoords) = 0.0d0
    outputCLS(Ncoords+1:Ncoords+Ninterpolation) = 0.0d0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      Dropoff Part 1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call getVarsMaxMin(coords1,Natoms,vals,Nvar,BOND_LABELLING_DATA)

    if (force_NoLabels) then
        do i = 1, Natoms
            coords2(:,i) = coords1(:,i)
            gradient2(:,i) = gradient1(:,i)
        end do
    else
        do i = 1, Natoms
            coords2(:,i) = coords1(:,BOND_LABELLING_DATA(i))
            gradient2(:,i) = gradient1(:,BOND_LABELLING_DATA(i))
        end do
    end if

    !Now lets add lots of bad points (but still in threshold)
!   do steps = 2, Ntest
    open(6666,file=gridpath5//"tmp_A.dat")
    min_rmsd_prime = 1.0d9
    max_rmsd_prime = 0.0d0
    do step = 1, Ninterpolation
!           do n = 1, Natoms
!           do i = 1, 3
!               delta_coords(i,n) = rand() - 0.5d0
!           end do
!           end do
!   
!           delta_length = sqrt(sum(delta_coords**2))
!   
!           if (delta_length >= 1.0d0) cycle
!           if (delta_length == 0.0d0) cycle

!           coords = coords1 + delta_coords * &
!                   (0.0d0 + (threshold_rmsd - 0.0d0)*rand())

        coords = libcoords(step,:,:)
        gradient_var = libgradients(step,:,:)

        call rmsd_dp(Natoms,coords,coords2,1,&
                     candidate_U,x_center,y_center,min_rmsd)

        if (min_rmsd < min_rmsd_prime) then
            min_rmsd_prime = min_rmsd
            min_rmsd_index = step
        end if

        if (min_rmsd > max_rmsd_prime) then
            max_rmsd_prime = min_rmsd
            max_rmsd_index = step
        end if

!       call Acceleration(vals,coords,gradient_var)
        gradient_var = matmul(candidate_U,gradient_var)
        gradient_steps(:,:,step) = gradient_var

        do i = 1, 3
            coords(i,:) = &
            coords(i,:) - x_center(i)
        end do

        coords = matmul(candidate_U,coords)

        do i = 1, 3
            coords(i,:) = &
            coords(i,:) + y_center(i)
        end do

        inputCLS2(1:Ncoords,step) = reshape(coords - coords2,(/Ncoords/))
        inputCLS2(Ncoords+step,:) = 0.0d0
        inputCLS2(Ncoords+step,step) = alpha_ratio * &
                sum(inputCLS2(1:Ncoords,step)**2)/Natoms
    
        call CLS2(inputCLS2(1:Ncoords+step,&
                  1:step),Ncoords+step,step,&
                  restraints,1,restraint_values,&
                  outputCLS(1:Ncoords+step),frame_weights(1:step))
    
        approx_gradient = 0.0d0
        do n = 1, step
            approx_gradient = approx_gradient + frame_weights(n) *&
                    gradient_steps(:,:,n)
        end do
    
        gradient_var_min = gradient_steps(:,:,min_rmsd_index)
        gradient_var_max = gradient_steps(:,:,max_rmsd_index)
    
        !Interpolated Error
        error1 = sqrt(sum((gradient2-approx_gradient)**2)/Natoms)
        !Accept Best Error
        error2 = sqrt(sum((gradient2-gradient_var_min)**2)/Natoms)
        !Accept Current Error
        error3 = sqrt(sum((gradient2-gradient_steps(:,:,step))**2)/Natoms)
        !Accept Worst Error
        error4 = sqrt(sum((gradient2-gradient_var_max)**2)/Natoms)
        !Current RMSD
        rmsd1 = sqrt(sum(inputCLS2(1:Ncoords,step)**2)/Natoms)

        write(6666,FMT=*) step, error1, error2, error3, error4, rmsd1

    end do
    close(6666)


    !Linear binning!
    Nbins = 20
!   min_rmsd_binwidth = (max_rmsd_prime - min_rmsd_prime)/Nbins
    min_rmsd_binwidth = (threshold_rmsd)/Nbins
    allocate(min_rmsd_binning(Nbins))

    min_rmsd_binning = 0
    do step = 1, Ninterpolation
        min_rmsd = sqrt(sum(inputCLS2(1:Ncoords,step)**2)/Natoms)
!       min_rmsd_bin = (min_rmsd - min_rmsd_prime) / &
!               min_rmsd_binwidth
        min_rmsd_bin = floor((min_rmsd - 0.0d0) / &
                min_rmsd_binwidth) + 1
        if (min_rmsd_bin < 1) min_rmsd_bin = 1
        if (min_rmsd_bin > Nbins) min_rmsd_bin = Nbins

        min_rmsd_binning(min_rmsd_bin) = &
                min_rmsd_binning(min_rmsd_bin) + 1

    end do

    open(6666,file=gridpath5//"tmp_B.dat")

!   write(6666,FMT=*) (1.0 - 0.5)*min_rmsd_binwidth + &
!           min_rmsd_prime, 0.0

    do n = 1, Nbins
!       write(6666,FMT=*) (n - 0.5)*min_rmsd_binwidth + &
!               min_rmsd_prime, min_rmsd_binning(n)
        write(6666,FMT=*) (n - 0.5)*min_rmsd_binwidth + &
                0.0d0, min_rmsd_binning(n)
    end do

!   write(6666,FMT=*) (Nbins + 1.0 - 0.5)*min_rmsd_binwidth + &
!           min_rmsd_prime, 0.0

    close(6666)

    deallocate(min_rmsd_binning)

    open(6666,file=gridpath5//"tmp_C.dat")

    do step = 1, Ninterpolation
        min_rmsd = sqrt(sum(inputCLS2(1:Ncoords,step)**2)/Natoms)
!       write(6666,FMT=*) min_rmsd, frame_weights(step), libNtraj(step)
        write(6666,FMT=*) min_rmsd, frame_weights(step)
    end do

    close(6666)

!    Nunique = 0
!    do step = 1, Ninterpolation
!        unique_flag = .true.
!        do n = 1, Nunique
!            if (libNtraj(step) == uniqueNtraj(n)) unique_flag = .false.
!        end do
!
!        if (.not.(unique_flag)) cycle
!
!        Nunique = Nunique + 1
!        uniqueNtraj(Nunique) = libNtraj(step)
!    end do
!
!    do n = 1, Nunique
!        meancoords = 0.0d0
!        unique_counter = 0
!        do step = 1, Ninterpolation
!            if (libNtraj(step) == uniqueNtraj(n)) then
!!               coords = coords + &
!!                       libcoords(step,:,:)
!!               coords = coords + &
!!                       reshape(inputCLS2(1:Ncoords,step),&
!!                       (/3,Natoms/))
!                meancoords = meancoords + &
!                        inputCLS2(1:Ncoords,step)
!                unique_counter = unique_counter + 1
!            end if
!        end do
!
!!       meancoords = reshape(coords,(/Ncoords/)) / unique_counter
!        meancoords = meancoords / unique_counter
!!       mu(n) = sqrt(sum((meancoords-coords1)**2))
!        mu(n) = sqrt(sum((meancoords)**2)/Natoms)
!
!        sigma(n) = 0.0d0
!        do step = 1, Ninterpolation
!            if (libNtraj(step) == uniqueNtraj(n)) then
!!               sigma(n) = sigma(n) + &
!!                       sum((meancoords - libcoords(step,:,:))**2)
!                sigma(n) = sigma(n) + &
!                        sqrt(sum((meancoords - inputCLS2(1:Ncoords,step))**2)/Natoms)
!            end if
!        end do
!!       sigma(n) = sqrt(sigma(n)/unique_counter)
!        sigma(n) = sigma(n)/unique_counter
!    end do






    open(6666,file=gridpath5//"tmp_E.dat")

    write(6666,FMT="(E14.6,1x,I0.2,1x,'""',A,'""')")&
            error2, 1, "accept best"
    write(6666,FMT="(E14.6,1x,I0.2,1x,'""',A,'""')")&
            error3, 2, "accept worst"

    do n = 1, 3

    do step = 1, Ninterpolation
        inputCLS2(Ncoords+step,step) = 1.0d-4 * (10.0d0**(2*n)) * &
                sum(inputCLS2(1:Ncoords,step)**2)/Natoms
    end do
    
    call CLS2(inputCLS2(1:Ncoords+Ninterpolation,&
              1:Ninterpolation),Ncoords+Ninterpolation,Ninterpolation,&
              restraints,1,restraint_values,&
              outputCLS(1:Ncoords+Ninterpolation),&
              frame_weights(1:Ninterpolation))

    if (n == 1) &
            open(6667,file=gridpath5//"tmp_F1.dat")
    if (n == 2) &
            open(6667,file=gridpath5//"tmp_F2.dat")
    if (n == 3) &
            open(6667,file=gridpath5//"tmp_F3.dat")
    approx_gradient = 0.0d0
    do step = 1, Ninterpolation
        approx_gradient = approx_gradient + frame_weights(step) *&
                gradient_steps(:,:,step)

        min_rmsd = sqrt(sum(inputCLS2(1:Ncoords,step)**2)/Natoms)
        write(6667,FMT=*) min_rmsd, frame_weights(step), libNtraj(step)
    end do
    close(6667)

    if (n == 1) &
    write(6666,FMT="(E14.6,1x,I0.2,1x,'""',A,'""')")&
            sqrt(sum((gradient2-approx_gradient)**2)/Natoms),&
            n+2, "alpha = 10^-2"
    if (n == 2) &
    write(6666,FMT="(E14.6,1x,I0.2,1x,'""',A,'""')")&
            sqrt(sum((gradient2-approx_gradient)**2)/Natoms),&
            n+2, "alpha = 10^0"
    if (n == 3) &
    write(6666,FMT="(E14.6,1x,I0.2,1x,'""',A,'""')")&
            sqrt(sum((gradient2-approx_gradient)**2)/Natoms),&
            n+2, "alpha = 10^+2"

    end do
    close(6666)

!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!   print *, "     TEST START"
!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!   cost_final = matmul(inputCLS2(1:Ncoords+Ninterpolation,&
!                                 1:Ninterpolation),&
!                       reshape(frame_weights,(/Ninterpolation,1/)))
!   total_cost = sum(cost_final**2)

!   print *, ""
!   print *, " final weighting"
!   print *, "frame weights:", frame_weights(1), "...", frame_weights(Ninterpolation)
!   print *, "  cost vector:", cost_final(1,1), "...", cost_final(Ninterpolation,1)
!   print *, ""
!   print *, "   total cost:", total_cost

!   do n = 1, 5
!   do step = 1, Ninterpolation
!       frame_weights(step) = rand() - 0.2d0
!   end do
!   if (sum(frame_weights) == 0.0d0) &
!           frame_weights(1) = frame_weights(1) + 1.0d0
!   frame_weights = frame_weights / sum(frame_weights)

!   cost_final = matmul(inputCLS2(1:Ncoords+Ninterpolation,&
!                                 1:Ninterpolation),&
!                       reshape(frame_weights,(/Ninterpolation,1/)))

!   print *, ""
!   print *, "truly random weighting: ", n
!   print *, "frame weights:", frame_weights(1), "...", frame_weights(Ninterpolation)
!   print *, "  cost vector:", cost_final(1,1), "...", cost_final(Ninterpolation,1)
!   print *, ""
!   print *, "   total cost:", sum(cost_final**2)
!   print *, "      LARGER?:",&
!          "                                                                    ", &
!          sum(cost_final**2) > total_cost
!   end do

!   do step = 1, Ninterpolation
!   frame_weights(step) = frame_weights(step) * 1.1d0
!   if (sum(frame_weights) == 0.0d0) &
!           frame_weights(1) = frame_weights(1) + 1.0d0
!   frame_weights = frame_weights / sum(frame_weights)

!   cost_final = matmul(inputCLS2(1:Ncoords+Ninterpolation,&
!                                 1:Ninterpolation),&
!                       reshape(frame_weights,(/Ninterpolation,1/)))

!   print *, ""
!   print *, "scaled random weighting: ", step
!   print *, "frame weights:", frame_weights(1), "...", frame_weights(Ninterpolation)
!   print *, "  cost vector:", cost_final(1,1), "...", cost_final(Ninterpolation,1)
!   print *, ""
!   print *, "   total cost:", sum(cost_final**2)
!   print *, "      LARGER?:",&
!          "                                                                    ", &
!          sum(cost_final**2) > total_cost
!   end do

!   print *, ""
!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!   print *, "     TEST END"
!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!   call sleep(5)
!   print *, ""

!    print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!    print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!    print *, "     TEST START"
!    print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!
!    do n = 1, Natoms
!    do i = 1, 3
!        coords(i,n) = rand()
!    end do
!    end do
!
!    print *, ""
!    print *, "Randomly chosen point from origin picked!"
!    print *, "RMSD: ", sqrt(sum(coords**2))
!    print *, ""
!    print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!
!    do step = 1, 10
!    alpha_test = 10.0d0**(5-step)
!
!    inputCLS2(1:Ncoords,1) = reshape(coords,(/Ncoords/))
!!   inputCLS2(1:Ncoords,2) = reshape(-2.00d0*coords,(/Ncoords/))
!    inputCLS2(1:Ncoords,2) = -2.00d0*reshape(coords,(/Ncoords/))
!
!    inputCLS2(Ncoords+1:Ncoords+2,1:2) = 0.0d0
!!   inputCLS2(Ncoords+1,1) = 1.0d-5
!!   inputCLS2(Ncoords+2,2) = 1.0d-5
!    inputCLS2(Ncoords+1,1) = alpha_test*sqrt(sum(coords**2))
!    inputCLS2(Ncoords+2,2) = alpha_test*sqrt(sum((coords)**2))
!
!!   inputCLS2(3:Ncoords+2,1) = reshape(coords,(/Ncoords/))
!!   inputCLS2(3:Ncoords+2,2) = -2.00d0*reshape(coords,(/Ncoords/))
!
!!   inputCLS2(1:2,1:2) = 0.0d0
!!   inputCLS2(1,1) = alpha_test*sqrt(sum(coords**2))
!!   inputCLS2(2,2) = alpha_test*sqrt(sum((coords)**2))
!
!    call CLS2(inputCLS2(1:Ncoords+2,&
!              1:2),Ncoords+2,2,&
!              restraints(1,1:2),1,&
!              restraint_values,&
!              outputCLS(1:Ncoords+2),&
!              frame_weights(1:2))
!
!    cost_final(1:Ncoords+2,:) = matmul(inputCLS2(1:Ncoords+2,1:2),&
!                        reshape(frame_weights,(/2,1/)))
!    total_cost = sum(cost_final(1:Ncoords+2,1)**2)
!
!    print *, ""
!    write(6,FMT="(A,I0.2)") " alpha test: 10^",5-step
!    print *, "frame weights:", frame_weights(1:2)
!    print *, "  cost vector:", cost_final(1,1)
!    do n = 2, Ncoords+2
!        print *, "              ", cost_final(n,1)
!    end do
!    print *, ""
!    print *, "   total cost:", total_cost
!
!    cost_final(1:Ncoords+2,:) = matmul(inputCLS2(1:Ncoords+2,1:2),&
!                        reshape((/2.0d0/3.0d0, 1.0d0/3.0d0 /),(/2,1/)))
!
!    print *, ""
!    print *, " true weights:", (/2.0d0/3.0d0, 1.0d0/3.0d0/)
!    print *, "  cost vector:", cost_final(1,1)
!    do n = 2, Ncoords+2
!        print *, "              ", cost_final(n,1)
!    end do
!    print *, ""
!    print *, "   total cost:", sum(cost_final(1:Ncoords+2,1)**2)
!    print *, "   delta cost:", sum(cost_final(1:Ncoords+2,1)**2) - total_cost
!
!    end do
!
!    print *, ""
!    print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!    print *, "     TEST END"
!    print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!    print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!
!    call sleep(5)
!    print *, ""

!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!   print *, "     TEST START"
!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!   do n = 1, Natoms
!   do i = 1, 3
!       coords(i,n) = rand()
!       coords2(i,n) = rand()
!   end do
!   end do

!   print *, ""
!   print *, "Randomly chosen point from origin picked!"
!   print *, "RMSD1: ", sqrt(sum(coords**2))
!   print *, "RMSD2: ", sqrt(sum(coords2**2))
!   print *, ""
!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!   do step = 1, 10
!   alpha_test = 10.0d0**(5-step)

!   inputCLS2(1:Ncoords,1) = reshape(coords,(/Ncoords/))
!   inputCLS2(1:Ncoords,2) = -1.0d0*reshape(coords,(/Ncoords/))
!   inputCLS2(1:Ncoords,3) = -1.0d0*reshape(coords2,(/Ncoords/))

!   inputCLS2(Ncoords+1:Ncoords+3,1:3) = 0.0d0
!   inputCLS2(Ncoords+1,1) = alpha_test*sqrt(sum(coords**2))
!   inputCLS2(Ncoords+2,2) = alpha_test*sqrt(sum((coords)**2))
!   inputCLS2(Ncoords+3,3) = alpha_test*sqrt(sum((coords2)**2))

!   call CLS2(inputCLS2(1:Ncoords+3,&
!             1:3),Ncoords+3,3,&
!             restraints(1,1:3),1,&
!             restraint_values,&
!             outputCLS(1:Ncoords+3),&
!             frame_weights(1:3))

!   cost_final(1:Ncoords+3,:) = matmul(inputCLS2(1:Ncoords+3,1:3),&
!                       reshape(frame_weights,(/3,1/)))
!   total_cost = sum(cost_final(1:Ncoords+3,1)**2)

!   print *, ""
!   write(6,FMT="(A,I0.2)") " alpha test: 10^",5-step
!   print *, "frame weights:", frame_weights(1:3)
!   print *, "  cost vector:", cost_final(1,1)
!   do n = 2, Ncoords+3
!       print *, "              ", cost_final(n,1)
!   end do
!   print *, ""
!   print *, "   total cost:", total_cost

!   cost_final(1:Ncoords+3,:) = matmul(inputCLS2(1:Ncoords+3,1:3),&
!                       reshape((/0.5d0,0.5d0,0.0d0 /),(/3,1/)))

!   print *, ""
!   print *, " true weights:", (/0.5d0,0.5d0,0.0d0/)
!   print *, "  cost vector:", cost_final(1,1)
!   do n = 2, Ncoords+3
!       print *, "              ", cost_final(n,1)
!   end do
!   print *, ""
!   print *, "   total cost:", sum(cost_final(1:Ncoords+3,1)**2)
!   print *, "   delta cost:", sum(cost_final(1:Ncoords+3,1)**2) - total_cost

!   end do

!   print *, ""
!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!   print *, "     TEST END"
!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!   call sleep(5)
!   print *, ""






!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!   print *, "     TEST START"
!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!   do step = 1, Ninterpolation
!       min_rmsd = sqrt(sum(inputCLS2(1:Ncoords,step)**2)/Natoms)

!       print *, ""
!       print *, "min_rmsd  (in): ", min_rmsd
!       open(6666,file=gridpath5//"tmp_D.dat",&
!               form="unformatted",position="append")
!       write(6666) min_rmsd
!       close(6666)

!       open(6666,action="read",&
!               form="unformatted",&
!               file=gridpath5//"tmp_D.dat")
!       do n = 1, step
!       read(6666) min_rmsd_prime
!       end do
!       close(6666)

!       print *, "min_rmsd (out): ", min_rmsd_prime
!       print *, "         delta: ", min_rmsd - min_rmsd_prime

!   end do
!   print *, ""
!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!   print *, "     TEST END"
!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!   print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!   call sleep(5)
!   call system("rm "//gridpath5//"tmp_D.dat")
!   print *, ""






    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      Plotting
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call getVarsMaxMin(coords1,Natoms,vals,Nvar,BOND_LABELLING_DATA)

    open(gnuplotchannel,file=gridpath5//gnuplotfile)
    write(gnuplotchannel,*) "set term pngcairo enhanced size 1800,2400"
    write(gnuplotchannel,*) 'set encoding utf8'
    write(gnuplotchannel,FMT="(A,2F7.4,A,I0.5,A)") &
            'set output "'//gridpath4, vals,'.png"'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
write(gnuplotchannel,*) 'unset xtics'
write(gnuplotchannel,*) 'set multiplot layout 5'//&
                        ',1 columnsfirst margins 0.1,0.95,.1,.9 spacing 0.1,0 title '//&
                        '"Single Frame Interpolation'//&
                        '" font ",32" offset 0,-3'
write(gnuplotchannel,*) 'error_scaling = ', RU_energy / eV
write(gnuplotchannel,*) 'unset ylabel'
write(gnuplotchannel,*) 'set y2label "Error (eV/A)" font ",24"'
write(gnuplotchannel,*) 'set ytics font ",16" nomirror'
write(gnuplotchannel,*) 'unset y2tics'
write(gnuplotchannel,*) 'unset xlabel'
!   write(gnuplotchannel,*) 'set title "Error Convergence As More Points Interpolate"'
!   write(gnuplotchannel,*) 'set xlabel "Ninterpolation"'
!   write(gnuplotchannel,*) 'set ylabel "RMSD Between Interpolated and Target Gradient"'
    write(gnuplotchannel,FMT="(A,I5,A)") 'set label 1 "N = ',1, '" at screen 0.1,0.925'
    write(gnuplotchannel,FMT='(A,F7.4,",",F7.4,A)') 'set label 2 "(var1,var2) = ',vals, '" at screen 0.2,0.925'
    write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 3 "Threshhold = ',threshold_rmsd, ' A" at screen 0.5,0.925'
!   write(gnuplotchannel,FMT='(A,E16.8,A)') 'set label 5 "AlphaRatio = ',alpha_ratio, '" at screen 0.50,0.900'
    write(gnuplotchannel,FMT='(A,E16.8,A)') 'set label 5 "AlphaRatio = ',alpha_ratio, '" at screen 0.8,0.925'
    write(gnuplotchannel,*) 'xmax = ', Ninterpolation
    write(gnuplotchannel,*) 'set xrange [1:xmax]'
    write(gnuplotchannel,*) 'set logscale y'
!   write(gnuplotchannel,*) 'set ytics ('//&
!                                '"1e-11" .00000000001, '//&
!                                '"5e-11" .00000000005, '//&
!                                 '"1e-10" .0000000001, '//&
!                                 '"5e-10" .0000000005, '//&
!                                   '"1e-9" .000000001, '//&
!                                   '"5e-9" .000000005, '//&
!                                    '"1e-8" .00000001, '//&
!                                    '"5e-8" .00000005, '//&
!                                     '"1e-7" .0000001, '//&
!                                     '"5e-7" .0000005, '//&
!                                      '"1e-6" .000001, '//&
!                                      '"5e-6" .000005, '//&
!                                      '"1e-5"  .00001, '//&
!                                      '"5e-5"  .00005, '//&
!                                      '"1e-4"   .0001, '//&
!                                      '"5e-4"   .0005, '//&
!                                      '"1e-3"    .001, '//&
!                                      '"5e-3"    .005, '//&
!                                      '"1e-2"     .01, '//&
!                                      '"5e-2"     .05, '//&
!                                      '"1e-1"      .1, '//&
!                                      '"5e-1"      .5, '//&
!                                      ' "1.0"       1, '//&
!                              ')'
    write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_A.dat" u 1:(error_scaling*($2)) '//&
                            'w lines lw 3 lc "green" t "Interpolation",\'
    write(gnuplotchannel,*) '     "'//gridpath5//'tmp_A.dat" u 1:(error_scaling*($4)) '//&
                            'w lines lw 3 lc "black" t "Accept Current",\'
    write(gnuplotchannel,*) '     "'//gridpath5//'tmp_A.dat" u 1:(error_scaling*($3)) '//&
                            'w lines lw 3 lc "blue" t "Accept Best",\'
    write(gnuplotchannel,*) '     "'//gridpath5//'tmp_A.dat" u 1:(error_scaling*($5)) '//&
                            'w lines lw 3 lc "red" t "Accept Worst"'
    write(gnuplotchannel,*) 'set xlabel "Ninterpolation" font ",24"'
write(gnuplotchannel,*) 'unset y2label'
write(gnuplotchannel,*) 'set ylabel "RMSD (A)" font ",24"'
write(gnuplotchannel,*) 'set y2tics font ",16" nomirror'
write(gnuplotchannel,*) 'unset ytics'
    write(gnuplotchannel,*) 'set logscale y2'
write(gnuplotchannel,*) 'set xtics'
    write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_A.dat" u 1:6 axes x1y2 w lines lw 3 lc "black" t ""'
    write(gnuplotchannel,*) 'set multiplot next'
    write(gnuplotchannel,*) 'unset logscale y'
    write(gnuplotchannel,*) 'unset logscale y2'
    write(gnuplotchannel,*) 'set autoscale x'
    write(gnuplotchannel,*) 'set xrange [0:', threshold_rmsd, ']'
write(gnuplotchannel,*) 'unset y2label'
write(gnuplotchannel,*) 'set ylabel "Weighting" font ",24"'
    write(gnuplotchannel,*) 'unset xlabel'
    write(gnuplotchannel,*) 'unset xtics'
    write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set y2tics'
!   write(gnuplotchannel,*) 'ymin = ', minval(frame_weights)
!   write(gnuplotchannel,*) 'ymax = ', maxval(frame_weights)
!   write(gnuplotchannel,*) 'ydelta = (ymax - ymin)*0.05'
!   write(gnuplotchannel,*) 'set yrange [ymin-ydelta:ymax+ydelta]'
!   write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_C.dat" u '//&
!           '1:2:(sprintf(''%d'',$3)) w labels offset 0,1 point pt 6 lw 4 lc "black" t ""'
    write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_C.dat" u '//&
            '1:2 w p pt 6 lw 4 lc "black" t ""'

    if (.false.) then
    write(gnuplotchannel,FMT="(A,I0.2,A,I0.8,A)") &
            'set label ',1,' "Ntraj    RMSD of Mean       SD of Mean" '//&
            'at screen 0.1, character ',76, ' font ",14"'
    do n = 1, Nunique
    write(gnuplotchannel,FMT="(A,I0.2,A,I0.3,6x,E12.7,7x,E12.7,A,I0.8,A)") &
            'set label ',n+1,' "',uniqueNtraj(n),mu(n),sigma(n),'" '//&
            'at screen 0.1, character ',76 - (n)*1, ' font ",14"'
    end do
    end if

write(gnuplotchannel,*) 'unset ylabel'
write(gnuplotchannel,*) 'set y2label "Occurence" font ",24"'
write(gnuplotchannel,*) 'set ytics font ",16" nomirror'
    write(gnuplotchannel,*) 'set xtics font ",16"'
write(gnuplotchannel,*) 'unset y2tics'
    write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set style histogram clustered gap 1'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
    write(gnuplotchannel,*) 'set format x "%6.3f"'
    write(gnuplotchannel,*) 'set xlabel "RMSD (A) Between Current and Target Frame" font ",24"'
    write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_B.dat" u 1:2 w boxes t ""'
    close(gnuplotchannel)

    call system("gnuplot < "//gridpath5//gnuplotfile)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      Second Plotting
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call getVarsMaxMin(coords1,Natoms,vals,Nvar,BOND_LABELLING_DATA)

    open(gnuplotchannel,file=gridpath5//gnuplotfile)
    write(gnuplotchannel,*) "set term pngcairo enhanced size 1800,2400"
    write(gnuplotchannel,*) 'set encoding utf8'
    write(gnuplotchannel,FMT="(A,2F7.4,A,I0.5,A)") &
            'set output "'//gridpath4, vals,'_alphalook.png"'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
!write(gnuplotchannel,*) 'unset xtics'
        write(gnuplotchannel,*) 'set xtics ('//&
                                         '"1e-8" .00000001, '//&
                                         '"2e-8" .00000002, '//&
                                         '"5e-8" .00000005, '//&
                                          '"1e-7" .0000001, '//&
                                         '"2e-7" .00000002, '//&
                                          '"5e-7" .0000005, '//&
                                           '"1e-6" .000001, '//&
                                         '"2e-6" .00000002, '//&
                                           '"5e-6" .000005, '//&
                                           '"1e-5"  .00001, '//&
                                         '"2e-5" .00000002, '//&
                                           '"5e-5"  .00005, '//&
                                           '"1e-4"   .0001, '//&
                                         '"2e-4" .00000002, '//&
                                           '"5e-4"   .0005, '//&
                                           '"1e-3"    .001, '//&
                                         '"2e-3" .00000002, '//&
                                           '"5e-3"    .005, '//&
                                           '"1e-2"     .01, '//&
                                         '"2e-2" .00000002, '//&
                                   ')'
write(gnuplotchannel,*) 'error_scaling = ', RU_energy / eV
write(gnuplotchannel,*) 'set multiplot layout 6'//&
                        ',1 columnsfirst margins 0.1,0.95,.1,.9 spacing 0.1,0 title '//&
                        '"Single Frame Interpolation'//&
                        '" font ",32" offset 0,-3'
write(gnuplotchannel,*) 'unset ylabel'
write(gnuplotchannel,*) 'unset ytics'
    write(gnuplotchannel,FMT="(A,I5,A)") 'set label 1 "N = ',1, '" at screen 0.1,0.925'
    write(gnuplotchannel,FMT='(A,F7.4,",",F7.4,A)') 'set label 2 "(var1,var2) = ',vals, '" at screen 0.2,0.925'
    write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 3 "Threshhold = ',threshold_rmsd, ' A" at screen 0.5,0.925'
!   write(gnuplotchannel,FMT='(A,E16.8,A)') 'set label 5 "AlphaRatio = ',alpha_ratio, '" at screen 0.50,0.900'
    write(gnuplotchannel,FMT='(A,E16.8,A)') 'set label 5 "AlphaRatio = ',alpha_ratio, '" at screen 0.8,0.925'
    write(gnuplotchannel,*) 'xmax = ', Ninterpolation
write(gnuplotchannel,*) 'set xlabel "Error (eV/A)"'
    write(gnuplotchannel,*) 'set logscale x'
    write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_E.dat" u '//&
            '(error_scaling*($1)):2:(sprintf(''%s'',stringcolumn(3))) w labels offset 0,1 point pt 6 lw 4 t ""'

    if (.false.) then
    write(gnuplotchannel,FMT="(A,I0.2,A,I0.8,A)") &
            'set label ',1,' "Ntraj    RMSD of Mean       SD of Mean" '//&
            'at screen 0.1, character ',76, ' font ",14"'
    do n = 1, Nunique
    write(gnuplotchannel,FMT="(A,I0.2,A,I0.3,6x,E12.7,7x,E12.7,A,I0.8,A)") &
            'set label ',n+1,' "',uniqueNtraj(n),mu(n),sigma(n),'" '//&
            'at screen 0.1, character ',76 - (n)*1, ' font ",14"'
    end do
    end if


    write(gnuplotchannel,*) 'set multiplot next'
    write(gnuplotchannel,*) 'unset logscale x'
    write(gnuplotchannel,*) 'set autoscale x'
    write(gnuplotchannel,*) 'set xrange [0:', threshold_rmsd, ']'
write(gnuplotchannel,*) 'set ylabel "Weighting" font ",24"'
    write(gnuplotchannel,*) 'unset xlabel'
    write(gnuplotchannel,*) 'unset xtics'
    write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set y2tics'
!   write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_F1.dat" u '//&
!           '1:2:(sprintf(''%d'',$3)) w labels offset 0,1 point pt 6 lw 4 lc "black" t ""'
!   write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_F2.dat" u '//&
!           '1:2:(sprintf(''%d'',$3)) w labels offset 0,1 point pt 6 lw 4 lc "black" t ""'
!   write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_F3.dat" u '//&
!           '1:2:(sprintf(''%d'',$3)) w labels offset 0,1 point pt 6 lw 4 lc "black" t ""'
    write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_F1.dat" u '//&
            '1:2 w points pt 6 lw 4 lc "black" t ""'
    write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_F2.dat" u '//&
            '1:2 w points pt 6 lw 4 lc "black" t ""'
    write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_F3.dat" u '//&
            '1:2 w points pt 6 lw 4 lc "black" t ""'

write(gnuplotchannel,*) 'unset ylabel'
write(gnuplotchannel,*) 'set y2label "Occurence" font ",24"'
write(gnuplotchannel,*) 'set ytics font ",16" nomirror'
    write(gnuplotchannel,*) 'set xtics font ",16"'
write(gnuplotchannel,*) 'unset y2tics'
    write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set style histogram clustered gap 1'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
    write(gnuplotchannel,*) 'set format x "%6.3f"'
    write(gnuplotchannel,*) 'set xlabel "RMSD (A) Between Current and Target Frame" font ",24"'
    write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_B.dat" u 1:2 w boxes t ""'
    close(gnuplotchannel)

    call system("gnuplot < "//gridpath5//gnuplotfile)




!   print *, ""
!   print *, "Finished Error Check 5!"
!   print *, ""

end subroutine errorCheck5

subroutine errorCheck6(filechannels,coords1,gradient1,&
                       Ninterpolation,&
                       libcoords,libgradients,libNtraj)
    use FUNCTIONS
    use PARAMETERS
    use PHYSICS
    use VARIABLES
    use ANALYSIS
    use ls_rmsd_original
    implicit none

    !Coordinates, Velocities, and Variables
    real(dp), dimension(3,Natoms),intent(in) :: coords1,gradient1
    real(dp), dimension(3,Natoms) :: coords2,gradient2
    integer,intent(in) :: Ninterpolation
    real(dp), dimension(Ninterpolation,3,Natoms),intent(in) :: libcoords,libgradients
    integer,dimension(Ninterpolation) :: libNtraj
    real(dp), dimension(3,Natoms) :: coords_labelled,delta_coords
    real(dp), dimension(3,Natoms) :: gradient_var,gradient_labelled,gradient_var_labelled
    real(dp), dimension(3,Natoms) :: gradient_var_min, gradient_var_max
    real(dp), dimension(3,Natoms) :: approx_gradient,approx_gradient_prime
    real(dp), dimension(Nvar) :: vals
    real(dp), dimension(3,Natoms) :: coords,gradient
    integer :: bond_index1, bond_index2

!   real(dp), allocatable :: rmsd_x(:,:),rmsd_x_interpolated(:,:)
!   real(dp), allocatable :: rmsd_fx(:,:),rmsd_fx_interpolated(:,:)
!   real(dp), allocatable :: rmsd_weights(:,:,:)
!   real(dp),dimension(4) :: selected_means,selected_SDs
    real(dp) :: delta_length
    real(dp),dimension(Ninterpolation,3,Natoms) :: randcoords
    integer :: Ntest,Nsamples,Nsample
    integer :: Nanomaly,Nanomaly_index
    character(3) :: Nanomaly_text
    real(dp), allocatable :: outputCLS(:)
    real(dp), dimension(3,Natoms,Ninterpolation) :: gradient_steps
    real(dp), dimension(Ncoords+Ninterpolation,Ninterpolation) :: inputCLS2
    real(dp), allocatable :: frame_weights(:)
    real(dp), allocatable :: restraints(:,:),restraint_values(:)
    real(dp), allocatable :: minimized_differences2(:,:)
    real(dp) :: error1,error2,error3,error4,rmsd1
    real(dp) :: dropoff_cutoff,dropoff_mean,dropoff_SD
    real(dp) :: convergence_mean,convergence_SD

    !Various other variables
    real(dp) :: min_rmsd,min_rmsd_prime,max_rmsd_prime
    integer :: min_rmsd_index,max_rmsd_index
    integer :: number_of_frames,order,neighbor_check
    character(9) :: vals_interpolation_text

    integer :: Ntrials

    integer :: Nbins
    integer :: min_rmsd_bin
    real(dp) :: min_rmsd_binwidth
    real(dp),allocatable :: min_rmsd_binning(:)
    real(dp),dimension(Ncoords+Ninterpolation,1) :: cost_final
    real(dp) :: alpha_test,total_cost

    integer :: Nunique, unique_counter
    real(dp),dimension(Ncoords) :: meancoords
    real(dp),dimension(Ninterpolation) :: mu, sigma
    integer,dimension(Ninterpolation) :: uniqueNtraj
    logical :: unique_flag

    integer,dimension(1+Ngrid_total),intent(in) :: filechannels
    real(dp), dimension(3) :: x_center, y_center
    real(dp), allocatable :: g(:,:)
    real(dp),dimension(3,3) :: U
    real(dp),dimension(3,3) :: candidate_U

    !Incremental Integer
    integer :: i,n
    integer :: step

!   print *, ""
!   print *, "Started Error Check 6!"
!   print *, ""

    Nanomaly = 0
    Ntest = Ninterpolation
    Nsamples = 1
    Ntrials = 0

    Nanomaly_index = 0

    allocate(frame_weights(Ninterpolation),&
             outputCLS(Ncoords+Ninterpolation),&
             restraints(1,Ninterpolation),&
             restraint_values(1),&
             minimized_differences2(Ninterpolation,1))

    restraints = 1.0d0
    restraint_values = 1.0d0
    outputCLS(1:Ncoords) = 0.0d0
    outputCLS(Ncoords+1:Ncoords+Ninterpolation) = 0.0d0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      Dropoff Part 1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call getVarsMaxMin(coords1,Natoms,vals,Nvar,BOND_LABELLING_DATA)

    if (force_NoLabels) then
        do i = 1, Natoms
            coords2(:,i) = coords1(:,i)
            gradient2(:,i) = gradient1(:,i)
        end do
    else
        do i = 1, Natoms
            coords2(:,i) = coords1(:,BOND_LABELLING_DATA(i))
            gradient2(:,i) = gradient1(:,BOND_LABELLING_DATA(i))
        end do
    end if

    !Now lets add lots of bad points (but still in threshold)
    open(6666,file=gridpath5//"tmp_A.dat")
    min_rmsd_prime = 1.0d9
    max_rmsd_prime = 0.0d0
    do step = 1, Ninterpolation
        coords = libcoords(step,:,:)
        gradient_var = libgradients(step,:,:)

        call rmsd_dp(Natoms,coords,coords2,1,&
                     candidate_U,x_center,y_center,min_rmsd)

        if (min_rmsd < min_rmsd_prime) then
            min_rmsd_prime = min_rmsd
            min_rmsd_index = step
        end if

        if (min_rmsd > max_rmsd_prime) then
            max_rmsd_prime = min_rmsd
            max_rmsd_index = step
        end if

        gradient_var = matmul(candidate_U,gradient_var)
        gradient_steps(:,:,step) = gradient_var

        do i = 1, 3
            coords(i,:) = &
            coords(i,:) - x_center(i)
        end do

        coords = matmul(candidate_U,coords)

        do i = 1, 3
            coords(i,:) = &
            coords(i,:) + y_center(i)
        end do

        inputCLS2(1:Ncoords,step) = reshape(coords - coords2,(/Ncoords/))
        inputCLS2(Ncoords+step,:) = 0.0d0
        inputCLS2(Ncoords+step,step) = alpha_ratio * &
                sum(inputCLS2(1:Ncoords,step)**2)/Natoms

    end do
    close(6666)


    !Linear binning!
    Nbins = 20
    min_rmsd_binwidth = (threshold_rmsd)/Nbins
    allocate(min_rmsd_binning(Nbins))

    min_rmsd_binning = 0
    do step = 1, Ninterpolation
        min_rmsd = sqrt(sum(inputCLS2(1:Ncoords,step)**2)/Natoms)
        min_rmsd_bin = floor((min_rmsd - 0.0d0) / &
                min_rmsd_binwidth) + 1
        if (min_rmsd_bin < 1) min_rmsd_bin = 1
        if (min_rmsd_bin > Nbins) min_rmsd_bin = Nbins

        min_rmsd_binning(min_rmsd_bin) = &
                min_rmsd_binning(min_rmsd_bin) + 1

    end do

    open(6666,file=gridpath5//"tmp_B.dat")

    do n = 1, Nbins
        write(6666,FMT=*) (n - 0.5)*min_rmsd_binwidth + &
                0.0d0, min_rmsd_binning(n)
    end do

    close(6666)

    deallocate(min_rmsd_binning)

    Nunique = 0
    do step = 1, Ninterpolation
        unique_flag = .true.
        do n = 1, Nunique
            if (libNtraj(step) == uniqueNtraj(n)) unique_flag = .false.
        end do

        if (.not.(unique_flag)) cycle

        Nunique = Nunique + 1
        uniqueNtraj(Nunique) = libNtraj(step)
    end do

    do n = 1, Nunique
        meancoords = 0.0d0
        unique_counter = 0
        do step = 1, Ninterpolation
            if (libNtraj(step) == uniqueNtraj(n)) then
                meancoords = meancoords + &
                        inputCLS2(1:Ncoords,step)
                unique_counter = unique_counter + 1
            end if
        end do

        meancoords = meancoords / unique_counter
        mu(n) = sqrt(sum((meancoords)**2)/Natoms)

        sigma(n) = 0.0d0
        do step = 1, Ninterpolation
            if (libNtraj(step) == uniqueNtraj(n)) then
                sigma(n) = sigma(n) + &
                        sqrt(sum((meancoords - inputCLS2(1:Ncoords,step))**2)/Natoms)
            end if
        end do
        sigma(n) = sigma(n)/unique_counter
    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      Plotting
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call getVarsMaxMin(coords1,Natoms,vals,Nvar,BOND_LABELLING_DATA)

    open(gnuplotchannel,file=gridpath5//gnuplotfile)
    write(gnuplotchannel,*) "set term pngcairo enhanced size 1800,2400"
    write(gnuplotchannel,*) 'set encoding utf8'
    write(gnuplotchannel,FMT="(A,2F7.4,A,I0.5,A)") &
            'set output "'//gridpath4, vals,'.png"'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
write(gnuplotchannel,*) 'unset xtics'
write(gnuplotchannel,*) 'set multiplot layout 2'//&
                        ',1 columnsfirst margins 0.1,0.95,.1,.9 spacing 0.1,0 title '//&
                        '"Single Frame Interpolation'//&
                        '" font ",32" offset 0,-3'

    write(gnuplotchannel,FMT="(A,I0.2,A,I0.8,A)") &
            'set label ',1,' "Ntraj    RMSD of Mean       SD of Mean" '//&
            'at screen 0.1, character ',76, ' font ",14"'
    do n = 1, Nunique
    write(gnuplotchannel,FMT="(A,I0.2,A,I0.3,6x,E12.7,7x,E12.7,A,I0.8,A)") &
            'set label ',n+1,' "',uniqueNtraj(n),mu(n),sigma(n),'" '//&
            'at screen 0.1, character ',76 - (n)*1, ' font ",14"'
    end do


write(gnuplotchannel,*) 'unset ylabel'
write(gnuplotchannel,*) 'set y2label "Occurence" font ",24"'
write(gnuplotchannel,*) 'set ytics font ",16" nomirror'
    write(gnuplotchannel,*) 'set xtics font ",16"'
write(gnuplotchannel,*) 'unset y2tics'
    write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set style histogram clustered gap 1'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
    write(gnuplotchannel,*) 'set format x "%6.3f"'
    write(gnuplotchannel,*) 'set xlabel "RMSD(A) Between Current and Target Frame" font ",24"'
    write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_B.dat" u 1:2 w boxes t ""'
    close(gnuplotchannel)

    call system("gnuplot < "//gridpath5//gnuplotfile)


!   print *, ""
!   print *, "Finished Error Check 6!"
!   print *, ""

end subroutine errorCheck6

subroutine errorCheck7(filechannels,coords1,gradient1,&
                       Ninterpolation,&
                       libcoords,libgradients,libNtraj)
    use FUNCTIONS
    use PARAMETERS
    use PHYSICS
    use VARIABLES
    use ANALYSIS
    use ls_rmsd_original
    implicit none

    !Coordinates, Velocities, and Variables
    real(dp), dimension(3,Natoms),intent(in) :: coords1,gradient1
    real(dp), dimension(3,Natoms) :: coords2,gradient2
    integer,intent(in) :: Ninterpolation
    real(dp), dimension(Ninterpolation,3,Natoms),intent(in) :: libcoords,libgradients
    integer,dimension(Ninterpolation) :: libNtraj
    real(dp), dimension(3,Natoms) :: coords_labelled,delta_coords
    real(dp), dimension(3,Natoms) :: gradient_var,gradient_labelled,gradient_var_labelled
    real(dp), dimension(3,Natoms) :: gradient_var_min, gradient_var_max
    real(dp), dimension(3,Natoms) :: approx_gradient,approx_gradient_prime
    real(dp), dimension(Nvar) :: vals
    real(dp), dimension(3,Natoms) :: coords,gradient
    integer :: bond_index1, bond_index2

!   real(dp), allocatable :: rmsd_x(:,:),rmsd_x_interpolated(:,:)
!   real(dp), allocatable :: rmsd_fx(:,:),rmsd_fx_interpolated(:,:)
!   real(dp), allocatable :: rmsd_weights(:,:,:)
!   real(dp),dimension(4) :: selected_means,selected_SDs
    real(dp) :: delta_length
    real(dp),dimension(Ninterpolation,3,Natoms) :: randcoords
    integer :: Ntest,Nsamples,Nsample
    integer :: Nanomaly,Nanomaly_index
    character(3) :: Nanomaly_text
    real(dp), allocatable :: outputCLS(:)
    real(dp), dimension(3,Natoms,Ninterpolation) :: gradient_steps
    real(dp), dimension(Ncoords+Ninterpolation,Ninterpolation) :: inputCLS2
    real(dp), allocatable :: frame_weights(:)
    real(dp), allocatable :: restraints(:,:),restraint_values(:)
    real(dp), allocatable :: minimized_differences2(:,:)
    real(dp) :: error1,error2,error3,error4,rmsd1
    real(dp) :: dropoff_cutoff,dropoff_mean,dropoff_SD
    real(dp) :: convergence_mean,convergence_SD

    !Various other variables
    real(dp) :: min_rmsd,min_rmsd_prime,max_rmsd_prime
    integer :: min_rmsd_index,max_rmsd_index
    integer :: number_of_frames,order,neighbor_check
    character(9) :: vals_interpolation_text

    integer :: Ntrials

    integer :: Nbins
    integer :: min_rmsd_bin
    real(dp) :: min_rmsd_binwidth
    real(dp),allocatable :: min_rmsd_binning(:)
    real(dp),dimension(Ncoords+Ninterpolation,1) :: cost_final
    real(dp) :: alpha_test,total_cost

    integer :: Nunique, unique_counter
    real(dp),dimension(Ncoords) :: meancoords
    real(dp),dimension(Ninterpolation) :: mu, sigma
    integer,dimension(Ninterpolation) :: uniqueNtraj
    logical :: unique_flag

    integer,dimension(1+Ngrid_total),intent(in) :: filechannels
    real(dp), dimension(3) :: x_center, y_center
    real(dp), allocatable :: g(:,:)
    real(dp),dimension(3,3) :: U
    real(dp),dimension(3,3) :: candidate_U

    real(dp) :: min_rmsd1, min_rmsd2, max_rmsd2
    real(dp),dimension(Ninterpolation,3) :: rmsd_catalog
    integer,dimension(3) :: min_rmsd_bins
    integer,allocatable :: rmsd_binning(:,:)

    !Incremental Integer
    integer :: i,n
    integer :: step

!   print *, ""
!   print *, "Started Error Check 7!"
!   print *, ""

    do step = 1, Ninterpolation

        call rmsd_dp(Natoms,libcoords(step,:,:),coords1,1,&
                     candidate_U,x_center,y_center,min_rmsd1)
        
        min_rmsd2 = min_rmsd1
        max_rmsd2 = min_rmsd1

        do n = 1, Nindistinguishables
            BOND_LABELLING_DATA = &
                INDISTINGUISHABLES(n,:)
            do i = 1, Natoms
                coords(:,i) = libcoords(step,:,&
                        BOND_LABELLING_DATA(i))
            end do

            call rmsd_dp(Natoms,coords,coords1,1,&
                         candidate_U,x_center,&
                         y_center,min_rmsd)

            min_rmsd2 = min(min_rmsd2,min_rmsd)
            max_rmsd2 = max(max_rmsd2,min_rmsd)
        end do

        rmsd_catalog(step,:) = &
                (/min_rmsd1, min_rmsd2, max_rmsd2/)
    end do

    Nbins = 20
    min_rmsd2 = minval(rmsd_catalog(:,2))
    max_rmsd2 = maxval(rmsd_catalog(:,1))
    min_rmsd_binwidth = (max_rmsd2 - min_rmsd2) / Nbins

    allocate(rmsd_binning(3,Nbins))
    rmsd_binning = 0

    do step = 1, Ninterpolation

        min_rmsd_bins = floor((rmsd_catalog(step,:) - &
                min_rmsd2) / min_rmsd_binwidth)
        do i = 1, 2
            if (min_rmsd_bins(i) > Nbins) &
                min_rmsd_bins(i) = Nbins
            if (min_rmsd_bins(i) < 1) &
                min_rmsd_bins(i) = 1

            rmsd_binning(i,min_rmsd_bins(i)) =&
            rmsd_binning(i,min_rmsd_bins(i)) + 1
        end do
    end do

    open(6666,file=gridpath5//"tmp_A.dat")
    do n = 1, Nbins
        write(6666,FMT=*) min_rmsd2 + (n-0.5)*min_rmsd_binwidth,&
                rmsd_binning(1:2,n)
    end do
    close(6666)

    deallocate(rmsd_binning)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      Plotting
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call getVarsMaxMin(coords1,Natoms,vals,Nvar,BOND_LABELLING_DATA)

    open(gnuplotchannel,file=gridpath5//gnuplotfile)
    write(gnuplotchannel,*) "set term pngcairo enhanced size 1800,2400"
    write(gnuplotchannel,*) 'set encoding utf8'
    write(gnuplotchannel,FMT="(A,2F7.4,A,I0.5,A)") &
            'set output "'//gridpath4, vals,'.png"'
write(gnuplotchannel,*) 'set tmargin 0'
write(gnuplotchannel,*) 'set bmargin 0'
write(gnuplotchannel,*) 'set lmargin 1'
write(gnuplotchannel,*) 'set rmargin 1'
write(gnuplotchannel,*) 'unset xtics'
write(gnuplotchannel,*) 'set multiplot layout 2'//&
                        ',1 columnsfirst margins 0.1,0.95,.1,.9 spacing 0.1,0 title '//&
                        '"Single Frame Interpolation'//&
                        '" font ",32" offset 0,-3'

write(gnuplotchannel,*) 'set ylabel "Occurence" font ",24"'
write(gnuplotchannel,*) 'set ytics font ",16" nomirror'
    write(gnuplotchannel,*) 'set xtics font ",16"'
    write(gnuplotchannel,*) 'set autoscale y'
write(gnuplotchannel,*) 'set yrange [0:]'
write(gnuplotchannel,*) 'set style histogram clustered gap 1'
write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
    write(gnuplotchannel,*) 'set format x ""'
    write(gnuplotchannel,*) 'unset xlabel'
write(gnuplotchannel,FMT="(A,F9.7,':',F9.7,A)") 'set xrange [',&
                                        min_rmsd2 - min_rmsd_binwidth,&
                                        max_rmsd2 + min_rmsd_binwidth, ']'
    write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_A.dat" u 1:2 w boxes t ""'
!   write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_A.dat" u 1:3 w boxes t ""'
    write(gnuplotchannel,*) 'set format x "%6.3f"'
    write(gnuplotchannel,*) 'set xlabel "RMSD(A) Between Current and Target Frame" font ",24"'
    write(gnuplotchannel,*) 'plot "'//gridpath5//'tmp_A.dat" u 1:3 w boxes t ""'
    close(gnuplotchannel)

    call system("gnuplot < "//gridpath5//gnuplotfile)


!   print *, ""
!   print *, "Finished Error Check 7!"
!   print *, ""

end subroutine errorCheck7

subroutine errorCheck8(filechannels,coords1,gradient1,&
                       Ninterpolation,&
                       libcoords,libgradients,&
                       alpha_flagging)
    use FUNCTIONS
    use PARAMETERS
    use PHYSICS
    use VARIABLES
    use ANALYSIS
    use ls_rmsd_original
    implicit none

    !Coordinates, Velocities, and Variables
    real(dp), dimension(3,Natoms),intent(in) :: coords1,gradient1
    integer,intent(in) :: Ninterpolation
    real(dp), dimension(Ninterpolation,3,Natoms),intent(in) :: libcoords,libgradients
    integer,dimension(Nalpha),intent(out) :: alpha_flagging
    real(dp), dimension(Nvar) :: vals
    real(dp), dimension(3,Natoms) :: coords,gradient,approx_gradient

    real(dp), dimension(Ncoords+Ninterpolation) :: outputCLS
    real(dp), dimension(3,Natoms,Ninterpolation) :: gradient_steps
    real(dp), dimension(Ncoords+Ninterpolation,Ninterpolation) :: inputCLS2
    real(dp), dimension(Ninterpolation) :: frame_weights
    real(dp), dimension(1,Ninterpolation) :: restraints
    real(dp), dimension(1) :: restraint_values

    !Various other variables
    real(dp) :: error1, max_error1, min_error1
    real(dp) :: best_error
    real(dp),dimension(Ninterpolation) :: min_rmsd

    real(dp),dimension(Nalpha) :: alpha_array

    integer,dimension(1+Ngrid_total),intent(in) :: filechannels
    real(dp), dimension(3) :: x_center, y_center
    real(dp), allocatable :: g(:,:)
    real(dp),dimension(3,3) :: U
    real(dp),dimension(3,3) :: candidate_U

    character(6) :: steps_text

    !Incremental Integer
    integer :: i,n,iostate
    integer :: step

    write(steps_text,FMT="(I0.6)") steps

!   print *, ""
!   print *, "Started Error Check 8!"
!   print *, ""

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      First Graph
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do step = 1, Ninterpolation

        coords = libcoords(step,:,:)

        call rmsd_dp(Natoms,coords,coords1,1,&
                     candidate_U,x_center,y_center,min_rmsd(step))

        gradient_steps(:,:,step) = matmul(candidate_U,libgradients(step,:,:))

        do i = 1, 3
            coords(i,:) = &
            coords(i,:) - x_center(i)
        end do

        coords = matmul(candidate_U,coords)

        do i = 1, 3
            coords(i,:) = &
            coords(i,:) + y_center(i)
        end do

        inputCLS2(1:Ncoords,step) = reshape(coords - coords1,(/Ncoords/))

        inputCLS2(Ncoords+step,:) = 0.0d0

    end do

    do n = 1, Nalpha-1
        alpha_array(n) = alpha_start + &
            (n-1)*(alpha_end-alpha_start)/&
            (1.0d0*(Nalpha-1))
    end do
    alpha_array(Nalpha) = alpha_end

    if (logarithmic_alpha_flag) then
        do n = 1, Nalpha
            alpha_array(n) = &
                10.0d0 ** alpha_array(n)
        end do
    end if

    !The frame with the lowest RMSD should be the first
    !in the library because of the sorting
    best_error = sqrt(sum((gradient1-gradient_steps(:,:,1))**2)/Natoms)
    alpha_flagging = 0

    max_error1 = 0.0d0
    min_error1 = 1.0d9

    open(6667,file=gridpath5//steps_text//"A.dat")
    do n = 1, Nalpha

        do step = 1, Ninterpolation
            inputCLS2(Ncoords+step,step) = alpha_array(n) * &
                    sum(inputCLS2(1:Ncoords,step)**2)/Natoms
        end do
         
        restraints = 1.0d0
        restraint_values = 1.0d0
        outputCLS(1:Ncoords) = 0.0d0
        outputCLS(Ncoords+1:Ncoords+Ninterpolation) = 0.0d0
        
        call CLS2(inputCLS2(1:Ncoords+Ninterpolation,&
                  1:Ninterpolation),Ncoords+Ninterpolation,Ninterpolation,&
                  restraints,1,restraint_values,&
                  outputCLS(1:Ncoords+Ninterpolation),&
                  frame_weights(1:Ninterpolation))

        iostate = 0
        if (all(frame_weights == 0.0d0)) then
            frame_weights(1) = 1.0d0
            iostate = 1
        end if
    
        approx_gradient = 0.0d0
        do step = 1, Ninterpolation
            approx_gradient = approx_gradient + frame_weights(step) *&
                    gradient_steps(:,:,step)
        end do
    
        error1 = sqrt(sum((gradient1-approx_gradient)**2)/Natoms)
        min_error1 = min(min_error1,error1)
        max_error1 = max(max_error1,error1)

        write(6667,FMT=*) alpha_array(n), &
                 error1, iostate
!                error1 * RU_energy / eV, iostate

        if (error1 > best_error) alpha_flagging(n) =&
                alpha_flagging(n) + 1

    end do
    close(6667)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      Second Graph
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! reset the bounds for readability
    max_error1 = 0.0d0
    min_error1 = 1.0d9

    open(6667,file=gridpath5//steps_text//"B.dat")
    do step = 1, Ninterpolation

        error1 = sqrt(sum((gradient1-gradient_steps(:,:,step))**2)/Natoms)
        min_error1 = min(min_error1,error1)
        max_error1 = max(max_error1,error1)

!       write(6667,FMT=*) n, min_rmsd(step), error1 * RU_energy / eV
        write(6667,FMT=*) n, min_rmsd(step), error1

    end do
    close(6667)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      Plotting
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   min_error1 = min_error1 * RU_energy / eV
!   max_error1 = max_error1 * RU_energy / eV

    call getVarsMaxMin(coords1,Natoms,vals,Nvar,BOND_LABELLING_DATA)

    open(gnuplotchannel,file=gridpath5//gnuplotfile)
    write(gnuplotchannel,*) "set term pngcairo enhanced size 2400,1800"
    write(gnuplotchannel,*) 'set encoding utf8'
    write(gnuplotchannel,FMT="(A,2F7.4,A)") &
            'set output "'//gridpath4//steps_text//"_", vals,'.png"'
    write(gnuplotchannel,*) 'set tmargin 0'
    write(gnuplotchannel,*) 'set bmargin 0'
    write(gnuplotchannel,*) 'set lmargin 1'
    write(gnuplotchannel,*) 'set rmargin 1'
    write(gnuplotchannel,*) 'set multiplot layout 1'//&
                            ',2 rowsfirst margins 0.1,0.95,.1,.9 spacing 0,0.1 title '//&
                            '"Single Frame Interpolation, Variable Alpha'//&
                            '" font ",32" offset 0,-3'

    write(gnuplotchannel,*) 'set ylabel "Error (eV/A)" font ",24"'
    write(gnuplotchannel,*) 'set ytics font ",16" nomirror'
    write(gnuplotchannel,*) 'set xtics font ",16"'
    write(gnuplotchannel,*) 'set autoscale y'
!   write(gnuplotchannel,FMT="(A,E16.6,':',E16.6,A)") &
!           'set yrange [',&
!           (min_error1 - (max_error1-min_error1)/Nalpha),&
!           (max_error1 + (max_error1-min_error1)/Nalpha),']'
    ! reset the bounds for readability
    write(gnuplotchannel,*) 'maxE = ', (max_error1 + (max_error1-min_error1)/Nalpha)
    write(gnuplotchannel,FMT="(A,E16.6,':',E16.6,A)") &
            'set yrange [0:maxE]'
    if (logarithmic_alpha_flag) then
        write(gnuplotchannel,*) 'set logscale x'
        write(gnuplotchannel,FMT="(A,E16.6,':',E16.6,A)") &
                'set xrange [',&
                10.0d0 ** (alpha_start - (alpha_end-alpha_start)/Nalpha),&
                10.0d0 ** (alpha_end + (alpha_end-alpha_start)/Nalpha),']'
    else
        write(gnuplotchannel,FMT="(A,E16.6,':',E16.6,A)") &
                'set xrange [',&
                (alpha_start - (alpha_end-alpha_start)/Nalpha),&
                (alpha_end + (alpha_end-alpha_start)/Nalpha),']'
    end if
    write(gnuplotchannel,*) 'set xlabel "Alpha Ratio" font ",24"'
    write(gnuplotchannel,*) 'plot "'//gridpath5//steps_text//'A.dat" u 1:($3==0?$2:1/0) w points ps 3 pt 7 lc "black" t "",\'
    write(gnuplotchannel,*) '     "'//gridpath5//steps_text//'A.dat" u 1:($2>maxE?maxE:1/0) w points ps 3 pt 7 lc "purple" t "",\'
    write(gnuplotchannel,*) '     "'//gridpath5//steps_text//'A.dat" u 1:($3==1?$2:1/0) w points ps 3 pt 7 lc "red" t ""'

    write(gnuplotchannel,*) 'unset ytics'
!   write(gnuplotchannel,FMT="(A,E16.6,':',E16.6,A)") &
!           'set yrange [',&
!           (min_error1 - (max_error1-min_error1)/Nalpha),&
!           (max_error1 + (max_error1-min_error1)/Nalpha),']'
    ! reset the bounds for readability
    write(gnuplotchannel,FMT="(A,E16.6,':',E16.6,A)") &
            'set yrange [0:maxE]'
    write(gnuplotchannel,*) 'set format y ""'
    write(gnuplotchannel,*) 'unset ylabel'
    if (logarithmic_alpha_flag) write(gnuplotchannel,*) 'unset logscale x'
    write(gnuplotchannel,*) 'set autoscale x'
    write(gnuplotchannel,FMT="(A,E16.6,':',E16.6,A)") &
            'set xrange [',&
            (minval(min_rmsd) - (maxval(min_rmsd)-minval(min_rmsd))/Nalpha),&
            (maxval(min_rmsd) + (maxval(min_rmsd)-minval(min_rmsd))/Nalpha),']'
    write(gnuplotchannel,*) 'set xlabel "RMSD(A) Between Current and Target Frame" font ",24"'
    write(gnuplotchannel,*) 'plot "'//gridpath5//steps_text//'B.dat" u 2:3 w points ps 3 pt 7 t ""'
    close(gnuplotchannel)

    call system("gnuplot < "//gridpath5//gnuplotfile)


!   print *, ""
!   print *, "Finished Error Check 8!"
!   print *, ""

end subroutine errorCheck8





subroutine plotErrorCheck1(vals,Ntrials,Nsamples,&
        dropoff_mean,dropoff_SD,convergence_mean,convergence_SD)
        use FUNCTIONS
        use PARAMETERS
        use PHYSICS
        use VARIABLES
        use ANALYSIS
        implicit none

        real(dp), dimension(Nvar), intent(in) :: vals
        integer, intent(in) :: Ntrials
        integer, intent(in) :: Nsamples
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
        write(gnuplotchannel,FMT="(A,2F7.4,A,I0.5,A)") &
                'set output "'//gridpath4, vals, '_', Ntrials, '.png"'
        write(gnuplotchannel,*) 'set title "Error Convergence As More Points Interpolate"'
        write(gnuplotchannel,*) 'set xlabel "Ninterpolation"'
        write(gnuplotchannel,*) 'set ylabel "RMSD Between Interpolated and Target Gradient"'
        write(gnuplotchannel,FMT="(A,I5,A)") 'set label 1 "N = ', Nsamples, '" at screen 0.1,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,",",F7.4,A)') 'set label 2 "(var1,var2) = ',vals, '" at screen 0.2,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 3 "Threshhold = ',threshold_rmsd, ' A" at screen 0.5,0.925'
        write(gnuplotchannel,FMT='(A,F7.2,A,F7.4,A)') &
                'set label 4 "log(Convergence) = ',convergence_mean,&
                             '     Dropoff = ',dropoff_mean,'" at screen 0.1,0.900'
        write(gnuplotchannel,FMT='(A,E16.8,A)') 'set label 5 "AlphaRatio = ',alpha_ratio, '" at screen 0.50,0.900'
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

!       max_rmsd = 0.0d0
!       min_rmsd = 1.0d9
        max_rmsd = 1.0d0
        min_rmsd = 1.0d0
        Ninterpolation = 0

        open(filechannel2,file=gridpath5//"relative"//errorcheckfile)
        do
                read(filechannel2,iostat=iostate,FMT=*) n, rmsd_x,dum1,&
                        rmsd_x_interpolated,dum2,&
                        rmsd_fx_interpolated,dum3,rmsd_fx,dum4
                if (iostate /= 0) exit

!               max_rmsd = max(max_rmsd,rmsd_fx_interpolated+dum3,rmsd_fx+dum4)
!               min_rmsd = min(min_rmsd,rmsd_fx_interpolated,rmsd_fx)
                max_rmsd = max(max_rmsd,rmsd_fx_interpolated+dum3)
                min_rmsd = min(min_rmsd,rmsd_fx_interpolated-dum3)
                Ninterpolation = n
        end do
        close(filechannel2)

        open(gnuplotchannel,file=gridpath5//gnuplotfile)
        write(gnuplotchannel,*) "set term pngcairo size 1200,1200"
        write(gnuplotchannel,FMT="(A,2F7.4,A,I0.5,A)") &
                'set output "'//gridpath4//"relative_", vals, '_', Ntrials, '.png"'
        write(gnuplotchannel,*) 'set title "Error Convergence As More Points Interpolate"'
        write(gnuplotchannel,*) 'set xlabel "Ninterpolation"'
        write(gnuplotchannel,*) 'set ylabel "Relative RMSD Between Interpolated and Target Gradient"'
        write(gnuplotchannel,FMT="(A,I5,A)") 'set label 1 "N = ', Nsamples, '" at screen 0.1,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,",",F7.4,A)') 'set label 2 "(var1,var2) = ',vals, '" at screen 0.2,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 3 "Threshhold = ',threshold_rmsd, ' A" at screen 0.5,0.925'
        write(gnuplotchannel,FMT='(A,E16.8,A)') 'set label 4 "AlphaRatio = ',alpha_ratio, '" at screen 0.50,0.895'
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
!               '" u 1:(1.0) w lines lc "red" t "Accept Current",\'
                '" u 1:(1.0) w lines lc "blue" t "Accept Best",\'
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
        write(gnuplotchannel,FMT="(A,2F7.4,A,I0.5,A)") &
                'set output "'//gridpath4//"weight1_", vals, '_', Ntrials, '.png"'
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
        write(gnuplotchannel,FMT='(A,E16.8,A)') 'set label 3 "AlphaRatio = ',alpha_ratio, '" at screen 0.50,0.900'
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
        write(gnuplotchannel,FMT="(A,2F7.4,A,I0.5,A)") &
                'set output "'//gridpath4//"new", vals, '_', Ntrials, '.png"'
        write(gnuplotchannel,*) 'set title "Error Convergence As Frames Get Closer"'
        write(gnuplotchannel,*) 'unset key'
        write(gnuplotchannel,*) 'set xlabel "RMSD Between Nearby and Target Frame"'
        write(gnuplotchannel,*) 'set ylabel "Error Between Nearby and Target Gradient"'
        write(gnuplotchannel,FMT="(A,I5,A)") 'set label 1 "N = ', Nsamples, '" at screen 0.1,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,",",F7.4,A)') 'set label 2 "(var1,var2) = ',vals, '" at screen 0.2,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 3 "Threshhold = ',threshold_rmsd, ' A" at screen 0.5,0.925'
        write(gnuplotchannel,*) 'set xrange [:', threshold_rmsd, ']'
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
        write(gnuplotchannel,FMT="(A,2F7.4,A,I0.5,A)") &
                'set output "'//gridpath4//"linear", vals, '_', Ntrials, '.png"'
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







        !Make a heatmap?
        if (.false.) then

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
                'set output "'//gridpath4//"heatmap_linear", vals, '.png"'
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
        end if







        !Make a heatmap gif?
        if (.false.) then
        print *, "Started making the gif!"

        call system("rm -r "//gridpath4//"png/")
        call system("mkdir "//gridpath4//"png")
        Ngif = Ninterpolation
        do n = 1, Ngif
        write(ntext,FMT="(I0.3)") n

        open(gnuplotchannel,file=gridpath5//gnuplotfile)
        write(gnuplotchannel,*) 'set term pngcairo size 2400,1200'
        write(gnuplotchannel,FMT="(A)") &
                'set output "'//gridpath4//"png/"//ntext//'.png"'
        write(gnuplotchannel,*) 'set multiplot'
        write(gnuplotchannel,*) 'set size 0.5, 1.0'
        write(gnuplotchannel,*) 'set origin 0.0, 0.0'
        write(gnuplotchannel,*) 'set title "Interpolation Error Comparison by RMSD Variables 1 and 2"'
        write(gnuplotchannel,*) 'set pm3d map'
        write(gnuplotchannel,*) 'unset key'
        write(gnuplotchannel,FMT="(A,I5,A)") 'set label 1 "N = ', Nsamples, '" at screen 0.055,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,",",F7.4,A)') 'set label 2 "(var1,var2) = ',vals, '" at screen 0.1,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 3 "Threshhold = ',threshold_rmsd, ' A" at screen 0.25,0.925'
        write(gnuplotchannel,FMT='(A,E16.8,A)') 'set label 4 "AlphaRatio = ',alpha_ratio, '" at screen 0.25,0.900'
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

        write(gnuplotchannel,*) "set lmargin at screen 0.05"
        write(gnuplotchannel,*) "set rmargin at screen 0.45"
        write(gnuplotchannel,*) "set bmargin at screen 0.10"
        write(gnuplotchannel,*) "set tmargin at screen 0.95"
        write(gnuplotchannel,*) 'splot "'//gridpath5//'heatmap_linear'//errorcheckfile//&
                '" u 1:2:3 w image palette'

        write(gnuplotchannel,*) "unset title"
        write(gnuplotchannel,*) "unset xlabel"
        write(gnuplotchannel,*) "unset ylabel"
        write(gnuplotchannel,*) 'unset xtics'
        write(gnuplotchannel,*) 'unset ytics'
        write(gnuplotchannel,*) 'plot "<(sed -n ''1,'//trim(adjustl(ntext))//&
                'p'' '//gridpath5//"heatmapline"//errorcheckfile//&
                ')" u (log10($2)-xmin):(log10($3)-ymin) w l lw 3 lc "black"'
        write(gnuplotchannel,*) 'set samples 100000'
        write(gnuplotchannel,*) 'alpha=',alpha_ratio**2
        write(gnuplotchannel,*) 'cA = 10**((xmin + (xmax-xmin)*1/5)*2) + '//&
                                     'alpha*10**((ymin + (ymax-ymin)*1/5)*2)'
        write(gnuplotchannel,*) 'cB = 10**((xmin + (xmax-xmin)*2/5)*2) + '//&
                                     'alpha*10**((ymin + (ymax-ymin)*2/5)*2)'
        write(gnuplotchannel,*) 'cC = 10**((xmin + (xmax-xmin)*3/5)*2) + '//&
                                     'alpha*10**((ymin + (ymax-ymin)*3/5)*2)'
        write(gnuplotchannel,*) 'cD = 10**((xmin + (xmax-xmin)*4/5)*2) + '//&
                                     'alpha*10**((ymin + (ymax-ymin)*4/5)*2)'
!       write(gnuplotchannel,FMT="(A,F15.11,A)") 'f(x) = (ymax-ymin)/2 + x/(', alpha_ratio, '*(10**(xmin/ymin)))'
!       write(gnuplotchannel,FMT="(A,F15.11,A)") 'f(x) = 0.5*log10(cA - 10**(2*(x+xmin))) - ymin'
        write(gnuplotchannel,FMT="(A,F15.11,A)") 'fA(x) = 0.5*log10((cA - 10**(2*(x+xmin)))/alpha) - ymin'
        write(gnuplotchannel,FMT="(A,F15.11,A)") 'fB(x) = 0.5*log10((cB - 10**(2*(x+xmin)))/alpha) - ymin'
        write(gnuplotchannel,FMT="(A,F15.11,A)") 'fC(x) = 0.5*log10((cC - 10**(2*(x+xmin)))/alpha) - ymin'
        write(gnuplotchannel,FMT="(A,F15.11,A)") 'fD(x) = 0.5*log10((cD - 10**(2*(x+xmin)))/alpha) - ymin'
!       write(gnuplotchannel,*) 'plot f(x) w l lw 2 lc "black"'
        write(gnuplotchannel,*) 'plot fA(x) w l lw 2 lc "black"'
        write(gnuplotchannel,*) 'plot fB(x) w l lw 2 lc "black"'
        write(gnuplotchannel,*) 'plot fC(x) w l lw 2 lc "black"'
        write(gnuplotchannel,*) 'plot fD(x) w l lw 2 lc "black"'

        write(gnuplotchannel,*) 'set xrange [1:',Ninterpolation,']'
        write(gnuplotchannel,*) 'set logscale y'
        write(gnuplotchannel,*) "unset lmargin"
        write(gnuplotchannel,*) "unset rmargin"
        write(gnuplotchannel,*) "unset bmargin"
        write(gnuplotchannel,*) "unset tmargin"
        write(gnuplotchannel,*) 'set xtics'
        write(gnuplotchannel,*) 'set ytics'
        write(gnuplotchannel,*) 'set xlabel "Ninterpolation"'

        write(gnuplotchannel,*) 'set size 0.5, 0.5'
        write(gnuplotchannel,*) 'set origin 0.5, 0.0'
        write(gnuplotchannel,*) 'set title "Title 2"'
        write(gnuplotchannel,*) 'set yrange [',min_rmsd1,':',max_rmsd1,']'
        write(gnuplotchannel,*) 'set ylabel "RMSD Variable 1"'
        write(gnuplotchannel,*) 'plot "<(sed -n ''1,'//trim(adjustl(ntext))//&
                'p'' '//gridpath5//"heatmapline"//errorcheckfile//&
                ')" u 1:2 w l lc "black" lw 3'

        write(gnuplotchannel,*) 'set size 0.5, 0.5'
        write(gnuplotchannel,*) 'set origin 0.5, 0.5'
        write(gnuplotchannel,*) 'set title "Title 3"'
        write(gnuplotchannel,*) 'set yrange [',min_rmsd2,':',max_rmsd2,']'
        write(gnuplotchannel,*) 'set ylabel "RMSD Variable 2"'
        write(gnuplotchannel,*) 'plot "<(sed -n ''1,'//trim(adjustl(ntext))//&
                'p'' '//gridpath5//"heatmapline"//errorcheckfile//&
                ')" u 1:3 w l lc "black" lw 3'
        close(gnuplotchannel)

        call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)
        end do

        call system("convert -background white -alpha remove -layers OptimizePlus "//&
                "-delay 20 -loop 0 "//gridpath4//"png/*.png "//&
                gridpath4//"heatmap_trajectory.gif")

        print *, "Finished making the gif!"
        call sleep(5)
        end if







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
        write(gnuplotchannel,FMT="(A,2F7.4,A,I0.5,A)") &
                'set output "'//gridpath4//"ratio", vals, '_', Ntrials, '.png"'
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
                'set output "'//gridpath4//'LocalErrorConvergence.png"'
        write(gnuplotchannel,*) 'set title "Local Error Convergence over a Trajectory"'
        write(gnuplotchannel,*) 'unset key'
        write(gnuplotchannel,FMT="(A,I5,A)") 'set label 1 "N = ', Nsamples, '" at screen 0.05,0.925'
        write(gnuplotchannel,FMT='(A,F7.4,A)') 'set label 2 "Threshhold = ',threshold_rmsd, &
                ' A" at screen 0.15,0.925'
        write(gnuplotchannel,FMT='(A,E16.8,A)') 'set label 3 "AlphaRatio = ',alpha_ratio, &
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

subroutine plotCheckstateEnergy(PNGname)
        use FUNCTIONS
        use PARAMETERS
        use PHYSICS
        use VARIABLES
        use ANALYSIS
        implicit none

        character(*),intent(in) :: PNGname
        integer :: step, frames
        integer :: Nsaved,Nwasted,Nrejected

        real(dp) :: min_rmsd, min_rmsd_prime
        real(dp) :: min_rmsd_previous
        real,dimension(Nvar) :: vals
        real(dp),dimension(6) :: rsv

        real(dp) :: U, KE
        real(dp) :: E, E_previous
        real(dp) :: avg_DE,max_absDE
        real(dp) :: avg_absDE,sd_absDE

        real(dp) :: DE_binwidth
        integer :: Nbins,DE_bin
        integer,allocatable :: DE_binning1(:), DE_binning2(:)

        integer :: i1, i2, i3
        integer :: n, m, iostate

        open(filechannel2,file=gridpath5//&
                "truncated"//checkstatefile)
        read(filechannel2,FMT="(1x,4(I8,1x))") &
                Nsaved,Nwasted,steps,Nrejected

        max_absDE = 0.0d0
        avg_absDE = 0.0d0
        avg_DE = 0.0d0
        frames = 1
        do
            read(filechannel2,iostat=iostate,FMT=*) &
                    i1,i2,i3,step,&
                    min_rmsd,min_rmsd_prime,&
                    vals(1),vals(2),U,KE
            if (iostate /= 0) exit

            if (step == 1) then
                E_previous = U + KE
                min_rmsd_previous = min_rmsd_prime
                cycle
            end if

            E = U + KE

            max_absDE = max(abs(E-E_previous),max_absDE)
            avg_absDE = avg_absDE + abs(E-E_previous)
            avg_DE = avg_DE + E - E_previous

            E_previous = E
            frames = frames + 1
        end do
        close(filechannel2)

        avg_absDE = avg_absDE / frames
        avg_DE = avg_DE / frames


        Nbins = 40
        DE_binwidth = 2 * max_absDE / Nbins
        allocate(DE_binning1(Nbins),DE_binning2(Nbins))
        DE_binning1 = 0
        DE_binning2 = 0

        if (gather_interpolation_flag) &
                open(filechannel1,file=gridpath5//interpolationfile)
        open(filechannel2,file=gridpath5//&
                "truncated"//checkstatefile)
        read(filechannel2,FMT="(1x,4(I8,1x))") &
                Nsaved,Nwasted,steps,Nrejected

        open(filechannel3,file=gridpath5//temporaryfile2)

        min_rmsd_previous = default_rmsd
        sd_absDE = 0.0d0
        do
            read(filechannel2,iostat=iostate,FMT=*) &
                    i1,i2,i3,step,&
                    min_rmsd,min_rmsd_prime,&
                    vals(1),vals(2),U,KE
            if (iostate /= 0) exit

            if (step == 1) then
                if (min_rmsd_previous < default_rmsd*0.9) then
                    if (gather_interpolation_flag) then
                        read(filechannel1,FMT=*) vals, i1, rsv
                    end if
                end if

                E_previous = U + KE
                min_rmsd_previous = min_rmsd_prime
                cycle
            end if

            E = U + KE
            sd_absDE = sd_absDE + (abs(E-E_previous) -&
                    avg_absDE)**2

            DE_bin = floor((max_absDE + E-E_previous)/&
                    DE_binwidth)
            if (DE_bin == 0) DE_bin = 1
            if (DE_bin > Nbins) DE_bin = Nbins

            if (min_rmsd_previous < default_rmsd*0.9) then
                DE_binning1(DE_bin) = DE_binning1(DE_bin) + 1

!               write(filechannel3,FMT=*) vals(1),vals(2),&
!                       Ninterpolation,largest_weighted_rmsd2,&
!                       largest_weighted_rmsd,candidate_rmsd,&
!                       min_rmsd,error1,error2
                if (gather_interpolation_flag) then
                    read(filechannel1,FMT=*) vals, i1, rsv
    
                    write(filechannel3,FMT=*) &
                            (E-E_previous) * RU_energy / eV, &
                            min_rmsd_prime, &
                            rsv(6) * RU_energy / eV
                end if
            else
                DE_binning2(DE_bin) = DE_binning2(DE_bin) + 1
            end if

            E_previous = E
            min_rmsd_previous = min_rmsd_prime
        end do

        if (gather_interpolation_flag) &
                close(filechannel1)
        close(filechannel2)
        close(filechannel3)

        sd_absDE = sqrt(sd_absDE / frames)

        open(filechannel2,file=gridpath5//temporaryfile1)
        do n = 1, Nbins
            write(filechannel2,FMT=*) &
                    (DE_binwidth * (n + 0.5) - max_absDE) *&
                    RU_energy / eV,&
                    DE_binning1(n), DE_binning2(n)
        end do
        close(filechannel2)
        deallocate(DE_binning1, DE_binning2)



        open(filechannel2,file=gridpath5//&
                "energyconservation.dat",position="append")
        write(filechannel2,FMT=*) avg_DE, &
                max_absDE, avg_absDE, sd_absDE
        close(filechannel2)



        if (.not.(gather_interpolation_flag)) then
        open(gnuplotchannel,file=gridpath5//gnuplotfile)
        write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
        write(gnuplotchannel,FMT="(A)") &
                'set output "'//gridpath4//PNGname//'_rewind.png"'
        write(gnuplotchannel,*) 'set title "Energy Conservation over a Trajectory"'
        write(gnuplotchannel,*) 'set xlabel "Time (fs)"'
        write(gnuplotchannel,*) 'xscale = ', dt
        write(gnuplotchannel,*) 'set ylabel "Energy (eV)"'
        write(gnuplotchannel,*) 'yscale = ', RU_energy / eV
        write(gnuplotchannel,*) 'plot "'//gridpath5//checkstatefile//&
                '" u (($4)*xscale):(($9)*yscale) w l lc "blue" t "U",\'
        write(gnuplotchannel,*) '     "'//gridpath5//checkstatefile//&
                '" u (($4)*xscale):(($10)*yscale) w l lc "red" t "KE",\'
        write(gnuplotchannel,*) '     "'//gridpath5//checkstatefile//&
                '" u (($4)*xscale):((($9)+($10))*yscale) w l lc "black" t "Total",\'
        write(gnuplotchannel,*) '     "'//gridpath5//checkstatefile//&
                '" u (($4)*xscale):($6>1.0?1/0:((($9)+($10))*yscale)) w p pt 7 ps 1 lc "red" t ""'
        close(gnuplotchannel)

        call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

        open(gnuplotchannel,file=gridpath5//gnuplotfile)
        write(gnuplotchannel,*) 'set term pngcairo size 1200,1200'
        write(gnuplotchannel,FMT="(A)") &
                'set output "'//gridpath4//PNGname//'_final.png"'
        write(gnuplotchannel,*) 'set title "Energy Conservation over a Trajectory"'
        write(gnuplotchannel,*) 'set tmargin at screen 0.90'
        write(gnuplotchannel,FMT="(A,I8,A,I8,A,I8,A,I8,A)") &
           'set label 1 "Nsaved:',Nsaved,&
           ' Nwasted:',Nwasted,&
           ' Ncomputed(rejected:isolated):',Nrejected,':',steps-Nsaved-Nrejected,&
           '" at screen 0.3,0.975'

        write(gnuplotchannel,*) 'set xlabel "Time (fs)"'
        write(gnuplotchannel,*) 'xscale = ', dt
        write(gnuplotchannel,*) 'set ylabel "Energy (eV)"'
        write(gnuplotchannel,*) 'yscale = ', RU_energy / eV
        write(gnuplotchannel,*) 'plot "'//gridpath5//"truncated"//checkstatefile//&
                '" u (($4)*xscale):(($9)*yscale) w l lc "blue" t "U",\'
        write(gnuplotchannel,*) '     "'//gridpath5//"truncated"//checkstatefile//&
                '" u (($4)*xscale):(($10)*yscale) w l lc "red" t "KE",\'
        write(gnuplotchannel,*) '     "'//gridpath5//"truncated"//checkstatefile//&
                '" u (($4)*xscale):((($9)+($10))*yscale) w l lc "black" t "Total",\'
        write(gnuplotchannel,*) '     "'//gridpath5//"truncated"//checkstatefile//&
                '" u (($4)*xscale):($6>1.0?1/0:((($9)+($10))*yscale)) w p pt 7 ps 1 lc "red" t ""'
        close(gnuplotchannel)

        call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

        else

        open(gnuplotchannel,file=gridpath5//gnuplotfile)
        write(gnuplotchannel,*) 'set term pngcairo size 1200,1800'
        write(gnuplotchannel,FMT="(A)") &
                'set output "'//gridpath4//PNGname//'_Distribution.png"'
        write(gnuplotchannel,*) 'set tmargin 0'
        write(gnuplotchannel,*) 'set bmargin 0'
        write(gnuplotchannel,*) 'set lmargin 1'
        write(gnuplotchannel,*) 'set rmargin 1'
        write(gnuplotchannel,*) 'unset title'
        write(gnuplotchannel,FMT="(A)") 'set multiplot layout '//&
                                '4,1 columnsfirst margins 0.1,0.95,.1,.9 spacing 0.1,0 '//&
                                'title "Energy Jump Distributions" font ",32"'
        write(gnuplotchannel,FMT='(A,E16.6,":",E16.6,A)') &
                'set xrange [',&
                (-max_absDE - 0.5 * DE_binwidth) * RU_energy / eV,&
                ( max_absDE + 0.5 * DE_binwidth) * RU_energy / eV,&
                ']'
        write(gnuplotchannel,*) 'set ytics font ",16"'
        write(gnuplotchannel,*) 'set xtics font ",16"'
        write(gnuplotchannel,*) 'set format x ""'
        write(gnuplotchannel,*) 'unset xlabel'

        write(gnuplotchannel,*) 'set ylabel "Error (eV/A)" font ",24"'
!       write(gnuplotchannel,FMT='(A,E16.6,":",E16.6,A)') &
!               'set yrange [',&
!               inner_threshold,&
!               threshold_rmsd,&
!               ']'
        write(gnuplotchannel,*) 'plot "'//gridpath5//temporaryfile2//&
                '" u 1:3 w points pt 1 ps 1 t "Accept Best"'
 
        write(gnuplotchannel,*) 'set ylabel "RMSD (A)" font ",24"'
        write(gnuplotchannel,FMT='(A,E16.6,":",E16.6,A)') &
                'set yrange [',&
                inner_threshold,&
                threshold_rmsd,&
                ']'
        write(gnuplotchannel,*) 'plot "'//gridpath5//temporaryfile2//&
                '" u 1:2 w points pt 1 ps 1 t "Accept Best"'

        write(gnuplotchannel,*) 'set style histogram clustered gap 1'
        write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'
        write(gnuplotchannel,*) 'set ylabel "Occurence" font ",24"'
        write(gnuplotchannel,*) 'set autoscale y'
        write(gnuplotchannel,*) 'set yrange [0:]'
        
        write(gnuplotchannel,*) 'plot "'//gridpath5//temporaryfile1//&
                '" u 1:2 w boxes t "Accept Best"'
        write(gnuplotchannel,*) 'set format x "%.1e"'
        write(gnuplotchannel,*) 'set xlabel "Energy Jump (eV)"'
        write(gnuplotchannel,*) 'plot "'//gridpath5//temporaryfile1//&
                '" u 1:3 w boxes t "Reject"'
        close(gnuplotchannel)

        call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

        end if


end subroutine plotCheckstateEnergy

subroutine plotEnergyConservationInformatics()
    use FUNCTIONS
    use PARAMETERS
    use PHYSICS
    use VARIABLES
    use ANALYSIS
    implicit none

    integer :: frames

!   real(dp) :: avg_DE,max_absDE
!   real(dp) :: avg_absDE,sd_absDE

    real(dp),dimension(4) :: DEs
    real(dp),dimension(4) :: max_DEs,min_DEs
    real(dp),dimension(4) :: DE_binwidths
    integer,dimension(4) :: DE_bins

    real(dp) :: avg_absDE1
    real(dp) :: max_DE2
    real(dp) :: avg_DE3
    real(dp) :: avg_DE4

    integer :: Nbins
    integer,allocatable :: DE_binning(:,:)

    integer :: n, m, iostate

    frames = 0
    max_DEs = -1.0d19
    min_DEs = 1.0d19

    avg_absDE1 = 0.0d0
    max_DE2 = -1.0d19
    avg_DE3 = 0.0d0
    avg_DE4 = 0.0d0

    open(filechannel2,file=gridpath5//&
            "energyconservation.dat")
    do 
!       read(filechannel2,iostat=iostate,FMT=*) avg_DE, &
!               max_absDE, avg_absDE, sd_absDE
        read(filechannel2,iostat=iostate,FMT=*) DEs
        if (iostate /= 0) exit

        do n = 1, 4
            max_DEs(n) = max(max_DEs(n),DEs(n))
            min_DEs(n) = min(min_DEs(n),DEs(n))
        end do

        avg_absDE1 = avg_absDE1 + abs(DEs(1))
        max_DE2 = max(max_DE2,DEs(2))
        avg_DE3 = avg_DE3 + DEs(3)
        avg_DE4 = avg_DE4 + DEs(4)

        frames = frames + 1
    end do
    close(filechannel2)

    avg_absDE1 = (avg_absDE1 / frames)  * RU_energy / eV
    max_DE2 = max_DE2 * RU_energy / eV
    avg_DE3 = (avg_DE3 / frames) * RU_energy / eV
    avg_DE4 = (avg_DE4 / frames) * RU_energy / eV

    Nbins = 40
    allocate(DE_binning(4,Nbins))
    DE_binning = 0

    DE_binwidths = (max_DEs - min_DEs) / Nbins

    open(filechannel2,file=gridpath5//&
            "energyconservation.dat")
    do 
        read(filechannel2,iostat=iostate,FMT=*) DEs
        if (iostate /= 0) exit

        DE_bins = floor((DEs - min_DEs)/DE_binwidths)

        do n = 1, 4
            if (DE_bins(n) == 0) DE_bins(n) = 1
            if (DE_bins(n) > Nbins) DE_bins = Nbins

            DE_binning(n,DE_bins(n)) = &
                DE_binning(n,DE_bins(n)) + 1
        end do
    end do
    close(filechannel2)

    open(filechannel2,file=gridpath5//temporaryfile1)
    do n = 1, Nbins
        write(filechannel2,FMT=*) &
                (min_DEs + (n-0.5) * DE_binwidths) * RU_energy / eV,&
                DE_binning(:,n)
    end do
    close(filechannel2)
    

    open(gnuplotchannel,file=gridpath5//gnuplotfile)
    write(gnuplotchannel,*) 'set term pngcairo size 4200,1200'
    write(gnuplotchannel,FMT="(A)") &
            'set output "'//gridpath4//'EnergyConservationInformatics.png"'
    write(gnuplotchannel,*) 'unset title'
    write(gnuplotchannel,*) 'set multiplot layout 1,4 '//&
            'title "Energy Conservation Informatics" font ",32"'
    write(gnuplotchannel,*) 'set ylabel "Occurence" font ",24"'
    write(gnuplotchannel,*) 'set yrange [0:]'
    write(gnuplotchannel,*) 'set ytics font ",16"'
    write(gnuplotchannel,*) 'set xtics font ",16"'
!   write(gnuplotchannel,*) 'set format x "%e"'
!   write(gnuplotchannel,*) 'set format x "%.1t*10^%+03T"'
    write(gnuplotchannel,*) 'set format x "%.1e"'
    write(gnuplotchannel,*) 'unset key'
    write(gnuplotchannel,*) 'set style histogram clustered gap 1'
    write(gnuplotchannel,*) 'set style fill solid 1.0 noborder'

    write(gnuplotchannel,*) 'set xlabel "Total Energy Drift (eV)" font ",24"'
    write(gnuplotchannel,FMT="(A)") 'unset arrow'
    write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
            'set arrow from ', avg_absDE1,&
            ',graph 0 to ', avg_absDE1, ', graph 1 nohead front '//&
            'lw 2 lc rgb "black"'
    write(gnuplotchannel,FMT='(A,E16.6,":",E16.6,A)') &
            'set xrange [',&
            (min_DEs(1) - 0.5 * DE_binwidths(1)) * RU_energy / eV,&
            (max_DEs(1) + 0.5 * DE_binwidths(1)) * RU_energy / eV,&
            ']'
    write(gnuplotchannel,*) 'plot "'//gridpath5//temporaryfile1//&
            '" u 1:5 w boxes'

    write(gnuplotchannel,*) 'set xlabel "Maximum Energy Jump (eV)" font ",24"'
    write(gnuplotchannel,FMT="(A)") 'unset arrow'
    write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
            'set arrow from ', max_DE2,&
            ',graph 0 to ', max_DE2, ', graph 1 nohead front '//&
            'lw 2 lc rgb "black"'
    write(gnuplotchannel,FMT='(A,E16.6,":",E16.6,A)') &
            'set xrange [',&
            (min_DEs(2) - 0.5 * DE_binwidths(2)) * RU_energy / eV,&
            (max_DEs(2) + 0.5 * DE_binwidths(2)) * RU_energy / eV,&
            ']'
    write(gnuplotchannel,*) 'plot "'//gridpath5//temporaryfile1//&
            '" u 2:6 w boxes'

    write(gnuplotchannel,*) 'set xlabel "Average Energy Jump (eV)" font ",24"'
    write(gnuplotchannel,FMT="(A)") 'unset arrow'
    write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
            'set arrow from ', avg_DE3,&
            ',graph 0 to ', avg_DE3, ', graph 1 nohead front '//&
            'lw 2 lc rgb "black"'
    write(gnuplotchannel,FMT='(A,E16.6,":",E16.6,A)') &
            'set xrange [',&
            (min_DEs(3) - 0.5 * DE_binwidths(3)) * RU_energy / eV,&
            (max_DEs(3) + 0.5 * DE_binwidths(3)) * RU_energy / eV,&
            ']'
    write(gnuplotchannel,*) 'plot "'//gridpath5//temporaryfile1//&
            '" u 3:7 w boxes'

    write(gnuplotchannel,*) 'set xlabel "Standard Deviation of Energy Jumps (eV)" font ",24"'
    write(gnuplotchannel,FMT="(A)") 'unset arrow'
    write(gnuplotchannel,FMT="(A,E16.6,A,E16.6,A)")&
            'set arrow from ', avg_DE4,&
            ',graph 0 to ', avg_DE4, ', graph 1 nohead front '//&
            'lw 2 lc rgb "black"'
    write(gnuplotchannel,FMT='(A,E16.6,":",E16.6,A)') &
            'set xrange [',&
            (min_DEs(4) - 0.5 * DE_binwidths(4)) * RU_energy / eV,&
            (max_DEs(4) + 0.5 * DE_binwidths(4)) * RU_energy / eV,&
            ']'
    write(gnuplotchannel,*) 'plot "'//gridpath5//temporaryfile1//&
            '" u 4:8 w boxes'

    call system(path_to_gnuplot//"gnuplot < "//gridpath5//gnuplotfile)

end subroutine plotEnergyConservationInformatics



end module runTrajectory
