!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PROGRAM
!               runTrajectory
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!		This program simulates an MD collision of two molecules
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
!		getVar3				VARIABLES
!		getVar4				VARIABLES
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
!       FILES	                          FILETYPE	SUBROUTINE
!
!		gridpath0//trajectory		XYZ		checkMultipleTrajectories
!		gridpath1//trajectory		XYZ		checkTrajectory
!		gridpath1//checkstatefile	DAT		checkTrajectory
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
!	INPUT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	OUTPUT				KIND				DESCRIPTION
!
!		velocityH1			REAL(DP)			The velocity of the incoming H/H2 molecule
!										at the start of the trajectory
!		velocityH2			REAL(DP)			The velocity of the incoming H/H2 molecule
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
!	OUTPUT				KIND				DESCRIPTION
!
!		velocityH1			REAL(DP)			The velocity of the incoming H/H2 molecule
!										at the start of the trajectory
!		velocityH2			REAL(DP)			The velocity of the incoming H/H2 molecule
!										at the end of the trajectory (post-collision)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES	                          FILETYPE			SUBROUTINE
!
!		gridpath0//trajectory		XYZ				checkMultipleTrajectories
!		gridpath1//trajectory		XYZ				checkTrajectory
!		gridpath1//checkstatefile	DAT				checkTrajectory
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine addTrajectory(velocityH1,velocityH2)
	use VARIABLES
	use PHYSICS
	use PARAMETERS
	use interactSingleGrid
        implicit none

	!Coordinates, Velocities, and Variables
	real(dp), dimension(3,Natoms) :: coords,gradient,velocities,approx_gradient
	real(dp), dimension(Nvar) :: vals
	real(dp), dimension(3),intent(out) :: velocityH1, velocityH2

	call system_clock(c1,count_rate=cr)
	system_clock_rate = 1.0/real(cr)

        !Initialize the scene
        call InitialSetup3(coords,velocities)
	Norder1 = 0

	!velocityH is just the velocity from the get-go
        !velocityH1 = velocities(:,1) - (velocities(:,1)+velocities(:,2)+velocities(:,3))/3
	velocityH1 = velocities(:,1)

	call getVar3(coords,Natoms,vals(1))
	call getVar4(coords,Natoms,vals(2))

        !Accelerate the velcocities for a half step (verlet)
        call Acceleration(vals,coords,gradient)
	call addState(vals,coords,gradient)

	velocities = velocities + 0.5d0 * gradient

        do steps = 1, Nsteps


		!Just for bug-testing
                if (.false.) then !(modulo(steps,50) == 0) then
                        open(filechannel1,file=gridpath0//trajectoryfile,position="append")
                        write(filechannel1,'(I1)') 3
                        write(filechannel1,*) ""
                        write(filechannel1,'(A1,3F10.6)') 'H',&
                                coords(1,1), coords(2,1), coords(3,1)
                        write(filechannel1,'(A1,3F10.6)') 'H',&
                                coords(1,2), coords(2,2), coords(3,2)
                        write(filechannel1,'(A1,3F10.6)') 'H',&
                                coords(1,3), coords(2,3), coords(3,3)
			close(filechannel1)
                end if
 
                !Just to see progress, print something out every 500 steps
                if (modulo(steps,500) == 1) then      
 			if ((vals(1)>max_var1) .or. (vals(2)>max_var2)) then
 				exit
 			end if
                endif

                !Update the coordinates with the velocities
		coords = coords + dt * velocities

		call getVar3(coords,Natoms,vals(1))
		call getVar4(coords,Natoms,vals(2))


                !Accelerate and update velocities
                call Acceleration(vals,coords,gradient)
        	call addState(vals,coords,gradient)

		!If there are too many subdivisions and counter0 will get out-of-bounds
		!we would have to call this to exit
                if (header_max_flag) exit

		velocities = velocities + gradient

        end do

        !velocityH2 = velocities(:,1) - (velocities(:,1)+velocities(:,2)+velocities(:,3))/3
	velocityH2 = velocities(:,1)

end subroutine addTrajectory



subroutine checkTrajectory(velocityH1,velocityH2)
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
        real(dp), dimension(3),intent(out) :: velocityH1,velocityH2

        !Various other variables
        real(dp) :: U, KE
        real(dp) :: min_rmsd,min_rmsd_prime
        integer :: TimeA,TimeB
        integer :: number_of_frames,order,neighbor_check

        !Initialize the scene
        call InitialSetup3(coords,velocities)

        !velocityH1 = velocities(:,1) - (velocities(:,1)+velocities(:,2)+velocities(:,3))/3
        velocityH1 = velocities(:,1)

        !Accelerate the velcocities for a half step (verlet)

        call getVar3(coords,Natoms,vals(1))
        call getVar4(coords,Natoms,vals(2))

        call Acceleration(vals,coords,gradient)

         open(filechannel1,file=gridpath1//trajectoryfile)
         write(filechannel1,'(I1)') 3
         write(filechannel1,*) ""
         write(filechannel1,'(A1,3F10.6)') 'H',&
               coords(1,1), coords(2,1), coords(3,1)
         write(filechannel1,'(A1,3F10.6)') 'H',&
               coords(1,2), coords(2,2), coords(3,2)
         write(filechannel1,'(A1,3F10.6)') 'H',&
               coords(1,3), coords(2,3), coords(3,3)
        close(filechannel1)

        velocities = velocities + 0.5d0 * gradient

        !Keep track of the time
        TimeA = time()
open(filechannel2,file=gridpath1//checkstatefile)
        do steps = 1, Nsteps

                !Every 50 frames, print to an xyz file for visualization
                 if (modulo(steps,10) == 0) then
                         open(filechannel1,file=gridpath1//trajectoryfile,position="append")
                         write(filechannel1,'(I1)') 3
                         write(filechannel1,*) ""
                         write(filechannel1,'(A1,3F10.6)') 'H',&
                                 coords(1,1), coords(2,1), coords(3,1)
                         write(filechannel1,'(A1,3F10.6)') 'H',&
                                 coords(1,2), coords(2,2), coords(3,2)
                         write(filechannel1,'(A1,3F10.6)') 'H',&
                                 coords(1,3), coords(2,3), coords(3,3)
                        close(filechannel1)
                 end if

                !Just to see progress, print something out every 500 steps
                if (modulo(steps,500) == 1) then
                        if ((vals(1)>max_var1) .or.&
                            (vals(2)>max_var2)) then
                                exit
                        end if
                endif

                !Upate the coordinates with the velocities
                coords = coords + dt * velocities

                !Get the variables
                call getVar3(coords,Natoms,vals(1))
                call getVar4(coords,Natoms,vals(2))

                !Check for similar frames
                min_rmsd_prime = .200100d0
                call checkState(coords,approx_gradient,min_rmsd_prime,&
                                number_of_frames,order,neighbor_check)
                min_rmsd = min_rmsd_prime

                if (force_Neighbors) then
                call checkState(coords,approx_gradient,min_rmsd)
                end if

                !Update the gradient
                if ((min_rmsd .ge. threshold_rmsd) .or. (reject_flag)) then
                        call Acceleration(vals,coords,gradient)
                else
                        gradient = approx_gradient
                end if

	        U = MorsePotential(coords(:,1),coords(:,2))
	        U = U + MorsePotential(coords(:,1),coords(:,3))
	        U = U + HOPotential(coords(:,2),coords(:,3))
	        KE = KineticEnergy(velocities(:,1))
	        KE = KE + KineticEnergy(velocities(:,2))
	        KE = KE + KineticEnergy(velocities(:,3))
		write(filechannel2,FMT=*) number_of_frames,order,neighbor_check,steps,&
                                          min_rmsd,min_rmsd_prime,vals(1),vals(2),U,KE

                !Update the velocities
                velocities = velocities + gradient

        end do
close(filechannel2)

        !velocityH2 = velocities(:,1) - (velocities(:,1)+velocities(:,2)+velocities(:,3))/3
        velocityH2 = velocities(:,1)

end subroutine checkTrajectory




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       SUBROTUINE checkMultipleTrajectories
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!INPUT: real :: initial_bond_distance           "the length of the H2 bond initially"
!               initial_rotational_speed        "the rotational speed of the H2"
!               initial_rotation_angle          "the direction the H2 is spinning (relative)"
!
!       real :: initial_bond_angle1             "the angle the H2 bond makes with the x-axis"
!               initial_bond_angle2             "the angle the H2 bond makes with the y-z plane"
!
!       logical            :: force_Neighbors           "the choice to check adjacent cells"
!       integer            :: Ngrid_total                       "the number of grids to check"
!       int,dim(Ngrid_total) :: filechannels            "the files that we write RMSD calls to"
!       character(*)       :: path_to_directory         "the path to the directory with the grids"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Simulates a H --> H2 collision with the above parameters
!       Uses physical constants and parameters as supplied by f2_physics_parameters.f90
!       Checks a folder full of grids of files formatted by f2_parameters.f90
!       Checks each individual frame according to subroutines in checkCells5.f90
!       Makes use of a cross product which is supplied by f1_functions.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine checkMultipleTrajectories(filechannels,velocityH1,velocityH2)
        use PHYSICS
        use PARAMETERS
        use VARIABLES
	use ANALYSIS
        use interactMultipleGrids
        implicit none

        !Coordinates, Velocities, and Variables
        real(dp),dimension(3,Natoms) :: coords,velocities
        real(dp),dimension(3,Natoms) :: gradient,approx_gradient
        real(dp),dimension(Nvar) :: vals
        real(dp),dimension(3),intent(out) :: velocityH1, velocityH2
        real(dp),dimension(3) :: velocity_rotation,velocity_vibration,velocity_translation

        !Grid Parameters
        integer,dimension(Ngrid_total),intent(in) :: filechannels

        !Various other variables
        real(dp) :: U, KE
        real(dp) :: min_rmsd,min_rmsd_prime
        integer :: TimeA,TimeB
        integer :: number_of_frames,order,neighbor_check

        !Initialize the scene
        call InitialSetup3(coords,velocities)

        !velocityH1 = velocities(:,1) - (velocities(:,1)+velocities(:,2)+velocities(:,3))/3
        velocityH1 = velocities(:,1)

        open(filechannel1,file=gridpath0//trajectoryfile)
        write(filechannel1,'(I1)') 3
        write(filechannel1,*) ""
        write(filechannel1,'(A1,3F10.6)') 'H',&
              coords(1,1), coords(2,1), coords(3,1)
        write(filechannel1,'(A1,3F10.6)') 'H',&
              coords(1,2), coords(2,2), coords(3,2)
        write(filechannel1,'(A1,3F10.6)') 'H',&
              coords(1,3), coords(2,3), coords(3,3)
        close(filechannel1)

        !Calculate the variables for use in acceleration
        call getVar3(coords,Natoms,vals(1))
        call getVar4(coords,Natoms,vals(2))

        !Accelerate the velcocities for a half step (verlet)
        call Acceleration(vals,coords,gradient)

        velocities = velocities + 0.5d0 * gradient

        !Keep track of the time
        TimeA = time()
        do steps = 1, Nsteps

                if (modulo(steps,50) == 1) then
                        !Just to see progress, print something out every 500 steps
                        !If we are too far out, stop the simulation
                
                        open(filechannel1,file=gridpath0//trajectoryfile,position="append")
                        write(filechannel1,'(I1)') 3
                        write(filechannel1,*) ""
                        write(filechannel1,'(A1,3F10.6)') 'H',&
                              coords(1,1), coords(2,1), coords(3,1)
                        write(filechannel1,'(A1,3F10.6)') 'H',&
                              coords(1,2), coords(2,2), coords(3,2)
                        write(filechannel1,'(A1,3F10.6)') 'H',&
                              coords(1,3), coords(2,3), coords(3,3)
                        close(filechannel1)
                        if (modulo(steps,500) == 1) then
        
                                if ((vals(1)>max_var1) .or.&
                                    (vals(2)>max_var2)) then
                                        exit
                                end if
                        endif
                end if

                !Upate the coordinates with the velocities
                coords = coords + dt * velocities

                !Get the variables
                call getVar3(coords,Natoms,vals(1))
                call getVar4(coords,Natoms,vals(2))

                min_rmsd = .200100d0
                !Check for similar frames
                call checkState(vals,coords,approx_gradient,min_rmsd,&
                         filechannels,number_of_frames,order,neighbor_check)

                !Update the gradient
                if ((min_rmsd .ge. threshold_RMSD).or.(reject_flag)) then
                        call Acceleration(vals,coords,gradient)
                else
                        gradient = approx_gradient
                end if

                !Update the velocities
                velocities = velocities + gradient

        end do

        !velocityH2 = velocities(:,1) - (velocities(:,1)+velocities(:,2)+velocities(:,3))/3
        velocityH2 = velocities(:,1)

end subroutine checkMultipleTrajectories





subroutine Acceleration(vals,coords,gradient)
	use PARAMETERS
        use PHYSICS
        implicit none
        real(dp), dimension(Nvar), intent(in) :: vals
        real(dp), dimension(3,Natoms), intent(in) :: coords
        real(dp), dimension(3,Natoms), intent(out) :: gradient
        integer :: i, start_index1, start_index2
        integer :: index1, index2, bond_index1, bond_index2

        gradient = 0.0d0

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                 NON-BONDED INTERACTIONS
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        start_index1 = 1
        start_index2 = 2
        do i = 1, Nbonds
                bond_index1 = BONDING_DATA(i,1)
                do index1 = start_index1, bond_index1
                        if (index1 == bond_index1) then
                                bond_index2 = BONDING_DATA(i,2) - 1
                        else
                                bond_index2 = Natoms
                        end if

                        do index2 = start_index2, bond_index2
                                !Remark: the optional 5th argument is the distance between
                                !        coord1 and coord2 (so it doesn't have to recalculate)
                                call NonBondedForce(coords(:,index1),coords(:,index2),&
                                                gradient(:,index1),gradient(:,index2),vals(index2-1))
                        end do

                        start_index2 = index1 + 2
                end do

                start_index1 = bond_index1
                start_index2 = bond_index2 + 2
        end do
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

        do i = 1, Nbonds
                index1 = BONDING_DATA(i,1)
                index2 = BONDING_DATA(i,2)
                call BondedForce(coords(:,index1),coords(:,index2),&
                             gradient(:,index1),gradient(:,index2))
        end do

	return

end subroutine Acceleration


end module runTrajectory
