
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       MODULE
!               interactMultipleGrids
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This module governs how a frame interacts with a grid
!
!               A frame can be added to the grid or it can be
!               checked alongside the grid (for the closest frame)
!
!               The grid has a particular file formatting system
!               which stores frames according to the variables they are
!               associated with (ex. r1, r2); this formatting is
!               initialized in PARAMETERS
!
!               The grid also has an internal file counting system
!               which keeps track of the size of files, whether they
!               have children or not, and which children have a file
!               associated with them; this system is initialized in
!               PARAMETERS and reset whenever a grid is made 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILECHANNELS                    ACTION
!
!               FILECHANNELS(2:Ngrid_total+1)   WRITE
!               FILECHANNELS(1)                 OPEN, READ, CLOSE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINES                     ARGUMENTS               KIND
!
!		checkState                      vals                            intent(in),REAL(DP),DIM(Nvar)
!                                               coords                          intent(in),REAL(DP),DIM(3,Natoms)
!                                               gradient                        intent(inout),REAL(DP),DIM(3,Natoms)
!                                               min_rmsd                        intent(inout),REAL(DP)
!                                               filechannels                    intent(in),INTEGER,DIM(Ngrid_total+1)
!
!                                               number_of_frames                intent(out),INTEGER
!                                               order                           intent(out),INTEGER
!                                               neighbor_check                  intent(out),INTEGER
!
!		getNeighbors                    scaling                         intent(in),REAL
!                                               multiplier                      intent(in),REAL
!                                               FMTorder                        intent(in),CHARACTER(6)
!                                               index                           intent(in),INTEGER
!                                               filechannel                     intent(in),INTEGER
!
!                                               coords                          intent(in),REAL(DP),DIM(3,Natoms)
!                                               min_rmsd                        intent(inout),REAL(DP)
!                                               gradient                        intent(inout),REAL(DP),DIM(3,Natoms)
!                                               U                               intent(inout),REAL(DP),DIM(3,3)
!                                               number_of_frames                intent(inout),INTEGER
!
!		getRMSD_dp                      filechannel_thread              intent(in),INTEGER
!                                               subcell                         intent(in),CHARACTER(*)
!
!                                               coords                          intent(in),REAL(DP),DIM(3,Natoms)
!                                               population                      intent(out),INTEGER
!                                               min_rmsd                        intent(inout),REAL(DP)
!                                               gradient                        intent(inout),REAL(DP),DIM(3,Natoms)
!                                               U                               intent(inout),REAL(DP),DIM(3,3)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       CALLS                           MODULE
!
!               getVarMaxMin                    VARIABLES
!
!               getNeighbors                    interactMultipleGrids
!               getRMSD_dp                      interactMultipleGrids
!
!               rmsd                            ls_rmsd_original
!               matmul                          INTRINSIC
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE      SUBROUTINE                      FMT
!
!               gridpath2//                     DAT             getRMSD_dp              FMT7,FMT3
!                 trim(subcell).dat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module interactMultipleGrids
use ANALYSIS
use VARIABLES
use PARAMETERS
implicit none

!Whereas the buffer1 (the array holding information on
!cells of order = 1) may be large

real(dp),allocatable :: valsbuffer1(:,:)
real(dp),allocatable :: coordsbuffer1(:,:,:)
real(dp),allocatable :: gradientbuffer1(:,:,:)
real(dp),allocatable :: potEbuffer1(:)
real(dp),allocatable :: Ubuffer1(:,:,:)
integer,allocatable :: Ntrajbuffer1(:)

!For the memory buffer

integer,dimension(Nvar) :: previous_var_index
!integer,allocatable :: vals_hash(:,:)

integer,allocatable :: populationbuffer2(:)
real(dp),allocatable :: valsbuffer2(:,:,:)
integer,allocatable :: Ntrajbuffer2(:,:)
real(dp),allocatable :: coordsbuffer2(:,:,:,:)
real(dp),allocatable :: gradientbuffer2(:,:,:,:)
real(dp),allocatable :: potEbuffer2(:,:)

integer,allocatable :: temppopulationbuffer2(:)
real(dp),allocatable :: tempvalsbuffer2(:,:,:)
integer,allocatable :: tempNtrajbuffer2(:,:)
real(dp),allocatable :: tempcoordsbuffer2(:,:,:,:)
real(dp),allocatable :: tempgradientbuffer2(:,:,:,:)
real(dp),allocatable :: temppotEbuffer2(:,:)

!Arrays for interpolation

real(dp),allocatable :: inputCLS(:,:)

!In the latter case, we need to keep track of how large
!the buffer gets, and increase it if it gets too large

integer :: buffer1_size
integer :: buffer2_size

integer,dimension(Norder_max+1) :: local_frame_count
integer :: number_of_cells

!Variables related to interpolation

integer :: Ninterpolation

!Other global variables to clean things up

integer :: Totalnumber_of_frames
integer :: Norder
integer :: var_filechannel
real(dp) :: candidate_rmsd
real(dp),dimension(3,Natoms) :: candidate_gradient

integer :: Naccept




contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      checkState_PCM FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       INPUT: real(dp) dim(Nvar), vals
!              real(dp) dim(3,Natoms), coords
!              integer dim(Ngrid_total+1), filechannels
!      IN/OUT: real(dp) min_rmsd
!      OUTPUT: real(dp) dim(3,Natoms), gradient
!              integer number_of_frames (optional)
!              integer order (optional)
!              integer neighbor_check (optional)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This is a function that, given a frame, can
!       look through the library and retrieve its
!       neighboring frames; it then uses these frames
!       to interpolate an energy gradient. Meta-data
!       is also retrieved, such as the number of
!       frames searched, the order of the frames, and
!       how many cells were searced. There are also
!       several global variables that are updated
!       which relate to the error in interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine checkState_PCM(vals,coords,gradient,&
                min_rmsd,filechannels,&
                number_of_frames,order,neighbor_check)
use ls_rmsd_original
use ANALYSIS
use VARIABLES
use PARAMETERS
use FUNCTIONS
use SIMILARITY
implicit none
integer :: i,j,k,l
integer,intent(out),optional :: order,number_of_frames,neighbor_check
integer :: population
integer :: chosen_index
integer :: OMP_GET_THREAD_NUM
logical :: stop_flag
real(dp), dimension(Nvar), intent(in) :: vals
real(dp), dimension(3,Natoms), intent(in) :: coords
real(dp), dimension(3,Natoms) :: new_coords
real(dp), dimension(3,Natoms), intent(out) :: gradient
real(dp), intent(inout) :: min_rmsd
integer, dimension(1+Ngrid_total),intent(in) :: filechannels
character(Ngrid_text_length) :: Ngrid_text
character(5) :: variable_length_text
real(dp), dimension(3,3) :: new_U
real(dp),dimension(NSIs) :: new_SI

real(dp),dimension(Nvar) :: var_cell
integer,dimension(Nvar,Norder_max+1) :: var_index
integer,dimension(Nvar,Norder_max+1) :: var_index_diff

integer  :: largest_rmsd_error_index
real(dp) :: largest_rmsd_error
real(dp),dimension(Natoms,Natoms) :: temp_CM

! Initializing the trajectory numbers to all be zero
! (assuming the trajectory number can't naturally be
!  zero) turns out to be efficient for the diversity
! checking
Ntrajbuffer1 = 0

! The filechannel which does most of the heavy duty
! (reading and writing to files) is supplied as the
! first filechannel from the arguments
var_filechannel = filechannels(1)

! To calculate the similarity identifiers we must
! set the target coordinates
call setTarget(coords)

! We start off with zero frames having been checked
Totalnumber_of_frames = 0

! We start off assuming we are checking a first-
! order cell (order = 0)
order = 0

! We also start off assuming we have checked zero
! neighbors
neighbor_check = 0

! We initialize with zero non-empty cells having
! been checked
number_of_cells = 0

! The local_frame_count keeps track of the number
! of frames in the cell of the target (in case we
! need to divy-up)
local_frame_count = 0

! We start off with having zero frames suitable
! for interpolation
Ninterpolation = 0

! We may need to stop checking cells before the
! full search ends
stop_flag = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 SUBCELL TARGETING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Retrieve the index of each variable with
! respect to the grid and the real number
! (rounded) that represents that index
do i = 1, Nvar
    ! Repeat this for however many orders of cells
    ! deep we have been instructed to go
    do j = 1, Norder_max + 1
        var_index(i,j) = int(vals(i) * divisor(i,j))
    end do
end do

! To access the data in the memory buffer loaded
! from the previous cell search, we must subtract
! the var_indexes of the current and previous
! frames and supply that to shiftMemory
! (Note: we are using the first order (parent)
!        by default)
call shiftMemory((var_index(:,Norder_order(1)+1) -&
                 previous_var_index))

! And the previous var_index for the next frame
! will be the current var_index
previous_var_index = var_index(:,Norder_order(1)+1)

! Now, we start iterating over the grids;
! usually we are only looking at one grid
do k = 1, Ngrid_total
    
    ! Streamline the process by storing the
    ! path to the grid in a string gridpath3
    write(variable_length_text,FMT=FMT5_variable)&
            Ngrid_text_length
    write(Ngrid_text,FMT="(I0."//&
            trim(adjustl(variable_length_text))//&
            ")") k
    gridpath3 = gridpath0//Ngrid_text//"/grid/"
    
    ! The way we check, we check by order first
    do l = 1, Norder_max+1
    
        ! The user can specify in what order to
        ! check the orders through this array
        ! (so we can look at children cells
        !  before parent cells, for instance)
        Norder = Norder_order(l)
        
        ! Read the frames in the cells and
        ! process their RMSDs

        ! If we are adding to grid no. K then the
        ! user can specify that through the
        ! variable grid_addition
        if (k == grid_addition) then

            ! If this is our first time searching
            ! then of course we must look through
            ! the cell; we also look through the
            ! cell if we have not yet found any
            ! frames
            if ((l==1) .or. &
                (Totalnumber_of_frames == 0)) then
                call getRMSD_1_PCM(&
                        var_index(:,Norder+1),&
                        population)

            ! To add a frame to the grid, we need
            ! population data which we do through
            ! a "dummy" search (getRMSD_2)
            else
                call getRMSD_2(&
                        var_index(:,Norder+1),&
                        population)
            end if

            ! And here we can specify how many
            ! frames we found in this cell (in
            ! case we need to divy it up)
            local_frame_count(Norder+1) = &
                    population

        ! If we have no grid addition:
        else

            ! We follow a similar process
            if ((l==1) .or. &
                (Totalnumber_of_frames == 0)) then
                call getRMSD_1_PCM(&
                        var_index(:,Norder+1),&
                        population)

            ! But no need to do a dummy search
            ! if we don't need to
            else
            end if
        end if
        
        ! If the cell we are looking at is
        ! overcrowded then that means it has
        ! already been subdivided (ex. parent
        ! cell divies up to children cell)
        if (population >= var_overcrowd(Norder+1)) then

            ! In this case, we want to go deeper
            ! (ex. look at the children cells,
            !      not the parent cells) so long
            ! as we can go deeper
            if (Norder < Norder_max) cycle
        end if

        ! We keep track of how many frames we
        ! have had to look through
        Totalnumber_of_frames = &
            Totalnumber_of_frames + population
        
        ! If the cell is unpopulated or a certain flag is true,
        ! then we go ahead and look at the neighbors of this cell
        if ((force_Neighbors) .or. (population == 0)) then
        
            ! Integer i keeps track of how far away from the original
            ! cell we are; we look at cells on the 'diamond'
            ! surrounding the original cell
            do i = 1, subcellsearch_max(Norder+1)
        
                ! Initialize this to zero
                var_index_diff = 0
        
                ! And search through the required cells
                call getRelativeIndex_PCM(1,var_index(:,Norder+1),i,&
                                      var_index_diff,stop_flag)
                
                ! If we have found a frame, then stop_flag
                ! is true and we can exit. However, if we
                ! are forced to we can continue to search
                if ((stop_flag) .and. (.not. force_Neighbors)) exit
            end do
        end if
        

        ! If we found a non-empty cell then our search has
        ! been over this particular order and does not need
        ! to go over other orders
        ! (ex. if we have found at least one frame in the
        !      parent cells, then no need to look at the
        !      children cells)
        if (Totalnumber_of_frames > 0) then
                order = order + Norder
                exit
        end if
    
    end do

    ! Record the SI encountered here for later analysis:
    ! particularly, percent-RMSD graphs
    trajRMSDbuffer(k,Naccept+1) = SIbuffer1(Nsort,1)
    
end do


! Keep track of how many frames we have
! searched
number_of_frames = Totalnumber_of_frames

! Keep track of how many non-empty cells
! we have searched
neighbor_check = number_of_cells

! Keep track of what order (ex. parent or
! child) we have ultimately searched
Norder_total(1+order/Ngrid_total) = Norder_total(1+order/Ngrid_total) + 1

! If we have found no suitable frames
! then just return; the global variable
! candidate_rmsd is given the default
! original rmsd called min_rmsd
if (Ninterpolation == 0) then
    candidate_rmsd = min_rmsd
    return
end if

! If we have more than one frame then we
! can interpolate
if ((interpolation_flag).and.(Ninterpolation > 1)) then

    ! This subroutine actually does the
    ! interpolation
    call getInterpolatedGradient(new_coords,gradient)

    ! Calculate the SIs of the interpolated
    ! frame and energy gradient
    call getSIs(coords+new_coords,&
                new_coords,new_U,new_SI)

    ! Choose the minimum of the SIs; this
    ! assumes SIbuffer1 is sorted so that
    ! the minimum value is first
    min_SIs = SIbuffer1(:,1)

    ! What we call the interpolated_SIs
    ! are the SIs we just calculated
    interpolated_SIs = new_SI

else

    ! If we are not interpolating then
    ! we take the "best" frame; this also
    ! assumes all of buffer1 is sorted so
    ! that the best frame is first
    gradient = matmul(Ubuffer1(:,:,1),&
                      gradientbuffer1(:,:,1))

!   largest_weighted_rmsd = min_rmsd
!   largest_weighted_rmsd2 = min_rmsd**2

!   interpolated_CMdiff = SIbuffer1(2,1)

    ! Choose the minimum of the SIs; this
    ! assumes SIbuffer1 is sorted so that
    ! the minimum value is first
    min_SIs = SIbuffer1(:,1)

    ! What we call the interpolated_SIs
    ! would be those of the best frame
    interpolated_SIs = SIbuffer1(:,1)

    ! We must manually calculate these
    ! global variables (this is usually
    ! done behind-the-scenes in
    ! getInterpolatedGradient)
    largest_weighted_SIs = SIbuffer1(:,1)
    largest_weighted_SIs2 = SIbuffer1(:,1)**2
    
end if

! What we call min_rmsd is the the
! SI of the interpolated frame; if we
! sort by RMSD then it is the minRMSD
min_rmsd = interpolated_SIs(Nsort)

! What we call the candidate_rmsd
! and candidate_gradient are the SI
! and energy gradient of the best
! frame (which, again, should be the
! first frame in buffer1)
candidate_rmsd = min_SIs(Nsort)
candidate_gradient = &
        matmul(Ubuffer1(:,:,1),&
        gradientbuffer1(:,:,1))

return

end subroutine checkState_PCM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      getRelativeIndex_PCM FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       INPUT: integer currentVar
!              integer dim(Nvar), var_index
!              integer N
!      IN/OUT: integer dim(Nvar), var_index_diff
!              logical stop_flag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       This is a recursive function that should
!       technically be able search through all
!       neighboring cells with var_indexes that
!       have a Manhattan distance of less than N
!       to the target var_index; but it has only
!       been tested for Nvar = 1 and 2.
!
!       This is (on paper) a little more efficient
!       than just looking at all cells with an
!       index between -N and +N and then filtering
!       out anything with a Manhattan distance
!       larger than the search radius.
!
!       Each new cell it visits it calls
!       getRMSD_1_PCM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine getRelativeIndex_PCM(currentVar,var_index,&
                                      N,var_index_diff,stop_flag)
use ANALYSIS
use PARAMETERS
implicit none

integer,intent(in) :: currentVar
integer,dimension(Nvar),intent(in) :: var_index
integer,intent(in) :: N

integer :: nextVar
integer :: j, k
integer :: population

integer,dimension(Nvar),intent(inout) :: var_index_diff
logical,intent(inout) :: stop_flag

! If we are just starting (indicated by the
! fact that currentVar = 1 so we are on the
! first index):
if (currentVar == 1) then
!   stop_flag = .false.

    ! We have not yet branched out from the
    ! origin yet, so set the distance in each
    ! index to be zero
    var_index_diff = 0

    ! And so the sum of the indexes, k, must
    ! also be zero
    k = 0
else

    ! The sum of the indexes, k, is the
    ! current Manhattan distance of our
    ! branch from the origin; it will let
    ! us know how much more distance we can
    ! distribute for the current index
    k = sum(abs(var_index_diff(1:currentVar-1)))
end if

! If the internal distance tracking has
! reached the maximum Manhattan distance (N):
if (k == N) then

    ! Then we there is no more "distance" to
    ! cover so all the next indexes must be
    ! zero
    do nextVar = currentVar, Nvar
        var_index_diff(nextVar) = 0
    end do

    ! This branch has reached its "end" so we
    ! read the cell corresponding to it
    call getRMSD_1_PCM(var_index + var_index_diff,population)

    ! We have a global variable called
    ! totalnumber_of_frames to keep track
    ! of how many frames have been read
    if (population > 0) then
        Totalnumber_of_frames = &
        Totalnumber_of_frames + population

        ! (Deprecated) The logical stop_flag
        ! indicates that we have found a
        ! non-empty cell so we can technically
        ! exit out early for "accept_first"
        stop_flag = .true.
    end if

! If the current index is at the
! maximum index (which is Nvar) then:
else if (currentVar == Nvar) then

    ! This branch has reached its "end" so we
    ! set the remainder of the Manhattan
    ! distance in this last index
    var_index_diff(currentVar) = N - k

    ! Read the frames in
    call getRMSD_1_PCM(var_index + var_index_diff,population)

    ! We have a global variable called
    ! totalnumber_of_frames to keep track
    ! of how many frames have been read
    if (population > 0) then
        Totalnumber_of_frames = &
        Totalnumber_of_frames + population

        ! (Deprecated) The logical stop_flag
        ! indicates that we have found a
        ! non-empty cell so we can technically
        ! exit out early for "accept_first"
        stop_flag = .true.
    end if

    ! Because this is a nonzero number
    ! (otherwise we would be in the previous
    !  if statement) we must also consider
    ! the negative value as well
    var_index_diff(currentVar) = -N + k

    ! Read the frames in
    call getRMSD_1_PCM(var_index + var_index_diff,population)

    ! We have a global variable called
    ! totalnumber_of_frames to keep track
    ! of how many frames have been read
    if (population > 0) then
        Totalnumber_of_frames = &
        Totalnumber_of_frames + population

        ! (Deprecated) The logical stop_flag
        ! indicates that we have found a
        ! non-empty cell so we can technically
        ! exit out early for "accept_first"
        stop_flag = .true.
    end if

! If neither of the previous cases
! occurred, then we have not reached
! the "end" of a branch
else
    ! We must consider each branch out which
    ! corresponds to a different distance
    ! from the origin for the current index
    do j = -N + k, N - k

        ! This next branch will have a distance
        ! of j from the origin in this index
        var_index_diff(currentVar) = j

        ! And then we recursively call this
        ! same subroutine but now on the next
        ! index
        nextVar = currentVar + 1
        call getRelativeIndex_PCM(nextVar,var_index,&
                              N,var_index_diff,stop_flag)
    end do
end if

end subroutine getRelativeIndex_PCM




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      getInterpolatedGradient FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      OUTPUT: real(dp) dim(3,Natoms), weighted_coords
!              real(dp) dim(3,Natoms), weighted_gradient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      This can be called after reading in all
!      potential interpolating frames and preparing
!      the upper half of the matrix inputCLS; it
!      outputs the interpolated frame and
!      interpolated energy gradient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getInterpolatedGradient(weighted_coords,weighted_gradient)
use ls_rmsd_original
use FUNCTIONS
use ANALYSIS
use PARAMETERS
use SIMILARITY
implicit none

!Information of the frame is held here
real(dp), dimension(3,Natoms), intent(out) :: weighted_coords
real(dp), dimension(3,Natoms), intent(out) :: weighted_gradient
real(dp) :: current_rmsd

real(dp), allocatable :: outputCLS(:)
real(dp), allocatable :: restraints(:,:)

real(dp), allocatable :: frame_weights(:)
real(dp), allocatable :: restraint_values(:)

real(dp),allocatable :: minimized_differences2(:,:)

real(dp),allocatable :: temp_valsbuffer1(:,:)
real(dp),allocatable :: temp_coordsbuffer1(:,:,:)
real(dp),allocatable :: temp_gradientbuffer1(:,:,:)
real(dp),allocatable :: temp_Ubuffer1(:,:,:)
real(dp),allocatable :: temp_inputCLS(:,:)

integer,allocatable :: RMSD_indexes(:)

integer :: Ninterpolation_true
real(dp) :: weight
real(dp) :: total_weight

!Incremental integers
integer :: i, j, k

! We have a user-defined cap to the number of
! frames we can interpolate with; change
! Ninterpolation to go below this cap
Ninterpolation = min(Ninterpolation,&
        Ninterpolation_cap)

! Allocate arrays accordingly
allocate(frame_weights(Ninterpolation),&
         outputCLS(Ncoords+Ninterpolation),&
         restraints(1,Ninterpolation),&
         restraint_values(1),&
         minimized_differences2(Ninterpolation,1))

! Now let's start filling out the bottom half
! of the matrix inputCLS
do i = 1, Ninterpolation

    ! We want all non-diagonal terms to be
    ! zero
    inputCLS(Ncoords+i,:) = 0.0d0

    ! And all diagonal terms to be of size
    ! hyperparameter alpha * RMSD^2
    ! (or whatever SI we are using)
    inputCLS(Ncoords+i,i) = alpha_ratio * &
        SIbuffer1(Nsort,i)**2
end do

! In the default case we are using only one
! restraint; here, all weights must sum to one
! so set the constraint to be a row of 1s
restraints = 1.0d0

! And when we multiply this with the weights
! we should get a sum equalling one
restraint_values = 1.0d0

! The target vector for the top half of the
! matrix is just the zero vector
outputCLS(1:Ncoords) = 0.0d0

! And same for the bottom half of the matrix
outputCLS(Ncoords+1:Ncoords+Ninterpolation) = 0.0d0

! We can now call the constrained least
! squares algorithm on this input, output,
! and constraint matrices to retrive a
! set of weights
call CLS2(inputCLS(1:Ncoords+Ninterpolation,&
                   1:Ninterpolation),&
          Ncoords+Ninterpolation,&
          Ninterpolation,&
          restraints,1,restraint_values,&
          outputCLS,frame_weights)

! If there is a numerical instability (such
! as the matrix being close to linearly
! dependent) then the weights are output as
! being all zero
if (all(frame_weights == 0.0d0)) then

    ! In this, pretend we are in the "accept
    ! best" regime and only weight the best
    ! frame for the final interpolation
    frame_weights(minloc(outputCLS(1:Ninterpolation))) = 1.0d0
end if

! Initialize the interpolations to be zero
weighted_coords = 0.0d0
weighted_gradient = 0.0d0

! Keep track of the total weight for
! debugging purposes
total_weight = 0.0d0

! Also keep track of the largest weighted
! SIs for interpolation error information
largest_SIs = 0.0d0
largest_weighted_SIs = 0.0d0
largest_weighted_SIs2 = 0.0d0

! And keep track of the minimzed differences
! (the sum of the R1 and R2 errors) for
! additional interpolation error info
!minimized_differences2 = &
!        matmul(inputCLS(Ncoords+1:&
!                 Ncoords+Ninterpolation,1:&
!                 Ninterpolation),&
!               reshape(frame_weights,&
!                       (/Ninterpolation,1/)))

! For each frame:
do i = 1, Ninterpolation

    ! Retrieve the weight
    weight = frame_weights(i)
    total_weight = total_weight + weight

    ! And for each SI
    do j = 1, NSIs

        ! Keep track of the largest value,
        largest_SIs(j) = max(largest_SIs(j),&
                             SIbuffer1(j,i))

        ! The largest weighted value
        largest_weighted_SIs(j) = max(&
                largest_weighted_SIs(j),&
                abs(weight)*SIbuffer1(j,i))

        ! And the sum of the weighted
        ! values squared
        largest_weighted_SIs2(j) = &
                largest_weighted_SIs2(j) +&
                (weight*(SIbuffer1(j,i)**2))**2
    end do

    ! Add the coordinates to the interpolation
    ! according to its weight
    weighted_coords = weighted_coords + &
             weight * reshape(inputCLS(1:Ncoords,i),&
                                       (/3,Natoms/))
    
    ! Add the energy gradinet to the
    ! interpolation according to its weight
    weighted_gradient = weighted_gradient + &
             weight * matmul(Ubuffer1(:,:,i),&
                             gradientbuffer1(:,:,i))
end do

! Take the square root of the sum of the
! weighted values squared (alike to RMSD)
largest_weighted_SIs2 = &
        sqrt(largest_weighted_SIs2 / Natoms)

! And finally deallocate!
deallocate(frame_weights,outputCLS,&
           restraints,restraint_values,&
           minimized_differences2)

return

end subroutine getInterpolatedGradient



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      getRMSD_1_PCM FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       INPUT: integer dim(Nvar), var_index
!      OUTPUT: integer population
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Reads through a cell (which has a physical file
!       representation) and deposits the data (like
!       the coords, energy gradient, etc) into the 
!       appropriate buffers (global variables).
!       Returns how many frames it found into the
!       variable population
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getRMSD_1_PCM(var_index,population)
use PARAMETERS
use PHYSICS
use ANALYSIS
use FUNCTIONS
use SIMILARITY
implicit none

!Inputs for file reading
!character(*),intent(in) :: subcell
integer,dimension(Nvar),intent(in) :: var_index
integer :: single_index, index_switch

character(50) :: var_filename
character(150) :: subcell
logical :: subcell_existence
logical :: stop_flag

!Outputs from file reading
integer,intent(out) :: population

!Stores values temporarily
real(dp),dimension(Nvar) :: current_vals
integer :: current_Ntraj
real(dp) :: current_potE
real(dp),dimension(3,Natoms) :: temp_coords
real(dp),dimension(3,Natoms) :: current_coords,current_gradient
logical :: rejectNtraj_flag

!Outputs from the similarity identifier
real(dp),dimension(3,Natoms) :: new_coords
real(dp),dimension(3,3) :: new_U
real(dp),dimension(NSIs) :: new_SI

!Incremental integers and iostate checking
integer :: n,i,j,k,iostate
integer :: endpoint

! If we are in the correct order (in the default case
! any parent cell) then we may have an active memory
! buffer for this data
if (Norder == Norder_order(1)) then

    ! To access the correct index in the memory buffer
    ! we need to convert the Nvar-dimensional index
    ! called var_index into the one-dimensional index
    ! called single_index
    call getFlattened(Nvar,&
            (var_index-previous_var_index),&
            single_index)

    ! And we can access the population from this
    population = populationbuffer2(single_index)

! Otherwise, there should be no index in the memory
! buffer that corresponds to this cell
else

    ! To signify this, we set the single_index to
    ! be out-of-range
    single_index = single_index_max + 1

    ! We set the population to be -1 to signify that
    ! we have not yet searched it
    population = -1
end if

! If the population is -1, then that means we have not
! yet searched the cell according to the memory buffer
if (population < 0) then

    ! Given the var_index (integers describing the cell index)
    ! and the order (parent, cell, etc) we can construct the
    ! filename of the cell with this special format
    write(var_filename,FMT=var_multipleFMT&
          (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
            var_index * multiplier(:,Norder+1)
    
    ! Construct the subcell filename
    subcell = gridpath3//trim(var_filename)
    
    ! See whether this subcell exists
    inquire(file=trim(subcell),exist=subcell_existence)
    
    ! If not, then there is population zero
    if (.not. subcell_existence) then
        population = 0
        populationbuffer2(single_index) = 0
        return

    ! If so, count this to the number of non-empty cells
    ! we have searched, number_of_cells
    else
        number_of_cells = number_of_cells + 1
    end if

    !Open the file corresponding to the cell
    if (unreadable_flag) then
        open(var_filechannel,action="read",form="unformatted",&
             file=trim(subcell))
    else
        open(var_filechannel,action="read",&
             file=trim(subcell))
    end if
    
    ! Initialize the population
    population = 0
    stop_flag = .true.
    do
    
        !Read the candidate frame
        if (unreadable_flag) then
    
            !In unformatted files, the first line are the variables
            !which do not need to be stored
            read(var_filechannel,iostat=iostate) &
                    (current_vals(i),i=1,Nvar)
 
            !If there are no more lines, stop; the population of the cell should be
            !one less the number of times that this portion of the loop was called
            if (iostate /= 0) exit

            population = population + 1

            valsbuffer2(:,population,single_index) = current_vals

            !The next line is trajectory number
            read(var_filechannel) &
                   Ntrajbuffer2(population,single_index)

            !The next line is the potential energy
            read(var_filechannel) &
                   potEbuffer2(population,single_index)

            !The next line is the coordinates
            read(var_filechannel) &
                   ((coordsbuffer2(i,j,population,single_index),i=1,3),j=1,Natoms)
    
            !In unformatted files, the last line is the gradient
            read(var_filechannel) &
                    ((gradientbuffer2(i,j,population,single_index),i=1,3),j=1,Natoms)
        else
    
            !In formatted files, everything (variables, coordinates, and gradient)
            !are stored in one line; FMT1 reads the variables
            read(var_filechannel,FMT=FMT1,advance="no",iostat=iostate) &
                    (current_vals(i),i=1,Nvar)

            !If there are no more lines, stop; the population of the cell should be
            !one less the number of times that this portion of the loop was called
            if (iostate /= 0) exit

            population = population + 1

            valsbuffer2(:,population,single_index) = current_vals

            read(var_filechannel,FMT="(I5)",advance="no") Ntrajbuffer2(population,single_index)
            read(var_filechannel,FMT="(ES10.2)",advance="no") potEbuffer2(population,single_index)
    
            !In formatted files, FMT3 reads the coordinates
            read(var_filechannel,FMT=FMT3,advance="no") &
                   ((coordsbuffer2(i,j,population,single_index),i=1,3),j=1,Natoms)
    
            !In formatted files, FMT3 reads the gradient as well
            read(var_filechannel,FMT=FMT3) &
                    ((gradientbuffer2(i,j,population,single_index),i=1,3),j=1,Natoms)
        end if

        !For diversity, accept first
!       rejectNtraj_flag = .false.
!       do i = 1, Ninterpolation
!           if (Ntrajbuffer2(population,single_index) == Ntrajbuffer1(i)) then
!               rejectNtraj_flag = .true.
!           end if
!       end do
!       if (rejectNtraj_flag) cycle

        ! For permutations, we must check the RMSD of the
        ! target compared to each of the frame's possible
        ! permutations
        do n = 1, Nindistinguishables
    
            ! We store the new order of indexes in the
            ! variable BOND_LABELLING_DATA
            BOND_LABELLING_DATA = INDISTINGUISHABLES(n,:)
        
            ! Permute the new indexes into the variable
            ! temp_coords
            do j = 1, Natoms
                 temp_coords(:,j) = &
                     coordsbuffer2(:,BOND_LABELLING_DATA(j),population,single_index)
            end do

            ! Get the RMSD (and other similarity identifiers)
            call getSIs(temp_coords,new_coords,new_U,new_SI)

            ! If in the "accept best" method and the relevant
            ! SI is within the acceptable thresholds:
            if ((new_SI(Nsort) < outer_threshold_SI).and.&
                (new_SI(Nsort) >= inner_threshold_SI)) then

            ! For diversity, accept best:
            if (diversity_flag) then

                ! Initialize by assuming we are not going
                ! to reject the frame and set index_switch
                ! to 1 to indicate we are looking at the
                ! first frame in the buffer
                rejectNtraj_flag = .false.
                index_switch = 1
                do
                    ! If we have no more frames to compare
                    ! then exit out
                    if (index_switch > Ninterpolation) exit

                    ! If the frame we are looking at has the
                    ! same trajectory number as the one in
                    ! the buffer right now that we want to
                    ! add, then:
                    if (Ntrajbuffer2(population,single_index) &
                            == Ntrajbuffer1(index_switch)) then

                        ! If the SI of interest is larger,
                        ! then we reject the frame we might
                        ! have added
                        if (new_SI(Nsort) >= SIbuffer1(Nsort,index_switch)) then
                            rejectNtraj_flag = .true.

                        ! Otherwise, we add the frame and decrement
                        ! Ninterpolation because we are replacing
                        ! the other frame
                        else
                            Ninterpolation = Ninterpolation - 1
                        end if
        
!                       index_switch = index_switch + 1

                        ! Either way, we can now exit out because
                        ! we are either not adding it or we are
                        ! going to replace the frame we are looking
                        ! at now
                        exit
                    end if

                    ! Increment index_switch to indicate that
                    ! we now look at the next frame in the
                    ! buffer
                    index_switch = index_switch + 1
                end do

                ! If we have rejected a trajectory then there
                ! is no need to replace any data
                if (rejectNtraj_flag) cycle
                
            ! For non-diverse situations:
            else

                ! Always choose to move ALL data over in a
                ! buffer since there is no replacement
                index_switch = Ninterpolation + 1
            end if
!           index_switch = index_switch - 1
        
                ! If we are in "accept first", then we can
                ! stop reading frames after accepting this
                ! first frame
                if (accept_first) iostate = 1
        
                ! For "accept worst" we are going to be
                ! sorting by the WORST SI
                ! Currently deprecated because the diversity
                ! flag does not work with this
                if (accept_worst) then
!                   do i = 1, index_switch+1
                    do i = 1, index_switch
                        if (new_SI(Nsort) > SIbuffer1(Nsort,i)) exit
                    end do

                ! For "accept best" we are going to be
                ! sorting by the BEST SI
                else
!                   do i = 1, index_switch+1
                    do i = 1, index_switch
                        if (new_SI(Nsort) < SIbuffer1(Nsort,i)) exit
                    end do
                end if
        
                ! This frame is considered to be accepted
                ! so increment Ninterpolation
                if (Ninterpolation < Ninterpolation_max) &
                    Ninterpolation = Ninterpolation + 1
        
                ! We shift over everything from i (where we
                ! determined the next best frame is) to
                ! index_switch (where we determined the last
                ! or to-be-replaced frame is)
                call shiftBuffer(i+1,index_switch)
        
                ! Now, index i in the buffer should be empty
                ! so place the new data in that spot

                ! The coords, energy gradient, potential
                ! energy, vals, and trajectory number should
                ! all be in the memory buffer
                
                do j = 1, Natoms
                    coordsbuffer1(:,j,i) = &
                        coordsbuffer2(:,BOND_LABELLING_DATA(j),population,single_index)
                    gradientbuffer1(:,j,i) = &
                        gradientbuffer2(:,BOND_LABELLING_DATA(j),population,single_index)
                end do
        
                potEbuffer1(i) = potEbuffer2(population,single_index)
                valsbuffer1(:,i) = valsbuffer2(:,population,single_index)
                Ntrajbuffer1(i) = Ntrajbuffer2(population,single_index)

                ! The rotation matrix U and the similarity
                ! identifiers SI are not in memory but are
                ! put in the buffer
                Ubuffer1(:,:,i) = new_U
                SIbuffer1(:,i) = new_SI
        
                ! Also put the vectorized coordinates into
                ! the constrained least squares input
                ! matrix called inputCLS, after substracting
                ! the target coordinates
                inputCLS(1:Ncoords,i) =&
                           reshape(new_coords - &
                           coords_target,(/Ncoords/))
            end if
    
        end do

        ! The memory buffer has a maximum capacity such that
        ! afterwards we can only read in frames that are
        ! then not stored in memory
        if (population == buffer2_size) then
            stop_flag = .false.
            exit
        end if
    
    end do

    ! The memory buffer can now be updated to
    ! the population of the cell
    populationbuffer2(single_index) = population

    ! Close the file corresponding to the cell
    ! if we are done, then return
    if (stop_flag) then
        close(var_filechannel)
        return
    end if

! If the population is 0, then that means we have 
! searched the cell according to the memory buffer
! but it is empty. Thus, we can just return
else if (population == 0) then
    return

! If the population is greater than zero, then that
! means we have searched the cell according to the
! and we can read them in
else if (population < buffer2_size) then

    ! Iterate over each frame
    do k = 1, population

        !For diversity, accept first
!       rejectNtraj_flag = .false.
!       do i = 1, Ninterpolation
!           if (Ntrajbuffer2(k,single_index) == Ntrajbuffer1(i)) then
!               rejectNtraj_flag = .true.
!           end if
!       end do
!       if (rejectNtraj_flag) cycle
    
        ! For permutations, we must check the RMSD of the
        ! target compared to each of the frame's possible
        ! permutations
        do n = 1, Nindistinguishables

            ! We store the new order of indexes in the
            ! variable BOND_LABELLING_DATA
            BOND_LABELLING_DATA = INDISTINGUISHABLES(n,:)
        
            ! Permute the new indexes into the variable
            ! temp_coords
            do j = 1, Natoms
                 temp_coords(:,j) = &
                     coordsbuffer2(:,BOND_LABELLING_DATA(j),k,single_index)
            end do

            ! Get the RMSD (and other similarity identifiers)
            call getSIs(temp_coords,new_coords,new_U,new_SI)

            ! If in the "accept best" method and the relevant
            ! SI is within the acceptable thresholds:
            if ((new_SI(Nsort) < outer_threshold_SI).and.&
                (new_SI(Nsort) >= inner_threshold_SI)) then

            ! For diversity, accept best:
            if (diversity_flag) then

                ! Initialize by assuming we are not going
                ! to reject the frame and set index_switch
                ! to 1 to indicate we are looking at the
                ! first frame in the buffer
                rejectNtraj_flag = .false.
                index_switch = 1

                do
                    ! If we have no more frames to compare
                    ! then exit out
                    if (index_switch > Ninterpolation) exit

                    ! If the frame we are looking at has the
                    ! same trajectory number as the one in
                    ! the buffer right now that we want to
                    ! add, then:
                    if (Ntrajbuffer2(k,single_index) &
                            == Ntrajbuffer1(index_switch)) then

                        ! If the SI of interest is larger,
                        ! then we reject the frame we might
                        ! have added
                        if (new_SI(Nsort) >= SIbuffer1(Nsort,index_switch)) then
                            rejectNtraj_flag = .true.

                        ! Otherwise, we add the frame and decrement
                        ! Ninterpolation because we are replacing
                        ! the other frame
                        else
                            Ninterpolation = Ninterpolation - 1
                        end if
        
!                       index_switch = index_switch + 1

                        ! Either way, we can now exit out because
                        ! we are either not adding it or we are
                        ! going to replace the frame we are looking
                        ! at now
                        exit
                    end if

                    ! Increment index_switch to indicate that
                    ! we now look at the next frame in the
                    ! buffer
                    index_switch = index_switch + 1
                end do

                ! If we have rejected a trajectory then there
                ! is no need to replace any data
                if (rejectNtraj_flag) cycle

            ! For non-diverse situations:
            else

                ! Always choose to move ALL data over in a
                ! buffer since there is no replacement
                index_switch = Ninterpolation + 1
            end if
!           index_switch = index_switch - 1
        
                ! If we are in "accept first", then we can
                ! stop reading frames after accepting this
                ! first frame
                if (accept_first) iostate = 1
        
                ! For "accept worst" we are going to be
                ! sorting by the WORST SI
                ! Currently deprecated because the diversity
                ! flag does not work with this
                if (accept_worst) then
!                   do i = 1, index_switch+1
                    do i = 1, index_switch
                        if (new_SI(Nsort) > SIbuffer1(Nsort,i)) exit
                    end do

                ! For "accept best" we are going to be
                ! sorting by the BEST SI
                else
!                   do i = 1, index_switch+1
                    do i = 1, index_switch
                        if (new_SI(Nsort) < SIbuffer1(Nsort,i)) exit
                    end do
                end if
        
                ! This frame is considered to be accepted
                ! so increment Ninterpolation
                if (Ninterpolation < Ninterpolation_max) &
                    Ninterpolation = Ninterpolation + 1
        
                ! We shift over everything from i (where we
                ! determined the next best frame is) to
                ! index_switch (where we determined the last
                ! or to-be-replaced frame is)
                call shiftBuffer(i+1,index_switch)

                ! Now, index i in the buffer should be empty
                ! so place the new data in that spot

                ! The coords, energy gradient, potential
                ! energy, vals, and trajectory number should
                ! all be in the memory buffer
        
                do j = 1, Natoms
                    coordsbuffer1(:,j,i) = &
                        coordsbuffer2(:,BOND_LABELLING_DATA(j),k,single_index)
                    gradientbuffer1(:,j,i) = &
                        gradientbuffer2(:,BOND_LABELLING_DATA(j),k,single_index)
                end do

                potEbuffer1(i) = potEbuffer2(k,single_index)
                valsbuffer1(:,i) = valsbuffer2(:,k,single_index)
                Ntrajbuffer1(i) = Ntrajbuffer2(k,single_index)

                ! The rotation matrix U and the similarity
                ! identifiers SI are not in memory but are
                ! put in the buffer
                Ubuffer1(:,:,i) = new_U
                SIbuffer1(:,i) = new_SI
        
                ! Also put the vectorized coordinates into
                ! the constrained least squares input
                ! matrix called inputCLS, after substracting
                ! the target coordinates
                inputCLS(1:Ncoords,i) =&
                           reshape(new_coords - &
                           coords_target,(/Ncoords/))
            end if
    
        end do
    
    end do

    ! Count this to the number of non-empty cells
    ! we have searched, number_of_cells
    number_of_cells = number_of_cells + 1

! If the population of the memory buffer is equal
! in size to the memory buffer's maximum capacity
! then there may be MORE frames to read in in the
! cell so we should not stop
else

    ! Given the var_index (integers describing the cell index)
    ! and the order (parent, cell, etc) we can construct the
    ! filename of the cell with this special format
    write(var_filename,FMT=var_multipleFMT&
          (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
            var_index * multiplier(:,Norder+1)
    
    ! Construct the subcell filename
    subcell = gridpath3//trim(var_filename)

    ! Count this to the number of non-empty cells
    ! we have searched, number_of_cells
    number_of_cells = number_of_cells + 1

    ! Initialize the population
    population = 0

    ! Open the file corresponding to the cell
    if (unreadable_flag) then
        open(var_filechannel,action="read",form="unformatted",&
             file=trim(subcell))
    else
        open(var_filechannel,action="read",&
             file=trim(subcell))
    end if
end if

! When we have not read all of the frames (because we have
! filled up the memory buffer to the maximum capacity) then
! we read the excess frames in a final do loop

do

    ! Read the candidate frame; these are not stored in
    ! memory so they are stored in a variable called
    ! current_coords, etc.
    if (unreadable_flag) then

        ! In unformatted files, the first line are the
        ! variables which do not need to be stored
        read(var_filechannel,iostat=iostate) &
                (current_vals(i),i=1,Nvar)
 
        ! If there are no more lines, stop; the population
        ! of the cell should be one less the number of times
        ! that this portion of the loop was called
        if (iostate /= 0) exit

        ! For a successful frame reading, increment the
        ! population
        population = population + 1

        ! The next line is the trajectory number
        read(var_filechannel) &
               current_Ntraj

        ! The next line is the potential energy
        read(var_filechannel) &
               current_potE

        ! The next line is the coordinates
        read(var_filechannel) &
               ((current_coords(i,j),i=1,3),j=1,Natoms)

        ! In unformatted files, the last line is the gradient
        read(var_filechannel) &
                ((current_gradient(i,j),i=1,3),j=1,Natoms)
    else

        ! In formatted files, everything (variables,
        ! coordinates, and gradient) are stored in one line;
        ! FMT1 reads the variables
        read(var_filechannel,FMT=FMT1,advance="no",iostat=iostate) &
                (current_vals(i),i=1,Nvar)

        ! If there are no more lines, stop; the population
        ! of the cell should be one less the number of times
        ! that this portion of the loop was called
        if (iostate /= 0) exit

        ! For a successful frame reading, increment the
        ! population
        population = population + 1

        ! In formatted files, the trajectory number has
        ! format I5
        read(var_filechannel,FMT="(I5)",advance="no") current_Ntraj

        ! In formatted files, the potential energy has
        ! format ES10.2
        read(var_filechannel,FMT="(ES10.2)",advance="no") current_potE

        ! In formatted files, FMT3 reads the coordinates
        read(var_filechannel,FMT=FMT3,advance="no") &
               ((current_coords(i,j),i=1,3),j=1,Natoms)

        ! In formatted files, FMT3 reads the energy gradient
        ! as well
        read(var_filechannel,FMT=FMT3) &
                ((current_gradient(i,j),i=1,3),j=1,Natoms)

    end if

    !For diversity, accept first
!   rejectNtraj_flag = .false.
!   do i = 1, Ninterpolation
!       if (current_Ntraj == Ntrajbuffer1(i)) then
!           rejectNtraj_flag = .true.
!       end if
!   end do
!   if (rejectNtraj_flag) cycle

    ! For permutations, we must check the RMSD of the
    ! target compared to each of the frame's possible
    ! permutations
    do n = 1, Nindistinguishables

        ! We store the new order of indexes in the
        ! variable BOND_LABELLING_DATA
        BOND_LABELLING_DATA = INDISTINGUISHABLES(n,:)
    
        ! Permute the new indexes into the variable
        ! temp_coords
        do j = 1, Natoms
             temp_coords(:,j) = &
                 current_coords(:,BOND_LABELLING_DATA(j))
        end do
    
        ! Get the RMSD (and other similarity identifiers)
        call getSIs(temp_coords,new_coords,new_U,new_SI)

        ! If in the "accept best" method and the relevant
        ! SI is within the acceptable thresholds:
        if ((new_SI(Nsort) < outer_threshold_SI).and.&
            (new_SI(Nsort) >= inner_threshold_SI)) then

        ! For diversity, accept best
        if (diversity_flag) then

            ! Initialize by assuming we are not going
            ! to reject the frame and set index_switch
            ! to 1 to indicate we are looking at the
            ! first frame in the buffer
            rejectNtraj_flag = .false.
            index_switch = 1

            do
                ! If we have no more frames to compare
                ! then exit out
                if (index_switch > Ninterpolation) exit

                ! If the frame we are looking at has the
                ! same trajectory number as the one in
                ! the buffer right now that we want to
                ! add, then:
                if (current_Ntraj == Ntrajbuffer1(index_switch)) then

                    ! If the SI of interest is larger,
                    ! then we reject the frame we might
                    ! have added
                    if (new_SI(Nsort) >= SIbuffer1(Nsort,index_switch)) then
                        rejectNtraj_flag = .true.

                    ! Otherwise, we add the frame and decrement
                    ! Ninterpolation because we are replacing
                    ! the other frame
                    else
                        Ninterpolation = Ninterpolation - 1
                    end if
    
!                   index_switch = index_switch + 1

                    ! Either way, we can now exit out because
                    ! we are either not adding it or we are
                    ! going to replace the frame we are looking
                    ! at now
                    exit
                end if

                ! Increment index_switch to indicate that
                ! we now look at the next frame in the
                ! buffer
                index_switch = index_switch + 1
            end do

            ! If we have rejected a trajectory then there
            ! is no need to replace any data
            if (rejectNtraj_flag) cycle

        ! For non-diverse situations:
        else

            ! Always choose to move ALL data over in a
            ! buffer since there is no replacement
            index_switch = Ninterpolation + 1
        end if
!       index_switch = index_switch - 1
    
            ! If we are in "accept first", then we can
            ! stop reading frames after accepting this
            ! first frame
            if (accept_first) iostate = 1
    
            ! For "accept worst" we are going to be
            ! sorting by the WORST SI
            ! Currently deprecated because the diversity
            ! flag does not work with this
            if (accept_worst) then
!               do i = 1, index_switch+1
                do i = 1, index_switch
                    if (new_SI(Nsort) > SIbuffer1(Nsort,i)) exit
                end do

            ! For "accept best" we are going to be
            ! sorting by the BEST SI
            else
!               do i = 1, index_switch+1
                do i = 1, index_switch
                    if (new_SI(Nsort) < SIbuffer1(Nsort,i)) exit
                end do
            end if
    
            ! This frame is considered to be accepted
            ! so increment Ninterpolation
            if (Ninterpolation < Ninterpolation_max) &
                Ninterpolation = Ninterpolation + 1
    
            ! We shift over everything from i (where we
            ! determined the next best frame is) to
            ! index_switch (where we determined the last
            ! or to-be-replaced frame is)
            call shiftBuffer(i+1,index_switch)

            ! Now, index i in the buffer should be empty
            ! so place the new data in that spot

            ! The coords, energy gradient, potential
            ! energy, vals, and trajectory number should
            ! all be in the memory buffer
    
            do j = 1, Natoms
                coordsbuffer1(:,j,i) = &
                    current_coords(:,BOND_LABELLING_DATA(j))
                gradientbuffer1(:,j,i) = &
                    current_gradient(:,BOND_LABELLING_DATA(j))
            end do
    
            potEbuffer1(i) = current_potE
            valsbuffer1(:,i) = current_vals
            Ntrajbuffer1(i) = current_Ntraj

            ! The rotation matrix U and the similarity
            ! identifiers SI are not in memory but are
            ! put in the buffer
            Ubuffer1(:,:,i) = new_U
            SIbuffer1(:,i) = new_SI
    
            ! Also put the vectorized coordinates into
            ! the constrained least squares input
            ! matrix called inputCLS, after substracting
            ! the target coordinates
            inputCLS(1:Ncoords,i) =&
                       reshape(new_coords - &
                       coords_target,(/Ncoords/))
        end if

    end do

end do

! We should be finally done, so close this file
! and return
close(var_filechannel)
return

end subroutine getRMSD_1_PCM



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      getRMSD_2 FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       INPUT: integer dim(Nvar), var_index
!      OUTPUT: integer population
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Reads through a cell (which has a physical file
!       representation) and dumps the data;
!       Returns how many frames it found into the
!       variable population
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getRMSD_2(var_index,population)
use ANALYSIS
use PARAMETERS
implicit none

!Inputs for file reading
!character(*),intent(in) :: subcell
integer,dimension(Nvar),intent(in) :: var_index

character(50) :: var_filename
character(150) :: subcell
logical :: subcell_existence

!Outputs from file reading
integer,intent(out) :: population

real(dp),dimension(Nvar) :: dummy_vals
real(dp),dimension(3,Natoms) :: dummy_coords
integer :: dummy_Ntraj
real(dp) :: dummy_potE

!Incremental integers and iostate checking
integer :: i,j,k,iostate
integer :: endpoint

! Given the var_index (integers describing the cell index)
! and the order (parent, cell, etc) we can construct the
! filename of the cell with this special format
write(var_filename,FMT=var_multipleFMT&
      (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
        var_index * multiplier(:,Norder+1)

! Construct the subcell filename
subcell = gridpath3//trim(var_filename)

! See whether this subcell exists
inquire(file=trim(subcell),exist=subcell_existence)

! If not, then there is population zero
if (.not. subcell_existence) then
    population = 0
    return

! If so, count this to the number of non-empty cells
! we have searched, number_of_cells
else
    number_of_cells = number_of_cells + 1
end if

! Open the file corresponding to the cell; there are different
! options for opening depending on whether the grid is
! specified to be formatted or not
if (unreadable_flag) then
    open(var_filechannel,action="read",form="unformatted",&
         file=trim(subcell))
else
    open(var_filechannel,action="read",&
         file=trim(subcell))
end if

! Initialize the population
population = 1

!Because the buffer may not be large enough to hold all the frames
!and, theoretically, we don't know how many times we need to
!increase the buffer, we need an overarching do loop that is capable
!of increasing the buffer without bounds
do
    !Read the candidate frame
    if (unreadable_flag) then

        !In unformatted files, the first line are the variables
        !which do not need to be stored
        read(var_filechannel,iostat=iostate) &
                (dummy_vals(i),i=1,Nvar)

        !If there are no more lines, stop; the population of the cell should be
        !one less the number of times that this portion of the loop was called
        if (iostate /= 0) then
            population = population - 1
            exit
        end if

        !The next line is the trajectory number
        read(var_filechannel) &
               dummy_Ntraj

        !The next line is the potential energy
        read(var_filechannel) &
               dummy_potE

        !The next line is the coordinates
        read(var_filechannel) &
               ((dummy_coords(i,j),i=1,3),j=1,Natoms)

        !In unformatted files, the last line is the gradient
        read(var_filechannel) &
                ((dummy_coords(i,j),i=1,3),j=1,Natoms)
    else

        !In formatted files, everything (variables, coordinates, and gradient)
        !are stored in one line; FMT1 reads the variables
        read(var_filechannel,FMT=FMT1,advance="no",iostat=iostate) &
                (dummy_vals(i),i=1,Nvar)

        !If there are no more lines, stop; the population of the cell should be
        !one less the number of times that this portion of the loop was called
        if (iostate /= 0) then
            population = population - 1
            exit
        end if

        read(var_filechannel,FMT="(I5)",advance="no") dummy_Ntraj
        read(var_filechannel,FMT="(ES10.2)",advance="no") dummy_potE

        !In formatted files, FMT3 reads the coordinates
        read(var_filechannel,FMT=FMT3,advance="no") &
               ((dummy_coords(i,j),i=1,3),j=1,Natoms)

        !In formatted files, FMT3 reads the gradient as well
        read(var_filechannel,FMT=FMT3) &
                ((dummy_coords(i,j),i=1,3),j=1,Natoms)
    end if

    !Increment the number of frames visited
    population = population + 1
end do
close(var_filechannel)

end subroutine getRMSD_2
!end subroutine getRMSD_1_PCM
!
!subroutine getRMSD_2(var_index,population)
!use ls_rmsd_original
!use ANALYSIS
!use PARAMETERS
!implicit none
!
!!Inputs for file reading
!!character(*),intent(in) :: subcell
!integer,dimension(Nvar),intent(in) :: var_index
!
!character(50) :: var_filename
!character(150) :: subcell
!logical :: subcell_existence
!
!!Outputs from file reading
!integer,intent(out) :: population
!
!real(dp),dimension(Nvar) :: dummy_vals
!real(dp),dimension(3,Natoms) :: dummy_coords
!
!!Incremental integers and iostate checking
!integer :: i,j,k,iostate
!integer :: endpoint
!
!write(var_filename,FMT=var_multipleFMT&
!      (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
!        var_index * multiplier(:,Norder+1)
!
!!Construct the subcell filename
!subcell = gridpath3//trim(var_filename)
!
!!See whether this subcell exists
!inquire(file=trim(subcell),exist=subcell_existence)
!
!if (.not. subcell_existence) then
!        population = 0
!        return
!else
!    number_of_cells = number_of_cells + 1
!end if
!
!!Open the file corresponding to the cell
!if (unreadable_flag) then
!        open(var_filechannel,action="read",form="unformatted",&
!             file=trim(subcell))
!else
!        open(var_filechannel,action="read",&
!             file=trim(subcell))
!end if
!
!!Initialize a variable
!population = 1
!
!!Because the buffer may not be large enough to hold all the frames
!!and, theoretically, we don't know how many times we need to
!!increase the buffer, we need an overarching do loop that is capable
!!of increasing the buffer without bounds
!do
!        !Read the candidate frame
!        if (unreadable_flag) then
!
!                !In unformatted files, the first line are the variables
!                !which do not need to be stored
!                read(var_filechannel,iostat=iostate) &
!                        (dummy_vals(i),i=1,Nvar)
!
!                !If there are no more lines, stop; the population of the cell should be
!                !one less the number of times that this portion of the loop was called
!	        if (iostate /= 0) then
!                        population = population - 1
!                        exit
!                end if
!
!                !The next line is the coordinates
!                read(var_filechannel) &
!                       ((dummy_coords(i,j),i=1,3),j=1,Natoms)
!
!                !In unformatted files, the last line is the gradient
!                read(var_filechannel) &
!                        ((dummy_coords(i,j),i=1,3),j=1,Natoms)
!        else
!
!                !In formatted files, everything (variables, coordinates, and gradient)
!                !are stored in one line; FMT1 reads the variables
!                read(var_filechannel,FMT=FMT1,advance="no",iostat=iostate) &
!                        (dummy_vals(i),i=1,Nvar)
!
!                !If there are no more lines, stop; the population of the cell should be
!                !one less the number of times that this portion of the loop was called
!	        if (iostate /= 0) then
!                        population = population - 1
!                        exit
!                end if
!
!                !In formatted files, FMT3 reads the coordinates
!                read(var_filechannel,FMT=FMT7,advance="no") &
!                       ((dummy_coords(i,j),i=1,3),j=1,Natoms)
!
!                !In formatted files, FMT3 reads the gradient as well
!                read(var_filechannel,FMT=FMT3) &
!                        ((dummy_coords(i,j),i=1,3),j=1,Natoms)
!        end if
!
!        !Increment the number of frames visited
!        population = population + 1
!end do
!close(var_filechannel)
!
!end subroutine getRMSD_2





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      addState_new FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       INPUT: real(dp) dim(Nvar), vals
!              real(dp) dim(3,Natoms), coords
!              real(dp) dim(3,Natoms), gradient
!              real(dp) potE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Keeps internal arrays up-to-date on how
!       many new files have been made, what frames
!       to add, what cells to divy-up, all AFTER
!       the cells have been checked by checkState
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine addState_new(vals,coords,gradient,potE)
use PHYSICS
use PARAMETERS
implicit none

!Inputs for the file
real(dp),dimension(Nvar),intent(in) :: vals
real(dp),dimension(3,Natoms),intent(in) :: coords
real(dp),dimension(3,Natoms),intent(in) :: gradient
real(dp),intent(in) :: potE

!Keeps track of how many frames are in a cell
integer :: population

!An index used to uniquely identify the cell
integer,dimension(Nvar,Norder_max+1) :: var_index

!Character strings used to identify the file
character(50) :: var_filename

!Incremental integers
integer :: i,j, k,l

!Retreive the index of each variable with respect to the grid
!and the real number (rounded) that represents that index
do i = 1, Nvar
    !Repeat this for however many orders of cells deep
    !we have been instructed to go
    do j = 1, Norder_max + 1
        var_index(i,j) = int(vals(i) * divisor(i,j))
    end do
end do

do l = 1, Norder_max+1

    Norder = Norder_order(l)

    population = local_frame_count(Norder+1)

    if (population == var_overcrowd(Norder+1)) then

        call frameAddition(vals,coords,gradient,&
                           potE,&
                           var_index(:,Norder+1))

        valsbuffer1(:,population+1) = vals
        coordsbuffer1(:,:,population+1) = coords
        gradientbuffer1(:,:,population+1) = gradient

        headers(Norder+1) = headers(Norder+1) + 1

        if (Norder < Norder_max) then

            call divyUp(population+1)

        end if

        exit

    else if (population < var_overcrowd(Norder+1)) then

        if (population == 0) then

            if ((Norder == 1) .and. &
                (local_frame_count(1) <= &
                 var_overcrowd(1))) cycle

            Nfile = Nfile + 1

        end if

        call frameAddition(vals,coords,gradient,&
                           potE,&
                           var_index(:,Norder+1))

        exit

    else
        cycle
    end if

end do

end subroutine addState_new


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      frameAddition FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       INPUT: real(dp) dim(Nvar), vals
!              real(dp) dim(3,Natoms), coords
!              real(dp) dim(3,Natoms), gradient
!              real(dp) potE
!              integer dim(Nvar), var_index
!              logical nolabel_flag (deprecated)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Write the data given to the cell specified
!       by var_index and Norder (global variable)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine frameAddition(vals,coords,gradient,&
                         potE,&
                         var_index,nolabel_flag)
use PARAMETERS
use PHYSICS
implicit none

!Inputs for file writing
real(dp),dimension(Nvar),intent(in) :: vals
real(dp),dimension(3,Natoms),intent(in) :: coords,gradient
real(dp),intent(in) :: potE
integer,dimension(Nvar),intent(in) :: var_index

!Character strings used to identify the file
character(50) :: var_filename

logical,optional :: nolabel_flag

!Incremental integers
integer :: i,j,k

! Given the var_index (integers describing the cell index)
! and the order (parent, cell, etc) we can construct the
! filename of the cell with this special format
write(var_filename,FMT=var_multipleFMT&
      (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
        var_index * multiplier(:,Norder+1)

! If the files in the grid are unreadable (unformatted) then:
if (unreadable_flag) then

    ! Open them as being unformatted
    open(filechannel1,file=gridpath2//trim(var_filename),position="append",form="unformatted")

    ! If we are "labelling" the frames then we need to
    ! reorder them before writing to the file
    ! (currently deprecated because the theory is not
    !  not sound)
    if ((force_NoLabels).or.(present(nolabel_flag).and.(nolabel_flag))) then

        ! In the current version of the code, we are writing:
        !    1) The vals
        !    2) The trajectory number
        !    3) The potential energy
        !    4) The coordinates
        !    5) The energy gradient
        write(filechannel1) (vals(j),j=1,Nvar)
        write(filechannel1) Ntraj
        write(filechannel1) potE
        write(filechannel1) ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1) ((gradient(i,j),i=1,3),j=1,Natoms)
    else
        write(filechannel1) (vals(j),j=1,Nvar)
        write(filechannel1) Ntraj
        write(filechannel1) potE
        write(filechannel1) ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1) ((gradient(i,j),i=1,3),j=1,Natoms)
    end if

! Otherwise:
else

    ! Open them as being formatted (normal)
    open(filechannel1,file=gridpath2//trim(var_filename),position="append")
    if ((force_NoLabels).or.(present(nolabel_flag).and.(nolabel_flag))) then
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT="(I5)",advance="no") Ntraj
        write(filechannel1,FMT="(ES10.2)",advance="no") potE
        write(filechannel1,FMT=FMT3,advance="no") ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
    else
        write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
        write(filechannel1,FMT="(I5)",advance="no") Ntraj
        write(filechannel1,FMT="(ES10.2)",advance="no") potE
        write(filechannel1,FMT=FMT3,advance="no") ((coords(i,j),i=1,3),j=1,Natoms)
        write(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
    end if
end if
close(filechannel1)

end subroutine frameAddition


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      divyUp FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       INPUT:  integer frames
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Add frames from buffer1 (the current cell)
!       to the next order of cells (usually from
!       parent to child cell)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine divyUp(frames)
use PARAMETERS
use FUNCTIONS
implicit none

!The number of frames that are to be distributed
integer,intent(in) :: frames

!An index used to uniquely identify the cell
integer,dimension(Nvar) :: var_index

logical :: nolabel_flag = .true.

integer :: i, j

! Increase the order
! (ex. Parent = 0, Child = 1)

Norder = Norder + 1

do i = 1, frames

    ! Calculate the new var_index of each
    ! frame; this depends on divisor from
    ! PARAMETERS

    do j = 1, Nvar
        var_index(j) = int(valsbuffer1(j,i) * &
                           divisor(j,Norder+1))
    end do

    ! Add this frame to the cell; this means
    ! its vals, coordinates, gradient,
    ! potential energy, and trajectory number
    ! (the last variable is implicit)

    call frameAddition(valsbuffer1(:,i),&
                       coordsbuffer1(:,:,i),&
                       gradientbuffer1(:,:,i),&
                       potEbuffer1(i),&
                       var_index,nolabel_flag)
end do

! Go back to the original order

Norder = Norder - 1

end subroutine divyUp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      setAllocations FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Allocate buffer1 and buffer2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine setAllocations()
    use PARAMETERS
    use ANALYSIS
    use SIMILARITY
    implicit none

    ! Buffer1 is for frames that could be
    ! used for interpolation; just good to
    ! have in general, sometimes referred
    ! to as the "buffer"

    allocate(valsbuffer1(Nvar,buffer1_size),&
             coordsbuffer1(3,Natoms,buffer1_size),&
             gradientbuffer1(3,Natoms,buffer1_size),&
             potEbuffer1(buffer1_size),&
             Ubuffer1(3,3,buffer1_size),&
             SIbuffer1(NSIs,buffer1_size))

    allocate(Ntrajbuffer1(buffer1_size))
    allocate(inputCLS(Ncoords+buffer1_size,buffer1_size))

    ! Buffer2 is the "memory buffer" or
    ! sometimes referred to as the "memory"
    ! and keeps track of the data in each
    ! cell according to their single_index
    if (memory_flag) then

!       ! (Deprecated)
!       allocate(vals_hash(single_index_max,Nvar))

        allocate(populationbuffer2(&
                    single_index_max+1),&
                 valsbuffer2(Nvar,&
                    buffer2_size,&
                    single_index_max+1),&
                 Ntrajbuffer2(&
                    buffer2_size,&
                    single_index_max+1),&
                 coordsbuffer2(3,Natoms,&
                    buffer2_size,&
                    single_index_max+1),&
                 gradientbuffer2(3,Natoms,&
                    buffer2_size,&
                    single_index_max+1),&
                 potEbuffer2(&
                    buffer2_size,&
                    single_index_max+1))

        allocate(temppopulationbuffer2(&
                    single_index_max+1),&
                 tempvalsbuffer2(Nvar,&
                    buffer2_size,&
                    single_index_max+1),&
                 tempNtrajbuffer2(&
                    buffer2_size,&
                    single_index_max+1),&
                 tempcoordsbuffer2(3,Natoms,&
                    buffer2_size,&
                    single_index_max+1),&
                 tempgradientbuffer2(3,Natoms,&
                    buffer2_size,&
                    single_index_max+1),&
                 temppotEbuffer2(&
                    buffer2_size,&
                    single_index_max+1))
    end if

    return

end subroutine setAllocations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      unsetAllocations FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Deallocate buffer1 and buffer2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine unsetAllocations()
    use PARAMETERS
    use ANALYSIS
    use SIMILARITY
    implicit none

    deallocate(valsbuffer1,&
               coordsbuffer1,&
               gradientbuffer1,&
               potEbuffer1,&
               Ubuffer1,&
               SIbuffer1)

    deallocate(Ntrajbuffer1)

    deallocate(inputCLS)

    if (memory_flag) then
!       deallocate(vals_hash)

        deallocate(populationbuffer2,&
                   valsbuffer2,&
                   Ntrajbuffer2,&
                   coordsbuffer2,&
                   gradientbuffer2,&
                   potEbuffer2)

        deallocate(temppopulationbuffer2,&
                   tempvalsbuffer2,&
                   tempNtrajbuffer2,&
                   tempcoordsbuffer2,&
                   tempgradientbuffer2,&
                   temppotEbuffer2)
    end if

    return

end subroutine unsetAllocations


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      setSubcellSearchMax FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Calculate how many cells are going to be
!      searched for the subcellsearch given in the
!      ANALYSIS file
!
!      Set single_index_max accordingly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine setSubcellSearchMax()
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

integer :: i,j
integer :: single_index
integer,dimension(Nvar) :: var_index

! The variables ssm1 and ssm2 are specified
! in ANALYSIS where they are transfered to
! subcellsearch_max here.

! They are separate because the dimension
! is specified in PARAMETERS (Norder_max)
! and is, by principle, separate from
! ANALYSIS

do i = 1, min(Norder_max+1,ssm_length)
    subcellsearch_max(i) = ssm1(i)
    subcellsearch_max1(i) = ssm1(i)
    subcellsearch_max2(i) = ssm2(i)
end do

if (memory_flag) then

    ! Initialize the memory buffer assuming
    ! only subcellsearch_max1 and the first
    ! order search are used
    j = subcellsearch_max(&
            Norder_order(1)+1)
    single_index_max = 0

    ! The single_index gets larger as we add
    ! more and more dimensions, so we must
    ! iteratively increase the maximum

    do i = 1, Nvar

        ! The largest single_index has to
        ! occur at the extrema (+/- j) so
        ! only check those points

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      shiftBuffer FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       INPUT:  integer first_index
!               integer last_index
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Moves the first index of all buffers in buffer1
!       so to remove only data corresponding to
!       last_index and free data corresponding to
!       first_index-1. Analagous to a "right shift".
!
!       Ex.
!           [ 1, 2, 3, 4, 5, 6, 7, 8, 9 ]
!
!         first_index = 3; last_index = 7:
!           [ 1, 2, 2, 3, 4, 5, 6, 8, 9 ]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine shiftBuffer(first_index,last_index)
    use PARAMETERS
    use ANALYSIS
    use SIMILARITY
    use FUNCTIONS
    implicit none
    integer,intent(in) :: first_index,last_index
    integer :: j

    ! We have to go in reverse order otherwise
    ! we will overwrite
 
    do j = last_index, first_index, -1
        valsbuffer1(:,j) = valsbuffer1(:,j-1)
        Ntrajbuffer1(j) = Ntrajbuffer1(j-1)
        coordsbuffer1(:,:,j) = coordsbuffer1(:,:,j-1)
        gradientbuffer1(:,:,j) = gradientbuffer1(:,:,j-1)
        potEbuffer1(j) = potEbuffer1(j-1)
        Ubuffer1(:,:,j) = Ubuffer1(:,:,j-1)
        SIbuffer1(:,j) = SIbuffer1(:,j-1)
        inputCLS(:,j) = inputCLS(:,j-1)
    end do

    return
        
end subroutine shiftBuffer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      shiftMemory FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       INPUT:  integer dim (Nvar), delta_var_index
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Move the first index of all buffers in buffer2
!       according to how we moved in the grid
!
!       Ex.
!            Nvar = 2;  single_index_max = 5
!
!        populationbuffer2          single_index (example)
!              + - +                      + - +
!              |-1 |                      | 4 |
!          + - + - + - +              + - + - + - +
!          | 8 | 6 | 0 |              | 2 | 1 | 3 |
!          + - + - + - +              + - + - + - +
!              | 2 |                      | 5 |
!              + - +                      + - +
!
!                |    delta_var_index       |
!                V         = (1,0)          V
!
!              + - +                      + - +
!              |-1 |                      | 9 |
!          + - + - + - +              + - + - + - +
!          | 6 | 0 |-1 |              | 1 | 3 | 7 |
!          + - + - + - +              + - + - + - +
!              |-1 |                      | 8 |
!              + - +                      + - +
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine shiftMemory(delta_var_index)
    use PARAMETERS
    use ANALYSIS
    use FUNCTIONS
    implicit none
    integer,dimension(Nvar),intent(in) :: delta_var_index
    integer :: population
    integer :: single_index,xflattened
    integer,dimension(Nvar) :: xexpanded

    ! If we are in the same cell, then the first index
    ! should be the same, so we can just do nothing

    if (all(delta_var_index == 0)) return

    ! Save the memory buffer to a temporary buffer

    temppopulationbuffer2 = populationbuffer2
    tempvalsbuffer2 = valsbuffer2
    tempNtrajbuffer2 = Ntrajbuffer2
    tempcoordsbuffer2 = coordsbuffer2
    tempgradientbuffer2 = gradientbuffer2
    temppotEbuffer2 = potEbuffer2

    ! Any cell that has a population of -1 will be
    ! considered "undiscovered"

    populationbuffer2 = -1

    ! Traverse through each cell by their now
    ! one-dimensional index

    do single_index = 1, single_index_max

        ! Access the cell's population

        population = temppopulationbuffer2(single_index)

        ! If the cell is undiscovered, there is no data to move
        ! (Note: even if the cell is empty (population = 0) we
        !        must move the knowledge that the cell is empty)

        if (population == -1) cycle

        ! Otherwise, we must calculate the NEW one-dimensional
        ! index of the cell. We do this by:

        ! First, find the Nvar-dimensional index of the cell

        call getExpanded(Nvar,single_index,xexpanded)

        ! Second, get the NEW Nvar-dimensional index of the cell
        ! (simply by adding the change in index)

        ! Third, flatten this to get the one-dimensional index

        call getFlattened(Nvar,xexpanded-delta_var_index,xflattened)

        ! If the index is inside the potential search
        ! radius, then we move the data

        if (xflattened <= single_index_max) then
            populationbuffer2(xflattened) = population
            valsbuffer2(:,1:population,xflattened) = &
                tempvalsbuffer2(:,1:population,single_index)
            Ntrajbuffer2(1:population,xflattened) = &
                tempNtrajbuffer2(1:population,single_index)
            coordsbuffer2(:,:,1:population,xflattened) = &
                tempcoordsbuffer2(:,:,1:population,single_index)
            gradientbuffer2(:,:,1:population,xflattened) = &
                tempgradientbuffer2(:,:,1:population,single_index)
            potEbuffer2(1:population,xflattened) = &
                temppotEbuffer2(1:population,single_index)
        end if

    end do

    ! Comment this line out if you do not want
    ! a memory buffer working (the memory buffer
    ! in this case ALWAYS increases performance)

!   populationbuffer2 = -1

    return

end subroutine shiftMemory


end module interactMultipleGrids
