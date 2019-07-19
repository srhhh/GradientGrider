
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
real(dp),allocatable :: Ubuffer1(:,:,:)
real(dp),allocatable :: RMSDbuffer1(:)

real(dp),allocatable :: CMdiffbuffer1(:)
integer,allocatable :: Ntrajbuffer1(:)
integer,allocatable :: varindexbuffer1(:)

!For the memory buffer

integer,dimension(Nvar) :: previous_var_index
integer,allocatable :: vals_hash(:,:)

integer,allocatable :: populationbuffer2(:)
real(dp),allocatable :: valsbuffer2(:,:,:)
real(dp),allocatable :: coordsbuffer2(:,:,:,:)
real(dp),allocatable :: gradientbuffer2(:,:,:,:)

integer,allocatable :: temppopulationbuffer2(:)
real(dp),allocatable :: tempvalsbuffer2(:,:,:)
real(dp),allocatable :: tempcoordsbuffer2(:,:,:,:)
real(dp),allocatable :: tempgradientbuffer2(:,:,:,:)

!Arrays for interpolation

logical, allocatable :: acceptable_frame_mask(:)
real(dp),allocatable :: inputCLS(:,:)

!In the latter case, we need to keep track of how large
!the buffer gets, and increase it if it gets too large

integer :: buffer1_size
integer :: buffer2_size

!It is useful to estimate how many cells we may be
!looking at at one time

!integer,dimension(Norder_max+1) :: subcellsearch_max1 = (/ 0, 3 /)
!integer,dimension(Norder_max+1) :: subcellsearch_max2 = (/ 0, 3 /)

!integer,dimension(Norder_max+1) :: subcellsearch_max = (/ 0, 3 /)
integer,dimension(Norder_max+1) :: local_frame_count

!Assuming that each subcellsearch_max < 5
!This only approximates the number of cells in the outer
!shell of a Nvar-dimensional cube (not diamond, like
!used in the program)
integer,parameter :: number_of_cells_max = &
        2 * ((1 + 2 * (5))**(Nvar) - (1 + 2 * (5 - 1))**(Nvar))
integer,parameter :: number_of_frames_max = &
        number_of_cells_max * var_overcrowd(2)

!Store the index of the best candidate in approximation_index

integer :: number_of_cells
integer,allocatable :: approximation_index(:)

!Variables related to interpolation

integer :: Ninterpolation
real(dp) :: largest_rmsd
real(dp) :: largest_weighted_rmsd
integer  :: largest_rmsd_index
real(dp) :: largest_weighted_rmsd2
real(dp) :: interpolated_CMdiff

real(dp), allocatable :: temp_frame_weights(:)
real(dp), allocatable :: temp_rmsd_weights(:)

!Other global variables to clean things up

integer :: Totalnumber_of_frames
integer :: Norder
integer :: var_filechannel
real(dp),dimension(3,Natoms) :: var_coords
real(dp),dimension(Natoms,Natoms) :: var_CM
real(dp) :: candidate_rmsd
real(dp),dimension(3,Natoms) :: candidate_gradient

integer :: Naccept




contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               checkState
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine bins a frame into a subcell then checks the file associated
!               with that subcell for frames that are close in RMSD to the original frame
!
!               Only the subcell that is most granular and which is not overcrowded
!               is checked; if this subcell is empty then the subcells closest to it
!               in terms of variable1, variable2 are checked
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               vals                            REAL(DP),DIM(Nvar)              The collective variables used to bin a frame
!               coords                          REAL(DP),DIM(3,Natoms)          The coordinates defining the reference frame
!		filechannels			INTEGER,DIM(Ngrid_total)	The filechannels of each grid to write data
!										such as minimum rmsd
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!               min_rmsd                        REAL(DP)                        The minimum rmsd found between the refence frame and
!                                                                               all frames in the subcell associated with it
!               gradient                        REAL(DP),DIM(3,Natoms)          The gradient associated with a frame
!
!               number_of_frames                INTEGER                         The total number of frames checked for the refernce frame
!               order                           INTEGER                         The order of the subcell that was checked
!               neighbor_check                  INTEGER                         May either be 0 or 1 depending on whether the
!                                                                               subcells surrounding the original subcell were checked
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               subcell                         CHARACTER                       The name of the file associated with a subcell
!               var1_filename                   CHARACTER                       The portion of the subcell string that is
!                                                                               associated with variable1
!               var2_filename                   CHARACTER                       The portion of the subcell string that is
!                                                                               associated with variable2
!
!               var                             REAL                            The original value of the variable
!               var_cell                        REAL                            The value of the variable minus all rounding;
!                                                                               can be considered the "round off"
!               var_roundN                      REAL                            The rounded value of the variable, rounded down
!                                                                               to order N (specified in PARAMETERS)
!               var_index                       INTEGER                         The modulo of the variable with respect to
!                                                                               a subcell spacing
!
!		subcell_existence		LOGICAL				A flag indicating whethere a file associated with
!										a subcell exists
!               population                      INTEGER                         The number of frames checked in a cell
!               subcellsearch_max               INTEGER                         The maximum number of cells radii to check
!                                                                               before quitting getNeighbors
!
!		old_min_rmsd			REAL(DP)			Temporary buffer for the value of min_rmsd
!										when comparing the min_rmsd output of two grids
!
!               stop_flag                       LOGICAL                         A flag to indicate when a nonempty cell was found
!                                                                               in the search for neighbors, so the search can stop
!               index1_1                        INTEGER                         The first index to check for variable1 (neighbors)
!               index1_2                        INTEGER                         The second index to check for variable1 (neighbors)
!               index2_1                        INTEGER                         The first index to check for variable2 (neighbors)
!               index2_2                        INTEGER                         The second index to check for variable2 (neighbors)
!
!               U                               REAL(DP),DIM(3,3)               The rotation matrix used to rotate a frame
!                                                                               and minimize its RMSD with respect to the reference frame
!
!               key                             INTEGER                         An index used to uniquely associate files in traversal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!subroutine checkState(vals,coords,gradient,min_rmsd,filechannels,&
!                      number_of_frames,order,neighbor_check)
!use ANALYSIS
!use VARIABLES
!use PARAMETERS
!implicit none
!integer :: i,j,k
!integer,intent(out),optional :: order,number_of_frames,neighbor_check
!integer :: var1_index0,var2_index0,var1_index1,var2_index1
!integer :: index1_1,index1_2,index2_1,index2_2
!integer :: population
!integer :: OMP_GET_THREAD_NUM
!logical :: stop_flag
!real :: var1,var2,var1_cell,var2_cell,var1_round,var2_round
!real :: var1_round0,var2_round0,var1_round1,var2_round1,var1_round2,var2_round2,var1_round3,var2_round3
!real(dp), dimension(Nvar), intent(in) :: vals
!real(dp), dimension(3,Natoms), intent(in) :: coords
!real(dp), dimension(3,Natoms), intent(inout) :: gradient
!real(dp), dimension(3,3) :: U, old_U
!real(dp), intent(inout) :: min_rmsd
!real(dp)  :: old_min_rmsd
!real(dp), dimension(3,Natoms) :: old_gradient
!logical :: subcell_existence
!character(100) :: subcell
!character(FMTlength) ::  var1_filename0, var2_filename0, var1_filename1, var2_filename1
!character(FMTlength) :: var1_filename,var2_filename
!integer, dimension(1+Ngrid_total),intent(in) :: filechannels
!integer :: subcell0search_max,subcell1search_max
!character(Ngrid_text_length) :: Ngrid_text
!character(5) :: variable_length_text
!integer :: key0, key1
!
!!In order to decrease the number of frames checked
!!If we need to look at adjacent cells for a frame
!!(for instance, if the target cell is empty)
!!Then we stop looking after we are subcellsearch_max cells away
!subcell0search_max = 1
!subcell1search_max = 1
!
!!We start off with zero frames having been checked
!number_of_frames = 0
!neighbor_check = 0
!
!!If we do not find good candidates, then we restore our defaults
!!that are stored in old_...
!old_min_rmsd = min_rmsd
!old_gradient = gradient
!old_U = 0.0d0
!stop_flag = .false.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                 SUBCELL TARGETING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!var_cell keeps track of the child subcell portion of the decimal
!!and throws away the parent subcell portion of the decimal
!var1_cell = vals(1)
!var2_cell = vals(2)
!
!!var_index keeps track of the index of the child subcell
!!inside of its parent subcell
!var1_index0 = int(var1_cell * divisor1_0)
!var2_index0 = int(var2_cell * divisor2_0)
!
!!var_round is the variable rounded to the subcell length multiple
!!This may have more or even less digits than required
!var1_round0 = multiplier1_0 * var1_index0
!var2_round0 = multiplier2_0 * var2_index0
!
!!Repeat for the child subcell portion
!var1_cell = var1_cell - var1_round0
!var2_cell = var2_cell - var2_round0
!
!var1_index1 = int(var1_cell * divisor1_1)
!var2_index1 = int(var2_cell * divisor2_1)
!
!var1_round1 = multiplier1_1 * var1_index1
!var2_round1 = multiplier2_1 * var2_index1
!
!var1_round = var1_round0 + var1_round1
!var2_round = var2_round0 + var2_round1
!
!!The child-level subcells will use var_round
!write(var1_filename1,FMT=FMTorder1) var1_round
!write(var2_filename1,FMT=FMTorder1) var2_round
!
!!The parent-level subcells will use var_round0
!write(var1_filename0,FMT=FMTorder0) var1_round0
!write(var2_filename0,FMT=FMTorder0) var2_round0
!
!!Now, we start iterating over the grids
!do Ngrid = 1, Ngrid_total
!
!order = 1
!neighbor_check = 0
!
!write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
!write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
!gridpath2 = gridpath0//Ngrid_text//"/grid/"
!
!!Key creation (for traversal checking)
!if (traversal_flag) then
!        key0 = bounds1*var2_index0 + var1_index0 + 1
!        inquire(file=gridpath2//trim(adjustl(var1_filename0))//"_"//&
!                                trim(adjustl(var2_filename0))//".dat",&
!                exist=subcell_existence)
!        if (subcell_existence) then
!                inquire(file=gridpath2//trim(adjustl(var1_filename1))//"_"//&
!                                        trim(adjustl(var2_filename1))//".dat",&
!                        exist=subcell_existence)
!                if (subcell_existence) then
!                        key1 = modulo(traversal0(Ngrid,key0),key_start) + &
!                               bounds2*var2_index1 + var1_index1 + 1
!                        traversal1(Ngrid,key1) = traversal1(Ngrid,key1) + 1
!                else
!                        traversal0(Ngrid,key0) = traversal0(Ngrid,key0) + 1  
!                end if
!        else
!                header1 = header1 + 1
!                traversal0(Ngrid,key0) = header1*key_start + traversal0(Ngrid,key0) + 1
!        end if
!end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                 ORDER 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!Construct the subcell filename
!subcell = gridpath2//trim(adjustl(var1_filename1))//"_"//trim(adjustl(var2_filename1))
!
!!See whether this subcell exists
!inquire(file=trim(subcell)//".dat",exist=subcell_existence)
!
!!For bug-testing
!!if (.false.) then
!!print *, ""
!!print *, ""
!!print *, "inquiring for: ", subcell
!!print *, "exists? ", subcell_existence
!!print *, ""
!!end if
!
!!If the subcell exists, we will look in here for a closeby frame
!if (subcell_existence) then
!        !getRMSD reads off the coordinates and calculated the RMSD with
!        !ls_rmsd module
!	call getRMSD_dp(filechannels(1),subcell,coords,population,min_rmsd,gradient,U)
!
!	!However many frames are in the subcell we increment to number_of_frames
!	number_of_frames = number_of_frames + population
!
!	if (min_rmsd /= old_min_rmsd) then
!		old_min_rmsd = min_rmsd
!		old_gradient = gradient
!		old_U = U
!	end if
!
!	stop_flag = .true.
!end if
!
!!If there are no frames in this cell, we can look at one of its
!!neighbors (better than having to look at the parent)
!if ((.not.(subcell_existence)).or.(force_Neighbors)) then
!        neighbor_check = 1
!
!        !Integer i keeps track of how far away from the original
!        !subcell we are; we look at cells on the 'diamond' surrounding
!        !the original subcell
!        do i = 1, subcell1search_max
!                neighbor_check = i
!
!                index1_1 = var1_index1 + i
!                index1_2 = var1_index1 - i
!
!                index2_1 = var2_index1 + i
!                index2_2 = var2_index1 - i
!                
!                if (index1_1 < scaling1_0) then
!                
!                        !Make the name of the subcell
!                        var1_round = var1_round0 + index1_1 * multiplier1_1
!                        var2_round = var2_round0 + var2_round1
!                        write(var1_filename,FMT=FMTorder1) var1_round
!                        write(var2_filename,FMT=FMTorder1) var2_round
!                        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
!                
!                        !Check if it exists
!                        inquire(file=trim(subcell)//".dat",exist=subcell_existence)
!                        if (subcell_existence) then
!                                call getRMSD_dp(filechannels(1),subcell,coords,population,min_rmsd,gradient,U)
!                        end if
!                
!                        number_of_frames = number_of_frames + population
!                
!                end if
!                !Rinse and repeat
!                if (index1_2 .ge. 0) then
!                
!                        var1_round = var1_round0 + index1_2 * multiplier1_1
!                        var2_round = var2_round0 + var2_round1
!                        write(var1_filename,FMT=FMTorder1) var1_round
!                        write(var2_filename,FMT=FMTorder1) var2_round
!                        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
!                
!                        inquire(file=trim(subcell)//".dat",exist=subcell_existence)
!                        if (subcell_existence) then
!                                call getRMSD_dp(filechannels(1),subcell,coords,population,min_rmsd,gradient,U)
!                        end if
!                
!                        number_of_frames = number_of_frames + population
!                end if
!                !Rinse and repeat
!                if (index2_1 < scaling2_0) then
!                
!                        var1_round = var1_round0 + var1_round1
!                        var2_round = var2_round0 + index2_1 * multiplier2_1
!                        write(var1_filename,FMT=FMTorder1) var1_round
!                        write(var2_filename,FMT=FMTorder1) var2_round
!                        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
!                
!                        inquire(file=trim(subcell)//".dat",exist=subcell_existence)
!                        if (subcell_existence) then
!                                call getRMSD_dp(filechannels(1),subcell,coords,population,min_rmsd,gradient,U)
!                        end if
!                
!                        number_of_frames = number_of_frames + population
!                end if
!                !Rinse and repeat
!                if (index2_2 .ge. 0) then
!                
!                        var1_round = var1_round0 + var1_round1
!                        var2_round = var2_round0 + index2_2 * multiplier2_1
!                        write(var1_filename,FMT=FMTorder1) var1_round
!                        write(var2_filename,FMT=FMTorder1) var2_round
!                        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
!                
!                        inquire(file=trim(subcell)//".dat",exist=subcell_existence)
!                        if (subcell_existence) then
!                                call getRMSD_dp(filechannels(1),subcell,coords,population,min_rmsd,gradient,U)
!                        end if
!                
!                        number_of_frames = number_of_frames + population
!                end if
!
!                if (min_rmsd /= old_min_rmsd) then
!                        stop_flag = .true.
!			old_min_rmsd = min_rmsd
!			old_gradient = gradient
!			old_U = U
!                end if
!        
!                !And integer j keeps track of where on the circumfrence
!                !Of the diamond surrounding the subcell we are at
!                do j = 1, i-1
!                
!                        !Yes, there are redundancies for j = 0 and j = i
!                        !This is the first thing we can improve upon
!                        !(maybe with an if statement?)
!                        index1_1 = var1_index1 + j
!                        index1_2 = var1_index1 - j
!        
!                        index2_1 = var2_index1 + i - j
!                        index2_2 = var2_index1 - i + j
!        
!                        !This does the nitty-gritty of checking to see
!                        !If any of the indexes are out of bounds and 
!                        !obtaining their rmsds and coordinates
!                        call getNeighbors(scaling1_0,scaling2_0,multiplier1_1,multiplier2_1,FMTorder1,&
!                                          index1_1,index1_2,index2_1,index2_2,&
!                                          var1_round0,var2_round0,filechannels(1),&
!                                          coords,min_rmsd,gradient,U,number_of_frames)
!        
!                        if (min_rmsd /= old_min_rmsd) then
!				old_min_rmsd = min_rmsd
!				old_gradient = gradient
!				old_U = U
!                        end if
!                end do
!        
!                if (stop_flag) exit
!        end do
!end if
!
!if (stop_flag) then
!	write(filechannels(1+Ngrid),FMT=FMT6) old_min_rmsd
!	cycle
!end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                 ORDER 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!Rinse and repeat
!subcell = gridpath2//trim(adjustl(var1_filename0))//"_"//trim(adjustl(var2_filename0))
!
!inquire(file=trim(subcell)//".dat",exist=subcell_existence)
!
!order = 0
!
!!For bug-testing
!!if (.false.) then
!!print *, ""
!!print *, ""
!!print *, "inquiring for: ", trim(subcell)//".dat" 
!!print *, "exists? ", subcell_existence
!!print *, ""
!!end if
!
!if (subcell_existence) then
!	call getRMSD_dp(filechannels(1),subcell,coords,population,min_rmsd,gradient,U)
!
!	number_of_frames = number_of_frames + population
!
!	if (min_rmsd /= old_min_rmsd) then
!		old_min_rmsd = min_rmsd
!		old_gradient = gradient
!		old_U = U
!	end if
!end if
!
!if ((.not.(subcell_existence)).or.(force_Neighbors)) then
!        stop_flag = .false.
!
!        do i = 1, subcell0search_max
!                neighbor_check = i
!
!                index1_1 = var1_index0 + i
!                index1_2 = var1_index0 - i
!
!                index2_1 = var2_index0 + i
!                index2_2 = var2_index0 - i
!                
!                if (index1_1 < bounds1) then
!                
!                        !Make the name of the subcell
!                        var1_round = index1_1 * multiplier1_0
!                        var2_round = var2_round0
!                        write(var1_filename,FMT=FMTorder0) var1_round
!                        write(var2_filename,FMT=FMTorder0) var2_round
!                        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
!                
!                        !Check if it exists
!                        inquire(file=trim(subcell)//".dat",exist=subcell_existence)
!                        if (subcell_existence) then
!                                call getRMSD_dp(filechannels(1),subcell,coords,population,min_rmsd,gradient,U)
!                        end if
!                
!                        number_of_frames = number_of_frames + population
!                end if
!                
!                !Rinse and repeat
!                if (index1_2 .ge. 0) then
!                
!                        var1_round = index1_2 * multiplier1_0
!                        var2_round = var2_round0
!                        write(var1_filename,FMT=FMTorder0) var1_round
!                        write(var2_filename,FMT=FMTorder0) var2_round
!                        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
!                
!                        inquire(file=trim(subcell)//".dat",exist=subcell_existence)
!                        if (subcell_existence) then
!                                call getRMSD_dp(filechannels(1),subcell,coords,population,min_rmsd,gradient,U)
!                        end if
!                
!                        number_of_frames = number_of_frames + population
!                end if
!                !Rinse and repeat
!                if (index2_1 < bounds2) then
!                
!                        var1_round = var1_round0
!                        var2_round = index2_1 * multiplier2_0
!                        write(var1_filename,FMT=FMTorder0) var1_round
!                        write(var2_filename,FMT=FMTorder0) var2_round
!                        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
!                
!                        inquire(file=trim(subcell)//".dat",exist=subcell_existence)
!                        if (subcell_existence) then
!                                call getRMSD_dp(filechannels(1),subcell,coords,population,min_rmsd,gradient,U)
!                        end if
!                
!                        number_of_frames = number_of_frames + population
!                end if
!                
!                !Rinse and repeat
!                if (index2_2 .ge. 0) then
!                
!                        var1_round = var1_round0
!                        var2_round = index2_2 * multiplier2_0
!                        write(var1_filename,FMT=FMTorder0) var1_round
!                        write(var2_filename,FMT=FMTorder0) var2_round
!                        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
!                
!                        inquire(file=trim(subcell)//".dat",exist=subcell_existence)
!                        if (subcell_existence) then
!                                call getRMSD_dp(filechannels(1),subcell,coords,population,min_rmsd,gradient,U)
!                        end if
!                
!                        number_of_frames = number_of_frames + population
!                end if
!
!                if (min_rmsd /= old_min_rmsd) then
!                        stop_flag = .true.
!			old_min_rmsd = min_rmsd
!			old_gradient = gradient
!			old_U = U
!                end if
!
!        do j = 1, i-1
!        
!                index1_1 = var1_index0 + j
!                index1_2 = var1_index0 - j
!
!                index2_1 = var2_index0 + i - j
!                index2_2 = var2_index0 - i + j
!
!                call getNeighbors(bounds1,bounds2,multiplier1_0,multiplier2_0,FMTorder1,&
!                                  index1_1,index1_2,index2_1,index2_2,&
!                                  0.0,0.0,filechannels(1),&
!                                  coords,min_rmsd,gradient,U,number_of_frames)
!
!                if (min_rmsd /= old_min_rmsd) then
!                        stop_flag = .true.
!			old_min_rmsd = min_rmsd
!			old_gradient = gradient
!			old_U = U
!                end if
!
!        end do
!
!        if (stop_flag) exit
!        end do
!
!end if
!
!write(filechannels(1+Ngrid),FMT=FMT6) old_min_rmsd
!end do
!
!!This is to see whether old_U has its default value (all zeroes) or not
!if (any(abs(old_U) > 1.0d-20)) gradient = matmul(old_U,old_gradient)
!min_rmsd = old_min_rmsd
!
!end subroutine checkState
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!       SUBROUTINE
!!               getNeighbors
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!       PURPOSE
!!               This subroutine checks the four subcells in the cartesian product
!!               of relative indexes (index1_1,index1_2) x (index2_1, index2_2), with knowledge
!!               of the scaling, the variables rounded down, and the counter
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!       INPUT                           KIND                            DESCRIPTION
!!
!!               scaling1                        INTEGER                         The ratio between the sizes of parent vs child
!!                                                                               subcells with respect to variable1
!!               scaling2                        INTEGER                         The ratio between the sizes of parent vs child
!!                                                                               subcells with respect to variable2
!!               multiplier1                     REAL                            The actual size of the child subcell
!!                                                                               with respect to variable1
!!               multiplier2                     REAL                            The actual size of the child subcell
!!                                                                               with respect to variable2
!!               FMTorder                        CHRACTER(6)                     The formatting for filenames associated with the
!!                                                                               child subcell
!!               var_round                       REAL                            The variables associated with the reference frame
!!                                                                               after rounding to order N
!!
!!               index1_1                        INTEGER                         The first index to check for variable1
!!               index1_2                        INTEGER                         The second index to check for variable1
!!               index2_1                        INTEGER                         The first index to check for variable2
!!               index2_2                        INTEGER                         The second index to check for variable2
!!
!!               coords                          REAL(DP),DIM(3,Natoms)          The coordinates defining the reference frame
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!       OUTPUT                          KIND                            DESCRIPTION
!!
!!               min_rmsd                        REAL(DP)                        The minimum rmsd found between the reference frame and
!!                                                                               all frames in the subcell associated with it
!!               gradient                        REAL(DP),DIM(3,Natoms)          The gradient associated with a frame
!!               U                               REAL(DP),DIM(3,3)               The rotation matrix used to rotate a frame
!!                                                                               and minimize its RMSD with respect to the reference frame
!!
!!               number_of_frames                INTEGER                         The total number of frames checked for the reference frame
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!!
!!               subcell                         CHARACTER                       The name of the file associated with a subcell
!!               var1_filename                   CHARACTER                       The portion of the subcell string that is
!!                                                                               associated with variable1
!!               var2_filename                   CHARACTER                       The portion of the subcell string that is
!!                                                                               associated with variable2
!!
!!               var_round                       REAL                            The variable associated with the subcell in question,
!!                                                                               rounded down to the order specified by multiplier
!!
!!               flag1                           LOGICAL                         Is relative index (index1_1) not out-of-bounds?
!!               flag2                           LOGICAL                         Is relative index (index1_2) not out-of-bounds?
!!               flag3                           LOGICAL                         Is relative index (index2_1) not out-of-bounds?
!!               flag4                           LOGICAL                         Is relative index (index2_2) not out-of-bounds?
!!
!!               population                      INTEGER                         The number of frames in a subcell
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!       FILES                             FILETYPE                      DESCRIPTION
!!
!!               gridpath2//                     DAT                     The file houses all the information (var1,var2,coords,
!!                 trim(subcell).dat                                     gradient) of each frame that is in the cell associated
!!                                                                       with this file
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!subroutine getNeighbors(scaling1,scaling2,multiplier1,multiplier2,FMTorder,&
!                        index1_1,index1_2,index2_1,index2_2,&
!                        var1_round0,var2_round0,filechannel_thread,&
!                        coords,min_rmsd,gradient,U,number_of_frames)
!use ANALYSIS
!use PARAMETERS
!implicit none
!logical :: flag1,flag2,flag3,flag4,subcell_existence
!integer,intent(in) :: scaling1, scaling2
!integer,intent(in) :: index1_1,index1_2,index2_1,index2_2
!real,intent(in) :: multiplier1,multiplier2
!real,intent(in) :: var1_round0,var2_round0
!real :: var1_round,var2_round
!integer,intent(inout) :: number_of_frames
!character(6),intent(in) :: FMTorder
!integer,intent(in) :: filechannel_thread
!character(FMTlength) :: var1_filename,var2_filename
!character(100) :: subcell
!real(dp),dimension(3,3),intent(inout) :: U
!real(dp), intent(inout) :: min_rmsd
!real(dp), dimension(3,Natoms), intent(inout) :: gradient
!integer :: population
!real(dp),intent(in), dimension(3,Natoms) :: coords
!
!!We need to check if any of the indexes (neighbors!) are out of bounds
!flag1 = index1_1 < scaling1
!flag2 = index1_2 .ge. 0
!flag3 = index2_1 < scaling2
!flag4 = index2_2 .ge. 0
!
!!We check indexes in pairs (var1, var2)
!if ((flag1) .and. (flag3)) then
!
!        !Make the name of the subcell
!	var1_round = var1_round0 + index1_1 * multiplier1
!	var2_round = var2_round0 + index2_1 * multiplier2
!        write(var1_filename,FMT=FMTorder) var1_round
!        write(var2_filename,FMT=FMTorder) var2_round
!        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
!
!	!Check if it exists
!	inquire(file=trim(subcell)//".dat",exist=subcell_existence)
!
!	!If it exists, get the rmsd of each frame in it
!	if (subcell_existence) then
!		call getRMSD_dp(filechannel_thread,subcell,coords,population,min_rmsd,gradient,U)
!	end if
!
!	!Increment however many frames we have
!	number_of_frames = number_of_frames + population
!end if
!
!!Rinse and repeat
!if ((flag1) .and. (flag4)) then
!
!	var1_round = var1_round0 + index1_1 * multiplier1
!	var2_round = var2_round0 + index2_2 * multiplier2
!        write(var1_filename,FMT=FMTorder) var1_round
!        write(var2_filename,FMT=FMTorder) var2_round
!        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
!
!	inquire(file=trim(subcell)//".dat",exist=subcell_existence)
!	if (subcell_existence) then
!		call getRMSD_dp(filechannel_thread,subcell,coords,population,min_rmsd,gradient,U)
!	end if
!
!	number_of_frames = number_of_frames + population
!end if
!
!!Rinse and repeat
!if ((flag2) .and. (flag3)) then
!
!	var1_round = var1_round0 + index1_2 * multiplier1
!	var2_round = var2_round0 + index2_1 * multiplier2
!        write(var1_filename,FMT=FMTorder) var1_round
!        write(var2_filename,FMT=FMTorder) var2_round
!        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
!
!	inquire(file=trim(subcell)//".dat",exist=subcell_existence)
!	if (subcell_existence) then
!		call getRMSD_dp(filechannel_thread,subcell,coords,population,min_rmsd,gradient,U)
!	end if
!
!	number_of_frames = number_of_frames + population
!end if
!
!!Rinse and repeat
!if ((flag2) .and. (flag4)) then
!
!	var1_round = var1_round0 + index1_2 * multiplier1
!	var2_round = var2_round0 + index2_2 * multiplier2
!        write(var1_filename,FMT=FMTorder) var1_round
!        write(var2_filename,FMT=FMTorder) var2_round
!        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
!
!	inquire(file=trim(subcell)//".dat",exist=subcell_existence)
!	if (subcell_existence) then
!		call getRMSD_dp(filechannel_thread,subcell,coords,population,min_rmsd,gradient,U)
!	end if
!
!	number_of_frames = number_of_frames + population
!end if
!
!
!
!end subroutine getNeighbors
!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!               FUNCTION GETRMSD_DP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  actually pretty self-explanatory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!subroutine getRMSD_dp(filechannel_thread,subcell,coords,population,min_rmsd,gradient,U)
!
!use ls_rmsd_original
!use ANALYSIS
!use PARAMETERS
!implicit none
!integer,intent(out) :: population
!real(dp),intent(inout), dimension(3,Natoms) :: gradient
!real(dp),dimension(3,Natoms) :: candidate_gradient
!real(dp),intent(inout) :: min_rmsd
!integer,intent(in) :: filechannel_thread
!integer :: OMP_GET_THREAD_NUM
!integer :: i,j,k,iostate
!character(*),intent(in) :: subcell
!real(dp),intent(in), dimension(3,Natoms) :: coords
!real(dp) :: candidate_min_rmsd
!real(dp), dimension(3) :: x_center, y_center
!real(dp), allocatable :: g(:,:)
!real(dp),dimension(3,3),intent(inout) :: U
!real(dp),dimension(3,3) :: candidate_U
!real(dp), dimension(3,Natoms) :: coords2
!
!!Open the file corresponding to the cell
!if (unreadable_flag) then
!        open(filechannel_thread,action="read",file=trim(subcell)//".dat",form="unformatted")
!else
!        open(filechannel_thread,action="read",file=trim(subcell)//".dat")
!end if
!
!population = 0
!do
!        !Read the candidate frame
!        if (unreadable_flag) then
!                read(filechannel_thread,iostat=iostate)
!	        if (iostate /= 0) exit
!                read(filechannel_thread) ((coords2(i,j),i=1,3),j=1,Natoms)
!        else
!                read(filechannel_thread,FMT=FMT7,advance="no",iostat=iostate) ((coords2(i,j),i=1,3),j=1,Natoms)
!	        if (iostate /= 0) exit
!        end if
!
!        call rmsd_dp(Natoms,coords2,coords,1,candidate_U,x_center,y_center,candidate_min_rmsd)!,.false.,g)
!
!        !Increment the number of frames visited
!        population = population + 1
!
!        !If in the "accept worst" method:
!        if (accept_worst) then
!
!        !If it is too low, reject the candidate frame
!        if (candidate_min_rmsd < min_rmsd) then
!                if (unreadable_flag) then
!                        read(filechannel_thread) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
!                else
!                        read(filechannel_thread,FMT=FMT3) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
!                end if
!
!        !Otherwise if it is low enough, accept the candidate frame (and save its rotation matrix!)
!        else if (candidate_min_rmsd < threshold_rmsd) then
!                min_rmsd = candidate_min_rmsd
!                U = candidate_U
!
!                if (unreadable_flag) then
!                        read(filechannel_thread) ((gradient(j,k),j=1,3),k=1,Natoms)
!                else
!                        read(filechannel_thread,FMT=FMT3) ((gradient(j,k),j=1,3),k=1,Natoms)
!                end if
!
!                !In special cases, we exit immediately afterwards
!                if (accept_first) exit
!
!        !If it is too high, also reject the candidate frame
!        else
!                if (unreadable_flag) then
!                        read(filechannel_thread) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
!                else
!                        read(filechannel_thread,FMT=FMT3) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
!                end if
!        end if
!
!        !If in the "accept best" method:
!        else
!
!        !If it is low enough, accept the candidate frame (and save its rotation matrix!)
!        if (candidate_min_rmsd < min_rmsd) then
!                min_rmsd = candidate_min_rmsd
!                U = candidate_U
!
!                if (unreadable_flag) then
!                        read(filechannel_thread) ((gradient(j,k),j=1,3),k=1,Natoms)
!                else
!                        read(filechannel_thread,FMT=FMT3) ((gradient(j,k),j=1,3),k=1,Natoms)
!                end if
!
!                !In special cases, we exit immediately afterwards
!                if (accept_first) exit
!
!        !Otherwise, just don't record the gradient
!        else
!                if (unreadable_flag) then
!                        read(filechannel_thread) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
!                else
!                        read(filechannel_thread,FMT=FMT3) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
!                end if
!        end if
!
!        end if
!
!end do
!close(filechannel_thread)
!
!end subroutine getRMSD_dp
!







subroutine checkState_new(vals,coords,gradient,min_rmsd,filechannels,&
                          number_of_frames,order,neighbor_check)
use ANALYSIS
use VARIABLES
use PARAMETERS
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
real(dp), dimension(3) :: x_center, y_center
real(dp), allocatable :: g(:,:)
real(dp), dimension(3,3) :: new_U

real(dp),dimension(Nvar) :: var_cell
integer,dimension(Nvar,Norder_max+1) :: var_index
integer,dimension(Nvar,Norder_max+1) :: var_index_diff

integer  :: largest_rmsd_error_index
real(dp) :: largest_rmsd_error

var_filechannel = filechannels(1)
var_coords = coords

!We start off with zero frames having been checked
local_frame_count = 0
Totalnumber_of_frames = 0
order = 0
number_of_cells = 0
neighbor_check = 0
approximation_index = 0
Ninterpolation = 0
stop_flag = .false.

acceptable_frame_mask = .false.

candidate_rmsd = min_rmsd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 SUBCELL TARGETING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Retrieve the index of each variable with
!respect to the grid and the real number
!(rounded) that represents that index
do i = 1, Nvar
    !Repeat this for however many orders of cells deep
    !we have been instructed to go
    do j = 1, Norder_max + 1
        var_index(i,j) = int(vals(i) * divisor(i,j))
    end do
end do

!Now, we start iterating over the grids
do k = 1, Ngrid_total
    
    !Streamline the process by storing the
    !path to the grid in a string gridpath3
    write(variable_length_text,FMT=FMT5_variable)&
            Ngrid_text_length
    write(Ngrid_text,FMT="(I0."//&
            trim(adjustl(variable_length_text))//&
            ")") k
    gridpath3 = gridpath0//Ngrid_text//"/grid/"
    
    !The way we check, we check by order first
    do l = 1, Norder_max+1
    
        !The user can specify in what order to
        !check the orders through this array
        Norder = Norder_order(l)
        
        !Read the frames in the cells and
        !process their RMSDs

        !If 
        if (k == grid_addition) then
            if ((l==1) .or. &
                (Totalnumber_of_frames == 0)) then
                call getRMSD_1(&
                        var_index(:,Norder+1),&
                        population)
            else
                call getRMSD_2(&
                        var_index(:,Norder+1),&
                        population)
            end if
            local_frame_count(Norder+1) = &
                    population
        else
            if ((l==1) .or. &
                (Totalnumber_of_frames == 0)) then
                call getRMSD_1(&
                        var_index(:,Norder+1),&
                        population)
            else
            end if
        end if
        
!        !If the cell is populated...
!        if (population > 0) then
!            if ((l==1).or.(Totalnumber_of_frames == 0)) then
!        
!            !However many frames are in the subcell we
!            !increment to number_of_frames
!            Totalnumber_of_frames = &
!            Totalnumber_of_frames + population
!            
!!           if ((population > var_overcrowd(Norder+1)) &
!!               .and. (Norder < Norder_max)) then
!            if ((population >= var_overcrowd(Norder+1)) &
!                .and. (Norder < Norder_max)) then
!                print *, "hello"
!                cycle
!            end if
!
!            !This is good enough, no need to look at more frames
!            stop_flag = .true.
!    
!            end if
!        end if

        Totalnumber_of_frames = &
            Totalnumber_of_frames + population

        if (population >= var_overcrowd(Norder+1)) then
            if (Norder < Norder_max) cycle
        end if
        
        !If the cell is unpopulated or a certain flag is true,
        !then we go ahead and look at the neighbors of this cell
        if ((force_Neighbors) .or. (population == 0)) then
        
            !Integer i keeps track of how far away from the original
            !subcell we are; we look at cells on the 'diamond' surrounding
            !the original subcell
            do i = 1, subcellsearch_max(Norder+1)
        
                var_index_diff = 0
        
                call getRelativeIndex(1,var_index(:,Norder+1),i,&
                                      var_index_diff,stop_flag)
                
                if ((stop_flag) .and. (.not. force_Neighbors)) exit
            end do
        end if
        
        !Psyche!
!       if (stop_flag) exit
        if (population > 0) exit
    
    end do
    
    order = order + Norder
    
    if (testtraj_flag) write(filechannels(1+k),FMT=FMT6) candidate_rmsd
    
end do

number_of_frames = Totalnumber_of_frames
neighbor_check = number_of_cells
Norder_total(1+order/Ngrid_total) = Norder_total(1+order/Ngrid_total) + 1

if ((number_of_cells == 0)) then
    return
else
    chosen_index = maxval(approximation_index(1:number_of_cells),DIM=1)
end if

if (chosen_index > 0) then
    if ((interpolation_flag).and.(Ninterpolation > 1)) then
        call getInterpolatedGradient(new_coords,gradient)

        !Remark: Ken Dill's code uses a version of the RMSD
        !        that divides by N, not N - 1

        min_rmsd = sqrt(sum((new_coords)**2)/(Natoms))
    else
        gradient = matmul(Ubuffer1(:,:,chosen_index),&
                          gradientbuffer1(:,:,chosen_index))

        min_rmsd = RMSDbuffer1(chosen_index)

        largest_weighted_rmsd = min_rmsd
        largest_weighted_rmsd2 = min_rmsd**2
    end if

    candidate_rmsd = RMSDbuffer1(chosen_index)
    candidate_gradient = matmul(Ubuffer1(:,:,chosen_index),&
            gradientbuffer1(:,:,chosen_index))
else
    min_rmsd = default_rmsd
end if

return

end subroutine checkState_new

subroutine checkState_new_permute(vals,coords,gradient,&
                min_rmsd,filechannels,&
                number_of_frames,order,neighbor_check)
use ANALYSIS
use VARIABLES
use PARAMETERS
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
real(dp), dimension(3) :: x_center, y_center
real(dp), allocatable :: g(:,:)
real(dp), dimension(3,3) :: new_U

real(dp),dimension(Nvar) :: var_cell
integer,dimension(Nvar,Norder_max+1) :: var_index
integer,dimension(Nvar,Norder_max+1) :: var_index_diff

integer  :: largest_rmsd_error_index
real(dp) :: largest_rmsd_error

var_filechannel = filechannels(1)
var_coords = coords

!We start off with zero frames having been checked
local_frame_count = 0
Totalnumber_of_frames = 0
order = 0
number_of_cells = 0
neighbor_check = 0
approximation_index = 0
Ninterpolation = 0
stop_flag = .false.

acceptable_frame_mask = .false.

candidate_rmsd = min_rmsd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 SUBCELL TARGETING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Retrieve the index of each variable with
!respect to the grid and the real number
!(rounded) that represents that index
do i = 1, Nvar
    !Repeat this for however many orders of cells deep
    !we have been instructed to go
    do j = 1, Norder_max + 1
        var_index(i,j) = int(vals(i) * divisor(i,j))
    end do
end do

!Now, we start iterating over the grids
do k = 1, Ngrid_total
    
    !Streamline the process by storing the
    !path to the grid in a string gridpath3
    write(variable_length_text,FMT=FMT5_variable)&
            Ngrid_text_length
    write(Ngrid_text,FMT="(I0."//&
            trim(adjustl(variable_length_text))//&
            ")") k
    gridpath3 = gridpath0//Ngrid_text//"/grid/"
    
    !The way we check, we check by order first
    do l = 1, Norder_max+1
    
        !The user can specify in what order to
        !check the orders through this array
        Norder = Norder_order(l)
        
        !Read the frames in the cells and
        !process their RMSDs

        !If 
        if (k == grid_addition) then
            if ((l==1) .or. &
                (Totalnumber_of_frames == 0)) then
                call getRMSD_1_permute(&
                        var_index(:,Norder+1),&
                        population)
            else
                call getRMSD_2(&
                        var_index(:,Norder+1),&
                        population)
            end if
            local_frame_count(Norder+1) = &
                    population
        else
            if ((l==1) .or. &
                (Totalnumber_of_frames == 0)) then
                call getRMSD_1_permute(&
                        var_index(:,Norder+1),&
                        population)
            else
            end if
        end if
        
!       !If the cell is populated...
!       if (population > 0) then
!           if ((l==1).or.(Totalnumber_of_frames == 0)) then
!       
!           !However many frames are in the subcell we
!           !increment to number_of_frames
!           Totalnumber_of_frames = &
!           Totalnumber_of_frames + population
!           
!           !This is good enough, no need to look at more frames
!           stop_flag = .true.
!   
!           end if
!       end if

        if (population >= var_overcrowd(Norder+1)) then
            if (Norder < Norder_max) cycle
        else
            Totalnumber_of_frames = &
                Totalnumber_of_frames + population
        end if
        
        !If the cell is unpopulated or a certain flag is true,
        !then we go ahead and look at the neighbors of this cell
        if ((force_Neighbors) .or. (population == 0)) then
        
            !Integer i keeps track of how far away from the original
            !subcell we are; we look at cells on the 'diamond' surrounding
            !the original subcell
            do i = 1, subcellsearch_max(Norder+1)
        
                var_index_diff = 0
        
                call getRelativeIndex(1,var_index(:,Norder+1),i,&
                                      var_index_diff,stop_flag)
                
                if ((stop_flag) .and. (.not. force_Neighbors)) exit
            end do
        end if
        
        !Psyche!
!       if (stop_flag) exit

        !If we found a non-empty cell then our search has
        !been over this particular order and does not need
        !to go over other orders
!       if (population > 0) then
        if (Totalnumber_of_frames > 0) then
                order = order + Norder
                exit
        end if
    
    end do
    
    !Record the RMSD encountered here for later analysis:
    !particularly, percent-RMSD graphs
    if (testtraj_flag) write(filechannels(1+k),FMT=FMT6) candidate_rmsd
    
end do

number_of_frames = Totalnumber_of_frames
neighbor_check = number_of_cells
Norder_total(1+order/Ngrid_total) = Norder_total(1+order/Ngrid_total) + 1

if ((number_of_cells == 0)) then
    return
else
    chosen_index = maxval(approximation_index(1:number_of_cells),DIM=1)
end if

if (chosen_index > 0) then
    if ((interpolation_flag).and.(Ninterpolation > 1)) then
        call getInterpolatedGradient(new_coords,gradient)

        !Remark: Ken Dill's code uses a version of the RMSD
        !        that divides by N, not N - 1

        min_rmsd = sqrt(sum((new_coords)**2)/(Natoms))
    else
        gradient = matmul(Ubuffer1(:,:,chosen_index),&
                          gradientbuffer1(:,:,chosen_index))

        min_rmsd = RMSDbuffer1(chosen_index)

        largest_weighted_rmsd = min_rmsd
        largest_weighted_rmsd2 = min_rmsd**2
    end if

    candidate_rmsd = RMSDbuffer1(chosen_index)
    candidate_gradient = matmul(Ubuffer1(:,:,chosen_index),&
            gradientbuffer1(:,:,chosen_index))

else
    min_rmsd = default_rmsd
end if

return

end subroutine checkState_new_permute

subroutine checkState_new_permute_cap(vals,coords,gradient,&
                min_rmsd,filechannels,&
                number_of_frames,order,neighbor_check)
use ls_rmsd_original
use ANALYSIS
use VARIABLES
use PARAMETERS
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
real(dp), dimension(3) :: x_center, y_center
real(dp), allocatable :: g(:,:)
real(dp), dimension(3,3) :: new_U

real(dp),dimension(Nvar) :: var_cell
integer,dimension(Nvar,Norder_max+1) :: var_index
integer,dimension(Nvar,Norder_max+1) :: var_index_diff

integer  :: largest_rmsd_error_index
real(dp) :: largest_rmsd_error
real(dp),dimension(Natoms,Natoms) :: temp_CM

var_filechannel = filechannels(1)
var_coords = coords

call getCoulombMatrix(Natoms,var_coords,charges,var_CM)

!We start off with zero frames having been checked
local_frame_count = 0
Totalnumber_of_frames = 0
order = 0
number_of_cells = 0
neighbor_check = 0
Ninterpolation = 0
stop_flag = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 SUBCELL TARGETING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Retrieve the index of each variable with
!respect to the grid and the real number
!(rounded) that represents that index
do i = 1, Nvar
    !Repeat this for however many orders of cells deep
    !we have been instructed to go
    do j = 1, Norder_max + 1
        var_index(i,j) = int(vals(i) * divisor(i,j))
    end do
end do

!Now, we start iterating over the grids
do k = 1, Ngrid_total
    
    !Streamline the process by storing the
    !path to the grid in a string gridpath3
    write(variable_length_text,FMT=FMT5_variable)&
            Ngrid_text_length
    write(Ngrid_text,FMT="(I0."//&
            trim(adjustl(variable_length_text))//&
            ")") k
    gridpath3 = gridpath0//Ngrid_text//"/grid/"
    
    !The way we check, we check by order first
    do l = 1, Norder_max+1
    
        !The user can specify in what order to
        !check the orders through this array
        Norder = Norder_order(l)
        
        !Read the frames in the cells and
        !process their RMSDs

        !If 
        if (k == grid_addition) then
            if ((l==1) .or. &
                (Totalnumber_of_frames == 0)) then
                call getRMSD_1_permute_cap(&
                        var_index(:,Norder+1),&
                        population)
            else
                call getRMSD_2(&
                        var_index(:,Norder+1),&
                        population)
            end if
            local_frame_count(Norder+1) = &
                    population
        else
            if ((l==1) .or. &
                (Totalnumber_of_frames == 0)) then
                call getRMSD_1_permute_cap(&
                        var_index(:,Norder+1),&
                        population)
            else
            end if
        end if
        
!       !If the cell is populated...
!       if (population > 0) then
!           if ((l==1).or.(Totalnumber_of_frames == 0)) then
!       
!           !However many frames are in the subcell we
!           !increment to number_of_frames
!           Totalnumber_of_frames = &
!           Totalnumber_of_frames + population
!           
!           !This is good enough, no need to look at more frames
!           stop_flag = .true.
!   
!           end if
!       end if

        if (population >= var_overcrowd(Norder+1)) then
            if (Norder < Norder_max) cycle
        else
            Totalnumber_of_frames = &
                Totalnumber_of_frames + population
        end if
        
        !If the cell is unpopulated or a certain flag is true,
        !then we go ahead and look at the neighbors of this cell
        if ((force_Neighbors) .or. (population == 0)) then
        
            !Integer i keeps track of how far away from the original
            !subcell we are; we look at cells on the 'diamond' surrounding
            !the original subcell
            do i = 1, subcellsearch_max(Norder+1)
        
                var_index_diff = 0
        
                call getRelativeIndex_permute_cap(1,var_index(:,Norder+1),i,&
                                      var_index_diff,stop_flag)
                
                if ((stop_flag) .and. (.not. force_Neighbors)) exit
            end do
        end if
        

        !If we found a non-empty cell then our search has
        !been over this particular order and does not need
        !to go over other orders
!       if (population > 0) then
        if (Totalnumber_of_frames > 0) then
                order = order + Norder
                exit
        end if
    
    end do

    !Record the RMSD encountered here for later analysis:
    !particularly, percent-RMSD graphs
!   if (testtraj_flag) &
         trajRMSDbuffer(k,Naccept+1) = RMSDbuffer1(1)
!        write(filechannels(1+k),FMT=FMT6) RMSDbuffer1(1)
    
end do

number_of_frames = Totalnumber_of_frames
neighbor_check = number_of_cells
Norder_total(1+order/Ngrid_total) = Norder_total(1+order/Ngrid_total) + 1

if (Ninterpolation == 0) then
    candidate_rmsd = min_rmsd
    return
end if

if ((interpolation_flag).and.(Ninterpolation > 1)) then
    call getInterpolatedGradient(new_coords,gradient)

    call getCoulombMatrix(Natoms,var_coords+new_coords,&
            charges,temp_CM)
    call getCoulombMatrixDiff(Natoms,temp_CM,var_CM,&
            interpolated_CMdiff)

    !Remark: Ken Dill's code uses a version of the RMSD
    !        that divides by N, not N - 1

    min_rmsd = sqrt(sum((new_coords)**2)/(Natoms))
else
    gradient = matmul(Ubuffer1(:,:,1),&
                      gradientbuffer1(:,:,1))

    min_rmsd = RMSDbuffer1(1)

    largest_weighted_rmsd = min_rmsd
    largest_weighted_rmsd2 = min_rmsd**2

    interpolated_CMdiff = CMdiffbuffer1(1)
    
!   !!!!!!!!!!!!!
!   ! TEST START
!   !!!!!!!!!!!!!

!   new_coords = coordsbuffer1(:,:,1)

!   call rmsd_dp(Natoms,new_coords,var_coords,&
!                1,new_U,x_center,y_center,&
!                min_rmsd)

!   do i = 1, 3
!       new_coords(i,:) = new_coords(i,:) - x_center(i)
!   end do

!   new_coords = matmul(Ubuffer1(:,:,1),new_coords)

!   do i = 1, 3
!       new_coords(i,:) = new_coords(i,:) + y_center(i)
!   end do

!   print *, ""
!   print *, ""
!   print *, RMSDbuffer1(1)
!   print *, min_rmsd
!   print *, sqrt(sum((new_coords-var_coords)**2)/Natoms)
!   print *, ""
!   print *, min_rmsd - sqrt(sum((new_coords-var_coords)**2)/Natoms)

!   min_rmsd = RMSDbuffer1(1)

!   !!!!!!!!!!!!!
!   ! TEST END
!   !!!!!!!!!!!!!

end if

candidate_rmsd = RMSDbuffer1(1)
candidate_gradient = matmul(Ubuffer1(:,:,1),&
        gradientbuffer1(:,:,1))

return

end subroutine checkState_new_permute_cap

subroutine checkState_new_cap(vals,coords,gradient,&
                min_rmsd,filechannels,&
                number_of_frames,order,neighbor_check)
use ANALYSIS
use VARIABLES
use PARAMETERS
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
real(dp), dimension(3) :: x_center, y_center
real(dp), allocatable :: g(:,:)
real(dp), dimension(3,3) :: new_U

real(dp),dimension(Nvar) :: var_cell
integer,dimension(Nvar,Norder_max+1) :: var_index
integer,dimension(Nvar,Norder_max+1) :: var_index_diff

integer  :: largest_rmsd_error_index
real(dp) :: largest_rmsd_error

var_filechannel = filechannels(1)
var_coords = coords

!We start off with zero frames having been checked
local_frame_count = 0
Totalnumber_of_frames = 0
order = 0
number_of_cells = 0
neighbor_check = 0
Ninterpolation = 0
stop_flag = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 SUBCELL TARGETING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Retrieve the index of each variable with
!respect to the grid and the real number
!(rounded) that represents that index
do i = 1, Nvar
    !Repeat this for however many orders of cells deep
    !we have been instructed to go
    do j = 1, Norder_max + 1
        var_index(i,j) = int(vals(i) * divisor(i,j))
    end do
end do

!Now, we start iterating over the grids
do k = 1, Ngrid_total
    
    !Streamline the process by storing the
    !path to the grid in a string gridpath3
    write(variable_length_text,FMT=FMT5_variable)&
            Ngrid_text_length
    write(Ngrid_text,FMT="(I0."//&
            trim(adjustl(variable_length_text))//&
            ")") k
    gridpath3 = gridpath0//Ngrid_text//"/grid/"
    
    !The way we check, we check by order first
    do l = 1, Norder_max+1
    
        !The user can specify in what order to
        !check the orders through this array
        Norder = Norder_order(l)
        
        !Read the frames in the cells and
        !process their RMSDs

        !If 
        if (k == grid_addition) then
            if ((l==1) .or. &
                (Totalnumber_of_frames == 0)) then
                call getRMSD_1_cap(&
                        var_index(:,Norder+1),&
                        population)
            else
                call getRMSD_2(&
                        var_index(:,Norder+1),&
                        population)
            end if
            local_frame_count(Norder+1) = &
                    population
        else
            if ((l==1) .or. &
                (Totalnumber_of_frames == 0)) then
                call getRMSD_1_cap(&
                        var_index(:,Norder+1),&
                        population)
            else
            end if
        end if
        
        if (population >= var_overcrowd(Norder+1)) then
            if (Norder < Norder_max) cycle
        else
            Totalnumber_of_frames = &
                Totalnumber_of_frames + population
        end if

        !If the cell is unpopulated or a certain flag is true,
        !then we go ahead and look at the neighbors of this cell
        if ((force_Neighbors) .or. (population == 0)) then
        
            !Integer i keeps track of how far away from the original
            !subcell we are; we look at cells on the 'diamond' surrounding
            !the original subcell
            do i = 1, subcellsearch_max(Norder+1)
        
                var_index_diff = 0
        
                call getRelativeIndex_cap(1,var_index(:,Norder+1),i,&
                                      var_index_diff,stop_flag)
                
                if ((stop_flag) .and. (.not. force_Neighbors)) exit
            end do
        end if

        !If we found a non-empty cell then our search has
        !been over this particular order and does not need
        !to go over other orders
!       if (population > 0) then
        if (Totalnumber_of_frames > 0) then
                order = order + Norder
                exit
        end if
    
    end do
    
    !Record the RMSD encountered here for later analysis:
    !particularly, percent-RMSD graphs
    if ((testtraj_flag).and.(.not.(&
         testtrajDetailedRMSD_flag)))&
         trajRMSDbuffer(k,Naccept+1) = RMSDbuffer1(1)
!        write(filechannels(1+k),FMT=FMT6) RMSDbuffer1(1)
    
end do

number_of_frames = Totalnumber_of_frames
neighbor_check = number_of_cells
Norder_total(1+order/Ngrid_total) = Norder_total(1+order/Ngrid_total) + 1

if (Ninterpolation == 0) then
    candidate_rmsd = min_rmsd
    return
end if

if ((interpolation_flag).and.(Ninterpolation > 1)) then
    call getInterpolatedGradient(new_coords,gradient)

    !Remark: Ken Dill's code uses a version of the RMSD
    !        that divides by N, not N - 1

    min_rmsd = sqrt(sum((new_coords)**2)/(Natoms))
else
    gradient = matmul(Ubuffer1(:,:,1),&
                      gradientbuffer1(:,:,1))

    min_rmsd = RMSDbuffer1(1)

    largest_weighted_rmsd = min_rmsd
    largest_weighted_rmsd2 = min_rmsd**2
end if

candidate_rmsd = RMSDbuffer1(1)
candidate_gradient = matmul(Ubuffer1(:,:,1),&
        gradientbuffer1(:,:,1))

return

end subroutine checkState_new_cap


subroutine checkState_PCM(vals,coords,gradient,&
                min_rmsd,filechannels,&
                number_of_frames,order,neighbor_check)
use ls_rmsd_original
use ANALYSIS
use VARIABLES
use PARAMETERS
use FUNCTIONS
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
real(dp), dimension(3) :: x_center, y_center
real(dp), allocatable :: g(:,:)
real(dp), dimension(3,3) :: new_U

real(dp),dimension(Nvar) :: var_cell
integer,dimension(Nvar,Norder_max+1) :: var_index
integer,dimension(Nvar,Norder_max+1) :: var_index_diff

integer  :: largest_rmsd_error_index
real(dp) :: largest_rmsd_error
real(dp),dimension(Natoms,Natoms) :: temp_CM

var_filechannel = filechannels(1)
var_coords = coords

call getCoulombMatrix(Natoms,var_coords,charges,var_CM)

!We start off with zero frames having been checked
local_frame_count = 0
Totalnumber_of_frames = 0
order = 0
number_of_cells = 0
neighbor_check = 0
Ninterpolation = 0
stop_flag = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 SUBCELL TARGETING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Retrieve the index of each variable with
!respect to the grid and the real number
!(rounded) that represents that index
do i = 1, Nvar
    !Repeat this for however many orders of cells deep
    !we have been instructed to go
    do j = 1, Norder_max + 1
        var_index(i,j) = int(vals(i) * divisor(i,j))
    end do
end do

call shiftMemory((var_index(:,Norder_order(1)+1) -&
                 previous_var_index))
previous_var_index = var_index(:,Norder_order(1)+1)

!Now, we start iterating over the grids
do k = 1, Ngrid_total
    
    !Streamline the process by storing the
    !path to the grid in a string gridpath3
    write(variable_length_text,FMT=FMT5_variable)&
            Ngrid_text_length
    write(Ngrid_text,FMT="(I0."//&
            trim(adjustl(variable_length_text))//&
            ")") k
    gridpath3 = gridpath0//Ngrid_text//"/grid/"
    
    !The way we check, we check by order first
    do l = 1, Norder_max+1
    
        !The user can specify in what order to
        !check the orders through this array
        Norder = Norder_order(l)
        
        !Read the frames in the cells and
        !process their RMSDs

        !If 
        if (k == grid_addition) then
            if ((l==1) .or. &
                (Totalnumber_of_frames == 0)) then
                call getRMSD_1_PCM(&
                        var_index(:,Norder+1),&
                        population)
            else
                call getRMSD_2(&
                        var_index(:,Norder+1),&
                        population)
            end if
            local_frame_count(Norder+1) = &
                    population
        else
            if ((l==1) .or. &
                (Totalnumber_of_frames == 0)) then
                call getRMSD_1_PCM(&
                        var_index(:,Norder+1),&
                        population)
            else
            end if
        end if
        
        if (population >= var_overcrowd(Norder+1)) then
            if (Norder < Norder_max) cycle
        else
            Totalnumber_of_frames = &
                Totalnumber_of_frames + population
        end if
        
        !If the cell is unpopulated or a certain flag is true,
        !then we go ahead and look at the neighbors of this cell
        if ((force_Neighbors) .or. (population == 0)) then
        
            !Integer i keeps track of how far away from the original
            !subcell we are; we look at cells on the 'diamond' surrounding
            !the original subcell
            do i = 1, subcellsearch_max(Norder+1)
        
                var_index_diff = 0
        
                call getRelativeIndex_PCM(1,var_index(:,Norder+1),i,&
                                      var_index_diff,stop_flag)
                
                if ((stop_flag) .and. (.not. force_Neighbors)) exit
            end do
        end if
        

        !If we found a non-empty cell then our search has
        !been over this particular order and does not need
        !to go over other orders
!       if (population > 0) then
        if (Totalnumber_of_frames > 0) then
                order = order + Norder
                exit
        end if
    
    end do

    !Record the RMSD encountered here for later analysis:
    !particularly, percent-RMSD graphs
!   if (testtraj_flag) &
         trajRMSDbuffer(k,Naccept+1) = RMSDbuffer1(1)
!        write(filechannels(1+k),FMT=FMT6) RMSDbuffer1(1)
    
end do

number_of_frames = Totalnumber_of_frames
neighbor_check = number_of_cells
Norder_total(1+order/Ngrid_total) = Norder_total(1+order/Ngrid_total) + 1

if (Ninterpolation == 0) then
    candidate_rmsd = min_rmsd
    return
end if

if ((interpolation_flag).and.(Ninterpolation > 1)) then
    call getInterpolatedGradient(new_coords,gradient)

    call getCoulombMatrix(Natoms,var_coords+new_coords,&
            charges,temp_CM)
    call getCoulombMatrixDiff(Natoms,temp_CM,var_CM,&
            interpolated_CMdiff)

    !Remark: Ken Dill's code uses a version of the RMSD
    !        that divides by N, not N - 1

    min_rmsd = sqrt(sum((new_coords)**2)/(Natoms))
else
    gradient = matmul(Ubuffer1(:,:,1),&
                      gradientbuffer1(:,:,1))

    min_rmsd = RMSDbuffer1(1)

    largest_weighted_rmsd = min_rmsd
    largest_weighted_rmsd2 = min_rmsd**2

    interpolated_CMdiff = CMdiffbuffer1(1)
    
end if

candidate_rmsd = RMSDbuffer1(1)
candidate_gradient = matmul(Ubuffer1(:,:,1),&
        gradientbuffer1(:,:,1))

return

end subroutine checkState_PCM



recursive subroutine getRelativeIndex(currentVar,var_index,&
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

if (currentVar == 1) then
!      stop_flag = .false.
       var_index_diff = 0
       k = 0
else
       k = sum(abs(var_index_diff(1:currentVar-1)))
end if

if (k == N) then
        do nextVar = currentVar, Nvar
                var_index_diff(nextVar) = 0
        end do

        call getRMSD_1(var_index + var_index_diff,population)

        if (population > 0) then
                Totalnumber_of_frames = &
                Totalnumber_of_frames + population

                stop_flag = .true.
        end if

else if (currentVar == Nvar) then

        var_index_diff(currentVar) = N - k

        call getRMSD_1(var_index + var_index_diff,population)

        if (population > 0) then
                Totalnumber_of_frames = &
                Totalnumber_of_frames + population

                stop_flag = .true.
        end if

        var_index_diff(currentVar) = -N + k

        call getRMSD_1(var_index + var_index_diff,population)

        if (population > 0) then
                Totalnumber_of_frames = &
                Totalnumber_of_frames + population

                stop_flag = .true.
        end if

else
        do j = -N + k, N - k
                var_index_diff(currentVar) = j
                nextVar = currentVar + 1
                call getRelativeIndex(nextVar,var_index,&
                                      N,var_index_diff,stop_flag)
        end do
end if

end subroutine getRelativeIndex

recursive subroutine getRelativeIndex_cap(currentVar,var_index,&
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

if (currentVar == 1) then
!      stop_flag = .false.
       var_index_diff = 0
       k = 0
else
       k = sum(abs(var_index_diff(1:currentVar-1)))
end if

if (k == N) then
        do nextVar = currentVar, Nvar
                var_index_diff(nextVar) = 0
        end do

        call getRMSD_1_cap(var_index + var_index_diff,population)

        if (population > 0) then
                Totalnumber_of_frames = &
                Totalnumber_of_frames + population

                stop_flag = .true.
        end if

else if (currentVar == Nvar) then

        var_index_diff(currentVar) = N - k

        call getRMSD_1_cap(var_index + var_index_diff,population)

        if (population > 0) then
                Totalnumber_of_frames = &
                Totalnumber_of_frames + population

                stop_flag = .true.
        end if

        var_index_diff(currentVar) = -N + k

        call getRMSD_1_cap(var_index + var_index_diff,population)

        if (population > 0) then
                Totalnumber_of_frames = &
                Totalnumber_of_frames + population

                stop_flag = .true.
        end if

else
        do j = -N + k, N - k
                var_index_diff(currentVar) = j
                nextVar = currentVar + 1
                call getRelativeIndex_cap(nextVar,var_index,&
                                      N,var_index_diff,stop_flag)
        end do
end if

end subroutine getRelativeIndex_cap

recursive subroutine getRelativeIndex_permute_cap(currentVar,var_index,&
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

if (currentVar == 1) then
!      stop_flag = .false.
       var_index_diff = 0
       k = 0
else
       k = sum(abs(var_index_diff(1:currentVar-1)))
end if

if (k == N) then
        do nextVar = currentVar, Nvar
                var_index_diff(nextVar) = 0
        end do

        call getRMSD_1_permute_cap(var_index + var_index_diff,population)

        if (population > 0) then
                Totalnumber_of_frames = &
                Totalnumber_of_frames + population

                stop_flag = .true.
        end if

else if (currentVar == Nvar) then

        var_index_diff(currentVar) = N - k

        call getRMSD_1_permute_cap(var_index + var_index_diff,population)

        if (population > 0) then
                Totalnumber_of_frames = &
                Totalnumber_of_frames + population

                stop_flag = .true.
        end if

        var_index_diff(currentVar) = -N + k

        call getRMSD_1_permute_cap(var_index + var_index_diff,population)

        if (population > 0) then
                Totalnumber_of_frames = &
                Totalnumber_of_frames + population

                stop_flag = .true.
        end if

else
        do j = -N + k, N - k
                var_index_diff(currentVar) = j
                nextVar = currentVar + 1
                call getRelativeIndex_permute_cap(nextVar,var_index,&
                                      N,var_index_diff,stop_flag)
        end do
end if

end subroutine getRelativeIndex_permute_cap

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

if (currentVar == 1) then
!      stop_flag = .false.
       var_index_diff = 0
       k = 0
else
       k = sum(abs(var_index_diff(1:currentVar-1)))
end if

if (k == N) then
        do nextVar = currentVar, Nvar
                var_index_diff(nextVar) = 0
        end do

        call getRMSD_1_PCM(var_index + var_index_diff,population)

        if (population > 0) then
                Totalnumber_of_frames = &
                Totalnumber_of_frames + population

                stop_flag = .true.
        end if

else if (currentVar == Nvar) then

        var_index_diff(currentVar) = N - k

        call getRMSD_1_PCM(var_index + var_index_diff,population)

        if (population > 0) then
                Totalnumber_of_frames = &
                Totalnumber_of_frames + population

                stop_flag = .true.
        end if

        var_index_diff(currentVar) = -N + k

        call getRMSD_1_PCM(var_index + var_index_diff,population)

        if (population > 0) then
                Totalnumber_of_frames = &
                Totalnumber_of_frames + population

                stop_flag = .true.
        end if

else
        do j = -N + k, N - k
                var_index_diff(currentVar) = j
                nextVar = currentVar + 1
                call getRelativeIndex_PCM(nextVar,var_index,&
                                      N,var_index_diff,stop_flag)
        end do
end if

end subroutine getRelativeIndex_PCM





subroutine getInterpolatedGradient(weighted_coords,weighted_gradient)
use ls_rmsd_original
use FUNCTIONS
use ANALYSIS
use PARAMETERS
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
integer,allocatable :: temp_approximation_index(:)
logical ,allocatable :: temp_acceptable_frame_mask(:)
real(dp),allocatable :: temp_inputCLS(:,:)

integer,allocatable :: RMSD_indexes(:)

integer :: Ninterpolation_true
real(dp) :: weight_threshold = 1.0d-6
real(dp) :: weight
real(dp) :: total_weight

!Incremental integers
integer :: i, j, k
!
!allocate(RMSD_indexes(buffer1_size))
!do i = 1, Ninterpolation
!    RMSD_indexes(i) = i
!end do
!
!call qsort2(RMSDbuffer1,RMSD_indexes,buffer1_size,&
!        1,1,Ninterpolation,1)
!        
!!Create a temporary buffer
!allocate(temp_valsbuffer1(Nvar,buffer1_size),&
!         temp_coordsbuffer1(3,Natoms,buffer1_size),&
!         temp_gradientbuffer1(3,Natoms,buffer1_size),&
!         temp_Ubuffer1(3,3,buffer1_size),&
!         temp_approximation_index(buffer1_size),&
!         temp_acceptable_frame_mask(buffer1_size),&
!         temp_inputCLS(Ncoords+buffer1_size,buffer1_size))
!
!!Store the buffer in the temporary buffer
!temp_valsbuffer1 = valsbuffer1
!temp_coordsbuffer1 = coordsbuffer1
!temp_gradientbuffer1 = gradientbuffer1
!temp_Ubuffer1 = Ubuffer1
!temp_approximation_index = approximation_index
!temp_acceptable_frame_mask = acceptable_frame_mask
!temp_inputCLS = inputCLS
!
!do i = 1, Ninterpolation
!    j = RMSD_indexes(i)
!    valsbuffer1(:,i) = temp_valsbuffer1(:,j)
!    coordsbuffer1(:,:,i) = temp_coordsbuffer1(:,:,j)
!    gradientbuffer1(:,:,i) = temp_gradientbuffer1(:,:,j)
!    Ubuffer1(:,:,i) = temp_Ubuffer1(:,:,j)
!    approximation_index(i) = temp_approximation_index(j)
!    acceptable_frame_mask(i) = temp_acceptable_frame_mask(j)
!    inputCLS(:,i) = temp_inputCLS(:,j)
!end do
!
!deallocate(temp_valsbuffer1,temp_coordsbuffer1,&
!        temp_gradientbuffer1,temp_Ubuffer1,&
!        temp_approximation_index,&
!        temp_acceptable_frame_mask,temp_inputCLS)
!deallocate(RMSD_indexes)
!


Ninterpolation = min(Ninterpolation,&
        Ninterpolation_cap)

allocate(frame_weights(Ninterpolation),&
         outputCLS(Ncoords+Ninterpolation),&
         restraints(1,Ninterpolation),&
         restraint_values(1),&
         minimized_differences2(Ninterpolation,1))

do i = 1, Ninterpolation
        inputCLS(Ncoords+i,:) = 0.0d0
!       inputCLS(Ncoords+i,i) = alpha_ratio*&
!               maxval(abs(inputCLS(1:Ncoords,i)))**2
        inputCLS(Ncoords+i,i) = alpha_ratio*&
                sum(inputCLS(1:Ncoords,i)**2)/Natoms
end do

restraints = 1.0d0
restraint_values = 1.0d0
outputCLS(1:Ncoords) = 0.0d0
outputCLS(Ncoords+1:Ncoords+Ninterpolation) = 0.0d0

call CLS2(inputCLS(1:Ncoords+Ninterpolation,&
          1:Ninterpolation),Ncoords+Ninterpolation,Ninterpolation,&
          restraints,1,restraint_values,&
          outputCLS,frame_weights)

!print *, frame_weights
if (all(frame_weights == 0.0d0)) then
        do i = 1, Ninterpolation
                outputCLS(i) = inputCLS(Ncoords+i,i)
        end do
        frame_weights(minloc(outputCLS(1:Ninterpolation))) = 1.0d0
!       print*, outputCLS(1:Ninterpolation)
!       print*, minloc(outputCLS(1:Ninterpolation))
end if

weighted_coords = 0.0d0
weighted_gradient = 0.0d0
total_weight = 0.0d0
largest_rmsd = 0.0d0
largest_weighted_rmsd = 0.0d0

minimized_differences2 = matmul(inputCLS(Ncoords+1:&
                         Ncoords+Ninterpolation,1:&
                         Ninterpolation),reshape(&
                         frame_weights,(/Ninterpolation,1/)))
!largest_weighted_rmsd2 = maxval(abs(minimized_differences2))**2
largest_weighted_rmsd2 = &
        sqrt(sum((minimized_differences2)**2)/Natoms)/alpha_ratio

do i = 1, Ninterpolation

        weight = frame_weights(i)
!       if (abs(weight) < 1.0d-9) cycle

        total_weight = total_weight + weight
        largest_rmsd = max(largest_rmsd,RMSDbuffer1(i))
        largest_weighted_rmsd = max(largest_weighted_rmsd,&
                abs(weight)*RMSDbuffer1(i))

        weighted_coords = weighted_coords + weight * reshape(&
                          inputCLS(1:Ncoords,i),(/3,Natoms/))
        
        weighted_gradient = weighted_gradient + &
                            weight * matmul(Ubuffer1(:,:,i),&
                                     gradientbuffer1(:,:,i))
end do

!print *, ""
!print *, "           weight threshold: ", weight_threshold
!print *, "total weight of next coords: ", total_weight
!print *, ""

!temp_frame_weights(1:Ninterpolation) = frame_weights

deallocate(frame_weights,outputCLS,&
           restraints,restraint_values,&
           minimized_differences2)

return

end subroutine getInterpolatedGradient


subroutine getWeightedGradient(weighted_coords,weighted_gradient)
use ls_rmsd_original
use ANALYSIS
use PARAMETERS
implicit none

!Information of the frame is held here
real(dp), dimension(3,Natoms), intent(out) :: weighted_coords
real(dp), dimension(3,Natoms), intent(out) :: weighted_gradient

real(dp) :: weight
real(dp) :: total_weight

!Incremental integers
integer :: i, j, k

weighted_coords = 0.0d0
weighted_gradient = 0.0d0
total_weight = 0.0d0

do i = 1, Totalnumber_of_frames

        if (RMSDbuffer1(i) > threshold_rmsd) cycle

        Ninterpolation = Ninterpolation + 1
        
        weight = (RMSDbuffer1(i))**(-interpolation_alpha1)
        
        weighted_coords = weighted_coords + &
                          weight * matmul(Ubuffer1(:,:,i),&
                                     coordsbuffer1(:,:,i))
        
        weighted_gradient = weighted_gradient + &
                            weight * matmul(Ubuffer1(:,:,i),&
                                     gradientbuffer1(:,:,i))
        
        total_weight = total_weight + weight

end do

weighted_coords = (weighted_coords) / total_weight
weighted_gradient = (weighted_gradient) / total_weight

return

end subroutine getWeightedGradient


subroutine getRMSD_1(var_index,population)
use ls_rmsd_original
use ANALYSIS
use PARAMETERS
implicit none

!Inputs for file reading
!character(*),intent(in) :: subcell
integer,dimension(Nvar),intent(in) :: var_index

character(50) :: var_filename
character(150) :: subcell
logical :: subcell_existence

!Variables used in RMSD calculations
real(dp), dimension(3) :: x_center, y_center
real(dp), allocatable :: g(:,:)

!Outputs from file reading
integer,intent(out) :: population

!In case the buffer overfills
real(dp),allocatable :: temp_valsbuffer1(:,:)
real(dp),allocatable :: temp_coordsbuffer1(:,:,:)
real(dp),allocatable :: temp_gradientbuffer1(:,:,:)
real(dp),allocatable :: temp_Ubuffer1(:,:,:)
real(dp),allocatable :: temp_RMSDbuffer1(:)
integer,allocatable :: temp_approximation_index(:)

integer,allocatable :: temp_Ntrajbuffer1(:)

logical ,allocatable :: temp_acceptable_frame_mask(:)
real(dp),allocatable :: temp_inputCLS(:,:)

!Stores values temporarily
real(dp) :: current_rmsd
real(dp),dimension(3,Natoms) :: current_coords

!Incremental integers and iostate checking
integer :: i,j,k,iostate
integer :: endpoint

write(var_filename,FMT=var_multipleFMT&
      (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
        var_index * multiplier(:,Norder+1)

!Construct the subcell filename
subcell = gridpath3//trim(var_filename)

!See whether this subcell exists
inquire(file=trim(subcell),exist=subcell_existence)

if (.not. subcell_existence) then
        population = 0
        return
else
        number_of_cells = number_of_cells + 1
        approximation_index(number_of_cells) = 0
end if

!Open the file corresponding to the cell
if (unreadable_flag) then
        open(var_filechannel,action="read",form="unformatted",&
             file=trim(subcell))
else
        open(var_filechannel,action="read",&
             file=trim(subcell))
end if

!Initialize a variable
population = 1

!Because the buffer may not be large enough to hold all the frames
!and, theoretically, we don't know how many times we need to
!increase the buffer, we need an overarching do loop that is capable
!of increasing the buffer without bounds
do
        !If we have just increased the buffer size, then that means that there
        !are already some frames in the buffer; in this case, do not start at
        !1, but at the current line number (population)
        do k = Totalnumber_of_frames + population, buffer1_size

                !Read the candidate frame
                if (unreadable_flag) then
        
                        !In unformatted files, the first line are the variables
                        !which do not need to be stored
                        read(var_filechannel,iostat=iostate) &
                                (valsbuffer1(i,k),i=1,Nvar)
 
                        !If there are no more lines, stop; the population of the cell should be
                        !one less the number of times that this portion of the loop was called
                        if (iostate /= 0) then
                                population = population - 1
                                exit
                        end if

                        read(var_filechannel) &
                                Ntrajbuffer1(k)
        
                        !The next line is the coordinates
                        read(var_filechannel) &
                               ((current_coords(i,j),i=1,3),j=1,Natoms)
        
                        !In unformatted files, the last line is the gradient
                        read(var_filechannel) &
                                ((gradientbuffer1(i,j,k),i=1,3),j=1,Natoms)
                else
        
                        !In formatted files, everything (variables, coordinates, and gradient)
                        !are stored in one line; FMT1 reads the variables
                        read(var_filechannel,FMT=FMT1,advance="no",iostat=iostate) &
                                (valsbuffer1(i,k),i=1,Nvar)

                        !If there are no more lines, stop; the population of the cell should be
                        !one less the number of times that this portion of the loop was called
                        if (iostate /= 0) then
                                population = population - 1
                                exit
                        end if
        
                        !In formatted files, FMT3 reads the coordinates
                        read(var_filechannel,FMT=FMT7,advance="no") &
                               ((current_coords(i,j),i=1,3),j=1,Natoms)
        
                        !In formatted files, FMT3 reads the gradient as well
                        read(var_filechannel,FMT=FMT3) &
                                ((gradientbuffer1(i,j,k),i=1,3),j=1,Natoms)
                end if

                coordsbuffer1(:,:,k) = current_coords
        
                !Calculate the RMSD between this frame and the incoming frame
                call rmsd_dp(Natoms,current_coords,var_coords,&
                             1,Ubuffer1(:,:,k),x_center,y_center,&
                             current_rmsd)                               !,.false.,g)
!               call rmsd_dp(Natoms,var_coords,current_coords,&
!                            1,Ubuffer1(:,:,k),y_center,x_center,&
!                            current_rmsd)                               !,.false.,g)
                RMSDbuffer1(k) = current_rmsd
        
                !If in the "accept worst" method:
                if (accept_worst) then
        
                        !If the RMSD is too low, reject the candidate frame for it is "too good"
                        if (current_rmsd < candidate_rmsd) then
                
                        !And if the RMSD if lower than the threshold, designate this frame as the "worst"
!                       else if (current_rmsd < threshold_rmsd) then
                        else if ((current_rmsd < threshold_rmsd).and.&
                                 (current_rmsd >= inner_threshold)) then
                
                                approximation_index(number_of_cells) = k
                                candidate_rmsd = current_rmsd
                
                                !In special cases, we exit immediately afterwards
                                if (accept_first) exit
                
                        !If the RMSD is too high, also reject the candidate frame
                        else
                
                        end if
        
                !If in the "accept best" method:
                else
                        if ((interpolation_flag).and.&
                            (current_rmsd < threshold_rmsd)) then

                                acceptable_frame_mask(k) = .true.

                                Ninterpolation = Ninterpolation + 1

                                do i = 1, 3
                                        current_coords(i,:) = &
                                        current_coords(i,:) - x_center(i)
                                end do

                                current_coords = matmul(&
                                        Ubuffer1(:,:,k),current_coords)

                                do i = 1, 3
                                        current_coords(i,:) = &
                                        current_coords(i,:) + y_center(i)
                                end do

                                inputCLS(1:Ncoords,Ninterpolation) =&
                                           reshape(current_coords-var_coords,(/Ncoords/))

!                               temp_rmsd_weights(Ninterpolation) = current_rmsd

                        end if

                        !If the RMSD is low enough, designate this frame as the "best"
!                       if (current_rmsd < candidate_rmsd) then
                        if ((current_rmsd < candidate_rmsd).and.&
                            (current_rmsd >= inner_threshold)) then
                        
                                approximation_index(number_of_cells) = k
                                candidate_rmsd = current_rmsd
                
                                !In special cases, we exit immediately afterwards
                                if (accept_first) exit
                
                        !Otherwise, do nothing
                        else
                
                        end if
        
                end if
        
                !Increment the number of frames visited
                population = population + 1
        end do
        
        !If the last line of the file has been read, we can exit out of the loop
        if (iostate /= 0) exit

        !Otherwise, the buffer has run out of room, and we need to first increase
        !it, then continuereading in frames
        
        !Create a temporary buffer
        allocate(temp_valsbuffer1(Nvar,buffer1_size),&
                 temp_coordsbuffer1(3,Natoms,buffer1_size),&
                 temp_gradientbuffer1(3,Natoms,buffer1_size),&
                 temp_Ubuffer1(3,3,buffer1_size),&
                 temp_RMSDbuffer1(buffer1_size),&
                 temp_approximation_index(buffer1_size),&
                 temp_acceptable_frame_mask(buffer1_size),&
                 temp_inputCLS(Ncoords+buffer1_size,buffer1_size))

        allocate(temp_Ntrajbuffer1(buffer1_size))
        
        !Store the buffer in the temporary buffer
        temp_valsbuffer1 = valsbuffer1
        temp_coordsbuffer1 = coordsbuffer1
        temp_gradientbuffer1 = gradientbuffer1
        temp_Ubuffer1 = Ubuffer1
        temp_RMSDbuffer1 = RMSDbuffer1
        temp_approximation_index = approximation_index
        temp_acceptable_frame_mask = acceptable_frame_mask
        temp_inputCLS = inputCLS

        temp_Ntrajbuffer1 = Ntrajbuffer1
        
        !For now, we simply double the buffer size each time it overfills
        deallocate(valsbuffer1,coordsbuffer1,gradientbuffer1,Ubuffer1,RMSDbuffer1,&
                   approximation_index,acceptable_frame_mask,inputCLS)

        deallocate(Ntrajbuffer1)
        allocate(valsbuffer1(Nvar,buffer1_size*2),&
                 coordsbuffer1(3,Natoms,buffer1_size*2),&
                 gradientbuffer1(3,Natoms,buffer1_size*2),&
                 Ubuffer1(3,3,buffer1_size*2),&
                 RMSDbuffer1(buffer1_size*2),&
                 approximation_index(buffer1_size*2),&
                 acceptable_frame_mask(buffer1_size*2),&
                 inputCLS(Ncoords+buffer1_size*2,buffer1_size*2))

        allocate(Ntrajbuffer1(buffer1_size*2))
        
        acceptable_frame_mask = .false.

        !And reload all the frames back into the buffer
        valsbuffer1(:,1:buffer1_size) = temp_valsbuffer1
        coordsbuffer1(:,:,1:buffer1_size) = temp_coordsbuffer1
        gradientbuffer1(:,:,1:buffer1_size) = temp_gradientbuffer1
        Ubuffer1(:,:,1:buffer1_size) = temp_Ubuffer1
        RMSDbuffer1(1:buffer1_size) = temp_RMSDbuffer1
        approximation_index(1:buffer1_size) = temp_approximation_index
        acceptable_frame_mask(1:buffer1_size) = temp_acceptable_frame_mask
        inputCLS(1:Ncoords+buffer1_size,1:buffer1_size) = temp_inputCLS

        Ntrajbuffer1(1:buffer1_size) = temp_Ntrajbuffer1

        !Destroy the temporary buffer
        deallocate(temp_valsbuffer1,temp_coordsbuffer1,temp_gradientbuffer1,&
                   temp_Ubuffer1,temp_RMSDbuffer1,temp_approximation_index,&
                   temp_acceptable_frame_mask,temp_inputCLS)

        deallocate(temp_Ntrajbuffer1)

        !Permanently increase the buffer size so this will not
        !have to happen next time
        buffer1_size = buffer1_size*2
end do
close(var_filechannel)

end subroutine getRMSD_1

subroutine getRMSD_1_permute(var_index,population)
use ls_rmsd_original
use PARAMETERS
use PHYSICS
use ANALYSIS
implicit none

!Inputs for file reading
!character(*),intent(in) :: subcell
integer,dimension(Nvar),intent(in) :: var_index

character(50) :: var_filename
character(150) :: subcell
logical :: subcell_existence

!Variables used in RMSD calculations
real(dp), dimension(3) :: x_center, y_center
real(dp), allocatable :: g(:,:)

!Outputs from file reading
integer,intent(out) :: population

!In case the buffer overfills
real(dp),allocatable :: temp_valsbuffer1(:,:)
real(dp),allocatable :: temp_coordsbuffer1(:,:,:)
real(dp),allocatable :: temp_gradientbuffer1(:,:,:)
real(dp),allocatable :: temp_Ubuffer1(:,:,:)
real(dp),allocatable :: temp_RMSDbuffer1(:)
integer,allocatable :: temp_approximation_index(:)

integer,allocatable :: temp_Ntrajbuffer1(:)

logical ,allocatable :: temp_acceptable_frame_mask(:)
real(dp),allocatable :: temp_inputCLS(:,:)

!Stores values temporarily
real(dp) :: current_rmsd
real(dp),dimension(3,Natoms) :: current_coords,temp_coords
real(dp),dimension(3,Natoms) :: current_gradient
real(dp),dimension(Nvar) :: current_vals
integer :: current_Ntraj

!Incremental integers and iostate checking
integer :: n,i,j,k,iostate
integer :: endpoint

write(var_filename,FMT=var_multipleFMT&
      (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
        var_index * multiplier(:,Norder+1)

!Construct the subcell filename
subcell = gridpath3//trim(var_filename)

!See whether this subcell exists
inquire(file=trim(subcell),exist=subcell_existence)

if (.not. subcell_existence) then
        population = 0
        return
else
        number_of_cells = number_of_cells + 1
        approximation_index(number_of_cells) = 0
end if

!Open the file corresponding to the cell
if (unreadable_flag) then
        open(var_filechannel,action="read",form="unformatted",&
             file=trim(subcell))
else
        open(var_filechannel,action="read",&
             file=trim(subcell))
end if

!Initialize a variable
!population = 1
population = 0
n = 1

!Because the buffer may not be large enough to hold all the frames
!and, theoretically, we don't know how many times we need to
!increase the buffer, we need an overarching do loop that is capable
!of increasing the buffer without bounds
do
        !If we have just increased the buffer size, then that means that there
        !are already some frames in the buffer; in this case, do not start at
        !1, but at the current line number (population)
!       do k = Totalnumber_of_frames + population, buffer1_size
!       do k = (Totalnumber_of_frames + population) * &
!              Nindistinguishables + n, buffer1_size
!       do k = Totalnumber_of_frames + Ninterpolation + 1, buffer1_size
        do k = Ninterpolation + 1, buffer1_size

        do

!       n = modulo(population,Nindistinguishables) + 1

        if (n == 1) then

                !Read the candidate frame
                if (unreadable_flag) then
        
                        !In unformatted files, the first line are the variables
                        !which do not need to be stored
                        read(var_filechannel,iostat=iostate) &
                                (current_vals(i),i=1,Nvar)
!                               (valsbuffer1(i,k),i=1,Nvar)
 
                        !If there are no more lines, stop; the population of the cell should be
                        !one less the number of times that this portion of the loop was called
                        if (iostate /= 0) then
!                               population = population - 1
                                exit
                        end if

!                       read(var_filechannel) &
!                               current_Ntraj
        
                        !The next line is the coordinates
                        read(var_filechannel) &
                               ((current_coords(i,j),i=1,3),j=1,Natoms)
        
                        !In unformatted files, the last line is the gradient
                        read(var_filechannel) &
                                ((current_gradient(i,j),i=1,3),j=1,Natoms)
!                               ((gradientbuffer1(i,j,k),i=1,3),j=1,Natoms)
                else
        
                        !In formatted files, everything (variables, coordinates, and gradient)
                        !are stored in one line; FMT1 reads the variables
                        read(var_filechannel,FMT=FMT1,advance="no",iostat=iostate) &
                                (current_vals(i),i=1,Nvar)
!                               (valsbuffer1(i,k),i=1,Nvar)

                        !If there are no more lines, stop; the population of the cell should be
                        !one less the number of times that this portion of the loop was called
                        if (iostate /= 0) then
!                               population = population - 1
                                exit
                        end if
        
                        !In formatted files, FMT3 reads the coordinates
                        read(var_filechannel,FMT=FMT7,advance="no") &
                               ((current_coords(i,j),i=1,3),j=1,Natoms)
        
                        !In formatted files, FMT3 reads the gradient as well
                        read(var_filechannel,FMT=FMT3) &
                                ((current_gradient(i,j),i=1,3),j=1,Natoms)
!                               ((gradientbuffer1(i,j,k),i=1,3),j=1,Natoms)
                end if

        end if

        BOND_LABELLING_DATA = INDISTINGUISHABLES(n,:)

                n = n + 1
                if (n > Nindistinguishables) then

                        !Increment the number of frames visited
                        population = population + 1
                        n = 1
        
                end if

!               valsbuffer1(:,k) = current_vals
!               Ntrajbuffer1(k) = current_Ntraj

!               do j = 1, Natoms
!                    coordsbuffer1(:,j,k) = &
!                        current_coords(:,BOND_LABELLING_DATA(j))
!                    gradientbuffer1(:,j,k) = &
!                        current_gradient(:,BOND_LABELLING_DATA(j))
!               end do

                do j = 1, Natoms
                     temp_coords(:,j) = &
                         current_coords(:,BOND_LABELLING_DATA(j))
                end do

!               coordsbuffer1(:,:,k) = current_coords
!               gradientbuffer1(:,:,k) = current_gradient
        
                !Calculate the RMSD between this frame and the incoming frame
                call rmsd_dp(Natoms,temp_coords,var_coords,&
                             1,Ubuffer1(:,:,k),x_center,y_center,&
                             current_rmsd)                               !,.false.,g)
!               call rmsd_dp(Natoms,coordsbuffer1(:,:,k),var_coords,&
!                            1,Ubuffer1(:,:,k),x_center,y_center,&
!                            current_rmsd)                               !,.false.,g)
!               call rmsd_dp(Natoms,current_coords,var_coords,&
!                            1,Ubuffer1(:,:,k),x_center,y_center,&
!                            current_rmsd)                               !,.false.,g)
!               call rmsd_dp(Natoms,var_coords,current_coords,&
!                            1,Ubuffer1(:,:,k),y_center,x_center,&
!                            current_rmsd)                               !,.false.,g)
!               RMSDbuffer1(k) = current_rmsd
        
                !If in the "accept worst" method:
                if (accept_worst) then
        
                        !If the RMSD is too low, reject the candidate frame for it is "too good"
                        if (current_rmsd < candidate_rmsd) then
                
                        !And if the RMSD if lower than the threshold, designate this frame as the "worst"
!                       else if (current_rmsd < threshold_rmsd) then
                        else if ((current_rmsd < threshold_rmsd).and.&
                                 (current_rmsd >= inner_threshold)) then
                
                                approximation_index(number_of_cells) = k
                                candidate_rmsd = current_rmsd
                
                                !In special cases, we exit immediately afterwards
                                if (accept_first) exit
                
                        !If the RMSD is too high, also reject the candidate frame
                        else
                
                        end if
        
                !If in the "accept best" method
                !and the RMSD is low enough:
                else if (current_rmsd < threshold_rmsd) then

                        !If the RMSD is low enough, designate this frame as the "best"
!                       if (current_rmsd < candidate_rmsd) then
                        if ((current_rmsd < candidate_rmsd).and.&
                            (current_rmsd >= inner_threshold)) then
                        
                                approximation_index(number_of_cells) = k
                                candidate_rmsd = current_rmsd
                
                                !In special cases, we exit immediately afterwards
!                               if (accept_first) exit
                                if (accept_first) iostate = 1
                
                        end if

!                       if ((interpolation_flag).and.&
!                           (current_rmsd < threshold_rmsd)) then
!                       if (&
!                           (current_rmsd < threshold_rmsd)) then
!                       if (((interpolation_flag).or.&
!                            (gather_interpolation_flag)).and.&
!                           (current_rmsd < threshold_rmsd)) then
                        if ((interpolation_flag).or.&
                             (gather_interpolation_flag)) then

                        RMSDbuffer1(k) = current_rmsd
                        valsbuffer1(:,k) = current_vals
!                       Ntrajbuffer1(k) = current_Ntraj
        
                        do j = 1, Natoms
                             coordsbuffer1(:,j,k) = &
                                 current_coords(:,BOND_LABELLING_DATA(j))
                             gradientbuffer1(:,j,k) = &
                                 current_gradient(:,BOND_LABELLING_DATA(j))
                        end do

                        acceptable_frame_mask(k) = .true.

                        Ninterpolation = Ninterpolation + 1

!                       temp_coords = coordsbuffer1(:,:,k)

                        do i = 1, 3
                                temp_coords(i,:) = &
                                temp_coords(i,:) - x_center(i)
                        end do

                        temp_coords(:,:) = matmul(&
                                Ubuffer1(:,:,k),&
                                temp_coords(:,:))

                        do i = 1, 3
                                temp_coords(i,:) = &
                                temp_coords(i,:) + y_center(i)
                        end do

                        inputCLS(1:Ncoords,Ninterpolation) =&
                                   reshape(temp_coords(:,:) -&
                                   var_coords,(/Ncoords/))

!                       temp_rmsd_weights(Ninterpolation) = current_rmsd

                        exit

                        end if
        
                end if

        end do

        !If the last line of the file has been read, we can exit out of the loop
        if (iostate /= 0) exit

        end do
        
        !If the last line of the file has been read, we can exit out of the loop
        if (iostate /= 0) exit

        print *, buffer1_size

        !Otherwise, the buffer has run out of room, and we need to first increase
        !it, then continuereading in frames
        
        !Create a temporary buffer
        allocate(temp_valsbuffer1(Nvar,buffer1_size),&
                 temp_coordsbuffer1(3,Natoms,buffer1_size),&
                 temp_gradientbuffer1(3,Natoms,buffer1_size),&
                 temp_Ubuffer1(3,3,buffer1_size),&
                 temp_RMSDbuffer1(buffer1_size),&
                 temp_approximation_index(buffer1_size),&
                 temp_acceptable_frame_mask(buffer1_size),&
                 temp_inputCLS(Ncoords+buffer1_size,buffer1_size))

!       allocate(temp_Ntrajbuffer1(buffer1_size))
        
        !Store the buffer in the temporary buffer
        temp_valsbuffer1 = valsbuffer1
        temp_coordsbuffer1 = coordsbuffer1
        temp_gradientbuffer1 = gradientbuffer1
        temp_Ubuffer1 = Ubuffer1
        temp_RMSDbuffer1 = RMSDbuffer1
        temp_approximation_index = approximation_index
        temp_acceptable_frame_mask = acceptable_frame_mask
        temp_inputCLS = inputCLS

!       temp_Ntrajbuffer1 = Ntrajbuffer1
        
        !For now, we simply double the buffer size each time it overfills
        deallocate(valsbuffer1,coordsbuffer1,gradientbuffer1,Ubuffer1,RMSDbuffer1,&
                   approximation_index,acceptable_frame_mask,inputCLS)

!       deallocate(Ntrajbuffer1)
        allocate(valsbuffer1(Nvar,buffer1_size*2),&
                 coordsbuffer1(3,Natoms,buffer1_size*2),&
                 gradientbuffer1(3,Natoms,buffer1_size*2),&
                 Ubuffer1(3,3,buffer1_size*2),&
                 RMSDbuffer1(buffer1_size*2),&
                 approximation_index(buffer1_size*2),&
                 acceptable_frame_mask(buffer1_size*2),&
                 inputCLS(Ncoords+buffer1_size*2,buffer1_size*2))

!       allocate(Ntrajbuffer1(buffer1_size*2))
        
        acceptable_frame_mask = .false.

        !And reload all the frames back into the buffer
        valsbuffer1(:,1:buffer1_size) = temp_valsbuffer1
        coordsbuffer1(:,:,1:buffer1_size) = temp_coordsbuffer1
        gradientbuffer1(:,:,1:buffer1_size) = temp_gradientbuffer1
        Ubuffer1(:,:,1:buffer1_size) = temp_Ubuffer1
        RMSDbuffer1(1:buffer1_size) = temp_RMSDbuffer1
        approximation_index(1:buffer1_size) = temp_approximation_index
        acceptable_frame_mask(1:buffer1_size) = temp_acceptable_frame_mask
        inputCLS(1:Ncoords+buffer1_size,1:buffer1_size) = temp_inputCLS

!       Ntrajbuffer1(1:buffer1_size) = temp_Ntrajbuffer1

        !Destroy the temporary buffer
        deallocate(temp_valsbuffer1,temp_coordsbuffer1,temp_gradientbuffer1,&
                   temp_Ubuffer1,temp_RMSDbuffer1,temp_approximation_index,&
                   temp_acceptable_frame_mask,temp_inputCLS)

!       deallocate(temp_Ntrajbuffer1)

        !Permanently increase the buffer size so this will not
        !have to happen next time
        buffer1_size = buffer1_size*2
end do
close(var_filechannel)


end subroutine getRMSD_1_permute

subroutine getRMSD_1_permute_cap(var_index,population)
use ls_rmsd_original
use PARAMETERS
use PHYSICS
use ANALYSIS
implicit none

!Inputs for file reading
!character(*),intent(in) :: subcell
integer,dimension(Nvar),intent(in) :: var_index

character(50) :: var_filename
character(150) :: subcell
logical :: subcell_existence

!Variables used in RMSD calculations
real(dp), dimension(3) :: x_center, y_center
real(dp), allocatable :: g(:,:)

!Outputs from file reading
integer,intent(out) :: population

!Stores values temporarily
real(dp) :: current_rmsd
real(dp),dimension(3,Natoms) :: current_coords,temp_coords
real(dp),dimension(3,Natoms) :: current_gradient
real(dp),dimension(Nvar) :: current_vals
real(dp),dimension(3,3) :: current_U

integer :: current_Ntraj
real(dp),dimension(Natoms,Natoms) :: temp_CM
real(dp) :: current_CMdiff

!Incremental integers and iostate checking
integer :: n,i,j,k,iostate
integer :: endpoint

write(var_filename,FMT=var_multipleFMT&
      (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
        var_index * multiplier(:,Norder+1)

!Construct the subcell filename
subcell = gridpath3//trim(var_filename)

!See whether this subcell exists
inquire(file=trim(subcell),exist=subcell_existence)

if (.not. subcell_existence) then
    population = 0
    return
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

!Initialize a variable
population = 0
n = Nindistinguishables

do

    if (n == Nindistinguishables) then

        !Read the candidate frame
        if (unreadable_flag) then
    
            !In unformatted files, the first line are the variables
            !which do not need to be stored
            read(var_filechannel,iostat=iostate) &
                    (current_vals(i),i=1,Nvar)
 
            !If there are no more lines, stop; the population of the cell should be
            !one less the number of times that this portion of the loop was called
            if (iostate /= 0) exit

!                       read(var_filechannel) &
!                               current_Ntraj
    
            !The next line is the coordinates
            read(var_filechannel) &
                   ((current_coords(i,j),i=1,3),j=1,Natoms)
    
            !In unformatted files, the last line is the gradient
            read(var_filechannel) &
                    ((current_gradient(i,j),i=1,3),j=1,Natoms)
        else
    
            !In formatted files, everything (variables, coordinates, and gradient)
            !are stored in one line; FMT1 reads the variables
            read(var_filechannel,FMT=FMT1,advance="no",iostat=iostate) &
                    (current_vals(i),i=1,Nvar)

            !If there are no more lines, stop; the population of the cell should be
            !one less the number of times that this portion of the loop was called
            if (iostate /= 0) exit
    
            !In formatted files, FMT3 reads the coordinates
            read(var_filechannel,FMT=FMT7,advance="no") &
                   ((current_coords(i,j),i=1,3),j=1,Natoms)
    
            !In formatted files, FMT3 reads the gradient as well
            read(var_filechannel,FMT=FMT3) &
                    ((current_gradient(i,j),i=1,3),j=1,Natoms)
        end if

        n = 0
        population = population + 1

    end if

    n = n + 1

    BOND_LABELLING_DATA = INDISTINGUISHABLES(n,:)

    do j = 1, Natoms
         temp_coords(:,j) = &
             current_coords(:,BOND_LABELLING_DATA(j))
    end do

    !Calculate the RMSD between this frame and the incoming frame
    call rmsd_dp(Natoms,temp_coords,var_coords,&
                 1,current_U,x_center,y_center,&
                 current_rmsd)                               !,.false.,g)

    call getCoulombMatrix(Natoms,temp_coords,charges,temp_CM)

    call getCoulombMatrixDiff(Natoms,temp_CM,var_CM,current_CMdiff)

    !If in the "accept best" method
    !and the RMSD is low enough:
    if ((current_rmsd < threshold_rmsd).and.&
        (current_rmsd >= inner_threshold)) then

        if (accept_first) iostate = 1

        if (accept_worst) then
            do i = 1, Ninterpolation+1
                if (current_rmsd > RMSDbuffer1(i)) exit
            end do
        else
            do i = 1, Ninterpolation+1
                if (current_rmsd < RMSDbuffer1(i)) exit
            end do
        end if

        if (Ninterpolation < Ninterpolation_max) &
            Ninterpolation = Ninterpolation + 1

        do j = Ninterpolation, i+1, -1
            valsbuffer1(:,j) = valsbuffer1(:,j-1)
            coordsbuffer1(:,:,j) = coordsbuffer1(:,:,j-1)
            gradientbuffer1(:,:,j) = gradientbuffer1(:,:,j-1)
            Ubuffer1(:,:,j) = Ubuffer1(:,:,j-1)
            RMSDbuffer1(j) = RMSDbuffer1(j-1)
            CMdiffbuffer1(j) = CMdiffbuffer1(j-1)
            inputCLS(:,j) = inputCLS(:,j-1)
        end do

        do k = 1, Natoms
            coordsbuffer1(:,k,i) = &
                current_coords(:,BOND_LABELLING_DATA(k))
            gradientbuffer1(:,k,i) = &
                current_gradient(:,BOND_LABELLING_DATA(k))
        end do

        valsbuffer1(:,i) = current_vals
        Ubuffer1(:,:,i) = current_U
        RMSDbuffer1(i) = current_rmsd
        CMdiffbuffer1(i) = current_CMdiff

        temp_coords = coordsbuffer1(:,:,i)

        do j = 1, 3
            temp_coords(j,:) = &
            temp_coords(j,:) - x_center(j)
        end do

        temp_coords(:,:) = matmul(&
                Ubuffer1(:,:,i),&
                temp_coords(:,:))

        do j = 1, 3
            temp_coords(j,:) = &
            temp_coords(j,:) + y_center(j)
        end do

        inputCLS(1:Ncoords,i) =&
                   reshape(temp_coords(:,:) - &
                   var_coords,(/Ncoords/))

    end if

    !If the last line of the file has been read, we can exit out of the loop
    if (iostate /= 0) exit

end do

close(var_filechannel)


end subroutine getRMSD_1_permute_cap

subroutine getRMSD_1_cap(var_index,population)
use ls_rmsd_original
use PARAMETERS
use PHYSICS
use ANALYSIS
implicit none

!Inputs for file reading
!character(*),intent(in) :: subcell
integer,dimension(Nvar),intent(in) :: var_index

character(50) :: var_filename
character(150) :: subcell
logical :: subcell_existence

!Variables used in RMSD calculations
real(dp), dimension(3) :: x_center, y_center
real(dp), allocatable :: g(:,:)

!Outputs from file reading
integer,intent(out) :: population

!Stores values temporarily
real(dp) :: current_rmsd
real(dp),dimension(3,Natoms) :: current_coords
real(dp),dimension(3,Natoms) :: current_gradient
real(dp),dimension(Nvar) :: current_vals
real(dp),dimension(3,3) :: current_U
integer :: current_Ntraj

!Incremental integers and iostate checking
integer :: i,j,k,iostate
integer :: endpoint

write(var_filename,FMT=var_multipleFMT&
      (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
        var_index * multiplier(:,Norder+1)

!Construct the subcell filename
subcell = gridpath3//trim(var_filename)

!See whether this subcell exists
inquire(file=trim(subcell),exist=subcell_existence)

if (.not. subcell_existence) then
    population = 0
    return
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

!Initialize a variable
population = 0

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

!                       read(var_filechannel) &
!                               current_Ntraj

        !The next line is the coordinates
        read(var_filechannel) &
               ((current_coords(i,j),i=1,3),j=1,Natoms)

        !In unformatted files, the last line is the gradient
        read(var_filechannel) &
                ((current_gradient(i,j),i=1,3),j=1,Natoms)
    else

        !In formatted files, everything (variables, coordinates, and gradient)
        !are stored in one line; FMT1 reads the variables
        read(var_filechannel,FMT=FMT1,advance="no",iostat=iostate) &
                (current_vals(i),i=1,Nvar)

        !If there are no more lines, stop; the population of the cell should be
        !one less the number of times that this portion of the loop was called
        if (iostate /= 0) exit

        !In formatted files, FMT3 reads the coordinates
        read(var_filechannel,FMT=FMT7,advance="no") &
               ((current_coords(i,j),i=1,3),j=1,Natoms)

        !In formatted files, FMT3 reads the gradient as well
        read(var_filechannel,FMT=FMT3) &
                ((current_gradient(i,j),i=1,3),j=1,Natoms)
    end if

    population = population + 1

    !Calculate the RMSD between this frame and the incoming frame
    call rmsd_dp(Natoms,current_coords,var_coords,&
                 1,current_U,x_center,y_center,&
                 current_rmsd)                               !,.false.,g)

    !If in the "accept best" method
    !and the RMSD is low enough:
    if ((current_rmsd < threshold_rmsd).and.&
        (current_rmsd >= inner_threshold)) then

        if (accept_first) iostate = 1

        if (accept_worst) then
            do i = 1, Ninterpolation+1
                if (current_rmsd > RMSDbuffer1(i)) exit
            end do
        else
            do i = 1, Ninterpolation+1
                if (current_rmsd < RMSDbuffer1(i)) exit
            end do
        end if

        if (Ninterpolation < Ninterpolation_max) &
            Ninterpolation = Ninterpolation + 1

        do j = Ninterpolation, i+1, -1
            valsbuffer1(:,j) = valsbuffer1(:,j-1)
            coordsbuffer1(:,:,j) = coordsbuffer1(:,:,j-1)
            gradientbuffer1(:,:,j) = gradientbuffer1(:,:,j-1)
            Ubuffer1(:,:,j) = Ubuffer1(:,:,j-1)
            RMSDbuffer1(j) = RMSDbuffer1(j-1)
            inputCLS(:,j) = inputCLS(:,j-1)
        end do

        do k = 1, Natoms
            coordsbuffer1(:,k,i) = &
                current_coords(:,k)
            gradientbuffer1(:,k,i) = &
                current_gradient(:,k)
        end do

        valsbuffer1(:,i) = current_vals
        Ubuffer1(:,:,i) = current_U
        RMSDbuffer1(i) = current_rmsd

        do j = 1, 3
            current_coords(j,:) = &
            current_coords(j,:) - x_center(j)
        end do

        current_coords(:,:) = matmul(&
                Ubuffer1(:,:,i),&
                current_coords(:,:))

        do j = 1, 3
            current_coords(j,:) = &
            current_coords(j,:) + y_center(j)
        end do

        inputCLS(1:Ncoords,i) =&
                   reshape(current_coords(:,:) - &
                   var_coords,(/Ncoords/))

    end if

end do

close(var_filechannel)


end subroutine getRMSD_1_cap

subroutine getRMSD_1_PCM(var_index,population)
use ls_rmsd_original
use PARAMETERS
use PHYSICS
use ANALYSIS
use FUNCTIONS
implicit none

!Inputs for file reading
!character(*),intent(in) :: subcell
integer,dimension(Nvar),intent(in) :: var_index
integer :: single_index

character(50) :: var_filename
character(150) :: subcell
logical :: subcell_existence
logical :: stop_flag

!Variables used in RMSD calculations
real(dp), dimension(3) :: x_center, y_center
real(dp), allocatable :: g(:,:)

!Outputs from file reading
integer,intent(out) :: population

!Stores values temporarily
real(dp),dimension(Nvar) :: current_vals
real(dp) :: temp_rmsd
real(dp),dimension(3,Natoms) :: temp_coords
real(dp),dimension(3,Natoms) :: current_coords,current_gradient
real(dp),dimension(3,3) :: temp_U

real(dp),dimension(Natoms,Natoms) :: temp_CM
real(dp) :: temp_CMdiff

!Incremental integers and iostate checking
integer :: n,i,j,k,iostate
integer :: endpoint

call getFlattened(Nvar,&
        (var_index-previous_var_index),&
        single_index)

population = populationbuffer2(single_index)

if (population < 0) then

    write(var_filename,FMT=var_multipleFMT&
          (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
            var_index * multiplier(:,Norder+1)
    
    !Construct the subcell filename
    subcell = gridpath3//trim(var_filename)
    
    !See whether this subcell exists
    inquire(file=trim(subcell),exist=subcell_existence)
    
    if (.not. subcell_existence) then
        population = 0
        populationbuffer2(single_index) = 0
        return
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
    
    !Initialize a variable
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
    
            !In formatted files, FMT3 reads the coordinates
            read(var_filechannel,FMT=FMT7,advance="no") &
                   ((coordsbuffer2(i,j,population,single_index),i=1,3),j=1,Natoms)
    
            !In formatted files, FMT3 reads the gradient as well
            read(var_filechannel,FMT=FMT3) &
                    ((gradientbuffer2(i,j,population,single_index),i=1,3),j=1,Natoms)
        end if
 
        do n = 1, Nindistinguishables
    
            BOND_LABELLING_DATA = INDISTINGUISHABLES(n,:)
        
            do j = 1, Natoms
                 temp_coords(:,j) = &
                     coordsbuffer2(:,BOND_LABELLING_DATA(j),population,single_index)
            end do
 
            !Calculate the RMSD between this frame and the incoming frame
            call rmsd_dp(Natoms,temp_coords,var_coords,&
                         1,temp_U,x_center,y_center,&
                         temp_rmsd)                               !,.false.,g)

            call getCoulombMatrix(Natoms,temp_coords,charges,temp_CM)
        
            call getCoulombMatrixDiff(Natoms,temp_CM,var_CM,temp_CMdiff)

!           do i = 1, 3
!               temp_U(i,:) = 0.0d0
!               x_center(i) = 0.0d0
!               y_center(i) = 0.0d0
!           end do
!           temp_rmsd = rand() * threshold_rmsd * 4
!           temp_CMdiff = 0.0d0
        
            !If in the "accept best" method
            !and the RMSD is low enough:
            if ((temp_rmsd < threshold_rmsd).and.&
                (temp_rmsd >= inner_threshold)) then
        
                if (accept_first) iostate = 1
        
                if (accept_worst) then
                    do i = 1, Ninterpolation+1
                        if (temp_rmsd > RMSDbuffer1(i)) exit
                    end do
                else
!                   do i = 1, Ninterpolation+1
!                       if (temp_rmsd < RMSDbuffer1(i)) exit
!                   end do
                    do i = 1, Ninterpolation+1
                        if (temp_CMdiff < CMdiffbuffer1(i)) exit
                    end do
                end if
        
                if (Ninterpolation < Ninterpolation_max) &
                    Ninterpolation = Ninterpolation + 1
        
                do j = Ninterpolation, i+1, -1
                    valsbuffer1(:,j) = valsbuffer1(:,j-1)
                    coordsbuffer1(:,:,j) = coordsbuffer1(:,:,j-1)
                    gradientbuffer1(:,:,j) = gradientbuffer1(:,:,j-1)
                    Ubuffer1(:,:,j) = Ubuffer1(:,:,j-1)
                    RMSDbuffer1(j) = RMSDbuffer1(j-1)
                    CMdiffbuffer1(j) = CMdiffbuffer1(j-1)
                    inputCLS(:,j) = inputCLS(:,j-1)
                end do
        
                do j = 1, Natoms
                    coordsbuffer1(:,j,i) = &
                        coordsbuffer2(:,BOND_LABELLING_DATA(j),population,single_index)
                    gradientbuffer1(:,j,i) = &
                        gradientbuffer2(:,BOND_LABELLING_DATA(j),population,single_index)
                end do
        
                valsbuffer1(:,i) = valsbuffer2(:,population,single_index)
                Ubuffer1(:,:,i) = temp_U
                RMSDbuffer1(i) = temp_rmsd
                CMdiffbuffer1(i) = temp_CMdiff
        
                temp_coords = coordsbuffer1(:,:,i)
        
                do j = 1, 3
                    temp_coords(j,:) = &
                    temp_coords(j,:) - x_center(j)
                end do
        
                temp_coords(:,:) = matmul(&
                        Ubuffer1(:,:,i),&
                        temp_coords(:,:))
        
                do j = 1, 3
                    temp_coords(j,:) = &
                    temp_coords(j,:) + y_center(j)
                end do
        
                inputCLS(1:Ncoords,i) =&
                           reshape(temp_coords(:,:) - &
                           var_coords,(/Ncoords/))
            end if
    
        end do
    
        if (population == buffer2_size) then
            stop_flag = .false.
            exit
        end if
    
    end do

    populationbuffer2(single_index) = population

    if (stop_flag) then
        close(var_filechannel)
        return
    end if

else if (population == 0) then
    return

else if (population < buffer2_size) then

    do k = 1, population
    
        do n = 1, Nindistinguishables
    
            BOND_LABELLING_DATA = INDISTINGUISHABLES(n,:)
        
            do j = 1, Natoms
                 temp_coords(:,j) = &
                     coordsbuffer2(:,BOND_LABELLING_DATA(j),k,single_index)
            end do
        
            !Calculate the RMSD between this frame and the incoming frame
            call rmsd_dp(Natoms,temp_coords,var_coords,&
                         1,temp_U,x_center,y_center,&
                         temp_rmsd)                               !,.false.,g)

            call getCoulombMatrix(Natoms,temp_coords,charges,temp_CM)
        
            call getCoulombMatrixDiff(Natoms,temp_CM,var_CM,temp_CMdiff)

!           do i = 1, 3
!               temp_U(i,:) = 0.0d0
!               x_center(i) = 0.0d0
!               y_center(i) = 0.0d0
!           end do
!           temp_rmsd = rand() * threshold_rmsd * 4
!           temp_CMdiff = 0.0d0
        
            !If in the "accept best" method
            !and the RMSD is low enough:
            if ((temp_rmsd < threshold_rmsd).and.&
                (temp_rmsd >= inner_threshold)) then
        
                if (accept_first) iostate = 1
        
                if (accept_worst) then
                    do i = 1, Ninterpolation+1
                        if (temp_rmsd > RMSDbuffer1(i)) exit
                    end do
                else
!                   do i = 1, Ninterpolation+1
!                       if (temp_rmsd < RMSDbuffer1(i)) exit
!                   end do
                    do i = 1, Ninterpolation+1
                        if (temp_CMdiff < CMdiffbuffer1(i)) exit
                    end do
                end if
        
                if (Ninterpolation < Ninterpolation_max) &
                    Ninterpolation = Ninterpolation + 1
        
                do j = Ninterpolation, i+1, -1
                    valsbuffer1(:,j) = valsbuffer1(:,j-1)
                    coordsbuffer1(:,:,j) = coordsbuffer1(:,:,j-1)
                    gradientbuffer1(:,:,j) = gradientbuffer1(:,:,j-1)
                    Ubuffer1(:,:,j) = Ubuffer1(:,:,j-1)
                    RMSDbuffer1(j) = RMSDbuffer1(j-1)
                    CMdiffbuffer1(j) = CMdiffbuffer1(j-1)
                    inputCLS(:,j) = inputCLS(:,j-1)
                end do
        
                do j = 1, Natoms
                    coordsbuffer1(:,j,i) = &
                        coordsbuffer2(:,BOND_LABELLING_DATA(j),k,single_index)
                    gradientbuffer1(:,j,i) = &
                        gradientbuffer2(:,BOND_LABELLING_DATA(j),k,single_index)
                end do
        
                valsbuffer1(:,i) = valsbuffer2(:,k,single_index)
                Ubuffer1(:,:,i) = temp_U
                RMSDbuffer1(i) = temp_rmsd
                CMdiffbuffer1(i) = temp_CMdiff
        
                temp_coords = coordsbuffer1(:,:,i)
        
                do j = 1, 3
                    temp_coords(j,:) = &
                    temp_coords(j,:) - x_center(j)
                end do
        
                temp_coords(:,:) = matmul(&
                        Ubuffer1(:,:,i),&
                        temp_coords(:,:))
        
                do j = 1, 3
                    temp_coords(j,:) = &
                    temp_coords(j,:) + y_center(j)
                end do
        
                inputCLS(1:Ncoords,i) =&
                           reshape(temp_coords(:,:) - &
                           var_coords,(/Ncoords/))
            end if
    
        end do
    
    end do

    number_of_cells = number_of_cells + 1
else
    write(var_filename,FMT=var_multipleFMT&
          (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
            var_index * multiplier(:,Norder+1)
    
    subcell = gridpath3//trim(var_filename)

    number_of_cells = number_of_cells + 1
    population = 0

    !Open the file corresponding to the cell
    if (unreadable_flag) then
        open(var_filechannel,action="read",form="unformatted",&
             file=trim(subcell))
    else
        open(var_filechannel,action="read",&
             file=trim(subcell))
    end if
end if

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

        !The next line is the coordinates
        read(var_filechannel) &
               ((current_coords(i,j),i=1,3),j=1,Natoms)

        !In unformatted files, the last line is the gradient
        read(var_filechannel) &
                ((current_gradient(i,j),i=1,3),j=1,Natoms)
    else

        !In formatted files, everything (variables, coordinates, and gradient)
        !are stored in one line; FMT1 reads the variables
        read(var_filechannel,FMT=FMT1,advance="no",iostat=iostate) &
                (current_vals(i),i=1,Nvar)

        !If there are no more lines, stop; the population of the cell should be
        !one less the number of times that this portion of the loop was called
        if (iostate /= 0) exit

        population = population + 1

        !In formatted files, FMT3 reads the coordinates
        read(var_filechannel,FMT=FMT7,advance="no") &
               ((current_coords(i,j),i=1,3),j=1,Natoms)

        !In formatted files, FMT3 reads the gradient as well
        read(var_filechannel,FMT=FMT3) &
                ((current_gradient(i,j),i=1,3),j=1,Natoms)

    end if

    do n = 1, Nindistinguishables

        BOND_LABELLING_DATA = INDISTINGUISHABLES(n,:)
    
        do j = 1, Natoms
             temp_coords(:,j) = &
                 current_coords(:,BOND_LABELLING_DATA(j))
        end do
    
        !Calculate the RMSD between this frame and the incoming frame
        call rmsd_dp(Natoms,temp_coords,var_coords,&
                     1,temp_U,x_center,y_center,&
                     temp_rmsd)                               !,.false.,g)

        call getCoulombMatrix(Natoms,temp_coords,charges,temp_CM)
    
        call getCoulombMatrixDiff(Natoms,temp_CM,var_CM,temp_CMdiff)

!       do i = 1, 3
!           temp_U(i,:) = 0.0d0
!           x_center(i) = 0.0d0
!           y_center(i) = 0.0d0
!       end do
!       temp_rmsd = rand() * threshold_rmsd * 4
!       temp_CMdiff = 0.0d0
    
        !If in the "accept best" method
        !and the RMSD is low enough:
        if ((temp_rmsd < threshold_rmsd).and.&
            (temp_rmsd >= inner_threshold)) then
    
            if (accept_first) iostate = 1
    
            if (accept_worst) then
                do i = 1, Ninterpolation+1
                    if (temp_rmsd > RMSDbuffer1(i)) exit
                end do
            else
!               do i = 1, Ninterpolation+1
!                   if (temp_rmsd < RMSDbuffer1(i)) exit
!               end do
                do i = 1, Ninterpolation+1
                    if (temp_CMdiff < CMdiffbuffer1(i)) exit
                end do
            end if
    
            if (Ninterpolation < Ninterpolation_max) &
                Ninterpolation = Ninterpolation + 1
    
            do j = Ninterpolation, i+1, -1
                valsbuffer1(:,j) = valsbuffer1(:,j-1)
                coordsbuffer1(:,:,j) = coordsbuffer1(:,:,j-1)
                gradientbuffer1(:,:,j) = gradientbuffer1(:,:,j-1)
                Ubuffer1(:,:,j) = Ubuffer1(:,:,j-1)
                RMSDbuffer1(j) = RMSDbuffer1(j-1)
                CMdiffbuffer1(j) = CMdiffbuffer1(j-1)
                inputCLS(:,j) = inputCLS(:,j-1)
            end do
    
            do j = 1, Natoms
                coordsbuffer1(:,j,i) = &
                    current_coords(:,BOND_LABELLING_DATA(j))
                gradientbuffer1(:,j,i) = &
                    current_gradient(:,BOND_LABELLING_DATA(j))
            end do
    
            valsbuffer1(:,i) = current_vals
            Ubuffer1(:,:,i) = temp_U
            RMSDbuffer1(i) = temp_rmsd
            CMdiffbuffer1(i) = temp_CMdiff
    
            temp_coords = coordsbuffer1(:,:,i)
    
            do j = 1, 3
                temp_coords(j,:) = &
                temp_coords(j,:) - x_center(j)
            end do
    
            temp_coords(:,:) = matmul(&
                    Ubuffer1(:,:,i),&
                    temp_coords(:,:))
    
            do j = 1, 3
                temp_coords(j,:) = &
                temp_coords(j,:) + y_center(j)
            end do
    
            inputCLS(1:Ncoords,i) =&
                       reshape(temp_coords(:,:) - &
                       var_coords,(/Ncoords/))
        end if

    end do

end do
close(var_filechannel)

end subroutine getRMSD_1_PCM




subroutine getRMSD_2(var_index,population)
use ls_rmsd_original
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

!Incremental integers and iostate checking
integer :: i,j,k,iostate
integer :: endpoint

write(var_filename,FMT=var_multipleFMT&
      (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
        var_index * multiplier(:,Norder+1)

!Construct the subcell filename
subcell = gridpath3//trim(var_filename)

!See whether this subcell exists
inquire(file=trim(subcell),exist=subcell_existence)

if (.not. subcell_existence) then
        population = 0
        return
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

!Initialize a variable
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

                !In formatted files, FMT3 reads the coordinates
                read(var_filechannel,FMT=FMT7,advance="no") &
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






subroutine addState_new(vals,coords,gradient)
use PHYSICS
use PARAMETERS
implicit none

!Inputs for the file
real(dp),dimension(Nvar) :: vals
real(dp),dimension(3,Natoms) :: coords
real(dp),dimension(3,Natoms) :: gradient

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
                                   var_index(:,Norder+1))

                exit
        
        else
                cycle
        end if

end do

end subroutine addState_new


subroutine frameAddition(vals,coords,gradient,var_index,nolabel_flag)
use PARAMETERS
use PHYSICS
implicit none

!Inputs for file writing
real(dp),dimension(Nvar),intent(in) :: vals
real(dp),dimension(3,Natoms),intent(in) :: coords,gradient
integer,dimension(Nvar),intent(in) :: var_index

!Character strings used to identify the file
character(50) :: var_filename

logical,optional :: nolabel_flag

!Incremental integers
integer :: i,j,k

write(var_filename,FMT=var_multipleFMT&
      (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
        var_index * multiplier(:,Norder+1)

!Pretty self-explanatory
if (unreadable_flag) then
        open(filechannel1,file=gridpath2//trim(var_filename),position="append",form="unformatted")
        if ((force_NoLabels).or.(present(nolabel_flag).and.(nolabel_flag))) then
                write(filechannel1) (vals(j),j=1,Nvar)
!               write(filechannel1) Ntraj
                write(filechannel1) ((coords(i,j),i=1,3),j=1,Natoms)
                write(filechannel1) ((gradient(i,j),i=1,3),j=1,Natoms)
        else
                write(filechannel1) (vals(j),j=1,Nvar)
!               write(filechannel1) Ntraj
!               write(filechannel1) ((coords(i,BOND_LABELLING_DATA(j)),i=1,3),j=1,Natoms)
!               write(filechannel1) ((gradient(i,BOND_LABELLING_DATA(j)),i=1,3),j=1,Natoms)
                write(filechannel1) ((coords(i,j),i=1,3),j=1,Natoms)
                write(filechannel1) ((gradient(i,j),i=1,3),j=1,Natoms)
        end if
else
        open(filechannel1,file=gridpath2//trim(var_filename),position="append")
        if ((force_NoLabels).or.(present(nolabel_flag).and.(nolabel_flag))) then
                write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
                write(filechannel1,FMT=FMT3,advance="no") ((coords(i,j),i=1,3),j=1,Natoms)
                write(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
        else
                write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
!               write(filechannel1,FMT=FMT3,advance="no") ((coords(i,BOND_LABELLING_DATA(j)),i=1,3),j=1,Natoms)
!               write(filechannel1,FMT=FMT3) ((gradient(i,BOND_LABELLING_DATA(j)),i=1,3),j=1,Natoms)
                write(filechannel1,FMT=FMT3,advance="no") ((coords(i,j),i=1,3),j=1,Natoms)
                write(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
        end if
end if
close(filechannel1)

end subroutine frameAddition



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

Norder = Norder + 1

do i = 1, frames

        do j = 1, Nvar
                var_index(j) = int(valsbuffer1(j,i) * divisor(j,Norder+1))
        end do

        call frameAddition(valsbuffer1(:,i),coordsbuffer1(:,:,i),&
                           gradientbuffer1(:,:,i),var_index,nolabel_flag)
end do

Norder = Norder - 1

end subroutine divyUp



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               FUNCTION GETRMSD_DP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  actually pretty self-explanatory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getRMSD_dp(filechannel_thread,subcell,coords,population,min_rmsd,gradient,U)

use ls_rmsd_original
use ANALYSIS
use PARAMETERS
implicit none
integer,intent(out) :: population
real(dp),intent(inout), dimension(3,Natoms) :: gradient
!real(dp),dimension(3,Natoms) :: candidate_gradient
real(dp),intent(inout) :: min_rmsd
integer,intent(in) :: filechannel_thread
integer :: OMP_GET_THREAD_NUM
integer :: i,j,k,iostate
character(*),intent(in) :: subcell
real(dp),intent(in), dimension(3,Natoms) :: coords
real(dp) :: candidate_min_rmsd
real(dp), dimension(3) :: x_center, y_center
real(dp), allocatable :: g(:,:)
real(dp),dimension(3,3),intent(inout) :: U
real(dp),dimension(3,3) :: candidate_U
real(dp), dimension(3,Natoms) :: coords2

!Open the file corresponding to the cell
if (unreadable_flag) then
        open(filechannel_thread,action="read",file=trim(subcell)//".dat",form="unformatted")
else
        open(filechannel_thread,action="read",file=trim(subcell)//".dat")
end if

population = 0
do
        !Read the candidate frame
        if (unreadable_flag) then
                read(filechannel_thread,iostat=iostate)
	        if (iostate /= 0) exit
                read(filechannel_thread) ((coords2(i,j),i=1,3),j=1,Natoms)
        else
                read(filechannel_thread,FMT=FMT7,advance="no",iostat=iostate) ((coords2(i,j),i=1,3),j=1,Natoms)
	        if (iostate /= 0) exit
        end if

        call rmsd_dp(Natoms,coords2,coords,1,candidate_U,x_center,y_center,candidate_min_rmsd)!,.false.,g)

        !Increment the number of frames visited
        population = population + 1

        !If in the "accept worst" method:
        if (accept_worst) then

        !If it is too low, reject the candidate frame
        if (candidate_min_rmsd < min_rmsd) then
                if (unreadable_flag) then
                        read(filechannel_thread) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
                else
                        read(filechannel_thread,FMT=FMT3) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
                end if

        !Otherwise if it is low enough, accept the candidate frame (and save its rotation matrix!)
        else if (candidate_min_rmsd < threshold_rmsd) then
                min_rmsd = candidate_min_rmsd
                U = candidate_U

                if (unreadable_flag) then
                        read(filechannel_thread) ((gradient(j,k),j=1,3),k=1,Natoms)
                else
                        read(filechannel_thread,FMT=FMT3) ((gradient(j,k),j=1,3),k=1,Natoms)
                end if

                !In special cases, we exit immediately afterwards
                if (accept_first) exit

        !If it is too high, also reject the candidate frame
        else
                if (unreadable_flag) then
                        read(filechannel_thread) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
                else
                        read(filechannel_thread,FMT=FMT3) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
                end if
        end if

        !If in the "accept best" method:
        else

        !If it is low enough, accept the candidate frame (and save its rotation matrix!)
        if (candidate_min_rmsd < min_rmsd) then
                min_rmsd = candidate_min_rmsd
                U = candidate_U

                if (unreadable_flag) then
                        read(filechannel_thread) ((gradient(j,k),j=1,3),k=1,Natoms)
                else
                        read(filechannel_thread,FMT=FMT3) ((gradient(j,k),j=1,3),k=1,Natoms)
                end if

                !In special cases, we exit immediately afterwards
                if (accept_first) exit

        !Otherwise, just don't record the gradient
        else
                if (unreadable_flag) then
                        read(filechannel_thread) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
                else
                        read(filechannel_thread,FMT=FMT3) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
                end if
        end if

        end if

end do
close(filechannel_thread)

end subroutine getRMSD_dp



subroutine setAllocations()
    use PARAMETERS
    use ANALYSIS
    implicit none

    allocate(valsbuffer1(Nvar,buffer1_size),&
             coordsbuffer1(3,Natoms,buffer1_size),&
             gradientbuffer1(3,Natoms,buffer1_size),&
             Ubuffer1(3,3,buffer1_size),&
             RMSDbuffer1(buffer1_size))

    allocate(CMdiffbuffer1(buffer1_size),&
             Ntrajbuffer1(buffer1_size),&
             varindexbuffer1(buffer1_size))

    allocate(acceptable_frame_mask(buffer1_size),&
             inputCLS(Ncoords+buffer1_size,buffer1_size))

    if (memory_flag) then
        allocate(vals_hash(single_index_max,Nvar))

        allocate(populationbuffer2(&
                    single_index_max),&
                 valsbuffer2(Nvar,&
                    buffer2_size,&
                    single_index_max),&
                 coordsbuffer2(3,Natoms,&
                    buffer2_size,&
                    single_index_max),&
                 gradientbuffer2(3,Natoms,&
                    buffer2_size,&
                    single_index_max))

        allocate(temppopulationbuffer2(&
                    single_index_max),&
                 tempvalsbuffer2(Nvar,&
                    buffer2_size,&
                    single_index_max),&
                 tempcoordsbuffer2(3,Natoms,&
                    buffer2_size,&
                    single_index_max),&
                 tempgradientbuffer2(3,Natoms,&
                    buffer2_size,&
                    single_index_max))
    end if

    return

end subroutine setAllocations

subroutine unsetAllocations()
    use PARAMETERS
    use ANALYSIS
    implicit none

    deallocate(valsbuffer1,&
               coordsbuffer1,&
               gradientbuffer1,&
               Ubuffer1,&
               RMSDbuffer1)

    deallocate(CMdiffbuffer1,&
               Ntrajbuffer1,&
               varindexbuffer1)

    deallocate(acceptable_frame_mask,&
               inputCLS)

    if (memory_flag) then
        deallocate(vals_hash)

        deallocate(populationbuffer2,&
                   valsbuffer2,&
                   coordsbuffer2,&
                   gradientbuffer2)

        deallocate(temppopulationbuffer2,&
                   tempvalsbuffer2,&
                   tempcoordsbuffer2,&
                   tempgradientbuffer2)
    end if

    return

end subroutine unsetAllocations

subroutine shiftMemory(delta_var_index)
    use PARAMETERS
    use ANALYSIS
    use FUNCTIONS
    implicit none
    integer,dimension(Nvar),intent(in) :: delta_var_index
    integer :: population
    integer :: single_index,xflattened
    integer,dimension(Nvar) :: xexpanded

    if (all(delta_var_index == 0)) return

    temppopulationbuffer2 = populationbuffer2
    tempvalsbuffer2 = valsbuffer2
    tempcoordsbuffer2 = coordsbuffer2
    tempgradientbuffer2 = gradientbuffer2

    populationbuffer2 = -1

    do single_index = 1, single_index_max

        population = temppopulationbuffer2(single_index)
        if (population == -1) cycle

        call getExpanded(Nvar,single_index,xexpanded)
        call getFlattened(Nvar,xexpanded-delta_var_index,xflattened)

        if (xflattened <= single_index_max) then
            populationbuffer2(xflattened) = population
            valsbuffer2(:,1:population,xflattened) = &
                tempvalsbuffer2(:,1:population,single_index)
            coordsbuffer2(:,:,1:population,xflattened) = &
                tempcoordsbuffer2(:,:,1:population,single_index)
            gradientbuffer2(:,:,1:population,xflattened) = &
                tempgradientbuffer2(:,:,1:population,single_index)
        end if

    end do

!   populationbuffer2 = -1

    return

end subroutine shiftMemory


end module interactMultipleGrids
