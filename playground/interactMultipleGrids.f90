!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!In the latter case, we need to keep track of how large
!the buffer gets, and increase it if it gets too large

integer :: buffer1_size

!It is useful to estimate how many cells we may be
!looking at at one time

integer,dimension(Norder_max+1) :: subcellsearch_max1 = (/ 1, 5 /)
integer,dimension(Norder_max+1) :: subcellsearch_max2 = (/ 2, 2 /)

integer,dimension(Norder_max+1) :: subcellsearch_max = (/ 1, 5 /)
integer,dimension(Norder_max+1) :: local_frame_count

!Assuming that each subcellsearch_max < 5
!This only approximates the number of cells in the outer
!shell of a Nvar-dimensional cube (not diamond, like
!used in the program)
integer,parameter :: number_of_cells_max = &
        2*((1 + 2 * (5))**(Nvar) - (1 + 2 * (5 - 1))**(Nvar))
integer,parameter :: number_of_frames_max = &
        number_of_cells_max * var_overcrowd(2)

!Store the index of the best candidate in approximation_index

integer :: number_of_cells
integer :: approximation_index(number_of_cells_max)

!Other global variables to clean things up

integer :: Totalnumber_of_frames
integer :: Norder
integer :: var_filechannel
real(dp),dimension(3,Natoms) :: var_coords
real(dp) :: candidate_rmsd




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
real(dp), dimension(3,Natoms), intent(out) :: gradient
real(dp), intent(inout) :: min_rmsd
integer, dimension(1+Ngrid_total),intent(in) :: filechannels
character(Ngrid_text_length) :: Ngrid_text
character(5) :: variable_length_text

real(dp),dimension(Nvar) :: var_cell
integer,dimension(Nvar,Norder_max+1) :: var_index
integer,dimension(Nvar,Norder_max+1) :: var_index_diff

var_filechannel = filechannels(1)
var_coords = coords

!We start off with zero frames having been checked
local_frame_count = 0
Totalnumber_of_frames = 0
number_of_cells = 0
neighbor_check = 0
approximation_index = 0
stop_flag = .false.

candidate_rmsd = min_rmsd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 SUBCELL TARGETING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Retrive the index of each variable with respect to the grid
!and the real number (rounded) that represents that index
do i = 1, Nvar
       !Repeat this for however many orders of cells deep
       !we have been instructed to go
       do j = 1, Norder_max + 1
              var_index(i,j) = int(vals(i) * divisor(i,j))
       end do
end do

!print *, ""
!print *, ""
!print *, ""
!print *, "vals:", vals

!print *, "value: ", vals
!print *, "var_index: ", var_index(:,1)
!print *, "var_round: ", var_index(:,1) * multiplier(:,1)
!print *, "var_index: ", var_index(:,2)
!print *, "var_round: ", var_index(:,2) * multiplier(:,2)

!Now, we start iterating over the grids
do k = 1, Ngrid_total

write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") k
gridpath3 = gridpath0//Ngrid_text//"/grid/"

do l = 1, Norder_max+1

        Norder = Norder_order(l)
        
        !Read the frames in the cells and process their RMSDs
        call getRMSD_1(var_index(:,Norder+1),population)
        local_frame_count(Norder+1) = population
        
        !If the cell is populated...
        if (population > 0) then
        
                !However many frames are in the subcell we
                !increment to number_of_frames
                Totalnumber_of_frames = &
                Totalnumber_of_frames + population
                
                !This is good enough, no need to look at more frames
                stop_flag = .true.
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
        
        if (stop_flag) exit

end do

if (testtraj_flag) write(filechannels(1+k),FMT=FMT6) candidate_rmsd

end do

number_of_frames = Totalnumber_of_frames
order = Norder
neighbor_check = number_of_cells

!print *, ""
!print *, "step: ", steps
!print *, "Norder: ", Norder
!print *, "value: ", vals
!print *, "coords: ", coords
!print *, "var_index: ", var_index(:,1)
!print *, "var_index: ", var_index(:,2)
!print *, "approx_index: ", approximation_index(1:number_of_cells)
!print *, "RMSDbuffer1: ", RMSDbuffer1(1:Totalnumber_of_frames)
!print *, "candidatermsd: ", candidate_rmsd

if (number_of_cells == 0) then
        return
else
        Norder_total(Norder+1) = Norder_total(Norder+1) + 1

!if (Norder == 1) then
!print *, ""
!print *, "step: ", steps
!print *, "value: ", vals
!print *, "coords: ", coords
!print *, "var_index: ", var_index(:,1)
!print *, "var_index: ", var_index(:,2)
!print *, "approx_index: ", approximation_index(1:number_of_cells)
!print *, "RMSDbuffer1: ", RMSDbuffer1(1:Totalnumber_of_frames)
!print *, "candidatermsd: ", candidate_rmsd
!end if

        chosen_index = maxval(approximation_index(1:number_of_cells),DIM=1)
end if

if (chosen_index > 0) then
        min_rmsd = RMSDbuffer1(chosen_index)

        gradient = matmul(Ubuffer1(:,:,chosen_index),&
                          gradientbuffer1(:,:,chosen_index))
else
        min_rmsd = default_rmsd
end if

return

end subroutine checkState_new



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
       stop_flag = .false.
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
        
                        !The next line is the coordinates
                        read(var_filechannel) &
                               ((coordsbuffer1(i,j,k),i=1,3),j=1,Natoms)
        
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
                               ((coordsbuffer1(i,j,k),i=1,3),j=1,Natoms)
        
                        !In formatted files, FMT3 reads the gradient as well
                        read(var_filechannel,FMT=FMT3) &
                                ((gradientbuffer1(i,j,k),i=1,3),j=1,Natoms)
                end if
        
                !Calculate the RMSD between this frame and the incoming frame
                call rmsd_dp(Natoms,coordsbuffer1(:,:,k),var_coords,&
                             1,Ubuffer1(:,:,k),x_center,y_center,&
                             RMSDbuffer1(k))                               !,.false.,g)
        
                !If in the "accept worst" method:
                if (accept_worst) then
        
                        !If the RMSD is too low, reject the candidate frame for it is "too good"
                        if (RMSDbuffer1(k) < candidate_rmsd) then
                
                        !And if the RMSD if lower than the threshold, designate this frame as the "worst"
                        else if (RMSDbuffer1(k) < threshold_rmsd) then
                
                                approximation_index(number_of_cells) = k
                                candidate_rmsd = RMSDbuffer1(k)
                
                                !In special cases, we exit immediately afterwards
                                if (accept_first) exit
                
                        !If the RMSD is too high, also reject the candidate frame
                        else
                
                        end if
        
                !If in the "accept best" method:
                else

                        !If the RMSD is low enough, designate this frame as the "best"
                        if (RMSDbuffer1(k) < candidate_rmsd) then
                        
                                approximation_index(number_of_cells) = k
                                candidate_rmsd = RMSDbuffer1(k)
                
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
                 temp_RMSDbuffer1(buffer1_size))
        
        !Store the buffer in the temporary buffer
        temp_valsbuffer1 = valsbuffer1
        temp_coordsbuffer1 = coordsbuffer1
        temp_gradientbuffer1 = gradientbuffer1
        temp_Ubuffer1 = Ubuffer1
        temp_RMSDbuffer1 = RMSDbuffer1
        
        !For now, we simply double the buffer size each time it overfills
        deallocate(valsbuffer1,coordsbuffer1,gradientbuffer1,Ubuffer1,RMSDbuffer1)
        allocate(valsbuffer1(Nvar,buffer1_size*2),&
                 coordsbuffer1(3,Natoms,buffer1_size*2),&
                 gradientbuffer1(3,Natoms,buffer1_size*2),&
                 Ubuffer1(3,3,buffer1_size*2),&
                 RMSDbuffer1(buffer1_size*2))
        
        !And reload all the frames back into the buffer
        valsbuffer1(:,1:buffer1_size) = temp_valsbuffer1
        coordsbuffer1(:,:,1:buffer1_size) = temp_coordsbuffer1
        gradientbuffer1(:,:,1:buffer1_size) = temp_gradientbuffer1
        Ubuffer1(:,:,1:buffer1_size) = temp_Ubuffer1
        RMSDbuffer1(1:buffer1_size) = temp_RMSDbuffer1

        !Destroy the temporary buffer
        deallocate(temp_valsbuffer1,temp_coordsbuffer1,temp_gradientbuffer1,&
                   temp_Ubuffer1,temp_RMSDbuffer1)

        !Permanently increase the buffer size so this will not
        !have to happen next time
        buffer1_size = buffer1_size*2
end do
close(var_filechannel)

end subroutine getRMSD_1




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
integer :: i,j

!Retrive the index of each variable with respect to the grid
!and the real number (rounded) that represents that index
do i = 1, Nvar
       !Repeat this for however many orders of cells deep
       !we have been instructed to go
       do j = 1, Norder_max + 1
              var_index(i,j) = int(vals(i) * divisor(i,j))
       end do
end do

do Norder = 0, Norder_max

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

                if (population == 0) Nfile = Nfile + 1
        
                call frameAddition(vals,coords,gradient,&
                                   var_index(:,Norder+1))

                exit
        
        else
                cycle
        end if

end do

end subroutine addState_new


subroutine frameAddition(vals,coords,gradient,var_index)
use PARAMETERS
use PHYSICS
implicit none

!Inputs for file writing
real(dp),dimension(Nvar),intent(in) :: vals
real(dp),dimension(3,Natoms),intent(in) :: coords,gradient
integer,dimension(Nvar),intent(in) :: var_index

!Character strings used to identify the file
character(50) :: var_filename

!Incremental integers
integer :: i,j,k

write(var_filename,FMT=var_multipleFMT&
      (1+Norder*multipleFMT_length:(Norder+1)*multipleFMT_length) )&
        var_index * multiplier(:,Norder+1)

!Pretty self-explanatory
if (unreadable_flag) then
        open(filechannel1,file=gridpath3//trim(var_filename),position="append",form="unformatted")
        if (force_NoLabels) then
                write(filechannel1) (vals(j),j=1,Nvar)
                write(filechannel1) ((coords(i,j),i=1,3),j=1,Natoms)
                write(filechannel1) ((gradient(i,j),i=1,3),j=1,Natoms)
        else
                write(filechannel1) (vals(j),j=1,Nvar)
                write(filechannel1) ((coords(i,BOND_LABELLING_DATA(j)),i=1,3),j=1,Natoms)
                write(filechannel1) ((gradient(i,BOND_LABELLING_DATA(j)),i=1,3),j=1,Natoms)
        end if
else
        open(filechannel1,file=gridpath3//trim(var_filename),position="append")
        if (force_NoLabels) then
                write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
                write(filechannel1,FMT=FMT3,advance="no") ((coords(i,j),i=1,3),j=1,Natoms)
                write(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
        else
                write(filechannel1,FMT=FMT1,advance="no") (vals(j),j=1,Nvar)
                write(filechannel1,FMT=FMT3,advance="no") ((coords(i,BOND_LABELLING_DATA(j)),i=1,3),j=1,Natoms)
                write(filechannel1,FMT=FMT3) ((gradient(i,BOND_LABELLING_DATA(j)),i=1,3),j=1,Natoms)
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

integer :: i, j

Norder = Norder + 1

do i = 1, frames

        do j = 1, Nvar
                var_index(j) = int(valsbuffer1(j,i) * divisor(j,Norder+1))
        end do

        call frameAddition(valsbuffer1(:,i),coordsbuffer1(:,:,i),&
                           gradientbuffer1(:,:,i),var_index)
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
real(dp),dimension(3,Natoms) :: candidate_gradient
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



end module interactMultipleGrids
