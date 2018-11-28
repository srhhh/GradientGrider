!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       MODULE
!               interactSingleGrid
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
!               FILECHANNEL1                    OPEN, WRITE, CLOSE
!               FILECHANNEL1                    OPEN, READ, CLOSE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINES                     ARGUMENTS               KIND
!
!		checkState
!
!		getNeighbors
!
!		getRMSD_dp
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       CALLS                           MODULE
!
!               getVar3                         VARIABLES
!               getVar4                         VARIABLES
!
!               getNeighbors                    interactMultipleGrids
!               getRMSD_dp                      interactMultipleGrids
!
!               rmsd                            ls_rmsd_original
!               matmul                          PHYSICS
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
implicit none


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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine checkState(vals,coords,gradient,min_rmsd,filechannels,&
                      number_of_frames,order,neighbor_check)
use ANALYSIS
use VARIABLES
use PARAMETERS
implicit none
integer :: i,j,k
integer,intent(out),optional :: order,number_of_frames,neighbor_check
integer :: var1_index0,var2_index0,var1_index1,var2_index1
integer :: index1_1,index1_2,index2_1,index2_2
integer :: population
logical :: stop_flag
real :: var1,var2,var1_cell,var2_cell,var1_round,var2_round
real :: var1_round0,var2_round0,var1_round1,var2_round1,var1_round2,var2_round2,var1_round3,var2_round3
real(dp), dimension(Nvar), intent(in) :: vals
real(dp), dimension(3,Natoms), intent(in) :: coords
real(dp), dimension(3,Natoms), intent(inout) :: gradient
real(dp), dimension(3,3) :: U, old_U
real(dp), intent(inout) :: min_rmsd
real(dp)  :: old_min_rmsd
real(dp), dimension(3,Natoms) :: old_gradient
logical :: subcell_existence
character(100) :: subcell
character(FMTlength) ::  var1_filename0, var2_filename0, var1_filename1, var2_filename1
character(FMTlength) :: var1_filename,var2_filename
integer, dimension(Ngrid_total),intent(in) :: filechannels
integer :: subcell0search_max,subcell1search_max
character(Ngrid_text_length) :: Ngrid_text
character(5) :: variable_length_text

!In order to decrease the number of frames checked
!If we need to look at adjacent cells for a frame
!(for instance, if the target cell is empty)
!Then we stop looking after we are subcellsearch_max cells away
subcell0search_max = 1
subcell1search_max = 1

!We start off with zero frames having been checked
number_of_frames = 0
neighbor_check = 0

old_min_rmsd = min_rmsd
old_gradient = gradient
old_U = 0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 SUBCELL TARGETING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!var_cell keeps track of the child subcell portion of the decimal
!and throws away the parent subcell portion of the decimal
var1_cell = vals(1)
var2_cell = vals(2)

!var_index keeps track of the index of the child subcell
!inside of its parent subcell
var1_index0 = int(var1_cell * divisor1_0)
var2_index0 = int(var2_cell * divisor2_0)

!var_round is the variable rounded to the subcell length multiple
!This may have more or even less digits than required
var1_round0 = multiplier1_0 * var1_index0
var2_round0 = multiplier2_0 * var2_index0

!Repeat for the child subcell portion
var1_cell = var1_cell - var1_round0
var2_cell = var2_cell - var2_round0

var1_index1 = int(var1_cell * divisor1_1)
var2_index1 = int(var2_cell * divisor2_1)

var1_round1 = multiplier1_1 * var1_index1
var2_round1 = multiplier2_1 * var2_index1

var1_round = var1_round0 + var1_round1
var2_round = var2_round0 + var2_round1

!The child-level subcells will use var_round
write(var1_filename1,FMT=FMTorder1) var1_round
write(var2_filename1,FMT=FMTorder1) var2_round

!The parent-level subcells will use var_round0
write(var1_filename0,FMT=FMTorder0) var1_round0
write(var2_filename0,FMT=FMTorder0) var2_round0

!Now, we start iterating over the grids
do Ngrid = 1, Ngrid_total

order = 1
neighbor_check = 0

write(variable_length_text,FMT=FMT5_variable) Ngrid_text_length
write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
gridpath2 = gridpath0//Ngrid_text//"/grid/"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 ORDER 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Construct the subcell filename
subcell = gridpath2//trim(adjustl(var1_filename1))//"_"//trim(adjustl(var2_filename1))

!See whether this subcell exists
inquire(file=trim(subcell)//".dat",exist=subcell_existence)
stop_flag = .false.

!For bug-testing
if (.false.) then
print *, ""
print *, ""
print *, "inquiring for: ", subcell
print *, "exists? ", subcell_existence
print *, ""
end if

!If the subcell exists, we will look in here for a closeby frame
if (subcell_existence) then
        !getRMSD reads off the coordinates and calculated the RMSD with
        !ls_rmsd module
	call getRMSD_dp(subcell,coords,population,min_rmsd,gradient,U)

	!However many frames are in the subcell we increment to number_of_frames
	number_of_frames = number_of_frames + population

	if (min_rmsd /= old_min_rmsd) then
		old_min_rmsd = min_rmsd
		old_gradient = gradient
		old_U = U
	end if

	stop_flag = .true.
end if

!If there are no frames in this cell, we can look at one of its
!neighbors (better than having to look at the parent)
if ((.not.(subcell_existence)).or.(force_Neighbors)) then
        neighbor_check = 1

        !Integer i keeps track of how far away from the original
        !subcell we are; we look at cells on the 'diamond' surrounding
        !the original subcell
        do i = 1, subcell1search_max
                neighbor_check = i

                index1_1 = var1_index1 + i
                index1_2 = var1_index1 - i

                index2_1 = var2_index1 + i
                index2_2 = var2_index1 - i
                
                if (index1_1 < scaling1_0) then
                
                        !Make the name of the subcell
                        var1_round = var1_round0 + index1_1 * multiplier1_1
                        var2_round = var2_round0 + var2_round1
                        write(var1_filename,FMT=FMTorder1) var1_round
                        write(var2_filename,FMT=FMTorder1) var2_round
                        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
                
                        !Check if it exists
                        inquire(file=trim(subcell)//".dat",exist=subcell_existence)
                        if (subcell_existence) then
                                call getRMSD_dp(subcell,coords,population,min_rmsd,gradient,U)
                        end if
                
                        number_of_frames = number_of_frames + population
                
                end if
                !Rinse and repeat
                if (index1_2 .ge. 0) then
                
                        var1_round = var1_round0 + index1_2 * multiplier1_1
                        var2_round = var2_round0 + var2_round1
                        write(var1_filename,FMT=FMTorder1) var1_round
                        write(var2_filename,FMT=FMTorder1) var2_round
                        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
                
                        inquire(file=trim(subcell)//".dat",exist=subcell_existence)
                        if (subcell_existence) then
                                call getRMSD_dp(subcell,coords,population,min_rmsd,gradient,U)
                        end if
                
                        number_of_frames = number_of_frames + population
                end if
                !Rinse and repeat
                if (index2_1 < scaling2_0) then
                
                        var1_round = var1_round0 + var1_round1
                        var2_round = var2_round0 + index2_1 * multiplier2_1
                        write(var1_filename,FMT=FMTorder1) var1_round
                        write(var2_filename,FMT=FMTorder1) var2_round
                        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
                
                        inquire(file=trim(subcell)//".dat",exist=subcell_existence)
                        if (subcell_existence) then
                                call getRMSD_dp(subcell,coords,population,min_rmsd,gradient,U)
                        end if
                
                        number_of_frames = number_of_frames + population
                end if
                !Rinse and repeat
                if (index2_2 .ge. 0) then
                
                        var1_round = var1_round0 + var1_round1
                        var2_round = var2_round0 + index2_2 * multiplier2_1
                        write(var1_filename,FMT=FMTorder1) var1_round
                        write(var2_filename,FMT=FMTorder1) var2_round
                        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
                
                        inquire(file=trim(subcell)//".dat",exist=subcell_existence)
                        if (subcell_existence) then
                                call getRMSD_dp(subcell,coords,population,min_rmsd,gradient,U)
                        end if
                
                        number_of_frames = number_of_frames + population
                end if

                if (min_rmsd /= old_min_rmsd) then
                        stop_flag = .true.
			old_min_rmsd = min_rmsd
			old_gradient = gradient
			old_U = U
                end if
        
                !And integer j keeps track of where on the circumfrence
                !Of the diamond surrounding the subcell we are at
                do j = 1, i-1
                
                        !Yes, there are redundancies for j = 0 and j = i
                        !This is the first thing we can improve upon
                        !(maybe with an if statement?)
                        index1_1 = var1_index1 + j
                        index1_2 = var1_index1 - j
        
                        index2_1 = var2_index1 + i - j
                        index2_2 = var2_index1 - i + j
        
                        !This does the nitty-gritty of checking to see
                        !If any of the indexes are out of bounds and 
                        !obtaining their rmsds and coordinates
                        call getNeighbors(scaling1_0,scaling2_0,multiplier1_1,multiplier2_1,FMTorder1,&
                                          index1_1,index1_2,index2_1,index2_2,&
                                          var1_round0,var2_round0,&
                                          coords,min_rmsd,gradient,U,number_of_frames)
        
                        if (min_rmsd /= old_min_rmsd) then
				old_min_rmsd = min_rmsd
				old_gradient = gradient
				old_U = U
                        end if
                end do
        
                if (stop_flag) exit
        end do
end if

if (stop_flag) then
	write(filechannels(Ngrid),FMT=FMT6) old_min_rmsd
	cycle
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 ORDER 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Rinse and repeat
subcell = gridpath2//trim(adjustl(var1_filename0))//"_"//trim(adjustl(var2_filename0))

inquire(file=trim(subcell)//".dat",exist=subcell_existence)

order = 0

!For bug-testing
if (.false.) then
print *, ""
print *, ""
print *, "inquiring for: ", trim(subcell)//".dat" 
print *, "exists? ", subcell_existence
print *, ""
end if

if (subcell_existence) then
	call getRMSD_dp(subcell,coords,population,min_rmsd,gradient,U)

	number_of_frames = number_of_frames + population

	if (min_rmsd /= old_min_rmsd) then
		old_min_rmsd = min_rmsd
		old_gradient = gradient
		old_U = U
	end if
end if

if ((.not.(subcell_existence)).or.(force_Neighbors)) then
        stop_flag = .false.

        do i = 1, subcell0search_max
                neighbor_check = i

                index1_1 = var1_index0 + i
                index1_2 = var1_index0 - i

                index2_1 = var2_index0 + i
                index2_2 = var2_index0 - i
                
                if (index1_1 < bounds1) then
                
                        !Make the name of the subcell
                        var1_round = index1_1 * multiplier1_0
                        var2_round = var2_round0
                        write(var1_filename,FMT=FMTorder0) var1_round
                        write(var2_filename,FMT=FMTorder0) var2_round
                        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
                
                        !Check if it exists
                        inquire(file=trim(subcell)//".dat",exist=subcell_existence)
                        if (subcell_existence) then
                                call getRMSD_dp(subcell,coords,population,min_rmsd,gradient,U)
                        end if
                
                        number_of_frames = number_of_frames + population
                end if
                
                !Rinse and repeat
                if (index1_2 .ge. 0) then
                
                        var1_round = index1_2 * multiplier1_0
                        var2_round = var2_round0
                        write(var1_filename,FMT=FMTorder0) var1_round
                        write(var2_filename,FMT=FMTorder0) var2_round
                        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
                
                        inquire(file=trim(subcell)//".dat",exist=subcell_existence)
                        if (subcell_existence) then
                                call getRMSD_dp(subcell,coords,population,min_rmsd,gradient,U)
                        end if
                
                        number_of_frames = number_of_frames + population
                end if
                !Rinse and repeat
                if (index2_1 < bounds2) then
                
                        var1_round = var1_round0
                        var2_round = index2_1 * multiplier2_0
                        write(var1_filename,FMT=FMTorder0) var1_round
                        write(var2_filename,FMT=FMTorder0) var2_round
                        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
                
                        inquire(file=trim(subcell)//".dat",exist=subcell_existence)
                        if (subcell_existence) then
                                call getRMSD_dp(subcell,coords,population,min_rmsd,gradient,U)
                        end if
                
                        number_of_frames = number_of_frames + population
                end if
                
                !Rinse and repeat
                if (index2_2 .ge. 0) then
                
                        var1_round = var1_round0
                        var2_round = index2_2 * multiplier2_0
                        write(var1_filename,FMT=FMTorder0) var1_round
                        write(var2_filename,FMT=FMTorder0) var2_round
                        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))
                
                        inquire(file=trim(subcell)//".dat",exist=subcell_existence)
                        if (subcell_existence) then
                                call getRMSD_dp(subcell,coords,population,min_rmsd,gradient,U)
                        end if
                
                        number_of_frames = number_of_frames + population
                end if

                if (min_rmsd /= old_min_rmsd) then
                        stop_flag = .true.
			old_min_rmsd = min_rmsd
			old_gradient = gradient
			old_U = U
                end if

        do j = 1, i-1
        
                index1_1 = var1_index0 + j
                index1_2 = var1_index0 - j

                index2_1 = var2_index0 + i - j
                index2_2 = var2_index0 - i + j

                call getNeighbors(bounds1,bounds2,multiplier1_0,multiplier2_0,FMTorder1,&
                                  index1_1,index1_2,index2_1,index2_2,&
                                  0.0,0.0,&
                                  coords,min_rmsd,gradient,U,number_of_frames)

                if (min_rmsd /= old_min_rmsd) then
                        stop_flag = .true.
			old_min_rmsd = min_rmsd
			old_gradient = gradient
			old_U = U
                end if

        end do

        if (stop_flag) exit
        end do

end if

write(filechannels(Ngrid),FMT=FMT6) old_min_rmsd
end do

!This is to see whether old_U has its default value (all zeroes) or not
if (any(abs(old_U) > 1.0d-20)) gradient = matmul(old_U,old_gradient)
min_rmsd = old_min_rmsd

end subroutine checkState

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               getNeighbors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine checks the four subcells in the cartesian product
!               of relative indexes (index1_1,index1_2) x (index2_1, index2_2), with knowledge
!               of the scaling, the variables rounded down, and the counter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               scaling1                        INTEGER                         The ratio between the sizes of parent vs child
!                                                                               subcells with respect to variable1
!               scaling2                        INTEGER                         The ratio between the sizes of parent vs child
!                                                                               subcells with respect to variable2
!               multiplier1                     REAL                            The actual size of the child subcell
!                                                                               with respect to variable1
!               multiplier2                     REAL                            The actual size of the child subcell
!                                                                               with respect to variable2
!               FMTorder                        CHRACTER(6)                     The formatting for filenames associated with the
!                                                                               child subcell
!               var_round                       REAL                            The variables associated with the reference frame
!                                                                               after rounding to order N
!
!               index1_1                        INTEGER                         The first index to check for variable1
!               index1_2                        INTEGER                         The second index to check for variable1
!               index2_1                        INTEGER                         The first index to check for variable2
!               index2_2                        INTEGER                         The second index to check for variable2
!
!               coords                          REAL(DP),DIM(3,Natoms)          The coordinates defining the reference frame
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!               min_rmsd                        REAL(DP)                        The minimum rmsd found between the reference frame and
!                                                                               all frames in the subcell associated with it
!               gradient                        REAL(DP),DIM(3,Natoms)          The gradient associated with a frame
!               U                               REAL(DP),DIM(3,3)               The rotation matrix used to rotate a frame
!                                                                               and minimize its RMSD with respect to the reference frame
!
!               number_of_frames                INTEGER                         The total number of frames checked for the reference frame
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
!               var_round                       REAL                            The variable associated with the subcell in question,
!                                                                               rounded down to the order specified by multiplier
!
!               flag1                           LOGICAL                         Is relative index (index1_1) not out-of-bounds?
!               flag2                           LOGICAL                         Is relative index (index1_2) not out-of-bounds?
!               flag3                           LOGICAL                         Is relative index (index2_1) not out-of-bounds?
!               flag4                           LOGICAL                         Is relative index (index2_2) not out-of-bounds?
!
!               population                      INTEGER                         The number of frames in a subcell
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!               gridpath2//                     DAT                     The file houses all the information (var1,var2,coords,
!                 trim(subcell).dat                                     gradient) of each frame that is in the cell associated
!                                                                       with this file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getNeighbors(scaling1,scaling2,multiplier1,multiplier2,FMTorder,&
                        index1_1,index1_2,index2_1,index2_2,&
                        var1_round0,var2_round0,&
                        coords,min_rmsd,gradient,U,number_of_frames)
use ANALYSIS
use PARAMETERS
implicit none
logical :: flag1,flag2,flag3,flag4,subcell_existence
integer,intent(in) :: scaling1, scaling2
integer,intent(in) :: index1_1,index1_2,index2_1,index2_2
real,intent(in) :: multiplier1,multiplier2
real,intent(in) :: var1_round0,var2_round0
real :: var1_round,var2_round
!integer,intent(in) :: decimals1, decimals2
integer,intent(inout) :: number_of_frames
character(6),intent(in) :: FMTorder
character(FMTlength) :: var1_filename,var2_filename
character(100) :: subcell
real(dp),dimension(3,3),intent(inout) :: U
real(dp), intent(inout) :: min_rmsd
real(dp), dimension(3,Natoms), intent(inout) :: gradient
integer :: population
real(dp),intent(in), dimension(3,Natoms) :: coords

!We need to check if any of the indexes (neighbors!) are out of bounds
flag1 = index1_1 < scaling1
flag2 = index1_2 .ge. 0
flag3 = index2_1 < scaling2
flag4 = index2_2 .ge. 0

!We check indexes in pairs (var1, var2)
if ((flag1) .and. (flag3)) then

        !Make the name of the subcell
	var1_round = var1_round0 + index1_1 * multiplier1
	var2_round = var2_round0 + index2_1 * multiplier2
        write(var1_filename,FMT=FMTorder) var1_round
        write(var2_filename,FMT=FMTorder) var2_round
        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

	!Check if it exists
	inquire(file=trim(subcell)//".dat",exist=subcell_existence)

	!If it exists, get the rmsd of each frame in it
	if (subcell_existence) then
		call getRMSD_dp(subcell,coords,population,min_rmsd,gradient,U)
	end if

	!Increment however many frames we have
	number_of_frames = number_of_frames + population
end if

!Rinse and repeat
if ((flag1) .and. (flag4)) then

	var1_round = var1_round0 + index1_1 * multiplier1
	var2_round = var2_round0 + index2_2 * multiplier2
        write(var1_filename,FMT=FMTorder) var1_round
        write(var2_filename,FMT=FMTorder) var2_round
        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

	inquire(file=trim(subcell)//".dat",exist=subcell_existence)
	if (subcell_existence) then
		call getRMSD_dp(subcell,coords,population,min_rmsd,gradient,U)
	end if

	number_of_frames = number_of_frames + population
end if

!Rinse and repeat
if ((flag2) .and. (flag3)) then

	var1_round = var1_round0 + index1_2 * multiplier1
	var2_round = var2_round0 + index2_1 * multiplier2
        write(var1_filename,FMT=FMTorder) var1_round
        write(var2_filename,FMT=FMTorder) var2_round
        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

	inquire(file=trim(subcell)//".dat",exist=subcell_existence)
	if (subcell_existence) then
		call getRMSD_dp(subcell,coords,population,min_rmsd,gradient,U)
	end if

	number_of_frames = number_of_frames + population
end if

!Rinse and repeat
if ((flag2) .and. (flag4)) then

	var1_round = var1_round0 + index1_2 * multiplier1
	var2_round = var2_round0 + index2_2 * multiplier2
        write(var1_filename,FMT=FMTorder) var1_round
        write(var2_filename,FMT=FMTorder) var2_round
        subcell = gridpath2//trim(adjustl(var1_filename))//"_"//trim(adjustl(var2_filename))

	inquire(file=trim(subcell)//".dat",exist=subcell_existence)
	if (subcell_existence) then
		call getRMSD_dp(subcell,coords,population,min_rmsd,gradient,U)
	end if

	number_of_frames = number_of_frames + population
end if



end subroutine getNeighbors




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               FUNCTION GETRMSD_DP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  actually pretty self-explanatory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getRMSD_dp(subcell,coords,population,min_rmsd,gradient,U)

use ls_rmsd_original
use ANALYSIS
use PARAMETERS
implicit none
integer,intent(out) :: population
real(dp),intent(inout), dimension(3,Natoms) :: gradient
real(dp),dimension(3,Natoms) :: candidate_gradient
real(dp),intent(inout) :: min_rmsd
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
        open(filechannel1,file=trim(subcell)//".dat",form="unformatted")
else
        open(filechannel1,file=trim(subcell)//".dat")
end if

population = 0
do
        !Read the candidate frame
        if (unreadable_flag) then
                read(filechannel1,iostat=iostate)
	        if (iostate /= 0) exit
                read(filechannel1) ((coords2(i,j),i=1,3),j=1,Natoms)
        else
                read(filechannel1,FMT=FMT7,advance="no",iostat=iostate) ((coords2(i,j),i=1,3),j=1,Natoms)
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
                        read(filechannel1) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
                else
                        read(filechannel1,FMT=FMT3) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
                end if

        !Otherwise if it is low enough, accept the candidate frame (and save its rotation matrix!)
        else if (candidate_min_rmsd < threshold_rmsd) then
                min_rmsd = candidate_min_rmsd
                U = candidate_U

                if (unreadable_flag) then
                        read(filechannel1) ((gradient(j,k),j=1,3),k=1,Natoms)
                else
                        read(filechannel1,FMT=FMT3) ((gradient(j,k),j=1,3),k=1,Natoms)
                end if

                !In special cases, we exit immediately afterwards
                if (accept_first) exit

        !If it is too high, also reject the candidate frame
        else
                if (unreadable_flag) then
                        read(filechannel1) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
                else
                        read(filechannel1,FMT=FMT3) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
                end if
        end if

        !If in the "accept best" method:
        else

        !If it is low enough, accept the candidate frame (and save its rotation matrix!)
        if (candidate_min_rmsd < min_rmsd) then
                min_rmsd = candidate_min_rmsd
                U = candidate_U

                if (unreadable_flag) then
                        read(filechannel1) ((gradient(j,k),j=1,3),k=1,Natoms)
                else
                        read(filechannel1,FMT=FMT3) ((gradient(j,k),j=1,3),k=1,Natoms)
                end if

                !In special cases, we exit immediately afterwards
                if (accept_first) exit

        !Otherwise, just don't record the gradient
        else
                if (unreadable_flag) then
                        read(filechannel1) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
                else
                        read(filechannel1,FMT=FMT3) ((candidate_gradient(j,k),j=1,3),k=1,Natoms)
                end if
        end if

        end if

!	if (candidate_min_rmsd < min_rmsd) then
!		min_rmsd = candidate_min_rmsd
!		U = candidate_U
!
!                if (unreadable_flag) then
!                        read(filechannel1) ((gradient(i,j),i=1,3),j=1,Natoms)
!                else
!                        read(filechannel1,FMT=FMT3) ((gradient(i,j),i=1,3),j=1,Natoms)
!                end if
!
!                if (accept_first) exit
!	else
!                if (unreadable_flag) then
!                        read(filechannel1) ((candidate_gradient(i,j),i=1,3),j=1,Natoms)
!                else
!                        read(filechannel1,FMT=FMT3) ((candidate_gradient(i,j),i=1,3),j=1,Natoms)
!                end if
!	end if
end do
close(filechannel1)

end subroutine getRMSD_dp




end module interactMultipleGrids
