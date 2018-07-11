module FUNCTIONS
implicit none

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      PARTITION FUNCTION (QuickSort Alogrithm by Tony Hoare)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      INPUT:   matrix A, dim (rows,cols)       "values to be sorted"
!               matrix B, dim (rows,1)          "pointers to orignal indexes"
!               integer start_index             "in case of sorting only a part of the martix A, 
!                                                the starting index (before sorting) of the to-be-sorted part"
!               integer end_index               "in case of sorting only a part of the martix A, 
!                                                the final index (before sorting) of the to-be-sorted part"
!               integer var                     "the index of the column of the martix A that is used as sorting criteria"
!      OUTPUT:  integer pivot_index             "the index of the row that partitions matrix A according to var"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Column "var" of matrix A is partitioned
!       but only for the sub-matrix A (start_index:end_index)
!       Any position changes in A are mimicked in B.
!       The value of end_index is moved to pivot_index so that
!       all values of A(start:pivot) are <= A(pivot) and
!       all values of A(pivot:end) are >= A(pivot).
!       Pivot_index is then returned
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine partition(A,B,rows,cols,start_index,end_index,var,pivot_index)
implicit none
integer, intent(in) :: start_index,end_index,rows,cols,var
integer, intent(out) :: pivot_index
integer, dimension(rows,cols), intent(out) :: A
integer, dimension(rows,1), intent(out) :: B
integer :: i, j
integer :: test_value,pivot_value

pivot_value = A(end_index,var)
j = start_index-1
do i = start_index, end_index-1
        test_value = A(i,var)
        if (test_value.le.pivot_value) then
                j = j + 1
                call swapI(A,rows,cols,i,j)
                call swapI(B,rows,1,i,j)
        end if
end do
pivot_index = j+1
call swapI(A,rows,cols,pivot_index,end_index)
call swapI(B,rows,1,pivot_index,end_index)

end subroutine partition

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      SWAP FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IN/OUT:   matrix A, dim (rows,cols)       "values"
!      INPUT:    integer i1                      "spot1"
!                integer i2                      "spot2"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Spots 1 and 2 are swapped in matrix A
!       
!       swapR is for matrices of type real (vals)
!       swapI is for matrices of type integer (indexer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine swapR(A,rows,cols,i1, i2)
implicit none
integer, intent(in) :: rows, cols, i1, i2
real, dimension(rows,cols), intent(out) :: A
integer :: i
real, dimension(cols) :: val1, val2

do i = 1, cols
        val1(i) = A(i1,i)
        val2(i) = A(i2,i)
end do
do i = 1, cols
        A(i2,i) = val1(i)
        A(i1,i) = val2(i)
end do

end subroutine swapR

subroutine swapI(A,rows,cols,i1, i2)
implicit none
integer, intent(in) :: rows, cols, i1, i2
integer, dimension(rows,cols), intent(out) :: A
integer :: i
integer, dimension(cols) :: val1, val2

do i = 1, cols
        val1(i) = A(i1,i)
        val2(i) = A(i2,i)
end do
do i = 1, cols
        A(i2,i) = val1(i)
        A(i1,i) = val2(i)
end do

end subroutine swapI


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      QSORT2 FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IN/OUT:  matrix A, dim (rows,cols)       "values to be sorted"
!               matrix B, dim (rows,1)          "pointers to orignal indexes"
!  The following inputs control the sorting
!      INPUT:   integer start_index             "in case of sorting only a part of the martix A, 
!                                                the starting index (before sorting) of the to-be-sorted part"
!               integer end_index               "in case of sorting only a part of the martix A, 
!                                                the final index (before sorting) of the to-be-sorted part"
!               integer var                     "the index of the column of the martix A that is used as sorting criteria"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Column "var" of matrix A is sorted but only for the sub-matrix A (start_index:end_index).
!       Any position changes in A are mimicked in B.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ex.  A = [ 1.45  2.39 /       B = [ 1 /
!             0.98  7.45 /             2 /
!             2.11  3.00   ]           3   ]
!
!       QSORT2(A,B,3,2,   1,3,1) produces:
!       A = [ 0.98  7.45 /       B = [ 2 /
!             1.45  2.39 /             1 /
!             2.11  3.00   ]           3   ]
!
!       QSORT(A,B,3,2,    1,3,2) produces:
!       A = [ 1.45  2.39 /       B = [ 1 /
!             2.11  3.00 /             3 /
!             0.98  7.45   ]           2   ]
!
!       QSORT2(A,B,3,2,   1,2,1) produces:
!       A = [ 0.98  7.45 /       B = [ 2 /
!             1.45  2.39 /             1 /
!             2.11  3.00   ]           3   ]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine qsort2(A,B,rows,cols,start_index,end_index,var)
implicit none
integer, intent(in) :: rows, cols, start_index, end_index, var
integer :: pivot_index
integer, dimension(rows,cols), intent(out) :: A
integer, dimension(rows,1), intent(out) :: B

call partition(A,B,rows,cols,start_index,end_index,var,pivot_index)
if (start_index /= pivot_index) &
                call qsort2(A,B,rows,cols,start_index,pivot_index-1,var)
if (end_index /= pivot_index) &
                call qsort2(A,B,rows,cols,pivot_index+1,end_index,var)

end subroutine qsort2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      GRIDER FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      INPUT:   matrix A, dim (rows,cols)       "already-sorted values to be grided"
!               real gridline_spacing           "length of grid spacing"
!               real gridline_start             "value of first gridline"
!               integer max_gridlines           "the number of gridlines"
!               integer start_index             "start"
!               integer end_index               "end"
!               integer var                     "which column"
!      OUTPUT:  array grid, dim(Ngridmax)       "indexes of A for gridlines"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Column var of matrix A is assumed to be sorted
!       This functions grids submatrix A(start:end) into grid
!       There is a plain test description in addCells4 as well
!       
!       As a visual example, let's grid
!               A = [1.1, 2.7, 4.1, 4.2, 8.9, 9.6]
!
!       grider(grid1,A, 1.0, 0.0, 10, 6,1, 1, 6, 1)
!       grider(grid2,A, 2.0, 0.0, 10, 6,1, 1, 6, 1)
!       grider(grid3,A, 1.0, 2.0, 10, 6,1, 1, 6, 1)
!       grider(grid4,A, 1.0, 0.0,  5, 6,1, 1, 6, 1)
!       grider(grid5,A, 1.0, 0.0, 10, 6,1, 4, 6, 1)
!
!           A = | | 1.1 | 2.7 | | 4.1 4.2 | | | | 8.9 | 9.6 |
!       grid1 = 1,1,    2,    3,3,        5,5,5,5     6    
!
!           A = | 1.1 | 2.7 | 4.1 4.2 | | 8.9 9.6 | | | | | |
!       grid2 = 1,    2,    3,        5,5,        7,7,7,7,7
!
!           A = 1.1 | 2.7 | | 4.1 4.2 | | | | 8.9 | 9.6 | | |
!       grid3 =     2,    3,3,        5,5,5,5,    6,    7,7
!
!           A = | | 1.1 | 2.6 | | 4.1 4.2 |
!       grid4 = 1,1,    2,    3,3,        
!
!           A = | | 1.1 | 2.7 | | 4.1 4.2 | | | | 8.9 | 9.6 |
!       grid1 = 4,4,    4,    4,4,        5,5,5,5     6    
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine grider(grid,A,gridline_spacing,gridline_start,max_gridlines,&
                  rows,cols,start_index,end_index,var)
implicit none
integer, intent(in) :: rows, cols, var, max_gridlines
integer, intent(in) :: start_index, end_index
integer :: A_index, grid_index, i
real, intent(in) :: gridline_spacing, gridline_start
real, dimension(rows,cols), intent(in) :: A
real :: gridline
integer, dimension(max_gridlines), intent(out) :: grid

A_index = start_index
grid_index = 1
gridline = gridline_start + gridline_spacing
do
! RS: Question
! RS: gridline is the lower boundary of the grid 
! RS: b.c. grid_index-1 = 0 at the beginning
! RS: Is this correct? Shouldn't you use the upper boundary?
! RS: the first round until (A_index == end_index) is pointless
! RS: Since A is already sorted, I think we can use a much algorithm to grid
! RS: Like the "finding the root" homework
! RS: Let's talk about this tomorrow

! RS: Do you plan to change this later?
!                       KF: yes, it will be the next push, most likely
!       gridline = gridline_start + (grid_index-1)*gridline_spacing
        if (gridline < A(A_index,var)) then
                grid(grid_index) = A_index
                grid_index = grid_index + 1
                gridline = gridline + gridline_spacing
        else if (A_index == end_index) then
!               do
!                       if (grid_index > max_gridlines) exit
!                       grid(grid_index) = A_index + 1
!                       grid_index = grid_index + 1
!               end do
!               exit
                A_index = A_index + 1
                do i = grid_index, max_gridlines
                        grid(i) = A_index
                end do
                exit
        else
                A_index = A_index + 1
        end if
end do

end subroutine grider


end module FUNCTIONS

