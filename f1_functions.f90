


module f1_functions
implicit none

contains


!This partitions A (rows by cols) and B (rows) by the first column of A
!Used in quicksort
subroutine partition(A,B,rows,cols,p,r,m,var)
!The arguments are (vals,indexer,overcrowdN,Nvar,1,overcrowdN,???,1) from divyup

implicit none
integer, intent(in) :: p, r, rows, cols, var
integer, intent(out) :: m
real, dimension(rows,cols), intent(out) :: A
integer, dimension(rows,1), intent(out) :: B
! RS: The name of the variables should be systematic
! RS: for example, val and var has been used in the code throughoutly for something else
! RS: use "_" to connect terms to make the name of the variable more intuitive
! RS: Even thought this is not a priority, please be mindful
integer :: i, j, i1, i2
! RS: Did you use i1, i2, and val1 in this subroutine at all?
real :: val1,val2,val3

! RS: I see you are doing a bubble sort but don't you need to loop through val3 as well?
! RS: https://www.youtube.com/watch?v=nmhjrI-aW5o
val3 = A(r,var)
j = p-1
do i = p, r-1
        val2 = A(i,var)
        if (val2.le.val3) then
                j = j + 1
                call swapR(A,rows,cols,i,j)
                call swapI(B,rows,1,i,j)
        end if
end do
m = j+1
call swapR(A,rows,cols,m,r)
call swapI(B,rows,1,m,r)

end subroutine partition

!This swaps two rows (i1 and i2) of a matrix A (rows by cols)
!Used in partition
subroutine swapR(A,rows,cols,i1, i2)
!The arguments are (vals,overcrowdN,Nvar,i,j) 
!first three args are from from divyup
!last two args are from partition
implicit none
integer, intent(in) :: rows, cols, i1, i2
real, dimension(rows,cols), intent(out) :: A
integer :: i, j
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

!This is the same as above but for integers
! RS: for indexer??
!Used in partition
subroutine swapI(A,rows,cols,i1, i2)
!The arguments are (indexer,overcrowdN,1,i,j) 
!first three args are from from divyup
!last two args are from partition
implicit none
integer, intent(in) :: rows, cols, i1, i2
integer, dimension(rows,cols), intent(out) :: A
integer :: i, j
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

!rename please!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This is the quicksort algorithm
!It's pretty interesting
!It sorts matrices A and B by the column of A indicated by var
recursive subroutine qsort(A,B,rows,cols,p,r,var)
!The arguments are (vals,indexer,overcrowdN,Nvar,1,overcrowdN,1) from divyup
implicit none
integer, intent(in) :: rows, cols, p, r, var
integer :: m
real, dimension(rows,cols), intent(out) :: A
integer, dimension(rows,1), intent(out) :: B

call partition(A,B,rows,cols,p,r,m,var)
if (p /= m) call qsort(A,B,rows,cols,p,m-1,var)
if (r /= m) call qsort(A,B,rows,cols,m+1,r,var)

end subroutine qsort








!This goes through the values of A in column var, rows [p1, p2)
!(A NEEDS TO BE PRE-ORDERED) and marks in grid rows [r1, r2]
!It indicates which state p is directly following gridline r
subroutine grider(grid,A,spacings,r0,Ngridmax,rows,cols,p1,p2,r1,r2,var)
implicit none
integer, intent(in) :: rows, cols, var, Ngridmax
integer, intent(in) :: p1, p2, r1, r2
integer :: p, r
real, intent(in) :: spacings, r0
real, dimension(rows,cols), intent(in) :: A
real :: grid_line
integer, dimension(Ngridmax), intent(out) :: grid

p = p1
r = 0
do
        grid_line = r0 + r*spacings
        if (grid_line < A(p,var)) then
                grid(r1+r) = p
                r = r + 1
        else if (p == p2) then
                do
                        if (r1+r > r2) exit
                        grid(r1+r) = p + 1
                        r = r + 1
                end do
                exit
        else
                p = p + 1
        end if
end do

end subroutine grider








!Used to increase the array length of a too-full array
subroutine growI(A,B,rows1,cols1,rows2,cols2)
implicit none
integer, intent(in) :: rows1,cols1,rows2,cols2
integer :: i,j
integer, dimension(rows1,cols1), intent(in) :: A
integer, dimension(rows2,cols2), intent(out) :: B

do i = 1, rows1
        do j = 1, cols1
                B(i,j) = A(i,j)
        end do
end do

end subroutine growI

!Same as above but for real numbers
subroutine growR(A,B,rows1,cols1,rows2,cols2)
implicit none
integer, intent(in) :: rows1,cols1,rows2,cols2
integer :: i,j
real, dimension(rows1,cols1),intent(in) :: A
real, dimension(rows2,cols2), intent(out) :: B

do i = 1, rows1
        do j = 1, cols1
                B(i,j) = A(i,j)
        end do
end do

end subroutine growR



end module f1_functions




