module FUNCTIONS
use DOUBLE
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
real(dp), dimension(rows,cols), intent(out) :: A
integer, dimension(rows,1), intent(out) :: B
integer :: i, j
integer :: test_value,pivot_value

pivot_value = A(end_index,var)
j = start_index-1
do i = start_index, end_index-1
        test_value = A(i,var)
        if (test_value.le.pivot_value) then
                j = j + 1
                call swapD(A,rows,cols,i,j)
                call swapI(B,rows,1,i,j)
        end if
end do
pivot_index = j+1
call swapD(A,rows,cols,pivot_index,end_index)
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
!       swapD is for matrices of type real(dp) (vals)
!       swapR is for matrices of type real (vals)
!       swapI is for matrices of type integer (indexer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine swapD(A,rows,cols,i1, i2)
implicit none
integer, intent(in) :: rows, cols, i1, i2
real(dp), dimension(rows,cols), intent(out) :: A
integer :: i
real(dp), dimension(cols) :: val1, val2

do i = 1, cols
        val1(i) = A(i1,i)
        val2(i) = A(i2,i)
end do
do i = 1, cols
        A(i2,i) = val1(i)
        A(i1,i) = val2(i)
end do

end subroutine swapD

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
real(dp), dimension(rows,cols), intent(out) :: A
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



subroutine chooseINT(Nintegers,lowerBound,upperBound,randomIntegers)
implicit none
integer,intent(in) :: Nintegers, lowerBound, upperBound
integer,dimension(Nintegers),intent(out) :: randomIntegers
real :: rand_num
integer :: rand_int
integer :: Nrange
logical :: stop_flag
integer :: i,j

Nrange = upperBound - lowerBound + 1

do i = 1, Nintegers
        do
                stop_flag = .false.
                rand_num = rand()
                rand_int = floor(rand_num*Nrange)
                if ((rand_int < 0).or.(rand_int == Nrange)) cycle
                do j = 1, i-1
                        if (rand_int == randomIntegers(j)) stop_flag = .true.
                end do
                if (stop_flag) cycle
                randomIntegers(i) = rand_int
                exit
        end do
end do

randomIntegers = randomIntegers + lowerBound

end subroutine chooseINT

subroutine QRDECOMP(A,m,n,Q,R)
    implicit none
    integer,intent(in) :: m,n
    real(dp), dimension(m,n), intent(in) :: A
    real(dp), dimension(m,m), intent(out) :: Q
    real(dp), dimension(m,n), intent(out) :: R
    real(dp), dimension(m,m) :: Qn
    real(dp), dimension(m,n) :: Rn

    real(dp), dimension(m) :: x
    real(dp) :: u, xnorm, tau
    integer :: i, j, t

    do i = 1, m
        do j = 1, m
            if (i == j) then
                Q(i,j) = 1
            else
                Q(i,j) = 0
            end if
        end do
    end do

    Rn = A

    do t = 1, min(m-1,n)

        do i = 1, m
            do j = 1, m
                if (i == j) then
                    Qn(i,j) = 1
                else
                    Qn(i,j) = 0
                end if
            end do
        end do

        x = Rn(:,t)
        xnorm = sqrt(sum(x(t:m)**2,dim=1))
        x(t) = x(t) - sign(xnorm,Rn(t,t))
        xnorm = sqrt(sum(x(t:m)**2,dim=1))
        x = x / xnorm

        Qn(t:m,t:m) = Qn(t:m,t:m) - 2*matmul(&
                reshape(x(t:m),(/m-t+1,1/)),&
                reshape(x(t:m),(/1,m-t+1/)))

        Rn = matmul(transpose(Qn),Rn)

        Q = matmul(Q,transpose(Qn))

    end do

    R = matmul(transpose(Q),A)

end subroutine QRDECOMP

subroutine testQRDECOMP()
    implicit none
    integer,parameter :: m = 3
    integer,parameter :: n = 3
    character(1) :: ntext
    character(1) :: mtext
    real(dp), dimension(m,n) :: A, B, R
    real(dp), dimension(m,m) :: Qtotal, Q

    integer i,j,k

    write(mtext,FMT="(I1)") m
    write(ntext,FMT="(I1)") n

    do i = 1, m
        do j = 1, n
            A(i,j) = rand()
            R(i,j) = A(i,j)

            if (i == j) then
                Qtotal(i,j) = 1
            else
                Qtotal(i,j) = 0
            end if
        end do
    end do

    A(1,1) = 1.0d0
    A(1,2) = 2.0d0
    A(1,3) = 3.0d0

    A(2,1) = 2.0d0
    A(2,2) = 4.0d0
    A(2,3) = 6.0d0

    print *, ""
    print *, "A"
    print *, ""

    do i = 1, n
        write(6,FMT="("//mtext//"(F14.8))") A(:,i)
    end do

    call QRDECOMP(A,m,n,Q,R)

    A = matmul(transpose(Q),A)
    B = matmul(Q,R)

    print *, ""
    print *, "Q"
    print *, ""

    do j = 1, m
        write(6,FMT="("//mtext//"(F14.8))") Q(:,j)
    end do

    print *, ""
    print *, "R"
    print *, ""

    do j = 1, n
        write(6,FMT="("//mtext//"(F14.8))") R(:,j)
    end do

    print *, ""
    print *, "QR"
    print *, ""

    do j = 1, n
        write(6,FMT="("//mtext//"(F14.8))") B(:,j)
    end do

    print *, ""
    print *, "QtA"
    print *, ""

    do j = 1, n
        write(6,FMT="("//mtext//"(F14.8))") A(:,j)
    end do

end subroutine testQRDECOMP

subroutine CLS(A,m,n,C,p,d,b,x)
    implicit none
    integer,intent(in) :: m, n, p
    character(1) :: ntext
    character(1) :: ptext
    character(1) :: mptext
    real(dp), dimension(m,n), intent(in) :: A
    real(dp), dimension(p,n), intent(in) :: C
    real(dp), dimension(m), intent(in) :: b
    real(dp), dimension(p), intent(in) :: d
    real(dp), dimension(n), intent(out) :: x

    real(dp), dimension(m+p,n) :: S,R
    real(dp), dimension(m,m+p) :: Q1,Q1R
    real(dp), dimension(p,m+p) :: Q2,Q2R
    real(dp), dimension(m+p,m+p) :: Q2Tq
    real(dp), dimension(m+p,p) :: Q2T,Q2Tr
    real(dp), dimension(m+p,m+p) :: Q,Qt,QQ
    real(dp), dimension(m+p,n) :: QR
    real(dp), dimension(m+p) :: u,ce
    real(dp), dimension(p) :: w
    real(dp), dimension(m+p) :: y
    real(dp), dimension(m) :: bout
    real(dp), dimension(p,1) :: dout
    real(dp), dimension(m+p,1) :: ceout

    integer :: i, j

    write(ntext,FMT="(I1)") n
    write(ptext,FMT="(I1)") p
    write(mptext,FMT="(I1)") m+p

    do i = 1, m
        do j = 1, n
            S(i,j) = A(i,j)
        end do
        y(i) = b(i)
    end do

    do i = 1, p
        do j = 1, n
            S(m+i,j) = C(i,j)
        end do
        y(m+i) = d(i)
    end do

    call QRDECOMP(S,m+p,n,Q,R)

    print *, ""
    print *, "S ( 'stacked' [ A^T C^T ]^T )"
    print *, ""

    do i = 1, m+p
        write(6,FMT="("//ntext//"(F10.6))") S(i,:)
    end do

    do i = 1, n
        if (all(abs(R(:,n)) < 1.0d-9)) then
            print *, ""
            print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print *, ""
            print *, "ERROR: S is not left-invertible"
            print *, " (linearly independent columns)"
            print *, ""
            print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print *, ""
        end if
    end do

    Qt = transpose(Q)

    print *, ""
    print *, "Q"
    print *, ""

    do i = 1, m+p
        write(6,FMT="("//mptext//"(F10.6))") Q(i,:)
    end do

    print *, ""
    print *, "R"
    print *, ""

    do i = 1, m+p
        write(6,FMT="("//ntext//"(F10.6))") R(i,:)
    end do

    QQ = matmul(Qt,Q)

    print *, ""
    print *, "Q^T Q"
    print *, ""

    do i = 1, m+p
        write(6,FMT="("//mptext//"(F10.6))") QQ(i,:)
    end do

    QR = matmul(Q,R)

    print *, ""
    print *, "QR"
    print *, ""

    do i = 1, m+p
        write(6,FMT="("//ntext//"(F10.6))") QR(i,:)
    end do

    Q1 = Q(1:m,:)

    print *, ""
    print *, "Q1"
    print *, ""

    do i = 1, m
        write(6,FMT="("//mptext//"(F10.6))") Q1(i,:)
    end do

    Q1R = matmul(Q1,R)

    print *, ""
    print *, "Q1R"
    print *, ""

    do i = 1, m
        write(6,FMT="("//ntext//"(F10.6))") Q1R(i,:)
    end do

    Q2 = Q(m+1:m+p,:)

    print *, ""
    print *, "Q2"
    print *, ""

    do i = 1, p
        write(6,FMT="("//mptext//"(F10.6))") Q2(i,:)
    end do

    Q2R = matmul(Q2,R)

    print *, ""
    print *, "Q2R"
    print *, ""

    do i = 1, p
        write(6,FMT="("//ntext//"(F10.6))") Q2R(i,:)
    end do

    call QRDECOMP(transpose(Q2),m+p,p,Q2Tq,Q2Tr)

    print *, ""
    print *, "Q2Tq"
    print *, ""

    do i = 1, m+p
        write(6,FMT="("//mptext//"(F10.6))") Q2Tq(i,:)
    end do

    print *, ""
    print *, "Q2Tr"
    print *, ""

    do i = 1, m+p
        write(6,FMT="("//ptext//"(F10.6))") Q2Tr(i,:)
    end do

    Q2T = matmul(Q2Tq,Q2Tr)

    print *, ""
    print *, "Q2Tq Q2Tr"
    print *, ""

    do i = 1, m+p
        write(6,FMT="("//ptext//"(F10.6))") Q2T(i,:)
    end do

    if (abs(R(m+p,m+p-1)) < 1.0d-10) then
        print *, ""
        print *, "!!!!!!!!!!!!!!!!!!!!!!!"
        print *, ""
        print *, "Exact solution found!"
        print *, ""
        print *, "!!!!!!!!!!!!!!!!!!!!!!!"
        print *, ""

        print *, ""
        print *, "d"
        print *, ""

        do i = 1, p
            write(6,FMT="(F10.6)") d(i)
        end do

        call BackSubstitute(&
            transpose(Q2Tr),p,m+p,u,d)

        dout = matmul(transpose(Q2Tr),&
                  reshape(u,(/m+p,1/)))

        print *, ""
        print *, "dout"
        print *, ""

        do i = 1, p
            write(6,FMT="(F10.6)") dout(i,1)
        end do

!       ce = reshape(&
!            matmul(transpose(Q2Tq),&
!            matmul(transpose(Q1),&
!            reshape(b,(/m,1/)))),&
!            (/m+p/)) - u

        ce  = -u

        print *, ""
        print *, "ce"
        print *, ""

        do i = 1, m+p
            write(6,FMT="(F10.6)") ce(i)
        end do

        call ForwardSubstitute(Q2Tr,m+p,p,&
        w,reshape(ce,(/m+p/)))

!       w = -1.0d0

        ceout = matmul(Q2Tr,&
                   reshape(w,(/p,1/)))

        print *, ""
        print *, "ceout"
        print *, ""

        do i = 1, m+p
            write(6,FMT="(F10.6)") ceout(i,1)
        end do

        y = reshape(&
            matmul(transpose(Q1),&
                reshape(b,(/m,1/))) - &
            matmul(transpose(Q2),&
                reshape(w,(/p,1/))),&
            (/m+p/))

        print *, ""
        print *, "y"
        print *, ""

        do i = 1, m+p
            write(6,FMT="(F10.6)") y(i)
        end do

        call ForwardSubstitute(R,m+p,n,x,y)

        print *, ""
        print *, "x"
        print *, ""

        do i = 1, n
            write(6,FMT="(F10.6)") x(i)
        end do

        y = reshape(matmul(R,&
                reshape(x,(/n,1/))),&
                (/m+p/))

        print *, ""
        print *, "R x"
        print *, ""

        do i = 1, m+p
            write(6,FMT="(F10.6)") y(i)
        end do



        bout = matmul(A,x) - b

        print *, ""
        print *, "A x - b"
        print *, ""

        do i = 1, m
            write(6,FMT="(F10.6)") bout(i)
        end do

        x = 1.0d0 / n

        bout = matmul(A,x) - b

        print *, ""
        print *, "suppose x = 1/n"
        print *, ""
        print *, "A x - b"
        print *, ""

        do i = 1, m
            write(6,FMT="(F10.6)") bout(i)
        end do

    end if

end subroutine CLS

subroutine CLS2(A,m,n,C,p,d,b,x)
    implicit none
    integer,intent(in) :: m, n, p
    character(2) :: ntext
    character(2) :: ptext
    character(2) :: mptext
    character(2) :: nptext
    real(dp), dimension(m,n), intent(in) :: A
    real(dp), dimension(p,n), intent(in) :: C
    real(dp), dimension(m), intent(in) :: b
    real(dp), dimension(p), intent(in) :: d
    real(dp), dimension(n), intent(out) :: x

    real(dp), dimension(n+p,n+p) :: S,R
    real(dp), dimension(n+p,n+p) :: Q,Qt,QQ
    real(dp), dimension(n+p,n+p) :: QR
    real(dp), dimension(n,n) :: H
    real(dp), dimension(n+p,n+p) :: E
    real(dp), dimension(n,1) :: z
    real(dp), dimension(n+p,1) :: w
    real(dp), dimension(n+p) :: y, u
    real(dp), dimension(m,1) :: bout

    integer :: i, j

!   write(ntext,FMT="(I2)") n
!   write(ptext,FMT="(I2)") p
!   write(mptext,FMT="(I2)") m+p
!   write(nptext,FMT="(I2)") n+p

!   print *, ""
!   print *, "A"
!   print *, ""

!   do i = 1, m
!       write(6,FMT="("//ntext//"(F14.10))") A(i,:)
!   end do


!   print *, ""
!   print *, "b"
!   print *, ""

!   do i = 1, m
!       write(6,FMT="(F14.10)") b(i)
!   end do

    H = matmul(transpose(A),A)

    E(1:n,1:n) = H
    E(1:n,n+1:n+p) = transpose(C)
    E(n+1:n+p,1:n) = C
    E(n+1:n+p,n+1:n+p) = 0.0d0

    call QRDECOMP(E,n+p,n+p,Q,R)

!   print *, ""
!   print *, "E"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="("//nptext//"(F10.6))") E(i,:)
!   end do

!    do i = 1, n
!        if (all(abs(R(:,n)) < 1.0d-9)) then
!            print *, ""
!            print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!            print *, ""
!            print *, "ERROR: S is not left-invertible"
!            print *, " (linearly independent columns)"
!            print *, ""
!            print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!            print *, ""
!        end if
!    end do

    Qt = transpose(Q)

!   print *, ""
!   print *, "Q"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="("//nptext//"(F10.6))") Q(i,:)
!   end do

!   print *, ""
!   print *, "R"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="("//nptext//"(F10.6))") R(i,:)
!   end do

!   QQ = matmul(Qt,Q)

!   print *, ""
!   print *, "Q^T Q"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="("//nptext//"(F10.6))") QQ(i,:)
!   end do

!   QR = matmul(Q,R)

!   print *, ""
!   print *, "QR"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="("//nptext//"(F10.6))") QR(i,:)
!   end do

    z = matmul(transpose(A),&
           reshape(b,(/m,1/)))

    w(1:n,1) = z(1:n,1)
    w(n+1:p,1) = d(:)

!   print *, ""
!   print *, "w"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="(F10.6)") w(i,1)
!   end do

    w = matmul(Qt,w)
    y = reshape(w,(/n+p/))

!   print *, ""
!   print *, "y = Qt w"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="(F10.6)") y(i)
!   end do

    call ForwardSubstitute(R,n+p,n+p,u,y)

!   print *, ""
!   print *, "u"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="(F10.6)") u(i)
!   end do

!   w = matmul(R,reshape(u,(/n+p,1/)))
!   y = reshape(w,(/n+p/))

!   print *, ""
!   print *, "y = R u"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="(F10.6)") y(i)
!   end do

    x = u(1:n)

!   bout = matmul(A,reshape(x,(/n,1/)))

!   print *, ""
!   print *, "bout = A x"
!   print *, ""

!   do i = 1, m
!       write(6,FMT="(F14.10)") bout(i,1)
!   end do

!   print *, ""
!   print *, "||bout-b|| = ", sqrt(sum((bout(:,1)-b)**2))
!   print *, "||bout-b|| = ", sqrt(sum((bout(1:12,1)-b(1:12))**2)/4)
!   print *, ""

!   x = 1.0d0 / n

!   bout = matmul(A,reshape(x,(/n,1/)))

!   print *, ""
!   print *, "suppose x = 1/n"
!   print *, ""
!   print *, "bout = A x"
!   print *, ""

!   do i = 1, m
!       write(6,FMT="(F10.6)") bout(i,1)
!   end do

!   print *, ""
!   print *, "||bout-b|| = ", sqrt(sum((bout(:,1)-b)**2))
!   print *, "||bout-b|| = ", sqrt(sum((bout(1:12,1)-b(1:12))**2)/4)
!   print *, ""

!   x = u(1:n)
!   x(1) = x(1) + 0.01d0 * u(1)
!   x = x / (sum(x))

!   bout = matmul(A,reshape(x,(/n,1/)))

!   print *, ""
!   print *, "suppose x = normalized(x + 0.01 * x1 * e1)"
!   print *, ""
!   print *, "bout = A x"
!   print *, ""

!   do i = 1, m
!       write(6,FMT="(F10.6)") bout(i,1)
!   end do

!   print *, ""
!   print *, "||bout-b|| = ", sqrt(sum((bout(:,1)-b)**2))
!   print *, "||bout-b|| = ", sqrt(sum((bout(1:12,1)-b(1:12))**2)/4)
!   print *, ""

!   x = u(1:n)

end subroutine CLS2

subroutine LS(A,m,n,b,x)
    implicit none
    integer,intent(in) :: m, n
    character(2) :: ntext
    real(dp), dimension(m,n), intent(in) :: A
    real(dp), dimension(m), intent(in) :: b
    real(dp), dimension(n), intent(out) :: x

    real(dp), dimension(n,n) :: H, Q, R, Qt
    real(dp), dimension(n,1) :: z
    real(dp), dimension(n,1) :: w
    real(dp), dimension(n) :: y, u

    integer :: i, j

    H = matmul(transpose(A),A)

    call QRDECOMP(H,n,n,Q,R)

    Qt = transpose(Q)

    z = matmul(transpose(A),&
           reshape(b,(/m,1/)))

    w(1:n,1) = z(1:n,1)

    w = matmul(Qt,w)
    y = reshape(w,(/n/))

    call ForwardSubstitute(R,n,n,u,y)

    x = u(1:n)

end subroutine LS

subroutine testCLS()
    implicit none
    integer, parameter :: m = 4
    integer, parameter :: p = 1
    integer, parameter :: n = 4
    character(1) :: mtext
    character(1) :: ptext
    character(1) :: ntext
    real(dp), dimension(m,n) :: A
    real(dp), dimension(p,n) :: C
    real(dp), dimension(m) :: b
    real(dp), dimension(p) :: d
    real(dp), dimension(n) :: x

    real(dp), dimension(m+p,n) :: S,R
    real(dp), dimension(m+p,m+p) :: Q

    integer :: i, j

    write(mtext,FMT="(I1)") m
    write(ptext,FMT="(I1)") p
    write(ntext,FMT="(I1)") n

    do i = 1, m
        do j = 1, n
            A(i,j) = rand()
        end do
    end do

    do i = 1, p
        do j = 1, n
            C(i,j) = rand()
        end do
    end do

    C = 1.0d0

    b = 0.0d0
    d = 1.0d0

    print *, ""
    print *, "A"
    print *, ""

    do i = 1, m
        write(6,FMT="("//ntext//"(F10.6))") A(i,:)
    end do

    print *, ""
    print *, "C"
    print *, ""

    do i = 1, p
        write(6,FMT="("//ntext//"(F10.6))") C(i,:)
    end do

    call CLS2(A,m,n,C,p,d,b,x)

    print *, ""
    print *, "x"
    print *, ""

    do i = 1, n
        write(6,FMT="(F10.6)") x(i)
    end do

    b = reshape(matmul(A,reshape(x,(/n,1/))),(/m/))

    print *, ""
    print *, "b = A x"
    print *, ""

    do i = 1, m
        write(6,FMT="(F10.6)") b(i)
    end do

    print *, ""
    print *, "||b|| = ", sqrt(sum(b**2))
    print *, ""

end subroutine testCLS

!This is for LOWER (LEFT) TRIANGULAR matrix A
!with equation of form Ax = b
subroutine BackSubstitute(A,m,n,x,b)
    implicit none
    integer,intent(in) :: m, n
    real(dp), dimension(m,n), intent(in) :: A
    real(dp), dimension(n), intent(out) :: x
    real(dp), dimension(m), intent(in) :: b
    
    real(dp) :: other_terms

    integer :: i, j

    x = 0.0d0

    do i = 1, m
        if (abs(A(i,i)) < 1.0d-15) then
            print *, ""
            print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print *, ""
            print *, "ERROR in back substitution... EXITING"
            print *, ""
            print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print *, ""

            exit
        end if

        other_terms = 0.0d0
        do j = 1, i-1
            other_terms = other_terms + x(j)*A(i,j)
        end do

        x(i) = (b(i) - other_terms) / A(i,i)
    end do

end subroutine BackSubstitute

!This is for UPPER (RIGHT) TRIANGULAR matrix A
!with equation of form Ax = b
subroutine ForwardSubstitute(A,m,n,x,b)
    implicit none
    integer,intent(in) :: m, n
    real(dp), dimension(m,n), intent(in) :: A
    real(dp), dimension(n), intent(out) :: x
    real(dp), dimension(m), intent(in) :: b
    
    real(dp) :: other_terms, dividor
    integer :: extra_pivots

    integer :: i, j

    extra_pivots = m+1
    x = 0.0d0

    do i = 1, m
        dividor = A(m-i+1,m-i+1)

        if (abs(dividor) < 1.0d-15) then
            if (extra_pivots <= n) then
                !need to add this
                extra_pivots = extra_pivots + 1
            end if

            other_terms = 0.0d0
            do j = (m-i+1), m-1
                other_terms = other_terms + x(j)*A(m-i+1,j)
            end do

            if (abs(other_terms - b(m-i+1)) < 1.d-10) then
                dividor = 1.0d0
            else
                print *, ""
                print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                print *, ""
                print *, "ERROR in forward substitution... EXITING"
                print *, ""
                print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                print *, ""

                x = 0.0d0
                exit
            end if
        end if

        other_terms = 0.0d0
        do j = (m-i+1), m
            other_terms = other_terms + x(j)*A(m-i+1,j)
        end do

        x(m-i+1) = (b(m-i+1) - other_terms) / dividor
    end do

end subroutine ForwardSubstitute

end module FUNCTIONS

