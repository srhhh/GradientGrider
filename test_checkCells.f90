program test_checkCells
use checkCells4
use f1_variables
use f1_parameters

implicit none
real, dimension(Nvar) :: vals
real, dimension(3*Natoms) :: coords, gradient
real, dimension(6*Natoms) :: coords_Cells
double precision :: min_rmsd
integer :: i, j, k, rand_subcell_position, rand_file_position, state1
integer,dimension(counter0_max) :: counter0
integer,dimension(counter1_max) :: counter1
integer,dimension(counter2_max) :: counter2
integer,dimension(counter3_max) :: counter3
logical :: flag1
character(50) :: subcell

!Instead of printing to the terminal, we print to the progress files
!First we check if it exists and just wipe
!We want to start from a clean slate
inquire(file=trim(path4)//trim(progressfile),exist=flag1)
if (flag1) call system("rm "//trim(path4)//trim(progressfile))
open(70,file=trim(path4)//trim(progressfile),status="new")
write(70,FMT="(A50)") "Let's start"
write(70,FMT="(A50)") ""
close(70)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          COUNTER ARRAYS RETRIEVAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!remark: to TEST the grid, the grid has to actually have been MADE beforehand
!        if you want to quickly make the grid, do into getCells.f90
!        and make the program quit after maybe twenty trajectories.
!        This program will work regardless

open(70,file=trim(path4)//"counter0.txt",status="old")
do i = 1, counter0_max
        read(70,FMT="(I8)") counter0(i)
end do
close(70)

open(70,file=trim(path4)//"counter1.txt",status="old")
do i = 1, counter1_max
        read(70,FMT="(I8)") counter1(i)
end do
close(70)

open(70,file=trim(path4)//"counter2.txt",status="old")
do i = 1, counter2_max
        read(70,FMT="(I8)") counter2(i)
end do
close(70)

open(70,file=trim(path4)//"counter3.txt",status="old")
do i = 1, counter3_max
        read(70,FMT="(I8)") counter3(i)
end do
close(70)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          RANDOM FRAME RETRIEVAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rand_subcell_position = 29
rand_file_position = 12

!Get a list of all subcells; pick a random one
call system("ls -p "//trim(path3)//" | grep '.dat' > "//trim(path4)//trim(temporaryfile1))
open(70,file=trim(path4)//trim(temporaryfile1))
do i = 1, rand_subcell_position
        read(70,FMT="(A50)") subcell
end do
close(70)

!Open up that subcell; pick a random frame
open(70,file=trim(path4)//trim(progressfile),position="append")
write(70,FMT="(A50)") "For the random frame, looking in ", trim(subcell)
write(70,FMT="(A50)") ""
close(70)

open(70,file=trim(path3)//trim(subcell))
i = 0
do
        i = i + 1
        if (i < rand_file_position) then
                read(70,FMT="(A50)",iostat=state1) subcell
                if (state1 /= 0) then
                        i = i - 1
                        rewind(70)
                end if
        else
                read(70,FMT=FMT1,advance="no",iostat=state1) (vals(j),j=1,Nvar)
                if (state1 /= 0) then
                        i = i - 1
                        close(70)
                        open(70,file=trim(path3)//trim(subcell))
                        cycle
                end if
                read(70,FMT=FMT2) (coords(j),j=1,3*Natoms)
                exit 
        end if
end do
close(70)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               FRAME MODIFICATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(70,file=trim(path4)//trim(progressfile),position="append")
write(70,*) "The original coordinates are: "
write(70,*) coords
write(70,*) ""
write(70,*) "With variables: ", vals(1), " A and ", vals(2), " A"
write(70,*) ""
close(70)

!This is how I choose to change the coordinates
coords(3) = coords(3) + 0.1
coords(4) = coords(4) - 0.000
coords(18) = coords(18) + .000

open(70,file=trim(path4)//trim(progressfile),position="append")
write(70,FMT="(A50)") "Making these changes to it:"
write(70,FMT="(A50)") ""
write(70,*) "coords(3) = coords(3) + 0.1"
write(70,*) "coords(4) = coords(4) - 0.000"
write(70,*) "coords(18) = coords(18) + .000"
write(70,FMT="(A50)") ""
close(70)

call getVar1(coords,Natoms,vals(1))
call getVar2(coords,Natoms,vals(2))

open(70,file=trim(path4)//trim(progressfile),position="append")
write(70,*) "The modified coordinates are: "
write(70,*) coords
write(70,*) ""
write(70,*) "With variables: ", vals(1), " A and ", vals(2), " A"
write(70,*) ""
close(70)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               GRID CHECKING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call checkState(coords,coords_Cells,min_rmsd,&
                counter0,counter1,counter2,counter3)

call getVar1(coords_Cells,Natoms,vals(1))
call getVar2(coords_Cells,Natoms,vals(2))

open(70,file=trim(path4)//trim(progressfile),position="append")
write(70,*) "The closests coordinates are: "
write(70,*) coords_Cells(1:3*Natoms)
write(70,*) ""
write(70,*) "With variables: ", vals(1), " A and ", vals(2), " A"
write(70,*) ""
write(70,*) ""
write(70,*) "This has an rmsd of: ", min_rmsd
write(70,*) ""
write(70,*) ""
write(70,*) "The closest frame has a gradient as follows: "
write(70,*) coords_Cells(3*Natoms+1:6*Natoms)
close(70)




end program test_checkCells
