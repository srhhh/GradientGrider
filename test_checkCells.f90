program test_checkCells
use checkCells3
use addCells3
use ls_rmsd
use f1_variables
use f1_parameters

implicit none
real, dimension(Nvar) :: vals
real, dimension(3*Natoms) :: coords, gradient
real, dimension(6*Natoms) :: coords_Cells
double precision :: min_rmsd
integer :: i, j, k, rand_subcell_position, rand_file_position, state1
logical :: flag1
character(50) :: subcell

rand_subcell_position = 29
rand_file_position = 12


!Instead of printing to the terminal, we print to the progress files
!First we check if it exists and just wipe
!We want ot start from a clean slate
inquire(file=trim(path4)//trim(progressfile),exist=flag1)
if (flag1) call system("rm "//trim(path4)//trim(progressfile))
open(70,file=trim(path4)//trim(progressfile),status="new")
write(70,FMT="(A50)") "Let's start"
write(70,FMT="(A50)") ""
close(70)

!We assume that the grid already exists
!We take a random frame in the grid, change it slightly, and
!see what frame is closest to it
call system("ls -p "//trim(path3)//" | grep '.dat' > "//trim(path4)//trim(temporaryfile1))
open(70,file=trim(path4)//trim(temporaryfile1))
do i = 1, rand_subcell_position
        read(70,FMT="(A50)") subcell
end do
close(70)

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



!This is how I choose to change the coordinates
coords(1) = coords(1) + 0.01
coords(4) = coords(4) - 0.004
coords(18) = coords(18) + .005

open(70,file=trim(path4)//trim(progressfile),position="append")
write(70,FMT="(A50)") "Making these changes to it:"
write(70,FMT="(A50)") ""
write(70,*) "coords(1) = coords(1) + 0.01"
write(70,*) "coords(4) = coords(4) - 0.004"
write(70,*) "coords(18) = coords(18) + .005"
write(70,FMT="(A50)") ""
close(70)




!Let's see how it works out
call checkState(coords,coords_Cells,min_rmsd)
open(70,file=trim(path4)//trim(progressfile),position="append")
write(70,*) "The original coordinates are: "
write(70,*) coords
write(70,*) ""
write(70,*) "The closests coordinates are: "
write(70,*) coords_Cells(1:3*Natoms)
write(70,*) ""
write(70,*) "With variables: ", vals(1), " A and ", vals(2), " A"
write(70,*) ""
write(70,*) "This has an rmsd of: ", min_rmsd
write(70,*) ""
write(70,*) "The gradient is: "
write(70,*) coords_Cells(3*Natoms+1:6*Natoms)
close(70)




end program test_checkCells
