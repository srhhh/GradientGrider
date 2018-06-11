program test_checkCells
use checkCells4
use f1_variables
use f1_parameters

implicit none
real, dimension(Nvar) :: vals
! RS: I suggest we stay with the same set of data for reading and storing
! RS: i.e. be consistant with (3*Natoms) :: coords, gradient
! RS:      or (6*Natoms) :: coords (as you have in addCells)
!                               KF: in this particular case, we may not have the
!                               gradient for coords, so I use only 3*Natoms;
!                               but we could easily increase it to 6*Natoms if
!                               it may be more convenient later
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

! RS: Why would we want to rid of this progress file? 
! RS: I think we should treat this as a two-step process
! RS: 1. construct the data (when there is not much data it does not make much sense to check in the database)
! RS: 2. checking the data and maybe write more into the data
!                       KF: I just copy-pasted those line from getCells....
!                       (this was when I was just starting the project)
!                       in reality, I could replace all those lines with:

!Instead of printing to the terminal, we print to the progress files
open(progresschannel,file=trim(path4)//trim(progressfile))
write(progresschannel,FMT="(A50)") "Let's start"
write(progresschannel,FMT="(A50)") ""
close(progresschannel)

!                       KF: this is more compact
!                       remark: this is just so that I don't have to print
!                               and then look at the terminal output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          COUNTER ARRAYS RETRIEVAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!remark: to TEST the grid, some amount of grid (the folder?) has to actually have been MADE beforehand
!        if you want to quickly make the grid, do into getCells.f90
!        and make the program quit after maybe twenty trajectories.
!        This program will work regardless

! counter'X' is the array storing the counter (# of states in a cell/subcell) and pointer. 
! Pionter is used in case where a cell/subcell get overcrowded, pointer points to the 'address' 
!    of the counter and pointer in the children-level array counter'X+1'. 
! X is the order of the cell/subcell, with 0 corresponding to the parent-level cell. 
!    e.g. 1 corresponding to the first parent-level cell, etc.
! 
! In the counter'X' array, every element has the type of integer*8, in which the first four digits 
!    represent the 'pointer', and the last four digits represent the 'counter' (# of states in a cell/subcell)
! For example,
!    counter0(i) = aaaa|bbbb   (the pipleline is imaginary) 
!    if aaaa = 0, it indicates the ith cell (X=0) has bbbb states stored
!    if aaaa > 0, the ith cell is overcrowded and the cell has been divided (through divyup in addCells)
!        and the data has been stored in subcell. The counter and pointer can be found at counter1 (X=1)
!
!                       KF: yes, I should add something like this to the README

open(filechannel1,file=trim(path4)//"counter0.txt",status="old")
do i = 1, counter0_max
        read(filechannel1,FMT="(I8)") counter0(i)
end do
close(filechannel1)

open(filechannel1,file=trim(path4)//"counter1.txt",status="old")
do i = 1, counter1_max
        read(filechannel1,FMT="(I8)") counter1(i)
end do
close(filechannel1)

open(filechannel1,file=trim(path4)//"counter2.txt",status="old")
do i = 1, counter2_max
        read(filechannel1,FMT="(I8)") counter2(i)
end do
close(filechannel1)

open(filechannel1,file=trim(path4)//"counter3.txt",status="old")
do i = 1, counter3_max
        read(filechannel1,FMT="(I8)") counter3(i)
end do
close(filechannel1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          RANDOM FRAME RETRIEVAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rand_subcell_position = 29
rand_file_position = 12

! Get a list of all cells and subcells; pick a random one
call system("ls -p "//trim(path3)//" | grep '.dat' > "//trim(path4)//trim(temporaryfile1))
open(trajectorieschannel,file=trim(path4)//trim(temporaryfile1))
! RS: There is an option in OPEN call RECL that allows you to directly read from a specific line
! RS: https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnaf/index.html
!                               KF: I don't believe you just yet
do i = 1, rand_subcell_position
        read(trajectorieschannel,FMT="(A50)") subcell
end do
close(trajectorieschannel)

!writing the checking information into the progress file
open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
write(progresschannel,FMT="(A50)") "For the random frame, looking in ", trim(subcell)
write(progresschannel,FMT="(A50)") ""
close(progresschannel)

! RS: One good practice of file I/O is to use different indexer for the file
! RS: i.e. you keep using progresschannel for various files, which has two potential disadvantages:
! RS:    a. Hard to keep track of what progresschannel if the operation is long
! RS:    b. Hard to debug if you accidential forgot to close the progresschannel before
!                               KF: maybe set a variable. for example:
!                               integer, parameter :: progresschannel = progresschannel
!                               integer, parameter :: filechannel1 = 71
!                                       etc

!Open up that subcell; pick a random frame
open(frameschannel,file=trim(path3)//trim(subcell))
i = 0
do
        i = i + 1
        if (i < rand_file_position) then
                read(frameschannel,FMT="(A50)",iostat=state1) subcell
                if (state1 /= 0) then
                        i = i - 1
! RS: Did you use rewind here so that you won't get an error if rand_file_position is larger
! RS:    than the cell/subcell that you are in?
! RS: I think we should be considerate about that case and just add an exit with some writing
! RS:    statement, e.g. "the frame that you are trying to check does not exist"
!                                       KF: yes... it is inefficient but it
!                                       doesn't seem that important (just a
!                                       test); I wanted as little friction as
!                                       possible from number input to frame
!                                       output
                        rewind(frameschannel)
                end if
        else
                read(frameschannel,FMT=FMT1,advance="no",iostat=state1) (vals(j),j=1,Nvar)
                if (state1 /= 0) then
                        i = i - 1
                        close(frameschannel)
                        open(frameschannel,file=trim(path3)//trim(subcell))
                        cycle
                end if
                read(frameschannel,FMT=FMT2) (coords(j),j=1,3*Natoms)
                exit 
        end if
end do
close(frameschannel)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               FRAME MODIFICATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
write(progresschannel,*) "The original coordinates are: "
write(progresschannel,*) coords
write(progresschannel,*) ""
write(progresschannel,*) "With variables: ", vals(1), " A and ", vals(2), " A"
write(progresschannel,*) ""
close(progresschannel)

!This is how I choose to change the coordinates

! RS: It is fine that you choose the following but why we are doing this at the first place?
! RS: We do we want to chagne the coord?
!                       KF: This part is just to play with;
!                       Because the frame is already part of the grid, if we
!                       don't make any changes, we'll just get the same frame
!                       back (boring). This change demonstrates a new frame will
!                       produce a frame from the grid that approximates it
coords(3) = coords(3) + 0.1
coords(4) = coords(4) - 0.000
coords(18) = coords(18) + .000

open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
write(progresschannel,FMT="(A50)") "Making these changes to it:"
write(progresschannel,FMT="(A50)") ""
write(progresschannel,*) "coords(3) = coords(3) + 0.1"
write(progresschannel,*) "coords(4) = coords(4) - 0.000"
write(progresschannel,*) "coords(18) = coords(18) + .000"
write(progresschannel,FMT="(A50)") ""
close(progresschannel)

call getVar1(coords,Natoms,vals(1))
call getVar2(coords,Natoms,vals(2))

open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
write(progresschannel,*) "The modified coordinates are: "
write(progresschannel,*) coords
write(progresschannel,*) ""
write(progresschannel,*) "With variables: ", vals(1), " A and ", vals(2), " A"
write(progresschannel,*) ""
close(progresschannel)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               GRID CHECKING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RS: To check weather the modified grid exsit?
!                                       KF: Everything up until here was just to
!                                       get a new frame;
!                                       This subroutine call actually goes into
!                                       the grid and checks for a frame that is
!                                       close to it
call checkState(coords,coords_Cells,min_rmsd,&
                counter0,counter1,counter2,counter3)

call getVar1(coords_Cells,Natoms,vals(1))
call getVar2(coords_Cells,Natoms,vals(2))

open(progresschannel,file=trim(path4)//trim(progressfile),position="append")
write(progresschannel,*) "The closests coordinates are: "
write(progresschannel,*) coords_Cells(1:3*Natoms)
write(progresschannel,*) ""
write(progresschannel,*) "With variables: ", vals(1), " A and ", vals(2), " A"
write(progresschannel,*) ""
write(progresschannel,*) ""
write(progresschannel,*) "This has an rmsd of: ", min_rmsd
write(progresschannel,*) ""
write(progresschannel,*) ""
write(progresschannel,*) "The closest frame has a gradient as follows: "
write(progresschannel,*) coords_Cells(3*Natoms+1:6*Natoms)
close(progresschannel)




end program test_checkCells
