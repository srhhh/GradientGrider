program getCells
use f1_parameters
use f1_variables
use f1_functions
use addCells3
implicit none
integer :: Ntraj,Nstates,state1,state2,state3,line_num,skips,i,j,k
logical :: flag1, flag2
real :: var1,var2,var3
integer :: var1_int, var2_int, var3_int
integer :: header1 = 1
integer :: header2 = 1
integer :: header3 = 1
integer, dimension(ceiling(max_var1*max_var2)) :: counter0 = 0
integer, dimension(counter1_max) :: counter1 = 0
integer, dimension(counter2_max) :: counter2 = 0
integer, dimension(counter3_max) :: counter3 = 0
integer, allocatable :: number_of_states(:,:)
real :: t1,t2
real,allocatable :: coords(:),vals(:)
character(50) :: current_path,line_data
character(20) :: subcell,trajectory_path,descriptor1,descriptor2


!Instead of printing to the terminal, we print to the progressfile
!First we check if it already exists and just wipe it
!We want to start from a clean slate
inquire(file=trim(path4)//trim(progressfile),exist=flag1)
if (flag1) call system("rm "//trim(path4)//trim(progressfile))
open(70,file=trim(path4)//trim(progressfile),status="new")
write(70,FMT="(A50)") "Let's start"
write(70,FMT="(A50)") ""
close(70)

!All trajectory folders are formatted as a number
!So search for these numbered folders and read them
call system("ls -p "//trim(path1)//" | grep '[0123456789]/' > "//trim(path4)//trajectories)
open(70,file=trim(path4)//trim(trajectories),action="read")

!coords  ---  this array keeps the coordinates and gradients of a state
!vals    ---  this array keeps var1, var2, var3 of a state
allocate(coords(6*Natoms))
allocate(vals(Nvar))

!The way the files are formatted, every seventh line is a new state
skips = Natoms+1

!We will keep track of how many trajectories and states we encounter
call CPU_time(t1)
Ntraj = 0
Nstates = 0
do

        if (.not. (start_from_scratch)) exit
!       if (Ntraj == 20) exit       !(if we wanted to end data collection prematurely)

        !Fetch the name of one folder (a trajectory)
        !Format its contents with sed into tmp.txt
        !If there are no more trajectories, iostat returns nonzero
        read(70,FMT="(A20)",iostat=state1) trajectory_path
        if (state1 /= 0) exit
        call system(trim(path2)//"trajectory_sed.sh "//trim(path1)//&
                trim(trajectory_path)//"kkk.out "//trim(path4)//trim(temporaryfile1))

        !A successful folder opening is a successful trajectory reading
        Ntraj = Ntraj+1



!To track the progress, open up 80
open(80,file=trim(path4)//trim(progressfile),position="append")
write(80,*) "Now accessing folder:", trajectory_path
close(80)



        !Open the now-formatted trajectory
        open(71,file=trim(path4)//trim(temporaryfile1))

        do
                !Get rid of the first line (just characters)
                !If it doesn't exists, we're at the end of the file
                read(71,FMT="(A50)", iostat=state2) line_data
                if (state2 /= 0) exit

                !Read the six lines of coordinates
                do line_num = 1, 6
                        i = 3*line_num
                        j = 3*Natoms
                        read(71,FMT="(10x,3(1x,F10.6),1x,3(1x,F10.6))") coords(i-2), coords(i-1), &
                                coords(i), coords(j+i-2), coords(j+i-1), coords(j+i)
                end do

                        !With the fully described state, calculate the
                        !variables wanted
                        call getVar1(coords(1:3*Natoms),Natoms,vals(1))
                        call getVar2(coords(1:3*Natoms),Natoms,vals(2))
                        call getVar3(coords(1:3*Natoms),Natoms,vals(3))

                        !If they are outliers, just skip this cycle
                        if ((vals(1) > max_var1).or.(vals(2) > max_var2)) then
                                Nstates = Nstates - 1
                                cycle
                        end if

                        !Write to the file the variables, coordinates, and gradients
                        call addState(vals,coords,&
                                header1,header2,header3,&
                                counter0,counter1,counter2,counter3)

                        !And this is one successful state/frame
                        Nstates = Nstates + 1
 
        end do              
        close(71)
end do
deallocate(vals,coords)
close(70)

call CPU_time(t2)
open(70,file=trim(path4)//trim(progressfile),position="append")
write(70,*), ""
write(70,*) "The reading and collecting took:", t2-t1, "seconds"
write(70,*) "The number of states is:", Nstates
close(70)
t1 = t2



open(70,file=trim(path4)//trim(counter0file))
do i = 1, ceiling(max_var1*max_var2)
        write(70,FMT="(I8)") counter0(i)
end do
close(70)

open(70,file=trim(path4)//trim(counter1file))
do i = 1, 250*resolution
        write(70,FMT="(I8)") counter1(i)
end do
close(70)

open(70,file=trim(path4)//trim(counter2file))
do i = 1, 500*resolution
        write(70,FMT="(I8)") counter2(i)
end do
close(70)

open(70,file=trim(path4)//trim(counter3file))
do i = 1, 250*resolution
        write(70,FMT="(I8)") counter3(i)
end do
close(70)




call CPU_time(t2)
open(70,file=trim(path4)//trim(progressfile),position="append")
write(70,*), ""
write(70,*) "The writing of counters took:", t2-t1, "seconds"
close(70)

end program getCells


