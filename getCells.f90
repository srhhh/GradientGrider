program getCells
use f1_parameters
use f1_variables
use f1_functions
use addCells4
implicit none
integer :: Ntraj,Nstates,Total_States,state1,state2,line_num,skips,i,gradient_index
logical :: flag1, flag2
real :: var1,var2,var3
integer :: var1_int, var2_int, var3_int
integer :: header1 = 1
integer :: header2 = 1
integer :: header3 = 1
integer, dimension(counter0_max) :: counter0 = 0
integer, dimension(counter1_max) :: counter1 = 0
integer, dimension(counter2_max) :: counter2 = 0
integer, dimension(counter3_max) :: counter3 = 0
integer, allocatable :: number_of_states(:,:)
real :: t1,t2,system_clock_rate
integer :: total_time1,total_time2,c1,c2,cr,cm
real,allocatable :: coords(:),vals(:)
character(50) :: current_path,line_data
character(20) :: subcell,trajectory_path,descriptor1,descriptor2


!Instead of printing to the terminal, we print to the progressfile
open(progresschannel,file=path4//progressfile)
write(progresschannel,FMT="(A50)") "Let's start"
write(progresschannel,FMT="(A50)") ""
close(progresschannel)

!All trajectory folders are formatted as a number
!So search for these numbered folders and read them
call system("ls -p "//path1//" | grep '[0123456789]/' > "//path4//trajectories)
open(trajectorieschannel,file=path4//trajectories,action="read")

!coords  ---  this array keeps the coordinates and gradients of a state
!vals    ---  this array keeps var1, var2, var3 of a state
allocate(coords(6*Natoms))
allocate(vals(Nvar))

gradient_index = 3*Natoms

!We want to know how fast this program takes
call system_clock(count_rate=cr)
call system_clock(count_max=cm)
system_clock_rate = real(cr)
call system_clock(total_time1)

!We will keep track of how many trajectories and states we encounter
Ntraj = 0
Total_States = 0
do

call system_clock(c1)

!       if (.not. (start_from_scratch)) exit
!       if (Ntraj == 10) exit       !(if we wanted to end data collection prematurely)

        !Fetch the name of one folder (a trajectory)
        !Format its contents with sed into tmp.txt
        !If there are no more trajectories, iostat returns nonzero
        read(trajectorieschannel,FMT="(A20)",iostat=state1) trajectory_path
        if (state1 /= 0) exit
        call system(path2//"trajectory_sed.sh "//path1//&
                trim(trajectory_path)//"kkk.out "//path4//temporaryfile1)

        !A successful folder opening is a successful trajectory reading
        Ntraj = Ntraj+1
        Nstates = 0

call system_clock(c2)
!To track the progress, open up progresschannel
open(progresschannel,file=path4//progressfile,position="append")
write(progresschannel,*) "Folder:", trajectory_path, " Formatting Time: ",(c2-c1)/system_clock_rate
close(progresschannel)



        !Open the now-formatted trajectory
        open(frameschannel,file=path4//temporaryfile1)

        do
                !Get rid of the first line (just characters)
                !If it doesn't exists, we're at the end of the file
                read(frameschannel,FMT="(A50)", iostat=state2) line_data
                if (state2 /= 0) exit

                !Read the six lines of coordinates
                do line_num = 1, 6
                        i = 3*line_num
                        read(frameschannel,FMT="(10x,3(1x,F10.6),1x,3(1x,F10.6))") coords(i-2), coords(i-1), &
                                coords(i), coords(gradient_index+i-2), &
                                coords(gradient_index+i-1),coords(gradient_index+i)
                end do

                        !With the fully described state, calculate the
                        !variables wanted
                        call getVar1(coords(1:3*Natoms),Natoms,vals(1))
                        call getVar2(coords(1:3*Natoms),Natoms,vals(2))
                        call getVar3(coords(1:3*Natoms),Natoms,vals(3))

                        !If they are outliers, just skip this cycle
                        if ((vals(1) > max_var1).or.(vals(2) > max_var2)) then
                                cycle
                        end if

                        !Write to the file the variables, coordinates, and gradients
                        call addState(vals,coords,&
                                header1,header2,header3,&
                                counter0,counter1,counter2,counter3)

                        !And this is one successful state/frame
                        Nstates = Nstates + 1
 
        end do              
        close(frameschannel)

call system_clock(c2)
open(progresschannel,file=path4//progressfile,position="append")
write(progresschannel,*) "Total Time:", (c2-c1)/system_clock_rate, "   Frames: ",Nstates
write(progresschannel,*) ""
close(progresschannel)
Total_States = Total_States + Nstates

end do
deallocate(vals,coords)
close(trajectorieschannel)

call system_clock(total_time2)

open(progresschannel,file=path4//progressfile,position="append")
write(progresschannel,*), ""
write(progresschannel,*), ""
write(progresschannel,*) "The reading and collecting took: ", &
                         (total_time2-total_time1)/system_clock_rate, " System time"
write(progresschannel,*) "Total number of frames: ", Total_States
close(progresschannel)
total_time1 = total_time2


open(filechannel1,file=path4//counter0file)
do i = 1, counter0_max
        write(filechannel1,FMT="(I8)") counter0(i)
end do
close(filechannel1)

open(filechannel1,file=path4//counter1file)
do i = 1, counter1_max
        write(filechannel1,FMT="(I8)") counter1(i)
end do
close(filechannel1)

open(filechannel1,file=path4//counter2file)
do i = 1, counter2_max
        write(filechannel1,FMT="(I8)") counter2(i)
end do
close(filechannel1)

open(filechannel1,file=path4//counter3file)
do i = 1, counter3_max
        write(filechannel1,FMT="(I8)") counter3(i)
end do
close(filechannel1)





call system_clock(total_time2)

open(progresschannel,file=path4//progressfile,position="append")
write(progresschannel,*), ""
write(progresschannel,*) "The writing of counters took: ", &
                         (total_time2-total_time1)/system_clock_rate, " System time"
close(progresschannel)

end program getCells


