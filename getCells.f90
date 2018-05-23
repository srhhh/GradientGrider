
program getCells
use f1_parameters
use f1_variables
use f1_functions
use addCells
implicit none
integer :: Ntraj,Nstates,state1,state2,state3,line_num,skips,i,j,k
logical :: flag1, flag2
real :: var1,var2,var3
integer :: var1_int, var2_int, var3_int
real :: t1,t2
real,allocatable :: coords(:),vals(:)
character(50) :: line_data
character(20) :: file2,file3,some_folder,descriptor1,descriptor2


file3 = "progress.txt"
inquire(file=trim(path2)//trim(file3),exist=flag1)
if (flag1) call system("rm "//trim(path2)//trim(file3))
open(70,file=trim(path2)//trim(file3),status="new")
close(70)
call system("ls -p " // trim(path1) // " | grep '[0123456789]/' > " // file1)
open(70,file=trim(file1),action="read")




!coords  ---  this array keeps the coordinates and gradients of a state
!vals    ---  this array keeps var1, var2, var3 of a state
allocate(coords(6*Natoms))
allocate(vals(Nvar))
skips = Natoms+1


!We will keep track of how many trajectories and states we encounter
call CPU_time(t1)
Ntraj = 0
Nstates = 0
do

 exit    !Comment this out later

        !Fetch the name of one folder (a trajectory)
        read(70,FMT="(A20)",iostat=state1) some_folder
        if (state1 /= 0) exit
open(80,file=trim(path2)//trim(file3),position="append")
write(80,*) "Now accessing folder:", some_folder
close(80)
        Ntraj = Ntraj+1
        
        !Open the file in the folder which has the trajectories (the fort.8
        !file)
        open(71,file=trim(path1)//trim(some_folder)//"kkk.out")

        !flag2 indicates when coordinates and gradients are NOT on the next line
        !All following lines will be in a pattern that can be read
        flag2 = .true.
        do
                if (flag2) then
                        !If the line contains "     x         y    ", this
                        !indicates a state coming up
                        read(71,FMT="(A50)", iostat=state2) line_data
                        if ((flag2).and.(index(line_data,"x          y")/= 0)) then
                                flag2 = .false.
                                Nstates = Nstates + 1
                                line_num = 0
                        end if
                else
                        !We only have a fully-described state when all atoms
                        !have their coordinate and gradient written
                        line_num = line_num + 1
                        if (line_num == skips) then
                                !With the fully described state, calculate the
                                !variables wanted
                                call getVar1(coords(1:3*Natoms),Natoms,var1)
                                call getVar2(coords(1:3*Natoms),Natoms,var2)
                                call getVar3(coords(1:3*Natoms),Natoms,var3)

                                !If they are outliers, just skip this cycle
                                if ((var1 > max_var1).or.(var2 > max_var2)) then
                                        Nstates = Nstates - 1
                                        flag2 = .true.
                                        cycle
                                end if

                                vals(1) = var1
                                vals(2) = var2
                                vals(3) = var3

                                !Find which appropriate cell the state is in
                                var1_int = floor(var1*spacing1)
                                var2_int = floor(var2*spacing2)
                                var3_int = floor(var3*spacing3)

                                !Make the filename for the cell
                                write(descriptor1,FMT="(I4)") var1_int
                                write(descriptor2,FMT="(I4)") var2_int
                                file2 = trim(adjustl(descriptor1))//"_"//trim(adjustl(descriptor2))//".dat"
                                inquire(file=trim(path3)//trim(file2),exist=flag1)


                                !Write to the file the variables, coordinates,
                                !and gradients
                                if (flag1) then
                                        open(72,file=trim(path3)//trim(file2),position="append")
                                        write(72,FMT="(3(1x,F11.6))",advance="no") (vals(i),i=1,Nvar)
                                        write(72,FMT="(36(1x,F11.6))")(coords(i),i=1,6*Natoms)
                                        close(72)
                                else
                                        open(72,file=trim(path3)//trim(file2),position="append",status="new")
                                        write(72,FMT="(3(1x,F11.6))",advance="no") (vals(i),i=1,Nvar)
                                        write(72,FMT="(36(1x,F11.6))")(coords(i),i=1,6*Natoms)
                                        close(72)
                                end if

                                flag2 = .true.
                                cycle
                        end if

                        !If the state is not yet fully described, continue
                        !adding coordinates; (line_num keeps track of atom #)
                        state3 = 1
                        call read_coords(coords,Natoms,line_num,state3)

                        !If there is a z coordinate equal to ******, skip it
                        if (state3 /= 0) then
                                read(71,FMT="(A50)") line_data
                                Nstates = Nstates - 1
                                flag2 = .true.
                        end if
                end if
                

                !At the end of the file, quit
                if (state2 /= 0) exit
        end do              
        close(71)
!       if (Ntraj == 20) exit
end do
call CPU_time(t2)
close(70)

open(70,file=trim(path2)//trim(file3),position="append")
write(70,*) "The reading and collecting took:", t2-t1, "seconds"
write(70,*) "The number of states is:", Nstates
close(70)
t1 = t2

deallocate(vals,coords)







!Now we just need to organize every folder
!Note: probably more efficient do this inside the loop
!maybe by using the addState subroutine instead
line_data = ""
do i = 1, floor(max_var1/spacing1)-1
do j = 1, floor(max_var2/spacing1)-1

open(70,file=trim(path2)//trim(file3),position="append")
write(70,*) "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
write(70,*) ""
write(70,*) "divying up   var1:", i, "var2:", j
write(70,*) ""
close(70)

call divyUp(i,j,1,0,line_data)


end do

end do







end program getCells







!This exists because some z coordinates are ******
subroutine read_coords(coords,Natoms,line_num,stat)
implicit none
integer, intent(in) :: Natoms, line_num
integer, intent(out) :: stat
character(11) :: cvar
real, dimension(6*Natoms), intent(out) :: coords

read(71,FMT="(10x,2(1x,F10.6),(A11),1x,3(1x,F10.6)))") coords(3*line_num-2), &
& coords(3*line_num-1),cvar,coords(3*Natoms+3*line_num-2), &
& coords(3*Natoms+3*line_num-1),coords(3*Natoms+3*line_num)

if (cvar /= " **********") then
        read (cvar, FMT="(1x,F10.6)") coords(3*line_num)
        stat = 0
end if
end subroutine






