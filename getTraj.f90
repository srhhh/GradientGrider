
program getTraj
use f1_functions
implicit none
integer :: Ntraj,Nlength,Nlength1,Nvar,Nstates,state1,state2,state3,line_num,Natoms,m_i,skips,i,j,k,i1,i2
integer :: Ngrid1,Ngrid2,Ngrid3,j1,j2
logical :: flag1, flag2
real :: x,y,z,val1,val2,val3
real :: var1,var2,var3,spacing1,spacing2,spacing3,max_var1,max_var2,max_var3
real :: t1,t2
real,allocatable :: coords(:,:),vals(:,:),temp_vals(:,:),temp_coords(:,:)
integer,allocatable :: indexer(:,:),temp_indexer(:,:),grid1(:),grid2(:),grid3(:),grid_total(:)
character(50) :: path,line_data
character(20) :: file1,file2,some_folder



!Path to the files
path = "/home/ruisun/proj/ruisun/B0/"

!File that will keep the trajectory folder names
file1 = "f1_trajectories.txt"

!File that will keep the var1, var2, var3, coordinates, gradients
file2 = "f1_data_temp.dat"

call system("ls -p " // trim(path) // " | grep '[0123456789]/' > " // file1)
open(70,file=trim(file1),action="read")
open(72,file=trim(file2),action="write")

!Some parameters, Nlength is the initial guess for how many states; this can
!grow later
Natoms = 6
Nlength = 10000
Nvar = 3

!The spacing is the spacing for the gridlines; currently var3 is not grided
spacing1 = 1.0
spacing2 = 1.0
spacing3 = 0.02

!The scaling is the amount that is resolved for a variable when there is
!overcrowding (0.1,0.1,0.1) = x1000 magnification
scaling1 = 0.10
scaling2 = 0.10
scaling3 = 0.10

!There are some outliers; making a maximum throws these away
max_var1 = 50.0
max_var2 = 50.0



!coords  ---  this array keeps the coordinates and gradients, and is UNSORTED
!vals    ---  this array keeps var1, var2, var3 and GETS SORTED
!indexer ---  this array points towards elemnts of coords and also GETS SORTED
allocate(coords(Nlength,6*Natoms))
allocate(vals(Nlength,Nvar),indexer(Nlength,1))
skips = Natoms+1
Nlength1 = Nlength - 1


!We will keep track of how many trajectories and states we encounter
call CPU_time(t1)
Ntraj = 0
Nstates = 0
do
        !Fetch the name of one folder (a trajectory)
        read(70,FMT="(A20)",iostat=state1) some_folder
        if (state1 /= 0) exit
        print *, "Now accessing folder:", some_folder
        Ntraj = Ntraj+1
        
        !Open the file in the folder which has the trajectories (the fort.8
        !file)
        open(71,file=trim(path)//trim(some_folder)//"kkk.out")

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
                                call getVar1(coords(Nstates,:),Natoms,var1)
                                call getVar2(coords(Nstates,:),Natoms,var2)
                                call getVar3(coords(Nstates,:),Natoms,var3)

                                !This increases the array length when it is
                                !about to overflow (seems inefficient to me)
                                !Also: I could not figure out how to put this in
                                !a subroutine... some allocation error
                                if (Nstates > Nlength1) then
                                        print *, "   Reallocating...", Nlength
                                        Nlength = Nlength*2
                                        Nlength1 = Nlength-1
                                        allocate(temp_vals(Nlength,Nvar),temp_indexer(Nlength,1),temp_coords(Nlength,6*Natoms))
                                        call growR(vals,temp_vals,Nlength/2,Nvar,Nlength,Nvar)
                                        call growR(coords,temp_coords,Nlength/2,6*Natoms,Nlength,6*Natoms)
                                        call growI(indexer,temp_indexer,Nlength/2,1,Nlength,1)
                                        deallocate(vals,indexer,coords)
                                        allocate(vals(Nlength,Nvar),indexer(Nlength,1),coords(Nlength,6*Natoms))
                                        vals = temp_vals
                                        indexer = temp_indexer
                                        coords = temp_coords
                                        deallocate(temp_vals,temp_indexer,temp_coords)
                                        print *, "   Reallocated to:", Nlength
                                end if

                                !We put the variables in vals and indices that
                                !point to their corresponding state in indexer
                                !In the beginning, they are unordered
                                vals(Nstates,1) = var1
                                vals(Nstates,2) = var2
                                vals(Nstates,3) = var3

                                !If they are outliers, just skip this cycle
                                if ((var1 > max_var1).or.(var2 > max_var2)) then
                                        Nstates = Nstates - 1
                                else
                                        indexer(Nstates,1) = Nstates
                                end if

                                !And repeat
                                flag2 = .true.
                                cycle
                        end if

                        !If the state is not yet fully described, continue
                        !adding coordinates; (line_num keeps track of atom #)
                        state3 = 1
                        call read_coords(coords,Nlength,Natoms,Nstates,line_num,state3)

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
        if (Ntraj == 2) exit
end do
call CPU_time(t2)
print *, "The reading and collecting took:", t2-t1, "seconds"
print *, "The number of states is:", Nstates
t1 = t2


print *, ""
print *, "minimum r1:", minval(vals(1:Nstates,1)), "at", minloc(vals(1:Nstates,1))
!print *, coords(minloc(vals(1:Nstates,1)),:)
print *, "minimum r2:", minval(vals(1:Nstates,2)), "at", minloc(vals(1:Nstates,2))
!print *, coords(minloc(vals(1:Nstates,2)),:)
print *, ""
print *, "maximum r1:", maxval(vals(1:Nstates,1)), "at", maxloc(vals(1:Nstates,1))
print *, "maximum r2:", maxval(vals(1:Nstates,2)), "at", maxloc(vals(1:Nstates,2))
print *, ""







!The first sorting doesn't require a do loop
!We can also gauge how much time it takes to sort a list of this size
call qsort(vals,indexer,Nlength,Nvar,1,Nstates,1)
call CPU_time(t2)
print *, "The broad sorting took:", t2-t1, "seconds"
print *, ""
t1 = t2






!This is the griding for the first variable
Ngrid1 = ceiling(maxval(vals(:,1))/spacing1)
print *, "The number of gridlines for var1:", Ngrid1
allocate(grid1(Ngrid1))
call grider(grid1,vals,spacing1,Ngrid1,Nlength,Nvar,1,Nstates,1,Ngrid1,1)
call CPU_time(t2)
print *, "The first griding took:", t2-t1, "seconds"
print *, ""
t1 = t2





!This is the griding for the second variable
Ngrid2 = ceiling(maxval(vals(:,2))/spacing2)
print *, "The number of gridlines for var2:", Ngrid2
allocate(grid2(Ngrid1*Ngrid2))
do i1 = 1, Ngrid1
        j1 = grid1(i1)
        if (i1 == Ngrid1) then
                j2 = Nstates+1
        else
                j2 = grid1(i1+1)
        end if
        if (j1 == j2) then
                do i = 1, Ngrid2
                        grid2(Ngrid2*(i1-1)+i) = j1
                end do
        else
                call qsort(vals,indexer,Nlength,Nvar,j1,j2-1,2)
                call grider(grid2,vals,spacing2,Ngrid1*Ngrid2,Nlength,Nvar,j1,j2-1,Ngrid2*(i1-1)+1,Ngrid2*i1,2)
        end if
                if (.false.) then
                        print *, grid2(Ngrid2*(i1-1)+1:Ngrid2*i1)
                end if

end do
print *, "The number of states is:", Nstates
call CPU_time(t2)
print *, "The second-level griding took:", t2-t1, "seconds"
print *, ""
t1 = t2
!I did not yet implement griding for the third variable







!Finally, we write to the files the organized list
do i = 1, Nstates
        j = indexer(i,1)
        write(72,FMT="(3(1x,F11.6))",advance="no") (vals(i,m_i),m_i=1,Nvar)
        write(72,FMT="(36(1x,F11.6))")(coords(j,m_i),m_i=1,6*Natoms)
end do
call CPU_time(t2)
print *, "The final writing took:", t2-t1, "seconds"






close(70)
close(71)
close(72)
deallocate(vals,indexer,coords)

end program getTraj








!Variable one is the squared distance between atom 1 (fluorine) and 2 (carbon)
subroutine getVar1(coords,Natoms,var1)
implicit none
integer, intent(in) :: Natoms
real, dimension(6*Natoms), intent(in) :: coords
real, intent(out) :: var1

var1 =sqrt( (coords(7)-coords(1))**2 + (coords(8)-coords(2))**2 + (coords(9)-coords(3))**2)
end subroutine getVar1

!Variable two is the squared distance between atom 2 (carbon) and 6 (iodine)
subroutine getVar2(coords,Natoms,var2)
implicit none
integer, intent(in) :: Natoms
real, dimension(6*Natoms), intent(in) :: coords
real, intent(out) :: var2

var2 = sqrt((coords(31)-coords(7))**2 + (coords(32)-coords(8))**2 + (coords(33)-coords(9))**2)
end subroutine getVar2

!Variable three is the cosine of the angle between atoms 1 - 2 - 6
subroutine getVar3(coords,Natoms,var3)
implicit none
integer, intent(in) :: Natoms
real, dimension(6*Natoms), intent(in) :: coords
real, intent(out) :: var3
real :: x1,y1,z1,x2,y2,z2

x1 = coords(31)-coords(7)
y1 = coords(32)-coords(8)
z1 = coords(33)-coords(9)
x2 = coords(1)-coords(7)
y2 = coords(2)-coords(8)
z2 = coords(3)-coords(9)
var3 = (x1*x2+y1*y2+z1*z2)/sqrt((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))
end subroutine getVar3








!This exists because some z coordinates are ******
subroutine read_coords(coords,Nlength,Natoms,Nstates,line_num,stat)
implicit none
integer, intent(in) :: Nlength,Natoms, Nstates, line_num
integer, intent(out) :: stat
character(11) :: cvar
real, dimension(Nlength,6*Natoms), intent(out) :: coords

read(71,FMT="(10x,2(1x,F10.6),(A11),1x,3(1x,F10.6)))") coords(Nstates,6*line_num-5), &
& coords(Nstates,6*line_num-4),cvar,coords(Nstates,6*line_num-2), &
& coords(Nstates,6*line_num-1),coords(Nstates,6*line_num)

if (cvar /= " **********") then
        read (cvar, FMT="(1x,F10.6)") coords(Nstates,6*line_num-3)
        stat = 0
end if
end subroutine






end program getTraj
