
program getGrid
use f1_parameters
use f1_functions
use addCells
implicit none
integer :: Ntraj,Nlength,Nlength1,Nvar,Nstates,state1,state2,line_num,skips,i,j,k
logical :: flag1
real :: var1,var2,var3
real :: t1,t2
real,allocatable :: coords(:,:),vals(:,:)!,temp_vals(:,:),temp_coords(:,:)
integer,allocatable :: indexer(:,:)!,temp_indexer(:,:),grid1(:),grid2(:),grid3(:),grid_total(:)
character(50) :: line_data
character(20) :: some_folder



call system("ls -p "//trim(path1)//" | grep '[0123456789]/' > "//file1)
open(70,file=trim(path2)//trim(file1),action="read")

!Nlength is the initial guess for how many states; this can grow later
Nlength = 10000


!coords  ---  this array keeps the coordinates and gradients, and is SORTED
!vals    ---  this array keeps var1, var2, var3 and is ASSUMED SORTED
allocate(coords(Nlength,6*Natoms))
allocate(vals(Nlength,Nvar))
allocate(indexer(Nlength,1))
Nlength1 = Nlength - 1
skips = Natoms+1

!We will keep track of how many states and trajectories we encounter
call CPU_time(t1)
Ntraj = 0
Nstates = 0
do
        !Fetch the name of one folder (a trajectory)
        read(70,FMT="(A50)",iostat=state1) some_folder
        if (state1 /= 0) exit
        print *, "Now accessing folder:", some_folder
        Ntraj = Ntraj+1
        open(71,file=trim(path1)//trim(some_folder)//"kkk.out"
        flag1 = .true.

do
if (flag1) then
        read(71,FMT="(A50)", iostat=state2) line_data
        if ((flag1).and.(index(line_data,"x          y")/=0) then
                flag1 = .false.
                Nstates = Nstates + 1
                line_num = 0
        end if
else
line_num = line_num + 1
if (line_num == skips) then
        call getVar1(coords,Natoms,var1)
        call getVar2(coords,Natoms,var2)
        call getVar3(coords,Natoms,var3)

        !This increases the array length when it is
        !about to overflow (seems inefficient to me)
        !Also: I could not figure out how to put this in
        !a subroutine... some allocation error
!       if (Nstates > Nlength1) then
!               print *, "   Reallocating...", Nlength
!               Nlength = Nlength*2
!               Nlength1 = Nlength-1
!               allocate(temp_vals(Nlength,Nvar),temp_coords(Nlength,6*Natoms),temp_indexer(Nlength,1))
!               call growR(vals,temp_vals,Nlength/2,Nvar,Nlength,Nvar)
!               call growR(coords,temp_coords,Nlength/2,6*Natoms,Nlength,6*Natoms)
!               call growI(indexer,temp_indexer,Nlength/2,1,Nlength,1)
!               deallocate(vals,indexer,coords)
!               allocate(vals(Nlength,Nvar),coords(Nlength,6*Natoms),indexer(Nlength,1))
!               vals = temp_vals
!               coords = temp_coords
!               indexer = temp_indexer
!               deallocate(temp_vals,temp_indexer,temp_coords)
!               print *, "   Reallocated to:", Nlength
!       end if

        !If they are outliers, just skip this cycle
        if ((var1 > max_var1).or.(var2 > max_var2)) then
                Nstates = Nstates - 1
                flag1 = .true.
                cycle
        end if


        !We put the variables in vals and indices that
        !point to their corresponding state in indexer
        !In the beginning, they are unordered
        vals(Nstates,1) = var1
        vals(Nstates,2) = var2
        vals(Nstates,3) = var3

!       indexer(Nstates,1) = Nstates
        
        !Add the state
        call addState(var1,var2,var3,vals,coords)

        !And repeat
        flag2 = .true.
        cyc1e
end if



!If the state is not yet fully described, continue
!adding coordinates; (line_num keeps track of atom #)
state3 = 1
call read_coords(coords,Nlength,Natoms,Nstates,line_num,state3)

!If there is a z coordinate equal to ******, skip it
if (state3 /= 0) then
        read(71,FMT="(A50)") line_data
        Nstates = Nstates - 1
        flag1 = .true.
end if

!At the end of the file, quit
if (state2 /= 0) exit
end do              
close(71)



!If we want to end this program early, uncomment this
if (Ntraj == 2) exit
end do

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
!call qsort(vals,indexer,Nlength,Nvar,1,Nstates,1)
!call CPU_time(t2)
!print *, "The broad sorting took:", t2-t1, "seconds"
!print *, ""
!t1 = t2






!This is the griding for the first variable
!Ngrid1 = ceiling(maxval(vals(:,1))/spacing1)
!print *, "The number of gridlines for var1:", Ngrid1
!allocate(grid1(Ngrid1))
!call grider(grid1,vals,spacing1,Ngrid1,Nlength,Nvar,1,Nstates,1,Ngrid1,1)
!call CPU_time(t2)
!print *, "The first griding took:", t2-t1, "seconds"
!print *, ""
!t1 = t2





!This is the griding for the second variable
!Ngrid2 = ceiling(maxval(vals(:,2))/spacing2)
!print *, "The number of gridlines for var2:", Ngrid2
!allocate(grid2(Ngrid1*Ngrid2))
!do i1 = 1, Ngrid1
!        j1 = grid1(i1)
!        if (i1 == Ngrid1) then
!                j2 = Nstates+1
!        else
!                j2 = grid1(i1+1)
!        end if
!        if (j1 == j2) then
!                do i = 1, Ngrid2
!                        grid2(Ngrid2*(i1-1)+i) = j1
!                end do
!        else
!                call qsort(vals,indexer,Nlength,Nvar,j1,j2-1,2)
!                call grider(grid2,vals,spacing2,Ngrid1*Ngrid2,Nlength,Nvar,j1,j2-1,Ngrid2*(i1-1)+1,Ngrid2*i1,2)
!        end if
!                if (.false.) then
!                        print *, grid2(Ngrid2*(i1-1)+1:Ngrid2*i1)
!                end if
!
!end do
!print *, "The number of states is:", Nstates
!call CPU_time(t2)
!print *, "The second-level griding took:", t2-t1, "seconds"
!print *, ""
!t1 = t2
!I did not yet implement griding for the third variable













close(70)
close(71)
deallocate(vals,indexer,coords)








end program getGrid
