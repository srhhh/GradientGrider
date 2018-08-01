module analyzeHeatMapswithMultipleGrids
implicit none



contains

subroutine analyzeHeatMaps1(Ngrid_cap)
use PARAMETERS
use mapCellData
implicit none

!CAP TO NUMBER OF GRIDS ANALYZED
integer, intent(in) :: Ngrid_cap
integer :: Ngrid_total

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(6) :: Ntraj_text

!I/O HANDLING
integer :: iostate

!INTEGER INCREMENTALS
integer :: n


!All trajectory folders are formatted as I0.3 (3-digit integer)
!So search for these numbered folders and read them
call system("ls -p "//gridpath0//" | grep '[0123456789]/' > "//gridpath0//trajectories)

!This is a slight 'hack'; for the 'true' scattering angle plots based on
!the trajectories FROM THE GRID, we only need to concatenate
!all the trajectories and read off a specific column
!To do this, we make a long string " file1 file2 file3 ... fileN"
!And only take whatever portion out we want (ex. " file1 file2")
open(trajectorieschannel,file=gridpath0//trajectories,action="read")
do Ngrid_total = 1, Ngrid_max
        read(trajectorieschannel,FMT="(A4)",iostat=iostate) folder_text
        if (iostate /= 0) exit
end do
close(trajectorieschannel)



!We increment before we read, so we need to subtract out one increment here
if (Ngrid_total < 2) return
Ngrid_total = Ngrid_total - 1

!For now, we are putting a soft cap on how many trajectories we are using
!Each grid has Ntraj_max trajectories so we will use a maximum of 4 * Ntraj_max trajectories
Ngrid_total = min(4, Ngrid_total)

!First, we do the 'true' scattering angle plots
!This data was made during creation (or should have been!) so all we need to do
!is merge, read, and plot them
write(variable_length_text,FMT="(I5)") Ngrid_text_length
do Ngrid = 1, Ngrid_total

        !The folders are named starting from 001 by increments of 1
        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid

        !We will produce a basic parent-level heat map for the grid
        open(filechannel1,file=gridpath0//Ngrid_text//"/"//counter0file)
        do n = 1, counter0_max
        read(filechannel1,FMT="(I8)") counter0(n)
        end do
        close(filechannel1)

        call mapCell(0,counter0,counter0_max,overcrowd0,bounds1,bounds2,multiplier1_0,multiplier2_0,gridpath0//Ngrid_text//"/")
end do

end subroutine analyzeHeatMaps1

end module analyzeHeatMapswithMultipleGrids
