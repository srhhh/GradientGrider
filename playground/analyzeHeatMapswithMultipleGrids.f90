module analyzeHeatMapswithMultipleGrids
implicit none



contains

subroutine analyzeHeatMaps1()
use PARAMETERS
use mapCellData
use ANALYSIS
implicit none

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(6) :: Ntraj_text

!I/O HANDLING
integer :: iostate

!INTEGER INCREMENTALS
integer :: n


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
