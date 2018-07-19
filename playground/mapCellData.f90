module mapCellData
implicit none

contains

subroutine mapCell(indexN,counterN,counterN_max,overcrowd,&
                   scaling1,scaling2,multiplier1,multiplier2,path_to_grid)
use PARAMETERS
implicit none
logical :: flag1
integer,intent(in) :: indexN,counterN_max,overcrowd
integer,intent(in) :: scaling1, scaling2
real,intent(in) :: multiplier1, multiplier2
integer,intent(in),dimension(counterN_max) :: counterN
integer :: i, j, indexer, population
real :: var1, var2
character(*),intent(in) :: path_to_grid
character(20) :: descriptor1

!Iterate through these cells
open(filechannel1,file=path_to_grid//temporaryfile3)
do i = 1, scaling1
var1 = i * multiplier1

do j = 1, scaling2
var2 = j * multiplier2

        indexer = indexN + scaling2*(i-1) + j
        population = counterN(indexer)

        write(filechannel1,FMT="(F5.2,1x,F5.2,1x,I8)") var1, var2, population

end do
        write(filechannel1,*) ""
end do
close(filechannel1)

!Write to the plot
write(descriptor1,FMT="(I8)") overcrowd
open(filechannel1,file=path_to_grid//gnuplotfile)
write(filechannel1,*) 'set terminal jpeg size 1600,1300'
write(filechannel1,*) 'set output "'//path_to_grid//'HeatMap_TopLevel.jpg"'
write(filechannel1,*) 'unset key'
write(filechannel1,*) 'set palette rgbformula -7,2,-7'
write(filechannel1,*) 'set cbrange [0:'//trim(adjustl(descriptor1))//']'
write(filechannel1,*) 'set cblabel "Population"'
write(filechannel1,*) 'set xlabel "Var1 (A)"'
write(filechannel1,*) 'set ylabel "Var2 (A)"'
write(filechannel1,*) 'set cbtics'
write(filechannel1,*) 'set view map'
write(filechannel1,*) 'set pm3d interpolate 1,1'
write(filechannel1,*) 'splot "'//path_to_grid//temporaryfile3//'" u 1:2:3 w pm3d'
close(filechannel1)
call system("gnuplot < "//path_to_grid//gnuplotfile)

!Then call gnuplot to make the image
call system("gnuplot < makeHeatMap")

end subroutine mapCell

end module mapCellData
