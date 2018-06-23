module mapCellData
implicit none

contains

subroutine mapCell(indexN,counterN,counterN_max,&
                   scaling1,scaling2)
use f2_parameters
implicit none
logical :: flag1
integer,intent(in) :: indexN,counterN_max
integer,intent(in) :: scaling1, scaling2
integer,intent(in),dimension(counterN_max) :: counterN
integer :: i, j, indexer, population
integer,allocatable :: cells(:,:)
character(20) :: descriptor1
!character(20) :: descriptor2, descriptor2, subcell, line_data

!How many cells are in the heatmap grid?
allocate(cells(scaling1,scaling2))

!Iterate through these cells
do i = 1, scaling1
do j = 1, scaling2

        indexer = indexN + scaling1*(j-1) + i - 1
        population = counterN(indexer)

        cells(i,j) = population

end do
end do



!Then print the data onto temporaryfile3
open(72,file=trim(path4)//trim(temporaryfile3))
write(descriptor1,FMT=FMT4) scaling1
do j = 1, scaling2
write(72,FMT="("//trim(adjustl(descriptor1))//"I9)") cells(:,j)
end do




open(filechannel1,file=path4//temporaryfile1)
write(filechannel1,*) 'set terminal jpeg'
write(filechannel1,*) 'set output "heatmap_2.jpg"'
write(filechannel1,*) 'unset key'
write(filechannel1,*) 'set tic scale 0'
write(filechannel1,*) 'set palette rgbformula -7,2,-7'
write(filechannel1,*) 'set cbrange [0:100]'
write(filechannel1,*) 'set cblabel "Population"'
write(filechannel1,*) 'unset cbtics'
write(filechannel1,*) 'set view map'
write(filechannel1,*) 'splot "'//path4//temporaryfile3//'" matrix w image'
close(filechannel1)
call system("gnuplot < "//path4//temporaryfile1)















!Then call gnuplot to make the image
call system("gnuplot < makeHeatMap")

end subroutine mapCell

end module mapCellData
