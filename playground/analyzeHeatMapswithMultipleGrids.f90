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






subroutine analyzeHeatMaps2(counter0_optional,counter1_optional,filename_optional)
use PARAMETERS
use mapCellData
use ANALYSIS
implicit none

!FORMATTING OF JPG FILES
character(5) :: variable_length_text
character(Ngrid_text_length) :: Ngrid_text
character(Ngrid_text_length+1) :: folder_text
character(6) :: Ntraj_text

integer,dimension(counter0_max),optional :: counter0_optional
integer,dimension(counter1_max),optional :: counter1_optional
character(*),optional :: filename_optional

integer :: occurence_max
integer :: indexer0, indexer1
integer :: population
real :: var1, var2

!I/O HANDLING
integer :: iostate

!INTEGER INCREMENTALS
integer :: i, j, n



do Ngrid = 1, max(Ngrid_total,1)

        !The folders are named starting from 001 by increments of 1
        write(variable_length_text,FMT="(I5)") Ngrid_text_length
        write(Ngrid_text,FMT="(I0."//trim(adjustl(variable_length_text))//")") Ngrid
	gridpath1 = gridpath0//Ngrid_text//"/"

	!This heat map will look only at counters 0 and 1 to make a heat map

	!If the counters are already open, then no need to read them off again
	if ((present(counter0_optional)).and.(present(counter1_optional)).and.(present(filename_optional))) then
		counter0 = counter0_optional
		counter1 = counter1_optional

	else
	        open(filechannel1,file=gridpath1//counter0file)
	        do n = 1, counter0_max
		        read(filechannel1,FMT="(I8)") counter0(n)
	        end do
	        close(filechannel1)
	
	        open(filechannel1,file=gridpath1//counter1file)
	        do n = 1, counter1_max
		        read(filechannel1,FMT="(I8)") counter1(n)
	        end do
	        close(filechannel1)
	end if

        !Iterate through these cells
        open(filechannel1,file=gridpath0//temporaryfile3)
	occurence_max = overcrowd0

        do i = 1, bounds2
        var2 = i * multiplier2_0
        indexer0 = bounds1*(i-1)
        
        do j = 1, bounds1
        var1 = j * multiplier1_0
        
                population = counter0(indexer0 + j)
		if (population < overcrowd0) then
                	write(filechannel1,FMT="(F5.2,1x,F5.2,1x,I8)") var1, var2, population
			cycle
        	else
			indexer1 = population / key_start
			population = modulo(sum(counter1(indexer1+1:indexer1+1+resolution_0)),key_start)
                	write(filechannel1,FMT="(F5.2,1x,F5.2,1x,I8)") var1, var2, population

			occurence_max = max(occurence_max,population)
		end if
        end do
                write(filechannel1,*) ""
        end do
        close(filechannel1)
        
!
!        do i = 1, bounds1
!        var1 = i * multiplier1_0
!        indexer0 = bounds2*(i-1)
!        
!        do j = 1, bounds2
!        var2 = j * multiplier2_0
!        
!                population = counter0(indexer0 + j)
!		if (population < overcrowd0) then
!                	write(filechannel1,FMT="(F5.2,1x,F5.2,1x,I8)") var1, var2, population
!			cycle
!        	else
!			indexer1 = population / key_start
!			population = sum(counter1(indexer1+1:indexer1+17))
!                	write(filechannel1,FMT="(F5.2,1x,F5.2,1x,I8)") var1, var2, population
!
!			occurence_max = max(occurence_max,population)
!		end if
!        end do
!                write(filechannel1,*) ""
!        end do
!        close(filechannel1)
        
        !Write to the plot
        write(variable_length_text,FMT="(I5)") floor(occurence_max*.9)
        open(filechannel1,file=gridpath0//gnuplotfile)
        write(filechannel1,*) 'set terminal pngcairo size 1600,1300'
	if ((present(counter0_optional)).and.(present(counter1_optional)).and.(present(filename_optional))) then
        	write(filechannel1,*) 'set output "'//gridpath1//filename_optional//'"'
	else
	        write(filechannel1,*) 'set output "'//gridpath1//'HeatMap_TopLevel_Detailed1.png"'
	end if
        write(filechannel1,*) 'unset key'
        write(filechannel1,*) 'set palette defined ( 0 0 1 0, 0.3333 0 0 1, 0.6667 1 0 0,\'
        write(filechannel1,*) '     1 1 0.6471 0 )'
        write(filechannel1,*) 'set cbrange [0:'//trim(adjustl(variable_length_text))//']'
        write(filechannel1,*) 'set cblabel "Population"'
        write(filechannel1,*) 'set title "Configurational Heatmap of an H - H_2 System" font ",32" offset 0,3'
        write(filechannel1,*) 'set xlabel "Var1 (A)" font ",28" offset 0,-2'
        write(filechannel1,*) 'set xtics font ",24"'
        write(filechannel1,*) 'set xrange [2:11]'
        write(filechannel1,*) 'set ylabel "Var2 (A)" font ",28" offset -5,0'
        write(filechannel1,*) 'set ytics font ",24"'
        write(filechannel1,*) 'set yrange [2:12]'
        write(filechannel1,*) 'set cbtics'
        write(filechannel1,*) 'set view map'
        write(filechannel1,*) 'set pm3d interpolate 1,1'
        write(filechannel1,*) 'splot "'//gridpath0//temporaryfile3//'" u 1:2:3 w pm3d'
        close(filechannel1)
        call system(path_to_gnuplot//"gnuplot < "//gridpath0//gnuplotfile)

	!If we only have one counter (that was provided) just exit early
	if ((present(counter0_optional)).and.(present(counter1_optional)).and.(present(filename_optional))) exit
        
end do

end subroutine analyzeHeatMaps2


end module analyzeHeatMapswithMultipleGrids
