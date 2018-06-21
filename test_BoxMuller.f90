program test_BoxMuller
implicit none
real :: Z1, Z2
real,parameter :: pi2 = 2.0*3.14159
integer :: seed,sample_size,i
real,allocatable :: U1(:), U2(:), X(:), Y(:), XY(:)


seed = 699
seed = rand(seed)

sample_size = 100000
allocate(U1(sample_size),U2(sample_size),X(sample_size),Y(sample_size))
allocate(XY(sample_size*2))

do i = 1, sample_size
     U1(i) = rand()
     U2(i) = rand()

     Z1 = sqrt(-2.0*log(U1(i)))
     Z2 = pi2*U2(i)

     X(i) = Z1*cos(Z2)
     Y(i) = Z1*sin(Z2)
end do

XY(1:sample_size) = X
XY(sample_size+1:sample_size*2) = Y



open(69,file="temporaryfile1.dat")
do i = 1, sample_size
	write(69,*) U1(i), U2(i), X(i), Y(i)
end do
close(69)

open(69,file="temporaryfile2.dat")
do i = 1, sample_size*2
	write(69,*) XY(i)
end do
close(69)

open(69,file="temporaryfile3.txt")
write(69,*) 'set terminal jpeg'
write(69,*) 'set output "test_BoxMuller.jpg"'
write(69,*) 'set multiplot layout 2,3 rowsfirst'
write(69,*) 'unset key'
write(69,*) 'set border 3'
write(69,*) 'set boxwidth 0.05 absolute'
write(69,*) 'bin_width = 0.1'
write(69,*) 'bin_number(x) = floor(x/bin_width)'
write(69,*) 'rounded(x) = bin_width * (bin_number(x) + 0.5)'
write(69,*) 'set style fill solid 1.0 noborder'
write(69,*) 'plot "temporaryfile1.dat" u (rounded($1)):(1) '//&
            'smooth frequency with boxes'
write(69,*) 'plot "temporaryfile1.dat" u (rounded($2)):(2) '//&
            'smooth frequency with boxes'
write(69,*) 'plot "temporaryfile1.dat" u (rounded($3)):(3) '//&
            'smooth frequency with boxes'
write(69,*) 'plot "temporaryfile1.dat" u (rounded($4)):(4) '//&
            'smooth frequency with boxes'
write(69,*) 'plot "temporaryfile2.dat" u (rounded($1)):(1) '//&
            'smooth frequency with boxes'
close(69)

call system("gnuplot < temporaryfile3.txt")
call system("rm temporaryfile1.dat temporaryfile2.dat")
call system("rm temporaryfile3.txt")

end program test_BoxMuller
