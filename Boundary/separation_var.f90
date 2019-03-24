program separation_var
implicit none
	real(8):: Bi, Fo, Fon, eps, h
	real(8), allocatable:: Td(:, :), xhd(:)
	integer:: i, n
	
	!parameter
	Bi = 10
	Fon = 3
	n = 30
	h = 1.0/ float(n - 1)
	

	n = 10
	allocate(Td(n, n), xhd(n))
	eps = 1
	write(*,*) ('-',i=1,79)
	Bi = series_T(eps, Bi, Fo, xhd(1), n)
	
	
	
	
		
		
			
	
	
	
	

contains
real(8) function series_T(eps, Bi, Fo, xhd, n)
	integer:: n
	real(8):: eps, Bi, Fo, xhd
	real(8):: pi, Td, a, b, mu
	integer:: i
	pi = 4* atan(1.0)
	
	Td = 0
	i = 1
	do while (.TRUE.)
		a = pi* (n - 1.0) + eps
		b = pi* (n - 0.5)
		mu = bisec(a, b, eps, Bi)
		a = 2.0* exp(-mu**2 *Fo)* cos(mu* xhd)
		Td = Td + a

		i = i + 1
		if (abs(a/ Td) < eps ) then 
			exit
		end if
	end do
	series_T = Td	
end function

real(8) function bisec(a, b, eps, Bi)
	real(8):: a, b, eps, Bi
	real(8):: fa, fb, c, fc
	fa = f(a, Bi)
	fb = f(b, Bi)
	if (fa* fb > 0) then
		write(*,*) 'Error: fa* fb > 0 ...'
		stop
	end if
	do while (abs(b - a) > eps)
		c = (b + a)/ 2
		fc = f(c, Bi)
		if (fc .eq. 0) then
			exit		
		end if
		if (fa* fc < 0 ) then
			b = c
		else 
			a = c
			fa = fc
		end if
	end do
	bisec = c
end function

real(8) function f(x, Bi)
	real(8):: x, Bi
	f = 1.0/tan(x) - x/ Bi
end function
subroutine export(Td, Bi, Fo, Fon, n, nt)
	real(8):: Td(nt, n), Bi, Fo(nt), Fon
	integer:: n, nt
	real(8):: T(nt, n), time(nt)	
	integer:: i
	


	open(1, file='file_out.txt')
	write(1,*) 'Created from separation_var program'
	write(1,*) 'Parameter: \n'
	write(1,*) 'Bi = ', Bi
	write(1,*) 'Fo = ', Fon
	write(1,*) 'Coordinate step: ', n
	write(1,*) 'Time step: '
	write(1,*) ('-',i=1,79)
	do i = 1, nt
		write(1,*) time(i)
		write(1,*) ' | '
		write(1,*) T(i, :)
	end do
	close(1)	


end subroutine
end program










