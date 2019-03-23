program razed_per
implicit none
real(8):: Bi, Fo, Fon, eps
real(8), allocatable:: T(:), xh(:), xhd(:)
integer:: i, n
	
	!parameter
	Bi = 10
	Fon = 3

	n = 10
	allocate(T(n), xh(n), xhd(n))
	eps = 1
	Bi = rad_t(eps, Bi, Fo, xhd, n)
	
	
	
	
		
		
			
	
	
	
	

contains
real(8) function rad_t(eps, Bi, Fo, xhd, n)
	integer:: n
	real(8):: eps, Bi, Fo, xhd
	real(8):: pi, Td, a, b
	integer:: i
	pi = 4* atan(1.0)
	
	Td = 0
	i = 1
	do while (.TRUE.)
		write(*,*)

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
	
end function

real(8) function bisec(a, b, eps, Bi)
	real(8):: a, b, eps, Bi
	real(8):: fa, fb
	fa = f(a, Bi)
	fb = f(b, Bi)
	if (fa* fb > 0) then
		write(*,*) 'Error: fa* fb > 0 ...'
		stop
	end if
	do while (abs(b - a) > eps)
		c = (b + a)/ 2
		fc = f(c, Bi)
		if (fc = 0) then
			exit		
		end if
		if (fa* fc < 0 ) then
			b = c
		else 
			a = c
			fa = fc
		end if
	bisec = c
end function

real(8) function f(x, Bi)
	real(8):: x, Bi
	f = 1.0/tan(x) - x/ Bi
end function
end program
