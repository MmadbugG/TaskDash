program separation_var
implicit none
	real(8):: Bi, Fon, eps, h, ht
	real(8), allocatable:: Td(:, :), xh(:), Fo(:)
	integer:: i, j, n, nt
	real(8):: xsl !(/0.0, 0.5, 1.0/)
	
	Bi = 1
	Fon = 1
	
	n = 41 !nodes length, the length of interval is equal to n - 1
	nt = 35
	
	allocate(xh(n), Fo(nt), Td(n,nt))
	ht = Fon/ nt
	h = 1.0/ (n - 1)
	
	eps = 1.0e-2

        do i = 1, n
                xh(i) = (i-1)* h
        end do
	do i = 1, nt
	        Fo(i) = i* ht
	        do j = 1, n
	                Td(j,i) = series_T(eps, Bi, Fo(i), xh(j))
	        end do
	end do
	call export_txt(Td, Bi, Fo, xh, Fon, n, nt)
	xsl = 0.5
	call slice(xh, xsl, n)

contains
real(8) function series_T(eps, Bi, Fo, xh)
	integer:: n
	real(8):: eps, Bi, Fo, xh
	real(8):: pi, Td, a, b, mu
	integer:: i
	pi = 4* atan(1.0)
	
	Td = 0
	do i = 1, 100
		a = pi* (i - 1.0) + eps
		b = pi* (i - 0.5)
		mu = bisec(a, b, eps, Bi)
		a = 2.0* sin(mu)/ (mu + sin(mu)* cos(mu))
		a = a* exp(-mu**2 *Fo)* cos(mu* xh)
		Td = Td + a
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
subroutine export_txt(Td, Bi, Fo, xh, Fon, n, nt)
	real(8):: Td(n, nt), Bi, Fo(nt), Fon, xh(n)
	integer:: n, nt
	real(8):: T(nt, n), time(nt)	
	integer:: i
	character(len=15):: str
	
	open(1, file='file_out.txt')
	write(1,'(a40/)') 'Created from separation_var program'
	write(1,'(a20/)') 'Parameter: '
	write(1,'(a6,f12.1)') 'Bi = ', Bi
	write(1,'(a6,f12.1/)') 'Fo = ', Fon
	write(1,'(a20,i4)') 'Coordinate step: ', n
	write(1,'(a20,i4/)') 'Time step: ', nt
	write(1,*) ('-', i=1,30 )
	write(str , *) n
	write(1,'(16x,'//str//'f11.3/)') xh(:)
	do i = 1, nt
		write(1,'(i5,f11.4,a3,'//str//'f11.7)') &
		i, Fo(i), ' | ', Td(:, i)
	end do
	close(1)	
end subroutine
subroutine slice(xh, xsl, n)
        real(8):: xh(n), xsl
        integer:: isl, n
        real(8):: dist, buff
        integer:: i, j
        buff = 1.0
        do i = 1, n
                dist = abs(xh(i) - xsl)
                if (dist > buff) then
                        exit
                end if
                buff = dist
        end do
        isl = i - 1
        write(*,*) xh(isl)        
end subroutine
end program




