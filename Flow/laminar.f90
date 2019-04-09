program laminar
implicit none
	integer:: n
	real:: k
	real(8):: f0, y0, yn
	real(8), allocatable:: h(:), A(:), B(:), C(:), F(:), y(:), x(:)
	integer:: i
	
	n = 100	
	if (mod(n,2) == 1) then
		write(*,*) 'n+1'
		n = n + 1
	end if
	k = 0.9
	allocate(h(n), A(n), B(n), C(n), F(n), y(n), x(n))	
	
	call grid(n, k, h)
	f0 = 50
	y0 = 0
	yn = 1
	call matrix(h, f0, y0, yn, n, A, B, C, F)
	!subroutine find n-1 parameter and change n to n-1
	call sweep(A, B, C, F, n, y)

	n = n + 1
	x(n/2) = 0.5
	do i = n/2+1, n-1
		x(i) = h(n/2)* ((1 - k** (i-n/2))/ (1 - k)) + 0.5
		x(n-i) = 1 - x(i)
	end do
		
	!export
	open(1, file='out_data.txt')
	write(1,*) n + 1
	write(1,2) 'i', 'x', 'y'
	write(1,*) ('-', i=1,50)
	write(1,1) 1, 0.0, y0
	do i = 1, n-1
		write(1,1) i+1, x(i), y(i)
	end do
	write(1,1) n+1, 1.0, yn
1 format(i4,2f23.10)
2 format(a4,2a23)
contains
subroutine grid(n, k, h)
	integer:: n
	real:: k
	real(8):: h(n)
	integer:: i

	h(n/2+1) = 0.5* (1 - k)/ (1 - k** (n/2))
	h(n/2) = h(n/2+1)
	do i = n/2+1, n-1
		h(i+1) = h(i)* k
		h(n-i) = h(i+1)
	end do
end subroutine
subroutine matrix(h, f0, y0, yn, n, A, B, C, F)
	real(8):: h(n), f0, y0, yn
	integer:: n
	real(8):: A(n), B(n), C(n), F(n)
	integer:: i

	do i = 1, n-1
		C(i) = - 2* (1/ h(i) + 1/ h(i+1))/ (h(i) + h(i+1))
		A(i) = 2/ (h(i)* (h(i) + h(i+1)) )
		B(i) = 2/ (h(i+1)* (h(i) + h(i+1)) )
		F(i) = f0
	end do
	F(1) = f0 - y0* 2/ (h(1)* (h(1) + h(2)) )
	F(n-1) = f0 - yn* 2/ (h(n)* (h(n) + h(n-1)) )
	A(1) = 0
	B(n-1) = 0
end subroutine
subroutine sweep(A, B, C, F, n, y)
	real(8):: A(n), B(n), C(n), F(n)
	integer:: n
	real(8):: y(n)
	real(8):: alp(n), bet(n)
	integer:: i
	
	n = n - 1
	alp(2) = - B(1)/ C(1)
	bet(2) = F(1)/ C(1)

	do i = 2, n-1
		alp(i+1) = -B(i)/ (A(i)* alp(i) + C(i))
		bet(i+1) = (F(i) - A(i)* bet(i))/&
        	(A(i)* alp(i) + C(i))
	end do

	y(n) = (F(n) - A(n)* bet(n))/ (C(n) + A(n)* alp(n))
	do i = n-1 , 1, -1
		y(i) = alp(i+1)* y(i+1) + bet(i+1)
	end do
end subroutine
end program




