program turb
implicit none
        integer:: i, j, n
        real(8), allocatable:: x(:), v(:), h(:), dv(:), mu(:)
        real(8), allocatable:: A(:), B(:), C(:), F(:)
        real(8):: nu, rho, k, p, v0, vn
        
        rho = 1.293
        nu = 13.28
        n = 101
	k = 0.9
	mu(:) = rho* nu
	
	allocate(h(n-1), A(n-2), B(n-2), C(n-2), F(n-2), v(n), &
	x(n), dv(n))
		
	call grid(n-1, k, h)
	stop
	n = n - 1
	x(n/2) = 0.5
	do i = n/2+1, n-1
		x(i) = h(n/2)* ((1 - k** (i-n/2))/ (1 - k)) + 0.5
		x(n-i) = 1 - x(i)
	end do
	n = n + 1
	write(*,*) n
	x(1) = 0.0
	x(n) = 1.0
	stop
	do i = 1, n
	write(*,*) x(i)
	end do
	stop
	
	
	p = 500
	v0 = 0
	vn = 1
	
	call matrix(h, p, mu, v0, vn, n-2, A, B, C, F)
	call sweep(A, B, C, F, n-2, v)
	v(1) = v0
	v(n) = vn

	
        call deriv(v, h, dv, n)
        call calc_mu(nu, rho, x, dv, v, n, mu)

contains
subroutine deriv(v, h, dv, n)
        integer:: n, i
        real(8):: v(n), h(n-1), dv(n)
        
        dv(1) = (v(2) - v(1))/ h(1)
        dv(n) = (v(n) - v(n-1))/ h(n-1)
        do i = 2, n-1
                dv(i) = (v(i+1) - v(i))/ h(i)
        end do      
end subroutine
subroutine calc_mu(nu, rho, x, dv, v, n, mu)
        integer:: n, i
        real(8):: rho, nu, muo, mut, k, demf, l, nut
        real(8):: x(n), dv(n), v(n), mu(n)
        
        muo = rho* nu
        k = 0.41
        demf = 1
        do i = 1, n
                l = k* x(i)
                nut = demf* l** 2* dv(i)
                mut = nut* rho
                mu(i) = muo + mut
        end do
end subroutine
!_______________________________

subroutine grid(n, k, h)
	integer:: n
	real(8):: k
	real(8):: h(n)
	integer:: i

	h(n/2+1) = 0.5* (1 - k)/ (1 - k** (n/2))
	h(n/2) = h(n/2+1)
	do i = n/2+1, n-1
		h(i+1) = h(i)* k
		h(n-i) = h(i+1)
	end do
end subroutine
subroutine matrix(h, p, mu, y0, yn, n, A, B, C, F)
	real(8):: h(n+1), mu(n+2), y0, yn, p
	integer:: n
	real(8):: A(n), B(n), C(n), F(n)
	integer:: i

	do i = 1, n
		C(i) = - 2* (1/ h(i) + 1/ h(i+1))/ (h(i) + h(i+1))
		A(i) = 2/ (h(i)* (h(i) + h(i+1)) )
		B(i) = 2/ (h(i+1)* (h(i) + h(i+1)) )
		F(i) = p/ mu(i+1)
	end do
	F(1) = p/ mu(2) - y0* 2/ (h(1)* (h(1) + h(2)) )
	F(n) = p/ mu(n+1) - yn* 2/ (h(n)* (h(n) + h(n-1)) )
	A(1) = 0
	B(n) = 0
end subroutine
subroutine sweep(A, B, C, F, n, y)
	real(8):: A(n), B(n), C(n), F(n)
	integer:: n
	real(8):: y(n+2)
	real(8):: alp(n), bet(n)
	integer:: i
	
	alp(2) = - B(1)/ C(1)
	bet(2) = F(1)/ C(1)

	do i = 2, n-1
		alp(i+1) = -B(i)/ (A(i)* alp(i) + C(i))
		bet(i+1) = (F(i) - A(i)* bet(i))/&
        	(A(i)* alp(i) + C(i))
	end do

	y(n+1) = (F(n) - A(n)* bet(n))/ (C(n) + A(n)* alp(n))
	do i = n-1 , 1, -1
		y(i+1) = alp(i+1)* y(i+1) + bet(i+1)
	end do
end subroutine
end program
