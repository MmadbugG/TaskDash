program turb
implicit none
        integer:: i, j, n
        real(8), allocatable:: x(:), v(:), h(:), dv(:), mu(:)
        real(8), allocatable:: A(:), B(:), C(:), F(:), vplus(:), xplus(:)
        real(8):: nu, rho, p, v0, vn, q, h1
        
        rho = 1.293
        nu = 13.28
	q = 1.1 ! q > 1
	h1 = 0.01

	n = (int(log((h1* 2 - 1 + q)/ (2.0* h1))/ log(q)) + 1)* 2 + 1
	
	allocate(h(n-1), A(n-2), B(n-2), C(n-2), F(n-2), v(n), &
	x(n), dv(n), mu(n), vplus(n), xplus(n))
	
	call grid(n, q, h1, h, x)
        
	p = 45
	v0 = 0
	vn = 1
	mu(:) = rho* nu
	v(1) = v0
	v(n) = vn
	
	do i = 1, 100000
        	call matrix(h, p, mu, v0, vn, n-2, A, B, C, F)
        	call sweep(A, B, C, F, n-2, v)
                call deriv(v, h, dv, n)
                call calc_mu(nu, rho, x, dv, v, n, mu, vplus, xplus)
        end do
        open(1, file="turb_out_plus.txt")
        do i = 1, n
        write(1,*) xplus(i), vplus(i)
        end do
        close(1)
        open(1, file="turb_out.txt")
        do i = 1, n
        write(1,*) x(i), v(i)
        end do
        close(1)

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
subroutine calc_mu(nu, rho, x, dv, v, n, mu, vplus, xplus)
        integer:: n, i
        real(8):: rho, nu, muo, mut, k, dempf, l, nut, &
        vt, tau, vplus(n), xplus(n), aa        
        real(8):: x(n), dv(n), v(n), mu(n)
        
        
        muo = rho* nu
        k = 0.41
        aa = 26
        do i = 1, n
                tau = nu* abs(dv(i))
                vt = sqrt(tau/ rho)
                vplus(i) = v(i)/ vt
                xplus(i) = x(i)* vt/ nu
                dempf = (1 - exp(xplus(i)/aa))**2
                l = k* x(i)
                nut = dempf* l** 2* abs(dv(i))
                mut = nut* rho
                mu(i) = muo + mut
        end do
end subroutine
!_______________________________

subroutine grid(n, q, h1, h, x)
	integer:: n
	real(8):: q, h1
	real(8):: h(n-1), x(n)
	integer:: i
	
	i = 1
	x(1) = 0
	x(2) = h1
	h(1) = h1
	do while(x(i+1) <= 0.5)
	        i = i + 1 
	        h(i) = h(i-1)* q
	        x(i+1) = x(i) + h(i)
	end do
	q = 0.5/x(i+1)
	do i = 2, n/2
	        h(i) = h(i)* q
	        x(i+1) = x(i) + h(i)
	end do
	do i = 1, n/2
	        h(n-i) = h(i)
	        x(n-i+1) = 1 - x(i)
	end do
	x(n/2 + 1) = 0.5
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
	F(n) = p/ mu(n+1) - yn* 2/ (h(n+1)* (h(n+1) + h(n)) )
	A(1) = 0
	B(n) = 0
end subroutine
subroutine sweep(A, B, C, F, n, v)
	real(8):: A(n), B(n), C(n), F(n)
	integer:: n
	real(8):: v(n+2)
	real(8):: alp(n), bet(n)
	integer:: i

	alp(2) = - B(1)/ C(1)
	bet(2) = F(1)/ C(1)

	do i = 2, n-1
		alp(i+1) = -B(i)/ (A(i)* alp(i) + C(i))
		bet(i+1) = (F(i) - A(i)* bet(i))/&
        	(A(i)* alp(i) + C(i))
	end do

	v(n+1) = (F(n) - A(n)* bet(n))/ (C(n) + A(n)* alp(n))
	do i = n-1 , 1, -1
		v(i+1) = alp(i+1)* v(i+2) + bet(i+1)
	end do
end subroutine
end program
