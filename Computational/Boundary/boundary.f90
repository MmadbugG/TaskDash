program boundary
implicit none
    real(8):: xa, xb, ya, yb, h
    integer:: i, n
    real(8), allocatable:: xh(:), yh(:), A(:), B(:), C(:), F(:), y_exact(:)

    xa = 0.2
    xb = 1.0
    ya = exact_solution(xa)
    yb = exact_solution(xb)
    open(1, file='out.txt')
    n = 3
    do while(n < 10**6) 
        h = (xb - xa)/ (n - 1)
        allocate(xh(n), yh(n), A(n-2), B(n-2), C(n-2), F(n-2), y_exact(n))
        do i = 1, n
            xh(i) = h*(i-1) + xa
            y_exact(i) = exact_solution(xh(i))
        end do
        
        do i = 1, n-2
            f(i) = 3* exp(xh(i+1))
            a(i) = 1/ h **2 - (1 + sin(xh(i+1)) **2)/ h
            b(i) = -2/ h **2 + (1 + sin(xh(i+1)) **2)/ h + cos(xh(i+1))**2
            c(i) = 1/ h **2
        end do
        f(1) = f(1) + ya* ((1 + sin(xh(2)) **2)/ h - 1/ h **2)
        f(n-2) = f(n-2) - yb/ h**2
        a(1) = 0
        c(n-2) = 0

        call sweep(a, c, b, f, n-2, yh)
        yh(1) = ya
        yh(n) = yb
        
        write(1,*) n, max(abs(yh - y_exact), n)
        deallocate(xh, yh, A, B, C, F, y_exact)
        n = n* 2
    end do
    close(1)
contains
real(8) function max(x, n)
    real(8) x(n)
    integer n
    max = x(1)
    do i =2, n
        if( x(i) > max) max = x(i)
    end do
end function
real(8) function exact_solution(x)
    real(8) x
    exact_solution =  exp(x)
end function
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
        y(i+1) = alp(i+1)* y(i+2) + bet(i+1)
    end do
end subroutine
end program
