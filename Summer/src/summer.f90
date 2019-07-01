program main
implicit none
    integer:: n, i, numIter
    real(8), allocatable:: h(:), x(:), T(:), TNext(:), lambda(:), A(:), B(:), C(:), D(:)
    real(8):: k, L, error, eps, lambda0, kt, Qv, alp1, alp2, Te1, Te2, ht

    n = 31
    k = 1.3
    L = 50e-2 
    lambda0 = 110
    kt = 6e-3
    Qv = 3e+3
    Te1 = -20
    Te2 = 10
    alp1 = 120
    alp2 = 20

    allocate(h(n), x(n), T(n), TNext(n), lambda(n))
    call makeMesh(n, h, k)
    h = h* L
    
    x(1) = 0
    do i = 2, n
        x(i) = x(i - 1) + h(i)
    end do
    
    error = 1
    eps = 1e-9
    T = 0
    numIter = 0
    ht = 0.001
    allocate(A(n), B(n), C(n), D(n))
    do while(error > eps .and. numIter < 10000)
        lambda = lambda0* (1 + kt* T)
        
        do i = 2, n-1
            A(i) = -4/ (h(i) + h(i+1))* 1/(h(i+1)*(1/lambda(i+1) + 1/lambda(i)))
            B(i) = 4/ (h(i) + h(i+1))*&
                (1/(h(i+1)*(1/lambda(i+1) + 1/lambda(i))) + 1/(h(i)*(1/lambda(i-1) + 1/lambda(i))))&
                + 1/ht
            C(i) = -4/ (h(i) + h(i+1))* 1/(h(i)*(1/lambda(i-1) + 1/lambda(i)))
            D(i) = Qv + T(i)/ht
        end do
        D(1) = - alp1* Te1
        D(n) = - alp2* Te2
        A(1) = lambda(1)/h(2)
        B(1) = -lambda(1)/h(2) -alp1
        B(n) = -lambda(n)/h(n) -alp2
        C(n) = lambda(n)/h(n)
        
        call sweep(A, B, C, D, n, TNext)
        error = maxError(abs(TNext - T), n)
        T = TNext
        numIter = numIter + 1
    end do
    open(1, file='out.txt')
    write(1,*) x(:)
    write(1,*) TNext(:)
    close(1)
    write(*,*) numIter
contains
subroutine makeMesh(n, h, k)
    integer:: n, i
    real(8):: h(n), k, hMin

    hMin = (k - 1)/ (k**((n+1)/2 - 1) - 1) /2
    h(2) = hMin
    h(n) = hMin
    do i = 3, (n + 1)/2
        h(i) = h(i-1)* k
        h(n - i + 2) = h(i)
    end do
end subroutine
subroutine sweep(A, B, C, D, n, TNext)
	real(8):: A(n), B(n), C(n), D(n)
    integer:: n
	real(8):: TNext(n)
	real(8):: alp(n), bet(n)
    integer:: i
    
    alp(2) = - A(1)/ B(1)
    bet(2) = D(1)/ B(1)
    do i = 2, n-1
        alp(i+1) = -A(i)/ (C(i)* alp(i) + B(i))
        bet(i+1) = (D(i) - C(i)* bet(i))/&
            (C(i)* alp(i) + B(i))
    end do
    TNext(n) = (D(n) - C(n)* bet(n))/ (B(n) + C(n)* alp(n))
    do i = n-1 , 1, -1
        TNext(i) = alp(i+1)* TNext(i+1) + bet(i+1)
    end do
end subroutine
real(8) function maxError(error, n)
    integer:: n
    real(8):: error(n)
    maxError = 0
    do i = 1, n
        if(error(i) > maxError) then
            maxError = error(i)
        end if
    end do
end function
end program
