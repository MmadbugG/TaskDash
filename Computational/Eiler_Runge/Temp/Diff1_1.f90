program Diff
    
    implicit none
    
    integer :: n, i, mMaxRG, mMinRG, mMaxEu, mMinEu, m, k
    real(16) :: xa, xb, ya, yiEu, yiRG, y1, y2, h, epsEu, epsRG, eps
    real(16), allocatable :: X(:), Yan(:)
    
    xa = 0.0
    xb = 4.0
    ya = 1.0
    
    !eps = 1D-5
    
    n = 20
    h = (xb - xa) / n
    
    allocate(X(n + 1))
    allocate(Yan(n + 1))
    
    do i = 1, n + 1
        X(i) = xa + (i - 1) * h
        Yan(i) = analyticSolution(X(i))
    end do
    
    do k = 1, 6
        
        eps = 10**(-1.0 * k)
    
    yiEu = ya
    yiRG = ya
    epsEu = 0
    epsRG = 0
    do i = 1, n
        m = 2
        mMaxEu = m
        y1 = Yan(i)
        y2 = y1
        do while (.true.)
           
            call Euler(Yan(i), X(i), y1, h, m, dev)
            if (abs(y1 - y2) < eps) exit            
            y2 = y1            
            m = 2 * m
            
        end do
        yiEu = y1
        if (abs(yiEu - Yan(i + 1)) > epsEu) epsEu = abs(yiEu - Yan(i + 1))
        if (i .eq. 1) mMinEu = m
        if (m > mMaxEu) mMaxEu = m
        if (m < mMinEu) mMinEu = m
    end do
    
    do i = 1, n
        m = 2
        mMaxRG = m
        y1 = Yan(i)
        y2 = y1
        do while (.true.)
           
            call RungeKutta3(Yan(i), X(i), y1, h, m, dev)
            if (abs(y1 - y2) < eps) exit            
            y2 = y1            
            m = 2 * m
            
        end do
        yiRG = y1
        if (abs(yiRG - Yan(i + 1)) > epsRG) epsRG = abs(yiRG - Yan(i + 1))
        if (i .eq. 1) mMinRG = m
        if (m > mMaxRG) mMaxRG = m
        if (m < mMinRG) mMinRG = m
    end do
    
    write(*,'(a)') '-------------------------------------------------------------------------------'
    write(*,'(a,f3.1,a,f3.1,a)') 'Interval: [', xa, ';', xb,  ']'
    write(*,'(a,es7.1)') 'Error: ', eps
    write(*,'(a,es7.1)') 'Step: ', h
    write(*,'(a)') '            Euler'
    write(*,'(a,i7/,a,i7/,a,es7.1)') 'Min number of divides: ', mMinEu, 'Max number of divides: ', mMaxEu, 'Error: ', epsEu
    write(*,'(a)') '         Runge-Kutta'
    write(*,'(a,i7/,a,i7/,a,es7.1)') 'Min number of divides: ', mMinRG, 'Max number of divides: ', mMaxRG, 'Error: ', epsRG
    
    end do

    write(*,'(a)') '-------------------------------------------------------------------------------'
    
    deallocate(X)
    deallocate(Yan)
    
    contains
    
    real(16) function analyticSolution(x)
        real(16) :: x
    
        analyticSolution = exp(x) * (log(abs(x)) + 1)
    end function analyticSolution
    
    real(16) function dev(x, y)
        real(16) :: x, y
    
        dev = (4*x + 2*y)/(2*x + 1) 
    end function dev
    
    subroutine Euler(yi, xi, yRes, h, m, dev)
    
        real(16) :: yi, xi, yRes, h
        integer :: m
        
        interface
            real(16) function dev(x, y)
                real(16) :: x, y
            end function
        end interface
        
        real(16) :: hi, x1, y1, y2
        integer :: i
        
        hi = h / m
        
        y1 = yi
        do i = 2, m + 1
            x1 = xi + (i - 1) * hi
            y2 = y1 + hi * dev(x1, y1) 
            y1 = y2
        end do
        
        yRes = y2
        
    end subroutine
    
    subroutine RungeKutta3(yi, xi, yRes, h, m, dev)
    
        real(16) :: yi, xi, yRes, h
        integer :: m
        
        interface
            real(16) function dev(x, y)
                real(16) :: x, y
            end function
        end interface
        
        real(16) :: k1, k2, k3, hi, x1, y1, y2
        integer :: i
        
        hi = h / m
        
        y1 = yi
        do i = 2, m + 1
            x1 = xi + (i - 1) * hi
            k1 = hi * dev(x1, y1)
            k2 = hi * dev(x1 + hi / 2, y1 + k1 / 2)
            k3 = hi * dev(x1 + hi, y1 - k1 + 2 * k2)
            y2 = y1 + 1.0 / 6 * (k1 + 4 * k2 + k3)
            y1 = y2
        end do
        
        yRes = y2
        
    end subroutine
    
end program
