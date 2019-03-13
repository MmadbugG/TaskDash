program skalar
implicit none
real(8), allocatable:: x1(:), x2(:), A(:,:)
integer i, n, k, iter
real(8) lambd, lam0, eps

open(1, file='in.txt')
read(1,*) n
allocate(x1(n),x2(n), A(n,n))
write(*,*) 'Matrix'
do i = 1, n
        read(1,*) A(i,:)
        write(*,1) A(i,:)
end do
read(1,*) x1(:)
write(*,*)
write(*,*) 'Zero approximation of x1'
write(*,1) x1(:)
write(*,*)
lam0 = 1
eps = 1.0e-10
do while (abs(lambd - lam0) > eps)
        if (iter > 1000) then
                write(*,*) 'Ooops..., number of iteration is over 1000!'
                exit
        end if
        lam0 = lambd
        iter = iter + 1
        x1 = x1/sqrt(skal(x1,x1,n))
        do k = 1, n
                x2(k) = skal(A(k,:),x1,n)
        end do
        lambd = skal(x1,x2,n)
        x1 = x2
end do

write(*,*) "Iter: ", iter
write(*,*)
write(*,*) "Result lambda"      
write(*,*) lambd  
write(*,*) 
write(*,*) "Result x1"
write(*,*) x1(:)
write(*,*)
write(*,*) "||A x1 - lambda x1||"
write(*,*) norm(A, x1, n, lambd*x1)
1 format(10f15.5)
contains
real(8) function skal(x1, x2, n)
        integer i, n
        real(8) x1(n), x2(n)
        skal = 0
        do i = 1, n
                skal = skal + x1(i)*x2(i) 
        end do
end
real(8) function norm(A, x, n, b)
        integer n
        real(8) A(n,n), x(n), b(n)
        integer k, i, j
        real(8) s
        
        k = 2
        norm = 0
        do i = 1, n
          s = 0
          do j = 1, n
            s = s + A(i, j)*x(j)
          end do
          norm = norm + (abs(s - b(i) ))**k
        end do
        norm = norm**(1.0/k)
end
end program





