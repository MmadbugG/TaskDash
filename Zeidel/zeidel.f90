  program zeidel
      implicit none
      real(8), allocatable:: x(:), b(:), A(:,:)
      integer:: i, j, n, iter
      real(8):: s, k, nor, eps, xt(5)=(/1, 2, 3, 4, 5/)

      open(1, file='matrix.txt')
      write(*,*) "Sustem matrix"
      read(1,*) n
      
      allocate(A(n,n), x(n), b(n))

      do i = 1, n
        read(1,*) A(i,:)
        write(*,1) A(i,:)
      end do

      write(*,*)
      write(*,*) "Right side of equation"
      read(1,*) b(:)
      write(*,1) b(:)
      write(*,*)
      write(*,*) "Zero approximation of x"
      read(1,*) x(:)
      write(*,1) x(:)
      write(*,*)

      do i = 1, n
        do j = 1, n
                if (i.ne.j) then 
                  A(i,j) = A(i,j)/A(i,i)
                end if
        end do
        
        b(i) = b(i)/A(i,i)
        
        A(i,i) = 1
      end do
      
      nor = 1
      write(*,*) "Eps"
      eps = 10.0**(-15)
      write(*,"(d10.2)") eps
      write(*,*)
      
      iter = 0
      do while(nor .ge. eps)
        iter = iter + 1
        do i = 1, n
          s = 0
          do j = 1, n
            if (i.ne.j) then 
              s = s + x(j)*A(i,j)
            end if
          end do
          
          x(i) = b(i) - s
        end do
       ! write(*,*) x(:)  
        
        k = 3
        nor = 0
        do i = 1, n
          s = 0
          do j = 1, n
            s = s + A(i, j)*x(j)
          end do
          nor = nor + (abs(s - b(i) ))**k
        end do
        nor = nor**(1.0/k)
        
      end do
      
      
      write(*,*) "Iter: ", iter
      write(*,*)
      write(*,*) "Result"      
      write(*,*) x(:)   
      write(*,*)
      b = x - xt
      write(*,*) sqrt(skal(b, b, n))
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
      
  end program
  
  
  
  
  
  
  
  
