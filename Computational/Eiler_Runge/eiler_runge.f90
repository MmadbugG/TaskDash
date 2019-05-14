program pract2
implicit none
real, allocatable:: x(:), y(:), y_old(:)
integer:: N, i, m, iter, p
real:: k1, k2, k3, eps, a, b, y0, h

!Открывам файл с исходными данными и создаем файл для выходных данных
!Open the file with the original values and create a file for the output values

open(1, file='vhod.txt', status='old')
open(2, file='out2.txt', status='replace')

!Ввод параметров задачи
!input of task parameters

read(1,*) a, b, y0, eps

!Находим число шагов и создаем массив нужного размера
!find the number of steps and create an array of the correct size
h = 0.01

!Находим значения x и t по формуле явного метода Рунге-Кутты 2-го порядка или Эйлера
!find the values of x and t by the formula of the explicit Runge-Kutta method of the second order or Euler method

iter = 1
allocate(y_old(1), y(1), x(1))
y_old(1) = 0
write(*,*)'Use Euler method - press 2, use Runge-Kutta method - press 1'
read(*,*)m
if(m.ne.1 .and. m.ne. 2) then
        write(*,*) "Choose 1 or 2!!!"
        stop
end if

if(m==1) then 
       p = 3
      else
       p = 1
end if 
do while(1>0)
n = 2**iter
deallocate(x, y)
allocate(x(n), y(n))

x(1) = a
y(1) = y0
h = abs(b-a)/(n-1)
if(m==1)then
    do i = 1,N-1
          k1 = f(x(i), y(i))
          k2 = f(x(i) + h/2.0, y(i) + h/2.0 *k1)
          k3 = f(x(i) + h, y(i) - h* k1 + 2.0* h* k2)
          y(i+1) = y(i) + h/6.0* (k1 + 4*k2 + k3) 
          x(i+1) = i* h + a
    end do
end if
if(m==2)then
    do i = 2,N
    y(i) = y(i-1) + h* f(x(i-1), y(i-1))
    x(i) = (i-1) * h + a
    end do
end if
iter = iter + 1
if(infty_norm(y(::2) - y_old(:), n/2)/(2**p - 1) < eps) then
        exit
end if
deallocate(y_old)
allocate(y_old(n))
y_old = y
end do
!Описание функции f
!Description of the function f
write(*,*) iter
do i = 1,n
        write(2,*) i, x(i), y(i)
end do
contains
real function infty_norm(x, n)
        real x(n), res
        integer n
        res = 0
        do i = 1, n
                if( abs(x(i)) > res) then
                        res = abs(x(i))
                end if
        end do
        infty_norm = res
end function
real function f(x, y)
real, intent(in):: x, y
   f = (4*x + 2*y)/(2*x + 1)
end function
end program
