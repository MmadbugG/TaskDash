program pract2
implicit none
real, allocatable:: x(:), y(:)
integer:: N, i, m
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
N = aint((b - a)/ h + 2.0)
allocate(x(N), y(N))

!Находим значения x и t по формуле явного метода Рунге-Кутты 2-го порядка или Эйлера
!find the values of x and t by the formula of the explicit Runge-Kutta method of the second order or Euler method

x(1) = a
y(1) = y0

write(*,*)'Use Euler method - press 2, use Runge-Kutta method - press 1'
read(*,*)m
if(m==1)then
    do i = 1,N-1
          k1 = f(x(i), y(i))
          k2 = f(x(i) + h/2.0, y(i) + h/2.0 *k1)
          k3 = f(x(i) + h, y(i) - h* k1 + 2.0* h* k2)
          y(i+1) = y(i) + h/6.0* (k1 + 4*k2 + k3) 
          x(i+1) = i* h
          write(2,*) i, y(i+1), x(i+1)
    end do
endif
if(m==2)then
    do i = 2,N
    y(i) = y(i-1) + h* f(x(i-1), y(i-1))
    x(i) = (i-1) * h
    write(2,*) i, y(i), x(i)
    end do
end if

!Описание функции f
!Description of the function f

contains
real function f(x, y)
real, intent(in):: x, y
   f = (4*x + 2*y)/(2*x + 1)
end function
end
