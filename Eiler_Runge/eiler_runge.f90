program pract2
implicit none
real:: x0, a, tmax, h
real, allocatable:: x(:), t(:)
integer:: N,i,m
real:: k

!Открывам файл с исходными данными и создаем файл для выходных данных
!Open the file with the original values and create a file for the output values

open(1, file='vhod.txt', status='old')
open(2, file='out2.txt', status='replace')

!Ввод параметров задачи
!input of task parameters

read(1,*) a, x0, tmax, h

!Находим число шагов и создаем массив нужного размера
!find the number of steps and create an array of the correct size

N = aint(tmax/h + 2.0)
allocate(x(N), t(N))

!Находим значения x и t по формуле явного метода Рунге-Кутты 2-го порядка или Эйлера
!find the values ​​of x and t by the formula of the explicit Runge-Kutta method of the second order or Euler method

x(1) = x0
t(1) = 0
write(*,*)'Use Euler method - press 2, use Runge-Kutta method - press 1'
read(*,*)m
if(m==1)then
    do i = 2,N
    k = x(i-1) + h/2 * f(x(i-1), a)
    x(i) = x(i-1) + h* f(k, a)
    t(i) = (i-1) * h
    write(2,*) i,t(i), x(i)
    end do
endif
if(m==2)then
    do i = 2,N
    x(i) = x(i-1) + h* f(x(i-1), a)
    t(i) = (i-1) * h
    write(2,*) i,t(i),x(i)
    end do
end if

!Описание функции f
!Description of the function f

contains
real function f(x, a)
real, intent(in):: x, a
   f = a * sqrt(x)
end function
end
