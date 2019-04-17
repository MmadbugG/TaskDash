program fourier
implicit none
real(8), allocatable:: T(:,:), X(:), Fo(:)
real(8):: Bi, dFo, h, ht, a1, a2, f1, f2, a3, f3, an, mu, pi
integer:: i, j, n, nt, k
real(8):: T0, Te, alp, l, c, rho, lam, a, Tn

open(1, file='out.txt')
open(2, file='temp.txt')

!Ввод исходных данных
	l = 1
	alp = 110000.0
	dFo = 0.2 !шаг по времени
	lam = 0.74
	c = 670.0
	rho = 2500.0
bi = alp*(l/2)/lam !число Био
a = lam/(c*rho) !коэффициент температуропроводности
pi = 4* atan(1.0)



n = 41 !число узлов по координате
nt = 3200 !число узлов по времени
	allocate(T(n,nt), X(n), Fo(nt))
	h = 1./(n-1) !шаг по координате
	!write(1,*) h

	do j = 1, n !координатная сетка
	X(j) = (j - 1)* h
	write(1,*) X(j)
	end do

	!do i = 2, nt !временная сетка
	Fo(1) = 0
	!Fo(i) = (i - 1)* dFo
	!write(1,*) Fo(i)
	!end do

	T(:,1) = 1.0 !НУ
	ht =  h**2 /4 
	write(1, *) 0.0
		
	do i = 2, nt!внешний цикл по времени
	        Fo(i) = (i - 1)* ht
		write(1,*) Fo(i)
	        do j = 2, n-1 !Внутренний цикл по координате
	                T(j, i) = (1.0 - 2.0* ht/ (h** 2))* T(j, i-1) + &
	                (ht/ (h** 2))* (T(j+1, i-1) + T(j-1, i-1))
	        end do
	        T(1, i) = T(2, i)
	        T(n, i) = T(n-1, i)/ (1.0 + h* bi)
	end do
	
	do i = 1, nt 
	do j = 1, n 
	write(2,*) T(j,i)
	end do
	end do

end program




