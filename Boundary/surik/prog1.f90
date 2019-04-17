program fourier
implicit none
real(8), allocatable:: T(:,:), X(:), Fo(:)
real(8):: Bi, dFo, eps, h, ht, a1, a2, f1, f2, a3, f3, an, mu, pi
integer:: i, j, n, nt, k
real(8):: T0, Te, alp, l, c, rho, lam, a, Tn

open(1, file='out.txt')
open(2, file='temp.txt')

!Ввод исходных данных
	l = 1
	alp = 110000.0
	eps = 10e-11
	dFo = 0.0005 !шаг по времени
	lam = 0.74
	c = 670.0
	rho = 2500.0
bi = alp*(l/2)/lam !число Био
a = lam/(c*rho) !коэффициент температуропроводности
pi = 4* atan(1.0)



n = 31 !число узлов по координате
nt = 1000 !число узлов по времени
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
	write(2,*) T(:,1)
	ht = dFo
	write(1, *) 0.0
	do i = 2, nt !внешний цикл по времени
	Fo(i) = (i-1)*ht
	do j = 1, n !Внутренний цикл по координате
	X(j) = (j-1)*h
	Tn = 0
		do k = 1, 10
		a1 = pi* (k - 1.0) + eps
		a2 = pi* (k - 0.5)
		mu = bisec(a1, a2, eps, bi)
		an = 2.0* sin(mu)/ (mu + sin(mu)* cos(mu))
		an = an* exp(-mu**2 *Fo(i))* cos(mu* X(j))
		Tn = Tn + an
		if (abs(an) < eps ) then 
			exit
		end if
		end do
	T(j,i) = Tn
	
	write(2,*) T(j,i)
		
	end do
	write(1,*) Fo(i)
	end do

contains 
real(8) function bisec(a1, a2, eps, Bi) !Метод половинного деления
	real(8):: a1, a2, eps, Bi
	real(8):: f1, f2, a3, f3
	f1 = f(a1, Bi)
	f2 = f(a2, Bi)
	if (f1* f2 > 0) then
		write(*,*) 'Error: f1* f2 > 0 ...'
		stop
	end if
	do while (abs(a2 - a1) > eps)
		a3 = (a2 + a1)/ 2
		f3 = f(a3, Bi)
		if (f3 .eq. 0) then
			exit		
		end if
		if (f1* f3 < 0 ) then
			a2 = a3
		else 
			a1 = a3
			f1 = f3
		end if
	end do
	bisec = a3
end function
real(8) function f(x, Bi)
	real(8):: x, Bi
	f = 1.0/tan(x) - x/ Bi
end function

end program




