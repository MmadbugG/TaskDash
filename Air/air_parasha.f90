module air_curve
implicit none
        type circle !окружность с центром в точке о(x,y) и радиусом r
                real(8):: r
                real(8):: o(2)
        end type
contains
subroutine calc_curve(r1, r2, alp, bet, gam, o1o2, curve, n, k) 
	! функция, строящая половину (верхнюю или нижнюю) профиля
        ! расстояние между точками на профиле одинаковое
        real(8), intent(in):: r1, r2, alp, bet, o1o2
	! входные параметры: 
	! радиусы передней и задней кромок, углы касательных в точках прикреления к кромкам, 
	! угол, задающий положение точки сопряжения, 
	! расстояние между центрами передней и задней кромок
        integer:: i, j, n, k(4), offset !i,j - индексы для циклов, n - количество точек на профиле
        ! k - массив кол-ва точек на каждой окружности профиля, offset - количество уже записанных точек
        real(8):: curve(n, 2), offtet, rr(4), y, gam
        ! curve - массив точек профиля, offtet - пойденная дуга окружности, rr - радиусы окружностей
        real(8):: tet(4), dtet(4), p(4, 2), po(4, 2), l, s(2) ! tet - угол, на оторый опирается дуга каждой окружности
        ! dtet - шаг угла для каждой окружности
        real(8):: m (2,2), fi, rs, Rr1, Rr2, ab, star, a(2), b(2), q(2)
        real(8):: pi = 3.14159265, zero = 0.0
        type(circle):: o1, o2, oR1, oR2! передняя кромка, задняя, первая и вторая сопряженные окружности
        ! радиусы и центры окружностей кромок
	o1%r = r1
	o1%o = (/ 0.0, 0.0 /)
	o2%r = r2
	o2%o = (/ o1o2, zero /)
	! вычисление раидусов сопряженных окружностей
        star = r1* cos(alp) - r2* cos(bet)
	fi = acos(star/ o1o2)
	a = fr(o1, fi + alp)
	b = fr(o2, fi - bet)
	ab = norm(b - a, 2)
        rs = ab/ (2* sin((alp + bet)/ 2))
	! координата y центра вспомогательной окружности
	! От знака зависит направление отсчета гамма
	y = (a(2) * b(2) * o2%o(1))/(a(1) * b(2) - a(2) * b(1) + a(2) * o2%o(1))
	if (y < 0) then
                i = 1
        else
                i = -1
        end if
	gam = gam* i

        Rr1 = rs* sin(gam/ 2)/ sin((alp - bet + gam)/ 2)
        Rr2 = rs* sin((alp + bet - gam)/ 2)/ sin(bet - gam/ 2)
        
        ! радиусы и центры сопряженных окружностей
        oR1%o = a + Rr1* (o1%o - a)/ norm(o1%o - a, 2)
        oR2%o = b + Rr2* (o2%o - b)/ norm(o2%o - b, 2)
        oR1%r = Rr1
        oR2%r = Rr2
        ! условие для нажодения точки сопряжения: сложение или вычетание векторов
        if (Rr1 > Rr2) then
                i = 1
        else
                i = -1
        end if
        ! точка сопряжения
        q = oR1%o + (oR2%o - oR1%o)/ norm(oR2%o - oR1%o, 2)* i* Rr1
        
	! р - массив векторов из центров окружностей, которые будут поворачиваться матричей поворота
	! fr - функция нахождения координат точки на окружности (заданы окружность и угол)
        p(1, :) = fr(o1, pi) - o1%o ! fr(o1, pi) - передняя точка профиля
        p(2, :) = a - oR1%o
        p(3, :) = q - oR2%o
        p(4, :) = b - o2%o
	! точки, вокруг которых производится поворот (центры окружностей)
        po(1, :) = o1%o
        po(2, :) = oR1%o
        po(3, :) = oR2%o
        po(4, :) = o2%o
	! вычисление углов, на которые опираются дуги окружностей
	! vectorangle - функция вычисления угла между векторами
        tet(1) = vectorangle(p(1,:), a, 2)
        tet(2) = vectorangle(p(2,:), q - oR1%o, 2)
        tet(3) = vectorangle(p(3,:), b - oR2%o ,2)
        tet(4) = vectorangle(p(4,:), fr(o2, zero), 2)
	! длина половины профиля (сумма четырех дуг)
	! знак радиуса показывает вогнутость или выпуклость участка профиля
        l = o1%r* tet(1) + abs(oR1%r* tet(2)) + abs(oR2%r* tet(3)) + o2%r* tet(4)
	! нахождение dtet для каждой дуги 
	! Для обеспечения одинакового расстояния между точками профиля 
	! значение берется обратно пропорционально радиусу
	rr = (/ r1, Rr1, Rr2, r2 /)
        dtet = (l/ (n-1)) / rr 
	
        offset = 0
	offtet = 0.0
	! цикл по 4-м окружностям
	! вычисление точек на профиле поворотом векторов
        do i = 1, 4
		j = 1	
		offtet = offtet/ rr(i) ! начальный угол
		m(1,:) = (/ cos(offtet), -sin(offtet) /)
                m(2,:) = (/ sin(offtet), cos(offtet) /)
		s = dot( p(i, :), m, 2) ! поворот начального вектора на угол offtet
		curve(offset + j, :) = s(:) + po(i, :) ! запись первой точки
		j = j + 1 ! увеличение индекса
		! матрица поворота 
                m(1,:) = (/ cos(dtet(i)), -sin(dtet(i)) /)
                m(2,:) = (/ sin(dtet(i)), cos(dtet(i)) /)
		
		! до тех пор пока не превысим угол
                do while(abs(dtet(i)*(j-1) + offtet) < tet(i) - 1e-9) ! шаг* кол-во пройденных точек + начальный угол < tet(i)
			! если равно tet(i), то выходим из цикла (поправка -1е-9)
                        s = dot(s, m, 2) 
                        curve(offset + j, :) = s(:) + po(i, :)
			j = j + 1
                end do
		! вычисление начального сдвига для следующей окружности, 
		! связанного с тем, что на предыдущей не уложилось целое число точек контура
		offtet = abs(dtet(i)*(j-1) + offtet) - tet(i) ! начальный сдвиг по углу
		offtet = offtet* abs(rr(i)) ! начальный сдвиг по дуге
		k(i) = j - 1 ! количество точек, отложенных на текущей окружности
                offset = offset + k(i) ! кол-во записанных точек контура
        end do
	! добавляем последнюю крайнюю точку контура
	curve(n,:) = po(4,:) + (/ o2%r, zero /)
	k(4) = k(4) + 1 ! число точек на задней кромке с учетом крайней

end subroutine
function dot(v, m, n) result(v2)
        real(8):: v(n), v2(n), m(n,n), s
        integer:: i, j, n
        do i = 1, n
                s = 0
                do j = 1, n
                        s = s +  m(j, i)* v(j)
                end do
                v2(i) = s
        end do
end function
function fr(circ, fi) result(v)
	real(8):: fi, v(2)
	type(circle):: circ
	v(1) = circ%r* cos(fi)
	v(2) = circ%r* sin(fi)
	v = v + circ%o
end function	
real(8) function vectorangle(v1, v2, n)
        real(8):: v1(n), v2(n)
	integer:: n
	vectorangle = acos(scal(v1, v2, n)&
	/ (norm(v1, n)* norm(v2, n)))       
end function
real(8) function scal(v1, v2, n)
	real(8):: v1(n), v2(n), s
	integer:: n, i
	s = 0
	do i = 1, n
		s = s + v1(i)* v2(i)
	end do
	scal = s
end function
real(8) function norm(v, n)
	integer:: n, i
	real(8):: s, v(n)
	s = 0
	do i = 1, n
		s = s + v(i)**2
	end do
	norm = sqrt(s)
end function
end module


program air
use air_curve
implicit none
        real(8):: r1, r2, alp1, bet1, gam1, alp2, bet2, gam2, o1o2
        real(8), allocatable:: curve(:, :)
        integer:: n, i, j, k(4) 
        real(8):: pi = 3.14159265       
        r1 = 1.0
        r2 = 0.5
        o1o2 = 10.0
        alp1 = pi/ 8
        bet1 = pi/ 6
        gam1 = 0.07

        alp2 = alp1
        bet2 = -bet1
        gam2 = 0.11
        n = 70
        
	!read(*,*) r1, r2, o1o2 
	!read(*,*) alp1, bet1, gam1, alp2, bet2, gam2
        open(1, file="out_data.txt")
        allocate(curve(n, 2))
        call calc_curve(r1, r2, alp1, bet1, gam1, o1o2, curve, n, k)
        do i = 1, n  - k(4)
                write(1,*) curve(i,:)
        end do

        call calc_curve(r1, r2, alp2, bet2, gam2, o1o2, curve, n, k)
        do i = n - k(4), 1, -1
                write(1,*) curve(i, 1), -curve(i,2)
        end do
end program


