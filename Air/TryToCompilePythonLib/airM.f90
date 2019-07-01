        type circle
                real r
                real o(2)
        end type


subroutine airM(r1, r2, alp, bet, gam, o1o2, curve, n)
implicit none
        real, intent(in):: r1, r2, alp, bet, gam, o1o2
        integer, intent(in):: n
	integer:: i, j, n, k(4), offset
        real, intent(out):: curve(n, 2)
        real:: tet(4), dtet(4), p(4, 2), po(4, 2), l, s(2)
        real:: m (2,2), fi, rs, Rr1, Rr2, ab, star, a(2), b(2), q(2)
        real:: pi = 3.14159265
	!f2py    intent(in)    r1
	!f2py    intent(in)    r2
	!f2py    intent(in)    alp
	!f2py    intent(in)    bet
	!f2py    intent(in)    gam
	!f2py    intent(in)    o1o2
	!f2py    intent(in)    n
	!f2py    intent(out)   curve 
        type(circle):: o1, o2, oR1, oR2
	o1%r = r1
	o1%o = (/ 0.0, 0.0 /)
	o2%r = r2
	o2%o = (/ o1o2, 0.0 /)
	
        star = r1* cos(alp) - r2* cos(bet)
	fi = acos(star/ o1o2)
	a = fr(o1, fi + alp)
	b = fr(o2, fi - bet)
	ab = norm(b - a, 2)
        rs = ab/ (2* sin((alp + bet)/ 2))
        Rr1 = rs* sin(gam/ 2)/ sin((alp - bet + gam)/ 2)
        Rr2 = rs* sin((alp + bet - gam)/ 2)/ sin(bet - gam/ 2)

        oR1%o = a + Rr1* (o1%o - a)/ norm(o1%o - a, 2)
        oR2%o = b + Rr2* (o2%o - b)/ norm(o2%o - b, 2)
        oR1%r = Rr1
        oR2%r = Rr2
        if (Rr1 > Rr2) then
                i = 1
        else
                i = -1
        end if
        q = oR1%o + (oR2%o - oR1%o)/ norm(oR2%o - oR1%o, 2)* i* Rr1
        
        p(1, :) = fr(o1, pi) - o1%o
        p(2, :) = a - oR1%o
        p(3, :) = q - oR2%o
        p(4, :) = b - o2%o
        po(1, :) = o1%o
        po(2, :) = oR1%o
        po(3, :) = oR2%o
        po(4, :) = o2%o
        tet(1) = vectorangle(p(1,:), a, 2)
        tet(2) = vectorangle(p(2,:), q - oR1%o, 2)
        tet(3) = vectorangle(p(3,:), b - oR2%o ,2)
        tet(4) = vectorangle(p(4,:), fr(o2, 0.0), 2)
        l = o1%r* tet(1) + abs(oR1%r* tet(2)) + abs(oR2%r* tet(3)) + o2%r* tet(4)
        dtet = (l/ (n-1))* 1/ (/ r1, Rr1, Rr2, r2 /)
        k = nint(tet/dtet)
        dtet = tet/ k
        k = abs(k)
        curve(1, :) = p(1, :) + po(1, :)
        offset = 1
        do i = 1, 4
                m(1,:) = (/ cos(dtet(i)), -sin(dtet(i)) /)
                m(2,:) = (/ sin(dtet(i)), cos(dtet(i)) /)
                s = p(i, :)
                do j = 1, k(i)
                        s = dot(s, m, 2)
                        curve(offset + j, :) = s(:) + po(i, :)
                end do
                offset = offset + k(i)
        end do
end subroutine
function dot(v, m, n) result(v2)
        real:: v(n), v2(n), m(n,n), s
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
	real:: fi, v(2)
	type(circle):: circ
	v(1) = circ%r* cos(fi)
	v(2) = circ%r* sin(fi)
	v = v + circ%o
end function	
real function vectorangle(v1, v2, n)
        real:: v1(n), v2(n)
	integer:: n
	vectorangle = acos(scal(v1, v2, n)&
	/ (norm(v1, n)* norm(v2, n)))       
end function
real function scal(v1, v2, n)
	real:: v1(n), v2(n), s
	integer:: n, i
	s = 0
	do i = 1, n
		s = s + v1(i)* v2(i)
	end do
	scal = s
end function
real function norm(v, n)
	integer:: n, i
	real:: s, v(n)
	s = 0
	do i = 1, n
		s = s + v(i)**2
	end do
	norm = sqrt(s)
end function

