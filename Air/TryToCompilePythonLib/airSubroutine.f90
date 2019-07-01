subroutine calc_curve(r1, r2, alp, bet, gam, o1o2, curve, n)
implicit none
type circle
                real r
                real o(2)
end type

        real, intent(in):: r1, r2, alp, bet, gam, o1o2
        integer:: i, j, n, k(4), offset
        real, intent(out):: curve(n, 2)
        real:: tet(4), dtet(4), p(4, 2), po(4, 2), l, s(2)
        real:: m (2,2), fi, rs, Rr1, Rr2, ab, star, a(2), b(2), q(2)
        real:: pi = 3.14159265
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
