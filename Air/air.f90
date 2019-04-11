module air_curve
implicit none
                     
contains
subroutine calc_curve(r1, r2, alp, bet, gam, o1o2, curve, n)
        real, intent(in):: r1, r2, alp, bet, gam, o1o2
        integer:: i, j, n
        real:: curve(n)
        real:: tet(4), p(4), po(4), m (2,2), fi, rs, Rr1, Rr2, &
	ab, star, a(2), b(2)
        type circle
                real r
                real o(2)
        end type
        type(circle):: o1, o2
	o1%r = 4
	o1%o = 0
	a = 2
	b = func(o1, a(1))
	write(*,*) b
        star = r1* cos(alp) - r2* cos(bet)
	fi = acos(star/ o1o2)
        rs = ab/ (2* sin((alp + bet)/ 2))
        Rr1 = rs* sin(gam/ 2)/ sin((alp - bet + gam)/ 2)
        Rr2 = rs* sin((alp + bet - gam)/ 2)/ sin(bet - gam/ 2)

end subroutine
function func(circ, fi) result(v)
	real:: fi
	type circle
                real r
                real o(2)
        end type
	type(circle):: circ
	real:: v(2)
	v(1) = fi
	v(2) = circ%r

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
end module

program air
use air_curve
implicit none
        real:: r1, r2, alp1, bet1, gam1, alp2, bet2, gam2, o1o2
        real, allocatable:: curve(:)
        integer:: n, i, j        
        
        n = 100
        allocate(curve(n))
        call calc_curve(r1, r2, alp1, bet1, gam1, o1o2, curve, n)
end program

