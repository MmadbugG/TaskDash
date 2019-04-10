module air_curve
implicit none
                     
contains
subroutine calc_curve(r1, r2, alp, bet, gam, o1o2, curve, n)
        real:: r1, r2, alp, bet, gam, o1o2
        integer:: i, j, n
        real:: curve(n)
        real:: tet(4), p(4), po(4), m(2,2), ab, fi, rs, Rr1, Rr2
        type circle
                real r
                real o(2)
        end type
        type(circle):: c
        
        rs = ab/ (2* sin((alp + bet)/ 2))
        Rr1 = rs* sin(gam/ 2)/ sin((alp - bet + gam)/ 2)
        Rr2 = rs* sin((alp + bet - gam)/ 2)/ sin(bet - gam/ 2)

end subroutine
real function vectorangle(v1,v2)
        real:: v1(2), v2(2)
        
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

