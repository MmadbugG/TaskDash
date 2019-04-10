program air
implicit none
        real:: r1, r2, alp1, bet1, gam1, alp2, bet2, gam2, o1o2
        integer:: i, j, n
        real, allocatable:: curve(:)
        
        n = 100
        allocate(curve(n))
        call calc_curve(r1, r2, alp1, bet1, gam1, o1o2, curve, n)

contains
subroutine calc_curve(r1, r2, alp, bet, gam, o1o2, curve, n)
        real:: r1, r2, alp, bet, gam, o1o2
        integer:: i, j, n
        real:: curve(n)
        real:: tet(4), p(4), po(4), m(2,2), ab, fi, rs, Rr2, Rr2
        type circle
                real r
                real ox
                real oy
        end type
        type(circle):: ci1
        
        rs = ab/ (2* sin((alp + bet)/ 2)
        Rr1 = rs* sin(gam/ 2)/ sin((alp - bet + gam)/ 2)
        Rr2 = rs* sin((alp + bet - gam)/ 2)/ sin(bet - gam/ 2)
        

end subroutine
end program
