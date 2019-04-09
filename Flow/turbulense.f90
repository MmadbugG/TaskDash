program turb
implicit none
        integer:: i, j, n
        real(8), allocatable:: x(:), v(:), h(:), dv(:), mu(:)
        character*1 c
        
        open(1, file="out_data.txt", status="old")
        read(1,*) n
        allocate(x(n), v(n), h(n-1), dv(n), mu(n))
        
        read(1,*)
        read(1,*)
        do j = 1, n
                read(1,*) i, x(j), v(j)
        end do
        close(1)
        do j = 2, n
                h(j-1) = x(j) - x(j-1)
        end do
        call deriv(v, h, dv, n)
        call calc_mu(x, dv, v, n, mu)

contains
subroutine deriv(v, h, dv, n)
        integer:: n, i
        real(8):: v(n), h(n-1), dv(n)
        
        dv(1) = (v(2) - v(1))/ h(1)
        dv(n) = (v(n) - v(n-1))/ h(n-1)
        do i = 2, n-1
                dv(i) = (v(i+1) - v(i))/ h(i)
        end do      
end subroutine
subroutine calc_mu(x, dv, v, n, mu)
        integer:: n, i
        real(8):: rho, nu, muo, mut, k, demf, l, nut
        real(8):: x(n), dv(n), v(n), mu(n)
        rho = 1.293
        nu = 13.28
        muo = rho* nu
        k = 0.41
        demf = 1
        do i = 1, n
                l = k* x(i)
                nut = demf* l** 2* dv(i)
                mut = nut* rho
                mu(i) = muo + mut
        end do
end subroutine
end program
