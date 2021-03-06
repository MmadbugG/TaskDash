program separation_var
implicit none
	real(8):: Bi, Fon, eps, h, ht
	real(8), allocatable:: Td(:, :), xh(:), Fo(:)
	integer:: i, j, n, nt
        real(8):: Ta, Te, alp, delt, c, rho, lam, a
        real(8):: fsl(3)=(/0.05, 0.1, 0.5/),&
                  xsl(3)=(/0.0, 0.5, 1.0/)
	real(8), allocatable:: Am(:), Bm(:), Cm(:), Fm(:), ym(:)
	
	alp = 11000
	c = 670
	rho = 2500.0
	lam = 0.74
	Fon = 0.5
	Ta = 20
	Te = 100
	delt = 0.005
	
	Bi = alp* delt/ lam
	a = lam/ (c* rho)
	
	n = 1000 !nodes length, the length of interval is equal to n - 1
	nt = 1000
	h = 1.0/ (n - 1)
        ht = Fon/ (nt - 1)
	allocate(xh(n), Fo(nt), Td(n,nt))
	allocate(Am(n), Bm(n), Cm(n), Fm(n), ym(n))
	eps = 1.0e-6
        do i = 1, n
                xh(i) = (i - 1)* h
        end do
        do i = 1, nt
                Fo(i) = (i - 1)* ht
        end do
        
        Td(:, 1) = 1
        do i = 2, nt
                call matrix(h, ht, n, Td(:, i-1), Am, Bm, Cm, Fm) 
                call sweep(Am, Bm, Cm, Fm, n, ym)
                Td(:, i) = ym
        end do
                
	call units(Td, Fo, xh, n, nt, a, delt, Ta, Te)
	do i = 1, 3
	        fsl(i) = fsl(i)* delt**2 / a
	        xsl(i) = xsl(i)* delt
	end do 
	call export_txt(Td, Bi, Fo, xh, Fon, n, nt)
        call export_fo(xh, Fo, Td, n, nt, fsl)
        call export_xh(xh, Fo, Td, n, nt, xsl)
               

contains
subroutine matrix(h, ht, n, tet, A, B, C, F)
	real(8):: h, ht, tet(n)
	integer:: n
	real(8):: A(n), B(n), C(n), F(n)
	integer:: i

	do i = 2, n-1
		C(i) = - 2/ h** 2 - 1/ ht
		A(i) = 1/ h** 2
		B(i) = 1/ h** 2
		F(i) = - tet(i)/ ht
	end do
	F(1) = 0
	F(n) = 0
	C(1) = - 1/ h
	C(n) = - Bi - 1/ h
	A(1) = 0
	A(n) = 1/ h
	B(1) = 1/ h 
	B(n) = 0
end subroutine
subroutine sweep(A, B, C, F, n, y)
	real(8):: A(n), B(n), C(n), F(n)
	integer:: n
	real(8):: y(n)
	real(8):: alp(n), bet(n)
	integer:: i
	
	alp(2) = - B(1)/ C(1)
	bet(2) = F(1)/ C(1)
	do i = 2, n-1
		alp(i+1) = -B(i)/ (A(i)* alp(i) + C(i))
		bet(i+1) = (F(i) - A(i)* bet(i))/&
        	(A(i)* alp(i) + C(i))
	end do
	y(n) = (F(n) - A(n)* bet(n))/ (C(n) + A(n)* alp(n))
	do i = n-1 , 1, -1
		y(i) = alp(i+1)* y(i+1) + bet(i+1)
	end do
end subroutine
subroutine units(Td, Fo, xh, n, nt, a, delt, Ta, Te)
        real(8):: Ta, Te, a, delt
        real(8):: Td(n, nt), xh(n), Fo(nt)
        integer:: n, nt, i, j
        
        do j = 1, nt
                Fo(j) = delt**2 / a* Fo(j)
        end do
        do i = 1, n
                do j = 1, nt
                        Td(i,j) = Td(i,j)* (Ta - Te) + Te
                end do
                xh(i) = xh(i)* delt
        end do
end subroutine
subroutine export_txt(Td, Bi, Fo, xh, Fon, n, nt)
	real(8):: Td(n, nt), Bi, Fo(nt), Fon, xh(n)
	integer:: n, nt
	real(8):: T(nt, n), time(nt)	
	integer:: i
	character(len=15):: str
	open(1, file='file_out.txt')
	write(1,'(a40/)') 'Created from separation_var program'
	write(1,'(a20/)') 'Parameter: '
	write(1,'(a6,f12.1)') 'Bi = ', Bi
	write(1,'(a6,f12.1/)') 'Fo = ', Fon
	write(1,'(a20,i4)') 'Coordinate step: ', n
	write(1,'(a20,i4/)') 'Time step: ', nt
	write(1,*) ('-', i=1,30 )
	write(str , *) n
	write(1,'(16x,'//str//'f15.7/)') xh(:)
	do i = 1, nt
		write(1,'(i5,f15.7,a3,'//str//'f15.7)') &
		i, Fo(i), ' | ', Td(:, i)
	end do
	close(1)	
end subroutine
integer function slice(xh, xsl, n)
        real(8):: xh(n)
        integer:: isl, n
        real(8):: dist, buff, xsl
        integer:: i
        buff = abs(xh(n)-xh(1)) + 1
        do i = 1, n
                dist = abs(xh(i) - xsl)
                if (dist > buff) then
                        if (i == 1) then
                                write(*,*) 'slice bad...'
                                stop
                        end if
                        exit
                end if
                buff = dist
        end do
        isl = i - 1
        slice = isl
end function
subroutine export_fo(xh, Fo, Td, n, nt, fsl)
        real(8):: xh(n), Fo(nt), Td(n, nt), fsl(3)
        integer:: n, nt
        integer:: isl(3), i
        open(1, file='exp_fo.txt')
        write(1,'(a40/)') 'Created from separation_var program'
	write(1,'(a20/)') 'Parameter: '
	write(1,'(a6,f12.1)') 'Bi = ', Bi
	write(1,'(a6,f12.1/)') 'Fo = ', Fon
	write(1,'(a20,i4)') 'Coordinate step: ', n
	write(1,'(a20,i4/)') 'Time step: ', nt
	write(1,*) ('-', i=1,30 )
	
	write(1, *) xh(:)
        do i = 1, 3
                isl(i) = slice(Fo, fsl(i), nt)
                write(1, *) Td(:, isl(i))
        end do
        close(1)
end subroutine
subroutine export_xh(xh, Fo, Td, n, nt, xsl)
        real(8):: xh(n), Fo(nt), Td(n, nt), xsl(3)
        integer:: n, nt
        integer:: isl(3), i
        open(1, file='exp_xh.txt')
        write(1,'(a40/)') 'Created from separation_var program'
	write(1,'(a20/)') 'Parameter: '
	write(1,'(a6,f12.1)') 'Bi = ', Bi
	write(1,'(a6,f12.1/)') 'Fo = ', Fon
	write(1,'(a20,i4)') 'Coordinate step: ', n
	write(1,'(a20,i4/)') 'Time step: ', nt
	write(1,*) ('-', i=1,30 )
	
        write(1, *) Fo(:)
        do i = 1, 3
                isl(i) = slice(xh, xsl(i), n)
                write(1, *) Td(isl(i), :)
        end do
        
        close(1)
end subroutine
end program
