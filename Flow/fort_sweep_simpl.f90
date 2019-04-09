program main
implicit none

integer i
integer, parameter:: n=10
real(8):: aa,bb,y(n+1),x(n+1),h,alf(n+1),&
       bet(n+1),a(n+1),b(n+1),c(n+1),f(n+1)

open(1, file="out.txt")
aa=0.0
bb=1.0
h=(bb-aa)/n

open(2, file="mat.txt")
read(2,*) a(:)
read(2,*) b(:)
read(2,*) c(:)
read(2,*) f(:)


alf(1)=-c(1)/b(1)
bet(1)=f(1)/b(1)

do i=2,n
alf(i)=-c(i)/(a(i)*alf(i-1)+b(i))
bet(i)=(f(i)-a(i)*bet(i-1))/&
        (a(i)*alf(i-1)+b(i))
end do

y(n+1)=1

do i=n,2,-1
y(i)=alf(i)*y(i+1)+bet(i)
end do

x(1)=0
do i=2,n+1
x(i)=x(i-1)+h
end do
do i=1,n+1
write(1,*) x(i), y(i)
end do
close(1)
close(2)
contains
!ne pashet prosto sohranil
subroutine deriv(v, h, dv, n)
        integer:: n, i
        real(8):: v(n), h(n-1), dv(n), delt
        
        delt = h(2)/ h(1)
        dv(1) = 1/ (h(1) + h(2))* &
        (-(2 + delt)* v(1) + (1 + delt)** 2 / delt * v(2) - v(3)/ delt)
        
        do i = 2, n-1
                delt = h(i)/ h(i-1)
                dv(i) = 1/ (h(i+1) + h(i))* &
                (-delt* v(i-1) + (delt**2 - 1)/ delt* v(i) + 1/ delt* v(i+1))
        end do  
               
        delt = h(n-1)/ h(n-2)
        dv(n) = 1/ (h(n) + h(n-1))* &
        (delt* v(n-2) - (1 + delt)** 2/ delt* v(n-1) + (2 + delt)/ delt* v(n))
end subroutine
end program main
