	subroutine evaKer(G,x,y,n,m,LDG)
c	/////////////////////////////////////////////////////////////
c	Compute kernel K(x,y) evaluate at x and y, where x and y are
c	vectors of lenth n and m.
c	
c	n,m: Input. Scalar. Length of x and y.
c	x, y: Input. Vector
c	LDG: Input. Scalar. Leading dimension of G. LDG >= n
c	G: n*m matrix. Kernel evaluate at (x,y).
c
c
c
c	
c
c	/////////////////////////////////////////////////////////////
	integer n, m,LDG
	real *8 x(n), y(m), G(LDG,m)
	integer i,j

	do j = 1,m
		do i = 1,n
			call ker1(G(i,j),x(i),y(j))
		enddo !i
	enddo !j

	return
	
	end



	subroutine ker1(G,x,y)
c	/////////////////////////////////////////////////////////////
c	Compute kernel K(x,y) evaluate at x and y.
c	
c	x, y: scalar. Input
c	G: scalar. Kernel evaluate at (x,y).
c
c
c
c	
c
c	/////////////////////////////////////////////////////////////

	real *8 G,x,y

c	G = exp((x-y)*(x-y))
c	G = cos(x/3)*sin(2*y-1)
c	G = 3
c	G = (3*x+2)*(4*y-4)
c	G = (2*x*x+4*x+3)*(3*y*y+4)
	G = (x+2.5)*(x-1.5)*(x-2)*(y+2)*(y-2)*(y-2)
c	G = (x-4)*(x-2)*(x-3)*(x+4)*(y+3)*(y-2)*(y+2.5)*(y-2)
	return

	end

c	subroutine quickSortD(a,n)
c	implicit none
c	integer n
c	real *8 a(n),swap

c	integer i,j
c
c	do i = 1,n-1
c		do j = i+1,n
c			if (a(i)>a(j)) then
c				swap = a(i)
c				a(i) = a(j)
c				a(j) = swap
c			endif
c		enddo
c	enddo
c
c	return
c	end

