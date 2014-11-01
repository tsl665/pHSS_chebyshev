	subroutine evaKer(G,x,y,n,m)
c	/////////////////////////////////////////////////////////////
c	Compute kernel K(x,y) evaluate at x and y, where x and y are
c	vectors of lenth n and m.
c	
c	n,m: scalar. Length of x and y.
c	x, y: vector. Input
c	G: n*m matrix. Kernel evaluate at (x,y).
c
c
c
c	
c
c	/////////////////////////////////////////////////////////////
	integer n, m
	real *8 x(n), y(m), G(n,m)
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

	G = exp((x-y)*(x-y))

	return

	end
