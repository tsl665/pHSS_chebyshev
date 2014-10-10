	subroutine cheb_interp2D(n,z,x,y,f)




	integer n
	real *8 x,y,z,f(n)
	integer j
	

	










	end
	
	
	
	
	subroutine cheb_interp1D(n,y,x,f)
c	/////////////////////////////////////////////////////////////
c	n: scalar. Number of points using in interpolation
c	f: vector of length n. Function value at chebyshev points
c	x: input
c	y: output value of interpolated function
c
c
c
c	
c
c	/////////////////////////////////////////////////////////////
	
	integer n
	real *8 x,y,f(n)
	integer j
	real *8 z(n),w(n),deno

	call std_cheb_nodes(n-1,z,w)
	deno = 0
	do j = 1,n
		y = y + w(j)*f(j)/(x - z(j))
		deno = deno + w(j)/(x - z(j))
	enddo
	
	y = y/deno


	return

	end




	subroutine cheb_nodes(n,x,w,a,b)
c	/////////////////////////////////////////////////////////////
c	generate n+1 Chebyshev nodes and weights in interval [a,b]
c	/////////////////////////////////////////////////////////////
	integer n
	real *8 x(n+1), w(n+1)
	integer j

	call std_cheb_nodes(n,x,w)
	do j = 1,n+1
		x(j) = ((b-a)*x(j)+(b+a))/2
	enddo

	return

	end




	subroutine std_cheb_nodes(n,x,w)
c	/////////////////////////////////////////////////////////////
c	generate n+1 Chebyshev nodes and weights in interval [-1,1]
c	n+1: number of nodes (by conversion)
c	x: Chebyshev nodes of second kind, which is defined as
c			x_j = cos((j-1)*PI/n)	for j = 1,..,n+1
c	w: Chebyshev weights, which is defined as
c			w_j = (-1)^(j-1)/2	for j = 1 or j = n+1
c			w_j = (-1)^(j-1)		for j = 2,..,n
c	
c	/////////////////////////////////////////////////////////////
	real *8 PI
	parameter (PI = 3.141592653589793)
	integer n
	real *8 x(n+1), w(n+1)
	integer j,s
	s = -1
	do j = 1,n+1
		x(j) = cos((j-1)*PI/n)
	enddo

	w(1) = 1/2

	do j = 2,n
		w(j) = s
		s = s*(-1)
	enddo

	w(n+1) = s/2
	
	return
	end

