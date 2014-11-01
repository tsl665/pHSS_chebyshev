	subroutine cheb_interp1D(x,y,n,m,p,U,VT,LDU,LDV)
c	/////////////////////////////////////////////////////////////
c	n, m: Input. Scalar. Length of x and y
c	x, y: Input. Vector of length n and m. Points where kernel is
c		evaluated.
c	p: Input. Scalar. Rank
c	LDU, LDV: Input. Scalar. Leading dimension of U and VT. 
c		LDU >= p, LDV >= m.
c	U, VT: Output. Matrix of size (n * p) and (p * m). Low rank
c		approximation of kernel evaluated at x and y.
c
c	 
c
c
c
c	
c
c	/////////////////////////////////////////////////////////////
	
	integer n,m,p,LDU, LDV
	real *8 x(n), y(m), U(n,LDU), VT(p, LDV)
	
	integer i,j,k
	real *8 xcheb(p), ycheb(p),wx(p),wy(p)
	real *8 xmax, ymax, xmin, ymin

	xmin = x(1)
	xmax = x(n)
	ymin = y(1)
	ymax = y(m)
c	Obtain p chebyshev nodes for x
	call cheb_nodes(p-1,xmin,xmax,xcheb,wx)

c	Obtain p chebyshev nodes for y
	call cheb_nodes(p-1,ymin,ymax,ycheb,wy)

c	Evaluate kernel at chebyshev nodes
	call evaKer(G,xcheb,ycheb,p,p)

c	Compute U
c	First, precompute the denominator in Barycentric Formula
	do i = 1,n
		denoX(i) = 0
		do k = 1,p
			denoX(i) = denoX(i) + wx(k)/(x(i) - xcheb(k))
		enddo !k
	enddo !i

	do k = 1,p
		do i = 1,n
			U(i,k) = wx(i)/(x(i)-xcheb(k))/denoX(i)
		enddo
	enddo

c	Compute VT
c	Precompute denominator
	do j = 1,m
		denoY(j) = 0
		do k = 1,p
			denoY(j) = denoY(j) + wy(k)/(x(j) - xcheb(k))
		enddo !k
	enddo !i

	do j = 1,m
		do k = 1,p
			VT(k,j) = 0
			do l = 1,p
			VT(k,j) = VT(k,j)+G(k,l)*wy(j)/(y(j)-ycheb(l))
			enddo !l
			VT(k,j) = VT(k,j)/denoY(j)
		enddo !k
	enddo !j
	
	return

	end



	subroutine cheb_nodes(n,a,b,x,w)
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
c	n+1: number of nodes (by convention)
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


