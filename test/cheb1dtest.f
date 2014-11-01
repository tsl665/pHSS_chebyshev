	implicit none

	integer n,m,i,j,k,p

	real *8, dimension (:), allocatable :: x,y,sigma,f_true
	real *8, dimension (:), allocatable :: f_app,tempVec
	real *8, dimension (:,:), allocatable :: G, U, VT
	real *8 error
	integer AllocateStatus



c	print *, 'Enter n (Length of x)'
c	read *, n
c	print *, 'Enter m (Length of y)'
c	read *, m
	n = 100
	m = 100
	print *, 'Enter p (Rank of Approximation)'
	read *, p

	allocate (x(n), stat = AllocateStatus)
	IF (AllocateStatus /= 0) stop "*** Not Enough Memory"
	allocate (y(m), stat = AllocateStatus)
	IF (AllocateStatus /= 0) stop "*** Not Enough Memory"
	allocate (G(n,m), stat = AllocateStatus)
	IF (AllocateStatus /= 0) stop "*** Not Enough Memory"
	allocate (sigma(m), stat = AllocateStatus)
	IF (AllocateStatus /= 0) stop "*** Not Enough Memory"
	allocate (f_true(n), stat = AllocateStatus)
	IF (AllocateStatus /= 0) stop "*** Not Enough Memory"
	allocate (f_app(n), stat = AllocateStatus)
	IF (AllocateStatus /= 0) stop "*** Not Enough Memory"
	allocate (U(n,p), stat = AllocateStatus)
	IF (AllocateStatus /= 0) stop "*** Not Enough Memory"
	allocate (VT(p,m), stat = AllocateStatus)
	IF (AllocateStatus /= 0) stop "*** Not Enough Memory"
	allocate (tempVec(p), stat = AllocateStatus)
	IF (AllocateStatus /= 0) stop "*** Not Enough Memory"
	
	call prini(6,13)
	
	
	call RANDOM_SEED
	call RANDOM_NUMBER(x)
	call RANDOM_NUMBER(y)
	call RANDOM_NUMBER(sigma)
c	call quickSortD(x,n)
c	call quickSortD(y,m)
c	do i = 1,n
c		x(i) = i/(n+1)
c	enddo
c	do j = 1,m
c		y(j) = j/(m+1)
c	enddo
	call evaKer(G,x,y,n,m)
	f_true = matmul(G,sigma)

	call cheb_interp1D(x,y,n,m,p,U,VT,p,m)
	tempVec = matmul(VT,sigma)
	f_app = matmul(U,tempVec)

	call prin2('f_true = *',f_true,n)

c	call prin2('U = *',U,n*p)
c	call prin2('VT = *',VT,p*m)
c	call prin2('tempVec = *',tempVec,p)
	call prin2('f_app = *',f_app,n)

	error = 0
	do i = 1,n
		error = error + (f_app(i)-f_true(i))*(f_app(i)-f_true(i))
	enddo
	error = sqrt(error)

	call prin2('Error = *',error,1)
	

	end	
