	implicit none

	integer n,m,i,j,p

	real *8, dimension (:), allocatable :: x,y,sigma,f_true
	real *8, dimension (:), allocatable :: f_app,tempVec
	real *8, dimension (:,:), allocatable :: K, U, G, VT
	real *8 error, dm(4)
	integer AllocateStatus



c	print *, 'Enter n (Length of x)'
c	read *, n
c	print *, 'Enter m (Length of y)'
c	read *, m
	n = 1000
	m = 1000
	print *, 'Enter p (Rank of Approximation)'
	read *, p

	allocate (x(n), stat = AllocateStatus)
	IF (AllocateStatus /= 0) stop "*** Not Enough Memory"

	allocate (y(m), stat = AllocateStatus)
	IF (AllocateStatus /= 0) stop "*** Not Enough Memory"

	allocate (K(n,m), stat = AllocateStatus)
	IF (AllocateStatus /= 0) stop "*** Not Enough Memory"

	allocate (G(p,p), stat = AllocateStatus)
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

	dm = [0,1,0,1]
	
	
	call RANDOM_SEED
	call RANDOM_NUMBER(x)
	do i = 1,n
		x(i) = x(i)*(dm(2)-dm(1)) + dm(1)
	enddo
	call RANDOM_NUMBER(y)
	do j = 1,m
		y(j) = y(j)*(dm(4)-dm(3)) + dm(3)
	enddo
	call RANDOM_NUMBER(sigma)
c	call quickSortD(x,n)
c	call quickSortD(y,m)
c	do i = 1,n
c		x(i) = i/(n+1)
c	enddo
c	do j = 1,m
c		y(j) = j/(m+1)
c	enddo

	call evaKer(K,x,y,n,m,n)
	f_true = matmul(K,sigma)

	call cheb_interp1D(x,y,n,m,dm,p,U,VT,G,n,p,p)
	tempVec = matmul(VT,sigma)
	tempVec = matmul(G,tempVec)
	f_app = matmul(U,tempVec)

c	call prin2('f_true = *',f_true,n)

c	call prin2('U = *',U,n*p)
c	call prin2('VT = *',VT,p*m)
c	call prin2('tempVec = *',tempVec,p)
c	call prin2('f_app = *',f_app,n)

	error = 0
	do i = 1,n
		error = error + (f_app(i)-f_true(i))*(f_app(i)-f_true(i))
	enddo
	error = sqrt(error)

	call prin2('Error = *',error,1)


	deallocate (x,y,K,G,U,VT,sigma,tempVec,f_true,f_app)
	

	end	
