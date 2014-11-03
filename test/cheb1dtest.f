	implicit none

	integer n,m,i,j,p

	real *8, dimension (:), allocatable :: x,y,sigma,f_true
	real *8, dimension (:), allocatable :: f_app,tempVec
	real *8, dimension (:,:), allocatable :: K, U, G, VT, Kapp
	real *8 error, dm(4), errorM, normf, normM,errms
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

	allocate (U(n,p), stat = AllocateStatus)
	IF (AllocateStatus /= 0) stop "*** Not Enough Memory"

	allocate (VT(p,m), stat = AllocateStatus)
	IF (AllocateStatus /= 0) stop "*** Not Enough Memory"

	allocate (Kapp(n,m), stat = AllocateStatus)
	IF (AllocateStatus /= 0) stop "*** Not Enough Memory"

	
	call prini(6,13)

	dm = [-1,3,2,5]
	
	
	call RANDOM_SEED
	call RANDOM_NUMBER(x)
	do i = 1,n
		x(i) = x(i)*(dm(2)-dm(1)) + dm(1)
	enddo
	call RANDOM_NUMBER(y)
	do j = 1,m
		y(j) = y(j)*(dm(4)-dm(3)) + dm(3)
	enddo

c	Compute the original matrix
	call evaKer(K,x,y,n,m,n)

c	Compute Approximation matrix Kapp = U*G*VT
	call cheb1d(x,y,n,m,dm,p,U,VT,G,n,p,p)
	Kapp = matmul(U,matmul(G,VT))

	call materr(n,m,K,Kapp,error,errms)
	call prin2('errmax = *',error,1)
	call prin2('err_root-mean-square = *',errms,1)


	deallocate (x,y,K,G,U,VT,Kapp)
	

	end	
