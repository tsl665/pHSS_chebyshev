	subroutine sparseIndex(work,Ap,Ai,Ax,lp)

	implicit none
	integer work(*), Ap(*), Ai(*), lp
	real *8 Ax(*)

	integer i,j
	integer, dimension (:), allocatable :: Aps,Ais
	real *8, dimension (:), allocatable :: Axs
	integer lps

	kappa = work(21)
	sns = work(22)
	n = work(24)
	p = work(25)
	nn = work(26)

c	First, Construct the diagonal dense part

	pw = (2**kappa-1)*sns+work(1) !the first leaf

	do k = 1, 2**kappa
c		irow and icol are starting indices of diagonal
c		dense matrices in HODLR matrix A. They happen to
c		be the same as starting indicies of near field(dense)
c		blocks in extended matrix.
		irow = work(pw + work(3))
		icol = irow
		nrow = work(pw + work(4))
		if (k == 1) then
			call indexConstruct(A(irow,icol),nrow,nrow,LDA,
     &			Ap,Ai,Ax,lp)
		else
			allocate(Aps(nrow))
			allocate(Ais(nrow*nrow))
			allocate(Axs(nrow*nrow))
			call indexConstruct(A(irow,icol),nrow,nrow,LDA,
     &			Aps,Ais,Axs,lps)
			call indexMerge(Ap,Ai,Ax,lp,Aps,Ais,Axs,lps,
     &			irow,icol)
			deallocate(Aps,Ais,Axs)
		endif
		
		pw = pw + sns
	enddo !k

c	Then, Construct the rest. For each level (from fine to 
c	coarse), first construct S2M, then construct M2L
c	
	do l = kappa, 0,-1
		pw = (2**l-1)*sns+work(1) ! First node in level l

		do k = 1,2**l
			odsri = work(pw + work(3))
			odsci = work(pw + work(5))
ccccc	I need to check with Mike here... Is there any how to deal
c	with the domain of x and y on each subblocks. It's given
c	or we calculate its value by max and min
			call cheb1d()
			
		


	return
	
	end

	subroutine indexIdentity(n,Ap,Ai,Ax,lp)
c	/////////////////////////////////////////////////////////////
c	Construct a input to SuiteSparse representing n*n identity
c	matrix.
c
c	n: Input. Scalar. Size of identity matrix.
c	Ap,Ai,Ax: Output. Vector of Length (lp, li, li). Input for
c		SuiteSparse
c	lp: Output. Scalar. Length of Ap.
c	
c	Remark: length of Ai = Ap(lp)
c
c	/////////////////////////////////////////////////////////////

	implicit none
	integer n, Ap(*), Ai(*), lp
	real *8 Ax(*)
	integer i

	do i = 1,n
		Ap(i) = i
		Ai(i) = i
		Ax(i) = 1.0
	enddo
	Ap(n+1) = n+1
	lp = n+1

	return
	
	end
	

		
	


	subroutine indexConstruct(A,m,n,LDA,Ap,Ai,Ax,lp)
c	/////////////////////////////////////////////////////////////
c	Construct a input to SuiteSparse representing matrix block A.
c
c	A: Input. Matrix of size (n*m). The Matrix block to fit in
c		SuiteSparse.
c	n, m: Input. Scalar. Dimensions of A.
c	LDA: Input. Scalar. Leading dimension of A.
c	Ap,Ai,Ax: Output. Vector of Length (lp, li, li). Input for
c		SuiteSparse
c	lp: Output. Scalar. Length of Ap.
c	
c	Remark: length of Ai = Ap(lp)
c
c	/////////////////////////////////////////////////////////////

	implicit none
	integer n,m,LDA, istart
	real *8 A(LDA,m)
	integer Ap(1), Ai(1),lp
	real *8 Ax(1)

	integer i,j,k
	real *8 eps
	parameter (eps = 1e-15)
	integer pp, pi

	pp = 1
	pi = 1

	do j = 1,n
		Ap(pp) = pi - 1
		pp = pp + 1
		do i = 1,m
			if (abs(A(i,j)>eps)) then
				Ai(pi) = i - 1
				Ax(pi) = A(i,j)
				pi = pi + 1
			endif
		enddo !i
	enddo !j

	Ap(pp) = pi
	lp = pp

	return

	end

	subroutine indexWrap(work,Ap,Ai,lp,srow,scol)
c	/////////////////////////////////////////////////////////////
c	Wrap the index information of submatrix A into a big vector
c
c	work: Output. Vector. Contain all index info of A.
c		work(1): length of work. work(1) = 5 + lp + li
c		work(2): lp
c		work(3): srow -- starting row index
c		work(4): scol -- starting column index
c		work(5:4+lp): Ap
c		work(5+lp:work(1)): Ai
c	n, m: Input. Scalar. Dimensions of A.
c	LDA: Input. Scalar. Leading dimension of A.
c	Ap,Ai: Input. Vector of Length (lp, Ap(lp)). Input for
c		SuiteSparse
c	lp: Input. Scalar. Length of Ap.
c
c	/////////////////////////////////////////////////////////////
	implicit none
	integer work(1), Ap(1), Ai(1), lp,li,srow,scol
	integer pw,i

	work(1) = 5 + lp + Ap(lp)
	work(2) = lp
	work(3) = srow
	work(4) = scol

	pw = 5
	do i = 1,lp
		work(pw) = Ap(i)
		pw = pw + 1
	enddo

	do i = 1,li
		work(pw) = Ai(i)
		pw = pw + 1
	enddo

	IF (pw - 1 /= work(1)) stop "*** Index Mistake in indexWrap"	

	return

	end

	subroutine indexUnwrap(Ap,Ai,lp,srow,scol,work)

c	/////////////////////////////////////////////////////////////
c	Unwrap information from work. Definition of variable is the 
c	same as in indexWrap, with input and output exchanged. 
c
c	/////////////////////////////////////////////////////////////
	implicit none  
	integer work(1), Ap(1), Ai(1), lp,srow,scol
	integer pw,i

	lp = work(2)
	srow = work(3) 
	scol = work(4)  
	pw = 5
	do i = 1,lp
		Ap(i) = work(pw)
		pw = pw + 1
	enddo

	do i = 1,li
		Ai(i) = work(pw)
		pw = pw + 1
	enddo
	
	IF (pw - 1 /= work(1)) stop "*** Index Mistake in indexUnwrap"	

	return

	end

	subroutine indexMerge(Ap,Ai,Ax,lp,Aps,Ais,Axs,lps,srow,scol)

c	/////////////////////////////////////////////////////////////
c	
c	Given total Sparse set Ap, Ai, Ax, and information about a
c	submatrix (subm, Axs). Merge the information into big set.
c	In-place merge. Time complexity is O(n)
c	
c	/////////////////////////////////////////////////////////////


	implicit none
	integer Ap(*),Ai(*),Ax(*),Aps(*),Ais(*),Axs(*)
c	integer, dimension (:), allocatable :: Aps,Ais
	integer lps,srow,scol,lp
	integer i,j,pp,pi,pps,pis,gap,li

c	allocate (Aps(subm(2)))
c	allocate (Ais(subm(3)))
	
c	call indexUnwrap(Aps,Ais,lps,srow,scol,subm)

	pp = scol
	pps = 1
	li = Ap(lp) - 1
	do i = 1,lps
c		Update pointers for pi, by Ap
		pi = Ap(pp)
		do while ((srow<Ai(pi)).AND.(pi<Ap(pp+1)))
			pi = pi + 1
		enddo

c		How many non-zero elements in submatrix in this column
		gap = Aps(pps+1)-Aps(pps)

c		Move Space for Ais(Axs) to merge in
		do j = li,pi,-1
			Ai(j + gap) = Ai(j)
			Ax(j + gap) = Ax(j)
		enddo
c		Update li and Ap accordingly
		li = li + gap
		do j = pp+1,lp
			Ap(j) = Ap(j) + gap
		enddo

c		Fill in Ais(Axs)
		do pis = Aps(pps), Aps(pps+1)-1
			Ai(pi) = Ais(pis) + srow - 1
			Ax(pi) = Axs(pis)
			pi = pi + 1
		enddo

c		Update pointers for Ap and Aps
		pp = pp + 1
		pps = pps + 1
		IF (li /= Ap(lp)-1) stop "*** Error in indexMerge"
	enddo

	return

	end

	subroutine indexF2C(Ap,Ai,lp)

c	/////////////////////////////////////////////////////////////
c	
c	Transfer Fortran style indices (start from 1) to C style
c	indices (start from 0)
c
c	/////////////////////////////////////////////////////////////
	implicit none
	integer Ap(*), Ai(*),lp
	integer li, i

	li = Ap(lp)-1
	do i = 1,lp
		Ap(i) = Ap(i)-1
	enddo
	
	do i = 1,li
		Ai(i) = Ai(i)-1
	enddo

	return

	end
