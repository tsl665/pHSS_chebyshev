	subroutine treeCIR(A,LDA,n,p,work)
c	////////////////////////////////////////////////////////////
c	
c	treeCIR = (construct) TREE (using) Chebyshev Interpolation
c		(to a specific) Rank
c
c	Construct Tree Structure for HODLR matrix to a specific
c	rank, using Chebyshev interpolation of continuous kernel.
c
c	A: Input. Matrix (n*n). HODLR matrix
c	LDA: Input. Scalar. Leading dimension of A.
c	n: Input. Scalar. Size of A.
c	p: Input. Scalar. Desired rank.
c
c	work: Output. Data structure stores the tree.
c	Ap,Ai,Ax: Output. Vector.
c	lp: Output. Scalar. Length of Ap.
c
c	////////////////////////////////////////////////////////////

	implicit none
	integer n,p,LDA
	real *8 A(LDA,*), Ax(*)
	integer Ap(*), Ai(*), work(*)
c	work: the object stores the tree structure of HODLR matrix.
c	work(1): starting index of tree structure. The first 1 to 
c		work(1)-1 memories are for constants. Here work(1) = 31.
c	
c	work(2) - work(20) is constants about attribution of each
c	single node. Let "tag" be the information we want. "pw"
c	is pointer for current node. Then we access that information
c	by work(pw+tag). "tag" is one of the following:
c
c	work(2): [lvl] Level of current node
c	work(3): [mstart] The starting column(row) of current node
c		in A
c	work(4): [msize] Size of current node
c	work(5): [odsci] Off-Diagonal Starting Column Index
c		Notice that the [odsri] is the same as [mstart]
c
c	work(6): [prt] Pointer of parent of current node. e.g. parent
c		of node specified by pointer "pw" is 
c		work(work(pw + parent))
c		If this value is 0, the current node is root.
c	work(7),work(8): [lc], [rc]: Pointer of left/right child of
c		current node. If both these two value is 0, the current
c		node is leaf.
c		If lc = -1, the left child of current node is
c		not yet specified.
c	work(9), [leaf]: work(pw+leaf) = 0 means not leaf
c		work(pw+leaf) = 1 means leaf
c
c
c	work(10) - work(20) is reserved.
c
c	work(21) - work(30) is other constants.
c
c	work(21): [kappa]: finest(max) level. Given by
c		kappa = INT(log_2(n/2/p))
c		In practice, kappa = INT(log(n*done/2/p)/log(2.0)+eps)
c	work(22): [sns]: Single Node Size. Number of memories need to
c		store a single node.
c		sns = number of memories used from work(2) to work(30)
c	work(23): [WE]: End pointer for "work".
c	work(24): [n]
c	work(25): [p]
c	work(26): [nn]: Size of extended matrix
c


	integer pw
	
	call initTree(work,n,p)

	pw = work(1)

	do
		if (pw >= work(23)) exit

		if (work(pw+work(7)) == -1) then
			call newNode(A,LDA,work, pw, 'L')
		endif

		if (work(pw+work(8)) == -1) then
			call newNode(A, LDA, work, pw, 'R')
		endif

		pw = pw + work(22) ! pw = pw + sns

	enddo

	return
	end
	
	subroutine newNode(A,LDA,work,pw,child)
	implicit none
	integer LDA, pw, work(*)
	real *8 A(LDA,*)
	character child

	integer pw1, pw2, i,j,k, nL, nR, nS
	integer kappa, n, p

	pw1 = pw !pointer of parent
	pw2 = work(23) !WE
	kappa = work(21)
	n = work(24)
	p = work(25)
	nL = INT(work(pw1 + work(4))/2) ! left child size
	nR = work(pw1 + work(4)) - nL   ! right child size
	nS = work(pw1 + work(3))        ! parent starting index

c	Update independent attribution
	work(pw2 + work(2)) = work(pw1 + work(2)) + 1 !lvl
	work(pw2 + work(6)) = pw1 !prt

c	Update attribution depending on level
	if (work(pw2 + work(2)) < kappa) then
		work(pw2 + work(7)) = -1 !lc
		work(pw2 + work(8)) = -1 !rc
		work(pw2 + work(9)) = 0  !leaf
	else
		work(pw2 + work(7)) = 0
		work(pw2 + work(8)) = 0
		work(pw2 + work(9)) = 1
	endif
	
c	Update attribution depending on child type
	if ( child == 'L' ) then
		work(pw2 + work(3)) = nS      !mstart
		work(pw2 + work(4)) = nL      !msize
		work(pw2 + work(5)) = nS + nL !odsci
	
		work(pw1 + work(7)) = pw2	!lc of parent
	elseif (child == 'R') then
		work(pw2 + work(3)) = nS + nL !mstart
		work(pw2 + work(4)) = nR	!msize
		work(pw2 + work(5)) = nS	!odsci

		work(pw1 + work(8)) = pw2	!rc of parent
	else
		print *, 'Wrong input "child": Must be "L" or "R"!'
	endif


c	Update pointers

	work(23) = work(23) + work(22) !WE = WE + sns

	return

	end
	

	

	subroutine initTree(work,n,p)

	implicit none
	integer work(:)
	integer i, root, n,p

	real *8 done, eps
	parameter (done = 1.0, eps = 1E-8)

c	Initialize Constants
	work(1) = 31
	do i = 2,9
		work(i) = i - 2
	enddo
	work(21) = INT (log(n*done/2/p)/log(2.0)+eps) !kappa
	work(22) = 8          !sns
	work(23) = work(1) - 1 !WE
	work(24) = n
	work(25) = p
	work(26) = n + 4 * (2**work(21) - 1) * p

c	Initialize root node
	root = work(1)
	work(root + work(2)) = 0  !lvl
	work(root + work(3)) = 1  !mstart
	work(root + work(4)) = n  !msize
	work(root + work(5)) = 0  !odsci(no low-rank part)
	work(root + work(6)) = 0  !prt(no parent)
	work(root + work(7)) = -1 !lc
	work(root + work(8)) = -1 !rc
	work(root + work(9)) = 0  !leaf
	work(root + work(10)) = 0 !iL
	work(root + work(11)) = 0 !iP

	work(23) = root + work(22) ! WE = root + sns
	
	return

	end

	
	
	
	
	

	
	


