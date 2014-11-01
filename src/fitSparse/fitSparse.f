	subroutine fitSparse(A,n,m,LDA,istart,Ap,Ai,Ax)
	implicit none
	integer n,m,LDA, istart
	real *8 A(LDA,m)
	integer Ap(1), Ai(1)
	real *8 Ax(1)

	integer i,j,k
	real *8 eps
	parameter (eps = 1e-16)
	integer pp, pi, px

	

	
	

	

