      implicit real *8 (a-h,o-z)
      implicit integer (i-n)

      real *8 dm(4)
      real *8, dimension (:), allocatable :: x,y
      real *8, dimension (:,:), allocatable :: K,Kapp,G,U,VT
      
      integer AllocateStatus
      n = 10
      m = 10
      print *, 'Enter p (Rank of Approximation)'
      read *, krank
      
      allocate (x(n), stat = AllocateStatus)
      IF (AllocateStatus /= 0) stop "*** Not Enough Memory"

      allocate (y(m), stat = AllocateStatus)
      IF (AllocateStatus /= 0) stop "*** Not Enough Memory"

      allocate (K(n,m), stat = AllocateStatus)
      IF (AllocateStatus /= 0) stop "*** Not Enough Memory"

      allocate (G(krank,krank), stat = AllocateStatus)
      IF (AllocateStatus /= 0) stop "*** Not Enough Memory"

      allocate (U(n,krank), stat = AllocateStatus)
      IF (AllocateStatus /= 0) stop "*** Not Enough Memory"

      allocate (VT(krank,m), stat = AllocateStatus)
      IF (AllocateStatus /= 0) stop "*** Not Enough Memory"

      allocate (Kapp(n,m), stat = AllocateStatus)
      IF (AllocateStatus /= 0) stop "*** Not Enough Memory"
      

      call prini(6,13)

      dm = [-1,1,-1,1]

      print *, 'Predicted Error: '
      call checkRank(krank,dm)
  
  
      call RANDOM_SEED
      call RANDOM_NUMBER(x)
      do i = 1,n
        x(i) = x(i)*(dm(2)-dm(1)) + dm(1)
      enddo
      call RANDOM_NUMBER(y)
      do j = 1,m
        y(j) = y(j)*(dm(4)-dm(3)) + dm(3)
      enddo

c     Compute the original matrix
      call evaKer(K,x,y,n,m,n)

c     Compute Approximation matrix Kapp = U*G*VT
      call cheb1d(x,y,n,m,dm,krank,U,VT,G,n,krank,krank)
      Kapp = matmul(U,matmul(G,VT))

      call materr(n,m,K,Kapp,error,errms)


      call prin2('errmax = *',error,1)
      call prin2('err_root-mean-square = *',errms,1)


      deallocate (x,y,K,G,U,VT,Kapp)
  

      end
      

