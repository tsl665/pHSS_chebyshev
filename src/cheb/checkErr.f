      subroutine getRank1(krank,tol,dm)
      implicit real *8 (a-h,o-z)
      implicit integer (i-n)

      real *8 dm(4)

      krank = 1
      do
        call checkAprxErr1(errMax,krank,dm)
        if (errMax<tol) exit
        krank = krank * 2
      end do

      return

      end


      subroutine checkAprxErr1(errMax,krank,dm)

      implicit real *8 (a-h,o-z)
      implicit integer (i-n)

      real *8 dm(4)
      real *8, dimension (:), allocatable :: x,y
      real *8, dimension (:,:), allocatable :: K,Kapp,G,U,VT
      
      integer AlcSta
      n = krank * 2
      
      allocate (x(n), stat = AlcSta)
      IF ( AlcSta/= 0) stop "*** Not Enough Memory"

      allocate (y(n), stat = AlcSta)
      IF ( AlcSta/= 0) stop "*** Not Enough Memory"

      allocate (K(n,n), stat = AlcSta)
      IF ( AlcSta/= 0) stop "*** Not Enough Memory"

      allocate (G(krank,krank), stat = AlcSta)
      IF ( AlcSta/= 0) stop "*** Not Enough Memory"

      allocate (U(n,krank), stat = AlcSta)
      IF ( AlcSta/= 0) stop "*** Not Enough Memory"

      allocate (VT(krank,n), stat = AlcSta)
      IF ( AlcSta/= 0) stop "*** Not Enough Memory"

      allocate (Kapp(n,n), stat = AlcSta)
      IF ( AlcSta/= 0) stop "*** Not Enough Memory"
      



      do i = 1,n
        x(i) = (i-1)/(n-1)*(dm(2)-dm(1)) + dm(1)
      enddo

      do j = 1,n
        y(j) = (j-1)/(n-1)*(dm(4)-dm(3)) + dm(3)
      enddo

c     Compute the original matrix
      call evaKer(K,x,y,n,n,n)

c     Compute Approximation matrix Kapp = U*G*VT
      call cheb1d(x,y,n,n,dm,krank,U,VT,G,n,krank,krank)
      Kapp = matmul(U,matmul(G,VT))
      
      errMax = 0
      valMax = 0
      do i = 1,n
        do j = 1,n
          if (abs(K(i,j))>valMax) valMax = abs(K(i,j))
          if (abs(K(i,j)-Kapp(i,j))>errMax) errMax = abs(K(i,j)
     &        -Kapp(i,j))
        enddo
      enddo



c      call materr(n,n,K,Kapp,error,errms)

      deallocate (x,y,K,G,U,VT,Kapp)
  

      end




      subroutine checkRank2(krank,dm)
      implicit real *8 (a-h,o-z)
      implicit integer (i-n)

      integer krank
      real *8 dm(4),a,b,c,d
      real *8 x(krank),y(krank),w(krank)
      real *8 cheCM(krank,krank),G(krank,krank)
      
      PI = 4.d0*atan(1.d0)

      a = dm(1)
      b = dm(2)
      c = dm(3)
      d = dm(4)

      call cheb_nodes(krank-1,a,b,x,w,2)
      
      call cheb_nodes(krank-1,c,d,y,w,2)

c     compute diagonal elements


      call evaKer(G,x,y,krank,krank,krank)
      do i = 1,krank
        do j = 1,krank
          cheCM(i,j) = 0
          do k = 1,krank
            if ((k == 1).OR.(k==krank)) then
              s1 = 0.5d0
            else
              s1 = 1.d0
            endif
            do l = 1,krank
              if ((k == 1).OR.(k==krank)) then
                s2 = 0.5d0
              else
                s2 = 1.d0
              endif
              cheCM(i,j) = cheCM(i,j)+G(i,j)*cos(PI*k*i/krank)
     &          *cos(PI*l*j/krank)*s1*s2
            enddo !l
          enddo !k
          cheCM(i,j) = cheCM(i,j)/krank/krank
        enddo !j
      enddo !i

      dcv = 0
      do i = 1,krank
        j = krank+1-i
        dcv = dcv + cheCM(i,j)
      enddo
      dcv = dcv / krank

      hfcv = 0
      k =0
      do i = krank/2+1,krank
        do j = krank/2+1,krank
          hfcv = hfcv + cheCM(i,j)
          k = k + 1
        enddo
      enddo
      hfcv = hfcv / k

      do i = 1,krank
        call prinf('i = *',i,1)
        call prin2('i-th column = *', cheCM(1,i),krank)
      enddo

      call prin2('dcv = *',dcv,1)
      call prin2('hfcv = *',hfcv,1)

      return

      end



        
        

      
