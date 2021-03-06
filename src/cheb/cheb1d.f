    
      subroutine cheb1d(x,y,n,m,dm,p,U,VT,G,LDU,LDV,LDG)
c     /////////////////////////////////////////////////////////////
c     n, m: Input. Scalar. Length of x and y
c     x, y: Input. Vector of length n and m. Points where kernel is
c           evaluated.
c     p: Input. Scalar. Rank
c     dm: Input. Vector of length 4. Represent the domain
c           where x and y lies in. dm has the structure
c                       [ax,bx,ay,by],
c           where ax <= x <= bx and ay <= y <= by.
c     LDU, LDV,LDG: Input. Scalar. Leading dimension of U and VT. 
c           LDU >= n; LDV >= p; LDG >=p.
c     U: Output. Matrix of size (n * p). L2L matrix.
c     VT: Output. Matrix of size (p * m). M2M matrix.
c     G: Output. Matrix of size (p * p).  M2L matrix.
c     
c     U, G, and VT construct a low rank approximation of 
c                       K ~ U*G*VT
c           where K is the kernel matrix. i.e. kernel evaluated at
c           x and y.
c
c      
c
c
c
c     
c
c     /////////////////////////////////////////////////////////////
      implicit none
      integer n,m,p,LDU, LDV,LDG
      real *8 x(n), y(m), U(LDU,p), VT(LDV, m), G(LDG, p)
      real *8 dm(4)
      
      integer i,j,k,l,ctype
      real *8 xcheb(p), ycheb(p),wx(p),wy(p)
      real *8 denoX(n),denoY(m)     

      ctype = 1

c     Obtain p chebyshev nodes for x
      call cheb_nodes(p-1,dm(1),dm(2),xcheb,wx,ctype)

c     Obtain p chebyshev nodes for y
      call cheb_nodes(p-1,dm(3),dm(4),ycheb,wy,ctype)

c     Evaluate kernel at chebyshev nodes
      call evaKer(G,xcheb,ycheb,p,p,p)

c     Compute U
c     First, precompute the denominator in Barycentric Formula
      do i = 1,n
            denoX(i) = 0
            do k = 1,p
                  denoX(i) = denoX(i) + wx(k)/(x(i) - xcheb(k))
            enddo !k
      enddo !i

      do k = 1,p
            do i = 1,n
                  U(i,k) = wx(k)/(x(i)-xcheb(k))/denoX(i)
            enddo
      enddo
c     Compute VT
c     Precompute denominator
      do j = 1,m
            denoY(j) = 0
            do k = 1,p
                  denoY(j) = denoY(j) + wy(k)/(y(j) - ycheb(k))
            enddo !k
      enddo !i
      
      do j = 1,m
            do k = 1,p
                  VT(k,j) = wy(k)/(y(j)-ycheb(k))/denoY(j)
            enddo
      enddo


      return

      end



      subroutine cheb_nodes(n,a,b,x,w,ctype)
c     /////////////////////////////////////////////////////////////
c     generate n+1 Chebyshev nodes and weights in interval [a,b]
c     /////////////////////////////////////////////////////////////
      implicit none
      integer n,ctype
      real *8 x(n+1), w(n+1),a,b
      integer j
      
      if (ctype == 1) then
        call std_cheb_nodes1(n,x,w)
      else
        call std_cheb_nodes2(n,x,w)
      endif
      
      do j = 1,n+1
            x(j) = ((b-a)*x(j)+(b+a))/2
      enddo

      return

      end




      subroutine std_cheb_nodes2(n,x,w)
c     /////////////////////////////////////////////////////////////
c     generate n+1 Chebyshev nodes of 2nd kind and corresponding
c     weights in interval [-1,1]
c     n+1: number of nodes (by convention)
c     x: Chebyshev nodes of second kind, which is defined as
c                 x_j = cos((j-1)*PI/n)   for j = 1,..,n+1
c     w: Chebyshev weights, which is defined as
c                 w_j = (-1)^(j-1)/2      for j = 1 or j = n+1
c                 w_j = (-1)^(j-1)        for j = 2,..,n
c     
c     /////////////////////////////////////////////////////////////
      implicit none

      real *8 PI, done
      integer n
      real *8 x(n+1), w(n+1)
      integer j,s

      done = 1.0d0

      PI = 4.d0*atan(done)
      
      do j = 1,n+1
            x(j) = cos((j-1)*PI/n)
      enddo

      s = -1
      w(1) = 0.5

      do j = 2,n
            w(j) = s*1.0
            s = -s
      enddo

      w(n+1) = s*0.5

      return
      end

      subroutine std_cheb_nodes1(n,x,w)
c     /////////////////////////////////////////////////////////////
c     generate n+1 Chebyshev nodes of 1st kind and corresponding
c     weights in interval [-1,1]
c     n+1: number of nodes (by convention)
c     x: Chebyshev nodes of second kind, which is defined as
c                 x_j = cos((j-1)*PI/n)   for j = 1,..,n+1
c     w: Chebyshev weights, which is defined as
c                 w_j = (-1)^(j-1)/2      for j = 1 or j = n+1
c                 w_j = (-1)^(j-1)        for j = 2,..,n
c     /////////////////////////////////////////////////////////////
      implicit real *8 (a-h,o-z)
      implicit integer (i-n)

      real *8 x(n+1), w(n+1)

      done = 1.d0
      PI = 4.d0*atan(done)

      s = 1.d0
      do j = 1,n+1
        x(j) = cos((2*j-1)*PI/(2*n+2))
        w(j) = s*sin((2*j-1)*PI/(2*n+2))
        s = -s
      enddo

      return
      end


