      subroutine secant(rfun,t,x,y,u0,um1,g,alpha,u,fu,tol,msi,
     1           isi,ndamp,ier)
C
C     Damped secant method for solving
C
C           u = g + alpha * f(u,x,t)
C
C     INPUT:
C
C     A subroutine of the form
C           subroutine rfun(t,x,u,f)
c           real *8 t,x(2)
C           real *8 u,f
C     which computes the value f(u,x,t) 
C
C     u0,um1 = two initial guesses
C     alpha = coefficient of nonlinear term above
C     tol = tolerance factor for solver
C     msi = maximum number of iterates
C
C     OUTPUT:
C
C     u   = solution
C     fu  = f(u,x,t)
C     isi  = number of iterates used
C     ier  = error code.
C          ier = 0 => solution achieved within max
C                     number of secant steps 
C          ier = 1 => max number of secant steps exceeded
C
      implicit none
      integer nnodes,i,msi,nhalf,isi,j,ndamp,ier
      real *8 x,y,pi,done,tol,tes,tes2,alpha
      real *8 u0,sol,resid,hstep,u,fu,t
      real *8 fprime,sol0,solm1,um1,fum1,s0,sm1,g
      external rfun,rjac
C
      done = 1.0d0
      pi=4*datan(done)
      nhalf = ndamp+1
      ier = 0
C
C     sol -> current guess, solm1 -> previous guess
C
      sol = u0
      solm1 = um1
c      write(*,*) 'sol,solm1=',sol,solm1
      call rfun(t,x,y,solm1,fum1)
      do 200 i = 1,msi
	 isi = i
	 sol0 = sol
c       write(*,*) ' '
c       write(*,*) 'sol0=', sol0
	 call rfun(t,x,y,sol0,fu)
c       write(*,*) 'fu=', fu
c       write(*,*) t, x, y, sol0, alpha, fu
	 resid = g - sol0 + alpha*fu
	 tes = abs(resid)
c         write(*,*) ' inside secant, tol = ',tol
c         write(*,*) ' inside secant, tes = ',tes
C
C-----   if tolerance achieved, we are done.
C
	 if (tes .le. tol) goto 201
C
C-----   compute Newton iterate
C
         fprime = (fu-fum1)/(sol0-solm1)
C         
C        now that we have fprime, set current guess to 
C        old guess      
C
         fum1 = fu
         solm1 = sol0
C
ccc         write(6,*)  ' inside secant, fprime = ',fprime
	 hstep =  resid/(1 - alpha*fprime)
ccc         write(6,*)  ' inside secant, hstep = ',hstep
c--------
c       some debugging code
c       if(((x+0.43810045998739905d0) .lt. 1d-14) .and.
c     1    ((y+0.37439954001260095d0) .lt. 1d-14)) then
c         write(*,*) fu, fum1, sol0, solm1, hstep
c       endif
c--------
c
c-----  damping step (check that residual has decreased )
c
	 do 100 j = 1,nhalf
ccc            write(6,*)  ' damping step j = ',j
	     sol = sol0 + hstep
c
c-----      check that residual has decreased
c
	    call rfun(t,x,y,sol,fu)
	    resid = g - sol + alpha*fu
	    tes2 = abs(resid)
ccc            write(6,*)  ' tes,tes2 = ',tes,tes2
	    if (tes2 .lt. tol) then
	       goto 201
	    else if (tes2 .lt. tes) then
	       goto 101
	    else
	       hstep = 0.5D0*hstep
	    endif
100      continue
ccc	 write(6,*) ' exceeded max damping steps ',nhalf
ccc	 write(6,*) ' tes2 = ',tes2
101      continue
ccc         write(6,*)  ' number of damping steps = ',j
200   continue
c
c     looped msi times without success
c
      ier = 1
ccc      write(6,*) ' exceeded max secant steps *',msi
ccc      write(6,*) ' tes2  ',tes2
201   continue
ccc      write(6,*)  ' number of secant steps = ',it
      u = sol
c
c      write(127,*) u0, u
c
      return
      end
C






c------------------------------------------------
c     a vanilla newton's method for the solution
c     of nonlinear systems
c------------------------------------------------
c
      subroutine newton0(nsys, fun, dfun, x0, rtol, atol,
     1           kmax, ipars, dpars, zpars, x1, f1, res, 
     2           ktot, ierr, istat)
      implicit real*8 (a-h, o-z)
      integer nsys, kmax, ktot, ierr
      integer ipars(1)
      real*8 x0(nsys), x1(nsys), dx(nsys), f1(nsys)
      real*8 df(nsys,nsys)
      real*8 cond, tol, dpars(1)
c     local vars      
      integer ipvt(nsys), job
      real*8 rcond, z(nsys)
      complex*16 zpars(1)
      external fun, dfun
c
      done = 1.0d0
      pi=4*datan(done)
c
c     initialization
      ktot = 0
      istat = 0
c      write(*,*) 'nsys=', nsys
c      write(*,*) 'initial value:'
      do i=1, nsys
        x1(i)=x0(i)
c        write(*,*) i, x1(i)
      enddo
c
      do while((ktot .le. kmax) .and. (istat .eq. 0))
c        write(*,*) ' '
c        write(*,*) 'k=', ktot
        ktot = ktot + 1
        call fun(nsys, x1, ipars, dpars, zpars, f1)
        call dfun(nsys, x1, ipars, dpars, zpars, df)
c        write(*,*) 'f1='
c        write(*,*) f1(1)
c        write(*,*) f1(2)
c
c        write(*,*) 'df='
c        write(*,*) df(1,1), df(1,2)
c        write(*,*) df(2,1), df(2,2)
c
        dnorm = 0.0d0
        do i=1, nsys
          dnorm = dnorm + f1(i)**2
        enddo
        dnorm = sqrt(dnorm)
        res = dnorm
c        write(*,*) 'residual='
c        write(*,*) ktot, res
c
c       now call dgeco followed by dgesl to solve 
c       a linear system to iterate
c       dx = df^{-1} f(x1)
c       x1 = x1 - dx
        call dgeco(df,nsys,nsys,ipvt,rcond,z)
        call dgesl(df,nsys,nsys,ipvt,f1,job)
c       the solution saved in f1
c        write(*,*) 'rcond=', rcond
c
c       decide if it's enough to terminate
        dnorm = 0.0d0
        xnorm = 0.0d0
        do i=1, nsys
          dnorm = dnorm + f1(i)**2
          xnorm = xnorm + x1(i)**2
        enddo
        dnorm = sqrt(dnorm)
        xnorm = sqrt(xnorm)
        if(xnorm .gt. 1.0d-10) then
          drel = dnorm/xnorm
        endif
c
c       istat = 1: absolute tolerance reached
c       istat = 2: relative tolerance reached
c       istat = 3: jacobian is ill-conditioned
        if(dnorm .lt. atol) istat = 1
        if(drel .lt. rtol) istat = 2
        if(rcond .lt. 1.0d-9) istat = 3
c
        if(istat .ne. 3) then
c          write(*,*) 'x1='
          do i=1, nsys
            x1(i) = x1(i) - f1(i)
c            write(*,*) i, x1(i)
          enddo
        endif
c       one step of newton's iteration done
c
      enddo
c
      if(res .gt. atol) ierr=1
c
      end subroutine





********************************
* dgeco, dgeco -> dgesl
* zgeco, zgeco -> zgesl
* all come from linpack
* except dasum.f from blas
********************************
      subroutine dgeco(a,lda,n,ipvt,rcond,z)
      integer lda,n,ipvt(1)
      double precision a(lda,1),z(1)
      double precision rcond
c
c     dgeco factors a double precision matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, dgefa is slightly faster.
c     to solve  a*x = b , follow dgeco by dgesl.
c     to compute  inverse(a)*c , follow dgeco by dgesl.
c     to compute  determinant(a) , follow dgeco by dgedi.
c     to compute  inverse(a) , follow dgeco by dgedi.
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       double precision(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack dgefa
c     blas daxpy,ddot,dscal,dasum
c     fortran dabs,dmax1,dsign
c
c     internal variables
c
      double precision ddot,ek,t,wk,wkm
      double precision anorm,s,dasum,sm,ynorm
      integer info,j,k,kb,kp1,l
c
c
c     compute 1-norm of a
c
      anorm = 0.0d0
      do 10 j = 1, n
         anorm = dmax1(anorm,dasum(n,a(1,j),1))
   10 continue
c
c     factor
c
      call dgefa(a,lda,n,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
c     trans(a)  is the transpose of a .  the components of  e  are
c     chosen to cause maximum local growth in the elements of w  where
c     trans(u)*w = e .  the vectors are frequently rescaled to avoid
c     overflow.
c
c     solve trans(u)*w = e
c
      ek = 1.0d0
      do 20 j = 1, n
         z(j) = 0.0d0
   20 continue
      do 100 k = 1, n
         if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
         if (dabs(ek-z(k)) .le. dabs(a(k,k))) go to 30
            s = dabs(a(k,k))/dabs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         if (a(k,k) .eq. 0.0d0) go to 40
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
         go to 50
   40    continue
            wk = 1.0d0
            wkm = 1.0d0
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
               sm = sm + dabs(z(j)+wkm*a(k,j))
               z(j) = z(j) + wk*a(k,j)
               s = s + dabs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               do 70 j = kp1, n
                  z(j) = z(j) + t*a(k,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
c     solve trans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + ddot(n-k,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 110
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
      ynorm = 1.0d0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call daxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 130
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = v
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (dabs(z(k)) .le. dabs(a(k,k))) go to 150
            s = dabs(a(k,k))/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (a(k,k) .ne. 0.0d0) z(k) = z(k)/a(k,k)
         if (a(k,k) .eq. 0.0d0) z(k) = 1.0d0
         t = -z(k)
         call daxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
c     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end
      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision da,dx(1)
      integer i,incx,m,mp1,n,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n .lt. 1 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end


      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),b(1)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end


      double precision function dasum(n,dx,incx)
c
c     takes the sum of the absolute values.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision dx(1),dtemp
      integer i,incx,m,mp1,n,nincx
c
      dasum = 0.0d0
      dtemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 continue
      dasum = dtemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))
     *  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp
      return
      end





**************************************************************
c     Complex version
**************************************************************
c
c subroutines and functions from linpack,
c except zdotc, dzasum, izamax from BLAS
c double complex version
c
c zgeco -> zgesl
c-----------------------------------------

      subroutine zgeco(a,lda,n,ipvt,rcond,z)
      integer lda,n,ipvt(1)
      complex*16 a(lda,1),z(1)
      double precision rcond
c
c     zgeco factors a complex*16 matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, zgefa is slightly faster.
c     to solve  a*x = b , follow zgeco by zgesl.
c     to compute  inverse(a)*c , follow zgeco by zgesl.
c     to compute  determinant(a) , follow zgeco by zgedi.
c     to compute  inverse(a) , follow zgeco by zgedi.
c
c     on entry
c
c        a       complex*16(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       complex*16(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack zgefa
c     blas zaxpy,zdotc,zdscal,dzasum
c     fortran dabs,dmax1,dcmplx,dconjg
c
c     internal variables
c
      complex*16 zdotc,ek,t,wk,wkm
      double precision anorm,s,dzasum,sm,ynorm
      integer info,j,k,kb,kp1,l
c
      complex*16 zdum,zdum1,zdum2,csign1
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
      csign1(zdum1,zdum2) = cabs1(zdum1)*(zdum2/cabs1(zdum2))
c
c     compute 1-norm of a
c
      anorm = 0.0d0
      do 10 j = 1, n
         anorm = dmax1(anorm,dzasum(n,a(1,j),1))
   10 continue
c
c     factor
c
      call zgefa(a,lda,n,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  ctrans(a)*y = e .
c     ctrans(a)  is the conjugate transpose of a .
c     the components of  e  are chosen to cause maximum local
c     growth in the elements of w  where  ctrans(u)*w = e .
c     the vectors are frequently rescaled to avoid overflow.
c
c     solve ctrans(u)*w = e
c
      ek = (1.0d0,0.0d0)
      do 20 j = 1, n
         z(j) = (0.0d0,0.0d0)
   20 continue
      do 100 k = 1, n
         if (cabs1(z(k)) .ne. 0.0d0) ek = csign1(ek,-z(k))
         if (cabs1(ek-z(k)) .le. cabs1(a(k,k))) go to 30
            s = cabs1(a(k,k))/cabs1(ek-z(k))
            call zdscal(n,s,z,1)
            ek = dcmplx(s,0.0d0)*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = cabs1(wk)
         sm = cabs1(wkm)
         if (cabs1(a(k,k)) .eq. 0.0d0) go to 40
            wk = wk/dconjg(a(k,k))
            wkm = wkm/dconjg(a(k,k))
         go to 50
   40    continue
            wk = (1.0d0,0.0d0)
            wkm = (1.0d0,0.0d0)
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
               sm = sm + cabs1(z(j)+wkm*dconjg(a(k,j)))
               z(j) = z(j) + wk*dconjg(a(k,j))
               s = s + cabs1(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               do 70 j = kp1, n
                  z(j) = z(j) + t*dconjg(a(k,j))
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dzasum(n,z,1)
      call zdscal(n,s,z,1)
c
c     solve ctrans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + zdotc(n-k,a(k+1,k),1,z(k+1),1)
         if (cabs1(z(k)) .le. 1.0d0) go to 110
            s = 1.0d0/cabs1(z(k))
            call zdscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0d0/dzasum(n,z,1)
      call zdscal(n,s,z,1)
c
      ynorm = 1.0d0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call zaxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (cabs1(z(k)) .le. 1.0d0) go to 130
            s = 1.0d0/cabs1(z(k))
            call zdscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0d0/dzasum(n,z,1)
      call zdscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = v
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (cabs1(z(k)) .le. cabs1(a(k,k))) go to 150
            s = cabs1(a(k,k))/cabs1(z(k))
            call zdscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (cabs1(a(k,k)) .ne. 0.0d0) z(k) = z(k)/a(k,k)
         if (cabs1(a(k,k)) .eq. 0.0d0) z(k) = (1.0d0,0.0d0)
         t = -z(k)
         call zaxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
c     make znorm = 1.0
      s = 1.0d0/dzasum(n,z,1)
      call zdscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end





c-------------------------------------------
      subroutine zgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      complex*16 a(lda,1)
c
c     zgefa factors a complex*16 matrix by gaussian elimination.
c
c     zgefa is usually called by zgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for zgeco) = (1 + 9/n)*(time for zgefa) .
c
c     on entry
c
c        a       complex*16(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that zgesl or zgedi will divide by zero
c                     if called.  use  rcond  in zgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas zaxpy,zscal,izamax
c     fortran dabs
c
c     internal variables
c
      complex*16 t
      integer izamax,j,k,kp1,l,nm1
c
      complex*16 zdum
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = izamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (cabs1(a(l,k)) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -(1.0d0,0.0d0)/a(k,k)
            call zscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call zaxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (cabs1(a(n,n)) .eq. 0.0d0) info = n
      return
      end







c-------------------------------------------
      subroutine zgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job
      complex*16 a(lda,1),b(1)
c
c     zgesl solves the complex*16 system
c     a * x = b  or  ctrans(a) * x = b
c     using the factors computed by zgeco or zgefa.
c
c     on entry
c
c        a       complex*16(lda, n)
c                the output from zgeco or zgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from zgeco or zgefa.
c
c        b       complex*16(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  ctrans(a)*x = b  where
c                            ctrans(a)  is the conjugate transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if zgeco has set rcond .gt. 0.0
c        or zgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call zgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call zgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas zaxpy,zdotc
c     fortran dconjg
c
c     internal variables
c
      complex*16 zdotc,t
      integer k,kb,l,nm1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call zaxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call zaxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  ctrans(a) * x = b
c        first solve  ctrans(u)*y = b
c
         do 60 k = 1, n
            t = zdotc(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/dconjg(a(k,k))
   60    continue
c
c        now solve ctrans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + zdotc(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end







c-------------------------------------------
      subroutine zdscal(n,da,zx,incx)
*     .. scalar arguments ..
      double precision da
      integer incx,n
*     ..
*     .. array arguments ..
      double complex zx(*)
*     ..
*
*  purpose
*  =======
*
*     zdscal scales a vector by a constant.
*
*  further details
*  ===============
*
*     jack dongarra, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. local scalars ..
      integer i,nincx
*     ..
*     .. intrinsic functions ..
      intrinsic dcmplx
*     ..
      if (n.le.0 .or. incx.le.0) return
      if (incx.eq.1) then
*
*        code for increment equal to 1
*
         do i = 1,n
            zx(i) = dcmplx(da,0.0d0)*zx(i)
         end do
      else
*
*        code for increment not equal to 1
*
         nincx = n*incx
         do i = 1,nincx,incx
            zx(i) = dcmplx(da,0.0d0)*zx(i)
         end do
      end if
      return
      end






c-------------------------------------------
      subroutine zaxpy(n,za,zx,incx,zy,incy)
*     .. scalar arguments ..
      double complex za
      integer incx,incy,n
*     ..
*     .. array arguments ..
      double complex zx(*),zy(*)
*     ..
*
*  purpose
*  =======
*
*     zaxpy constant times a vector plus a vector.
*
*  further details
*  ===============
*
*     jack dongarra, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. local scalars ..
      integer i,ix,iy
*     ..
*     .. external functions ..
      double precision dcabs1
      external dcabs1
*     ..
      if (n.le.0) return
      if (dcabs1(za).eq.0.0d0) return
      if (incx.eq.1 .and. incy.eq.1) then
*
*        code for both increments equal to 1
*
         do i = 1,n
            zy(i) = zy(i) + za*zx(i)
         end do
      else
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         ix = 1
         iy = 1
         if (incx.lt.0) ix = (-n+1)*incx + 1
         if (incy.lt.0) iy = (-n+1)*incy + 1
         do i = 1,n
            zy(iy) = zy(iy) + za*zx(ix)
            ix = ix + incx
            iy = iy + incy
         end do
      end if
*
      return
      end






c-------------------------------------------
      double complex function zdotc(n,zx,incx,zy,incy)
*     .. scalar arguments ..
      integer incx,incy,n
*     ..
*     .. array arguments ..
      double complex zx(*),zy(*)
*     ..
*
*  purpose
*  =======
*
*  zdotc forms the dot product of a vector.
*
*  further details
*  ===============
*
*     jack dongarra, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. local scalars ..
      double complex ztemp
      integer i,ix,iy
*     ..
*     .. intrinsic functions ..
      intrinsic dconjg
*     ..
      ztemp = (0.0d0,0.0d0)
      zdotc = (0.0d0,0.0d0)
      if (n.le.0) return
      if (incx.eq.1 .and. incy.eq.1) then
*
*        code for both increments equal to 1
*
         do i = 1,n
            ztemp = ztemp + dconjg(zx(i))*zy(i)
         end do
      else
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         ix = 1
         iy = 1
         if (incx.lt.0) ix = (-n+1)*incx + 1
         if (incy.lt.0) iy = (-n+1)*incy + 1
         do i = 1,n
            ztemp = ztemp + dconjg(zx(ix))*zy(iy)
            ix = ix + incx
            iy = iy + incy
         end do
      end if
      zdotc = ztemp
      return
      end







c-------------------------------------------
      double precision function dzasum(n,zx,incx)
*     .. scalar arguments ..
      integer incx,n
*     ..
*     .. array arguments ..
      double complex zx(*)
*     ..
*
*  purpose
*  =======
*
*     dzasum takes the sum of the absolute values.
*
*  further details
*  ===============
*
*     jack dongarra, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. local scalars ..
      double precision stemp
      integer i,nincx
*     ..
*     .. external functions ..
      double precision dcabs1
      external dcabs1
*     ..
      dzasum = 0.0d0
      stemp = 0.0d0
      if (n.le.0 .or. incx.le.0) return
      if (incx.eq.1) then
*
*        code for increment equal to 1
*
         do i = 1,n
            stemp = stemp + dcabs1(zx(i))
         end do
      else
*
*        code for increment not equal to 1
*
         nincx = n*incx
         do i = 1,nincx,incx
            stemp = stemp + dcabs1(zx(i))
         end do
      end if
      dzasum = stemp
      return
      end





c-------------------------------------------
      double precision function dcabs1(z)
*     .. scalar arguments ..
      double complex z
*     ..
*     ..
*  purpose
*  =======
*
*  dcabs1 computes absolute value of a double complex number 
*
*  =====================================================================
*
*     .. intrinsic functions ..
      intrinsic abs,dble,dimag
*
      dcabs1 = abs(dble(z)) + abs(dimag(z))
      return
      end






c-------------------------------------------
      integer function izamax(n,zx,incx)
*     .. scalar arguments ..
      integer incx,n
*     ..
*     .. array arguments ..
      double complex zx(*)
*     ..
*
*  purpose
*  =======
*
*     izamax finds the index of element having max. absolute value.
*
*  further details
*  ===============
*
*     jack dongarra, 1/15/85.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. local scalars ..
      double precision dmax
      integer i,ix
*     ..
*     .. external functions ..
      double precision dcabs1
      external dcabs1
*     ..
      izamax = 0
      if (n.lt.1 .or. incx.le.0) return
      izamax = 1
      if (n.eq.1) return
      if (incx.eq.1) then
*
*        code for increment equal to 1
*
         dmax = dcabs1(zx(1))
         do i = 2,n
            if (dcabs1(zx(i)).gt.dmax) then
               izamax = i
               dmax = dcabs1(zx(i))
            end if
         end do
      else
*
*        code for increment not equal to 1
*
         ix = 1
         dmax = dcabs1(zx(1))
         ix = ix + incx
         do i = 2,n
            if (dcabs1(zx(ix)).gt.dmax) then
               izamax = i
               dmax = dcabs1(zx(ix))
            end if
            ix = ix + incx
         end do
      end if
      return
      end







c-------------------------------------------
      subroutine zscal(n,za,zx,incx)
*     .. scalar arguments ..
      double complex za
      integer incx,n
*     ..
*     .. array arguments ..
      double complex zx(*)
*     ..
*
*  purpose
*  =======
*
*     zscal scales a vector by a constant.
*
*  further details
*  ===============
*
*     jack dongarra, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. local scalars ..
      integer i,nincx
*     ..
      if (n.le.0 .or. incx.le.0) return
      if (incx.eq.1) then
*
*        code for increment equal to 1
*
         do i = 1,n
            zx(i) = za*zx(i)
         end do
      else
*
*        code for increment not equal to 1
*
         nincx = n*incx
         do i = 1,nincx,incx
            zx(i) = za*zx(i)
         end do
      end if
      return
      end



c--------------------------------------------------------------
c     newton's method for the solution of nonlinear systems
c     of the form
c
c     u = g+alpha*f(t,x,u)
c
c     where u, f, g are vectors of length nsys
c     x is a vec os length ndim
c
c     input:
c     fevaln: the subroutine of f
c     devaln: jacobian, or the inverse of jacobian, of f
c     u0: the initial guess
c     rtol, atol: relative and absolute tolerance
c     kmax: maximum number of iterations
c
c     output:
c     u1: the solution 
c     f1: f(t,x,u1)
c
c--------------------------------------------------------------

      subroutine newton(nsys, ndim, fevaln, devaln, u0, 
     1           t, xy, g, alpha, ipars, dpars, zpars,
     2           rtol, atol, kmax, u1, f1, res, ktot, 
     3           ierr, istat)
      implicit real*8 (a-h, o-z)
      integer ndim, nsys, kmax, ktot, ierr, ifinv
      integer ipars(1)
      real*8 u0(nsys), u1(nsys), du(nsys), f1(nsys)
      real*8 df(nsys,nsys), alpha, g(nsys)
      real*8 ff(nsys), dff(nsys,nsys), xy(ndim)
      real*8 cond, tol, dpars(1)
c     local vars      
      integer ipvt(nsys), job
      real*8 rcond, z(nsys)
      complex*16 zpars(1)
      external fevaln, devaln
c
c
      done = 1.0d0
      pi=4*datan(done)
c
c     initialization
      ktot = 0
      istat = 0
      ierr = 0
c      write(*,*) 'nsys=', nsys
c      write(*,*) 'initial value:'
      do i=1, nsys
        u1(i)=u0(i)
c        write(*,*) i, u1(i)
      enddo
c
      do while((ktot .le. kmax) .and. (istat .eq. 0))
c        write(*,*) ' '
c        write(*,*) 'k=', ktot
        ktot = ktot + 1
c       I suppose it is yx instead of x1
        call fevaln(nsys,t,xy,u1,dpars,zpars,ipars,
     1       f1)
        call devaln(nsys,t,xy,u1,dpars,zpars,ipars,
     1       df)
c       the capital F
        do i=1, nsys
          ff(i)=u1(i)-alpha*f1(i)-g(i)
        enddo
c       the jacobian of F
        do j=1, nsys
        do i=1, nsys
          dff(i,j)=-alpha*df(i,j)
        enddo
        enddo
c
c        write(301,*) xy(1), xy(2), u1(1), u1(2),
c     1               df(1,1), df(1,2),
c     1               df(2,1), df(2,2)
c
        do i=1, nsys
          dff(i,i)=dff(i,i)+1.0d0
        enddo
c
c        write(302,*) dff(1,1), dff(1,2),
c     1               dff(2,1), dff(2,2)
c
c       decide if it's enough to terminate
        dnorm = 0.0d0
        xnorm = 0.0d0
        do i=1, nsys
          dnorm = dnorm + ff(i)**2
          xnorm = xnorm + u1(i)**2
        enddo
        dnorm = sqrt(dnorm)
        xnorm = sqrt(xnorm)
c
        res = dnorm
c
        if(xnorm .gt. 1.0d-10) then
          drel = dnorm/xnorm
        else
          drel = 1.0d0
        endif
c
        call dgeco(dff,nsys,nsys,ipvt,rcond,z)
        call dgesl(dff,nsys,nsys,ipvt,ff,job)
c
c       istat = 1: absolute tolerance reached
c       istat = 2: relative tolerance reached
c       istat = 3: jacobian is ill-conditioned
        if(dnorm .lt. atol) istat = 1
c        if(drel .lt. rtol) istat = 2
        if(rcond .lt. 1.0d-9) istat = 2
c
       if(istat .ne. 3) then
c          write(*,*) 'x1='
          do i=1, nsys
            u1(i) = u1(i) - ff(i)
c            write(*,*) i, x1(i)
          enddo
        endif
c       one step of newton's iteration done
c
      enddo
c
      if(res .gt. atol) ierr=1
c
c      write(125,*) u0(1), u1(1), res, rcond
c      write(126,*) u0(2), u1(2), res

      end subroutine





c
c     a testing newton solver
c     with the same interface as secant
c     some inputs are redundant
c
c     a newton solver that solves
c     u = g +alpha * f(t,x,u)
c
      subroutine newton_scal(rfun,rjac,t,x,y,u0,um1,g,alpha,u,fu,tol,
     1           msi,ktot,istat,ier)
      implicit real*8 (a-h,o-z)
      integer nnodes,i,msi,nhalf,isi,j,ndamp,ier
      integer ktot, istat
      real *8 x,y,pi,done,tol,tes,tes2,alpha
      real *8 u0,sol,resid,hstep,u,fu,t, ff
      real *8 fprime,sol0,solm1,um1,fum1,s0,sm1,g
      external rfun, rjac
C
      done = 1.0d0
      pi = 4*datan(done)
c
c     u0: initial guess
c     u: the solution, fu: f(t,x,u)
      ktot = 0
      istat = 0
      ier = 0
c
      u = u0
      do while((ktot .le. msi) .and. (istat .eq. 0))
        ktot = ktot + 1
c       call rfun -> fu
        call rfun(t,x,y,u,fu)
c       call rjac -> df
        call rjac(t,x,y,u,df)
c        write(311,*) x,y,u,df
c
        ff = u - alpha*fu -g 
        dff = 1.0d0 -alpha* df
        res = abs(ff)
c
        dnorm = abs(ff)
        drel = abs(ff)/abs(u)
        rcond = 1.0d0/abs(dff)
c
c       istat = 1: absolute tolerance reached
c       istat = 2: relative tolerance reached
c       istat = 3: jacobian is ill-conditioned
        if(dnorm .lt. tol) istat = 1
c        if(drel .lt. tol) istat = 2
        if(rcond .lt. 1.0d-9) istat = 2
c
c       attention: u and fu inconsistent here
c       modify later
        ff = ff/dff
        if(istat .ne. 3) then
          u = u -ff
        endif
c
      enddo
c
      if(res .gt. tol) ier = 1
c      write(128,*) u0, u, res
 


      end subroutine

