c     a testing driver
c     which tests the correctness of
c     subroutine heat2dperfm
c
c     this version tests the correctness
c     of the solver for an adaptive case
c     with a unknown solution
c
c     see subroutines finit, fforce, uexact
c     for more details
c
c     self convergence
c
      program test
      implicit real*8 (a-h,o-z)
      integer norder, ndim, iprec,nx,ny
      integer ntot, mxltree, mxboxes, mxlevels,ntarg
      parameter(ndim=2)
      parameter(norder=8)
      parameter(mxboxes=100000)
      parameter(mxltree=200000)
      parameter(mxlevels=30)
      parameter(ntarg=1)
      parameter(nx = 1)
      parameter(ny = 1)
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 usol(norder**ndim,mxboxes)
      real*8 usolp(ntarg)
      real*8 uextp(ntarg)
      real*8 rintl(0:200), dpars(1000)
      real *8, allocatable :: xref(:,:),wts(:)
      real*8, allocatable:: targ(:)
      real*8, allocatable:: targs(:)
      real*8 hx,hy
      complex*16 zpars(10), zk
      real *8 daprs(10)
      
      real*8 errpe
      character*1 type
      external finit, fforce, uexact
c
      open(17,file='usolp11.data')

      pi=3.1415926535897932d0
c
      write(*,*)'*********',ntot
      tf=ntot*dt
c
      eps=1.0d-3
      npbox=norder**ndim
c
c      nordert=2
      nordert=4
      ifnewtree=1
c
      allocate(targs(ndim))

c
      do i = 5,5
        
         dt= 1.0d-1/(2**i)
         ntot = 2**i

      call heat2dperfm3(norder, nordert,finit,
     1      fforce,dt, ntot, eps, mxltree, mxboxes, mxlevels,
     2     ltree, itree, iptr, centers, nlevels,
     3     boxsize, nboxes, ifnewtree, usol,ntarg,nx,ny,usolp)
c
       do k = 1,ntarg
           write(17,*)usolp(k)
       enddo

      enddo
c


      end program





c------------------------------------------
c     should have the same calling seq
c     as rhsfun in test_boxfgt*.f
      subroutine finit(nd,xyz,dpars,zpars,
     1           ipars,f)
c     the initial condition
c     defined to be a single Fourier mode
      implicit real *8 (a-h,o-z)
      integer nd, ndim, ipars(*)
      real *8 dpars(*),f(nd),xyz(*)
      real *8 s1, s2, res
      complex *16 zpars(*)
c
      nd=1
      ndim=2
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
c
      f(1)=0.0d0
      
      end subroutine






c------------------------------------------
c     the forcing term
      subroutine fforce(nd,xyz,dpars,zpars,ipars,u)
      implicit real*8 (a-h,o-z)
      integer nd, ndim, ipars(*)
      real *8 dpars(*), xyz(*), t, u(nd)
      complex *16 zpars(*)
      real *8  c1(2),c2(2)
      real *8 del,u1,u2
      
      nd=1
c
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
      t=dpars(1)
      del = 1.0d-2
c
      c1(1)=0.25d0*dcos(20.0d0*pi*t)
      c1(2)=0.25d0*dsin(20.0d0*pi*t)

      c2(1)=0.25d0*dcos(40.0d0*pi*t+pi)
      c2(2)=0.25d0*dsin(40.0d0*pi*t+pi)
c
      u(1)=0.0d0
      do i = 0,4
        do j = 0,4
            u1=(x-c1(1)-i+2)**2+(y-c1(2)-j+2)**2
            u2=(x-c2(1)-i+2)**2+(y-c2(2)-j+2)**2
            u(1)=u(1)+dexp(-u1/del)-0.5d0*dexp(-u2/del)
        enddo
      enddo


      

      end subroutine




c---------------------------------------
c     should have the same calling seq
c     as rhsfun in test_boxfgt*.f
ccccc just for the program to run properly
ccccc not true exact solution
ccccc this test uses self convergence to get L2 error
      subroutine uexact(t,xyz,dpars,zpars,
     1           ipars,u)
      implicit real *8 (a-h,o-z)
      integer nd, ndim, ipars(*)
      real *8 dpars(*), xyz(*), t, u
      complex *16 zpars(*)
c
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
c
      u=dsin(2.0d0*pi*x+dexp(t))*dcos(2.0d0*pi*y)
c
      end subroutine
