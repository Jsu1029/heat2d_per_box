c     a testing driver test11
c     which tests the correctness of
c     subroutine ns2dperfm4(2d periodic forcing term order4)
c
c     this version tests the correctness
c     of the solver for an adaptive case
c     with a known solution
c
c     see subroutines finit, fforce, uexact
c     for more details
c
c
      program test
      implicit real*8 (a-h,o-z)
      integer norder, ndim, iprec,nx,ny
      integer ntot, mxltree, mxboxes, mxlevels,ntarg
      parameter(ndim=2)
      parameter(nd=2)
      parameter(norder=8)
      parameter(mxboxes=100000)
      parameter(mxltree=200000)
      parameter(mxlevels=30)
      parameter(ntarg=1000000)
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
c
      real*8,allocatable :: velocity(:,:,:)
      real*8,allocatable :: vorticity(:,:)
      real*8 velocityp(nd,ntarg)
      real*8,allocatable :: velocity_ex(:,:,:)
      real*8 rintl(0:200), dpars(1000)
      real *8 visc
      complex*16 zpars(10)
      real *8 daprs(10)
      
      real*8 errpe
      character*1 type
      external finit, fforce, uexact
c
      pi=3.1415926535897932d0
c
      dt=1.2d0/1200
      ntot=1200
      ipoly=1

      allocate(velocity_ex(nd,norder**ndim,mxboxes))
      allocate(velocity(nd,norder**ndim,mxboxes))
      allocate(vorticity(norder**ndim,mxboxes))
c
      write(*,*)'*********',ntot
      tf=ntot*dt
c
      eps=1.0d-9

      visc=1.0d-3
      nordert=4
c
      ifnewtree=1
c
      irhsfun=0
c
      call ns2dperfm4(nd,norder,nordert,finit,
     1      fforce,dt, ntot, eps, mxltree, mxboxes, mxlevels,
     2     ltree, itree, iptr, centers, nlevels,
     3     boxsize, nboxes, ifnewtree,velocity,vorticity,
     4     ntarg,velocityp,visc,irhsfun)

      end program





c------------------------------------------
c     should have the same calling seq
c     as rhsfun in test_boxfgt*.f
      subroutine finit(nd,xyz,dpars,zpars,
     1           ipars,f)
c     the initial condition
c     defined to be a single Fourier mode
      implicit real *8 (a-h,o-z)
      integer nd, ipars(*)
      real *8 dpars(*),f(nd),xyz(*)
      complex *16 zpars(*)
c
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
      rho=30.0d0
      delta=5.0d-2
c
      if (y .le. 0.0)then
        f(1) = tanh(rho*(y+0.25d0))
      else
        f(1) = tanh(rho*(0.25d0-y))
      endif

      f(2)=delta*sin(2.0d0*pi*x)
      
      end subroutine






c------------------------------------------
c     the forcing term
      subroutine fforce(nd,xyz,dpars,zpars,ipars,u)
      implicit real*8 (a-h,o-z)
      integer nd, ipars(*)
      real *8 dpars(*), xyz(*), t, u(nd)
      complex *16 zpars(*)
c
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
      t=dpars(1)
c
        u(1) =0.0d0
c
        u(2) =0.0d0

      end subroutine






c     as rhsfun in test_boxfgt*.f
      subroutine funs(nd,xyz,dpars,zpars,
     1           ipars,f)
c     the initial condition
c     defined to be a single Fourier mode
      implicit real *8 (a-h,o-z)
      integer nd, ipars(*)
      real *8 dpars(*),f(nd),xyz(*)
      complex *16 zpars(*)
c
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
      t=dpars(1)
      
      f(1)=0.0d0

      f(2)= 0.0d0
      end subroutine


c---------------------------------------
c     should have the same calling seq
c     as rhsfun in test_boxfgt*.f
ccccc just for the program to run properly
ccccc not true exact solution
ccccc this test uses self convergence to get L2 error
      subroutine uexact(nd,xyz,dpars,zpars,
     1           ipars,u)
      implicit real *8 (a-h,o-z)
      integer nd, ipars(*)
      real *8 dpars(*), xyz(*), t, u(nd)
      complex *16 zpars(*)
c
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
      t=dpars(1)
      u(1)=dsin(2.0d0*pi*x+dexp(t))*dcos(2.0d0*pi*y)
      u(2)=-dcos(2.0d0*pi*x+dexp(t))*dsin(2.0d0*pi*y)
c
      end subroutine
