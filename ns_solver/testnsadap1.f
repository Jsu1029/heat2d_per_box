c     a testing driver testns7
c     which tests the correctness of
c     subroutine heat2dperfm
c
c     this version tests the correctness
c     of the solver for an adaptive case
c     with a known solution
c
c     see subroutines finit, fforce, uexact
c     for more details
c
      program test
      implicit real*8 (a-h,o-z)
      integer norder, ndim, iprec,nx,ny
      integer ntot, mxltree, mxboxes, mxlevels,ntarg
      integer nordert
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
c
      real*8 errpe
      character*1 type
      external finit, fforce, uexact
      open(220,file='usoltargetx2.txt')
c
      pi=3.1415926535897932d0
c
      dt=1.0d0/100
      ntot=100
      ipoly=1

      allocate(velocity_ex(nd,norder**ndim,mxboxes))
      allocate(velocity(nd,norder**ndim,mxboxes))
      allocate(vorticity(norder**ndim,mxboxes))
c
      write(*,*)'*********',ntot
      tf=ntot*dt
c
      eps=1.0d-9

      visc=1.0d0
      nordert=4
c
      ifnewtree=1
c
      irhsfun=1
c
      call ns2dperfm4(nd,norder,nordert,finit,
     1      fforce,dt, ntot, eps, mxltree, mxboxes, mxlevels,
     2     ltree, itree, iptr, centers, nlevels,
     3     boxsize, nboxes, ifnewtree,velocity,vorticity,
     4     ntarg,velocityp,visc,irhsfun)

c    velocity at target points
      do i=1,ntarg
         write(220,*)velocityp(1,i)
      enddo

      end program





c---------------------------------------
c     should have the same calling seq
c     as rhsfun in test_boxfgt*.f
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


c----------------------------------



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
        u(1) =dexp(t)*dcos(2.0d0*pi*x+dexp(t))*dcos(2.0d0*pi*y)+
     1      8.0d0*pi**2*dsin(2.0d0*pi*x+dexp(t))*dcos(2.0d0*pi*y)+
     2     2.0d0*pi*dcos(2.0d0*pi*x)*dcos(2.0d0*pi*y)*dexp(t)+
     3     2.0d0*pi*dsin(2.0d0*pi*x+dexp(t))*dcos(2.0d0*pi*y)*
     4     dcos(2.0d0*pi*x+dexp(t))*dcos(2.0d0*pi*y)+2.0d0*pi*
     5     dcos(2.0d0*pi*x+dexp(t))*dsin(2.0d0*pi*y)*
     6     dsin(2.0d0*pi*x+dexp(t))*dsin(2.0d0*pi*y)

c
       u(2) =dexp(t)*dsin(2.0d0*pi*x+dexp(t))*dsin(2.0d0*pi*y)
     1   -8.0d0*pi**2*dcos(2.0d0*pi*x+dexp(t))*dsin(2.0d0*pi*y)-
     2   2.0d0*pi*dsin(2.0d0*pi*x)*dsin(2.0d0*pi*y)*dexp(t)+
     3   2.0d0*pi*dsin(2.0d0*pi*x+dexp(t))*dcos(2.0d0*pi*y)*
     4   dsin(2.0d0*pi*x+dexp(t))*dsin(2.0d0*pi*y)+2.0d0*pi*
     5   dcos(2.0d0*pi*x+dexp(t))*dsin(2.0d0*pi*y)*
     6   dcos(2.0d0*pi*x+dexp(t))*dcos(2.0d0*pi*y)

      end subroutine
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
      
      f(1)= dsin(2.0d0*pi*x+1.0d0)*dcos(2.0d0*pi*y)
      f(2)=-dcos(2.0d0*pi*x+1.0d0)*dsin(2.0d0*pi*y)
      
      end subroutine

c------------------------------------------
c     should have the same calling seq
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
      
      f(1)=dexp(t)*dcos(2.0d0*pi*x+dexp(t))*dcos(2.0d0*pi*y)+
     1      8.0d0*pi**2*dsin(2.0d0*pi*x+dexp(t))*dcos(2.0d0*pi*y)

      f(2)= dexp(t)*dsin(2.0d0*pi*x+dexp(t))*dsin(2.0d0*pi*y)
     1   -8.0d0*pi**2*dcos(2.0d0*pi*x+dexp(t))*dsin(2.0d0*pi*y)
      end subroutine

    
