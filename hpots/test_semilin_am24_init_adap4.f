c     a testing driver
c     which tests the correctness of
c     subroutine heat2d_am24
c
c
      program test
      implicit real*8 (a-h,o-z)
      integer norder, ndim, iprec
      integer ntot, mxltree, mxboxes, mxlevels
      parameter(ndim=2)
      parameter(norder=8)
      parameter(mxboxes=100000)
      parameter(mxltree=200000)
      parameter(mxlevels=30)
      parameter(ntarg=1000000)
      parameter(nx=1000)
      parameter(ny=1000)
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 u0(norder**ndim,mxboxes)
      real*8 u1(norder**ndim,mxboxes)
      real*8 un(norder**ndim,mxboxes)
      real*8 usol(norder**ndim,mxboxes)
      real*8 uext(norder**ndim,mxboxes)
      real*8 usolp(ntarg)
      real*8 rintl(0:200), dpars(1000)
      real *8, allocatable :: xref(:,:),wts(:)
      real*8, allocatable:: targ(:)
      complex*16 zpars(10), zk
      character*1 type
      external finit, feval, uexact
      open(166,file='self_convergence_test.data')
c
c      dt=1.0d-5
      dt=1.0d-3
      dt=1.0d0/2**10
c     ntot >=2 to begin with
c     fix possible problems later
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0

      ntot=1
      ntot=2**10
c      ntot=ntot*2
c      ntot=ntot*2
c      ntot=ntot*2
c      ntot=ntot*2
c      ntot=ntot*2
c      ntot=ntot*2
c      ntot=ntot*2
c      ntot=ntot*2
c      ntot=ntot*2
c
      tf=ntot*dt
c
      eps=1.0d-9
      npbox=norder**ndim
      nd=1
c
      nordert=4
c      nordert=3
c      nordert=4
c     don't change the tree within fgt for now
      ifnewtree=0

      ifexact = 0

      iadap=1
c
      do i = 3,7
     
      dt=1.0d-2/2**i
      ntot=2**i

      call heat2d_am24_self(norder, nordert, finit, feval,
     1     dt, ntot, eps, mxltree, mxboxes, mxlevels,
     2     ltree, itree, iptr, centers, nlevels,
     3     boxsize, nboxes, ifnewtree, ifexact,iadap,
     4     usol,ntarg,nx,ny,usolp)
c
       do k = 1,ntarg
           write(166,*)usolp(k)
       enddo

      enddo



      end program

c------------------------------------------
c     the forcing term
      subroutine feval(t,x,y,u,fu)
      implicit real*8 (a-h,o-z)
      real*8 x, y, t, u
      real*8 uex,uex1,uex2
      real *8  c1(2),c2(2)
      real *8 del,u1,u2,uu
c
      pi=3.1415926535897932d0
ccc      fu=dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)*dexp(-t)
ccc      fu=(-1.0d0+8.0d0*pi**2)*fu
ccc      fu=(-1.0d0+8.0d0*pi**2)*u
c
      del=1.0d-2
      c1(1)=0.25d0*dcos(20.0d0*pi*t)
      c1(2)=0.25d0*dsin(20.0d0*pi*t)

      c2(1)=0.25d0*dcos(40.0d0*pi*t+pi)
      c2(2)=0.25d0*dsin(40.0d0*pi*t+pi)
c
      uu=0.0d0
      do i = 0,4
        do j = 0,4
            u1=(x-c1(1)-i+2)**2+(y-c1(2)-j+2)**2
            u2=(x-c2(1)-i+2)**2+(y-c2(2)-j+2)**2
            uu=uu+dexp(-u1/del)-0.5d0*dexp(-u2/del)
        enddo
      enddo
      
      
      fu=uu

      end subroutine





c     as rhsfun in test_boxfgt*.f
      subroutine uexact(nd,xyz,dpars,zpars,
     1           ipars,u)
      implicit real *8 (a-h,o-z)
      integer nd, ndim, ipars(*)
      real *8 dpars(*), xyz(*), t, u(nd)
      complex *16 zpars(*)
c
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
c
      nd=1
      t=dpars(1)
c
      u=0.0d0
c
      end subroutine

c     as rhsfun in test_boxfgt*.f
      subroutine finit(nd,xyz,dpars,zpars,
     1           ipars,f)
c     the initial condition
c     defined to be a single Fourier mode
      implicit real *8 (a-h,o-z)
      integer nd, ndim, ipars(*)
      real *8 dpars(*),f(nd),xyz(*)
      real *8 s1, s2, res,del,u1
      complex *16 zpars(*)
c
      nd=1
      ndim=2
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
      del=1.0d-1

      f(1)=dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)
      
      end subroutine
