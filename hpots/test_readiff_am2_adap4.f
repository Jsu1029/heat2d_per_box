c     a testing driver 
c     which tests the correctness of 
c     subroutine heatd_sys_am2
c
c
      program test
      implicit real*8 (a-h,o-z)
      integer norder, ndim, iprec, nd
      integer ntot, mxltree, mxboxes, mxlevels
c     nd: dimension of the system
c     nd: dimension of the xyz param
      parameter(nd=2)
      parameter(ndim=2)
      parameter(norder=8)
      parameter(mxboxes=100000)
      parameter(mxltree=1000000)
      parameter(mxlevels=30)
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8, allocatable :: usol(:,:,:)
c      real*8 usol(2,64,mxboxes)
c      real*8 uext(nd,norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000), visc(2)
      real*8 uu(2)
      
      complex*16 zpars(10), zk
      character*1 type
      external finit, fevaln, devaln, uexact
      open(166,file='test3_am2.data')
c
c      dt=1.0d-5
c      dt=1.0d-3
      dt=10.0d0/1000
      ntot=1000
c
      tf=ntot*dt
c
      eps=1.0d-9
      npbox=norder**ndim
c
      nordert=2
c      nordert=4
c     don't change the tree within fgt for now
      ifnewtree=0
c
      iadap = 1
      visc(1) = 2.0d-3
      visc(2) = 2.0d-3
c     target point
      allocate(usol(nd,norder**ndim,mxboxes))
c
      call heat2d_sys_am2_video(nd, norder, nordert, finit,
     1     fevaln, devaln, visc, dt, ntot, eps, mxltree,
     2     mxboxes, mxlevels, ltree, itree, iptr, 
     3     centers, nlevels, boxsize, nboxes, iadap, usol)



      end program
c
c
c
c
c
c     be careful about the calling seqs
c     of the following functions
c     in the solver
c
c     fevaln is a little different
c     (for the use of secant.f/secant)
c
c-----------------------------------------
c
c     Initial condition
c     u(x,t)=sin(2*pi*x)*cos(2*pi*y)*exp(-t)
c     v(x,t)=sin(4*pi*x)*cos(2*pi*y)exp(-2*t)
c
c     u(x,0)=sin(2*pi*x)*cos(2*pi*y)
c     v(x,0)=sin(4*pi*x)*cos(2*pi*y)
c
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
      ndim=2
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)


c
      if(xyz(1) .lt. 0.5d0 .and. xyz(1) .ge. -0.5d0 .and.
     1    xyz(2) .lt. 0.5d0 .and. xyz(2) .ge. -0.5d0) then
c(0,0)
c
      f(1)=1.0d0+y
      f(2)=3.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. 1.5d0 .and. xyz(1) .ge. 0.5d0 .and.
     1    xyz(2) .lt. 0.5d0 .and. xyz(2) .ge. -0.5d0) then
c(1,0)
c
      f(1)=1.0d0+y
      f(2)=-1.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. 2.5d0 .and. xyz(1) .ge. 1.5d0 .and.
     1    xyz(2) .lt. 0.5d0 .and. xyz(2) .ge. -0.5d0) then
c(2,0)
c
      f(1)=1.0d0+y
      f(2)=-6.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. -0.5d0 .and. xyz(1) .ge. -1.5d0 .and.
     1    xyz(2) .lt. 0.5d0 .and. xyz(2) .ge. -0.5d0) then
c(-1,0)
c
      f(1)=1.0d0+y
      f(2)=8.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. -1.5 .and. xyz(1) .ge. -2.5d0 .and.
     1    xyz(2) .lt. 0.5d0 .and. xyz(2) .ge. -0.5d0) then
c(-2,0)
c
      f(1)=1.0d0+y
      f(2)=13.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. 0.5d0 .and. xyz(1) .ge. -0.5d0 .and.
     1    xyz(2) .lt. 1.5d0 .and. xyz(2) .ge. 0.5d0) then
c(0,1)
c
      f(1)=y
      f(2)=3.5d0+5.0d0*x

      elseif(xyz(1) .lt. 1.5d0 .and. xyz(1) .ge. 0.5d0 .and.
     1    xyz(2) .lt. 1.5d0 .and. xyz(2) .ge. 0.5d0) then
c(1,1)
c
      f(1)=y
      f(2)=-1.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. 2.5d0 .and. xyz(1) .ge. 1.5d0 .and.
     1    xyz(2) .lt. 1.5d0 .and. xyz(2) .ge. 0.5d0) then
c(2,1)
c
      f(1)=y
      f(2)=-6.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. -0.5d0 .and. xyz(1) .ge. -1.5d0 .and.
     1    xyz(2) .lt. 1.5d0 .and. xyz(2) .ge. 0.5d0) then
c(-1,1)
c
      f(1)=y
      f(2)=8.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. -1.5 .and. xyz(1) .ge. -2.5d0 .and.
     1    xyz(2) .lt. 1.5d0 .and. xyz(2) .ge. 0.5d0) then
c（-2,1)
c
      f(1)=y
      f(2)=13.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. 0.5d0 .and. xyz(1) .ge. -0.5d0 .and.
     1    xyz(2) .lt. 2.5d0 .and. xyz(2) .ge. 1.5d0) then
c(0,2)
c
      f(1)=-1.0d0+y
      f(2)=3.5d0+5.0d0*x

      elseif(xyz(1) .lt. 1.5d0 .and. xyz(1) .ge. 0.5d0 .and.
     1    xyz(2) .lt. 2.5d0 .and. xyz(2) .ge. 1.5d0) then
c(1,2)
c
      f(1)=-1.0d0+y
      f(2)=-1.5d0+5.0d0*x2
c
      elseif(xyz(1) .lt. 2.5d0 .and. xyz(1) .ge. 1.5d0 .and.
     1    xyz(2) .lt. 2.5d0 .and. xyz(2) .ge. 1.5d0) then
c(2,2)
c
      f(1)=-1.0d0+y
      f(2)=-6.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. -0.5d0 .and. xyz(1) .ge. -1.5d0 .and.
     1    xyz(2) .lt. 2.5d0 .and. xyz(2) .ge. 1.5d0) then
c(-1,2)
c
      f(1)=-1.0d0+y
      f(2)=8.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. -1.5 .and. xyz(1) .ge. -2.5d0 .and.
     1    xyz(2) .lt. 2.5d0 .and. xyz(2) .ge. 1.5d0) then
c（-2,2)
c
      f(1)=-1.0d0+y
      f(2)=13.5d0+5.0d0*x

c
      elseif(xyz(1) .lt. 0.5d0 .and. xyz(1) .ge. -0.5d0 .and.
     1    xyz(2) .lt. -0.5d0 .and. xyz(2) .ge. -1.5d0) then
c(0,-1)
c
      f(1)=2.0d0+y
      f(2)=3.5d0+5.0d0*x

      elseif(xyz(1) .lt. 1.5d0 .and. xyz(1) .ge. 0.5d0 .and.
     1    xyz(2) .lt. -0.5d0 .and. xyz(2) .ge. -1.5d0) then
c(1,-1)
c
      f(1)=2.0d0+y
      f(2)=-1.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. 2.5d0 .and. xyz(1) .ge. 1.5d0 .and.
     1    xyz(2) .lt. -0.5d0 .and. xyz(2) .ge. -1.5d0) then
c(2,-1)
c
      f(1)=2.0d0+y
      f(2)=-6.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. -0.5d0 .and. xyz(1) .ge. -1.5d0 .and.
     1    xyz(2) .lt. -0.5d0 .and. xyz(2) .ge. -1.5d0) then
c(-1,-1)
c
      f(1)=2.0d0+y
      f(2)=8.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. -1.5 .and. xyz(1) .ge. -2.5d0 .and.
     1    xyz(2) .lt. -0.5d0 .and. xyz(2) .ge. -1.5d0) then
c（-2,-1)
c
      f(1)=2.0d0+y
      f(2)=13.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. 0.5d0 .and. xyz(1) .ge. -0.5d0 .and.
     1    xyz(2) .lt. -1.5d0 .and. xyz(2) .ge. -2.5d0) then
c(0,-2)
c
      f(1)=3.0d0+y
      f(2)=3.5d0+5.0d0*x

      elseif(xyz(1) .lt. 1.5d0 .and. xyz(1) .ge. 0.5d0 .and.
     1    xyz(2) .lt. -1.5d0 .and. xyz(2) .ge. -2.5d0) then
c(1,-2)
c
      f(1)=3.0d0+y
      f(2)=-1.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. 2.5d0 .and. xyz(1) .ge. 1.5d0 .and.
     1    xyz(2) .lt. -1.5d0 .and. xyz(2) .ge. -2.5d0) then
c(2,-2)
c
      f(1)=3.0d0+y
      f(2)=-6.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. -0.5d0 .and. xyz(1) .ge. -1.5d0 .and.
     1    xyz(2) .lt. -1.5d0 .and. xyz(2) .ge. -2.5d0) then
c(-1,-2)
c
      f(1)=3.0d0+y
      f(2)=8.5d0+5.0d0*x
c
      elseif(xyz(1) .lt. -1.5 .and. xyz(1) .ge. -2.5d0 .and.
     1    xyz(2) .lt. -1.5d0 .and. xyz(2) .ge. -2.5d0) then
c（-2,-2)
c
      f(1)=3.0d0+y
      f(2)=13.5d0+5.0d0*x
       
      else
c
      f(1)=0.0d0
      f(2)=0.0d0

      endif


      
      end subroutine






c------------------------------------------
c     the forcing term
c     don't forget to pass visc in dpars to it
      subroutine fevaln(nd,t,xy,u,dpars,zpars,ipars,fu)
      implicit real*8 (a-h,o-z)
      integer nd, ndim, ipars(*)
      parameter(ndim=2)
      real*8 xy(ndim), t, u(nd), fu(nd)
      real*8 dpars(*), visc(nd)
      complex*16 zpars(*)
      real *8 acoeff, bcoeff
c
      x=xy(1)
      y=xy(2)
      visc(1)=dpars(1)
      visc(2)=dpars(2)
c
      pi=3.1415926535897932d0
c
      acoeff=1.0d0
      bcoeff=3.4d0

      fu(1)= aceoff+u(1)*u(2)**2 - (bcoeff+1.0d0)*u(1)

      fu(2)=bcoeff*u(1)-u(1)*u(2)**2
c
c     simple ODE test
c
c      fu=-10.0d0*u
      end subroutine




c---------------------------------------
c
c     exact solution:
c     u(x,t)=sin(2*pi*x)*cos(2*pi*y)*exp(-t)
c     v(x,t)=sin(2*pi*x)*cos(2*pi*y)exp(-t)
c---------------------------------------
c     should have the same calling seq
c     as rhsfun in test_boxfgt*.f
      subroutine uexact(nd,xyz,t,dpars,zpars,
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
      u(1)=dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)*dexp(-t)
      u(2)=dsin(4.0d0*pi*x)*dcos(2.0d0*pi*y)*dexp(-2.0d0*t)
c
      end subroutine





c------------------------------------------
c     jacobian of the forcing term
      subroutine devaln(nd,t,xy,u,dpars,zpars,ipars,du)
      implicit real*8 (a-h,o-z)
      integer nd, ndim, ipars(*)
      parameter(ndim=2)
      real*8 xy(ndim), t, u(nd), du(nd,nd)
      real*8 dpars(*), visc(nd)
      complex*16 zpars(*)
      real *8 acoeff,bcoeff
c
      x=xy(1)
      y=xy(2)
      visc(1)=dpars(1)
      visc(2)=dpars(2)
c
      acoeff=1.0d0
      bcoeff=3.4d0
c
      pi=3.1415926535897932d0
ccc      fu=dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)*dexp(-t)
ccc      fu=(-1.0d0+8.0d0*pi**2)*fu
ccc      fu=(-1.0d0+8.0d0*pi**2)*u
      du(1,1)= u(2)**2-(bcoeff+1.0d0)
      du(1,2)= 2.0d0*u(1)*u(2)
      du(2,1)= bcoeff-u(2)**2
      du(2,2)= -2.0d0*u(1)*u(2)
c
c     simple ODE test
c
c      fu=-10.0d0*u
      end subroutine




