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
c      parameter(mxboxes=100000)
c      parameter(mxltree=200000)
      parameter(mxboxes=50000)
      parameter(mxltree=100000)
      parameter(mxlevels=30)
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 usol(nd,norder**ndim,mxboxes)
c      real*8 usol(2,64,mxboxes)
c      real*8 uext(nd,norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000), visc(2)
      real*8 uu(2)
      real *8, allocatable :: xref(:,:),wts(:)
      real*8, allocatable:: targ(:)
      complex*16 zpars(10), zk
      character*1 type
      external finit, fevaln, devaln, uexact
c
c      dt=1.0d-5
c      dt=1.0d-3
      dt=500.0d0/2000
c
      ntot=2000
c      ntot=ntot*2
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
c
      nordert=2
c      nordert=4
c     don't change the tree within fgt for now
      ifnewtree=0
c
      iadap = 1
      visc(1) = 8.0d-5
      visc(2) = 4.0d-5
c     target point

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
      if(xyz(1) .lt. 0.125d0 .and. xyz(1) .gt. -0.125d0 .and.
     1    xyz(2) .lt. 0.125d0 .and. xyz(2) .gt. -0.125d0) then
c(0,0)
         f(2)=1.0d0/4.0d0*(dsin(8.0d0*pi*
     1        x))**2*(dsin(8.0d0*pi*y))**2
c
      elseif(xyz(1) .lt. 1.125d0 .and. xyz(1) .gt. 0.875d0 .and.
     1    xyz(2) .lt. 0.125d0 .and. xyz(2) .gt. -0.125d0)then
c(1,0)
         f(2)=1.0d0/4.0d0*(dsin(8.0d0*pi*
     1        x))**2*(dsin(8.0d0*pi*y))**2
c
      elseif(xyz(1) .lt. 2.125d0 .and. xyz(1) .gt. 1.875d0 .and.
     1    xyz(2) .lt. 0.125d0 .and. xyz(2) .gt. -0.125d0) then
c(2,0)
         f(2)=1.0d0/4.0d0*(dsin(8.0d0*pi*
     1        x))**2*(dsin(8.0d0*pi*y))**2
c
      elseif(xyz(1) .lt. -0.875d0 .and. xyz(1) .gt. -1.125d0 .and.
     1    xyz(2) .lt. 0.125d0 .and. xyz(2) .gt. -0.125d0) then
c(-1,0)
         f(2)=1.0d0/4.0d0*(dsin(8.0d0*pi*
     1        x))**2*(dsin(8.0d0*pi*y))**2
c
      elseif(xyz(1) .lt. -1.875 .and. xyz(1) .gt. -2.125d0 .and.
     1    xyz(2) .lt. 0.125d0 .and. xyz(2) .gt. -0.125d0) then
c(-2,0)
         f(2)=1.0d0/4.0d0*(dsin(8.0d0*pi*
     1        x))**2*(dsin(8.0d0*pi*y))**2
c
      elseif(xyz(1) .lt. 0.125d0 .and. xyz(1) .gt. -0.125d0 .and.
     1    xyz(2) .lt. 1.125d0 .and. xyz(2) .gt. 0.875d0) then
c(0,1)
         f(2)=1.0d0/4.0d0*(dsin(8.0d0*pi*
     1        x))**2*(dsin(8.0d0*pi*y))**2

      elseif(xyz(1) .lt. 1.125d0 .and. xyz(1) .gt. 0.875d0 .and.
     1    xyz(2) .lt. 1.125d0 .and. xyz(2) .gt. 0.875d0)then
c(1,1)
         f(2)=1.0d0/4.0d0*(dsin(8.0d0*pi*
     1        x))**2*(dsin(8.0d0*pi*y))**2
c
      elseif(xyz(1) .lt. 2.125d0 .and. xyz(1) .gt. 1.875d0 .and.
     1    xyz(2) .lt. 1.125d0 .and. xyz(2) .gt. 0.875d0) then
c(2,1)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2
c
      elseif(xyz(1) .lt. -0.875d0 .and. xyz(1) .gt. -1.125d0 .and.
     1    xyz(2) .lt. 1.125d0 .and. xyz(2) .gt. 0.875d0) then
c(-1,1)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2
c
      elseif(xyz(1) .lt. -1.875 .and. xyz(1) .gt. -2.125d0 .and.
     1    xyz(2) .lt. 1.125d0 .and. xyz(2) .gt. 0.875d0) then
c（-2,1)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2

c
      elseif(xyz(1) .lt. 0.125d0 .and. xyz(1) .gt. -0.125d0 .and.
     1    xyz(2) .lt. 2.125d0 .and. xyz(2) .gt. 1.875d0) then
c(0,2)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2

      elseif(xyz(1) .lt. 1.125d0 .and. xyz(1) .gt. 0.875d0 .and.
     1    xyz(2) .lt. 2.125d0 .and. xyz(2) .gt. 1.875d0)then
c(1,2)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2
c
      elseif(xyz(1) .lt. 2.125d0 .and. xyz(1) .gt. 1.875d0 .and.
     1    xyz(2) .lt. 2.125d0 .and. xyz(2) .gt. 1.875d0) then
c(2,2)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2
c
      elseif(xyz(1) .lt. -0.875d0 .and. xyz(1) .gt. -1.125d0 .and.
     1    xyz(2) .lt. 2.125d0 .and. xyz(2) .gt. 1.875d0) then
c(-1,2)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2
c
      elseif(xyz(1) .lt. -1.875 .and. xyz(1) .gt. -2.125d0 .and.
     1    xyz(2) .lt. 2.125d0 .and. xyz(2) .gt. 1.875d0) then
c（-2,2)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2

c
      elseif(xyz(1) .lt. 0.125d0 .and. xyz(1) .gt. -0.125d0 .and.
     1    xyz(2) .lt. -0.875d0 .and. xyz(2) .gt. -1.125d0) then
c(0,-1)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2

      elseif(xyz(1) .lt. 1.125d0 .and. xyz(1) .gt. 0.875d0 .and.
     1    xyz(2) .lt. -0.875d0 .and. xyz(2) .gt. -1.125d0)then
c(1,-1)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2
c
      elseif(xyz(1) .lt. 2.125d0 .and. xyz(1) .gt. 1.875d0 .and.
     1    xyz(2) .lt. -0.875d0 .and. xyz(2) .gt. -1.125d0) then
c(2,-1)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2
c
      elseif(xyz(1) .lt. -0.875d0 .and. xyz(1) .gt. -1.125d0 .and.
     1    xyz(2) .lt. -0.875d0 .and. xyz(2) .gt. -1.125d0) then
c(-1,-1)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2
c
      elseif(xyz(1) .lt. -1.875 .and. xyz(1) .gt. -2.125d0 .and.
     1    xyz(2) .lt. -0.875d0 .and. xyz(2) .gt. -1.125d0) then
c（-2,-1)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2
c
      elseif(xyz(1) .lt. 0.125d0 .and. xyz(1) .gt. -0.125d0 .and.
     1    xyz(2) .lt. -1.875d0 .and. xyz(2) .gt. -2.125d0) then
c(0,-2)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2

      elseif(xyz(1) .lt. 1.125d0 .and. xyz(1) .gt. 0.875d0 .and.
     1    xyz(2) .lt. -1.875d0 .and. xyz(2) .gt. -2.125d0)then
c(1,-2)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2
c
      elseif(xyz(1) .lt. 2.125d0 .and. xyz(1) .gt. 1.875d0 .and.
     1    xyz(2) .lt. -1.875d0 .and. xyz(2) .gt. -2.125d0) then
c(2,-2)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2
c
      elseif(xyz(1) .lt. -0.875d0 .and. xyz(1) .gt. -1.125d0 .and.
     1    xyz(2) .lt. -1.875d0 .and. xyz(2) .gt. -2.125d0) then
c(-1,-2)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2
c
      elseif(xyz(1) .lt. -1.875 .and. xyz(1) .gt. -2.125d0 .and.
     1    xyz(2) .lt. -1.875d0 .and. xyz(2) .gt. -2.125d0) then
c（-2,-2)
         f(2)=1.0d0/4.0d0*dsin(8.0d0*pi*
     1        x)**2*dsin(8.0d0*pi*y)**2


      else
         f(2)=0.0d0
      endif

c
      f(1)=1.0d0-2.0d0*f(2)

c      f(1)=dsin(4.0d0*pi*x)*dsin(4.0d0*pi*y)
c      f(2)=f(1)

      
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
      real *8 gamma, kappa
c
      x=xy(1)
      y=xy(2)
      visc(1)=dpars(1)
      visc(2)=dpars(2)
c
      pi=3.1415926535897932d0
c
      gamma=0.024d0
      kappa=0.06d0

      fu(1)= -u(1)*u(2)**2 + gamma*(1-u(1))

      fu(2)=u(1)*u(2)**2 - (gamma+kappa)*u(2)
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
      real *8 gamma, kappa
c
      x=xy(1)
      y=xy(2)
      visc(1)=dpars(1)
      visc(2)=dpars(2)
c
      gamma=0.024d0
      kappa=0.06d0
c
      pi=3.1415926535897932d0
ccc      fu=dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)*dexp(-t)
ccc      fu=(-1.0d0+8.0d0*pi**2)*fu
ccc      fu=(-1.0d0+8.0d0*pi**2)*u
      du(1,1)= -u(2)**2-gamma
      du(1,2)= -2.0d0*u(1)*u(2)
      du(2,1)= u(2)**2
      du(2,2)= 2.0d0*u(1)*u(2)-(kappa+gamma)
c
c     simple ODE test
c
c      fu=-10.0d0*u
      end subroutine




