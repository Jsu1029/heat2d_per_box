c     a testing driver tests4
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
      external finit, fforce, uexact,funs


      pi=3.1415926535897932d0
      visc=1.0d0
  
c
      dt=1.0d-1/100

      ntot=100


    
      allocate(velocity_ex(nd,norder**ndim,mxboxes))
      allocate(velocity(nd,norder**ndim,mxboxes))
      allocate(vorticity(norder**ndim,mxboxes))
c
      write(*,*)'time steps=',ntot
      tf=ntot*dt
c
      eps=1.0d-9
      npbox=norder**ndim
c
      ifnewtree=0

      nordert=4

      write(*,*)'eps =',eps
      call unsteady_stokes2dperfm4_video(nd,norder,nordert,finit,
     1      fforce,dt, ntot, eps, mxltree, mxboxes, mxlevels,
     2     ltree, itree, iptr, centers, nlevels,
     3     boxsize, nboxes, ifnewtree,velocity,vorticity,
     4     ntarg,velocityp,visc,funs)



   





      end program



c---------------------------------------
c     should have the same calling seq
c     as rhsfun in test_boxfgt*.f
      subroutine uexact(nd,xyz,dpars,zpars,
     1           ipars,f)
      implicit real *8 (a-h,o-z)
      integer nd, ipars(*)
      real *8 dpars(*), xyz(*), t, f(nd)
      complex *16 zpars(*)
      real *8  c1(2),c2(2)
c
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
      t=dpars(1)
c
      del = 1.0d-2
c
      c1(1)=0.25d0*dcos(20.0d0*pi*t)
      c1(2)=0.25d0*dsin(20.0d0*pi*t)

c
      f(1)=0.0d0
      f(2)=0.0d0
      do i = 0,4
            u1=(y-c1(1)-i+2)**2
            u2=(x-c1(2)-i+2)**2
            f(1)=f(1)+dexp(-u1/del)+dexp(-u2/del)
            f(2)=f(2)+(2.0d0*(x-c1(2)-i+2)/del)*
     1           dexp(-u2/del)*(y-i+2)
      enddo
c
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
      real *8  c1(2),c2(2)
c
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
      t=0.0d0
c
      del = 1.0d-2
c
      c1(1)=0.25d0*dcos(20.0d0*pi*t)
      c1(2)=0.25d0*dsin(20.0d0*pi*t)

c
      f(1)=0.0d0
      f(2)=0.0d0
      do i = 0,4
            u1=(y-c1(1)-i+2)**2
            u2=(x-c1(2)-i+2)**2
            f(1)=f(1)+dexp(-u1/del)+dexp(-u2/del)
            f(2)=f(2)+(2.0d0*(x-c1(2)-i+2)/del)*
     1           dexp(-u2/del)*(y-i+2)
      enddo
      
      end subroutine






c------------------------------------------
c     the forcing term
      subroutine funs(nd,xyz,dpars,zpars,ipars,u)
      implicit real*8 (a-h,o-z)
      integer nd, ipars(*)
      real *8 dpars(*), xyz(*), t, u(nd)
      complex *16 zpars(*)
      real *8  c1(2),c2(2)
c
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
      t=dpars(1)
c
      c1(1)=0.25d0*dcos(20.0d0*pi*t)
      c1(2)=0.25d0*dsin(20.0d0*pi*t)
c
      del = 1.0d-2
c
      u(1)=0.0d0
      u(2)=0.0d0
      do i = 0,4
            u1=(y-c1(1)-i+2)**2
            u2=(x-c1(2)-i+2)**2
            usqrt1=(y-c1(1)-i+2)
            usqrt2=(x-c1(2)-i+2)
            u(1)=u(1)+
     1         dexp(-u1/del)*(-2.0d0/del*usqrt1*0.25d0*
     2          dsin(20.0d0*pi*t)*20.0d0*pi)
     3       -(4.0d0/(del**2)*u1)*dexp(-u1/del)+(2.0d0/del)*
     4        dexp(-u1/del)+dexp(-u2/del)*(2.0d0/del*usqrt2*0.25d0*
     2         dcos(20.0d0*pi*t)*20.0d0*pi)
     2        -(4.0d0/(del**2)*u2)*dexp(-u2/del)+(2.0d0/del)*
     3        dexp(-u2/del)

           u(2)=u(2)-(2.0d0/del)*0.25d0*
     1         dcos(20.0d0*pi*t)*20.0d0*pi*dexp(-u2/del)*(y-i+2)+
     2         (4.0d0/(del**2)*u2)*dexp(-u2/del)*0.25d0*
     3         dcos(20.0d0*pi*t)*20.0d0*pi*(y-i+2)+
     4         ((2.0d0/del-4.0d0/(del**2)*u2)*2.0d0/del*usqrt2)*
     5         dexp(-u2/del)*(y-i+2)+(8.0d0/(del**2)*usqrt2)*
     6         dexp(-u2/del)*(y-i+2)
      enddo

      end subroutine






c     as rhsfun in test_boxfgt*.f
      subroutine fforce(nd,xyz,dpars,zpars,
     1           ipars,f)
c     the initial condition
c     defined to be a single Fourier mode
      implicit real *8 (a-h,o-z)
      integer nd, ipars(*)
      real *8 dpars(*),f(nd),xyz(*)
      complex *16 zpars(*)
      real *8  c1(2),c2(2)
c
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
      t=dpars(1)
c
      del = 1.0d-2
c
      c1(1)=0.25d0*dcos(20.0d0*pi*t)
      c1(2)=0.25d0*dsin(20.0d0*pi*t)

c
      f(1)=0.0d0
      f(2)=0.0d0
      do i = 0,4
            u1=(y-c1(1)-i+2)**2
            u2=(x-c1(2)-i+2)**2
            usqrt1=(y-c1(1)-i+2)
            usqrt2=(x-c1(2)-i+2)
            f(1)=f(1)+
     1         dexp(-u1/del)*(-2.0d0/del*usqrt1*0.25d0*
     2          dsin(20.0d0*pi*t)*20.0d0*pi)
     3       -(4.0d0/(del**2)*u1)*dexp(-u1/del)+(2.0d0/del)*
     4        dexp(-u1/del)+dexp(-u2/del)*(2.0d0/del*usqrt2*0.25d0*
     2         dcos(20.0d0*pi*t)*20.0d0*pi)
     2        -(4.0d0/(del**2)*u2)*dexp(-u2/del)+(2.0d0/del)*
     3        dexp(-u2/del)

            f(2)=f(2)-(2.0d0/del)*0.25d0*
     1         dcos(20.0d0*pi*t)*20.0d0*pi*dexp(-u2/del)*(y-i+2)+
     2         (4.0d0/(del**2)*u2)*dexp(-u2/del)*0.25d0*
     3         dcos(20.0d0*pi*t)*20.0d0*pi*(y-i+2)+
     4         ((2.0d0/del-4.0d0/(del**2)*u2)*2.0d0/del*usqrt2)*
     5         dexp(-u2/del)*(y-i+2)+(8.0d0/(del**2)*usqrt2)*
     6         dexp(-u2/del)*(y-i+2)
      enddo

      end subroutine


