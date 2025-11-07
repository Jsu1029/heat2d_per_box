c     a testing driver tests3
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
      parameter(mxboxes=1000000)
      parameter(mxltree=2000000)
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
c
      open(21,file='step64.txt')
      open(221,file='usoltarget2.data')

      pi=3.1415926535897932d0
      visc=1.0d0


    
      allocate(velocity_ex(nd,norder**ndim,mxboxes))
      allocate(velocity(nd,norder**ndim,mxboxes))
      allocate(vorticity(norder**ndim,mxboxes))
c
c      write(*,*)'time steps=',ntot
c      tf=ntot*dt
c
      eps=1.0d-9
      npbox=norder**ndim
c
      ifnewtree=0

      nordert=4

c
c      do k=6,6

      dt=1.0d-1/100
      ntot=100
      write(*,*)'eps =',eps
      call unsteady_stokes2dperfm4(nd,norder,nordert,finit,
     1      fforce,dt, ntot, eps, mxltree, mxboxes, mxlevels,
     2     ltree, itree, iptr, centers, nlevels,
     3     boxsize, nboxes, ifnewtree,velocity,vorticity,
     4     ntarg,velocityp,visc,funs)

c    velocity at target points
      do i=1,ntarg
         write(21,*)velocityp(1,i)
         write(221,*)velocityp(2,i)
      enddo

c      enddo





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
      f(1)=0.0d0  ! u component
      f(2)=0.0d0  ! v component
      do i = 0,4
            ! Calculate the distance components
            dx = x - (c1(2) + i - 2)
            dy = y - (c1(1) + i - 2)
            
            ! Calculate the Gaussian term
            gaussian = dexp(-(dx**2 + dy**2)/del)
            
            ! Calculate the derivatives for velocity components
            ! u = -∂ψ/∂y = -∑[2(y - c1(1) - i + 2)/del * exp(-(...)/del)]
            ! v = ∂ψ/∂x = ∑[2(x - c1(2) - i + 2)/del * exp(-(...)/del)]
            f(1) = f(1) - (2.0d0*dy/del) * gaussian  ! u component
            f(2) = f(2) + (2.0d0*dx/del) * gaussian  ! v component
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
      real *8  c1(2),c2(2),dx,dy,guassian
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
      f(1)=0.0d0  ! u component
      f(2)=0.0d0  ! v component
      do i = 0,4
            ! Calculate the distance components
            dx = x - (c1(2) + i - 2)
            dy = y - (c1(1) + i - 2)
            
            ! Calculate the Gaussian term
            gaussian = dexp(-(dx**2 + dy**2)/del)
            
            ! Calculate the derivatives for velocity components
            ! u = -∂ψ/∂y = -∑[2(y - c1(1) - i + 2)/del * exp(-(...)/del)]
            ! v = ∂ψ/∂x = ∑[2(x - c1(2) - i + 2)/del * exp(-(...)/del)]
            f(1) = f(1) - (2.0d0*dy/del) * gaussian  ! u component
            f(2) = f(2) + (2.0d0*dx/del) * gaussian  ! v component
      enddo
c
      end subroutine






c------------------------------------------
c     the forcing term
      subroutine funs(nd,xyz,dpars,zpars,ipars,u)
      implicit real*8 (a-h,o-z)
      integer nd, ipars(*)
      real *8 dpars(*), xyz(*), t, u(nd)
      complex *16 zpars(*)
      real *8  c1(2),c2(2),c1_t(2)
c
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
      t=dpars(1)
c
      c1(1)=0.25d0*dcos(20.0d0*pi*t)
      c1(2)=0.25d0*dsin(20.0d0*pi*t)
      
      c1_t(1) = -0.25d0 * 20.0d0 * pi * dsin(20.0d0*pi*t)  ! dc1(1)/dt
      c1_t(2) =  0.25d0 * 20.0d0 * pi * dcos(20.0d0*pi*t)  ! dc1(2)/dt
c
      del = 1.0d-2
c
      u(1)=0.0d0
      u(2)=0.0d0
      do i = 0,4
            dx = x - (c1(2) + i - 2)
            dy = y - (c1(1) + i - 2)
            r2 = dx**2 + dy**2
            g = dexp(-r2/del)
            
            ! Calculate u and v components
            u_comp = - (2.0d0*dy/del) * g
            v_comp =   (2.0d0*dx/del) * g
            
            ! Calculate time derivatives u_t and v_t
            u_t = -(2.0d0/del) * g * (
     1           dy * (2.0d0/del) * (dx*c1_t(2) + dy*c1_t(1))
     2           - c1_t(1))
            
            v_t = -(2.0d0/del) * g * (
     1           -dx * (2.0d0/del) * (dx*c1_t(2) + dy*c1_t(1))
     2           + c1_t(2))
c
            ! Calculate Laplace u and Laplace v
            laplace_u = (4.0d0/(del**2)) * g * (
     1           dy * (1.0d0-2.0d0/del*dx**2)
     2           - 2.0d0*dy**3/(del) + 3.0d0*dy)
            
            laplace_v = (4.0d0/(del**2)) * g * (
     1           -dx * (1.0d0-2.0d0/del*dy**2)
     2           +2.0d0*dx**3/(del) - 3.0d0*dx)
            
            ! Sum up contributions
            u(1) = u(1) + u_t - laplace_u
            u(2) = u(2) + v_t - laplace_v
      enddo

      end subroutine






c---------------------------------------
c     should have the same calling seq
c     as rhsfun in test_boxfgt*.f
      subroutine fforce(nd,xyz,dpars,zpars,
     1           ipars,f)
c     the initial condition
c     defined to be a single Fourier mode
      implicit real *8 (a-h,o-z)
      integer nd, ipars(*)
      real *8 dpars(*),f(nd),xyz(*)
      complex *16 zpars(*)
      real *8  c1(2), c1_t(2)
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
      
      c1_t(1) = -0.25d0 * 20.0d0 * pi * dsin(20.0d0*pi*t)
      c1_t(2) =  0.25d0 * 20.0d0 * pi * dcos(20.0d0*pi*t)

c
      f(1)=0.0d0
      f(2)=0.0d0
      do i = 0,4
            dx = x - (c1(2) + i - 2)
            dy = y - (c1(1) + i - 2)
            r2 = dx**2 + dy**2
            g = dexp(-r2/del)
            
            ! Calculate u and v components
            u_comp = - (2.0d0*dy/del) * g
            v_comp =   (2.0d0*dx/del) * g
            
            ! Calculate time derivatives u_t and v_t
            u_t = -(2.0d0/del) * g * (
     1           dy * (2.0d0/del) * (dx*c1_t(2) + dy*c1_t(1))
     2           - c1_t(1))
            
            v_t = -(2.0d0/del) * g * (
     1           -dx * (2.0d0/del) * (dx*c1_t(2) + dy*c1_t(1))
     2           + c1_t(2))
c
            ! Calculate Laplace u and Laplace v
            laplace_u = (4.0d0/(del**2)) * g * (
     1           dy * (1.0d0-2.0d0/del*dx**2)
     2           - 2.0d0*dy**3/(del) + 3.0d0*dy)
            
            laplace_v = (4.0d0/(del**2)) * g * (
     1           -dx * (1.0d0-2.0d0/del*dy**2)
     2           +2.0d0*dx**3/(del) - 3.0d0*dx)
            
            ! Sum up contributions
            f(1) = f(1) + u_t - laplace_u
            f(2) = f(2) + v_t - laplace_v
      enddo

      end subroutine
