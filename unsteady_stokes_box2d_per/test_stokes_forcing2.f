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
c      open(21,file='usoltarget.txt')
c      open(221,file='usoltarget2.data')

      pi=3.1415926535897932d0
      visc=1.0d0
  
c
      dt=1.0d0/100
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

c
c      do k=3,7

c      dt=1.0d0/(2**k)
c      ntot=2**k
c      write(*,*)'eps =',eps
      call unsteady_stokes2dperfm4_video(nd,norder,nordert,finit,
     1      fforce,dt, ntot, eps, mxltree, mxboxes, mxlevels,
     2     ltree, itree, iptr, centers, nlevels,
     3     boxsize, nboxes, ifnewtree,velocity,vorticity,
     4     ntarg,velocityp,visc,funs)

c    velocity at target points
c      do i=1,ntarg
c         write(21,*)velocityp(1,i)
c         write(221,*)velocityp(2,i)
c      enddo

c      enddo
c    exact solution
      ipoly=1
      dpars(1)=tf
      write(*,*)'final time =', dpars(1)
      call evalf(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,
     2       dpars,zpars,ipars,uexact,velocity_ex)

c
      sum_err=0.0d0
      sum_u=0.0d0
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               sum_u=sum_u+velocity_ex(1,j,ibox)**2
               sum_err=sum_err+abs(velocity(1,j,ibox)-
     1                   velocity_ex(1,j,ibox))**2
             enddo
          endif
        enddo
      enddo
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'usol err2_rel=',err2_rel




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
      u(1)=pi*dexp(dsin(2.0d0*pi*x+dexp(t)))*dexp(dsin(2.0d0*pi*y
     1       +dexp(2.0d0*t)))*dcos(2.0d0*pi*y+dexp(2.0d0*t))
      u(2)=-pi*dexp(dsin(2.0d0*pi*x+dexp(t)))*dexp(dsin(2.0d0*pi*y
     1         +dexp(2.0d0*t)))*dcos(2.0d0*pi*x+dexp(t))
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

c
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
      f(1)=pi*dexp(dsin(2.0d0*pi*x+1.0d0))*dexp(dsin(2.0d0*pi*y
     1       +1.0d0))*dcos(2.0d0*pi*y+1.0d0)
      f(2)=-pi*dexp(dsin(2.0d0*pi*x+1.0d0))*dexp(dsin(2.0d0*pi*y
     1         +1.0d0))*dcos(2.0d0*pi*x+1.0d0)
      
      end subroutine






c------------------------------------------
c     the forcing term
      subroutine fforce(nd,xyz,dpars,zpars,ipars,u)
      implicit real*8 (a-h,o-z)
      integer nd, ipars(*)
      real *8 dpars(*), xyz(*), t, u(nd)
      complex *16 zpars(*)
c
      pi = 3.1415926535897932d0
      x = xyz(1)
      y = xyz(2)
      t = dpars(1)
c
      theta = 2.0d0 * pi * x + dexp(t)
      phi = 2.0d0 * pi * y + dexp(2.0d0*t)
c
      cos_theta = dcos(theta)
      sin_theta = dsin(theta)
      cos_phi = dcos(phi)
      sin_phi = dsin(phi)
c
      A = dexp(dsin(theta))
      B = dexp(dsin(phi))
c
      ! Compute u(1) and u(2) for convenience
      u1 = pi * A * B * cos_phi
      u2 = -pi * A * B * cos_theta
c
      ! Calculate f(1) = u1_t - laplace(u1)
      term1 = u1 * cos_theta * dexp(t)
      term2 = 2.0d0 * dexp(2.0d0*t) * pi * A * B *
     1        (cos_phi**2 - sin_phi)
      term3 = (2.0d0*pi)**2 * u1 *
     1        (cos_theta**2 - sin_theta +
     2         cos_phi**2 - 3.0d0*sin_phi - 1.0d0)
      u(1) = term1 + term2 - term3
c
      ! Calculate f(2) = u2_t - laplace(u2)
c
      u2_t = -pi * A * B * (
     1        (cos_theta**2 - sin_theta) * exp(t) +
     2        cos_theta * cos_phi * 2.0d0 * exp(2.0d0*t))
      term23 = u2 * (2.0d0*pi)**2 *
     1             (cos_theta**2 - 3.0d0*sin_theta - 1.0d0 +
     2              cos_phi**2 - sin_phi)
      u(2) = u2_t - term23
      end subroutine






c---------------------------------------
c     should have the same calling seq
c     as rhsfun in test_boxfgt*.f
      subroutine funs(nd,xyz,dpars,zpars,
     1           ipars,f)
c     the source term f = u_t - laplace(u)
      implicit real *8 (a-h,o-z)
      integer nd, ipars(*)
      real *8 dpars(*),f(nd),xyz(*)
      complex *16 zpars(*)
c
      pi = 3.1415926535897932d0
      x = xyz(1)
      y = xyz(2)
      t = dpars(1)
c
      theta = 2.0d0 * pi * x + dexp(t)
      phi = 2.0d0 * pi * y + dexp(2.0d0*t)
c
      cos_theta = dcos(theta)
      sin_theta = dsin(theta)
      cos_phi = dcos(phi)
      sin_phi = dsin(phi)
c
      A = dexp(dsin(theta))
      B = dexp(dsin(phi))
c
      ! Compute u1_y (partial derivative of u1 with respect to y)
      u1_y = 2.0d0 * pi * pi * A * B *
     1       (cos_phi**2 - sin_phi )
c
      ! Compute u2_x (partial derivative of u2 with respect to x)
      u2_x = -2.0d0 * pi * pi * A * B *
     1       (cos_theta**2 - sin_theta )
c
      ! Compute vorticity = u2_x - u1_y
      f(1) = u2_x - u1_y
c

      f(2) = 0.0
c
      end subroutine
