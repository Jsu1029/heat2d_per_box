c     a testing driver tests1
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
      open(21,file='usoltarget.txt')
      open(234,file='exactgird18.txt')

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
      write(*,*)'eps =',eps
      call unsteady_stokes2dperfm4(nd,norder,nordert,finit,
     1      fforce,dt, ntot, eps, mxltree, mxboxes, mxlevels,
     2     ltree, itree, iptr, centers, nlevels,
     3     boxsize, nboxes, ifnewtree,velocity,vorticity,
     4     ntarg,velocityp,visc,funs)

c    velocity at target points
      do i=1,ntarg
         write(21,*)i,velocityp(1,i),velocityp(2,i)
      enddo

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
     1      8.0d0*pi**2*dsin(2.0d0*pi*x+dexp(t))*dcos(2.0d0*pi*y)
     2     +2.0d0*pi*dcos(2.0d0*pi*x)*dcos(2.0d0*pi*y)*dexp(t)

c
       u(2) = dexp(t)*dsin(2.0d0*pi*x+dexp(t))*dsin(2.0d0*pi*y)
     1   -8.0d0*pi**2*dcos(2.0d0*pi*x+dexp(t))*dsin(2.0d0*pi*y)
     2   -2.0d0*pi*dsin(2.0d0*pi*x)*dsin(2.0d0*pi*y)*dexp(t)

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

c       f(1)=2.0d0*pi*dcos(2.0d0*pi*x)*dcos(2.0d0*pi*y)*dexp(t)
c       f(2)=-2.0d0*pi*dsin(2.0d0*pi*x)*dsin(2.0d0*pi*y)*dexp(t)

      end subroutine
    
