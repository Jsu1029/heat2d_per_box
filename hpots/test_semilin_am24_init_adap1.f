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
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 u0(norder**ndim,mxboxes)
      real*8 u1(norder**ndim,mxboxes)
      real*8 un(norder**ndim,mxboxes)
      real*8 usol(norder**ndim,mxboxes)
      real*8 uext(norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
      real *8, allocatable :: xref(:,:),wts(:)
      real*8, allocatable:: targ(:)
      complex*16 zpars(10), zk
      character*1 type
      external finit, feval, uexact
c
c      dt=1.0d-5
c      dt=1.0d-3
      dt=1.0d-1/2**5
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
      ntot=2**5
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
      nordert=2
c      nordert=3
      nordert=4
c     don't change the tree within fgt for now
      ifnewtree=0

      ifexact = 0

      iadap=1
      call heat2d_am24(norder, nordert, finit, feval,
     1     dt, ntot, eps, mxltree, mxboxes, mxlevels,
     2     ltree, itree, iptr, centers, nlevels,
     3     boxsize, nboxes, ifnewtree, ifexact,iadap,
     4     usol)
c
c     need to create the tree and stuff
c     use the tree of the previous solver
c
      allocate(xref(ndim,npbox),wts(npbox))
      allocate(targ(ndim))
      itype = 0
      type='f'
c
      ipoly=1
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c
      sum_err=0.0d0
      sum_u=0.0d0
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               do k=1,ndim
                 targ(k)=centers(k,ibox) + xref(k,j)*bs
               enddo
c               write(*,*) j, ibox, targ(1), targ(2), tf
               dpars(1)=tf
               call uexact(nd,targ,dpars,zpars,ipars,uu)
               sum_u=sum_u+uu**2
               sum_err=sum_err+abs(usol(j,ibox)-uu)**2
c
               write(234,*) targ(1),targ(2), uu, usol(j,ibox),
     1                     abs(usol(j,ibox)-uu)
             enddo
          endif
        enddo
      enddo
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'err2_rel=',err2_rel

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
c     v(x,t)=sin(2*pi*x)*cos(2*pi*y)*exp(-t)
c     v(x,0)=sin(2*pi*x)*cos(2*pi*y)
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
      nd=1
      ndim=2
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
c
      f(1)=dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)
      
      end subroutine






c------------------------------------------
c     the forcing term
      subroutine feval(t,x,y,u,fu)
      implicit real*8 (a-h,o-z)
      real*8 x, y, t, u
c
      pi=3.1415926535897932d0
ccc      fu=dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)*dexp(-t)
ccc      fu=(-1.0d0+8.0d0*pi**2)*fu
ccc      fu=(-1.0d0+8.0d0*pi**2)*u
      uex=dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)*dexp(-t)
      fu= u*u + (-1.0d0+8.0d0*pi**2)*uex - uex*uex
c
c     simple ODE test
c
c      fu=-10.0d0*u
      end subroutine




c---------------------------------------
c
c     exact solution:
c     v(x,t)=sin(2*pi*x)*cos(2*pi*y)*exp(-t)
c---------------------------------------
c     should have the same calling seq
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
      u=dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)*dexp(-t)
c
      end subroutine

