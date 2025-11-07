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
      parameter(ntarg=100)
      parameter(nx=10)
      parameter(ny=10)
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 usol(nd,norder**ndim,mxboxes)
      real*8 usolp(nd,ntarg)
c      real*8 usol(2,64,mxboxes)
c      real*8 uext(nd,norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000), visc(2)
      real*8 uu(2)
      real *8, allocatable :: xref(:,:),wts(:)
      real*8, allocatable:: targ(:)
      complex*16 zpars(10), zk
      integer itarg
      character*1 type
      external finit, fevaln, devaln, uexact
c
c      dt=1.0d-5
c      dt=1.0d-3
      dt=1.0d-2/2**2
      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c
      ntot=1
      ntot=2**2
      ntot=ntot*2
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
c      ntot=3
      tf=ntot*dt
c
      eps=1.0d-9
      npbox=norder**ndim
c
      nordert=4
c      nordert=4
c     don't change the tree within fgt for now
      ifnewtree=0
c
      iadap = 1
      visc(1) = 2.0d0
      visc(2) = 1.0d0

      itarg=0
c
      call heat2d_sys_am24(nd, norder, nordert, finit,
     1     fevaln, devaln, visc, dt, ntot, eps, mxltree,
     2     mxboxes, mxlevels, ltree, itree, iptr, 
     3     centers, nlevels, boxsize, nboxes, iadap,
     4     ntarg, nx, ny, usol,itarg, usolp)
c
      write(*,*) 'done calling heat2d_sys_am2'
      write(*,*) 'nboxes, nlevels=', nboxes, nlevels
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
      write(*,*) ''
      write(*,*) 'dim 1:'
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
c               dpars(1)=tf
c               call uexact(nd,targ,dpars,zpars,ipars,uu)
               call uexact(nd, targ, tf, dpars, zpars,
     1              ipars, uu)
               sum_u=sum_u+uu(1)**2
               sum_err=sum_err+abs(usol(1,j,ibox)-uu(1))**2
c
               write(234,*) targ(1),targ(2), uu(1), usol(1,j,ibox),
     1                     abs(usol(1,j,ibox)-uu(1))
             enddo
          endif
        enddo
      enddo
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'err2_rel=',err2_rel

      write(*,*) ''
      write(*,*) 'dim 2:'
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
c               dpars(1)=tf
c               call uexact(nd,targ,dpars,zpars,ipars,uu)
               call uexact(nd, targ, tf, dpars, zpars,
     1              ipars, uu)
               sum_u=sum_u+uu(2)**2
               sum_err=sum_err+abs(usol(2,j,ibox)-uu(2))**2
c
               write(235,*) targ(1),targ(2), uu(2), usol(2,j,ibox),
     1                     abs(usol(2,j,ibox)-uu(2))
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
      f(1)=dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)
      f(2)=dsin(4.0d0*pi*x)*dcos(2.0d0*pi*y)
      
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
c
      x=xy(1)
      y=xy(2)
      visc(1)=dpars(1)
      visc(2)=dpars(2)
c
      pi=3.1415926535897932d0
ccc      fu=dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)*dexp(-t)
ccc      fu=(-1.0d0+8.0d0*pi**2)*fu
ccc      fu=(-1.0d0+8.0d0*pi**2)*u
      uex=dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)*dexp(-t)
      fu(1)= u(1)*u(1) + (-1.0d0+visc(1)*8.0d0*pi**2)*uex - uex*uex
c
      vex=dsin(4.0d0*pi*x)*dcos(2.0d0*pi*y)*dexp(-2.0d0*t)
      fu(2)= u(2)*u(2) + (-2.0d0+visc(2)*20.0d0*pi**2)*vex - vex*vex
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
c
      x=xy(1)
      y=xy(2)
      visc(1)=dpars(1)
      visc(2)=dpars(2)
c
      pi=3.1415926535897932d0
ccc      fu=dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)*dexp(-t)
ccc      fu=(-1.0d0+8.0d0*pi**2)*fu
ccc      fu=(-1.0d0+8.0d0*pi**2)*u
      du(1,1)= 2.0d0*u(1) 
      du(1,2)= 0.0d0
      du(2,1)= 0.0d0
      du(2,2)= 2.0d0*u(2)
c
c     simple ODE test
c
c      fu=-10.0d0*u
      end subroutine




