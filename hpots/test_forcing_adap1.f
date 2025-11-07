c     a testing driver 
c     which tests the correctness of 
c     subroutine heat2dperfm
c
c     this version tests the correctness
c     of the solver for an essentially nonadaptive case
c     with a known solution
c
c     see subroutines finit, fforce, uexact 
c     for more details
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
      real*8 usol(norder**ndim,mxboxes)
      real*8 uext(norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
      real *8, allocatable :: xref(:,:),wts(:)
      real*8, allocatable:: targ(:)
      complex*16 zpars(10), zk
      character*1 type
      external finit, fforce, uexact
c
      open(21,file='exacttarget18.txt')
      open(234,file='exactgird18.txt')
c
      dt=1.0d-2/(2**9)
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0
c      dt=dt/2.0d0

      ntot=2**9
c
      tf=ntot*dt
c
      eps=1.0d-9
      npbox=norder**ndim
c
      nordert=2
c      nordert=4
      ifnewtree=1
c
c     now call the solver 
      call heat2dperfm(norder, nordert, finit, fforce,
     1     dt, ntot, eps, mxltree, mxboxes, mxlevels,
     2     ltree, itree, iptr, centers, nlevels, 
     3     boxsize, nboxes, ifnewtree, usol)
c
c
c      call vol_tree_unpack(iptr, iladdr, ilevel, iparent,     
c     1     nchild, ichild, ncoll, icoll, ltree)
c      write(*,*) 'iladdr =', iladdr
c      write(*,*) 'ilevel =', ilevel
c      write(*,*) 'iparent =', iparent
c      write(*,*) 'nchild =', nchild
c      write(*,*) 'ichild =', ichild
c      write(*,*) 'ncoll =', ncoll
c      write(*,*) 'coll =', icoll
c      write(*,*) 'ltree =', ltree+1
c     why is ltree increased by 1? I'd like to fix this.
c      ltree=ltree-1
c
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
c      do i=1, npbox
c        write(122,*) xref(1,i), xref(2,i)
c      enddo
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
               call uexact(tf,targ,dpars,zpars,ipars,uu)
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
      write(*,*) 'solution and error written in fort.234'



      end program





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
      subroutine fforce(nd,xyz,dpars,zpars,ipars,u)
      implicit real*8 (a-h,o-z)
      integer nd, ndim, ipars(*)
      real *8 dpars(*), xyz(*), t, u(nd)
      complex *16 zpars(*)
c
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
      t=dpars(1)
c
      u(1)=dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)*dexp(-t)
      u(1)=(-1.0d0+8.0d0*pi**2)*u(1)

      end subroutine





c---------------------------------------
c     should have the same calling seq
c     as rhsfun in test_boxfgt*.f
      subroutine uexact(t,xyz,dpars,zpars,
     1           ipars,u)
      implicit real *8 (a-h,o-z)
      integer nd, ndim, ipars(*)
      real *8 dpars(*), xyz(*), t, u
      complex *16 zpars(*)
c
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
c
      u=dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)*dexp(-t)
c
      end subroutine

