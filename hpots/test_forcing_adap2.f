c     a testing driver 
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
      parameter(ndim=2)
      parameter(norder=8)
      parameter(mxboxes=100000)
      parameter(mxltree=200000)
      parameter(mxlevels=30)
      parameter(ntarg=1000000)
      parameter(nx = 1000)
      parameter(ny = 1000)
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 usol(norder**ndim,mxboxes)
      real*8 usolp(ntarg)
      real*8 uextp(ntarg)
      real*8 rintl(0:200), dpars(1000)
      real *8, allocatable :: xref(:,:),wts(:)
      real*8, allocatable:: targ(:)
      real*8, allocatable:: targs(:)
      real*8 hx,hy
      complex*16 zpars(10), zk
      real *8 daprs(10)
      
      real*8 errpe
      character*1 type
      external finit, fforce, uexact
c
      open(21,file='exacttarget18.txt')
      open(234,file='exactgird18.txt')

      pi=3.1415926535897932d0
c
      dt=1.0d-2/(2**14)


      ntot=2**14

c
      write(*,*)'*********',ntot
      tf=ntot*dt
c
      eps=1.0d-9
      npbox=norder**ndim
c
      nordert=2
c      nordert=4
      ifnewtree=1
c
      allocate(targs(ndim))


c     now call the solver
c      call heat2dperfm2(norder, nordert,ntarg, nx,ny,finit,
c     1      fforce,dt, ntot, eps, mxltree, mxboxes, mxlevels,
c     2     ltree, itree, iptr, centers, nlevels,
c     3     boxsize, nboxes, ifnewtree, usol,usolp)
      call heat2dperfm3(norder, nordert,finit,
     1      fforce,dt, ntot, eps, mxltree, mxboxes, mxlevels,
     2     ltree, itree, iptr, centers, nlevels,
     3     boxsize, nboxes, ifnewtree, usol,ntarg,nx,ny,usolp)

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
cc     compute exact solutions on arbitrary targets
c
      do i = 1,ntarg
         uextp(i)=0.0d0
      enddo
c
      ra = 0
      erra = 0

      hx=1.0d0/(nx-1)
      hy=1.0d0/(ny-1)

      do i=1,nx
      do j=1,ny
        ii=(i-1)*ny+j
        targs(1)=-0.5d0+(i-1)*hx
        targs(2)=-0.5d0+(j-1)*hy
        call uexact(tf,targs,dpars,zpars,ipars,uu)
        ra = ra + uu**2
        erra = erra + (uu-usolp(ii))**2
        write(21,*)ii,targs(1),uu,usolp(ii),abs(uu-usolp(ii))
      enddo
      enddo

         erra = sqrt(erra/ra)
 
      write(*,*)'relative pottarg l2 error=*',erra

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
      f(1)=dsin(2.0d0*pi*x+1.0d0)*dcos(2.0d0*pi*y)
      
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
      u(1)= dexp(t)*dcos(2.0d0*pi*x+dexp(t))*dcos(2.0d0*pi*y)
     1   +8.0d0*pi**2*dsin(2.0d0*pi*x+dexp(t))*dcos(2.0d0*pi*y)

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
      u=dsin(2.0d0*pi*x+dexp(t))*dcos(2.0d0*pi*y)
c
      end subroutine



c----------------------------------



