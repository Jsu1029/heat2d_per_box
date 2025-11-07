c test3
c test the subroutine evalt
c subroutine defines the array that represents the
c     right hand side of the equation.

c
      program test
      implicit real *8 (a-h,o-z)
      integer ipars(100), ltree
      integer iptr(8), icoll,ipoly
c     what is iptr? add an unpack subroutine?
      real*8 done, pi, rintl(0:200)
      real*8 dpars(1000)
      complex*16 zpars(10)
      complex*16 zk, ztmp, ima, zz
      real *8 tt
      parameter(nd=2)
      parameter(norder=8)
      parameter(ndim=2)
      integer ntarg
      integer nx,ny
      real*8, allocatable:: targs(:,:)
      real *8 hx,hy
      real *8 targ(ndim)
      real *8 uu(nd)
c
      integer, allocatable :: itree(:)
      integer, allocatable :: isrc(:), isrcse(:,:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8, allocatable :: src(:,:)
      real *8, allocatable:: pote(:,:)

      integer nlevels
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      external uexact
  
      open(21,file='test3.data')

c     set parameters
      call prini(6,13)
      done = 1
      pi = atan(done)*4
c
      eps=1.0d-9
      ipoly=1
      iperiod=1
c
      ifnewtree=0
c     when ifnewtree=1, prepare to refine for
c     one more level? I guess so.
c
      zk=1.0d0
      eta=1.0d0
c     not sure about the above two parameters
c
      epstree=1.0d-1*eps
      boxlen=1.0d0
      iptype=2
      npbox= norder**ndim
      ntarg=1000000
      nx=1000
      ny=1000

      tt=0.0d0
      dpars(1)=0.0d0
c
      allocate(targs(ndim,ntarg))
      hx=1.0d0/(nx-1)
      hy=1.0d0/(ny-1)

      do i=1,nx
      do j=1,ny
        ii=(i-1)*ny+j
        targs(1,ii)=-0.5d0+(i-1)*hx
        targs(2,ii)=-0.5d0+(j-1)*hy
      enddo
      enddo
c
c     subroutine vol_tree_mem returns:
c     nboxes, nlevels, ltree, rintl(0:200)
      call cpu_time(t1)
      call vol_tree_mem(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,uexact,nd,dpars,zpars,ipars,ifnewtree,
     2    nboxes,nlevels,ltree,rintl)
      write(*,*) 'nboxes =', nboxes
      write(*,*) 'nlevels=', nlevels
      write(*,*) 'ltree=', ltree
c     still not sure what is rintl, but okay
      call cpu_time(t2)
c
c     2. allocate memory and then call vol_tree_build
c     (itree, iptr) defines the logic structure
c     while (centers, boxsize) defines the geometry
c     and fvals defines the function values on the vol grid
      allocate(fvals(nd,npbox,nboxes),centers(ndim,nboxes))
      allocate(boxsize(0:nlevels),itree(ltree))
c
c     Is level restriction enforced? Yes.
      call vol_tree_build(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,uexact,nd,dpars,zpars,ipars,rintl,
     2    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals)

c

      allocate(pote(nd,ntarg))

      call evalt(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,fvals,
     2       ntarg,targs,pote)

c
      sum_err=0.0d0
      sum_u=0.0d0
      do i=1,ntarg
         do j =1,ndim
            targ(j)=targs(j,i)
         enddo
         call uexact(nd,targ,dpars,zpars,
     1           ipars,uu)
         write(21,*)uu(1),pote(1,i),uu(2),pote(2,i)
         sum_u=sum_u+uu(1)**2
         sum_err=sum_err+abs(pote(1,i)-uu(1))**2
      enddo
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'err2_rel=',err2_rel


      end program




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
      t=dpars(1)
c
      u(1)=dsin(2.0d0*pi*x+dexp(t))*dcos(2.0d0*pi*y)
      u(2)=1.0d0
c
      end subroutine

