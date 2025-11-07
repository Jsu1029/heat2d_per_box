c test2
c test the subroutine evalf
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
c
      parameter(mxboxes=100000)
      parameter(mxltree=200000)
      parameter(mxlevels=30)
      integer, allocatable :: itree(:)
      integer, allocatable :: isrc(:), isrcse(:,:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8, allocatable :: src(:,:)
      real *8, allocatable :: fright(:,:,:)

      integer nlevels
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      external fforce
      external uexact
  
      open(21,file='test2.data')

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

      tt=0.0d0
      dpars(1)=0.0d0

c
c     subroutine vol_tree_mem returns:
c     nboxes, nlevels, ltree, rintl(0:200)
      call cpu_time(t1)
      call vol_tree_mem(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,fforce,nd,dpars,zpars,ipars,ifnewtree,
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
     1    norder,iptype,eta,fforce,nd,dpars,zpars,ipars,rintl,
     2    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals)

      allocate(fright(nd,npbox,nboxes))
c
      call evalf(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,
     2       dpars,zpars,ipars,fforce,fright)


c
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               write(21,*)fright(1,j,ibox),fright(2,j,ibox)
             enddo
          endif
        enddo
      enddo

      end program




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
      u(2)=10.0d0

      end subroutine

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
