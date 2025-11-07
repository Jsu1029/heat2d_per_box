c test1
c
c test the subroutine newtree2oldtree2d
c plot old tree and new tree for comparison
c
c     tests pts_tree_sort
c     is the map from point to box needed?
c
      program test
      implicit real *8 (a-h,o-z)
      integer ipars(100), ltree
      integer iptr(8), icoll
c     what is iptr? add an unpack subroutine?
      real*8 done, pi, rintl(0:200)
      real*8 dpars(1000)
      complex*16 zpars(10)
      complex*16 zk, ztmp, ima, zz
c
      parameter(mxboxes=100000)
      parameter(mxltree=200000)
      parameter(mxlevels=30)
      integer, allocatable :: itree(:)
      integer, allocatable :: isrc(:), isrcse(:,:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8, allocatable :: src(:,:)
c

      integer   levelbox(mxboxes)
      integer   nlev, nboxes
      integer   icolbox(mxboxes), irowbox(mxboxes)
      integer   iparentbox(mxboxes), ichildbox(4,mxboxes)
      integer   nblevel(0:mxlevels), iboxlev(mxboxes)
      integer   istartlev(0:mxlevels)
      real *8   cent0(2),xsize0

      integer nlevels
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      external rhsfun
c     set parameters
      call prini(6,13)
      done = 1
      pi = atan(done)*4
c
      eps=1.0d-9
      ipoly=1
      iperiod=1
      norder=8
      nd=1
      ndim=2
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
c
c     subroutine vol_tree_mem returns:
c     nboxes, nlevels, ltree, rintl(0:200)
      call cpu_time(t1)
      call vol_tree_mem(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,rhsfun,nd,dpars,zpars,ipars,ifnewtree,
     2    nboxes,nlevels,ltree,rintl)
      write(*,*) 'nboxes =', nboxes
      write(*,*) 'nlevels=', nlevels 
      write(*,*) 'ltree=', ltree
      do i=1, nlevels
        write(*,*) i, rintl(i)
      enddo
c     still not sure what is rintl, but okay
      call cpu_time(t2)
c
c     2. allocate memory and then call vol_tree_build
c     (itree, iptr) defines the logic structure
c     while (centers, boxsize) defines the geometry
c     and fvals defines the function values on the vol grid
      write(*,*) 'ltree=', ltree
      allocate(fvals(nd,npbox,nboxes),centers(ndim,nboxes))
      allocate(boxsize(0:nlevels),itree(ltree))
c
c     Is level restriction enforced? Yes.
      call vol_tree_build(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,rhsfun,nd,dpars,zpars,ipars,rintl,
     2    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals)

      nlfbox = 0
      do il=1,nlevels
        do ibox=itree(iladdr+2*il),itree(iladdr+2*il+1)
          if(itree(iptr(4)+ibox-1).eq.0) then
            nlfbox = nlfbox+1
          endif
        enddo
      enddo
      call prinf('nlfbox=*',nlfbox,1)

c     3. call print_tree_matlab for visualization
c     make sure you can use it

      ifplot=1
      if (ifplot.eq.1) then
         fname1 = 'trgtree.data'
         fname2 = 'src.data'
         fname3 = 'targ.data'
      
         nt=0
         ns=0
c        given ns and nt
c        generate random sources and targets 
c        and sort into the tree?
         call print_tree_matlab(ndim,itree,ltree,nboxes,centers,
     1       boxsize,nlevels,iptr,ns,src,nt,targs,fname1,fname2,fname3)
      endif
      write(*,*) 'tree printed out for visualization'


c  oldtree to newtree
      call newtree2oldtree2d(nlev,levelbox,iparentbox,
     1    ichildbox,icolbox,irowbox,nboxes,nblevel,
     2    iboxlev,istartlev,cent0,xsize0,iperiod,
     3    ltree,nlevels,itree,iptr,centers,boxsize)

c
       iprint=101
       call printtree(levelbox,icolbox,irowbox,nboxes,
     1       nlev,iparentbox,ichildbox,nblevel,iboxlev,
     2       istartlev,iprint,nleaves,cent0,xsize0)
      write(*,*) 'nleaves=',nleaves


      end program






      subroutine rhsfun(nd,xyz,dpars,zpars,ipars,f)
c     right-hand-side function
c     f has dimension nd
      implicit real *8 (a-h,o-z)
      integer nd,ndim,ipars(*)
      complex *16 zpars(*)
      real *8 dpars(*),f(nd),xyz(*)
      real *8 w(5),c1(5),c2(5)
c
      pi=3.1415926535897932d0
c      f(1)=dsin(4.0d0*pi*xyz(1))*dsin(2.0d0*pi*xyz(2))
c      f(1)=dsin(4.0d0*pi*xyz(1))
      w(1)=100.0d0
      w(2)=200.0d0
      w(3)=300.0d0
      w(4)=400.0d0
      w(5)=500.0d0

      c1(1)=-0.3d0
      c1(2)=-0.188d0
      c1(3)=0.178d0
      c1(4)=-0.088d0
      c1(5)=-0.38d0

      c2(1)=-0.4d0
      c2(2)=0.0d0
      c2(3)=-0.1d0
      c2(4)=0.3d0
      c2(5)=-0.05d0

      res=0.0d0
      do i=1,5
        res=res+dexp(-((xyz(1)-c1(i))**2
     1                +(xyz(2)-c2(i))**2)*w(i))
      enddo
      f(1)=res
c     slightly weird

      end subroutine


c**********************************************************
c
      subroutine posbox(xll,yll,dx,dy,icol,irow,level,
     1           cent0,xsize0)
      implicit real*8 (a-h,o-z)
      integer icol, irow, level
      real*8 cent0(2), xsize0
      
      xlength = xsize0 / dble(2**level)
      x0 = cent0(1)-xsize0/2.0d0
      y0 = cent0(2)-xsize0/2.0d0
c
      dx=xlength
      dy=xlength

      xll=x0+(icol-1)*xlength
      yll=y0+(irow-1)*xlength

      end subroutine


