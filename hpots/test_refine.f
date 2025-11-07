c     a testing driver
c     which tests subroutines:
c     vol_tree_find_box_refine and 
c     vol_tree_refine_boxes
c
c     by calling vol_tree_build to resolve
c     a simple function rhsfun1
c     and calling subsequent subroutines
c     to refine more towards function rhsfun2
c
      program test
      implicit real *8 (a-h,o-z)
      integer ipars(100), ltree
      integer iptr(12), icoll
c     what is iptr? add an unpack subroutine?
      parameter(mxboxes=100000)
      parameter(mxltree=200000)
      parameter(mxlevels=30)
      real*8 done, pi, rintl(0:200)
      real*8 dpars(1000)
      complex*16 zpars(10)
      complex*16 zk, ztmp, ima, zz
      integer, allocatable :: itree(:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8, allocatable :: fval1(:,:,:,:),centerstmp(:,:,:)
      integer, allocatable :: irefinebox(:),ifrefine(:)
c     what are these things...
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8, allocatable :: rmask(:)
c
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      external rhsfun1, rhsfun2
c
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
c     one more level?
c
      zk=1.0d0
      eta=0.0d0
c     not sure about the above two parameters
c
      epstree=1.0d-1*eps
      boxlen=1.0d0
      iptype=2
      npbox= norder**ndim
      npols = norder**ndim
c
      mc=2**ndim
      mnbors=3**ndim
c
c     subroutine vol_tree_mem returns:
c     nboxes, nlevels, ltree, rintl(0:200)
      call cpu_time(t1)
      call vol_tree_mem(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,rhsfun1,nd,dpars,zpars,ipars,ifnewtree,
     2    nboxes,nlevels,ltree,rintl)
      write(*,*) 'nboxes =', nboxes
      write(*,*) 'nlevels=', nlevels 
      write(*,*) 'ltree=', ltree 
      do i=0, nlevels
        write(*,*) i, rintl(i)
      enddo
c     still not sure what is rintl, but okay
      call cpu_time(t2)
c
c     2. allocate memory and then call vol_tree_build
c     (itree, iptr) defines the logic structure
c     while (centers, boxsize) defines the geometry
c     and fvals defines the function values on the vol grid
c
c     attention: need to increase the sizes 
c     since the tree has to be adjusted later
      write(*,*) 'ltree=', ltree
      allocate(fvals(nd,npbox,mxboxes),centers(ndim,mxboxes))
      allocate(boxsize(0:mxlevels),itree(mxltree))
c
c     Is level restriction enforced? Yes.
      call vol_tree_build(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,rhsfun1,nd,dpars,zpars,ipars,rintl,
     2    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals)
      call prinf('laddr=*',itree,2*(nlevels+1))
c     laddr is within itree, starting from 2*(nlevels+1),
c     I suppose?
c     not familiar with laddr, is this iboxlev?
c
      nlfbox = 0
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nlfbox = nlfbox+1
      enddo
      call prinf('nlfbox=*',nlfbox,1)
c
c     3. call print_tree_matlab for visualization
c     make sure you can use it
c
      ifplot=0
      if (ifplot.eq.1) then
         fname1 = 'trgtree.data'
         fname2 = 'src.data'
         fname3 = 'targ.data'
      
         ns=0
         nt=0
c        given ns and nt
c        generate random sources and targets 
c        and sort into the tree?
         call print_tree_matlab(ndim,itree,ltree,nboxes,centers,
     1       boxsize,nlevels,iptr,ns,src,nt,targs,fname1,fname2,fname3)
         write(*,*) 'tree printed out for visualization'
      endif
c
c------------------------------------------------------------
      write(*,*) '-----------'
      write(*,*) 'the pointers:'
      write(*,*) iptr(1)
      write(*,*) iptr(2), itree(iptr(2))
      write(*,*) iptr(3)
      write(*,*) iptr(4)
      write(*,*) iptr(5)
      write(*,*) iptr(6)
      write(*,*) iptr(7)
      write(*,*) iptr(8)
      write(*,*) '-----------'
c     that's the problem
c     not enough for refinement
c
c------------------------------------------------------------
c
c     now a tree (with only one node) is created 
c     call refinement related routines to resolve
c     the other function rhsfun2
cc
c     let's first do a uniform refinement of several levels
c     shall we?
c
c     allocate: rmask
c     grid, xq, wts, umat, vmat
c
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
      allocate(grid(ndim,npbox))
c
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)
      do i=1,norder
        xq(i) = xq(i)/2
      enddo

c      if(1 .eq. 2) then
      if (ndim.eq.1) then
         do i=1,norder
            grid(1,i)=xq(i)
         enddo
      elseif (ndim.eq.2) then
         call mesh2d(xq,norder,xq,norder,grid)
      elseif (ndim.eq.3) then
         call mesh3d(xq,norder,xq,norder,xq,norder,grid)
      endif
c     what is xq... 
c     grid is the grid on the root box
c
c     print out grid
c      do j=1, npbox
c        write(*,*) j, grid(1,j), grid(2,j)
c      enddo 
c
      allocate(rmask(npbox))
      call tens_prod_get_rmask(ndim,iptype,norder,npbox,
     1    rmask,rsum)
c
cccccc     
c      do j=1, npbox
c        write(*,*) j, rmask(j)
c      enddo
c
c
c      if(1 .eq. 2) then
c------------
      do iref = 1, 1
        write(*,*) 'iref=', iref
        nbctr = nboxes
c       get all the leaf boxes, in this situation: 
c       the finest level
        ilev = nlevels
c       the finest level
        ifirstbox = itree(2*ilev+1) 
        ilastbox = itree(2*ilev+2)
        write(*,*) 'ilev, ifirstbox, ilastbox=', ilev,
     1              ifirstbox, ilastbox

        nbloc = ilastbox-ifirstbox+1
        write(*,*) 'nbloc=', nbloc
c
        allocate(fval1(nd,npols,mc,nbloc))
        allocate(centerstmp(ndim,mc,nbloc))
        allocate(irefinebox(nbloc))
c
        rsc = rsc*rintl(ilev)
c       whatever this means
c
        do ii=1, nbloc
          irefinebox(ii)=1
        enddo
c       set flag irefinebox: refine them all
        irefine = 1
c       what does this mean?
c       how to determine if a new level is created? 
c-----------------------------------
c
        if(irefine .eq. 1) then
c         set a few parameters outside
          boxsize(ilev+1) = boxsize(ilev)/2
c          itree(2*ilev+3) = nbctr+1
c         what is 2*ilev+3?
c
c        call vol_tree_refine_boxes
c        to do the main work
c         needs a grid
          write(*,*) 'before vol_tree_refine_boxes'
          write(*,*) 2*ilev+3, iptr(2), itree(iptr(2))
          call vol_tree_refine_boxes(ndim,irefinebox,nd,npbox,fvals,
     1      rhsfun1,dpars,zpars,ipars,grid,nboxes,ifirstbox,nbloc,
     2      centers,boxsize(ilev+1),nbctr,ilev+1,itree(iptr(2)),
     3      itree(iptr(3)),itree(iptr(4)),itree(iptr(5)))
c        is the geometry related stuff modified inside?
          
c          itree(2*ilev+4) = nbctr
c         what is 2*ilev+4?
          write(*,*) 'nbctr=', nbctr
c         nboxes need to be adjusted outside?
          
        else
          exit
        endif
c
        deallocate(fval1,irefinebox,centerstmp)
c
c       complete the tree structure?
        nboxes = nbctr
        nlevels = ilev+1
c
        write(*,*) 'nboxes, nlevels=', nboxes, nlevels

        do i=1,nboxes
           itree(iptr(6)+i-1) = 0
           do j=1,mnbors
              itree(iptr(7)+mnbors*(i-1)+j-1) = -1
           enddo
        enddo
c
        call computecoll(ndim,nlevels,nboxes,itree(iptr(1)),
     1      boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2      itree(iptr(5)),iperiod,itree(iptr(6)),itree(iptr(7)))

c        if(nlevels.ge.2) then
c           call vol_tree_fix_lr(ndim,iperiod,fun,nd,dpars,zpars,ipars,
c     1         norder,npbox,fvals,grid,centers,nlevels,nboxes,boxsize,
c     2         nboxes,nlevels,itree(iptr(1)),itree(iptr(2)),
c     3         itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
c     4         itree(iptr(6)),itree(iptr(7)))
c        endif
c
      enddo
cccccc
cccccc
c      endif
      ltree = iptr(8)-1
c     iptr(8) is not modified either      
c
c     as it seems: vol_tree_refine_boxes 
c     doesn't change nboxes, nlevels, ltree???
      write(*,*) 'nboxes =', nboxes
      write(*,*) 'nlevels=', nlevels 
      write(*,*) 'ltree=', ltree 
c
      write(*,*) ' '
      nlfbox = 0
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nlfbox = nlfbox+1
      enddo
      call prinf('nlfbox=*',nlfbox,1)
c
c     not quite right yet
c
c------------
c      endif
c
c------------------------------------------------------------
c     before printing the geometry
c     check the self-consistency of the tree
c
      if(2 .eq. 2) then
      do ilevel=0, nlevels
        do ibox = itree(2*ilevel+1), itree(2*ilevel+2)
          write(*,*) ilevel,itree(iptr(2)+ibox-1),
     1               ibox
          write(*,*) 'parent:', ibox, '->', 
     1               itree(iptr(3)+ibox-1)
           write(*,*) 'num of children of', ibox, 
     1                itree(iptr(4)+ibox-1)
           write(*,*) 'children of', ibox, ':'
           write(*,*) itree(iptr(5)+4*ibox-4)
           write(*,*) itree(iptr(5)+4*ibox-3)
           write(*,*) itree(iptr(5)+4*ibox-2)
           write(*,*) itree(iptr(5)+4*ibox-1)
           write(*,*) 'num of colls of', ibox, 
     1                itree(iptr(6)+ibox-1)
ccc         periodic copies included
ccc         (only for leaf nodes)
ccc         (ordering of children boxes?)
ccc         (ref. vol_tree_refine_boxes)
           write(*,*) itree(iptr(7)+9*ibox-9)
           write(*,*) itree(iptr(7)+9*ibox-8)
           write(*,*) itree(iptr(7)+9*ibox-7)
           write(*,*) itree(iptr(7)+9*ibox-6)
           write(*,*) itree(iptr(7)+9*ibox-5)
           write(*,*) itree(iptr(7)+9*ibox-4)
           write(*,*) itree(iptr(7)+9*ibox-3)
           write(*,*) itree(iptr(7)+9*ibox-2)
           write(*,*) itree(iptr(7)+9*ibox-1)
        enddo
      enddo
      write(*,*) 'iptr(8)=',iptr(8)
      endif
c
c
c------------------------------------------------------------
c
c     call print_tree_matlab for visualization
c     make sure you can use it
c
      ifplot=1
      if (ifplot.eq.1) then
         fname1 = 'trgtree.data'
         fname2 = 'src.data'
         fname3 = 'targ.data'
      
         ns=0
         nt=0
c        given ns and nt
c        generate random sources and targets 
c        and sort into the tree?
         call print_tree_matlab(ndim,itree,ltree,nboxes,centers,
     1       boxsize,nlevels,iptr,ns,src,nt,targs,fname1,fname2,fname3)
      endif
      write(*,*) 'tree printed out for visualization'
c
c------------------------------------------------------------
c
c     what is rsc though? a scaling factor
c
c     what is rintbs?
c     laddr: one field of itree,
c     laddr(2,0:nlmax)
c     ifirstbox, ilastbox of each level
c
c     call vol_tree_find_box_refine
c     
c     rmask computed above
      if(1 .eq. 2) then
c
      do ilev=0, nlevels-1
        irefine = 0

        ifirstbox = itree(2*ilev+1) 
        ilastbox = itree(2*ilev+2)

        nbloc = ilastbox-ifirstbox+1

        allocate(fval1(nd,npols,mc,nbloc))
        allocate(centerstmp(ndim,mc,nbloc))
        allocate(irefinebox(nbloc))

        rsc = rsc*rintl(ilev)
        call vol_tree_find_box_refine(ndim,nd,iptype,eta,eps,zk,
     1      norder,npbox,fvals,npols,umat,rmask,rsum,boxsize(ilev),
     2      nboxes,ifirstbox,nbloc,rsc,irefinebox,irefine)
c 
        if(irefine.eq.1) then
          boxsize(ilev+1) = boxsize(ilev)/2
          itree(2*ilev+3) = nbctr+1

          call vol_tree_refine_boxes(ndim,irefinebox,nd,npbox,fvals,
     1      rhsfun2,dpars,zpars,ipars,grid,nboxes,ifirstbox,nbloc,
     2      centers, boxsize(ilev+1),nbctr,ilev+1,itree(iptr(2)),
     3      itree(iptr(3)),itree(iptr(4)),itree(iptr(5)))
          
          itree(2*ilev+4) = nbctr
        else
          exit
        endif

        deallocate(fval1,irefinebox,centerstmp)
      enddo
c
      nboxes0 = nbctr
      nlevels = ilev

      do i=1,nboxes0
         itree(iptr(6)+i-1) = 0
         do j=1,mnbors
            itree(iptr(7)+mnbors*(i-1)+j-1) = -1
         enddo
      enddo

      call computecoll(ndim,nlevels,nboxes0,itree(iptr(1)),
     1    boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2    itree(iptr(5)),iperiod,itree(iptr(6)),itree(iptr(7)))

      if(nlevels.ge.2) then
         call vol_tree_fix_lr(ndim,iperiod,fun,nd,dpars,zpars,ipars,
     1       norder,npbox,fvals,grid,centers,nlevels,nboxes0,boxsize,
     2       nboxes,nlevels,itree(iptr(1)),itree(iptr(2)),
     3       itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     4       itree(iptr(6)),itree(iptr(7)))
      endif
c
c------------------------------------------------------------
      endif
c


      


      end program







c------------------------------------------------------------
c
      subroutine rhsfun1(nd,xyz,dpars,zpars,ipars,f)
c     right-hand-side function
c     f has dimension nd
      implicit real *8 (a-h,o-z)
      integer nd,ndim,ipars(*)
      complex *16 zpars(*)
      real *8 dpars(*),f(nd),xyz(*)
      real *8 w(5),c1(5),c2(5)
c
      pi=3.1415926535897932d0
      f(1)=1.0d0
c
      end subroutine





c------------------------------------------------------------
c
      subroutine rhsfun2(nd,xyz,dpars,zpars,ipars,f)
c     right-hand-side function
c     f has dimension nd
c
c     one single gaussian 
c     centered at the origin
c     pretty sharp
      implicit real *8 (a-h,o-z)
      integer nd,ndim,ipars(*)
      complex *16 zpars(*)
      real *8 dpars(*),f(nd),xyz(*)
      real *8 w(5),c1(5),c2(5)
c
      pi=3.1415926535897932d0
      w(1)=500.0d0

      c1(1)=0.0d0

      c2(1)=0.0d0

      res=0.0d0
      do i=1,1
        res=res+dexp(-((xyz(1)-c1(i))**2
     1                +(xyz(2)-c2(i))**2)*w(i))
      enddo
c
      f(1)=res

      end subroutine







