c

      subroutine vol_tree_adapt(ndim,ipoly,iperiod,eps,zk,boxlen,
     1    norder,iptype,eta,fun,nd,dpars,zpars,ipars,rintl,nboxes,
     2    nlevels,ltree,itree,iptr,centers,boxsize,fvals)
c
c      compute the tree
c
c      input parameters:
c        eps - double precision
c           precision requested
c        zk - double complex
c           Helmholtz parameter
c        boxlen - double precision
c           length of box in which volume is contained, 
c           if boxlen = L, then volume is [-L/2,L/2]^2
c        norder - integer
c           order of discretization
c        iptype - integer
c           error norm
c           iptype = 0 - linf
c           iptype = 1 - l1
c           iptype = 2 - l2
c        eta - double precision
c           scaling parameter for error
c        fun - function handle
c           function to evalute it everywhere in the volume
c        nd - integer
c           number of real functions returned by fun
c        dpars - double precision
c           real parameters for function evaluation
c        zpars - double complex
c           complex parameters for function evaluation
c        ipars - integer
c           integer parameters for function evaluation
c        nlevels - integer
c          number of levels
c        nboxes - integer
c          number of boxes
c        ltree - integer
c          length of tree = 2*(nlevels+1)+17*nboxes
c        rintl - real *8 (0:nlevels)
c          estimate of lp norm for scaling the errors
c          at various levels. 
c          We require the estimate at each level to make sure
c          that the memory estimate code is consitent
c          with the build code else there could be potential
c          memory issues 
c         fvals - function values at discretization nodes
c
c      output:
c        itree - integer(ltree)
c          tree info
c        iptr - integer(8)
c          iptr(1) - laddr
c          iptr(2) - ilevel
c          iptr(3) - iparent
c          iptr(4) - nchild
c          iptr(5) - ichild
c          iptr(6) - ncoll
c          iptr(7) - coll
c          iptr(8) - ltree
c        fvals - double precision (nd,norder**ndim,nboxes)
c          function values at new discretization nodes
c        centers - double precision (ndim,nboxes)
c          xyz coordinates of box centers in the oct tree
c        boxsize - double precision (0:nlevels)
c          size of box at each of the levels
c

      implicit none
      real *8 eps,boxlen,eta,dpars(*)
      complex *16 zk,zpars(*)
      integer ndim,ipoly,iperiod,nd,ipars(*),iptype
      integer nlevels,nboxes,norder
      integer iptr(8),ltree
      integer itree(ltree),ier
      real *8 fvals(nd,norder**ndim,nboxes),centers(ndim,nboxes)
      real *8, allocatable :: fval1(:,:,:,:),centerstmp(:,:,:)
      integer, allocatable :: irefinebox(:),ifrefine(:)
      real *8 boxsize(0:nlevels)
c      real *8 xq(norder),wts(norder),umat(norder,norder)
c      real *8 vmat(norder,norder)
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8, allocatable :: rmask(:)
      real *8 rintl(0:nlevels)
      real *8 xy(ndim)

      integer i,ilev,irefine,itype,nbmax,nlmax,npbox,npc,ii
      integer ifirstbox,ilastbox,nbctr,nbloc
      real *8 rsc,rsum

      real *8 ra
      integer j,nboxes0,npols,mc,mnbors,nrefine

      external fun

      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
c
      mc=2**ndim
      mnbors=3**ndim
      
      iptr(1) = 1
      iptr(2) = 2*(nlevels+1)+1
      iptr(3) = iptr(2) + nboxes
      iptr(4) = iptr(3) + nboxes
      iptr(5) = iptr(4) + nboxes
      iptr(6) = iptr(5) + mc*nboxes
      iptr(7) = iptr(6) + nboxes
      iptr(8) = iptr(7) + mnbors*nboxes



      boxsize(0) = boxlen

      do i=1,ndim
         centers(i,1) = 0
      enddo

c
c      set tree info for level 0
c
      itree(1) = 1
      itree(2) = 1
      itree(iptr(2)) = 0
      itree(iptr(3)) = -1
      itree(iptr(4)) = 0
      do i=1,mc
        itree(iptr(5)+i-1) = -1
      enddo
c
c
      npbox = norder**ndim
      allocate(grid(ndim,npbox))
c
c     Generate a grid on the box [-1/2,1/2]^2
c
c
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)
      do i=1,norder
        xq(i) = xq(i)/2
      enddo

      if (ndim.eq.1) then
         do i=1,norder
            grid(1,i)=xq(i)
         enddo
      elseif (ndim.eq.2) then
         call mesh2d(xq,norder,xq,norder,grid)
      elseif (ndim.eq.3) then
         call mesh3d(xq,norder,xq,norder,xq,norder,grid)
      endif
      
      npols = norder**ndim

      do i=1,npbox
         do j=1,ndim
            xy(j) = grid(j,i)*boxlen
         enddo
        call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,1))
      enddo

c     compute rmask and rsum for estimating whether the function is resolved
      allocate(rmask(npbox))

      call tens_prod_get_rmask(ndim,iptype,norder,npbox,
     1    rmask,rsum)
      

c
c       Reset nlevels, nboxes
c
      nbctr = 1

      do ilev=0,nlevels-1
        irefine = 0

        ifirstbox = itree(2*ilev+1) 
        ilastbox = itree(2*ilev+2)

        nbloc = ilastbox-ifirstbox+1

        allocate(fval1(nd,npols,mc,nbloc))
        allocate(centerstmp(ndim,mc,nbloc))
        allocate(irefinebox(nbloc))

        
        if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**ndim)
        if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**ndim)
        if(iptype.eq.0) rsc = 1.0d0
        rsc = rsc*rintl(ilev)
        call vol_tree_find_box_refine(ndim,nd,iptype,eta,eps,zk,
     1      norder,npbox,fvals,npols,umat,rmask,rsum,boxsize(ilev),
     2      nboxes,ifirstbox,nbloc,rsc,irefinebox,irefine)
        

        if(irefine.eq.1) then
          boxsize(ilev+1) = boxsize(ilev)/2
          itree(2*ilev+3) = nbctr+1

          call vol_tree_refine_boxes(ndim,irefinebox,nd,npbox,fvals,
     1      fun,dpars,zpars,ipars,grid,nboxes,ifirstbox,nbloc,centers,
     2      boxsize(ilev+1),nbctr,ilev+1,itree(iptr(2)),itree(iptr(3)),
     3      itree(iptr(4)),itree(iptr(5)))
          
          itree(2*ilev+4) = nbctr
        else
          exit
        endif

        deallocate(fval1,irefinebox,centerstmp)
      enddo

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

      call prinf('nboxes in tree_build=*',nboxes,1)

cccc      call prinf('nboxes0=*',nboxes0,1)
cccc      call prinf('nlevels=*',nlevels,1)
      return
      end
c
