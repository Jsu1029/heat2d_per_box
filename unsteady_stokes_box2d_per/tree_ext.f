      subroutine vol_tree_unpack(iptr, iladdr, ilevel, iparent,
     1           nchild, ichild, ncoll, coll, ltree)
      implicit none
      integer iptr(8) 
      integer iladdr, ilevel, iparent, nchild, ichild
      integer ncoll, coll, ltree
c
      iladdr = iptr(1)
      ilevel = iptr(2)
      iparent = iptr(3)
      nchild = iptr(4)
      ichild = iptr(5)
      ncoll = iptr(6)
      coll = iptr(7)
      ltree = iptr(8)-1
c     the above entries define the logic structure
c     the geometry is defined by two more arrays
c     centers(ndim, nboxes) and boxlen(0:nlevels)
c
      end subroutine





c----------------------------------------------------------------
c     prints the tree structure to matlab files for visualization
c----------------------------------------------------------------
      subroutine print_tree_matlab(ndim,itree,ltree,nboxes,centers,
     1   boxsize,nlevels,iptr,ns,src,nt,targ,fname1,fname2,fname3)
c
c        this subroutine writes the tree info to a file
c
c        input arguments:
c          itree - integer (ltree)
c             packed array containing tree info
c          ltree - integer
c            length of itree
c          nboxes - integer
c             number of boxes
c          centers - real *8 (2,nboxes)
c             xy coordinates of box centers in tree hierarchy
c          boxsize - real *8 (0:nlevels)
c             size of box at various levels
c          nlevels - integer
c             number of levels
c          iptr - integer(12)
c            pointer to various arrays inside itree
c          ns - integer
c            number of sources
c          src - real *8 (2,ns)
c            xy coorindates of source locations
c          nt - integer
c            number of targets
c          targ - real *8 (2,nt)
c            xy cooridnates of target locations
c          fname1 - character *
c            file name to which tree info is to be written
c          fname1 - character *
c            file name to which source points are to be written
c          fname3 - character *
c            file name to which target points are to be written
c 
c          output
c          files with name fname1, fname2, fname3,
c            which contains the tree info, source points, target points
c            file can be plotted using the matlab script
c            tree_plot.m

      implicit real *8 (a-h,o-z)
      integer itree(ltree),ltree,nboxes,nlevels,iptr(*),ns,nt
      real *8 centers(ndim,nboxes),boxsize(0:nlevels)
      real *8 x(5),y(5),src(ndim,ns),targ(ndim,nt)
      character (len=*) fname1,fname2,fname3

      open(unit=33,file=trim(fname1))
c      open(55,file='targetp.txt')
      nleafbox = 0
      
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
      enddo
      write(*,*) 'nleafbox=',nleafbox

c 1111 format(10(2x,e11.5))      

      do ibox=1,nboxes
         if(itree(iptr(4)+ibox-1).eq.0) then
           ilev = itree(iptr(2)+ibox-1)
           bs = boxsize(ilev)
           x1 = centers(1,ibox) - bs/2
           x2 = centers(1,ibox) + bs/2
c          if a leaf node,
c          get level and box size
c          get xrange [x1,x2]
           if (ndim.eq.2) then
              y1 = centers(2,ibox) - bs/2
              y2 = centers(2,ibox) + bs/2
           endif
c           write(55,*)centers(1,ibox),centers(2,ibox),ibox,bs,ilev,x1,x2
c          if ndim == 2
c          get yrange [x1,x2]

           if (ndim.eq.2) then
              write(33,*) x1,x2,x2,x1,x1,y1,y1,y2,y2,y1
c              write(34,*) x1, y1, x2, y2, ibox, ilev 
cccc         change to the format I'm familiar with?
           else
              write(33,*) x1,x2
           endif
         endif
      enddo
      close(33)
    
 2222 format(2(2x,e11.5))
     
      open(unit=33,file=trim(fname2))
      if (ns .gt. 0) then
         do i=1,ns
            write(33,2222) src(1,i),src(2,i)
         enddo

      endif
      
      close(33)
      open(unit=33,file=trim(fname3))
      if (nt .gt. 0) then
         do i=1,nt
            write(33,2222) targ(1,i),targ(2,i)
         enddo
      endif

      close(33)
     
      return
      end
c
c
c
c
c
c------------------------------------------------
c     given a function handle: fun
c     and an initial tree: (itree, iptr)
c     adapts the tree
c     maintaining arrays of function values: hvals
c
c     Attention: the initial tree should be 
c     big enough to hold newly added entries
c
c     (it should be proceeded by vol_tree_mem)
c     (be a little more generous in vol_tree_mem)
c     (at least add some length in iptr(2))
c
c     Attention: remember to maintain 
c     actual number of boxes and levels 
c     and max number of boxes and levels
c
c     Attention: iperiod added to the calling seq
c
c------------------------------------------------
c
      subroutine vol_tree_adap_fun(nd,ndim,eps,ipoly,
     1    iptype,norder,npbox,nboxes,mboxes,nlevels,
     2    mlevels,ltree,itree,iptr,centers,boxsize,
     3    rintl,dpars,zpars,ipars,iperiod,fun,fvals)
c
c    Input:
c     nd - integer
c            number of right hand sides
c     ndim - integer
c            dimension of the underlying space
c     eps - double precision
c            scaled desired precision
c     ipoly - integer
c            0: Legendre polynomials
c            1: Chebyshev polynomials
c     iptype - integer
c           error norm
c           iptype = 0 - linf
c           iptype = 1 - l1
c           iptype = 2 - l2
c     norder - integer
c           order of expansions for input function value array
c     npbox - integer
c           number of points per box where potential is to be dumped = (norder**ndim)
c     nboxes - integer
c            number of boxes
c
c     input and output:
c     nlevels - integer
c            number of levels
c     ltree - integer
c            length of array containing the tree structure
c     itree - integer(ltree)
c            array containing the tree structure
c     iptr - integer(8)
c            pointer to various parts of the tree structure
c           iptr(1) - laddr
c           iptr(2) - ilevel
c           iptr(3) - iparent
c           iptr(4) - nchild
c           iptr(5) - ichild
c           iptr(6) - ncoll
c           iptr(7) - coll
c           iptr(8) - ltree
c     centers - double precision (ndim,mboxes)
c           xyz coordintes of boxes in the tree structure
c     boxsize - double precision (0:nlevels)
c           size of boxes at each of the levels
c
c     dpars - double precision
c           real parameters for function evaluation
c     zpars - double complex
c           complex parameters for function evaluation
c     ipars - integer
c           integer parameters for function evaluation
c     fun - a function handle, the funciton to be
c           resolved on the tree
c
c----------------
c     Output:
c      (itree, iptr): adjusted
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
c      (centers, boxsize): adjusted
c        centers - double precision (ndim,mboxes)
c          xyz coordinates of box centers in the oct tree
c        boxsize - double precision (0:nlevels)
c          size of box at each of the levels
c
c      fvals: double precision (nd,norder**ndim,mboxes)
c        array of function values, corresponding to fun
c
c     attention:
c     mboxes: here it means the max number of boxes
c     nlevels: here it means the actual number of levels on input
c     
c------------
c     we use nboxes1 and nlevels1 to mean the new numbers
c     of boxes and levels inside this subroutine
c     decide on whether to overwrite nboxes and nlevels later
c------------
c
      implicit real*8 (a-h,o-z) 
      integer nd,ndim,mboxes,nlevels
      integer iptr(8),ltree,ipoly
      integer itree(ltree),norder,npbox
      real *8 eps, eta
      real *8 fvals(nd,npbox,mboxes),fcoefs(nd,npbox,mboxes)
      real *8 centers(ndim,mboxes),boxsize(0:mlevels)
      external fun
      complex *16 zk
c     allocatables
c     grid related stuff
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8, allocatable :: rmask(:), umat_nd(:,:,:)
      real *8 rintl(0:mlevels)
      real *8 xy(ndim)
      real *8, allocatable:: coefsp(:,:,:)
c     flags for refinement and coarsening
      integer isgn(ndim,2**ndim)
c      integer nblock(0:nlevels), ilevstart(0:nlevels)
      integer, allocatable:: nblock(:)
      integer, allocatable:: ilevstart(:)
      integer, allocatable:: nboxid(:), ifdelete(:)
      integer, allocatable:: ifrefine(:)
      integer, allocatable:: irefinelev(:), imaxrefinelev(:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
c
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
      allocate(umat_nd(norder,norder,ndim))
c
      mc=2**ndim
      mnbors=3**ndim
      npols = norder**ndim
c
c      write(*,*) '***'
c      write(*,*) nd, ndim, eps
c      write(*,*) ipoly, iptype, norder
c      write(*,*) npbox, nboxes, mboxes
c      write(*,*) nlevels, mlevels, ltree
c
c      write(*,*) 'centers:'
c      do ibox=1, nboxes
c        write(*,*) centers(1,ibox), centers(2,ibox)
c      enddo
c
c      write(*,*) 'boxsize:'
c      do ilev=0, nlevels
c        write(*,*) ilev, boxsize(ilev)
c      enddo
c     done with sanity check
c      write(*,*) '***'
c      write(*,*) ''
c
      call get_child_box_sign(ndim,isgn)
c      do j=1, 2**ndim
c      do k=1, ndim
c        write(*,*) k,j,isgn(k,j)
c      enddo
c      enddo
c
      allocate(grid(ndim,npbox))
      allocate(grad(nd,ndim,npbox,1))
      allocate(hess(nd,ndim*(ndim+1)/2,npbox,1))
c
c     initialize the grid related stuff
c
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)
c
      norder2=norder*norder
      do i=1,ndim
         call dcopy_f77(norder2,umat,1,umat_nd(1,1,i),1)
      enddo
c
      do i=1,norder
        xq(i) = xq(i)/2
      enddo
c
      if (ndim.eq.1) then
         do i=1,norder
            grid(1,i)=xq(i)
         enddo
      elseif (ndim.eq.2) then
         call mesh2d(xq,norder,xq,norder,grid)
      elseif (ndim.eq.3) then
         call mesh3d(xq,norder,xq,norder,xq,norder,grid)
      endif
c
c      write(*,*) 'grid:'
c      do i=1, ndim
c      do j=1, npbox
c        write(*,*) i,j,grid(i,j)
c        write(201,*) i,j,grid(i,j)
c      enddo
c      enddo
c     grid looks correct
c
c--------------------------------------------------
c     on the currect tree
c     eval fun on the grid of leaf nodes
      do ilev=0, nlevels
        do ibox = itree(2*ilev+1), itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild .eq. 0) then
            do i=1,npbox
              do j=1,ndim
                xy(j) = grid(j,i)*boxsize(ilev)+centers(j,ibox)
c                write(*,*) i, j, xy(j)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,ibox))
c              write(*,*) i, ibox, fvals(1,i,ibox)
            enddo
          endif
        enddo
      enddo
c
c     compute rmask and rsum for estimating whether the function is resolved
      allocate(rmask(npbox))

      call tens_prod_get_rmask(ndim,iptype,norder,npbox,
     1    rmask,rsum)
c
c      write(*,*) 'rsum=',rsum
c
      if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**ndim)
      if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**ndim)
      if(iptype.eq.0) rsc = 1.0d0
c
c     compute the Lp norm of pot??
c
c     rsc = rsc*rnorm
c     problematic?
      reps = eps*rsc
c      write(*,*) 'reps=',reps
c
c     save the original number of boxes
ccc   no need to do this now 
ccc   since nboxes is an input
      nboxes0 = nboxes
      nlevels0 = nlevels
      write(*,*) 'nboxes0=', nboxes0
      write(*,*) 'nboxes=', nboxes
      write(*,*) 'mboxes=', mboxes
c
      write(*,*) '***'
c
c     get rid of the max refine levels and stuff
      allocate(ifrefine(mboxes))
c
c     mboxes here means the max of nboxes
      do ibox = 1, mboxes
        ifrefine(ibox)=0
      enddo
c
c     convert potential to its polynomial expansion coefficients
c     (for all the leaf nodes)
c
      allocate(coefsp(nd,npbox,mboxes))
      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,fvals,coefsp,umat_nd)
c
c      do j=1, npbox
c        write(*,*) coefsp(1,j,5)
c      enddo
c
c     now go down the original tree level by level
c     mark boxes that need to be refined
c
      eta=1.0d0
      zk=1.0d0
c
      do ilev=0, nlevels
        sc = boxsize(ilev)/2
        rscale=sc**eta
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
c         loop over boxes of level ilev
          if (itree(iptr(4)+ibox-1).eq.0) then
c           if a leaf box, check the resolution
             call fun_err(nd,npbox,coefsp(1,1,ibox),
     1            rmask,iptype,rscale,erra)
             erra = erra/rsum
c             write(*,*) ibox, erra, reps
c            need to scale by rsum
             if(erra.gt.reps*rintl(ilev)) then
                ifrefine(ibox)=1
c                write(123,*) ibox, ifrefine(ibox)
             endif
          endif
        enddo
      enddo
c
c     leaf boxes that need to be refined
c     are flaged in ifrefine
c
c-----------------------------------------
c
c      write(*,*) nlevels, mlevels
c
      do ilev = nlevels+1, mlevels
        boxsize(ilev)=boxsize(ilev-1)/2.0d0
c        write(*,*) ilev, boxsize(ilev)
      enddo
c
      allocate(nblock(0:mlevels))
      allocate(ilevstart(0:mlevels+1))
c
      do i=0,mlevels
         nblock(i)=0
      enddo
c
c-----------------------------------------
c
      nnewboxes=0
      do ibox = 1, mboxes
c       change this loop to a level-by-level loop?
c       if no refinement is done on the finest level, exit
c       loop over all boxes, including newly added ones
        write(124,*) 'inspecting box:', ibox
        write(124,*) 'ifrefine=', ifrefine(ibox)
c
        if(ifrefine(ibox) .eq. 1) then
          write(124,*) 'refining', ibox
          ilev=itree(iptr(2)+ibox-1)
          bsh = boxsize(ilev)/4
c         half boxsize on the children's level?
          rscale=bsh**eta
c         set up nchild of ibox
          itree(iptr(4)+ibox-1)=mc
          write(124,*) 'nchild=', itree(iptr(4)+ibox-1)
c         setup for each child
          do j=1, mc
            nnewboxes = nnewboxes+1
            jbox = nboxes0+nnewboxes
            write(124,*) 'jbox=', jbox
c
c           set up centers for the new children
            write(124,*) 'centers of the new boxes:'
            do k=1,ndim
              centers(k,jbox) = centers(k,ibox)+isgn(k,j)*bsh
              write(124,*) k,j,jbox,centers(k,jbox), isgn(k,j), bsh
            enddo
            write(124,*) 'centers of the parent box'
            do k=1, ndim
              write(124,*) k, ibox, centers(k,ibox)
            enddo
c           set up parent for the new children
            itree(iptr(3)+jbox-1)=ibox
            write(124,*) 'parent:'
            write(124,*) jbox, '->', ibox
c
c           set up nchild for the new children
            itree(iptr(4)+jbox-1)=0
c
c           set up children for the new children
            do l=1,mc
              itree(iptr(5)+mc*(jbox-1)+l-1) = -1
            enddo
c
c           set up children for the parent
            itree(iptr(5)+mc*(ibox-1)+j-1)=jbox
            write(124,*) 'children:'
            write(124,*) ibox, '->', jbox
c
c           set up level for the new children
            jlev=ilev+1
            itree(iptr(2)+jbox-1) = jlev
c           count the number of boxes on the level jlev
            nblock(jlev)=nblock(jlev)+1
c            write(124,*) jlev, nblock(jlev)
            write(124,*) 'new box level:'
            write(124,*) jbox, itree(iptr(2)+jbox-1)
c            write(*,*) 'nblock', jlev, nblock(jlev)
c           update nlevels along the way
            if(jlev .gt. nlevels) then
              nlevels = jlev
            endif
c
c           now compute fvals for
c           the newly generated boxes
c           clearly not the right grid
            do ii=1, npbox
              do jj=1, ndim
                xy(jj) = grid(jj,ii)*boxsize(jlev)+centers(jj,jbox)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,ii,jbox))
c              write(*,*) ii,jbox,xy(1),xy(2),fvals(1,ii,jbox)
c              write(*,*) ii,jbox,xy(1),xy(2),boxsize(jlev)
c
c             okay, need to check the error
            enddo
c
c           does it go over all the boxes
c           including newly generated ones?
c
            call ortho_trans_nd(ndim,nd,0,norder,fvals(1,1,jbox),
     1           coefsp(1,1,jbox),umat_nd)

            call fun_err(nd,npbox,coefsp(1,1,jbox),rmask,
     1            iptype,rscale,erra)
            erra = erra/rsum
            write(124,*) 'error on the new box:'
            write(124,*) jbox, erra, reps*rintl(jlev)
c         error obtained for the newly added
c         box jbox
c
            if(erra .gt. reps*rintl(jlev)) then
c             flag the new box
              ifrefine(jbox) = 1
            endif
c         end of the loop over mc children
          enddo
c       end refining box ibox
        endif
      enddo
c
c------------------------------------------------
c
c      write(*,*) 'after refinement, nlevels=',nlevels
c      write(*,*) 'nnewboxes=', nnewboxes
c
c     nblock: number of new boxes on each level
c      write(*,*) 'nblock:'
c      do jj=0, mlevels
c        write(*,*) jj, nblock(jj)
c      enddo
c
c     old id of the boxes??
      allocate(nboxid(mboxes))
      ilevstart(0)=0
      call cumsum(nlevels+1,nblock(0),ilevstart(1))
c      write(*,*) 'after cumsum'
c      write(*,*) 'ilevstart='
c      do jj=0, mlevels
c        write(*,*) jj, ilevstart(jj)
c      enddo
c
      do j=1,nnewboxes
         jbox=nboxes0+j
         jlev=itree(iptr(2)+jbox-1)
         ilevstart(jlev)=ilevstart(jlev)+1
         nboxid(ilevstart(jlev))=jbox
c         write(*,*) j, jlev, jbox, nboxid(ilevstart(jlev))
      enddo
c
c      write(*,*) 'ilevstart='
c      do jj=0, mlevels
c        write(*,*) jj, ilevstart(jj)
c      enddo
c
c     nboxid: the current id of each new box
      write(125,*) 'nboxid='
      do j=1, mboxes
        write(125,*) j, nboxid(j)
      enddo
c
c     nlevels updated along the way, retrieve nnewleverls
c     update nboxes, since nnewboxes is avaiable
c
      nnewlevels = nlevels - nlevels0
      nboxes = nboxes0 + nnewboxes
c
      write(*,*) ' '
      write(*,*) 'mboxes=',mboxes
      write(*,*) 'nboxes=',nboxes
      write(*,*) 'nboxes0=',nboxes0
      write(*,*) 'nnewboxes=',nnewboxes
      write(*,*) ' '
      write(*,*) 'mlevels=',mlevels
      write(*,*) 'nlevels=',nlevels
      write(*,*) 'nlevels0=',nlevels0
      write(*,*) 'nnewlevels=',nnewlevels
      write(*,*) ' '
c
      if (nnewboxes.eq.0) goto 4600
c
      ifpgh = 1
      call vol_tree_reorg_laddr(ndim,mboxes,nd,npbox,
     1    nblock,nboxid,nnewboxes,nboxes0,mlevels,
     2    nnewlevels,nlevels0,centers,itree(iptr(1)),
     3    itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),ifpgh,fvals,coefsp,grad,hess)

c
 4600 continue
c
c
c      write(*,*) 'after reorg, nlevels=', nlevels
c      write(*,*) 'nboxes=', nboxes
c
c     coarsening part, hopefylly done by
c     vol_tree_coarsen and
c     fgt_vol_tree_reorg_after_coarsen 
c
      ifpgh = 1
      allocate(ifdelete(mboxes))
c
c     pot -> fvals
c     grad and hess unavailable yet
c     go over vol_tree_coarsen
      write(*,*) 'calling coarsening routine:'
      call vol_tree_coarsen(nd,ndim,reps,ipoly,norder,npbox,
     1    mboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    ifpgh,fvals,coefsp,grad,hess,
     3    nblock,nboxid,ndelboxes,ifdelete)
c
c      write(*,*) 'after calling coarsening:'
c      write(*,*) 'ndelboxes=', ndelboxes
c
      if (ndelboxes.gt.0) then
         call fgt_vol_tree_reorg_after_coarsen(ndim,mboxes,nd,npbox,
     1       nblock,nboxid,ndelboxes,ifdelete,
     2       centers,nlevels,itree(iptr(1)),itree(iptr(2)),
     3       itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     4       ifpgh,fvals,coefsp,grad,hess)
      endif
c
      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
      nboxes = nboxes - ndelboxes
      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
c------------------------------------------------
c
c     now make it level-restricted again
      call computecoll(ndim,nlevels,nboxes,itree(iptr(1)),
     1    boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2    itree(iptr(5)),iperiod,itree(iptr(6)),itree(iptr(7)))
c
c     call vol_tree_fix_lr_interp
c      if(1 .eq. 2) then
      call vol_tree_fix_lr_interp(ndim,nd,ipoly,iperiod,norder,npbox,
     1    ifpgh,fvals,coefsp,grad,hess,
     2    mboxes,mlevels,centers,boxsize,nboxes,nlevels,
     3    itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),itree(iptr(6)),itree(iptr(7)))
c
      write(*,*) 'after fixing lr:'
      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c      endif


      end subroutine







c-------------------------------------------------
c     subroutine vol_tree_reorg_laddr
c-------------------------------------------------
c
c     re-label boxes of a given tree after refinement
c     so that the box ids are in the ascending order
c
c     this is a modified version of fgt_vol_tree_reorg
c     which is more general, allowing for a trivial
c     tree on input
c
c    INPUT/OUTPUT arguments
c    ndim           in: dimension
c
c    nboxes         in: integer
c                   max number of boxes, could be bigger 
c                   than the actual number of boxes
c
c    nd             in: integer
c                   number of real value functions
c
c    npbox          in: integer
c                   number of grid points per function
c
c    nblock         in: integer(0:nlevels)
c                   number of new boxes on each level
c
c    nboxid:        in: integer(nboxes)
c                   current id of the new boxes
c
c    nnewboxes:     in: integer
c                   number of new boxes 
c    
c    nboxes0:       in: integer
c                   number of old boxes
c
c    centers        in/out: double precision(ndim,nboxes)
c                   x and y coordinates of the center of boxes
c
c    nlevels        in: integer
c                   Number of levels in the tree
c
c    laddr          in/out: integer(2,0:nlevels)
c                   boxes at level i are numbered between
c                   laddr(1,i) to laddr(2,i)
c
c    ilevel      in/out: integer(nboxes)
c                ilevel(i) is the level of box i
c
c    iparent     in/out: integer(nboxes)
c                 iparent(i) is the parent of box i
c
c    nchild      in/out: integer(nboxes)
c                nchild(i) is the number of children 
c                of box i
c
c    ichild       in/out: integer(8,nboxes)
c                 ichild(j,i) is the jth child of box i
c
c ------------------------------------------------
c     RMK1: centers and ilevel - ichild
c                are available for both new and old boxes,
c                while nblock and nboxid are for the new boxes 
c                and laddr is for the old boxes only. 
c
c     RMK2: nlevels is the actual number of levels
c           while nboxes is the max number of boxes
c ------------------------------------------------
c
c   ifpgh:       in: integer
c                if the pot, grad and hess arrays are needed
c   pot, coefsp, grad, hess: arrays of function values,
c                coefficients, gradients and hessians 
c
c
c   On output:
c
c   the tree is reorganized so that boxes are arranged in 
c   ascending order of levels
c   and the arrays of function values, coefficients, gradients
c   and hessians are reorganized accordingly
c      
c-------------------------------------------------
c
      subroutine vol_tree_reorg_laddr(ndim,mboxes,nd,npbox,
     1    nblock,nboxid,nnewboxes,nboxes0,mlevels,nnewlevels,
     2    nlevels0,centers,laddr,ilevel,iparent,
     3    nchild,ichild,ifpgh,pot,coefsp,grad,hess)
      implicit real *8 (a-h,o-z)
c     Calling sequence variables and temporary variables
      integer ndim,mboxes,mlevels,npbox,nd
      double precision centers(ndim,mboxes)
      integer laddr(2,0:mlevels), tladdr(2,0:mlevels)
      integer ilevel(mboxes)
      integer iparent(mboxes)
      integer nchild(mboxes)
      integer ichild(2**ndim,mboxes)
      integer nblock(0:mlevels),nboxid(mboxes)
      double precision pot(nd,npbox,*)
      double precision coefsp(nd,npbox,*)
      double precision grad(nd,ndim,npbox,*)
      double precision hess(nd,ndim*(ndim+1)/2,npbox,*)
      
      integer, allocatable :: tilevel(:),tiparent(:),tnchild(:)
      integer, allocatable :: tichild(:,:)
      integer, allocatable :: iboxtocurbox(:)

      double precision, allocatable :: tcenters(:,:)
      double precision, allocatable :: tpot(:,:,:)
      double precision, allocatable :: tcoefsp(:,:,:)
      double precision, allocatable :: tgrad(:,:,:,:)
      double precision, allocatable :: thess(:,:,:,:)

c     Temporary variables
      integer i,j,k,l,ilevstart(0:mlevels+1)
      integer ibox,ilev, curbox,idim,nblev,mc
c
c      write(*,*) 'inside vol_tree_reorg_laddr'
c
c      t1=second()
      nlevels=nlevels0+nnewlevels
      nboxes=nboxes0+nnewboxes
c
c     retrieve ilevstart
      ilevstart(0)=0
      call cumsum(nlevels+1,nblock(0),ilevstart(1))
      mc=2**ndim
      nhess=ndim*(ndim+1)/2
c
c
c     allocate temp arrays      
      allocate(tilevel(mboxes),tiparent(mboxes),tnchild(mboxes))
      allocate(tichild(mc,mboxes),iboxtocurbox(mboxes))
      allocate(tpot(nd,npbox,mboxes),tcenters(ndim,mboxes))
      allocate(tcoefsp(nd,npbox,mboxes))

      if (ifpgh.ge.2) allocate(tgrad(nd,ndim,npbox,mboxes))
      if (ifpgh.ge.3) allocate(thess(nd,nhess,npbox,mboxes))
c
c      write(*,*) ' '
c      write(*,*) 'within subroutine vol_tree_reorg_laddr:'
c      write(*,*) '... mlevels=',mlevels
c      write(*,*) '... nlevels=',nlevels
c      write(*,*) '... nlevels0=',nlevels0
c
c     copy over laddr for the old boxes
c      write(*,*) 'tladdr='
      do ilev = 0,nlevels0
         tladdr(1,ilev) = laddr(1,ilev)
         tladdr(2,ilev) = laddr(2,ilev)
c         write(*,*) ilev, tladdr(1,ilev), tladdr(2,ilev)
      enddo

c      t2=second()
c      write(*,*) ' '
c
c     copy different fields of the tree:
c     tcenter, tilevel, tiparent, tnchild, tichild, tpot
c     be careful about mboxes v.s. nboxes
c     I believe it should be mboxes here

      call vol_tree_copy(ndim,nd,nboxes,npbox,centers,ilevel,
     1    iparent,nchild,ichild,pot,tcenters,tilevel,
     2    tiparent,tnchild,tichild,tpot)
c      write(*,*) 'done with vol_tree_copy'
c
c     copy coefsp
c     copy grad and hess if asked
      ntot=nd*npbox*nboxes
      call dcopy_f77(ntot,coefsp,1,tcoefsp,1)
      ntotg=ntot*ndim
      if (ifpgh.ge.2) call dcopy_f77(ntotg,grad,1,tgrad,1)
      ntoth=ntot*nhess
      if (ifpgh.ge.3) call dcopy_f77(ntoth,hess,1,thess,1)

c      t3=second()
c      write(*,*) 'done with dcopy_f77'
c      
c     Rearrange old arrays now
c
c     initialize for level 0
      do ilev = 0,0
         do ibox = laddr(1,ilev),laddr(2,ilev)
           iboxtocurbox(ibox) = ibox
c           write(*,*) ibox, ilev
         enddo
      enddo

      np=nd*npbox
      ng=np*ndim
      nh=np*nhess
c
c
c      t4=second()
c     attention: laddr(1,1) unavailable in this case
      if(nboxes .gt. 1) then
c     start from the second box, 
c     the first one besides the root box
      curbox = 2
c      write(*,*) 'curbox=',curbox
c     curbox: the current box, with the correct new id
      do ilev=1,nlevels
         laddr(1,ilev) = curbox
c         write(*,*) ilev, laddr(1,ilev)
c        if ilev <= nlevels0, copy over the old arrays first
         if(ilev .le. nlevels0) then
         do ibox = tladdr(1,ilev),tladdr(2,ilev)
            if((ibox .ge. 1).and.(ibox .le. nboxes0) )then
c            write(*,*) '...',ibox
            ilevel(curbox) = tilevel(ibox)
            nchild(curbox) = tnchild(ibox)
            do k=1,ndim
               centers(k,curbox) = tcenters(k,ibox)
            enddo
            call dcopy_f77(np,tpot(1,1,ibox),1,pot(1,1,curbox),1)
            call dcopy_f77(np,tcoefsp(1,1,ibox),1,coefsp(1,1,curbox),1)
            if (ifpgh.ge.2) call dcopy_f77(ng,tgrad(1,1,1,ibox),1,
     1          grad(1,1,1,curbox),1)
            if (ifpgh.ge.3) call dcopy_f77(nh,thess(1,1,1,ibox),1,
     1          hess(1,1,1,curbox),1)
               
            iboxtocurbox(ibox) = curbox

            curbox = curbox + 1
            endif
         enddo
         endif
c
c     now add new boxes to level ilev
         do i=1,nblock(ilev)
c           loop over all new boxes on level ilev
            ibox=nboxid(ilevstart(ilev)+i)
            ilevel(curbox) = tilevel(ibox)
            nchild(curbox) = tnchild(ibox)
            do k=1,ndim
               centers(k,curbox) = tcenters(k,ibox)
            enddo

            call dcopy_f77(np,tpot(1,1,ibox),1,pot(1,1,curbox),1)
            call dcopy_f77(np,tcoefsp(1,1,ibox),1,coefsp(1,1,curbox),1)
            if (ifpgh.ge.2) call dcopy_f77(ng,tgrad(1,1,1,ibox),1,
     1          grad(1,1,1,curbox),1)
            if (ifpgh.ge.3) call dcopy_f77(nh,thess(1,1,1,ibox),1,
     1          hess(1,1,1,curbox),1)

            iboxtocurbox(ibox) = curbox

            curbox = curbox + 1
         enddo
c         write(*,*) ilev, curbox-1, laddr(1,ilev), laddr(2,ilev)
         laddr(2,ilev) = curbox-1
      enddo
      endif

c
c     Handle the parent children part of the tree 
c     using the mapping iboxtocurbox
c
c     I believe nlevels should be the updated 
c     actual number of levels here
c
      do ibox=1,laddr(2,nlevels)
         curbox=iboxtocurbox(ibox)
         if(tiparent(ibox).eq.-1) iparent(curbox) = -1
         if(tiparent(ibox).gt.0) 
     1       iparent(curbox) = iboxtocurbox(tiparent(ibox))
         do i=1,mc
            if(tichild(i,ibox).eq.-1) ichild(i,curbox) = -1
            if(tichild(i,ibox).gt.0) 
     1      ichild(i,curbox) = iboxtocurbox(tichild(i,ibox))
         enddo
      enddo
 

c      write(*,*)'time of 3 process',t2-t1,t3-t2,t4-t3,t4-t1
      return

      end subroutine





c     a modified version of vol_tree_mem
c     which works better for the adaptive case
c     with mlevels, mboxes
c     (mlevels = 2* nlevels +1, shall we?)
c     instead of the actual number of boxes
c     and actual number of levels
c
      subroutine vol_tree_mem2(ndim,ipoly,iperiod,eps,zk,boxlen,norder,
     1    iptype,eta,fun,nd,dpars,zpars,ipars,ifnewtree,mboxes,mlevels,
     2    ltree,rintl)
c
c      get memory requirements for the tree
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
c        ifnewtree - interger
c           0: tree is unchanged; 
c           1: tree will be changed and thus needs extra memory
c        output:
c           nlevels - integer
c             number of levels
c           mboxes - integer
c             max number of boxes
c           nlboxes - integer
c             number of leaf boxes
c             not returned in this version
c           ltree - integer
c             length of tree
c           rintl(0:nlevels) - real *8
c             lp norm to scale the functions by
c             (on input rintl should be of size(0:200)
c
c      
c

      implicit none
      real *8 eps,boxlen,eta,dpars(*),delta
      real *8, allocatable :: fvals(:,:,:)
      complex *16 zpars(*),zk
      integer nd,ipars(*),iptype
      integer ltree,ifnewtree
      integer ndim,ipoly,iperiod,nlevels,mboxes,norder
      integer mlevels

      integer, allocatable :: laddr(:,:),ilevel(:),iparent(:),nchild(:)
      integer, allocatable :: ichild(:,:),ncoll(:),icoll(:,:)
      real *8, allocatable :: centers(:,:)
      integer, allocatable :: nbors(:,:),nnbors(:),ifrefine(:)

      integer, allocatable :: ilevel2(:),iparent2(:),nchild2(:),
     1    ichild2(:,:)
      real *8, allocatable :: centers2(:,:),fvals2(:,:,:)

      integer nbmax,nlmax,npbox,npc,mc,mnbors
      real *8, allocatable :: grid(:,:),qwts(:)
      real *8, allocatable:: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8 xy(ndim)
      real *8, allocatable :: wts2(:),xref2(:,:)
      real *8 rintl(0:200)
      real *8 rint
      real *8, allocatable :: rintbs(:),rintbs2(:)
      integer i,itype,j,npols

      real *8, allocatable :: fval1(:,:,:,:),centerstmp(:,:,:)
      real *8, allocatable :: boxsize(:)
      real *8, allocatable :: rmask(:)
      integer, allocatable :: irefinebox(:)

      real *8 rsc,rsum,ra,utmp,vtmp
      integer nbloc,nbctr,nbadd,irefine,ilev,ifirstbox,ilastbox
      integer nbtot,iii,idim,iper,isep,nrefine

      external fun
      nbmax = 100000
      nbmax = 10 000
      nlmax = 200

      allocate(boxsize(0:nlmax))

      mc=2**ndim
      mnbors=3**ndim
      npbox = norder**ndim
      
      allocate(laddr(2,0:nlmax),ilevel(nbmax),iparent(nbmax))
      allocate(nchild(nbmax),ichild(mc,nbmax))

      allocate(fvals(nd,npbox,nbmax),centers(ndim,nbmax))

      allocate(rintbs(nbmax))

c
c      set tree info for level 0
c
      laddr(1,0) = 1
      laddr(2,0) = 1
      ilevel(1) = 0
      iparent(1) = -1
      nchild(1) = 0
      do i=1,mc
        ichild(i,1) = -1
      enddo

      do i=1,ndim
         centers(i,1) = 0
      enddo
c
c     compute rmask and rsum for estimating whether the function is resolved
      allocate(rmask(npbox))

      call tens_prod_get_rmask(ndim,iptype,norder,npbox,
     1    rmask,rsum)
      
c
      allocate(grid(ndim,npbox))
c
c     Generate a grid on the box [-1/2,1/2]^2
c
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
      
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
      allocate(wts2(npbox),xref2(ndim,npbox))

      itype = 1
      call polytens_exps_nd(ndim,ipoly,itype,norder,'f',xref2,
     1    utmp,1,vtmp,1,wts2)
c
c       compute fvals at the grid
c

      boxsize(0) =boxlen

      rint = 0
      rintbs(1) = 0


c
c   note extra factor of 4 sincee wts2 are on [-1,1]^2 
c   as opposed to [-1/2,1/2]^2
c
      rsc = boxlen**2/mc

      do i=1,npbox
        do j=1,ndim
           xy(j) = grid(j,i)*boxlen
        enddo
        call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,1))
        if(iptype.eq.0) then
          do idim=1,nd
            if(abs(fvals(idim,i,1)).gt.rintbs(1)) rintbs(1) = 
     1          abs(fvals(idim,i,1))
          enddo
        endif

        if(iptype.eq.1) then
          do idim=1,nd
            rintbs(1) = rintbs(1) + abs(fvals(idim,i,1))*wts2(i)*rsc
          enddo
        endif

        if(iptype.eq.2) then
          do idim=1,nd
            rintbs(1) = rintbs(1) + fvals(idim,i,1)**2*wts2(i)*rsc
          enddo
        endif
      enddo

      if(iptype.eq.0.or.iptype.eq.1) rint = rintbs(1)
      if(iptype.eq.2) rint = sqrt(rintbs(1))

      rintl(0) = rint
      
      nbctr = 1



     

      do ilev=0,nlmax-1
        irefine = 0

        ifirstbox = laddr(1,ilev) 
        ilastbox = laddr(2,ilev)

        nbloc = ilastbox-ifirstbox+1

        allocate(fval1(nd,npbox,mc,nbloc))
        allocate(centerstmp(ndim,mc,nbloc))
        allocate(irefinebox(nbloc))

        
        if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**ndim)
        if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**ndim)
        if(iptype.eq.0) rsc = 1.0d0

        rsc = rsc*rint
ccc     what are these outputs?
ccc        print *, ilev, rint,rsc

        call vol_tree_find_box_refine(ndim,nd,iptype,eta,eps,zk,
     1      norder,npbox,fvals,npols,umat,rmask,rsum,boxsize(ilev),
     2      nbmax,ifirstbox,nbloc,rsc,irefinebox,irefine)
     

c
c
c          figure out if current set of boxes is sufficient
c

        nbadd = 0 
        do i=1,nbloc
          if(irefinebox(i).eq.1) nbadd = nbadd+mc
        enddo

        nbtot = nbctr+nbadd

c
c         if current memory is not sufficient reallocate
c
        if(nbtot.gt.nbmax) then
c          print *, "Reallocating"
          allocate(centers2(ndim,nbmax),ilevel2(nbmax),iparent2(nbmax))
          allocate(nchild2(nbmax),ichild2(mc,nbmax))
          allocate(fvals2(nd,npbox,nbmax),rintbs2(nbmax))

          call vol_tree_copy(ndim,nd,nbctr,npbox,centers,ilevel,iparent,
     1            nchild,ichild,fvals,centers2,ilevel2,iparent2,nchild2,
     2            ichild2,fvals2)
          call dcopy_f77(nbctr,rintbs,1,rintbs2,1)

          deallocate(centers,ilevel,iparent,nchild,ichild,fvals,rintbs)

          nbmax = nbtot
          allocate(centers(ndim,nbmax),ilevel(nbmax),iparent(nbmax))
          allocate(nchild(nbmax),ichild(mc,nbmax),fvals(nd,npbox,nbmax))
          allocate(rintbs(nbmax))

          call vol_tree_copy(ndim,nd,nbctr,npbox,centers2,ilevel2,
     1        iparent2,nchild2,ichild2,fvals2,centers,ilevel,iparent,
     2        nchild,ichild,fvals)
          call dcopy_f77(nbctr,rintbs2,1,rintbs,1)

          deallocate(centers2,ilevel2,iparent2,nchild2,ichild2,fvals2)
          deallocate(rintbs2)
        endif

        if(irefine.eq.1) then
          boxsize(ilev+1) = boxsize(ilev)/2
          laddr(1,ilev+1) = nbctr+1

          call vol_tree_refine_boxes(ndim,irefinebox,nd,npbox,fvals,
     1      fun,dpars,zpars,ipars,grid,nbmax,ifirstbox,nbloc,centers,
     2      boxsize(ilev+1),nbctr,ilev+1,ilevel,iparent,nchild,ichild)

          rsc = boxsize(ilev+1)**ndim/mc
          call update_rints(ndim,nd,npbox,nbmax,fvals,ifirstbox,nbloc,
     1       iptype,nchild,ichild,wts2,rsc,rintbs,rint)
          rintl(ilev+1) = rint
          
          laddr(2,ilev+1) = nbctr
        else
          exit
        endif

        deallocate(fval1,irefinebox,centerstmp)
      enddo

      mboxes = nbctr
      nlevels = ilev

      if(nlevels.ge.2) then

        nbtot = 2*mc*mboxes
        if(nbtot.gt.nbmax) then
c          print *, "Reallocating"
          allocate(centers2(ndim,nbmax),ilevel2(nbmax),iparent2(nbmax))
          allocate(nchild2(nbmax),ichild2(mc,nbmax))
          allocate(fvals2(nd,npbox,nbmax),rintbs2(nbmax))

          call vol_tree_copy(ndim,nd,mboxes,npbox,centers,ilevel,
     1       iparent,nchild,ichild,fvals,centers2,ilevel2,
     2       iparent2,nchild2,ichild2,fvals2)
          call dcopy_f77(mboxes,rintbs,1,rintbs2,1)

          deallocate(centers,ilevel,iparent,nchild,ichild,fvals,rintbs)

          nbmax = nbtot
          allocate(centers(ndim,nbmax),ilevel(nbmax),iparent(nbmax))
          allocate(nchild(nbmax),ichild(mc,nbmax),fvals(nd,npbox,nbmax))
          allocate(rintbs(nbmax))

          call vol_tree_copy(ndim,nd,mboxes,npbox,centers2,ilevel2,
     1        iparent2,nchild2,ichild2,fvals2,centers,ilevel,
     2        iparent,nchild,ichild,fvals)
          call dcopy_f77(mboxes,rintbs2,1,rintbs,1)

          deallocate(centers2,ilevel2,iparent2,nchild2,ichild2,fvals2)
          deallocate(rintbs2)
        endif

        allocate(nnbors(nbmax))
        allocate(nbors(mnbors,nbmax))

        call computecoll(ndim,nlevels,mboxes,laddr,boxsize,centers,
     1        iparent,nchild,ichild,iperiod,nnbors,nbors)

        if(nlevels.ge.2) then
          call vol_tree_fix_lr(ndim,iperiod,fun,nd,dpars,zpars,ipars,
     1      norder,npbox,fvals,grid,centers,nlevels,mboxes,boxsize,
     2      nbmax,nlmax,laddr,ilevel,iparent,nchild,ichild,nnbors,nbors)
        endif
      endif

c     triple the number of boxes for possible refinement
c     added by Shidong Jiang 04/13/2022
      if (ifnewtree.eq.1) mboxes = mboxes*3

      ltree = (4+mc+mnbors)*mboxes + 2*(nlevels+1)
c
      mlevels = 2*nlevels + 1
ccc   add in the following
      do ilev=nlevels+1, 200
        rintl(ilev)=rint
      enddo

      return
      end
c
c
c
c
c
c     a modified version of vol_tree_build
c     which works better for the adaptive case
c     mboxes and mlevels are inputs
c     nboxes and nlevels are outputs
      subroutine vol_tree_build2(ndim,ipoly,iperiod,eps,zk,boxlen,
     1    norder,iptype,eta,fun,nd,dpars,zpars,ipars,rintl,nboxes,
     2    mboxes,nlevels,mlevels,ltree,itree,iptr,centers,boxsize,
     3    fvals)
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
c         
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
c          function values at discretization nodes
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
      integer mlevels,mboxes
      integer iptr(8),ltree
      integer itree(ltree),ier
      real *8 fvals(nd,norder**ndim,mboxes),centers(ndim,mboxes)
      real *8, allocatable :: fval1(:,:,:,:),centerstmp(:,:,:)
      integer, allocatable :: irefinebox(:),ifrefine(:)
      real *8 boxsize(0:nlevels)
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8, allocatable :: rmask(:)
      real *8 rintl(0:mlevels)
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
c     set nlevels here?
      nlevels=(mlevels-1)/2
c      
c-----------------------------------------
c     like this?
      iptr(1) = 1
      iptr(2) = 2*(mlevels+1)+1
      iptr(3) = iptr(2) + mboxes
      iptr(4) = iptr(3) + mboxes
      iptr(5) = iptr(4) + mboxes
      iptr(6) = iptr(5) + mc*mboxes
      iptr(7) = iptr(6) + mboxes
      iptr(8) = iptr(7) + mnbors*mboxes
c-----------------------------------------
c
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
c     what is rmask and rsum?      
c     chose indices based on iptype

c
c       Reset nlevels, nboxes
c
      nbctr = 1

c     nlevels correct here?
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
ccccccc     why this line??
ccccccc     lost this line in my subroutine:
ccccccc     vol_tree_adap_fun
        rsc = rsc*rintl(ilev)
c       nboxes or mboxes??
        call vol_tree_find_box_refine(ndim,nd,iptype,eta,eps,zk,
     1      norder,npbox,fvals,npols,umat,rmask,rsum,boxsize(ilev),
     2      mboxes,ifirstbox,nbloc,rsc,irefinebox,irefine)
        

        if(irefine.eq.1) then
          boxsize(ilev+1) = boxsize(ilev)/2
          itree(2*ilev+3) = nbctr+1

c         nboxes or mboxes?
          call vol_tree_refine_boxes(ndim,irefinebox,nd,npbox,fvals,
     1      fun,dpars,zpars,ipars,grid,mboxes,ifirstbox,nbloc,centers,
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

c      if(1 .eq. 2) then
      if(nlevels.ge.2) then
         call vol_tree_fix_lr(ndim,iperiod,fun,nd,dpars,zpars,ipars,
     1       norder,npbox,fvals,grid,centers,nlevels,nboxes0,boxsize,
     2       mboxes,nlevels,itree(iptr(1)),itree(iptr(2)),
     3       itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     4       itree(iptr(6)),itree(iptr(7)))
      endif
c      endif
c
      nboxes = itree(2*nlevels+2)
      if(nboxes .ne. nboxes0) then
        write(*,*) 'inconsistent nboxes and nboxes0'
      endif
c
c      call prinf('mboxes in tree_build=*',mboxes,1)
c      call prinf('nboxes=*',nboxes,1)
c      call prinf('nboxes0=*',nboxes0,1)
cccc      call prinf('nlevels=*',nlevels,1)
      return
      end
c
c
c      
c
c
c-------------------------------------------------
c     subroutine vol_tree_set_vals
c
c     sample a function on a given tree
c
c-------------------------------------------------
c
      subroutine vol_tree_set_vals(nd,ndim,ipoly,
     1    norder,npbox,nboxes,mboxes,nlevels,
     2    mlevels,ltree,itree,iptr,centers,boxsize,
     3    dpars,zpars,ipars,fun,fvals)
c
c    Input:
c     nd - integer
c            number of right hand sides
c     ndim - integer
c            dimension of the underlying space
c     ipoly - integer
c            0: Legendre polynomials
c            1: Chebyshev polynomials
c     norder - integer
c           order of expansions for input function value array
c     npbox - integer
c           number of points per box where potential is to be dumped = (norder**ndim)
c     nboxes - integer
c            number of boxes
c
c     input and output:
c     nlevels - integer
c            number of levels
c     ltree - integer
c            length of array containing the tree structure
c     itree - integer(ltree)
c            array containing the tree structure
c     iptr - integer(8)
c            pointer to various parts of the tree structure
c           iptr(1) - laddr
c           iptr(2) - ilevel
c           iptr(3) - iparent
c           iptr(4) - nchild
c           iptr(5) - ichild
c           iptr(6) - ncoll
c           iptr(7) - coll
c           iptr(8) - ltree
c     centers - double precision (ndim,mboxes)
c           xyz coordintes of boxes in the tree structure
c     boxsize - double precision (0:nlevels)
c           size of boxes at each of the levels
c
c     dpars - double precision
c           real parameters for function evaluation
c     zpars - double complex
c           complex parameters for function evaluation
c     ipars - integer
c           integer parameters for function evaluation
c     fun - a function handle, the funciton to be
c           resolved on the tree
c
c----------------
c     Output:
c
c     fvals: the function sampled on the tree
c
c     attention:
c     mboxes: here it means the max number of boxes
c     nlevels: here it means the actual number of levels on input
c     
      implicit real*8 (a-h,o-z) 
      integer nd,ndim,mboxes,nlevels
      integer iptr(8),ltree,ipoly
      integer itree(ltree),norder,npbox
      real *8 eps, eta
      real *8 fvals(nd,npbox,mboxes),fcoefs(nd,npbox,mboxes)
      real *8 centers(ndim,mboxes),boxsize(0:mlevels)
      external fun
      complex *16 zk
c     allocatables
c     grid related stuff
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8, allocatable :: rmask(:), umat_nd(:,:,:)
      real *8 xy(ndim)
      real *8, allocatable:: coefsp(:,:,:)
      integer isgn(ndim,2**ndim)
      integer, allocatable:: nblock(:)
      integer, allocatable:: ilevstart(:)
      integer, allocatable:: nboxid(:), ifdelete(:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
c
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
      allocate(umat_nd(norder,norder,ndim))
c
      call get_child_box_sign(ndim,isgn)
c
      allocate(grid(ndim,npbox))
      allocate(grad(nd,ndim,npbox,1))
      allocate(hess(nd,ndim*(ndim+1)/2,npbox,1))
c
c     initialize the grid related stuff
c
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)
c
      norder2=norder*norder
      do i=1,ndim
         call dcopy_f77(norder2,umat,1,umat_nd(1,1,i),1)
      enddo
c
      do i=1,norder
        xq(i) = xq(i)/2
      enddo
c
      if (ndim.eq.1) then
         do i=1,norder
            grid(1,i)=xq(i)
         enddo
      elseif (ndim.eq.2) then
         call mesh2d(xq,norder,xq,norder,grid)
      elseif (ndim.eq.3) then
         call mesh3d(xq,norder,xq,norder,xq,norder,grid)
      endif
c
c      write(*,*) 'grid:'
c      do i=1, ndim
c      do j=1, npbox
c        write(*,*) i,j,grid(i,j)
c        write(201,*) i,j,grid(i,j)
c      enddo
c      enddo
c     grid looks correct
c
c--------------------------------------------------
c     on the currect tree
c     eval fun on the grid of leaf nodes
      do ilev=0, nlevels
        do ibox = itree(2*ilev+1), itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild .eq. 0) then
            do i=1,npbox
              do j=1,ndim
                xy(j) = grid(j,i)*boxsize(ilev)+centers(j,ibox)
c                write(*,*) i, j, xy(j)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,ibox))
c              write(*,*) i, ibox, fvals(1,i,ibox)
            enddo
          endif
        enddo
      enddo
c


      end subroutine
c
c
c      
c
c
c------------------------------------------------------
c     a modified version of vol_tree_coarsen
c     which checks the resolution of two functions
c     one given by a function handle: fun
c     another given by sampled function values: uvals
c     coarsen if both functions are over-resolved
c
c     fvals, fcoefs -> uvals, ucoefs, fun, grid
c     interp u and eval f
c
c------------------------------------------------------
c
      subroutine vol_tree_coarsen_fv(nd,ndim,reps,ipoly,norder,npbox,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    ifpgh,uvals,ucoefs,fun,grid,grad,hess,
     3    nblock,nboxid,ndelboxes,ifdelete)

c    This subroutine marks all boxes that can be coarsened
c    and computes the expansion coefficients for parent boxes
c    when its children are to be deleted in a subsequent call
c    to tree_reorg.
c
c    input:
c     nd - integer
c            number of right hand sides
c     ndim - integer
c            dimension of the underlying space
c     reps - double precision
c            scaled desired precision
c     ipoly - integer
c            0: Legendre polynomials
c            1: Chebyshev polynomials
c     norder - integer
c           order of expansions for input function value array
c     npbox - integer
c           number of points per box where potential is to be dumped = (norder**ndim)
c     nboxes - integer
c            number of boxes
c
c     input and output:
c     nlevels - integer
c            number of levels
c     ltree - integer
c            length of array containing the tree structure
c     itree - integer(ltree)
c            array containing the tree structure
c     iptr - integer(8)
c            pointer to various parts of the tree structure
c           iptr(1) - laddr
c           iptr(2) - ilevel
c           iptr(3) - iparent
c           iptr(4) - nchild
c           iptr(5) - ichild
c           iptr(6) - ncoll
c           iptr(7) - coll
c           iptr(8) - ltree
c     centers - double precision (ndim,nboxes)
c           xyz coordintes of boxes in the tree structure
c     boxsize - double precision (0:nlevels)
c           size of boxes at each of the levels
c
c     input/output:
c     fvals - double precision (nd,npbox,nboxes)
c           function values of the given tree data on each leaf box
c     fcoefs - double precision (nd,npbox,nboxes)
c           expansion coefficients of the given tree data on each leaf box
c     grad - double precision (nd,ndim,npbox,nboxes)
c           gradient 
c     hess - double precision (nd,npbox,nboxes)
c           hessian
c
c     output:
c     nblock - number of boxes to be removed on each level
c     nboxid - indices of the boxes to be removed
c     ndelboxes - total number of removable boxes
c     ifdelete - integer nboxes
c                    1 - box can be deleted
c                    0 - box should be kept
c
c
c
      implicit real *8 (a-h,o-z)
      real *8 reps
      integer nboxes,nlevels
      integer iptr(8),ltree
      integer itree(ltree),norder,npbox
      integer nblock(0:nlevels),nboxid(nboxes),ifdelete(nboxes)

      real *8 uvals(nd,npbox,nboxes),ucoefs(nd,npbox,nboxes)
      real *8 grad(nd,ndim,npbox,*),hess(nd,ndim*(ndim+1)/2,npbox,*)

      real *8 centers(ndim,nboxes),boxsize(0:nlevels)

      real *8 uvals2(nd,npbox/2**ndim), grid(ndim,npbox)
      real*8 xy(ndim)
      external fun
      real *8, allocatable :: polyv(:,:,:,:) 
      real *8, allocatable :: gvals2(:,:,:) 
      real *8, allocatable :: hvals2(:,:,:) 
      real *8, allocatable :: gcoefs(:,:,:) 
      real *8, allocatable :: hcoefs(:,:,:) 
      real *8, allocatable :: fvals(:,:,:)
      real *8, allocatable :: fcoefs(:,:,:)

      real *8 umat_nd(norder,norder,ndim)
      real *8 rmask(npbox)
      integer isgn(ndim,2**ndim)

      mc=2**ndim
      n2=norder/2

      np2=npbox/2**ndim
      nhess=ndim*(ndim+1)/2

      ng=nd*ndim
      nh=nd*nhess
c

      allocate(fvals(nd,npbox,nboxes),fcoefs(nd,npbox,nboxes))
      
      if (ifpgh.ge.2) then
         allocate(gvals2(nd,ndim,np2))
         allocate(gcoefs(nd,ndim,npbox))
      endif
      if (ifpgh.ge.3) then
         allocate(hvals2(nd,nhess,np2))
         allocate(hcoefs(nd,nhess,npbox))
      endif
      
      call get_child_box_sign(ndim,isgn)
      allocate(polyv(norder,n2,ndim,mc))
      call get_c2p_interp_matrices(ndim,ipoly,norder,isgn,polyv)
      call get_val2coefs_matrices(ndim,ipoly,norder,umat_nd)

      do i=1,nboxes
         ifdelete(i)=0
      enddo

      do i=0,nlevels
         nblock(i)=0
      enddo
      
      iptype=2
cccc      eta=1.0d0
      eta=0.0d0

      call tens_prod_get_rmask(ndim,iptype,norder,npbox,
     1    rmask,rsum)
      
      ndelboxes=0
      
      itype=0
      do ilev = nlevels-1,0,-1
        sc = boxsize(ilev)/2
        rscale=sc**eta
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.mc) then
c           first, make sure all children are leaf boxes
            nleafkids=0
            do j=1,mc
               jbox=itree(iptr(5)+mc*(ibox-1)+j-1)
               if (itree(iptr(4) + jbox-1).eq.0) then
                  nleafkids=nleafkids+1
               endif
            enddo
c           second, check whether the function can be merged
            if (nleafkids.eq.mc) then
               do j=1,mc
                  jbox=itree(iptr(5)+mc*(ibox-1)+j-1)
                  
                  call ortho_eval_nd(ndim,nd,norder,ucoefs(1,1,jbox),
     1                n2,uvals2,polyv(1,1,1,j))
                  call tree_data_c2p_copy_nd(ndim,nd,norder,uvals2,
     1                uvals(1,1,ibox),j)
               enddo
               call ortho_trans_nd(ndim,nd,itype,norder,uvals(1,1,ibox),
     1             ucoefs(1,1,ibox),umat_nd)
               call fun_err(nd,npbox,ucoefs(1,1,ibox),
     1            rmask,iptype,rscale,erra)
cccc               print *, 'erra=',erra, eps*rsc
cccc           add: sample the function: fun
c                   on the parent box: ibox
               do i=1, npbox
                 do j=1, ndim
                   xy(j) = grid(j,i)*boxsize(ilev)+centers(j,ibox)
                 enddo
                 call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,ibox))
               enddo
cccc           convert to coefficients: fcoefs
               call ortho_trans_nd(ndim,nd,itype,norder,fvals(1,1,ibox),
     1             fcoefs(1,1,ibox),umat_nd)
cccc           check the error: errb
               call fun_err(nd,npbox,fcoefs(1,1,ibox),
     1            rmask,iptype,rscale,errb)
c               
               erra = erra/rsum
               errb = errb/rsum
c              erra: u, errb: f
               if((erra.le.reps).and.(errb.le.reps)) then
c                  write(225,*) 'coarsening', ibox, erra, errb, reps
                  if (ifpgh.ge.2) then
                     do j=1,mc
                        jbox=itree(iptr(5)+mc*(ibox-1)+j-1)
                  
                        call ortho_trans_nd(ndim,ng,itype,norder,
     1                      grad(1,1,1,jbox),gcoefs,umat_nd)
                        call ortho_eval_nd(ndim,ng,norder,
     1                      gcoefs,n2,gvals2,polyv(1,1,1,j))
                        call tree_data_c2p_copy_nd(ndim,ng,norder,
     1                      gvals2,grad(1,1,1,ibox),j)
                     enddo
                  endif
                  if (ifpgh.ge.3) then
                     do j=1,mc
                        jbox=itree(iptr(5)+mc*(ibox-1)+j-1)
                  
                        call ortho_trans_nd(ndim,nh,itype,norder,
     1                      hess(1,1,1,jbox),hcoefs,umat_nd)
                        call ortho_eval_nd(ndim,nh,norder,
     1                      hcoefs,n2,hvals2,polyv(1,1,1,j))
                        call tree_data_c2p_copy_nd(ndim,nh,norder,
     1                      hvals2,hess(1,1,1,ibox),j)
                     enddo
                  endif
                  
                  do j=1,mc
                     jbox=itree(iptr(5)+mc*(ibox-1)+j-1)

                     ndelboxes=ndelboxes+1
                     nboxid(ndelboxes)=jbox
                     nblock(ilev+1)=nblock(ilev+1)+1
                     
c                    parent
                     itree(iptr(3)+jbox-1)=-1
c                    nchild of jbox
                     itree(iptr(4)+jbox-1)=0
c                    jbox's children
                     do l=1,mc
                        itree(iptr(5)+mc*(jbox-1)+l-1) = -1
                     enddo
                     
                     ifdelete(jbox)=1
                  enddo
c                 nchild of ibox
                  itree(iptr(4)+ibox-1)=0
                  do j=1,mc
                     itree(iptr(5)+mc*(ibox-1)+j-1)=-1
                  enddo
               endif
             endif
          endif
        enddo
      enddo
      
      return
      end
c
c
c
c
c------------------------------------------------------
c     another modified version of vol_tree_coarsen
c     which checks the resolution of two functions
c     both given by sampled function values
c
c     fvals, fcoefs -> fvals, fcoefs, uvals, ucoefs,
c     interp both to check the resolution
c
c------------------------------------------------------
c
      subroutine vol_tree_coarsen_f2(nd,ndim,reps,ipoly,
     1    norder,npbox,nboxes,nlevels,ltree,itree,iptr,
     2    centers,boxsize,ifpgh,fvals,fcoefs,uvals,ucoefs,
     3    grad,hess,nblock,nboxid,ndelboxes,ifdelete)

c    This subroutine marks all boxes that can be coarsened
c    and computes the expansion coefficients for parent boxes
c    when its children are to be deleted in a subsequent call
c    to tree_reorg.
c
c    input:
c     nd - integer
c            number of right hand sides
c     ndim - integer
c            dimension of the underlying space
c     reps - double precision
c            scaled desired precision
c     ipoly - integer
c            0: Legendre polynomials
c            1: Chebyshev polynomials
c     norder - integer
c           order of expansions for input function value array
c     npbox - integer
c           number of points per box where potential is to be dumped = (norder**ndim)
c     nboxes - integer
c            number of boxes
c
c     input and output:
c     nlevels - integer
c            number of levels
c     ltree - integer
c            length of array containing the tree structure
c     itree - integer(ltree)
c            array containing the tree structure
c     iptr - integer(8)
c            pointer to various parts of the tree structure
c           iptr(1) - laddr
c           iptr(2) - ilevel
c           iptr(3) - iparent
c           iptr(4) - nchild
c           iptr(5) - ichild
c           iptr(6) - ncoll
c           iptr(7) - coll
c           iptr(8) - ltree
c     centers - double precision (ndim,nboxes)
c           xyz coordintes of boxes in the tree structure
c     boxsize - double precision (0:nlevels)
c           size of boxes at each of the levels
c
c     input/output:
c     fvals - double precision (nd,npbox,nboxes)
c           function values of the given tree data on each leaf box
c     fcoefs - double precision (nd,npbox,nboxes)
c           expansion coefficients of the given tree data on each leaf box
c     grad - double precision (nd,ndim,npbox,nboxes)
c           gradient 
c     hess - double precision (nd,npbox,nboxes)
c           hessian
c
c     output:
c     nblock - number of boxes to be removed on each level
c     nboxid - indices of the boxes to be removed
c     ndelboxes - total number of removable boxes
c     ifdelete - integer nboxes
c                    1 - box can be deleted
c                    0 - box should be kept
c
c
c
      implicit real *8 (a-h,o-z)
      real *8 reps
      integer nboxes,nlevels
      integer iptr(8),ltree
      integer itree(ltree),norder,npbox
      integer nblock(0:nlevels),nboxid(nboxes),ifdelete(nboxes)

      real *8 fvals(nd,npbox,nboxes),fcoefs(nd,npbox,nboxes)
      real *8 uvals(nd,npbox,nboxes),ucoefs(nd,npbox,nboxes)
      real *8 grad(nd,ndim,npbox,*),hess(nd,ndim*(ndim+1)/2,npbox,*)

      real *8 centers(ndim,nboxes),boxsize(0:nlevels)

      real *8 fvals2(nd,npbox/2**ndim)
      real *8 uvals2(nd,npbox/2**ndim)

      real *8, allocatable :: polyv(:,:,:,:) 
      real *8, allocatable :: gvals2(:,:,:) 
      real *8, allocatable :: hvals2(:,:,:) 
      real *8, allocatable :: gcoefs(:,:,:) 
      real *8, allocatable :: hcoefs(:,:,:) 

      real *8 umat_nd(norder,norder,ndim)
      real *8 rmask(npbox)
      integer isgn(ndim,2**ndim)

      mc=2**ndim
      n2=norder/2

      np2=npbox/2**ndim
      nhess=ndim*(ndim+1)/2

      ng=nd*ndim
      nh=nd*nhess
      
      if (ifpgh.ge.2) then
         allocate(gvals2(nd,ndim,np2))
         allocate(gcoefs(nd,ndim,npbox))
      endif
      if (ifpgh.ge.3) then
         allocate(hvals2(nd,nhess,np2))
         allocate(hcoefs(nd,nhess,npbox))
      endif
      
      call get_child_box_sign(ndim,isgn)
      allocate(polyv(norder,n2,ndim,mc))
      call get_c2p_interp_matrices(ndim,ipoly,norder,isgn,polyv)
      call get_val2coefs_matrices(ndim,ipoly,norder,umat_nd)

      do i=1,nboxes
         ifdelete(i)=0
      enddo

      do i=0,nlevels
         nblock(i)=0
      enddo
      
      iptype=2
cccc      eta=1.0d0
      eta=0.0d0

      call tens_prod_get_rmask(ndim,iptype,norder,npbox,
     1    rmask,rsum)
      
      ndelboxes=0
      
      itype=0
      do ilev = nlevels-1,0,-1
        sc = boxsize(ilev)/2
        rscale=sc**eta
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.mc) then
c           first, make sure all children are leaf boxes
            nleafkids=0
            do j=1,mc
               jbox=itree(iptr(5)+mc*(ibox-1)+j-1)
               if (itree(iptr(4) + jbox-1).eq.0) then
                  nleafkids=nleafkids+1
               endif
            enddo
c           second, check whether the function can be merged
            if (nleafkids.eq.mc) then
               do j=1,mc
                  jbox=itree(iptr(5)+mc*(ibox-1)+j-1)
                  
                  call ortho_eval_nd(ndim,nd,norder,fcoefs(1,1,jbox),
     1                n2,fvals2,polyv(1,1,1,j))
                  call tree_data_c2p_copy_nd(ndim,nd,norder,fvals2,
     1                fvals(1,1,ibox),j)
c
                  call ortho_eval_nd(ndim,nd,norder,ucoefs(1,1,jbox),
     1                n2,uvals2,polyv(1,1,1,j))
                  call tree_data_c2p_copy_nd(ndim,nd,norder,uvals2,
     1                uvals(1,1,ibox),j)
               enddo
               call ortho_trans_nd(ndim,nd,itype,norder,fvals(1,1,ibox),
     1             fcoefs(1,1,ibox),umat_nd)
               call fun_err(nd,npbox,fcoefs(1,1,ibox),
     1            rmask,iptype,rscale,erra)
c
               call ortho_trans_nd(ndim,nd,itype,norder,uvals(1,1,ibox),
     1             ucoefs(1,1,ibox),umat_nd)
               call fun_err(nd,npbox,ucoefs(1,1,ibox),
     1            rmask,iptype,rscale,errb)

cccc               print *, 'erra=',erra, eps*rsc
c               
               erra = erra/rsum
               errb = errb/rsum
c
c               if((ibox .eq. 63).or.(ibox .eq. 66)) then
c                 write(325,*) 'no-coarsening', ibox, erra, reps
c               endif 
               if((erra.le.reps).and.(errb.le.reps) )then
c                  write(225,*) 'coarsening', ibox, erra, errb,reps
                  if (ifpgh.ge.2) then
                     do j=1,mc
                        jbox=itree(iptr(5)+mc*(ibox-1)+j-1)
                  
                        call ortho_trans_nd(ndim,ng,itype,norder,
     1                      grad(1,1,1,jbox),gcoefs,umat_nd)
                        call ortho_eval_nd(ndim,ng,norder,
     1                      gcoefs,n2,gvals2,polyv(1,1,1,j))
                        call tree_data_c2p_copy_nd(ndim,ng,norder,
     1                      gvals2,grad(1,1,1,ibox),j)
                     enddo
                  endif
                  if (ifpgh.ge.3) then
                     do j=1,mc
                        jbox=itree(iptr(5)+mc*(ibox-1)+j-1)
                  
                        call ortho_trans_nd(ndim,nh,itype,norder,
     1                      hess(1,1,1,jbox),hcoefs,umat_nd)
                        call ortho_eval_nd(ndim,nh,norder,
     1                      hcoefs,n2,hvals2,polyv(1,1,1,j))
                        call tree_data_c2p_copy_nd(ndim,nh,norder,
     1                      hvals2,hess(1,1,1,ibox),j)
                     enddo
                  endif
                  
                  do j=1,mc
                     jbox=itree(iptr(5)+mc*(ibox-1)+j-1)

                     ndelboxes=ndelboxes+1
                     nboxid(ndelboxes)=jbox
                     nblock(ilev+1)=nblock(ilev+1)+1
                     
c                    parent
                     itree(iptr(3)+jbox-1)=-1
c                    nchild of jbox
                     itree(iptr(4)+jbox-1)=0
c                    jbox's children
                     do l=1,mc
                        itree(iptr(5)+mc*(jbox-1)+l-1) = -1
                     enddo
                     
                     ifdelete(jbox)=1
                  enddo
c                 nchild of ibox
                  itree(iptr(4)+ibox-1)=0
                  do j=1,mc
                     itree(iptr(5)+mc*(ibox-1)+j-1)=-1
                  enddo
               endif
             endif
          endif
        enddo
      enddo
      
      return
      end
c
c
c
c
c
c
c-------------------------------------------------
c     subroutine vol_tree_adap_fv
c-------------------------------------------------
c
c     a modified version of vol_tree_adap_fun
c
c     given a tree, and a function resolved on the 
c     leaf nodes, and another function (given by a 
c     function handle), reifne and coarsen the tree
c     so that both functions are resolved but not
c     over-resolved on the new tree
c
c     the I/O is pretty much the same as
c     subroutine vol_tree_adap_fun
c
c     except for:
c     uvals: (in/out) a function sampled on the 
c     tensor product grid of each leaf node
c
c     on output, the tree is modified, and uvals
c     is interpolated and rearranged to match 
c     the new tree
c
c     rmk: 1. after refinement, interpolate to the children
c          2. when coarsening, check the resolution of uvals too
c
c-------------------------------------------------
c
      subroutine vol_tree_adap_fv(nd,ndim,eps,ipoly,
     1    iptype,norder,npbox,nboxes,mboxes,nlevels,
     2    mlevels,ltree,itree,iptr,centers,boxsize,
     3    rintl,dpars,zpars,ipars,iperiod,fun,uvals)
      implicit real*8 (a-h,o-z) 
      integer nd,ndim,mboxes,nlevels
      integer iptr(8),ltree,ipoly
      integer itree(ltree),norder,npbox
      integer ipars(1)
      real *8 eps, eta
      real *8 uvals(nd,npbox,mboxes),fcoefs(nd,npbox,mboxes)
      real *8 centers(ndim,mboxes),boxsize(0:mlevels)
      real *8 dpars(1)
      external fun
      complex *16 zk
      complex *16 zpars(1)
c     allocatables
      real*8, allocatable :: fvals(:,:,:)
c     grid related stuff
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8, allocatable :: rmask(:), umat_nd(:,:,:)
      real *8 rintl(0:mlevels)
      real *8 xy(ndim)
      real *8, allocatable:: coefsp(:,:,:), ucoefs(:,:,:)
      real *8, allocatable:: polyvc(:,:,:,:)
c     flags for refinement and coarsening
      integer isgn(ndim,2**ndim)
c      integer nblock(0:nlevels), ilevstart(0:nlevels)
      integer, allocatable:: nblock(:)
      integer, allocatable:: ilevstart(:)
      integer, allocatable:: nboxid(:), ifdelete(:)
      integer, allocatable:: ifrefine(:)
      integer, allocatable:: irefinelev(:), imaxrefinelev(:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
c
c      write(*,*) nd, ndim, eps, ipoly, iptype
c      write(*,*) norder, npbox
c      write(*,*) nboxes, mboxes, nlevels, mlevels, ltree
c      write(*,*) dpars(1), iperiod
c
      mc=2**ndim
      mnbors=3**ndim
      npols = norder**ndim
c
      allocate(fvals(nd,npbox,mboxes))
c
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
      allocate(umat_nd(norder,norder,ndim))
      allocate(polyvc(norder,norder,ndim,mc))
c
      call get_child_box_sign(ndim,isgn)
      call get_p2c_interp_matrices(ndim,ipoly,norder,isgn,polyvc)
c
      allocate(grid(ndim,npbox))
      allocate(grad(nd,ndim,npbox,1))
      allocate(hess(nd,ndim*(ndim+1)/2,npbox,1))
    
c
c     initialize the grid related stuff
c
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)
c
      norder2=norder*norder
      do i=1,ndim
         call dcopy_f77(norder2,umat,1,umat_nd(1,1,i),1)
      enddo

c
      do i=1,norder
        xq(i) = xq(i)/2
      enddo
c
      if (ndim.eq.1) then
         do i=1,norder
            grid(1,i)=xq(i)
         enddo
      elseif (ndim.eq.2) then
         call mesh2d(xq,norder,xq,norder,grid)
      elseif (ndim.eq.3) then
         call mesh3d(xq,norder,xq,norder,xq,norder,grid)
      endif
 
c
c--------------------------------------------------
c     on the currect tree
c     eval fun on the grid of leaf nodes
      do ilev=0, nlevels
        do ibox = itree(2*ilev+1), itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild .eq. 0) then
            do i=1,npbox
              do j=1,ndim
                xy(j) = grid(j,i)*boxsize(ilev)+centers(j,ibox)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,ibox))
c              write(*,*)i,ibox,fvals(1,i,ibox)
            enddo
          endif
        enddo
      enddo

c
c     compute rmask and rsum for estimating whether the function is resolved
      allocate(rmask(npbox))

      call tens_prod_get_rmask(ndim,iptype,norder,npbox,
     1    rmask,rsum)
c
      if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**ndim)
      if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**ndim)
      if(iptype.eq.0) rsc = 1.0d0
c
c     rsc = rsc*rnorm
c     problematic?
      reps = eps*rsc

c
c     save the original number of boxes
ccc   no need to do this now 
ccc   since nboxes is an input
      nboxes0 = nboxes
      nlevels0 = nlevels
c      write(*,*) 'nboxes0=', nboxes0
c      write(*,*) 'nboxes=', nboxes
c      write(*,*) 'mboxes=', mboxes
c
c      write(*,*) '***'
c
c     get rid of the max refine levels and stuff
      allocate(ifrefine(mboxes))
c
c     mboxes here means the max of nboxes
      do ibox = 1, mboxes
        ifrefine(ibox)=0
      enddo
c
c     convert potential to its polynomial expansion coefficients
c     (for all the leaf nodes)
c
c     add: ucoefs, transform uvals -> ucoefs
c
      allocate(coefsp(nd,npbox,mboxes))
      allocate(ucoefs(nd,npbox,mboxes))
      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,fvals,coefsp,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,uvals,ucoefs,umat_nd)
c
c     now go down the original tree level by level
c     mark boxes that need to be refined
c
      eta=1.0d0
      zk=1.0d0
c
      do ilev=0, nlevels
        sc = boxsize(ilev)/2
        rscale=sc**eta
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
c         loop over boxes of level ilev
          if (itree(iptr(4)+ibox-1).eq.0) then
c           if a leaf box, check the resolution
             call fun_err(nd,npbox,coefsp(1,1,ibox),
     1            rmask,iptype,rscale,erra)
             erra = erra/rsum
c             write(*,*) ibox, ifrefine(ibox), erra,
c     1                    reps*rintl(ilev), rintl(ilev), reps
             if(erra.gt.reps*rintl(ilev)) then
                ifrefine(ibox)=1
             endif
          endif
        enddo
      enddo
c
c     leaf boxes that need to be refined
c     are flaged in ifrefine
c
c-----------------------------------------
c
c      write(*,*) nlevels, mlevels
c
      do ilev = nlevels+1, mlevels
        boxsize(ilev)=boxsize(ilev-1)/2.0d0
c        write(*,*) ilev, boxsize(ilev)
      enddo
c
      allocate(nblock(0:mlevels))
      allocate(ilevstart(0:mlevels+1))
c
      do i=0,mlevels
         nblock(i)=0
      enddo
c
c-----------------------------------------
c
      nnewboxes=0
      do ibox = 1, mboxes
c       change this loop to a level-by-level loop?
c       if no refinement is done on the finest level, exit
c       loop over all boxes, including newly added ones
c
c        write(*,*) ibox, ifrefine(ibox)
        if(ifrefine(ibox) .eq. 1) then
c          write(123,*) ibox, ifrefine(ibox)
          ilev=itree(iptr(2)+ibox-1)
          bsh = boxsize(ilev)/4
c         half boxsize on the children's level?
          rscale=bsh**eta
c         set up nchild of ibox
          itree(iptr(4)+ibox-1)=mc
c         setup for each child
          do j=1, mc
            nnewboxes = nnewboxes+1
            jbox = nboxes0+nnewboxes
c
c           set up centers for the new children
            do k=1,ndim
              centers(k,jbox) = centers(k,ibox)+isgn(k,j)*bsh
            enddo
c           set up parent for the new children
            itree(iptr(3)+jbox-1)=ibox
c
c           set up nchild for the new children
            itree(iptr(4)+jbox-1)=0
c
c           set up children for the new children
            do l=1,mc
              itree(iptr(5)+mc*(jbox-1)+l-1) = -1
            enddo
c
c           set up children for the parent
            itree(iptr(5)+mc*(ibox-1)+j-1)=jbox
c
c           set up level for the new children
            jlev=ilev+1
            itree(iptr(2)+jbox-1) = jlev
c           count the number of boxes on the level jlev
            nblock(jlev)=nblock(jlev)+1
c            write(*,*) 'nblock', jlev, nblock(jlev)
c           update nlevels along the way
            if(jlev .gt. nlevels) then
              nlevels = jlev
            endif
c
c           now compute fvals for
c           the newly generated boxes
c           clearly not the right grid
            do ii=1, npbox
              do jj=1, ndim
                xy(jj) = grid(jj,ii)*boxsize(jlev)+centers(jj,jbox)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,ii,jbox))
c             okay, need to check the error
            enddo
c     
c-------------------
c           add here: interpolate uvals to the new boxes
            call ortho_eval_nd(ndim,nd,norder,
     1           ucoefs(1,1,ibox),norder,uvals(1,1,jbox),
     2           polyvc(1,1,1,j))
c
            call ortho_trans_nd(ndim,nd,0,norder,fvals(1,1,jbox),
     1           coefsp(1,1,jbox),umat_nd)
c          add here: convert uvals to ucoefs on jbox
c
            call ortho_trans_nd(ndim,nd,0,norder,uvals(1,1,jbox),
     1           ucoefs(1,1,jbox),umat_nd)
c
            call fun_err(nd,npbox,coefsp(1,1,jbox),rmask,
     1            iptype,rscale,erra)
            erra = erra/rsum
c           error obtained for the newly added box jbox
c
            if(erra .gt. reps*rintl(jlev)) then
c             flag the new box
              ifrefine(jbox) = 1
            endif
c            write(124,*) erra, reps*rintl(jlev), ifrefine(jbox)
c            write(*,*) '...', erra,reps, rintl(jlev)
c         end of the loop over mc children
          enddo
c       end refining box ibox
        endif
      enddo
c
c------------------------------------------------
c
c      write(*,*) 'after refinement, nlevels=',nlevels
c      write(*,*) 'nnewboxes=', nnewboxes
c      stop
c
c     old id of the boxes
      allocate(nboxid(mboxes))
      ilevstart(0)=0
      call cumsum(nlevels+1,nblock(0),ilevstart(1))
c
      do j=1,nnewboxes
         jbox=nboxes0+j
         jlev=itree(iptr(2)+jbox-1)
         ilevstart(jlev)=ilevstart(jlev)+1
         nboxid(ilevstart(jlev))=jbox
      enddo
c
c     nlevels updated along the way, retrieve nnewleverls
c     update nboxes, since nnewboxes is avaiable
c
      nnewlevels = nlevels - nlevels0
      nboxes = nboxes0 + nnewboxes
c
c      write(*,*) ' '
c      write(*,*) 'mboxes=',mboxes
c      write(*,*) 'nboxes=',nboxes
c      write(*,*) 'nboxes0=',nboxes0
c      write(*,*) 'nnewboxes=',nnewboxes
c      write(*,*) ' '
c      write(*,*) 'mlevels=',mlevels
c      write(*,*) 'nlevels=',nlevels
c      write(*,*) 'nlevels0=',nlevels0
c      write(*,*) 'nnewlevels=',nnewlevels
c      write(*,*) ' '
c
      if (nnewboxes.eq.0) goto 4600
c
      ifpgh = 1
c     modification fvals -> uvals
      call vol_tree_reorg_laddr(ndim,mboxes,nd,npbox,
     1    nblock,nboxid,nnewboxes,nboxes0,mlevels,
     2    nnewlevels,nlevels0,centers,itree(iptr(1)),
     3    itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),ifpgh,uvals,ucoefs,grad,hess)
c
c
 4600 continue
c
c      write(*,*) 'after reorg, nlevels=', nlevels
c      write(*,*) 'nboxes=', nboxes
c
c
c------------------------------------------------
c     latest version:
c     set function values on the tree      
c     call a modified version of coarsening
c
      call vol_tree_set_vals(nd,ndim,ipoly,
     1     norder,npbox,nboxes,mboxes,nlevels,
     2     mlevels,ltree,itree,iptr,centers,boxsize,
     3     dpars,zpars,ipars,fun,fvals)
c
      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1     iptr,boxsize,norder,fvals,coefsp,umat_nd)
c
c------------------------------------------------
c
c     coarsening part, hopefylly done by
c     vol_tree_coarsen and
c     fgt_vol_tree_reorg_after_coarsen 
c
      ifpgh = 1
      allocate(ifdelete(mboxes))
c
c     pot -> fvals
c     grad and hess unavailable yet
c     go over vol_tree_coarsen
c
c      write(*,*) 'calling coarsening routine:'
      call vol_tree_coarsen_f2(nd,ndim,reps,ipoly,norder,npbox,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    ifpgh,fvals,coefsp,uvals,ucoefs,grad,hess,
     3    nblock,nboxid,ndelboxes,ifdelete)
c
c      write(*,*) 'after calling coarsening:'
c      write(*,*) 'ndelboxes=', ndelboxes
c
      if (ndelboxes.gt.0) then
c        modification: fvals -> uvals
         call fgt_vol_tree_reorg_after_coarsen(ndim,nboxes,nd,npbox,
     1       nblock,nboxid,ndelboxes,ifdelete,
     2       centers,nlevels,itree(iptr(1)),itree(iptr(2)),
     3       itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     4       ifpgh,uvals,ucoefs,grad,hess)
      endif
c
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
      nboxes = nboxes - ndelboxes
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
c------------------------------------------------
c
c     now make it level-restricted again
      call computecoll(ndim,nlevels,nboxes,itree(iptr(1)),
     1    boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2    itree(iptr(5)),iperiod,itree(iptr(6)),itree(iptr(7)))
c
c     call vol_tree_fix_lr_interp
c     fvals -> uvals, coefsp -> ucoefs
      call vol_tree_fix_lr_interp(ndim,nd,ipoly,iperiod,norder,npbox,
     1    ifpgh,uvals,ucoefs,grad,hess,
     2    mboxes,mlevels,centers,boxsize,nboxes,nlevels,
     3    itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),itree(iptr(6)),itree(iptr(7)))
c
c      write(*,*) 'after fixing lr:'
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
      end subroutine




c newtree to oldtree
c for old adaptfmm2d
c newtree-->oldtree--> fmm -->newtree
c reorder???
c
c
c     convert an new tree to the old tree
c
c     input:
c     nlevels - integer
c          number of levels
c     itree - integer(ltree)
c          tree info
c     iptr - integer(8)
c          iptr(1) - laddr
c          iptr(2) - ilevel
c          iptr(3) - iparent
c          iptr(4) - nchild
c          iptr(5) - ichild
c          iptr(6) - ncoll
c          iptr(7) - coll
c          iptr(8) - ltree
c     centers - double precision (2,nboxes)
c          xy coordinates of box centers in the oct tree
c     boxsize - double precision (0:nlevels)
c          size of box at each of the levels
c     ltree - integer
c          length of tree = 2*(nlevels+1)+17*nboxes

c     output parameters:
c     nlev - total number of  levels
c     levelbox - an array determining the level of each box
c     iparentbox - the parent of each box
c     ichildbox - the four children of each box
c     icolbox - the column of each box
c     irowbox - the row of each box
c     nboxes - integer
c          number of boxes
c     nblevel - the total number of boxes per level
c     iboxlev - the array in which the boxes are arranged
c     istartlev - the pointer to where each level
c               begins in the iboxlev array
c     cent0 - center of the root box
c     xsize0 - size of the root box
c     iperiod = 0 : free space
c               1: periodic

c

      subroutine newtree2oldtree2d(nlev,levelbox,iparentbox,
     1    ichildbox,icolbox,irowbox,nboxes,nblevel,
     2    iboxlev,istartlev,cent0,xsize0,iperiod,
     3    ltree,nlevels,itree,iptr,centers,boxsize)

      implicit real *8 (a-h,o-z)

      integer  levelbox(1)
      integer  nlev, nboxes
      integer  icolbox(nboxes), irowbox(nboxes)
      integer  iparentbox(1), ichildbox(4,nboxes)
      integer  nblevel(0:1), iboxlev(nboxes), istartlev(0:1)
      integer  nlevels
      integer iptr(8),ltree
      integer itree(ltree)
      real *8 cent0(2),xsize0
      real *8 centers(2,*),boxsize(0:1)
      integer iboxlevinv(nboxes)
      integer ibox,i
      
      nlev = nlevels
c
      do ilev=0,nlevels
        do i=itree(2*ilev+1),itree(2*ilev+2)
          ibox=i

c          print *, ilev, ibox,nlevels,iptr(2)
         enddo
      enddo

      

      xsize0=boxsize(0)

      cent0(1)=centers(1,1)
      cent0(2)=centers(2,1)


      do ilev=0,nlevels
cccc  istartlev  nblevel
         istartlev(ilev)=itree(iptr(1)+2*ilev)
         nblevel(ilev)=itree(iptr(1)+2*ilev+1)
     1              -itree(iptr(1)+2*ilev)+1
      enddo
      
        
c
      do ilev=0,nlevels
        do i=itree(2*ilev+1),itree(2*ilev+2)
          ibox=i
cccc  iboxlev(i)
          iboxlev(ibox)=i
cccc  levelbox(ibox)
          levelbox(ibox)=itree(iptr(2)+i-1)
cccc          write(*,*)ibox,levelbox(ibox)
cccc  iparent(ibox)
          jbox=itree(iptr(3)+i-1)
c
          if (jbox.gt.0) then
             iparentbox(ibox)=itree(iptr(3)+i-1)
          else
             iparentbox(ibox) = -1
          endif
cccc ichild(ibox)
          do j=1,4
             ichildbox(j,ibox)=-1
          enddo

          
          do j=1,4
c
            ichildbox(1,ibox)=itree(iptr(5)+4*(i-1)+2)
c
            ichildbox(2,ibox)=itree(iptr(5)+4*(i-1)+3)
c
            ichildbox(3,ibox)=itree(iptr(5)+4*(i-1)+1)
c
            ichildbox(4,ibox)=itree(iptr(5)+4*(i-1))
          enddo
ccccc
c          if (itree(iptr(4)+i-1) .eq. ichild)then
c             write(*,*)'iptr(4) right!'
c          else
c             write(*,*)'iptr(4) 404!'
c          endif
ccccc      icolbox - the column of each box
ccccc      irowbox - the row of each box

            bs = boxsize(ilev)
c
            icol=(centers(1,i)+xsize0/2 - cent0(1))/bs +0.5d0
            irow=(centers(2,i) - cent0(2)+xsize0/2)/bs +0.5d0
c
            icolbox(ibox)=icol
            irowbox(ibox)=irow


        enddo
      enddo
      
      return
      end subroutine

c-------------------------------------------------
c     subroutine vol_tree_adap_fv4(cancellation)
c-------------------------------------------------
c
c     a modified version of vol_tree_adap_fun
c
c     given a tree, and four functions resolved on the
c     leaf nodes, and another function (given by a
c     function handle), reifne and coarsen the tree
c     so that both functions are resolved but not
c     over-resolved on the new tree
c
c     the I/O is pretty much the same as
c     subroutine vol_tree_adap_fun
c
c     except for:
c      uvals1:
c      uvals2:
c      uvals3:
c      uvals4:
c             (in/out) a function sampled on the
c     tensor product grid of each leaf node
c
c     on output, the tree is modified, and uvals
c     is interpolated and rearranged to match
c     the new tree
c
c     rmk: 1. after refinement, interpolate to the children
c          2. when coarsening, check the resolution of uvals too
c
c-------------------------------------------------
c
      subroutine vol_tree_adap_fv4(nd,ndim,eps,ipoly,
     1    iptype,norder,npbox,nboxes,mboxes,nlevels,
     2    mlevels,ltree,itree,iptr,centers,boxsize,
     3    rintl,dpars,zpars,ipars,iperiod,fun,
     4    uvals1,uvals2,uvals3,uvals4)
      implicit real*8 (a-h,o-z)
      integer nd,ndim,mboxes,nlevels
      integer iptr(8),ltree,ipoly
      integer itree(ltree),norder,npbox
      integer ipars(1)
      real *8 eps, eta
      real *8 uvals1(nd,npbox,mboxes),uvals2(nd,npbox,mboxes)
      real *8 uvals3(nd,npbox,mboxes),uvals4(nd,npbox,mboxes)
      real *8 fcoefs(nd,npbox,mboxes)
      real *8 centers(ndim,mboxes),boxsize(0:mlevels)
      real *8 dpars(1)
      external fun
      complex *16 zk
      complex *16 zpars(1)
c     allocatables
      real*8, allocatable :: fvals(:,:,:)
c     grid related stuff
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8, allocatable :: rmask(:), umat_nd(:,:,:)
      real *8 rintl(0:mlevels)
      real *8 xy(ndim)
      real *8, allocatable:: coefsp(:,:,:), ucoefs1(:,:,:)
      real *8, allocatable:: ucoefs2(:,:,:),ucoefs3(:,:,:)
      real *8, allocatable:: ucoefs4(:,:,:)
      real *8, allocatable:: polyvc1(:,:,:,:)
      real *8, allocatable:: polyvc2(:,:,:,:)
      real *8, allocatable:: polyvc3(:,:,:,:)
      real *8, allocatable:: polyvc4(:,:,:,:)

c     flags for refinement and coarsening
      integer isgn(ndim,2**ndim)
c      integer nblock(0:nlevels), ilevstart(0:nlevels)
      integer, allocatable:: nblock(:)
      integer, allocatable:: ilevstart(:)
      integer, allocatable:: nboxid(:), ifdelete(:)
      integer, allocatable:: ifrefine(:)
      integer, allocatable:: irefinelev(:), imaxrefinelev(:)
c
      real*8, allocatable:: grad1(:,:,:,:)
      real*8, allocatable:: hess1(:,:,:,:)
c
      real*8, allocatable:: grad2(:,:,:,:)
      real*8, allocatable:: hess2(:,:,:,:)
c
      real*8, allocatable:: grad3(:,:,:,:)
      real*8, allocatable:: hess3(:,:,:,:)
c
      real*8, allocatable:: grad4(:,:,:,:)
      real*8, allocatable:: hess4(:,:,:,:)
c
c      write(*,*) nd, ndim, eps, ipoly, iptype
c      write(*,*) norder, npbox
c      write(*,*) nboxes, mboxes, nlevels, mlevels, ltree
c      write(*,*) dpars(1), iperiod
c
      mc=2**ndim
      mnbors=3**ndim
      npols = norder**ndim
c
      allocate(fvals(nd,npbox,mboxes))
c
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
      allocate(umat_nd(norder,norder,ndim))
      allocate(polyvc1(norder,norder,ndim,mc))
      allocate(polyvc2(norder,norder,ndim,mc))
      allocate(polyvc3(norder,norder,ndim,mc))
      allocate(polyvc4(norder,norder,ndim,mc))
c
      call get_child_box_sign(ndim,isgn)
      call get_p2c_interp_matrices(ndim,ipoly,norder,isgn,polyvc1)
      call get_p2c_interp_matrices(ndim,ipoly,norder,isgn,polyvc2)
      call get_p2c_interp_matrices(ndim,ipoly,norder,isgn,polyvc3)
      call get_p2c_interp_matrices(ndim,ipoly,norder,isgn,polyvc4)
c
      allocate(grid(ndim,npbox))
c
      allocate(grad1(nd,ndim,npbox,1))
      allocate(hess1(nd,ndim*(ndim+1)/2,npbox,1))
c
      allocate(grad2(nd,ndim,npbox,1))
      allocate(hess2(nd,ndim*(ndim+1)/2,npbox,1))
c
      allocate(grad3(nd,ndim,npbox,1))
      allocate(hess3(nd,ndim*(ndim+1)/2,npbox,1))
c
      allocate(grad4(nd,ndim,npbox,1))
      allocate(hess4(nd,ndim*(ndim+1)/2,npbox,1))

c
c     initialize the grid related stuff
c
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)
c
      norder2=norder*norder
      do i=1,ndim
         call dcopy_f77(norder2,umat,1,umat_nd(1,1,i),1)
      enddo
c
      do i=1,norder
        xq(i) = xq(i)/2
      enddo
c
      if (ndim.eq.1) then
         do i=1,norder
            grid(1,i)=xq(i)
         enddo
      elseif (ndim.eq.2) then
         call mesh2d(xq,norder,xq,norder,grid)
      elseif (ndim.eq.3) then
         call mesh3d(xq,norder,xq,norder,xq,norder,grid)
      endif
c
c--------------------------------------------------
c     on the currect tree
c     eval fun on the grid of leaf nodes
      do ilev=0, nlevels
        do ibox = itree(2*ilev+1), itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild .eq. 0) then
            do i=1,npbox
              do j=1,ndim
                xy(j) = grid(j,i)*boxsize(ilev)+centers(j,ibox)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,ibox))
            enddo
          endif
        enddo
      enddo
c      write(*,*)'1111111'
c
c     compute rmask and rsum for estimating whether the function is resolved
      allocate(rmask(npbox))

      call tens_prod_get_rmask(ndim,iptype,norder,npbox,
     1    rmask,rsum)
c
      if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**ndim)
      if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**ndim)
      if(iptype.eq.0) rsc = 1.0d0
c
c     rsc = rsc*rnorm
c     problematic?
      reps = eps*rsc
c
c     save the original number of boxes
ccc   no need to do this now
ccc   since nboxes is an input
      nboxes0 = nboxes
      nlevels0 = nlevels
c      write(*,*) 'nboxes0=', nboxes0
c      write(*,*) 'nboxes=', nboxes
c      write(*,*) 'mboxes=', mboxes
c
c      write(*,*) '***'
c
c     get rid of the max refine levels and stuff
      allocate(ifrefine(mboxes))
c
c     mboxes here means the max of nboxes
      do ibox = 1, mboxes
        ifrefine(ibox)=0
      enddo
c
c     convert potential to its polynomial expansion coefficients
c     (for all the leaf nodes)
c
c     add: ucoefs, transform uvals -> ucoefs
c
      allocate(coefsp(nd,npbox,mboxes))
      allocate(ucoefs1(nd,npbox,mboxes))
      allocate(ucoefs2(nd,npbox,mboxes))
      allocate(ucoefs3(nd,npbox,mboxes))
      allocate(ucoefs4(nd,npbox,mboxes))
      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,fvals,coefsp,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,uvals1,ucoefs1,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,uvals2,ucoefs2,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,uvals3,ucoefs3,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,uvals4,ucoefs4,umat_nd)

c      write(*,*)'222222222'
c
c     now go down the original tree level by level
c     mark boxes that need to be refined
c
      eta=1.0d0
      zk=1.0d0
c
      do ilev=0, nlevels
        sc = boxsize(ilev)/2
        rscale=sc**eta
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
c         loop over boxes of level ilev
          if (itree(iptr(4)+ibox-1).eq.0) then
c           if a leaf box, check the resolution
             call fun_err(nd,npbox,coefsp(1,1,ibox),
     1            rmask,iptype,rscale,erra)
             erra = erra/rsum
c             write(223,*) ibox, ifrefine(ibox), erra,
c     1                    reps*rintl(ilev), rintl(ilev), reps
             if(erra.gt.reps*rintl(ilev)) then
                ifrefine(ibox)=1
             endif
          endif
        enddo
      enddo
c
c     leaf boxes that need to be refined
c     are flaged in ifrefine
c
c-----------------------------------------
c
c      write(*,*) nlevels, mlevels
c
      do ilev = nlevels+1, mlevels
        boxsize(ilev)=boxsize(ilev-1)/2.0d0
c        write(*,*) ilev, boxsize(ilev)
      enddo
c
      allocate(nblock(0:mlevels))
      allocate(ilevstart(0:mlevels+1))
c
      do i=0,mlevels
         nblock(i)=0
      enddo
c
c-----------------------------------------
c
c      write(*,*)'33333'
      nnewboxes=0
      do ibox = 1, mboxes
c       change this loop to a level-by-level loop?
c       if no refinement is done on the finest level, exit
c       loop over all boxes, including newly added ones
c
c        write(*,*) ibox, ifrefine(ibox)
        if(ifrefine(ibox) .eq. 1) then
c          write(123,*) ibox, ifrefine(ibox)
          ilev=itree(iptr(2)+ibox-1)
          bsh = boxsize(ilev)/4
c         half boxsize on the children's level?
          rscale=bsh**eta
c         set up nchild of ibox
          itree(iptr(4)+ibox-1)=mc
c         setup for each child
          do j=1, mc
            nnewboxes = nnewboxes+1
            jbox = nboxes0+nnewboxes
c
c           set up centers for the new children
            do k=1,ndim
              centers(k,jbox) = centers(k,ibox)+isgn(k,j)*bsh
            enddo
c           set up parent for the new children
            itree(iptr(3)+jbox-1)=ibox
c
c           set up nchild for the new children
            itree(iptr(4)+jbox-1)=0
c
c           set up children for the new children
            do l=1,mc
              itree(iptr(5)+mc*(jbox-1)+l-1) = -1
            enddo
c
c           set up children for the parent
            itree(iptr(5)+mc*(ibox-1)+j-1)=jbox
c
c           set up level for the new children
            jlev=ilev+1
            itree(iptr(2)+jbox-1) = jlev
c           count the number of boxes on the level jlev
            nblock(jlev)=nblock(jlev)+1
c            write(*,*) 'nblock', jlev, nblock(jlev)
c           update nlevels along the way
            if(jlev .gt. nlevels) then
              nlevels = jlev
            endif
c
c           now compute fvals for
c           the newly generated boxes
c           clearly not the right grid
            do ii=1, npbox
              do jj=1, ndim
                xy(jj) = grid(jj,ii)*boxsize(jlev)+centers(jj,jbox)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,ii,jbox))
c             okay, need to check the error
            enddo
c
c-------------------
c           add here: interpolate uvals to the new boxes
c
            call ortho_eval_nd(ndim,nd,norder,
     1           ucoefs1(1,1,ibox),norder,uvals1(1,1,jbox),
     2           polyvc1(1,1,1,j))
c
            call ortho_eval_nd(ndim,nd,norder,
     1           ucoefs2(1,1,ibox),norder,uvals2(1,1,jbox),
     2           polyvc2(1,1,1,j))
c
            call ortho_eval_nd(ndim,nd,norder,
     1           ucoefs3(1,1,ibox),norder,uvals3(1,1,jbox),
     2           polyvc3(1,1,1,j))
c
            call ortho_eval_nd(ndim,nd,norder,
     1           ucoefs4(1,1,ibox),norder,uvals4(1,1,jbox),
     2           polyvc4(1,1,1,j))
c            write(*,*)'dddddd'
c
            call ortho_trans_nd(ndim,nd,0,norder,fvals(1,1,jbox),
     1           coefsp(1,1,jbox),umat_nd)

c          add here: convert uvals to ucoefs on jbox
c
            call ortho_trans_nd(ndim,nd,0,norder,uvals1(1,1,jbox),
     1           ucoefs1(1,1,jbox),umat_nd)
c
  
            call ortho_trans_nd(ndim,nd,0,norder,uvals2(1,1,jbox),
     1           ucoefs2(1,1,jbox),umat_nd)
c
 
            call ortho_trans_nd(ndim,nd,0,norder,uvals3(1,1,jbox),
     1           ucoefs3(1,1,jbox),umat_nd)
c
          
            call ortho_trans_nd(ndim,nd,0,norder,uvals4(1,1,jbox),
     1           ucoefs4(1,1,jbox),umat_nd)


c
            call fun_err(nd,npbox,coefsp(1,1,jbox),rmask,
     1            iptype,rscale,erra)
            erra = erra/rsum
c           error obtained for the newly added box jbox
c
            if(erra .gt. reps*rintl(jlev)) then
c             flag the new box
              ifrefine(jbox) = 1
            endif
c            write(124,*) erra, reps*rintl(jlev), ifrefine(jbox)
c            write(*,*) '...', erra,reps, rintl(jlev)
c         end of the loop over mc children
          enddo
c       end refining box ibox
        endif
      enddo


c
c------------------------------------------------
c
c      write(*,*) 'after refinement, nlevels=',nlevels
c      write(*,*) 'nnewboxes=', nnewboxes
c      stop
c
c     old id of the boxes
      allocate(nboxid(mboxes))
      ilevstart(0)=0
      call cumsum(nlevels+1,nblock(0),ilevstart(1))
c
      do j=1,nnewboxes
         jbox=nboxes0+j
         jlev=itree(iptr(2)+jbox-1)
         ilevstart(jlev)=ilevstart(jlev)+1
         nboxid(ilevstart(jlev))=jbox
      enddo
c
c     nlevels updated along the way, retrieve nnewleverls
c     update nboxes, since nnewboxes is avaiable
c
      nnewlevels = nlevels - nlevels0
      nboxes = nboxes0 + nnewboxes
c
c      write(*,*) ' '
c      write(*,*) 'mboxes=',mboxes
c      write(*,*) 'nboxes=',nboxes
c      write(*,*) 'nboxes0=',nboxes0
c      write(*,*) 'nnewboxes=',nnewboxes
c      write(*,*) ' '
c      write(*,*) 'mlevels=',mlevels
c      write(*,*) 'nlevels=',nlevels
c      write(*,*) 'nlevels0=',nlevels0
c      write(*,*) 'nnewlevels=',nnewlevels
c      write(*,*) ' '
c
      if (nnewboxes.eq.0) goto 4600
c
      ifpgh = 1
c     modification fvals -> uvals
c
      call vol_tree_reorg_laddr(ndim,mboxes,nd,npbox,
     1    nblock,nboxid,nnewboxes,nboxes0,mlevels,
     2    nnewlevels,nlevels0,centers,itree(iptr(1)),
     3    itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),ifpgh,uvals1,ucoefs1,grad1,hess1)
c
      call vol_tree_reorg_laddr(ndim,mboxes,nd,npbox,
     1    nblock,nboxid,nnewboxes,nboxes0,mlevels,
     2    nnewlevels,nlevels0,centers,itree(iptr(1)),
     3    itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),ifpgh,uvals2,ucoefs2,grad2,hess2)
c
      call vol_tree_reorg_laddr(ndim,mboxes,nd,npbox,
     1    nblock,nboxid,nnewboxes,nboxes0,mlevels,
     2    nnewlevels,nlevels0,centers,itree(iptr(1)),
     3    itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),ifpgh,uvals3,ucoefs3,grad3,hess3)
c
      call vol_tree_reorg_laddr(ndim,mboxes,nd,npbox,
     1    nblock,nboxid,nnewboxes,nboxes0,mlevels,
     2    nnewlevels,nlevels0,centers,itree(iptr(1)),
     3    itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),ifpgh,uvals4,ucoefs4,grad4,hess4)
c
c
 4600 continue
c
c      write(*,*) 'after reorg, nlevels=', nlevels
c      write(*,*) 'nboxes=', nboxes
c
c
c------------------------------------------------
c     latest version:
c     set function values on the tree
c     call a modified version of coarsening
c
      call vol_tree_set_vals(nd,ndim,ipoly,
     1     norder,npbox,nboxes,mboxes,nlevels,
     2     mlevels,ltree,itree,iptr,centers,boxsize,
     3     dpars,zpars,ipars,fun,fvals)
c
      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1     iptr,boxsize,norder,fvals,coefsp,umat_nd)
c
c------------------------------------------------
c
c     coarsening part, hopefylly done by
c     vol_tree_coarsen and
c     fgt_vol_tree_reorg_after_coarsen
c
      ifpgh = 1
      allocate(ifdelete(mboxes))
c
c     pot -> fvals
c     grad and hess unavailable yet
c     go over vol_tree_coarsen
c
c      write(*,*) 'calling coarsening routine:'
      call vol_tree_coarsen_f2(nd,ndim,reps,ipoly,norder,npbox,
     1    mboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    ifpgh,fvals,coefsp,uvals1,ucoefs1,grad1,hess1,
     3    nblock,nboxid,ndelboxes,ifdelete)

      call vol_tree_coarsen_f2(nd,ndim,reps,ipoly,norder,npbox,
     1    mboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    ifpgh,fvals,coefsp,uvals2,ucoefs2,grad2,hess2,
     3    nblock,nboxid,ndelboxes,ifdelete)

      call vol_tree_coarsen_f2(nd,ndim,reps,ipoly,norder,npbox,
     1    mboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    ifpgh,fvals,coefsp,uvals3,ucoefs3,grad3,hess3,
     3    nblock,nboxid,ndelboxes,ifdelete)

      call vol_tree_coarsen_f2(nd,ndim,reps,ipoly,norder,npbox,
     1    mboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    ifpgh,fvals,coefsp,uvals4,ucoefs4,grad4,hess4,
     3    nblock,nboxid,ndelboxes,ifdelete)
c
c      write(*,*) 'after calling coarsening:'
c      write(*,*) 'ndelboxes=', ndelboxes
c
      if (ndelboxes.gt.0) then
c        modification: fvals -> uvals
         call fgt_vol_tree_reorg_after_coarsen(ndim,mboxes,nd,npbox,
     1       nblock,nboxid,ndelboxes,ifdelete,
     2       centers,nlevels,itree(iptr(1)),itree(iptr(2)),
     3       itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     4       ifpgh,uvals1,ucoefs1,grad1,hess1)

         call fgt_vol_tree_reorg_after_coarsen(ndim,mboxes,nd,npbox,
     1       nblock,nboxid,ndelboxes,ifdelete,
     2       centers,nlevels,itree(iptr(1)),itree(iptr(2)),
     3       itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     4       ifpgh,uvals2,ucoefs2,grad2,hess2)

         call fgt_vol_tree_reorg_after_coarsen(ndim,mboxes,nd,npbox,
     1       nblock,nboxid,ndelboxes,ifdelete,
     2       centers,nlevels,itree(iptr(1)),itree(iptr(2)),
     3       itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     4       ifpgh,uvals3,ucoefs3,grad3,hess3)

         call fgt_vol_tree_reorg_after_coarsen(ndim,mboxes,nd,npbox,
     1       nblock,nboxid,ndelboxes,ifdelete,
     2       centers,nlevels,itree(iptr(1)),itree(iptr(2)),
     3       itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     4       ifpgh,uvals4,ucoefs4,grad4,hess4)
      endif
c
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
      nboxes = nboxes - ndelboxes
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
c------------------------------------------------
c
c     now make it level-restricted again
      call computecoll(ndim,nlevels,nboxes,itree(iptr(1)),
     1    boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2    itree(iptr(5)),iperiod,itree(iptr(6)),itree(iptr(7)))
c
c     call vol_tree_fix_lr_interp
c     fvals -> uvals, coefsp -> ucoefs
      call vol_tree_fix_lr_interp(ndim,nd,ipoly,iperiod,norder,npbox,
     1    ifpgh,uvals1,ucoefs1,grad1,hess1,
     2    mboxes,mlevels,centers,boxsize,nboxes,nlevels,
     3    itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),itree(iptr(6)),itree(iptr(7)))

      call vol_tree_fix_lr_interp(ndim,nd,ipoly,iperiod,norder,npbox,
     1    ifpgh,uvals2,ucoefs2,grad2,hess2,
     2    mboxes,mlevels,centers,boxsize,nboxes,nlevels,
     3    itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),itree(iptr(6)),itree(iptr(7)))

      call vol_tree_fix_lr_interp(ndim,nd,ipoly,iperiod,norder,npbox,
     1    ifpgh,uvals3,ucoefs3,grad3,hess3,
     2    mboxes,mlevels,centers,boxsize,nboxes,nlevels,
     3    itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),itree(iptr(6)),itree(iptr(7)))

      call vol_tree_fix_lr_interp(ndim,nd,ipoly,iperiod,norder,npbox,
     1    ifpgh,uvals4,ucoefs4,grad4,hess4,
     2    mboxes,mlevels,centers,boxsize,nboxes,nlevels,
     3    itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),itree(iptr(6)),itree(iptr(7)))
c
c      write(*,*) 'after fixing lr:'
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
      end subroutine

      
c-------------------------------------------------
c     subroutine vol_tree_adap_fv4ns
c-------------------------------------------------
c
c     a modified version of vol_tree_adap_fun
c
c     given a tree, and a function resolved on the
c     leaf nodes, and another function (given by a
c     function handle), reifne and coarsen the tree
c     so that both functions are resolved but not
c     over-resolved on the new tree
c
c     the I/O is pretty much the same as
c     subroutine vol_tree_adap_fun
c
c     except for:
c     uvals: (in/out) a function sampled on the
c     tensor product grid of each leaf node
c
c     on output, the tree is modified, and uvals
c     is interpolated and rearranged to match
c     the new tree
c
c     rmk: 1. after refinement, interpolate to the children
c          2. when coarsening, check the resolution of uvals too
c
c-------------------------------------------------
c
      subroutine vol_tree_adap_fv4ns(ndim,eps,ipoly,
     1    iptype,norder,npbox,nboxes,mboxes,nlevels,
     2    mlevels,ltree,itree,iptr,centers,boxsize,
     3    rintl,dpars,zpars,ipars,iperiod,fun,uvals1,
     4    uvals2,uvals3,uvals4)
      implicit real*8 (a-h,o-z)
      parameter(nd=8)
      integer ndim,mboxes,nlevels
      integer iptr(8),ltree,ipoly
      integer itree(ltree),norder,npbox
      integer ipars(1)
      real *8 eps, eta
      real *8 uvals1(2,npbox,mboxes)
      real *8 uvals2(2,npbox,mboxes)
      real *8 uvals3(2,npbox,mboxes)
      real *8 uvals4(2,npbox,mboxes)

      real *8 uvals(nd,npbox,mboxes),fcoefs(nd,npbox,mboxes)
      real *8 centers(ndim,mboxes),boxsize(0:mlevels)
      real *8 dpars(1)
      external fun
      complex *16 zk
      complex *16 zpars(1)
c     allocatables
      real*8, allocatable :: fvals(:,:,:)
c     grid related stuff
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8, allocatable :: rmask(:), umat_nd(:,:,:)
      real *8 rintl(0:mlevels)
      real *8 xy(ndim)
      real *8, allocatable:: coefsp(:,:,:), ucoefs(:,:,:)
      real *8, allocatable:: polyvc(:,:,:,:)
c     flags for refinement and coarsening
      integer isgn(ndim,2**ndim)
c      integer nblock(0:nlevels), ilevstart(0:nlevels)
      integer, allocatable:: nblock(:)
      integer, allocatable:: ilevstart(:)
      integer, allocatable:: nboxid(:), ifdelete(:)
      integer, allocatable:: ifrefine(:)
      integer, allocatable:: irefinelev(:), imaxrefinelev(:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
c
c      write(*,*) nd, ndim, eps, ipoly, iptype
c      write(*,*) norder, npbox
c      write(*,*) nboxes, mboxes, nlevels, mlevels, ltree
c      write(*,*) dpars(1), iperiod
c
c  Initialization
      do i = 1,nboxes
         do j =1,norder**ndim
c
            uvals(1,j,i)=uvals1(1,j,i)
c
            uvals(2,j,i)=uvals1(2,j,i)
c
            uvals(3,j,i)=uvals2(1,j,i)
c
            uvals(4,j,i)=uvals2(2,j,i)
c
            uvals(5,j,i)=uvals3(1,j,i)
c
            uvals(6,j,i)=uvals3(2,j,i)
c
            uvals(7,j,i)=uvals4(1,j,i)
c
            uvals(8,j,i)=uvals4(2,j,i)
   
         enddo
      enddo



      mc=2**ndim
      mnbors=3**ndim
      npols = norder**ndim
c
      allocate(fvals(nd,npbox,mboxes))
c
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
      allocate(umat_nd(norder,norder,ndim))
      allocate(polyvc(norder,norder,ndim,mc))
c
      call get_child_box_sign(ndim,isgn)
      call get_p2c_interp_matrices(ndim,ipoly,norder,isgn,polyvc)
c
      allocate(grid(ndim,npbox))
      allocate(grad(nd,ndim,npbox,1))
      allocate(hess(nd,ndim*(ndim+1)/2,npbox,1))
c
c     initialize the grid related stuff
c
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)
c
      norder2=norder*norder
      do i=1,ndim
         call dcopy_f77(norder2,umat,1,umat_nd(1,1,i),1)
      enddo
c
      do i=1,norder
        xq(i) = xq(i)/2
      enddo
c
      if (ndim.eq.1) then
         do i=1,norder
            grid(1,i)=xq(i)
         enddo
      elseif (ndim.eq.2) then
         call mesh2d(xq,norder,xq,norder,grid)
      elseif (ndim.eq.3) then
         call mesh3d(xq,norder,xq,norder,xq,norder,grid)
      endif
c
c--------------------------------------------------
c     on the currect tree
c     eval fun on the grid of leaf nodes
      do ilev=0, nlevels
        do ibox = itree(2*ilev+1), itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild .eq. 0) then
            do i=1,npbox
              do j=1,ndim
                xy(j) = grid(j,i)*boxsize(ilev)+centers(j,ibox)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,ibox))
c              write(*,*)ibox,fvals(3,i,ibox)
            enddo
          endif
        enddo
      enddo

c
c     compute rmask and rsum for estimating whether the function is resolved
      allocate(rmask(npbox))

      call tens_prod_get_rmask(ndim,iptype,norder,npbox,
     1    rmask,rsum)
c
      if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**ndim)
      if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**ndim)
      if(iptype.eq.0) rsc = 1.0d0
c
c     rsc = rsc*rnorm
c     problematic?
      reps = eps*rsc
c
c     save the original number of boxes
ccc   no need to do this now
ccc   since nboxes is an input
      nboxes0 = nboxes
      nlevels0 = nlevels
c      write(*,*) 'nboxes0=', nboxes0
c      write(*,*) 'nboxes=', nboxes
c      write(*,*) 'mboxes=', mboxes
c
c      write(*,*) '***'
c
c     get rid of the max refine levels and stuff
      allocate(ifrefine(mboxes))
c
c     mboxes here means the max of nboxes
      do ibox = 1, mboxes
        ifrefine(ibox)=0
      enddo
c
c     convert potential to its polynomial expansion coefficients
c     (for all the leaf nodes)
c
c     add: ucoefs, transform uvals -> ucoefs
c
      allocate(coefsp(nd,npbox,mboxes))
      allocate(ucoefs(nd,npbox,mboxes))

      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,fvals,coefsp,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,uvals,ucoefs,umat_nd)

c
c     now go down the original tree level by level
c     mark boxes that need to be refined
c
      eta=1.0d0
      zk=1.0d0
c
      do ilev=0, nlevels
        sc = boxsize(ilev)/2
        rscale=sc**eta
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
c         loop over boxes of level ilev
          if (itree(iptr(4)+ibox-1).eq.0) then
c           if a leaf box, check the resolution
             call fun_err(nd,npbox,coefsp(1,1,ibox),
     1            rmask,iptype,rscale,erra)
             erra = erra/rsum
c             write(223,*) ibox, ifrefine(ibox), erra,
c     1                    reps*rintl(ilev), rintl(ilev), reps
             if(erra.gt.reps*rintl(ilev)) then
                ifrefine(ibox)=1
             endif
          endif
        enddo
      enddo

c
c     leaf boxes that need to be refined
c     are flaged in ifrefine
c
c-----------------------------------------
c
c      write(*,*) nlevels, mlevels
c
      do ilev = nlevels+1, mlevels
        boxsize(ilev)=boxsize(ilev-1)/2.0d0
c        write(*,*) ilev, boxsize(ilev)
      enddo
c
      allocate(nblock(0:mlevels))
      allocate(ilevstart(0:mlevels+1))
c
      do i=0,mlevels
         nblock(i)=0
      enddo
c
c-----------------------------------------
c
      nnewboxes=0
      do ibox = 1, mboxes
c       change this loop to a level-by-level loop?
c       if no refinement is done on the finest level, exit
c       loop over all boxes, including newly added ones
c
c        write(*,*) ibox, ifrefine(ibox)
        if(ifrefine(ibox) .eq. 1) then
c          write(123,*) ibox, ifrefine(ibox)
          ilev=itree(iptr(2)+ibox-1)
          bsh = boxsize(ilev)/4
c         half boxsize on the children's level?
          rscale=bsh**eta
c         set up nchild of ibox
          itree(iptr(4)+ibox-1)=mc
c         setup for each child
          do j=1, mc
            nnewboxes = nnewboxes+1
            jbox = nboxes0+nnewboxes
c
c           set up centers for the new children

            do k=1,ndim
              centers(k,jbox) = centers(k,ibox)+isgn(k,j)*bsh
            enddo
c           set up parent for the new children
            itree(iptr(3)+jbox-1)=ibox
c
c           set up nchild for the new children
            itree(iptr(4)+jbox-1)=0
c
c           set up children for the new children

            do l=1,mc
              itree(iptr(5)+mc*(jbox-1)+l-1) = -1
            enddo
c
c           set up children for the parent
            itree(iptr(5)+mc*(ibox-1)+j-1)=jbox
c
c           set up level for the new children
            jlev=ilev+1
            itree(iptr(2)+jbox-1) = jlev
c           count the number of boxes on the level jlev
            nblock(jlev)=nblock(jlev)+1
c            write(*,*) 'nblock', jlev, nblock(jlev)
c           update nlevels along the way
            if(jlev .gt. nlevels) then
              nlevels = jlev
            endif
c
c           now compute fvals for
c           the newly generated boxes
c           clearly not the right grid
            do ii=1, npbox
              do jj=1, ndim
                xy(jj) = grid(jj,ii)*boxsize(jlev)+centers(jj,jbox)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,ii,jbox))
c             okay, need to check the error
            enddo
c
c-------------------
c           add here: interpolate uvals to the new boxes

            call ortho_eval_nd(ndim,nd,norder,
     1           ucoefs(1,1,ibox),norder,uvals(1,1,jbox),
     2           polyvc(1,1,1,j))

c
            call ortho_trans_nd(ndim,nd,0,norder,fvals(1,1,jbox),
     1           coefsp(1,1,jbox),umat_nd)
c          add here: convert uvals to ucoefs on jbox
c
        
            call ortho_trans_nd(ndim,nd,0,norder,uvals(1,1,jbox),
     1           ucoefs(1,1,jbox),umat_nd)
c
   
            call fun_err(nd,npbox,coefsp(1,1,jbox),rmask,
     1            iptype,rscale,erra)
            erra = erra/rsum
c           error obtained for the newly added box jbox
c
            if(erra .gt. reps*rintl(jlev)) then
c             flag the new box
              ifrefine(jbox) = 1
            endif
c            write(124,*) erra, reps*rintl(jlev), ifrefine(jbox)
c            write(*,*) '...', erra,reps, rintl(jlev)
c         end of the loop over mc children
          enddo
      
c       end refining box ibox
        endif
      enddo
  
c
c------------------------------------------------
c
c      write(*,*) 'after refinement, nlevels=',nlevels
c      write(*,*) 'nnewboxes=', nnewboxes
c      stop
c
c     old id of the boxes
      allocate(nboxid(mboxes))
      ilevstart(0)=0
      call cumsum(nlevels+1,nblock(0),ilevstart(1))
c
      do j=1,nnewboxes
         jbox=nboxes0+j
         jlev=itree(iptr(2)+jbox-1)
         ilevstart(jlev)=ilevstart(jlev)+1
         nboxid(ilevstart(jlev))=jbox
      enddo
c
c     nlevels updated along the way, retrieve nnewleverls
c     update nboxes, since nnewboxes is avaiable
c
      nnewlevels = nlevels - nlevels0
      nboxes = nboxes0 + nnewboxes
c
c      write(*,*) ' '
c      write(*,*) 'mboxes=',mboxes
c      write(*,*) 'nboxes=',nboxes
c      write(*,*) 'nboxes0=',nboxes0
c      write(*,*) 'nnewboxes=',nnewboxes
c      write(*,*) ' '
c      write(*,*) 'mlevels=',mlevels
c      write(*,*) 'nlevels=',nlevels
c      write(*,*) 'nlevels0=',nlevels0
c      write(*,*) 'nnewlevels=',nnewlevels
c      write(*,*) ' '
c
      if (nnewboxes.eq.0) goto 4600
c
      ifpgh = 1
c     modification fvals -> uvals
      call vol_tree_reorg_laddr(ndim,mboxes,nd,npbox,
     1    nblock,nboxid,nnewboxes,nboxes0,mlevels,
     2    nnewlevels,nlevels0,centers,itree(iptr(1)),
     3    itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),ifpgh,uvals,ucoefs,grad,hess)
c
c
 4600 continue
c
c      write(*,*) 'after reorg, nlevels=', nlevels
c      write(*,*) 'nboxes=', nboxes
c
c
c------------------------------------------------
c     latest version:
c     set function values on the tree
c     call a modified version of coarsening
c
      call vol_tree_set_vals(nd,ndim,ipoly,
     1     norder,npbox,nboxes,mboxes,nlevels,
     2     mlevels,ltree,itree,iptr,centers,boxsize,
     3     dpars,zpars,ipars,fun,fvals)
c
      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1     iptr,boxsize,norder,fvals,coefsp,umat_nd)
c
c------------------------------------------------
c
c     coarsening part, hopefylly done by
c     vol_tree_coarsen and
c     fgt_vol_tree_reorg_after_coarsen
c
      ifpgh = 1
      allocate(ifdelete(mboxes))
c
c     pot -> fvals
c     grad and hess unavailable yet
c     go over vol_tree_coarsen
c
c      write(*,*) 'calling coarsening routine:'
      call vol_tree_coarsen_f2(nd,ndim,reps,ipoly,norder,npbox,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    ifpgh,fvals,coefsp,uvals,ucoefs,grad,hess,
     3    nblock,nboxid,ndelboxes,ifdelete)
c
c      write(*,*) 'after calling coarsening:'
c      write(*,*) 'ndelboxes=', ndelboxes
c
      if (ndelboxes.gt.0) then
c        modification: fvals -> uvals
         call fgt_vol_tree_reorg_after_coarsen(ndim,nboxes,nd,npbox,
     1       nblock,nboxid,ndelboxes,ifdelete,
     2       centers,nlevels,itree(iptr(1)),itree(iptr(2)),
     3       itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     4       ifpgh,uvals,ucoefs,grad,hess)
      endif
c
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
      nboxes = nboxes - ndelboxes
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
c------------------------------------------------
c
c     now make it level-restricted again
      call computecoll(ndim,nlevels,nboxes,itree(iptr(1)),
     1    boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2    itree(iptr(5)),iperiod,itree(iptr(6)),itree(iptr(7)))
c
c     call vol_tree_fix_lr_interp
c     fvals -> uvals, coefsp -> ucoefs
      call vol_tree_fix_lr_interp(ndim,nd,ipoly,iperiod,norder,npbox,
     1    ifpgh,uvals,ucoefs,grad,hess,
     2    mboxes,mlevels,centers,boxsize,nboxes,nlevels,
     3    itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),itree(iptr(6)),itree(iptr(7)))
c
c      write(*,*) 'after fixing lr:'
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
c  Initialization
      do i = 1,nboxes
         do j =1,norder**ndim
c
            uvals1(1,j,i)=uvals(1,j,i)
c
            uvals1(2,j,i)=uvals(2,j,i)
c
            uvals2(1,j,i)=uvals(3,j,i)
c
            uvals2(2,j,i)=uvals(4,j,i)
c
            uvals3(1,j,i)=uvals(5,j,i)
c
            uvals3(2,j,i)=uvals(6,j,i)
c
            uvals4(1,j,i)=uvals(7,j,i)
c
            uvals4(2,j,i)=uvals(8,j,i)

   
         enddo
      enddo
      end subroutine


c-------------------------------------------------
c     subroutine vol_tree_adap_fv2ns
c-------------------------------------------------
c
c     a modified version of vol_tree_adap_fun
c
c     given a tree, and a function resolved on the
c     leaf nodes, and another function (given by a
c     function handle), reifne and coarsen the tree
c     so that both functions are resolved but not
c     over-resolved on the new tree
c
c     the I/O is pretty much the same as
c     subroutine vol_tree_adap_fun
c
c     except for:
c     uvals: (in/out) a function sampled on the
c     tensor product grid of each leaf node
c
c     on output, the tree is modified, and uvals
c     is interpolated and rearranged to match
c     the new tree
c
c     rmk: 1. after refinement, interpolate to the children
c          2. when coarsening, check the resolution of uvals too
c
c-------------------------------------------------
c
      subroutine vol_tree_adap_fv2ns(ndim,eps,ipoly,
     1    iptype,norder,npbox,nboxes,mboxes,nlevels,
     2    mlevels,ltree,itree,iptr,centers,boxsize,
     3    rintl,dpars,zpars,ipars,iperiod,fun,uvals1,
     4    uvals2)
      implicit real*8 (a-h,o-z)
      parameter(nd=4)
      integer ndim,mboxes,nlevels
      integer iptr(8),ltree,ipoly
      integer itree(ltree),norder,npbox
      integer ipars(1)
      real *8 eps, eta
      real *8 uvals1(2,npbox,mboxes)
      real *8 uvals2(2,npbox,mboxes)

      real *8 uvals(nd,npbox,mboxes),fcoefs(nd,npbox,mboxes)
      real *8 centers(ndim,mboxes),boxsize(0:mlevels)
      real *8 dpars(1)
      external fun
      complex *16 zk
      complex *16 zpars(1)
c     allocatables
      real*8, allocatable :: fvals(:,:,:)
c     grid related stuff
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8, allocatable :: rmask(:), umat_nd(:,:,:)
      real *8 rintl(0:mlevels)
      real *8 xy(ndim)
      real *8, allocatable:: coefsp(:,:,:), ucoefs(:,:,:)
      real *8, allocatable:: polyvc(:,:,:,:)
c     flags for refinement and coarsening
      integer isgn(ndim,2**ndim)
c      integer nblock(0:nlevels), ilevstart(0:nlevels)
      integer, allocatable:: nblock(:)
      integer, allocatable:: ilevstart(:)
      integer, allocatable:: nboxid(:), ifdelete(:)
      integer, allocatable:: ifrefine(:)
      integer, allocatable:: irefinelev(:), imaxrefinelev(:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
c
c      write(*,*) nd, ndim, eps, ipoly, iptype
c      write(*,*) norder, npbox
c      write(*,*) nboxes, mboxes, nlevels, mlevels, ltree
c      write(*,*) dpars(1), iperiod
c
c  Initialization
      do i = 1,nboxes
         do j =1,norder**ndim
c
            uvals(1,j,i)=uvals1(1,j,i)
c
            uvals(2,j,i)=uvals1(2,j,i)
c
            uvals(3,j,i)=uvals2(1,j,i)
c
            uvals(4,j,i)=uvals2(2,j,i)

   
         enddo
      enddo



      mc=2**ndim
      mnbors=3**ndim
      npols = norder**ndim
c
      allocate(fvals(nd,npbox,mboxes))
c
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
      allocate(umat_nd(norder,norder,ndim))
      allocate(polyvc(norder,norder,ndim,mc))
c
      call get_child_box_sign(ndim,isgn)
      call get_p2c_interp_matrices(ndim,ipoly,norder,isgn,polyvc)
c
      allocate(grid(ndim,npbox))
      allocate(grad(nd,ndim,npbox,1))
      allocate(hess(nd,ndim*(ndim+1)/2,npbox,1))
c
c     initialize the grid related stuff
c
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)
c
      norder2=norder*norder
      do i=1,ndim
         call dcopy_f77(norder2,umat,1,umat_nd(1,1,i),1)
      enddo
c
      do i=1,norder
        xq(i) = xq(i)/2
      enddo
c
      if (ndim.eq.1) then
         do i=1,norder
            grid(1,i)=xq(i)
         enddo
      elseif (ndim.eq.2) then
         call mesh2d(xq,norder,xq,norder,grid)
      elseif (ndim.eq.3) then
         call mesh3d(xq,norder,xq,norder,xq,norder,grid)
      endif
c
c--------------------------------------------------
c     on the currect tree
c     eval fun on the grid of leaf nodes
      do ilev=0, nlevels
        do ibox = itree(2*ilev+1), itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild .eq. 0) then
            do i=1,npbox
              do j=1,ndim
                xy(j) = grid(j,i)*boxsize(ilev)+centers(j,ibox)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,ibox))
c              write(*,*)ibox,fvals(3,i,ibox)
            enddo
          endif
        enddo
      enddo

c
c     compute rmask and rsum for estimating whether the function is resolved
      allocate(rmask(npbox))

      call tens_prod_get_rmask(ndim,iptype,norder,npbox,
     1    rmask,rsum)
c
      if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**ndim)
      if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**ndim)
      if(iptype.eq.0) rsc = 1.0d0
c
c     rsc = rsc*rnorm
c     problematic?
      reps = eps*rsc
c
c     save the original number of boxes
ccc   no need to do this now
ccc   since nboxes is an input
      nboxes0 = nboxes
      nlevels0 = nlevels
c      write(*,*) 'nboxes0=', nboxes0
c      write(*,*) 'nboxes=', nboxes
c      write(*,*) 'mboxes=', mboxes
c
c      write(*,*) '***'
c
c     get rid of the max refine levels and stuff
      allocate(ifrefine(mboxes))
c
c     mboxes here means the max of nboxes
      do ibox = 1, mboxes
        ifrefine(ibox)=0
      enddo
c
c     convert potential to its polynomial expansion coefficients
c     (for all the leaf nodes)
c
c     add: ucoefs, transform uvals -> ucoefs
c
      allocate(coefsp(nd,npbox,mboxes))
      allocate(ucoefs(nd,npbox,mboxes))

      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,fvals,coefsp,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,uvals,ucoefs,umat_nd)

c
c     now go down the original tree level by level
c     mark boxes that need to be refined
c
      eta=1.0d0
      zk=1.0d0
c
      do ilev=0, nlevels
        sc = boxsize(ilev)/2
        rscale=sc**eta
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
c         loop over boxes of level ilev
          if (itree(iptr(4)+ibox-1).eq.0) then
c           if a leaf box, check the resolution
             call fun_err(nd,npbox,coefsp(1,1,ibox),
     1            rmask,iptype,rscale,erra)
             erra = erra/rsum
c             write(223,*) ibox, ifrefine(ibox), erra,
c     1                    reps*rintl(ilev), rintl(ilev), reps
             if(erra.gt.reps*rintl(ilev)) then
                ifrefine(ibox)=1
             endif
          endif
        enddo
      enddo

c
c     leaf boxes that need to be refined
c     are flaged in ifrefine
c
c-----------------------------------------
c
c      write(*,*) nlevels, mlevels
c
      do ilev = nlevels+1, mlevels
        boxsize(ilev)=boxsize(ilev-1)/2.0d0
c        write(*,*) ilev, boxsize(ilev)
      enddo
c
      allocate(nblock(0:mlevels))
      allocate(ilevstart(0:mlevels+1))
c
      do i=0,mlevels
         nblock(i)=0
      enddo
c
c-----------------------------------------
c
      nnewboxes=0
      do ibox = 1, mboxes
c       change this loop to a level-by-level loop?
c       if no refinement is done on the finest level, exit
c       loop over all boxes, including newly added ones
c
c        write(*,*) ibox, ifrefine(ibox)
        if(ifrefine(ibox) .eq. 1) then
c          write(123,*) ibox, ifrefine(ibox)
          ilev=itree(iptr(2)+ibox-1)
          bsh = boxsize(ilev)/4
c         half boxsize on the children's level?
          rscale=bsh**eta
c         set up nchild of ibox
          itree(iptr(4)+ibox-1)=mc
c         setup for each child
          do j=1, mc
            nnewboxes = nnewboxes+1
            jbox = nboxes0+nnewboxes
c
c           set up centers for the new children

            do k=1,ndim
              centers(k,jbox) = centers(k,ibox)+isgn(k,j)*bsh
            enddo
c           set up parent for the new children
            itree(iptr(3)+jbox-1)=ibox
c
c           set up nchild for the new children
            itree(iptr(4)+jbox-1)=0
c
c           set up children for the new children

            do l=1,mc
              itree(iptr(5)+mc*(jbox-1)+l-1) = -1
            enddo
c
c           set up children for the parent
            itree(iptr(5)+mc*(ibox-1)+j-1)=jbox
c
c           set up level for the new children
            jlev=ilev+1
            itree(iptr(2)+jbox-1) = jlev
c           count the number of boxes on the level jlev
            nblock(jlev)=nblock(jlev)+1
c            write(*,*) 'nblock', jlev, nblock(jlev)
c           update nlevels along the way
            if(jlev .gt. nlevels) then
              nlevels = jlev
            endif
c
c           now compute fvals for
c           the newly generated boxes
c           clearly not the right grid
            do ii=1, npbox
              do jj=1, ndim
                xy(jj) = grid(jj,ii)*boxsize(jlev)+centers(jj,jbox)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,ii,jbox))
c             okay, need to check the error
            enddo
c
c-------------------
c           add here: interpolate uvals to the new boxes

            call ortho_eval_nd(ndim,nd,norder,
     1           ucoefs(1,1,ibox),norder,uvals(1,1,jbox),
     2           polyvc(1,1,1,j))

c
            call ortho_trans_nd(ndim,nd,0,norder,fvals(1,1,jbox),
     1           coefsp(1,1,jbox),umat_nd)
c          add here: convert uvals to ucoefs on jbox
c
        
            call ortho_trans_nd(ndim,nd,0,norder,uvals(1,1,jbox),
     1           ucoefs(1,1,jbox),umat_nd)
c
   
            call fun_err(nd,npbox,coefsp(1,1,jbox),rmask,
     1            iptype,rscale,erra)
            erra = erra/rsum
c           error obtained for the newly added box jbox
c
            if(erra .gt. reps*rintl(jlev)) then
c             flag the new box
              ifrefine(jbox) = 1
            endif
c            write(124,*) erra, reps*rintl(jlev), ifrefine(jbox)
c            write(*,*) '...', erra,reps, rintl(jlev)
c         end of the loop over mc children
          enddo
      
c       end refining box ibox
        endif
      enddo
  
c
c------------------------------------------------
c
c      write(*,*) 'after refinement, nlevels=',nlevels
c      write(*,*) 'nnewboxes=', nnewboxes
c      stop
c
c     old id of the boxes
      allocate(nboxid(mboxes))
      ilevstart(0)=0
      call cumsum(nlevels+1,nblock(0),ilevstart(1))
c
      do j=1,nnewboxes
         jbox=nboxes0+j
         jlev=itree(iptr(2)+jbox-1)
         ilevstart(jlev)=ilevstart(jlev)+1
         nboxid(ilevstart(jlev))=jbox
      enddo
c
c     nlevels updated along the way, retrieve nnewleverls
c     update nboxes, since nnewboxes is avaiable
c
      nnewlevels = nlevels - nlevels0
      nboxes = nboxes0 + nnewboxes
c
c      write(*,*) ' '
c      write(*,*) 'mboxes=',mboxes
c      write(*,*) 'nboxes=',nboxes
c      write(*,*) 'nboxes0=',nboxes0
c      write(*,*) 'nnewboxes=',nnewboxes
c      write(*,*) ' '
c      write(*,*) 'mlevels=',mlevels
c      write(*,*) 'nlevels=',nlevels
c      write(*,*) 'nlevels0=',nlevels0
c      write(*,*) 'nnewlevels=',nnewlevels
c      write(*,*) ' '
c
      if (nnewboxes.eq.0) goto 4600
c
      ifpgh = 1
c     modification fvals -> uvals
      call vol_tree_reorg_laddr(ndim,mboxes,nd,npbox,
     1    nblock,nboxid,nnewboxes,nboxes0,mlevels,
     2    nnewlevels,nlevels0,centers,itree(iptr(1)),
     3    itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),ifpgh,uvals,ucoefs,grad,hess)
c
c
 4600 continue
c
c      write(*,*) 'after reorg, nlevels=', nlevels
c      write(*,*) 'nboxes=', nboxes
c
c
c------------------------------------------------
c     latest version:
c     set function values on the tree
c     call a modified version of coarsening
c
      call vol_tree_set_vals(nd,ndim,ipoly,
     1     norder,npbox,nboxes,mboxes,nlevels,
     2     mlevels,ltree,itree,iptr,centers,boxsize,
     3     dpars,zpars,ipars,fun,fvals)
c
      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1     iptr,boxsize,norder,fvals,coefsp,umat_nd)
c
c------------------------------------------------
c
c     coarsening part, hopefylly done by
c     vol_tree_coarsen and
c     fgt_vol_tree_reorg_after_coarsen
c
      ifpgh = 1
      allocate(ifdelete(mboxes))
c
c     pot -> fvals
c     grad and hess unavailable yet
c     go over vol_tree_coarsen
c
c      write(*,*) 'calling coarsening routine:'
      call vol_tree_coarsen_f2(nd,ndim,reps,ipoly,norder,npbox,
     1    mboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    ifpgh,fvals,coefsp,uvals,ucoefs,grad,hess,
     3    nblock,nboxid,ndelboxes,ifdelete)
c
c      write(*,*) 'after calling coarsening:'
c      write(*,*) 'ndelboxes=', ndelboxes
c
      if (ndelboxes.gt.0) then
c        modification: fvals -> uvals
         call fgt_vol_tree_reorg_after_coarsen(ndim,mboxes,nd,npbox,
     1       nblock,nboxid,ndelboxes,ifdelete,
     2       centers,nlevels,itree(iptr(1)),itree(iptr(2)),
     3       itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     4       ifpgh,uvals,ucoefs,grad,hess)
      endif
c
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
      nboxes = nboxes - ndelboxes
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
c------------------------------------------------
c
c     now make it level-restricted again
      call computecoll(ndim,nlevels,nboxes,itree(iptr(1)),
     1    boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2    itree(iptr(5)),iperiod,itree(iptr(6)),itree(iptr(7)))
c
c     call vol_tree_fix_lr_interp
c     fvals -> uvals, coefsp -> ucoefs
      call vol_tree_fix_lr_interp(ndim,nd,ipoly,iperiod,norder,npbox,
     1    ifpgh,uvals,ucoefs,grad,hess,
     2    mboxes,mlevels,centers,boxsize,nboxes,nlevels,
     3    itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),itree(iptr(6)),itree(iptr(7)))
c
c      write(*,*) 'after fixing lr:'
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
c  Initialization
      do i = 1,nboxes
         do j =1,norder**ndim
c
            uvals1(1,j,i)=uvals(1,j,i)
c
            uvals1(2,j,i)=uvals(2,j,i)
c
            uvals2(1,j,i)=uvals(3,j,i)
c
            uvals2(2,j,i)=uvals(4,j,i)

   
         enddo
      enddo
      end subroutine


c-------------------------------------------------
c
      subroutine vol_tree_adap_fv4nsnew(ndim,eps,ipoly,
     1    iptype,norder,npbox,nboxes,mboxes,nlevels,
     2    mlevels,ltree,itree,iptr,centers,boxsize,
     3    rintl,dpars,zpars,ipars,iperiod,fun,uvals1,
     4    uvals2,uvals3,uvals4)
      implicit real*8 (a-h,o-z)
      parameter(nd=8)
      integer ndim,mboxes,nlevels
      integer iptr(8),ltree,ipoly
      integer itree(ltree),norder,npbox
      integer ipars(1)
      real *8 eps, eta
c      real*8, allocatable :: uvals1(:,:,:)
c      real*8, allocatable :: uvals2(:,:,:)
c      real*8, allocatable :: uvals3(:,:,:)
c      real*8, allocatable :: uvals4(:,:,:)

      real *8 uvals1(2,npbox,mboxes)
      real *8 uvals2(2,npbox,mboxes)
      real *8 uvals3(2,npbox,mboxes)
      real *8 uvals4(2,npbox,mboxes)
      real*8, allocatable :: uvals(:,:,:)
      real*8, allocatable :: fcoefs(:,:,:)
c      real *8 uvals(nd,npbox,mboxes),fcoefs(nd,npbox,mboxes)
      real *8 centers(ndim,mboxes),boxsize(0:mlevels)
      real *8 dpars(1)
      external fun
      complex *16 zk
      complex *16 zpars(1)
c     allocatables
      real*8, allocatable :: fvals(:,:,:)
c     grid related stuff
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8, allocatable :: rmask(:), umat_nd(:,:,:)
      real *8 rintl(0:mlevels)
      real *8 xy(ndim)
      real *8, allocatable:: coefsp(:,:,:), ucoefs(:,:,:)
      real *8, allocatable:: polyvc(:,:,:,:)
c     flags for refinement and coarsening
      integer isgn(ndim,2**ndim)
c      integer nblock(0:nlevels), ilevstart(0:nlevels)
      integer, allocatable:: nblock(:)
      integer, allocatable:: ilevstart(:)
      integer, allocatable:: nboxid(:), ifdelete(:)
      integer, allocatable:: ifrefine(:)
      integer, allocatable:: irefinelev(:), imaxrefinelev(:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
c
c      write(*,*) nd, ndim, eps, ipoly, iptype
c      write(*,*) norder, npbox
c      write(*,*) nboxes, mboxes, nlevels, mlevels, ltree
c      write(*,*) dpars(1), iperiod
c


      mc=2**ndim
      mnbors=3**ndim
      npols = norder**ndim
c
c      allocate(uvals1(2,npbox,mboxes))
c
c      allocate(uvals2(2,npbox,mboxes))
c
c      allocate(uvals3(2,npbox,mboxes))
cc
c      allocate(uvals4(2,npbox,mboxes))
c
      allocate(uvals(nd,npbox,mboxes))
c
c
      allocate(fcoefs(nd,npbox,mboxes))
c
      allocate(fvals(nd,npbox,mboxes))
c
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
      allocate(umat_nd(norder,norder,ndim))
      allocate(polyvc(norder,norder,ndim,mc))
c
      call get_child_box_sign(ndim,isgn)
      call get_p2c_interp_matrices(ndim,ipoly,norder,isgn,polyvc)
c
      allocate(grid(ndim,npbox))
      allocate(grad(nd,ndim,npbox,1))
      allocate(hess(nd,ndim*(ndim+1)/2,npbox,1))

c  Initialization
      do i = 1,nboxes
         do j =1,norder**ndim
c
            uvals(1,j,i)=uvals1(1,j,i)
c
            uvals(2,j,i)=uvals1(2,j,i)
c
            uvals(3,j,i)=uvals2(1,j,i)
c
            uvals(4,j,i)=uvals2(2,j,i)
c
            uvals(5,j,i)=uvals3(1,j,i)
c
            uvals(6,j,i)=uvals3(2,j,i)
c
            uvals(7,j,i)=uvals4(1,j,i)
c
            uvals(8,j,i)=uvals4(2,j,i)
   
         enddo
      enddo

c
c     initialize the grid related stuff
c
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)
c
      norder2=norder*norder
      do i=1,ndim
         call dcopy_f77(norder2,umat,1,umat_nd(1,1,i),1)
      enddo
c
      do i=1,norder
        xq(i) = xq(i)/2
      enddo
c
      if (ndim.eq.1) then
         do i=1,norder
            grid(1,i)=xq(i)
         enddo
      elseif (ndim.eq.2) then
         call mesh2d(xq,norder,xq,norder,grid)
      elseif (ndim.eq.3) then
         call mesh3d(xq,norder,xq,norder,xq,norder,grid)
      endif
c
c--------------------------------------------------
c     on the currect tree
c     eval fun on the grid of leaf nodes
      do ilev=0, nlevels
        do ibox = itree(2*ilev+1), itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild .eq. 0) then
            do i=1,npbox
              do j=1,ndim
                xy(j) = grid(j,i)*boxsize(ilev)+centers(j,ibox)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,ibox))
            enddo
c          write(*,*)ibox,fvals(1,1,ibox),fvals(2,1,ibox)
          endif
        enddo
      enddo

c
c     compute rmask and rsum for estimating whether the function is resolved
      allocate(rmask(npbox))

      call tens_prod_get_rmask(ndim,iptype,norder,npbox,
     1    rmask,rsum)
c
      if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**ndim)
      if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**ndim)
      if(iptype.eq.0) rsc = 1.0d0
c
c     rsc = rsc*rnorm
c     problematic?
      reps = eps*rsc
c
c     save the original number of boxes
ccc   no need to do this now
ccc   since nboxes is an input
      nboxes0 = nboxes
      nlevels0 = nlevels
c      write(*,*) 'nboxes0=', nboxes0
c      write(*,*) 'nboxes=', nboxes
c      write(*,*) 'mboxes=', mboxes
c
c      write(*,*) '***'
c
c     get rid of the max refine levels and stuff
      allocate(ifrefine(mboxes))
c
c     mboxes here means the max of nboxes
      do ibox = 1, mboxes
        ifrefine(ibox)=0
      enddo
c
c     convert potential to its polynomial expansion coefficients
c     (for all the leaf nodes)
c
c     add: ucoefs, transform uvals -> ucoefs
c
      allocate(coefsp(nd,npbox,mboxes))
      allocate(ucoefs(nd,npbox,mboxes))
     
      do i=1,nboxes
         do j =1,npbox
            do k=1,nd
               coefsp(k,j,i)=0.0d0
               ucoefs(k,j,i)=0.0d0
            enddo
         enddo
      enddo


      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,fvals,coefsp,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,uvals,ucoefs,umat_nd)

c
c     now go down the original tree level by level
c     mark boxes that need to be refined
c
      eta=1.0d0
      zk=1.0d0
c
      nnewbox=0
      do ilev=0, nlevels
        sc = boxsize(ilev)/2
        rscale=sc**eta
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
c         loop over boxes of level ilev
          if (itree(iptr(4)+ibox-1).eq.0) then
c           if a leaf box, check the resolution
             call fun_err(nd,npbox,coefsp(1,1,ibox),
     1            rmask,iptype,rscale,erra)
c           if a leaf box, check the resolution
             call fun_err(nd,npbox,ucoefs(1,1,ibox),
     1            rmask,iptype,rscale,errb)

             erra = erra/rsum
             errb = errb/rsum
c             write(*,*) ibox, ifrefine(ibox), erra,errb,
c     1                    reps*rintl(ilev), rintl(ilev), reps
             if((erra.gt.reps*rintl(ilev)) .or.
     1              (errb.gt.reps*rintl(ilev)) ) then
                ifrefine(ibox)=1
              nnewbox=nnewbox+1
             endif
c           write(*,*)ibox,erra,errb,reps*rintl(ilev),ifrefine(ibox)
          endif
        enddo
      enddo
c      write(*,*)'number of boxes needed refine=',nnewbox
c
c     leaf boxes that need to be refined
c     are flaged in ifrefine
c
c-----------------------------------------
c
c      write(*,*) nlevels, mlevels
c
      do ilev = nlevels+1, mlevels
        boxsize(ilev)=boxsize(ilev-1)/2.0d0
c        write(*,*) ilev, boxsize(ilev)
      enddo
c
      allocate(nblock(0:mlevels))
      allocate(ilevstart(0:mlevels+1))
c
      do i=0,mlevels
         nblock(i)=0
      enddo
c
c-----------------------------------------
c
      nnewboxes=0
      do ibox = 1, mboxes
c       change this loop to a level-by-level loop?
c       if no refinement is done on the finest level, exit
c       loop over all boxes, including newly added ones
c
c        write(*,*) ibox, ifrefine(ibox)
        if(ifrefine(ibox) .eq. 1) then
c          write(123,*) ibox, ifrefine(ibox)
          ilev=itree(iptr(2)+ibox-1)
          bsh = boxsize(ilev)/4
c         half boxsize on the children's level?
          rscale=bsh**eta
c         set up nchild of ibox
          itree(iptr(4)+ibox-1)=mc
c         setup for each child
          do j=1, mc
            nnewboxes = nnewboxes+1
            jbox = nboxes0+nnewboxes
c
c           set up centers for the new children

            do k=1,ndim
              centers(k,jbox) = centers(k,ibox)+isgn(k,j)*bsh
            enddo
c           set up parent for the new children
            itree(iptr(3)+jbox-1)=ibox
c
c           set up nchild for the new children
            itree(iptr(4)+jbox-1)=0
c
c           set up children for the new children

            do l=1,mc
              itree(iptr(5)+mc*(jbox-1)+l-1) = -1
            enddo
c
c           set up children for the parent
            itree(iptr(5)+mc*(ibox-1)+j-1)=jbox
c
c           set up level for the new children
            jlev=ilev+1
            itree(iptr(2)+jbox-1) = jlev
c           count the number of boxes on the level jlev
            nblock(jlev)=nblock(jlev)+1
c            write(*,*) 'nblock', jlev, nblock(jlev)
c           update nlevels along the way
            if(jlev .gt. nlevels) then
              nlevels = jlev
            endif
c
c           now compute fvals for
c           the newly generated boxes
c           clearly not the right grid
            do ii=1, npbox
              do jj=1, ndim
                xy(jj) = grid(jj,ii)*boxsize(jlev)+centers(jj,jbox)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,ii,jbox))
c             okay, need to check the error
            enddo
c
c-------------------
c           add here: interpolate uvals to the new boxes

            call ortho_eval_nd(ndim,nd,norder,
     1           ucoefs(1,1,ibox),norder,uvals(1,1,jbox),
     2           polyvc(1,1,1,j))

c
            call ortho_trans_nd(ndim,nd,0,norder,fvals(1,1,jbox),
     1           coefsp(1,1,jbox),umat_nd)
c          add here: convert uvals to ucoefs on jbox
c
        
            call ortho_trans_nd(ndim,nd,0,norder,uvals(1,1,jbox),
     1           ucoefs(1,1,jbox),umat_nd)
c
   
            call fun_err(nd,npbox,coefsp(1,1,jbox),rmask,
     1            iptype,rscale,erra)

            erra = erra/rsum
c
   
            call fun_err(nd,npbox,ucoefs(1,1,jbox),rmask,
     1            iptype,rscale,errb)

            errb = errb/rsum

c           error obtained for the newly added box jbox
c
            if((erra .gt. reps*rintl(jlev)) .or.
     1       (errb .gt. reps*rintl(jlev)) ) then
c             flag the new box
              ifrefine(jbox) = 1
            endif
c          write(*,*)jbox,erra,errb,reps*rintl(ilev),ifrefine(jbox)
c            write(124,*) erra, reps*rintl(jlev), ifrefine(jbox)
c            write(*,*) '...', erra,reps, rintl(jlev)
c         end of the loop over mc children
          enddo
      
c       end refining box ibox
        endif
      enddo
  
c
c------------------------------------------------
c
c      write(*,*) 'after refinement, nlevels=',nlevels
c      write(*,*) 'nnewboxes=', nnewboxes
c      stop
c
c     old id of the boxes
      allocate(nboxid(mboxes))
      ilevstart(0)=0
      call cumsum(nlevels+1,nblock(0),ilevstart(1))
c
      do j=1,nnewboxes
         jbox=nboxes0+j
         jlev=itree(iptr(2)+jbox-1)
         ilevstart(jlev)=ilevstart(jlev)+1
         nboxid(ilevstart(jlev))=jbox
      enddo
c
c     nlevels updated along the way, retrieve nnewleverls
c     update nboxes, since nnewboxes is avaiable
c
      nnewlevels = nlevels - nlevels0
      nboxes = nboxes0 + nnewboxes
c
c      write(*,*) ' '
c      write(*,*) 'mboxes=',mboxes
c      write(*,*) 'nboxes=',nboxes
c      write(*,*) 'nboxes0=',nboxes0
c      write(*,*) 'nnewboxes=',nnewboxes
c      write(*,*) ' '
c      write(*,*) 'mlevels=',mlevels
c      write(*,*) 'nlevels=',nlevels
c      write(*,*) 'nlevels0=',nlevels0
c      write(*,*) 'nnewlevels=',nnewlevels
c      write(*,*) ' '
c
      if (nnewboxes.eq.0) goto 4600
c
      ifpgh = 1
c     modification fvals -> uvals
      call vol_tree_reorg_laddr(ndim,mboxes,nd,npbox,
     1    nblock,nboxid,nnewboxes,nboxes0,mlevels,
     2    nnewlevels,nlevels0,centers,itree(iptr(1)),
     3    itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),ifpgh,uvals,ucoefs,grad,hess)
c
c
 4600 continue
c
c      write(*,*) 'after reorg, nlevels=', nlevels
c      write(*,*) 'nboxes=', nboxes
c
c
c------------------------------------------------
c     latest version:
c     set function values on the tree
c     call a modified version of coarsening
c
      call vol_tree_set_vals(nd,ndim,ipoly,
     1     norder,npbox,nboxes,mboxes,nlevels,
     2     mlevels,ltree,itree,iptr,centers,boxsize,
     3     dpars,zpars,ipars,fun,fvals)
c
      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1     iptr,boxsize,norder,fvals,coefsp,umat_nd)
c
c------------------------------------------------
c
c     coarsening part, hopefylly done by
c     vol_tree_coarsen and
c     fgt_vol_tree_reorg_after_coarsen
c
      ifpgh = 1
      allocate(ifdelete(mboxes))
c
c     pot -> fvals
c     grad and hess unavailable yet
c     go over vol_tree_coarsen
c
c      write(*,*) 'calling coarsening routine:'
      call vol_tree_coarsen_f2(nd,ndim,reps,ipoly,norder,npbox,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    ifpgh,fvals,coefsp,uvals,ucoefs,grad,hess,
     3    nblock,nboxid,ndelboxes,ifdelete)
c
c      write(*,*) 'after calling coarsening:'
c      write(*,*) 'ndelboxes=', ndelboxes
c
      if (ndelboxes.gt.0) then
c        modification: fvals -> uvals
         call fgt_vol_tree_reorg_after_coarsen(ndim,nboxes,nd,npbox,
     1       nblock,nboxid,ndelboxes,ifdelete,
     2       centers,nlevels,itree(iptr(1)),itree(iptr(2)),
     3       itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     4       ifpgh,uvals,ucoefs,grad,hess)
      endif
c
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
      nboxes = nboxes - ndelboxes
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
c------------------------------------------------
c
c     now make it level-restricted again
      call computecoll(ndim,nlevels,nboxes,itree(iptr(1)),
     1    boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2    itree(iptr(5)),iperiod,itree(iptr(6)),itree(iptr(7)))
c
c     call vol_tree_fix_lr_interp
c     fvals -> uvals, coefsp -> ucoefs
      call vol_tree_fix_lr_interp(ndim,nd,ipoly,iperiod,norder,npbox,
     1    ifpgh,uvals,ucoefs,grad,hess,
     2    mboxes,mlevels,centers,boxsize,nboxes,nlevels,
     3    itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),itree(iptr(6)),itree(iptr(7)))
c
c      write(*,*) 'after fixing lr:'
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
c  Initialization
      do i = 1,nboxes
         do j =1,norder**ndim
c
            uvals1(1,j,i)=uvals(1,j,i)
c
            uvals1(2,j,i)=uvals(2,j,i)
c
            uvals2(1,j,i)=uvals(3,j,i)
c
            uvals2(2,j,i)=uvals(4,j,i)
c
            uvals3(1,j,i)=uvals(5,j,i)
c
            uvals3(2,j,i)=uvals(6,j,i)
c
            uvals4(1,j,i)=uvals(7,j,i)
c
            uvals4(2,j,i)=uvals(8,j,i)

   
         enddo
      enddo
      end subroutine

c-------------------------------------------------
c     subroutine vol_tree_adap_fv
c-------------------------------------------------
c
c     a modified version of vol_tree_adap_fun
c
c     given a tree, and a function resolved on the
c     leaf nodes, and another function (given by a
c     function handle), reifne and coarsen the tree
c     so that both functions are resolved but not
c     over-resolved on the new tree
c
c     the I/O is pretty much the same as
c     subroutine vol_tree_adap_fun
c
c     except for:
c     uvals: (in/out) a function sampled on the
c     tensor product grid of each leaf node
c
c     on output, the tree is modified, and uvals
c     is interpolated and rearranged to match
c     the new tree
c
c     rmk: 1. after refinement, interpolate to the children
c          2. when coarsening, check the resolution of uvals too
c
c-------------------------------------------------
c
      subroutine vol_tree_adap_fvnew(nd,ndim,eps,ipoly,
     1    iptype,norder,npbox,nboxes,mboxes,nlevels,
     2    mlevels,ltree,itree,iptr,centers,boxsize,
     3    rintl,dpars,zpars,ipars,iperiod,fun,uvals)
      implicit real*8 (a-h,o-z)
      integer nd,ndim,mboxes,nlevels
      integer iptr(8),ltree,ipoly
      integer itree(ltree),norder,npbox
      integer ipars(1)
      real *8 eps, eta
      real *8 uvals(nd,npbox,mboxes),fcoefs(nd,npbox,mboxes)
      real *8 centers(ndim,mboxes),boxsize(0:mlevels)
      real *8 dpars(1)
      external fun
      complex *16 zk
      complex *16 zpars(1)
c     allocatables
      real*8, allocatable :: fvals(:,:,:)
c     grid related stuff
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8, allocatable :: rmask(:), umat_nd(:,:,:)
      real *8 rintl(0:mlevels)
      real *8 xy(ndim)
      real *8, allocatable:: coefsp(:,:,:), ucoefs(:,:,:)
      real *8, allocatable:: polyvc(:,:,:,:)
c     flags for refinement and coarsening
      integer isgn(ndim,2**ndim)
c      integer nblock(0:nlevels), ilevstart(0:nlevels)
      integer, allocatable:: nblock(:)
      integer, allocatable:: ilevstart(:)
      integer, allocatable:: nboxid(:), ifdelete(:)
      integer, allocatable:: ifrefine(:)
      integer, allocatable:: irefinelev(:), imaxrefinelev(:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
c
c      write(*,*) nd, ndim, eps, ipoly, iptype
c      write(*,*) norder, npbox
c      write(*,*) nboxes, mboxes, nlevels, mlevels, ltree
c      write(*,*) dpars(1), iperiod
c
      mc=2**ndim
      mnbors=3**ndim
      npols = norder**ndim
c
      allocate(fvals(nd,npbox,mboxes))
c
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
      allocate(umat_nd(norder,norder,ndim))
      allocate(polyvc(norder,norder,ndim,mc))
c
      call get_child_box_sign(ndim,isgn)
      call get_p2c_interp_matrices(ndim,ipoly,norder,isgn,polyvc)
c
      allocate(grid(ndim,npbox))
      allocate(grad(nd,ndim,npbox,1))
      allocate(hess(nd,ndim*(ndim+1)/2,npbox,1))
    
c
c     initialize the grid related stuff
c
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)
c
      norder2=norder*norder
      do i=1,ndim
         call dcopy_f77(norder2,umat,1,umat_nd(1,1,i),1)
      enddo

c
      do i=1,norder
        xq(i) = xq(i)/2
      enddo
c
      if (ndim.eq.1) then
         do i=1,norder
            grid(1,i)=xq(i)
         enddo
      elseif (ndim.eq.2) then
         call mesh2d(xq,norder,xq,norder,grid)
      elseif (ndim.eq.3) then
         call mesh3d(xq,norder,xq,norder,xq,norder,grid)
      endif
 
c
c--------------------------------------------------
c     on the currect tree
c     eval fun on the grid of leaf nodes
      do ilev=0, nlevels
        do ibox = itree(2*ilev+1), itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild .eq. 0) then
            do i=1,npbox
              do j=1,ndim
                xy(j) = grid(j,i)*boxsize(ilev)+centers(j,ibox)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,ibox))
            enddo
          endif
        enddo
      enddo

c
c     compute rmask and rsum for estimating whether the function is resolved
      allocate(rmask(npbox))

      call tens_prod_get_rmask(ndim,iptype,norder,npbox,
     1    rmask,rsum)
c
      if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**ndim)
      if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**ndim)
      if(iptype.eq.0) rsc = 1.0d0
c
c     rsc = rsc*rnorm
c     problematic?
      reps = eps*rsc
c
c     save the original number of boxes
ccc   no need to do this now
ccc   since nboxes is an input
      nboxes0 = nboxes
      nlevels0 = nlevels
c      write(*,*) 'nboxes0=', nboxes0
c      write(*,*) 'nboxes=', nboxes
c      write(*,*) 'mboxes=', mboxes
c
c      write(*,*) '***'
c
c     get rid of the max refine levels and stuff
      allocate(ifrefine(mboxes))
c
c     mboxes here means the max of nboxes
      do ibox = 1, mboxes
        ifrefine(ibox)=0
      enddo
c
c     convert potential to its polynomial expansion coefficients
c     (for all the leaf nodes)
c
c     add: ucoefs, transform uvals -> ucoefs
c
      allocate(coefsp(nd,npbox,mboxes))
      allocate(ucoefs(nd,npbox,mboxes))
      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,fvals,coefsp,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,uvals,ucoefs,umat_nd)
c
c     now go down the original tree level by level
c     mark boxes that need to be refined
c
      eta=1.0d0
      zk=1.0d0
c
      do ilev=0, nlevels
        sc = boxsize(ilev)/2
        rscale=sc**eta
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
c         loop over boxes of level ilev
          if (itree(iptr(4)+ibox-1).eq.0) then
c           if a leaf box, check the resolution
             call fun_err(nd,npbox,coefsp(1,1,ibox),
     1            rmask,iptype,rscale,erra)
             call fun_err(nd,npbox,ucoefs(1,1,ibox),
     1            rmask,iptype,rscale,errb)
             erra = erra/rsum
             errb = errb/rsum
c             write(*,*) ibox, ifrefine(ibox), erra,errb,
c     1                    reps*rintl(ilev), rintl(ilev), reps
             if((erra.gt.reps*rintl(ilev)) .or.
     1       (errb.gt.reps*rintl(ilev)) ) then
                ifrefine(ibox)=1
             endif
          endif
        enddo
      enddo
c
c     leaf boxes that need to be refined
c     are flaged in ifrefine
c
c-----------------------------------------
c
c      write(*,*) nlevels, mlevels
c
      do ilev = nlevels+1, mlevels
        boxsize(ilev)=boxsize(ilev-1)/2.0d0
c        write(*,*) ilev, boxsize(ilev)
      enddo
c
      allocate(nblock(0:mlevels))
      allocate(ilevstart(0:mlevels+1))
c
      do i=0,mlevels
         nblock(i)=0
      enddo
c
c-----------------------------------------
c
      nnewboxes=0
      do ibox = 1, mboxes
c       change this loop to a level-by-level loop?
c       if no refinement is done on the finest level, exit
c       loop over all boxes, including newly added ones
c
c        write(*,*) ibox, ifrefine(ibox)
        if(ifrefine(ibox) .eq. 1) then
c          write(123,*) ibox, ifrefine(ibox)
          ilev=itree(iptr(2)+ibox-1)
          bsh = boxsize(ilev)/4
c         half boxsize on the children's level?
          rscale=bsh**eta
c         set up nchild of ibox
          itree(iptr(4)+ibox-1)=mc
c         setup for each child
          do j=1, mc
            nnewboxes = nnewboxes+1
            jbox = nboxes0+nnewboxes
c
c           set up centers for the new children
            do k=1,ndim
              centers(k,jbox) = centers(k,ibox)+isgn(k,j)*bsh
            enddo
c           set up parent for the new children
            itree(iptr(3)+jbox-1)=ibox
c
c           set up nchild for the new children
            itree(iptr(4)+jbox-1)=0
c
c           set up children for the new children
            do l=1,mc
              itree(iptr(5)+mc*(jbox-1)+l-1) = -1
            enddo
c
c           set up children for the parent
            itree(iptr(5)+mc*(ibox-1)+j-1)=jbox
c
c           set up level for the new children
            jlev=ilev+1
            itree(iptr(2)+jbox-1) = jlev
c           count the number of boxes on the level jlev
            nblock(jlev)=nblock(jlev)+1
c            write(*,*) 'nblock', jlev, nblock(jlev)
c           update nlevels along the way
            if(jlev .gt. nlevels) then
              nlevels = jlev
            endif
c
c           now compute fvals for
c           the newly generated boxes
c           clearly not the right grid
            do ii=1, npbox
              do jj=1, ndim
                xy(jj) = grid(jj,ii)*boxsize(jlev)+centers(jj,jbox)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,ii,jbox))
c             okay, need to check the error
            enddo
c
c-------------------
c           add here: interpolate uvals to the new boxes
            call ortho_eval_nd(ndim,nd,norder,
     1           ucoefs(1,1,ibox),norder,uvals(1,1,jbox),
     2           polyvc(1,1,1,j))
c
            call ortho_trans_nd(ndim,nd,0,norder,fvals(1,1,jbox),
     1           coefsp(1,1,jbox),umat_nd)
c          add here: convert uvals to ucoefs on jbox
c
            call ortho_trans_nd(ndim,nd,0,norder,uvals(1,1,jbox),
     1           ucoefs(1,1,jbox),umat_nd)
c
            call fun_err(nd,npbox,coefsp(1,1,jbox),rmask,
     1            iptype,rscale,erra)
c
            call fun_err(nd,npbox,ucoefs(1,1,jbox),rmask,
     1            iptype,rscale,errb)

            erra = erra/rsum
            errb = errb/rsum
c           error obtained for the newly added box jbox
c
            if((erra .gt. reps*rintl(jlev)) .or .
     1    (errb .gt. reps*rintl(jlev))) then
c             flag the new box
              ifrefine(jbox) = 1
            endif
c            write(124,*) erra, reps*rintl(jlev), ifrefine(jbox)
c            write(*,*) '...', erra,reps, rintl(jlev)
c         end of the loop over mc children
          enddo
c       end refining box ibox
        endif
      enddo
c
c------------------------------------------------
c
c      write(*,*) 'after refinement, nlevels=',nlevels
c      write(*,*) 'nnewboxes=', nnewboxes
c      stop
c
c     old id of the boxes
      allocate(nboxid(mboxes))
      ilevstart(0)=0
      call cumsum(nlevels+1,nblock(0),ilevstart(1))
c
      do j=1,nnewboxes
         jbox=nboxes0+j
         jlev=itree(iptr(2)+jbox-1)
         ilevstart(jlev)=ilevstart(jlev)+1
         nboxid(ilevstart(jlev))=jbox
      enddo
c
c     nlevels updated along the way, retrieve nnewleverls
c     update nboxes, since nnewboxes is avaiable
c
      nnewlevels = nlevels - nlevels0
      nboxes = nboxes0 + nnewboxes
c
c      write(*,*) ' '
c      write(*,*) 'mboxes=',mboxes
c      write(*,*) 'nboxes=',nboxes
c      write(*,*) 'nboxes0=',nboxes0
c      write(*,*) 'nnewboxes=',nnewboxes
c      write(*,*) ' '
c      write(*,*) 'mlevels=',mlevels
c      write(*,*) 'nlevels=',nlevels
c      write(*,*) 'nlevels0=',nlevels0
c      write(*,*) 'nnewlevels=',nnewlevels
c      write(*,*) ' '
c
      if (nnewboxes.eq.0) goto 4600
c
      ifpgh = 1
c     modification fvals -> uvals
      call vol_tree_reorg_laddr(ndim,mboxes,nd,npbox,
     1    nblock,nboxid,nnewboxes,nboxes0,mlevels,
     2    nnewlevels,nlevels0,centers,itree(iptr(1)),
     3    itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),ifpgh,uvals,ucoefs,grad,hess)
c
c
 4600 continue
c
c      write(*,*) 'after reorg, nlevels=', nlevels
c      write(*,*) 'nboxes=', nboxes
c
c
c------------------------------------------------
c     latest version:
c     set function values on the tree
c     call a modified version of coarsening
c
      call vol_tree_set_vals(nd,ndim,ipoly,
     1     norder,npbox,nboxes,mboxes,nlevels,
     2     mlevels,ltree,itree,iptr,centers,boxsize,
     3     dpars,zpars,ipars,fun,fvals)
c
      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1     iptr,boxsize,norder,fvals,coefsp,umat_nd)
c
c------------------------------------------------
c
c     coarsening part, hopefylly done by
c     vol_tree_coarsen and
c     fgt_vol_tree_reorg_after_coarsen
c
      ifpgh = 1
      allocate(ifdelete(mboxes))
c
c     pot -> fvals
c     grad and hess unavailable yet
c     go over vol_tree_coarsen
c
c      write(*,*) 'calling coarsening routine:'
      call vol_tree_coarsen_f2(nd,ndim,reps,ipoly,norder,npbox,
     1    mboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    ifpgh,fvals,coefsp,uvals,ucoefs,grad,hess,
     3    nblock,nboxid,ndelboxes,ifdelete)
c
c      write(*,*) 'after calling coarsening:'
c      write(*,*) 'ndelboxes=', ndelboxes
c
      if (ndelboxes.gt.0) then
c        modification: fvals -> uvals
         call fgt_vol_tree_reorg_after_coarsen(ndim,mboxes,nd,npbox,
     1       nblock,nboxid,ndelboxes,ifdelete,
     2       centers,nlevels,itree(iptr(1)),itree(iptr(2)),
     3       itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     4       ifpgh,uvals,ucoefs,grad,hess)
      endif
c
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
      nboxes = nboxes - ndelboxes
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
c------------------------------------------------
c
c     now make it level-restricted again
      call computecoll(ndim,nlevels,nboxes,itree(iptr(1)),
     1    boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2    itree(iptr(5)),iperiod,itree(iptr(6)),itree(iptr(7)))
c
c     call vol_tree_fix_lr_interp
c     fvals -> uvals, coefsp -> ucoefs
      call vol_tree_fix_lr_interp(ndim,nd,ipoly,iperiod,norder,npbox,
     1    ifpgh,uvals,ucoefs,grad,hess,
     2    mboxes,mlevels,centers,boxsize,nboxes,nlevels,
     3    itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),itree(iptr(6)),itree(iptr(7)))
c
c      write(*,*) 'after fixing lr:'
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
      end subroutine

c-------------------------------------------------
c     subroutine vol_tree_adap_fv2ns
c-------------------------------------------------
c
c     a modified version of vol_tree_adap_fun
c
c     given a tree, and a function resolved on the
c     leaf nodes, and another function (given by a
c     function handle), reifne and coarsen the tree
c     so that both functions are resolved but not
c     over-resolved on the new tree
c
c     the I/O is pretty much the same as
c     subroutine vol_tree_adap_fun
c
c     except for:
c     uvals: (in/out) a function sampled on the
c     tensor product grid of each leaf node
c
c     on output, the tree is modified, and uvals
c     is interpolated and rearranged to match
c     the new tree
c
c     rmk: 1. after refinement, interpolate to the children
c          2. when coarsening, check the resolution of uvals too
c
c-------------------------------------------------
c
      subroutine vol_tree_adap_fv2nsnew(ndim,eps,ipoly,
     1    iptype,norder,npbox,nboxes,mboxes,nlevels,
     2    mlevels,ltree,itree,iptr,centers,boxsize,
     3    rintl,dpars,zpars,ipars,iperiod,fun,uvals1,
     4    uvals2)
      implicit real*8 (a-h,o-z)
      parameter(nd=4)
      integer ndim,mboxes,nlevels
      integer iptr(8),ltree,ipoly
      integer itree(ltree),norder,npbox
      integer ipars(1)
      real *8 eps, eta
      real *8 uvals1(2,npbox,mboxes)
      real *8 uvals2(2,npbox,mboxes)

      real *8 uvals(nd,npbox,mboxes),fcoefs(nd,npbox,mboxes)
      real *8 centers(ndim,mboxes),boxsize(0:mlevels)
      real *8 dpars(1)
      external fun
      complex *16 zk
      complex *16 zpars(1)
c     allocatables
      real*8, allocatable :: fvals(:,:,:)
c     grid related stuff
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8, allocatable :: rmask(:), umat_nd(:,:,:)
      real *8 rintl(0:mlevels)
      real *8 xy(ndim)
      real *8, allocatable:: coefsp(:,:,:), ucoefs(:,:,:)
      real *8, allocatable:: polyvc(:,:,:,:)
c     flags for refinement and coarsening
      integer isgn(ndim,2**ndim)
c      integer nblock(0:nlevels), ilevstart(0:nlevels)
      integer, allocatable:: nblock(:)
      integer, allocatable:: ilevstart(:)
      integer, allocatable:: nboxid(:), ifdelete(:)
      integer, allocatable:: ifrefine(:)
      integer, allocatable:: irefinelev(:), imaxrefinelev(:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
c
c      write(*,*) nd, ndim, eps, ipoly, iptype
c      write(*,*) norder, npbox
c      write(*,*) nboxes, mboxes, nlevels, mlevels, ltree
c      write(*,*) dpars(1), iperiod
c
c  Initialization
      do i = 1,nboxes
         do j =1,norder**ndim
c
            uvals(1,j,i)=uvals1(1,j,i)
c
            uvals(2,j,i)=uvals1(2,j,i)
c
            uvals(3,j,i)=uvals2(1,j,i)
c
            uvals(4,j,i)=uvals2(2,j,i)

   
         enddo
      enddo



      mc=2**ndim
      mnbors=3**ndim
      npols = norder**ndim
c
      allocate(fvals(nd,npbox,mboxes))
c
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
      allocate(umat_nd(norder,norder,ndim))
      allocate(polyvc(norder,norder,ndim,mc))
c
      call get_child_box_sign(ndim,isgn)
      call get_p2c_interp_matrices(ndim,ipoly,norder,isgn,polyvc)
c
      allocate(grid(ndim,npbox))
      allocate(grad(nd,ndim,npbox,1))
      allocate(hess(nd,ndim*(ndim+1)/2,npbox,1))
c
c     initialize the grid related stuff
c
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)
c
      norder2=norder*norder
      do i=1,ndim
         call dcopy_f77(norder2,umat,1,umat_nd(1,1,i),1)
      enddo
c
      do i=1,norder
        xq(i) = xq(i)/2
      enddo
c
      if (ndim.eq.1) then
         do i=1,norder
            grid(1,i)=xq(i)
         enddo
      elseif (ndim.eq.2) then
         call mesh2d(xq,norder,xq,norder,grid)
      elseif (ndim.eq.3) then
         call mesh3d(xq,norder,xq,norder,xq,norder,grid)
      endif
c
c--------------------------------------------------
c     on the currect tree
c     eval fun on the grid of leaf nodes
      do ilev=0, nlevels
        do ibox = itree(2*ilev+1), itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild .eq. 0) then
            do i=1,npbox
              do j=1,ndim
                xy(j) = grid(j,i)*boxsize(ilev)+centers(j,ibox)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,ibox))
c              write(*,*)ibox,fvals(3,i,ibox)
            enddo
          endif
        enddo
      enddo

c
c     compute rmask and rsum for estimating whether the function is resolved
      allocate(rmask(npbox))

      call tens_prod_get_rmask(ndim,iptype,norder,npbox,
     1    rmask,rsum)
c
      if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**ndim)
      if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**ndim)
      if(iptype.eq.0) rsc = 1.0d0
c
c     rsc = rsc*rnorm
c     problematic?
      reps = eps*rsc
c
c     save the original number of boxes
ccc   no need to do this now
ccc   since nboxes is an input
      nboxes0 = nboxes
      nlevels0 = nlevels
c      write(*,*) 'nboxes0=', nboxes0
c      write(*,*) 'nboxes=', nboxes
c      write(*,*) 'mboxes=', mboxes
c
c      write(*,*) '***'
c
c     get rid of the max refine levels and stuff
      allocate(ifrefine(mboxes))
c
c     mboxes here means the max of nboxes
      do ibox = 1, mboxes
        ifrefine(ibox)=0
      enddo
c
c     convert potential to its polynomial expansion coefficients
c     (for all the leaf nodes)
c
c     add: ucoefs, transform uvals -> ucoefs
c
      allocate(coefsp(nd,npbox,mboxes))
      allocate(ucoefs(nd,npbox,mboxes))

      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,fvals,coefsp,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,uvals,ucoefs,umat_nd)

c
c     now go down the original tree level by level
c     mark boxes that need to be refined
c
      eta=1.0d0
      zk=1.0d0
c
      do ilev=0, nlevels
        sc = boxsize(ilev)/2
        rscale=sc**eta
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
c         loop over boxes of level ilev
          if (itree(iptr(4)+ibox-1).eq.0) then
c           if a leaf box, check the resolution
             call fun_err(nd,npbox,coefsp(1,1,ibox),
     1            rmask,iptype,rscale,erra)
c           if a leaf box, check the resolution
             call fun_err(nd,npbox,ucoefs(1,1,ibox),
     1            rmask,iptype,rscale,errb)

             erra = erra/rsum
             errb = errb/rsum
c             write(223,*) ibox, ifrefine(ibox), erra,
c     1                    reps*rintl(ilev), rintl(ilev), reps
             if((erra.gt.reps*rintl(ilev)) .or.
     1     (errb.gt.reps*rintl(ilev)) ) then
                ifrefine(ibox)=1
             endif
          endif
        enddo
      enddo

c
c     leaf boxes that need to be refined
c     are flaged in ifrefine
c
c-----------------------------------------
c
c      write(*,*) nlevels, mlevels
c
      do ilev = nlevels+1, mlevels
        boxsize(ilev)=boxsize(ilev-1)/2.0d0
c        write(*,*) ilev, boxsize(ilev)
      enddo
c
      allocate(nblock(0:mlevels))
      allocate(ilevstart(0:mlevels+1))
c
      do i=0,mlevels
         nblock(i)=0
      enddo
c
c-----------------------------------------
c
      nnewboxes=0
      do ibox = 1, mboxes
c       change this loop to a level-by-level loop?
c       if no refinement is done on the finest level, exit
c       loop over all boxes, including newly added ones
c
c        write(*,*) ibox, ifrefine(ibox)
        if(ifrefine(ibox) .eq. 1) then
c          write(123,*) ibox, ifrefine(ibox)
          ilev=itree(iptr(2)+ibox-1)
          bsh = boxsize(ilev)/4
c         half boxsize on the children's level?
          rscale=bsh**eta
c         set up nchild of ibox
          itree(iptr(4)+ibox-1)=mc
c         setup for each child
          do j=1, mc
            nnewboxes = nnewboxes+1
            jbox = nboxes0+nnewboxes
c
c           set up centers for the new children

            do k=1,ndim
              centers(k,jbox) = centers(k,ibox)+isgn(k,j)*bsh
            enddo
c           set up parent for the new children
            itree(iptr(3)+jbox-1)=ibox
c
c           set up nchild for the new children
            itree(iptr(4)+jbox-1)=0
c
c           set up children for the new children

            do l=1,mc
              itree(iptr(5)+mc*(jbox-1)+l-1) = -1
            enddo
c
c           set up children for the parent
            itree(iptr(5)+mc*(ibox-1)+j-1)=jbox
c
c           set up level for the new children
            jlev=ilev+1
            itree(iptr(2)+jbox-1) = jlev
c           count the number of boxes on the level jlev
            nblock(jlev)=nblock(jlev)+1
c            write(*,*) 'nblock', jlev, nblock(jlev)
c           update nlevels along the way
            if(jlev .gt. nlevels) then
              nlevels = jlev
            endif
c
c           now compute fvals for
c           the newly generated boxes
c           clearly not the right grid
            do ii=1, npbox
              do jj=1, ndim
                xy(jj) = grid(jj,ii)*boxsize(jlev)+centers(jj,jbox)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,ii,jbox))
c             okay, need to check the error
            enddo
c
c-------------------
c           add here: interpolate uvals to the new boxes

            call ortho_eval_nd(ndim,nd,norder,
     1           ucoefs(1,1,ibox),norder,uvals(1,1,jbox),
     2           polyvc(1,1,1,j))

c
            call ortho_trans_nd(ndim,nd,0,norder,fvals(1,1,jbox),
     1           coefsp(1,1,jbox),umat_nd)
c          add here: convert uvals to ucoefs on jbox
c
        
            call ortho_trans_nd(ndim,nd,0,norder,uvals(1,1,jbox),
     1           ucoefs(1,1,jbox),umat_nd)
c
   
            call fun_err(nd,npbox,coefsp(1,1,jbox),rmask,
     1            iptype,rscale,erra)
c
   
            call fun_err(nd,npbox,ucoefs(1,1,jbox),rmask,
     1            iptype,rscale,errb)

            erra = erra/rsum
            errb = errb/rsum
c           error obtained for the newly added box jbox
c
            if((erra .gt. reps*rintl(jlev)) .or.
     1         (errb .gt. reps*rintl(jlev))   ) then
c             flag the new box
              ifrefine(jbox) = 1
            endif
c            write(124,*) erra, reps*rintl(jlev), ifrefine(jbox)
c            write(*,*) '...', erra,reps, rintl(jlev)
c         end of the loop over mc children
          enddo
      
c       end refining box ibox
        endif
      enddo
  
c
c------------------------------------------------
c
c      write(*,*) 'after refinement, nlevels=',nlevels
c      write(*,*) 'nnewboxes=', nnewboxes
c      stop
c
c     old id of the boxes
      allocate(nboxid(mboxes))
      ilevstart(0)=0
      call cumsum(nlevels+1,nblock(0),ilevstart(1))
c
      do j=1,nnewboxes
         jbox=nboxes0+j
         jlev=itree(iptr(2)+jbox-1)
         ilevstart(jlev)=ilevstart(jlev)+1
         nboxid(ilevstart(jlev))=jbox
      enddo
c
c     nlevels updated along the way, retrieve nnewleverls
c     update nboxes, since nnewboxes is avaiable
c
      nnewlevels = nlevels - nlevels0
      nboxes = nboxes0 + nnewboxes
c
c      write(*,*) ' '
c      write(*,*) 'mboxes=',mboxes
c      write(*,*) 'nboxes=',nboxes
c      write(*,*) 'nboxes0=',nboxes0
c      write(*,*) 'nnewboxes=',nnewboxes
c      write(*,*) ' '
c      write(*,*) 'mlevels=',mlevels
c      write(*,*) 'nlevels=',nlevels
c      write(*,*) 'nlevels0=',nlevels0
c      write(*,*) 'nnewlevels=',nnewlevels
c      write(*,*) ' '
c
      if (nnewboxes.eq.0) goto 4600
c
      ifpgh = 1
c     modification fvals -> uvals
      call vol_tree_reorg_laddr(ndim,mboxes,nd,npbox,
     1    nblock,nboxid,nnewboxes,nboxes0,mlevels,
     2    nnewlevels,nlevels0,centers,itree(iptr(1)),
     3    itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),ifpgh,uvals,ucoefs,grad,hess)
c
c
 4600 continue
c
c      write(*,*) 'after reorg, nlevels=', nlevels
c      write(*,*) 'nboxes=', nboxes
c
c
c------------------------------------------------
c     latest version:
c     set function values on the tree
c     call a modified version of coarsening
c
      call vol_tree_set_vals(nd,ndim,ipoly,
     1     norder,npbox,nboxes,mboxes,nlevels,
     2     mlevels,ltree,itree,iptr,centers,boxsize,
     3     dpars,zpars,ipars,fun,fvals)
c
      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1     iptr,boxsize,norder,fvals,coefsp,umat_nd)
c
c------------------------------------------------
c
c     coarsening part, hopefylly done by
c     vol_tree_coarsen and
c     fgt_vol_tree_reorg_after_coarsen
c
      ifpgh = 1
      allocate(ifdelete(mboxes))
c
c     pot -> fvals
c     grad and hess unavailable yet
c     go over vol_tree_coarsen
c
c      write(*,*) 'calling coarsening routine:'
      call vol_tree_coarsen_f2(nd,ndim,reps,ipoly,norder,npbox,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    ifpgh,fvals,coefsp,uvals,ucoefs,grad,hess,
     3    nblock,nboxid,ndelboxes,ifdelete)
c
c      write(*,*) 'after calling coarsening:'
c      write(*,*) 'ndelboxes=', ndelboxes
c
      if (ndelboxes.gt.0) then
c        modification: fvals -> uvals
         call fgt_vol_tree_reorg_after_coarsen(ndim,nboxes,nd,npbox,
     1       nblock,nboxid,ndelboxes,ifdelete,
     2       centers,nlevels,itree(iptr(1)),itree(iptr(2)),
     3       itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     4       ifpgh,uvals,ucoefs,grad,hess)
      endif
c
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
      nboxes = nboxes - ndelboxes
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
c------------------------------------------------
c
c     now make it level-restricted again
      call computecoll(ndim,nlevels,nboxes,itree(iptr(1)),
     1    boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2    itree(iptr(5)),iperiod,itree(iptr(6)),itree(iptr(7)))
c
c     call vol_tree_fix_lr_interp
c     fvals -> uvals, coefsp -> ucoefs
      call vol_tree_fix_lr_interp(ndim,nd,ipoly,iperiod,norder,npbox,
     1    ifpgh,uvals,ucoefs,grad,hess,
     2    mboxes,mlevels,centers,boxsize,nboxes,nlevels,
     3    itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),itree(iptr(6)),itree(iptr(7)))
c
c      write(*,*) 'after fixing lr:'
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
c  Initialization
      do i = 1,nboxes
         do j =1,norder**ndim
c
            uvals1(1,j,i)=uvals(1,j,i)
c
            uvals1(2,j,i)=uvals(2,j,i)
c
            uvals2(1,j,i)=uvals(3,j,i)
c
            uvals2(2,j,i)=uvals(4,j,i)

   
         enddo
      enddo
      end subroutine

c
      subroutine vol_tree_adap_fv_refine(nd,ndim,eps,ipoly,
     1    iptype,norder,npbox,nboxes,mboxes,nlevels,
     2    mlevels,ltree,itree,iptr,centers,boxsize,
     3    rintl,dpars,zpars,ipars,iperiod,fun,uvals)
      implicit real*8 (a-h,o-z)
      integer nd,ndim,mboxes,nlevels
      integer iptr(8),ltree,ipoly
      integer itree(ltree),norder,npbox
      integer ipars(1)
      real *8 eps, eta
      real *8 uvals(nd,npbox,mboxes),fcoefs(nd,npbox,mboxes)
      real *8 centers(ndim,mboxes),boxsize(0:mlevels)
      real *8 dpars(1)
      external fun
      complex *16 zk
      complex *16 zpars(1)
c     allocatables
      real*8, allocatable :: fvals(:,:,:)
c     grid related stuff
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8, allocatable :: rmask(:), umat_nd(:,:,:)
      real *8 rintl(0:mlevels)
      real *8 xy(ndim)
      real *8, allocatable:: coefsp(:,:,:), ucoefs(:,:,:)
      real *8, allocatable:: polyvc(:,:,:,:)
c     flags for refinement and coarsening
      integer isgn(ndim,2**ndim)
c      integer nblock(0:nlevels), ilevstart(0:nlevels)
      integer, allocatable:: nblock(:)
      integer, allocatable:: ilevstart(:)
      integer, allocatable:: nboxid(:), ifdelete(:)
      integer, allocatable:: ifrefine(:)
      integer, allocatable:: irefinelev(:), imaxrefinelev(:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
c
c      write(*,*) nd, ndim, eps, ipoly, iptype
c      write(*,*) norder, npbox
c      write(*,*) nboxes, mboxes, nlevels, mlevels, ltree
c      write(*,*) dpars(1), iperiod
c
      mc=2**ndim
      mnbors=3**ndim
      npols = norder**ndim
c
      allocate(fvals(nd,npbox,mboxes))
c
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
      allocate(umat_nd(norder,norder,ndim))
      allocate(polyvc(norder,norder,ndim,mc))
c
      call get_child_box_sign(ndim,isgn)
      call get_p2c_interp_matrices(ndim,ipoly,norder,isgn,polyvc)
c
      allocate(grid(ndim,npbox))
      allocate(grad(nd,ndim,npbox,1))
      allocate(hess(nd,ndim*(ndim+1)/2,npbox,1))
    
c
c     initialize the grid related stuff
c
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)
c
      norder2=norder*norder
      do i=1,ndim
         call dcopy_f77(norder2,umat,1,umat_nd(1,1,i),1)
      enddo

c
      do i=1,norder
        xq(i) = xq(i)/2
      enddo
c
      if (ndim.eq.1) then
         do i=1,norder
            grid(1,i)=xq(i)
         enddo
      elseif (ndim.eq.2) then
         call mesh2d(xq,norder,xq,norder,grid)
      elseif (ndim.eq.3) then
         call mesh3d(xq,norder,xq,norder,xq,norder,grid)
      endif
 
c
c--------------------------------------------------
c     on the currect tree
c     eval fun on the grid of leaf nodes
      do ilev=0, nlevels
        do ibox = itree(2*ilev+1), itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild .eq. 0) then
            do i=1,npbox
              do j=1,ndim
                xy(j) = grid(j,i)*boxsize(ilev)+centers(j,ibox)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,ibox))
            enddo
          endif
        enddo
      enddo

c
c     compute rmask and rsum for estimating whether the function is resolved
      allocate(rmask(npbox))

      call tens_prod_get_rmask(ndim,iptype,norder,npbox,
     1    rmask,rsum)
c
      if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**ndim)
      if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**ndim)
      if(iptype.eq.0) rsc = 1.0d0
c
c     rsc = rsc*rnorm
c     problematic?
      reps = eps*rsc
c
c     save the original number of boxes
ccc   no need to do this now
ccc   since nboxes is an input
      nboxes0 = nboxes
      nlevels0 = nlevels
c      write(*,*) 'nboxes0=', nboxes0
c      write(*,*) 'nboxes=', nboxes
c      write(*,*) 'mboxes=', mboxes
c
c      write(*,*) '***'
c
c     get rid of the max refine levels and stuff
      allocate(ifrefine(mboxes))
c
c     mboxes here means the max of nboxes
      do ibox = 1, mboxes
        ifrefine(ibox)=0
      enddo
c
c     convert potential to its polynomial expansion coefficients
c     (for all the leaf nodes)
c
c     add: ucoefs, transform uvals -> ucoefs
c
      allocate(coefsp(nd,npbox,mboxes))
      allocate(ucoefs(nd,npbox,mboxes))
      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,fvals,coefsp,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,uvals,ucoefs,umat_nd)
c
c     now go down the original tree level by level
c     mark boxes that need to be refined
c
      eta=1.0d0
      zk=1.0d0
c
      do ilev=0, nlevels
        sc = boxsize(ilev)/2
        rscale=sc**eta
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
c         loop over boxes of level ilev
          if (itree(iptr(4)+ibox-1).eq.0) then
c           if a leaf box, check the resolution
             call fun_err(nd,npbox,coefsp(1,1,ibox),
     1            rmask,iptype,rscale,erra)
             erra = erra/rsum
c             write(223,*) ibox, ifrefine(ibox), erra,
c     1                    reps*rintl(ilev), rintl(ilev), reps
             if(erra.gt.reps*rintl(ilev)) then
                ifrefine(ibox)=1
             endif
          endif
        enddo
      enddo
c
c     leaf boxes that need to be refined
c     are flaged in ifrefine
c
c-----------------------------------------
c
c      write(*,*) nlevels, mlevels
c
      do ilev = nlevels+1, mlevels
        boxsize(ilev)=boxsize(ilev-1)/2.0d0
c        write(*,*) ilev, boxsize(ilev)
      enddo
c
      allocate(nblock(0:mlevels))
      allocate(ilevstart(0:mlevels+1))
c
      do i=0,mlevels
         nblock(i)=0
      enddo
c
c-----------------------------------------
c
      nnewboxes=0
      do ibox = 1, mboxes
c       change this loop to a level-by-level loop?
c       if no refinement is done on the finest level, exit
c       loop over all boxes, including newly added ones
c
c        write(*,*) ibox, ifrefine(ibox)
        if(ifrefine(ibox) .eq. 1) then
c          write(123,*) ibox, ifrefine(ibox)
          ilev=itree(iptr(2)+ibox-1)
          bsh = boxsize(ilev)/4
c         half boxsize on the children's level?
          rscale=bsh**eta
c         set up nchild of ibox
          itree(iptr(4)+ibox-1)=mc
c         setup for each child
          do j=1, mc
            nnewboxes = nnewboxes+1
            jbox = nboxes0+nnewboxes
c
c           set up centers for the new children
            do k=1,ndim
              centers(k,jbox) = centers(k,ibox)+isgn(k,j)*bsh
            enddo
c           set up parent for the new children
            itree(iptr(3)+jbox-1)=ibox
c
c           set up nchild for the new children
            itree(iptr(4)+jbox-1)=0
c
c           set up children for the new children
            do l=1,mc
              itree(iptr(5)+mc*(jbox-1)+l-1) = -1
            enddo
c
c           set up children for the parent
            itree(iptr(5)+mc*(ibox-1)+j-1)=jbox
c
c           set up level for the new children
            jlev=ilev+1
            itree(iptr(2)+jbox-1) = jlev
c           count the number of boxes on the level jlev
            nblock(jlev)=nblock(jlev)+1
c            write(*,*) 'nblock', jlev, nblock(jlev)
c           update nlevels along the way
            if(jlev .gt. nlevels) then
              nlevels = jlev
            endif
c
c           now compute fvals for
c           the newly generated boxes
c           clearly not the right grid
            do ii=1, npbox
              do jj=1, ndim
                xy(jj) = grid(jj,ii)*boxsize(jlev)+centers(jj,jbox)
              enddo
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,ii,jbox))
c             okay, need to check the error
            enddo
c
c-------------------
c           add here: interpolate uvals to the new boxes
            call ortho_eval_nd(ndim,nd,norder,
     1           ucoefs(1,1,ibox),norder,uvals(1,1,jbox),
     2           polyvc(1,1,1,j))
c
            call ortho_trans_nd(ndim,nd,0,norder,fvals(1,1,jbox),
     1           coefsp(1,1,jbox),umat_nd)
c          add here: convert uvals to ucoefs on jbox
c
            call ortho_trans_nd(ndim,nd,0,norder,uvals(1,1,jbox),
     1           ucoefs(1,1,jbox),umat_nd)
c
            call fun_err(nd,npbox,coefsp(1,1,jbox),rmask,
     1            iptype,rscale,erra)
            erra = erra/rsum
c           error obtained for the newly added box jbox
c
            if(erra .gt. reps*rintl(jlev)) then
c             flag the new box
              ifrefine(jbox) = 1
            endif
c            write(124,*) erra, reps*rintl(jlev), ifrefine(jbox)
c            write(*,*) '...', erra,reps, rintl(jlev)
c         end of the loop over mc children
          enddo
c       end refining box ibox
        endif
      enddo
c
c------------------------------------------------
c
c      write(*,*) 'after refinement, nlevels=',nlevels
c      write(*,*) 'nnewboxes=', nnewboxes
c      stop
c
c     old id of the boxes
      allocate(nboxid(mboxes))
      ilevstart(0)=0
      call cumsum(nlevels+1,nblock(0),ilevstart(1))
c
      do j=1,nnewboxes
         jbox=nboxes0+j
         jlev=itree(iptr(2)+jbox-1)
         ilevstart(jlev)=ilevstart(jlev)+1
         nboxid(ilevstart(jlev))=jbox
      enddo
c
c     nlevels updated along the way, retrieve nnewleverls
c     update nboxes, since nnewboxes is avaiable
c
      nnewlevels = nlevels - nlevels0
      nboxes = nboxes0 + nnewboxes
c
c      write(*,*) ' '
c      write(*,*) 'mboxes=',mboxes
c      write(*,*) 'nboxes=',nboxes
c      write(*,*) 'nboxes0=',nboxes0
c      write(*,*) 'nnewboxes=',nnewboxes
c      write(*,*) ' '
c      write(*,*) 'mlevels=',mlevels
c      write(*,*) 'nlevels=',nlevels
c      write(*,*) 'nlevels0=',nlevels0
c      write(*,*) 'nnewlevels=',nnewlevels
c      write(*,*) ' '
c
      if (nnewboxes.eq.0) goto 4600
c
      ifpgh = 1
c     modification fvals -> uvals
      call vol_tree_reorg_laddr(ndim,mboxes,nd,npbox,
     1    nblock,nboxid,nnewboxes,nboxes0,mlevels,
     2    nnewlevels,nlevels0,centers,itree(iptr(1)),
     3    itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),ifpgh,uvals,ucoefs,grad,hess)
c
c
 4600 continue

c------------------------------------------------
c
c     now make it level-restricted again
      call computecoll(ndim,nlevels,nboxes,itree(iptr(1)),
     1    boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2    itree(iptr(5)),iperiod,itree(iptr(6)),itree(iptr(7)))
c
c     call vol_tree_fix_lr_interp
c     fvals -> uvals, coefsp -> ucoefs
      call vol_tree_fix_lr_interp(ndim,nd,ipoly,iperiod,norder,npbox,
     1    ifpgh,uvals,ucoefs,grad,hess,
     2    mboxes,mlevels,centers,boxsize,nboxes,nlevels,
     3    itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),itree(iptr(6)),itree(iptr(7)))
c
c      write(*,*) 'after fixing lr:'
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
      end subroutine
