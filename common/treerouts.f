c**********************************************************
c
c     Quadtree related subroutines
c
c**********************************************************
c
c     Primary routines:
c     settree: initializes the tree structure
c     unitree: generates a uniform tree
c     mktreef: generates a tree to resolve a given function
c              to a given accuracy
c     mktreep: given a set of points, this routine 
c              generates a tree so that each leaf box 
c              contains no more than a given number of pts
c              (points are reordered on output)
c     refinetree: given a quad-tree and a set of points,
c                 refine the tree until each leaf box
c                 contains no more than a given number 
c                 of points. sort the point along the way.
c     spreadtree1: given a quad-tree, the center and side 
c                 length of the root box, spread out the 
c                 tree by one level, so that it covers a
c                 bigger region (intermediate routine, 
c                 used by spreadtree only) 
c     spreadtree: given a quad-tree, the center and side 
c                 length of the root box, spread out the 
c                 tree by one or two levels, so that it covers a
c                 bigger region 
c     restriction: checks if the tree is level restricted
c     fixtree: adapts the tree into a level restricted one
c     fixtreenf: same as fixtree, except that this version
c                deals with array of function values at the
c                same time (after fixing, interpolate to 
c                newly generated boxes)
c     fixtreenf2: same as fixtreef, except this maintains 
c                 two arrays of function values
c     mkcolls: generates colleagues 
c              (including the periodic case)
c     mkchild: gets children of each box
c     mkcoeffs: get the coefficients of chebyshev poly 
c               approx of a given function
c     treesort: sort points into leaf boxes of a
c               given tree. on output, the array
c               of points is overwritten so that
c               points of the same box are adjacent
c     printtree: print out the tree
c     tree2targs: assemble grid points on the leaf nodes
c                 of the tree as target points, return
c                 the array of points, as well as maps
c                 between the two kinds of indices
c
c------------------------
c     Adaptive mesh related routines:
c
c     movemesh: given a set of points, a tree, and a function
c               resolved on the old tree, adjust the tree
c               so that each leaf box contains no more than a
c               given number of points. Interpolate the function
c               and sort the points onto the new tree.
c               (calling adaptree)
c
c     subdivide: divides a leaf box into four children
c     subdivide1: same as subdivide, except this version
c                 doesn't consider colleague boxes
c     coarsen1: merges children (that are leaf boxes) of
c             a given parent box (which introduces
c             redudant box numbers, should be followed by
c             subroutine reorder before using the tree)
c     reorder: to clean up the mess caused by merging
c              boxes (all the box numbers might be
c              reassigned)
c     reorderf: reorder the array fval, according to the
c               reordering of the tree
c     adaptree: (the main routine that does the
c              adaptive tree refinement and coarsening,
c              by calling subdivide1, coarsen1, reorder, and
c              chebyk2p, chebyp2k)
c     adaptreef: another version of adaptive mesh refinement
c                and coarsening, which is forcing term driven
c     (RMK: colleagues are always regenerated after tree
c      adjustment, don't do it inside refinement or
c      coarsening)
c
c------------------------
c     Utility routines:
c     posbox: returns the coords of lower left corner
c             and side length in x and y directions
c     mkgrid: generates the tensor product chebyshev grid on 
c             a given box of the tree
c     setf: evaluates a given function on the chebyshev
c           grid on each leaf node
c     getcoeff: 2D Chebyshev transform (move to chebex2d.f)
c     geterror: returns an estimated error for a given 
c               Chebyshev series (sum of tails)
c     sortboxes: sorts boxes into one long array in 
c                ascending order of level
c     chebyk2p: given function values on the k*k grids on 
c               four children boxes, interpolate to get the
c               chebyshev coefficients on the parent box
c     chebyp2k: given function values on the k*k grid on 
c               parent box, interpolate to get the function
c               values on the k*k grids on children boxes
c     chebyvol2p: interpolate a function from the chebyshev
c                 grid on the leaf nodes of a tree to a
c                 given set of points (to be added)
c     outputpot: print out the function values of a function
c                sampled on the n-by-n grid of leaf nodes
c                of a given tree
c     interppot: interpolate the function sampled on 
c                the leaf nodes at a regular tensor 
c                product grid 
c     interpv2p: interpolate the function sampled on the leaf
c                node of a given tree to a given array of points
c     rvecold2new: reorder points from old order to new order
c              given permutation ipold that maps new indices 
c              to the old ones
c     rvecnew2old: reorder points from new order to old order
c              given permutation ipold that maps new indices 
c              to the old ones
c
c------------------
c
c     get_cent_box: returns the center of a given box in the
c                   tree (assuming that the unit box is centered
c                   at the origin)
c     get_rel_coords: given a tree and an array of points sorted
c                     in the tree (leaf nodes), return the relative
c                     coordinates w.r.t. box centers
c     get_dist_small: given a box (with a certain side length) in
c                    the tree, compute the distance from its center
c                    to the center of a colleague box's children
c     get_dist_col: given a box (with a certain side length) in 
c                   the tree, compute the distance from its center
c                   to the center of a colleague box
c
c*****************************************************
c
c     the following subroutine is just set up to initialize the
c     tree structure.  for each box, the level, parent, row, column,
c     and children must be specified.  if any of these things are 
c     not present, then that quantity is set to -1.
c     (this is used as a flag in later routines.)
c     nothing is set upon input to this routine.  all of the above
c     arrays are just set in the subroutine below.
c
c
c     input:
c     
c     nothing is defined on input as this is an initialization routine.
c
c
c     output:
c
c     levelbox is an array defining the level of each box
c
c     icolbox denotes the column of each box
c
c     irowbox denotes the row of each box
c
c     nboxes is the total number of boxes
c
c     nlev is the finest level
c
c**********************************************************
c
      subroutine settree(levelbox,icolbox,irowbox,nboxes,nlev)
      implicit none
c-----global variables
      integer *4 levelbox(1)
      integer *4 icolbox(1), irowbox(1)
      integer *4 nboxes, nlev

c-----local variables
      integer *4 ibox

c     check to see if we are in the adaptive
c     or the uniform case:

ccc      goto 400
c     mock tree number one:
      ibox = 1
      levelbox(ibox) = 0
      icolbox(ibox) = 1
      irowbox(ibox) = 1

      ibox = 2
      levelbox(ibox) = 1
      icolbox(ibox) = 1
      irowbox(ibox) = 1

      ibox = 3
      levelbox(ibox) = 1
      icolbox(ibox) = 2
      irowbox(ibox) = 1

      ibox = 4
      levelbox(ibox) = 1
      icolbox(ibox) = 1
      irowbox(ibox) = 2

      ibox = 5
      levelbox(ibox) = 1
      icolbox(ibox) = 2
      irowbox(ibox) = 2

      ibox = 6
      levelbox(ibox) = 2
      icolbox(ibox) = 1
      irowbox(ibox) = 1

      ibox = 7
      levelbox(ibox) = 2
      icolbox(ibox) = 2
      irowbox(ibox) = 1

      ibox = 8
      levelbox(ibox) = 2
      icolbox(ibox) = 1
      irowbox(ibox) = 2

      ibox = 9
      levelbox(ibox) = 2
      icolbox(ibox) = 2
      irowbox(ibox) = 2

      nboxes = 9
      nlev = 2
      return

300   continue
      ibox = 1
      levelbox(ibox) = 0
      icolbox(ibox) = 1
      irowbox(ibox) = 1
 
      ibox = 2
      levelbox(ibox) = 1
      icolbox(ibox) = 1
      irowbox(ibox) = 1

      ibox = 3
      levelbox(ibox) = 1
      icolbox(ibox) = 2
      irowbox(ibox) = 1
 
      ibox = 4
      levelbox(ibox) = 1
      icolbox(ibox) = 1
      irowbox(ibox) = 2
 
      ibox = 5
      levelbox(ibox) = 1
      icolbox(ibox) = 2
      irowbox(ibox) = 2

      ibox = 6
      levelbox(ibox) = 2
      icolbox(ibox) = 1
      irowbox(ibox) = 1

      ibox = 7
      levelbox(ibox) = 2
      icolbox(ibox) = 2
      irowbox(ibox) = 1
 
      ibox = 8
      levelbox(ibox) = 2
      icolbox(ibox) = 3
      irowbox(ibox) = 1

      ibox = 9
      levelbox(ibox) = 2
      icolbox(ibox) = 4
      irowbox(ibox) = 1

      nboxes = 9
      nlev = 2
      return

400   continue
      ibox = 1
      levelbox(ibox) = 0
      icolbox(ibox) = 1
      irowbox(ibox) = 1

      ibox = 2
      levelbox(ibox) = 1
      icolbox(ibox) = 1
      irowbox(ibox) = 1

      ibox = 3
      levelbox(ibox) = 1
      icolbox(ibox) = 2
      irowbox(ibox) = 1

      ibox = 4
      levelbox(ibox) = 1
      icolbox(ibox) = 1
      irowbox(ibox) = 2

      ibox = 5
      levelbox(ibox) = 1
      icolbox(ibox) = 2
      irowbox(ibox) = 2

      ibox = 6
      levelbox(ibox) = 2
      icolbox(ibox) = 1
      irowbox(ibox) = 1

      ibox = 7
      levelbox(ibox) = 2
      icolbox(ibox) = 2
      irowbox(ibox) = 1

      ibox = 8
      levelbox(ibox) = 2
      icolbox(ibox) = 1
      irowbox(ibox) = 2

      ibox = 9
      levelbox(ibox) = 2
      icolbox(ibox) = 2
      irowbox(ibox) = 2

      ibox = 10
      levelbox(ibox) = 2
      icolbox(ibox) = 1
      irowbox(ibox) = 3

      ibox = 11
      levelbox(ibox) = 2
      icolbox(ibox) = 2
      irowbox(ibox) = 3

      ibox = 12
      levelbox(ibox) = 2
      icolbox(ibox) = 1
      irowbox(ibox) = 4

      ibox = 13
      levelbox(ibox) = 2
      icolbox(ibox) = 2
      irowbox(ibox) = 4

      ibox = 14
      levelbox(ibox) = 2
      icolbox(ibox) = 3
      irowbox(ibox) = 1

      ibox = 15
      levelbox(ibox) = 2
      icolbox(ibox) = 4
      irowbox(ibox) = 1

      ibox = 16
      levelbox(ibox) = 2
      icolbox(ibox) = 3
      irowbox(ibox) = 2

      ibox = 17
      levelbox(ibox) = 2
      icolbox(ibox) = 4
      irowbox(ibox) = 2

      ibox = 18
      levelbox(ibox) = 2
      icolbox(ibox) = 3
      irowbox(ibox) = 3

      ibox = 19
      levelbox(ibox) = 2
      icolbox(ibox) = 4
      irowbox(ibox) = 3

      ibox = 20
      levelbox(ibox) = 2
      icolbox(ibox) = 3
      irowbox(ibox) = 4

      ibox = 21
      levelbox(ibox) = 2
      icolbox(ibox) = 4
      irowbox(ibox) = 4

      ibox = 22
      levelbox(ibox) = 3
      icolbox(ibox) = 5
      irowbox(ibox) = 7

      ibox = 23
      levelbox(ibox) = 3
      icolbox(ibox) = 5
      irowbox(ibox) = 8

      ibox = 24
      levelbox(ibox) = 3
      icolbox(ibox) = 6
      irowbox(ibox) = 7

      ibox = 25
      levelbox(ibox) = 3
      icolbox(ibox) = 6
      irowbox(ibox) = 8

      ibox = 26
      levelbox(ibox) = 4
      icolbox(ibox) = 9
      irowbox(ibox) = 15

      ibox = 27
      levelbox(ibox) = 4
      icolbox(ibox) = 10
      irowbox(ibox) = 15

      ibox = 28
      levelbox(ibox) = 4
      icolbox(ibox) = 9
      irowbox(ibox) = 16

      ibox = 29
      levelbox(ibox) = 4
      icolbox(ibox) = 10
      irowbox(ibox) = 16


      nboxes = 29
      nlev = 4
      return
      end subroutine





c**********************************************************
c
c     the following subroutine is used to generate a uniform tree.
c     this is only for testing purposes.
c
c
c     input:
c
c     nlev is the finest level
c     ladder: an array of dim (0:max) where max>=nlev
c
c     output:
c
c     levelbox is an array determining the level of each box
c
c     icolbox denotes the column of each box
c
c     irowbox denotes the row of each box
c
c     nboxes is the total number of boxes
c
c
c**********************************************************
c
      subroutine unitree(levelbox,icolbox,irowbox,nboxes,nlev, 
     1           iparentbox, ichildbox, nblevel, iboxlev,
     2           istartlev, ladder)
      implicit none
c-----global variables
      integer *4  nboxes, nlev
      integer *4  levelbox(1)
      integer *4  icolbox(1), irowbox(1)
      integer *4  ichildbox(4,1), iparentbox(1)
      integer *4  nblevel(0:nlev)
      integer *4  iboxlev(1),istartlev(0:nlev)
      integer *4  ladder(0:1)
c-----local variables
      integer *4  i, ibox, j, k, l, nside
      integer *4  icol, irow
      integer *4  ic1, ic2, ic3, ic4
      integer *4  ilength,  ilev, istart, idim
      integer *4  mybox, istartc

c     just create a uniform tree, given nlev as the input:
c     (in this loop, let's initialize the parents
c     and children arrays to be -1.)
      ibox = 1
      nside = 1
      do i = 0, nlev
        do j = 1, nside
          do k = 1, nside
            levelbox(ibox) = i
            icolbox(ibox)  = j
            irowbox(ibox)  = k
            iparentbox(ibox) = -1
            do l = 1, 4
              ichildbox(l,ibox) = -1
            end do
            ibox = ibox + 1
          end do
        end do
       nside = 2 * nside
      end do
      nboxes = ibox - 1

c     now initialize the ladder array:
      j = 1
      ladder(0) = 0
      ladder(1) = 1
      do i = 2, nlev + 1
        j = 4*j
        ladder(i) = ladder(i-1) +  j
      end do


c     now set all of the parents and children:
      do i = 0, nlev - 1
        istart = ladder(nlev-i-1)
        ilength = ladder(nlev-i) - istart
        ilev = nlev-i-1
        istartc = ladder(nlev-i)
        idim = 2**ilev
        nside = 2*idim

        do j = 1, ilength
          icol = 1 + mod(j-1,idim )
          irow = 1 + (j-1)/idim

          mybox = istart+j
          ic1 = istartc + (2*icol) + (2*irow-2)*nside
          ic2 = istartc + (2*icol) + (2*irow-1)*nside
          ic3 = istartc + (2*icol-1) + (2*irow-1)*nside
          ic4 = istartc + (2*icol-1) + (2*irow-2)*nside

          ichildbox(1,mybox) = ic1
          ichildbox(2,mybox) = ic2
          ichildbox(3,mybox) = ic3
          ichildbox(4,mybox) = ic4

          iparentbox(ic1) = mybox
          iparentbox(ic2) = mybox
          iparentbox(ic3) = mybox
          iparentbox(ic4) = mybox

        end do
      end do

      nboxes=(4**(nlev+1)-1)/3

      do l=0,nlev
        nblevel(l)=4**l
        istartlev(l)=(4**l-1)/3+1
        do i=istartlev(l),istartlev(l)+nblevel(l)-1
          iboxlev(i)=i
        enddo
      enddo

      return
      end subroutine




c**********************************************************
c
c     the following subroutine is used to generate a tree given a
c     right hand side, fright.  
c     the algorithm works by testing the function values vs. an approximating 
c     polynomial and dividing if necessary.  it is important to note that
c     the tree generated by this algorithm may not satisfy the level 
c     restriction.  it may be necessary to call the routine fixtree 
c     after this to make sure that the tree can be handled properly by the
c     adapfmm6  and boundfmm6 routines.
c
c
c     input:
c
c     maxboxes denotes the maximum number of boxes allowed
c
c     maxlevel denotes the deepest level allowed
c
c     itemparray is just a dummy array
c
c     eps denotes the error bound that determines 
c         when refinement take place
c
c     h is the real function that is the right hand side
c       of the poisson equation
c
c     cent0: center of the region covered
c     xsize0: side length of the unit box 
c     ndeg: spatial order
c
c     output:
c
c     istartlev is the pointer to where each level
c               begins in the iboxlev array
c
c     levelbox is an array determining the level of each box
c
c     nboxes is the total number of boxes
c
c     nlev is the finest level
c
c     icolbox denotes the column of each box
c
c     irowbox denotes the row of each box
c
c     iparentbox denotes the parent of each box
c
c     ichildbox denotes the four children of each box
c
c     nblevel is the total number of boxes per level
c
c     iboxlev is the array in which the boxes are arranged
c
c**********************************************************
c
      subroutine mktreef(levelbox, icolbox, irowbox, 
     1           nboxes, nlev, iparentbox, ichildbox,
     2           nblevel, iboxlev, istartlev,
     3           maxboxes, itemparray, maxlevel, eps, h,
     4           ndeg, cent0, xsize0)
c-----global variables
      implicit none
      integer *4  levelbox(1), maxboxes
      integer *4  nlev, nboxes,  maxlevel
      integer *4  icolbox(1), irowbox(1)
      integer *4  iparentbox(1), ichildbox(4,1)
      integer *4  nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer *4  itemparray(1), ndeg
      real *8  eps, h
      real *8 cent0(2), xsize0
c-----local variables
      integer *4  i, ibox, iflag
      integer *4  j, istart, iend, jj
      integer *4  levflag, ii
      real *8  coefftemp(ndeg,ndeg)
      real *8  epscaled
      real *8  error, hh
      real *8  xf(ndeg), yf(ndeg)
      real *8  ftemp(ndeg,ndeg)
      real *8  wsave2(1000)

      do i = 0, maxlevel
        nblevel(i) = 0
        istartlev(i) = 0
      end do
      do i = 1, maxboxes
        iboxlev(i) = 0
      end do
c
c     first set the big parent box to the 
c     appropriate settings:
c     (initially, there is just one box and
c     it is the big parent box at level 0)
c
      ibox = 1
      levelbox(ibox) = 0
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = -1
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1

      nboxes = 1
      nlev = 0

c     we also need to initialize the adaptive 'ladder'
c     structures to the correct initial values:
      nblevel(0) = 1
      istartlev(0) = 1
      iboxlev(1) = 1

      call chxcin(ndeg,wsave2)
c
      do i = 0, maxlevel - 1
      iflag = 0
      istart = istartlev(i)
      iend = istart + nblevel(i) - 1
      do j = istart, iend
       ibox = iboxlev(j)

       call mkgrid(xf,yf,icolbox(ibox),irowbox(ibox),
     1      levelbox(ibox),ndeg,cent0,xsize0)

       do jj = 1, ndeg
         do ii = 1, ndeg
           ftemp(jj,ii) = h(xf(jj),yf(ii))
         end do
       end do


c      compute chebyshev transforms
       call getcoeff(ndeg,ndeg,ftemp,coefftemp,wsave2)

       call geterror(ndeg,coefftemp,error)

       hh = dble(4**levelbox(ibox))
       epscaled = eps * hh

c       write(*,*) error, epscaled

       if(error .ge. epscaled)then
c        call subdivide
         call subdivide1(ibox,iparentbox,ichildbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev,itemparray)
         iflag = 1
       endif
      end do
      if(iflag .eq. 0)then
c      nothing was divided at the
c      last level, so exit the loop.
       return
      endif
      end do
      return
      end subroutine




c******************************************************
c     subroutine mktreep
c******************************************************
c
c     This subroutine is used to generate a tree given
c     a discrete set of points. The algorithm works by
c     refining until there's no more than a given
c     number of points in each leaf box. Points are
c     reordered on output. The original id is saved.
c
c     INPUT:
c     maxboxes: the maximum number of boxes allowed
c     maxlevel: deepest level allowed
c     maxppl: the maximum number of boundary points
c             allowed in leaf box
c     npts: number of points
c     xpts: coordinates of points
c     itemparray: a dummy array to be used in 
c             subroutine subdivide1, which is called within
c     cent0: center of the root box
c     xsize0: side length of the root box 
c     
c     OUTPUT: 
c     (tree structure:)
c     istartlev is the pointer to where each level
c               begins in the iboxlev array
c     levelbox is an array determining the level of each box
c     nboxes is the total number of boxes
c     nlev is the finest level
c     icolbox denotes the column of each box
c     irowbox denotes the row of each box
c     iparentbox denotes the parent of each box
c     ichildbox denotes the four children of each box
c     nblevel is the total number of boxes per level
c     iboxlev is the array in which the boxes are arranged
c
c     (tree-point relation:)
c     npbox: number of sources in each box 
c                     (including non-leaf boxes)
c     ipold: ids of points on input
c            (before the overwriting)
c     istartbox: starting point in the array xpts
c                for each box
c     ibofp: indices of leaf boxes that the point
c            belongs to (after overwriting)
c
c     (xpts: overwritten, reordered so that points in
c            a box are stored in adjacent places)
c
c
c*****************************************************
c
      subroutine mktreep(levelbox,icolbox,irowbox,nboxes,
     1           nlev,iparentbox,ichildbox,nblevel,iboxlev,
     2           istartlev,npbox,ipold,istartbox,
     3           ibofp,maxboxes,itemparray,maxlevel,
     4           maxppl,npts,xpts,cent0,xsize0)
      implicit real*8 (a-h,o-z)
c-----global variables
      integer nboxes, nlev, maxboxes, maxlevel
      integer maxppl, npts
      integer levelbox(1)
      integer icolbox(1), irowbox(1)
      integer iparentbox(1), ichildbox(4,1)
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer itemparray(1)
      integer npbox(1), ipold(1), istartbox(1), ibofp(1)
      real*8 xpts(2,npts)
      real*8 cent0(2), xsize0
c-----local variables
      integer nsinchild(4)
      integer, allocatable:: isinchild(:,:)
      integer, allocatable:: itmp(:)
      real*8, allocatable:: xtmp(:,:)
c
c-----------------------------------------------------
c
      allocate(isinchild(npts,4))
      allocate(itmp(npts))
      allocate(xtmp(2,npts))
c
      do i = 0, maxlevel
        nblevel(i) = 0
        istartlev(i) = 0
      end do
      do i = 1, maxboxes
        iboxlev(i) = 0
      end do

c     first set the big parent box to the 
c     appropriate settings:
c     (initially, there is just one box and
c     it is the big parent box at level 0)
      ibox = 1
      levelbox(ibox) = 0
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = -1
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1

      nboxes = 1
      nlev = 0

      nblevel(0) = 1
      istartlev(0) = 1
      iboxlev(1) = 1
c
c-----------------------
c
      do i=1,nboxes
        npbox(i)=0
        istartbox(i)=-1
      enddo 
      
      do j=1,npts
        ipold(j)=j
      enddo

      npbox(1)=npts
      istartbox(1)=1
c
c         the root box contains all the points
c
c------------------------------------
c
      ix=1
      iy=1
      xlength=xsize0
      x0=cent0(1)-xsize0/2.0d0
      y0=cent0(2)-xsize0/2.0d0
c
      nside=1
c
      do i=0,maxlevel-1
        iflag=0
c           iflag: if the last level has been divided at all
        xlength=xlength/2.0d0
        nside=nside*2
c
        istart=istartlev(i)
        iend=istart+nblevel(i)-1
        do j=istart, iend
          ibox=iboxlev(j)
          if(npbox(ibox) .gt. maxppl) then
c            write(*,*) 'box:', ibox
c            write(*,*) 'number of sources in box:',nsinbox(ibox)
c--------------
c           call subdivide1 to subdivide the current box
            call subdivide1(ibox,iparentbox,ichildbox,
     1           nboxes,irowbox,icolbox,levelbox,nlev,
     2           istartlev, nblevel, iboxlev, itemparray)
c
c           sort points in ibox into its newborn babies
c           attention: this time physically rewrite
c           the array xpts so that points in the same
c           box are stored in adjacent places
c
            iss=istartbox(ibox)
            ies=istartbox(ibox)+npbox(ibox)-1
c
            do kk=1,4
              nsinchild(kk)=0
c             initialize the number of points in each kid
c             to be zero
              ic=ichildbox(kk,ibox)
              npbox(ic)=0
            enddo
c
            do jj=iss,ies
c             if we don't reorder along the way
c             this part is ugly
              xx=xpts(1,jj)              
              yy=xpts(2,jj)              
c
              xtmp(1,jj)=xx
              xtmp(2,jj)=yy
              itmp(jj)=ipold(jj)
c                 iss-th to ies-th entries of xpts
c                 and ipold copied to xtmp and itmp at
c                 corresponding places
c
ccccccccc----------------------------------
c              ix=ceiling((xx+0.5d0)/xlength)
c              iy=ceiling((yy+0.5d0)/xlength)
              ix=ceiling((xx-x0)/xlength)
              iy=ceiling((yy-y0)/xlength)
ccccccccc----------------------------------
c              nstart=nside/2+1
c              ix=nstart+ceiling(xx/xlength)
c              iy=nstart+ceiling(yy/xlength)
ccccccccc----------------------------------

              if(ix.le.0) then
                ix=1
              endif

              if(ix.gt.nside) then
                ix=nside
              endif

              if(iy.le.0) then
                iy=1
              endif

              if(iy.gt.nside) then
                iy=nside
              endif
c               x and y index of the point 
c               on the children's level
c
              do kk=1,4
                ic=ichildbox(kk,ibox)
                if((icolbox(ic).eq.ix).and.(irowbox(ic).eq.iy)) then
                  nsinchild(kk)=nsinchild(kk)+1
                  npbox(ic)=npbox(ic)+1
                  isinchild(nsinchild(kk),kk)=jj
                endif
              enddo
            enddo
cccccc------
c
c         all points in box ibox sorted into isinchild,
c         indicating which child the points belong to. 
c         now overwrite the iss-th to ies-th spots in
c         xpts and ipold together
c
cccccc------
            jpt=iss
            do jc=1,4
              ic=ichildbox(jc,ibox)
              istartbox(ic)=jpt
              do js=1,nsinchild(jc)
c                issorted(jpt)=isinchild(js,jc)
                idnow=isinchild(js,jc)
                xpts(1,jpt)=xtmp(1,idnow)  
                xpts(2,jpt)=xtmp(2,idnow)
                ipold(jpt)=itmp(idnow)
                jpt=jpt+1
              enddo
            enddo
c           done with box: ibox
c--------------
            iflag=1
          endif      
        enddo 

        if(iflag .eq. 0) then
c         nothing was divides at the last level,
c         exit the loop (on levels)
          goto 123
        endif
      enddo
c
c     now get information for ibofp
123   do l=0,nlev
        istart=istartlev(l)
        iend=istartlev(l)+nblevel(l)-1
        do ii=istart,iend
          ibox=iboxlev(ii)
          if((ichildbox(1,ibox).lt.0).and.(npbox(ibox).gt.0)) then
            iss=istartbox(ibox)
            ies=iss+npbox(ibox)-1
            do kk=iss,ies
              ibofp(kk)=ibox
            enddo
          endif
        enddo
      enddo
c
      deallocate(isinchild)
      deallocate(itmp)
      deallocate(xtmp)

      return 
      end subroutine




c**********************************************************
c     subroutine mkcolls
c**********************************************************
c     the following subroutine is used to generate the colleagues
c     for all of the boxes in the tree structure.  if a colleague
c     doesn't exist it is set to -1.  each box has nine colleagues
c     and they are ordered as follows:
c
c                        7     8     9
c               
c                        4     5     6
c
c                        1     2     3
c
c     you are your own colleague number 5.
c     the algorithm used here is recursive and takes advantage of
c     the fact that your colleagues can only be the children of 
c     your parents colleagues.  there is no need to scan all of the
c     boxes.  iperiod denotes whether or not we are in a periodic
c     or free space case.  the basic algorithm is the same, but in the
c     periodic case we have to look for boxes that are 'outside' of
c     the standard size box.
c
c
c     input:
c
c     levelbox is an array determining the level of each box
c
c     icolbox denotes the column of each box
c
c     irowbox denotes the row of each box
c
c     nboxes is the total number of boxes
c
c     nlev is the finest level
c
c     iparentbox denotes the parent of each box
c
c     ichildbox denotes the four children of each box
c
c     nblevel is the total number of boxes per level
c
c     iboxlev is the array in which the boxes are arranged
c
c     istartlev is the pointer to where each level begins in the
c               iboxlev array
c
c     iperiod denotes what kind of colleagues are to be generated
c             iperiod = 0 : free space
c             iperiod = 1 or 2 : periodic
c             iperiod = 3 : periodic up/down and free space left/right
c
c     output:
c
c     icolleagbox denotes the colleagues of a given box
c
c**********************************************************
c
      subroutine mkcolls(icolbox,
     1      irowbox, icolleagbox, nboxes, nlev,
     2      iparentbox, ichildbox, nblevel,
     3      iboxlev, istartlev, iperiod)
      implicit none
c-----global variables
      integer *4 icolleagbox(9,1)
      integer *4 icolbox(1), irowbox(1)
      integer *4 nboxes, nlev, iparentbox(1)
      integer *4 ichildbox(4,1)
      integer *4 nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer *4 iperiod
c-----local variables
      integer *4 colleague, partemp
      integer *4 jcntr, ibox, itest
      integer *4 icntr, ilev, j, l, nside
      integer *4 irowtemp, icoltemp
      integer *4 irowtest, icoltest

c     initialize colleague number 5 to
c     yourself and all other colleagues to
c     -1.  -1 is the flag for the case when
c     the colleagues don't exist.  it can 
c     be overwritten below. 
      do ibox = 1, nboxes
       icolleagbox(5,ibox) = ibox
       do j = 1, 4
         icolleagbox(j,ibox) = -1
       end do
       do j = 6, 9
         icolleagbox(j,ibox) = -1
       end do
      end do


c     scan through all of the levels except the coarsest level.
c     the one box at this level cannot have any colleagues.
c     do the uniform case first:
      if(iperiod .eq. 0)then
      do ilev = 1, nlev
c      scan through all of the boxes on each level.  for each test
c      box, scan the parent's colleagues and test to see if 
c      their children are in contact with the box being tested.
c      each colleague is placed in the correct order, so there is
c      no need to 'shuffle' them later on.
       do l = istartlev(ilev), istartlev(ilev) + nblevel(ilev) - 1
         ibox    = iboxlev(l)
         partemp = iparentbox(ibox)

c        irowtemp and icoltemp denote the row and column of
c        the test box.
         irowtemp = irowbox(ibox)
         icoltemp = icolbox(ibox)
      
         do 100 jcntr = 1, 9
c          colleague denotes the colleague of the parent box.
           colleague = icolleagbox(jcntr,partemp)
c          if the colleague doesn't exist
c          or is childless, skip it:
           if (colleague .lt. 0)goto 100
           if (ichildbox(1,colleague) .lt. 0)goto 100
           do icntr = 1, 4
             j = ichildbox(icntr,colleague)
c            irowtest and icoltest denote the row and column of
c            the box being compared to the test box.
             irowtest = irowbox(j)
             icoltest = icolbox(j)

             if(irowtemp .eq. irowtest+1)then
               if(icoltemp .eq. icoltest+1)then
                 icolleagbox(1,ibox) = j
               elseif(icoltemp .eq. icoltest)then
                 icolleagbox(2,ibox) = j
               elseif(icoltemp .eq. icoltest-1)then
                 icolleagbox(3,ibox) = j
               endif
             elseif(irowtemp .eq. irowtest)then
               if(icoltemp .eq. icoltest+1)then
                 icolleagbox(4,ibox) = j
               elseif(icoltemp .eq. icoltest-1)then
                 icolleagbox(6,ibox) = j
               endif
             elseif(irowtemp .eq. irowtest-1)then
               if(icoltemp .eq. icoltest+1)then
                 icolleagbox(7,ibox) = j
               elseif(icoltemp .eq. icoltest)then
                 icolleagbox(8,ibox) = j
               elseif(icoltemp .eq. icoltest-1)then
                 icolleagbox(9,ibox) = j
               endif
             endif
          end do 
100      continue
       end do
      end do 

c     now compute the colleagues in
c     the periodic case, if necessary:
      elseif(iperiod .eq. 1 .or. iperiod .eq. 2)then
c     initialize the first box (level 0) so
c     that it has its own colleagues.
c     this is necessary, because at deeper 
c     levels, the algorithm works by scanning 
c     the parent boxes colleagues.
      ibox = iboxlev(istartlev(0))
      icolleagbox(1,ibox) = ibox
      icolleagbox(2,ibox) = ibox
      icolleagbox(3,ibox) = ibox
      icolleagbox(4,ibox) = ibox
      icolleagbox(5,ibox) = ibox
      icolleagbox(6,ibox) = ibox
      icolleagbox(7,ibox) = ibox
      icolleagbox(8,ibox) = ibox
      icolleagbox(9,ibox) = ibox

      do ilev = 1, nlev
      nside = 2**ilev
c      scan through all of the boxes on each level.  for each test
c      box, scan the parent's colleagues and test to see if
c      their children are in contact with the box being tested.
c      each colleague is placed in the correct order, so there is
c      no need to 'shuffle' them later on.
       do l = istartlev(ilev), istartlev(ilev) + nblevel(ilev) - 1
        ibox = iboxlev(l)

c       irowtemp and icoltemp denote the
c       row and column of the test box.
        irowtemp = irowbox(ibox)
        icoltemp = icolbox(ibox)

c       irowtest and icoltest denote the row and column of
c       the box being compared to the test box.

        do 300 jcntr = 1, 9
c         first determine the column and row numbers
c         of all of the potential colleagues:
          if(jcntr .eq. 5)goto 300
          if(jcntr .eq. 1)then
            icoltest = icoltemp - 1
            irowtest = irowtemp - 1
          elseif(jcntr .eq. 2)then
            icoltest = icoltemp
            irowtest = irowtemp - 1
          elseif(jcntr .eq. 3)then
            icoltest = icoltemp + 1
            irowtest = irowtemp - 1
          elseif(jcntr .eq. 4)then
            icoltest = icoltemp - 1
            irowtest = irowtemp
          elseif(jcntr .eq. 6)then
            icoltest = icoltemp + 1
            irowtest = irowtemp
          elseif(jcntr .eq. 7)then
            icoltest = icoltemp - 1
            irowtest = irowtemp + 1
          elseif(jcntr .eq. 8)then
            icoltest = icoltemp
            irowtest = irowtemp + 1
          elseif(jcntr .eq. 9)then
            icoltest = icoltemp + 1
            irowtest = irowtemp + 1
          endif

c         now test to see if the test parameters 
c         lie in the domain:
c         (if they are outside of the domain, just
c         add or subtract the appropriate number so
c         that the boxes 'wrap around.')
          if(icoltest .lt. 1)then
            icoltest = icoltest + nside
          elseif(icoltest .gt. nside)then
            icoltest = icoltest - nside
          endif
          if(irowtest .lt. 1)then
            irowtest = irowtest + nside
          elseif(irowtest .gt. nside)then
            irowtest = irowtest - nside
          endif


       do 200 j = 1, 9
        if(icolleagbox(j,iparentbox(ibox)) .lt. 0)goto 200
        if(ichildbox(1,icolleagbox(j,iparentbox(ibox))) .lt. 0)goto 200
          do icntr = 1, 4
            itest = ichildbox(icntr,icolleagbox(j,iparentbox(ibox)))
            if(irowbox(itest) .eq. irowtest
     1         .and. icolbox(itest) .eq. icoltest)then
               icolleagbox(jcntr,ibox) = itest
c              assign ishcoll(1,jcntr,ibox)=ishiftx
c              assign ishcoll(2,jcntr,ibox)=ishifty
            endif
          end do
200    continue
300    continue
      end do
      end do
      elseif(iperiod .eq. 3)then
c     initialize the first box (level 0) so
c     that it has its own colleagues.
c     this is necessary, because at deeper 
c     levels, the algorithm works by scanning 
c     the parent boxes colleagues.
      ibox = iboxlev(istartlev(0))
      icolleagbox(1,ibox) = ibox
      icolleagbox(2,ibox) = ibox
      icolleagbox(3,ibox) = ibox
      icolleagbox(4,ibox) = ibox
      icolleagbox(5,ibox) = ibox
      icolleagbox(6,ibox) = ibox
      icolleagbox(7,ibox) = ibox
      icolleagbox(8,ibox) = ibox
      icolleagbox(9,ibox) = ibox

      do ilev = 1, nlev
      nside = 2**ilev
c      scan through all of the boxes on each level.  for each test
c      box, scan the parent's colleagues and test to see if
c      their children are in contact with the box being tested.
c      each colleague is placed in the correct order, so there is
c      no need to 'shuffle' them later on.
       do l = istartlev(ilev), istartlev(ilev) + nblevel(ilev) - 1
        ibox = iboxlev(l)

c       irowtemp and icoltemp denote the
c       row and column of the test box.
        irowtemp = irowbox(ibox)
        icoltemp = icolbox(ibox)

c       irowtest and icoltest denote the row and column of
c       the box being compared to the test box.

        do 500 jcntr = 1, 9
c         first determine the column and row numbers
c         of all of the potential colleagues:
          if(jcntr .eq. 5)goto 500
          if(jcntr .eq. 1)then
            icoltest = icoltemp - 1
            irowtest = irowtemp - 1
          elseif(jcntr .eq. 2)then
            icoltest = icoltemp
            irowtest = irowtemp - 1
          elseif(jcntr .eq. 3)then
            icoltest = icoltemp + 1
            irowtest = irowtemp - 1
          elseif(jcntr .eq. 4)then
            icoltest = icoltemp - 1
            irowtest = irowtemp
          elseif(jcntr .eq. 6)then
            icoltest = icoltemp + 1
            irowtest = irowtemp
          elseif(jcntr .eq. 7)then
            icoltest = icoltemp - 1
            irowtest = irowtemp + 1
          elseif(jcntr .eq. 8)then
            icoltest = icoltemp
            irowtest = irowtemp + 1
          elseif(jcntr .eq. 9)then
            icoltest = icoltemp + 1
            irowtest = irowtemp + 1
          endif

c         now test to see if the test parameters 
c         lie in the domain:
c         (if they are outside of the domain, just
c         add or subtract the appropriate number so
c         that the boxes 'wrap around.')
          if(irowtest .lt. 1)then
            irowtest = irowtest + nside
          elseif(irowtest .gt. nside)then
            irowtest = irowtest - nside
          endif


       do 400 j = 1, 9
        if(icolleagbox(j,iparentbox(ibox)) .lt. 0)goto 400
        if(ichildbox(1,icolleagbox(j,iparentbox(ibox))) .lt. 0)goto 400
          do icntr = 1, 4
            itest = ichildbox(icntr,icolleagbox(j,iparentbox(ibox)))
            if(irowbox(itest) .eq. irowtest
     1         .and. icolbox(itest) .eq. icoltest)then
               icolleagbox(jcntr,ibox) = itest
            endif
          end do
400    continue
500    continue
      end do
      end do
      endif
      return
      end subroutine




c**********************************************************
c
c     this subroutine will determine whether or not a given tree satisfies 
c     the level restriction.  if it doesn't, call fixtree to fix the
c     tree.
c
c**********************************************************
c
      subroutine restriction(levelbox,iparentbox,ichildbox,icolbox, 
     1             irowbox,icolleagbox,nboxes,nlev,
     2             nblevel,iboxlev,istartlev,iperiod,ifixflag)
      implicit none
c-----global variables
      integer *4 levelbox(1), icolleagbox(9,1)
      integer *4 iparentbox(1), ichildbox(4,1)
      integer *4 icolbox(1), irowbox(1)
      integer *4 nboxes, nlev
      integer *4 nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer *4 iperiod, ifixflag
c-----local variables
      integer *4 ichild(4),icoll(9), ibox
      integer *4 i, ipar, itest, j, nb
      integer *4 itemp


c     let's sort all of the boxes by level.
c     this takes the place of the ladder structure
c     in the uniform case.  all boxes are sorted
c     into the array and placed in their proper places.
      call sortboxes(levelbox,nboxes,nlev,
     1           nblevel,iboxlev,istartlev)


c     first let's call a subroutine that will
c     generate all of the colleagues for each
c     box.  the colleagues are generated in the
c     correct order so there is no need to 'shuffle'
c     them later on.
      call mkcolls(icolbox,
     1     irowbox,icolleagbox,nboxes,nlev,
     2     iparentbox,ichildbox,nblevel,
     3     iboxlev, istartlev,iperiod)


      do i = nlev, 2, -1
        do j = istartlev(i), istartlev(i) + nblevel(i) - 1
          ibox = iboxlev(j)
          ipar  = iparentbox(ibox)
          itest = iparentbox(ipar)

            icoll(1) = icolleagbox(1,itest)
            icoll(2) = icolleagbox(2,itest)
            icoll(3) = icolleagbox(3,itest)
            icoll(4) = icolleagbox(4,itest)
            icoll(5) = icolleagbox(5,itest)
            icoll(6) = icolleagbox(6,itest)
            icoll(7) = icolleagbox(7,itest)
            icoll(8) = icolleagbox(8,itest)
            icoll(9) = icolleagbox(9,itest)

            ichild(1) = ichildbox(1,itest)
            ichild(2) = ichildbox(2,itest)
            ichild(3) = ichildbox(3,itest)
            ichild(4) = ichildbox(4,itest)


          do nb = 1, 9
            itemp = icoll(nb)
            if(ichildbox(1,itemp) .lt. 0)then
c             the neighboring box is not divided
c             we could have problems.
              if (nb .eq. 1)then
                if(ipar .eq. ichild(4))then
                    ifixflag = 1
                end if
              elseif (nb .eq. 2)then
                if(ipar .eq. ichild(3) .or. ipar .eq. ichild(4))then
                    ifixflag = 1
                end if
              elseif (nb .eq. 3)then
                if(ipar .eq. ichild(3))then
                    ifixflag = 1
                end if
              elseif (nb .eq. 4)then
                if(ipar .eq. ichild(4) .or. ipar .eq. ichild(1))then
                    ifixflag = 1
                end if
              elseif (nb .eq. 6)then
                if(ipar .eq. ichild(2) .or. ipar .eq. ichild(3))then
                    ifixflag = 1
                end if
              elseif (nb .eq. 7)then 
                if(ipar .eq. ichild(1))then
                    ifixflag = 1
                end if
              elseif (nb .eq. 8)then
                if(ipar .eq. ichild(1) .or. ipar .eq. ichild(2))then
                    ifixflag = 1
                end if
              elseif (nb .eq. 9)then
                if(ipar .eq. ichild(2))then
                    ifixflag = 1
                end if
              endif
            endif
          end do
        end do
      end do

      return
      end subroutine




c**********************************************************
c     subroutine fixtree
c**********************************************************
c
c     the following subroutine is designed to take a correctly defined
c     tree and alter it so that no two boxes that contact each other
c     are more than one level apart.  this is corrected only by adding
c     boxes.  the procedure involves flagging down bigger boxes and
c     dividing them and their children as necessary.
c     this routine also produces an array of colleagues for each box
c     that is numbered in the correct order.  all of the children are set
c     so that they satisfy our ordering convention.
c     the algorithm in the periodic case works the same way, it is just 
c     that upon subdivision the new colleagues must be put down to 
c     account for the periodicity.
c
c     input:
c
c     levelbox is an array determining the level of each box
c
c     iparentbox denotes the parent of each box
c
c     ichildbox denotes the four children of each box
c
c     icolbox denotes the column of each box
c
c     irowbox denotes the row of each box
c
c     icolleagbox denotes the colleagues of a given box
c
c     nboxes is the total number of boxes
c
c     nlev is the finest level
c
c     nblevel is the total number of boxes per level
c
c     iboxlev is the array in which the boxes are arranged
c
c     istartlev is the pointer to where each level begins in the
c               iboxlev array
c
c     iperiod denotes what kind of colleagues are to be generated
c             iperiod = 0 : free space
c             iperiod = 1 or 2 : periodic
c             iperiod = 3 : periodic up/down and free space left/right
c
c     iflag is just a dummy array
c
c     maxboxes is the maximum number of boxes we have storage for
c
c     output:
c
c     icolbox, irowbox, icolleagbox, nboxes, and all other
c     arrays describing the tree may be change on output
c
c**********************************************************
c
      subroutine fixtree(levelbox,iparentbox,ichildbox,icolbox, 
     1           irowbox,icolleagbox,nboxes,nlev,
     2           nblevel,iboxlev,istartlev,iperiod,
     3           iflag, maxboxes,itemparray)
      implicit none
c-----global variables
      integer *4 levelbox(1), icolleagbox(9,1)
      integer *4 maxboxes
      integer *4 iparentbox(1), ichildbox(4,1)
      integer *4 icolbox(1), irowbox(1)
      integer *4 nboxes, nlev
      integer *4 nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer *4 iperiod
      integer *4 itemparray(1)
c-----local variables
      integer *4 iflag(maxboxes)
      integer *4 ichild(4),icoll(9), ibox
      integer *4 i, ipar, itest, j, nb
      integer *4 itemp, ntemp, jcntr, icntr
      integer *4 start, istop


c     let's sort all of the boxes by level.
c     this takes the place of the ladder structure
c     in the uniform case.  all boxes are sorted
c     into the array and placed in their proper places.
      call sortboxes(levelbox,nboxes,nlev,
     1           nblevel,iboxlev,istartlev)


c     first let's call a subroutine that will
c     generate all of the colleagues for each
c     box.  the colleagues are generated in the
c     correct order so there is no need to 'shuffle'
c     them later on.
      call mkcolls(icolbox,
     1       irowbox,icolleagbox,nboxes,nlev,
     2       iparentbox,ichildbox,nblevel,
     3       iboxlev, istartlev,iperiod)


c     let's initialize all of the flags to zero.
      do i = 1, nboxes
        iflag(i) = 0
      end do

c     find all of the boxes that need to be
c     flagged.  a flagged box will be denoted by 
c     setting iflag(box) = 1.
c     this refers to any box that is directly touching 
c     a box that is more than one level smaller than
c     it.  it is found by performing an upward pass
c     and looking a box's parents parents and seeing
c     if they are childless and contact the given box.
c     note that we only need to get up to level two, as
c     we will not find a violation at a coarser level
c     than that.
      do i = nlev, 2, -1
        do j = istartlev(i), istartlev(i) + nblevel(i) - 1
          ibox = iboxlev(j)
          ipar  = iparentbox(ibox)
          itest = iparentbox(ipar)

            icoll(1) = icolleagbox(1,itest)
            icoll(2) = icolleagbox(2,itest)
            icoll(3) = icolleagbox(3,itest)
            icoll(4) = icolleagbox(4,itest)
            icoll(5) = icolleagbox(5,itest)
            icoll(6) = icolleagbox(6,itest)
            icoll(7) = icolleagbox(7,itest)
            icoll(8) = icolleagbox(8,itest)
            icoll(9) = icolleagbox(9,itest)

            ichild(1) = ichildbox(1,itest)
            ichild(2) = ichildbox(2,itest)
            ichild(3) = ichildbox(3,itest)
            ichild(4) = ichildbox(4,itest)


          do nb = 1, 9
            itemp = icoll(nb)
            if(ichildbox(1,itemp) .lt. 0)then
c             the neighboring box is not divided
c             we could have problems.
              if (nb .eq. 1)then
                if(ipar .eq. ichild(4))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 2)then
                if(ipar .eq. ichild(3) .or. ipar .eq. ichild(4))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 3)then
                if(ipar .eq. ichild(3))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 4)then
                if(ipar .eq. ichild(4) .or. ipar .eq. ichild(1))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 6)then
                if(ipar .eq. ichild(2) .or. ipar .eq. ichild(3))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 7)then 
                if(ipar .eq. ichild(1))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 8)then
                if(ipar .eq. ichild(1) .or. ipar .eq. ichild(2))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 9)then
                if(ipar .eq. ichild(2))then
                    iflag(itemp) = 1
                end if
              endif
            endif
          end do
        end do
      end do


c     find all of the boxes that need to be
c     given a flag+.  a flag+ box will be denoted by 
c     setting iflag(box) = 2.
c     this refers to any box that is not already flagged
c     and is bigger than and is contacting a flagged box
c     or another box that has already been given a flag+.
c     it is found by performing an upward pass
c     and looking at a flagged box's parents colleagues
c     and a flag+ box's parents colleagues and seeing if
c     they are childless and present the case where a 
c     bigger box is contacting a flagged or a flag+ box.
      do i = nlev, 2, -1
        do j = istartlev(i), istartlev(i) + nblevel(i) - 1
         ibox = iboxlev(j)
          if(iflag(ibox) .eq. 1 .or. iflag(ibox) .eq. 2)then

          ipar  = iparentbox(ibox)
 
            icoll(1) = icolleagbox(1,ipar)
            icoll(2) = icolleagbox(2,ipar)
            icoll(3) = icolleagbox(3,ipar)
            icoll(4) = icolleagbox(4,ipar)
            icoll(5) = icolleagbox(5,ipar)
            icoll(6) = icolleagbox(6,ipar)
            icoll(7) = icolleagbox(7,ipar)
            icoll(8) = icolleagbox(8,ipar)
            icoll(9) = icolleagbox(9,ipar)

            ichild(1) = ichildbox(1,ipar)
            ichild(2) = ichildbox(2,ipar)
            ichild(3) = ichildbox(3,ipar)
            ichild(4) = ichildbox(4,ipar)
          

          do nb = 1, 9
            itemp = icoll(nb)
c           let's check using the same criteria as above, but noting that
c           a flag will take precedence over a flag+.
            if(ichildbox(1,itemp) .lt. 0 
     1          .and. iflag(itemp) .ne. 1)then
c             the neighboring box is not divided
c             we could have problems.
              if (nb .eq. 1)then
                if(ibox .eq. ichild(4))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 2)then
                if(ibox .eq. ichild(3) .or. ibox .eq. ichild(4))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 3)then
                if(ibox .eq. ichild(3))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 4)then
                if(ibox .eq. ichild(4) .or. ibox .eq. ichild(1))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 6)then
                if(ibox .eq. ichild(2) .or. ibox .eq. ichild(3))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 7)then
                if(ibox .eq. ichild(1))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 8)then
                if(ibox .eq. ichild(1) .or. ibox .eq. ichild(2))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 9)then
                if(ibox .eq. ichild(2))then
                    iflag(itemp) = 2
                end if
              endif
            endif
           end do
          endif
        end do
      end do


c     now let's divide the boxes that need to be immediately
c     divided up.  all of the flagged and flag+ boxes need to
c     be divided one time.  the distinction lies in the fact
c     that the children of a flag+ box will never need to be
c     divided but the children of a flagged box may need to 
c     be divided further.
c     below, all flagged and flag+ boxes are divided once.  the
c     children of a flag+ box are left unflagged while those of
c     the flagged boxes are given a flag++ (denoted by setting
c     iflag(box) = 3) which will be needed in the downward pass.     
      ntemp = nboxes
      do i = 1, ntemp
c      divide flagged boxes:
       if (iflag(i) .eq. 1)then

         if(ichildbox(1,i) .lt. 0)then
         call subdivide(i,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
         endif


c        give flag++ to children of flagged boxes.
         itemp = ichildbox(1,i)
         iflag(itemp) = 3

         itemp = ichildbox(2,i)
         iflag(itemp) = 3

         itemp = ichildbox(3,i)
         iflag(itemp) = 3

         itemp = ichildbox(4,i)
         iflag(itemp) = 3

c      divide flag+ boxes.
       elseif (iflag(i) .eq. 2)then
  
         if(ichildbox(1,i) .lt. 0)then
         call subdivide(i,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
         endif

       endif
      end do 


c     now we need to do a downward pass.
c     we will concern ourselves only with the children of
c     flagged boxes and their children.  at each level,
c     for each flag++ box, test colleagues children and see
c     if they have children that are contacting you.  if so,
c     divide and flag++ all children that are created.     

      do i = 0, nlev
      ntemp = nboxes
      start = istartlev(i)
      istop  = istartlev(i) + nblevel(i) - 1
      do 500 j = start, istop
       ibox = iboxlev(j)
c      only be concerned with boxes on this level and
c      boxes that are given a flag++:
       if(iflag(ibox) .ne. 3)goto 500

         icoll(1) = icolleagbox(1,ibox)
         icoll(2) = icolleagbox(2,ibox)
         icoll(3) = icolleagbox(3,ibox)
         icoll(4) = icolleagbox(4,ibox)
         icoll(5) = icolleagbox(5,ibox)
         icoll(6) = icolleagbox(6,ibox)
         icoll(7) = icolleagbox(7,ibox)
         icoll(8) = icolleagbox(8,ibox)
         icoll(9) = icolleagbox(9,ibox)


c       scan colleagues.
        do 400 jcntr = 1, 9
        if(icoll(jcntr) .lt. 0)goto 400
        if(ichildbox(1,icoll(jcntr)) .lt. 0)goto 400

         ichild(1) = ichildbox(1,icoll(jcntr))
         ichild(2) = ichildbox(2,icoll(jcntr))
         ichild(3) = ichildbox(3,icoll(jcntr))
         ichild(4) = ichildbox(4,icoll(jcntr))


c          scan colleague's children.
           do 300 icntr = 1, 4
           if (ichildbox(1,ichild(icntr)) .lt. 0)goto 300 

           if(jcntr .eq. 1 .and. icntr .eq. 2)then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
         endif



c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 2 .and. 
     1        (icntr .eq. 1 .or. icntr .eq. 2))then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
         endif



c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 3 .and. icntr .eq. 1)then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
         endif



c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 4 .and. 
     1        (icntr .eq. 2 .or. icntr .eq. 3))then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
         endif


c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 6 .and. 
     1        (icntr .eq. 1 .or. icntr .eq. 4))then
             
             
c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
         endif


c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 7 .and. icntr .eq. 3)then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
         endif


c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 8 .and. 
     1        (icntr .eq. 3 .or. icntr .eq. 4))then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
         endif


c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 9 .and. icntr .eq. 4)then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
         endif


c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3
           endif
300     continue   
400     continue
500     continue
      end do
      return
      end subroutine




c**********************************************************
c
c     This subroutine gives the position of a box in the tree
c     by returning the coordinates of the lower left vertex, and
c     side lengths in x and y directions
c     (assume that the tree is built on [-0.5,0.5]^2)
c     
c     INPUT:
c     icol: number of column (within the current level)
c     irow: number of row (within the current level)
c     level: the current level
c
c     OUTPUT:
c     (xll,yll): coordinates of the lower left vertex
c     dx: side length in x direction
c     dy: side length in y direction
c           
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




c**********************************************************
c
c     This subroutine generates the ndeg x ndeg tensor 
c     product Chebyshev grid on a given box
c
c     input:
c
c     icol denotes the column of the given box
c
c     irow denotes the row of the given box
c
c     level denotes the level of the given box
c
c     ndeg denotes the spatial order
c
c     cent0: center of the root box
c
c     xsize0: side length of the root box
c
c     output:
c
c     xf denotes the x values 
c     yf denotes the y values
c
c**********************************************************
c
      subroutine mkgrid(xf,yf,icol,irow,level,ndeg,
     1           cent0,xsize0) 
      implicit none
c-----global variables
      integer *4  level
      integer *4  icol, irow, ndeg
      real *8  xf(ndeg), yf(ndeg)
      real*8 cent0(2), xsize0
c-----local variables
      integer *4  i, j
      real *8  xstart
      real *8  xshift, yshift
      real *8  temp1, pi16, pi2n
      real *8  xx(ndeg), xscale(ndeg)

c       pi16 = dacos(-1.0d0) / 16.0d0
c       xx(1) = dcos(15.0d0*pi16) / 2.0d0
c       xx(2) = dcos(13.0d0*pi16) / 2.0d0
c       xx(3) = dcos(11.0d0*pi16) / 2.0d0
c       xx(4) = dcos( 9.0d0*pi16) / 2.0d0
c       xx(5) = dcos( 7.0d0*pi16) / 2.0d0
c       xx(6) = dcos( 5.0d0*pi16) / 2.0d0
c       xx(7) = dcos( 3.0d0*pi16) / 2.0d0
c       xx(8) = dcos( 1.0d0*pi16) / 2.0d0
       pi2n =dacos(-1.0d0) /dble(2.0d0*ndeg)
       do i=1, ndeg
         xx(i) =dcos((2.0d0*(ndeg-i)+1)*pi2n)/2.0d0
       enddo

       temp1 = xsize0 / dble(2**level)

c       xscale(1) = xx(1) * temp1 
c       xscale(2) = xx(2) * temp1 
c       xscale(3) = xx(3) * temp1
c       xscale(4) = xx(4) * temp1
c       xscale(5) = xx(5) * temp1
c       xscale(6) = xx(6) * temp1
c       xscale(7) = xx(7) * temp1 
c       xscale(8) = xx(8) * temp1
        do i=1, ndeg
          xscale(i) = xx(i) * temp1
        enddo


       xshift = cent0(1)-xsize0/2.0d0 + (icol-0.5d0) * temp1
       yshift = cent0(2)-xsize0/2.0d0 + (irow-0.5d0) * temp1

       do j = 1, ndeg
        xf(j) = xscale(j) + xshift
        yf(j) = xscale(j) + yshift
       end do

      return
      end subroutine





c**********************************************************
c
      subroutine getcoeff(n,m,fdat,coeff,wsave2)
      implicit real *8 (a-h,o-z)
      real *8 fdat(n,m),coeff(n,m)
      real *8 wsave2(1000)
      real *8 work(1000)
      real *8 f(1000),texp(1000)
c
c     transform rows
c
      do j = 1,m
         do i = 1,n
            f(i) = fdat(i,j)
         enddo
         call chexfc(f,n,texp,wsave2,work)
         do i = 1,n
            coeff(i,j) = texp(i)
         enddo
      enddo
c
c     transform columns
c
      do i = 1,n
         do j = 1,m
            f(j) = coeff(i,j)
         enddo
         call chexfc(f,n,texp,wsave2,work)
         do j = 1,m
            coeff(i,j) = texp(j)
         enddo
      enddo
      return
      end subroutine





c**********************************************************
c
c     the following subroutine calculates the l2 error between two vectors
c     that is needed in order to generate the tree structure.  this is used 
c     to determine how accurate the approximating polynomial is.
c
c     input:
c
c     xf denotes the x values of 36 cell centered points in the box
c
c     yf denotes the y values of 36 cell centered points in the box
c
c     xf2 denotes the x values of 144 cell centered points in the box
c
c     yf2 denotes the y values of 144 cell centered points in the box
c
c     level denotes the level of the given box
c
c     a is the matrix that maps 36 function values
c       onto the 21 basis function coefficients
c
c     output:
c
c     error is the l2 error between the exact solution evaluated
c     at the 144 points and the approximation given by the polynomial
c     (coefficients found from the 36 points)
c
c**********************************************************
c
      subroutine geterror(ndeg,coefftemp,error)
      implicit none
      integer ndeg, i, j
      real *8  error, coefftemp(ndeg,ndeg)
c-----local variables

c      just sum over the antidiagonal: (all (ndeg-1)th degree coefficients)
        error = 0.0d0
c        error = error + dabs(coefftemp(8,1))
c        error = error + dabs(coefftemp(7,2))
c        error = error + dabs(coefftemp(6,3))
c        error = error + dabs(coefftemp(5,4))
c        error = error + dabs(coefftemp(4,5))
c        error = error + dabs(coefftemp(3,6))
c        error = error + dabs(coefftemp(2,7))
c        error = error + dabs(coefftemp(1,8))
        do i=1, ndeg
          error=error+dabs(coefftemp(ndeg-i+1,i))
        enddo

c      with a weighted version of the off diagonal terms
c        error = error + dabs(coefftemp(7,1))/4.0d0
c        error = error + dabs(coefftemp(6,2))/4.0d0
c        error = error + dabs(coefftemp(5,3))/4.0d0
c        error = error + dabs(coefftemp(4,4))/4.0d0
c        error = error + dabs(coefftemp(3,5))/4.0d0
c        error = error + dabs(coefftemp(2,6))/4.0d0
c        error = error + dabs(coefftemp(1,7))/4.0d0
        do i=1, ndeg-1
         error=error+dabs(coefftemp(ndeg-i,i))/4.0d0
        enddo
c
c       modification: change to the relative error
c       by dividing abs(coefftemp(1,1)+coefftemp(1,2)
c       +coefftemp(2,1))
c         write(*,*) dabs(coefftemp(1,1)+coefftemp(1,2)+coefftemp(2,1))
        error=error/dabs(coefftemp(1,1))
c        error=error/dabs(coefftemp(1,1)+coefftemp(1,2)
c     1                                 +coefftemp(2,1))

      return
      end subroutine




c**********************************************************
c
c     the following subroutine sets up a structure that is analogous
c     to the 'ladder' structure in the nonadaptive case.
c     it is just a way of organizing the boxes by level in one long
c     array and denoting where in the array the levels change.
c
c     input:
c
c     levelbox is an array determining the level of each box
c
c     nboxes is the total number of boxes
c
c     nlev is the finest level
c
c     ichildbox denotes the four children of each box
c
c     output:
c
c     nblevel is the total number of boxes per level
c
c     iboxlev is the array in which the boxes are arranged
c
c     istartlev is the pointer to where each level begins in the
c               iboxlev array
c
c**********************************************************
c
      subroutine sortboxes(levelbox,nboxes,nlev,
     1           nblevel,iboxlev,istartlev)
      implicit none
c-----global variables
      integer *4 levelbox(1), nboxes, nlev
      integer *4 nblevel(0:1), iboxlev(1), istartlev(0:1)
c-----local variables
      integer *4 ncntr, i, j

        ncntr = 1
        do i = 0, nlev
          nblevel(i) = 0
          istartlev(i) = ncntr
          do j = 1, nboxes
            if(levelbox(j) .eq. i)then
             iboxlev(ncntr) = j
             ncntr = ncntr + 1
             nblevel(i) = nblevel(i) + 1
            endif
          end do
        end do
      return
      end subroutine




c**********************************************************
c     subroutine fixtreenf
c**********************************************************
c
c     the following subroutine is designed to take a correctly defined
c     tree and alter it so that no two boxes that contact each other
c     are more than one level apart.  this is corrected only by adding
c     boxes.  the procedure involves flagging down bigger boxes and
c     dividing them and their children as necessary.
c     this routine also produces an array of colleagues for each box
c     that is numbered in the correct order.  all of the children are set
c     so that they satisfy our ordering convention.
c     the algorithm in the periodic case works the same way, it is just 
c     that upon subdivision the new colleagues must be put down to 
c     account for the periodicity.
c
c     input:
c
c     levelbox is an array determining the level of each box
c
c     iparentbox denotes the parent of each box
c
c     ichildbox denotes the four children of each box
c
c     icolbox denotes the column of each box
c
c     irowbox denotes the row of each box
c
c     icolleagbox denotes the colleagues of a given box
c
c     nboxes is the total number of boxes
c
c     nlev is the finest level
c
c     nblevel is the total number of boxes per level
c
c     iboxlev is the array in which the boxes are arranged
c
c     istartlev is the pointer to where each level begins in the
c               iboxlev array
c
c     iperiod denotes what kind of colleagues are to be generated
c             iperiod = 0 : free space
c             iperiod = 1 or 2 : periodic
c             iperiod = 3 : periodic up/down and free space left/right
c
c     iflag is just a dummy array
c
c     maxboxes is the maximum number of boxes we have storage for
c
c     output:
c
c     icolbox, irowbox, icolleagbox, nboxes, and all other
c     arrays describing the tree may be change on output
c
c     fval may be overwritten too (interpolated to the new grid)
c
c**********************************************************
      subroutine fixtreenf(levelbox,iparentbox,ichildbox,
     1           icolbox,irowbox,icolleagbox,nboxes,nlev,
     2           nblevel,iboxlev,istartlev,iperiod,
     3           iflag, maxboxes,itemparray,ndeg,fval)
      implicit none
c-----global variables
      integer *4 levelbox(1), icolleagbox(9,1)
      integer *4 maxboxes, ndeg
      integer *4 iparentbox(1), ichildbox(4,1)
      integer *4 icolbox(1), irowbox(1)
      integer *4 nboxes, nlev
      integer *4 nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer *4 iperiod
      integer *4 iflag(maxboxes)
      integer *4 itemparray(1)
      real *8 fval(ndeg*ndeg, 1)
c-----local variables
      integer *4 ichild(4),icoll(9), ibox
      integer *4 i, ipar, itest, j, nb
      integer *4 itemp, ntemp, jcntr, icntr
      integer *4 start, istop
      integer *4 ik1, ik2, ik3, ik4


c     let's sort all of the boxes by level.
c     this takes the place of the ladder structure
c     in the uniform case.  all boxes are sorted
c     into the array and placed in their proper places.
      call sortboxes(levelbox,nboxes,nlev,
     1           nblevel,iboxlev,istartlev)


c     first let's call a subroutine that will
c     generate all of the colleagues for each
c     box.  the colleagues are generated in the
c     correct order so there is no need to 'shuffle'
c     them later on.
      call mkcolls(icolbox,
     1       irowbox,icolleagbox,nboxes,nlev,
     2       iparentbox,ichildbox,nblevel,
     3       iboxlev, istartlev,iperiod)


c     let's initialize all of the flags to zero.
      do i = 1, nboxes
        iflag(i) = 0
      end do

c     find all of the boxes that need to be
c     flagged.  a flagged box will be denoted by 
c     setting iflag(box) = 1.
c     this refers to any box that is directly touching 
c     a box that is more than one level smaller than
c     it.  it is found by performing an upward pass
c     and looking a box's parents parents and seeing
c     if they are childless and contact the given box.
c     note that we only need to get up to level two, as
c     we will not find a violation at a coarser level
c     than that.
      do i = nlev, 2, -1
        do j = istartlev(i), istartlev(i) + nblevel(i) - 1
          ibox = iboxlev(j)
          ipar  = iparentbox(ibox)
          itest = iparentbox(ipar)

            icoll(1) = icolleagbox(1,itest)
            icoll(2) = icolleagbox(2,itest)
            icoll(3) = icolleagbox(3,itest)
            icoll(4) = icolleagbox(4,itest)
            icoll(5) = icolleagbox(5,itest)
            icoll(6) = icolleagbox(6,itest)
            icoll(7) = icolleagbox(7,itest)
            icoll(8) = icolleagbox(8,itest)
            icoll(9) = icolleagbox(9,itest)

            ichild(1) = ichildbox(1,itest)
            ichild(2) = ichildbox(2,itest)
            ichild(3) = ichildbox(3,itest)
            ichild(4) = ichildbox(4,itest)


          do nb = 1, 9
            itemp = icoll(nb)
            if(ichildbox(1,itemp) .lt. 0)then
c             the neighboring box is not divided
c             we could have problems.
              if (nb .eq. 1)then
                if(ipar .eq. ichild(4))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 2)then
                if(ipar .eq. ichild(3) .or. ipar .eq. ichild(4))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 3)then
                if(ipar .eq. ichild(3))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 4)then
                if(ipar .eq. ichild(4) .or. ipar .eq. ichild(1))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 6)then
                if(ipar .eq. ichild(2) .or. ipar .eq. ichild(3))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 7)then 
                if(ipar .eq. ichild(1))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 8)then
                if(ipar .eq. ichild(1) .or. ipar .eq. ichild(2))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 9)then
                if(ipar .eq. ichild(2))then
                    iflag(itemp) = 1
                end if
              endif
            endif
          end do
        end do
      end do


c     find all of the boxes that need to be
c     given a flag+.  a flag+ box will be denoted by 
c     setting iflag(box) = 2.
c     this refers to any box that is not already flagged
c     and is bigger than and is contacting a flagged box
c     or another box that has already been given a flag+.
c     it is found by performing an upward pass
c     and looking at a flagged box's parents colleagues
c     and a flag+ box's parents colleagues and seeing if
c     they are childless and present the case where a 
c     bigger box is contacting a flagged or a flag+ box.
      do i = nlev, 2, -1
        do j = istartlev(i), istartlev(i) + nblevel(i) - 1
         ibox = iboxlev(j)
          if(iflag(ibox) .eq. 1 .or. iflag(ibox) .eq. 2)then

          ipar  = iparentbox(ibox)
 
            icoll(1) = icolleagbox(1,ipar)
            icoll(2) = icolleagbox(2,ipar)
            icoll(3) = icolleagbox(3,ipar)
            icoll(4) = icolleagbox(4,ipar)
            icoll(5) = icolleagbox(5,ipar)
            icoll(6) = icolleagbox(6,ipar)
            icoll(7) = icolleagbox(7,ipar)
            icoll(8) = icolleagbox(8,ipar)
            icoll(9) = icolleagbox(9,ipar)

            ichild(1) = ichildbox(1,ipar)
            ichild(2) = ichildbox(2,ipar)
            ichild(3) = ichildbox(3,ipar)
            ichild(4) = ichildbox(4,ipar)
          

          do nb = 1, 9
            itemp = icoll(nb)
c           let's check using the same criteria as above, but noting that
c           a flag will take precedence over a flag+.
            if(ichildbox(1,itemp) .lt. 0 
     1          .and. iflag(itemp) .ne. 1)then
c             the neighboring box is not divided
c             we could have problems.
              if (nb .eq. 1)then
                if(ibox .eq. ichild(4))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 2)then
                if(ibox .eq. ichild(3) .or. ibox .eq. ichild(4))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 3)then
                if(ibox .eq. ichild(3))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 4)then
                if(ibox .eq. ichild(4) .or. ibox .eq. ichild(1))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 6)then
                if(ibox .eq. ichild(2) .or. ibox .eq. ichild(3))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 7)then
                if(ibox .eq. ichild(1))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 8)then
                if(ibox .eq. ichild(1) .or. ibox .eq. ichild(2))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 9)then
                if(ibox .eq. ichild(2))then
                    iflag(itemp) = 2
                end if
              endif
            endif
           end do
          endif
        end do
      end do


c     now let's divide the boxes that need to be immediately
c     divided up.  all of the flagged and flag+ boxes need to
c     be divided one time.  the distinction lies in the fact
c     that the children of a flag+ box will never need to be
c     divided but the children of a flagged box may need to 
c     be divided further.
c     below, all flagged and flag+ boxes are divided once.  the
c     children of a flag+ box are left unflagged while those of
c     the flagged boxes are given a flag++ (denoted by setting
c     iflag(box) = 3) which will be needed in the downward pass.     
      ntemp = nboxes
      do i = 1, ntemp
c      divide flagged boxes:
       if (iflag(i) .eq. 1)then

         if(ichildbox(1,i) .lt. 0)then
         call subdivide(i,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c        after calling subdivide or subdivide1, assign fval
         ik1=ichildbox(1,i)
         ik2=ichildbox(2,i)
         ik3=ichildbox(3,i)
         ik4=ichildbox(4,i)
         call chebyp2k(ndeg,fval(1,i),fval(1,ik1),fval(1,ik2),
     1                 fval(1,ik3),fval(1,ik4))
c------------------------
         endif


c        give flag++ to children of flagged boxes.
         itemp = ichildbox(1,i)
         iflag(itemp) = 3

         itemp = ichildbox(2,i)
         iflag(itemp) = 3

         itemp = ichildbox(3,i)
         iflag(itemp) = 3

         itemp = ichildbox(4,i)
         iflag(itemp) = 3

c      divide flag+ boxes.
       elseif (iflag(i) .eq. 2)then
  
         if(ichildbox(1,i) .lt. 0)then
         call subdivide(i,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c        after calling subdivide or subdivide1, assign fval
         ik1=ichildbox(1,i)
         ik2=ichildbox(2,i)
         ik3=ichildbox(3,i)
         ik4=ichildbox(4,i)
         call chebyp2k(ndeg,fval(1,i),fval(1,ik1),fval(1,ik2),
     1                 fval(1,ik3),fval(1,ik4))
c------------------------
         endif

       endif
      end do 


c     now we need to do a downward pass.
c     we will concern ourselves only with the children of
c     flagged boxes and their children.  at each level,
c     for each flag++ box, test colleagues children and see
c     if they have children that are contacting you.  if so,
c     divide and flag++ all children that are created.     

      do i = 0, nlev
      ntemp = nboxes
      start = istartlev(i)
      istop  = istartlev(i) + nblevel(i) - 1
      do 500 j = start, istop
       ibox = iboxlev(j)
c      only be concerned with boxes on this level and
c      boxes that are given a flag++:
       if(iflag(ibox) .ne. 3)goto 500

         icoll(1) = icolleagbox(1,ibox)
         icoll(2) = icolleagbox(2,ibox)
         icoll(3) = icolleagbox(3,ibox)
         icoll(4) = icolleagbox(4,ibox)
         icoll(5) = icolleagbox(5,ibox)
         icoll(6) = icolleagbox(6,ibox)
         icoll(7) = icolleagbox(7,ibox)
         icoll(8) = icolleagbox(8,ibox)
         icoll(9) = icolleagbox(9,ibox)


c       scan colleagues.
        do 400 jcntr = 1, 9
        if(icoll(jcntr) .lt. 0)goto 400
        if(ichildbox(1,icoll(jcntr)) .lt. 0)goto 400

         ichild(1) = ichildbox(1,icoll(jcntr))
         ichild(2) = ichildbox(2,icoll(jcntr))
         ichild(3) = ichildbox(3,icoll(jcntr))
         ichild(4) = ichildbox(4,icoll(jcntr))


c          scan colleague's children.
           do 300 icntr = 1, 4
           if (ichildbox(1,ichild(icntr)) .lt. 0)goto 300 

           if(jcntr .eq. 1 .and. icntr .eq. 2)then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c          after calling subdivide or subdivide1, assign fval
           ik1=ichildbox(1,ibox)
           ik2=ichildbox(2,ibox)
           ik3=ichildbox(3,ibox)
           ik4=ichildbox(4,ibox)
           call chebyp2k(ndeg,fval(1,ibox),fval(1,ik1),fval(1,ik2),
     1                 fval(1,ik3),fval(1,ik4))
c------------------------
         endif



c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 2 .and. 
     1        (icntr .eq. 1 .or. icntr .eq. 2))then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c          after calling subdivide or subdivide1, assign fval
           ik1=ichildbox(1,ibox)
           ik2=ichildbox(2,ibox)
           ik3=ichildbox(3,ibox)
           ik4=ichildbox(4,ibox)
           call chebyp2k(ndeg,fval(1,ibox),fval(1,ik1),fval(1,ik2),
     1                 fval(1,ik3),fval(1,ik4))
c------------------------

         endif



c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 3 .and. icntr .eq. 1)then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c           after calling subdivide or subdivide1, assign fval
            ik1=ichildbox(1,ibox)
            ik2=ichildbox(2,ibox)
            ik3=ichildbox(3,ibox)
            ik4=ichildbox(4,ibox)
            call chebyp2k(ndeg,fval(1,ibox),fval(1,ik1),fval(1,ik2),
     1                 fval(1,ik3),fval(1,ik4))
c------------------------

         endif



c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 4 .and. 
     1        (icntr .eq. 2 .or. icntr .eq. 3))then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c           after calling subdivide or subdivide1, assign fval
            ik1=ichildbox(1,ibox)
            ik2=ichildbox(2,ibox)
            ik3=ichildbox(3,ibox)
            ik4=ichildbox(4,ibox)
            call chebyp2k(ndeg,fval(1,ibox),fval(1,ik1),fval(1,ik2),
     1                 fval(1,ik3),fval(1,ik4))
c------------------------

         endif


c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 6 .and. 
     1        (icntr .eq. 1 .or. icntr .eq. 4))then
             
             
c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c           after calling subdivide or subdivide1, assign fval
            ik1=ichildbox(1,ibox)
            ik2=ichildbox(2,ibox)
            ik3=ichildbox(3,ibox)
            ik4=ichildbox(4,ibox)
            call chebyp2k(ndeg,fval(1,ibox),fval(1,ik1),fval(1,ik2),
     1                 fval(1,ik3),fval(1,ik4))
c------------------------

         endif


c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 7 .and. icntr .eq. 3)then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c           after calling subdivide or subdivide1, assign fval
            ik1=ichildbox(1,ibox)
            ik2=ichildbox(2,ibox)
            ik3=ichildbox(3,ibox)
            ik4=ichildbox(4,ibox)
            call chebyp2k(ndeg,fval(1,ibox),fval(1,ik1),fval(1,ik2),
     1                 fval(1,ik3),fval(1,ik4))
c------------------------

         endif


c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 8 .and. 
     1        (icntr .eq. 3 .or. icntr .eq. 4))then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c           after calling subdivide or subdivide1, assign fval
            ik1=ichildbox(1,ibox)
            ik2=ichildbox(2,ibox)
            ik3=ichildbox(3,ibox)
            ik4=ichildbox(4,ibox)
            call chebyp2k(ndeg,fval(1,ibox),fval(1,ik1),fval(1,ik2),
     1                 fval(1,ik3),fval(1,ik4))
c------------------------

         endif


c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 9 .and. icntr .eq. 4)then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c           after calling subdivide or subdivide1, assign fval
            ik1=ichildbox(1,ibox)
            ik2=ichildbox(2,ibox)
            ik3=ichildbox(3,ibox)
            ik4=ichildbox(4,ibox)
            call chebyp2k(ndeg,fval(1,ibox),fval(1,ik1),fval(1,ik2),
     1                 fval(1,ik3),fval(1,ik4))
c------------------------

         endif


c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3
           endif
300     continue   
400     continue
500     continue
      end do
      return
      end subroutine





c----------------------------------------------------------
c     subroutine chebyk2p
c----------------------------------------------------------
c
c     This subroutine generates the k*k Chebyshev coefficients on a 
c     parent box, given function values on the k*k grids of children
c     boxes. (to be used to determine the resolution of the function
c     on the parent box)
c
c     INPUT:
c     ndeg: order of Chebyshev approx
c     valk1: function values on kid 1
c     valk2: function values on kid 2
c     valk3: function values on kid 3
c     valk4: function values on kid 4
c
c     OUTPUT:
c     valp: function values on parent box
c     coeffp: chebyshev coefficients on parent box
c
c     (remember to save the function values on parent box, when
c      doing a merging)
c
c------------------------------------------------------------------
c
      subroutine chebyk2p(ndeg, valk1, valk2, valk3, valk4, 
     1           valp, coeffp)
      implicit real*8 (a-h,o-z)
      integer ndeg
      real*8 valk1(ndeg,ndeg), valk2(ndeg,ndeg)
      real*8 valk3(ndeg,ndeg), valk4(ndeg,ndeg)
      real*8 valp(ndeg,ndeg), coeffp(ndeg,ndeg)
c     local vars
      real*8 wsave(1000), work(1000)
      real*8 xf(ndeg), yf(ndeg)
      real*8 coeffk1(ndeg,ndeg), coeffk2(ndeg,ndeg)
      real*8 coeffk3(ndeg,ndeg), coeffk4(ndeg,ndeg)
      real*8 cent0(2)
c
      cent0(1)=0.0d0
      cent0(2)=0.0d0
      call chxcin(ndeg, wsave)
c
      call chexfc2d(valk1,ndeg,coeffk1,wsave,work)
      call chexfc2d(valk2,ndeg,coeffk2,wsave,work)
      call chexfc2d(valk3,ndeg,coeffk3,wsave,work)
      call chexfc2d(valk4,ndeg,coeffk4,wsave,work)
c          got the chebyshev coefficients of each kid
c
      call mkgrid(xf, yf, 1, 1, 0, ndeg, cent0, 1.0d0)
c          pretend that the parent box is the root box
c          (since scale doesn't matter here)
c
      do k2=1,ndeg
        do k1=1,ndeg
          x=xf(k1)
          y=yf(k2)
          if((x.gt.0.0d0) .and. (y.gt.0.0d0)) then
            a=0.0d0
            b=0.5d0
            c=0.0d0
            d=0.5d0
            call cheval2d(x,y,valp(k1,k2),coeffk2,ndeg,a,b,c,d)
          elseif((x.lt.0.0d0) .and. (y.gt.0.0d0)) then
            a=-0.5d0
            b=0.0d0
            c=0.0d0
            d=0.5d0
            call cheval2d(x,y,valp(k1,k2),coeffk1,ndeg,a,b,c,d)
          elseif((x.lt.0.0d0) .and. (y.lt.0.0d0)) then
            a=-0.5d0
            b=0.0d0
            c=-0.5d0
            d=0.0d0
            call cheval2d(x,y,valp(k1,k2),coeffk4,ndeg,a,b,c,d)
          elseif((x.gt.0.0d0) .and. (y.lt.0.0d0)) then
            a=0.0d0
            b=0.5d0
            c=-0.5d0
            d=0.0d0
            call cheval2d(x,y,valp(k1,k2),coeffk3,ndeg,a,b,c,d)
          endif
        enddo
      enddo
c       interpolated to the k*k grid on the parent box
c
      call chexfc2d(valp,ndeg,coeffp,wsave,work) 
c

      end subroutine





c--------------------------------------------------------
c     subroutine chebyp2k
c--------------------------------------------------------
c
c     This subroutine interpolates function values on a 
c     k*k chebyshev grid to the k*k chebyshev grids of
c     the four children boxes 
c
c     INPUT:
c     ndeg: order of Chebyshev approx
c     valp: function values on parent box
c
c     OUTPUT:
c     valk1: function values on kid 1
c     valk2: function values on kid 2
c     valk3: function values on kid 3
c     valk4: function values on kid 4
c
c     (remember to update the function values on children
c      boxes when subdividing, to maintain the function
c      value array that makes sense for the current tree)
c
c--------------------------------------------------------
c
      subroutine chebyp2k(ndeg,valp,valk1,valk2,valk3,valk4)
      implicit real*8 (a-h,o-z)
      integer ndeg
      real*8 valp(ndeg,ndeg)
      real*8 valk1(ndeg,ndeg), valk2(ndeg,ndeg)
      real*8 valk3(ndeg,ndeg), valk4(ndeg,ndeg)
c     local vars
      real*8 wsave(1000), work(1000), valtmp(ndeg,ndeg)
      real*8 xf(ndeg), yf(ndeg), coeffp(ndeg,ndeg)
      real*8 cent0(2)
c
      cent0(1)=0.0d0
      cent0(2)=0.0d0
c
      do k2=1,ndeg
      do k1=1,ndeg
        valtmp(k1,k2)=valp(k1,k2)
      enddo
      enddo
c       not sure if the following calls would destroy
c       the original values, so I copy it over to be safe
c
      call chxcin(ndeg, wsave)
      call chexfc2d(valtmp,ndeg,coeffp,wsave,work)
c          got the chebyshev coefficients of the 
c          parent box. now it remains to evaluate
c          the chebyshev series on the grids of children
c
      a=-0.5d0
      b=0.5d0
      c=-0.5d0
      d=0.5d0
c     
      call mkgrid(xf, yf, 1, 2, 1, ndeg, cent0, 1.0d0)
      do k2=1,ndeg
      do k1=1,ndeg
        x=xf(k1)
        y=yf(k2)
        k=(k2-1)*ndeg+k1
        call cheval2d(x,y,valk1(k1,k2),coeffp,ndeg,a,b,c,d)
      enddo
      enddo
c          child 1
c
      call mkgrid(xf, yf, 2, 2, 1, ndeg, cent0, 1.0d0)
      do k2=1,ndeg
      do k1=1,ndeg
        x=xf(k1)
        y=yf(k2)
        k=(k2-1)*ndeg+k1
        call cheval2d(x,y,valk2(k1,k2),coeffp,ndeg,a,b,c,d)
      enddo
      enddo
c          child 2
c
      call mkgrid(xf, yf, 2, 1, 1, ndeg, cent0, 1.0d0)
      do k2=1,ndeg
      do k1=1,ndeg
        x=xf(k1)
        y=yf(k2)
        k=(k2-1)*ndeg+k1
        call cheval2d(x,y,valk3(k1,k2),coeffp,ndeg,a,b,c,d)
      enddo
      enddo
c          child 3
c
      call mkgrid(xf, yf, 1, 1, 1, ndeg, cent0, 1.0d0)
      do k2=1,ndeg
      do k1=1,ndeg
        x=xf(k1)
        y=yf(k2)
        k=(k2-1)*ndeg+k1
        call cheval2d(x,y,valk4(k1,k2),coeffp,ndeg,a,b,c,d)
      enddo
      enddo
c          child 4

      end subroutine




c**********************************************************
c
c     the following subroutine is used to set up the arrays of parents
c     and children for each box.  if a box is childless, all of its
c     children are set to -1.
c     each box that has children has four.  they are numbered
c     as follows:
c
c                        1     2
c
c                        4     3
c
c
c     this is done merely by scanning through all of the boxes and
c     comparing the column and row numbers of boxes that are one level
c     apart.
c
c     input:
c
c     levelbox is an array determining the level of each box
c
c     iparentbox denotes the parent of each box
c
c     icolbox denotes the column of each box
c
c     irowbox denotes the row of each box
c
c     nboxes is the total number of boxes
c
c     nlev is the finest level
c
c     output:
c
c     ichildbox denotes the four children of each box
c
c**********************************************************
c
      subroutine mkchild(levelbox,iparentbox,ichildbox,icolbox,
     1         irowbox,nboxes,nlev)
      implicit none
c-----global variables
      integer *4 levelbox(1)
      integer *4 iparentbox(1), ichildbox(4,1)
      integer *4 icolbox(1), irowbox(1)
      integer *4 nboxes, nlev
c-----local variables
      integer *4 i, j, l
      integer *4 icol, irow
      integer *4 icoltest, irowtest

c     first let's initialize all of the ichildbox
c     arrays to zero:
      do j = 1, nboxes
        iparentbox(j) = -1
        do i = 1, 4
          ichildbox(i,j) = -1
        end do
      end do

      do i = 0, nlev
        do 200 j = 1, nboxes
        if (levelbox(j) .ne. i)goto 200
c         figure out where the first child should be
c         by looking at the parents row and column.
          icoltest = 2*(icolbox(j) - 1) + 1
          irowtest = 2*(irowbox(j) - 1) + 1
c         now scan all of the boxes on the next
c         finest level.
          do 100 l = 1, nboxes
          if(levelbox(l) .ne. i+1)goto 100
            icol = icolbox(l)
            irow = irowbox(l)
            if(icol .eq. icoltest .and. irow .eq. irowtest)then
              ichildbox(4,j) = l
              iparentbox(l) = j
            elseif(icol .eq. icoltest + 1 
     1                        .and. irow .eq. irowtest + 1)then
              ichildbox(2,j) = l
              iparentbox(l) = j
            elseif(icol .eq. icoltest 
     1                        .and. irow .eq. irowtest + 1)then
              ichildbox(1,j) = l
              iparentbox(l) = j
            elseif(icol .eq. icoltest + 1 
     1                        .and. irow .eq. irowtest)then
              ichildbox(3,j) = l
              iparentbox(l) = j
            endif
100       continue
200     continue
      end do
      return
      end




c**********************************************************
c     subroutine mkcoeffs
c**********************************************************
c     the following subroutine defines the array coeffs that contains
c     the polynomial coefficients for the polynomial that approximates
c     the right hand side.
c
c     input:
c
c     ndeg is the order of approximation
c     fright is the right hand side defined on the old grid
c     levelbox is an array determining the level of each box
c     nboxes is the total number of boxes
c     nlev is the finest level
c     ichildbox denotes the four children of each box
c     nblevel is the total number of boxes per level
c     iboxlev is the array in which the boxes are arranged
c     istartlev is the pointer to where each level begins
c               in the iboxlev array
c
c     output:
c
c     coeffs is the array of coefficients for the basis functions
c
c**********************************************************
c
      subroutine mkcoeffs(ndeg,coeffs,fright,
     1         nlev, ichildbox, nblevel, 
     2         iboxlev, istartlev)
      implicit none
c-----global variables
      integer *4  ndeg,nlev
      integer *4  ichildbox(4,1)
      integer *4  nblevel(0:1), iboxlev(1), istartlev(0:1)
      real *8  coeffs(0:ndeg-1,0:ndeg-1,1), fright(ndeg**2,1)
c-----local variables
      integer *4  i, ibox, j, l, ii
      real *8  ftemp(ndeg,ndeg)
      real *8  coefftemp(ndeg,ndeg)
      real *8  wsave2(1000)

      do l = 0, nlev
       do i = istartlev(l), istartlev(l) + nblevel(l) - 1
       ibox = iboxlev(i)
        if (ichildbox(1,ibox) .lt.0)then


c       initialize wsave2 array (needed by getcoeff).
        call chxcin(ndeg,wsave2)

        do j = 1, ndeg
          do ii = 1, ndeg
            ftemp(j,ii) = fright(ndeg*(ii-1)+j,ibox)
          end do
        end do


c       compute chebyshev transforms
        call getcoeff(ndeg,ndeg,ftemp,coefftemp,wsave2)

c       now set these values to out coefficients
        do j = 0, ndeg-1
          do ii = 0, ndeg-1
            coeffs(j,ii,ibox) = coefftemp(j+1,ii+1)
          end do
        end do
        endif
       end do
      end do
      return
      end




c*****************************************************
c     subroutine treesort
c*****************************************************
c
c     Given the tree structure and a set of discrete
c     points, sort the points into boxes of the tree
c
c     INPUT:
c     nlev-istartlev: the tree structure, see
c          subroutine mktree for details
c     xpts: coordinates of the points
c     npts: number of points
c     
c     OUTPUT:
c     npbox: number of sources in each box 
c                     (including non-leaf boxes)
c     ipold: ids of points on input
c            (before the overwriting)
c     istartbox: starting point in the array xpts
c                for each box
c     ibofp: indices of leaf boxes that the point
c            belongs to (after overwriting)
c
c     (xpts: overwritten, reordered so that points in
c            a box are stored in adjacent places)
c
c*************************************************************
c
      subroutine treesort(nlev,levelbox,iparentbox,ichildbox,
     1           icolbox,irowbox,nboxes,nblevel,iboxlev,istartlev,
     2           xpts,npts,npbox,ipold,istartbox,ibofp,
     3           cent0,xsize0)
      implicit real*8 (a-h,o-z)
c     global vars
      integer nlev, nboxes, npts
      integer levelbox(nboxes), iparentbox(nboxes)
      integer ichildbox(4,nboxes)
      integer icolbox(nboxes), irowbox(nboxes)
      integer nblevel(0:nlev), iboxlev(nboxes)
      integer istartlev(0:nlev)
      integer npbox(nboxes), ipold(npts)
      integer istartbox(nboxes), ibofp(npts)
      real*8 xpts(2,npts), cent0(2), xsize0
c     local vars
      integer nsinchild(4)
      integer, allocatable:: isinchild(:,:)
      integer, allocatable:: itmp(:)
      real*8, allocatable:: xtmp(:,:)
c
c
      allocate(isinchild(npts,4))
      allocate(itmp(npts))
      allocate(xtmp(2,npts))
c
      do i=1,nboxes
        npbox(i)=0
        istartbox(i)=-1
      enddo 
      
      do j=1,npts
        ipold(j)=j
      enddo
c
      iroot=iboxlev(istartlev(0))

      npbox(iroot)=npts
      istartbox(iroot)=1
c
c    the root box contains all the points
c    RMK: in the current version, root box doesn't
c         necessarily have index 1
c         use iroot=iboxlev(istartlev(0))
c         instead of 1
c
c----------------------------------
c
      ix=1
      iy=1
c
      xlength=xsize0
      x0=cent0(1)-xsize0/2.0d0
      y0=cent0(2)-xsize0/2.0d0
      nside=1
c
      do l=0,nlev-1
c       go through all the levels
        xlength=xlength/2.0d0
        nside=nside*2

        istart=istartlev(l)
        iend=istartlev(l)+nblevel(l)-1

        do ii=istart,iend
          ibox=iboxlev(ii)
c         go through all non-empty boxes that have children
c
          if((npbox(ibox) .gt. 0) .and. 
     1               (ichildbox(1,ibox)).gt.0) then
c-------------------
            iss=istartbox(ibox)
            ies=istartbox(ibox)+npbox(ibox)-1
c
            do kk=1,4
              nsinchild(kk)=0
c             initialize the number of points in each
c             child to be zero
              ic=ichildbox(kk,ibox)
              npbox(ic)=0
            enddo
c
c             go through all the sources points in the box
c             (which occupy adjacent spots in issorted)      
c
            do j=iss,ies
              xx=xpts(1,j)
              yy=xpts(2,j) 
c
              xtmp(1,j)=xx
              xtmp(2,j)=yy
              itmp(j)=ipold(j)
c                 iss-th to ies-th entries of xpts
c                 and ipold copied to xtmp and itmp at
c                 corresponding places
              ix=ceiling((xx-x0)/xlength)
              iy=ceiling((yy-y0)/xlength)

              if(ix.le.0) then
                ix=1
              endif

              if(ix.gt.nside) then
                ix=nside
              endif

              if(iy.le.0) then
                iy=1
              endif

              if(iy.gt.nside) then
                iy=nside
              endif
c               x and y index of the point 
c               on the children's level
c
              do kk=1,4
                ic=ichildbox(kk,ibox)
                if((icolbox(ic).eq.ix).and.(irowbox(ic).eq.iy)) then
                  nsinchild(kk)=nsinchild(kk)+1
                  npbox(ic)=npbox(ic)+1
                  isinchild(nsinchild(kk),kk)=j
                endif
              enddo
            enddo
c
c           all points in box ibox sorted into isinchild,
c           indicating which child the points belong to. 
c           now overwrite the iss-th to ies-th spot in issorted
c
            jpt=iss
            do jc=1,4
              ic=ichildbox(jc,ibox)
              istartbox(ic)=jpt
              do js=1,nsinchild(jc)
                idnow=isinchild(js,jc)
                xpts(1,jpt)=xtmp(1,idnow)  
                xpts(2,jpt)=xtmp(2,idnow)
                ipold(jpt)=itmp(idnow)
                jpt=jpt+1
              enddo
            enddo
c-------------------
          endif
        enddo
c       go to the next level
      enddo

c     now get information for ibofp
      do l=0,nlev
        istart=istartlev(l)
        iend=istartlev(l)+nblevel(l)-1
        do ii=istart,iend
          ibox=iboxlev(ii)
          if((ichildbox(1,ibox).lt.0).and.(npbox(ibox).gt.0)) then
            iss=istartbox(ibox)
            ies=iss+npbox(ibox)-1
            do kk=iss,ies
              ibofp(kk)=ibox
            enddo
          endif
        enddo
      enddo
c
c
      deallocate(isinchild)
      deallocate(itmp)
      deallocate(xtmp)

      return
      end subroutine




c*****************************************************
c     subroutine refinetree
c*****************************************************
c
c     Given a discrete set of points and a quad-tree,
c     this subroutine refines the tree until there's no
c     more than a given number of points in each leaf
c     box.Points are reordered on output. The original
c     id is saved.
c
c     INPUT:
c     maxboxes: the maximum number of boxes allowed
c     maxlevel: deepest level allowed
c     maxppl: the maximum number of boundary points
c             allowed in leaf box
c     npts: number of points
c     xpts: coordinates of points
c     itemparray: a dummy array to be used in 
c             subroutine subdivide1, which is called within
c     cent0: center of the root box
c     xsize0: side length of the root box 
c     
c     (tree structure:)
c     istartlev is the pointer to where each level
c               begins in the iboxlev array
c     levelbox is an array determining the level of each box
c     nboxes is the total number of boxes
c     nlev is the finest level
c     icolbox denotes the column of each box
c     irowbox denotes the row of each box
c     iparentbox denotes the parent of each box
c     ichildbox denotes the four children of each box
c     nblevel is the total number of boxes per level
c     iboxlev is the array in which the boxes are arranged
c
c     (tree-point relation:)
c     npbox: number of sources in each box 
c                     (including non-leaf boxes)
c     ipold: ids of points on input
c            (before the overwriting)
c     istartbox: starting point in the array xpts
c                for each box
c     ibofp: indices of leaf boxes that the point
c            belongs to (after overwriting)
c
c     (xpts: overwritten, reordered so that points in
c            a box are stored in adjacent places)
c
c     OUTPUT:
c     Tree structure and point-tree relation are
c     rewritten to reflect the changes
c
c******************************************************
c
      subroutine refinetree(levelbox,icolbox,irowbox,
     1           nboxes,nlev,iparentbox,ichildbox,nblevel,
     2           iboxlev,istartlev,npbox,ipold,istartbox,
     3           ibofp,maxboxes,itemparray,maxlevel,
     4           maxppl,npts,xpts,cent0,xsize0)
      implicit real*8 (a-h,o-z)
c-----global variables
      integer nboxes, nlev, maxboxes, maxlevel
      integer maxppl, npts
      integer levelbox(1)
      integer icolbox(1), irowbox(1)
      integer iparentbox(1), ichildbox(4,1)
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer itemparray(1)
      integer npbox(1), ipold(1), istartbox(1), ibofp(1)
      real*8 xpts(2,npts)
      real*8 cent0(2), xsize0
c-----local variables
      integer nsinchild(4)
      integer, allocatable:: isinchild(:,:)
      integer, allocatable:: itmp(:)
      real*8, allocatable:: xtmp(:,:)
c
c-----------------------------------------------------
c
      allocate(isinchild(npts,4))
      allocate(itmp(npts))
      allocate(xtmp(2,npts))
c
      ix=1
      iy=1
      xlength=xsize0
      x0=cent0(1)-xsize0/2.0d0
      y0=cent0(2)-xsize0/2.0d0
      nside=1
c
      do i=0,maxlevel-1
        iflag=0
c           iflag: if the last level has been divided at all
        xlength=xlength/2.0d0
        nside=nside*2
c
        istart=istartlev(i)
        iend=istart+nblevel(i)-1
        do j=istart, iend
          ibox=iboxlev(j)
          if((ichildbox(1,ibox) .lt. 0).and.
     1       (npbox(ibox) .gt. maxppl)) then
c--------------
c           call subdivide1 to subdivide the current box
            call subdivide1(ibox,iparentbox,ichildbox,
     1           nboxes,irowbox,icolbox,levelbox,nlev,
     2           istartlev, nblevel, iboxlev, itemparray)
c
c           sort points in ibox into its newborn babies
c           attention: this time physically rewrite
c           the array xpts so that points in the same
c           box are stored in adjacent places
c
            iss=istartbox(ibox)
            ies=istartbox(ibox)+npbox(ibox)-1
c
            do kk=1,4
              nsinchild(kk)=0
c             initialize the number of points in each kid
c             to be zero
              ic=ichildbox(kk,ibox)
              npbox(ic)=0
            enddo
c
            do jj=iss,ies
              xx=xpts(1,jj)              
              yy=xpts(2,jj)              
c
              xtmp(1,jj)=xx
              xtmp(2,jj)=yy
              itmp(jj)=ipold(jj)
c                 iss-th to ies-th entries of xpts
c                 and ipold copied to xtmp and itmp at
c                 corresponding places
c
              ix=ceiling((xx-x0)/xlength)
              iy=ceiling((yy-y0)/xlength)

              if(ix.le.0) then
                ix=1
              endif

              if(ix.gt.nside) then
                ix=nside
              endif

              if(iy.le.0) then
                iy=1
              endif

              if(iy.gt.nside) then
                iy=nside
              endif
c               x and y index of the point 
c               on the children's level
c
              do kk=1,4
                ic=ichildbox(kk,ibox)
                if((icolbox(ic).eq.ix).and.(irowbox(ic).eq.iy)) then
                  nsinchild(kk)=nsinchild(kk)+1
                  npbox(ic)=npbox(ic)+1
                  isinchild(nsinchild(kk),kk)=jj
                endif
              enddo
            enddo
cccccc------
c
c         all points in box ibox sorted into isinchild,
c         indicating which child the points belong to. 
c         now overwrite the iss-th to ies-th spots in
c         xpts and ipold together
c
cccccc------
            jpt=iss
            do jc=1,4
              ic=ichildbox(jc,ibox)
              istartbox(ic)=jpt
              do js=1,nsinchild(jc)
                idnow=isinchild(js,jc)
                xpts(1,jpt)=xtmp(1,idnow)  
                xpts(2,jpt)=xtmp(2,idnow)
                ipold(jpt)=itmp(idnow)
                jpt=jpt+1
              enddo
            enddo
c           done with box: ibox
c            write(*,*) ibox
c--------------
            iflag=1
          endif      
        enddo 

        if((iflag .eq. 0).and.(i.eq.nlev))then
c         nothing was divides at the last level,
c         exit the loop (on levels)
          goto 456
        endif
      enddo
c
c     now get information for ibofp
456   do l=0,nlev
        istart=istartlev(l)
        iend=istartlev(l)+nblevel(l)-1
        do ii=istart,iend
          ibox=iboxlev(ii)
          if((ichildbox(1,ibox).lt.0).and.(npbox(ibox).gt.0)) then
            iss=istartbox(ibox)
            ies=iss+npbox(ibox)-1
            do kk=iss,ies
              ibofp(kk)=ibox
            enddo
          endif
        enddo
      enddo
c
      deallocate(isinchild)
      deallocate(itmp)
      deallocate(xtmp)

      return 
      end subroutine





c**********************************************************
c     subroutine subdivide1
c**********************************************************
c
c     this routine is identical to the subdivide routine except that
c     it does not concern itself at all with generating colleagues.
c     this routine is used only within the mktree8 routine.
c
c     input:
c
c     iparbox denotes the box being divided
c
c     iparentbox denotes the parent of each box
c
c     ichildbox denotes the four children of each box
c
c     nboxes is the total number of boxes
c
c     irowbox denotes the row of each box
c
c     icolbox denotes the column of each box
c
c     levelbox is an array determining the level of each box
c
c     nlev is the finest level of the current tree
c
c     istartlev is the pointer to where each level begins in the
c               iboxlev array
c
c     nblevel is the total number of boxes per level
c
c     iboxlev is the array in which the boxes are arranged
c
c     itemparray is just a dummy array
c
c     output:
c
c     The tree structure altered to reflect the change
c
c**********************************************************
c
      subroutine subdivide1(iparbox,iparentbox,ichildbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev,itemparray)
      implicit real*8 (a-h,o-z)
c-----global variables
      integer  nboxes, nlev
      integer  levelbox(1)
      integer  icolbox(1), irowbox(1)
      integer  ichildbox(4,1), iparentbox(1)
      integer  nblevel(0:nlev)
      integer  iboxlev(1),istartlev(0:nlev)
      integer itemparray(1)
c-----local variables
c
c     initialize the temp array itemparray to zero
      do i = 1, nboxes + 4
        itemparray(i) = 0
      end do
c
c     level, icolumn, and irow refer to the level, column,
c     and row of the parent box, respectively.
      level   = levelbox(iparbox)
      icolumn = icolbox(iparbox)
      irow    = irowbox(iparbox)
c
c     here are the new boxes placed in the
c     correct positions.  they are all childless.
c     there columns and rows are determined from
c     the parents columns and rows.  the level is
c     obviously one level finer than the parent.
      ibox = nboxes + 1
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 1
      irowbox(ibox) = 2*(irow-1) + 1

      ibox = nboxes + 2
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 2
      irowbox(ibox) = 2*(irow-1) + 1

      ibox = nboxes + 3
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 1
      irowbox(ibox) = 2*(irow-1) + 2

      ibox = nboxes + 4
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 2
      irowbox(ibox) = 2*(irow-1) + 2

      ichildbox(1,iparbox) = nboxes + 3
      ichildbox(2,iparbox) = nboxes + 4
      ichildbox(3,iparbox) = nboxes + 2
      ichildbox(4,iparbox) = nboxes + 1

c     set up a temporary array to store the old one in:
      do i = 1, nboxes
        itemparray(i) = iboxlev(i)
      end do

c     update the level info:
      if(level .eq. nlev) then
c      subdividing a box on the finest level
        istartlev(level+1)=nboxes+1
        iboxlev(istartlev(level+1))  =nboxes+1
        iboxlev(istartlev(level+1)+1)=nboxes+2
        iboxlev(istartlev(level+1)+2)=nboxes+3
        iboxlev(istartlev(level+1)+3)=nboxes+4

        nblevel(level+1)=4
        nlev=nlev+1
      elseif(level .eq. nlev-1) then
c      subdividing a box on the 2nd fines level
c      (i.e. the new boxes are on the fines level)
c      (istartlev stays the same)
        iboxlev(nboxes+1)=nboxes+1
        iboxlev(nboxes+2)=nboxes+2
        iboxlev(nboxes+3)=nboxes+3
        iboxlev(nboxes+4)=nboxes+4
c
        nblevel(level+1)=nblevel(level+1)+4
c       (nlev stays the same)
      else
c      subdividing a box not on the finest two levels
c      (level+2 exists)
c       put the new boxes in the right place in iboxlev
        iboxlev(istartlev(level+2))   = nboxes + 1
        iboxlev(istartlev(level+2)+1) = nboxes + 2
        iboxlev(istartlev(level+2)+2) = nboxes + 3
        iboxlev(istartlev(level+2)+3) = nboxes + 4
c
c       copy over the old iboxlev info
        do i = istartlev(level+2) + 4, nboxes + 4
          iboxlev(i) = itemparray(i - 4)
        end do 
c
c       update istartlev for finer levels
        do i = level + 2, nlev
          istartlev(i) = istartlev(i) + 4
        end do
c       
c       update nblevel
        nblevel(level+1)=nblevel(level+1)+4
      endif

      nboxes=nboxes+4
     
      end subroutine





c**********************************************************
c     subroutine subdivide
c**********************************************************
c
c     the following subroutine is designed to divide up a childless
c     box into four children.  
c     the children are placed in correct order (clockwise starting
c     from the upper left corner) so there is no need to 'shuffle' 
c     the child order later on. in the periodic case, the colleagues
c     must be obtained by looking at the potential colleague numbers
c     and their row and column and seeing if they lie outside of 
c     the domain. if they do it must be readjusted to account for the 
c     periodicity.
c
c     input:
c
c     iparbox denotes the box being divided
c
c     iparentbox denotes the parent of each box
c
c     ichildbox denotes the four children of each box
c
c     icolleagbox denotes the colleagues of a given box
c
c     nboxes is the total number of boxes
c
c     irowbox denotes the row of each box
c
c     icolbox denotes the column of each box
c
c     levelbox is an array determining the level of each box
c
c     nlev is the finest level of the current tree
c
c     istartlev is the pointer to where each level begins in the
c               iboxlev array
c
c     nblevel is the total number of boxes per level
c
c     iboxlev is the array in which the boxes are arranged
c
c     iperiod denotes what kind of colleagues are to be generated
c             iperiod = 1 : free space
c             iperiod = 2 : periodic
c             iperiod = 3 : periodic up/down and free space left/right
c
c     itemparray is just a dummy array
c
c     output:
c     
c     The tree structure altered to reflect the change
c
c**********************************************************
c
      subroutine subdivide(iparbox,iparentbox,ichildbox,
     1           icolleagbox,nboxes,irowbox,icolbox,
     2           levelbox,nlev,istartlev, nblevel, iboxlev,
     3           iperiod,itemparray)
      implicit real*8 (a-h,o-z)
c-----global variables
      integer *4  iparentbox(1), ichildbox(4,1)
      integer *4  icolbox(1), irowbox(1)
      integer *4  levelbox(1), nboxes
      integer *4  iparbox, nlev, iperiod
      integer *4  nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer *4  icolleagbox(9,1)
c-----local variables
      integer *4  level, ibox, itemparray(1)
      integer *4  icolumn, irow, i, j
      integer *4  icntr, jcntr, isister
      integer *4  icoltemp, icoltest
      integer *4  irowtemp, irowtest
      integer *4  partemp, colleague
      integer *4  nside, ilev, itest, l

c     let's initialize the array itemparray to zero:
      do i = 1, nboxes + 4
        itemparray(i) = 0
      end do

c     level, icolumn, and irow refer to the level, column,
c     and row of the parent box, respectively.
      level   = levelbox(iparbox)
      icolumn = icolbox(iparbox)
      irow    = irowbox(iparbox)

c     here are the new boxes placed in the
c     correct positions.  they are all childless.
c     there columns and rows are determined from
c     the parents columns and rows.  the level is
c     obviously one level finer than the parent.
      ibox = nboxes + 1
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 1
      irowbox(ibox) = 2*(irow-1) + 1

      ibox = nboxes + 2
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 2
      irowbox(ibox) = 2*(irow-1) + 1

      ibox = nboxes + 3
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 1
      irowbox(ibox) = 2*(irow-1) + 2

      ibox = nboxes + 4
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 2
      irowbox(ibox) = 2*(irow-1) + 2

      ichildbox(1,iparbox) = nboxes + 3
      ichildbox(2,iparbox) = nboxes + 4
      ichildbox(3,iparbox) = nboxes + 2
      ichildbox(4,iparbox) = nboxes + 1

c     set up a temporary array to store the old one in:
      do i = 1, nboxes
        itemparray(i) = iboxlev(i)
      end do

c     update the level info:
      if(level .eq. nlev) then
c      subdividing a box on the finest level
        istartlev(level+1)=nboxes+1
        iboxlev(istartlev(level+1))  =nboxes+1
        iboxlev(istartlev(level+1)+1)=nboxes+2
        iboxlev(istartlev(level+1)+2)=nboxes+3
        iboxlev(istartlev(level+1)+3)=nboxes+4

        nblevel(level+1)=4
        nlev=nlev+1
      elseif(level .eq. nlev-1) then
c      subdividing a box on the 2nd fines level
c      (i.e. the new boxes are on the fines level)
c      (istartlev stays the same)
        iboxlev(nboxes+1)=nboxes+1
        iboxlev(nboxes+2)=nboxes+2
        iboxlev(nboxes+3)=nboxes+3
        iboxlev(nboxes+4)=nboxes+4
c
        nblevel(level+1)=nblevel(level+1)+4
c       (nlev stays the same)
      else
c      subdividing a box not on the finest two levels
c      (level+2 exists)
c       put the new boxes in the right place in iboxlev
        iboxlev(istartlev(level+2))   = nboxes + 1
        iboxlev(istartlev(level+2)+1) = nboxes + 2
        iboxlev(istartlev(level+2)+2) = nboxes + 3
        iboxlev(istartlev(level+2)+3) = nboxes + 4
c
c       copy over the old iboxlev info
        do i = istartlev(level+2) + 4, nboxes + 4
          iboxlev(i) = itemparray(i - 4)
        end do 
c
c       update istartlev for finer levels
        do i = level + 2, nlev
          istartlev(i) = istartlev(i) + 4
        end do
c       
c       update nblevel
        nblevel(level+1)=nblevel(level+1)+4
      endif

      nboxes=nboxes+4

c --------------------------------
c     don't change anything below
c --------------------------------
c      now let's go through the process of reforming any 
c      necessary colleagues.  for each of the child boxes 
c      that we just formed, all we need to do is scan through
c      the boxes that are children of the above parent boxes 
c      colleagues and test the column and row numbers.  we can 
c      also take advantage of the fact that for every one of 
c      the newly formed boxes colleagues, that box will list 
c      the newly formed box as one of its colleagues.  
c      the colleague numbers can be found easily if we think 
c      of a 'reflection.'  colleague 1 and 9 are opposites, 
c      3 and 7 etc.
c      first do the free space case:
       if(iperiod .eq. 0)then
       do 200 i = 1, 4
         if(ichildbox(1,iparbox) .lt. 0)goto 200
         ibox = ichildbox(i,iparbox)
         icolleagbox(5,ibox) = ibox
         do j = 1, 4
           icolleagbox(j,ibox) = -1
         end do
         do j = 6, 9
           icolleagbox(j,ibox) = -1
         end do

         partemp = iparentbox(ibox)

c        irowtemp and icoltemp denote the
c        row and column of the test box.
         irowtemp = irowbox(ibox)
         icoltemp = icolbox(ibox)

         do 100 jcntr = 1, 9
c          colleague denotes the colleague of the parent box.
           colleague = icolleagbox(jcntr,partemp)
c          if the colleague doesn't exist
c          or is childless, skip it:
           if (colleague .lt. 0)goto 100
           if (ichildbox(1,colleague) .lt. 0)goto 100
c          otherwise scan the four children:
           do icntr = 1, 4
             j = ichildbox(icntr,colleague)
c            irowtest and icoltest denote the row and column of
c            the box being compared to the test box.
             irowtest = irowbox(j)
             icoltest = icolbox(j)

             if(irowtemp .eq. irowtest+1)then
               if(icoltemp .eq. icoltest+1)then
                 icolleagbox(1,ibox) = j
                 icolleagbox(9,j) = ibox
               elseif(icoltemp .eq. icoltest)then
                 icolleagbox(2,ibox) = j
                 icolleagbox(8,j) = ibox
               elseif(icoltemp .eq. icoltest-1)then
                 icolleagbox(3,ibox) = j
                 icolleagbox(7,j) = ibox
               endif
             elseif(irowtemp .eq. irowtest)then
               if(icoltemp .eq. icoltest+1)then
                 icolleagbox(4,ibox) = j
                 icolleagbox(6,j) = ibox
               elseif(icoltemp .eq. icoltest-1)then
                 icolleagbox(6,ibox) = j
                 icolleagbox(4,j) = ibox
               endif
             elseif(irowtemp .eq. irowtest-1)then
               if(icoltemp .eq. icoltest+1)then
                 icolleagbox(7,ibox) = j
                 icolleagbox(3,j) = ibox
               elseif(icoltemp .eq. icoltest)then
                 icolleagbox(8,ibox) = j
                 icolleagbox(2,j) = ibox
               elseif(icoltemp .eq. icoltest-1)then
                 icolleagbox(9,ibox) = j
                 icolleagbox(1,j) = ibox
               endif
             endif
          end do
100      continue
200   continue
c     next do the periodic case:
      elseif(iperiod .eq. 1 .or. iperiod .eq. 2)then
       ilev = level + 1
       nside = 2**ilev

       do l = 1, 4
        ibox = ichildbox(l,iparbox)
c       initialize colleague number 5 to
c       yourself and all other colleagues to
c       -1.  -1 is the flag for the case when
c       the colleagues don't exist.  this is 
c       for all of the newly created boxes.
         icolleagbox(5,ibox) = ibox
         do j = 1, 4
           icolleagbox(j,ibox) = -1
         end do
         do j = 6, 9
           icolleagbox(j,ibox) = -1
         end do

c       irowtemp and icoltemp denote the row and column of
c       the test box.
        irowtemp = irowbox(ibox)
        icoltemp = icolbox(ibox)

c       irowtest and icoltest denote the row and column of
c       the box being compared to the test box.
             do 400 jcntr = 1, 9
               if(jcntr .eq. 5)goto 400
               if(jcntr .eq. 1)then
                 icoltest = icoltemp - 1
                 irowtest = irowtemp - 1
                 isister = 9
               elseif(jcntr .eq. 2)then 
                 icoltest = icoltemp
                 irowtest = irowtemp - 1
                 isister = 8
               elseif(jcntr .eq. 3)then
                 icoltest = icoltemp + 1
                 irowtest = irowtemp - 1
                 isister = 7
               elseif(jcntr .eq. 4)then
                 icoltest = icoltemp - 1
                 irowtest = irowtemp
                 isister = 6
               elseif(jcntr .eq. 6)then
                 icoltest = icoltemp + 1
                 irowtest = irowtemp
                 isister = 4
               elseif(jcntr .eq. 7)then
                 icoltest = icoltemp - 1
                 irowtest = irowtemp + 1
                 isister = 3
               elseif(jcntr .eq. 8)then
                 icoltest = icoltemp
                 irowtest = irowtemp + 1
                 isister = 2
               elseif(jcntr .eq. 9)then
                 icoltest = icoltemp + 1
                 irowtest = irowtemp + 1
                 isister = 1
               endif

c         now test to see if the test parameters
c         lie in the domain.
            if(icoltest .lt. 1)then
               icoltest = icoltest + nside
            elseif(icoltest .gt. nside)then
               icoltest = icoltest - nside
            endif
            if(irowtest .lt. 1)then
               irowtest = irowtest + nside
            elseif(irowtest .gt. nside)then
               irowtest = irowtest - nside
            endif


       do 300 j = 1, 9
        if(icolleagbox(j,iparentbox(ibox)) .lt. 0)goto 300
        if(ichildbox(1,icolleagbox(j,iparentbox(ibox))) .lt. 0)goto 300
          do icntr = 1, 4
            itest = ichildbox(icntr,icolleagbox(j,iparentbox(ibox)))
            if(irowbox(itest) .eq. irowtest
     1         .and. icolbox(itest) .eq. icoltest)then
               icolleagbox(jcntr,ibox) = itest
               icolleagbox(isister,itest) = ibox
            endif
          end do
300    continue
400    continue
      end do
      elseif(iperiod .eq. 3)then
       ilev = level + 1
       nside = 2**ilev

       do l = 1, 4
        ibox = ichildbox(l,iparbox)
c       initialize colleague number 5 to
c       yourself and all other colleagues to
c       -1.  -1 is the flag for the case when
c       the colleagues don't exist.  this is 
c       for all of the newly created boxes.
         icolleagbox(5,ibox) = ibox
         do j = 1, 4
           icolleagbox(j,ibox) = -1
         end do
         do j = 6, 9
           icolleagbox(j,ibox) = -1
         end do

c       irowtemp and icoltemp denote the row and column of
c       the test box.
        irowtemp = irowbox(ibox)
        icoltemp = icolbox(ibox)

c       irowtest and icoltest denote the row and column of
c       the box being compared to the test box.
             do 600 jcntr = 1, 9
               if(jcntr .eq. 5)goto 600
               if(jcntr .eq. 1)then
                 icoltest = icoltemp - 1
                 irowtest = irowtemp - 1
                 isister = 9
               elseif(jcntr .eq. 2)then 
                 icoltest = icoltemp
                 irowtest = irowtemp - 1
                 isister = 8
               elseif(jcntr .eq. 3)then
                 icoltest = icoltemp + 1
                 irowtest = irowtemp - 1
                 isister = 7
               elseif(jcntr .eq. 4)then
                 icoltest = icoltemp - 1
                 irowtest = irowtemp
                 isister = 6
               elseif(jcntr .eq. 6)then
                 icoltest = icoltemp + 1
                 irowtest = irowtemp
                 isister = 4
               elseif(jcntr .eq. 7)then
                 icoltest = icoltemp - 1
                 irowtest = irowtemp + 1
                 isister = 3
               elseif(jcntr .eq. 8)then
                 icoltest = icoltemp
                 irowtest = irowtemp + 1
                 isister = 2
               elseif(jcntr .eq. 9)then
                 icoltest = icoltemp + 1
                 irowtest = irowtemp + 1
                 isister = 1
               endif

c         now test to see if the test parameters
c         lie in the domain.
            if(irowtest .lt. 1)then
               irowtest = irowtest + nside
            elseif(irowtest .gt. nside)then
               irowtest = irowtest - nside
            endif


       do 500 j = 1, 9
        if(icolleagbox(j,iparentbox(ibox)) .lt. 0)goto 500
        if(ichildbox(1,icolleagbox(j,iparentbox(ibox))) .lt. 0)goto 500
          do icntr = 1, 4
            itest = ichildbox(icntr,icolleagbox(j,iparentbox(ibox)))
            if(irowbox(itest) .eq. irowtest
     1         .and. icolbox(itest) .eq. icoltest)then
               icolleagbox(jcntr,ibox) = itest
               icolleagbox(isister,itest) = ibox
            endif
          end do
500    continue
600    continue
      end do
      endif

      return
      end subroutine





c-----------------------------------------------------
c     subroutine coarsen1
c-----------------------------------------------------
c
c     this subroutine merges four kids of a given box
c     where all kids are leaf nodes
c     (Attention: this version doesn't alter icolleaguebox!)
c
c     INPUT:
c     iparbox denotes the box whose kids are being merged
c     iparentbox denotes the parent of each box
c     ichildbox denotes the four children of each box
c     nboxes is the total number of boxes
c     irowbox denotes the row of each box
c     icolbox denotes the column of each box
c     levelbox is an array determining the level of each box
c     nlev is the finest level
c     istartlev is the pointer to where each level begins in the
c               iboxlev array
c     nblevel is the total number of boxes per level
c     iboxlev is the array in which the boxes are arranged
c---------------------------
c     iempty: array of length nboxes, binary indicator:
c             iempty = 1: empty box
c             iempty = 0: not empty
c
c     OUTPUT:
c     iparentbox, ichildbox, and iempty are altered
c     to reflect the change
c
c-----------------------------------------------------
c
      subroutine coarsen1(iparbox,iparentbox,ichildbox,
     1           nboxes,irowbox,icolbox,levelbox,nlev,
     2           istartlev,nblevel,iboxlev,iempty)
      implicit real*8 (a-h,o-z)
c     global vars
      integer nboxes, nlev, iparbox
      integer iparentbox(1), ichildbox(4,1)
      integer irowbox(1), icolbox(1), levelbox(1)
      integer istartlev(0:1), nblevel(0:1), iboxlev(1)
      integer iempty(1)
c
c      write(*,*)"test if the subroutine coarsen1 run"
c

      kid1=ichildbox(1,iparbox)
      kid2=ichildbox(2,iparbox)
      kid3=ichildbox(3,iparbox)
      kid4=ichildbox(4,iparbox)
c
      ichildbox(1,iparbox)=-1
      ichildbox(2,iparbox)=-1
      ichildbox(3,iparbox)=-1
      ichildbox(4,iparbox)=-1
c
      iparentbox(kid1)=-1 
      iparentbox(kid2)=-1 
      iparentbox(kid3)=-1 
      iparentbox(kid4)=-1 
c
      iempty(kid1)=1
      iempty(kid2)=1
      iempty(kid3)=1
      iempty(kid4)=1

      end subroutine




c-----------------------------------------------------
c     subroutine reorder
c-----------------------------------------------------
c
c     This subroutine cleans up an intermediate tree, to get rid of
c     empty boxes and reassigns box numbers
c
c     INPUT:
c     iparentbox denotes the parent of each box
c     ichildbox denotes the four children of each box
c     nboxes is the total number of boxes
c     irowbox denotes the row of each box
c     icolbox denotes the column of each box
c     levelbox is an array determining the level of each box
c     nlev is the finest level
c     istartlev is the pointer to where each level begins in the
c               iboxlev array
c     nblevel is the total number of boxes per level
c     iboxlev is the array in which the boxes are arranged
c---------------------------
c     iempty: array of length nboxes, binary indicator:
c             iempty = 1: empty box
c             iempty = 0: not empty
c     itemp: temporary working array, length>= nboxes
c     itemp4: temporary working array, dimension(4,nboxes)
c
c     OUTPUT:
c     all the tree structure may be altered to reflect the change
c     idold2new: maps old box id to new box id
c      
c-----------------------------------------------------
c
      subroutine reorder(iparentbox, ichildbox, nboxes, 
     1           irowbox, icolbox, levelbox, nlev,
     2           istartlev, nblevel, iboxlev, iempty, 
     3           idold2new, nboxesold, itemp, itemp4)
      implicit real*8 (a-h,o-z)
c     global vars
      integer nboxes, nlev, iparbox
      integer iparentbox(1), ichildbox(4,1)
      integer irowbox(1), icolbox(1), levelbox(1)
      integer istartlev(0:1), nblevel(0:1), iboxlev(1)
      integer iempty(1), itemp(1), itemp4(4,1)
      integer idold2new(nboxes)
c
c------------------------------
c
c     First, go through the boxes, generate the new id for
c     each box, store them in array idold2new. at the same 
c     time, determine the actual number of boxes on each 
c     level, rewrite nblevel
c
      nshift=0
      do i=1,nboxes
        if(iempty(i) .eq. 0) then
c         box i is not empty
          idold2new(i)=i-nshift
        elseif(iempty(i) .eq. 1) then
c         box i is empty
          nshift=nshift+1
          idold2new(i)=-1 
          lev=levelbox(i)
          nblevel(lev)=nblevel(lev)-1                        
c           update the number of boxes on the current level
        endif
      enddo
      nempty=nshift
c      write(*,*) 'nempty=',nempty
c       got the new number of boxes
c       and updated nblevel
c
c---------------------------------
c
c     second, update the tree structure, using the new numbering      
c
c     1. update iparentbox
c
      do i=1,nboxes
        itemp(i)=iparentbox(i)
      enddo
c
      do iold=1, nboxes
        if(iempty(iold) .eq. 0) then
          jold=itemp(iold)
          i=idold2new(iold)
          if(jold .gt. 0) then
            j=idold2new(jold)
          else
            j=-1
          endif
          iparentbox(i)=j
        endif
      enddo      
c
c---------------------------------
c
c     2. update ichildbox
c
      do i=1,nboxes
        itemp4(1,i)=ichildbox(1,i)
        itemp4(2,i)=ichildbox(2,i)
        itemp4(3,i)=ichildbox(3,i)
        itemp4(4,i)=ichildbox(4,i)
      enddo
c
      do iold=1, nboxes
        if(iempty(iold) .eq. 0) then
          jold1=itemp4(1,iold)
          jold2=itemp4(2,iold)
          jold3=itemp4(3,iold)
          jold4=itemp4(4,iold)
c
          i=idold2new(iold) 
          if(jold1 .gt. 0) then
            j1=idold2new(jold1)
            j2=idold2new(jold2)
            j3=idold2new(jold3)
            j4=idold2new(jold4)
          else
            j1=-1
            j2=-1
            j3=-1
            j4=-1
          endif
c
          ichildbox(1,i)=j1
          ichildbox(2,i)=j2
          ichildbox(3,i)=j3
          ichildbox(4,i)=j4
        endif
      enddo 
c
c---------------------------------
c
c     3. update irowbox, icolbox, and levelbox
c
      do i=1,nboxes
        itemp(i)=irowbox(i)
      enddo  
c
      do iold=1,nboxes
        if(iempty(iold) .eq. 0) then
          i=idold2new(iold)  
          irowbox(i)=itemp(iold)
        endif
      enddo
c
c----------------------------
c
      do i=1,nboxes
        itemp(i)=icolbox(i)
      enddo 
c
      do iold=1,nboxes
        if(iempty(iold) .eq. 0) then
          i=idold2new(iold) 
          icolbox(i)=itemp(iold)
        endif
      enddo
c
c----------------------------
c
      do i=1,nboxes
        itemp(i)=levelbox(i)
      enddo
c
      do iold=1,nboxes
        if(iempty(iold) .eq. 0) then
          i=idold2new(iold)
          levelbox(i)=itemp(iold)
        endif
      enddo
c
c---------------------------------
c
c     4. update istartlev, nblevel, iboxlev 
c
c     4.1 use the new numbering of boxes in: iboxlev
c     
      do i=1,nboxes
        iold=iboxlev(i)
        if(iempty(iold) .eq. 0) then  
          inew=idold2new(iold)
          itemp(i)=inew
        else
          itemp(i)=-1
        endif
      enddo
c
c     4.2 update istartlev, nblevel, and iboxlev.
c         shift boxes in iboxlev first, and then
c         count istartlev and nblevel
c
      nshift=0
      do i=1, nboxes
        if(itemp(i) .gt. 0) then
          iboxlev(i-nshift)=itemp(i)
        else
          nshift=nshift+1
        endif
      enddo
c
c     now count and overwrite istartlev and nblevel      
c
cccccc
c      write(*,*) nboxes, nempty, nboxes-nempty
c      do ii=1,nboxes
c        write(*,*) ii, idold2new(ii)
c      enddo
cccccc
c
      nboxesold=nboxes
      nboxes=nboxes-nempty
      do i=1, nboxes
        iempty(i)=0
      enddo
c           now that almost everything has been overwritten
c           let's finally change nboxes
c
      levpre=0
      istartlev(0)=1
c
      do i=1,nboxes
        ibox=iboxlev(i)
        lev=levelbox(ibox)
        if(lev .ne. levpre) then
          istartlev(lev)=i
          levpre=lev
        endif 
      enddo
c       this is enough since the 'ladder' ordering is preserved
c       
c     one last thing, update nlev
c     (there could be levels at the end disappearing)
c
      nlevnew=nlev
      do l=0,nlev
        if(nblevel(l) .le. 0) then
          nlevnew=nlevnew-1 
        endif
      enddo
      nlev=nlevnew
c

      end subroutine




c-----------------------------------------------------
c     subroutine reorderf
c-----------------------------------------------------
c
c     Given an array fval, and a reordering of box ids
c     reorder the array accordingly
c
c-----------------------------------------------------
c
      subroutine reorderf(na, nold, fval, idold2new)
      implicit real*8 (a-h,o-z)
      integer na, nboxes
      integer idold2new(nold)
      real*8 fval(na, nold)
      real*8, allocatable:: ftemp(:,:)
c
      allocate(ftemp(na,nold))

c     save the old array in ftemp
      do j=1,nold
        do i=1,na 
          ftemp(i,j)=fval(i,j)
        enddo
      enddo
c
      do j=1,nold
        jnew=idold2new(j)
        if(jnew .gt. 0) then
          do i=1,na
            fval(i,jnew)=ftemp(i,j)
          enddo
        endif
      enddo

      deallocate(ftemp)

      end subroutine




c-----------------------------------------------------
c     subroutine adaptree
c-----------------------------------------------------
c
c     This subroutine deals with adaptive mesh refining 
c     and coarsening of a quad-tree, based on two criteria
c     1) the underlying function should be resolved
c     2) each leaf box contains no more than O(1) points 
c
c     This subroutine consists of two 'sweeps':
c     1) a downward sweep to refine (until the two
c        criteria are met)
c     2) an upward sweep to merge when over-resolved
c        (consider empty boxes only)
c
c     RMK: 1) level-restriction may be violated on output
c             call 'fixtree' after this 
c          2) all box ids may be changed due to coarsening
c          3) this subroutine is not optimal, but
c             sufficient for now. optimize when needed.
c             in particular, the 'merge' part is not 
c             stringent. (when an object moves across a
c             region, it may leave a bit of a tail, but
c             not for long)
c
c     INPUT:
c     maxboxes: the maximum number of boxes allowed
c     maxlevel: deepest level allowed
c     maxppl: the maximum number of boundary points allowed
c             in leaf box
c     levelbox - istartlev: the quad-tree structure
c     npts: number of points
c     xpts: coordinates of the points
c     npbox: number of sources in each box 
c                     (including non-leaf boxes)
c     ipold: ids of points on input
c            (before the overwriting)
c     istartbox: starting point in the array xpts
c                for each box
c     ibofp: indices of leaf boxes that the point
c            belongs to (after overwriting)
c     ndeg: degree of chebyshev approx 
c     fval: function values on the ndeg x ndeg grid of
c           the leaf nodes
c     eps: error tolerance of function approx
c     isort: indicator. If isort = 0, do not sort the 
c            sources after altering the tree (the 
c            relation might be wrong, but we choose to 
c            sort them outside since there might be more 
c            changes to the tree)
c            If isort =1, sort the sources
c
c     OUTPUT:
c     the tree overwritten to reflect the changes
c     the points may be reordered, and point-tree relation
c         may be overwritten
c     the function values interpolated to the new 
c         leaf nodes of the tree
c
c-----------------------------------------------------
c
      subroutine adaptree(maxboxes, maxlevel, maxppl,
     1           levelbox, icolbox, irowbox, nboxes,
     2           nlev, iparentbox, ichildbox, nblevel,
     3           iboxlev, istartlev, npts, xpts, 
     4           npbox, ipold, istartbox, ibofp,
     5           ndeg, fval, eps, isort, cent0, xsize0)
      implicit real*8 (a-h,o-z)
c-----global variables
      integer maxboxes, maxlevel, maxppl
      integer nboxes, nlev, npts, ndeg, isort
      integer levelbox(1)
      integer icolbox(1), irowbox(1)
      integer iparentbox(1), ichildbox(4,1)
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer npbox(1), ipold(1), istartbox(1), ibofp(1)
      real*8 xpts(2,1), fval(ndeg*ndeg,1), eps
      real*8 cent0(2), xsize0
c-----local variables
      integer kids(4), nsinchild(4)
      integer, allocatable:: iempty(:), itemparray(:)
      integer, allocatable:: itemp(:), itemp4(:,:)
      integer, allocatable:: isinchild(:,:), itmp(:)
      integer, allocatable:: idold2new(:)
      real*8, allocatable:: ftemp(:,:), coefftemp(:,:)
      real*8, allocatable:: xtmp(:,:)
      real*8, allocatable:: valk1(:,:), valk2(:,:)
      real*8, allocatable:: valk3(:,:), valk4(:,:)
      real*8 wsave(1000), chwork(1000)
c
      allocate(ftemp(ndeg,ndeg))
      allocate(coefftemp(ndeg,ndeg))
      allocate(itemparray(maxboxes))
c
      allocate(isinchild(npts,4))
      allocate(itmp(npts))
      allocate(xtmp(2,npts))
c
      allocate(iempty(maxboxes))
      allocate(itemp(maxboxes))
      allocate(itemp4(4,maxboxes))
      allocate(idold2new(maxboxes))
c
      allocate(valk1(ndeg,ndeg))
      allocate(valk2(ndeg,ndeg))
      allocate(valk3(ndeg,ndeg))
      allocate(valk4(ndeg,ndeg))
c
c---------------------------------------------------------
c     1) refinement sweep. go down the tree level by level
c        refine leaf boxes if one of the criteria is 
c        violated, sort sources along the way
c---------------------------------------------------------
c
      ix=1
      iy=1
c      xlength=1.0d0
      xlength=xsize0
      x0=cent0(1)-xsize0/2.0d0
      y0=cent0(2)-xsize0/2.0d0
c
      nside=1
c
      call chxcin(ndeg,wsave)
      do i=0, maxlevel-1
        iflag=0
c          iflag: if the last level has been divided at all
        xlength=xlength/2.0d0
        nside=nside*2
c
        istart=istartlev(i)
        iend=istart+nblevel(i)-1
        do j=istart, iend
          ibox=iboxlev(j)
          if(ichildbox(1,ibox) .lt. 0) then
            npb=npbox(ibox)
c           1) get the number of points in a leaf box 
c
            do k2=1,ndeg
            do k1=1,ndeg
              k=(k2-1)*ndeg+k1
              ftemp(k1,k2)=fval(k,ibox)
            enddo
            enddo
c
            call getcoeff(ndeg,ndeg,ftemp,coefftemp,wsave)
            call geterror(ndeg,coefftemp,error)
            hh = dble(4**levelbox(ibox))
            epscaled = eps * hh
c           2) check resolution of the underlying function
c
            
            if((error .ge. epscaled).or.(npb .gt. maxppl)) then
              call subdivide1(ibox,iparentbox,ichildbox,
     1           nboxes,irowbox,icolbox,levelbox, nlev,
     2           istartlev, nblevel, iboxlev, itemparray)
c              after calling subdivide1, remember to assign
c              interpolated function values to children
              ik1=ichildbox(1,ibox)
              ik2=ichildbox(2,ibox)
              ik3=ichildbox(3,ibox)
              ik4=ichildbox(4,ibox)
              call chebyp2k(ndeg,fval(1,ibox),fval(1,ik1),
     1             fval(1,ik2),fval(1,ik3),fval(1,ik4))
c             now sort points in ibox into its children
c             (as in 'refinetree')
c
              iss=istartbox(ibox)
              ies=istartbox(ibox)+npbox(ibox)-1
c
              do kk=1,4
                nsinchild(kk)=0
c               initialize the number of points in each kid
c               to be zero
                ic=ichildbox(kk,ibox)
                npbox(ic)=0
              enddo
c
              do jj=iss,ies
                xx=xpts(1,jj)              
                yy=xpts(2,jj)              
c
                xtmp(1,jj)=xx
                xtmp(2,jj)=yy
                itmp(jj)=ipold(jj)
c                 iss-th to ies-th entries of xpts
c                 and ipold copied to xtmp and itmp at
c                 corresponding places
c
                ix=ceiling((xx-x0)/xlength)
                iy=ceiling((yy-y0)/xlength)

                if(ix.le.0) then
                  ix=1
                endif
  
                if(ix.gt.nside) then
                  ix=nside
                endif

                if(iy.le.0) then
                  iy=1
                endif

                if(iy.gt.nside) then
                  iy=nside
                endif
c               x and y index of the point 
c               on the children's level
c
                do kk=1,4
                  ic=ichildbox(kk,ibox)
                  if((icolbox(ic).eq.ix).and.(irowbox(ic).eq.iy)) then
                    nsinchild(kk)=nsinchild(kk)+1
                    npbox(ic)=npbox(ic)+1
                    isinchild(nsinchild(kk),kk)=jj
                  endif
                enddo
              enddo
c
c             all points in box ibox sorted into isinchild,
c             indicating which child the points belong to. 
c             now overwrite the iss-th to ies-th spots in
c             xpts and ipold together
c
              jpt=iss
              do jc=1,4
                ic=ichildbox(jc,ibox)
                istartbox(ic)=jpt
                do js=1,nsinchild(jc)
                  idnow=isinchild(js,jc)
                  xpts(1,jpt)=xtmp(1,idnow)  
                  xpts(2,jpt)=xtmp(2,idnow)
                  ipold(jpt)=itmp(idnow)
                  jpt=jpt+1
                enddo
              enddo
c             done with box: ibox
              iflag=1
            endif
          endif
        enddo
c
        if((iflag .eq. 0) .and. (i .eq. nlev)) then
c         nothing was divides at the last level,
c         exit the loop (on levels)
          goto 234
        endif
c     the end of the loop on levels
      enddo
c
c----------------------------------
c     2) an upward sweep to coarsen 
c----------------------------------
c
234   do ii=1,nboxes
        iempty(ii)=0
      enddo
c       before merging, there is no non-existant box
c       (with a slight abuse of language, 
c        we call them empty boxes)
c
      do i=nlev,0,-1
c       to merge, go up the tree level by level
        istart=istartlev(i)
        iend=istart+nblevel(i)-1
        do j=istart, iend
          ibox=iboxlev(j)
c         two criteria: 1) no points in it
c                       2) children are all leaves
          nss=npbox(ibox)
          nleafkids=0
          if(ichildbox(1,ibox) .gt. 0) then
            kids(1)=ichildbox(1,ibox)      
            kids(2)=ichildbox(2,ibox)      
            kids(3)=ichildbox(3,ibox)      
            kids(4)=ichildbox(4,ibox)      
c
            do k=1,4
              if(ichildbox(1,kids(k)).lt.0) then
                nleafkids=nleafkids+1
              endif
            enddo
          endif
c         got nss and nleafkids
c
          if((nss.le.0).and.(nleafkids.eq.4)) then
c           a regular box containing no boundary points
c           with four leaf kids, do the following:
            do k2=1,ndeg
            do k1=1,ndeg
              k=(k2-1)*ndeg+k1
              valk1(k1,k2)=fval(k,kids(1))
              valk2(k1,k2)=fval(k,kids(2))
              valk3(k1,k2)=fval(k,kids(3))
              valk4(k1,k2)=fval(k,kids(4))
            enddo
            enddo
c
            call chebyk2p(ndeg, valk1, valk2, valk3, valk4,
     1           ftemp, coefftemp)
c
            call geterror(ndeg,coefftemp,error)
            hh = dble(4**levelbox(ibox))
            epscaled = eps * hh
c
            if(error .lt. epscaled) then
              call coarsen1(ibox,iparentbox,ichildbox,
     1             nboxes,irowbox,icolbox,levelbox,nlev,
     2             istartlev,nblevel,iboxlev,iempty)
c
              do k2=1,ndeg
              do k1=1,ndeg
                k=(k2-1)*ndeg+k1
                fval(k,ibox)=ftemp(k1,k2)
              enddo
              enddo
c             remember to assign interpolated function values
c             to the box that has become a leaf box
            endif            
          endif
        enddo
      enddo
c
c------------------------------------------------------
c     3) clean up the intermediate tree with deleted
c        boxes, reorder fval at the same time
c        update the point-tree relation when isort>0
c------------------------------------------------------
c
      call reorder(iparentbox, ichildbox, nboxes, 
     1     irowbox, icolbox, levelbox, nlev,
     2     istartlev, nblevel, iboxlev, iempty, 
     3     idold2new, nold, itemp, itemp4)
c
      call reorderf(ndeg*ndeg, nold, fval, idold2new)

      if(isort .gt. 0) then
        call treesort(nlev,levelbox,iparentbox,ichildbox,
     1       icolbox,irowbox,nboxes,nblevel,iboxlev,
     2       istartlev,xpts,npts,npbox,ipold,
     3       istartbox,ibofp,cent0,xsize0)
      endif
c
c
      deallocate(ftemp)
      deallocate(coefftemp)
      deallocate(itemparray)
c
      deallocate(isinchild)
      deallocate(itmp)
      deallocate(xtmp)
c
      deallocate(iempty)
      deallocate(itemp)
      deallocate(itemp4)
      deallocate(idold2new)
c
      deallocate(valk1)
      deallocate(valk2)
      deallocate(valk3)
      deallocate(valk4)
c

      end subroutine




c-----------------------------------------------------
c     subroutine movemesh
c-----------------------------------------------------
c
c     Given a set of points, a tree, and a function
c     resolved on the old tree, this subroutine adjusts
c     the tree so that each leaf box contains no more 
c     than a given number of points, and is level-
c     restricted.
c     It then interpolates the function to the new grid
c     and sorts the points onto the new tree.
c
c     INPUT:
c     maxboxes: the maximum number of boxes allowed
c     maxlevel: deepest level allowed
c     maxppl: the maximum number of boundary points allowed
c             in leaf box
c     levelbox - istartlev: the quad-tree structure
c       (RMK: including colleague list!)
c     npts: number of points
c     xpts: coordinates of the points

c     ndeg: degree of chebyshev approx 
c     fval: function values on the ndeg x ndeg grid of
c           the leaf nodes
c     eps: error tolerance of function approx
c     iperiod: iperiod=0 => non-periodic
c              iperiod=1 => periodic
c
c     OUTPUT:
c     the tree overwritten to reflect the changes
c     the function values interpolated to the new grid
c     the points reordered, and point-tree relation 
c         written in the following array:
c
c     npbox: number of sources in each box 
c                     (including non-leaf boxes)
c     ipold: ids of points on input
c            (before the overwriting)
c     istartbox: starting point in the array xpts
c                for each box
c     ibofp: indices of leaf boxes that the point
c            belongs to (after overwriting)
c
c-----------------------------------------------------
c
      subroutine movemesh(maxboxes, maxlevel, maxppl,
     1           levelbox, icolbox, irowbox, icolleagbox,
     2           nboxes, nlev, iparentbox, ichildbox, 
     3           nblevel, iboxlev, istartlev, npts, xpts, 
     4           npbox, ipold, istartbox, ibofp, ndeg, 
     5           fval, eps, iperiod, ilevrestr, ifcolls,
     6           cent0, xsize0)
      implicit real*8 (a-h,o-z)
c-----global variables
      integer maxboxes, maxlevel, maxppl
      integer nboxes, nlev, npts, ndeg
      integer iperiod, ilevrestr
      integer levelbox(1)
      integer icolbox(1), irowbox(1), icolleagbox(9,1)
      integer iparentbox(1), ichildbox(4,1)
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer npbox(1), ipold(1), istartbox(1), ibofp(1)
      real*8 xpts(2,1), fval(ndeg*ndeg,1), eps
      real*8 cent0(2), xsize0
c-----local variables
      integer, allocatable:: itemparray(:)
      integer, allocatable:: iflag(:)
c
      allocate(itemparray(maxboxes))
      allocate(iflag(maxboxes))
c
c     1) sort points into the tree
      call treesort(nlev,levelbox,iparentbox,ichildbox,
     1     icolbox,irowbox,nboxes,nblevel,iboxlev,
     2     istartlev,xpts,npts,npbox,ipold,istartbox,ibofp,
     3     cent0,xsize0)
c
c     2) adjust the tree
      isort=0
      call adaptree(maxboxes, maxlevel, maxppl,
     1     levelbox, icolbox, irowbox, nboxes,
     2     nlev, iparentbox, ichildbox, nblevel,
     3     iboxlev, istartlev, npts, xpts, 
     4     npbox, ipold, istartbox, ibofp,
     5     ndeg, fval, eps, isort, cent0, xsize0)
c
c     3) check level-restriction
      if(ilevrestr .gt. 0) then
        ifixflag=0
        call restriction(levelbox, iparentbox, ichildbox,
     1       icolbox, irowbox, icolleagbox, nboxes, nlev,
     2       nblevel, iboxlev, istartlev, iperiod, ifixflag)
c
c       fix the tree if needed
        if(ifixflag .eq. 1) then
          call fixtreenf(levelbox,iparentbox,ichildbox,
     1         icolbox,irowbox,icolleagbox,nboxes,nlev,
     2         nblevel,iboxlev,istartlev,iperiod,iflag,
     3         maxboxes,itemparray,ndeg,fval)
        endif     
      endif
c
c     4) sort the points into the new tree
c     reverse the order of the point before sorting
      call rvecnew2old(2,npts,xpts,ipold)
      call treesort(nlev,levelbox,iparentbox,ichildbox,
     1     icolbox,irowbox,nboxes,nblevel,iboxlev,
     2     istartlev,xpts,npts,npbox,ipold,istartbox,ibofp,
     3     cent0,xsize0)
c
c     5) generate colleague list
      if(ifcolls .gt. 0) then
        call mkcolls(icolbox, irowbox, icolleagbox, nboxes,
     1       nlev, iparentbox, ichildbox, nblevel, iboxlev,
     2       istartlev, iperiod)
      endif

      deallocate(itemparray)
      deallocate(iflag) 

      end subroutine




c********************************************************
c     subroutine setf
c********************************************************
c
c     the following subroutine defines the array that represents the 
c     right hand side of the poisson equation. 
c
c     input:
c
c     levelbox is an array determining the level of each box
c     nboxes is the total number of boxes
c     icolbox denotes the column of each box
c     irowbox denotes the row of each box
c     ichildbox denotes the four children of each box
c     nlev is the finest level
c     iboxlev is the array in which the boxes are arranged
c     nblevel is the total number of boxes per level
c     istartlev is the pointer to where each level begins in the
c               iboxlev array
c     h is the real function that is the right hand side
c       of the poisson equation
c
c     output:
c
c     fright is the right hand side defined on the tree
c
c********************************************************
c
      subroutine setf(ndeg,fright,
     1       icolbox, irowbox, ichildbox,nlev,
     2       nblevel, iboxlev, istartlev, h,
     3       cent0, xsize0)
      implicit none
c-----global variables
      integer *4  ndeg
      integer *4  icolbox(1), irowbox(1), nlev
      integer *4  ichildbox(4,1)
      integer *4  nblevel(0:1), iboxlev(1), istartlev(0:1)
      real *8  fright(ndeg*ndeg,1), xf(ndeg*ndeg), yf(ndeg*ndeg)
      real*8 cent0(2), xsize0, x0, y0
c-----local variables
      integer *4  i, ibox, j, k, l
      real *8  temp1, pi2n
      real *8  xstart
      real *8  xshift, xx(ndeg)
      real *8  xscale(ndeg)
      real *8  yshift
c-----external functions
      real *8  h

c     xx: Chebyshev nodes on [-0.5,0.5]
      pi2n=dacos(-1.0d0)/(2.0d0*ndeg)
      do i=1,ndeg
         xx(i)=dcos((2*(ndeg-i)+1)*pi2n)/2.0d0
      enddo 

c      temp1 = 1.0d0
      temp1 = xsize0
      x0 = cent0(1)-xsize0/2.0d0
      y0 = cent0(2)-xsize0/2.0d0
      do k = 0, nlev
c      xstart = (1.0d0 - temp1) / 2.0d0

c     Chebyshev nodes scaled and shifted to the first box 
      do i=1,ndeg
c        xscale(i)=xx(i)*temp1-xstart
        xscale(i)=xx(i)*temp1
      enddo

      do 100 i = istartlev(k), istartlev(k) + nblevel(k) - 1
        ibox = iboxlev(i)
        if(ichildbox(1,ibox) .gt. 0)goto 100

        xshift  =  x0+dble(icolbox(ibox)-0.5d0) * temp1
        yshift  =  y0+dble(irowbox(ibox)-0.5d0) * temp1

        do j = 1, ndeg
          do l = 1, ndeg
            xf(ndeg*(l-1)+j) = xscale(j) + xshift
            yf(ndeg*(j-1)+l) = xscale(j) + yshift
          end do
        end do
       
        do j = 1, ndeg*ndeg
          fright(j,ibox) = h(xf(j),yf(j))
        end do
100   continue
      temp1 = temp1 / 2.0d0
      end do

      return
      end subroutine





c*************************************************************
c     subroutine sortsrc
c*************************************************************
c
c     Given the tree structure and the chunk discretization of
c     the boundary, sort the boundary points into leaf boxes
c     of the tree
c
c     INPUT:
c     nlev-istartlev: the tree structure, see subroutine mktree
c                     for details
c     src: points on the boundary
c     nsrc: number of sources
c     
c     OUTPUT:
c     nsinbox: number of sources in each box 
c                     (including non-leaf boxes)
c     issorted: indices of the sources in each leaf box, sorted
c              in the order of boxes
c              (at least, include info of sources in leaf boxes)
c     istartbox: starting point in issorted for each box
c     ibofs: indices of leaf boxes that the source point belongs to
c
c
c*************************************************************
c
      subroutine sortsrc(nlev,levelbox,iparentbox,ichildbox,
     1           icolbox,irowbox,nboxes,nblevel,iboxlev,istartlev,
     2           src,nsrc,nsinbox,issorted,istartbox,ibofs,
     3           cent0,xsize0)
      implicit real*8 (a-h,o-z)
c     global vars
      integer nlev, nboxes, nsrc
      integer levelbox(nboxes), iparentbox(nboxes)
      integer ichildbox(4,nboxes)
      integer icolbox(nboxes), irowbox(nboxes)
      integer nblevel(0:nlev), iboxlev(nboxes)
      integer istartlev(0:nlev)
      integer nsinbox(nboxes), issorted(nsrc)
      integer istartbox(nboxes), ibofs(nsrc)
      real*8 src(2,nsrc), cent0(2), xsize0
c     local vars
      integer nsinchild(4), isinchild(nsrc,4)
c
      do i=1,nboxes
        nsinbox(i)=0
        istartbox(i)=-1
      enddo 
      
      do j=1,nsrc
        issorted(j)=j
      enddo
c
      iroot=iboxlev(istartlev(0))
c
      nsinbox(iroot)=nsrc
      istartbox(iroot)=1
c
c         the root box contains all the sources
c
c----------------------------------
c
      ix=1
      iy=1
      xlength=xsize0
      x0=cent0(1)-xsize0/2.0d0
      y0=cent0(2)-xsize0/2.0d0

      nside=1

      do l=0,nlev-1
c       go through all the levels
        xlength=xlength/2.0d0
        nside=nside*2

        istart=istartlev(l)
        iend=istartlev(l)+nblevel(l)-1

        do ii=istart,iend
          ibox=iboxlev(ii)
c         go through all non-empty boxes that have children
          if((nsinbox(ibox) .gt. 0) .and. 
     1               (ichildbox(1,ibox)).gt.0) then
c-------------------
            iss=istartbox(ibox)
            ies=istartbox(ibox)+nsinbox(ibox)-1
c
            do kk=1,4
              nsinchild(kk)=0
c             initialize the number of points in each
c             child to be zero
            enddo
c
c             go through all the sources points in the box
c             (which occupy adjacent spots in issorted)      
c
            do j=iss,ies
              isrc=issorted(j)
              xpt=src(1,isrc)
              ypt=src(2,isrc)

c              ix=ceiling((xpt+0.5d0)/xlength)
c              iy=ceiling((ypt+0.5d0)/xlength)
              ix=ceiling((xx-x0)/xlength)
              iy=ceiling((yy-y0)/xlength)

              if(ix.le.0) then
                ix=1
              endif

              if(ix.gt.nside) then
                ix=nside
              endif

              if(iy.le.0) then
                iy=1
              endif

              if(iy.gt.nside) then
                iy=nside
              endif
c               x and y index of the point 
c               on the children's level
c
              do kk=1,4
                ic=ichildbox(kk,ibox)
                if((icolbox(ic).eq.ix).and.(irowbox(ic).eq.iy)) then
                  nsinchild(kk)=nsinchild(kk)+1
                  nsinbox(ic)=nsinbox(ic)+1
                  isinchild(nsinchild(kk),kk)=isrc
                endif
              enddo
            enddo
c
c           all points in box ibox sorted into isinchild,
c           indicating which child the points belong to. 
c           now overwrite the iss-th to ies-th spot in issorted
c
            jpt=iss
            do jc=1,4
              ic=ichildbox(jc,ibox)
              istartbox(ic)=jpt
              do js=1,nsinchild(jc)
                issorted(jpt)=isinchild(js,jc)
                jpt=jpt+1
              enddo
            enddo
c-------------------
          endif
        enddo
c       go to the next level
      enddo

c     now get information for ibofs
      do l=0,nlev
        istart=istartlev(l)
        iend=istartlev(l)+nblevel(l)-1
        do ii=istart,iend
          ibox=iboxlev(ii)
          if((ichildbox(1,ibox).lt.0).and.(nsinbox(ibox).gt.0)) then
            iss=istartbox(ibox)
            ies=iss+nsinbox(ibox)-1
            do kk=iss,ies
              isrc=issorted(kk)
              ibofs(isrc)=ibox
            enddo
          endif
        enddo
      enddo

      end subroutine





c-----------------------------------------------------
c     subroutine adaptreef
c-----------------------------------------------------
c
c     Given a tree, and a function resolved on the 
c     leaf nodes, and another function (by function
c     evaluation routine), refine and coarsen the tree
c     so that both functions are resolved but not
c     over-resolved on the new tree
c
c     INPUT:
c     maxboxes: the maximum number of boxes allowed
c     maxlevel: deepest level allowed
c     levelbox - istartlev: the quad-tree structure
c     ndeg: degree of chebyshev approx on each leaf node
c     uval: array of function values (resolved already)
c     feval: function evaluation routine of another
c            function
c     eps: error tolerance
c     t: parameter (e.g. time) needed in feval
c
c     OUTPUT:
c     the tree: overwritten to reflect the change
c     uval: interpolated to the new tree
c
c-----------------------------------------------------
c
      subroutine adaptreef(maxboxes, maxlevel, levelbox,
     1           icolbox, irowbox, nboxes, nlev, iparentbox,
     2           ichildbox, nblevel, iboxlev, istartlev,
     3           ndeg, uval, feval, t, eps, cent0, xsize0)
      implicit real*8 (a-h,o-z)
c-----global variables
      integer maxboxes, maxlevel
      integer nboxes, nlev, ndeg
      integer levelbox(maxboxes)
      integer icolbox(maxboxes), irowbox(maxboxes)
      integer iparentbox(maxboxes), ichildbox(4,maxboxes)
      integer nblevel(0:maxlevel), iboxlev(maxboxes)
      integer istartlev(0:maxlevel)
      real*8 uval(ndeg*ndeg,maxboxes), eps
      real*8 cent0(2), xsize0
      external feval
c-----local variables
      integer kids(4)
      integer, allocatable:: iempty(:)
      integer, allocatable:: itemp(:), itemp4(:,:)
      integer, allocatable:: idold2new(:)
      integer, allocatable:: itemparray(:)
      real*8 wsave(1000), chwork(1000)
      real*8 xf(ndeg), yf(ndeg), ftemp(ndeg,ndeg)
      real*8 coefftemp(ndeg,ndeg)
      real*8, allocatable:: valk1(:,:), valk2(:,:)
      real*8, allocatable:: valk3(:,:), valk4(:,:)
      allocate(itemparray(maxboxes))
      allocate(iempty(maxboxes))
      allocate(itemp(maxboxes))
      allocate(itemp4(4,maxboxes))
      allocate(idold2new(maxboxes))
c
      allocate(valk1(ndeg,ndeg))
      allocate(valk2(ndeg,ndeg))
      allocate(valk3(ndeg,ndeg))
      allocate(valk4(ndeg,ndeg))
c
c      write(*,*) ndeg, t, eps
c
c---------------------------------------------------------
c     1) refinement sweep. go down the tree level by level
c        refine leaf boxes if feval(:,t) is unresolved
c        no need to worry about uval, since by assumption
c        it is resolved
c---------------------------------------------------------
c
      ix=1
      iy=1
c      xlength=1.0d0
      xlength=xsize0
      x0=cent0(1)-xsize0/2.0d0
      y0=cent0(2)-xsize0/2.0d0
      nside=1
c
      call chxcin(ndeg,wsave)
      do i=0, maxlevel-1
        iflag=0
c          iflag: if the last level has been divided at all
        xlength=xlength/2.0d0
        nside=nside*2
c
        istart=istartlev(i)
        iend=istart+nblevel(i)-1
        do j=istart, iend
          ibox=iboxlev(j)
          if(ichildbox(1,ibox) .lt. 0) then
            call mkgrid(xf,yf,icolbox(ibox),irowbox(ibox),
     1           levelbox(ibox),ndeg,cent0,xsize0)
c
            do jj = 1, ndeg
              do ii = 1, ndeg
                call feval(xf(jj), yf(ii), t, ftemp(jj,ii))
              enddo
            enddo
c
            call getcoeff(ndeg,ndeg,ftemp,coefftemp,wsave)
            call geterror(ndeg,coefftemp,error)
            hh = dble(4**levelbox(ibox))
            epscaled = eps * hh


            if(error .ge. epscaled)then
c             call subdivide
c              write(*,*) ibox, error, epscaled,hh
              call subdivide1(ibox,iparentbox,ichildbox,
     1             nboxes,irowbox,icolbox,levelbox,nlev,
     2             istartlev, nblevel, iboxlev,itemparray)
c             need to assign function values to 
              ik1=ichildbox(1,ibox)
              ik2=ichildbox(2,ibox)
              ik3=ichildbox(3,ibox)
              ik4=ichildbox(4,ibox)
              call chebyp2k(ndeg,uval(1,ibox),uval(1,ik1),
     1             uval(1,ik2),uval(1,ik3),uval(1,ik4))
c
              iflag = 1
            endif
          endif
        enddo
c
        if((iflag .eq. 0) .and. (i .eq. nlev))then
          goto 110
        endif
      enddo

c------------------------------------------
c 2. go up the tree to coarsen:
c    if on a certain parent box, both functions
c    are resolved, merge its children
c------------------------------------------
c
110   do ii=1,nboxes
        iempty(ii)=0
      enddo   

      do i=nlev,0,-1
c       to merge, go up the tree level by level
        istart=istartlev(i)
        iend=istart+nblevel(i)-1
        do j=istart, iend
          ibox=iboxlev(j)
          nleafkids=0
          if(ichildbox(1,ibox) .gt. 0) then
            kids(1)=ichildbox(1,ibox)      
            kids(2)=ichildbox(2,ibox)      
            kids(3)=ichildbox(3,ibox)      
            kids(4)=ichildbox(4,ibox)      
c
            do k=1,4
              if(ichildbox(1,kids(k)).lt.0) then
                nleafkids=nleafkids+1
              endif
            enddo
          endif
c         got nleafkids
            
          if(nleafkids .eq. 4) then
            call mkgrid(xf,yf,icolbox(ibox),irowbox(ibox),
     1           levelbox(ibox),ndeg,cent0,xsize0)
c
            do jj = 1, ndeg
              do ii = 1, ndeg
                call feval(xf(jj), yf(ii), t, ftemp(jj,ii))
              enddo
            enddo
c
            call getcoeff(ndeg,ndeg,ftemp,coefftemp,wsave)
            call geterror(ndeg,coefftemp,ferror)
c            ferror=1.0d0
c                got ferror

            do k2=1,ndeg
            do k1=1,ndeg
              k=(k2-1)*ndeg+k1
              valk1(k1,k2)=uval(k,kids(1))
              valk2(k1,k2)=uval(k,kids(2))
              valk3(k1,k2)=uval(k,kids(3))
              valk4(k1,k2)=uval(k,kids(4))
            enddo
            enddo
c
            call chebyk2p(ndeg, valk1, valk2, valk3, valk4,
     1           ftemp, coefftemp)
            call geterror(ndeg,coefftemp,uerror)
c                got uerror
c
            hh = dble(4**levelbox(ibox))
            epscaled = eps * hh 
c                the error tolerance scaled

c            write(*,*) ibox, ferror, uerror, epscaled
c
            if((uerror.lt.epscaled) .and. 
     1         (ferror.lt.epscaled)) then
c             the condition is correct here!
              call coarsen1(ibox,iparentbox,ichildbox,
     1             nboxes,irowbox,icolbox,levelbox,nlev,
     2             istartlev,nblevel,iboxlev,iempty)
              do k2=1,ndeg
              do k1=1,ndeg
                k=(k2-1)*ndeg+k1
                uval(k,ibox)=ftemp(k1,k2)
              enddo
              enddo
c             assign interpolated function values
c             to the box that has become a leaf box
            endif
          endif
        enddo
      enddo

c------------------------------------------
c     3) clean up the intermediate tree with deleted
c        boxes, reorder uval at the same time
c------------------------------------------
c
      call reorder(iparentbox, ichildbox, nboxes, 
     1     irowbox, icolbox, levelbox, nlev,
     2     istartlev, nblevel, iboxlev, iempty, 
     3     idold2new, nold, itemp, itemp4)
c

      call reorderf(ndeg*ndeg, nold, uval, idold2new)


      deallocate(itemparray)
      deallocate(iempty)
      deallocate(itemp)
      deallocate(itemp4)
      deallocate(idold2new)
c
      deallocate(valk1)
      deallocate(valk2)
      deallocate(valk3)
      deallocate(valk4)
c
c      write(*,*) nboxes
c
      end subroutine




      
c-----------------------------------------------------
c     suboroutine outputpot
c-----------------------------------------------------
c
c     this subroutine writes function values of a given
c     function sampled on the leaf nodes of a given
c     quadtree
c
c-----------------------------------------------------
c
      subroutine outputpot(levelbox, icolbox, irowbox, 
     1           nboxes, nlev, iparentbox, ichildbox,
     2           nblevel, iboxlev, istartlev, ndeg, pot,
     3           iprint, cent0, xsize0)
      implicit real*8 (a-h,o-z)
      integer nlev, nboxes, ndeg
      integer levelbox(1)
      integer icolbox(1), irowbox(1)
      integer iparentbox(1), ichildbox(4,1)
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
      real*8 pot(ndeg*ndeg,nboxes)
      real*8, allocatable:: xf(:,:), yf(:,:)
      real*8, allocatable:: xx(:), yy(:)
      real*8 cent0(2), xsize0
c
      allocate(xf(ndeg**2,nboxes))
      allocate(yf(ndeg**2,nboxes))

      allocate(xx(ndeg))
      allocate(yy(ndeg))
c
c      call getxyclassical8(xf, yf,
c     1     icolbox, irowbox, ichildbox,nlev,
c     2     nblevel, iboxlev, istartlev,
c     3     cent0, xsize0)
      do l=0, nlev
        istart=istartlev(l)
        iend=istart+nblevel(l)-1
        do ii=istart, iend
          ibox=iboxlev(ii)
          if(ichildbox(1,ibox) .le. 0) then
            call mkgrid(xx,yy,icolbox(ibox),
     1           irowbox(ibox),levelbox(ibox),ndeg,
     2           cent0,xsize0)
            do k2=1, ndeg
            do k1=1, ndeg
              xf(ndeg*(k2-1)+k1,ibox)=xx(k1)
              yf(ndeg*(k2-1)+k1,ibox)=yy(k2)
            enddo
            enddo
          endif
        enddo 
      enddo
c     change to arbitary order, delete the 
c     obsolete routine getxyclassical8
c

      do i = 0, nlev
      do 100 j = istartlev(i), istartlev(i) + nblevel(i) - 1
       ibox = iboxlev(j)
       if(ichildbox(1,ibox) .gt. 0)goto 100
       do l = 1, ndeg
       write(iprint,*)' y(:,',l,') = ['
        do jcntr = 1, ndeg
         write(iprint,*)yf(ndeg*(jcntr-1)+l,ibox)
        end do
       write(iprint,*)' ];'
       write(iprint,*)' x(:,',l,') = ['
        do jcntr = 1, ndeg
         write(iprint,*)xf(ndeg*(jcntr-1)+l,ibox)
        end do
       write(iprint,*)' ];'
       write(iprint,*)' pot(:,',l,') = ['
        do jcntr = 1, ndeg
         write(iprint,*)pot(ndeg*(jcntr-1)+l,ibox)
        end do
       write(iprint,*)' ];' 
      end do
       write(iprint,*)'surf(x,y,pot);'
c       write(iprint,*)'view(2);'
       write(iprint,*)'colormap(jet);'
       write(iprint,*)'hold on;'
       write(iprint,*)'axis equal;'
100   continue
      end do

      write(iprint,*) 'view(2);'

      deallocate(xf)
      deallocate(yf)

      deallocate(xx)
      deallocate(yy)
c
      end subroutine





c-----------------------------------------------------
c     subroutine interppot
c-----------------------------------------------------
c
c     This subroutine interpolates the function sampled on 
c     the Chebyshev grids of leaf nodes of an adaptive
c     quadtree to a regular tensor product grid 
c
c     INPUT: 
c     levelbox - istartlev: the tree structure
c     ndeg: degree of function approx on leaf nodes
c     pot: the array of function values
c     nx, ny: number of grid points in x and y direction
c     xg, yg: the mesh grid
c     fg: the function values on the mesh grid
c
c-----------------------------------------------------
c
      subroutine interppot(levelbox, icolbox, irowbox, 
     1           nboxes, nlev, iparentbox, ichildbox,
     2           nblevel, iboxlev, istartlev, ndeg, pot,
     3           ng, xg, yg, fg, cent0, xsize0)
      implicit real*8 (a-h,o-z)
      integer nlev, nboxes, ndeg, nx
      integer levelbox(1)
      integer icolbox(1), irowbox(1)
      integer iparentbox(1), ichildbox(4,1)
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
      real*8 pot(ndeg*ndeg,nboxes)
      real*8 xg(ng), yg(ng), fg(ng)
      real*8 cent0(2), xsize0
c     local vars
      integer npts
      integer,allocatable:: npbox(:), ipold(:)
      integer,allocatable:: istartbox(:), ibofp(:)
      real*8 wsave(1000), chwork(1000)
      real*8, allocatable:: xpts(:,:), chpot(:,:)
c
      npts=ng
      allocate(xpts(2,npts))
      allocate(chpot(ndeg,ndeg))
      allocate(npbox(nboxes))
      allocate(ipold(npts))
      allocate(istartbox(nboxes))
      allocate(ibofp(npts))
c
      do i=1,npts
        xpts(1,i)=xg(i)
        xpts(2,i)=yg(i)
      enddo
c
c     sort the grid points into the tree
      t1=second()
      call treesort(nlev,levelbox,iparentbox,ichildbox,
     1    icolbox,irowbox,nboxes,nblevel,iboxlev,istartlev,
     2    xpts,npts,npbox,ipold,istartbox,ibofp,
     3    cent0,xsize0)
      t2=second()
c      write(*,*) 't_treesort=',t2-t1
c
c     go through leaf nodes, 
c     if nonempty, do chebyshev interpolation
c
      t1=second()
      call chxcin(ndeg, wsave)
      do ib=1, nboxes
        npb=npbox(ib)
c        write(*,*) ib, npb
        if((ichildbox(1,ib).lt.0) .and. (npb.gt.0)) then
          call chexfc2d(pot(1,ib),ndeg,chpot,wsave,chwork)
c
          lev=levelbox(ib)
          icol=icolbox(ib)
          irow=irowbox(ib)
  
          xlength=0.5d0**lev*xsize0
          x0=cent0(1)-xsize0/2.0d0
          y0=cent0(2)-xsize0/2.0d0
c
          a=dble(icol-1)*xlength+x0
          b=a+xlength
          c=dble(irow-1)*xlength+y0
          d=c+xlength
c           basic parameters (of current box)
c           RMK: for a generic unit box, treesort needs to 
c                change, and these lines need to change
         
          istart=istartbox(ib) 
          iend=istart+npb-1
          do ip=istart, iend
            xx=xpts(1,ip) 
            yy=xpts(2,ip)
            call cheval2d(xx,yy,val,chpot,ndeg,a,b,c,d)
            iold=ipold(ip)
            fg(iold)=val
          enddo
        endif
      enddo
      t2=second()
c      write(*,*) 't_interp=',t2-t1

      deallocate(xpts)
      deallocate(npbox)
      deallocate(ipold)
      deallocate(istartbox)
      deallocate(ibofp)
      deallocate(chpot)

      end subroutine







c-----------------------------------------------------
c     subroutine rvecnew2old
c-----------------------------------------------------
c
c Reorder a vector from the new order to the old order
c given permutation ipold that maps new indices 
c to the old ones
c
c INPUT:
c ndim: dimension of the vector
c nvec: length of the vector
c rvec: the real vector  
c ipold: permutation vector that maps the new indices
c        to the old ones
c
c OUTPUT:
c rvec: rewritten in the new order
c
c-----------------------------------------------------
c
      subroutine rvecnew2old(ndim, nvec, rvec, ipold)
      implicit real*8 (a-h,o-z)
      integer ndim, nvec
      integer ipold(nvec)
      real*8 rvec(ndim,nvec)    
      real*8, allocatable:: rnew(:,:)
c
      allocate(rnew(ndim,nvec))
c
      do i=1,nvec
      do j=1,ndim
        rnew(j,i)=rvec(j,i)
      enddo
      enddo
c
      do i=1,nvec
        io=ipold(i)
        do j=1,ndim
          rvec(j,io)=rnew(j,i)
        enddo
      enddo
c
      deallocate(rnew)
      
      end subroutine






c-----------------------------------------------------
c     subroutine rvecold2new
c-----------------------------------------------------
c
c Reorder a vector from the old order to the new order
c given permutation ipold that maps new indices 
c to the old ones
c
c INPUT:
c ndim: dimension of the vector
c nvec: length of the vector
c rvec: the real vector  
c ipold: permutation vector that maps the new indices
c        to the old ones
c
c OUTPUT:
c rvec: rewritten in the old order
c
c-----------------------------------------------------
c
      subroutine rvecold2new(ndim, nvec, rvec, ipold)
      implicit real*8 (a-h,o-z)
      integer ndim, nvec
      integer ipold(nvec)
      real*8 rvec(ndim,nvec)    
      real*8, allocatable:: rold(:,:)
c
      allocate(rold(ndim,nvec))
c
      do i=1,nvec
      do j=1,ndim
        rold(j,i)=rvec(j,i)
      enddo
      enddo
c
      do i=1,nvec
        io=ipold(i)
        do j=1,ndim
          rvec(j,i)=rold(j,io)
        enddo
      enddo
c
      deallocate(rold)
      
      end subroutine




c*****************************************************
c     subroutine printtree
c*****************************************************
      subroutine printtree(levelbox,icolbox,irowbox,nboxes,
     1           nlev,iparentbox,ichildbox,nblevel,iboxlev,
     2           istartlev,iprint,nleaves,cent0,xsize0)
      implicit real*8 (a-h,o-z)
      integer levelbox(1)
      integer nlev, nboxes
      integer icolbox(1), irowbox(1)
      integer iparentbox(1), ichildbox(4,1)
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
      real*8 cent0(2), xsize0 
c
      nleaves=0
      do l=0,nlev
        istart=istartlev(l)
        iend=istart+nblevel(l)-1
        do i=istart,iend
          ibox=iboxlev(i)
c
          if(ichildbox(1,ibox) .eq. 0) then
            write(*,*) 'ichildbox(1,ibox)=0!'
            write(*,*) 'ibox=',ibox
          endif
c
          if(ichildbox(1,ibox).lt.0) then
            nleaves=nleaves+1
            if(iprint .gt. 0) then
c             if a leaf node, print out end points
              call posbox(xll,yll,dx,dy,icolbox(ibox),
     1             irowbox(ibox),l,cent0,xsize0)
              write(iprint,*) xll,yll,dx,dy
            endif
          endif
        enddo
      enddo
c
c      write(*,*) 'tree data saved in file',iprint
      close(iprint)

      end subroutine






c**********************************************************
c     subroutine fixtreenf2
c**********************************************************
c
c     the following subroutine is designed to take a correctly defined
c     tree and alter it so that no two boxes that contact each other
c     are more than one level apart.  this is corrected only by adding
c     boxes.  the procedure involves flagging down bigger boxes and
c     dividing them and their children as necessary.
c     this routine also produces an array of colleagues for each box
c     that is numbered in the correct order.  all of the children are set
c     so that they satisfy our ordering convention.
c     the algorithm in the periodic case works the same way, it is just 
c     that upon subdivision the new colleagues must be put down to 
c     account for the periodicity.
c
c     input:
c
c     levelbox is an array determining the level of each box
c
c     iparentbox denotes the parent of each box
c
c     ichildbox denotes the four children of each box
c
c     icolbox denotes the column of each box
c
c     irowbox denotes the row of each box
c
c     icolleagbox denotes the colleagues of a given box
c
c     nboxes is the total number of boxes
c
c     nlev is the finest level
c
c     nblevel is the total number of boxes per level
c
c     iboxlev is the array in which the boxes are arranged
c
c     istartlev is the pointer to where each level begins in the
c               iboxlev array
c
c     iperiod denotes what kind of colleagues are to be generated
c             iperiod = 0 : free space
c             iperiod = 1 or 2 : periodic
c             iperiod = 3 : periodic up/down and free space left/right
c
c     iflag is just a dummy array
c
c     maxboxes is the maximum number of boxes we have storage for
c
c     output:
c
c     icolbox, irowbox, icolleagbox, nboxes, and all other
c     arrays describing the tree may be change on output
c
c     fval1 and fval2 may be overwritten too 
c     (interpolated to the new grid)
c
c**********************************************************
      subroutine fixtreenf2(levelbox,iparentbox,ichildbox,
     1           icolbox,irowbox,icolleagbox,nboxes,nlev,
     2           nblevel,iboxlev,istartlev,iperiod,
     3           iflag, maxboxes,itemparray,ndeg,fval1,fval2)
      implicit none
c-----global variables
      integer *4 levelbox(1), icolleagbox(9,1)
      integer *4 maxboxes, ndeg
      integer *4 iparentbox(1), ichildbox(4,1)
      integer *4 icolbox(1), irowbox(1)
      integer *4 nboxes, nlev
      integer *4 nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer *4 iperiod
      integer *4 iflag(maxboxes)
      integer *4 itemparray(1)
      real *8 fval1(ndeg*ndeg, 1)
      real *8 fval2(ndeg*ndeg, 1)
c-----local variables
      integer *4 ichild(4),icoll(9), ibox
      integer *4 i, ipar, itest, j, nb
      integer *4 itemp, ntemp, jcntr, icntr
      integer *4 start, istop
      integer *4 ik1, ik2, ik3, ik4


c     let's sort all of the boxes by level.
c     this takes the place of the ladder structure
c     in the uniform case.  all boxes are sorted
c     into the array and placed in their proper places.
      call sortboxes(levelbox,nboxes,nlev,
     1           nblevel,iboxlev,istartlev)


c     first let's call a subroutine that will
c     generate all of the colleagues for each
c     box.  the colleagues are generated in the
c     correct order so there is no need to 'shuffle'
c     them later on.
      call mkcolls(icolbox,
     1       irowbox,icolleagbox,nboxes,nlev,
     2       iparentbox,ichildbox,nblevel,
     3       iboxlev, istartlev,iperiod)


c     let's initialize all of the flags to zero.
      do i = 1, nboxes
        iflag(i) = 0
      end do

c     find all of the boxes that need to be
c     flagged.  a flagged box will be denoted by 
c     setting iflag(box) = 1.
c     this refers to any box that is directly touching 
c     a box that is more than one level smaller than
c     it.  it is found by performing an upward pass
c     and looking a box's parents parents and seeing
c     if they are childless and contact the given box.
c     note that we only need to get up to level two, as
c     we will not find a violation at a coarser level
c     than that.
      do i = nlev, 2, -1
        do j = istartlev(i), istartlev(i) + nblevel(i) - 1
          ibox = iboxlev(j)
          ipar  = iparentbox(ibox)
          itest = iparentbox(ipar)

            icoll(1) = icolleagbox(1,itest)
            icoll(2) = icolleagbox(2,itest)
            icoll(3) = icolleagbox(3,itest)
            icoll(4) = icolleagbox(4,itest)
            icoll(5) = icolleagbox(5,itest)
            icoll(6) = icolleagbox(6,itest)
            icoll(7) = icolleagbox(7,itest)
            icoll(8) = icolleagbox(8,itest)
            icoll(9) = icolleagbox(9,itest)

            ichild(1) = ichildbox(1,itest)
            ichild(2) = ichildbox(2,itest)
            ichild(3) = ichildbox(3,itest)
            ichild(4) = ichildbox(4,itest)


          do nb = 1, 9
            itemp = icoll(nb)
            if(ichildbox(1,itemp) .lt. 0)then
c             the neighboring box is not divided
c             we could have problems.
              if (nb .eq. 1)then
                if(ipar .eq. ichild(4))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 2)then
                if(ipar .eq. ichild(3) .or. ipar .eq. ichild(4))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 3)then
                if(ipar .eq. ichild(3))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 4)then
                if(ipar .eq. ichild(4) .or. ipar .eq. ichild(1))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 6)then
                if(ipar .eq. ichild(2) .or. ipar .eq. ichild(3))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 7)then 
                if(ipar .eq. ichild(1))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 8)then
                if(ipar .eq. ichild(1) .or. ipar .eq. ichild(2))then
                    iflag(itemp) = 1
                end if
              elseif (nb .eq. 9)then
                if(ipar .eq. ichild(2))then
                    iflag(itemp) = 1
                end if
              endif
            endif
          end do
        end do
      end do


c     find all of the boxes that need to be
c     given a flag+.  a flag+ box will be denoted by 
c     setting iflag(box) = 2.
c     this refers to any box that is not already flagged
c     and is bigger than and is contacting a flagged box
c     or another box that has already been given a flag+.
c     it is found by performing an upward pass
c     and looking at a flagged box's parents colleagues
c     and a flag+ box's parents colleagues and seeing if
c     they are childless and present the case where a 
c     bigger box is contacting a flagged or a flag+ box.
      do i = nlev, 2, -1
        do j = istartlev(i), istartlev(i) + nblevel(i) - 1
         ibox = iboxlev(j)
          if(iflag(ibox) .eq. 1 .or. iflag(ibox) .eq. 2)then

          ipar  = iparentbox(ibox)
 
            icoll(1) = icolleagbox(1,ipar)
            icoll(2) = icolleagbox(2,ipar)
            icoll(3) = icolleagbox(3,ipar)
            icoll(4) = icolleagbox(4,ipar)
            icoll(5) = icolleagbox(5,ipar)
            icoll(6) = icolleagbox(6,ipar)
            icoll(7) = icolleagbox(7,ipar)
            icoll(8) = icolleagbox(8,ipar)
            icoll(9) = icolleagbox(9,ipar)

            ichild(1) = ichildbox(1,ipar)
            ichild(2) = ichildbox(2,ipar)
            ichild(3) = ichildbox(3,ipar)
            ichild(4) = ichildbox(4,ipar)
          

          do nb = 1, 9
            itemp = icoll(nb)
c           let's check using the same criteria as above, but noting that
c           a flag will take precedence over a flag+.
            if(ichildbox(1,itemp) .lt. 0 
     1          .and. iflag(itemp) .ne. 1)then
c             the neighboring box is not divided
c             we could have problems.
              if (nb .eq. 1)then
                if(ibox .eq. ichild(4))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 2)then
                if(ibox .eq. ichild(3) .or. ibox .eq. ichild(4))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 3)then
                if(ibox .eq. ichild(3))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 4)then
                if(ibox .eq. ichild(4) .or. ibox .eq. ichild(1))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 6)then
                if(ibox .eq. ichild(2) .or. ibox .eq. ichild(3))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 7)then
                if(ibox .eq. ichild(1))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 8)then
                if(ibox .eq. ichild(1) .or. ibox .eq. ichild(2))then
                    iflag(itemp) = 2
                end if
              elseif (nb .eq. 9)then
                if(ibox .eq. ichild(2))then
                    iflag(itemp) = 2
                end if
              endif
            endif
           end do
          endif
        end do
      end do


c     now let's divide the boxes that need to be immediately
c     divided up.  all of the flagged and flag+ boxes need to
c     be divided one time.  the distinction lies in the fact
c     that the children of a flag+ box will never need to be
c     divided but the children of a flagged box may need to 
c     be divided further.
c     below, all flagged and flag+ boxes are divided once.  the
c     children of a flag+ box are left unflagged while those of
c     the flagged boxes are given a flag++ (denoted by setting
c     iflag(box) = 3) which will be needed in the downward pass.     
      ntemp = nboxes
      do i = 1, ntemp
c      divide flagged boxes:
       if (iflag(i) .eq. 1)then

         if(ichildbox(1,i) .lt. 0)then
         call subdivide(i,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c        after calling subdivide or subdivide1, assign fval
         ik1=ichildbox(1,i)
         ik2=ichildbox(2,i)
         ik3=ichildbox(3,i)
         ik4=ichildbox(4,i)
         call chebyp2k(ndeg,fval1(1,i),fval1(1,ik1),fval1(1,ik2),
     1                 fval1(1,ik3),fval1(1,ik4))
         call chebyp2k(ndeg,fval2(1,i),fval2(1,ik1),fval2(1,ik2),
     1                 fval2(1,ik3),fval2(1,ik4))
c------------------------
         endif


c        give flag++ to children of flagged boxes.
         itemp = ichildbox(1,i)
         iflag(itemp) = 3

         itemp = ichildbox(2,i)
         iflag(itemp) = 3

         itemp = ichildbox(3,i)
         iflag(itemp) = 3

         itemp = ichildbox(4,i)
         iflag(itemp) = 3

c      divide flag+ boxes.
       elseif (iflag(i) .eq. 2)then
  
         if(ichildbox(1,i) .lt. 0)then
         call subdivide(i,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c        after calling subdivide or subdivide1, assign fval
         ik1=ichildbox(1,i)
         ik2=ichildbox(2,i)
         ik3=ichildbox(3,i)
         ik4=ichildbox(4,i)
         call chebyp2k(ndeg,fval1(1,i),fval1(1,ik1),fval1(1,ik2),
     1                 fval1(1,ik3),fval1(1,ik4))
         call chebyp2k(ndeg,fval2(1,i),fval2(1,ik1),fval2(1,ik2),
     1                 fval2(1,ik3),fval2(1,ik4))
c------------------------
         endif

       endif
      end do 


c     now we need to do a downward pass.
c     we will concern ourselves only with the children of
c     flagged boxes and their children.  at each level,
c     for each flag++ box, test colleagues children and see
c     if they have children that are contacting you.  if so,
c     divide and flag++ all children that are created.     

      do i = 0, nlev
      ntemp = nboxes
      start = istartlev(i)
      istop  = istartlev(i) + nblevel(i) - 1
      do 500 j = start, istop
       ibox = iboxlev(j)
c      only be concerned with boxes on this level and
c      boxes that are given a flag++:
       if(iflag(ibox) .ne. 3)goto 500

         icoll(1) = icolleagbox(1,ibox)
         icoll(2) = icolleagbox(2,ibox)
         icoll(3) = icolleagbox(3,ibox)
         icoll(4) = icolleagbox(4,ibox)
         icoll(5) = icolleagbox(5,ibox)
         icoll(6) = icolleagbox(6,ibox)
         icoll(7) = icolleagbox(7,ibox)
         icoll(8) = icolleagbox(8,ibox)
         icoll(9) = icolleagbox(9,ibox)


c       scan colleagues.
        do 400 jcntr = 1, 9
        if(icoll(jcntr) .lt. 0)goto 400
        if(ichildbox(1,icoll(jcntr)) .lt. 0)goto 400

         ichild(1) = ichildbox(1,icoll(jcntr))
         ichild(2) = ichildbox(2,icoll(jcntr))
         ichild(3) = ichildbox(3,icoll(jcntr))
         ichild(4) = ichildbox(4,icoll(jcntr))


c          scan colleague's children.
           do 300 icntr = 1, 4
           if (ichildbox(1,ichild(icntr)) .lt. 0)goto 300 

           if(jcntr .eq. 1 .and. icntr .eq. 2)then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c          after calling subdivide or subdivide1, assign fval
           ik1=ichildbox(1,ibox)
           ik2=ichildbox(2,ibox)
           ik3=ichildbox(3,ibox)
           ik4=ichildbox(4,ibox)
           call chebyp2k(ndeg,fval1(1,ibox),fval1(1,ik1),fval1(1,ik2),
     1                 fval1(1,ik3),fval1(1,ik4))
           call chebyp2k(ndeg,fval2(1,ibox),fval2(1,ik1),fval2(1,ik2),
     1                 fval2(1,ik3),fval2(1,ik4))
c------------------------
         endif



c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 2 .and. 
     1        (icntr .eq. 1 .or. icntr .eq. 2))then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c          after calling subdivide or subdivide1, assign fval
           ik1=ichildbox(1,ibox)
           ik2=ichildbox(2,ibox)
           ik3=ichildbox(3,ibox)
           ik4=ichildbox(4,ibox)
           call chebyp2k(ndeg,fval1(1,ibox),fval1(1,ik1),fval1(1,ik2),
     1                 fval1(1,ik3),fval1(1,ik4))
c------------------------

         endif



c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 3 .and. icntr .eq. 1)then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c           after calling subdivide or subdivide1, assign fval
            ik1=ichildbox(1,ibox)
            ik2=ichildbox(2,ibox)
            ik3=ichildbox(3,ibox)
            ik4=ichildbox(4,ibox)
            call chebyp2k(ndeg,fval1(1,ibox),fval1(1,ik1),fval1(1,ik2),
     1                 fval1(1,ik3),fval1(1,ik4))
            call chebyp2k(ndeg,fval2(1,ibox),fval2(1,ik1),fval2(1,ik2),
     1                 fval2(1,ik3),fval2(1,ik4))
c------------------------

         endif



c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 4 .and. 
     1        (icntr .eq. 2 .or. icntr .eq. 3))then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c           after calling subdivide or subdivide1, assign fval
            ik1=ichildbox(1,ibox)
            ik2=ichildbox(2,ibox)
            ik3=ichildbox(3,ibox)
            ik4=ichildbox(4,ibox)
            call chebyp2k(ndeg,fval1(1,ibox),fval1(1,ik1),fval1(1,ik2),
     1                 fval1(1,ik3),fval1(1,ik4))
            call chebyp2k(ndeg,fval2(1,ibox),fval2(1,ik1),fval2(1,ik2),
     1                 fval2(1,ik3),fval2(1,ik4))
c------------------------

         endif


c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 6 .and. 
     1        (icntr .eq. 1 .or. icntr .eq. 4))then
             
             
c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c           after calling subdivide or subdivide1, assign fval
            ik1=ichildbox(1,ibox)
            ik2=ichildbox(2,ibox)
            ik3=ichildbox(3,ibox)
            ik4=ichildbox(4,ibox)
            call chebyp2k(ndeg,fval1(1,ibox),fval1(1,ik1),fval1(1,ik2),
     1                 fval1(1,ik3),fval1(1,ik4))
            call chebyp2k(ndeg,fval2(1,ibox),fval2(1,ik1),fval2(1,ik2),
     1                 fval2(1,ik3),fval2(1,ik4))
c------------------------

         endif


c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 7 .and. icntr .eq. 3)then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c           after calling subdivide or subdivide1, assign fval
            ik1=ichildbox(1,ibox)
            ik2=ichildbox(2,ibox)
            ik3=ichildbox(3,ibox)
            ik4=ichildbox(4,ibox)
            call chebyp2k(ndeg,fval1(1,ibox),fval1(1,ik1),fval1(1,ik2),
     1                 fval1(1,ik3),fval1(1,ik4))
            call chebyp2k(ndeg,fval2(1,ibox),fval2(1,ik1),fval2(1,ik2),
     1                 fval2(1,ik3),fval2(1,ik4))
c------------------------

         endif


c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 8 .and. 
     1        (icntr .eq. 3 .or. icntr .eq. 4))then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c           after calling subdivide or subdivide1, assign fval
            ik1=ichildbox(1,ibox)
            ik2=ichildbox(2,ibox)
            ik3=ichildbox(3,ibox)
            ik4=ichildbox(4,ibox)
            call chebyp2k(ndeg,fval1(1,ibox),fval1(1,ik1),fval1(1,ik2),
     1                 fval1(1,ik3),fval1(1,ik4))
            call chebyp2k(ndeg,fval2(1,ibox),fval2(1,ik1),fval2(1,ik2),
     1                 fval2(1,ik3),fval2(1,ik4))
c------------------------

         endif


c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3

           elseif(jcntr .eq. 9 .and. icntr .eq. 4)then

c           call subdivide
         if(ichildbox(1,ibox) .lt. 0)then
            call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev, iperiod,
     3         itemparray)
c------------------------
c           after calling subdivide or subdivide1, assign fval
            ik1=ichildbox(1,ibox)
            ik2=ichildbox(2,ibox)
            ik3=ichildbox(3,ibox)
            ik4=ichildbox(4,ibox)
            call chebyp2k(ndeg,fval1(1,ibox),fval1(1,ik1),fval1(1,ik2),
     1                 fval1(1,ik3),fval1(1,ik4))
            call chebyp2k(ndeg,fval2(1,ibox),fval2(1,ik1),fval2(1,ik2),
     1                 fval2(1,ik3),fval2(1,ik4))
c------------------------

         endif


c           flag++ all children created
            itemp = ichildbox(1,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(2,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(3,ibox)
            iflag(itemp) = 3
            itemp = ichildbox(4,ibox)
            iflag(itemp) = 3
           endif
300     continue   
400     continue
500     continue
      end do
      return
      end subroutine





      subroutine refine1(levelbox,icolbox,irowbox,nboxes,nlev,
     1         iparentbox,ichildbox,nblevel,iboxlev,istartlev)
      implicit real*8 (a-h,o-z)
      integer  nboxes, nlev
      integer  levelbox(1)
      integer  icolbox(1), irowbox(1)
      integer  ichildbox(4,1), iparentbox(1)
      integer  nblevel(0:nlev)
      integer  iboxlev(1),istartlev(0:nlev)
      integer, allocatable:: ichildold(:,:)
      integer, allocatable:: itemparray(:)
c
      allocate(itemparray(nboxes*3))
      allocate(ichildold(4,nboxes*3))
c
      nold=nboxes
      do i=1,nold
        do k=1,4
          ichildold(k,i)=ichildbox(k,i)
        enddo
      enddo
c
      do i=1,nold
        if(ichildold(1,i) .lt. 0) then
          iparbox=i
          call subdivide1(iparbox,iparentbox,ichildbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev,itemparray)
        endif
      enddo
c
      deallocate(itemparray)
      deallocate(ichildold)

      end subroutine






c**********************************************************
c     This subroutine returns the center of a given box in 
c     the tree 
c
c     INPUT:
c     xsize0: side length of the root box
c     icol: column index of the box
c     irow: row index of the box
c     level: level of the box
c
c     OUTPUT:
c     (cx, cy): coordinates of the center of the box
c
c**********************************************************
c
      subroutine get_cent_box(cx,cy,cent0,xsize0,icol,irow,level)
      implicit real*8 (a-h,o-z)
      integer icol, irow, level
      real*8 cx, cy, xsize0, cent0(2)
c
      xlength = xsize0 / dble(2**level)
      x0 = cent0(1)-xsize0/2.0d0
      y0 = cent0(2)-xsize0/2.0d0
      cx = x0 + (icol-0.5d0)*xlength
      cy = y0 + (irow-0.5d0)*xlength

      end subroutine








c**********************************************************
c
c     Given a tree and an array of points sorted
c     in the tree (leaf nodes), this subroutine returns the 
c     relative coordinates w.r.t. box centers
c
c     INPUT:
c     nlev - istartlev: the tree
c     xpts - ibofp: the points and point-tree relation
c
c     OUTPUT:
c     xrel: relative coordinate of each point w.r.t. the 
c           center of the leaf box it belongs to
c
c**********************************************************
c
      subroutine get_rel_coords(nlev,levelbox,iparentbox,
     1           ichildbox,icolbox,irowbox,nboxes,nblevel,
     2           iboxlev,istartlev,xpts,npts,npbox,ipold,
     3           istartbox,xrel,cent0,xsize0)
      implicit real*8 (a-h,o-z)
      integer nlev, nboxes, npts
      integer levelbox(nboxes), iparentbox(nboxes)
      integer ichildbox(4,nboxes)
      integer icolbox(nboxes), irowbox(nboxes)
      integer nblevel(0:nlev), iboxlev(nboxes)
      integer istartlev(0:nlev)
      integer npbox(nboxes), ipold(npts)
      integer istartbox(nboxes)
      real*8 xpts(2,npts), xrel(2,npts)
      real*8 cent0(2), xsize0
c
      do i=1,npts
        xrel(1,i)=0.0d0
        xrel(2,i)=0.0d0
      enddo
c     initialize to zero
c
      do i=1,nboxes
        ic=ichildbox(1,i)
        np=npbox(i)
        if((ic .lt. 0) .and. (np .gt. 0)) then
c       go through non-empty leaf boxes
          call get_cent_box(cx,cy,cent0,xsize0,icolbox(i),irowbox(i),
     1         levelbox(i)) 
          istart=istartbox(i)
          iend=istart+npbox(i)-1
          do ii=istart, iend
            xrel(1,ii)=xpts(1,ii)-cx
            xrel(2,ii)=xpts(2,ii)-cy
          enddo 
        endif
      enddo

      end subroutine







c**********************************************************
c
c     This subroutine computes the distance from the center
c     of a box (with a given side length) to the center of 
c     a colleague box's children
c
c     INPUT:
c     xsize: side length of the box
c     nb: neighbor index (1-9, nb=5 => itself)
c     ic: child index (1-4), numbered as follows
c
c                        1     2
c
c                        4     3
c
c     OUTPUT:
c     dx, dy: coords(center of neighbor's child) 
c            -coords(the box)
c
c**********************************************************
c
      subroutine get_dist_small(xsize0, level, nb, ic, dx, dy)
      implicit real*8 (a-h,o-z) 
      integer nb, ic, level
      real*8 xsize, dx, dy
      real*8 dist(2,4,9)
c
      xsize = xsize0 / dble(2**level)
c
      dist(1,1,1)=-1.250d0
      dist(2,1,1)=-0.750d0
      dist(1,2,1)=-0.750d0
      dist(2,2,1)=-0.750d0
      dist(1,3,1)=-0.750d0
      dist(2,3,1)=-1.250d0
      dist(1,4,1)=-1.250d0
      dist(2,4,1)=-1.250d0
      dist(1,1,2)=-0.250d0
      dist(2,1,2)=-0.750d0
      dist(1,2,2)=0.250d0
      dist(2,2,2)=-0.750d0
      dist(1,3,2)=0.250d0
      dist(2,3,2)=-1.250d0
      dist(1,4,2)=-0.250d0
      dist(2,4,2)=-1.250d0
      dist(1,1,3)=0.750d0
      dist(2,1,3)=-0.750d0
      dist(1,2,3)=1.250d0
      dist(2,2,3)=-0.750d0
      dist(1,3,3)=1.250d0
      dist(2,3,3)=-1.250d0
      dist(1,4,3)=0.750d0
      dist(2,4,3)=-1.250d0
      dist(1,1,4)=-1.250d0
      dist(2,1,4)=0.250d0
      dist(1,2,4)=-0.750d0
      dist(2,2,4)=0.250d0
      dist(1,3,4)=-0.750d0
      dist(2,3,4)=-0.250d0
      dist(1,4,4)=-1.250d0
      dist(2,4,4)=-0.250d0
      dist(1,1,5)=-0.250d0
      dist(2,1,5)=0.250d0
      dist(1,2,5)=0.250d0
      dist(2,2,5)=0.250d0
      dist(1,3,5)=0.250d0
      dist(2,3,5)=-0.250d0
      dist(1,4,5)=-0.250d0
      dist(2,4,5)=-0.250d0
      dist(1,1,6)=0.750d0
      dist(2,1,6)=0.250d0
      dist(1,2,6)=1.250d0
      dist(2,2,6)=0.250d0
      dist(1,3,6)=1.250d0
      dist(2,3,6)=-0.250d0
      dist(1,4,6)=0.750d0
      dist(2,4,6)=-0.250d0
      dist(1,1,7)=-1.250d0
      dist(2,1,7)=1.250d0
      dist(1,2,7)=-0.750d0
      dist(2,2,7)=1.250d0
      dist(1,3,7)=-0.750d0
      dist(2,3,7)=0.750d0
      dist(1,4,7)=-1.250d0
      dist(2,4,7)=0.750d0
      dist(1,1,8)=-0.250d0
      dist(2,1,8)=1.250d0
      dist(1,2,8)=0.250d0
      dist(2,2,8)=1.250d0
      dist(1,3,8)=0.250d0
      dist(2,3,8)=0.750d0
      dist(1,4,8)=-0.250d0
      dist(2,4,8)=0.750d0
      dist(1,1,9)=0.750d0
      dist(2,1,9)=1.250d0
      dist(1,2,9)=1.250d0
      dist(2,2,9)=1.250d0
      dist(1,3,9)=1.250d0
      dist(2,3,9)=0.750d0
      dist(1,4,9)=0.750d0
      dist(2,4,9)=0.750d0    
c
      dx = dist(1,ic,nb)*xsize
      dy = dist(2,ic,nb)*xsize
c
c
      end subroutine






c**********************************************************
c
c     This subroutine computes the distance from the center
c     of a box (with a given side length) to the center of 
c     a colleague box
c
c     the colleague boxes are indexed as follows
c
c                        7     8     9
c               
c                        4     5     6
c
c                        1     2     3
c
c
c**********************************************************
c
      subroutine get_dist_col(xsize0, level, nb, dx, dy)
      implicit real*8 (a-h,o-z)
      integer nb, ic, level
      real*8 xsize0, xsize, dx, dy
      real*8 dist(2,9)
c
      xsize = xsize0 / dble(2**level)
c
      dist(1,1)=-1.0d0
      dist(2,1)=-1.0d0
c
      dist(1,2)=0.0d0
      dist(2,2)=-1.0d0
c
      dist(1,3)=1.0d0
      dist(2,3)=-1.0d0
c
      dist(1,4)=-1.0d0
      dist(2,4)=0.0d0
c
      dist(1,5)=0.0d0
      dist(2,5)=0.0d0
c
      dist(1,6)=1.0d0
      dist(2,6)=0.0d0
c
      dist(1,7)=-1.0d0
      dist(2,7)=1.0d0
c
      dist(1,8)=0.0d0
      dist(2,8)=1.0d0
c
      dist(1,9)=1.0d0
      dist(2,9)=1.0d0
c
      dx=dist(1,nb)*xsize
      dy=dist(2,nb)*xsize
c

      end subroutine





c**********************************************************
c     subroutine mkcolls_shift
c**********************************************************
c     the following subroutine is used to generate the colleagues
c     for all of the boxes in the tree structure.  if a colleague
c     doesn't exist it is set to -1.  each box has nine colleagues
c     and they are ordered as follows:
c
c                        7     8     9
c               
c                        4     5     6
c
c                        1     2     3
c
c     you are your own colleague number 5.
c
c---------------
c     In the periodic case (iperiod = 1,2 for now)
c     if a box is in the colleague list in the sense
c     of a periodic image,
c     it also returns the 'shift' (in the unit of one
c     root box size) from the original box to the 
c     periodic image
c     update: 02-02-2019, J.W.
c---------------
c
c     the algorithm used here is recursive and takes advantage of
c     the fact that your colleagues can only be the children of 
c     your parents colleagues.  there is no need to scan all of the
c     boxes.  iperiod denotes whether or not we are in a periodic
c     or free space case.  the basic algorithm is the same, but in the
c     periodic case we have to look for boxes that are 'outside' of
c     the standard size box.
c
c
c     input:
c
c     levelbox is an array determining the level of each box
c
c     icolbox denotes the column of each box
c
c     irowbox denotes the row of each box
c
c     nboxes is the total number of boxes
c
c     nlev is the finest level
c
c     iparentbox denotes the parent of each box
c
c     ichildbox denotes the four children of each box
c
c     nblevel is the total number of boxes per level
c
c     iboxlev is the array in which the boxes are arranged
c
c     istartlev is the pointer to where each level begins in the
c               iboxlev array
c
c     iperiod denotes what kind of colleagues are to be generated
c             iperiod = 0 : free space
c             iperiod = 1 or 2 : periodic
c             iperiod = 3 : periodic up/down and free space left/right
c
c     output:
c
c     icolleagbox denotes the colleagues of a given box
c
c**********************************************************
c
      subroutine mkcolls_shift(icolbox,
     1      irowbox, icolleagbox, nboxes, nlev,
     2      iparentbox, ichildbox, nblevel,
     3      iboxlev, istartlev, iperiod, ishift)
      implicit none
c-----global variables
      integer icolleagbox(9,1)
      integer icolbox(1), irowbox(1)
      integer nboxes, nlev, iparentbox(1)
      integer ichildbox(4,1)
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer iperiod, ishift(2,9,1)
c-----local variables
      integer colleague, partemp
      integer jcntr, ibox, itest
      integer icntr, ilev, j, l, nside
      integer irowtemp, icoltemp
      integer irowtest, icoltest
      integer ishiftx, ishifty


c     initialize colleague number 5 to
c     yourself and all other colleagues to
c     -1.  -1 is the flag for the case when
c     the colleagues don't exist.  it can 
c     be overwritten below. 
      do ibox = 1, nboxes
       icolleagbox(5,ibox) = ibox
       do j = 1, 4
         icolleagbox(j,ibox) = -1
       end do
       do j = 6, 9
         icolleagbox(j,ibox) = -1
       end do
      end do
c
c     assign ishift to zeros
c     for the non-periodic case
      do ibox = 1, nboxes
        do j=1,9
          ishift(1,j,ibox)=0
          ishift(2,j,ibox)=0
        enddo
      enddo
c
c     scan through all of the levels except the coarsest level.
c     the one box at this level cannot have any colleagues.
c     do the uniform case first:
      if(iperiod .eq. 0)then
      do ilev = 1, nlev
c      scan through all of the boxes on each level.  for each test
c      box, scan the parent's colleagues and test to see if 
c      their children are in contact with the box being tested.
c      each colleague is placed in the correct order, so there is
c      no need to 'shuffle' them later on.
       do l = istartlev(ilev), istartlev(ilev) + nblevel(ilev) - 1
         ibox    = iboxlev(l)
         partemp = iparentbox(ibox)

c        irowtemp and icoltemp denote the row and column of
c        the test box.
         irowtemp = irowbox(ibox)
         icoltemp = icolbox(ibox)
      
         do 100 jcntr = 1, 9
c          colleague denotes the colleague of the parent box.
           colleague = icolleagbox(jcntr,partemp)
c          if the colleague doesn't exist
c          or is childless, skip it:
           if (colleague .lt. 0)goto 100
           if (ichildbox(1,colleague) .lt. 0)goto 100
           do icntr = 1, 4
             j = ichildbox(icntr,colleague)
c            irowtest and icoltest denote the row and column of
c            the box being compared to the test box.
             irowtest = irowbox(j)
             icoltest = icolbox(j)

             if(irowtemp .eq. irowtest+1)then
               if(icoltemp .eq. icoltest+1)then
                 icolleagbox(1,ibox) = j
               elseif(icoltemp .eq. icoltest)then
                 icolleagbox(2,ibox) = j
               elseif(icoltemp .eq. icoltest-1)then
                 icolleagbox(3,ibox) = j
               endif
             elseif(irowtemp .eq. irowtest)then
               if(icoltemp .eq. icoltest+1)then
                 icolleagbox(4,ibox) = j
               elseif(icoltemp .eq. icoltest-1)then
                 icolleagbox(6,ibox) = j
               endif
             elseif(irowtemp .eq. irowtest-1)then
               if(icoltemp .eq. icoltest+1)then
                 icolleagbox(7,ibox) = j
               elseif(icoltemp .eq. icoltest)then
                 icolleagbox(8,ibox) = j
               elseif(icoltemp .eq. icoltest-1)then
                 icolleagbox(9,ibox) = j
               endif
             endif
          end do 
100      continue
       end do
      end do 

c     now compute the colleagues in
c     the periodic case, if necessary:
      elseif(iperiod .eq. 1 .or. iperiod .eq. 2)then
c     initialize the first box (level 0) so
c     that it has its own colleagues.
c     this is necessary, because at deeper 
c     levels, the algorithm works by scanning 
c     the parent boxes colleagues.
      ibox = iboxlev(istartlev(0))
      icolleagbox(1,ibox) = ibox
      icolleagbox(2,ibox) = ibox
      icolleagbox(3,ibox) = ibox
      icolleagbox(4,ibox) = ibox
      icolleagbox(5,ibox) = ibox
      icolleagbox(6,ibox) = ibox
      icolleagbox(7,ibox) = ibox
      icolleagbox(8,ibox) = ibox
      icolleagbox(9,ibox) = ibox
c
      ishift(1,1,ibox)=-1
      ishift(2,1,ibox)=-1
c
      ishift(1,2,ibox)=0
      ishift(2,2,ibox)=-1
c
      ishift(1,3,ibox)=1
      ishift(2,3,ibox)=-1
c
      ishift(1,4,ibox)=-1
      ishift(2,4,ibox)=0
c
      ishift(1,5,ibox)=0
      ishift(2,5,ibox)=0
c
      ishift(1,6,ibox)=1
      ishift(2,6,ibox)=0
c
      ishift(1,7,ibox)=-1
      ishift(2,7,ibox)=1
c
      ishift(1,8,ibox)=0
      ishift(2,8,ibox)=1
c
      ishift(1,9,ibox)=1
      ishift(2,9,ibox)=1
c
      do ilev = 1, nlev
      nside = 2**ilev
c      scan through all of the boxes on each level.  for each test
c      box, scan the parent's colleagues and test to see if
c      their children are in contact with the box being tested.
c      each colleague is placed in the correct order, so there is
c      no need to 'shuffle' them later on.
       do l = istartlev(ilev), istartlev(ilev) + nblevel(ilev) - 1
        ibox = iboxlev(l)

c       irowtemp and icoltemp denote the
c       row and column of the test box.
        irowtemp = irowbox(ibox)
        icoltemp = icolbox(ibox)

c       irowtest and icoltest denote the row and column of
c       the box being compared to the test box.

        do 300 jcntr = 1, 9
c         first determine the column and row numbers
c         of all of the potential colleagues:
          if(jcntr .eq. 5)goto 300
          if(jcntr .eq. 1)then
            icoltest = icoltemp - 1
            irowtest = irowtemp - 1
          elseif(jcntr .eq. 2)then
            icoltest = icoltemp
            irowtest = irowtemp - 1
          elseif(jcntr .eq. 3)then
            icoltest = icoltemp + 1
            irowtest = irowtemp - 1
          elseif(jcntr .eq. 4)then
            icoltest = icoltemp - 1
            irowtest = irowtemp
          elseif(jcntr .eq. 6)then
            icoltest = icoltemp + 1
            irowtest = irowtemp
          elseif(jcntr .eq. 7)then
            icoltest = icoltemp - 1
            irowtest = irowtemp + 1
          elseif(jcntr .eq. 8)then
            icoltest = icoltemp
            irowtest = irowtemp + 1
          elseif(jcntr .eq. 9)then
            icoltest = icoltemp + 1
            irowtest = irowtemp + 1
          endif

c         now test to see if the test parameters 
c         lie in the domain:
c         (if they are outside of the domain, just
c         add or subtract the appropriate number so
c         that the boxes 'wrap around.')
          ishiftx=0
          ishifty=0
          if(icoltest .lt. 1)then
            icoltest = icoltest + nside
            ishiftx = -1
          elseif(icoltest .gt. nside)then
            icoltest = icoltest - nside
            ishiftx = 1
          endif
          if(irowtest .lt. 1)then
            irowtest = irowtest + nside
            ishifty = -1
          elseif(irowtest .gt. nside)then
            irowtest = irowtest - nside
            ishifty = 1
          endif


       do 200 j = 1, 9
        if(icolleagbox(j,iparentbox(ibox)) .lt. 0)goto 200
        if(ichildbox(1,icolleagbox(j,iparentbox(ibox))) .lt. 0)goto 200
          do icntr = 1, 4
            itest = ichildbox(icntr,icolleagbox(j,iparentbox(ibox)))
            if(irowbox(itest) .eq. irowtest
     1         .and. icolbox(itest) .eq. icoltest)then
               icolleagbox(jcntr,ibox) = itest
              ishift(1,jcntr,ibox)=ishiftx
              ishift(2,jcntr,ibox)=ishifty
            endif
          end do
200    continue
300    continue
      end do
      end do
      elseif(iperiod .eq. 3)then
c     initialize the first box (level 0) so
c     that it has its own colleagues.
c     this is necessary, because at deeper 
c     levels, the algorithm works by scanning 
c     the parent boxes colleagues.
      ibox = iboxlev(istartlev(0))
      icolleagbox(1,ibox) = ibox
      icolleagbox(2,ibox) = ibox
      icolleagbox(3,ibox) = ibox
      icolleagbox(4,ibox) = ibox
      icolleagbox(5,ibox) = ibox
      icolleagbox(6,ibox) = ibox
      icolleagbox(7,ibox) = ibox
      icolleagbox(8,ibox) = ibox
      icolleagbox(9,ibox) = ibox

      do ilev = 1, nlev
      nside = 2**ilev
c      scan through all of the boxes on each level.  for each test
c      box, scan the parent's colleagues and test to see if
c      their children are in contact with the box being tested.
c      each colleague is placed in the correct order, so there is
c      no need to 'shuffle' them later on.
       do l = istartlev(ilev), istartlev(ilev) + nblevel(ilev) - 1
        ibox = iboxlev(l)

c       irowtemp and icoltemp denote the
c       row and column of the test box.
        irowtemp = irowbox(ibox)
        icoltemp = icolbox(ibox)

c       irowtest and icoltest denote the row and column of
c       the box being compared to the test box.

        do 500 jcntr = 1, 9
c         first determine the column and row numbers
c         of all of the potential colleagues:
          if(jcntr .eq. 5)goto 500
          if(jcntr .eq. 1)then
            icoltest = icoltemp - 1
            irowtest = irowtemp - 1
          elseif(jcntr .eq. 2)then
            icoltest = icoltemp
            irowtest = irowtemp - 1
          elseif(jcntr .eq. 3)then
            icoltest = icoltemp + 1
            irowtest = irowtemp - 1
          elseif(jcntr .eq. 4)then
            icoltest = icoltemp - 1
            irowtest = irowtemp
          elseif(jcntr .eq. 6)then
            icoltest = icoltemp + 1
            irowtest = irowtemp
          elseif(jcntr .eq. 7)then
            icoltest = icoltemp - 1
            irowtest = irowtemp + 1
          elseif(jcntr .eq. 8)then
            icoltest = icoltemp
            irowtest = irowtemp + 1
          elseif(jcntr .eq. 9)then
            icoltest = icoltemp + 1
            irowtest = irowtemp + 1
          endif

c         now test to see if the test parameters 
c         lie in the domain:
c         (if they are outside of the domain, just
c         add or subtract the appropriate number so
c         that the boxes 'wrap around.')
          if(irowtest .lt. 1)then
            irowtest = irowtest + nside
          elseif(irowtest .gt. nside)then
            irowtest = irowtest - nside
          endif


       do 400 j = 1, 9
        if(icolleagbox(j,iparentbox(ibox)) .lt. 0)goto 400
        if(ichildbox(1,icolleagbox(j,iparentbox(ibox))) .lt. 0)goto 400
          do icntr = 1, 4
            itest = ichildbox(icntr,icolleagbox(j,iparentbox(ibox)))
            if(irowbox(itest) .eq. irowtest
     1         .and. icolbox(itest) .eq. icoltest)then
               icolleagbox(jcntr,ibox) = itest
            endif
          end do
400    continue
500    continue
      end do
      end do
      endif
      return
      end subroutine







c-----------------------------------------------------
c     subroutine interpv2p
c-----------------------------------------------------
c
c     This subroutine interpolates the function sampled on 
c     the Chebyshev grids of leaf nodes of an adaptive
c     quadtree to a given array of points
c
c     INPUT: 
c     levelbox - istartlev: the tree structure
c     ndeg: degree of function approx on leaf nodes
c     pot: the array of function values
c     npts: number of point targets
c     xpts: xy-coordinates of the points targets
c
c     OUTPUT:
c     potp: the array of function values at xpts
c
c-----------------------------------------------------
c
      subroutine interpv2b(levelbox, icolbox, irowbox, 
     1           nboxes, nlev, iparentbox, ichildbox,
     2           nblevel, iboxlev, istartlev, ndeg, pot,
     3           npts, xpts, potp, cent0, xsize0)
      implicit real*8 (a-h,o-z)
      integer nlev, nboxes, ndeg, npts
      integer levelbox(1)
      integer icolbox(1), irowbox(1)
      integer iparentbox(1), ichildbox(4,1)
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
      real*8 pot(ndeg*ndeg,nboxes)
      real*8 xpts(2,npts), potp(npts)
      real*8 cent0(2), xsize0
c     local vars
      integer,allocatable:: npbox(:), ipold(:)
      integer,allocatable:: istartbox(:), ibofp(:)
      real*8 wsave(1000), chwork(1000)
      real*8, allocatable:: chpot(:,:)
c
      allocate(chpot(ndeg,ndeg))
      allocate(npbox(nboxes))
      allocate(ipold(npts))
      allocate(istartbox(nboxes))
      allocate(ibofp(npts))
c
c     sort the grid points into the tree
      call treesort(nlev,levelbox,iparentbox,ichildbox,
     1    icolbox,irowbox,nboxes,nblevel,iboxlev,istartlev,
     2    xpts,npts,npbox,ipold,istartbox,ibofp,
     3    cent0,xsize0)
c
c     go through leaf nodes, 
c     if nonempty, do chebyshev interpolation
c
      call chxcin(ndeg, wsave)
      do ib=1, nboxes
        npb=npbox(ib)
        if((ichildbox(1,ib).lt.0) .and. (npb.gt.0)) then
          call chexfc2d(pot(1,ib),ndeg,chpot,wsave,chwork)
c
          lev=levelbox(ib)
          icol=icolbox(ib)
          irow=irowbox(ib)
  
          xlength=0.5d0**lev*xsize0
          x0=cent0(1)-xsize0/2.0d0
          y0=cent0(2)-xsize0/2.0d0
          a=dble(icol-1)*xlength+x0
          b=a+xlength
          c=dble(irow-1)*xlength+y0
          d=c+xlength
c           basic parameters
         
          istart=istartbox(ib) 
          iend=istart+npb-1
          do ip=istart, iend
            xx=xpts(1,ip) 
            yy=xpts(2,ip)
            call cheval2d(xx,yy,val,chpot,ndeg,a,b,c,d)
            iold=ipold(ip)
            potp(iold)=val
          enddo
        endif
      enddo
c
      call rvecnew2old(2, npts, xpts, ipold)
c     xpts reversed to the original order

      deallocate(npbox)
      deallocate(ipold)
      deallocate(istartbox)
      deallocate(ibofp)
      deallocate(chpot)

      end subroutine





c-----------------------------------------------------
c     subroutine spreadtree1
c-----------------------------------------------------
c     
c     given a quad-tree, the center and side 
c     length of the root box, spread out the 
c     tree by one level, so that it covers a
c     bigger region 
c
c-----------------------------------------------------
c
      subroutine spreadtree1(levelbox,icolbox,irowbox,
     1           nboxes,nlev,iparentbox,ichildbox,nblevel,
     2           iboxlev,istartlev,xsize0,cent0,idir,
     3           ndeg,fval)
      implicit real*8 (a-h,o-z)
c-----global variables
      integer nboxes, nlev, maxboxes, maxlevel
      integer maxppl, npts, ndeg
      integer levelbox(1)
      integer icolbox(1), irowbox(1)
      integer iparentbox(1), ichildbox(4,1)
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer idir(2)
      real*8 xsize0, cent0(2)
      real*8 fval(ndeg*ndeg,1)
      real*8, allocatable:: coeffp(:,:)
c
      allocate(coeffp(ndeg,ndeg))

      iroot=iboxlev(istartlev(0))
c     current root box
c
      if(abs(idir(1)*idir(2)) .gt. 1.0d-12) then

c     1. deal with level info: existing boxes
        do i=nboxes, 1, -1
          iboxlev(i+4)=iboxlev(i)
c         move back by 4 entries
          levelbox(i)=levelbox(i)+1
c         move down by 1 level
        enddo

        do l=nlev, 1, -1
          nblevel(l+1)=nblevel(l)
c         move back by 4 entries
          istartlev(l+1)=istartlev(l)+4
c         move down by 1 level
        enddo
c
c       2. deal with level info: new boxes
        do ibox=nboxes+1, nboxes+3
          levelbox(ibox)=1
        enddo
        levelbox(nboxes+4)=0
c       levels for the 4 new boxes
c
        nblevel(0)=1
        nblevel(1)=4
c
        istartlev(0)=1
        istartlev(1)=2
      
        do i=1,4
          iboxlev(i)=nboxes+5-i
        enddo
c
c       done with level structure, now deal with parent-child
c       and location info      
c
        icolbox(nboxes+4)=1
        irowbox(nboxes+4)=1
c
        iparentbox(nboxes+4)=-1
c
        iparentbox(iroot)=nboxes+4
        iparentbox(nboxes+3)=nboxes+4
        iparentbox(nboxes+2)=nboxes+4
        iparentbox(nboxes+1)=nboxes+4
c
        do i=1,4
          ichildbox(i,nboxes+1)=-1
          ichildbox(i,nboxes+2)=-1
          ichildbox(i,nboxes+3)=-1
        enddo
c
c  hand-rolled case select on idir, for the sake of readability
c
        if((idir(1) .gt. 0) .and. (idir(2).gt.0)) then
          cent0(1)=cent0(1)+xsize0/2.0d0
          cent0(2)=cent0(2)+xsize0/2.0d0
          icolbox(iroot)=1
          icolbox(nboxes+1)=1
          icolbox(nboxes+2)=2
          icolbox(nboxes+3)=2
c
          irowbox(iroot)=1
          irowbox(nboxes+1)=2
          irowbox(nboxes+2)=2
          irowbox(nboxes+3)=1
c
c       in this case, icolbox and irowbox remain the same
c       for all the other boxes
c
          ichildbox(1,nboxes+4)=nboxes+1
          ichildbox(2,nboxes+4)=nboxes+2
          ichildbox(3,nboxes+4)=nboxes+3
          ichildbox(4,nboxes+4)=iroot
c
        elseif((idir(1) .gt. 0).and.(idir(2).lt.0)) then
          cent0(1)=cent0(1)+xsize0/2.0d0
          cent0(2)=cent0(2)-xsize0/2.0d0
          icolbox(iroot)=1
          icolbox(nboxes+1)=2
          icolbox(nboxes+2)=2
          icolbox(nboxes+3)=1
c
          irowbox(iroot)=2
          irowbox(nboxes+1)=2
          irowbox(nboxes+2)=1
          irowbox(nboxes+3)=1
c
c       in this case, icolbox remains the same, while irowbox
c       is increased
c
          do l=2, nlev+1
            istart=istartlev(l)
            iend=istart+nblevel(l)-1
            ntot=2**(l-1)
            do ii=istart,iend
              ibox=iboxlev(ii)
              irowbox(ibox)=irowbox(ibox)+ntot
            enddo
          enddo
c
          ichildbox(1,nboxes+4)=iroot
          ichildbox(2,nboxes+4)=nboxes+1
          ichildbox(3,nboxes+4)=nboxes+2
          ichildbox(4,nboxes+4)=nboxes+3
        elseif((idir(1) .lt. 0).and.(idir(2).gt.0)) then
          cent0(1)=cent0(1)-xsize0/2.0d0
          cent0(2)=cent0(2)+xsize0/2.0d0
          icolbox(iroot)=2
          icolbox(nboxes+1)=1
          icolbox(nboxes+2)=2
          icolbox(nboxes+3)=1
c
          irowbox(iroot)=1
          irowbox(nboxes+1)=2
          irowbox(nboxes+2)=2
          irowbox(nboxes+3)=1
c
c       in this case, irowbox remains the same, while icolbox
c       is increased
c
          do l=2, nlev+1
            istart=istartlev(l)
            iend=istart+nblevel(l)-1
            ntot=2**(l-1)
            do ii=istart,iend
              ibox=iboxlev(ii)
              icolbox(ibox)=icolbox(ibox)+ntot
            enddo
          enddo
c
          ichildbox(1,nboxes+4)=nboxes+1
          ichildbox(2,nboxes+4)=nboxes+2
          ichildbox(3,nboxes+4)=iroot
          ichildbox(4,nboxes+4)=nboxes+3
        elseif((idir(1) .lt. 0).and.(idir(2).lt.0)) then
          cent0(1)=cent0(1)-xsize0/2.0d0
          cent0(2)=cent0(2)-xsize0/2.0d0
          icolbox(iroot)=2
          icolbox(nboxes+1)=1
          icolbox(nboxes+2)=2
          icolbox(nboxes+3)=1
c
          irowbox(iroot)=2
          irowbox(nboxes+1)=2
          irowbox(nboxes+2)=1
          irowbox(nboxes+3)=1
c
c       in this case, both irowbox and icolbox are increased
c
          do l=2, nlev+1
            istart=istartlev(l)
            iend=istart+nblevel(l)-1
            ntot=2**(l-1)
            do ii=istart,iend
              ibox=iboxlev(ii)
              irowbox(ibox)=irowbox(ibox)+ntot
              icolbox(ibox)=icolbox(ibox)+ntot
            enddo
          enddo
c
          ichildbox(1,nboxes+4)=nboxes+1
          ichildbox(2,nboxes+4)=iroot
          ichildbox(3,nboxes+4)=nboxes+2
          ichildbox(4,nboxes+4)=nboxes+3
        endif

c
c     done spreading the tree structure by 1 level
c     now interpolate the function value array fval
c
        do ibox=nboxes+1, nboxes+3
          do k=1,ndeg**2
            fval(k,ibox)=0.0d0
          enddo
        enddo
c
        k1=ichildbox(1,nboxes+4)
        k2=ichildbox(2,nboxes+4)
        k3=ichildbox(3,nboxes+4)
        k4=ichildbox(4,nboxes+4)
        call chebyk2p(ndeg, fval(1,k1), fval(1,k2), fval(1,k3),
     1        fval(1,k4), fval(1,nboxes+4), coeffp)
c
c
        xsize0=xsize0*2.0d0
        nlev=nlev+1
        nboxes=nboxes+4     
c
      endif
c
      deallocate(coeffp)

      end subroutine





c-----------------------------------------------------
c     subroutine spreadtree
c-----------------------------------------------------
c     
c     given a quad-tree, the center and side 
c     length of the root box, spread out the 
c     tree by two levels, so that it covers a
c     bigger region 
c
c     RMK: first spread out in directions idir1,
c          then spread out in the opposite directions
c
c-----------------------------------------------------
c
      subroutine spreadtree(levelbox,icolbox,irowbox,
     1           nboxes,nlev,iparentbox,ichildbox,nblevel,
     2           iboxlev,istartlev,xsize0,cent0,idir1,
     3           idir2,ndeg,fval)
      implicit real*8 (a-h,o-z)
c-----global variables
      integer nboxes, nlev, maxboxes, maxlevel
      integer maxppl, npts, ndeg
      integer levelbox(1)
      integer icolbox(1), irowbox(1)
      integer iparentbox(1), ichildbox(4,1)
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer idir1(2), idir2(2)
      real*8 xsize0, cent0(2)
      real*8 fval(ndeg*ndeg,1)
c
      if(idir1(1)*idir1(2) .ne. 0) then
        call spreadtree1(levelbox,icolbox,irowbox,
     1       nboxes,nlev,iparentbox,ichildbox,nblevel,
     2       iboxlev,istartlev,xsize0,cent0,idir1,
     3       ndeg,fval)
      endif
c
      if(idir2(1)*idir2(2) .ne. 0) then
        call spreadtree1(levelbox,icolbox,irowbox,
     1       nboxes,nlev,iparentbox,ichildbox,nblevel,
     2       iboxlev,istartlev,xsize0,cent0,idir2,
     3       ndeg,fval)
      endif
c
c
      end subroutine




c-----------------------------------------------------
c     subroutine checkedges
c-----------------------------------------------------
c
c     given an adapt tree, and function values given on 
c     the Chebyshev grids of leaf nodes, check function
c     values along the edges, and return the inf norm
c     on each of the edges
c
c     output:
c     nspread: number of spreads needed
c     idir1: direction of the first spreading
c     idir2: direction of the second spreading
c
c-----------------------------------------------------
c
      subroutine checkedges(levelbox, icolbox, irowbox, 
     1           nboxes, nlev, iparentbox, ichildbox,
     2           nblevel, iboxlev, istartlev, ndeg, fval,
     3           xsize0, cent0, nspread, idir1, idir2)
      implicit real*8 (a-h,o-z)
      integer nlev, nboxes, ndeg
      integer levelbox(1)
      integer icolbox(1), irowbox(1)
      integer iparentbox(1), ichildbox(4,1)
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer nspread, idir1(2), idir2(2)
      real*8 fval(ndeg*ndeg,nboxes), xsize0, cent0(2)
      real*8 fmaxedges(4)
      real*8, allocatable:: xg(:), yg(:), fg(:)
c
      nx=100
      ng=4*nx
      allocate(xg(ng))
      allocate(yg(ng))
      allocate(fg(ng))
c
      h=xsize0/(nx-1)
c
c     right edge
      do i=1,nx
        xg(i)=cent0(1)+xsize0/2.0d0
        yg(i)=cent0(2)-xsize0/2.0d0+h*(i-1)
      enddo
c
c     left edge
      do i=(nx+1),2*nx
        xg(i)=cent0(1)-xsize0/2.0d0
        yg(i)=cent0(2)-xsize0/2.0d0+h*(i-nx-1)
      enddo
c
c     top edge
      do i=(2*nx+1),3*nx
        xg(i)=cent0(1)-xsize0/2.0d0+h*(i-2*nx-1)
        yg(i)=cent0(2)+xsize0/2.0d0
      enddo
c
c     bottom edge
      do i=(3*nx+1),4*nx
        xg(i)=cent0(1)-xsize0/2.0d0+h*(i-3*nx-1)
        yg(i)=cent0(2)-xsize0/2.0d0
      enddo
c
c      do i=1,ng
c        write(23,*) xg(i), yg(i)
c      enddo
c
      call interppot(levelbox, icolbox, irowbox, 
     1     nboxes, nlev, iparentbox, ichildbox,
     2     nblevel, iboxlev, istartlev, ndeg, fval,
     3     ng, xg, yg, fg, cent0, xsize0)
c
      do i=1, 4
        istart=(i-1)*nx+1
        iend=i*nx
        fmax=0.0d0
        do ii=istart, iend
          if(abs(fg(ii)) .gt. fmax) then
            fmax=abs(fg(ii))
          endif
c          write(24,*) i, xg(ii), yg(ii), fg(ii)
        enddo
        fmaxedges(i)=fmax
c        write(*,*) i, fmaxedges(i)
      enddo
c
c     now that we've got fmaxedges, we can decide
c     the tree spreadings needed based on that
c
      if(fmaxedges(1) .gt. fmaxedges(2)) then
        idir1(1)=-1
        idir2(1)=1
      else
        idir1(1)=1
        idir2(1)=-1
      endif
c
c
      if(fmaxedges(3) .gt. fmaxedges(4)) then
        idir1(2)=-1
        idir2(2)=1
      else
        idir1(2)=1
        idir2(2)=-1
      endif
c
c----------------------------------------
c     
      tol=1.0d-11
c
      nzero=0
      do i=1,4
        if(fmaxedges(i) .lt. tol) then
          nzero=nzero+1
        endif
      enddo
c
      if(nzero .eq. 0) then
        nspread=2
      elseif(nzero .eq. 1) then
        nspread=2
      elseif(nzero .eq. 2) then
        fx=max(fmaxedges(1),fmaxedges(2))
        fy=max(fmaxedges(3),fmaxedges(4))
        if((fx .lt. tol) .or. (fy .lt. tol)) then
          nspread=2
        else
          nspread=1
          idir1(1)=idir2(1)
          idir1(2)=idir2(2)
          idir2(1)=0
          idir2(2)=0
        endif
      elseif(nzero .eq. 3) then
        nspread=1
        idir1(1)=idir2(1)
        idir1(2)=idir2(2)
        idir2(1)=0
        idir2(2)=0
      elseif(nzero .eq. 4) then
        nspread=0
        idir1(1)=0
        idir1(2)=0
        idir2(1)=0
        idir2(2)=0
      endif
c
c
      deallocate(xg)
      deallocate(yg)
      deallocate(fg)

      end subroutine




c-------------------------------------------------
c     a simple subroutine to get the number 
c     of leaf boxes in a tree
c-------------------------------------------------
      subroutine getnleaf(levelbox,icolbox,irowbox,
     1           nboxes,nlev,iparentbox,ichildbox, 
     2           nblevel,iboxlev,istartlev,nleaf)
      implicit none
      integer nboxes, nlev, nleaf
      integer levelbox(1)
      integer icolbox(1), irowbox(1)
      integer ichildbox(4,1), iparentbox(1)
      integer nblevel(0:nlev)
      integer iboxlev(1),istartlev(0:nlev)
      integer l, istart, iend, ii, ibox
c
      nleaf=0
c
      do l=0, nlev
        istart=istartlev(l)
        iend=istart+nblevel(l)-1
        do ii=istart, iend
          ibox=iboxlev(ii)
          if(ichildbox(1,ibox).lt.0) then
            nleaf=nleaf+1
          endif
        enddo
      enddo

      end subroutine






c-------------------------------------------------
c     subroutine tree2targs
c-------------------------------------------------
c
c     this subroutine assembles grid points on the 
c     leaf nodes of the tree as target points, 
c     return the array of points, as well as maps 
c     between the two kinds of indices
c
c     input:
c     levelbox - xsize0: the tree structure
c     
c     output:
c     ntarg: number of targets
c     xtarg: coordinates of the targets
c     itarg2tree: the map from target index to tree
c                 index
c     itree2targ: the map from tree index to target
c                 index 
c
c-------------------------------------------------
c
      subroutine tree2targs(levelbox, icolbox, irowbox, 
     1           nboxes, nlev, iparentbox, ichildbox,
     2           nblevel, iboxlev, istartlev, maxboxes,
     3           maxlevel, cent0, xsize0, ndeg, 
     4           ntarg, xtarg, itarg2tree, itree2targ)
      implicit real*8 (a-h,o-z)
      integer levelbox(1), maxboxes
      integer nlev, nboxes, maxlevel, ndeg
      integer icolbox(1), irowbox(1)
      integer iparentbox(1), ichildbox(4,1)
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer ntarg, itarg2tree(2,1)
      integer itree2targ(ndeg*ndeg,maxboxes)
      real *8 cent0(2), xsize0, xf(ndeg), yf(ndeg)
      real*8 xtarg(2,1)
c
      ntarg=0
      do ibox=1,nboxes
        if(ichildbox(1,ibox) .lt. 0) then
          call mkgrid(xf,yf,icolbox(ibox),irowbox(ibox),
     1         levelbox(ibox),ndeg,cent0,xsize0)
          do k2=1,ndeg
          do k1=1,ndeg
            ntarg=ntarg+1
            kk=(k2-1)*ndeg+k1
            xtarg(1,ntarg)=xf(k1)
            xtarg(2,ntarg)=yf(k2)
            itree2targ(kk,ibox)=ntarg
            itarg2tree(1,ntarg)=kk
            itarg2tree(2,ntarg)=ibox
          enddo
          enddo
        endif
      enddo


      end subroutine






c********************************************************
c     subroutine chebyvol2p
c********************************************************
c
c     chebyshev interpolation of a given function from
c     the tensor product chebyshev grid on the leaf nodes
c     of a given tree to a given set of points
c
c     INPUT:
c     (tree structure:)
c     istartlev is the pointer to where each level
c               begins in the iboxlev array
c     levelbox is an array determining the level of each box
c     nboxes is the total number of boxes
c     nlev is the finest level
c     icolbox denotes the column of each box
c     irowbox denotes the row of each box
c     iparentbox denotes the parent of each box
c     ichildbox denotes the four children of each box
c     nblevel is the total number of boxes per level
c     iboxlev is the array in which the boxes are arranged
c
c     ndeg: degree of chebyshev approx      
c     fval: array of dim(ndeg*ndeg,nboxes), function values
c           on the leaf boxes
c     npts: number of points
c     xpts: coordinates of points
c
c     (tree-point relation:)
c     npbox: number of sources in each box 
c                     (including non-leaf boxes)
c     ipold: ids of points on input
c            (before the overwriting)
c     istartbox: starting point in the array xpts
c                for each box
c     ibofp: indices of leaf boxes that the point
c            belongs to (after overwriting)
c
c     OUTPUT:
c     fpts: array of length npts, function values fval
c           interpolated to the points xpts
c
c     RMK: not tested yet
c     RMK: which order is better for the output?
c          the order of xpts or ipold? or allow for 
c          both?
c     RMK: ever called?
c
c********************************************************
c
      subroutine chebyvol2p(levelbox,icolbox,irowbox,
     1           nboxes,nlev,iparentbox,ichildbox,nblevel,
     2           iboxlev,istartlev,ndeg,fval,npts,xpts,
     3           npbox,ipold,istartbox,ibofp,fpts)
      implicit real*8 (a-h,o-z)
c-----global variables
      integer nboxes, nlev, npts, ndeg
      integer levelbox(1)
      integer icolbox(1), irowbox(1)
      integer iparentbox(1), ichildbox(4,1)
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer npbox(1), ipold(1), istartbox(1), ibofp(1)
      real*8 fval(ndeg*ndeg,nboxes)
      real*8 xpts(2,npts), fpts(npts)
c-----global variables
      real*8, allocatable:: wsave(:), chwork(:)
      real*8, allocatable:: chcoef(:,:)
c
c     go through boxes
c     if (nonempty leaf box)
c        interpolate to the points
c     endif      
c
      allocate(wsave(1000))
      allocate(chwork(1000))
      allocate(chcoef(ndeg,ndeg))
c
      call chxcin(ndeg, wsave)            
      do l=0, nlev
        istart=istartlev(l)
        iend=istart+nblevel(l)-1
        do i=istart, iend
          ibox=iboxlev(i)
          npb=npbox(ibox)
          if((ichildbox(1,ibox).lt.0) .and. (npb .gt. 0)) then
            call chexfc2d(fval(1,ibox),ndeg,chcoef,wsave,chwork) 
c           convert the point values on the grid to 
c           Chebyshev coefficients 
c
            lev=levelbox(ibox)
            icol=icolbox(ibox)
            irow=irowbox(ibox)
c
            xlength=0.5d0**lev
            a=dble(icol-1)*xlength-0.5d0
            b=a+xlength
            c=dble(irow-1)*xlength-0.5d0
            d=c+xlength
c
            jstart=istartbox(ibox)
            jend=jstart+npbox(ibox)-1
            do j=jstart, jend
              call cheval2d(xpts(1,j),xpts(2,j),val,chcoef,
     1             ndeg,a,b,c,d)
              fpts(j)=val
c             RMK: which order is better?
c                  the order of xpts or ipold?
            enddo

          endif
        enddo 
      enddo
c
c
      deallocate(wsave)
      deallocate(chwork)
      deallocate(chcoef)
c
      end subroutine






c********************************************************
c     suborutine restore
c********************************************************
c     a hacking routine, to be eliminated later
c
c     input:
c     npts, xpts, ibofs, ipold
c
c     output:
c     issorted
c     xpts and ibofs rewritten
c********************************************************
c
      subroutine restore(npts,xpts,ibofs,ipold,issorted)
      implicit real*8 (a-h,o-z)
      integer npts, ibofs(npts)
      integer ipold(npts), issorted(npts)
      real*8 xpts(2,npts)
      integer, allocatable:: ibofs0(:)
      real*8, allocatable:: xpts0(:,:)
c
      allocate(ibofs0(npts))
      allocate(xpts0(2,npts))
c
      do i=1,npts
        ibofs0(i)=ibofs(i)
        xpts0(1,i)=xpts(1,i)
        xpts0(2,i)=xpts(2,i)
      enddo
c     store the new order
c
      do i=1,npts
        iold=ipold(i)
        issorted(i)=iold
        ibofs(iold)=ibofs0(i)
        xpts(1,iold)=xpts0(1,i)
        xpts(2,iold)=xpts0(2,i)
      enddo
c     recover the original order
c
      deallocate(ibofs0)
      deallocate(xpts0)
c
      end subroutine





c********************************************************
c     the following subroutine generates the x and y values on which a 
c     function is defined, given the tree structure as input.
c
c     input:
c
c     levelbox is an array determining the level of each box
c
c     nboxes is the total number of boxes
c
c     icolbox denotes the column of each box
c
c     irowbox denotes the row of each box
c
c     ichildbox denotes the four children of each box
c
c     nlev is the finest level
c
c     nblevel is the total number of boxes per level
c
c     iboxlev is the array in which the boxes are arranged
c
c     istartlev is the pointer to where each level begins in the
c               iboxlev array
c
c     output:
c
c     xf denotes the x values on the leaf nodes
c
c     yf denotes the y values on the leaf nodes
c
c     note the arrays xf and yf are defined on childless boxes
c     exactly the same way as the right hand side.  the 36 points
c     on each node are numbered the same way.
c
c**************************************************************************
      subroutine getxyclassical8(xf, yf,
     1       icolbox, irowbox, ichildbox,nlev,
     2       nblevel, iboxlev, istartlev,
     3       cent0, xsize0)
      implicit none
c-----global variables
      integer *4  icolbox(1), irowbox(1), nlev
      integer *4  ichildbox(4,1)
      integer *4  nblevel(0:1), iboxlev(1), istartlev(0:1)
      real *8  xf(64,1), yf(64,1)
      real *8 cent0(2), xsize0, x0, y0
c-----local variables
      integer *4  i, ibox, j, k, l
      real *8  temp1
      real *8  xstart, pi16
      real *8  xx(8), xshift, yshift
      real *8  xscale(8)

      pi16 = dacos(-1.0d0) / 16.0d0
      xx(1) = dcos(15.0d0*pi16) / 2.0d0
      xx(2) = dcos(13.0d0*pi16) / 2.0d0
      xx(3) = dcos(11.0d0*pi16) / 2.0d0
      xx(4) = dcos( 9.0d0*pi16) / 2.0d0
      xx(5) = dcos( 7.0d0*pi16) / 2.0d0
      xx(6) = dcos( 5.0d0*pi16) / 2.0d0
      xx(7) = dcos( 3.0d0*pi16) / 2.0d0
      xx(8) = dcos( 1.0d0*pi16) / 2.0d0


      temp1 = xsize0 
      x0 = cent0(1)-xsize0/2.0d0
      y0 = cent0(2)-xsize0/2.0d0
c
      do k = 0, nlev
      xscale(1) = xx(1) * temp1 
      xscale(2) = xx(2) * temp1 
      xscale(3) = xx(3) * temp1 
      xscale(4) = xx(4) * temp1 
      xscale(5) = xx(5) * temp1 
      xscale(6) = xx(6) * temp1 
      xscale(7) = xx(7) * temp1 
      xscale(8) = xx(8) * temp1 

      do 100 i = istartlev(k), istartlev(k) + nblevel(k) - 1
        ibox = iboxlev(i)
        if(ichildbox(1,ibox) .gt. 0)goto 100

        xshift  =  x0+dble(icolbox(ibox)-0.5d0) * temp1
        yshift  =  y0+dble(irowbox(ibox)-0.5d0) * temp1

        do j = 1, 8
          do l = 1, 8
            xf(8*(l-1)+j,ibox) = xscale(j) + xshift
            yf(8*(j-1)+l,ibox) = xscale(j) + yshift
          end do
        end do
100   continue
      temp1 = temp1/2.0d0
      end do
      return
      end


      subroutine getxypractical8(xfp, yfp,
     1       icolbox, irowbox, ichildbox,nlev,
     2       nblevel, iboxlev, istartlev,
     3       cent0, xsize0)
      implicit none
c-----global variables
      integer *4  icolbox(1), irowbox(1), nlev
      integer *4  ichildbox(4,1)
      integer *4  nblevel(0:1), iboxlev(1), istartlev(0:1)
      real *8  xfp(81,1), yfp(81,1)
      real *8 cent0(2), xsize0, x0, y0
c-----local variables
      integer *4  i, ibox, j, k, l
      real *8  temp1
      real *8  xstart, pi16
      real *8  xx(9), xshift, yshift
      real *8  xscale(9)

      pi16 = dacos(-1.0d0) / 16.0d0
      xx(1) = -1.0d0 / 2.0d0
      xx(2) = dcos(14.0d0*pi16) / 2.0d0
      xx(3) = dcos(12.0d0*pi16) / 2.0d0
      xx(4) = dcos(10.0d0*pi16) / 2.0d0
      xx(5) = 0.0d0
      xx(6) = dcos( 6.0d0*pi16) / 2.0d0
      xx(7) = dcos( 4.0d0*pi16) / 2.0d0
      xx(8) = dcos( 2.0d0*pi16) / 2.0d0
      xx(9) = 1.0d0 / 2.0d0


      temp1 = xsize0
      do k = 0, nlev

      xscale(1) = xx(1) * temp1 
      xscale(2) = xx(2) * temp1 
      xscale(3) = xx(3) * temp1 
      xscale(4) = xx(4) * temp1 
      xscale(5) = xx(5) * temp1 
      xscale(6) = xx(6) * temp1 
      xscale(7) = xx(7) * temp1 
      xscale(8) = xx(8) * temp1 
      xscale(9) = xx(9) * temp1 

      do 100 i = istartlev(k), istartlev(k) + nblevel(k) - 1
        ibox = iboxlev(i)
        if(ichildbox(1,ibox) .gt. 0)goto 100

        xshift  =  x0+dble(icolbox(ibox)-0.5d0) * temp1
        yshift  =  y0+dble(irowbox(ibox)-0.5d0) * temp1

        do j = 1, 9
          do l = 1, 9
            xfp(9*(l-1)+j,ibox) = xscale(j) + xshift
            yfp(9*(j-1)+l,ibox) = xscale(j) + yshift
          end do
        end do
100   continue
      temp1 = temp1/2.0d0
      end do
      return
      end

c -----------------------------------------------------
c     subroutine adaptreef2d
c-----------------------------------------------------
c
c     Given a tree, and a function resolved on the 
c     leaf nodes, and another function (by function
c     evaluation routine), refine and coarsen the tree
c     so that both functions are resolved but not
c     over-resolved on the new tree
c
c     INPUT:
c     maxboxes: the maximum number of boxes allowed
c     maxlevel: deepest level allowed
c     levelbox - istartlev: the quad-tree structure
c     ndeg: degree of chebyshev approx on each leaf node
c     uval: array of function values (resolved already)
c     feval1,feval2: function evaluation routine of another
c            function
c     eps: error tolerance
c     t: parameter (e.g. time) needed in feval
c
c     OUTPUT:
c     the tree: overwritten to reflect the change
c     uval1,uval2: interpolated to the new tree
c     
c     local variables
c     valk1: function values on kid 1
c     valk2: function values on kid 2
c     valk3: function values on kid 3
c     valk4: function values on kid 4

      subroutine adaptreef2d(maxboxes, maxlevel, levelbox,
     1           icolbox, irowbox, nboxes, nlev, iparentbox,
     2           ichildbox, nblevel, iboxlev, istartlev,
     3           ndeg, uval1,uval2,feval1,feval2, t,
     4           eps, cent0, xsize0)

      implicit real*8 (a-h,o-z)
c-----global variables
      integer maxboxes, maxlevel
      integer nboxes, nlev, ndeg
      integer levelbox(1)
      integer icolbox(1), irowbox(1)
      integer iparentbox(1), ichildbox(4,1)
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
      real*8 uval1(ndeg*ndeg,1),uval2(ndeg*ndeg,1) 
      real*8 eps
      real*8 cent0(2), xsize0
      external feval1
      external feval2
c-----local variables
      integer kids(4)
      integer, allocatable:: iempty(:)
      integer, allocatable:: itemp(:), itemp4(:,:)
      integer, allocatable:: idold2new(:)
      integer, allocatable:: itemparray(:)
      real*8 wsave(1000), chwork(1000)
      real*8 xf(ndeg), yf(ndeg)
      real*8 ftemp1(ndeg,ndeg),ftemp2(ndeg,ndeg)
      real*8 coefftemp1(ndeg,ndeg),coefftemp2(ndeg,ndeg)
      real*8, allocatable:: valk1(:,:), valk2(:,:)
      real*8, allocatable:: valk3(:,:), valk4(:,:)
      real*8, allocatable:: valk1b(:,:), valk2b(:,:)
      real*8, allocatable:: valk3b(:,:), valk4b(:,:)
      
      allocate(itemparray(maxboxes))
      allocate(iempty(maxboxes))
      allocate(itemp(maxboxes))
      allocate(itemp4(4,maxboxes))
      allocate(idold2new(maxboxes))
c    demension1
      allocate(valk1(ndeg,ndeg))
      allocate(valk2(ndeg,ndeg))
      allocate(valk3(ndeg,ndeg))
      allocate(valk4(ndeg,ndeg))
c    demension2
      allocate(valk1b(ndeg,ndeg))
      allocate(valk2b(ndeg,ndeg))
      allocate(valk3b(ndeg,ndeg))
      allocate(valk4b(ndeg,ndeg))

c---------------------------------------------------------
c     1) refinement sweep. go down the tree level by level
c        refine leaf boxes if feval(:,t) is unresolved
c        no need to worry about uval, since by assumption
c        it is resolved
c---------------------------------------------------------
c
      ix=1
      iy=1
c      xlength=1.0d0
      xlength=xsize0
      x0=cent0(1)-xsize0/2.0d0
      y0=cent0(2)-xsize0/2.0d0
      nside=1
c
      call chxcin(ndeg,wsave)
      do i=0, maxlevel-1
        iflag=0
c          iflag: if the last level has been divided at all
        xlength=xlength/2.0d0
        nside=nside*2
c
        istart=istartlev(i)
        iend=istart+nblevel(i)-1

        do j=istart, iend
          ibox=iboxlev(j)
          if(ichildbox(1,ibox) .lt. 0) then

            call mkgrid(xf,yf,icolbox(ibox),irowbox(ibox),
     1           levelbox(ibox),ndeg,cent0,xsize0)

            do jj = 1, ndeg
              do ii = 1, ndeg
                 call feval1(xf(jj), yf(ii), t, ftemp1(jj,ii))
                 call feval2(xf(jj), yf(ii), t, ftemp2(jj,ii))
              enddo
            enddo
c demension1
            call getcoeff(ndeg,ndeg,ftemp1,coefftemp1,wsave)
            call geterror(ndeg,coefftemp1,error1)
c demension2
            call getcoeff(ndeg,ndeg,ftemp2,coefftemp2,wsave)
            call geterror(ndeg,coefftemp2,error2)

            hh = dble(4**levelbox(ibox))
            epscaled = eps * hh
            if(error1 .ge. epscaled .or. error2 .ge. epscaled)then
c             call subdivide 
              call subdivide1(ibox,iparentbox,ichildbox,
     1             nboxes,irowbox,icolbox,levelbox,nlev,
     2             istartlev, nblevel, iboxlev,itemparray)
              
c             need to assign function values to 
              ik1=ichildbox(1,ibox)
              ik2=ichildbox(2,ibox)
              ik3=ichildbox(3,ibox)
              ik4=ichildbox(4,ibox)
              call chebyp2k(ndeg,uval1(1,ibox),uval1(1,ik1),
     1             uval1(1,ik2),uval1(1,ik3),uval1(1,ik4))
              call chebyp2k(ndeg,uval2(1,ibox),uval2(1,ik1),
     1             uval2(1,ik2),uval2(1,ik3),uval2(1,ik4))
c
              iflag = 1
            endif
          endif

        enddo
c
        if((iflag .eq. 0) .and. (i .eq. nlev))then
          goto 110
        endif
      enddo

c------------------------------------------
c 2. go up the tree to coarsen:
c    if on a certain parent box, both functions
c    are resolved, merge its children
c------------------------------------------
c

110   do ii=1,nboxes
        iempty(ii)=0
      enddo   
c
      do i=nlev,0,-1
c       to merge, go up the tree level by level
        istart=istartlev(i)
        iend=istart+nblevel(i)-1
        do j=istart, iend
          ibox=iboxlev(j)
          nleafkids=0
          if(ichildbox(1,ibox) .gt. 0) then
            kids(1)=ichildbox(1,ibox)
            kids(2)=ichildbox(2,ibox)
            kids(3)=ichildbox(3,ibox)
            kids(4)=ichildbox(4,ibox)
c     get nleafkids
            do k=1,4
              if(ichildbox(1,kids(k)).lt.0) then
                nleafkids=nleafkids+1
              endif
            enddo
          endif
c
          if(nleafkids .eq. 4) then
            call mkgrid(xf,yf,icolbox(ibox),irowbox(ibox),
     1           levelbox(ibox),ndeg,cent0,xsize0)
            do jj = 1, ndeg
              do ii = 1, ndeg
                  call feval1(xf(jj), yf(ii), t, ftemp1(jj,ii))
                  call feval2(xf(jj), yf(ii), t, ftemp2(jj,ii))
              enddo
            enddo
c     get ferror
            call getcoeff(ndeg,ndeg,ftemp1,coefftemp1,wsave)
            call geterror(ndeg,coefftemp1,ferror1)

            call getcoeff(ndeg,ndeg,ftemp2,coefftemp2,wsave)
            call geterror(ndeg,coefftemp2,ferror2)

            do k2=1,ndeg
            do k1=1,ndeg
              k=(k2-1)*ndeg+k1

              valk1(k1,k2)=uval1(k,kids(1))
              valk2(k1,k2)=uval1(k,kids(2))
              valk3(k1,k2)=uval1(k,kids(3))
              valk4(k1,k2)=uval1(k,kids(4))

              valk1b(k1,k2)=uval2(k,kids(1))
              valk2b(k1,k2)=uval2(k,kids(2))
              valk3b(k1,k2)=uval2(k,kids(3))
              valk4b(k1,k2)=uval2(k,kids(4))
            enddo
            enddo

            call chebyk2p(ndeg, valk1, valk2, valk3, valk4,
     1           ftemp1, coefftemp1)
            call geterror(ndeg,coefftemp1,uerror1)

            call chebyk2p(ndeg, valk1b, valk2b, valk3b, valk4b,
     1           ftemp2, coefftemp2)
            call geterror(ndeg,coefftemp2,uerror2)

            hh = dble(4**levelbox(ibox))
            epscaled = eps * hh

            if((uerror1 .lt. epscaled .and.uerror2 .lt. epscaled) .and.
     1         (ferror1 .lt. epscaled .and. ferror2 .lt. epscaled)) then
c             the condition is correct here!
              call coarsen1(ibox,iparentbox,ichildbox,
     1             nboxes,irowbox,icolbox,levelbox,nlev,
     2             istartlev,nblevel,iboxlev,iempty)
              do k2=1,ndeg
              do k1=1,ndeg
                k=(k2-1)*ndeg+k1
                uval1(k,ibox)=ftemp1(k1,k2)
                uval2(k,ibox)=ftemp2(k1,k2)
              enddo
              enddo
c             assign interpolated function values
c             to the box that has become a leaf box
            endif
          endif
        enddo
      enddo



c------------------------------------------
c     3) clean up the intermediate tree with deleted
c        boxes, reorder uval at the same time
c------------------------------------------
c

      call reorder(iparentbox, ichildbox, nboxes,
     1     irowbox, icolbox, levelbox, nlev,
     2     istartlev, nblevel, iboxlev, iempty,
     3     idold2new, nold, itemp, itemp4)
c
      call reorderf(ndeg*ndeg, nold, uval1, idold2new)
      call reorderf(ndeg*ndeg, nold, uval2, idold2new)

      deallocate(itemparray)
      deallocate(iempty)
      deallocate(itemp)
      deallocate(itemp4)
      deallocate(idold2new)
c
      deallocate(valk1)
      deallocate(valk2)
      deallocate(valk3)
      deallocate(valk4)

      deallocate(valk1b)
      deallocate(valk2b)
      deallocate(valk3b)
      deallocate(valk4b)



      end subroutine


c
      subroutine mktreef2d(levelbox, icolbox, irowbox,
     1           nboxes, nlev, iparentbox, ichildbox,
     2           nblevel, iboxlev, istartlev,
     3           maxboxes, itemparray, maxlevel, eps, h1,h2,
     4           ndeg, cent0, xsize0,t)
c-----global variables
      implicit none
      integer *4  levelbox(1), maxboxes
      integer *4  nlev, nboxes,  maxlevel
      integer *4  icolbox(1), irowbox(1)
      integer *4  iparentbox(1), ichildbox(4,1)
      integer *4  nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer *4  itemparray(1), ndeg
      real *8  eps, h1,h2
      real *8 cent0(2), xsize0
c-----local variables
      integer *4  i, ibox, iflag
      integer *4  j, istart, iend, jj
      integer *4  levflag, ii
      real *8  coefftemp1(ndeg,ndeg)
      real *8  coefftemp2(ndeg,ndeg)
      real *8  epscaled
      real *8  error1,error2, hh
      real *8  xf(ndeg), yf(ndeg)
      real *8  ftemp1(ndeg,ndeg)
      real *8  ftemp2(ndeg,ndeg)
      real *8  wsave2(1000)
      real *8  wsave3(1000)
      real *8  t

      do i = 0, maxlevel
        nblevel(i) = 0
        istartlev(i) = 0
      end do
      do i = 1, maxboxes
        iboxlev(i) = 0
      end do
c
c     first set the big parent box to the
c     appropriate settings:
c     (initially, there is just one box and
c     it is the big parent box at level 0)
c
      ibox = 1
      levelbox(ibox) = 0
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = -1
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1

      nboxes = 1
      nlev = 0

c     we also need to initialize the adaptive 'ladder'
c     structures to the correct initial values:
      nblevel(0) = 1
      istartlev(0) = 1
      iboxlev(1) = 1

      call chxcin(ndeg,wsave2)
      call chxcin(ndeg,wsave3)
c
      do i = 0, maxlevel - 1
      iflag = 0
      istart = istartlev(i)
      iend = istart + nblevel(i) - 1
      do j = istart, iend
       ibox = iboxlev(j)

       call mkgrid(xf,yf,icolbox(ibox),irowbox(ibox),
     1      levelbox(ibox),ndeg,cent0,xsize0)

       do jj = 1, ndeg
         do ii = 1, ndeg
           ftemp1(jj,ii) = h1(xf(jj),yf(ii),t)
           ftemp2(jj,ii) = h2(xf(jj),yf(ii),t)
         end do
       end do


c      compute chebyshev transforms
       call getcoeff(ndeg,ndeg,ftemp1,coefftemp1,wsave2)
       call getcoeff(ndeg,ndeg,ftemp2,coefftemp2,wsave3)
       call geterror(ndeg,coefftemp1,error1)
       call geterror(ndeg,coefftemp2,error2)

       hh = dble(4**levelbox(ibox))
       epscaled = eps * hh

c       write(*,*) error, i

       if(error1 .ge. epscaled .or. error2 .ge. epscaled)then
c        call subdivide
         call subdivide1(ibox,iparentbox,ichildbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev,itemparray)
         iflag = 1
       endif
      end do
      if(iflag .eq. 0)then
c      nothing was divided at the
c      last level, so exit the loop.
       return
      endif
      end do
      return
      end subroutine


c -----------------------------------------------------
c     subroutine adaptree2
c-----------------------------------------------------
c
c     Given a tree, and a function resolved on the
c     leaf nodes, and another function (by function
c     evaluation routine), refine and coarsen the tree
c     so that both functions are resolved but not
c     over-resolved on the new tree
c
c     INPUT:
c     maxboxes: the maximum number of boxes allowed
c     maxlevel: deepest level allowed
c     levelbox - istartlev: the quad-tree structure
c     ndeg: degree of chebyshev approx on each leaf node
c     uval: array of function values (resolved already)
c     feval1,feval2: function evaluation routine of another
c            function
c     eps: error tolerance
c     t: parameter (e.g. time) needed in feval
c
c     OUTPUT:
c     the tree: overwritten to reflect the change
c     uval1,uval2: interpolated to the new tree
c
c     local variables
c     valk1: function values on kid 1
c     valk2: function values on kid 2
c     valk3: function values on kid 3
c     valk4: function values on kid 4

      subroutine adaptree2(maxboxes, maxlevel, levelbox,
     1           icolbox, irowbox, nboxes, nlev, iparentbox,
     2           ichildbox, nblevel, iboxlev, istartlev,
     3           ndeg, uval1,uval2,
     4           eps, cent0, xsize0)

      implicit real*8 (a-h,o-z)
c-----global variables
      integer maxboxes, maxlevel
      integer nboxes, nlev, ndeg
      integer levelbox(1)
      integer icolbox(1), irowbox(1)
      integer iparentbox(1), ichildbox(4,1)
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
      real*8 uval1(ndeg*ndeg,1),uval2(ndeg*ndeg,1)
      real*8 eps
      real*8 cent0(2), xsize0
      external feval1
      external feval2
c-----local variables
      integer kids(4)
      integer, allocatable:: iempty(:)
      integer, allocatable:: itemp(:), itemp4(:,:)
      integer, allocatable:: idold2new(:)
      integer, allocatable:: itemparray(:)
      real*8 wsave(1000), chwork(1000)
      real*8 xf(ndeg), yf(ndeg)
      real*8 ftemp1(ndeg,ndeg),ftemp2(ndeg,ndeg)
      real*8 coefftemp1(ndeg,ndeg),coefftemp2(ndeg,ndeg)
      real*8, allocatable:: valk1(:,:), valk2(:,:)
      real*8, allocatable:: valk3(:,:), valk4(:,:)
      real*8, allocatable:: valk1b(:,:), valk2b(:,:)
      real*8, allocatable:: valk3b(:,:), valk4b(:,:)
      
      allocate(itemparray(maxboxes))
      allocate(iempty(maxboxes))
      allocate(itemp(maxboxes))
      allocate(itemp4(4,maxboxes))
      allocate(idold2new(maxboxes))
c    demension1
      allocate(valk1(ndeg,ndeg))
      allocate(valk2(ndeg,ndeg))
      allocate(valk3(ndeg,ndeg))
      allocate(valk4(ndeg,ndeg))
c    demension2
      allocate(valk1b(ndeg,ndeg))
      allocate(valk2b(ndeg,ndeg))
      allocate(valk3b(ndeg,ndeg))
      allocate(valk4b(ndeg,ndeg))

c---------------------------------------------------------
c     1) refinement sweep. go down the tree level by level
c        refine leaf boxes if feval(:,t) is unresolved
c        no need to worry about uval, since by assumption
c        it is resolved
c---------------------------------------------------------
c
      ix=1
      iy=1
c      xlength=1.0d0
      xlength=xsize0
      x0=cent0(1)-xsize0/2.0d0
      y0=cent0(2)-xsize0/2.0d0
      nside=1
c
      call chxcin(ndeg,wsave)
      do i=0, maxlevel-1
        iflag=0
c          iflag: if the last level has been divided at all
        xlength=xlength/2.0d0
        nside=nside*2
c
        istart=istartlev(i)
        iend=istart+nblevel(i)-1

        do j=istart, iend
          ibox=iboxlev(j)
          if(ichildbox(1,ibox) .lt. 0) then

            do jj = 1, ndeg
              do ii = 1, ndeg
                 ftemp1(jj,ii) = uval1(jj+(ii-1)*8,ibox)
                 ftemp2(jj,ii) = uval2(jj+(ii-1)*8,ibox)
              enddo
            enddo
c demension1
            call getcoeff(ndeg,ndeg,ftemp1,coefftemp1,wsave)
            call geterror(ndeg,coefftemp1,error1)
c demension2
            call getcoeff(ndeg,ndeg,ftemp2,coefftemp2,wsave)
            call geterror(ndeg,coefftemp2,error2)

            hh = dble(4**levelbox(ibox))
            epscaled = eps * hh
            if(error1 .ge. epscaled .or. error2 .ge. epscaled)then
c             call subdivide
              call subdivide1(ibox,iparentbox,ichildbox,
     1             nboxes,irowbox,icolbox,levelbox,nlev,
     2             istartlev, nblevel, iboxlev,itemparray)
              
c             need to assign function values to
              ik1=ichildbox(1,ibox)
              ik2=ichildbox(2,ibox)
              ik3=ichildbox(3,ibox)
              ik4=ichildbox(4,ibox)
              call chebyp2k(ndeg,uval1(1,ibox),uval1(1,ik1),
     1             uval1(1,ik2),uval1(1,ik3),uval1(1,ik4))
              call chebyp2k(ndeg,uval2(1,ibox),uval2(1,ik1),
     1             uval2(1,ik2),uval2(1,ik3),uval2(1,ik4))
c
              iflag = 1
            endif
          endif

        enddo
c
        if((iflag .eq. 0) .and. (i .eq. nlev))then
          goto 110
        endif
      enddo

c------------------------------------------
c 2. go up the tree to coarsen:
c    if on a certain parent box, both functions
c    are resolved, merge its children
c------------------------------------------
c

110   do ii=1,nboxes
        iempty(ii)=0
      enddo
c
      do i=nlev,0,-1
c       to merge, go up the tree level by level
        istart=istartlev(i)
        iend=istart+nblevel(i)-1
        do j=istart, iend
          ibox=iboxlev(j)
          nleafkids=0
          if(ichildbox(1,ibox) .gt. 0) then
            kids(1)=ichildbox(1,ibox)
            kids(2)=ichildbox(2,ibox)
            kids(3)=ichildbox(3,ibox)
            kids(4)=ichildbox(4,ibox)
c     get nleafkids
            do k=1,4
              if(ichildbox(1,kids(k)).lt.0) then
                nleafkids=nleafkids+1
              endif
            enddo
          endif
c
          if(nleafkids .eq. 4) then
            do k2=1,ndeg
            do k1=1,ndeg
              k=(k2-1)*ndeg+k1

              valk1(k1,k2)=uval1(k,kids(1))
              valk2(k1,k2)=uval1(k,kids(2))
              valk3(k1,k2)=uval1(k,kids(3))
              valk4(k1,k2)=uval1(k,kids(4))

              valk1b(k1,k2)=uval2(k,kids(1))
              valk2b(k1,k2)=uval2(k,kids(2))
              valk3b(k1,k2)=uval2(k,kids(3))
              valk4b(k1,k2)=uval2(k,kids(4))
            enddo
            enddo

            call chebyk2p(ndeg, valk1, valk2, valk3, valk4,
     1           ftemp1, coefftemp1)
            call geterror(ndeg,coefftemp1,uerror1)

            call chebyk2p(ndeg, valk1b, valk2b, valk3b, valk4b,
     1           ftemp2, coefftemp2)
            call geterror(ndeg,coefftemp2,uerror2)

            hh = dble(4**levelbox(ibox))
            epscaled = eps * hh
c            write(*,*)uerror1,uerror2,epscaled
            if(
     1         uerror1 .lt.epscaled .and. uerror2 .lt.epscaled) then
c             the condition is correct here!
              call coarsen1(ibox,iparentbox,ichildbox,
     1             nboxes,irowbox,icolbox,levelbox,nlev,
     2             istartlev,nblevel,iboxlev,iempty)
              do k2=1,ndeg
              do k1=1,ndeg
                k=(k2-1)*ndeg+k1
                uval1(k,ibox)=ftemp1(k1,k2)
                uval2(k,ibox)=ftemp2(k1,k2)
              enddo
              enddo
c             assign interpolated function values
c             to the box that has become a leaf box
            endif
          endif
        enddo
      enddo



c------------------------------------------
c     3) clean up the intermediate tree with deleted
c        boxes, reorder uval at the same time
c------------------------------------------
c

      call reorder(iparentbox, ichildbox, nboxes,
     1     irowbox, icolbox, levelbox, nlev,
     2     istartlev, nblevel, iboxlev, iempty,
     3     idold2new, nold, itemp, itemp4)
c
      call reorderf(ndeg*ndeg, nold, uval1, idold2new)
      call reorderf(ndeg*ndeg, nold, uval2, idold2new)

      deallocate(itemparray)
      deallocate(iempty)
      deallocate(itemp)
      deallocate(itemp4)
      deallocate(idold2new)
c
      deallocate(valk1)
      deallocate(valk2)
      deallocate(valk3)
      deallocate(valk4)

      deallocate(valk1b)
      deallocate(valk2b)
      deallocate(valk3b)
      deallocate(valk4b)



      end subroutine

c --------------------------------------------------------
c subtoutine mktreef4
c
c Input:
c h1,h2:source term
c f1,f2:initial functions(independent of time)
c
      subroutine mktreef4(levelbox, icolbox, irowbox,
     1           nboxes, nlev, iparentbox, ichildbox,
     2           nblevel, iboxlev, istartlev,
     3           maxboxes, itemparray, maxlevel, eps, h1,h2,
     4           f1,f2,
     5           ndeg, cent0, xsize0,t)
c-----global variables
      implicit none
      integer *4  levelbox(1), maxboxes
      integer *4  nlev, nboxes,  maxlevel
      integer *4  icolbox(1), irowbox(1)
      integer *4  iparentbox(1), ichildbox(4,1)
      integer *4  nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer *4  itemparray(1), ndeg
      real *8  eps, h1,h2,f1,f2
      real *8 cent0(2), xsize0
c-----local variables
      integer *4  i, ibox, iflag
      integer *4  j, istart, iend, jj
      integer *4  levflag, ii

      real *8  coefftemp1(ndeg,ndeg)
      real *8  coefftemp2(ndeg,ndeg)
      real *8  coefftemp3(ndeg,ndeg)
      real *8  coefftemp4(ndeg,ndeg)

      real *8  epscaled
      real *8  error1,error2, error3,error4
      real *8  hh
      real *8  xf(ndeg), yf(ndeg)

      real *8  ftemp1(ndeg,ndeg)
      real *8  ftemp2(ndeg,ndeg)
      real *8  ftemp3(ndeg,ndeg)
      real *8  ftemp4(ndeg,ndeg)

      real *8  wsave2(1000)
      real *8  wsave3(1000)
      real *8  wsave4(1000)
      real *8  wsave5(1000)
      real *8  t

      do i = 0, maxlevel
        nblevel(i) = 0
        istartlev(i) = 0
      end do
      do i = 1, maxboxes
        iboxlev(i) = 0
      end do
c
c     first set the big parent box to the
c     appropriate settings:
c     (initially, there is just one box and
c     it is the big parent box at level 0)
c
      ibox = 1
      levelbox(ibox) = 0
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = -1
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1

      nboxes = 1
      nlev = 0

c     we also need to initialize the adaptive 'ladder'
c     structures to the correct initial values:
      nblevel(0) = 1
      istartlev(0) = 1
      iboxlev(1) = 1

      call chxcin(ndeg,wsave2)
      call chxcin(ndeg,wsave3)      
      call chxcin(ndeg,wsave4)
      call chxcin(ndeg,wsave5)
c
      do i = 0, maxlevel - 1
      iflag = 0
      istart = istartlev(i)
      iend = istart + nblevel(i) - 1
      do j = istart, iend
       ibox = iboxlev(j)

       call mkgrid(xf,yf,icolbox(ibox),irowbox(ibox),
     1      levelbox(ibox),ndeg,cent0,xsize0)

       do jj = 1, ndeg
         do ii = 1, ndeg
           ftemp1(jj,ii) = h1(xf(jj),yf(ii),t)
           ftemp2(jj,ii) = h2(xf(jj),yf(ii),t)
           ftemp3(jj,ii) = f1(xf(jj),yf(ii))
           ftemp4(jj,ii) = f2(xf(jj),yf(ii))
         end do
       end do

c      compute chebyshev transforms
       call getcoeff(ndeg,ndeg,ftemp1,coefftemp1,wsave2)
       call getcoeff(ndeg,ndeg,ftemp2,coefftemp2,wsave3)
       call getcoeff(ndeg,ndeg,ftemp3,coefftemp3,wsave4)
       call getcoeff(ndeg,ndeg,ftemp4,coefftemp4,wsave5)

       call geterror(ndeg,coefftemp1,error1)
       call geterror(ndeg,coefftemp2,error2)
       call geterror(ndeg,coefftemp3,error3)
       call geterror(ndeg,coefftemp4,error4)

       hh = dble(4**levelbox(ibox))
       epscaled = eps * hh


c       write(*,*)error1,error2,error3,error4,epscaled
       if(error1 .ge. epscaled .or. error2 .ge. epscaled
     1     .or. error3 .ge. epscaled .or.
     2      error4 .ge. epscaled )then
c        call subdivide
         call subdivide1(ibox,iparentbox,ichildbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev,itemparray)
         iflag = 1
       endif
      end do
      if(iflag .eq. 0)then
c      nothing was divided at the
c      last level, so exit the loop.
       return
      endif
      end do
      return
      end subroutine
