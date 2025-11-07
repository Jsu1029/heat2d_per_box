c fmm2d(Frank version)
c convert old tree to new tree
c periodic boundary condition
c
c
c     This code computes the volume potential on a box for densities
c     defined on a tensor product grid of each leaf node in an adaptive tree.
c
c     input
c     nd - integer
c          number of right hand sides
c     ndim - integer
c           dimension of the underlying space
c     eps - double precision
c           precision requested
c     ikernel - integer
c            0: the Yukawa kernel; 1: the Laplace kernel; 2: the square-root Laplace kernel.
c     beta - double precision
c            either the parameter in the Yukawa kernel or the exponent of the power
c            function kernel
c     ipoly - integer
c            0: Legendre polynomials
c            1: Chebyshev polynomials
c     norder - integer
c           order of expansions for input function value array
c     npbox - integer
c           number of points per box where potential is to be dumped = (norder**ndim)
c     fvals - double precision (nd,npbox,nboxes)
c           function values tabulated on a tensor grid in each leaf node
c     ifpgh   : flag for computing pot/grad/hess
c                   ifpgh = 1, only potential is computed
c                   ifpgh = 2, potential and gradient are computed
c                   ifpgh = 3, potential, gradient, and hessian
c                   are computed
c     ifpghtarg: flag for computing pottarg/gradtarg/hesstarg
c                    ifpghtarg = 1, only potential is computed at targets
c                    ifpghtarg = 2, potential and gradient are
c                    computed at targets
c                    ifpghtarg = 3, potential, gradient, and hessian are
c                    computed at targets
c     nboxes - integer
c            number of boxes
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
c     output:
c     pot - double precision (nd,npbox,nboxes)
c            volume potential on the tree structure (note that
c            the potential is non-zero only in the leaf boxes of the new tree
c     grad - double precision (nd,ndim,npbox,nboxes)
c            gradient of the volume potential on the tree structure
c     hess - double precision (nd,ndim*(ndim+1)/2,npbox,nboxes)
c            hessian of the volume potential on the tree structure
c            in 2d, the order is xx, xy, yy
c            in 3d, the order is xx, yy, zz, xy, xz, yz
c     pote - double precision (nd,ntarg)
c            volume potential at targets
c     grade - double precision (nd,ndim,ntarg)
c            gradient of the volume potential at targets
c     hesse - double precision (nd,ndim*(ndim+1)/2,ntarg)
c            hessian of the volume potential at targets
c
      subroutine bdmk_old_per(nd,ndim,eps,ikernel,beta,ipoly,norder,
     1    npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    fvals,ifpgh,pot,grad,hess,ntarg,targs,
     3    ifpghtarg,pote,grade,hesse,tottimeinfo)

      implicit real*8 (a-h,o-z)
      real *8 eps,beta
      integer nd,ndim
      integer ikernel,nboxes,nlevels,ntarg,ifpgh,ifpghtarg
      integer iptr(8),ltree
      integer itree(ltree),norder,npbox
      real *8 targs(ndim,ntarg)
      real *8 fvals(nd,npbox,nboxes)
      integer ipoly
      
      real *8 pot(nd,npbox,nboxes)
      real *8 grad(nd,ndim,npbox,nboxes)
      real *8 hess(nd,ndim*(ndim+1)/2,npbox,nboxes)

      real *8 pote(nd,ntarg)
      real *8 grade(nd,ndim,nboxes)
      real *8 hesse(nd,ndim*(ndim+1)/2,nboxes)

      real *8 centers(ndim,nboxes)
      real *8 boxsize(0:nlevels)
      real *8 timeinfo(20,-100:20)
c old tree
      integer   levelbox(nboxes)
      integer   nlev
      integer   icolbox(nboxes), irowbox(nboxes)
      integer   iparentbox(nboxes), ichildbox(4,nboxes)
      integer   nblevel(0:nlevels), iboxlev(nboxes)
      integer   istartlev(0:nlevels)
      real *8   cent0(2),xsize0
c    old fmm version related
      integer iperiod
      integer ndeg
      integer iprec
      real *8 fright1(npbox,nboxes)
      real *8 fright2(npbox,nboxes)
      real *8 pot1(npbox,nboxes),pot2(npbox,nboxes)
      real *8 ttotal1,ttotal2

      ndeg = norder
c
      if (eps .lt. 1.0d-8)then

          iprec=2
      else 
          iprec=3
      end if

c 1. newtree to oldtree
c      t1=omp_get_wtime()
      call newtree2oldtree2d(nlev,levelbox,iparentbox,
     1    ichildbox,icolbox,irowbox,nboxes,nblevel,
     2    iboxlev,istartlev,cent0,xsize0,iperiod,
     3    ltree,nlevels,itree,iptr,centers,boxsize)
c      t2=omp_get_wtime()
c      write(*,*)'time for newtree2oldtree2d',t2-t1

        do i=1, nboxes
          do j=1, npbox
            fright1(j,i)=fvals(1,j,i)
            fright2(j,i)=fvals(2,j,i)
c            write(23,*)i,fright1(j,i),fright2(j,i)
            pot1(j,i)=0.0d0
            pot2(j,i)=0.0d0
          enddo
        enddo
        iperiod=1

c 3.   call adapfmm2d8

       t1=omp_get_wtime()
      call adaptfmm2d8(ndeg, nlev, levelbox, iparentbox,
     1     ichildbox, icolbox, irowbox, nboxes, nblevel,
     2     iboxlev, istartlev, fright1, iprec, iperiod,
     3     pot1, ttotal1)
       t2=omp_get_wtime()
       write(*,*)'time for dmk=',t2-t1
       write(*,*)'speed for dmk=',npbox*nboxes/(t2-t1)
c
      call adaptfmm2d8(ndeg, nlev, levelbox, iparentbox,
     1     ichildbox, icolbox, irowbox, nboxes, nblevel,
     2     iboxlev, istartlev, fright2, iprec, iperiod,
     3     pot2, ttotal2)


c 4. output on new tree
c
      do i = 0,nlev
         do j = istartlev(i),istartlev(i)+nblevel(i)-1
           ibox = iboxlev(j)
           if(ichildbox(1,ibox) .lt. 0)then
              do l = 1,64
            pot(1,l,ibox)=pot1(l,ibox)
            pot(2,l,ibox)=pot2(l,ibox)
              end do
           end if
         end do
       end do

c  ifpgh
c       t1=omp_get_wtime()
      if (ifpgh .eq. 3)then
c
      call evalgh(ndim,nd,ipoly,nlevels,ltree,itree,iptr,
     1       boxsize,norder,nboxes,pot,grad,hess)

      endif
c ifpghtarg
      if (ifpghtarg .eq. 1)then
      call evalt(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,pot,
     2       ntarg,targs,pote)
      endif
c       t2=omp_get_wtime()
c       write(*,*)'time for eval=',t2-t1

      end subroutine

