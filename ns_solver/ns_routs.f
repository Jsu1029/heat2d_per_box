c     This code evaluates and hessian at the tensor
c     grid on an adaptive tree given the orthogonal polynomial expansion coefficients
c     at each leaf box.
c
c     input:
c     ndim - dimension of the underlying space
c     nd - integer,   number of functions
c     ipoly - 0: Legendre polynomials
c                 1: Chebyshev polynomials
c     nlevels - integer
c            number of levels
c     nboxes - integer
c            number of boxes
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
c     norder - integer
c           order of expansions for input coefficients array
c     fvals  - function values on the tree
c
c     itype - itype = i>0 - transformation is only along ith coordinates,
c                         for example, calculating the expansion coefficients of the x-derivative
c             itype = 0 - transformation is done on both directions, for example,
c                         val2coefs or coef2vals
c
c     output:
c     grad - double precision (nd,ndim,norder**ndim,nboxes)
c            gradient values on tensor grid on each leaf box
c     hess - double precision (nd,nhess,norder**ndim,nboxes)
c            hessian values on tensor grid on each leaf box
c            in the order of uxx, uxy, uyy in 2d
c            and uxx,uyy,uzz,uxy,uxz,uyz in 3d. nhess=ndim*(ndim+1)/2
      subroutine evalgh(ndim,nd,ipoly,nlevels,ltree,itree,iptr,
     1       boxsize,norder,nboxes,fvals,grad,hess)
c
      implicit real *8 (a-h,o-z)
      integer nd,ndim,ltree
      integer norder,ipoly
      integer nlevels,nboxes
      integer itree(ltree),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 fvals(nd,norder**ndim,nboxes)
  
      real *8 umat_nd(norder,norder,ndim)

      real *8 coefs(nd,norder**ndim,nboxes)
      real *8 grad(nd,ndim,norder**ndim,nboxes)
      real *8 hess(nd,ndim*(ndim+1)/2,norder**ndim,nboxes)

      real *8 umat(norder,norder)
      real *8 vmat(norder,norder)
      real *8 wts(norder),xq(norder)

      
c constructs Legendre nodes, and  corresponding Gaussian
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)

c
      norder2=norder*norder
      do i=1,ndim
         call dcopy_f77(norder2,umat,1,umat_nd(1,1,i),1)
      enddo

      itype=0
      
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1       iptr,boxsize,norder,fvals,coefs,umat_nd)
      call treedata_evalgh_nd(ndim,nd,ipoly,nlevels,itree,iptr,
     1       boxsize,norder,coefs,grad,hess)
     
      end subroutine


c     This code evaluates function values at nt target points where
c     the function is given
c
      subroutine evalt(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,pot,
     2       ntarg,targs,pote)
c
      implicit real *8 (a-h,o-z)
      integer nd,ndim,ltree
      integer norder,ipoly
      integer ntarg
      integer nlevels,nboxes
      integer itree(ltree),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 centers(ndim,nboxes)
      real *8 pot(nd,norder**ndim,nboxes)
      real *8 coefsp(nd,norder**ndim,nboxes)
      real *8 umat_nd(norder,norder,ndim)
      real *8 targs(ndim,ntarg)
      real *8 pote(nd,ntarg)

      real *8 umat(norder,norder)
      real *8 vmat(norder,norder)
      real *8 wts(norder),xq(norder)

c constructs Legendre nodes, and  corresponding Gaussian
      itype = 2
      if (ipoly.eq.0) call legeexps(itype,norder,xq,umat,vmat,wts)
      if (ipoly.eq.1) call chebexps(itype,norder,xq,umat,vmat,wts)


      norder2=norder*norder
      do i=1,ndim
         call dcopy_f77(norder2,umat,1,umat_nd(1,1,i),1)
      enddo

      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,pot,coefsp,umat_nd)
c
      call treedata_evalt_nd(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,coefsp,
     2       ntarg,targs,pote)
 
      end subroutine


c     the following subroutine defines the array that represents the
c     right hand side of the equation.

      subroutine evalf(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,
     2       dpars,zpars,ipars,fforce,fright)
c
      implicit real *8 (a-h,o-z)
      integer nd,ndim,ltree
      integer norder
      integer nlevels,nboxes,itype
      integer itree(ltree),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 centers(ndim,nboxes)
      real *8 u(nd)
      real *8 fright(nd,norder**ndim,nboxes)
      integer ipoly
c
      integer ipars(100)
      real *8 dpars(1000)
      complex*16 zpars(10)
      character*1 type
c
      real *8 umat(norder,norder)
      real *8 vmat(norder,norder)
c
      real *8 xref(ndim,norder**ndim),wts(norder**ndim)
      real *8 targ(ndim)

      external fforce
  
      npbox=norder**ndim
c
      itype = 0
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               do k=1,ndim
                 targ(k)=centers(k,ibox) + xref(k,j)*bs
               enddo
               call fforce(nd,targ,dpars,zpars,ipars,u)
               fright(1,j,ibox)=u(1)
               fright(2,j,ibox)=u(2)
             enddo
          endif
        enddo
      enddo
  

      end subroutine
