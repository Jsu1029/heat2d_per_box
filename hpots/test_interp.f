c     a testing driver
c     which tests the interpolation related subroutines
c     get_c2p_interp_matrices
c     get_p2c_interp_matrices
c     get_val2coefs_matrices
c     get_coefs2val_matrices
c
c     along with
c     ortho_trans_nd and ortho_eval_nd
c
c------------------------------------------------------
c
      program test
      implicit real *8 (a-h,o-z)
      integer ipars(100), ltree
      integer iptr(12), icoll
      integer, allocatable:: isgn(:,:)
      real*8 done, pi, rintl(0:200)
      real*8 dpars(1000), eta
      complex*16 zpars(10)
      complex*16 zk, ztmp, ima, zz
      integer, allocatable :: itree(:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8, allocatable :: fvals1(:,:,:), fvals2(:,:)
      real *8, allocatable :: fcoefs(:,:,:)
      real *8, allocatable :: fvalsc(:,:,:), fvalsp(:,:,:)
c     interp related variables
      real*8, allocatable:: polyv(:,:,:,:)  , polyvc(:,:,:,:)
      real*8, allocatable:: umat_nd(:,:,:)  
      real*8, allocatable:: grid(:,:), xy(:)
c     cheby grid related
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      character *14 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      external rhsfun1, rhsfun2, rhsfun3
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
c     what is n2
      n2=norder/2
      nd=1
      ndim=2
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
      eta=1.0d0
      zk=1.0d0
c
      ifnewtree=1
c
      call cpu_time(t1)
c     attention: this subroutine initializes rintl too
c     this version is better for adaptivity
      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,rhsfun1,nd,dpars,zpars,ipars,ifnewtree,
     2    mboxes,mlevels,ltree,rintl)
      write(*,*) 'mboxes =', mboxes
      write(*,*) 'mlevels=', mlevels 
      write(*,*) 'ltree=', ltree 
      call cpu_time(t2)
c
c
      allocate(fvals(nd,npbox,mboxes),centers(ndim,mboxes))
      allocate(boxsize(0:mlevels),itree(ltree))
      allocate(fvals1(nd,npbox,mboxes))
      allocate(fvalsp(nd,npbox,mboxes))
      allocate(fvalsc(nd,npbox,mboxes))
      allocate(fcoefs(nd,npbox,mboxes))
      allocate(fvals2(nd,npbox/2**ndim))
c
      call vol_tree_build2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,rhsfun1,nd,dpars,zpars,ipars,rintl,
     2    nboxes,mboxes,nlevels,mlevels,ltree,itree,iptr,centers,
     3    boxsize,fvals)
c
      write(*,*) ' '
      write(*,*) '******'
      write(*,*) 'after calling vol_tree_build'
      write(*,*) 'nboxes =', nboxes
      write(*,*) 'mboxes =', mboxes
c
      write(*,*) 'nlevels=', nlevels
      write(*,*) 'mlevels=', mlevels
c
      write(*,*) 'ltree=', ltree 
c
      write(*,*) ''
      write(*,*) 'laddr='
      do i=1, 2*nlevels+1, 2
        write(*,*) itree(i), itree(i+1)
      enddo
      write(*,*) '******'
c
c-----------------------------------------------
c
c     print out the tree for visualization
c
      ifplot = 1
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
c-----------------------------------------------
c
      ibox=84
      ipar=itree(ibox+iptr(3)-1) 
      ic1=itree(iptr(5)+mc*(ipar-1)+0)
      ic2=itree(iptr(5)+mc*(ipar-1)+1)
      ic3=itree(iptr(5)+mc*(ipar-1)+2)
      ic4=itree(iptr(5)+mc*(ipar-1)+3)
      write(*,*) 'ibox=', ibox
      write(*,*) 'ipar=', ipar
      write(*,*) 'ic1, ic2, ic3, ic4=', 
     1            ic1, ic2, ic3, ic4
c
      allocate(isgn(ndim,2**ndim))
      allocate(polyv(norder,n2,ndim,mc))
      allocate(polyvc(norder,norder,ndim,mc))
      allocate(umat_nd(norder,norder,ndim))
c
c     precomp the interp matrices
c     make sure you understand what they are and how to use them
      call get_child_box_sign(ndim,isgn)
      call get_c2p_interp_matrices(ndim,ipoly,norder,isgn,polyv)
      call get_p2c_interp_matrices(ndim,ipoly,norder,isgn,polyvc)
      call get_val2coefs_matrices(ndim,ipoly,norder,umat_nd)
c
c     sample the function value on the box ibox
c     and interp to the parent box ipar
c
c     grid: the classical grid
c     xy: the xy coord of the current point
      allocate(grid(ndim,npbox))
      allocate(xy(ndim))
c     local cheby grid
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
c
c-----------------------------------------------
c     initialize the grid
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
c
c-----------------------------------------------
c     sample the function on the leaf box-> fvals1
c     compare with fvals
c
c     sample the function on ipar
      ilev = itree(iptr(2)+ipar-1) 
      do i=1,npbox
        do j=1,ndim
          xy(j) = grid(j,i)*boxsize(ilev)+centers(j,ipar)
c          write(*,*) i, j, xy(j)
        enddo
        call rhsfun1(nd,xy,dpars,zpars,ipars,fvals1(1,i,ipar))
c        write(*,*) i, ipar, fvals(1,i,ipar), fvals1(1,i,ipar), 
c     1      abs(fvals(1,i,ipar)-fvals1(1,i,ipar))
      enddo
c
c
c-----------------------------------------------
c     convert all the function vals fvals
c     to cheby coeffs: fcoefs
c     by calling tree_data_trans_nd
c     ref: tree_ext.f/vol_tree_adap_fun
      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,fvals,fcoefs,umat_nd)
c
c     on each individual box: ortho_trans_nd
c
c      do j=1, npbox
c        write(*,*) fcoefs(1,j,ic1), fcoefs(1,j,ic2),
c     1             fcoefs(1,j,ic3), fcoefs(1,j,ic4)
c      enddo
c
c-----------------------------------------------
c     now do the interp to the parent
c     ref: coarsen
c
c      subroutine ortho_eval_nd(ndim,nd,n1,fin,n2,fout,umat)
c                  call ortho_eval_nd(ndim,nd,norder,fcoefs(1,1,jbox),
c     1                n2,fvals2,polyv(1,1,1,j))
c                  call tree_data_c2p_copy_nd(ndim,nd,norder,fvals2,
c     1                fvals(1,1,ibox),j)
c
c     fcoefs -> fvals2(nd,npbox/2**ndim)
c     eval cheby coeffs of each child by multiplying
c     polyv
      do j=1, mc
        jbox=itree(iptr(5)+mc*(ipar-1)+j-1)
        write(*,*) 'jbox=', jbox
        call ortho_eval_nd(ndim,nd,norder,fcoefs(1,1,jbox),
     1       n2,fvals2,polyv(1,1,1,j))
c       copy fvals2 -> fvalsp (in the right place)
        call tree_data_c2p_copy_nd(ndim,nd,norder,fvals2,
     1       fvalsp(1,1,ipar),j)
      enddo
c
c     now print out and check
      do i=1,npbox
c        write(*,*) i, ipar, fvals(1,i,ipar), fvalsp(1,i,ipar), 
c     1      abs(fvals(1,i,ipar)-fvalsp(1,i,ipar))
      enddo
c     looking good, proceed on
c
c-----------------------------------------------
c     now interp from parent ipar to child ibox
c
      call ortho_trans_nd(ndim,nd,itype,norder,
     1     fvals(1,1,ipar),fcoefs(1,1,ipar),umat_nd)
c
c     use p2c matrix polyvc to eval the chey exp
c     at children boxes
c
      call ortho_eval_nd(ndim,nd,norder,fcoefs(1,1,ipar),
     1     norder,fvalsc(1,1,ic1),polyvc(1,1,1,1))

      call ortho_eval_nd(ndim,nd,norder,fcoefs(1,1,ipar),
     1     norder,fvalsc(1,1,ic2),polyvc(1,1,1,2))

      call ortho_eval_nd(ndim,nd,norder,fcoefs(1,1,ipar),
     1     norder,fvalsc(1,1,ic3),polyvc(1,1,1,3))

      call ortho_eval_nd(ndim,nd,norder,fcoefs(1,1,ipar),
     1     norder,fvalsc(1,1,ic4),polyvc(1,1,1,4))
c
c     now print out and check the error
      do i=1,npbox
        write(*,*) i, ic1, fvals(1,i,ic1), fvalsc(1,i,ic1), 
     1      abs(fvals(1,i,ic1)-fvalsc(1,i,ic1))
      enddo
c
      do i=1,npbox
        write(*,*) i, ic2, fvals(1,i,ic2), fvalsc(1,i,ic2), 
     1      abs(fvals(1,i,ic2)-fvalsc(1,i,ic2))
      enddo
c
      do i=1,npbox
        write(*,*) i, ic3, fvals(1,i,ic3), fvalsc(1,i,ic3), 
     1      abs(fvals(1,i,ic3)-fvalsc(1,i,ic3))
      enddo
c
      do i=1,npbox
        write(*,*) i, ic4, fvals(1,i,ic4), fvalsc(1,i,ic4), 
     1      abs(fvals(1,i,ic4)-fvalsc(1,i,ic4))
      enddo
c
      end program




c------------------------------------------------------
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
      f(1)=sin(2.0d0*pi*xyz(1))*cos(2.0d0*pi*xyz(2))
c
      end subroutine



