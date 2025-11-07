c**************************************************************
c
c     The Heat Equation:
c     heat2dper0: homogeneous heat eqn solver, nonadaptive version 
c     heat2dperf: inhomogeneous heat eqn solver, nonadaptive version
c     heat2dperfm: inhomogeneous heat eqn solver, adaptive version
c
c     Semilinear Heat Equation:
c
c     heat2d_am24: semilinear heat eqn solver, 2nd -4th order adams-moulton,
c                 nonadaptive, for now. 
c
c     RMK: eval at arbitrary user-specified targets can be done 
c          by calling tree_data_routs_nd.f/treedata_evalt_nd
c
c       (all with periodic boundary condition)
c
c     heat2d_sys_am2: reaction-diffusion system, 2nd order adams-moulton,
c                     adaptive
c     heat2d_sys_am24
c
c     modified on 2025.10.11
c     heat2d_sys_am6: reaction-diffusion system, 2nd 3rd,4th,5th and
c                    6th order adams-moulton,adaptive
c
c
c**************************************************************
c
c
c
c
c
      subroutine heat2dper0(norder, finit, dt, ntot,
     1           eps, mxltree, mxboxes, mxlevels,
     2           ltree, itree, iptr, centers, nlevels, 
     3           boxsize, nboxes, usol)
      implicit real*8 (a-h,o-z)
      integer norder, ndim 
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels) 
      real*8 usol(norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
      complex*16 zpars(10), zk
c     allocatables
c     the tree structure
      real*8, allocatable:: fvals(:,:,:)
      real*8, allocatable:: fright(:,:)
      real*8, allocatable:: potp(:)
c     for tests only
      real*8, allocatable:: targs(:,:)
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      external finit
c
c     1. create a tree to resolve the initial function
c        given by finit
      call prini(6,13)
      pi=3.1415926535897932d0
      ifnewtree=0
      iptype=2
c
      ipoly=1
      iperiod=1
c
      ifpgh=1
      ifpghtarg=0
c
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      boxlen=1.0d0
      nd=1
c
      epstree=1.0d-1*eps
      eta=1.0d0
      zk=1.0d0
c
      ntarg=1
c
c     call vol_tree_mem
c     rhsfun -> finit
      call vol_tree_mem(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1     norder,iptype,eta,finit,nd,dpars,zpars,ipars,
     2     ifnewtree,nboxes,nlevels,ltree,rintl)
      write(*,*) 'nboxes=', nboxes
      write(*,*) 'nlevels=', nlevels
      write(*,*) 'ltree=', ltree
c
c     allocate mem
      allocate(fvals(nd,npbox,nboxes))
      allocate(fright(npbox,nboxes))
c
      allocate(targs(ndim,ntarg))
      allocate(potp(ntarg))
c
c     call vol_tree_build
      call vol_tree_build(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,finit,nd,dpars,zpars,ipars,rintl,
     2    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals)
c
c--------------------------------------------------------------
c     2. the initial condition resolved
c     on the tree, function vals saved
c     in the array fvals
c     now call subroutine initpot
c     to perform the actual computation
c
      do j=1, nboxes
        do i=1, npbox
          usol(i,j)=fvals(1,i,j)
        enddo
      enddo
c
      do it=1, ntot
        tt=it*dt
        write(*,*) '*********************************'
        write(*,*) 'time step:', it
        write(*,*) '*********************************'
c
        write(*,*) 'calling initpot'
        do j=1, nboxes
          do i=1, npbox
            fright(i,j)=usol(i,j)
          enddo
        enddo
c
      ifvtarg=1
      ifptarg=0
c
      call initpot(norder, ltree, itree, iptr, centers, 
     1     nlevels, boxsize, nboxes, fright, dt,
     2     iperiod, ntarg, targs, ifvtarg, ifptarg, 
     3     eps, usol, potp)
c
      enddo

c
      end subroutine





c--------------------------------------------------
c     IO list is pretty much the same as heat2dper0
c
c     except for:
c     finit: the handle for the initial condition 
c     fforce: the handle for the forcing term
c
c--------------------------------------------------
c
      subroutine heat2dperf(norder, nordert, finit, fforce,
     1           dt, ntot, eps, mxltree, mxboxes, mxlevels,
     2           ltree, itree, iptr, centers, nlevels, 
     3           boxsize, nboxes, ifnewtree, usol)
      implicit real*8 (a-h,o-z)
      integer norder, ndim 
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels) 
      real*8 usol(norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
      complex*16 zpars(10), zk
c     allocatables
c     the tree structure
      real*8, allocatable:: fvals(:,:,:)
      real*8, allocatable:: fright(:,:)
      real*8, allocatable:: pot(:,:), vpot(:,:)
      real*8, allocatable:: potp(:), vpotp(:)
c     for tests only
      real*8, allocatable:: targs(:,:)
      real*8, allocatable:: targ(:)
      real *8, allocatable :: xref(:,:),wts(:)
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      character*1 type
      external finit, fforce
c
c     1. create a tree to resolve the initial function
c        given by finit
      call prini(6,13)
      pi=3.1415926535897932d0
      ifnewtree=0
c     ifnewtree: use the adaptive mesh or not
c     fix it to be zero temporarily
      iptype=2
c
      ipoly=1
      iperiod=1
c
      ifpgh=1
      ifpghtarg=0
c
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      boxlen=1.0d0
      nd=1
c
      epstree=1.0d-1*eps
c      epstree=1.0d-4*eps
ccc   tree refined for debugging
      eta=1.0d0
      zk=1.0d0
c
      ntarg=1
c
c     call vol_tree_mem
c     rhsfun -> finit
      call vol_tree_mem(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1     norder,iptype,eta,finit,nd,dpars,zpars,ipars,
     2     ifnewtree,nboxes,nlevels,ltree,rintl)
      write(*,*) 'nboxes=', nboxes
      write(*,*) 'nlevels=', nlevels
      write(*,*) 'ltree=', ltree
c
c     allocate mem
      allocate(fvals(nd,npbox,nboxes))
      allocate(fright(npbox,nboxes))
c     the initial potential and volume potential
      allocate(pot(npbox,mxboxes))
      allocate(vpot(npbox,mxboxes))
c
      allocate(targs(ndim,ntarg))
      allocate(potp(ntarg))
      allocate(vpotp(ntarg))
c
      allocate(targ(ndim))
      allocate(xref(ndim,npbox),wts(npbox))
c
c     call vol_tree_build
      call vol_tree_build(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,finit,nd,dpars,zpars,ipars,rintl,
     2    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals)
c
c-------------------------------------
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
      endif
      write(*,*) 'tree printed out for visualization'
c
c--------------------------------------------------------------
c     2. the initial condition resolved
c     on the tree, function vals saved
c     in the array fvals
c     now call subroutine initpot
c     to perform the actual computation
c
      itype = 0
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c
      do j=1, nboxes
        do i=1, npbox
          usol(i,j)=fvals(1,i,j)
        enddo
      enddo
c
      do it=1, ntot
        tt=it*dt
        tpre=(it-1)*dt
        write(*,*) '*********************************'
        write(*,*) 'time step:', it
        write(*,*) '*********************************'
c
        if(ifnewtree .gt. 0) then
c         do something to resolve things 
c         for the current time step
        endif
c
        write(*,*) 'calling initpot'
        do j=1, nboxes
          do i=1, npbox
            fright(i,j)=usol(i,j)
          enddo
        enddo
c
        do i=1, nboxes
        do j=1, npbox
          pot(j,i)=0.0d0
          vpot(j,i)=0.0d0
        enddo
        enddo
c
c     i. compute the initial potential
        ifvtarg=1
        ifptarg=0
        call initpot(norder, ltree, itree, iptr, centers, 
     1       nlevels, boxsize, nboxes, fright, dt,
     2       iperiod, ntarg, targs, ifvtarg, ifptarg, 
     3       eps, pot, potp)
c
c     ii. compute the volume potential
c
      write(*,*) 'calling volpot'
      call volpot(norder, nordert, ltree, itree, iptr, 
     4     centers, nlevels, boxsize, nboxes, fforce, tpre,
     2     dt, iperiod, ntarg, targs, ifvtarg, ifptarg, 
     3     eps, vpot, vpotp)
c
c     iii. add the potentials
        do i=1, nboxes
          do j=1, npbox
            usol(j,i)=pot(j,i)+vpot(j,i)
          enddo
        enddo
c
c     print out for debugging
c     pot is correct 
c     it must be vpot
c      do ilevel=0,nlevels
c        bs = boxsize(ilevel)/2.0d0
c        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
c          if(itree(iptr(4)+ibox-1).eq.0) then
c             do j=1,npbox
c               do k=1,ndim
c                 targ(k)=centers(k,ibox) + xref(k,j)*bs
c               enddo
c
c               write(34,*) targ(1),targ(2), pot(j,ibox), 
c     2                     vpot(j,ibox), usol(j,ibox)
c             enddo
c          endif
c        enddo
c      enddo
c
c
c       interpolate to an user specified grid if needed      
      enddo
c
c     at the final time step, interpolate to the 
c     user-specified grid

c
      end subroutine





c--------------------------------------------------
c     IO list is pretty much the same as heat2dper0
c
c     except for:
c     finit: the handle for the initial condition 
c     fforce: the handle for the forcing term
c
c     attention: the input list of finit/fforce is 
c                slightly different from that of heat2dperf
c
c--------------------------------------------------
c
      subroutine heat2dperfm(norder, nordert, finit, fforce,
     1           dt, ntot, eps, mxltree, mxboxes, mxlevels,
     2           ltree, itree, iptr, centers, nlevels, 
     3           boxsize, nboxes, ifnewtree, usol)
      implicit real*8 (a-h,o-z)
      integer norder, ndim 
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels) 
      real*8 usol(norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
      complex*16 zpars(10), zk
c     allocatables
c     the tree structure
      real*8, allocatable:: fvals(:,:,:)
      real*8, allocatable:: fright(:,:)
      real*8, allocatable:: pot(:,:), vpot(:,:)
      real*8, allocatable:: potp(:), vpotp(:)
c     for tests only
      real*8, allocatable:: targs(:,:)
      real*8, allocatable:: targ(:)
      real *8, allocatable :: xref(:,:),wts(:)
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      character*1 type
      external finit, fforce
c
c     1. create a tree to resolve the initial function
c        given by finit
      call prini(6,13)
      pi=3.1415926535897932d0
      ifnewtree=1
c     ifnewtree: use the adaptive mesh or not
c     fix this to be 1 momentarily in this case
c
      iptype=2
c
      ipoly=1
      iperiod=1
c
      ifpgh=1
      ifpghtarg=0
c
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      boxlen=1.0d0
      nd=1
c
      epstree=1.0d-1*eps
c      epstree=1.0d-4*eps
      eta=1.0d0
      zk=1.0d0
c
      ntarg=1
c
c     call vol_tree_mem2, the version that works better
c     for the adaptive case
c     be careful about the interface
c     not necessarily the same as vol_tree_mem

      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1     norder,iptype,eta,finit,nd,dpars,zpars,ipars,
     2     ifnewtree,mboxes,mlevels,ltree,rintl)
c
      dpars(1)=dt
      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1     norder,iptype,eta,fforce,nd,dpars,zpars,ipars,
     2     ifnewtree,mboxes,mlevels,ltree,rintl)
c
      mboxes=mboxes*2
      mlevels=mlevels+5
      ltree=ltree*2
c
c     not enough. need to increase these numbers
c     a better memory estimate?
c
c------------------------------
c     rintl initialized in mem2
c     different for different input of functions
c     leads to different trees
c     seems okay for now, 
c     different by one uniform refinement
c------------------------------
c
      write(*,*) 'mboxes=', mboxes
      write(*,*) 'mlevels=', mlevels
      write(*,*) 'ltree=', ltree
c     not sure if we need mxboxes, mxlevels, mxltree
c     as inputs
c
c     allocate mem
      allocate(fvals(nd,npbox,mboxes))
      allocate(fright(npbox,mboxes))
c     the initial potential and volume potential
      allocate(pot(npbox,mboxes))
      allocate(vpot(npbox,mboxes))
c
c     what is ntarg though
      allocate(targs(ndim,ntarg))
      allocate(potp(ntarg))
      allocate(vpotp(ntarg))
c
c     why do we need these?
      allocate(targ(ndim))
      allocate(xref(ndim,npbox),wts(npbox))
c
c     call vol_tree_build2 to resolve finit
c     again, the adaptive version
      call vol_tree_build2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,finit,nd,dpars,zpars,ipars,rintl,
     2    nboxes,mboxes,nlevels,mlevels,ltree,itree,iptr,centers,
     3    boxsize,fvals)
      write(*,*) 'nboxes, mboxes=', nboxes, mboxes
      write(*,*) 'nlevels, mlevels=', nlevels, mlevels
      write(*,*) 'ltree=', ltree
c
      do ilev=1, mlevels
        write(*,*) '---',ilev,rintl(ilev)
      enddo
c
c---------------------------------------------------------------
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
      endif
      write(*,*) 'tree printed out for visualization'
c
c--------------------------------------------------------------
c     2. the initial condition resolved
c     on the tree, function vals saved
c     in the array fvals
c     now call subroutine initpot
c     to perform the actual computation
c
      itype = 0
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c     what is this though...
c     shouldn't it be changed to another grid setting routine?
c
      do j=1, nboxes
        do i=1, npbox
          usol(i,j)=fvals(1,i,j)
        enddo
      enddo
c

      do it=1, ntot
        tt=it*dt
        tpre=(it-1)*dt
        write(*,*) '*********************************'
        write(*,*) 'time step:', it
        write(*,*) 'tt=', tt
        write(*,*) '*********************************'
c
        if(ifnewtree .gt. 0) then
c         adapt the tree so that fforce is resolved at
c         the current time step
c         maintain the function value array: usol
          write(*,*) '--------------------'
          write(*,*) 'calling vol_tree_adap_fv'
          dpars(1)=tt
          write(*,*) 'dpars(1)=',dpars(1)
          call vol_tree_adap_fv(nd,ndim,epstree,ipoly,
     1         iptype,norder,npbox,nboxes,mboxes,nlevels,
     2         mlevels,ltree,itree,iptr,centers,boxsize,
     3         rintl,dpars,zpars,ipars,iperiod,fforce,usol)
c         debug: send the testing functions of test16
c         to this interface, disable all the other things
c         to debug the seg fault
c
c         take out the rest of the code
c         debug this subroutine alone
          write(*,*) 'done calling adap_fv'
          write(*,*) 'nboxes, mboxes=', nboxes, mboxes
          write(*,*) '--------------------'

        endif
c
        do j=1, nboxes
          do i=1, npbox
            fright(i,j)=usol(i,j)
          enddo
        enddo
c
        do i=1, nboxes
        do j=1, npbox
          pot(j,i)=0.0d0
          vpot(j,i)=0.0d0
        enddo
        enddo
c
c     i. compute the initial potential
        ifvtarg=1
        ifptarg=0
c       might need to modify this interface
c       add in mboxes, mlevels?
c       I believe it's not a problem
        write(*,*) 'calling initpot'
        call initpot(norder, ltree, itree, iptr, centers, 
     1       nlevels, boxsize, nboxes, fright, dt,
     2       iperiod, ntarg, targs, ifvtarg, ifptarg, 
     3       eps, pot, potp)
c
c     ii. compute the volume potential
c       might need to modify this interface
c       add in mboxes, mlevels
        write(*,*) 'calling volpot'
        call volpot(norder, nordert, ltree, itree, iptr, 
     4       centers, nlevels, boxsize, nboxes, fforce, tpre,
     2       dt, iperiod, ntarg, targs, ifvtarg, ifptarg, 
     3       eps, vpot, vpotp)
c
c     iii. add the potentials
        do i=1, nboxes
          do j=1, npbox
            usol(j,i)=pot(j,i)+vpot(j,i)
          enddo
        enddo
c
c     print out for debugging
c     pot is correct 
c     it must be vpot
c      do ilevel=0,nlevels
c        bs = boxsize(ilevel)/2.0d0
c        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
c          if(itree(iptr(4)+ibox-1).eq.0) then
c             do j=1,npbox
c               do k=1,ndim
c                 targ(k)=centers(k,ibox) + xref(k,j)*bs
c               enddo
c
c               write(34,*) targ(1),targ(2), pot(j,ibox), 
c     2                     vpot(j,ibox), usol(j,ibox)
c             enddo
c          endif
c        enddo
c      enddo
c
c
c       interpolate to an user specified grid if needed      
      enddo
c
c     at the final time step, interpolate to the 
c     user-specified grid

c
      end subroutine

c--------------------------------------------------
c     IO list is pretty much the same as heat2dper0
c
c     except for:
c     finit: the handle for the initial condition
c     fforce: the handle for the forcing term
c
c     attention: the input list of finit/fforce is
c                slightly different from that of heat2dperf
c
c--------------------------------------------------
c
      subroutine heat2dperfm3(norder, nordert, finit,
     1           fforce,dt, ntot, eps, mxltree, mxboxes,
     2           mxlevels,ltree, itree, iptr, centers,
     3           nlevels,boxsize, nboxes, ifnewtree, usol,ntarg,
     4           nx,ny,usolp)
      implicit real*8 (a-h,o-z)
      integer norder, ndim
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes,ntarg
      integer nx,ny
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps,hx,hy
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 usol(norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
      complex*16 zpars(10), zk
c     allocatables
c     the tree structure
      real*8, allocatable:: fvals(:,:,:)
      real*8, allocatable:: fright(:,:)
      real*8, allocatable:: pot(:,:), vpot(:,:)
      real*8, allocatable:: potp(:), vpotp(:)
c     for tests only
      real*8, allocatable:: targs(:,:)
      real*8, allocatable:: targ(:)
      real *8, allocatable :: xref(:,:),wts(:)
c     user-specified points
c     solution at targets
      real*8 usolp(ntarg)
c     for extra targets
      real *8, allocatable :: coefsp(:,:,:)
      real *8, allocatable :: umat_nd(:,:,:)
      real *8, allocatable :: xq(:),umat(:,:),vmat(:,:)
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      character*1 type
      external finit, fforce
      open(77,file='target.txt')
      open(367,file='speedadap.txt')
      open(368,file='speedtotal.txt')
      open(369,file='nleaf.txt')
c
c     1. create a tree to resolve the initial function
c        given by finit
      call prini(6,13)
      pi=3.1415926535897932d0
      ifnewtree=1
c     ifnewtree: use the adaptive mesh or not
c     fix this to be 1 momentarily in this case
c
      iptype=2
c
      ipoly=1
      iperiod=1
c
      ifpgh=1
      ifpghtarg=0
c
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      boxlen=1.0d0
      nd=1
c
      epstree=1.0d-1*eps
c      epstree=1.0d-4*eps
      eta=1.0d0
      zk=1.0d0
c
      allocate(umat_nd(norder,norder,ndim))
      allocate(xq(norder),umat(norder,norder),
     1    vmat(norder,norder))
c
c      ntarg=1
c
c     call vol_tree_mem2, the version that works better
c     for the adaptive case
c     be careful about the interface
c     not necessarily the same as vol_tree_mem

      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1     norder,iptype,eta,finit,nd,dpars,zpars,ipars,
     2     ifnewtree,mboxes,mlevels,ltree,rintl)
c
      dpars(1)=dt
      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1     norder,iptype,eta,fforce,nd,dpars,zpars,ipars,
     2     ifnewtree,mboxes,mlevels,ltree,rintl)
c
      mboxes=mboxes*2
      mlevels=mlevels+5
      ltree=ltree*2
c
c     not enough. need to increase these numbers
c     a better memory estimate?
c
c------------------------------
c     rintl initialized in mem2
c     different for different input of functions
c     leads to different trees
c     seems okay for now,
c     different by one uniform refinement
c------------------------------
c
      write(*,*) 'mboxes=', mboxes
      write(*,*) 'mlevels=', mlevels
      write(*,*) 'ltree=', ltree
c     not sure if we need mxboxes, mxlevels, mxltree
c     as inputs
c
c     allocate mem
      allocate(fvals(nd,npbox,mboxes))
      allocate(fright(npbox,mboxes))
c     the initial potential and volume potential
      allocate(pot(npbox,mboxes))
      allocate(vpot(npbox,mboxes))
c
c     what is ntarg though
      allocate(targs(ndim,ntarg))
      allocate(potp(ntarg))
      allocate(vpotp(ntarg))
c
c     why do we need these?
      allocate(targ(ndim))
      allocate(xref(ndim,npbox),wts(npbox))
c     extra targets
      allocate(coefsp(nd,npbox,mboxes))
c

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
c     call vol_tree_build2 to resolve finit
c     again, the adaptive version
      call vol_tree_build2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,finit,nd,dpars,zpars,ipars,rintl,
     2    nboxes,mboxes,nlevels,mlevels,ltree,itree,iptr,centers,
     3    boxsize,fvals)
      write(*,*) 'nboxes, mboxes=', nboxes, mboxes
      write(*,*) 'nlevels, mlevels=', nlevels, mlevels
      write(*,*) 'ltree=', ltree
c
      do ilev=1, mlevels
        write(*,*) '---',ilev,rintl(ilev)
      enddo
c
c---------------------------------------------------------------
c
c      ifplot=1
c      if (ifplot.eq.1) then
c         fname1 = 'trgtree.data'
c         fname2 = 'src.data'
c         fname3 = 'targ.data'
      
c         ns=0
c         nt=0
c        given ns and nt
c        generate random sources and targets
c        and sort into the tree?
c         call print_tree_matlab(ndim,itree,ltree,nboxes,centers,
c     1       boxsize,nlevels,iptr,ns,src,nt,targs,fname1,fname2,fname3)
c      endif
c      write(*,*) 'tree printed out for visualization'

c
c--------------------------------------------------------------
c     2. the initial condition resolved
c     on the tree, function vals saved
c     in the array fvals
c     now call subroutine initpot
c     to perform the actual computation
c
      itype = 0
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c     what is this though...
c     shouldn't it be changed to another grid setting routine?
c


      do j=1, nboxes
        do i=1, npbox
          usol(i,j)=fvals(1,i,j)
        enddo
      enddo
c

      do it=1, ntot
        tt=it*dt
        tpre=(it-1)*dt
        write(*,*) '*********************************'
        write(*,*) 'time step:', it
        write(*,*) 'tt=', tt
        write(*,*) '*********************************'
        
        t3=omp_get_wtime()
c
        if(ifnewtree .gt. 0) then
c         adapt the tree so that fforce is resolved at
c         the current time step
c         maintain the function value array: usol
          write(*,*) '--------------------'
          write(*,*) 'calling vol_tree_adap_fv'
          dpars(1)=tt
          write(*,*) 'dpars(1)=',dpars(1)
          t1=omp_get_wtime()
          call vol_tree_adap_fv(nd,ndim,epstree,ipoly,
     1         iptype,norder,npbox,nboxes,mboxes,nlevels,
     2         mlevels,ltree,itree,iptr,centers,boxsize,
     3         rintl,dpars,zpars,ipars,iperiod,fforce,usol)
          t2=omp_get_wtime()
          write(367,*)t2-t1
c         debug: send the testing functions of test16
c         to this interface, disable all the other things
c         to debug the seg fault
c
c         take out the rest of the code
c         debug this subroutine alone
          write(*,*) 'done calling adap_fv'
          write(*,*) 'nboxes, mboxes=', nboxes, mboxes
          write(*,*) '--------------------'

        endif
c
        do j=1, nboxes
          do i=1, npbox
            fright(i,j)=usol(i,j)
          enddo
        enddo
c
        do i=1, nboxes
        do j=1, npbox
          pot(j,i)=0.0d0
          vpot(j,i)=0.0d0
        enddo
        enddo
c
c     i. compute the initial potential
        ifvtarg=1
        ifptarg=1
c       might need to modify this interface
c       add in mboxes, mlevels?
c       I believe it's not a problem
        write(*,*) 'calling initpot'
        call initpot(norder, ltree, itree, iptr, centers,
     1       nlevels, boxsize, nboxes, fright, dt,
     2       iperiod, ntarg, targs, ifvtarg, ifptarg,
     3       eps, pot, potp)
c
c     ii. compute the volume potential
c       might need to modify this interface
c       add in mboxes, mlevels
        write(*,*) 'calling volpot'
c        t5=omp_get_wtime()
        call volpot(norder, nordert, ltree, itree, iptr,
     4       centers, nlevels, boxsize, nboxes, fforce, tpre,
     2       dt, iperiod, ntarg, targs, ifvtarg, ifptarg,
     3       eps, vpot, vpotp)
c        t6=omp_get_wtime()
c        write(*,*)t6-t5
c
c     iii. add the potentials
        do i=1, nboxes
          do j=1, npbox
            usol(j,i)=pot(j,i)+vpot(j,i)
          enddo
        enddo
c
       do i = 1,ntarg
         usolp(i)=potp(i)+vpotp(i)
       enddo
c
c
      nleafbox = 0
      
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
      enddo
      write(*,*)'number of leaf box = ',nleafbox


      t4=omp_get_wtime()
      write(368,*)t4-t3
      write(369,*)nleafbox

c
c       interpolate to an user specified grid if needed
      enddo

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
      endif
      write(*,*) 'tree printed out for visualization'
c
c
c     at the final time step, interpolate to the
c     user-specified grid

c
      end subroutine



c--------------------------------------------------
c     the user interface for
c     2nd and 4th order semilinear heat eqn
c     which allocates memory
c     creates the tree
c     initializes the first few steps
c     and calls heat2d_am24_main
c
c     rmk:
c     ifexact: use the exact solution to initialize or not
c     iadap:
c     iadap = 0: no refinement or coarsening
c     iadap = 1: refinement and coarsening
c     iadap = 2: refinement only
c     ifnewtree: use the adaptive mesh or not
c     uexact: function handle of the exact solution, for tests only
c     nordert: 2 or 4
c--------------------------------------------------
      subroutine heat2d_am24(norder, nordert, finit, fevaln,
     1           dt, ntot, eps, mxltree, mxboxes, mxlevels,
     2           ltree, itree, iptr, centers, nlevels,
     3           boxsize, nboxes, ifnewtree, ifexact,
     4           iadap, usol)
      implicit real*8 (a-h,o-z)
      integer norder, ndim
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes
      integer mboxes,mlevels
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 usol(norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
      real*8 xtarg(ndim)
      complex*16 zpars(10), zk
c     allocatables
c     the tree structure
      real*8, allocatable:: fvals(:,:,:)
      real*8, allocatable:: fright(:,:)
      real*8, allocatable:: u0(:,:), u1(:,:), un(:,:)
      real*8, allocatable:: pot(:,:), upre(:,:)
      real*8, allocatable:: potp(:), vpotp(:)
c     for tests only
      real*8, allocatable:: targs(:,:)
      real*8, allocatable:: targ(:)
      real *8, allocatable :: xref(:,:),wts(:)
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      character*1 type
      external finit, fevaln, uexact
c
c     1. create a tree to resolve the initial function
c        given by finit
c
c      open(21,file='ifexact0.txt')
c      open(22,file='ifexact1.txt')
c      open(666,file='init_value.txt')
      call prini(6,13)
      pi=3.1415926535897932d0
      ifnewtree=0
c     ifnewtree: use the adaptive mesh or not
c     fix it to be zero temporarily
      iptype=2
      ipoly=1
      iperiod=1
c
      ifpgh=1
      ifpghtarg=0
c
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      boxlen=1.0d0
      nd=1
c
      epstree=1.0d-1*eps
c      epstree=1.0d-2*eps
c      epstree=1.0d-4*eps
ccc   tree refined for debugging
      eta=1.0d0
      zk=1.0d0
c
      ntarg=1

c     debugging
      write(*,*) norder, nordert, dt, ntot, eps
      write(*,*) mxltree, mxboxes, mxlevels
      write(*,*) nlevels, nboxes
c
c     call vol_tree_mem
c     rhsfun -> finit
c     iadapt
      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,norder,
     1    iptype,eta,finit,nd,dpars,zpars,ipars,iadap,
     2    mboxes,mlevels,ltree,rintl)


      mboxes=mboxes*10
      mlevels=mlevels+5
c     debugging
      write(*,*) norder, nordert, dt, ntot, eps
      write(*,*) mxltree, mxboxes,mboxes,mxlevels,mlevels
      write(*,*) nlevels, nboxes

c
c     allocate mem
      allocate(fvals(nd,npbox,mxboxes))
      allocate(fright(npbox,mxboxes))
c     the initial potential and volume potential
      allocate(pot(npbox,mxboxes))
      allocate(u0(npbox,mxboxes))
      allocate(u1(npbox,mxboxes))
      allocate(un(npbox,mxboxes))
c
      allocate(upre(npbox,mxboxes))
c
      allocate(targs(ndim,ntarg))
      allocate(potp(ntarg))
      allocate(vpotp(ntarg))
c
      allocate(targ(ndim))
      allocate(xref(ndim,npbox),wts(npbox))
c
c
      call vol_tree_build2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,finit,nd,dpars,zpars,ipars,rintl,nboxes,
     2    mboxes,nlevels,mlevels,ltree,itree,iptr,centers,boxsize,
     3    fvals)

c
      write(*,*) 'nboxes, mboxes=', nboxes,mboxes
      write(*,*) 'nlevels,mlevels=', nlevels,mlevels
      write(*,*) 'ltree=', ltree
c
c-------------------------------------
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
c-------------------------------------------------
c     2. the initial condition resolved
c     on the tree, function vals saved
c     in the array fvals
c     now call subroutine heat2d_am24_main
c     to perform the actual computation
c
      itype = 0
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c
      do j=1, nboxes
        do i=1, npbox
          u0(i,j)=fvals(1,i,j)
         write(666,*)j,i,u0(i,j)
        enddo
      enddo
c     assign u0 = fvals
c
c      print out usol for debugging
c      do ilevel=0,nlevels
c        bs = boxsize(ilevel)/2.0d0
c        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
c          if(itree(iptr(4)+ibox-1).eq.0) then
c            do j=1, npbox
c              write(201,*) usol(j,ibox)
c            enddo
c          endif
c        enddo
c      enddo
c
c
c-------------------------------------------------
c
c     needs allocatables: u0, u1, un
c
c     assign initial values using exact values
c     if allowed
      if(ifexact .gt. 0) then
        do ilevel=0, nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
c             a leaf node
c             get the grid and eval the forcing term (using un)
              do j=1, npbox
                do k=1, ndim
                  xtarg(k)=centers(k,ibox) + xref(k,j)*bs
                enddo
ccc
                if(nordert .eq. 4) then
                  dpars(1)=dt
                  call uexact(nd,xtarg,dpars,zpars,ipars,u1(j,ibox))
                  dpars(1)=2.0d0*dt
                  call uexact(nd,xtarg,dpars,zpars,ipars,un(j,ibox))
                elseif(nordert .eq. 3) then
                  dpars(1)=dt
                  call uexact(nd,xtarg,dpars,zpars,ipars,un(j,ibox))
                endif
c             end of the grid pt loop
              enddo
            endif
c         end of the ibox loop
          enddo
c       end of the ilevel loop
        enddo
c
      do i=1, nboxes
        do j=1, npbox
          write(22,*)i,j,u0(j,i),u1(j,i),un(j,i)
        enddo
      enddo

      else
c       ifexact = 0
c       use smaller steps
c       w/ richardson extrapolation
        write(*,*)'calling heat2d_init'
c
        call heat2d_init(norder, nordert, fevaln,
     1       dt, eps, mxltree, mxboxes, mxlevels,rintl,
     2       ltree, itree, iptr, centers, nlevels,
     3       boxsize, nboxes, ifnewtree, u0, u1, un)
c       be careful about the calling sequence
c       right now it seems okay
c
      do i=1, nboxes
        do j=1, npbox
          write(21,*)i,j,u0(j,i),u1(j,i),un(j,i)
        enddo
      enddo

      endif



c
c     call heat2d_am24_main
c      write(*,*) norder, nordert, dt, ntot, eps
c      write(*,*) nlevels, nboxes
c
      call heat2d_am24_main(norder, nordert, fevaln,
     1     dt, ntot, eps, mxltree,mxboxes,mxlevels,rintl,
     2     ltree, itree, iptr, centers, nlevels,
     3     boxsize, nboxes, ifnewtree,iadap, u0, u1, un)
c
      do i=1, nboxes
        do j=1, npbox
          usol(j,i)=un(j,i)
        enddo
      enddo
c
      end subroutine

c********************************************************
c     subroutine heat2d_init
c********************************************************
c
c     interface is the same as heat2d_am24_main,
c     except ntot is not needed
c
c     input:
c     u0: the initial condition resolved on a tree
c     ltree - nboxes: the tree
c     fevaln: the eval routine of the forcing term
c
c     output:
c     u1: the solution at t_1, only needed if nordert ==4
c     un: the solution at t_{nordert-2}, only needed if nordert >=3
c
c********************************************************
c
      subroutine heat2d_init(norder, nordert, fevaln,
     1           dt0, eps, mxltree, mxboxes,mxlevels,rintl,
     2           ltree, itree, iptr, centers, nlevels,
     3           boxsize, nboxes, ifnewtree, u0, u1, un)
      implicit real*8 (a-h,o-z)
      integer norder, ndim
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes
      integer nboxes,nlevels
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt0, eps
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
c
      real*8 u0(norder**ndim,mxboxes)
      real*8 u1(norder**ndim,mxboxes)
      real*8 un(norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
c     richardson related
      integer nr(3), ntots(3)
      real*8 rwhts(3), dts(3)
      complex*16 zpars(10), zk
c     allocatables, not sure yet
      real*8, allocatable::usol(:,:,:), targs(:,:)
      real*8, allocatable::usave(:,:,:)
      external fevaln
c      open(225,file='test_init.txt')
c      open(667,file='init_time1.txt')
c
      pi=3.1415926535897932d0
      ifnewtree=0
c
c     ifnewtree: use the adaptive mesh or not
c     fix it to be zero
c     i.e. do not mess around with the tree inside
c     the initialization routine
c
      iptype=2
      ipoly=1
      iperiod=1
c
      ifpgh=1
      ifpghtarg=0
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      boxlen=1.0d0
      nd=1
c
      epstree=1.0d-1*eps
      eta=1.0d0
      zk=1.0d0

c
      allocate(usol(npbox,mxboxes,3))
      allocate(usave(npbox,mxboxes,2))
      allocate(targs(ndim,ntarg))
c
c
c     fix the total num of steps to be 2
      ntot0=2
c
c     richardson parameters
c     do three passes of the 2nd order solve
c     set parameters for each
      nr(1)=1
      nr(2)=2
      nr(3)=4
c     number of refinements for each pass
c
      dts(1)=dt0/nr(1)
      dts(2)=dt0/nr(2)
      dts(3)=dt0/nr(3)
c     smaller delta_t values
c
      ntots(1)=ntot0*nr(1)
      ntots(2)=ntot0*nr(2)
      ntots(3)=ntot0*nr(3)
c     total num of (smaller) steps
c
      rwhts(1)=1.0d0/21.0d0
      rwhts(2)=-4.0d0/7.0d0
      rwhts(3)=32.0d0/21.0d0
c     weights for each Richardson pass
c
c     initialize usave to be zero
      do jj=1, 2
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
              do j=1, npbox
                usave(j,ibox,jj)=0.0d0
              enddo
            endif
          enddo
        enddo
      enddo
c
c------------------------------------------------------
c
c     do the three passes of Richarson solve
c     2nd order -> 3rd order -> 4th order
      do jp=1, 3
        dt=dts(jp)
        ntot=ntots(jp)
c
c       set dt and ntot, march in time using
c       the 2nd order adams-moulton
c
c       initialize each pass
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
              do j=1, npbox
                usol(j,ibox,jp)=u0(j,ibox)
              enddo
            endif
          enddo
        enddo

c
        do it=1, ntot
          tt=it*dt
          tpre=(it-1)*dt
          write(*,*) '*********************************'
          write(*,*) 'Richardson step:', it, jp
          write(*,*) '*********************************'
c
          ifvtarg = 1
          ifptarg = 0
          nordert0 = 2
c
c         remember to allocate usol
c         and save things in the right place
c          write(*,*) norder, nordert0, nboxes
c          write(*,*) tpre, dt, iperiod, ntarg
c          write(*,*) ifvtarg, ifptarg, eps

          iadap = 0
  
          call vpot_semilin_am2(norder, nordert0, ltree,
     1         itree, iptr, centers, nlevels,mxlevels, boxsize,
     2         nboxes, mxboxes,rintl,fevaln, tpre, dt, iperiod, ntarg,
     3         targs, ifvtarg, ifptarg,iadap,eps,usol(1,1,jp))
c       initialize each pass
        if(jp .eq. 1 .and. it .eq. 1)then
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
              do j=1, npbox
                write(667,*)ibox,j,usol(j,ibox,1)
              enddo
            endif
          enddo
        enddo
        endif
c
c         solution obtained for each Richardson pass
c         do the linear combination of them all
c         to get 3rd and 4th order initial values
c
c         save the solution if it is one big step
          if(mod(it,nr(jp)) .eq. 0) then
            ii=it/nr(jp)
            do ibox=1, nboxes
            do i=1, npbox
              usave(i,ibox,ii)=usave(i,ibox,ii)+rwhts(jp)*
     1                                    usol(i,ibox,jp)
            enddo
            enddo
          endif
c        test
          if(mod(it,nr(jp)) .eq. 0) then
            ii=it/nr(jp)
            do ibox=1, nboxes
            do i=1, npbox
              write(225,*)jp,it,ii,ibox,i,usol(i,ibox,jp)
            enddo
            enddo
          endif
c
c       end of the time stepping it-loop
        enddo
c
c     end of the Richardson jp-loop
      enddo
c
c     update: un, u1, copy over from usave
      if(nordert .eq. 3) then
        do ibox=1, nboxes
        do i=1, npbox
          un(i,ibox)=usave(i,ibox,1)
        enddo
        enddo
      elseif(nordert .eq. 4) then
        do ibox=1, nboxes
        do i=1, npbox
          u1(i,ibox)=usave(i,ibox,1)
          un(i,ibox)=usave(i,ibox,2)
c        write(*,*)ibox,i,usave(i,ibox,1),usave(i,ibox,1)
        enddo
        enddo
      endif
c     like this?
c     yeah, like this
c
      end subroutine

c--------------------------------------------------
c     subroutine heat2d_am24_main
c--------------------------------------------------
c
c     semilinear heat eqn solver, nonadaptive version
c     periodic boundary condition is imposed
c     2nd - 4th order adam-moulton method implemented
c
c     u_t=Lapace(u) + f(t,x,u)
c     u periodic on [-0.5,0.5]^2
c
c     interface is practically the same as others
c     except for the user-specified subroutine fevaln
c
c     fevaln has the calling sequence
c     fevaln(t,x,u,f)
c     where x is a dim 2 vector
c     secant.f/secant is called for the
c     solution of the resulting nonlinear eqn
c
c     remember to interpolate the output
c     to a given grid
c
c     input:
c     the same as heat2d_am2, except for
c     the tree is given on input
c     u0: initial cond, sampled on the 2D tensor product grid on the leaf nodes
c     u1: the solution at t_1, only needed if nordert ==4
c     un: the solution at t_{nordert-2}, only needed if nordert >=3
c
c--------------------------------------------------
c
      subroutine heat2d_am24_main(norder, nordert, fevaln,
     1        dt, ntot, eps, mxltree, mxboxes, mxlevels,rintl,
     2           ltree, itree, iptr, centers, nlevels,
     3           boxsize, nboxes, ifnewtree,iadap,u0, u1, un)
      implicit real*8 (a-h,o-z)
      integer norder, ndim
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 u0(norder**ndim,mxboxes)
      real*8 u1(norder**ndim,mxboxes)
      real*8 un(norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
      complex*16 zpars(10), zk
c     allocatables
c     the tree structure
      real*8, allocatable:: fvals(:,:,:)
      real*8, allocatable:: fright(:,:)
      real*8, allocatable:: pot(:,:), upre(:,:)
      real*8, allocatable:: potp(:), vpotp(:)
c     for tests only
      real*8, allocatable:: targs(:,:)
      real*8, allocatable:: targ(:)
      real *8, allocatable :: xref(:,:),wts(:)
      real*8, allocatable:: fnm1(:,:),fnm2(:,:)
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      character*1 type
      external fevaln,uexact
c
c     set parameters and initialize things
      call prini(6,13)
      pi=3.1415926535897932d0
      ifnewtree=0
c     ifnewtree: use the adaptive mesh or not
c     fix it to be zero temporarily
      iptype=2
c
      ipoly=1
      iperiod=1
c
      ifpgh=1
      ifpghtarg=0
c
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      boxlen=1.0d0
      nd=1
c
      epstree=1.0d-1*eps
c      epstree=1.0d-2*eps
c      epstree=1.0d-4*eps
ccc   tree refined for debugging
      eta=1.0d0
      zk=1.0d0
c

c     attention: tree is built outside
c
      allocate(fnm1(npbox,mxboxes))
      allocate(fnm2(npbox,mxboxes))
c
      allocate(targ(ndim))
      allocate(xref(ndim,npbox),wts(npbox))
      itype = 0
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c
c     initialize un, fnm1, fnm2
      if(nordert .eq. 2) then
c       2nd order:
c       assign un = u0
        do j=1, nboxes
          if(itree(iptr(4)+j-1).eq.0) then
            do i=1, npbox
              un(i,j)=u0(i,j)
            enddo
          else
            do i=1, npbox
              un(i,j)=0.0d0
            enddo
          endif
        enddo
      elseif(nordert .eq. 3) then
c       3rd order:
c       retrieve the grid on the leaves
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
               do j=1,npbox
                 do k=1,ndim
                   targ(k)=centers(k,ibox) + xref(k,j)*bs
                 enddo
c                eval f(tt=0, u=u0)-> fnm1
                 tt=0.0d0
                 call fevaln(tt,targ(1),targ(2),
     1                u0(j,ibox),fnm1(j,ibox))
               enddo
            else
              do j=1, npbox
                fnm1(j,ibox)=0.0d0
              enddo
            endif
c         end of the loop over boxes
          enddo
c       end of the loop over levels
        enddo
c
      elseif(nordert .eq. 4) then
c       4th order:
c       retrieve the grid on the leaves
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
               do j=1,npbox
                 do k=1,ndim
                   targ(k)=centers(k,ibox) + xref(k,j)*bs
                 enddo
c              eval f(tt=0, u=u0)-> fnm2
               tt=0.0d0
               call fevaln(tt,targ(1),targ(2),
     1              u0(j,ibox),fnm2(j,ibox))
c              eval f(tt=dt, u=u1)-> fnm1
               tt=dt
               call fevaln(tt,targ(1),targ(2),
     1              u1(j,ibox),fnm1(j,ibox))
               enddo
            else
              do j=1, npbox
                fnm2(j,ibox)=0.0d0
                fnm1(j,ibox)=0.0d0
              enddo
            endif
c         end of the loop over boxes
          enddo
c       end of the loop over levels
        enddo
      endif
c
c     time stepping:
      do it=nordert-1, ntot
        tt=it*dt
        tpre=(it-nordert+1)*dt
c       correct? correct.
c       adams-moulton, remember
        write(*,*) '*********************************'
        write(*,*) 'time step:', it
        write(*,*) '*********************************'
c
c       call vpot_semilin_am24
c       which contains one fgt and a pt-wise nonlinear solve
        ifvtarg = 1
        ifptarg = 0
        write(*,*) 'calling vpot---'
c
c        write(*,*) nd, ndim, delta, eps
c        write(*,*) ipoly, iperiod, norder, npbox
c        write(*,*) nboxes, nlevels, ifvtarg
c        write(*,*) ntarg, ifptarg,mxltree,mxboxes,mxlevels
c        write(*,*)nordert,ltree,tpre,dt,iadap,ifnewtree

        call vpot_semilin_am24(norder, nordert, ltree, itree,
     1           iptr, centers,mxltree, nlevels,mxlevels, boxsize,
     2         nboxes,mxboxes,fevaln, tpre, dt, iperiod, ntarg, targs,
     3     rintl,ifvtarg, ifptarg,ifnewtree,iadap, eps, fnm2, fnm1,un)

c
      sum_err=0.0d0
      sum_u=0.0d0
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               do k=1,ndim
                 targ(k)=centers(k,ibox) + xref(k,j)*bs
               enddo
c               write(*,*) j, ibox, targ(1), targ(2), tf
               dpars(1)=tt
               call uexact(nd,targ,dpars,zpars,ipars,uu)
               sum_u=sum_u+uu**2
               sum_err=sum_err+abs(un(j,ibox)-uu)**2
             enddo
          endif
        enddo
      enddo
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'err2_rel=',err2_rel


c
      enddo
c

      end subroutine

c--------------------------------------------------
c self convergence version
      subroutine heat2d_am24_self(norder, nordert, finit, fevaln,
     1           dt, ntot, eps, mxltree, mxboxes, mxlevels,
     2           ltree, itree, iptr, centers, nlevels,
     3           boxsize, nboxes, ifnewtree, ifexact,
     4           iadap, usol,ntarg,nx,ny,usolp)
      implicit real*8 (a-h,o-z)
      integer norder, ndim
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes
      integer mboxes,mlevels
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 usol(norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
      real*8 xtarg(ndim)
      integer ntarg,nx,ny
      real*8 hx,hy
      real*8 usolp(ntarg)

      complex*16 zpars(10), zk
c     allocatables
c     the tree structure
      real*8, allocatable:: fvals(:,:,:)
      real*8, allocatable:: fright(:,:)
      real*8, allocatable:: u0(:,:), u1(:,:), un(:,:)
      real*8, allocatable:: pot(:,:), upre(:,:)
      real*8, allocatable:: potp(:), vpotp(:)
c     for tests only
      real*8, allocatable:: targs(:,:)
      real*8, allocatable:: targ(:)
      real *8, allocatable :: xref(:,:),wts(:)
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      character*1 type
      external finit, fevaln, uexact
c
c     1. create a tree to resolve the initial function
c        given by finit
c
c      open(21,file='ifexact0.txt')
c      open(22,file='ifexact1.txt')
      call prini(6,13)
      pi=3.1415926535897932d0
      ifnewtree=0
c     ifnewtree: use the adaptive mesh or not
c     fix it to be zero temporarily
      iptype=2
      ipoly=1
      iperiod=1
c
      ifpgh=1
      ifpghtarg=0
c
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      boxlen=1.0d0
      nd=1
c
      epstree=1.0d-1*eps
c      epstree=1.0d-2*eps
c      epstree=1.0d-4*eps
ccc   tree refined for debugging
      eta=1.0d0
      zk=1.0d0


c     debugging
      write(*,*) norder, nordert, dt, ntot, eps
      write(*,*) mxltree, mxboxes, mxlevels
      write(*,*) nlevels, nboxes
c
c     call vol_tree_mem
c     rhsfun -> finit
c     iadapt
      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,norder,
     1    iptype,eta,finit,nd,dpars,zpars,ipars,iadap,
     2    mboxes,mlevels,ltree,rintl)

      mboxes=mboxes*10
      mlevels=mlevels+5
c     debugging
      write(*,*) norder, nordert, dt, ntot, eps
      write(*,*) mxltree, mxboxes,mboxes,mxlevels,mlevels
      write(*,*) nlevels, nboxes

c
c     allocate mem
      allocate(fvals(nd,npbox,mxboxes))
      allocate(fright(npbox,mxboxes))
c     the initial potential and volume potential
      allocate(pot(npbox,mxboxes))
      allocate(u0(npbox,mxboxes))
      allocate(u1(npbox,mxboxes))
      allocate(un(npbox,mxboxes))
c
      allocate(upre(npbox,mxboxes))
c
      allocate(targs(ndim,ntarg))
      allocate(potp(ntarg))
      allocate(vpotp(ntarg))
c
      allocate(targ(ndim))
      allocate(xref(ndim,npbox),wts(npbox))
c

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
c
      call vol_tree_build2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,finit,nd,dpars,zpars,ipars,rintl,nboxes,
     2    mboxes,nlevels,mlevels,ltree,itree,iptr,centers,boxsize,
     3    fvals)

c
      write(*,*) 'nboxes, mboxes=', nboxes,mboxes
      write(*,*) 'nlevels,mlevels=', nlevels,mlevels
      write(*,*) 'ltree=', ltree

c       do i=1,nboxes
c          do j =1,npbox
c            fvals(1,j,i)=0.0d0
c          enddo
c       enddo
c
c-------------------------------------
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
c-------------------------------------------------
c     2. the initial condition resolved
c     on the tree, function vals saved
c     in the array fvals
c     now call subroutine heat2d_am24_main
c     to perform the actual computation
c
      itype = 0
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c
      do j=1, nboxes
        do i=1, npbox
          u0(i,j)=fvals(1,i,j)
c         write(*,*)j,i,u0(i,j)
        enddo
      enddo
c     assign u0 = fvals
c
c      print out usol for debugging
c      do ilevel=0,nlevels
c        bs = boxsize(ilevel)/2.0d0
c        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
c          if(itree(iptr(4)+ibox-1).eq.0) then
c            do j=1, npbox
c              write(201,*) usol(j,ibox)
c            enddo
c          endif
c        enddo
c      enddo
c
c
c-------------------------------------------------
c
c     needs allocatables: u0, u1, un
c
c     assign initial values using exact values
c     if allowed
      if(ifexact .gt. 0) then
        do ilevel=0, nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
c             a leaf node
c             get the grid and eval the forcing term (using un)
              do j=1, npbox
                do k=1, ndim
                  xtarg(k)=centers(k,ibox) + xref(k,j)*bs
                enddo
ccc
                if(nordert .eq. 4) then
                  dpars(1)=dt
                  call uexact(nd,xtarg,dpars,zpars,ipars,u1(j,ibox))
                  dpars(1)=2.0d0*dt
                  call uexact(nd,xtarg,dpars,zpars,ipars,un(j,ibox))
                elseif(nordert .eq. 3) then
                  dpars(1)=dt
                  call uexact(nd,xtarg,dpars,zpars,ipars,un(j,ibox))
                endif
c             end of the grid pt loop
              enddo
            endif
c         end of the ibox loop
          enddo
c       end of the ilevel loop
        enddo
c
      do i=1, nboxes
        do j=1, npbox
          write(22,*)i,j,u0(j,i),u1(j,i),un(j,i)
        enddo
      enddo

      else
c       ifexact = 0
c       use smaller steps
c       w/ richardson extrapolation
        write(*,*)'calling heat2d_init'

        call heat2d_init(norder, nordert, fevaln,
     1       dt, eps, mxltree, mxboxes, mxlevels,rintl,
     2       ltree, itree, iptr, centers, nlevels,
     3       boxsize, nboxes, ifnewtree, u0, u1, un)
c       be careful about the calling sequence
c       right now it seems okay
c
      do i=1, nboxes
        do j=1, npbox
          write(21,*)i,j,u0(j,i),u1(j,i),un(j,i)
        enddo
      enddo

      endif



c
c     call heat2d_am24_main
c      write(*,*) norder, nordert, dt, ntot, eps
c      write(*,*) nlevels, nboxes
c
      call heat2d_am24_main(norder, nordert, fevaln,
     1     dt, ntot, eps, mxltree,mxboxes,mxlevels,rintl,
     2     ltree, itree, iptr, centers, nlevels,
     3     boxsize, nboxes, ifnewtree,iadap, u0, u1, un)
c
      do i=1, nboxes
        do j=1, npbox
          usol(j,i)=un(j,i)
        enddo
      enddo

c

      call evalt(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,usol,
     2       ntarg,targs,usolp)
c
      end subroutine



c--------------------------------------------------
c     subroutine heat2d_sys_am2
c--------------------------------------------------
c
c     semilinear heat eqn solver, nonadaptive version
c     periodic boundary condition is imposed
c     2nd order adam-moulton method implemented
c     adam-moulton of arbitrary order will be added soon
c
c     u_t=vis1*Lapace(u) + f_1(t,x,u,v)
c     v_t=vis2*Lapace(v) + f_2(t,x,u,v)
c     u, v: periodic on [-0.5,0.5]^2
c
c     be careful about the calling seq of fevaln
c     which needs to be consistent with that of
c     secant.f/newton
c
c     nd: the dimension of the system 
c         not to be confused with ndim, which
c         is the dimension of the variable
c
c--------------------------------------------------
c
      subroutine heat2d_sys_am2(nd, norder, nordert, finit, 
     1           fevaln, devaln, visc, dt, ntot, eps, mxltree, 
     2           mxboxes, mxlevels, ltree, itree, iptr, 
     3           centers, nlevels, boxsize, nboxes, iadap, usol)
      implicit real*8 (a-h,o-z)
      integer nd, norder, ndim 
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps, visc(nd)
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels) 
      real*8 usol(nd,norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
      complex*16 zpars(10), zk
c     allocatables
c     the tree structure
      real*8, allocatable:: fvals(:,:,:)
      real*8, allocatable:: fright(:,:,:)
      real*8, allocatable:: pot(:,:,:), upre(:,:,:)
      real*8, allocatable:: potp(:,:), vpotp(:,:)
c     for tests only
      real*8, allocatable:: targs(:,:)
      real*8, allocatable:: targ(:)
      real *8, allocatable :: xref(:,:),wts(:)
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      character*1 type
      external finit, fevaln, devaln
c
c     debugging
c      write(*,*) 'input parameters:'
c      write(*,*) nd, norder, nordert
c      write(*,*) visc(1), visc(2)
c      write(*,*) dt, ntot, eps
c      write(*,*) mxltree, mxboxes, mxlevels
c      write(*,*) nlevels, nboxes
c      write(*,*) iadap
c      write(*,*) ''
c
c     1. create a tree to resolve the initial function
c        given by finit
c
      call prini(6,13)
      pi=3.1415926535897932d0
      ifnewtree=0
c     ifnewtree: use the adaptive mesh or not
c     fix it to be zero temporarily
      iptype=2
c
      ipoly=1
      iperiod=1
c
      ifpgh=1
      ifpghtarg=0
c
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      boxlen=1.0d0
c
c      epstree=1.0d-1*eps
      epstree=eps
c      epstree=1.0d-2*eps
c      epstree=1.0d-4*eps
ccc   tree refined for debugging
      eta=1.0d0
      zk=1.0d0
c
      ntarg=1
c
c     iadap ---> ifnewtree
      ifnewtree=0
      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,norder,
     1    iptype,eta,finit,nd,dpars,zpars,ipars,ifnewtree,
     2    mboxes,mlevels,ltree,rintl)
c
      mboxes = mboxes*10
      mlevels = mlevels +10
c
      write(*,*) 'mboxes,mxboxes=',
     1            mboxes,mxboxes
      write(*,*) 'mlevels,mxlevels=', 
     1            mlevels,mxlevels
      write(*,*) 'ltree,mxltree=', ltree,mxltree
      write(*,*) ''
c
      if((mxboxes .lt. mboxes).or.
     1   (mxlevels .lt. mlevels) .or. 
     1   (mxltree .lt. mltree)) then
        write(*,*) 'not enough mem for the tree'
        write(*,*) 'please increase the mem'
        stop
      endif
c
c     allocate mem
c     be careful: all the arrays are nd-vecs
      allocate(fvals(nd,npbox,mxboxes))
      allocate(fright(nd,npbox,mxboxes))
c     the initial potential and volume potential
      allocate(pot(nd,npbox,mxboxes))
      allocate(upre(nd,npbox,mxboxes))
c
      allocate(targs(ndim,ntarg))
      allocate(potp(nd,ntarg))
      allocate(vpotp(nd,ntarg))
c
      allocate(targ(ndim))
      allocate(xref(ndim,npbox),wts(npbox))
c
      call vol_tree_build2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,finit,nd,dpars,zpars,ipars,rintl,nboxes,
     2    mboxes,nlevels,mlevels,ltree,itree,iptr,centers,boxsize,
     3    fvals)
      write(*,*) 'nboxes, nlevels=', nboxes, nlevels
c     
c     the above routines are theoretically nd-veced
c     but be careful about it anyway...
c
c-------------------------------------
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
      endif
      write(*,*) 'tree printed out for visualization'
c
c
c-------------------------------------
c     2. the initial condition resolved
c     on the tree, function vals saved
c     in the array fvals
c     now call subroutine vpot_semilin_am2
c     to perform the actual computation
c
      itype = 0
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c
c      write(*,*) nd, npbox, nboxes
c
      do j=1, nboxes
        do i=1, npbox
          do k=1, nd
            usol(k,i,j)=fvals(k,i,j)
c            write(201,*) k,i,j,fvals(k,i,j)
          enddo
        enddo
      enddo
c     assign usol
c
c
c      print out usol for debugging
c      do ilevel=0,nlevels
c        bs = boxsize(ilevel)/2.0d0
c        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
c          if(itree(iptr(4)+ibox-1).eq.0) then
c            do j=1, npbox
c              write(201,*) usol(j,ibox)
c            enddo
c          endif
c        enddo
c      enddo
c
c
c----------------------------------------------------
c
c     time stepping
      do it=nordert-1, ntot
        tt=it*dt
        tpre=(it-nordert+1)*dt
        write(*,*) '*********************************'
        write(*,*) 'time step:', it
        write(*,*) '*********************************'
c
c-----------------------------------------------
c
c        write(*,*) 'calling vpot_readiff_am2:'
c        write(*,*) 'nboxes=', nboxes
c        write(*,*) nd, npbox, nboxes, ltree
c
        ifvtarg = 1
        ifptarg = 0
c
c        write(*,*) norder, nordert, nboxes
c        write(*,*) tpre, dt, iperiod, ntarg
c        write(*,*) ifvtarg, ifptarg, eps
c
c       add iadap
c       ltree -> mxltree
        ltree = iptr(8)
        write(*,*) 'ltree=', ltree
c       why is it inconsistent with the ltree returned above?
c       something is wrong with the tree builder 
c       hope it doesn't matter
c
        call vpot_readiff_am2(nd, norder, nordert, ltree, 
     1       itree, iptr, centers, nlevels, mxlevels, boxsize, 
     2       nboxes, mxboxes, rintl, fevaln, devaln, tpre, dt, 
     3       iperiod, ntarg, targs, ifvtarg, ifptarg, visc, 
     4       iadap, eps, usol)
c
        write(*,*) ' '
        write(*,*) 'after calling vpot_readiff_am2'
        write(*,*) 'mxltree=',mxltree
        write(*,*) 'nlevels, mxlevels', nlevels, mxlevels
        write(*,*) 'nboxes, mxboxes', nboxes, mxboxes
c        write(*,*) 'tpre, dt', tpre, dt
c        write(*,*) 'iperiod, ntarg, ifvtarg, ifptarg, iadap, eps', 
c     1             iperiod, ntarg, ifvtarg, ifptarg, iadap, eps
c
      enddo
c
c
      end subroutine





c--------------------------------------------------
c     subroutine heat2d_sys_am2
c--------------------------------------------------
c
c     semilinear heat eqn solver, nonadaptive version
c     periodic boundary condition is imposed
c     2nd order adam-moulton method implemented
c     adam-moulton of arbitrary order will be added soon
c
c     u_t=vis1*Lapace(u) + f_1(t,x,u,v)
c     v_t=vis2*Lapace(v) + f_2(t,x,u,v)
c     u, v: periodic on [-0.5,0.5]^2
c
c     be careful about the calling seq of fevaln
c     which needs to be consistent with that of
c     secant.f/newton
c
c     nd: the dimension of the system
c         not to be confused with ndim, which
c         is the dimension of the variable
c
c--------------------------------------------------
c
      subroutine heat2d_sys_am2_video(nd, norder, nordert, finit,
     1           fevaln, devaln, visc, dt, ntot, eps, mxltree,
     2           mxboxes, mxlevels, ltree, itree, iptr,
     3           centers, nlevels, boxsize, nboxes, iadap, usol)
      implicit real*8 (a-h,o-z)
      integer nd, norder, ndim
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps, visc(nd)
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 usol(nd,norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
      complex*16 zpars(10), zk
c     allocatables
c     the tree structure
      real*8, allocatable:: fvals(:,:,:)
      real*8, allocatable:: fright(:,:,:)
      real*8, allocatable:: pot(:,:,:), upre(:,:,:)
      real*8, allocatable:: potp(:,:), vpotp(:,:)
c     for tests only
      real*8, allocatable:: targs(:,:)
      real*8, allocatable:: targ(:)
      real *8, allocatable :: xref(:,:),wts(:)
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      character*1 type

      real *8 xf(8),yf(8)
      external finit, fevaln, devaln
c------open output file:
      open(21,file='xf2.dat')
      open(22,file='yf2.dat')
      open(190,file='nleaves2.dat')
      open(26,file='velocity2.dat')
      open(367,file='speed2.txt')
c
c     debugging
c      write(*,*) 'input parameters:'
c      write(*,*) nd, norder, nordert
c      write(*,*) visc(1), visc(2)
c      write(*,*) dt, ntot, eps
c      write(*,*) mxltree, mxboxes, mxlevels
c      write(*,*) nlevels, nboxes
c      write(*,*) iadap
c      write(*,*) ''
c
c     1. create a tree to resolve the initial function
c        given by finit
c
      call prini(6,13)
      pi=3.1415926535897932d0
      ifnewtree=0
c     ifnewtree: use the adaptive mesh or not
c     fix it to be zero temporarily
      iptype=2
c
      ipoly=1
      iperiod=1
c
      ifpgh=1
      ifpghtarg=0
c
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      boxlen=1.0d0
c
      epstree=1.0d-1*eps
c      epstree=1.0d-2*eps
c      epstree=1.0d-4*eps
ccc   tree refined for debugging
      eta=1.0d0
      zk=1.0d0
c
      ntarg=1
c
c     iadap ---> ifnewtree
      ifnewtree=0
      t1=second()
      write(*,*)'iadap=',iadap
      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,norder,
     1    iptype,eta,finit,nd,dpars,zpars,ipars,iadap,
     2    mboxes,mlevels,ltree,rintl)
c
      mboxes = mboxes*10
      mlevels = mlevels +10
c
      write(*,*) 'mboxes,mxboxes=',
     1            mboxes,mxboxes
      write(*,*) 'mlevels,mxlevels=',
     1            mlevels,mxlevels
      write(*,*) 'ltree,mxltree=', ltree,mxltree
      write(*,*) ''
c
      if((mxboxes .lt. mboxes).or.
     1   (mxlevels .lt. mlevels) .or.
     1   (mxltree .lt. mltree)) then
        write(*,*) 'not enough mem for the tree'
        write(*,*) 'please increase the mem'
        stop
      endif
c
c     allocate mem
c     be careful: all the arrays are nd-vecs
      allocate(fvals(nd,npbox,mxboxes))
      allocate(fright(nd,npbox,mxboxes))
c     the initial potential and volume potential
      allocate(pot(nd,npbox,mxboxes))
      allocate(upre(nd,npbox,mxboxes))
c
      allocate(targs(ndim,ntarg))
      allocate(potp(nd,ntarg))
      allocate(vpotp(nd,ntarg))
c
      allocate(targ(ndim))
      allocate(xref(ndim,npbox),wts(npbox))
c
      call vol_tree_build2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,finit,nd,dpars,zpars,ipars,rintl,nboxes,
     2    mboxes,nlevels,mlevels,ltree,itree,iptr,centers,boxsize,
     3    fvals)
       t2=second()
      write(367,*)'time for building tree=',t2-t1
      write(367,*)'speed for building tree=',npbox*nboxes/(t2-t1)

      write(*,*) 'nboxes, nlevels=', nboxes, nlevels
c
c     the above routines are theoretically nd-veced
c     but be careful about it anyway...
c
c-------------------------------------
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
      endif
      write(*,*) 'tree printed out for visualization'
c
c
c-------------------------------------
c     2. the initial condition resolved
c     on the tree, function vals saved
c     in the array fvals
c     now call subroutine vpot_semilin_am2
c     to perform the actual computation
c
      itype = 0
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c
c      write(*,*) nd, npbox, nboxes
c
      do j=1, nboxes
        do i=1, npbox
          do k=1, nd
            usol(k,i,j)=fvals(k,i,j)
c            write(201,*) k,i,j,fvals(k,i,j)
          enddo
        enddo
      enddo
c     assign usol
c
c
c      print out usol for debugging
c      do ilevel=0,nlevels
c        bs = boxsize(ilevel)/2.0d0
c        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
c          if(itree(iptr(4)+ibox-1).eq.0) then
c            do j=1, npbox
c              write(201,*) usol(j,ibox)
c            enddo
c          endif
c        enddo
c      enddo
c
c
c----------------------------------------------------
c
c     time stepping
      do it=nordert-1, ntot
        tt=it*dt
        tpre=(it-nordert+1)*dt
        write(*,*) '*********************************'
        write(*,*) 'time step:', it
        write(*,*) '*********************************'
c
c-----------------------------------------------
c
c        write(*,*) 'calling vpot_readiff_am2:'
c        write(*,*) 'nboxes=', nboxes
c        write(*,*) nd, npbox, nboxes, ltree
c
        ifvtarg = 1
        ifptarg = 0
c
c        write(*,*) norder, nordert, nboxes
c        write(*,*) tpre, dt, iperiod, ntarg
c        write(*,*) ifvtarg, ifptarg, eps
c
c       add iadap
c       ltree -> mxltree
        ltree = iptr(8)
        write(*,*) 'ltree=', ltree
c       why is it inconsistent with the ltree returned above?
c       something is wrong with the tree builder
c       hope it doesn't matter
c
        t1=second()
c        call vpot_readiff_am2(nd, norder, nordert, ltree,
c     1       itree, iptr, centers, nlevels, mlevels, boxsize,
c     2       nboxes, mboxes, rintl, fevaln, devaln, tpre, dt,
c     3       iperiod, ntarg, targs, ifvtarg, ifptarg, visc,
c     4       iadap, eps, usol)
c
        call vpot_readiff_am2(nd, norder, nordert, ltree,
     1       itree, iptr, centers, nlevels, mxlevels, boxsize,
     2       nboxes, mxboxes, rintl, fevaln, devaln, tpre, dt,
     3       iperiod, ntarg, targs, ifvtarg, ifptarg, visc,
     4       iadap, eps, usol)

        t2=second()
c
        write(367,*)'time step=',it
        write(367,*)'time for one step=',t2-t1
        write(367,*)'speed for one step=',npbox*nboxes/(t2-t1)
c
        write(*,*) ' '
        write(*,*) 'after calling vpot_readiff_am2'
        write(*,*) 'mxltree=',mxltree
        write(*,*) 'nlevels, mxlevels', nlevels, mxlevels
        write(*,*) 'nboxes, mxboxes', nboxes, mxboxes
c        write(*,*) 'tpre, dt', tpre, dt
c        write(*,*) 'iperiod, ntarg, ifvtarg, ifptarg, iadap, eps',
c     1             iperiod, ntarg, ifvtarg, ifptarg, iadap, eps
c

c
      iplot=1
      if(iplot .eq. 1)then

      if(mod(it,100) .eq. 0) then

      itype = 0
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)

      write(*,*)'output data files to plot the tree',it
c
      nleafbox = 0
      
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
      enddo
      write(190,*)nleafbox

       
c
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               do k=1,ndim
                 targ(k)=centers(k,ibox) + xref(k,j)*bs
               enddo
c
              if(j .eq. 1)xf(1)=targ(1)
c
              if(j .eq. 2)xf(2)=targ(1)
c
              if(j .eq. 3)xf(3)=targ(1)
c
              if(j .eq. 4)xf(4)=targ(1)
c
              if(j .eq. 5)xf(5)=targ(1)
c
              if(j .eq. 6)xf(6)=targ(1)
c
              if(j .eq. 7)xf(7)=targ(1)
c
              if(j .eq. 8)xf(8)=targ(1)
c
              if(j .eq. 1)yf(1)=targ(2)
c
              if(j .eq. 9)yf(2)=targ(2)
c
              if(j .eq. 17)yf(3)=targ(2)
c
              if(j .eq. 25)yf(4)=targ(2)
c
              if(j .eq. 33)yf(5)=targ(2)
c
              if(j .eq. 41)yf(6)=targ(2)
c
              if(j .eq. 49)yf(7)=targ(2)
c
              if(j .eq. 57)yf(8)=targ(2)
             enddo
            do l = 1, 8
            write(21,*)xf(1),xf(2),xf(3),xf(4),xf(5),xf(6),xf(7),xf(8)
            write(22,*)yf(l),yf(l),yf(l),yf(l),yf(l),yf(l),yf(l),yf(l)
            end do
c
       do k = 1, 8
         write(26,*)usol(2,8*(k-1)+1,ibox),usol(2,8*(k-1)+2,ibox),
     1 usol(2,8*(k-1)+3,ibox),usol(2,8*(k-1)+4,ibox),
     2 usol(2,8*(k-1)+5,ibox),usol(2,8*(k-1)+6,ibox),
     3 usol(2,8*(k-1)+7,ibox),usol(2,8*(k-1)+8,ibox)
       end do
          endif
        enddo
      enddo

      endif

      endif

      enddo



c
c
      end subroutine




c--------------------------------------------------
c     subroutine heat2d_sys_am24
c--------------------------------------------------
c
c     semilinear heat eqn solver, nonadaptive version
c     periodic boundary condition is imposed
c     2nd order adam-moulton method implemented
c     adam-moulton of arbitrary order will be added soon
c
c     u_t=vis1*Lapace(u) + f_1(t,x,u,v)
c     v_t=vis2*Lapace(v) + f_2(t,x,u,v)
c     u, v: periodic on [-0.5,0.5]^2
c
c     be careful about the calling seq of fevaln
c     which needs to be consistent with that of
c     secant.f/newton
c
c     nd: the dimension of the system
c         not to be confused with ndim, which
c         is the dimension of the variable
c
c--------------------------------------------------
c
      subroutine heat2d_sys_am24(nd, norder, nordert, finit,
     1           fevaln, devaln, visc, dt, ntot, eps, mxltree,
     2           mxboxes, mxlevels, ltree, itree, iptr,
     3           centers, nlevels, boxsize, nboxes, iadap,
     4           ntarg, nx, ny, usol, itarg, usolp)
      implicit real*8 (a-h,o-z)
      integer nd, norder, ndim
      integer ntarg,nx,ny
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps, visc(nd)
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 usol(nd,norder**ndim,mxboxes)
      real*8 usolp(nd,ntarg)
      real*8 rintl(0:200), dpars(1000)
      complex*16 zpars(10), zk
c     for test
c      real*8 ustep1(nd,norder**ndim,mxboxes)
   
      real*8 hx,hy
      integer itarg
c     allocatables
c     the tree structure
      real*8, allocatable:: fvals(:,:,:)
      real*8, allocatable:: fright(:,:,:)
      real*8, allocatable:: u0(:,:,:), u1(:,:,:), un(:,:,:)
      real*8, allocatable:: pot(:,:,:), upre(:,:,:)
      real*8, allocatable:: potp(:,:), vpotp(:,:)

c     for tests only
      real*8, allocatable:: targs(:,:)
      real*8, allocatable:: targ(:)
      real *8, allocatable :: xref(:,:),wts(:)
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      character*1 type
      external finit, fevaln, devaln
c      open(901,file='init_step1.txt')
c
c     debugging
c      write(*,*) 'input parameters:'
c      write(*,*) nd, norder, nordert
c      write(*,*) visc(1), visc(2)
c      write(*,*) dt, ntot, eps
c      write(*,*) mxltree, mxboxes, mxlevels
c      write(*,*) nlevels, nboxes
c      write(*,*) iadap
c      write(*,*) ''
c
c     1. create a tree to resolve the initial function
c        given by finit
c
      call prini(6,13)
      pi=3.1415926535897932d0
      ifnewtree=0
c     ifnewtree: use the adaptive mesh or not
c     fix it to be zero temporarily
      iptype=2
c
      ipoly=1
      iperiod=1
c
      ifpgh=1
      ifpghtarg=0
c
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      boxlen=1.0d0
c
      epstree=1.0d-1*eps
c      epstree=1.0d-2*eps
c      epstree=1.0d-4*eps
ccc   tree refined for debugging
      eta=1.0d0
      zk=1.0d0

c
c     iadap ---> ifnewtree
      ifnewtree=0
      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,norder,
     1    iptype,eta,finit,nd,dpars,zpars,ipars,ifnewtree,
     2    mboxes,mlevels,ltree,rintl)
c
      mboxes = mboxes*10
      mlevels = mlevels +10
c
      write(*,*) 'mboxes,mxboxes=',
     1            mboxes,mxboxes
      write(*,*) 'mlevels,mxlevels=',
     1            mlevels,mxlevels
      write(*,*) 'ltree,mxltree=', ltree,mxltree
      write(*,*) ''
c
      if((mxboxes .lt. mboxes).or.
     1   (mxlevels .lt. mlevels) .or.
     1   (mxltree .lt. mltree)) then
        write(*,*) 'not enough mem for the tree'
        write(*,*) 'please increase the mem'
        stop
      endif
c
c     allocate mem
c     be careful: all the arrays are nd-vecs
      allocate(fvals(nd,npbox,mxboxes))
      allocate(fright(nd,npbox,mxboxes))
c     the initial potential and volume potential
      allocate(pot(nd,npbox,mxboxes))
      allocate(upre(nd,npbox,mxboxes))
c
      allocate(targs(ndim,ntarg))
      allocate(potp(nd,ntarg))
      allocate(vpotp(nd,ntarg))
c
      allocate(targ(ndim))
      allocate(xref(ndim,npbox),wts(npbox))
c
      allocate(u0(nd,npbox,mxboxes))
      allocate(u1(nd,npbox,mxboxes))
      allocate(un(nd,npbox,mxboxes))
c
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
      call vol_tree_build2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,finit,nd,dpars,zpars,ipars,rintl,nboxes,
     2    mboxes,nlevels,mlevels,ltree,itree,iptr,centers,boxsize,
     3    fvals)
      write(*,*) 'nboxes, nlevels=', nboxes, nlevels
c
c     the above routines are theoretically nd-veced
c     but be careful about it anyway...

c
c
c-------------------------------------
c     2. the initial condition resolved
c     on the tree, function vals saved
c     in the array fvals
c     now call subroutine vpot_semilin_am2
c     to perform the actual computation
c
      itype = 0
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c
c      write(*,*) nd, npbox, nboxes
c
      do j=1, nboxes
        do i=1, npbox
          do k=1, nd
            usol(k,i,j)=fvals(k,i,j)
            u0(k,i,j)=fvals(k,i,j)
          enddo
        enddo
      enddo

c       ifexact = 0
c       use smaller steps
c       w/ richardson extrapolation
        write(*,*)'calling heatsys2d_init'

        ifvtarg=1
        ifptarg=0

c        write(*,*)nd,norder,nordert,iperiod
c        write(*,*)ntarg,dt,eps,mxltree, mxboxes,mxlevels
c        write(*,*)nlevels,nboxes, ifnewtree,visc(1)
c
   
      call heat2d_sys_init(nd, norder, nordert,
     1           fevaln, devaln, iperiod, ntarg, targs,
     2           dt, eps, mxltree, mxboxes,mxlevels,rintl,
     3           ltree, itree, iptr, centers, nlevels,
     4           boxsize, nboxes, ifnewtree, visc,
     5           u0, u1, un)
       
c       test sys_init ---- pass
c        if(it .eq. 1)then
c        do ilevel=0,nlevels
c          bs = boxsize(ilevel)/2.0d0
c          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
c            if(itree(iptr(4)+ibox-1).eq.0) then
c              do j=1, npbox
c                do k=1,nd
c                write(*,*)ibox,j,u1(k,j,ibox),usol(k,j,ibox)
c                enddo
c              enddo
c            endif
c          enddo
c        enddo
c        endif
c
      call heat2d_sys_am24_main(nd,norder, nordert, fevaln,
     1     devaln, iperiod, ntarg, targs,
     2     dt, ntot, eps, mxltree,mxboxes,mxlevels,rintl,
     3     ltree, itree, iptr, centers, nlevels,
     4     boxsize, nboxes, ifnewtree,iadap,visc, u0, u1, un)
c
      do i=1, nboxes
        do j=1, npbox
          do k =1,nd
          usol(k,j,i)=un(k,j,i)
          enddo
        enddo
      enddo

c
      if(itarg .eq. 1) then
c
      call evalt(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,usol,
     2       ntarg,targs,usolp)

      endif
c
      end subroutine










c********************************************************
c     subroutine heat2d_sys_init
c********************************************************
c
c     interface is the same as heat2d_am24_main,
c     except ntot is not needed
c
c     input:
c     u0: the initial condition resolved on a tree
c     ltree - nboxes: the tree
c     fevaln: the eval routine of the forcing term
c
c     output:
c     u1: the solution at t_1, only needed if nordert ==4
c     un: the solution at t_{nordert-2}, only needed if nordert >=3
c
c********************************************************
c
      subroutine heat2d_sys_init(nd, norder, nordert,
     1           fevaln, devaln, iperiod, ntarg, targs,
     2           dt0, eps, mxltree, mxboxes,mxlevels,rintl,
     3           ltree, itree, iptr, centers, nlevels,
     4           boxsize, nboxes, ifnewtree, visc,
     5           u0, u1, un)
      implicit real*8 (a-h,o-z)
      integer norder, ndim, nd, nordert
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes,mxlevels
      integer nboxes,nlevels,ntarg
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt0, eps,visc(nd)
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
c
      real*8 u0(nd,norder**ndim,mxboxes)
      real*8 u1(nd,norder**ndim,mxboxes)
      real*8 un(nd,norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
c     richardson related
      integer nr(3), ntots(3)
      real*8 rwhts(3), dts(3)
      complex*16 zpars(10), zk

      real*8 targs(ndim,ntarg)
c     allocatables, not sure yet
      real*8, allocatable::usol(:,:,:,:)
      real*8, allocatable::usave(:,:,:,:)
      external fevaln, devaln
c      open(225,file='test_init.txt')
c      open(667,file='init_time1.txt')
c
      pi=3.1415926535897932d0
      
c
c     ifnewtree: use the adaptive mesh or not
c     fix it to be zero
c     i.e. do not mess around with the tree inside
c     the initialization routine
c
      iptype=2
      ipoly=1
c
      ifpgh=1
      ifpghtarg=0
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      boxlen=1.0d0
c
      epstree=1.0d-1*eps
      eta=1.0d0
      zk=1.0d0

     
c
      allocate(usol(nd,npbox,mxboxes,3))
      allocate(usave(nd,npbox,mxboxes,2))

c
c
c     fix the total num of steps to be 2
      ntot0=2
c
c     richardson parameters
c     do three passes of the 2nd order solve
c     set parameters for each
      nr(1)=1
      nr(2)=2
      nr(3)=4
c     number of refinements for each pass
c
      dts(1)=dt0/nr(1)
      dts(2)=dt0/nr(2)
      dts(3)=dt0/nr(3)
c     smaller delta_t values
c
      ntots(1)=ntot0*nr(1)
      ntots(2)=ntot0*nr(2)
      ntots(3)=ntot0*nr(3)
c     total num of (smaller) steps
c
      rwhts(1)=1.0d0/21.0d0
      rwhts(2)=-4.0d0/7.0d0
      rwhts(3)=32.0d0/21.0d0
c     weights for each Richardson pass
c
c     initialize usave to be zero
      do jj=1, 2
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
              do j=1, npbox
                do k=1,nd
                usave(k,j,ibox,jj)=0.0d0
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
c
c------------------------------------------------------
c
c     do the three passes of Richarson solve
c     2nd order -> 3rd order -> 4th order
      do jp=1, 3
        dt=dts(jp)
        ntot=ntots(jp)
c
c       set dt and ntot, march in time using
c       the 2nd order adams-moulton
c
c       initialize each pass
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
              do j=1, npbox
                do k=1,nd
                usol(k,j,ibox,jp)=u0(k,j,ibox)
                enddo
              enddo
            endif
          enddo
        enddo

c
        do it=1, ntot
          tt=it*dt
          tpre=(it-1)*dt
          write(*,*) '*********************************'
          write(*,*) 'Richardson step:', it, jp
          write(*,*) '*********************************'
c
          ifvtarg = 1
          ifptarg = 0
          nordert0 = 2
c
c         remember to allocate usol
c         and save things in the right place
c          write(*,*) norder, nordert0, nboxes
c          write(*,*) tpre, dt, iperiod, ntarg
c          write(*,*) ifvtarg, ifptarg, eps

          iadap = 0
c
        call vpot_readiff_am2(nd, norder, nordert0, ltree,
     1       itree, iptr, centers, nlevels, mxlevels, boxsize,
     2       nboxes, mxboxes, rintl, fevaln, devaln, tpre, dt,
     3       iperiod, ntarg, targs, ifvtarg, ifptarg, visc,
     4       iadap, eps, usol(1,1,1,jp))

c
c         solution obtained for each Richardson pass
c         do the linear combination of them all
c         to get 3rd and 4th order initial values
c
c         save the solution if it is one big step
          if(mod(it,nr(jp)) .eq. 0) then
            ii=it/nr(jp)
            do ibox=1, nboxes
            do i=1, npbox
              do k=1,nd
              usave(k,i,ibox,ii)=usave(k,i,ibox,ii)+rwhts(jp)*
     1                                    usol(k,i,ibox,jp)
              enddo
            enddo
            enddo
          endif
c
c       end of the time stepping it-loop
        enddo
c
c     end of the Richardson jp-loop
      enddo
c
c     update: un, u1, copy over from usave
      if(nordert .eq. 3) then
        do ibox=1, nboxes
        do i=1, npbox
         do k=1,nd
          un(k,i,ibox)=usave(k,i,ibox,1)
         enddo
        enddo
        enddo

      elseif(nordert .eq. 4) then
        do ibox=1, nboxes
        do i=1, npbox
          do k=1,nd
          u1(k,i,ibox)=usave(k,i,ibox,1)
          un(k,i,ibox)=usave(k,i,ibox,2)
          enddo
        enddo
        enddo
      endif


      end subroutine

c--------------------------------------------------
c
      subroutine heat2d_sys_am24_main(nd,norder, nordert,
     1        fevaln, devaln, iperiod, ntarg, targs,
     2        dt, ntot, eps, mxltree, mxboxes, mxlevels,rintl,
     3        ltree, itree, iptr, centers, nlevels,
     4        boxsize, nboxes, ifnewtree,iadap,visc,u0, u1, un)
      implicit real*8 (a-h,o-z)
      integer norder, ndim, nd
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps, visc(nd)
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 u0(nd,norder**ndim,mxboxes)
      real*8 u1(nd,norder**ndim,mxboxes)
      real*8 un(nd,norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
      complex*16 zpars(10), zk
      integer ntarg
      real*8 targs(ndim,ntarg)
c     allocatables
      real*8, allocatable:: targ(:)
      real *8, allocatable :: xref(:,:),wts(:)
      real*8, allocatable:: fnm1(:,:,:),fnm2(:,:,:)
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      character*1 type
      external fevaln,devaln
c
c     set parameters and initialize things
      call prini(6,13)
      pi=3.1415926535897932d0
      ifnewtree=0
c     ifnewtree: use the adaptive mesh or not
c     fix it to be zero temporarily
      iptype=2
c
      ipoly=1
      iperiod=1
c
      ifpgh=1
      ifpghtarg=0
c
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      boxlen=1.0d0
c
      epstree=1.0d-1*eps
c      epstree=1.0d-2*eps
c      epstree=1.0d-4*eps
ccc   tree refined for debugging
      eta=1.0d0
      zk=1.0d0
c

c     attention: tree is built outside
c
      allocate(fnm1(nd,npbox,mxboxes))
      allocate(fnm2(nd,npbox,mxboxes))
c
      allocate(targ(ndim))
      allocate(xref(ndim,npbox),wts(npbox))
      itype = 0
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c
c     initialize un, fnm1, fnm2
      if(nordert .eq. 2) then
c       2nd order:
c       assign un = u0
        do j=1, nboxes
          if(itree(iptr(4)+j-1).eq.0) then
            do i=1, npbox
              do k=1,nd
              un(k,i,j)=u0(k,i,j)
              enddo
            enddo
          else
            do i=1, npbox
              do k=1,nd
              un(k,i,j)=0.0d0
              enddo
            enddo
          endif
        enddo
      elseif(nordert .eq. 3) then
c       3rd order:
c       retrieve the grid on the leaves
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
               do j=1,npbox
                 do k=1,ndim
                   targ(k)=centers(k,ibox) + xref(k,j)*bs
                 enddo
c                eval f(tt=0, u=u0)-> fnm1
                 tt=0.0d0
c
                   dpars(1)=visc(1)
                   dpars(2)=visc(2)
                   call fevaln(nd, tt, targ, u0(1,j,ibox),
     1                  dpars, zpars, ipars, fnm1(1,j,ibox))
               enddo
            else
              do j=1, npbox
                do k=1,nd
                fnm1(k,j,ibox)=0.0d0
                enddo
              enddo
            endif
c         end of the loop over boxes
          enddo
c       end of the loop over levels
        enddo
c
      elseif(nordert .eq. 4) then
c       4th order:
c       retrieve the grid on the leaves
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
               do j=1,npbox
                 do k=1,ndim
                   targ(k)=centers(k,ibox) + xref(k,j)*bs
                 enddo
c              eval f(tt=0, u=u0)-> fnm2
               tt=0.0d0
c
                   dpars(1)=visc(1)
                   dpars(2)=visc(2)
                   call fevaln(nd, tt, targ, u0(1,j,ibox),
     1                  dpars, zpars, ipars, fnm2(1,j,ibox))
c              eval f(tt=dt, u=u1)-> fnm1
               tt=dt
c
                   dpars(1)=visc(1)
                   dpars(2)=visc(2)
                   call fevaln(nd, tt, targ, u1(1,j,ibox),
     1                  dpars, zpars, ipars, fnm1(1,j,ibox))
               enddo
            else
              do j=1, npbox
                do k=1,nd
                fnm2(k,j,ibox)=0.0d0
                fnm1(k,j,ibox)=0.0d0
                enddo
              enddo
            endif
c         end of the loop over boxes
          enddo
c       end of the loop over levels
        enddo
      endif
c
c     time stepping:
      do it=nordert-1, ntot
        tt=it*dt
        tpre=(it-nordert+1)*dt
c       correct? correct.
c       adams-moulton, remember
        write(*,*) '*********************************'
        write(*,*) 'time step:', it
        write(*,*) '*********************************'
c
c       call vpot_semilin_am24
c       which contains one fgt and a pt-wise nonlinear solve
        ifvtarg = 1
        ifptarg = 0
        write(*,*) 'calling vpot---'
c
c        write(*,*) nd, ndim, delta, eps
c        write(*,*) ipoly, iperiod, norder, npbox
c        write(*,*) nboxes, nlevels, ifvtarg
c        write(*,*) ntarg, ifptarg,mxltree,mxboxes,mxlevels
c        write(*,*)nordert,ltree,tpre,dt,iadap,ifnewtree

        call vpot_readiff_am24(nd, norder, nordert, ltree,
     1     itree, iptr, centers, nlevels, mxlevels, boxsize,
     2     nboxes, mxboxes, rintl, fevaln, devaln, tpre, dt,
     3     iperiod, ntarg, targs, ifvtarg, ifptarg, visc,
     4     iadap, eps, fnm2, fnm1, un)

c
      enddo
c

      end subroutine

c--------------------------------------------------
c     subroutine heat2d_sys_am6
c--------------------------------------------------
c
c     semilinear heat eqn solver, nonadaptive version
c     periodic boundary condition is imposed
c     2nd order adam-moulton method implemented
c     adam-moulton of arbitrary order will be added soon
c
c     u_t=vis1*Lapace(u) + f_1(t,x,u,v)
c     v_t=vis2*Lapace(v) + f_2(t,x,u,v)
c     u, v: periodic on [-0.5,0.5]^2
c
c     be careful about the calling seq of fevaln
c     which needs to be consistent with that of
c     secant.f/newton
c
c     nd: the dimension of the system
c         not to be confused with ndim, which
c         is the dimension of the variable
c
c--------------------------------------------------
c
      subroutine heat2d_sys_am6(nd, norder, nordert, finit,
     1           fevaln, devaln, visc, dt, ntot, eps, mxltree,
     2           mxboxes, mxlevels, ltree, itree, iptr,
     3           centers, nlevels, boxsize, nboxes, iadap,
     4           ntarg, nx, ny, usol, itarg, usolp)
      implicit real*8 (a-h,o-z)
      integer nd, norder, ndim
      integer ntarg,nx,ny
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps, visc(nd)
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 usol(nd,norder**ndim,mxboxes)
      real*8 usolp(nd,ntarg)
      real*8 rintl(0:200), dpars(1000)
      complex*16 zpars(10), zk
c     for test
c      real*8 ustep1(nd,norder**ndim,mxboxes)
   
      real*8 hx,hy
      integer itarg
c     allocatables
c     the tree structure
      real*8, allocatable:: fvals(:,:,:)
      real*8, allocatable:: fright(:,:,:)
      real*8, allocatable:: u0(:,:,:), u1(:,:,:), un(:,:,:)
      real*8, allocatable:: u2(:,:,:),u3(:,:,:)
      real*8, allocatable:: pot(:,:,:), upre(:,:,:)
      real*8, allocatable:: potp(:,:), vpotp(:,:)

c     for tests only
      real*8, allocatable:: targs(:,:)
      real*8, allocatable:: targ(:)
      real *8, allocatable :: xref(:,:),wts(:)
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      character*1 type
      external finit, fevaln, devaln
c      open(901,file='init_step1.txt')
c
c     debugging
c      write(*,*) 'input parameters:'
c      write(*,*) nd, norder, nordert
c      write(*,*) visc(1), visc(2)
c      write(*,*) dt, ntot, eps
c      write(*,*) mxltree, mxboxes, mxlevels
c      write(*,*) nlevels, nboxes
c      write(*,*) iadap
c      write(*,*) ''
c
c     1. create a tree to resolve the initial function
c        given by finit
c
      call prini(6,13)
      pi=3.1415926535897932d0
      ifnewtree=0
c     ifnewtree: use the adaptive mesh or not
c     fix it to be zero temporarily
      iptype=2
c
      ipoly=1
      iperiod=1
c
      ifpgh=1
      ifpghtarg=0
c
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      boxlen=1.0d0
c
      epstree=1.0d-1*eps
c      epstree=1.0d-2*eps
c      epstree=1.0d-4*eps
ccc   tree refined for debugging
      eta=1.0d0
      zk=1.0d0

c
c     iadap ---> ifnewtree
      ifnewtree=0
      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,norder,
     1    iptype,eta,finit,nd,dpars,zpars,ipars,ifnewtree,
     2    mboxes,mlevels,ltree,rintl)
c
      mboxes = mboxes*10
      mlevels = mlevels +10
c
      write(*,*) 'mboxes,mxboxes=',
     1            mboxes,mxboxes
      write(*,*) 'mlevels,mxlevels=',
     1            mlevels,mxlevels
      write(*,*) 'ltree,mxltree=', ltree,mxltree
      write(*,*) ''
c
      if((mxboxes .lt. mboxes).or.
     1   (mxlevels .lt. mlevels) .or.
     1   (mxltree .lt. mltree)) then
        write(*,*) 'not enough mem for the tree'
        write(*,*) 'please increase the mem'
        stop
      endif
c
c     allocate mem
c     be careful: all the arrays are nd-vecs
      allocate(fvals(nd,npbox,mxboxes))
      allocate(fright(nd,npbox,mxboxes))
c     the initial potential and volume potential
      allocate(pot(nd,npbox,mxboxes))
      allocate(upre(nd,npbox,mxboxes))
c
      allocate(targs(ndim,ntarg))
      allocate(potp(nd,ntarg))
      allocate(vpotp(nd,ntarg))
c
      allocate(targ(ndim))
      allocate(xref(ndim,npbox),wts(npbox))
c
      allocate(u0(nd,npbox,mxboxes))
      allocate(u1(nd,npbox,mxboxes))
      allocate(un(nd,npbox,mxboxes))
      allocate(u2(nd,npbox,mxboxes))
      allocate(u3(nd,npbox,mxboxes))
c
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
      call vol_tree_build2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,finit,nd,dpars,zpars,ipars,rintl,nboxes,
     2    mboxes,nlevels,mlevels,ltree,itree,iptr,centers,boxsize,
     3    fvals)
      write(*,*) 'nboxes, nlevels=', nboxes, nlevels
c
c     the above routines are theoretically nd-veced
c     but be careful about it anyway...

c
c
c-------------------------------------
c     2. the initial condition resolved
c     on the tree, function vals saved
c     in the array fvals
c     now call subroutine vpot_semilin_am2
c     to perform the actual computation
c
      itype = 0
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c
c      write(*,*) nd, npbox, nboxes
c
      do j=1, nboxes
        do i=1, npbox
          do k=1, nd
            usol(k,i,j)=fvals(k,i,j)
            u0(k,i,j)=fvals(k,i,j)
          enddo
        enddo
      enddo

c       ifexact = 0
c       use smaller steps
c       w/ richardson extrapolation
        write(*,*)'calling heatsys2d_init'

        ifvtarg=1
        ifptarg=0

c        write(*,*)nd,norder,nordert,iperiod
c        write(*,*)ntarg,dt,eps,mxltree, mxboxes,mxlevels
c        write(*,*)nlevels,nboxes, ifnewtree,visc(1)
c
   
      call heat2d_sys_init6(nd, norder, nordert,
     1           fevaln, devaln, iperiod, ntarg, targs,
     2           dt, eps, mxltree, mxboxes,mxlevels,rintl,
     3           ltree, itree, iptr, centers, nlevels,
     4           boxsize, nboxes, ifnewtree, visc,
     5           u0, u1, u2, u3, un)
       
c       test sys_init ---- pass
c        if(it .eq. 1)then
c        do ilevel=0,nlevels
c          bs = boxsize(ilevel)/2.0d0
c          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
c            if(itree(iptr(4)+ibox-1).eq.0) then
c              do j=1, npbox
c                do k=1,nd
c                write(*,*)ibox,j,un(k,j,ibox),usol(k,j,ibox)
c                enddo
c              enddo
c            endif
c          enddo
c        enddo
c        endif
c
      call heat2d_sys_am6_main(nd,norder, nordert, fevaln,
     1     devaln, iperiod, ntarg, targs,
     2     dt, ntot, eps, mxltree,mxboxes,mxlevels,rintl,
     3     ltree, itree, iptr, centers, nlevels,
     4     boxsize, nboxes, ifnewtree,iadap,visc, u0, u1,
     5     u2, u3, un)
c
      do i=1, nboxes
        do j=1, npbox
          do k =1,nd
          usol(k,j,i)=un(k,j,i)
          enddo
        enddo
      enddo

c
      if(itarg .eq. 1) then
c
      call evalt(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,usol,
     2       ntarg,targs,usolp)

      endif
c
      end subroutine

c********************************************************
c     subroutine heat2d_sys_init6
c********************************************************
c
c     interface is the same as heat2d_am24_main,
c     except ntot is not needed
c
c     input:
c     u0: the initial condition resolved on a tree
c     ltree - nboxes: the tree
c     fevaln: the eval routine of the forcing term
c
c     output:
c     u1: the solution at t_1, only needed if nordert ==4
c     un: the solution at t_{nordert-2}, only needed if nordert >=3
c
c********************************************************
c
      subroutine heat2d_sys_init6(nd, norder, nordert,
     1           fevaln, devaln, iperiod, ntarg, targs,
     2           dt0, eps, mxltree, mxboxes,mxlevels,rintl,
     3           ltree, itree, iptr, centers, nlevels,
     4           boxsize, nboxes, ifnewtree, visc,
     5           u0, u1, u2, u3, un)
      implicit real*8 (a-h,o-z)
      integer norder, ndim, nd, nordert
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes,mxlevels
      integer nboxes,nlevels,ntarg
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt0, eps,visc(nd)
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
c
      real*8 u0(nd,norder**ndim,mxboxes)
      real*8 u1(nd,norder**ndim,mxboxes)
      real*8 un(nd,norder**ndim,mxboxes)
      real*8 u2(nd,norder**ndim,mxboxes)
      real*8 u3(nd,norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
c     richardson related
      integer nr(5), ntots(5)
      real*8 rwhts(5), dts(5)
      complex*16 zpars(10), zk

      real*8 targs(ndim,ntarg)
c     allocatables, not sure yet
      real*8, allocatable::usol(:,:,:,:)
      real*8, allocatable::usave(:,:,:,:)
      external fevaln, devaln
c      open(225,file='test_init.txt')
c      open(667,file='init_time1.txt')
c
      pi=3.1415926535897932d0
      
c
c     ifnewtree: use the adaptive mesh or not
c     fix it to be zero
c     i.e. do not mess around with the tree inside
c     the initialization routine
c
      iptype=2
      ipoly=1
c
      ifpgh=1
      ifpghtarg=0
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      boxlen=1.0d0
c
      epstree=1.0d-1*eps
      eta=1.0d0
      zk=1.0d0

     
c
      allocate(usol(nd,npbox,mxboxes,5))
      allocate(usave(nd,npbox,mxboxes,4))

c
c
c     fix the total num of steps to be 2
      ntot0=4
c
c     richardson parameters
c     do three passes of the 2nd order solve
c     set parameters for each
      nr(1)=1
      nr(2)=2
      nr(3)=4
      nr(4)=8
      nr(5)=16
c     number of refinements for each pass
c
      dts(1)=dt0/nr(1)
      dts(2)=dt0/nr(2)
      dts(3)=dt0/nr(3)
      dts(4)=dt0/nr(4)
      dts(5)=dt0/nr(5)

c     smaller delta_t values
c
      ntots(1)=ntot0*nr(1)
      ntots(2)=ntot0*nr(2)
      ntots(3)=ntot0*nr(3)
      ntots(4)=ntot0*nr(4)
      ntots(5)=ntot0*nr(5)

c     total num of (smaller) steps
c
      rwhts(1)=1.0d0/9765.0d0
      rwhts(2)=-60.0d0/9765.0d0
      rwhts(3)=1120.0d0/9765.0d0
      rwhts(4)=-7680.0d0/9765.0d0
      rwhts(5)=16384.0d0/9765.0d0
c     weights for each Richardson pass
c
c     initialize usave to be zero
      do jj=1, 4
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
              do j=1, npbox
                do k=1,nd
                usave(k,j,ibox,jj)=0.0d0
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
c
c------------------------------------------------------
c
c     do the three passes of Richarson solve
c     2nd order -> 3rd order -> 4th order -> 5th order -> 6th order
      do jp=1, 5
        dt=dts(jp)
        ntot=ntots(jp)
c
c       set dt and ntot, march in time using
c       the 2nd order adams-moulton
c
c       initialize each pass
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
              do j=1, npbox
                do k=1,nd
                usol(k,j,ibox,jp)=u0(k,j,ibox)
                enddo
              enddo
            endif
          enddo
        enddo

c
        do it=1, ntot
          tt=it*dt
          tpre=(it-1)*dt
          write(*,*) '*********************************'
          write(*,*) 'Richardson step:', it, jp
          write(*,*) '*********************************'
c
          ifvtarg = 1
          ifptarg = 0
          nordert0 = 2
c
c         remember to allocate usol
c         and save things in the right place
c          write(*,*) norder, nordert0, nboxes
c          write(*,*) tpre, dt, iperiod, ntarg
c          write(*,*) ifvtarg, ifptarg, eps

          iadap = 0
c
        call vpot_readiff_am2(nd, norder, nordert0, ltree,
     1       itree, iptr, centers, nlevels, mxlevels, boxsize,
     2       nboxes, mxboxes, rintl, fevaln, devaln, tpre, dt,
     3       iperiod, ntarg, targs, ifvtarg, ifptarg, visc,
     4       iadap, eps, usol(1,1,1,jp))

c
c         solution obtained for each Richardson pass
c         do the linear combination of them all
c         to get 3rd and 4th order initial values
c
c         save the solution if it is one big step
          if(mod(it,nr(jp)) .eq. 0) then
            ii=it/nr(jp)
            do ibox=1, nboxes
            do i=1, npbox
              do k=1,nd
              usave(k,i,ibox,ii)=usave(k,i,ibox,ii)+rwhts(jp)*
     1                                    usol(k,i,ibox,jp)
              enddo
            enddo
            enddo
          endif
c
c       end of the time stepping it-loop
        enddo
c
c     end of the Richardson jp-loop
      enddo
c
c     update: un, u1, copy over from usave
      if(nordert .eq. 3) then
        do ibox=1, nboxes
        do i=1, npbox
         do k=1,nd
          un(k,i,ibox)=usave(k,i,ibox,1)
         enddo
        enddo
        enddo
c
      elseif(nordert .eq. 4) then
        do ibox=1, nboxes
        do i=1, npbox
          do k=1,nd
          u1(k,i,ibox)=usave(k,i,ibox,1)
          un(k,i,ibox)=usave(k,i,ibox,2)
          enddo
        enddo
        enddo
c
      elseif(nordert .eq. 5) then
        do ibox=1, nboxes
        do i=1, npbox
          do k=1,nd
          u1(k,i,ibox)=usave(k,i,ibox,1)
          u2(k,i,ibox)=usave(k,i,ibox,2)
          un(k,i,ibox)=usave(k,i,ibox,3)
          enddo
        enddo
        enddo
c
      elseif(nordert .eq. 6) then
        do ibox=1, nboxes
        do i=1, npbox
          do k=1,nd
          u1(k,i,ibox)=usave(k,i,ibox,1)
          u2(k,i,ibox)=usave(k,i,ibox,2)
          u3(k,i,ibox)=usave(k,i,ibox,3)
          un(k,i,ibox)=usave(k,i,ibox,4)
          enddo
        enddo
        enddo
      endif


      end subroutine


c--------------------------------------------------
c
      subroutine heat2d_sys_am6_main(nd,norder, nordert,
     1        fevaln, devaln, iperiod, ntarg, targs,
     2        dt, ntot, eps, mxltree, mxboxes, mxlevels,rintl,
     3        ltree, itree, iptr, centers, nlevels,
     4        boxsize, nboxes, ifnewtree,iadap,visc,u0, u1,
     5        u2, u3, un)
      implicit real*8 (a-h,o-z)
      integer norder, ndim, nd
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps, visc(nd)
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 u0(nd,norder**ndim,mxboxes)
      real*8 u1(nd,norder**ndim,mxboxes)
      real*8 u2(nd,norder**ndim,mxboxes)
      real*8 u3(nd,norder**ndim,mxboxes)
      real*8 un(nd,norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
      complex*16 zpars(10), zk
      integer ntarg
      real*8 targs(ndim,ntarg)
c     allocatables
      real*8, allocatable:: targ(:)
      real *8, allocatable :: xref(:,:),wts(:)
      real*8, allocatable:: fnm1(:,:,:),fnm2(:,:,:)
      real*8, allocatable:: fnm3(:,:,:),fnm4(:,:,:)
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      character*1 type
      external fevaln,devaln
c
c     set parameters and initialize things
      call prini(6,13)
      pi=3.1415926535897932d0
      ifnewtree=0
c     ifnewtree: use the adaptive mesh or not
c     fix it to be zero temporarily
      iptype=2
c
      ipoly=1
      iperiod=1
c
      ifpgh=1
      ifpghtarg=0
c
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      boxlen=1.0d0
c
      epstree=1.0d-1*eps
c      epstree=1.0d-2*eps
c      epstree=1.0d-4*eps
ccc   tree refined for debugging
      eta=1.0d0
      zk=1.0d0
c

c     attention: tree is built outside
c
      allocate(fnm1(nd,npbox,mxboxes))
      allocate(fnm2(nd,npbox,mxboxes))
      allocate(fnm3(nd,npbox,mxboxes))
      allocate(fnm4(nd,npbox,mxboxes))
c
      allocate(targ(ndim))
      allocate(xref(ndim,npbox),wts(npbox))
      itype = 0
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c
c     initialize un, fnm1, fnm2
      if(nordert .eq. 2) then
c       2nd order:
c       assign un = u0
        do j=1, nboxes
          if(itree(iptr(4)+j-1).eq.0) then
            do i=1, npbox
              do k=1,nd
              un(k,i,j)=u0(k,i,j)
              enddo
            enddo
          else
            do i=1, npbox
              do k=1,nd
              un(k,i,j)=0.0d0
              enddo
            enddo
          endif
        enddo
      elseif(nordert .eq. 3) then
c       3rd order:
c       retrieve the grid on the leaves
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
               do j=1,npbox
                 do k=1,ndim
                   targ(k)=centers(k,ibox) + xref(k,j)*bs
                 enddo
c                eval f(tt=0, u=u0)-> fnm1
                 tt=0.0d0
c
                   dpars(1)=visc(1)
                   dpars(2)=visc(2)
                   call fevaln(nd, tt, targ, u0(1,j,ibox),
     1                  dpars, zpars, ipars, fnm1(1,j,ibox))
               enddo
            else
              do j=1, npbox
                do k=1,nd
                fnm1(k,j,ibox)=0.0d0
                enddo
              enddo
            endif
c         end of the loop over boxes
          enddo
c       end of the loop over levels
        enddo
c
      elseif(nordert .eq. 4) then
c       4th order:
c       retrieve the grid on the leaves
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
               do j=1,npbox
                 do k=1,ndim
                   targ(k)=centers(k,ibox) + xref(k,j)*bs
                 enddo
c              eval f(tt=0, u=u0)-> fnm2
               tt=0.0d0
c
                   dpars(1)=visc(1)
                   dpars(2)=visc(2)
                   call fevaln(nd, tt, targ, u0(1,j,ibox),
     1                  dpars, zpars, ipars, fnm2(1,j,ibox))
c              eval f(tt=dt, u=u1)-> fnm1
               tt=dt
c
                   dpars(1)=visc(1)
                   dpars(2)=visc(2)
                   call fevaln(nd, tt, targ, u1(1,j,ibox),
     1                  dpars, zpars, ipars, fnm1(1,j,ibox))
               enddo
            else
              do j=1, npbox
                do k=1,nd
                fnm2(k,j,ibox)=0.0d0
                fnm1(k,j,ibox)=0.0d0
                enddo
              enddo
            endif
c         end of the loop over boxes
          enddo
c       end of the loop over levels
        enddo
c
      elseif(nordert .eq. 5) then
c       4th order:
c       retrieve the grid on the leaves
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
               do j=1,npbox
                 do k=1,ndim
                   targ(k)=centers(k,ibox) + xref(k,j)*bs
                 enddo
c              eval f(tt=0, u=u0)-> fnm3
               tt=0.0d0
c
                   dpars(1)=visc(1)
                   dpars(2)=visc(2)
                   call fevaln(nd, tt, targ, u0(1,j,ibox),
     1                  dpars, zpars, ipars, fnm3(1,j,ibox))
c              eval f(tt=dt, u=u1)-> fnm2
               tt=dt
c
                   dpars(1)=visc(1)
                   dpars(2)=visc(2)
                   call fevaln(nd, tt, targ, u1(1,j,ibox),
     1                  dpars, zpars, ipars, fnm2(1,j,ibox))
c              eval f(tt=dt, u=u2)-> fnm1
               tt=2.0d0*dt
c
                   dpars(1)=visc(1)
                   dpars(2)=visc(2)
                   call fevaln(nd, tt, targ, u2(1,j,ibox),
     1                  dpars, zpars, ipars, fnm1(1,j,ibox))
               enddo
            else
              do j=1, npbox
                do k=1,nd
                fnm3(k,j,ibox)=0.0d0
                fnm2(k,j,ibox)=0.0d0
                fnm1(k,j,ibox)=0.0d0
                enddo
              enddo
            endif
c         end of the loop over boxes
          enddo
c       end of the loop over levels
        enddo
c
      elseif(nordert .eq. 6) then
c       4th order:
c       retrieve the grid on the leaves
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
               do j=1,npbox
                 do k=1,ndim
                   targ(k)=centers(k,ibox) + xref(k,j)*bs
                 enddo
c              eval f(tt=0, u=u0)-> fnm4
               tt=0.0d0
c
                   dpars(1)=visc(1)
                   dpars(2)=visc(2)
                   call fevaln(nd, tt, targ, u0(1,j,ibox),
     1                  dpars, zpars, ipars, fnm4(1,j,ibox))
c              eval f(tt=dt, u=u1)-> fnm3
               tt=dt
c
                   dpars(1)=visc(1)
                   dpars(2)=visc(2)
                   call fevaln(nd, tt, targ, u1(1,j,ibox),
     1                  dpars, zpars, ipars, fnm3(1,j,ibox))
c              eval f(tt=dt, u=u2)-> fnm2
               tt=2.0d0*dt
c
                   dpars(1)=visc(1)
                   dpars(2)=visc(2)
                   call fevaln(nd, tt, targ, u2(1,j,ibox),
     1                  dpars, zpars, ipars, fnm2(1,j,ibox))
c              eval f(tt=dt, u=u3)-> fnm1
               tt=3.0d0*dt
c
                   dpars(1)=visc(1)
                   dpars(2)=visc(2)
                   call fevaln(nd, tt, targ, u3(1,j,ibox),
     1                  dpars, zpars, ipars, fnm1(1,j,ibox))
               enddo
            else
              do j=1, npbox
                do k=1,nd
                fnm4(k,j,ibox)=0.0d0
                fnm3(k,j,ibox)=0.0d0
                fnm2(k,j,ibox)=0.0d0
                fnm1(k,j,ibox)=0.0d0
                enddo
              enddo
            endif
c         end of the loop over boxes
          enddo
c       end of the loop over levels
        enddo
      endif
c
c     time stepping:
      do it=nordert-1, ntot
        tt=it*dt
        tpre=(it-nordert+1)*dt
c       correct? correct.
c       adams-moulton, remember
        write(*,*) '*********************************'
        write(*,*) 'time step:', it
        write(*,*) '*********************************'
c
c       call vpot_semilin_am24
c       which contains one fgt and a pt-wise nonlinear solve
        ifvtarg = 1
        ifptarg = 0
        write(*,*) 'calling vpot---'
c
c        write(*,*) nd, ndim, delta, eps
c        write(*,*) ipoly, iperiod, norder, npbox
c        write(*,*) nboxes, nlevels, ifvtarg
c        write(*,*) ntarg, ifptarg,mxltree,mxboxes,mxlevels
c        write(*,*)nordert,ltree,tpre,dt,iadap,ifnewtree

        call vpot_readiff_am6(nd, norder, nordert, ltree,
     1     itree, iptr, centers, nlevels, mxlevels, boxsize,
     2     nboxes, mxboxes, rintl, fevaln, devaln, tpre, dt,
     3     iperiod, ntarg, targs, ifvtarg, ifptarg, visc,
     4     iadap, eps, fnm4, fnm3, fnm2, fnm1, un)

c
      enddo
c

      end subroutine


