c--------------------------------------------------
c
      subroutine unsteady_stokes2dperfm4(nd,norder,nordert, finit,
     1           fforce,dt, ntot, eps, mxltree, mxboxes,
     2           mxlevels,ltree, itree, iptr, centers,
     3           nlevels,boxsize, nboxes, ifnewtree,velocity,vorticity,
     4           ntarg,velocityp,visc,funs)
      implicit real*8 (a-h,o-z)
      integer norder, ndim
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes,ntarg
      integer nordert
      integer nx,ny
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps
      real*8 visc
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 velocity(nd,norder**ndim,mxboxes)
      real*8 vorticity(norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
      complex*16 zpars(10), zk
c     allocatables
c     the tree structure
      real*8, allocatable:: fvals(:,:,:)
      real*8, allocatable:: fright(:,:,:)
      real*8, allocatable:: pot(:,:,:), vpot(:,:,:)
      real*8, allocatable:: potp(:,:), vpotp(:,:)
c     for tests only
      real*8, allocatable:: targs(:,:)
      real*8, allocatable:: targ(:)
      real *8, allocatable :: xref(:,:),wts(:)
c     user-specified points
c     solution at targets
      real*8 velocityp(nd,ntarg)
c     for extra targets
      real *8, allocatable :: coefsp(:,:,:)
      real *8, allocatable :: umat_nd(:,:,:)
      real *8, allocatable :: xq(:),umat(:,:),vmat(:,:)
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      character*1 type
c initialization
c
      real *8,allocatable :: M1(:,:,:)
      real *8,allocatable :: M2(:,:,:)
      real *8,allocatable :: M3(:,:,:)
      real *8,allocatable :: M4(:,:,:)
      real *8,allocatable :: usolp(:,:)
c
      real *8,allocatable :: usol1(:,:,:)
      real *8,allocatable :: usol2(:,:,:)
      real *8,allocatable :: usol3(:,:,:)
      real *8,allocatable :: usol4(:,:,:)

      real *8 dt_init
      integer ntot_init
      real *8,allocatable :: grad(:,:,:,:)
      real *8,allocatable :: hess(:,:,:,:)
      real *8,allocatable :: velocity_ex(:,:,:)

      external finit, fforce,uexact,funs
c      open(77,file='target.txt')
c      open(78,file='usol.txt')
c      open(19,file='velocity.txt')
c      open(555,file='initfun_value.txt')
c       open(232,file='vpot.txt')
       open(999,file='testtttt.txt')
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
c
      write(*,*)'eps =',eps
c      epstree=1.0d-1*eps
      epstree=eps*1000
c      epstree=eps
      write(*,*)'epstree =',epstree
c      epstree=1.0d-4*eps
      eta=1.0d0
      zk=1.0d0
c      zk=60.0d0
c
      allocate(umat_nd(norder,norder,ndim))
      allocate(xq(norder),umat(norder,norder),
     1    vmat(norder,norder))
c
c     call vol_tree_mem2, the version that works better
c     for the adaptive case
c
      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1     norder,iptype,eta,finit,nd,dpars,zpars,ipars,
     2     ifnewtree,mboxes,mlevels,ltree,rintl)
c
      dpars(1)=dt
      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1     norder,iptype,eta,fforce,nd,dpars,zpars,ipars,
     2     ifnewtree,mboxes,mlevels,ltree,rintl)
      
      mboxes=mboxes*10
      mlevels=mlevels+5
      ltree=ltree*2

      write(*,*) 'mboxes=', mboxes
      write(*,*) 'mlevels=', mlevels
      write(*,*) 'ltree=', ltree
c     not sure if we need mxboxes, mxlevels, mxltree
c     as inputs
c
c     allocate mem
      allocate(fvals(nd,npbox,mxboxes))
      allocate(fright(nd,npbox,mxboxes))
c     the initial potential and volume potential
      allocate(pot(nd,npbox,mxboxes))
      allocate(vpot(nd,npbox,mxboxes))
c
c     what is ntarg though
      allocate(targs(ndim,ntarg))
      allocate(potp(nd,ntarg))
      allocate(vpotp(nd,ntarg))
c
c     why do we need these?
      allocate(targ(ndim))
      allocate(xref(ndim,npbox),wts(npbox))
c     extra targets
      allocate(coefsp(nd,npbox,mxboxes))
c
      allocate(M1(nd,norder**ndim,mxboxes))
c
      allocate(M2(nd,norder**ndim,mxboxes))
c
      allocate(M3(nd,norder**ndim,mxboxes))
c
      allocate(M4(nd,norder**ndim,mxboxes))
      allocate(usolp(nd,ntarg))
c
      allocate(usol1(nd,norder**ndim,mxboxes))
c
      allocate(usol2(nd,norder**ndim,mxboxes))
c
      allocate(usol3(nd,norder**ndim,mxboxes))
c
      allocate(usol4(nd,norder**ndim,mxboxes))
      
      allocate(grad(nd,ndim,npbox,mxboxes))
      allocate(hess(nd,nhess,npbox,mxboxes))
      allocate(velocity_ex(nd,npbox,mxboxes))

c
      nx=1000
      ny=1000
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
c     call vol_tree_build2 to resolve finit
c     again, the adaptive version
      write(*,*)'resolve initial function'
      call vol_tree_build2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,finit,nd,dpars,zpars,ipars,rintl,
     2    nboxes,mboxes,nlevels,mlevels,ltree,itree,iptr,centers,
     3    boxsize,fvals)
      write(*,*) 'nboxes, mboxes=', nboxes, mboxes
      write(*,*) 'nlevels, mlevels=', nlevels, mlevels
      write(*,*) 'ltree=', ltree
c
c      do ilevel=0,mlevels
c         rintl(ilevel)=0.70710678118654757
c      enddo
c
      do j=1, nboxes
        do i=1, npbox
          velocity(1,i,j)=fvals(1,i,j)
          velocity(2,i,j)=fvals(2,i,j)
c         write(*,*)fvals(1,i,j),fvals(2,i,j)
        enddo
      enddo



c
       do it = 1,ntot

        tt=it*dt
        tpre=(it-1)*dt

        write(*,*) '*********************************'
        write(*,*) 'initial time step:', it
        write(*,*) 'time=', tt
        write(*,*) '*********************************'
c
          t1=omp_get_wtime()
          write(*,*) '--------------------'
          dpars(1)=tt
          
          call vol_tree_adap_fv(nd,ndim,epstree,ipoly,
     1         iptype,norder,npbox,nboxes,mboxes,nlevels,
     2         mlevels,ltree,itree,iptr,centers,boxsize,
     3         rintl,dpars,zpars,ipars,iperiod,fforce,velocity)
          t3=omp_get_wtime()
          write(*,*) 'nboxes, mboxes=', nboxes, mboxes
c
          nleafbox = 0
      
          do i=1,nboxes
            if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
          enddo
          write(*,*)'number of leaf box = ',nleafbox

          write(*,*) '--------------------'
c
        write(*,*)'time for adap=',t3-t1
        write(*,*)'speed for adap=',npbox*nleafbox/(t3-t1)
c
        do j=1, nboxes
          do i=1, npbox
            fright(1,i,j)=velocity(1,i,j)
            fright(2,i,j)=velocity(2,i,j)
          enddo
        enddo

c
        do i=1, nboxes
        do j=1, npbox
          pot(1,j,i)=0.0d0
          pot(2,j,i)=0.0d0
          vpot(1,j,i)=0.0d0
          vpot(2,j,i)=0.0d0
        enddo
        enddo

c    compute volume potential
        write(*,*)'compute volume potential'
        ifvtarg=1
        ifptarg=0
c
c        write(*,*)norder,nordert,ltree
c        write(*,*)nlevels,mxlevels, nboxes,mxboxes
c        write(*,*)tpre,dt,ifvtarg,ifptarg,eps

        call volpot(norder,nordert, ltree, itree, iptr,
     1       centers, nlevels,mxlevels, boxsize, nboxes,mxboxes,
     2       fforce, tpre,dt, iperiod, ntarg, targs, ifvtarg,
     3       ifptarg,eps,vpot,vpotp,visc,funs)


c    compute initial potential
        write(*,*)'compute initial potential'
      call compute_initpot(norder,ltree, itree, iptr, centers,
     1           nlevels,mxlevels, boxsize, nboxes,mxboxes, fright,
     2           dt,iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           eps, pot, potp,visc)

c
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
                velocity(1,j,ibox)=vpot(1,j,ibox) + pot(1,j,ibox)
                velocity(2,j,ibox)=vpot(2,j,ibox) + pot(2,j,ibox)
c                write(232,*)j,ibox,vpot(1,j,ibox),vpot(2,j,ibox),
c     1          pot(1,j,ibox),pot(2,j,ibox)
             enddo
           endif
        enddo
      enddo
        t2=omp_get_wtime()
c
        write(*,*)'time for one step=',t2-t1
        write(*,*)'speed for one step=',npbox*nleafbox/(t2-t1)
c    exact solution
      ipoly=1
      dpars(1)=tt
      write(*,*)'final time =', dpars(1)
      call evalf(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,
     2       dpars,zpars,ipars,uexact,velocity_ex)

c
      sum_err=0.0d0
      sum_u=0.0d0
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               sum_u=sum_u+velocity_ex(1,j,ibox)**2
               sum_err=sum_err+abs(velocity(1,j,ibox)-
     1                   velocity_ex(1,j,ibox))**2
             enddo
          endif
        enddo
      enddo
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'usol err2_rel dimension1_pre=',err2_rel


c   end of it
       enddo
c   output vorticity
      call evalgh(ndim,nd,ipoly,nlevels,ltree,itree,iptr,
     1       boxsize,norder,nboxes,velocity,grad,hess)
    

c   output velocity at target points
      call evalt(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,velocity,
     2       ntarg,targs,velocityp)

c
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
              vorticity(j,ibox)=grad(2,1,j,ibox)-
     1              grad(1,2,j,ibox)
              write(999,*)ibox,j,grad(1,1,j,ibox),grad(1,2,j,ibox),
     1    grad(2,1,j,ibox), grad(2,2,j,ibox),velocity(1,j,ibox),
     2     velocity(2,j,ibox)
             enddo
           endif
        enddo
      enddo


c
      end subroutine

c**************************************************************
c     subroutine compute_initpot
c**************************************************************
c     This subroutine computes the initial potential
c
c  INPUT:
c  norder: degree of poly approx on each leaf node
c  ltree - nboxes: the tree structure
c  fright: the initial condition, resolved on the leaf nodes
c          of the tree, by norder-th order poly
c  dt: time step
c  iperiod: periodic BC on the unit box (or not)
c           iperiod = 0 => free space BC on the unit box
c           iperiod = 1 => periodic BC on the unit box
c  ntarg: number of extra targets
c  targs: coords of extra targets
c  ifvtarg: consider volume target or not
c  ifptarg: consider point target or not
c
c  OUTPUT:
c  pot: the initial potential on the volume grid
c  potp: the initial potential at extra targets
c
c  RMK: initpot = \int_{unit box} G(x-y,dt)f(y) dy
c       (here we consider the volume integral on the whole
c        unit box)
c
c**************************************************************
c
      subroutine compute_initpot(norder, ltree, itree, iptr, centers,
     1  nlevels,mxlevels, boxsize, nboxes,mxboxes, fright, dt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           eps, pot, potp,visc)
      implicit real*8 (a-h,o-z)
c     global vars
      parameter(nd=2)
      integer norder, ltree, nboxes, ndim, nlevels
      integer mxlevels,mxboxes
      integer iperiod, ntarg
      integer ifvtarg, ifptarg
      parameter(ndim=2)
      integer itree(ltree), iptr(12)
      real*8 centers(ndim, mxboxes), boxsize(0:mxlevels)
      real*8 fright(nd,norder**ndim,mxboxes), dt
      real*8 targs(ndim,ntarg), eps
      real*8 pot(nd,norder**ndim,mxboxes), potp(nd,ntarg)
      real*8 visc
c     local vars
      real*8, allocatable:: timeinfo(:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
      real*8, allocatable:: gradp(:,:,:)
      real*8, allocatable:: hessp(:,:,:)
c      open(31,file='fvals.txt')
c
      pi=3.1415926535897932d0
      nhess = ndim*(ndim+1)/2
      npbox=norder**ndim
c
      allocate(timeinfo(100))
c
      allocate(grad(nd,ndim,npbox,mxboxes))
      allocate(hess(nd,nhess,npbox,mxboxes))
      allocate(gradp(nd,ndim,ntarg))
      allocate(hessp(nd,nhess,ntarg))
c
c     allocate certain vars
c
c     simply call boxfgt and divide the
c     output by dt
c
      delta=4.0d0*dt*visc
      ifnewtree=0
      ipoly=1
c
      if(ifvtarg .gt. 0) then
        ifpgh=1
      else
        ifpgh=0
      endif
c
      if(ifptarg .gt. 0) then
        ifpghtarg=1
      else
        ifpghtarg=0
      endif
c
      nleafbox = 0
      
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
      enddo
c
c---------------------------------------
c
      do j=1, nboxes
        do i=1, npbox
          pot(1,i,j)=0.0d0
          pot(2,i,j)=0.0d0
        enddo
      enddo
   
c
      t1=omp_get_wtime()

      call boxfgt(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fright,
     2    ifpgh,pot,grad,hess,ifnewtree,ntarg,targs,
     3    ifpghtarg,potp,gradp,hessp,timeinfo)
c
       t2=omp_get_wtime()
       write(*,*)'time for fgt=',t2-t1
       write(*,*)'speed for fgt=',npbox*nleafbox/(t2-t1)

c
c---------------------------------------
c
c     divide the results by (4.0d0*pi*dt)
      if(ifvtarg .gt. 0) then
      do ilevel=0,nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               pot(1,j,ibox)=pot(1,j,ibox)/(pi*delta)
               pot(2,j,ibox)=pot(2,j,ibox)/(pi*delta)
             enddo
          endif
        enddo
      enddo
      endif
c
c----------------------------------------
c
c     divide the results by (4.0d0*pi*dt)
c
      if(ifptarg .gt. 0) then
        do i=1, ntarg
           potp(1,i)=potp(1,i)/(pi*delta)
           potp(2,i)=potp(2,i)/(pi*delta)
        enddo
      endif
c
c
      deallocate(timeinfo)
      deallocate(grad)
      deallocate(hess)
      deallocate(gradp)
      deallocate(hessp)
c
c
      end subroutine


c**************************************************************
c
      subroutine compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           initvalue,eps,frightS,frightG)
      implicit real*8 (a-h,o-z)
c     global vars
      parameter(nd=2)
      integer norder, ltree, nboxes, ndim, nlevels
      integer iperiod, ntarg
      integer ifvtarg, ifptarg
      parameter(ndim=2)
      real*8 tt
      integer itree(ltree), iptr(12)
      real*8 centers(ndim, nboxes), boxsize(0:nlevels)
      real*8 initvalue(nd,norder**ndim,nboxes)
      real*8 targs(ndim,ntarg), eps, dt
      real*8 targ(ndim)
      real *8 frightS(nd,norder**ndim,nboxes)
      real *8 frightG(nd,norder**ndim,nboxes)
c     local vars
      real *8, allocatable:: fright(:,:,:)
      real *8, allocatable:: convective(:,:,:)
c      real *8, allocatable:: frightG(:,:,:)
      integer, allocatable:: ipars(:)
      real*8, allocatable:: tottimeinfo(:)
      real*8, allocatable:: pot(:,:,:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
      real*8, allocatable:: potp(:,:)
      real*8, allocatable:: gradp(:,:,:)
      real*8, allocatable:: hessp(:,:,:)
      real*8, allocatable:: dpars(:)
      complex*16, allocatable:: zpars(:)
      real *8, allocatable:: grad_init(:,:,:,:)
      real *8, allocatable:: hess_init(:,:,:,:)
c     time quad
      character*1 type
      real*8, allocatable:: vpotinc(:,:,:), vpotpinc(:,:)
      external fforce
c
      pi=3.1415926535897932d0
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      if(ifvtarg .gt. 0) then
        ifpgh=1
      else
        ifpgh=0
      endif
c
      if(ifptarg .gt. 0) then
        ifpghtarg=1
      else
        ifpghtarg=0
      endif
c
      nleafbox = 0
      
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
      enddo
c
      allocate(fright(nd,npbox,nboxes))
      allocate(convective(nd,npbox,nboxes))
c      allocate(frightG(nd,npbox,nboxes))
      allocate(vpotinc(nd,npbox,nboxes))
      allocate(vpotpinc(nd,ntarg))

      allocate(tottimeinfo(100))
c
      allocate(pot(nd,npbox,ntarg))
      allocate(grad(nd,ndim,npbox,nboxes))
      allocate(hess(nd,nhess,npbox,nboxes))
      allocate(potp(nd,ntarg))
      allocate(gradp(nd,ndim,ntarg))
      allocate(hessp(nd,nhess,ntarg))
c
      allocate(ipars(100))
      allocate(dpars(1000))
      allocate(zpars(10))

      allocate(grad_init(nd,ndim,norder**ndim,nboxes))
      allocate(hess_init(nd,ndim*(ndim+1)/2,norder**ndim,nboxes))
c
      ipoly=1

c   1.   call fmm need hess
c
      ifpghdmk=3
      ifpghtargdmk=0
      ikernel=1
      beta=0.1d-5
c
      call bdmk_old_per(nd,ndim,eps,ikernel,beta,ipoly,norder,
     1    npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    initvalue,ifpghdmk,pot,grad,hess,ntarg,targs,
     3    ifpghtargdmk,potp,gradp,hessp,tottimeinfo)
       
    
c   3.  compute fs
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
              frightG(1,j,ibox) = hess(1,1,j,ibox)+hess(2,2,j,ibox)
              frightG(2,j,ibox) = hess(1,2,j,ibox)+hess(2,3,j,ibox)
              frightS(1,j,ibox) = initvalue(1,j,ibox)
     1               - frightG(1,j,ibox)
              frightS(2,j,ibox) = initvalue(2,j,ibox)
     1               - frightG(2,j,ibox)
             enddo
           endif
        enddo
      enddo


c
      end subroutine

c
      subroutine volpot(norder, nordert, ltree, itree, iptr,
     1           centers, nlevels,mxlevels, boxsize, nboxes,
     2           mxboxes, fforce, tpre,
     3           dt, iperiod, ntarg, targs, ifvtarg, ifptarg,
     4           eps, vpot, vpotp,visc,funs)
      implicit real*8 (a-h,o-z)
c     global vars
      integer norder, ltree, nboxes, ndim, nlevels
      integer nordert
      integer mxboxes,mxlevels
      integer iperiod, ntarg
      integer ifvtarg, ifptarg
      parameter(ndim=2)
      parameter(nd=2)
      integer itree(ltree), iptr(12)
      real*8 centers(ndim, mxboxes), boxsize(0:mxlevels)
      real*8 targs(ndim,ntarg), eps, dt
      real*8 targ(ndim)
      real*8 vpot(nd,norder**ndim,mxboxes), vpotp(nd,ntarg)
      real*8 visc
c     local vars
      integer, allocatable:: ipars(:)
      real*8, allocatable:: timeinfo(:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
      real*8, allocatable:: gradp(:,:,:)
      real*8, allocatable:: hessp(:,:,:)
      real*8, allocatable:: dpars(:)
      complex*16, allocatable:: zpars(:)
c     time quad
      character*1 type
      real*8, allocatable:: xtarg(:)
      real*8, allocatable:: tnodes(:), twhts(:)
      real*8, allocatable:: fright(:,:,:)
      real*8, allocatable:: frightS(:,:,:)
      real*8, allocatable:: frightG(:,:,:)
      real*8, allocatable:: vpotinc(:,:,:), vpotpinc(:,:)
      real*8, allocatable:: xref(:,:), wts(:)
      real *8,allocatable :: velocity_ex(:,:,:)
      external fforce funs

c      open(231,file='fs2.txt')
c
      pi=3.1415926535897932d0
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      if(ifvtarg .gt. 0) then
        ifpgh=1
      else
        ifpgh=0
      endif
c
      if(ifptarg .gt. 0) then
        ifpghtarg=1
      else
        ifpghtarg=0
      endif
c
      nleafbox = 0
      
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
      enddo
c
      allocate(tnodes(nordert))
      allocate(twhts(nordert))
      allocate(fright(nd,npbox,mxboxes))
      allocate(frightS(nd,npbox,mxboxes))
      allocate(frightG(nd,npbox,mxboxes))
      allocate(vpotinc(nd,npbox,mxboxes))
      allocate(vpotpinc(nd,ntarg))
c
      allocate(xref(ndim,npbox))
      allocate(wts(npbox))
c
      allocate(xtarg(ndim))
      allocate(timeinfo(100))
c
      allocate(grad(nd,ndim,npbox,mxboxes))
      allocate(hess(nd,nhess,npbox,mxboxes))
      allocate(gradp(nd,ndim,ntarg))
      allocate(hessp(nd,nhess,ntarg))
c
      allocate(ipars(100))
      allocate(dpars(1000))
      allocate(zpars(10))
      allocate(velocity_ex(nd,npbox,mxboxes))
c
      itype = 0
      ipoly=1
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c
      if(nordert .eq. 2) then
        nnodes=2
        tnodes(1)=0.0d0
        tnodes(2)=dt
        twhts(1)=1.0d0/2.0d0
        twhts(2)=1.0d0/2.0d0
      elseif(nordert .eq. 4) then
        nnodes=3
        tnodes(1)=0.0d0
        tnodes(2)=dt/2.0d0
        tnodes(3)=dt
        twhts(1)=1.0d0/6.0d0
        twhts(2)=4.0d0/6.0d0
        twhts(3)=1.0d0/6.0d0
      else
        norder=2
        nnodes=2
        tnodes(1)=0.0d0
        tnodes(2)=dt
        twhts(1)=1.0d0/2.0d0
        twhts(2)=1.0d0/2.0d0
        write(*,*) 'volpot: only norder=2,4 allowed!'
        write(*,*) 'norder reset to 2!'
      endif
c
c       twhts divided by dt
c
c--------------------------------------------------
c     initialize vpot to zero
c
      do i=1,nboxes
        do k=1,npbox
          vpot(1,k,i)=0.0d0
          vpot(2,k,i)=0.0d0
        enddo
      enddo
     
      do i = 1,ntarg
        vpotp(1,i)=0.0d0
        vpotp(2,i)=0.0d0
      enddo
c
c
c        write(*,*)norder,nordert,ltree
c        write(*,*)nlevels,mxlevels, nboxes,mxboxes
c        write(*,*)tpre,dt,ifvtarg,ifptarg,eps,npbox
c
c--------------------------------------------------
c     contribution from t(jj): an FGT from last step
      do jj=1,nnodes-1
        tt=tpre+tnodes(jj)
   
        scal=dt
c
c       assign fright: fforce(targ,tt,val) val->fright
        dpars(1)=tt
        ipoly=1
c
        call evalf(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,
     2       dpars,zpars,ipars,fforce,fright)
c       compute fs
        write(*,*)'compute helmholtz decomposition of f'
c
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           fright,eps,frightS,frightG)
c
c
c   -----------test fs ---------------------------------
c    exact solution for debuging
      ipoly=1
      dpars(1)=tt
      call evalf(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,
     2       dpars,zpars,ipars,funs,velocity_ex)

c
      sum_err=0.0d0
      sum_u=0.0d0
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               sum_u=sum_u+velocity_ex(1,j,ibox)**2
               sum_err=sum_err+abs(frightS(1,j,ibox)-
     1                   velocity_ex(1,j,ibox))**2
c             write(231,*)ibox,j,frightS(1,j,ibox),velocity_ex(1,j,ibox)
             enddo
          endif
        enddo
      enddo
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'frightS err2_rel=',err2_rel
c     ----------------------------------------------------------------

c       fright assigned
c
        do i=1,nboxes
          do j=1, npbox
            vpotinc(1,j,i)=0.0d0
            vpotinc(2,j,i)=0.0d0
          enddo
        enddo

        do i = 1,ntarg
          vpotpinc(1,i)=0.0d0
          vpotpinc(2,i)=0.0d0
        enddo
c
c       call boxfgt: fright -> vpotinc
c       be careful about the scaling factor

        delta=4.0d0*(dt-tnodes(jj))*visc
        ifnewtree = 0

        eps0=eps

c
      t1=omp_get_wtime()
        call boxfgt(nd,ndim,delta,eps0,ipoly,iperiod,norder,npbox,
     1       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2       frightS,ifpgh,vpotinc,grad,hess,ifnewtree,ntarg,targs,
     3       ifpghtarg,vpotpinc,gradp,hessp,timeinfo)
c
       t2=omp_get_wtime()
       write(*,*)'time for fgt=',t2-t1
       write(*,*)'speed for fgt=',npbox*nleafbox/(t2-t1)
c
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do j=1, npbox
               vpot(1,j,i)=vpot(1,j,i)+vpotinc(1,j,i)/(pi*delta)
     1            *dt*twhts(jj)

               vpot(2,j,i)=vpot(2,j,i)+vpotinc(2,j,i)/(pi*delta)
     1            *dt*twhts(jj)
            enddo
          endif
        enddo
c
      enddo
c
c--------------------------------------------------
c
c     contribution from t(nnodes): function value at
c     the current time step, no transform needed
      jj=nnodes
      tt=tpre+tnodes(jj)
c       assign fright: fforce(targ,tt,val) val->fright
        dpars(1)=tt
        ipoly=1
c
        call evalf(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,
     2       dpars,zpars,ipars,fforce,fright)
c
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           fright,eps,frightS,frightG)


      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               vpot(1,j,ibox)=vpot(1,j,ibox)+
     1                  frightS(1,j,ibox)*dt*twhts(jj)
               vpot(2,j,ibox)=vpot(2,j,ibox)+
     1                  frightS(2,j,ibox)*dt*twhts(jj)
             enddo
          endif
        enddo
      enddo
c
      if(ifptarg .gt. 0) then
c
      call evalt(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,vpot,
     2       ntarg,targs,vpotp)
      endif

c
      end subroutine


c--------------------------------------------------
c
      subroutine unsteady_stokes2dperfm4_video(nd,norder,nordert, finit,
     1           fforce,dt, ntot, eps, mxltree, mxboxes,
     2           mxlevels,ltree, itree, iptr, centers,
     3           nlevels,boxsize, nboxes, ifnewtree,velocity,vorticity,
     4           ntarg,velocityp,visc,funs)
      implicit real*8 (a-h,o-z)
      integer norder, ndim
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes,ntarg
      integer nordert
      integer nx,ny
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps
      real*8 visc
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 velocity(nd,norder**ndim,mxboxes)
      real*8 vorticity(norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
      complex*16 zpars(10), zk
c     allocatables
c     the tree structure
      real*8, allocatable:: fvals(:,:,:)
      real*8, allocatable:: fright(:,:,:)
      real*8, allocatable:: pot(:,:,:), vpot(:,:,:)
      real*8, allocatable:: potp(:,:), vpotp(:,:)
c     for tests only
      real*8, allocatable:: targs(:,:)
      real*8, allocatable:: targ(:)
      real *8, allocatable :: xref(:,:),wts(:)
c     user-specified points
c     solution at targets
      real*8 velocityp(nd,ntarg)
c     for extra targets
      real *8, allocatable :: coefsp(:,:,:)
      real *8, allocatable :: umat_nd(:,:,:)
      real *8, allocatable :: xq(:),umat(:,:),vmat(:,:)
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      character*1 type
c initialization
c
      real *8,allocatable :: M1(:,:,:)
      real *8,allocatable :: M2(:,:,:)
      real *8,allocatable :: M3(:,:,:)
      real *8,allocatable :: M4(:,:,:)
      real *8,allocatable :: usolp(:,:)
c
      real *8,allocatable :: usol1(:,:,:)
      real *8,allocatable :: usol2(:,:,:)
      real *8,allocatable :: usol3(:,:,:)
      real *8,allocatable :: usol4(:,:,:)

      real *8 dt_init
      integer ntot_init
      real *8,allocatable :: grad(:,:,:,:)
      real *8,allocatable :: hess(:,:,:,:)
      real *8,allocatable :: velocity_ex(:,:,:)

      real *8 xf(8),yf(8)

      external finit, fforce,uexact,funs
c------open output file:
      open(21,file='xf.dat')
      open(22,file='yf.dat')
      open(190,file='nleaves.dat')
      open(26,file='velocity.dat')
      open(27,file='gradientvx.dat')
      open(28,file='gradientuy.dat')
      open(29,file='gradientvy.dat')
      open(30,file='gradientux.dat')
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
c
c      write(*,*)'eps2 =',eps
c       epstree=1.0d-1*eps
c      epstree=eps*10
      epstree=eps*1000
      write(*,*)'epstree =',epstree
c      epstree=1.0d-4*eps
      eta=1.0d0
      zk=1.0d0
c      zk=60.0d0
c
      allocate(umat_nd(norder,norder,ndim))
      allocate(xq(norder),umat(norder,norder),
     1    vmat(norder,norder))
c
c     call vol_tree_mem2, the version that works better
c     for the adaptive case
c
      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1     norder,iptype,eta,finit,nd,dpars,zpars,ipars,
     2     ifnewtree,mboxes,mlevels,ltree,rintl)
c
      dpars(1)=dt
      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1     norder,iptype,eta,fforce,nd,dpars,zpars,ipars,
     2     ifnewtree,mboxes,mlevels,ltree,rintl)
      
      mboxes=mboxes*10
      mlevels=mlevels+5
      ltree=ltree*2

      write(*,*) 'mboxes=', mboxes
      write(*,*) 'mlevels=', mlevels
      write(*,*) 'ltree=', ltree
c     not sure if we need mxboxes, mxlevels, mxltree
c     as inputs
c
c     allocate mem
      allocate(fvals(nd,npbox,mxboxes))
      allocate(fright(nd,npbox,mxboxes))
c     the initial potential and volume potential
      allocate(pot(nd,npbox,mxboxes))
      allocate(vpot(nd,npbox,mxboxes))
c
c     what is ntarg though
      allocate(targs(ndim,ntarg))
      allocate(potp(nd,ntarg))
      allocate(vpotp(nd,ntarg))
c
c     why do we need these?
      allocate(targ(ndim))
      allocate(xref(ndim,npbox),wts(npbox))
c     extra targets
      allocate(coefsp(nd,npbox,mxboxes))
c
      allocate(M1(nd,norder**ndim,mxboxes))
c
      allocate(M2(nd,norder**ndim,mxboxes))
c
      allocate(M3(nd,norder**ndim,mxboxes))
c
      allocate(M4(nd,norder**ndim,mxboxes))
      allocate(usolp(nd,ntarg))
c
      allocate(usol1(nd,norder**ndim,mxboxes))
c
      allocate(usol2(nd,norder**ndim,mxboxes))
c
      allocate(usol3(nd,norder**ndim,mxboxes))
c
      allocate(usol4(nd,norder**ndim,mxboxes))
      
      allocate(grad(nd,ndim,npbox,mxboxes))
      allocate(hess(nd,nhess,npbox,mxboxes))
      allocate(velocity_ex(nd,npbox,mxboxes))

c
      nx=1000
      ny=1000
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
c     call vol_tree_build2 to resolve finit
c     again, the adaptive version
      write(*,*)'resolve initial function'
      call vol_tree_build2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,finit,nd,dpars,zpars,ipars,rintl,
     2    nboxes,mboxes,nlevels,mlevels,ltree,itree,iptr,centers,
     3    boxsize,fvals)
      write(*,*) 'nboxes, mboxes=', nboxes, mboxes
      write(*,*) 'nlevels, mlevels=', nlevels, mlevels
      write(*,*) 'ltree=', ltree
c
c      do ilevel=0,mlevels
c         rintl(ilevel)=0.70710678118654757
c      enddo
c
      do j=1, nboxes
        do i=1, npbox
          velocity(1,i,j)=fvals(1,i,j)
          velocity(2,i,j)=fvals(2,i,j)

        enddo
      enddo
 


c
       do it = 1,ntot

        tt=it*dt
        tpre=(it-1)*dt

        write(*,*) '*********************************'
        write(*,*) 'initial time step:', it
        write(*,*) 'time=', tt
        write(*,*) '*********************************'
c
          write(*,*) '--------------------'
          dpars(1)=tt
          call vol_tree_adap_fv(nd,ndim,epstree,ipoly,
     1         iptype,norder,npbox,nboxes,mboxes,nlevels,
     2         mlevels,ltree,itree,iptr,centers,boxsize,
     3         rintl,dpars,zpars,ipars,iperiod,fforce,velocity)

          write(*,*) 'nboxes, mboxes=', nboxes, mboxes
c
          nleafbox = 0
      
          do i=1,nboxes
            if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
          enddo
          write(*,*)'number of leaf box = ',nleafbox

          write(*,*) '--------------------'

c
        do j=1, nboxes
          do i=1, npbox
            fright(1,i,j)=velocity(1,i,j)
            fright(2,i,j)=velocity(2,i,j)
          enddo
        enddo

c
        do i=1, nboxes
        do j=1, npbox
          pot(1,j,i)=0.0d0
          pot(2,j,i)=0.0d0
          vpot(1,j,i)=0.0d0
          vpot(2,j,i)=0.0d0
        enddo
        enddo

c    compute volume potential
        write(*,*)'compute volume potential'
        ifvtarg=1
        ifptarg=0
c
        write(*,*)norder,nordert,ltree
        write(*,*)nlevels,mxlevels, nboxes,mxboxes
        write(*,*)tpre,dt,ifvtarg,ifptarg,eps

        call volpot(norder,nordert, ltree, itree, iptr,
     1       centers, nlevels,mxlevels, boxsize, nboxes,mxboxes,
     2       fforce, tpre,dt, iperiod, ntarg, targs, ifvtarg,
     3       ifptarg,eps,vpot,vpotp,visc,funs)


c    compute initial potential
        write(*,*)'compute initial potential'
      call compute_initpot(norder,ltree, itree, iptr, centers,
     1           nlevels,mxlevels, boxsize, nboxes,mxboxes, fright,
     2           dt,iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           eps, pot, potp,visc)

c
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
                velocity(1,j,ibox)=vpot(1,j,ibox) + pot(1,j,ibox)
                velocity(2,j,ibox)=vpot(2,j,ibox) + pot(2,j,ibox)
c                write(232,*)j,ibox,vpot(1,j,ibox),vpot(2,j,ibox),
c     1          pot(1,j,ibox),pot(2,j,ibox)
             enddo
           endif
        enddo
      enddo

c    exact solution
      ipoly=1
      dpars(1)=tt
      write(*,*)'final time =', dpars(1)
      call evalf(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,
     2       dpars,zpars,ipars,uexact,velocity_ex)

c
      sum_err=0.0d0
      sum_u=0.0d0
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               sum_u=sum_u+velocity_ex(2,j,ibox)**2
               sum_err=sum_err+abs(velocity(2,j,ibox)-
     1                   velocity_ex(2,j,ibox))**2
             enddo
          endif
        enddo
      enddo
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'usol err2_rel dimension2_pre=',err2_rel

c   output vorticity
      call evalgh(ndim,nd,ipoly,nlevels,ltree,itree,iptr,
     1       boxsize,norder,nboxes,velocity,grad,hess)
    

c   output velocity at target points
      call evalt(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,velocity,
     2       ntarg,targs,velocityp)

c
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
              vorticity(j,ibox)=grad(2,1,j,ibox)-
     1              grad(1,2,j,ibox)
             enddo
           endif
        enddo
      enddo
c    exact solution
      ipoly=1
      dpars(1)=tt
      write(*,*)'final time =', dpars(1)
      call evalf(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,
     2       dpars,zpars,ipars,funs,velocity_ex)

c
      sum_err=0.0d0
      sum_u=0.0d0
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               sum_u=sum_u+velocity_ex(1,j,ibox)**2
               sum_err=sum_err+abs(vorticity(j,ibox)-
     1                   velocity_ex(1,j,ibox))**2
             enddo
          endif
        enddo
      enddo
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'vorticity err2_rel =',err2_rel
c
      iplot=1
      if(iplot .eq. 1)then

      if(mod(it,1) .eq. 0) then

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
         write(26,*)vorticity(8*(k-1)+1,ibox),vorticity(8*(k-1)+2,ibox),
     1 vorticity(8*(k-1)+3,ibox),vorticity(8*(k-1)+4,ibox),
     2 vorticity(8*(k-1)+5,ibox),vorticity(8*(k-1)+6,ibox),
     3 vorticity(8*(k-1)+7,ibox),vorticity(8*(k-1)+8,ibox)
       end do
c
       do k = 1, 8
         write(27,*)grad(2,1,8*(k-1)+1,ibox),grad(2,1,8*(k-1)+2,ibox),
     1 grad(2,1,8*(k-1)+3,ibox),grad(2,1,8*(k-1)+4,ibox),
     2 grad(2,1,8*(k-1)+5,ibox),grad(2,1,8*(k-1)+6,ibox),
     3 grad(2,1,8*(k-1)+7,ibox),grad(2,1,8*(k-1)+8,ibox)
       end do
c
       do k = 1, 8
         write(28,*)grad(1,2,8*(k-1)+1,ibox),grad(1,2,8*(k-1)+2,ibox),
     1 grad(1,2,8*(k-1)+3,ibox),grad(1,2,8*(k-1)+4,ibox),
     2 grad(1,2,8*(k-1)+5,ibox),grad(1,2,8*(k-1)+6,ibox),
     3 grad(1,2,8*(k-1)+7,ibox),grad(1,2,8*(k-1)+8,ibox)
       end do
c
       do k = 1, 8
         write(29,*)grad(2,2,8*(k-1)+1,ibox),grad(2,2,8*(k-1)+2,ibox),
     1 grad(2,2,8*(k-1)+3,ibox),grad(2,2,8*(k-1)+4,ibox),
     2 grad(2,2,8*(k-1)+5,ibox),grad(2,2,8*(k-1)+6,ibox),
     3 grad(2,2,8*(k-1)+7,ibox),grad(2,2,8*(k-1)+8,ibox)
       end do
c
       do k = 1, 8
         write(30,*)grad(1,1,8*(k-1)+1,ibox),grad(1,1,8*(k-1)+2,ibox),
     1 grad(1,1,8*(k-1)+3,ibox),grad(1,1,8*(k-1)+4,ibox),
     2 grad(1,1,8*(k-1)+5,ibox),grad(1,1,8*(k-1)+6,ibox),
     3 grad(1,1,8*(k-1)+7,ibox),grad(1,1,8*(k-1)+8,ibox)
       end do

          endif
        enddo
      enddo

      endif

      endif
c   end of it
       enddo

c
      close(21)
      close(22)
      close(190)
      close(26)


c
      end subroutine
