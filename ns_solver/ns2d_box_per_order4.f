c u0: the solution at t_1
c u1: the solution at t_1

      subroutine ns2dperfm4(nd,norder,nordert, finit,
     1           fforce,dt, ntot, eps, mxltree, mxboxes,
     2           mxlevels,ltree, itree, iptr, centers,
     3           nlevels,boxsize, nboxes, ifnewtree,velocity,vorticity,
     4           ntarg,velocityp,visc,irhsfun)
      implicit real*8 (a-h,o-z)
      integer norder, ndim,nordert
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes,ntarg
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
      real *8,allocatable :: u0(:,:,:)
      real *8,allocatable :: u1(:,:,:)
      real *8,allocatable :: u2(:,:,:)
      real *8,allocatable :: u3(:,:,:)
      real *8,allocatable :: usolp(:,:)


      real *8 dt_init
      integer ntot_init
      real *8,allocatable :: grad(:,:,:,:)
      real *8,allocatable :: hess(:,:,:,:)
      real *8,allocatable :: velocity_ex(:,:,:)

      external finit, fforce,uexact
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
      ifvtarg=1
      ifptarg=0
c
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
c
      boxlen=1.0d0
c
c
      epstree=1.0d-1*eps
c      epstree=eps
      write(*,*)'epstree=',epstree
c      epstree=1.0d-4*eps
      eta=1.0d0
      zk=1.0d0
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
      if(irhsfun .gt. 0)then
       dpars(1)=dt
      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1     norder,iptype,eta,fforce,nd,dpars,zpars,ipars,
     2     ifnewtree,mboxes,mlevels,ltree,rintl)
      endif

      mboxes=mboxes*6
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
      allocate(u0(nd,norder**ndim,mxboxes))
c
      allocate(u1(nd,norder**ndim,mxboxes))
c
      allocate(u2(nd,norder**ndim,mxboxes))
c
      allocate(u3(nd,norder**ndim,mxboxes))
      allocate(usolp(nd,ntarg))
      
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
      dt0=dt/8.0d0
      ntot0=24
      call build_init_tree(nd,norder, ltree, itree, iptr,
     1 centers,nlevels,mlevels, boxsize, nboxes,mboxes,fforce,
     2   dt0,ntot0,iperiod, ntarg, targs, ifvtarg, ifptarg,
     3   eps,rintl, fvals,visc)

c    exact solution
      ipoly=1
      dpars(1)=0.0d0
      call evalf(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,
     2       dpars,zpars,ipars,finit,fvals)

c
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               u0(1,j,ibox)=fvals(1,j,ibox)
               u0(2,j,ibox)=fvals(2,j,ibox)
             enddo
          endif
        enddo
      enddo

c     initialization of u1,u2,u3
      call ns2d_init_usol(nd,norder,nordert,
     1           fforce,dt, eps, mxltree, mboxes,
     2           mlevels,rintl,ltree, itree, iptr, centers,
     3           nlevels,boxsize, nboxes, ifnewtree,ntarg,targs,
     4           u0,u1,u2,u3,visc)





ccccccccccccc call volpot_main
c
      iadap=1
      call ns2d_am24_main(nd,norder, nordert, fforce,
     1     dt, ntot, eps, mxltree,mboxes,mlevels,rintl,
     2     ltree, itree, iptr, centers, nlevels,
     3     boxsize, nboxes, ifnewtree,iadap,ntarg,targs,
     4     u0, u1,u2,u3,visc)
c
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               velocity(1,j,ibox)=u3(1,j,ibox)
               velocity(2,j,ibox)=u3(2,j,ibox)
             enddo
          endif
        enddo
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
             enddo
           endif
        enddo
      enddo

       end subroutine



c
      subroutine build_init_tree(nd,norder, ltree, itree, iptr,
     1 centers,nlevels,mlevels, boxsize, nboxes,mboxes,fforce,
     2   dt,ntot,iperiod, ntarg, targs, ifvtarg, ifptarg,
     3   eps,rintl, initvalue,visc)
      implicit real*8 (a-h,o-z)
c     global vars
      integer nd
      integer norder, ltree, nboxes, ndim, nlevels
      integer mboxes,mlevels
      integer iperiod, ntarg
      integer ifvtarg, ifptarg
      parameter(ndim=2)
      real*8 rintl(0:200)
      integer itree(ltree), iptr(12)
      real*8 centers(ndim, mboxes), boxsize(0:mlevels)
      real*8 dt
      real*8,allocatable :: fright(:,:,:)
      real*8 targs(ndim,ntarg), eps
      real *8 initvalue(nd,norder**ndim,mboxes)
      real *8 usol(nd,norder**ndim,mboxes)
      real *8 usolp(nd,ntarg)
      real *8 visc
c     local vars
      real*8, allocatable:: timeinfo(:)
      real*8, allocatable:: pot(:,:,:)
      real *8, allocatable:: potp(:,:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
      real*8, allocatable:: gradp(:,:,:)
      real*8, allocatable:: hessp(:,:,:)
      real*8, allocatable:: vpot(:,:,:)
      real*8, allocatable:: vpotp(:,:)
c     initial potential
c      real *8 initpot(nd,norder**ndim,nboxes)
c     real *8 initpotp(nd,ntarg)
      real *8, allocatable:: initpot(:,:,:)
      real *8, allocatable:: initpotp(:,:)
      real *8, allocatable:: frightS1(:,:,:)
      real *8, allocatable:: frightS2(:,:,:)

      integer, allocatable:: ipars(:)
      real*8, allocatable:: dpars(:)
      complex*16, allocatable:: zpars(:)
      external fforce
c
      pi=3.1415926535897932d0
      nhess = ndim*(ndim+1)/2
      npbox=norder**ndim
c
      allocate(timeinfo(100))
c
c      allocate(initvalue(nd,norder**ndim,nboxes))

c
      allocate(pot(nd,npbox,mboxes))
      allocate(potp(nd,ntarg))
      allocate(vpot(nd,npbox,mboxes))
      allocate(vpotp(nd,ntarg))
      allocate(grad(nd,ndim,npbox,mboxes))
      allocate(hess(nd,nhess,npbox,mboxes))
      allocate(gradp(nd,ndim,ntarg))
      allocate(hessp(nd,nhess,ntarg))
c
      allocate(fright(nd,npbox,mboxes))
      allocate(initpot(nd,npbox,mboxes))
      allocate(initpotp(nd,ntarg))
      allocate(frightS1(nd,npbox,mboxes))
      allocate(frightS2(nd,npbox,mboxes))
c
      allocate(ipars(100))
      allocate(dpars(1000))
      allocate(zpars(10))

c
c     allocate certain vars
c
c     simply call boxfgt and divide the
c     output by dt
c
      epstree=1.0d-1*eps
c      epstree=eps
      write(*,*)epstree
      ifnewtree=0
      ipoly=1
      iptype=2
c
        do j=1, nboxes
          do i=1, npbox
            fright(1,i,j)=initvalue(1,i,j)
            fright(2,i,j)=initvalue(2,i,j)
          enddo
        enddo
c
      do it=1, ntot
        tt=it*dt
        tpre=(it-1)*dt
c
          dpars(1)=dt
          call vol_tree_adap_fv(nd,ndim,epstree,ipoly,
     1         iptype,norder,npbox,nboxes,mboxes,nlevels,
     2         mlevels,ltree,itree,iptr,centers,boxsize,
     3         rintl,dpars,zpars,ipars,iperiod,fforce,fright)

         write(*,*) 'nboxes, mboxes=', nboxes, mboxes
c
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tpre,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           fright,eps,frightS1)
c
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           fright,eps,frightS2)


c    compute volume potential
c
c      write(*,*)'call volpot1'
      call volpot1(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes,
     2           dt, iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           frightS1,frightS2,eps, vpot, vpotp,visc)

c      write(*,*)'call compute_initpot'
c    compute initial potential
      call compute_initpot(norder, ltree, itree, iptr, centers,
     1  nlevels,mxlevels, boxsize, nboxes,mxboxes, fright, dt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           eps, initpot, initpotp,visc)
c
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
                usol(1,j,ibox)=vpot(1,j,ibox) + initpot(1,j,ibox)
                usol(2,j,ibox)=vpot(2,j,ibox) + initpot(2,j,ibox)
                fright(1,j,ibox)=usol(1,j,ibox)
                fright(2,j,ibox)=usol(2,j,ibox)
             enddo
           endif
        enddo
      enddo


c  end of it
      enddo


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
     1           centers, nlevels, boxsize, nboxes, fforce,tt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           initvalue,eps,frightS)
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
c     local vars
      real *8, allocatable:: fright(:,:,:)
      real *8, allocatable:: convective(:,:,:)
      real *8, allocatable:: frightG(:,:,:)
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
      allocate(fright(nd,npbox,nboxes))
      allocate(convective(nd,npbox,nboxes))
      allocate(frightG(nd,npbox,nboxes))
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

c     compute convective term at t_{n-1}
      call evalgh(ndim,nd,ipoly,nlevels,ltree,itree,iptr,
     1       boxsize,norder,nboxes,initvalue,grad_init,hess_init)

      dpars(1)=tt
c
      call evalf(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,
     2       dpars,zpars,ipars,fforce,fright)

c
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
                convective(1,j,ibox) = fright(1,j,ibox)-
     1     initvalue(1,j,ibox)*grad_init(1,1,j,ibox)-
     2     initvalue(2,j,ibox)*grad_init(1,2,j,ibox)


                convective(2,j,ibox) = fright(2,j,ibox)-
     1     initvalue(1,j,ibox)*grad_init(2,1,j,ibox)-
     2     initvalue(2,j,ibox)*grad_init(2,2,j,ibox)
             enddo
          endif
        enddo
      enddo
c   2.   call fmm need hess
c
      ifpghdmk=3
      ifpghtargdmk=0
      ikernel=1
      beta=0.1d-5
c
c      write(*,*)'dmkdmk'
      call bdmk_old_per(nd,ndim,eps,ikernel,beta,ipoly,norder,
     1    npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    convective,ifpghdmk,pot,grad,hess,ntarg,targs,
     3    ifpghtargdmk,potp,gradp,hessp,tottimeinfo)
c       write(*,*)'koko'
    
c   3.  compute fs
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
              frightG(1,j,ibox) = hess(1,1,j,ibox)+hess(2,2,j,ibox)
              frightG(2,j,ibox) = hess(1,2,j,ibox)+hess(2,3,j,ibox)

              frightS(1,j,ibox) = convective(1,j,ibox)
     1               - frightG(1,j,ibox)
              frightS(2,j,ibox) = convective(2,j,ibox)
     1               - frightG(2,j,ibox)
             enddo
           endif
        enddo
      enddo


c
      end subroutine

c**************************************************************
c     subroutine volpot1
c**************************************************************
c     This subroutine computes the volume potential
c
c     trapezoidal rule and simpson's rule for the time quad used
c     which gives 2nd and 4th order in time
c
c  INPUT:
c  norder: degree of poly approx on each leaf node
c  ltree - nboxes: the tree structure
c  fforce: the handle for the forcing term
c          (can be replaced later)
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
c  vpot: the initial potential on the volume grid
c  vpotp: the initial potential at extra targets
c
c  attention: vpotp not assigned yet in this version
c
c**************************************************************
c
      subroutine volpot1(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes,
     2           dt, iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           fright1,fright2,eps, vpot, vpotp,visc)
      implicit real*8 (a-h,o-z)
c     global vars
      parameter(nd=2)
      integer norder, ltree, nboxes, ndim, nlevels
      integer iperiod, ntarg
      integer ifvtarg, ifptarg
      parameter(ndim=2)
      integer itree(ltree), iptr(12)
      real*8 centers(ndim, nboxes), boxsize(0:nlevels)
      real*8 fright1(nd,norder**ndim,nboxes)
      real*8 fright2(nd,norder**ndim,nboxes)
      real*8 targs(ndim,ntarg), eps, dt
      real*8 targ(ndim)
      real*8 vpot(nd,norder**ndim,nboxes), vpotp(nd,ntarg)
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
      allocate(vpotinc(nd,npbox,nboxes))
      allocate(vpotpinc(nd,ntarg))

      allocate(timeinfo(100))
c
      allocate(grad(nd,ndim,npbox,nboxes))
      allocate(hess(nd,nhess,npbox,nboxes))
      allocate(gradp(nd,ndim,ntarg))
      allocate(hessp(nd,nhess,ntarg))
c
      allocate(ipars(100))
      allocate(dpars(1000))
      allocate(zpars(10))
c
      ipoly=1
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
c--------------------------------------------------
c     contribution from t(jj): an FGT from last step
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
        delta=4.0d0*dt*visc
        ifnewtree = 0
        t1=omp_get_wtime()
        call boxfgt(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2       fright1,ifpgh,vpotinc,grad,hess,ifnewtree,ntarg,targs,
     3       ifpghtarg,vpotpinc,gradp,hessp,timeinfo)
        t2=omp_get_wtime()
        write(*,*)'time for fgt=',t2-t1

c
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do j=1, npbox
               vpot(1,j,i)=vpot(1,j,i)+vpotinc(1,j,i)/(pi*delta)
     1                  *(dt/2.0d0)
               vpot(2,j,i)=vpot(2,j,i)+vpotinc(2,j,i)/(pi*delta)
     1                  *(dt/2.0d0)
c               write(*,*)vpot(1,j,i),vpot(2,j,i)
            enddo
          endif
        enddo
c     contribution from t(nnodes): function value at
c     the current time step, no transform needed
        do i=1, nboxes
            do j=1, npbox
               vpot(1,j,i)=vpot(1,j,i)+fright2(1,j,i)*(dt/2.0d0)
               vpot(2,j,i)=vpot(2,j,i)+fright2(2,j,i)*(dt/2.0d0)
c            write(*,*)vpot(1,j,i),fright2(1,j,i)*(dt/2.0d0)
            enddo
        enddo
       
c      for target point
c
      if(ifptarg .gt. 0) then
c
      call evalt(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,vpot,
     2       ntarg,targs,vpotp)
      endif
 
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
       t2=omp_get_wtime()
       write(*,*)'time for fgt=',t2-t1
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

c     initialization of u1,u2,u3
      subroutine ns2d_init_usol(nd,norder,nordert,
     1           fforce,dt0, eps, mxltree, mxboxes,
     2           mxlevels,rintl,ltree, itree, iptr, centers,
     3           nlevels,boxsize, nboxes, ifnewtree,ntarg,
     4           targs,u0,u1,u2,u3,visc)
      implicit real*8 (a-h,o-z)
      integer norder,nordert
      integer nd,ntarg
      parameter(ndim=2)
      integer  mxltree, mxboxes
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt0, eps,dt
      integer ntot
      real*8 visc
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 targs(ndim,ntarg)
      real*8 rintl(0:200), dpars(1000)
      complex*16 zpars(10), zk
      real*8 u0(nd,norder**ndim,mxboxes)
      real*8 u1(nd,norder**ndim,mxboxes)
      real*8 u2(nd,norder**ndim,mxboxes)
      real*8 u3(nd,norder**ndim,mxboxes)
c     richardson related
      integer nr(4), ntots(4)
      real*8 rwhts(4), dts(4)
c     allocatables, not sure yet
      real*8, allocatable::usave(:,:,:,:)
      real*8,allocatable :: fright(:,:,:)
c     local vars
      real*8, allocatable:: timeinfo(:)
      real*8, allocatable:: pot(:,:,:)
      real *8, allocatable:: potp(:,:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
      real*8, allocatable:: gradp(:,:,:)
      real*8, allocatable:: hessp(:,:,:)
      real*8, allocatable:: vpot(:,:,:)
      real*8, allocatable:: vpotp(:,:)
c     initial potential
c      real *8 initpot(nd,norder**ndim,nboxes)
c     real *8 initpotp(nd,ntarg)
      real *8, allocatable:: initpot(:,:,:)
      real *8, allocatable:: initpotp(:,:)
      real *8, allocatable:: frightS1(:,:,:)
      real *8, allocatable:: frightS2(:,:,:)
       
      external  fforce
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
      epstree=1.0d-1*eps
c      epstree=1.0d-4*eps
      eta=1.0d0
      zk=1.0d0


c
      allocate(usave(nd,npbox,mxboxes,3))
      allocate(fright(nd,npbox,mxboxes))
c
      allocate(timeinfo(100))
c
c      allocate(initvalue(nd,norder**ndim,nboxes))

c
      allocate(pot(nd,npbox,mxboxes))
      allocate(potp(nd,ntarg))
      allocate(vpot(nd,npbox,mxboxes))
      allocate(vpotp(nd,ntarg))
      allocate(grad(nd,ndim,npbox,mxboxes))
      allocate(hess(nd,nhess,npbox,mxboxes))
      allocate(gradp(nd,ndim,ntarg))
      allocate(hessp(nd,nhess,ntarg))
c
      allocate(initpot(nd,npbox,mxboxes))
      allocate(initpotp(nd,ntarg))
      allocate(frightS1(nd,npbox,mxboxes))
      allocate(frightS2(nd,npbox,mxboxes))

c     richardson parameters
c     do three passes of the 2nd order solve
c     set parameters for each
      nr(1)=1
      nr(2)=2
      nr(3)=4
      nr(4)=8
c     number of refinements for each pass
c
      dts(1)=dt0/nr(1)
      dts(2)=dt0/nr(2)
      dts(3)=dt0/nr(3)
      dts(4)=dt0/nr(4)
c     smaller delta_t values
c
      ntot0=3
      ntots(1)=ntot0*nr(1)
      ntots(2)=ntot0*nr(2)
      ntots(3)=ntot0*nr(3)
      ntots(4)=ntot0*nr(4)
c     total num of (smaller) steps
c
      rwhts(1)=-1.0d0/21.0d0
      rwhts(2)=2.0d0/3.0d0
      rwhts(3)=-8.0d0/3.0d0
      rwhts(4)=64.0d0/21.0d0
c     weights for each Richardson pass

c     call ns2d_init_tree


c     ---------------------------

c     initialize usave to be zero
      do jj=1, 3
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
              do j=1, npbox
                usave(1,j,ibox,jj)=0.0d0
                usave(2,j,ibox,jj)=0.0d0
              enddo
            endif
          enddo
        enddo
      enddo

    


c------------------------------------------------------
c
c     do the three passes of Richarson solve
c     2nd order -> 3rd order -> 4th order
      do jp=1, 4
        dt=dts(jp)
        ntot=ntots(jp)
c     initialize fright to u0
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
              do j=1, npbox
                fright(1,j,ibox)=u0(1,j,ibox)
                fright(2,j,ibox)=u0(2,j,ibox)
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
c      write(*,*)'call volpot1'
c
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tpre,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           fright,eps,frightS1)
c
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           fright,eps,frightS2)

      call volpot1(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes,
     2           dt, iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           frightS1,frightS2,eps, vpot, vpotp,visc)

c    compute initial potential
      call compute_initpot(norder, ltree, itree, iptr, centers,
     1  nlevels,mxlevels, boxsize, nboxes,mxboxes, fright, dt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           eps, initpot, initpotp,visc)

c
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
                fright(1,j,ibox)=vpot(1,j,ibox) + initpot(1,j,ibox)
                fright(2,j,ibox)=vpot(2,j,ibox) + initpot(2,j,ibox)
             enddo
           endif
        enddo
      enddo
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
              usave(1,i,ibox,ii)=usave(1,i,ibox,ii)+rwhts(jp)*
     1                                    fright(1,i,ibox)
              usave(2,i,ibox,ii)=usave(2,i,ibox,ii)+rwhts(jp)*
     1                                    fright(2,i,ibox)
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
      if(nordert .eq. 4) then
        do ibox=1, nboxes
        do i=1, npbox
          do k =1,nd
          u1(k,i,ibox)=usave(k,i,ibox,1)
          u2(k,i,ibox)=usave(k,i,ibox,2)
          u3(k,i,ibox)=usave(k,i,ibox,3)
          enddo
        enddo
        enddo
      else if(nordert .eq. 2) then
        do ibox=1, nboxes
        do i=1, npbox
          do k =1,nd
          u1(k,i,ibox)=usave(k,i,ibox,1)
          u2(k,i,ibox)=0.0d0
          u3(k,i,ibox)=0.0d0
          enddo
        enddo
        enddo
      endif
c     like this?
c     yeah, like this
c
      end subroutine

c

c
      subroutine ns2d_am24_main(nd,norder, nordert, fforce,
     1     dt, ntot, eps, mxltree,mxboxes,mxlevels,rintl,
     2     ltree, itree, iptr, centers, nlevels,
     3     boxsize, nboxes, ifnewtree,iadap,ntarg,targs,
     4     u0, u1,u2,u3,visc)
      implicit real*8 (a-h,o-z)
      integer norder, ndim,nd,ntarg
      parameter(ndim=2)
      integer ntot, mxltree, mxboxes,mxlevels
      integer ipars(100), iptr(12)
      integer itree(mxltree)
      real*8 dt, eps,visc
      real*8 targs(ndim,ntarg)
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 u0(nd,norder**ndim,mxboxes)
      real*8 u1(nd,norder**ndim,mxboxes)
      real*8 u2(nd,norder**ndim,mxboxes)
      real*8 u3(nd,norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
      complex*16 zpars(10), zk
c     allocatables
c     the tree structure
      real *8, allocatable :: xref(:,:),wts(:)
      real *8 xf(8),yf(8)
      real*8,allocatable :: vorticity(:,:)
      real *8,allocatable :: grad(:,:,:,:)
      real *8,allocatable :: hess(:,:,:,:)
      real *8,allocatable :: velocity_ex(:,:,:)
      character*1 type
      real*8, allocatable:: targ(:)
c     for tests only

      external fforce,uexact
c------open output file:
      open(21,file='xf.dat')
      open(22,file='yf.dat')
      open(190,file='nleaves.dat')
      open(26,file='velocity.dat')
c      open(299,file='u3.txt')
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
c      epstree=1.0d-1*eps
c      epstree=1.0d-4*eps
ccc   tree refined for debugging
      eta=1.0d0
      zk=1.0d0
      allocate(xref(ndim,npbox),wts(npbox))
c
      allocate(targ(ndim))
      allocate(grad(nd,ndim,npbox,mxboxes))
      allocate(hess(nd,nhess,npbox,mxboxes))
      allocate(vorticity(npbox,mxboxes))
      allocate(velocity_ex(nd,npbox,mxboxes ))

c
c     time stepping:
      do it=nordert, ntot
        tt=it*dt
        tpre=(it-nordert)*dt
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
c        write(*,*) nd, ndim, eps,visc
c        write(*,*) ipoly, iperiod, norder, npbox
c        write(*,*) nboxes, nlevels, ifvtarg
c        write(*,*) ntarg, ifptarg,mxltree,mxboxes,mxlevels
c        write(*,*)nordert,ltree,tpre,dt,iadap,ifnewtree
        t1=omp_get_wtime()
        call vpot_ns_am24(nd,norder, nordert, ltree, itree,
     1           iptr, centers,mxltree, nlevels,mxlevels, boxsize,
     2         nboxes,mxboxes,fforce, tpre, dt, iperiod, ntarg, targs,
     3     rintl,ifvtarg, ifptarg,ifnewtree,iadap, eps,visc,
     4     u0,u1,u2,u3)
         t2=omp_get_wtime()
c
      nleafbox = 0
      
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
      enddo
c
        write(*,*)'time for one step=',t2-t1
        write(*,*)'speed for one step=',npbox*nleafbox/(t2-t1)
c
      nleafbox = 0
      
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
      enddo
      write(*,*)'nboxes=',nboxes
      write(*,*)'nleaves=',nleafbox
      write(*,*)'nlevels=',nlevels

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
               sum_err=sum_err+abs(u3(1,j,ibox)-
     1                   velocity_ex(1,j,ibox))**2
             enddo
          endif
        enddo
      enddo
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'usol err2_rel=',err2_rel,sum_err,sum_u
c
c      do ilevel=0,nlevels
c        bs = boxsize(ilevel)/2.0d0
c        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
c          if(itree(iptr(4)+ibox-1).eq.0) then
c             do j=1,npbox
c               write(299,*)it,ibox,j,u3(1,j,ibox),u3(2,j,ibox)
c             enddo
c           endif
c        enddo
c      enddo


      iplot=1
      if(iplot .eq. 1) then
c
      if(mod(it,15) .eq. 0) then
c   output vorticity
      write(*,*)'plot for vorticity',it
      call evalgh(ndim,nd,ipoly,nlevels,ltree,itree,iptr,
     1       boxsize,norder,nboxes,u3,grad,hess)

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

      itype = 0
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)

      write(*,*)'output data files to plot the tree'

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
          endif
        enddo
      enddo

      endif

      endif

      enddo
c

      end subroutine






c**************************************************************
      subroutine vpot_ns_am24(nd,norder, nordert, ltree, itree,
     1       iptr, centers,mxltree, nlevels,mxlevels, boxsize,
     2       nboxes,mxboxes,fforce, tpre, dt, iperiod, ntarg, targs,
     3   rintl, ifvtarg, ifptarg,ifnewtree,iadap, eps,visc,
     4    u0,u1,u2,u3)
      implicit real*8 (a-h,o-z)
c     global vars
      integer nd
      integer norder, ltree, nboxes, ndim, nlevels
      integer mxltree,mxlevels,mxboxes
      integer iperiod, ntarg
      integer ifvtarg, ifptarg
      integer iadap
      parameter(ndim=2)
      integer itree(ltree), iptr(12)
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 targs(ndim,ntarg), eps, dt
      real*8 u0(nd,norder**ndim,mxboxes)
      real*8 u1(nd,norder**ndim,mxboxes)
      real*8 u2(nd,norder**ndim,mxboxes)
      real*8 u3(nd,norder**ndim,mxboxes)
      real*8 rintl(0:200)
      real*8 visc

c     local vars
      character*1 type
      real*8, allocatable:: vpotinc(:,:,:), vpotpinc(:,:)
      real*8, allocatable:: tnodes(:), twhts(:)

      integer, allocatable:: ipars(:)
      real*8, allocatable:: timeinfo(:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
      real*8, allocatable:: gradp(:,:,:)
      real*8, allocatable:: hessp(:,:,:)
      real*8, allocatable:: dpars(:),unm1(:,:,:)
      complex*16, allocatable:: zpars(:)
c
      real*8, allocatable:: fnm1(:,:,:)
      real*8, allocatable:: fnm2(:,:,:)
      real*8, allocatable:: fnm3(:,:,:)
      real*8, allocatable:: fnm4(:,:,:)
      real*8, allocatable:: fright(:,:,:)
c
      real*8, allocatable:: pot(:,:,:), vpot(:,:,:)
      real*8, allocatable:: potp(:,:), vpotp(:,:)
      real *8, allocatable:: velocity_ex(:,:,:)
c
      real*8, allocatable:: vpot4(:,:,:,:)
      real*8, allocatable:: vpot_correction(:,:,:)
      real*8, allocatable:: initpot(:,:,:)
      external fforce,funs

c
      pi=3.1415926535897932d0
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
      iptype=2
      ipoly=1
c
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
      allocate(tnodes(nordert))
      allocate(twhts(nordert))
c
      allocate(vpotinc(nd,npbox,mxboxes))
      allocate(vpotpinc(nd,ntarg))

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
      allocate(unm1(nd,npbox,mxboxes))
c
      allocate(fnm1(nd,npbox,mxboxes))
      allocate(fnm2(nd,npbox,mxboxes))
      allocate(fnm3(nd,npbox,mxboxes))
      allocate(fnm4(nd,npbox,mxboxes))
      allocate(fright(nd,npbox,mxboxes))
c     the initial potential and volume potential
      allocate(pot(nd,npbox,mxboxes))
      allocate(vpot(nd,npbox,mxboxes))
      allocate(potp(nd,ntarg))
      allocate(vpotp(nd,ntarg))
      allocate(velocity_ex(nd,npbox,mxboxes))
      allocate(vpot4(nd,npbox,mxboxes,4))
      allocate(vpot_correction(nd,npbox,mxboxes))
      allocate(initpot(nd,npbox,mxboxes))

c
      if(nordert .eq. 2) then
        nnodes=2
        tnodes(1)=0.0d0
        tnodes(2)=dt
        twhts(1)=-1.0d0/2.0d0
        twhts(2)=3.0d0/2.0d0
      elseif(nordert .eq. 4) then
        nnodes=4
        tnodes(1)=0.0d0
        tnodes(2)=dt
        tnodes(3)=2*dt
        tnodes(4)=3*dt
        twhts(1)= -3.0d0/8.0d0
        twhts(2)=37.0d0/24.0d0
        twhts(3)=-59.0d0/24.0d0
        twhts(4)= 55.0d0/24.0d0
      else
        nordert=2
        nnodes=2
        tnodes(1)=0.0d0
        tnodes(2)=dt
        twhts(1)=-1.0d0/2.0d0
        twhts(2)=3.0d0/2.0d0
        write(*,*) 'volpot: only nordert=2,3,4 allowed!'
        write(*,*) 'nordert reset to 2!'
      endif

c     initialize vpot, fn to zero
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
      do i=1,nboxes
         do j=1,npbox
            do k=1, nd
              fnm1(k,j,i)=0.0d0
              fnm2(k,j,i)=0.0d0
              fnm3(k,j,i)=0.0d0
              fnm4(k,j,i)=0.0d0
              unm1(k,j,i)=u3(k,j,i)
              fright(k,j,i)=0.0d0
            enddo
         enddo
      enddo
c
      nleafbox = 0
      
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
      enddo

c--------------------------------------------------
c     contribution from t(jj): an FGT from last step
      scal=dt
      do jj=1,nnodes
c       loop over nnodes-1 nodes (in time)
        tt=tpre+tnodes(jj)
c        write(*,*) tpre, tnodes(jj), tt
        delta=(nordert*dt-tnodes(jj))*visc*4.0d0
c        write(*,*) 'jj, nnodes', jj, nnodes
c
c
      do i=1,nboxes
        do k=1,npbox
          fright(1,k,i)=0.0d0
          fright(2,k,i)=0.0d0
        enddo
      enddo

      if(jj .eq. 1) then
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           u0,eps,fnm1)
      
c
          do i=1,nboxes
            do k=1,npbox
              fright(1,k,i)=fnm1(1,k,i)
              fright(2,k,i)=fnm1(2,k,i)
c              write(203,*)i,k,fnm1(1,k,i),fnm1(2,k,i)
            enddo
          enddo

c
      elseif(jj .eq. 2) then
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           u1,eps,fnm2)
      
c
          do i=1,nboxes
            do k=1,npbox
              fright(1,k,i)=fnm2(1,k,i)
              fright(2,k,i)=fnm2(2,k,i)
c              write(202,*)i,k,fnm2(1,k,i),fnm2(2,k,i)
            enddo
          enddo
c
      elseif(jj .eq. 3) then
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           u2,eps,fnm3)
      
c
          do i=1,nboxes
            do k=1,npbox
              fright(1,k,i)=fnm3(1,k,i)
              fright(2,k,i)=fnm3(2,k,i)
c              write(201,*)i,k,fnm3(1,k,i),fnm3(2,k,i)
            enddo
          enddo
c
      elseif(jj .eq. 4) then
c
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           u3,eps,fnm4)


c
          do i=1,nboxes
            do k=1,npbox
              fright(1,k,i)=fnm4(1,k,i)
              fright(2,k,i)=fnm4(2,k,i)
c              write(204,*)i,k,fnm4(1,k,i),fnm4(2,k,i)
            enddo
          enddo
c    end of if jj
       endif
c
        do i=1, nboxes
          do k=1, npbox
            vpotinc(1,k,i)=0.0d0
            vpotinc(2,k,i)=0.0d0
          enddo
        enddo

c        delta = 4.0d0*delta
        ifnewtree = 0
c        write(*,*) 'calling boxfgt'

c
c        write(*,*) nd, ndim, delta, eps
c        write(*,*) ipoly, iperiod, norder, npbox
c        write(*,*) nboxes, nlevels, ifpgh, ifnewtree
c        write(*,*) ntarg,  ifpghtarg
c
c
c      do i=1,nboxes
c        do k=1,npbox
c          fright(1,k,i)=0.0d0
c          fright(2,k,i)=0.0d0
c        enddo
c      enddo
        t1=omp_get_wtime()
        call boxfgt(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2       fright,ifpgh,vpotinc,grad,hess,ifnewtree,ntarg,targs,
     3       ifpghtarg,vpotpinc,gradp,hessp,timeinfo)
        t2=omp_get_wtime()
c
        write(*,*)'time for fgt=',jj,t2-t1
        write(*,*)'speed for fgt=',npbox*nleafbox/(t2-t1)
c
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do j=1, npbox
               vpot4(1,j,i,jj)=vpotinc(1,j,i)
               vpot4(2,j,i,jj)=vpotinc(2,j,i)
               vpotinc(1,j,i)=vpotinc(1,j,i)*scal*twhts(jj)/(pi*delta)
               vpotinc(2,j,i)=vpotinc(2,j,i)*scal*twhts(jj)/(pi*delta)
               vpot(1,j,i)=vpot(1,j,i)+vpotinc(1,j,i)
               vpot(2,j,i)=vpot(2,j,i)+vpotinc(2,j,i)
c               write(205,*)jj,i,j,vpot(1,j,i),vpot(2,j,i)
            enddo
          endif
        enddo
c
c     end of the jj loop
      enddo

c     compute initial potential
c
        do i=1, nboxes
          do k=1, npbox
            initpot(1,k,i)=0.0d0
            initpot(2,k,i)=0.0d0
          enddo
        enddo

        delta = 4.0d0*dt*visc
        ifnewtree = 0
c        write(*,*) 'calling boxfgt'
        t1=omp_get_wtime()
        call boxfgt(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2       u3,ifpgh,initpot,grad,hess,ifnewtree,ntarg,targs,
     3       ifpghtarg,vpotpinc,gradp,hessp,timeinfo)
        t2=omp_get_wtime()
c
        write(*,*)'time for fgt=',t2-t1
        write(*,*)'speed for fgt=',npbox*nleafbox/(t2-t1)
c
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do j=1, npbox
               vpot(1,j,i)=vpot(1,j,i)+initpot(1,j,i)/(pi*delta)
               vpot(2,j,i)=vpot(2,j,i)+initpot(2,j,i)/(pi*delta)
            enddo
          endif
        enddo
c
c    correction process
c
      if(nordert .eq. 2) then
        nnodes=2
        tnodes(1)=0.0d0
        tnodes(2)=dt
        twhts(1)=1.0d0/2.0d0
        twhts(2)=1.0d0/2.0d0
      elseif(nordert .eq. 4) then
        nnodes=4
        tnodes(1)=0.0d0
        tnodes(2)=dt
        tnodes(3)=2*dt
        tnodes(4)=3*dt
        twhts(1)= 1.0d0/24.0d0
        twhts(2)=-5.0d0/24.0d0
        twhts(3)=19.0d0/24.0d0
        twhts(4)= 9.0d0/24.0d0
      else
        nordert=2
        nnodes=2
        tnodes(1)=0.0d0
        tnodes(2)=dt
        twhts(1)=1.0d0/2.0d0
        twhts(2)=1.0d0/2.0d0
        write(*,*) 'volpot: only nordert=2,3,4 allowed!'
        write(*,*) 'nordert reset to 2!'
      endif
c
c  times of correction process
      ncorrection=1

      do ll=1,ncorrection

      do i=1,nboxes
        do k=1,npbox
           vpot_correction(1,k,i)=0.0d0
           vpot_correction(2,k,i)=0.0d0
        enddo
      enddo



      scal=dt
      do jj=1,nnodes-1
       kk=jj+1
       delta=((nordert-1)*dt-tnodes(jj))*visc*4.0d0
c
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do j=1, npbox
               vpot_correction(1,j,i)=vpot_correction(1,j,i)+
     1              vpot4(1,j,i,kk)*scal*twhts(jj)/(pi*delta)
               vpot_correction(2,j,i)=vpot_correction(2,j,i)+
     1              vpot4(2,j,i,kk)*scal*twhts(jj)/(pi*delta)
            enddo
          endif
        enddo
c
c     end of the jj loop
      enddo
c
      tt=tpre+nordert*dt
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           vpot,eps,fnm4)

c
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do j=1, npbox
               vpot_correction(1,j,i)=vpot_correction(1,j,i)+
     1              fnm4(1,j,i)*scal*twhts(4)+
     2              initpot(1,j,i)/(pi*delta)
               vpot_correction(2,j,i)=vpot_correction(2,j,i)+
     1              fnm4(2,j,i)*scal*twhts(4)+
     2              initpot(2,j,i)/(pi*delta)
            enddo
          endif
        enddo
c
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do j=1, npbox
               vpot(1,j,i) = vpot_correction(1,j,i)
               vpot(2,j,i) = vpot_correction(2,j,i)
            enddo
          endif
        enddo

c   end of ll
        enddo

c

c
      if(nordert .eq. 4) then
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do k=1, npbox
              u0(1,k,i)=u1(1,k,i)
              u0(2,k,i)=u1(2,k,i)
              u1(1,k,i)=u2(1,k,i)
              u1(2,k,i)=u2(2,k,i)
              u2(1,k,i)=u3(1,k,i)
              u2(2,k,i)=u3(2,k,i)
              u3(1,k,i)=vpot(1,k,i)
              u3(2,k,i)=vpot(2,k,i)
            enddo
          endif
        enddo
      endif
c
      if(nordert .eq. 2) then
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do k=1, npbox
              u0(1,k,i)=u1(1,k,i)
              u0(2,k,i)=u1(2,k,i)
              u1(1,k,i)=vpot(1,k,i)
              u1(2,k,i)=vpot(2,k,i)
            enddo
          endif
        enddo
      endif

c
      iadap=1
      if(iadap .gt. 0 ) then
          write(*,*)'calling vol_tree_adap_fv '
          dpars(1)=tt+dt
          epstree=1.0d-1*eps
c          epstree=eps
          write(*,*)'epstree=',epstree
c          write(*,*)rintl(0),rintl(1),rintl(2),rintl(3)
c          write(*,*)rintl(4),rintl(5),rintl(6),rintl(7)
          t1=omp_get_wtime()
          call vol_tree_adap_fv4nsnew(ndim,epstree,ipoly,
     1         iptype,norder,npbox,nboxes,mxboxes,nlevels,
     2         mxlevels,ltree,itree,iptr,centers,boxsize,
     3         rintl,dpars,zpars,ipars,iperiod,fforce,
     4         u0,u1,u2,u3)
          t2=omp_get_wtime()
c
      nleafbox = 0
      
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
      enddo
c
        write(*,*)'time for adap=',t2-t1
        write(*,*)'speed for adap=',npbox*nleafbox/(t2-t1)
      endif
c
c
      end subroutine
