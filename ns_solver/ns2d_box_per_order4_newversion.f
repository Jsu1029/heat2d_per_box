c u0: the solution at t_1
c u1: the solution at t_1

      subroutine ns2dperfm4(nd,norder,nordert, finit,
     1           fforce,dt, ntot, eps, mxltree, mxboxes,
     2           mxlevels,ltree, itree, iptr, centers,
     3           nlevels,boxsize, nboxes, ifnewtree,velocity,vorticity,
     4           ntarg,velocityp,visc)
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
      real *8,allocatable :: un(:,:,:)
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
       dpars(1)=dt
      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1     norder,iptype,eta,fforce,nd,dpars,zpars,ipars,
     2     ifnewtree,mboxes,mlevels,ltree,rintl)

      
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
      allocate(un(nd,norder**ndim,mxboxes))
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
     1 centers,nlevels,mlevels,boxsize, nboxes,mboxes,fforce,
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
     4           u0,u1,u2,un,visc)



c   test init
c
      ipoly=1
      dpars(1)=dt
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
               sum_err=sum_err+abs(u1(1,j,ibox)-
     1                   velocity_ex(1,j,ibox))**2
c             write(*,*)u0(2,j,ibox),u1(2,j,ibox),u2(2,j,ibox)
             enddo
          endif
        enddo
      enddo
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'usol1 err2_rel=',err2_rel
c
      ipoly=1
      dpars(1)=dt*2
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
               sum_err=sum_err+abs(u2(1,j,ibox)-
     1                   velocity_ex(1,j,ibox))**2
             enddo
          endif
        enddo
      enddo
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'usol2 err2_rel=',err2_rel
c
      ipoly=1
      dpars(1)=dt*3
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
               sum_err=sum_err+abs(un(1,j,ibox)-
     1                   velocity_ex(1,j,ibox))**2
             enddo
          endif
        enddo
      enddo
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'usol3 err2_rel=',err2_rel


c
      iadap=1
      call ns2d_am24_main(nd,norder, nordert, fforce,
     1     dt, ntot, eps, mxltree,mboxes,mlevels,rintl,
     2     ltree, itree, iptr, centers, nlevels,
     3     boxsize, nboxes, ifnewtree,iadap,ntarg,targs,
     4     u0, u1,u2,un,visc)

c
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               velocity(1,j,ibox)=un(1,j,ibox)
               velocity(2,j,ibox)=un(2,j,ibox)
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
        write(*,*)'*****************************'
c
        t3=second()
          dpars(1)=tt
          t1=second()
          call vol_tree_adap_fv(nd,ndim,epstree,ipoly,
     1         iptype,norder,npbox,nboxes,mboxes,nlevels,
     2         mlevels,ltree,itree,iptr,centers,boxsize,
     3         rintl,dpars,zpars,ipars,iperiod,fforce,fright)
          t2=second()
c
      nleafbox = 0
      
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
      enddo
c
        write(*,*)'******'
        write(*,*)'time for adap_fv=',t2-t1
        write(*,*)'speed for adap_fv',npbox*nleafbox/(t2-t1)

         write(*,*) 'nboxes, mboxes=', nboxes, mboxes
c
c
        write(*,*)'******'
c      t1=second()
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tpre,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           fright,eps,frightS1)
c
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           fright,eps,frightS2)
c      t2=second()
c
c        write(*,*)'time for compute fs=',t2-t1
c        write(*,*)'speed for compute fs',npbox*nleafbox/(t2-t1)

c    compute volume potential
c
c      write(*,*)'call volpot1'
c
        write(*,*)'******'
      t1=second()
      call volpot1(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes,
     2           dt, iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           frightS1,frightS2,eps, vpot, vpotp,visc)
       t2=second()

        write(*,*)'time for compute vol1=',t2-t1
        write(*,*)'speed for compute vol1=',npbox*nleafbox/(t2-t1)

c      write(*,*)'call compute_initpot'
c    compute initial potential
c
        write(*,*)'******'
       t1=second()
      call compute_initpot(norder, ltree, itree, iptr, centers,
     1  nlevels,mxlevels, boxsize, nboxes,mxboxes, fright, dt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           eps, initpot, initpotp,visc)
        t2=second()

        write(*,*)'time for compute Init=',t2-t1
        write(*,*)'speed for compute Init=',npbox*nleafbox/(t2-t1)
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

      t4=second()
      write(*,*)'******'
      write(*,*)'time for one step of init_tree',t4-t3
      write(*,*)'speed for one step of init_tree',
     1           npbox*nleafbox/(t4-t3)
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
      do i = 1, nboxes
         do j = 1, npbox
            do k = 1, nd
               frightS(k,j,i) = 0.0d0
            enddo
         enddo
      enddo
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
c      t1=second()
      call evalgh(ndim,nd,ipoly,nlevels,ltree,itree,iptr,
     1       boxsize,norder,nboxes,initvalue,grad_init,hess_init)

      dpars(1)=tt
c
      call evalf(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,
     2       dpars,zpars,ipars,fforce,fright)
c
c      t2=second()
c      write(*,*)'time for eval',t2-t1
c      write(*,*)'speed for eval',npbox*nleafbox/(t2-t1)

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
c      t1=second()
      call bdmk_old_per(nd,ndim,eps,ikernel,beta,ipoly,norder,
     1    npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2    convective,ifpghdmk,pot,grad,hess,ntarg,targs,
     3    ifpghtargdmk,potp,gradp,hessp,tottimeinfo)
c      t2=second()
c      write(*,*)'time for bdmk_old_per =',t2-t1
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
      nleafbox = 0
      
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
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

        t1=second()
        call boxfgt(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2       fright1,ifpgh,vpotinc,grad,hess,ifnewtree,ntarg,targs,
     3       ifpghtarg,vpotpinc,gradp,hessp,timeinfo)
        t2=second()
c
        write(*,*)'time for fgt=',t2-t1
        write(*,*)'speed for fgt=',npbox*nleafbox/(t2-t1)
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
c      if(ifptarg .gt. 0) then
c
c      call evalt(ndim,nd,ipoly,norder,nboxes,nlevels,
c     1       ltree,itree,iptr,centers,boxsize,vpot,
c     2       ntarg,targs,vpotp)
c      endif
 
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
      nleafbox = 0
      
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
      enddo
c
      t1=second()
      call boxfgt(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fright,
     2    ifpgh,pot,grad,hess,ifnewtree,ntarg,targs,
     3    ifpghtarg,potp,gradp,hessp,timeinfo)
      t2=second()
c
        write(*,*)'time for fgt=',t2-t1
        write(*,*)'speed for fgt',npbox*nleafbox/(t2-t1)
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
c     output:
c     u1: the solution at t_1, only needed if nordert ==4
c     u2: the solution at t_2, only needed if nordert ==4
c     un: the solution at t_{nordert-1}
      subroutine ns2d_init_usol(nd,norder,nordert,
     1           fforce,dt0, eps, mxltree, mxboxes,
     2           mxlevels,rintl,ltree, itree, iptr, centers,
     3           nlevels,boxsize, nboxes, ifnewtree,ntarg,
     4           targs,u0,u1,u2,un,visc)
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
      real*8 un(nd,norder**ndim,mxboxes)
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
      nleafbox = 0
      
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
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
      t1=second()
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
          un(k,i,ibox)=usave(k,i,ibox,3)
c          write(*,*)k,i,ibox,un(k,i,ibox)
          enddo
        enddo
        enddo
      else if(nordert .eq. 2) then
        do ibox=1, nboxes
        do i=1, npbox
          do k =1,nd
          u1(k,i,ibox)=0.0d0
          u2(k,i,ibox)=0.0d0
          un(k,i,ibox)=usave(k,i,ibox,1)
          enddo
        enddo
        enddo
      endif
c     like this?
c     yeah, like this
c
      t2=second()
      write(*,*)'******'
      write(*,*)'time for Richardson on step',t2-t1
      write(*,*)'speed for Richardson on step',
     1     npbox*nleafbox/(t2-t1)
      end subroutine

c

c
      subroutine ns2d_am24_main(nd,norder, nordert, fforce,
     1     dt, ntot, eps, mxltree,mxboxes,mxlevels,rintl,
     2     ltree, itree, iptr, centers, nlevels,
     3     boxsize, nboxes, ifnewtree,iadap,ntarg,targs,
     4     u0, u1,u2,un,visc)
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
      real*8 un(nd,norder**ndim,mxboxes)
      real*8 rintl(0:200), dpars(1000)
      complex*16 zpars(10), zk
c     allocatables
c     the tree structure
      real *8, allocatable :: xref(:,:),wts(:)
      real*8, allocatable:: fnm1(:,:,:),fnm2(:,:,:)
      real*8, allocatable:: fnm3(:,:,:)
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
      open(21,file='xf3.dat')
      open(22,file='yf3.dat')
      open(190,file='nleaves3.dat')
      open(26,file='velocity3.dat')
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
c      epstree=1.0d-2*eps
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
      allocate(fnm1(nd,npbox,mxboxes))
      allocate(fnm2(nd,npbox,mxboxes))
      allocate(fnm3(nd,npbox,mxboxes))

c     initialize fnm1, fnm2, fnm3
      if(nordert .eq. 2) then
c       2nd order:
c
      tt=0.0d0
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           u0,eps,fnm1)

c
      elseif(nordert .eq. 4) then
c       4th order:
c
      tt=0.0d0
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           u0,eps,fnm3)
c
      tt=dt
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           u1,eps,fnm2)
c
      tt=2.0d0*dt
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           u2,eps,fnm1)
      endif

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
        t1=second()
        call vpot_ns_am24(nd,norder, nordert, ltree, itree,
     1           iptr, centers,mxltree, nlevels,mxlevels, boxsize,
     2         nboxes,mxboxes,fforce, tpre, dt, iperiod, ntarg, targs,
     3     rintl,ifvtarg, ifptarg,ifnewtree,iadap, eps,visc,
     4     fnm3,fnm2,fnm1,un)
        t2=second()
c
      nleafbox = 0
      
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
      enddo
c
        write(*,*)'time for one step=',t2-t1
        write(*,*)'speed for one step=',npbox*nleafbox/(t2-t1)
      


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


      sum_err=0.0d0
      sum_u=0.0d0
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               sum_u=sum_u+velocity_ex(1,j,ibox)**2
               sum_err=sum_err+abs(un(1,j,ibox)-
     1                   velocity_ex(1,j,ibox))**2
             enddo
          endif
        enddo
      enddo
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'usol err2_rel=',err2_rel,sum_err,sum_u

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
     4    fnm3,fnm2,fnm1,un)
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
      real*8 fnm3(nd,norder**ndim,mxboxes)
      real*8 fnm2(nd,norder**ndim,mxboxes)
      real*8 fnm1(nd,norder**ndim,mxboxes)
      real*8 un(nd,norder**ndim,mxboxes)
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
      real*8, allocatable:: fn(:,:,:)
      real*8, allocatable:: fnm(:,:,:)
      real*8, allocatable:: fright(:,:,:)
c
      real*8, allocatable:: pot(:,:,:), vpot(:,:,:)
      real*8, allocatable:: potp(:,:), vpotp(:,:)
      real *8, allocatable:: velocity_ex(:,:,:)
c
      real*8, allocatable:: vpot4(:,:,:,:)
      real*8, allocatable:: vpot_correction(:,:,:)
      real*8, allocatable:: initpot(:,:,:)
      real*8, allocatable:: g(:,:,:)
      external fforce,funs
c      open(91,file='test_fnm3.txt')
c
      pi=3.1415926535897932d0
      npbox=norder**ndim
      nhess = ndim*(ndim+1)/2
      iptype=2
      ipoly=1
c
c
      t1=second()
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

      ifpgh=1
      ifpghtarg=0
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
      allocate(fn(nd,npbox,mxboxes))
      allocate(fnm(nd,npbox,mxboxes))
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
      allocate(g(nd,npbox,mxboxes))
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
              fn(k,j,i)=0.0d0
              unm1(k,j,i)=un(k,j,i)
              fright(k,j,i)=0.0d0
            enddo
            write(101,*)i,j,fnm1(1,j,i),fnm1(2,j,i)
            write(102,*)i,j,fnm2(1,j,i),fnm2(2,j,i)
            write(103,*)i,j,fnm3(1,j,i),fnm3(2,j,i)
         enddo
      enddo
c
      nleafbox = 0
      
      do i=1,nboxes
        if(itree(iptr(4)+i-1).eq.0) nleafbox = nleafbox+1
      enddo
      t2=second()
      write(*,*)'time for set=',t2-t1

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
      if(jj .eq. nnodes-3) then
c
          do i=1,nboxes
            do k=1,npbox
              do j=1,nd
              fright(j,k,i)=fnm3(j,k,i)
              enddo
            enddo
          enddo

c
      elseif(jj .eq. nnodes-2) then
c
          do i=1,nboxes
            do k=1,npbox
              do j=1,nd
              fright(j,k,i)=fnm2(j,k,i)
              enddo
c            write(91,*)i,k,fnm2(1,k,i),fnm2(2,k,i)
            enddo
          enddo
c
      elseif(jj .eq. nnodes-1) then
c
          do i=1,nboxes
            do k=1,npbox
              do j=1,nd
              fright(j,k,i)=fnm1(j,k,i)
c            write(91,*)i,k,fnm1(1,k,i),fnm1(2,k,i)
              enddo
            enddo
          enddo
c
      elseif(jj .eq. nnodes) then
c
      t1=second()
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           un,eps,fn)
      t2=second()
      write(*,*)'time for fs=',t2-t1
          do i=1,nboxes
            do k=1,npbox
              do j=1,nd
              fright(j,k,i)=fn(j,k,i)
              enddo
              write(104,*)i,k,fn(1,k,i),fn(2,k,i)
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
c
c        do i=1, ntarg
c          do k=1, nd
c            vpotpinc(k,i)=0.0d0
c            vpotpinc(k,i)=0.0d0
c          enddo
c        enddo

        ifnewtree = 0

c
        t1=second()
        call boxfgt(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2       fright,ifpgh,vpotinc,grad,hess,ifnewtree,ntarg,targs,
     3       ifpghtarg,vpotpinc,gradp,hessp,timeinfo)
        t2=second()
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
               write(105,*)jj,i,j,vpot(1,j,i),vpot(2,j,i)
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
        t1=second()
        call boxfgt(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2       un,ifpgh,initpot,grad,hess,ifnewtree,ntarg,targs,
     3       ifpghtarg,vpotpinc,gradp,hessp,timeinfo)
        t2=second()
c
        write(*,*)'time for fgt=',t2-t1
        write(*,*)'speed for fgt=',npbox*nleafbox/(t2-t1)
c
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do j=1, npbox
               vpot(1,j,i)=vpot(1,j,i)+initpot(1,j,i)/(pi*delta)
               vpot(2,j,i)=vpot(2,j,i)+initpot(2,j,i)/(pi*delta)
c            write(91,*)j,i,un(2,j,i),initpot(1,j,i),initpot(2,j,i)
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
c
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
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do j=1, npbox
               g(1,j,i) = vpot_correction(1,j,i)+
     1                     initpot(1,j,i)/(pi*delta)
               g(2,j,i) = vpot_correction(2,j,i)+
     2                     initpot(1,j,i)/(pi*delta)
            enddo
          endif
        enddo
c
      tt=tpre+nordert*dt
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           vpot,eps,fnm)
       
        jj=nnodes
c
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do j=1, npbox
               vpot_correction(1,j,i)=vpot_correction(1,j,i)+
     1              fnm(1,j,i)*scal*twhts(jj)
               vpot_correction(2,j,i)=vpot_correction(2,j,i)+
     1              fnm(2,j,i)*scal*twhts(jj)
            enddo
          endif
        enddo
c
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do j=1, npbox
               un(1,j,i) = vpot_correction(1,j,i)
               un(2,j,i) = vpot_correction(2,j,i)
            enddo
          endif
        enddo

c   end of ll
        enddo

c

c
      if(nordert .eq. 4) then
        do i=1, nboxes
c          if(itree(iptr(4)+i-1).eq.0) then
            do k=1, npbox
              do j=1,nd
              fnm3(j,k,i)=fnm2(j,k,i)
              fnm2(j,k,i)=fnm1(j,k,i)
              fnm1(j,k,i)=fnm(j,k,i)
              enddo
            enddo
c          endif
        enddo
      endif
c
      if(nordert .eq. 2) then
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do k=1, npbox
               do j=1,nd
                  fnm1(j,k,i)=fnm(j,k,i)
               enddo
            enddo
          endif
        enddo
      endif

c
      iadap=0
      if(iadap .gt. 0 ) then
          write(*,*)'calling vol_tree_adap_fv '
          dpars(1)=visc
          tt=tpre+nordert*dt
c          dpars(1)=tt
          epstree=1.0d0*eps
c          epstree=eps*10
          write(*,*)'epstree=',epstree
          t1=second()

          call am24_adap_tree_nd(nd,norder, ltree, itree,
     1     iptr, rintl, centers, nlevels, mxlevels, boxsize,
     2     nboxes, mxboxes, fforce, dpars, zpars, ipars,
     3     tt, dt, iperiod, eps, ntarg, targs, nordert,
     4     g,fnm3,fnm2,fnm1,un,initpot,iadap)
        
          t2=second() 
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


c
      subroutine am24_adap_tree_nd(nd, norder, ltree, itree,
     1     iptr, rintl, centers, nlevels, mlevels, boxsize,
     2     nboxes, mboxes,fforce, dpars, zpars, ipars,
     3     tt, dt, iperiod, eps, ntarg, targs, nordert,
     4     vpot, fnm3, fnm2, fnm1, un, unm1,iadap)
      implicit real*8 (a-h,o-z)
c     global vars
      integer norder, ltree, nboxes, ndim, nlevels
      integer iadap
      integer mboxes, mlevels, nd, ntarg, nordert
      integer iperiod, ipars(100)
      parameter(ndim=2)
      integer itree(ltree), iptr(12)
      real*8 centers(ndim, mboxes), boxsize(0:mlevels)
      real*8 tt, dt, eps, alpha, f0(1)
      real*8 vpot(nd,norder**ndim,mboxes)
      real*8 unm1(nd,norder**ndim,mboxes)
      real*8 fnm1(nd,norder**ndim,mboxes)
      real*8 fnm2(nd,norder**ndim,mboxes)
      real*8 fnm3(nd,norder**ndim,mboxes)
      real*8 un(nd,norder**ndim,mboxes)
      real*8 dpars(1000)
      complex *16 zk, zpars(10)
c     local vars
c     flags for refinement and coarsening
      integer isgn(ndim,2**ndim)
      integer, allocatable:: nblock(:)
      integer, allocatable:: ilevstart(:)
      integer, allocatable:: nboxid(:), ifdelete(:)
      integer, allocatable:: ifrefine(:)
      integer, allocatable:: irefinelev(:), imaxrefinelev(:)
c     real allocatables
c     need more of these fvals thing?
      real*8, allocatable :: fn(:,:,:)
      real*8, allocatable :: uallvals(:,:,:), uallcoefs(:,:,:)
      real*8, allocatable :: fcoefs(:,:,:), ucoefs(:,:,:)
      real*8, allocatable :: fcoefs1(:,:,:), fcoefs2(:,:,:)
      real*8, allocatable :: fcoefs3(:,:,:)
      real*8, allocatable :: ucoefs1(:,:,:), vcoefs(:,:,:)
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8, allocatable :: rmask(:), umat_nd(:,:,:)
      real*8, allocatable:: tnodes(:), twhts(:)
      real *8 rintl(0:mlevels)
      real *8 xy(ndim)
      real*8 targs(ndim,ntarg)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
c     change the following names
      real *8, allocatable:: polyvc(:,:,:,:)
c     nd-vecs
      real*8, allocatable:: u0(:), u1(:), f1(:), g(:)
c     complex allocatables
      external fevaln, devaln
c
c     do down the tree
c     check resolution of un and f(tt,un)
c     on each leaf node
c     subdivide if needed
c     interpolate vpot, unm1 (as the initial guess)
c     solve for un again
c
c-------------------------------------------------------
c
c      write(*,*) norder, ltree
c      write(*,*) nlevels, mlevels, nboxes, mboxes
c      write(*,*) tt, iperiod, eps, alpha
c
c
      allocate(tnodes(nordert))
      allocate(twhts(nordert))
      mc=2**ndim
      mnbors=3**ndim
      npols = norder**ndim
      npbox = norder**ndim
      ipoly=1
      iptype=2
      pi=3.1415926535897932d0
      visc=dpars(1)
      delta = 4.0d0*dt*visc
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
c      write(*,*) npols, npbox, ndim, nd,tt, iperiod, eps
c      write(*,*)ltree,nlevels, mlevels,nboxes, mboxes
c
      allocate(fn(nd,npbox,mboxes))
c     allocate coeffs for them all:
c     un, fn, unm1, vpot
      allocate(ucoefs(nd,npbox,mboxes))
      allocate(fcoefs(nd,npbox,mboxes))
      allocate(fcoefs1(nd,npbox,mboxes))
      allocate(fcoefs2(nd,npbox,mboxes))
      allocate(fcoefs3(nd,npbox,mboxes))
      allocate(ucoefs1(nd,npbox,mboxes))
      allocate(vcoefs(nd,npbox,mboxes))
c
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
      allocate(umat_nd(norder,norder,ndim))
      allocate(polyvc(norder,norder,ndim,mc))
c
      allocate(u0(nd))
      allocate(u1(nd))
      allocate(f1(nd))
      allocate(g(nd))
c
      call get_child_box_sign(ndim,isgn)
      call get_p2c_interp_matrices(ndim,ipoly,norder,isgn,polyvc)
      allocate(grid(ndim,npbox))
c
c------------------------------------------------------
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
c     eval fun on the grid of leaf nodes -> fvals
c
      ifvtarg=1
      ifptarg=0
c      t1=second()
      call compute_fs(norder, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes, fforce,tt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           un,eps,fn)
c      t2=second()
c      write(*,*)'time for fs=',t2-t1
c
c     compute rmask and rsum for estimating whether the function is resolved
      allocate(rmask(npbox))
c
      call tens_prod_get_rmask(ndim,iptype,norder,npbox,
     1    rmask,rsum)
c
      if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**ndim)
      if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**ndim)
      if(iptype.eq.0) rsc = 1.0d0
c
      reps = eps*rsc
      nboxes0 = nboxes
      nlevels0 = nlevels
c
      allocate(ifrefine(mboxes))
c
c     mboxes here means the max of nboxes
 
      do ibox = 1, mboxes
        ifrefine(ibox)=0
      enddo


c
c--------------------------------------------------
c     convert potential to its polynomial expansion coefficients
c     (for all the leaf nodes)
c
c
c     be careful: which to which
c     convert them all to coeffs, I believe
      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,fn,fcoefs,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,fnm1,fcoefs1,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,fnm2,fcoefs2,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,fnm3,fcoefs3,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,un,ucoefs,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,unm1,ucoefs1,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,vpot,vcoefs,umat_nd)
c
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
             call fun_err(nd,npbox,fcoefs(1,1,ibox),
     1            rmask,iptype,rscale,erra)
             call fun_err(nd,npbox,ucoefs(1,1,ibox),
     1            rmask,iptype,rscale,errb)
             erra = erra/rsum
             errb = errb/rsum
             erri = max(erra, errb)
c             write(*,*)ibox,erra,errb
             if(erri.gt.reps*rintl(ilev)) then
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
      t1=second()
      nnewboxes=0
c      write(*,*)'nnewboxes=',nnewboxes
      if(iadap .gt. 0) then
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
c-------------------
c           modification here:
c           solve for u again on the new box
c           calling the nonlinear solver
c           might need to interpolate vpot to the new box first
c           but the assumption is vpot is resolved
c           interpolate vpot, fn, unm1
c           and use unm1 as initial value to solve
c
c           interpolate vpot from ibox to jbox, which is
c           the j-th child of ibox
c
c           interpolate unm1 too
c           as an initial guess for the nonlinear solve
            call ortho_eval_nd(ndim,nd,norder,
     1           ucoefs1(1,1,ibox),norder,unm1(1,1,jbox),
     2           polyvc(1,1,1,j))
c
            call ortho_trans_nd(ndim,nd,0,norder,unm1(1,1,jbox),
     1           ucoefs1(1,1,jbox),umat_nd)
c           interpolate vpot from ibox to jbox, which is
c           the j-th child of ibox
            call ortho_eval_nd(ndim,nd,norder,
     1           vcoefs(1,1,ibox),norder,vpot(1,1,jbox),
     2           polyvc(1,1,1,j))
c
            call ortho_trans_nd(ndim,nd,0,norder,vpot(1,1,jbox),
     1           vcoefs(1,1,jbox),umat_nd)
c           interpolate fn too
c           as an initial guess for the nonlinear solve
            call ortho_eval_nd(ndim,nd,norder,
     1           fcoefs(1,1,ibox),norder,fn(1,1,jbox),
     2           polyvc(1,1,1,j))
c
            call ortho_trans_nd(ndim,nd,0,norder,fn(1,1,jbox),
     1           fcoefs(1,1,jbox),umat_nd)
c
c           interpolate fnm1 too
c           for subsequent steps
            call ortho_eval_nd(ndim,nd,norder,
     1           fcoefs1(1,1,ibox),norder,fnm1(1,1,jbox),
     2           polyvc(1,1,1,j))
c
            call ortho_trans_nd(ndim,nd,0,norder,fnm1(1,1,jbox),
     1           fcoefs1(1,1,jbox),umat_nd)
c
c           interpolate fnm2 too
c           for subsequent steps
            call ortho_eval_nd(ndim,nd,norder,
     1           fcoefs2(1,1,ibox),norder,fnm2(1,1,jbox),
     2           polyvc(1,1,1,j))
c
            call ortho_trans_nd(ndim,nd,0,norder,fnm2(1,1,jbox),
     1           fcoefs2(1,1,jbox),umat_nd)
c
c           interpolate fnm3 too
c           for subsequent steps
            call ortho_eval_nd(ndim,nd,norder,
     1           fcoefs3(1,1,ibox),norder,fnm3(1,1,jbox),
     2           polyvc(1,1,1,j))
c
            call ortho_trans_nd(ndim,nd,0,norder,fnm3(1,1,jbox),
     1           fcoefs3(1,1,jbox),umat_nd)

c             I suppose the indices are correct?
              do ii=1,npbox
              do kk=1, nd
                un(kk,ii,jbox)=vpot(kk,ii,jbox)+unm1(kk,ii,jbox)/
     1     (pi*delta)+fn(kk,ii,jbox)*dt*twhts(nordert)
              enddo
              enddo

c
c-------------------
c           now instead of interpolating un and fn
c           we have the function values re-evaluated
c           convert them to coefficients and estimate errors

c
            call ortho_trans_nd(ndim,nd,0,norder,un(1,1,jbox),
     1           ucoefs(1,1,jbox),umat_nd)

            call fun_err(nd,npbox,ucoefs(1,1,jbox),rmask,
     1           iptype,rscale,erra)
c
            erra = erra/rsum
c
            if(erra .gt. reps*rintl(jlev)) then
c             flag the new box
              ifrefine(jbox) = 1
            endif

          enddo
c       end refining box ibox
        endif
      enddo
c     end of the iadap if
      endif
      t2=second()
c      write(*,*)'time for refine=',t2-t1
c
c      stop
c------------------------------------------------
c
      write(*,*) 'after refinement, nlevels=',nlevels
      write(*,*) 'nnewboxes=', nnewboxes
c      stop
c
      t1=second()
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
c      write(*,*) 'nnewlevels=',nnewlevels
c
c     attention: what if nnewboxes = 0?
c     move the allocation and stuff upward?
c      nd0 = 2
      nd0 = 6*nd
      allocate(uallvals(nd0,npbox,mboxes))
      allocate(uallcoefs(nd0,npbox,mboxes))
      allocate(grad(nd0,ndim,npbox,1))
      allocate(hess(nd0,ndim*(ndim+1)/2,npbox,1))
c
c      allocate(grad(nd,ndim,npbox,1))
c      allocate(hess(nd,ndim*(ndim+1)/2,npbox,1))
c
c     do we need fn? not sure. do it anyway.
c     yes we do need fn
c     we need to make sure both fn and un are resolved
      do ibox = 1, nboxes
        do i=1, npbox
          do kk=1, nd
            uallvals(kk,i,ibox) = un(kk,i,ibox)
            uallvals(2+kk,i,ibox) = unm1(kk,i,ibox)
            uallvals(4+kk,i,ibox) = fnm1(kk,i,ibox)
            uallvals(6+kk,i,ibox) = fnm2(kk,i,ibox)
            uallvals(8+kk,i,ibox) = fnm3(kk,i,ibox)
            uallvals(10+kk,i,ibox) = vpot(kk,i,ibox)
c
          uallcoefs(kk,i,ibox) = ucoefs(kk,i,ibox)
          uallcoefs(2+kk,i,ibox) = ucoefs1(kk,i,ibox)
          uallcoefs(4+kk,i,ibox) = fcoefs1(kk,i,ibox)
          uallcoefs(6+kk,i,ibox) = fcoefs2(kk,i,ibox)
          uallcoefs(8+kk,i,ibox) = fcoefs3(kk,i,ibox)
          uallcoefs(10+kk,i,ibox) = vcoefs(kk,i,ibox)
          enddo
        enddo
      enddo
c
      if (nnewboxes.eq.0) goto 5600
c
      ifpgh = 1
c     after refinement, reorg un and fn together
c
c

c      stop
c
c     reorg the tree and the packed arrays: uallvals, uallcoefs
      call vol_tree_reorg_laddr(ndim, mboxes, nd0, npbox,
     1    nblock, nboxid, nnewboxes, nboxes0, mlevels,
     2    nnewlevels, nlevels0, centers, itree(iptr(1)),
     3    itree(iptr(2)), itree(iptr(3)), itree(iptr(4)),
     4    itree(iptr(5)), ifpgh, uallvals, uallcoefs, grad, hess)
c
5600  continue
c
      write(*,*) 'after reorg, nlevels=', nlevels
      write(*,*) 'nboxes=', nboxes
      t2=second()
c      write(*,*)'time for reorg=',t2-t1
c      stop
c
c------------------------------------------------
c
      t1=second()
      if(iadap .eq. 1) then
c     when iadap = 1, do the coarsening
c
      itype=0
c      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
c     1     iptr,boxsize,norder,uallvals,uallcoefs,umat_nd)
c
c     nd v.s. nd0? I believe it's nd0
c      write(*,*)'koko'
      call treedata_trans_nd(ndim,nd0,itype,nlevels,itree,
     1     iptr,boxsize,norder,uallvals,uallcoefs,umat_nd)
c     no need to call this routine?
c     uallcoefs already available?
c
      ifpgh = 1
      allocate(ifdelete(mboxes))
c
c     vol_tree_coarsen needs to check
c     the resolution of both un and fn
c     so...
c
c      write(*,*)'here??'
      call vol_tree_coarsen(nd0,ndim,reps,ipoly,norder,npbox,
     1     nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2     ifpgh,uallvals,uallcoefs,grad,hess,
     3     nblock,nboxid,ndelboxes,ifdelete)
c     I believe the coarsening routine does the
c     interpolation from children to parent
c     not the fault of vol_tree_coarsen as it seems
c
      write(*,*) 'ndelboxes=', ndelboxes
c
      if(ndelboxes .gt. 0) then
        call fgt_vol_tree_reorg_after_coarsen(ndim,nboxes,nd0,npbox,
     1       nblock,nboxid,ndelboxes,ifdelete,
     2       centers,nlevels,itree(iptr(1)),itree(iptr(2)),
     3       itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     4       ifpgh,uallvals,uallcoefs,grad,hess)

      endif
c
c      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
      nboxes = nboxes - ndelboxes
      write(*,*) 'after vol_tree_coarsen:'
      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
c     end of the iadap if
      endif
      t2=second()
c      write(*,*)'time for coarsen=',t2-t1
c
c------------------------------------------------
c
c     now compute the colls and make it
c     level-restricted again
c
      t1=second()
      call computecoll(ndim,nlevels,nboxes,itree(iptr(1)),
     1    boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2    itree(iptr(5)),iperiod,itree(iptr(6)),itree(iptr(7)))
c
      call vol_tree_fix_lr_interp(ndim,nd0,ipoly,iperiod,norder,npbox,
     1    ifpgh,uallvals,uallcoefs,grad,hess,
     2    mboxes,mlevels,centers,boxsize,nboxes,nlevels,
     3    itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),itree(iptr(6)),itree(iptr(7)))
c
      write(*,*) 'after fixing lr:'
      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
      t2=second()
c      write(*,*)'time for level-res=',t2-t1
c
c
c     now unpack the packed arrays
c
      do ibox = 1, nboxes
        do i=1, npbox
          do kk=1, nd
            un(kk,i,ibox) = uallvals(kk,i,ibox)
            unm1(kk,i,ibox) = uallvals(2+kk,i,ibox)
            fnm1(kk,i,ibox) = uallvals(4+kk,i,ibox)
            fnm2(kk,i,ibox) = uallvals(6+kk,i,ibox)
            fnm3(kk,i,ibox) = uallvals(8+kk,i,ibox)
            vpot(kk,i,ibox) = uallvals(10+kk,i,ibox)
          enddo
        enddo
      enddo

c
c------------------------------------------------
c
      end subroutine
