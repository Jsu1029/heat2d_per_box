c test fmm2dmk

c
      program test
      implicit real*8 (a-h,o-z)
c
      parameter(nd=2)
      parameter(ndim=2)
      parameter(norder=8)
      parameter(mxboxes=100000)
      parameter(mxltree=200000)
      parameter(mxlevels=30)
      parameter(ntarg=1000000)
      parameter(nx = 1000)
      parameter(ny = 1000)
      integer ipars(100), iptr(8)
      integer itree(mxltree)
      real*8 dt, eps
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      real*8 usol(norder**ndim,mxboxes)
      real*8 usolp(ntarg)
      real*8 uextp(ntarg)
      real*8 rintl(0:200), dpars(1000)
      integer nlevels

      real *8,allocatable :: fvals(:,:,:)
      real *8 targs(nd,ntarg)
      real*8 hx,hy
      complex*16 zpars(10), zk

      integer ikernel
      integer ifpgh,ifpghtarg
      real *8 beta

c
      real *8 timeinfo(100)

      real *8, allocatable :: pot(:,:,:)
      real *8, allocatable :: potp(:,:)

      real *8, allocatable :: potex(:,:,:)
      real *8, allocatable :: potpex(:,:)

      character *12 fname1
      character *8 fname2
      character *9 fname3
      
      real*8 errpe
      external fforce, uexact

      open(21,file='testns9_pot.data')
      open(22,file='fvals.data')
      open(26,file='potex.data')
      npbox=norder**ndim

c
      allocate(fvals(nd,npbox,mxboxes))
     
      allocate(pot(nd,npbox,mxboxes))
      allocate(potp(nd,ntarg))
c
      allocate(potex(nd,npbox,mxboxes))
      allocate(potpex(nd,ntarg))
    
      hx=1.0d0/(nx-1)
      hy=1.0d0/(ny-1)

      do i=1,nx
      do j=1,ny
        ii=(i-1)*ny+j
        targs(1,ii)=-0.5d0+(i-1)*hx
        targs(2,ii)=-0.5d0+(j-1)*hy
      enddo
      enddo

      done = 1
      pi = atan(done)*4
c
      eps=1.0d-9
      ipoly=1
      iperiod=1

c
      ifnewtree=0
c     when ifnewtree=1, prepare to refine for
c     one more level? I guess so.
c
      zk=1.0d0
      eta=1.0d0
c     not sure about the above two parameters
c
      epstree=1.0d-1*eps
      boxlen=1.0d0
      iptype=2
      npbox= norder**ndim

  
      ifpgh=1
      ifpghtarg=0
     
      dpars(1)=0.0d0
      
      write(*,*)eps,epstree
      call vol_tree_mem2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1     norder,iptype,eta,fforce,nd,dpars,zpars,ipars,
     2     ifnewtree,mboxes,mlevels,ltree,rintl)
c
      mboxes=mboxes*2
      mlevels=mlevels+5
      ltree=ltree*2
c
      write(*,*) 'mboxes=', mboxes
      write(*,*) 'mlevels=', mlevels
      write(*,*) 'ltree=', ltree

c
      call vol_tree_build2(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,uexact,nd,dpars,zpars,ipars,rintl,
     2    nboxes,mboxes,nlevels,mlevels,ltree,itree,iptr,centers,
     3    boxsize,fvals)
c
      write(*,*) ' '
      write(*,*) '******'
      write(*,*) 'after calling vol_tree_build2'
      write(*,*) 'nboxes =', nboxes
      write(*,*) 'mboxes =', mboxes
c
      write(*,*) 'nlevels=', nlevels
      write(*,*) 'mlevels=', mlevels
c
      write(*,*) 'ltree=', ltree
c
      write(*,*) ''

c
      dt=1.0d-5

      
      call compute_initpot(nd,norder, ltree, itree, iptr,
     1           centers,nlevels, boxsize, nboxes, fvals, dt,
     2           iperiod, ntarg, targs, ifpgh, ifpghtarg,
     3           eps, pot, potp)

c
      dpars(1)=dt
      call evalf(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,
     2       dpars,zpars,ipars,uexact,potex)

c
      sum_err=0.0d0
      sum_u=0.0d0
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               sum_u=sum_u+potex(1,j,ibox)**2
               sum_err=sum_err+abs(pot(1,j,ibox)
     1                -potex(1,j,ibox))**2
               write(21,*)pot(1,j,ibox),potex(1,j,ibox),
     1    abs(pot(1,j,ibox) -potex(1,j,ibox))
             enddo
          endif
        enddo
      enddo
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'err2_rel=',err2_rel
      write(*,*) 'solution and error written in fort.234'

      end program



c     as rhsfun in test_boxfgt*.f
      subroutine fforce(nd,xyz,dpars,zpars,
     1           ipars,u)
      implicit real *8 (a-h,o-z)
      integer nd, ndim, ipars(*)
      real *8 dpars(*), xyz(nd), t, u(nd)
      complex *16 zpars(*)
      real*8 s1, s2, res
c
      one=1.0d0
      pi=3.1415926535897932d0
c
      s1=2.0d0
      s2=2.0d0
c
      res=dsin(2*s1*pi*xyz(1))*dcos(2*s2*pi*xyz(2))
      u(1)=res
      u(2)=u(1)

      end subroutine

c     as rhsfun in test_boxfgt*.f
      subroutine uexact(nd,xyz,dpars,zpars,
     1           ipars,u)
      implicit real *8 (a-h,o-z)
      integer nd, ndim, ipars(*)
      real *8 dpars(*), xyz(*), t, u(nd)
      complex *16 zpars(*)
      real*8 s1, s2,x1,x2
c
      one=1.0d0
      pi=3.1415926535897932d0
c
      s1=2.0d0
      s2=2.0d0
c
      d=dpars(1)

      const1=1.0d0
      const2=4.0d0
c
      x1=xyz(1)
      x2=xyz(2)
c
      u(1)=const1*dexp(-(s1**2+s2**2)*pi**2*d*const2)*
     1              dsin(2*s1*pi*x1)*dcos(2*s2*pi*x2)

      u(2)=u(1)
c
      end subroutine


c**************************************************************
c     subroutine initpot
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
      subroutine compute_initpot(nd,norder, ltree, itree, iptr,
     1           centers,nlevels, boxsize, nboxes, fright, dt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           eps, pot, potp)
      implicit real*8 (a-h,o-z)
c     global vars
      integer norder, ltree, nboxes, ndim, nlevels
      integer iperiod, ntarg,nd
      integer ifvtarg, ifptarg
      parameter(ndim=2)
      integer itree(ltree), iptr(12)
      real*8 centers(ndim, nboxes), boxsize(0:nlevels)
      real*8 fright(nd,norder**ndim,nboxes), dt
      real*8 targs(ndim,ntarg), eps
c      real*8, allocatable:: pot(:,:,:)
c      real*8, allocatable:: potp(:,:)
      real*8 pot(nd,norder**ndim,nboxes), potp(nd,ntarg)
c     local vars
      real*8, allocatable:: timeinfo(:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
      real*8, allocatable:: gradp(:,:,:)
      real*8, allocatable:: hessp(:,:,:)
      open(31,file='fvlas.txt')
c
      pi=3.1415926535897932d0
      nhess = ndim*(ndim+1)/2
      npbox=norder**ndim
c
      allocate(timeinfo(100))
c
      
c      allocate(pot(nd,norder**ndim,nboxes))
c      allocate(potp(nd,ntarg))
      allocate(grad(nd,ndim,npbox,nboxes))
      allocate(hess(nd,nhess,npbox,nboxes))
      allocate(gradp(nd,ndim,ntarg))
      allocate(hessp(nd,nhess,ntarg))
c
c     allocate certain vars
c
c     simply call boxfgt and divide the
c     output by dt
c
      delta=4.0d0*dt
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
      call boxfgt(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fright,
     2    ifpgh,pot,grad,hess,ifnewtree,ntarg,targs,
     3    ifpghtarg,potp,gradp,hessp,timeinfo)

c
c---------------------------------------
c
c     divide the results by (4.0d0*pi*dt)
      if(ifvtarg .gt. 0) then
 
      do ilevel=0,nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               pot(1,j,ibox)=pot(1,j,ibox)/(4.0d0*pi*dt)
               pot(2,j,ibox)=pot(2,j,ibox)/(4.0d0*pi*dt)
               write(31,*)pot(1,j,ibox)
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
          potp(1,i)=potp(1,i)/(4.0d0*pi*dt)
          potp(2,i)=potp(2,i)/(4.0d0*pi*dt)
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

