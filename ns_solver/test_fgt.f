c test fmm2dmk

c
      program test
      implicit real*8 (a-h,o-z)
      integer norder, ndim, iprec,nx,ny
      integer ntot, mxltree, mxboxes, mxlevels,ntarg
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
c      real *8 fvals(nd,norder**ndim,mxboxes)
       real *8, allocatable :: fvals(:,:,:)
       real *8 targs(nd,ntarg)
     
      real*8, allocatable:: targ(:)
c      real*8, allocatable:: targs(:,:)
      real*8 hx,hy
      complex*16 zpars(10), zk

      integer ikernel
      integer ifpgh,ifpghtarg
      real *8 beta

c
      real *8 timeinfo(100)

      real *8, allocatable :: pot(:,:,:), potex(:,:,:)
      real *8, allocatable :: grad(:,:,:,:), gradex(:,:,:,:)
      real *8, allocatable :: hess(:,:,:,:), hessex(:,:,:,:)

      real *8, allocatable :: pote(:,:)
      real *8, allocatable :: grade(:,:,:)
      real *8, allocatable :: hesse(:,:,:)

      real *8, allocatable :: src(:,:)
      character *12 fname1
      character *8 fname2
      character *9 fname3
      
      real*8 errpe
      character*1 type
      external fforce, uexact

      open(21,file='test5_pot.data')
      open(22,file='fvals.data')
      open(26,file='potex.data')

c     set parameters
      call prini(6,13)
c
c      allocate(targs(ndim,ntarg))
      allocate(targ(ndim))
    
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
      eps=1.0d-6
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
c     nboxes, nlevels, ltree, rintl(0:200)
c     ikernel - 0: Yukawa kernel; 1: Laplace kernel; 2: square root Laplace kernel
      ikernel = 1
c     beta - the parameter in the Yukawa kernel or the exponent of the
c     power function kernel
      beta=0.1d-5
  
      ifpgh=1
      ifpghtarg=1
c      allocate(pot(nd,npbox,mxboxes))
      call cpu_time(t1)
      call vol_tree_mem(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,fforce,nd,dpars,zpars,ipars,ifnewtree,
     2    nboxes,nlevels,ltree,rintl)
      write(*,*) 'nboxes =', nboxes
      write(*,*) 'nlevels=', nlevels
      write(*,*) 'ltree=', ltree
c     still not sure what is rintl, but okay
      call cpu_time(t2)
c
c     2. allocate memory and then call vol_tree_build
c     (itree, iptr) defines the logic structure
c     while (centers, boxsize) defines the geometry
c     and fvals defines the function values on the vol grid
      allocate(fvals(nd,npbox,mxboxes))
c
c     Is level restriction enforced? Yes.
      call vol_tree_build(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,fforce,nd,dpars,zpars,ipars,rintl,
     2    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals)


c     3. call print_tree_matlab for visualization
c     make sure you can use it

      ifplot=1
      if (ifplot.eq.1) then
         fname1 = 'trgtree.data'
         fname2 = 'src.data'
         fname3 = 'targ.data'
      
         nt=0
         ns=0
c        given ns and nt
c        generate random sources and targets
c        and sort into the tree?
         call print_tree_matlab(ndim,itree,ltree,nboxes,centers,
     1       boxsize,nlevels,iptr,ns,src,nt,targs,fname1,fname2,fname3)
      endif
      write(*,*) 'tree printed out for visualization'

c
      nhess = ndim*(ndim+1)/2
      allocate(pot(nd,npbox,mxboxes))
      allocate(grad(nd,ndim,npbox,mboxes))
      allocate(hess(nd,nhess,npbox,mxboxes))
      allocate(pote(nd,ntarg))
      allocate(grade(nd,ndim,ntarg),hesse(nd,nhess,ntarg))
c
      allocate(potex(nd,npbox,mxboxes))
      allocate(gradex(nd,ndim,npbox,mxboxes))
      allocate(hessex(nd,nhess,npbox,mxboxes))
c
       dt=10-2
        call boxfgt(nd,ndim,dt,eps,ipoly,iperiod,norder,npbox,
     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals,
     2    ifpgh,pot,grad,hess,ifnewtree,ntarg,targs,
     3    ifpghtarg,pote,grade,hesse,timeinfo)

c
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
               sum_err=sum_err+abs(potex(1,j,ibox)
     1                    -pot(1,j,ibox))**2
               sum_u=sum_u+potex(1,j,ibox)**2
             write(21,*)ibox,j,pot(1,j,ibox),potex(1,j,ibox)
             enddo
          endif
        enddo
      enddo
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'err2_rel=',err2_rel

      end


c     as rhsfun in test_boxfgt*.f
      subroutine fforce(nd,xyz,dpars,zpars,
     1           ipars,u)
      implicit real *8 (a-h,o-z)
      integer nd, ndim, ipars(*)
      real *8 dpars(*), xyz(nd), t, u(nd)
      complex *16 zpars(*)
c
      pi=3.1415926535897932d0
c
      r=(xyz(1))**2+(xyz(2))**2
      r2=(xyz(1)-0.1)**2+(xyz(2)-0.1)**2
      r3=(xyz(1)+0.15)**2+(xyz(2)-0.1)**2
      w=250
      w2=w*w
      u(1)=(-4*w+4*w2*r)*dexp(-r*w)+
     1       (-4*w+4*w2*r2)*dexp(-r2*w)+
     2       (-4*w+4*w2*r3)*dexp(-r3*w)
      u(2)=-8.0d0*pi**2*dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)
      end subroutine

c     as rhsfun in test_boxfgt*.f
      subroutine uexact(nd,xyz,dpars,zpars,
     1           ipars,u)
      implicit real *8 (a-h,o-z)
      integer nd, ndim, ipars(*)
      real *8 dpars(*), xyz(*), t, u(nd)
      complex *16 zpars(*)
c
      pi=3.1415926535897932d0
      x=xyz(1)
      y=xyz(2)
      t=dpars(1)
c
      r=(xyz(1))**2+(xyz(2))**2
      r2=(xyz(1)-0.1)**2+(xyz(2)-0.1)**2
      r3=(xyz(1)+0.15)**2+(xyz(2)-0.1)**2
      w=250
      u(1)=dexp(-r*w)+dexp(-r2*w)+
     1          dexp(-r3*w)
      u(2)=dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)
c
      end subroutine
