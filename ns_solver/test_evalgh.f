c test2
c test the subroutine evalf
c subroutine defines the array that represents the
c     right hand side of the equation.

c
      program test
      implicit real *8 (a-h,o-z)
      integer ipars(100), ltree
      integer iptr(8), icoll,ipoly
c     what is iptr? add an unpack subroutine?
      real*8 done, pi, rintl(0:200)
      real*8 dpars(1000)
      complex*16 zpars(10)
      complex*16 zk, ztmp, ima, zz
      real *8 tt
      parameter(nd=2)
      parameter(norder=8)
      parameter(ndim=2)
c
      parameter(mxboxes=100000)
      parameter(mxltree=200000)
      parameter(mxlevels=30)
      integer, allocatable :: itree(:)
      integer, allocatable :: isrc(:), isrcse(:,:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8, allocatable :: src(:,:)
      real *8, allocatable :: fright(:,:,:)
      real *8, allocatable :: grad(:,:,:,:)
      real *8, allocatable :: hess(:,:,:,:)

      real *8, allocatable :: grad_ex1(:,:,:)
      real *8, allocatable :: grad_ex2(:,:,:)
      real *8, allocatable :: hess_ex1(:,:,:)
      real *8, allocatable :: hess_ex2(:,:,:)

      integer nlevels
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      external fforce
      external fforce_grad1
      external fforce_grad2
      external fforce_hess1
      external fforce_hess2
  
      open(21,file='test4_grad.data')
      open(22,file='test4_hess.data')

c     set parameters
      call prini(6,13)
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

      tt=0.0d0
      dpars(1)=0.0d0

c
c     subroutine vol_tree_mem returns:
c     nboxes, nlevels, ltree, rintl(0:200)
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
      allocate(fvals(nd,npbox,nboxes),centers(ndim,nboxes))
      allocate(boxsize(0:nlevels),itree(ltree))
      allocate(grad(nd,ndim,norder**ndim,nboxes))
      allocate(hess(nd,ndim*(ndim+1)/2,norder**ndim,nboxes))
      allocate(grad_ex1(nd,npbox,nboxes))
      allocate(grad_ex2(nd,npbox,nboxes))
      allocate(hess_ex1(nd,npbox,nboxes))
      allocate(hess_ex2(nd,npbox,nboxes))
c
c     Is level restriction enforced? Yes.
      call vol_tree_build(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,fforce,nd,dpars,zpars,ipars,rintl,
     2    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals)


c
      call evalgh(ndim,nd,ipoly,nlevels,ltree,itree,iptr,
     1       boxsize,norder,nboxes,fvals,grad,hess)

c
      call evalf(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,
     2       dpars,zpars,ipars,fforce_grad1,grad_ex1)
c
      call evalf(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,
     2       dpars,zpars,ipars,fforce_grad2,grad_ex2)
c
      call evalf(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,
     2       dpars,zpars,ipars,fforce_hess1,hess_ex1)
c
      call evalf(ndim,nd,ipoly,norder,nboxes,nlevels,
     1       ltree,itree,iptr,centers,boxsize,
     2       dpars,zpars,ipars,fforce_hess2,hess_ex2)


c
      sum_err=0.0d0
      sum_u=0.0d0

      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               sum_err=sum_err+abs(grad_ex1(1,j,ibox)
     1                    -grad(1,1,j,ibox))**2
               sum_u=sum_u+grad_ex1(1,j,ibox)**2
               write(21,*)grad(1,1,j,ibox),grad(1,2,j,ibox),
     1        grad(2,1,j,ibox),grad(2,2,j,ibox)
             enddo
          endif
        enddo
      enddo
   
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'grad x direction err2_rel=',err2_rel
c
      sum_err=0.0d0
      sum_u=0.0d0

      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               sum_err=sum_err+abs(grad_ex2(1,j,ibox)
     1                    -grad(1,2,j,ibox))**2
               sum_u=sum_u+grad_ex2(1,j,ibox)**2
             enddo
          endif
        enddo
      enddo
   
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'grad y direction err2_rel=',err2_rel

c
      sum_err=0.0d0
      sum_u=0.0d0

      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               sum_err=sum_err+abs(hess_ex1(1,j,ibox)
     1                    -hess(1,1,j,ibox))**2
               sum_u=sum_u+hess_ex1(1,j,ibox)**2
             write(22,*)hess(1,1,j,ibox),hess(1,2,j,ibox),
     1    hess(2,1,j,ibox)
             enddo
          endif
        enddo
      enddo
   
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'hess xx direction err2_rel=',err2_rel

      

c
      sum_err=0.0d0
      sum_u=0.0d0

      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               sum_err=sum_err+abs(hess_ex2(1,j,ibox)
     1                    -hess(1,2,j,ibox))**2
               sum_u=sum_u+hess_ex2(1,j,ibox)**2
             write(22,*)hess(1,1,j,ibox),hess(1,2,j,ibox),
     1    hess(2,1,j,ibox)
             enddo
          endif
        enddo
      enddo
   
c
      err2_rel=dsqrt(sum_err/sum_u)
      write(*,*) 'hess xy direction err2_rel=',err2_rel

      end program



c     as rhsfun in test_boxfgt*.f
      subroutine fforce(nd,xyz,dpars,zpars,
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
      u(1)=dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)
      u(2)=1.0d0
c
      end subroutine

c     as rhsfun in test_boxfgt*.f
      subroutine fforce_grad1(nd,xyz,dpars,zpars,
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
      u(1)=2.0d0*pi*dcos(2.0d0*pi*x)*dcos(2.0d0*pi*y)
      u(2)=0.0d0
c
      end subroutine


c     as rhsfun in test_boxfgt*.f
      subroutine fforce_grad2(nd,xyz,dpars,zpars,
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
      u(1)=-2.0d0*pi*dsin(2.0d0*pi*x)*dsin(2.0d0*pi*y)
      u(2)=0.0d0
c
      end subroutine


c     as rhsfun in test_boxfgt*.f
      subroutine fforce_hess1(nd,xyz,dpars,zpars,
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
      u(1)=-4.0d0*pi**2*dsin(2.0d0*pi*x)*dcos(2.0d0*pi*y)
      u(2)=0.0d0
c
      end subroutine
c     as rhsfun in test_boxfgt*.f
      subroutine fforce_hess2(nd,xyz,dpars,zpars,
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
      u(1)=-4.0d0*pi**2*dcos(2.0d0*pi*x)*dsin(2.0d0*pi*y)
      u(2)=0.0d0
c
      end subroutine

