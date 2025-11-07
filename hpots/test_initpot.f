      program test
      implicit real*8 (a-h,o-z)
      integer iptr(12), ipars(100)
      real*8 epsvals(5), deltas(20)
      real*8 pps(20,5), dpars(1000)
      real *8 rintl(0:200)
      real *8 timeinfo(100)
      complex*16 ima,zz,ztmp,zk
      complex*16 zpars(10)
c     allocatables
c     targets and outputs
      real*8, allocatable:: targs(:,:)
      real*8, allocatable:: src(:,:)
c
c      real*8, allocatable:: pot(:,:,:)
      real*8, allocatable:: pot(:,:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
      real*8, allocatable:: potex(:,:,:)
      real*8, allocatable:: gradex(:,:,:,:)
      real*8, allocatable:: hessex(:,:,:,:)
c
c      real*8, allocatable:: pote(:,:)
      real*8, allocatable:: pote(:)
      real*8, allocatable:: grade(:,:,:)
      real*8, allocatable:: hesse(:,:,:)
      real*8, allocatable:: potexe(:,:)
      real*8, allocatable:: gradexe(:,:,:)
      real*8, allocatable:: hessexe(:,:,:)
c
c     what is whts?
      real *8, allocatable :: xref(:,:),wts(:)
      real *8, allocatable :: targ(:)
c
c     the tree structure
      integer, allocatable:: itree(:)
      real*8, allocatable:: fvals(:,:,:)
      real*8, allocatable:: centers(:,:)
      real*8, allocatable:: boxsize(:)
      character*1 type
      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      external rhsfun, uexact
c
c-------------------------------------
c
      call prini(6,13)
      pi=3.1415926535897932d0
c
      epsvals(1) = 1d-3
      epsvals(2) = 1d-6
      epsvals(3) = 1d-9
      epsvals(4) = 1d-12
c
      do i=1,20
         deltas(i)=10.0d0**(-i+1)
      enddo
c
      ipoly=1
      iperiod=1
      norder=8
c      norder=16
c     nd: number of different densities
      nd=1
c     ndim: dimension of the underlying space
      ndim=2
c
      ifpgh=1
      ifpghtarg=0
c      ntarg=1000000
      ntarg=0
c
c     eps = 1e-12
      eps=epsvals(3)
c     delta=1e-3
c      delta=deltas(4)
      delta=1.0d-5
c     eps and delta defined
c     start the test the boxfgt
      
      ifnewtree=0
c     Lp norm
      iptype=2
c     what is eta?
c
      npbox=norder**ndim
      type='f'
c
      boxlen=1.0d0
c
c     allocate targs, pote, grade, hesse
      nsrc=ntarg
c      allocate(targs(ndim,ntarg),src(ndim,nsrc),pote(nd,ntarg))
      allocate(targs(ndim,ntarg),src(ndim,nsrc),pote(ntarg))
      allocate(grade(nd,ndim,ntarg),hesse(nd,nhess,ntarg))
c
      do i=1,ntarg
         do j=1,ndim
            targs(j,i) = hkrand(0)-0.5d0
         enddo
        write(101,*) targs(1,i), targs(2,i)
      enddo

      call prin2('delta=*',delta,1)

c
c-------------------------------------
c
      epstree=1.0d-10
c      epstree=1.0d-11
c      epstree=1.0d-12
c     sometimes not quite robust, creating
c     a tree that is too small
c     but there is not much to be done about it
c
      eta=1.0d0
      zk=1.0d0
c     start the refining from which level...
c
      call cpu_time(t1)
c     iptype: choose an Lp norm
c     tree mem: (nboxes,nlevels,ltree,rintl)
c     what is rintl?
      call vol_tree_mem(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1     norder,iptype,eta,rhsfun,nd,dpars,zpars,ipars,
     2     ifnewtree,nboxes,nlevels,ltree,rintl)
      call cpu_time(t2)
      write(*,*) 'time of vol_tree_mem=', t2-t1
      write(*,*) 'parameters:'
      write(*,*) 'nboxes =', nboxes
      write(*,*) 'nlevels=', nlevels 
      write(*,*) 'ltree=', ltree 
c
c     mem alloc for the tree structure:
c     logic structure: (itree, iptr)
c     geometry: (centers, boxsize)
c     function values: fvals
c
      allocate(fvals(nd,npbox,nboxes),centers(ndim,nboxes))
      allocate(boxsize(0:nlevels),itree(ltree))
c
c     Is level restriction enforced? Yes.
      call vol_tree_build(ndim,ipoly,iperiod,epstree,zk,boxlen,
     1    norder,iptype,eta,rhsfun,nd,dpars,zpars,ipars,rintl,
     2    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals)
      call prinf('laddr=*',itree,2*(nlevels+1))
c
c      call vol_tree_unpack(iptr, iladdr, ilevel, iparent,     
c     1     nchild, ichild, ncoll, icoll, ltree)
c      write(*,*) 'iladdr =', iladdr
c      write(*,*) 'ilevel =', ilevel
c      write(*,*) 'iparent =', iparent
c      write(*,*) 'nchild =', nchild
c      write(*,*) 'ichild =', ichild
c      write(*,*) 'ncoll =', ncoll
c      write(*,*) 'coll =', icoll
c      write(*,*) 'ltree =', ltree+1
c     why is ltree increased by 1? I'd like to fix this.
c      ltree=ltree-1
c
c     count the number of leaf boxes
c      nlfbox = 0
c      do il=1,nlevels
c        do ibox=itree(iladdr+2*il),itree(iladdr+2*il+1)
c          write(*,*) iladdr+2*il, iladdr+2*il+1
c          write(*,*) ibox, itree(nchild+ibox-1)
c          if(itree(nchild+ibox-1).eq.0) then
c            nlfbox = nlfbox+1
c          endif
c        enddo
c      enddo
c      call prinf('nlfbox=*',nlfbox,1)
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
c-------------------------------------
c
c     alloc mem for output
      nhess = ndim*(ndim+1)/2
c      allocate(pot(nd,npbox,nboxes))
      allocate(pot(npbox,nboxes))
      allocate(grad(nd,ndim,npbox,nboxes))
      allocate(hess(nd,nhess,npbox,nboxes))
c
c     initialize output
      do i=1,nboxes
      do j=1,npbox
      do ind=1,nd
         pot(j,i) = 0
         if (ifpgh.ge.2) then
            do k=1,ndim
               grad(ind,k,j,i) = 0
            enddo
         endif
         if (ifpgh.ge.3) then
            do k=1,nhess
               hess(ind,k,j,i) = 0
            enddo
         endif
      enddo
      enddo
      enddo
c
c     the main task: call the boxfgt
      call cpu_time(t1) 
c      dt=4.0d0*delta
c
      do ll=0,nlevels
        write(*,*) ll, boxsize(ll)
      enddo
c
c------------------------------------------------------------
c
c      call boxfgt(nd,ndim,dt,eps,ipoly,iperiod,norder,npbox,
c     1    nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals,
c     2    ifpgh,pot,grad,hess,ifnewtree,ntarg,targs,
c     3    ifpghtarg,pote,grade,hesse,timeinfo)
c
c     output of initpot: pot, potp
c
      dt=delta
      ifvtarg=1
      ifptarg=0
      write(*,*) 'before calling initpot'
      write(*,*) norder, ltree
      write(*,*) nlevels, nboxes, dt
      write(*,*) iperiod, ntarg, ifvtarg, ifptarg
      write(*,*) eps

      call initpot(norder, ltree, itree, iptr, centers, 
     1     nlevels, boxsize, nboxes, fvals, dt, 
     2     iperiod, ntarg, targs, ifvtarg, ifptarg, 
     3     eps, pot, pote)
c
c------------------------------------------------------------
c
      call cpu_time(t2) 
      call prin2('time taken in fgt=*',t2-t1,1)
      ntot=npbox*nlfbox+ntarg
      pps=(npbox*nlfbox+ntarg+0.0d0)/(t2-t1)
      call prinf('ntotal=*',ntot,1)
      call prin2('speed in pps=*',pps,1)
c
c     allocate mem for exact values
      allocate(potex(nd,npbox,nboxes))
      allocate(gradex(nd,ndim,npbox,nboxes))
      allocate(hessex(nd,nhess,npbox,nboxes))

      allocate(xref(ndim,npbox),wts(npbox))
      itype = 0
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c     itype = 0 -> only xref returned
c     utmp, vtmp, whts: the transf matrices and weights
c
      do i=1, npbox
        write(112,*) xref(1,i), xref(2,i)
      enddo

      allocate(targ(ndim))
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
             do j=1,npbox
               do k=1,ndim
                 targ(k)=centers(k,ibox) + xref(k,j)*bs
               enddo
c
c               pot(j,ibox)=pot(j,ibox)/(4.0d0*pi*delta)
c
c              cent0, xsize0?, d=delta
               dpars(1)=delta
               dpars(2)=-boxlen/2.0d0
               dpars(3)=boxlen/2.0d0
               dpars(4)=-boxlen/2.0d0
               dpars(5)=boxlen/2.0d0
               call uexact(nd,targ,dpars,zpars,ipars,
     1             potex(1,j,ibox),gradex(1,1,j,ibox),
     2             hessex(1,1,j,ibox))
cccc               call rhsfun(nd,targ,dpars,zpars,ipars,fval)

cccc               write(34,*) targ(1), fval, potex(1,j,ibox), pot(1,j,ibox)
               err = abs(potex(1,j,ibox)-pot(j,ibox))
               write(134,*) targ(1), targ(2), 
     1                     potex(1,j,ibox), pot(j,ibox), err
               write(132,*) j,ibox, fvals(1,j,ibox), pot(j,ibox)
             enddo
          endif
        enddo
      enddo



      end program




c-----------------------------------------
      subroutine rhsfun(nd,xyz,dpars,zpars,
     1           ipars,f)
c     the right hand side function
c     defined as the sum of five gaussians
      implicit real *8 (a-h,o-z)
      integer nd, ndim, ipars(*)
      real *8 dpars(*),f(nd),xyz(*)
      real*8 s1, s2, res
      complex *16 zpars(*)
c
      one=1.0d0
      pi=3.1415926535897932d0
c
      s1=2.0d0
      s2=2.0d0
c
      res=dsin(2*s1*pi*xyz(1))*dcos(2*s2*pi*xyz(2))
      f(1)=res

      end subroutine




c
c-------------------------------------
c
      subroutine uexact(nd,targ,dpars,zpars,ipars,
     1    pot,grad,hess)
c     exact solution for the given rhs function
      implicit real*8 (a-h,o-z)
      real*8 targ(*),pot(nd),grad(nd,*),hess(nd,*)
      real*8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)
      real*8 s1, s2
      
      one=1.0d0
      pi=3.1415926535897932d0
c
      s1=2.0d0
      s2=2.0d0
c
      d=dpars(1)
c
c      const1=pi*d
c      const2=1.0d0
      const1=1.0d0
      const2=4.0d0
c
      x1=targ(1)
      x2=targ(2)
c
      pot(1)=const1*dexp(-(s1**2+s2**2)*pi**2*d*const2)*
     1              dsin(2*s1*pi*x1)*dcos(2*s2*pi*x2)
c
c
c
      end subroutine
