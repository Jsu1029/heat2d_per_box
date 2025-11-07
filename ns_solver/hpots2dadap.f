**************************************************************
c                                                
c     Subroutines for the evaluation of heat potentials
c     (in a possibly moving geometry)
c
**************************************************************
c
c
c
c
c
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
      subroutine initpot(norder, ltree, itree, iptr, centers, 
     1           nlevels, boxsize, nboxes, fright, dt,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg, 
     3           eps, pot, potp)
      implicit real*8 (a-h,o-z)
c     global vars
      integer norder, ltree, nboxes, ndim, nlevels
      integer iperiod, ntarg
      integer ifvtarg, ifptarg
      parameter(ndim=2)
      integer itree(ltree), iptr(12)
      real*8 centers(ndim, nboxes), boxsize(0:nlevels)
      real*8 fright(norder**ndim,nboxes), dt
      real*8 targs(ndim,ntarg), eps
      real*8 pot(norder**ndim,nboxes), potp(ntarg)
c     local vars      
      real*8, allocatable:: timeinfo(:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
      real*8, allocatable:: gradp(:,:,:)
      real*8, allocatable:: hessp(:,:,:)
c
      pi=3.1415926535897932d0
      nd=1
      nhess = ndim*(ndim+1)/2
      npbox=norder**ndim
c
      allocate(timeinfo(100))
c
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
          pot(i,j)=0.0d0
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
               pot(j,ibox)=pot(j,ibox)/(4.0d0*pi*dt)
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
          potp(i)=potp(i)/(4.0d0*pi*dt)
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
c     subroutine volpot
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
      subroutine volpot(norder, nordert, ltree, itree, iptr, 
     4           centers, nlevels, boxsize, nboxes, fforce, tpre, 
     1           dt, iperiod, ntarg, targs, ifvtarg, ifptarg, 
     3           eps, vpot, vpotp)
      implicit real*8 (a-h,o-z)
c     global vars
      integer norder, ltree, nboxes, ndim, nlevels
      integer iperiod, ntarg
      integer ifvtarg, ifptarg
      parameter(ndim=2)
      integer itree(ltree), iptr(12)
      real*8 centers(ndim, nboxes), boxsize(0:nlevels)
      real*8 targs(ndim,ntarg), eps, dt
      real*8 vpot(norder**ndim,nboxes), vpotp(ntarg)
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
      real*8, allocatable:: fright(:,:)
      real*8, allocatable:: vpotinc(:,:), vpotpinc(:)
      real*8, allocatable:: xref(:,:), wts(:)
      external fforce
c
      pi=3.1415926535897932d0
      npbox=norder**ndim
      nd=1
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
      allocate(fright(npbox,nboxes))
      allocate(vpotinc(npbox,nboxes))
      allocate(vpotpinc(ntarg))
c
      allocate(xref(ndim,npbox))
      allocate(wts(npbox))
c
      allocate(xtarg(ndim))
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
          vpot(k,i)=0.0d0
        enddo
      enddo
c
c--------------------------------------------------
c     contribution from t(jj): an FGT from last step
      do jj=1,nnodes-1
        tt=tpre+tnodes(jj)
        delta=dt-tnodes(jj)
        scal=dt
c
c       assign fright: fforce(targ,tt,val) val->fright
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
c              a leaf node
c              get the grid and compute the forcing term 
               do j=1,npbox
                 do k=1,ndim
                   xtarg(k)=centers(k,ibox) + xref(k,j)*bs
                 enddo
                 dpars(1)=tt
                 call fforce(nd,xtarg,dpars,zpars,ipars,val)
                 fright(j,ibox)=val
               enddo
            else
c             otherwise assign fright to be zero
              do j=1, npbox
                fright(j,ibox)=0.0d0
              enddo
            endif
          enddo
        enddo
c       fright assigned
c
        do i=1,nboxes
          do j=1, npbox
            vpotinc(j,i)=0.0d0
          enddo
        enddo
c
c       call boxfgt: fright -> vpotinc
c       be careful about the scaling factor
        delta=4.0d0*(dt-tnodes(jj))
        ifnewtree = 0
ccc -------------------------------------------------
ccc     I believe it's the problem of the boxfgt?
ccc     dividing by a small number? unlikely.
ccc     the error estimate of the boxfgt not enough?
ccc     number of terms on the edge?
ccc -------------------------------------------------
c        eps0=eps*1e3
c        eps0=eps*1e2
c        eps0=eps*1e1
        eps0=eps
c        eps0=eps*1e-2
c        eps0=eps*1e-3
        call boxfgt(nd,ndim,delta,eps0,ipoly,iperiod,norder,npbox,
     1       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2       fright,ifpgh,vpotinc,grad,hess,ifnewtree,ntarg,targs,
     3       ifpghtarg,vpotpinc,gradp,hessp,timeinfo)
c
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do j=1, npbox
               vpot(j,i)=vpot(j,i)+vpotinc(j,i)/(4.0d0*pi)
     1                  /(dt-tnodes(jj))*dt*twhts(jj)
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
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
c            a leaf node
c            get the grid and compute the forcing term 
             do j=1,npbox
               do k=1,ndim
                 xtarg(k)=centers(k,ibox) + xref(k,j)*bs
               enddo
               dpars(1)=tt
               call fforce(nd,xtarg,dpars,zpars,ipars,val)
               vpot(j,ibox)=vpot(j,ibox)+val*dt*twhts(jj)
             enddo
          endif
        enddo
      enddo
c



      end subroutine




c**************************************************************
c     subroutine volpotam
c**************************************************************
c     This subroutine computes the volume potential
c     by 2nd and 4th order Adams-Moulton methods
c
c     where the forcing term = f(t,x,u)
c
c  INPUT:
c  norder: degree of poly approx on each leaf node
c  ltree - nboxes: the tree structure
c  fevaln: the handle for the forcing term
c          (can be replaced later)
c  dt: time step
c  iperiod: periodic BC on the unit box (or not)
c           iperiod = 0 => free space BC on the unit box
c           iperiod = 1 => periodic BC on the unit box 
c  ntarg: number of extra targets
c  targs: coords of extra targets
c  ifvtarg: consider volume target or not
c  ifptarg: consider point target or not
c  upre: the solution of the heat eqn at the previous step
c  gpot: the initial potential of the current step
c
c  OUTPUT:
c  usol: solution of the current step on the volume grid
c  usolp: solution of the current step at targ
c         (not implemented yet)
c
c  attention: vpotp not assigned yet in this version
c  fforce has the form f(t,x,u) and the calling sequence
c
c  this version implemented 2nd order Adams-Moulton method only
c  and will include 4th AM method soon
c  the resulting nonlinear equation is solved by secant method
c
c
c**************************************************************
      subroutine volpotam(norder, nordert, ltree, itree, 
     1           iptr, centers, nlevels, boxsize, nboxes, 
     2           fevaln, tpre, dt, iperiod, ntarg, targs, 
     3           ifvtarg, ifptarg, upre, gpot, eps, usol)
      implicit real*8 (a-h,o-z)
c     global vars
      integer norder, ltree, nboxes, ndim, nlevels
      integer iperiod, ntarg
      integer ifvtarg, ifptarg
      parameter(ndim=2)
      integer itree(ltree), iptr(12)
      real*8 centers(ndim, nboxes), boxsize(0:nlevels)
      real*8 targs(ndim,ntarg), eps, dt
      real*8 gpot(norder**ndim,nboxes)
      real*8 upre(norder**ndim,nboxes)
      real*8 usol(norder**ndim,nboxes)
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
      real*8, allocatable:: fright(:,:)
      real*8, allocatable:: fpre(:,:), vpot(:,:)
      real*8, allocatable:: vpotinc(:,:), vpotpinc(:)
      real*8, allocatable:: xref(:,:), wts(:)
      external fevaln
c
      pi=3.1415926535897932d0
      npbox=norder**ndim
      nd=1
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
      allocate(fright(npbox,nboxes))
      allocate(fpre(npbox,nboxes))
      allocate(vpot(npbox,nboxes))
      allocate(vpotinc(npbox,nboxes))
      allocate(vpotpinc(ntarg))
c
      allocate(xref(ndim,npbox))
      allocate(wts(npbox))
c
      allocate(xtarg(ndim))
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
      itype = 0
      ipoly=1
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c
c     only nordert .eq. 2 implemented
c     for now
c     AM2, a.k.a. Crank-Nicolson, 
c     a.k.a. trapezoidal rule
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
        nordert=2
        nnodes=2
        tnodes(1)=0.0d0
        tnodes(2)=dt
        twhts(1)=1.0d0/2.0d0
        twhts(2)=1.0d0/2.0d0
        write(*,*) 'volpotam: only norder=2,4 allowed!'
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
          vpot(k,i)=0.0d0
          fpre(k,i)=0.0d0
        enddo
      enddo
c
c--------------------------------------------------
c     contribution from t(jj): an FGT from last step
      do jj=1,nnodes-1
        tt=tpre+tnodes(jj)
        delta=dt-tnodes(jj)
        scal=dt
        write(*,*) 'jj, nnodes', jj, nnodes
c
c       assign fright: fevaln(targ,tt,val) val->fright
        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
c              a leaf node
c              get the grid and compute the forcing term 
               do j=1,npbox
                 do k=1,ndim
                   xtarg(k)=centers(k,ibox) + xref(k,j)*bs
                 enddo
                 dpars(1)=tt
ccc              modify the calling sequence here
c                 call fforce(nd,xtarg,upre(j,ibox),dpars,zpars,
c     1                ipars,fpre(j,ibox))
c                fevaln(t,x,y,u,f)
                 call fevaln(tt,xtarg(1),xtarg(2),
     1                upre(j,ibox),fpre(j,ibox))
                 fright(j,ibox)=fpre(j,ibox)
               enddo
            else
c             otherwise assign fright to be zero
              do j=1, npbox
                fright(j,ibox)=0.0d0
              enddo
            endif
          enddo
        enddo
c       fright assigned
c
        do i=1,nboxes
          do j=1, npbox
            vpotinc(j,i)=0.0d0
          enddo
        enddo
c
c       call boxfgt: fright -> vpotinc
c       be careful about the scaling factor
        delta=4.0d0*(dt-tnodes(jj))
        ifnewtree = 0
c
        write(*,*) 'calling boxfgt'
        call boxfgt(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2       fright,ifpgh,vpotinc,grad,hess,ifnewtree,ntarg,targs,
     3       ifpghtarg,vpotpinc,gradp,hessp,timeinfo)
c
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do j=1, npbox
c               vpot(j,i)=vpot(j,i)+vpotinc(j,i)/(4.0d0*pi)
c     1                  /(dt-tnodes(jj))*dt*twhts(jj)
               vpotinc(j,i)=vpotinc(j,i)/(4.0d0*pi)
     1                  /(dt-tnodes(jj))*dt*twhts(jj)
               vpot(j,i)=vpot(j,i)+vpotinc(j,i)
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
c
c     however the function value is unknown
c     a nonlinear solver is needed
c     (initial guess: usol of the last step
c     and the initial potential of the current step)
c
c     the current version is nonadaptive
c     in the adaptive version, an adaptive stage
c     is required after the nonlinear solve, 
c     a.k.a. the update of the volpot
c
      jj=nnodes
      tt=tpre+tnodes(jj)
      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
c            a leaf node
c            get the grid and compute the forcing term 
             do j=1,npbox
               do k=1,ndim
                 xtarg(k)=centers(k,ibox) + xref(k,j)*bs
               enddo
c               dpars(1)=tt
c               call fforce(nd,xtarg,dpars,zpars,ipars,val)
c               vpot(j,ibox)=vpot(j,ibox)+val*dt*twhts(jj)
c--------------------------
c              now call secant to solve the scalar nonlinear eqn
               msi = 10
               ndamp = 1
               alpha= dt*twhts(jj)
               g = vpot(j,ibox)+gpot(j,ibox)
               u0 = upre(j,ibox)
               u1 = g+alpha*fpre(j,ibox)
               tol = 1.0d-12
c              be careful about the 
c              calling seq of fevaln for secant
c              which solves
c              u = g +alpha * f(u,x,t)
c              where r is given by a subroutine rfun
c              with calling sequence
c              rfun(t,x,u,f)
               call secant(fevaln,tt,xtarg(1),xtarg(2),
     1              u0,u1,g,alpha,u,fu,tol,msi,isi,
     2              ndamp,ier)
               if(ier .ne. 0) then
                 write(*,*) 'ier=',ier
               endif
c--------------------------
c              assign u to usol: solution of the 
c              current time step
               usol(j,ibox)=u
             enddo
          endif
        enddo
      enddo
c
c
c
      end subroutine





c**************************************************************
c     subroutine vpot_semilin_am2
c**************************************************************
c
c     This subroutine computes the volume potential
c     by 2nd and 4th order Adams-Moulton methods
c
c     where the forcing term = f(t,x,u)
c
c  INPUT:
c  norder: degree of poly approx on each leaf node
c  ltree - nboxes: the tree structure
c  fevaln: the handle for the forcing term
c          (can be replaced later)
c  dt: time step
c  iperiod: periodic BC on the unit box (or not)
c           iperiod = 0 => free space BC on the unit box
c           iperiod = 1 => periodic BC on the unit box 
c  ntarg: number of extra targets
c  targs: coords of extra targets
c  ifvtarg: consider volume target or not
c  ifptarg: consider point target or not
c  un: the solution of the heat eqn at the previous step
c
c  OUTPUT:
c  un: solution of the current step on the volume grid
c  unp: solution of the current step at targ
c         (not implemented yet)
c
c  attention: vpotp not assigned yet in this version
c  fforce has the form f(t,x,u) and the calling sequence
c
c  this version implemented 2nd order Adams-Moulton method only
c  and will include 4th AM method soon
c  the resulting nonlinear equation is solved by secant method
c
c     rintl and iadap added
c     mboxes and mlevels added
c
c**************************************************************
      subroutine vpot_semilin_am2(norder, nordert, ltree, itree, 
     1           iptr, centers, nlevels, mlevels, boxsize, nboxes, 
     2           mboxes, rintl, fevaln, tpre, dt, iperiod, ntarg, 
     3           targs, ifvtarg, ifptarg, iadap, eps, un)
      implicit real*8 (a-h,o-z)
c     global vars
      integer norder, ltree, nboxes, ndim, nlevels
      integer iperiod, ntarg, mboxes, mlevels
      integer ifvtarg, ifptarg
      parameter(ndim=2)
      integer itree(ltree), iptr(12)
      real*8 centers(ndim, mboxes), boxsize(0:mlevels)
      real*8 targs(ndim,ntarg), eps, dt
      real*8 un(norder**ndim,mboxes)
      real *8 rintl(0:mlevels)
c     local vars      
      integer, allocatable:: ipars(:)
      real*8, allocatable:: timeinfo(:)
      real*8, allocatable:: unm1(:,:)
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
      real*8, allocatable:: fright(:,:)
      real*8, allocatable:: fn(:,:), vpot(:,:)
      real*8, allocatable:: vpotinc(:,:), vpotpinc(:)
      real*8, allocatable:: xref(:,:), wts(:)
      external fevaln
c      open(333,file='vpot_am2.txt')
c      open(334,file='vpot_am22.txt')
c      open(335,file='vpot_am23.txt')
c
      pi=3.1415926535897932d0
      npbox=norder**ndim
      nd=1
      ierrr=0
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
      allocate(fright(npbox,mboxes))
      allocate(unm1(npbox,mboxes))
      allocate(fn(npbox,mboxes))
      allocate(vpot(npbox,mboxes))
      allocate(vpotinc(npbox,mboxes))
      allocate(vpotpinc(ntarg))
c
      allocate(xref(ndim,npbox))
      allocate(wts(npbox))
c
      allocate(xtarg(ndim))
      allocate(timeinfo(100))
c
      allocate(grad(nd,ndim,npbox,mboxes))
      allocate(hess(nd,nhess,npbox,mboxes))
      allocate(gradp(nd,ndim,ntarg))
      allocate(hessp(nd,nhess,ntarg))
c
      allocate(ipars(100))
      allocate(dpars(1000))
      allocate(zpars(10))



c ---------------debugging------------------------------
c      do ilevel=0, nlevels
c        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
c          if(itree(iptr(4)+ibox-1).eq.0) then
c             a leaf node
c            do i=1, npbox
c             write(333,*)ibox,i,un(i,ibox)
c            enddo
c          endif
c        enddo
c      enddo
c-------------------------------------------------------
c
      itype = 0
      ipoly=1
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c
c     only nordert .eq. 2 implemented
c     for now
c     AM2, a.k.a. Crank-Nicolson, 
c     a.k.a. trapezoidal rule
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
        write(*,*) 'volpot: only norder=2,4 allowed!'
        write(*,*) 'norder reset to 2!'
      endif
c
c
c--------------------------------------------------
c     debugging
c      print out usol for debugging
c      do ilevel=0,nlevels
c        bs = boxsize(ilevel)/2.0d0
c        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
c          if(itree(iptr(4)+ibox-1).eq.0) then
c            do j=1, npbox
c              write(205,*) un(j,ibox)
c            enddo
c          endif
c        enddo
c      enddo
c
c--------------------------------------------------
c     initialize vpot, fn to zero
c     save un -> unm1
c
      do i=1,nboxes
        do k=1,npbox
          vpot(k,i)=0.0d0
        enddo
      enddo
      do i=1,nboxes
         do k=1,npbox
            fn(k,i)=0.0d0
            unm1(k,i)=un(k,i)
         enddo
      enddo
c
      write(*,*) 'fine here,mboxes=', 1,mboxes
c
c--------------------------------------------------
c     contribution from t(jj): an FGT from last step
      scal=dt
      do jj=1,nnodes-1
c       loop over nnodes-1 nodes (in time)
        tt=tpre+tnodes(jj)
c        write(*,*) tpre, tnodes(jj), tt
        delta=(norder-1)*dt-tnodes(jj)
c        write(*,*) 'jj, nnodes', jj, nnodes
c
        if(jj .eq. nnodes-3) then
        elseif(jj .eq. nnodes-2) then
        elseif(jj .eq. nnodes-1) then
c         evaluate fevaln at the grid pts
c         using un, solution of the last step
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
                  call fevaln(tt,xtarg(1),xtarg(2),
     1                 un(j,ibox),fn(j,ibox))
                  if(abs(tt-0.0d-3).lt.1.0d-10) then
                    write(206,*) xtarg(1), xtarg(2),tt, 
     1                           un(j,ibox),fn(j,ibox), j, ibox
                  else
                    write(207,*) xtarg(1), xtarg(2),tt, 
     1                           un(j,ibox),fn(j,ibox), j, ibox
                  endif
cccccc            problem found: un and fn are different each time
cccccc            the program is run
                enddo
c
              endif
            enddo
          enddo

c
c         adding the contribution from un
c         and fn (initpot+volpot)
          do i=1, nboxes
            do k=1, npbox
              fright(k,i)=fn(k,i)*scal*twhts(jj)+un(k,i)
c              if(itree(iptr(4)+i-1).eq.0) then
c                write(335,*)jj,i,k, fn(k,i), un(k,i), scal, twhts(jj)
c              endif
            enddo
          enddo
c
        endif
c
        do i=1, nboxes
          do k=1, npbox
            vpotinc(k,i)=0.0d0
          enddo
        enddo
c
c       call boxfgt
c       fright -> vpotinc
c       remember to scale
c
        delta=4.0d0*(dt-tnodes(jj))
        ifnewtree = 0
        write(*,*) 'calling boxfgt'
c       be careful: dt*twhts multiplied
c       above (in fright)
c
        write(*,*) nd, ndim, delta, eps
        write(*,*) ipoly, iperiod, norder, npbox
        write(*,*) nboxes, nlevels, ifpgh, ifnewtree
        write(*,*) ntarg,  ifpghtarg
c
c ---------------debugging------------------------------
c      do ilevel=0, nlevels
c        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
c          if(itree(iptr(4)+ibox-1).eq.0) then
c             a leaf node
c            do i=1, npbox
c             write(334,*)ibox,i,un(i,ibox)
c            enddo
c          endif
c        enddo
c      enddo
c-------------------------------------------------------

        call boxfgt(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2       fright,ifpgh,vpotinc,grad,hess,ifnewtree,ntarg,targs,
     3       ifpghtarg,vpotpinc,gradp,hessp,timeinfo)
c
      write(*,*) 'fine here', 2
c
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do j=1, npbox
               vpotinc(j,i)=vpotinc(j,i)/(pi*delta)
               vpot(j,i)=vpot(j,i)+vpotinc(j,i)
c               write(334,*)jj, i,j, vpot(j,i), fright(j,i)
            enddo
          endif
        enddo
c
c     end of the jj loop 
      enddo
c
      write(*,*) 'fine here', 3
c
cccccc TBD:
c     update fnm2, fnm1 for higher order
c      if(nordert .eq. 4) then
c      endif
c
c--------------------------------------------------
c
c     contribution from t(nnodes): function value at 
c     the current time step, no transform but a nonlinear solve (scalar)
c
      write(*,*) 'doing the nonlinear solve'
      jj=nnodes
      tt=tpre+tnodes(jj)
      do ilevel=0, nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
c         pt by pt, get the coords
c         and do a pointwise nonlinear solve
            do j=1, npbox
              do k=1, ndim
                xtarg(k)=centers(k,ibox) + xref(k,j)*bs
              enddo
c
              msi=20
              ndamp=1
              g=vpot(j,ibox)
              alpha=scal*twhts(jj)
c
c             initial values
              u0=un(j,ibox)
              u1=g+alpha*fn(j,ibox)
              tol=1.0d-12
c
              call secant(fevaln,tt,xtarg(1),xtarg(2),u0,u1,g,alpha,u,
     1             fu,tol,msi,isi,ndamp,ier)
              if(ier .gt. 0) then
                write(111,*) 'secant, ier=1', j, ibox
                ierrr = 1
              endif
c
              un(j,ibox)=u
c             solution updated to that of the new time step
cccccc        print out for debugging
c              write(202,*) j, ibox, u0, u1, g, u 
c           end of the pointwise loop
            enddo
c
          else
            do j=1, npbox
              un(j,ibox)=0.0d0
            enddo
          endif
c       end of the loop over boxes on a level
        enddo
c     end of the loop over levels 
      enddo


c
      if(ierrr .gt. 0) then
        write(*,*) 'ierrr =1, failure of secant'
        write(*,*) 'see fort.111 for details'
      endif
c
      if(iadap .gt. 0 ) then
        write(*,*) ' '
        write(*,*) 'calling am2_adap_tree'
c
        write(*,*) nd, ndim, delta, eps
        write(*,*) ipoly, iperiod, norder, npbox
        write(*,*) nboxes, nlevels, ifpgh, ifnewtree
        write(*,*) ntarg,  ifpghtarg
c ---------------debugging------------------------------
      do ilevel=0, nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
c             a leaf node
            do i=1, npbox
             write(333,*)ibox,i,un(i,ibox)
            enddo
          endif
        enddo
      enddo
c-------------------------------------------------------
        call am2_adap_tree(norder, ltree, itree,
     1       iptr, rintl, centers, nlevels, mlevels, boxsize, 
     2       nboxes, mboxes, fevaln, tt, iperiod, eps, alpha,
     3       vpot, unm1, un, iadap)
c ---------------debugging------------------------------
      do ilevel=0, nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
c             a leaf node
            do i=1, npbox
             write(334,*)ibox,i,un(i,ibox)
            enddo
          endif
        enddo
      enddo
c-------------------------------------------------------

      endif
c
      nleaf=0
      do ilevel=0, nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
            nleaf=nleaf+1
          else
            do i=1, npbox
               un(i,ibox) =0.0d0
            enddo
          endif
        enddo
      enddo
c      write(*,*)'nleaf=',nleaf
c

      end subroutine





c**************************************************************
c     subroutine vpot_semilin_am24
c**************************************************************
c
c     This subroutine computes the volume potential
c     by 2nd to 4th order Adams-Moulton methods
c
c     where the forcing term = f(t,x,u)
c
c  INPUT:
c  norder: degree of poly approx on each leaf node
c  ltree - nboxes: the tree structure
c  fevaln: the handle for the forcing term
c          (can be replaced later)
c  dt: time step
c  iperiod: periodic BC on the unit box (or not)
c           iperiod = 0 => free space BC on the unit box
c           iperiod = 1 => periodic BC on the unit box 
c  ntarg: number of extra targets
c  targs: coords of extra targets
c  ifvtarg: consider volume target or not
c  ifptarg: consider point target or not
c     fnm2: the inhomogeneous term f(u,x,t) at t_{n-2}, 
c           only needed if norder == 4
c     fnm1: the inhomogeneous term f(u,x,t) at t_{n-1}, 
c           only needed if norder >=3
c  un: the solution of the heat eqn at the previous step
c
c  OUTPUT:
c  un: solution of the current step on the volume grid
c  unp: solution of the current step at targ
c         (not implemented yet)
c     fnm2: the inhomogeneous term f(u,x,t) at t_{n-1}, 
c           only needed if norder == 4
c     fnm1: the inhomogeneous term f(u,x,t) at t_{n}, 
c           only needed if norder >=3
c
c  attention: vpotp not assigned yet in this version
c  fforce has the form f(t,x,u) and the calling sequence
c
c  rmk: this version has the input of ifnewtree, which 
c       controls the adaptivity
c       it has adaptivity implemented inside
c
c**************************************************************
      subroutine vpot_semilin_am24(norder, nordert, ltree, itree, 
     1           iptr, centers, nlevels, boxsize, nboxes, 
     2           fevaln, tpre, dt, iperiod, ntarg, targs, 
     3           ifnewtree, ifvtarg, ifptarg, eps, fnm2, fnm1,un)
      implicit real*8 (a-h,o-z)
c     global vars
      integer norder, ltree, nboxes, ndim, nlevels
      integer iperiod, ntarg
      integer ifvtarg, ifptarg
      parameter(ndim=2)
      integer itree(ltree), iptr(12)
      real*8 centers(ndim, nboxes), boxsize(0:nlevels)
      real*8 targs(ndim,ntarg), eps, dt
      real*8 fnm2(norder**ndim,nboxes)
      real*8 fnm1(norder**ndim,nboxes)
      real*8 un(norder**ndim,nboxes)
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
      real*8, allocatable:: fright(:,:)
      real*8, allocatable:: fn(:,:), vpot(:,:)
      real*8, allocatable:: vpotinc(:,:), vpotpinc(:)
      real*8, allocatable:: xref(:,:), wts(:)
      external fevaln
c
      pi=3.1415926535897932d0
      npbox=norder**ndim
      nd=1
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
      allocate(fright(npbox,nboxes))
      allocate(fn(npbox,nboxes))
      allocate(vpot(npbox,nboxes))
      allocate(vpotinc(npbox,nboxes))
      allocate(vpotpinc(ntarg))
c
      allocate(xref(ndim,npbox))
      allocate(wts(npbox))
c
      allocate(xtarg(ndim))
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
      elseif(nordert .eq. 3) then
        nnodes=3
        tnodes(1)=0.0d0
        tnodes(2)=dt
        tnodes(3)=2*dt
        twhts(1)=-1.0d0/12.0d0
        twhts(2)= 8.0d0/12.0d0
        twhts(3)= 5.0d0/12.0d0
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
c
c--------------------------------------------------
c     debugging
c      print out usol for debugging
c      do ilevel=0,nlevels
c        bs = boxsize(ilevel)/2.0d0
c        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
c          if(itree(iptr(4)+ibox-1).eq.0) then
c            do j=1, npbox
c              write(205,*) un(j,ibox)
c            enddo
c          endif
c        enddo
c      enddo
c
c--------------------------------------------------
c     initialize vpot, fn to zero
c
      do i=1,nboxes
        do k=1,npbox
          vpot(k,i)=0.0d0
        enddo
      enddo
      do i=1,nboxes
         do k=1,npbox
            fn(k,i)=0.0d0
         enddo
      enddo
c
c--------------------------------------------------
c     contribution from t(jj): an FGT from last step
      scal=dt
      do jj=1,nnodes-1
c       loop over nnodes-1 nodes (in time)
        tt=tpre+tnodes(jj)
c        write(*,*) tpre, tnodes(jj), tt
        delta=(nordert-1)*dt-tnodes(jj)
c        write(*,*) 'jj, nnodes', jj, nnodes
c
        if(jj .eq. nnodes-3) then
          do i=1,nboxes
            do k=1,npbox
              fright(k,i)=fnm2(k,i)*scal*twhts(jj)
            enddo
          enddo
        elseif(jj .eq. nnodes-2) then
          do i=1,nboxes
            do k=1,npbox
              fright(k,i)=fnm1(k,i)*scal*twhts(jj)
            enddo
          enddo
        elseif(jj .eq. nnodes-1) then
c         evaluate fevaln at the grid pts
c         using un, solution of the last step
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
                  call fevaln(tt,xtarg(1),xtarg(2),
     1                 un(j,ibox),fn(j,ibox))
c                  write(206,*) xtarg(1), xtarg(2),tt, 
c     1            un(j,ibox),fn(j,ibox)
                enddo
c
              endif
            enddo
          enddo
c
c         adding the contribution from un
c         and fn (initpot+volpot)
          do i=1, nboxes
            do k=1, npbox
              fright(k,i)=fn(k,i)*scal*twhts(jj)+un(k,i)
c              if(itree(iptr(4)+i-1).eq.0) then
c                write(204,*) fn(k,i), un(k,i), scal, twhts(jj)
c              endif
            enddo
          enddo
c       end of the jj if
        endif
c
        do i=1, nboxes
          do k=1, npbox
            vpotinc(k,i)=0.0d0
          enddo
        enddo
c
c       call boxfgt
c       fright -> vpotinc
c       remember to scale
c
cccccc        delta=4.0d0*(dt-tnodes(jj))
        delta = 4.0d0*delta
        ifnewtree = 0
        write(*,*) 'calling boxfgt'
c       be careful: dt*twhts multiplied
c       above (in fright)
c
c        write(*,*) nd, ndim, delta, eps
c        write(*,*) ipoly, iperiod, norder, npbox
c        write(*,*) nboxes, nlevels, ifpgh, ifnewtree
c        write(*,*) ntarg,  ifpghtarg
c
        call boxfgt(nd,ndim,delta,eps,ipoly,iperiod,norder,npbox,
     1       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2       fright,ifpgh,vpotinc,grad,hess,ifnewtree,ntarg,targs,
     3       ifpghtarg,vpotpinc,gradp,hessp,timeinfo)
c
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do j=1, npbox
               vpotinc(j,i)=vpotinc(j,i)/(pi*delta)
               vpot(j,i)=vpot(j,i)+vpotinc(j,i)
c               write(203,*) j, i, vpot(j,i), fright(j,i)
            enddo
          endif
        enddo
c
c     end of the jj loop 
      enddo
c
cccccc
c     done with the boxfgt's, contributions saved in vpot
cccccc
c
c     update fnm2, fnm1 for higher order
      if(nordert .eq. 4) then
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do k=1, npbox
              fnm2(k,i)=fnm1(k,i)
            enddo
          endif
        enddo
      endif
c
      if(nordert .ge. 3) then
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do k=1, npbox
              fnm1(k,i)=fn(k,i)
            enddo
          endif
        enddo
      endif
c
c--------------------------------------------------
c
c     contribution from t(nnodes): function value at 
c     the current time step, no transform but a nonlinear solve (scalar)
c
      jj=nnodes
      tt=tpre+tnodes(jj)
      do ilevel=0, nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
c         pt by pt, get the coords
c         and do a pointwise nonlinear solve
            do j=1, npbox
              do k=1, ndim
                xtarg(k)=centers(k,ibox) + xref(k,j)*bs
              enddo
c
              msi=20
              ndamp=1
              g=vpot(j,ibox)
              alpha=scal*twhts(jj)
c
c             initial values
              u0=un(j,ibox)
              u1=g+alpha*fn(j,ibox)
              tol=1.0d-12
c
              call secant(fevaln,tt,xtarg(1),xtarg(2),u0,u1,g,alpha,u,
     1             fu,tol,msi,isi,ndamp,ier)
              if(ier .gt. 0) then
                write(*,*) 'secant, ier=1', j, ibox
              endif
c
              un(j,ibox)=u
c             solution updated to that of the new time step
cccccc        print out for debugging
c              write(202,*) j, ibox, u0, u1, g, u 
c           end of the pointwise loop
            enddo
c
          else
            do j=1, npbox
              un(j,ibox)=0.0d0
            enddo
          endif
c       end of the loop over boxes on a level
        enddo
c     end of the loop over levels 
      enddo
c
c----------------------------------------------
c     lastly, adjust the tree by calling
c     subroutine am24_adap_tree
c     which adjusts the tree so that u and f(t,x,u)
c     are both resolved at the current step
c     and interpolates the history data
c     which is fnm1 and fnm2
c
c      call am24_adap_tree(norder, nordert, ltree, itree, 
c     1     iptr, centers, nlevels, boxsize, nboxes, 
c     2     fevaln, tt, iperiod, eps, fnm2, fnm1, un)
c
c     checks the resolution of un and fevaln(tt,un)
c     refine or coarsen if needed
c     interpolate fnm2, fnm1, solve for un again,
c     since it's a pointwise nonlinear solve
c
c
      end subroutine




c----------------------------------------------
c     subroutine am24_adap_tree
c----------------------------------------------
c
c     adapt the tree so that un and f(tt,un)
c     are both resolved but not over-resolved
c
c     maintain the functions fnm2 and fnm1
c     for the solution of the next step
c
c     needs the input of alpha, the coefficient
c     in front of the nonlinear term
c     for the re-solve step
c
c     eventually it's the un, fn, fnm1, fnm2
c     arrays that need to be resolved
c
c     while vpot and unm1 are needed for the re-solve
c
c     output:
c     fnm2, fnm1: interpolated to the new tree
c     un: solved again on the new tree
c     the tree: updated so that the above arrays
c               are resolved but not over-resolved
c
c----------------------------------------------
c     presumably the vector version would be 
c     am24_adap_tree_nd
c
      subroutine am24_adap_tree(norder, nordert, ltree, 
     1     itree, iptr, rintl, centers, nlevels, mlevels, boxsize, 
     2     nboxes, mboxes, fevaln, tt, iperiod, eps, alpha,
     3     vpot, unm1, fnm2, fnm1, un, rhsfun1, rhsfun2, uexact)
      implicit real*8 (a-h,o-z)
c     global vars
      integer norder, ltree, nboxes, ndim, nlevels
      integer mboxes, mlevels, nd
      integer iperiod, ntarg
      parameter(ndim=2)
      parameter(nd=1)
      integer itree(ltree), iptr(12)
      real*8 centers(ndim, mboxes), boxsize(0:mlevels)
      real*8 tt, eps, alpha, f0(1)
      real*8 fnm2(nd,norder**ndim,mboxes)
      real*8 fnm1(nd,norder**ndim,mboxes)
      real*8 vpot(nd,norder**ndim,mboxes)
      real*8 unm1(nd,norder**ndim,mboxes)
      real*8 un(nd,norder**ndim,mboxes)
      complex *16 zk
c     local vars      
      integer, allocatable:: ipars(:)
c     flags for refinement and coarsening
      integer isgn(ndim,2**ndim)
      integer, allocatable:: nblock(:)
      integer, allocatable:: ilevstart(:)
      integer, allocatable:: nboxid(:), ifdelete(:)
      integer, allocatable:: ifrefine(:)
      integer, allocatable:: irefinelev(:), imaxrefinelev(:)
c     real allocatables
      real*8, allocatable:: dpars(:)
c     need more of these fvals thing?
      real*8, allocatable :: fn(:,:,:)
      real*8, allocatable :: uallvals(:,:,:), uallcoefs(:,:,:)
      real*8, allocatable :: fcoefs(:,:,:), ucoefs(:,:,:)
      real*8, allocatable :: fcoefs1(:,:,:), fcoefs2(:,:,:)
      real*8, allocatable :: ucoefs1(:,:,:), vcoefs(:,:,:)
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8, allocatable :: rmask(:), umat_nd(:,:,:)
      real *8 rintl(0:mlevels)
      real *8 xy(ndim)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
c     change the following names
      real *8, allocatable:: polyvc(:,:,:,:)
c     complex allocatables      
      complex*16, allocatable:: zpars(:)
      external fevaln, rhsfun1, rhsfun2, uexact
c
c     do down the tree
c     check resolution of un and f(tt,un)
c     on each leaf node
c     subdivide if needed
c     interpolate fnm2, fnm1, vpot, unm1
c     solve for un again
c
c-------------------------------------------------------
c
c      write(*,*) norder, nordert, ltree
c      write(*,*) nlevels, mlevels, nboxes, mboxes
c      write(*,*) tt, iperiod, eps, alpha
c
      mc=2**ndim
      mnbors=3**ndim
      npols = norder**ndim
      npbox = norder**ndim
      ipoly=1
      iptype=2
c
c      write(*,*) npols, npbox, ndim, nd
c
      allocate(dpars(1000))
      allocate(zpars(1000))
      allocate(ipars(1000))
      allocate(fn(nd,npbox,mboxes))
c     allocate coeffs for them all:
c     un, fn, unm1, fnm1, fnm2, vpot
      allocate(ucoefs(nd,npbox,mboxes))
      allocate(fcoefs(nd,npbox,mboxes))
      allocate(fcoefs1(nd,npbox,mboxes))
      allocate(fcoefs2(nd,npbox,mboxes))
      allocate(ucoefs1(nd,npbox,mboxes))
      allocate(vcoefs(nd,npbox,mboxes))
c
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
      allocate(umat_nd(norder,norder,ndim))
      allocate(polyvc(norder,norder,ndim,mc))
c
      call get_child_box_sign(ndim,isgn)
      call get_p2c_interp_matrices(ndim,ipoly,norder,isgn,polyvc)
c     I don't think grad/hess are needed tho
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
      do ilev=0, nlevels
        do ibox = itree(2*ilev+1), itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild .eq. 0) then
            do i=1,npbox
              do j=1,ndim
                xy(j) = grid(j,i)*boxsize(ilev)+centers(j,ibox)
              enddo
c             attention: check the calling sequence
              call fevaln(tt,xy(1),xy(2),un(1,i,ibox),
     1             fn(1,i,ibox))
            enddo
          endif
        enddo
      enddo
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
     1    iptr,boxsize,norder,un,ucoefs,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,unm1,ucoefs1,umat_nd)
c
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
     1    iptr,boxsize,norder,vpot,vcoefs,umat_nd)
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
             do i=1, npbox
               write(301,*) i, ibox, fcoefs(1,i,ibox),ucoefs(1,i,ibox),
     1                      fn(1,i,ibox),un(1,i,ibox)
             enddo
             erra = erra/rsum
             errb = errb/rsum
c            refine if one of them is greater than the tol
             erri = max(erra, errb)
             if(erri.gt.reps*rintl(ilev)) then
                ifrefine(ibox)=1
             endif
c
c             write(223,*) ibox, erra, errb, erri, reps*rintl(ilev), 
c     1                    ifrefine(ibox)
c     2                    rintl(ilev), reps
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
      nnewboxes=0
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
            call ortho_eval_nd(ndim,nd,norder,
     1           vcoefs(1,1,ibox),norder,vpot(1,1,jbox),
     2           polyvc(1,1,1,j))
c
            call ortho_trans_nd(ndim,nd,0,norder,vpot(1,1,jbox),
     1           vcoefs(1,1,jbox),umat_nd)
c
c           interpolate unm1 too
c           as an initial guess for the nonlinear solve
            call ortho_eval_nd(ndim,nd,norder,
     1           ucoefs1(1,1,ibox),norder,unm1(1,1,jbox),
     2           polyvc(1,1,1,j))
c
            call ortho_trans_nd(ndim,nd,0,norder,unm1(1,1,jbox),
     1           ucoefs1(1,1,jbox),umat_nd)
c
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
c           now get the grid on the box jbox
c           and do the nonlinear solve again for un
            do ii=1, npbox
              do jj=1, ndim
                xy(jj) = grid(jj,ii)*boxsize(jlev)+centers(jj,jbox)
              enddo
c
c           set parameters
c           secant solves:
c           u = g+alpha*f(u,x,t)
c
              msi = 40
              ndamp = 1
              g = vpot(1,ii,jbox)
              u0 = unm1(1,ii,jbox)
              u1 = g+alpha*fn(1,ii,jbox)
              tol = 1.0d-14
c             fix a problem: u0 and u1 are too close
              if(abs(u0-u1) .lt. 10.0d0*tol) then
                u1 = u0+10.0d0*tol
              endif
              call secant(fevaln,tt,xy(1),xy(2),u0,u1,g,alpha,u,fu,tol,
     1             msi,isi,ndamp,ier)
c
c              dpars(1) = tt
c              call rhsfun1(nd,xy,dpars,zpars,ipars,f0)
c              write(227,*) xy(1), xy(2), dpars(1),  
c     1                     f0(1), fn(1,ii,jbox), 
c     1                     abs(f0(1)-fn(1,ii,jbox))
c              call rhsfun2(nd,xy,dpars,zpars,ipars,g0)
c              write(227,*) xy(1), xy(2), dpars(1), g0, g,
c     1                     abs(g0-g)
c              dpars(1)=0.2d0
c              call uexact(nd,xy,dpars,zpars,ipars,u01)
c----------------
c             fn, vpot, unm1 are all interpolated correctly
c             why is this solve still problematic??
c----------------
c              write(227,*) xy(1), xy(2), dpars(1), u01, unm1(1,ii,jbox),
c     1                     abs(u01-unm1(1,ii,jbox))
c
c             I suppose the indices are correct?
              un(1,ii,jbox)=u
ccc           fn = fu? check the calling sequence
ccc           of secant
              fn(1,ii,jbox)=fu
c              write(225,*) tt, xy(1), xy(2), u0, u1, g, alpha,
c     1            ii, jbox
c              write(226,*) xy(1), xy(2), u, fu, isi, ier
            enddo
c
c-------------------
c           now instead of interpolating un and fn
c           we have the function values re-evaluated
c           convert them to coefficients and estimate errors
c
            call ortho_trans_nd(ndim,nd,0,norder,fn(1,1,jbox),
     1           fcoefs(1,1,jbox),umat_nd)
c
            call ortho_trans_nd(ndim,nd,0,norder,un(1,1,jbox),
     1           ucoefs(1,1,jbox),umat_nd)
c
            call fun_err(nd,npbox,fcoefs(1,1,jbox),rmask,
     1           iptype,rscale,erra)
c
            call fun_err(nd,npbox,ucoefs(1,1,jbox),rmask,
     1           iptype,rscale,errb)
c
            erra = erra/rsum
            errb = errb/rsum
c           error obtained for the newly added box jbox
c
            erri = max(erra, errb)
c
            if(erri .gt. reps*rintl(jlev)) then
c             flag the new box
              ifrefine(jbox) = 1
            endif
c            write(224,*) erri, reps*rintl(jlev), ifrefine(jbox)
ccc         refine a few more boxes to trigger coarsening
c            if(jbox .eq. 325) then
c              ifrefine(jbox)=1
c            endif
c            write(*,*) '...', erra,reps, rintl(jlev)
c         end of the loop over mc children
          enddo
c       end refining box ibox
        endif
      enddo
c
c------------------------------------------------
c
      write(*,*) 'after refinement, nlevels=',nlevels
      write(*,*) 'nnewboxes=', nnewboxes
c      stop
c
c     old id of the boxes
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
c
      if (nnewboxes.eq.0) goto 4600
c
      ifpgh = 1
c     modification uvals -> uallvals
c     remember to reorg them all:
c     un, fn, fnm1, fnm2,
c     unm1 and vpot are no longer needed 
ccc   attention: be careful about which arrays are needed!!
c
c     pack together un, fn, fnm1, fnm2
      nd0 = 5
      allocate(uallvals(nd0,npbox,mboxes))
      allocate(uallcoefs(nd0,npbox,mboxes))
      allocate(grad(nd0,ndim,npbox,1))
      allocate(hess(nd0,ndim*(ndim+1)/2,npbox,1))
c
      do ibox = 1, mboxes
        do i=1, npbox
          uallvals(1,i,ibox) = un(1,i,ibox)
          uallvals(2,i,ibox) = unm1(1,i,ibox)
          uallvals(3,i,ibox) = fnm1(1,i,ibox)
          uallvals(4,i,ibox) = fnm2(1,i,ibox)
          uallvals(5,i,ibox) = vpot(1,i,ibox)
c
          uallcoefs(1,i,ibox) = ucoefs(1,i,ibox)
          uallcoefs(2,i,ibox) = ucoefs1(1,i,ibox)
          uallcoefs(3,i,ibox) = fcoefs1(1,i,ibox)
          uallcoefs(4,i,ibox) = fcoefs2(1,i,ibox)
          uallcoefs(5,i,ibox) = vcoefs(1,i,ibox)
        enddo
      enddo
c
c     reorg the tree and the packed arrays: uallvals, uallcoefs
c     with nd0 = 5
      call vol_tree_reorg_laddr(ndim, mboxes, nd0, npbox,
     1    nblock, nboxid, nnewboxes, nboxes0, mlevels,
     2    nnewlevels, nlevels0, centers, itree(iptr(1)),
     3    itree(iptr(2)), itree(iptr(3)), itree(iptr(4)),
     4    itree(iptr(5)), ifpgh, uallvals, uallcoefs, grad, hess)
c
c
 4600 continue
c
      write(*,*) 'after reorg, nlevels=', nlevels
      write(*,*) 'nboxes=', nboxes
c      stop
c
c
c------------------------------------------------
c     latest version:
c     set function values on the tree      
c     call a modified version of coarsening
c     subroutine vol_tree_coarsen_f4?
c     stick with the many copies of the same subroutine thing 
c
c     no need to call set_vals again since the arrays
c     are maintained
c     now we can call the coarsening and reorganizing 
c     routines with the packed arrays
c
      itype=0
      call treedata_trans_nd(ndim,nd0,itype,nlevels,itree,
     1     iptr,boxsize,norder,uallvals,uallcoefs,umat_nd)
c     no need to call this routine?
c     uallcoefs already available?
c
      ifpgh = 1
      allocate(ifdelete(mboxes))
c
c     yet another modified version of vol_tree_coarsen
c     vol_tree_coarsen_f4?
c     now call the original version
c     vol_tree_coarsen
      call vol_tree_coarsen(nd0,ndim,reps,ipoly,norder,npbox,
     1     mboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2     ifpgh,uallvals,uallcoefs,grad,hess,
     3     nblock,nboxid,ndelboxes,ifdelete)
c     I believe the coarsening routine does the
c     interpolation from children to parent
c     not the fault of vol_tree_coarsen as it seems
c     
c
      if(ndelboxes .gt. 0) then
        call fgt_vol_tree_reorg_after_coarsen(ndim,mboxes,nd0,npbox,
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
c------------------------------------------------
c
c     now make it level-restricted again
      call computecoll(ndim,nlevels,nboxes,itree(iptr(1)),
     1    boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2    itree(iptr(5)),iperiod,itree(iptr(6)),itree(iptr(7)))
c
c     modification: we need a modified version again
c     interpolating all the arrays:
c     un, fnm1, fnm2
c
c     call vol_tree_fix_lr_interp
      call vol_tree_fix_lr_interp(ndim,nd0,ipoly,iperiod,norder,npbox,
     1    ifpgh,uallvals,uallcoefs,grad,hess,
     2    mboxes,mlevels,centers,boxsize,nboxes,nlevels,
     3    itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4    itree(iptr(5)),itree(iptr(6)),itree(iptr(7)))
c
      write(*,*) 'after fixing lr:'
      write(*,*) 'nboxes=', nboxes, itree(2*nlevels+2)
c
c     now unpack the packed arrays
c
      do ibox = 1, mboxes
        do i=1, npbox
          un(1,i,ibox) = uallvals(1,i,ibox)
          unm1(1,i,ibox) = uallvals(2,i,ibox)
          fnm1(1,i,ibox) = uallvals(3,i,ibox)
          fnm2(1,i,ibox) = uallvals(4,i,ibox)
          vpot(1,i,ibox) = uallvals(5,i,ibox)
        enddo
      enddo
c
c
      end subroutine





c----------------------------------------------
c     subroutine am2_adap_tree
c----------------------------------------------
c
c     2nd order version of am24_adap_tree
c     which checks the resolution of the grid
c     after solving the nonlinear eqn 
c
c     the interface is more or less the same
c     as that of am24_adap_tree
c     except that it's 2nd order only and 
c     needs to resolve fewer arrays
c
c     iadap = 0: no refinement or coarsening
c     iadap = 1: refinement and coarsening
c     iadap = 2: refinement only 
c
c     when used in the initialization stage
c     we can be more generous and use iadap =2
c     so that we can make sure the grid resolves
c     everything in the initialization stage
c
c----------------------------------------------
c
      subroutine am2_adap_tree(norder, ltree, itree,
     1     iptr, rintl, centers, nlevels, mlevels, boxsize, 
     2     nboxes, mboxes, fevaln, tt, iperiod, eps, alpha,
     3     vpot, unm1, un, iadap)
      implicit real*8 (a-h,o-z)
c     global vars
      integer norder, ltree, nboxes, ndim, nlevels
      integer mboxes, mlevels, nd
      integer iperiod, ntarg
      parameter(ndim=2)
      parameter(nd=1)
      integer itree(ltree), iptr(12)
      real*8 centers(ndim, mboxes), boxsize(0:mlevels)
      real*8 tt, eps, alpha, f0(1)
      real*8 vpot(nd,norder**ndim,mboxes)
      real*8 unm1(nd,norder**ndim,mboxes)
      real*8 un(nd,norder**ndim,mboxes)
      complex *16 zk
c     local vars      
      integer, allocatable:: ipars(:)
c     flags for refinement and coarsening
      integer isgn(ndim,2**ndim)
      integer, allocatable:: nblock(:)
      integer, allocatable:: ilevstart(:)
      integer, allocatable:: nboxid(:), ifdelete(:)
      integer, allocatable:: ifrefine(:)
      integer, allocatable:: irefinelev(:), imaxrefinelev(:)
c     real allocatables
      real*8, allocatable:: dpars(:)
c     need more of these fvals thing?
      real*8, allocatable :: fn(:,:,:)
      real*8, allocatable :: uallvals(:,:,:), uallcoefs(:,:,:)
      real*8, allocatable :: fcoefs(:,:,:), ucoefs(:,:,:)
      real*8, allocatable :: ucoefs1(:,:,:), vcoefs(:,:,:)
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8, allocatable :: rmask(:), umat_nd(:,:,:)
      real *8 rintl(0:mlevels)
      real *8 xy(ndim)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
c     change the following names
      real *8, allocatable:: polyvc(:,:,:,:)
c     complex allocatables      
      complex*16, allocatable:: zpars(:)
      external fevaln
c      open(337,file='del1.txt')
c      open(223,file='del2.txt')
c      open(338,file='del3.txt')
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
      mc=2**ndim
      mnbors=3**ndim
      npols = norder**ndim
      npbox = norder**ndim
      ipoly=1
      iptype=2
c
c      write(*,*) npols, npbox, ndim, nd
c
      allocate(dpars(1000))
      allocate(zpars(1000))
      allocate(ipars(1000))
      allocate(fn(nd,npbox,mboxes))
c     allocate coeffs for them all:
c     un, fn, unm1, vpot
      allocate(ucoefs(nd,npbox,mboxes))
      allocate(fcoefs(nd,npbox,mboxes))
      allocate(ucoefs1(nd,npbox,mboxes))
      allocate(vcoefs(nd,npbox,mboxes))
c
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
      allocate(umat_nd(norder,norder,ndim))
      allocate(polyvc(norder,norder,ndim,mc))
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
c ---------------debugging------------------------------
      do ilevel=0, nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
c             a leaf node
            do i=1, npbox
             write(337,*)ibox,i,un(1,i,ibox),unm1(1,i,ibox)
            enddo
          endif
        enddo
      enddo
c-------------------------------------------------------

c
c
c--------------------------------------------------
c     on the currect tree
c     eval fun on the grid of leaf nodes -> fvals
      do ilev=0, nlevels
        do ibox = itree(2*ilev+1), itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild .eq. 0) then
            do i=1,npbox
              do j=1,ndim
                xy(j) = grid(j,i)*boxsize(ilev)+centers(j,ibox)
              enddo
c             attention: check the calling sequence
              call fevaln(tt,xy(1),xy(2),un(1,i,ibox),
     1             fn(1,i,ibox))
            enddo
          endif
        enddo
      enddo
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
             do i=1, npbox
               write(301,*) i, ibox, fcoefs(1,i,ibox),ucoefs(1,i,ibox),
     1                      fn(1,i,ibox),un(1,i,ibox)
             enddo
             erra = erra/rsum
             errb = errb/rsum
c            refine if one of them is greater than the tol
             erri = max(erra, errb)
             if(erri.gt.reps*rintl(ilev)) then
                ifrefine(ibox)=1
             endif
c
             write(223,*) ibox, erra, errb, erri, reps*rintl(ilev), 
     1                    ifrefine(ibox)
c     2                    rintl(ilev), reps
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
      nnewboxes=0
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
            call ortho_eval_nd(ndim,nd,norder,
     1           vcoefs(1,1,ibox),norder,vpot(1,1,jbox),
     2           polyvc(1,1,1,j))
c
            call ortho_trans_nd(ndim,nd,0,norder,vpot(1,1,jbox),
     1           vcoefs(1,1,jbox),umat_nd)
c
c           interpolate unm1 too
c           as an initial guess for the nonlinear solve
            call ortho_eval_nd(ndim,nd,norder,
     1           ucoefs1(1,1,ibox),norder,unm1(1,1,jbox),
     2           polyvc(1,1,1,j))
c
            call ortho_trans_nd(ndim,nd,0,norder,unm1(1,1,jbox),
     1           ucoefs1(1,1,jbox),umat_nd)
c
c           interpolate fn too
c           as an initial guess for the nonlinear solve
            call ortho_eval_nd(ndim,nd,norder,
     1           fcoefs(1,1,ibox),norder,fn(1,1,jbox),
     2           polyvc(1,1,1,j))
c
            call ortho_trans_nd(ndim,nd,0,norder,fn(1,1,jbox),
     1           fcoefs(1,1,jbox),umat_nd)
c
c           now get the grid on the box jbox
c           and do the nonlinear solve again for un
            do ii=1, npbox
              do jj=1, ndim
                xy(jj) = grid(jj,ii)*boxsize(jlev)+centers(jj,jbox)
              enddo
c
c           set parameters
c           secant solves:
c           u = g+alpha*f(u,x,t)
c
              msi = 40
              ndamp = 1
              g = vpot(1,ii,jbox)
              u0 = unm1(1,ii,jbox)
              u1 = g+alpha*fn(1,ii,jbox)
              tol = 1.0d-14
c             fix a problem: u0 and u1 are too close
              if(abs(u0-u1) .lt. 10.0d0*tol) then
                u1 = u0+10.0d0*tol
              endif
              call secant(fevaln,tt,xy(1),xy(2),u0,u1,g,alpha,u,fu,tol,
     1             msi,isi,ndamp,ier)
c
c              dpars(1) = tt
c              call rhsfun1(nd,xy,dpars,zpars,ipars,f0)
c              write(227,*) xy(1), xy(2), dpars(1),  
c     1                     f0(1), fn(1,ii,jbox), 
c     1                     abs(f0(1)-fn(1,ii,jbox))
c              call rhsfun2(nd,xy,dpars,zpars,ipars,g0)
c              write(227,*) xy(1), xy(2), dpars(1), g0, g,
c     1                     abs(g0-g)
c              dpars(1)=0.2d0
c              call uexact(nd,xy,dpars,zpars,ipars,u01)
c----------------
c             fn, vpot, unm1 are all interpolated correctly
c             why is this solve still problematic??
c----------------
c              write(227,*) xy(1), xy(2), dpars(1), u01, unm1(1,ii,jbox),
c     1                     abs(u01-unm1(1,ii,jbox))
c
c             I suppose the indices are correct?
              un(1,ii,jbox)=u
ccc           fn = fu? check the calling sequence
ccc           of secant
              fn(1,ii,jbox)=fu
c              write(225,*) tt, xy(1), xy(2), u0, u1, g, alpha,
c     1            ii, jbox
c              write(226,*) xy(1), xy(2), u, fu, isi, ier
            enddo
c
c-------------------
c           now instead of interpolating un and fn
c           we have the function values re-evaluated
c           convert them to coefficients and estimate errors
c
            call ortho_trans_nd(ndim,nd,0,norder,fn(1,1,jbox),
     1           fcoefs(1,1,jbox),umat_nd)
c
            call ortho_trans_nd(ndim,nd,0,norder,un(1,1,jbox),
     1           ucoefs(1,1,jbox),umat_nd)
c
            call fun_err(nd,npbox,fcoefs(1,1,jbox),rmask,
     1           iptype,rscale,erra)
c
            call fun_err(nd,npbox,ucoefs(1,1,jbox),rmask,
     1           iptype,rscale,errb)
c
            erra = erra/rsum
            errb = errb/rsum
c           error obtained for the newly added box jbox
c
            erri = max(erra, errb)
c
            if(erri .gt. reps*rintl(jlev)) then
c             flag the new box
              ifrefine(jbox) = 1
            endif
            write(224,*) erra, errb, reps*rintl(jlev), ifrefine(jbox)
c           still the problem of rintl
ccc         refine a few more boxes to trigger coarsening
c            if(jbox .eq. 325) then
c              ifrefine(jbox)=1
c            endif
c            write(*,*) '...', erra,reps, rintl(jlev)
c         end of the loop over mc children
          enddo
c       end refining box ibox
        endif
      enddo
c     end of the iadap if
      endif
c
c      stop
c------------------------------------------------
c
      write(*,*) 'after refinement, nlevels=',nlevels
      write(*,*) 'nnewboxes=', nnewboxes
c      stop
c


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
      nd0 = 2
      allocate(uallvals(nd0,npbox,mboxes))
      allocate(uallcoefs(nd0,npbox,mboxes))
      allocate(grad(nd0,ndim,npbox,1))
      allocate(hess(nd0,ndim*(ndim+1)/2,npbox,1))
c
c      allocate(grad(nd,ndim,npbox,1))
c      allocate(hess(nd,ndim*(ndim+1)/2,npbox,1))
c
c
c     do we need fn? not sure. do it anyway.
c     yes we do need fn
c     we need to make sure both fn and un are resolved
      do ibox = 1, mboxes
        do i=1, npbox
          uallvals(1,i,ibox) = un(1,i,ibox)
          uallvals(2,i,ibox) = fn(1,i,ibox)
c
          uallcoefs(1,i,ibox) = ucoefs(1,i,ibox)
          uallcoefs(2,i,ibox) = fcoefs(1,i,ibox)
        enddo
      enddo

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
c ---------------debugging------------------------------
      do ilevel=0, nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
c             a leaf node
            do i=1, npbox
             write(338,*)ibox,i,un(1,i,ibox),unm1(1,i,ibox)
            enddo
          endif
        enddo
      enddo
c-------------------------------------------------------

c      stop
c
c------------------------------------------------
c
      if(iadap .eq. 1) then
c     when iadap = 1, do the coarsening
c
      itype=0
      call treedata_trans_nd(ndim,nd,itype,nlevels,itree,
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
      call vol_tree_coarsen(nd0,ndim,reps,ipoly,norder,npbox,
     1     mboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2     ifpgh,uallvals,uallcoefs,grad,hess,
     3     nblock,nboxid,ndelboxes,ifdelete)
c     I believe the coarsening routine does the
c     interpolation from children to parent
c     not the fault of vol_tree_coarsen as it seems
c
      write(*,*) 'ndelboxes=', ndelboxes
c     
      if(ndelboxes .gt. 0) then
        call fgt_vol_tree_reorg_after_coarsen(ndim,mboxes,nd0,npbox,
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
c
c------------------------------------------------
c
c     now compute the colls and make it 
c     level-restricted again
c
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
c
c
c     now unpack the packed arrays
c
      do ibox = 1, mboxes
        do i=1, npbox
          un(1,i,ibox) = uallvals(1,i,ibox)
          fn(1,i,ibox) = uallvals(2,i,ibox)
        enddo
      enddo
c
c------------------------------------------------
c
      end subroutine
