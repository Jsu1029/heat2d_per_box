c     Navier-Stokes solver adaptive version
c     periodic boundary condition
c     impressive + viscocity
c     boxcode [-0.5,0.5]^2
c     2d
c
c     Input
c     fforce: right hand funtion
c     finit: initial function
c
c     nordert: convergence order
c        2: trapezoidal rule for the integral about time
c           2nd-order BDF for the linear terms
c        4: simpson's rule
c           4th-order BDF for the linear terms
c     visc: viscocity
c     dt: delta t
c     ntot: time steps
c     eps: error tolerance
c
c     output:
c     velocity: solution of the equation vector
c     vorticity:vorticity of velocity
c
c
c     Iniatialization
c     M_1(deltat) 1st method
c     u_0,u_1,u_2,u_3 Richardson extrapolation
c     order2 -- 2 M_1(deltat/2) - M_1(deltat)
c     order4 -- 8. M_1(deltat/8) - 8/3 M_1(deltat/4)
c                + 2/3 M_1(deltat/2) - 1/21 M_1(deltat)
c
c--------------------------------------------------
c
      subroutine ns2dperfm4(nd,norder, finit,
     1           fforce,dt,ntot, eps, mxltree,
     2           mxboxes,mxlevels,ltree, itree, iptr, centers,
     3           nlevels,boxsize, nboxes, ifnewtree,ntarg,
     4           nx,ny,velocity,vorticity,velocityp)
      implicit real*8 (a-h,o-z)
c     input:
      integer nd,norder,tot
      integer mxltree, mxboxes,mxlevels
      integer ntot,ntarg,nx,ny,ifnewtree
      real *8 dt,eps
c     output
      integer ltree
      real*8 centers(ndim,mxltree), boxsize(0:mxlevels)
      integer itree(mxltree)
      integer nlevels, iptr(8),nboxes
      real*8
      
c
      integer ipars(100)
      real *8 dpars(1000)
      complex*16 zpars(10), zk


      character *12 fname1
      character *8 fname2
      character *9 fname3
      data ima/(0.0d0,1.0d0)/
      character*1 type
      external finit,fforce
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
c      nd=1
c
      epstree=1.0d-1*eps
c      epstree=1.0d-4*eps
      eta=1.0d0
      zk=1.0d0
c
      xsize0=1.0d0
      cent0(1)=0.0d0
      cent0(2)=0.0d0
c
      if(nordert .eq. 2) then
        nnodes=2
      elseif(nordert .eq. 4) then
        nnodes=3
      endif
c

c
      allocate(umat_nd(norder,norder,ndim))
      allocate(xq(norder),umat(norder,norder),
     1    vmat(norder,norder))
      allocate(usol_pre(nordert*2,npbox,mxboxes))
c
c
      call mktreef2d(levelbox,icolbox,irowbox,nboxes,nlev,
     1          iparentbox, ichildbox,
     2          nblevel, iboxlev, istartlev,
     3          mxboxes, itemparray, mxlevels, eps,
     4          finit1,finit2,ndeg,cent0,xsize0)
c
      if(nboxes .gt. mxboxes)then
       write(*,*)'there are more boxes needed than are allotted'
       write(*,*)'for in the memory space.'
       write(*,*)'readjust the memory.'
       write(*,*)'i am stopping'
       stop
      endif
c
      ifixflag=0
      call restriction(levelbox,iparentbox,ichildbox,icolbox,
     1     irowbox,icolleagbox,nboxes,nlev,
     2     nblevel,iboxlev,istartlev,iperiod,ifixflag)
c
      if(ifixflag .eq. 1)then
c     the correction routine is used only in the adaptive case:
         write(*,*)'correcting tree'
         call fixtree(levelbox,iparentbox,ichildbox,icolbox,
     1        irowbox,icolleagbox,nboxes,nlev,
     2        nblevel, iboxlev, istartlev,iperiod,
     3        iflag, mxboxes,itemparray)
      
      endif
c
c     test to see if the total number of boxes
c     is over the allotted memory space.
      if(nboxes .gt. mxboxes)then
       write(*,*)'there are more boxes needed than are allotted'
       write(*,*)'for in the memory space.'
       write(*,*)'readjust the memory.'
       write(*,*)'i am stopping'
       stop
      endif
       
      write(*,*)'nboxes = ',nboxes
      write(*,*)'nlev= ',nlev
c
      iprint=-101
       call printtree(levelbox,icolbox,irowbox,nboxes,
     1       nlev,iparentbox,ichildbox,nblevel,iboxlev,
     2       istartlev,iprint,nleaves,cent0,xsize0)
      write(*,*) 'nleaves=',nleaves
c
c
c     allocate mem
      allocate(fvals(nd,npbox,mxboxes))
c     the initial potential and volume potential
      allocate(vpot(nd,npbox,mxboxes))
c
      allocate(valueinit(nd,npbox,mxboxes))
      allocate(potin(nd,npbox,mxboxes))
      allocate(potinp(nd,ntarg))
c
c     what is ntarg though
      allocate(targs(ndim,ntarg))
      allocate(potp(ntarg))
      allocate(vpotp(nd,ntarg))
c
c     why do we need these?
      allocate(targ(ndim))
      allocate(xref(ndim,npbox),wts(npbox))
c     extra targets
      allocate(coefsp(nd,npbox,mxboxes))
c
      allocate(fs(nnodes*2,npbox,mxboxes))
      
c
      allocate(timeinfo(100))
c
      allocate(grad(nd,ndim,npbox,mxboxes))
      allocate(hess(nd,nhess,npbox,mxboxes))
      allocate(gradp(nd,ndim,ntarg))
      allocate(hessp(nd,nhess,ntarg))


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
      itype = 0
      type='f'
      call polytens_exps_nd(ndim,ipoly,itype,norder,type,xref,
     1    utmp,1,vtmp,1,wts)
c
      call setf(ndeg,valueinit1,icolbox,irowbox,ichildbox,
     1    nlev,nblevel,iboxlev,istartlev,finit1,
     2    cent0,xsize0)
      
      call setf(ndeg,valueinit2,icolbox,irowbox,ichildbox,
     1    nlev,nblevel,iboxlev,istartlev,finit2,
     2    cent0,xsize0)

c
        do j=1, nboxes
        do i=1, npbox
           usol_pre(1,i,j)=valueinit1(i,j)
           usol_pre(2,i,j)=valueinit2(i,j)
        enddo
        enddo

c Initialization using Richradson extrapolation

c

c     1-order solver

      dt_init = dt/100.0d0
      ntot_init = (nordert-1)*100
c     i. compute the initial potential
        ifvtarg=1
        ifptarg=1
        ifpghtarg=1

      do it = 1,ntot_init

        tt=it*dt_init
        tpre=(it-1)*dt_init
c
      write(*,*)'runrunrun???'
      call setf8(fright,icolbox, irowbox, ichildbox,
     1          nlev, nblevel, iboxlev, istartlev,fforce1,tpre)
c
      call setf8(fright2,icolbox,irowbox,ichildbox,
     1           nlev,nblevel,iboxlev,istartlev,fforce2,tpre)

c compute fs(t_{n-1})
       do i = 0,nlev
         do j = istartlev(i),istartlev(i)+nblevel(i)-1
           ibox = iboxlev(j)
           if(ichildbox(1,ibox) .lt. 0)then
              do l = 1,64
                convective1(l,ibox) = fright(l,ibox)-
     1     valueinit1(l,ibox)*frightinitialx1(l,ibox)-
     2     valueinit2(l,ibox)*frightinitialy1(l,ibox)

                convective2(l,ibox) = fright2(l,ibox)-
     1     valueinit1(l,ibox)*frightinitialx2(l,ibox)-
     2     valueinit2(l,ibox)*frightinitialy2(l,ibox)

              valueinit(1,l,ibox) = valueinit1(l,ibox)
              valueinit(2,l,ibox) = valueinit2(l,ibox)
              end do
           end if
         end do
       end do
c
      call adaptfmm2d8(ndeg, nlev, levelbox, iparentbox,
     1     ichildbox, icolbox, irowbox, nboxes, nblevel,
     2     iboxlev, istartlev, convective1, iprec, iperiod,
     3     pot, ttotal1)

      call getderivatives8(pot,
     1      nlev, ichildbox, nblevel, iboxlev, istartlev,
     2      potx, poty, potxy,potxx,potyy, lap, coeffsu5)


      call adaptfmm2d8(ndeg,nlev,levelbox,iparentbox,
     1        ichildbox,icolbox, irowbox, nboxes, nblevel,
     2        iboxlev, istartlev, convective2,iprec,iperiod,
     3        pot2,ttotal2)
  
      call getderivatives8(pot2,
     1     nlev,ichildbox,nblevel,iboxlev,istartlev,
     2     potx2,poty2,potxy2,potxx2,potyy2,lap2,coeffsu6)
c
      do i = 0,nlev
        do j = istartlev(i),istartlev(i)+nblevel(i)-1
          ibox = iboxlev(j)
          if (ichildbox(1,ibox) .gt. 0)goto 2200
          do l = 1,64
            frightG1(l,ibox) = potxx(l,ibox)+potxy2(l,ibox)
            frightG2(l,ibox) = potxy(l,ibox)+potyy2(l,ibox)

            frightS(1,l,ibox) = convective1(l,ibox) - frightG1(l,ibox)
            frightS(1,l,ibox) = convective2(l,ibox) - frightG2(l,ibox)
          end do
 2200  continue
       end do
      end do

c  oletree to newtree
      call oldtree2newtree2d(nlev,levelbox,iparentbox,
     1    ichildbox,icolbox,irowbox,nboxes,nblevel,
     2    iboxlev,istartlev,cent0,xsize0,iperiod,
     3    ltree,nlevels,itree,iptr,centers,boxsize)
c     initialize vpot to zero
c
      do i=1,nboxes
        do k=1,npbox
          do l = 1,nd
          vpot(l,k,i)=0.0d0
          enddo
        enddo
      enddo

c   Volume potential
        delta=dt_init

        call boxfgt(nd,ndim,delta,eps0,ipoly,iperiod,norder,npbox,
     1       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,
     2       frightS,ifpgh,vpot,grad,hess,ifnewtree,ntarg,targs,
     3       ifpghtarg,vpotp,gradp,hessp,timeinfo)
 
c
        do i=1, nboxes
          if(itree(iptr(4)+i-1).eq.0) then
            do j=1, npbox
               do l=1,nd
               vpot(l,j,i)=vpot(l,j,i)/(4.0d0*pi)
     1                  /(dt_init)*dt_init/2.0d0
               enddo
            enddo
          endif
        enddo

c

      call setf8(fright,icolbox, irowbox, ichildbox,
     1          nlev, nblevel, iboxlev, istartlev,fforce1,tt)
c
      call setf8(fright2,icolbox,irowbox,ichildbox,
     1           nlev,nblevel,iboxlev,istartlev,fforce2,tt)

c compute fs(t_{n})
       do i = 0,nlev
         do j = istartlev(i),istartlev(i)+nblevel(i)-1
           ibox = iboxlev(j)
           if(ichildbox(1,ibox) .lt. 0)then
              do l = 1,64
                convective1(l,ibox) = fright(l,ibox)-
     1     valueinit1(l,ibox)*frightinitialx1(l,ibox)-
     2     valueinit2(l,ibox)*frightinitialy1(l,ibox)

                convective2(l,ibox) = fright2(l,ibox)-
     1     valueinit1(l,ibox)*frightinitialx2(l,ibox)-
     2     valueinit2(l,ibox)*frightinitialy2(l,ibox)
              end do
           end if
         end do
       end do
c
      call adaptfmm2d8(ndeg, nlev, levelbox, iparentbox,
     1     ichildbox, icolbox, irowbox, nboxes, nblevel,
     2     iboxlev, istartlev, convective1, iprec, iperiod,
     3     pot, ttotal1)

      call getderivatives8(pot,
     1      nlev, ichildbox, nblevel, iboxlev, istartlev,
     2      potx, poty, potxy,potxx,potyy, lap, coeffsu5)


      call adaptfmm2d8(ndeg,nlev,levelbox,iparentbox,
     1        ichildbox,icolbox, irowbox, nboxes, nblevel,
     2        iboxlev, istartlev, convective2,iprec,iperiod,
     3        pot2,ttotal2)
  
      call getderivatives8(pot2,
     1     nlev,ichildbox,nblevel,iboxlev,istartlev,
     2     potx2,poty2,potxy2,potxx2,potyy2,lap2,coeffsu6)
c
      do i = 0,nlev
        do j = istartlev(i),istartlev(i)+nblevel(i)-1
          ibox = iboxlev(j)
          if (ichildbox(1,ibox) .gt. 0)goto 2300
          do l = 1,64
            frightG1(l,ibox) = potxx(l,ibox)+potxy2(l,ibox)
            frightG2(l,ibox) = potxy(l,ibox)+potyy2(l,ibox)

            frightS(1,l,ibox) = convective1(l,ibox) - frightG1(l,ibox)
            frightS(1,l,ibox) = convective2(l,ibox) - frightG2(l,ibox)
          end do
 2300  continue
       end do
      end do

      do ilevel=0,nlevels
        bs = boxsize(ilevel)/2.0d0
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
c            a leaf node
c            get the grid and compute the forcing term
             do j=1,npbox
               vpot(1,j,ibox)=vpot(1,j,ibox)+
     1                 frightS(1,j,ibox)*dt_init/2.0d0
               vpot(2,j,ibox)=vpot(2,j,ibox)+
     1                 frightS(2,j,ibox)*dt_init/2.0d0
             enddo
          endif
        enddo
      enddo

c
      call initpot2(nd,norder, ltree, itree, iptr, centers,
     1           nlevels, boxsize, nboxes, valueinit, dt_init,
     2           iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           eps, potin, potinp)


c
c     iii. add the potentials
        do i=1, nboxes
          do j=1, npbox
            valueinit1(j,i)=potin(1,j,i)+vpot(1,j,i)
            valueinit2(j,i)=potin(2,j,i)+vpot(2,j,i)
          enddo
        enddo
c

c    Initilization
c
        if (it .eq. 100)then
        do j=1, nboxes
        do i=1, npbox
           usol_pre(3,i,j)=valueinit1(i,j)
           usol_pre(4,i,j)=valueinit2(i,j)
        enddo
        enddo
        endif
c
        if (it .eq. 200)then
        do j=1, nboxes
        do i=1, npbox
           usol_pre(5,i,j)=valueinit1(i,j)
           usol_pre(6,i,j)=valueinit2(i,j)
        enddo
        enddo
        endif
c
        if (it .eq. 300)then
        do j=1, nboxes
        do i=1, npbox
           usol_pre(7,i,j)=valueinit1(i,j)
           usol_pre(8,i,j)=valueinit2(i,j)
        enddo
        enddo
        endif

c
      enddo


c
      do it=nordert, ntot
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
          dpars(1)=tt
c
         call adaptreef2d(mxboxes, mxlevel, levelbox,
     1           icolbox, irowbox, nboxes, nlev, iparentbox,
     2           ichildbox, nblevel, iboxlev, istartlev,
     3           ndeg,valueinit1,valueinit2,fforce1,
     4           fforce2, tt,
     5           eps, cent0, xsize0)
c
        ifixflag=0
        call restriction(levelbox, iparentbox, ichildbox,
     1       icolbox, irowbox, icolleagbox, nboxes, nlev,
     2       nblevel, iboxlev, istartlev, iperiod, ifixflag)

c       fix the tree if needed
        if(ifixflag .eq. 1) then
          write(*,*)'correcting tree'
c
          call fixtreenf2(levelbox,iparentbox,ichildbox,
     1         icolbox,irowbox,icolleagbox,nboxes,nlev,
     2         nblevel,iboxlev,istartlev,iperiod,iflag,
     3         mxboxes,itemparray,ndeg,valueinit1,valueinit2)

        endif

        write(*,*) 'nboxes=',nboxes
        write(*,*)'nlev= ',nlev
c
       iprint=-100
       call printtree(levelbox,icolbox,irowbox,nboxes,
     1       nlev,iparentbox,ichildbox,nblevel,iboxlev,
     2       istartlev,iprint,nleaves,cent0,xsize0)
      
        write(*,*) 'nleaves=',nleaves


        endif

        if (nordert .eq. 4 )then
         call compute_fs4(ndeg,nordert,nlev,levelbox,iparentbox,
     1        ichildbox,icolbox, irowbox, nboxes, nblevel,
     2        iboxlev, istartlev,iprec,iperiod, nnodes,
     3        fforce1,fforce2,dt,tt,nnoeds,usol_pre,
     4        fs,fsp_now1,fsp_now2,norder,ndim)
        elseif (nordert .eq. 2)then
          call compute_fs2(ndeg,nordert,nlev,levelbox,iparentbox,
     1        ichildbox,icolbox, irowbox, nboxes, nblevel,
     2        iboxlev, istartlev,iprec,iperiod, nnodes,
     3        fforce1,fforce2,dt,tt,nnoeds,usol_pre,ntarg
     4        fs,fsp_now1,fsp_now2,norder,ndim)

        endif
 
        do j=1,ntarg
           potinp(1,j)=0.0d0
           potinp(2,j)=0.0d0
        enddo
c
c     i. compute the initial potential
        ifvtarg=1
        ifptarg=1
c       might need to modify this interface
c       add in mxboxes, mlevels?
c       I believe it's not a problem
        write(*,*) 'calling initpot'
c     Notice dt!
c     I[u(t_n)]
        call initpot2(nd,norder, ltree, itree, iptr, centers,
     1       nlevels, boxsize, nboxes, valueinit, dt,
     2       iperiod, ntarg, targs, ifvtarg, ifptarg,
     3       eps, potin, potinp)
c
c     ii. compute the volume potential
c       might need to modify this interface
c       add in mxboxes, mlevels
        write(*,*) 'calling volpot'
        call volpot2(norder, nordert, ltree, itree, iptr,
     4       centers, nlevels, boxsize, nboxes, frightS, tpre,
     2       dt, iperiod, ntarg, targs, ifvtarg, ifptarg,
     3       eps, vpot, vpotp)
c
c     iii. add the potentials
        do i=1, nboxes
          do j=1, npbox
            usol(1,j,i)=potin(1,j,i)+vpot(1,j,i)
            usol(2,j,i)=potin(2,j,i)+vpot(2,j,i)
            velocity1(j,i)=usol(1,j,i)
            velocity2(j,i)=usol(2,j,i)
          enddo
        enddo
c
       do i = 1,ntarg
         usolp(1,i)=potinp(1,i)+vpotp(1,i)
         usolp(2,i)=potinp(2,i)+vpotp(2,i)
         velocityp1(i)=usolp(1,i)
         velocityp2(i)=usolp(2,i)
       enddo
c    Initilization
c
        do j=1, nboxes
        do i=1, npbox
           do k = 1,nordert-1
            usol_pre(2*k-1,i,j) = usol_pre(2*k+1,i,j)
            usol_pre(2*k,i,j) = usol_pre(2*k+2,i,j)
           enddo
        enddo
        enddo
c
        do j=1, nboxes
        do i=1, npbox
           usol_pre(nordert*2-1,i,j)=velocity1(i,j)
           usol_pre(nordert*2,i,j)=velocity2(i,j)
        enddo
        enddo
c
      enddo

c
      end subroutine









c  compute helomholtz decomposition FS (nordert = 4)
c  we need FS(t_n) FS(t_{n+1/2}) FS(t_{n+1})
c  then need u_{tn} u_{t_{n+1/2}} u_{t_{n+1}}
c  the latter two term using BDF method
      subroutine compute_fs4(ndeg,nordert,nlev,levelbox,iparentbox,
     1        ichildbox,icolbox, irowbox, nboxes, nblevel,
     2        iboxlev, istartlev,iprec,iperiod, nnodes,
     3        fforce1,fforce2,dt,tt,nnoeds,usol_pre,
     4        fs,fsp_now1,fsp_now2,ntarg,targs,norder,ndim)
      implicit real*8 (a-h,o-z)
c
      integer nlev,ndeg,nordert,ntarg,norder,ndim
      integer levelbox(nboxes), iparentbox(nboxes)
      integer ichildbox(4,nboxes), icolleagbox(9,nboxes)
      integer icolbox(nboxes), irowbox(nboxes)
      integer nblevel(0:nlev), iboxlev(nboxes)
      integer itemparray(nboxes)
      integer istartlev(0:nlev)
      integer iflag(nboxes)
      integer noodes
      real *8 targs(nd,ntarg)
      real *8 xg(ntarg)
      real *8 yg(ntarg)
      real *8 tt

      real *8 fright(64,nboxes),fright2(64,nboxes)
      real *8 usol_pre(nordert*2,ndeg**2,nboxes)
c
      real *8 usol_prex1(ndeg**2,nboxes)
c
      real *8 usol_prey1(ndeg**2,nboxes)
c
      real *8 usol_prex2(ndeg**2,nboxes)
c
      real *8 usol_prey2(ndeg**2,nboxes)
c
      real *8 usol_prex3(ndeg**2,nboxes)
c
      real *8 usol_prey3(ndeg**2,nboxes)
c
      real *8 usol_prex4(ndeg**2,nboxes)
c
      real *8 usol_prey4(ndeg**2,nboxes)
c
      real *8 usol_half1(ndeg**2,nboxes)
      real *8 usol_now1(ndeg**2,nboxes)
c
      real *8 usol_half2(ndeg**2,nboxes)
      real *8 usol_now2(ndeg**2,nboxes)

      real*8 convective1(64,nboxes),convective2(64,nboxes)

      real *8 cent0(2), xsize0
c     output
      real *8 fs(nnodes*2,ndeg**2,nboxes)
c
c     initial value's derivates
      real *8  frightinitialx1(64,nboxes)
      real *8  frightinitialy1(64,nboxes)
      real *8  frightinitialx2(64,nboxes)
      real *8  frightinitialy2(64,nboxes)
      real *8  frightinitialxy1(64,nboxes)
      real *8  frightinitialxx1(64,nboxes)
      real *8  frightinitialyy1(64,nboxes)
      real *8  frightinitiallap1(64,nboxes)
      real *8  frightinitialxy2(64,nboxes)
      real *8  frightinitialxx2(64,nboxes)
      real *8  frightinitialyy2(64,nboxes)
      real *8  frightinitiallap2(64,nboxes)
c
      real *8  coeffsu(0:7,0:7,nboxes)
      real *8  coeffsu2(0:7,0:7,nboxes)
      real *8  coeffsu3(0:7,0:7,nboxes)
      real *8  coeffsu4(0:7,0:7,nboxes)
      real *8  coeffsu5(0:7,0:7,nboxes)
      real *8  coeffsu6(0:7,0:7,nboxes)
      real *8  coeffsu7(0:7,0:7,nboxes)
      real *8  coeffsu8(0:7,0:7,nboxes)

c
      real *8  fmmpot(64,nboxes)
      real *8  fmmpotx(64,nboxes), fmmpoty(64,nboxes)
      real *8  fmmpotxy(64,nboxes)
      real *8  fmmpotxx(64,nboxes),fmmpotyy(64,nboxes),fmmlap(64,nboxes)
      real *8  fmmpot2(64,nboxes)
      real *8  fmmpotx2(64,nboxes),fmmpoty2(64,nboxes)
      real *8  fmmpotxy2(64,nboxes)
      real *8  fmmpotxx2(64,nboxes),fmmpotyy2(64,nboxes)
      real *8  fmmlap2(64,nboxes)

c    helmholtz
      real *8  frightG1(norder**ndim,nboxes)
      real *8  frightG2(norder**ndim,nboxes)
      real *8  frightS1(norder**ndim,nboxes)
      real *8  frightS2(norder**ndim,nboxes)
c    fs at targets
      real *8  fsp_now1(ntarg)
      real *8  fsp_now2(ntarg)
c
      external fforce1,fforce2
c
      if (nordert .eq. 4)then
        do j = 1,nboxes
           do i = 1,64
           usol_prex1(i,j) = usol_pre(1,i,j)
           usol_prey1(i,j) = usol_pre(2,i,j)
           usol_prex2(i,j) = usol_pre(3,i,j)
           usol_prey2(i,j) = usol_pre(4,i,j)
           usol_prex3(i,j) = usol_pre(5,i,j)
           usol_prey3(i,j) = usol_pre(6,i,j)
           usol_prex4(i,j) = usol_pre(7,i,j)
           usol_prey4(i,j) = usol_pre(8,i,j)
c  compute u_{n+1/2}  nordert=4
           usol_half1(i,j) = 35.0d0/16.0d0*usol_prex1(i,j)-
     1    35.0d0/16.0d0*usol_prex2(i,j)+21.0d0/16.0d0*
     2     usol_prex3(i,j)-5.0d0/16.0d0*usol_prex4(i,j)

           usol_half2(i,j) = 35.0d0/16.0d0*usol_prey1(i,j)-
     1    35.0d0/16.0d0*usol_prey2(i,j)+21.0d0/16.0d0*
     2     usol_prey3(i,j)-5.0d0/16.0d0*usol_prey4(i,j)
c  compute u_{n+1} nordert=4
           usol_now1(i,j) = 4.0d0*usol_prex1(i,j)-
     1    6.0d0*usol_prex2(i,j)+4.0d0*
     2     usol_prex3(i,j)-1.0d0*usol_prex4(i,j)

           usol_now2(i,j) = 4.0d0*usol_prey1(i,j)-
     1    6.0d0*usol_prey2(i,j)+4.0d0*
     2     usol_prey3(i,j)-1.0d0*usol_prey4(i,j)
    
           enddo
        enddo
      endif


      t_pre=tt-dt
      t_half = tt-dt/2.0d0
      t_now = tt
c
      do j = 1,nboxes
        do i = 1,64
          fright(i,j) = 0.0d0
          fright2(i,j) = 0.0d0
          convective1(i,j)=0.0d0
          convective2(i,j)=0.0d0
          fmmpot(i,j) = 0.0d0
          fmmpot2(i,j) = 0.0d0
        enddo
       enddo

c     compute fs(t_n)
c
      call getderivatives8(usol_prex1,
     1      nlev, ichildbox, nblevel, iboxlev, istartlev,
     2      frightinitialx1, frightinitialy1, frightinitialxy1,
     3      frightinitialxx1,frightinitialyy1,frightinitiallap1,
     4       coeffsu)
c
      call getderivatives8(usol_prey1,
     1      nlev, ichildbox, nblevel, iboxlev, istartlev,
     2      frightinitialx2, frightinitialy2, frightinitialxy2,
     3      frightinitialxx2,frightinitialyy2,frightinitiallap2,
     4       coeffsu2)
c
      call setf8(fright,icolbox, irowbox, ichildbox,
     1          nlev, nblevel, iboxlev, istartlev,fforce1,t_pre)
c
      call setf8(fright2,icolbox,irowbox,ichildbox,
     1           nlev,nblevel,iboxlev,istartlev,fforce2,t_pre)
c     convective term of ns eq
c
        do i=1, nboxes
        do j=1, npbox
                convective1(j,i) = fright(j,i)-
     1     usol_prex1(j,i)*frightinitialx1(j,i)-
     2     usol_prey1(j,i)*frightinitialy1(j,i)

                convective2(j,i) = fright2(j,i)-
     1     usol_prex1(j,i)*frightinitialx2(j,i)-
     2     usol_prey1(j,i)*frightinitialy2(j,i)
        enddo
        enddo
c
      call adaptfmm2d8(ndeg, nlev, levelbox, iparentbox,
     1     ichildbox, icolbox, irowbox, nboxes, nblevel,
     2     iboxlev, istartlev, convective1, iprec, iperiod,
     3     fmmpot, ttotal1)

      call getderivatives8(pot,
     1      nlev, ichildbox, nblevel, iboxlev, istartlev,
     2      fmmpotx, fmmpoty, fmmpotxy,fmmpotxx,
     3      fmmpotyy, fmmlap, coeffsu5)


      call adaptfmm2d8(ndeg,nlev,levelbox,iparentbox,
     1        ichildbox,icolbox, irowbox, nboxes, nblevel,
     2        iboxlev, istartlev, convective2,iprec,iperiod,
     3        fmmpot2,ttotal2)
  
      call getderivatives8(pot2,
     1     nlev,ichildbox,nblevel,iboxlev,istartlev,
     2     fmmpotx2,fmmpoty2,fmmpotxy2,fmmpotxx2,fmmpotyy2,
     3     fmmlap2,coeffsu6)
c
c
        do i=1, nboxes
        do j=1, npbox
            frightG1(j,i) = fmmpotxx(j,i)+fmmpotxy2(j,i)
            frightG2(j,i) = fmmpotxy(j,i)+fmmpotyy2(j,i)

            frightS1(j,i) = convective1(j,i) - frightG1(j,i)
            frightS2(j,i) = convective2(j,i) - frightG2(j,i)
        enddo
        enddo
c

        do i=1, nboxes
        do j=1, npbox
          fs(1,j,i)=frightS1(j,i)
          fs(2,j,i)=frightS2(j,i)
        enddo
        enddo
c
      do j = 1,nboxes
        do i = 1,64
          fright(i,j) = 0.0d0
          fright(i,j) = 0.0d0
          convective1(i,j)=0.0d0
          convective2(i,j)=0.0d0
          fmmpot(i,j) = 0.0d0
          fmmpot2(i,j) = 0.0d0
        enddo
      enddo

c     compute fs(t_{n+1/2})
c
      call getderivatives8(usol_half1,
     1      nlev, ichildbox, nblevel, iboxlev, istartlev,
     2      frightinitialx1, frightinitialy1, frightinitialxy1,
     3      frightinitialxx1,frightinitialyy1,frightinitiallap1,
     4       coeffsu)
c
      call getderivatives8(usol_half2,
     1      nlev, ichildbox, nblevel, iboxlev, istartlev,
     2      frightinitialx2, frightinitialy2, frightinitialxy2,
     3      frightinitialxx2,frightinitialyy2,frightinitiallap2,
     4       coeffsu2)
c
      call setf8(fright,icolbox, irowbox, ichildbox,
     1          nlev, nblevel, iboxlev, istartlev,fforce1,t_half)
c
      call setf8(fright2,icolbox,irowbox,ichildbox,
     1           nlev,nblevel,iboxlev,istartlev,fforce2,t_half)
c     convective term of ns eq
c
        do i=1, nboxes
        do j=1, npbox
                convective1(j,i) = fright(j,i)-
     1     usol_half1(j,i)*frightinitialx1(j,i)-
     2     usol_half2(j,i)*frightinitialy1(j,i)

                convective2(j,i) = fright2(j,i)-
     1     usol_half1(j,i)*frightinitialx2(j,i)-
     2     usol_half2(j,i)*frightinitialy2(j,i)
        enddo
        enddo
c
      call adaptfmm2d8(ndeg, nlev, levelbox, iparentbox,
     1     ichildbox, icolbox, irowbox, nboxes, nblevel,
     2     iboxlev, istartlev, convective1, iprec, iperiod,
     3     fmmpot, ttotal1)

      call getderivatives8(fmmpot,
     1      nlev, ichildbox, nblevel, iboxlev, istartlev,
     2      fmmpotx, fmmpoty, fmmpotxy,fmmpotxx,
     3      fmmpotyy, fmmlap, coeffsu5)


      call adaptfmm2d8(ndeg,nlev,levelbox,iparentbox,
     1        ichildbox,icolbox, irowbox, nboxes, nblevel,
     2        iboxlev, istartlev, convective2,iprec,iperiod,
     3        fmmpot2,ttotal2)
  
      call getderivatives8(fmmpot2,
     1     nlev,ichildbox,nblevel,iboxlev,istartlev,
     2     fmmpotx2,fmmpoty2,fmmpotxy2,fmmpotxx2,fmmpotyy2,
     3     fmmlap2,coeffsu6)
c
c
        do i=1, nboxes
        do j=1, npbox
            frightG1(j,i) = fmmpotxx(j,i)+fmmpotxy2(j,i)
            frightG2(j,i) = fmmpotxy(j,i)+fmmpotyy2(j,i)

            frightS1(j,i) = convective1(j,i) - frightG1(j,i)
            frightS2(j,i) = convective2(j,i) - frightG2(j,i)
        enddo
        enddo
c

        do i=1, nboxes
        do j=1, npbox
          fs(3,j,i)=frightS1(j,i)
          fs(4,j,i)=frightS2(j,i)
        enddo
        enddo
c
      do j = 1,nboxes
        do i = 1,64
          fright(i,j) = 0.0d0
          fright(i,j) = 0.0d0
          convective1(i,j)=0.0d0
          convective2(i,j)=0.0d0
          fmmpot(i,j) = 0.0d0
          fmmpot2(i,j) = 0.0d0
        enddo
      enddo

c     compute fs(t_{n+1})
c
      call getderivatives8(usol_now1,
     1      nlev, ichildbox, nblevel, iboxlev, istartlev,
     2      frightinitialx1, frightinitialy1, frightinitialxy1,
     3      frightinitialxx1,frightinitialyy1,frightinitiallap1,
     4       coeffsu)
c
      call getderivatives8(usol_now2,
     1      nlev, ichildbox, nblevel, iboxlev, istartlev,
     2      frightinitialx2, frightinitialy2, frightinitialxy2,
     3      frightinitialxx2,frightinitialyy2,frightinitiallap2,
     4       coeffsu2)
c
      call setf8(fright,icolbox, irowbox, ichildbox,
     1          nlev, nblevel, iboxlev, istartlev,fforce1,t_now)
c
      call setf8(fright2,icolbox,irowbox,ichildbox,
     1           nlev,nblevel,iboxlev,istartlev,fforce2,t_now)
c     convective term of ns eq
c
        do i=1, nboxes
        do j=1, npbox
                convective1(j,i) = fright(j,i)-
     1     usol_half1(j,i)*frightinitialx1(j,i)-
     2     usol_half2(j,i)*frightinitialy1(j,i)

                convective2(j,i) = fright2(j,i)-
     1     usol_half1(j,i)*frightinitialx2(j,i)-
     2     usol_half2(j,i)*frightinitialy2(j,i)
        enddo
        enddo
c
      call adaptfmm2d8(ndeg, nlev, levelbox, iparentbox,
     1     ichildbox, icolbox, irowbox, nboxes, nblevel,
     2     iboxlev, istartlev, convective1, iprec, iperiod,
     3     fmmpot, ttotal1)

      call getderivatives8(fmmpot,
     1      nlev, ichildbox, nblevel, iboxlev, istartlev,
     2      fmmpotx, fmmpoty, fmmpotxy,fmmpotxx,
     3      fmmpotyy, fmmlap, coeffsu5)


      call adaptfmm2d8(ndeg,nlev,levelbox,iparentbox,
     1        ichildbox,icolbox, irowbox, nboxes, nblevel,
     2        iboxlev, istartlev, convective2,iprec,iperiod,
     3        fmmpot2,ttotal2)
  
      call getderivatives8(fmmpot2,
     1     nlev,ichildbox,nblevel,iboxlev,istartlev,
     2     fmmpotx2,fmmpoty2,fmmpotxy2,fmmpotxx2,fmmpotyy2,
     3     fmmlap2,coeffsu6)
c
c
        do i=1, nboxes
        do j=1, npbox
            frightG1(j,i) = fmmpotxx(j,i)+fmmpotxy2(j,i)
            frightG2(j,i) = fmmpotxy(j,i)+fmmpotyy2(j,i)

            frightS1(j,i) = convective1(j,i) - frightG1(j,i)
            frightS2(j,i) = convective2(j,i) - frightG2(j,i)
        enddo
        enddo
c

        do i=1, nboxes
        do j=1, npbox
          fs(5,j,i)=frightS1(j,i)
          fs(6,j,i)=frightS2(j,i)
        enddo
        enddo
c
      do i =1,ntarg
         xg(i)=targs(1,i)
         yg(i)=targs(2,i)
      enddo
c
      call interppot(levelbox, icolbox, irowbox,
     1     nboxes, nlev, iparentbox, ichildbox,
     2     nblevel, iboxlev, istartlev, ndeg,frightS1,
     3     ntarg, xg, yg, fsp_now1, cent0, xsize0)
c
      call interppot(levelbox, icolbox, irowbox,
     1     nboxes, nlev, iparentbox, ichildbox,
     2     nblevel, iboxlev, istartlev, ndeg,frightS2,
     3     ntarg, xg, yg, fsp_now2, cent0, xsize0)


      end subroutine


c  compute helomholtz decomposition FS (nordert = 2)
c  we need FS(t_n) FS(t_{n+1/2}) FS(t_{n+1})
c  then need u_{tn} u_{t_{n+1/2}} u_{t_{n+1}}
c  the latter two term using BDF method
      subroutine compute_fs2(ndeg,nordert,nlev,levelbox,iparentbox,
     1        ichildbox,icolbox, irowbox, nboxes, nblevel,
     2        iboxlev, istartlev,iprec,iperiod, nnodes,
     3        fforce1,fforce2,dt,tt,nnoeds,usol_pre,
     4        fs,fsp_now1,fsp_now2,ntarg,targs,norder,ndim)
      implicit real*8 (a-h,o-z)
c
      integer nlev,ndeg,nordert,ntarg,norder,ndim
      integer levelbox(nboxes), iparentbox(nboxes)
      integer ichildbox(4,nboxes), icolleagbox(9,nboxes)
      integer icolbox(nboxes), irowbox(nboxes)
      integer nblevel(0:nlev), iboxlev(nboxes)
      integer itemparray(nboxes)
      integer istartlev(0:nlev)
      integer iflag(nboxes)
      integer noodes
      real *8 targs(nd,ntarg)

      real *8 fright(64,nboxes),fright2(64,nboxes)
      real *8 usol_pre(nordert*2,ndeg**2,nboxes)
c
      real *8 usol_prex1(ndeg**2,nboxes)
c
      real *8 usol_prey1(ndeg**2,nboxes)
c
      real *8 usol_prex2(ndeg**2,nboxes)
c
      real *8 usol_prey2(ndeg**2,nboxes)
c

      real *8 usol_now1(ndeg**2,nboxes)

      real *8 usol_now2(ndeg**2,nboxes)

      real*8 convective1(64,nboxes),convective2(64,nboxes)
      

      real *8 cent0(2), xsize0
c     output
      real *8 fs(nnodes*2,ndeg**2,nboxes)
c
c     initial value's derivates
      real *8  frightinitialx1(64,nboxes)
      real *8  frightinitialy1(64,nboxes)
      real *8  frightinitialx2(64,nboxes)
      real *8  frightinitialy2(64,nboxes)
      real *8  frightinitialxy1(64,nboxes)
      real *8  frightinitialxx1(64,nboxes)
      real *8  frightinitialyy1(64,nboxes)
      real *8  frightinitiallap1(64,nboxes)
      real *8  frightinitialxy2(64,nboxes)
      real *8  frightinitialxx2(64,nboxes)
      real *8  frightinitialyy2(64,nboxes)
      real *8  frightinitiallap2(64,nboxes)
c
      real *8  coeffsu(0:7,0:7,nboxes)
      real *8  coeffsu2(0:7,0:7,nboxes)
      real *8  coeffsu3(0:7,0:7,nboxes)
      real *8  coeffsu4(0:7,0:7,nboxes)
      real *8  coeffsu5(0:7,0:7,nboxes)
      real *8  coeffsu6(0:7,0:7,nboxes)
      real *8  coeffsu7(0:7,0:7,nboxes)
      real *8  coeffsu8(0:7,0:7,nboxes)
c
      
      real *8  fmmpot(64,nboxes)
      real *8  fmmpotx(64,nboxes), fmmpoty(64,nboxes)
      real *8  fmmpotxy(64,nboxes)
      real *8  fmmpotxx(64,nboxes),fmmpotyy(64,nboxes),fmmlap(64,nboxes)
      real *8  fmmpot2(64,nboxes)
      real *8  fmmpotx2(64,nboxes),fmmpoty2(64,nboxes)
      real *8  fmmpotxy2(64,nboxes)
      real *8  fmmpotxx2(64,nboxes),fmmpotyy2(64,nboxes)
      real *8  fmmlap2(64,nboxes)

c    helmholtz
      real *8  frightG1(norder**ndim,nboxes)
      real *8  frightG2(norder**ndim,nboxes)
      real *8  frightS1(norder**ndim,nboxes)
      real *8  frightS2(norder**ndim,nboxes)
c
      real *8  xg(ntarg)
      real *8  yg(ntarg)
c
      real *8  fsp_now1(ntarg)
      real *8  fsp_now2(ntarg)


      external fforce1,fforce2
c
c
      if (nordert .eq. 2)then
        do j = 1,nboxes
           do i = 1,64
           usol_prex1(i,j) = usol_pre(1,i,j)
           usol_prey1(i,j) = usol_pre(2,i,j)
           usol_prex2(i,j) = usol_pre(3,i,j)
           usol_prey2(i,j) = usol_pre(4,i,j)

c  compute u_{n+1} nordert=2
           usol_now1(i,j) = 2.0d0*usol_prex1(i,j)-
     1    usol_prex2(i,j)+4.0d0

           usol_now2(i,j) = 2.0d0*usol_prey1(i,j)-
     1    usol_prey2(i,j)
    
           enddo
        enddo
      endif

      t_pre=tt-dt
      t_now = tt
c
      do j = 1,nboxes
        do i = 1,64
          fright(i,j) = 0.0d0
          fright2(i,j) = 0.0d0
          convective1(i,j)=0.0d0
          convective2(i,j)=0.0d0
          fmmpot(i,j) = 0.0d0
          fmmpot2(i,j) = 0.0d0
        enddo
      enddo

c     compute fs(t_n)
c
      call getderivatives8(usol_prex1,
     1      nlev, ichildbox, nblevel, iboxlev, istartlev,
     2      frightinitialx1, frightinitialy1, frightinitialxy1,
     3      frightinitialxx1,frightinitialyy1,frightinitiallap1,
     4       coeffsu)
c
      call getderivatives8(usol_prey1,
     1      nlev, ichildbox, nblevel, iboxlev, istartlev,
     2      frightinitialx2, frightinitialy2, frightinitialxy2,
     3      frightinitialxx2,frightinitialyy2,frightinitiallap2,
     4       coeffsu2)
c
      call setf8(fright,icolbox, irowbox, ichildbox,
     1          nlev, nblevel, iboxlev, istartlev,fforce1,t_pre)
c
      call setf8(fright2,icolbox,irowbox,ichildbox,
     1           nlev,nblevel,iboxlev,istartlev,fforce2,t_pre)
c     convective term of ns eq
c
        do i=1, nboxes
        do j=1, npbox
                convective1(j,i) = fright(j,i)-
     1     usol_prex1(j,i)*frightinitialx1(j,i)-
     2     usol_prey1(j,i)*frightinitialy1(j,i)

                convective2(j,i) = fright2(j,i)-
     1     usol_prex1(j,i)*frightinitialx2(j,i)-
     2     usol_prey1(j,i)*frightinitialy2(j,i)
        enddo
        enddo
c
      call adaptfmm2d8(ndeg, nlev, levelbox, iparentbox,
     1     ichildbox, icolbox, irowbox, nboxes, nblevel,
     2     iboxlev, istartlev, convective1, iprec, iperiod,
     3     fmmpot, ttotal1)

      call getderivatives8(pot,
     1      nlev, ichildbox, nblevel, iboxlev, istartlev,
     2      fmmpotx, fmmpoty, fmmpotxy,fmmpotxx,
     3      fmmpotyy, fmmlap, coeffsu5)


      call adaptfmm2d8(ndeg,nlev,levelbox,iparentbox,
     1        ichildbox,icolbox, irowbox, nboxes, nblevel,
     2        iboxlev, istartlev, convective2,iprec,iperiod,
     3        fmmpot2,ttotal2)
  
      call getderivatives8(pot2,
     1     nlev,ichildbox,nblevel,iboxlev,istartlev,
     2     fmmpotx2,fmmpoty2,fmmpotxy2,fmmpotxx2,fmmpotyy2,
     3     fmmlap2,coeffsu6)
c
c
        do i=1, nboxes
        do j=1, npbox
            frightG1(j,i) = fmmpotxx(j,i)+fmmpotxy2(j,i)
            frightG2(j,i) = fmmpotxy(j,i)+fmmpotyy2(j,i)

            frightS1(j,i) = convective1(j,i) - frightG1(j,i)
            frightS2(j,i) = convective2(j,i) - frightG2(j,i)
        enddo
        enddo
c

        do i=1, nboxes
        do j=1, npbox
          fs(1,j,i)=frightS1(j,i)
          fs(2,j,i)=frightS2(j,i)
        enddo
        enddo

c
      do j = 1,nboxes
        do i = 1,64
          fright(i,j) = 0.0d0
          fright(i,j) = 0.0d0
          convective1(i,j)=0.0d0
          convective2(i,j)=0.0d0
          fmmpot(i,j) = 0.0d0
          fmmpot2(i,j) = 0.0d0
        enddo
      enddo

c     compute fs(t_{n+1})
c
      call getderivatives8(usol_now1,
     1      nlev, ichildbox, nblevel, iboxlev, istartlev,
     2      frightinitialx1, frightinitialy1, frightinitialxy1,
     3      frightinitialxx1,frightinitialyy1,frightinitiallap1,
     4       coeffsu)
c
      call getderivatives8(usol_now2,
     1      nlev, ichildbox, nblevel, iboxlev, istartlev,
     2      frightinitialx2, frightinitialy2, frightinitialxy2,
     3      frightinitialxx2,frightinitialyy2,frightinitiallap2,
     4       coeffsu2)
c
      call setf8(fright,icolbox, irowbox, ichildbox,
     1          nlev, nblevel, iboxlev, istartlev,fforce1,t_now)
c
      call setf8(fright2,icolbox,irowbox,ichildbox,
     1           nlev,nblevel,iboxlev,istartlev,fforce2,t_now)
c     convective term of ns eq
c
        do i=1, nboxes
        do j=1, npbox
                convective1(j,i) = fright(j,i)-
     1     usol_half1(j,i)*frightinitialx1(j,i)-
     2     usol_half2(j,i)*frightinitialy1(j,i)

                convective2(j,i) = fright2(j,i)-
     1     usol_half1(j,i)*frightinitialx2(j,i)-
     2     usol_half2(j,i)*frightinitialy2(j,i)
        enddo
        enddo
c
      call adaptfmm2d8(ndeg, nlev, levelbox, iparentbox,
     1     ichildbox, icolbox, irowbox, nboxes, nblevel,
     2     iboxlev, istartlev, convective1, iprec, iperiod,
     3     fmmpot, ttotal1)

      call getderivatives8(fmmpot,
     1      nlev, ichildbox, nblevel, iboxlev, istartlev,
     2      fmmpotx, fmmpoty, fmmpotxy,fmmpotxx,
     3      fmmpotyy, fmmlap, coeffsu5)


      call adaptfmm2d8(ndeg,nlev,levelbox,iparentbox,
     1        ichildbox,icolbox, irowbox, nboxes, nblevel,
     2        iboxlev, istartlev, convective2,iprec,iperiod,
     3        fmmpot2,ttotal2)
  
      call getderivatives8(fmmpot2,
     1     nlev,ichildbox,nblevel,iboxlev,istartlev,
     2     fmmpotx2,fmmpoty2,fmmpotxy2,fmmpotxx2,fmmpotyy2,
     3     fmmlap2,coeffsu6)
c
c
        do i=1, nboxes
        do j=1, npbox
            frightG1(j,i) = fmmpotxx(j,i)+fmmpotxy2(j,i)
            frightG2(j,i) = fmmpotxy(j,i)+fmmpotyy2(j,i)

            frightS1(j,i) = convective1(j,i) - frightG1(j,i)
            frightS2(j,i) = convective2(j,i) - frightG2(j,i)
        enddo
        enddo
c

        do i=1, nboxes
        do j=1, npbox
          fs(3,j,i)=frightS1(j,i)
          fs(4,j,i)=frightS2(j,i)
        enddo
        enddo

c
c
      do i =1,ntarg
         xg(i)=targs(1,i)
         yg(i)=targs(2,i)
      enddo

      call interppot(levelbox, icolbox, irowbox,
     1     nboxes, nlev, iparentbox, ichildbox,
     2     nblevel, iboxlev, istartlev, ndeg,frightS1,
     3     ntarg, xg, yg, fsp_now1, cent0, xsize0)
c
      call interppot(levelbox, icolbox, irowbox,
     1     nboxes, nlev, iparentbox, ichildbox,
     2     nblevel, iboxlev, istartlev, ndeg,frightS2,
     3     ntarg, xg, yg, fsp_now2, cent0, xsize0)

      end subroutine

c**************************************************************
c     subroutine initpot2
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
c  nd=2!
c
c**************************************************************
c
      subroutine initpot2(nd,norder, ltree, itree, iptr, centers,
     1           nlevels, boxsize, nboxes, fright, dt,
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
      real*8 pot(nd,norder**ndim,nboxes), potp(nd,ntarg)
c     local vars
      real*8, allocatable:: timeinfo(:)
      real*8, allocatable:: grad(:,:,:,:)
      real*8, allocatable:: hess(:,:,:,:)
      real*8, allocatable:: gradp(:,:,:)
      real*8, allocatable:: hessp(:,:,:)
c
      pi=3.1415926535897932d0
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
          do l =1,nd
          pot(l,i,j)=0.0d0
          enddo
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

c**************************************************************
c     subroutine volpot2
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
      subroutine volpot2(nd,norder, nordert, ltree, itree, iptr,
     1           centers, nlevels, boxsize, nboxes,tpre,
     2           dt, iperiod, ntarg, targs, ifvtarg, ifptarg,
     3           eps, vpot, vpotp,fs,fsp_now1,fsp_now2,nnodes)
      implicit real*8 (a-h,o-z)
c     global vars
      integer norder, ltree, nboxes, ndim, nlevels
      integer iperiod, ntarg,nd,nnodes
      integer ifvtarg, ifptarg
      parameter(ndim=2)
      parameter(npbox=64)
      integer itree(ltree), iptr(12)
      real*8 centers(ndim, nboxes), boxsize(0:nlevels)
      real*8 targs(ndim,ntarg), eps, dt
      real*8 targ(ndim)
      real*8 fs(nnodes*2,64,nboxes)
      real*8 vpot(nd,norder**ndim,nboxes), vpotp(nd,ntarg)
      real*8 fright(nd,norder**ndim,nboxes)
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
      real*8, allocatable:: vpotinc(:,:,:), vpotpinc(:,:)
      real*8, allocatable:: xref(:,:), wts(:)
      external fforce
c
      pi=3.1415926535897932d0
      
    
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
      allocate(vpotinc(nd,npbox,nboxes))
      allocate(vpotpinc(nd,ntarg))
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
      do jj=1,nnodes-1
        tt=tpre+tnodes(jj)
        delta=dt-tnodes(jj)
        scal=dt
c       fright assigned
c       assign fright: fforce(targ,tt,val) val->fright
        do i =1,nboxes
           do j = 1,64
              fright(1,j,i)=0.0d0
              fright(2,j,i)=0.0d0
           enddo
        enddo

        do ilevel=0,nlevels
          bs = boxsize(ilevel)/2.0d0
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
c              a leaf node
c              get the grid and compute the forcing term
               do j=1,npbox
                 fright(1,j,ibox)=fs(jj*2-1,j,ibox)
                 fright(2,j,ibox)=fs(jj*2,j,ibox)
               enddo
            else
c             otherwise assign fright to be zero
              do j=1, npbox
                fright(1,j,ibox)=0.0d0
                fright(2,j,ibox)=0.0d0
              enddo
            endif
          enddo
        enddo
c
        do i=1,nboxes
          do j=1, npbox
            vpotinc(1,j,i)=0.0d0
            vpotinc(2,j,i)=0.0d0
          enddo
        enddo

c
        do i=1,nboxes
          do j=1, npbox
            fright(1,j,i)=0.0d0
            fright(2,j,i)=0.0d0
          enddo
        enddo

        do i = 1,ntarg
          vpotpinc(1,i)=0.0d0
          vpotpinc(2,i)=0.0d0
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
               vpot(1,j,i)=vpot(1,j,i)+vpotinc(1,j,i)/(4.0d0*pi)
     1                  /(dt-tnodes(jj))*dt*twhts(jj)
               vpot(2,j,i)=vpot(2,j,i)+vpotinc(2,j,i)/(4.0d0*pi)
     1                  /(dt-tnodes(jj))*dt*twhts(jj)

            enddo
          endif
        enddo
       
c      for target point
c
      if(ifptarg .gt. 0) then
        do i=1, ntarg
              vpotp(1,i)=vpotp(1,i)+vpotpinc(1,i)/(4.0d0*pi)
     1                  /(dt-tnodes(jj))*dt*twhts(jj)
              vpotp(2,i)=vpotp(2,i)+vpotpinc(2,i)/(4.0d0*pi)
     1                  /(dt-tnodes(jj))*dt*twhts(jj)
        enddo
      endif
 
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
c
               vpot(1,j,ibox)=vpot(1,j,ibox)+
     1        fs(nnodes*2-1,j,ibox)*dt*twhts(jj)
c
               vpot(2,j,ibox)=vpot(2,j,ibox)+
     1        fs(nnodes*2,j,ibox)*dt*twhts(jj)

             enddo
          endif
        enddo
      enddo
c     contribution from t(nnodes): function value at
c     the current time step, no transform needed
c     targets
      jj=nnodes
      tt=tpre+tnodes(jj)
      do i = 1,ntarg
        vpotp(1,i)=vpotp(1,i)+fsp_now1(i)*dt*twhts(jj)
        vpotp(2,i)=vpotp(2,i)+fsp_now2(i)*dt*twhts(jj)
      enddo
c
      end subroutine
