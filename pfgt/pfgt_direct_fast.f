c     This file contains nine direction evaluation subroutines
c     for point Gauss transform in 1,2,3 dimensions.
c
c     Naming convention for subroutines:
c     1. all subroutines start with gnd_direct, indicating that 
c        they are for direct evaluation of Gauss transform in general dimensions;
c     2. they are followed by c (charge input only), d (dipole input only),
c        cd (charge+dipole input), p (potential output only), g (potential+gradient output),
c        h (potential+gradient+hessian output).
c
c     pot(ii,i)  = \sum_j charge(ii,j)*exp(-(t_i-s_j)^2/delta)
c                + \sum_j dipstr(ii,j)*(grad_s(exp(-(t_i-s_j)^2/delta) \cdot rnormal(s_j))
c
c     grad  = d(pot)
c     hess  = dd(pot)
c     Note: in 2d, hessian has the order dxx, dxy, dyy
c           in 3d, hessian has the order dxx, dyy, dzz, dxy, dxz, dyz
c     
c     Note: all output variables are incremented. So proper initialization are needed, but
c     the subroutines can be called many times to calculate the output variable due to all
c     sources.
c
c     The union of input and output arguments are as follows.
c
c     Input parameters:
c     nd: number of input vectors (charge, rnormal, dipstr) and output vectors (pot, grad, hess)
c     dim: dimension of the underlying space
c     delta: the Gaussian variance
c     dmax: the squared distance outside which the Gaussian kernel is regarded as 0
c     ns: number of sources
c     sources: (dim,ns) source coordinates
c     ntarg: number of targets
c     targ: (dim,ntarg) target coordinates
c     charge: (nd,ns) charge strengths
c     rnormal: (dim,ns) dipole orientation vectors
c     dipstr: (nd,ns) dipole strengths
c
c     Output parameters:
c     pot: (nd,ntarg) incremented potential at targets
c     grad: (nd,dim,ntarg) incremented gradient at targets
c     hess: (nd,dim*(dim+1)/2,ntarg) incremented hessian at targets
c***********************************************************************
c
c     charge to potential
c
c**********************************************************************
      subroutine g3d_directcp(nd,dim,eps,delta,dmax,sources,ns,charge,
     $           xtarg,ytarg,ztarg,ntarg,pot)
      implicit none
c**********************************************************************
      integer dim,ntarg
      integer ns,nd,digits
      real *8 sources(dim,ns),xtarg(*),ytarg(*),ztarg(*)
      real *8 dr(dim)
      real *8 rtmp,delta,dmax,deltainvm,eps
      real *8 pot(nd,ntarg)
      real *8 charge(nd,ns)
c
      deltainvm=-1.0d0/delta
      digits = log10(1.0d0/eps)
      call g3d_directcp_cpp(nd,dim,digits,deltainvm,
     $    sources,ns,charge,xtarg,ytarg,ztarg,ntarg,pot)
      
      return
      end
c
c
c
