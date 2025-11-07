c     a simple test for the subroutine secant
c     which solves a nonlinear equation
c     of the form:
c     u = g +alpha * f(u,x,t)
c
      program test
      implicit real*8 (a-h,o-z)
      real*8 t, x, y, uext, usol
      real*8 alpha, f, g
      external rfun
c
      x=1.0d0
      y=2.0d0
      t = 0.5d0
      alpha = 0.1d0
c
      uext = 3.0d0*y**2+x
c     uext = 13.0
      f = y*sin(10.0d0*t*uext)+t
      g = uext - alpha*f
c     eqn set up
c      write(*,*) x, y, t
c      write(*,*) alpha, g 
c      write(*,*) uext
c      write(*,*) ''
c
c     now set parameters and call secant to solve
      msi = 10
      ndamp = 1
c     since uext =5.0d0
      u0 = 11.9d0
      u1 = 11.8d0
c
      tol =1.0d-12
      call secant(rfun,t,x,y,
     1     u0,u1,g,alpha,usol,fusol,tol,msi,isi,
     2     ndamp,ier)
c
      write(*,*) 'ier=', ier
      write(*,*) 'isi=', isi
      write(*,*) 'usol=', usol
      write(*,*) 'error=', abs(usol-uext)

c     test example:
c     solution:
c     u(t,x)= 3*x(2)^2+x(1)
c
c     equation:
c     f(t,x,u) = x(2)*sin(10*t*u)+t
c     so that it has some nontrivial adaptivity in t
c     fix t <= T=1, so that it's not crazy
c
c     alpha = 0.1
c     u = g+alpha*f(t,x,u)
c     i.e. g = u-alpha*f(t,x,u)
c
c     f given by subroutine rfun
c     with calling seq
c     rfun(t,x,u,f)

      end program



      subroutine rfun(t,x,y,u,f)
      implicit real*8 (a-h,o-z)
      real*8 t, x, y, u, f
c
      f = y*sin(10.0d0*t*u)+t
c
      end subroutine
