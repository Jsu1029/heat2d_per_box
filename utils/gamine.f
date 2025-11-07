c*****************************************************************
c     subroutines pulled from netlib to compute incomplete
c     gamma function: dgamic and dependencies
c
c     with one function added: 
c          ei32: exponential integral of order 3/2
c                which is based on incomplete gamma function
c
c*****************************************************************
c
      function ei32(x)
c     a function that computes one specific special function
c     that I need: exponential integral of order 3/2
c     namely: ei32(x)=int_1^{inf} exp(-t*x)/t^{3/2} dt
c     input: x>=0
c     output: the function value
c-------------------------------------
c
      implicit real*8 (a-h,o-z)
      real*8 x, ei32

      ei32=2.0d0*(dexp(-x)-dsqrt(x)*dgamic(0.5d0,x))

      end function





c------------------------------------------------------
c
        subroutine expint1(x,e1)
c
c       ============================================
c       purpose: compute exponential integral expint1(x)
c       input :  x  --- argument of expint1(x)
c       output:  e1 --- expint1(x)  ( x > 0 )
c       ============================================
c
        implicit double precision (a-h,o-z)
        if (x.eq.0.0) then
           e1=1.0d+300
        else if (x.le.1.0) then
           e1=1.0d0
           r=1.0d0
           do 10 k=1,25
              r=-r*k*x/(k+1.0d0)**2
              e1=e1+r
              if (dabs(r).le.dabs(e1)*1.0d-15) go to 15
10         continue
15         ga=0.5772156649015328d0
           e1=-ga-dlog(x)+x*e1
        else
           m=20+int(80.0/x)
           t0=0.0d0
           do 20 k=m,1,-1
              t0=k/(1.0d0+k/(x+t0))
20         continue
           t=1.0d0/(x+t0)
           e1=dexp(-x)*t
        endif

        end






c----------------------------------------------------------------
*deck dgamic
      double precision function dgamic (a, x)
c***begin prologue  dgamic
c***purpose  calculate the complementary incomplete gamma function.
c***library   slatec (fnlib)
c***category  c7e
c***type      double precision (gamic-s, dgamic-d)
c***keywords  complementary incomplete gamma function, fnlib,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c   evaluate the complementary incomplete gamma function
c
c   dgamic = integral from x to infinity of exp(-t) * t**(a-1.)  .
c
c   dgamic is evaluated for arbitrary real values of a and for non-
c   negative values of x (even though dgamic is defined for x .lt.
c   0.0), except that for x = 0 and a .le. 0.0, dgamic is undefined.
c
c   dgamic, a, and x are double precision.
c
c   a slight deterioration of 2 or 3 digits accuracy will occur when
c   dgamic is very large or very small in absolute value, because log-
c   arithmic variables are used.  also, if the parameter a is very close
c   to a negative integer (but not a negative integer), there is a loss
c   of accuracy, which is reported if the result is less than half
c   machine precision.
c
c***references  w. gautschi, a computational procedure for incomplete
c                 gamma functions, acm transactions on mathematical
c                 software 5, 4 (december 1979), pp. 466-481.
c               w. gautschi, incomplete gamma functions, algorithm 542,
c                 acm transactions on mathematical software 5, 4
c                 (december 1979), pp. 482-489.
c***routines called  d1mach, d9gmic, d9gmit, d9lgic, d9lgit, dlgams,
c                    dlngam, xerclr, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920528  description and references sections revised.  (wrb)
c***end prologue  dgamic
      double precision a, x, aeps, ainta, algap1, alneps, alngs, alx,
     1  bot, e, eps, gstar, h, sga, sgng, sgngam, sgngs, sqeps, t,
     2  d1mach, dlngam, d9gmic, d9gmit, d9lgic, d9lgit
      logical first
      save eps, sqeps, alneps, bot, first
      data first /.true./
c***first executable statement  dgamic
      if (first) then
         eps = 0.5d0*d1mach(3)
         sqeps = sqrt(d1mach(4))
         alneps = -log (d1mach(3))
         bot = log (d1mach(1))
      endif
      first = .false.
c
      if (x .lt. 0.d0) call xermsg ('slatec', 'dgamic', 'x is negative'
     +   , 2, 2)
c
      if (x.gt.0.d0) go to 20
      if (a .le. 0.d0) call xermsg ('slatec', 'dgamic',
     +   'x = 0 and a le 0 so dgamic is undefined', 3, 2)
c
      dgamic = exp (dlngam(a+1.d0) - log(a))
      return
c
 20   alx = log (x)
      sga = 1.0d0
      if (a.ne.0.d0) sga = sign (1.0d0, a)
      ainta = aint (a + 0.5d0*sga)
      aeps = a - ainta
c
      izero = 0
      if (x.ge.1.0d0) go to 40
c
      if (a.gt.0.5d0 .or. abs(aeps).gt.0.001d0) go to 30
      e = 2.0d0
      if (-ainta.gt.1.d0) e = 2.d0*(-ainta+2.d0)/(ainta*ainta-1.0d0)
      e = e - alx * x**(-0.001d0)
      if (e*abs(aeps).gt.eps) go to 30
c
      dgamic = d9gmic (a, x, alx)
      return
c
 30   call dlgams (a+1.0d0, algap1, sgngam)
      gstar = d9gmit (a, x, algap1, sgngam, alx)
      if (gstar.eq.0.d0) izero = 1
      if (gstar.ne.0.d0) alngs = log (abs(gstar))
      if (gstar.ne.0.d0) sgngs = sign (1.0d0, gstar)
      go to 50
c
 40   if (a.lt.x) dgamic = exp (d9lgic(a, x, alx))
      if (a.lt.x) return
c
      sgngam = 1.0d0
      algap1 = dlngam (a+1.0d0)
      sgngs = 1.0d0
      alngs = d9lgit (a, x, algap1)
c
c evaluation of dgamic(a,x) in terms of tricomi-s incomplete gamma fn.
c
 50   h = 1.d0
      if (izero.eq.1) go to 60
c
      t = a*alx + alngs
      if (t.gt.alneps) go to 70
      if (t.gt.(-alneps)) h = 1.0d0 - sgngs*exp(t)
c
      if (abs(h).lt.sqeps) call xerclr
      if (abs(h) .lt. sqeps) call xermsg ('slatec', 'dgamic',
     +   'result lt half precision', 1, 1)
c
 60   sgng = sign (1.0d0, h) * sga * sgngam
      t = log(abs(h)) + algap1 - log(abs(a))
      if (t.lt.bot) call xerclr
      dgamic = sgng * exp(t)
      return
c
 70   sgng = -sgngs * sga * sgngam
      t = t + algap1 - log(abs(a))
      if (t.lt.bot) call xerclr
      dgamic = sgng * exp(t)
      return
c
      end






c------------------------------------------------------------------
*deck d9gmic
      double precision function d9gmic (a, x, alx)
c***begin prologue  d9gmic
c***subsidiary
c***purpose  compute the complementary incomplete gamma function for a
c            near a negative integer and x small.
c***library   slatec (fnlib)
c***category  c7e
c***type      double precision (r9gmic-s, d9gmic-d)
c***keywords  complementary incomplete gamma function, fnlib, small x,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c compute the complementary incomplete gamma function for a near
c a negative integer and for small x.
c
c***references  (none)
c***routines called  d1mach, dlngam, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   890911  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900720  routine changed from user-callable to subsidiary.  (wrb)
c***end prologue  d9gmic
      double precision a, x, alx, alng, bot, eps, euler, fk, fkp1, fm,
     1  s, sgng, t, te, d1mach, dlngam
      logical first
      save euler, eps, bot, first
      data euler / 0.5772156649 0153286060 6512090082 40 d0 /
      data first /.true./
c***first executable statement  d9gmic
      if (first) then
         eps = 0.5d0*d1mach(3)
         bot = log (d1mach(1))
      endif
      first = .false.
c
      if (a .gt. 0.d0) call xermsg ('slatec', 'd9gmic',
     +   'a must be near a negative integer', 2, 2)
      if (x .le. 0.d0) call xermsg ('slatec', 'd9gmic',
     +   'x must be gt zero', 3, 2)
c
      m = -(a - 0.5d0)
      fm = m
c
      te = 1.0d0
      t = 1.0d0
      s = t
      do 20 k=1,200
        fkp1 = k + 1
        te = -x*te/(fm+fkp1)
        t = te/fkp1
        s = s + t
        if (abs(t).lt.eps*s) go to 30
 20   continue
      call xermsg ('slatec', 'd9gmic',
     +   'no convergence in 200 terms of continued fraction', 4, 2)
c
 30   d9gmic = -alx - euler + x*s/(fm+1.0d0)
      if (m.eq.0) return
c
      if (m.eq.1) d9gmic = -d9gmic - 1.d0 + 1.d0/x
      if (m.eq.1) return
c
      te = fm
      t = 1.d0
      s = t
      mm1 = m - 1
      do 40 k=1,mm1
        fk = k
        te = -x*te/fk
        t = te/(fm-fk)
        s = s + t
        if (abs(t).lt.eps*abs(s)) go to 50
 40   continue
c
 50   do 60 k=1,m
        d9gmic = d9gmic + 1.0d0/k
 60   continue
c
      sgng = 1.0d0
      if (mod(m,2).eq.1) sgng = -1.0d0
      alng = log(d9gmic) - dlngam(fm+1.d0)
c
      d9gmic = 0.d0
      if (alng.gt.bot) d9gmic = sgng * exp(alng)
      if (s.ne.0.d0) d9gmic = d9gmic +
     1  sign (exp(-fm*alx+log(abs(s)/fm)), s)
c
      if (d9gmic .eq. 0.d0 .and. s .eq. 0.d0) call xermsg ('slatec',
     +   'd9gmic', 'result underflows', 1, 1)
      return
c
      end







c------------------------------------------------------------------
*deck d9gmit
      double precision function d9gmit (a, x, algap1, sgngam, alx)
c***begin prologue  d9gmit
c***subsidiary
c***purpose  compute tricomi's incomplete gamma function for small
c            arguments.
c***library   slatec (fnlib)
c***category  c7e
c***type      double precision (r9gmit-s, d9gmit-d)
c***keywords  complementary incomplete gamma function, fnlib, small x,
c             special functions, tricomi
c***author  fullerton, w., (lanl)
c***description
c
c compute tricomi's incomplete gamma function for small x.
c
c***references  (none)
c***routines called  d1mach, dlngam, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   890911  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900720  routine changed from user-callable to subsidiary.  (wrb)
c***end prologue  d9gmit
      double precision a, x, algap1, sgngam, alx, ae, aeps, algs, alg2,
     1  bot, eps, fk, s, sgng2, t, te, d1mach, dlngam
      logical first
      save eps, bot, first
      data first /.true./
c***first executable statement  d9gmit
      if (first) then
         eps = 0.5d0*d1mach(3)
         bot = log (d1mach(1))
      endif
      first = .false.
c
      if (x .le. 0.d0) call xermsg ('slatec', 'd9gmit',
     +   'x should be gt 0', 1, 2)
c
      ma = a + 0.5d0
      if (a.lt.0.d0) ma = a - 0.5d0
      aeps = a - ma
c
      ae = a
      if (a.lt.(-0.5d0)) ae = aeps
c
      t = 1.d0
      te = ae
      s = t
      do 20 k=1,200
        fk = k
        te = -x*te/fk
        t = te/(ae+fk)
        s = s + t
        if (abs(t).lt.eps*abs(s)) go to 30
 20   continue
      call xermsg ('slatec', 'd9gmit',
     +   'no convergence in 200 terms of taylor-s series', 2, 2)
c
 30   if (a.ge.(-0.5d0)) algs = -algap1 + log(s)
      if (a.ge.(-0.5d0)) go to 60
c
      algs = -dlngam(1.d0+aeps) + log(s)
      s = 1.0d0
      m = -ma - 1
      if (m.eq.0) go to 50
      t = 1.0d0
      do 40 k=1,m
        t = x*t/(aeps-(m+1-k))
        s = s + t
        if (abs(t).lt.eps*abs(s)) go to 50
 40   continue
c
 50   d9gmit = 0.0d0
      algs = -ma*log(x) + algs
      if (s.eq.0.d0 .or. aeps.eq.0.d0) go to 60
c
      sgng2 = sgngam * sign (1.0d0, s)
      alg2 = -x - algap1 + log(abs(s))
c
      if (alg2.gt.bot) d9gmit = sgng2 * exp(alg2)
      if (algs.gt.bot) d9gmit = d9gmit + exp(algs)
      return
c
 60   d9gmit = exp (algs)
      return
c
      end






c------------------------------------------------------------------
*deck d9lgic
      double precision function d9lgic (a, x, alx)
c***begin prologue  d9lgic
c***subsidiary
c***purpose  compute the log complementary incomplete gamma function
c            for large x and for a .le. x.
c***library   slatec (fnlib)
c***category  c7e
c***type      double precision (r9lgic-s, d9lgic-d)
c***keywords  complementary incomplete gamma function, fnlib, large x,
c             logarithm, special functions
c***author  fullerton, w., (lanl)
c***description
c
c compute the log complementary incomplete gamma function for large x
c and for a .le. x.
c
c***references  (none)
c***routines called  d1mach, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900720  routine changed from user-callable to subsidiary.  (wrb)
c***end prologue  d9lgic
      double precision a, x, alx, eps, fk, p, r, s, t, xma, xpa, d1mach
      save eps
      data eps / 0.d0 /
c***first executable statement  d9lgic
      if (eps.eq.0.d0) eps = 0.5d0*d1mach(3)
c
      xpa = x + 1.0d0 - a
      xma = x - 1.d0 - a
c
      r = 0.d0
      p = 1.d0
      s = p
      do 10 k=1,300
        fk = k
        t = fk*(a-fk)*(1.d0+r)
        r = -t/((xma+2.d0*fk)*(xpa+2.d0*fk)+t)
        p = r*p
        s = s + p
        if (abs(p).lt.eps*s) go to 20
 10   continue
      call xermsg ('slatec', 'd9lgic',
     +   'no convergence in 300 terms of continued fraction', 1, 2)
c
 20   d9lgic = a*alx - x + log(s/xpa)
c
      return
      end







c------------------------------------------------------------------
*deck d9lgit
      double precision function d9lgit (a, x, algap1)
c***begin prologue  d9lgit
c***subsidiary
c***purpose  compute the logarithm of tricomi's incomplete gamma
c            function with perron's continued fraction for large x and
c            a .ge. x.
c***library   slatec (fnlib)
c***category  c7e
c***type      double precision (r9lgit-s, d9lgit-d)
c***keywords  fnlib, incomplete gamma function, logarithm,
c             perron's continued fraction, special functions, tricomi
c***author  fullerton, w., (lanl)
c***description
c
c compute the log of tricomi's incomplete gamma function with perron's
c continued fraction for large x and for a .ge. x.
c
c***references  (none)
c***routines called  d1mach, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900720  routine changed from user-callable to subsidiary.  (wrb)
c***end prologue  d9lgit
      double precision a, x, algap1, ax, a1x, eps, fk, hstar, p, r, s,
     1  sqeps, t, d1mach
      logical first
      save eps, sqeps, first
      data first /.true./
c***first executable statement  d9lgit
      if (first) then
         eps = 0.5d0*d1mach(3)
         sqeps = sqrt(d1mach(4))
      endif
      first = .false.
c
      if (x .le. 0.d0 .or. a .lt. x) call xermsg ('slatec', 'd9lgit',
     +   'x should be gt 0.0 and le a', 2, 2)
c
      ax = a + x
      a1x = ax + 1.0d0
      r = 0.d0
      p = 1.d0
      s = p
      do 20 k=1,200
        fk = k
        t = (a+fk)*x*(1.d0+r)
        r = t/((ax+fk)*(a1x+fk)-t)
        p = r*p
        s = s + p
        if (abs(p).lt.eps*s) go to 30
 20   continue
      call xermsg ('slatec', 'd9lgit',
     +   'no convergence in 200 terms of continued fraction', 3, 2)
c
 30   hstar = 1.0d0 - x*s/a1x
      if (hstar .lt. sqeps) call xermsg ('slatec', 'd9lgit',
     +   'result less than half precision', 1, 1)
c
      d9lgit = -x - algap1 - log(hstar)
      return
c
      end







c------------------------------------------------------------------
*deck d9lgmc
      double precision function d9lgmc (x)
c***begin prologue  d9lgmc
c***subsidiary
c***purpose  compute the log gamma correction factor so that
c            log(dgamma(x)) = log(sqrt(2*pi)) + (x-5.)*log(x) - x
c            + d9lgmc(x).
c***library   slatec (fnlib)
c***category  c7e
c***type      double precision (r9lgmc-s, d9lgmc-d, c9lgmc-c)
c***keywords  complete gamma function, correction term, fnlib,
c             log gamma, logarithm, special functions
c***author  fullerton, w., (lanl)
c***description
c
c compute the log gamma correction factor for x .ge. 10. so that
c log (dgamma(x)) = log(sqrt(2*pi)) + (x-.5)*log(x) - x + d9lgmc(x)
c
c series for algm       on the interval  0.          to  1.00000e-02
c                                        with weighted error   1.28e-31
c                                         log weighted error  30.89
c                               significant figures required  29.81
c                                    decimal places required  31.48
c
c***references  (none)
c***routines called  d1mach, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900720  routine changed from user-callable to subsidiary.  (wrb)
c***end prologue  d9lgmc
      double precision x, algmcs(15), xbig, xmax, dcsevl, d1mach
      logical first
      save algmcs, nalgm, xbig, xmax, first
      data algmcs(  1) / +.1666389480 4518632472 0572965082 2 d+0      /
      data algmcs(  2) / -.1384948176 0675638407 3298605913 5 d-4      /
      data algmcs(  3) / +.9810825646 9247294261 5717154748 7 d-8      /
      data algmcs(  4) / -.1809129475 5724941942 6330626671 9 d-10     /
      data algmcs(  5) / +.6221098041 8926052271 2601554341 6 d-13     /
      data algmcs(  6) / -.3399615005 4177219443 0333059966 6 d-15     /
      data algmcs(  7) / +.2683181998 4826987489 5753884666 6 d-17     /
      data algmcs(  8) / -.2868042435 3346432841 4462239999 9 d-19     /
      data algmcs(  9) / +.3962837061 0464348036 7930666666 6 d-21     /
      data algmcs( 10) / -.6831888753 9857668701 1199999999 9 d-23     /
      data algmcs( 11) / +.1429227355 9424981475 7333333333 3 d-24     /
      data algmcs( 12) / -.3547598158 1010705471 9999999999 9 d-26     /
      data algmcs( 13) / +.1025680058 0104709120 0000000000 0 d-27     /
      data algmcs( 14) / -.3401102254 3167487999 9999999999 9 d-29     /
      data algmcs( 15) / +.1276642195 6300629333 3333333333 3 d-30     /
      data first /.true./
c***first executable statement  d9lgmc
      if (first) then
         nalgm = initds (algmcs, 15, real(d1mach(3)) )
         xbig = 1.0d0/sqrt(d1mach(3))
         xmax = exp (min(log(d1mach(2)/12.d0), -log(12.d0*d1mach(1))))
      endif
      first = .false.
c
      if (x .lt. 10.d0) call xermsg ('slatec', 'd9lgmc',
     +   'x must be ge 10', 1, 2)
      if (x.ge.xmax) go to 20
c
      d9lgmc = 1.d0/(12.d0*x)
      if (x.lt.xbig) d9lgmc = dcsevl (2.0d0*(10.d0/x)**2-1.d0, algmcs,
     1  nalgm) / x
      return
c
 20   d9lgmc = 0.d0
      call xermsg ('slatec', 'd9lgmc', 'x so big d9lgmc underflows', 2,
     +   1)
      return
c
      end







c------------------------------------------------------------------
*deck dcsevl
      double precision function dcsevl (x, cs, n)
c***begin prologue  dcsevl
c***purpose  evaluate a chebyshev series.
c***library   slatec (fnlib)
c***category  c3a2
c***type      double precision (csevl-s, dcsevl-d)
c***keywords  chebyshev series, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c  evaluate the n-term chebyshev series cs at x.  adapted from
c  a method presented in the paper by broucke referenced below.
c
c       input arguments --
c  x    value at which the series is to be evaluated.
c  cs   array of n terms of a chebyshev series.  in evaluating
c       cs, only half the first coefficient is summed.
c  n    number of terms in array cs.
c
c***references  r. broucke, ten subroutines for the manipulation of
c                 chebyshev series, algorithm 446, communications of
c                 the a.c.m. 16, (1973) pp. 254-256.
c               l. fox and i. b. parker, chebyshev polynomials in
c                 numerical analysis, oxford university press, 1968,
c                 page 56.
c***routines called  d1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900329  prologued revised extensively and code rewritten to allow
c           x to be slightly outside interval (-1,+1).  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dcsevl
      double precision b0, b1, b2, cs(*), onepl, twox, x, d1mach
      logical first
      save first, onepl
      data first /.true./
c***first executable statement  dcsevl
      if (first) onepl = 1.0d0 + d1mach(4)
      first = .false.
      if (n .lt. 1) call xermsg ('slatec', 'dcsevl',
     +   'number of terms .le. 0', 2, 2)
      if (n .gt. 1000) call xermsg ('slatec', 'dcsevl',
     +   'number of terms .gt. 1000', 3, 2)
      if (abs(x) .gt. onepl) call xermsg ('slatec', 'dcsevl',
     +   'x outside the interval (-1,+1)', 1, 1)
c
      b1 = 0.0d0
      b0 = 0.0d0
      twox = 2.0d0*x
      do 10 i = 1,n
         b2 = b1
         b1 = b0
         ni = n + 1 - i
         b0 = twox*b1 - b2 + cs(ni)
   10 continue
c
      dcsevl = 0.5d0*(b0-b2)
c
      return
      end







c------------------------------------------------------------------
*deck dgamlm
      subroutine dgamlm (xmin, xmax)
c***begin prologue  dgamlm
c***purpose  compute the minimum and maximum bounds for the argument in
c            the gamma function.
c***library   slatec (fnlib)
c***category  c7a, r2
c***type      double precision (gamlim-s, dgamlm-d)
c***keywords  complete gamma function, fnlib, limits, special functions
c***author  fullerton, w., (lanl)
c***description
c
c calculate the minimum and maximum legal bounds for x in gamma(x).
c xmin and xmax are not the only bounds, but they are the only non-
c trivial ones to calculate.
c
c             output arguments --
c xmin   double precision minimum legal value of x in gamma(x).  any
c        smaller value of x might result in underflow.
c xmax   double precision maximum legal value of x in gamma(x).  any
c        larger value of x might cause overflow.
c
c***references  (none)
c***routines called  d1mach, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dgamlm
      double precision xmin, xmax, alnbig, alnsml, xln, xold, d1mach
c***first executable statement  dgamlm
      alnsml = log(d1mach(1))
      xmin = -alnsml
      do 10 i=1,10
        xold = xmin
        xln = log(xmin)
        xmin = xmin - xmin*((xmin+0.5d0)*xln - xmin - 0.2258d0 + alnsml)
     1    / (xmin*xln+0.5d0)
        if (abs(xmin-xold).lt.0.005d0) go to 20
 10   continue
      call xermsg ('slatec', 'dgamlm', 'unable to find xmin', 1, 2)
c
 20   xmin = -xmin + 0.01d0
c
      alnbig = log (d1mach(2))
      xmax = alnbig
      do 30 i=1,10
        xold = xmax
        xln = log(xmax)
        xmax = xmax - xmax*((xmax-0.5d0)*xln - xmax + 0.9189d0 - alnbig)
     1    / (xmax*xln-0.5d0)
        if (abs(xmax-xold).lt.0.005d0) go to 40
 30   continue
      call xermsg ('slatec', 'dgamlm', 'unable to find xmax', 2, 2)
c
 40   xmax = xmax - 0.01d0
      xmin = max (xmin, -xmax+1.d0)
c
      return
      end







c------------------------------------------------------------------
*deck dgamma
      double precision function dgamma (x)
c***begin prologue  dgamma
c***purpose  compute the complete gamma function.
c***library   slatec (fnlib)
c***category  c7a
c***type      double precision (gamma-s, dgamma-d, cgamma-c)
c***keywords  complete gamma function, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c dgamma(x) calculates the double precision complete gamma function
c for double precision argument x.
c
c series for gam        on the interval  0.          to  1.00000e+00
c                                        with weighted error   5.79e-32
c                                         log weighted error  31.24
c                               significant figures required  30.00
c                                    decimal places required  32.05
c
c***references  (none)
c***routines called  d1mach, d9lgmc, dcsevl, dgamlm, initds, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   890911  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920618  removed space from variable name.  (rwc, wrb)
c***end prologue  dgamma
      double precision x, gamcs(42), dxrel, pi, sinpiy, sq2pil, xmax,
     1  xmin, y, d9lgmc, dcsevl, d1mach
      logical first
c
      save gamcs, pi, sq2pil, ngam, xmin, xmax, dxrel, first
      data gamcs(  1) / +.8571195590 9893314219 2006239994 2 d-2      /
      data gamcs(  2) / +.4415381324 8410067571 9131577165 2 d-2      /
      data gamcs(  3) / +.5685043681 5993633786 3266458878 9 d-1      /
      data gamcs(  4) / -.4219835396 4185605010 1250018662 4 d-2      /
      data gamcs(  5) / +.1326808181 2124602205 8400679635 2 d-2      /
      data gamcs(  6) / -.1893024529 7988804325 2394702388 6 d-3      /
      data gamcs(  7) / +.3606925327 4412452565 7808221722 5 d-4      /
      data gamcs(  8) / -.6056761904 4608642184 8554829036 5 d-5      /
      data gamcs(  9) / +.1055829546 3022833447 3182350909 3 d-5      /
      data gamcs( 10) / -.1811967365 5423840482 9185589116 6 d-6      /
      data gamcs( 11) / +.3117724964 7153222777 9025459316 9 d-7      /
      data gamcs( 12) / -.5354219639 0196871408 7408102434 7 d-8      /
      data gamcs( 13) / +.9193275519 8595889468 8778682594 0 d-9      /
      data gamcs( 14) / -.1577941280 2883397617 6742327395 3 d-9      /
      data gamcs( 15) / +.2707980622 9349545432 6654043308 9 d-10     /
      data gamcs( 16) / -.4646818653 8257301440 8166105893 3 d-11     /
      data gamcs( 17) / +.7973350192 0074196564 6076717535 9 d-12     /
      data gamcs( 18) / -.1368078209 8309160257 9949917230 9 d-12     /
      data gamcs( 19) / +.2347319486 5638006572 3347177168 8 d-13     /
      data gamcs( 20) / -.4027432614 9490669327 6657053469 9 d-14     /
      data gamcs( 21) / +.6910051747 3721009121 3833697525 7 d-15     /
      data gamcs( 22) / -.1185584500 2219929070 5238712619 2 d-15     /
      data gamcs( 23) / +.2034148542 4963739552 0102605193 2 d-16     /
      data gamcs( 24) / -.3490054341 7174058492 7401294910 8 d-17     /
      data gamcs( 25) / +.5987993856 4853055671 3505106602 6 d-18     /
      data gamcs( 26) / -.1027378057 8722280744 9006977843 1 d-18     /
      data gamcs( 27) / +.1762702816 0605298249 4275966074 8 d-19     /
      data gamcs( 28) / -.3024320653 7353062609 5877211204 2 d-20     /
      data gamcs( 29) / +.5188914660 2183978397 1783355050 6 d-21     /
      data gamcs( 30) / -.8902770842 4565766924 4925160106 6 d-22     /
      data gamcs( 31) / +.1527474068 4933426022 7459689130 6 d-22     /
      data gamcs( 32) / -.2620731256 1873629002 5732833279 9 d-23     /
      data gamcs( 33) / +.4496464047 8305386703 3104657066 6 d-24     /
      data gamcs( 34) / -.7714712731 3368779117 0390152533 3 d-25     /
      data gamcs( 35) / +.1323635453 1260440364 8657271466 6 d-25     /
      data gamcs( 36) / -.2270999412 9429288167 0231381333 3 d-26     /
      data gamcs( 37) / +.3896418998 0039914493 2081663999 9 d-27     /
      data gamcs( 38) / -.6685198115 1259533277 9212799999 9 d-28     /
      data gamcs( 39) / +.1146998663 1400243843 4761386666 6 d-28     /
      data gamcs( 40) / -.1967938586 3451346772 9510399999 9 d-29     /
      data gamcs( 41) / +.3376448816 5853380903 3489066666 6 d-30     /
      data gamcs( 42) / -.5793070335 7821357846 2549333333 3 d-31     /
      data pi / 3.1415926535 8979323846 2643383279 50 d0 /
      data sq2pil / 0.9189385332 0467274178 0329736405 62 d0 /
      data first /.true./
c***first executable statement  dgamma
      if (first) then
         ngam = initds (gamcs, 42, 0.1*real(d1mach(3)) )
c
         call dgamlm (xmin, xmax)
         dxrel = sqrt(d1mach(4))
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.10.d0) go to 50
c
c compute gamma(x) for -xbnd .le. x .le. xbnd.  reduce interval and find
c gamma(1+y) for 0.0 .le. y .lt. 1.0 first of all.
c
      n = x
      if (x.lt.0.d0) n = n - 1
      y = x - n
      n = n - 1
      dgamma = 0.9375d0 + dcsevl (2.d0*y-1.d0, gamcs, ngam)
      if (n.eq.0) return
c
      if (n.gt.0) go to 30
c
c compute gamma(x) for x .lt. 1.0
c
      n = -n
      if (x .eq. 0.d0) call xermsg ('slatec', 'dgamma', 'x is 0', 4, 2)
      if (x .lt. 0.0 .and. x+n-2 .eq. 0.d0) call xermsg ('slatec',
     +   'dgamma', 'x is a negative integer', 4, 2)
      if (x .lt. (-0.5d0) .and. abs((x-aint(x-0.5d0))/x) .lt. dxrel)
     +   call xermsg ('slatec', 'dgamma',
     +   'answer lt half precision because x too near negative integer',
     +   1, 1)
c
      do 20 i=1,n
        dgamma = dgamma/(x+i-1 )
 20   continue
      return
c
c gamma(x) for x .ge. 2.0 and x .le. 10.0
c
 30   do 40 i=1,n
        dgamma = (y+i) * dgamma
 40   continue
      return
c
c gamma(x) for abs(x) .gt. 10.0.  recall y = abs(x).
c
 50   if (x .gt. xmax) call xermsg ('slatec', 'dgamma',
     +   'x so big gamma overflows', 3, 2)
c
      dgamma = 0.d0
      if (x .lt. xmin) call xermsg ('slatec', 'dgamma',
     +   'x so small gamma underflows', 2, 1)
      if (x.lt.xmin) return
c
      dgamma = exp ((y-0.5d0)*log(y) - y + sq2pil + d9lgmc(y) )
      if (x.gt.0.d0) return
c
      if (abs((x-aint(x-0.5d0))/x) .lt. dxrel) call xermsg ('slatec',
     +   'dgamma',
     +   'answer lt half precision, x too near negative integer', 1, 1)
c
      sinpiy = sin (pi*y)
      if (sinpiy .eq. 0.d0) call xermsg ('slatec', 'dgamma',
     +   'x is a negative integer', 4, 2)
c
      dgamma = -pi/(y*sinpiy*dgamma)
c
      return
      end







c------------------------------------------------------------------
*deck dlgams
      subroutine dlgams (x, dlgam, sgngam)
c***begin prologue  dlgams
c***purpose  compute the logarithm of the absolute value of the gamma
c            function.
c***library   slatec (fnlib)
c***category  c7a
c***type      double precision (algams-s, dlgams-d)
c***keywords  absolute value of the logarithm of the gamma function,
c             fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c dlgams(x,dlgam,sgngam) calculates the double precision natural
c logarithm of the absolute value of the gamma function for
c double precision argument x and stores the result in double
c precision argument dlgam.
c
c***references  (none)
c***routines called  dlngam
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  dlgams
      double precision x, dlgam, sgngam, dlngam
c***first executable statement  dlgams
      dlgam = dlngam(x)
      sgngam = 1.0d0
      if (x.gt.0.d0) return
c
      int = mod (-aint(x), 2.0d0) + 0.1d0
      if (int.eq.0) sgngam = -1.0d0
c
      return
      end







c------------------------------------------------------------------
*deck dlngam
      double precision function dlngam (x)
c***begin prologue  dlngam
c***purpose  compute the logarithm of the absolute value of the gamma
c            function.
c***library   slatec (fnlib)
c***category  c7a
c***type      double precision (alngam-s, dlngam-d, clngam-c)
c***keywords  absolute value, complete gamma function, fnlib, logarithm,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c dlngam(x) calculates the double precision logarithm of the
c absolute value of the gamma function for double precision
c argument x.
c
c***references  (none)
c***routines called  d1mach, d9lgmc, dgamma, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900727  added external statement.  (wrb)
c***end prologue  dlngam
      double precision x, dxrel, pi, sinpiy, sqpi2l, sq2pil, xmax,
     1  y, dgamma, d9lgmc, d1mach, temp
      logical first
      external dgamma
      save sq2pil, sqpi2l, pi, xmax, dxrel, first
      data sq2pil / 0.9189385332 0467274178 0329736405 62 d0 /
      data sqpi2l / +.2257913526 4472743236 3097614947 441 d+0    /
      data pi / 3.1415926535 8979323846 2643383279 50 d0 /
      data first /.true./
c***first executable statement  dlngam
      if (first) then
         temp = 1.d0/log(d1mach(2))
         xmax = temp*d1mach(2)
         dxrel = sqrt(d1mach(4))
      endif
      first = .false.
c
      y = abs (x)
      if (y.gt.10.d0) go to 20
c
c log (abs (dgamma(x)) ) for abs(x) .le. 10.0
c
      dlngam = log (abs (dgamma(x)) )
      return
c
c log ( abs (dgamma(x)) ) for abs(x) .gt. 10.0
c
 20   if (y .gt. xmax) call xermsg ('slatec', 'dlngam',
     +   'abs(x) so big dlngam overflows', 2, 2)
c
      if (x.gt.0.d0) dlngam = sq2pil + (x-0.5d0)*log(x) - x + d9lgmc(y)
      if (x.gt.0.d0) return
c
      sinpiy = abs (sin(pi*y))
      if (sinpiy .eq. 0.d0) call xermsg ('slatec', 'dlngam',
     +   'x is a negative integer', 3, 2)
c
      if (abs((x-aint(x-0.5d0))/x) .lt. dxrel) call xermsg ('slatec',
     +   'dlngam',
     +   'answer lt half precision because x too near negative integer',
     +   1, 1)
c
      dlngam = sqpi2l + (x-0.5d0)*log(y) - x - log(sinpiy) - d9lgmc(y)
      return
c
      end







c------------------------------------------------------------------
*deck initds
      function initds (os, nos, eta)
c***begin prologue  initds
c***purpose  determine the number of terms needed in an orthogonal
c            polynomial series so that it meets a specified accuracy.
c***library   slatec (fnlib)
c***category  c3a2
c***type      double precision (inits-s, initds-d)
c***keywords  chebyshev, fnlib, initialize, orthogonal polynomial,
c             orthogonal series, special functions
c***author  fullerton, w., (lanl)
c***description
c
c  initialize the orthogonal series, represented by the array os, so
c  that initds is the number of terms needed to insure the error is no
c  larger than eta.  ordinarily, eta will be chosen to be one-tenth
c  machine precision.
c
c             input arguments --
c   os     double precision array of nos coefficients in an orthogonal
c          series.
c   nos    number of coefficients in os.
c   eta    single precision scalar containing requested accuracy of
c          series.
c
c***references  (none)
c***routines called  xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891115  modified error message.  (wrb)
c   891115  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  initds
      double precision os(*)
c***first executable statement  initds
      if (nos .lt. 1) call xermsg ('slatec', 'initds',
     +   'number of coefficients is less than 1', 2, 1)
c
      err = 0.
      do 10 ii = 1,nos
        i = nos + 1 - ii
        err = err + abs(real(os(i)))
        if (err.gt.eta) go to 20
   10 continue
c
   20 if (i .eq. nos) call xermsg ('slatec', 'initds',
     +   'chebyshev series too short for specified accuracy', 1, 1)
      initds = i
c
      return
      end







c------------------------------------------------------------------
c     secondary dependency subroutines and god knows what
c------------------------------------------------------------------
*deck xermsg
      subroutine xermsg (librar, subrou, messg, nerr, level)
c***begin prologue  xermsg
c***purpose  process error messages for slatec and other libraries.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xermsg-a)
c***keywords  error message, xerror
c***author  fong, kirby, (nmfecc at llnl)
c***description
c
c   xermsg processes a diagnostic message in a manner determined by the
c   value of level and the current value of the library error control
c   flag, kontrl.  see subroutine xsetf for details.
c
c    librar   a character constant (or character variable) with the name
c             of the library.  this will be 'slatec' for the slatec
c             common math library.  the error handling package is
c             general enough to be used by many libraries
c             simultaneously, so it is desirable for the routine that
c             detects and reports an error to identify the library name
c             as well as the routine name.
c
c    subrou   a character constant (or character variable) with the name
c             of the routine that detected the error.  usually it is the
c             name of the routine that is calling xermsg.  there are
c             some instances where a user callable library routine calls
c             lower level subsidiary routines where the error is
c             detected.  in such cases it may be more informative to
c             supply the name of the routine the user called rather than
c             the name of the subsidiary routine that detected the
c             error.
c
c    messg    a character constant (or character variable) with the text
c             of the error or warning message.  in the example below,
c             the message is a character constant that contains a
c             generic message.
c
c                   call xermsg ('slatec', 'mmpy',
c                  *'the order of the matrix exceeds the row dimension',
c                  *3, 1)
c
c             it is possible (and is sometimes desirable) to generate a
c             specific message--e.g., one that contains actual numeric
c             values.  specific numeric values can be converted into
c             character strings using formatted write statements into
c             character variables.  this is called standard fortran
c             internal file i/o and is exemplified in the first three
c             lines of the following example.  you can also catenate
c             substrings of characters to construct the error message.
c             here is an example showing the use of both writing to
c             an internal file and catenating character strings.
c
c                   character*5 charn, charl
c                   write (charn,10) n
c                   write (charl,10) lda
c                10 format(i5)
c                   call xermsg ('slatec', 'mmpy', 'the order'//charn//
c                  *   ' of the matrix exceeds its row dimension of'//
c                  *   charl, 3, 1)
c
c             there are two subtleties worth mentioning.  one is that
c             the // for character catenation is used to construct the
c             error message so that no single character constant is
c             continued to the next line.  this avoids confusion as to
c             whether there are trailing blanks at the end of the line.
c             the second is that by catenating the parts of the message
c             as an actual argument rather than encoding the entire
c             message into one large character variable, we avoid
c             having to know how long the message will be in order to
c             declare an adequate length for that large character
c             variable.  xermsg calls xerprn to print the message using
c             multiple lines if necessary.  if the message is very long,
c             xerprn will break it into pieces of 72 characters (as
c             requested by xermsg) for printing on multiple lines.
c             also, xermsg asks xerprn to prefix each line with ' *  '
c             so that the total line length could be 76 characters.
c             note also that xerprn scans the error message backwards
c             to ignore trailing blanks.  another feature is that
c             the substring '$$' is treated as a new line sentinel
c             by xerprn.  if you want to construct a multiline
c             message without having to count out multiples of 72
c             characters, just use '$$' as a separator.  '$$'
c             obviously must occur within 72 characters of the
c             start of each line to have its intended effect since
c             xerprn is asked to wrap around at 72 characters in
c             addition to looking for '$$'.
c
c    nerr     an integer value that is chosen by the library routine's
c             author.  it must be in the range -99 to 999 (three
c             printable digits).  each distinct error should have its
c             own error number.  these error numbers should be described
c             in the machine readable documentation for the routine.
c             the error numbers need be unique only within each routine,
c             so it is reasonable for each routine to start enumerating
c             errors from 1 and proceeding to the next integer.
c
c    level    an integer value in the range 0 to 2 that indicates the
c             level (severity) of the error.  their meanings are
c
c            -1  a warning message.  this is used if it is not clear
c                that there really is an error, but the user's attention
c                may be needed.  an attempt is made to only print this
c                message once.
c
c             0  a warning message.  this is used if it is not clear
c                that there really is an error, but the user's attention
c                may be needed.
c
c             1  a recoverable error.  this is used even if the error is
c                so serious that the routine cannot return any useful
c                answer.  if the user has told the error package to
c                return after recoverable errors, then xermsg will
c                return to the library routine which can then return to
c                the user's routine.  the user may also permit the error
c                package to terminate the program upon encountering a
c                recoverable error.
c
c             2  a fatal error.  xermsg will not return to its caller
c                after it receives a fatal error.  this level should
c                hardly ever be used; it is much better to allow the
c                user a chance to recover.  an example of one of the few
c                cases in which it is permissible to declare a level 2
c                error is a reverse communication library routine that
c                is likely to be called repeatedly until it integrates
c                across some interval.  if there is a serious error in
c                the input such that another step cannot be taken and
c                the library routine is called again without the input
c                error having been corrected by the caller, the library
c                routine will probably be called forever with improper
c                input.  in this case, it is reasonable to declare the
c                error to be fatal.
c
c    each of the arguments to xermsg is input; none will be modified by
c    xermsg.  a routine may make multiple calls to xermsg with warning
c    level messages; however, after a call to xermsg with a recoverable
c    error, the routine should return to the user.  do not try to call
c    xermsg with a second recoverable error after the first recoverable
c    error because the error package saves the error number.  the user
c    can retrieve this error number by calling another entry point in
c    the error handling package and then clear the error number when
c    recovering from the error.  calling xermsg in succession causes the
c    old error number to be overwritten by the latest error number.
c    this is considered harmless for error numbers associated with
c    warning messages but must not be done for error numbers of serious
c    errors.  after a call to xermsg with a recoverable error, the user
c    must be given a chance to call numxer or xerclr to retrieve or
c    clear the error number.
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  fdump, j4save, xercnt, xerhlt, xerprn, xersve
c***revision history  (yymmdd)
c   880101  date written
c   880621  revised as directed at slatec cml meeting of february 1988.
c           there are two basic changes.
c           1.  a new routine, xerprn, is used instead of xerprt to
c               print messages.  this routine will break long messages
c               into pieces for printing on multiple lines.  '$$' is
c               accepted as a new line sentinel.  a prefix can be
c               added to each line to be printed.  xermsg uses either
c               ' ***' or ' *  ' and long messages are broken every
c               72 characters (at most) so that the maximum line
c               length output can now be as great as 76.
c           2.  the text of all messages is now in upper case since the
c               fortran standard document does not admit the existence
c               of lower case.
c   880708  revised after the slatec cml meeting of june 29 and 30.
c           the principal changes are
c           1.  clarify comments in the prologues
c           2.  rename xrprnt to xerprn
c           3.  rework handling of '$$' in xerprn to handle blank lines
c               similar to the way format statements handle the /
c               character for new records.
c   890706  revised with the help of fred fritsch and reg clemens to
c           clean up the coding.
c   890721  revised to use new feature in xerprn to count characters in
c           prefix.
c   891013  revised to correct comments.
c   891214  prologue converted to version 4.0 format.  (wrb)
c   900510  changed test on nerr to be -9999999 < nerr < 99999999, but
c           nerr .ne. 0, and on level to be -2 < level < 3.  added
c           level=-1 logic, changed calls to xersav to xersve, and
c           xerctl to xercnt.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xermsg
      character*(*) librar, subrou, messg
      character*8 xlibr, xsubr
      character*72  temp
      character*20  lfirst
c***first executable statement  xermsg
      lkntrl = j4save (2, 0, .false.)
      maxmes = j4save (4, 0, .false.)
c
c       lkntrl is a local copy of the control flag kontrl.
c       maxmes is the maximum number of times any particular message
c          should be printed.
c
c       we print a fatal error message and terminate for an error in
c          calling xermsg.  the error number should be positive,
c          and the level should be between 0 and 2.
c
      if (nerr.lt.-9999999 .or. nerr.gt.99999999 .or. nerr.eq.0 .or.
     *   level.lt.-1 .or. level.gt.2) then
         call xerprn (' ***', -1, 'fatal error in...$$ ' //
     *      'xermsg -- invalid error number or level$$ '//
     *      'job abort due to fatal error.', 72)
         call xersve (' ', ' ', ' ', 0, 0, 0, kdummy)
         call xerhlt (' ***xermsg -- invalid input')
         return
      endif
c
c       record the message.
c
      i = j4save (1, nerr, .true.)
      call xersve (librar, subrou, messg, 1, nerr, level, kount)
c
c       handle print-once warning messages.
c
      if (level.eq.-1 .and. kount.gt.1) return
c
c       allow temporary user override of the control flag.
c
      xlibr  = librar
      xsubr  = subrou
      lfirst = messg
      lerr   = nerr
      llevel = level
      call xercnt (xlibr, xsubr, lfirst, lerr, llevel, lkntrl)
c
      lkntrl = max(-2, min(2,lkntrl))
      mkntrl = abs(lkntrl)
c
c       skip printing if the control flag value as reset in xercnt is
c       zero and the error is not fatal.
c
      if (level.lt.2 .and. lkntrl.eq.0) go to 30
      if (level.eq.0 .and. kount.gt.maxmes) go to 30
      if (level.eq.1 .and. kount.gt.maxmes .and. mkntrl.eq.1) go to 30
      if (level.eq.2 .and. kount.gt.max(1,maxmes)) go to 30
c
c       announce the names of the library and subroutine by building a
c       message in character variable temp (not exceeding 66 characters)
c       and sending it out via xerprn.  print only if control flag
c       is not zero.
c
      if (lkntrl .ne. 0) then
         temp(1:21) = 'message from routine '
         i = min(len(subrou), 16)
         temp(22:21+i) = subrou(1:i)
         temp(22+i:33+i) = ' in library '
         ltemp = 33 + i
         i = min(len(librar), 16)
         temp(ltemp+1:ltemp+i) = librar (1:i)
         temp(ltemp+i+1:ltemp+i+1) = '.'
         ltemp = ltemp + i + 1
         call xerprn (' ***', -1, temp(1:ltemp), 72)
      endif
c
c       if lkntrl is positive, print an introductory line before
c       printing the message.  the introductory line tells the choice
c       from each of the following three options.
c       1.  level of the message
c              'informative message'
c              'potentially recoverable error'
c              'fatal error'
c       2.  whether control flag will allow program to continue
c              'prog continues'
c              'prog aborted'
c       3.  whether or not a traceback was requested.  (the traceback
c           may not be implemented at some sites, so this only tells
c           what was requested, not what was delivered.)
c              'traceback requested'
c              'traceback not requested'
c       notice that the line including four prefix characters will not
c       exceed 74 characters.
c       we skip the next block if the introductory line is not needed.
c
      if (lkntrl .gt. 0) then
c
c       the first part of the message tells about the level.
c
         if (level .le. 0) then
            temp(1:20) = 'informative message,'
            ltemp = 20
         elseif (level .eq. 1) then
            temp(1:30) = 'potentially recoverable error,'
            ltemp = 30
         else
            temp(1:12) = 'fatal error,'
            ltemp = 12
         endif
c
c       then whether the program will continue.
c
         if ((mkntrl.eq.2 .and. level.ge.1) .or.
     *       (mkntrl.eq.1 .and. level.eq.2)) then
            temp(ltemp+1:ltemp+14) = ' prog aborted,'
            ltemp = ltemp + 14
         else
            temp(ltemp+1:ltemp+16) = ' prog continues,'
            ltemp = ltemp + 16
         endif
c
c       finally tell whether there should be a traceback.
c
         if (lkntrl .gt. 0) then
            temp(ltemp+1:ltemp+20) = ' traceback requested'
            ltemp = ltemp + 20
         else
            temp(ltemp+1:ltemp+24) = ' traceback not requested'
            ltemp = ltemp + 24
         endif
         call xerprn (' ***', -1, temp(1:ltemp), 72)
      endif
c
c       now send out the message.
c
      call xerprn (' *  ', -1, messg, 72)
c
c       if lkntrl is positive, write the error number and request a
c          traceback.
c
      if (lkntrl .gt. 0) then
         write (temp, '(''error number = '', i8)') nerr
         do 10 i=16,22
            if (temp(i:i) .ne. ' ') go to 20
   10    continue
c
   20    call xerprn (' *  ', -1, temp(1:15) // temp(i:23), 72)
         call fdump
      endif
c
c       if lkntrl is not zero, print a blank line and an end of message.
c
      if (lkntrl .ne. 0) then
         call xerprn (' *  ', -1, ' ', 72)
         call xerprn (' ***', -1, 'end of message', 72)
         call xerprn ('    ',  0, ' ', 72)
      endif
c
c       if the error is not fatal or the error is recoverable and the
c       control flag is set for recovery, then return.
c
   30 if (level.le.0 .or. (level.eq.1 .and. mkntrl.le.1)) return
c
c       the program will be stopped due to an unrecovered error or a
c       fatal error.  print the reason for the abort and the error
c       summary if the control flag and the maximum error count permit.
c
      if (lkntrl.gt.0 .and. kount.lt.max(1,maxmes)) then
         if (level .eq. 1) then
            call xerprn
     *         (' ***', -1, 'job abort due to unrecovered error.', 72)
         else
            call xerprn(' ***', -1, 'job abort due to fatal error.', 72)
         endif
         call xersve (' ', ' ', ' ', -1, 0, 0, kdummy)
         call xerhlt (' ')
      else
         call xerhlt (messg)
      endif
      return
      end







c-------------------------------------------------------------------
*deck xerclr
      subroutine xerclr
c***begin prologue  xerclr
c***purpose  reset current error number to zero.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xerclr-a)
c***keywords  error, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        this routine simply resets the current error number to zero.
c        this may be necessary in order to determine that a certain
c        error has occurred again since the last time numxer was
c        referenced.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  j4save
c***revision history  (yymmdd)
c   790801  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xerclr
c***first executable statement  xerclr
      junk = j4save(1,0,.true.)
      return
      end







c-------------------------------------------------------------------
*deck fdump
      subroutine fdump
c***begin prologue  fdump
c***purpose  symbolic dump (should be locally written).
c***library   slatec (xerror)
c***category  r3
c***type      all (fdump-a)
c***keywords  error, xermsg
c***author  jones, r. e., (snla)
c***description
c
c        ***note*** machine dependent routine
c        fdump is intended to be replaced by a locally written
c        version which produces a symbolic dump.  failing this,
c        it should be replaced by a version which prints the
c        subprogram nesting list.  note that this dump must be
c        printed on each of up to five files, as indicated by the
c        xgetua routine.  see xsetua and xgetua for details.
c
c     written by ron jones, with slatec common math library subcommittee
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   790801  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  fdump
c***first executable statement  fdump
      return
      end






c-------------------------------------------------------------------
*deck j4save
      function j4save (iwhich, ivalue, iset)
c***begin prologue  j4save
c***subsidiary
c***purpose  save or recall global variables needed by error
c            handling routines.
c***library   slatec (xerror)
c***type      integer (j4save-i)
c***keywords  error messages, error number, recall, save, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        j4save saves and recalls several global variables needed
c        by the library error handling routines.
c
c     description of parameters
c      --input--
c        iwhich - index of item desired.
c                = 1 refers to current error number.
c                = 2 refers to current error control flag.
c                = 3 refers to current unit number to which error
c                    messages are to be sent.  (0 means use standard.)
c                = 4 refers to the maximum number of times any
c                     message is to be printed (as set by xermax).
c                = 5 refers to the total number of units to which
c                     each error message is to be written.
c                = 6 refers to the 2nd unit for error messages
c                = 7 refers to the 3rd unit for error messages
c                = 8 refers to the 4th unit for error messages
c                = 9 refers to the 5th unit for error messages
c        ivalue - the value to be set for the iwhich-th parameter,
c                 if iset is .true. .
c        iset   - if iset=.true., the iwhich-th parameter will be
c                 given the value, ivalue.  if iset=.false., the
c                 iwhich-th parameter will be unchanged, and ivalue
c                 is a dummy parameter.
c      --output--
c        the (old) value of the iwhich-th parameter will be returned
c        in the function value, j4save.
c
c***see also  xermsg
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  (none)
c***revision history  (yymmdd)
c   790801  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900205  minor modifications to prologue.  (wrb)
c   900402  added type section.  (wrb)
c   910411  added keywords section.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  j4save
      logical iset
      integer iparam(9)
      save iparam
      data iparam(1),iparam(2),iparam(3),iparam(4)/0,2,0,10/
      data iparam(5)/1/
      data iparam(6),iparam(7),iparam(8),iparam(9)/0,0,0,0/
c***first executable statement  j4save
      j4save = iparam(iwhich)
      if (iset) iparam(iwhich) = ivalue
      return
      end







c-------------------------------------------------------------------
*deck xercnt
      subroutine xercnt (librar, subrou, messg, nerr, level, kontrl)
c***begin prologue  xercnt
c***subsidiary
c***purpose  allow user control over handling of errors.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xercnt-a)
c***keywords  error, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        allows user control over handling of individual errors.
c        just after each message is recorded, but before it is
c        processed any further (i.e., before it is printed or
c        a decision to abort is made), a call is made to xercnt.
c        if the user has provided his own version of xercnt, he
c        can then override the value of kontrol used in processing
c        this message by redefining its value.
c        kontrl may be set to any value from -2 to 2.
c        the meanings for kontrl are the same as in xsetf, except
c        that the value of kontrl changes only for this message.
c        if kontrl is set to a value outside the range from -2 to 2,
c        it will be moved back into that range.
c
c     description of parameters
c
c      --input--
c        librar - the library that the routine is in.
c        subrou - the subroutine that xermsg is being called from
c        messg  - the first 20 characters of the error message.
c        nerr   - same as in the call to xermsg.
c        level  - same as in the call to xermsg.
c        kontrl - the current value of the control flag as set
c                 by a call to xsetf.
c
c      --output--
c        kontrl - the new value of kontrl.  if kontrl is not
c                 defined, it will remain at its original value.
c                 this changed value of control affects only
c                 the current occurrence of the current message.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  (none)
c***revision history  (yymmdd)
c   790801  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900206  routine changed from user-callable to subsidiary.  (wrb)
c   900510  changed calling sequence to include library and subroutine
c           names, changed routine name from xerctl to xercnt.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xercnt
      character*(*) librar, subrou, messg
c***first executable statement  xercnt
      return
      end







c-------------------------------------------------------------------
*deck xerhlt
      subroutine xerhlt (messg)
c***begin prologue  xerhlt
c***subsidiary
c***purpose  abort program execution and print error message.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xerhlt-a)
c***keywords  abort program execution, error, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        ***note*** machine dependent routine
c        xerhlt aborts the execution of the program.
c        the error message causing the abort is given in the calling
c        sequence, in case one needs it for printing on a dayfile,
c        for example.
c
c     description of parameters
c        messg is as in xermsg.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  (none)
c***revision history  (yymmdd)
c   790801  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900206  routine changed from user-callable to subsidiary.  (wrb)
c   900510  changed calling sequence to delete length of character
c           and changed routine name from xerabt to xerhlt.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xerhlt
      character*(*) messg
c***first executable statement  xerhlt
      stop
      end







c-------------------------------------------------------------------
*deck xerprn
      subroutine xerprn (prefix, npref, messg, nwrap)
c***begin prologue  xerprn
c***subsidiary
c***purpose  print error messages processed by xermsg.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xerprn-a)
c***keywords  error messages, printing, xerror
c***author  fong, kirby, (nmfecc at llnl)
c***description
c
c this routine sends one or more lines to each of the (up to five)
c logical units to which error messages are to be sent.  this routine
c is called several times by xermsg, sometimes with a single line to
c print and sometimes with a (potentially very long) message that may
c wrap around into multiple lines.
c
c prefix  input argument of type character.  this argument contains
c         characters to be put at the beginning of each line before
c         the body of the message.  no more than 16 characters of
c         prefix will be used.
c
c npref   input argument of type integer.  this argument is the number
c         of characters to use from prefix.  if it is negative, the
c         intrinsic function len is used to determine its length.  if
c         it is zero, prefix is not used.  if it exceeds 16 or if
c         len(prefix) exceeds 16, only the first 16 characters will be
c         used.  if npref is positive and the length of prefix is less
c         than npref, a copy of prefix extended with blanks to length
c         npref will be used.
c
c messg   input argument of type character.  this is the text of a
c         message to be printed.  if it is a long message, it will be
c         broken into pieces for printing on multiple lines.  each line
c         will start with the appropriate prefix and be followed by a
c         piece of the message.  nwrap is the number of characters per
c         piece; that is, after each nwrap characters, we break and
c         start a new line.  in addition the characters '$$' embedded
c         in messg are a sentinel for a new line.  the counting of
c         characters up to nwrap starts over for each new line.  the
c         value of nwrap typically used by xermsg is 72 since many
c         older error messages in the slatec library are laid out to
c         rely on wrap-around every 72 characters.
c
c nwrap   input argument of type integer.  this gives the maximum size
c         piece into which to break messg for printing on multiple
c         lines.  an embedded '$$' ends a line, and the count restarts
c         at the following character.  if a line break does not occur
c         on a blank (it would split a word) that word is moved to the
c         next line.  values of nwrap less than 16 will be treated as
c         16.  values of nwrap greater than 132 will be treated as 132.
c         the actual line length will be npref + nwrap after npref has
c         been adjusted to fall between 0 and 16 and nwrap has been
c         adjusted to fall between 16 and 132.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  i1mach, xgetua
c***revision history  (yymmdd)
c   880621  date written
c   880708  revised after the slatec cml subcommittee meeting of
c           june 29 and 30 to change the name to xerprn and to rework
c           the handling of the new line sentinel to behave like the
c           slash character in format statements.
c   890706  revised with the help of fred fritsch and reg clemens to
c           streamline the coding and fix a bug that caused extra blank
c           lines to be printed.
c   890721  revised to add a new feature.  a negative value of npref
c           causes len(prefix) to be used as the length.
c   891013  revised to correct error in calculating prefix length.
c   891214  prologue converted to version 4.0 format.  (wrb)
c   900510  added code to break messages between words.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xerprn
      character*(*) prefix, messg
      integer npref, nwrap
      character*148 cbuff
      integer iu(5), nunit
      character*2 newlin
      parameter (newlin = '$$')
c***first executable statement  xerprn
      call xgetua(iu,nunit)
c
c       a zero value for a logical unit number means to use the standard
c       error message unit instead.  i1mach(4) retrieves the standard
c       error message unit.
c
      n = i1mach(4)
      do 10 i=1,nunit
         if (iu(i) .eq. 0) iu(i) = n
   10 continue
c
c       lpref is the length of the prefix.  the prefix is placed at the
c       beginning of cbuff, the character buffer, and kept there during
c       the rest of this routine.
c
      if ( npref .lt. 0 ) then
         lpref = len(prefix)
      else
         lpref = npref
      endif
      lpref = min(16, lpref)
      if (lpref .ne. 0) cbuff(1:lpref) = prefix
c
c       lwrap is the maximum number of characters we want to take at one
c       time from messg to print on one line.
c
      lwrap = max(16, min(132, nwrap))
c
c       set lenmsg to the length of messg, ignore any trailing blanks.
c
      lenmsg = len(messg)
      n = lenmsg
      do 20 i=1,n
         if (messg(lenmsg:lenmsg) .ne. ' ') go to 30
         lenmsg = lenmsg - 1
   20 continue
   30 continue
c
c       if the message is all blanks, then print one blank line.
c
      if (lenmsg .eq. 0) then
         cbuff(lpref+1:lpref+1) = ' '
         do 40 i=1,nunit
            write(iu(i), '(a)') cbuff(1:lpref+1)
   40    continue
         return
      endif
c
c       set nextc to the position in messg where the next substring
c       starts.  from this position we scan for the new line sentinel.
c       when nextc exceeds lenmsg, there is no more to print.
c       we loop back to label 50 until all pieces have been printed.
c
c       we look for the next occurrence of the new line sentinel.  the
c       index intrinsic function returns zero if there is no occurrence
c       or if the length of the first argument is less than the length
c       of the second argument.
c
c       there are several cases which should be checked for in the
c       following order.  we are attempting to set lpiece to the number
c       of characters that should be taken from messg starting at
c       position nextc.
c
c       lpiece .eq. 0   the new line sentinel does not occur in the
c                       remainder of the character string.  lpiece
c                       should be set to lwrap or lenmsg+1-nextc,
c                       whichever is less.
c
c       lpiece .eq. 1   the new line sentinel starts at messg(nextc:
c                       nextc).  lpiece is effectively zero, and we
c                       print nothing to avoid producing unnecessary
c                       blank lines.  this takes care of the situation
c                       where the library routine has a message of
c                       exactly 72 characters followed by a new line
c                       sentinel followed by more characters.  nextc
c                       should be incremented by 2.
c
c       lpiece .gt. lwrap+1  reduce lpiece to lwrap.
c
c       else            this last case means 2 .le. lpiece .le. lwrap+1
c                       reset lpiece = lpiece-1.  note that this
c                       properly handles the end case where lpiece .eq.
c                       lwrap+1.  that is, the sentinel falls exactly
c                       at the end of a line.
c
      nextc = 1
   50 lpiece = index(messg(nextc:lenmsg), newlin)
      if (lpiece .eq. 0) then
c
c       there was no new line sentinel found.
c
         idelta = 0
         lpiece = min(lwrap, lenmsg+1-nextc)
         if (lpiece .lt. lenmsg+1-nextc) then
            do 52 i=lpiece+1,2,-1
               if (messg(nextc+i-1:nextc+i-1) .eq. ' ') then
                  lpiece = i-1
                  idelta = 1
                  goto 54
               endif
   52       continue
         endif
   54    cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
         nextc = nextc + lpiece + idelta
      elseif (lpiece .eq. 1) then
c
c       we have a new line sentinel at messg(nextc:nextc+1).
c       don't print a blank line.
c
         nextc = nextc + 2
         go to 50
      elseif (lpiece .gt. lwrap+1) then
c
c       lpiece should be set down to lwrap.
c
         idelta = 0
         lpiece = lwrap
         do 56 i=lpiece+1,2,-1
            if (messg(nextc+i-1:nextc+i-1) .eq. ' ') then
               lpiece = i-1
               idelta = 1
               goto 58
            endif
   56    continue
   58    cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
         nextc = nextc + lpiece + idelta
      else
c
c       if we arrive here, it means 2 .le. lpiece .le. lwrap+1.
c       we should decrement lpiece by one.
c
         lpiece = lpiece - 1
         cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
         nextc  = nextc + lpiece + 2
      endif
c
c       print
c
      do 60 i=1,nunit
         write(iu(i), '(a)') cbuff(1:lpref+lpiece)
   60 continue
c
      if (nextc .le. lenmsg) go to 50
      return
      end







c-------------------------------------------------------------------
*deck xersve
      subroutine xersve (librar, subrou, messg, kflag, nerr, level,
     +   icount)
c***begin prologue  xersve
c***subsidiary
c***purpose  record that an error has occurred.
c***library   slatec (xerror)
c***category  r3
c***type      all (xersve-a)
c***keywords  error, xerror
c***author  jones, r. e., (snla)
c***description
c
c *usage:
c
c        integer  kflag, nerr, level, icount
c        character * (len) librar, subrou, messg
c
c        call xersve (librar, subrou, messg, kflag, nerr, level, icount)
c
c *arguments:
c
c        librar :in    is the library that the message is from.
c        subrou :in    is the subroutine that the message is from.
c        messg  :in    is the message to be saved.
c        kflag  :in    indicates the action to be performed.
c                      when kflag > 0, the message in messg is saved.
c                      when kflag=0 the tables will be dumped and
c                      cleared.
c                      when kflag < 0, the tables will be dumped and
c                      not cleared.
c        nerr   :in    is the error number.
c        level  :in    is the error severity.
c        icount :out   the number of times this message has been seen,
c                      or zero if the table has overflowed and does not
c                      contain this message specifically.  when kflag=0,
c                      icount will not be altered.
c
c *description:
c
c   record that this error occurred and possibly dump and clear the
c   tables.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  i1mach, xgetua
c***revision history  (yymmdd)
c   800319  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900413  routine modified to remove reference to kflag.  (wrb)
c   900510  changed to add library name and subroutine to calling
c           sequence, use if-then-else, make number of saved entries
c           easily changeable, changed routine name from xersav to
c           xersve.  (rwc)
c   910626  added libtab and subtab to save statement.  (bks)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xersve
      parameter (lentab=10)
      integer lun(5)
      character*(*) librar, subrou, messg
      character*8  libtab(lentab), subtab(lentab), lib, sub
      character*20 mestab(lentab), mes
      dimension nertab(lentab), levtab(lentab), kount(lentab)
      save libtab, subtab, mestab, nertab, levtab, kount, kountx, nmsg
      data kountx/0/, nmsg/0/
c***first executable statement  xersve
c
      if (kflag.le.0) then
c
c        dump the table.
c
         if (nmsg.eq.0) return
c
c        print to each unit.
c
         call xgetua (lun, nunit)
         do 20 kunit = 1,nunit
            iunit = lun(kunit)
            if (iunit.eq.0) iunit = i1mach(4)
c
c           print the table header.
c
            write (iunit,9000)
c
c           print body of table.
c
            do 10 i = 1,nmsg
               write (iunit,9010) libtab(i), subtab(i), mestab(i),
     *            nertab(i),levtab(i),kount(i)
   10       continue
c
c           print number of other errors.
c
            if (kountx.ne.0) write (iunit,9020) kountx
            write (iunit,9030)
   20    continue
c
c        clear the error tables.
c
         if (kflag.eq.0) then
            nmsg = 0
            kountx = 0
         endif
      else
c
c        process a message...
c        search for this messg, or else an empty slot for this messg,
c        or else determine that the error table is full.
c
         lib = librar
         sub = subrou
         mes = messg
         do 30 i = 1,nmsg
            if (lib.eq.libtab(i) .and. sub.eq.subtab(i) .and.
     *         mes.eq.mestab(i) .and. nerr.eq.nertab(i) .and.
     *         level.eq.levtab(i)) then
                  kount(i) = kount(i) + 1
                  icount = kount(i)
                  return
            endif
   30    continue
c
         if (nmsg.lt.lentab) then
c
c           empty slot found for new message.
c
            nmsg = nmsg + 1
            libtab(i) = lib
            subtab(i) = sub
            mestab(i) = mes
            nertab(i) = nerr
            levtab(i) = level
            kount (i) = 1
            icount    = 1
         else
c
c           table is full.
c
            kountx = kountx+1
            icount = 0
         endif
      endif
      return
c
c     formats.
c
 9000 format ('0          error message summary' /
     +   ' library    subroutine message start             nerr',
     +   '     level     count')
 9010 format (1x,a,3x,a,3x,a,3i10)
 9020 format ('0other errors not individually tabulated = ', i10)
 9030 format (1x)
      end







c-------------------------------------------------------------------
*deck xgetua
      subroutine xgetua (iunita, n)
c***begin prologue  xgetua
c***purpose  return unit number(s) to which error messages are being
c            sent.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xgetua-a)
c***keywords  error, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        xgetua may be called to determine the unit number or numbers
c        to which error messages are being sent.
c        these unit numbers may have been set by a call to xsetun,
c        or a call to xsetua, or may be a default value.
c
c     description of parameters
c      --output--
c        iunit - an array of one to five unit numbers, depending
c                on the value of n.  a value of zero refers to the
c                default unit, as defined by the i1mach machine
c                constant routine.  only iunit(1),...,iunit(n) are
c                defined by xgetua.  the values of iunit(n+1),...,
c                iunit(5) are not defined (for n .lt. 5) or altered
c                in any way by xgetua.
c        n     - the number of units to which copies of the
c                error messages are being sent.  n will be in the
c                range from 1 to 5.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  j4save
c***revision history  (yymmdd)
c   790801  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xgetua
      dimension iunita(5)
c***first executable statement  xgetua
      n = j4save(5,0,.false.)
      do 30 i=1,n
         index = i+4
         if (i.eq.1) index = 3
         iunita(i) = j4save(index,0,.false.)
   30 continue
      return
      end







c-------------------------------------------------------------------
      integer function i1mach(i)
      integer i
c
c    i1mach( 1) = the standard input unit.
c    i1mach( 2) = the standard output unit.
c    i1mach( 3) = the standard punch unit.
c    i1mach( 4) = the standard error message unit.
c    i1mach( 5) = the number of bits per integer storage unit.
c    i1mach( 6) = the number of characters per character storage unit.
c    integers have form sign ( x(s-1)*a**(s-1) + ... + x(1)*a + x(0) )
c    i1mach( 7) = a, the base.
c    i1mach( 8) = s, the number of base-a digits.
c    i1mach( 9) = a**s - 1, the largest magnitude.
c    floats have form  sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )
c               where  emin .le. e .le. emax.
c    i1mach(10) = b, the base.
c  single-precision
c    i1mach(11) = t, the number of base-b digits.
c    i1mach(12) = emin, the smallest exponent e.
c    i1mach(13) = emax, the largest exponent e.
c  double-precision
c    i1mach(14) = t, the number of base-b digits.
c    i1mach(15) = emin, the smallest exponent e.
c    i1mach(16) = emax, the largest exponent e.
c
      integer imach(16), output, sc, small(2)
      save imach, sc
      real rmach
      equivalence (imach(4),output), (rmach,small(1))
      integer i3, j, k, t3e(3)
      data t3e(1) / 9777664 /
      data t3e(2) / 5323660 /
      data t3e(3) / 46980 /
c  this version adapts automatically to most current machines,
c  including auto-double compilers.
c  to compile on older machines, add a c in column 1
c  on the next line
      data sc/0/
c  and remove the c from column 1 in one of the sections below.
c  constants for even older machines can be obtained by
c          mail netlib@research.bell-labs.com
c          send old1mach from blas
c  please send corrections to dmg or ehg@bell-labs.com.
c
c     machine constants for the honeywell dps 8/70 series.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /   43 /
c      data imach( 4) /    6 /
c      data imach( 5) /   36 /
c      data imach( 6) /    4 /
c      data imach( 7) /    2 /
c      data imach( 8) /   35 /
c      data imach( 9) / o377777777777 /
c      data imach(10) /    2 /
c      data imach(11) /   27 /
c      data imach(12) / -127 /
c      data imach(13) /  127 /
c      data imach(14) /   63 /
c      data imach(15) / -127 /
c      data imach(16) /  127 /, sc/987/
c
c     machine constants for pdp-11 fortrans supporting
c     32-bit integer arithmetic.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   32 /
c      data imach( 6) /    4 /
c      data imach( 7) /    2 /
c      data imach( 8) /   31 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /    2 /
c      data imach(11) /   24 /
c      data imach(12) / -127 /
c      data imach(13) /  127 /
c      data imach(14) /   56 /
c      data imach(15) / -127 /
c      data imach(16) /  127 /, sc/987/
c
c     machine constants for the univac 1100 series.
c
c     note that the punch unit, i1mach(3), has been set to 7
c     which is appropriate for the univac-for system.
c     if you have the univac-ftn system, set it to 1.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   36 /
c      data imach( 6) /    6 /
c      data imach( 7) /    2 /
c      data imach( 8) /   35 /
c      data imach( 9) / o377777777777 /
c      data imach(10) /    2 /
c      data imach(11) /   27 /
c      data imach(12) / -128 /
c      data imach(13) /  127 /
c      data imach(14) /   60 /
c      data imach(15) /-1024 /
c      data imach(16) / 1023 /, sc/987/
c
      if (sc .ne. 987) then
*        *** check for autodouble ***
         small(2) = 0
         rmach = 1e13
         if (small(2) .ne. 0) then
*           *** autodoubled ***
            if (      (small(1) .eq. 1117925532
     *           .and. small(2) .eq. -448790528)
     *       .or.     (small(2) .eq. 1117925532
     *           .and. small(1) .eq. -448790528)) then
*               *** ieee ***
               imach(10) = 2
               imach(14) = 53
               imach(15) = -1021
               imach(16) = 1024
            else if ( small(1) .eq. -2065213935
     *          .and. small(2) .eq. 10752) then
*               *** vax with d_floating ***
               imach(10) = 2
               imach(14) = 56
               imach(15) = -127
               imach(16) = 127
            else if ( small(1) .eq. 1267827943
     *          .and. small(2) .eq. 704643072) then
*               *** ibm mainframe ***
               imach(10) = 16
               imach(14) = 14
               imach(15) = -64
               imach(16) = 63
            else
               write(*,9010)
               stop 777
               end if
            imach(11) = imach(14)
            imach(12) = imach(15)
            imach(13) = imach(16)
         else
            rmach = 1234567.
            if (small(1) .eq. 1234613304) then
*               *** ieee ***
               imach(10) = 2
               imach(11) = 24
               imach(12) = -125
               imach(13) = 128
               imach(14) = 53
               imach(15) = -1021
               imach(16) = 1024
               sc = 987
            else if (small(1) .eq. -1271379306) then
*               *** vax ***
               imach(10) = 2
               imach(11) = 24
               imach(12) = -127
               imach(13) = 127
               imach(14) = 56
               imach(15) = -127
               imach(16) = 127
               sc = 987
            else if (small(1) .eq. 1175639687) then
*               *** ibm mainframe ***
               imach(10) = 16
               imach(11) = 6
               imach(12) = -64
               imach(13) = 63
               imach(14) = 14
               imach(15) = -64
               imach(16) = 63
               sc = 987
            else if (small(1) .eq. 1251390520) then
*              *** convex c-1 ***
               imach(10) = 2
               imach(11) = 24
               imach(12) = -128
               imach(13) = 127
               imach(14) = 53
               imach(15) = -1024
               imach(16) = 1023
            else
               do 10 i3 = 1, 3
                  j = small(1) / 10000000
                  k = small(1) - 10000000*j
                  if (k .ne. t3e(i3)) go to 20
                  small(1) = j
 10               continue
*              *** cray t3e ***
               imach( 1) = 5
               imach( 2) = 6
               imach( 3) = 0
               imach( 4) = 0
               imach( 5) = 64
               imach( 6) = 8
               imach( 7) = 2
               imach( 8) = 63
               call i1mcr1(imach(9), k, 32767, 16777215, 16777215)
               imach(10) = 2
               imach(11) = 53
               imach(12) = -1021
               imach(13) = 1024
               imach(14) = 53
               imach(15) = -1021
               imach(16) = 1024
               go to 35
 20            call i1mcr1(j, k, 16405, 9876536, 0)
               if (small(1) .ne. j) then
                  write(*,9020)
                  stop 777
                  end if
*              *** cray 1, xmp, 2, and 3 ***
               imach(1) = 5
               imach(2) = 6
               imach(3) = 102
               imach(4) = 6
               imach(5) = 46
               imach(6) = 8
               imach(7) = 2
               imach(8) = 45
               call i1mcr1(imach(9), k, 0, 4194303, 16777215)
               imach(10) = 2
               imach(11) = 47
               imach(12) = -8188
               imach(13) = 8189
               imach(14) = 94
               imach(15) = -8141
               imach(16) = 8189
               go to 35
               end if
            end if
         imach( 1) = 5
         imach( 2) = 6
         imach( 3) = 7
         imach( 4) = 6
         imach( 5) = 32
         imach( 6) = 4
         imach( 7) = 2
         imach( 8) = 31
         imach( 9) = 2147483647
 35      sc = 987
         end if
 9010 format(/' adjust autodoubled i1mach by uncommenting data'/
     * ' statements appropriate for your machine and setting'/
     * ' imach(i) = imach(i+3) for i = 11, 12, and 13.')
 9020 format(/' adjust i1mach by uncommenting data statements'/
     * ' appropriate for your machine.')
      if (i .lt. 1  .or.  i .gt. 16) go to 40
      i1mach = imach(i)
      return
 40   write(*,*) 'i1mach(i): i =',i,' is out of bounds.'
      stop
* /* c source for i1mach -- remove the * in column 1 */
* /* note that some values may need changing. */
*#include <stdio.h>
*#include <float.h>
*#include <limits.h>
*#include <math.h>
*
*long i1mach_(long *i)
*{
*	switch(*i){
*	  case 1:  return 5;	/* standard input */
*	  case 2:  return 6;	/* standard output */
*	  case 3:  return 7;	/* standard punch */
*	  case 4:  return 0;	/* standard error */
*	  case 5:  return 32;	/* bits per integer */
*	  case 6:  return sizeof(int);
*	  case 7:  return 2;	/* base for integers */
*	  case 8:  return 31;	/* digits of integer base */
*	  case 9:  return long_max;
*	  case 10: return flt_radix;
*	  case 11: return flt_mant_dig;
*	  case 12: return flt_min_exp;
*	  case 13: return flt_max_exp;
*	  case 14: return dbl_mant_dig;
*	  case 15: return dbl_min_exp;
*	  case 16: return dbl_max_exp;
*	  }
*	fprintf(stderr, "invalid argument: i1mach(%ld)\n", *i);
*	exit(1);return 0; /* some compilers demand return values */
*}
      end
      subroutine i1mcr1(a, a1, b, c, d)
**** special computation for old cray machines ****
      integer a, a1, b, c, d
      a1 = 16777216*b + c
      a = 16777216*a1 + d
      end







c-------------------------------------------------------------------
      real function r1mach(i)
      integer i
c
c  single-precision machine constants
c  r1mach(1) = b**(emin-1), the smallest positive magnitude.
c  r1mach(2) = b**emax*(1 - b**(-t)), the largest magnitude.
c  r1mach(3) = b**(-t), the smallest relative spacing.
c  r1mach(4) = b**(1-t), the largest relative spacing.
c  r1mach(5) = log10(b)
c
      integer small(2)
      integer large(2)
      integer right(2)
      integer diver(2)
      integer log10(2)
c     needs to be (2) for autodouble, harris slash 6, ...
      integer sc
      save small, large, right, diver, log10, sc
      real rmach(5)
      equivalence (rmach(1),small(1))
      equivalence (rmach(2),large(1))
      equivalence (rmach(3),right(1))
      equivalence (rmach(4),diver(1))
      equivalence (rmach(5),log10(1))
      integer j, k, l, t3e(3)
      data t3e(1) / 9777664 /
      data t3e(2) / 5323660 /
      data t3e(3) / 46980 /
c  this version adapts automatically to most current machines,
c  including auto-double compilers.
c  to compile on older machines, add a c in column 1
c  on the next line
      data sc/0/
c  and remove the c from column 1 in one of the sections below.
c  constants for even older machines can be obtained by
c          mail netlib@research.bell-labs.com
c          send old1mach from blas
c  please send corrections to dmg or ehg@bell-labs.com.
c
c     machine constants for the honeywell dps 8/70 series.
c      data rmach(1) / o402400000000 /
c      data rmach(2) / o376777777777 /
c      data rmach(3) / o714400000000 /
c      data rmach(4) / o716400000000 /
c      data rmach(5) / o776464202324 /, sc/987/
c
c     machine constants for pdp-11 fortrans supporting
c     32-bit integers (expressed in integer and octal).
c      data small(1) /    8388608 /
c      data large(1) / 2147483647 /
c      data right(1) /  880803840 /
c      data diver(1) /  889192448 /
c      data log10(1) / 1067065499 /, sc/987/
c      data rmach(1) / o00040000000 /
c      data rmach(2) / o17777777777 /
c      data rmach(3) / o06440000000 /
c      data rmach(4) / o06500000000 /
c      data rmach(5) / o07746420233 /, sc/987/
c
c     machine constants for the univac 1100 series.
c      data rmach(1) / o000400000000 /
c      data rmach(2) / o377777777777 /
c      data rmach(3) / o146400000000 /
c      data rmach(4) / o147400000000 /
c      data rmach(5) / o177464202324 /, sc/987/
c
      if (sc .ne. 987) then
*        *** check for autodouble ***
         small(2) = 0
         rmach(1) = 1e13
         if (small(2) .ne. 0) then
*           *** autodoubled ***
            if (      small(1) .eq. 1117925532
     *          .and. small(2) .eq. -448790528) then
*              *** ieee big endian ***
               small(1) = 1048576
               small(2) = 0
               large(1) = 2146435071
               large(2) = -1
               right(1) = 1017118720
               right(2) = 0
               diver(1) = 1018167296
               diver(2) = 0
               log10(1) = 1070810131
               log10(2) = 1352628735
            else if ( small(2) .eq. 1117925532
     *          .and. small(1) .eq. -448790528) then
*              *** ieee little endian ***
               small(2) = 1048576
               small(1) = 0
               large(2) = 2146435071
               large(1) = -1
               right(2) = 1017118720
               right(1) = 0
               diver(2) = 1018167296
               diver(1) = 0
               log10(2) = 1070810131
               log10(1) = 1352628735
            else if ( small(1) .eq. -2065213935
     *          .and. small(2) .eq. 10752) then
*              *** vax with d_floating ***
               small(1) = 128
               small(2) = 0
               large(1) = -32769
               large(2) = -1
               right(1) = 9344
               right(2) = 0
               diver(1) = 9472
               diver(2) = 0
               log10(1) = 546979738
               log10(2) = -805796613
            else if ( small(1) .eq. 1267827943
     *          .and. small(2) .eq. 704643072) then
*              *** ibm mainframe ***
               small(1) = 1048576
               small(2) = 0
               large(1) = 2147483647
               large(2) = -1
               right(1) = 856686592
               right(2) = 0
               diver(1) = 873463808
               diver(2) = 0
               log10(1) = 1091781651
               log10(2) = 1352628735
            else
               write(*,9010)
               stop 777
               end if
         else
            rmach(1) = 1234567.
            if (small(1) .eq. 1234613304) then
*              *** ieee ***
               small(1) = 8388608
               large(1) = 2139095039
               right(1) = 864026624
               diver(1) = 872415232
               log10(1) = 1050288283
            else if (small(1) .eq. -1271379306) then
*              *** vax ***
               small(1) = 128
               large(1) = -32769
               right(1) = 13440
               diver(1) = 13568
               log10(1) = 547045274
            else if (small(1) .eq. 1175639687) then
*              *** ibm mainframe ***
               small(1) = 1048576
               large(1) = 2147483647
               right(1) = 990904320
               diver(1) = 1007681536
               log10(1) = 1091781651
            else if (small(1) .eq. 1251390520) then
*              *** convex c-1 ***
               small(1) = 8388608
               large(1) = 2147483647
               right(1) = 880803840
               diver(1) = 889192448
               log10(1) = 1067065499
            else
               do 10 l = 1, 3
                  j = small(1) / 10000000
                  k = small(1) - 10000000*j
                  if (k .ne. t3e(l)) go to 20
                  small(1) = j
 10               continue
*              *** cray t3e ***
               call i1mcra(small, k, 16, 0, 0)
               call i1mcra(large, k, 32751, 16777215, 16777215)
               call i1mcra(right, k, 15520, 0, 0)
               call i1mcra(diver, k, 15536, 0, 0)
               call i1mcra(log10, k, 16339, 4461392, 10451455)
               go to 30
 20            call i1mcra(j, k, 16405, 9876536, 0)
               if (small(1) .ne. j) then
                  write(*,9020)
                  stop 777
                  end if
*              *** cray 1, xmp, 2, and 3 ***
               call i1mcra(small(1), k, 8195, 8388608, 1)
               call i1mcra(large(1), k, 24574, 16777215, 16777214)
               call i1mcra(right(1), k, 16338, 8388608, 0)
               call i1mcra(diver(1), k, 16339, 8388608, 0)
               call i1mcra(log10(1), k, 16383, 10100890, 8715216)
               end if
            end if
 30      sc = 987
         end if
*     sanity check
      if (rmach(4) .ge. 1.0) stop 776
      if (i .lt. 1 .or. i .gt. 5) then
         write(*,*) 'r1mach(i): i =',i,' is out of bounds.'
         stop
         end if
      r1mach = rmach(i)
      return
 9010 format(/' adjust autodoubled r1mach by getting data'/
     *' appropriate for your machine from d1mach.')
 9020 format(/' adjust r1mach by uncommenting data statements'/
     *' appropriate for your machine.')
* /* c source for r1mach -- remove the * in column 1 */
*#include <stdio.h>
*#include <float.h>
*#include <math.h>
*float r1mach_(long *i)
*{
*	switch(*i){
*	  case 1: return flt_min;
*	  case 2: return flt_max;
*	  case 3: return flt_epsilon/flt_radix;
*	  case 4: return flt_epsilon;
*	  case 5: return log10((double)flt_radix);
*	  }
*	fprintf(stderr, "invalid argument: r1mach(%ld)\n", *i);
*	exit(1); return 0; /* else complaint of missing return value */
*}
      end
      subroutine i1mcra(a, a1, b, c, d)
**** special computation for cray machines ****
      integer a, a1, b, c, d
      a1 = 16777216*b + c
      a = 16777216*a1 + d
      end







c-------------------------------------------------------------------
      double precision function d1mach(i)
      integer i
c
c  double-precision machine constants
c  d1mach( 1) = b**(emin-1), the smallest positive magnitude.
c  d1mach( 2) = b**emax*(1 - b**(-t)), the largest magnitude.
c  d1mach( 3) = b**(-t), the smallest relative spacing.
c  d1mach( 4) = b**(1-t), the largest relative spacing.
c  d1mach( 5) = log10(b)
c
      integer small(2)
      integer large(2)
      integer right(2)
      integer diver(2)
      integer log10(2)
      integer sc, cray1(38), j
      common /d9mach/ cray1
      save small, large, right, diver, log10, sc
      double precision dmach(5)
      equivalence (dmach(1),small(1))
      equivalence (dmach(2),large(1))
      equivalence (dmach(3),right(1))
      equivalence (dmach(4),diver(1))
      equivalence (dmach(5),log10(1))
c  this version adapts automatically to most current machines.
c  r1mach can handle auto-double compiling, but this version of
c  d1mach does not, because we do not have quad constants for
c  many machines yet.
c  to compile on older machines, add a c in column 1
c  on the next line
      data sc/0/
c  and remove the c from column 1 in one of the sections below.
c  constants for even older machines can be obtained by
c          mail netlib@research.bell-labs.com
c          send old1mach from blas
c  please send corrections to dmg or ehg@bell-labs.com.
c
c     machine constants for the honeywell dps 8/70 series.
c      data small(1),small(2) / o402400000000, o000000000000 /
c      data large(1),large(2) / o376777777777, o777777777777 /
c      data right(1),right(2) / o604400000000, o000000000000 /
c      data diver(1),diver(2) / o606400000000, o000000000000 /
c      data log10(1),log10(2) / o776464202324, o117571775714 /, sc/987/
c
c     machine constants for pdp-11 fortrans supporting
c     32-bit integers.
c      data small(1),small(2) /    8388608,           0 /
c      data large(1),large(2) / 2147483647,          -1 /
c      data right(1),right(2) /  612368384,           0 /
c      data diver(1),diver(2) /  620756992,           0 /
c      data log10(1),log10(2) / 1067065498, -2063872008 /, sc/987/
c
c     machine constants for the univac 1100 series.
c      data small(1),small(2) / o000040000000, o000000000000 /
c      data large(1),large(2) / o377777777777, o777777777777 /
c      data right(1),right(2) / o170540000000, o000000000000 /
c      data diver(1),diver(2) / o170640000000, o000000000000 /
c      data log10(1),log10(2) / o177746420232, o411757177572 /, sc/987/
c
c     on first call, if no data uncommented, test machine types.
      if (sc .ne. 987) then
         dmach(1) = 1.d13
         if (      small(1) .eq. 1117925532
     *       .and. small(2) .eq. -448790528) then
*           *** ieee big endian ***
            small(1) = 1048576
            small(2) = 0
            large(1) = 2146435071
            large(2) = -1
            right(1) = 1017118720
            right(2) = 0
            diver(1) = 1018167296
            diver(2) = 0
            log10(1) = 1070810131
            log10(2) = 1352628735
         else if ( small(2) .eq. 1117925532
     *       .and. small(1) .eq. -448790528) then
*           *** ieee little endian ***
            small(2) = 1048576
            small(1) = 0
            large(2) = 2146435071
            large(1) = -1
            right(2) = 1017118720
            right(1) = 0
            diver(2) = 1018167296
            diver(1) = 0
            log10(2) = 1070810131
            log10(1) = 1352628735
         else if ( small(1) .eq. -2065213935
     *       .and. small(2) .eq. 10752) then
*               *** vax with d_floating ***
            small(1) = 128
            small(2) = 0
            large(1) = -32769
            large(2) = -1
            right(1) = 9344
            right(2) = 0
            diver(1) = 9472
            diver(2) = 0
            log10(1) = 546979738
            log10(2) = -805796613
         else if ( small(1) .eq. 1267827943
     *       .and. small(2) .eq. 704643072) then
*               *** ibm mainframe ***
            small(1) = 1048576
            small(2) = 0
            large(1) = 2147483647
            large(2) = -1
            right(1) = 856686592
            right(2) = 0
            diver(1) = 873463808
            diver(2) = 0
            log10(1) = 1091781651
            log10(2) = 1352628735
         else if ( small(1) .eq. 1120022684
     *       .and. small(2) .eq. -448790528) then
*           *** convex c-1 ***
            small(1) = 1048576
            small(2) = 0
            large(1) = 2147483647
            large(2) = -1
            right(1) = 1019215872
            right(2) = 0
            diver(1) = 1020264448
            diver(2) = 0
            log10(1) = 1072907283
            log10(2) = 1352628735
         else if ( small(1) .eq. 815547074
     *       .and. small(2) .eq. 58688) then
*           *** vax g-floating ***
            small(1) = 16
            small(2) = 0
            large(1) = -32769
            large(2) = -1
            right(1) = 15552
            right(2) = 0
            diver(1) = 15568
            diver(2) = 0
            log10(1) = 1142112243
            log10(2) = 2046775455
         else
            dmach(2) = 1.d27 + 1
            dmach(3) = 1.d27
            large(2) = large(2) - right(2)
            if (large(2) .eq. 64 .and. small(2) .eq. 0) then
               cray1(1) = 67291416
               do 10 j = 1, 20
                  cray1(j+1) = cray1(j) + cray1(j)
 10               continue
               cray1(22) = cray1(21) + 321322
               do 20 j = 22, 37
                  cray1(j+1) = cray1(j) + cray1(j)
 20               continue
               if (cray1(38) .eq. small(1)) then
*                  *** cray ***
                  call i1mcry(small(1), j, 8285, 8388608, 0)
                  small(2) = 0
                  call i1mcry(large(1), j, 24574, 16777215, 16777215)
                  call i1mcry(large(2), j, 0, 16777215, 16777214)
                  call i1mcry(right(1), j, 16291, 8388608, 0)
                  right(2) = 0
                  call i1mcry(diver(1), j, 16292, 8388608, 0)
                  diver(2) = 0
                  call i1mcry(log10(1), j, 16383, 10100890, 8715215)
                  call i1mcry(log10(2), j, 0, 16226447, 9001388)
               else
                  write(*,9000)
                  stop 779
                  end if
            else
               write(*,9000)
               stop 779
               end if
            end if
         sc = 987
         end if
*    sanity check
      if (dmach(4) .ge. 1.0d0) stop 778
      if (i .lt. 1 .or. i .gt. 5) then
         write(*,*) 'd1mach(i): i =',i,' is out of bounds.'
         stop
         end if
      d1mach = dmach(i)
      return
 9000 format(/' adjust d1mach by uncommenting data statements'/
     *' appropriate for your machine.')
* /* standard c source for d1mach -- remove the * in column 1 */
*#include <stdio.h>
*#include <float.h>
*#include <math.h>
*double d1mach_(long *i)
*{
*	switch(*i){
*	  case 1: return dbl_min;
*	  case 2: return dbl_max;
*	  case 3: return dbl_epsilon/flt_radix;
*	  case 4: return dbl_epsilon;
*	  case 5: return log10((double)flt_radix);
*	  }
*	fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
*	exit(1); return 0; /* some compilers demand return values */
*}
      end
      subroutine i1mcry(a, a1, b, c, d)
**** special computation for old cray machines ****
      integer a, a1, b, c, d
      a1 = 16777216*b + c
      a = 16777216*a1 + d
      end







