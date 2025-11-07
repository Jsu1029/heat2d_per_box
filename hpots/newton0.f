      subroutine newrap(f,f1,p0,delta,epsilon,max,p1,dp,y1,cond,k)
      parameter(small=1e-20)
      integer cond,k,max
      real delta,epsilon,df,dp,p0,p1,y0,y1,relerr
      external f,f1
      k=0
      cond=0
      y0=f(p0)
      p1=p0+1
      while ((k.le.max).and.(cond.eq.0))
        df=f1(p0)
        if (df.eq.0) then
          cond=1
          dp=p1-p0
          p1=p0
        else
          dp=y0/df
          p1=p0-dp
        endif
        y1=f(p1)
        relerr=abs(dp)/(abs(p1)+small)
        if (relerr.lt.delta) cond=2
        if (abs(y1).lt.epsilon) cond=3
        if (cond.eq.2).and.(abs(y1).lt.epsilon) cond=4
        p0=p1
        y0=y1
        k=k+1
        write(9,1000) k,p1,y1
      repeat
c      pause
      return
1000  format(i2,4x,f15.7,4x,f15.7)
      end subroutine

