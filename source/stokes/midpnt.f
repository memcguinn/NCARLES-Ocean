      subroutine midpnt(a,b,s,n)
c
      integer it,j

c
c       write(6,1200) a,b,s,n
 1200 format(' 1200 a = ',e15.6,' b = ',e15.6,/,
     +       '      s = ',e15.6,' n = ',i5)
C
      if(n .eq. 1) then
         s = (b - a)*stokes_ker(0.5*(a+b))
      else
         it   = 3**(n-2)
c         write(6,6000) n, it
 6000    format(' n = ',i4,' it = ',i4)
         tnm  = float(it)
         del  = (b - a)/(3.0*tnm)
         ddel = del + del
         x    = a + 0.5*del
         sum = 0.0
         do j=1,it
            sum = sum + stokes_ker(x)
            x   = x + ddel
            sum = sum + stokes_ker(x)
            x   = x + del
         enddo
         s = (s + (b - a)*sum/tnm)/3.0
      endif
      return
      end
