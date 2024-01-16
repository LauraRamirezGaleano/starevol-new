
      SUBROUTINE spline_interp(xa,ya,y2a,n,x,y)

************************************************************************
*                                                                      *
* Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a   *
* function (with the xai’s in order), and given the array y2a(1:n),    *
* which is the output from spline above, and given a value of x, this  *
* routine returns a cubic-spline interpolated value y.                 *
* From Press, W.H. et al. (1992). Numerical recipes in FORTRAN. The    *
* art of scientific computing. ISBN : 0-521-43064-X                    * 
*                                                                      * 
************************************************************************
      
      implicit none
      
      integer n
      double precision x,y,xa(n),y2a(n),ya(n)
      integer k,khi,klo
      double precision a,b,h
      
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x) then
            khi=k
         else
            klo=k
         endif
         goto 1

      endif
      h=xa(khi)-xa(klo)
      a=(xa(khi)-x)/h          
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     &     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      
      return
      END
