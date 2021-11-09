

************************************************************************

      DOUBLE PRECISION FUNCTION PCONV (a,V,B)

************************************************************************
* Search for the convective solution in case V <= 0.1                  *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      double precision a,V,B
      double precision a2,a3,V2,p,q,p1,f,g,dp,tol

      V2 = V*V
      a2 = a*a
      a3 = a*a2
      p = 1.d0/a
      p1 = 1.d0/a3
      tol = max((p-p1)*1.d-8,1.d-12)
      pconv = p
 10   q = (p+p1)*0.5d0
      f =  B*(1.d0-p*a3)*(1.d0-p*a3)**2*p*(1.d0+p**2*a)+
     &     (1+p)*(p**2*a2-1.d0)*((1.d0-p*a**3)*a2*p*(p-a)-
     &     (p**2*a2-1)**2/V2)
      g =  B*(1.d0-q*a3)*(1.d0-q*a3)**2*q*(1.d0+q**2*a)+
     &     (1+q)*(q**2*a2-1.d0)*((1.d0-q*a**3)*a2*q*(q-a)-
     &     (q**2*a2-1)**2/V2)
      if ((f.gt.0.d0.and.g.lt.0.d0).or.(f.lt.0.d0.and.g.gt.0.d0))
     &     then
         p1 = q
      else
         p = q
      endif
      pconv = p
      dp = abs(p-p1)/p
      if (dp.lt.tol) goto 30
      goto 10

 30   return
      end
