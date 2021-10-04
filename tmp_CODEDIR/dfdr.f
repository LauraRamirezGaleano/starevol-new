

************************************************************************

      DOUBLE PRECISION FUNCTION dfdr (i,h,f,schwarr,idet)

************************************************************************
* Calculate the derivative of the penetration function                 *
* for the accretion process                                            *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.grad'
      include 'evolcom.ion'
      include 'evolcom.rot'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer i,idet,isup

      double precision h,f,schwarr
      double precision rr,muacc,ri,vf,b0,gradmu,b1,c,det,fthacc

      isup = nmod
      dfdr = 1.d-21
      if (idet.eq.0) return
      rr = r(i)+h
      fthacc = f*thacc
      if (iaccr.eq.1.or.iaccr.eq.3) muacc = mueinv(i)+fthacc*muiacc+
     &     (1.d0-fthacc)*muiinv(i)
      if (iaccr.eq.2.or.iaccr.eq.4) muacc = mueinv(i)+(fthacc*
     &     muiacc+muiinv(i))/(1.d0+fthacc)
      if (abrad(i).gt.abad(i)) then
         schwarr = -abs(schwarr)
         ri = ric
      else
         ri = rir
      endif

      vf = vomega(i)*r(i)*r(i)/sigma0
      if (iaccr.eq.1.or.iaccr.eq.3) then
         b0 = -thacc*(muiacc-muiinv(i))*phiKS(i)/muacc
         gradmu = abmu(i)*(1.d0-thacc*f)
      endif
      if (iaccr.eq.2.or.iaccr.eq.4) then
         b0 = -thacc*(muiacc-muiinv(i))/(1.d0+thacc*f)**2*phiKS(i)/
     &        muacc
         gradmu = abmu(i)/(1.d0+fthacc)
      endif
      if (ri.eq.ric) then
         gradmu = 0.d0
         abmu(i) = 0.d0
      endif
      b1 = m(i)/(thacc*m(isup)*xiaccr*xiaccr)-(f-vf)**2*r(isup)/rr
      b0 = b1*b0-2.d0*(f-vf)*r(isup)/rr
      c = (-1.d0/(mu(i)*muacc)*gradmu*phiKS(i)-prrc*deltaKS(i)*
     &     schwarr)/hp(i)*b1
      det = b0*b0-4.d0*r(isup)*ri*c
      if (det.lt.0.d0) then
         idet = 0
         return
      endif

      dfdr = abs(b0/(2.d0*r(isup)*ri)-abs(dsqrt(det)/(2.d0*r(isup)*ri)))
      if (dfdr.lt.1.d-20) dfdr = 1.d-21

      return
      end
