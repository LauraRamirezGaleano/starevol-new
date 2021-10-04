

************************************************************************

      SUBROUTINE centeq

************************************************************************
* Give the central boundary conditions for stellar structure           *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.eng'
      include 'evolcom.eos'
      include 'evolcom.grad'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.opa'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.var'

      logical vlconv

      integer k,l

      double precision dtninv,dmkb,dlnt1,dt4,sph,gmrua1,plo
      double precision ablainv
      double precision lconv,lrad
      double precision ft,fp

      common /resolution/ vlconv
      common /ftfp/ ft(nsh),fp(nsh)

      dtninv = 1.d0/dtn

      sph = 3.d0/(pim4*ro(1))
      dmkb = (dm(2)+dm(1))*0.5d0
      dlnt1 = log(t(2)/t(1))

      forall (k=1:neq,l=1:neq1) eq(k,l) = 0.d0

***   central velocity
      if (hydro) then
         eq(1,16) = u(1)
      endif
      eq(1,1) = 1.d0

***   central radius
c      eq(2,16) = r(1)
      eq(2,2) = r(1)

***   mass conservation
      eq(3,16) = r(2)**3/dm(1)-sph
      eq(3,3) = sph*drodf(1)/ro(1)
      eq(3,4) = sph*drodt(1)/ro(1)
      eq(3,7) = 3.d0*r(2)**3/dm(1)

***   central luminosity
      eq(4,16) = lum(1)
      eq(4,5) = 1.d0

***   equation of transport
      dt4 = t(2)**4-t(1)**4
      if (abs(dt4).lt.1.d-5) dt4 = dt4+1.d-5
      if (vlconv) then
         lrad = -64.d0*sig*pi**2*r(2)**4*dt4/(3.d0*kapm(2)*dmkb)
         lconv = vhconv(2)*rom(2)*tm(2)*r(2)

         eq(5,16) = lum(2)-lrad-lconv
         eq(5,3) = wi(2)*(lrad*dkapdf(1)/kap(1)-lconv*drodf(1)/ro(1))
         eq(5,4) = lrad*(wi(2)*dkapdt(1)/kap(1)+4.d0*t(1)**4/dt4)-
     &        lconv*wi(2)*(1.d0+drodt(1)/ro(1))
         eq(5,7) =  -(4.d0*lrad+lconv)
         eq(5,8) = wj(2)*(lrad*dkapdf(2)/kap(2)-lconv*drodf(2)/ro(2))
         eq(5,9) = lrad*(wj(2)*dkapdt(2)/kap(2)-4.d0*t(2)**4/dt4)-
     &        lconv*wj(2)*(1.d0+drodt(2)/ro(2))
         eq(5,10) = 1.d0
      else
         if (crz(2).lt.-1) then
            gmrua1 = gmr(2)*fp(2)+accel(2)
            plo = gmrua1*abla(2)*dmkb/(pim4*r(2)*r(2)*pm(2))
            ablainv = 1.d0/abla(2)
            eq(5,16) = dlnt1+plo
            eq(5,1) = plo*abdu1(2)*ablainv
            eq(5,3) = plo*(abdf1(2)*ablainv-wi(2)*dpdf(1)/p(1))
            eq(5,4) = plo*(abdt1(2)*ablainv-wi(2)*dpdt(1)/p(1))-1.d0
            eq(5,6) = plo*(abdu2(2)*ablainv+dynfac*dtninv/gmrua1)
            eq(5,7) = plo*(abdr(2)*ablainv*r(2)-2.d0*(1.d0+gmr(2)/
     &           gmrua1))
            eq(5,8) = plo*(abdf2(2)*ablainv-wj(2)*dpdf(2)/p(2))
            eq(5,9) = plo*(abdt2(2)*ablainv-wj(2)*dpdt(2)/p(2))+1.d0
            eq(5,10) = plo*abdl(2)*ablainv
         else
            plo = 3.d0*dmkb/(256.d0*sig*pi**2)*lum(2)*kapm(2)/
     &           (r(2)*tm(2))**4*ft(2)
            eq(5,16) = dlnt1+plo
            eq(5,3) = plo*wi(2)*dkapdf(1)/kap(1)
            eq(5,4) = plo*wi(2)*(dkapdt(1)/kap(1)-4.d0)-1.d0
            eq(5,7) = -4.d0*plo
            eq(5,8) = plo*wj(2)*dkapdf(2)/kap(2)
            eq(5,9) = plo*wj(2)*(dkapdt(2)/kap(2)-4.d0)+1.d0
            eq(5,10) = plo/lum(2)
         endif
      endif

      return
      end
