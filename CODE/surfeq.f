      SUBROUTINE surfeq

************************************************************************
* Give the outer boundary conditions for stellar structure             *
* When ntprof = 2, the outer boundary conditions for T and rho depend  *
* only on T(Eddington) and L(Eddington) and not on tau (tau0 = 2/3 in  *
* this case)                                                           *
*                                                                      *
* $LastChangedDate:: 2016-05-13 11:39:14 +0200 (Ven, 13 mai 2016)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 75                                                          $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.eng'
      include 'evolcom.eos'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.opa'
      include 'evolcom.rot'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer k,l

      logical vlconv,ntp2

      double precision dtninv,dtiomi,dtipsi
      double precision dp,gmrua,dtau,tsf,vfac,rosurf1,rosurf2,rosurf,
     &     vrm,vgmr,vdp,ddtau,sigvrm,sigma1
      double precision fac1,rrk,lrad,drvfacc
      double precision tsurfr,kiinv,kiminv
      double precision ledd
c      double precision fac2,dvfacc,efacc
      double precision spgp,sfequi

      double precision rae
      double precision disctime
!      double precision rhoatm(nmod),Patm(nmod),convatm(nmod)
      double precision rhoat
      double precision Tatm

      common /resolution/ vlconv
      common /spgp/ spgp(nsh),sfequi(nsh)
      common /radi/ rae(intmax)
      common /disclocking/ disctime
      common /interpatm/ rhoat(nsh),Tatm(nsh),ntp2


      dtninv = 1.d0/dtn
      dp = p(nmod)-p(nmod1)
      vdp = vp(nmod)-vp(nmod1)
      gmrua = gmr(nmod)+accel(nmod)
      if (ntprof.le.1.or.ntprof.ge.3) then
         dtau = (qtau0-tau0-qtaus)
         if (hydrorot) then
*** Eq. A.40 Appendix in Meynet & Maeder 1997, A&A 321,465
	    tsurfr = teff
            tsf = (pw34*(sfequi(nmod)/(pim4*r(nmod)**2)*ft(nmod)*
     &           tau0+qtaus*gmrua
     $           *sfequi(nmod)/(pim4*spgp(nmod))))**0.25d0
         else
            tsurfr = teff
            tsf = (pw34*(tau0+qtaus))**0.25d0
         endif
         vfac = 1.d0/(rk*tsurfr*tsf*muinv(nmod))
         rosurf1 = vfac*tau0*gmrua/kap(nmod)
         rosurf2 = vfac*sig/c*(tsurfr*tsf)**4*dtau
         rosurf = rosurf1+rosurf2
         if ((ntprof.eq.4.and..not.ntp2).or.ntprof.eq.6.or.
     &        ntprof.eq.7) then
            rosurf = rhoat(nmod)
         endif

***   Surface boundary condition defined with Eddington values
***   No atmosphere treated (in this case, tau0 = 2/3) (see Langer)
      else if (ntprof.eq.2) then
         tsurfr = (abs(lum(nmod))/(pim4*sig*r(nmod)**2))**0.25d0
         vfac = 1.d0/(rk*tsurfr*muinv(nmod))
         ledd = pim4*c*g*m(nmod)/kap(nmod)
         rosurf2 = -pw23*vfac*sig/c*tsurfr**4
         rosurf = vfac*pw23*(ledd-lum(nmod))*gmrua/(kap(nmod)*ledd)
c     &        +rosurf2
      endif
      if (ntprof.eq.1.or.(ntprof.ge.3.and..not.ntp2).or.ntprof.eq.6
     &     .or.ntprof.eq.7) then
         ddtau = 0.d0
         if (abs(dqdtau(nmod)).gt.1.d-6) ddtau = ddqtau(nmod)/
     &        dqdtau(nmod)
      endif
      vrm = pim8*r(nmod)*r(nmod)/dm(nmod1)
      vgmr = g*m(nmod)/(vr(nmod)*vr(nmod))
      sigma1 = 1.d0-sigma
      sigvrm = sigma*vrm
      dtiomi = dtninv*omi(nmod)
      dtipsi = dtninv*psi(nmod)
c      dtiomi = 0.d0
c      dtipsi = 0.d0

      forall (k=1:neq,l=1:neq1) eq(k,l) = 0.d0


***   equation of motion
************************************************************************

      if (hydrorot) then
         eq(1,16) = accel(nmod)+sigma*(vrm*dp+gmr(nmod)*fp(nmod))
         eq(1,7) = 2.d0*sigma*(vrm*dp-gmr(nmod)*fp(nmod))
      else
         eq(1,16) = accel(nmod)+sigma*(vrm*dp+gmr(nmod)*fp(nmod))
     &        +sigma1*(vrm*vdp+vgmr)
         eq(1,7) = 2.d0*sigma*(vrm*dp-gmr(nmod)*fp(nmod))
      endif
      eq(1,1) = dynfac*dtipsi
      eq(1,3) = -sigvrm*dpdf(nmod1)
      eq(1,4) = -sigvrm*dpdt(nmod1)
      eq(1,6) = dynfac*(dtninv-dtipsi)
      eq(1,8) = sigvrm*dpdf(nmod)
      eq(1,9) = sigvrm*dpdt(nmod)
      if (ntprof.eq.1.or.(ntprof.ge.3.and..not.ntp2).or.ntprof.eq.6
     &     .or.ntprof.eq.7) then
         if (numeric.eq.2.or.numeric.eq.4) then
            kiinv = kapm(nmod)/kap(nmod)**2
            kiminv = kapm(nmod)/kap(nmod1)**2
         else
            kiinv = 1.d0/kap(nmod)
            kiminv = 1.d0/kap(nmod1)
         endif
         eq(1,16) = eq(1,16)+dptau(nmod)
         eq(1,3) = eq(1,3)+dptau(nmod)*dkapdf(nmod1)*kiminv
         eq(1,4) = eq(1,4)+dptau(nmod)*dkapdt(nmod1)*kiminv
         eq(1,7) = eq(1,7)
     &        -2.d0*sigma*gmr(nmod)*xsitau(nmod)
     &        +gmr(nmod)*dxsitaudr(nmod)
     &        !+pim4*r(nmod)**3*(gmr(nmod)+accel(nmod))*
     &        !(1+xsitau(nmod))*ro(nmod)
         eq(1,8) = eq(1,8)+gmr(nmod)*dxsitaudf(nmod)
         eq(1,9) = eq(1,9)+gmr(nmod)*dxsitaudt(nmod)
         eq(1,10) = eq(1,10)+gmr(nmod)*xsitau(nmod)/lum(nmod)
      endif

***   lagrangian velocity
************************************************************************

      if (hydro) then
         eq(2,16) = (lnr(nmod)-vlnr(nmod)-psi(nmod)*(lnr(nmod)-
     &        lnr(nmod1)))*dtninv-sigma*u(nmod)/r(nmod)-sigma1*
     &        vu(nmod)/vr(nmod)
         eq(2,2) = dtipsi
         eq(2,6) = -sigma/r(nmod)
         eq(2,7) = dtninv-dtipsi+sigma*u(nmod)/r(nmod)
      else
         eq(2,6) = 1.d0
      endif

***   surface density
************************************************************************

      if
     $     (ntprof.le.1.or.ntprof.eq.3.or.ntprof.eq.5.or.ntp2) then

         eq(3,16) = ro(nmod)-rosurf
         eq(3,1) = -rosurf1/gmrua*dynfac*dtipsi
         eq(3,6) = -rosurf1/gmrua*dynfac*(dtninv-dtipsi)
         eq(3,7) = 2.d0*rosurf1/gmrua*gmr(nmod)
         eq(3,8) = drodf(nmod)-rosurf*dmudf(nmod)*muinv(nmod)+rosurf1*
     &        dkapdf(nmod)/kap(nmod)
         eq(3,9) = drodt(nmod)-rosurf*dmudt(nmod)*muinv(nmod)+rosurf1*
     &        dkapdt(nmod)/kap(nmod)
         if (tsurfr.ne.teff) then
            eq(3,7) = eq(3,7)-0.5d0*rosurf1+1.5d0*rosurf2
            eq(3,10) = (0.25d0*rosurf1-0.75d0*rosurf2)/lum(nmod)
         endif

      else if (ntprof.eq.2) then

         eq(3,16) = ro(nmod)-rosurf
         eq(3,1) = -rosurf/gmrua*dynfac*dtipsi
         eq(3,6) = -rosurf/gmrua*dynfac*(dtninv-dtipsi)
         eq(3,7) = 2.d0*rosurf*gmr(nmod)/gmrua-0.5d0*rosurf
         eq(3,8) = drodf(nmod)-rosurf*dmudf(nmod)*muinv(nmod)+rosurf*
     &        dkapdf(nmod)/kap(nmod)*ledd/(ledd-lum(nmod))
         eq(3,9) = drodt(nmod)-rosurf*dmudt(nmod)*muinv(nmod)+rosurf*
     &        dkapdt(nmod)/kap(nmod)*ledd/(ledd-lum(nmod))
         eq(3,10) = rosurf/(ledd-lum(nmod))+0.25d0*rosurf/lum(nmod)

      elseif (ntprof.eq.4.or.ntprof.eq.6.or.ntprof.eq.7) then
         eq(3,16) = ro(nmod)-rosurf
         !eq(3,7) = 
     &   !     -3*lum(nmod)*kapm(nmod)*(1+phitau(nmod))/
     &   !     (64*pi*sig*t(nmod)**2*r(nmod))*ro(nmod)
         eq(3,8) = drodf(nmod)
         eq(3,9) = drodt(nmod)
      endif

***   energy conservation
************************************************************************

      if (lgrav.le.3) then
         eq(4,16) = -egrav(nmod1)-evisc(nmod1)+sigma*((lum(nmod)-
     &        lum(nmod1))/dm(nmod1)-enucl(nmod1))+sigma1*((vlum(nmod)-
     &        vlum(nmod1))/dm(nmod1)-venucl(nmod1))
         eq(4,1) = -devdu1(nmod1)
         eq(4,2) = -devdr1(nmod1)*r(nmod1)
         eq(4,3) = -degdf1(nmod1)-devdf1(nmod1)-sigma*denucldf(nmod1)
         eq(4,4) = -degdt1(nmod1)-devdt1(nmod1)-sigma*denucldt(nmod1)
         eq(4,5) = -sigma/dm(nmod1)
         eq(4,6) = -devdu2(nmod1)
         eq(4,7) = -devdr2(nmod1)*r(nmod)
         eq(4,8) = -degdf2(nmod1)-devdf2(nmod1)
         eq(4,9) = -degdt2(nmod1)-devdt2(nmod1)
         eq(4,10) = sigma/dm(nmod1)
      endif

      if (lgrav.ge.4) then
         drvfacc = facc(nmod1)*vro(nmod1)*dtninv/ro(nmod1)**2
***   case egrav = -dE/dt - 4*pi(P+Q)*d(r^2v)/dm
         eq(4,16) = (e(nmod1)-ve(nmod1)-omi(nmod)*(e(nmod)-e(nmod1)))*
     &        dtninv+sigma*(pim4*drvdm(nmod1)*p(nmod1)+(lum(nmod)-
     &        lum(nmod1))/dm(nmod1)-enucl(nmod1))+sigma1*(pim4*
     &        vdrvdm(nmod1)*vp(nmod1)+(vlum(nmod)-vlum(nmod1))/
     &        dm(nmod1)-venucl(nmod1))
     &        -(p(nmod1)+pvisc(nmod1))*drvfacc
         eq(4,1) = sigma*pim4*(drvdm(nmod)*dpvdu1(nmod1)-2.d0*(p(nmod1)+
     &        pvisc(nmod1))*r(nmod1)**2/dm(nmod1))
     &        -dpvdu1(nmod1)*drvfacc
         eq(4,2) = sigma*pim4*(drvdm(nmod)*dpvdr1(nmod1)-(p(nmod1)+
     &        pvisc(nmod1))*2.d0*r(nmod1)*u(nmod1)/dm(nmod1))
     &        -dpvdr1(nmod1)*drvfacc
         eq(4,2) = eq(4,2)*r(nmod1)
         eq(4,3) = dedf(nmod1)*(dtninv+dtiomi)+sigma*(pim4*drvdm(nmod)*
     &        (dpdf(nmod1)+dpvdf1(nmod1))-denucldf(nmod1))
     &        -(dpdf(nmod1)+dpvdf1(nmod1))*drvfacc+2.d0*(p(nmod1)+
     &        pvisc(nmod1))*drodf(nmod1)/ro(nmod1)
         eq(4,4) = dedt(nmod1)*(dtninv+dtiomi)+sigma*(pim4*drvdm(nmod)*
     &        (dpdt(nmod1)+dpvdt1(nmod1))-denucldt(nmod1))
     &        -(dpdt(nmod1)+dpvdt1(nmod1))*drvfacc+2.d0*(p(nmod1)+
     &        pvisc(nmod1))*drodt(nmod1)/ro(nmod1)
         eq(4,5) = -sigma/dm(nmod1)
         eq(4,6) = sigma*pim4*(p(nmod1)+pvisc(nmod1)*r(nmod)**2/
     &        dm(nmod1)+drvdm(nmod)*dpvdu2(nmod1))
     &        -dpvdu2(nmod1)*drvfacc
         eq(4,7) = sigma*pim4*(drvdm(nmod)*dpvdr2(nmod1)+(p(nmod1)+
     &        pvisc(nmod1))*2.d0*r(nmod)*u(nmod)/dm(nmod1))
     &        -dpvdr2(nmod1)*drvfacc
         eq(4,7) = eq(4,7)*r(nmod)
         eq(4,8) = dedf(nmod)*dtiomi+sigma*pim4*drvdm(nmod)*
     &        dpvdf2(nmod1)
     &        -dpvdf2(nmod1)*drvfacc
         eq(4,9) = dedt(nmod)*dtiomi+sigma*pim4*drvdm(nmod)*
     &        dpvdt2(nmod1)
     &        -dpvdt2(nmod1)*drvfacc
         eq(4,10) = sigma/dm(nmod1)
      endif

***   surface temperature
************************************************************************

      if (vlconv) then
         rrk = dsqrt(r(nmod)*r(nmod1))
         fac1 = 0.25d0*kap(nmod)*dm(nmod1)/(pim4*rrk*r(nmod))
         lrad = 16.d0*sig*pi*pw13*r(nmod)**2*t(nmod)**4/(pw23+fac1)

         eq(5,16) = lum(nmod)-lrad
         eq(5,2) = -0.5d0*lrad/(pw23+fac1)*fac1/r(nmod1)
         eq(5,2) = eq(5,2)*r(nmod1)
         eq(5,7) = -2.d0*lrad/r(nmod)-1.5d0*lrad/(pw23+fac1)*fac1/
     &        r(nmod)
         eq(5,7) = eq(5,7)*r(nmod)
         eq(5,8) = lrad/(pw23+fac1)*fac1*dkapdf(nmod)/kap(nmod)
         eq(5,9) = -lrad*(4.d0-fac1/(pw23+fac1)*dkapdt(nmod)/kap(nmod))
         eq(5,10) = 1.d0

c         lrad = pim4*sig*r(nmod)**2*teff**4
c         eq(5,16) = lum(nmod)-lrad
c         eq(5,7) = -2.d0*lrad
c         eq(5,10) = 1.d0

c         lrad = leff
c         eq(5,16) = lum(nmod)-lrad
c         eq(5,10) = 1.d0

c         fac1 = dsqrt(kap(nmod)*kap(nmod1))*0.5d0*dm(nmod1)/(pim4*
c     &        r(nmod)**2)
c         fac2 = tau0 + pw23
c         lrad = 16.d0*pi*sig*pw13*r(nmod)**2*t(nmod1)**4/(fac1+fac2)
c         eq(5,16) = lum(nmod)-lrad
c         eq(5,2) = 0.d0
c         eq(5,3) = lrad*fac1/(fac1+fac2)*dkapdf(nmod1)/(2.d0*kap(nmod1))
c         eq(5,4) = lrad*(fac1/(fac1+fac2)*dkapdt(nmod1)/(2.d0*
c     &        kap(nmod1))-4.d0)
c         eq(5,7) = -2.d0*lrad/r(nmod)*(2.d0*fac1+fac2)/(fac1+fac2)
c         eq(5,8) = lrad*fac1/(fac1+fac2)*dkapdf(nmod)/(2.d0*kap(nmod))
c         eq(5,9) = lrad*fac1/(fac1+fac2)*dkapdt(nmod)/(2.d0*kap(nmod))
c         eq(5,10) = 1.d0

      else if (ntprof.le.1.or.ntprof.ge.3) then
         eq(5,9) = t(nmod)
         eq(5,16) = t(nmod)-tsurfr*tsf
         if (tsurfr.ne.teff) then
            eq(5,7) = tsurfr*0.5d0*tsf
            eq(5,10) = -tsurfr*tsf*0.25d0/lum(nmod)
         endif

***   spherical case
c      eq(5,7) = tsurfr*0.5d0
c      eq(5,10) = -tsurfr*0.25d0/lum(nmod)

      else if (ntprof.eq.2) then
         tsurfr = teff
         eq(5,16) = t(nmod)-tsurfr
         eq(5,7) = 0.5d0*tsurfr
         eq(5,9) = t(nmod)
         eq(5,10) = -tsurfr/(4.d0*lum(nmod))
      endif

      return
      end
