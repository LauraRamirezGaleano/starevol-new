      SUBROUTINE inteq

************************************************************************
* Write the difference equations (interior) for stellar structure      *
* Modif CC -ST energie ondes soleil (5/10/07 --> 23/11/07)             *
* $LastChangedDate:: 2016-05-11 17:20:46 +0200 (Mer, 11 mai 2016)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 62                                                          $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
    !  include 'evolcom.conv'
      include 'evolcom.therm2'
      include 'evolcom.eng'
      include 'evolcom.eos'
      include 'evolcom.grad'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.opa'
      include 'evolcom.conv'
      include 'evolcom.rot'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer i,im,ip
      integer k,l

      logical vlconv,ntp2

      double precision dtninv,dmk,dp,sph,vrm,lrad,lconv
      double precision vdp,vvrm,vgmr,gmrua1,plo,dlnt1,dpvis,
     &     vdpvis,dmkb,dxsi,ddtau,sigma1,dt4,ablainv,dtiomi,
     &     dtipsi,dtipsi1,sigvrm,drvfacc
c      double precision dvfacc,efacc
      double precision wip,wjp,kiinv,kipinv,kiminv,sqrp
      double precision rhoat,Tatm

      common /resolution/ vlconv
      common /interpatm/ rhoat(nsh),Tatm(nsh),ntp2

      i = ish
      im = i-1
      ip = i+1
      sigma1 = 1.d0-sigma
      dtninv = 1.d0/dtn
      dmk = (dm(i)+dm(im))*0.5d0
      dmkb = (dm(i)+dm(ip))*0.5d0
      dp = p(i)-p(im)
      vdp = vp(i)-vp(im)
      sph = 3.d0/(pim4*ro(i))
      vrm = pim4*r(i)*r(i)/dmk
      vvrm = pim4*vr(i)*vr(i)/dmk
      vgmr = g*m(i)/(vr(i)*vr(i))
      dlnt1 = log(t(ip)/t(i))

      forall (k=1:neq,l=1:neq1) eq(k,l) = 0.d0
      wip = wi(ip)
      wjp = wj(ip)

      dpvis = pvisc(i)-pvisc(im)
      vdpvis = vpvisc(i)-vpvisc(im)
      sigvrm = sigma*vrm
      dtiomi = omi(i)*dtninv
      dtipsi = psi(i)*dtninv
      dtipsi1 = psi(ip)*dtninv
c      dtiomi = 0.d0
c      dtipsi = 0.d0
c      dtipsi1 = 0.d0

***   equation of motion
************************************************************************
!!! xsitau corrections not implemented !!!!

c      dvfacc = u(i)*facc(i)*dtninv*vro(i)/ro(i)
      eq(1,16) = accel(i)+sigma*(vrm*(dp+dpvis)+gmr(i)*fp(i))+sigma1*
     &     (vvrm*(vdp+vdpvis)+vgmr)
c     &     +dvfacc
      eq(1,1) = dynfac*dtipsi-sigvrm*dpvdu1(im)
      eq(1,2) = -sigvrm*dpvdr1(im)*r(im)
      eq(1,3) = -sigvrm*(dpdf(im)+dpvdf1(im))
      eq(1,4) = -sigvrm*(dpdt(im)+dpvdt1(im))
      eq(1,6) = dynfac*(dtninv-dtipsi)+sigvrm*(dpvdu1(i)-dpvdu2(im))
c     &     +dvfacc/u(i)
      eq(1,7) = 2.d0*sigma*(vrm*(dpvis+dp)-gmr(i)*fp(i))+sigvrm*
     &     (dpvdr1(i)-dpvdr2(im))*r(i)
      eq(1,8) = sigvrm*(dpdf(i)+dpvdf1(i)-dpvdf2(im))
c     &     -dvfacc*drodf(i)/ro(i)
      eq(1,9) = sigvrm*(dpdt(i)+dpvdt1(i)-dpvdt2(im))
c     &     -dvfacc*drodt(i)/ro(i)
      eq(1,11) = sigvrm*dpvdu2(i)
      eq(1,12) = sigvrm*dpvdr2(i)*r(ip)
      eq(1,13) = sigvrm*dpvdf2(i)
      eq(1,14) = sigvrm*dpvdt2(i)
      if (tau(i).lt.taulim.and.((ntprof.eq.1.or.ntprof.ge.3.and.
     &     .not.ntp2).or.ntprof.eq.6.or.ntprof.eq.7)) then
         if (crz(im)*crz(i).lt.0.and.crz(i).lt.-1.and.nretry.eq.4) then
            kiminv = 0.d0
            kiinv = 1.d0/kap(i)
         elseif (crz(i)*crz(ip).lt.0.and.crz(i).lt.-1.and.nretry.eq.4)
     &           then
            kiminv = 1.d0/kap(im)
            kiinv = 0.d0
         else
            if (numeric.eq.2.or.numeric.eq.4) then
               kiinv = kapm(i)/kap(i)**2
               kiminv = kapm(i)/kap(im)**2
            else
               kiinv = 1.d0/kap(i)
               kiminv = 1.d0/kap(im)
            endif
         endif
!     .. spherical corrections
         if ((ntprof.eq.4.or.ntprof.eq.6.or.ntprof.eq.7)
     &        .and.abs(phitau(i)).gt.1.d-6) then
            eq(1,16) = eq(1,16)+dptau(i)
            eq(1,3) = eq(1,3)+dptau(i)*dkapdf(im)*kiminv
            eq(1,4) = eq(1,4)+dptau(i)*dkapdt(im)*kiminv
            eq(1,7) = eq(1,7)
     &           -2.d0*gmr(i)*xsitau(i)
     &           +gmr(i)*dxsitaudr(i)
     &           !+pim4*r(i)**3*(gmr(i)+accel(i))*
     &           !(1+xsitau(i))/dmk*ro(i)
            eq(1,8) = eq(1,8)+gmr(i)*dxsitaudf(i)
            eq(1,9) = eq(1,9)+gmr(i)*dxsitaudt(i)
            eq(1,10) = eq(1,10)+gmr(i)*xsitau(i)/lum(i)
         end if
      endif

***   definition of lagrangian velocity
************************************************************************

      if (hydro) then
         eq(2,16) = (lnr(i)-vlnr(i)-psi(i)*(lnr(i)-lnr(i-1)))*dtninv-
     &        sigma*u(i)/r(i)-sigma1*vu(i)/vr(i)
         eq(2,2) = dtipsi
         eq(2,6) = -sigma/r(i)
         eq(2,7) = dtninv-dtipsi+sigma*u(i)/r(i)
      else
         eq(2,6) = 1.d0
      endif

***   mass conservation
************************************************************************

      eq(3,16) = (r(ip)**3-r(i)**3)/dm(i)-sph
      eq(3,7) = -3.d0*r(i)*r(i)*r(i)/dm(i)
c$$$      if (tau(i).lt.taulim.and.ntprof.eq.6) then
c$$$         eq(3,7) = eq(3,7)
c$$$     &        -sph*3*lum(i)*kapm(i)*(1+phitau(i))/(64*pi*sig*t(i)**2*
c$$$     &        r(i))
c$$$      end if 
      eq(3,8) = sph*drodf(i)/ro(i)
      eq(3,9) = sph*drodt(i)/ro(i)
      eq(3,12) = 3.d0*r(ip)*r(ip)*r(ip)/dm(i)

***   energy conservation
************************************************************************

      if (lgrav.le.3) then
c         efacc = u(im)**2*facc(im)*vro(im)*dtninv/ro(im)
         eq(4,16) = -egrav(im)-evisc(im)+sigma*((lum(i)-lum(im))/
     &        dm(im)-enucl(im))+sigma1*((vlum(i)-vlum(im))/dm(im)-
     &        venucl(im))-eloc(im)-egconv(im)
         eq(4,1) = -devdu1(im)
c     &        +2.d0*efacc/u(im)
         eq(4,2) = -devdr1(im)
         eq(4,2) = eq(4,2)*r(im)
         eq(4,3) = -degdf1(im)-devdf1(im)-sigma*denucldf(im)-delocdf(im)
     &        -degconvdf(im)
c     &        -efacc*drodf(im)/ro(im)
         eq(4,4) = -degdt1(im)-devdt1(im)-sigma*denucldt(im)-delocdt(im)
     &        -degconvdt(im)
c     &        -efacc*drodt(im)/ro(im)
         eq(4,5) = -sigma/dm(im)
         eq(4,6) = -devdu2(im)
         eq(4,7) = -devdr2(im)
         eq(4,7) = eq(4,7)*r(i)
         eq(4,8) = -degdf2(im)-devdf2(im)
         eq(4,9) = -degdt2(im)-devdt2(im)
         eq(4,10) = sigma/dm(im)
      endif

      if (lgrav.ge.4) then
c         efacc = u(im)**2*facc(im)*vro(im)*dtninv/ro(im)
         drvfacc = facc(im)*vro(im)*dtninv/ro(im)**2

***   case egrav = -dE/dt - 4*pi*(P+Q)*d(r^2v)/dm
         eq(4,16) = (e(im)-ve(im)-omi(i)*(e(i)-e(im)))*dtninv+
     &        sigma*((lum(i)-lum(im))/dm(im)+pim4*drvdm(i)*(p(im)+
     &        pvisc(im))-enucl(im))+sigma1*((vlum(i)-vlum(im))/dm(im)-
     &        venucl(im)+pim4*vdrvdm(i)*(vp(im)+vpvisc(im)))
     &        -(p(im)+pvisc(im))*drvfacc
c     &        +efacc
         eq(4,1) = sigma*pim4*(drvdm(i)*dpvdu1(im)-(p(im)+pvisc(im))*
     &        r(im)**2/dm(im))
     &        -dpvdu1(im)*drvfacc
c     &        +2.d0*efacc/u(im)
         eq(4,2) = sigma*pim4*(drvdm(i)*dpvdr1(im)-(p(im)+pvisc(im))*
     &        2.d0*r(im)*u(im)/dm(im))
     &        -dpvdr1(im)*drvfacc
         eq(4,2) = eq(4,2)*r(im)
         eq(4,3) = dedf(im)*(dtninv+dtiomi)+sigma*(pim4*drvdm(i)*
     &        (dpdf(im)+dpvdf1(im))-denucldf(im))
     &        -(dpdf(im)+dpvdf1(im))*drvfacc+2.d0*(p(im)+pvisc(im))*
     &        drodf(im)*drvfacc/ro(im)
c     &        -efacc*drodf(im)/ro(im)
         eq(4,4) = dedt(im)*(dtninv+dtiomi)+sigma*(pim4*drvdm(i)*
     &        (dpdt(im)+dpvdt1(im))-denucldt(im))
     &        -(dpdt(im)+dpvdt1(im))*drvfacc+2.d0*(p(im)+pvisc(im))*
     &        drodt(im)*drvfacc/ro(im)
c     &        -efacc*drodt(im)/ro(im)
         eq(4,5) = -sigma/dm(im)
         eq(4,6) = sigma*pim4*(drvdm(i)*dpvdu2(im)+(p(im)+pvisc(im))*
     &        r(i)**2/dm(im))
     &        -dpvdu2(im)*drvfacc
         eq(4,7) = sigma*pim4*(drvdm(i)*dpvdr2(im)+(p(im)+pvisc(im))*
     &        2.d0*r(i)*u(i)/dm(im))
     &        -dpvdr2(im)*drvfacc
         eq(4,7) = eq(4,7)*r(i)
         eq(4,8) = -dedf(i)*dtiomi+sigma*pim4*drvdm(i)*dpvdf2(im)
     &        -dpvdf2(im)*drvfacc
         eq(4,9) = -dedt(i)*dtiomi+sigma*pim4*drvdm(i)*dpvdt2(im)
     &        -dpvdt2(im)*drvfacc
         eq(4,10) = sigma/dm(im)
      endif

***   equation of transport
************************************************************************
      
      if (crz(i)*crz(ip).lt.0.and.crz(ip).lt.-1) then
         kiinv = 0.d0
         kipinv = 1.d0/kap(ip)/wjp
      elseif (crz(ip)*crz(min(ip+1,nmod)).lt.0.and.crz(ip).lt.-1) then
         kiinv = 1.d0/kap(i)/wip
         kipinv = 0.d0
      else
         if (numeric.eq.2.or.numeric.eq.4) then
            kiinv = kapm(ip)/kap(i)**2
            kipinv = kapm(ip)/kap(ip)**2
         else
            kiinv = 1.d0/kap(i)
            kipinv = 1.d0/kap(ip)
         endif
      endif
      dt4 = t(ip)**4-t(i)**4
      if (abs(dt4).lt.1.d-5) dt4 = dt4+1.d-5
      if (vlconv) then
         lrad = -64.d0*sig*pi**2*r(ip)**4*dt4/(3.d0*kapm(ip)*dmkb)
         lconv = vhconv(ip)*rom(ip)*tm(ip)*r(ip)

         eq(5,16) = lum(ip)-lrad-lconv
         eq(5,8) = wip*(lrad*dkapdf(i)*kiinv-lconv*drodf(i)/ro(i))
         eq(5,9) = lrad*(wip*dkapdt(i)*kiinv+4.d0*t(i)**4/dt4)-
     &        lconv*wip*(1.d0+drodt(i)/ro(i))
         eq(5,12) = -(4.d0*lrad+lconv)
         eq(5,13) = wjp*(lrad*dkapdf(ip)*kipinv-lconv*drodf(ip)/
     &        ro(ip))
         eq(5,14) = lrad*(wjp*dkapdt(ip)*kipinv-4.d0*t(ip)**4/dt4)-
     &        lconv*wjp*(1.d0+drodt(ip)/ro(ip))
         eq(5,15) = 1.d0
      else

***   version with xsitau and without psi and omi
!!!  TO BE CHECKED !!!!!!!
         if (tau(ip).lt.taulim.and.(ntprof.eq.1.or.(ntprof.ge.3.and.
     &        .not.ntp2))) then
            
 !           if (abs(dqdtau(ip)).gt.1d-3) then
               dxsi = -xsitau(i)/(g*m(ip)+r(ip)*r(ip)*accel(ip))
               ddtau = 0.d0
               if (abs(dqdtau(i)).gt.1.d-6) ddtau = ddqtau(i)/dqdtau(i)
c     plo = 3.d0*lum(ip)*kapm(ip)*dmkb/(64.d0*sig*pi*pim4*
c     &           (r(ip)*tm(ip))**4)*ft(ip)!/fp(ip)
!               write(*,*) ip, tau(ip), dqdtau(ip)
               gmrua1 = (gmr(ip)+accel(ip))*fp(ip)
               
               ablainv = 1.d0/abla(ip)
               sqrp = sqrt(p(i)*p(ip))
               
               plo = gmrua1*abla(ip)*dmkb/(pim4*r(ip)*r(ip)*sqrp)

               if ((ntprof.eq.4.or.ntprof.eq.6.or.ntprof.eq.7)
     &              .and.abs(phitau(i)).gt.1.d-6) then
                  
                  eq(5,16) = dlnt1+plo*(1.d0+xsitau(ip))
                  eq(5,6) = plo*(abdu1(ip)*ablainv+dynfac*dtninv
     &                 *psi(ip)/gmrua1)*(1.d0+xsitau(ip))+plo*dxsi
     &                 *r(ip)*r(ip)*dtninv*psi(ip)
                  !eq(5,7) = plo*r(i)/p(i)*
     &             !    gmrua1*(1+xsitau(i))*(1+xsitau(ip))*ro(i)
                  eq(5,8) = plo*(abdf1(ip)*ablainv-wip*dpdf(i)/p(i))*
     &                 (1.d0+xsitau(ip))
                  eq(5,9) = plo*(abdt1(ip)*ablainv-wip*dpdt(i)/p(i))*
     &                 (1.d0+xsitau(ip))-1.d0 
                  eq(5,11) = plo*(abdu2(ip)*ablainv+dynfac*(dtninv-
     &                 dtipsi1)/gmrua1)*(1.d0+xsitau(ip))
     &                 +plo*dxsi*r(ip)*
     &                 r(ip)*dtninv*(1.d0-psi(ip))
                  eq(5,12) = plo*(1.d0+xsitau(ip))*(abdr(ip)*ablainv
     &                 -2.d0*(1.d0+gmr(ip)/gmrua1))!+r(ip)/p(ip)*
     &                 !gmrua1*(1+xsitau(ip)*ro(ip)))*plo*(1+xsitau(ip))
     &                 +plo*dxsitaudr(ip)
                  eq(5,13) = plo*(abdf2(ip)*ablainv-wjp*dpdf(ip)/p(ip))*
     &                 (1.d0+xsitau(ip))+plo*dxsitaudf(ip)
                  eq(5,14) = plo*(abdt2(ip)*ablainv-wjp*dpdt(ip)/p(ip))*
     &                 (1.d0+xsitau(ip))+plo*dxsitaudt(ip)+1.d0
                  eq(5,15) = plo*(abdl(ip)*ablainv*(1.d0+xsitau(ip))+
     &                 xsitau(ip)/lum(ip))
                  
               else
                  
                  eq(5,16) = dlnt1+plo*(1.d0+xsitau(ip))
                  eq(5,6) = plo*(abdu1(ip)*ablainv+dynfac*dtninv*
     &                 psi(ip)/gmrua1)*(1.d0+xsitau(ip))+plo*dxsi*
     &                 r(ip)*r(ip)*dtninv*psi(ip)
                  eq(5,7) = plo*xsitau(ip)*ddtau*kap(i)*ro(i)
                  eq(5,8) = plo*(abdf1(ip)*ablainv-wip*dpdf(i)/p(i))*
     &                 (1.d0+xsitau(ip))+plo*xsitau(ip)*(wip*dkapdf(i)*
     &                 kiinv)
c     &           xsitau(ip))+plo*xsitau(ip)*(wip*dkapdf(i)*kiinv+ddtau*
c     &           dtaudf(i))
                  eq(5,9) = plo*(abdt1(ip)/abla(ip)-wip*dpdt(i)/p(i))*
     &                 (1.d0+xsitau(ip))+plo*xsitau(ip)*(wip*dkapdt(i)*
     &                 kiinv)  !+ddtau*
c     &           xsitau(ip))+plo*xsitau(ip)*(wip*dkapdt(i)*kiinv+ddtau*
c     &           dtaudt(i))
     &                 -1.d0
                  eq(5,11) = plo*(abdu2(ip)*ablainv+dynfac*(dtninv-
     &                 dtipsi1)/gmrua1)*(1.d0+xsitau(ip))+plo*dxsi
     &                 *r(ip)*r(ip)*dtninv*(1.d0-psi(ip))
                  eq(5,12) = plo*(abdr(ip)*ablainv-2.d0*(1.d0+gmr(ip)/
     &                 gmrua1)/r(ip))*(1.d0+xsitau(i))+plo*dxsi*2.d0*
     &                 r(ip)*accel(ip)
                  eq(5,13) = plo*(abdf2(ip)*ablainv-wjp*dpdf(ip)/p(ip))*
     &                 (1.d0+xsitau(ip))+plo*xsitau(ip)*wjp*dkapdf(ip)
     &                 *kipinv
                  eq(5,14) = plo*(abdt2(ip)/abla(ip)-wjp*dpdt(ip)/p(ip))
     &                 *(1.d0+xsitau(ip))+plo*xsitau(ip)*wjp*dkapdt(ip)
     &                 *kipinv+1.d0
                  eq(5,15) = plo*(abdl(ip)*ablainv*(1.d0+xsitau(ip))+
     &                 xsitau(ip)/lum(ip))
                  
               end if

***   version without atmospheric corrections (xsitau)
         else
c..  convective zones
            if (crz(ip).lt.-1) then
               gmrua1 = gmr(ip)*fp(ip)+accel(ip)
c               plo = gmrua1*abla(ip)*dmkb/(pim4*r(ip)*r(ip)*pm(ip))*
c     &              ft(ip)/fp(ip)
               plo = gmrua1*abla(ip)*dmkb/(pim4*r(ip)*r(ip)*pm(ip))
               ablainv = 1.d0/abla(ip)
               eq(5,16) = dlnt1+plo
               eq(5,6) = plo*(abdu1(ip)*ablainv+dynfac*dtipsi1/gmrua1)
               eq(5,8) = plo*(abdf1(ip)*ablainv-wip*dpdf(i)/p(i))
               eq(5,9) = plo*(abdt1(ip)*ablainv-wip*dpdt(i)/p(i))-1.d0
               eq(5,11) = plo*(abdu2(ip)*ablainv+dynfac*
     &              (dtninv-dtipsi1)/gmrua1)
               eq(5,12) = plo*(abdr(ip)*ablainv*r(ip)-2.d0*(1.d0+
     &              gmr(ip)/gmrua1))
               eq(5,13) = plo*(abdf2(ip)*ablainv-wjp*dpdf(ip)/p(ip))
               eq(5,14) = plo*(abdt2(ip)*ablainv-wjp*dpdt(ip)/p(ip))+
     &              1.d0
               eq(5,15) = plo*abdl(ip)*ablainv

            else
c..  radiative zones
c               plo = 3.d0*lum(ip)*kapm(ip)*dmkb/(64.d0*pi*sig*pim4*
c     &              (r(ip)*tm(ip))**4)*ft(ip)/fp(ip)
               plo = 3.d0*lum(ip)*kapm(ip)*dmkb/(64.d0*pi*sig*pim4*
     &              (r(ip)*tm(ip))**4)*ft(ip)
               eq(5,16) = dlnt1+plo
               eq(5,8) = plo*wip*dkapdf(i)*kiinv
               eq(5,9) = plo*wip*(dkapdt(i)*kiinv-4.d0)-1.d0
               eq(5,12) = -4.d0*plo
               eq(5,13) = plo*wjp*dkapdf(ip)*kipinv
               eq(5,14) = plo*wjp*(dkapdt(ip)*kipinv-4.d0)+1.d0
               eq(5,15) = plo/lum(ip)
c               plo = 3.d0*lum(ip)*kapm(ip)*dmkb/(16.d0*pi*sig*pim4*
c     &              r(ip)**4)*ft(ip)/fp(ip)
c               eq(5,16) = t(ip)**4-t(i)**4+plo
c               eq(5,8) = plo*wip*dkapdf(i)*kiinv
c               eq(5,9) = plo*wip*dkapdt(i)*kiinv-4.d0*t(i)**4
c               eq(5,12) = -4.d0*plo
c               eq(5,13) = plo*wjp*dkapdf(ip)*kipinv
c               eq(5,14) = plo*wjp*dkapdt(ip)*kipinv+4.d0*t(ip)**4
c               eq(5,15) = plo/lum(ip)
            endif
         endif
      endif

      return
      end
