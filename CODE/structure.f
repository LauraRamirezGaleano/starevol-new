      SUBROUTINE structure (error)

************************************************************************
*  Calculate egrav, viscous pressure, their dervatives and             *
*  related quantities related to the equation of state (eos)           *
*  Modif CC - ST energie ondes soleil (5/10/07 --> 23/11/07)           *
* $LastChangedDate:: 2014-02-04 15:45:03 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 11                                                          $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.eos'
      include 'evolcom.grad'
      include 'evolcom.ion'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.opa'
      include 'evolcom.rot'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.var'
      include 'evolcom.igw' ! LR 20250616


      integer error,jpass
      integer imin,imax
      integer i,j,k,l,ip,kl
      integer itmax,ishockb,ishockt
      integer kc,kc1,kc2,iedge,imid
      integer fcs         ! first convective shell
C Modif CC - ST energie ondes soleil (5 octobre 2007)
      integer nenvconv
      integer klenvstruc,klcorestruc,ndbstruc,ndtstruc
C <--

      double precision domdr(nsh),d2omdr(nsh),dmdr(nsh),djdt(nsh),
     &     sigm(nsh),vsigm(nsh),djdtr(nsh)
      double precision eshrdro1,eshrdro2,tvt,pvp,dtiomi,dtipsi
      double precision kiinv,kipinv,piinv,pipinv,pro2inv,
     &     pro2ipinv,dro1,devisc,devisc1,dpvis,dpvis1
      double precision dtninv,gmrr2,deltamu,dlogmu
      double precision cortau1,cortau2
      double precision Lmax,tmax,Mackmax,enucmax
      double precision q1,rim,dvr,dpvis2
C Modif CC - ST energie ondes soleil (5 octobre 2007)
      double precision renvconv,den
      double precision eondestot,eondestotzr,eondestotec
C <--
      double precision fe(nsh),efermi(nsh),efermim(nsh),defermi(nsh)
      double precision mueconv,efs,mconv
      double precision dvaria,varia

      double precision rhoat,Tatm
      double precision geffrot,Rp2
      double precision disctime   ! ajout TD Fev 2019

      logical vlconv,lvisc,lshock
      logical partialmix,neutrality,bp
      logical ntp2
c      logical ecap

      common /disclocking/ disctime ! ajout TD Fev 2019
      common /hydrodyn/ Lmax,tmax,Mackmax,enucmax,itmax,ishockb,ishockt
      common /resolution/ vlconv
      common /geffective/ geffrot(nsh),Rp2(nsh)

      common /interpatm/ rhoat(nsh),Tatm(nsh),ntp2


*__________________________________________________________________
***   Calculate the atmospheric temperature profile and derivatives
*------------------------------------------------------------------

      call atmtp

*____________________________________________________________________
***   calculation of the mean thermodynamic values, the gravitational
***   energy production rate, the shear energy production rate
***   (eventually) and the radiative gradient
*--------------------------------------------------------------------
      
      bp = nphase.eq.4.and.xsp(1,ihe4).lt.0.5d0.and.xsp(1
     $     ,ihe4).gt.0.003d0

      lvisc = hydro.and.q0.gt.0.d0.and.ivisc.gt.0
      dmdr(1) = 0.d0
      domdr(1) = 0.d0
      d2omdr(1) = 0.d0
      if (iaccr.gt.0) then
         sigm(1) = facc(1)*dm(1)*sigma0
         vsigm(1) = vfacc(1)*dm(1)*sigma0
      endif
      omi(1) = 0.d0
      psi(1) = 0.d0
      accel(1) = 0.d0
      dtninv = 1.d0/dtn

      if (hydro) then
         do i = 2,nmod
            accel(i) = (u(i)-vu(i)-psi(i)*(u(i)-u(i-1)))*dtninv
         enddo
      else
         accel(2:nmod) = 0.d0
      endif

c..   compute mean thermodynamic values
      tm(1) = t(1)
      pm(1) = p(1)
      rom(1) = ro(1)
      kapm(1) = kap(1)
      cpm(1) = cp(1)
      deltaKSm(1) = deltaKS(1)
      abm(1) = abad(1)
      drvdm(1) = 0.d0
cx!$OMP PARALLEL
cx!$OMP DO PRIVATE(ip)
      do i = 1,nmod1
         ip = i+1
         rom(ip) = ro(i)**wi(ip)*ro(ip)**wj(ip)
         tm(ip) = exp(wi(ip)*lnt(i)+wj(ip)*lnt(ip))
         pm(ip) = p(i)**wi(ip)*p(ip)**wj(ip)
         cpm(ip) = cp(i)**wi(ip)*cp(ip)**wj(ip)
         deltaKSm(ip) = deltaKS(i)**wi(ip)*deltaKS(ip)**wj(ip)
         abm(ip) = abad(i)**wi(ip)*abad(ip)**wj(ip)
         if (numeric.eq.2.or.numeric.eq.4) then
            kapm(ip) = kap(i)*kap(ip)/(wj(ip)*kap(i)+wi(ip)*kap(ip))
         else
            kapm(ip) = kap(i)**wi(ip)*kap(ip)**wj(ip)
         endif
c         if (.not.hydrorot) gmr(ip) = g*m(ip)/(r(ip)*r(ip))
         gmr(ip) = g*m(ip)/(r(ip)*r(ip))
         hp(ip) = pm(ip)/rom(ip)/abs(gmr(ip)+accel(ip))
c         hp(ip) = pm(ip)/rom(ip)/gmr(ip)
      enddo
cx!$OMP END DO
cx!$OMP END PARALLEL

c..   redefine thermo at convective borders
      if (nsconv.gt.0.and.nretry.eq.4) then
         do kl = 3,4
            if (kl.eq.3) then
               iedge = 1
               imid = 0
            else
               iedge = -1
               imid = -1
            endif
            do k = 1,nsconv
               kc = novlim(k,kl)
               if (kc.gt.2) then
                  kc1 = kc+iedge
                  kc2 = kc1+iedge
                  if (iedge.eq.1) then
                     dvaria = -dm(kc)/(dm(kc)+dm(kc1))
                  else
                     dvaria = -dm(kc1)/(dm(kc2)+dm(kc1))
                  endif
                  varia = dvaria*(kap(kc1+imid)-kap(kc+imid))+
     &                 kap(kc+imid)
                  if (varia.gt.0.d0) kapm(kc) = varia
                  varia = dvaria*(p(kc1+imid)-p(kc+imid))+
     &                 p(kc+imid)
                  if (varia.gt.0.d0) pm(kc) = varia
                  varia = dvaria*(t(kc1+imid)-t(kc+imid))+
     &                 t(kc+imid)
                  if (varia.gt.0.d0) tm(kc) = varia
                  varia = dvaria*(ro(kc1+imid)-ro(kc+imid))+
     &                 ro(kc+imid)
                  if (varia.gt.0.d0) rom(kc) = varia
                  varia = dvaria*(abad(kc1+imid)-abad(kc+imid))+
     &                 abad(kc+imid)
                  if (varia.gt.0.d0) abm(kc) = varia
                  varia = dvaria*(cp(kc1+imid)-cp(kc+imid))+
     &                 cp(kc+imid)
                  if (varia.gt.0.d0) cpm(kc) = varia
                  varia = dvaria*(deltaKS(kc1+imid)-deltaKS(kc+imid))+
     &                 deltaKS(kc+imid)
                  if (varia.gt.0.d0) deltaKSm(kc) = varia
                  hp(kc) = pm(kc)/rom(kc)/abs(gmr(kc)+accel(kc))
               endif
            enddo
         enddo
      endif

c..   URCA terms
      if (urca) then
         efs = dsqrt(1.0d0+1.018d-4*(ro(1)*mueinv(1))**pw23)
         efermi(1) = 8.187578d-7*(efs-1.d0)
         efermim(1) = efermi(1)
         defermi(1) = 2.778318d-11*efs*mueinv(1)**pw23/ro(1)**pw13
         fe(1) = 0.d0
         do k = 2,nmod
            fe(k) = 0.d0
            efs = dsqrt(1.d0+1.018d-4*(ro(k)*mueinv(k))**pw23)
            efermi(k) = 8.187578d-7*(efs-1.d0)+1.d-90
            defermi(k) = 2.778318d-11*efs*mueinv(k)**pw23/ro(k)**pw13
            efermim(k) = wi(k)*efermi(k-1)+wj(k)*efermi(k)
            eloc(k) = -avn*efermi(k)*(mueinv(k)-vmueinv(k))*dtninv
            delocdf(k) = eloc(k)/efermi(k)*defermi(k)*drodf(k)
            delocdt(k) = eloc(k)/efermi(k)*defermi(k)*drodt(k)
            if (crz(k).ge.0) then
               egconv(k) = 0.d0
               degconvdt(k) = 0.d0
               degconvdf(k) = 0.d0
            endif
c            mueinvm = mueinv(i-1)**wi(i)+mueinv(i)**wj(i)
c            efsm = dsqrt(1.0d0+1.018d-4*(rom(k)*mueinvm)**pw23)
c            efermim(k) = 8.187578d-7*(efsm-1.d0)
         enddo
         if (nsconv.gt.0) then
            do kl = 1,nsconv
c               ecap = ro(novlim(kl,3)+1).gt.1.d8
c               if (ecap) then
               mueconv = 0.d0
               mconv = 0.d0
               do k = novlim(kl,3)+1,novlim(kl,4)
                  mueconv = mueconv+(mueinv(k)-vmueinv(k))*dm(k)
                  mconv = mconv+dm(k)
               enddo
               mueconv = mueconv/mconv
               do i = novlim(kl,3)+1,novlim(kl,4)
                  ip = i+1
                  fe(i) = fe(i-1)+(mueinv(i)-vmueinv(i)-mueconv)*
     &                 dtninv*dm(i)
                  efs =  avn*fe(i)/dm(i)
                  egconv(i) = efs*(efermim(ip)-efermim(i))
                  degconvdf(i) = efs*(wj(ip)*defermi(ip)*drodf(ip)+
     &                 defermi(i)*drodf(i)*(wi(ip)-wj(i))-
     &                 wi(i)*defermi(i-1)*drodf(i-1))
                  degconvdt(i) = efs*(wj(ip)*defermi(ip)*drodt(ip)+
     &                 defermi(i)*(wi(ip)*drodt(ip)-wj(i)*drodt(i))
     &                 -wi(i)*defermi(i-1)*drodt(i-1))
c                degconvdt0(i) = -efs*wi(i)*defermi(i-1)*drodt(i-1)
c                degconvdf0(i) = -efs*wi(i)*defermi(i-1)*drodf(i-1)
c                degconvdt1(i) = efs*(wi(ip)-wj(i))*defermi(i)*drodt(i)
c                degconvdf1(i) = efs*(wi(ip)-wj(i))*defermi(i)*drodf(i)
c                degconvdt2(i) = efs*wj(ip)*defermi(ip)*drodt(ip)
c                degconvdf2(i) = efs*wj(ip)*defermi(ip)*drodf(ip)
               enddo
            enddo
         endif
      else
         forall (i=1:nmod)
            eloc(i) = 0.d0
            delocdf(i) = 0.d0
            delocdt(i) = 0.d0
            egconv(i) = 0.d0
            degconvdt(i) = 0.d0
            degconvdf(i) = 0.d0
         end forall
      endif

!$OMP PARALLEL
!$OMP DO PRIVATE(ip) SCHEDULE(dynamic,10)
      do i = 1,nmod1
         ip = i+1
         dtiomi = omi(ip)*dtninv
         dtipsi = psi(ip)*dtninv
c         dtiomi = 0.d0
c         dtipsi = 0.d0
         piinv = 1.d0/p(i)
         pipinv = 1.d0/p(ip)

***   compute D(1/ro)/Dt
         if (lgrav.le.1.or.(hydro.and.ivisc.gt.0.and.q0.gt.0.d0))
     &        dvdti(i) = (1.d0/ro(i)-1.d0/vro(i)-omi(ip)*(1.d0/ro(ip)
     &        -1.d0/ro(i)))*dtninv

***   Egrav = -P D(1/ro)/Dt-De/Dt
!         if (lgrav.le.1) then
         if (lgrav.le.1.and..not.bp) then
            pvp = p(i)
c1            pvp = 0.5d0*(vp(i)+p(i))
c2            pvp = dsqrt(vp(i)*p(i))
            pro2inv = pvp/(ro(i)*ro(i))
            pro2ipinv = pvp/(ro(ip)*ro(ip))
            egrav(i) = -(pvp*dvdti(i)+(e(i)-ve(i)-omi(ip)*(e(ip)-e(i)))*
     &           dtninv)
            degdf1(i) = -(dedf(i)-pro2inv*drodf(i))*(dtiomi+dtninv)
     &           -dpdf(i)*dvdti(i)
c1     &           -0.5d0*dpdf(i)*dvdti(i)
c2     &          -0.5d0*pvp*piinv*dpdf(i)*dvdti(i)
            degdt1(i) = -(dedt(i)-pro2inv*drodt(i))*(dtiomi+dtninv)
     &           -dpdt(i)*dvdti(i)
c1     &           -0.5d0*dpdt(i)*dvdti(i)
c2     &          -0.5d0*pvp*piinv*dpdt(i)*dvdti(i)
            degdf2(i) = dtiomi*(dedf(ip)-pro2ipinv*drodf(ip))
            degdt2(i) = dtiomi*(dedt(ip)-pro2ipinv*drodt(ip))
            
         elseif (lgrav.le.1.and.bp) then
            egrav(i) = 0.d0
            degdf1(i) = 0.d0
            degdt1(i) = 0.d0
            degdf2(i) = 0.d0
            degdt2(i) = 0.d0
            
***   Egrav = -cp DT/Dt + (cp*T*abad)/P * DP/Dt
         elseif (lgrav.eq.2) then
            dpdti(i) = (p(i)-vp(i)-omi(ip)*(p(ip)-p(i)))*dtninv
            dtdti(i) = (t(i)-vt(i)-omi(ip)*(t(ip)-t(i)))*dtninv
            egrav(i) = -cp(i)*(dtdti(i)-dpdti(i)*abad(i)*t(i)*piinv)
            degdf1(i) = egrav(i)*dcpdf(i)/cp(i)+cp(i)*t(i)*piinv*
     &           (dpdti(i)*dabadf(i)+abad(i)*dpdf(i)*(dtninv+dtiomi-
     &           dpdti(i)*piinv))
            degdt1(i) = egrav(i)*(dcpdt(i)/cp(i)+1.d0)-cp(i)*t(i)*
     &           (dtninv+dtiomi-dtdti(i)/t(i)-(dpdti(i)*dabadt(i)+
     &           abad(i)*dpdt(i)*(dtninv+dtiomi-dpdti(i)*piinv))*piinv)
            degdf2(i) = -dtiomi*cp(i)*t(i)*dpdf(ip)*abad(i)*piinv
            degdt2(i) = -dtiomi*cp(i)*(t(i)*dpdt(ip)*abad(i)*piinv-
     &           t(ip))

***   Egrav = -T DS/Dt
         elseif (lgrav.eq.3) then
c0            tvt = t(i)
            tvt = 0.5d0*(t(i)+vt(i))
c2            tvt = dsqrt(t(i)*vt(i))
            egrav(i) = -tvt*(s(i)-vs(i)-omi(ip)*(s(ip)-s(i)))*dtninv
            degdf1(i) = -tvt*dsdf(i)*(dtninv+dtiomi)
c0            degdt1(i) = egrav(i)-tvt*dsdt(i)*(dtninv+dtiomi)
            degdt1(i) = 0.5d0*egrav(i)-tvt*dsdt(i)*(dtninv+dtiomi)
c2            degdt1(i) = 0.5d0*egrav(i)-tvt*dsdt(i)*(dtninv+dtiomi)
            degdf2(i) = dtiomi*tvt*dsdf(ip)
            degdt2(i) = dtiomi*tvt*dsdt(ip)

***   Egrav = -(P+Q)/ro*div(v)-De/Dt = -4*pi*(P+Q)*d(r^2v)/dm-De/Dt
         elseif (lgrav.ge.4) then
            egrav(i) = -pim4*p(i)*(r(ip)**2*u(ip)-r(i)**2*u(i))/
     &           dm(i)-(e(i)-ve(i)-omi(ip)*(e(ip)-e(i)))*dtninv
     &           -pim4*p(i)*facc(i)*vro(i)/(ro(i)*ro(i))*dtninv
         endif

         if ((ivisc.eq.2.and.q0.gt.0.d0).or.lgrav.ge.4) then
            drvdm(ip) = (r(ip)**2*u(ip)-r(i)**2*u(i))/dm(i)
            vdrvdm(ip) = (vr(ip)**2*vu(ip)-vr(i)**2*vu(i))/dm(i)
         else
            drvdm(ip) = 1.d0
            vdrvdm(ip) = 1.d0
         endif


***   treatment of shocks
         lshock = lvisc.and.i.ge.ishockb-nq0.and.i.le.ishockt+nq0
c         lshock = .true.
c         icut = lvisc.and.abs(i-ishockb).le.10.and.abs(i-ishockt).le.10

         if (lshock) then
            jpass = 0

***   viscous pressure

***   version  Q = q0*l^2*rho*(d ln(rho)/dt)**2
            if (ivisc.eq.1.and.dvdti(i).lt.0.d0) then
               jpass = 1
c                dro1 = -dvdti(i)*ro(i)
c                pvisc(i) = q0*ro(i)*((r(ip)-r(i))*dro1)**2
c                dpvis = pvisc(i)*(3.d0+2.d0*(dtninv+dtiomi)/dro1)/
c     &               ro(i)
c                dpvis1 = -2.d0*dtiomi*pvisc(i)*ro(i)/(dro1*ro(ip)*
c     &               ro(ip))
               dro1 = (log(ro(i)/vro(i))-omi(ip)*log(ro(ip)/ro(i)))*
     &              dtninv
               pvisc(i) = q0*ro(i)*((r(ip)-r(i))*dro1)**2
               dpvis = pvisc(i)*(1.d0+2.d0*(dtninv+dtiomi)/dro1)/ro(i)
               dpvis1 = -2.d0*dtiomi*pvisc(i)/(dro1*ro(ip))
               dpvdu1(i) = 0.d0
               dpvdr1(i) = -2.d0*pvisc(i)/(r(ip)-r(i))
               dpvdf1(i) = dpvis*drodf(i)
               dpvdf2(i) = dpvis1*drodf(ip)
               dpvdt1(i) = dpvis*drodt(i)
               dpvdt2(i) = dpvis1*drodt(ip)
               dpvdr2(i) = -dpvdr1(i)
               dpvdu2(i) = 0.d0

***   version Q = q0*l^2*rho*(div U)**2
            elseif (ivisc.eq.2.and.drvdm(ip).lt.0.d0) then
               jpass = 1
               pvisc(i) =q0*ro(i)*(pim4*ro(i)*(r(ip)-r(i))*drvdm(ip))**2
               dpvis = 3.d0*pvisc(i)/ro(i)
               dpvis1 = 2.d0*pvisc(i)/(drvdm(ip)*dm(i))
               dpvdu1(i) = -dpvis1*r(i)**2
               dpvdr1(i) = -2.d0*pvisc(i)/(r(ip)-r(i))-2.d0*dpvis1*
     &              r(i)*u(i)
               dpvdf1(i) = dpvis*drodf(i)
               dpvdt1(i) = dpvis*drodt(i)
               dpvdf2(i) = 0.d0
               dpvdt2(i) = 0.d0
               dpvdu2(i) = dpvis1*r(ip)**2
               dpvdr2(i) = 2.d0*pvisc(i)/(r(ip)-r(i))+2.d0*dpvis1*
     &              r(ip)*u(ip)

***   version Q = q0*l^2*rho*|div U|*(dU/dr-1/3*div U)
            elseif (ivisc.eq.3.and.drvdm(ip).lt.0.d0) then
               jpass = 1
               q1 = pw23*q0*pim4*pim4
               rim = (r(ip)**3+r(i)**3)*0.5d0
               dvr = (u(ip)/r(ip)-u(i)/r(i))/dm(i)
               pvisc(i) = q1*ro(i)**3*rim*(r(ip)-r(i))**2*dvr*drvdm(ip)
               if (pvisc(i).lt.0.d0) write (nout,10) pvisc(i),i
               dpvis = 3.d0*pvisc(i)/ro(i)
               dpvis1 = pvisc(i)/(drvdm(ip)*dm(i))
               dpvis2 = pvisc(i)/(dvr*dm(i))
               dpvdu1(i) = -dpvis1*r(i)**2-dpvis2/r(i)
               dpvdr1(i) = -2.d0*pvisc(i)/(r(ip)-r(i))-2.d0*dpvis1*r(i)*
     &              u(i)+1.5d0*pvisc(i)*r(i)**2/rim+dpvis2*u(i)/r(i)**2
               dpvdf1(i) = dpvis*drodf(i)
               dpvdt1(i) = dpvis*drodt(i)
               dpvdf2(i) = 0.d0
               dpvdt2(i) = 0.d0
               dpvdu2(i) = dpvis1*r(ip)**2+dpvis2/r(ip)
               dpvdr2(i) = 2.d0*pvisc(i)/(r(ip)-r(i))+2.d0*dpvis1*r(ip)*
     &              u(ip)+1.5d0*pvisc(i)*r(ip)**2/rim-dpvis2*u(ip)/
     &              r(ip)**2
            endif
            if (jpass.eq.0.or.pvisc(i).lt.0.d0) then
               pvisc(i) = 0.d0
               dpvdu1(i) = 0.d0
               dpvdr1(i) = 0.d0
               dpvdf1(i) = 0.d0
               dpvdt1(i) = 0.d0
               dpvdf2(i) = 0.d0
               dpvdt2(i) = 0.d0
               dpvdr2(i) = 0.d0
               dpvdu2(i) = 0.d0
            endif

***   computation of Evisc (and derivatives)
            if (lgrav.le.3) then
               devisc = pvisc(i)*(dtninv+dtiomi)/(ro(i)*ro(i))
               devisc1 = pvisc(i)*dtiomi/(ro(ip)*ro(ip))
               evisc(i) = -pvisc(i)*dvdti(i)
               devdu1(i) = -dpvdu1(i)*dvdti(i)
               devdr1(i) = -dpvdr1(i)*dvdti(i)
               devdf1(i) = -dpvdf1(i)*dvdti(i)+devisc*drodf(i)
               devdf2(i) = -dpvdf2(i)*dvdti(i)-devisc1*drodf(ip)
               devdt1(i) = -dpvdt1(i)*dvdti(i)+devisc*drodt(i)
               devdt2(i) = -dpvdt2(i)*dvdti(i)-devisc1*drodt(ip)
               devdu2(i) = -dpvdu2(i)*dvdti(i)
               devdr2(i) = -dpvdr2(i)*dvdti(i)
            else
               evisc(i) = -pim4*drvdm(i)*pvisc(i)
            endif
         else
            pvisc(i) = 0.d0
            evisc(i) = 0.d0
            devdu1(i) = 0.d0
            devdr1(i) = 0.d0
            devdf1(i) = 0.d0
            devdt1(i) = 0.d0
            devdu2(i) = 0.d0
            devdr2(i) = 0.d0
            devdf2(i) = 0.d0
            devdt2(i) = 0.d0
            dpvdu1(i) = 0.d0
            dpvdr1(i) = 0.d0
            dpvdf1(i) = 0.d0
            dpvdt1(i) = 0.d0
            dpvdr2(i) = 0.d0
            dpvdu2(i) = 0.d0
            dpvdf2(i) = 0.d0
            dpvdt2(i) = 0.d0
         endif


***   treatment of shear energy (!!! not checked !!!)
***   discretization wrong !!!!!!!!
         eshr(i) = 0.d0
         eshrdr1(i) = 0.d0
         eshrdr2(i) = 0.d0
         eshrdf1(i) = 0.d0
         eshrdf2(i) = 0.d0

         if (iaccr.ge.3.and.model.ge.1.and.accphase.eq.0) then
            dmdr(ip) = pim4*r(ip)*r(ip)*ro(ip)
            sigm(ip) = sigm(i)+dmdr(ip)*r(ip)*r(ip)*(r(ip)-r(i))*
     &           omega(ip)
            vsigm(ip) = vsigm(i)+pim4*vr(ip)**4*(vr(ip)-vr(i))*
     &           vro(ip)*vomega(ip)
            if (i.gt.1.and.crz(i).gt.0.and.crz(ip).gt.0) then
               domdr(i) = (omega(i)-omega(i-1))/(r(i)-r(i-1))
c               d2omdr(i) = (domdr(ip)-domdr(i))/(r(ip)-r(i))
               d2omdr(i) = 0.d0
               dmdr(ip) = pim4*r(ip)*r(ip)*ro(ip)
               sigm(ip) = sigm(i)+dmdr(ip)*r(ip)*r(ip)*
     &              (r(ip)-r(i))*omega(ip)
               vsigm(ip) = vsigm(i)+pim4*vr(ip)**4*(vr(ip)-vr(i))*
     &              vro(ip)*vomega(ip)
               djdt(i) = (sigm(i)-vsigm(i)+omi(ip)*(sigm(ip)-
     &              sigm(i)))*dtninv*pw23
               djdtr(i) = (dmdr(i)*r(i)*((4.d0*omega(i)+r(i)*domdr(i))*
     &              (r(i)-r(i-1))+omega(i)*r(i))-omi(ip)*pim4*
     &              ro(ip)*r(ip)**4*omega(ip))*dtninv
               eshr(i) = domdr(i)*djdt(i)/dmdr(i)
               eshrdr1(i) = (djdtr(i)*domdr(i)+djdt(i)*d2omdr(i))/
     &              dmdr(i)-djdt(i)*domdr(i)*2.d0/r(i)/dmdr(i)
               eshrdr2(i) = pim4*ro(ip)*r(ip)**3*((5.d0*r(ip)-
     &              4.d0*r(i))*omega(ip)+r(ip)*domdr(ip)*(r(ip)-
     &              r(i)))*omi(ip)*dtninv*domdr(i)/dmdr(i)
               eshrdro1 = (pim4*r(i)**4*(r(i)-r(i-1))*omega(i)*
     &              dtninv-djdt(i)/ro(i))*domdr(i)/dmdr(i)
               eshrdro2 = pim4*r(ip)**4*(r(ip)-r(i))*omega(ip)*
     &              domdr(i)/dmdr(i)*omi(ip)*dtninv
               eshrdf1(i) = eshrdro1*drodf(i)
               eshrdf2(i) = eshrdro2*drodf(ip)
               eshrdt1(i) = eshrdro1*drodt(i)
               eshrdt2(i) = eshrdro2*drodt(ip)
            endif
         endif

***   define gradients (interface variables)

         dlogmu = log(mui(ip))-log(mui(i))
         deltamu = (mu(i)*muiinv(i))**wi(ip)*(mu(ip)*muiinv(ip))**wj(ip)
     &        *dlogmu
         abmu(ip) = deltamu/(log(p(ip)/p(i))+1.d-20) !d ln mu/ d ln P
c         abmuKS(ip) = abmu(ip)*sign(abs(phiKS(i)/deltaKS(i))**wi(ip)*
c     &        abs(phiKS(ip)/deltaKS(ip))**wj(ip),phiKS(i)/deltaKS(i))
         abmuKS(ip) = abmu(ip)*(phiKS(i)/deltaKS(i)*wi(ip)+
     &        phiKS(ip)/deltaKS(ip)*wj(ip))

         abled(ip) = abm(ip)+abmuKS(ip)

         gmrr2 = g*m(ip)+r(ip)*r(ip)*accel(ip)
         abrad(ip) = 3.d0*kapm(ip)*lum(ip)*pm(ip)/(64.d0*pi*sig*
     &        tm(ip)**4*gmrr2)*ft(ip)/fp(ip)
         abla(ip) = abrad(ip)
         if (numeric.eq.2.or.numeric.eq.4) then
            kiinv = kapm(ip)/kap(i)**2
            kipinv = kapm(ip)/kap(ip)**2
         else
            kiinv = 1.d0/kap(i)
            kipinv = 1.d0/kap(ip)
         endif

*________________________________________________________________________
***   Corrections to the radiative gradient in atmospheric layers
***   when plane-parallel, grey and Eddington approximations are ruled out
*------------------------------------------------------------------------

         xsitau(ip) = 0.d0
         phitau(ip) = dqdtau(ip)+2.d0*(qtau(ip)+
     &        tau(ip))/(r(ip)*kapm(ip)*rom(ip))
         
         if (tau(ip).lt.taulim.and.((ntprof.eq.1.or.ntprof.ge.3.and.
     &        .not.ntp2).or.ntprof.eq.6.or.ntprof.eq.7)) then
            
***   Plane-parallel corrections
*-------------------------------
            
c$$$            dptau(ip) = -kapm(ip)*lum(ip)/(pim4*c*r(ip)**2)*dqdtau(ip)
c$$$            xsitau(ip) = -dptau(ip)*r(ip)**2/gmrr2
c$$$  xsitau(ip) = 0.d0   ! TEST 31/03/2020
!            print*, dqdtau(ip), xsitau(ip)

***   
c$$$            if (abs(dqdtau(ip)).gt.1d-1) then
c$$$               abrad(ip) = abrad(ip)*(1.d0+dqdtau(ip))/(1.d0+xsitau(ip))
c$$$               dxsitaudt(ip) = xsitau(ip)*(kipinv*dkapdt(ip)
c$$$     &              + ddqdtaudt(ip)/dqdtau(ip))
c$$$               abdt1(ip) = abrad(ip)*wi(ip)*(kiinv*dkapdt(i)+
c$$$     &              piinv*dpdt(i)-4.d0 + ddqdtaudt(ip)/(1.d0+dqdtau(ip))
c$$$     &              - dxsitaudt(ip)/(1.d0+xsitau(ip)))
c$$$               abdt2(ip) = abrad(ip)*wj(ip)*(kipinv*dkapdt(ip)+pipinv*
c$$$     &              dpdt(ip)-4.d0 + ddqdtaudt(ip)/(1.d0+dqdtau(ip))
c$$$     &              - dxsitaudt(ip)/(1.d0+xsitau(ip)))
c$$$            else
c$$$               abdf1(ip) = abrad(ip)*wi(ip)*(kiinv*dkapdf(i)+
c$$$     &              piinv*dpdf(i))
c$$$               abdf2(ip) = abrad(ip)*wj(ip)*(kipinv*dkapdf(ip)+pipinv*
c$$$     &              dpdf(ip))
c$$$               abdt1(ip) = abrad(ip)*wi(ip)*(kiinv*dkapdt(i)+
c$$$     &              piinv*dpdt(i)-4.d0)
c$$$               abdt2(ip) = abrad(ip)*wj(ip)*(kipinv*dkapdt(ip)+pipinv*
c$$$     &              dpdt(ip)-4.d0)
c$$$            end if
c$$$            abdl(ip) = abrad(ip)/lum(ip)
c            write(555,'(2(1x,i4),7(1x,1pe11.4))'),i,crz(i)
c     $        ,dqdtau(i),ddqtau(i),xsitau(i),abrad(i),tau(i),t(i)

***   
            if ((ntprof.eq.4.or.ntprof.eq.6.or.ntprof.eq.7)
     &        .and.dqdtau(ip).gt.1.d-6) then

               phitau(ip) = dqdtau(ip)
               dptau(ip) = -kapm(ip)*lum(ip)/(pim4*c*r(ip)**2)*
     &              phitau(ip)
               xsitau(ip) = -dptau(ip)*r(ip)**2/gmrr2
               dphitaudr(ip) = 0.d0 
               dphitaudf(ip) = 0.d0 
               dphitaudt(ip) = 0.d0 
               dxsitaudr(ip) = xsitau(ip)*(dphitaudr(ip)/phitau(ip))
               dxsitaudf(ip) = xsitau(ip)*(kipinv*dkapdf(ip)
     &              +dphitaudf(ip)/phitau(ip))
               dxsitaudt(ip) = xsitau(ip)*(kipinv*dkapdt(ip)
     &              +dphitaudt(ip)/phitau(ip))
            
               abrad(ip) = abrad(ip)*(1.d0+phitau(ip))/(1.d0+xsitau(ip))
               abdf1(ip) = abrad(ip)*wi(ip)*(kiinv*dkapdf(i)+
     &              piinv*dpdf(i)+dphitaudf(ip)/(1.d0+phitau(ip))
     &              -dxsitaudf(ip)/(1.d0+xsitau(ip)))
               abdf2(ip) = abrad(ip)*wj(ip)*(kipinv*dkapdf(ip)+pipinv*
     &              dpdf(ip)+dphitaudf(ip)/(1.d0+phitau(ip))
     &              -dxsitaudf(ip)/(1.d0+xsitau(ip)))
               abdt1(ip) = abrad(ip)*wi(ip)*(kiinv*dkapdt(i)+
     &              piinv*dpdt(i)-4.d0+dphitaudt(ip)/(1.d0+phitau(ip))
     &              -dxsitaudt(ip)/(1.d0+xsitau(ip)))
               abdt2(ip) = abrad(ip)*wj(ip)*(kipinv*dkapdt(ip)+pipinv*
     &              dpdt(ip)-4.d0+dphitaudt(ip)/(1.d0+phitau(ip))
     &              -dxsitaudt(ip)/(1.d0+xsitau(ip)))
               abdl(ip) = abrad(ip)/lum(ip)*(1.d0-xsitau(ip)/
     &              (1.d0+xsitau(ip)))
               abla(ip) = abrad(ip)
           
            
***   Spherical corrections
*--------------------------

c$$$            if ((ntprof.eq.4.or.ntprof.eq.6.or.ntprof.eq.7)
c$$$     &           .and.abs((dqdtau(ip)+2.d0*(qtau(ip)+tau(ip))
c$$$     &           /(r(ip)*kapm(ip)*rom(ip)))).gt.1.d-6) then
c$$$
c$$$               phitau(ip) = dqdtau(ip)+2.d0*(qtau(ip)+
c$$$     &              tau(ip))/(r(ip)*kapm(ip)*rom(ip))
c$$$               dptau(ip) = -kapm(ip)*lum(ip)/(pim4*c*r(ip)**2)*
c$$$     &              phitau(ip)
c$$$               xsitau(ip) = -dptau(ip)*r(ip)**2/gmrr2
c$$$               dphitaudr(ip) = 
c$$$     &              -2.d0*(qtau(ip)+tau(ip))/(r(ip)*kapm(ip)*rom(ip))
c$$$     &              -2.d0*(1.d0+dqdtau(ip))
c$$$               dphitaudf(ip) =
c$$$     &              +2.d0/(r(ip)*kapm(ip)*rom(ip))*dtaudf(ip)*
c$$$     &              (1.d0+dqdtau(ip))
c$$$     &              -2.d0*(qtau(ip)+tau(ip))/(r(ip)*kapm(ip)*rom(ip))*
c$$$     &              (kipinv*dkapdf(ip)+drodf(ip)/rom(ip))
c$$$               dphitaudt(ip) = 
c$$$     &              +2.d0/(r(ip)*kapm(ip)*rom(ip))*dtaudt(ip)*
c$$$     &              (1.d0+dqdtau(ip))
c$$$     &              -2.d0*(qtau(ip)+tau(ip))/(r(ip)*kapm(ip)*rom(ip))*
c$$$     &              (kipinv*dkapdt(ip)+drodt(ip)/rom(ip))
c$$$               dxsitaudr(ip) = xsitau(ip)*(dphitaudr(ip)/phitau(ip))
c$$$               dxsitaudf(ip) = xsitau(ip)*(kipinv*dkapdf(ip)
c$$$     &              +dphitaudf(ip)/phitau(ip))
c$$$               dxsitaudt(ip) = xsitau(ip)*(kipinv*dkapdt(ip)
c$$$     &              +dphitaudt(ip)/phitau(ip))
c$$$            
c$$$               abrad(ip) = abrad(ip)*(1.d0+phitau(ip))/(1.d0+xsitau(ip))
c$$$               abdf1(ip) = abrad(ip)*wi(ip)*(kiinv*dkapdf(i)+
c$$$     &              piinv*dpdf(i)+dphitaudf(ip)/(1.d0+phitau(ip))
c$$$     &              -dxsitaudf(ip)/(1.d0+xsitau(ip)))
c$$$               abdf2(ip) = abrad(ip)*wj(ip)*(kipinv*dkapdf(ip)+pipinv*
c$$$     &              dpdf(ip)+dphitaudf(ip)/(1.d0+phitau(ip))
c$$$     &              -dxsitaudf(ip)/(1.d0+xsitau(ip)))
c$$$               abdt1(ip) = abrad(ip)*wi(ip)*(kiinv*dkapdt(i)+
c$$$     &              piinv*dpdt(i)-4.d0+dphitaudt(ip)/(1.d0+phitau(ip))
c$$$     &              -dxsitaudt(ip)/(1.d0+xsitau(ip)))
c$$$               abdt2(ip) = abrad(ip)*wj(ip)*(kipinv*dkapdt(ip)+pipinv*
c$$$     &              dpdt(ip)-4.d0+dphitaudt(ip)/(1.d0+phitau(ip))
c$$$     &              -dxsitaudt(ip)/(1.d0+xsitau(ip)))
c$$$               abdl(ip) = abrad(ip)/lum(ip)*(1.d0-xsitau(ip)/
c$$$     &              (1.d0+xsitau(ip)))
c$$$               abla(ip) = abrad(ip)
               
c            xsitau(ip) = lum(ip)/(pim2*gmrr2*c)*(kapm(ip)*dqdtau(ip)*
c     &           0.5d0+(tau(ip)+qtau(ip))/(rom(ip)*r(ip)))
c            abrad(ip) = abrad(ip)*(1.d0+dqdtau(ip)+2.d0*(tau(ip)+
c     &           qtau(ip))/(kapm(ip)*rom(ip)*r(ip)))/(1.d0+xsitau(ip))
            
***   P/T corrections
*--------------------
               
c            cortau1 = ((tau(ip)+qtau(ip))/(tau(ip)+pw23))**0.25d0
c            cortau2 = lum(ip)*(qtau(ip)-pw23)/(pim4*c*r(ip)*r(ip)*
c     &           pm(ip))
c            cortau2 = 1.d0/(1.d0-cortau2)
c     abrad(ip) = abrad(ip)*cortau2/cortau1

***   If not in the atmosphere, classical case
*---------------------------------------------
               
            else
               abdf1(ip) = abrad(ip)*wi(ip)*(kiinv*dkapdf(i)+
     &              piinv*dpdf(i))
               abdf2(ip) = abrad(ip)*wj(ip)*(kipinv*dkapdf(ip)+pipinv*
     &              dpdf(ip))
               abdt1(ip) = abrad(ip)*wi(ip)*(kiinv*dkapdt(i)+
     &              piinv*dpdt(i)-4.d0)
               abdt2(ip) = abrad(ip)*wj(ip)*(kipinv*dkapdt(ip)+pipinv*
     &              dpdt(ip)-4.d0)
               abdl(ip) = abrad(ip)/lum(ip)
            end if
         else
            abla(ip) = abrad(ip)
            abdf1(ip) = abrad(ip)*wi(ip)*(kiinv*dkapdf(i)+piinv*dpdf(i))
            abdf2(ip) = abrad(ip)*wj(ip)*(kipinv*dkapdf(ip)+pipinv*
     &           dpdf(ip))
            abdt1(ip) = abrad(ip)*wi(ip)*(kiinv*dkapdt(i)+piinv*dpdt(i)-
     &           4.d0)
            abdt2(ip) = abrad(ip)*wj(ip)*(kipinv*dkapdt(ip)+pipinv*
     &           dpdt(ip)-4.d0)
            abdl(ip) = abrad(ip)/lum(ip)
         end if

***   Derivatives with respect to u and to lnr
*---------------------------------------------
         
         if (hydro) then
            abdr(ip) = -abrad(ip)/gmrr2*2.d0*r(ip)*accel(ip)
            abdu1(ip) = -abrad(ip)/gmrr2*r(ip)*r(ip)*dtipsi
            abdu2(ip) = -abrad(ip)/gmrr2*r(ip)*r(ip)*(dtninv-dtipsi)
            if (tau(ip).lt.taulim) then
               abdr(ip) = abdr(ip)
     &              + abrad(ip)*(dphitaudr(ip)/(1.d0+phitau(ip))
     &              - dxsitaudr(ip)/(1.d0+xsitau(ip)))
            end if
         else
            if (ntprof.eq.4.or.ntprof.eq.6.or.ntprof.eq.7) then
               if (tau(ip).lt.taulim) then
                  abdr(ip) = abrad(ip)*(dphitaudr(ip)/(1.d0+phitau(ip))
     &                 - dxsitaudr(ip)/(1.d0+xsitau(ip)))
                  abdu1(ip) = 0.d0
                  abdu2(ip) = 0.d0
               else
                  abdr(ip) = 0.d0
                  abdu1(ip) = 0.d0
                  abdu2(ip) = 0.d0
               end if
            else
               abdr(ip) = 0.d0
               abdu1(ip) = 0.d0
               abdu2(ip) = 0.d0
            endif
         endif
      enddo


c$$$         abla(ip) = abrad(ip)
c$$$         if (numeric.eq.2.or.numeric.eq.4) then
c$$$            kiinv = kapm(ip)/kap(i)**2
c$$$            kipinv = kapm(ip)/kap(ip)**2
c$$$         else
c$$$            kiinv = 1.d0/kap(i)
c$$$            kipinv = 1.d0/kap(ip)
c$$$         endif
c$$$
c$$$         abdf1(ip) = abrad(ip)*wi(ip)*(kiinv*dkapdf(i)+piinv*dpdf(i))
c$$$         abdf2(ip) = abrad(ip)*wj(ip)*(kipinv*dkapdf(ip)+pipinv*
c$$$     &        dpdf(ip))
c$$$         abdt1(ip) = abrad(ip)*wi(ip)*(kiinv*dkapdt(i)+piinv*dpdt(i)-
c$$$     &        4.d0)
c$$$         abdt2(ip) = abrad(ip)*wj(ip)*(kipinv*dkapdt(ip)+pipinv*dpdt(ip)
c$$$     &        -4.d0)
c$$$         abdl(ip) = abrad(ip)/lum(ip)
c$$$         if (hydro) then
c$$$            abdr(ip) = -abrad(ip)/gmrr2*2.d0*r(ip)*accel(ip)
c$$$            abdu1(ip) = -abrad(ip)/gmrr2*r(ip)*r(ip)*dtipsi
c$$$            abdu2(ip) = -abrad(ip)/gmrr2*r(ip)*r(ip)*(dtninv-dtipsi)
c$$$         else
c$$$            abdr(ip) = 0.d0
c$$$            abdu1(ip) = 0.d0
c$$$            abdu2(ip) = 0.d0
c$$$         endif


!$OMP END DO

      if (icorr.eq.'a'.or.icorr.eq.'b') then
         do j = nmod-1,2,-1
            if (tau(j).gt.taulim) exit
c               if (abs(mue(j)-mue(j-1)).gt.1.d-1) exit
         enddo
         forall (i=j:nmod)
            egrav(i) = 0.d0
            degdf1(i) = 0.d0
            degdf2(i) = 0.d0
            degdt1(i) = 0.d0
            degdt2(i) = 0.d0
         end forall
      endif

      gmr(1) = g*(4.188787d0*ro(1))**pw23*(0.5d0*(m(1)+m(2)))**pw13
      eshr(1) = eshr(2)
      tconv(nmod) = tconv(nmod1)
      egrav(nmod) = egrav(nmod1)
      abrad(1) = abrad(2)
      abla(1) = abla(2)
      abmu(1) = abmu(2)
      abmuKS(1) = abmuKS(2)
      abled(1) = abled(2)
      hp(1) = hp(2)
      pvisc(nmod) = pvisc(nmod1)
!$OMP END PARALLEL


*____________________________________________________________________
***   calculation of the location and structure of convective regions
***   (also application of overshooting, if required)
*--------------------------------------------------------------------

      lmix = .false.
      if (vlconv) return

      if (.not.((idup.eq.2.or.idup.eq.3).and.iter.gt.abs(itermix))) then
         sconv(1:nmod) = 0.d0
         Dconv(1:nmod) = 0.d0
         tconv(1:nmod) = 1.d99
         call convzone
      endif

      do kl = 1,nsconv
         imin = novlim(kl,3)
         imax = novlim(kl,4)
         if (imin.gt.0) then
            call mlt (imin,imax,error)
            if (error.gt.0) return
         endif
      enddo
      
***   compute parametric overshoot (alpha * Hp)
      if (novopt.ge.1) call oversh
            
*_______________________________________________________________________
***   mixing of convective (+overshoot) regions at each iteration
***   if nmixd.gt.0, diffusion only, no nucleosynthesis
*-----------------------------------------------------------------------

      neutrality = idup.eq.1.or.idup.eq.3
      if (nsconv.gt.0.and..not.vlconv.and.(neutrality.or.(mixopt.and.
     &     (nucreg.eq.1.or.nucreg.eq.2)))) then
         if (iter.le.abs(itermix)) then
c..   restore initial abundances
            do l = 1,nsp
               do k = 1,nmod
                  xsp(k,l) = vxsp(k,l)
                  ysp(k,l) = xsp(k,l)/anuc(l)
               enddo
            enddo
c$$$  !     Tests activation diffusion atomique ou non (modif Td Fev.2019)
            if (microdiffus.and.nphase.gt.5) then
               microdiffus = .false.
               print *, 'phase after 5 microdiffus false'
            else if (nphase.eq.1.and.time*seci.lt.2e7) then ! Ajout test TD Juillet.2019
               microdiffus = .false.
               print *,'no diffusion this young star'
            else if (nphase.eq.1.and.novlim(1,3).eq.1.and.nsconv.eq.1)
     $              then   
               microdiffus = .false.
               print *,'microdiffus false - no radiative core'
               print *, 'novlim(1,3)',novlim(1,3),'nsconv',nsconv
            endif
            
! Fin test activation diffusion atomique 
            if (diffzc.or.nmixd.gt.0) then
               call diffusion (0,error)
               if (error.gt.0) return
            else
               do kl = 1,nsconv
                  imin = novlim(kl,3)
                  imax = novlim(kl,4)
                  call mix (dm,imin,imax,kl,0,partialmix)
               enddo
               if (idiffcc.and..not.microdiffus) then  ! pour gain de temps
                  call diffusion (0,error)
                  if (error.gt.0) return
               endif
            endif
         elseif (iter.eq.abs(itermix)+1.and.itermix.gt.0) then
c..   restore abundances from the previous model
            do l = 1,nsp
               do k = 1,nmod
                  xsp(k,l) = vxsp(k,l)
                  ysp(k,l) = xsp(k,l)/anuc(l)
               enddo
            enddo
         endif
      else
         lmix = .false.
      endif

*__________________________________________________________
***   mixing atmosphere abundances with convective envelope
*----------------------------------------------------------

      if (nsconv.gt.0) then
         if (tau(novlim(nsconv,4)).lt.1.d2) then
            imin = novlim(nsconv,4)
            do i = 1,nsp
               do j = imin,nmod
                  xsp(j,i) = xsp(imin-1,i)
                  ysp(j,i) = xsp(j,i)/anuc(i)
               enddo
            enddo
         endif
      endif

 10   format(2x,'WARNING : pvisc (',1pe10.3,') < 0 #',i4)

      return
      end
