

************************************************************************

      SUBROUTINE mlt (imin,imax,error)

************************************************************************
* Calculate the MLT solution with or without compression               *
* with hp: Kippenhahn analytic formulation                             *
* with hro (with or without turbulence): Pfenniger numerical solution  *
* V2.76 : lambda limited to the distance to the upper conv. boundary   *
************************************************************************
c Modification3 changement de dependance de la convection par rapport au
c temps. Wood 1974 ApJ 190, 609


      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.eng'
      include 'evolcom.eos'
      include 'evolcom.grad'
      include 'evolcom.mod'
      include 'evolcom.opa'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.var'

      integer imin,imax,imin0
      integer error,ibiff
      integer it,iross
      integer i,j,im,ip

      double precision ca,cb,cg,y,sqrz,x4,x
      double precision lambdac,lambdaa,a0
      double precision dfzconv,fzconv,faccor,taumix,betacor
c     double precision dcuarr,dcuau1,dcuau2,dcuaf1,dcuaf2,dcuat1,dcuat2
      double precision delroa,cua,adrad,cc,kiminv,kiinv
      double precision con1,con2,con3,alphac2,ut,ut2,ze,al,
     &     epsf,we,we2,dr1,dr2,dr3,dr4,dr5,dr6,dr7,dutdf1,dutdf2,
     &     dutdt1,dutdt2,dutrr,daldf1,daldf2,daldt1,daldt2,dalrr,daldl,
     &     dalu1,dalu2,wedf1,wedf2,wedt1,wedt2,werr,wedl,weu1,weu2,
     &     zedf1,zedf2,zedt1,zedt2,zerr,zedl,zeu1,zeu2,cua2,tol,z0,
     &     dzadrad,ziter,dznew,znew,ddz,zn,yn,xn4,xn,cu
c     double precision dtni,hpdf1,hpdf2,hpdt1,hpdt2,hpdrr,hpdu1,hpdu2
      double precision xconv,sconv0,xsiconv1,xsiconv
      double precision a1conv,a2conv,a3conv,rconv,qconv,rqconv,
     &     rsign,rqexp,yconv,ztm,csm,abel,a1,A

c      double precision hpdf1,hpdf2,hpdt1,hpdt2,hpdrr,hpdu1,hpdu2,
c     &     c1df1,c1df2,c1dt1,c1dt2,c1drr,c1du1,c1du2,c5df1,c5df2,c5dt1,
c     &     c5dt2,c5drr,c5du1,c5du2,dc1,c3df1,c3df2,c3dt1,c3dt2,c3drr,
c     &     c3du1,c3du2,vdf1,vdf2,vdt1,vdt2,vdrr,vdu1,vdu2,dkvc87
      double precision adraddf1,adraddf2,adraddt1,adraddt2,dkvc8
      double precision alphab,alphat,alphab2,ralpha,ralpha2,ralpha3,dtni

      double precision fconvb(nsh),fconvt(nsh),efft,effb,sconvb,sconvt,
     &     tconvt,tconvb

      double precision om_b,om_t,g0_b,g0_t,c_b,c_t,V_b,yy,c1,B,a00,
     &     pconv,profil
      double precision lambdat(nsh),lambdab(nsh)
      double precision tconv_midCE
 
      dimension lambdac(nsh),lambdaa(nsh),taumix(nsh),abel(nsh)

      common / rossby_number/ tconv_midCE

C.. Modif TD 11/2019 - Ajout Penetration convective Kyle A.       
      double precision zz, omega
      common /rotvar/ omega(nsh)
C.. Fin Modif 

      external pconv

      dtni = 1.d0/dtn
      a0 = 2.25d0 ! form factor for the convective globules
      con1 = 4672.d0/19683.d0
      con2 = 368.d0/729.d0
      con3 = 19.d0/27.d0
      fkcr(1:nmod1) = 0.d0
      alphac2 = alphac*alphac

      if (imin.gt.imax) return
      imin0 = max(imin,2)

*_____________________________
***   MLT prescription with Hp
*-----------------------------

c      print *,'mrBCZ',mr(novlim(nsconv,3)),'mrTCZ',mr(novlim(nsconv,4)),
c    &     'imin0',imin0,'imax',imax,'T(imax)',tm(imax)


      if (.not.ihro) then

***   with analytic (exact) root of the third order polynomial equation
***   following the Cox's formalism (1984)

         if (hpmlt.eq.1.or.hpmlt.eq.2) then
            do 10 i = imin0,imax
c               if (abm(i).gt.abrad(i)) goto 10
               im = i-1
               lambdac(i) = hp(i)*alphac
c               if (r(imax+1).gt.r(i)) lambdac(i) = min(lambdac(i),
c     &              r(imax+1)-r(i))
               A = cpm(i)*kapm(i)*alphac*lambdac(i)*dsqrt(0.5d0*
     &              pm(i)*deltaKSm(i)*rom(i)**3)/(48.d0*sig*tm(i)**3)
               xconv = (A**2/a0*dabs(abrad(i)-abm(i)))**pw13
               a1conv = 1.d0/(a0*xconv)
               a2conv = a1conv/xconv
               a3conv = -1.d0
               rconv = (2.d0*a1conv**3-9.d0*a1conv*a2conv+27.d0*a3conv)/
     &              54.d0
               qconv = (a1conv*a1conv-3.d0*a2conv)/9.d0
               rqconv = rconv*rconv-qconv**3
               rsign = 1.d0
               if (rconv.gt.0.d0) rsign = -1.d0
               rqexp = (dsqrt(rqconv)+dabs(rconv))**pw13
               yconv = rsign*(rqexp+qconv/rqexp)-a1conv/3.d0
               xsiconv = yconv**3
               xsiconv1 = 1.d0-xsiconv
               abla(i) = max(abm(i),xsiconv1*abrad(i)+xsiconv*abm(i))
               eff(i) = abs(xconv*yconv)
               fconv(i) = xsiconv*(abrad(i)-abm(i))/abrad(i)*lum(i)
     &              /(pim4*r(i)**2)
               sconv(i) = sign(alphac*dsqrt(deltaKSm(i)*pm(i)/
     &              (8.d0*rom(i)))*eff(i)/A,fconv(i))
               sconv(i) = dabs(sconv(i))
               abel(i) = (abla(i)+eff(i)*abm(i))/(1.d0+eff(i))

***   other formulations for the flux
c               f1 = 0.5d0*rom(i)*cpm(i)*tm(i)*sconv(i)*alphac*(eff(i)
c     &              /A)**2
c               f2 = 4.d0*rom(i)**2*cpm(i)*tm(i)/(deltaKSm(i)*alphac*
c     &              pm(i))*sconv(i)**3
c               f3 = 0.25d0*cpm(i)*tm(i)*dsqrt(0.5d0*deltaKSm(i)*rom(i)*
c     &              pm(i))*alphac2*(eff(i)/A)**3
c               faccor = 2.25d0*eff(i)**2/(1.d0+eff(i))
c               f4 = (abrad(i)-abm(i))/abrad(i)*lum(i)/(pim4*r(i)**2)*
c     &              faccor/(1.d0+faccor)
c               f6 = (abrad(i)+faccor*abm(i))/(1.d0+faccor)
c               f5 =  (abrad(i)-f6)/abrad(i)*lum(i)/(pim4*r(i)*r(i))
c               write (nout,*) 'flux',i,lum(i)/(pim4*r(i)**2*fconv(i)),f1
c     &           /fconv(i),f2/fconv(i),f3/fconv(i),f4/fconv(i),f5
c     &           /fconv(i)
c               write (nout,*) 'effi',i,(abla(i)-abel(i))/(abel(i)-abm(i)),
c     &              cpm(i) *kapm(i)*rom(i)**2*lambdac(i)*sconv(i)/tm(i)
c     &              **3/(24.d0 *sig),eff(i)

               frad(i) = lum(i)/(pim4*r(i)*r(i))-fconv(i)
               fkin(i) = 0.d0
               fkcr(i) = 0.d0
               Dconv(i) = sconv(i)*lambdac(i)/3.d0
c              tconv(i) = (r(imax)-r(imin))**2/Dconv(i)
               tconv(i) = lambdac(i)/sconv(i)
               taumix(i) = 2.d0*dtn*(sconv(i)+vsconv(i))/lambdac(i)
               if (taumix(i).gt.1.d0.or.hpmlt.eq.1) then
                  abdf1(i) = xsiconv1*abdf1(i)+xsiconv*wi(i)*abm(i)*
     &                 dabadf(im)/abad(im)
                  abdf2(i) = xsiconv1*abdf2(i)+xsiconv*wj(i)*abm(i)*
     &                 dabadf(i)/abad(i)
                  abdt1(i) = xsiconv1*abdt1(i)+xsiconv*wi(i)*abm(i)*
     &                 dabadt(im)/abad(im)
                  abdt2(i) = xsiconv1*abdt2(i)+xsiconv*wj(i)*abm(i)*
     &                 dabadt(i)/abad(i)
                  abdr(i) = xsiconv1*abdr(i)
                  abdl(i) = xsiconv1*abdl(i)
                  abdu1(i) = xsiconv1*abdu1(i)
                  abdu2(i) = xsiconv1*abdu2(i)
               endif
 10         continue
         endif

!     Rotationally modified convection e.g. Augustson & Mathis 2019
         if (hpmlt.eq.6) then
            do 11 i = imin0,imax
               im = i-1
               lambdac(i) = hp(i)*alphac
               A = cpm(i)*kapm(i)*alphac*lambdac(i)*dsqrt(0.5d0*
     &              pm(i)*deltaKSm(i)*rom(i)**3)/(48.d0*sig*tm(i)**3)
               xconv = (A**2/a0*dabs(abrad(i)-abm(i)))**pw13
               a1conv = 1.d0/(a0*xconv)
               a2conv = a1conv/xconv
               a3conv = -1.d0
               rconv = (2.d0*a1conv**3-9.d0*a1conv*a2conv+27.d0*a3conv)/
     &              54.d0
               qconv = (a1conv*a1conv-3.d0*a2conv)/9.d0
               rqconv = rconv*rconv-qconv**3
               rsign = 1.d0
               if (rconv.gt.0.d0) rsign = -1.d0
               rqexp = (dsqrt(rqconv)+dabs(rconv))**pw13
               yconv = rsign*(rqexp+qconv/rqexp)-a1conv/3.d0
               xsiconv = yconv**3
               xsiconv1 = 1.d0-xsiconv
               eff(i) = abs(xconv*yconv)
               fconv(i) = xsiconv*(abrad(i)-abm(i))/abrad(i)*lum(i)
     &              /(pim4*r(i)**2)
               sconv(i) = sign(alphac*dsqrt(deltaKSm(i)*pm(i)/
     &              (8.d0*rom(i)))*eff(i)/A,fconv(i))
               sconv(i) = dabs(sconv(i))

!     Rotational modification
               Call Quintic(omega(i),sconv(i),alphac*hp(i),zz)

               yconv = ((2.5d0)**(1d0/6d0))*yconv/dsqrt(zz)
               xsiconv = yconv**3
               xsiconv1 = 1.d0-xsiconv
               eff(i) = abs(xconv*yconv)
!     It is assumed that the flux is unchanged, but that the gradients
!     do change to compensate for a lower velocity
!     fconv(i) = xsiconv*(abrad(i)-abm(i))/abrad(i)*lum(i)
!     &              /(pim4*r(i)**2)
               sconv(i) = sign(alphac*dsqrt(deltaKSm(i)*pm(i)/
     &              (8.d0*rom(i)))*eff(i)/A,fconv(i))
               sconv(i) = dabs(sconv(i))
               abla(i) = max(abm(i),xsiconv1*abrad(i)+xsiconv*abm(i))
               abel(i) = (abla(i)+eff(i)*abm(i))/(1.d0+eff(i))
               
               frad(i) = lum(i)/(pim4*r(i)*r(i))-fconv(i)
               fkin(i) = 0.d0
               fkcr(i) = 0.d0
               Dconv(i) = sconv(i)*lambdac(i)/3.d0
               tconv(i) = lambdac(i)/sconv(i)
               taumix(i) = 2.d0*dtn*(sconv(i)+vsconv(i))/lambdac(i)
               
               if (taumix(i).gt.1.d0) then
                  abdf1(i) = xsiconv1*abdf1(i)+xsiconv*wi(i)*abm(i)*
     &                 dabadf(im)/abad(im)
                  abdf2(i) = xsiconv1*abdf2(i)+xsiconv*wj(i)*abm(i)*
     &                 dabadf(i)/abad(i)
                  abdt1(i) = xsiconv1*abdt1(i)+xsiconv*wi(i)*abm(i)*
     &                 dabadt(im)/abad(im)
                  abdt2(i) = xsiconv1*abdt2(i)+xsiconv*wj(i)*abm(i)*
     &                 dabadt(i)/abad(i)
                  abdr(i) = xsiconv1*abdr(i)
                  abdl(i) = xsiconv1*abdl(i)
                  abdu1(i) = xsiconv1*abdu1(i)
                  abdu2(i) = xsiconv1*abdu2(i)
               endif
 11         continue
         endif
         
***   with analytic (approximate) root of the third order polynomial
***   equation following Kippenhahn's formalism (1991)

         if (hpmlt.eq.3.or.hpmlt.eq.4) then
            do 20 i = imin0,imax
               if (abm(i).gt.abrad(i)) goto 20
               im = i-1
               lambdac(i) = hp(i)*alphac
c               if (lambdac(i).gt.(r(imax+1)-r(i))) then
c                  lambdac(i) = min(lambdac(i),r(imax+1)-r(i))
c               endif
               ut = 24.d0*sig*tm(i)**3/(kapm(i)*cpm(i)*lambdac(i)*
     &              alphac*dsqrt(0.5d0*deltaKSm(i)*pm(i)*rom(i)**3))
               if (ut.lt.1.d-12) then
                  ut = 0.d0
                  ze = 0.d0
                  zedf1 = 0.d0
                  zedf2 = 0.d0
                  zedt1 = 0.d0
                  zedt2 = 0.d0
                  dutdf1 = 0.d0
                  dutdf2 = 0.d0
                  dutdt1 = 0.d0
                  dutdt2 = 0.d0
                  zerr = 0.d0
                  dutrr = 0.d0
                  zedl = 0.d0
                  zeu1 = 0.d0
                  zeu2 = 0.d0
                  abla(i) = abm(i)
                  fconv(i) = (abrad(i)-abla(i))/abrad(i)*lum(i)/(pim4*
     &                 r(i)*r(i))
                  sconv(i) = sign((0.25d0*alphac*dabs(fconv(i))*pm(i)*
     &                 deltaKSm(i)/(cpm(i)*rom(i)**2*tm(i)))**pw13,
     &                 fconv(i))
                  sconv(i) = dabs(sconv(i))
                  eff(i) = cpm(i)*kapm(i)*rom(i)**2*lambdac(i)*sconv(i)/
     &                 (24.d0*sig*tm(i)**3)
                  abel(i) = (abla(i)+eff(i)*abm(i))/(1.d0+eff(i))
                  frad(i) = lum(i)/(pim4*r(i)*r(i))-fconv(i)
                  fkin(i) = 0.d0
                  fkcr(i) = 0.d0
                  Dconv(i) = sconv(i)*lambdac(i)/3.d0
c                 tconv(i) = (r(imax)-r(imin))**2/Dconv(i)
                  tconv(i) = lambdac(i)/sconv(i)
                  goto 25
               endif
               ut2 = ut*ut
               al = ((abrad(i)-abm(i))/a0+con1*ut2)*ut
               epsf = dsqrt(al*al+(con2*ut2)**3)
               we = (dabs(al+epsf))**pw13
               ze = we+con3*ut-con2*ut2/we
               A1 = alphac*dsqrt(deltaKSm(i)*pm(i)/(8.d0*rom(i)))
               sconv(i) = A1*(ze-ut)+1.d-10
               sconv(i) = dabs(sconv(i))
c..   impose convective cells to be sub-sonic
c               if (sconv(i).gt.0.8d0*cs(i)) then
c                  sconv(i) = 0.8d0*cs(i)
c                  ze = sconv(i)/A1+ut
c               endif
               fconv(i) = 0.25d0*cpm(i)*tm(i)*dsqrt(0.5d0*deltaKSm(i)*
     &              rom(i)*pm(i))*alphac2*(ze-ut)**3
               eff(i) = 0.5d0*(ze/ut-1.d0)
               abla(i) = max(abm(i),abm(i)+ze*ze-ut2)

               abel(i) = abla(i)-(ze-ut)**2
               frad(i) = lum(i)/(pim4*r(i)*r(i))-fconv(i)
               fkin(i) = 0.d0
               fkcr(i) = 0.d0
               Dconv(i) = sconv(i)*lambdac(i)/3.d0
c              tconv(i) = (r(imax)-r(imin))**2/Dconv(i)
               tconv(i) = lambdac(i)/sconv(i)
               we2 = we*we
               dr1 = al/ut+2.d0*con1*ut2
               dr2 = 1.d0+con2*ut2/we2
               dr3 = con3-2.d0*con2*ut/we
               dr4 = 1.d0/(3.d0*we2)
               dr5 = 3.d0*con2**3*ut**5
               dr6 = ut/a0
               dr7 = (1.d0+al/epsf)*dr4
               if (numeric.eq.2.or.numeric.eq.4) then
                  kiminv = kapm(i)/kap(im)**2
                  kiinv = kapm(i)/kap(i)**2
               else
                  kiminv = 1.d0/kap(im)
                  kiinv = 1.d0/kap(i)
               endif
               dutdf1 = -wi(i)*ut*(1.5d0*dpdf(im)/p(im)+dcpdf(im)/
     &              cp(im)+dkapdf(im)*kiminv+0.5d0*
     &              deltaKSdf(im)/deltaKS(im)+0.5d0*drodf(im)/ro(im))
               dutdf2 = -wj(i)*ut*(1.5d0*dpdf(i)/p(i)+dcpdf(i)/cp(i)+
     &              dkapdf(i)*kiinv+0.5d0*deltaKSdf(i)/deltaKS(i)+
     &              0.5d0*drodf(i)/ro(i))
               dutdt1 = wi(i)*ut*(3.d0-dcpdt(im)/cp(im)-
     &              dkapdt(im)*kiminv-1.5d0*dpdt(im)/p(im)-0.5d0*
     &              deltaKSdt(im)/deltaKS(im)-0.5d0*drodt(im)/ro(im))
               dutdt2 = wj(i)*ut*(3.d0-dcpdt(i)/cp(i)-
     &              dkapdt(i)*kiinv-1.5d0*dpdt(i)/p(i)-
     &              0.5d0*deltaKSdt(i)/deltaKS(i)-0.5d0*drodt(i)/ro(i))
               dutrr = -2.d0*ut/r(i)
               daldf1 = dr6*(abdf1(i)-wi(i)*abm(i)*dabadf(im)/abad(im))+
     &              dr1*dutdf1
               daldf2 = dr6*(abdf2(i)-wj(i)*abm(i)*dabadf(i)/abad(i))+
     &              dr1*dutdf2
               daldt1 = dr6*(abdt1(i)-wi(i)*abm(i)*dabadt(im)/abad(im))+
     &              dr1*dutdt1
               daldt2 = dr6*(abdt2(i)-wj(i)*abm(i)*dabadt(i)/abad(i))+
     &              dr1*dutdt2
               dalrr = dr6*abdr(i)+dr1*dutrr
               daldl = dr6*abdl(i)
               dalu1 = dr6*abdu1(i)
               dalu2 = dr6*abdu2(i)
               wedf1 = dr4*(daldf1+(al*daldf1+dr5*dutdf1)/epsf)
               wedf2 = dr4*(daldf2+(al*daldf2+dr5*dutdf2)/epsf)
               wedt1 = dr4*(daldt1+(al*daldt1+dr5*dutdt1)/epsf)
               wedt2 = dr4*(daldt2+(al*daldt2+dr5*dutdt2)/epsf)
               werr = dr4*(dalrr+(al*dalrr+dr5*dutrr)/epsf)
               wedl = dr7*daldl
               weu1 = dr7*dalu1
               weu2 = dr7*dalu2
               zedf1 = wedf1*dr2+dutdf1*dr3
               zedf2 = wedf2*dr2+dutdf2*dr3
               zedt1 = wedt1*dr2+dutdt1*dr3
               zedt2 = wedt2*dr2+dutdt2*dr3
               zerr = werr*dr2+dutrr*dr3
               zedl = wedl*dr2
               zeu1 = weu1*dr2
               zeu2 = weu2*dr2
 25            taumix(i) = 2.d0*dtn*(sconv(i)+vsconv(i))/lambdac(i)
               if (taumix(i).gt.1.d0.or.hpmlt.eq.3) then
                  abdf1(i) = wi(i)*abm(i)*dabadf(im)/abad(im)+
     &                 2.d0*(ze*zedf1-ut*dutdf1)
                  abdf2(i) = wj(i)*abm(i)*dabadf(i)/abad(i)+
     &                 2.d0*(ze*zedf2-ut*dutdf2)
                  abdt1(i) = wi(i)*abm(i)*dabadt(im)/abad(im)+
     &                 2.d0*(ze*zedt1-ut*dutdt1)
                  abdt2(i) = wj(i)*abm(i)*dabadt(i)/abad(i)+
     &                 2.d0*(ze*zedt2-ut*dutdt2)
                  abdr(i) = 2.d0*(ze*zerr-ut*dutrr)
                  abdl(i) = 2.d0*ze*zedl
                  abdu1(i) = 2.d0*ze*zeu1
                  abdu2(i) = 2.d0*ze*zeu2
               endif
 20         continue
         endif

*____________________________________________________________
***   convection model with compression effects (and with hp)
*     Forestini, Lumer, Arnould 1991, A&A 252, 127
*------------------------------------------------------------

         if (hpmlt.eq.5) then

            alphat = alphac
            alphab = alphatu
            alphab2 = alphab*alphab
            ralpha = alphab/alphat
            ralpha2 = ralpha*ralpha
            ralpha3 = ralpha*ralpha2
            a00 = etaturb

            do i = imin0,imax
	       im = i - 1
               adrad = dabs(abrad(i)-abm(i))
               lambdat(i) = hp(i)*alphat
               lambdab(i) = hp(i)*alphab
               lambdac(i) = 0.5d0 * (lambdat(i) + lambdab(i))

               om_b = kapm(i)*lambdab(i)*rom(i)
               om_t = kapm(i)*lambdat(i)*rom(i)
               c1 = cpm(i)*rom(i)/(8.d0*sig*tm(i)**3)
               g0_b = c1*(1.d0+pw13*om_b**2)/om_b
               g0_t = c1*(1.d0+pw13*om_t**2)/om_t

               c_b = alphab2*gmr(i)*hp(i)*deltaKSm(i)/(8.d0*a00)
               c_t = c_b/ralpha2
               V_b = 1.d0/(g0_b*sqrt(c_b*adrad))
               B = 0.75d0*om_b*c1/g0_b
               B = a0

               profil = pconv (ralpha,V_b,B)
               yy = V_b*(ralpha3*profil-1.d0)/(1.d0-ralpha2*profil**2)

               sconvb = yy/(g0_b*V_b) 
               sconvt = profil*sconvb
               abla(i) = max(abm(i),abm(i) + adrad*yy*(yy+V_b))
c               abelb = abla(i) - sconvb*sconvb/c_b
c               abelt = abla(i) - sconvt*sconvt/c_t
               effb = g0_b*sconvb
               efft = g0_t*sconvt

               fconvb(i) = 0.5d0*profil*cpm(i)*rom(i)*alphab*sconvb**3*
     &              tm(i)/(c_b*(1.d0+profil))

               fconvt(i) = 0.5d0*cpm(i)*rom(i)*alphat*sconvt**3*
     &              tm(i)/(c_t*(1.d0+profil))
               fconv(i) = fconvb(i)+fconvt(i)
               fkin(i) = rom(i)*sconvb**3*profil*(profil-1.d0)/2.d0
               frad(i) = lum(i)/(pim4*r(i)*r(i)) - fconv(i)
               fkcr(i) = fkin(i)/fconv(i)
               !frad(i) = lum(i)/(pim4*r(i)*r(i)) - fconv(i) -fkin(i)
c               write (*,'(i4,3(1x,1pe13.6))') i,fconv(i),fkin(i),fkcr(i)

               tconvb = lambdab(i)/sconvb
               tconvt = lambdat(i)/sconvt

               sconv(i) = 0.5d0*(sconvb+sconvt)
               eff(i) = 0.5d0*(efft+effb)
               tconv(i) = (tconvt+tconvb)*0.5d0
               taumix(i) = 2.d0*dtn*(sconv(i)+vsconv(i))/lambdac(i)
               Dconv(i) = (sconvb*lambdab(i)+sconvt*lambdat(i))/6.d0

c                hpdf1 = wi(i)*(dpdf(im)/pm(i)-drodf(im)/rom(i))
c                hpdf2 = wj(i)*(dpdf(i)/pm(i)-drodf(i)/rom(i))
c                hpdt1 = wi(i)*(dpdt(im)/pm(i)-drodt(im)/rom(i))
c                hpdt2 = wj(i)*(dpdt(i)/pm(i)-drodt(i)/rom(i))
c                hpdrr = 2.d0*gmr(i)/(gmr(i)+accel(i))
c                hpdu1 = -dtni*dynfac*psi(i)/(gmr(i)+accel(i))
c                hpdu2 = -dtni*dynfac*(1.d0-psi(i))/(gmr(i)+accel(i))

c                if (numeric.eq.2.or.numeric.eq.4) then
c                   kiminv = kapm(i)/kap(im)**2
c                   kiinv = kapm(i)/kap(i)**2
c                else
c                   kiminv = 1.d0/kap(im)
c                   kiinv = 1.d0/kap(i)
c                endif
c                c1df1 = wi(i)*(dkapdf(im)*kiminv+drodf(im)/rom(i))+hpdf1 
c                c1df2 = wj(i)*(dkapdf(i) *kiinv +drodf(i) /rom(i))+hpdf2     
c                c1dt1 = wi(i)*(dkapdt(im)*kiminv+drodt(im)/rom(i))+hpdt1
c                c1dt2 = wj(i)*(dkapdt(i) *kiinv +drodt(i) /rom(i))+hpdt2 
c                c1drr = hpdrr
c                c1du1 = hpdu1
c                c1du2 = hpdu2

c                c5df1 = wi(i)*deltaKSdf(im)/deltaKSm(i)+hpdf1
c                c5df2 = wj(i)*deltaKSdf(i) /deltaKSm(i)+hpdf2
c                c5dt1 = wi(i)*deltaKSdt(im)/deltaKSm(i)+hpdt1
c                c5dt2 = wj(i)*deltaKSdt(i) /deltaKSm(i)+hpdt2
c                c5drr = hpdrr-2.d0
c                c5du1 = hpdu1
c                c5du2 = hpdu2

c                dc1 = (om_b*om_b-3.d0)/(3.d0+om_b*om_b)
c                c3df1 = wi(i)*(dcpdf(im)/cpm(i)+drodt(im)/rom(i))
c      &              + dc1*c1df1/om_b
c                c3df2 = wj(i)*(dcpdf(i)/cpm(i)+drodt(i)/rom(i))
c      &              + dc1*c1df2/om_b
c                c3dt1 = wi(i)*(dcpdt(im)/cpm(i)+drodt(im)/rom(i)-3.d0)
c      &              + dc1*c1dt1/om_b
c                c3dt2 = wj(i)*(dcpdt(i)/cpm(i) + drodt(i)/rom(i) - 3.d0)
c      &              + dc1*c1dt2/om_b
c                c3drr = dc1*c1drr/om_b
c                c3du1 = dc1*c1du1*dynfac
c                c3du2 = dc1*c1du2*dynfac

               adraddf1 = abdf1(i)-wi(i)*abm(i)*dabadf(im)/abad(im)
               adraddf2 = abdf2(i)-wj(i)*abm(i)*dabadf(i)/abad(i)
               adraddt1 = abdt1(i)-wi(i)*abm(i)*dabadt(im)/abad(im)
               adraddt2 = abdt2(i)-wj(i)*abm(i)*dabadt(i)/abad(i)

c                vdf1 = -c3df1-0.5d0*(c5df1+adraddf1/adrad)
c                vdf2 = -c3df2-0.5d0*(c5df2+adraddf2/adrad)
c                vdt1 = -c3dt1-0.5d0*(c5dt1+adraddt1/adrad)
c                vdt2 = -c3dt2-0.5d0*(c5dt2+adraddt2/adrad)
c                vdrr = -c3drr-0.5d0*(c5drr+abdr(i)/adrad)
c                vdu1 = -c3du1-0.5d0*(c5du1+abdu1(i)/adrad)
c                vdu2 = -c3du2-0.5d0*(c5du2+abdu2(i)/adrad)

                dkvc8 = yy*(yy+V_b)
c                dkvc87 = adrad*yy*V_b

               abdf1(i) = wi(i)*abm(i)*dabadf(im)/abad(im)+dkvc8*
     &              adraddf1 !+dkvc87*vdf1
               abdf2(i) = wj(i)*abm(i)*dabadf(i)/abad(i)+dkvc8*adraddf2
     &              !+dkvc87*vdf2
               abdt1(i) = wi(i)*abm(i)*dabadt(im)/abad(im)+dkvc8*
     &              adraddt1 !+dkvc87*vdt1
               abdt2(i) = wj(i)*abm(i)*dabadt(i)/abad(i)+dkvc8*adraddt2
     &              !+dkvc87*vdt2
               abdl(i) = dkvc8*abdl(i)
               abdr(i) = dkvc8*abdr(i) !+dkvc87*vdrr
               abdu1(i) = dkvc8*abdu1(i) !+dkvc87*vdu1
               abdu2(i) = dkvc8*abdu2(i) !+dkvc87*vdu2
            enddo

         endif


         if (hpmlt.eq.0) then
            do i = imin0,imax
               fconv(i) = 0.d0
               frad(i) = lum(i)/(pim4*r(i)*r(i))-fconv(i)
               fkin(i) = 0.d0
               fkcr(i) = 0.d0
               Dconv(i) = 0.d0
               tconv(i) = 1.d99
            enddo
         endif
         if (t(imax).lt.8.d3) then
            do i = imax,imin,-1
               if (t(i).le.8.d3) exit
            enddo
            do j = i,imax
               Dconv(j) = Dconv(i)
               sconv(j) = 3.d0*Dconv(j)/lambdac(j)
               tconv(j) = lambdac(j)/sconv(j)
            enddo
         endif

      else

*_______________________________________________________________
***   MLT prescription with hro (and with or without turbulence)
*---------------------------------------------------------------

         do i = imin0,imax
            im = i-1
            if (zt(i)*zt(i-1).lt.0.d0) write (nout,*) 'WARNING : ',
     &           'zt < 0 in MLT, shell',i
            ztm = dabs(zt(i-1))**wi(i)*dabs(zt(i))**wj(i)
            csm = dabs(cs(i-1))**wi(i)*dabs(cs(i))**wj(i)
            delroa = 1.d0/ztm-abm(i)
            adrad = abrad(i)-abm(i)
            lambdaa(i) = alphatu*hp(i)/(deltaKSm(i)*delroa)
c            if (r(imax+1).gt.r(i))  lambdaa(i) = min(lambdaa(i),
c     &           r(imax+1)-r(i))
            cua = 12.d0*sig*tm(i)**3/(cpm(i)*rom(i)*rom(i)*kapm(i)*
     &           lambdaa(i)*lambdaa(i))*dsqrt(8.d0*hp(i)/(gmr(i)*
     &           deltaKSm(i)))
            cua2 = cua*cua
            cc = etaturb/9.d0*(12.d0*sig*tm(i)**3/(cpm(i)*rom(i)*
     &           rom(i)*kapm(i)))**3*(gmr(i)*deltaKSm(i))**2/(cpm(i)*
     &           tm(i)*hp(i)*csm**5)
            ca = delroa/(cua*cua)
            cb = 9.d0/16.d0*cua*cua/adrad
            cg = adrad*cc

            it = 0
            tol = 1.d-8
            z0 = 1.d0
 30         z0 = z0*5.d0
            y = z0**3+8.d0/9.d0*z0*(z0+2.d0)
            sqrz = dsqrt(y*y+cg*z0**8)
            x4 = cb*(y+sqrz)
            x = x4**0.25d0
            dzadrad = -2.d0*(z0+1.d0)+ca*cb*(1.d0-3.d0/(4.d0*x))/sqrz*
     &           (4.d0*cg*z0**7+(y+sqrz)*(3.d0*z0*z0+16.d0/9.d0*
     &           (z0+1.d0)))
            if (dzadrad.gt.0.d0) goto 40
            goto 30
 40         ziter = z0
 50         it = it+1
            y = ziter**3+8.d0/9.d0*ziter*(ziter+2.d0)
            sqrz = dsqrt(y*y+cg*ziter**8)
            x4 = cb*(y+sqrz)
            x = x4**0.25d0
            fzconv = ca*x**3*(x-1.d0)-ziter*(ziter+2.d0)
            dfzconv = -2.d0*(ziter+1.d0)+ca*cb*(1.d0-3.d0/(4.d0*x))/
     &           sqrz*(4.d0*cg*ziter**7+(y+sqrz)*(3.d0*ziter*ziter+
     &           16.d0/9.d0*(ziter+1.d0)))
            dznew = fzconv/dfzconv
            znew = ziter-dznew
            ddz = dabs(dznew)/znew
            if (ddz.lt.tol.and.znew.gt.0.d0) goto 60
            if (it.gt.50) then
               write (90,1000) i,ziter,znew
               error = 2
               return
            endif
            ziter = znew
            goto 50
 60         zn = znew
            yn = zn**3+8.d0/9.d0*zn*(zn+2.d0)
            xn4 = 9.d0/16.d0*cua2/adrad*(yn+dsqrt(yn*yn+adrad*cc*zn**8))
            xn = xn4**0.25d0
            cu = cua/(xn*xn)
            lambdac(i) = xn*lambdaa(i)
            abla(i) = max(abm(i),abm(i)+cu*cu*zn*(zn+2.d0))
            eff(i) = zn*0.5d0
            sconv(i) = dsqrt(gmr(i)*lambdac(i)*lambdac(i)*deltaKSm(i)/
     &           (8.d0*hp(i)))*cu*zn
            sconv(i) = dabs(sconv(i))
            abel(i) = (abla(i)+eff(i)*abm(i))/(eff(i)+1.d0)
            fconv(i) = lambdac(i)*cpm(i)*rom(i)*tm(i)*sconv(i)/
     &           (2.d0*hp(i))*(abla(i)-abel(i))
            fkin(i) = etaturb*rom(i)*sconv(i)**8/csm**5
            fkcr(i) = fkin(i)/fconv(i)
            frad(i) = lum(i)/(pim4*r(i)*r(i))-fconv(i)-fkin(i)
            Dconv(i) = sconv(i)*lambdac(i)/3.d0
c     tconv(i) = (r(imax)-r(imin))**2/Dconv(i)
            tconv(i) = lambdac(i)/sconv(i)
c      hpdf1 = wi(i)*(dpdf(im)/p(im)-drodf(im)/ro(im))
c      hpdf2 = wj(i)*(dpdf(i)/p(i)-drodf(i)/ro(i))
c      hpdt1 = wi(i)*(dpdt(im)/p(im)-drodt(im)/ro(im))
c      hpdt2 = wj(i)*(dpdt(i)/p(i)-drodt(i)/ro(i))
c      hpdrr = 2.d0*gmr(i)/(gmr(i)+accel(i))
c      hpdu1 = -dtni*psi(i)/(gmr(i)+accel(i))*dynfac
c      hpdu2 = -dtni*(1.d0-psi(i))/(gmr(i)+accel(i))*dynfac
c      dcuarr = -3.d0/hp(i)*hpdrr+2.d0/(gmr(i)*r(i))
c      dcuau1 = -3.d0/hp(i)*hpdu1
c      dcuau2 = -3.d0/hp(i)*hpdu2
c      if (numeric.eq.2.or.numeric.eq.4) then
c         kiminv = kapm(i)/kap(im)**2
c         kiinv = kapm(i)/kap(i)**2
c      else
c         kiminv = 1.d0/kap(im)
c         kiinv = 1.d0/kap(i)
c      endif
c      dcuaf1 = -4.d0/delroa/(zt(im)*zt(im))*dztdf(im)-2.d0/delroa*
c     &           wi(i)*abm(i)*dabadf(im)/abad(im)+1.5d0/deltaKS(im)*
c     &           deltaKSdf(im)-dcpdf(im)/cp(im)-dkapdf(im)*kininv-
c     &           2.d0*drodf(im)/ro(im)
c      dcuaf2 = -2.d0/delroa*wj(i)*abm(i)*dabadf(i)/abad(i)+
c     &           1.5d0/deltaKS(i)*deltaKSdf(i)-dcpdf(i)/cp(i)-
c     &           dkapdf(i)*kiinv-2.d0*drodf(i)/ro(i)
c      dcuat1 = -4.d0/delroa/(zt(i)*zt(i))*dztdt(i)-2.d0/delroa*
c     &           wi(i)*abm(i)*dabadt(im)/abad(im)+1.5d0/deltaKS(im)*
c     &           deltaKSdt(im)-dcpdt(im)/cp(im)-dkapdt(im)*kininv+3.d0-
c     &           2.d0*drodt(im)/ro(im)
c      dcuat2 = -2.d0/delroa*wj(i)*abm(i)*dabadt(i)/abad(i)+
c     &           1.5d0/deltaKS(i)*deltaKSdt(i)-dcpdt(i)/cp(i)-
c     &           dkapdt(i)*kiinv+3.d0-2.d0*drodt(i)/ro(i)
c      abdf1(i) = wi(i)*abm(i)*dabadf(im)/abad(im)+cua2*dcuaf1
c      abdf2(i) = wj(i)*abm(i)*dabadf(i)/abad(i)+cua2*dcuaf2
c      abdt1(i) = wi(i)*abm(i)*dabadt(im)/abad(im)+cua2*dcuat1
c      abdt2(i) = wj(i)*abm(i)*dabadt(i)/abad(i)+cua2*dcuat2
c      abdr(i) = cua2*dcuarr
c      abdl(i) = 0.d0
c      abdu1(i) = cua2*dcuau1
c      abdu2(i) = cua2*dcuau2            
            abdf1(i) = wi(i)*abm(i)*dabadf(im)/abad(im)
            abdf2(i) = wj(i)*abm(i)*dabadf(i)/abad(i)
            abdt1(i) = wi(i)*abm(i)*dabadt(im)/abad(im)
            abdt2(i) = wj(i)*abm(i)*dabadt(i)/abad(i)
            abdr(i) = 0.d0
            abdl(i) = 0.d0
            abdu1(i) = 0.d0
            abdu2(i) = 0.d0
         enddo
      endif

      if (hpmlt.le.1.or.hpmlt.eq.3.or.hpmlt.eq.5) goto 80

*____________________________________________________
***   treatment of time-dependent convection
***   Sparks, Starrfield, Truran, 1978, ApJ 220, 1063
***   Wood 1974 ApJ 190, 609
*----------------------------------------------------

      ibiff = 0
      do i = imin0,imax
         if (taumix(i).le.1.d0) then
            if (taumix(i).lt.0.d0) then
               ibiff = ibiff + 1
            endif
            error = -9
            sconv0 = sconv(i)
            sconv(i) = vsconv(i)+taumix(i)*(sconv(i)-vsconv(i))
            im = max(i-1,imin0)
            ip = min(i+1,imax)
            con1 = max(0.d0,pw13*(1.d0-(vr(i)-vr(im))/lambdac(i)))
            con2 = max(0.d0,pw13*(1.d0-(vr(ip)-vr(i))/lambdac(i)))
            sconv(i) = con1*vsconv(im)+(1.d0-con1-con2)*sconv(i)+con2*
     &           vsconv(ip)

            eff(i) = eff(i)*sconv(i)/sconv0
            faccor = 2.25d0*eff(i)**2/(1.d0+eff(i))
            fconv(i) = 4.d0*rom(i)**2*cpm(i)*tm(i)/(deltaKSm(i)*alphac*
     &           pm(i))*sconv(i)**3
            abla(i) = (abrad(i)+faccor*abm(i))/(1.d0+faccor)

c            A = cpm(i)*kapm(i)*alphac*lambdac(i)*dsqrt(0.5d0*
c     &           pm(i)*deltaKSm(i)*rom(i)**3)/(48.d0*sig*tm(i)**3)
c            f2 = 0.5d0*rom(i)*cpm(i)*tm(i)*sconv(i)*alphac*(eff(i)
c     &           /A)**2
c            f3 = sign(0.25d0*cpm(i)*tm(i)*dsqrt(0.5d0*deltaKSm(i)*
c     &           rom(i)*pm(i))*alphac2*(eff(i)/A)**3,fconv(i))
c            faccor = 2.25d0*eff(i)**2/(1.d0+eff(i))
c            f4 = (abrad(i)-abm(i))/abrad(i)*lum(i)/(pim4*r(i)**2)*
c     &           faccor/(1.d0+faccor)
c            f5 =  (abrad(i)-abla(i))/abrad(i)*lum(i)/(pim4*r(i)*r(i))
c            faccor = eff(i)/(1.d0+eff(i))
c            f6 = 0.5d0*rom(i)*cpm(i)*tm(i)*sconv(i)*alphac*(abla(i)
c     &           -abm(i))*faccor
c            write (nout,*) i,f1/fconv(i),f2/fconv(i),f3/fconv(i),
c                 f4/fconv(i),f5/fconv(i),f6/fconv(i),sconv(i)/sconv0

            betacor = 1.d0/(1.d0+faccor)
            abdr(i) = abdr(i)*betacor
            abdl(i) = abdl(i)*betacor
            abdu1(i) = abdu1(i)*betacor
            abdu2(i) = abdu2(i)*betacor
            abdf1(i) = betacor*(abdf1(i)+faccor*wi(i)*abm(i)*
     &           dabadf(im)/abad(im))
            abdf2(i) = betacor*(abdf2(i)+faccor*wj(i)*abm(i)*
     &           dabadf(i)/abad(i))
            abdt1(i) = betacor*(abdt1(i)+faccor*wi(i)*abm(i)*
     &           dabadt(im)/abad(im))
            abdt2(i) = betacor*(abdt2(i)+faccor*wj(i)*abm(i)*
     &           dabadt(i)/abad(i))
         endif
      enddo

      if (ibiff.ge.1) then
         write (nout,*) 'WARNING: ',ibiff,' shell(s) with taumix < 0 !'
      endif

 80   if (imin.eq.1) then
         fconv(1) = 0.d0
         frad(1) = 0.d0
         fkin(1) = 0.d0
         sconv(1) = sconv(2)
         abel(1) = abel(2)
         Dconv(1) = Dconv(2)
         eff(1) = eff(2)
         abla(1) = abla(2)
         taumix(1) = 2.d0*dtn*(sconv(1)+vsconv(1))/(hp(1)*alphac)
         if (sconv(1).ne.0.d0) then
            tconv(1) = hp(1)*alphac/dabs(sconv(1))
	 else
            tconv(1) = 1.d99 
         endif
      endif
c..   to prevent "abnormal" convective timescales at the upper boundary
      tconv(imax) = tconv(imax-1)
      tconv_midCE = 0.d0
c...  Rossby number
      iross = imin0
      do i = imin0,imax-1
         if (hp(i+1).le.0.1d0*hp(imin).and.hp(i).gt.0.1d0*hp(imin))
     &        iross = i+1
      enddo
c      print *,'imin0,imax,iross,imin0',imax,iross
      if (abs(iross).le.imax) tconv_midCE = tconv(iross)
   
 1000 format (5x,'50 iterations reached for MLT with hro, at shell ',
     &     i4,', with ziter = ',1pe10.4,' and znew = ',1pe10.4)

      return
      end
