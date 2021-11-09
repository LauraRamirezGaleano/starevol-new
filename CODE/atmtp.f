      SUBROUTINE atmtp

************************************************************************
* Calculate the atmospheric temperature profile and derivatives        *
*                                                                      *
* $LastChangedDate:: 2018-01-22 10:12:58 +0000 (Mon, 22 Jan 2018)    $ *
* $Author:: amard                                                    $ *
* $Rev:: 110                                                         $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.atm'
      include 'evolcom.mod'
      include 'evolcom.mass'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.var'
      include 'evolcom.therm'

      include 'evolcom.grad'

      integer i,j,k
      integer hs,hhs

      double precision exv
      double precision loggs,ts,gs,gs1,ms,ms1
      double precision fa1(4),fa2(4),fa3(8),fa4(4),fa5(8)
      double precision a1,a2,a3,a4,a5
      double precision expq1,expq2
      double precision alpha20,beta0,gamma0,delta0,pwa1,pwa2,pwa3,pwb1,
     &     pwb2,pwb3,pwc1,pwc2,pwc3,pwd1,pwd2,pwd3
      double precision alphap,betap,gammap,deltap
      double precision qq1,qq2
      double precision Tatm,sig_av(nmod)
      double precision qtauKS(nmod),dqdtauKS(nmod),qtauLS(nmod)
     $     ,dqdtauLS(nmod),ddqtauKS(nmod),ddqtauLS(nmod),qtaubis(nmod)
      double precision Tatmd,rhoatmd,Ratmd
      double precision rhoat,valmax,taumin,taumax,taulim0,taue
      double precision rhoatm(nmod),convatm(nmod)
      double precision rtau,rT,rhotab,Fconvtab
      double precision Ttabeff,gtabeff,Ztabeff
      double precision FeH,taueffb
      double precision dqtautemp(nmod),ddqtautemp(nmod)
      double precision Tatms(nmod),taumod(nmod)

      integer filesize

      logical ntp,ntp2

      external exv

      data (fa1(j),j = 1,4) / 3.497886d0,-6.24d-4,-2.5416d-2,9.095d-3 /
      data (fa2(j),j = 1,4) / -1.837334d0,3.78d-4,1.7968d-2,-8.235d-3 /
      data (fa3(j),j = 1,8) / 2.817937d0,3.27984d-1,3.264645d0,
     &     -2.704504d0,-1.576d-3,-1.08d-4,-1.156d-3,9.22d-4 /
      data (fa4(j),j = 1,4) / -1.030175d0,2.26d-4,4.521d-3,2.103d-3 /
      data (fa5(j),j = 1,8) / 6.382663d0,4.130736d0,-1.7722d-2,
     &     -1.316d-3,8.359776d0,1.30637215d2,-3.04d-3,-4.2022d-2 /

      common /interpatm/ rhoat(nsh),Tatm(nsh),ntp2
!      common /atmospheres/ rtau(127,34,10),rT(127,34,10),
!     &     rhotab(127,34,10),Fconvtab(127,34,10),Ttabeff(34),
!     &     gtabeff(10),nb_geff,nb_Teff,filesize
      common /atmospheres/ rtau(127,10,59,6),rT(127,10,59,6),
     &     rhotab(127,10,59,6),Fconvtab(127,10,59,6),Ttabeff(59),
     &     gtabeff(10),Ztabeff(6),filesize,taumin,taumax

c$$$      common /atmospheres/ rtau(127,60,10,7),rT(127,60,10,7),
c$$$     &     rhotab(127,60,10,7),Fconvtab(127,60,10,7),Ttabeff(60),
c$$$     &     gtabeff(10),Ztabeff(7),filesize

      common /metal/ FeH


      qtau(1:nmod) = 0.d0
      dqdtau(1:nmod) = 0.d0
      ddqtau(1:nmod) = 0.d0

      loggs = log10(g*m(neff)/(reff*reff))
      ntp = teff.gt.2750.d0.and.teff.lt.3500.d0.and.loggs.gt.-0.5d0.
     &     and.loggs.lt.1.5d0

      ntp2 = (teff.gt.7.2d3).and.(ntprof.eq.4)
!!      ntp2 = .false.

*____________________________________________________________
***   Grey Atmosphere - Eddington approximation
*------------------------------------------------------------

      if (ntprof.eq.0.or.(ntprof.eq.1.and..not.ntp).or.ntp2) then
c         print *,'entering ?'
         taueff = max(pw23,tau0)
         qtau0 = max(pw23,tau0)
         do i = nmod,1,-1
            if (tau(i).lt.taulim) then
               qtau(i) = max(pw23,tau0)
            endif
         enddo
         qtaus = max(pw23,tau0)
      endif

c$$$*____________________
c$$$***      Decoupling
c$$$*--------------------
c$$$
c$$$      if (ntprof.eq.2) then
c$$$         call atmos(Tatmd,rhoatmd,Ratmd,Patmd,Reff,Leff)
c$$$      endif
c$$$***************************

      if (ntprof.ne.1.and.ntprof.ne.3.and.ntprof.ne.4.and.ntprof.ne.5.
     &     or.ntp2)
     &     return

*____________________________________________________________
***   Atmosphere Models (Plez): analytic fits for cool giants
*------------------------------------------------------------

      if (ntprof.eq.1.and.ntp) then

         ts = teff
         gs = g*m(neff)/(reff*reff)
         gs1 = 1.d0/gs
         ms = m(neff)/msun
         ms1 = 1.d0/ms

         a1 = fa1(1)+fa1(2)*ts+fa1(3)*gs1+fa1(4)*ms
         a2 = fa2(1)+fa2(2)*ts+fa2(3)*gs1+fa2(4)*ms
         a3 = fa3(1)+fa3(2)*ms1+(fa3(3)+fa3(4)*ms1)*gs1+(fa3(5)+fa3(6)*
     &        ms1+(fa3(7)+fa3(8)*ms1)*gs1)*ts
         a4 = fa4(1)+fa4(2)*ts+fa4(3)*gs1+fa4(4)*ms
         a5 = fa5(1)+fa5(2)*ms1+(fa5(3)+fa5(4)*ms1)*ts+(fa5(5)+fa5(6)*
     &        ms1+(fa5(7)+fa5(8)*ms1)*ts)*gs1

         do i = nmod,1,-1
            if (tau(i).lt.taulim) then
               expq1 = a2*exv(a3*tau(i))
               expq2 = a4*exv(a5*tau(i))
               qtau(i) = a1+expq1+expq2
               dqdtau(i) = a3*expq1+a5*expq2
               ddqtau(i) = a3*a3*expq1+a5*a5*expq2
               if (abs(qtau(i)).lt.1.d-5) then
                  qtau(i) = 0.d0
                  dqdtau(i) = 0.d0
                  ddqtau(i) = 0.d0
               endif
            else
               goto 10
            endif
         enddo
 10      qtaus = qtau(nmod)
         qtau0 = 1.d0/dsqrt(3.d0)

      endif


*________________________________________________________________
***   Atmosphere Models (Kurucz + Eriksson + Plez): analytic fits
***   models used for PMS grid
*----------------------------------------------------------------

      if (ntprof.eq.3) then

         ts = teff
         gs = g*m(neff)/(reff*reff)

         alpha20 = 0.065d0
         beta0 = 5.7d2
         gamma0 = 0.0943d0
         delta0 = 2.7d0
         pwa1 = -1.1d0
         pwa2 = -1.d-2
         pwa3 = 1.d-2
         pwb1 = -4.1d0
         pwb2 = -4.d-2
         pwb3 = -8.d-2
         pwc1 = -0.7d0
         pwc2 = -5.d-2
c         pwc3 = -0.12d0
         pwc3 = 0.12d0
         pwd1 = -3.d0
         pwd2 = -0.17d0
         pwd3 = 0.24d0

         alphap = alpha20*(ts*1.d-4)**pwa1*(gs*1.d-4)**pwa2*
     &        (zeff/0.02d0)**pwa3
         betap = beta0*(ts*1.d-4)**pwb1*(gs*1.d-4)**pwb2*
     &        (zeff/0.02d0)**pwb3
         gammap = gamma0*(ts*1.d-4)**pwc1*(gs*1.d-4)**pwc2*
     &        (zeff/0.02d0)**pwc3
         deltap = delta0*(ts*1.d-4)**pwd1*(gs*1.d-4)**pwd2*
     &        (zeff/0.02d0)**pwd3

         qtau0 = -alphap-gammap

         do i = nmod,1,-1
            if (tau(i).lt.taulim) then
               qq1 = alphap*exp(-betap*tau(i))
               qq2 = gammap*exp(-deltap*tau(i))
               qtau(i) = -qq1-qq2
               dqdtau(i) = qq1*betap+qq2*deltap
               ddqtau(i) = -qq1*betap*betap-qq2*deltap*deltap
               if (abs(qtau(i)).lt.1.d-5) then
                  qtau(i) = 0.d0
                  dqdtau(i) = 0.d0
                  ddqtau(i) = 0.d0
               endif
!               dqdtau(i) = dqdtau(i)/1.d5
               qtau(i) = qtau(i)+pw23
            else
               goto 20
            endif
         enddo

c         write (850,2000) (qtau(i),tau(i),dqdtau(i),abrad(i), i=1,nmod)

! 20      qtaus = max(pw23,tau0)
 20      qtaus = qtau(nmod)
!         if (no.eq.maxmod) then
!            do i=1,nmod
!               write (2000,2010) i,tau(i),qtau(i),dqdtau(i)
!            enddo
!         endif

      endif

*********************************************************************************************
***   Atmosphere Models (PHOENIX): interpolation
***   models used for PMS grid 2016
*----------------------------------------------------------------

      if (ntprof==4.and..not.ntp2) then
         gs = g*m(neff)/(reff*reff)
         ts = teff
         qtau = 0.d0
         Tatm = 0.d0
         convatm = 0.d0
         rhoatm = 0.d0
!     print *,'## ATMOS : teff =',ts,' geff =',gs
!     &           ,' Zeff=',FeH
!     !            print *,tau(nmod-3),gs,taulim,nmod,ts,Tatm(nmod-3),
!     !     &           rhoatm(nmod-3),convatm(nmod-3)
         call select_table(tau,gs,taulim
     &        ,nmod,ts,Tatm,rhoatm,convatm)
!     !            print *,tau(nmod-3),gs,taulim,nmod,ts,Tatm(nmod-3),
!     !     &           rhoatm(nmod-3),convatm(nmod-3)
         
         k=nmod
         do while (tau(k+10)<=taulim)
            qtau(k) = pw43*(Tatm(k)/ts)**4.d0
     &           -tau(k)
            k=k-1
!!!   print *,tau(k),qtau(k)
         enddo
         valmax = 0.d0
         
         dqtautemp = 0.d0
         ddqtautemp = 0.d0
c     print *,(k,tau(k),Tatm(k),k=1,nsh)
         do i=1,k
            Tatms(i)=t(i)
            taumod(i)=tau(i)
         enddo
         do i=k+1,nmod
            Tatms(i)=Tatm(i)
            taumod(i)=tau(i)
         enddo
c     print *,(k,taumod(k),Tatms(k),k=1,nmod)
         taue = 0.d0
         
         call splineatmb (taumod,Tatms,nmod,1.d50,1.d50,ddqtautemp,
     &        dqtautemp)
         call splintatmb (taumod,Tatms,ddqtautemp,nmod,pw23,taue)
         
!!!         print *,'qtaue=',taue
         
         call splineatmb (Tatms,taumod,nmod,1.d50,1.d50,ddqtautemp,
     &        dqtautemp)
         
         call splintatmb (Tatms,taumod,ddqtautemp,nmod,ts,taueffb)
!!!         print *,'taueffb=',taueffb,ts
         call splineatm (tau,qtau,nmod,1.d50,1.d50,ddqtau,dqdtau)
         
         qtaus = qtau(nmod)
         qtau0 = 1.d0/dsqrt(3.d0)
         
         qtausv(1:nmod) = qtau(1:nmod)
         ddqtausv(1:nmod) = ddqtau(1:nmod)
         dqdtausv(1:nmod) = dqdtau(1:nmod)
         rhoatmsv(1:nmod) = rhoatm(1:nmod)
         
!     print *,ts,t(neff),tatm(neff)
         
         ts = teff
c     ts = (lum(nmod)/(pi*sig*r(nmod)*r(nmod)))**0.25d0
         gs = g*m(neff)/(r(nmod)*r(nmod))
!         gs = g*m(neff)/(reff*reff)

         rhoat(1:nmod) = rhoatm(1:nmod)


         alpha20 = 0.065d0
         beta0 = 5.7d2
         gamma0 = 0.0943d0
         delta0 = 2.7d0
         pwa1 = -1.1d0
         pwa2 = -1.d-2
         pwa3 = 1.d-2
         pwb1 = -4.1d0
         pwb2 = -4.d-2
         pwb3 = -8.d-2
         pwc1 = -0.7d0
         pwc2 = -5.d-2
c         pwc3 = -0.12d0
         pwc3 = 0.12d0
         pwd1 = -3.d0
         pwd2 = -0.17d0
         pwd3 = 0.24d0

         alphap = alpha20*(ts*1.d-4)**pwa1*(gs*1.d-4)**pwa2*
     &        (zeff/0.02d0)**pwa3
         betap = beta0*(ts*1.d-4)**pwb1*(gs*1.d-4)**pwb2*
     &        (zeff/0.02d0)**pwb3
         gammap = gamma0*(ts*1.d-4)**pwc1*(gs*1.d-4)**pwc2*
     &        (zeff/0.02d0)**pwc3
         deltap = delta0*(ts*1.d-4)**pwd1*(gs*1.d-4)**pwd2*
     &        (zeff/0.02d0)**pwd3

         qtau0 = -alphap-gammap

         do i = nmod,1,-1
            if (tau(i).lt.taulim) then
               qq1 = alphap*exp(-betap*tau(i))
               qq2 = gammap*exp(-deltap*tau(i))
               qtauLS(i) = -qq1-qq2
               dqdtauLS(i) = qq1*betap+qq2*deltap
               ddqtauLS(i) = -qq1*betap*betap-qq2*deltap*deltap
               if (abs(qtau(i)).lt.1.d-5) then
                  qtauLS(i) = 0.d0
                  dqdtauLS(i) = 0.d0
                  ddqtauLS(i) = 0.d0
               endif
               qtauLS(i) = qtauLS(i)+pw23
            endif
         enddo




c$$$         qtauKS(taulim:nmod) = 1.39d0 - 0.815d0 *exp(-2.54d0
c$$$     $        *tau(taulim:nmod)) -0.025d0 * exp(-30.d0
c$$$     $        *tau(taulim:nmod))
c$$$         dqdtauKS(taulim:nmod) = 2.0701d0*exp(-2.54d0*tau(taulim:nmod))
c$$$     &        + 0.75d0 * exp(-30.d0 *tau(taulim:nmod))
c$$$         ddqtauKS(taulim:nmod) = -5.258054d0*exp(-2.54d0*
c$$$     &        tau(taulim:nmod)) - 22.5d0 * exp(-30.d0 *tau(taulim:nmod))

         k=nmod
         do while (tau(k).lt.1.d3)
            k = k-1
         enddo

C         rewind (2000)
C         do i=k,nmod
C            write (2000,2000) i,tau(i),T(i),Tatm(i),qtau(i),dqdtau(i)
C!!!     &           ,qtauKS(i),dqdtauKS(i),qtauLS(i)
C         enddo
      endif

**************************************************
c.. Atmosphere from Krishna Swamy (1966)
*-------------------------------------------------
      if (ntprof.eq.5) then
         ts = teff
         a1 = 1.39d0
         a2 = 0.815d0
         a3 = -2.54d0
         a4 = 0.025d0
         a5 = -30.d0
         qtau(taulim:nmod) = a1-a2*exp(a3*tau(taulim:nmod))-a4*
     $        exp(a5*tau(taulim:nmod))
         dqdtau(taulim:nmod) = -a2*a3*exp(a3*tau(taulim:nmod))-a4*a5*
     &        exp(a5*tau(taulim:nmod))
         ddqtau(taulim:nmod) = -a2*a3*a3*exp(a3*tau(taulim:nmod))-
     &        a4*a5*a5*exp(a5*tau(taulim:nmod))
         qtaus = qtau(nmod)
      endif


 2000 format (i4,2x,5(1x,0pe12.5))
 2010 format (i4,2x,3(1x,0pe12.5))
 2020 format (i4,2x,8(1x,0pe12.5))

      return
      end
