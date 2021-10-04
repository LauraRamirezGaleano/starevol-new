      SUBROUTINE init

************************************************************************
* Initialize various quantities and determine the time-step            *
* 3003 constrain on nuclear luminosities
*
* ifail .le. 0 convergence successful or initial model
* ifail .gt. 0 failed model : retry
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
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.eos'
      include 'evolcom.flame'
      include 'evolcom.lum'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'


      integer icourant,itmax,ishockt,ishockb,imack
      integer neff1,nn,iburning
      integer idtnt,idtnf,idtnro,idtnr,idtnu,idtnl,idtnx
      integer nconvsh,ilum,ienvsh
      integer inshock,outshock
      integer klenv,klpulse,klcore
      integer i,j,k,kl,nq
      integer idiffvr

      character cmax*2,cburn*2

      double precision dtnflame
      double precision fatsh,xmaxH,tflash
      double precision mack(nsh),mackmax,mackmin,tshock
      double precision rne1l,rnel,lne1l,lnel,deltnel,deltes,dfact,syr,
     &     taune1l,taunel,tne1l,tnel,taueffl,facdtlum,faclc
      double precision dtnt,dtnf,dtnro,dtnr,dtnu,dtnl,dtnx,dtkh
      double precision dtiter,dtalt,dtsav,dtkhpms,dtnbrake
      double precision dtcour,dtcour0,dtnhydro,vdtn,dthdro,dvlum,
     &     dvlmax,dtlum,dtpulse,ftslum,slmax,Lmax,tmax,enucmax
      double precision menvconv,dmenvconv,mcore
      double precision TburnC,TburnNe
c      double precision stopage
      double precision diffvr,breaktime,disctime

      character*3 cyr
      character*10 stopage

      logical shockdetect
      logical endhb,endheb,ibrake

      common /burning/ TburnC,TburnNe
      common /hydrodyn/ Lmax,tmax,mackmax,enucmax,itmax,ishockb,ishockt
      common /meshconv/ menvconv,dmenvconv,mcore,nconvsh,ienvsh
      common /overshoot/ klcore,klenv,klpulse
      common /rotlaw/ diffvr,breaktime,idiffvr
      common /disclocking/ disctime
      
*_______________________________________________________
***   initializations for the computation of a new model
***   in case of failure or for initial model
*-------------------------------------------------------
      print *,'In init',xsp(1,2),xsp(10,2),xsp(20,2),xsp(30,2),xsp(40,2)
     $     ,xsp(50,2),xsp(60,2),xsp(70,2),xsp(80,2),xsp(90,2),xsp(100,2)

      ibrake = .false.
      ilum = 0
      if (no.eq.0.or.ifail.gt.0) then
***   surface quantities
         neff1 = neff-1
         taueffl = log(abs(taueff))
         taune1l = log(abs(tau(neff1)))
         taunel = log(abs(tau(neff)))
         tne1l = lnt(neff1)
         tnel = lnt(neff)
         deltnel = tne1l-tnel
         tphot = tne1l-(taune1l-taueffl)*deltnel/(taune1l-taunel)
         tphot = exp(tphot)
c         teff = tne1l-(taune1l-taueffl)*deltnel/(taune1l-taunel)
c         teff = exp(teff)
         taudnuc = 1.d99
         rne1l = lnr(neff1)
         rnel = lnr(neff)
         lne1l = log(abs(lum(neff1)))
         lnel = log(abs(lum(neff)))
         deltes = lnt(neff1)-log(tphot)
c         deltes = lnt(neff1)-log(teff)
         reff = rne1l-deltes*(rne1l-rnel)/deltnel
         leff = lne1l-deltes*(lne1l-lnel)/deltnel
         reff = exp(reff)
         leff = exp(leff)
         teff = (leff/(pim4*sig*reff*reff))**0.25d0

!!         teff = tphot

         solr = reff/rsun
         soll = leff/lsun
         teffy = teff
      endif
      dtalt = dtn

c..   in case of failure goto 30
      if (ifail.gt.0) goto 30

      call cpu_time (timebeg)
      vdtn = dtn


***   initialize variables
      Lmax = lum(nmod)
      enucmax = abs(venucl(1))
      tmax = t(1)
      itmax = 1
      do k = 2,nmod1
         if (t(k).gt.tmax) then
            tmax = t(k)
            itmax = k
         endif
         if (abs(lum(k)).gt.Lmax) Lmax = abs(lum(k))
         if (abs(venucl(k)).gt.enucmax) enucmax = abs(venucl(k))
      enddo

*______________________________________________________________
***   Define core and enveloppe mass coordinates and convective
***   boundaries where resolution needs to be increased
*--------------------------------------------------------------

      klcore = 0
      klpulse = 0
      klenv = 0
      nconvsh = 0
      mcore = 0.d0
      menvconv = 0.d0
      dmenvconv = 0.d0
      i = 1
      if (nsconv.gt.0) then
c..   define core
         if (novlim(1,3).eq.1) then
            klcore = 1
            if (mr(novlim(1,4)).lt.0.9d0) then
               mcore = m(novlim(1,4))
               i = 2
            endif
         endif
c..   define envelope
         do kl = i,nsconv
            if ((tau(novlim(kl,4)).lt.1.d2.or.t(novlim(kl,4)).lt.1.d6)
     &           .and.klenv.eq.0.and.mr(novlim(kl,4)).gt.0.7d0) then
               klenv = kl
               menvconv = m(novlim(klenv,3))
               dmenvconv = m(novlim(klenv,4))-menvconv
               ienvsh = novlim(klenv,3)
            else
c..   define pulse
               if (t(novlim(kl,3)).gt.2.d8.and.novlim(kl,3).gt.1.
     &              and.nphase.eq.5.and.klpulse.eq.0) klpulse = kl
c..   define regions where resolution needs to be increased
               if (novlim(kl,4)-novlim(kl,3).gt.20.or.(t(novlim(kl,3))
     &              .gt.4.d8.and.novlim(kl,4)-novlim(kl,3).gt.10)) then
                  nconvsh = nconvsh+1
                  novlim(nconvsh,7) = novlim(kl,3)
                  novlim(nconvsh,8) = novlim(kl,4)
               endif
            endif
         enddo
         if (nconvsh.gt.1) then
            do k = 1,int(nconvsh/2)
               j = novlim(k,7)
               novlim(k,7) = novlim(nconvsh+1-k,7)
               novlim(nconvsh+1-k,7) = j
               j = novlim(k,8)
               novlim(k,8) = novlim(nconvsh+1-k,8)
               novlim(nconvsh+1-k,8) = j
            enddo
         endif
      endif
      rgbphase = nphase.eq.3.and.teff.lt.teffagb.and.totm.le.2.5d0
      agbphase = nphase.eq.5.and.teff.lt.teffagb.and.totm.le.10.d0
      superagb = nphase.eq.5.and.teff.lt.teffagb.and.totm.ge.7.d0.and.
     &     totm.le.12.d0.and.xsp(1,ic12).lt.2d-2
c      postagb = nphase.eq.5.and..not.agbphase.and..not.agbphase.and.
c     &     dmenvconv/menvconv.lt.0.1d0

      if (klenv.gt.0) then
         hbb = t(novlim(klenv,3)).gt.4.d7.and.nuclopt.ne.'h'
      else
         hbb = .false.
      endif
      if (nsconv.ge.2) then
         flame = xsp(1,ine20).lt.0.1d0.and.t(novlim(1,3)).gt.5.d8.and.
     &        nphase.eq.6
         thermalpulse = lhe.gt.1.d3.and.t(novlim(1,3)).gt.1.5d8.and.
     &        agbphase.and.novlim(1,3).gt.1.and.m(novlim(1,4)).lt.
     &        1.4d0*msun
      else
         flame = .false.
         thermalpulse = .false.
      endif
      relaxpulse = .false.
      if (totnucl.gt.0.d0.and.klenv.gt.0) relaxpulse = agbphase.and.
     &     totgrav.lt.0.d0.and.m(novlim(klenv,3))/msun.lt.1.4d0.and.
     &     lhe.gt.1.d4

c..   envelope overshooting only during agbphase (if idup = 5)
      if (.not.agbphase.and.lover.eq.23.and.idup.eq.5) then
         diffover = .false.
         diffusconv = diffzc.or.semiconv.or.diffover
      endif

c..   determine deepest envelope penetration (DUPs)
      if (nphase.eq.3.or.(nphase.eq.5.and..not.thermalpulse)) then
         mdredgeup = min(mdredgeup,menvconv)
      else
         mdredgeup = menvconv
      endif

***   print headers
      if (no.eq.0) then
         if (time.lt.1.d0) then
            write (nout,1000) phase(nphase),model,min(soll,9999999.d0),
     &           min(9999.d0,solr),teffy,totm,dtn,time
         else
            write (nout,1001) phase(nphase),model,min(soll,9999999.d0),
     &           min(9999.d0,solr),teffy,totm,dtn*seci,time*seci
         endif
         modeli = model
      endif
      write (90,1100) model,nphase,phase(nphase)


*____________________________
***
***   TIMESTEP DETERMINATION
***
*----------------------------


***   duration of central nuclear burning phase
*----------------------------------------------

c..   H burning phase (Eggleton's book p36)
      dtmax = (1532.d0+2740.d0*mtini**4+146.d0*mtini**(5.5d0)+mtini**7)/
     &     (0.0397d0*mtini**2+0.3432d0*mtini**7)*1.d6*sec
      dtmax = min(dtmax*2.5d-3,5.d7*sec)
      xmaxH = 3.d-1
      TburnC = 5.8d8
      TburnNe = 9.d8
      if (nphase.eq.2.and.xsp(1,ih1).le.xmaxH) dtmax = dtmax*(0.1d0+
     &     0.9d0*xsp(1,ih1)/xmaxH)

      if (totm.le.15.d0.and.nphase.eq.4) dtmax = dtmax*0.3d0

c..   for massive stars (M > 15Mo)
      if (totm.gt.15.d0) then
c..   He burning phase : 0.200 dtmax --> 0.005 dtmax
         if (nphase.ge.4.and.xsp(1,ihe4).lt.5.d-1)
     &        dtmax = dtmax*(0.25d0+xsp(1,ihe4)*0.1d0)
c..   transition phase between He and C burning
c..   dt : 0.005 dtmax --> 50 yr & T < TburnC
         cburn = ' C'
         if (nphase.ge.5.and.tmax.le.TburnC) dtmax = dtmax*0.05d0*
     &        (TburnC-tmax)/tmax+50.d0*sec
c..   C burning  phase : 50 yr --> 0.1 yr & TburnC < T < 1.2d9
         if (nphase.ge.5.and.tmax.gt.TburnC.and.tmax.le.1.2d9) then
            dtmax = 50.d0*sec*(1.2d9-tmax)/tmax+0.1d0*sec
         endif
c..   Ne burning  phase
         if (tmax.gt.1.2d9.and.tmax.le.2.d9) then
            cburn = 'Ne'
         endif
c..   O burning  phase
         if (tmax.gt.2.d9) then
            dtmax = 1.d-3*sec
            cburn = ' O'
         endif
      endif

      if (dtmax.lt.dtmin) then
         dtmax = dtmin*facdt**(9-abs(icrash0))
         write (90,800) dtmax*seci
         write (nout,800) dtmax*seci
      endif

      if (hydro.and.ivisc.gt.0) dtmax = dtmax0

      ireset = 0
!      modatmos = .true.
      write (nout,1200) model+1
      dtalt = dtn
      dfact = min(dble(itermax)/dble(iter),facdt)
      dtiter = dtn*dfact


***   Shock fronts, Mack number and Courant time-step determinations
*-------------------------------------------------------------------

      shockdetect = .false.
      mackmax = 0.d0
      mackmin = prrc
      dtcour = 1.d-20
      tshock = 5.d8
      ishockb = nmod1
      ishockt = 1
      icourant = itmax
      imack = itmax
      if (model.eq.modeli.and.iter.eq.1.and.ifail.eq.0.and.hydro.and.
     &     ivisc.gt.0) write (90,2000) mackmin
      if (hydro.and.ivisc.gt.0.and.m(nmod1).gt.mcut.and.tmax.gt.tshock)
     &     then
***   determine largest Mack number in the high temperature region
         do i = 2,nmod1
            mack(i) =  abs(u(i))/dsqrt(pw53*p(i)/ro(i))
            if (mack(i).gt.mackmax.and.t(i).gt.1.d8) then
               mackmax = mack(i)
               imack = i
            endif
         enddo

c..   determine shock fronts (T > Tshock and Mack > Mackmin)
         inshock = nmod1
         outshock = nmod1
         do i = nmod1,2,-1
            if (m(i).lt.mcut) goto 20
c.. ori      outshock = i
            if (mack(i).gt.mackmin) outshock = i
            if (t(i+1).gt.tshock.or.t(i-1).gt.tshock) then
c.. ori         if (inshock.eq.nmod1) inshock = i
c.. ori         if (t(i+1).lt.tshock.and.t(i).ge.tshock.and.ishockt
c.. ori     &              .eq.nmod1) ishockt = i+1
c.. ori         if (t(i+1).gt.tshock.and.t(i).le.tshock) ishockb = i
               if (inshock.eq.nmod1.and.mack(i).gt.mackmin) inshock = i
               if (t(i+1).lt.tshock.and.t(i).ge.tshock.and.ishockt
     &              .eq.1.and.mack(i).gt.mackmin) ishockt = i+1
               if (t(i+1).gt.tshock.and.t(i).le.tshock.and.mack(i).gt.
     &              mackmin) ishockb = i
            endif
         enddo
 20      if (ishockt.ne.1) then
            ishockt = min(ishockt,inshock)
            if (ishockb.ne.nmod1) then
               ishockb = max(ishockb,outshock)
            else
               ishockb = i
            endif
         endif
         if (tmax.lt.tshock) then
            ishockt = min(ishockt,itmax+10)
            ishockb = max(ishockb,itmax-10)
         endif

         if (ishockt.ge.ishockb) shockdetect = .true.

***   determine minimum Courant time-steps in the shock region
         dtcour = 1.d99
         if (sigma.lt.0.4d0.and.mackmax.gt.0.3d0) then
            nq = 1
            do k = ishockb-nq,ishockt+nq
               dtcour0 = (r(k+1)-r(k))/dsqrt(pw53*p(k)/ro(k))
               if (dtcour0.lt.dtcour) then
                  dtcour = dtcour0
                  icourant = k
               endif
            enddo
         else
            dtcour = (r(ishockt+1)-r(ishockt))/dsqrt(pw53*p(ishockt)/
     &           ro(ishockt))
            dtcour0 = (r(ishockb+1)-r(ishockb))/dsqrt(pw53*p(ishockb)/
     &           ro(ishockb))
            icourant = ishockt
            if (dtcour0.lt.dtcour) then
               dtcour = dtcour0
               icourant = ishockb
            endif
         endif
      endif


***   evolution timescale of the structure
*------------------------------------------

      nn = nmod1
      if (ishtest.eq.'t'.or.ishtest.eq.'s') then
         do i = nmod1,2,-1
            if (nphase.gt.1) then
               if (t(i).gt.2.d6) goto 10
c               if (abs(enucl(i)).gt.shlim.and.t(i).gt.1.d6) goto 10
            else
               if (abs(enucl(i)).gt.1.d0.and.t(i).gt.1.d6) goto 10
            endif
         enddo
 10      nn = i+1
      endif


***   structural changes
*-----------------------

      if (accphase.eq.0) then
         iacctop = nmod
         iaccbot = 1
      endif
      call tstep (nn,dtnt,dtnf,dtnro,dtnr,dtnu,dtnl,idtnt,idtnf,
     &     idtnro,idtnr,idtnu,idtnl)
c      dtn = min(dtiter,dtnf,dtnt,dtnro,dtnacc)
      if (totm.lt.1.4d0.and.agbphase) then
         dtn = min(dtiter,dtnf,dtnacc)
      else
         dtn = min(dtiter,dtnf,dtnt,dtnacc)
      endif
c      if (rotation) dtn = min(dtn,dtnrot)
      if (ishtest.eq.'r') dtn = min(dtn,dtnr)


***   nuclear burning timescales
*-------------------------------

      dtnx = 1.d99
      idtnx = 1
      call tstepx (mtini,totm,nmod1,dtnx,idtnx,iburning)
      if (nphase.gt.2) then
         dtn = min(dtn,dtnx*10.d0)
      else
         dtn = min(dtn,dtnx)
      endif

***   Kelvin-Helmholtz timescale for contraction phases
*------------------------------------------------------

      dtkh = 1.d99
      endhb = xsp(1,ih1).lt.1.d-1.and.totm.ge.1.1d0.and.nphase.lt.3
      endheb = xsp(1,ihe4).lt.1.d-1.and.nphase.gt.3
      if (nphase.eq.1.or.nphase.eq.3.or.nphase.ge.5.or.endhb.or.
     &     endheb.or.(iaccr.gt.0.and.massrate.gt.0.d0)) then
         if (iaccr.gt.0.and.massrate.gt.0.d0) then
            if (accphase.eq.6.or.accphase.eq.7) then
               dtkh = g*menv*(m(nmod)-menv)/abs(r(nmod)*lum(nmod))
            else
               if (accphase.eq.8) then
                  dtkh = g*(1.d0-menv)*menv*m(nmod)**2/abs(r(nmod)*
     &                 lum(nmod))
               else
                  dtkh = g*m(iacctop)*(m(iacctop)-m(iaccbot))/
     &                 abs(r(nmod)*lum(nmod))
               endif
            endif
            dtn = min(dtkh*fkhdt,dtn)
            write (90,1450) vdtn*seci,dtnt*seci,idtnt,dtnf*seci,idtnf,
     &           dtnro*seci,idtnro,dtnacc*seci,cphase(iburning),
     &           dtnx*seci,idtnx,dtkh*seci,dtn*seci
         else if (nphase.eq.1) then
c            if (igw) then
c.. Timestep control on PMS for low-mass stars
c            if (novlim(nsconv,3).eq.1) then
               dtkhpms = m(nmod)*m(nmod)*g/abs(r(nmod)*lum(nmod))
c            else
               if (mtini.le.1.4d0) then             ! commenter en cas de modele trop lent sur pms
!                  dtkhpms = 1.d99
                  do k = 10,nmod
                     if (dtkhpms.gt.m(k)*m(k)*g/abs(r(k)*lum(k)).and.
     &                    m(k).gt.3.d-6*msun)
     &                    dtkhpms = m(k)*m(k)*g/abs(r(k)*lum(k))
                     dtkhpms = max(dtkhpms,1.d3*sec)
                  enddo
               else
                  dtkhpms = m(nmod)*m(nmod)*g/abs(r(nmod)*lum(nmod))
                  dtkhpms = min(dtkhpms,1.d4*sec)
               endif
!!            endif

!!            else
c$$$               if (egrav(1).lt.1.d0.and.novlim(nsconv,3).gt.1) then
c$$$c                  dtkh = dtkhpms
c$$$c                  dtn = min(dtkh*fkhdt,dtn)
c$$$                  dtkhpms = 1.d99
c$$$                  do k = 10,nmod
c$$$                     if (dtkhpms.gt.m(k)*m(k)*g/abs(r(k)*lum(k)).and.
c$$$     &                    m(k).gt.3.d-6*msun) then
c$$$                        dtkhpms = m(k)*m(k)*g/abs(r(k)*lum(k))
c$$$                     endif
c$$$                  enddo
c$$$                  dtkhpms = max(dtkhpms,1.d5*sec)
c$$$               else
c$$$                  dtkhpms = m(nmod)*m(nmod)*g/abs(r(nmod)*lum(nmod))
c$$$               endif
!!            endif
c            else
c               dtkhpms = m(nmod)*m(nmod)*g/abs(r(nmod)*lum(nmod))
c            endif
            dtkh = dtkhpms
            dtn = min(dtkh*fkhdt,dtn)

         else
            dtkh = m(nmod)*m(nmod)*g/abs(r(nmod)*lum(nmod))

c$$$c..   time step reduced during CNO buring on PMS
c$$$            if (nphase.eq.1.and.mtini.ge.0.6d0.and.novlim(1,3).eq.1
c$$$     $           .and.nmaxconv.gt.1) then
c$$$               dtkh = dtkh*0.1d0
            dtn = min(dtkh*fkhdt,dtn)
c$$$            endif
c$$$c..   time step reduced when the star approaches the ZAMS
c$$$            if (nphase.eq.1.and.novlim(1,3).gt.0) then
c$$$               dtkh = dtkh*5.d-2
c$$$               dtn = min(dtkh*fkhdt,dtn)
c$$$            endif
c---
c..   end of MS for rotating low-mass stars
            if (endhb.and.rotation.and.fkhdt.gt.0.5d0) dtn = min(dtkh*
     &           fkhdt*0.1d0,dtn)
            write (90,1400) vdtn*seci,dtnt*seci,idtnt,dtnf*seci,idtnf,
     &           dtnro*seci,idtnro,cphase(iburning),dtnx*seci,idtnx,
     &           dtkh*seci,dtn*seci
         endif
      else
         if (iaccr.gt.0) then
            write (90,1350) vdtn*seci,dtnt*seci,idtnt,dtnf*seci,idtnf,
     &           dtnro*seci,idtnro,dtnacc*seci,dtnx*seci,idtnx,
     &           cphase(iburning),dtn*seci
         else
            write (90,1300) vdtn*seci,dtnt*seci,idtnt,dtnf*seci,idtnf,
     &           dtnro*seci,idtnro,cphase(iburning),dtnx*seci,idtnx,
     &           dtn*seci
         endif
      endif


***   constraints based on nuclear evolutionary phase
*------------------------------------------------------

      dtsav = dtn
      if (nphase.le.2.and.totm.gt.2.d0) dtn = min(dtmax,dtn)
      if (nphase.eq.4) dtn = min(dtmax,dtn)
      if (nphase.ge.6) dtn = min(dtmax,dtn)
      if (dtn.lt.dtsav.and.nphase.eq.2) endhb = .false.
      if (dtn.lt.dtsav.and.nphase.eq.4) endheb = .false.
      if (dtn.lt.dtsav.and.nphase.gt.4) ilum = 3


***   hydrodynamical constraints
*-------------------------------

c..   Magnetic braking at the ZAMS (rotating case)
      if (mtini.ge.0.6d0.and.rotation.and.((idiffvr.eq.5.or.
     $     idiffvr.ge.8).and.xsp(1,ih1)/xsp(nmod,ih1).le.0.99d0).and.
     $     xsp(1,ih1)/xsp(nmod,ih1).gt.0.9979d0) then
         ibrake = .true.
      endif


      if (iaccr.gt.0.and.accphase.gt.4.and.tmax.gt.2.d8) iaccr = 0
      if (shockdetect.and.tmax.gt.2.d9.and.time.gt.1.d10.and.dtn.lt.
     &     dtcour*1.d2) time = 0.d0
c..   use Courant time-step
      tflash = 1.d8
      if (shockdetect.and.tmax.gt.tflash) then
         dthdro = 1.d-8
         dtnhydro = max(0.d0,(2.5d9-tmax)/(2.5d9-tflash))
         dtnhydro = dtnhydro**4
         if (vdtn.gt.dtcour) facdt = facdt0*2.d0
         dtnhydro = dtnhydro*dtmax
         fatsh = 1.d0
         if (tmax.lt.2.d9) then
            fatsh = 1.d2
            dtnhydro = dtnhydro+dthdro
         endif
c         dtcour = max(dtcour,5.d-8/ftsh)
         if (dtcour.gt.dthdro) then
            dtn = min(vdtn*facdt,(dtnhydro+dtcour*ftsh*fatsh))
         else
            dtn = min(vdtn*facdt,dthdro)
         endif
      endif


***   constraints based on nuclear luminosities
***   (works for no > 1)
*----------------------------------------------

      if (ftnuc.lt.1.d10) then
         ftslum = ftst
         dvlmax = 0.d0
         dvlum = 0.d0
         slmax = Lmax/lsun
c         if (vlpp.gt.0.d0.and.lpp.gt.0.d0.and.lpp/slmax.gt.1.d-1) then
c            dvlmax = abs(log10(lpp/vlpp))
c            cmax = 'pp'
c         endif
         if (vlh.gt.0.d0.and.lh.gt.0.d0.and.lh/slmax.gt.1.d-1) then
            dvlum = abs(log10(lh/vlh))
            if (dvlum.gt.dvlmax) then
               dvlmax = dvlum
               cmax = 'H '
            endif
         endif
         if (vlhe.gt.0.d0.and.lhe.gt.0.d0.and.lhe/slmax.gt.1.d-1) then
            dvlum = abs(log10(lhe/vlhe))
            if (dvlum.gt.dvlmax) then
               dvlmax = dvlum
               cmax = 'He'
            endif
         endif
c..   if flame detected increase accuracy
         faclc = 1.d-1
c      if (flame) then
c         faclc = 5.d-2
c         ftslum = faclc
c      endif
         if (vlc.gt.0.d0.and.lc.gt.0.d0.and.lc/slmax.gt.faclc) then
            dvlum = abs(log10(lc/vlc))
            if (dvlum.gt.dvlmax) then
               dvlmax = dvlum
               cmax = 'C '
            endif
         endif
         if (abs(vlne).gt.0.d0.and.abs(lne).gt.0.d0.and.abs(lne)/slmax
     &        .gt.1.d-1) then
            dvlum = abs(log10(lne/vlne))
            if (dvlum.gt.dvlmax) then
               dvlmax = dvlum
               cmax = 'Ne'
            endif
         endif
         if (dvlmax.gt.log10(1.d0+ftslum)) then
            facdtlum = facdt0
            if (dvlmax.gt.0.3d0) facdtlum = 10.d0**dvlmax
c..init236         dtlum = dtalt/min(10.d0**dvlmax,facdt0)
c..init231         dtlum = dtalt/max(10.d0**dvlmax,facdt0)
            dtlum = dtalt/facdtlum
            if (dtlum.lt.dtn) then
               ilum = 1
               dtn = dtlum
            endif
         endif
         dvlmax = min(dvlmax,99999.d0)
      endif

***   constraints based on the presence of a thermal pulse
*---------------------------------------------------------

      if ((thermalpulse.or.relaxpulse).and.dtn.gt.sec) then
c..   thibaut
c         dtpulse = min(facdt*dtalt,dtn)
         dtpulse = max(1.d0*sec,dtalt/facdt0,dtmin*sec)
         if (dtpulse.lt.dtn) then
            ilum = 2
            dtn = dtpulse
         endif
      endif


*____________________________________________________________
***   Final time-step constraints from the previous time-step
*------------------------------------------------------------

      dtn = min(facdt*dtalt,dtn)

c      if (hydro.and.(q0.gt.0.d0.or.ivisc.gt.0).and.tmax.gt.2.d8)
c     &        dtn = min(dtn,dtcour*fts)


***   Treatment of Flames
*------------------------

      if (imodpr.eq.11.or.imodpr.eq.12) dtn = min(dtn,dtnflame(dtn))

      if (nphase.gt.1.and.nphase.lt.6.and..not.dup3) then
         dtn = max(dtalt/facdt**3,dtn)
      else
         dtn = max(dtalt/facdt**2,dtn)
      endif

c..   first model : dtn = dtin
      if (model.eq.modeli.and.ireset.eq.0.and.dtin.ne.0.d0) dtn = dtin

c..   luminosity variation too large. Do not increase dtn
      if (ifail.lt.0) then
         dtn = min(dtn,dtalt)
         write (nout,*) 'luminosity variation too large, dtn not ',
     &        'increased'
      endif

*____________________________________________________________
***   Treatment of core hydrogen burning exhaustion (modif Louis 09/2018)
*------------------------------------------------------------

      xmaxH = 1.d-1
      if (nphase.eq.2.and.xsp(1,ih1).le.xmaxH) dtn = min(dtn,dtmax0*
     &     (0.05d0+0.95d0*xsp(1,ih1)/xmaxH))
c$$$      if (nphase.eq.2.and.xsp(1,ih1).le.xmaxH) dtn = min(dtn,dtmax0*
c$$$     &     (0.005d0+0.95d0*xsp(1,ih1)/xmaxH))                                  ! In case of atomic diffusion
      
*____________________________________________________________
***   Treatment of centrifugal forces in rotating models
*------------------------------------------------------------

 30   dtn = max(dtn,dtmin)
      if (hydrorot.and.nphase.eq.1.and.novlim(nsconv,3).eq.1) then
         dtn = min(dtn,dtmax0*5.d-3)
      else
         dtn = min(dtn,dtmax0)
      endif
      if (ibrake) dtn=1.d4/seci
      


*_________________
***   print banner
*-----------------

      if (ifail.le.0) then
         dtn0 = dtn
         if (dtn.lt.1.d2) then
            cyr = ' s '
            syr = 1.d0
         else
            cyr = ' yr'
            syr = seci
         endif
         if (nphase.eq.1.or.nphase.eq.3.or.nphase.ge.5.or.endhb.or.
     &        endheb) then
            if (iaccr.eq.0) then
               write (nout,1700) dtalt*syr,cyr,dtnt*syr,cyr,idtnt,dtnf*
     &              syr,cyr,idtnf,cphase(iburning),dtnx*syr,cyr,idtnx,
     &              dtnro*syr,cyr,idtnro,turnenv*syr,cyr,dtnu*syr,cyr,
     &              idtnu,dtnr*syr,cyr,idtnr,dtnl*seci,idtnl,dtkh*seci,
     &              tmax,mr(itmax)*totm,itmax
            else
               write (nout,1750) dtalt*syr,cyr,dtnt*syr,cyr,idtnt,dtnf*
     &              syr,cyr,idtnf,cphase(iburning),dtnx*syr,cyr,idtnx,
     &              dtnro*syr,cyr,idtnro,turnenv*syr,cyr,dtnacc*syr,cyr,
     &              dtnu* syr,cyr,idtnu,dtnr*syr,cyr,idtnr,dtnl*syr,cyr,
     &              idtnl,dtkh*seci,tmax,mr(itmax)*totm,itmax
            endif
         else
            if (iaccr.eq.0) then
               write (nout,1800) dtalt*syr,cyr,dtnt*syr,cyr,idtnt,dtnf*
     &              syr,cyr,idtnf,cphase(iburning),dtnx*syr,cyr,idtnx,
     &              dtnro*syr,cyr,idtnro,turnenv*syr,cyr,dtnu*syr,cyr,
     &              idtnu,dtnr*syr,cyr,idtnr,dtnl*syr,cyr,idtnl,dtmax*
     &              seci,tmax,mr(itmax)*totm,itmax
            else
               write (nout,1850) dtalt*syr,cyr,dtnt*syr,cyr,idtnt,dtnf*
     &              syr,cyr,idtnf,cphase(iburning),dtnx*syr,cyr,idtnx,
     &              dtnro*syr,cyr,idtnro,turnenv*syr,cyr,dtnacc*syr,cyr,
     &              dtnu*syr,cyr,idtnu,dtnr*syr,cyr,idtnr,dtnl*syr,cyr,
     &              idtnl,dtmax*seci,tmax,mr(itmax)*totm,itmax
            endif
         endif
         if (hydro.and.ivisc.gt.0.and.shockdetect) then
            write (nout,1580) dtcour,icourant,t(icourant),ishockb,
     &           ishockt,mr(ishockb)*totm,mr(ishockt)*totm,mack(imack)
            write (90,1580) dtcour,icourant,t(icourant),ishockb,ishockt,
     &           mr(ishockb)*totm,mr(ishockt)*totm,mack(imack)
         endif
c         if (irotbin.eq.1) write (nout,1860) dtnrot*syr,cyr,idtnrot
         if (ilum.eq.1) then
            write (nout,900) cmax,(10.d0**dvlmax-1.d0)*1.d2,ftslum*1.d2,
     &           dtlum*seci
         elseif (ilum.eq.2) then
            write (nout,1950) dtn*syr,cyr
         elseif (ilum.eq.3) then
            write (nout,1970) cburn,dtn*syr,cyr
         else
            write (nout,1900) dtn*syr,cyr
         endif
      endif

      phi = min(phi0*dtn/dtalt,1.d0)
      time = time+dtn

*____________________________________________________________________
***   Modification of the time-step in case it will lead to overpass
***   the adopted age of the Sun (Ana 04/11/2005)
***   To be used only if calibrating a solar model
*--------------------------------------------------------------------

      if (mtini.eq.1.0d0.and.nphase.eq.2) then
         if (time.gt.4.57d9*sec.and.(time-dtn).lt.(4.57d9*sec)) then
            time = time-dtn
            dtn = 4.5700000d9*sec-time
            time = time+dtn
            maxmod = no
            write(*,1980) dtn*syr
         endif
      endif

*____________________________________________________________________
***   Modification of the time-step in case it will lead to overpass
***   the disk-locking time (Louis 16/01/2013)
***   To be used only if idiffvr = 5, 8 or 9
*--------------------------------------------------------------------
      if (nphase.eq.1.and.(idiffvr.eq.5.or.idiffvr.ge.8.or.
     &     idiffvr.eq.4)) then
         if (time.gt.disctime*sec.and.(time-dtn).lt.(disctime*sec)) then
            time = time-dtn
            dtn = disctime*sec-time
            time = time+dtn
            maxmod = no
            write(*,1980) dtn*syr
         endif
      endif

c      call getenv('STOPAGE',stopage)
c      print *,'stopage',stopage
c      if (time.gt.stopage*sec) then
c         time = time-dtn
c         dtn = stopage*sec-time
c         time = time+dtn
c         maxmod = no
c         write(*,1980) dtn*syr
c      endif

      do kl = 1,nmaxconv
         turnover(kl) = 0.d0
      enddo

 800  format ('WARNING dtmin < dtmax : dtmax increased to ',1pe9.3,'yr')
 900  format (/,' too large increase in ',a2,' luminosity : ',0pf6.0,
     &     '% > ',0pf3.0,', time step reduced to ',1pe9.3,/)
 1000 format (/,132('_'),//,16x,a6,' phase:',/,1x,'initial model:',4x,
     &     'model',5x,'L',8x,'R',6x,'Teff',7x,'M',9x,'dt (s)',8x,
     &     't (s)',/,16x,75('-'),//,16x,i8,1x,f10.2,1x,f7.2,1x,f7.0,
     &     1x,f11.8,1x,1pe10.4,1x,1pe16.10,/)
 1001 format (/,132('_'),//,16x,a6,' phase:',/,1x,'initial model:',4x,
     &     'model',5x,'L',8x,'R',6x,'Teff',7x,'M',8x,'dt (yr)',7x,
     &     't (yr)',/,16x,75('-'),//,16x,i8,1x,f10.2,1x,f7.2,1x,f7.0,
     &     1x,f11.8,1x,1pe10.4,1x,1pe16.10,/)
 1100 format ('model #',i8,', Phase ',i1,': ',a6)
 1200 format (132('_'),//,1x,'model #',i8,':')
 1300 format (/,1x,'dtni =',1pe9.3,': T =',1pe9.3,' #',i4,'; lnf =',
     &     1pe9.3,' #',i4,'; rho =',1pe9.3,' #',i4,'; ',a4,'= ',1pe9.3,
     &     ' #',i4,'; dtn =',1pe9.3)
 1350 format (/,1x,'dtni =',1pe9.3,': T =',1pe9.3,' #',i4,'; lnf =',
     &     1pe9.3,' #',i4,'; rho =',1pe9.3,' #',i4,'; acc =',1pe9.3,
     &     '; ',a4,'=',1pe9.3,' #',i4,'; dtn =',1pe9.3)
 1400 format (/,1x,'dtni =',1pe9.3,': T =',1pe9.3,' #',i4,'; lnf =',
     &     1pe9.3,' #',i4,'; rho =',1pe9.3,' #',i4,'; ',a4,'=',1pe9.3,
     &     ' #',i4,'; KH =',1pe9.3,'; dtn =',1pe9.3)
 1450 format (/,1x,'dtni =',1pe9.3,': T =',1pe9.3,' #',i4,'; lnf =',
     &     1pe9.3,' #',i4,'; rho =',1pe9.3,' #',i4,'; acc =',1pe9.3,
     &     '; ',a4,'=',1pe9.3,' #',i4,'; KHacc =',1pe9.3,'; dtn =',
     &     1pe9.3)
 1580 format (1x,'shock : Courant time = ',1pe9.3,' s (#',i4,', T = ',
     &     1pe9.3,'), location (',i4,',',i4,'), Mr (',0pf7.4,',',0pf7.4,
     &     '), Mack_Tm = ',0pf7.4)
 1700 format (/,1x,'dtni = ',1pe12.6,a3,'; from T: ',1pe10.4,a3,' (#',
     &     i4,'), from lnf: ',1pe10.4,a3,' (#',i4,'), from ',a4,': ',
     &     1pe10.4,a3,' (#',i4,')',/,25x,'from rho: ',1pe10.4,a3,' (#',
     &     i4,'), from envelope turnover: ',1pe10.4,a3,/,25x,
     &     '[ from u: ',1pe10.4,a3,' (#',i4,'), r: ',1pe10.4,a3,
     &     ' (#',i4,'), L: ',1pe10.4,' (#',i4,') ]',/,1x,
     &     'relevant Kelvin-Helmholtz time-scale = ',1pe10.4,' yr',
     &     ', Tmax = ',1pe9.3,' (Mr = ',0pf6.4,', #',i4,')')
 1750 format (/,1x,'dtni = ',1pe12.6,a3,'; from T: ',1pe10.4,a3,' (#',
     &     i4,'), from lnf: ',1pe10.4,a3,' (#',i4,'), from ',a4,': ',
     &     1pe10.4,a3,' (#',i4,')',/,25x,'from rho: ',1pe10.4,a3,' (#',
     &     i4,'), from envelope turnover: ',1pe10.4,a3,/,16x,
     &     '[ from Macc: ',1pe10.4,a3,', u: ',1pe10.4,a3,' (#',i4,
     &     '), r: ',1pe10.4,a3,' (#',i4,'), L: ',1pe10.4,a3,' (#',i4,
     &     ') ]',/,1x,'relevant Kelvin-Helmholtz time-scale in the ',

     &     'accretion region = ',1pe10.4,' yr',
     &     ', Tmax = ',1pe9.3,' (Mr = ',0pf6.4,', #',i4,')')
 1800 format (/,1x,'dtni = ',1pe12.6,a3,'; from T: ',1pe10.4,a3,' (#',
     &     i4,'), from lnf: ',1pe10.4,a3,' (#',i4,'), from ',a4,': ',
     &     1pe10.4,a3,' (#',i4,')',/,25x,'from rho: ',1pe10.4,a3,' (#',
     &     i4,'), from envelope turnover: ',1pe10.4,a3,/,25x,
     &     '[ from u: ',1pe10.4,a3,' (#',i4,'), r: ',1pe10.4,a3,
     &     ' (#',i4,'), L: ',1pe10.4,a3,' (#',i4,') ]',/,1x,
     &     'characteristic core burning timescale = ',1pe10.4,' yr',
     &     ', Tmax = ',1pe9.3,' (Mr = ',0pf6.4,', #',i4,')')
 1850 format (/,1x,'dtni = ',1pe12.6,a3,'; from T: ',1pe12.6,a3,' (#',
     &     i4,'), from lnf: ',1pe12.6,a3,' (#',i4,'), from ',a4,': ',
     &     1pe12.6,a3,' (#',i4,')',/,25x,'from rho: ',1pe12.6,a3,' (#',
     &     i4,') from envelope turnover: ',1pe10.4,a3,/,16x,
     &     '[ from Macc: ',1pe12.6,a3,', u: ',1pe12.6,a3,' (#',i4,
     &     '), r: ',1pe12.6,a3,' (#',i4,'), L: ',1pe12.6,a3,' (#',i4,
     &     ') ]',/,1x,'characteristic core burning timescale = ',
     &     1pe10.4,' yr, Tmax = ',1pe9.3,' (Mr = ',0pf6.4,', #',i4,')')
c 1860 format (/,16x,'[ from xpsi (rotation): ',1pe10.4,a3,' (#',i4,
c     &     ') ]')
 1900 format (/,1x,'first try: dtn = ',1pe12.6,a3,/)
 1950 format (/,1x,'time step constrained by pulse/3DUP, first try: ',
     &     'dtn = ',1pe12.6,a3,/)
 1970 format (/,1x,'time step constrained from ',a2,' burning, first ',
     &     'try: dtn = ',1pe12.6,a3,/)
 1980 format (/,1x,'Computation of solar model : timestep modified for
     &     the last model computation',1pe12.6,/)
 2000 format (2x,'** Shock declare if Mack number > Mackmin = ',0pf5.3)

      return
      end
