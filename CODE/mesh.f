      SUBROUTINE mesh

************************************************************************
*     setup the mesh grid-point                                        *
*     WARNING : never use t, use lnt !                                 *
*     In case of rotational mixing, add a criterion on shell mass      *
*     (in mass) for grid points removal.                               *
*                                                                      *
* $LastChangedDate:: 2016-05-11 17:20:46 +0200 (Mer, 11 mai 2016)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 62                                                          $ *
*                                                                      *
************************************************************************
c..   14/08 : mesh add test on enucmax
c..   modifs : freeze mesh if lnt < 4.d8 and nphase >= 6
c..   adaptative mesh resolution
c..   17/09 : in case of rotation, mesh frozen for the first model
c..   following the CC disappearance
c..   20/02 : increase resolution at convective boundaries+ bugs
c..   08/05 : add test on flame phase for denn
c..   05/10 : rewrite generalized logic

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.flame'
      include 'evolcom.ion'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.rot'
      include 'evolcom.conv'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.transp'
      include 'evolcom.var'


      integer icrzc(nsh),icz,icz1,icz0,kdup,icount
      integer idwn,iup,kp,km,is,imesh,ienvsh,ndiff_ov
      integer ibov,nresol,nrconv,nconvmin,ndup
      integer nshdel(nsh),nshadd(nsh),nshad(nsh),kshell(nsh)
      integer natm,nshock
      integer ie,io,mm0,mm,mm1,mm2
      integer itmax,ishockb,ishockt,kshockb,kshockt
      integer iadd,iremove,ir,il,ic,nmod0,nconvsh
      integer i,j,k,ij,ik
      integer iHBSb,iHBSt,iCEb,iCEt,iionH,iionHe
      integer iIGWb,iIGWt
      integer fcs         ! first convective shell

      logical meshfrozen,zonetest,meshcrit,pass,azer
      logical duptest,shellcenter
      logical structest,tautest,hydrotest,preshocktest,pvistest,
     &     courantest,shockdetect,center,edge,chimi,edge_bis
      logical flametest,incrHBS,incrIGW,incrATM
      logical addition,suppression,debug
      logical solmod, diffumod
      logical maskGamma(nmod)

      double precision tshock,Lmax,tmax,Mackmax
      double precision lntfreeze,lntmesh,meshcut
      double precision lntc,lnfc,lnroc,lnpc,venuclc
      double precision facshmin,lntmax,massfactor
      double precision dfnf,dror,dtnt,dlnl,dpnp,dmnm,denn,domega,ddmnm
      double precision dnucmax,dvmax,dvmin,dlummax,drvmax,fdm
      double precision dtaur,dtaumesh,drr,dxmue,dxnx,dtaumeshbis
      double precision tausup,tauinf,mrmax,mshockb,mshockt,
     &     dmrmax,dmrmin,dmrcore,dmrenv,lumlim,enucmax,enuclim
c      double precision dpvis,dmshock
      double precision dtcourf,dtcourb,dtcour,facHBS
      double precision dlmnm,xmresol,xmrconv
      double precision tdust

      double precision facIGW,dtauaccu,facION

      double precision rhp
      double precision facCE

      double precision menvconv,dmenvconv,mcore,dmcore,ddmcore,dmenv,
     &     ddmenv,dmcz,dmcz0,ddmcz,dmrcz

c     Add by T.Dumont Jan. 2018 -> variations allowed between two successives shells 
      double precision drmax,dHmax,dHe4max,dO16max,dFemax,delemmax,
     &     dPtotmax,dTempmax

      double precision times,chronos ! Ajout chronos pour mesure du tps de calcul
      double precision disctime     ! modif TD Fev.2019

      common /meshconv/ menvconv,dmenvconv,mcore,nconvsh,ienvsh
      common /hydrodyn/ Lmax,tmax,Mackmax,enucmax,itmax,ishockb,ishockt
      common /disclocking/ disctime ! modif TD Fev.2019

      debug = .false.
      icount = 0
      pass = .false.
      massfactor = 1.63d0
      nconvmin = 20
      dlummax = 1.d-2
      dmrmax = dmrma
      dmrmin = dmrmi
c      lntmax = log(tmax*ftacc)
      lntmax = log(tmax)

      if (imodpr.eq.11.or.imodpr.eq.12) then
         do i = 1,nbcz
            dlummax = min(dlummax,(lum(itCOcz(i))/Lmax)*0.5d0)
         enddo
      endif
      lumlim = Lmax*dlummax
      enuclim = enucmax*1.d-3
      if (agbphase.and.enuclim.gt.1.d2) enuclim = 1.d2
      print *,'debut', 'dmrmax=',dmrmax, 'dmrmin=',dmrmin
c      print *, 'tmax=', tmax, 'lntmax', lntmax
      print *, 'lumlim', lumlim, 'Lmax', Lmax, 'dlummax', dlummax
c     print *, 'nshock',nshock,'dtcourf',dtcourf,'dmnm',dmnm
c     print *, 'enuclim', enuclim
c      print *,'crz', crz
      DO i = 1,nsh
         IF (crz(i+1).LT.0.AND.crz(i).GT.1) THEN
c            PRINT *, 'last radiative shell ',i
c           last_rad = i 
c            PRINT *, 'first convective shell ', i+1
            fcs = i+1
c            first_conv = i+1
         ENDIF
      ENDDO
*-------------------
***   3DUP detection
*-------------------

      ndiff_ov = 0
      if (agbphase) then
         ibov = ienvsh
         if (diffover.and.crz(ienvsh-1).eq.1) then
            ibov = ibov-1
            do while (crz(ibov).eq.1.and.crz(ibov-1).eq.1)
               ibov = ibov-1
            enddo
         endif
         ndiff_ov = ienvsh-ibov
c..     ndup : number of shells below the envelope where the composition
c..     is checked
         ndup = 20
         if (dup3) ndup = 50
         dup3 = xsp(ienvsh,ih1)/xsp(ibov-ndup,ih1).gt.2.d1.and.
     &        xsp(ibov-ndup,ic12).gt.1.d-1
      else
         dup3 = .false.
      endif
*----------------------------
***   HYDRODYNAMICS of SHOCKS
*----------------------------

      io = 0
      ie = 0
      kshockb = ishockb
      kshockt = ishockt
      forall (i=1:nsh)
         kshell(i) = 0
         nshdel(i) = 0
         nshadd(i) = 0
         nshad(i) = 0
      end forall
      shockdetect = .false.
      tshock = 5.d8
      nshock = nretry           ! number of shells in the shock neighborhood
      hydrotest = hydro.and.q0.gt.0.d0.and.ivisc.gt.0.and.tmax.gt.
     &     tshock
      if (hydrotest.and.maxsh.ne.0.and.ishockb.le.ishockt) then
         shockdetect = .true.
         dvmax = dlnvma
         drvmax = dlnvma
         dnucmax = dlnenuc
c..   inward shockfront
         do k = nshock,1,-1
            ij = ishockb-k+1
            dtcourb = (r(ij+1)-r(ij))/dsqrt(pw53*p(ij)/ro(ij))
            if ((io+ie).le.maxsh) then
               if (dtcourb.gt.1.d-4) then
                  io = io+1
                  kshell(ij) = kshell(ij)+1
                  nshadd(io) = ij
                  write (nout,*) ij,'add1',dtcourb,', iadd:',io,
     &                 ', idel:',ie,exp(lnt(ij))
               elseif (dtcourb.lt.1.d-7) then
                  ie = ie+1
                  kshell(ij) = kshell(ij)-1
                  nshdel(ie) = ij
                  write (nout,*) ij,'del1',dtcourb,', iadd:',io,
     &                 ', idel:',ie,exp(lnt(ij))
               endif
            endif
         enddo
c..   outward shockfront
         do k = 1,nshock
            ik = ishockt+k-1
            dtcourf = (r(ik+1)-r(ik))/dsqrt(pw53*p(ik)/ro(ik))
            if ((io+ie).le.maxsh) then
               if (dtcourf.gt.1.d-3) then
                  io = io+1
                  kshell(ik) =  kshell(ik)+1
                  nshadd(io) = ik
                  write (nout,*) ik,'add2',dtcourb,', iadd:',io,
     &                 ', idel:',ie,exp(lnt(ik))
               elseif (dtcourf.lt.1.d-7) then
                  ie = ie+1
                  kshell(ik) =  kshell(ik)-1
                  nshdel(ie) = ik
                  write (nout,*) ik,'del2',dtcourb,', iadd:',io,
     &                 ', idel:',ie,exp(lnt(ik))
               endif
            endif
         enddo
c..   delete shells inside shock fronts
         do k = ishockb+nq0,ishockt-nq0
            dtcourf = (r(k+1)-r(k))/dsqrt(pw53*p(k)/ro(k))
            dmnm = (dm(k-1)+dm(k))/m(nmod)
c            print *,'dmnm', dmnm, 'dmrmax', dmrmax, 'dm',dm
            if (dtcourf.lt.1.d-5.and.(io+ie).le.maxsh.and.
     &           dmnm.lt.dmrmax) then
               ie = ie+1
               kshell(k) =  kshell(k)-1
               nshdel(ie) = k
               write (nout,*) k,dtcourf,', iadd:',io,', idel:',ie,dmnm
     &              ,m(k+1)/msun-m(k)/msun
            endif
         enddo
*     add shells in the shock region (peak temperature)
         mrmax = mr(itmax)
         mshockb = mrmax*0.97d0
         mshockt = mrmax*1.07d0
c         dmshock = (mshockt-mshockb)*1.d-2
         iremove = ie
         iadd = io+1
c..   do not add shells
c         maxsh = ie
c         io = 0
         if (ishtest.eq.'t'.or.ishtest.eq.'s') then
            if (ie.ne.0) then
               print *, 'goto 500'
               goto 500
            else
               nmod0 = nmod
               print *,'goto 700'
               goto 700
            endif
         endif
      else
         mshockb = 0.d0
         mshockt = mr(nmod)
      endif


*--------------------
*** ROTATIONAL MIXING
*--------------------

c..   in case of rotational, increase resolution in HBS
      incrHBS = nphase.gt.2.and.nphase.lt.5.and.rotation
      if (incrHBS) then
         facHBS = 2.d0
         iHBSb = 0
         iHBSt = 0
         do k = 1,nmod
            if (xsp(k,ih1).gt.1.d-5.and.xsp(k,ihe4).lt.xsp(k,ih1)
     &           .and.enucl(k).gt.1.d2) then
               iHBSb = k
               exit
            endif
         enddo
         do k = iHBSb,nmod
            if (xsp(k,ih1).gt.4.d-1) then
               iHBSt = k
               exit
            endif
         enddo
         if (iHBSt.eq.nmod) incrHBS = .false.
         write (nout,'(" HBS location : [",i4,",",i4,"]")') iHBSb,iHBSt
      endif
c..  modify mesh in case of rotational mixing at the end of H burning
c      if (nphase.ge.2.and.xsp(1,ih1).lt.1.d-2.and.totm.lt.2.d0.and.
c     &     irotbin.eq.1) maxsh = 0


c.. In case of rotation + igw : increase resolution in the core.

      incrIGW = nphase.eq.2.and.igwrot
      if (incrIGW) then
         facIGW = 4.d0
         k=1
         do while (enucl(k).gt.0.6d0*enucl(1))
            k = k+1
         enddo
         iIGWt = k
         iIGWb = 1
      endif
      
c.. In case of atmospheric model, increase resolution in ionization region

      incrATM = ntprof.ne.0
!      incrATM = .false.      ! Modif 17/04/2020
      if (incrATM) then
         facION = 3.d0
         maskGamma(1:nmod) = .false.
         do k=1,nmod
            if (gamma1(k).lt.1.6d0) maskGamma(k) = .true.
         enddo
      endif


*-----------------------------
***   SOLAR MODEL
*_____________________________
      solmod = totm.eq.1.d0.and.nphase.eq.2
      iCEb = novlim(nsconv,3)
      iCEt = novlim(nsconv,4)
      k = 1
      do while (t(k).gt.5.d5)
         k = k+1
      enddo
      iionHe = k-1
      do while (t(k).gt.1.d5)
         k = k+1
      enddo
      iionH = k-1
c      print *,'iCEb',iCEb,'iCEt',iCEt
c      print *,'iionHe',iionHe,'iionH',iionH
      
c     facCE = 3.d0
c      facCE = 5.d1
       facCE = 1.d1

*----------------------------------------------------------------------------------------------------------
***   Microscopic diffusion case
***   Implementation of the criteria for the case of microscopic diffusion without rotation or atmosphere
***   Part 1 : definition of the case       
***   January 2018 - Thibaut Dumont (Genève)
*----------------------------------------------------------------------------------------------------------
  
       diffumod = microdiffus.AND.nphase.EQ.2

       IF (diffumod) PRINT *,
     $      "mesh mode diffusion microscopique actif phase 2"
       

*-----------------------------*
***   NO REZONING : RETURN  ***
*-----------------------------*

      if (maxsh.eq.0) then
         if (dup3) then
            write (nout,1030) nmod
            write (90,1030) nmod
         else
            write (nout,1020) nmod
            write (90,1020) nmod
         endif
c     see parameter card for the value of lnucl
         if (.not.lnucl) then
            forall (i=1:nmod)
               enucl(i) = venucl(i)
               denucldf(i) = 0.d0
               denucldt(i) = 0.d0
               denucldro(i) = 0.d0
            end forall
         endif
         return
      endif



***   initialize variables
*-------------------------

      nmod1 = nmod-1
      meshcut = 0.d0*msun
c..   central values
      lntc = lnt(1)+(m(2)/m(3))*(lnt(1)-lnt(2))
      lnfc = lnf(1)+(m(2)/m(3))*(lnf(1)-lnf(2))
      lnroc = log(ro(1))+(m(2)/m(3))*(log(ro(1)/ro(2)))
      lnpc = log(p(1))+(m(2)/m(3))*(log(p(1)/p(2)))
      venuclc = venucl(1)+(m(2)/m(3))*(venucl(1)-venucl(2))

c..   width of shell 1 is contrained to be no more than 1/10 of the mass
c..   between the center and the base of CB convective zone
      facshmin = 0.05d0

c..   variable needed for zonetest if not diffusion
c      if (nresconv.gt.0.or.(nphase.eq.4.and..not.rotation))
      if (nresconv.gt.0)
     &     then
         do i = 1,nmod
            if (crz(i).lt.0) then
               icrzc(i) = 0
            else
               icrzc(i) = 1
            endif
         enddo
      endif


***   treatment of convective boundaries
*---------------------------------------

c     initialisation rayon coeur, rayon enveloppe conv et zone conv ? (com TD)      
      dmrcore = 1.d-4
      dmrenv = 1.d-4
      dmrcz = 1.d-4
c...  factor by which the resolution is increased at the convective front
      xmresol = 3.d-1
c...  fraction of the mass of the conv. zone outside which the resolution
c     is increased
      xmrconv = 2.d-2
c..   thibaut
c      dmcore = max(0.1d0*mcore,0.2d0*msun)
c      ddmcore = 0.03d0*mcore
c      dmcore = min(0.1d0*mcore,0.2d0*msun)    avant modif Ana
c      ddmcore = min(0.1d0*mcore,0.2d0*msun)   avant modif Ana
C...  ANA 11/09/18 - Test
      dmcore = min(0.01d0*mcore,0.2d0*msun)
      ddmcore = min(0.1d0*mcore,0.2d0*msun)
      nresol = abs(nresconv)
      if (nresol.ne.0) then
c..   at least max(5,nresol) shells on each side of the conv. core
         if (mcore.gt.0.d0) then
c         dmcore = m(novlim(1,4)+nresol)-m(novlim(1,4))
            dmcore = max(dmcore,m(novlim(1,4)+nresol)-m(novlim(1,4)))
            if (novlim(1,4).gt.nresol) ddmcore = max(m(novlim(1,4))
     &           -m(novlim(1,4)-nresol),ddmcore)
         endif
c..  resolution increased in 1xHp about the convective envelope boundary
         if (dup3) then
            k = ienvsh
            rhp = r(k)-max(1.1d0,2.d0*etaturb)*p(k)*r(k)**2/
     &           (ro(k)*g*m(k))
            dmenv = 0.d0
            do while (r(k).gt.rhp.and.k.gt.1)
                k = k-1
               dmenv = dmenv+dm(k)
            enddo
            kdup = k
            dmenv = min(dmenv,0.1d0*dmenvconv)
c            dmenv = min(1.d-5*dtn*msun*seci,0.1d0*dmenvconv)
         else
            dmenv = min(1.d-7*dtn*msun*seci,0.1d0*dmenvconv)
         endif
         ddmenv = dmenv
         if (ndiff_ov.gt.0) dmenv = max(menvconv-m(ienvsh-ndiff_ov-10),
     &        dmenv)
c..  final phase of AGB evolution
         if (agbphase.and.ro(ienvsh).lt.1.d-5) xmresol = 1.d-1
      endif


***   treatment of atmosphere
*----------------------------
c..   at least "natm" shells at the surface between taulim > tau > tau0
      if (ntprof.eq.2) then
         natm = 10
      else
         natm = 100
      endif

      
      if (nretry.eq.1) then
         tdust = 1700.d0
c..   verify if T is outside molecular opacity range T(tau0)>1700K
         if (t(nmod).lt.tdust) then
            write (nout,850) tdust
            do i = nmod,1,-1
               if (t(i).gt.tdust) goto 30
            enddo
 30         tau0 = tau0+0.3d0*(tau(i)-tau0)
         else
            if (tau0.gt.0.005d0.and.t(nmod).gt.2.d3)
     &           tau0 = tau0-0.3d0*(tau(nmod)-0.005d0)
         endif
      endif
      if (ntprof.eq.1.or.ntprof.ge.3) natm = 150
      
      tauinf = max(tau0,1.d-1)
      tausup = max(2.d0*taulim,tauinf,tau(neff-1),1.d1)
      dtaumesh = log(tausup/tauinf)/dble(natm)


***   mesh freeze-out
*--------------------
      lntfreeze = 0.d0
      lntmesh = log(5.d9)
c..   freeze H and He BS after C burning except in super AGB stars
c      if (nphase.ge.6.and.mtini.gt.13.d0) lntfreeze = log(4.d8)
c..   freeze CBS
c      if (lntfreeze.gt.0.d0.and.tmax.gt.1.4d9) lntfreeze = log(1.d9)
c..   freeze NeBS
c      if (lntfreeze.gt.0.d0.and.tmax.gt.1.8d9) lntfreeze = log(1.5d9)


*_____________________________________________________
***
***  DEFINE WHERE THE SHELLS MUST BE ADDED AND DELETED
***
*-----------------------------------------------------
      
c..   during core He burning stabilizes mesh by increasing dnucmax
      if (nphase.eq.4) then
         print *, 'phase 4 reached -> increase dnucmax'
         dnucmax = 1.d50
      else
         dnucmax = dlnenuc
      endif

      do imesh = 1,2
         icount = icount+1
         if (imesh.eq.1) then
            addition = .false.
         else
            addition = .true.
c           myflamespeed concerne le stade evolutif où l'étoile brule le carbone au coeur 
c           print *, 'call myflamespeed'
            call myflamespeed (.true.,imodpr)
         endif
         suppression = .not.addition

c         if (maxsh0.gt.0.and.((suppression.and.ie.ge.maxsh).or.
c     &        (addition.and.io.ge.maxsh))) goto 300

         nmod0 = nmod
         mm = 0

 100     mm = mm+1
         mm1 = mm+1
         mm2 = mm+2
         mm0 = max(mm-1,1)
         meshfrozen = .false.

c         IF (nphase.EQ.2) THEN
c            print *,'mm =',mm, 'm(mm)',m(mm),'dm(mm)',dm(mm)
c            print *, 'ie =',ie, 'io =', io
c            print *, 'addition', addition
c            print *, 'suppression', suppression
c         ENDIF
         
         if (suppression) then
            is = ie
            kp = mm2
            k = mm
            km = mm0
         else
            is = io
            kp = mm1
            k = mm
            km = mm
         endif

***   Freeze grid points
c..   if lnt < lntfreeze
c      if (lnt(mm).lt.lntfreeze) meshfrozen = .true.
c..   if m < meshcut
c     if (m(mm).lt.meshcut) meshfrozen = .true.
         
         IF (meshfrozen) print *, 'goto 300'
         if (meshfrozen) goto 300

         if (mm.ge.ishockb-nq0.and.mm1.le.ishockt+nq0.and.shockdetect)
     &        then
            mm = ishockt+nq0
c            print *, 'goto 100 1'
            goto 100
         endif

***   initialize logicals
         zonetest = .true.
         pvistest = .true.
         courantest = .true.
         duptest = .true.
         tautest = .true.
         structest = .true.
         flametest = .true.
         shellcenter = .true.
***   initialize tolerances
         dvmax = dlnvma
         drvmax = dvmax
         dvmin = dlnvmi
         
c increase resolution for H-gap and bottom of RGB
         if (nphase.eq.3) then
c            print *, 'increase resolution for H-gap and bot of RGB'
            dvmax = dlog(0.5d0)+dlog(dexp(dvmax)+1.d0)
            dvmin = dlog(0.5d0)+dlog(dexp(dvmin)+1.d0)
         endif
         dmrmax = dmrma
         dmrmin = dmrmi

         if (mm.eq.1) then
            print *, 'mm = 1 so dvmax = dlnvmi and dmrmax = dmrmi'
            dvmax = dlnvmi
            drvmax = dvmax
            dmrmax = dmrmi
            print *, 'mm=1 dmrmax =',dmrmax
         endif
       
***   increase resolution in pre-sn stage
         if (lnt(mm).ge.lntmesh) then
            print *, 'increase resolution : pre-sn stage'
            dmrmax = dmrma*xmresol
            dmrmin = dmrmi*xmresol
         endif

***   increase resolution at core/DUP fronts (DUP -> phase RGB ?)
         if (mcore.gt.0.d0.and.nresol.ne.0) then
            if (m(mm).ge.mcore-ddmcore.and.m(mm).le.mcore+dmcore) then
c               print *, 'dmrmax l.569', dmrmax
               dmrmax = min(dmrma*xmresol,dmrcore)
               dmrmin = min(dmrmi*xmresol,dmrcore)
            endif
         endif
         if (menvconv.gt.0.d0.and.nresol.ne.0) then
            if (m(mm).ge.menvconv-dmenv.and.m(mm).le.menvconv+ddmenv)
     &           then
c               drvmax = min(dvmax,log(1.05d0))
c               dvmax = min(dvmax,0.0953d0)
               dmrmax = min(dmrma*xmresol,dmrenv)
               dmrmin = min(dmrmi*xmresol,dmrenv)
            endif
         endif
         if (nconvsh.gt.0.and.nresol.ne.0) then
            do i = 1,nconvsh
               iup = novlim(i,8)
               idwn = novlim(i,7)
***   increase resolution if convective zone large enough
c               print *, 'increase because of conv zone large enough'
               if (iup-idwn.gt.nconvmin) then
                  nrconv = min(int((iup-idwn)/4),10)
                  dmcz = 0.d0
                  dmcz0 = (m(iup)-m(idwn))*xmrconv
                  dmcz = min(m(iup+nrconv)-m(iup),dmcz0)
                  dmcz = max(dmcz,m(iup+2)-m(iup-2))
                  ddmcz = dmcz
                  if (m(mm).ge.m(iup)-ddmcz.and.m(mm).le.m(iup)+dmcz
     &                 .and.mm.lt.iup+nresol) then
                     dmrmax = min(dmrma*xmresol,dmrcz)
                     dmrmin = min(dmrmi*xmresol,dmrcz)
                     goto 200
                  endif
                  if (idwn.gt.nrconv) then
                     dmcz = min(m(idwn)-m(idwn-nrconv),dmcz0)
                  else
                     dmcz = dmcz0
                  endif
                  dmcz = max(dmcz,m(idwn+2)-m(idwn-2))
                  ddmcz = dmcz
                  if (m(mm).ge.m(idwn)-dmcz.and.m(mm).le.m(idwn)+ddmcz
     &                 .and.mm.gt.idwn-nresol) then
                     dmrmax = min(dmrma*xmresol,dmrcz)
                     dmrmin = min(dmrmi*xmresol,dmrcz)
                     goto 200
                  endif
               endif
            enddo
         endif

***   increase resolution in HBS in case of rotational mixing
         if (incrHBS) then
            if (mm.ge.iHBSb.and.mm.le.iHBSt) then
               dmrmax = dmrmax/facHBS
               dmrmin = dmrmin/facHBS
               dvmax = dvmax/facHBS
               drvmax = dvmax
               dvmin = dvmin/facHBS
            endif
         endif

***   increase resolution in core in case of rotational mixing + IGW
         if (incrIGW) then
            if (mm.ge.iIGWb.and.mm.le.iIGWt) then
               dmrmax = dmrmax/facIGW
               dmrmin = dmrmin/facIGW
               dvmax = dvmax/facIGW
               drvmax = dvmax
               dvmin = dvmin/facIGW
            endif
         endif

***   increase resolution in ionization region in case of ATM
         if (incrATM) then
c            print *, 'case of atm'
            if (maskGamma(mm)) then
               dmrmax = dmrmax/facION
               dmrmin = dmrmin/facION
               dvmax = dvmax/facION
               drvmax = dvmax
               dvmin = dvmin/facION
            endif
         endif


***   increase resolution in CE in case of solar model
         if (solmod.AND..NOT.diffumod) then
c            print *,'cas modèle solaire'
c            print *,'solmod in newmesh',solmod,iCEb,iCEt,mm,iionHe-10
c     &           ,iionHe+10,iionH-10,iionH+10
            if (mm.ge.iCEb.and.mm.le.iCEt.and.((mm.gt.iionHe-10.and.mm
     &           .lt.iionHe+10).or.(mm.gt.iionH-10.and.mm.lt.iionH+10)))
     &           then
c            if (mm.ge.iCEb.and.mm.le.iCEt) then
               dmrmax = dmrmax/facCE
               dmrmin = dmrmin/facCE
               dvmax = dvmax/facCE
               dvmin = dvmin/facCE
            endif
         endif

*----------------------------------------------------------------------------------------------------------
***   Microscopic diffusion case
***   Implementation of the criteria for the case of microscopic diffusion without rotation or atmosphere
***   Part 2 : Criterias for the case       
***   January 2018 - Thibaut Dumont (Genève)
*----------------------------------------------------------------------------------------------------------
         IF (diffumod) THEN
c     PRINT *,'cas diffusion atomique'
               IF (mm.EQ.1) print *,'dmrmin ini',dmrmin
               dmrmin = 5.d-3          ! Difference minimale de masse entre deux couches
               drmax = 1.d-2           ! Difference maximale de rayon entre deux couches
               dPtotmax = 2.d-2        ! Difference maximale de pression totale entre deux couches
               dlummax = 2.d-3         ! Difference maximale de luminosité entre deux couches
               dHmax = 5.d-2           ! Difference maximale d'abondance en H1 entre deux couches
               dHe4max = 2.5d-2       ! Difference maximale d'abondance en He4 entre deux couches
               dO16max = 1.d-1  ! Difference maximale d'abondance en O16 entre deux couches
c               dO16max = 1.d-6         ! Difference maximale d'abondance en O16 entre deux couches
               dFemax = 1.d-1          ! Difference maximale d'abondance en Fer entre deux couches
               delemmax = 3.d-1        ! list of element considered here : Ca, Ti, Cr, Mn, Ni -> Pb
               IF (mm.EQ.1.AND.xsp(mm,ih1).LT.1.d-3) dmrmin = 1.d-5 ! condition pour XH_centre 
         ENDIF      
c         dTempmax = 10.d0   ! Difference minimale de température entre deux couches (TD 03/2020)
         
***   treatment of shock regions : impose 1.d-6 =< tcourant =< 1.d-4
 200     preshocktest = hydrotest.and.mr(mm).gt.mshockb.and.mr(mm).lt.
     &        mshockt
         if (preshocktest) then
c         dpvis = abs(log((vpvisc(mm)+1.d0)/(vpvisc(mm1)+1.d0)))
c         pvistest = dpvis.lt.dvmin
c         dmrmax = dmshock
c         dmrmin = min(0.9d0*dmshock,dmrmi)
c      else
c         dmrmax = dmrma
c         dmrmin = dmrmi
            dtcour = (r(mm1)-r(mm))/dsqrt(pw53*p(mm)/ro(mm))
            if (addition) then
               courantest = dtcour.gt.1.d-4
            else
               courantest = dtcour.lt.1.d-6
            endif
         endif

***   increase resolution in ionization regions
         if (mm.eq.1) then
            dxmue = 0.d0
         else
            dxmue = abs(log(vmueinv(kp)/vmueinv(km)))
         endif

***   increase resolution in s-process region (if ishtest = 's')
         dxnx = 0.d0
         if (ishtest.eq.'s'.and.m(mm).ge.mdredgeup-1.d-3*msun.and.
     &        m(mm).le.mdredgeup+1.d-4*msun.and..not.thermalpulse) then
            if (max(xsp(kp,ic12),xsp(km,ic12)).gt.1.d-6.or.
     &           max(xsp(kp,iH1),xsp(km,iH1)).gt.1.d-8) then
               dxnx = max(abs(log(xsp(kp,ic12)/xsp(km,ic12))),
     &              abs(log(xsp(kp,ih1)/xsp(km,ih1))))
            else
               dxnx = 0.d0
            endif
         endif

***   treatment of the structure
         dror = 0.d0                ! density
         dpnp = 0.d0                ! pressure
         dlnl = 0.d0                ! luminosity
         denn = 0.d0                ! nuclear energy
         drr = 0.d0                 ! radius
         dfnf = 0.d0                ! e- degen (taux)
         dtnt = 0.d0                ! temps (?)
         dlmnm = 1.d-10             ! mass of a shell
         dtaur = -1.d99
         ddmnm = 0.d0


***   treatment of centered variables
*-------------------------------------
c         print *, 'treatment of centered variables'
         if (mm.eq.1) then
            ddmnm = abs(log(dm(kp)/dm(1)))
            dfnf = abs(lnf(kp)-lnfc)
            dtnt = abs(lnt(kp)-lntc)
            dror = abs(log(ro(kp))-lnroc)
            dpnp = abs(log(p(kp))-lnpc)
            if (abs(venuclc).lt.enuclim.and..not.flame) then
               denn = 0.d0
            else
               denn = abs(log(abs(venucl(kp)/venuclc)))
            endif
         else
	    if (addition) then
               ddmnm = abs(log(2.d0*min(dm(min(kp,nmod-1)),dm(km-1))/
     &              dm(km)))
            else
               ddmnm = abs(log(min(dm(min(kp,nmod-1)),dm(km))/
     &              (dm(k)+dm(k+1))))
            endif
            dfnf = abs(lnf(kp)-lnf(km))
            dtnt = abs(lnt(kp)-lnt(km))
            dror = abs(log(ro(kp)/ro(km)))
            dpnp = abs(log(p(kp)/p(km)))
c            print *, 'kp',kp,'km',km
c         if (max(abs(venucl(kp)),abs(venucl(km))).lt.enuclim.or.
c     &        dtnt.lt.0.01d0*dvmin) then
            if (max(abs(venucl(kp)),abs(venucl(km))).lt.enuclim) then
               denn = 0.d0
            else
               denn = abs(log(abs(venucl(kp)/venucl(km))))
            endif
         endif


***   treatment of interface variables
*-------------------------------------
c         print *, 'treatment of interface variables'
         dmnm = mr(kp)-mr(k)
c         if (addition) dmnm = 2.d0*dmnm
c..   width of central shell is constrained based on the distance from the
c..   center to the first convective zone
         if (mm.eq.1) then
            if (nconvsh.gt.0.and.novlim(1,7).gt.2) then
               if (m(3).gt.(m(novlim(1,7))*facshmin))
     &              shellcenter = .false.
            endif
c..   impose central shell to be less massive than shell number 2
            if (suppression.and.m(2).gt.1.d-6*msun)
     &           shellcenter = .false.
         else
            drr = lnr(kp)-lnr(k)
            dlmnm = abs(log(m(kp)/m(k)))
         endif
c..   the luminosity profile is constrained only in regions where the
c..   temperature profile also shows substantial variations
cc257   if ((mm.eq.1.or.(abs(lum(mm1)).lt.lumlim.and.denn.eq.0.d0).or.
cc     &     dtnt.lt.1.d-2*dvmin.or.lnt(mm).lt.lntmax)) then
cc      if (mm.eq.1.or.abs(lum(mm1)).lt.lumlim.or.dtnt.lt.1.d-2*dvmin)
cc     &     then
c276         if (mm.eq.1.or.(denn.eq.0.d0.and.(max(abs(lum(kp)),abs(lum(k)))
c276     &        .lt.lumlim.or.dtnt.lt.1.d-2*dvmin))) then
         if (mm.eq.1.or.max(abs(lum(kp)),abs(lum(k))).lt.lumlim.or.
     &        dtnt.lt.1.d-2*dvmin) then
            dlnl = 0.d0
         else
            dlnl = abs(log(abs(lum(kp)/lum(k))))
         endif

         if (k.gt.100) dtaur = log(tau(k)/tau(kp))


***   treatment of rotation
*--------------------------
         if (rotation) then
            domega = abs(log(abs(vomega(kp)/vomega(k))))
         else
            domega = 0.d0
         endif


***   define zonetest (adjacent shells identical = .true.)
*---------------------------------------------------------

c..   addition : because split most massive shell check, k-1,k,k+1 and k+2
c         if (mm.gt.1.and.mm.lt.nmod0-2.and.(nresconv.gt.0.or.
c     &        (nphase.eq.4.and..not.rotation))) then
         if (mm.gt.1.and.mm.lt.nmod0-2.and.nresconv.gt.0) then
            if (addition) then
               icz = abs(icrzc(mm)+icrzc(mm1)-1)
               icz1 = abs(icrzc(mm+1)+icrzc(mm+2)-1)
               icz0 = abs(icrzc(mm)+icrzc(mm-1)-1)
               zonetest = icz.eq.1.and.icz1.eq.1.and.icz0.eq.1
c               print *,'zonetest definie'
            else
               icz = abs(icrzc(mm)+icrzc(mm1)-1)
               zonetest = icz.eq.1
            endif
c            zonetest = zonetest.or.idiffcc
            pass = .true.
         endif


***   treatment of flames
*------------------------

         if (imodpr.eq.11.or.imodpr.eq.12) then
c..   shell suppression in the precursor flame region forbidden
            if (suppression) then
               do i = 1,nbcz
                  if (mm.le.(ibCOcz(i)+2).and.mm.ge.iminlumprof(i).and.
     &                 iminlumprof(i).ne.0) then
                     flametest = .false.
                     exit
                  endif
               enddo
            else
               do i = 1,nbcz
c..   skip this shell so that we can add one under CO convective zone
                  if (mm.eq.(ibCOcz(i)-2)) print *, 'goto 100 2'
                  if (mm.eq.(ibCOcz(i)-2)) goto 100
c..   add one under CO convective zone if necessary
                  if (mm.eq.(ibCOcz(i)-1)) then
                     if ((r(ibCOcz(i))-r(ibCOcz(i)-1)).gt.
     &                    (0.5d0*thurff(i)*dtn)) then
                        flametest = .true.
                     else
                        flametest = .false.
                     endif
                     exit
                  endif
               enddo
            endif
         endif


*_________________________________________________
***
***  PERFORM FINAL TEST FOR SHELL ADDITION/REMOVAL
***
*-------------------------------------------------

         dlmnm = 0.d0
         dxmue = 0.d0
         fdm = 2.4d0
c..  do not check mass distribution in uppermost atmospheric layers
         if (tau(mm).lt.tausup.and.tau(mm).gt.tau0) ddmnm = 0.d0
         center = dfnf.lt.dvmax.and.dpnp.lt.dvmax.and.dtnt.lt.dvmax
     &        .and.dror.lt.dvmax.and.denn.lt.dnucmax.and.pvistest
     &        .and.ddmnm.lt.fdm.and.dxnx.lt.log(1.1d0)
     &        .and.dxmue.lt.dvmax*0.4d0
   
         
c..  ddmnm impose that the ratio of the mass of 2 consecutive shells
c..  does not differ by a more than a factor f=4 (1.61=ln(1+f))
         edge = dmnm.lt.dmrmax.and.dlnl.lt.dvmax.and.drr.lt.drvmax
     &        .and.dlmnm.lt.dvmax.and.domega.lt.dvmax.and.shellcenter
     &        .and.courantest

         IF (diffumod.AND.crz(mm).EQ.4) THEN

            chimi = ((xsp((mm+1),ih1)-xsp(mm,ih1)).LT.dHmax
     &           .AND.(xsp((mm+1),ihe4)-xsp(mm,ihe4)).LT.dHe4max
     &           .AND.(xsp((mm+1),20)-xsp(mm,20)).LT.dO16max
     &           .AND.(xsp((mm+1),53)-xsp(mm,53)).LT.dFemax)
         ENDIF


         IF (diffumod) structest = center.AND.edge.AND.chimi
         IF (.NOT.diffumod) structest = center.AND.edge

        
         if (tau(mm).lt.tausup.and.tau(mm).ge.tauinf)
     &        tautest = dtaur.lt.dtaumesh
         meshcrit = (structest.and.tautest.and.duptest.and.zonetest
     &        .and..not.meshfrozen.and.flametest).or.dmnm.lt.1.d-12

c     Si suppression, test contraintes et si respect on supprime une couche
c     Si addition, test contraintes et si non respect on ajoute une couche
         if (addition) meshcrit = .not.meshcrit                  ! Inverse meshcrit (devient faux si était vraie et devient vraie si était faux)
         
         if (debug.and.mm.gt.960.and.mm.lt.980) then
            print *,mm,meshcrit,addition,structest,tautest,duptest
     &           ,zonetest
            if (.not.center) print *,'c',center,dfnf.lt.dvmax,dpnp.lt
     &           .dvmax,dtnt.lt.dvmax,dror.lt.dvmax,denn.lt.dnucmax
     &           ,ddmnm.lt.1.61d0,ddmnm,tau(mm),tauinf,dvmax
            if (.not.edge) print *,'e',mm,edge,dmnm.lt.dmrmax
     &           ,dlnl.lt.dvmax
              print *,'e',edge,dmnm.lt.dmrmax,dlnl.lt.dvmax
     &           ,drr.lt.dvmax,dlmnm.lt.dvmax,domega.lt.dvmax
     &           ,shellcenter,dlmnm,dvmax
         endif
         if (meshcrit) then
c..   split the most massive shell
            if (addition.and.dm(kp).gt.dm(k).and..not.center
     &           .and.mm.lt.nmod1) mm = mm+1

            if (suppression) then
               ie = ie+1
               kshell(mm) = kshell(mm)-1
               nshdel(ie) = mm
            else
               kshell(mm) = kshell(mm)+1
               if (mm.gt.1) then
                  if (kshell(mm-1).eq.-1.and.dm(mm-1).eq.dm(mm)) then
                     kshell(mm-1) = 0
                     kshell(mm) = 0
                     goto 300
                  endif
               endif
	       if (nmod0-ie+io+2.gt.nsh) then
		  mm = nsh 
		  goto 300
 	       endif
               io = io+1
               nshadd(io) = mm
            endif
c..   avoid the removal/addition of 2 consecutive shells
c            if (tau(mm).gt.taulim) mm = mm+2
            mm = mm+1
         endif

c         if ((suppression.and.mm+2.le.nmod1).or.(addition.and.mm+1.lt.
c     &        nmod0).and.mm+2.lt.nsh) print *, 'goto 100 3'
         
 300     if ((suppression.and.mm+2.le.nmod1).or.(addition.and.mm+1.lt.
     &        nmod0).and.mm+2.lt.nsh) goto 100

      enddo
      if (icount.ne.2) print *,'WARNING : PROBLEM WITH IMESH'

      print *, 'ie',ie, 'io', io

***____________________________________
***
***          SHELL SUPPRESSION
***
***         remove shells il+1 = nshdel(m)+1
***     [il,il+1,il+2] --> [il,il+2]
***     [nshdel(m),nshdel(m)+1,nshdel(m)+2] --> [nshdel(m),nshdel(m)+2]
***   only centered variable re-defined
***------------------------------------


 500  iremove = ie
      i = 1
*   discard shells that were both assigned for addition and suppression
      do while (i.le.iremove)
         if (kshell(nshdel(i)).ne.-1) then
            do k = i,ie-1
               nshdel(k) = nshdel(k+1)
            enddo
            nshdel(iremove) = 0
            iremove = iremove-1
         else
            i = i+1
         endif
      enddo
      ie = min(iremove,maxsh)
c      ie = iremove
      iadd = io
      i = 1
      do while (i.le.iadd)    
         if (kshell(nshadd(i)).ne.1) then
            do k = i,io-1
               nshadd(k) = nshadd(k+1)
            enddo
            nshadd(iadd) = 0
            iadd = iadd-1
         else
            i = i+1
         endif
      enddo
      io = min(iadd,maxsh)
c      io = iadd
      do k = 1,io
         nshad(k) = nshadd(k)
      enddo
      do i = 1,ie
         if (kshell(nshdel(i)).ne.-1) write (nout,900) nshdel(i)
         il = nshdel(i)-i+1
         ir = il+1
*     interpolate variables in merged shells
         call interpmesh (il,ir,1)
*     redefinition of the shell mass
         dm(il) = dm(il)+dm(ir)   ! la couche restante obtient la masse de la couche supprimée
         m(ir+1) = m(il)+dm(il)
         mr(ir+1) = m(ir+1)/(totm*msun)
         nmod = nmod-1            ! le nombre de couche diminue d'une unité
*     shell shifting
         do k = 1,io
            if (nshad(k).gt.il) nshad(k) =  nshad(k)-1
         enddo
         do k = ir,nmod
            call interpmesh (k,k+1,0)
         enddo
         if (il.lt.kshockb) kshockb = kshockb-1
         if (il.lt.kshockt) kshockt = kshockt-1
*     redefine convective borders
         nconvsh = nsconv
         do k = 1,nconvsh
            if (il.le.novlim(k,3)) then
               do j = k,nconvsh
                  novlim(j,7) = novlim(j,7)-1
                  novlim(j,8) = novlim(j,8)-1
                  novlim(j,3) = novlim(j,3)-1
                  novlim(j,4) = novlim(j,4)-1
               enddo
               goto 550
            endif
            if (il.le.novlim(k,4)) then
               novlim(k,8) = novlim(k,8)-1
               novlim(k,4) = novlim(k,4)-1
               if (k.lt.nconvsh) then
                  do j = k+1,nconvsh
                     novlim(j,7) = novlim(j,7)-1
                     novlim(j,8) = novlim(j,8)-1
                     novlim(j,3) = novlim(j,3)-1
                     novlim(j,4) = novlim(j,4)-1
                  enddo
               endif
               goto 550
            endif
         enddo
 550     mr(nmod) = 1.d0
         dm(nmod) = 0.d0
      enddo
      nmod1 = nmod-1     
      
***---------------------------
***
***      SHELL ADDITION
***
***     max(dm[il],dm[il+1])
***        divided by 2
*-----------------------------

 700  iadd = min((io+1),maxsh+1)
c 700  iadd = io+1                          
c      if (io.EQ.0) print *, 'goto 800'
      if (io.eq.0) goto 800
      if (nshad(io).eq.0) then
         do k = 1,io
            nshad(k) = nshadd(k)
         enddo
      endif
      nshad(iadd) = nmod
*     redefinition of the mass grid : shell shifting
      do i = iadd,2,-1
         il = nshad(i-1)
         ir = nshad(i-1)+i-1
         if (il.lt.kshockb) kshockb = kshockb+1
         if (il.lt.kshockt) kshockt = kshockt+1
*     shell shifting
         do k = nshad(i)+i-1,nshad(i-1)+i-1,-1
            call interpmesh (k,k-i+1,0)
         enddo
         nmod = nmod+1          ! le nombre de couche augmente d'une unité
         dm(il) = dm(il)*0.5d0
         dm(ir) = dm(il)           ! partage de la masse de la couche initiale en 2 entre celle-ci et la nouvelle
         m(ir) = m(il)+dm(il)
         mr(ir) = m(ir)/(totm*msun)
      enddo
      nmod1 = nmod-1
*     interpolate variables in added shells
 800  if (shockdetect) then
         do i = 1,iadd-1
            il = nshad(i)+i-1
            ic = il+1
            ir = ic+1
            if (il.eq.1) then
               call interpmesh (il,ir,7)
            else
               if (nshadd(i).eq.ishockb) call interpmesh (il,ir,5)
               if (nshadd(i).lt.ishockb) call interpmesh (il,ir,2)
               if (nshadd(i).eq.ishockt) call interpmesh (il,ir,6)
               if (nshadd(i).gt.ishockt) call interpmesh (il,ir,2)
            endif
         enddo
      else
         do i = 1,iadd-1
            il = nshad(i)+i-1
            ic = il+1
            ir = ic+1
            if (tau(ir).lt.taulim.or.ir.ge.nmod1) then
               call interpmesh (il,ir,5)
            else
               if (il.eq.1) then
                  call interpmesh (il,ir,7)
               else
                  call interpmesh (il,ir,2)
               endif
            endif
         enddo
      endif
      ishockb = kshockb
      ishockt = kshockt

      pvisc(1) = 0.d0
      pvisc(nmod) = 0.d0
      r(1) = 1.d-99
      vr(1) = 1.d-99
      forall (i=1:nmod)
         lnr(i) = log(r(i))
         vlnr(i) = log(vr(i))
      end forall
      if (.not.lnucl) then
         forall (i=1:nmod)
            enucl(i) = venucl(i)
            denucldf(i) = 0.d0
            denucldt(i) = 0.d0
            denucldro(i) = 0.d0
         end forall
      endif

      if (ireset.ne.0) return


***  screen outputs
*------------------

      if (model.eq.modeli.and.iter.eq.1.and.ifail.eq.0) then
         if (pass) then
            write (90,*) '    *** MESH : zonetest activated'
         else
            write (90,*) '    *** MESH : zonetest NOT activated'
         endif
         if (dup3) then
            print *, 'dup3 active'
            write (90,1010) ie,io,nmod
         else
            write (90,1000) ie,io,nmod
         endif
         if (ie.ge.1) then
            write (90,1100) (nshdel(k),k = 1,min(ie,16))
            write (90,1200) (nshdel(k),k = 17,ie)
         endif
         if (io.ge.1) then
            write (90,1300) (nshadd(k),k = 1,min(io,16))
            write (90,1400) (nshadd(k),k = 17,io)
         endif
         write (90,*)
      endif
      if (dup3) then
         write (nout,1010) ie,io,nmod
      else
         write (nout,1000) ie,io,nmod
      endif
      if (ie.ge.1) then
         write (nout,1100) (nshdel(k)+1,k = 1,min(ie,16))
         write (nout,1200) (nshdel(k)+1,k = 17,ie)
      endif
      if (io.ge.1) then
         write (nout,1300) (nshadd(k),k = 1,min(io,16))
         write (nout,1400) (nshadd(k),k = 17,io)
      endif
      write (nout,*)

c..   test mesh
c      call resulpr (0)
c      stop
      
 850  format (3x,'WARNING : T(tau0) < Tdust = ',0pf5.0)
 900  format (3x,'WARNING : mesh logic problem at shell # ',i4)
 1000 format (5x,'suppressed shells: - ',i3,'; added shells: + ',i3,
     &     ' [shell number = ',i4,']')
 1010 format (5x,'suppressed shells: - ',i3,'; added shells: + ',i3,
     &     ' [shell number = ',i4,']; 3DUP phase')
 1020 format (1x,'no rezoning applied, nshell = ',i4,/)
 1030 format (1x,'3DUP phase, no rezoning applied, nshell = ',i4,/)
 1100 format (5x,'suppressed shells: ',16(1x,i4,1x))
 1200 format (20(1x,i4,1x))
 1300 format (5x,'added shells: ',5x,16(1x,i4,1x))
 1400 format (20(1x,i4,1x))


      return
      end
