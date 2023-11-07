      SUBROUTINE rmodpar

************************************************************************
*     Read parameter card and models properties                        *
*     Modifs CC ondes (23/11/07)                                       *
* $LastChangedDate:: 2018-01-22 10:25:40 +0000 (Mon, 22 Jan 2018)    $ *
* $Author:: amard                                                    $ *
* $Rev:: 114                                                         $ *
*                                                                      *
************************************************************************
    
      implicit none

************************************************************************
    
      include 'evolpar.star'
    
      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.eos'
      include 'evolcom.ion'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.opa'
      include 'evolcom.rot'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.transp'
      include 'evolcom.var'
C     Modifs CC ondes (2/04/07)
      include 'evolcom.igw'

************************************************************************
      
      integer imerk,icrasht,ldum,ihydro
      integer i,j,k,l

      double precision disctime,dlnvma0,dlnvmi0,dlnenuc0
      double precision dm_ext,theta
      double precision etaturby,dtiny,dtminy,dtmaxy,massend0,massrate0
      double precision paq_c,paq_d
      double precision evolpar_version

      character chydro*1
      character cdiff*18
      character cover*36,machine*20
      character month*36,GetDateTimeStr*15,fmt*26
      parameter (month = 'JanFebMarAprMayJunJulAugSepOctNovDec')
      parameter (fmt = '(I2.2,A1,I2.2,I3,A3,I4)')
      integer itime(8), n_turbul(1)  ! modif TD AP
      double precision Tfix, om_turbul, PM ! modif TD AP
      double precision tol_u,tol_lnr,tol_lnf,tol_lnT,tol_l,tol_ang

      logical vlconv
      logical locconsAM

      common /conserv/ locconsAM
      common /disclocking/ disctime
      common /fitcd/ paq_c(3,8,50),paq_d(8,50)
      common /resolution/ vlconv
      common /convergence/ imerk,icrasht
      common /turbulparam/ n_turbul,Tfix,om_turbul,PM

      namelist /RICHER/ Tfix, om_turbul, n_turbul
      namelist /PROFITTMICHAUD/ PM
      namelist /MassiveMloss / clumpfac,zscaling

      namelist /starevol_parameters/ evolpar_version,
     &     maxmod,imodpr,mtini,zkint,
     &     addH2,addHm,tmaxioH,tmaxioHe,
     &     ionizHe,ionizC,ionizN,ionizO,ionizNe,
     &     lgrav,lnucl,
     &     alpha_mlt_hp,alpha_mlt_hrho,dtnmix,mixopt,hpmlt,nczm,
     &     ihro,iturb,etaturby,alphatu,
     &     novopt,aup,adwn,idup,
     &     tau0,ntprof,taulim,
     &     mlp,etapar,dmlinc,
     &     nucreg,nuclopt,tolnuc,ftnuc,znetmax,
     &     idiffcc,idiffnuc,idiffty,diffzc,diffst,
     &     zgradmu,grad_crit,del_zc,ledouxmlt,
     &     omegaconv,Dh_prescr,dm_ext,
     &     Dv_prescr,om_sat,D_zc,disctime,thermal_equilibrium,
     &     diffvr,idiffvr,idiffex,diffcst,breaktime,
     &     nmixd,itmind,itmaxd,rprecd,nretry,
     &     lover,lthal,ltach,lmicro,
     &     iaccr,accphase,itacc,massrate,massend,menv,
     &     ric,prrc,xiaccr,fdtacc,alphaccr,
     &     chydro,ivisc,q0,mcut,
     &     maxsh,nresconv,
     &     dlnvma,dlnvmi,dlnenuc,dmrma,dmrmi,
     &     dtiny,dtminy,dtmaxy,facdt,fkhdt,
     &     ishtest,fts,ftsh,ftshe,ftsc,ftacc,ftst,
     &     itmin,itermax,itermix,icrash,numeric,icorr,
     &     tol_u,tol_lnr,tol_lnf,tol_lnT,tol_l,tol_ang,
     &     phi,alpha,iacc,alphmin,alphmax,sigma,vlconv

************************************************************************
      
***   Read parameter card
      
      open(UNIT=79, FILE='starevol.par', STATUS='OLD')
      rewind(79)
      read(79, NML=starevol_parameters) 
      close(79)

      mlp0 = mlp
      massrate0 = massrate
      massend0 = massend
      maxsh0 = maxsh
      dlnvma0 = dlnvma
      dlnvmi0 = dlnvmi
      dlnenuc0 = dlnenuc
      icrash0 = icrash
      eps(1) = tol_u
      eps(2) = tol_lnr
      eps(3) = tol_lnf
      eps(4) = tol_lnT
      eps(5) = tol_l
      epsa = tol_ang

***   Checks
      
      if (maxmod.eq.0) then
         write (nout,*) 'maxmod = 0 : program stopped'
         stop ' stop maxmod = 0'
      endif

      if (ishtest.eq.'T') ishtest = 't'
      if (ishtest.eq.'S') ishtest = 's'
      if (omegaconv.eq.'T') omegaconv = 't'
      if (omegaconv.eq.'F') omegaconv = 'f'
      if (omegaconv.eq.'S') omegaconv = 's'
      if (omegaconv.eq.'M') omegaconv = 'm'
      omegaconv0 = omegaconv
      mlp = abs(mlp0)
      if (mlp.gt.49) mlp = mlp-50
      om_sat=dble(om_sat)

***   Check versions compatibility

      if(code_version.lt.evolpar_version.and.evolpar_version.lt.2.1d0)
     &     then
         write (nout,*) 'starevol.par not compatible with executable :',
     &        ' network problem ',code_version,evolpar_version
         stop ' stop rmodpar'
      endif

***   Check parameters

      if (lgrav.lt.0.or.lgrav.gt.5) then
         write (nout,*) 'choose a correct value for lgrav (0-4) !'
         stop 'rmodpar : bad lgrav'
      endif
      if (mlp.ne.0.and.iaccr.ne.0) then
         write (nout,*) 'choose between accretion and mass loss !'
         stop 'rmodpar : bad mlp/iaccr'
      endif
      if (hpmlt.lt.0.or.hpmlt.gt.6) then
         write (nout,*) 'choose a correct value for hpmlt (1-5) !'
         stop 'rmodpar : bad hpmlt'
      endif
      if (hpmlt.eq.5.and.alphatu.lt.alpha_mlt_hp) then
         write (nout,*) 'when HMLT = 5, alphatu MUST be > '//
     &        'alpha_mlt_hp !'
         stop 'rmodpar : bad alphatu and alpha_mlt_hp'
      endif
      if (.not.(mlp.le.28.or.mlp.ne.55.or.mlp.ne.56)) then
         write (nout,*) 'choose a correct value for mlp (0-22) !'
         stop 'rmodpar : bad mlp'
      endif
      if (ivisc.lt.0.or.ivisc.gt.3) then
         write (nout,*) 'choose a correct value for ivisc (0-3) !'
         stop 'rmodpar : bad ivisc'
      endif
      if (nucreg.lt.0.or.nucreg.gt.3) then
         write (nout,*) 'choose a correct value for nucreg (0-3) !'
         stop 'rmodpar : bad nucreg'
      endif
      if (nmixd.lt.0.or.nmixd.gt.5) then
         write (nout,*) 'choose a correct value for nmixd (0-5) !'
         stop 'rmodpar : bad nmixd'
      endif
      if (dtnmix.ne.'f'.and.dtnmix.ne.'F'.and.dtnmix.ne.'t'.and.
     &     dtnmix.ne.'T'.and.dtnmix.ne.'g'.and.dtnmix.ne.'G'.and.
     &     dtnmix.ne.'u'.and.dtnmix.ne.'U') then
         write (nout,*) 'choose correct value for dtnmix (f,t,g or u)'
         stop 'rmodpar : bad dtnmix'
      endif
c        if (idiffvr.eq.6..or.idiffvr.eq.7.and.omegaconv.ne.'f') then
      if (idiffvr.eq.6.and.omegaconv.ne.'f') then
         write (nout,*) 'imcompatible choice of idiffvr and omegaconv'
         stop 'rmodpar : bad idiffvr'
      endif

***   Read massloss.par

      if (mlp.eq.15.or.mlp.eq.16) then
        open(UNIT=79, FILE='massloss.par', STATUS='OLD')
        rewind(79)
        read(79, NML=MassiveMloss) 
        close(79)
      endif

      
***   define mixing type

c..   if 0 < nmixd < 5, you take control of the treatment of the convective
c     zones. Therefore you impose  diffzc = .false.

c..   if rotational mixing activated, microscopic diffusion for stars only
c..   in the mass range 1.4 < M < 2 (Li dip thing)
c      if (lmicro.eq.3.and.idiffcc.and.(idiffty.eq.11.or.idiffty.eq.13
c     &     .or.idiffty.eq.14).and.(totm.gt.1.4d0.or.totm.gt.2.d0))
c     &     lmicro = 0
      if (nmixd.gt.0.and.nmixd.lt.5) diffzc = .false.
      diffover = idiffcc.and.((lover.ge.22.and.lover.le.39).or.
     &     (lover.ge.70.and.lover.le.73))
      if (novopt.gt.0.and.diffover) then
         write (nout,*) 'INCOMPATIBLE PARAMETERS : overshooting ',
     &        'treated EITHER by novopt OR by idiffty, make a choice !'
         stop
      endif
      if (diffover.and.(lover.eq.60.or.lover.eq.61.or.lover.eq.70.or.
     $     lover.eq.71)) then
         open(UNIT=70, FILE='transport_param.par', STATUS='OLD')
         rewind(70)
         read(70, NML=Richer) 
         read(70, NML=ProfittMichaud) 
         close(70)
      endif
      
      semiconv = idiffcc.and.ledouxmlt
      if (ledouxmlt.and..not.diffzc) then
         write (nout,7300)
         write (90,7300)
      endif
      rotation = idiffcc.and.(idiffty.ge.8.and.idiffty.le.17)                 ! TD ajout 15 (03/2020)
      difftacho = idiffcc.and.(ltach.ge.41.and.ltach.le.44)
      turbulence = idiffcc.and.((lover.ge.60.and.lover.le.61).or.
     &     (lover.ge.70.and.lover.le.73))
      diffusconv = diffzc.or.diffover.or.semiconv
      microdiffus = idiffcc.and.(lmicro.ge.2.and.lmicro.le.6) ! modif thoul = 4
      viscosityadd = idiffcc.and.idiffty.eq.15 ! TD pour nu_add (03/2020)
      if (viscosityadd.and.grad_crit.lt.1) then ! TD pour nu_add (03/2020)
         write (nout,*)
     $        'choose a correct value for additional viscosity'
         write (nout,*) 'nuadd = 3.5e4 is an exemple of possible value'
         stop
      endif
      if (idiffnuc.eq.3.and..not.diffusconv) then
         write (nout,*) 'choose a correct value for idiffnuc/diffusconv'
         write (nout,*) 'idiffnuc.eq.3 not compatible with diffzc = f',
     &        ', semiconv = f and diffover = f'
         stop
      endif
      if (nuclopt.eq.'m'.and.diffusconv) then
         write (nout,*) 'INCOMPATIBLE PARAMETERS : mixing in ',
     &        'convective zones'
         write (nout,*) 'if you set nuclopt = m, diffzc must be false !'
         stop
      endif
      locconsAM = idiffcc.and.idiffty.eq.16
      if (locconsAM) rotation = .true.

      if (rotation) thermal_equilibrium = .true.
      thermohaline = idiffcc.and.(lthal.eq.1.or.lthal.eq.2.
     &   or.lthal.eq.3.or.lthal.eq.4)
C Modifs CC ondes (12/10/07) -->
C logicalondes : test logique pour le calcul du transport de moment
C                cinetique par les ondes
C logicalondeschim : test logique pour le calcul du transport des
C                especes chimiques par les ondes
C      logicalondes = .true.
C      logicalondeschim = .true.
      igwsurfchim = idiffcc.and.
     &              ((lover.eq.41).or.(lover.eq.43))
      igwcorechim = idiffcc.and.
     &              ((lover.eq.42).or.(lover.eq.43)) ! pas de prescription pour l'instant
      igwsurfrot = idiffcc.and.rotation.and.
     &              ((lover.eq.41).or.(lover.eq.43))
      igwcorerot = idiffcc.and.rotation.and.
     &              ((lover.eq.42).or.(lover.eq.43))
      igwrot = igwsurfrot.or.igwcorerot
      igwsurfenerg = .false.
      igwcoreenerg = .false.
      igw = igwsurfchim.or.igwcorechim.or.igwrot
c      print *,'igw,igwsurfchim,igwrot,igwsurfenerg, = ',igw,igwsurfchim,
c     &         igwrot,igwsurfenerg
      if (igw.and.omegaconv.eq.'t') then
         write (nout,*) "Bad parameter choice: conflict omegaconv - igw"
c         stop "Bad parameter choice: conflict omegaconv - igw"
      endif
C <--

C Check value for Dh (rotation mixing)
      if (.not.(Dh_prescr.eq.'Zahn1992'.or.Dh_prescr.eq.'MPZ_2004'.or.
     &     Dh_prescr.eq.'Maeder03'.or.Dh_prescr.eq.'Maeder06'.or.
     &     Dh_prescr.eq.'Mathis16'.or.Dh_prescr.eq.'Mathis02'.or.
     &     Dh_prescr.eq.'Mathiepi'))
     &     stop 'Bad value for Dh'

***   Initialize/setup parameters

      if (q0.le.0.d0) ivisc = 0
      if (ivisc.eq.0) mcut = 0.d0
      icrasht = 0
      mcut = mcut*msun
c      etaturb = etaturby*dble(iturb)
      etaturb = etaturby
      dtin = dtiny*sec
      dtmax = dtmaxy*sec
      dtmax0 = dtmax
      dtmin = dtminy*sec
      phi0 = phi
      nq0 = int(q0)+imodpr
      maxmod0 = maxmod
      alpha0 = alpha
      imodpr0 = imodpr
      maxsh = abs(maxsh0)
      nmixd0 = nmixd
      mixopt0 = mixopt
      idiffcc0 = idiffcc
      nuclopt0 = nuclopt
      if (dtnmix.eq.'g'.or.dtnmix.eq.'u') then
         chkmix = .true.
      else
         chkmix = .false.
      endif
      if (icorr.eq.'T') icorr = 't'
      if (icorr.eq.'F') icorr = 'f'
      if (icorr.eq.'I') icorr = 'i'
      if (icorr.eq.'M') icorr = 'm'
      if (icorr.eq.'H') icorr = 'h'
      if (icorr.eq.'A') icorr = 'a'
      if (nuclopt.eq.'u') then
         urca = .true.
      else
         urca = .false.
      endif
      icorr0 = icorr
      dlnvma = dlog(dlnvma0+1.d0)
      dlnvmi = dlog(dlnvmi0+1.d0)
      dlnenuc = dlog(dlnenuc0+1.d0)
      ireverse = 1
      if (icrash0.lt.0) ireverse = -1
      icrash = abs(icrash0)
      massend = massend0*msun
      ric = -ric
      if (iaccr.eq.0.and.mlp.ne.17.and.mlp.ne.18) then
         massrate = 0.d0
      else
         massrate = massrate0
      endif
      irotbin = 0
      iacc0 = iacc
      if (rotation) irotbin = 1
c     if (irotbin.eq.1) vsurf = diffvr*1.d5
      if (vlconv) mixopt = .false.
      do l = 1,5
         eps0(l) = eps(l)
      enddo
      call date_and_time (values=itime)
      write (GetDateTimeStr,FMT) itime(5),':',itime(6),itime(3),
     &     month(3*itime(2)-2:3*itime(2)),itime(1)
      call getenv('HOST',machine)

      write (nout,4100) GetDateTimeStr,machine,code_version,code_version
      write (90,4100) GetDateTimeStr,machine,code_version,code_version
      write (90,4200) maxmod,imodpr,mtini,zkint
      write (90,4300) addH2,addHm,tmaxioH,tmaxioHe
      write (90,4400) ionizHe,ionizC,ionizN,ionizO,ionizNe
      write (90,4500) lgrav,lnucl
      write (90,4600) alpha_mlt_hp,dtnmix,mixopt,hpmlt,nczm
      write (90,4700) ihro,iturb,etaturby,alphatu
      write (90,4800) novopt,aup,adwn,idup
      write (90,4900) tau0,ntprof,taulim
      write (90,5000) mlp0,etapar,dmlinc
      write (90,5100) nucreg,nuclopt,tolnuc,ftnuc,znetmax
      write (90,5200) idiffcc,idiffnuc,idiffty,diffzc,diffst
      write (90,5300) zgradmu,grad_crit,del_zc,ledouxmlt
      write (90,5310) omegaconv,Dh_prescr,dm_ext
      write (90,5320) Dv_prescr,om_sat,D_zc,disctime,thermal_equilibrium
      write (90,5330) diffvr,idiffvr,idiffex,diffcst,breaktime
      write (90,5400) nmixd,itmind,itmaxd,rprecd,nretry
      write (90,5500) lover,lthal,ltach,lmicro
      write (90,5600) iaccr,accphase,itacc,massrate0,massend0,menv
      write (90,5700) abs(ric),prrc,xiaccr,fdtacc,alphaccr
      write (90,5800) chydro,ivisc,q0,mcut/msun
      write (90,5900) maxsh0,nresconv
      write (90,6000) dlnvma0,dlnvmi0,dlnenuc0,dmrma,dmrmi
      write (90,6100) dtiny,dtminy,dtmaxy,facdt,fkhdt
      write (90,6200) ishtest,fts,ftsh,ftshe,ftsc,ftacc,ftst
      write (90,6300) itmin,itermax,itermix,icrash0,numeric,icorr
      write (90,6400) (eps(l),l = 1,5),epsa
      write (90,6500) phi,alpha,iacc,alphmin,alphmax,sigma,vlconv

      write (90,8000)
      write (90,8100) diffusconv
      write (90,8200) semiconv
      write (90,8300) difftacho
      if (microdiffus) then
         if (lmicro.eq.2) then
            cdiff = 'Chapman & Cowling'
         else
            cdiff = 'Paquette'
         endif
         write (90,8400) microdiffus,cdiff
      else
         write (90,8450) microdiffus
      endif
      write (90,8500) rotation
      if (diffover) then
         if (lover.eq.22) cover = 'below convective envelope'
         if (lover.eq.23) cover = 'below convective envelope'
         if (lover.eq.24) cover = 'above core'
         if (lover.eq.25) cover = 'below pulse'
         if (lover.eq.26) cover = 'above pulse'
         if (lover.eq.27) cover = 'above core and below envelope'
         if (lover.eq.28) cover = 'below and above pulse'
         if (lover.eq.29) cover = 'core+pulse+envelope (up+down)'
         if (lover.eq.30) cover = 'overshoot below all conv. zones'
         if (lover.eq.31) cover = 'overshoot above all conv. zones'
         if (lover.eq.32) cover = 'overshoot everywhere'
         if (lover.eq.22) then
            write (90,8550)
         else
            write (90,8600) diffover,etaturb,cover
         endif
      else
         write (90,8650) diffover
      endif
      write (90,8700) thermohaline
      if (idup.eq.1.or.idup.eq.3) write (90,8800)
      write (90,'(/)')

*____________________________________
***   set values for fixed parameters
*------------------------------------

      modeli = 0
      model_old = 0

***   conserv the initial value of idiffvr (needed in case of disc-locking)

      idiffvr0 = idiffvr

***   equation of state

      Ztotioni = .true.
c      Ztotioni = .false.          ! determine if elements are totaly ionised or neutral (partial ionization is indep)
      addfitZ = .false.
      numderiv = .false.
      addprio = .true.
c      addprio = .false.

***   atmosphere

      if (ntprof.eq.2) tau0 = pw23
      taueff = max(pw23,tau0)
      if (ntprof.eq.1) then
         taueff = 1.d0
         taulim = 0.5d0
      endif
      teffagb = 5500.d0

***   convective zone definition

      delmcz = 5.d-5
      dup3 = .false.
      if (dtnmix.eq.'T'.or.dtnmix.eq.'u'.or.dtnmix.eq.'U') dtnmix = 't'
      if (dtnmix.eq.'F'.or.dtnmix.eq.'g'.or.dtnmix.eq.'G') dtnmix = 'f'
c      nczm = max(1,nczm)

***   nucleosynthesis

      tolnuc0 = tolnuc
      tnucmin = 5.d5
      shlim = 1.d1

***   accretion

      rir = 0.25d0
      accrw = 0.5d0
      if (accphase.ne.8) menv = menv*msun

***   rotation (define angular grid)

      dtheta = 0.5d0*pi/dble(intmax-1)
      dtheta2 = 0.5d0*dtheta
      do i = 1,intmax
         theta = dble(i-1)*dtheta
         cost(i) = cos(theta)
         sint(i) = sin(theta)
         sint2(i) = sint(i)*sint(i)
         plegendre(i) = 1.5d0*cost(i)*cost(i)-0.5d0
      enddo

***   hydro
*     0 : hydrostatic, 1 : hydrodynamic (2,3 idem but rotation included)

      dynfac = 0.d0
      ihydro = 0
      hydro = .false.
      hydrorot = .false.
      if (chydro.eq.'t'.or.chydro.eq.'T'.or.chydro.eq.'1') ihydro = 1
      if (chydro.eq.'2') ihydro = 2
      if (chydro.eq.'3') ihydro = 3
      if (ihydro.eq.1.or.ihydro.eq.3) hydro = .true.
      if (hydro) dynfac = 1.d0
      if (rotation.and.hydro) hydrorot = .true.            ! Correction 17/04/2020
      if (lgrav.eq.4.and..not.hydro) then
         write (nout,*) 'if you use lgrav = 4, hydro MUST be on !'
         stop 'rmodpar : bad combination of lgrav-hydro'
      endif
      hydro0 = hydro
      if (hydro.and.ivisc.gt.0) write (90,6700) nq0,nretry

***   opacity

      tlhkap = 8.d3

***   mesh law

      dmdup = 2.5d-6/mtini

***   maximum allowed corrections for Newton-Raphson

      iresmax = 9
      if (icrash.gt.iresmax) icrash = iresmax
      epsmax(1) = 1.d20
      do l = 2,4
         epsmax(l) = 1.d6
      enddo
      epsmax(5) = 1.d10
      if (numeric.eq.1) write (90,7000)
      if (numeric.eq.2.or.numeric.eq.3) write (90,7100)
      if (numeric.eq.2.or.numeric.eq.4) write (90,7200)
      write (90,'(1x)')

***   time steps

      if (fts.eq.0.d0) fts = 1.d50
      if (ftsh.eq.0.d0) ftsh = 1.d50
      if (ftshe.eq.0.d0) ftshe = 1.d50
      if (ftsc.eq.0.d0) ftsc = 1.d50
      if (fkhdt.eq.0.d0) fkhdt = 1.d50
      if (ftacc.eq.0.d0) ftacc = 1.d50
      if (ftst.eq.0.d0) ftst = 1.d50
      if (ftnuc.eq.0.d0) ftnuc = 1.d50
      if (dlnenuc0.eq.0.d0) dlnenuc = 1.d50
      facdt0 = facdt
      fkhdt0 = fkhdt
      fts0 = fts

***   open specific ascii and binary files in case of
***   angular momentum transport

      if (irotbin.eq.1.or.irotbin.eq.2) then
C     Modifs CC ondes (2/04/07)
         open (unit = 93,file = 'modang.bin',form = 'unformatted',
     &        status = 'unknown')
c    &        status = 'unknown',convert = 'big_endian')
         open (unit = 94,file = 'nextang.bin',form = 'unformatted',
     &        status = 'unknown')
c    &        status = 'unknown',convert = 'big_endian')
      endif

* Read the Paquette coefficients for radiative diffusion
      if (microdiffus.or.thermohaline.or.igw) then
***   open specific ascii file in case of microscopic diffusion
         open (unit = 85,file = 'evolvarang',status = 'unknown',form
     $        ='unformatted')
         open (unit = 71,file='paq_rp1.dat',status='old',action='read')
         open (unit = 72,file='paq_rp2.dat',status='old',action='read')
         open (unit = 73,file='paq_rp3.dat',status='old',action='read')
         open (unit = 74,file='paq_rp4.dat',status='old',action='read')
         open (unit = 75,file='paq_ap1.dat',status='old',action='read')
         open (unit = 76,file='paq_ap2.dat',status='old',action='read')
         open (unit = 77,file='paq_ap3.dat',status='old',action='read')
         open (unit = 78,file='paq_ap4.dat',status='old',action='read')
         do k = 1,50
            read (71,*) ldum,(paq_c(1,j,k),j = 1,4)
            read (75,*) ldum,(paq_c(1,j,k),j = 5,8)
            read (72,*) ldum,(paq_c(2,j,k),j = 1,4)
            read (76,*) ldum,(paq_c(2,j,k),j = 5,8)
            read (73,*) ldum,(paq_c(3,j,k),j = 1,4)
            read (77,*) ldum,(paq_c(3,j,k),j = 5,8)
         enddo

         do k = 1,50
            read (74,*) ldum,(paq_d(j,k),j = 1,4)
            read (78,*) ldum,(paq_d(j,k),j = 5,8)
         enddo
         close (71)
         close (72)
         close (73)
         close (74)
         close (75)
         close (76)
         close (77)
         close (78)
      endif

************************************************************************

 900  format (/,43x,f4.2)
 1000 format (/,25x,i4,11x,i3,10x,f6.2,10x,f8.6)
 1100 format (25x,l1,10x,l1,12x,1pd8.2,13x,1pd8.2)
 1200 format (27x,l1,11x,l1,11x,l1,11x,l1,12x,l1)
 1300 format (25x,i1,10x,l1)
 1400 format (26x,f6.4,11x,a1,11x,l1,10x,i1,8x,i3)
 1500 format (24x,l1,10x,i1,12x,d9.3,12x,f6.4)
 1600 format (26x,i1,8x,f4.2,9x,f4.2,9x,i1)
 1700 format (24x,f7.5,11x,i1,11x,d8.2)
 1800 format (22x,i3,11x,f5.2,11x,d7.1)
 1900 format (26x,i1,12x,a1,11x,d7.1,10x,1pd8.2,12x,i3)
 2000 format (27x,l1,13x,i1,12x,i2,11x,l1,11x,d9.3)
 2100 format (27x,l1,14x,d10.4,11x,d8.2,14x,l1)
 2110 format (29x,a1,14x,a8,11x,f5.2)
 2120 format (23x,a4,11x,f4.1,9x,d8.2,13x,d8.2,17x,l1)
 2130 format (26x,d8.2,12x,i1,12x,i1,12x,d9.3,14x,d8.2)
 2200 format (25x,i1,11x,i2,11x,i2,11x,d8.2,11x,i2)
 2300 format (33x,i2,2x,i2,2x,i2,2x,i2)
 2400 format (25x,i1,13x,i1,10x,i1,13x,d8.2,12x,f5.2,9x,d10.4)
 2500 format (23x,d8.2,9x,f5.3,11x,f5.3,11x,f4.2,13x,d8.2)
 2600 format (25x,a1,10x,i1,7x,f6.3,9x,d10.4)
 2700 format (25x,i4,13x,i3,11x,d9.3)
 2800 format (26x,d8.2,11x,d8.2,12x,d8.2,10x,d8.2,10x,d8.2)
 2900 format (24x,d9.3,10x,d9.3,10x,d9.3,10x,f4.2,10x,d8.2)
 3000 format (27x,a1,8x,f5.3,9x,f5.2,10x,f4.2,9x,f4.2,10x,f5.3,9x,f5.3)
 3100 format (25x,i2,12x,i3,12x,i3,11x,i2,12x,i1,10x,a1)
 3200 format (21x,d8.2,6x,d8.2,8x,d8.2,8x,d8.2,6x,d8.2,8x,d8.2)
 3300 format (23x,f4.2,10x,f4.2,9x,l1,12x,f4.2,12x,f4.2,10x,f5.3,11x,
     &     l1,///)
 4100 format (/,'Stellar Evolution Code - computation listing - ',a15,
     &     ' - on host : ',a20,//,' AGB network : version ',0pf4.2,//,
     &     ' PARAMETER CARD for STAREVOL : version AGB ',0pf4.2,/)
 4200 format (1x,'evol. sequence: maxmod =',i4,', imodpr = ',i3,
     &     ', mtini = ',0pf6.2,', zkint = ',0pf8.6)
 4300 format (1x,'eos           : addH2 = ',l1,', addHm = ',l1,
     &     ', tmaxioH = ',1pd8.2,', tmaxioHe = ',1pd8.2)
 4400 format (17x,'ionizHe = ',l1,', ionizC = ',l1,', ionizN = ',l1,
     &     ', ionizO = ',l1,', ionizNe = ',l1)
 4500 format (1x,'gravitation   : lgrav = ',i1,', lnucl = ',l1)
 4600 format (1x,'convection    : alpha_mlt_hp = ',0pf6.4,', dtnmix = ',
     &     a1,', mixopt = ',l1,', hpmlt = ',i1,', nczm =',i3)
 4700 format (17x,'ihro = ',l1,', iturb = ',i1,', etaturb = ',1pd9.3,
     &     ', alphatu = ',0pf6.4)
 4800 format (17x,'novopt = ',i1,', aup = ',f4.2,', adwn = ',f4.2,
     &     ', idup = ',i1)
 4900 format (1x,'atmosphere    : tau0 = ',0pf7.5,', ntprof = ',i1,
     &     ', taulim = ',1pd8.2)
 5000 format (1x,'mass loss     : mlp =',i3,', etapar = ',0pf5.2,
     &     ', dmlinc = ',1pd7.1)
 5100 format (1x,'nuclear       : nucreg = ',i1,', nuclopt = ',a1,
     &     ', tolnuc = ',1pd7.1,', ftnuc = ',1pd8.2,', znetmax = ',i3)
 5200 format (1x,'diffusion     : idiffcc = ',l1,', idiffnuc = ',i1,
     &     ', idiffty = ',i2,', diffzc = ',l1,', diffst = ',1pd9.3)
 5300 format (17x,'zgradmu = ',l1,', grad_crit = ',1pd10.4,
     &     ', del_zc = ',1pd8.2,', semiconve = ',l1)
 5310 format (1x,'rotation      : omegaconv = ',a1,', Dh_prescr = ',a8,
     &     ', dm_ext = ',0pf5.2)
 5320 format (17x,'Dv =  ',a4,', om_sat = ',0pf4.1,', D_zc = ',
     &     1pd8.2,', disctime = ',1pd8.2,', thermal_equilibrium = ',l1)
 5330 format (17x,'diffvr = ',1pd8.2,', idiffvr = ',i1,', idiffex = ',
     &     i1,', diffcst = ',1pd9.3,', breaktime = ',1pd8.2)
 5400 format (1x,'time-dep conv.: nmixd = ',i1,', itmind = ',i2,
     &     ', itmaxd = ',i2,', rprecd = ',1pd8.2,', nretry = ',i2)
 5500 format (1x,'extra mixing  : Ov Thl Tac Mic  ',i2.2,', ',i2.2,', ',
     &     i2.2,', ',i2.2)
 5600 format (1x,'accretion     : iaccr = ',i1,', accphase = ',i1,
     &     ', itacc = ',i1,', massrate = ',1pd8.2,', massend = ',0pf5.2,
     &     ', menv = ',1pd10.4)
 5700 format (17x,'ric = ',1pd8.2,', prrc = ',0pf5.3,', xiaccr = ',
     &     0pf5.3,', fdtacc = ',0pf4.2,', alphaccr = ',1pd8.2)
 5800 format (1x,'hydrodynamics : hydro = ',a1,', ivisc = ',i1,
     &     ', q0 = ',f6.3,', mcut = ',1pd10.4)
 5900 format (1x,'shell masses  : maxsh = ',i4,', nresconv = ',i3)
 6000 format (17x,'dvarma = ',1pd8.2,', dvarmi = ',1pd8.2,
     &     ', denucma = ',1pd8.2,', dmrma = ',1pd8.2,', dmrmi = ',
     &     1pd8.2)
 6100 format (1x,'time-step     : dtin = ',1pd9.3,', dtmin = ',1pd9.3,
     &     ', dtmax = ',1pd9.3,', facdt = ',0pf4.2,', fkhdt = ',1pd8.2)
 6200 format (17x,'ishtest = ',a1,', fts = ',0pf5.3,', ftsh = ',0pf5.2,
     &     ', ftshe = ',0pf4.2,', ftsc = ',0pf4.2,', ftacc = ',0pf5.3,
     &     ', ftst = ',0pf5.3)
 6300 format (1x,'iterations    : itmin = ',i2,', itermax = ',i3,
     &     ', itermix = ',i3,', icrash = ',i2,', numeric = ',i1,
     &     ', icorr = ',a1)
 6400 format (1x,'tolerances    : u = ',1pd8.2,', r = ',1pd8.2,
     &     ', lnf = ',1pd8.2,', lnT = ',1pd8.2,', l = ',1pd8.2,
     &     ', ang = ',1pd8.2)
 6500 format (1x,'convergence   : phi = ',0pf4.2,', alpha = ',0pf4.2,
     &     ', iacc = ',l1,', alphmin = ',0pf4.2,', alphmax = ',
     &     0pf4.2,', sigma = ',0pf5.3,', vlconv = ',l1,/)
 6700 format (' ** viscous pressure allowed to spread over nq0 =',i2,
     &     ' shells on each side of the shock fronts',/,' ** mesh ',
     &     'applied nretry =',i2,' shells on each side of the shock ',
     &     'fronts',/)
 7000 format (2x,'WARNING : convergence not followed in the surface ',
     &     'layers (T < tmaxioHe)')
 7100 format (2x,'WARNING : use first order accuracy in spatial ',
     &     'derivatives')
 7200 format (2x,'WARNING : use arithmetic mean for the opacity')
 7300 format (2x,'WARNING : Ledoux criterion used without ',
     &     'semiconvective mixing')
c 7400 format (2x,'WARNING : diffusive processes in radiative zones but',
c     &     ' convective zones treated as homogeneous (change diffzc ?)')
 8000 format (1x,'MIXING PARAMETERS')
 8100 format (3x,'diffusive mixing in convective zones : ',l1)
 8200 format (3x,'semiconvection         : ',l1)
 8300 format (3x,'tachocline mixing      : ',l1)
 8400 format (3x,'microscopic diffusion  : ',l1,', ',a18)
 8450 format (3x,'microscopic diffusion  : ',l1)
 8500 format (3x,'rotational mixing      : ',l1)
 8550 format (3x,'Internal Gravity Waves :  T')
 8600 format (3x,'overshooting           : ',l1,', fover = ',0pf6.4,
     &     ', ',a36)
 8650 format (3x,'overshooting           : ',l1)
 8700 format (3x,'thermohaline mixing    : ',l1)
 8800 format (3x,'gradiant neutrality    : T')
      
************************************************************************

      return
      end

