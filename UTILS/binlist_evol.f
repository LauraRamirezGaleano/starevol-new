
      PROGRAM BINLIST_EVOL

************************************************************************
* Create the ascii output (the listing) from binary output evol file   *
* (structure and chemical profiles)                                    *
*                                                                      *
* $LastChangedDate:: 2017-11-10 13:22:58 +0100 (Ven, 10 nov 2017)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 109                                                         $ *
*                                                                      *
************************************************************************

      implicit none

************************************************************************

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.atm'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.data'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.eos'
      include 'evolcom.grad'
      include 'evolcom.ion'
      include 'evolcom.lum'
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
      include 'evolcom.transp'
      include 'evolcom.var'
C Modif CC ondes (11/07/07)
      include 'evolcom.igw'

************************************************************************

      character name*150,nameang*20,namediff*150,atmlib*200
      character*100 dirmongo,diresults,dirbatch,opalib
      character cmodel*7
      character*1 cdum

      integer n1,n2,n3,iscr,rcrz
c      double precision rklth,tklth,opaclth

      integer itmax,ishockt,ishockb,imack
      integer imod,error,pass
      integer lenci,lenca,sm
      integer nsequence
      integer i,j,k,l,ip,j1,ij
      integer crzi
      integer klenv,klpulse,klcore
      integer inshock,outshock
C Modif CC (02/02/07) 
      integer id,nmodmdix
C Modif CC ondes (12/04/07)
      integer ndb,ndt
C Modif CC ondes (11/07/07)
      integer ii

      integer BCE,k_He,kconv

      double precision dturb,dmic,vdmic,dtinv,vvro,vvmue,deltamu
      double precision sr
      double precision ebind
      double precision version,FeH
      double precision Dhold
      double precision abmurj,abmuj
C Modif CC (02/02/07) 
      double precision abmuliss
C Modif CC (02/02/07) 
      double precision abmax,abmin
      double precision xspr,y
      double precision epot,eint,ebindenv,epotenv,eintenv
      double precision comp,heat,work,etherm,vve,vvp,flux
      double precision renucl,enucl_save
      double precision Dtot
      double precision mack,Mackmax,Mackmin,tshock
C Modif CC ondes (11/07/07)
C      double precision dtcour,lmax,tmax,enucmax
      double precision dtcour,tmax,enucmax
      double precision zbar,ztilde,zx,zy,corr
      double precision dift,dife,rpvisc,rsconv,dlumdt
      double precision vrray,egravj1,rkonvmj
      double precision xpsim,urm,auxm,omm,dholdm,vxpsim,vom,vxpsi_save
      double precision tamp1a,tamp1b,raprho,aalpha,divis,epsmoym,rapeps
      double precision fepsi,bbeta,omm2,tamp1,tamp2a,tamp2b,tamp2,tamp3
      double precision tamp,tamp4,tamp6a,dt,tamp6,rhmoy,epsmoy
      double precision StAdv,StB,StTh,StNG,StNS,xlambdam,xxnorm
      double precision StTh1,StTh2,StB1,StB2,StB3,StB4
      double precision om,ur,vxpsi,tamp5a,tamp5b,tamp5
      double precision sum,Dnu,Dnuech,errr
      double precision sum_Dp,rBCE,sum_BCE
      double precision addvar,addsum
      double precision T_tot,t_BCE,Dp,nu_max
      double precision t_He,sum_He,r_He
      double precision hpbce,tc,tg
      double precision hptce,rtce,rc,tc_hp,rc_hp,rce,tc_r,rc_r
      double precision tc_m,rc_m,tc_max,rc_max,secd,mce,factdc
      double precision lgdmic,lgdturb,vpsi
      double precision Ri_n,S_n,Y_n !Added LR Jan 2025
      logical shockdetect

      integer kbce,ktce
c      double precision xspsol,zsol
c      character*5 refsolar

      double precision system
c      double precision sigc
      double precision tautime,Ereac,taureac
      character*120 cmd
      integer dummy

c      parameter(mlth = 19,nlth = 85,itlth = 120)
      common /difcirc/ dift(nsh),dife(nsh)
      common /coefdiff/ dturb(nsh),dmic(nsh),vdmic(nsh)
      common /calcDh/ Dhold(nsh)
      common /gradmuj/ abmurj(nsh),abmuj(nsh)
c      common /lthopac/ rklth(19),tklth(85),opaclth(19,85,120)
c      common /lthopac/ rklth(19),tklth(63),opaclth(19,63,104)
      common /overshoot/ klcore,klenv,klpulse
C Modif CC ondes (11/07/07)
C      common /hydrodyn/ lmax,tmax,Mackmax,enucmax,itmax,ishockb,ishockt
      common /hydrodyn/ tmax,Mackmax,enucmax,itmax,ishockb,ishockt
      common /metal/ FeH
      common /moy/ rhmoy(nsh),epsmoy(nsh)
      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)
      common /nucleaire/ tautime(nsh,nreac),Ereac(nsh,nreac),taureac(nsh
     $     ,nreac)
!!      common /Dvprescription / Dv_prescrp,thetac(nsh)

!!      character*11 Dv_prescrp
!!      double precision thetac

      integer mlth,nlth,itlth
      double precision rklth,tklth,opaclth
      double precision rklth_as09,tklth_as09,opaclth_as09
      double precision rklth_gn93,tklth_gn93,opaclth_gn93
      double precision rklth_grid,tklth_grid,opaclth_grid
      double precision rklth_ay18,tklth_ay18,opaclth_ay18  ! TD 12/2019

c      common /lthopac/ opaclth(19,85,155),tklth(85),rklth(19),
c     &     mlth,nlth,itlth
      common /lthopac/ opaclth(19,81,155),tklth(81),rklth(19),   ! TD 03/2020
     &     mlth,nlth,itlth
      common /lthopac_as09/ rklth_as09(19),tklth_as09(85),
     &     opaclth_as09(19,85,155)
      common /lthopac_gn93/ rklth_gn93(19),tklth_gn93(85),
     &     opaclth_gn93(19,85,120)
      common /lthopac_grid/ rklth_grid(19),tklth_grid(63),
     &     opaclth_grid(19,63,104)
c      common /lthopac_ay18/ rklth_ay18(19),tklth_ay18(85),  ! TD 12/2019
c     &     opaclth_ay18(19,85,155)
      common /lthopac_ay18/ rklth_ay18(19),tklth_ay18(81),  ! TD 03/2020
     &     opaclth_ay18(19,81,155)

c Ajout par TD (22/02/2018) + modif (25/04/2019)
      common /coefdiffpaqhe/ xvdiffhe(nsh),xvdiffhenoc(nsh)
      double precision xvdiffhe,xvdiffhenoc
      common /coefdiffpaqli/ xvdiffli(nsh),xvdli(nsh),
     &     xvthermli(nsh),xdmicli(nsh)
      double precision xvdiffli,xvdli,xvthermli,xdmicli
      common /coefdiffpaqc12/ xvdiffc12(nsh),xvdiffc12noc(nsh)
     $     ,xvdiffc13(nsh),xvdiffc13noc(nsh)
      double precision xvdiffc12,xvdiffc12noc
     $     ,xvdiffc13,xvdiffc13noc
      common /coefdiffpaqn14/ xvdiffn14(nsh),xvdiffn14noc(nsh)
     $     ,xvdiffn15(nsh),xvdiffn15noc(nsh)
      double precision xvdiffn14,xvdiffn14noc
     $     ,xvdiffn15,xvdiffn15noc      
      common /coefdiffpaqo16/ xvdiffo16(nsh),xvdiffo16noc(nsh)
     $     ,xvdiffo17(nsh),xvdiffo17noc(nsh),xvdiffo18(nsh)
     $     ,xvdiffo18noc(nsh)
      double precision xvdiffo16,xvdiffo16noc
     $     ,xvdiffo17,xvdiffo17noc,xvdiffo18,xvdiffo18noc
      common /coefdiffpaqne20/ xvdiffne20(nsh),xvdiffne20noc(nsh)
      double precision xvdiffne20,xvdiffne20noc
      common /coefdiffpaqna23/ xvdiffna23(nsh),xvdiffna23noc(nsh)
      double precision xvdiffna23,xvdiffna23noc
c Fin ajout
   
      integer bip1
      integer filesize
      double precision rtau,rT,rhotab,Fconvtab
      double precision Ttabeff,gtabeff,Ztabeff
      double precision T_Zeff(6),rho_Zeff(6),Fconv_Zeff(6)
      double precision tableTZ(128,59,10),tablerhoZ(128,59,10),
     &     tableFcZ(128,59,10)
      double precision sTZ,srhotabZ,sFconvtabZ
      double precision bip,taumin,taumax
      character(len=89) ref,nom_fichier
      character(len=len(ref)) DATUM
      character(len=64) fmt_init
      character(len=64) fmt_lect1,fmt_lect2,fmt_lect3
      character(len=5) charT
      character(len=4) charg,charZ
      character(len=21) nameatm,nameg,nameZ,nameinit
      character atmtype*7
      common /atmospheres/ rtau(127,10,59,6),rT(127,10,59,6),
     &     rhotab(127,10,59,6),Fconvtab(127,10,59,6),Ttabeff(59),
     &     gtabeff(10),Ztabeff(6),taumin,taumax,filesize


c      common /solar/ xspsol(nis+6,2),zsol,refsolar

      dimension rpvisc(nsh),rsconv(nsh),vvro(nsh),vvp(nsh),vvmue(nsh)
      dimension Dtot(nsh)
      dimension mack(nsh),sr(nsh)

      dimension crzi(nsh),rcrz(nsh),vve(nsh)
      dimension xspr(nsh,nsp),y(nsp)
      dimension comp(nsh),heat(nsh),work(nsh),etherm(nsh),dlumdt(nsh),
     &     renucl(nsh),corr(nsh),enucl_save(nsh)
      dimension ur(nsh),vxpsi_save(nsh)
      dimension StAdv(nsh),StB(nsh),StTh(nsh),StNG(nsh),StNS(nsh)
      dimension StTh1(nsh),StTh2(nsh),StB1(nsh),StB2(nsh),StB3(nsh),
     &     StB4(nsh)
C Modif CC (02/02/07) 
      dimension abmuliss(nsh)
      dimension Ri_n(nsh), S_n(nsh) !LR 21 Jan 2025

      external system
      
C Modif Nadège (16/02/10)
      double precision Nu2,Nt2,bruntV2,bruntV
      
      dimension Nu2(nsh),Nt2(nsh),bruntV2(nsh),bruntV(nsh)
      dimension rBCE(nsh),addvar(nsh)

      dimension lgdmic(nsh),lgdturb(nsh),vpsi(nsh)

C     Modif 24/10/2023 -------------------------------------------------
                                                                       
      character hydrodynamics*1
      double precision evolpar_version
      double precision dm_ext
      integer n_turbul(1)
      double precision Tfix,om_turbul,PM_turbul
      double precision disctime,dlnvma0,dlnvmi0,dlnenuc0
      double precision massend0,massrate0
      double precision tol_u,tol_lnr,tol_lnf,tol_lnT,tol_l,tol_ang
      logical vlconv
      namelist /starevol_parameters/ evolpar_version,
     &     maxmod,imodpr,mtini,zkint,
     &     addH2,addHm,tmaxioH,tmaxioHe,
     &     ionizHe,ionizC,ionizN,ionizO,ionizNe,
     &     lgrav,lnucl,
     &     alphac,dtnmix,mixopt,hpmlt,nczm,
     &     ihro,iturb,etaturb,fover,alphatu,
     &     novopt,aup,adwn,idup,
     &     tau0,ntprof,taulim,
     &     mlp,etapar,dmlinc,
     &     clumpfac,zscaling,
     &     nucreg,nuclopt,tolnuc,ftnuc,znetmax,
     &     idiffcc,idiffnuc,idiffty,diffzc,diffst,
     &     zgradmu,nu_add,del_zc,ledouxmlt,
     &     Tfix,om_turbul,n_turbul,PM_turbul,
     &     omegaconv,Dh_prescr,dm_ext,
     &     Dv_prescr,om_sat,D_zc,disctime,thermal_equilibrium,
     &     diffvr,idiffvr,idiffex,diffcst,breaktime,
     &     nmixd,itmind,itmaxd,rprecd,nretry,
     &     lover,lthal,ltach,lmicro,
     &     iaccr,accphase,itacc,massrate,massend,menv,
     &     ric,prrc,xiaccr,fdtacc,alphaccr,
     &     hydrodynamics,ivisc,q0,mcut,
     &     maxsh,nresconv,
     &     dlnvma,dlnvmi,dlnenuc,dmrma,dmrmi,
     &     dtin,dtmin,dtmax,facdt,fkhdt,
     &     ishtest,fts,ftsh,ftshe,ftsc,ftacc,ftst,
     &     itermin,itermax,itermix,icrash,numeric,icorr,
     &     tol_u,tol_lnr,tol_lnf,tol_lnT,tol_l,tol_ang,
     &     phi,alpha,iacc,alphmin,alphmax,sigma,vlconv

************************************************************************
      
      pi = 3.1415926535d0
      sec = 3.1557807d7
      seci = 1.d0/sec
      c = 2.99792458d10
      g = 6.67259d-8
      boltz = 1.380658d-16
      h = 6.6260755d-27
      sig = 2.d0*pi**5*boltz**4/(15.d0*h**3*c*c)
      sigc = 4.d0*sig/c
      avn = 6.0221367d23
      rk = avn*boltz
      mprot = 1.6726231d-24
      ech = 1.60217733d-19
      econv = ech*avn*1.d13
      msun = 1.9891d33
      rsun = 6.9599d10
      lsun = 3.846d33
      cmodel = ' '
      code_version = 9.15d0

      pim2 = 2.d0*pi
      pim4 = 4.d0*pi
      pim8 = 8.d0*pi
      pw14 = 1.d0/4.d0
      pw13 = 1.d0/3.d0
      pw23 = 2.d0/3.d0
      pw34 = 3.d0/4.d0
      pw43 = 4.d0/3.d0
      pw53 = 5.d0/3.d0
      vlog10 = log(10.d0)
      error = 0
      saha = (sqrt(pim2*boltz*9.1093897d-28)/h)**3 !2.4147031864d15
      

***   define reference solar composition 
      call getenv('SOLCOMP',refsolar)
      refsolar = trim(refsolar)
      print *,'refsolar',refsolar
      if (refsolar.eq.'AGS05') then
         do i = 1,nis+6
            xspsol(i) = xspref(2,i)
         enddo
** Modif Asplund 2009 - AP feb 2012 -->
      elseif (refsolar.eq.'AGSS09') then
         do i = 1,nis+6
            xspsol(i) = xspref(3,i)
         enddo
**  Modif Asplund 2009 - AP feb 2012 <--
      elseif (refsolar.eq.'GRID') then
         do i = 1,nis+6
            xspsol(i) = xspref(4,i)
         enddo
**  Modif Young 2018 - TD Dec 2019 <--
      elseif (refsolar.eq.'AY18') then
         do i = 1,nis+6
            xspsol(i) = xspref(5,i)
         enddo         
      elseif (refsolar.eq.'GN93') then
         refsolar = 'GN93'
         do i = 1,nis+6
            xspsol(i) = xspref(1,i)
         enddo
      else
         stop 'Wrong solar composition, check $SOLCOMP variable'
      endif

*** Define model atmospheres
      call getenv ('STAR_ATM',atmtype)
      atmtype = trim(atmtype)
      if (atmtype.eq.'none') write(6,*) ' No model atmospheres used'
        print *,atmtype

*------------------
***  initialisation
*------------------

      do i = 1,nsh
         Dsc(i) = 0.d0          ! semiconvection
         Dthc(i) = 0.d0         ! thermohaline
         Dherw(i) = 0.d0        ! overshoot
         Dhd(i) = 0.d0          ! rotation
c         coefDtacho(i) = 0.d0   ! tachocline - commenter (TD 09/2019)
         dmic(i) = 0.d0         ! microscopic diffusion
         vdmic(i) = 0.d0        ! microscopic velocity
c         Dturbul(i) = 0.d0      ! ad oc turbulence (TD 09/2019)
      enddo
      
c..   imod : 9 : nextini or modini, else filename_xxxx[abd]
c..   sm = 0 : generate kippenhahn file
c..   sm = 1 : generate smfile : filename.p[i]
c..   ex : 1 1 shock2007_0200b

      nsequence = 9999
      read (5,999) imod,sm,name
      write (*,5) trim(name),imod
 5    format (/,' filename : ',A,', imod = ',i6)
      if (imod.ne.9) then
         i = lenci(name)
         cmodel = name(i+1:i+6)
c         call make_number (cmodel,nsequence)
         read (cmodel,'(i7)') nsequence
         print *,cmodel,nsequence
      endif

      call getenv ('DIR_SMONGO',dirmongo)
      call getenv ('DIR_RESULTS',diresults)
      call getenv ('DIR_BATCH',dirbatch)

***   standard output unit
      nout = 6

      if (imod.eq.9) then
         open (unit = 91,file = trim(diresults) // trim(name),
     &        form = 'unformatted',status = 'old', action = 'read')
         open (unit = 99,file = trim(diresults) // '/starevol.par',
     &        status = 'old', action = 'read')
         
         read(99, NML=starevol_parameters) 
         version = evolpar_version
         rewind(99)
      else
         open (unit = 91,file = trim(diresults) // trim(name) // 'b',
     &        form = 'unformatted', status = 'old', action = 'read')
         open (unit = 99,file = trim(diresults) // trim(name) // 'd',
     &        status = 'old', action = 'read')
c         if (idiffcc.and.nphase.ge.2) then
         if (idiffcc.and.nphase.ge.1) then   ! Modif by TD Nov.2018
         open (unit = 98,file = trim(diresults) // trim(name) // 'o',         ! Add by TD Fev.2018
     &        form = 'unformatted',status = 'old', action = 'read')
         endif
            
         cdum = ' '
         do while (cdum.ne.'S')
            read (99,'(a1)') cdum
         enddo
         read (99,20) version
 20      format (///,43x,0pf4.2)
         backspace(99)
         backspace(99)
      endif
      write (*,30) version
 30   format (' STAREVOL version : ',f4.2)
c      open (unit = 92,file = '/home/decressin/tmp/nextini.bin',
c     &     form = 'unformatted',
c     &     status = 'unknown')
     

*********************************************************
*  Set the reference solar chemical composition adopted *
*********************************************************

c      call getenv ('SOLCOMP',refsolar)
c      refsolar = trim(refsolar)
c      if (refsolar.eq.'GN93') refsolchem = 1
c      if (refsolar.eq.'AGS05') refsolchem = 2
      
c	print *,refsolchem,refsolar

*-----------------------
***  read opacity tables
*-----------------------
c#ifdef GRID
c      call getenv ('DIR_OPAGRID',opalib)
c#elif AS09
c      call getenv ('DIR_OPAAS09',opalib)
c#else
c      call getenv ('DIR_OPA',opalib)
c#endif
c      call set_opal_dir(opalib)

      if (refsolar.eq.'AGS05') then
         call getenv('DIR_OPA2',opalib)
         stop 'No opacity files: check with Ana'
 
**  Modif Asplund 2009 - AP feb 2012 -->      
      elseif (refsolar.eq.'AGSS09') then
         call getenv('DIR_OPAAGSS09',opalib) !LIB_AGSS09

         mlth = 19
         nlth = 85
         itlth = 155
         do i = 1,mlth
            rklth(i) = rklth_as09(i)
         enddo
         do i = 1,nlth
            tklth(i) = tklth_as09(i)
         enddo
         do i = 1,mlth
            do j = 1,nlth
               do k = 1,itlth
                  opaclth(i,j,k) = opaclth_as09(i,j,k)
               enddo
            enddo
         enddo

**  Modif Asplund 2009 - TD Dec 2019 -->      
      elseif (refsolar.eq.'AY18') then
         call getenv('DIR_OPAAY18',opalib) !LIB_AY18

         mlth = 19
c         nlth = 85
         nlth = 81      ! Modif pour new opa Ferguson
         itlth = 155
         do i = 1,mlth
            rklth(i) = rklth_ay18(i)
         enddo
         do i = 1,nlth
            tklth(i) = tklth_ay18(i)
         enddo
         do i = 1,mlth
            do j = 1,nlth
               do k = 1,itlth
                  opaclth(i,j,k) = opaclth_ay18(i,j,k)
               enddo
            enddo
         enddo         
       
      elseif (refsolar.eq.'GRID') then
         call getenv('DIR_OPAGRID',opalib) !LIB_GRID

         mlth = 19
         nlth = 63
         itlth = 104
         do i = 1,mlth
            rklth(i) = rklth_grid(i)
         enddo
         do i = 1,nlth
            tklth(i) = tklth_grid(i)
         enddo
         do i = 1,mlth
            do j = 1,nlth
               do k = 1,itlth
                  opaclth(i,j,k) = opaclth_grid(i,j,k)
               enddo
            enddo
         enddo

      else
c..   Opacities for the Grevesse & Noels 1993 solar chemical composition
         call getenv('DIR_OPA',opalib)     !LIB_EXT

         mlth = 19
         nlth = 85
         itlth = 120
         do i = 1,mlth
            rklth(i) = rklth_gn93(i)
         enddo
         do i = 1,nlth
            tklth(i) = tklth_gn93(i)
         enddo
         do i = 1,mlth
            do j = 1,nlth
               do k = 1,itlth
                  opaclth(i,j,k) = opaclth_gn93(i,j,k)
               enddo
            enddo
         enddo         

      endif

!      print *,'opalib',opalib
      call set_opal_dir(opalib)

      
************************************************************************
*     Atmospheric table                                                *
************************************************************************

      nameinit = "Ms02500_g+2.00_z+0.00"


      fmt_lect1 = "(2x,1pe11.5,2x,1pe11.5,2x,1pe11.5,1x,1pe12.5)"
      fmt_lect2 = "(1x,1pe11.5,1x,1pe11.5,1x,1pe11.5,1x,1pe12.5)"
      fmt_init = "(2(/))"
      nameatm = nameinit


*** PHOENIX ATMOSPHERES
      if (atmtype.eq.'PHOENIX') then
         nameinit = "Ms02500_g+2.00_z-0.00"
         nameatm = nameinit
         call getenv ('DIR_ATMP',atmlib)
         nameatm(1:1) = "P"
         nameg = nameatm
         nb_geff = 10
!         nb_Teff = 76
!         nb_Teff = 60
         nb_Teff = 59
         nb_Zeff = 6
         filesize = 127
         
         gtabeff(1:nb_geff) = [(i,i=1,nb_geff)]
         gtabeff=(gtabeff)/2.
         print *,gtabeff
         Ttabeff(1:45) = [(i,i=2600,7000,100)]
!         Ttabeff(1:25) = [(i,i=2600,5000,100)]
         Ttabeff(46:nb_Teff) = [(i,i=7200,9800,200)]
!         Ttabeff(46:60) = [(i,i=7200,12000,200)]
!         Ttabeff(71:nb_Teff) = [(i,i=12500,15000,500)]

!         Ztabeff(1:nb_Zeff) = [-4.d0,-3.5d0,-2.d0,-1.5d0,-1.d0,0.d0
!     &        ,0.3d0]

         Ztabeff(1:nb_Zeff) = [-2.d0,-1.5d0,-1.d0,0.d0
     &        ,0.3d0,0.5d0]
         do ij = 1,nb_Zeff
            nameatm = nameg
            write(charZ,'(f3.1)') abs(Ztabeff(ij))            
            if (ij==nb_Zeff.or.(ij==(nb_Zeff-1))) nameatm(17:17)='+'
            nameatm(18:20)=charZ
            nameZ = nameatm
            do j = 1,nb_geff
               nameatm = nameZ
               write(charg,'(f3.1)') gtabeff(j)
               print *,charg,gtabeff(j)
               nameatm(11:13)=charg
               do i = 1,nb_Teff
                  if ((Ttabeff(i)==10000.and.gtabeff(j)==2.0
     &                 .and.Ztabeff(ij)==0.0).or.
     &                 (Ttabeff(i)==10000.and.gtabeff(j)==2.0
     &                 .and.Ztabeff(ij)==-2.0).or.
     &                 (Ttabeff(i)==10000.and.gtabeff(j)==4.0
     &                 .and.Ztabeff(ij)==-2.0).or.
     &                 (Ttabeff(i)==3000.and.gtabeff(j)==3.5
     &                 .and.Ztabeff(ij)==-3.5).or.
     &                 (Ttabeff(i)==4000.and.gtabeff(j)==4.5
     &                 .and.Ztabeff(ij)==-3.5).or.
     &                 (Ttabeff(i)==3900.and.gtabeff(j)==4.5
     &                 .and.Ztabeff(ij)==-3.5).or.
     &                 (Ttabeff(i)==10000.and.gtabeff(j)==2.0
     &                 .and.Ztabeff(ij)==-4.0).or.
     &                 (Ttabeff(i)==10000.and.gtabeff(j)==4.0
     &                 .and.Ztabeff(ij)==-4.0))  then 
                     rtau(:,j,i,ij) = 0.d0
                     rT(:,j,i,ij) = 0.d0
                     rhotab(:,j,i,ij) = 0.d0
                     Fconvtab(:,j,i,ij) = 0.d0
                     stop 'banane mauvaise table'
                  else
                     if (Ttabeff(i).lt.10000) then 
                        write(charT,'(i4)') int(Ttabeff(i))
                        nameatm(4:7)=charT
                     else 
                        write(charT,'(i5)') int(Ttabeff(i))
                        nameatm(3:7)=charT
                     endif
                     ref = trim(atmlib)//trim(nameatm)
                     open(unit = 701,file=ref,
     &                    form="formatted", action="read")
                     read (unit = 701,fmt=fmt_init)
                     do k=1,filesize
                        if (ij == nb_Zeff-1) then 
                           read (unit = 701,fmt=fmt_lect1) 
     &                          rtau(k,j,i,ij),rT(k,j,i,ij),
     &                          rhotab(k,j,i,ij),Fconvtab(k,j,i,ij)
                        else
                           read (unit = 701,fmt=fmt_lect2) 
     &                          rtau(k,j,i,ij),rT(k,j,i,ij),
     &                          rhotab(k,j,i,ij),Fconvtab(k,j,i,ij)
                        endif
                     enddo
                     close(unit = 701)
                  endif
               enddo
            enddo             
         enddo

*** MARCS ATMOSPHERES

      else if (atmtype.eq.'MARCS') then
         call getenv ('DIR_ATMM',atmlib)
         nameatm(1:1) = "M"
         nameg = nameatm

         nb_geff = 8
!         nb_Teff = 32
         nb_Teff = 23
         filesize = 55
        
         gtabeff(1:nb_geff) = [(i,i=20,55,5)] ! 8 values of log g
         gtabeff=gtabeff*0.1d0
         do i = 1,nb_Teff ! 32 values of Teff
            if (i.le.16) then
               Ttabeff(i) = 2500+(i-1)*100
            else
               Ttabeff(i) = Ttabeff(i-1)+250
            endif
         enddo            
  

       do j = 1,nb_geff
            nameatm = nameg
            if (gtabeff(j).gt.3.5) nameatm(2:2) = "p"
            do i = 1,nb_Teff
               if ((Ttabeff(i).ge.6000.and.gtabeff(j).lt.3.00).or.
     &              (Ttabeff(i).ge.7000.and.gtabeff(j).lt.4.00).or.
     &              (Ttabeff(i).lt.3000.and.gtabeff(j).ge.5.50)) exit
               write(charg,'(f3.1)') gtabeff(j)
               if (Ttabeff(i).lt.10000) then 
                  write(charT,'(i4)') int(Ttabeff(i))
                  nameatm(4:7)=charT
               else 
                  write(charT,'(i5)') int(Ttabeff(i))
                  nameatm(3:7)=charT
               endif
               nameatm(11:13)=charg
               ref = trim(atmlib)//trim(nameatm)
               open(unit = 701,file=ref,
     &              form="formatted", action="read")
               read (unit = 701,fmt=fmt_init)
               do k=1,filesize
                  read (unit = 701,fmt=fmt_lect1) rtau(k,j,i,ij)
     &                 ,rT(k,j,i,ij),rhotab(k,j,i,ij),Fconvtab(k,j,i,ij)
               enddo
               close(unit = 701)
            enddo
         enddo             
c$$$       
c$$$         do j=1,nb_geff
c$$$!            if (gtabeff(i).gt.3.5) then
c$$$!               ref(53:53) = "p"
c$$$!               ref(66:66) = "0"
c$$$!            endif
c$$$            do i=1,nb_Teff
c$$$               write(charg,'(f3.1)') gtabeff(j)
c$$$               write(charT,'(i4)') int(Ttabeff(i))
c$$$               ref(54:57)=charT
c$$$               ref(61:63)=charg
c$$$               open(unit = 701,file=ref,
c$$$     &              form="formatted", action="read")
c$$$               read (unit = 701,fmt=fmt_init)
c$$$               do k=1,filesize
c$$$                  read (unit = 701,fmt=fmt_lect) rtau(k,i,j)
c$$$     &                 ,rT(k,i,j),rhotab(k,i,j),Fconvtab(k,i,j)
c$$$               enddo
c$$$               close(unit = 701)
c$$$            enddo
c$$$         enddo             

      endif

!      do j=1,nb_geff
!         do i=1,nb_Teff
!            close(unit = unite(j,i))
!         enddo
!      enddo



*_______________________________
***   reading of the binary file
*-------------------------------
      
      call rmodpar
      if (idup.eq.1) idup = 0
      print *,'rotation        :',rotation
      print *,'overshooting    :',diffover
      print *,'semiconvection  :',semiconv
      print *,'thermohaline    :',thermohaline
      print *,'micro-diffusion :',microdiffus
      print *,'tachocline      :',difftacho

      if (ntprof.eq.4) call interpZ

      close (99)
      mixopt = .false.
      lnucl = .true.
      nuclopt = 'u'
      urca = .true.
      diffusconv = .true.
      if (rotation) then
         if (imod.eq.9) then
            if (trim(name).eq.'modini.bin') then
               nameang = 'modang.bin'
            else
               nameang = 'nextang.bin'
            endif
            print *,'open ',trim(diresults) // nameang(1:lenca(nameang))
            open (unit = 93,file = trim(diresults) // 
     &           nameang(1:lenca(nameang)),form = 'unformatted',
     &           status = 'unknown')
         else
            print *,'open ',trim(name) // 'a'
            open (unit = 93,file = trim(diresults) // trim(name) // 'a',
     &           form = 'unformatted', status = 'old', action = 'read')
         endif
      endif

      if (version.ge.2.1d0.and.version.lt.2.5d0) then
         open (unit = 99, file = trim(dirbatch) // 
     &        '/starevolnuc_evol.par', status = 'old', action = 'read')
      else if (version.ge.2.5d0) then
         open (unit = 99, file = trim(dirbatch) //
     &        '/starevolnuc_evol300.par', status = 'old', 
     &        action = 'read')
      else
         open (unit = 99, file = trim(dirbatch) // 
     &        '/starevolnuc_old_evol.par', status = 'old', 
     &        action = 'read')
      endif
      if (version.lt.2.5d0) then
         read (99,40)
 40      format (1x,//)
      end if
      call rinimod (91,93,92,94,0)
      print *,'nmod after rinimod',nmod

c.. Smoothing the soundspeed profile

      call lissage (cs,1,nmod-1,5)
      call lissage (p,1,nmod-1,5)

      if (microdiffus.or.thermohaline.or.igw) then
         if (igw) then
            open (unit=85, file = trim(diresults) // trim(name) // 'o',
     &           form = 'unformatted', status = 'old', action = 'read', 
     &           err=42)
            !print *,'LR open ',trim(diresults) // trim(name) // 'o'

            read (85, end=42, err=42) (lgdmic(k),lgdturb(k),vom(k),
     &           depottot(k),depottot_surf(k),
     &           depottot_core(k),
     &           depotondestot(k),Dondes(k),Dondeschim(k),
     &           lumwave(k),lumwave_core(k),
     &           lumondes_surf(k),lumondes_core(k),lumondestot(k),
     &           brunt_o(k),vr(k),vxpsi(k),
     &           Dmicro(k),vmicro(k),Dthc(k),phiKS(k),          
     &           deltaKS(k),tautime(k,ippg),Ereac(k,ippg),tautime(k
     &           ,ipdg), Ereac(k,ipdg),tautime(k,i2he3),Ereac(k
     &           ,i2he3),tautime(k,ihe3ag), Ereac(k,ihe3ag),tautime(k
     &           ,ibe7beta),Ereac(k,ibe7beta), tautime(k,ili7pa)
     &           ,Ereac(k,ili7pa),tautime(k,ibe7pg), Ereac(k,ibe7pg)
     &           ,tautime(k,ib8beta),Ereac(k,ib8beta), tautime(k
     &           ,ic13pg),Ereac(k,ic13pg),tautime(k,in14pg), Ereac(k
     &           ,in14pg),tautime(k,icpg),Ereac(k,icpg), k = 1,nmod)
            do k =1,nmod
               print *,'LR read bin', lumwave(k), deltaKS(k)
            enddo
         else
            open (unit=85, file = trim(diresults) // trim(name) // 'o',
     &           form = 'unformatted', status = 'old', action = 'read', 
     &           err=42)
            read (85, end=42, err=42) (Dmicro(k),vmicro(k),Dthc(k),
     $           phiKS(k),Dturbul(k),Dbar(k),Dkyle(k),coefDtacho(k), ! Ajout TD Fev.2018 + dturbul 10/2019
     $           xvdiffhe(k),xvdiffhenoc(k),xvdiffc12(k),
     $           xvdiffc12noc(k),xvdiffc13(k),xvdiffc13noc(k),
     $           xvdiffn14(k),xvdiffn14noc(k),xvdiffn15(k),
     $           xvdiffn15noc(k),xvdiffo16(k),xvdiffo16noc(k),
     $           xvdiffo17(k),xvdiffo17noc(k),xvdiffo18(k)
     $           ,xvdiffo18noc(k), xvdiffne20(k),xvdiffne20noc(k),
     $           xvdiffna23(k),xvdiffna23noc(k),
     $           deltaKS(k),tautime(k,ippg),Ereac(k ! Fin ajout TD Fev.2018
     $           ,ippg), tautime(k,ipdg), Ereac(k,ipdg),tautime(k
     $           ,i2he3), Ereac(k,i2he3),tautime(k,ihe3ag),Ereac(k
     $           ,ihe3ag), tautime(k,ibe7beta),Ereac(k,ibe7beta),
     $           tautime(k,ili7pa),Ereac(k,ili7pa),tautime(k,ibe7pg),
     $           Ereac(k,ibe7pg),tautime(k,ib8beta),Ereac(k,ib8beta),
     $           tautime(k,ic13pg),Ereac(k,ic13pg),tautime(k,in14pg),
     $           Ereac(k,in14pg),tautime(k,icpg),Ereac(k,icpg), k = 1
     $           ,nmod)
         endif
 42      close (85)
      endif

      !print *,'DEPOTTOT',depottot(:)

      write (*,45) nmod,model,time*seci
 45   format (/,' Total number of shells :',i5,', model number :',i8,/
     &     ' Age = ',1pe12.6,' yr',/)

*----------------------------------------
***   save file for color structure files
*----------------------------------------
      if (sm.eq.0) then
         open (unit=10,file=trim(dirmongo) // 'DATA/' //
     &        name(1:lenci(name)-1) // '.p99')
         rewind (99)
c     &        name(1:lenci(name)-1) // '.p99',position = 'append')
         print *,'generate file',trim(dirmongo) // 'DATA/' //
     &        name(1:lenci(name)-1) // '.p99'
         write (10,*) nmod
         r(1) = 0.d0
         do i = 1,nmod
            write (10,200) i,max(1e-10,r(i)/rsun),m(i)/msun,enucl(i),
     &           max(1.d-36,xsp(i,ih1)),max(1.d-36,xsp(i,ihe4)),
     &           max(1.d-36,xsp(i,ic12)),max(1.d-36,xsp(i,io16)),t(i),
     &           vomega(i)
         enddo
 200     format (1x,i5,2x,1pe14.8,1x,0pf14.9,1x,1pe13.6,6(1x,1pe12.6))
         close (10)
         stop
      endif

*--------------------
***   initializations
*--------------------

      do i = 1,nmod
         enucl_save(i) = enucl(i)
      enddo

      no = 0
      dtinv = 1.d0/dtn

      if (numeric.eq.2.or.numeric.eq.3) then
c.. first order accuracy in the spatial derivatives
         do i = 1,nmod
            wi(i) = 0.5d0
            wj(i) = 0.5d0
         enddo
      else
c.. second order accuracy in the spatial deriatives
         wi(1) = 0.5d0
         wj(1) = 0.5d0
         do i = 1,nmod1-1
            ip = i+1
            wi(ip) = dm(ip)/(dm(i)+dm(ip))
            wj(ip) = 1.d0-wi(ip)
         enddo
         wi(nmod) = 0.5d0
         wj(nmod) = 0.5d0
      endif
      do i = 1,nmod
         do l = 1,nsp
            xspr(i,l) = xsp(i,l)
            vxsp(i,l) = xsp(i,l)
         enddo
         fconv(i) = 0.d0
         frad(i) = 1.d0
         omi(i) = 0.d0
         psi(i) = 0.d0
         rcrz(i) = crz(i)
         renucl(i) = venucl(i)
         rpvisc(i) = vpvisc(i)
         rsconv(i) = vsconv(i)
         t(i) = exp(lnt(i))
         vt(i) = exp(vlnt(i))
         vvro(i) = ro(i)
         vve(i) = e(i)
         vvmue(i) = vmueinv(i)
         vvp(i) = p(i)
         sr(i) = s(i)
         if (vfacc(i).gt.0) iaccr = 2
         if (rotation) then
            Dhd(i) = dift(i)+dife(i)
            Dtot(i) = Dhd(i)+Dconv(i)
         endif
         comp(i) = (ro(i)-vro(i))*dtinv/ro(i) !d ln rho/dt
         heat(i) = (lnt(i)-vlnt(i))*dtinv     !d ln T/dt
         work(i) = -p(i)*(1.d0/ro(i)-1.d0/vro(i))*dtinv  !-Pd(1/rho)/dt
         etherm(i) = (ve(i)-e(i))*dtinv       !-de/dt
c                  egrav = -Pd(1/rho)/dt-de/dt)
         dlumdt(i) = (lum(i)-vlum(i))*dtinv/(lum(i)+1.d-30)
      enddo

*-------------------------------------------------------------------
***   Shock fronts, Mack number and Courant time-step determinations
*-------------------------------------------------------------------

      shockdetect = .false.
      mackmax = 0.d0
      mackmin = prrc
      dtcour = 1.d-20
      tshock = 5.d8
      ishockb = nmod1
      ishockt = 1
      imack = 1

      if (hydro.and.ivisc.gt.0.and.m(nmod1).gt.mcut.and.tmax.gt.tshock)
     &     then

c..   determine largest Mack number in the high temperature region
         do i = 2,nmod1
            mack(i) =  abs(u(i))/sqrt(pw53*p(i)/ro(i))
            if (mack(i).gt.Mackmax.and.t(i).gt.1.d8) then
               mackmax = mack(i)
               imack = i
            endif
         enddo

c..   determine shock fronts (T > Tshock and Mack > Mackmin)
         inshock = nmod1
         outshock = nmod1
         do i = nmod1,2,-1
            if (m(i).lt.mcut) goto 50
            if (mack(i).gt.mackmin) outshock = i
            if (t(i+1).gt.tshock.or.t(i-1).gt.tshock) then
               if (inshock.eq.nmod1.and.mack(i).gt.mackmin) inshock = i
               if (t(i+1).lt.tshock.and.t(i).ge.tshock.and.ishockt
     &              .eq.1.and.mack(i).gt.mackmin) ishockt = i+1
               if (t(i+1).gt.tshock.and.t(i).le.tshock.and.mack(i).gt.
     &              mackmin) ishockb = i
            endif
         enddo
 50      if (ishockt.ne.1) then
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
      endif


*__________________________________________________________________
***   compute of the factor ft and fp for the treatment of rotation
*------------------------------------------------------------------

      if (hydrorot) then
         call potential
      else
         do k = 1,nmod
            ft(k) = 1.d0
            fp(k) = 1.d0
         enddo
      endif

      if (rcrz(1).gt.0.and.rotation) Dtot(1) = Dtot(2)

*--------------------------
***   compute the structure
*--------------------------
      ireset = 0
      call thermo (error)

*____________________________________________________
***   Compute diffusion coefficients
***   coefficients that are independent of abundances
*----------------------------------------------------
c$$$      if (microdiffus.or.difftacho) then
c$$$         call diffusion (2,error) ! comment by Thibaut Dumont Jan.2019
c$$$      endif
      do j = 1,nmod
         dturb(j) = Dconv(j)+Dsc(j)+Dherw(j)+Dthc(j)+coefDtacho(j)
     $        +Dturbul(j)+Dbar(j)+Dkyle(j)           ! add Dturbul TD 09/2019
      enddo
      
C Modif CC (01/02/07) Thermohaline
*______________________________________________________
***   Compute atomic diffusion coefficient and velocity
***   coefficients that are independent of abundances
*------------------------------------------------------

c      if (xsp(1,ih1).gt.1.d-5) then 
c          call diffpaq(ihe4,nmod,1,vdmic,dmic,xsp,ratneh)
c      endif

C Modif CC ondes (12/04/07) -->
*______________________________________________________
***   Compute the diffusion coefficient
***   due to internal gravity waves
*------------------------------------------------------
      ndt = novlim(klenv,3)-1
      ndb = novlim(klcore,4)+1

      call coefondes_surf(ndb,ndt)
C <--

*________________________________________________
***   compute of the gravitational binding energy
*------------------------------------------------

      ebind = 0.d0
      epot = 0.d0
      eint = 0.d0
      ebindenv = 0.d0
      epotenv = 0.d0
      eintenv = 0.d0
      do k = 1,nmod-1
         dm(k) = mr(k+1)-mr(k)
         eint = eint+e(k)*dm(k)
         epot = epot+g*m(k)/r(k)*dm(k)
         if (nsconv.gt.0.and.k.ge.novlim(nsconv,3)) then
            eintenv = eintenv+e(k)*dm(k)
            epotenv = epotenv+g*m(k)/r(k)*dm(k)
         endif
      enddo
      ebind = eint-epot
      ebindenv = eintenv-epotenv
      write (*,110) ebind,eint,epot
      write (*,111) ebindenv,eintenv,epotenv
 110  format (/,' STAR     : Ebind =',1pe11.4,' erg, Eint =',1pe11.4,
     &     ' erg, Epot =',1pe11.4,' erg')
 111  format (' envelope : Ebind =',1pe11.4,' erg, Eint =',1pe11.4,
     &     ' erg, Epot =',1pe11.4,' erg')


*--------------------------------------
***   compute of the correlation factor
*--------------------------------------

      do k = 1,nmod
         do l = 1,nsp
            y(l) = max(1.d-50,xsp(k,l))/anuc(l)
         enddo
         zx = 0.d0
         zbar = 0.d0
         ztilde = 0.d0            
         do i = 1,nis
            zx = zx+y(i)
            zy = znuc(i)*y(i)
            zbar = zbar+zy
            ztilde = ztilde+znuc(i)*zy
         enddo
c..   note zx = muiinv(k)
         zbar = zbar/zx
         ztilde = sqrt(ztilde/zx+rhpsi(k)*zbar)
         corr(k) = 1.88d8*sqrt(ro(k)*zx/t(k)**3)
      enddo


*-------------------------------------------
***   useful quantities for angular momentum
*-------------------------------------------

      if (rotation) then
c         print *, 'calc rot'
c remove too little value for xlambda
         do i = 1,nmod
            vxpsi_save(i) = vxpsi(i)
            if (dabs(xlambdas(i)).lt.1d-99) xlambdas(i) = 0.d0
         enddo

c compute terms for Ur
         call diffinit (1)
c.. not done in diffinit
         do i = 1,nmod1
            rro(i) = pw34*dm(i)*msun*totm/(pi*(r(i+1)**3-r(i)**3))
         enddo

! Richardson number Added LR 21 Jan 2025
         rkonv_o(1) = abla(1)-abad(1)
         do i = 2,nmod
            rkonv_o(i) = abla(i)-abm(i)
         enddo
         do i = 1,nmod
            xNt_o(i) = -grav(i)*gs/hp(i)*rkonv_o(i)*deltaKSm(i)
            xNmu_o = grav(i)*gs/hp(i)*abmu(i)*phiKSm(i)
            xN_o(i) = xNt_o(i)+xNmu_o(i)

         enddo
         do i = 1,nmod-1 
            print *,'lumwave =',lumwave(i)
            tamp = (vomega(i+1)-vomega(i))/(r(i+1)-r(i))
            !Ri_n(i) = xN_o(i)/(r(i)*5.d1*tamp)**2 !Previus Version, Why 5.d1, Should be 0.5?
            Ri_n(i) = xN_o(i)/(r(i)*tamp)**2
c            Ri(i) = 0.d0
            if (Ri_n(i).gt.1.d30) Ri_n(i) = 0.d0
            if (Ri_n(i).lt.-1.d30) Ri_n(i) = 0.d0     
cc            print *,i,Ri(i),xN_o(i),tamp,vomega(i),r(i)
         enddo
         Ri_n(nmod) = Ri_n(nmod-1)

! S_NUMBER Added LR 21 Jan 2025
         
         do i = 1,nmod-1
            if (xN_o(i) .gt. 0.d0) then
                S_n(i) = 2*vomega(i)/sqrt(xN_o(i))
            else
               S_n(i) = 0.d0
            endif

            if (S_n(i).gt.1.d30) S_n(i) = 0.d0
            if (S_n(i).lt.-1.d30) S_n(i) = 0.d0   

cc           
         enddo
         S_n(nmod) = S_n(nmod-1)



         vomega(1) = vomega(2)
         omega_S = vomega(nmod)
         dt = dtn/rtot
         xnorm = rtot*omega_S*omega_S/(gs*gs)
         call val_moy
         do j = 2,nmod     

            j1 = j-1
c            print *,j,vxpsi(j),xpsis(j),vxpsi(j)-xpsis(j)
            vxpsim = 0.5d0*(vxpsi_save(j)+vxpsi_save(j1))
            xpsim = 0.5d0*(xpsis(j)+xpsis(j1))
            urm = 0.5d0*(urs(j)+urs(j1))
            auxm = 0.5d0*(auxs(j)+auxs(j1))
            omm = 0.5d0*(vomega(j)+vomega(j1))/omega_S
            Dholdm = 0.5d0*(Dhold(j)+Dhold(j1)) 
            omm2 = omm*omm
            tamp = hnenn(j)


            if (idiffty.eq.13.or.idiffty.eq.15) then
               tamp6a = rmassm(j)*cp(j1)*t(j1)/(dt*rtot*xlumlm(j))
               tamp6 = tamp6a*(xpsim-vxpsim)
               egravj1 = egrav(j1)
               rkonvmj = rkonvLm(j)
            else
               tamp6a = 0.d0
               tamp6 = 0.d0
               egravj1 = 0.d0
               rkonvmj = rkonvm(j)
           endif
               
            if (rcrz(j).lt.0) then
               aalpha = 0.d0
            else
               divis = -gravm(j)*rmassm(j)*cp(j1)*t(j1)*rkonvmj
               aalpha = xlumlm(j)*p(j1)*xnorm/divis
            endif

            raprho = rro(j1)/rhmoy(j)
            if (abs(enucl_save(j1)).gt.1.d-40) then
               epsmoym = 0.5d0*(epsmoy(j)+epsmoy(j1))
               rapeps = (enucl_save(j1)+egravj1)/epsmoym
               fepsi = enucl_save(j1)/(enucl_save(j1)+egravj1)
            else
               rapeps = 0.d0
               fepsi = 0.d0
            endif

            bbeta = 2.d0*(pw43-raprho)

            tamp1a = omm2*omega_S*omega_S*2.384989d6/rro(j1)
            tamp1a = tamp1a-vomega(1)*vomega(1)*2.384989d6/rro(1)

            tamp1b = pw23/raprho*deltaKS(j1)*omega_S*omega_S
            tamp1 = 1.d0-tamp1a+tamp1b*xpsim-rapeps
            tamp2a = rapeps*(fepsi*(epsit(j1)-deltaKS(j1))+deltaKS(j1))
            tamp2b = (2.d0*hhtm(j)/rraym(j)*(1.d0+Dholdm/xKt(j1))-pw23*
     &           deltaKS(j1))/raprho
            tamp2 = tamp2a-tamp2b
            tamp3 = (pw13/raprho)*rraym(j)*tamp
            tamp4 = bbeta*rraym(j)*omm/gravm(j)

            StAdv(j) = -urm*ur_Ss
            StB(j) = aalpha/rro(j)*(tamp4*omm*tamp1+
     &           pw23*deltaKS(j1)*xpsim/raprho) 
            StTh(j) = aalpha/rro(j)*(tamp3*(auxs(j1)-auxs(j))-
     &           2.d0*hhtm(j)/(rraym(j)*rtot)*(1.d0+Dholdm/xKt(j1))
     &           *xpsim/raprho)
            StNG(j) = aalpha/rro(j)*(rapeps*auxm+tamp2a*xpsim)
            StNS(j) = aalpha/rro(j)*(tamp6)

c    

            StB1(j) = aalpha/rro(j)*(tamp4*omm*(1.d0-tamp1a))
            StB2(j) = aalpha/rro(j)*(tamp4*omm*(tamp1b*xpsim))
            StB3(j) = aalpha/rro(j)*(tamp4*omm*(-rapeps))
            StB4(j) = aalpha/rro(j)*(pw23*deltaKS(j1)*xpsim/raprho)

            StTh1(j) = aalpha/rro(j)*(tamp3*(auxs(j1)-auxs(j))
     &           -2.d0*hhtm(j)/(rraym(j)*rtot)*xpsim/raprho)
            StTh2(j) = -aalpha/rro(j)*(2.d0*hhtm(j)/(rraym(j)*rtot)*
     &           (Dholdm/xKt(j1))*xpsim/raprho)

            if (zgradmu) then
               xlambdam = 0.5d0*(xlambdas(j)+xlambdas(j1))
               xxnorm = xnorm*gs

               tamp4 = bbeta*rraym(j)*omm/(omega_S*gravm(j))

* variable Ur
               tamp5a = -pw23*phiKS(j1)/raprho
               tamp5b = rapeps*(fepsi*(epsimu(j1)+phiKS(j1))-phiKS(j1))
               tamp5 = (tamp5a+tamp5b)*xlambdam
               tamp6 = -pw23*phiKS(j1)/raprho*omega_S*omega_S

               StB(j) = StB(j)+aalpha/rro(j)*(tamp4*omm*tamp6
     &              *xlambdam+tamp5a*xlambdam)
               StNG(j) = StNG(j)+aalpha/rro(j)*(tamp5b*xlambdam)

               StB4(j) = StB4(j)+aalpha/rro(j)*(tamp5a*xlambdam)
               StB2(j) =StB2(j)+aalpha/rro(j)*(tamp4*omm*tamp6*xlambdam)

            endif

            aalpha = rro(j1)*t(j1)*cp(j1)*gravm(j)*deltaKS(j)*rkonvm(j)/
     &           hhp(j)
            StAdv(j) = aalpha*StAdv(j)
            StB(j) = aalpha*StB(j)
            StTh(j) = aalpha*StTh(j)
            StNG(j) = aalpha*StNG(j)
            StNS(j) = aalpha*StNS(j)
            
            StB1(j) = aalpha*StB1(j)
            StB2(j) = aalpha*StB2(j)
            StB3(j) = aalpha*StB3(j)
            StB4(j) = aalpha*StB4(j)
            StTh1(j) = aalpha*StTh1(j)
            StTh2(j) = aalpha*StTh2(j)


C Normalisatuon du psi
            xpsis(j) = xpsis(j)/gs*(rtot*omega_S**2)

         enddo
         StAdv(1) = 0.d0
         StB(1) = 0.d0
         StTh(1) = 0.d0
         StNG(1) = 0.d0
         StNS(1) = 0.d0
         StB1(1) = 0.d0
         StB2(1) = 0.d0
         StB3(1) = 0.d0
         StB4(1) = 0.d0
         StTh1(1) = 0.d0
         StTh2(1) = 0.d0
      endif

*----------------------------------------
***   compute turnover-timescales for CE
*----------------------------------------

      if (novlim(klenv,3).ne.0) then
         kbce = novlim(klenv,3)
         hpbce = hp(kbce)

         if (0.05d0*hp(novlim(klenv,3)).lt.r(nmod)) then
            do while (r(kbce).lt.r(novlim(klenv,3))+
     &           0.05d0*hp(novlim(klenv,3)).and.kbce.lt.nmod)
               kbce = kbce+1
            enddo
         endif
         rbce = r(kbce)
         
         ktce = novlim(klenv,4)
         hptce = hp(ktce)
         if (nphase.gt.2) then
            factdc = 0.05d0
         else 
            factdc = 0.1d0
         endif
         do while (r(ktce).gt.r(novlim(klenv,4))-
     &        factdc*hp(novlim(klenv,4))
     &        .and.kbce.gt.1)
            ktce = ktce-1
         enddo
      
C turnover time at Hp/2 over bottom pf CE
         k = novlim(klenv,3)
         do while (r(k).lt.r(novlim(klenv,3))+0.5d0*hp(novlim(klenv,3))
     &        .and.k.lt.nmod)
            k = k+1
         enddo
         print *,'tttt',k,novlim(klenv,3),novlim(klenv,4)
         if (k.eq.0) k = 1
         if (k.ge.novlim(klenv,4).or.sconv(k).eq.0.d0) then
            tc = 0.d0
            rc = 0.d0
         else
            tc = alphac*hp(k)/sconv(k)
            rc = r(k)
         endif

c turnover time at Hp over bottom pf CE
         k = novlim(klenv,3)
         do while (r(k).lt.r(novlim(klenv,3))+hp(novlim(klenv,3))
     &        .and.k.lt.nmod)
            k = k+1
         enddo
         if (k.eq.0) k = 1
         if (k.ge.novlim(klenv,4).or.sconv(k).eq.0.d0) then
            tc_hp = 0.d0
            rc_Hp = 0.d0
         else
            tc_hp = alphac*hp(k)/sconv(k)
            rc_hp = r(k)
         endif
         
c turnover time at 1/2R_ce
         rce = r(novlim(klenv,4))-r(novlim(klenv,3))
         k = novlim(klenv,3)
         do while (r(k).lt.r(novlim(klenv,3))+0.5d0*rce.and.
     &        k.lt.novlim(klenv,4))
            k = k+1
         enddo
         if (k.eq.0) k = 1
         if (k.eq.novlim(klenv,4)) k = novlim(klenv,4)-1
         if (sconv(k).eq.0.d0) then
            tc_hp = 0.d0
            rc_Hp = 0.d0
         else
            tc_r = alphac*hp(k)/sconv(k)
            rc_r = r(k)
         endif

c turnover time at 1/2M_ce
         mce =  m(novlim(klenv,4))-m(novlim(klenv,3))
         k = novlim(klenv,3)
         do while (m(k).lt.m(novlim(klenv,3))+0.5d0*mce.and.
     &        k.lt.novlim(klenv,4))
            k = k+1
         enddo
         if (k.eq.0) k = 1
         if (k.eq.novlim(klenv,4)) k = novlim(klenv,4)-1
         if (sconv(k).eq.0.d0) then
            tc_hp = 0.d0
            rc_Hp = 0.d0
         else
            tc_m = alphac*hp(k)/sconv(k)
            rc_m = r(k)
         endif

c turnover time with maximal value (bottom and top layers are not taken into account)
         tc_max = -1.d0
         do k = kbce,ktce
            if (sconv(k).gt.0.d0.and.
     &           tc_max.lt.alphac*hp(k)/sconv(k)) then
               tc_max = alphac*hp(k)/sconv(k)
               rc_max = r(k)
            endif
         enddo
         if (rc_max.lt.1.d-30) rc_max = 1.d-30
         if (tc_max.eq.-1.d0) then
            tc_max = 0.d0
            rc_max = 0.d0
         endif
      

c integreted turmover time (bottom and top layers are not taken into account)
         tg = 0.d0
         do k = kbce,ktce
            if (sconv(k).gt.0.d0) tg = tg+(r(k+1)-r(k))/sconv(k)
         enddo
         
         secd = 3600.d0*24.d0
         write (*,6666) model,time*seci,tc/secd,rc/rsun,tc_hp/secd,
     &        rc_hp/rsun,tc_r/secd,rc_r/rsun,tc_m/secd,rc_m/rsun,
     &        tc_max/secd,rc_max/rsun,tg/secd
 6666    format('CE-TURNOVER ',i8,1x,1pe16.10,11(1x,1pe10.4))
      else
          write (*,6666) model,time*seci,0.d0,0.d0,0.d0,
     &        0.d0,0.d0,0.d0,0.d0,0.d0,
     &        0.d0,0.d0,0.d0
      endif
      
*--------
* fréquence de brunt Vaisalla
*--------
      
	Nt2=0.d0
	Nu2=0.d0
	bruntV2=0.d0
	bruntV=0.d0
c 	rkonv=0.d0
	
c	rkonv(1)=abla(nmod)-abad(nmod)
c      	do i = 2,nmod
c         j=nmod-i+1
c         rkonv(j) = abla(i)-dsqrt(abad(i)*abad(i-1))
c     	enddo
      
	do i=1,nmod-1
c		gmr_oo(i) = g*m(i)/(r(i)*r(i))
		grav(i) = gmr(i)/gmr(nmod)
		gs = gmr(nmod)
c		hp(i) = pm(i)/rom(i)/abs(gmr(i)+accel(i))
		
		Nt2(i) = gs*grav(i)*deltaKS(i)*(abad(i)-abla(i))/(hp(i))
		Nu2(i) = gs*grav(i)*phiKS(i)*abmu(i)/(hp(i))
		
c		write(*,*) 'Hp'
c		write(*,*) hp(i)
	
		bruntV2(i) = Nt2(i)+Nu2(i)
c		write(*,*) bruntV2(i)

	enddo
		
	do i=1,nmod-1	
           if (bruntV2(i).gt.0.d0.and.crz(i).ne.-2) then
              
              bruntV(i) = dsqrt(bruntV2(i))
              
           elseif (bruntV2(i).le.0.d0.and.bruntV2(i-1).gt.0.d0
     &             .and.bruntV2(i+1).gt.0.d0
     &             .and.crz(i).ne.-2) then
              
              bruntV(i) = dsqrt(-bruntV2(i))
              bruntV2(i) = -bruntV2(i)
           else  		
              bruntV(i) = 0.d0
              bruntV2(i) = 0.d0
           endif 
	enddo
	
      
C...fin modif Nadège
*----------------------------------------
***   compute asteroseimo-stuff to check
*----------------------------------------

      sum = 0.d0
      do k = 1,nmod-1
         sum = sum + (r(k+1)-r(k))/cs(k)
      enddo

c.... grande séparation avec relation asymptotique
      Dnu = 1.d0/(2.d0*sum)

c.... grande séparation avec les relations d'échelle
      Dnuech = 134.9d-6*(m(nmod)/msun)**0.5d0*(reff/rsun)**(-1.5d0)
      nu_max = 3150.d-6*(m(nmod)/msun)*(reff/rsun)**(-2.d0)*
     &  (teff/5780)**(-0.5d0)

c      print *, nu_max

c.... erreur entre les relations asymptotiques et d'échelle

      errr = (Dnu-Dnuech)/Dnu

c....rayon acoustique totale 
      T_tot=1/(2*Dnu)
      
c      write(*,*) 'toto'
      
c... plaçons nous à la base de l'env. convective      
      sum_BCE = 0.d0

      BCE = novlim(klenv,3)
      
      do k = 1,BCE
      sum_BCE=sum_BCE+(r(k+1)-r(k))/cs(k)
      enddo
      
c... rayon acoustique à la base de l'env.conv
      t_BCE=sum_BCE  
      
c... Rayon acoustique a la base de la zone d'ionisation de l'He ++
 
	sum_He = 0.d0
	
	do k= 1,nmod-1
	 
	 if (gamma1(k).lt.1.55.and.t(k).le.1e5.and.
     &      gamma1(k).lt.gamma1(k+1)) then 
     		k_He=k
		goto 23
	 endif
	enddo
	
23	print *, 'kHe', k_He

  	do k = 1,k_He
      	sum_He=sum_He+(r(k+1)-r(k))/cs(k)
      	enddo
	
	r_He=r(k_He)
	t_He=sum_He
	write(*,*) r_He/rsun
      
c.... séparation en période entre les modes g pour l=1   
    
      sum_Dp=0.d0
      Dp=0.d0
c      addvar=0.d0
c      addsum=0.d0
     
c      write(*,*) 'floup'

c      do k=2,nmod-1
c      do k = nmod-1,2,-1 
c      do k=594,nmod-1
c	write(*,*) 'tot2222'

c	if ((2.d0*pi*nu_max)**(2.d0).lt.bruntV2(k).and.
c     &   (2.d0*pi*nu_max)**(2.d0).lt.(2.d0*(cs(k)/r(k))**(2.d0))) then 
     
c     	addvar(k) = 1.d0

c	else 
	
c	addvar(k) = -1.d0
     
c     	endif
		
c	if (addvar(k)*addvar(k+1).lt.0.d0
c     &       .and.(crz(k+1).ne.-2
c     &       .or.k.eq.novlim(klenv,3)
c     & 	     .or.k.eq.novlim(klenv-1,3)
c     &       .or.bruntV2(k+1).eq.0.d0)) then
	
c		addsum=addsum+1
		
c	elseif (addvar(k)*addvar(k+1).lt.0.d0
c     &       .and.crz(k+1).eq.-2.and.
c     &      (k.ne.novlim(klenv,3)
c     &     .or.k.ne.novlim(klenv-1,3)))then
     
c      write(*,*) 'zone conv'	
c	write(*,*) k
	
c		addsum=-abs(addsum)
c	stop
	
c	elseif (addvar(k)*addvar(k+1).lt.0.d0
c     &          .and.Dconv(k+1).le.0.d0) then
     
c	write(*,*) 'totooooooooo'	
	
c		addsum=addsum+1
	
c	endif
c	
c	if (addvar(k).eq.1.d0.and.addsum.gt.0.d0) then 
c	write(*,*) 'toto'
c	write(*,*) 'k=',k

	
      do k = 2,nmod
      
      	if (Dconv(k-1).eq.0.d0.and.Dconv(k).gt.0.d0) then
	
	kconv=k
	goto 25
	
	endif
      enddo  
     
25    if (kconv.eq.novlim(klenv,3)) then
	
	kconv=2
	
      endif
	
c	write(*,*) kconv,novlim(klenv,3)	
      pass = 0
	
26    do k=kconv,nmod-1

          if ((2.d0*pi*nu_max)**(2.d0).lt.bruntV2(k).and.
     &   (2.d0*pi*nu_max)**(2.d0).lt.(2.d0*(cs(k)/r(k))**(2.d0))) then

	    sum_Dp=sum_Dp+((r(k+1)-r(k))*bruntV(k)*(r(k))**(-1.d0))
	    
	  endif
	
c	write(*,*) k,addvar(k),addvar(k+1),addsum
	
      enddo
      
c      write(*,*) sum_Dp
	
      if (sum_Dp.lt.1.d-3
c      .and.sum_Dp.gt.0.d0
     &    .and.pass.eq.0
     &     .and.kconv.ne.2) then
c      	print *,"redo sum"
		sum_Dp = 0.d0
		kconv = 2
		pass = 1
		goto 26
      endif	 
	
	write(*,*) sum_Dp	
		
c      write(*,*) sum_Dp
		
      if (sum_Dp.gt.0.d0) then
	  Dp=((2.0d0)**(0.5)*pi**(2))/(sum_Dp)
	else
	  Dp=-1.0d0
      endif
     

      write (*,6667) model,nphase,time*seci,
c     & lum(nmod)/lsun,
     & Dnu,Dnuech,errr,T_tot,
     & t_BCE,t_He,nu_max,Dp
c     & xspr(1,2),xspr(1,5)
 6667 format ('ASTERO ',i6,1x,i1,1x,1pe16.10,
c     &   1x,1pe11.4,
     &   1x,1pe10.4,1x,1pe11.4,1x,1pe11.4,1x,1pe11.4,
     &   1x,1pe11.4,1x,1pe11.4,1x,1pe11.4,1x,1pe11.4)
c     &   1x,1pe12.6,1x,1pe12.6)
     
     

CCC
C  A commenter
CCC
c      stop 'generation of grid files'
CCC
C  A commenter
CCC

*----------------------------------------------
***   generation of ascii files for Super Mongo
*----------------------------------------------
      write (cmodel,'(i7.7)') model
c      call make_char (model,cmodel)
      open (unit=10,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p0',status = 'unknown')
      open (unit=11,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p1',status = 'unknown')
      open (unit=12,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p2',status = 'unknown')
      open (unit=13,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p3',status = 'unknown')
      open (unit=14,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p4',status = 'unknown')
      open (unit=15,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p5',status = 'unknown')
      open (unit=16,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p6',status = 'unknown')
      open (unit=17,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p7',status = 'unknown')
      open (unit=18,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p8',status = 'unknown')
      open (unit=19,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p9',status = 'unknown')

c.. file.p0
      write (10,1205)
      write (10,1210) nmod,totm,lum(nmod)/lsun,r(nmod)/rsun,time,
     &     iaccr,irotbin,nphase,nmod,nsp-1,model,dtn,nsequence,
     &     int(bin_version*1.d2),FeH

c..   banners
      write (11,1211)
      write (12,1212)
      if (microdiffus) then
         write(13,1224)
      else
         write (13,1213)
      endif
      if (version.lt.2.1d0) then
         write (*,*) 'Old network (version < 2.10) detected'
         write (14,2214)
         write (15,2215)
         write (16,2216)
         write (17,2217)
         write (18,2218)
      else
         write (*,*) 'URCA network (version >= 2.10) detected'
         write (14,1214)
         write (15,1215)
         write (16,1216)
         write (17,1217)
         write (18,1218)
      endif
      write (19,1219)
C... Modif Nadège 15/02/10
      write (25,1225)
C...fin modif Nadège 

c.. file.p1
      do i = 1,nmod
            if (vr(i).gt.1.d30) vr(i) = 0.d0
            if (vr(i).lt.-1.d30) vr(i) = 0.d0    
         write (11,2810) i,rcrz(i),r(i)/rsun,t(i),vvro(i),vvp(i),beta(i)
     &        ,eta(i),lnf(i),sr(i),lum(i)/lsun,u(i),mr(i),accel(i),vr(i)
      enddo
      
c.. file.p2
      do i = 1,nmod1
         ip = i+1
         deltamu = log(mui(ip))-log(mui(i))
         deltamu = (mu(i)*muiinv(i))**wi(ip)*(mu(ip)*muiinv(ip))**wj(ip)
     &        *deltamu
         abmu(ip) = deltamu/(log(vvp(ip))-log(vvp(i))) !d ln mu/ d ln P
      enddo
      abmu(1) = abmu(2)
C Modif CC (02/02/07) Thermohaline --> 
      abmax = -1.d99
      abmin = 1.d99
      do i = 1, nmod
         abmax = max(abmax,abmu(i))
         abmin = min(abmin,abmu(i))
         abmuliss(i) = abmu(i)
      enddo

c remove noise
      abmax = abmax*1.d-2
      abmin = abmin*1.d-2
      do i = 1,nmod
         if (abmu(i).gt.0.d0.and.abmu(i).le.abmax) then
             abmuliss(i) = 0.d0
         endif
         if (abmu(i).lt.0.d0.and.abmu(i).ge.abmin) then
             abmuliss(i) = 0.d0
         endif
      enddo

c smooth
      id = 5
      nmodmdix = nmod-10
      do k = 1,10
         call lissage (abmuliss,nmodmdix,10,id)
      enddo 
      do i = 1,nmod
         if (abrad(i).gt.1.d35) abrad(i) = 1.d35
         write (12,3010) i,tau(i),kap(i),mu(i),1.d0/vvmue(i),
     &        abad(i),abmu(i),abmuliss(i),abrad(i),abla(i),cs(i),
     &        cp(i),gamma1(i),vve(i),rpvisc(i)+1.d-20
      enddo
      
C Fin Modif CC (02/02/07) Thermohaline
c..   file.p3
      do i = 1,nmod
         iscr = 0
            if (rcrz(i).eq.-2) iscr = iscr+1
            if (Dsc(i).gt.0.d0.or.abs(rcrz(i)).eq.3) iscr = iscr+10
            if (Dthc(i).gt.0.d0.or.rcrz(i).eq.2) iscr = iscr+100
            if (Dherw(i).gt.0.d0.or.rcrz(i).eq.1) iscr = iscr+1000
            if (tconv(i).gt.1.d37) tconv(i) = 1.d37
c$$$            if (microdiffus) then ! Ajout TD Fev.2018
               if (i.eq.1) then
                  flux = 1.d0
               else
                  flux = frad(i)+fconv(i)
               endif
               write(13,3213) i,xvdiffhe(i),xvdiffhenoc(i) ! Fin ajout TD Fev.2018 + modif TD Avr.2019 
     $              ,xvdiffc12(i),xvdiffc12noc(i)
     $              ,xvdiffc13(i),xvdiffc13noc(i),xvdiffn14(i)
     $              ,xvdiffn14noc(i),xvdiffn15(i)
     $              ,xvdiffn15noc(i),xvdiffo16(i) ,xvdiffo16noc(i)
     $              ,xvdiffo17(i),xvdiffo17noc(i)
     $              ,xvdiffo18(i),xvdiffo18noc(i),xvdiffne20(i)
     $              ,xvdiffne20noc(i),xvdiffna23(i)
     $              ,xvdiffna23noc(i),tconv(i)
     $              ,abs(rsconv(i)),fconv(i)/flux,frad(i)/flux,renucl(i)
     $              ,enupla(i),egrav(i),enunucl(i),Dconv(i),dturb(i)
     $              ,dherw(i),Dsc(i),Dthc(i),Dondes(i) ,Dondes(i)
     $              ,Dturbul(i),Dbar(i),Dkyle(i),coefDtacho(i),iscr
c$$$         else
c$$$            if (i.eq.1) then
c$$$               flux = 1.d0
c$$$            else
c$$$               flux = frad(i)+fconv(i)
c$$$            endif
c$$$C Modif CC ondes (11/07/07)
c$$$            ii = nmod-i+1
c$$$            write (13,3210) i,tconv(i),abs(rsconv(i)),fconv(i)/flux,
c$$$     &           frad(i)/flux,renucl(i),enupla(i),egrav(i),enunucl(i),
c$$$C Modif CC ondes (12/04/07) -->
c$$$     &           Dconv(i),dturb(i),dherw(i),Dsc(i),dmic(i),vdmic(i),
c$$$C Modif CC ondes (11/07/07)
c$$$C     &        Dthc(i),Dondes(i),Dondeschim(i),iscr
c$$$     &           Dthc(i),Dondes(i),Dondes(i),lumwave(ii),iscr
c$$$     &           
c$$$         endif
      enddo

c.. chemical profiles
c.. file.p4-9
      do i = 1,nmod
         do l = 1,nsp
            xspr(i,l) = max(xspr(i,l),1.d-35)
         enddo
      enddo
      do i = 1,nmod
         write (14,3710) i,(xspr(i,l),l = 1,10)
      enddo
      do i = 1,nmod
         write (15,3710) i,(xspr(i,l),l = 11,20)
      enddo
      do i = 1,nmod
         write (16,3710) i,(xspr(i,l),l = 21,30)
      enddo
      do i = 1,nmod
         write (17,3710) i,(xspr(i,l),l = 31,40)
      enddo
      do i = 1,nmod
         write (18,3710) i,(xspr(i,l),l = 41,50)
      enddo
      do i = 1,nmod
         write (19,4110) i,(xspr(i,l),l = 51,nsp-1)
      enddo
      
c.. file.p10
      if (rotation) then
c      call lissage (Ur,1,nmod-1,5)
c      call lissage (Dhold,1,nmod-1,5)
c      call lissage (dift,1,nmod-1,5)
c      call lissage (dife,1,nmod-1,5)
      if (microdiffus) Dtot(1:nsh) = Dtot(1:nsh)+Dmicro(1:nsh)
         open (unit=20,file=trim(dirmongo) // 'DATA/' // 
     &        name(1:lenci(name)) //cmodel//'.p10',status = 'unknown')
         write (20,1220)
         do i = 1,nmod
            write (20,3460) i,vomega(i),urs(i)*ur_Ss,dift(i),dife(i),
     &           Dtot(i),Dhold(i),xpsis(i),xlambdas(i),abmuj(i),
     &           xnuvv(i),xnum(i),xKt(i),auxs(i),StB(i),StTh(i),
     &           StNG(i),StNS(i),StAdv(i),StTh1(i),StTh2(i),
     &           StB1(i),StB2(i),StB3(i),StB4(i),omega(i),Ri_n(i),
     &           S_n(i)!!,thetac(i)
         enddo
      endif

c.. file.p12
      open (unit=22,file=trim(dirmongo) // 'DATA/' // 
     &        name(1:lenci(name)) //cmodel//'.p12',status = 'unknown')
      write (22,1223)
      do k = 1,nmod
         write (22,3463) k,tautime(k,ippg),Ereac(k,ippg),tautime(k
     $        ,ipdg), Ereac(k,ipdg),tautime(k,i2he3),Ereac(k
     $        ,i2he3),tautime(k,ihe3ag), Ereac(k,ihe3ag),tautime(k
     $        ,ibe7beta),Ereac(k,ibe7beta), tautime(k,ili7pa)
     $        ,Ereac(k,ili7pa),tautime(k,ibe7pg), Ereac(k,ibe7pg)
     $        ,tautime(k,ib8beta),Ereac(k,ib8beta), tautime(k
     $        ,ic13pg),Ereac(k,ic13pg),tautime(k,in14pg), Ereac(k
     $        ,in14pg),tautime(k,icpg),Ereac(k,icpg)
      enddo

c.. file.p11
      open (unit=21,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel//'.p11',status = 'unknown')
      eloc(nmod) = 0.d0
      egconv(nmod) = 0.d0
      if (iaccr.gt.0) then
         write (21,1221)
         do i = 1,nmod
            write (21,3461) i,comp(i),heat(i),work(i),etherm(i),
     &           dlumdt(i),eloc(i),egconv(i),corr(i),vfacc(i),
     &           macc(i),eacc(i),eshr(i),tacc(i)
         enddo
      else
         write (21,1222)
         do i = 1,nmod
            write (21,3462) i,comp(i),heat(i),work(i),etherm(i),
     &           dlumdt(i),eloc(i),egconv(i),corr(i)
         enddo
      endif
      
C... Modif Nadège 15/02/10** 
c  fichier sortie des valeurs frequence de Brunt-Vaisala .p13 
c ------------------------------------------------------- 

      open (unit=25,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p13',status = 'unknown')
      write(25,1225)
      do i=1,nmod
         write (25,2222) i,hp(i),gmr(i),gs,grav(i),phiKS(i),deltaKS(i),
     &        abad(i),abla(i),abmu(i),Nt2(i),Nu2(i),bruntV2(i),bruntV(i)
      enddo

C...fin modif Nadège
           

c.. fichier onde .p14
      if (igw) then
         open (unit=26,file=trim(dirmongo) // 'DATA/' // 
     &        name(1:lenci(name)) //cmodel// '.p14',status = 'unknown')
         write(26,2626)
         do k = 1,nmod
            write (26,2627) k,lgdmic(k),lgdturb(k),vom(k),
     &           depottot(k),depottot_surf(k),
     &           depottot_core(k),
     &           depotondestot(k),Dondes(k),Dondeschim(k),
     &           lumwave(k),lumwave_core(k),
     &           lumondes_surf(k),lumondes_core(k),lumondestot(k),
     &           brunt_o(k),vr(k),vxpsi(k),
     &           Dmicro(k),vmicro(k),Dthc(k),phiKS(k),
     &           deltaKS(k)
         enddo
         close (26)
      endif
 
        
      close (10)
      close (11)
      close (12)
      close (13)
      close (14)
      close (15)
      close (16)
      close (17)
      close (18)
      close (19)
      if (rotation) close (20)
      close (21)
      close(25)

      write (*,*)
      write (*,*) 'Files saved as : ',trim(dirmongo) // 'DATA/' //
     &     name(1:lenci(name)) //cmodel
      write (*,*) 'readsi ',name(1:lenci(name)) //cmodel

      write (*,*) 'sm files generated !'

 999  format (i1,1x,i1,1x,a150)

c..   sm files formats

c.. *.p0
 1205 format ('#  nsh',8x,'MM',13x,'LL',10x,'RR',11x,'times',13x,
     &     'iacc irot npha  shell  nxsp   modl   tstep',7x,
     &     'nseq version bintype  FeH')
 1210 format (1x,i5,2x,f14.10,2x,1pe11.4,2x,1pe10.3,2x,1pe21.15,3(4x,i1)
     &     ,2x,i5,2x,i4,2x,i7,2x,1pe9.3,1x,i6,3x,i4,5x,'AGB',2x,0pf6.3)

c.. *.p1
 1211 format ('# nsh yzi',9x,'r',11x,'T',9x,'rho',11x,'P',7x,'beta',4x,
     &     'eta',6x,'lnf',8x,'s',11x,'Lr',9x,'u',11x,'xmr',14x,'accel',
     &     10x,'vr')
 2810 format (1x,i4,2x,i2,1x,0pf14.9,3(1x,1pe11.4),1x,0pf6.4,1x,0pf8.3,
     &     1x,0pf8.3,1x,1pe11.4,1x,1pe11.4,1x,1pe10.3,2(1x,0pf15.13)
     &     ,1x,1pe11.4)

c.. *.p2
C Modif CC (02/02/07) -->
 1212 format ('# nsh',5x,'tau',7x,'kap',8x,'mu',7x,'mue',5x,'abadd',5x,
     &     'abmu',7x,'abmuliss',7x,'abrad',8x,'abla',8x,'cs',9x,'cp',6x,
     &     'gamma1',4x,'eint',7x,'pvisc')
c 3010 format (1x,i4,1x,1pe10.3,1x,1pe10.4,1x,0pf7.4,1x,1pe10.4,1x,
c     &     0pf7.5,1x,1pe10.3,1x,1pe11.4,1x,1pe11.4,1x,1x,1pe11.4,1x,
c     &     1pe10.4,1x,
c     &     1pe10.4,1x,0pf7.4,1x,1pe10.4,1x,1pe10.3)
 3010 format (1x,i4,1x,1pe10.3,1x,1pe10.4,1x,0pf10.6,1x,1pf12.6,1x,           ! Modif TD Mars 2018 (modif pour sortie mue pf10 -> pf12)
     &     0pf7.5,1x,1pe10.3,1x,1pe11.4,1x,1pe11.4,1x,1x,1pe11.4,1x,
     &     1pe10.4,1x,
     &     1pe10.4,1x,0pf7.4,1x,1pe10.4,1x,1pe10.3)                           ! Fin Modif TD Mars 2018
C Fin Modif CC (02/02/07) -->

c.. *.p3
 1213 format ('# nsh',3x,'tconv',6x,'Vconv',6x,'rfconv',5x,'rfrad',6x,
     &     'enucl',7x,'enupla',7x,'egrav',5x,'enunucl',6x,'Dconv',7x,
     &     'dturb',6x,
     &     'Dherw ',6x,' Dsc  ',6x,'Dmicro',6x,'Vmicro',7x,'Dthc',6x,
!C Modif CC ondes (12/04/07) 
!C     &     'Dondes',6x,'Dondeschim',3x,'scz')
!C Modif CC ondes (11/07/07) 
     &     'Dondes',6x,'Dondeschim',6x,'lumwave',3x,'scz')
c Ajout TD Fev.2018      
 1224 format ('# nsh',3x,'VmicroHe',6x,'VmicroHeTP',6x ,'VmicroC12',6x
     $     ,'VmicroC12TP',6x,'VmicroC13' ,6x ,'VmicroC13TP',6x
     $     ,'VmicroN14',6x,'VmicroN14TP',6x ,'VmicroN15',6x
     $     ,'VmicroN15TP',6x, 'VmicroO16' ,6x,'VmicroO16TP',6x
     $     ,'VmicroO17' ,6x,'VmicroO17TP',6x,'VmicroO18' ,6x
     $     ,'VmicroO18TP'6x, 'VmicroNe',6x ,'VmicroNeTP',6x,'VmicroNa'
     $     ,6x,'VmicroNaTP' ,3x,'tconv',6x ,'Vconv',6x,'rfconv',5x
     $     ,'rfrad' ,6x, 'enucl',7x,'enupla',7x ,'egrav',5x,'enunucl',6x
     $     ,'Dconv' ,7x, 'dturb',6x, 'Dherw ',6x ,' Dsc  ',7x,'Dthc',6x
     $     ,'Dondes' ,6x ,'Dondeschim',6x,'Dturbulence',6x,'Dbar',6x
     $     ,'Dkyle',6x,'Dtacho',3x,'scz')
cFin ajout TD Fev.2018      
 3210 format (1x,i4,1x,1pe10.4,1x,1pe10.4,1x,1pe12.3,1x,
!C Modif CC ondes (12/04/07) 
!C    &     1pe9.3,10(1x,1pe11.4),1x,i4)
!C Modif CC ondes (11/07/07) 
     &     1pe12.3,14(1x,1pe11.4),1x,i4)
 3213 format (1x,i4,20(1x,1pe11.4),1x,1pe10.4,1x,1pe10.4,1x,1pe12.3,1x
     $     ,1pe12.3,15(1x,1pe11.4),1x,i4)   ! Add by TD Fev.2018        
c.. *.p10
 1220 format ('# nsh',4x,'omega',8x,'Ucirc',7x,'Dshear',5x,'Dcirc',
     &     6x,'Dtot',8x,'Dh',9 x,' psi ',7x,'xlambda',5x,
     &     'abmuj',6x,'Nuturb',5x,'Numol',8x,'Kt',9x,'aux',8x,
     &     ' SB ',8x,' STh ',7x,' SNG ',7x,' SNS ',7x,' Sadv ',
     &     7x,' STh1',7x,' STh2',7x,' SB1 ',7x,' SB2 ',7x,
     &     ' SB3 ',7x,' SB4 ',7x,'vomega',7x,'Ri',7x,'S')
 3460 format (1x,i4,1x,1pe12.5,1x,1pe12.5,4(1x,1pe10.4),3(1x,
     &     1pe11.4),3(1x,1pe10.4),12(1x,1pe11.4),1x,4(1pe12.5,1x)) ! Add by LR 21 Jan 2025     


c.. *.p11
 1221 format ('# nsh',5x,'comp',9x,'heat',9x,'work',9x,'etherm',6x,
     &     'dlumdt',8x,'eloc',8x,'egconv',8x,'corr',9x,'facc',
     &     7x,'macc',7x,'eacc',7x,'eshr',7x,'tacc')
 3461 format (1x,i4,8(1x,1pe12.5),1x,0pf9.7,1x,1pe12.5,1x,1pe11.4,1x,
     &     1pe11.4,1x,1pe11.5)
 1222 format ('# nsh',5x,'comp',9x,'heat',9x,'work',9x,'etherm',6x,
     &     'dlumdt',8x,'eloc',8x,'egconv',7x,'corr')
 3462 format (1x,i4,7(1x,1pe12.5),2(1x,1pe11.5))
C... Modif nadège 15/02/10 .p13
 1225 format ('# nsh',12x,'Hp',12x,'gmr',12x,'gs',12x,'grav',12x,'phi',
     &    12x,'delta',12x,'abad',12x,'abla',12x,'abmu',12x,'Nt2',12x,
     &     'Nmu2',15x,'N2',15x,'N') 
 2222 format (i8,1x,e15.6,1x,e15.6,1x,e15.6,1x,e15.6,1x,e15.6,
     &   1x,e15.6,1x,e15.6,1x,e15.6,1x,e15.6,1x,e15.6,1x,
     &   e15.6,1x,e15.6,1x,e15.6)
C...fin modif Nadège 


c..   chemical profiles
c..   new version >= 2.10
 1214 format ('# nsh',7x,'n',11x,'H1',11x,'H2',10x,'He3',10x,'He4',10x,
     &     'Li6',10x,'Li7',10x,'Be7',10x,'B8',11x,'Be9')
 1215 format ('# nsh',5x,'B10',10x,'B11',10x,'C12',10x,'C13',10x,'N13',
     &     10x,'C14',10x,'N14',10x,'N15',10x,'O15',10x,'O16')
 1216 format ('# nsh',6x,'O17',10x,'O18',10x,'F18',10x,'F19',10x,'F20',
     &     9x,'Ne20',9x,'Ne21',9x,'Ne22',9x,'Na22',9x,'Ne23')
 1217 format ('# nsh',5x,'Na23',9x,'Na24',9x,'Mg24',9x,'Na25',9x,
     &     'Mg25',9x,'Mg26',8x,'Al26m',8x,'Al26g',9x,'Mg27',9x,'Al27')
 1218 format ('# nsh',5x,'Si28',9x,'Si29',9x,'Si30',10x,'P31',10x,
     &     'S32',10x,'S33',10x,'S34',10x,'S35',9x,'Cl35',10x,'S36')
 1219 format ('# nsh',4x,'Cl36',7x,'Cl37',6x,'heavy',9x,'sumX')
 3710 format (1x,i4,10(1x,1pe12.6))
 4110 format (1x,i4,3(1x,1pe10.4),1x,1pe11.4)

c..   old version < 2.10
 2214 format ('# nsh',5x,'n',10x,'H1',9x,'H2',8x,'He3',8x,'He4',8x,
     &     'Li6',8x,'Li7',8x,'Be7',9x,'B8',8x,'Be9')
 2215 format ('# nsh',4x,'B10',8x,'B11',8x,'C11',8x,'C12',8x,'C13',8x,
     &     'N13',8x,'C14',8x,'N14',8x,'N15',8x,'O15')
 2216 format ('# nsh',4x,'O16',8x,'O17',8x,'O18',8x,'F18',8x,'F19',7x,
     &     'Ne20',7x,'Ne21',7x,'Ne22',7x,'Na22',7x,'Na23')
 2217 format ('# nsh',4x,'Mg24',7x,'Mg25',7x,'Mg26',6x,'Alm26',6x,
     &     'Alg26',7x,'Al27',7x,'Si28',7x,'Si29',7x,'Si30',7x,'Si31')
 2218 format ('# nsh',5x,'P31',7x,'Si32',8x,'P32',8x,'S32',8x,'P33',
     &     8x,'S33',8x,'S34',8x,'S35',7x,'Cl35',8x,'S36')
 1223 format('# nsh',4x,'Tpp',6x,'Epp',6x,'Tpd',6x,'Epd',6x,'T2he3',6x
     $     ,'E2he3',6x,'The3he4',6x,'Ehe3he4',6x,'Tbe7b',6x,'Ebe7b',6x
     $     ,'Tli7p' ,6x,'Eli7p' ,6x,'Tbe7p',6x,'Ebe7p',6x,'Tb8b',6x
     $     ,'Eb8b',6x ,'Tc13p',6x ,'Ec13p',6x ,'Tn14p',6x,'En14p',6x
     $     ,'Tc12p',6x ,'Ec12p')
 3463 format (1x,i4,22(1x,e15.4E3))
 8100 format(4x,'i',4x,'Dmicro',6x,'Vmicro',5x,'Dthc',6x,'phiKS',6x
     $     ,'deltaKS',6x,'Tpp',6x,'Epp',6x,'Tpd',6x,'Epd',6x,'T2he3',6x
     $     ,'E2he3',6x,'The3he4',6x,'Ehe3he4',6x,'Tbe7b',6x,'Ebe7b',6x
     $     ,'Tli7p' ,6x,'Eli7p' ,6x,'Tbe7p',6x,'Ebe7p',6x,'Tb8b',6x
     $     ,'Eb8b',6x ,'Tc13p',6x ,'Ec13p',6x ,'Tn14p',6x,'En14p',6x
     $     ,'Tc12p',6x ,'Ec12p')

 2626 format('# nsh',5x,'lgdmic',6x,'lgdturb',5x,'vom',9x,'depottot',4x
     &     ,'depottot_surf',1x,'depottot_core',1x,'depotondestot',1x
     &     ,'Dondes',6x,'Dondeschim',2x,'lumwave',5x,'lumwave_core',1x
     &     ,'lumondes_surf',1x,'lumondes_core',1x,'lumondestot',1x
     &     ,'brunt_o',5x,'vr',10x,'vxpsi',7x,'Dmicro2',5x,'vmicro2',5x
     &     ,'Dthc',8x,'phiKS',7x,'deltaKS')
 2627 format(i4,22(1x,1pe11.4))

 8101 format (2x,i4,5(1x,1pe11.4),22(1x,e15.4E2))

      
      end


************************************************************************

      INTEGER FUNCTION LENCI (chain)

************************************************************************
*   Last non-blank character position in string a                      *
************************************************************************

      character*(*) chain

      integer i,n

      n = len(chain)
      do i = n,1,-1
         if (chain(i:i).eq.'_') then
            lenci = i
            return
         endif
      enddo

      lenci = 1

      return
      end


************************************************************************

      INTEGER FUNCTION LENCA (chain)

************************************************************************
*   Last non-blank character position in string a                      *
************************************************************************

      character*(*) chain

 
      integer i,n

      n = len(chain)
      do i = 1,n
         if (chain(i:i).eq.' ') then
            lenca = i-1
            return
         endif
      enddo

      lenca = 0

      return
      end
