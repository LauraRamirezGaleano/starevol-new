                        SUBROUTINE   starevol_init
*                                                                      *
*                       Version 3.30   (Mars 2012)                     *
*                                                                      *
*                   by  Lionel Siess                                   *
*                       Ana Palacios                                   *
*                       Thibaut Decressin                              *
*                       Corinne Charbonnel                             *
*                                                                      *
*   originally developed by  Manuel Forestini                          *
*                                                                      *
*                       Institut d'Astronomie et d'Astrophysique       *
*                       I A A - ULB Brussels                           *
*                                                                      *
*                                                                      *
*                ************************************                  *
*                * this code is NOT A FREEWARE and  *                  *
*                * MUST NOT BE DISTRIBUTED unless   *                  *
*                * the developers have given their  *                  *
*                * unanimous authorization          *                  *
*                ************************************                  *
*                                                                      *
************************************************************************

*************  definitive modifications from Version 2.00  *************
*  incorporation of time-dependent convective mixing        (LS) : V2.01
*  version ln(r) -                                   09/02  (LS)
*  treatment of mesh in case of shocks+reformualtion 10/02  (LS) : V2.02
*  reformulation of structure,inteq,surfeq,centeq    11/02  (LS) : V2.03
*  bugs + new starevol.par, jtranspc, diffusc        01/03  (LS) : V2.05
*  change network : URCA pairs included              04/03  (LS) : V2.10
*  change binary files, change u(i)+psi(i), check compatibility
*  vit (Li7), mesh (lmax), network, prvar, surfeq    06/03  (LS) : V2.13
*  mesh (conservation of J) + file evolvar12         06/03  (LS) : V2.13
*  mesh : enucmax
*  calcevo and rinimod: simplification of initial
*  model reading management                          09/03  (LS) : V2.13
*  new logical for mixing                                        : V2.15
*  partialmix+diffusion+tridiag+mixopt+dtnmix        16/10  (LS) : V2.16
*  change nwtraf, define ledoux regions, writings in init, diffusion
*  removal of mixe, rpaq, demix, abundinit, mixatm, case hpmlt=2
*  change logics of mix, bug in diffusc, surfeq      24/10  (LS) : V2.20
****  synchronization with starevolc_sn_2.20 ****
*  finite difference and new version for diffusion   13/01  (LS) : V2.22
*  mesh improvement (treatment of convective boundaries), new vitsub
*  with additionnel points at low temperature, kapm  20.03  (LS) : V2.30
*  tstepx, mchange + init (dtn < nuclear luminosity) 02.04  (LS) : V2.31
*  correction in convzone (l126 : merging of cz)     13.04  (LS) : V2.32
*  time dep. conv. burning (chedif,nmixd) init.f, dfzconv, totenrg,
*  fzconv, tnuci, diffusion, init                    27.04  (LS) : V2.35
*  semiconvection-ledoux, procconv, convzone         05.05  (LS) : V2.36
*  adaptation for AGB (dtn,maxmod), difledoux        03.06  (LS)
*  generalized version of jtranspc                   11.06  (AP) : V2.37
*  newton raphson for nucleosynthesis, prinit        06.08  (LS) : V2.40
*  effect of rotation on the structure               06.08  (LS) : V2.50
*  bug (mesh, netdiff) fixed 			     05.10  (LS) : V2.51
*  molecular opacity included			     22.10  (TD) : V2.52
*  improve angular momentum conservation (+diffinit) 09.11  (AP) : V2.53
*  2nd order derivatives and new Eq.5               25.01.05(LS) : V2.55
*  new vitsub, starevol.par, neutron burning, mesh   28.01  (LS) : V2.57
*  new opacity tables + numeric + screening factor   21.02  (LS) : V2.58
*  special treatment of unstable nuclei & mixing     28.02  (TD) : V2.60
*  maxshxsp, treatment on neutrons, renormalization  08.03  (LS) : V2.61
*  treatment of convective URCA process, reorganization of evolcom-par
*  bug in chedif, zonetest, convdiff -> dconv,ledoux 30.03  (LS) : V2.65
*  update Fergusson opacity + thermohaline mixing    14.04  (TD) : V2.70
*  update network : include 16O+12C reactions + idup 05.05  (LS) : V2.75
*  fuse evolpar, modif eosc (P ioni), kiinv, nuclopt 02.08  (LS) : V2.76
*  change network B8(g,p), new molecular opacity, generalized
*  jtranspc, new mesh, new logic lmix                22.11 (LSTD): V2.80
*  treatment of atmosphere (xsitau...)               24.12  (LS) : V2.81
*  output nout                                       23.02  (LS) : V2.82
*  network and coupling corrections                  23.02  (LS) : V2.83
*  default nucleosynthesis = one-zone, bug mesh (zonetest), lmix,
*  dtnmnix, HBB adapt., ddmnm                        28.04  (LS) : V2.85
*  new vitsub (al26)                                 08.06  (LS) : V2.90
*  new vitsub (new T grid), jtransp (stabilized), screen (Mitler)
*  opacity (transition), diffusion (omegaconv)       09.11 (ALL) : V2.91
*  change crzc, modif in jtranspc                    09.12  (LS) : V2.92
*  binary modif : irrad + IGW + hpmlt=5 + mesh       17.02  (LS) : V2.95
*  choice between Grevesse/asplund opacities         22.03  (LS) : V2.97
*  bug thermohaline + convzone                       26.09 (LSCC): V2.98
*  new network (D+p)                                 05.05  (LS) : V3.00
*  add IGW transport of AM for enveloppe             21.03  (TD) : V3.10
*  opacity from Asplund et al. 2009                  23.04  (AP) : V3.20
*  asteroseismic data outputs                        04.02  (NL) : V3.30
*  atmospheres, NACRE2 reaction rates, PMS AM transport (torques, disc),
*  new prescriptions for turbulent transport, mass loss for massive
*  stars, several bug corrections                    13.02 (LAAP): V3.40
*                                                                      *
* $LastChangedDate:: 2017-05-23 11:17:15 +0100 (Tue, 23 May 2017)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 94                                                          $ *
*                                                                      *
************************************************************************



      implicit none

c#ifdef _OPENMP
c      include 'omp_lib.h'       !needed for OMP_GET_NUM_THREADS()
c#endif

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.atm'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.data'
      include 'evolcom.mod'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.diff'
      include 'evolcom.igw'
      include 'evolcom.transp'

      character(len=1) calib
      integer n1,n2,n3,i,j,k,ij
      integer mlth,nlth,itlth
      integer filesize
c      integer Ttabeff

      double precision rklth,tklth,opaclth
      double precision rklth_as05,tklth_as05,opaclth_as05
      double precision rklth_as09,tklth_as09,opaclth_as09
      double precision rklth_gn93,tklth_gn93,opaclth_gn93
      double precision rklth_grid,tklth_grid,opaclth_grid
      double precision rklth_ay18,tklth_ay18,opaclth_ay18  ! TD 12/2019
      double precision dturb,dmic,vdmic,lgdturb
     &     ,lgdmic
      double precision rtau,rT,rhotab,Fconvtab
      double precision Ttabeff,gtabeff,Ztabeff
c      double precision gtabeff,Ztabeff
      double precision FeH
      double precision T_Zeff(6),rho_Zeff(6),Fconv_Zeff(6)
      double precision tableTZ(128,59,12),tablerhoZ(128,59,12),
     &     tableFcZ(128,59,12)
      double precision sTZ,srhotabZ,sFconvtabZ
      double precision vom,disctime
      double precision alpharefcalib      ! TD 07/2021

      integer bip1
      double precision bip,taumin,taumax
      character(len=89) ref,nom_fichier
      character(len=len(ref)) DATUM
      character(len=64) fmt_init
      character(len=64) fmt_lect1,fmt_lect2,fmt_lect3
      character(len=5) charT
      character(len=4) charg,charZ
      character(len=21) name,nameg,nameZ,nameinit


      character opalib*200,atmlib*200,stdout*2,file_output*50
      character month*36,GetDateTimeStr*15,fmt*26
      character atmtype*7
       parameter (month = 'JanFebMarAprMayJunJulAugSepOctNovDec')
      parameter (fmt = '(I2.2,A1,I2.2,I3,A3,I4)')
      integer itime(8)

      common /oldstar/ vom(nsh)
      common /coefdiff/ dturb(nsh),dmic(nsh),vdmic(nsh),lgdturb(nsh)
     &     ,lgdmic(nsh)
      common /metal/ FeH

      common /lthopac/ opaclth(19,85,155),tklth(85),rklth(19), ! AGSS09
     &     mlth,nlth,itlth
   !    common /lthopac/ opaclth(19,81,155),tklth(81),rklth(19),
   !   &     mlth,nlth,itlth                                       ! TD pour Opacité Fergu 01/2020
      common /lthopac_as05/ rklth_as05(19),tklth_as05(85),
     &     opaclth_as05(19,85,155)
      common /lthopac_as09/ rklth_as09(19),tklth_as09(85),
     &     opaclth_as09(19,85,155)
      common /lthopac_gn93/ rklth_gn93(19),tklth_gn93(85),
     &     opaclth_gn93(19,85,120)
      common /lthopac_grid/ rklth_grid(19),tklth_grid(63),
     &     opaclth_grid(19,63,104)
   !    common /lthopac_ay18/ rklth_ay18(19),tklth_ay18(85),  ! TD 12/2019
   !   &     opaclth_ay18(19,85,155)
      common /lthopac_ay18/ rklth_ay18(19),tklth_ay18(81),  ! TD pour Opacité Fergu 01/2020
     &     opaclth_ay18(19,81,155)
      common /disclocking/ disctime
!      common /atmospheres/ rtau(127,76,10),rT(127,76,10),
!     &     rhotab(127,76,10),Fconvtab(127,76,10),Ttabeff(76),
!     &     gtabeff(10),Ztabeff(1),nb_geff,nb_Teff,nb_Zeff,filesize

c      common /atmospheres/ rtau(127,10,60,6),rT(127,10,60,6),
c     &     rhotab(127,10,60,6),Fconvtab(127,10,60,6),Ttabeff(60),
c     &     gtabeff(10),Ztbeff(6),taumin,taumax,filesize
      common /atmospheres/ rtau(127,12,59,6),rT(127,12,59,6),
     &     rhotab(127,12,59,6),Fconvtab(127,12,59,6),Ttabeff(59),
     &     gtabeff(12),Ztabeff(6),filesize,taumin,taumax
c$$$      common /atmospheres/ rtau(120,3,17,3),rT(120,3,17,3),
c$$$     &     rhotab(120,3,17,3),Fconvtab(120,3,17,3),Ttabeff(17),
c$$$     &     gtabeff(3),Ztabeff(3),filesize,taumin,taumax
c      common /atmospheres/ rtau(127,60,10,6),rT(127,60,10,6),
c     &     rhotab(127,60,10,6),Fconvtab(127,60,10,6),Ttabeff(60),
c     &     gtabeff(10),Ztabeff(6),filesize


***   define STAREVOL version number
      code_version = 3.30d0
      file_output='00_run.out'

*** Define model atmospheres
      call getenv ('STAR_ATM',atmtype)
      atmtype = trim(atmtype)
      if (atmtype.eq.'none') write(6,*) ' No model atmospheres used'
      

***   standard output unit
c..   6 screen else file
      call getenv('STAR_STDOUT',stdout)
      if (stdout.ne.'96') then
         nout = 6
      else
         call getenv('STAR_FILEOUT',file_output)
         nout = 96
         write (6,*) 'output re-directed to file ',file_output
      endif

      open (unit = 99,file = 'network.par',status = 'old',
     &     action='read')
      open (unit = 90,file = 'evoldisp',status = 'unknown')
      open (unit = 91,file = 'modini.bin',form = 'unformatted',
     &     status = 'unknown')
      open (unit = 92,file = 'nextini.bin',form = 'unformatted',
     &     status = 'unknown')
      if (nout.ne.6) open (unit = nout,file = file_output,
     &     status = 'unknown')

*** units 93 and 94 used for modang and nextang
*** units between 50 and 55 used by opacity tables
*** units between 70 and 78 used by paquette data files
      open (unit = 10,file = 'evolhr',status = 'unknown')
      open (unit = 11,file = 'evolvar1',status = 'unknown')
      open (unit = 12,file = 'evolvar2',status = 'unknown')
      open (unit = 13,file = 'evolvar3',status = 'unknown')
      open (unit = 14,file = 'evolvar4',status = 'unknown')
      open (unit = 15,file = 'evolvar5',status = 'unknown')
      open (unit = 16,file = 'evolvar6',status = 'unknown')
      open (unit = 17,file = 'evolvar7',status = 'unknown')
      open (unit = 18,file = 'evolvar8',status = 'unknown')
      open (unit = 27,file = 'evolvar11',status = 'unknown')
      open (unit = 28,file = 'evolvar12',status = 'unknown')
      open (unit = 29,file = 'evolvar13',status = 'unknown')
      open (unit = 21,file = 'evolchc1',status = 'unknown')
      open (unit = 22,file = 'evolchc2',status = 'unknown')
      open (unit = 123,file = 'evolchc3',status = 'unknown')
      open (unit = 124,file = 'evolchc4',status = 'unknown')
c      open (unit = 25,file = 'evolchc5',status = 'unknown')
      open (unit = 31,file = 'evolchs1',status = 'unknown')
      open (unit = 32,file = 'evolchs2',status = 'unknown')
      open (unit = 33,file = 'evolchs3',status = 'unknown')
      open (unit = 34,file = 'evolchs4',status = 'unknown')
      open (unit = 40,file = 'evolflame',status = 'unknown')
      open (unit = 41,file = 'evoltc1',status = 'unknown')
      open (unit = 42,file = 'evoltc2',status = 'unknown')
      open (unit = 43,file = 'evolas',status = 'unknown')



**************************
*  definition of constants
**************************

      pi = 3.1415926535d0
      pim2 = 2.d0*pi
      pim4 = 4.d0*pi
      pim8 = 8.d0*pi
      sec = 3.1557807d7
      seci = 1.d0/sec
      c = 2.99792458d10
      g = 6.67259d-8
      boltz = 1.380658d-16
      h = 6.6260755d-27
      sig = 2.d0*pi**5*boltz**4/(15.d0*h**3*c*c) ! 5.67050854d-5
c      arad = 4.d0*sig/c  ! radiative density cst : a=7.5659122d-15
      avn = 6.0221367d23
      rk = avn*boltz ! 8.3145112d7
      ech = 4.8032068d-10
      econv = 1.60217733d-6*avn ! 9.648530899d17
      mprot = 1.6726231d-24
      saha = (dsqrt(pim2*boltz*9.1093897d-28)/h)**3 !2.4147031864d15

      msun = 1.9891d33
      rsun = 6.9599d10
      lsun = 3.846d33
      omsun = 2.86d-6

      pw17 = 1.d0/7.d0
      pw14 = 1.d0/4.d0
      pw13 = 1.d0/3.d0
      pw16 = 1.d0/6.d0
      pw23 = 2.d0/3.d0
      pw34 = 3.d0/4.d0
      pw43 = 4.d0/3.d0
      pw53 = 5.d0/3.d0
      pw235 = 2.0d0/35.0d0
      pw635 = 3.0d0*pw235
      vlog10 = log(10.d0)

      phase(1) = "PreHB "
      phase(2) = "HB    "
      phase(3) = "PstHB "
      phase(4) = "HeB   "
      phase(5) = "PstHeB"
      phase(6) = "CB    "
      phase(7) = "NeB   "
      phase(8) = "OB    "
      phase(9) = "SiB   "
      cphase(1) = "D   "
      cphase(2) = "pp  "
      cphase(3) = "CN  "
      cphase(4) = "CNO "
      cphase(5) = "He  "
      cphase(6) = "C   "
      cphase(7) = "Ne  "
      cphase(8) = "O   "
      cphase(9) = "Si  "

************************************************************************
*  read parameter card                                                 *
************************************************************************

      call rmodpar

************************************************************************
*     choice of solar reference and input physics (blockdata evodat.f)
************************************************************************

      call getenv('SOLCOMP',refsolar)
      print *, 'refsolar starevol.f : ', refsolar
** Modif Calib sun - read xspsol from zini.out of previous model - LA 0919
c$$$      call getenv('CALSOL',calib)
c$$$      if (calib.eq.'Y') then
c$$$         print *,'SOLAR CALIBRATION : XSPSOL FROM ZINI.OUT'
c$$$         open (unit = 31,file = 'zini.out',status = 'old')
c$$$         read (31,9002) xspsol(1)
c$$$         do i=2,nis
c$$$            read (31,9003) xspsol(i)
c$$$         enddo
c$$$         xspsol(54) = 5.652d-2*xspsol(53) ! Ar
c$$$         xspsol(55) = 4.197d-2* xspsol(53) ! Ca
c$$$         xspsol(56) = 0.204d-2*xspsol(53) ! Ti
c$$$         xspsol(57) = 1.087d-2*xspsol(53) ! Cr
c$$$         xspsol(58) = 0.708d-2*xspsol(53) ! Mn
c$$$         xspsol(59) = 8.4143d-1*xspsol(53) ! Fe
c$$$         xspsol(60) = 4.664d-2*xspsol(53) ! Ni
c$$$         close (31)
c$$$ 9002    format (/////////,28x,1pe11.5)
c$$$ 9003    format (28x,1pe11.5)
c$$$      else
***   define reference solar composition
c$$$cAP      call getenv('STAR_COMPO',refsolar)
      refsolar = trim(refsolar)
      print *,'refsolar : ',refsolar
      if (refsolar.eq.'GN93') then
         do i = 1,nis+6
            xspsol(i) = xspref(1,i)
         enddo
cAP      if (refsolar.eq.'AC06') then
      elseif (refsolar.eq.'AGS05') then
         alpharefcalib = 1.6046d0
         if (alphac.ne.alpharefcalib) then
            write (nout,*) 'choose the correct value
     $for alphac in parameter card !'
            stop 'starevol.f : bad alpha'
         endif
         do i = 1,nis+6
            xspsol(i) = xspref(2,i)
         enddo
** Modif Asplund 2009 - AP feb 2012 - TD July 2021 -->
      elseif (refsolar.eq.'AGSS09'.and..not.idiffcc.and.ntprof.eq.0
     $        ) then
         alpharefcalib = 1.6602d0
         if (alphac.ne.alpharefcalib) alpharefcalib = 1.6267
c         if (alphac.ne.alpharefcalib) alpharefcalib = 1.7020
c         if (alphac.ne.alpharefcalib) then
c            write (nout,*) 'choose the correct value
c     $for alphac in parameter card !'
c            stop 'starevol.f : bad alpha'
c         endif
         if (alphac.ne.alpharefcalib) alpharefcalib = 1.6304
         do i = 1,nis+6
            xspsol(i) = xspref(3,i)
         enddo
      elseif (refsolar.eq.'AGSS09'.and..not.idiffcc.and.ntprof.eq.4
     $        .and.idiffty.eq.0) then
         alpharefcalib = 1.9730d0
         if (alphac.ne.alpharefcalib) then
            write (nout,*) 'choose the correct value
     $for alphac in parameter card !'
            stop 'starevol.f : bad alpha'
         endif
         do i = 1,nis+6
            xspsol(i) = xspref(5,i)
         enddo
      elseif (refsolar.eq.'AGSS09'.and..not.idiffcc.and.ntprof.eq.5
     $        ) then
         alpharefcalib = 2.1069d0
         if (alphac.ne.alpharefcalib) then
            write (nout,*) 'choose the correct value
     $for alphac in parameter card !'
            stop 'starevol.f : bad alpha'
         endif
         do i = 1,nis+6
            xspsol(i) = xspref(6,i)
         enddo
      
      elseif (refsolar.eq.'AGSS09'.and.ntprof.eq.6
     $        ) then
         alpharefcalib = 1.9730d0
         if (alphac.ne.alpharefcalib) then
            write (nout,*) 'choose the correct value
     $for alphac in parameter card !'
            stop 'starevol.f : bad alpha'
         endif
         do i = 1,nis+6
            xspsol(i) = xspref(6,i)
         enddo
      
      elseif (refsolar.eq.'AGSS09'.and..not.idiffcc.and.ntprof.eq.7
     $        ) then
         alpharefcalib = 1.9730d0
         if (alphac.ne.alpharefcalib) then
            write (nout,*) 'choose the correct value
     $for alphac in parameter card !'
            stop 'starevol.f : bad alpha'
         endif
         do i = 1,nis+6
            xspsol(i) = xspref(6,i)
         enddo
         
      elseif (refsolar.eq.'AGSS09'.and.microdiffus.and.ntprof.
     $        eq.0) then
         alpharefcalib = 1.7936d0
         if (alphac.ne.alpharefcalib) alpharefcalib = 1.7727d0
         if (alphac.ne.alpharefcalib) then
            write (nout,*) 'choose the correct value
     $for alphac in parameter card !'
            stop 'starevol.f : bad alpha'
         endif
         do i = 1,nis+6
            xspsol(i) = xspref(7,i)
         enddo
      elseif (refsolar.eq.'AGSS09'.and.microdiffus.and.ntprof
     $        .eq.5.and.idiffty.eq.0) then
         alpharefcalib = 2.2573d0
         if (alphac.ne.alpharefcalib) then
            write (nout,*) 'choose the correct value
     $for alphac in parameter card !'
            stop 'starevol.f : bad alpha'
         endif
         do i = 1,nis+6
            xspsol(i) = xspref(8,i)
         enddo
      elseif (refsolar.eq.'AGSS09'.and.microdiffus.and.ntprof
     $        .eq.0.and.idiffty.gt.0) then
         print *, 'Warning ! This block data does not exist.'
         stop
c            do i = 1,nis+6
c               xspsol(i) = xspref(9,i)
c            enddo
      elseif (refsolar.eq.'AGSS09'.and.microdiffus.and.ntprof
     $        .eq.4.and.idiffty.gt.0) then
         alpharefcalib = 2.1954d0
         if (alphac.ne.alpharefcalib) then
            write (nout,*) 'choose the correct value
     $for alphac in parameter card !'
            stop 'starevol.f : bad alpha'
         endif
         do i = 1,nis+6
            xspsol(i) = xspref(10,i)
         enddo
      elseif (refsolar.eq.'AGSS09'.and.microdiffus.and.ntprof
     $        .eq.5.and.idiffty.gt.0) then
         alpharefcalib = 2.1987d0
         if (alphac.ne.alpharefcalib) then
            write (nout,*) 'choose the correct value
     $for alphac in parameter card !'
            stop 'starevol.f : bad alpha'
         endif
         do i = 1,nis+6
            xspsol(i) = xspref(11,i)
         enddo
**  Modif Asplund 2009 - AP feb 2012 <--
      elseif (refsolar.eq.'GRID') then
         do i = 1,nis+6
            xspsol(i) = xspref(4,i)
         enddo
**  Modif Young 2018 - TD Dec 2019 - TD July 2021 <--
c      elseif (refsolar.eq.'AY18'.and..not.idiffcc.and.ntprof.eq.5
c     $        ) then
      elseif (refsolar.eq.'AY18'.and.idiffty.eq.0.and.ntprof.eq.5
     $        ) then
         alpharefcalib = 2.1100d0
         if (alphac.ne.alpharefcalib) then
            write (nout,*) 'choose the correct value
     $for alphac in parameter card !'
            stop 'starevol.f : bad alpha'
         endif
         do i = 1,nis+6
            xspsol(i) = xspref(12,i)
         enddo
      elseif (refsolar.eq.'AY18'.and.microdiffus.and.ntprof.eq.5
     $        .and.idiffty.eq.0) then
         alpharefcalib = 2.2600d0
         if (alphac.ne.alpharefcalib) then
            write (nout,*) 'choose the correct value
     $for alphac in parameter card !'
            stop 'starevol.f : bad alpha'
         endif
         do i = 1,nis+6
            xspsol(i) = xspref(13,i)
         enddo
      elseif (refsolar.eq.'AY18'.and.microdiffus.and.ntprof.eq.5
     $        .and.idiffty.gt.0) then
         alpharefcalib = 2.2236d0
         if (alphac.ne.alpharefcalib) then
            write (nout,*) 'choose the correct value
     $for alphac in parameter card !'
            stop 'starevol.f : bad alpha'
         endif
         do i = 1,nis+6
            xspsol(i) = xspref(14,i)
         enddo
      elseif (refsolar.eq.'AY18'.and.ntprof.eq.0) then
         do i = 1,nis+6
            xspsol(i) = xspref(15,i)
         enddo
      elseif (refsolar.eq.'AGSS09') then  !LAURA ADDED 2024/11/19
         alpharefcalib = 2.1954d0
         do i = 1,nis+6
            xspsol(i) = xspref(3,i)
         enddo
         
      else
         stop 'Wrong solar composition, check $SOLCOMP variable'
      endif
c$$$     endif


************************************************************************
*  setting the directory where the opacity tables are to be read       *
*  and reading low-temperature opacity tables                          *
************************************************************************

      if (refsolar.eq.'AGS05') then
         call getenv('DIR_OPAAGS05',opalib)
         mlth = 19
         nlth = 85
         itlth = 155
         print *,'opalib',opalib,mlth,rklth_as05(1)
         do i = 1,mlth
            rklth(i) = rklth_as05(i)
         enddo
         do i = 1,nlth
            tklth(i) = tklth_as05(i)
         enddo
         do i = 1,mlth
            do j = 1,nlth
               do k = 1,itlth
                  opaclth(i,j,k) = opaclth_as05(i,j,k)
               enddo
            enddo
         enddo

**  Modif Asplund 2009 - AP feb 2012 -->
      elseif (refsolar.eq.'AGSS09') then
         call getenv('DIR_OPAAGSS09',opalib) !LIB_AGSS09
         print *,'use new opacity'
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

**  Modif Young 2018 - TD Dec 2019 -->
      elseif (refsolar.eq.'AY18') then
         call getenv('DIR_OPAAY18',opalib) !LIB_AY18
         print *,'use new opacity'
         mlth = 19
          nlth = 85
c        nlth = 81   ! modif TD 30/01/2020 pour opacite Ferguson
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

      print *,'opalib',opalib
      call set_opal_dir(opalib)

************************************************************************
*     Atmospheric tables                                               *
************************************************************************

      nameinit = "Ms02500_g+2.00_z+0.00"
      fmt_lect1 = "(2x,1pe11.5,2x,1pe11.5,2x,1pe11.5,1x,1pe12.5)"
      fmt_lect2 = "(1x,1pe11.5,1x,1pe11.5,1x,1pe11.5,1x,1pe12.5)"
      fmt_init = "(2(/))"
      name = nameinit
      print *, "Atmosphere tables : ", atmtype
      
***   CMFGEN atmospheres
************************************************************************
      
      if (atmtype.eq.'CMFGEN') then
         call getenv ('DIR_ATMC',atmlib)
         call load_CMFGEN(trim(atmlib))
      end if
      
***   TLUSTY atmospheres
************************************************************************

      if (atmtype.eq.'TLUSTY') then
         call getenv ('DIR_ATMT',atmlib)
         call load_TLUSTY(trim(atmlib))
      end if
      
***   PHOENIX atmospheres
************************************************************************

      if (atmtype.eq.'PHOENIX') then
         nameinit = "Ms02500_g+2.00_z-0.00"
         name = nameinit
         call getenv ('DIR_ATMP',atmlib)
         name(1:1) = "P"
         nameg = name
         nb_geff = 12
         nb_Teff = 59
         nb_Zeff = 6
         filesize = 127

         gtabeff(1:nb_geff) = [(i,i=1,nb_geff)]
         gtabeff = gtabeff/2. - 1.
         Ttabeff(1:45) = [(i,i=2600,7000,100)]
         Ttabeff(46:nb_Teff) = [(i,i=7200,9800,200)]

         Ztabeff(1:nb_Zeff) = [-2.d0,-1.5d0,-1.d0,0.d0
     &        ,0.3d0,0.5d0]
         do ij = 1,nb_Zeff
            name = nameg
            write(charZ,'(f3.1)') abs(Ztabeff(ij))
            if (ij.eq.nb_Zeff.or.(ij.eq.(nb_Zeff-1))) name(17:17)='+'
            name(18:20)=charZ
            nameZ = name
            do j = 1,nb_geff
               name = nameZ
               if (j.eq.1) name(10:10)='-'
               write(charg,'(f3.1)') abs(gtabeff(j))
               name(11:13)=charg
               do i = 1,nb_Teff
                  if (i.gt.45.and.gtabeff(j).eq.-0.5d0) then
                     rtau(:,j,i,ij) = 0.d0
                     rT(:,j,i,ij) = 0.d0
                     rhotab(:,j,i,ij) = 0.d0
                     Fconvtab(:,j,i,ij) = 0.d0
                     exit
                  endif
                  if ((Ttabeff(i).eq.10000.d0.and.gtabeff(j).eq.2.d0
     &                 .and.Ztabeff(ij).eq.0.d0).or.
     &                 (Ttabeff(i).eq.10000.d0.and.gtabeff(j).eq.2.d0
     &                 .and.Ztabeff(ij).eq.-2.d0).or.
     &                 (Ttabeff(i).eq.10000.d0.and.gtabeff(j).eq.4.d0
     &                 .and.Ztabeff(ij).eq.-2.d0).or.
     &                 (Ttabeff(i).eq.3000.d0.and.gtabeff(j).eq.3.5d0
     &                 .and.Ztabeff(ij).eq.-3.5d0).or.
     &                 (Ttabeff(i).eq.4000.d0.and.gtabeff(j).eq.4.5d0
     &                 .and.Ztabeff(ij).eq.-3.5d0).or.
     &                 (Ttabeff(i).eq.3900.d0.and.gtabeff(j).eq.4.5d0
     &                 .and.Ztabeff(ij).eq.-3.5d0).or.
     &                 (Ttabeff(i).eq.10000.d0.and.gtabeff(j).eq.2.d0
     &                 .and.Ztabeff(ij).eq.-4.d0).or.
     &                 (Ttabeff(i).eq.10000.d0.and.gtabeff(j).eq.4.d0
     &                 .and.Ztabeff(ij).eq.-4.d0))  then
                     rtau(:,j,i,ij) = 0.d0
                     rT(:,j,i,ij) = 0.d0
                     rhotab(:,j,i,ij) = 0.d0
                     Fconvtab(:,j,i,ij) = 0.d0
                     stop 'Reading wrong tables'
                  else
                     if (Ttabeff(i).lt.10000.d0) then
                        write(charT,'(i4)') int(Ttabeff(i))
                        name(4:7)=charT
                     else
                        write(charT,'(i5)') int(Ttabeff(i))
                        name(3:7)=charT
                     endif
                     ref = trim(atmlib)//trim(name)
                     open(unit = 701,file=ref,
     &                    form="formatted", action="read")
                     read (unit = 701,fmt=fmt_init)
                     do k=1,filesize
                        if (ij.eq.nb_Zeff-1) then
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
      end if
         
***   MARCS atmospheres
************************************************************************
      
      if (atmtype.eq.'MARCS') then
         call getenv ('DIR_ATMM',atmlib)
         name(1:1) = "M"
         nameg = name

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
            name = nameg
            if (gtabeff(j).gt.3.5) name(2:2) = "p"
            do i = 1,nb_Teff
               if ((Ttabeff(i).ge.6000.and.gtabeff(j).lt.3.00).or.
     &              (Ttabeff(i).ge.7000.and.gtabeff(j).lt.4.00).or.
     &              (Ttabeff(i).lt.3000.and.gtabeff(j).ge.5.50)) exit
               write(charg,'(f3.1)') gtabeff(j)
               if (Ttabeff(i).lt.10000) then
                  write(charT,'(i4)') int(Ttabeff(i))
                  name(4:7)=charT
               else
                  write(charT,'(i5)') int(Ttabeff(i))
                  name(3:7)=charT
               endif
               name(11:13)=charg
               ref = trim(atmlib)//trim(name)
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

************************************************************************
*     read initial model                                               *
************************************************************************

      call rinimod (91,93,92,94,0)

************************************************************************
*     Selection of metallicity table
************************************************************************

      if (ntprof.eq.4) call interpZ

************************************************************************
*     determination of time-step, mesh law and mass change             *
************************************************************************

      if (imodpr.eq.11.or.imodpr.eq.12) call myflamespeed(.true.,imodpr)

   10 call init

c.. accretion
      if (iaccr.gt.0) call accret

c.. rezoning
      call mesh
c.. mass loss
      call mchange

************************************************************************
*      calculation of the stellar structure at the next time-step      *
************************************************************************
      
      if (imodpr.eq.11.or.imodpr.eq.12) call myflamespeed(.true.,imodpr)

      call calcevo (*10)

      call date_and_time (values=itime)
      write (GetDateTimeStr,FMT) itime(5),':',itime(6),itime(3),
     &     month(3*itime(2)-2:3*itime(2)),itime(1)
      write (90,'("job ended  : ",a15)') GetDateTimeStr
c      print *,'Dondeschim = ',Dondeschim(nmod)
      close (90)  !  evoldisp
      close (91)  !  modini.bin
      close (92)  !  nextini.bin
      if (irotbin.eq.1) then
         close (93)  !  modang.bin
         close (94)  !  nextang.bin
      endif
      if (nout.ne.6) close (nout)  !  00.log
      close (10)  !  evolhr
      close (11)  !  evolvar1
      close (12)  !  evolvar2
      close (13)  !  evolvar3
      close (14)  !  evolvar4
      close (15)  !  evolvar5
      close (16)  !  evolvar6
      close (17)  !  evolvar7
      close (18)  !  evolvar8
      close (27)  !  evolvar11
      close (28)  !  evolvar12
      close (29)  !  evolvar13
      close (21)  !  evolchc1
      close (22)  !  evolchc2
      close (123)  !  evolchc3
      close (124)  !  evolchc4
c      close (25)  !  evolchc5
      close (31)  !  evolchs1
      close (32)  !  evolchs2
      close (33)  !  evolchs3
      close (34)  !  evolchs4
      close (40)  !  evolflame
      close (41)  !  evoltc1
      close (42)  !  evoltc2
      close (43)  !  evolas

!      if (astero) close (86)
      
      if (irotbin.ne.0.or.microdiffus.or.thermohaline.or.igw) close (85) ! evolvarang

c      close(81)                 ! evoldiff
c      close (85) ! evolvarang

      return
      end
