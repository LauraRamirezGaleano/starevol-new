      PROGRAM runevol

#ifdef GYRE
      use evolstell_gyre
#endif
      implicit none 

      integer i, i_mod, ibeg, iend, filesize
      character (len=5) :: istr, iendstr, ibegstr, bintype, phase
      character(len=200) :: filename, phdir, subinis
      character(len=24) :: starevolog, stopage
      character(len=25) :: filex, hostname, root_name, v_file
      character(len=100) :: stardate
      character(len=10) :: date
      integer error, size_error
      character(len=30) :: seqa, seqaz, seqb, seqbz, seqd, seqdz
      integer sizetot, sizea, sizeaz, sizeb, sizebz, sized, sizedz

      ! Initialize GYRE
#ifdef GYRE
      call evolstell_startup('gyre.in','evolgyre')
#endif
      ! If runevol.e launched by itself
      if ( COMMAND_ARGUMENT_COUNT() == 0 ) then

	  call starevol_init

      ! If it is launched by batchesun_3.40
      else if ( COMMAND_ARGUMENT_COUNT() == 10 ) then
      
         ! Collect arguments
         call getarg(1,filename)
         starevolog=trim(filename)//'.rec'
         call getarg(2,phdir)   ! character string
         call getarg(3,subinis) ! character string         
         call getarg(4,ibegstr) ! convert to integer
         read(ibegstr,'(I5)') ibeg 
         call getarg(5,iendstr) ! convert to integer
         read(iendstr,'(I5)') iend
         call getarg(6,phase)   ! character string
         call getarg(7, stopage) ! character string         
         call getarg(8, bintype) ! get it in char
         call getarg(9, root_name)
         call getarg(10, v_file)
      
         i=ibeg
         do while ( i <= iend )
         
            i_mod=i + 10000

            write(istr,'(I5)') i_mod 
            istr=istr(2:5)
            
            ! Run bash_before         
            call system('$DIR_BATCH/bash_before '//trim(filename)
     &        //' '//trim(ibegstr)//' '//trim(iendstr)//' '//trim(istr))

            ! Retrieve the error from bash_before if there was one
            inquire(file='error.log', size=size_error)
            if ( size_error /=0 ) then
               open(unit=1, file='error.log', action='read') 
               read(1,*) error
               close(1)
               if ( error /= 0) then 
                 call exit(error)
               endif
            endif
            
            ! If the computation is run for the first time, tables must be read
            if ( i == ibeg ) then
               call starevol_init
               
            ! If it is not, tables are already in memory 
            else
               call starevol
      
            end if

            ! Read execution date to use it in bash_after_sun
            open(unit=1, file='date.log', action='read')
            read(1, *) stardate
            close(1)
            
            ! Run bash_after_sun      
            call system('$DIR_BATCH/bash_after_sun '//trim(filename)//
     &        ' '//trim(phdir)//' '//trim(subinis)//' '//trim(ibegstr)//
     &        ' '//istr//' '//trim(iendstr)//' '//trim(stardate)//' '
     &        //trim(phase)//' '//trim(stopage))

            ! Retrieve the error from bash_after_sun if there was one
            inquire(file='error.log', size=size_error)
            if ( size_error /=0 ) then
               open(unit=1, file='error.log', action='read') 
               read(1,*) error
               close(1)
               if ( error /= 0) then
                  call exit(error)
               endif
            endif
            
            ! Run bash_evolfile           
            call system('$DIR_BATCH/bash_evolfile '//
     &        trim(filename)//' '//trim(phdir)//' ' //trim(subinis)//
     &        ' '//trim(root_name)//' '//trim(v_file))
            
            ! Retrieve the error from bash_after_sun if there was one
            inquire(file='error.log', size=size_error)
            if ( size_error /=0 ) then
               open(unit=1, file='error.log', action='read') 
               read(1,*) error
               close(1)
               if ( error /= 0) then
                  call exit(error)
               endif
            endif
            
            ! Go to next iteration
            i=i+1
         end do
         
      ! If it is launched by batche_3.40
      else 
      
         ! Collect arguments
         call getarg(1,filename)
         starevolog=trim(filename)//'.rec'
         call getarg(2,phdir)   ! character string
         call getarg(3,subinis) ! character string
         call getarg(4,ibegstr) ! convert to integer
         read(ibegstr,'(I5)') ibeg
         call getarg(5,iendstr) ! convert to integer
         read(iendstr,'(I5)') iend
         call getarg(6,phase)   ! character string
         call getarg(7, bintype) ! get it in char
         call getarg(8, root_name)
         call getarg(9, v_file)
      
         i=ibeg
         do while ( i <= iend )
         
            i_mod=i + 10000

            write(istr,'(I5)') i_mod
            istr=istr(2:5)
         
            ! Run bash_before
            call system('$DIR_BATCH/bash_before '//trim(filename)
     &        //' '//trim(ibegstr)//' '//trim(iendstr)//' '//trim(istr))

            ! Retrieve the error from bash_before if there was one
            inquire(file='error.log', size=size_error)
            if ( size_error /=0 ) then
               open(unit=1, file='error.log', action='read') 
               read(1,*) error
               close(1)
               if ( error /= 0) then
                  call exit(error)
               endif
            endif
            
            ! If the computation is run for the first time, tables must be read
            if ( i == ibeg ) then
               call starevol_init
            
            ! If it is not, tables are already in memory 
            else
               call starevol
         
            end if
            
            ! Read execution date to use it in bash_after_sun
            open(unit=1, file='date.log', action='read')
            read(1, *) stardate
            close(1)
            
            ! Run bash_after      
            call system('$DIR_BATCH/bash_after '//trim(filename)//
     &        ' '//trim(phdir)//' '//trim(subinis)//' '//trim(ibegstr)//
     &        ' '//istr//' '//trim(iendstr)//' '//trim(stardate)//' '
     &        //trim(phase))
     
            ! Retrieve the error from bash_after if there was one
            inquire(file='error.log', size=size_error)
            if ( size_error /=0 ) then
               open(unit=1, file='error.log', action='read') 
               read(1,*) error
               close(1)
               if ( error /= 0) then
                  call exit(error)
               endif
            endif
       
            ! Run bash_evolfile        
            call system('$DIR_BATCH/bash_evolfile '//
     &        trim(filename)//' '//trim(phdir)//' ' //trim(subinis)//
     &        ' '//trim(root_name)//' '//trim(v_file))
     
            ! Retrieve the error from bash_evolfile if there was one
            inquire(file='error.log', size=size_error)
            if ( size_error /=0 ) then
               open(unit=1, file='error.log', action='read') 
               read(1,*) error
               close(1)
               if ( error /= 0) then
                  call exit(error)
               endif
            endif
            
            ! Go to next iteration
            i=i+1
         end do
      
      end if
      
      ! Finalize GYRE
#ifdef GYRE      
      call evolstell_final
#endif
      end program runevol
      
      

      
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
* modification and correction of atomic diffusion : inclusion of Thoul
* formalism, modification of EoS to account
* for partial ionization                             2018 (TD)   : V3.40
* new solar reference chemical composition (AY18)    2019  (TD)  : V3.40
* inclusion of namelists, modification of evodat     2021  (TD)  : V4.00
* modification of opacity tables call to speed up code 2021 (AF) : V4.00
* coupling with the GYRE oscillation code            2021 (AF,AP): V4.00
* modification of mass loss for massive stars        2021  (AP)  : V4.00
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


       common /lthopac/ opaclth(19,85,155),tklth(85),rklth(19),
     &     mlth,nlth,itlth
c     common /lthopac/ opaclth(19,81,155),tklth(81),rklth(19),
c    &     mlth,nlth,itlth                                       ! TD pour Opacité Fergu 01/2020
      common /lthopac_as05/ rklth_as05(19),tklth_as05(85),
     &     opaclth_as05(19,85,155)
      common /lthopac_as09/ rklth_as09(19),tklth_as09(85),
     &     opaclth_as09(19,85,155)
      common /lthopac_gn93/ rklth_gn93(19),tklth_gn93(85),
     &     opaclth_gn93(19,85,120)
      common /lthopac_grid/ rklth_grid(19),tklth_grid(63),
     &     opaclth_grid(19,63,104)
      common /lthopac_ay18/ rklth_ay18(19),tklth_ay18(85),  ! TD 12/2019
     &     opaclth_ay18(19,85,155)
c      common /lthopac_ay18/ rklth_ay18(19),tklth_ay18(81),  ! TD pour Opacité Fergu 01/2020
c     &     opaclth_ay18(19,81,155)
      common /disclocking/ disctime

      common /atmospheres/ rtau(127,12,59,6),rT(127,12,59,6),
     &     rhotab(127,12,59,6),Fconvtab(127,12,59,6),Ttabeff(59),
     &     gtabeff(12),Ztabeff(6),filesize,taumin,taumax




***   define STAREVOL version number
      code_version = 3.30d0
      file_output='00_run.out'

*** Define model atmospheres
      call getenv ('STAR_ATM',atmtype)
      atmtype = trim(atmtype)
      if (atmtype.eq.'none') write(6,*) ' No model atmospheres used'
        print *,atmtype

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

      open (unit = 99,file = 'starevol.par',status = 'old',
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
      print *, 'refsolar starevol.f', refsolar
      refsolar = trim(refsolar)
      print *,'refsolar',refsolar
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
*     Atmospheric table                                                *
************************************************************************

      nameinit = "Ms02500_g+2.00_z+0.00"
      fmt_lect1 = "(2x,1pe11.5,2x,1pe11.5,2x,1pe11.5,1x,1pe12.5)"
      fmt_lect2 = "(1x,1pe11.5,1x,1pe11.5,1x,1pe11.5,1x,1pe12.5)"
      fmt_init = "(2(/))"
      name = nameinit

*** PHOENIX ATMOSPHERES
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

*** MARCS ATMOSPHERES

      else if (atmtype.eq.'MARCS') then
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

      endif


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

      if (irotbin.ne.0.or.microdiffus.or.thermohaline.or.igw) close (85) ! evolvarang

      return
      end
                        SUBROUTINE   starevol
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
      double precision FeH
      double precision T_Zeff(6),rho_Zeff(6),Fconv_Zeff(6)
      double precision tableTZ(128,59,12),tablerhoZ(128,59,12),
     &     tableFcZ(128,59,12)
      double precision sTZ,srhotabZ,sFconvtabZ
      double precision vom,disctime

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


       common /lthopac/ opaclth(19,85,155),tklth(85),rklth(19),
     &     mlth,nlth,itlth
c     common /lthopac/ opaclth(19,81,155),tklth(81),rklth(19),
c    &     mlth,nlth,itlth                                      ! TD pour Opacité Fergu 01/2020
      common /lthopac_as05/ rklth_as05(19),tklth_as05(85),
     &     opaclth_as05(19,85,155)
      common /lthopac_as09/ rklth_as09(19),tklth_as09(85),
     &     opaclth_as09(19,85,155)
      common /lthopac_gn93/ rklth_gn93(19),tklth_gn93(85),
     &     opaclth_gn93(19,85,120)
      common /lthopac_grid/ rklth_grid(19),tklth_grid(63),
     &     opaclth_grid(19,63,104)
c      common /lthopac_ay18/ rklth_ay18(19),tklth_ay18(85),  ! TD 12/2019
c     &     opaclth_ay18(19,85,155)
      common /lthopac_ay18/ rklth_ay18(19),tklth_ay18(81),  !  TD pour Opacité Fergu 01/2020
     &     opaclth_ay18(19,81,155)
      common /disclocking/ disctime
      common /atmospheres/ rtau(127,12,59,6),rT(127,12,59,6),
     &     rhotab(127,12,59,6),Fconvtab(127,12,59,6),Ttabeff(59),
     &     gtabeff(12),Ztabeff(6),filesize,taumin,taumax




***   define STAREVOL version number
      code_version = 4.00d0
      file_output='00_run.out'

*** Define model atmospheres
      call getenv ('STAR_ATM',atmtype)
      atmtype = trim(atmtype)
      if (atmtype.eq.'none') write(6,*) ' No model atmospheres used'
        print *,atmtype

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

      open (unit = 99,file = 'starevol.par',status = 'old',
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
*  setting the directory where the opacity tables are to be read       *
*  and reading low-temperature opacity tables -> moved to starevol_init*
************************************************************************

************************************************************************
*  read parameter card                                                 *
************************************************************************

      call rmodpar


************************************************************************
*     choice of solar reference and input physics (blockdata evodat.f)
************************************************************************

      call getenv('SOLCOMP',refsolar)
      print *, 'refsolar starevol.f', refsolar
      
      refsolar = trim(refsolar)
      print *,'refsolar',refsolar
      if (refsolar.eq.'GN93') then
         do i = 1,nis+6
            xspsol(i) = xspref(1,i)
         enddo
cAP      if (refsolar.eq.'AC06') then
      elseif (refsolar.eq.'AGS05') then
         do i = 1,nis+6
            xspsol(i) = xspref(2,i)
         enddo
c$$$** Modif Asplund 2009 - AP feb 2012 - TD July 2021 -->
      elseif (refsolar.eq.'AGSS09'.and..not.idiffcc.and.ntprof.eq.0
     $        ) then
         do i = 1,nis+6
            xspsol(i) = xspref(3,i)
         enddo
      elseif (refsolar.eq.'AGSS09'.and..not.idiffcc.and.ntprof.eq.4
     $        .and.idiffty.eq.0) then
         do i = 1,nis+6
            xspsol(i) = xspref(5,i)
         enddo
      elseif (refsolar.eq.'AGSS09'.and..not.idiffcc.and.ntprof.eq.5
     $        ) then
         do i = 1,nis+6
            xspsol(i) = xspref(6,i)
         enddo
      elseif (refsolar.eq.'AGSS09'.and.microdiffus.and.ntprof.
     $        eq.0) then
         do i = 1,nis+6
            xspsol(i) = xspref(7,i)
         enddo
      elseif (refsolar.eq.'AGSS09'.and.microdiffus.and.ntprof
     $        .eq.5.and.idiffty.eq.0) then
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
         do i = 1,nis+6
            xspsol(i) = xspref(10,i)
         enddo
      elseif (refsolar.eq.'AGSS09'.and.microdiffus.and.ntprof
     $        .eq.5.and.idiffty.gt.0) then
         do i = 1,nis+6
            xspsol(i) = xspref(11,i)
         enddo
**  Modif Asplund 2009 - AP feb 2012 <--
      elseif (refsolar.eq.'GRID') then
         do i = 1,nis+6
            xspsol(i) = xspref(4,i)
         enddo
**  Modif Young 2018 - TD Dec 2019 - TD July 2021 <--
      elseif (refsolar.eq.'AY18'.and..not.idiffcc.and.ntprof.eq.5
     $        ) then
         do i = 1,nis+6
            xspsol(i) = xspref(12,i)
         enddo
      elseif (refsolar.eq.'AY18'.and.microdiffus.and.ntprof.eq.5
     $        .and.idiffty.eq.0) then
         do i = 1,nis+6
            xspsol(i) = xspref(13,i)
         enddo
      elseif (refsolar.eq.'AY18'.and.microdiffus.and.ntprof.eq.5
     $        .and.idiffty.gt.0) then
         do i = 1,nis+6
            xspsol(i) = xspref(14,i)
         enddo  
      else
         stop 'Wrong solar composition, check $SOLCOMP variable'
      endif
c$$$      endif
         
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

      if (irotbin.ne.0.or.microdiffus.or.thermohaline.or.igw) close (85) ! evolvarang

c$$$      if (ifail.gt.0) then
c$$$         stop 'failed'
c$$$      else
c$$$         stop 'normal'
c$$$      endif

      return
      end


************************************************************************

      SUBROUTINE accret

************************************************************************
* calculate the mass deposition and angular velocity in each shell due *
* to accretion                                                         *
*                                                                      *
*               PARAMETRIC ACCRETION : accphase values                 *
*                                                                      *
*  0  : accretion model, including D burning (suited for PMS phase)    *
* 1,4 : accretion inside the star (planet accretion)                   *
*  5  : uniform accretion from the surface with facc=ric (menv unknown)*
*  6  : uniform accretion from the surface M* to menv (facc unknown)   *
*  7  : uniform accretion from the surface M* to M*-menv (facc unknown)*
*  8  : uniform accretion from the surface M* to menv*M  (facc unknown)*
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.eng'
      include 'evolcom.grad'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer itnr,itime
      integer iexit,icrasha,icenter,idet
      integer i,j,k
      integer iaccmax,itmax,ishockb,ishockt

      integer ibotC,itopC,ibotH,itopH,itopHe,ibotHe
      integer lb,lt,ideb

      double precision pts,shlim0,Lmax,tmax,Mackmax,enucmax
      double precision schwar(nsh)
      double precision dfdr
      double precision totmv,sinaccr,fdtacc0,dmaccrg,fmax,fmin0,
     &     maccr1,maccr,hx,y,fmin,rmid,dtacc0,dlim,fmean,flinea,fmx,
     &     faccmax,dminf

c     parameter (itnr = 10,pts = 1.d-2)
      parameter (itnr = 4,pts = 1.d-2)

      common /hydrodyn/ Lmax,tmax,Mackmax,enucmax,itmax,ishockb,ishockt

      external dfdr

      if (iaccr.eq.0) goto 30

      faccmax = 1.d0
      totmv = m(nmod)
      dlim = 1.d0
      lacc = 0.d0
      itime = 0
      iaccbot = 2
      iaccmax = 50
      iexit = 0
      idet = 0
      icrasha = 0
      fdtacc0 = fdtacc
      dmaccr = 0.d0
      sigma0 = xiaccr*dsqrt(g*r(nmod)*totmv)

      if (dtn*seci.lt.1.d-2.or.totmv.gt.massend.or.(tmax.gt.1.d8.and.
     &     xsp(itmax,ihe4).gt.0.5d0)) then
         iaccr = 0
         write (nout,*) 'accretion stopped'
         goto 30
      endif

***   initialization of accretion variables

      if (model.lt.0) goto 30

      dtacc0 = dtn
      dmaccr = massrate*dtn*seci
      dmaccrg = dmaccr*msun

*---------------------
***   Planet accretion
*---------------------

      if (accphase.gt.0.and.accphase.le.4) then

***   determination of the nuclear burning region boundaries

         sinthac = 1.d0
         shlim0 = shlim-1.d0
 11      shlim0 = shlim0+1.d0
         itopH = 0
         itopHe = 0
         ibotH = 0
         ibotHe = 0
         itopC = 0
         ibotC = 0
         nnulim(1:nmaxenu,1:3) = 0
         k = 1
         if (enucl(1).ge.shlim0) then
            nnulim(1,1) = 1
            nnulim(1,2) = 1
            do i = 2,nmod
               if (enucl(i).lt.shlim0) goto 101
            enddo
 101        nnulim(1,3) = i
            k = 2
         endif
         do 201 i = 2,nmod
            if (enucl(i-1).lt.shlim0.and.enucl(i).ge.shlim0) then
               nnulim(k,1) = i
               nnulim(k,2) = i
            endif
            if (nnulim(k,2).eq.0) goto 201
            if (enucl(i).gt.enucl(i-1).and.enucl(i).gt.
     &           enucl(i+1).and.enucl(i).gt.enucl(nnulim(k,2)))
     &           nnulim(k,2) = i
            if (enucl(i).lt.shlim0.and.enucl(i-1).ge.shlim0) then
               nnulim(k,3) = i
               k = k+1
               if (k.gt.nmaxenu) goto 301
            endif
 201     continue
 301     ideb = 1
         if (nnulim(1,1).eq.1) ideb = 2
         do 401 i = ideb,nmaxenu
            if (nnulim(i,1).eq.0) goto 501
            if (enucl(nnulim(i,2)).lt.shlim0) goto 401
            lb = nnulim(i,1)
            lt = nnulim(i,3)
            if (xsp(lb,ih1).lt.xsp(lt,ih1).and.xsp(lb,ihe4).gt.
     &           xsp(lt,ihe4)) then
               itopH = nnulim(i,3)
               ibotH = nnulim(i,1)
            endif
            if (xsp(lb,ihe4).lt.xsp(lt,ihe4).and.itopHe.eq.0) then
               ibotHe = nnulim(i,1)
               itopHe = nnulim(i,3)
            endif
 401     continue
 501     if (agbphase.and.k.eq.2.and.itopH.eq.0) goto 11
         do i = nsconv,1,-1
            if (novlim(i,3).ne.0) goto 601
         enddo
 601     if (i.lt.1) i = 1
         ibotC = novlim(i,3)
         itopC = novlim(i,4)
         if (nsconv.gt.1.and.mr(ibotC).gt.0.95d0) then
            ibotC = novlim(i-1,3)
            iacctop = ibotC+5
         endif
         iacctop = ibotC+5
         if (itopH+1.gt.iacctop) iacctop = itopH
         if (accphase.eq.1) then
            if (itopH+1.lt.iacctop) then
               iaccbot = itopH+1
            else
               iaccbot = ibotH
c              iaccbot = ibotC
            endif
            if (itopH.eq.0) then
               iaccbot = ibotC-1
               iacctop = iaccbot+10
            endif
         endif
         if (accphase.eq.2) iaccbot = ibotH
         if (accphase.eq.3) then
            if (itopH.ne.0) then
               iaccbot = itopHe+1
            else
               iaccbot = ibotHe
            endif
            if (iaccbot.gt.iacctop) iaccbot = ibotHe
         endif
         facc(1:nmod) = 0.d0
         macc(1:nmod) = 0.d0
c        if (itacc.eq.3) then
c           fmean = dmaccrg/(sinthac*(m(iacctop+1)-m(iaccbot)))
c           do k = iaccbot,iacctop
c              facc(k) = fmean
c            enddo
c        endif
 701     fmean = dmaccrg/(sinthac*(m(iacctop+1)-m(iaccbot)))
         flinea = 2.d0*fmean/(m(iacctop+1)-m(iaccbot))
         fmx = flinea*((m(iacctop+1)+m(iacctop))*0.5d0-m(iaccbot))
         if (fmx.gt.faccmax.and.mr(iacctop).lt.0.9d0) then
            iacctop = iacctop+10
            goto 701
         endif
         do k = iaccbot,iacctop
            facc(k) = flinea*((m(k+1)+m(k))*0.5d0-m(iaccbot))
            macc(k+1) = macc(k)+dm(k)*facc(k)*sinthac
         enddo
         icenter = 0
         iaccbot = iaccbot-1
         dmaccr = dmaccrg/msun
         raccbot = r(iaccbot)/r(nmod)
         maccbot = m(iaccbot)/m(nmod)
         if (facc(iacctop).gt.1.d0/sinthac) then
            if (iaccr.eq.1) iaccr = 2
            if (iaccr.eq.3) iaccr = 4
         endif
         goto 25
      else
         iacctop = nmod-2
      endif

      if (totmv.gt.6.d0) faccmax = 3.d0

*--------------------------------------------------
***   Novae explosion : pill-up mass on top of star
*--------------------------------------------------

      if (accphase.ge.5) then
         sinthac = 1.d0
***   accrete with f = cste
         if (accphase.eq.5) then
            fmean = abs(ric)
            dminf = m(nmod)-dmaccrg/(sinthac*fmean)
         endif
***   accrete from a fixed mass coordinate (Menv)
         if (accphase.eq.6) then
            if (menv.gt.m(nmod)) then
               write (nout,100)
               stop 'accret : accphase = 6'
            endif
            dminf = menv
            fmean = (m(nmod)-menv)/(sinthac*dmaccrg)
         endif
***   accretion with a constant mass depth
         if (accphase.eq.7) then
            if (menv.gt.m(nmod)) then
               write (nout,100)
               stop 'accret : accphase = 7'
            endif
            dminf = m(nmod)-menv
            fmean = menv/(sinthac*dmaccrg)
         endif
***   accretion with a constant relative mass depth
         if (accphase.eq.8) then
            if (menv.gt.1.d0) then
               write (nout,200)
               stop 'accret : accphase = 8'
            endif
            dminf = m(nmod)*menv
            fmean = m(nmod)*(1.d0-menv)/(sinthac*dmaccrg)
         endif
         do k = nmod,2,-1
            if (m(k).le.dminf) then
               iaccbot = k
               iacctop = nmod
               goto 1
            endif
         enddo
 1       if (accphase.eq.5) iaccbot = min(iaccbot,nmod-10)
         do k = 1,iaccbot-1
            facc(k) = 0.d0
            macc(k+1) = 0.d0
         enddo
         fmean = dmaccrg/(sinthac*(m(nmod)-m(iaccbot)))
c        flinea = 2.d0*fmean/(m(nmod)-m(iaccbot))
c        fmx = flinea*((m(nmod)+m(nmod1))*0.5d0-m(iaccbot))
         do k = iaccbot,nmod1
c           facc(k) = flinea*((m(k+1)+m(k))*0.5d0-m(iaccbot))
            facc(k) = fmean
            macc(k+1) = macc(k)+dm(k)*facc(k)*sinthac
         enddo
         macc(1) = 0.d0
         facc(nmod) = facc(nmod1)
         icenter = 0
         iaccbot = iaccbot-1
         raccbot = r(iaccbot)/r(nmod)
         maccbot = m(iaccbot)/m(nmod)
         dmaccr = dmaccrg/msun
         maccr = macc(nmod)
         time = time-dtn
         dtn = maccr*sec/(massrate*msun)
         time = time+dtn
         dmaccr = dmaccrg/msun
         if (iaccr.eq.1) iaccr = 2
         if (iaccr.eq.3) iaccr = 4
         goto 25
      endif


*-----------------------------
***   Accretion onto PMS stars
*-----------------------------

      sinaccr = dsqrt(1.d0-accrw**rir)
      thacc = 1.d0-sinaccr*sinaccr*pw13
      sinthac = sinaccr*thacc

***   Initializations

      do j = 2,nmod1
         schwar(j) = abad(j)-abla(j)
      enddo
      schwar(1) = schwar(2)
      schwar(nmod) = schwar(nmod1)
      do j = 1,itnr
         do i = 2,nmod1
            if (abs((schwar(i+1)-schwar(i))/schwar(i)).gt.pts.or.
     &           abs((schwar(i)-schwar(i-1))/schwar(i)).gt.pts)
     &           schwar(i) = abs(schwar(i+1)+schwar(i-1))*0.5d0
         enddo
      enddo

 5    if (icrasha.gt.iaccmax) then
         if (iexit.gt.10) stop 'accret'
         write (nout,300) iaccbot,dtn*seci,maccr
         write (90,300) iaccbot,dtn*seci,maccr
         iexit = iexit+1
         time = time-dtn
         if (iexit.eq.6) then
            dtn = dtacc0
            fdtacc0 = 1.d0/fdtacc
         endif
         dtn = dtn/fdtacc0
         time = time+dtn
         dmaccr = massrate*dtn*seci
         dmaccrg = dmaccr*msun
      endif

      iaccbot = 1
      fmax = dmaccrg/(totmv*sinthac)
      fmin0 = fmax/fdtacc
      facc(iaccbot) = 0.d0
      macc(iaccbot) = 0.d0
      icenter = 0
      maccr1 = 0.d0

 10   iaccbot = iaccbot+1
      facc(iaccbot) = 0.d0
      macc(iaccbot) = 0.d0
      macc(iaccbot+1) = 0.d0
      icrasha = 0

 20   facc(1) = facc(iaccbot)
      if (icenter.eq.1) then
         if (iaccr.eq.2.or.iaccr.eq.4) maccr = m(2)*sinthac*facc(1)
         if (iaccr.eq.1.or.iaccr.eq.3) maccr = m(2)*sinthac*facc(1)/
     &        (1.d0-sinthac*facc(1))
      else
         maccr = 0.d0
      endif
      macc(2) = maccr
      idet = 1
      icrasha = icrasha+1
      if (icrasha.gt.iaccmax) goto 5
      do i = iaccbot+1,iacctop
         hx = r(i)-r(i-1)
         y = hx*dfdr (i-1,0.d0,facc(i-1),schwar(i-1),idet)
         facc(i) = facc(i-1)+y/6.d0
         y = hx*dfdr (i-1,hx*0.5d0,facc(i-1)+y*0.5d0,schwar(i-1),idet)
         facc(i) = facc(i)+y/3.d0
         y = hx*dfdr (i-1,hx*0.5d0,facc(i-1)+y*0.5d0,schwar(i-1),idet)
         facc(i) = facc(i)+y/3.d0
         y = hx*dfdr (i-1,hx,facc(i-1)+y,schwar(i-1),idet)
         facc(i) = facc(i)+y/6.d0

         if (((iaccr.eq.1.or.iaccr.eq.3).and.facc(i).gt.1.d0).or.
     &        iaccbot.ge.(nmod-10)) then
            icrasha = 101
            goto 5
         endif
         if (facc(i).lt.0.d0.or.idet.eq.0) then
            goto 10
         endif
         if ((iaccr.eq.2.or.iaccr.eq.4).and.facc(i).gt.faccmax) then
            facc(i) = faccmax
         endif
         if (iaccr.eq.1.or.iaccr.eq.3) maccr = maccr+dm(i)*sinthac*
     &        facc(i)/(1.d0-sinthac*facc(i))
         if (iaccr.eq.2.or.iaccr.eq.4) maccr = maccr+dm(i)*sinthac*
     &        facc(i)
         macc(i+1) = maccr
         if (itime.eq.0.and.maccr.gt.(1.05d0*dmaccrg)) then
            if (icenter.eq.0) then
               maccr1 = maccr
               goto 10
            else
               if (icrasha.eq.1.and.facc(iaccbot).eq.fmin0) then
                  icrasha = 0
                  facc(iaccbot) = facc(iaccbot)/fdtacc
                  fmin0 = facc(iaccbot)
               else
                  fmax = facc(iaccbot)
                  facc(iaccbot) = facc(iaccbot)-(fmax-fmin)*0.5d0
               endif
               goto 20
            endif
         endif
      enddo

      if (maccr1.eq.0.d0.and.icenter.eq.0) then
         icenter = 1
         iaccbot = 2
         icrasha = 0
         fmin = facc(iaccbot)
         facc(iaccbot) = fmin0
         goto 20
      endif

      if (maccr.lt.(0.95d0*dmaccrg).and.icenter.eq.1) then
         fmin = facc(iaccbot)
         facc(iaccbot) = facc(iaccbot)+(fmax-fmin)*0.5d0
         goto 20
      endif

      time = time-dtn
      dtn = maccr*sec/(massrate*msun)
      dmaccr = maccr/msun
      time = time+dtn
      raccbot = r(iaccbot)/r(nmod)
      maccbot = m(iaccbot)/m(nmod)
      lacc = alphaccr*g*totmv*thacc/(r(nmod)*dtn)
      dlim = dtacc0/dtn
      macc(nmod) = macc(nmod1)

      if (itime.eq.0.and.icrasha.eq.1.and.icenter.eq.0.and.dlim.gt.
     &     facdt) then
         itime = 1
         iaccbot = iaccbot-2
         goto 10
      endif

***   end accretion : define accretion variables

 25   do i = nmod1,1,-1
         if (i.gt.iaccbot.and.i.lt.iacctop) then
            rmid = (r(i)-r(i-1))*0.5d0
            dfaccdr(i) = dfdr (i-1,rmid,facc(i),schwar(i),idet)
            if (i.eq.iaccbot+1.and.idet.eq.0) dfaccdr(i) = 0.d0
            tacc(i) = (macc(nmod)-macc(i))/dmaccr
            if (accphase.eq.0.and.tacc(i).gt.0.d0) then
               vacc(i) = alphac*hp(i)/tacc(i)
            else
               vacc(i) = 0.d0
            endif
         else
            if (icenter.eq.0) facc(i) = 0.d0
            dfaccdr(i) = 0.d0
            tacc(i) = macc(nmod)/dmaccr
            vacc(i) = 1.d-37
         endif
      enddo
      tacc(nmod1) = tacc(nmod1-1)
      tacc(nmod) = tacc(nmod1)
      vacc(nmod1) = vacc(nmod1-1)
      vacc(nmod) = vacc(nmod1)
      dfaccdr(nmod) = dfaccdr(nmod1)
      if (accphase.eq.0.or.accphase.ge.5) then
         write (nout,400) iaccbot+1,iacctop,facc(max(iacctop-1,1))*
     &        sinthac,t(max(iaccbot-1,1)),m(max(iaccbot-1,1))/msun
      else
         write (nout,500) iaccbot+1,iacctop,facc(iacctop)*sinthac,
     &        ibotHe,itopHe,ibotH,itopH,ibotC,itopC,shlim0
      endif

      return

***   no accretion, initializations

 30   massrate = 0.d0
      dmaccr = 0.d0
      maccbot = 1.d0
      raccbot = 1.d0
      lacc = 0.d0
      facc(1:nmod) = 0.d0
      vfacc(1:nmod) = 0.d0
      macc(1:nmod) = 0.d0
      dfaccdr(1:nmod) = 0.d0
      tacc(1:nmod) = 1.d37
      vacc(1:nmod) = 1.d-37


 100  format ('** CHANGE PARAMETER CARD : menv > m(nmod)')
 200  format ('** CHANGE PARAMETER CARD : menv > 1.d0')
 300  format (5x,'non-convergence of f [',i3,']  -->  dtn = ',1pe11.5,
     &     ' yr, Macc =',1pe11.5,/)
 400  format (' accretion : ',i4,' - ',i4,' ; fmax = ',1pe8.2,
     &     ', Tbot = ',1pe8.2,', Mbot = ',0pf7.4,/)
 500  format (' accretion : ',i4,' - ',i4,' ; fmax = ',1pe8.2,
     &     ' | He shell : ',i4,' - ',i4,' | H shell : ',i4,
     &     ' - ',i4,' | envelope ',i4,' - ',i4,' | limit : ',0pf5.0,/)

      return
      end
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
      common /atmospheres/ rtau(127,12,59,6),rT(127,12,59,6),
     &     rhotab(127,12,59,6),Fconvtab(127,12,59,6),Ttabeff(59),
     &     gtabeff(12),Ztabeff(6),filesize,taumin,taumax


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


! 20      qtaus = max(pw23,tau0)
 20      qtaus = qtau(nmod)

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
         call select_table(tau,gs,taulim
     &        ,nmod,ts,Tatm,rhoatm,convatm)
         k=nmod
         do while (tau(k+10)<=taulim)
            qtau(k) = pw43*(Tatm(k)/ts)**4.d0
     &           -tau(k)
            k=k-1
         enddo
         valmax = 0.d0
         
         dqtautemp = 0.d0
         ddqtautemp = 0.d0
         do i=1,k
            Tatms(i)=t(i)
            taumod(i)=tau(i)
         enddo
         do i=k+1,nmod
            Tatms(i)=Tatm(i)
            taumod(i)=tau(i)
         enddo
         taue = 0.d0
         
         call splineatmb (taumod,Tatms,nmod,1.d50,1.d50,ddqtautemp,
     &        dqtautemp)
         call splintatmb (taumod,Tatms,ddqtautemp,nmod,pw23,taue)
         
         
         call splineatmb (Tatms,taumod,nmod,1.d50,1.d50,ddqtautemp,
     &        dqtautemp)
         
         call splintatmb (Tatms,taumod,ddqtautemp,nmod,ts,taueffb)
         call splineatm (tau,qtau,nmod,1.d50,1.d50,ddqtau,dqdtau)
         
         qtaus = qtau(nmod)
         qtau0 = 1.d0/dsqrt(3.d0)
         
         qtausv(1:nmod) = qtau(1:nmod)
         ddqtausv(1:nmod) = ddqtau(1:nmod)
         dqdtausv(1:nmod) = dqdtau(1:nmod)
         rhoatmsv(1:nmod) = rhoatm(1:nmod)
         
         
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


         k=nmod
         do while (tau(k).lt.1.d3)
            k = k-1
         enddo

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
      SUBROUTINE BOREAS (Mstar, Rstar, Lstar, Protday,grav,FeH,
     $     Mdot,Teff,icrazy)

************************************************************************
***   Fortran adaptation of the IDL procedure boreas.pro from S. Cranmer
***
***    boreas.pro:  Stand-alone IDL procedure to compute the mass loss rate of
***   a cool, late-type star using the models of Cranmer & Saar (2011).
***
***   Inputs (all scalars):
***
***   Mstar:  mass of star (in units of solar mass)
***   Rstar:  radius of star (in units of solar radius)
***   Lstar:  bolometric luminosity of star (in units of solar luminosity)
***   Protday:  rotation period of star (in days)
***   FeH:  metallicity (log of iron/hydrogen abundance ratio, in solar units)
***
*** Outputs: (all scalars):
***
***   Mdot:  mass loss rate (in units of solar masses per year)

***   THESE OUTPUTS BELOW ARE NOT USED IN STAREVOL
***
***   Rossby:  dimensionless Rossby number for this star
***   fstar:  dimensionless filling factor of open flux tubes at photosphere
***   Mdot_hot:  mass loss rate from just hot coronal (gas pressure) model
***   Mdot_cold:  mass loss rate from just cold (wave pressure) model
***   MATR:  hot wind's Alfvenic Mach number at the transition region
***
*** Version 1.0 : A. Palacios March 2012
***
***   This mass loss predictions are valid in the follwing domain:
***   2500 K <= Teff <= 6000 K
***   -2 dex <= log g <= 6 dex
***   0.05 <= L*/Lsun <= 8000

************************************************************************

      implicit none

      integer niter,icrazy,i,i1,iter,iok
      integer irotbin,model,modeli,model_old

      double precision Gconst,xMsun,xRsun,xLsun,boltzK,xmHyd,stefan
     $     ,xMdotyr
      double precision alphaturb,heightfac,theta,ffloor,ellperpSUN
      double precision xMcgs,xRcgs,xLcgs,Teff,grav,logg,ZoZsun,Vesc
     $     ,Vesc2,TeffSUN, ProtdaySUN, gravSUN
      double precision xteff
      double precision tauc, taucSUN, gratio,Rossby,RossbySUN,Ronorm
     $     ,fstar,pi,Mstar,Rstar,Lstar,FeH
      double precision rhophoto
      double precision C0fit(18),C1fit(18),Teffgrid(18),C0now
     $     ,C1now

      double precision xmuavg, gadia, Pphoto,Bequi,Bphoto,cs2,Hphoto
      double precision xmuavgSUN, csSUN, HphotoSUN
      double precision alpha_MU02, T0_MU02,F0_MU02
      double precision FluxAphoto,Valfphoto,vperpPHOT,ellperpphoto
     $     ,Qphoto
      double precision fTR,QTR,BTR,TTR,uinf,Lammax,LammaxSUN
      double precision cs2SUN,ReflCoef,quench,rhoTR,PTR,ValfTR
      double precision ReflCoefnew,sqbrac
      double precision FluxTR0,FluxA_TR,velfacTR,FcondTR,fraccond,fccmax
     $     ,FluxTR,AreaTR,Mdotcgs_hot,Mdot_hot,uTR_hot,MATR
      double precision enn, beta, eee36,rcrit,ucrit,Bcrit,Areaphoto
     $     ,Areacrit,action_photo,vperpcrit,rhocrit,Valfcrit,MAcrit
     $     ,macfac
      double precision Protday
      double precision Mdotcgs_cold,Mdot_cold,Mdot
      double precision fmin,fmax

      data (C0fit(i),i=1,18) / -9.0329342d0,-8.8057005d0,-8.6187019d0,
     $     -8.5343736d0,-8.5792239d0,-8.6915425d0,-8.8033701d0,
     $     -8.8791533d0, -8.9288180d0,-8.9604793d0,-8.9954977d0,
     $     -9.0593624d0, -9.1566287d0, -9.2743908d0,-9.4120161d0,
     $     -9.5781877d0, -9.7290674d0, -9.8972636d0 /

      data (C1fit(i),i=1,18) / 0.78081830d0,0.68284416d0,0.60471646d0
     $     ,0.56629124d0,0.57113263d0,0.59584083d0,0.61352883d0
     $     ,0.61218030d0,0.59646729d0,0.57949132d0,0.56963417d0
     $     ,0.57219919d0,0.58595944d0,0.60671743d0,0.63103575d0
     $     ,0.65574884d0,0.67753323d0,0.69808401d0 /

      common /outp/ irotbin,model,modeli,model_old
      common /boreas_Mloss/ Bcrit,fstar,Bequi

      Gconst  = 6.6732d-08
      xMsun   = 1.989d+33
      xRsun   = 6.96d+10
      xLsun   = 3.826d+33
      boltzK  = 1.380622d-16
      xmHyd   = 1.67333d-24
      stefan  = 5.66961d-5
      xMdotyr = 6.30276d+25
      pi =  3.1415926535d0

*** Some of the key parameters for the models are defined here

      alphaturb  = 0.5d0
      heightfac  = 0.5d0
      theta      = 0.333d0
      ffloor     = 1.0d-4
      ellperpSUN = 3.0d7
      niter      = 50

***  Set up basic stellar parameters in cgs units

      xMcgs  = Mstar * xMsun
      xRcgs  = Rstar * xRsun
      xLcgs  = Lstar * xLsun
!      Teff   = (xLcgs/(4.d0*pi*xRcgs*xRcgs*stefan))**0.25d0
      grav   = Gconst*xMcgs/(xRcgs*xRcgs)
      logg   = log10(grav)
      ZoZsun = 10.d0**(FeH)

c      Vesc2  = 2.d0*grav*xRcgs
      Vesc2  = 2.d0*Gconst*xMcgs/xRcgs
      Vesc   = sqrt(Vesc2)

      TeffSUN    = 5770.2d0
      ProtdaySUN = 25.3d0
      gravSUN    = Gconst*xMsun/(xRsun*xRsun)

***    Make sure parameters aren't too crazy

      icrazy = 0
      if (Teff.lt.1.5d3.or.Teff.gt.1.2d4)  icrazy = 1
c      if (Teff.lt.2.5d3.or.Teff.gt.6.d3)  icrazy = 1
      if (logg.lt.-4.d0.or.logg.gt.7.d0) icrazy = 1
      if (Lstar.lt.1.d-4.or.Lstar.gt.1.d6)  icrazy = 1
c      if (Mstar.lt.0.001d0.or.Mstar.gt.1.1d0)  icrazy = 1
      if (FeH.lt.-5.d0.or.FeH.gt.2.d0)  icrazy = 1

      if (icrazy.ne.0) then
         write(*,*)' Input parameters seem to be out of bounds! '
         return
      endif

***   Estimate Rossby number and open-flux filling factor

      tauc    = 314.241d0*exp(-Teff/1952.5d0)*exp(-(Teff/6250.d0)**18)
     $     +0.002d0
      taucSUN = 314.241d0*exp(-TeffSUN/1952.5d0)*exp(-(TeffSUN
     $     /6250.d0)**18)+0.002d0

      gratio = gravSUN/grav
      if (gratio.gt.1.d0) tauc = tauc * (gratio**0.18)

      Rossby    = Protday / tauc
      RossbySUN = ProtdaySUN / taucSUN
      Ronorm    = Rossby/RossbySUN
c..   Modification as in Gallet & Bouvier 2013 (F. Gallet private communication)
      fmin = 0.5d0 / (1.d0 + (Ronorm/0.16d0)**2.6d0)**1.3d0
      fmax =  1.d0 / (1.d0 + (Ronorm/0.31d0)**2.5d0)
c..   Original expression from G&B 2013
cc      fstar = 0.55d0 /(1.d0 + (Ronorm/0.16d0)**2.3d0)**1.22d0
c..   original expression in the Boreas routine
c..      fstar     = 0.5d0 / (1.d0 + (Ronorm/0.16d0)**2.6d0)**1.3d0
c..   fstar modified for STAREVOL
      fstar     = 0.4d0 / (1.d0 + (Ronorm/0.16d0)**2.1d0)**1.22d0

      if (fstar.lt.ffloor)  fstar = ffloor


***   Compute photospheric mass density
***   linera interpolation for C0 and C1

      Teffgrid(18) = 6000.d0
      do i = 1,17
         Teffgrid(i) = 2.5d3 + (6.d3 - 2.5d3)*i1/17.d0
      enddo
      do i = 1,17
         i1 = i-1
         if (Teff.ge.Teffgrid(i).and.Teff.le.Teffgrid(i+1)) then
            C0now = C0fit(i)+(C0fit(i+1)-C0fit(i))*(Teff-Teffgrid(i))
     $           /(Teffgrid(i+1)-Teffgrid(i))
            C1now = C1fit(i)+(C1fit(i+1)-C1fit(i))*(Teff-Teffgrid(i))
     $           /(Teffgrid(i+1)-Teffgrid(i))
         endif
      enddo
      if (Teff.gt.Teffgrid(18)) then
           C0now = C0fit(17)+(C0fit(18)-C0fit(17))*(Teff-Teffgrid(17))
     $           /(Teffgrid(18)-Teffgrid(17))
            C1now = C1fit(17)+(C1fit(18)-C1fit(17))*(Teff-Teffgrid(17))
     $           /(Teffgrid(18)-Teffgrid(17))
      endif


**** Fits given by S. Cranmer to F. Gallet
       xteff = Teff / 1000.
       C0now=4.0872049-48.979961*xteff+47.135345*(xteff**2.)-
     *   20.204336*(xteff**3.)+4.4110828*(xteff**4.)-
     *   0.48112223*(xteff**5.)+0.020825121*(xteff**6.)
       C1now=24.562656-25.713078*xteff+9.8731818*(xteff**2.)-
     *   1.3825986*(xteff**3.)-0.055190235*(xteff**4.)
     *   + 0.031616983*(xteff**5.)-0.0021585992*(xteff**6.)
****

      rhophoto = 10.d0**(C0now + C1now*logg)

***   Compute the photospheric equation of state (from fit to OPAL models),
***   magnetic field, and scale height

      xmuavg = 1.75d0 + 0.5d0*tanh((3.5d3-Teff)/6.d2)
      gadia  = 5.d0/3.d0

      Pphoto = rhophoto*boltzK*Teff/(xmuavg*xmHyd)
      Bequi  = sqrt(8.*pi*Pphoto)
      Bphoto = 1.13d0*Bequi

      cs2    = gadia*boltzK*Teff/(xmuavg*xmHyd)
      Hphoto = cs2 / (gadia*grav)

      xmuavgSUN = 1.75d0 + 0.5d0*tanh((3.5d3-TeffSUN)/6.d2)
      cs2SUN    = gadia*boltzK*TeffSUN/(xmuavgSUN*xmHyd)
      HphotoSUN = cs2SUN / (gadia*gravSUN)

***   Estimate surface flux of Alfven waves using a parameterized fit to the
***   Musielak et al. (2002a) kink-mode flux models

      alpha_MU02 = 6.774d0 + 0.5057d0*logg
      T0_MU02    = 5624.d0 + 600.2d0*logg
      F0_MU02    = exp(22.468d0 - 0.0871d0*logg)

c      FluxAphoto = F0_MU02 * ((Teff/T0_MU02)**alpha_MU02) * exp(-(Teff
c     $     /T0_MU02)**25)
c.. Modification according to Gallet & Bouvier 2013
      FluxAphoto = F0_MU02 * ((Teff/T0_MU02)**alpha_MU02) * exp(-(Teff
     $     /T0_MU02)**25)/3.2d0
      if (FluxAphoto.lt.1.d-10) FluxAphoto = 1.d-10

***   Set up MHD turbulence parameters at the photosphere

      Valfphoto    = Bphoto / sqrt(4.d0*pi*rhophoto)
      vperpPHOT    = sqrt(FluxAphoto/rhophoto/Valfphoto)
      ellperpphoto = ellperpSUN * (Hphoto/HphotoSUN)

      Qphoto       = alphaturb*rhophoto*(vperpPHOT**3)/ellperpphoto

***   Extrapolate MHD turbulence parameters up to the transition region (TR)

      fTR  = fstar**theta
      BTR  = Bphoto * (fstar/fTR)
      uinf = Vesc

      TTR       = 2.0d5
      Lammax    = 7.4d-23 + 4.2d-22*(ZoZsun**1.13)
      LammaxSUN = 7.4d-23 + 4.2d-22

***   Iterate to find self-consistent solution for density and heating rate
***   at the TR, assuming that the non-WKB reflection at the TR is imperfect.
***   The reflection coefficient is given by an approximate version of the
***   low-frequency Cranmer (2010) limit.

      ReflCoef = 0.5d0

      do iter=1,niter
         quench  = ReflCoef*(1.d0+ReflCoef)/(1.d0+ReflCoef**2)**1.5d0 *
     $        sqrt(2.d0)
         sqbrac  = quench*Qphoto*xmHyd*xmHyd/(rhophoto**0.25)/Lammax
         rhoTR   = (sqbrac**(4.d0/7.d0)) * (fstar**(2.d0*(1.d0-theta)
     $        /7.d0))
         QTR     = Qphoto*quench*((rhoTR/rhophoto)**0.25)*sqrt(BTR
     $        /Bphoto)
         ValfTR  = BTR / sqrt(4.d0*pi*rhoTR)
         PTR     = 2.d0*rhoTR*boltzk*TTR/xmHyd
         ReflCoefnew = abs((ValfTR-uinf)/(ValfTR+uinf))
         ReflCoef    = sqrt(ReflCoefnew*ReflCoef)
      enddo

***   Does the heating-related energy flux at the TR exceed the flux in
***   "passive propagation" of the Alfven waves?  If so, cap it!

      FluxTR0  = heightfac*QTR*xRcgs
      FluxA_TR = FluxAphoto * fstar/fTR

      if (FluxTR0.gt.FluxA_TR) FluxTR0 = FluxA_TR

***   Estimate the mass loss rate for a hot coronal wind, using the
***   Hansteen et al. (1995) energy balance approximation.

      velfacTR = 1.4d6 * sqrt(Lammax/LammaxSUN)
      FcondTR  = PTR * velfacTR
      fraccond = FcondTR / FluxTR0
      fccmax   = 0.9d0
      if (fraccond.gt.fccmax) fraccond = fccmax
      FluxTR   = FluxTR0 * (1.d0-fraccond)

      AreaTR      = fTR * (4.d0*pi*xRcgs*xRcgs)
      Mdotcgs_hot = AreaTR*FluxTR/ (0.5d0*(Vesc2+uinf*uinf))
      Mdot_hot    = Mdotcgs_hot / xMdotyr

      uTR_hot = Mdotcgs_hot / (rhoTR*AreaTR)
      MATR    = uTR_hot / ValfTR

***   For the Holzer et al. (1983) cold wave-driven wind, first estimate the
***   radius, speed, and magnetic field at the wave-modified critical point.

      enn   = 2.0d0       ! B \propto r^{-enn} at crit point
      beta  = 0.5d0*enn
      eee36 = (3.d0/(7.d0*beta)+4.d0/7.d0)/(1.d0+(vperpPHOT/Vesc)**2)

      rcrit = 1.75*eee36*xRcgs
      ucrit = sqrt(Gconst*xMcgs/(enn*rcrit))

      Bcrit = Bphoto * ((xRcgs/rcrit)**2) * fstar
      print *,'Bcrit,Bphoto',Bcrit,Bphoto

***   At the critical point, iterate on the definition of u_crit and the
***   condition of wave action conservation to get density and wave amplitude.

      Areaphoto    = fstar * (4.d0*pi*xRcgs*xRcgs)
      Areacrit     = 4.d0*pi*rcrit*rcrit
      action_photo = rhophoto*(vperpPHOT**2)*Valfphoto*Areaphoto

      vperpcrit = 2.d0*ucrit
      rhocrit   = 4.d0*pi*
     $     (action_photo/((vperpcrit**2)*Bcrit*Areacrit))**2

      do iter=1,niter
         Valfcrit  = Bcrit / sqrt(4.d0*pi*rhocrit)
         MAcrit    = ucrit / Valfcrit
         macfac    = (1.d0 + 3.d0*MAcrit)/(1.d0 + MAcrit)
         vperpcrit = 2.d0*ucrit / sqrt(macfac)
         rhocrit   = action_photo /((vperpcrit**2)*Valfcrit*Areacrit*
     $        (1.d0+MAcrit)**2)
      enddo

      Mdotcgs_cold = rhocrit*ucrit*Areacrit
      Mdot_cold    = Mdotcgs_cold / xMdotyr

***   Estimate the actual mass loss from both hot and cold processes.

      Mdot = Mdot_cold + (Mdot_hot*exp(-4.d0*MATR**2))

      return
      end
      SUBROUTINE calcevo (*)

************************************************************************
*     Compute the next model and, in case of success, save the results     *
*     If convergence fails, the procedure is repeated by changing the      *
*     current time-step                                                    *
*     09/03: Simplification of the procedure of re-reading initial model   *
*     in case of crashes ( see also modification of rinimod)        *
*     Modifs CC ondes (2/04/07 --> 23/11/07)                               *
*     *
*     $LastChangedDate:: 2016-05-13 11:10:12 +0200 (Ven, 13 mai 2016)    $ *
*     $Author:: palacios                                                 $ *
*     $Rev:: 73                                                          $ *
*     *
************************************************************************
#ifdef GYRE
c........................Rajoute par Alice
      use evolstell_gyre
#endif
      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.atm'
      include 'evolcom.cons'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.ion'
      include 'evolcom.lum'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.opa'
      include 'evolcom.rot'
      include 'evolcom.conv'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.var'

#ifdef GYRE     
c......................Rajoute par Alice
      include 'evolcom.grad'
      include 'omp_lib.h'
#endif      
      integer imerk,icrasht,iflag
      integer ntprof0,itermax0
      integer neff1, nconvsh, ienvsh
      integer error,iexp
      integer i,j,k,l,kl,ksh
      integer klenv,klpulse,klcore,iHecore
      integer imin,imax,itmaxf,itmax,ishockb,ishockt
      integer fcs               ! first convective shell

      double precision alpha01,alpha00,alphmin0
      double precision ledd,lnuc,vlnuc,dabslum,dloglum,dfnuc
      double precision sum1,taueffl,taune1l,taunel,rne1l,rnel,tne1l,
     &     tnel,rone1l,ronel,lne1l,lnel,deltnel,deltes,dtna,tmaxf,tmax,
     &     Lmax,Mackmax,enucmax,renucl,vi
      double precision zzmean,amean,xamean,yeq,y0,x,y,irrad
      double precision eng,engdt,engdro,engnup,engnupdt,engnupdro,engro
     &     ,engrocap
      double precision elemsave,mconv
      double precision menvconv,dmenvconv,mcore
      double precision omegacritic,Rpolinv
      double precision mHecore
      double precision geffrot,Rp2
      double precision times,chronos ! Ajout chronos pour mesure du tps de calcul

      logical partialmix,pass
      logical vlconv,adapt

      character month*36,GetDateTimeStr*15,fmt*26
      parameter (month = 'JanFebMarAprMayJunJulAugSepOctNovDec')
      parameter (fmt = '(I2.2,A1,I2.2,I3,A3,I4)')
      integer itime(8)

c................Rajoute par Alice
      logical :: res, resdet
      integer, dimension(nsh) :: normx
      logical, dimension(nsh) :: mask
      double precision brunt
      double precision bruntv(nsh), bruntv2(nsh)      
      common /brunt/ brunt(nsh),bruntv2,bruntv
      double precision xmHb, xmHt
      common /hburn/ xmHb, xmHt
      logical rgbphase_gyre
c.......................Fin rajout
      
      common /saveMS/ elemsave
      common /convergence/ imerk,icrasht
      common /grid/ alpha00,alpha01,alphmin0,ntprof0,itermax0
      common /overshoot/ klcore,klenv,klpulse
      common /resolution/ vlconv
      common /hydrodyn/ Lmax,tmax,Mackmax,enucmax,itmax,ishockb,ishockt
      common /meshconv/ menvconv,dmenvconv,mcore,nconvsh,ienvsh
      common /geffective/ geffrot(nsh),Rp2(nsh)
c     common /eddington/ ledd

      dimension renucl(nsh),vi(nreac),y0(nsp),x(nsp),y(nsp)


      
      error = 0
c       if (nphase.eq.4.and.novlim(1,3).eq.1) then
c        iHecore = novlim(1,4)
c        mHecore = m(iHecore)
c       else
c        iHecore = 1
c        mHecore = 0.d0
c       endif
      model = model+1
      maxsh = abs(maxsh0)
      if (model.le.0) time = 0.d0
      if (no.eq.0) then
         itermax0 = itermax
         alpha00 = alpha0
         alpha01 = alpha
         ntprof0 = ntprof
         alphmax0 = alphmax
         alphmin0 = alphmin
      endif
      if (model.le.3.and.no.eq.2) then
         itermax = itermax0
         alphmax = alphmax0
         alphmin = alphmin0
         alpha = alpha01
         alpha0 = alpha00
         ntprof = ntprof0
      endif
      
      if (imodpr.eq.11.and.ireset.eq.0) call mycompo

*------------------------------
***   some specific adaptations
*------------------------------

      if (icrasht.ne.0) then
         icrash = 9
         ireverse = 1
         write (nout,300)
         write (90,300)
      else
         icrash = abs(icrash0)
         if (icrash0.lt.0) ireverse = -1
      endif

c..   before AGB phase
      if (nphase.le.4.and..not.rgbphase.and..not.igwrot)
     &     maxmod = min(maxmod,300) ! modif Fev 2019

c..   AGB phase

c..   increase structural time step & constrain the number of saved models
c..   during the pulse
      if (ifail.eq.0) then
         if (thermalpulse.and.imodpr.ne.0) then
            fts = min(fts0,2.d-2)
c..   pulse starts, save immediately
            if (maxmod.eq.maxmod0) then
               if (model-modeli.gt.10) then
                  maxmod = no+1
               else
                  maxmod = min(11-mod(model-modeli,10)+no,maxmod0)
               endif
               write (nout,600) fts,maxmod
            else
               write (nout,650) fts
            endif
         else
            fts = fts0
         endif
c..   further constrain structural time step during 3DUP
         if ((icorr.eq.'t'.or.icorr.eq.'a').and.dup3) then
c     if (dup3) then
            fts = min(fts0,2.d-2)
            if (imodpr.eq.33.and.maxmod.eq.maxmod0) then
               maxmod = min(41-mod(model-modeli,40)+no,maxmod0)
               write (nout,700) fts,maxmod
            else
               write (nout,750) fts
            endif
         else
            fts = fts0
         endif
      endif
c..   thibaut
c     if (diffherwig.and.agbphase.and.totm.lt.10.1d0.and.dup3.and.
c     &     dtn*seci.gt.0.1d0) then
c     dtn = 0.1d0*sec
c     write (nout,108) dtn*seci
c     write (90,108) dtn*seci
c     endif

      ifail = 0
      no = no+1
      
*______________________________
***   abundance renormalization
*------------------------------

c..   index of maximum for each row
      normx(:) = maxloc(xsp, 2)
c..   loop over rows of xsp
      do k = 1,nmod
c..   compute sum the row      
         sum1 = sum(xsp(k, 1:nsp))
c..   for the element of the row which is the maximum,
c     remove 1+(sum of the row) so that sum of the row = 1 now         
         xsp(k, normx(k)) = maxval(xsp(k,1:nsp))+(1.d0-sum1)
      enddo

*__________________________________________________________________
***   extrapolation of the next model and variable initializations
*------------------------------------------------------------------

      call guess
      
*______________________________________________________________________
***   calculation of the factor ft and fp for the treatment of rotation
*----------------------------------------------------------------------

      if (hydrorot) then
         call potential
      else
         do k = 1,nmod
            ft(k) = 1.d0
            fp(k) = 1.d0
         enddo
      endif
      
*___________________________________________________
***   calculation of the stellar structure evolution
*---------------------------------------------------

      if (ireset.eq.0) alpha = alpha0
      
      call nwtraf (*20,error)
      
*________________________
***   test He-core growth
*------------------------
c     if (nphase.eq.4.and.xsp(novlim(1,4),ihe4).gt.
c     &     vxsp(novlim(1,4),ihe4)) then
c     write(nout,*)'Breathing pulse developing : central He4',
c     &        ' abundance increased during this timestep'
c     error = 35
c     goto 20
c     endif

c     if (nphase.eq.4.and.mHecore.gt.0.d0.and.xsp(1,ihe4).lt.0.5d0.and.
c     &     xsp(1,ihe4).gt.0.03d0) then
c     if ((mHecore-m(novlim(1,4))).gt.0.1d0*msun.and.
c     &        (novlim(2,3)-novlim(1,4)).gt.5) then
c     c         if (abs(m(novlim(1,4))-mHecore).gt.0.1d0*msun) then
c     write (nout,75) m(novlim(1,4))/msun,novlim(1,4),
c     &           mHecore/msun,iHecore
c     if (xsp(1,ihe4)/vxsp(1,ihe4).gt.0.001d0) then
c     write (nout,76) vxsp(1,ihe4),novlim(1,3),xsp(1,ihe4)
c     $           ,novlim(1,3)
c     error = 35
c     goto 20
c     endif
c     endif

*______________________________
***   test temperature increase
*------------------------------

      if (ftst.gt.0.d0.and.model.gt.10) then
         mask(:)=m.gt.mcut     
c..   Apply mask for t corresponding to m(i)<mcut except for t(1)
c     with merge function   
c..   Find tmax for elements if m(index_element)<mcut
         tmaxf=maxval([t(1),merge(t(2:), 0.d0, mask(2:))])
c..   Find index of tmax     
         itmaxf=maxloc([t(1),merge(t(2:), 0.d0, mask(2:))],1)
         if (abs(log10(tmax/tmaxf)).gt.log10(1.d0+ftst).and.max(tmax,
     &        tmaxf).gt.2.d7) then
            icrasht = icrasht+1
            if (icrasht.le.iresmax) then
               write (nout,80) int(ftst*100.d0),tmax,itmax,m(itmax)/msun
     &              ,tmaxf,itmaxf,m(itmaxf)/msun
               error = 40
               ireset = max(ireset-1,0)
               goto 20
            else
               stop 'Temperature increase too large : modify ftst'
            endif
         endif
      endif

      call thermo (error)


c     if (imodpr.eq.-1) then
c     call initondes
c     write (*,'(" mv sortie_ondes sortie_ondes_",i5)') model
c     stop
c     endif
      if (error.gt.0) goto 20



*____________________________________________________________________
***   
***   beyond this point the STRUCTURE has converged otherwise goto 20
***   
*--------------------------------------------------------------------


      icrasht = 0
      maxsh = abs(maxsh0)
      mixopt = mixopt0
c     print *,coefssv(:)
c     coefssv(:) = coefsv(:)
c     print *,coefssv(:)

*____________________________________________________________
***   calculation of the energy loss by plasma neutrinos
***   nuclear energy production rates and neutron irradiation
***   in case lnucl = .false.
*------------------------------------------------------------

      if (.not.lnucl) then
         do k = 1,nmod
            do l = 1,nis-1
               x(l) = max(ymin*anuc(l),xsp(k,l))
               y(l) = ysp(k,l)
            enddo
            x(nis) = 1.d-50
            x(nsp) = 1.d-50
            y(nis) = 1.d-50
            y(nsp) = 1.d-50
            tconv(k) = 1.d99
            if (t(k).gt.1.d7) then
               zzmean = 0.d0
               xamean = 0.d0
               do j = 1,nsp
                  zzmean = zzmean+znuc(j)*y(j)
                  xamean = xamean+y(j)
               enddo
               amean = 1.d0/xamean
               zzmean = zzmean*amean
               call nulosc (t(k),ro(k),amean,zzmean,engnup,engnupdt,
     &              engnupdro)
            else
               engnup = 0.d0
               engnupdt = 0.d0
               engnupdro = 0.d0
            endif
            enupla(k) = engnup

c..   calculation of the nuclear energy production rates
c..   neutrino energy losses and neutron irradiation

            irrad = 7.68771025912336d0*ysp(k,1)*ro(k)*dtn*dsqrt(t(k))
            if (irrad.le.1.72d0) then
               signt(k) = 43.8372d0*irrad+10.4d0
            else
               if (irrad.lt.2.92d0) then
                  signt(k) = -69.8417d0*irrad+205.93d0
               else
                  signt(k) = 1.99d0
               endif
            endif
            ksh = k
            call denucl (t(k),ro(k),mueinv(k),y,ksh,eng,engdt,engdro)
            enucl(k) = eng-enupla(k)

         enddo
      endif
      
*_______________________________________________________________
***   computation of vHconv and of the new convective boundaries
*---------------------------------------------------------------

      vhconv(1) = 0.d0
      do i = 2,nmod
         vhconv(i) = 0.d0
         tm(i) = t(i-1)**wi(i)*t(i)**wj(i)
      enddo
      tm(1) = t(1)
      rom(1) = ro(1)
      if (vlconv) call convzone
      do kl = 1,nsconv
         imin = novlim(kl,3)
         imax = novlim(kl,4)
         if (vlconv) call mlt (imin,imax,error)
         do i = imin,imax
            vhconv(i) = pim4*r(i)*fconv(i)/(tm(i)*rom(i))
         enddo
      enddo

***   compute parametric overshoot (alpha * Hp)
      if (novopt.ge.1) call oversh

*________________________________________________________________
***   calculation of the associated nucleosynthesis and diffusion
*----------------------------------------------------------------

c     if (ondes) call ondes
      if (idiffcc.and.agbphase.and..not.idiffcc0) then
         write (nout,100) idiffcc,idiffcc0
         idiffcc = idiffcc0
      endif

      if (nucreg.ne.3) then
         klcore = 0
         klpulse = 0
         klenv = 0

***   Compute mixing coefficients and nucleosynthesis
         pass = .true.
c     x         if (flame.and.diffzc) pass = .false.

c..   CASE 1 : diffusion (if nmixd > 0, diffusion treated in netdiff)
c..   
C     Modif CC ondes (19/10/07)
         if ((microdiffus.or.rotation.or.diffusconv.or.difftacho.or.
     &        thermohaline.or.igw).and.nmixd.eq.0.and.pass) then
C     &        thermohaline).and.nmixd.eq.0.and.pass) then
C     <--
!     Tests activation diffusion atomique ou non (modif Td Fev.2019)
            if (microdiffus.and.nphase.gt.5) then
               microdiffus = .false.
c     print *, 'phase after 5 microdiffus false'
            else if (nphase.eq.1.and.time*seci.lt.2e7) then ! Ajout test TD Juillet.2019 si rotation
c     else if (nphase.eq.1.and.time*seci.lt.2e6) then ! Ajout test TD AP Fev.2019
               microdiffus = .false.
c     print *,'no diffusion this young star'
            else if (nphase.eq.1.and.novlim(1,3).eq.1.and.nsconv.eq.1)
     $              then   
               microdiffus = .false.
c     print *,'microdiffus false - no radiative core'
            endif
!     Fin test activation diffusion atomique
            if (idiffnuc.eq.2) then
               if (nucreg.ne.0) call chedif (error)
               if (error.gt.0) goto 20
               call diffusion (1,error)
               if (error.gt.0) goto 20
            else
               call diffusion (1,error)
               if (error.gt.0) goto 20
               if (nucreg.ne.0) call chedif (error)
               if (error.gt.0) goto 20
            endif
            if (idiffnuc.eq.3) then
               call diffusion (1,error)
               if (error.gt.0) goto 20
            endif

         else

c..   CASE 2 : no diffusion or full coupling
c..   
            if (nucreg.ne.0) then
               if (nsconv.gt.0) then
                  if (novlim(1,3).eq.1) klcore = 1
                  do kl = 1,nsconv
                     if ((tau(novlim(kl,4)).lt.1.d2.or.t(novlim(kl,4))
     &                    .lt.1.d6).and.mr(novlim(kl,4)).gt.0.7d0.and.
     &                    klenv.eq.0) klenv = kl
                     if (t(novlim(kl,3)).gt.2.d8.and.novlim(kl,3).gt.1.
     &                    and.nphase.eq.5.and.klpulse.eq.0) klpulse = kl
                  enddo
                  if (nsconv.eq.1.and.klcore.eq.1) then
                     klenv = 0
                     klpulse = 0
                  endif
               endif
               call chedif (error)
               if (error.gt.0) goto 20
            else

c..   CASE 3 : no nucleosynthesis  (nucreg = 0)
c..   restore initial abundances and do mixing only
               do l = 1,nsp
                  do k = 1,nmod
                     xsp(k,l) = vxsp(k,l)
                     ysp(k,l) = xsp(k,l)/anuc(l)
                  enddo
               enddo
               do kl = 1,nsconv
                  imin = novlim(kl,3)
                  imax = novlim(kl,4)
                  call mix (dm,imin,imax,kl,0,partialmix)
               enddo
            endif
         endif
      else
c..   nucleosynthesis computed during convergence process
c..   simply update composition
         do l = 1,nsp
            do k = 1,nmod
               vxsp(k,l) = xsp(k,l)
            enddo
         enddo
      endif


***   check that diffusion and nucleosynthesis have a small energetic impact
***   and compute neutron equilibrium abundances

      lnuc = 0.d0
      vlnuc = 0.d0
      if (no.eq.1.and.nuclopt.eq.'n') write (90,*) '*** Compute final ',
     &     'neutron equilibrium abundance after nucleosynthesis'
      do k = 1,nmod
         do l = 1,nsp
            x(l) = xsp(k,l)
            y(l) = xsp(k,l)/anuc(l)
         enddo
         vlnuc = vlnuc+enucl(k)*dm(k)
         call vit (t(k),ro(k),mueinv(k),k,vi,0)
         if (nuclopt.eq.'n') then
            do l = 1,nsp
               y0(l) = x(l)/anuc(l)
            enddo
            call yneutron (y0,vi,yeq)
            xsp(k,1) = yeq
            vxsp(k,1) = yeq
            x(1) = yeq
         endif
         ksh = k
         call nuceng (-1,y,vi,ksh,eng,engro,engrocap)
         renucl(k) = eng-enupla(k)
         lnuc = lnuc+renucl(k)*dm(k)
      enddo
      if (nuclopt.eq.'m') then
         do kl = 1,nsconv
            imin = novlim(kl,3)
            imax = novlim(kl,4)
            yeq = 0.d0
            mconv = 0.d0
            do i = imin,imax
               mconv = mconv+dm(i)
               yeq = yeq+xsp(i,1)*dm(i)
            enddo
            yeq = yeq/mconv
            do i = imin,imax
               xsp(i,1) = yeq
            enddo
         enddo
      endif

      dloglum = abs(log(abs(vlnuc/lnuc)))
      dabslum = min(exp(dloglum)-1.d0,1.d6)
      dfnuc = log(1.d0+ftnuc)
      if (nphase.gt.1.and.model.gt.10.and.
     &     dloglum.gt.min(log(2.d0),dfnuc)) then
         if (dloglum.gt.dfnuc) then
            write (90,1810) model,min(9999,int(dabslum*1.d2)),
     &           int(ftnuc*1.d2)
            write (nout,1810) model,min(9999,int(dabslum*1.d2)),
     &           int(ftnuc*1.d2)
            error = 30
c     mixopt = .true.
         else
            error = -30
         endif
         if (dloglum.gt.1.d0.and.nmixd.eq.0.and.nphase.gt.2.and..not.
     &        thermalpulse.and.(nuclopt.eq.'c'.or.nuclopt.eq.'u')) then
            error = 30
            tolnuc = max(tolnuc0*1.d2,1.d-6)
c..   increase number of saved models in case nmixd = 5
            if (imodpr.eq.33) then
               maxmod = min(11-mod(model-modeli,10)+no,maxmod0)
            else
               maxmod = min(51-mod(model-modeli,50)+no,maxmod0)
            endif
            write (90,1850) model,nmixd0,maxmod0,maxmod
            write (nout,1850) model,nmixd0,maxmod0,maxmod
            nmixd = 5
         endif
         if (error.gt.0) goto 20
      endif
c..   save additional models if H is modified by more than 0.05 during MS
c..   at least 10 models computed
      if (xsp(1,ih1).lt.elemsave.and.nphase.eq.2) then
         maxmod = min(model-modeli,maxmod0)
         maxmod = max(maxmod,10)
         write (nout,2000)
         write (90,2000)
      endif
c..   idem for He during core He burning
      if (xsp(1,ihe4).lt.elemsave.and.nphase.eq.4) then
         maxmod = min(model-modeli,maxmod0)
         maxmod = max(maxmod,10)
         write (nout,2100)
         write (90,2100)
      endif

c..   save more models in the Hertzsprung gap
c      if (nphase.eq.3.and..not.rgbphase) maxmod = min(maxmod,25) ! couper TD 03/2020

*_______________________________________________________________________
***   determination of the photosphere location, temperature and density
*-----------------------------------------------------------------------

      neff = 0
      do k = nmod-2,2,-1
         if (tau(k).gt.taueff.and.tau(k+1).ge.taueff) exit
      enddo
      neff = k+2


      print *,'taueff=',taueff

      neff = min(neff,nmod)
      neff1 = neff-1
      taueffl = log(abs(taueff))
      taune1l = log(abs(tau(neff1)))
      taunel = log(abs(tau(neff)))
      rne1l = lnr(neff1)
      rnel = lnr(neff)
      tne1l = lnt(neff1)
      tnel = lnt(neff)
      rone1l = log(abs(ro(neff1)))
      ronel = log(abs(ro(neff)))
      lne1l = log(abs(lum(neff1)))
      lnel = log(abs(lum(neff)))
      deltnel = tne1l-tnel
      tphot = tne1l-(taune1l-taueffl)*deltnel/(taune1l-taunel)
      deltes = tne1l-tphot
      reff = rne1l-deltes*(rne1l-rnel)/deltnel
      roeff = rone1l-deltes*(rone1l-ronel)/deltnel
      roeff = min(roeff,1.d0)
      roeff = max(roeff,-20.d0)
      leff = lne1l-deltes*(lne1l-lnel)/deltnel
      leff = exp(leff)
      tphot = exp(tphot)
      reff = exp(reff)
      roeff = exp(roeff)

c..   TEST ROTATION DECEMBRE 2014 For a rotating star, the expression of
c..   the Eddington luminosity is modified to account for the change in
c..   the gravitational potential that affects the mass

      if (hydrorot) then
         ledd = pim4*c*g*m(neff)/kap(neff)*(1.d0-pw23*vomega(neff)**2
     $        *Rp2(neff)**3/(g*m(neff)))
      else
         ledd = pim4*c*g*m(neff)/kap(neff)
      endif
      print *,'Ledd',ledd, pim4*c*g*m(neff)/kap(neff)
      if (leff.gt.ledd) then
         write (nout,90) model,leff/lsun,ledd/lsun
         error = 66
         goto 20
c      nmod = neff
      endif

*______________________________________
***   calculation of surface quantities
*--------------------------------------
      if (hydrorot) then
         if (ntprof.eq.2) then
            geff = log10(geffrot(neff-1))
         else
            geff = log10(geffrot(neff))
         endif

      else
         geff = log10(g*m(neff)/(reff*reff))
      endif
      soll = leff/lsun
      solr = reff/rsun
c     teffy = tsurf
      teffy = teff
      if (time.lt.1.d0) then
         write (nout,1100) model,min(soll,9999999.d0),min(9999.d0,solr),
     &        teffy,totm,dtn,time
      else
         write (nout,1200) model,min(soll,9999999.d0),min(9999.d0,solr),
     &        teffy,totm,dtn*seci,time*seci
      endif

*__________________________________________________________
***   mixing atmosphere abundances with convective envelope
*----------------------------------------------------------

      if (nsconv.gt.0.and.tau(novlim(nsconv,4)).lt.1.d2) then
         imin = novlim(nsconv,4)
         do i = 1,nsp
            do j = imin,nmod
               xsp(j,i) = xsp(imin-1,i)
               ysp(j,i) = xsp(j,i)/anuc(i)
            enddo
         enddo
      endif


*___________________________
***   storage of the results
*---------------------------

      call cpu_time (timeend)
      cputime = min(9999.d0,(timeend-timebeg))

      call prvar

*___________________________
***   Rajoute par Alice
*---------------------------

#ifdef GYRE

      ! For first iteration
      if (model .eq. 1) then
         gyreprec = .false.
         tpgy = 0.d0
         lpgy = 0.d0
         deltat = 0.d0
      endif

      dtgy = dtn*seci
      mgy = totm
      tgy = t(nmod)
      lgy = lum(nmod)/lsun

      xspgy = xsp(1,2)

      rgbphase_gyre = rgbphase.and.xmHb/xmHt.ge.0.8d0
      
      iter_gyre = imodpr  ! by default
      n_time_gyre = model

      call evo_iter("criteria.in", res, resdet, nphase, rgbphase_gyre)
   
      if (res .or. resdet) then ! if gyre must run in summary or detail mode     

         !msun = 1.989d33
         !rsun = 6.9599d10
         !lsun = 3.86d33
         nversion_gyre = 101
         n_gyre = nmod
         mstar_gyre = totm*msun
         rstar_gyre = r(nmod)
         lstar_gyre = lum(nmod)
         age_gyre = time*seci

         if (allocated(omega_gyre)) then
            deallocate(omega_gyre)
         endif
         if (allocated(kapkt_gyre)) then
            deallocate(kapkt_gyre)
         endif
         if (allocated(kapkro_gyre)) then
            deallocate(kapkro_gyre)
         endif
         if (allocated(enuclt_gyre)) then
            deallocate(enuclt_gyre)
         endif
         if (allocated(enuclro_gyre)) then
            deallocate(enuclro_gyre)
         endif
         if (allocated(r_gyre)) then
            deallocate(r_gyre)
         endif
         if (allocated(mr_gyre)) then
            deallocate(mr_gyre)
         endif
         if (allocated(lum_gyre)) then
            deallocate(lum_gyre)
         endif
         if (allocated(rho_gyre)) then
            deallocate(rho_gyre)
         endif
         if (allocated(t_gyre)) then
            deallocate(t_gyre)
         endif
         if (allocated(p_gyre)) then
            deallocate(p_gyre)
         endif
         if (allocated(abla_gyre)) then
            deallocate(abla_gyre)
         endif
         if (allocated(bruntv2_gyre)) then
            deallocate(bruntv2_gyre)
         endif
         if (allocated(gamma1_gyre)) then
            deallocate(gamma1_gyre)
         endif
         if (allocated(deltaks_gyre)) then
            deallocate(deltaks_gyre)
         endif
         if (allocated(kap_gyre)) then
            deallocate(kap_gyre)
         endif
         if (allocated(enucl_gyre)) then
            deallocate(enucl_gyre)
         endif
         if (allocated(abad_gyre)) then
            deallocate(abad_gyre)
         endif


         allocate(omega_gyre(n_gyre))
         allocate(kapkt_gyre(n_gyre))
         allocate(kapkro_gyre(n_gyre))
         allocate(enuclt_gyre(n_gyre))
         allocate(enuclro_gyre(n_gyre))
         allocate(r_gyre(n_gyre))
         allocate(mr_gyre(n_gyre))
         allocate(lum_gyre(n_gyre))
         allocate(rho_gyre(n_gyre))
         allocate(t_gyre(n_gyre))
         allocate(p_gyre(n_gyre))
         allocate(abla_gyre(n_gyre))
         allocate(bruntv2_gyre(n_gyre))
         allocate(gamma1_gyre(n_gyre))
         allocate(deltaks_gyre(n_gyre))
         allocate(kap_gyre(n_gyre))
         allocate(enucl_gyre(n_gyre))
         allocate(abad_gyre(n_gyre))

         
         do i = 1,nmod
            r_gyre(i) = r(i)
            mr_gyre(i) = mr(i)*msun
            lum_gyre(i) = lum(i)
            p_gyre(i) = p(i)
            t_gyre(i) = t(i)
            rho_gyre(i) = ro(i)
            abla_gyre(i) = abla(i)
            bruntv2_gyre(i) = bruntv2(i)
            gamma1_gyre(i) = gamma1(i)
            abad_gyre(i) = abad(i)
            deltaks_gyre(i) = deltaks(i)
            kap_gyre(i) = kap(i)
            kapkt_gyre(i) = t(i)*dkapdt(i)
            kapkro_gyre(i) = ro(i)*dkapdro(i)
            enucl_gyre(i) = enucl(i)
            enuclt_gyre(i) = t(i)*denucldt(i)
            enuclro_gyre(i) = ro(i)*denucldro(i)
            omega_gyre(i) = omega(i)
         enddo

!     open(unit=32, file="verif.log", position=
!     & "append", status ="unknown")
         
!     write(32, *)
!     write(32, *) "Time iteration :", n_time_gyre
!     write(32, 5100) "mr_gyre", "lum_gyre",
!     & "bruntv2_gyre", "enucl_gyre", "omega_gyre",
!     & "dkapdt_gyre", 
!     & "dkapdro_gyre", "denucldt_gyre", "denucldro_gyre", 
!     & "deltaks_gyre"
!     write_loop : do j = 1, size(r_gyre)
!     write(32, 6000) mr_gyre(j), lum_gyre(j),
!     & bruntv2_gyre(j), enucl_gyre(j), omega_gyre(j), 
!     & kapkt_gyre(j), 
!     &   kapkro_gyre(j), enuclt_gyre(j),         
!     &   enuclro_gyre(j), deltaks_gyre(j)
!     end do write_loop

!     close(unit=32) 

!     5100     format(10(A25))
!     6000     format(10(E25.16E3))

         call evolstell_build_data
         call evolstell_run('evolgyre','detail', resdet)

         tpgy = tgy
         lpgy = lgy

      endif
#endif      
c....................................Fin rajout
      
      call resulpr (error)
      print *, 'error',error
      if (error.eq.-9) then
         write (90,*) '   time-dependent convection activated'
         write (nout,*) '   time-dependent convection activated'
      endif
      
      iacc = iacc0
      hydro = hydro0
      ifail = 0
      if (error.eq.-30) then
         write (90,1800) model,min(9999,int(dabslum*1.d2))
         write (nout,1800) model,min(9999,int(dabslum*1.d2))
         if (icorr.eq.'t') ifail = -1
      endif
      write (90,1900)

      if (no.lt.maxmod) then
         call flush (90)
         return 1
      endif

      call date_and_time (values=itime)
      write (GetDateTimeStr,FMT) itime(5),':',itime(6),itime(3),
     &     month(3*itime(2)-2:3*itime(2)),itime(1)


      write (90,*)
      write (90,*)
      write (90,'("stop normal : ",a15)') GetDateTimeStr
      write (90,*)
      write (nout,*)
      write (nout,*) "stop 'normal'"


      return



*____________________________________
***   Managing crash situations
*------------------------------------

 20   if (error.eq.18) then
         write (90,1400)
         write (nout,1400)
      else
         if (error.eq.66) then
            write (90,1666) model
            write (nout,1666) model
         endif
         if (error.eq.50) then
            write (90,1555) model
            write (nout,1555) model
         endif
         if (error.eq.99) then
            write (90,1300) model
            write (nout,1300) model
         endif
         if (error.eq.1) then
            write (90,3001)
            write (nout,3001)
         endif
         if (error.eq.2) then
            write (90,3002)
            write (nout,3002)
         endif
         if (error.eq.5) then
            write (90,3005)
            write (nout,3005)
         endif
         if (error.eq.6) then
            write (90,3006)
            write (nout,3006)
         endif
c..   OPACITY problem
         if (error.eq.7) then
            write (90,3007)
            write (nout,3007)
         endif
         if (error.eq.80) then
            write (90,3080)
            write (nout,3080)
         endif
c     if (error.eq.81) then
c     write (90,3081)
c     write (nout,3081)
c     endif
         if (error.eq.82) then
            write (90,3082)
            write (nout,3082)
         endif
         if (error.eq.83) then
            write (90,3083)
            write (nout,3083)
         endif
         if (error.eq.84) then
            write (90,3084)
            write (nout,3084)
         endif
c..   EOS problem
         if (error.eq.10) then
            write (90,3010)
            write (nout,3010)
         endif
         if (error.eq.11) then
            write (90,3011)
            write (nout,3011)
         endif
         if (error.eq.12) then
            write (90,3012)
            write (nout,3012)
         endif
c..   DIFFUSION problem
         if (error.eq.20) then
            write (90,3020)
            write (nout,3020)
         endif
         if (error.eq.21) then
            write (90,3021)
            write (nout,3021)
         endif
         if (error.eq.22) then
            write (90,3022)
            write (nout,3022)
         endif
         if (error.eq.23) then
            write (90,3023)
            write (nout,3023)
         endif
         if (error.eq.15) then
            write (90,3015)
            write (nout,3015)
         endif
         if (error.eq.16) then
            write (90,3016)
            write (nout,3016)
         endif
         if (error.eq.9) then
            write (90,3009)
            write (nout,3009)
         endif
c..   NUCLEAR problem
         if (error.eq.13) then
            write (90,3013)
            write (nout,3013)
         endif
         if (error.eq.25) then
            write (90,3025)
            write (nout,3025)
         endif
         if (error.eq.26) then
            write (90,3026)
            write (nout,3026)
         endif
         if (error.eq.27) then
            write (90,3027)
            write (nout,3027)
         endif
         if (error.eq.14) then
            write (90,3014)
            write (nout,3014)
         endif
         if (error.eq.3) then
            write (90,3003)
            write (nout,3003)
            maxsh = 0
         endif
         if (error.eq.4) then
            write (90,3004)
            write (nout,3004)
         endif
      endif

      ifail = 1
      iexp = 2
      if (error.eq.40.or.error.eq.30) then
         dtn = dtn/(facdt**4)
      else
         if (ireverse.ge.0) then
            if (ireset.lt.abs(icrash)) dtn = dtn/(facdt**iexp)
            if (ireset.eq.abs(icrash)) dtn = dtn0*facdt
            if (ireset.gt.abs(icrash).and.ireset.lt.iresmax) dtn = dtn*
     &           facdt
         else
            if (ireset.lt.abs(icrash)) dtn = dtn*facdt
            if (ireset.eq.abs(icrash)) dtn = dtn0/(facdt**iexp)
            if (ireset.gt.abs(icrash).and.ireset.lt.iresmax) dtn = dtn/
     &           (facdt**iexp)
         endif
      endif
      ireset = ireset+1

      if (ireset.eq.iresmax) goto 60
      if (dtn.lt.dtmin) goto 50


***   in case of failure, adaptative parameter card changes

      if (icorr.ne.'f') then
         if (icorr.eq.'m'.or.icorr.eq.'i'.or.icorr.eq.'h'.or.
     &        icorr.eq.'a') then
            adapt = .false.
         else
            adapt = .true.
         endif

c..   after 2 crashes, mesh is frozen
         if (ireset.ge.1.and.(icorr.eq.'m'.or.icorr.eq.'b'.or.adapt))
     &        then
            maxsh = 0
c     if (ireset.eq.1) dtn = dtn0
            write (nout,400)
            write (90,400)
            ireverse = 1
         else
            maxsh = abs(maxsh0)
            if (maxsh0.lt.0) ireverse = -1
         endif

c..   after 2 crashes, convergence acceleration disabled
         if (ireset.ge.2.and.iacc0.and.iacc.and.(icorr.eq.'h'.or.
     &        adapt)) then
            iacc = .false.
            write (nout,500)
            write (90,500)
         endif

c..   after 2 crashes if end AGB, change tau0
c     if (ireset.eq.3.and.nphase.eq.5.and.ntprof.ne.2.and.(adapt.or.
c     &        icorr.eq.'h').and.(dmenvconv.lt.0.15d0*msun.or.
c     &        dmenvconv/mtini/msun.lt.0.1d0)) then
c     if (tau0.lt.5.d-4) then
c     tau0 = 1.d-3
c     else
c     tau0 = 1.d-4
c     endif
c     dtn = dtn0
c     write (nout,520) tau0
c     write (90,520) tau0
c     endif

c..   after 2 crashes, tolerence on acceleration reduced
         if (hydro.and.(icorr.eq.'i'.or.adapt)) then
            if (ireset.ge.2) then
               eps(1) = eps(1)*1.d1
               write (nout,550) eps(1)
               write (90,550) eps(1)
            else
               eps(1) = eps0(1)
            endif
            if (eps(1).gt.1.d2) then
               eps(1) = eps0(1)
               hydro = .false.
            endif
         endif
      endif

      dtna = dtn
      if (dtn.gt.1.d2) then
         write (90,1500) model,dtn*seci,alpha
         write (nout,1500) model,dtn*seci,alpha
      else
         write (90,1510) model,dtn,alpha
         write (nout,1510) model,dtn,alpha
      endif

      error = 0
      

c..   restore saved model and current time step
      rewind (92)
      rewind (94)
      iflag = ifail
      call rinimod (92,94,92,94,iflag)

      dtn = dtna
      no = no-1
      
      call flush (90)
      return 1

 50   write (nout,1600)
      write (90,1600)
      goto 70
 60   write (nout,1700) model
      write (90,1700) model

 70   write (90,*)
      write (90,*) "failed model:"
      write (90,*)

      write (90,*)
      write (90,*) "stop 'failed'"
      write (90,*)
      write (nout,*)

      return


 75   format (/,' He-core mass changes too much : ',0pf6.3,' [#',i4,
     &     '] --> ',0pf6.3,' [#',i4,']')
 76   format (/,' He-core abundance changes too much : ',0pf6.3,' [#',i4
     &     ,'] --> ',0pf6.3,' [#',i4,']')
 80   format (/,' Temperature increase exceeds ftst = ',i3,'% : T = ',
     &     1pe9.3,' [#',i4,' Mr = ',0pf5.3,'] --> ',1pe9.3,' [#',i4,
     &     ' Mr = ',0pf5.3,']')
 90   format (1x,'WARNING : model ',i8,', Lsurf (',1pe9.3,') > ',
     &     'Eddington limit (Ledd = ',1pe9.3,')')
 100  format (1x,'WARNING : parameter card change: idiffcc = ',l1,
     &     'set to',l1)
 300  format (1x,'WARNING : parameter card change: icrash = 9')
 400  format (1x,'WARNING : parameter card change: maxsh = 0')
 500  format (1x,'WARNING : parameter card change: iacc = .false.')
 520  format (1x,'WARNING : parameter card change: tau0 = ',0pf7.5,
     &     ' and time step reset')
 550  format (1x,'WARNING : parameter card change: u = ',1pe8.2)
 600  format (1x,' pulse : fts changed to ',1pe8.2,' and maxmod ',
     &     'decreased to ',i3)
 650  format (1x,'WARNING :  pulse : fts changed to ',1pe8.2)
 700  format (1x,'WARNING :  3DUP  : fts changed to ',1pe8.2,
     &     ' and maxmod decreased to ',i3)
 750  format (1x,'WARNING :  3DUP  : fts changed to ',1pe8.2)
 1100 format (/,3x,'model',5x,'L',8x,'R',6x,'Teff',7x,'M',9x,'dt (s)',
     &     8x,'t (s)',/,75('-'),/,i8,1x,f10.2,1x,f7.2,1x,f7.0,1x,
     &     f11.8,1x,1pe10.4,1x,1pe16.10,/)
 1200 format (/,3x,'model',5x,'L',8x,'R',6x,'Teff',7x,'M',8x,'dt (yr)',
     &     7x,'t (yr)',/,75('-'),/,i8,1x,f10.2,1x,f7.2,1x,f7.0,1x,
     &     f8.3,1x,1pe10.4,1x,1pe16.10,/)
 1300 format (/,10x,'model ',i8,' failed : too many iterations',/)
 1400 format (/,10x,'new model computation with a reduced time-step,',
     &     ' due to a convergence problem for angular momentum',
     &     ' transport',/)
 1500 format (/,1x,'retry model ',i8,': dtn = ',1pe11.5,
     &     ' yr, with alpha = ',0pf6.4,/)
 1510 format (/,1x,'retry model ',i8,': dtn = ',1pe11.5,
     &     ' s, with alpha = ',0pf6.4,/)
 1555 format (/,1x,'retry model ',i8,':Too large variation of Vsurf',/)
 1600 format (//,5x,'model ',i8,' failed : time-step too low !',/)
 1666 format  (/,1x,'retry model ',i8,'because L > Ledd',/)
 1700 format (//,5x,'model ',i8,' failed : maximum number of crashes ',
     &     ' reached',/)
 1800 format (/,1x,'WARNING : model ',i8,', large variation of the ',
     &     'TOTAL nuclear luminosity (',i4,'%) : time step should be ',
     &     'reduced                  !')
 1810 format (/,1x,'WARNING : model ',i8,', too large variation of ',
     &     'the TOTAL nuclear luminosity (',i4,'% > ftnuc = ',i3,'%),')
c     &     ' mixopt changed : t')
 1850 format (/,1x,'WARNING : model ',i8,', too large variation of ',
     &     'the TOTAL nuclear luminosity (> 500%), nmixd changed : ',
     &     i1,' --> 5, nmaxmod changed : ',i3,' --> ',i3)
 1900 format (132('='))
 2000 format ('save model because of too large H depletion')
 2100 format ('save model because of too large He depletion')

c..   error labels

 3001 format (3x,'MLT : too many convective zones')
 3002 format (3x,'max iterations reached for MLT with hro')
 3005 format (3x,'Pb in NWRMAT')
 3006 format (3x,'Corrections too big')
 3007 format (3x,'KAPPA : Mass fractions exceed unity')
 3080 format (3x,'KAPPA : pb table 2')
c     3081 format (3x,'KAPPA : pb table 3 - O rich')
 3082 format (3x,'KAPPA : pb table 3 - N rich')
 3083 format (3x,'KAPPA : pb table 4')
 3084 format (3x,'KAPPA : pb molecular opacity')
 3010 format (3x,'EOS : rho outside reasonable limits')
 3011 format (3x,'EOS : pb with ionization of H')
 3012 format (3x,'EOS : pb with ionization of heavy elements')
 3020 format (3x,'DIFFSOLVE : division by zero at center !')
 3021 format (3x,'DIFFSOLVE : division by zero in the interior')
 3022 format (3x,'DIFFSOLVE : division by zero at surface !')
 3023 format (3x,'DIFFSOLVE : too many iterations')
 3015 format (3x,'ELEMENTS_FINIS : no solution found')
 3016 format (3x,'DIFFUSION : tstep too large')
 3009 format (3x,'NWRMAT : pb with convective diffusion')
 3013 format (3x,'NETDIFF : maximum iterations and/or time-step ',
     &     'too low !')
 3025 format (3x,'NETDIFF : negative abundance')
 3026 format (3x,'NETDIFF : renormalisation too large')
 3027 format (3x,'NETDIFF : maximum iterations and/or non-convergence')
 3014 format (3x,'NUCSOLVE : maximum iterations and/or non-convergence')
 3003 format (3x,'Pb in NETWORK, maxsh set to ZERO')
 3004 format (3x,'Pb in CHEDIF')

      end


************************************************************************

      SUBROUTINE centeq

************************************************************************
* Give the central boundary conditions for stellar structure           *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.eng'
      include 'evolcom.eos'
      include 'evolcom.grad'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.opa'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.var'

      logical vlconv

      integer k,l

      double precision dtninv,dmkb,dlnt1,dt4,sph,gmrua1,plo
      double precision ablainv
      double precision lconv,lrad
      double precision ft,fp

      common /resolution/ vlconv
      common /ftfp/ ft(nsh),fp(nsh)

      dtninv = 1.d0/dtn

      sph = 3.d0/(pim4*ro(1))
      dmkb = (dm(2)+dm(1))*0.5d0
      dlnt1 = log(t(2)/t(1))

      forall (k=1:neq,l=1:neq1) eq(k,l) = 0.d0

***   central velocity
      if (hydro) then
         eq(1,16) = u(1)
      endif
      eq(1,1) = 1.d0

***   central radius
c      eq(2,16) = r(1)
      eq(2,2) = r(1)

***   mass conservation
      eq(3,16) = r(2)**3/dm(1)-sph
      eq(3,3) = sph*drodf(1)/ro(1)
      eq(3,4) = sph*drodt(1)/ro(1)
      eq(3,7) = 3.d0*r(2)**3/dm(1)

***   central luminosity
      eq(4,16) = lum(1)
      eq(4,5) = 1.d0

***   equation of transport
      dt4 = t(2)**4-t(1)**4
      if (abs(dt4).lt.1.d-5) dt4 = dt4+1.d-5
      if (vlconv) then
         lrad = -64.d0*sig*pi**2*r(2)**4*dt4/(3.d0*kapm(2)*dmkb)
         lconv = vhconv(2)*rom(2)*tm(2)*r(2)

         eq(5,16) = lum(2)-lrad-lconv
         eq(5,3) = wi(2)*(lrad*dkapdf(1)/kap(1)-lconv*drodf(1)/ro(1))
         eq(5,4) = lrad*(wi(2)*dkapdt(1)/kap(1)+4.d0*t(1)**4/dt4)-
     &        lconv*wi(2)*(1.d0+drodt(1)/ro(1))
         eq(5,7) =  -(4.d0*lrad+lconv)
         eq(5,8) = wj(2)*(lrad*dkapdf(2)/kap(2)-lconv*drodf(2)/ro(2))
         eq(5,9) = lrad*(wj(2)*dkapdt(2)/kap(2)-4.d0*t(2)**4/dt4)-
     &        lconv*wj(2)*(1.d0+drodt(2)/ro(2))
         eq(5,10) = 1.d0
      else
         if (crz(2).lt.-1) then
            gmrua1 = gmr(2)*fp(2)+accel(2)
            plo = gmrua1*abla(2)*dmkb/(pim4*r(2)*r(2)*pm(2))
            ablainv = 1.d0/abla(2)
            eq(5,16) = dlnt1+plo
            eq(5,1) = plo*abdu1(2)*ablainv
            eq(5,3) = plo*(abdf1(2)*ablainv-wi(2)*dpdf(1)/p(1))
            eq(5,4) = plo*(abdt1(2)*ablainv-wi(2)*dpdt(1)/p(1))-1.d0
            eq(5,6) = plo*(abdu2(2)*ablainv+dynfac*dtninv/gmrua1)
            eq(5,7) = plo*(abdr(2)*ablainv*r(2)-2.d0*(1.d0+gmr(2)/
     &           gmrua1))
            eq(5,8) = plo*(abdf2(2)*ablainv-wj(2)*dpdf(2)/p(2))
            eq(5,9) = plo*(abdt2(2)*ablainv-wj(2)*dpdt(2)/p(2))+1.d0
            eq(5,10) = plo*abdl(2)*ablainv
         else
            plo = 3.d0*dmkb/(256.d0*sig*pi**2)*lum(2)*kapm(2)/
     &           (r(2)*tm(2))**4*ft(2)
            eq(5,16) = dlnt1+plo
            eq(5,3) = plo*wi(2)*dkapdf(1)/kap(1)
            eq(5,4) = plo*wi(2)*(dkapdt(1)/kap(1)-4.d0)-1.d0
            eq(5,7) = -4.d0*plo
            eq(5,8) = plo*wj(2)*dkapdf(2)/kap(2)
            eq(5,9) = plo*wj(2)*(dkapdt(2)/kap(2)-4.d0)+1.d0
            eq(5,10) = plo/lum(2)
         endif
      endif

      return
      end


***********************************************************************

      SUBROUTINE chedif (error)

************************************************************************
* Calculate the chemical evolution due to nuclear reactions, mixing    *
* and slow transport processes                                         *
* 17/09/03 In case of rotation, nucleosynthesis is computed in a       *
* radiative way in the ENTIRE star. Homogeneisation of CZ is done      *
* afterwards in the DIFFUSION routine                                  *
*  Coupling nucleosynthesis-diffusion equations                        *
*  nmixd = 1 : core only                                               *
*  nmixd = 2 : enveloppe only                                          *
*  nmixd = 3 : core+enveloppe                                          *
*  nmixd = 4 : all convective zones                                    *
*  nmixd = 5 : entire star                                             *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eos'
      include 'evolcom.ion'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer imin,imax,mmean
      integer kl,i,l,k,ij
      integer error
      integer klenv,klcore,klpulse

      double precision xi(nsp),xx(nsh,nsp)
      double precision dturb,dmic,vdmic,lgdturb,lgdmic

      logical nucl(nsh),ipmix(nmaxconv),pass,partialmix

      common /coefdiff/ dturb(nsh),dmic(nsh),vdmic(nsh),lgdturb(nsh)
     &     ,lgdmic(nsh)
      common /overshoot/ klcore,klenv,klpulse


      if (t(1).eq.t(2)) then
         error = 4
         return
      endif

*___________________________________________________________________
***   various abundance initializations (in case of call to ydequil)
*-------------------------------------------------------------------

      if (dmaccr.gt.0.d0.and.xspacc(ih2).gt.1.d-10.and.accphase.eq.0)
     &     then
         do k = iaccbot,iacctop
            vxsp(k,ih1) = xsp(k,ih1)
            vxsp(k,ih2) = xsp(k,ih2)
            vxsp(k,ihe3) = xsp(k,ihe3)
         enddo
      endif
      do k = 1,nmod
         nucl(k) = .false.
      enddo

*-----------------------------------------------------------------------*
***  coupling nucleosynthesis and time-dependent mixing in the all star *
*-----------------------------------------------------------------------*

c..  compute diffusion coefficients
      if (nmixd.gt.0) then
         call diffusion (1,error)
         if (error.gt.0) return
      endif
      if (nmixd.eq.5) then
         kcentd = 1
         ksurfd = 1
         write (nout,*) 'Coupling diffusion/nucleosynthesis equations ',
     &        'throughout the entire structure'
         call netdiff (1,nmod1,dturb,dtn,error)
         return
      endif


*_________________________________________
***   nucleosynthesis in convective shells
*-----------------------------------------


***   Instantaneous mixing in convective zones if NO DIFFUSION
      pass = .not.diffzc
cx      if (flame.and.diffzc) pass = .true.
      if (nsconv.gt.0.and.nmixd.lt.4.and.pass) then
c..   restore initial abundances
         do l = 1,nsp
            do k = 1,nmod
               xsp(k,l) = vxsp(k,l)
               ysp(k,l) = xsp(k,l)/anuc(l)
            enddo
         enddo
         do 10 kl = 1,nsconv
            if ((nmixd.eq.1.and.kl.eq.klcore).or.
     &           (nmixd.eq.2.and.kl.eq.klenv).or.
     &           (nmixd.eq.3.and.(kl.eq.klcore.or.kl.eq.klenv))) goto 10
            imin = novlim(kl,3)
            imax = novlim(kl,4)
            am(kl) = 0.d0
            do i = imin,imax
               am(kl) = am(kl)+dm(i)
            enddo
            am(kl) = 1.d0/am(kl)
            if (nucreg.eq.2.and..not.flame) then
               call mix (dm,imin,imax,kl,5,partialmix)
               ipmix(kl) = .true.
            else
               call mix (dm,imin,imax,kl,2,partialmix)
               ipmix(kl) = partialmix
            endif
            if (.not.ipmix(kl)) then
               if (no.eq.1) write (90,100) imin,imax
               write (nout,100) imin,imax
               mmean = int((imin+imax)/2)
               do l = 1,nsp
                  xi(l) = xsp(mmean,l)
               enddo
               call network ('c',xi,am(kl),dtn,tnucmin,imin,imax,error)
               if (error.gt.0) return
               do k = imin,imax
                  nucl(k) = .true.
                  do l = 1,nsp
                     xsp(k,l) = xi(l)
                  enddo
c..  update composition
                  if (nucreg.ne.3) then
                     do l = 1,nsp
                        vxsp(k,l) = xi(l)
                     enddo
                  endif
               enddo
            else
***   in case of partial mixing, convective zone treated as radiative
               if (no.eq.1) write (90,200) imin,imax
               write (nout,200) imin,imax
               do k = imin,imax
                  ij = k
                  do l = 1,nsp
                     xi(l) = xsp(k,l)
                  enddo
                  call network ('r',xi,0.d0,dtn,tnucmin,ij,ij,error)
                  if (error.gt.0) return
                  nucl(k) = .true.
                  if (nucreg.eq.3) then
                     do l = 1,nsp
                        xx(k,l) = vxsp(k,l)
                     enddo
                  endif
c..  update composition
                  do l = 1,nsp
                     xsp(k,l) = xi(l)
                     vxsp(k,l) = xi(l)
                  enddo
               enddo
               call mix (dm,imin,imax,kl,3,partialmix)
               if (nucreg.eq.3) then
                  do l = 1,nsp
                     do k = imin,imax
                        vxsp(k,l) = xx(k,l)
                     enddo
                  enddo
               endif
             endif
 10      continue
      endif


*-----------------------------------------------------------------------*
***  nucleosynthesis calculation with time-dependent convective mixing  *
*-----------------------------------------------------------------------*


      if (nmixd.gt.0.and.nsconv.gt.0) then
         kcentd = 1
         ksurfd = 1
         if ((nmixd.eq.1.or.nmixd.eq.3).and.klcore.eq.1) then
            write (nout,400) novlim(1,8)
            imax = novlim(1,8)
            call netdiff (1,imax,dturb,dtn,error)
            do k = 1,imax
               nucl(k) = .true.
            enddo
         endif
         if ((nmixd.eq.2.or.nmixd.eq.3).and.klenv.gt.0) then
            write (nout,500) novlim(klenv,7),nmod1
            imin = novlim(klenv,7)
            call netdiff (imin,nmod1,dturb,dtn,error)
            do k = imin,nmod
               nucl(k) = .true.
            enddo
         endif
         if (nmixd.eq.4) then
            do kl = 1,nsconv
               imin = max(1,novlim(kl,7))
               imax = min(nmod1,novlim(kl,8)+1)
               if (kl.eq.klenv) imax = nmod1
               write (nout,600) imin,imax
               call netdiff (imin,imax,dturb,dtn,error)
               do k = imin,imax
                  nucl(k) = .true.
               enddo
            enddo
         endif
      endif

*_________________________________________________________
***   nucleosynthesis in radiative and thermohaline shells
*---------------------------------------------------------

      if (no.eq.1) write (90,*) 'Processing radiative nucleosynthesis'
      write (nout,*) 'Processing radiative nucleosynthesis'
      do k = 1,nmod
         if (.not.nucl(k)) then
            ij = k
            do l = 1,nsp
               xi(l) = xsp(k,l)
            enddo
            call network ('r',xi,0.d0,dtn,tnucmin,ij,ij,error)
            if (error.gt.0) return
            nucl(k) = .true.
            do l = 1,nsp
               xsp(k,l) = xi(l)
            enddo
c..  update composition
            if (nucreg.ne.3) then
               do l = 1,nsp
                  vxsp(k,l) = xi(l)
               enddo
            endif
         endif
      enddo


 100  format (' Processing one-zone convective nucleosynthesis in [',
     &     i4,',',i4,']')
 200  format (' Processing convective zone : [',i4,',',i4,
     &     '] as a radiative one')
 400  format (' Coupling diffusion/nucleosynthesis equations in the ',
     &     'core + overshoot regions : [ 1,',i4,']')
 500  format (' Coupling diffusion/nucleosynthesis equations in the ',
     &     'envelope + overshoot regions : [',i4,',',i4,']')
 600  format (' Coupling diffusion/nucleosynthesis equations in ',
     &     'convective zone + overshoot regions : [',i4,',',i4,']')

      return
      end


************************************************************************

      SUBROUTINE convzone

************************************************************************
*                     Delineate the convective boundaries
* novlim(kl,1-2)   -->  Schwarzschild criterion
* novlim(kl,3-4)   -->  CURRENT LIMITS and parametric overshoot   ! note td : 3 - base ZC et 4 - haut ZC (nsconv : nbre ZC)
* novlim(kl,5-6)   -->  Ledoux criterion
* novlim(kl,7-8)   -->  diffusive overshoot and mesh
* novlim(kl,9-10)  -->  semiconvective zones
* novlim(kl,11-12) -->  thermohaline zones
*
* crzc = c  -2  : convective shell      abled < abrad or abadd < abrad
* crzc = a  -1  : atmospheric shell     tau < 10
* crzc = s  -3  : semi-convective shell abadd < abrad < abled  (Ledoux)
*                 shell mixed : semiconv+diffzc
*        S   3  : if Ledoux criterion used, shell *not* mixed (semiconv)
* crzc = t   2  : thermohaline shell    abled < abrad < abadd
* crzc = o   1  : overshoot layers
* crzc = r   4  : radiative shell       abrad < abled or abrad < abadd
*
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
C... modif NL 
      include 'evolcom.diff'
      include 'evolcom.therm2'
      include 'evolcom.transp'
      include 'evolcom.nuc'
      include 'evolcom.spec'
C...fin modif 
      include 'evolcom.grad'
      include 'evolcom.mod'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer kl,nconv1,nconv2,nrad,nshdel
      integer delnconv,nc,ncz
      integer nsconv0,novlim0
      integer i,j,k,jl,it
      integer ishift,kc,kc0,kc1,imin
      integer klenv,klpulse,klcore
C... modif NL 
      integer nmodmdix,id,ip,nshell,ibth,itth
C.. fin modif 
      

      double precision extabrad
      double precision delabrad
      double precision delmconv
C... modif NL
      double precision abmuliss,abmax,abmin,gradledoux
C...fin modif 

      logical isconv,isconvS,isconvT,isconvTnuc,isconvL,
     &       tbot,ttop,neutrality

      dimension isconv(nsh),isconvL(nsh),isconvS(nsh),isconvT(nsh),
     &     novlim0(nmaxconv,ntypeconv),isconvTnuc(nsh)
C... modif NL
      dimension abmuliss(nsh),gradledoux(nsh)
C...fin modif      
      
      common /overshoot/ klcore,klenv,klpulse

***   define gradient differences
      do i = 1,nmod1
         gradledoux(i) = abled(i)-abrad(i)
         isconv(i) = abrad(i).gt.abm(i)   !Schwarzchild
         isconvL(i) = abrad(i).gt.abled(i).and.isconv(i)  !Ledoux
         isconvS(i) = isconv(i).and.abmu(i).gt.0.d0  !semiconvection
	 isconvT(i) = .not.isconv(i).and.abmu(i).lt.0.d0.
     &                  and.gradledoux(i).gt.0.d0  !thermohaline
	 isconvTnuc(i) = .not.isconv(i).and.abmu(i).lt.0.d0.
     &                   and.gradledoux(i).gt.0.d0.
     &                   and.xsp(i,ihe3).ge.1.d-08 !thermohaline
     
     
      enddo
      isconv(nmod) = .false.
      isconvS(nmod) = .false.
      isconvT(nmod) = .false.
      isconvTnuc(nmod) = .false.
      isconvL(nmod) = .false. 

      scr = 'r'
      do i = 1,nmod
         crz(i) = 4
      enddo

***   define novlim
      nsconv0 = nsconv
      do kl = 1,ntypeconv
         do j = 1,nmaxconv
            novlim0(j,kl) = novlim(j,kl)
            novlim(j,kl) = 0
         enddo
      enddo
      nsconv = 0
      nschwar = 0
      nledoux = 0
      nsemiconv = 0
      nthermha = 0

***   define Schwarzschild and Ledoux regions
      do i = 1,nmod
         if (isconv(i)) crz(i) = -2
      enddo
      it = -3
      if (ledouxmlt) it = 3
      do j = 1,2
         kl = 1
         jl = 1
         if (j.eq.2) then
            do i = 1,nmod
               if (isconvS(i)) crz(i) = it
            enddo
            jl = 5
         endif
         if (isconv(1)) novlim(kl,jl) = 1
         do k = 2,nmod1
            tbot = crz(k).ne.-2.and.crz(k+1).eq.-2
            ttop = crz(k).eq.-2.and.crz(k+1).ne.-2
            if (tbot) novlim(kl,jl) = k+1
            if (ttop.and.k-novlim(kl,jl).gt.nczm) then
               novlim(kl,jl+1) = k
               kl = kl+1
               if (kl.gt.nmaxconv) then
                  write (nout,*) 'too many convective zones !',j
                  goto 10
               endif
            endif
         enddo
 10      if (j.eq.1) nschwar = kl-1
         if (j.eq.2) nledoux = kl-1
      enddo

***   define semiconvective and thermohaline zones
      do i = 1,nmod
         if (isconvT(i)) crz(i) = 2
      enddo
      do j = 1,2
         kl = 1
         jl = 2*j+7
         if (j.eq.2) it = 2
         if (crz(1).eq.it) novlim(kl,jl) = 1
         do k = 2,nmod1
            tbot = crz(k).ne.it.and.crz(k+1).eq.it
            ttop = crz(k).eq.it.and.crz(k+1).ne.it
            if (tbot) novlim(kl,jl) = k+1
            if (ttop.and.k-novlim(kl,jl).ge.nczm) then
               novlim(kl,jl+1) = k
               kl = kl+1
               if (kl.gt.nmaxconv) then
                  write (nout,*) kl,'too many unstable zones !',j
                  goto 20
               endif
            endif
         enddo
 20      if (j.eq.1) nsemiconv = kl-1
         if (j.eq.2) nthermha = kl-1
      enddo

      if (ledouxmlt) then
         nsconv = nledoux
         do kl = 1,nsconv
            novlim(kl,3) = novlim(kl,5)
            novlim(kl,4) = novlim(kl,6)
         enddo
      else
         nsconv = nschwar
         do kl = 1,nsconv
            novlim(kl,3) = novlim(kl,1)
            novlim(kl,4) = novlim(kl,2)
         enddo
      endif

C... modif NL 
      do i= 2,nmod 
         if (.not.isconvTnuc(i-1).and.isconvTnuc(i)) then 
	     novlim(1,11)= i
         goto 22 
	 endif 
      enddo
 22   do i= 1,nmod-1 
	 if (.not.isconvTnuc(i+1).and.isconvTnuc(i)) then 
	     novlim(1,12)= i
         goto 25 
	 endif
      enddo       
 

*___________________________________________________________________
***   suppress convective zones whose width is less than nczm shells
***   additional suppression during agbphase if shell mass too small
*-------------------------------------------------------------------

c      if (nphase.ge.5.and.nsconv.ge.2.and.nczm.gt.0) then
 25   ncz = abs(nczm)
      if (nsconv.ge.2.and.nczm.ne.0) then
         k = 0
 50      k = k+1
         delmconv = 1.d99
         if (agbphase) delmconv = (m(novlim(k,4))-m(novlim(k,3)))/msun
         if (nczm.gt.0.or.t(novlim(k,4)).lt.1.d6) then
            delnconv = novlim(k,4)-novlim(k,3)
         else
            delnconv = ncz+1
         endif
         if (delmconv.lt.delmcz.or.delnconv.le.ncz) then
c..   the suppressed convective shells are assumed to be radiative
            do i = novlim(k,3),novlim(k,4)
c               crz(i) = 4
               if (crz(i).eq.-2) crz(i) = 4
            enddo
            do kl = k,nsconv-1
               novlim(kl,3) = novlim(kl+1,3)
               novlim(kl,4) = novlim(kl+1,4)
            enddo
            novlim(nsconv,4) = 0
            novlim(nsconv,3) = 0
            k = k-1
            nsconv = nsconv-1
         endif
         if (k.lt.nsconv) goto 50
      endif

*____________________________________________________
***   merging of adjacent convective zones,
***   if separated by less than nczm radiative shells
*----------------------------------------------------

      if (nphase.gt.1.and.nsconv.gt.1.and.nczm.ne.0) then
         nc = nsconv
         do kl = 1,nsconv-1
            nconv1 = novlim(kl,4)-novlim(kl,3)+1
            nconv2 = novlim(kl+1,4)-novlim(kl+1,3)+1
            if (nczm.gt.0.or.t(novlim(kl,4)).lt.1.d6) then
               nrad = novlim(kl+1,3)-novlim(kl,4)-1
            else
               nrad = ncz+1
            endif
            if (nrad.lt.nconv1.and.nrad.le.ncz.and.nrad.le.nconv2)
     &           then
               delabrad = abrad(novlim(kl+1,3))-abrad(novlim(kl,4))
               nshdel = novlim(kl+1,3)-novlim(kl,4)
               do k = novlim(kl,4)+1,novlim(kl+1,3)-1
                  abrad(k) = abrad(novlim(kl,4))+delabrad/dble(nshdel)*
     &                 dble(k-novlim(kl,4))
c                  abrad(k) = max(abrad(k),abm(k)*1.001d0)
                  crz(k) = 4
               enddo
               novlim(kl,4) = novlim(kl+1,4)
               do j = kl+1,nc-1
                  novlim(j,3) = novlim(j+1,3)
                  novlim(j,4) = novlim(j+1,4)
               enddo
               novlim(nc,4) = 0
               novlim(nc,3) = 0
               nc = nc-1
               if (nc.eq.nsconv-1) goto 60
            endif
         enddo
 60      nsconv = nc
      endif

      do kl = 1,nsconv
         novlim(kl,7) = novlim(kl,3)
         novlim(kl,8) = novlim(kl,4)
      enddo
      if (nsconv.ge.1) then
         if (ledouxmlt) then
            scr = 'l'
         else
            scr = 'c'
         endif
      endif

***   define atmospheric layers
      if (nsconv.gt.0) then
         if (t(novlim(nsconv,4)).lt.1.d6) then
            imin = novlim(nsconv,4)
             do i = imin,nmod
               crz(i) = -1
               Dconv(i) = Dconv(imin-1)
            enddo
         endif
      endif

*____________________________________________________________________
***   search for gradient neutrality (abrad = abad) at all convective 
***   boundaries during the first itermix iterations
*--------------------------------------------------------------------

c      neutrality = idup.eq.1.or.idup.eq.3  
      if (nphase.eq.2) neutrality = .true. ! modif Ana Palacios 09/2018
      if (nsconv.gt.0.and.neutrality.and.iter.le.abs(itermix)) then
         do kl = 3,4
            if (kl.eq.3) then
               ishift = 1
            else
               ishift = -1
            endif
            do k = 1,nsconv
               kc = novlim(k,kl)
c..  do not treat upper part of conv. envelope
               if (k.eq.klenv.and.kl.eq.4) kc = 0
               if (kc.gt.1) then
                  kc1 = kc+ishift
                  kc0 = kc-ishift
                  if (ishift.eq.1) then
                     delabrad = dm(kc0)/dm(kc)
                  else
                     delabrad = dm(kc)/dm(kc1)
                  endif
                  extabrad = abrad(kc)-delabrad*(abrad(kc1)-abrad(kc))
c                  if (extabrad.gt.abm(kc0).and.extabrad.lt.abrad(kc))
c     &                 then
	          if (extabrad.gt.abm(kc0)) then
                     novlim(k,kl) = kc0
                     novlim(k,kl-2) = novlim(k,kl)
                     abrad(kc0) = extabrad
                     crz(kc0) = -2
                  endif
               endif
            enddo
         enddo
      endif

c..   check if convective boundaries have changed
      if (chkmix) then
         if (nsconv.ne.nsconv0) then
            lmix = .true.
         else
            do j = 1,nsconv
               do kl = 3,4
                  if (novlim0(j,kl).ne.novlim(j,kl)) then
                     lmix = .true.
                     return
                  endif
               enddo
            enddo
         endif
      endif

 2000 format (i8,2x,10(1x,0pf10.7),2x,1x,i2,a1,a2)

      return
      end



************************************************************************

      SUBROUTINE denucl (tk,rok,mueinvk,x,ksh,eng,engdt,engdro)

************************************************************************
* Compute the nuclear energy production rates and its derivatives      *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      integer nflu,ksh

      double precision tk,rok,mueinvk,x
      double precision tkp,rokp,eng,engro,engdt,engdro,engrocap
      double precision v
      double precision dtk,drok,engro1,eng1,engrocap1,engro2,eng2,
     &     engrocap2

      dimension x(nsp),v(nreac)

      if (tk.ge.1.d5) then
***   rho,T+dT
         dtk = tk*1.d-4
         tkp = tk+dtk
         nflu = 1
         call vit (tkp,rok,mueinvk,ksh,v,0)
         call nuceng (nflu,x,v,ksh,eng1,engro1,engrocap1)
***   rho+drho,T
         nflu = 2
         drok = rok*1.d-4
         rokp = rok+drok
         call vit (tk,rokp,mueinvk,ksh,v,2)
         call nuceng (nflu,x,v,ksh,eng2,engro2,engrocap2)
***   rho,T
         nflu = 0
         call vit (tk,rok,mueinvk,ksh,v,0)
         call nuceng (nflu,x,v,ksh,eng,engro,engrocap)
         engdt = (eng1-eng)/dtk
         engdro = engro/rok+(engrocap2-engrocap)/drok
      else
***   rho,T
         nflu = 0
         call vit (tk,rok,mueinvk,ksh,v,0)
         call nuceng (nflu,x,v,ksh,eng,engro,engrocap)
         engdt = 0.d0
         engdro = 0.d0
      endif

      return
      end


************************************************************************

      DOUBLE PRECISION FUNCTION dfdr (i,h,f,schwarr,idet)

************************************************************************
* Calculate the derivative of the penetration function                 *
* for the accretion process                                            *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.grad'
      include 'evolcom.ion'
      include 'evolcom.rot'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer i,idet,isup

      double precision h,f,schwarr
      double precision rr,muacc,ri,vf,b0,gradmu,b1,c,det,fthacc

      isup = nmod
      dfdr = 1.d-21
      if (idet.eq.0) return
      rr = r(i)+h
      fthacc = f*thacc
      if (iaccr.eq.1.or.iaccr.eq.3) muacc = mueinv(i)+fthacc*muiacc+
     &     (1.d0-fthacc)*muiinv(i)
      if (iaccr.eq.2.or.iaccr.eq.4) muacc = mueinv(i)+(fthacc*
     &     muiacc+muiinv(i))/(1.d0+fthacc)
      if (abrad(i).gt.abad(i)) then
         schwarr = -abs(schwarr)
         ri = ric
      else
         ri = rir
      endif

      vf = vomega(i)*r(i)*r(i)/sigma0
      if (iaccr.eq.1.or.iaccr.eq.3) then
         b0 = -thacc*(muiacc-muiinv(i))*phiKS(i)/muacc
         gradmu = abmu(i)*(1.d0-thacc*f)
      endif
      if (iaccr.eq.2.or.iaccr.eq.4) then
         b0 = -thacc*(muiacc-muiinv(i))/(1.d0+thacc*f)**2*phiKS(i)/
     &        muacc
         gradmu = abmu(i)/(1.d0+fthacc)
      endif
      if (ri.eq.ric) then
         gradmu = 0.d0
         abmu(i) = 0.d0
      endif
      b1 = m(i)/(thacc*m(isup)*xiaccr*xiaccr)-(f-vf)**2*r(isup)/rr
      b0 = b1*b0-2.d0*(f-vf)*r(isup)/rr
      c = (-1.d0/(mu(i)*muacc)*gradmu*phiKS(i)-prrc*deltaKS(i)*
     &     schwarr)/hp(i)*b1
      det = b0*b0-4.d0*r(isup)*ri*c
      if (det.lt.0.d0) then
         idet = 0
         return
      endif

      dfdr = abs(b0/(2.d0*r(isup)*ri)-abs(dsqrt(det)/(2.d0*r(isup)*ri)))
      if (dfdr.lt.1.d-20) dfdr = 1.d-21

      return
      end
      SUBROUTINE diffinit (icall)

************************************************************************
* Initialisation of variables for radiative diffusion computations     *
*
* icall = 0 : call during convergence (mixopt=t), do not treat rotational
*             mixing and do not update vxsp
* icall = 1 : call after convergence
* icall = 2 : for binlist_evol.f
*                                                                      *
* $LastChangedDate:: 2016-05-11 17:20:46 +0200 (Mer, 11 mai 2016)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 62                                                          $ *
*                                                                      *
************************************************************************
c.. 20/10: Correction des expressions de khi_mu et epsimu

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.grad'
      include 'evolcom.ion'
      include 'evolcom.igw'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.opa'
      include 'evolcom.conv'
      include 'evolcom.spec'
      include 'evolcom.rot'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.transp'
      include 'evolcom.var'

      integer i,j,k
      integer nmodmdix,icall,nczm0
      integer id

      double precision delradm
      double precision vrray,vom,vxpsi
      double precision xnurad,xnumol,ddmax
      double precision abmuj,abmurj,abmax,abmin
      double precision tamp

      common /gradmuj/ abmurj(nsh),abmuj(nsh)
      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)
      common /diffvisc/ xnumol(nsh),xnurad(nsh)

*___________________
***   initialisation
*-------------------

      do i = 1,nmod
         muizero(i) = 0.d0
         epsimu(i) = 0.d0
         do j = 1,nis   
            muizero(i) = muizero(i)+xsp(i,j)/anuc(j)  ! 1 / mui physique
         enddo        
         muizero(i) = 1.d0/muizero(i)    ! mui
         khi(i) = 3.02427122d-4*t(i)**3/(kap(i)*ro(i))
      enddo
      do i = 2,nmod
         phiKSm(i) = phiKS(i-1)*wi(i)+phiKS(i)*wj(i)
         khim(i) = khi(i-1)**wi(i)*khi(i)**wj(i)
      enddo
      phiKSm(1) = phiKS(1)
      khim(1) = khi(1)

      do i = 2,nmod
         hhp(i) = hp(i)
         hht(i) = hhp(i)/abrad(i)
      enddo
      hhp(1) = hp(1)
      delradm = 0.5d0*(abrad(2)+abrad(1))
      hht(1) = hhp(1)/delradm

      if (.not.(icall.eq.1.and.rotation).or.icall.eq.2) return

c..   remove/merge small convective zones
      nczm0 = nczm
      nczm = max(nczm,10)
      call convzone
      nczm = nczm0

*_______________________________
***   density redefinition
*-------------------------------
      do i = 1,nmod1
         rro(i) = pw34*dm(i)/(pi*(r(i+1)**3-r(i)**3))
         vrro(i) = pw34*dm(i)/(pi*(vr(i+1)**3-vr(i)**3))
      enddo
      rro(nmod) = ro(nmod)
      vrro(nmod) = vro(nmod)


*__________________________________________________________
*** smooth mean molecular weight gradient needed for lambda
*----------------------------------------------------------

      rtot = r(nmod)
      gs = gmr(nmod)

      abmax = -1.d99
      abmin = 1.d99
      do i = 1,nmod
         abmax = max(abmax,abmu(i))
         abmin = min(abmin,abmu(i))
         abmuj(i) = abmu(i)
         abmurj(i) = -abmuj(i)/hp(i)
      enddo

c..  remove noise
      abmax = abmax*1.d-8
      abmin = abmin*1.d-8
      do i = 1,nmod
         if (abmu(i).gt.0.d0.and.abmu(i).le.abmax) then
            abmuj(i) = 0.d0
            abmurj(i) = 0.d0
         endif
         if (abmu(i).lt.0.d0.and.abmu(i).ge.abmin) then
            abmuj(i) = 0.d0
            abmurj(i) = 0.d0
         endif
      enddo

c.. smooth
      id = 5
      nmodmdix = nmod-10
      do k = 1,10
         call lissage (abmuj,nmodmdix,10,id)
         call lissage (abmurj,nmodmdix,10,id)
      enddo

c..   after smoothing, ensure conv. zones are chemically homogeneous
      if (nsconv.gt.0) then
         do i = 1,nsconv
            do j = novlim(i,3),novlim(i,4)
               abmurj(j) = 0.d0
               abmuj(j) = 0.d0
            enddo
         enddo
      endif
c..   Add a fetch factor fmu to abmu and abmurj as in Chieffi & Limongi (2013)
c      fmu = 1.d0
c      do i = 1,nmod
c         abmurj(i) = fmu * abmurj(i)
c         abmu(i) = fmu * abmu(i)
c      enddo


*___________________________________________________________
***   compute smoothed normalized variables for AM transport
*-----------------------------------------------------------

      rray(1) = 0.d0
      vrray(1) = 0.d0
      grav(1) = gmr(1)/gmr(nmod)
      rkonv(1) = abrad(1)-abad(1)
      rkonvL(1) = rkonv(1)-abmuj(1)*phiKS(1)/deltaKS(1)
      vom(1) = vomega(1)

      do i = 2,nmod
         rray(i) = r(i)/rtot
         vrray(i) = vr(i)/rtot
         grav(i) = gmr(i)/gmr(nmod)
         rkonv(i) = abrad(i)-abm(i) ! abrad = abla in RZ
         rkonvL(i) = rkonv(i)-abmuj(i)*phiKSm(i)/deltaKSm(i)
         vom(i) = vomega(i)
      enddo


*______________________________________________________________________
*** compute logarithmic derivatives of nuclear energy and thermal
*** conductivity with respect to mean molecular weight and temperature
*-----------------------------------------------------------------------

      do i = 1,nmod
         epsit(i) = denucldt(i)/enucl(i)
         khit(i) = 3.d0-dkapdt(i)/kap(i)+deltaKS(i)
         epsimu(i) = ro(i)/enucl(i)*denucldro(i)*phiKS(i)
         khi_mu(i) = -phiKS(i)*(1.d0+ro(i)*dkapdro(i)/kap(i))
      enddo

*_____________________________________________________________
***   computation of quantities for angular momentum transport
*-------------------------------------------------------------

      if (rotation) then
         do i = 2,nmod
            vxpsi(i) = xpsis(i)
         enddo
         vxpsi(1) = xpsis(1)
      endif

*_______________________________________________________
***   Computation of radiative and molecular viscosities
***   numol -> Schatzman 1977, A&A,56
***   nurad -> Kippenhahn & Weigert
*-------------------------------------------------------

      do i = 1,nmod
         xnurad(i) = khim(i)*tm(i)/rom(i)/4.4937758d21
c         xnurad(i) = 6.72991d-26*t(i)**4/(kap(i)*rro(i)**2)
*       val num : 5*c**2
         ddmax = 1.5964182d-8*tm(i)**1.5d0/dsqrt(rom(i))
*       val num : 1.5/e**3*(mp*k**3/pi)**0.5
         xnumol(i) = 2.1688529d-15*tm(i)**2.5d0/rom(i)/log(ddmax)
*       val num : 0.4*mp**0.5*k**2.5/e**4
         xnum(i) = xnumol(i)+xnurad(i)
      enddo

      return
      end
*************************************************************************

      SUBROUTINE diffsolve (abond,ratneh,anuc,znuc,dmix,f1,f2,f3,dmdr,
     &     ielem,nshell,icall,error)

*************************************************************************
*     solve diffusion equation : dX/dt = d/dm(DdX/dt-vdiff*V)
*                                                                      *
* $LastChangedDate:: 2014-02-04 15:45:03 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 11                                                          $ *
*                                                                      *
*************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.diff'
      include 'evolcom.mass'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer, intent (in) :: ielem
      integer ij,ndt,ndb,error
      integer iiter,iitermin,iitermax,ip,im
      integer isolve,icall
      integer i,ijk,nshell,icor,icormin,jshell
      integer mspec             ! Ajout pour form. thoul
      integer elemthoul
      parameter (mspec = 32)    ! Ajout pour form. thoul

      integer ineutron,ih1,ih2,ihe3,ihe4,ili6,ili7,ibe7,ibe9,ib8,ib10, ! Ajout pour form. thoul
     &     ib11,ic12,ic13,ic14,in13,in14,in15,io15,io16,io17,io18,if19,
     &     ine20,ine21,ine22,ina23,img24,img26,ial27,isi28,ip31,is32,
     &	   icl35,idxspc,idxsps,idxspci
      integer ina24,ina25,ine23,ial26g,ial26m,if18,if20,ina22,img27,
     &     is35,icl36,img25

      double precision abond,ratneh,anuc,znuc,tsh,rhosh !Ajout des deux derniers pour form. thoul
      double precision vdiffp,diff13,ration  ! M&M ioni part 09/18
      double precision rray,rkonv,rkonvL,muizero
      double precision, save :: cd(nsh),cdm(nsh)
      double precision xa,d12,dmix,vdiff,vv  !,vdiffthoul !Ajout des deux derniers pour form. thoul
      double precision dmdr,x1,x2,x3,x
      double precision v,w,f1,f2,f3
      double precision tol,alpha,corm,cormin,dx,a1,eqd,dv,cdi,cdip
      double precision vdiffthoul, dmicthoul
      
      logical indiffus

      common /star/ rray(nsh),rkonv(nsh),rkonvL(nsh),muizero(nsh)
      common /thoulindex/ elemthoul(mspec-1)                        ! Ajout TD AP
      common /specindex/ ineutron,ih1,ih2,ihe3,ihe4,ili6,ili7,ibe7,ibe9,                    ! Ajout pour form. thoul
     &     ib8,ib10,ib11,ic12,ic13,ic14,in13,in14,in15,io15,io16,io17
     &     ,io18,if19,ine20,ine21,ine22,ina23,img24,ial27,isi28,ip31,
     &     is32,icl35,ina24,ina25,ine23,img26,ial26g,ial26m,if18,if20,
     &     ina22,img27,is35,icl36,img25,idxspc(nprintc),idxsps(nprint)
     &	   ,idxspci(9)
      common /dmicthoul/ vdiffthoul(nsh,mspec), dmicthoul(nsh,mspec) ! modif

      dimension abond(nsh,nis),xa(nshell),x(nsh),ratneh(nsh),
     &     anuc(nsp),znuc(nsp)      
      dimension d12(nsh),dmix(nsh),eqd(4)
      dimension vdiff(nsh),vv(nshell)
      dimension v(nshell),w(nshell),f1(nsh),f2(nsh),f3(nsh),dmdr(nsh)

      isolve = 1
      ij = ielem
      tol = 1.d-14
      alpha = 0.75d0
      iitermax = 200

c..   initialisations
      if (ielem.eq.2.or.(ielem.eq.3.and.microdiffus)) then
         dmix(nmod) = dmix(nmod1)
         do i = 1,nmod1
            ip = i+1
            cd(i) = dmix(i)
            cdm(ip) = dmix(i)**wi(ip)*dmix(ip)**wj(ip)
         enddo
      endif
         
 100  iiter = 0
      iitermin = 3
      ndb = 1
      ndt = nmod-1
      indiffus = .false.
      icormin = 1
      cormin = 1.d99
      do i = 1,nmod1
         vv(i) = 0.d0
         xa(i) = abond(i,ielem)
         x(i) = xa(i)
      enddo

*_____________________________
***   Newton Raphson procedure
*-----------------------------

  200 continue
      icor = 1
      corm = 0.d0
      iiter=iiter+1

*_____________________________________________
***   Compute microscopic diffusion velocities
*---------------------------------------------

      if (microdiffus.and.icall.eq.1) then
         indiffus = .true.
c..   Diffusion "Chapman & Cowling"
         if (lmicro.eq.2) then
            call diffmic (ij,ndt,ndb,vdiff,d12,anuc,znuc)
c..   Diffusion "Paquette et al. (1986)"
         else if (lmicro.eq.3.or.lmicro.eq.5) then
            call diffpaq (ij,ndt,ndb,vdiff,d12,abond,ratneh)
c..   Diffusion "Thoul et al. (1994)"
         else if (lmicro.eq.4.or.lmicro.eq.6) then
            do i = 2,nmod1            
               if (crz(i).gt.0) then
                  cd(i) = dmix(i)
               endif
               do ijk = 1,mspec-1
                  if (ij.eq.elemthoul(ijk)) then
                     vv(i)=vdiffthoul(i,ijk)*dmdr(i)
                  endif
               enddo
               cdm(i) = cd(i-1)**wi(i)*cd(i)**wj(i)
            enddo
         endif
         Dmicro(1:nsh) = d12(1:nsh)
         vmicro(1:nsh) = vdiff(1:nsh)
         if (lmicro.eq.2.or.lmicro.eq.3.or.lmicro.eq.5) then
            do i = 2,nmod1
               if (crz(i).gt.0) then
                  cd(i) = dmix(i)+d12(i)*dmdr(i)*dmdr(i)
                  vv(i) = vdiff(i)*dmdr(i)
               endif
               cdm(i) = cd(i-1)**wi(i)*cd(i)**wj(i)
            enddo
         endif
      endif

c..   center
      cdip = cdm(2)
      x1 = f1(1)*vv(2)
      x2 = f2(1)*cdip
      eqd(4) = x(1)-xa(1)-x2*(x(2)-x(1))+x1*x(2)
      eqd(1) = 0.d0
      eqd(2) = 1.d0+x2
      eqd(3) = -x2+x1
      if (eqd(2).eq.0.d0) then
         error = 20
         return
      endif
      w(1) = eqd(4)/eqd(2)
      v(1) = eqd(3)/eqd(2)

c..   interior
      do i = 2,nmod1-1
         ip = i+1
         im = i-1
         cdi = cdm(i)
         cdip = cdm(ip)
         x3 = f3(i)*cdi
         x2 = f2(i)*cdip
         eqd(4) = x(i)-xa(i)-x2*(x(ip)-x(i))+x3*(x(i)-x(im))+
     &        f1(i)*(vv(ip)*x(ip)-vv(i)*x(i))
         eqd(1) = -x3
         eqd(2) = 1.d0+x3+x2-f1(i)*vv(i)
         eqd(3) = -x2+f1(i)*vv(ip)
         if (isolve.eq.0) then
            a1 = eqd(2)-eqd(1)*v(im)
            if (a1.eq.0.d0) a1 = -eqd(3)-eqd(1)*(1.d0+v(im))+1.d0
         elseif (isolve.eq.1) then
            a1 = -eqd(3)-eqd(1)*(1.d0+v(im))+1.d0
            if (a1.eq.0.d0) a1 = eqd(2)-eqd(1)*v(im)
         elseif (isolve.eq.2) then
            a1 = eqd(2)-eqd(1)*v(im)
            a1 = 0.5d0*(a1-eqd(3)-eqd(1)*(1.d0+v(im)))+0.5d0
         endif
         if (a1.eq.0.d0) then
            iiter = iitermax
            goto 300
         endif
         dv = eqd(1)*w(im)-eqd(4)
         w(i) = -dv/a1
         v(i) = eqd(3)/a1
      enddo

c..   surface
      i = nmod1
      cdi = cdm(i)
      x1 = f1(i)*vv(i)
      x3 = f3(i)*cdi
      eqd(4) = x(i)-xa(i)+x3*(x(i)-x(i-1))-x1*x(i)
      eqd(1) = -x3
      eqd(2) = 1.d0+x3-x1
      eqd(3) = 0.d0
      a1 = eqd(2)-eqd(1)*v(i-1)
      if (a1.eq.0.d0) a1 = -eqd(3)-eqd(1)*(1.d0+v(i-1))+1.d0
      if (a1.eq.0.d0) then
         write (*,'("a1 = 0 in diffsolve at surface #",i4," for species"
     &        ,i3,", nshell =",i4," ["i2"]")') i,ielem,NSHELL,crz(i)
         write (*,'("shell",5x,"X",8x,"Xold",8x,"VV",9x,"cd",9x,"f1",
     &        9x,"f2",9x,"f3")')
         error = 22
         stop
      endif
      dx = (eqd(1)*w(i-1)-eqd(4))/a1

      x(i) = x(i)+alpha*dx
      if (x(i).lt.0.d0) x(i) = 0.5d0*(x(i)-alpha*dx)

      if (x(i).gt.10.d0*tol) then
         corm = dx/x(i)
      else
         corm = dx
      endif

      do i = nmod1-1,1,-1
         dx = -w(i)-v(i)*dx
         x(i) = x(i)+alpha*dx
         if (x(i).lt.0.d0) x(i) = 0.5d0*(x(i)-alpha*dx)

         if (x(i).gt.10.d0*tol) then
            if (abs(dx/x(i)).gt.abs(corm)) then
               icor = i
               corm = dx/x(i)
            endif
         else
            if (abs(dx).gt.abs(corm)) corm = dx
         endif
      enddo
      if (abs(corm).le.cormin) then
         cormin = abs(corm)
         icormin = icor
      endif

      if (indiffus) then
         x(nmod) = x(nmod1)
         do i = 1,nmod1
            abond(i,ielem) = x(i)
         enddo
      endif

 300  if (iiter.ge.iitermax) then
c.. change convergence acceleration
         if (isolve.eq.1) then
            alpha = 0.2d0
            iitermax = 2000
            isolve = 0
            goto 100
         endif
c.. change resolution of equation ? don't remember the origin
         if (isolve.eq.0) then
            isolve = 2
            goto 100
         endif
c.. no convergence
         write (nout,'(3x,"Convergence problem species #",i2)') ielem
         error = 23
         return
      endif

      if (crz(icor).ne.4.and.(abs(corm).gt.tol.or.iiter.lt.iitermin))
     &     goto 200

      if (.not.indiffus) then
         x(nmod) = x(nmod1)
         do i = 1,nmod1
            abond(i,ielem) = x(i)
         enddo
      endif


      return
      end
      SUBROUTINE diffusion (icall,error)

************************************************************************
*     Compute diffusion coefficients and chemical mixing
*
*       The different types of mixing are :
* lmicro = 2   : microscopic diffusion, Chapman & Cowling
* lmicro = 3   : microscopic diffusion, Montmerle&Michaud + Paquette 86 (ionidation partielle)
* lmicro = 4   : atomic diffusion, Thoul et al. 94 + Paquette 86 (ionidation partielle)
* lmicro = 5   : microscopic diffusion, Montmerle&Michaud + Paquette 86 (ionidation totale)     
* lmicro = 6   : atomic diffusion, Thoul et al. 94 + Paquette 86 (ionidation totale)    
* idiffty = 4   : mixing recipe used in Charbonnel ApJ 1995
* idiffty = 8   : chemical mixing + AM transport Talon 1997
* idiffty = 9   : AM transport only, Talon 1997
* idiffty = 10  : (8)+(2)
* idiffty = 11  : (3) + chemical mixing +  AM transport Maeder,Zahn 1998
* idiffty = 13  : (11) + non stationnary terms in AM transport eqs
*                 more complete treatment, recommended for giants
* idiffty = 14  : Omega evolves as a results of structural changes. Angular
*     momentum assumed to be preserved in each shell
* idiffty = 15  : (13) + addition of an additional viscosity      
* idiffty = 17  : Use simplified expression for Ur (case of massive stars models)
* lthal = 01    : Thermohaline mixing
* ltach = 41-44 : Tachocline mixing (see Tachocline diffusc_ondes.f)
* lover = 23-32 : diffusive convective overshoot (Herwig)
*    lover = 23 : below convective envelope
*    lover = 24 : above core
*    lover = 25 : below pulse
*    lover = 26 : above pulse
*    lover = 27 : above core and below enveloppe (=23+24)
*    lover = 28 : below and above pulse (=25+26)
*    lover = 29 : treat all cases (=27+28)
*    lover = 30 : overshoot below all convective zones
*    lover = 31 : overshoot above all convective zones
*     lover = 32 : overshoot everywhere (=30+31)
C Modif TD 2019-2020      
*    lover = 33 : overshoot below EC : Baraffe et al 2017 penetrative convection (see diffbaraffe.f)
*    lover = 34 : overshoot below a convective zone : Augustson & Mathis 2019 rotating convection model - Baraffe and Pratt type diffusion (see diffkyle.f)
*    lover = 35 : overshoot above a convection zone : Augustson & Mathis 2019 rotating convection model - Baraffe and Pratt type diffusion (see diffkyle.f)
*    lover = 36 : overshoot everywhere (=34+35) : Augustson & Mathis 2019 rotating convection model (see diffkyle.f)
*    lover = 37 : overshoot below a convective zone : Korre et al. 2019 type diffusion (see diffkyle.f) 2
*    lover = 38 : overshoot above a convective zone : Korre et al. 2019 (see diffkyle.f) 2
*    lover = 39 : overshoot everywhere (=37+38) : Korre et al. 2019 (see diffkyle.f) 2
*    lover = 70 : overshoot below a convective zone : Augustson & Mathis 2019 rotating convection model - Baraffe and Pratt type diffusion (see diffkyle.f)   + Dturbul T6.4xx
*    lover = 71 : overshoot below a convective zone : Augustson & Mathis 2019 rotating convection model - Baraffe and Pratt type diffusion (see diffkyle.f)  + Dturbul PM10000
*    lover = 72 : overshoot below a convective zone : Korre et al. 2019 (see diffkyle.f) 2 + Dturbul T6.4xx
*    lover = 73 : overshoot below a convective zone : Korre et al. 2019 (see diffkyle.f) 2 + Dturbul PM10000
*
*    lover = 60 : compute and add turbulence coefficient (temperature fixation T.6xx - see diffturbul.f) - only if atomic diffusion active
*    lover = 61 : compute and add turbulence coefficient (BCZ fixation PM10000 - see diffturbul.f) - only if atomic diffusion active     
*     
C Modif CC ondes (5/12/07)
* lover = 41-43 : Internal gravity waves
*    lover = 41 : igw below the convective envelope
*    lover = 42 : igw above the convective core
*     lover = 43 : igw in both cases (=41+42)    
* icall = 0 : call during convergence (mixopt=t), do not treat rotational
*             mixing and do not update vxsp
* icall = 1 : call after convergence
* icall = 2 : for binlist_evol.f
*
* NOTE : if nmixd.gt.0 then diffzc = .FALSE.
*                                                                      *
* $LastChangedDate:: 2016-05-13 11:39:14 +0200 (Fri, 13 May 2016)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 75                                                          $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
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
      include 'evolcom.transp'
      include 'evolcom.var'
      include 'evolcom.igw'
c      include 'evolcom.transpondesexcit'


      integer imerk,icrasht,error
      integer ndt,ndb,ndtenv,imin,imax
      integer i,ielem,j,k,l,ik,kl
      integer nshell,ndbmel,nxmax,jshell        ! le dernier modif TD
      integer icall,kiter,inorm
      integer neqj,ncls,nclc
      integer klenv,klpulse,klcore,klenv0,klpulse0,klcore0,istart,iend
      integer nmodmdix,id,ip
      integer ii,mtu,npasr
      

      integer BCE,EMB

      double precision rho_bcz   ! dturbul TD 09/19
      double precision times,chronos      ! Ajout chronos pour mesure du tps de calcul
      double precision fover,pwi,f1,f2,f3
      double precision abond,summx,xm_sum,anucni
      double precision cd,cd0,dd,dd0,dd1,dmdr
      double precision Dmicro_He,Dmicro_O,Vmicro_He,Vmicro_O ! maj TD Jan.2018
      double precision Dmicro_C,Vmicro_C                     ! maj TD Fev.2018
      double precision dturb,dmic,vdmic,lgdturb,lgdmic
      double precision dlntdr,dlnpdr,dlnxdr,dxdr ! maj thoul
      double precision dri1,dri,a,b,dtdr,dpdr
      double precision xmoment_tot
      double precision vom,vrray,vtheta
      double precision alphasc
      double precision xlambdaold
      double precision dift,dife
      double precision mueinvk,vi
      double precision fLi7,fLi6,fB11,xt,fdecay
      double precision disctime
      double precision abmuliss,abmax,abmin,gradledoux
      double precision THBS,TEMb,F_BS99
      double precision xnurad,xnumol,ddmax
      double precision brunt(nsh)
      double precision dtacho
      double precision Dhold
      double precision rhosh,tsh,muesh ! MODIF THOUL + modif Fev.2019 muesh
      double precision tdiff,tempo,dtkhpms                        ! diffusion PMS
      double precision xdmiche!,Dturbul ! dturbul TD 09/19
      double precision xvdiffhe,xvdiffhenoc
      
!!      double precision thetac
!!      character*11 Dv_prescrp

      logical partialmix,envdef
      logical ipass,istop

      common /envel/ ndtenv
      common /coefdiff/ dturb(nsh),dmic(nsh),vdmic(nsh),lgdturb(nsh)
     &     ,lgdmic(nsh)
      common /coefdiffpaqhe/ xvdiffhe(nsh),xvdiffhenoc(nsh),xdmiche(nsh) ! dturbul TD 09/19
      common /convergence/ imerk,icrasht
      common /dlntdrray/ dlntdr(nsh),dlnpdr(nsh),dlnxdr(nsh,nsp),       ! modif Thoul
     &     dxdr(nsh,nsp)
      common /lambdaold/ xlambdaold(nsh)
      common /moment/ xmoment_tot
      common /oldstar/ vom(nsh),vrray(nsh),vtheta(nsh)
      common /overshoot/ klcore,klenv,klpulse
      common /difcirc/ dift(nsh),dife(nsh)
      common /disclocking/ disctime
      common /calcDh/ Dhold(nsh)
!      common /rotvar/ omega(nsh)    ! add by TD to compute Dkyle
!!      common /Dvprescription / Dv_prescrp,thetac(nsh)


      dimension abond(nsh,nis)
      dimension cd(nsh),dd(nsh),dd0(nsh,3),dd1(nsh)
      dimension Dmicro_He(nsh),Dmicro_O(nsh),Dmicro_C(nsh)                   ! maj TD Jan.2018
      dimension Vmicro_He(nsh),Vmicro_O(nsh),Vmicro_C(nsh) ! maj TD Jan.2018
      dimension dmdr(nsh),f1(nsh),f2(nsh),f3(nsh)
      dimension istop(nis)
      dimension vi(nreac)
      dimension abmuliss(nsh),gradledoux(nsh)
C Modif thermohaline (24/11/2010)
      double precision xKmu,xtau_NL
      double precision xR0_NL,xr_NL,xR0inv_NL

      dimension xR0_NL(nsh),xr_NL(nsh),xtau_NL(nsh)
      dimension xR0inv_NL(nsh)
      dimension  xnurad(nsh),xnumol(nsh)

      common /brunt/ brunt

*____________________
***   Initializations
*--------------------

c..   stop microscopic diffusion after the TO
c      if (microdiffus.and.nphase.gt.2) then
c$$$      if (microdiffus.and.xsp(1,ih1).lt.1.d-5) then
c$$$         microdiffus = .false.
c$$$         write (nout,50)
c$$$         write (90,50)
c$$$  endif
c$$$      if (nphase.ge.3) then
c$$$         lmicro=0               ! coupe diffu phase 3
c$$$         print *, 'no atomic diffusion at phase 3 or later'
c$$$  endif
      
      if (.not.microdiffus.and..not.rotation.and..not.diffusconv
     $     .and..not.difftacho.and..not.thermohaline.and..not.igw.and.
     $     .not.diffover.and..not.turbulence) then
         print *,'no diffusion process active'
         return
      endif
       
      if (idiffty.eq.15) print *,
     $     'Additional viscosity for the transport of angular momentum a
     $ctive'

      ndt = 0
      ndb = 1
      times = time*seci
      call cpu_time(chronos)
      call diffinit (icall)
      if (igw) call diffinitondesexcit
      do k = 1,nmod
         brunt(k)=brunt_o(nmod-k+1)
      enddo

!!      thetac=0.d0

      if (icall.ne.2) then
         do i = 1,nmod
            cd(i) = 0.d0
            Dsc(i) = 0.d0        ! semiconvection
            Dthc(i) = 0.d0       ! thermohaline
            Dherw(i) = 0.d0      ! overshoot
            Dbar(i) = 0.d0       ! overshoot
            Dkyle(i) = 0.d0      ! Convective penetration
            Dhd(i) = 0.d0        ! rotation
            coefDtacho(i) = 0.d0 ! tachocline
            dmic(i) = 0.d0       ! microscopic diffusion
            vdmic(i) = 0.d0      ! microscopic velocity
            Dondeschim(i) = 0.d0 ! Dwaves
            Dturbul(i) = 0.d0    ! ad oc turbulence (TD 09/2019)
         enddo
         if (microdiffus) then
            do i = 1,nmod
               do j = 1,nis
                  abond(i,j) = xsp(i,j)*muizero(i)/anuc(j)
               enddo
            enddo
         else
            do i = 1,nmod
               do j = 1,nis
                  abond(i,j) = xsp(i,j)
               enddo
            enddo
         endif
      endif

      klenv = 0
      klcore = 0
      klpulse = 0
      envdef = .true.
      if (novlim(1,3).eq.1.and.nsconv.gt.1) klcore = 1
      do kl = 1,nsconv
         if (t(novlim(kl,3)).gt.2.d8.and.novlim(kl,3).gt.1.and.
     &        nphase.eq.5.and.klpulse.eq.0) klpulse = kl
         if ((tau(novlim(kl,4)).lt.1.d2.or.t(novlim(kl,4)).lt.1.d6)
     &        .and.mr(novlim(kl,4)).gt.0.5d0.and.klenv.eq.0.and.envdef)
     &        then
            klenv = kl
            envdef = .false.
         endif
         imin = novlim(kl,3)
         imax = novlim(kl,4)
!!         if (kl==klenv) print *,'imin',imin,imax,sconv(imin+1)
         call mix (dm,imin,imax,kl,5,partialmix)
      enddo

c      if (.not.diffzc.and.microdiffus.and.klenv.eq.0) then
c         print *,'out of diffusion.f'
c         return
c      else
c         print *, 'diffusion.f active'
c      endif
      
c..   Initialisations for Rotational mixing (only after convergence)
      if (icall.eq.1.and.(rotation.or.microdiffus.or.igw)) then
         if (klenv.gt.0) then
            ndtenv = novlim(klenv,3)-1
         else
            ndtenv = nmod
         endif
c.. case full coupling
         if (omegaconv.eq.'t') then
            ndb = 1
            ndt = nmod
            if (nphase.lt.2.and.(idiffvr.eq.5.or.idiffvr.ge.8.or.
     $           idiffvr.eq.4).and.novlim(klenv,3)-1.le.1) ndt = 0
         else
c.. case envelope boundary condition
c.. define upper integration boundaries for AM transport
            if (klcore.eq.0) then
               ndb = 1
            else
               ndb = novlim(klcore,4)+1
            endif
            if (klenv.gt.0) then
               ndt = novlim(klenv,3)-1
            else
               ndt = nmod
            endif
         endif


c..   define the number of equations NEQJ to be solved for
c..   the transport of angular momentum, as well as the number
c..   of outer and inner boundary conditions (NCLS and NCLC)

         if (zgradmu.and..not.inits) then
            neqj = 5
            ncls = 3
            nclc = neqj-ncls
         else
            neqj = 4
            ncls = 2
            nclc = neqj-ncls
         endif
      endif

c..   Initialisations for diffusion (Thoul with lmicro = 4)
      if (microdiffus.or.difftacho) then
         do i = 2,nmod1
            dri = dsqrt((r(i+1)**2 + r(i+1)*r(i)+r(i)**2)/3.d0) -
     &           dsqrt((r(i)**2 + r(i)*r(i-1) + r(i-1)**2)/3.d0)
            dtdr = (t(i) - t(i-1))/dri
            dlntdr(i) = 2.d0*dtdr/(t(i)+t(i-1))
            if (lmicro.eq.4) dlntdr(i) = dlntdr(i)*rsun ! modif Thoul

            dpdr = (p(i) - p(i-1))/dri
            dlnpdr(i) = 2.d0*dpdr/(p(i)+p(i-1))
            if (lmicro.eq.4) dlnpdr(i) = dlnpdr(i)*rsun ! modif Thoul
            
            do j = 1,nis        ! modif Thoul (fraction masse)
               dxdr(i,j) = (xsp(i,j) - xsp(i-1,j))/dri
               dlnxdr(i,j) = 2.d0*dxdr(i,j)/(xsp(i,j)+xsp(i-1,j))
               if (lmicro.eq.4) dlnxdr(i,j) = dlnxdr(i,j)*rsun ! modif Thoul
            enddo
         enddo
         dlntdr(1) = dlntdr(2)
      endif

*____________________________________________________
***   Compute diffusion coefficients
***   coefficients that are independent of abundances
*----------------------------------------------------

***   tachocline diffusion
cc      if (difftacho.and.nphase.eq.2.and.klenv.gt.1.and.icall.eq.1) then
cc         ndt = novlim(klenv,3)-1
cc         call tachocline (ndt,1,coefDtacho)
c         stop 'phase 2 ?'
cc      endif

***   Semiconvective diffusion (Langer et al 1985, A&A, 145, 179)

      if (semiconv) then
         alphasc = 4.d-2
         do k = 2,nmod1
            if (crz(k).eq.-3) then
               Dsc(k) = alphasc*khim(k)/(6.d0*rom(k)*cpm(k))*
     &              (abla(k)-abm(k))/(abled(k)-abla(k))
               if (Dconv(k).gt.0.d0) then
                  Dconv(k) = 0.d0
                  write (nout,60) k
               endif
            endif
         enddo
         if (crz(1).eq.3) Dsc(1) = Dsc(2)
      endif

      if (thermohaline) then

         TEMb=0.d0
         THBS=0.d0
         EMB=0

         do k = 2,nmod
            if (xsp(k-1,ih1).lt.1.d-10.and.xsp(k,ih1).ge.1.d-10) then
               THBS = t(k)
               goto 33
            endif
         enddo

 33      TEMb=10**(log10(THBS)-0.262d0)


         do k = 1,nmod
            if (t(k).lt.TEMb) then
	        EMB = k
		goto 331
             endif
        enddo

 331    do k = 1,nmod
           if (Dconv(k).ge.1d5)then
              BCE=k
              goto 221
           endif
        enddo
	

 221    do k = 10,nmod1
           gradledoux(k) = abled(k)-abrad(k)
C Version utilisee pour RGB :
           if ((abmu(k).lt.0.d0).and.(gradledoux(k).gt.0.d0)
     &          .and.(lthal.eq.1)) then
              Dthc(k) = -khim(k)*phiKSm(k)*abmu(k)/(rom(k)*
     &             cpm(k)*deltaKSm(k)*(abm(k)-abla(k)))
              Dthc(k) = Dthc(k)*1.0d3
            endif

******************************************************************************
*... fixer l'efficacite du melange thermohaline jusqu'a une profondeur donnee
*... Test pour comparer avec Boothroyd & Sackmann 99 utilise un parametre fixe
* ...deltalogT = log THBb - log TEMb = 0.262
* ... A commenter pour un calcul normal
******************************************************************************
            if ((abmu(k).lt.0.d0).and.(gradledoux(k).gt.0.d0)
     &           .and.lthal.eq.4) then
               F_BS99=1.0d-4*1.989d33/3.1557d7
               Dthc(k)= F_BS99*(r(BCE)-r(EMB))/(4*pi*rom(EMB)*r(EMB)**2)
            endif

********************************************************************************
*****..... Coefficient pour la diffusion thermohaline donne par Traxler et al
********   ApJ 728:L29 2011
********************************************************************************
	    if ((abmu(k).lt.0.d0).and.(gradledoux(k).gt.0.d0)
     &           .and.(lthal.eq.2)) then
               xnurad(k) = khim(k)*tm(k)/rom(k)/4.4937758d21
               ddmax = 1.5964182d-8*tm(k)**1.5d0/dsqrt(rom(k))
               xnumol(k) = 2.1688529d-15*tm(k)**2.5d0/rom(k)/dlog(ddmax)

               xKt(k)=khim(k)/(rom(k)*cpm(k))
	
               xtau_NL(k)=xnumol(k)/xKt(k)
               xR0_NL(k)=(abm(k)-abla(k))/abmu(k)
               xR0inv_NL(k)=1.d0/xR0_NL(k)
               xr_NL(k)=(xR0inv_NL(k)-1.d0)/((1/xtau_NL(k))-1.d0)
	
	
               Dthc(k) = 101*(xnumol(k)*(xnumol(k)+xnurad(k)))**(0.5d0)*
     &              dexp(-3.6*xr_NL(k))*(1-xr_NL(k))**(1.1d0)
            endif
	
********************************************************************************
*****..... Coefficient pour la diffusion thermohaline donne par Denissenkov 2010
********   ApJ 723:563-579
********************************************************************************

	    if ((abmu(k).lt.0.d0).and.(gradledoux(k).gt.0.d0)
     &           .and.(lthal.eq.3)) then
               Dthc(k)=abmu(k)/(abrad(k)-abm(k))
               Dthc(k)=Dthc(k)*(0.5d0)**(2.d0)*khim(k)/(rom(k)*cpm(k))

               Dthc(k)=19.7d0*Dthc(k)
	    endif
         enddo
	
         if (iter.eq.1) then
            if (lthal.eq.1) then
               write (*,*)'thermohaline Charbonnel & Zahn 07'
            elseif (lthal.eq.2) then
               write(*,*)'thermohaline Traxler et al 2011'
            elseif (lthal.eq.3) then
               write(*,*)'thermohaline Denissenkov 2010'
            elseif (lthal.eq.4) then
               write(*,*)'thermohaline BS 1999'
            endif
         endif
         if (crz(1).eq.2) Dthc(1) = Dthc(2)
      endif

***************************************************
*--   Diffusive overshoot : various prescriptions
***************************************************
      if (diffover.and.nsconv.gt.0) then
         pwi = 1.d0
         iend = 0
         istart = 0

***   diffusive overshoot : Baraffe et al 2017 formalism, low-mass stars (AP Oct.2019)
         if (lover.eq.33) then
            klenv0 = klenv
            istart = novlim(klenv,3)
            if (klenv0.ne.0.and.istart.gt.1) then
               fover = etaturb
               cd0 = Dconv(istart)
               call diffbaraffe (istart,iend,cd0,Dbar,fover)
            endif
            novlim(klenv,7) = iend
         endif

***   diffusive overshoot : Augustson & Mathis 2019 formalism (code KA, adapt. TD Nov. 2019)
c..   Diffusion below a convection zone         
         if (lover.eq.34.or.lover.eq.36.or.lover.eq.37.or.lover.eq.39
     $        .or.(lover.ge.70.and.lover.lt.73)) then
            klenv0 = klenv
            istart = novlim(klenv,3)
            if (klenv0.ne.0.and.istart.gt.1) then
               fover = etaturb
               cd0 = Dconv(istart)
               if (lover.eq.34.or.lover.eq.36.or.lover.eq.70.or.
     $              lover.eq.71) call diffkyle (istart,1,iend,cd0,fover
     $              ,Dkyle,1)
               if (lover.eq.37.or.lover.eq.39.or.lover.eq.72.or.
     $              lover.eq.73)
     $              call diffkyle (istart,1,iend,cd0,fover,Dkyle,2)
            endif
            novlim(klenv,7) = iend
            do i = 1,nmod
               if (Dkyle(i).eq.0.d0) then
                  Dkyle(i) = 1.d0
               endif
            enddo
         endif

c..   Diffusion above a convection zone
         if (lover.eq.35.or.lover.eq.36.or.lover.eq.38.or.lover.eq.39)
     $        then
            klenv0 = klenv
            istart = novlim(klenv,4)
            if (klenv0.ne.0.and.istart.gt.1) then
               fover = etaturb
               cd0 = Dconv(istart)
               if (lover.eq.35.or.lover.eq.36) 
     $              call diffkyle (istart,-1,iend,cd0,fover,Dkyle,1)
               if (lover.eq.38.or.lover.eq.39) 
     $              call diffkyle (istart,-1,iend,cd0,fover,Dkyle,2)
            endif
            novlim(klenv,8) = iend
         endif

***   diffusive overshoot : Herwig formalism, AGB phase
c..   overshooting below ALL convective zones
         if (lover.ge.30.and.lover.le.32) then
            if (lover.eq.30.or.lover.eq.32) then
               do k = 1,nsconv
                  istart = novlim(k,3)
                  if (istart.gt.1) then
                     fover = etaturb
                     cd0 = Dconv(istart)
                     call diffherwig (istart,iend,cd0,Dherw,fover,pwi,1)
                     novlim(k,7) = iend
                     if (icall.eq.1.and.istart.ne.0.and.iend.ne.0) then
                        write (nout,200) istart,iend,t(iend),
     &                       abs(m(istart)-m(iend))/msun,fover
                     endif
                  endif
               enddo
            endif
c..   overshooting above ALL convective zones
            if (lover.eq.31.or.lover.eq.32) then
               do k = 1,nsconv
                  istart = novlim(k,4)
                  fover = etaturb
                  cd0 = Dconv(istart)
                  call diffherwig (istart,iend,cd0,Dherw,fover,pwi,-1)
                  novlim(k,8) = iend
                  if (icall.eq.1.and.istart.ne.0.and.iend.ne.0) then
                     write (nout,200) istart,iend,t(iend),abs(m(istart)-
     &                    m(iend))/msun,fover
                  endif
               enddo
            endif
            goto 10
         endif
         klenv0 = klenv
         klcore0 = klcore
         klpulse0 = klpulse
         if (lover.ne.23.and.lover.ne.27.and.lover.ne.29) klenv0=0
         if (lover.ne.24.and.lover.ne.27.and.lover.ne.29) klcore0=0
         if (lover.eq.23.or.lover.eq.24.or.lover.eq.27) klpulse0=0
         if (lover.eq.22) klenv0 = klenv

c..   overshooting below convective envelope
         if (klenv0.ne.0) then
            istart = novlim(klenv,3)
            fover = etaturb
            cd0 = Dconv(istart)
            if (lover.eq.22) then
               call diffherwig (istart,iend,cd0,Dherw,fover,pwi,2)
            else
               call diffherwig (istart,iend,cd0,Dherw,fover,pwi,1)
            endif
            novlim(klenv,7) = iend
            if (icall.eq.1.and.istart.ne.0.and.iend.ne.0) then
               write (nout,200) istart,iend,t(iend),
     &              abs(m(istart)-m(iend))/msun,fover
            endif
         endif
c..   overshooting above pulse
         if (klpulse0.gt.0) then
            if (lover.eq.26.or.lover.eq.28.or.lover.eq.29) then
               istart = novlim(klpulse,4)
               fover = etaturb
               cd0 = Dconv(istart)
               call diffherwig (istart,iend,cd0,Dherw,fover,pwi,-1)
               novlim(klpulse,8) = iend
               if (icall.eq.1.and.istart.ne.0.and.iend.ne.0) then
                  write (nout,200) istart,iend,t(iend),abs(m(istart)-
     &                 m(iend))/msun,fover
               endif
            endif
c..   overshooting below pulse
            if (lover.eq.25.or.lover.eq.28.or.lover.eq.29) then
               istart = novlim(klpulse,3)
               fover = 1.d-4
               cd0 = Dconv(istart)
               call diffherwig (istart,iend,cd0,Dherw,fover,pwi,1)
               novlim(klpulse,7) = iend
               if (icall.eq.1.and.istart.ne.0.and.iend.ne.0) then
                  write (nout,200) istart,iend,t(iend),abs(m(istart)-
     &                 m(iend))/msun,fover
               endif
            endif
         endif
c..   overshooting above core
         if (klcore0.ne.0) then
            istart = novlim(klcore,4)
            fover = etaturb
            cd0 = Dconv(istart)
            call diffherwig (istart,iend,cd0,Dherw,fover,pwi,-1)
            novlim(klcore,8) = iend
            if (icall.eq.1.and.istart.ne.0.and.iend.ne.0) then
               write (nout,200) istart,iend,t(iend),
     &              abs(m(istart)-m(iend))/msun,fover
            endif
         endif
      endif

 10   if (icall.eq.2) return

*_____________________________________________
***   rotational mixing : angular momentum and
***   chemical species (Maeder & Zahn 1998)
*---------------------------------------------
      
      if (icall.eq.1.and.(rotation.or.difftacho)) then

         vsurf = vomega(nmod)*vr(nmod)*1.d-5
         omega_S = vomega(nmod)

c..   compute angular momentum evolution + diffusion coefficients
c..   imerk = 0 means the structure has converged
c..   routines not called during convergence process
         if (imerk.eq.0) then

            if (zgradmu) then
               xlambdaold(ndb:ndt) = xlambdas(ndb:ndt)
               xlambdaold(1:ndb-1) = 0.d0
               xlambdaold(ndt+1:nmod) = 0.d0
c               call sgsmooth (xlambdaold,ndt,ndb,17,21)
            endif

C Adaptation for IGW in case of fully convective star
C Use omegaconv.eq.true
            if (igwrot.and.ndt.eq.0) then
               ndt = nmod
               ndb = 1
               if (nphase.lt.2.and.(idiffvr.eq.5.or.idiffvr.ge.8.or.
     $              idiffvr.eq.4).and.novlim(klenv,3)-1.le.1) ndtenv = 0
            endif
            if (idiffty.eq.14) then
               ndt = nmod1
               call assymptDV (ndt,xmoment_tot)
            elseif (idiffty.eq.9) then
               call assympt (ndt,neqj,ncls,nclc,times,error)
            else


c..   Solid body rotation during the fully convective PMS phase
               if (nphase.lt.2.and.(idiffvr.eq.5.or.idiffvr.ge.8.or.
     $              idiffvr.eq.4).and.(time/sec.le.disctime)) then
cAP : ajut 10 dec 2012
                  if (idiffvr.eq.5.or.idiffvr.ge.8) idiffvr = 4
                  omega_S = breaktime
                  if (omegaconv.eq.'t') then
                     omegaconv = 'f'
                  endif
                  if (ndt.eq.0.or.(ndtenv.eq.0.and.igwrot)) then
                     omega(1:nmod) = omega_S
                     oms(1:nmod) = 1.d0
                     auxs(1) = 0.d0
                     urs(1) = 0.d0
                     xpsis(1) = 0.d0
                     if (zgradmu) xlambdas(1:nmod) = 0.d0
                     vsurf = omega(nmod)*vr(nmod)*1.d-5
                  else
                     call dif_omega (ndb,ndt,neqj,ncls,nclc,times,error)
                  endif
               else
                  if (idiffvr.ne.idiffvr0) idiffvr = idiffvr0
                  if (omegaconv.ne.omegaconv0) omegaconv = omegaconv0
c                  if (ndt.eq.0) ndt = nmod1
                  if (ndt.eq.0) omegaconv = 's'
                  print *,'ndt=',ndt
                  if (mtini.ge.15.d0.and.xsp(1,ih1).le.5.d-6) omegaconv
     $                 = 'm'
                  call dif_omega (ndb,ndt,neqj,ncls,nclc,times,error)
               endif
            endif

c..   Check that surface velocity has not change too much
c..   (in case of AM transport by IGW)
            do i = ndb+1,ndt-1
               if (omega(i)*omega(i-1).lt.0.d0.and.omega(i)*omega(i+1)
     &              .lt.0.d0) error = 50
            enddo
            if (error.gt.0) return
         endif
          
c.. Tachocline diffusion -> Modif. TD 30/09/2019 - no restriction on phase - specify the choice of the ZC         
c     if (difftacho.and.nphase.eq.2.and.klenv.gt.1) then
         if (difftacho) then
            if (nsconv.eq.1.and.novlim(1,3).ne.1) then
               print *, 'One ZC', 'Dtacho active'
c     ndt = novlim(klenv,3)-1
               ndt = novlim(1,3)-1
               if  (ltach.eq.41) call tachocline (ndt,1,coefDtacho,1)
               if  (ltach.eq.42) call tachocline (ndt,1,coefDtacho,2)
               if  (ltach.eq.43) call tachocline (ndt,1,coefDtacho,3)
               print *,'ndt-1=',ndt-1
            else if (nsconv.gt.1.and.novlim(2,3).ne.1) then
               print *, 'Nbr ZC=',nsconv, 'Dtacho active'
               ndt = novlim(2,3)-1
               if  (ltach.eq.41) call tachocline (ndt,1,coefDtacho,1)
               if  (ltach.eq.42) call tachocline (ndt,1,coefDtacho,2)
               if  (ltach.eq.43) call tachocline (ndt,1,coefDtacho,3)
               print *,'ndt-1=',ndt-1
            endif
            if (igwrot) call coefondes_surf(ndb,ndtenv)
         endif
         do jshell = 1,nmod
            if (coefDtacho(jshell).eq.0.d0) coefDtacho(jshell) = 1.d0 ! set a minimum different of 0
         enddo
      endif
      
c.. diffusion at the center
      Dhd(1) = Dhd(2)
      if (rotation) then
         xnuvv(1) = xnuvv(2)
         V_circ(1) = V_circ(2)
         dift(1) = dift(2)
         dife(1) = dife(2)
      endif

***   Parametric mixing for giant stars

      if (idiffcc.and.ndt.gt.1) then
         if (idiffty.eq.21) then
            call zonedeltaxrgb (ndt,ndb,ndbmel)
            call dturburgb (ndbmel,ndt,diffst,cd)
         endif

         if (zgradmu.and.idiffty.eq.31) then
c            call zonemucarbon (ndb,ibottom)
            call zonemelange (ndt,ndb)
         endif

         if (idiffty.eq.31) call dturbudenissenkov (ndt,ndb,cd)
      endif

      if (microdiffus.and.lmicro.eq.4) then
         do jshell = 1,nmod
            tsh = 0.5*(t(jshell)+t(jshell-1))     ! remise à la couche i, 01/2019
            rhosh = 0.5*(ro(jshell)+ro(jshell-1)) ! remise à la couche i, 01/2019
            muesh = mue(jshell) ! ajout TD AP Fev.2019
c$$$ Note : if convective core, nsconv and novlim(1,4) have to be used to control diffusion            
            call diffthoul (jshell,xsp,anuc,znuc,tsh,rhosh,muesh,abond,
     $           nmod,nsconv,novlim(1,4))              
         enddo
      endif

C Ajout TD 09/2019 - Compute turbulence coefficient     
***   compute the diffusion coefficient of turbulence (Richer et al. 2000; Richard et al.2005)
      if ((lover.ge.60.and.lover.le.61).or.(lover.ge.70
     $     .and.lover.le.73)) then
         print *, 'dturb active'
         if (nphase.eq.1.and.novlim(1,3).eq.1.and.nsconv.eq.1) then
            print *, 'convective core, no ad oc turbulence'
         else if (nphase.eq.1.and.novlim(1,3).ne.1.and.nsconv.eq.1) then
            print *, 'convective zone, ad oc turbulence'
            print *, 'radius bcz',r(novlim(1,3))/rsun
            rho_bcz = ro(novlim(1,3))
            call diffturbul (lover,rho_bcz,xdmiche,Dturbul)
         else if (nphase.eq.1.and.novlim(1,3).eq.1.and.nsconv.ne.1) then
            print *, 'convective zones, ad oc turbulence'
            print *, 'radius bcz',r(novlim(2,3))/rsun
            rho_bcz = ro(novlim(2,3))
            call diffturbul (lover,rho_bcz,xdmiche,Dturbul)
         else if (nphase.gt.1.and.nsconv.eq.1) then
            rho_bcz = ro(novlim(1,3))
            print *, 'phase',nphase
            print *, 'convective zone','radius bcz',r(novlim(1,3))/rsun
            call diffturbul (lover,rho_bcz,xdmiche,Dturbul)
         else if (nphase.gt.1.and.nsconv.gt.1) then
            rho_bcz = ro(novlim(2,3))
            print *, 'more than 1 convective zone'
            print *,'radius bcz',r(novlim(2,3))
     $           /rsun
            call diffturbul (lover,rho_bcz,xdmiche,Dturbul)
         endif
         do jshell = 1,nmod
            if (r(jshell)/rsun.le.0.1) Dturbul(jshell) = 1.d0 ! set a limit for some depth
         enddo
      endif
      
C     Fin ajout TD 09/2019      

***   compute total diffusion coefficient for chemical species

 20   forall (i = 1:nmod) dturb(i) = cd(i)+Dconv(i)+Dsc(i)+Dherw(i)
     &     +coefDtacho(i)+Dhd(i)+Dthc(i)+Dturbul(i)+Dbar(i)+Dkyle(i)

      if (nmixd.gt.0.and.icall.eq.1) return

      if (icall.eq.1) write (nout,500)
c$$$      if (omegaconv.eq.'s'.and.hydrorot) goto 45    ! commenter 17/04/2020

c      if (.not.diffzc.and.microdiffus.and.klenv.eq.0) return   ! commenter
      if (.not.diffzc.and.rotation.and.novlim(nsconv,3)-1.le.1.and.
     &     t(1).lt.1.8d7) goto 45
c     &     nphase.eq.1) goto 45

*______________________________________________
***   Computation of the new abundance profiles
***   Solve the diffusion equation
*----------------------------------------------

c..   define matrix elements for diffusion   
      ipass = .false.
      kiter = 1
      dmdr(1) = 0.d0
      f1(1) = dtn/dm(1)/dble(kiter)
      f2(1) = 2.d0*f1(1)/(dm(1)+dm(2))
      f3(1) = 0.d0
      do i = 2,nmod1
         dmdr(i) = pim4*r(i)*r(i)*rom(i)
         f1(i) = dtn/dm(i)/dble(kiter)
         f2(i) = 2.d0*f1(i)/(dm(i)+dm(i+1))
         f3(i) = 2.d0*f1(i)/(dm(i)+dm(i-1))
         dd(i) = dturb(i)*dmdr(i)*dmdr(i)
         dd1(i) = dd(i)
      enddo
      dd(1) = dd(2)
      dd1(1) = dd1(2)
      
c..   No diffusion for NEUTRONS and for instable isotopes if
c..   tconv >> tdecay & dtn >> tconv (mixing efficient if dtn > dtconv)
      istop(1) = .true.
      istop(2:nis) = .false.

*________________________________________
***   short lived elements do not diffuse
*----------------------------------------

c.. during core H and He burning : instable elements do not diffuse
c..   grid adaptation
      
      if (nphase.le.2.or.nphase.eq.4) goto 30
      
      if (nuclopt.ne.'i'.and.nuclopt.ne.'j') then
         do kl = 1,nsconv
            imin = novlim(kl,3)
            do i = 1,ndecay
               j = kdecay(i)
               if (.not.istop(j)) then
                  if (tdecay(j).lt.1.d-2*turnover(kl)) then
                     if (icall.eq.1) then
                        write (nout,300) elem(j),imin,abond(imin,j),kl
                        if (no.eq.1) write (90,300) elem(j),imin,
     &                       abond(imin,j),kl
                     endif
                     istop(j) = .true.
                  endif
               endif
            enddo
         enddo
      endif
      
C Modif CC Thermohaline (28/03/07) -->
      if (thermohaline.and.nphase.eq.3) then
          istop(ibe7) = .true.
          write(*,*) 'istop(7Be) = ',istop(ibe7)
      endif
C Modif CC <--

*__________________________________________
***   special treatment of LiBeB during HBB
*------------------------------------------

c      if (hbb.and.nuclopt.ne.'j') then
      if ((hbb.and.nuclopt.ne.'j').or.(thermohaline.and.nphase.eq.3))
     &     then
         imin = novlim(klenv,3)
         imax = novlim(klenv,4)
         do k = imin,imax
            mueinvk = mueinv(k)
            call vit (t(k),ro(k),mueinvk,k,vi,1)
            xt = turnover(klenv)*vxsp(k,ih1)
            fLi7 = vi(ili7pa)*xt
            fB11 = vi(ib11pa)*xt
            fLi6 = vi(ili6pa)*xt
c.. dd0 = dd*tnuc/(tnuc+tconv)
            dd0(k,1) = dd(k)/(1.d0+fLi7)
            dd0(k,2) = dd(k)/(1.d0+fB11)
            dd0(k,3) = dd(k)/(1.d0+fLi6)
         enddo
         istop(ili7) = .true.
         istop(ib11) = .true.
         istop(ili6) = .true.
         if (icall.eq.1) then
            write (nout,350) imin,imax
            if (no.eq.1) write (90,350) imin,imax
         endif
      endif
        
c..   compute diffusion
 30   nshell = nmod
      do ik = 1,kiter
c..   loop starts at index 2 because neutrons do not diffuse !!
         if (microdiffus) then
            ielem = 3
         else
            ielem = 2
         endif
         do kl = ielem,nis
            if (istop(kl)) then
               i = 0
               if (kl.eq.ili7) i = 1
               if (kl.eq.ib11) i = 2
               if (kl.eq.ili6) i = 3
               if (i.ne.0) then
                  do k = novlim(klenv,3),novlim(klenv,4)
                     dd1(k) = dd0(k,i)
                  enddo
               else
                  fdecay = tdecay(kl)/(tdecay(kl)+turnover(klenv))
                  do k = novlim(klenv,3),novlim(klenv,4)
                     dd1(k) = dd(k)*fdecay
                  enddo
               endif
               
C Modif CC Thermohaline (28/03/07)
               if (kl.ne.ibe7)
     &          call diffsolve(abond,ratneh,anuc,znuc,dd1,f1,f2,f3,dmdr,
     &              kl,nshell,icall,error)
            else
               call diffsolve(abond,ratneh,anuc,znuc,dd,f1,f2,f3,dmdr,
     &              kl,nshell,icall,error)
            endif
            if (error.gt.0) return
         enddo

c..   normalize and update abundances for nuclear calculation and
c..   for binary output
         summx = 0.d0
         nxmax = 1         
         do j = 1,nmod
            xm_sum = 0.d0
            inorm = 2
c..   loop starts at index 2 because neutrons did not diffuse !!
            if (microdiffus) then
               anucni = anuc(inorm)
               do k = 2,nis
                  anucni = anucni+(anuc(k)-anuc(inorm))*abond(j,k)
               enddo
               do k = 2,nis
                  xsp(j,k) = abond(j,k)*anuc(k)/anucni
                  xm_sum = xm_sum+xsp(j,k)
                  if (xsp(j,k).gt.xsp(j,inorm)) inorm = k ! changement d elem ref si H pas le plus abondant
               enddo
               xsp(j,inorm) = xsp(j,inorm)+(1.d0-xm_sum) ! element ref (H en general) absorbe la difference pour retomber à 1
               if (icall.eq.1) vxsp(j,inorm) = xsp(j,inorm)
               abond(j,inorm) = xsp(j,inorm)*anucni/anuc(inorm)  !modif +
               summx = abs(1.d0-xm_sum)
            else
               do k = 2,nis
                  xsp(j,k) = abond(j,k)
                  xm_sum = xm_sum+xsp(j,k)
                  if (xsp(j,k).gt.xsp(j,inorm)) inorm = k
               enddo
               if (abs(1.d0-xm_sum).gt.summx) then
                  summx = abs(1.d0-xm_sum)
                  nxmax = j
               endif
               xsp(j,inorm) = xsp(j,inorm)+(1.d0-xm_sum)
               if (icall.eq.1)  vxsp(j,inorm) = xsp(j,inorm)
               abond(j,inorm) = xsp(j,inorm)
               if (summx.gt.1.d-8) then
                  if (.not.ipass) then
                     write (nout,400) nxmax,summx,int(diffst)
                     write (90,400) nxmax,summx,int(diffst)
                     ipass = .true.
                     kiter = int(diffst)
                     do i = 1,nmod1
                        f1(i) = f1(i)/diffst
                        f2(i) = f2(i)/diffst
                        f3(i) = f3(i)/diffst
                        do l = 1,nis
c                           abond(i,l) = xsp(i,l)*muizero(i)/anuc(l)
                           abond(i,l) = xsp(i,l)
                        enddo
                     enddo
                     goto 30
                  else
                     error = 17
                     return
                  endif
               endif
            endif
         enddo
      enddo
      if (icall.eq.1.and.summx.gt.1.d-10.and..not.microdiffus) then
         write (nout,100) nxmax,kiter,summx
         write (90,100) nxmax,kiter,summx
      endif
      
c..   if mixing done before nucleosynthesis, update abundances
      if (idiffnuc.eq.1.and.icall.eq.1) then
         vxsp(1:nmod,1:nis) = xsp(1:nmod,1:nis)
      endif


*______________________________________________________
***   Storage of the angular momentum transport results
*------------------------------------------------------

 45   if (inits) inits = .false.
c      if (icall.eq.1.and.rotation.and.ndt.gt.1) then
!!      print *,'icall =',icall,rotation,omega_S,idiffvr

      if (icall.eq.1.and.rotation) then
         call cons_l (oms,xmoment_tot,2,ndb,ndt)
         vom(1:nmod) = vomega(1:nmod)
         vomega(1:nmod) = omega(1:nmod)
         xmom_tots = xmoment_tot
         omega_S = omega(nmod)
      endif


c..   During H burning, if the convective "tongue" is close to the core,
c..   homogeneize between both convective zones
c      if (rotation.and.nphase.eq.2.and.totm.gt.9.d0.and.ndb.eq.
c     &     novlim(2,4)+1) then
c         imin = novlim(1,4)
c         write (nout,*) 'Warning : merging convective zone to the core'
c         write (90,*) 'Warning : merging convective zone to the core'
c         do j = 1,nsp
c            do i = imin,novlim(2,4)
c               xsp(i,j) = xsp(i,imin-1)
c               ysp(i,j) = xsp(i,j)/anuc(j)
c            enddo
c         enddo
c     endif
     

 50   format (' Microscopic diffusion stopped after TO')
 60   format (' Problem ? semiconvective shell ',i4,' has changed !!')
 100  format (' diffusion : MAX CORRECTION [',i4,'], kiter = ',i2,
     &     ', |dX| = ',1pe9.3)
 200  format (1x,'Overshoot from [',i4,'] --> [',i4,', T = ',1pe9.3,']',
     &     ', dm = ',1pe8.2,', fover = ',0pf8.6)
 300  format ('  partial diffusion for ',a5,' : X[',i4,'] = ',1pe9.3,
     &     ', CZ # ',i2)
 350  format ('  HBB : partial mixing for Li6, Li7 and B11 in the ',
     &     'envelope [',i4,',',i4,']')
 400  format (' Maximum correction dX[',i4,'] = ',1pe10.3,
     &     ' > 1.d-10 --> diffusion time step reduced by ',i3)
 500  format (1x,'Processing chemical transport')

      return
      end


************************************************************************

      BLOCK DATA evodat

*                                                                      *
* $LastChangedDate:: 2016-05-30 11:41:18 +0200 (Lun, 30 mai 2016)    $ *
*     $Author:: amard                                                $ *
* Note TD: xspref(3,i) -> AGSS09 - xspref(5,i) -> AY18                 *
* $Rev:: 85                                                          $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.data'

      integer i,j

c..   parameters for the de Jager mass loss prescription
      data ((amlos(i,j),i = 1,6),j = 1,5)/6.34916d0,3.41678d0,
     &   -1.08683d0,0.13095d0,0.22427d0,0.11968d0,-5.0424d0,0.15629d0,
     &   0.41952d0,-0.09825d0,0.46591d0,0.d0,-0.83426d0,2.96244d0,
     &   -1.37272d0,0.13025d0,0.d0,0.d0,-1.13925d0,0.33659d0,-1.07493d0,
     &   0.d0,0.d0,0.d0,-0.12202d0,0.57576d0,0.d0,0.d0,0.d0,0.d0 / 

!!!!!!!!!
!
!Grevesse & Noels (1993)      
!
!Classical (no Dmic, no rotation, grey atm)
!!!!!!!!!
      
c..   solar abundances : Grevesse & Noels (1993)
c..   WARNING : the sum of xspsol must = 1, that's why X(H) has many digits
      data (xspref(1,i),i = 1,nis+6) /1.d-50,0.706437027647287d+00,
     & 4.80396d-05,2.93896d-05,2.76114d-01,6.48936d-10,9.33982d-09,
     & 1.00000d-50,1.00000d-50,1.67237d-10,1.06576d-09,4.72266d-09,
     & 2.97490d-03,3.58090d-05,1.00000d-50,1.00000d-50,9.19646d-04,
     & 3.63188d-06,1.00000d-50,8.35937d-03,3.38761d-06,1.88884d-05,
     & 1.00000d-50,4.87353d-07,1.00000d-50,1.57939d-03,4.02677d-06,
     & 1.27042d-04,1.00000d-50,1.00000d-50,3.47391d-05,1.00000d-50,
     & 5.09446d-04,1.00000d-50,6.69598d-05,7.67973d-05,1.00000d-50,
     & 1.00000d-50,1.00000d-50,5.62930d-05,6.47306d-04,3.39575d-05,
     & 2.33187d-05,6.17237d-06,3.48358d-04,2.83577d-06,1.64250d-05,
     & 1.00000d-50,8.25353d-08,5.92704d-06,1.00000d-50,2.00021d-06,
     & 1.51058d-03,3.74987d-06,6.23352d-05,2.89829d-06,1.76933d-05,
     & 1.32341d-05,1.28476d-03,7.42195d-05 /

!!!!!!!!!
!
!AGS05 (2005)
!      
!Classical (no Dmic, no rotation, grey atm)
!!!!!!!!!      

C...	solar abundances AGS05 calibrated for solar model
C...	without diffusion and alpha_MLT = 1.6046
      data (xspref(2,i),i = 1,nis+6) / 1.D-50,0.726189670660267D+00,
     & 5.01012D-05,2.67030D-05,2.61517D-01,5.89481D-10,8.48412D-09,
     & 1.00000D-50,1.00000D-50,1.59074D-10,8.23997D-10,3.65136D-09,
     & 2.14655D-03,2.58381D-05,1.00000D-50,1.00000D-50,6.19283D-04,
     & 2.44568D-06,1.00000D-50,5.37579D-03,2.17852D-06,1.21468D-05,
     & 1.00000D-50,3.76802D-07,1.00000D-50,9.47888D-04,2.41672D-06,
     & 7.62457D-05,1.00000D-50,1.00000D-50,3.15563D-05,1.00000D-50,
     & 4.73550D-04,1.00000D-50,6.22419D-05,7.13862D-05,1.00000D-50,
     & 1.00000D-50,1.00000D-50,5.35454D-05,6.15712D-04,3.23001D-05,
     & 2.21805D-05,5.73747D-06,3.23814D-04,2.63596D-06,1.52678D-05,
     & 1.00000D-50,7.67200D-08,3.31976D-06,1.00000D-50,1.12032D-06,
     & 1.36771D-03,3.74987D-06,6.23352D-05,2.89829D-06,1.76933D-05,
     & 1.32341D-05,1.28476D-03,7.42195D-05 /

!!!!!!!!!
!
! Asplund 2009      
!
! Classical (no Dmic, no rotation, grey atm)      
!!!!!!!!!      
      
c$$$c      solar Abundances Asplund et al. 2009 calibrated for solar model
c$$$c      without diffusion and alpha_MLT = 1.7020  data refdolar /'AGSS09'/
c$$$c	correct isotopic ratio and Y = 0.26857  and Z = 0.01340 dr=1e-05 dl=0
c$$$c      from zini_4_1 in ~/STAREVOL/MODELS/Sun_std_0.0134460.2485001.7020
c$$$c      on regor2 GENEVA ... Nadege
c$$$      data (xspref(3,i),i = 1,nis+6) / 1.D-50,0.7179541322222622D+00,
c$$$     & 5.00478E-05,2.66687D-05,2.68569D-01,7.09808D-10,9.68678D-09,
c$$$     & 1.00000D-50,1.00000D-50,1.58346D-10,7.31503D-10,3.23882D-09,
c$$$     & 2.34268D-03,2.83883D-05,1.00000D-50,1.00000D-50,6.92619D-04,
c$$$     & 1.70329D-06,1.00000D-50,5.73387D-03,2.31446D-06,1.29319D-05,
c$$$     & 1.00000D-50,5.05960D-07,1.00000D-50,1.16040D-03,2.92077D-06,
c$$$     & 9.39517D-05,1.00000D-50,1.00000D-50,2.93150D-05,1.00000D-50,
c$$$     & 5.53536D-04,1.00000D-50,7.29966D-05,8.35840D-05,1.00000D-50,
c$$$     & 1.00000D-50,1.00000D-50,5.58118D-05,6.12902D-04,3.22332D-05,
c$$$     & 2.19811D-05,5.84419D-06,2.93708D-04,2.42488D-06,1.41026D-05,
c$$$     & 1.00000D-50,6.96137D-08,6.15156D-06,1.00000D-50,2.07844D-06,
c$$$     & 1.54110D-03,8.71029D-05,6.46799D-05,3.14384D-06,1.67516D-05,
c$$$  & 1.09110D-05,1.29673D-03,7.18769D-05 /
      

c      solar Abundances Asplund et al. 2009 calibrated for solar model
c      without diffusion and alpha_MLT = 1.6267  data refdolar /'AGSS09'/
c correct isotopic ratio and Y = 0.268885  and Z = 0.013446 dr=1e-05 dl=0
c      from zini.out_2_1 in ~/MODELS/Sun_std_0.0134460.2700001.7500
c      on port-palacios2 MONTPELLIER ... Ana
c$$$      data (xspref(3,i),i = 1,nis+6) / 1.D-50,0.717591601961690E+00,
c$$$     & 4.88113D-05,2.85992D-05,2.68885D-01,7.12237D-10,9.71992D-09,
c$$$     & 1.00000D-50,1.00000D-50,1.58888D-10,7.34006D-10,3.24990D-09,
c$$$     & 2.35070D-03,2.84854D-05,1.00000D-50,1.00000D-50,6.94989D-04,
c$$$     & 1.70912D-06,1.00000D-50,5.75349D-03,2.32238D-06,1.29761D-05,
c$$$     & 1.00000D-50,5.07691D-07,1.00000D-50,1.16437D-03,2.93076D-06,
c$$$     & 9.42732D-05,1.00000D-50,1.00000D-50,2.94153D-05,1.00000D-50,
c$$$     & 5.55430D-04,1.00000D-50,7.32464D-05,8.38700D-05,1.00000D-50,
c$$$     & 1.00000D-50,1.00000D-50,5.60028D-05,6.14999D-04,3.23435D-05,
c$$$     & 2.20563D-05,5.86419D-06,2.94713D-04,2.43318D-06,1.41509D-05,
c$$$     & 1.00000D-50,6.98519D-08,6.17261D-06,1.00000D-50,2.08555D-06,
c$$$     & 1.54637D-03,8.74008D-05,6.49011D-05,3.15459D-06,1.68090D-05,
c$$$  & 1.09483D-05,1.30116D-03,7.21227D-05 /

      
! Solar abundances Asplund et al. 2009 calibrated for solar model without diffusion
! or rotation but grey atmosphere ;
! alpha_MLT = 1.6602.
! Isotopic ratio Y = 0.268986, Z = 0.013446 ; calib : dr = 1.44e-5, dl = 8.05e-6
! from zini.out_4_1 in ~/STAREVOL/MODELS/Sun_std_0.0134460.2688851.6267_greyatm/ on port-dumont. Thibaut (11/07/2019). 
c$$$      data (xspref(3,i),i = 1,nis+6) / 1.D-50,0.717491000000000D+00,
c$$$     & 4.87946D-05,2.86253D-05,2.68986D-01,7.12237D-10,9.71992D-09,
c$$$     & 1.00000D-50,1.00000D-50,1.58888D-10,7.34006D-10,3.24990D-09, 
c$$$     & 2.35070D-03,2.84854D-05,1.00000D-50,1.00000D-50,6.94989D-04,
c$$$     & 1.70912D-06,1.00000D-50,5.75349D-03,2.32238D-06,1.29761D-05,
c$$$     & 1.00000D-50,5.07691D-07,1.00000D-50,1.16437D-03,2.93076D-06, 
c$$$     & 9.42732D-05,1.00000D-50,1.00000D-50,2.94153D-05,1.00000D-50,
c$$$     & 5.55430D-04,1.00000D-50,7.32464D-05,8.38700D-05,1.00000D-50,
c$$$     & 1.00000D-50,1.00000D-50,5.60028D-05,6.14999D-04,3.23435D-05, 
c$$$     & 2.20563D-05,5.86419D-06,2.94713D-04,2.43318D-06,1.41509D-05,
c$$$     & 1.00000D-50,6.98519D-08,6.17261D-06,1.00000D-50,2.08555D-06,
c$$$     & 1.54637D-03,1.00000D-50,1.00000D-50 /      

C...	solar abundances Asplund et al. 2009 calibrated for solar model
C...	without diffusion and alpha_MLT = 1.6304
C...  from zini.out_3_1 iin ~/CALCULS/MODELS/Sun_std_0.0134460.2700001.7500
C...  on kamina (@LUPM)
      data (xspref(3,i),i = 1,nis+6) / 1.E-50,0.7177592127594249E+00,
     & 5.00427D-05,2.66718D-05,2.68713D-01,6.02506D-10,8.67158D-09,
     & 1.00000D-50,1.00000D-50,1.58889D-10,7.33530D-10,3.25047D-09,
     & 2.35089D-03,2.82978D-05,1.00000D-50,1.00000D-50,6.94034D-04,
     & 2.74089D-06,1.00000D-50,5.75352D-03,2.33160D-06,1.30004D-05,
     & 1.00000D-50,5.07696D-07,1.00000D-50,1.16479D-03,2.96973D-06,
     & 9.36931D-05,1.00000D-50,1.00000D-50,2.94157D-05,1.00000D-50,
     & 5.55722D-04,1.00000D-50,7.30423D-05,8.37733D-05,1.00000D-50,
     & 1.00000D-50,1.00000D-50,5.60034D-05,6.14993D-04,3.22623D-05,
     & 2.21545D-05,5.86424D-06,2.94976D-04,2.40122D-06,1.39081D-05,
     & 1.00000D-50,6.98876D-08,6.17445D-06,1.00000D-50,2.08371D-06,
     & 1.54644D-03,8.74069D-05,6.49057D-05,3.15482D-06,1.68102D-05,
     & 1.09491D-05,1.30125D-03,7.21277D-05 /

      
!!!!!!!!!
!
! Grid version
!      
!!!!!!!!!       

C Warning: Solar composition changed in GRID version
c..   WARNING : the sum of xspsol must = 1, that's why X(H) has many digits
c..   solar abundaces: Asplund 2009 with enhanced Ne (GVA GRID)
      data refsolar /'GRID'/
      data (xspref(4,i),i = 1,nis+6) /1.00000D-50,7.199079447D-1,
     & 4.89617D-05,4.30936D-05,2.66000D-01,6.34807D-10,8.99628D-09,
     & 1.00000D-50,1.00000D-50,1.68736D-10,8.09058D-10,3.94194D-09,
     & 2.26531D-03,3.63117D-05,1.00000D-50,1.00000D-50,6.56299D-04,
     & 2.34183D-06,1.00000D-50,5.69064D-03,3.82024D-06,1.28412D-05,
     & 1.00000D-50,5.38319D-07,1.00000D-50,1.78666D-03,5.69874D-06,
     & 2.39665D-04,1.00000D-50,1.00000D-50,2.65373D-05,1.00000D-50,
     & 4.99158D-04,1.00000D-50,6.69308D-05,7.67465D-05,1.00000D-50,
     & 1.00000D-50,1.00000D-50,4.93616D-05,5.97360D-04,6.65170D-05,
     & 4.81417D-05,5.53754D-06,3.24069D-04,3.48086D-06,1.80401D-05,
     & 1.00000D-50,6.54058D-06,1.00000D-50,1.00000D-50,2.14896D-06,
     & 1.46620D-03,
     & 3.74987d-06,6.23352d-05,2.89829d-06,1.76933d-05,1.32341d-05,
     & 1.28476d-03,7.42195d-05 /


!!!!!!!!!
!
! Asplund 2009      
!
! Classical (no Dmic, no rotation, Phoenix)      
!!!!!!!!!  

! Solar Abundances Asplund et al. 2009 calibrated for solar model without diffusion,
! with atmospheric interpolation in PHOENIX atmosphere table and alpha_MLT = 1.9730,
! Isotopic ratio Y = 0.269076, Z = 0.0134460 dr=1e-5, dl=1e-6
! from zini.out_31_1 in ~/STAREVOL/MODELS/Sun_std_0.0134460.2690621.9730 on port-amard
! Louis
      data (xspref(5,i),i = 1,nis+6) / 1.E-50,0.717400000000000E+00,
     & 4.88113D-05,2.85992D-05,2.69076D-01,7.12237D-10,9.71992D-09,
     & 1.00000D-50,1.00000D-50,1.58888D-10,7.34006D-10,3.24990D-09,
     & 2.35070D-03,2.84854D-05,1.00000D-50,1.00000D-50,6.94989D-04,
     & 1.70912D-06,1.00000D-50,5.75349D-03,2.32238D-06,1.29761D-05,
     & 1.00000D-50,5.07691D-07,1.00000D-50,1.16437D-03,2.93076D-06,
     & 9.42732D-05,1.00000D-50,1.00000D-50,2.94153D-05,1.00000D-50,
     & 5.55430D-04,1.00000D-50,7.32464D-05,8.38700D-05,1.00000D-50,
     & 1.00000D-50,1.00000D-50,5.60028D-05,6.14999D-04,3.23435D-05,
     & 2.20563D-05,5.86419D-06,2.94713D-04,2.43318D-06,1.41509D-05,
     & 1.00000D-50,6.98519D-08,6.17261D-06,1.00000D-50,2.08555D-06,
     & 1.54637D-03,8.74008D-05,6.49011D-05,3.15459D-06,1.68090D-05,
     & 1.09483D-05,1.30116D-03,7.21227D-05 /

!!!!!!!!!
!
! Asplund 2009      
!
! Classical (no Dmic, no rotation, KS66 atm)      
!!!!!!!!!      
      
! Solar abundances Asplund et al.2009 calibrated for solar model without diffusion
! and rotation but with KS66 atmosphere ;
! alpha_MLT = 2.1069.
! Isotopic ratio Y = 0.268558, Z = 0.013446 ; calib : dr = 4.66e-5, dl = 2.86e-5 
! from zini.out_5_1 in ~/STAREVOL/MODELS/Sun_std_0.0134460.2688851.6267  / on port-dumont.       
! Thibaut (22/07/2019).      
      data (xspref(6,i),i = 1,nis+6) / 1.D-50,0.717919000000000D+00,
     & 4.87946D-05,2.86253D-05,2.68558D-01,7.12237D-10,9.71992D-09,
     & 1.00000D-50,1.00000D-50,1.58888D-10,7.34006D-10,3.24990D-09, 
     & 2.35070D-03,2.84854D-05,1.00000D-50,1.00000D-50,6.94989D-04,
     & 1.70912D-06,1.00000D-50,5.75349D-03,2.32238D-06,1.29761D-05,
     & 1.00000D-50,5.07691D-07,1.00000D-50,1.16437D-03,2.93076D-06, 
     & 9.42732D-05,1.00000D-50,1.00000D-50,2.94153D-05,1.00000D-50,
     & 5.55430D-04,1.00000D-50,7.32464D-05,8.38700D-05,1.00000D-50,
     & 1.00000D-50,1.00000D-50,5.60028D-05,6.14999D-04,3.23435D-05, 
     & 2.20563D-05,5.86419D-06,2.94713D-04,2.43318D-06,1.41509D-05,
     & 1.00000D-50,6.98519D-08,6.17261D-06,1.00000D-50,2.08555D-06,
     & 1.54637D-03,8.74008D-05,6.49012D-05,3.15460D-06,1.68090D-05,
     & 1.09483D-05,1.30116D-03,7.21227D-05 /

!!!!!!!!!
!
! Asplund 2009      
!
! Standard (Dmic, no rotation, grey atm)      
!!!!!!!!!      

! Solar abundances Asplund et al. 2009 calibrated for solar model with diffusion
! (microscopic diffusion, M&M + Paquette) and without rotation and atmosphere ;
! alpha_MLT = 1.7936.
! Isotopic ratio Y = 0.273093, Z = 0.01443014 ; calib : dr = 6.27e-6, dl = 1.41e-4, dZsx = 5.18e-5
! from zini.out_24_1 in ~/STAREVOL/MODELS/Sun_diff_0.0147910.2757801.8003/ on port-dumont. 
! Thibaut (24/04/2018).
c$$$      data (xspref(7,i),i = 1,nis+6) / 1.D-50,0.712039000000000D+00,
c$$$     & 4.82344D-05,2.93589D-05,2.73092D-01,7.83474D-10,1.06921D-08,
c$$$     & 1.00000D-50,1.00000D-50,1.74779D-10,8.07421D-10,3.57496D-09, 
c$$$     & 2.58582D-03,3.13346D-05,1.00000D-50,1.00000D-50,7.64502D-04,
c$$$     & 1.88007D-06,1.00000D-50,6.32895D-03,2.55467D-06,1.42741D-05,
c$$$     & 1.00000D-50,5.58470D-07,1.00000D-50,1.28083D-03,3.22390D-06, 
c$$$     & 1.03702D-04,1.00000D-50,1.00000D-50,3.23575D-05,1.00000D-50,
c$$$     & 6.10983D-04,1.00000D-50,8.05725D-05,9.22586D-05,1.00000D-50,
c$$$     & 1.00000D-50,1.00000D-50,6.16041D-05,6.76511D-04,3.55785D-05, 
c$$$     & 2.42623D-05,6.45071D-06,3.24191D-04,2.67654D-06,1.55663D-05,
c$$$     & 1.00000D-50,7.68386D-08,6.78999D-06,1.00000D-50,2.29414D-06,
c$$$     & 1.70105D-03,1.00000D-50,1.00000D-50 /

! Solar abundances Asplund et al.2009 calibrated for solar model with diffusion
! (microscopic diffusion, Thoul + coeff. Paquette + ionisation partielle) and without rotation and atmosphere ;
! alpha_MLT = 1.7727.
! Isotopic ratio Y = 0.273325, Z = 0.014395 ; calib : dr = -1.03e-5, dl = -3.37e-4, dZsx = -8.55e-5
! from zini.out_1_1 in ~/STAREVOL/MODELS/Sun_diff_0.0143280.2728391.7715  / on port-dumont. 
! Thibaut (18/04/2019).
      data (xspref(7,i),i = 1,nis+6) / 1.D-50,0.712270000000000D+00,
     & 4.85023D-05,2.89888D-05,2.73325D-01,7.58928D-10,1.03571D-08,
     & 1.00000D-50,1.00000D-50,1.69303D-10,7.82125D-10,3.46296D-09, 
     & 2.50481D-03,3.03529D-05,1.00000D-50,1.00000D-50,7.40551D-04,
     & 1.82117D-06,1.00000D-50,6.13067D-03,2.47463D-06,1.38269D-05,
     & 1.00000D-50,5.40973D-07,1.00000D-50,1.24070D-03,3.12290D-06, 
     & 1.00453D-04,1.00000D-50,1.00000D-50,3.13438D-05,1.00000D-50,
     & 5.91841D-04,1.00000D-50,7.80482D-05,8.93682D-05,1.00000D-50,
     & 1.00000D-50,1.00000D-50,5.96741D-05,6.55316D-04,3.44638D-05, 
     & 2.35022D-05,6.24861D-06,3.14034D-04,2.59269D-06,1.50786D-05,
     & 1.00000D-50,7.44313D-08,6.57726D-06,1.00000D-50,2.22227D-06,
     & 1.64776D-03,9.31314D-05,7.11668D-05,3.88887D-06,1.88504D-05,
     & 1.04633D-05,1.36934D-03,8.09215D-05 /

!!!!!!!!!
!
! Asplund 2009      
!
! Standard (Dmic, no rotation, KS66 atm)      
!!!!!!!!!      
      
! Solar abundances Asplund et al.2009 calibrated for solar model with atomic diffusion
! and  KS66 atmosphere but without rotation ;
! alpha_MLT = 2.2573.
! Isotopic ratio Y = 0.271557, Z = 0.014140 ; calib : dr = 2.60e-5, dl = 4.31e-5, dZsx = -3.26e-5 
! from zini.out_4_1 in ~/STAREVOL/MODELS/Sun_diff_0.0144050.2728251.7127_(Vgood)  / on port-dumont. 
! Thibaut (28/06/2019).
      data (xspref(8,i),i = 1,nis+6) / 1.D-50,0.713961000000000D+00,
     & 4.86296D-05,2.87814D-05,2.71557D-01,7.63035D-10,1.04132D-08,
     & 1.00000D-50,1.00000D-50,1.70220D-10,7.86357D-10,3.48169D-09, 
     & 2.51836D-03,3.05170D-05,1.00000D-50,1.00000D-50,7.44557D-04,
     & 1.83102D-06,1.00000D-50,6.16384D-03,2.48802D-06,1.39016D-05,
     & 1.00000D-50,5.43900D-07,1.00000D-50,1.24742D-03,3.13979D-06, 
     & 1.00997D-04,1.00000D-50,1.00000D-50,3.15133D-05,1.00000D-50,
     & 5.95044D-04,1.00000D-50,7.84705D-05,8.98518D-05,1.00000D-50,
     & 1.00000D-50,1.00000D-50,5.99970D-05,6.58862D-04,3.46503D-05, 
     & 2.36294D-05,6.28243D-06,3.15732D-04,2.60672D-06,1.51602D-05,
     & 1.00000D-50,7.48339D-08,6.61285D-06,1.00000D-50,2.23430D-06,
     & 1.65666D-03,9.36344D-05,6.95300D-05,3.37959D-06,1.80079D-05,
     & 1.17292D-05,1.39396D-03,7.72666D-05 /

!!!!!!!!!
!
! Asplund 2009      
!
! Non-Standard (Dmic, Rotation, grey atm)      
!!!!!!!!!

!!!!!!!!!
!
! Asplund 2009      
!
! Non-Standard (Dmic, Rotation, Phoenix atm)      
!!!!!!!!!
      
! Solar abundances Asplund et al. 2009 calibrated for solar model with diffusion and rotation,
! (prescrp = Dh from Mathis et al 2016, Dv from Zahn 1992, magnetic braking from Matt 2015)
! with atmospheric interpolation in PHOENIX atmosphere table and alpha_MLT = 2.1954,
! Isotopic ratio Y = 0.273294, Z = 0.01441608, dr = 1e-4, dl = 2e-4
! from zini.out_8_4 in ~/STAREVOL/MODELS/Sun_rot_0.02095595852250.2690002.19862.2kms on port-amard
! Louis.
      data (xspref(10,i),i = 1,nis+6) / 1.E-50,0.712213127849917E+00,
     & 4.86114D-05,2.88088D-05,2.73294D-01,7.63590D-10,1.04207D-08,
     & 1.00000D-50,1.00000D-50,1.70344D-10,7.86929D-10,3.48422D-09,
     & 2.52019D-03,3.05392D-05,1.00000D-50,1.00000D-50,7.45098D-04,
     & 1.83235D-06,1.00000D-50,6.16832D-03,2.48983D-06,1.39117D-05,
     & 1.00000D-50,5.44296D-07,1.00000D-50,1.24832D-03,3.14207D-06,
     & 1.01070D-04,1.00000D-50,1.00000D-50,3.15362D-05,1.00000D-50,
     & 5.95477D-04,1.00000D-50,7.85275D-05,8.99171D-05,1.00000D-50,
     & 1.00000D-50,1.00000D-50,6.00407D-05,6.59341D-04,3.46755D-05,
     & 2.36466D-05,6.28700D-06,3.15962D-04,2.60861D-06,1.51712D-05,
     & 1.00000D-50,7.48883D-08,6.61766D-06,1.00000D-50,2.23592D-06,
     & 1.65786D-03,9.37022D-05,6.95804D-05,3.38203D-06,1.80209D-05,
     & 1.17376D-05,1.39497D-03,7.73226D-05 /      

!!!!!!!!!
!
! Asplund 2009      
!
! Non-Standard (Dmic, Rotation, KS66 atm)      
!!!!!!!!!      

! Solar abundances Asplund et al.2009 calibrated for solar model with atomic diffusion
! and  KS66 atmosphere and rotation (prescription = Dh from Mathis et al 2018,
! Dv from Zahn 1992, magnetic braking from Matt 2015) ; alpha_MLT = 2.1987.
! Isotopic ratio Y = 0.270127, Z = 0.013850 ; calib : dr = 8.76e-5, dl = 5.05e-4, dZsx = 3.13e-4 
! from zini.out_2_1 in ~/STAREVOL/MODELS/Sun_rot_0.0141400.2699712.17462.2kms  / on port-dumont. 
! Thibaut (13/08/2019).
      data (xspref(11,i),i = 1,nis+6) / 1.D-50,0.715655000000000D+00,
     & 4.85657D-05,2.89096D-05,2.70127D-01,7.48990D-10,1.02215D-08,
     & 1.00000D-50,1.00000D-50,1.67087D-10,7.71883D-10,3.41761D-09, 
     & 2.47201D-03,2.99554D-05,1.00000D-50,1.00000D-50,7.30853D-04,
     & 1.79732D-06,1.00000D-50,6.05039D-03,2.44223D-06,1.36458D-05,
     & 1.00000D-50,5.33889D-07,1.00000D-50,1.22446D-03,3.08200D-06, 
     & 9.91380D-05,1.00000D-50,1.00000D-50,3.09333D-05,1.00000D-50,
     & 5.84092D-04,1.00000D-50,7.70262D-05,8.81980D-05,1.00000D-50,
     & 1.00000D-50,1.00000D-50,5.88927D-05,6.46735D-04,3.40125D-05, 
     & 2.31945D-05,6.16679D-06,3.09921D-04,2.55874D-06,1.48811D-05,
     & 1.00000D-50,7.34565D-08,6.49113D-06,1.00000D-50,2.19318D-06,
     & 1.62617D-03,9.19111D-05,6.82504D-05,3.31739D-06,1.76765D-05,
     & 1.15133D-05,1.36831D-03,7.58445D-05 /

!!!!!!!!!
!
! Young 2018      
!
! Classical (no Dmic, no rotation, KS66 atm)      
!!!!!!!!!

! Solar abundances Asplund et al.2009 + Ne enhancement (Young+18) calibrated for solar model without diffusion
! and rotation but with KS66 atmosphere ;
! alpha_MLT = 2.1100.
! Isotopic ratio Y = 0.268507, Z = 0.013446 ; calib : dr = 1.00e-6, dl = 1.00e-6 
! from zini.out_2_1 in ~/STAREVOL/MODELS/Sun_std_0.0134460.2678612.1237  / on port-dumont.       
! Thibaut (17/04/2020). Maj 10/06/2020 with mass loss.     
      data (xspref(12,i),i = 1,nis+6) / 1.D-50,0.717969000000000D+00,
     & 4.88825D-05,2.84877D-05,2.68507D-01,5.80068D-10,8.34864D-09,
     & 1.00000D-50,1.00000D-50,1.52972D-10,7.06213D-10,3.12942D-09, 
     & 2.26334D-03,2.72439D-05,1.00000D-50,1.00000D-50,6.68188D-04,
     & 2.63881D-06,1.00000D-50,5.53926D-03,2.24478D-06,1.25163D-05,
     & 1.00000D-50,4.88789D-07,1.00000D-50,1.58404D-03,4.03864D-06, 
     & 1.27415D-04,1.00000D-50,1.00000D-50,2.83202D-05,1.00000D-50,
     & 5.35027D-04,1.00000D-50,7.03221D-05,8.06536D-05,1.00000D-50,
     & 1.00000D-50,1.00000D-50,5.39177D-05,5.92090D-04,3.10608D-05, 
     & 2.13295D-05,5.64586D-06,2.83991D-04,2.31180D-06,1.33902D-05,
     & 1.00000D-50,6.72850D-08,5.94449D-06,1.00000D-50,2.00610D-06,
     & 1.48884D-03,8.41475D-05,6.24854D-05,3.03717D-06,1.61834D-05,
     & 1.05408D-05,1.25273D-03,6.94381D-05 / 
      
      
!!!!!!!!!
!
! Young 2018      
!
! classical : not calibrated      
!!!!!!!!!
      data (xspref(15,i),i = 1,nis+6) / 1.E-50,0.735583029735419E+00,
     & 2.94239D-05,3.11795D-05,2.50397D-01,6.34481D-10,8.65879D-09,
     & 1.00000D-50,1.00000D-50,1.58812D-10,7.33658D-10,3.24837D-09,
     & 2.34958D-03,2.84720D-05,1.00000D-50,1.00000D-50,6.94661D-04,
     & 1.70831D-06,1.00000D-50,5.75077D-03,2.32129D-06,1.29701D-05,
     & 1.00000D-50,5.07451D-07,1.00000D-50,1.64395D-03,4.13785D-06,
     & 1.33101D-04,1.00000D-50,1.00000D-50,2.94014D-05,1.00000D-50,
     & 5.55167D-04,1.00000D-50,7.32117D-05,8.38303D-05,1.00000D-50,
     & 1.00000D-50,1.00000D-50,5.59763D-05,6.14709D-04,3.23282D-05,
     & 2.20458D-05,5.86141D-06,2.94574D-04,2.43202D-06,1.41441D-05,
     & 1.00000D-50,6.98189D-08,6.16969D-06,1.00000D-50,2.08457D-06,
     & 1.54565D-03,0.00000D+00,0.00000D+00,0.00000D+00,0.00000D+00,
     & 1.05408D-50,1.25273D-50,6.94381D-50 / 

      
      
!!!!!!!!!
!
! Young 2018      
!
! Standard (Dmic, no rotation, KS66 atm)      
!!!!!!!!!       
      

! Solar abundances Asplund et al.2009 + Ne enhancement (Young+18) calibrated for solar model with atomic diffusion
! and  KS66 atmosphere but without rotation ;
! alpha_MLT = 2.2357.
! Isotopic ratio Y = 0.267852, Z = 0.013628 ; calib : dr = 1.87e-3, dl = 2.18e-4, dZsx = 1.73e-3 
! from zini.out_0_1 in ~/STAREVOL/MODELS/Sun_diff_0.0136280.2678522.2357  / on port-dumont. 
! Thibaut (20/04/2020). 
c$$$      data (xspref(13,i),i = 1,nis+6) / 1.D-50,0.718443000000000D+00,
c$$$     & 4.89741D-05,2.83250D-05,2.67852D-01,5.87907D-10,8.46146D-09,
c$$$     & 1.00000D-50,1.00000D-50,1.55039D-10,7.15756D-10,3.17171D-09, 
c$$$     & 2.29393D-03,2.76121D-05,1.00000D-50,1.00000D-50,6.77217D-04,
c$$$     & 2.67447D-06,1.00000D-50,5.61411D-03,2.27511D-06,1.26854D-05,
c$$$     & 1.00000D-50,4.95394D-07,1.00000D-50,1.60545D-03,4.09322D-06, 
c$$$     & 1.29137D-04,1.00000D-50,1.00000D-50,2.87029D-05,1.00000D-50,
c$$$     & 5.42257D-04,1.00000D-50,7.12724D-05,8.17435D-05,1.00000D-50,
c$$$     & 1.00000D-50,1.00000D-50,5.46463D-05,6.00091D-04,3.14805D-05, 
c$$$     & 2.16177D-05,5.72215D-06,2.87829D-04,2.34304D-06,1.35711D-05,
c$$$     & 1.00000D-50,6.81942D-08,6.02482D-06,1.00000D-50,2.03321D-06,
c$$$     & 1.50896D-03,8.52864D-05,6.33311D-05,3.07828D-06,1.64024D-05,
c$$$     &     1.06834D-05,1.26968D-03,7.03779D-05 /

! Isotopic ratio Y = 0.271557, Z = 0.014140, alpha_MLT = 2.2600      
calib : dr = 2.10e-3, dl = 1.3e-3, dZsx = 8.20e-5      
      data (xspref(13,i),i = 1,nis+6) / 1.D-50,0.713961000000000D+00,
     & 4.86296D-05,2.87814D-05,2.71557D-01,7.63035D-10,1.04132D-08,
     & 1.00000D-50,1.00000D-50,1.70220D-10,7.86357D-10,3.48169D-09, 
     & 2.51836D-03,3.05170D-05,1.00000D-50,1.00000D-50,7.44557D-04,
     & 1.83102D-06,1.00000D-50,6.16384D-03,2.48802D-06,1.39016D-05,
     & 1.00000D-50,5.43900D-07,1.00000D-50,1.24742D-03,3.13979D-06, 
     & 1.00997D-04,1.00000D-50,1.00000D-50,3.15133D-05,1.00000D-50,
     & 5.95044D-04,1.00000D-50,7.84705D-05,8.98518D-05,1.00000D-50,
     & 1.00000D-50,1.00000D-50,5.99970D-05,6.58862D-04,3.46503D-05, 
     & 2.36294D-05,6.28243D-06,3.15732D-04,2.60672D-06,1.51602D-05,
     & 1.00000D-50,7.48339D-08,6.61285D-06,1.00000D-50,2.23430D-06,
     & 1.65666D-03,9.36344D-05,6.95300D-05,3.37959D-06,1.80079D-05,
     & 1.17292D-05,1.39396D-03,7.72666D-05 /

!!!!!!!!!
!
! Young 2018      
!
! Non-Standard (Dmic, Rotation, KS66 atm)      
!!!!!!!!!        

! Solar abundances Asplund et al.2009 + Ne enhancement (Young+18) calibrated for solar model with atomic diffusion
! and KS66 atmosphere and rotation (prescription = Dh from Mathis et al 2018,
! Dv from Zahn 1992, magnetic braking from Matt 2015) ; alpha_MLT = 2.2236.
! Isotopic ratio Y = 0.271809, Z = 0.014227 ; calib : dr = 6.29e-4, dl = 2.99e-4, dZsx = 5.34e-5 
! from zini.out_1_3 in ~/STAREVOL/MODELS/Sun_rot_0.0133490.2663522.20572.2kms  / on port-dumont. 
! Thibaut (20/04/2020).
      data (xspref(14,i),i = 1,nis+6) / 1.D-50,0.713887000000000D+00,
     & 4.88804D-05,2.84079D-05,2.71809D-01,6.13729D-10,8.83310D-09,
     & 1.00000D-50,1.00000D-50,1.61848D-10,7.47192D-10,3.31101D-09, 
     & 2.39468D-03,2.88248D-05,1.00000D-50,1.00000D-50,7.06961D-04,
     & 2.79194D-06,1.00000D-50,5.86068D-03,2.37503D-06,1.32425D-05,
     & 1.00000D-50,5.17153D-07,1.00000D-50,1.67596D-03,4.27299D-06, 
     & 1.34809D-04,1.00000D-50,1.00000D-50,2.99635D-05,1.00000D-50,
     & 5.66073D-04,1.00000D-50,7.44028D-05,8.53338D-05,1.00000D-50,
     & 1.00000D-50,1.00000D-50,5.70465D-05,6.26447D-04,3.28632D-05, 
     & 2.25672D-05,5.97348D-06,3.00470D-04,2.44595D-06,1.41672D-05,
     & 1.00000D-50,7.11894D-08,6.28944D-06,1.00000D-50,2.12251D-06,
     & 1.57524D-03,8.90326D-05,6.61128D-05,3.21349D-06,1.71229D-05,
     & 1.11527D-05,1.32545D-03,7.34692D-05 /
      
c ##### Montréal / Montpellier #############      
! Abundances from Montreal code provide by O.Richard (Montpellier)
! calibrated for solar model with diffusion and without rotation and atmosphere ; 
! alpha_MLT = 1.6562.
! Isotopic ratio Y = 0.268169, Z = 0.01501039 ; calib : dr = ?, dl = ?, dZsx = ?
! from fort.21 in ~/Bureau/Thèse/Diffusion_atomique/Calibration_solaire/OR/Data/AGSS09_Diff_Ed
! on port-dumont. Thibaut (30/11/2017). (+ ~/ port-richard. Montpellier)
c$$$      data (xspref(7,i),i = 1,nis+6) / 1.D-50,0.71677700000000D+00,
c$$$     & 4.36207D-05,4.36207D-05,2.68169D-01,7.00000D-10,8.30000D-09,
c$$$     & 1.00000D-50,1.00000D-50,1.30000D-10,9.50000D-10,3.75000D-09, 
c$$$     & 2.65217D-03,2.95320D-05,1.00000D-50,1.00000D-50,7.75154D-04,
c$$$     & 1.00000D-50,1.00000D-50,6.41982D-03,1.00000D-50,1.00000D-50,
c$$$     & 1.00000D-50,1.00000D-50,1.00000D-50,1.41170D-03,1.00000D-50, 
c$$$     & 1.00000D-50,1.00000D-50,1.00000D-50,3.27783D-05,1.00000D-50,
c$$$     & 7.96248D-04,1.00000D-50,1.00000D-50,1.00000D-50,1.00000D-50,
c$$$     & 1.00000D-50,1.00000D-50,6.25099D-05,7.45911D-04,1.00000D-50, 
c$$$     & 1.00000D-50,6.52989D-06,3.45978D-04,1.00000D-50,1.00000D-50,
c$$$     & 1.00000D-50,1.00000D-50,1.00000D-50,1.00000D-50,9.17006D-06,
c$$$     & 1.62100D-03,1.00000D-50,1.00000D-50 / 

! Abundances from Montreal code provide by O.Richard (Montpellier)
! calibrated for pop 2 model without diffusion and without rotation and atmosphere ; 
! alpha_MLT = 1.5162.
!     Isotopic ratio Y = 0.248031, Z = 0.00039898 ; calib : dr = ?, dl = ?, dZsx = ?
! Fe/H = -1.75 and alpha/Fe = 0.3      
! from fort.21 in ~/Bureau/These/Data/Modeles_PopII/M0.8_FeH-1.75_Y_MoMo_NoDiff     
! on port-dumont. Thibaut (18/06/2019). (+ ~/ port-richard. Montpellier)
c$$$      data (xspref(3,i),i = 1,nis+6) / 1.D-50,0.75152660000000D+00,
c$$$     & 1.00000D-10,4.36207D-05,2.48031D-01,1.40000D-11,2.90000D-09,
c$$$     & 1.00000D-50,1.00000D-50,3.20000D-12,1.40000D-11,6.50000D-11, 
c$$$     & 4.28653D-05,4.77305D-07,1.00000D-50,1.00000D-50,1.25283D-05,
c$$$     & 1.00000D-50,1.00000D-50,2.07027D-04,1.00000D-50,1.00000D-50,
c$$$     & 1.00000D-50,1.00000D-50,1.00000D-50,4.55245D-05,1.00000D-50, 
c$$$     & 1.00000D-50,1.00000D-50,1.00000D-50,1.05704D-06,1.00000D-50,
c$$$     & 2.56775D-05,1.00000D-50,1.00000D-50,1.00000D-50,1.00000D-50,
c$$$     & 1.00000D-50,1.00000D-50,5.06352D-07,2.40542D-05,1.00000D-50, 
c$$$     & 1.00000D-50,2.10576D-07,1.11571D-05,1.00000D-50,1.00000D-50,
c$$$     & 1.00000D-50,1.00000D-50,1.00000D-50,1.00000D-50,2.95717D-07,
c$$$     & 3.03140D-05,1.00000D-50,1.00000D-50 / 
      
! Abundances from Montreal code provide by O.Richard (Montpellier)
! calibrated for pop 2 model with diffusion and without rotation and atmosphere ; 
! alpha_MLT = 1.6562.
!     Isotopic ratio Y = 0.248031, Z = 0.00039898 ; calib : dr = ?, dl = ?, dZsx = ?
! Fe/H = -1.75 and alpha/Fe = 0.3      
! from abondance.out in ~/Bureau/Thèse/Diffusion_atomique/Modeles_PopII/M0.8_FeH-1.75_Y_MoMo
! firts model at 30Myrs, first line (surface)      
! on port-dumont. Thibaut (03/05/2019). (+ ~/ port-richard. Montpellier)
c$$$      data (xspref(7,i),i = 1,nis+6) / 1.D-50,0.75152660000000D+00,
c$$$     & 1.00000D-10,4.36207D-05,2.48031D-01,1.40000D-11,2.90000D-09,
c$$$     & 1.00000D-50,1.00000D-50,3.20000D-12,1.40000D-11,6.50000D-11, 
c$$$     & 4.28653D-05,4.77305D-07,1.00000D-50,1.00000D-50,1.25283D-05,
c$$$     & 1.00000D-50,1.00000D-50,2.07027D-04,1.00000D-50,1.00000D-50,
c$$$     & 1.00000D-50,1.00000D-50,1.00000D-50,4.55245D-05,1.00000D-50, 
c$$$     & 1.00000D-50,1.00000D-50,1.00000D-50,1.05704D-06,1.00000D-50,
c$$$     & 2.56775D-05,1.00000D-50,1.00000D-50,1.00000D-50,1.00000D-50,
c$$$     & 1.00000D-50,1.00000D-50,5.06352D-07,2.40542D-05,1.00000D-50, 
c$$$     & 1.00000D-50,2.10576D-07,1.11571D-05,1.00000D-50,1.00000D-50,
c$$$     & 1.00000D-50,1.00000D-50,1.00000D-50,1.00000D-50,2.95717D-07,
c$$$     & 3.03140D-05,1.00000D-50,1.00000D-50 / 

c ##### Montréal / Montpellier ############# 
      
c..   atomic mass for heavy elements (> Cl) used for molecular opacity
c..   calculations (K, Ca, Ti, Cr, Mn, Fe, Ni)
      data (aheavy(i),i = 1,7) / 19.d0,20.d0,22.d0,24.d0,25.d0,26.d0,
     &     28.d0/

c...  coeff for K (dissociation constants, Table 2 Rossi & Maciel 83)
      data ((akap(i,j),j = 1,5),i = 1,11)/
     &     1.3590d1,-1.1523d1,6.6519d-2,-6.3343d-3,2.3781d-4, !CO
     &     1.2602d1,-4.9711d0,7.2280d-2,-6.1443d-3,2.0995d-4, !H2
     &     1.2198d1,-4.8641d0,6.6674d-2,-5.6697d-3,1.9546d-4, !OH
     &     2.5315d1,-1.0343d1,1.0807d-1,-8.7408d-3,2.9465d-4, !H2O
     &     1.2699d1,-8.0228d0,1.8920d-2,-6.1584d-4,1.1931d-6, !CN
     &     1.2654d1,-6.4287d0,3.4061d-2,-3.0069d-3,1.0475d-4, !C2
     &     1.3236d1,-1.0177d1,6.3664d-2,-5.9114d-3,2.1845d-4, !N2
     &     1.1793d1,-3.7508d0,2.1569d-2,-8.0102d-4,5.5938d-6, !CH
     &     1.2243d1,-6.8750d0,5.6616d-2,-5.2959d-3,1.9710d-4, !NO
     &     1.3114d1,-5.4035d0,2.9003d-1,-2.0968d-3,6.5775d-5, !O2
     &     1.1821d1,-3.8755d0,5.2603d-2,-4.1713d-3,1.3709d-4/ !NH
      data (molname(i),i=1,nmole+4)/'H','C','N','O','H2','H2O','OH',
     &     'CO','CN','C2','N2'/

c      data  iih,iic,iin,iio,iih2,iih2o,iioh,iico,iicn,iic2,iin2/1,2,3,4,
c     &     5,6,7,8,9,10,11/
      data iiH,iiC,iiN,iiO,iiCO,iiH2,iiOH,iiH2O,iiCN,iiC2,iiN2,iiCH,
     &     iiNO,iiO2,iiNH/001,002,003,004,005,006,007,008,009,010,011,
     &     012,013,014,015/
      data neqchim,nspecies,natoms/011,015,004/

      data ((ireac(i,j),j = 1,6),i = 1,11)/
     &     02,02,01,04,01,05,   !C + 0 = CO
     &     01,01,02,06,00,00,   !2*H = H2
     &     02,01,01,04,01,07,   !H + O = OH
     &     02,01,02,04,01,08,   !2*H + 0 = H2O
     &     02,02,01,03,01,09,   !C + N = CN
     &     01,02,02,10,00,00,   !2*C = C2
     &     01,03,02,11,00,00,   !2*N = N2
     &     02,01,01,02,01,12,   !C + H = CH
     &     02,02,01,04,01,13,   !N + O = NO
     &     01,04,02,14,00,00,   !2O = O2
     &     02,01,01,03,01,15/   !H + N = NH

c...  partition functions from Irwin (1981)
      data ((aIrwin(i,j),j = 1,12),i = 1,7)/
     &     -5.07315371d1,3.11367765d1,-7.50713553d0,9.14197587d-1, !F
     &     -5.59801362d-2,1.37482919d-3,
     &     4.70347846d1,-2.80260758d1,6.70365937d0,-7.72451873d-1, !F+
     &     4.28627164d-2,-9.09734601d-4,
     &     -5.21521085d2,4.10712446d2,-1.26375869d2,1.90976557d1,  !Ca
     &     -1.42269952d0,4.19085267d-2,
     &     1.65874025d3,-1.04392185d3,2.61702809d2,-3.26369424d1,  !Ca+
     &     2.02350114d0,-4.98594460d-2,
     &     -5.42849573d2,3.46635141d2,-8.90099156d1,1.15603978d1,  !Ti
     &     -7.60309452d-1,2.02863075d-2,
     &     7.41156498d0,-1.62989642d1,6.88299761d0,-1.14483565d0,  !Ti+
     &     8.56663842d-2,-2.39509300d-3,
     &     -8.80838409d2,5.72933536d2,-1.49289839d2,1.95455757d1,  !Cr
     &     -1.28737652d0,3.41709398d-2,
     &     3.32026608d3,-2.07680692d3,5.17436644d2,-6.4121353d1,   !Cr+
     &     3.94957079d0,-9.66571288d-2,
     &     -1.54913173d3,1.04762114d3,-2.82777415d2,3.81262914d1,  !Mn
     &     -2.56771825d0,6.91046934d-2,
     &     1.36807626d3,-8.39953845d2,2.05165234d2,-2.48590226d1,  !Mn+
     &     1.49207266d0,-3.54240524d-2,
     &     -1.15609527d3,7.46597652d2,-1.92865672d2,2.49658410d1,  !Fe
     &     -1.61934455d0,4.21182087d-2,
     &     2.71692895d2,-1.52697440d2,3.36119665d1,-3.56415427d0,  !Fe+
     &     1.80193259d-1,-3.38654879d-3,
     &     -1.02017263d3,6.68758438d2,-1.74841425d2,2.28263910d1,  !Ni
     &     -1.48693754d0,3.86538338d-2,
     &     -2.01279351d2,1.01910058d2,-1.92092108d1,1.62959379d0,  !Ni+
     &     -5.55991468d-2,3.32146945d-4/

c...  partition functions from Sauval & Tatum (1984)
      data ((aSauval(i,j),j = 1,10),i =1,nioniz)/
     &     0.30103,-0.00001,0.d0,0.d0,0.d0,                       !H
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !H+
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !He
     &     0.30103d0,-0.00001d0,0.d0,0.d0,0.d0,                   !He+
     &     0.31804d0,-0.20616d0,0.91456d0,-1.66121d0,1.04159d0,   !Li
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !Li+
     &     0.00801,-0.17135d0,0.62921d0,-0.58945d0,0.d0,          !Be
     &     0.30389d0,-0.00819d0,0.d0,0.d0,0.d0,                   !Be+
     &     0.78028d0,-0.01622d0,0.d0,0.d0,0.d0,                   !B
     &     0.00349d0,-0.01035d0,0.d0,0.d0,0.d0,                   !B+
     &     0.96752d0,-0.09452d0,0.08055d0,0.d0,0.d0,              !C
     &     0.77239d0,-0.03540d0,0.d0,0.d0,0.d0,                   !C+
     &     0.60683d0,-0.08674d0,0.30565d0,-0.28114d0,0.d0,        !N
     &     0.94968d0,-0.06463d0,-0.01291d0,0.d0,0.d0,             !N+
     &     0.95033d0,-0.05703d0,0.d0,0.d0,0.d0,                   !O
     &     0.60405d0,-0.03025d0,0.04525d0,0.d0,0.d0,              !O+
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !F
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !F+ 
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !Ne
     &     0.74847d0,-0.06562d0,-0.07088d0,0.d0,0.d0,             !Ne+
     &     0.30955d0,-0.17778d0,1.10594d0,-2.42847d0,1.70721d0,   !Na
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !Na+
     &     0.00556d0,-0.12840d0,0.81506d0,-1.79635d0,1.26292d0,   !Mg
     &     0.30257d0,-0.00451d0,0.d0,0.d0,0.d0,                   !Mg+
     &     0.76786d0,-0.05207d0,0.14713d0,-0.21376d0,0.d0,        !Al
     &     0.00334d0,-0.00995d0,0.d0,0.d0,0.d0,                   !Al+
     &     0.97896d0,-0.19208d0,0.04753d0,0.d0,0.d0,              !Si
     &     0.75647d0,-0.05490d0,-0.10126d0,0.d0,0.d0,             !Si+
     &     0.64618d0,-0.31132d0,0.68633d0,-0.47505d0,0.d0,        !P
     &     0.93588d0,-0.18848d0,0.08921d0,-0.22447d0,0.d0,        !P+
     &     0.95254d0,-0.15166d0,0.02340d0,0.d0,0.d0,              !S
     &     0.61971d0,-0.17465d0,0.48283d0,-0.39157d0,0.d0,        !S+
     &     0.74465d0,-0.07389d0,-0.06965d0,0.d0,0.d0,             !Cl
     &     0.92728d0,-0.15913d0,-0.01983d0,0.d0,0.d0,             !Cl+
     &     0.34419d0,-0.48157d0,1.92563d0,-3.17826d0,1.83211d0,   !K
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !K+
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !Ca
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !Ca+ 
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !Ti
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !Ti+ 
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !Cr
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !Cr+ 
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !Mn
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !Mn+ 
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !Fe
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !Fe+ 
     &     0.d0,0.d0,0.d0,0.d0,0.d0,                              !Ni
     &     0.d0,0.d0,0.d0,0.d0,0.d0 /                             !Ni+ 

c..   first ionization energy
      data (Eion(i),i = 1,nioniz)/13.5984d0,24.5874d0,5.3917d0,9.3227d0,
     &     8.2980d0,11.2603d0,14.5341d0,13.6181d0,17.4228d0,21.5645d0,
     &     5.1391d0,7.6462d0,5.9858d0,8.1517d0,10.4867d0,10.3600d0,
     &     12.9676d0,4.3407d0,6.1132d0,6.8281d0,6.7665d0,7.4340d0,
     &     7.9024d0,7.6398d0 /



      end

************************************************************************

      SUBROUTINE fdiff (eqd,cd,x,f1,f2,f3,id,imin,imax)

************************************************************************
*     compute the diffusion matrix elements and derivatives
*     eq[i,neqd1] = dYi = dt*Fi(Yj)  (eq of nucleosynthesis at shell id)
*                 = dt * d (D dYi/dm)/dm
*     eq[i,j=1,neqd]        = [dFi/dYj]_{id-1} (jacobian)
*     eq[i,j=neqd+1,2neqd]  = [dFi/dYj]_{id}   (jacobian)
*     eq[i,j=2neqd+1,3neqd] = [dFi/dYj]_{id+1} (jacobian)
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.diff'
      include 'evolcom.conv'
      include 'evolcom.teq'

      integer j,id,neqd,neqd1,neqd2,imin,imax

      double precision x
      double precision cd,f1,f2,f3,eqd
      double precision xbcz,x2,x3

      dimension cd(nsh),f1(nsh),f2(nsh),f3(nsh),eqd(nsp,3*nsp+1),
     &     x(nsh,nsp)

      neqd = nsp
      neqd2 = 2*neqd
      neqd1 = 3*neqd+1
      kcentd = 1
      ksurfd = 1

c..   start filling matrix elements at j=2 because neutrons don't diffuse
***   interior
      if (id.gt.imin.and.id.lt.imax) then
         x3 = f3(id)*cd(id)
         x2 = f2(id)*cd(id+1)
         do j = 2,neqd
            eqd(j,neqd1) = eqd(j,neqd1)+x2*(x(id+1,j)-x(id,j))-
     &           x3*(x(id,j)-x(id-1,j))
            eqd(j,j) = eqd(j,j)+x3
            eqd(j,neqd+j) = -x2-x3+eqd(j,neqd+j)
            eqd(j,neqd2+j) = eqd(j,neqd2+j)+x2
         enddo
         return
      endif


***   inner boundary condition
      if (id.eq.imin) then
c..   center or continuity equation
         if (kcentd.eq.1.or.id.eq.1) then
            x2 = f2(id)*cd(id+1)
            do j = 2,neqd
               eqd(j,neqd1) = eqd(j,neqd1)+x2*(x(id+1,j)-x(id,j))
               eqd(j,j) = eqd(j,j)-x2
               eqd(j,neqd+j) = eqd(j,neqd+j)+x2
            enddo
            return
         endif
c..   reservoir
         if (kcentd.eq.2) then
            do j = 1,neqd
               xbcz = cd(id+1)*yconvlc(j)*f1(id)/delmc
               eqd(j,neqd1) = eqd(j,neqd1)-xbcz*(x(id+1,j)/x(id,j)-1.d0)
               eqd(j,j) = eqd(j,j)+xbcz*x(id+1,j)/(x(id,j)*x(id,j))
               eqd(j,neqd+j) = eqd(j,neqd+j)-xbcz/x(id,j)
            enddo
         endif
c..   Y = cste
c         if (kcentd.eq.3) then
c            do j = 2,neqd
c               eqd(j,neqd1) = -x(id,j)+x(id-1,j)
c               eqd(j,j) = -1.d0
c            enddo
c         endif
      endif

***   outer boundary condition
      if (id.eq.imax) then
c..   surface or continuity equation
         if (ksurfd.eq.1.or.id.eq.nmod1) then
            x3 = f3(id)*cd(id)
            do j = 2,neqd
               eqd(j,neqd1) = eqd(j,neqd1)-x3*(x(id,j)-x(id-1,j))
               eqd(j,j) = eqd(j,j)+x3
               eqd(j,neqd+j) = eqd(j,neqd+j)-x3
            enddo
            return
         endif
c..   reservoir
         if (ksurfd.eq.2) then
            do j = 2,neqd
               xbcz = cd(id)*yconvls(j)*f1(id)/delms
               eqd(j,neqd1) = eqd(j,neqd1)-xbcz*(1.d0-x(id-1,j)/x(id,j))
               eqd(j,j) = eqd(j,j)+xbcz/x(id,j)
               eqd(j,neqd+j) = eqd(j,neqd+j)-xbcz*x(id-1,j)/
     &              (x(id,j)*x(id,j))
            enddo
         endif
c..   Y = cste
c          if (ksurfd.eq.3) then
c             do j = 2,neqd
c                eqd(j,neqd1) = -x(id,j)+x(id+1,j)
c                eqd(j,neqd+j) = -1.d0
c             enddo
c          endif
      endif


      return
      end


************************************************************************

      SUBROUTINE freac (vv,eqd,y,dt,lflag)

************************************************************************
*  compute the nuclear matrix elements and derivatives
*  eq[i,neqd1] = dt*dYi/dt = dt*Fi(Yj)     (nucleosynthesis equation)
*  eq[i,j=neqd+1,2neqd] = dt*d(dYi/dt)/dYj      (jacobian)
*  lflag = .true.  : nucleosynthesis alone
*  lflag = .false. : coupling nucleosynthesis-mixing
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.nuc'

      integer neqd,neqd1
      integer kna,kia,knb,kib,knd,kid,knc,kic
      integer i,nn

      double precision dt,z,za,zb
      double precision vv(nreac),eqd(nsp,3*nsp+1),y(nsp)

      logical lflag

      neqd = nsp
      neqd1 = 3*neqd+1


c..  reaction : kna X(kia) + knb X(kib) --> knd X(kid) + knc X(kic)
c..  reaction : k2  X(k1)  + k4  X(k3)  --> k8  X(k7)  + k6  X(k5)
c          ex :  1    C12  +  1    He   -->  0  gamma  +  1   O16

      if (lflag) then
         nn = 0
      else
         nn = neqd
      endif
      do i = 1,nreac
         kna = k2(i)
         kia = k1(i)
         knb = k4(i)
         kib = k3(i)
         knd = k6(i)
         kid = k5(i)
         knc = k8(i)
         kic = k7(i)
         if (kna.eq.1.and.knb.eq.1) then
            z = vv(i)*y(kia)*y(kib)
            za = vv(i)*y(kib)
            zb = vv(i)*y(kia)
         else
            if (knb.eq.0) then
               if (kna.eq.1) then
                  z = vv(i)*y(kia)
                  za = vv(i)
               else
                  z = vv(i)*y(kia)**kna/fact(kna)
                  za = vv(i)*y(kia)**(kna-1)/fact(kna-1)
               endif
               zb = 0.d0
            else
               z = vv(i)*y(kia)**kna*y(kib)**knb/(fact(kna)*fact(knb))
               za = vv(i)*y(kia)**(kna-1)*y(kib)**knb/
     &              (fact(kna-1)*fact(knb))
               zb = vv(i)*y(kia)**kna*y(kib)**(knb-1)/
     &              (fact(kna)*fact(knb-1))
            endif
         endif
         z = z*dt
         za = za*dt
         zb = zb*dt

         eqd(kia,neqd1) = eqd(kia,neqd1)-kna*z
         eqd(kib,neqd1) = eqd(kib,neqd1)-knb*z
         eqd(kic,neqd1) = eqd(kic,neqd1)+knc*z
         eqd(kid,neqd1) = eqd(kid,neqd1)+knd*z

c..  jacobian
         eqd(kia,nn+kia) = eqd(kia,nn+kia)-kna*za
         eqd(kia,nn+kib) = eqd(kia,nn+kib)-kna*zb
         eqd(kib,nn+kia) = eqd(kib,nn+kia)-knb*za
         eqd(kib,nn+kib) = eqd(kib,nn+kib)-knb*zb
         eqd(kic,nn+kia) = eqd(kic,nn+kia)+knc*za
         eqd(kic,nn+kib) = eqd(kic,nn+kib)+knc*zb
         eqd(kid,nn+kia) = eqd(kid,nn+kia)+knd*za
         eqd(kid,nn+kib) = eqd(kid,nn+kib)+knd*zb

      enddo

      return
      end


************************************************************************

      SUBROUTINE fscreen (t,ro,xsp,rhpsi,anuc,znuc)

************************************************************************
* Calculation of the screening factors                                 *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.ion'
      include 'evolcom.teq'

      integer i,ii,k,l

      double precision t(nsh),ro(nsh),xsp(nsh,nsp),rhpsi(nsh),anuc(nsp),
     &     znuc(nsp)
      double precision fscr
      double precision screen
      double precision dzi
      double precision y(nsp),zx,zy,zbar,ztilde,l0
      double precision tburn,tnucmin,fnetdt,fnetdt0,tolnuc,tolnuc0
      double precision nelectronk,zeta,h0,kfact,hfact,zstar,kscreen

      common /funcscr/ fscr(nsh,nscr)
      common /nuc_param/ tburn,tnucmin,fnetdt,fnetdt0,tolnuc,tolnuc0

      external screen

      kfact = pim4*ech**2/boltz
      hfact = -0.4d0*pi*kfact

      fscr = 1.d0


***   Charged reaction factors
      do k = 1,nmod
         nelectronk = ro(k)*avn*mueinv(k)
         if (t(k).gt.tnucmin) then
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
            zstar = ztilde/zbar
            kscreen = dsqrt(kfact*nelectronk*(zstar+rhpsi(k))/t(k))
            zeta = 3.d0*kscreen**3/(pim4*nelectronk)
            h0 = hfact*nelectronk**2/(t(k)*kscreen**5)
            zbar = zbar/zx
            ztilde = dsqrt(ztilde/zx+rhpsi(k)*zbar)
            l0 = 1.88d8*dsqrt(ro(k)*zx/t(k)**3)
c..   nbz-1 because nbz = heavy
            do i = 1,nbz-1
               ii = nbz-1+i
               dzi = dble(i)
               fscr(k,i) = screen (y,zbar,ztilde,zeta,l0,h0,dzi,1.d0)
               fscr(k,ii) = screen (y,zbar,ztilde,zeta,l0,h0,dzi,2.d0)
            enddo
            fscr(k,2*nbz-1) = screen(y,zbar,ztilde,zeta,l0,h0,6.d0,6.d0)
            fscr(k,2*nbz) = screen (y,zbar,ztilde,zeta,l0,h0,8.d0,8.d0)
            fscr(k,2*nbz+1) = screen(y,zbar,ztilde,zeta,l0,h0,8.d0,6.d0)
         endif
      enddo

      return
      end


************************************************************************

      SUBROUTINE guess

************************************************************************
*  estimate the new stellar model at the first iteration               *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.mod'
      include 'evolcom.rot'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'

      double precision xmod,xa,eps
      double precision t,vt,r,vr

      common /varia/ xmod(nsh,neq),xa(nsh,neq),eps(neq)
      common /vartvt/ t(nsh),vt(nsh)
      common /varrvr/ r(nsh),vr(nsh)

      xa(1:nmod,1:neq) = xmod(1:nmod,1:neq)

c..   guess new structure
      xmod(1:nmod,1:neq) = xmod(1:nmod,1:neq)+(xmod(1:nmod,1:neq)-
     &     xa(1:nmod,1:neq))*phi

      t(1:nmod) = exp(xmod(1:nmod,4))
      vt(1:nmod) = exp(xa(1:nmod,4))
      r(1:nmod) = exp(xmod(1:nmod,2))
      vr(1:nmod) = exp(xa(1:nmod,2))
      vro(1:nmod) = ro(1:nmod)
      vp(1:nmod) = p(1:nmod)
      ve(1:nmod) = e(1:nmod)
      vs(1:nmod) = s(1:nmod)
      vxsp(1:nmod,1:nsp) = xsp(1:nmod,1:nsp)
      omega(1:nmod) = vomega(1:nmod)
      if (dmaccr.gt.0.d0) then
         eacc(1:nmod) = lacc*facc(1:nmod)/(1.d0+facc(1:nmod))
         eaccdr(1:nmod) = lacc*dfaccdr(1:nmod)/(1.d0+facc(1:nmod))**2
      else
         eacc(1:nmod) = 0.d0
         eaccdr(1:nmod) = 0.d0
      endif
      nnucacc(1:nmaxconv) = 0
      iprint(1:nmaxconv) = 0
      rnucacc = 0.d0

      return
      end
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
      print *,'dtin,dtn',dtin,dtn
c..   luminosity variation too large. Do not increase dtn
      if (ifail.lt.0.and.mtini.lt.300.d0) then
         dtn = min(dtn,dtalt)
         write (nout,*) 'luminosity variation too large, dtn not ',
     &        'increased',dtn,dtalt
      endif

*____________________________________________________________
***   Treatment of core hydrogen burning exhaustion (modif Louis 09/2018)
*------------------------------------------------------------

      xmaxH = 1.d-1

      if (nphase.eq.2.and.xsp(1,ih1).le.xmaxH) dtn = min(dtn,dtmax0*
     &     (0.05d0+0.95d0*xsp(1,ih1)/xmaxH))
      


*____________________________________________________________
***   Treatment of turn-off
*------------------------------------------------------------
c      if (ntprof.eq.5.and.nphase.eq.3.and..not.rgbphase) then !Test TD 03/2020
c         dtn = 1d4/seci
c      endif


      
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
      if (tau(i).lt.taulim.and.(ntprof.eq.1.or.ntprof.ge.3.and.
     &     .not.ntp2)) then
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
         eq(1,16) = eq(1,16)+dptau(i)
         eq(1,3) = eq(1,3)+dptau(i)*dkapdf(im)*kiminv
         eq(1,4) = eq(1,4)+dptau(i)*dkapdt(im)*kiminv
         eq(1,7) = eq(1,7)-2.d0*dptau(i)
         eq(1,8) = eq(1,8)+dptau(i)*dkapdf(i)*kiinv
         eq(1,9) = eq(1,9)+dptau(i)*dkapdt(i)*kiinv
         eq(1,10) = eq(1,10)+dptau(i)/lum(i)
      endif

***   definition of lagrangian velocity

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

      eq(3,16) = (r(ip)**3-r(i)**3)/dm(i)-sph
      eq(3,7) = -3.d0*r(i)*r(i)*r(i)/dm(i)
      eq(3,8) = sph*drodf(i)/ro(i)
      eq(3,9) = sph*drodt(i)/ro(i)
      eq(3,12) = 3.d0*r(ip)*r(ip)*r(ip)/dm(i)

***   energy conservation

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
         if (tau(ip).lt.taulim.and.(ntprof.eq.1.or.ntprof.ge.3.and.
     &        .not.ntp2)) then
            dxsi = -xsitau(i)/(g*m(ip)+r(ip)*r(ip)*accel(ip))
            ddtau = 0.d0
            if (abs(dqdtau(i)).gt.1.d-6) ddtau = ddqtau(i)/dqdtau(i)
c            plo = 3.d0*lum(ip)*kapm(ip)*dmkb/(64.d0*sig*pi*pim4*
c     &           (r(ip)*tm(ip))**4)*ft(ip)!/fp(ip)
            gmrua1 = (gmr(ip)+accel(ip))*fp(ip)
            
            ablainv = 1.d0/abla(ip)
            sqrp = sqrt(p(i)*p(ip))

            plo = gmrua1*abla(ip)*dmkb/(pim4*r(ip)*r(ip)*sqrp)

            eq(5,16) = dlnt1+plo*(1.d0+xsitau(ip))
            eq(5,6) = plo*(abdu1(ip)*ablainv+dynfac*dtninv*psi(ip)/
     &           gmrua1)*(1.d0+xsitau(ip))+plo*dxsi*r(ip)*r(ip)*dtninv*
     &           psi(ip)
            eq(5,7) = plo*xsitau(ip)*ddtau*kap(i)*ro(i)
            eq(5,8) = plo*(abdf1(ip)*ablainv-wip*dpdf(i)/p(i))*(1.d0+
     &           xsitau(ip))+plo*xsitau(ip)*(wip*dkapdf(i)*kiinv)
c     &           xsitau(ip))+plo*xsitau(ip)*(wip*dkapdf(i)*kiinv+ddtau*
c     &           dtaudf(i))
            eq(5,9) = plo*(abdt1(ip)/abla(ip)-wip*dpdt(i)/p(i))*(1.d0+
     &           xsitau(ip))+plo*xsitau(ip)*(wip*dkapdt(i)*kiinv)!+ddtau*
c     &           xsitau(ip))+plo*xsitau(ip)*(wip*dkapdt(i)*kiinv+ddtau*
c     &           dtaudt(i))
     &           -1.d0
            eq(5,11) = plo*(abdu2(ip)*ablainv+dynfac*(dtninv-dtipsi1)/
     &           gmrua1)*(1.d0+xsitau(ip))+plo*dxsi*r(ip)*r(ip)*dtninv*
     &           (1.d0-psi(ip))
            eq(5,12) = plo*(abdr(ip)*ablainv-2.d0*(1.d0+gmr(ip)/gmrua1)/
     &           r(ip))*(1.d0+xsitau(i))+plo*dxsi*2.d0*r(ip)*accel(ip)
            eq(5,13) = plo*(abdf2(ip)*ablainv-wjp*dpdf(ip)/p(ip))*
     &           (1.d0+xsitau(ip))+plo*xsitau(ip)*wjp*dkapdf(ip)*kipinv
            eq(5,14) = plo*(abdt2(ip)/abla(ip)-wjp*dpdt(ip)/p(ip))*
     &           (1.d0+xsitau(ip))+plo*xsitau(ip)*wjp*dkapdt(ip)*kipinv+
     &           1.d0
            eq(5,15) = plo*(abdl(ip)*ablainv*(1.d0+xsitau(ip))+
     &           xsitau(ip)/lum(ip))

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


************************************************************************

      SUBROUTINE interpmesh (il,ir,ls)

************************************************************************
* Interpolate variables in case of shell addition/removal              *
*                                                                      *
* delete shell                                                         *
* ------------                                                         *
* before    |      |     |    |        |                               *
*          il-1    il  il+1   ir     ir+1                              *
* after     |      |          |        |                               *
*          il-1    il         ir     ir+1                              *
*                                                                      *
* add shell                                                            *
* ---------                                                            *
* before    |      |          |        |                               *
*          il-1    il         ir     ir+1                              *
* after     |      |     |    |        |                               *
*          il-1    il  il+1   ir     ir+1                              *
*                                                                      *
* ls = 0 : shell shifting                                              *
* ls = 1 : shell suppression                                           *
* ls = 2 : shell addition, interpolate il+1 from il and ir (default)   *
* ls = 3 : shell addition, interpolate il+1 from il-1 and ir           *
* ls = 4 : shell addition, interpolate il+1 from il-1 and il           *
*   treatment of mass loss/accretion                                   *
* ls = 5 : rescaling the pseudo-lagrangian coordinate                  *
*   treatment of shock fronts                                          *
* ls = 6 : shell addition, interpolate il from ir and ir+1             *
* ls = 7 : shell addition, special treatment of central shell          *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.eng'
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

      integer ls,il,il0,il1,ir,ir1,ir2
      integer ia,ib,ic,id,iab,icd,ka,kb
      integer iqmrbot,iqmrtop
      integer k,l
      integer ja,jb,ixa,ixb,ixc,ixd

      double precision ai,bi,ci,di,del,ak,bk
      double precision midmr,vmidmr,qmr,vqmr
      double precision r2dm0,r2dm1,r2dm2,r2dm3,r2dm4
      double precision zj1,zj2,zj3,rm1,rm2,rm3
      double precision xai,xbi,xci,xdi
      double precision domj,vconv

      common /newmesh/ qmr(nsh),vqmr(nsh),iqmrbot,iqmrtop


      il0 = il-1
      il1 = il+1
      ir1 = ir+1
      ir2 = ir+2


*__________________________
*
***      SHELL SHIFTING
*
*--------------------------


      if (ls.eq.0) then
***   interface variables
         vu(il) = vu(ir)
         vlum(il) = vlum(ir)
         vr(il) = vr(ir)

         u(il) = u(ir)
         r(il) = r(ir)
         lum(il) = lum(ir)

         vsconv(il) = vsconv(ir)
         vhconv(il) = vhconv(ir)
         mr(il) = mr(ir)
         m(il) = m(ir)

***   centered variables
         lnf(il) = lnf(ir)
         lnt(il) = lnt(ir)
         ro(il) = ro(ir)
         p(il) = p(ir)
         e(il) = e(ir)
         s(il) = s(ir)

         vlnf(il) = vlnf(ir)
         vlnt(il) = vlnt(ir)
         vro(il) = vro(ir)
         vp(il) = vp(ir)
         ve(il) = ve(ir)
         vs(il) = vs(ir)

         crz(il) = crz(ir)
         dm(il) = dm(ir)
         tau(il) = tau(ir)
         vmueinv(il) = vmueinv(ir)
         venucl(il) = venucl(ir)
         vpvisc(il) = vpvisc(ir)
         exposure(il) = exposure(ir)

         do l = 1,nsp
            xsp(il,l) = xsp(ir,l)
         enddo
         vomega(il) = vomega(ir)
         if (irotbin.eq.1) then
            auxs(il) = auxs(ir)
            urs(il) = urs(ir)
            xpsis(il) = xpsis(ir)
            xlambdas(il) = xlambdas(ir)
         endif
         if (iaccr.gt.0) then
            facc(il) = facc(ir)
            macc(il) = macc(ir)
            vfacc(il) = vfacc(ir)
            dfaccdr(il) = dfaccdr(ir)
            vacc(il) = vacc(ir)
            tacc(il) = tacc(ir)
         endif
      endif



*___________________________________________________________
***
***                       SUPPRESSION
***
***   adjustment of shells adjacent to a suppressed one
***   only applies to centered variables (.i.e. not lum,r,u)
***   new version
*-----------------------------------------------------------


      if (ls.eq.1) then
         ai = dm(il)/(dm(il)+dm(ir))
         bi = 1.d0-ai

         lnf(il) = ai*lnf(il)+bi*lnf(ir)
         lnt(il) = ai*lnt(il)+bi*lnt(ir)
         p(il) = exp(ai*log(p(il))+bi*log(p(ir)))
         e(il) = exp(ai*log(e(il))+bi*log(e(ir)))
         if (s(il).gt.0.d0.and.s(ir).gt.0.d0) then
            s(il) = exp(ai*log(s(il))+bi*log(s(ir)))
         else
            s(il) = ai*s(il)+bi*s(ir)
         endif
         ro(il) = exp(ai*log(ro(il))+bi*log(ro(ir)))

         vlnf(il) = ai*vlnf(il)+bi*vlnf(ir)
         vlnt(il) = ai*vlnt(il)+bi*vlnt(ir)
         vp(il) = exp(ai*log(vp(il))+bi*log(vp(ir)))
         ve(il) = exp(ai*log(ve(il))+bi*log(ve(ir)))
         if (vs(il).gt.0.d0.and.vs(ir).gt.0.d0) then
            vs(il) = exp(ai*log(vs(il))+bi*log(vs(ir)))
         else
            vs(il) = ai*vs(il)+bi*vs(ir)
         endif
         vro(il) = exp(ai*log(vro(il))+bi*log(vro(ir)))

         tau(il) = ai*tau(il)+bi*tau(ir)
         vmueinv(il) = ai*vmueinv(il)+bi*vmueinv(ir)
         venucl(il) = ai*venucl(il)+bi*venucl(ir)
         vpvisc(il) = ai*vpvisc(il)+bi*vpvisc(ir)
c         exposure(il) = ai*exposure(il)+bi*exposure(ir)
         exposure(il) = max(exposure(il),exposure(ir))

         if (xsp(il,2).ne.xsp(ir,2)) then
            do l = 1,nsp
               xsp(il,l) = ai*xsp(il,l)+bi*xsp(ir,l)
            enddo
         endif

***   Interpolation ensures angular momentum conservation

         if (irotbin.eq.1) then
c           vomega(il) = ai*vomega(il)+bi*vomega(ir)
            if (il.gt.1) then 
               r2dm0 = dm(il0)*(r(il0)**2+r(il)*r(il0)+r(il)**2)/3.d0
            else
               r2dm0 = 0.d0
            endif
            r2dm3 = (dm(il)+dm(ir))*(r(il)**2+r(il)*r(ir1)+r(ir1)**2)
     &           /3.d0
            r2dm1 = dm(il)*(r(il)**2+r(il)*r(ir)+r(ir)**2)/3.d0
            r2dm2 = dm(ir)*(r(ir)**2+r(ir)*r(ir1)+r(ir1)**2)/3.d0
            r2dm4 = dm(ir1)*(r(ir2)**2+r(ir2)*r(ir1)+r(ir1)**2)/3.d0

c..   Increment domj is applied to omega at shells il and ir1 in order
c..   to allow angular momentum conservation when shell ir is suppressed
c            domj = (r2dm1*(vomega(il)+vomega(ir))+r2dm2*(vomega(ir)+
c     &           vomega(ir1))-r2dm3*(vomega(il)+vomega(ir1)))/
c     &           (r2dm0+r2dm4+2.d0*r2dm3)
c            vomega(il) = vomega(il)+domj
c            vomega(ir1) = vomega(ir1)+domj

c..   New formulation for the redistribution of omega
            domj = (vomega(ir)*(r2dm2+r2dm1)+vomega(il)*(r2dm1-r2dm3)+
     &           vomega(ir1)*(r2dm2-r2dm3))/(vomega(ir1)*(r2dm3+r2dm4)+
     &           vomega(il)*(r2dm3+r2dm0))

            vomega(il) = vomega(il)*(1.d0+domj)
            vomega(ir1) = vomega(ir1)*(1.d0+domj)

            auxs(il) = ai*auxs(il)+bi*auxs(ir)
            urs(il) = ai*urs(il)+bi*urs(ir)
            xpsis(il) = ai*xpsis(il)+bi*xpsis(ir)
            xlambdas(il) = ai*xlambdas(il)+bi*xlambdas(ir)
         endif
      endif



*______________________________________________________
***
***                       ADDITION
***
***   interpolation of the variables in the added shell
***   define interpolation coefficients
*------------------------------------------------------



      if (ls.eq.2.or.ls.eq.3.or.ls.eq.5.or.ls.eq.6.or.ls.eq.7) then

***   interpolation at il+1 with variables defined at il and ir
*     default
         if (ls.eq.2) then
            ka = il
            kb = ir
            ak = 0.5d0
            bk = 0.5d0
c..   centered variables
            iab = il1
            bi = dm(il1)/(dm(ir)+2.d0*dm(il1))
            ai = 1.d0-bi
            ia = il
            ib = ir
            icd = il
            di = dm(il)/(dm(il0)+2.d0*dm(il))
            ci = 1.d0-di
            ic = il
            id = il0
c..   chemical composition
            xai = ai
            xbi = bi
            xci = ci
            xdi = di
            ixa = ia
            ixb = ib
            ixc = ic
            ixd = id
         endif

***   interpolation at il+1 with variables defined at il-1 and ir
         if (ls.eq.3) then
            ka = il
            kb = ir
            ak = 0.5d0
            bk = 0.5d0
c..   centered variables
            del = dm(il0)+4.d0*dm(il)+dm(ir)
            del = 1.d0/del
            iab = il1
            ai = (3.d0*dm(il)+dm(ir))*del
            bi = 1.d0-ai
            ia = il0
            ib = ir
            icd = il
            ci = (dm(il)+dm(ir))*del
            di = 1.d0-ci
            ic = il0
            id = ir
c..   chemical composition
            xai = ai
            xbi = bi
            xci = ci
            xdi = di
            ixa = ia
            ixb = ib
            ixc = ic
            ixd = id
         endif

***   special treatment of inner shock-front
***   extrapolation at ir from variables defined at il-1 and il
         if (ls.eq.5) then
            ka = il0
            kb = il
            ak = -dm(il)/dm(il0)
            bk = 1.d0-ak
c..   centered variables
            iab = il
            bi = dm(il)/(dm(il0)+2.d0*dm(il))
            ai = 1.d0-bi
            ia = il
            ib = il0
            icd = il1
            ci = -2.d0*dm(il)/(dm(il)+dm(il0))
            di = 1.d0-ci
            ic = il0
            id = il
c..   chemical composition
            xai = ai
            xbi = bi
            xci = ci
            xdi = di
            ixa = ia
            ixb = ib
            ixc = ic
            ixd = id
         endif

***   special treatment of upper shock-front
***   extrapolation at il from variables defined at ir and ir+1
         if (ls.eq.6) then
            ka = ir
            kb = ir+1
            bk = -dm(il)/dm(ir)
            ak = 1.d0-bk
c..   centered variables
            iab = il1
            ai = dm(il1)/(dm(ir)+2.d0*dm(il1))
            bi = 1.d0-ai
            ia = ir
            ib = il
            icd = il
            ci = -dm(il1)/(2.d0*dm(il1)+dm(ir))
            di = 1.d0-ci
            ic = ir
            id = il
c..   chemical composition
            xai = ai
            xbi = bi
            xci = ci
            xdi = di
            ixa = ia
            ixb = ib
            ixc = ic
            ixd = id
         endif

***   special treatment of central quantities
         if (ls.eq.7) then
            ka = il
            kb = ir
            ak = 0.5d0
            bk = 0.5d0
c..   centered variables
            iab = il1
            ia = ir
            ib = il
            ai = dm(iab)/(dm(ia)+2.d0*dm(iab))
            bi = 1.d0-ai
            icd = il
            ic = il
            id = il1
            ci = 2.d0
            di = -1.d0
c..   chemical composition
            xai = ai
            xbi = bi
            xci = ci
            xdi = di
            ixa = ia
            ixb = ib
            ixc = ic
            ixd = id
         endif

*------------------------------------
*
***    INTERPOLATION (shell addition)
*
*------------------------------------

c..   interface variables
         vr(il1) = (0.5d0*(vr(il)**3+vr(ir)**3))**pw13
         vu(il1) = ak*vu(ka)+bk*vu(kb)
         vlum(il1) = ak*vlum(ka)+bk*vlum(kb)

         r(il1) = (0.5d0*(r(il)**3+r(ir)**3))**pw13
         u(il1) = ak*u(ka)+bk*u(kb)
         lum(il1) = ak*lum(ka)+bk*lum(kb)

         vconv = ak*vsconv(ka)+bk*vsconv(kb)
         if (vconv.gt.0.d0.or.ls.eq.2) then
            vsconv(il1) = vconv
            vhconv(il1) = ak*vhconv(ka)+bk*vhconv(kb)
         else
            vsconv(il1) = 0.5d0*(vsconv(il)+vsconv(ir))
            vhconv(il1) = 0.5d0*(vhconv(il)+vhconv(ir))
         endif
         crz(il1) = crz(il)

***   coefficients for omega interpolation with
***   angular momentum conservation
         if (irotbin.eq.1) then
            rm1 = (r(il)**2+r(il)*r(ir)+r(ir)**2)/3.d0
            rm2 = (r(il)**2+r(il)*r(il1)+r(il1)**2)/3.d0
            rm3 = (r(ir)**2+r(ir)*r(il1)+r(il1)**2)/3.d0
            zj1 = rm1*(dm(il)+dm(il1))
            zj2 = rm2*dm(il)
            zj3 = rm3*dm(il1)
            ja = il
            jb = ir
            vomega(il1) = (vomega(ja)*(zj1-zj2)+vomega(jb)*(zj1-zj3))/
     &           (zj2+zj3)
            auxs(il1) = ak*auxs(ka)+bk*auxs(kb)
            urs(il1) = ak*urs(ka)+bk*urs(kb)
            xpsis(il1) = ak*xpsis(ka)+bk*xpsis(kb)
            xlambdas(il1) = ak*xlambdas(ka)+bk*xlambdas(kb)
         else
            vomega(il1) = ak*vomega(ka)+bk*vomega(kb)
         endif

***   centered variables
***   WARNING : order matters, first interpolate at iab and then at icd
         lnf(iab) = ai*lnf(ia)+bi*lnf(ib)
         lnt(iab) = ai*lnt(ia)+bi*lnt(ib)
         p(iab) = exp(ai*log(p(ia))+bi*log(p(ib)))
         e(iab) = exp(ai*log(e(ia))+bi*log(e(ib)))
         s(iab) = exp(ai*log(s(ia))+bi*log(s(ib)))
         ro(iab) = exp(ai*log(ro(ia))+bi*log(ro(ib)))

         vlnf(iab) = ai*vlnf(ia)+bi*vlnf(ib)
         vlnt(iab) = ai*vlnt(ia)+bi*vlnt(ib)
         vp(iab) = exp(ai*log(vp(ia))+bi*log(vp(ib)))
         ve(iab) = exp(ai*log(ve(ia))+bi*log(ve(ib)))
         vs(iab) = exp(ai*log(vs(ia))+bi*log(vs(ib)))
         vro(iab) = exp(ai*log(vro(ia))+bi*log(vro(ib)))

         tau(iab) = ai*tau(ia)+bi*tau(ib)
         vmueinv(iab) = ai*vmueinv(ia)+bi*vmueinv(ib)
         venucl(iab) = ai*venucl(ia)+bi*venucl(ib)
         if (vpvisc(ia).gt.0.d0.and.vpvisc(ib).gt.0.d0) then
            vpvisc(iab) = exp(ai*log(vpvisc(ia))+bi*log(vpvisc(ib)))
         else
            vpvisc(iab) = ai*vpvisc(ia)+bi*vpvisc(ib)
         endif
         exposure(iab) = exposure(ia)

         lnf(icd) = ci*lnf(ic)+di*lnf(id)
         lnt(icd) = ci*lnt(ic)+di*lnt(id)
         p(icd) = exp(ci*log(p(ic))+di*log(p(id)))
         e(icd) = exp(ci*log(e(ic))+di*log(e(id)))
         if (s(ic).gt.0.d0.and.s(id).gt.0.d0) then
            s(icd) = exp(ci*log(s(ic))+di*log(s(id)))
         else
            s(icd) = ci*s(ic)+di*s(id)
         endif
         ro(icd) = exp(ci*log(ro(ic))+di*log(ro(id)))

         vlnf(icd) = ci*vlnf(ic)+di*vlnf(id)
         vlnt(icd) = ci*vlnt(ic)+di*vlnt(id)
         vp(icd) = exp(ci*log(vp(ic))+di*log(vp(id)))
         ve(icd) = exp(ci*log(ve(ic))+di*log(ve(id)))
         if (vs(ic).gt.0.d0.and.vs(id).gt.0.d0) then
            vs(icd) = exp(ci*log(vs(ic))+di*log(vs(id)))
         else
            s(icd) = ci*vs(ic)+di*vs(id)
         endif
         vro(icd) = exp(ci*log(vro(ic))+di*log(vro(id)))

         tau(icd) = ci*tau(ic)+di*tau(id)
         vmueinv(icd) = ci*vmueinv(ic)+di*vmueinv(id)
         venucl(icd) = ci*venucl(ic)+di*venucl(id)
         if (vpvisc(ic).gt.0.d0.and.vpvisc(id).gt.0.d0) then
            vpvisc(icd) = exp(ci*log(vpvisc(ic))+di*log(vpvisc(id)))
         else
            vpvisc(icd) = ci*vpvisc(ic)+di*vpvisc(id)
         endif
         exposure(icd) = exposure(iab)

***   conservative laws (mass fraction conserved)
         if (xsp(ixb,io16).ne.xsp(ixd,io16)) then
            do l = 1,nsp
c..   default
               xsp(iab,l) = xai*xsp(ixa,l)+xbi*xsp(ixb,l)
               xsp(icd,l) = xci*xsp(ixc,l)+xdi*xsp(ixd,l)
c..   mesh2
c               xsp(iab,l) = xsp(ixa,l)
c               xsp(icd,l) = xsp(iab,l)
c..   mesh4
c               xsp(iab,l) = xsp(ixb,l)
c               xsp(icd,l) = xsp(ixd,l)
c..   mesh3
c               if (crz(ia).eq.crz(ixb)) then
c                  xsp(iab,l) = xai*xsp(ixa,l)+xbi*xsp(ixb,l)
c               else
c                  xsp(iab,l) = xsp(ixa,l)
c               endif
c               if (crz(ic).eq.crz(ixd)) then
c                  xsp(icd,l) = xci*xsp(ixc,l)+xdi*xsp(ixd,l)
c               else
c                  xsp(icd,l) = xsp(ixa,l)
c               endif
c..   mesh1
c               if (crz(il0).eq.crz(ir)) then
c                  xsp(iab,l) = xai*xsp(ixa,l)+xbi*xsp(ixb,l)
c                  xsp(icd,l) = xci*xsp(ixc,l)+xdi*xsp(ixd,l)
c               else
c                  xsp(iab,l) = xsp(ixa,l)
c                  xsp(icd,l) = xsp(ixa,l)
c               endif
c..   mesh5
c              if (crz(il0).eq.crz(ir)) then
c                 xsp(iab,l) = xai*xsp(ixa,l)+xbi*xsp(ixb,l)
c                 xsp(icd,l) = xci*xsp(ixc,l)+xdi*xsp(ixd,l)
c              else
c                 xsp(iab,l) = xsp(ixb,l)
c                 xsp(icd,l) = xsp(ixd,l)
c              endif
            enddo
         endif

         if (iaccr.gt.0) then
            facc(iab) = ai*facc(ia)+bi*facc(ib)
            macc(iab) = ai*macc(ia)+bi*macc(ib)
            vfacc(iab) = ai*vfacc(ia)+bi*vfacc(ib)
            dfaccdr(iab) = ai*dfaccdr(ia)+bi*dfaccdr(ib)
            vacc(iab) = ai*vacc(ia)+bi*vacc(ib)
            tacc(iab) = ai*tacc(ia)+bi*tacc(ib)

            facc(icd) = ci*facc(ic)+di*facc(id)
            macc(icd) = ci*macc(ic)+di*macc(id)
            vfacc(icd) = ci*vfacc(ic)+di*vfacc(id)
            dfaccdr(icd) = ci*dfaccdr(ic)+di*dfaccdr(id)
            vacc(icd) = ci*vacc(ic)+di*vacc(id)
            tacc(icd) = ci*tacc(ic)+di*tacc(id)
         endif
      endif



*__________________________________________________________________
***
***               TREATMENT OF MASS LOSS/ACCRETION
***
***   change of independent variable in case of accretion/mass-loss
***   and interpolation of old variables at the new mesh point qmr
*------------------------------------------------------------------


      if (ls.eq.4) then
         ia = ir-1
         ib = ir
         ic = ir-1
         id = ir
         bi = (qmr(il)-vqmr(ir-1))/(vqmr(ir)-vqmr(ir-1))
         ai = 1.d0-bi
***   special treatment for centered variables
         midmr = 0.5d0*(qmr(il)+qmr(il+1))
         vmidmr = 0.5d0*(vqmr(ia)+vqmr(ir))
         if (midmr.lt.vmidmr) then
            ic = ia-1
            id = ir-1
            goto 20
         endif
         if (ir.lt.nmod) then
            vmidmr = 0.5d0*(vqmr(ir)+vqmr(ir+1))
            if (midmr.gt.vmidmr) then
               do k = ir+1,nmod1
                  vmidmr = 0.5d0*(vqmr(k)+vqmr(k+1))
                  if (vmidmr.gt.midmr) then
                     id = k
                     ic = id-1
                     goto 20
                  endif
               enddo
            endif
         else
            id = nmod1
            ic = nmod1
         endif
 20      di = (qmr(il)+qmr(il+1)-vqmr(ic)-vqmr(ic+1))/(vqmr(id+1)-
     &        vqmr(ic))
         ci = 1.d0-di

***   centered variable
         lnf(il) = (ci*lnf(ic)+di*lnf(id))
         lnt(il) = (ci*lnt(ic)+di*lnt(id))
         p(il) = exp(ci*log(p(ic))+di*log(p(id)))
         e(il) = exp(ci*log(e(ic))+di*log(e(id)))
         s(il) = exp(ci*log(s(ic))+di*log(s(id)))
         ro(il) = exp(ci*log(ro(ic))+di*log(ro(id)))

         vlnf(il) = (ci*vlnf(ic)+di*vlnf(id))
         vlnt(il) = (ci*vlnt(ic)+di*vlnt(id))
         vp(il) = exp(ci*log(vp(ic))+di*log(vp(id)))
         ve(il) = exp(ci*log(ve(ic))+di*log(ve(id)))
         vs(il) = exp(ci*log(vs(ic))+di*log(vs(id)))
         vro(il) = exp(ci*log(vro(ic))+di*log(vro(id)))

         tau(il) = (ci*tau(ic)+di*tau(id))
         vmueinv(il) = (ci*vmueinv(ic)+di*vmueinv(id))
         venucl(il) = (ci*venucl(ic)+di*venucl(id))
         if (vpvisc(ic).gt.0.d0.and.vpvisc(id).gt.0.d0) then
            vpvisc(il) = exp(ci*log(vpvisc(ic))+di*log(vpvisc(id)))
         else
            vpvisc(il) = ci*vpvisc(ic)+di*vpvisc(id)
         endif
         exposure(il) = (ci*exposure(ic)+di*exposure(id))

***   variables defined at interface
         u(il) = ai*u(ia)+bi*u(ib)
         r(il) = (ai*r(ia)**3+bi*r(ib)**3)**pw13
         lum(il) = ai*lum(ia)+bi*lum(ib)

         vu(il) = ai*vu(ia)+bi*vu(ib)
         vr(il) = (ai*vr(ia)**3+bi*vr(ib)**3)**pw13
         vlum(il) = ai*vlum(ia)+bi*vlum(ib)

         vsconv(il) = ai*vsconv(ia)+bi*vsconv(ib)
         vhconv(il) = ai*vhconv(ia)+bi*vhconv(ib)

         vqmr(ia) = qmr(il)

      endif


      return
      end


************************************************************************

      SUBROUTINE kappa (tk,rok,muiinvk,x,ksh,kap1,kapr1,kapt1,kapx1
     &     ,kapy1,opamol,iter,iopa,iopaf,error)

************************************************************************
* Determine the total opacity coefficient and its derivatives          *
*     kap1 = kappa
*     kapr1 = d kappa/dro
*     kapt1 = d kappa/dT
*     kapx1 = d kappa/dX  (X = H mass fraction)
*     kapy1 = d kappa/dY  (Y = He mass fraction)
*     zkint : metallicity of the initial model
*     (must be specified in the parameter card)
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.diff'
c      include 'evolcom.kap'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.opa'

      integer mlth,nlth,itlth
      
      double precision rklth,tklth,opaclth
C     Warning: dimensions of opacity tables change with solar comp.
c..   Young+18 (TD: 30/01/2020)
c     common /lthopac/ opaclth(19,81,155),tklth(81),rklth(19),
c    &     mlth,nlth,itlth
       common /lthopac/ opaclth(19,85,155),tklth(85),rklth(19),
     &     mlth,nlth,itlth
c#ifdef GRID
c      parameter (mlth = 19,nlth = 63,itlth = 104)
c      common /lthopac/ rklth(19),tklth(63),opaclth(19,63,104)
c#elif AS09
c      parameter (mlth = 19,nlth = 85,itlth = 155)
c      common /lthopac/ rklth(19),tklth(85),opaclth(19,85,155)      
c#else
c      parameter (mlth = 19,nlth = 85,itlth = 120)
c      common /lthopac/ rklth(19),tklth(85),opaclth(19,85,120)
c#endif

      integer mt,nt,kdel,kxdel,kk(4)
      integer ir,it,j,iter
      integer ii,iim,iin
      integer error,opamol,ksh,iopa,iopaf

      double precision tk,rok,muiinvk
      double precision x(nsp)
      double precision kap1,kapr1,kapt1,kapx1,kapy1
      double precision dxy,kap1x,kap1y
      double precision ck(4),ckr(4),ckt(4)
      double precision tik,rotik,til,rotil,xx,xy,xc,xn,xo,xno,xsolc,
     &     xsoln,xsolo,xz,topa,uopa,d1opa,d2opa,xsi1,xsi2,xsi11,xsi12,
     &     yyopa1,yyopa2,yyopa3,yyopa4,kap4,kapr4,kapt4
      double precision yopa(4),y1opa(4),y2opa(4),y12opa(4)

      double precision opact,dopact,dopacr,dopactd,ftredge,fedge,fzedge
      double precision xspsol,zsol

      character*6 refsolar

      common /extopac/ opact,dopact,dopacr,dopactd,fedge,ftredge,fzedge
      common /solar/ xspsol(nis+6),zsol,refsolar

      mt = 0
      nt = 0
      dxy = 1.d-2
      tik = tk*1.d-3
      rotik = log10(rok/((tik*1.d-3)**3))
      rotik = max(rotik,-9.d0)
      rotik = min(rotik,1.d1)
      til = tk*1.d-6
      rotil = rok/til**3
      rotil = max(rotil,1.d-9)
      rotil = min(rotil,1.d10)
      xx = x(ih1)+x(ih2)
      xy = x(ihe3)+x(ihe4)
      xsolc = zkint*(xspsol(ic12)+xspsol(ic13)+xspsol(ic14))/zsol
      xsoln = zkint*(xspsol(in13)+xspsol(in14)+xspsol(in15))/zsol
      xsolo = zkint*(xspsol(io16)+xspsol(io17)+xspsol(io18))/zsol
      xc = x(ic12)+x(ic13)+x(ic14)-xsolc
      xn = x(in13)+x(in14)+x(in15)-xsoln
      xo = x(io16)+x(io17)+x(io18)-xsolo
      xc = max(0.d0,xc)
      xn = max(0.d0,xn)
      xo = max(0.d0,xo)
***   If solar calibration with diffusion, then metallicity should correspond to the actual
***   metallicity of the model
c      if (totm.eq.1.d0) then
c         xz =  1.d0-xx-xy
c      else
         xz = zkint
c      endif

*------------------------------------------------------------------------
***   calculation of the opacity coefficient for low-T regions and
***   variable mixture of C and O (including molecules of H20, OH, H2, CO
***   C2 CN and N2)
***   Only if the C/O ratio (in number) is greater than 0.8 (i.e. C/O in
***   mass greater than 0.6)
***   WARNING : below 1700K, dust not taken into account
*------------------------------------------------------------------------

      if (tk.le.5.d3.and.opamol.eq.1.and.x(ic12).gt.0.6d0*x(io16)) then
         if (iopa.eq.0) then 
            iopa = ksh
            if (iter.eq.1) write (nout,100) iopa
         endif
         if (tk.lt.1.7d3) write (nout,200)
         call opa_co (tk,rok,muiinvk,x,ksh,kap1,kapr1,kapt1,error)
         kapx1 = 0.d0
         kapy1 = 0.d0

         return

      endif
      
*------------------------------------------------------------------------
***   calculation of the opacity coefficient for H-rich and low-T regions
***   low temperature opacities : Ferguson et al. (2005)
*------------------------------------------------------------------------
      
      if (tk.le.tlhkap) then
         if (iopaf.eq.0) then 
            iopaf = ksh
            if (iter.eq.1.and.model.eq.modeli) write (nout,250) iopaf
         endif
         do ir = 1,mlth
            if (rklth(ir).gt.rotik) goto 10
         enddo
 10      mt = ir-1
         do it = 1,nlth
            if (tklth(it).le.tk) goto 20
         enddo
 20      nt = it-1
         if (mt.lt.1) mt = 1
         if (nt.lt.1) nt = 1
         if (mt.gt.(mlth-1)) mt = mlth-1
         if (nt.gt.(nlth-1)) nt = nlth-1
         topa = (rotik-rklth(mt))/(rklth(mt+1)-rklth(mt))
         uopa = (tk-tklth(nt))/(tklth(nt+1)-tklth(nt))
         d1opa = rklth(mt+1)-rklth(mt)
         d2opa = tklth(nt+1)-tklth(nt)
         
         if (refsolar.eq.'GN93') then
            
            if (xz.le.1.d-5) then
               kdel = 14
               xsi1 = xz/1.d-5
            elseif (xz.le.3.d-5) then
               kdel = 13
               xsi1 = (xz-1.d-5)/3.d-5
            elseif (xz.le.1.d-4) then
               kdel = 12
               xsi1 = (xz-3.d-5)/7.d-5
            elseif (xz.le.3.d-4) then
               kdel = 11
               xsi1 = (xz-1.d-4)/3.d-4
            elseif (xz.le.1.d-3) then
               kdel = 10
               xsi1 = (xz-3.d-4)/7.d-4
            elseif (xz.le.2.d-3) then
               kdel = 9
               xsi1 = (xz-1.d-3)*1.d3
            elseif (xz.le.4.d-3) then
               kdel = 8
               xsi1 = (xz-2.d-3)*5.d2
            elseif (xz.le.1.d-2) then
               kdel = 7
               xsi1 = (xz-4.d-3)/6.d-3
            elseif (xz.le.2.d-2) then
               kdel = 6
               xsi1 = (xz-1.d-2)*1.d2
            elseif (xz.le.3.d-2) then
               kdel = 5
               xsi1 = (xz-2.d-2)*1.d2
            elseif (xz.le.4.d-2) then
               kdel = 4
               xsi1 = (xz-3.d-2)*1.d2
            elseif (xz.le.6.d-2) then
               kdel = 3
               xsi1 = (xz-4.d-2)*5.d1
            elseif (xz.le.8.d-2) then
               kdel = 2
               xsi1 = (xz-6.d-2)*5.d1
            elseif (xz.le.1.d-1) then
               kdel = 1
               xsi1 = (xz-8.d-2)*5.d1
            elseif (xz.gt.1.d-1) then
               kdel = 1
               xsi1 = 1.d0
            endif
            
            if (xx.le.0.1d0) then
               kxdel = 0
               xsi2 = 1.d0
            elseif (xx.le.0.2d0) then
               kxdel = 15
               xsi2 = (xx-0.1d0)*1.d1
            elseif (xx.le.0.35d0) then
               kxdel = 30
               xsi2 = (xx-0.2d0)/1.5d-1
            elseif (xx.le.0.5d0) then
               kxdel = 45
               xsi2 = (xx-0.35d0)/1.5d-1
            elseif (xx.le.0.7d0) then
               kxdel = 60
               xsi2 = (xx-0.5d0)*5.d0
            elseif (xx.le.0.8d0) then
               kxdel = 75
               xsi2 = (xx-0.7d0)*1.d1
            elseif (xx.le.0.9d0) then
               kxdel = 90
               xsi2 = (xx-0.8d0)*1.d1
            elseif (xx.gt.0.9d0) then
               kxdel = 105
               xsi2 = 0.d0
            endif

            kk(1) = kxdel+kdel
            kk(2) = kk(1)+1
            kk(3) = kk(1)+15
            kk(4) = kk(3)+1
         else if
     $           (refsolar.eq.'AGS05'.or.refsolar.eq.'AGSS09'.or.
     $           refsolar.eq.'AY18') then
            if (xz.le.1.d-5) then
               kdel = 15
               xsi1 = xz*1.d5
            elseif (xz.le.3.d-5) then
               kdel = 14
               xsi1 = (xz-1.d-5)/3.d-5
            elseif (xz.le.1.d-4) then
               kdel = 13
               xsi1 = (xz-3.d-5)/7.d-5
            elseif (xz.le.3.d-4) then
               kdel = 12
               xsi1 = (xz-1.d-4)/3.d-4
            elseif (xz.le.1.d-3) then
               kdel = 11
               xsi1 = (xz-3.d-4)/7.d-4
            elseif (xz.le.2.d-3) then
               kdel = 10
               xsi1 = (xz-1.d-3)*1.d3
            elseif (xz.le.4.d-3) then
               kdel = 9
               xsi1 = (xz-2.d-3)*5.d2
            elseif (xz.le.1.d-2) then
               kdel = 8
               xsi1 = (xz-4.d-3)/6.d-3
            elseif (xz.le.2.d-2) then
               kdel = 7
               xsi1 = (xz-1.d-2)*1.d2
            elseif (xz.le.3.d-2) then
               kdel = 6
               xsi1 = (xz-2.d-2)*1.d2
            elseif (xz.le.4.d-2) then
               kdel = 5
               xsi1 = (xz-3.d-2)*1.d2
            elseif (xz.le.5.d-2) then
               kdel = 4
               xsi1 = (xz-4.d-2)*1.d2
            elseif (xz.le.6.d-2) then
               kdel = 3
               xsi1 = (xz-5.d-2)*5.d1
            elseif (xz.le.8.d-2) then
               kdel = 2
               xsi1 = (xz-6.d-2)*5.d1
            elseif (xz.le.1.d-1) then
               kdel = 1
               xsi1 = (xz-8.d-2)*5.d1
            elseif (xz.gt.1.d-1) then
               kdel = 1
               xsi1 = 1.d0
            endif
          
            if (xx.le.0.1d0) then
               kxdel = 0
               xsi2 = 1.d0
            elseif (xx.le.0.2d0) then
               kxdel = 16
               xsi2 = (xx-0.1d0)*1.d1
            elseif (xx.le.0.35d0) then
               kxdel = 32
               xsi2 = (xx-0.2d0)/1.5d-1
            elseif (xx.le.0.5d0) then
               kxdel = 48
               xsi2 = (xx-0.35d0)/1.5d-1
            elseif (xx.le.0.7d0) then
               kxdel = 64
               xsi2 = (xx-0.5d0)*5.d0
            elseif (xx.le.0.8d0) then
               kxdel = 80
               xsi2 = (xx-0.7d0)*1.d1
            elseif (xx.le.0.9d0) then
               kxdel = 96
               xsi2 = (xx-0.8d0)*1.d1
            elseif (xx.gt.0.9d0) then
               kxdel = 112
               xsi2 = 0.d0
            endif

            kk(1) = kxdel+kdel
            kk(2) = kk(1)+1
            kk(3) = kk(1)+16
            kk(4) = kk(3)+1
            
         elseif (refsolar.eq.'GRID') then
C Warning: dimension of opacitiy tables changes for GRID version 
C Indices for Asplund/Cunha solar comp.

            if (xz.le.1.d-4) then
               kdel = 12
               xsi1 = xz*1.d4
            elseif (xz.le.3.d-4) then
               kdel = 11
               xsi1 = (xz-1.d-4)/2.d-4
            elseif (xz.le.1.d-3) then
               kdel = 10
               xsi1 = (xz-3.d-4)/7.d-4
            elseif (xz.le.2.d-3) then
               kdel = 9
               xsi1 = (xz-1.d-3)*1.d3
            elseif (xz.le.4.d-3) then
               kdel = 8
               xsi1 = (xz-2.d-3)*5.d2
            elseif (xz.le.1.d-2) then
               kdel = 7
               xsi1 = (xz-4.d-3)/6.d-3
            elseif (xz.le.2.d-2) then
               kdel = 6
               xsi1 = (xz-1.d-2)*1.d2
            elseif (xz.le.3.d-2) then
               kdel = 5
               xsi1 = (xz-2.d-2)*1.d2
            elseif (xz.le.4.d-2) then
               kdel = 4
               xsi1 = (xz-3.d-2)*1.d2
            elseif (xz.le.6.d-2) then
               kdel = 3
               xsi1 = (xz-4.d-2)*5.d1
            elseif (xz.le.8.d-2) then
               kdel = 2
               xsi1 = (xz-6.d-2)*5.d1
            elseif (xz.le.1.d-1) then
               kdel = 1
               xsi1 = (xz-8.d-2)*5.d1
            elseif (xz.gt.1.d-1) then
               kdel = 1
               xsi1 = 1.d0
            endif
            
            if (xx.le.0.1d0) then
               kxdel = 0
               xsi2 = 1.d0
            elseif (xx.le.0.2d0) then
               kxdel = 13
               xsi2 = (xx-0.1d0)*1.d1
            elseif (xx.le.0.35d0) then
               kxdel = 26
               xsi2 = (xx-0.2d0)/1.5d-1
            elseif (xx.le.0.5d0) then
               kxdel = 39
               xsi2 = (xx-0.35d0)/1.5d-1
            elseif (xx.le.0.7d0) then
               kxdel = 52
               xsi2 = (xx-0.5d0)*5.d0
            elseif (xx.le.0.8d0) then
               kxdel = 65
               xsi2 = (xx-0.7d0)*1.d1
            elseif (xx.le.0.9d0) then
               kxdel = 78
               xsi2 = (xx-0.8d0)*1.d1
            elseif (xx.gt.0.9d0) then
               kxdel = 91         
               xsi2 = 0.d0
            endif
            
            
            kk(1) = kxdel+kdel
            kk(2) = kk(1)+1
            kk(3) = kk(1)+13
            kk(4) = kk(3)+1
         endif

         do j = 1,4
            ck(j) = 0.d0
            ckr(j) = 0.d0
            ckt(j) = 0.d0
            if (kk(j).ne.0) then
               yyopa1 = opaclth(mt,nt,kk(j))
               yyopa2 = opaclth(mt+1,nt,kk(j))
               yyopa3 = opaclth(mt+1,nt+1,kk(j))
               yyopa4 = opaclth(mt,nt+1,kk(j))
               call interlin (topa,uopa,d1opa,d2opa,yyopa1,yyopa2,
     &              yyopa3,yyopa4,ck(j),ckr(j),ckt(j))
            endif
         enddo
         xsi12 = 1.d0-xsi2
         xsi11 = 1.d0-xsi1
         kap4 = xsi12*xsi11*ck(2)+xsi2*xsi11*ck(4)+xsi1*xsi2*ck(3)+
     &        xsi12*xsi1*ck(1)
         kapr4 = xsi12*xsi11*ckr(2)+xsi2*xsi11*ckr(4)+xsi1*xsi2*ckr(3)+
     &        xsi12*xsi1*ckr(1)
         kapt4 = xsi12*xsi11*ckt(2)+xsi2*xsi11*ckt(4)+xsi1*xsi2*ckt(3)+
     &        xsi12*xsi1*ckt(1)
         kap1 = 10.d0**kap4
         kapr1 = kap1*kapr4/rok
         kapt1 = kap1*(-3.d0*kapr4/tk+kapt4*2.302585093d0)
         kapx1 = 0.d0
         kapy1 = 0.d0
         return
         
*-------------------------------------------------
*** OPAL opacities for different Z and CO mixtures
*-------------------------------------------------

      else
         call opac (xz,xx,xc,xo,til,rotil,error)

         if (opact.gt.1.d15.or.error.gt.0) then
            write (nout,300) ksh,xz,xx,xc,xo,tk,rok,rotik
            error = 81
            return
         endif

         kap1 = 10.d0**opact
         kapr1 = kap1*dopacr/rok
         kapt1 = kap1*(-3.d0*dopacr+dopact)/tk
         kapx1 = 0.d0
         kapy1 = 0.d0

      endif

      kapx1 = 0.d0
      kapy1 = 0.d0
      if (rotation) then
         xno = xo
         call opac (xz,xx+xx*dxy,xc,xno,til,rotil,error)
         if (error.gt.0) then
            write (nout,400) ksh,xz,xx+xx*dxy,xc,xno,tk,rok,rotik
            error = 81
            return
         endif
         kap1x = 10.d0**opact
         kapx1 = (kap1x-kap1)/(kap1*dxy)
         call opac (xz,xx-xx*dxy,xc,xno,til,rotil,error)
         if (error.gt.0) then
            write (nout,400) ksh,xz,xx-xx*dxy,xc,xno,tk,rok,rotik
            error = 81
            return
         endif
         kap1y = 10.d0**opact
         kapy1 = xy*(kap1y-kap1)/(kap1*xx*dxy)
      endif

 100  format (' molecular opacity computed from shell #',i4)
 200  format (2x,'WARNING : T < 1700K ! dust not taken into account ',
     &        'in opacity')
 250  format (' Fergusson opacity computed from shell #',i4)
 300  format (3x,'KAPPA : pb table 3 - O rich : #',i4,', xZ =',1pe10.3,
     &     ', xH = ',1pe9.3,', xC = ',1pe9.3,', xO = ',1pe9.3,', T = ',
     &     1pe9.3,', rho = ',1pe9.3,', R = ',0pf8.3)
 400  format (3x,'KAPPA rotation : pb table 3 - O rich : #',i4,', xZ =',
     &     1pe9.3,', xH = ',1pe9.3,', xC = ',1pe9.3,', xO = ',1pe9.3,
     &     ', rho = ',1pe9.3,', R = ',0pf8.3)

      return
      end
      SUBROUTINE mchange

************************************************************************
*     Modify the mass distribution due to mass-loss and/or accretion   *
*     WARNING : never use t it is not defined yet, use lnt !
*     mlp = 1,2   : Reimers
*     mlp = 3,4   : de Jager
*     mlp = 5,6   : Vassiliadis & Wood without delaying the onset of super-wind
*     mlp = 55,56 : Vassiliadis & Wood original prescription
*     mlp = 7,8   : Blocker
*     mlp = 9,10  : Arndt
*     mlp = 11,12 : Schaller et al.
*     mlp = 13,14 : Chiosi
*     mlp = 15,16 : Vink et al.
*     mlp = 17,18 : massrate, specified in starevol.par
*     mlp = 19,20 : Crowther
*     mlp = 21,22 : Van loon et al (2005)
*     mlp = 23,24 : Cranmer & Saar (2011)
*     mlp = 25,26 : Graefener (2021)
*     mlp = 27,28 : Sanders & Vink (2022) -> VMS stars
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.data'
      include 'evolcom.eng'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.rot'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.transp'
      include 'evolcom.var'

      integer imenv,ii1,jj,jj1,nn1,nold,mm
      integer ibot,itop,iqmrbot,iqmrtop
      integer mlp1,nmod0,ibind,mlsh
      integer i,ii,ij,j,k,nn,j1
      integer itmax,ishockb,ishockt
      integer nconvsh,ienvsh
      integer icrazy

      logical arndt,blocker,wood,limit!,zscaling
      character prescrip*120,zdep*50

      double precision xmm,momlim,dmsmec,rapom,omegacrit
      double precision meini,Leini,omsurfini,dLrad,omegamax
      double precision mgammaedd
      double precision dlogm,dma
      double precision mominnew,mominold
      double precision plog,pcomp,mdlog,vesc,dmsmax,dmsreim
      double precision qmr,vqmr,xtot
      double precision rco,tteff,lleff,xmlos,ymlos,totmv,totmd,
     &     xxm,dxxm,ddm,vddm,vtotma
      double precision dmaccrg,fedd,ebind
      double precision sigmae,gammae,brackro,teffjump1,teffjump2,
     &     ratio,logdms,ff,fxm
      double precision pxsp(nsh,nsp),vdm(nsh),vm(nsh),vvr(nsh)
      double precision malpha,mexp,msigmae,mrom,mgamma,mfac,dms0,mfac1
      double precision Lmax,tmax,mackmax,enucmax
      double precision menvconv,dmenvconv,mcore
      double precision prot,FeH,logg
      double precision Bcrit,fstar,Bequi
      double precision qsig, HeH
      double precision ledd
      double precision geffrot,Rp2
      double precision kap,kapm
      double precision disctime
      double precision alfa,logdms1,logdms2,dms1,dms2,dT
      double precision tau1,tau2,taus,gammaeff,vinf,gammarot
     $     ,vesceff,dmsthick
      double precision Dclump,etawind
      double precision acoeff, cbd, logmdotoff, gammaeb,logdmsvms

      common /disclocking/ disctime
      common /hydrodyn/ Lmax,tmax,mackmax,enucmax,itmax,ishockb,ishockt
      common /meshconv/ menvconv,dmenvconv,mcore,nconvsh,ienvsh
      common /newmesh/ qmr(nsh),vqmr(nsh),iqmrbot,iqmrtop
      common /metal/ FeH
      common /boreas_Mloss/ Bcrit,fstar,Bequi
      common /geffective/ geffrot(nsh),Rp2(nsh)
      common /opa/ kap(nsh),kapm(nsh)

      if (abs(mlp0).gt.49) then
         limit = .true.
      else
         limit = .false.
      endif

*____________________
***   Initializations
*--------------------
      
      zeff = 1.d0-xsp(neff,ih1)-xsp(neff,ih2)-xsp(neff,ihe3)-
     &     xsp(neff,ihe4)
      zeff = max(zeff,5.d-5)
      ibind = 0
      iqmrbot = nmod
      iqmrtop = nmod
      iaccbot = nmod
      iacctop = 1
      mlsh = nmod
      m(1) = 0.d0
      dm(nmod) = 0.d0
      forall (i = 1:nmod)
         omi(i) = 0.d0
         psi(i) = 0.d0
         vdm(i) = dm(i)
         vm(i) = m(i)
      end forall
      dlogm = 0.d0
      dms = 0.d0
      Bequi = 0.d0

***   check if surface layers are still gravitationally bound
      ebind = 0.d0
      do i = nmod1,2,-1
         ebind = ebind+(e(i)-g*m(i)/r(i))*dm(i)
      enddo
      if (ebind.gt.0.d0) then
         write (nout,*) 'WARNING : star not gravitationnaly bound !!!',
     &        ebind
c         totmv = totm*msun
c         nold = nmod
c         dms = (m(nmod)-m(i))/dtn
c         totmd = m(i)
c         ibind = 1
c         goto 20
      endif


***   if no mass loss, define interpolation coefficients
      if (dmaccr.eq.0.d0.and.((nphase.lt.2.and.mlp.ne.23)
     $    .or.mlp.eq.0)) then
            if (numeric.eq.2.or.numeric.eq.3) then
c..   first order accuracy in the spatial derivatives
               forall (i = 1:nmod)
               wi(i) = 0.5d0
               wj(i) = 0.5d0
               end forall
            else
c..   second order accuracy in the spatial deriatives
               wi(1) = 0.5d0
               wj(1) = 0.5d0
               forall (i = 2:nmod1)
               wi(i) = dm(i)/(dm(i)+dm(i-1))
               wj(i) = 1.d0-wi(i)
               end forall
               wi(nmod) = 0.5d0
               wj(nmod) = 0.5d0
            endif
            return
!!         endif
      endif


*__________________________________
***   choice of the mass loss regim
*----------------------------------

      if (mlp.gt.0) then

         if (irotbin.eq.1.and.(mlp.eq.23.or.mlp.eq.24).and.
     &        (time/sec.gt.disctime)) then
            prot = pi/(vomega(nmod)*43200.d0)
            logg =  log10(g*m(neff)/(reff*reff))
            call boreas (totm,solr,soll,prot,gmr(nmod),0.d0,dms,teff,
     &           icrazy)
            if (icrazy.gt.0) mlp = 01
            print *,'Cranmer mass loss',dms,'solar masses / year '
            print *,'with gmr =',gmr(nmod)
         endif

         if (tmax.gt.5.d6.or.mlp.eq.23.or.mlp.eq.24) then
            nmod0 = 0
            if (mod(mlp,2).ne.0.and.dm(nmod1).le.1.d-15*m(nmod1))
     &           nmod0 = nmod
            do k = nmod1,1,-1
               if (dm(k).le.1.d-15*m(k).and.nmod0.ne.0.and.nmod0-1.eq.k)
     &              nmod0 = k
               if (lnt(k).gt.log(5.d6)) goto 10
            enddo
 10         imenv = max(k+1,2)
            
            mlsh = imenv
            mlp1 = mlp
            rco = (xsp(neff,ic12)+xsp(neff,ic13)+xsp(neff,ic14))/
     &           (xsp(neff,io16)+xsp(neff,io17)+xsp(neff,io18))
            tteff = log10(teff)
            lleff = log10(soll)
            if (tteff.le.3.2d0.or.tteff.ge.4.7d0.or.lleff.le.2.6d0.or.
     &           lleff.ge.6.6d0) then
               if (mlp.eq.3) mlp = 1
               if (mlp.eq.4) mlp = 2
            endif
***   change in mass-loss rate beyond the main sequence for B type stars
            
            if (mtini.ge.7.d0.and.mtini.le.12.d0.and.nphase.gt.2) then
               if (mlp.eq.15) mlp = 1
               if (mlp.eq.16) mlp = 2
            endif

***	switch to de Jager et al 1988 when Vink et al. 2000,2001 is not valid

c     if (tteff.le.3.9d0.and.(mlp.eq.15.or.mlp.eq.16)) then
            if (tteff.le.4.0969d0.and.(mlp.eq.15.or.mlp.eq.16)) then
               if (tteff.gt.3.7d0) then
                  write (nout,200)
                  if (mlp.eq.15) then
                     mlp = 3
                  else
                     mlp = 4
                  endif
               else
                  write (nout,300)
                  mlp = 19
               endif
            endif

***   change in mass-loss rate during AGB evolution

            if (agbphase.or.superagb) then
               arndt = totm.lt.1.8d0.and.rco.gt.1.2d0.and.teff.lt.3.d3
               blocker = totm.ge.1.d0.and.totm.le.7.d0
               if (mlp.eq.7.and..not.blocker) mlp = 1
               if (mlp.eq.8.and..not.blocker) mlp = 2
               if (mlp.eq.1.or.(mlp.eq.3.and..not.superagb)) mlp = 5
               if (mlp.eq.2.or.(mlp.eq.4.and..not.superagb)) mlp = 6
               if ((mlp.eq.5.or.mlp.eq.55.or.mlp.eq.7).and.arndt) then
                  mlp = 9
               endif
               if ((mlp.eq.6.or.mlp.eq.56.or.mlp.eq.8).and.arndt) then
                  mlp = 10
               endif
               if (mlp.ne.mlp1) then
                  write (nout,400) mlp
                  write (90,400) mlp
               endif
            endif
            wood = mlp.eq.5.or.mlp.eq.6.or.mlp.eq.55.or.mlp.eq.56


               
*_________________________________
***   mass-loss rate prescriptions
*---------------------------------

***   Reimers (1975) prescription (MS and RGB phases)

            if (mlp.eq.1.or.mlp.eq.2) then
               dms = min(etapar*3.98d-13*soll*solr/totm,1.d-4)
               print *,'Reimers mass loss',dms,'solar masses / year'
            endif


***   de Jager etal. (1988 A&AS,72,259) (RGB phase + massive stars)

            if (mlp.eq.3.or.mlp.eq.4) then
               if (tteff.gt.3.7d0) then
                  dms = 0.d0
                  xmlos = (tteff-4.05d0)/0.75d0
                  ymlos = (lleff-4.6d0)/2.1d0
                  xmlos = max(xmlos,-1.d0)
                  xmlos = min(xmlos,1.d0)
                  ymlos = max(ymlos,-1.d0)
                  ymlos = min(ymlos,1.d0)
                  do nn = 1,6
                     nn1 = nn-1
                     do ii = 1,nn
                        ii1 = ii-1
                        if (nn1-ii1.ne.5) then
                           jj1 = nn1-ii1
                           jj = jj1+1
                           dms = dms+amlos(ii,jj)*cos(dble(ii1)*
     &                          acos(xmlos))*cos(dble(jj1)*acos(ymlos))
                        endif
                     enddo
                  enddo
               else
*    Mass loss rate for RSG based on Sylvester+1998 and van Loon+99 (LMC)
                  dms = -(1.7d0*lleff-13.83d0)
               endif
               dms = 10.d0**(-dms)
               if (limit) dms = min(dms,1.d-4)
            endif

***   Vassiliadis & Wood (1993, ApJ 413, 641) prescription (AGB phase)

            if (wood) then
               plog = -2.07d0+1.94d0*log10(solr)-0.9d0*log10(totm)
               pcomp = 10.d0**plog
               if (totm.le.2.5d0.and.(mlp.eq.55.or.mlp.eq.56)) then
c..   modification for M > 2.5Mo to delay the onset of the superwind phase
                  mdlog = -11.4d0+0.0125d0*(pcomp-1.d2*(totm-2.5d0))
               else
                  mdlog = -11.4d0+0.0125d0*pcomp
               endif
c..   WARNING : limitation of mass loss rate to 1d-4 Mo/yr
               if (limit) mdlog = min(mdlog,-4.d0)
               dms = 10.d0**mdlog
               vesc = -13.5d0+0.056d0*pcomp
               vesc = vesc*1.d5
               vesc = max(3.d5,vesc)
               vesc = min(1.5d6,vesc)
c     dmsmax = lum(nmod)/(c*vesc)
               dmsmax = leff/(c*vesc)
               dmsmax = dmsmax*sec/msun
               dms = min(dms,dmsmax)
            endif

***   Blocker (1995, A&A 297, 727) prescription (AGB phase)

            if (mlp.eq.7.or.mlp.eq.8) then
               plog = -2.07d0+1.94d0*log10(solr)-0.9d0*log10(totm)
               pcomp = 10.d0**plog
               dmsreim = etapar*3.98d-13*soll*solr/totm
               if (pcomp.lt.1.d2) then
                  dms = dmsreim
               else
                  dms = 4.83d-9*soll**2.7d0*mtini**(-2.1d0)*dmsreim
               endif
c..   WARNING : limitation of mass loss rate to 1d-3 Mo/yr
               if (limit) dms = min(dms,1.d-3)
            endif

***   Arndt (1997, A&A 327, 614) prescription (low mass AGB)

            if (mlp.eq.9.or.mlp.eq.10) then
               dms = 17.158d0-8.26d0*tteff+1.53d0*lleff-2.88d0*
     &              log10(totm)
               dms = 10.d0**dms
            endif

***   Schaller et al. (1992) prescription (massive stars and WR phase)

            if (mlp.eq.11.or.mlp.eq.12) then
               dms = 4.d-8*totm**2.5d0
            endif

***   Chiosi (1981, A&A 93, 163) prescription (massive O stars)

            if (mlp.eq.13.or.mlp.eq.14) then
               fedd = abs(1.d0-1.31237d-5*lum(nmod)/m(nmod))
c     dms = 5.5d-14*soll**1.5d0*solr**2.25d0*fedd**(-1.75d0)*
c     &           totm**(-2.25d0)
               dms=2.6d-10*soll**0.72d0*(solr/(fedd*totm))**2.5d0
            endif

***   Vink et al. (2000) (massive stars with log Teff > 3.9) then
***   Nugis & Lamers (2000) for the WR phase , i.e. H_surf < 0.4 and
***   log Teff > 4.0

            if (mlp.eq.15.or.mlp.eq.16) then
               print *,'M loss for massive stars, Vink 2000'
               if (xsp(nmod1,ih1).le.0.4d0.and.tteff.gt.4.d0) then
                  print *,'M loss for massive stars : NL2000'
***   Nugis & Lamers (2000) for WR stars
                  write (nout,500)
                  logdms =-11.d0+1.29d0*log10(soll)+
     &                 1.73d0*log10(xsp(nmod,ihe4))+0.47d0*log10(zeff)
                  dms = 10**logdms
               else
***   Vink et al. 2001
c..   Lamers & Leitherer (1993), ApJ 412, p771, eq. 2
                  HeH = xsp(nmod1,ihe4)/(xsp(nmod1,ihe4)+xsp(nmod1,ih1))
                  if (teff.lt.3.d4) then
                     qsig = 0d0
                  elseif (teff.ge.3.d4.and.teff.lt.3.5d4) then
                     qsig = 0.5d0
                  else
                     qsig = 1.d0
                  endif
c                  sigmae = 0.401d0*(1.d0+qsig*HeH)/(1.d0+3.d0*HeH)
                  sigmae = 0.325d0
                  gammae = 7.66d-5*sigmae*soll/totm
                  brackro = -14.94d0+3.1857d0*gammae+0.85d0*
     &                 log10(zeff/0.019d0)
c     &                 log10(zeff/zsol)
                  teffjump1 = 1.d3*(61.2d0+2.59d0*brackro)
                  teffjump2 = 1.d3*(1.d2+6.d0*brackro)
                  print *,'teffjump1',teffjump1,'teffjump2',teffjump2
                  if (teffjump1.le.teffjump2) then
                     write (nout,'("unrealistic stellar parameters")')
                     stop 'mchange'
                  endif
c..   Smoothing the transitions following the eval_Vink_wind routine from MESA
                  if (teff.ge.27500.d0) then
                     ratio = 2.6d0
                     logdms1 = -6.697d0+2.194d0*log10(soll*1.d-5)
     &                    -1.313d0*log10(totm/30.d0)
     &                    -1.226d0*log10(ratio*0.5d0)
     &                    +9.33d-1*log10(teff/4.d4)
     &                    -10.92d0*(log10(teff/4.d4))**2
     &                    +0.85d0*log10(zeff/0.019d0)
                        dms = 10.d0**logdms1
                  else if (teff.le.22500.d0.and.teff.ge.12500.d0) then
                     if (teff.ge.teffjump2) then
                        ratio = 1.3d0
                        logdms2 = -6.688d0+2.210d0*log10(soll*1.d-5)
     &                       -1.339d0*log10(totm/30.d0)
     &                       -1.601d0*log10(ratio*0.5d0)
     &                       +1.07d0*log10(teff/2.d4)
     &                       +0.85d0*log10(zeff/0.019d0)
                        dms = 10.d0**logdms2
                     else
                        ratio = 0.7d0
                        logdms2 = -5.99d0+2.210d0*log10(soll*1.d-5)
     &                       -1.339d0*log10(totm/30.d0)
     &                       -1.601d0*log10(ratio*0.5d0)
     &                       +1.07d0*log10(teff/2.d4)
     &                       +0.85d0*log10(zeff/0.019d0)
                        dms = 10.d0**logdms2
                     endif
                  else if (teff.gt.22500.d0.and.teff.lt.27500.d0) then
                     alfa = (teff-22500.d0)/(5.d3)
                     logdms1 = -6.697d0+2.194d0*log10(soll*1.d-5)
     &                    -1.313d0*log10(totm/30.d0)
     &                    -1.226d0*log10(2.6d0*0.5d0)
     &                    +9.33d-1*log10(teff/4.d4)
     &                    -10.92d0*(log10(teff/4.d4))**2
     &                    +0.85d0*log10(zeff/0.019d0)
                     logdms2 = -6.688d0+2.210d0*log10(soll*1.d-5)
     &                    -1.339d0*log10(totm/30.d0)
     &                    -1.601d0*log10(1.3d0*0.5d0)
     &                    +1.07d0*log10(teff/2.d4)
     &                    +0.85d0*log10(zeff/0.019d0)
                     dms1 = 10.d0**logdms1
                     dms2 = 10.d0**logdms2
                     dms =  (1-alfa)*dms2 + alfa*dms1
                  endif
c$$$                   if (teff.le.teffjump1) then
c$$$                     if (teff.lt.teffjump2) then
c$$$                        ratio = 0.7d0
c$$$                        logdms = -5.99d0+2.210d0*log10(soll*1.d-5)
c$$$     &                       -1.339d0*log10(totm/30.d0)
c$$$     &                       -1.601d0*log10(ratio*0.5d0)
c$$$     &                       +1.07d0*log10(teff/2.d4)
c$$$     &                       +0.85d0*log10(zeff/0.019d0)
c$$$                     else
c$$$                        ratio = 1.3d0
c$$$                        logdms = -6.688d0+2.210d0*log10(soll*1.d-5)
c$$$     &                       -1.339d0*log10(totm/30.d0)
c$$$     &                       -1.601d0*log10(ratio*0.5d0)
c$$$     &                       +1.07d0*log10(teff/2.d4)
c$$$     &                       +0.85d0*log10(zeff/0.019d0)
c$$$                     endif
c$$$                  else
c$$$                     ratio = 2.6d0
c$$$                     logdms = -6.697d0+2.194d0*log10(soll*1.d-5)
c$$$     &                    -1.313d0*log10(totm/30.d0)
c$$$     &                    -1.226d0*log10(ratio*0.5d0)
c$$$     &                    +9.33d-1*log10(teff/4.d4)
c$$$     &                    -10.92d0*(log10(teff/4.d4))**2
c$$$     &                    +0.85d0*log10(zeff/zsol)
c$$$                  endif
c$$$                  
c$$$                  dms = 10.d0**logdms
c Mass loss rate from Vink et al 2000 modified to account for the
c effects of rotation (factor 0.85)
c                  dms = 0.8d0*dms
cc Mass loss rate modified to account for the effect of clumping!
cc (reduction by a factor of 3 ; see Smith 2014 ARA&A)
c     dms = dms*pw13
                  print*,'clumpfac in mchange',clumpfac
                  dms = dms*clumpfac
               endif
            endif

***   Crowther (2000) : for red-super giants
            if (mlp.eq.19.or.mlp.eq.20) then
               logdms = -1.7d0*log10(soll)+13.83d0
               dms = 10.d0**(-logdms)
            endif

***   van loon et al (2005, A&A, 438, 273) : for AGB & red supergiants
            if (mlp.eq.21.or.mlp.eq.22) then
               logdms = -5.65d0+1.05d0*log10(soll*1.d-4)-6.3d0*
     &              log10(teff/3500.d0)
               dms = 10.d0**logdms
            endif

            if (mlp1.eq.15.or.mlp1.eq.16) mlp = mlp1

            if (mlp.eq.17.or.mlp.eq.18) then
               dms = massrate
               dmlinc = 1.d0
            endif

***   Gräfener (2021, A&A 647, A13) : for Very Massive Stars (VMS)
            if (mlp.eq.25.or.mlp.eq.26) then
             if (mtini.le.100.d0) then
                 if (mlp.eq.25) mlp = 15
                 if (mlp.eq.26) mlp = 16
              else
               tau1 = pw23
               tau2 = 1.d0
               gammae = 10.d0**(-4.813)*(1.d0+xsp(nmod,ih1))*lum(nmod)
     $              /(lsun*totm)
               if (hydrorot) then
                  gammarot = vomega(neff)**2*Rp2(neff)**3/(g*m(neff))
                  gammaeff = gammae + 0.5d0*gammarot
               else
                  gammaeff = gammae
               endif
               Dclump = 10.d0
               vesceff = dsqrt(2.d0*g*m(nmod)*(1.d0-gammae)/reff)
               vinf = 2.51d0 * vesceff
               dmsthick = 5.22d0*log10(gammaeff)-0.5d0*log10(Dclump)
     $              -2.6d0
               dmsthick = 10.d0**dmsthick
               etawind = dmsthick*vinf/(lum(nmod)/c)
               taus = dmsthick*msun*seci*vinf*c/lum(nmod)*(1.d0+1.d0
     $              /(2.51d0**2))
***   switch to WNh type mass loss for VMS on main sequence according to
***   their sonic point optical depth - Graefener2021 and Bestenlehner+2014
c               ratio = 2.6d0
               ratio = 2.51d0
               logdms1 = -6.697d0+2.194d0*log10(soll*1.d-5)
     &              -1.313d0*log10(totm/30.d0)
     &              -1.226d0*log10(ratio*0.5d0)
     &              +9.33d-1*log10(teff/4.d4)
     &              -10.92d0*(log10(teff/4.d4))**2
     &              +0.85d0*log10(zeff/0.019d0)
c               dms1 = 0.8d0*10.d0**logdms1
               dms1 = pw13*10.d0**logdms1
               dms2 = dmsthick
               if (taus.lt.tau1) then
                  write(90,*)'Modified Vink+01 mass loss for OB VMS ' ,
     &                 'according to Graefener21'
                  dms = dms1
               else if (taus.gt.tau2) then
                write(90,*) 'WNh type mass loss for optically thick',
     $                 ' winds according to Graefener21'
                   dms = dms2
                else
                 write(90,*) 'Intermediate regime: interpolation btw',
     $                 ' thin and thick winds according to Graefener21'
                  
                  alfa = (taus-pw23)*3.d0
                  dms =  (1-alfa)*dms1 + alfa*dms2 
               endif
            endif
         endif
***   Mass loss recipe combining Vink for OB stars and Sander et Vink 2020 for VMS in a similar way as done by Grafener
***   Mass loss for VMS from Sander & Vink 2020, MNRAS 499, vol 1 p 873-892, Eqs (28) to (32)
            if (mlp.eq.27.or.mlp.eq.28) then
               if (mtini.le.100.d0) then
                  if (mlp.eq.27) mlp = 15
                  if (mlp.eq.28) mlp = 16
               else
                  tau1 = pw23
                  tau2 = 1.d0
                  gammae = 10.d0**(-4.813)*(1.d0+xsp(nmod,ih1))
     $                 *lum(nmod)/(lsun*totm)
                  if (hydrorot) then
                     gammarot = vomega(neff)**2*Rp2(neff)**3/(g*m(neff))
                     gammaeff = gammae + 0.5d0*gammarot
                  else
                     gammaeff = gammae
                  endif
                  Dclump = 10.d0
                  gammaeb = -0.324d0*log10(zeff/zsol)+0.244
                  logmdotoff = 0.23d0**log10(zeff/zsol)-2.61
                  cbd = -0.44d0*log10(zeff/zsol)+9.15
                  acoeff = 2.932d0
                  logdmsvms = acoeff*log10(-log10(1-gammae))-log10(2.d0)
     $                 *(gammaeb/gammae)**cbd+logmdotoff-0.5d0
     $                 *log10(Dclump)
                  vesceff = dsqrt(2.d0*g*m(nmod)*(1.d0-gammae)/reff)
                  vinf = 2.51d0 * vesceff
                  taus = 10.d0**(logdmsvms)*msun*seci*vinf*c/lum(nmod)
     $                 *(1.d0+1.d0/(2.51d0**2))
                  print *,'taus',taus, logdmsvms,vinf*c/lum(nmod)
     $                 ,gammae, gammaeb,logmdotoff
***   switch to WNh type mass loss for VMS on main sequence according to
***   their sonic point optical depth - Graefener2021 and Bestenlehner+2014
c     ratio = 2.6d0
                  ratio = 2.51d0
                  logdms1 = -6.697d0+2.194d0*log10(soll*1.d-5)
     &                 -1.313d0*log10(totm/30.d0)
     &                 -1.226d0*log10(ratio*0.5d0)
     &                 +9.33d-1*log10(teff/4.d4)
     &                 -10.92d0*(log10(teff/4.d4))**2
     &                 +0.85d0*log10(zeff/0.019d0)
c     dms1 = 0.8d0*10.d0**logdms1
                  dms1 = pw13*10.d0**logdms1
                  dms2 = 10**logdmsvms
                  if (taus.lt.tau1) then
                     write(90,*)'Modified Vink+01 mass loss for OB ',
     &                    'VMS according to Graefener21'
                     dms = dms1
                  else if (taus.gt.tau2) then
                     write(90,*) 'He VMS mass loss according to Sander & 
     &Vink 2020'
                     dms = dms2
                  else
                     write(90,*) 'Intermediate regime: interpolation',
     $' btw thin and He rich WR winds as in Graefener 2021'
                     
                     alfa = (taus-pw23)*3.d0
                     dms =  (1-alfa)*dms1 + alfa*dms2 
                  endif
            endif
            
         endif
         print *,'dms Sander and Vink 2020',dms,dms1,dms2,taus,tau1,tau2
        
            
***   Correction for rotating stars : Maeder & Meynet, 2001, A&A 373, 555
***   New version : Georgy et al. 2011 A&A 527, A52, Eq. (4)
            if (irotbin.eq.1) then
               malpha = 1.d0
               dms0 = dms
c...  mgamma = Eddington factor for electron scattering opacities
c            msigmae = 0.401d0*(zeff*0.25d0+xsp(neff,ih1)+xsp(neff,ih2)+
c     &           (xsp(neff,ihe3)+xsp(neff,ihe4))*0.5d0)
               msigmae = 0.401d0*(zeff*0.5d0+xsp(neff,ih1)+xsp(neff,ih2)
     &              +(xsp(neff,ihe3)+xsp(neff,ihe4))*0.5d0)
               mgamma = msigmae*lum(nmod)/(4*pi*g*c*m(nmod))
               if (hydrorot) then
                  ledd = pim4*c*g*m(neff)/kap(neff)*(1.d0-pw23
     &                 *vomega(neff)**2*Rp2(neff)**3/(g*m(neff)))
c               else
c                  ledd = pim4*c*g*m(neff)/kap(neff)
                  mgammaedd = lum(neff)/ledd
               else
                  mgammaedd = mgamma
               endif
               mrom = m(nmod)*pw34/(pi*r(nmod)**3)
               if (tteff.ge.4.35d0) malpha = 0.52d0
               if (tteff.ge.4.30d0.and.tteff.lt.4.35d0)
     &              malpha = 0.24d0
               if (tteff.ge.4.d0.and.tteff.lt.4.3d0)
     &              malpha = 0.17d0
               if (tteff.ge.3.9d0.and.tteff.lt.4.d0)
     &              malpha = 0.15d0
               mexp = 1.d0/malpha-1.d0
               vsurf = r(nmod)*vomega(nmod)
               mfac1 = vomega(nmod)**2*2.384989d6/mrom
c               omegacrit = dsqrt(g*m(nmod)/(1.5d0*r(nmod))**3)
c               mfac1 = 4.d0/27.d0*(vomega(nmod)/omegacrit)**2
c..   absolute value added on Jan 2013, 28th to solve NaN problem when
c..   (1.d0-mgamma-mfac1) < 0 in a 40 Msun rotating model at the end of the MS
C               mfac = dabs((1.d0-mgamma)/(1.d0-mgamma-mfac1))**mexp
               mfac = dabs((1.d0-mgamma)/((1.d0-mgammaedd)*
     &              (1.d0-mfac1)))**mexp
               print *,"mfac",mfac
               dms = dms*mfac
               if (abs(mfac-1.d0).gt.1.d-3) then
                  write (nout,600) mfac,dms0,dms
               endif
c..   In order to do as Ekstroem et al. 2012 (section 2.6.2), multiply the RSG
c..   mass loss by a factor of 3 whenever mgamma > 5 and M > 15 Msun
               if ((mlp.eq.19.or.mlp.eq.20).and.mtini.gt.15.d0.and.
     &              mgamma.gt.5.d0) then
                  dms = 3.d0*dms
               endif
            endif
         endif

         
*_____________________________________________________
***   change shells mass distribution due to mass-loss
*-----------------------------------------------------

***   Metallicity dependence
         zdep = ' '
         if (zscaling) then
C     Mokiem et al 2007
C            dms = dms*dsqrt(zeff/zsol)
            dms = dms*(zeff/zsol)**0.8d0
            zdep = ' with metallicity dependence activated'
         endif

         if (model.eq.modeli.and.iter.eq.1.and.ifail.eq.0) then
            if (mlp.eq.1.or.mlp.eq.2) prescrip = 'Reimers'
            if (mlp.eq.3.or.mlp.eq.4) prescrip = 'de Jager'
            if (mlp.eq.5.or.mlp.eq.6) prescrip = 'Vassiliadis & '//
     &           'Wood without delaying the onset of super-wind'
            if (mlp.eq.55.or.mlp.eq.56) prescrip = 'Vassiliadis & ' //
     &           'Wood original prescription'
            if (mlp.eq.7.or.mlp.eq.9) prescrip = 'Blocker'
            if (mlp.eq.9.or.mlp.eq.10) prescrip = 'Arndt'
            if (mlp.eq.11.or.mlp.eq.12) prescrip = 'Schaller et al.'
            if (mlp.eq.13.or.mlp.eq.14) prescrip = 'Chiosi'
            if (mlp.eq.15.or.mlp.eq.16) prescrip = 'Vink et al.'
            if (mlp.eq.17.or.mlp.eq.18) prescrip = 'massrate, '//
     &           'specified in starevol.par'
            if (mlp.eq.19.or.mlp.eq.20) prescrip = 'Crowther'
            prescrip = trim(prescrip) // zdep
            write (90,700) prescrip
         endif

 3       dms = dms*dmlinc
         dma = dms*dtn*seci
         totmv = totm*msun
         totm = totm-dma
         totmd = totm*msun
         if (dma*msun/(m(nmod)-m(ienvsh)).gt.1.d-3)
     &        mlsh = max(mlsh,min(ienvsh+100,nmod-100))


c...  Introduction of mechanical mass loss for massive stars at break-up
c...  (see Georgy et al. 2013, A&A 553, A24, section 2.2)
c...  Maximum authorized surface angular velocity = 99% of critical velocity

c$$$      omegamax = 0.99d0*dsqrt(g*m(nmod)/(1.5d0*r(nmod))**3)

c$$$c..   Inititial total angular momentum prior any mass loss is apply
c$$$      xmm = 0.d0
c$$$      do j = 2,nmod1
c$$$         j1 = j-1
c$$$         xmm = xmm+(r(j)**2+r(j)*r(j1)+r(j1)**2)/6.d0*
c$$$     &        (vomega(j)+vomega(j1))*dm(j1)
c$$$      enddo
c$$$      print *,""
c$$$      print *,'momentum before mass loss',xmm,omegamax,vomega(nmod)
c$$$
c$$$c..   Angular momentum lost by standard setllar winds
c$$$
c$$$      dLrad = pw23*dma*msun*vomega(nmod)*r(nmod)**2
c$$$
c$$$c..   Initial surface angular velocity
c$$$
c$$$      omsurfini = vomega(nmod)
c$$$
c$$$c..   Initial mass and angular momentum in region where mass loss is
c$$$c..   applied in case of uneven value of mlp
c$$$
c$$$      meini = m(nmod)-m(novlim(nsconv,3))
c$$$      Leini = 0.0d0
c$$$      do j = novlim(nsconv,3)+1,nmod
c$$$         j1 = j-1
c$$$         Leini = Leini+(r(j)**2+r(j)*r(j1)+r(j1)**2)/6.d0*
c$$$     &        (vomega(j)+vomega(j1))*dm(j1)
c$$$      enddo
c$$$


***   if mlp even : shells are removed
         if (mod(mlp,2).eq.0.and.dms.gt.0.d0) then
            nold = nmod
            do i = nmod,1,-1
               if (m(i).lt.totmd) goto 20
            enddo
 20         totmd = m(i)
            totm = totmd/msun
            dma = (totmv-totmd)/msun
            dms = dma*sec/dtn
            write (nout,800) nmod,i,totmv/msun,totm
            nmod = i
            nmod1 = nmod-1
            do k = 2,nmod1
               mr(k) = m(k)/totmd
            enddo
            mr(1) = 0.d0
            mr(nmod) = 1.d0
            dlogm = 0.d0
            do j = 1,nsp
               xtot = 0.d0
               do k = nmod,nold-1
                  xtot = xtot+xsp(k,j)*dm(k)
               enddo
               mtotlos(j) = mtotlos(j)+xtot/msun
            enddo
!***   if mlp odd : mass removed in shells above m = m(mlsh)

         else
            if (nmod0.ne.0) then
               if (dma*msun.gt.(m(nmod)-m(nmod0))) then
                  nmod = nmod0
                  nmod1 = nmod-1
               endif
            endif

            do j = 1,nsp
               do k = 1,nmod
                  pxsp(k,j) = xsp(k,j)
               enddo
            enddo

c..   mass uniformally removed in the (1-fxm) percent of the star
            fxm = 0.95d0
            xxm = m(mlsh)
            if (totmd.lt.xxm) then
               xxm = fxm*totmd
               do k = mlsh,1,-1
                  if (m(k).lt.xxm) goto 30
               enddo
 30            mlsh = k
               xxm = m(mlsh)
            endif

            dxxm = (totmd-xxm)/(totmv-xxm)
            do k = mlsh,nmod1
               dm(k) = vdm(k)*dxxm
               m(k+1) = m(k)+dm(k)
            enddo


c.. For massive stars, mass is only removed below the "envelope"
c
c            if (totm.gt.15.d0) then
c              k = nmod
c               menv = 0.985d0*totmv
c               do while (m(k).ge.menv)
c                  k = k-1
c               enddo
c               nenv = k
c
c               dxxm = (totmd-totmv)/(menv-xxm)+1.d0
c               do k = mlsh,nenv-1
c                  dm(k) = vdm(k)*dxxm
c               enddo
c               do k = nenv,nmod1
c                  dm(k) = vdm(k)
c               enddo
c               do k = mlsh,nmod1
c                  m(k+1) = m(k)+dm(k)
c               enddo
c            else
c... version Wagenhuber & Weiss, 1994, A&A, 286, 121
c..  not accurate for low mass los rates
c                 dxxm = dma*msun/m(nmod)
c                 xxm = 0.05d0/(dma*msun)
c                 m(1) = 0.d0
c                 mlsh = nmod
c                 do k = 2,nmod
c                    dmsreim = 0.d0
c                    do i = k,nmod
c                       dmsreim = dmsreim-dm(i)
c                    enddo
c                    dmsreim = dmsreim*xxm
c                    if (dmsreim.gt.-26.d0) then
c                       print *,k,exp(dmsreim),dxxm*dexp(dmsreim)
c                       if (mlsh.eq.nmod) mlsh = k
c                       m(k) = m(k)*(1.d0-dxxm*dexp(dmsreim))
c                       dm(k-1) = m(k)-m(k-1)
c                    endif
c                 enddo
c                print *,m(nmod),totmd,totmv,(totmv-totmd)/msun,(m(nmod)
c     &               -totmv)/msun,dma,(totmd-m(nmod))/m(nmod)
c            endif

            dm(nmod) = 0.d0
            do k = 2,nmod1
               mr(k) = m(k)/totmd
            enddo
            mr(1) = 0.d0
            mr(nmod) = 1.d0
            dlogm = totmd/totmv-1.d0
c            dlogm = -dma*msun/totmv
c            if (totm.gt.10.d0) nend = nenv
c            if (totm.le.10.d0) nend = nmod

c..   Interpolate new chemical composition
            do i = mlsh+1,nmod1
               mm = mlsh
 40            mm = mm+1
               if (vm(mm).gt.m(i)) then
                  ff = (m(i)-vm(mm-1))/(vm(mm)-vm(mm-1))
                  do j = 1,nsp
                     xsp(i,j) = pxsp(mm-1,j)+ff*(pxsp(mm,j)-
     &                    pxsp(mm-1,j))
                  enddo
               else
                  goto 40
               endif
            enddo
            do j = 1,nsp
               xtot = 0.d0
               do k = mlsh+1,nmod1
                  xtot = xtot+xsp(k,j)*(vdm(k)-dm(k))
               enddo
               mtotlos(j) = mtotlos(j)+xtot
            enddo
         endif
      endif

*___________________________
***   treatment of accretion
*---------------------------

      if (dmaccr.gt.0.d0) then

         iaccbot = 0
         iacctop = 0
         if (accphase.gt.0.and.accphase.lt.5) mixopt = .true.

*________________________________________________________
***   change in shells mass distribution due to accretion
*--------------------------------------------------------

         m(1) = 0.d0
         macc(1) = 0.d0
         vtotma = m(nmod)

         do k = 1,nmod1
            facc(k) = facc(k)*sinthac
            if (iaccr.eq.2) then
               macc(k+1) = macc(k)+facc(k)*dm(k)
               dm(k) = dm(k)*(1.d0+facc(k))
            else
               macc(k+1) = macc(k)+facc(k)/(1.d0+facc(k))*dm(k)
               dm(k) = dm(k)*(1.d0+facc(k)/(1.d0+facc(k)))
            endif
            vdm(k) = dm(k)
            do j = 1,nsp
               pxsp(k,j) = xsp(k,j)
            enddo
            m(k+1) = m(k)+dm(k)
            if (facc(k).gt.0.d0.and.iaccbot.eq.0) iaccbot = k
            if (facc(nmod+1-k).gt.0.d0.and.iacctop.eq.0) iacctop =
     &           nmod+1-k
         enddo
         do j = 1,nsp
            pxsp(nmod,j) = xsp(nmod,j)
         enddo
         facc(nmod) = facc(nmod)*sinthac
         dmaccrg = macc(nmod)
         iaccbot = max(iaccbot,2)

         dm(nmod) = 0.d0
         vdm(nmod) = 0.d0
         dmacc1 = m(nmod)-vtotma
         dmacc2 = dmaccrg
         ddmacc = dmacc1-dmacc2
         totmv = totm*msun
         totmd = m(nmod)
         totm = totmd/msun
         dmaccr = dmaccrg/msun
         time = time-dtn
         dtn = dmaccr*sec/massrate
         time = time+dtn

c         dlogm = dlogm+log(totmd/totmv)
         dlogm = dlogm+totmd/totmv-1.d0
         do k = 2,nmod1
            yd1(k) = xsp(k,ih2)
            yd0(k) = xsp(k,ih2)
            vyd0(k) = xsp(k,ih2)
            mr(k) = m(k)/totmd
         enddo
         vyd0(1) = xsp(1,ih2)
         yd0(1) = xsp(1,ih2)
         yd1(1) = xsp(1,ih2)
         mr(1) = 0.d0
         mr(nmod) = 1.d0
         vyd0(nmod) = xsp(nmod,ih2)
         yd0(nmod) = xsp(nmod,ih2)
         yd1(nmod) = xsp(nmod,ih2)

*_________________________________________________________________
***   change in chemical composition
*-----------------------------------------------------------------
* itacc = 0 : matter pills-up at the surface of the star
*             surface composition = composition of accreted matter
* itacc > 0 : mix accreted matter with stellar matter
*-----------------------------------------------------------------

c..  mix accreted matter with stellar matter
         if (accphase.le.4.or.itacc.gt.0) then
            if (facc(1).gt.0.d0) then
               if (iaccr.eq.1.or.iaccr.eq.3) then
                  do j = 1,nsp
                     xsp(1,j) = pxsp(1,j)*(1.d0-facc(1))+facc(1)*
     &                    xspacc(j)
                  enddo
               endif
               if (iaccr.eq.2.or.iaccr.eq.4) then
                  do j = 1,nsp
                     xsp(1,j) = (pxsp(1,j)+facc(1)*xspacc(j))/(1.d0+
     &                    facc(1))
                  enddo
               endif
            endif
            if (facc(nmod).gt.0.d0) then
               do j = 1,nsp
                  xsp(nmod,j) = xsp(nmod1,j)
               enddo
            endif
            if (iaccr.eq.1.or.iaccr.eq.3) then
               do j = 1,nsp
                  do i = iaccbot,min(iacctop,nmod-1)
                     xsp(i,j) = pxsp(i,j)*(1.d0-facc(i))+facc(i)*
     &                    xspacc(j)
                  enddo
               enddo
            endif
            if (iaccr.eq.2.or.iaccr.eq.4) then
               do j = 1,nsp
                  do i = iaccbot,min(iacctop,nmod-1)
                     xsp(i,j) = (pxsp(i,j)+facc(i)*xspacc(j))/(1.d0+
     &                    facc(i))
                  enddo
               enddo
            endif
         else
c..   Interpolate new chemical composition, pill-up mass on top of star
c..   surface composition = composition of accreted matter
            do i = nmod,iaccbot,-1
               mm = nmod
               if (m(i).gt.vm(mm)) then
                  do j = 1,nsp
                     xsp(i,j) = xspacc(j)
                  enddo
               else
 50               mm = mm-1
c..  if accretion inside the star (planet), interpolate composition
                  if (m(i).ge.vm(mm)) then
                     ff = (m(i)-vm(mm))/(vm(mm+1)-vm(mm))
                     do j = 1,nsp
                        xsp(i,j) = pxsp(mm,j)+ff*(pxsp(mm+1,j)-
     &                       pxsp(mm,j))
                     enddo
                  else
                     goto 50
                  endif
               endif
            enddo
         endif
      endif

      if (numeric.eq.2.or.numeric.eq.3) then
c.. first order accuracy in the spatial derivatives
         forall (i = 1:nmod)
            wi(i) = 0.5d0
            wj(i) = 0.5d0
         end forall
      else
c.. second order accuracy in the spatial derivatives
         forall (i = 2:nmod1)
            wi(i) = dm(i)/(dm(i)+dm(i-1))
            wj(i) = 1.d0-wi(i)
         end forall
         wi(1) = 0.5d0
         wj(1) = 0.5d0
         wi(nmod) = 0.5d0
         wj(nmod) = 0.5d0
      endif

***   if shells are removed no change of independent variable
      if (mod(mlp,2).eq.0.and.iaccr.eq.0) then
         if (irotbin.eq.1) then
            xmom_tots = 0.d0
            do j = 2,nmod1
               j1 = j-1
               xmom_tots = xmom_tots+(r(j)**2+r(j)*r(j1)+r(j1)**2)
     &              /6.d0*(vomega(j)+vomega(j1))*dm(j1)
            enddo
         endif

         return
      endif


*___________________________________________________________________
***               change of independent variable
***   define new independent variable qmr in the accretion/mass-loss
***   region and interpolate old variables at the new mesh point qmr
*-------------------------------------------------------------------


***   define new mesh point
      iqmrbot = min(mlsh,iaccbot)
      if (dup3) iqmrbot=max(iqmrbot,ienvsh+10)
      iqmrtop = nmod

      qmr(1) = 0.d0
      vqmr(1) = 0.d0
      if (iqmrbot.gt.1) then
         forall (i = 1:iqmrbot)
            qmr(i) = 0.d0
            vqmr(i) = 0.d0
         end forall
      endif
***   by construction vm(iqmrbot) = m(iqmrbot)
      ddm = 0.d0
      vddm = 0.d0
      do k = iqmrbot+1,nmod
         ddm = ddm+dm(k-1)
         vddm = vddm+vdm(k-1)
      enddo
      ddm = 1.d0/ddm
      vddm = 1.d0/vddm
      do k = iqmrbot+1,nmod
         qmr(k) = qmr(k-1)+dm(k-1)*ddm
         vqmr(k) = vqmr(k-1)+vdm(k-1)*vddm
      enddo
      qmr(iqmrbot-1) = -dm(iqmrbot-1)*ddm
      vqmr(iqmrbot-1) = -vdm(iqmrbot-1)*vddm
      qmr(iqmrbot) = 1.d-20

***   define omi & psi in the accretion/mass-loss region
      do k = iqmrbot+1,nmod1
         if (qmr(k).eq.qmr(k-1)) write (nout,900) k,dm(k)
         omi(k) = dlogm/log((qmr(k)+qmr(k+1))/(qmr(k-1)+qmr(k)))
         psi(k) = dlogm/log(qmr(k)/qmr(k-1))
      enddo
      omi(nmod) = omi(nmod-1)
      psi(nmod) = dlogm/log(qmr(nmod)/qmr(nmod1))

ccccccccccccccccccccccccccccccc
c      do k = iqmrbot+1,nmod
c         omi(k) = 0.d0
c         psi(k) = 0.d0
c      enddo
c      return
ccccccccccccccccccccccccccccccc

      if (irotbin.eq.1) forall (i=1:nmod) vvr(i) = r(i)

      ij = iqmrbot
      do j = iqmrbot,nmod-1
         itop = 0
         do k = ij,nmod
            if (vqmr(k).gt.qmr(j)) then
               itop = k
               goto 60
            endif
            ibot = k
         enddo
 60     if (itop-ibot.ne.1) then
            write (nout,*) ibot,itop,nmod,iqmrbot,nmod,dm(max(itop,1))
            stop 'mchange : shell problem, mass too small ?'
         endif
         if (itop.eq.0) stop 'mchange : variables interpolation failed'
         ij = itop-1
         ii = j
         call interpmesh (ii,itop,4)
      enddo


*________________________________________________________________
***   In order to ensure angular momentum conservation vomega has
***   to be rescaled according to the changes made on r.
***   Rescaling of vomega by propagation.
*----------------------------------------------------------------

      if (irotbin.eq.1) then
         xmom_tots = 0.d0
c..    Omega cst in the convective envelope
         if (idiffvr.le.5.or.idiffvr.eq.7.or.idiffvr.ge.8) then
            mominnew = 0.d0
            mominold = 0.d0
            do j = novlim(nsconv,3)+1,nmod
               j1 = j-1
               mominold = mominold+(vvr(j)**2+vvr(j)*vvr(j1)+
     &              vvr(j1)**2)*dm(j1)
               mominnew = mominnew+(r(j)**2+r(j1)**2+r(j)*r(j1))*
     &              dm(j1)
            enddo
            do j = novlim(nsconv,3)+1,nmod
               j1 = j-1
               vomega(j) = (vomega(j)+vomega(j1))*mominold/mominnew-
     &              vomega(j1)
            enddo
c..   j = cst in the convective envelope
         else if (idiffvr.eq.6) then
            do j = novlim(nsconv,3)+1,nmod
               j1 = j-1
               vomega(j) = (vomega(j)+vomega(j1))*(vvr(j)**2+vvr(j)*
     &              vvr(j1)+vvr(j1)**2)/(r(j)**2+r(j1)**2+r(j)*r(j1))-
     &              vomega(j1)
            enddo
         endif

         momlim = 0.d0
         do j = 2,nmod1
            j1 = j-1
            if (j.le.novlim(nsconv,3)) then
               momlim = momlim + (r(j)**2+r(j)*r(j1)+r(j1)**2)/6.d0*
     &              (vomega(j)+vomega(j1))*dm(j1)
            else
               momlim = momlim + (r(j)**2+r(j)*r(j1)+r(j1)**2)*pw13*
     &              omegamax*dm(j1)
            endif
         enddo


c..   Total AM after mass loss
         do j = 2,nmod1
            j1 = j-1
            xmom_tots = xmom_tots + pw13*(r(j)**2+r(j1)**2+r(j)*r(j1))
     &           *dm(j1)*0.5d0*(vomega(j)+vomega(j1))
         enddo

c$$$         rapom = omegamax/omsurfini
c$$$
c$$$         if (xmom_tots.gt.momlim) then
c$$$            dmsmec = (xmm*(1.d0-rapom)+Leini/meini*rapom*dma*msun-dLrad)
c$$$     &           /(omsurfini*1.5d0*r(nmod)**2-Leini/meini*rapom)/msun
c$$$c            dmsmec =  (-momlim
c$$$c     $        +xmom_tots)/(omsurfini*(1.5d0*r(nmod))**2)/msun
c$$$c            dms = dmsmec
c$$$c            goto 3
c$$$         endif
      endif

      print *,'xmom_tots in mchange',xmom_tots,momlim,xmom_tots/momlim

 100  format (5x,'WARNING : mchange mass shell [',i4,'] = ',1pe16.9,
     &     ' TOO SMALL !!!!')
 200  format (5x,'WARNING : Change mass loss --> de Jagger')
 300  format (5x,'WARNING : Change mass loss --> Crowther (2000)')
 400  format (5x,'WARNING : Change of Mass Loss Regim: mlp = ',i2,' !')
 500  format (5x,'Entering WR phase - Mass loss from Nugis & Lamers')
 600  format (5x,'Mass loss corrected by a factor of ',1pe10.4,' due ',
     &     'to rotation : Mloss = ',1pe10.4,' --> ',1pe10.4,' Mo/yr')
 700  format (/,'o MASS LOSS PRESCRIPTION',/,2x,A,/)
 800  format (5x,'shells removed : ',i4,'-->',i4,', M =',0pf8.5,'-->',
     &     0pf8.5)
 900  format (5x,'WARNING : qmr shell [',i4,'] = ',1pe16.9,
     &     ' TOO SMALL !!!!')

      return
      end
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
      incrATM = .false.      ! Modif 17/04/2020
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
c     lnucl = voir carte des paramètres
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
         
c..   during core He burning stabilizes mesh by increasing dnucmax
         if (nphase.eq.4) then
            print *, 'phase 4 reached -> increase dnucmax'
            dnucmax = 1.d50
         else
            dnucmax = dlnenuc
         endif
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


***********************************************************************

      SUBROUTINE mix (dm,imin,imax,kl,imix,partialmix)

************************************************************************
* Instantaneous mixing of species in convective zones                  *
* imix = 0 : default, do the mixing
* imix = 1 : special treatment of deuterium (accretion)
* imix = 2 : before nucleosynthesis, do the mixing
* imix = 3 : after nucleosynthesis, do the mixing
* imix = 5 : compute turnover timescales only
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.nuc'
      include 'evolcom.spec'
      include 'evolcom.teq'

c..   input
      double precision dm(nsh)

      integer imin,imax,imix
      integer ielabond
      integer i,ii,j,k,l,kl
      integer klenv,klpulse,klcore

      logical partialmix,irecede

      double precision mconv,drconv,tmix,dtni,tconvmix,dmeff
      double precision x,xk,xj
      double precision x0,x1
      double precision t,vt
      double precision r,vr

      dimension tmix(nsh)
      dimension ielabond(7)

      common /vartvt/ t(nsh),vt(nsh)
      common /varrvr/ r(nsh),vr(nsh)
      common /overshoot/ klcore,klenv,klpulse

      if (imin.ge.imax.or.dtn.lt.1.d-12) return

      tconvmix = 0.d0
      partialmix = .false.

*____________________________________________
***   compute convective turnover time scales
*--------------------------------------------

      if (imix.ne.3) then
         do j = imin,imax-1
            tconvmix = tconvmix+(r(j+1)-r(j))/sconv(j)
         enddo
         turnover(kl) = tconvmix
         if (nsconv.eq.1) then
            turnenv = turnover(kl)
         else
            if (t(imin).gt.1.d6.and.r(imax)/r(nmod).gt.0.9d0) then
               turnenv = turnover(kl)
            endif
         endif
      else
         tconvmix = turnover(kl)
      endif
      if (nuclopt.eq.'p'.and.tconvmix.gt.dtn) partialmix = .true.

      if (imix.eq.5) return

c..   do not mix before nucleosynthesis (imix=2) if dtnmix = t
c..   or tconvmix >> dtn (partialmix = .true.)
      if (imix.eq.2.and.(dtnmix.eq.'t'.or.partialmix)) return


*__________________________________
***   Instantaneous mixing
***   if dtnmix = f or tconv >> dtn
*----------------------------------

      if (dtnmix.eq.'f'.or..not.partialmix) then
         mconv = 0.d0
         do j = imin,imax
            mconv = mconv+dm(j)
         enddo
!$OMP PARALLEL REDUCTION(+:x)
!$OMP DO
         do l = 2,nis
            x = 0.d0
            do j = imin,imax
               x = x+vxsp(j,l)*dm(j)
            enddo
            xj = x/mconv
            do j = imin,imax
               xsp(j,l) = xj
               ysp(j,l) = xj/anuc(l)
            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL

*___________________________________________________
***   special treatment of LiBeB during HBB
***   restore old (equilibrium) abundances for LiBeB
*---------------------------------------------------

         if (imix.eq.2.or.imix.eq.3) write (nout,200) imin,imax
         if (hbb.and.kl.eq.klenv.and.partialmix) then
            do k = imin,imax
               xsp(k,ili7) = vxsp(k,ili7)
               xsp(k,ili6) = vxsp(k,ili6)
               xsp(k,ib11) = vxsp(k,ib11)
               ysp(k,ili7) = xsp(k,ili7)/anuc(ili7)
               ysp(k,ili6) = xsp(k,ili6)/anuc(ili6)
               ysp(k,ib11) = xsp(k,ib11)/anuc(ib11)
            enddo
            if (imix.eq.2.or.imix.eq.3) write (nout,250) imin,imax
         endif

      else

*_______________________________________________
***   time-dependent mixing when tauconv < dtn
***   if dtnmix = t and partialmix = .true
***
*** Sparks Endal, 1980, ApJ 237, 130
*** Chieffi et al. 1998, ApJ 502,737
*-----------------------------------------------

c..   check that the convective zone has moved into regions
c..   of different chemical composition
         if (iter.eq.1) then
            irecede = .false.
         elseif (mixopt) then
            irecede = .true.
            ielabond(1) = ih1
            ielabond(2) = ili7
            ielabond(3) = ic12
            ielabond(4) = io16
            ielabond(5) = ine20
            ielabond(6) = img24
            ielabond(7) = isi28
            do i = 1,7
               ii = ielabond(i)
               if (log10(abs(xsp(imin,ii)/vxsp(imin,ii))).gt.1.d-2.or.
     &              log10(abs(xsp(imax,ii)/vxsp(imax,ii))).gt.1.d-2) 
     &              then
                  irecede = .false.
                  goto 10
               endif
            enddo
         else
            irecede = .false.
         endif
 10      if (irecede.and.imix.ne.3) return

         write (nout,100) imin,imax,tconvmix*seci,dtn*seci
         dtni = 1.d0/dtn
         do k = imin,imax
            tmix(k) = 0.d0
            do j = k-1,imin,-1
               drconv = r(j+1)-r(j)
               tmix(j) = tmix(j+1)+drconv/sconv(j)
            enddo
            do j = k+1,imax
               drconv = r(j)-r(j-1)
               tmix(j) = tmix(j-1)+drconv/sconv(j-1)
            enddo
            do i = 2,nis
               xk = 0.d0
               mconv = 0.d0
               do l = imin,imax
                  dmeff = dm(l)*exp(-tmix(l)*dtni)
                  xk = xk+(vxsp(l,i)-vxsp(k,i))*dmeff
                  mconv = mconv+dmeff
               enddo
               xsp(k,i) = vxsp(k,i)+xk/mconv
            enddo
         enddo
         return
      endif

*_______________________________________________________________
***   specific mixing of light elements, in case of call ydequil
*---------------------------------------------------------------

      if (imix.eq.1) then
         imin = max(imin,idcur)
         x0 = 0.d0
         x1 = 0.d0
         mconv = 0.d0
         do j = imin,imax
            x0 = x0+vyd0(j)*dm(j)
            x1 = x1+yd0(j)*dm(j)
            mconv = mconv+dm(j)
         enddo
         x0 = x0/mconv
         x1 = x1/mconv
         do j = imin,imax
            vyd1(j) = x0
            yd1(j) = x1
         enddo
         do l = 2,4
            x = 0.d0
            do j = imin,imax
               x = x+vxsp(j,l)*dm(j)
            enddo
            x = x/mconv
            do j = imin,imax
               vxsp(j,l) = x
            enddo
         enddo
      endif

 100  format (' partial mixing in shells : [',i4,',',i4,'] , ',
     &     'tau_conv = ',1pe9.3,' > dtn = ',1pe9.3)
 200  format ('  instantaneous mixing in convective zone [',i4,',',
     &     i4,']')
 250  format ('   HBB : special mixing for Li6, Li7 and B11 in the ',
     &     'envelope [',i4,',',i4,']')

      return
      end


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
      double precision delroa,cua,adrad,cc,kiminv,kiinv
      double precision con1,con2,con3,alphac2,ut,ut2,ze,al,
     &     epsf,we,we2,dr1,dr2,dr3,dr4,dr5,dr6,dr7,dutdf1,dutdf2,
     &     dutdt1,dutdt2,dutrr,daldf1,daldf2,daldt1,daldt2,dalrr,daldl,
     &     dalu1,dalu2,wedf1,wedf2,wedt1,wedt2,werr,wedl,weu1,weu2,
     &     zedf1,zedf2,zedt1,zedt2,zerr,zedl,zeu1,zeu2,cua2,tol,z0,
     &     dzadrad,ziter,dznew,znew,ddz,zn,yn,xn4,xn,cu
      double precision xconv,sconv0,xsiconv1,xsiconv
      double precision a1conv,a2conv,a3conv,rconv,qconv,rqconv,
     &     rsign,rqexp,yconv,ztm,csm,abel,a1,A

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

               tconvb = lambdab(i)/sconvb
               tconvt = lambdat(i)/sconvt

               sconv(i) = 0.5d0*(sconvb+sconvt)
               eff(i) = 0.5d0*(efft+effb)
               tconv(i) = (tconvt+tconvb)*0.5d0
               taumix(i) = 2.d0*dtn*(sconv(i)+vsconv(i))/lambdac(i)
               Dconv(i) = (sconvb*lambdab(i)+sconvt*lambdat(i))/6.d0


               adraddf1 = abdf1(i)-wi(i)*abm(i)*dabadf(im)/abad(im)
               adraddf2 = abdf2(i)-wj(i)*abm(i)*dabadf(i)/abad(i)
               adraddt1 = abdt1(i)-wi(i)*abm(i)*dabadt(im)/abad(im)
               adraddt2 = abdt2(i)-wj(i)*abm(i)*dabadt(i)/abad(i)

                dkvc8 = yy*(yy+V_b)

               abdf1(i) = wi(i)*abm(i)*dabadf(im)/abad(im)+dkvc8*
     &              adraddf1 
               abdf2(i) = wj(i)*abm(i)*dabadf(i)/abad(i)+dkvc8*adraddf2
               abdt1(i) = wi(i)*abm(i)*dabadt(im)/abad(im)+dkvc8*
     &              adraddt1 
               abdt2(i) = wj(i)*abm(i)*dabadt(i)/abad(i)+dkvc8*adraddt2
               abdl(i) = dkvc8*abdl(i)
               abdr(i) = dkvc8*abdr(i) 
               abdu1(i) = dkvc8*abdu1(i) 
               abdu2(i) = dkvc8*abdu2(i) 
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
            ca = delroa/cua2
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
            fkin(i) = etaturb*rom(i)*sconv(i)**8/csm**5 ! Eq 1 Maeder 1987 A&A 173 : etaturb = 1e3 
            fkcr(i) = fkin(i)/fconv(i)
            frad(i) = lum(i)/(pim4*r(i)*r(i))-fconv(i)-fkin(i)
            Dconv(i) = sconv(i)*lambdac(i)/3.d0
c     tconv(i) = (r(imax)-r(imin))**2/Dconv(i)
            tconv(i) = lambdac(i)/sconv(i)
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
      if (abs(iross).le.imax) tconv_midCE = tconv(iross)
   
 1000 format (5x,'50 iterations reached for MLT with hro, at shell ',
     &     i4,', with ziter = ',1pe10.4,' and znew = ',1pe10.4)

      return
      end


************************************************************************

      SUBROUTINE netdiff (imin,imax,cd,dtn,error)

************************************************************************
* Compute the abundance evolution resulting from the coupling of the   *
* nucleosynthesis and diffusion equations                              *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eos'
      include 'evolcom.ion'
      include 'evolcom.therm2'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.spec'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer i,j,k,l,imin,imax,ielem,error
      integer iterd,neqd,neqd1,neqd2,mm,id,i1
      integer ncmmax,idmax,ntry,istep,itermax

      logical nequil,xequil

      double precision vv(nsh,nreac),vi(nreac)
      double precision f1(nsh),f2(nsh),f3(nsh),dd(nsh),cd(nsh),dmdr(nsh)
      double precision alphad,xx,xidk,xtol
      double precision eqd,xmod,xa,y,yeq
      double precision dx,dxmax,dv,w,a,um,v,alpha1,cormmax
      double precision xm_sum,dt
      double precision tk,rok,dtn,dtnuc,dtold,fac
      
      dimension eqd(nsp,3*nsp+1),xmod(nsh,nsp),xa(nsh,nsp),y(nsp)
      dimension dx(nsp),dv(nsp),w(nsh,nsp),a(nsp,2*nsp),um(nsp,nsp),
     &     v(nsh,nsp,nsp),alpha1(nsp)

      xequil = .false.
      nequil = .true.

      itermax = 1000
      icrash = 0
      ntry = 0
      istep = 0
      dt = 0.d0
      dtnuc = dtn*1.d-1

c..   computation nuclear reaction rates in all shells
      do k = imin,imax
         tk = 0.5d0*(vt(k)+t(k))
         rok = 0.5d0*(vro(k)+ro(k))
         call vit (tk,rok,mueinv(k),k,vi,0)
         do j = 1,nreac
            vv(k,j) = vi(j)
         enddo
      enddo

c..   initialisations
      alphad = 1.d0
      xtol = 1.d-8*tolnuc
      neqd = nsp
      neqd2 = 2*neqd
      neqd1 = 3*neqd+1

      do j = 1,nsp
         do i = imin,imax
            xmod(i,j) = xsp(i,j)/anuc(j)
            xa(i,j) = xmod(i,j)
         enddo
      enddo
      do i = 1,neqd
         alpha1(i) = alphad
      enddo

      do k = imin,imax
         if (k.eq.1) then
            dmdr(1) = 0.d0
            f1(1) = dtnuc/dm(1)
            f2(1) = 2.d0*f1(1)/(dm(1)+dm(2))
            f3(1) = 0.d0
         else
            dmdr(k) = pim4*r(k)*r(k)*rom(k)
            f1(k) = dtnuc/dm(k)
            f2(k) = 2.d0*f1(k)/(dm(k)+dm(k+1))
            f3(k) = 2.d0*f1(k)/(dm(k)+dm(k-1))
         endif
         dd(k) = cd(k)*dmdr(k)*dmdr(k)
      enddo
      dd(1) = dd(2)
      do k = imax-1,imin,-1
         dd(k+1) = dd(k)**wi(k+1)*dd(k+1)**wj(k+1)
      enddo
      dd(1) = dd(2)

*_______________________________________________________________________
***   calculation of the coupled nucleosynthesis and diffusion equations
***   Newton Raphson procedure
*-----------------------------------------------------------------------

      dtold = dtnuc
      iterd = 0
 10   iterd = iterd+1
      if (iterd.ge.itermax.or.dtnuc.lt.1.d-10) goto 20


c..   central boundary conditions
*--------------------------------

      id = imin
      do l = 1,neqd1
         do k = 1,neqd
            eqd(k,l) = 0.d0
         enddo
      enddo
      do j = 1,neqd
         y(j) = xmod(id,j)
      enddo
      do j = 1,nreac
         vi(j) = vv(id,j)
      enddo
c..   set neutrons to equilibrium abundance
      if (nequil) then
         call yneutron (y,vi,yeq)
         y(1) = yeq
      endif
c..   set  equilibrium abundance for elements with y < yeqmin
      if (xequil) call yequil (vi,y,ytminc,dtnuc,tconv(id))

      call freac (vi,eqd,y,dtnuc,.true.)
      if (kcentd.ne.3) call fdiff (eqd,dd,xmod,f1,f2,f3,id,imin,imax)

      do j = 1,neqd
         if (kcentd.ne.3) then
            eqd(j,neqd1) = eqd(j,neqd1)-xmod(id,j)+xa(id,j)
            eqd(j,j) = eqd(j,j)-1.d0
         else
            eqd(j,neqd1) = -xmod(id,j)+xmod(id-1,j)
            eqd(j,j) = -1.d0
         endif
         do l = 1,neqd
            a(l,j) = eqd(l,j)
            a(l,neqd+j) = 0.d0
         enddo
         a(j,neqd+j) = 1.d0
      enddo
      call nwrmat (a,um,neqd,neqd,error)
      if (error.gt.0) then
         error = 9
         return
      endif
      do mm = 1,neqd
         w(id,mm) = 0.d0
         do k = 1,neqd
            w(id,mm) = w(id,mm)+um(mm,k)*eqd(k,neqd1)
            v(id,mm,k) = 0.d0
            do j = 1,neqd
               v(id,mm,k) = v(id,mm,k)+um(mm,j)*eqd(j,neqd+k)
            enddo
         enddo
      enddo

      do l = 1,neqd1
         do k = 1,neqd
            eqd(k,l) = 0.d0
         enddo
      enddo

c..   internal shell diffusion equations
*---------------------------------------

      do id = imin+1,imax-1
         i1 = id-1

         do j = 1,neqd
            y(j) = xmod(id,j)
         enddo
         do j = 1,nreac
            vi(j) = vv(id,j)
         enddo
c..   set neutrons to equilibrium abundance
         if (nequil) then
            call yneutron (y,vi,yeq)
            y(1) = yeq
         endif
c..   set  equilibrium abundance for elements with y < yeqmin
         if (xequil) call yequil (vi,y,ytminc,dtnuc,tconv(id))

         call freac (vi,eqd,y,dtnuc,.false.)
         call fdiff (eqd,dd,xmod,f1,f2,f3,id,imin,imax)
         
         do j = 1,neqd
            eqd(j,neqd1) = eqd(j,neqd1)-xmod(id,j)+xa(id,j)
            eqd(j,neqd+j) = eqd(j,neqd+j)-1.d0
            do l = 1,neqd
               a(j,l) = eqd(j,neqd+l)
               do k = 1,neqd
                  a(j,l) = a(j,l)-eqd(j,k)*v(i1,k,l)
               enddo
               a(j,neqd+l) = 0.d0
            enddo
            a(j,neqd+j) = 1.d0
         enddo
         call nwrmat (a,um,neqd,neqd,error)
         if (error.gt.0) then
            error = 9
            return
         endif
         do l = 1,neqd
            dv(l) = -eqd(l,neqd1)
            do k = 1,neqd
               dv(l) = dv(l)+eqd(l,k)*w(i1,k)
            enddo
         enddo
         do k = 1,neqd
            w(id,k) = 0.d0
            do l = 1,neqd
               w(id,k) = w(id,k)-um(k,l)*dv(l)
               v(id,k,l) = 0.d0
               do j = 1,neqd
                  v(id,k,l) = v(id,k,l)+um(k,j)*eqd(j,neqd2+l)
               enddo
            enddo
         enddo
         do k = 1,neqd
            do l = neqd+1,neqd2
               eqd(k,l) = 0.d0
            enddo
            eqd(k,neqd1) = 0.d0
            eqd(k,k) = 0.d0
            eqd(k,k+neqd2) = 0.d0
         enddo
      enddo

c..   surface boundary conditions
*--------------------------------

      id = imax
      i1 = id-1

      do j = 1,neqd
         y(j) = xmod(id,j)
      enddo
      do j = 1,nreac
         vi(j) = vv(id,j)
      enddo

      call freac (vi,eqd,y,dtnuc,.false.)
      if (ksurfd.ne.3) call fdiff (eqd,dd,xmod,f1,f2,f3,id,imin,imax)

      do j = 1,neqd
         if (ksurfd.ne.3) then
            eqd(j,neqd1) = eqd(j,neqd1)-xmod(id,j)+xa(id,j)
            eqd(j,neqd+j) = eqd(j,neqd+j)-1.d0
         else
            eqd(j,neqd1) = -xmod(id,j)+xmod(id+1,j)
            eqd(j,neqd+j) = -1.d0
         endif
         do l = 1,neqd
            a(j,l) = eqd(j,neqd+l)
            do k = 1,neqd
               a(j,l) = a(j,l)-eqd(j,k)*v(i1,k,l)
            enddo
            a(j,neqd+l) = 0.d0
         enddo
         a(j,neqd+j) = 1.d0
      enddo
      call nwrmat (a,um,neqd,neqd,error)
      if (error.gt.0) then
         error = 9
         return
      endif
      do j = 1,neqd
         dv(j) = -eqd(j,neqd1)
         do k = 1,neqd
            dv(j) = dv(j)+eqd(j,k)*w(i1,k)
         enddo
      enddo

*______________________________________________________
***   determination of the abundance profile variations
***   (iterative process)
*------------------------------------------------------

      do mm = 1,neqd
         dx(mm) = 0.d0
         do j = 1,neqd
            dx(mm) = dx(mm)+um(mm,j)*dv(j)
         enddo
         xmod(id,mm) = xmod(id,mm)+alpha1(mm)*dx(mm)
      enddo

      cormmax = 0.d0
      ncmmax = 1
      idmax = 2

      do id = imax-1,imin,-1
         do k = 1,neqd
            dv(k) = -w(id,k)
            do l = 1,neqd
               dv(k) = dv(k)-v(id,k,l)*dx(l)
            enddo
         enddo
         do k = 1,neqd
            dx(k) = dv(k)
            xmod(id,k) = xmod(id,k)+alpha1(k)*dx(k)
            xidk = max(abs(xa(id,k)),abs(xmod(id,k)))
            xx = abs(dx(k))/xidk
            if (abs(xmod(id,k)).gt.xtol.and.abs(xa(id,k)).gt.xtol) then
               if (abs(xx).gt.abs(cormmax)) then
                  cormmax = xx
                  dxmax = dx(k)
                  ncmmax = id
                  idmax = k
               endif
            endif
         enddo
      enddo
cc..  check convergence
      if (abs(cormmax).gt.rprecd) then
         icrash = icrash+1
         if (mod(iterd,2).eq.0) then
            dtold = dtnuc
            dtnuc = dtnuc*0.33d0
         endif
         write(nout,'("shell",i4,", dt=",1pe10.4,1x,a5,", dY = ",
     &        1pe13.6,", Y =",1pe13.6,", dY/Y =",1pe13.6)') 
     &        ncmmax,dtnuc,elem(idmax),dxmax,xa(ncmmax,idmax),
     &        cormmax
         fac = dtnuc/dtold
         if (mod(iterd,2).eq.0) then
            do i = imin,imax
c                     do j = 1,nsp
c                        xmod(i,j) = xa(i,j)
c                     enddo
               f1(i) = f1(i)*fac
               f2(i) = f2(i)*fac
               f3(i) = f3(i)*fac
            enddo
         endif
         goto 10
      endif

c..   solution accepted, update abundances and increase time step 
      dt = dt+dtnuc
      if (dt.lt.dtn) then
         istep = istep+1
         dtold = dtnuc
         dtnuc = min(dtnuc*1.78d0,dtn-dt)
c         dtnuc = min(dtnuc*max(1.58d0,cormmax),dtn-dt)
         fac = dtnuc/dtold
         write (nout,100) min(istep,999),cormmax,elem(idmax),ncmmax,
     &        xa(ncmmax,idmax)*anuc(idmax),xmod(ncmmax,idmax)*
     &        anuc(idmax),icrash,dtnuc,1.d2*dt/dtn
         do i = imin,imax
            xmod(i,nis-1) = xmod(i,nis-1)+xmod(i,nis)*anuc(nis)/
     &           anuc(nis-1)
            xmod(i,nis) = 1.d-55
            xmod(i,nsp) = 1.d-55
            do j = 1,nsp
               xa(i,j) = xmod(i,j)
            enddo
            f1(i) = f1(i)*fac
            f2(i) = f2(i)*fac
            f3(i) = f3(i)*fac
         enddo
         iterd = 0
         icrash = 0
         goto 10
      else
c..   time step completed, return               
         do j = imin,imax
            xm_sum = 0.d0
            ielem = 2
            xmod(j,nis-1) = xmod(j,nis-1)+xmod(j,nis)*anuc(nis)/
     &           anuc(nis-1)
            xmod(j,nis) = 1.d-55
            xmod(j,nsp) = 1.d-55
            do k = 1,nsp
               xsp(j,k) = max(xmod(j,k)*anuc(k),1.d-50)
               vxsp(j,k) = xsp(j,k)
               xm_sum = xm_sum+xsp(j,k)
               if (xsp(j,k).gt.xsp(j,ielem)) ielem = k
            enddo
            if (xsp(j,ielem).lt.0.d0) then
               error = 25
               return
            endif            
            if (abs(1.d0-xm_sum).gt.xtol) then
               error = 26
               return
            endif            
            xsp(j,ielem) = xsp(j,ielem)+(1.d0-xm_sum)
            if (.not.nequil) then
               do l = 1,neqd
                  y(l) = xmod(j,l)
               enddo
               do i = 1,nreac
                  vi(i) = vv(j,i)
               enddo
               call yneutron (y,vi,yeq)
               xsp(j,1) = yeq
               vxsp(j,1) = yeq
            endif
            vxsp(j,ielem) = xsp(j,ielem)
         enddo
         return
      endif


***   Manage crash situations

 20   if (ntry.lt.2.and.error.eq.0) then
         ntry = ntry+1
c..   1st try : set equilibrium abundance
         xequil = .true.
c..   2nd try : set free neutron abundance
         if (ntry.eq.2) nequil = .not.nequil
         if (nequil) then
            error = 27
            return
         endif
         write (nout,*) iterd,dtnuc,xequil,nequil
         iterd = 0
         dtold = dtnuc
         dt = 0.d0
         dtnuc = dtn*0.125d0
         fac = dtnuc/dtold
         do k = imin,imax
            f1(k) = f1(k)*fac
            f2(k) = f2(k)*fac
            f3(k) = f3(k)*fac
            do j = 1,nsp
               xa(i,j) = xmod(i,j)
            enddo
         enddo
         goto 10
      else
         error = 27
         return
      endif


 100  format(i3,"| dY/Y=",1pe10.3,", for ",a5,", shell =",i4,", X0 =",
     &     1pe10.3,", X =",1pe10.3,", ifail=",i3,", dt=",1pe9.3,
     &     ", time=",0pf10.7,' %')

      return
      end


************************************************************************

      SUBROUTINE network (cshell,xi,mconv,dt,tnucmin,imin,imax,error)

************************************************************************
* Calculate the chemical evolution during the current model time-step  *
* cshell = c   : convective zone
* cshell = r   : radiative zone
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.ion'
      include 'evolcom.mod'
      include 'evolcom.nuc'
c      include 'evolcom.conv'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer error,imin,imax,ish
      integer normx,k
      integer i,j,l

      double precision tk,rok,mconv
      double precision xi,dt,y0,ysav,v,vi,irrad0,
     &     yeqmean,dmam,zlocal
      double precision dnucc,yeq,sumxsp,sen,tmean,tnucmin

      logical ipass,nsink

      character cshell*1

      dimension xi(nsp),y0(nsp),ysav(nsp),v(nreac),vi(nreac)

      ipass = .false.
      nsink = .true.

***   special treatment of neutrons for extremely metal poor stars
      if (zkint.lt.5.d-6) then
         zlocal = 0.d0
         do i = 2,ib11
            zlocal = zlocal+xi(i)
         enddo
         zlocal = 1.d0-zlocal
         if (zlocal.lt.1.d-10) then
            k7(iddn) = 2        !  2D(g,n)He3   --> 2D(g,p)He3
            k7(ili7d) = 2       !  Li7(D,n)2He4 --> Li7(D,p)2He4
            nsink = .false.
         else
            k7(iddn) = 1        !  2D(g,n)He3
            k7(ili7d) = 1       !  Li7(D,n)2He4
         endif
      endif

***   set up nuclear network tolerence in advanced evolutionary stages
      if (cshell.eq.'r') then
         ish = imin
         tmean = t(ish)
      else
         tmean = dsqrt(t(imin)*t(imax))
      endif
      if (tmean.gt.4.d8) then
         if (tmean.gt.1.2d9) then
            dnucc = min(1.d-5,dnuc*1.d4)
            ytminc = min(1.d-7,ytmin*1.d4)
         else
            dnucc = min(1.d-7,dnuc*1.d2)
            ytminc = min(1.d-9,ytmin*1.d2)
         endif
      else
         dnucc = dnuc
         ytminc = ytmin
      endif

***   initialisations

      irrad0 = 0.d0
      do l = 1,nsp
c         x(l) = xi(l)
         y0(l) = xi(l)/anuc(l)
         ysav(l) = y0(l)
      enddo

*_____________________________________________________________________
***   determine neutron equilibrium abundance, temperature and density
*---------------------------------------------------------------------

***   radiative shell
      if (cshell.eq.'r') then
         tk = (t(ish)+vt(ish))*0.5d0
         rok = (ro(ish)+vro(ish))*0.5d0
         call vit (tk,rok,mueinv(ish),ish,v,0)
         if (nsink) then 
            call yneutron (y0,v,yeq)
            y0(1) = yeq
            ysav(1) = yeq
         endif
         irrad0 = 7.68771025912336d0*yeq*rok*dsqrt(tk)
      endif

***   convective shell
      if (cshell.eq.'c') then
         do i = 1,nreac
            v(i) = 0.d0
         enddo
         yeqmean = 1.d-50
         do k = imin,imax
            tk = (vt(k)+t(k))*0.5d0
            rok = (vro(k)+ro(k))*0.5d0
            dmam = dm(k)*mconv
            call vit (tk,rok,mueinv(k),k,vi,0)
            if (nphase.le.2.or.tk.gt.tnucmin) then
               if (nsink) then 
                  call yneutron (y0,vi,yeq)
                  y0(1) = yeq
               endif
               irrad0 = irrad0+7.68771025912336d0*yeq*rok*dsqrt(tk)*dmam
               do j = 1,nreac
                  if (k1(j).eq.1.or.k3(j).eq.1) then
                     v(j) = v(j)+vi(j)*dmam*y0(1)
                  else
                     v(j) = v(j)+vi(j)*dmam
                  endif
               enddo
            else
               do i = 1,nbeta
                  j = jbeta(i)
                  v(j) = v(j)+vi(j)*dmam
               enddo
            endif
            yeqmean = yeqmean+y0(1)*dmam
         enddo
         y0(1) = yeqmean
         ysav(1) = y0(1)
         do j = 1,nreac
            if (k1(j).eq.1.or.k3(j).eq.1) v(j) = v(j)/yeqmean
         enddo
      endif
      irrad0 = irrad0*dt

*_______________________________________
***   calculation of the nucleosynthesis
*---------------------------------------

 10   call nucsolve (y0,v,dt,irrad0,nsink,imin,error)
      if (error.gt.0) return

      y0(nis-1) = y0(nis-1)+y0(nis)*anuc(nis)/anuc(nis-1)
      y0(nis) = 1.d-55
      y0(nsp) = 1.d-55
      sumxsp = 0.d0
      normx = 1
      do i = 1,nsp
         xi(i) = max(anuc(i)*y0(i),1.d-50)
         if (xi(i).gt.xi(normx)) normx = i
         sumxsp = sumxsp+xi(i)
      enddo
      sen = 1.d0-sumxsp
      if (abs(sen).gt.dnucc) then
         if (.not.ipass) then
            do l = 1,nsp
               y0(l) = ysav(l)
            enddo
            ytminc = ytminc*1.d2
            dnucc = dnucc*1.d1
            error = -1
            ipass = .true.
            goto 10
         else
            write (nout,1100) imin,sen,dnucc,ytminc,tk*1.d-9
            error = 3
         endif
      else
c..   convergence and normalization
         xi(normx) = xi(normx)+sen
c..   compute neutron irradiation
c..   irrad = y0(1)*rok*avn*dt*dsqrt(2.d0*boltz*tk*avn*0.98d0)*1.d-27
         if (cshell.eq.'c') then
            do k = imin,imax
               exposure(k) = exposure(k)+7.68771025912336d0*ro(k)*
     &              dsqrt(t(k))*irrad0
            enddo
         else
            exposure(imin) = exposure(ish)+7.68771025912336d0*ro(ish)*
     &           dsqrt(t(ish))*irrad0
         endif
      endif


 1000 format (/,5x,'shell #',i4,', 1-sumxsp = ',1pd13.6,' > dnucc = ',
     &     1pd12.6,' with ytminc = ',1pd10.3,', t9 = ',1pd9.3,/,5x,
     &     'Try again with smaller tolerance : ytminc = ',1pd10.3,/)
 1100 format (/,5x,'shell #',i4,', 1-sumxsp = ',1pd13.6,' > dnucc = ',
     &     1pd12.6,' with ytminc = ',1pd10.3,', t9 = ',1pd9.3)


      return
      end


************************************************************************

      SUBROUTINE neutri (tk,rok,mueinvk,x,ksh,q,dqt,dqro)

************************************************************************
* Calculate the neutrino loss rates due to the nuclear reactions       *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.eng'
      include 'evolcom.nuc'
      include 'evolcom.spec'

      integer i,j,i1,i2,i4
      integer ksh,iroe,ite

      double precision tk,rok,x,mueinvk
      double precision q,q2,qrho,qrho2,dqt,dqro
      double precision enu,enu2,v,tkp,dtk,rokp,drok
      double precision ee,ee2,eerho,eerho2,eecap,eecap2
      double precision a1,a2,b1,b2,y
      double precision vcapte,enuelec,rhoye,tcapte,qinuelec

      common /electcap/ rhoye(nroe),tcapte(nte),vcapte(nte,nroe,nrce),
     &     enuelec(nte,nroe,nrce)
      common /electcap_par/ ite,iroe,a1,a2,b1,b2

      dimension enu(nreac),enu2(nreac),qinuelec(nrce),v(nreac),x(nsp),
     &     y(nsp)

      ee = 0.d0
      ee2 = 0.d0
      eerho = 0.d0
      eerho2 = 0.d0
      eecap = 0.d0
      eecap2 = 0.d0
      do i = 1,nsp
         y(i) = x(i)/anuc(i)
      enddo
***   reactions independent of rho*Ye
      do j = 1,nre
         i2 = k2(j)
         i4 = k4(j)
         enu(j) = flu(ksh,j)*qinu(j)
         enu2(j) = flu2(j)*qinu(j)
         eerho = eerho+dble(i2+i4-1)*enu(j)
         ee = ee+enu(j)
         ee2 = ee2+enu2(j)
      enddo

***   rho*Ye-T dependent reactions (electron captures only)

***   rho,T+dT
      dtk = tk*1.d-4
      tkp = tk+dtk
      call vit (tkp,rok,mueinvk,ksh,v,2)
      do j = nre+1,nreac
         i = j-nre
         qinuelec(i) = a1*a2*enuelec(ite+1,iroe+1,i)+b1*b2*
     &        enuelec(ite,iroe,i)+a1*b2*enuelec(ite+1,iroe,i)+
     &        b1*a2*enuelec(ite,iroe+1,i)
      enddo
      do j = nre+1,nreac
         i1 = k1(j)
         i = j-nre
         flu2(j) = y(i1)*v(j)
         eecap2 = eecap2+flu2(j)*qinuelec(i)
      enddo

***   rho,T
      call vit (tk,rok,mueinvk,ksh,v,2)
      do j = nre+1,nreac
         i = j-nre
         qinuelec(i) = a1*a2*enuelec(ite+1,iroe+1,i)+b1*b2*
     &        enuelec(ite,iroe,i)+a1*b2*enuelec(ite+1,iroe,i)+
     &        b1*a2*enuelec(ite,iroe+1,i)
      enddo
      do j = nre+1,nreac
         i1 = k1(j)
         i = j-nre
         flu(ksh,j) = y(i1)*v(j)
         eecap = eecap+flu(ksh,j)*qinuelec(i)
      enddo

***   rho+drho,T
      drok = rok*1.d-4
      rokp = rok+drok
      call vit (tk,rokp,mueinvk,ksh,v,2)
      do j = nre+1,nreac
         i = j-nre
         qinuelec(i) = a1*a2*enuelec(ite+1,iroe+1,i)+b1*b2*
     &           enuelec(ite,iroe,i)+a1*b2*enuelec(ite+1,iroe,i)+
     &        b1*a2*enuelec(ite,iroe+1,i)
      enddo
      do j = nre+1,nreac
         i1 = k1(j)
         i = j-nre
         eerho2 = eerho2+y(i1)*v(j)*qinuelec(i)
      enddo
      ee = ee+eecap
      ee2 = ee2+eecap2
      q = ee*econv
      q2 = ee2*econv
      dqt = (q2-q)/dtk
      qrho = eecap*econv
      qrho2 = eerho2*econv
      dqro = eerho*econv/rok+(qrho2-qrho)/drok

      return
      end


************************************************************************

      SUBROUTINE nuceng (nflu,y,v,ksh,ee,eerho,eerhocap)

************************************************************************
* Calculate the nuclear energy production rate                         *
* nflu =< 0 : rhot,T
* nflu =  1 : rhot,T+dT
* nflu =  2 : rhot+drho,T
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.nuc'
      include 'evolcom.spec'

      integer nflu,ksh,ite,iroe
      integer i1,i2,i3,i4
      integer i,j,mm

      double precision ee,eerho,eerhocap,efn,flux,qinucap
      double precision v(nreac),vi(nreac),y(nsp)
      double precision a1,a2,b1,b2
      double precision vcapte,enuelec,rhoye,tcapte

      common /electcap/ rhoye(nroe),tcapte(nte),vcapte(nte,nroe,nrce),
     &     enuelec(nte,nroe,nrce)
      common /electcap_par/ ite,iroe,a1,a2,b1,b2

      common /engen/ efn(nreac),flux(nreac)

      ee = 0.d0
      eerho = 0.d0
      eerhocap = 0.d0
      do i = 1,nreac
         vi(i) = v(i)
      enddo

***   T dependence
*  reaction : k2  X(k1) + k4 X(k3) --> k8 X(k7) + k6 X(k5)
*        ex : 1    C12  + 1   He   --> 0  gamma + 1   O16

      if (nflu.lt.2) then
         do mm = 1,nreac
            i1 = k1(mm)
            i2 = k2(mm)
            i3 = k3(mm)
            i4 = k4(mm)
            if (i2.eq.1) then
               if (i4.eq.1) then
                  flux(mm) = vi(mm)*y(i1)*y(i3)
               else
                  flux(mm) = vi(mm)*y(i1)
               endif
            else
               flux(mm) = vi(mm)*y(i1)**i2/fact(i2)
            endif
            efn(mm) = flux(mm)*(qi(mm)-qinu(mm))
            ee = ee+efn(mm)
            if (nflu.le.0) then
               if (mm.le.nre) then
                  eerho = eerho+dble(i2+i4-1)*efn(mm)
               else
                  eerhocap = eerhocap+efn(mm)
               endif
            endif
         enddo
***   rho dependence
      else
         do j = nre+1,nreac
            i = j-nre
            qinucap = a1*a2*enuelec(ite+1,iroe+1,i)+b1*b2*
     &           enuelec(ite,iroe,i)+a1*b2*enuelec(ite+1,iroe,i)+
     &           b1*a2*enuelec(ite,iroe+1,i)
            i1 = k1(j)
            eerhocap = eerhocap+vi(j)*y(i1)*(qi(j)-qinucap)
         enddo
      endif
      ee = ee*econv
      eerho = eerho*econv
      eerhocap = eerhocap*econv

      return
      end


************************************************************************

      SUBROUTINE nucsolve (y,v,dtn,irrad0,nsink,ksh,error)

************************************************************************
*     Solve the nucleosynthesis equations
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.nuc'

      integer i,j,l
      integer istart,error,normx,negflag
      integer neqd,neqd1,ksh,ntry
      integer iterd,itermax
      integer indx(nis),info

      double precision d,suma,sumd,sumx
      double precision yeq,dtn,dtnuc,dt,irrad0
      double precision v(nreac),eqd(nsp,3*nsp+1)
      double precision y(nsp),x(nsp),y0(nsp),ysav(nsp)
      double precision p(nis),fjac(nis,nis)

      integer nz,lirn,ina,jna,licn
      parameter (nz = 370,lirn=nsp,licn=nz+1)

      logical nequil,xequil,nsink,normalize

      common /nuc_matrix/ jna(nz),ina(lirn)

      ntry = 0
      dt = 0.d0
      if (irrad0.lt.1.d-15) then
         dtnuc = dtn
      else
         dtnuc = dtn*5.d-2
      endif
      irrad0 = 0.d0
      itermax = 10000
      neqd1 = 3*nsp+1
      neqd = nis

      normalize = .false.
      xequil = .false.
      nequil = .true.
      yeq = 1.d-50

c..   for extremely metal poor stars (Z = 0), no neutrons sink
c      if (.not.nsink) nequil = .false.

      if (error.eq.-1) then
         nequil = .false.
         xequil = .true.
      endif
      ysav(1:nsp) = y(1:nsp)
      y0(1:nsp) = y(1:nsp)
      eqd(1:nsp,1:neqd1) = 0.d0
      if (nequil) then
         istart = 2
      else
         istart = 1
      endif

c..   set neutrons to equilibrium abundance
      if (nequil) then
         call yneutron (y,v,yeq)
         y(1) = yeq
      endif

      iterd = 0
 10   iterd = iterd+1
      if (iterd.ge.itermax.or.dtnuc.lt.1.d-10) goto 20

c..   set  equilibrium abundance for elements with y < yeqmin
      if (xequil) call yequil (v,y,ytminc,dtnuc,1.d99)
      call freac (v,eqd,y,dtnuc,.true.)
      do j = 1,neqd
c.withp         p(j) = y(j)-y0(j)-eqd(j,neqd1)
         p(j) = -eqd(j,neqd1)
         do l = 1,neqd
            fjac(l,j) = eqd(l,j)
         enddo
         fjac(j,j) = fjac(j,j)-1.d0
      enddo
c.. leqs
c      call leqs(fjac,p,neqd,neqd)
c.. numerical recipee
c      call ludcmp(fjac,neqd,neqd,indx,d)
c      call lubksb(fjac,neqd,neqd,indx,p)
c.. lapack
cx!$OMP PARALLEL PRIVATE(info)
      call dgetrf (neqd,neqd, fjac, neqd, indx, info)
      call dgetrs ('N',neqd,1,fjac,neqd,indx,p,neqd,info)
cx!$OMP END PARALLEL
c..  sparse matrix
c      p = 1.d0
c      j = 1
c      do i = 1,nz
c         if (i.eq.ina(j+1)) j = j+1
c         print *,i,j,jna(i)
c         fjac_s(i) = fjac(j,jna(i))
c      enddo
c      call ma28ad(neqd,nz,fjac_s,nz,ina,lirn,jna,d,indx_s,iw,w,info)
c      call ma28cd(neqd,fjac_s,nz,jna,indx_s,p,w,1)
c.. gift
c      do j = 1,neqd
c         fjac(j,neqd+1) = p(j)
c      enddo
c      call gift_test(fjac,neqd,neqd+1)
c      do j = 1,neqd
c         p(j) = fjac(j,neqd+1)
c      enddo
      sumd = 0.d0
      suma = 0.d0
      do i = istart,neqd
         sumd = sumd + p(i)*anuc(i) 
         suma = suma + abs(p(i))*anuc(i)
      enddo
      negflag = 0
      do i = istart,neqd
         y(i) = y(i)+p(i)
         if (y(i).lt.-2.d16) negflag = i
c..  check linearization approximation is valid : dY/dt * dtnuc << Y
         if ((abs(p(i)/y0(i)).gt.0.1d0.and.y(i).gt.ytminc.and.
     &        y0(i).gt.ytminc).or.negflag.gt.0.or.suma.ge.5.d-2.or.
     &        abs(sumd).ge.1.d-10) then
            dtnuc = dtnuc*0.23d0
            do l = 1,nsp
               y(l) = y0(l)
               do j = 1,nsp
                  eqd(j,l) = 0.d0
               enddo
               eqd(l,neqd1) = 0.d0
            enddo
            goto 10
         endif
      enddo

c..   solution accepted, update abundances and increase time step 
      dt = dt+dtnuc
      if (dt.lt.dtn) then
         do i = 1,nsp
            y0(i) = max(y(i),1.d-50)
            do j = 1,nsp
               eqd(j,i) = 0.d0
            enddo
            eqd(i,neqd1) = 0.d0
         enddo
         if (normalize) then
            sumx = 0.d0
            normx = 2
            do i = 1,nsp
               x(i) = max(anuc(i)*y(i),1.d-50)
               if (x(i).gt.x(normx)) normx = i
               sumx = sumx+x(i)
               y(normx) = (x(normx)+1.d0-sumx)/anuc(normx)
               y0(normx) = y(normx)
            enddo
         endif
c..   set neutrons to equilibrium abundance
         if (nequil) then
            call yneutron (y,v,yeq)
            y(1) = yeq
         endif
         irrad0 = irrad0+dtnuc*yeq
         dtnuc = min(dtnuc*1.58d0,dtn-dt)
         goto 10
      else
c..   time step completed, return
         if (nsink) then
            call yneutron (y,v,yeq)
            y(1) = yeq
         endif
         irrad0 = irrad0+dtnuc*yeq
         error = 0
         return
      endif


***   Manage crash situations

 20   if (ntry.lt.2.and.error.eq.0) then
         ntry = ntry+1
c..   1st try : set equilibrium abundance
         xequil = .true.
c..   2nd try : set free neutron abundance
         if (ntry.eq.2) nequil = .not.nequil
         if (.not.nsink.and.nequil) then
            error = 14
            return
         endif
         if (nequil) then
            istart = 2
         else
            istart = 1
         endif
         write (nout,100) ksh,iterd,dtnuc,xequil,nequil
         iterd = 0
         dt = 0.d0
         dtnuc = dtn*0.125d0
         y(1:nsp) = ysav(1:nsp)
         eqd(1:nsp,1:nsp) = 0.d0
         eqd(1:nsp,neqd1) = 0.d0
         goto 10
      else
         error = 14
         return
      endif

 100  format (' NUCSOLVE problem at shell #',i4,' : performed ',i5,
     &     ' iterations, final time step = ',1pe9.3,
     &     '(s), new try with Xequil = ',l1,', Nequil = ',l1)

      return
      end


************************************************************************

      SUBROUTINE nucprod(v,y,tk,epp,eCNO,e3a,eC12,eNe20,eO16)

************************************************************************
* Compute energy generation rate for specific reactions                *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.nuc'

      double precision epp,eCNO,e3a,eC12,eNe20,eO16
      double precision ylim,tk
      double precision v(nreac),y(nsp)

      ylim = 1.d-8
      epp = 0.d0
      eCNO = 0.d0
      e3a = 0.d0
      eC12 = 0.d0
      eNe20 = 0.d0
      eO16 = 0.d0

c PP chains
c.. exclude beta decay reactions

c      if (y(ih1).ge.ylim.or.nmixd.gt.0) then
      if (y(ih1).ge.ylim) then
         epp = y(ih1)**2*v(ippg)*qi(ippg)*0.5d0+
     &   y(ih1)*y(ih2)*v(ipdg)*qi(ipdg)+
     &   y(ihe3)**2*v(i2he3)*qi(i2he3)*0.5d0+
     &   y(ihe4)*y(ihe3)*v(ihe3ag)*qi(ihe3ag)+
c     &   y(ibe7)*v(ibe7beta)*qi(ibe7beta)+
     &   y(ili7)*y(ih1)*v(ili7pa)*qi(ili7pa)+
     &   y(ibe7)*y(ih1)*v(ibe7pg)*qi(ibe7pg)
c     &   y(ib8)*v(ib8beta)*qi(ib8beta)
         epp = epp*econv

c CNO bicycle

         eCNO = y(ic12)*y(ih1)*v(icpg)*qi(icpg)+
c     &   y(in13)*v(in13beta)*qi(in13beta)+
     &   y(ic13)*y(ih1)*v(ic13pg)*qi(ic13pg)+
     &   y(in14)*y(ih1)*v(in14pg)*qi(in14pg)+
c     &   y(io15)*v(io15beta)*qi(io15beta)+
     &   y(in15)*y(ih1)*v(in15pa)*qi(in15pa)+
     &   y(in15)*y(ih1)*v(in15pg)*qi(in15pg)+
     &   y(io16)*y(ih1)*v(io16pg)*qi(io16pg)+
     &   y(io17)*y(ih1)*v(io17pa)*qi(io17pa)
         eCNO = eCNO*econv
      endif

c 3a and C12(a,g)
      if (y(ihe4).gt.ylim) then
         e3a = y(ihe4)**3*v(i3a)*qi(i3a)/6.d0+
     &   y(ihe4)*y(ic12)*v(icag)*qi(icag)
         e3a = e3a*econv
      endif

c Carbon burning 12C(12C,g)24Mg is not taken into account in the network
      if (y(ic12).gt.ylim.and.tk.gt.3.7d8) then
         eC12 = y(ic12)**2*v(iccga)*qi(iccga)*0.5d0+
     &        y(ic12)**2*v(iccgp)*qi(iccgp)*0.5d0
         eC12 = eC12*econv
      endif

c Neon photodisintegration/Burning
      if (y(ine20).gt.ylim.and.tk.gt.7.d8) then
         eNe20 = y(ine20)*v(inega)*qi(inega)
c         if (eNe20.gt.shlim.or.nphase.eq.7) then
c            eNe20 = eNe20+y(ine20)*y(ihe4)*v(ine20ag)*qi(ine20ag)+
c     &      y(img24)*y(ihe4)*v(img24ag)*qi(img24ag)+
c     &      y(isi28)*y(ihe4)*v(isi28ag)*qi(isi28ag)
c         endif
         eNe20 = eNe20*econv
      endif

c Oxygen Burning
      if (y(io16).gt.ylim.and.tk.gt.1.d9) then
         eO16 = 0.5d0*y(io16)**2*(v(iooga)*qi(iooga)+v(ioogp)*qi(ioogp))
         eO16 = eO16*econv
      endif

      return
      end


************************************************************************

      SUBROUTINE nwracc (corm,iter,alpha1)

************************************************************************
* Accelerate the convergence of the Henyey's method                    *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.mod'

      integer iter
      integer niter,ncor,nav,nsame,nrecover,nsamemax
      integer k,l,ii

      double precision corm,alpha1
      double precision co,lnco,al
      double precision newalpha1
      double precision fac,f1,f2,facinv,ff1,f1f,ff2,f2f
      double precision f3,f4
      double precision eps,med
      double precision av,last,diflast,difmed,alp1,alp2,sum,a,b,c,d,e
      double precision factor1,factor2
      double precision step
      double precision s1,s2,s3,s4,s5
      double precision pic
      double precision p1,p2,p3,p4

      logical tochange,interpth
      logical nowait

      parameter (niter = 4, ncor = niter-1, nsamemax = 4)
c     parameter (f3 = 0.45d0, f4 = 1/dsqrt(f3))
      parameter (f3 = 0.45d0, f4 = 1.490712d0)
      parameter (fac = 3.8d0, f1 = 2.d0, f2 = 4.d0)
      parameter (facinv = 1.d0/fac, f1f = f1/fac, ff1 = fac/f1)
      parameter (f2f = f2/fac, ff2 = fac/f2)
      parameter (eps = 0.5d0)
      parameter (nowait = .true.)

      common /accelnwr/ co(neq,niter),lnco(neq,niter),
     &     al(neq,niter),nsame(neq),nrecover(neq)

      dimension corm(neq),alpha1(neq),newalpha1(neq)
      dimension tochange(neq)

      save tochange

      step (s1,s2,s3,s4,s5) = s1/(exp(s2*s5**3)+s3)+s4
      pic (p1,p2,p3,p4) = p1*exp(-p4*p4/p2)+p3

      if (iter.le.1) then
         do k = 1,neq
            nsame(k) = 1
            nrecover(k) = 0
            do l = 1,niter
               co(k,l) = 0.d0
               lnco(k,l) = 0.d0
               al(k,l) = 0.d0
            enddo
            tochange(k) = .false.
         enddo
      endif

***   Test for an interpulse period (AGB phase)
      interpth = (agbphase.and.(corm(5).lt.1.d0))

      ii = 2
      if (hydro) ii = 1
      do l = niter,2,-1
         do k = ii,neq
            co(k,l) = co(k,l-1)
            lnco(k,l) = lnco(k,l-1)
            al(k,l) = al(k,l-1)
         enddo
      enddo
      do k = ii,neq
         co(k,1) = corm(k)
         lnco(k,1) = log(abs(corm(k)))
         al(k,1) = alpha1(k)
      enddo

      if (iter.ge.ncor) then
         if (iter.ge.niter) then
            nav = niter
         else
            nav = ncor
         endif
         do k = ii,neq
            tochange(k) = .not.tochange(k)
            sum = 0.d0
            do l = 1,niter
               sum = sum+lnco(k,l)
            enddo
            med = lnco(k,1)-lnco(k,ncor)
            av = sum/nav
            last = lnco(k,1)
            diflast = av-last
            difmed = av-med
            if (difmed.le.(-abs(eps))) then
               c = (ff1-fac)/(facinv-ff1)
               a = (facinv-fac)*c
               d = fac
               b = 1.d0
            endif
            if (abs(difmed).lt.abs(eps)) then
               c = (1.d0-ff2)/(f2f-1.d0)
               a = (f2f-ff2)*c
               b = 1.d0
               d = ff2
            endif
            if (difmed.ge.abs(eps)) then
               c = (f1f-ff1)/(facinv-f1f)
               a = (facinv-ff1)*c
               b = 1.d0
               d = ff1
            endif
            e = min(diflast,3.d0)
            factor1 = step (a,b,c,d,e)
            if (diflast.le.(-abs(eps))) then
               if (difmed.gt.0.d0) then
                  a = f2f-facinv
                  b = 8.d0
                  c = facinv
                  d = 0.d0
                  factor2 = pic (a,b,c,difmed)
               endif
               if (difmed.le.0.d0) then
                  a = f2f-f1f
                  b = 8.d0
                  c = f1f
                  d = 0.d0
                  factor2 = pic (a,b,c,difmed)
               endif
            endif
            if (abs(diflast).lt.abs(eps)) then
               c = (1.d0-f1f)/(ff1-1.d0)
               a = (ff1-f1f)*c
               b = 1.d0
               d = f1f
               e = min(difmed,3.d0)
               factor2 = step (a,b,c,d,e)
            endif
            if (diflast.ge.abs(eps)) then
               if (difmed.gt.0.d0) then
                  a = ff2-ff1
                  b = 8.d0
                  c = ff1
                  d = 0.d0
                  factor2 = pic (a,b,c,difmed)
               endif
               if (difmed.le.0.d0) then
                  a = ff2-fac
                  b = 8.d0
                  c = fac
                  d = 0.d0
                  factor2 = pic (a,b,c,difmed)
               endif
            endif
            alp1 = alpha1(k)*factor1
            alp2 = alpha1(k)*factor2
            if (interpth) then
               newalpha1(k) = max(alp1,alp2)
            else
               newalpha1(k) = 0.5d0*(alp1+alp2)
            endif
            newalpha1(k) = min(newalpha1(k),alphmax)
            newalpha1(k) = max(newalpha1(k),alphmin)
** 3 consecutive decrease in the corrections
            if (lnco(k,1).lt.lnco(k,2).and.lnco(k,2).lt.lnco(k,3).and.
     &           newalpha1(k).le.alpha1(k))
     &           newalpha1(k) = min(1.2d0*alpha1(k),alphmax)
** 3 consecutive increase in the corrections
            if (lnco(k,1).gt.lnco(k,2).and.lnco(k,2).gt.lnco(k,3).and.
     &           newalpha1(k).ge.alpha1(k))
     &           newalpha1(k) = max(f3*alpha1(k),alphmin)
** Try to avoid a cycling in the luminosity
c           if (k.eq.5) then
               if (newalpha1(k).eq.alpha1(k)) then
                  nsame(k) = nsame(k)+1
               else
                  nsame(k) = 0
               endif
               if (nsame(k).eq.nsamemax) then
                  newalpha1(k) = max(f3*alpha1(k),alphmin)
                  nrecover(k) = 2
               endif
               if (nrecover(k).gt.0.and.nsame(k).ne.nsamemax) then
                  newalpha1(k) = min(f4*alpha1(k),alphmax)
                  nrecover(k) = nrecover(k)-1
               endif
c           endif
            if (nowait) then
               alpha1(k) = newalpha1(k)
            else
               if (tochange(k)) then
                  alpha1(k) = newalpha1(k)
               else
                  alpha1(k) = alpha1(k)
               endif
            endif
         enddo
      endif

      return
      end


************************************************************************

      SUBROUTINE nwtraf (*,error)

************************************************************************
* Apply the Newton-Raphson method for the five equations               *
* of the stellar structure evolution                                   *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
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
      include 'evolcom.ion'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'

      integer iiter
      integer ncm(neq)
      integer neq2
      integer iitmax,i1,imerk,icrasht
      integer klenv,klpulse,klcore
      integer i,j,k,kl,l,mm
      integer error
      integer itmax,ishockt,ishockb

      parameter (neq2 = 2*neq)

      double precision xmod,xa,eps
      double precision t,vt,r,vr
      double precision dx(neq),dv(neq),w(nsh,neq),a(neq,neq2),
     &     um(neq,neq),v(nsh,neq,neq),corm(neq),alpha1(neq)
      double precision xishk,xx
      double precision Lmax,tmax,mackmax,enucmax
      double precision dp,ag

      logical corrmax

      common /calcul/ iiter
      common /convergence/ imerk,icrasht
      common /overshoot/ klcore,klenv,klpulse
      common /varia/ xmod(nsh,neq),xa(nsh,neq),eps(neq)
      common /vartvt/ t(nsh),vt(nsh)
      common /varrvr/ r(nsh),vr(nsh)

      common /hydrodyn/ Lmax,tmax,mackmax,enucmax,itmax,ishockb,ishockt

      iitmax = itermax

      corrmax = .true.
      imerk = -1
      iter = 0
      error = 0

 10   iter = iter+1
      iiter = iter
      if (iter.le.1) then
         do i = 1,neq
            alpha1(i) = alpha*0.5d0
         enddo
      endif
      if (iter.eq.2) then
         do i = 1,neq
            alpha1(i) = alpha
         enddo
      endif
c      print *, 'dtn nwtraf',dtn*seci
*________________________________________________
***   calculation of all the constitutive physics
*------------------------------------------------
      
      call thermo (error)
      if (error.gt.0) return 1
      
*___________________________________________________
***   application of the central boundary conditions
*---------------------------------------------------

      ish = 1

      call centeq

      a(1:neq,1:neq) = eq(1:neq,1:neq)
      a(1:neq,neq+1:neq2) = 0.d0
      do j = 1,neq
c         do mm = 1,neq
c            a(mm,j) = eq(mm,j)
c            a(mm,neq+j) = 0.d0
c         enddo
         a(j,neq+j) = 1.d0
      enddo
      call nwrmat (a,um,neq,neq,error)
      if (error.gt.0) then
         write (nout,'("error nwtraf, central shell 1",i4)') crz(ish)
         do l = 1,neq
            do i = 1,neq
               write (nout,'(i2,i2,3(1x,1pe10.3))') i,l,eq(i,l),
     &              eq(i,neq+l),eq(i,2*neq+l)
            enddo
         enddo
         return 1
      endif
      do mm = 1,neq
         w(ish,mm) = 0.d0
         do k = 1,neq
            w(ish,mm) = w(ish,mm)+um(mm,k)*eq(k,neq1)
            v(ish,mm,k) = 0.d0
            do j = 1,neq
               v(ish,mm,k) = v(ish,mm,k)+um(mm,j)*eq(j,neq+k)
            enddo
         enddo
      enddo
      
*__________________________________________________
***   determination of the internal shell structure
*--------------------------------------------------

      do ish = 2,nmod1
         i1 = ish-1

         call inteq

         do j = 1,neq
            do l = 1,neq
               a(j,l) = eq(j,neq+l)
               do k = 1,neq
                  a(j,l) = a(j,l)-eq(j,k)*v(i1,k,l)
               enddo
               a(j,neq+l) = 0.d0
            enddo
            a(j,neq+j) = 1.d0
         enddo
         call nwrmat (a,um,neq,neq,error)
         if (error.gt.0) then
            write (nout,'("error nwtraf, shell #",i4,i3)') ish,crz(ish)
            do l = 1,neq
               do i = 1,neq
                  write (nout,'(i2,i2,3(1x,1pe10.3))') i,l,eq(i,l),
     &                 eq(i,neq+l),eq(i,2*neq+l)
               enddo
            enddo
            return 1
         endif
         do l = 1,neq
            dv(l) = -eq(l,neq1)
            do k = 1,neq
               dv(l) = dv(l)+eq(l,k)*w(i1,k)
            enddo
         enddo
         do k = 1,neq
            w(ish,k) = 0.d0
            do l = 1,neq
               w(ish,k) = w(ish,k)-um(k,l)*dv(l)
               v(ish,k,l) = 0.d0
               do j = 1,neq
                  v(ish,k,l) = v(ish,k,l)+um(k,j)*eq(j,neq2+l)
               enddo
            enddo
         enddo
      enddo

*___________________________________________________
***   application of the surface boundary conditions
*---------------------------------------------------

      ish = nmod
      i1 = ish-1
      
      call surfeq
      
      do j = 1,neq
         do l = 1,neq
            a(j,l) = eq(j,neq+l)
            do k = 1,neq
               a(j,l) = a(j,l)-eq(j,k)*v(i1,k,l)
            enddo
            a(j,neq+l) = 0.d0
         enddo
         a(j,neq+j) = 1.d0
      enddo
      call nwrmat (a,um,neq,neq,error)
      if (error.gt.0) then
         write (nout,'("error nwtraf, surface #",i4,a1)') ish,crz(ish)
         do l = 1,neq
            do i = 1,neq
               write (nout,'(i2,i2,3(1x,1pe10.3))') i,l,eq(i,l),
     &              eq(i,neq+l),eq(i,2*neq+l)
            enddo
         enddo
         return 1
      endif
      do j = 1,neq
         dv(j) = -eq(j,neq1)
         do k = 1,neq
            dv(j) = dv(j)+eq(j,k)*w(i1,k)
         enddo
      enddo
      
*_______________________________________________________________
***   determination of the evolved structure (iterative process)
*---------------------------------------------------------------

      do mm = 1,neq
         dx(mm) = 0.d0
         do j = 1,neq
            dx(mm) = dx(mm)+um(mm,j)*dv(j)
         enddo
      enddo
      do i = 1,neq
         corm(i) = 0.d0

         if (xmod(ish,i).ne.0.d0.and.((numeric.eq.1.and.t(ish).gt
     &        .tmaxioHe).or.numeric.ne.1)) corm(i) = dx(i)/xmod(ish,i)
c... norbert
c         if (i.ne.5.and.xmod(ish,i).ne.0.d0) corm(i) = dx(i)/xmod(ish,i)
c         if (i.eq.5) corm(i) = dx(i)
c         epsmax(5) = 1.d50
         ncm(i) = ish

cccc---------------------------------------------
ccc TEST : small corrections in the upper layers
c          if (t(ish).lt.1.d5) then
c             do l = 1,neq
c                alpha1(l) = 0.1d0
c             enddo
c          endif
cccc---------------------------------------------
         xmod(ish,i) = xmod(ish,i)+alpha1(i)*dx(i)
      enddo
      do ish = nmod1,1,-1
cccc---------------------------------------------
ccc TEST : small corrections in the upper layers
C          if (t(ish).lt.1.d5) then
C             do l = 1,neq
C                alpha1(l) = 0.1d0
C             enddo
C          else
C             do l = 1,neq
C                alpha1(l) = alpha2(l)
C             enddo
C          endif
cccc---------------------------------------------
         do k = 1,neq
            dv(k) = -w(ish,k)
            do l = 1,neq
               dv(k) = dv(k)-v(ish,k,l)*dx(l)
            enddo
         enddo
         do k = 1,neq
            dx(k) = dv(k)
            xishk = xmod(ish,k)
            if (ish.eq.1.or.xishk.eq.0.d0) goto 20
            xx = dx(k)/xmod(ish,k)
            if (k.eq.5.and.abs(xmod(ish,5)).lt.1.d-15*Lmax) xx = 0.d0
            if (abs(xx).le.abs(corm(k)).or.(numeric.eq.1.and.t(ish).lt.
     &           tmaxioHe)) goto 20
            corm(k) = xx
            ncm(k) = ish
 20         xmod(ish,k) = xmod(ish,k)+alpha1(k)*dx(k)
         enddo
      enddo
      
      xmod(1,2) = max(xmod(1,2),log(1.d-99))
      t(1) = exp(min(xmod(1,4),23.d0))
      r(1) = exp(xmod(1,2))
      do ish = 2,nmod
         t(ish) = exp(min(xmod(ish,4),23.d0))
         r(ish) = exp(min(xmod(ish,2),50.d0))
      enddo
      if (icorr.eq.'t') then
         do ish = neff,2,-1
c..  check if pressure gradient is a decreasing function of mass 
c..  except in the outer layers
            if (abs(corm(2)).lt.eps(2).or.(numeric.eq.1.and.t(ish).gt.
     &           tmaxioHe)) then
               dp = p(ish)-p(ish-1)
c            accel(ish) =(xmod(ish,1)-xa(ish,1)-psi(ish)*(xmod(ish,1)-
c     &           xmod(ish-1,1)))/dtn
               ag = accel(ish)+g*m(ish)/(r(ish)*r(ish))
               if (dp*ag.ge.0.d0) then
                  corm(2) = 9.99d0
                  ncm(2) = ish
                  exit
               endif
            endif
         enddo
      endif
     
      if (iter.eq.1) write (nout,1000)
      if (nsconv.le.1) write (nout,1001) iter,(alpha1(i),ncm(i),corm(i),
     &     i = 1,neq),(novlim(1,i),i = 3,4),nsconv,scr
      if (nsconv.eq.2) write (nout,1002) iter,(alpha1(i),ncm(i),corm(i),
     &     i = 1,neq),((novlim(kl,i),i = 3,4),kl = 1,2),max(nsconv,1),
     &     scr
      if (nsconv.eq.3) write (nout,1003) iter,(alpha1(i),ncm(i),corm(i),
     &     i = 1,neq),((novlim(kl,i),i = 3,4),kl = 1,3),max(nsconv,1),
     &     scr
      if (nsconv.ge.4) write (nout,1004) iter,(alpha1(i),ncm(i),corm(i),
     &     i = 1,neq),((novlim(kl,i),i = 3,4),kl = 1,4),max(nsconv,1),
     &     scr
      if (no.eq.1) then
         if (iter.eq.1) write (90,1000)
         if (nsconv.le.1) write (90,1001) iter,(alpha1(i),ncm(i),
     &        corm(i),i = 1,neq),(novlim(1,i),i = 3,4),nsconv,scr
         if (nsconv.eq.2) write (90,1002) iter,(alpha1(i),ncm(i),
     &        corm(i),i = 1,neq),((novlim(kl,i),i = 3,4),kl = 1,2),
     &        max(nsconv,1),scr
         if (nsconv.eq.3) write (90,1003) iter,(alpha1(i),ncm(i),
     &        corm(i),i = 1,neq),((novlim(kl,i),i = 3,4),kl = 1,3),
     &        max(nsconv,1),scr
         if (nsconv.ge.4) write (90,1004) iter,(alpha1(i),ncm(i),
     &        corm(i),i = 1,neq),((novlim(kl,i),i = 3,4),kl = 1,4),
     &        max(nsconv,1),scr
      endif
      
c.... NEW
      do i = 1,5
         corm(i) = abs(corm(i))
      enddo

      if (corm(1).gt.epsmax(1).or.corm(2).gt.epsmax(2).or.corm(3).gt.
     &     epsmax(3).or.corm(4).gt.epsmax(4).or.corm(5).gt.epsmax(5))
     &     then
         error = 6
         return 1
      endif
***   determine correction threshold in case nucreg = 3
      corrmax = abs(corm(2)).gt.eps(2)*1.d2.or.abs(corm(3)).gt.
     &     eps(3)*10.d0.or.abs(corm(4)).gt.eps(4)*10.d0.or.
     &     abs(corm(5)).gt.eps(5)*1.d5.or.abs(corm(1)).gt.eps(1)
      
*___________________________________________________________________
***   calculation of the associated nucleosynthesis and diffusion,
***   if required at each Newton-Raphson iteration for the structure
***   case : nucreg = 3
*-------------------------------------------------------------------
      
      if (nucreg.eq.3.and.iter.ge.(itmin-1).and..not.corrmax) then
         klcore = 0
         klpulse = 0
         klenv = 0
         if ((microdiffus.or.rotation.or.diffusconv.or.difftacho.or.
     &        thermohaline).and.nmixd.eq.0) then
c..   in case nmixd > 0, diffusion treated independently in netdiff
            if (idiffnuc.eq.2) then
               if (nucreg.ne.0) call chedif (error)
               if (error.gt.0) return 1
               call diffusion (1,error)
               if (error.gt.0) return 1
            else
               call diffusion (1,error)
               if (error.gt.0) return 1
               if (nucreg.ne.0) call chedif (error)
               if (error.gt.0) return 1
            endif
            if (idiffnuc.eq.3) then
               write (nout,*) 'Processing diffusive mixing'
               if (no.eq.1) write (90,*) 'Processing diffusive mixing'
               call diffusion (1,error)
               if (error.gt.0) return 1
            endif
         else
            if (nmixd.gt.0.and.nsconv.gt.0) then
               if (novlim(1,3).eq.1) klcore = 1
               do kl = 1,nsconv
                  if ((tau(novlim(kl,4)).lt.1.d2.or.t(novlim(kl,4))
     &                 .lt.1.d6).and.mr(novlim(kl,4)).gt.0.7d0.and.
     &                 klenv.eq.0) klenv = kl
                  if (t(novlim(kl,3)).gt.2.d8.and.novlim(kl,3).gt.1
     &                 .and.nphase.eq.5.and.klpulse.eq.0) klpulse = kl
               enddo
               if (nsconv.eq.1.and.klcore.eq.1) then
                  klenv = 0
                  klpulse = 0
               endif
               if (nmixd.eq.1) klenv = 0
               if (nmixd.eq.2) klcore = 0
               if (nmixd.le.3) klpulse = 0
            endif
            call chedif (error)
            if (error.gt.0) return 1
         endif
      endif
      if (iter.gt.iitmax) then
         error = 99
         return 1
      endif
      if (iacc) call nwracc (corm,iter,alpha1)
c..   if lmix (convective boundaries have moved and mixing has been performed) 
c..   do another iteration
      if (iter.lt.itmin.or.lmix) goto 10

      imerk = 0
      do i = 1,neq
         if (abs(corm(i)).gt.eps(i)) imerk = imerk+1
      enddo
      if (imerk.gt.0) goto 10

 1000 format (1x,'#','|',3x,'Variable u(m)',4x,'|',3x,'Variable r(m)',
     &     4x,'|',2x,'Variable lnf(m)',3x,'|',2x,'Variable lnT(m)',3x,
     &     '|',3x,'Variable l(m)',4x,'|',1x,'first conv. zones',/)
 1001 format (i2,5('|',0pf5.3,1x,i4,1x,1pe9.2),'|',(1x,i4,'-',i4),1x,
     &     i2,a1)
 1002 format (i2,5('|',0pf5.3,1x,i4,1x,1pe9.2),'|',2(1x,i4,'-',i4),1x,
     &     i2,a1)
 1003 format (i2,5('|',0pf5.3,1x,i4,1x,1pe9.2),'|',3(1x,i4,'-',i4),1x,
     &     i2,a1)
 1004 format (i2,5('|',0pf5.3,1x,i4,1x,1pe9.2),'|',4(1x,i4,'-',i4),1x,
     &     i2,a1)

      return
      end


************************************************************************

      SUBROUTINE opa_co (tk,rok,muiinvk,x,ksh,kap1,kapr1,kapt1,error)

************************************************************************
*     compute opacity and derivatives in cool surface layers           *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none
      
      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.data'
      include 'evolcom.ion'
      include 'evolcom.nuc'

      logical logicNR,debug

      integer i,j,k,ij,ksh
      integer order(neqchimmax),indx(neqchimmax+natomsmax)
      integer error,neqo,n1,info
      integer imol,im,ip,jcur,jmol
      integer ir1,nr1,ir2,nr2,ip1
      integer iterelec

C for test
c      double precision rapCO
C 
C for print
c      double precision nnh,nno2,nno,nch,nn2,no,nn,nc
C
      double precision x(nsp),tk,rok,muiinvk
      double precision xxh,xxhe,xxz
      double precision theta,rhompinv
      double precision xmol(neqchimmax+natomsmax),
     &     xmol0(neqchimmax+natomsmax)
      double precision logKp(neqchimmax),logK(neqchimmax)
      double precision kap1,kapr1,kapt1
      double precision logKt,tempnbr
      double precision kdiss,temp,tempbis,nel0im
      double precision errf,tolf,errx,tolx
      double precision fjac(neqchimmax+natomsmax,neqchimmax+natomsmax),
     &     fvec(neqchimmax+natomsmax),p(neqchimmax+natomsmax)
      double precision kdisso(neqchimmax)
      double precision nelec,zscale,lntk,logt,betaev,delta,nelecb
      double precision y,a,b,pel
      double precision nh,nh2,nh2o,noh,nco,ncn,nc2
      double precision Kion(nioniz),ntot(nioniz),lnQ(2),xnsp(nioniz)
      double precision expq
      double precision nhtot
      double precision tt,tm1,t2,tm2,t3,t4,tm4,t5,tm5,t6,tm6,t7,t10
      double precision tm7,t9,t05
      double precision rm1,r05,r025,rm075
      double precision temp1,temp2,temp3,temp4,temp5,temp6
      double precision temp1r,temp3r
      double precision temp1t,temp2t,temp3t,temp4t,temp5t,temp6t
      double precision temp3a
      double precision temp7,temp7r,temp7t
      double precision temp8,temp8r,temp9r,temp9t
      double precision temp10,temp10t,temp11
      double precision temp13r
      double precision temp11t,temp12,temp12r,temp12a,temp12t,temp13a
      double precision temp13b,temp13,temp9
      double precision num2,den2inv,dden2t
      double precision num3,num3inv,den3inv,dnum3r,dden3r,dnum3t,dden3t
      double precision den4inv,dden4t
      double precision num5,den5inv,dden5t
      double precision num6,den6inv,dden6t
      double precision num7,den7inv,dnum7t,dden7t
      double precision num9,den9inv,dden9t
      double precision num10,den10inv,dden10t
      double precision num11,den11inv,dden11t
      double precision expdeninv,expfrac
      double precision tkm1,thetatm1

c..  initial isotopic number abondances (number densities)
      rhompinv = rok/mprot
      do j = 1,nspecies
         xmol(j) = 0.d0
         xmol0(j) = 0.d0
      enddo  
      xmol(iih) = (x(ih1)/anuc(ih1)+x(ih2)/anuc(ih2))*rhompinv
      xmol(iic) = (x(ic12)/anuc(ic12)+x(ic13)/anuc(ic13)+
     &     x(ic14)/anuc(ic14))*rhompinv
      xmol(iin) = (x(in13)/anuc(in13)+x(in14)/anuc(in14)+
     &     x(in15)/anuc(in15))*rhompinv
      xmol(iio) = (x(io15)/anuc(io15)+x(io16)/anuc(io16)+
     &     x(io17)/anuc(io17)+x(io18)/anuc(io18))*rhompinv

C --> For test
c      RapCO = 0.99d0
c      xmol(iic) = RapCO*xmol(iio)
c      xmol(iin) = xmol(iic)*0.4d0
C <--

      xmol0(iih) = xmol(iih)
      xmol0(iic) = xmol(iic)
      xmol0(iin) = xmol(iin)
      xmol0(iio) = xmol(iio)

C need for opacity calculation
      xxh = x(ih1)+x(ih2)
      xxhe = x(ihe3)+x(ihe4)
      xxz = 1-xxh-xxhe
           
*-------------------------------------------------------------
***   calculation of the dissociation constants K'(T) and K(T)
*-------------------------------------------------------------

      theta = 5040.d0/tk
      logkT = log10(boltz*tk)

      do i = 1,natoms
         logKp(i) = 0.d0
      enddo
      do i = 1,neqchim
         logKp(i+natoms) = 0.d0
         do j = 1,5
            logKp(i+natoms) = logKp(i+natoms)+akap(i,j)*theta**(j-1)
         enddo
         if (ireac(i,1).eq.1) then
            tempnbr = ireac(i,3)-1.d0
         elseif (ireac(i,1).eq.2) then 
            tempnbr = ireac(i,3)+ireac(i,5)-1.d0
         endif
         logK(i+natoms) = logKp(i+natoms)-tempnbr*logkT
      enddo

c..   sort logK by increasing order of magnitude
      do i = 1,neqchim
         ij = natoms+1
         do j = natoms+1,natoms+neqchim
            if (logK(j).lt.logK(i+natoms)) ij = ij+1
         enddo
         order(ij-natoms) = i+natoms
         kdisso(i+natoms) = 10.d0**logK(i+natoms)
      enddo     

*----------------------------------------------------------------
***   compute molecular and atomic number densities : xmol (cm-3)
*----------------------------------------------------------------

C find a first solution
      do j = 1,neqchim
         jcur = order(j)
         kdiss = kdisso(jcur)
         jmol = jcur-natoms

C Molecules type X2
         if (ireac(jmol,1).eq.1) then 
            if (ireac(jmol,3).ne.2) stop 'Pb in opa_co 1'
            im = ireac(jmol,2)
            imol = ireac(jmol,4)
            temp = dsqrt(1.d0+8.d0*xmol(im)/kdiss)-1.d0
            xmol(im) = 0.25d0*kdiss*temp
            xmol(imol) = 0.25d0*xmol(im)*temp
C Molecules type XY
         elseif (ireac(jmol,1).eq.2) then 
c     find the most abundant species 
            if (xmol(ireac(jmol,2)).gt.xmol(ireac(jmol,4))) then
               ip = ireac(jmol,2)
               im = ireac(jmol,4)
            else
               ip = ireac(jmol,4)
               im = ireac(jmol,2)
            endif
            imol = ireac(jmol,6)
            if (ireac(jmol,3)+ireac(jmol,5).eq.2) then !molecule type XY
               temp = kdiss+xmol(ip)-xmol(im)
               nel0im = xmol(im)
               xmol(im) = 0.5d0*temp*(dsqrt(1.d0+4.d0*kdiss*xmol(im)/
     &              temp**2)-1.d0)
               xmol(imol) = abs(nel0im-xmol(im))
               xmol(ip) = xmol(ip)-xmol(imol)
            else                !molecule type H20
               if (ireac(jmol,3).eq.2) then
                  ip = ireac(jmol,2)
                  im = ireac(jmol,4)
               else
                  ip = ireac(jmol,4)
                  im = ireac(jmol,2)
               endif
                  imol = ireac(jmol,6)
               if (xmol(im).lt.xmol(ip)) then ! (nH>>n0)
                  xmol(im) = kdiss/(kdiss+xmol(ip)**2)*xmol(im)
                  xmol(imol) = xmol(ip)**2*xmol(im)/kdiss
                  xmol(ip) = xmol(ip)-2.d0*xmol(imol)
               else             ! (nH<<nO)
                  xmol(imol) = xmol(ip)
                  xmol(ip) = 0.d0
                  xmol(im) = xmol(im)-xmol(imol)
               endif
            endif
         else
            print *,jmol,j,nreac
            print *,ireac(jmol,1),ireac(jmol,2),ireac(jmol,3),
     &           ireac(jmol,4),ireac(jmol,5),ireac(jmol,6)
            stop 'Pb in opa_co 2'
         endif
      enddo

      logicNR = .true.
      debug = .false.
      neqo = neqchimmax+natomsmax
      n1 = natoms+neqchim

c.. Start newtom-raphson iterations

      if (logicNR) then
         tolf = 1.d-12
         tolx = 1.d-3
         do k = 1,100

C initialisation
            fvec(1:neqchimmax) = 0.d0
            fjac(1:neqchimmax,1:neqchimmax) = 0.d0


C --> Compute vectors and jacobien matrix
c     atoms
            do i = 1,natoms
               fvec(i) = xmol(i)-xmol0(i)
               fjac(i,i) = 1.d0
            enddo
            do i = 1,neqchim
               ir1 = ireac(i,2)
               nr1 = ireac(i,3)
c     monoatomic molecules
               if (ireac(i,1).eq.1) then
                  ip1 = ireac(i,4)

                  fvec(ir1) = fvec(ir1)+nr1*xmol(ip1)
                  fvec(ip1) = xmol(ir1)**nr1-kdisso(ip1)*xmol(ip1)

                  fjac(ir1,ip1) = fjac(ir1,ip1)+nr1
                  if (fjac(ip1,ip1).ne.0.d0) stop 'PB in usrfun 1'
                  fjac(ip1,ip1) = -kdisso(ip1)
                  fjac(ip1,ir1) = nr1*xmol(ir1)**(nr1-1)
c     diatomic molecules
               elseif (ireac(i,1).eq.2) then
                  ir1 = ireac(i,2)
                  nr1 = ireac(i,3)
                  ir2 = ireac(i,4)
                  nr2 = ireac(i,5)
                  ip1 = ireac(i,6)

                  fvec(ir1) = fvec(ir1)+nr1*xmol(ip1)
                  fvec(ir2) = fvec(ir2)+nr2*xmol(ip1)
                  fvec(ip1) = xmol(ir1)**nr1*xmol(ir2)**nr2-kdisso(ip1)*
     &                 xmol(ip1)

                  fjac(ir1,ip1) = fjac(ir1,ip1)+nr1
                  fjac(ir2,ip1) = fjac(ir2,ip1)+nr2
                  if (fjac(ip1,ip1).ne.0.d0) then
                     print *,ip1,fjac(ip1,ip1)
                     stop 'PB in usrfun 2'
                  endif
                  fjac(ip1,ip1) = -kdisso(ip1)
                  if (nr1.eq.1) then 
                     fjac(ip1,ir1) = xmol(ir2)**nr2
                  else
                     fjac(ip1,ir1) = nr1*xmol(ir1)**(nr1-1)*
     &                    xmol(ir2)**nr2
                  endif
                  if (nr2.eq.1) then 
                     fjac(ip1,ir2) = xmol(ir1)**nr1
                  else
                     fjac(ip1,ir2) = nr2*xmol(ir2)**(nr2-1)*
     &                    xmol(ir1)**nr1
                  endif
c     no triatomic molecules
               else
                  stop 'PB in usrfun 3'
               endif
            enddo

            if (debug) then
               write (nout,*) '/****  VECT  ****/'
               write (nout,100) (fvec(j),j=1,neqchim+natoms)
               write (nout,*) '/****  JACO  ****/'
               do i = 1,neqchim+natoms
                  write (nout,100) (fjac(i,j),j=1,neqchim+natoms)
 100              format (30(1x,1pe8.1))
               enddo
            endif

            errf = 0.d0
            do i = 1,nspecies     !Check function convergence.
               errf = errf+abs(fvec(i))
            enddo
            if(errf.le.tolf) then
               print *,'1',k
               goto 10
            endif
            do i = 1,natoms+neqchim !Right-hand side of linear equations
               p(i) = -fvec(i) 
            enddo

c --> Solve linear equations
c.. leqs
c            call leqs(fjac,p,n1,neqo)
c.. lapack
            call dgetrf (neqo,natoms+neqchim,fjac,neqo,indx,info)
            call dgetrs ('N',n1,1,fjac,neqo,indx,p,n1,info)

            errx = 0.d0           !Check root convergence.
            do i = 1,nspecies     !Update solution.
               if (xmol(i).ne.0.d0) errx = errx+abs(p(i)/xmol(i))
               if (.not.(p(i).le.0.d0).and..not.(p(i).ge.0.d0)) then
                  print *,i,p(i),ksh,tk,rok
                  stop 'PB !!!!!!!'
               endif
               xmol(i) = xmol(i)+p(i)
            enddo
            if(errx.le.tolx) then
               goto 10
            endif
         enddo

      endif
C verification de la solution converg�e :
 10   if (k.ge.100) then
         print *,'Too many iteration in NR procedure: T =',tk,
     &        ' rho =',rok
         error = 84
         return
c         stop "N-R failed"
      endif
      do i = 1,nmole+4
         if (xmol(i).lt.0.d0) then
            if (xmol(i).gt.-1.d-5) then
               xmol(i) = 0.d0
            else
               write (nout,200) imol,ksh
 200           format ('error in abundances determination for ',i2,
     &              ', Abund < 0 at shell :',i4)
               write (nout,*) 'abundance :',xmol(i)
               error = 84
               return
c               stop "Abundance < 0 in opa_co"
            endif
         endif
      enddo

*-------------------------------
***  determine electron pressure
*-------------------------------

      nelec = 0.d0
      xnsp(1:nioniz) = 1.d-50
c...  compute abundances according to molecules formation
      do i = 2,nbz-1
         if (i.lt.6.or.i.gt.8) then
            do j = 1,nis
               if (znuc(j).eq.i) xnsp(i) = xnsp(i)+x(j)/anuc(j)
            enddo
            xnsp(i) = xnsp(i)*rhompinv
         endif
      enddo
      xnsp(1) = xmol(1) !H
      xnsp(6) = xmol(2) !C
      xnsp(7) = xmol(3) !N
      xnsp(8) = xmol(4) !O
c..   scale from Grevesse 95 for elements > Cl
      zscale = zkint*rhompinv/zsol
      do k = 1,7
         xnsp(k+17) = xspsol(nis-1+k)/aheavy(k)*zscale
      enddo

c...  compute partition functions
      lntk = log(tk)
      logT = log10(5040.d0/tk)
      betaev = -1.d0/(tk*8.61173d-5)

      ij = 0
      do j = 1,nioniz
c..   Irwin table only considers F,Ca,Ti,Cr,Mn,Fe,Ni (i.e. j=9,19..24)
         if (j.eq.9.or.j.ge.19) then
            ij = ij+1
            lnQ(1) = aIrwin(ij,1)
            lnQ(2) = aIrwin(ij,7)
            do i = 2,6
               lnQ(1) = lnQ(1)+aIrwin(ij,i)*lntk**(i-1)
               lnQ(2) = lnQ(2)+aIrwin(ij,i+6)*lntk**(i-1)
            enddo
            expQ = exp(lnQ(1)-lnQ(2))

         else
            lnQ(1) = aSauval(j,1)
            lnQ(2) = aSauval(j,6)
            do i = 2,5
               lnQ(1) = lnQ(1)+aSauval(j,i)*logT**(i-1)
               lnQ(2) = lnQ(2)+aSauval(j,i+5)*logT**(i-1)
            enddo
            expQ = 10.d0**(lnQ(1)-lnQ(2))

         endif
         ntot(j) = xnsp(j)
         Kion(j) = 2.d0*expQ*saha*tk**1.5d0*exp(Eion(j)*betaev)
      enddo

      delta = 1.d-3
      nelecb = rhompinv*muiinvk*10**(-1.d1+tk*1.d-3)
      iterelec = 0
      y = 2.d0*delta
      do while (abs(y).gt.delta)
         iterelec = iterelec+1
         nelec = nelecb
         y = nelec
         a = 1
         do j = 1,nioniz
            temp = 1.d0/(nelec+Kion(j))
            tempbis = Kion(j)*ntot(j)*temp
            y = y-tempbis
            a = a+tempbis*temp
         enddo
         b = y-a*nelec
         nelecb = -b/a
c         if (iterelec.ge.1d2) print *,iterelec,y,delta,abs(y).gt.delta
         if (iterelec.eq.1d2) goto 20
      enddo

 20   pel = nelec*((1.d0-nelec/(muiinvk*rhompinv)))*boltz*tk

      nh = xmol(iih)
c      nc = xmol(iic)
c      nn = xmol(iin)
c      no = xmol(iio)
      nh2 = xmol(iih2)
      nh2o = xmol(iih2o)
      noh = xmol(iioh)
      nco = xmol(iico)
      ncn = xmol(iicn)
      nc2 = xmol(iic2)
c      nn2 = xmol(iin2)
c      nch = xmol(iich)
c      nno = xmol(iino)
c      no2 = xmol(iio2)
c      nnh = xmol(iinh)

Ccc... Verifications
c      write (111,11) ksh,tk,logK(natoms+1),logK(natoms+2),logK(natoms+3)
c     &     ,logK(natoms+4),logK(natoms+5),logK(natoms+6),logK(natoms+7),
c     &     nh,nc,nn,no,nh2,nh2o,noh,nco,ncn,nc2,rok,nn2,nch,nno,no2,nnh,
c     &     xmolhplus,nelec

c 11   format (1x,i4,40(1x,2pe11.4))


*------------------
***   compute kappa
*------------------

      nhtot = nh+2*nh2+2*nh2o+noh

      tt = tk*1.d-4
      tm1 = 1.d0/tt
      t2 = tt**2
      tm2 = tm1**2
      t3 = tt*t2
      t4 = t2**2
      tm4 = tm2**2
      t5 = t2*t3
      tm5 = tm1*tm4
      t6 = t3**2
      tm6 = tm1*tm5
      t7 = t3*t4
      tm7 = tm1*tm6
      t9 = t7*t2
      t10 = tt*t9
      t05 = dsqrt(tt)

      rm1 = 1.d0/rok
      r05 = dsqrt(rok)
      r025 = dsqrt(r05)
      rm075 = rm1**0.75d0

      kap1 = 0.d0
      kapr1 = 0.d0
      kapt1 = 0.d0


c Keeley 70
c electrons
      temp1 = 5.4d-13*rm1*tm1
      temp1t = -temp1*tm1
      temp1r = -temp1*rm1

      num2 = t05
      den2inv = 1.d0/(2.d6*tm4+2.1d0*t6)
      dden2t = -8.d6*tm5+12.6d0*t5
      temp2 = num2*den2inv
      temp2t= temp2*(0.5d0*tm1-dden2t*den2inv)

      temp3a = 4.d-3*r025+2.d-4*t4
      num3 = temp3a
      num3inv = 1.d0/num3
      den3inv = 1.d0/(4.5d0*t6*temp3a+r025*t3)
      dnum3r = 1.d-3*rm075
      dden3r = rm075*t3*(0.25d0+4.5d-3*t3)
      dnum3t = 8.d-4*t3
      dden3t = t2*(r025*(3.d0+0.108d0*t3)+9.d-3*t7)
      temp3 = (1.d0-2.d0*nh2/nhtot)*num3*den3inv
      temp3r = temp3*(dnum3r*num3inv-dden3r*den3inv)
      temp3t = temp3*(dnum3t*num3inv-dden3t*den3inv)

      den4inv = 1.d0/(1.4d3*tt+t6)
      dden4t = 1.4d3+6.d0*t5
      temp4 = den4inv
      temp4t = -temp4*dden4t*den4inv

      num5 = 1.5d0
      den5inv = 1.d0/(1.d6+0.1d0*t6)
      dden5t = 0.6d0*t5
      temp5 = num5*den5inv
      temp5t = -temp5*dden5t*den5inv

      num6 = t05
      den6inv = 1.d0/(20.d0*tt+5.d0*t4+t5)
      dden6t = 20.d0+20.d0*t3+5.d0*t4
      temp6 = num6*den6inv
      temp6t = temp6*(0.5d0*tm1-dden6t*den6inv)

      kap1 = pel*(temp1+xxh*(temp2+temp3)+xxhe*(temp4+temp5)+xxz*temp6)
      kapr1 = pel*(temp1r+xxh*temp3r)
      kapt1 = pel*(temp1t+xxh*(temp2t+temp3t)+xxhe*(temp4t+temp5t)+
     &     xxz*temp6t)

c H et H2
      num7 = 5.55d-27*t4
      den7inv = 1.d0/(1.d0+10.d0*t6+3.42d-5*tm6)
      dnum7t = 2.22d-26*t3
      dden7t = 60.d0*t5-2.052d-4*tm7
      temp7 = (nh+nh2)*rm1*num7*den7inv
      temp7r = -temp7*rm1
      temp7t = temp7*(dnum7t/num7-dden7t*den7inv)

c CO
      temp8 = 2.75d-26*nco*rm1
      temp8r = -temp8*rm1

c OH
      num9 = 1.4d-21*t6
      den9inv = 1.d0/(0.1d0+t6)
      dden9t = 6.d0*t5
      temp9 = noh*rm1*num9*den9inv
      temp9r = -temp9*rm1
      temp9t = temp9*(6.d0*tm1-dden9t*den9inv)

c H20
c Attention: changement du facteur temp11 suivant Marigo (2002)
      num10 = 2.6d-27
      den10inv = 1.d0/(4.23d-4+t4)
      dden10t = 4.d0*t3
      temp10 = num10*den10inv
      temp10t = -temp10*dden10t*den10inv

      expdeninv = 1.d0/(tt+0.37d0)
      expfrac = 3.2553d0*expdeninv
      num11 = 9.72d-21*exp(-expfrac)
      den11inv = 1.d0/(1.d0+3.78d3*t10)
      dden11t = 3.78d4*t9
      temp11 = num11*den11inv
      temp11t = temp11*(expfrac*expdeninv-dden11t*den11inv)

      temp12a = nh2o*rm1
      temp12 = temp12a*(temp10+temp11)
      temp12r = -temp12*rm1
      temp12t = temp12a*(temp10t+temp11t)

      kap1 = kap1+temp7+temp8+temp9+temp12
      kapr1 = kapr1+temp7r+temp8r+temp9r+temp12r
      kapt1 = 1.d-4*(kapt1+temp7t+temp9t+temp12t)

c Scalo & Ulrich 75 (CN) and  Querci et al 71 (C2)
      tkm1 = 1.d0/tk
      theta = 5040.d0*tkm1
      thetatm1 = theta*tkm1
      temp13a = -19.212d0+2.2479d0*theta-2.8069d0*theta**2+0.76d0*
     &     theta**3-0.078384d0*theta**4
      temp13b = -2.2479d0+2.d0*2.8069d0*theta-3.d0*
     &     0.76d0*theta**2+4.d0*0.078384d0*theta**3

      temp13 = 10.d0**temp13a*(ncn+nc2)*rm1
      temp13r = -temp13*rm1

      kap1 = kap1+temp13
      kapr1 = kapr1+temp13r
      kapt1 = kapt1+log(10.d0)*temp13*temp13b*thetatm1

      return
      end


************************************************************************

      SUBROUTINE oversh

************************************************************************
*   Calculate parametric overshooting : artificial extension of CZ     *
*   Also computes the diffusion coefficient associated with            *
*   convection in the overshooted layers (needed in case of rotation)  *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.conv'
      include 'evolcom.cons'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer imin,imax
      integer lmin,lmax
      integer kl,kl1,kdw,kup,nc
      integer i,j,k

      double precision dmconv,rext,rext1,rext2,dr,dxm
      double precision lambdac

      dimension lambdac(nsh)


      if (novopt.eq.3) novopt = 1
      if (novopt.eq.4) novopt = 2
      do kl = 1,nsconv
         imin = novlim(kl,1)
         imax = novlim(kl,2)
***   treat core overshoot only
         if (imin.gt.1.and.novopt.gt.2) goto 30

***   overshoot below CZ
         dmconv = mr(imax)-mr(imin)
         rext = adwn*hp(imin)
         if (t(imin).lt.1.d6) then
            i = imin-1
            goto 10
         endif
         do i = imin,1,-1
            if ((r(imin)-r(i)).gt.rext.or.(mr(imin)-mr(i)).gt.
     &           dmconv*adwn)  goto 10
         enddo
 10      if (i.eq.1) i = 0
         lmin = i+1

***   overshoot top CZ
         rext1 = aup*hp(imax)
         if (t(imax).lt.1.d6) then
            i = imax+1
            goto 20
         endif
         rext2 = aup*r(imax)
         do i = imax,nmod
            dr = r(i)-r(imax)
            dxm = mr(i)-mr(imax)
            if (novopt.eq.1.and.(dr.gt.rext1.or.
     &           dxm.gt.dmconv*aup)) goto 20
            if (novopt.eq.2.and.(dr.gt.min(rext1,rext2))) goto 20
         enddo
 20      lmax = i-1

         if (imin.gt.lmin) then
            do j = lmin,imin-1
               sconv(j) = sconv(imin)
               tconv(j) = tconv(imin)
               fconv(j) = lum(j)/(pim4*r(j)*r(j))
               lambdac(j) = hp(j)*alphac
               Dconv(j) = sconv(j)*lambdac(j)/3.d0
               crz(j) = 1
            enddo
         endif
         if (lmax.gt.imax) then
            do j = imax+1,lmax
               sconv(j) = sconv(imax)
               tconv(j) = tconv(imax)
               fconv(j) = lum(j)/(pim4*r(j)*r(j))
               lambdac(j) = hp(j)*alphac
               Dconv(j) = sconv(j)*lambdac(j)/3.d0
               crz(j) = 1
            enddo
         endif
         novlim(kl,3) = lmin
         novlim(kl,4) = lmax
         novlim(kl,7) = lmin
         novlim(kl,8) = lmax
      enddo


***   redefine novlim if zones overlap
 30   nc = nsconv
      do kl = 1,nsconv-1
         kl1 = kl+1
         if (novlim(kl,4).ge.novlim(kl1,3)) then
            kdw = kl
            if (novlim(kl1,3).lt.novlim(kl,3)) kdw = kl1
            kup = kl1
            if (novlim(kl,4).gt.novlim(kl1,4)) kup = kl
            novlim(kl,3) = novlim(kdw,3)
            novlim(kl,4) = novlim(kup,4)
            novlim(kl,7) = novlim(kdw,7)
            novlim(kl,8) = novlim(kup,8)
            do k = kl1+1,nc
               novlim(k-1,3) = novlim(k,3)
               novlim(k-1,4) = novlim(k,4)
               novlim(k-1,7) = novlim(k,3)
               novlim(k-1,8) = novlim(k,4)
            enddo
            novlim(nc,3) = 0
            novlim(nc,4) = 0
            novlim(nc,7) = 0
            novlim(nc,8) = 0
            nc = nc-1
            if (nc.eq.nsconv-1) goto 40
         endif
      enddo
 40   nsconv = nc

      return
      end


************************************************************************

      DOUBLE PRECISION FUNCTION PCONV (a,V,B)

************************************************************************
* Search for the convective solution in case V <= 0.1                  *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      double precision a,V,B
      double precision a2,a3,V2,p,q,p1,f,g,dp,tol

      V2 = V*V
      a2 = a*a
      a3 = a*a2
      p = 1.d0/a
      p1 = 1.d0/a3
      tol = max((p-p1)*1.d-8,1.d-12)
      pconv = p
 10   q = (p+p1)*0.5d0
      f =  B*(1.d0-p*a3)*(1.d0-p*a3)**2*p*(1.d0+p**2*a)+
     &     (1+p)*(p**2*a2-1.d0)*((1.d0-p*a**3)*a2*p*(p-a)-
     &     (p**2*a2-1)**2/V2)
      g =  B*(1.d0-q*a3)*(1.d0-q*a3)**2*q*(1.d0+q**2*a)+
     &     (1+q)*(q**2*a2-1.d0)*((1.d0-q*a**3)*a2*q*(q-a)-
     &     (q**2*a2-1)**2/V2)
      if ((f.gt.0.d0.and.g.lt.0.d0).or.(f.lt.0.d0.and.g.gt.0.d0))
     &     then
         p1 = q
      else
         p = q
      endif
      pconv = p
      dp = abs(p-p1)/p
      if (dp.lt.tol) goto 30
      goto 10

 30   return
      end
      SUBROUTINE prvar

************************************************************************
* Write results in evolhr, evolvar* and evolch* evolution files        *
* Modifs CC ondes (23/11/07)                                           *
* $LastChangedDate:: 2018-01-22 10:23:29 +0000 (Mon, 22 Jan 2018)     $ *
* $Author:: amard                                                    $ *
* $Rev:: 112                                                         $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.eng'
      include 'evolcom.eos'
      include 'evolcom.grad'
      include 'evolcom.lum'
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
      include 'evolcom.ion'

      include 'evolcom.igw'

      character omegaconv*1,omegaconv0*1

C modif NL
      integer nzoneth,itth,ibth
C... fin modif NL
      integer icore,ienv,ib,it,jtot,nnmax,nt,inuctop,jth,nnmaxth
      integer i,i1,j,k,l,kp1,kl,ii,ip
      integer ibC,imC,itC,ibNe,imNe,itNe
      integer ibH,imH,itH,ibHe,imHe,itHe,imax,nzone,idiffex,normx
      integer klenv,kbce,ktce,bce,k_he,kconv,pass,klpulse,klcore

      integer iromax,irohalf,ndtenv

      double precision xmm(nmaxconv,3),xrr(nmaxconv,3),xmmth(nmaxconv,3)
     $     ,xrrth(nmaxconv,3)
      double precision shlimH,shlimHe,shlimC,shlimNe,shlimO
      double precision tmax,xmtmax,xmHb,xtHb,xroHb,xrHb,xmHm,xmHt,xtHt
     &     ,xroHt,xrHt,xmHeb,xtHeb,xroHeb,xrHeb,xmHem,xmHet,xtHet
     &     ,xroHet,xrHet,emaxHe,emaxH,xmb,xtb,xrobenv,xrb,xmt,xtt
     &     ,xrotenv,xrt,xmenvb,xtenvb,xroenvb,xrenvb,xmCt,xtCt,xroCt
     &     ,xrCt,xmCb,xtCb,xroCb,xrCb,emaxC,enumaxC,etamC,xmNeb,xrNeb
     &     ,xtNeb,xroNeb,xmNet,xrNet,xtNet,xroNet,xmNem,emaxNe,enumaxNe
     &     ,etamNe,enumaxHe,etamHe,xmCm,xmxs,xrxs,summ,xnuc
      double precision xmxsth,xrxsth,xmthb,xmtht
     $     ,xrthb,xrtht
      double precision lpph,rhomax,irradmax,lnu,leloc,leconv,lcno
      double precision xmoment,momenvel,xmoment_vieux,varmom
      double precision  epp(nsh),eCNO(nsh),e3a(nsh),eC12(nsh),
     &     eNe20(nsh),eO16(nsh),eppcno(nsh)
      double precision depot,dekin,deint,eenutri,eephot,lnuc,etot,edyn,
     &     totequil,vlnuc,eenuc,xmrad
      double precision eepot,veepot,eekin,veekin,eeint,veeint,eevisc
      double precision y(nsp),v(nreac)
      double precision diffst,diffcst
      double precision tautime,Ereac,taureac,rok
      double precision rossby,tconv_midCE
      double precision jrad,jconv,jtrue,jobs

      double precision romax,rohalf
      double precision omegamax,Nmax,brunt,BVmax,Nu2max,somegam
      double precision omcore,rmid(nsh),rmax,rmin,rcav

C Modif LA Coupling timescale (01/15)
      integer ndt,nshmax
      double precision Dshear,Fcirc,Fshear,tautot,taucirc,taushear,tauJ

C Modif LA convective core tc (03/16)
      double precision tg_core,tc_max_core,rc_max_core,tc_m_core,
     &     rc_m_core,rc_r_core,tc_r_core,tc_hp_core,rc_hp_core,tc_core,
     &     rc_core,omega_core,mcc,rcc

C Modif CC energie ondes (11/10/07)
      double precision londes
      double precision hpbce,rbce,hptce,factdc,tc,rc,tc_hp,rc_hp,rce
      double precision tc_r,rc_r,mce,tc_m,rc_m,tc_max,rc_max,tg,secd
      double precision nt2(nsh),nu2(nsh),bruntv2(nsh),bruntv(nsh),sum
      double precision dnu,dnuech,nu_max,errr,t_tot,sum_bce,t_bce,sum_he
      double precision r_he,t_he,sum_dp,dp
      double precision omegacore,omegajrad
      double precision Bcrit,fstar,Bequi
      double precision tauc

      character*2 cseq
      logical diffzc,eq_thermique

      common / rossby_number / tconv_midCE
      common /omega_convection/ omegaconv,omegaconv0
      common /overshoot/ klcore,klenv,klpulse
      common /momreel/ xmoment,momenvel,xmoment_vieux
      common /turbulence/ diffst,diffcst,idiffex,diffzc,eq_thermique
      common /nucleaire/ tautime(nsh,nreac),Ereac(nsh,nreac),taureac(nsh
     $     ,nreac)
      common /envel/ ndtenv
      common /brunt/ brunt(nsh),bruntv2,bruntv
      common /boreas_Mloss/ Bcrit,fstar,Bequi
      common /tauc_hp/ tauc


*____________________________________
***   print banner of evolution files
*------------------------------------

      if (no.eq.1.and.maxmod.eq.1.and.imodpr.le.0) then
         write (10,100) code_version
         write (11,110) code_version
         write (12,120) code_version
         write (13,130) code_version
         write (14,140) code_version
         write (15,150) code_version !H
         write (16,160) code_version !He
         write (17,162) code_version !C
         write (18,164) code_version !Ne
         write (27,170) code_version
         write (28,180) code_version

         write (29,260) code_version
         write (30,260) code_version

         write (21,190) code_version,'Central '
         write (22,200) code_version,'Central '
         write (123,210) code_version,'Central '
         write (124,220) code_version,'Central '
         write (31,190) code_version,'Surface '
         write (32,200) code_version,'Surface '
         write (33,210) code_version,'Surface '
         write (34,220) code_version,'Surface '
         write (41,230) code_version
         write (42,240) code_version
         write (43,250) code_version
      endif
      cseq = '  '
      if (model.eq.modeli+1) cseq = ' @'

*__________________________________________________________________
***   restore the accretion equilibrium abundance of H1, H2 and He3
*------------------------------------------------------------------

      if (dmaccr.gt.0.d0.and.accphase.le.1) then
         do k = iaccbot,iacctop
            xsp(k,ih2) = vxsp(k,ih2)
            xsp(k,ihe3) = vxsp(k,ihe3)
            normx = 1
            summ = 0.d0
            do l = 1,nsp
               xsp(k,l) = max(xsp(k,l),1.d-50)
               if (xsp(k,l).gt.xsp(k,normx)) normx = l
               summ = summ+xsp(k,l)
            enddo
            xsp(k,normx) = xsp(k,normx)+(1.d0-summ)
         enddo
      endif


c..   evolhr
      if (model.eq.modeli+1) then
         if (soll.le.99.9999d0.and.solr.le.9.9999d0.and.r(nmod)/rsun.le.
     &        9.9999d0) then
            write (10,1000) model,nphase,soll,solr,r(nmod)/rsun,teffy,
     &           roeff,geff,dms,totm,dtn*seci,time*seci,min(iter,99),
     &           ireset,nmod,cputime,cseq
         else
            if (r(nmod)/rsun.lt.9999.9d0) then
               write (10,1100) model,nphase,soll,solr,r(nmod)/rsun,
     &              teffy,roeff,geff,dms,totm,dtn*seci,time*seci,
     &              min(iter,99),ireset,nmod,cputime,cseq
            else
!               write (10,1150) model,nphase,min(soll,9999999.d0),
               write (10,1150) model,nphase,soll,
     &              min(9999.d0,solr),r(nmod)/rsun,teffy,roeff,geff,
     &              dms,totm,dtn*seci,time*seci,min(iter,99),ireset,
     &              nmod,cputime,cseq
            endif
         endif
      else
         if (soll.le.99.9999d0.and.solr.le.9.9999d0.and.r(nmod)/rsun.le.
     &        9.9999d0) then
            write (10,1200) model,nphase,soll,solr,r(nmod)/rsun,teffy,
     &           roeff,geff,dms,totm,dtn*seci,time*seci,min(iter,99),
     &           ireset,nmod,cputime
         else
            if (r(nmod)/rsun.lt.9999.9d0) then
!               write (10,1300) model,nphase,min(soll,9999999.d0),solr,
               write (10,1300) model,nphase,soll,solr,
     &              r(nmod)/rsun,teffy,roeff,geff,dms,totm,dtn*seci,
     &              time*seci,min(iter,99),ireset,nmod,cputime
            else
!               write (10,1350) model,nphase,min(soll,9999999.d0),
               write (10,1350) model,nphase,soll,
     &              min(9999.d0,solr),r(nmod)/rsun,teffy,roeff,geff,
     &              dms,totm,dtn*seci,time*seci,min(iter,99),ireset,
     &              nmod,cputime
            endif
         endif
      endif

*__________________________________________
***   determination of the max. temperature
*------------------------------------------

      imax = 1
      tmax = t(1)
      xmtmax = 0.d0
      do i = 1,nmod1
         if (t(i).gt.tmax) then
            tmax = t(i)
            imax = i
         endif
      enddo
      xmtmax = m(imax)/msun
      rhomax = ro(imax)


c. Modif LA 15/01/2013
*______________________________________________
***   determination of the half-max. density
*----------------------------------------------

      romax = 0.d0
      do i = 1,nmod1
         if (ro(i).gt.romax) then
            romax = ro(i)
            iromax = i
         endif
      enddo
      ii = 0
      rohalf = romax
      do while (rohalf.gt.(romax/2.d0))
         ii = ii+1
         rohalf = ro(ii)
         irohalf = ii
      enddo
      if (irohalf.ge.ndtenv) irohalf = ndtenv

c. Fin Modif LA

*____________________________________________________________________
***   calculation of global energetic characteristics of the new star
*--------------------------------------------------------------------

      eepot = 0.d0
      veepot = 0.d0
      eekin = 0.d0
      veekin = 0.d0
      eeint = 0.d0
      veeint = 0.d0
      eevisc = 0.d0
      lnuc = 0.d0
C Modif CC energie ondes (11/10/07)
C*      londes = 0.d0
      vlnuc = 0.d0
      eenutri = 0.d0
      do k = 2,nmod1
         eepot = eepot-dm(k)*(m(k+1)/r(k+1)+m(k)/r(k))
         veepot = veepot-dm(k)*(m(k+1)/vr(k+1)+m(k)/vr(k))
         if (hydro) then
            eekin = eekin+(u(k)+u(k+1))**2*dm(k)
            veekin = veekin+(vu(k)+vu(k+1))**2*dm(k)
         endif
         eeint = eeint+e(k)*dm(k)
         veeint = veeint+ve(k)*dm(k)
         eevisc = eevisc+evisc(k)*dm(k)
         lnuc = lnuc+enucl(k)*dm(k)
         vlnuc = vlnuc+venucl(k)*dm(k)
         eenutri = eenutri+enupla(k)*dm(k)
C Modif CC energie ondes (11/10/07)
C*         londes = londes+eondes(k)*dm(k)
      enddo
      eepot = eepot*g*0.5d0
      veepot = veepot*g*0.5d0
      eekin = eekin*0.125d0
      veekin = veekin*0.125d0
      eenuc = 0.5d0*(lnuc+vlnuc)*dtn
      eenutri = eenutri*dtn
      eephot = 0.5d0*(lum(nmod)+vlum(nmod))*dtn
      depot = eepot-veepot
      dekin = eekin-veekin
      deint = eeint-veeint
      edyn = depot+dekin+deint
      etot = edyn+eephot
      totequil = abs(1.d0-eenuc/etot)
      if (abs(totequil).gt.99.99d0) totequil = 99.9999d0
c..   dEtot = d(Ekin+Epot+Eint+Ephot) = dEnuc
      write (90,1450) depot,dekin,deint,eephot,eenutri,eenuc,totequil

      tkh = m(nmod)*m(nmod)*g/abs(r(nmod)*lum(nmod))

*______________________________________________________________________
***   calculation of the momentum of inertia, of the gravothermal and
***   nuclear luminosities associated to H, He, C, Ne, O and Si-burning
***   and neutron irradiation (for He-burning zones)
*----------------------------------------------------------------------

      vlpp = lpp
      vlh = lh
      vlhe = lhe
      vlc = lc
      vlne = vlne
      vlo = vlo
      vlsi = vlsi

c     Determination of the nuclear energy contributions

c... modif NL 6/11/2008
	do k=1,nmod
           rok = ro(k)
            call vit (t(k),ro(k),mueinv(k),k,v,0)
            call nucprod (v,y,t(k),epp(k),eCNO(k),e3a(k),eC12(k),
     &           eNe20(k),eO16(k))
***   2 PROT  ( 0 OOOOO, 0 OOOOO)  1 DEUT
            if (v(2).ne.0.d0) then
               tautime(k,ippg)=(2.d0)/(v(2)*xsp(k,ih1)*sec)
            else
               tautime(k,ippg) = 1.d99
            endif
            tautime(k,ippg) = min(tautime(k,ippg),1.d20)
            Ereac(k,ippg)=xsp(k,ih1)**2*v(2)*qi(ippg)*avn/(2.d0)
            taureac(k,ippg)=avn*rok*xsp(k,ih1)/(tautime(k,ippg))
***   1 DEUT  ( 1 PROT , 0 OOOOO)  1 HE  3
            if (v(3).ne.0.d0) then
               tautime(k,ipdg)=2.d0/(v(3)*xsp(k,ih1)*sec)
            else
               tautime(k,ipdg) = 1.d99
            endif
            tautime(k,ipdg) = min(tautime(k,ipdg),1.d20)
            Ereac(k,ipdg)=xsp(k,ih1)*xsp(k,ih2)*v(3)*qi(ipdg)*avn/(2.d0)
            taureac(k,ipdg)=avn*rok*xsp(k,ih2)/(2*tautime(k,ipdg))
***   2 HE  3 ( 0 OOOOO, 2 PROT )  1 HE  4
            if (v(8).ne.0.d0) then
               tautime(k,i2he3)=3*2.d0/(v(8)*xsp(k,ihe3)*sec)
            else
               tautime(k,i2he3) = 1.d99
            endif
            tautime(k,i2he3) = min(tautime(k,i2he3),1.d20)
            Ereac(k,i2he3)=xsp(k,ihe3)**2*v(8)*qi(i2he3)*avn/(9*2.d0)
            taureac(k,i2he3)=avn*rok*xsp(k,ihe3)/(3*tautime(k,i2he3))
***   1 HE  4 ( 1 HE  3, 0 OOOOO)  1 BE  7
            if (v(10).ne.0.d0) then
               tautime(k,ihe3ag)=4.d0/(v(10)*xsp(k,ihe4)*sec)
            else
               tautime(k,ihe3ag) = 1.d99
            endif
             tautime(k,ihe3ag) = min(tautime(k,ihe3ag),1.d20)
            Ereac(k,ihe3ag)=xsp(k,ihe3)*xsp(k,ihe4)*v(10)*qi(ihe3ag)
     &                    *avn/(12.d0)
            taureac(k,i2he3)=avn*rok*xsp(k,ihe3)/(3*tautime(k,ihe3ag))
***   1 BE  7 ( 0 betap, 0 nutri)  1 LI  7
            if (v(ibe7beta).ne.0.d0) then
               tautime(k,ibe7beta)=1.d0/(v(ibe7beta)*sec)
            else
               tautime(k,ibe7beta) = 1.d99
            endif
            tautime(k,ibe7beta) = min(tautime(k,ibe7beta),1.d20)
            Ereac(k,ibe7beta)=xsp(k,ibe7)*v(ibe7beta)*qi(ibe7beta)
     &      *avn/(7.d0)
            taureac(k,ibe7beta)=avn*rok*xsp(k,ibe7)/(7*tautime(k
     $           ,ibe7beta))
***   1 LI  7 ( 1 PROT , 0 OOOOO)  2 HE  4
            if (v(14).ne.0.d0) then
               tautime(k,ili7pa)=1.d0/(v(14)*xsp(k,ih1)*sec)
            else
               tautime(k,ili7pa) = 1.d99
            endif
            tautime(k,ili7pa) = min(tautime(k,ili7pa),1.d20)
            Ereac(k,ili7pa)=xsp(k,ili7)*xsp(k,ih1)*v(14)*qi(ili7pa)*avn
     &                     /(7.d0)
            taureac(k,ili7pa)=avn*rok*xsp(k,ili7)/(7*tautime(k,ili7pa))
***   1 BE  7 ( 1 PROT , 0 OOOOO)  1 B   8
            if (v(17).ne.0.d0) then
            tautime(k,ibe7pg)=1.d0/(v(17)*xsp(k,ih1)*sec)
            else
               tautime(k,ibe7pg) = 1.d99
            endif
            tautime(k,ibe7pg) = min(tautime(k,ibe7pg),1.d20)
            Ereac(k,ibe7pg)=xsp(k,ibe7)*xsp(k,ih1)*v(17)*qi(ibe7pg)*avn
     &                     /(7.d0)
            taureac(k,ibe7pg)=avn*rok*xsp(k,ibe7)/(7*tautime(k,ibe7pg))
***   1 B   8 ( 0 OOOOO, 0 OOOOO)  2 HE  4
            if (v(ib8beta).ne.0.d0) then
               tautime(k,ib8beta)=1.d0/(v(ib8beta)*sec)
            else
               tautime(k,ib8beta) = 1.d99
            endif
            tautime(k,ib8beta) = min(tautime(k,ib8beta),1.d20)
            Ereac(k,ib8beta)=xsp(k,ib8)*v(ib8beta)*qi(ib8beta)*avn
     $           /(8.d0)
            taureac(k,ib8beta)=avn*rok*xsp(k,ib8)/(8*tautime(k,ib8beta))
***   1 C  13 ( 1 PROT , 0 OOOOO)  1 N  14
            if (v(35).ne.0.d0) then
               tautime(k,ic13pg)=1.d0/(v(35)*xsp(k,ih1)*sec)
            else
               tautime(k,ic13pg) = 1.d99
            endif
            tautime(k,ic13pg) = min(tautime(k,ic13pg),1.d20)
            Ereac(k,ic13pg)=xsp(k,ic13)*xsp(k,ih1)*v(35)*qi(ic13pg)*avn
     &                    /(13.d0)
            taureac(k,ic13pg)=avn*rok*xsp(k,in14)/(13*tautime(k,ic13pg))
***   1 N  14 ( 1 PROT , 0 OOOOO)  1 O  15
            if (v(48).ne.0.d0) then
               tautime(k,in14pg)=1.d0/(v(48)*xsp(k,ih1)*sec)
            else
               tautime(k,in14pg) = 1.d99
            endif
            tautime(k,in14pg) = min(tautime(k,in14pg),1.d20)
            Ereac(k,in14pg)=xsp(k,in14)*xsp(k,ih1)*v(48)*qi(in14pg)*avn
     &                    /(14.d0)
            taureac(k,in14pg)=avn*rok*xsp(k,in14)/(14*tautime(k,in14pg))
***   1 C  12 ( 1 PROT , 0 OOOOO)  1 N  13
            if (v(30).ne.0.d0) then
               tautime(k,icpg)=1.d0/(v(30)*xsp(k,ih1)*sec)
            else
               tautime(k,icpg) = 1.d99
            endif
            tautime(k,icpg) = min(tautime(k,icpg),1.d20)
            Ereac(k,icpg)=xsp(k,ic12)*xsp(k,ih1)*v(30)*qi(icpg)*avn
     &                    /(12.d0)
            taureac(k,icpg)=avn*rok*xsp(k,ic12)/(12*tautime(k,icpg))
c...fin de modif
	enddo

      inuctop = nmod1
      shlimH = 0.d0
      shlimHe = 0.d0
      shlimC = 0.d0
      shlimNe = 0.d0
      shlimO = 0.d0
      do k = 1,nmod1
         if (t(k).ge.tnucmin) then
            do l = 1,nsp
               y(l) = ysp(k,l)
            enddo
            call vit (t(k),ro(k),mueinv(k),k,v,0)
            call nucprod (v,y,t(k),epp(k),eCNO(k),e3a(k),eC12(k),
     &           eNe20(k),eO16(k))
            eppcno(k) = epp(k) + eCNO(k)
            shlimH = max(shlimH,eppcno(k))
            shlimHe = max(shlimHe,e3a(k))
            shlimC = max(shlimC,eC12(k))
            shlimNe = max(shlimNe,abs(eNe20(k)))
            shlimO = max(shlimO,eO16(k))
          else
            inuctop = k-1
            goto 111
         endif
      enddo

 111  irradmax = -1.d0
      do k = 1,inuctop
         if (irradmax.lt.exposure(k)) then
            irradmax = exposure(k)
            xmrad = m(k)/msun
         endif
      enddo

c..   evolvar1
      if (eta(1).lt.9999.9d0) then
         write (11,1400) model,t(1),tmax,xmtmax,ro(1),rhomax,p(1),
     &        beta(1),eta(1),degpe(1),enupla(1),enucl(1),egrav(1),cseq
      else
         write (11,1410) model,t(1),tmax,xmtmax,ro(1),rhomax,p(1),
     &        beta(1),eta(1),degpe(1),enupla(1),enucl(1),egrav(1),cseq
      endif


*___________________________________________________________
***   determination of the nuclear burning region boundaries
***   and properties
*-----------------------------------------------------------

      shlimH = min(shlimH,shlim)
      shlimHe = min(shlimHe,shlim)
      if (nphase.lt.2) shlimHe = 1.d10
      if (nphase.ge.5) then
         shlimC = min(shlimC,shlim)
         if (nphase.ge.6) then
            shlimNe = min(abs(shlimNe),shlim)
            shlimO = min(shlimO,shlim)
         endif
      else
         shlimC = shlim
         shlimNe = shlim
         shlimO = shlim
      endif
      shlimNe = max(shlimNe,1.d-10)

c Burning regions for PP, CNO, 3a, C12 and Ne20

      ibH = 0
      itH = 0
      imH = 0
      xmHm = 0.d0
      emaxH = 0.d0
c      enumaxH = 0.d0
c      etamH = 0.d0
      xmHb = 0.d0
      xrHb = 0.d0
      xtHb = 0.d0
      xroHb = 0.d0
      xmHt = 0.d0
      xrHt = 0.d0
      xtHt = 0.d0
      xroHt = 0.d0

      ibHe = 0
      itHe = 0
      imHe = 0
      xmHem = 0.d0
      emaxHe = 0.d0
      enumaxHe = 0.d0
      etamHe = 0.d0
      xmHeb = 0.d0
      xrHeb = 0.d0
      xtHeb = 0.d0
      xroHeb = 0.d0
      xmHet = 0.d0
      xrHet = 0.d0
      xtHet = 0.d0
      xroHet = 0.d0

      ibC = 0
      itC = 0
      imC = 0
      xmCm = 0.d0
      emaxC = 0.d0
      enumaxC = 0.d0
      etamC = 0.d0
      xmCb = 0.d0
      xrCb = 0.d0
      xtCb = 0.d0
      xroCb = 0.d0
      xmCt = 0.d0
      xrCt = 0.d0
      xtCt = 0.d0
      xroCt = 0.d0

      ibNe = 0
      itNe = 0
      imNe = 0
      xmNem = 0.d0
      emaxNe = 0.d0
      enumaxNe = 0.d0
      etamNe = 0.d0
      xmNeb = 0.d0
      xrNeb = 0.d0
      xtNeb = 0.d0
      xroNeb = 0.d0
      xmNet = 0.d0
      xrNet = 0.d0
      xtNet = 0.d0
      xroNet = 0.d0
      xmNet = 0.d0
      xrNet = 0.d0
      xtNet = 0.d0
      xroNet = 0.d0

      xnuc = 1.5d0

      do k = 1,inuctop-1
c..   location of maximum energy production rate
         if (eppcno(k).ge.shlimH.and.eppcno(k).gt.emaxH) then
            imH = k
            emaxH = eppcno(imH)
         endif
         if (e3a(k).ge.shlimHe.and.e3a(k).gt.emaxHe) then
            imHe = k
            emaxHe = e3a(imHe)
         endif
         if (nphase.ge.5) then
            if (eC12(k).ge.shlimC.and.eC12(k).gt.emaxC.and.eC12(k).gt.
     &           xnuc*abs(eNe20(k))) then
               imC = k
               emaxC = eC12(imC)
            endif
            if (abs(eNe20(k)).ge.shlimNe.and.abs(eNe20(k)).gt.emaxNe
     &           .and.((nphase.le.6.and.abs(eNe20(k)).gt.xnuc*eC12(k))
     &           .or.(abs(eNe20(k)).gt.xnuc*eO16(k).and.nphase.gt.6)))
     &           then
               imNe = k
               emaxNe = abs(eNe20(imNe))
            endif
         endif
c..      top of burning shells
         kp1 = k+1
         if (eppcno(k).ge.shlimH.and.eppcno(kp1).lt.shlimH) itH = k
         if (e3a(k).ge.shlimHe.and.e3a(kp1).lt.shlimHe) itHe = k
         if (eC12(k).ge.shlimC.and.eC12(kp1).lt.shlimC) itC = k
         if (abs(eNe20(k)).ge.shlimNe.and.abs(eNe20(k)).gt.xnuc*eC12(k))
     &        itNe = k
      enddo

      itH = max(itH,imH)
      itHe = max(itHe,imHe)
      itC = max(itC,imC)
      itNe = max(itNe,imNe)

      do k = 1,inuctop-1
        kp1 = k+1
c..     base of burning shells
        if (imH.gt.0) then
           if (eppcno(k).ge.shlimH.and.k.eq.1) then
              ibH = k
           else if (eppcno(k).lt.shlimH.and.eppcno(kp1).ge.shlimH.and.
     &             ibH.eq.0) then
              ibH = kp1
           endif
        endif
        if (e3a(k).ge.shlimHe.and.k.eq.1) then
           ibHe = k
        else if (e3a(k).lt.shlimHe.and.e3a(kp1).ge.shlimHe.and.
     &          ibHe.eq.0) then
           ibHe = kp1
        endif
        if (itC.gt.0.and.ibC.eq.0.and.eC12(k).ge.shlimC.and.
     &       eC12(k).ge.xnuc*abs(eNe20(k)).and.k.gt.itNe) then
           ibC = k
        endif
        if (itNe.gt.0.and.ibNe.eq.0.and.abs(eNe20(k)).ge.shlimNe.and.
     &       ((nphase.le.6.and.abs(eNe20(k)).gt.xnuc*eC12(k)).or.
     &       (abs(eNe20(k)).ge.xnuc*eO16(k).and.nphase.gt.6))) then
           ibNe = k
        endif
      enddo

      if (itHe.ne.0.and.ibH.ne.0) itHe = min(max(ibH,imHe),itHe)
c      if (itC.ne.0) itC = min(ibHe,itC)
      if (itNe.ne.0) itNe = min(ibC,itNe)

      if (ibH.eq.0.and.itH.ne.0) stop 'prvar : error HBS boundaries'
      if (ibHe.eq.0.and.itHe.ne.0.or.imHe.gt.itHe) then
         print *,'HeBS : [',ibHe,':',imHe,':',itHe,']'
         stop 'prvar : error HeBS boundaries'
      endif
      if (ibC.eq.0.and.itC.ne.0.or.imC.gt.itC) then
         stop 'prvar : error CBS boundaries'
      endif
      if (ibNe.eq.0.and.itNe.ne.0.or.imNe.gt.itNe)
     &     stop 'prvar : error NeBS boundaries'

      r(1) = 0.d0
      lh = 0.d0
      lpp = 0.d0
      lcno = 0.d0
      lhe = 0.d0
      lc = 0.d0
      lne = 0.d0
      lo = 0.d0
      lsi = 0.d0
      lnu = 0.d0

c hydrogen burning shell
      if (ibH.gt.0) then
         xmHb = m(ibH)/msun
         xrHb = r(ibH)/r(nmod)
         xtHb = t(ibH)
         xroHb = ro(ibH)
         xmHt = m(itH)/msun
         xrHt = r(itH)/r(nmod)
         xtHt = t(itH)
         xroHt = ro(itH)
         xmHm = m(imH)/msun
         emaxH = enucl(imH)
c         enumaxH = enupla(imH)
c         etamH = eta(imH)
         do i = ibH,itH
            lh = lh+(enucl(i)+enupla(i))*dm(i)
            lcno = lcno+eCNO(i)*dm(i)
            lpp = lpp+epp(i)*dm(i)
         enddo
      endif

c helium burning shell
      if (ibHe.gt.0) then
         xmHeb = m(ibHe)/msun
         xrHeb = r(ibHe)/r(nmod)
         xtHeb = t(ibHe)
         xroHeb = ro(ibHe)
         xmHet = m(itHe)/msun
         xrHet = r(itHe)/r(nmod)
         xtHet = t(itHe)
         xroHet = ro(itHe)
         xmHem = m(imHe)/msun
         emaxHe = enucl(imHe)
         enumaxHe = enupla(imHe)
         etamHe = eta(imHe)
         do i = ibHe,itHe
            lhe = lhe+(enucl(i)+enupla(i))*dm(i)
         enddo
      endif

c carbon burning shell
      if (ibC.gt.0) then
         xmCb = m(ibC)/msun
         xrCb = r(ibC)/r(nmod)
         xtCb = t(ibC)
         xroCb = ro(ibC)
         xmCt = m(itC)/msun
         xrCt = r(itC)/r(nmod)
         xtCt = t(itC)
         xroCt = ro(itC)
         xmCm = m(imC)/msun
         emaxC = enucl(imC)
         enumaxC = enupla(imC)
         etamC = eta(imC)
         do i = ibC,itC
            lc = lc+(enucl(i)+enupla(i))*dm(i)
         enddo
      endif

c neon burning shell
      if (ibNe.gt.0) then
         xmNeb = m(ibNe)/msun
         xrNeb = r(ibNe)/r(nmod)
         xtNeb = t(ibNe)
         xroNeb = ro(ibNe)
         xmNet = m(itNe)/msun
         xrNet = r(itNe)/r(nmod)
         xtNet = t(itNe)
         xroNet = ro(itNe)
         xmNem = m(imNe)/msun
         emaxNe = enucl(imNe)
         enumaxNe = enupla(imNe)
         etamNe = eta(imNe)
         do i = ibNe,itNe
            lne = lne+(enucl(i)+enupla(i))*dm(i)
         enddo
      endif

      totgrav = 0.d0
      totnucl = 0.d0
c      lshr = 0.d0

      do i = 1,nmod
         totnucl = totnucl+enucl(i)*dm(i)
         totgrav = totgrav+egrav(i)*dm(i)
         lnu = lnu+enupla(i)*dm(i)
c         lshr = lshr+eshr(i)*dm(i)*thacc
      enddo

      leloc = 0.d0
      leconv = 0.d0
      if (urca) then
         do i = 1,inuctop
            leloc = leloc+eloc(i)*dm(i)
            leconv = leconv+egconv(i)*dm(i)
         enddo
         leloc = leloc/lsun
         leconv = leconv/lsun
         leloc = sign(max(abs(leloc),1.d-30),leloc)
         leconv = sign(max(abs(leconv),1.d-30),leconv)
      endif

      totnucl = totnucl/lsun
      totgrav = totgrav/lsun
      lh = lh/lsun
      lpp = lpp/lsun
      lhe = lhe/lsun
      lc = lc/lsun
      lne = lne/lsun
c      lo = lo/lsun
c      lsi = lsi/lsun
      lnu = lnu/lsun
      lh = sign(max(abs(lh),1.d-30),lh)
      lhe = sign(max(abs(lhe),1.d-30),lhe)
      lc = sign(max(abs(lc),1.d-30),lc)
      lne = sign(max(abs(lne),1.d-30),lne)
c      lo = sign(max(abs(lo),1.d-30),lo)
c      lsi = sign(max(abs(lsi),1.d-30),lsi)
c      lshr = lshr/lsun

c..   evolvar2
      write (12,1500) model,lh,lhe,lc,lne,lo,xmrad,lnu,totnucl,
     &     totgrav,irradmax,cseq
c      write (12,1500) model,lh,lhe,lc,lne,lo,leloc,lnu,totnucl,
c     &     totgrav,leconv,cseq

c..   evolvar5
      if (lcno.eq.0.d0.and.lpp.eq.0.d0) then
         lpph = 0.d0
      else
         lpph = max(log(lpp/lcno),-99.9998d0)
      endif
      write (15,1600) model,xmHb,xrHb,xtHb,xroHb,xmHt,xrHt,xtHt,xroHt,
     &     xmHm,emaxH,lpph,cseq
c..   evolvar6
      write (16,1700) model,xmHeb,xrHeb,xtHeb,xroHeb,xmHet,xrHet,xtHet,
     &     xroHet,xmHem,emaxHe,enumaxHe,etamHe,cseq
c..   evolvar7
      write (17,1700) model,xmCb,xrCb,xtCb,xroCb,xmCt,xrCt,xtCt,
     &     xroCt,xmCm,emaxC,enumaxC,etamC,cseq
c..   evolvar8
      write (18,1700) model,xmNeb,xrNeb,xtNeb,xroNeb,xmNet,xrNet,xtNet,
     &     xroNet,xmNem,emaxNe,enumaxNe,etamNe,cseq
c..   evolvar9
c      write (19,1700) model,xmOb,xrOb,xtOb,xroOb,xmOt,xrOt,xtOt,
c     &     xroOt,xmOm,emaxO,enumaxO,etamO,cseq
cc..   evolvar10
c      write (20,1700) model,xmSib,xrSib,xtSib,xroSib,xmSit,xrSit,xtSit,
c     &     xroSit,xmSim,emaxSi,enumaxSi,etamSi,cseq



*_______________________________________________________
***   determination of the convective regions properties
*-------------------------------------------------------

      r(1) = 0.d0
      xmb = 0.d0
      xtb = 0.d0
      xrobenv = 0.d0
      xrb = 0.d0
      xmt = 0.d0
      xtt = 0.d0
      xrotenv = 0.d0
      xrt = 0.d0
      xmenvb = 0.d0
      xtenvb = 0.d0
      xroenvb = 0.d0
      xrenvb = 0.d0
      icore = 0
      ienv = 0
      if (nsconv.eq.0) goto 80
      do kl = 1,nsconv
c         if (novlim(kl,3).gt.1.and.mr(novlim(kl,4)).gt.0.7d0.and.
         if (novlim(kl,3).gt.1.and.
     &        (tau(novlim(kl,4)).lt.1.d2.or.t(novlim(kl,4)).lt.1.d6))
     &        then
            ienv = kl
            goto 60
         endif
      enddo
 60   if (ienv.gt.0) then
         ib = novlim(ienv,3)
         xmenvb = m(ib)/msun
         xtenvb = log10(t(ib))
         xroenvb = log10(ro(ib))
         xrenvb = r(ib)/r(nmod)
      endif
      if (nsconv.gt.0.and.ienv.ne.1) then
         ib = novlim(1,3)
         it = novlim(1,4)
         icore = 1
         xmb = m(ib)/msun
         xtb = log10(t(ib))
         xrobenv = log10(ro(ib))
         xrb = r(ib)/r(nmod)
         xmt = m(it+1)/msun
         xtt = log10(t(it))
         xrotenv = log10(ro(it))
         xrt = r(it)/r(nmod)
      endif

c..   evolvar3
 80   write (13,1900) model,xmb,xrb,xtb,xrobenv,xmt,xrt,xtt,xrotenv,
     &     xmenvb,xrenvb,xtenvb,xroenvb,cseq

*_____________________________________________________
***   selection of the remaining convective zones
***   and sorting as a function of their relative mass
*-----------------------------------------------------

      jtot = 0
      do i = 1,nsconv
         if (i.ne.icore.and.i.ne.ienv.and.novlim(i,3).ne.0) then
            jtot = jtot+1
            xmm(jtot,1) = m(novlim(i,3))/msun
            xmm(jtot,2) = m(novlim(i,4)+1)/msun
            xmm(jtot,3) = xmm(jtot,2)-xmm(jtot,1)
            xrr(jtot,1) = r(novlim(i,3))/rsun
            xrr(jtot,2) = r(novlim(i,4)+1)/rsun
            xrr(jtot,3) = xrr(jtot,2)-xrr(jtot,1)
         endif
      enddo
      do i = 1,jtot
         do j = i+1,jtot
            if (xmm(j,3).gt.xmm(i,3)) then
               xmxs = xmm(j,1)
               xmm(j,1) = xmm(i,1)
               xmm(i,1) = xmxs
               xmxs = xmm(j,2)
               xmm(j,2) = xmm(i,2)
               xmm(i,2) = xmxs
               xmxs = xmm(j,3)
               xmm(j,3) = xmm(i,3)
               xmm(i,3) = xmxs
               xrxs = xrr(j,1)
               xrr(j,1) = xrr(i,1)
               xrr(i,1) = xrxs
               xrxs = xrr(j,2)
               xrr(j,2) = xrr(i,2)
               xrr(i,2) = xrxs
               xrxs = xrr(j,3)
               xrr(j,3) = xrr(i,3)
               xrr(i,3) = xrxs
            endif
         enddo
      enddo

*_____________________________________________________________
***   selection of the next four most massive convective zones
*-------------------------------------------------------------

      nnmax = 5
      if (jtot.lt.nnmax) nnmax = jtot
      do i = 1,nnmax
         do j = i,nnmax
            if (xmm(j,1).lt.xmm(i,1)) then
               xmxs = xmm(j,1)
               xmm(j,1) = xmm(i,1)
               xmm(i,1) = xmxs
               xmxs = xmm(j,2)
               xmm(j,2) = xmm(i,2)
               xmm(i,2) = xmxs
               xmxs = xmm(j,3)
               xmm(j,3) = xmm(i,3)
               xmm(i,3) = xmxs
               xrxs = xrr(j,1)
               xrr(j,1) = xrr(i,1)
               xrr(i,1) = xrxs
               xrxs = xrr(j,2)
               xrr(j,2) = xrr(i,2)
               xrr(i,2) = xrxs
               xrxs = xrr(j,3)
               xrr(j,3) = xrr(i,3)
               xrr(i,3) = xrxs
            endif
         enddo
      enddo
      do i = nnmax+1,nmaxconv
         xmm(i,1) = 0.d0
         xmm(i,2) = 0.d0
         xmm(i,3) = 0.d0
         xrr(i,1) = 0.d0
         xrr(i,2) = 0.d0
         xrr(i,3) = 0.d0
      enddo
      if (nsconv.le.0) then
         nzone = 1
      else
         nzone = nsconv
      endif

c... modif NL  19/09/08

         xmthb=0.d0
         xmtht=0.d0
	 xrthb=0.d0
	 xrtht=0.d0

         if (novlim(1,11).ge.1) then
            xmthb=m(novlim(1,11))/msun
            xmtht=m(novlim(1,12)+1)/msun
            xrthb=r(novlim(1,11))/rsun
            xrtht=r(novlim(1,12)+1)/rsun
         endif
	 write(*,*) xmthb,xmtht
C	 stop

c..   evolvar4
      write (14,2000) model,xmm(1,1),xmm(1,2),xmm(2,1),xmm(2,2),
     &     xmm(3,1),xmm(3,2),xmm(4,1),xmm(4,2),xmthb,xmtht,
     &     min(nzone,99),scr,cseq
c.. evolvar12

       write (28,2050) model,xrr(1,1),xrr(1,2),xrr(2,1),xrr(2,2),
     &     xrr(3,1),xrr(3,2),xrr(4,1),xrr(4,2),xrthb,xrtht,
     &     cseq


*____________________________________________________
***   storage of the accretion and rotation variables
*----------------------------------------------------

      k2rad = 0.d0
      k2conv = 0.d0
      jrad = 0.d0
      jconv = 0.d0
      angconvr = 0.d0
      angradr = 0.d0
      jobs = 0.d0
      if (dmaccr.eq.0.d0) then
         raccbot = 1.d0
         maccbot = 1.d0
         lacc = 0.d0
         facc(1:nmod) = 0.d0
      endif
      if (irotbin.eq.0) then
         Fenerg = 0.d0
         vsurf = 0.d0
         ur_Ss = 0.d0
      endif

***   Computing the radius of gyration of the radiative zone (k2rad)
***   and of the convective envelope (k2conv) according to Rucinski 1988 (AJ 95)

c$$$      if (ienv.gt.0.and.(omegaconv.ne.'s')) then
      if (ienv.gt.0) then
         if (novlim(ienv,3).gt.0) then
            print *, 'ienv',ienv,'novlim(ienv,3)',novlim(ienv,3)
     $           ,'novlim(1,4)',novlim(1,4) 
            do i = novlim(ienv,3),novlim(ienv,4)
               i1 = i-1
               k2conv = k2conv+dm(i1)*(r(i)**2+r(i1)**2+r(i)*r(i1))
     &              /3.d0
            enddo

            if (nsconv.eq.1) then
               do i = 2,novlim(ienv,3)-1
                  i1 = i-1
                  k2rad = k2rad+dm(i1)*(r(i)**2+r(i1)**2+r(i)*r(i1))
     &                 /3.d0
                  jrad = jrad+dm(i1)*(r(i)**2+r(i1)**2+r(i)*r(i1))
     &                 /3.d0*0.5d0*(omega(i1)+omega(i))
                  omegacore = omegacore+omega(i1)+omega(i)
               enddo
               omegajrad = jrad/k2rad
               omegacore = omegacore/(2*novlim(ienv,3)-1)
               print *, 'nsconv',nsconv,'omegajrad',omegajrad
     $              ,'omegacore',omegacore
            else if (nsconv.ge.2) then
               do i = novlim(1,4)+1,novlim(ienv,3)-1 ! Test TD pour prendre en compte la présence d'un coeur convectif 
                  i1 = i-1
                  k2rad = k2rad+dm(i1)*(r(i)**2+r(i1)**2+r(i)*r(i1))
     &                 /3.d0
                  jrad = jrad+dm(i1)*(r(i)**2+r(i1)**2+r(i)*r(i1))
     &                 /3.d0*0.5d0*(omega(i1)+omega(i))
                  omegacore = omegacore+omega(i1)+omega(i)
               enddo
               omegajrad = jrad/k2rad
               omegacore = omegacore/(2*(novlim(ienv,3)-novlim(1,4)+1)
     $              -1)
               print *, 'nsconv',nsconv,'omegajrad',omegajrad
     $              ,'omegacore',omegacore
            endif
         else
            do i = 2,nmod
               i1 = i-1
               k2conv = k2conv+dm(i1)*(r(i)**2+r(i1)**2+r(i)*r(i1))
     &              /3.d0
            enddo
            omegacore = omega_S
         endif
      else
         do i = 2,nmod
            i1 = i-1
            k2conv = k2conv+dm(i1)*(r(i)**2+r(i1)**2+r(i)*r(i1))
     &           /3.d0
         enddo
         omegacore = omega_S
         omegajrad = omega_S
      endif

      jrad = pw23*jrad
      jtrue = (pw23*k2conv*omega_S+jrad)/(totm*msun)
      jobs = (pw23*(k2conv+k2rad))*omega_S/(totm*msun)

      k2conv = k2conv*2.d0/(3.d0*totm*msun*(solr*rsun)**2)
      k2rad = k2rad*2.d0/(3.d0*totm*msun*(solr*rsun)**2)
      varmom = 0.d0


      if (irotbin.eq.1) then
         if (.not.diffzc.and.novlim(nsconv,3)-1.le.1.and.nphase.eq.1)
     &        then
            varmom = 1.d0
         else
            momenvel = momenvel
            xmoment = xmoment
            varmom = 1.d0-xmoment_vieux/xmoment
         endif
c     Rossby number for rotating models : Prot / tconv
c     Ro = 2*pi*R/(vsurf*tc)*8.055 with R in units of Rsun, vsurf in units of km/s and tc in units of days
         rossby = pim2/(omega_S*tconv_midCE)
         if (tconv_midCE.eq.0.d0) rossby = 0.d0
c..   evolvar11
      else
         Fenerg = 0.d0
         omega_S = 0.d0
         rossby = 0.d0
         vsurf = 0.d0
         xmoment = 0.d0
         dmoment1 = 0.d0
      endif
      write (27,2100) model,massrate,raccbot,maccbot,lacc/lsun,k2conv
     &     ,k2rad,rossby,omega_S,vsurf,xmoment,Fenerg,dmoment1,cseq



!!! Mean the rotation rate below 0.025 Rtot for asteroseismic comparison

      Nt2=0.d0
      Nu2=0.d0
      do i = 1,nmod
         bruntV2(i) = 0.d0
         bruntV(i) = 0.d0
      enddo

      do i=1,nmod-1
         grav(i) = gmr(i)/gmr(nmod)
         gs = gmr(nmod)

         Nt2(i) = gs*grav(i)*deltaKS(i)*(abad(i)-abla(i))/(hp(i))
         Nu2(i) = gs*grav(i)*phiKS(i)*abmu(i)/(hp(i))

         bruntV2(i) = Nt2(i)+Nu2(i)
          
      enddo
     
      nshmax = 1
      do i=1,nmod
         if (i.lt.novlim(klenv,3).and.(Nu2(i).gt.Nu2max)) then
            Nu2max = Nu2(i)
            nshmax = i
         endif
      enddo
      
      somegam = 0.d0
      do i=1,nshmax
        if (omega(i).gt.0.d0) somegam = somegam + omega(i)*dm(i)
      enddo
      if (nshmax.ne.1) then 
         omegamax = somegam/m(nshmax)
      else 
         omegamax = 0.d0
      endif

c$$$      rmid(1)=r(1)
c$$$      omcore = 0.d0
c$$$      do i=2,nmod
c$$$         i1 = i-1
c$$$         rmid(i) = dsqrt(pw13*(r(i)**2+r(i1)**2+r(i)*r(i1)))
c$$$      enddo
c$$$      BVmax = 0.d0
c$$$      do i=1,novlim(klenv,3)
c$$$         if (bruntV2(i).gt.BVmax) then
c$$$            BVmax = bruntV2(i)
c$$$         endif
c$$$      enddo
c$$$
c$$$      i=0
c$$$      rmax = 0.d0
c$$$      rmin = r(nmod)
c$$$      do while (r(i).lt.0.10*r(nmod))
c$$$         if (bruntV2(i).ge.0.1d0*BVmax) then
c$$$            if (r(i).gt.rmax) rmax = r(i)
c$$$            if (r(i).lt.rmax) rmin = r(i)
c$$$            ip = i+1
c$$$            omcore=omcore+(rmid(ip)-rmid(i))*omega(i)
c$$$         endif      
c$$$         i = i+1
c$$$      enddo
c$$$      rcav = rmax-rmin
c$$$      omegamax = omcore/rcav
c$$$         
c$$$!!!
c$$$c..   evolvar13
c$$$      Nmax = 1.d-50
c$$$      omegamax = 1.d-50
c$$$      do i=1,nmod
c$$$         if (i.lt.novlim(klenv,3).and.(brunt(i).gt.Nmax))
c$$$     &        Nmax = brunt(i)
c$$$         if (omega(i).gt.omegamax) omegamax = omega(i)
c$$$      enddo
c$$$      print *,'omegamax2=',omegamax
c$$$
c$$$
c$$$      Nt2=0.d0
c$$$      Nu2=0.d0
c$$$      do i = 1,nmod
c$$$         bruntV2(i) = 0.d0
c$$$         bruntV(i) = 0.d0
c$$$      enddo
c$$$
c$$$      do i=1,nmod-1
c$$$         grav(i) = gmr(i)/gmr(nmod)
c$$$         gs = gmr(nmod)
c$$$
c$$$         Nt2(i) = gs*grav(i)*deltaKS(i)*(abad(i)-abla(i))/(hp(i))
c$$$         Nu2(i) = gs*grav(i)*phiKS(i)*abmu(i)/(hp(i))
c$$$
c$$$         bruntV2(i) = Nt2(i)+Nu2(i)
c$$$      enddo
c$$$      somegam = 0
c$$$      do i=1,nmod
c$$$         if (i.lt.novlim(klenv,3).and.(Nu2(i).gt.Nu2max)) then
c$$$            Nu2max = Nu2(i)
c$$$            nshmax = i
c$$$         endif
c$$$      enddo
c$$$      do i=1,nshmax
c$$$        somegam = somegam + omega(i)*dm(i)
c$$$      enddo
c$$$      omegamax = somegam/m(nshmax)
c$$$      print *,'omegamax1=',omegamax


      

cc Coupling between envelope and radiative core Timescale

      if ((novlim(ienv,3).gt.1).and.(novlim(klenv,3).gt.1)) then
         ndt = novlim(klenv,3)-1
         print *,'ndt=',ndt
         Dshear = 0.d0
         if (omega(ndt)/=omega(ndt-1))
     &        Dshear = 0.5d0*(xnuvv(ndt)+xnuvv(ndt-1))
         jconv = k2conv*omega_S*2.d0/(3.d0*totm*msun*(solr*rsun)**2)
         fcirc = -0.2d0*ro(ndt)*r(ndt)**4*omega(ndt)*urs(ndt)*ur_Ss
         print *,omega(ndt),omega(ndt-1)
         fshear = -ro(ndt)*r(ndt)**4*Dshear*
     &        (omega(ndt)-omega(ndt-1))/(r(ndt)-r(ndt-1))
         tauJ = seci*(k2conv*jrad-k2rad*jconv)/(k2conv+k2rad)
         if (fcirc.ne.0.d0.or.fshear.ne.0.d0) tautot = tauJ/
     &        (fcirc+fshear)
         if (fcirc.ne.0.d0) then
            taucirc = tauJ/fcirc
            if (fshear.ne.0.d0) then
               taushear = tauJ/fshear
            else
               taushear = 0.d0
            endif
         else
            taucirc = 0.d0
            if (fshear.eq.0.d0) then
               tautot = 0.d0
               taushear = 0.d0
            else
               taushear = tauJ/fshear
            endif
         endif
      else
         ndt = nmod
         fcirc = 0.d0
         fshear = 0.d0
         tauJ = 0.d0
         tautot = 0.d0
         taushear = 0.d0
         taucirc = 0.d0
      endif


c      if (irotbin.eq.1) then
c         omegacore = omega(irohalf)
c      else
c         omegacore = 0.d0
c      endif

c.. evolvar13
      write (29,2300) model,omegajrad,jtrue,jobs,jrad,Bequi,omegamax,
     &     Nmax,cseq


c$$$c.. evolvar14
c$$$      write (30,2600) model,omegacore,jtrue,jobs,jrad,Bequi,omegamax,
c$$$     &     Nmax,tautot,taucirc,taushear,cseq

*______________________________________
***   storage of the abundance profiles
*--------------------------------------

c..   evolchc*
      write (21,2200) model,(xsp(1,idxspc(i)),i = 1,11),cseq
      write (22,2200) model,(xsp(1,idxspc(i)),i = 12,22),cseq
      write (123,2200) model,(xsp(1,idxspc(i)),i = 23,33),cseq
      write (124,2200) model,(xsp(1,idxspc(i)),i = 34,44),cseq
c      write (25,2200) model,(xsp(1,idxspci(i)),i = 1,9),cseq
c      write (26) model,time,(xsp(1,i),i = 1,nsp)

      if (nsconv.gt.0) then
         nt = novlim(nsconv,4)
      else
         nt = neff
      endif
      if (tau(nt).ge.1.d2) nt = neff
c..   evolchs*
      write (31,2200) model,(xsp(nt,idxsps(i)),i = 1,11),cseq
      write (32,2200) model,(xsp(nt,idxsps(i)),i = 12,22),cseq
      write (33,2200) model,(xsp(nt,idxsps(i)),i = 23,33),cseq
      write (34,2200) model,(xsp(nt,idxsps(i)),i = 34,nprint),cseq

*__________________________________________________________________________
***   Storage of L,R,X and Y with maximum precision for solar calibration
*--------------------------------------------------------------------------

      if (time*seci.eq.4.5700000000d+09.and.totm.eq.1.d0) then
         open(unit=115,file='input_calibration',status='unknown')
         write(115,2400)soll,solr,xsp(nt,idxsps(2))+xsp(nt,idxsps(3)),
     &   xsp(nt,idxsps(4))+xsp(nt,idxsps(5))
         close(115)
      endif

*________________________________________
***   storage of convective turnover-time
*----------------------------------------

!!! Initialisation
      tg = 0.d0
      tc_max = 0.d0
      rc_max= 0.d0
      tc_m = 0.d0
      rc_m = 0.d0
      rc_r = 0.d0
      tc_r = 0.d0
      tc_hp = 0.d0
      rc_hp = 0.d0
      tc = 0.d0
      rc = 0.d0

      tg_core = 0.d0
      tc_max_core = 0.d0
      rc_max_core= 0.d0
      tc_m_core = 0.d0
      rc_m_core = 0.d0
      rc_r_core = 0.d0
      tc_r_core = 0.d0
      tc_hp_core = 0.d0
      rc_hp_core = 0.d0
      tc_core = 0.d0
      rc_core = 0.d0
      omega_core = 0.d0
!!!
      if (icore.eq.1) then
         print *,'convective core !'

C     turnover time at Hp/2 below top of CC
            k = it
            do while (r(k).lt.r(it)-0.5d0*hp(it)
     &           .and.k.gt.1)
               k = k-1
            enddo
            if (k.eq.0) k = 1
            if (k.lt.ib.or.sconv(k).eq.0.d0) then
               tc = 0.d0
               rc = 0.d0
            else
               tc_core = alphac*hp(k)/sconv(k)
               rc_core = r(k)
            endif

c     turnover time at Hp below top of CC
            k = it
            do while (r(k).gt.r(it)-hp(it).and.k.lt.nmod)
               k = k-1
            enddo
            if (k.eq.0) k = 1
            if (k.lt.ib.or.sconv(k).eq.0.d0) then
               tc_hp_core = 0.d0
               rc_hp_core = 0.d0
            else
               tc_hp_core = alphac*hp(k)/sconv(k)
               rc_hp_core = r(k)
            endif

c     turnover time at 1/2R_cc
            rcc = r(it)
            k = it
            do while (r(k).gt.r(it)-0.5d0*rcc.and.
     &           k.gt.ib)
               k = k-1
            enddo
            if (k.eq.0) k = 1
!            if (k.eq.ib) k = ib + 1
            if (sconv(k).eq.0.d0) then
               tc_r_core = 0.d0
               rc_r_core = 0.d0
            else
               tc_r_core = alphac*hp(k)/sconv(k)
               rc_r_core = r(k)
            endif

c     turnover time at 1/2M_cc
            mcc =  m(it)
            k = it
            do while (m(k).gt.0.5*mcc.and.
     &           k.gt.ib)
               k = k-1
            enddo
            if (k.eq.0) k = 1
!            if (k.eq.ib) k = ib + 1
            if (sconv(k).eq.0.d0) then
               tc_m_core = 0.d0
               rc_m_core = 0.d0
            else
               tc_m_core = alphac*hp(k)/sconv(k)
               rc_m_core = r(k)
            endif

c     turnover time with maximal value (bottom and top layers are not taken into account)
            tc_max_core = -1.d0
            do k = ib+1,it-1
               if (sconv(k).gt.0.d0.and.
     &              tc_max_core.lt.alphac*hp(k)/sconv(k)) then
                  tc_max_core = alphac*hp(k)/sconv(k)
                  rc_max_core = r(k)
               endif
            enddo
            if (rc_max_core.lt.1.d-30) rc_max = 1.d-30
            if (tc_max_core.eq.-1.d0) then
               tc_max_core = 0.d0
               rc_max_core = 0.d0
            endif

c     integreted turmover time (bottom and top layers are not taken into account)
            tg_core = 0.d0
            do k = ib+1,it-1
               if (sconv(k).gt.0.d0) then
                  tg_core = tg_core+(r(k+1)-r(k))/sconv(k)
               endif
            enddo
            if (irotbin.eq.1) omega_core = omega(it-1)
         endif

!!!! Envelope -----------------------------------------------
      if (ienv.gt.0.or.novlim(klenv,3)-1.le.1) then
         print *,'convective envelope !'
!         if (novlim(ienv,3).ne.0) then
            kbce = novlim(ienv,3)
            hpbce = hp(kbce)

            if (0.05d0*hp(novlim(ienv,3)).lt.r(nmod)) then
               do while (r(kbce).lt.r(novlim(ienv,3))+
     &              0.05d0*hp(novlim(ienv,3)).and.kbce.lt.nmod)
                  kbce = kbce+1
               enddo
            endif
            rbce = r(kbce)

            ktce = novlim(ienv,4)
            hptce = hp(ktce)
            if (nphase.gt.2) then
               factdc = 0.05d0
            else
               factdc = 0.1d0
            endif
            do while (r(ktce).gt.r(novlim(ienv,4))-
     &           factdc*hp(novlim(ienv,4))
     &           .and.kbce.gt.1)
               ktce = ktce-1
            enddo

C     turnover time at Hp/2 over bottom of CE
            k = novlim(ienv,3)
            do while (r(k).lt.r(novlim(ienv,3))+0.5d0*hp(novlim(ienv,3))
     &           .and.k.lt.nmod)
               k = k+1
            enddo
            if (k.eq.0) k = 1
            if (k.ge.novlim(ienv,4).or.sconv(k).eq.0.d0) then
               tc = 0.d0
               rc = 0.d0
            else
               tc = alphac*hp(k)/sconv(k)
               rc = r(k)
            endif

c     turnover time at Hp over bottom of CE
            k = novlim(ienv,3)
            do while (r(k).lt.r(novlim(ienv,3))+hp(novlim(ienv,3))
     &           .and.k.lt.nmod)
               k = k+1
            enddo
           if (k.eq.0) then
               k = 1
            else if (k.ge.novlim(ienv,4).or.sconv(k).eq.0.d0) then
               tc_hp = 0.d0
               rc_Hp = 0.d0
            else
               tc_hp = alphac*hp(k)/sconv(k)
               rc_hp = r(k)
            endif

c     turnover time at 1/2R_ce
            rce = r(novlim(ienv,4))-r(novlim(ienv,3))
            k = novlim(ienv,3)
            do while (r(k).lt.r(novlim(ienv,3))+0.5d0*rce.and.
     &           k.lt.novlim(klenv,4))
               k = k+1
            enddo
            if (k.eq.0) k = 1
            if (k.eq.novlim(ienv,4)) k = novlim(klenv,4)-1
            if (sconv(k).eq.0.d0) then
               tc_r = 0.d0
               rc_r = 0.d0
            else
               tc_r = alphac*hp(k)/sconv(k)
               rc_r = r(k)
            endif

c     turnover time at 1/2M_ce
            mce =  m(novlim(ienv,4))-m(novlim(ienv,3))
            k = novlim(ienv,3)
            do while (m(k).lt.m(novlim(ienv,3))+0.5d0*mce.and.
     &           k.lt.novlim(klenv,4))
               k = k+1
            enddo
            if (k.eq.0) k = 1
            if (k.eq.novlim(ienv,4)) k = novlim(klenv,4)-1
            if (sconv(k).eq.0.d0) then
               tc_m = 0.d0
               rc_m = 0.d0
            else
               tc_m = alphac*hp(k)/sconv(k)
               rc_m = r(k)
            endif

c     turnover time with maximal value (bottom and top layers are not taken into account)
            tc_max = -1.d0
            do k = kbce,ktce
               if (sconv(k).gt.0.d0.and.
     &              tc_max.lt.alphac*hp(k)/sconv(k)) then
                  tc_max = alphac*hp(k)/sconv(k)
                  rc_max = r(k)
               endif
            enddo
            if (rc_max.lt.1.d-30) rc_max = 1.d-30
            if (tc_max.eq.-1.d0) then
               tc_max = 0.d0
               rc_max = 0.d0
            endif


c     integrated turmover time (bottom and top layers are not taken into account)
            tg = 0.d0
            do k = kbce,ktce
               if (sconv(k).gt.0.d0) tg = tg+(r(k+1)-r(k))/sconv(k)
            enddo

!            write (41,2500) model,tc/secd,rc/rsun,tc_hp/secd,
!     &           rc_hp/rsun,tc_r/secd,rc_r/rsun,tauc,rc_m/rsun,
!     &           tc_max/secd,rc_max/rsun,tg/secd,cseq
!         else
!            write (41,2500) model,0.d0,0.d0,0.d0,
!     &           0.d0,0.d0,0.d0,0.d0,0.d0,
!     &           0.d0,0.d0,0.d0,cseq
!         endif
         endif
c..   evoltc1 : convective turnover in the envelope
            write (41,2500) model,tc*seci,rc/rsun,tc_hp*seci,
     &           rc_hp/rsun,tc_r*seci,rc_r/rsun,tauc,rc_m/rsun,
     &           tc_max*seci,rc_max/rsun,tg*seci,cseq
c..   evoltc2 : convective turnover in the core
            write (42,2510) model,tc_core*seci,rc_core/rsun,
     &           tc_hp_core*seci,rc_hp_core/rsun,tc_r_core*seci,
     &           rc_r_core/rsun,tauc,rc_m_core/rsun,tc_max_core*seci,
     &           rc_max_core/rsun,tg_core*seci,omega_core,cseq
!

*______________________________________
***   storage of asteroseimic relations
*--------------------------------------

*--------
* fréquence de brunt Vaisalla
*--------

      Nt2=0.d0
      Nu2=0.d0
      do i = 1,nmod
         bruntV2(i) = 0.d0
         bruntV(i) = 0.d0
      enddo

      do i=1,nmod-1
         grav(i) = gmr(i)/gmr(nmod)
         gs = gmr(nmod)

         Nt2(i) = gs*grav(i)*deltaKS(i)*(abad(i)-abla(i))/(hp(i))
         Nu2(i) = gs*grav(i)*phiKS(i)*abmu(i)/(hp(i))

         bruntV2(i) = Nt2(i)+Nu2(i)
      enddo

      do i=2,nmod-1	
         if (bruntV2(i).gt.0.d0.and.crz(i).ne.-2) then
            bruntV(i) = dsqrt(bruntV2(i))

         elseif (bruntV2(i).le.0.d0.and.bruntV2(i-1).gt.0.d0
     &           .and.bruntV2(i+1).gt.0.d0
     &           .and.crz(i).ne.-2) then
            bruntV(i) = dsqrt(-bruntV2(i))
            bruntV2(i) = -bruntV2(i)
         else  		
            bruntV(i) = 0.d0
            bruntV2(i) = 0.d0
         endif
      enddo

c.... grande séparation avec relation asymptotique
        sum = 0.d0
      do k = 1,nmod-1
         sum = sum+(r(k+1)-r(k))/cs(k)
      enddo
      Dnu = 1.d0/(2.d0*sum)

c.... grande séparation avec les relations d'échelle
c     Dnu_sun = 134.9d-6 microHz
      Dnuech = 134.9d-6*(m(nmod)/msun)**0.5d0*(reff/rsun)**(-1.5d0)

c.... frenquency at which oscilation modes are maximal
c     nu_max_sun = 3150.d-6 microHz
      nu_max = 3150.d-6*(m(nmod)/msun)*(reff/rsun)**(-2.d0)*
     &  (teff/5780)**(-0.5d0)

c.... difference entre les relations asymptotiques et d'échelle
      errr = (Dnu-Dnuech)/Dnu

c....rayon acoustique totale
      T_tot=1/(2*Dnu)

c... rayon acoustique à la base de l'env.conv
      sum_BCE = 0.d0
      if (ienv.gt.0) then
         BCE = novlim(ienv,3)
         do k = 1,BCE
            sum_BCE=sum_BCE+(r(k+1)-r(k))/cs(k)
         enddo
      endif
      t_BCE = sum_BCE

c... Rayon acoustique a la base de la zone d'ionisation de l'He++
      sum_He = 0.d0
      do k= 1,nmod-1
	 if (gamma1(k).lt.1.55.and.t(k).le.1e5.and.
     &        gamma1(k).lt.gamma1(k+1)) then
            k_He=k
            goto 23
         endif
      enddo
 23   do k = 1,k_He
         sum_He=sum_He+(r(k+1)-r(k))/cs(k)
      enddo
	
      r_He=r(k_He)
      t_He=sum_He

c.... séparation en période entre les modes g pour l=1
      sum_Dp=0.d0
      Dp=0.d0
      kconv = 2
      do k = 2,nmod
         if (Dconv(k-1).eq.0.d0.and.Dconv(k).gt.0.d0) then
            kconv=k
            goto 25
         endif
      enddo
 25   if (ienv.gt.0) then
         if (kconv.eq.novlim(ienv,3)) then
            kconv=2
         endif
      endif
	
      pass = 0
	
 26   do k=kconv,nmod-1
         if ((2.d0*pi*nu_max)**(2.d0).lt.bruntV2(k).and.
     &        (2.d0*pi*nu_max)**(2.d0).lt.(2.d0*(cs(k)/r(k))**(2.d0)))
     &        sum_Dp=sum_Dp+((r(k+1)-r(k))*bruntV(k)*(r(k))**(-1.d0))
      enddo
      if (sum_Dp.lt.1.d-3.and.pass.eq.0.and.kconv.ne.2) then
         sum_Dp = 0.d0
         kconv = 2
         pass = 1
         goto 26
      endif	
	
      if (sum_Dp.gt.0.d0) then
         Dp=((2.0d0)**(0.5)*pi**(2))/(sum_Dp)
      else
         Dp=-1.0d0
      endif

c...evolas
      write (43,2600) model,Dnu,Dnuech,errr,T_tot,t_BCE,t_He,nu_max,Dp
     &     ,cseq


 100  format (' AGB network : version ',0pf4.2,/,1x,
     &     'General informations:',//,
     &     1x,'model',1x,'phase',4x,'L',6x,'Reff',5x,'R',6x,'Teff',3x,
     &     'rhoeff',6x,'geff',6x,'mlos',10x,'M',9x,'dt',11x,'t',10x,
     &     'it',1x,'crash',1x,'nshell',1x,'cputime',/,132('-'),/)
 110  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Values at: center',//,2x,'model',7x,'Tc',6x,'Tmax',6x,
     &     'Tmax_Mr',4x,'rhoc',5x,'rho_max',4x,'Pc',6x,'betac',4x,
     &     'etac',4x,'degpec',5x,'Enu_Pla',6x,'Enucl',6x,
     &     'Egrav',/,130('-'),/)
 120  format (' AGB network : version ',0pf4.2,/,1x,'Energetics :',
     &     1x,'integrated values',//,
     &     2x,'model',7x,'LH',9x,'LHe',9x,
     &     'LC',9x,'LNe',9x,'LO',10x,'LSi',7x,'Lnu_Pla',7x,
     &     'Lnucl',7x,'Lgrav',3x,'irradmax',/,125('-'),/)
 130  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Convective zones I:',//,
     &     2x,'model',5x,'conv1_Mb',2x,'conv1_Rb',1x,'conv1_Tb',1x,
     &     'conv1_rob',1x,'conv1_Mt',2x,'conv1_Rt',1x,'conv1_Tt',1x,
     &     'conv1_rot',4x,'env_Mb',4x,'env_Rb',4x,'env_Tb',4x,'env_rob',
     &     /,128('-'),/)
 140  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Convective zones II:',//,
     &     2x,'model',5x,'conv2_Mb',3x,'conv2_Mt',3x,'conv3_Mb',3x,
     &     'conv3_Mt',3x,'conv4_Mb',3x,'conv4_Mt',3x,'conv5_Mb',3x,
     &     'conv5_Mt',3x,'thermoh_Mb',3x,'thermoh_Mt',2x,'nconvt',/,
c     &     'conv5_Mt',3x,'conv6_Mb',3x,'conv6_Mt',2x,'nconvt',/,
     &     127('-'),/)
 150  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Burning zones : H-burning',//
     &     2x,'model',4x,'Hburn_Mb',2x,'Hburn_Rb',2x,'Hburn_Tb',1x,
     &     'Hburn_rob',2x,'Hburn_Mt',2x,'Hburn_Rt',2x,'Hburn_Tt',1x,
     &     'Hburn_rot',2x,'Hburn_Mm',3x,'Hburn_em',5x,'Lpp',/,
     &     121('-'),/)
 160  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Burning zones : He-burning',//
     &     2x,'model',3x,'Heburn_Mb',1x,'Heburn_Rb',1x,'Heburn_Tb',1x,
     &     'Heburn_rob',1x,'Heburn_Mt',1x,'Heburn_Rt',1x,'Heburn_Tt',1x,
     &     'Heburn_rot',1x,'Heburn_Mm',1x,'Heburn_em',1x,'Heburn_enum',
     &     1x,'Heburn_etam',/,135('-'),/)
 162  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Burning zones : C-burning',//
     &     2x,'model',4x,'Cburn_Mb',2x,'Cburn_Rb',2x,'Cburn_Tb',1x,
     &     'Cburn_rob',2x,'Cburn_Mt',2x,'Cburn_Rt',2x,'Cburn_Tt',1x,
     &     'Cburn_rot',2x,'Cburn_Mm',2x,'Cburn_em',2x,'Cburn_enum',
     &     1x,'Cburn_etam',/,132('-'),/)
 164  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Burning zones : Ne-burning',//,
     &     2x,'model',2x,'Neburn_Mb',1x,'Neburn_Rb',1x,'Neburn_Tb',1x,
     &     'Neburn_rob',1x,'Neburn_Mt',1x,'Neburn_Rt',1x,'Neburn_Tt',1x,
     &     'Neburn_rot',1x,'Neburn_Mm',1x,'Neburn_em',1x,'Neburn_enum',
     &     1x,'Neburn_etam',/,134('-'),/)
 170  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Rotation and Accretion:',
     &     //,2x,'model',5x,'macc',5x,'Racc',5x,'Macc',4x,
     &     'Lshear',4x,'k2conv',3x,'k2rad',5x,'Jconv',5x,'Jrad',
     &     5x,'vsurf',7x,'Usurf',7x,'Fenerg',5x,'dmoment'/,117('-'),/)
 180  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Convective zones II:',//,
     &     2x,'model',5x,'conv2_Rb',4x,'conv2_Rt',4x,'conv3_Rb',4x,
     &     'conv3_Rt',4x,'conv4_Rb',4x,'conv4_Rt',4x,'conv5_Rb',4x,
     &     'conv5_Rt',4x,'therm_Rb',4x,'therm_Rt',/,129('-'),/)
c     &     'conv5_Rt',4x,'conv6_Rb',4x,'conv6_Rt',/,129('-'),/)
 190  format (' AGB network : version ',0pf4.2,/,1x,a7,
     &     ' chemical abundances:',
     &     //,2x,'model',6x,'n ',7x,' H1 ',7x,' H2 ',7x,' He3',7x,
     &     ' He4',7x,' Li6', 7x,' Li7',7x,' Be7',7x,' Be9',7x,' B10',
     &     7x,' B11', /,127('-'),/)
 200  format (' AGB network : version ',0pf4.2,/,1x,a7,
     &     ' chemical abundances:',
     &     //,2x,'model',4x,' C12',7x,' C13',7x,' C14',7x,' N14',7x,
     &     ' N15',7x,' O15', 7x,' O16',7x,' O17',7x,' O18',7x,' F19',
     &      7x, 'Ne20',/,127('-'),/)
 210  format (' AGB network : version ',0pf4.2,/,1x,a7,
     &     ' chemical abundances:',
     &     //,2x,'model',4x,'Ne21',7x,'Ne22',7x,'Na23',7x,'Mg24',7x,
     &     'Mg25',7x,'Mg26', 7x,'Al26m',6x,'Al26g',6x,'Al27',7x,'Si28',
     &     7x,'Si29',/,127('-'),/)
 220  format (' AGB network : version ',0pf4.2,/,1x,a7,
     &     ' chemical abundances:',
     &     //,2x,'model',4x, 'Si30',7x,' P31',7x,' S32',7x,' S33',7x,
     &     ' S34',7x,' S35',7x,'CL35',7x,' S36',7x,'CL36',7x,'CL37',6x,
     &     'heavy',/,127('-'),/)
 230  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Convective envelope turnover times:',//,
     &     2x,'model',4x,'tc',9x,'rc',9x,'tc_hp',6x,'rc_hp',6x,'tc_r',
     &     7x,'rc_r',7x,'tc_m',6x,'rc_m',6x,'tc_max',5x,'rc_max',5x,
     &     'tg',/,117('-'),/)
 240  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Convective core turnover times:',//,
     &     2x,'model',3x,'tc_cc',5x,'rc_cc',5x,'tc_hp_cc',3x,
     &     'rc_hp_cc',3x,'tc_r_cc',3x,'rc_r_cc',4x,'tc_m_cc',5x,
     &     'rc_m_cc',3x,'tc_max_cc',2x,'rc_max_cc',3x,'tg_cc',5x,
     &     'omegac',/,117('-'),/)
 250  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Asteroseimic relations:',//,
     &     2x,'model',5x,'Dnu',8x,'Dnu_ech',3x,'Dnu_error',5x,'Ttot',
     &     8x,'tbce',8x,'tHe',8x,'numax',7x,'Dpg ',/,117('-'),/)

 260  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Rotation and Accretion (bis):',
     &     //,2x,'model',5x,'omegacore',4x,'Jtrue',5x,'Jobs',7x,'Jrad',
     &     5x,'Bequi',3x,'omegamax',3x,'Nmax',4x,'tautot',3x,'taucirc',
     &     2x,'taushear',/,117('-'),/)

 270  format (' AGB network : version ',0pf4.2,/,1x,
     &     'F. Gallet data :',
     &     //,2x,'model',5x,'L/Lsun',4x,'Teff',5x,'time(yr)',7x,'R/Rsun',
     &     5x,'k2conv',3x,'k2rad',3x,'Mrad/Msun',4x,'Rrad/Rsun',/,
     &     117('-'),/)

 1000 format (i8,1x,i1,1x,0pf11.6,1x,0pf7.5,1x,0pf7.5,1x,0pf7.0,1x,
     &     1pe9.3,1x,1pe9.2,1x,1pe10.4,1x,0pf12.8,1x,1pe9.3,1x,
     &     1pe16.10,1x,i2,1x,i2,1x,i4,1x,0pe8.2,' ####',a2)
 1100 format (i8,1x,i1,1x,0pf11.2,1x,0pf7.2,1x,0pf7.2,1x,0pf7.0,1x,
     &     1pe9.3,1x,1pe9.2,1x,1pe10.4,1x,0pf12.8,1x,1pe9.3,1x,
     &     1pe16.10,1x,i2,1x,i2,1x,i4,1x,0pe8.2,' ####',a2)
 1150 format (i8,1x,i1,1x,0pf11.2,1x,0pf7.2,1x,0pf7.0,1x,0pf7.0,1x,
     &     1pe9.3,1x,1pe9.2,1x,1pe10.4,1x,0pf12.8,1x,1pe9.3,1x,
     &     1pe16.10,1x,i2,1x,i2,1x,i4,1x,0pe8.2,' ####',a2)
 1200 format (i8,1x,i1,1x,0pf11.6,1x,0pf7.5,1x,0pf7.5,1x,0pf7.0,1x,
     &     1pe9.3,1x,1pe9.2,1x,1pe10.4,1x,0pf12.8,1x,1pe9.3,1x,
     &     1pe16.10,1x,i2,1x,i2,1x,i4,1x,0pe8.2)
 1300 format (i8,1x,i1,1x,0pf11.2,1x,0pf7.2,1x,0pf7.2,1x,0pf7.0,1x,
     &     1pe9.3,1x,1pe9.2,1x,1pe10.4,1x,0pf12.8,1x,1pe9.3,1x,
     &     1pe16.10,1x,i2,1x,i2,1x,i4,1x,0pe8.2)
 1350 format (i8,1x,i1,1x,0pf11.2,1x,0pf7.2,1x,0pf7.0,1x,0pf7.0,1x,
     &     1pe9.3,1x,1pe9.2,1x,1pe10.4,1x,0pf12.8,1x,1pe9.3,1x,
     &     1pe16.10,1x,i2,1x,i2,1x,i4,1x,0pe8.2)
 1400 format (i8,2x,2(1x,1pe9.3),1x,0pf9.6,3(1x,1pe9.3),1x,
     &     0pf6.4,1x,0pf8.3,1x,0pf8.3,3x,1pe10.4,2(1x,1pe11.4),
     &     1x,a2)
 1410 format (i8,2x,2(1x,1pe9.3),1x,0pf9.6,3(1x,1pe9.3),1x,
     &     0pf6.4,1x,0pf8.1,1x,0pf8.3,3x,1pe10.4,2(1x,1pe11.4),
     &     1x,a2)
 1450 format (1x,'dEpot =',1pe10.3,', dEkin =',1pe10.3,', dEint =',
     &     1pe10.3,', dEphot =',1pe10.3,', dEnupla =',1pe10.3,
     &     ', dEnuc =',1pe10.3,', equil =',1pe9.3)
 1500 format (i8,2x,2(1x,1pe10.4),7(1x,1pe11.4),1x,1pe11.4,a2)
 1600 format (i8,2x,2(1x,0pf9.6,1x,1pe9.3,1x,1pe9.3,1x,1pe9.3),1x,
     &     0pf9.6,1x,1pe11.4,1x,0pf8.4,a2)
 1700 format (i8,2x,2(1x,0pf9.6,1x,1pe9.3,1x,1pe9.3,1x,1pe9.3),1x,
     &     0pf9.6,1x,1pe10.3,1x,1pe10.3,0pf8.3,a2)
 1900 format (i8,2x,1x,0pf10.7,1x,1pe9.3,1x,0pf7.4,1x,0pf8.4,1x,
     &     f10.6,1x,1pe9.3,1x,0pf7.4,1x,0pf8.4,2x,1x,0pf10.6,1x,
     &     1pe9.3,1x,0pf9.5,1x,f9.5,a2)
 2000 format (i8,2x,10(1x,0pf10.6),2x,1x,i2,a1,a2)
 2050 format (i8,2x,10(1x,1pe11.5),a2)
 2100 format (i8,2x,1x,1pe8.2,1x,1pe9.3,1x,0pf6.4,1x,1pe10.4,2(1x,
     &     1pe10.4),1x,1pe11.4,1x,1pe8.2,1x,1pe11.4,1x,1pe10.4,1x,
     &     1pe10.4,1x,1pe11.4,a2)
 2200 format (i8,11(1x,1pe10.4),a2)

 2300 format (i8,2x,1pe8.2,6(1x,1pe8.2),a2)

 2400 format(4(1x,0pf20.12))
 1211 format (i8,1x,22(e15.4E4,1x))
 2500 format(i8,11(1x,1pe10.4),a2)
 2510 format(i8,12(1x,1pe10.4),a2)
 2600 format (i8,1x,1pe10.4,1x,1pe11.4,1x,1pe11.4,1x,1pe11.4,
     &   1x,1pe11.4,1x,1pe11.4,1x,1pe11.4,1x,1pe11.4,a2)
      return
      end

************************************************************************

      SUBROUTINE Quintic(omkyle,vmlt,lmlt,val)

************************************************************************
*   Newton-Raphson solver for the solution to the quintic nondiffusive *
*   rotating convection model of Augustson & Mathis 2019,              *
*                                                                      *
*     Author: Kyle Augustson (CEA Saclay)                              *
*     Date of creation: 26 October 2019.                               *
*     Adaptation to Starevol : TD (04/11/2019)                         *
************************************************************************ 
!     omkyle is the local angular velocity
!     vmlt is the local mixing length velocity
!     lmlt is the local mixing length
!     val is the return value of the 3/2 root of the total wave vector
!
!     Note TD: routine to compute parameter z (see eq. 69, 70)      
************************************************************************
      
      Implicit None
      
      Real*8, Intent(In) :: omkyle, vmlt, lmlt
      Real*8, Intent(InOut) :: val
      Real*8 :: c, err, tol, x, guess, tmp, fun, funp, pi, Roc
      Integer :: n, nmax

      pi = 3.141592653589793238462643383279d0
      
!     Convective Rossby number (eq. 22)      
      Roc = vmlt/(2d0*omkyle*lmlt)
! c is the third term of equation (46)      
      c = 18d0/(5d0*pi*Roc)**2
      
      if (c.gt.20d0) then
         guess = 0.65d0*(c**(0.2d0))
      else
         guess = 1.1d0
      end if
      nmax = 100
      x = guess
!     Newton-Raphson algorithm
      tol = 1d-15
      err = 1d0
      n=0
      do while ((err>tol) .and. (n<nmax))
         fun = c+5d0*x**2-2d0*x**5
         funp = 10d0*x*(1d0-x**3)
         tmp = x - fun/funp
         err = abs(tmp-x)
         x = tmp
         n = n + 1
      end do
      val = tmp

      return
      end
c$$$      end subroutine Quintic
      SUBROUTINE resulpr (error)

************************************************************************
* Save computed model and write display file                           *
* 18/06 : adding some writings into evoldisp in case of rotation       *
* Modifs CC ondes (2/04/07 -> 23/11/07)                                *
* $LastChangedDate:: 2016-05-11 17:20:46 +0200 (Mer, 11 mai 2016)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 62                                                          $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eng'
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
      include 'evolcom.transp'
      include 'evolcom.var'
C Modifs CC ondes (2/04/07)
      include 'evolcom.igw'
c      include 'evolcom.transpondesexcit'

      integer error
      integer imin,imax
      integer j,k,kl,l,kl1
C Modifs CC ondes (2/04/07) -->
      integer kk
      double precision lgdmic,lgdturb
      double precision dturb,dmic,vdmic
C <--

      double precision nturb
      double precision Dhold
      double precision dift,dife
      double precision vom,vrray,vxpsi
      double precision abmurj,abmuj
      double precision FeH
      double precision tautime,Ereac,taureac
      double precision Ushear(nsh),ULj(nsh)
      double precision brunt,lumwaves(nsh),bruntV2,bruntV

      double precision depottotn,depottot_surfn,depottot_coren
      double precision depotondestotn,vrn

      double precision abundelec     ! Ajout pour creer donnees Montreal (Mars 2018)
      
      character*1 crzc(nsh)
      logical partialmix

      common /calcDh/ Dhold(nsh)
      common /difcirc/ dift(nsh),dife(nsh)
      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)
      common /gradmuj/ abmurj(nsh),abmuj(nsh)
      common /metal/ FeH

      common /brunt/ brunt(nsh),bruntV2(nsh),bruntV(nsh)
      common /abundelec/ abundelec(nsh)       ! Ajout pour creer donnees Montreal (Mars 2018)
      
C Modifs CC ondes (2/04/07) -->
      common /coefdiffpaqhe/ xvdiffhe(nsh),xvdiffhenoc(nsh)
      common /coefdiff/ dturb(nsh),dmic(nsh),vdmic(nsh),lgdturb(nsh)
     &     ,lgdmic(nsh)
      double precision xvdiffhe,xvdiffhenoc
C     <--
c Ajout TD (18/01/2019)
      common /diffutest/ zmeanO16(nsh),dens_elec(nsh) ,coulomb(nsh)
     $     ,Zpaq(nsh),Zpaq1(nsh),Zpaq2(nsh),Kpaq(nsh)
      double precision zmeanO16,dens_elec,coulomb
     $     ,Zpaq,Zpaq1,Zpaq2,Kpaq
c Fin ajout      
c Ajout par TD (22/02/2018) + modif 25/04/2019    
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
      
      common /dturbulente/ nturb(0:nsh)
      common /nucleaire/ tautime(nsh,nreac),Ereac(nsh,nreac),taureac(nsh
     $     ,nreac)
      
      dimension depottotn(nsh),depottot_surfn(nsh),depottot_coren(nsh)
      dimension depotondestotn(nsh),vrn(nsh)

      if (imodpr.eq.11.or.imodpr.eq.12)
     &     call myflamespeed (.false.,imodpr)
      
*________________________________________________________
***   storage of the calculated models in the binary file
*--------------------------------------------------------

      if (error.eq.-30) neff = -neff

      if (nsconv.gt.0.and.turnover(1).eq.0.d0) then
c..   determine turnover timescale
         do kl = 1,nsconv
            imin = novlim(kl,3)
            imax = novlim(kl,4)
            call mix (dm,imin,imax,kl,5,partialmix)
         enddo
      endif
      
      do k = 1,nmod
         if (crz(k).eq.-2) crzc(k) = 'c'
         if (crz(k).eq.-1) crzc(k) = 'a'
         if (crz(k).eq.-3) crzc(k) = 's'
         if (crz(k).eq.3) crzc(k) = 'S'
         if (crz(k).eq.2) crzc(k) = 't'
         if (crz(k).eq.1) crzc(k) = 'o'
         if (crz(k).eq.4) crzc(k) = 'r'
      enddo

      call system ('cp -f nextini.bin oldini.bin')
      rewind (92)
      write (92) model,nphase,totm,time,dtn,neff,code_version,
     &     nmod,nis-1,mdredgeup,FeH,turnenv,rayon0,
     &     (crzc(k),u(k),r(k),lnf(k),lnt(k),lum(k),
     &     vu(k),vr(k),vlnf(k),vlnt(k),vlum(k),facc(k),
     &     exposure(k),mueinv(k),dm(k),vhconv(k),omega(k),
     &     ro(k),vro(k),p(k),vp(k),e(k),ve(k),enucl(k),s(k),
     &     vs(k),sconv(k),pvisc(k),tau(k),
     &     (xsp(k,l),l = 1,nis-1),k = 1,nmod),
     &     (mtotlos(j),j = 1,nis-1)
      close (92)
      
      do k = 1,nmod
         venucl(k) = enucl(k)
         vsconv(k) = sconv(k)
         vpvisc(k) = pvisc(k)
         vmueinv(k) = mueinv(k)
      enddo
      if (irotbin.eq.1) then
c         call system ('cp -f nextang.bin modang.bin')
         rewind (94)
c         rewind(83)
         write (94) ndtold,ndbold,ur_Ss,xmom_tots,inits,
     &        (auxs(k),urs(k),vxpsi(k),xpsis(k),V_circ(k),
     &        xalpha(k),xlambdas(k),xKt(k),dift(k),dife(k),
     &        Dhold(k),Dconv(k),abmurj(k),xnuvv(k),xnum(k),
     &        k = 1,nmod),
     &        ((vxsp(k,l),l = 1,nis-1),k = 1,nmod)
C Modifs CC ondes (2/04/07) -->
c         print *,'urs=',urs(1:nmod)
         do k=1,nmod
         kk=nmod-k+1
            if (idiffcc) then
                lgdturb(k) = 0.d0
                lgdmic(k) = 0.d0
                if (dturb(k).gt.0.d0) lgdturb(k) = log10(dturb(k))
c                if (xdmiche(k).gt.0.d0) lgdmic(k) = log10(xdmiche(k))  ! modif Mai 2019
            else
                lgdturb(k) = 0.d0
                lgdmic(k) = 0.d0
            endif
         enddo
C <--
      endif
      open (unit = 92,file = 'nextini.bin',form = 'unformatted',
     &     status = 'unknown')



*_________________
***   print banner
*-----------------

      write (nout,*)
      if (nsconv.gt.0) then
         write (nout,1100)
      endif
      if (.not.ledouxmlt) then
         do j = 1,nsconv
            novlim(j,1) = novlim(j,3)
            novlim(j,2) = novlim(j,4)
         enddo
         nschwar = nsconv
      endif
      do kl = 1,nschwar
         if (novlim(kl,2)-novlim(kl,1).gt.5) then
            if ((diffover.or.novopt.gt.0).and..not.ledouxmlt.and.
     &           (novlim(kl,8).gt.novlim(kl,2).or.
     &           novlim(kl,1).gt.novlim(kl,7))) then
               write (90,1200) novlim(kl,1),novlim(kl,2),
     &              m(novlim(kl,1))/msun,m(novlim(kl,2))/msun,
     &              turnover(kl)*seci,novlim(kl,7),novlim(kl,8),
     &              (m(novlim(kl,1))-m(novlim(kl,7)))/msun,
     &              (m(novlim(kl,8))-m(novlim(kl,2)))/msun
               write (nout,1200) novlim(kl,1),novlim(kl,2),
     &              m(novlim(kl,1))/msun,m(novlim(kl,2))/msun,
     &              turnover(kl)*seci,novlim(kl,7),novlim(kl,8),
     &              (m(novlim(kl,1))-m(novlim(kl,7)))/msun,
     &              (m(novlim(kl,8))-m(novlim(kl,2)))/msun
            else
               write (90,1210) novlim(kl,1),novlim(kl,2),
     &              m(novlim(kl,1))/msun,m(novlim(kl,2))/msun,
     &              turnover(kl)*seci
               write (nout,1210) novlim(kl,1),novlim(kl,2),
     &              m(novlim(kl,1))/msun,m(novlim(kl,2))/msun,
     &              turnover(kl)*seci
            endif
         endif
      enddo
      if ((diffover.or.novopt.gt.0).and.ledouxmlt) then
         do kl = 1,nledoux
            kl1 = 0
            do l = 1,nschwar
               if (novlim(kl,8).gt.novlim(kl,6).or.
     &              novlim(kl,5).gt.novlim(kl,7)) then
                  kl1 = l
                  exit
               endif
            enddo
            if (kl1.ne.0) then
               write (90,1240) novlim(kl,5),novlim(kl,6),
     &              m(novlim(kl,5))/msun,m(novlim(kl,6))/msun,
     &              turnover(kl)*seci,novlim(kl1,7),novlim(kl1,8),
     &              (m(novlim(kl,5))-m(novlim(kl1,7)))/msun,
     &              (m(novlim(kl1,8))-m(novlim(kl,6)))/msun
               write (nout,1240) novlim(kl,5),novlim(kl,6),
     &              m(novlim(kl,5))/msun,m(novlim(kl,6))/msun,
     &              turnover(kl)*seci,novlim(kl1,7),novlim(kl1,8),
     &              (m(novlim(kl,5))-m(novlim(kl1,7)))/msun,
     &              (m(novlim(kl1,8))-m(novlim(kl,6)))/msun
            else
               write (90,1230) novlim(kl,5),novlim(kl,6),
     &              m(novlim(kl,5))/msun,m(novlim(kl,6))/msun,
     &              turnover(kl)*seci
               write (nout,1230) novlim(kl,5),novlim(kl,6),
     &              m(novlim(kl,5))/msun,m(novlim(kl,6))/msun,
     &              turnover(kl)*seci
            endif
         enddo
      endif
      do kl = 1,nsemiconv
         if (novlim(kl,10)-novlim(kl,9).gt.1.or.
     &        crz(novlim(kl,10)+1).lt.0) then
            write (90,1250) novlim(kl,9),novlim(kl,10),
     &           m(novlim(kl,9))/msun,m(novlim(kl,10))/msun
            write (nout,1250) novlim(kl,9),novlim(kl,10),
     &           m(novlim(kl,9))/msun,m(novlim(kl,10))/msun
         endif
      enddo
      do kl = 1,nthermha
         if (novlim(kl,12)-novlim(kl,11).gt.1.or.
     &        crz(novlim(kl,12)+1).lt.0) then
            write (90,1260) novlim(kl,11),novlim(kl,12),
     &           m(novlim(kl,11))/msun,m(novlim(kl,12))/msun
            write (nout,1260) novlim(kl,11),novlim(kl,12),
     &           m(novlim(kl,11))/msun,m(novlim(kl,12))/msun
         endif
      enddo

c..   Storage of critical Reynolds number, critical Richardson number
c..   and prescription used for horizontal turbulent diffusion coefficient

      if (model.eq.(modeli+1).and.rotation) then
         write (90,1270) Dh_prescr
      endif


      write (90,1000) totm,dms,time*seci,dtn*seci,teff,soll

c     if (no.eq.maxmod.or.time.eq.4.6000000d9*sec) then
      if (no.eq.maxmod) then
            if (igw) then

               do k = 1,nmod
                  if (dabs(lumondestot(k)).lt.1d-99)
     &                 lumondestot(k) = 0.d0
                  if (dabs(lumondes_surf(k)).lt.1d-99)
     &                 lumondes_surf(k) = 0.d0
                  if (dabs(lumondes_core(k)).lt.1d-99)
     &                 lumondes_core(k) = 0.d0
                  brunt(k)=brunt_o(nmod-k+1)
                  lumwaves(k)=lumwave(nmod-k+1)
               enddo
               rewind (85)

               write (85) (lgdmic(k),lgdturb(k),vom(k),
     &              depottot(k)/lsun,depottot_surf(k)/lsun,
     &              depottot_core(k)/lsun,
     &              depotwaves(k)/lsun,Dondes(k),Dondeschim(k),
     &              lumwaves(k),lumwave_core(k),
     &              lumondes_surf(k),lumondes_core(k),lumondestot(k),
     &              brunt(k),vr(k)/rsun,vxpsi(k),           
     &              Dmicro(k),vmicro(k),Dthc(k),phiKS(k),
     &              deltaKS(k),tautime(k,ippg),Ereac(k,ippg),tautime(k
     &              ,ipdg), Ereac(k,ipdg),tautime(k,i2he3),Ereac(k
     &              ,i2he3),tautime(k,ihe3ag), Ereac(k,ihe3ag),tautime(k
     &              ,ibe7beta),Ereac(k,ibe7beta), tautime(k,ili7pa)
     &              ,Ereac(k,ili7pa),tautime(k,ibe7pg), Ereac(k,ibe7pg)
     &              ,tautime(k,ib8beta),Ereac(k,ib8beta), tautime(k
     &              ,ic13pg),Ereac(k,ic13pg),tautime(k,in14pg), Ereac(k
     &              ,in14pg),tautime(k,icpg),Ereac(k,icpg), k = 1,nmod)

            else
               write (85) (Dmicro(k),vmicro(k),Dthc(k),phiKS(k),
     $              Dturbul(k),Dbar(k),Dkyle(k),coefDtacho(k),                                ! dturbul 10/2019 + dbar
     $              xvdiffhe(k),xvdiffhenoc(k),xvdiffc12(k)    ! Ajout TD Fev.2018
     $              ,xvdiffc12noc(k),xvdiffc13(k)
     $              ,xvdiffc13noc(k), xvdiffn14(k),xvdiffn14noc(k)
     $              ,xvdiffn15(k),xvdiffn15noc(k),
     $              xvdiffo16(k),xvdiffo16noc(k)
     $              ,xvdiffo17(k),xvdiffo17noc(k),xvdiffo18(k)
     $              ,xvdiffo18noc(k), xvdiffne20(k),xvdiffne20noc(k),
     $              xvdiffna23(k),xvdiffna23noc(k),
     $              deltaKS(k),tautime(k,ippg),Ereac(k        ! Fin ajout TD Fev.2018
     $              ,ippg), tautime(k,ipdg), Ereac(k,ipdg),tautime(k
     $              ,i2he3), Ereac(k,i2he3),tautime(k,ihe3ag),Ereac(k
     $              ,ihe3ag), tautime(k,ibe7beta),Ereac(k,ibe7beta),
     $              tautime(k,ili7pa),Ereac(k,ili7pa),tautime(k,ibe7pg),
     $              Ereac(k,ibe7pg),tautime(k,ib8beta),Ereac(k,ib8beta),
     $              tautime(k,ic13pg),Ereac(k,ic13pg),tautime(k,in14pg),
     $              Ereac(k,in14pg),tautime(k,icpg),Ereac(k,icpg), k = 1
     $              ,nmod)
            endif
      endif

     
c..   Storage of values for atomic diffusion (modif AP TD Jan.2019)
      if(model.eq.(modeli+1).and.microdiffus) then
         do k=1,nmod
            write(654,'(1x,i4,9(1x,1pe11.4))'),k,T(k),zmeanO16(k),xsp(k
     $           ,io16),dens_elec(k),coulomb(k),Zpaq(k),Zpaq1(k)
     $           ,Zpaq2(k),Kpaq(k)
         enddo
      endif
      
      
 1000 format (' Final mass = ',f11.8,', dm = ',f11.9,', age = ',1pe12.6,
     &     ' yr, dt = ',1pe13.7,' yr, Teff = ',0pf7.0,', L = ',1pe10.4)
 1100 format (5x,'convective limits:')
 1200 format (2x,'Schwar.: ',i4,' -->',i4,' (',1pe9.3,',',1pe9.3,
     &     ')   turnover: ',1pe9.3,'yr   Overshoot (dm): ',i4,' -->',i4,
     &     ' (',1pe9.3,',',1pe9.3,')')
 1210 format (2x,'Schwar.: ',i4,' -->',i4,' (',1pe9.3,',',1pe9.3,
     &     ')   turnover: ',1pe9.3,'yr')
 1230 format (2x,'Ledoux : ',i4,' -->',i4,' (',1pe9.3,',',1pe9.3,
     &     ')   turnover: ',1pe9.3,'yr')
 1240 format (2x,'Ledoux : ',i4,' -->',i4,' (',1pe9.3,',',1pe9.3,
     &     ')   turnover: ',1pe9.3,'yr   Overshoot (dm): ',i4,' -->',i4,
     &     ' (',1pe9.3,',',1pe9.3,')')
 1250 format (2x,'semiconvection : ',i4,' -->',i4,' (',1pe9.3,',',
     &     1pe9.3,')')
 1260 format (2x,'Thermohaline : ',i4,' -->',i4,' (',1pe9.3,',',
     &     1pe9.3,')')
 1270 format (//,'o PRESCRIPTIONS FOR ROTATIONNAL MIXING',/,3x,
     &     'Critical Richardson number Ric = 0.25',/,3x,
     &     'Prescription for horizontal turbulence diffusion',
     &     ' coefficient Dh : ',1x,a8,/)
 8100 format('# nsh lgdmic lgdturb vomega depottot depottot_surf',
     &     'depottot_core depotwaves Dondes Dondeschim lumwave',
     &     'lumwave_core brunt_o lumondes vr vpsi',
     &     'Dmicro',6x,'Vmicro',5x,'Dthc',6x,'phiKS',6x
     $     ,'deltaKS',6x,'Tpp',6x,'Epp',6x,'Tpd',6x,'Epd',6x,'T2he3',6x
     $     ,'E2he3',6x,'The3he4',6x,'Ehe3he4',6x,'Tbe7b',6x,'Ebe7b',6x
     $     ,'Tli7p' ,6x,'Eli7p' ,6x,'Tbe7p',6x,'Ebe7p',6x,'Tb8b',6x
     $     ,'Eb8b',6x ,'Tc13p',6x ,'Ec13p',6x ,'Tn14p',6x,'En14p',6x
     $     ,'Tc12p',6x ,'Ec12p')
 8101 format(4x,'i',4x,'Dmicro',6x,'Vmicro',5x,'Dthc',6x,'phiKS',6x
     $     ,'deltaKS',6x,'Tpp',6x,'Epp',6x,'Tpd',6x,'Epd',6x,'T2he3',6x
     $     ,'E2he3',6x,'The3he4',6x,'Ehe3he4',6x,'Tbe7b',6x,'Ebe7b',6x
     $     ,'Tli7p' ,6x,'Eli7p' ,6x,'Tbe7p',6x,'Ebe7p',6x,'Tb8b',6x
     $     ,'Eb8b',6x ,'Tc13p',6x ,'Ec13p',6x ,'Tn14p',6x,'En14p',6x
     $     ,'Tc12p',6x ,'Ec12p')


 8102 format (2x,i4,5(1x,1pe11.4),26(1x,1pe11.4),5(1x,1pe11.4),
     &     22(1x,e15.4E2))
 8103 format (2x,i4,5(1x,1pe11.4),22(1x,1pe14.6))

      return
      end


************************************************************************

      SUBROUTINE rinimod (uread,ureadrot,uwrite,uwriterot,flag)

************************************************************************
* Read the initial stellar model                                       *
* 09/03: mass fraction of photons & CAPTN not stored in binary anymore *
* 23/09: arguments added to the subroutine name to manage reading and  *
*        writing of binary files, including crash cases                *
*  flag  = 0 : first reading                                           *
*  flag <> 0 : failed computation, reload previous model               *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.ion'
      include 'evolcom.lum'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.opa'
      include 'evolcom.rot'
      include 'evolcom.conv'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.transp'
      include 'evolcom.var'

      integer nspr
      integer error,last_model
      integer i,j,k,l,kl,jl,it
      integer uread,ureadrot,uwrite,uwriterot,flag
      integer ishockb,ishockt,itmax
      integer klenv,klpulse,klcore

      character*1 crzc(nsh)
      logical tbot,ttop,ierr,test

      double precision vom,vrray,vxpsi
      double precision dift,dife,Dhold
      double precision xc11,xsi31,xsi32,xp32,xp33,xar36,FeH
      double precision abmurj,abmuj
      double precision zstart,y,xspr,vxspr,mtotlosr
      double precision Lmax,tmax,Mackmax,enucmax,elemsave
      double precision phim,deltam
      double precision vvar(4,nsh)
      double precision eta,degpe,rhpsi
      double precision fscr

      common /saveMS/ elemsave
      common /hydrodyn/ Lmax,tmax,Mackmax,enucmax,itmax,ishockb,ishockt
      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)
      common /difcirc/ dift(nsh),dife(nsh)
      common /calcDh/ Dhold(nsh)
      common /gradmuj/ abmurj(nsh),abmuj(nsh)
      common /metal/ FeH
      common /edegen/ eta(nsh),degpe(nsh),rhpsi(nsh)
      common /funcscr/ fscr(nsh,nscr)
      common /overshoot/ klcore,klenv,klpulse

      dimension y(nsp),xspr(nsh,100),vxspr(nsh,100),mtotlosr(100)

***   read binary models
      
      rewind (uread)

      read (uread) model,nphase,totm,time,dtn,neff,bin_version,nmod,
     &     nspr,mdredgeup,FeH,turnenv,rayon0,
     &     (crzc(k),u(k),r(k),lnf(k),lnt(k),lum(k),
     &     vu(k),vr(k),vlnf(k),vlnt(k),vlum(k),vfacc(k),
     &     exposure(k),vmueinv(k),dm(k),vhconv(k),vomega(k),
     &     ro(k),vro(k),p(k),vp(k),e(k),ve(k),venucl(k),s(k),vs(k),
     &     vsconv(k),vpvisc(k),tau(k),
     &     (xspr(k,l),l = 1,nspr),k = 1,nmod),
     &     (mtotlosr(j),j = 1,nspr)
      if (irotbin.eq.1) then
         read (ureadrot) ndtold,ndbold,ur_Ss,xmom_tots,inits,
     &        (auxs(k),urs(k),vxpsi(k),xpsis(k),V_circ(k),
     &        xalpha(k),xlambdas(k),xKt(k),dift(k),dife(k),
     &        Dhold(k),Dconv(k),abmurj(k),xnuvv(k),xnum(k),k = 1,
     &        nmod),((vxspr(k,l),l = 1,nspr),k = 1,nmod)
      endif
      ierr = .false.
      if (neff.lt.0) then
         ierr = .true.
         neff = abs(neff)
      endif
      
c..   compatibility check
      if (nphase.gt.1.and.mdredgeup.eq.0.d0) mdredgeup = totm*msun
      if (bin_version.lt.2.95d0) then
         write (*,*) '*** Initialize neutron exposure ***'
         forall (i = 1:nmod) exposure(i) = 0.d0
         mdredgeup = totm*msun
      endif
      do k = 1,nmod
         crz(k) = 0
         if (crzc(k).eq.'c') crz(k) = -2 
         if (crzc(k).eq.'a') crz(k) = -1 
         if (crzc(k).eq.'s') crz(k) = -3 
         if (crzc(k).eq.'l'.or.crzc(k).eq.'S') crz(k) = 3 
         if (crzc(k).eq.'t') crz(k) = 2 
         if (crzc(k).eq.'o') crz(k) = 1 
         if (crzc(k).eq.'r') crz(k) = 4 
         if (crz(k).eq.0) then
            write (*,'("Pb shell ",i4," crzc = ",a1," not recognized")')
     &           k,crzc(k)
            stop
         endif
      enddo
      
****************************************************
* check chemical composition of first sequence model
****************************************************

      if (flag.eq.0) then
         inquire (file = 'last_model', exist = test)
         if (test) then
            open (unit = 9, file = 'last_model')
            read (9,*) last_model
            close (9)
            if (model.ne.last_model) then
               write (*,*) 'The previous and new model numbers are ',
     &              'not consecutive : check binary/sequence name'
               write (90,*) 'The previous and new model numbers are ',
     &              'not consecutive : check binary/sequence name'
               stop 'rinimod'
            endif            
         endif
         write (nout,100) zkint,FeH
         write (90,100) zkint,FeH
 100     format (1x,' INITIAL COMPOSITION : Z = ',0pf8.6,' , [Fe/H] = ',
     &        0pf7.4)
      endif
      if (nspr+2.gt.nsp) then
         write (nout,*) 'binary conversion : MASSIVE --> EVOL format'
         do k = 1,nmod
            do l = 1,51
               xsp(k,l) = xspr(k,l)
            enddo
            xar36 = xspr(k,52)
            xsp(k,52) = xspr(k,53)
            xsp(k,53) = xar36
            do l = 54,nspr
               xsp(k,53) = xspr(k,53)+xspr(k,l)
            enddo
            if (irotbin.eq.1) then
               do l = 1,51
                  vxsp(k,l) = vxspr(k,l)
               enddo
               xar36 = vxspr(k,52)
               vxsp(k,52) = vxspr(k,53)
               vxsp(k,53) = xar36
               do l = 54,nspr
                  vxsp(k,53) = vxspr(k,53)+vxspr(k,l)
               enddo
            endif
         enddo
         do l = 1,51
            mtotlos(l) = mtotlosr(l)
         enddo
         xar36 = mtotlosr(52)
         mtotlosr(52) = mtotlosr(53)
         mtotlosr(53) = xar36
         do l = 1,nspr
            mtotlos(53) = mtotlosr(53)+mtotlosr(l)
         enddo
      else
         do l = 1,nspr
            mtotlos(l) = mtotlosr(l)
         enddo
         do l = 1,nis-1
            do k = 1,nmod
               xsp(k,l) = xspr(k,l)
            enddo
         enddo
         if (irotbin.eq.1) then
            do l = 1,nis-1
               do k = 1,nmod
                  vxsp(k,l) = vxspr(k,l)
               enddo
            enddo
         endif
         do l = 1,nis-1
            mtotlos(l) = mtotlosr(l)
         enddo
      endif

      mtotlos(nis) = 1.d-50
      mtotlos(nsp) = 1.d-50
      do k = 1,nmod
         xsp(k,nis) = 1.d-50
         xsp(k,nsp) = 1.d-50
      enddo


***   make correspondance between old and new binary files
      if (bin_version.gt.5.d2.or.bin_version.eq.0.d0) then
         write (nout,*) 'species conversion in binary : NETWORK change'
         do k = 1,nsp
            y(k) = xsp(1,k)
         enddo
         do k = 1,nmod
            xc11 = xsp(k,13)
            xsi31 = xsp(k,40)
            xsi32 = xsp(k,42)
            xp32 = xsp(k,43)
            xp33 = xsp(k,45)
            xsp(k,46) = xsp(k,46)+xp33
            xsp(k,45) = xsp(k,44)+xp32+xsi32
            xsp(k,44) = xsp(k,41)+xsi31
            xsp(k,43) = xsp(k,39) ! SI 30
            xsp(k,42) = xsp(k,38) ! SI 29
            xsp(k,41) = xsp(k,37) ! SI 28
            xsp(k,40) = xsp(k,36) ! AL 27
            xsp(k,39) = 1.d-50    ! Mg 27
            xsp(k,38) = xsp(k,35)
            xsp(k,37) = xsp(k,34)
            xsp(k,36) = xsp(k,33)
            xsp(k,35) = xsp(k,32)
            xsp(k,34) = 1.d-50    ! Na 25
            xsp(k,33) = xsp(k,31)
            xsp(k,32) = 1.d-50    ! Na 24
            xsp(k,31) = xsp(k,30)
            xsp(k,30) = 1.d-50    ! Ne 23
            do j = 13,24
               xsp(k,j) = xsp(k,j+1)
            enddo
            xsp(k,25) = 1.d-50    ! F  20
            xsp(k,12) = xsp(k,12)+xc11
         enddo
         do k = 1,nsp
            write (nout,200) k,k,y(k),k,xsp(1,k)
 200        format(' # :',i2,', old xsp(1,',i2,') = ',
     &           1pe9.3,' --> new xsp(1,',i2,') = ',1pe9.3)
         enddo
      endif
      
c..   initializations and variables definition

      r(1) = 1.d-99
      vr(1) = 1.d-99
      lnr(1) = -227.d0
      vlnr(1) = -227.d0
      enucl(1) = venucl(1)
      m(1) = 0.d0
      do k = 2,nmod
         m(k) = m(k-1)+dm(k-1)
         lnr(k) = log(r(k))
         vlnr(k) = log(vr(k))
         enucl(k) = venucl(k)
      enddo

      write(899,*)'IN RINIMOD'
      do i = 1,nmod
         write(899,*) i,vomega(i)
      enddo

      if (flag.eq.0) then
         rewind (uwrite)
         write (uwrite) model,nphase,totm,time,dtn,neff,bin_version,
     &        nmod,nis-1,mdredgeup,FeH,turnenv,rayon0,
     &        (crzc(k),u(k),r(k),lnf(k),lnt(k),lum(k),
     &        vu(k),vr(k),vlnf(k),vlnt(k),vlum(k),vfacc(k),
     &        exposure(k),vmueinv(k),dm(k),vhconv(k),vomega(k),
     &        ro(k),vro(k),p(k),vp(k),e(k),ve(k),venucl(k),s(k),vs(k),
     &        vsconv(k),vpvisc(k),tau(k),
     &        (xsp(k,l),l = 1,nis-1),k = 1,nmod),
     &        (mtotlos(j),j = 1,nis-1)
         if (irotbin.eq.1) then
            rewind (uwriterot)
            write (uwriterot) ndtold,ndbold,ur_Ss,xmom_tots,inits,
     &           (auxs(k),urs(k),vxpsi(k),xpsis(k),V_circ(k),
     &           xalpha(k),xlambdas(k),xKt(k),dift(k),dife(k),
     &           Dhold(k),Dconv(k),abmurj(k),xnuvv(k),xnum(k),
     &           k = 1,nmod),((vxsp(k,l),l = 1,nis-1),k = 1,nmod)
         endif
                  
c         if (time.gt.4.76d17) stop 'rinimod : age exceeds 15 Gyr !'
      endif

*____________________________
***   various initializations
*----------------------------

c..   stop microscopic diffusion if nphase > 4
      if (nphase.gt.4.and.lmicro.ne.0) then
         lmicro = 0
         microdiffus = .false.
         print *, 'microdiffus become false nphase > 4'
      endif

      nmod1 = nmod-1
      dmaccr = 0.d0
      error = 0
      mr(1) = 0.d0
      forall (k = 1:nmod)
         mr(k) = m(k)/m(nmod)
         t(k) = exp(lnt(k))
         vt(k) = exp(vlnt(k))
      end forall
      

***   initialize network parameters + thermodynamical variables
      if (flag.eq.0) then
         call rininet
         call eos (0,error)

c.. define initial H and He abundances for additional save
         if (nphase.le.2) then
            elemsave = max(xsp(1,ih1)-0.05d0,0.d0)
         elseif (nphase.le.4) then
            elemsave = max(xsp(1,ihe4)-0.05d0,0.d0)
         else
            elemsave = 0.d0
         endif

***   initialize model parameters
         no = 0      ! current model between [1,nmaxmod]
         if (ierr.and.icorr.eq.'t') then
            ifail = -1
         else
            ifail = 0
         endif

         iter = 1
c         if (nphase.eq.1) ishtest = .false.
         if (nphase.eq.2.or.nphase.eq.4.or.nphase.eq.6) ishtest = 't'
         if (nphase.lt.2) then
            if (abs(log(mtini/totm)).gt.log(1.02d0).and.iaccr.eq.0) then
               write (nout,300) totm
 300           format ('change mtini in starevol.par, does not ',
     &              'correspond to stellar mass model = ',0pf8.5)
               stop ' rinimod'
            endif
         endif
c..   check if metallicity table is correct
         Zstart = 1.d0-xsp(nmod,ih1)-xsp(nmod,ih2)-xsp(nmod,ihe3)-
     &        xsp(nmod,ihe4)
         
         if (((abs(log(zstart/zkint)).gt.log(1.15d0).and.Zstart.gt.
     &        1.d-6).or.(Zstart.lt.1.d-6.and.zkint.gt.1.d-6)).and.
     &        nphase.lt.3.and.(.not.microdiffus.or.crzc(1).eq.'c')) then
            write (*,400) zkint,zstart
 400        format (/,2x,'WARNING : Zkint in starevol.par (',f8.6,')'
     &           ' does not correspond to stellar composition ',f8.6,/)
c            if (nphase.le.2) stop ' rinimod'
         endif

***   initialize nuclear luminosities (needed for init)
         totgrav = 0.d0
         totnucl = 0.d0
         vlpp = 0.d0
         vlh = 0.d0
         vlhe = 0.d0
         vlc = 0.d0
         vlne = 0.d0
         vlo = 0.d0
         vlsi = 0.d0
         lpp = 0.d0
         lh = 0.d0
         lhe = 0.d0
         lc = 0.d0
         lne = 0.d0
         lo = 0.d0
         lsi = 0.d0
          
***   initialize thermo (needed for screening factors, rotation)
         forall (l = 1:nsp, k = 1:nmod) ysp(k,l) = xsp(k,l)/anuc(l)
         forall (k = 1:nmod )
            vvar(1,k) = e(k)
            vvar(2,k) = p(k)
            vvar(3,k) = s(k)
            vvar(4,k) = ro(k)
         end forall
         call eos (1,error)
         forall (k = 1:nmod )
            e(k) = vvar(1,k)
            p(k) = vvar(2,k)
            s(k) = vvar(3,k)
            ro(k) = vvar(4,k)
         end forall
c..   initialize screening factors
         call fscreen (t,ro,xsp,rhpsi,anuc,znuc)

      endif

***   initialize vmueinv for URCA process if old binary is used
      if (vmueinv(1).eq.0.d0) then
         write (nout,*) ' WARNING : binary change, compute mueinv'
         write (90,*) ' WARNING : binary change, compute mueinv'
         forall (k = 1:nmod) vmueinv(k) = mueinv(k)
      endif

***   initialize variables for ROTATION
      if (rotation.and.bin_version.lt.2.8d0) then
         write (nout,*) ' WARNING : binary modification, compute PSI'
         write (90,*) ' WARNING : binary modification, compute PSI'
         if (numeric.eq.2.or.numeric.eq.3) then
c.. first order accuracy in the spatial derivatives
            forall (i = 1:nmod)        
               wi(i) = 0.5d0
               wj(i) = 0.5d0
            end forall
         else
c.. second order accuracy in the spatial deriatives
            wi(1) = 0.5d0
            wj(1) = 0.5d0
            forall (i = 2:nmod1)
               wi(i) = dm(i)/(dm(i)+dm(i-1))
               wj(i) = 1.d0-wi(i)
            end forall
            wi(nmod) = 0.5d0
            wj(nmod) = 0.5d0
         endif
         do k = 2,nmod
            deltam = deltaKS(k-1)*wi(k)+deltaKS(k)*wj(k)
            phim = phiKS(k-1)*wi(k)+phiKS(k)*wj(k)
            xpsis(k) = (phim*xlambdas(k)-xpsis(k)*vomega(k)
     &           /vomega(nmod))/deltam
         enddo
         xpsis(1) = xpsis(2)
      endif

***   define novlim (needed for mesh)
      novlim(1:nmaxconv,1:ntypeconv) = 0
      nsconv = 0
      nschwar = 0
      nledoux = 0
      nsemiconv = 0
      nthermha = 0

***   define Schwarzschild and Ledoux boundaries
      kl = 1
      if (crz(1).lt.-1) novlim(kl,1) = 1
      do k = 2,nmod1
         tbot = crz(k).gt.0.and.crz(k+1).lt.-1
         ttop = crz(k).lt.-1.and.crz(k+1).gt.-2
         if (tbot) novlim(kl,1) = k+1
         if (ttop) then
            if (novlim(kl,1).lt.k) then
               novlim(kl,2) = k
               kl = kl+1
               if (kl.gt.nmaxconv) goto 20
            endif
         endif
      enddo
 20   nschwar = kl-1
      if (nschwar.gt.0) then
         if (novlim(nschwar,2).eq.0) then
            write (*,*)'Problem in determining Schwarzschild boundaries'
            stop 
         endif
      endif
      kl = 1
      it = -3
      if (ledouxmlt) it = 3
      if (crz(1).lt.-1) novlim(kl,5) = 1
      do k = 2,nmod1
         tbot = crz(k).ne.-2.and.crz(k+1).eq.-2
         ttop = crz(k).eq.-2.and.crz(k+1).ne.-2
         if (tbot) novlim(kl,5) = k+1
         if (ttop) then
            if (novlim(kl,5).lt.k) then
               novlim(kl,6) = k
               kl = kl+1
               if (kl.gt.nmaxconv) goto 30
            endif
         endif
      enddo
 30   nledoux = kl-1
      if (nledoux.gt.0) then
         if (novlim(nledoux,6).eq.0) then
            write (*,*) 'Problem in determining Ledoux boundaries'
            stop 
         endif
      endif
***   define semiconvective and thermohaline zones
      do j = 1,2
         kl = 1
         jl = 2*j+7
         if (j.eq.2) it = 2
         if (crz(1).eq.it) novlim(kl,jl) = 1
         do k = 2,nmod1
            tbot = crz(k).ne.it.and.crz(k+1).eq.it
            ttop = crz(k).eq.it.and.crz(k+1).ne.it
            if (tbot) novlim(kl,jl) = k+1
            if (ttop) then
               if (novlim(kl,jl).lt.k) then
                  novlim(kl,jl+1) = k
                  kl = kl+1
                  if (kl.gt.nmaxconv) goto 40
               endif
            endif
         enddo
 40      if (j.eq.1) nsemiconv = kl-1
         if (j.eq.2) nthermha = kl-1
      enddo

      if (ledouxmlt) then
         nsconv = nledoux
         novlim(1:nsconv,3) = novlim(1:nsconv,5)
         novlim(1:nsconv,4) = novlim(1:nsconv,6)
      else
         nsconv = nschwar
         novlim(1:nsconv,3) = novlim(1:nsconv,1)
         novlim(1:nsconv,4) = novlim(1:nsconv,2)
      endif
      scr = 'r'
      if (nsconv.gt.0) scr = 'c'

*** define phase
      teff = exp(lnt(neff))
      rgbphase = nphase.eq.3.and.teff.lt.teffagb.and.totm.le.2.5d0
      agbphase = nphase.eq.5.and.teff.lt.teffagb.and.totm.le.10.d0
      superagb = nphase.eq.5.and.teff.lt.teffagb.and.totm.ge.7.d0.and.
     &     totm.le.12.d0.and.xsp(1,ic12).lt.2d-2

***   initialize variables for accretion
      dmaccr = 0.d0
      if (iaccr.ge.1) then
         forall (l = 1:nsp, k = 1:nmod) vxsp(k,l) = xsp(k,l)
         if (numeric.eq.2.or.numeric.eq.3) then
c.. first order accuracy in the spatial derivatives
            forall (i = 1:nmod)        
               wi(i) = 0.5d0
               wj(i) = 0.5d0
            end forall
         else
c.. second order accuracy in the spatial deriatives
            forall (i = 2:nmod1)
               wi(i) = dm(i)/(dm(i)+dm(i-1))
               wj(i) = 1.d0-wi(i)
            end forall
            wi(1) = 0.5d0
            wj(1) = 0.5d0
            wi(nmod) = 0.5d0
            wj(nmod) = 0.5d0
         endif
         ishockb = nmod1
         ishockt = 1
         sigma0 = xiaccr*dsqrt(g*r(nmod)*m(nmod))
         iacctop = 1
         iaccbot = nmod
         do k = 1,nmod
            if (facc(k).gt.0.d0.and.iaccbot.eq.nmod) iaccbot = k
            if (facc(k).gt.0.d0) iacctop = max(iacctop,k)
         enddo
         forall (i = 1:nmod)
            facc(i) = vfacc(i)
            omi(i) = 0.d0
            psi(i) = 0.d0
            ft(i) = 1.d0
            fp(i) = 1.d0
         end forall
         hbb = .false.
         klenv = nsconv
         if (accphase.eq.0) call thermo (error)
      endif
      
      return
      end
      SUBROUTINE rininet

************************************************************************
* Read nuclear network and isotopes properties                         *
* 09/03: All nuclei with mass fractions greater than 1d-15 diffuse     *
*                                                                      *
* $LastChangedDate:: 2016-05-17 18:06:34 +0200 (Mar, 17 mai 2016)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 83                                                          $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.data'
      include 'evolcom.nuc'
      include 'evolcom.spec'

      integer normx,na,nb,nc,nd0,ia,ib,ic,id,iznuc,ianuc
      integer idum,kbeta,kbetaec,kneut
      integer ispec(0:nbz,2*nbz+3)
      integer i,j,k,l

      double precision sum1

      logical diffelem

      common /logicdiff/ diffelem(nsp)

      integer nz,lirn,ina,jna
      parameter (nz = 370,lirn=nsp)

c      double precision v(nreac),y(nsp),eqd(nsp,3*nsp+1)
c      common /nuc_matrix/ jna(nz),ina(lirn)

*___________________________________________
***   determination of numerical constraints
*-------------------------------------------

      ymin = 1.d-50
      ytmin = 1.d-12*tolnuc
      dnuc = 1.d-10*tolnuc

*_________________________________________________
***   determination of the nuclear matrix features
*-------------------------------------------------

      fact(0) = 1.d0
      fact(1) = 1.d0
      fact(2) = 2.d0
      fact(3) = 6.d0
      fact(4) = 24.d0
      fact(5) = 120.d0
      fact(6) = 720.d0

      do l = 1,2*nbz+3
         do j = 0,nbz
            ispec(j,l) = -1
         enddo
      enddo

*__________________
***   read nuclides
*------------------

      ial26g = -1
      ial26m = -1
      ndecay = 0
      do i = 1,nsp
         read (99,1000) idum,elem(i),anuc(i),znuc(i),xspacc(i),istate(i)
         if (istate(i).eq.'F'.and.i.lt.nis) then
            ndecay = ndecay+1
            kdecay(ndecay) = i
         endif
         if (i.lt.nsp) ispec(int(znuc(i)),int(anuc(i))) = i
         if (elem(i).eq.'AL26M'.or.elem(i).eq.'AL26m') ial26m = i
         if (elem(i).eq.'AL26G'.or.elem(i).eq.'AL26g') ial26g = i
      enddo
      if (ial26m.eq.-1) then
          write (nout,*) 'ial26m not defined'
         stop 'STOP rininet'
      endif
      if (ial26g.eq.-1) then
          write (nout,*) 'ial26g not defined'
         stop 'STOP rininet'
      endif
      if (ndecay.ne.10) then
         write (nout,*) ndecay
         stop 'rininet : wrong starevol.par, check istate'
      endif

*_______________________________________________________
***   initialize molecule index for opacity calculations
*-------------------------------------------------------

      idxmol(1,1) = iih2
      idxmol(1,2) = iih
      idxmol(1,3) = iih
      idxmol(2,1) = iih2o
      idxmol(2,2) = iih
      idxmol(2,3) = iio
      idxmol(3,1) = iioh
      idxmol(3,2) = iih
      idxmol(3,3) = iio
      idxmol(4,1) = iico
      idxmol(4,2) = iic
      idxmol(4,3) = iio
      idxmol(5,1) = iicn
      idxmol(5,2) = iic
      idxmol(5,3) = iin
      idxmol(6,1) = iic2
      idxmol(6,2) = iic
      idxmol(6,3) = iic
      idxmol(7,1) = iin2
      idxmol(7,2) = iin
      idxmol(7,3) = iin

*_____________________________
***   initialize species index
*-----------------------------

      ineutron = ispec(0,1)
      ih1 = ispec(1,1)
      ih2 = ispec(1,2)
      ihe3 = ispec(2,3)
      ihe4 = ispec(2,4)
      ili6 = ispec(3,6)
      ili7 = ispec(3,7)
      ibe7 = ispec(4,7)
      ibe9 = ispec(4,9)
      ib8 = ispec(5,8)
      ib11 = ispec(5,11)
      ic12 = ispec(6,12)
      ic13 = ispec(6,13)
      ic14 = ispec(6,14)
      in13 = ispec(7,13)
      in14 = ispec(7,14)
      in15 = ispec(7,15)
      io15 = ispec(8,15)
      io16 = ispec(8,16)
      io17 = ispec(8,17)
      io18 = ispec(8,18)
      if18 = ispec(9,18)
      if19 = ispec(9,19)
      if20 = ispec(9,20)
      ine20 = ispec(10,20)
      ine23 = ispec(10,23)
      ina22 = ispec(11,22)
      ina23 = ispec(11,23)
      ina24 = ispec(11,24)
      ina25 = ispec(11,25)
      img24 = ispec(12,24)
      img25 = ispec(12,25)
      img27 = ispec(12,27)
      ial27 = ispec(13,27)
      isi28 = ispec(14,28)
      ip31 = ispec(15,31)
      is32 = ispec(16,32)
      is35 = ispec(16,35)
      icl35 = ispec(17,35)
      icl36 = ispec(17,36)
c     Ajout pour diffusion Thoul                                 ! modif Thoul
      ine21 = ispec(10,21)
      ine22 = ispec(10,22)
      img26 = ispec(12,26)
      ib10 = ispec(5,10)
      
      zsol = 1.d0-xspsol(ih1)-xspsol(ih2)-xspsol(ihe3)-xspsol(ihe4)

*______________________________________________________
***   initialize desintegration half-lifes (in seconds)
*------------------------------------------------------

      do i = 1,nsp
         tdecay(i) = 1.d99
      enddo
      tdecay(ibe7)   = 4650048.d0
      tdecay(ib8)    = 0.770d0
      tdecay(in13)   = 598.2d0
      tdecay(ic14)   = 1.803478d11
      tdecay(io15)   = 122.2d0
      tdecay(if18)   = 6588.d0
      tdecay(if20)   = 11.d0
      tdecay(ina22)  = 82205792.d0
      tdecay(ine23)  = 37.2d0
      tdecay(ina24)  = 5.385d4
      tdecay(ina25)  = 59.3d0
      tdecay(ial26m) = 6.3452d0
      tdecay(ial26g) = 2.272098d13
      tdecay(img27)  = 567.d0
      tdecay(is35)   = 7534080.d0
      tdecay(icl36)  = 9.498634d14


*___________________________________________________
***   initialize index idxspc and idxsps for outputs
*---------------------------------------------------

      i = 0
      do k = 1,nsp-2
         iznuc = int(znuc(k))
         ianuc = int(anuc(k))
         if (.not.((iznuc.eq.5.and.ianuc.eq.8).or.k.eq.in13.or.
     &        (iznuc.eq.9.and.(ianuc.eq.18.or.ianuc.eq.20)).or.
     &        (iznuc.eq.10.and.ianuc.eq.23).or.(iznuc.eq.11.and.
     &        (ianuc.eq.22.or.ianuc.eq.24.or.ianuc.eq.25)).or.
     &        (iznuc.eq.12.and.ianuc.eq.27))) then
            i = i+1
            if (i.gt.nprint) then
               write (nout,*) 'nprint (',nprint,') must be set to',i,
     &              ' in evolpar.star'
               stop ' rininet : idxsp[c,s] dimension problem'
            endif
            idxspc(i) = k
            idxsps(i) = k
         endif
      enddo
c..   index for .c5 file
      idxspci(1) = in13
      idxspci(2) = if18
      idxspci(3) = if20
      idxspci(4) = ine23
      idxspci(5) = ina22
      idxspci(6) = ina24
      idxspci(7) = ina25
      idxspci(8) = img25
      idxspci(9) = img27

      do i = 1,nsp
         diffelem(i) = .true.
      enddo

*_______________________________________
***   accreted abundance renormalization
*---------------------------------------

      if (iaccr.ge.1) then
         normx = 1
         sum1 = 0.d0
         do l = 1,nsp
            if (xspacc(l).gt.xspacc(normx)) normx = l
            sum1 = sum1+xspacc(l)
         enddo
         if (abs(sum1-1.d0).gt.1.d-6) then
            write (nout,*) 'composition of accreted matter not ',
     &           'normalized : sumxsp-1 = ',sum1-1.d0
            stop 'rininet'
         endif
         xspacc(normx) = xspacc(normx)+(1.d0-sum1)

         muiacc = 0.d0
         do l = 1,nis
            muiacc = muiacc+xspacc(l)/anuc(l)
         enddo
      endif

      read (99,1200)

*________________________________________________________________
***   read nuclear reaction network
*  reaction : kna X(kia) + knb X(kib) --> knd X(kid) + knc X(kic)
***      ex :  1    C12  +  1    He   -->  0  gamma  +  1   O16
*----------------------------------------------------------------

      iddn = -1
      ippd = -1
      inpg = -1
      i3a = -1
      icag = -1
      ioag = -1
      iccga = -1
      iccgp = -1
      ioca = -1
      iocp = -1
      iocn = -1
      iooga = -1
      ioogp = -1
      inega = -1
      ippg = -1
      ipdg = -1
      i2he3 = -1
      ihe3ag = -1
      ili6pa = -1
      ili7pa = -1
      ili7d = -1
      ibe7beta = -1
      ibe7pg = -1
      ib8beta = -1
      ib11pa = -1
      icpg = -1
      ic13pg = -1
      in13beta = -1
      in14pg = -1
      in15pa = -1
      in15pg = -1
      io15beta = -1
      io16pg = -1
      io17pa = -1
      ine20ag = -1
      img24ag = -1
      isi28ag = -1

      kbeta = 0
      kbetaec = 0
      kneut = 0

      do l = 1,nreac
         read (99,1100) idum,na,nb,nd0,nc,qi(l),ia,ib,id,ic,qinu(l)
***   index of beta decay reactions
         if (ib.eq.0.and.id.eq.0) then
            kbeta = kbeta+1
            if (kbeta.gt.nbeta) goto 10
            jbeta(kbeta) = l
            if (l.gt.nre) kbetaec = kbetaec+1
         endif
***   index of neutron capture reactions
         if (ib.eq.1) then
            kneut = kneut+1
            if (kneut.gt.nneut) goto 10
            jneut(kneut) = l
         endif

         if (ib.eq.0) ib = nsp
         if (id.eq.0) id = nsp
         if (ia.eq.ih2.and.na.eq.2.and.ic.eq.ihe3) iddn = l
         if (ia.eq.ili7.and.ib.eq.ih2.and.ic.eq.ihe4) ili7d = l
         if (ia.eq.ih1.and.na.eq.2.and.ic.eq.ih2) ippd = l
         if (ia.eq.in14.and.ib.eq.ih1.and.ic.eq.io15) inpg = l
         if (ia.eq.ic12.and.ib.eq.ih1.and.ic.eq.in13) icpg = l
         if (ia.eq.ihe4.and.ic.eq.ic12) i3a = l
         if (ia.eq.ic12.and.ib.eq.ihe4.and.ic.eq.io16) icag = l
         if (ia.eq.io16.and.ib.eq.ihe4.and.ic.eq.ine20) ioag = l
         if (ia.eq.ic12.and.id.eq.ihe4.and.ic.eq.ine20) iccga = l
         if (ia.eq.ic12.and.id.eq.ih1.and.ic.eq.ina23) iccgp = l
         if (ia.eq.io16.and.id.eq.ihe4.and.ic.eq.isi28) iooga = l
         if (ia.eq.io16.and.id.eq.ih1.and.ic.eq.ip31) ioogp = l
         if (ia.eq.io16.and.id.eq.ihe4.and.ic.eq.img24) ioca = l
         if (ia.eq.io16.and.id.eq.ih1.and.ic.eq.ial27) iocp = l
         if (ia.eq.io16.and.id.eq.1.and.ic.eq.ial27) iocn = l
         if (ia.eq.ine20.and.id.eq.ihe4.and.ic.eq.io16) inega = l
         if (ia.eq.ih1.and.na.eq.2.and.ic.eq.ih2) ippg = l
         if (ia.eq.ih2.and.ib.eq.ih1.and.ic.eq.ihe3) ipdg = l
         if (ia.eq.ihe3.and.na.eq.2.and.ic.eq.ihe4) i2he3 = l
         if (ia.eq.ihe4.and.ib.eq.ihe3.and.ic.eq.ibe7) ihe3ag = l
         if (ia.eq.ili6.and.ib.eq.ih1.and.ic.eq.ihe4) ili6pa = l
         if (ia.eq.ili7.and.ib.eq.ih1.and.ic.eq.ihe4) ili7pa = l
         if (ia.eq.ibe7.and.ib.eq.nsp.and.ic.eq.ili7) ibe7beta = l
         if (ia.eq.ibe7.and.ib.eq.ih1.and.ic.eq.ib8) ibe7pg = l
         if (ia.eq.ib8.and.ib.eq.nsp.and.ic.eq.ihe4) ib8beta = l
         if (ia.eq.ib11.and.ib.eq.ih1.and.ic.eq.ihe4) ib11pa = l
         if (ia.eq.ic13.and.ib.eq.ih1.and.ic.eq.in14) ic13pg = l
         if (ia.eq.in13.and.ib.eq.nsp.and.ic.eq.ic13) in13beta = l
         if (ia.eq.in14.and.ib.eq.ih1.and.ic.eq.io15) in14pg = l
         if (ia.eq.in15.and.ib.eq.ih1.and.ic.eq.ic12) in15pa = l
         if (ia.eq.in15.and.ib.eq.ih1.and.ic.eq.io16) in15pg = l
         if (ia.eq.io15.and.ib.eq.nsp.and.ic.eq.in15) io15beta = l
         if (ia.eq.io16.and.ib.eq.ih1.and.ic.eq.io17) io16pg = l
         if (ia.eq.io17.and.ib.eq.ih1.and.id.eq.ihe4) io17pa = l
         if (ia.eq.ine20.and.ib.eq.ihe4.and.ic.eq.img24) ine20ag = l
         if (ia.eq.img24.and.ib.eq.ihe4.and.ic.eq.isi28) img24ag = l
         if (ia.eq.isi28.and.ib.eq.ihe4.and.ic.eq.is32) isi28ag = l

         k1(l) = ia
         k2(l) = na
         if (na.eq.0) then
            write (nout,*) ' Network error for reaction : ',idum
            write (nout,*) ' Wrong stoechiometric coefficients'
            stop 'rininet'
         endif
         k3(l) = ib
         k4(l) = nb
         k5(l) = ic
         k6(l) = nc
         k7(l) = id
         k8(l) = nd0
      enddo
      qdeut = qi(2)*econv
      if (kneut+1.ne.idneut) then
         write (nout,*) ' idneut = ',idneut,' not equal to kneut =',
     &        kneut+1
         write (nout,*) ' set idneut = ',kneut+1,' in evolpar.star'
         stop 'STOP rininet'
      endif
      if (ippd.eq.-1) then
         write (nout,*) 'ippd not defined'
         stop 'STOP rininet'
      endif
      if (iddn.eq.-1) then
         write (nout,*) 'iddn not defined'
         stop 'STOP rininet'
      endif
      if (inpg.eq.-1) then
         write (nout,*) 'inpg not defined'
         stop 'STOP rininet'
      endif
      if (i3a.eq.-1) then
         write (nout,*) 'i3a not defined'
         stop 'STOP rininet'
      endif
      if (icag.eq.-1) then
         write (nout,*) 'icag not defined'
         stop 'STOP rininet'
      endif
      if (ioag.eq.-1) then
         write (nout,*) 'ioag not defined'
         stop 'STOP rininet'
      endif
      if (iccga.eq.-1) then
         write (nout,*) 'iccga not defined'
         stop 'STOP rininet'
      endif
      if (iccgp.eq.-1) then
         write (nout,*) 'iccgp not defined'
         stop 'STOP rininet'
      endif
      if (iocp.eq.-1) then
         write (nout,*) 'iocp not defined'
         stop 'STOP rininet'
      endif
      if (ioca.eq.-1) then
         write (nout,*) 'ioca not defined'
         stop 'STOP rininet'
      endif
      if (iocn.eq.-1) then
         write (nout,*) 'iocn not defined'
         stop 'STOP rininet'
      endif
      if (inega.eq.-1) then
         write (nout,*) 'inega not defined'
         stop 'STOP rininet'
      endif
      if (ippg.eq.-1) then
         write (nout,*) 'ippg not defined'
         stop 'STOP rininet'
      endif
      if (ipdg.eq.-1) then
         write (nout,*) 'ipdg not defined'
         stop 'STOP rininet'
      endif
      if (i2he3.eq.-1) then
         write (nout,*) 'i2he3 not defined'
         stop 'STOP rininet'
      endif
      if (ihe3ag.eq.-1) then
         write (nout,*) 'ihe3ag not defined'
         stop 'STOP rininet'
      endif
      if (ili7d.eq.-1) then
         write (nout,*) 'ili7d not defined'
         stop 'STOP rininet'
      endif
      if (ili6pa.eq.-1) then
         write (nout,*) 'ili6pa not defined'
         stop 'STOP rininet'
      endif
      if (ili7pa.eq.-1) then
         write (nout,*) 'ili7pa not defined'
         stop 'STOP rininet'
      endif
      if (ibe7beta.eq.-1) then
         write (nout,*) 'ibe7beta not defined'
         stop 'STOP rininet'
      endif
      if (ibe7pg.eq.-1) then
         write (nout,*) 'ibe7pg not defined'
         stop 'STOP rininet'
      endif
      if (ib8beta.eq.-1) then
         write (nout,*) 'ib8beta not defined'
         stop 'STOP rininet'
      endif
      if (ib11pa.eq.-1) then
         write (nout,*) 'ib11pa not defined'
         stop 'STOP rininet'
      endif
      if (icpg.eq.-1) then
         write (nout,*) 'icpg not defined'
         stop 'STOP rininet'
      endif
      if (ic13pg.eq.-1) then
         write (nout,*) 'ic13pg not defined'
         stop 'STOP rininet'
      endif
      if (in13beta.eq.-1) then
         write (nout,*) 'in13beta not defined'
         stop 'STOP rininet'
      endif
      if (in14pg.eq.-1) then
         write (nout,*) 'in14pg not defined'
         stop 'STOP rininet'
      endif
      if (in15pa.eq.-1) then
         write (nout,*) 'in15pa not defined'
         stop 'STOP rininet'
      endif
      if (in15pg.eq.-1) then
         write (nout,*) 'in15pg not defined'
         stop 'STOP rininet'
      endif
      if (io15beta.eq.-1) then
         write (nout,*) 'io15beta not defined'
         stop 'STOP rininet'
      endif
      if (io16pg.eq.-1) then
         write (nout,*) 'ippd not defined'
         stop 'STOP rininet'
      endif
      if (io17pa.eq.-1) then
         write (nout,*) 'io17pa not defined'
         stop 'STOP rininet'
      endif
      if (ine20ag.eq.-1) then
         write (nout,*) 'ine20ag not defined'
         stop 'STOP rininet'
      endif
      if (img24ag.eq.-1) then
         write (nout,*) 'img24ag not defined'
         stop 'STOP rininet'
      endif
      if (isi28ag.eq.-1) then
         write (nout,*) 'isi28ag not defined'
         stop 'STOP rininet'
      endif

*________________________________
***   read nuclear reactions name
*--------------------------------

      rewind (99)
      do i = 1,33+nsp+1
         read (99,1200)
      enddo
      do j = 1,nreac
         read (99,1300) react(j)
      enddo

      close (99)  !  starevol.par

*______________________________________
***   initialize nuclear reaction types
*--------------------------------------

      do l = 1,nreac
***   index of neutron emission reactions
         if (k7(l).eq.1) then
            kneut = kneut+1
            jneut(kneut) = l
         endif
      enddo

***   LiBeB reactions
      jlibeb(1) = ili7pa
      jlibeb(2) = ib11pa
      jlibeb(3) = ili6pa

***   check dimension compatibility
 10   if (kneut.ne.nneut.or.kbeta.ne.nbeta.or.kbetaec.ne.nbetaec) then
         write (nout,*) 'Wrong dimension for vectors jneut/jbeta'
         if (kneut.gt.nneut) write (nout,*) 'The number of detected ',
     &        'reactions involving neutrons is kneut =',kneut,', it',
     &        ' must equal nneut =',nneut,' as defined in evolpar.star'
         if (kbeta.gt.nbeta) write (nout,*) 'The number of detected ',
     &        'beta decay reactions is kbeta =',kbeta,', it must equal',
     &        ' nbeta =',nbeta,' as defined in evolpar.star'
         if (kbetaec.gt.nbetaec) write (nout,*) 'The number of ',
     &        'detected electron capture reactions is kbetaec =',
     &        kbetaec,', it must equal nbetaec =',nbetaec,' as defined',
     &        ' in evolpar.star'
         write (nout,*) 'Check starevol.par (version >= 2.60 - 2.80)'
         stop 'STOP in rininet'
      endif


*_________________________________________________
***   initialize pointer for sparse nuclear matrix
*-------------------------------------------------

C$$$      v(1:nreac) = 1.d2
C$$$      y(1:nsp) = 1.d-2
C$$$      call freac (v,eqd,y,1.d-2,.true.)
C$$$      na = 0
C$$$      ia = 1
C$$$      ina(ia) = 1
C$$$      do  i = 1, nis
C$$$        do  j = 1, nis
C$$$           if (eqd(i,j) .ne. 0.d0) then
C$$$              na = na+1
C$$$              jna(na) = j
C$$$c              print *,na,j
C$$$           endif
C$$$        enddo
C$$$        ia = ia+1
C$$$        ina(ia) = na+1
C$$$      enddo
C$$$c      print *,ina

 1000 format (1x,i3,1x,a5,1x,f9.5,1x,f4.0,3x,1pe11.5,2x,a1)
 1100 format (i4,3x,i2,8x,i2,7x,i2,8x,i2,13x,f7.3,7x,4(i4),8x,f7.3)
 1200 format (1x)
 1300 format (7x,a37)

      return
      end
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

*__________________________________________
***   read values for adjustable parameters
*------------------------------------------

      read (99,900)  evolpar_version
      read (99,1000) maxmod,imodpr,mtini,zkint
      read (99,1100) addH2,addHm,tmaxioH,tmaxioHe
      read (99,1200) ionizHe,ionizC,ionizN,ionizO,ionizNe
      read (99,1300) lgrav,lnucl
      read (99,1400) alphac,dtnmix,mixopt,hpmlt,nczm
      read (99,1500) ihro,iturb,etaturby,alphatu
      read (99,1600) novopt,aup,adwn,idup
      read (99,1700) tau0,ntprof,taulim
      read (99,1800) mlp0,etapar,dmlinc
      read (99,1900) nucreg,nuclopt,tolnuc,ftnuc,znetmax
      read (99,2000) idiffcc,idiffnuc,idiffty,diffzc,diffst
      read (99,2100) zgradmu,grad_crit,del_zc,ledouxmlt
      read (99,2110) omegaconv,Dh_prescr,dm_ext
      read (99,2120) Dv_prescr,om_sat,D_zc,disctime,eq_thermique
      read (99,2130) diffvr,idiffvr,idiffex,diffcst,breaktime
      read (99,2200) nmixd,itmind,itmaxd,rprecd,nretry
      read (99,2300) lover,lthal,ltach,lmicro
      read (99,2400) iaccr,accphase,itacc,massrate0,massend0,menv
      read (99,2500) ric,prrc,xiaccr,fdtacc,alphaccr
      read (99,2600) chydro,ivisc,q0,mcut
      read (99,2700) maxsh0,nresconv
      read (99,2800) dlnvma0,dlnvmi0,dlnenuc0,dmrma,dmrmi
      read (99,2900) dtiny,dtminy,dtmaxy,facdt,fkhdt
      read (99,3000) ishtest,fts,ftsh,ftshe,ftsc,ftacc,ftst
      read (99,3100) itmin,itermax,itermix,icrash0,numeric,icorr
      read (99,3200) (eps(l),l = 1,5),epsa
      read (99,3300) phi,alpha,iacc,alphmin,alphmax,sigma,vlconv

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

     
***   check versions compatibility
      if (code_version.lt.evolpar_version.and.evolpar_version.lt.2.1d0)
     &     then
         write (nout,*) 'starevol.par not compatible with executable :',
     &        ' network problem ',code_version,evolpar_version
         stop ' stop rmodpar'
      endif

***   check parameters
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
      if (hpmlt.eq.5.and.alphatu.lt.alphac) then
         write (nout,*) 'when HMLT = 5, alphatu MUST be > alphac !'
         stop 'rmodpar : bad alphatu and alphac'
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


      if (mlp.eq.15.or.mlp.eq.16) then
         print *,'read namelist massloss.par'
         open(UNIT=79, FILE='massloss.par', STATUS='OLD')
         rewind(79)
         read(79, NML=MassiveMloss) 
         close(79)
         print *,'clumpfac',clumpfac
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
      write (nout,*)
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

      if (rotation) eq_thermique = .true.
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
      print *,'igw,igwsurfchim,igwrot,igwsurfenerg, = ',igw,igwsurfchim,
     &         igwrot,igwsurfenerg
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
      write (90,4600) alphac,dtnmix,mixopt,hpmlt,nczm
      write (90,4700) ihro,iturb,etaturby,alphatu
      write (90,4800) novopt,aup,adwn,idup
      write (90,4900) tau0,ntprof,taulim
      write (90,5000) mlp0,etapar,dmlinc
      write (90,5100) nucreg,nuclopt,tolnuc,ftnuc,znetmax
      write (90,5200) idiffcc,idiffnuc,idiffty,diffzc,diffst
      write (90,5300) zgradmu,grad_crit,del_zc,ledouxmlt
      write (90,5310) omegaconv,Dh_prescr,dm_ext
      write (90,5320) Dv_prescr,om_sat,D_zc,disctime,eq_thermique
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
 4600 format (1x,'convection    : alphac = ',0pf6.4,', dtnmix = ',a1,
     &     ', mixopt = ',l1,', hpmlt = ',i1,', nczm =',i3)
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
     &     1pd8.2,', disctime = ',1pd8.2,', eq_thermique = ',l1)
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


      return
      end


************************************************************************

      DOUBLE PRECISION FUNCTION screen (y,zbar,ztilde,zeta,l0,h0,z1,z2)

************************************************************************
* Calculate the screening factor for the thermonuclear reaction rates  *
* Input : nuclide charges z1 and z2                                    *
* Graboske etal 1973,ApJ,181,457 and Cox & Guili 17.15 for definitions *
* For intermediate screening :                                         *
* Mitler 1977, ApJ 212, 513 and Cox & Guili 17.15 for definitions      *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.nuc'

      double precision z1,z2,z12
      double precision zetafac,zeta1,zeta2
      double precision ztilde,zbar,h120,l0,l12,h120a
      double precision zeta,h0
      double precision y,zx,zmean,l086

      dimension y(nsp)

      h120 = 0.d0
      l12 = l0*z1*z2*ztilde

c... Graboske et al
c *___________________
c ***   weak screening
c *-------------------
c       if (l12.lt.0.1d0) then
c          h120 = l12

c *___________________________
c ***   intermediate screening
c *---------------------------
c       elseif (l12.ge.0.1d0.and.l12.le.5.d0) then
c          zmean = 0.d0
c          do i = 1,nis
c             zmean = zmean+(znuc(i)**1.58d0*y(i))
c          enddo
c          zmean = zmean/zx
c          l086 = 0.38d0*l0**0.86d0
c          h120 = zmean*ztilde**(-0.58d0)*zbar**(-0.28d0)
c          h120 = h120*((z1+z2)**1.86d0-z1**1.86d0-z2**1.86d0)*l086
c       endif


*____________________________________
***   weak and intermediate screening
***   Mitler (1977, ApJ, 212, 513) 
*------------------------------------
      if (l12.le.5.d0) then
         zeta1 = zeta*z1
         zeta2 = zeta*z2
         zetafac = (zeta1+zeta2+1.d0)**pw53-(zeta1+1.d0)**pw53-(zeta2+
     &        1.d0)**pw53+1.d0
         h120 = abs(h0*zetafac)
      endif

*_____________________
***   strong screening
*---------------------
      if (l12.ge.2.d0) then
         z12 = z1+z2
         h120a = (z12**pw53-z1**pw53-z2**pw53)
     &        +.316d0*zbar**pw13*(z12**pw43-z1**pw43-z2**pw43)
     &        +0.737d0/(zbar*l0**pw23)*(z12**pw23-z1**pw23-z2**pw23)
         h120a = h120a*0.624d0*zbar**pw13*l0**pw23
         if (l12.le.5.d0) h120 = min(h120,h120a)
         if (l12.gt.5.d0) h120 = h120a
         h120 = min(h120,h120a)
      endif
      screen = exp(min(220.d0,h120))

      return
      end
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
***   corrections to the radiative gradient in atmospheric layers
***   when plan-parallel, grey and Eddington approximations are ruled out
*------------------------------------------------------------------------

         xsitau(ip) = 0.d0

         if (tau(ip).lt.taulim.and.(ntprof.eq.1.or.ntprof.ge.3.and.
     &        .not.ntp2)) then
***   plan-parallel corrections
            dptau(ip) = -kapm(ip)*lum(ip)*dqdtau(ip)/(pim4*c*r(ip)**2)
            xsitau(ip) = -dptau(ip)*r(ip)**2/gmrr2
c$$$            xsitau(ip) = 0.d0   ! TEST 31/03/2020
            abrad(ip) = abrad(ip)*(1.d0+dqdtau(ip))/(1.d0+xsitau(ip))
            abla(ip) = abrad(ip)
            
c            write(555,'(2(1x,i4),7(1x,1pe11.4))'),i,crz(i)
c     $        ,dqdtau(i),ddqtau(i),xsitau(i),abrad(i),tau(i),t(i)   
***   spherical corrections
c            xsitau(ip) = lum(ip)/(pim2*gmrr2*c)*(kapm(ip)*dqdtau(ip)*
c     &           0.5d0+(tau(ip)+qtau(ip))/(rom(ip)*r(ip)))
c            abrad(ip) = abrad(ip)*(1.d0+dqdtau(ip)+2.d0*(tau(ip)+
c     &           qtau(ip))/(kapm(ip)*rom(ip)*r(ip)))/(1.d0+xsitau(ip))
***   P/T corrections
c            cortau1 = ((tau(ip)+qtau(ip))/(tau(ip)+pw23))**0.25d0
c            cortau2 = lum(ip)*(qtau(ip)-pw23)/(pim4*c*r(ip)*r(ip)*
c     &           pm(ip))
c            cortau2 = 1.d0/(1.d0-cortau2)
c            abrad(ip) = abrad(ip)*cortau2/cortau1
         endif
         abdf1(ip) = abrad(ip)*wi(ip)*(kiinv*dkapdf(i)+piinv*dpdf(i))
         abdf2(ip) = abrad(ip)*wj(ip)*(kipinv*dkapdf(ip)+pipinv*
     &        dpdf(ip))
         abdt1(ip) = abrad(ip)*wi(ip)*(kiinv*dkapdt(i)+piinv*dpdt(i)-
     &        4.d0)
         abdt2(ip) = abrad(ip)*wj(ip)*(kipinv*dkapdt(ip)+pipinv*dpdt(ip)
     &        -4.d0)
         abdl(ip) = abrad(ip)/lum(ip)
         if (hydro) then
            abdr(ip) = -abrad(ip)/gmrr2*2.d0*r(ip)*accel(ip)
            abdu1(ip) = -abrad(ip)/gmrr2*r(ip)*r(ip)*dtipsi
            abdu2(ip) = -abrad(ip)/gmrr2*r(ip)*r(ip)*(dtninv-dtipsi)
         else
            abdr(ip) = 0.d0
            abdu1(ip) = 0.d0
            abdu2(ip) = 0.d0
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
         if (ntprof.eq.4.and..not.ntp2) then
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
      if (ntprof.eq.1.or.(ntprof.ge.3.and..not.ntp2)) then
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
      if (ntprof.eq.1.or.ntprof.ge.3.and..not.ntp2) then
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
         eq(1,7) = eq(1,7)-2.d0*dptau(nmod)
         eq(1,8) = eq(1,8)+dptau(nmod)*dkapdf(nmod)*kiinv
         eq(1,9) = eq(1,9)+dptau(nmod)*dkapdt(nmod)*kiinv
         eq(1,10) = eq(1,10)+dptau(nmod)/lum(nmod)
      endif

***   lagrangian velocity

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

      elseif (ntprof.eq.4) then
         eq(3,16) = ro(nmod)-rosurf
         eq(3,8) = drodf(nmod)
         eq(3,9) = drodt(nmod)
      endif

***   energy conservation

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

!         if (ntprof.ne.4) then
            eq(5,16) = t(nmod)-tsurfr*tsf
!         else
!!            eq(5,16) = t(nmod)-tatm(nmod)
!         endif
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
      SUBROUTINE thermo (error)

************************************************************************
* Compute thermodynamical quantities and optical depth                 *
*                                                                      *
* $LastChangedDate:: 2016-05-11 17:20:46 +0200 (Mer, 11 mai 2016)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 62                                                          $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.eng'
      include 'evolcom.eos'
      include 'evolcom.ion'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.opa'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer error
      integer neff1
      integer j,k,l,iopa,iopaf

      integer i1,i2,i3,i4

      integer ktmax,ksh
      double precision tmax

      double precision kap1,kapt1,kapr1,kapx1,kapy1,engnup,engnupdt,
     &     engnupdro,eng,engdt,engdro,rst,rstdf,rstdt,taueffl,taune1l,
     &     taunel,rne1l,rnel,tne1l,tnel,lne1l,lnel,deltnel,deltes,drl
      double precision reffcorr,teffcorr,irrad
      double precision zzmean,amean,xamean
      double precision x(nsp),y(nsp)

*______________________________________________________________
***   calculation of the local metallicity and molar abundances
*--------------------------------------------------------------

      do l = 1,nsp
         do k = 1,nmod
            ysp(k,l) = xsp(k,l)/anuc(l)
         enddo
      enddo

*_____________________________________________________________
***   calculation of the equation of state, and of the related
***   thermodynamic quantities (+ derivatives)
*-------------------------------------------------------------
      
      call eos (1,error)
      if (error.gt.0) return
      
*_________________________________________
***   calculation of the screening factors
*-----------------------------------------

      if (iter.le.1) call fscreen (t,ro,xsp,rhpsi,anuc,znuc)

*_________________________________________________________________
***   calculation of H2 equilibrium abundance in case of accretion
*-----------------------------------------------------------------

      if (dmaccr.gt.0.d0.and.xspacc(ih2).gt.1.d-10.and.accphase.eq.0)
     &     call ydequil

*_______________________________
***   local variable definitions
*-------------------------------

      ktmax = 0
      tmax = -1.d0
      do k = 1,nmod
         if (t(k).gt.tmax) then
            ktmax = k
            tmax = t(ktmax)
         endif
      enddo

      iopa = 0
      iopaf = 0
      do k = 1,nmod
         ksh = k
         do l = 1,nis-1
            x(l) = max(ymin*anuc(l),xsp(k,l))
            y(l) = ysp(k,l)
         enddo
         x(nis) = 1.d-50
         x(nsp) = 1.d-50
         y(nis) = 1.d-50
         y(nsp) = 1.d-50
         tconv(k) = 1.d99

*_________________________________________________
***   calculation of the total opacity coefficient
*-------------------------------------------------


         call kappa (t(k),ro(k),muiinv(k),x,ksh,kap1,kapr1,kapt1,kapx1
     &        ,kapy1,nretry,iter,iopa,iopaf,error)


         if (error.gt.0) return
         kap(k) = kap1
         dkapdro(k) = kapr1
         dkapdf(k) = kapr1*drodf(k)
         dkapdt(k) = kapt1*t(k)+kapr1*drodt(k)
         dkapdx(k) = kapx1
         dkapdy(k) = kapy1

*_______________________________________________________
***   calculation of the energy loss by plasma neutrinos
*-------------------------------------------------------

         if (lnucl) then
            if (t(k).gt.1.d7) then
               zzmean = 0.d0
               xamean = 0.d0
               do j = 1,nsp
                  zzmean = zzmean+znuc(j)*y(j)
                  xamean = xamean+y(j)
               enddo
               amean = 1.d0/xamean
               zzmean = zzmean*amean
               call nulosc (t(k),ro(k),amean,zzmean,engnup,engnupdt,
     &              engnupdro)
            else
               engnup = 0.d0
               engnupdt = 0.d0
               engnupdro = 0.d0
            endif
            enupla(k) = engnup

*______________________________________________________________
***   calculation of the rate of produced energy by the nuclear
***   reactions and of the losses due to associated neutrinos
*--------------------------------------------------------------

            irrad = 7.68771025912336d0*ysp(k,1)*ro(k)*dtn*dsqrt(t(k))
            if (irrad.le.1.72d0) then
               signt(k) = 43.8372d0*irrad+10.4d0
            else
               if (irrad.lt.2.92d0) then
                  signt(k) = -69.8417d0*irrad+205.93d0
               else
                  signt(k) = 1.99d0
               endif
            endif

            call denucl (t(k),ro(k),mueinv(k),y,ksh,eng,engdt,engdro)

            if (iter.le.1.and.ktmax.gt.0.and.imodpr.eq.11.and.
     &           ireset.eq.0) call myengen(k,ktmax)

            denucldro(k) = engdro-engnupdro
            denucldf(k) = denucldro(k)*drodf(k)
            denucldt(k) = (engdt-engnupdt)*t(k)+denucldro(k)*drodt(k)
            enucl(k) = eng-enupla(k)
         endif
      enddo

*___________________________________________________
***   smooth opacity profile at molecular transition
*---------------------------------------------------

      if (iopa.gt.0) then
         do k = iopa-10,nmod
            kap(k) = log(kap(k))
         enddo
         call lissage (kap,min(nmod1-1,iopa+10),iopa-3,5)
         call lissage (dkapdf,min(nmod1-1,iopa+10),iopa-3,5)
         call lissage (dkapdt,min(nmod1-1,iopa+10),iopa-3,5)
         do k = iopa-10,nmod
            kap(k) = exp(kap(k))
         enddo
      endif

*______________________________________________________________
***   smooth opacity profile at OPAL <--> Fergusson transition
*--------------------------------------------------------------
* This smoothing leads to a step in temperature at this opacity
* transition
* Commented out on 04/02/2014

c$$$      if (iopaf.gt.0) then
c$$$         do k = iopaf-10,nmod
c$$$            kap(k) = log(kap(k))
c$$$         enddo
c$$$         call lissage (kap,min(nmod1-1,iopaf+10),iopaf-3,5)
c$$$         call lissage (dkapdf,min(nmod1-1,iopaf+10),iopaf-3,5)
c$$$         call lissage (dkapdt,min(nmod1-1,iopaf+10),iopaf-3,5)
c$$$         do k = iopaf-10,nmod
c$$$            kap(k) = exp(kap(k))
c$$$         enddo
c$$$      endif

*_____________________________________________
***   calculation of the optical depth profile
*---------------------------------------------
      
      tau(nmod) = tau0
      dtaudf(nmod) = 0.d0
      dtaudt(nmod) = 0.d0
      do l = nmod1,2,-1
         drl = abs(r(l+1)-r(l))
         rst = kap(l)*ro(l)*drl
         rstdf = (kap(l)*drodf(l)+ro(l)*dkapdf(l))*drl
         rstdt = (kap(l)*drodt(l)+ro(l)*dkapdt(l))*drl
         tau(l) = tau(l+1)+rst
         dtaudf(l) = dtaudf(l+1)+rstdf
         dtaudt(l) = dtaudt(l+1)+rstdt
      enddo
      tau(1) = tau(2)
      dtaudf(1) = dtaudf(2)
      dtaudt(1) = dtaudt(2)

*____________________________________________
***   calculation of the effective quantities
*--------------------------------------------

      neff = 0
      do k = nmod-2,2,-1
         if (tau(k).gt.taueff.and.tau(k+1).ge.taueff) exit
      enddo
      neff = k+2
      neff1 = neff-1
      taueffl = log(abs(taueff))
      taune1l = log(abs(tau(neff1)))
      taunel = log(abs(tau(neff)))
      rne1l = lnr(neff1)
      rnel = lnr(neff)
      tne1l = lnt(neff1)
      tnel = lnt(neff)
      lne1l = log(abs(lum(neff1)))
      lnel = log(abs(lum(neff)))
      deltnel = tne1l-tnel
***   interpolate quatities at tau = 2/3
      tphot = tne1l-(taune1l-taueffl)*deltnel/(taune1l-taunel)
      deltes = tne1l-tphot
      tphot = exp(tphot)
      if (deltnel.ne.0.d0) then
         reff = rne1l-deltes*(rne1l-rnel)/deltnel
         leff = lne1l-deltes*(lne1l-lnel)/deltnel
         reff = exp(reff)
         leff = exp(leff)
      else
         leff = lum(neff)
         reff = r(neff)
      endif
      teff = (leff/(pim4*sig*reff*reff))**0.25d0
      
*_____________________________________________________________
***   corrected atmospheric temperature profile for optically
***   thick winds (WR phase)
*-------------------------------------------------------------

c      if (xsp(nmod1,ih1).le.0.2d0.and.log10(teff).gt.4.d0.and.
c     &     m(nmod)/msun.gt.10.d0) then
c         call corrwind (teffcorr,reffcorr)
c         teff = 10.d0**teffcorr
c         reff = reffcorr
c      endif
      
*__________________________________________
***   calculation of the internal structure
*------------------------------------------

      call structure (error)
      
      return
      end


************************************************************************

      SUBROUTINE tstep (ng,dtnt,dtnf,dtnro,dtnr,dtnu,dtnl,idtnt,
     &     idtnf,idtnro,idtnr,idtnu,idtnl)

************************************************************************
* Find the smallest time-step according to the changes of              *
* the temperature and density profiles in the last two models          *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.transp'
      include 'evolcom.var'

      integer idtnt,idtnf,idtnro,idtnr,idtnu,idtnl
      integer k,ng,ni,nf

      double precision vom,vrray,vxpsi
      double precision dtnt,dtnf,dtnro,dtnr,dtnu,dtnl
      double precision dtt,dtf,dtro,dtr,dtu,dtl
      double precision ft

      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)

      dtnacc = 1.d99
      idtnu = 1
      dtnu = 1.d99
      idtnr = 1
      dtnr = 1.d99
      idtnf = 1
      dtnf = 1.d99
      idtnt = 1
      dtnt = 1.d99
      idtnro = 1
      dtnro = 1.d99
      idtnl = 1
      dtnl = 1.d99

      if (dtmin.eq.dtmax0) return

      ft = fts*dtn

      ni = 1
      nf = ng

***   Mass cut : only above Mcut is the time step constrained (hydro)
      if (mcut.lt.m(ng).and.hydro.and.ivisc.gt.0) then
         do k = ng,2,-1
            if (m(k).lt.mcut) goto 12
         enddo
 12      ni = k
      endif

*________________________
***   accretion time-step
*------------------------

      if (iaccr.gt.0.and.massrate.gt.0.d0) then
         if (accphase.eq.6) then
            dtnacc = ftacc*(m(nmod)-menv)*sec/(msun*massrate)
         elseif (accphase.eq.7) then
            dtnacc = ftacc*menv*sec/(msun*massrate)
         elseif (accphase.eq.8) then
            dtnacc = ftacc*m(nmod)*(1.d0-menv)*sec/(msun*massrate)
         else
            dtnacc = ftacc*(m(iacctop)-m(iaccbot))*sec/(msun*massrate)
         endif
         if (dtnacc.lt.0.d0) stop 'problem in mass accretion regim'
      endif

*________________________________________________________
***   timestep constrained in the shell region : ni to nf
*--------------------------------------------------------

      do k = ni+1,nf

         if (hydro.and.u(k).ne.vu(k)) then
            dtu = ft*abs(u(k)/(u(k)-vu(k)))
         else
            dtu = 1.d99
         endif
         dtr = ft*abs(r(k)/(r(k)-vr(k)))
         dtf = ft*abs(lnf(k)/(lnf(k)-vlnf(k)))
         dtt = ft*t(k)/abs(t(k)-vt(k))
         dtro = ft*ro(k)/abs(ro(k)-vro(k))
         dtl = ft*abs(lum(k)/(lum(k)-vlum(k)))

         if (dtu.lt.dtnu) then
            idtnu = k
            dtnu = dtu
         endif
         if (dtr.lt.dtnr) then
            idtnr = k
            dtnr = dtr
         endif
         if (dtf.lt.dtnf) then
            idtnf = k
            dtnf = dtf
         endif
         if (dtt.lt.dtnt) then
            idtnt = k
            dtnt = dtt
         endif
         if (dtro.lt.dtnro) then
            idtnro = k
            dtnro = dtro
         endif
         if (dtl.lt.dtnl) then
            idtnl = k
            dtnl = dtl
         endif
      enddo

      return
      end


************************************************************************

      SUBROUTINE tstepx (mtini,totm,nmod1,dtnx,idtnx,iburning)

************************************************************************
* Constraint the time-step by the abundances variations                *
* in the last two models                                               *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.diff'
      include 'evolcom.eos'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.conv'
      include 'evolcom.spec'
      include 'evolcom.surf'  ! TD 03/2020
      include 'evolcom.therm'
      include 'evolcom.var'

      integer nmod1,idtnx,iburning,inuctop
      integer ite,iroe,nphase0
      integer itmax,ishockb,ishockt
      integer i,k

      double precision tk,rok,mueinvk
      double precision mtini,totm,dtnx
      double precision ftsh0,ftshe0
      double precision v(nreac)
      double precision rcno
      double precision a1,a2,b1,b2
      double precision t9,t913,t923,vd,vdmean,dtnd,yp,yhe4,yc12,yne20,
     &     yn14,yo16,vpp,dtnH,vn14,v3a,vc12,vo16,vcc,vcca,vccp,voca,
     &     vocp,vocn,vne,vneinv,vHe4,dtHe4,dtcc,dtne,dtoo,dtsi
      double precision TburnC,TburnNe
      double precision Lmax,tmax,Mackmax,enucmax


      common /burning/TburnC,TburnNe
      common /electcap_par/ ite,iroe,a1,a2,b1,b2
      common /hydrodyn/ Lmax,tmax,Mackmax,enucmax,itmax,ishockb,ishockt

      nphase0 = nphase
      iburning = 1
      idtnx = 1
      dtnx = 1.d90
      dtnH = 1.d99
      dtHe4 = 1.d99
      dtcc = 1.d99
      dtne = 1.d99
      dtoo = 1.d99
      dtsi = 1.d99
      if (nphase.eq.2) then
         ftsh0 = min(ftsh,2.d-2)
      else
         ftsh0 = ftsh
      endif
      if (nphase.eq.4.or.thermalpulse) then
         ftshe0 = min(ftshe,5.d-2)
      else
         ftshe0 = ftshe
      endif

*_____________________________________________
***   determination of nuclearly active shells
*---------------------------------------------

      inuctop = nmod1
      do i = 2,nmod1
         if (t(i).lt.tnucmin) exit
      enddo
      inuctop = i-1

*________________________________________________________
***   time-step constrained by D-burning (PMS phase only)
*--------------------------------------------------------

      if (nphase.eq.1.and.xsp(1,ih2).gt.1.d-8.and.crz(1).lt.-1.and.
     &     iaccr.eq.0.and.totm.lt.3.d0.and.inuctop.gt.1) then
         vdmean = 0.d0
         do k = 1,inuctop
            t9 = t(k)*1.d-9
            t913 = t9**pw13
            t923 = t913*t913
            vd = (2.24d3/t923*exp(max(-220.d0,-3.72d0/t913))*(1.d0+
     &           0.112d0*t913+3.38d0*t923+2.65d0*t9))*ro(k)
            vdmean = vdmean+vd*dm(k)
         enddo
         vdmean = vdmean/m(inuctop)
         dtnd = 2.d0/(vdmean*xsp(1,ih1))
         dtnd = max(dtnd,1.d3*sec/mtini,1.d3)
         if (dtnd.le.dtnx) then
            dtnx = dtnd
            iburning = 1
         endif
      endif

*_____________________________
***   variable initializations
*-----------------------------

      do 20 k = 1,inuctop
         tk = t(k)
         rok = ro(k)
         mueinvk = 2.d0/(1.d0+xsp(inuctop,ih1))
         signt(k) = 10.4d0
         call vit (tk,rok,mueinvk,k,v,0)
         yp = xsp(k,ih1)
         yhe4 = xsp(k,ihe4)/anuc(ihe4)
         yc12 = xsp(k,ic12)/anuc(ic12)
         yn14 = xsp(k,in14)/anuc(in14)
         yo16 = xsp(k,io16)/anuc(io16)
         yne20 = xsp(k,ine20)/anuc(ine20)

*_______________________________________
***   time-step constrained by H-burning
*---------------------------------------

         if (tk.lt.5.d7.and.yhe4.gt.2.5d-4) then
***   2 proto ( 0 gamma, 0 nutri)  1 deutr
            vpp = v(ippd)*yp*0.5d0
            dtnH = ftsh0/vpp
            if (dtnH.lt.dtnx) then
               dtnx = dtnH
               idtnx = k
               iburning = 2
            endif
            goto 20
         endif
         if (tk.lt.1.5d8.and.yp.gt.1.d-3) then
***   1 C  12 ( 1 proto, 0 gamma)  1 N  13
            vc12 = v(icpg)*yc12
            dtnH = ftsh0/vc12
            if (dtnH.lt.dtnx) then
               dtnx = dtnH
               idtnx = k
               iburning = 3
            endif
***   1 N  14 ( 1 proto, 0 gamma)  1 O  15
            vn14 = v(inpg)*yn14
            dtnH = ftsh0/vn14
            if (dtnH.lt.dtnx) then
               dtnx = dtnH
               idtnx = k
               iburning = 4
            endif
            goto 20
         endif

*________________________________________
***   time-step constrained by He-burning
*----------------------------------------

         if (tk.gt.9.d7.and.tk.lt.TburnC.and.yhe4.gt.2.5d-4.and.
     &        yp.lt.1.d-10) then
***   3 alpha ( 0 gamma, 0 gamma)  1 C  12
            v3a = v(i3a)
***   1 C  12 ( 1 alpha, 0 gamma)  1 O  16
            vc12 = v(icag)
***   1 O  16 ( 1 alpha, 0 gamma)  1 Ne 20
            vo16 = v(ioag)
            vHe4 = (v3a*yhe4**2/6.d0+vc12*yc12+vo16*yo16)
            dtHe4 = ftshe0/vHe4
            if (dtHe4.le.dtnx) then
               dtnx = dtHe4
               idtnx = k
               iburning = 5
            endif
            goto 20
         endif

*_______________________________________
***   time-step constrained by C-burning
*---------------------------------------

         if (tk.ge.TburnC.and.yp.lt.1.d-10.and.yhe4.lt.2.5d-4) then
***   2 C  12 ( 0 gamma, 1 alpha)  1 Ne 20
            vcca = v(iccga)
***   2 C  12 ( 0 gamma, 1 proto)  1 Na 23
            vccp = v(iccgp)
***   1 O  16 ( 1 C  12, 1 NEUT )  1 AL 27
            vocn = v(iocn)
***   1 O  16 ( 1 C  12, 1 PROT )  1 AL 27
            vocp = v(iocp)
***   1 O  16 ( 1 C  12, 1 HE  4)  1 MG 24
            voca = v(ioca)
            vcc = abs((vcca+vccp)*yc12*0.5d0+(voca+vocp+vocn)*yo16)
            dtcc = ftsc/vcc
            if (dtcc.le.dtnx.and.yc12.gt.1.d-4) then
               dtnx = dtcc
               idtnx = k
               iburning = 6
            endif
         endif

*________________________________________
***   time-step constrained by Ne-burning
*----------------------------------------

         if (tk.gt.TburnNe) then
***   1 Ne 20 ( 0 gamma, 1 alpha)  1 O  16
            vne = v(inega)
***   1 O  16 ( 1 alpha, 0 gamma)  1 Ne 20
***   inverse only accounted for if T > 1.15 TburnNe     
            if (tk.gt.1.15d0*TburnNe) then
               vneinv = v(ioag)*yhe4*yo16/yne20
            else
               vneinv = 1.d-99
            endif
            dtne = ftsh0/abs(vne-vneinv)
            if (dtne.le.dtnx.and.yne20.gt.1.d-4) then
               dtnx = dtne
               idtnx = k
               iburning = 7
            endif
         endif

*__________________________________________
***   determine smallest burning time scale
*------------------------------------------

         if (k.eq.itmax) then
c           if (dtsi.lt.dtoo.and.dtsi.lt.dtne.and.dtsi.lt.dtcc)
c    &           nphase = 9
c           if (dtoo.lt.dtsi.and.dtoo.lt.dtne.and.dtoo.lt.dtcc)
c    &           nphase = 8
            if (dtne.lt.dtsi.and.dtne.lt.dtoo.and.dtne.lt.dtcc.and.
     &           tmax.gt.TburnNe) nphase = 7
            if (dtcc.lt.dtsi.and.dtcc.lt.dtoo.and.dtcc.lt.dtne.and.
     &           tmax.gt.TburnC) nphase = 6
         endif

 20   continue

*_________________________________________________________
***   determination of the current stellar evolution phase
*---------------------------------------------------------

      if (tmax.lt.TburnC) then
         if (xsp(1,ihe4).lt.1.d-4) then
            nphase = 5
         else
            rcno = xsp(1,ic12)+xsp(1,in14)+xsp(1,io16)
            rcno = rcno/(xsp(nmod1,ic12)+xsp(nmod1,in14)+
     &           xsp(nmod1,io16))
            if (xsp(1,ih1).lt.1.d-4.and.tmax.gt.8.d7.and.(rcno.gt.
     &           1.05d0.or.totm.gt.8.d0)) then
               nphase = 4
            else
               if (xsp(1,ih1).lt.1.d-7) then
                  nphase = 3
               else
                 if (xsp(1,ih1)/xsp(nmod1,ih1).lt.0.998d0.or.
     &                 t(1).gt.3.d7) then
                     nphase = 2
                  else
                     nphase = 1
                  endif
               endif
            endif
         endif
      endif
      if (nphase.lt.nphase0) then
         if (nphase.eq.5) then
c            nphase = nphase0
            write (nout,*) 'star entering the SUPER-AGB phase'
         else
            write (nout,'(" WARNING : PHASE decreasing !! : nphase =",
     &           i2," -->",i2)') nphase0,nphase
         endif
      endif
      if (mtini.gt.12.d0) nphase = max(nphase0,nphase,1)

      if (nphase.le.2.and.lmicro.ge.2.and.lmicro.le.6.and.idiffcc)         ! Modif AP TD 06/2018 puis TD 11/18 puis TD 07/21
     &     microdiffus = .true.
      
      return
      end

************************************************************************

      SUBROUTINE vit (tk,rok,mueinvk,ksh,v,iflag)

************************************************************************
* Calculate the nuclear reaction rates                                 *
* (rho*Navogadro*(sigma*vit)Maxwell, in 1/sec)                         *
*  iflag = 0 : defaut, compute all rates                               *
*  iflag = 1 : compute only rates associated to LiBeB (for HBB)        *
*  iflag = 2 : compute only electron capture rate (Ye dependence)      *
*                                                                      *
* $LastChangedDate:: 2014-02-04 15:45:03 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 11                                                          $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.mod'
      include 'evolcom.nuc'


      integer ivit,ite,iroe,iflag,idv0,iiv0
      integer i,j,jj,ksh

      double precision tburn,tnucmin,fnetdt,fnetdt0,tolnuc,tolnuc0
      double precision tk,rok,mueinvk,v
      double precision t9
      double precision tnuc9,vrate
      double precision fscr
      double precision dt9l,weight,t912,t913,t923,vsec,vmax,dh
      double precision a1,a2,b1,b2,vcapte,enuelec,rhoye,tcapte,romue

      common /funcscr/ fscr(nsh,nscr)
      common /nucvit/ tnuc9(nvit),vrate(nvit,nre)
      common /electcap/ rhoye(nroe),tcapte(nte),vcapte(nte,nroe,nrce),
     &     enuelec(nte,nroe,nrce)
      common /electcap_par/ ite,iroe,a1,a2,b1,b2
      common /nuc_param/ tburn,tnucmin,fnetdt,fnetdt0,tolnuc,tolnuc0

      dimension v(nreac),idv0(nreac)

      t9 = tk*1.d-9

*____________________________
***   Treatment of LiBeB only
*----------------------------

      if (iflag.eq.1) then
         do i = 1,nvit
            if (t9.lt.tnuc9(i)) goto 5
         enddo
 5       ivit = i-1
         if (ivit.eq.0) then
            v(ili7pa) = 1.d-99
            v(ib11pa) = 1.d-99
            v(ili6pa) = 1.d-99
            return
         else
            if (ivit.le.(nvit-1)) then
               dt9l = log10(t9/tnuc9(ivit))
               weight = dt9l/(log10(tnuc9(ivit+1)/tnuc9(ivit)))
               do i = 1,nlibeb
                  j = jlibeb(i)
                  v(j) = vrate(ivit,j)+weight*(vrate(ivit+1,j)-
     &                 vrate(ivit,j))
                  v(j) = 10.d0**v(j)
               enddo
            else
               do i = 1,nlibeb
                  j = jlibeb(i)
                  v(j) = 10.d0**vrate(nvit,j)
               enddo
            endif
         endif
c..   1 LI  6 ( 1 PROT , 1 HE  3)  1 HE  4
         v(ili6pa) = v(ili6pa)*rok*fscr(ksh,3)
c..   1 LI  7 ( 1 PROT , 0 OOOOO)  2 HE  4
         v(ili7pa) = v(ili7pa)*rok*fscr(ksh,3)
c..   1 B  11 ( 1 PROT , 0 OOOOO)  3 HE  4
         v(ib11pa) = v(ib11pa)*rok*fscr(ksh,5)
         return
      endif


*________________________________________
***   Treatment of electron captures only
*----------------------------------------
      if (iflag.eq.2) goto 25


*________________________________________________
***   interpolation of the nuclear reaction rates
*------------------------------------------------

      do i = 1,nvit
         if (t9.lt.tnuc9(i)) goto 10
      enddo
 10   ivit = i-1
      if (ivit.eq.0) then
         do j = 1,nre
            v(j) = 0.d0
         enddo
         do j = 1,nbeta-nbetaec
            jj = jbeta(j)
            v(jj) = vrate(1,jj)
            v(jj) = 10.d0**v(jj)
         enddo
***   1 BE  7 ( 0 betap, 0 nutri)  1 LI  7
         t912 = dsqrt(t9)
         t913 = t9**pw13
         t923 = t913*t913
         vsec = (1.34d-10/t912*(1.d0-0.537d0*t913+3.86d0*t923+2.7d-3/t9*
     &        exp(min(220.d0,2.515d-3/t9))))*rok*mueinvk
         if (t9.le.1.d-3) then
            vmax = 1.5032d-7
            vsec = min(vsec,vmax)
         endif
         v(ibe7beta) = vsec
         if (nphase.gt.2.and.tk.lt.tnucmin) goto 25
         do j = 1,idneut
            jj = jneut(j)
            v(jj) = vrate(1,jj)
            v(jj) = 10.d0**v(jj)
         enddo
         goto 20
      endif
      if (ivit.gt.(nvit-1)) then
         do j = 1,nre
            v(j) = vrate(nvit,j)
            v(j) = 10.d0**v(j)
         enddo
         goto 20
      endif
      dt9l = log10(t9/tnuc9(ivit))
      weight = dt9l/(log10(tnuc9(ivit+1)/tnuc9(ivit)))
      do j = 1,nre
         v(j) = vrate(ivit,j)+weight*(vrate(ivit+1,j)-vrate(ivit,j))
         v(j) = 10.d0**v(j)
      enddo

 20   continue

      iiv0 = 0
      do j = 1,nre
         if (v(j).lt.1.d-97) then
            iiv0 = iiv0+1
            idv0(iiv0) = j
         endif
      enddo

***   1 PROT  ( 1 NEUT , 0 OOOOO)  1 DEUT 
      v(1) = v(1)*rok
***   2 PROT  ( 0 OOOOO, 0 OOOOO)  1 DEUT 
      v(2) = v(2)*rok*fscr(ksh,1)
***   1 DEUT  ( 1 PROT , 0 OOOOO)  1 HE  3
      v(3) = v(3)*rok*fscr(ksh,1)
***   2 DEUT  ( 0 OOOOO, 0 OOOOO)  1 HE  4
      v(4) = v(4)*rok*fscr(ksh,1)
***   2 DEUT  ( 0 OOOOO, 1 NEUT )  1 HE  3
      v(5) = v(5)*rok*fscr(ksh,1)
***   1 HE  3 ( 1 NEUT , 0 OOOOO)  1 HE  4
      v(6) = v(6)*rok
***   1 HE  3 ( 1 DEUT , 1 PROT )  1 HE  4
      v(7) = v(7)*rok*fscr(ksh,2)
***   2 HE  3 ( 0 OOOOO, 2 PROT )  1 HE  4
      v(8) = v(8)*rok*fscr(ksh,19)
***   1 HE  4 ( 1 DEUT , 0 OOOOO)  1 LI  6
      v(9) = v(9)*rok*fscr(ksh,2)
***   1 HE  4 ( 1 HE  3, 0 OOOOO)  1 BE  7
      v(10) = v(10)*rok*fscr(ksh,19)
***   3 HE  4 ( 0 OOOOO, 0 OOOOO)  1 C  12
      v(11) = v(11)*rok*rok*fscr(ksh,19)*fscr(ksh,21)
***   1 LI  6 ( 1 PROT , 0 OOOOO)  1 BE  7
      v(12) = v(12)*rok*fscr(ksh,3)
***   1 LI  6 ( 1 PROT , 1 HE  3)  1 HE  4
      v(13) = v(13)*rok*fscr(ksh,3)
***   1 LI  7 ( 1 PROT , 0 OOOOO)  2 HE  4
      v(14) = v(14)*rok*fscr(ksh,3)
***   1 LI  7 ( 1 DEUT , 1 NEUT )  2 HE  4
      v(15) = v(15)*rok*fscr(ksh,3)
***   1 BE  7 ( 0 betap, 0 nutri)  1 LI  7
      if (ivit.ne.0) then
         vsec = v(ibe7beta)*rok*mueinvk
         if (t9.le.1.d-3) then
            vmax = 1.5032d-7
            vsec = min(vsec,vmax)
         endif
         v(ibe7beta) = vsec
      endif
***   1 BE  7 ( 1 PROT , 0 OOOOO)  1 B   8
      v(17) = v(17)*rok*fscr(ksh,4)
***   1 BE  7 ( 1 DEUT , 1 PROT )  2 HE  4
      v(18) = v(18)*rok*fscr(ksh,4)
***   1 B   8 ( 0 OOOOO, 0 OOOOO)  2 HE  4

***   1 B   8 ( 0 OOOOO, 1 PROT )  1 BE  7

***   1 BE  9 ( 1 PROT , 0 OOOOO)  1 B  10
      v(21) = v(21)*rok*fscr(ksh,4)
***   1 BE  9 ( 1 PROT , 1 DEUT )  2 HE  4
      v(22) = v(22)*rok*fscr(ksh,4)
***   1 BE  9 ( 1 PROT , 1 HE  4)  1 LI  6
      v(23) = v(23)*rok*fscr(ksh,4)
***   1 B  10 ( 1 PROT , 1 HE  4)  1 BE  7
      v(24) = v(24)*rok*fscr(ksh,5)
***   1 B  10 ( 1 PROT , 0 OOOOO)  1 B  11
      v(25) = v(25)*rok*fscr(ksh,5)
***   1 B  11 ( 1 PROT , 0 OOOOO)  1 C  12
      v(26) = v(26)*rok*fscr(ksh,5)
***   1 B  11 ( 1 PROT , 0 OOOOO)  3 HE  4
      v(27) = v(27)*rok*fscr(ksh,5)
***   1 B  11 ( 1 HE  4, 1 PROT )  1 C  14
      v(28) = v(28)*rok*fscr(ksh,22)
***   1 C  12 ( 1 NEUT , 0 OOOOO)  1 C  13
      v(29) = v(29)*rok
***   1 C  12 ( 1 PROT , 0 OOOOO)  1 N  13
      v(30) = v(30)*rok*fscr(ksh,6)
***   1 C  12 ( 1 HE  4, 0 OOOOO)  1 O  16
      v(31) = v(31)*rok*fscr(ksh,23)
***   2 C  12 ( 0 OOOOO, 1 HE  4)  1 NE 20
      v(32) = v(32)*rok*fscr(ksh,35)
***   2 C  12 ( 0 OOOOO, 1 PROT )  1 NA 23
      v(33) = v(33)*rok*fscr(ksh,35)
***   1 C  13 ( 1 NEUT , 0 OOOOO)  1 C  14
      v(34) = v(34)*rok
***   1 C  13 ( 1 PROT , 0 OOOOO)  1 N  14
      v(35) = v(35)*rok*fscr(ksh,6)
***   1 C  13 ( 1 HE  4, 1 NEUT )  1 O  16
      v(36) = v(36)*rok*fscr(ksh,23)
***   1 N  13 ( 0 OOOOO, 0 OOOOO)  1 C  13

***   1 N  13 ( 1 NEUT , 1 PROT )  1 C  13
      v(38) = v(38)*rok
***   1 N  13 ( 1 PROT , 0 OOOOO)  1 N  14
      v(39) = v(39)*rok*fscr(ksh,7)
***   1 C  14 ( 1 NEUT , 0 OOOOO)  1 N  15
      v(40) = v(40)*rok
***   1 C  14 ( 1 PROT , 0 OOOOO)  1 N  15
      v(41) = v(41)*rok*fscr(ksh,6)
***   1 C  14 ( 1 PROT , 1 NEUT )  1 N  14
      v(42) = v(42)*rok*fscr(ksh,6)
***   1 C  14 ( 1 PROT , 1 HE  4)  1 B  11
      v(43) = v(43)*rok*fscr(ksh,6)
***   1 C  14 ( 1 HE  4, 1 NEUT )  1 O  17
      v(44) = v(44)*rok*fscr(ksh,23)
***   1 C  14 ( 1 HE  4, 0 OOOOO)  1 O  18
      v(45) = v(45)*rok*fscr(ksh,23)
***   1 N  14 ( 1 NEUT , 0 OOOOO)  1 N  15
      v(46) = v(46)*rok
***   1 N  14 ( 1 NEUT , 1 PROT )  1 C  14
      v(47) = v(47)*rok
***   1 N  14 ( 1 PROT , 0 OOOOO)  1 O  15
      v(48) = v(48)*rok*fscr(ksh,7)
***   1 N  14 ( 1 HE  4, 0 OOOOO)  1 F  18
      v(49) = v(49)*rok*fscr(ksh,24)
***   1 N  15 ( 1 NEUT , 0 OOOOO)  1 O  16
      v(50) = v(50)*rok
***   1 N  15 ( 1 PROT , 1 HE  4)  1 C  12
      v(51) = v(51)*rok*fscr(ksh,7)
***   1 N  15 ( 1 PROT , 0 OOOOO)  1 O  16
      v(52) = v(52)*rok*fscr(ksh,7)
***   1 N  15 ( 1 HE  4, 0 OOOOO)  1 F  19
      v(53) = v(53)*rok*fscr(ksh,24)
***   1 O  15 ( 0 OOOOO, 0 OOOOO)  1 N  15

***   1 O  15 ( 1 NEUT , 1 HE  4)  1 C  12
      v(55) = v(55)*rok
***   1 O  15 ( 1 NEUT , 1 PROT )  1 N  15
      v(56) = v(56)*rok
***   1 O  16 ( 1 NEUT , 0 OOOOO)  1 O  17
      v(57) = v(57)*rok
***   1 O  16 ( 1 PROT , 0 OOOOO)  1 O  17
      v(58) = v(58)*rok*fscr(ksh,8)
***   1 O  16 ( 1 HE  4, 0 OOOOO)  1 NE 20
      v(59) = v(59)*rok*fscr(ksh,25)
***   1 O  16 ( 1 C  12, 1 NEUT )  1 AL 27
      v(60) = v(60)*rok*fscr(ksh,37)
***   1 O  16 ( 1 C  12, 1 PROT )  1 AL 27
      v(61) = v(61)*rok*fscr(ksh,37)
***   1 O  16 ( 1 C  12, 1 HE  4)  1 MG 24
      v(62) = v(62)*rok*fscr(ksh,37)
***   2 O  16 ( 0 OOOOO, 1 PROT )  1 P  31
      v(63) = v(63)*rok*fscr(ksh,36)
***   2 O  16 ( 0 OOOOO, 1 HE  4)  1 SI 28
      v(64) = v(64)*rok*fscr(ksh,36)
***   1 O  17 ( 1 NEUT , 1 HE  4)  1 C  14
      v(65) = v(65)*rok
***   1 O  17 ( 1 NEUT , 0 OOOOO)  1 O  18
      v(66) = v(66)*rok
***   1 O  17 ( 1 PROT , 0 OOOOO)  1 F  18
      v(67) = v(67)*rok*fscr(ksh,8)
***   1 O  17 ( 1 PROT , 1 HE  4)  1 N  14
      v(68) = v(68)*rok*fscr(ksh,8)
***   1 O  17 ( 1 HE  4, 1 NEUT )  1 NE 20
      v(69) = v(69)*rok*fscr(ksh,25)
***   1 O  17 ( 1 HE  4, 0 OOOOO)  1 NE 21
      v(70) = v(70)*rok*fscr(ksh,25)
***   1 O  18 ( 1 NEUT , 0 OOOOO)  1 F  19
      v(71) = v(71)*rok
***   1 O  18 ( 1 PROT , 0 OOOOO)  1 F  19
      v(72) = v(72)*rok*fscr(ksh,8)
***   1 O  18 ( 1 PROT , 1 HE  4)  1 N  15
      v(73) = v(73)*rok*fscr(ksh,8)
***   1 O  18 ( 1 HE  4, 0 OOOOO)  1 NE 22
      v(74) = v(74)*rok*fscr(ksh,25)
***   1 O  18 ( 1 HE  4, 1 NEUT )  1 NE 21
      v(75) = v(75)*rok*fscr(ksh,25)
***   1 F  18 ( 0 OOOOO, 0 OOOOO)  1 O  18

***   1 F  18 ( 1 NEUT , 1 PROT )  1 O  18
      v(77) = v(77)*rok
***   1 F  18 ( 1 NEUT , 1 HE  4)  1 N  15
      v(78) = v(78)*rok
***   1 F  18 ( 1 HE  4, 1 PROT )  1 NE 21
      v(79) = v(79)*rok*fscr(ksh,26)
***   1 F  19 ( 1 NEUT , 0 OOOOO)  1 NE 20
      v(80) = v(80)*rok
***   1 F  19 ( 1 PROT , 0 OOOOO)  1 NE 20
      v(81) = v(81)*rok*fscr(ksh,9)
***   1 F  19 ( 1 PROT , 1 HE  4)  1 O  16
      v(82) = v(82)*rok*fscr(ksh,9)
***   1 F  19 ( 1 HE  4, 1 PROT )  1 NE 22
      v(83) = v(83)*rok*fscr(ksh,26)
***   1 NE 20 ( 1 NEUT , 0 OOOOO)  1 NE 21
      v(84) = v(84)*rok
***   1 NE 20 ( 1 PROT , 0 OOOOO)  1 NE 21
      v(85) = v(85)*rok*fscr(ksh,10)
***   1 NE 20 ( 0 OOOOO, 1 HE  4)  1 O  16

***   1 NE 20 ( 1 HE  4, 0 OOOOO)  1 MG 24
      v(87) = v(87)*rok*fscr(ksh,27)
***   1 NE 20 ( 1 HE  4, 1 PROT )  1 NA 23
      v(88) = v(88)*rok*fscr(ksh,27)
***   1 NE 21 ( 1 NEUT , 0 OOOOO)  1 NE 22
      v(89) = v(89)*rok
***   1 NE 21 ( 1 NEUT , 1 HE  4)  1 O  18
      v(90) = v(90)*rok
***   1 NE 21 ( 1 PROT , 0 OOOOO)  1 NA 22
      v(91) = v(91)*rok*fscr(ksh,10)
***   1 NE 21 ( 1 HE  4, 0 OOOOO)  1 MG 25
      v(92) = v(92)*rok*fscr(ksh,27)
***   1 NE 21 ( 1 HE  4, 1 NEUT )  1 MG 24
      v(93) = v(93)*rok*fscr(ksh,27)
***   1 NE 22 ( 1 NEUT , 0 OOOOO)  1 NA 23
      v(94) = v(94)*rok
***   1 NE 22 ( 1 PROT , 0 OOOOO)  1 NA 23
      v(95) = v(95)*rok*fscr(ksh,10)
***   1 NE 22 ( 1 HE  4, 1 NEUT )  1 MG 25
      v(96) = v(96)*rok*fscr(ksh,27)
***   1 NE 22 ( 1 HE  4, 0 OOOOO)  1 MG 26
      v(97) = v(97)*rok*fscr(ksh,27)
***   1 NA 22 ( 0 OOOOO, 0 OOOOO)  1 NE 22

***   1 NA 22 ( 1 NEUT , 0 OOOOO)  1 NA 23
      v(99) = v(99)*rok
***   1 NA 22 ( 1 NEUT , 1 PROT )  1 NE 22
      v(100) = v(100)*rok
***   1 NA 22 ( 1 PROT , 0 OOOOO)  1 NA 23
      v(101) = v(101)*rok*fscr(ksh,11)
***   1 NA 23 ( 1 NEUT , 0 OOOOO)  1 MG 24
      v(102) = v(102)*rok
***   1 NA 23 ( 1 PROT , 1 HE  4)  1 NE 20
      v(103) = v(103)*rok*fscr(ksh,11)
***   1 NA 23 ( 1 PROT , 0 OOOOO)  1 MG 24
      v(104) = v(104)*rok*fscr(ksh,11)
***   1 NA 23 ( 1 HE  4, 1 PROT )  1 MG 26
      v(105) = v(105)*rok*fscr(ksh,28)
***   1 MG 24 ( 1 NEUT , 0 OOOOO)  1 MG 25
      v(106) = v(106)*rok
***   1 MG 24 ( 1 PROT , 0 OOOOO)  1 MG 25
      v(107) = v(107)*rok*fscr(ksh,12)
***   1 MG 24 ( 1 HE  4, 1 PROT )  1 AL 27
      v(108) = v(108)*rok*fscr(ksh,29)
***   1 MG 24 ( 1 HE  4, 0 OOOOO)  1 SI 28
      v(109) = v(109)*rok*fscr(ksh,29)
***   1 MG 25 ( 1 NEUT , 0 OOOOO)  1 MG 26
      v(110) = v(110)*rok
***   1 MG 25 ( 1 PROT , 0 OOOOO)  1 AL26m
      v(111) = v(111)*rok*fscr(ksh,12)
***   1 MG 25 ( 1 PROT , 0 OOOOO)  1 AL26g
      v(112) = v(112)*rok*fscr(ksh,12)
***   1 MG 25 ( 1 HE  4, 1 PROT )  1 SI 28
      v(113) = v(113)*rok*fscr(ksh,29)
***   1 MG 25 ( 1 HE  4, 1 NEUT )  1 SI 28
      v(114) = v(114)*rok*fscr(ksh,29)
***   1 MG 25 ( 1 HE  4, 0 OOOOO)  1 SI 29
      v(115) = v(115)*rok*fscr(ksh,29)
***   1 MG 26 ( 1 NEUT , 0 OOOOO)  1 AL 27
      v(116) = v(116)*rok
***   1 MG 26 ( 1 PROT , 0 OOOOO)  1 AL 27
      v(117) = v(117)*rok*fscr(ksh,12)
***   1 MG 26 ( 1 HE  4, 0 OOOOO)  1 SI 30
      v(118) = v(118)*rok*fscr(ksh,29)
***   1 MG 26 ( 1 HE  4, 1 PROT )  1 SI 29
      v(119) = v(119)*rok*fscr(ksh,29)
***   1 MG 26 ( 1 HE  4, 1 NEUT )  1 SI 29
      v(120) = v(120)*rok*fscr(ksh,29)
***   1 AL26m ( 0 OOOOO, 0 OOOOO)  1 MG 26

***   1 AL26m ( 1 NEUT , 1 PROT )  1 MG 26
      v(122) = v(122)*rok
***   1 AL26m ( 1 NEUT , 1 HE  4)  1 NA 23
      v(123) = v(123)*rok
***   1 AL26m ( 1 NEUT , 0 OOOOO)  1 AL 27
      v(124) = v(124)*rok
***   1 AL26m ( 1 PROT , 0 OOOOO)  1 AL 27
      v(125) = v(125)*rok*fscr(ksh,13)
***   1 AL26g ( 0 OOOOO, 0 OOOOO)  1 MG 26

***   1 AL26g ( 1 NEUT , 0 OOOOO)  1 AL 27
      v(127) = v(127)*rok
***   1 AL26g ( 1 NEUT , 1 PROT )  1 MG 26
      v(128) = v(128)*rok
***   1 AL26g ( 1 NEUT , 1 HE  4)  1 NA 23
      v(129) = v(129)*rok
***   1 AL26g ( 1 PROT , 0 OOOOO)  1 AL 27
      v(130) = v(130)*rok*fscr(ksh,13)
***   1 AL26g ( 0 OOOOO, 0 OOOOO)  1 AL26m

***   1 AL26m ( 0 OOOOO, 0 OOOOO)  1 AL26g

***   1 AL 27 ( 1 PROT , 0 OOOOO)  1 SI 28
      v(133) = v(133)*rok*fscr(ksh,13)
***   1 AL 27 ( 1 PROT , 1 HE  4)  1 MG 24
      v(134) = v(134)*rok*fscr(ksh,13)
***   1 AL 27 ( 1 HE  4, 1 NEUT )  1 SI 30
      v(135) = v(135)*rok*fscr(ksh,30)
***   1 AL 27 ( 1 HE  4, 1 PROT )  1 SI 30
      v(136) = v(136)*rok*fscr(ksh,30)
***   1 SI 28 ( 1 NEUT , 0 OOOOO)  1 SI 29
      v(137) = v(137)*rok
***   1 SI 28 ( 1 PROT , 0 OOOOO)  1 SI 29
      v(138) = v(138)*rok*fscr(ksh,14)
***   1 SI 28 ( 1 HE  4, 0 OOOOO)  1 S  32
      v(139) = v(139)*rok*fscr(ksh,31)
***   1 SI 29 ( 1 NEUT , 0 OOOOO)  1 SI 30
      v(140) = v(140)*rok
***   1 SI 29 ( 1 PROT , 0 OOOOO)  1 SI 30
      v(141) = v(141)*rok*fscr(ksh,14)
***   1 SI 30 ( 1 NEUT , 0 OOOOO)  1 P  31
      v(142) = v(142)*rok
***   1 SI 30 ( 1 PROT , 0 OOOOO)  1 P  31
      v(143) = v(143)*rok*fscr(ksh,14)
***   1 P  31 ( 1 NEUT , 0 OOOOO)  1 S  32
      v(144) = v(144)*rok
***   1 P  31 ( 1 PROT , 1 HE  4)  1 SI 28
      v(145) = v(145)*rok*fscr(ksh,15)
***   1 P  31 ( 1 PROT , 0 OOOOO)  1 S  32
      v(146) = v(146)*rok*fscr(ksh,15)
***   1 S  32 ( 1 NEUT , 0 OOOOO)  1 S  33
      v(147) = v(147)*rok
***   1 S  32 ( 1 NEUT , 1 HE  4)  1 SI 29
      v(148) = v(148)*rok
***   1 S  32 ( 1 PROT , 0 OOOOO)  1 S  33
      v(149) = v(149)*rok*fscr(ksh,16)
***   1 S  33 ( 1 NEUT , 0 OOOOO)  1 S  34
      v(150) = v(150)*rok
***   1 S  33 ( 1 NEUT , 1 HE  4)  1 SI 30
      v(151) = v(151)*rok
***   1 S  33 ( 1 PROT , 0 OOOOO)  1 S  34
      v(152) = v(152)*rok*fscr(ksh,16)
***   1 S  34 ( 1 NEUT , 0 OOOOO)  1 S  35
      v(153) = v(153)*rok
***   1 S  34 ( 1 PROT , 0 OOOOO)  1 CL 35
      v(154) = v(154)*rok*fscr(ksh,16)
***   1 S  35 ( 0 OOOOO, 0 OOOOO)  1 CL 35

***   1 S  35 ( 1 NEUT , 0 OOOOO)  1 S  36
      v(156) = v(156)*rok
***   1 S  35 ( 1 PROT , 0 OOOOO)  1 CL 36
      v(157) = v(157)*rok*fscr(ksh,16)
***   1 CL 35 ( 1 NEUT , 0 OOOOO)  1 CL 36
      v(158) = v(158)*rok
***   1 CL 35 ( 1 PROT , 0 OOOOO)  1 HEAVY
      v(159) = v(159)*rok*fscr(ksh,17)
***   1 S  36 ( 1 NEUT , 0 OOOOO)  1 CL 37
      v(160) = v(160)*rok
***   1 S  36 ( 1 PROT , 0 OOOOO)  1 CL 37
      v(161) = v(161)*rok*fscr(ksh,16)
***   1 CL 36 ( 0 OOOOO, 0 OOOOO)  1 HEAVY

***   1 CL 36 ( 1 NEUT , 0 OOOOO)  1 CL 37
      v(163) = v(163)*rok
***   1 CL 36 ( 1 PROT , 0 OOOOO)  1 CL 37
      v(164) = v(164)*rok*fscr(ksh,17)
***   1 CL 37 ( 1 NEUT , 0 OOOOO)  1 HEAVY
      v(165) = v(165)*rok
***   1 CL 37 ( 1 PROT , 0 OOOOO)  1 HEAVY
      v(166) = v(166)*rok*fscr(ksh,17)
***   1 HEAVY ( 1 NEUT , 0 OOOOO)  1 captn
      v(167) = (5.11404d0+signt(ksh))*1.44d5*rok

      do j = 1,iiv0
         v(idv0(j)) = 1.d-99
      enddo


*______________________________________________
***   electron capture rates (rho dependence)
* 166    1 N  14 ( 0 capte, 0 nutri)  1 C  14
* 167    1 C  14 ( 0 betam, 0 nutri)  1 N  14
* 168    1 NE 20 ( 0 capte, 0 nutri)  1 F  20
* 169    1 F  20 ( 0 betam, 0 nutri)  1 NE 20
* 170    1 NA 23 ( 0 capte, 0 nutri)  1 NE 23
* 171    1 NE 23 ( 0 betam, 0 nutri)  1 NA 23
* 172    1 MG 24 ( 0 capte, 0 nutri)  1 NA 24
* 173    1 NA 24 ( 0 betam, 0 nutri)  1 MG 24
* 174    1 MG 25 ( 0 capte, 0 nutri)  1 NA 25
* 175    1 NA 25 ( 0 betam, 0 nutri)  1 MG 25
* 176    1 AL 27 ( 0 capte, 0 nutri)  1 MG 27
* 177    1 MG 27 ( 0 betam, 0 nutri)  1 AL 27
*----------------------------------------------

 25   romue = log10(rok*mueinvk)
      jj = 0

***   extended table (high T and rho)
***   determine index "ite" for temperature
      if (t9.le.tcapte(1)) then
         ite = 1
         jj = 1
         a1 = 0.d0
         b1 = 1.d0
      endif
      if (t9.ge.tcapte(nte)) then
         ite = nte-1
         jj = 2
         a1 = 1.d0
         b1 = 0.d0
      endif
      if (jj.eq.0) then
         do i = 2,nte
            if (t9.lt.tcapte(i)) goto 50
         enddo
 50      ite = i-1
         dh = tcapte(ite+1)-tcapte(ite)
         a1 = (t9-tcapte(ite))/dh
         b1 = 1.d0-a1
      endif
***   determine index "iroe" for rho*Ye
      if (romue.le.rhoye(1)) then
         iroe = 1
         jj = 3
         a2 = 0.d0
         b2 = 1.d0
      endif
      if (romue.ge.rhoye(nroe)) then
         iroe = nroe-1
         jj = 4
         a2 = 1.d0
         b2 = 0.d0
      endif
      if (jj.lt.3) then
         do i = 2,nroe
            if (romue.lt.rhoye(i)) exit
         enddo
         iroe = i-1
         dh = rhoye(iroe+1)-rhoye(iroe)
         a2 = (romue-rhoye(iroe))/dh
         b2 = 1.d0-a2
      endif

***   determine reaction rates (bilinear interpolation)
      do j = nre+1,nreac
         i = j-nre
         v(j) = a1*a2*vcapte(ite+1,iroe+1,i)+b1*b2*
     &        vcapte(ite,iroe,i)+a1*b2*vcapte(ite+1,iroe,i)+
     &        b1*a2*vcapte(ite,iroe+1,i)
         v(j) = 10.d0**v(j)
      enddo

      return
      end

************************************************************************

      SUBROUTINE ydequil

************************************************************************
* Calculate the accretion equilibrium abundance of H1, H2 and He3      *
* * 1  :    2 proto ( 0 gamma, 0 nutri)  1 deutr   -->      used   +   *
* * 2  :    1 deutr ( 1 proto, 0 gamma)  1 he  3   -->      used   -   *
* * 3  :    2 deutr ( 0 gamma, 0 gamma)  1 alpha   --> not  used   -   *
* * 4  :    2 deutr ( 0 gamma, 1 neutr)  1 he  3   --> not  used   -   *
* * 5  :    1 he  3 ( 1 deutr, 1 proto)  1 alpha   --> not  used   -   *
* * 7  :    1 alpha ( 1 deutr, 0 gamma)  1 li  6   --> not  used   -   *
* * 13 :    1 li  7 ( 1 deutr, 1 neutr)  2 alpha   --> not  used   -   *
* * 16 :    1 be  7 ( 1 deutr, 1 proto)  2 alpha   --> not  used   -   *
* * 19 :    1 be  9 ( 1 proto, 1 deutr)  2 alpha   -->      used   +   *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.nuc'
      include 'evolcom.conv'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer mmean,ishell,iend,istart,nconv1,kl0
      integer i,j,k,l,kl,ideut,ivit,ndeut
      integer normx,imin,imax

      double precision sum1
      double precision exv
      double precision tdest,vdest,vdmean,prod,sprod
c     double precision fprod
      double precision xprot,t9,fscrd,fscrbe,xbe,taumean,dtdtaud,xdd,
     &     xdeq,xheq,xpeq,sdest,tauac,taud,xdconv,xdeut,xdburn,dxm,v
c    &     aa0,bb0,edx,edf,slope,v
      double precision dt9l,weight,tnuc9,vrate

      logical partialmix

      parameter (ndeut = 3)

      common /nucvit/ tnuc9(nvit),vrate(nvit,nre)

      dimension tdest(nsh),prod(nsh),taumean(nmaxconv),xdeq(nmaxconv),
     &     xheq(nmaxconv),xpeq(nmaxconv),xdd(nmaxconv),v(ndeut),
     &     ideut(ndeut)

      external exv

      data (ideut(k),k = 1,ndeut)/1,2,19/

      if (nsconv.eq.0) return

*______________________
***   convective region
*----------------------

      iprint(1) = 0
      tauac = 1.d99
      if (nsconv.ge.1) then
c..   restore initial abundances
         do l = 1,nsp
            do k = 1,nmod
               xsp(k,l) = vxsp(k,l)
               ysp(k,l) = xsp(k,l)/anuc(l)
            enddo
         enddo
         do kl = 1,nsconv
            mmin(kl) = novlim(kl,3)
            mmax(kl) = novlim(kl,4)
            tauaccr = 1.d99
            taud = 1.d99
            taumean(kl) = 1.d99
            itauac = -1
            mmean = int((mmin(kl)+mmax(kl))*0.5)
            vdmean = 0.d0
c           fprod = 0.d0
            sprod = 0.d0
            if (itacc.eq.2) then
               idcur = nnucacc(kl)
            else
               idcur = mmin(kl)
            endif
            imin = mmin(kl)
            imax = mmax(kl)
            call mix (dm,imin,imax,kl,1,partialmix)
            nnucacc(kl) = mmin(kl)-1
            if (t(imin).ge.tnucmin) then
               xprot = ysp(mmean,ih1)
               xbe = ysp(mmean,10)
               do 30 k = mmin(kl),mmax(kl)
                  if (t(k).gt.tnucmin) then

*________________________________________________
***   interpolation of the nuclear reaction rates
*------------------------------------------------

                     t9 = t(k)*1.d-9
                     do i = 1,nvit
                        if (t9.lt.tnuc9(i)) goto 10
                     enddo
 10                  ivit = i-1
                     if (ivit.eq.0) then
                        do j = 1,ndeut
                           v(j) = vrate(1,ideut(j))
                           v(j) = 10.d0**v(j)
                        enddo
                        goto 20
                     endif
                     if (ivit.gt.(nvit-1)) then
                        do j = 1,ndeut
                           v(j) = vrate(nvit,ideut(j))
                           v(j) = 10.d0**v(j)
                        enddo
                        goto 20
                     endif
                     dt9l = log10(t9/tnuc9(ivit))
                     weight = dt9l/(log10(tnuc9(ivit+1)/tnuc9(ivit)))
                     do j = 1,ndeut
                        v(j) = vrate(ivit,ideut(j))+(vrate(ivit+1,
     &                       ideut(j))-vrate(ivit,ideut(j)))*weight
                        v(j) = 10.d0**v(j)
                     enddo

 20                  fscrd = 1.d0
                     vdest = v(2)*ro(k)*fscrd
                     tdest(k) = 1.d0/(vdest*xprot)
                     fscrbe = 1.d0
                     prod(k) = (v(3)*xbe*fscrbe+v(1)*xprot*fscrd)*ro(k)
                     if (tdest(k).le.5.d0*tconv(k).and.itacc.lt.3)
     &                    then
                        nnucacc(kl) = k
                     else
                        sprod = sprod+prod(k)*dm(k)
                        vdmean = vdmean+vdest*dm(k)
                     endif
                  endif
 30            continue

               if (itacc.ne.1.or.itacc.ne.2.or.nnucacc(kl).lt.mmin(kl))
     &              then
                  am(kl) = 1.d0/(m(mmax(kl)+1)-m(mmin(kl)))
                  amacc(kl) = macc(mmax(kl)+1)-macc(mmin(kl))
                  if (amacc(kl).gt.0.d0) tauaccr = dtn/(am(kl)*
     &                 amacc(kl))
                  taumean(kl) = tauaccr
               else
                  am(kl) = 1.d0/(m(mmax(kl)+1)-m(nnucacc(kl)+1))
                  amacc(kl) = macc(mmax(kl)+1)-macc(nnucacc(kl)+1)
                  if (amacc(kl).gt.0.d0) tauaccr = dtn/(am(kl)*
     &                 amacc(kl))
                  if (macc(nnucacc(kl)+1).gt.0.d0) taumean(kl) = dtn*
     &                 (m(nnucacc(kl)+1)-m(mmin(kl)))/(macc(nnucacc(kl)+
     &                 1)-macc(mmin(kl)))
               endif
               if (nnucacc(kl).lt.mmax(kl)) then
                  vdmean = vdmean/(m(mmax(kl)+1)-m(nnucacc(kl)+1))
                  if (vdmean.gt.0.d0) taud = 1.d0/(vdmean*xprot)
                  taudnuc = min(taud,taudnuc)
                  dtdtaud = dtn/taud
c                 fprod = sprod/(m(mmax(kl)+1)-m(mmin(kl)))
                  xdd(kl) = vyd1(mmean)*(1.d0-am(kl)*amacc(kl))
                  if (dtdtaud.gt.1.d-9) then
                     xdeq(kl) = xdd(kl)*exv(-dtdtaud)+xspacc(ih2)*
     &                    (1.d0-exv(-dtdtaud))*taud/tauaccr
                  else
                     xdeq(kl) = xdd(kl)+xspacc(ih2)*dtn/tauaccr
                  endif
c                 xdeq(kl) = xprot*taud*fprod*anuc(ih2)+xdd(kl)*
c    &                 exv(-dtdtaud)+xspacc(ih2)*(1.d0-exv(-dtdtaud))*
c    &                 taud/tauaccr
                  xdburn = (xdd(kl)+xspacc(ih2)*taud/tauaccr)*(1.d0-
     &                 exv(-dtdtaud))
                  xpeq(kl) = vxsp(mmean,ih1)-xdburn*anuc(ih1)/
     &                 anuc(ih2)
                  xheq(kl) = vxsp(mmean,ihe3)+xdburn*anuc(ihe3)/
     &                 anuc(ih2)
                  if (itacc.eq.1) xdd(kl) = vyd1(mmean)*(1.d0-dtn/
     &                 taumean(kl))
                  if (nnucacc(kl).le.mmin(kl).or.itacc.ge.3) then
                     do k = mmin(kl),mmax(kl)
                        xsp(k,ih1) = xpeq(kl)
                        xsp(k,ih2) = xdeq(kl)
                        xsp(k,ihe3) = xheq(kl)
                     enddo
                  endif
               endif
            endif
         enddo
      endif

      ldacc = massrate*msun*qdeut*xspacc(ih2)/(anuc(ih2)*sec)

*_____________________
*     radiative region
*---------------------

      ishell = nsconv
      if (crz(1).gt.0.or.nnucacc(1).gt.0) ishell = nsconv+1
      do kl = 1,ishell
         iend = nmod
         if (crz(1).gt.0.or.nnucacc(1).gt.0) then
            if (kl.eq.1) then
               istart = 1
               kl0 = 1
               iend = nnucacc(kl0)
            else
               istart = mmax(kl-1)+1
               kl0 = kl
               if (kl.lt.ishell) iend = nnucacc(kl0)
            endif
         else
            istart = mmax(kl)+1
            kl0 = kl+1
            if (kl.lt.nsconv) iend = nnucacc(kl0)
         endif
         do k = istart,iend
            sdest = 0.d0
            sprod = 0.d0
            taud = 1.d99
            tauac = 1.d99
            if (t(k).gt.tnucmin) then
               if (crz(k).ge.-1.or.itacc.eq.2) then
                  if (facc(k).gt.0.d0) then
                     if (iaccr.eq.2.or.iaccr.eq.4) tauac = dtn*(facc(k)*
     &                    sinthac+1.d0)/(facc(k)*sinthac)
                     if (iaccr.eq.1.or.iaccr.eq.3) tauac = dtn/(facc(k)*
     &                    sinthac)
                  endif
                  if (tauac.lt.tauaccr) then
                     itauac = k
                     tauaccr = tauac
                  endif
               endif
               if (crz(k).ge.-1) then
                  xdeut = yd1(k)

*________________________________________________
***   interpolation of the nuclear reaction rates
*------------------------------------------------

                  t9 = t(k)*1.d-9
                  do i = 1,nvit
                     if (t9.lt.tnuc9(i)) goto 40
                  enddo
 40               ivit = i-1
                  if (ivit.eq.0) then
                     do j = 1,ndeut
                        v(j) = vrate(1,ideut(j))
                        v(j) = 10.d0**v(j)
                     enddo
                     goto 50
                  endif
                  if (ivit.gt.(nvit-1)) then
                     do j = 1,ndeut
                        v(j) = vrate(nvit,ideut(j))
                        v(j) = 10.d0**v(j)
                     enddo
                     goto 50
                  endif
                  dt9l = log10(t9/tnuc9(ivit))
                  weight = dt9l/(log10(tnuc9(ivit+1)/tnuc9(ivit)))
                  do j = 1,ndeut
                     v(j) = vrate(ivit,ideut(j))+(vrate(ivit+1,
     &                    ideut(j))-vrate(ivit,ideut(j)))*weight
                     v(j) = 10.d0**v(j)
                  enddo

 50               fscrd = 1.d0
                  sdest = v(2)*ro(k)*ysp(k,ih1)*fscrd
                  tdest(k) = 1.d0/sdest
                  taudnuc = min(taud,taudnuc)
                  fscrbe = 1.d0
                  xbe = ysp(k,10)
                  sprod = (v(3)*xbe*fscrbe+v(1)*ysp(k,ih1)*fscrd)*ro(k)
               else
                  if (itacc.le.1) then
                     tauac = taumean(kl0)
                     xdeut = xdd(kl0)
                  endif
                  if (itacc.eq.2) xdeut = yd1(k)
                  sprod = prod(k)
               endif
               taud = tdest(k)
               dtdtaud = dtn/taud
               xsp(k,ih2) = xspacc(ih2)*(1.d0-exv(-dtdtaud))*taud/tauac+
     &              xdeut*exv(-dtdtaud)
c    &              xdeut*exv(-dtdtaud)+xsp(k,ih1)*taud*sprod*anuc(ih2)/
c    &              anuc(ih1)
               xdburn = (xdeut+xspacc(ih2)*taud/tauac)*(1.d0-
     &              exv(-dtdtaud))
               xsp(k,ihe3) = vxsp(mmean,ihe3)+xdburn*anuc(ihe3)/
     &              anuc(ih2)
               xsp(k,ih1) = vxsp(mmean,ih1)-xdburn*anuc(ih1)/
     &              anuc(ih2)
            else
               tdest(k) = 1.d50
            endif
         enddo
      enddo

*__________________________
***   linear interpolation
*__________________________

      do 70 kl = 1,nsconv
         if (nnucacc(kl).gt.novlim(kl,3)) rnucacc = r(nnucacc(kl))/
     &        r(nmod)
         if (nnucacc(kl).eq.0.or.itacc.eq.3.or.(nnucacc(kl).lt.
     &        mmin(kl).and.itacc.lt.3).or.t(nnucacc(kl)).lt.
     &        tnucmin) goto 70
         if (facc(mmax(kl)).gt.0.d0.and.tdest(nnucacc(kl)+1).lt.
     &        1.d2*dtn) then
            nconv1 = min(nnucacc(kl)+1,mmax(kl))
            if ((nconv1-nnucacc(kl)).lt.2) goto 70
            dxm = m(nconv1)-m(nnucacc(kl))
            xdconv = (xdeq(kl)*(m(mmax(kl))-m(nnucacc(kl)))-0.5d0*
     &           dxm*xsp(nnucacc(kl),ih2))/(m(mmax(kl))-m(nconv1)+
     &           0.5d0*dxm)
            do k = nnucacc(kl),nconv1
               xsp(k,ih2) = (xdconv-xsp(nnucacc(kl),ih2))*(m(k)-
     &              m(nnucacc(kl)))/dxm+xsp(nnucacc(kl),ih2)
               xsp(k,ihe3) = vxsp(k,ihe3)+anuc(ihe3)/anuc(ih2)*
     &              (vxsp(k,ih2)-xsp(k,ih2))
               xsp(k,ih1) = vxsp(k,ih1)-anuc(ih1)/anuc(ih2)*
     &              (vxsp(k,ih2)-xsp(k,ih2))
            enddo
            if (mmax(kl).gt.nconv1) then
               do k = nconv1+1,mmax(kl)
                  xsp(k,ih2) = xsp(nconv1,ih2)
                  xsp(k,ihe3) = xsp(nconv1,ihe3)
                  xsp(k,ih1) = xsp(nconv1,ih1)
               enddo
            endif
         endif
 70   continue

      if (iter.eq.1) taudnuc0 = taudnuc

*______________________________
***   abundance renormalization
*------------------------------

      do k = 1,nmod
         normx = 1
         sum1 = 0.d0
         do j = 1,nsp
            if (xsp(k,j).gt.xsp(k,normx)) normx = j
            sum1 = sum1+xsp(k,j)
         enddo
         xsp(k,normx) = xsp(k,normx)+(1.d0-sum1)
      enddo

      return
      end


************************************************************************

      SUBROUTINE yequil (v,y,ytol,dtn,taumix)

************************************************************************
*   Set to low abundant elements (y<ytol) to their equilibrium abundance
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.nuc'

      integer i,j
      integer ii,ij,ik,il,nn,mm

      double precision rdest,sprod,am,dtdest,yeqmin,dtn
      double precision v,y,yeq,ytol,taumix

      dimension v(nreac),y(nsp)
      dimension nn(nsp)


      yeqmin = min(1.d-14,ytol*1.d-2)
      do i = 1,nis-1
         if (y(i).lt.yeqmin) then
            sprod = 0.d0
            rdest = 0.d0
            do mm = 1,nreac
               ii = k1(mm)
               nn(ii) = k2(mm)
               ij = k3(mm)
               nn(ij) = k4(mm)
               ik = k5(mm)
               il = k7(mm)
               if (i.eq.ii.or.i.eq.ij) then
                  if (i.eq.ij) then
                     j = ii
                  else
                     j = ij
                  endif
                  if (nn(i).eq.1.and.nn(j).eq.1) then
                     rdest = rdest+v(mm)*y(j)
                  else
                     am = v(mm)*nn(i)/(fact(nn(i))*fact(nn(j)))
                     rdest = rdest+am*y(j)**nn(j)
                  endif
               elseif (i.eq.ik.or.i.eq.il) then
                  if (nn(ii).eq.1.and.nn(ij).eq.1) then
                     sprod = sprod+v(mm)*y(ii)*y(ij)
                  else
                     am = v(mm)*nn(ii)/(fact(nn(ii))*fact(nn(ij)))
                     sprod = sprod+am*y(ii)**nn(ii)*y(ij)**nn(ij)
                  endif
               endif
            enddo
            if (rdest.gt.0.d0) then
               dtdest = 1.d0/rdest
               yeq = sprod/rdest
c..   set y to yeq if dtn > 10 dtdest & tau_mix > dtdest & y < yeqmin
               if (dtn.ge.10.d0*dtdest.and.yeq.lt.yeqmin.and.
     &              taumix.gt.dtdest) then
                  y(i) = yeq
               endif
            endif
         endif
      enddo


      return
      end


************************************************************************

      SUBROUTINE yneutron (y,v,yeq)

************************************************************************
*     Calculate the neutron equilibrium abundance yeq                  *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.nuc'

      integer mm
      integer ii,ij,j

      double precision y,v,yeq,zlocal
      double precision sprod,sdest

      dimension y(*),v(nreac)

***   no neutron captures for extremely metal poor stars (Z = 0)
      if (zkint.lt.5.d-6) then
         zlocal = 0.d0
         do j = 2,ib11
            zlocal = zlocal+y(j)*anuc(j)
         enddo
         zlocal = 1.d0-zlocal
         if (zlocal.lt.1.d-10) then
            yeq = 1.d-50
            return
         endif
      endif

      sdest = 0.d0
!$OMP PARALLEL
!$OMP DO PRIVATE(j) REDUCTION(+:sdest)
      do j = 1,idneut-1
         mm = jneut(j)
         ii = k1(mm)
         sdest = sdest+v(mm)*y(ii)
      enddo
!$OMP END DO

      sprod = 0.5d0*v(iddn)*y(ih2)*y(ih2)
!$OMP DO PRIVATE(j) REDUCTION(+:sprod)
      do j = idneut+1,nneut
         mm = jneut(j)
         ii = k1(mm)
         ij = k3(mm)
         sprod = sprod+v(mm)*y(ii)*y(ij)
      enddo
!$OMP END DO
!$OMP END PARALLEL

      if (sdest.gt.0.d0) then
         yeq = sprod/sdest
      else
         yeq = y(1)
      endif
      yeq = max(abs(yeq),1.d-50)
      return
      end
      Subroutine interpZ
*******************************************************************************************

      implicit none

      include 'evolpar.star'
      include 'evolcom.atm'

      integer i,j,k
      integer filesize
      double precision rtau,rT,rhotab,Fconvtab
      double precision T_Zeff(6),rho_Zeff(6),Fconv_Zeff(6)
      double precision Ttabeff,gtabeff,Ztabeff
      double precision FeH,DT,DDT,sTZ,srhotabZ,sFconvtabZ

      double precision tableTZ,tablerhoZ,tableFcZ

      dimension DT(512),DDT(512)

      common /atmospheres/ rtau(127,12,59,6),rT(127,12,59,6),
     &     rhotab(127,12,59,6),Fconvtab(127,12,59,6),Ttabeff(59),
     &     gtabeff(12),Ztabeff(6),filesize
      common /atmospheres2/ tableTZ(127,59,10),tablerhoZ(127,59,10),
     &     tableFcZ(127,59,10)
      common /metal/ FeH

      print *,FeH

      sTZ = 0.d0
      srhotabZ = 0.d0
      sFconvtabZ = 0.d0


      do i = 1,nb_Teff
         do j = 1,nb_geff
            do k = 1,filesize

               T_Zeff(:) = rT(k,j,i,:)
               call splineatm (Ztabeff,T_Zeff,nb_Zeff,1.d50,1.d50,
     &              DDT,DT)
               call splintatm (Ztabeff,T_Zeff,DDT,nb_Zeff,FeH,sTZ)
               tableTZ(k,i,j) = sTZ

               rho_Zeff(:) = rhotab(k,j,i,:)
               call splineatm (Ztabeff,rho_Zeff,nb_Zeff,1.d50,1.d50,
     &              DDT,DT)
               call splintatm (Ztabeff,rho_Zeff,DDT,nb_Zeff,FeH,
     &              srhotabZ)
               tablerhoZ(k,i,j) = srhotabZ

               Fconv_Zeff(:) = Fconvtab(k,j,i,:)
               call splineatm (Ztabeff,Fconv_Zeff,nb_Zeff,1.d50,1.d50,
     &              DDT,DT)
               call splintatm (Ztabeff,Fconv_Zeff,DDT,nb_Zeff,FeH,
     &              sFconvtabZ)
               tableFcZ(k,i,j) = sFconvtabZ

            enddo
         enddo
      enddo

      return
      end
      subroutine select_table(tau,gstruc,taulim,nmod,tsurf,Tatm
     &     ,rhoecatm,Fcatm)
***********************************************************************************************

      implicit none

      include 'evolpar.star'
      include 'evolcom.atm'
      include 'evolcom.conv'
      include 'evolcom.mass'
      include 'evolcom.var'

      integer i,j,k,l,ij
      integer nmod,nb_elmt
      integer klenv,klpulse,klcore
      integer filesize

      double precision tau
      double precision taulim,taumax,taumin
      double precision sum
      double precision DDT(nmod),DT(nmod)
      double precision DDT1(nmod),DDT2(nmod)
      double precision tsurf,Tatm(nmod),Fcatm(nmod),rhoecatm(nmod)
      double precision gstruc,lgstruc,FeH,FeHatm
      double precision rtau,rT,rhotab,Fconvtab
      double precision Ttabeff,gtabeff,Ztabeff
      double precision nouveauT(59),nouveauFc(59),nouveaurho(59)
      double precision teffinterp(nmod)
      double precision bip
      double precision tableT(nmod,59),tableFc(nmod,59),
     &     tablerho(nmod,59)
      double precision tableTZ,tableFcZ,tablerhoZ
      double precision Fc_atm,rho_atm,interpT
      double precision stau(127),sT(127),srhotab(127),sFconvtab(127)
      double precision srhotabZ,sFconvtabZ,sTZ
      double precision T_geff(10),rho_geff(10),Fconv_geff(10)
      double precision T_Zeff(6),rho_Zeff(6),Fconv_Zeff(6)

      double precision taucoupleatm,rayenv

      dimension sum(59,10)
      dimension tau(nsh)

      common /atmospheres/ rtau(127,12,59,6),rT(127,12,59,6),
     &     rhotab(127,12,59,6),Fconvtab(127,12,59,6),Ttabeff(59),
     &     gtabeff(12),Ztabeff(6),filesize,taumin,taumax
      common /atmospheres2/ tableTZ(127,59,10),tablerhoZ(127,59,10),
     &     tableFcZ(127,59,10)

      common /metal/ FeH
      common /overshoot/ klcore,klenv,klpulse
  

      lgstruc = LOG10(gstruc)
      stau(:) = rtau(:,1,1,1)

      taumax = 100!/mtini
      taumin = 10!/mtini

      taucoupleatm = taulim

      FeHatm = FeH
      DDT = 0.d0
      DT = 0.d0
      sTZ = 0.d0
      srhotabZ = 0.d0
      sFconvtabZ = 0.d0
      
c      print *,'entree dans tableatm',nb_Teff
!!!  2.  Interpolation in gravity
      do i = 1,nb_Teff
         do k=1,filesize
           
            T_geff(:) = tableTZ(k,i,:)
             
            call splineatm (gtabeff,T_geff,nb_geff,1.d50,1.d50,DDT,DT)
            call splintatm (gtabeff,T_geff,DDT,nb_geff,lgstruc,sT(k))

            rho_geff(:) = tablerhoZ(k,i,:)
            call splineatm (gtabeff,rho_geff,nb_geff,1.d50,1.d50,DDT,DT)
            call splintatm (gtabeff,rho_geff,DDT,nb_geff,lgstruc,
     &           srhotab(k))

            Fconv_geff(:) = tableFcZ(k,i,:)
            call splineatm (gtabeff,Fconv_geff,nb_geff,1.d50,1.d50,
     &           DDT,DT)
            call splintatm (gtabeff,Fconv_geff,DDT,nb_geff,lgstruc,
     &           sFconvtab(k))

         enddo

!!!  2.1  Set the values on the same grid as the rest of the code

         Fc_atm = 0.d0
         rho_atm = 0.d0
         DDT = 0.d0
         DT = 0.d0
         call splineatm (stau,sT,filesize,1.d50,1.d50,DDT,DT) !1.d50 to obtain 'natural' limit conditions : second derivative = 0 at 0 and N.
         call splineatm (stau,sFconvtab,filesize,1.d50,1.d50,DDT1,DT)
         call splineatm (stau,srhotab,filesize,1.d50,1.d50,DDT2,DT)
         k=nmod

         do while (tau(k)<=taumax)
            call splintatm (stau,sT,DDT,filesize,tau(k),interpT)
            call splintatm (stau,sFconvtab,DDT1,filesize,tau(k),
     &           Fc_atm)
            call splintatm (stau,srhotab,DDT2,filesize,tau(k),
     &           rho_atm)
            tableT(k,i) = interpT
            tablerho(k,i) = rho_atm
            tableFc(k,i) = Fc_atm
            k=k-1
         enddo
      enddo
      

!!!  3.  Selection of the closest temperature track for tau=[taumin;taumax] and mean it. => tsurf
      Teffinterp = 0.d0
      tsurf = 0.d0
      nb_elmt = 0
      do k=nmod,1,-1
         if (tau(k)>=taumin) exit
      enddo

      rayenv = r(novlim(klenv,4))-r(novlim(klenv,3))
!!      print *,'rayenv=',rayenv/r(nmod)
      if (rayenv.gt.r(nmod)*1.d-2) then
         do while (tau(k)<=taumax)
            nouveauT(:) = tableT(k,:)
            call splineatm (nouveauT,Ttabeff,nb_Teff,1.d50,1.d50,DDT,DT)
            call splintatm (nouveauT,Ttabeff,DDT,nb_Teff,t(k),
     &           Teffinterp(k))
            tsurf = tsurf + Teffinterp(k)
            k = k-1
            nb_elmt = nb_elmt + 1
         enddo
         tsurf = tsurf/(nb_elmt)
      else
         do k=nmod,1,-1
            if (tau(k)>=taucoupleatm) exit
         enddo
         nouveauT(:) = tableT(k,:)
         call splineatm (nouveauT,Ttabeff,nb_Teff,1.d50,1.d50,DDT,DT)
         call splintatm (nouveauT,Ttabeff,DDT,nb_Teff,t(k),
     &        Teffinterp(k))
         tsurf = Teffinterp(k)
      endif

!!!  4.  Creation of the effective T(tau) (and others) track(s) selected to determine q(tau)
      k = nmod
      do while (tau(k)<=taumax)
         nouveauT(:) = tableT(k,:)
         nouveauFc(:) = tableFc(k,:)
         nouveaurho(:) = tablerho(k,:)
         call splineatm (Ttabeff,nouveauT,nb_Teff,1.d50,1.d50,DDT,DT)
         call splintatm (Ttabeff,nouveauT,DDT,nb_Teff,tsurf,Tatm(k))
         call splineatm (Ttabeff,nouveauFc,nb_Teff,1.d50,1.d50,DDT,DT)
         call splintatm (Ttabeff,nouveauFc,DDT,nb_Teff,tsurf,Fcatm(k))
         call splineatm (Ttabeff,nouveaurho,nb_Teff,1.d50,1.d50,DDT,DT)
         call splintatm (Ttabeff,nouveaurho,DDT,nb_Teff,tsurf,
     &        rhoecatm(k))
         k = k-1
      enddo

      end subroutine select_table
