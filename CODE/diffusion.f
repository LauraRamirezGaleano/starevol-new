      SUBROUTINE diffusion (icall,error)

************************************************************************
*     Compute diffusion coefficients and chemical mixing
*
*       The different types of mixing are :
* lmicro = 2   : microscopic diffusion, Chapman & Cowling
* lmicro = 3   : microscopic diffusion, Montmerle&Michaud + Paquette
* lmicro = 4   : atomic diffusion, Thoul et al. 94 + Paquette 86 (ionisation partielle)
* lmicro = 5   : atomic diffusion, Thoul et al. 94 + Paquette 86 (ionisation totale) ne fonctionne pas en l'état !    
* idiffty = 4   : mixing recipe used in Charbonnel ApJ 1995
* idiffty = 8   : chemical mixing + AM transport Talon 1997
* idiffty = 9   : AM transport only, Talon 1997
* idiffty = 10  : (8)+(2)
* idiffty = 11  : (3) + chemical mixing +  AM transport Maeder,Zahn 1998
* idiffty = 13  : (11) + non stationnary terms in AM transport eqs
*                 more complete treatment, recommended for giants
* idiffty = 14  : Omega evolves as a results of structural changes. Angular
*                 momentum assumed to be preserved in each shell
* idiffty = 17  : Use simplified expression for Ur (case of massive stars models)
* lthal =   01  : Thermohaline mixing
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
*    lover = 32 : overshoot everywhere (=30+31)
*    lover = 33 : overshoot below EC : Baraffe et al 2017 penetrative convection (see diffbaraffe.f)
*    lover = 34 : overshoot below a convective zone : Augustson & Mathis 2019 rotating convection model - Baraffe and Pratt type diffusion (see diffkyle.f)
*    lover = 35 : overshoot above a convection zone : Augustson & Mathis 2019 rotating convection model - Baraffe and Pratt type diffusion (see diffkyle.f)
*     lover = 36 : overshoot everywhere (=34+35) : Augustson & Mathis 2019 rotating convection model (see diffkyle.f)
*      lover = 37 : overshoot below a convective zone : Augustson & Mathis 2019 rotating convection model - Baraffe and Pratt type diffusion (see diffkyle.f)   + Dturbul T6.45
*     lover = 38 : overshoot below a convective zone : Augustson & Mathis 2019 rotating convection model - Baraffe and Pratt type diffusion (see diffkyle.f)  + Dturbul PM12500
*     lover = 39 : overshoot below a convective zone : Augustson & Mathis 2019 rotating convection model - Korre type diffusion (see diffkyle.f) 2
C Modif CC ondes (5/12/07)
* lover = 41-43 : Internal gravity waves
*    lover = 41 : igw below the convective envelope
*    lover = 42 : igw above the convective core
*    lover = 43 : igw in both cases (=41+42)
C     Modif TD turbulence ad hoc (09/2019)
*    lover = 60 : compute and add turbulence coefficient (temperature fixation - see diffturbul.f) - only if atomic diffusion active
*    lover = 61 : compute and add turbulence coefficient (BCZ fixation - see diffturbul.f) - only if atomic diffusion active
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
c      double precision xeq(nsp),tnuc(nsp)
      double precision abmuliss,abmax,abmin,gradledoux
C.. test BS99
      double precision THBS,TEMb,F_BS99
      double precision xnurad,xnumol,ddmax
      double precision brunt(nsh)
      double precision dtacho
      double precision Dhold
      double precision rhosh,tsh,muesh ! MODIF THOUL + modif Fev.2019 muesh
c     double precision vdiffthoul
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
c      dimension Dturbul(nsh) ! dturbul TD 09/19
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
c$$$      endif
      if (.not.microdiffus.and..not.rotation.and..not.diffusconv
     $     .and..not.difftacho.and..not.thermohaline.and..not.igw.and.
     $     .not.diffover) then
         print *,'no diffusion process active'
         return
      endif
c      print *,'Entree dans diffusion'

      ndt = 0
      ndb = 1
      times = time*seci
c      call cpu_time(chronos)
c      print *, 'Chrono start diffusion', chronos, 'seconds'
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
C           Dondes(i) = 0.d0
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

      if (.not.diffzc.and.microdiffus.and.klenv.eq.0) then
c         print *,'out of diffusion.f'
         return
      else
c         print *, 'diffusion.f active'
      endif
      
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

c$$$C Modif CC ondes (19/10/07)
c$$$         if (icall.ne.1.and.igw) then
c$$$            print *,"diffusion icall ne 1 initialisation ndb,ndt"
c$$$            if (klenv.gt.0) then
c$$$               ndtenv = novlim(klenv,3)-1
c$$$            else
c$$$               ndtenv = nmod
c$$$            endif
c$$$c.. case envelope boundary condition
c$$$c.. define upper integration boundaries for AM transport
c$$$            if (klcore.eq.0) then
c$$$               ndb = 1
c$$$            else
c$$$               ndb = novlim(klcore,4)+1
c$$$            endif
c$$$            ndt = novlim(klenv,3)-1
c$$$         endif
c$$$C Fin modif CC

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
c$$$            dri1 = r(i)-r(i+1)
c$$$            dri = r(i-1)-r(i)
c$$$            a = -dri1/dri
c$$$            b = dri1+dri
c$$$            dtdr = ((-1.d0/a+a)*t(i)+t(i+1)/a-a*t(i-1))/b
c$$$  dlntdr(i) = dtdr/t(i)
            dri = dsqrt((r(i+1)**2 + r(i+1)*r(i)+r(i)**2)/3.d0) -
     &           dsqrt((r(i)**2 + r(i)*r(i-1) + r(i-1)**2)/3.d0)
            dtdr = (t(i) - t(i-1))/dri
            dlntdr(i) = 2.d0*dtdr/(t(i)+t(i-1))
            if (lmicro.eq.4) dlntdr(i) = dlntdr(i)*rsun ! modif Thoul

c            dlnpdr(i) = -ro(i)/p(i)*g*m(i)/(r(i)*r(i)) ! modif Thoul (pression)
            dpdr = (p(i) - p(i-1))/dri
            dlnpdr(i) = 2.d0*dpdr/(p(i)+p(i-1))
c           dlnpdr(i) = -(0.5*(ro(i)+ro(i-1)))/(0.5*(p(i)+p(i-1)))*g
c     $           *m(i)/(r(i)*r(i))                                ! remise à la couche i, 01/2019
            if (lmicro.eq.4) dlnpdr(i) = dlnpdr(i)*rsun ! modif Thoul
            
            do j = 1,nis        ! modif Thoul (fraction masse)
               dxdr(i,j) = (xsp(i,j) - xsp(i-1,j))/dri
               dlnxdr(i,j) = 2.d0*dxdr(i,j)/(xsp(i,j)+xsp(i-1,j))
               if (lmicro.eq.4) dlnxdr(i,j) = dlnxdr(i,j)*rsun ! modif Thoul
            enddo
c            print *, 'i',i,'dlntdr',dlntdr(i),'dlnpdr',dlnpdr(i)          
c$$$            write(69,'(1x,i4,1x,9(1x,1pe11.4))') i,r(i),t(i),dri,
c$$$     &           dtdr,dlntdr(i),dlnpdr(i),dlnxdr(i,5),m(i)
c$$$            write(65,'(1x,i4,1x,9(1x,1pe11.4))') i,r(i),t(i),p(i),
c$$$     &           m(i),ro(i)
         enddo
         dlntdr(1) = dlntdr(2)
      endif

*____________________________________________________
***   Compute diffusion coefficients
***   coefficients that are independent of abundances
*----------------------------------------------------

c$$$CC***   internal gravity waves
c$$$CC Modif CC (28/10/09)
c$$$C Calcul de Dondeschim dans le cas sans rotation
c$$$C Le calcul de Dondeschim dans le cas du calcul complet
c$$$C est deplace dans la partie "rotational mixing"
c$$$C de la routine diffusion
c$$$      if (igwsurfchim.and..not.igwrot) then
c$$$          mtu = nmod-ndt+1
c$$$          npasr = nmod-ndb+1
c$$$      print *,"diffusion, appel a calculflux_chimsurf,nmod,
c$$$     &         mtu,npasr",nmod,mtu,npasr
c$$$c          call calculflux_chimsurf(nmod,mtu,npasr)
c$$$c          call coefondes_surf(ndb,ndt)
c$$$      stop 'to be implemented'
c$$$      endif
c$$$C Fin modif CC

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
c               print *,'before baraffe',istart,r(istart),cd0
               call diffbaraffe (istart,iend,cd0,Dbar,fover)
            endif
            novlim(klenv,7) = iend
         endif

***   diffusive overshoot : Augustson & Mathis 2019 formalism (code KA, adapt. TD Nov. 2019)
c..   Diffusion below a convection zone         
         if (lover.eq.34.or.(lover.ge.36.and.lover.le.39))
     $        then
            klenv0 = klenv
            istart = novlim(klenv,3)
            if (klenv0.ne.0.and.istart.gt.1) then
               fover = etaturb
               cd0 = Dconv(istart)
c               print *,'before kyle',istart,r(istart),cd0
               if (lover.eq.34.or.(lover.ge.36.and.lover.le.38)) 
     $             call diffkyle (istart,1,iend,cd0,fover,Dkyle,1)
               if (lover.eq.39.) call diffkyle (istart,1,iend,cd0,fover
     $              ,Dkyle,2)
            endif
            novlim(klenv,7) = iend
c            print *, 'istart over',istart,'iend over',iend,'ndt',ndt
            do i = 1,nmod
               if (Dkyle(i).eq.0.d0) then
                  Dkyle(i) = 1.d0
               endif
            enddo
         endif

c..   Diffusion above a convection zone
         if (lover.eq.35.or.lover.eq.36) then
            klenv0 = klenv
            istart = novlim(klenv,4)
            if (klenv0.ne.0.and.istart.gt.1) then
               fover = etaturb
               cd0 = Dconv(istart)
c               print *,'before kyle',istart,r(istart),cd0
               call diffkyle (istart,-1,iend,cd0,fover,Dkyle)
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
c               write (90,200) istart,iend,t(iend),
c     &              abs(m(istart)-m(iend))/msun,fover
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
c                  write (90,200) istart,iend,t(iend),abs(m(istart)-
c     &                 m(iend))/msun,fover
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
c                  write (90,200) istart,iend,t(iend),abs(m(istart)-
c     &                 m(iend))/msun,fover
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
c               write (90,200) istart,iend,t(iend),abs(m(istart)-m(iend))
c     &              /msun,fover
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

C Modifs CC ondes (2/04/07)
c               if (igw) call diffinitondesexcit
C Fin modifs CC

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
c                  print *,'ndt=',ndt
                  if (mtini.ge.15.d0.and.xsp(1,ih1).le.5.d-6) omegaconv
     $                 = 'm'
                  call dif_omega (ndb,ndt,neqj,ncls,nclc,times,error)
               endif
            endif

c..   Check that surface velocity has not change too much
c..   (in case of AM transport by IGW)
c            print *, 'aie aie aie'
c            print *,omega(nmod),vomega(nmod),
c     &           dabs((omega(nmod)-vomega(nmod))/vomega(nmod))
c            if (igwrot.and.
c     &           dabs((omega(nmod)-vomega(nmod))/vomega(nmod)).gt.1d-2)
c     &           error = 50
            do i = ndb+1,ndt-1
               if (omega(i)*omega(i-1).lt.0.d0.and.omega(i)*omega(i+1)
     &              .lt.0.d0) error = 50
            enddo
            if (error.gt.0) return
         endif
         
c.. Tachocline diffusion -> Modif. TD 30/09/2019 - no restriction on phase - specify the choice of the ZC         
c         if (difftacho.and.nphase.eq.2.and.klenv.gt.1) then
         if (difftacho) then   
            if (nsconv.eq.1.and.novlim(1,3).ne.1) then
               print *, 'One ZC', 'Dtacho active'
c     ndt = novlim(klenv,3)-1
               ndt = novlim(1,3)-1
               call tachocline (ndt,1,coefDtacho)
c               print *,'ndt-1=',ndt-1
            else if (nsconv.gt.1.and.novlim(2,3).ne.1) then
c               print *, 'Nbr ZC=',nsconv, 'Dtacho active'
               ndt = novlim(2,3)-1
               call tachocline (ndt,1,coefDtacho)
c               print *,'ndt-1=',ndt-1
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
c         if (idiffty.eq.31) then
c             vturb = 1.d5
c             xomega = vturb/r(nmod)
c             call dturbucc95 (ndt,ndb,xomega,cd)
c         endif
      endif
c$$$      do jshell=1,nmod
c$$$         print *, mue(jshell),mu(jshell),muizero(i)
c$$$      enddo
c$$$  stop
c$$$      do j=1,nmod
c$$$         write(44,*) mue(j)
c$$$      enddo
c$$$      stop
      if (microdiffus.and.lmicro.eq.4) then
         do jshell = 1,nmod
c            tsh = t(jshell)
c            rhosh = ro(jshell)
            tsh = 0.5*(t(jshell)+t(jshell-1))     ! remise à la couche i, 01/2019
            rhosh = 0.5*(ro(jshell)+ro(jshell-1)) ! remise à la couche i, 01/2019
            muesh = mue(jshell) ! ajout TD AP Fev.2019
c     call diffthoul (jshell,xsp,anuc,znuc,tsh,rhosh)
c            print *, 'jshell', jshell
            call diffthoul (jshell,xsp,anuc,znuc,tsh,rhosh,muesh,abond,
     $           nmod)
         enddo
      endif

C Ajout TD 09/2019 - Compute turbulence coefficient     
***   compute the diffusion coefficient of turbulence (Richer et al. 2000; Richard et al.2005)
      if
     $     ((lover.ge.60.and.lover.le.61.and.microdiffus).or.(lover.eq.
     $     37.and.microdiffus).or.(lover.eq.38.and.microdiffus)) then
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
      
***   compute total diffusion coefficient for chemical species

 20   forall (i = 1:nmod) dturb(i) = cd(i)+Dconv(i)+Dsc(i)+Dherw(i)
     &     +coefDtacho(i)+Dhd(i)+Dthc(i)+Dturbul(i)+Dbar(i)+Dkyle(i)

      if (nmixd.gt.0.and.icall.eq.1) return

      if (icall.eq.1) write (nout,500)
      if (omegaconv.eq.'s'.and.hydrorot) goto 45

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
      
c      stop
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
c         if (ielem.eq.47) print *,'abond(1,47) avt diffsolve'
c     $        ,abond(1,47)
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
c               print *, 'abond helium av', kl,abond(1,5)
               call diffsolve(abond,ratneh,anuc,znuc,dd,f1,f2,f3,dmdr,
     &              kl,nshell,icall,error)
            endif
            if (error.gt.0) return
         enddo
c      print *, 'abond helium', abond(1,5)
c..   normalize and update abundances for nuclear calculation and
c..   for binary output
         summx = 0.d0
         nxmax = 1
c         print *, 'le rassemblement du corbeau !'         
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
                  if (xsp(j,k).gt.xsp(j,inorm)) inorm = k ! changement d'elem ref si H pas le plus abondant
               enddo
c$$$                if (icall.eq.1) write(777,*) j,abs(1.d0-xm_sum)  ! ecriture de summx
               xsp(j,inorm) = xsp(j,inorm)+(1.d0-xm_sum) ! element ref (H en general) absorbe la difference pour retomber à 1
               if (icall.eq.1) vxsp(j,inorm) = xsp(j,inorm)
c$$$               abond(j,inorm) = xsp(j,inorm)
               abond(j,inorm) = xsp(j,inorm)*anucni/anuc(inorm)  !modif +
c               if (inorm.ne.2) then
c                  write (nout,*) 'Atomic diffusion : most abundant ' //
c     &                 'species at shell ',j,'is not H but ',elem(inorm)
c     &                 , ' diffusion stopped'
c                  stop 'diffusion : microscopic'
c     endif
               summx = abs(1.d0-xm_sum)
            else
               do k = 2,nis
c                 if (j.le.10) print *,'xm_sum',xm_sum,
c     $                 xsp(j,k)
                  xsp(j,k) = abond(j,k)
c                  if (j.le.10) print *,'xm_sum',xm_sum,xsp(j,k),j,k
                  xm_sum = xm_sum+xsp(j,k)
                  if (xsp(j,k).gt.xsp(j,inorm)) inorm = k
               enddo
               if (abs(1.d0-xm_sum).gt.summx) then
                  summx = abs(1.d0-xm_sum)
                  nxmax = j
               endif
c               if (icall.eq.1) print *,'inorm,j,xsp(j,inorm)',inorm,j
c     $              ,xsp(j,inorm),abond(j,inorm),abond(j,ihe4),xm_sum
c     $              ,summx,xsp(j,1),xsp(j,2),xsp(j,3),xsp(j,4)
               xsp(j,inorm) = xsp(j,inorm)+(1.d0-xm_sum)
               if (icall.eq.1)  vxsp(j,inorm) = xsp(j,inorm)
               abond(j,inorm) = xsp(j,inorm)
c               print *,'ds diffusion, abond et xsp He et P',xsp(j,inorm)
c     $              ,xsp(j,ihe4),xsp(j,44)
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
c      write(398,*)model
c      do i = 1,nmod
c         write(398,*)i,xsp(i,2),crz(i)
c      enddo
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
c      do i = 1,ndt
c         print *,omega(i),vomega(i)
c      enddo
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
c$$$      if(icall.eq.1) then
c$$$         do i=ndt,0,-1
c$$$            write(554,'(2(1x,i4),4(1x,1pe11.4))'),i
c$$$     $           ,crz(i),r(i),t(i),ro(i)
c$$$         enddo
c$$$      endif
     
c      print *, 'Chrono end diffusion', chronos, 'seconds'

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
