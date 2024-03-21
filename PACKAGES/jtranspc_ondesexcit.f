
************************************************************************
* This part has been developed by:                                     *
*    Suzanne Talon and Ana Palacios                                    *
*    Version 2.92: Dec 2006                                            *
* Modifs ondes : Corinne Charbonnel (22/11/06)                         *
*                                                                      *
* $LastChangedDate:: 2018-01-22 10:28:22 +0000 (Mon, 22 Jan 2018)    $ *
* $Author:: amard                                                    $ *
* $Rev:: 115                                                         $ *
*                                                                      *
************************************************************************

      SUBROUTINE DIF_OMEGA (ndb,ndt,neqj,ncls,nclc,times,error)

*-----------------------------------------------------------------------
* Cette routine permet de faire le calcul de l'evolution du moment
* cinetique.  Le profil de rotation interne va aussi fournir
* le coefficient de transport turbulent.

* DIF_OMEGA est appelee par DIFFUSION,
* et appelle NORM, CONS_L, CI, LAMBDA, HEN_OMEGA, CDIFF, ROT_SOL,
*            ASSYMPT3, CAL_DH
*
*
* Arguments:
*
*     ndt = borne superieure de la zone de melange
*     neqj = nombre d'equations du systeme a resoudre
*     ncls = nombre de CL de surface associees au syst. a resoudre
*     nclc = nombre de CL au centre associees au syst. a resoudre
*     times = age du modele (en annees)
*-----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.data'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.grad'
      include 'evolcom.nuc'
      include 'evolcom.mod'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.transp'
      include 'evolcom.igw'
c      include 'evolcom.transpondesexcit'


      logical implicit_wave

      integer iiter,ivar,ipass,iiterondes
      integer ndt,ndb
      integer i,j,neqj,ncls,nclc
      integer klcore,klenv,klpulse

      integer pcol,error
      integer idiffvr_sav

      double precision omega,vomega,k2conv,k2rad,vsurf,
     &     angradr,angconvr
      double precision om,ur,vom,vrray,times,vxpsi
      double precision aux,xpsi,xlambda,xtheta
      double precision xjmod
      double precision dift,dife,Dhold
      double precision xflux_shear,xflux_cir
      double precision fluxshe,fluxcir
      double precision xmoment_tot
      double precision alph,dtom
      double precision rhmoy,epsmoy
      double precision nturb
      double precision r,vr
      double precision Deff
      double precision abmurj,abmuj
      integer ndtenv
      double precision dshear
c      double precision depotwaves

      logical init,MMM
      logical zgradmu_sav
      character omegaconv_sav*1


      common /moment/ xmoment_tot
      common /rotvar/ omega(nsh),vomega(nsh),k2conv,k2rad,vsurf,
     &     angradr,angconvr
      common /calcul/ iiter
      common /difcirc/ dift(nsh),dife(nsh)
      common /fluxomega/ fluxshe(nsh),fluxcir(nsh)
      common /moy/ rhmoy(nsh),epsmoy(nsh)
      common /dturbulente/ nturb(0:nsh)
      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)
      common /varrvr/ r(nsh),vr(nsh)
      common /calcDh/ Dhold(nsh)
      common /gradmuj/ abmurj(nsh),abmuj(nsh)
      common /init_rot/ init,omegaconv_sav
      common /envel/ ndtenv
      common /overshoot/ klcore,klenv,klpulse
      common /chemdiffshear/ dshear(nsh)

      double precision omdom(nsh)

      dimension om(nsh),ur(nsh),aux(nsh),xpsi(nsh),xlambda(nsh)
      dimension xjmod(nsh,neqj)
      dimension pcol(neqj)
      dimension Deff(nsh)
      double precision omprint(nmod)
c      dimension depotwaves(nsh)

*------------------------
***   Initialisations
*------------------------

      print *,'entree dans jtranspc'

      print *,idiffvr

      ipass = 0
      idiffvr_sav = idiffvr
      init = inits
      zgradmu_sav = zgradmu
      omegaconv_sav = omegaconv

      if (init.and.zgradmu) then
         zgradmu_sav = zgradmu
         zgradmu = .false.
      endif

      nturb(0) = 0.d0
      do i = 1,nmod
         nturb(i) = 0.d0
         Dhd(i) = 0.d0
         dift(i) = 0.d0
         xKt(i) = khi(i)/(cp(i)*rro(i))
c         if (Dconv(i).gt.1.d8.and.i.gt.klenv) 
c     &        Dconv(i) = max(Dconv(i),1.d17)
      enddo

c..   Initialisation du pointeur de colonnes pour ordonnancement
c..   des matrices

      if (zgradmu) then
         pcol(1) = 1
         pcol(2) = 2
         pcol(3) = 4
         pcol(4) = 5
         pcol(5) = 3
      else
         pcol(1) = 1
         pcol(2) = 2
         pcol(3) = 3
         pcol(4) = 4
      endif

      if (idiffvr.ge.6) init = .false.
      if (idiffvr.eq.7.and.nphase.le.2) then
         idiffvr = 0
      endif

*     Added in order to mimic the sequence adopted by MM at the end
*     of central H burning


      if (idiffcc.and.omegaconv.eq.'m') then
         MMM = .true.
          call assympt3 (ndt,om,ur)
          goto 10
      endif

      print *,'omegaconv=',omegaconv
*----------------------------------------------------------------
***   Rotation solide imposee dans la zone radiative de melange
*----------------------------------------------------------------

      if (omegaconv.eq.'s') then
         call rot_sol (ndt,om,aux,ur,xpsi,vom,times)
         goto 10
      endif


*--------------------------------------------------------------------
***   Calcul du tout premier modele avec rotation (si rotation solide
***   des ZC, calcul du profil de omega a l'eq. thermique)
*--------------------------------------------------------------------
      if (init) then
c$$$         if (nphase.eq.2) then
c$$$            vsurf = diffvr
c$$$            if (idiffvr.lt.6.or.idiffvr.eq.8) then
c$$$               omega_S = diffvr*1.d5/rtot
c$$$               if (idiffvr.eq.5.or.idiffvr.eq.8) then
c$$$                  idiffvr = 0
c$$$               endif
c$$$            else
c$$$               omega_S = vom(nmod)
c$$$            endif
c$$$         else
c$$$            if (idiffvr.eq.4) then
c$$$               omega_S = breaktime
c$$$            else
c$$$               omega_S = vsurf*1.d5/rtot
c$$$            endif
c$$$         endif

         if (nphase.eq.1.and.idiffvr.eq.4) then
            omega_S = breaktime
         endif

*     Calcule de la vitesse initiale de circulation et normalisation
*     des variables

         call ci (ndb,ndt,neqj,ncls,nclc,pcol,times,xjmod,error)

c..   Save profiles with distinction for massive hot stars when omegaconv = T

         if (omegaconv.eq.'f') then

            do i = 2,nmod-1
               if (crz(i-1).le.1.and.crz(i).le.1.and.crz(i+1).le.1) then
                  om(i) = (wi(i)*(xjmod(i-1,pcol(1))+xjmod(i,pcol(1)))+
     $                 wj(i)*(xjmod(i+1,pcol(1))+xjmod(i,pcol(1))))
     $                 *0.5d0
                  aux(i) = (wi(i)*(xjmod(i-1,pcol(2))+xjmod(i,pcol(2)))+
     $                 wj(i)*(xjmod(i+1,pcol(2))+xjmod(i,pcol(2))))
     $                 *0.5d0
                  ur(i) = (wi(i)*(xjmod(i-1,pcol(3))+xjmod(i,pcol(3)))+
     $                 wj(i)*(xjmod(i+1,pcol(3))+xjmod(i,pcol(3))))
     $                 *0.5d0
                  xpsi(i) = (wi(i)*(xjmod(i-1,pcol(4))+xjmod(i,pcol(4)))
     &                 +wj(i)*(xjmod(i+1,pcol(4))+xjmod(i,pcol(4))))*
     &                 0.5d0
                  if (idiffvr.eq.6.and.i.gt.ndt) then
                     xpsi(i) = (phiKSm(i)*xlambda(i)+pw43*om(nmod)*
     $                    rray(nmod)**2*om(i)/grav(i)/rray(i))
     $                    /deltaKSm(i)
                  else if (idiffvr.eq.7.and.i.gt.ndt) then
c     call ash_rgb_prof(ndt,i,om,omdom,xpsi)
                     call ash_rgb_prof(ndt,om,omdom,xpsi)
                     xpsi(i) = (phiKSm(i)*xlambda(i)-pw23*omdom(i)
     $                    *rray(i)**2/grav(i))/deltaKSm(i)
                  endif
               else
                  om(i) = xjmod(i,pcol(1))
                  aux(i) = xjmod(i,pcol(2))
                  ur(i) = xjmod(i,pcol(3))
c     if (idiffvr.eq.6) then
c     xpsi(i) = (phiKSm(i)*xlambda(i)+pw43*om(nmod)*
c     &                 rray(nmod)**2*om(i)/(grav(i)*rray(i)))
c     &                 /deltaKSm(i)
c     else if (idiffvr.eq.7) then
c     call ash_rgb_prof(ndt,om,omdom)
c     xpsi(i) = (phiKSm(i)*xlambda(i)-pw23*omdom(i)
c     $                 *rray(i)**2/grav(i))/deltaKSm(i)
c     else
                  xpsi(i) = xjmod(i,pcol(4))
c     endif
               endif
               if (zgradmu_sav) xlambda(i) = 0.d0
            enddo
            om(nmod) = xjmod(nmod,pcol(1))
            aux(nmod) = xjmod(nmod,pcol(2))
            xpsi(nmod) = xjmod(nmod,pcol(4))
            ur(nmod) = ur(nmod-1)

         else if (omegaconv.eq.'t') then
            do i = 2,nmod-1
c...  Treat radiation zone between surface ionisation regions as convective
               if (crz(i-1).le.1.and.crz(i).le.1.and.crz(i+1).le.1.and.
     &              i.lt.novlim(klenv,4)) then
                  om(i) = (wi(i)*(xjmod(i-1,pcol(1))+xjmod(i,pcol(1)))+
     $                 wj(i)*(xjmod(i+1,pcol(1))+xjmod(i,pcol(1))))
     $                 *0.5d0
                  aux(i) = (wi(i)*(xjmod(i-1,pcol(2))+xjmod(i,pcol(2)))+
     $                 wj(i)*(xjmod(i+1,pcol(2))+xjmod(i,pcol(2))))
     $                 *0.5d0
                  ur(i) = (wi(i)*(xjmod(i-1,pcol(3))+xjmod(i,pcol(3)))+
     $                 wj(i)*(xjmod(i+1,pcol(3))+xjmod(i,pcol(3))))
     $                 *0.5d0
                  xpsi(i) = (wi(i)*(xjmod(i-1,pcol(4))+xjmod(i,pcol(4)))
     $                 +wj(i)*(xjmod(i+1,pcol(4))+xjmod(i,pcol(4))))
     $                 *0.5d0
                  if (idiffvr.eq.6.and.i.gt.ndt) then
                     xpsi(i) = (phiKSm(i)*xlambda(i)+pw43*om(nmod)*
     &                    rray(nmod)**2*om(i)/grav(i)/rray(i))/
     &                    deltaKSm(i)
                  else if (idiffvr.eq.7.and.i.gt.ndt) then
c     call ash_rgb_prof(ndt,i,om,omdom,xpsi)
                     call ash_rgb_prof(ndt,om,omdom,xpsi)
                     xpsi(i) = (phiKSm(i)*xlambda(i)-pw23*omdom(i)
     $                    *rray(i)**2/grav(i))/deltaKSm(i)
                  endif
               else if (i.ge.novlim(klenv,4)) then
                  om(i) = xjmod(novlim(klenv,3)+1,pcol(1))
                  aux(i) = xjmod(novlim(klenv,3)+1,pcol(2))
                  ur(i) = xjmod(novlim(klenv,3)+1,pcol(3))
                  xpsi(i) = xjmod(novlim(klenv,3)+1,pcol(4))
               else
                  om(i) = xjmod(i,pcol(1))
                  aux(i) = xjmod(i,pcol(2))
                  ur(i) = xjmod(i,pcol(3))
c     if (idiffvr.eq.6) then
c     xpsi(i) = (phiKSm(i)*xlambda(i)+pw43*om(nmod)*
c     &                 rray(nmod)**2*om(i)/(grav(i)*rray(i)))
c     &                 /deltaKSm(i)
c     else if (idiffvr.eq.7) then
c     call ash_rgb_prof(ndt,om,omdom)
c     xpsi(i) = (phiKSm(i)*xlambda(i)-pw23*omdom(i)
c     $                 *rray(i)**2/grav(i))/deltaKSm(i)
c     else
                  xpsi(i) = xjmod(i,pcol(4))
c     endif
               endif
               if (zgradmu_sav) xlambda(i) = 0.d0
            enddo
            om(nmod) = xjmod(novlim(klenv,4),pcol(1))
            aux(nmod) = xjmod(novlim(klenv,4),pcol(2))
            xpsi(nmod) = xjmod(novlim(klenv,4),pcol(4))
            ur(nmod) = ur(nmod-1)
         endif
         om(1) = om(2)
         aux(1) = aux(2)
         ur(1) = xjmod(1,pcol(3))
         xpsi(1) = xpsi(2)
         if (zgradmu_sav) xlambda(1) = xlambda(2)

         if (error.gt.0) then
            zgradmu = zgradmu_sav
            model_old = model_old-1
            do j = 1,nmod
               aux(j) = 1.d-50
               xpsi(j) = 1.d-50
               xlambda(j) = 0.d0
               ur(j) = 0.d0
            enddo
            return
         endif

         do i = 1,nmod
            Dhd(i) = 0.d0
         enddo

      else

*-----------------------------------------------------------------------
***   Evolution du profil de omega
***   Determination des coefficients de diffusion associes a la rotation
*-----------------------------------------------------------------------

*     Normalisation des variables

c         if (idiffvr.ne.6) then
            ivar = 2
c         else
c            ivar = 3
c         endif
         alph = 0.1d0
         
         call norm (om,ur,aux,xpsi,xlambda,vom,vxpsi,vrray,ivar,ndb,ndt)
         


*     Calcul du coefficient horizontal de diffusion turbulente

 5       call cal_dh (ndb,ndt,ur,om,xpsi,xlambda,1)
         xmoment_tot = xmom_tots
c* Version pour laquelle le pas de temps d'evolution du moment cinetique
c* est identique au pas de temps d'evolution de l'etoile.

         dtom = dtn

*     Initialisation de la matrice a taille variable xjmod pour la 
*     resolution du syst. d'eq. par Newton-Raphson

C -----------------------
c Update tranport by IGW
C -----------------------


         implicit_wave = .false.
c         iiterondes = iterh
         if (igwrot) then
            if (implicit_wave.or.iiterondes.eq.1) 
     &           call transport_ondes (om,ndb,ndtenv,dtom,
     &           iiterondes)
c            implicit_wave = .false.           
         else
            do i = 1,nmod
               depotwaves(i) = 0.d0
            enddo
         endif
     

         do j = 1,nmod
            xjmod(j,pcol(1)) = om(j)
            xjmod(j,pcol(2)) = aux(j)
            xjmod(j,pcol(3)) = ur(j)
            xjmod(j,pcol(4)) = xpsi(j)
            if (zgradmu) xjmod(j,pcol(5)) = xlambda(j)
         enddo
         
         call hen_omega (ndb,ndt,neqj,ncls,nclc,pcol,dtom,alph,times,
c     &        xjmod,error,depotwaves)
     &        xjmod,error)
         

*     On n'a pas converge: on recalcule le modele avec un pas de temps 
*     d'evolution plus petit

         if (error.gt.0) then
            if (ipass.eq.0) then
               ipass = 1
               error = 0
               alph = 0.1d0
               nturb(0) = 0.d0
               do i = 1,nmod
                  nturb(i) = 0.d0
                  Dhd(i) = 0.d0
                  dift(i) = 0.d0
               enddo
               goto 5
            else
               model_old = model_old-1
               omegaconv = omegaconv_sav
            endif
            return
         endif

*     Mise a jour des differents vecteurs
         om(1) = xjmod(1,pcol(1))
         aux(1) = xjmod(1,pcol(2))
         ur(1) = xjmod(1,pcol(3))
         xpsi(1) = xjmod(1,pcol(4))
         if (zgradmu) xlambda(1) = xjmod(1,pcol(5))

         if (omegaconv.eq.'f') then
            do i = 2,ndt-1
c     if (crz(i-1).lt.0.and.crz(i).lt.0.and.crz(i+1).lt.0) then
               if ((crz(i-1).le.1.and.crz(i).le.1.and.crz(i+1).le.1)) 
     &              then
                  om(i) = (wi(i)*(xjmod(i-1,pcol(1))+xjmod(i,pcol(1)))+
     &                wj(i)*(xjmod(i+1,pcol(1))+xjmod(i,pcol(1))))*0.5d0
                  aux(i) = (wi(i)*(xjmod(i-1,pcol(2))+xjmod(i,pcol(2)))+
     &                wj(i)*(xjmod(i+1,pcol(2))+xjmod(i,pcol(2))))*0.5d0
                  ur(i) = (wi(i)*(xjmod(i-1,pcol(3))+xjmod(i,pcol(3)))+
     &                wj(i)*(xjmod(i+1,pcol(3))+xjmod(i,pcol(3))))*0.5d0
                 xpsi(i) = (wi(i)*(xjmod(i-1,pcol(4))+xjmod(i,pcol(4)))+
     &                wj(i)*(xjmod(i+1,pcol(4))+xjmod(i,pcol(4))))*0.5d0
              else
                 om(i) = xjmod(i,pcol(1))
                 aux(i) = xjmod(i,pcol(2))
                 ur(i) = xjmod(i,pcol(3))
                 xpsi(i) = xjmod(i,pcol(4))
              endif
              if (zgradmu) xlambda(i) = xjmod(i,pcol(5))
           enddo
           om(ndt) = xjmod(ndt,pcol(1))
           aux(ndt) = xjmod(ndt,pcol(2))
           xpsi(ndt) = xjmod(ndt,pcol(4))
           ur(ndt) = ur(ndt-1)
c         om(ndt) = om(ndt-1)
        else if (omegaconv.eq.'t') then
            do i = 2,ndt-1
               if ((crz(i-1).le.1.and.crz(i).le.1.and.crz(i+1).le.1)
     &              .and.i.lt.novlim(klenv,4)) then
                  om(i) = (wi(i)*(xjmod(i-1,pcol(1))+xjmod(i,pcol(1)))+
     &                wj(i)*(xjmod(i+1,pcol(1))+xjmod(i,pcol(1))))*0.5d0
                  aux(i) = (wi(i)*(xjmod(i-1,pcol(2))+xjmod(i,pcol(2)))+
     &                wj(i)*(xjmod(i+1,pcol(2))+xjmod(i,pcol(2))))*0.5d0
                  ur(i) = (wi(i)*(xjmod(i-1,pcol(3))+xjmod(i,pcol(3)))+
     &                wj(i)*(xjmod(i+1,pcol(3))+xjmod(i,pcol(3))))*0.5d0
                 xpsi(i) = (wi(i)*(xjmod(i-1,pcol(4))+xjmod(i,pcol(4)))+
     &                wj(i)*(xjmod(i+1,pcol(4))+xjmod(i,pcol(4))))*0.5d0
              else if (i.ge.novlim(klenv,4)) then
                 om(i) = xjmod(novlim(klenv,3)+1,pcol(1))
                 aux(i) = xjmod(novlim(klenv,3)+1,pcol(2))
                 ur(i) = xjmod(novlim(klenv,3)+1,pcol(3))
                 xpsi(i) = xjmod(novlim(klenv,3)+1,pcol(4))                 
              else
                 om(i) = xjmod(i,pcol(1))
                 aux(i) = xjmod(i,pcol(2))
                 ur(i) = xjmod(i,pcol(3))
                 xpsi(i) = xjmod(i,pcol(4))
              endif
              if (zgradmu) xlambda(i) = xjmod(i,pcol(5))
           enddo
           om(ndt) = xjmod(novlim(klenv,3)+1,pcol(1))
           aux(ndt) = xjmod(novlim(klenv,3)+1,pcol(2))
           xpsi(ndt) = xjmod(novlim(klenv,3)+1,pcol(4))
           ur(ndt) = ur(ndt-1)           
        endif
         V_circ(1) = 0.d0

c..  definition of integration variables in conv. zones (when not included)
         if (ndb.gt.1) then
            om(1) = om(2)
            xpsi(1) = xpsi(2)
            do j = 2,ndb-1
               om(j) = om(ndb)
               xpsi(j) = xpsi(ndb)
               ur(j) = 0.d0
               xlambda(j) = 0.d0
               aux(j) = 0.d0
               V_circ(j) = 0.d0
            enddo
         endif
         if (ndt.lt.nmod-1) then
            do j = ndt+1,nmod
               if (idiffvr.lt.6.or.idiffvr.ge.8) then
                  if (idiffvr.eq.4) then
                     om(j) = breaktime
                     xpsi(j) = 0.d0
                  endif
                  om(j) = om(ndt)
                  xpsi(j) = xpsi(ndt)
               else if (idiffvr.eq.6) then
                  xpsi(j) = (phiKSm(j)*xlambda(j)+pw43*om(ndt)*
     &                 rray(ndt)**2*om(j)/(grav(j)*rray(j)))/deltaKSm(j)
                  om(j) = om(ndt)*rray(ndt)**2/rray(j)**2
               else if (idiffvr.eq.7) then
c                  call ash_rgb_prof(ndt,j,om,omdom,xpsi)
                  call ash_rgb_prof(ndt,om,omdom,xpsi)
                  xpsi(j) = (phiKSm(j)*xlambda(j)-pw23*omdom(j)*
     &                 rray(ndt)**2/grav(j))/deltaKSm(j)
               endif
               xlambda(j) = 0.d0
               aux(j) = 0.d0
               ur(j) = 0.d0
               V_circ(j) = 0.d0
            enddo
         endif
      endif

*------------------------------------------------------
***   Computation of Deff for diffusion of chemicals
*------------------------------------------------------

 10   do j = 1,nmod
         Deff(j) = 0.d0
         dife(j) = 0.d0
      enddo
         
      call cal_dh (ndb,ndt,ur,om,xpsi,xlambda,2)
         
      do i = max(ndb,2),ndt
         if (Dhold(i).gt.0.d0) then
            Deff(i) = (rray(i)*ur(i)*ur_S*rtot)**2/(Dhold(i)*30.d0)   ! see eq. 26 Decressin et al. 2009
         else
            Deff(i) = 0.d0
         endif
         Dhd(i) = dshear(i)+Deff(i)
         dife(i) = Deff(i)
      enddo
      Dhd(ndb) = dshear(ndb)+Deff(ndb+1)
      dife(ndb) = Deff(ndb+1)
      Dhd(ndt) = dshear(ndt)+Deff(ndt)
      dife(ndt) = Deff(ndt)

c      call lissage (ur,ndt,1,5)

      inits = init

*     Initialisation de la variable xmom_tots a la fin du calcul du
*     tout premier modele avec rotation

      if (init) then
         ivar = 1
         call cons_l (om,xmoment_tot,ivar,ndb,ndt)
         xmom_tots = xmoment_tot
         zgradmu = zgradmu_sav
         idiffvr = idiffvr_sav
      endif

*     Stockage des profils pour sauvegarde et calcul au pas de temps suivant

      do i = 1,nmod
         omega(i) = om(i)*omega_S
         auxs(i) = aux(i)
         urs(i) = ur(i)
         xpsis(i) = xpsi(i)
         oms(i) = om(i)
         if (zgradmu.and..not.inits) xlambdas(i) = xlambda(i)
      enddo

      ur_Ss = ur_S
      ndtold = ndt
      do j = 2,nmod-1
         if (zgradmu) then
            xtheta = phiKSm(j)*xlambda(j)-deltaKSm(j)*xpsi(j)
         else
            xtheta = -deltaKSm(j)*xpsi(j)
         endif
         xflux_shear = 1.5d-36*grav(j)*rray(j)**2*dift(j)*xtheta/om(j)*
     &        rro(j)*omega_S/grav(nmod)
         xflux_cir = 2.d-37*rro(j)*rray(j)**4*om(j)*ur(j)*ur_S*
     &        omega_S*rtot**4
         fluxshe(j) = xflux_shear
         fluxcir(j) = xflux_cir
      enddo

      return
      end



** ----------------------------------------------------------------------------
*
      SUBROUTINE CI (ndb,ndt,neqj,ncls,nclc,pcol,times,xjmod,error)
*
** ----------------------------------------------------------------------------
* Initialisation des variables de rotation et calcul du premier modele
* avec rotation
*
* Calcul du profil de rotation qui assure l'equilibre thermique si idiffvr < 6
*
* Arguments:
*     ndt = borne superieure de la zone de melange
*     neqj = nbre d'eq du syst. a resoudre pour calculer le transport du moment
*            cinetique
*     ncls = nombre de CL de surface associees
*     nclc = nombre de CL au centre associees
*     times = age du modele (en annees)
*     xjmod(nsh,neqj) = matrice regroupant les vecteurs associes aux variables
*             independantes du syst.
*
*     NOTA BENE: Par defaut, on aura xjmod(i,1) = om(i)
*                                    xjmod(i,2) = aux(i)
*                                    xjmod(i,3) = ur(i)
*                                    xjmod(i,4) = xpsi(i)
*                TOUTE NOUVELLE VARIABLE INDEPENDANTE SERA RAJOUTEE A LA SUITE
*                exemple: xjmod(i,5) = xlambda(i) si zgradmu = T
*
*
* ATTENTION : Ne pas imposer de flux a la surface sinon il est
*             impossible d'obtenir om=vom (on tendra vers la
*             solution pour laquelle le flux est partout egal
*             au flux impose a la surface).
*
* 2.384989d6 = 1/(2*pi*G)
*
* Auteur: S.Talon
*
* Adaptation: A.Palacios (15/03/2000)
*
C Modifs CC ondes (22/11/06)
*
* Derniere version: 03 mai 2000
* -----------------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.grad'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.transp'
      include 'evolcom.igw'
      include 'evolcom.var'


      integer ndt,neqj,ncls,nclc,ndb
      integer i,j,nconv,idmax
      integer pcol,error
      integer klcore,klenv,klpulse

      logical convergence

      double precision xpsi,aux,om,ur,times,xlambda
      double precision vom,vrray,vxpsi
      double precision rhmoy,epsmoy
      double precision divis,dif
      double precision aalpha,bbeta,tamp4,raprho,rapeps
      double precision alph,deltat,ddmax
      double precision xjmod
      double precision ur_smax,ur_smin
      double precision tol
c      double precision depotwaves

      common /moy/ rhmoy(nsh),epsmoy(nsh)
      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)
      common /overshoot/ klcore,klenv,klpulse

      dimension om(nsh),aux(nsh),ur(nsh),xpsi(nsh),xlambda(nsh)
      dimension xjmod(nsh,neqj)
      dimension pcol(neqj)
c      dimension depotwaves(nsh)

*--------------------------------------
***   Normalisation des variables
*--------------------------------------

      tol = 1.d-5
      divis = 1.d0/rtot
      do j = 1,nmod
         aux(j) = 0.d0
         xpsi(j) = 0.d0
         xlambda(j) = 0.d0
         ur(j) = 0.d0
         hht(j) = hht(j)*divis
         hhp(j) = hhp(j)*divis
         if (idiffvr.lt.6) then
            om(j) = 1.d0
*     On initialise un omega deja normalise
            vom(j) = 1.d0
         else 
            vom(j) = vom(j)/vom(nmod)
            om(j) = vom(j)*vrray(j)*vrray(j)/(rray(j)*rray(j))
         endif
      enddo

*----------------------------------------------------
***   Calcul de la vitesse de circulation initiale
*----------------------------------------------------

      coeff = rtot*rtot
      xLtot = coeff*omega_S
      coeff = 4.d0*pi*rtot*coeff
      xLloc = coeff
      xnorm = rtot*omega_S*omega_S/(gs*gs)
      call val_moy

      do j =  max(ndb,2),ndt
c         if (crz(j).lt.0) then
          if (crz(j).le.1.or.j.ge.novlim(klenv,3)) then
            ur(j) = 0.d0
         else
            aalpha = lum(j)*pm(j)
            divis = -grav(j)*m(j)*cpm(j)*tm(j)*rkonv(j)
            aalpha = aalpha/(rdensm(j)*divis)*xnorm
            if (j.ne.2) then
               raprho = 0.5d0*(rro(j)/rhmoy(j)+rro(j-1)/rhmoy(j-1)) 
               rapeps = 0.5d0*(enucl(j)/epsmoy(j)+enucl(j-1)/
     &              epsmoy(j-1))
            else
               raprho = rro(j-1)/rhmoy(j)
               rapeps = enucl(j-1)/epsmoy(j)
            endif
            bbeta = 2.d0*(pw43-raprho) 
            tamp4 = 1.d0-omega_S*omega_S*2.384989d6/rdensm(j)-rapeps
            ur(j) = aalpha*bbeta*rray(j)/grav(j)*tamp4
         endif
      enddo
      if (ndb.eq.1.or.crz(ndb).le.1) ur(ndb) = 0.d0
c      if (ndb.eq.1.or.crz(ndb).lt.0) ur(ndb) = 0.d0
      if (ndt.eq.nmod1) ur(ndt) = 0.d0

*     Normalisation de la vitesse de circulation
      ur_smax = abs(urs(ndt))
      do j = ndt-1,1,-1
         ur_smax = max(ur_smax,abs(urs(j)))
      enddo
      ur_S = ur_smax*ur_Ss/2.d0
      ur_smin = 1.d-8
      ur_S = max(ur_smin,ur_S)
      if (abs(ur_S).gt.1.d-40) ur_Ss = ur_S
      do i = 1,ndt
         ur(i) = ur(i)/ur_S
      enddo

*     Initialisation de la matrice a taille variable xjmod pour la 
*     resolution du syst. d'eq. par Newton-Raphson

      do j = 1,nmod
         xjmod(j,pcol(1)) = om(j)
         xjmod(j,pcol(2)) = aux(j)
         xjmod(j,pcol(3)) = ur(j)
         xjmod(j,pcol(4)) = xpsi(j)
         if (zgradmu) xjmod(j,pcol(5)) = xlambda(j)
      enddo

*     Si on considere un moment specifique constant dans les ZC, on sort
      if (.not.thermal_equilibrium.or.idiffvr.ge.6) return


*------------------------------------------------------------------------
***   Calcul du profil de rotation initial a l'equilibre thermique
*------------------------------------------------------------------------

      do j = 1,nmod
        vrray(j) = rray(j)
        vrraym(j) = rraym(j)
        vrraym2(j) = rraym2(j)
      enddo

      convergence = .false.
      deltat = dtn
      alph = 0.15d0
      nconv = 0

      do i=1,nmod
         depotwaves(i)=0.d0
      enddo

      do while (.not.convergence)
         nconv = nconv+1

         call hen_omega (ndb,ndt,neqj,ncls,nclc,pcol,deltat,alph,times,
c     &        xjmod,error,depotwaves)
     &        xjmod,error)

         if (error.gt.0) return

         alph = 0.2d0
         if (ndb.gt.1) then
            do j = 1,ndb-1
               xjmod(j,pcol(1)) = xjmod(ndb,pcol(1))
            enddo
         endif
         if (ndt.lt.nmod) then
            do j = ndt+1,nmod
               xjmod(j,pcol(1)) = xjmod(ndt,pcol(1))
            enddo
         endif
         ddmax = 0.d0
         do j = 1,nmod
            dif = abs(xjmod(j,pcol(1))-vom(j))
            if (ddmax.lt.dif) then
               idmax = j
               ddmax = dif
            endif
            vom(j) = xjmod(j,pcol(1))
*     On garde ainsi vom(nmod) = 1-Vsurf = diffvr
         enddo
         if (ddmax.lt.tol) then
            convergence = .true.
         else
            deltat = deltat*1.7d0
            write (nout,100) ddmax,tol,idmax
         endif
c         if (ddmax.lt.1.d-1) convergence = .true.
c         if (ddmax.lt.1.d-5) convergence = .true.
      enddo

 100  format (' [CI] relaxation failed : ddmax = ',1pe10.4,' < ',
     &     1pe9.3,' [',i4,']')

      return
      end


** ---------------------------------------------------------------------
*
      SUBROUTINE CAL_DH (ndb,ndt,ur,om,xpsi,xlambda,ipar)
*
** ---------------------------------------------------------------------
* Calculation of Dh
*
* Author: S.Talon (Toulouse-Geneva stellar evolution code)
*
* La variable Dh_prescr est un ensemble de caracteres indiquant la prescription
* utilisee pour evaluer Dh
*
* ipar = 1 calcul Dhold a la premiere iteration
* ipar = 2 calcul Dh apres convergence rotation
*
* Adaptation: A.Palacios (15/03/2000)
*
* Last version: 15/03/2000
**----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.diff'
      include 'evolcom.grad'
      include 'evolcom.mod'
      include 'evolcom.rot'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.transp'

      integer ifin,j,j1,k
      integer ndt,ipar,ndb
      integer iiter

      double precision ur,om,xpsi
      double precision rtoturs
      double precision disctime,Dhold,Ch
      double precision tamp,tamp1,tamp2,tamp3,bbeta,coeffA,tamp4
      double precision rr2j1,rr2j
      double precision vgrav,vgs,V_circm
      double precision vom,vrray,vxpsi,vrdensm
      double precision r,vr
      double precision r_norm,om_norm
      double precision vhnenn
      double precision omm,omm1
!      double precision omega,vomega,k2conv,k2rad,vsurf,angradr,
!     &     angconvr
      double precision xlambda(nsh)
      
      double precision omgs,rays
      double precision Ric,Nt2(nsh),Nu2(nsh),BruntV2(nsh),xNt(nsh)
      double precision geffom(nsh),xNu(nsh),abmurj,abmuj,eps,coefomega
      double precision xnumol,xnurad,xNom,xnuvvm,Rec
      double precision intsin,term,term1,invS,sign,tautransp
      double precision terma,termb,termc,delta,x1,x2
      double precision term2,termden,termnum,termS,intsinsav
      real*4 Re

      common /disclocking/ disctime
      common /calcDh/ Dhold(nsh)
      common /calcul/ iiter
      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)
      common /varrvr/ r(nsh),vr(nsh)
!      common /rotvar/ omega(nsh),vomega(nsh),k2conv,k2rad,vsurf,
!     &     angradr,angconvr
      common /gradmuj/ abmurj(nsh),abmuj(nsh)
      common /diffvisc/ xnumol(nsh),xnurad(nsh)

      dimension ur(nsh),om(nsh),xpsi(nsh),vgrav(nsh),
     &     vhnenn(nsh),vrdensm(nsh)
      dimension omgs(nsh),rays(nsh)


*-------------------------------------------
***   Choix de la prescription pour Dh
*-------------------------------------------

      do j = 1,nmod
         V_circ(j) = 0.d0
         Dhold(j) = 0.d0
         xalpha(j) = 0.d0
         Nt2(j) = 0.d0
      enddo
         
      ifin = min(nmod1,ndt)

      if (ipar.eq.1) then

         rtoturs = vr(nmod)*ur_S
         r_norm = vr(nmod)
         om_norm = vomega(nmod)

         vgs = g*m(nmod)/(vr(nmod)*vr(nmod))
         vrdensm(1) = vrro(1)
         vgrav(1) = 0.d0
         do j = 2,nmod
            j1 = j-1
            vgrav(j) = g*m(j)/(vgs*vr(j)*vr(j))
            vrdensm(j) = wj(j)*vrro(j)+wi(j)*vrro(j1)
            vhnenn(j) = -1.d0/(vrray(j)-vrray(j1))
         enddo
         
         do j = max(2,ndb),ifin
            j1 = j-1
            rr2j1 = vrdensm(j1)*vrray(j1)**2
            rr2j = vrdensm(j)*vrray(j)**2
         
            V_circ(j) = (rr2j*ur(j)-rr2j1*ur(j1))*vhnenn(j)/(6.d0*
     &           vrro(j1)*vrraym(j))

c..  In estimating Dh, variations of mu (xlambda) neglected
            xalpha(j) = 1.d0-0.75d0*deltaKSm(j)*vxpsi(j)*vgrav(j)/
     &           (vom(j)**2*vrray(j))

            rays(j) = vrray(j)
            omgs(j) = vom(j)
         enddo

      else

         rtoturs = r(nmod)*ur_S
         r_norm = r(nmod)
         om_norm = omega_S

         do j = max(2,ndb),ifin
            j1 = j-1
            rr2j1 = rdensm(j1)*rray(j1)**2
            rr2j = rdensm(j)*rray(j)**2
         
            V_circ(j) = (rr2j*ur(j)-rr2j1*ur(j1))*hnenn(j)/(6.d0*
     &           rro(j1)*rraym(j))

c..  In estimating Dh, variations of mu (xlambda) neglected
            xalpha(j) = 1.d0-0.75d0*deltaKSm(j)*xpsi(j)*grav(j)/
     &           (om(j)**2*rray(j))

            rays(j) = rray(j)
            omgs(j) = om(j)
         enddo
      endif
c.. Marques & Goupil 2013 bbeta = 1.5d-6
c      bbeta = 1.5d-6*rtoturs*omega_S*r_norm**2
      bbeta = 2.d-6*rtoturs*omega_S*r_norm**2
      coeffA = (1.5d-3/pi)**pw13*rtot*rtoturs*om_norm*ur_S
c      Ric = 0.25d0
      Ric = pw16
      do j = max(2,ndb),ifin
         j1 = j-1
         V_circm = wj(j)*V_circ(j)+wi(j)*V_circ(j1)
         tamp = 2.d0*V_circm-xalpha(j)*ur(j)
         tamp1 = tamp*tamp+ur(j)*ur(j)
         tamp2 = dsqrt(tamp1)
         
         geffom(j) = gs*grav(j)-pw23*rray(j)*rtot*(om(j)*omega_S)**2
         if (rkonv(j).ge.-del_zc) then
            xNt(j) = geffom(j)/hhp(j)/rtot*del_zc*deltaKSm(j)
         else
            xNt(j) = -geffom(j)/hhp(j)/rtot*rkonv(j)*deltaKSm(j)
         endif
         xNu(j) = gs*grav(j)*abmurj(j)*phiKSm(j)
         
cc         coefomega = xKtm(j)*tamp**2
cc         Re = (coefomega)/(xNt(j)*xnumol(j))
         
         xnuvvm = 0.5d0*(xnuvv(j)+xnuvv(j1))
         Rec = 7.d0*xnum(j)

*     Prescription Zahn (1992)
         if (Dh_prescr.eq.'Zahn1992') then
            Ch = 1.d0
            Dhold(j) = rays(j)*tamp2/Ch*rtoturs

*     Prescription Mathis (2016)
c         else if (Dh_prescr.eq.'Mathis16') then
c            Nt2(j) = gs*grav(j)*deltaKS(j)*(abad(j)-abla(j))/(hp(j))
c            Dhold(j) = 0.75d0*Ric*xKtm(j)*(Nt2(j)/(4.d0*omega(j)**2))

            
**** M17+ - tau = 1/S
         else if (Dh_prescr.eq.'Mathis16') then
            eps = pw23

            tamp3 = omgs(j)*rays(j)*rays(j)
            tamp4 = sqrt(bbeta*tamp3*rays(j)*tamp2)
            if ((xnuvvm.le.Rec).and.(rkonvm(j).lt.0.d0).and.
     &           (xnuvv(j).ne.0.d0)) then
               Dhold(j) = tamp4
            else
               Dhold(j) = tamp4 + eps*Ric*xKtm(j)*(xNt(j)/(4.d0*
     &           (omgs(j)*omega_S)**2))
            endif
            
c$$$ Re = 1/RiPr = (Kt*S^2) / (\nu_m*N2)
c            tamp = (phiKSm(j)*xlambda(j)-deltaKSm(j)*xpsi(j))
c     &           /omgs(j)
c            coefomega = xKtm(j)/xNt(j)*9.d0/4.d0*(grav(j)*omega_S
c     $           /rray(j))**2*tamp**2
c            Re = (coefomega)/(xNt(j)*xnumol(j))
!!!            print *,'Reynold=',j,Re


*     Prescription Mathis (2017) bis
****       M17+ tau = 1/(2\Omega+S)
         else if (Dh_prescr.eq.'Mathis02') then              
            eps = 1.d0
c$$$  Re = 1/RiPr = (Kt*S^2) / (\nu_m*N2)
            omm = 0.5d0*(om(j)+om(j-1)) 
            omm1 = 0.5d0*(om(j)+om(j+1)) 
            termS = rray(j)*(omm1-omm)/(2.d0*omgs(j)*(rraym(j+1)
     &           -rraym(j)))   
            if (dabs(termS).le.1.d0) then
               term1 = 1.d0-termS**2
               term2 = termS/dsqrt(term1)
               termnum = ((termS-pi)*term1+termS)*dsqrt(term1)
               termnum = termnum+
     &              (3.d0*term1-1.d0)*(0.5d0*pi-datan(term2))
            else
               term1 = termS**2-1.d0
               term2 = dabs((termS-dsqrt(term1))/(termS+dsqrt(term1)))
               termnum = ((termS-pi)*term1-termS)*dsqrt(term1)
               termnum = termnum+
     &              0.5d0*(2.d0-3.d0*termS**2)*dlog(term2)
            endif
c$$$***   Ana coding -->
c$$$***   using the relation between omega derivative and phi*lambda - delta*psi
c$$$            
c$$$            termS = pw34*(phiKSm(j)*xlambda(j)-deltaKSm(j)*xpsi(j))
c$$$     &           *geffom(j)/(gs*rray(j)*omgs(j)**2)
c$$$c            print *,termS,sign*tamp,j
c$$$            if (dabs(termS).le.1.d0) then
c$$$               term1 = 1.d0-termS**2
c$$$               term2 = termS/dsqrt(term1)
c$$$               termnum = ((termS-pi)*term1+termS)*dsqrt(term1)
c$$$               termnum = termnum+
c$$$     &              (3.d0*term1-1.d0)*(0.5d0*pi-datan(term2))
c$$$            else
c$$$               term1 = termS**2-1.d0
c$$$               term2 = dabs((termS-dsqrt(term1))/(termS+dsqrt(term1)))
c$$$               termnum = ((termS-pi)*term1-termS)*dsqrt(term1)
c$$$               termnum = termnum+
c$$$     &              0.5d0*(-1.d0-3.d0*term1)*dlog(term2)
c$$$            endif
c$$$***   Ana coding <--
            termden = termS*term1**1.5d0
            if (termS.ne.0.d0.and.(termS.gt.-1.d0)) then
               intsin = termnum/termden
!     !               intsinsav = intsin
               intsin = min(1.d0,intsin)
            else
               print *,'termS < -1',termS
               intsin = 1.d0
c               intsin = intsinsav
            endif
            tamp3 = omgs(j)*rays(j)*rays(j)
            tamp4 = sqrt(bbeta*tamp3*rays(j)*tamp2)
            if ((xnuvvm.le.Rec).and.(rkonvm(j).lt.0.d0).and.
     &           (xnuvv(j).ne.0.d0)) then
               Dhold(j) = tamp4
            else
               Dhold(j) = tamp4 + eps*Ric*pw23*xKtm(j)*xNt(j)*intsin/
     &              ((2.d0*omgs(j)*omega_S)**2.d0)
            endif

!!            Dhold(j) = eps*Ric*pw23*xKtm(j)*xNt(j)*intsin/
!!     &           ((2.d0*omgs(j)*omega_S)**2.d0)
               
****     M17+ - tau = 1/N_Omega
         else if (Dh_prescr.eq.'Mathiepi') then
            eps = 1.d0
            intsin = 1.d0
c$$$ Re = 1/RiPr = (Kt*S^2) / (\nu_m*N2)
***   Ana coding -->
***   using the relation between omega derivative and phi*lambda - delta*psi
c            termS = pw34*(phiKSm(j)*xlambda(j)-deltaKSm(j)*xpsi(j))
c     &           *geffom(j)/(gs*rray(j)*omgs(j)**2)
            omm = 0.5d0*(om(j)+om(j-1)) 
            omm1 = 0.5d0*(om(j)+om(j+1))
            termS = rray(j)*(omm1-omm)/(2.d0*omgs(j)*(rraym(j+1)
     &           -rraym(j)))    ! gradient en j
            if (termS.ne.0.d0) term1 = sqrt(abs(1.d0+1.d0/termS))
            if (termS.lt.0.d0.and.termS.gt.-1.d0) then
               term2 = 1.d0/term1
               termnum = atan(term2)
            else if (termS.gt.0.d0) then
               term2 = 1.d0/(1.d0+term1)
               termnum = 0.5d0*log(abs((1.d0-term1)*term2))
            endif
***   Ana coding <--
               
            termden = termS*term1
            if (termS.ne.0.d0.and.(termS.gt.-1.d0)) then
               intsin = 1.d0+termnum/termden
               intsin = min(1.d0,intsin)
               intsin = max(-1.d0,intsin)
            else
               intsin = 1.d0
            endif
            
c$$$            if (intsin*termS.lt.0.d0) print *,'Dh < 0 ..?',intsin,
c$$$     &           termS*intsin

            tamp3 = omgs(j)*rays(j)*rays(j)
            tamp4 = sqrt(bbeta*tamp3*rays(j)*tamp2)
            if ((xnuvvm.le.Rec).and.(rkonvm(j).lt.0.d0).and.
     &           (xnuvv(j).ne.0.d0)) then
               Dhold(j) = tamp4
            else
               Dhold(j) = tamp4 + eps*Ric*pw23*xKtm(j)*xNt(j)*termS*
     &              intsin/((2.d0*omgs(j)*omega_S)**2.d0)           
            endif

c$$$            Dhold(j) = eps*Ric*pw23*xKtm(j)*xNt(j)*abs(termS*intsin)
c$$$     &           /((2.d0*omgs(j)*omega_S)**2.d0)
            
*     Prescription Mathis, Palacios & Zahn (2004)
         else if (Dh_prescr.eq.'MPZ_2004') then
            Ric=pw16

            tamp3 = omgs(j)*rays(j)*rays(j)
c            Dhold(j) = sqrt(bbeta*tamp3*rays(j)*dabs(tamp))
            Dhold(j) = sqrt(bbeta*tamp3*rays(j)*tamp2)
            print *,'omega*r^2,2V-alphaU',tamp3,tamp2,bbeta

*     Prescription Maeder (2003) ou Maeder (2006)
         else if (Dh_prescr.eq.'Maeder03'.or.Dh_prescr.eq.'Maeder06')
     &           then
            Dhold(j) = coeffA*rays(j)*abs(rays(j)*omgs(j)*
     &           V_circm*tamp2)**pw13
            if (Dh_prescr.eq.'Maeder06')
     &           Dhold(j) = 0.002d0*Dhold(j)
         endif

c         if (nphase.gt.1.and.Dhold(j).gt.Dh_max*xKtm(j))
         Dh_max = 10.d0
c         if (Dhold(j).gt.Dh_max*xKtm(j))
c     &        Dhold(j) = Dh_max*xKtm(j)
      enddo

      if (ndb.eq.1) Dhold(ndb) = Dhold(ndb+1)

      if (ipar.eq.1) write (nout,100) ndb,ndt,Dh_prescr
 100  format (1x,'Processing angular momentum transport from [',i4,',',
     &     i4,']',/,3x,'Type of Dh prescription : ',a8)


      return
      end






** ---------------------------------------------------------------------
*
      SUBROUTINE CALC_OMEGA (ttot,xomega,vitesse)
*
** ---------------------------------------------------------------------
* Calcule la vitesse angulaire de rotation a l'instant considere.
*
* Lois possibles :                      IDIFFVR :
*   Vitesse de rotation constante          1
*   Loi de puissance                       2
*   Loi de Skumanich                       3
*   Vitesse angulaire constante            4
*   Loi de Freinage en Omega^3             5
*   Moment angulaire spécifique constant   6
C   Periode de rotation constante          9
*
* Temps en annees, distances en km, vitesses en km/sec.
* diffcst est homogene a un temps en annees.
*
* OMEGA est appelee par DIFFUSION ou par DIF_OMEGA.
*
* Auteur : C.Charbonnel
* Modifs : S.Talon
*
* Adaptation: A.Palacios (14/03/2000)
*
* Derniere version : 2 octobre 2009
* ----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.diff'
      include 'evolcom.rot'
      include 'evolcom.teq'
      include 'evolcom.transp'

      integer psku1

      double precision psku2
      double precision xomega,vitesse,ttot
      double precision rr,xom
      double precision temps
      double precision xperiod

      common /omsurf/ xom

* Initialisations.
      if ((idiffty.ge.8.and.idiffty.le.13).or.idiffty.eq.15.or.
     &     idiffty.eq.17) then
         rr = rray(nmod)*1.d-5*rtot
         if (rayon0.lt.1.d-40.and.idiffvr.eq.3.or.idiffvr.eq.2) then
            rayon0 = rray(nmod)*1.d-5*rtot
         endif
c         if (rayon0.lt.1.d-40) rayon0 = rray(nmod)*1.d-5*rtot
      else
         if (rayon0.lt.1.d-40) rayon0 = rray(nmod)*1.d-5
         rr = rray(nmod)*1.d-5
      endif


      temps = ttot-breaktime
      if (idiffvr.eq.1) then
* Vitesse de rotation V constante
c         xomega = diffvr/rr
c         vitesse = diffvr
         xomega = vsurf/rr
         vitesse = vsurf
      else if (idiffvr.eq.2) then
* Loi de puissance
         if (temps.eq.0.d0) then
            xomega = diffvr/rayon0
         else
            xomega = diffvr/rayon0*(temps/diffcst)**idiffex
* idiffex ds evolcom.transp
         endif
         vitesse = rr*xomega
      else if (idiffvr.eq.3) then
* Loi de Skumanich
         psku1 = idiffex-1
         psku2 = 1.d0/dble(psku1)
         xomega = (diffvr)**psku1 + (diffcst*temps)
         xomega = (1.d0/xomega)**psku2
         vitesse = rr*xomega
      else if (idiffvr.eq.4) then
* Vitesse angulaire de rotation constante
cAP : Modification Ana March 2012 for PMS computations
c         xomega = breaktime
         omega_S = breaktime
         xomega = omega_S
         vitesse = xomega*rr

c$$$      else if (idiffvr.eq.9) then
c$$$C Periode de rotation constante
c$$$C diffcst est la periode de rotation en jours
c$$$         xperiod = diffcst*86400.d0
c$$$         vitesse = rr/xperiod
c$$$         xomega = vitesse/rr
      endif
      vsurf = vitesse
      xom = xomega
      
      return
      end
      

** ---------------------------------------------------------------------
*
      SUBROUTINE CONS_L (om,xmoment_tot,ivar,ndb,ndt)
*
** ---------------------------------------------------------------------
* Calcul du moment cinetique initial de l'etoile et verification
*   de sa conservation
*
* Pour ivar = 1 : Calcul de la quantite initiale
*      ivar = 2 : Verification de sa conservation
*
* Auteur: S.Talon
*
* Adaptation: A.Palacios (14/03/2000)
*
* Derniere version: 03/05/2000
* ----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.transp'

      integer idiffnuc
      integer i,i1,ivar
      integer ndt,ndb,ndtenv

      logical idiffcc,idiffcc0

      double precision om,xmoment_tot
      double precision xmoment,momenvel,xmoment_vieux
      double precision vomm,vdmom,omm,dmom
      double precision omega,vomega,k2conv,k2rad,vsurf,angradr,angconvr
      double precision rel
      double precision vom,vrray,vxpsi,dmomvieux
      double precision r,vr
      double precision momcore,momint,vmomenvel,vmomint

      common /rotvar/ omega(nsh),vomega(nsh),k2conv,k2rad,vsurf,
     &     angradr,angconvr
      common /momreel/ xmoment,momenvel,xmoment_vieux,momint
      common /diffcode/ idiffnuc,idiffcc,idiffcc0
      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)
      common /varrvr/ r(nsh),vr(nsh)
      common /envel/ ndtenv

      dimension om(nsh)


      xLtot = rtot*rtot*omega_S
      xmoment = 0.d0
      xmoment_vieux = 0.d0
      momenvel = 0.d0
      vmomint = 0.d0
      vmomenvel = 0.d0
      momcore = 0.d0
      momint = 0.d0
      do i = 2,nmod
         i1 = i-1
         omm = 0.5d0*(om(i)+om(i1))
         vomm = 0.5d0*(vom(i)+vom(i1))
         dmom = rraym2(i)*omm*dm(i1)
         vdmom = vrraym2(i)*vomm*dm(i1)
         xmoment = xmoment+dmom
         xmoment_vieux = xmoment_vieux+vdmom
         if (i.ge.ndtenv) then
            momenvel = momenvel+dmom
            vmomenvel = vmomenvel + vdmom
         else
            momint = momint+dmom
            vmomint = vmomint+vdmom
         endif
      enddo
      xmoment = xmoment*xLtot
      xmoment_vieux = xmoment_vieux*xLtot
      momenvel = momenvel*xLtot
      momcore = momcore*xLtot
      momint = momint*xLtot
      vmomenvel = vmomenvel*xLtot
      vmomint = vmomint*xLtot
      dmomvieux = (xmoment-xmoment_vieux)/xmoment_vieux
      if (ivar.eq.1) then
        xmoment_tot = xmoment
      else
        dmom = xmoment-xmoment_tot
        vsurf = om(nmod)*omega_S*rtot*1.d-5
        rel = dmom/xmoment
        write (nout,1000) xmoment,rel,dmomvieux
        write (nout,1001) vsurf,om(nmod),omega(2),omega_S,breaktime
        write (nout,1002) momenvel, momint
        write (nout,1003) vmomenvel, vmomint
      endif

 1000 format(/,' Total angular momentum: ',1pe12.5,/,
     &         ' Integrated conservation (dJ/dt)  : ',1pe12.5,/,
     &         ' Relative variation over the iteration: ',1pe12.5,/)
 1001 format(' Surface velocity =',1pe11.5,' (km/s)'/,' spin rate = ',
     &     1pe11.5,' (Hz), omegacore = ',E11.5,' omega_S = ',E11.5,
     &     ' breaktime =',E11.5)
 1002 format(' Current iteration',/
     $     ,' Angular momentum in the envelope : ',1pe12.5,/
     $     ,' Angular momentum in the core : ', 1pe12.5)
 1003 format(' Previous iteration',/
     $     ,' Angular momentum in the envelope : ',1pe12.5,/
     $     ,' Angular momentum in the core : ', 1pe12.5)

      return
      end



** ---------------------------------------------------------------------
*
      SUBROUTINE CENTOMEGA (ndb,neqj,nclc,pcol,dt,zj,xjmod)
*
** ---------------------------------------------------------------------
* Calcul des coefficients necessaires dans la methode de henyey
*   lors du calcul de l'evolution du moment cinetique ;
*     CONDITION LIMITE AU CENTRE
*
* Version generalisee incluant au moins 2 conditions limites :
* 1) CL sur xpsi
* 2) CL sur l'eq. d'evolution du moment cinetique
*
* CENTOMEGA est appelee par HEN_OMEGA
*
* Arguments:
*     neqj = nombre d'eq. du systeme a resoudre
*     nclc = nombre de CL au centre associees au syst. a resoudre
*     dt = dtom/rtot = pas de temps normalise au rayon
*                     (pour conditionnement de la matrice)
*
*     dtom = pas de temps de calcul
*     zj(neqj+nclc,neqj+nclc+1) = matrice pour le bloc limite central
*     xjmod(nsh,neqj) = matrice regroupant les vecteurs associes aux variables
*             independantes du syst.
*
*     NOTA BENE: Par defaut, on aura xjmod(i,pcol(1)) = om(i)
*                                    xjmod(i,2) = aux(i)
*                                    xjmod(i,3) = ur(i)
*                                    xjmod(i,4) = xpsi(i)
*                TOUTE NOUVELLE VARIABLE INDEPENDANTE SERA RAJOUTEE A LA SUITE
*                exemple: xjmod(i,5) = xlambda(i) si zgradmu = T
*
* Auteur: A.Palacios
*
* ----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.therm'
      include 'evolcom.transp'

      integer j,k,j1
      integer neqj,n1,n2,nclc,ndb
      integer pcol

      double precision vom,vrray,vxpsi,xmom_new,xmom_old,vomm,tamp1
      double precision zj,xjmod,om,ur,dt

      dimension zj(neqj+nclc,neqj+nclc+1),xjmod(nsh,neqj)
      dimension pcol(neqj)
      dimension om(nsh),ur(nsh)

      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)

*-----------------------
***   Initialisations
*-----------------------

      n1 = neqj+nclc-1
      n2 = neqj+nclc
      do j = neqj+1,neqj+nclc
         do k = 1,neqj+nclc+1
            zj(j,k) = 0.d0
         enddo
      enddo

      if (ndb.eq.1) then
c..    boundary condition : Ur(1) = 0
         zj(n1,neqj+nclc+1) = xjmod(1,pcol(3))
         zj(n1,nclc+pcol(3)) = 1.d0
      else
c..   Condition on the moment
         xmom_new = 0.d0
         xmom_old = 0.d0
         tamp1 = 0.2d0*rray(ndb)**4*rdensm(ndb)
         tamp1 = tamp1*xLloc         
         om(ndb) = xjmod(ndb,pcol(1))
         ur(ndb) = xjmod(ndb,pcol(3))
         do j = 2,ndb
            j1 = j-1
            vomm = 0.5d0*(vom(j)+vom(j1))
            xmom_new = xmom_new+rraym2(j)*dm(j1)
            xmom_old = xmom_old+vrraym2(j)*dm(j1)*vomm
         enddo
         zj(n1,neqj+nclc+1) = om(ndb)*ur(ndb)*ur_S*dt-(xmom_new
     &        *om(ndb)-xmom_old)/tamp1
         zj(n1,nclc+pcol(1)) = ur(ndb)*ur_S*dt-xmom_new/tamp1
         zj(n1,nclc+pcol(3)) = om(ndb)*ur_S*dt            
      endif

c..    boundary condition : xpsi(1) = 0 = xjmod(1,pcol(4))
      zj(n2,neqj+nclc+1) = 0.d0
      zj(n2,nclc+pcol(4)) = 1.d0

c..NEW
      zj(n1,neqj+nclc+1) = -zj(n1,neqj+nclc+1)
      zj(n2,neqj+nclc+1) = -zj(n2,neqj+nclc+1)

      return
      end



** ---------------------------------------------------------------------
*
      SUBROUTINE INTOMEGA (j,j1,iterh,ndb,ndt,neqj,ncls,nclc,pcol,dt,
     &     umat,xjmod,eqj,gj,depotgi)
*
** ---------------------------------------------------------------------
* Calcul des elements matriciels des blocs ventraux pour Newton-Raphson
*
*
* INTOMEGA est appelee par HEN_OMEGA
C CC : Ancienne SUBROUTINE GI_OMEGA
*
* Arguments:
*     j = numero de la couche courante
*     j1 = j-1
*     iterh = numero de l'iteration Newton-Raphson courante
*     ndt = borne superieure de la zone de melange
*     neqj = nombre d'eq. du syst. a resoudre pour calculer l'evolution
*            du moment cinetique
*     ncls = nombre de CL en surface associees au syst.
*     nclc = nombre de CL au centre associees au syst.
*     dt = dtom/rtot = pas de temps normalise au rayon total
*                     (necessaire pour conditionnement de la matrice)
*     umat(nsh,neqj,nclc+1) = matrice U,V,W permettant d'exprimer les neqj-nclc
*                      premieres variables en fonction des nclc variables
*                      restantes (lie a l'ecriture du bloc de surface)
*                      (c.f. pp 103-104 these A. Palacios)
*     xjmod(nsh,neqj) = matrice regroupant les vecteurs associes aux variables
*             independantes du syst.
***==========================================================================*
*     !! ATTENTION!!                                                         *
*     La matrice jacobienne associee au systeme a 5 equations peut etre      *
*     singuliere selon l'ordonnancement adopte.                              *
*     Pour eviter toute irregularite, on definit des pointeurs de colonnes   *
*     (pcol(neqj)) et des pointeurs de lignes (elin(neqj)) qui permettent    *
*     d'ordonner la matrice de la facon la plus appropriee.                  *
***==========================================================================*
*
*     pcol(neqj) = pointeur de colonnes pour ordonnancement des matrices
*     elin(neqj) = pointeur de lignes pour ordonnancement des matrices
*
*     NOTA BENE: Par defaut, on aura xjmod(i,pcol(1)) = om(i)
*                                    xjmod(i,pcol(2)) = aux(i)
*                                    xjmod(i,pcol(3)) = ur(i)
*                                    xjmod(i,pcol(4)) = xpsi(i)
*                TOUTE NOUVELLE VARIABLE INDEPENDANTE SERA RAJOUTEE A LA SUITE
*                exemple: xjmod(i,pcol(5)) = xlambda(i) si zgradmu = T
*
*     eqj(neqj,2*neqj+1) = elements de la matrice jacobienne Gi 
*                          (c.f. pp 102 these A. Palacios)
*                          derivees exprimees sur 2 couches :
*                              neqj premieres colonnes = derivees / aux variables en j
*                              neqj colonnes suivantes = derivees / aux variables en j1
*                              colonne 2*neqj+1 = equation a resoudre
*
*     Organisation des lignes:       eqj(1,x) = eq. pour xpsi
*                                    eqj(2,x) = eq. pour aux
*                                    eqj(3,x) = eq. pour Ur
*                                    eqj(4,x) = eq. d'evolution du moment cinetique
*                TOUTE NOUVELLE EQUATION SERA RAJOUTEE A LA SUITE
*                exemple: eqj(5,x) = eq. pour lambda si zgradmu = T
*
*     Organisation des colonnes:
* derivee / a  om(j)           aux(j)         ur(j)           xpsi(j)    xlambda(j)
*           eqj(x,pcol(1)) eqj(x,pcol(2)) eqj(x,pcol(3))  eqj(x,pcol(4))  eqj(x,pcol(5))
* derivee / a  om(j1)          aux(j1)             ur(j1)             xpsi(j1)          xlambda(j1)
*    eqj(x,neqj+pcol(1)) eqj(x,neqj+pcol(2)) eqj(x,neqj+pcol(3)) eqj(x,neqj+pcol(4)) eqj(x,neqj+pcol(5))
*
*     gj(neqj,neqj+nclc+1) = matrice "pseudo-jacobienne" effectivement inversee dans HEN_OMEGA
*                           (c.f. expression 5.53 pp 105 these A. Palacios)
*
* Dans cette version, on calcule au minimum la diffusion par viscosite
*   moleculaire
*
* Le vecteur nturb contient les coordonnees des couches pour lesquelles
*   xnuvv > xnum
*
*
* Auteur : A.Palacios (05/2004)
C Modifs CC ondes (22/11/06)
* ----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.transp'
      include 'evolcom.surf'


      integer j,j1,k,l,mm,i
      integer iterh
      integer ndt,ndtenv,ndb
      integer neqj,neqj1,ncls,nclc
      integer pcol,elin
      integer klcore,klenv,klpulse

      double precision dtom


      double precision om,aux,ur,xpsi,xlambda,dt,xjmod
      double precision vom,vrray,vxpsi
      double precision dift,dife,dj,dj1,diff
      double precision rhmoy,epsmoy,epsmoym
      double precision vomm,xpsim,auxm,urm,omm,omm2
      double precision tamp,tamp0,tamp1,tamp1a,tamp1b,tamp2,tamp3,tamp4
      double precision tamp5a,tamp5b,tamp5
      double precision aalpha,bbeta,divis,raprho,rapeps
      double precision dom,cte1,cte2
      double precision xnuvvm,xnuvv_j,xnuvv_j1,dnuvvo_j,dnuvvp_j,
     &     dnuvvl_j,dnuvvo_j1,dnuvvp_j1,dnuvvl_j1
      double precision umat,gj

      double precision vxpsim,fepsi
      double precision tamp2a,tamp2b,tamp6a,tamp6
      double precision Dhold,Dholdm,Rec
      double precision xnumol,xnurad,nturb
      double precision t,vt
      double precision eqj
      double precision abmurj,abmuj,xlambdam
      double precision xlambdaold,dlam
      double precision xxnorm
      double precision rmstar

      double precision omega_sat
      double precision totm
      double precision depotgi

      double precision dshear

      common /dturbulente/ nturb(0:nsh)
      common /difcirc/ dift(nsh),dife(nsh)
      common /moy/ rhmoy(nsh),epsmoy(nsh)
      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)

      common /calcDh/ Dhold(nsh)
      common /diffvisc/ xnumol(nsh),xnurad(nsh)
      common /vartvt/ t(nsh),vt(nsh)
      common /gradmuj/ abmurj(nsh),abmuj(nsh)
      common /lambdaold/ xlambdaold(nsh)
      common /envel/ ndtenv
      common /overshoot/ klcore,klenv,klpulse
      common /chemdiffshear/ dshear(nsh)
      
      dimension umat(nsh,neqj,nclc+1),gj(neqj,neqj+nclc+1),
     &     xjmod(nsh,neqj),eqj(neqj,2*neqj+1)
      dimension om(nsh),aux(nsh),ur(nsh),xpsi(nsh),xlambda(nsh)
      dimension pcol(neqj),elin(neqj)

      common /mass/ totm

*------------------------------------
***   Initialisations
*------------------------------------

      neqj1 = 2*neqj+1

      do k = j1,j
         om(k) = xjmod(k,pcol(1))
         aux(k) = xjmod(k,pcol(2))
         ur(k) = xjmod(k,pcol(3))
         xpsi(k) = xjmod(k,pcol(4))
         if (zgradmu) xlambda(k) = xjmod(k,pcol(5))
      enddo
      om(1) = xjmod(1,pcol(1))

      vomm = 0.5d0*(vom(j)+vom(j1))
      xpsim = 0.5d0*(xpsi(j)+xpsi(j1))
      urm = 0.5d0*(ur(j)+ur(j1))
      auxm = 0.5d0*(aux(j)+aux(j1))
      omm = 0.5d0*(om(j)+om(j1))
      Dholdm = 0.5d0*(Dhold(j)+Dhold(j1))

      tamp = hnenn(j)

      do i = 1,neqj
         do k = 1,neqj1
            eqj(i,k) = 0.d0
         enddo
         do k = 1,neqj+nclc+1
            gj(i,k) = 0.d0
         enddo
      enddo

      if (zgradmu) then
         elin(1) = 1
         elin(2) = 2
         elin(3) = 3
         elin(4) = 4
         elin(5) = 5
      else
         elin(1) = 1
         elin(2) = 2
         elin(3) = 3
         elin(4) = 4
      endif


*--------------------------------------------------------
***   Definition de xpsi ; eq. (5.24) these A. Palacios
***   Derniere ligne du bloc matriciel si neqj = 5
***   Premiere ligne du bloc matriciel si newj = 4
*--------------------------------------------------------

      tamp1 = hnenn(j)*pw23*rrog(j)
      tamp2 = tamp1*(om(j1)-om(j))
      tamp3 = deltaKS(j1)*xpsim/omm

      eqj(elin(1),neqj1) = tamp2+tamp3

      eqj(elin(1),pcol(1)) = -tamp1-0.5d0*tamp3/omm
      eqj(elin(1),neqj+pcol(1)) = tamp1-0.5d0*tamp3/omm
      eqj(elin(1),pcol(4)) = 0.5d0*deltaKS(j1)/omm
      eqj(elin(1),neqj+pcol(4)) = 0.5d0*deltaKS(j1)/omm



*--------------------------------------------------------
***   Definition de aux ; eq. (5.25) these A. Palacios
***   Seconde ligne du bloc matriciel si neqj = 4 ou 5
*--------------------------------------------------------

      tamp1 = hhtm(j)*hnenn(j)
      tamp2 = 1.d0-deltaKS(j1)+khit(j1)

      eqj(elin(2),neqj1) = auxm-tamp1*(xpsi(j1)-xpsi(j))+tamp2*xpsim

      eqj(elin(2),pcol(2)) = 0.5d0
      eqj(elin(2),neqj+pcol(2)) = 0.5d0
      eqj(elin(2),pcol(4)) = tamp1+0.5d0*tamp2
      eqj(elin(2),neqj+pcol(4)) = -tamp1+0.5d0*tamp2

*--------------------------------------------------------
***   Definition de la vitesse de circulation meridienne
***   eq. (5.26) these A. Palacios
***   Troisieme ligne du bloc matriciel si neqj = 4
***   Quatrieme ligne du bloc matriciel si neqj = 5
*--------------------------------------------------------

      omm2 = omm*omm
***   Nouveau cas avec expression de U simplifiée pour test 
***   dans le cas des modèles d'étoiles massives.
***   Utilisation de l'expression de Ur donnée dans la routine CI

      if (idiffty.eq.17) then
         rmstar = rmassm(j)*(1.d0-omm2*omega_S**2/(pim2*g*rhmoy(j)))
         if ((crz(j).lt.2.or.j.gt.novlim(klenv,3)).and.j.ne.ndt) then
c     if ((crz(j).lt.0.or.crz(j1).lt.0).and.j.ne.ndt) then
            aalpha = 0.d0
         else
            divis = -gravm(j)*rmstar*cp(j1)*t(j1)*rkonvm(j)
            aalpha = xlumlm(j)*p(j1)*xnorm/divis
         endif
         raprho = rro(j1)/rhmoy(j)
         bbeta = 2.d0*(pw43-raprho)
         tamp4 = bbeta*rraym(j)*omm/gravm(j) ! gtilde/gbar
         tamp1a = omm2*omega_S*omega_S*2.384989d6/rro(j1)
         tamp1 = 1.d0-tamp1a

         eqj(elin(3),neqj1) = -rro(j1)*urm*ur_S+aalpha*(tamp4*omm*tamp1)
         
         eqj(elin(3),pcol(1)) = aalpha*0.5d0*tamp4*2.d0*(tamp1-tamp1a)
         eqj(elin(3),neqj+pcol(1)) = eqj(elin(3),pcol(1))
         eqj(elin(3),pcol(3)) = -0.5d0*rro(j1)*ur_S
         eqj(elin(3),neqj+pcol(3)) = eqj(elin(3),pcol(3))        
                   
c         print *,'intomega',tamp4,tamp1,aalpha,divis,xnorm,p(j1),
c     &        xlumlm(j)
      else
         if (idiffty.eq.13) then

            vxpsim = 0.5d0*(vxpsi(j)+vxpsi(j1))
            rmstar = rmassm(j)*(1.d0-omm2*omega_S**2/(pim2*g*rhmoy(j)))
            tamp6a = rmstar*cp(j1)*t(j1)/(dt*rtot*xlumlm(j))
            tamp6 = tamp6a*(xpsim-vxpsim)
c     if (crz(j).lt.0.and.j.ne.ndt) then
            if ((crz(j).lt.2.or.j.gt.novlim(klenv,3)).and.j.ne.ndt) then
c     if ((crz(j).lt.0.or.crz(j1).lt.0).and.j.ne.ndt) then
               aalpha = 0.d0
            else
               divis = -gravm(j)*rmstar*cp(j1)*t(j1)*rkonvLm(j)
               aalpha = xlumlm(j)*p(j1)*xnorm/divis
            endif
            
            raprho = rro(j1)/rhmoy(j)
            if (abs(enucl(j1)).gt.1.d-40) then
               epsmoym = 0.5d0*(epsmoy(j)+epsmoy(j1))
               rapeps = (enucl(j1)+egrav(j1))/epsmoym
               fepsi = enucl(j1)/(enucl(j1)+egrav(j1))
            else
               rapeps = 0.d0
               fepsi = 0.d0
            endif
         else
            tamp6 = 0.d0
            tamp6a = 0.d0
            if ((crz(j).lt.2.or.j.gt.novlim(klenv,3)).and.j.ne.ndt) then
c     if ((crz(j).lt.0.or.crz(j1).lt.0).and.j.ne.ndt) then
               aalpha = 0.d0
            else
               divis = -gravm(j)*rmstar*cp(j1)*t(j1)*rkonvm(j)
               aalpha = xlumlm(j)*p(j1)*xnorm/divis
            endif
            
            raprho = rro(j1)/rhmoy(j)
            if (abs(enucl(j1)).gt.1.d-40) then
               epsmoym = 0.5d0*(epsmoy(j)+epsmoy(j1))
               rapeps = enucl(j1)/epsmoym
               fepsi = 1.d0
            else
               rapeps = 0.d0
               fepsi = 0.d0
            endif
         endif
         
         bbeta = 2.d0*(pw43-raprho)
         
         tamp1a = omm2*omega_S*omega_S*2.384989d6/rro(j1)
C     correction par les valeur centrales 
c..   A VERIFIER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         tamp1a = tamp1a-om(1)*om(1)*omega_S*omega_S*2.384989d6/rro(1)
         
         tamp1b = pw23/raprho*deltaKS(j1)*omega_S*omega_S
         tamp1 = 1.d0-tamp1a+tamp1b*xpsim-rapeps
         tamp2a = rapeps*(fepsi*(epsit(j1)-deltaKS(j1))+deltaKS(j1))
         tamp2b = (2.d0*hhtm(j)/rraym(j)*(1.d0+Dholdm/xKt(j1))-pw23*
     &        deltaKS(j1))/raprho
         tamp2 = tamp2a-tamp2b
         tamp3 = (pw13/raprho)*rraym(j)*tamp
         tamp4 = bbeta*rraym(j)*omm/gravm(j) ! gtilde/gbar
         
         eqj(elin(3),neqj1) = -rro(j1)*urm*ur_S+aalpha*(tamp4*omm*tamp1
     &        +tamp3*(aux(j1)-aux(j))+rapeps*auxm+tamp2*xpsim-tamp6)
         
         eqj(elin(3),pcol(1)) = aalpha*0.5d0*tamp4*2.d0*(tamp1-tamp1a)
         eqj(elin(3),neqj+pcol(1)) = eqj(elin(3),pcol(1))
         eqj(elin(3),pcol(2)) = (-tamp3+0.5d0*rapeps)*aalpha
         eqj(elin(3),neqj+pcol(2)) = (tamp3+0.5d0*rapeps)*aalpha
         eqj(elin(3),pcol(3)) = -0.5d0*rro(j1)*ur_S
         eqj(elin(3),neqj+pcol(3)) = eqj(elin(3),pcol(3))
         eqj(elin(3),pcol(4)) = 0.5d0*aalpha*(tamp2-tamp6a+tamp4*omm*
     &        tamp1b)
         eqj(elin(3),neqj+pcol(4)) = eqj(elin(3),pcol(4))
         
      endif

*--------------------------------------------------------
***   Equation d'evolution du moment cinetique
***   eq. (5.27) these A. Palacios
***   Quatrieme ligne du bloc matriciel si neqj = 4
***   Premiere ligne du bloc matriciel si neqj = 5
*--------------------------------------------------------

*     Terme en d/dt

      dom = rraym2(j)*omm-vrraym2(j)*vomm
      eqj(elin(4),neqj1) = dom
      eqj(elin(4),pcol(1)) = rraym2(j)*0.5d0
      eqj(elin(4),neqj+pcol(1)) = rraym2(j)*0.5d0

*     Partie advective de l'equation d'advection/diffusion

      tamp1 = rdensm(j1)*rray(j1)**4
      cte1 = tamp1*om(j1)*ur(j1)
      tamp2 = rdensm(j)*rray(j)**4
      cte2 = tamp2*om(j)*ur(j)
      tamp3 = -dt*ur_S*rtot**3*pim4/(5.d0*dm(j1))

      eqj(elin(4),neqj1) = eqj(elin(4),neqj1)-tamp3*(cte1-cte2)
     & -depotgi/dm(j1)

      eqj(elin(4),pcol(1)) = eqj(elin(4),pcol(1))+tamp3*tamp2*ur(j)
      eqj(elin(4),neqj+pcol(1)) = eqj(elin(4),neqj+pcol(1))-tamp3*tamp1
     &     *ur(j1)
      eqj(elin(4),pcol(3)) = tamp3*tamp2*om(j)
      eqj(elin(4),neqj+pcol(3))= -tamp3*tamp1*om(j1)

*     Partie diffusive de l'equation d'advection/diffusion
*     Attention a n'envoyer aucun flux d'une region turbulente vers
*     une region non-turbulente


      xnuvvm = 0.5d0*(xnuvv(j)+xnuvv(j1))
      Rec = 10.d0*xnum(j)
c      Rec = (3.d0+pw13)*xnum(j)
c.. A VERIFIER !!!!!!!!!!!!!!!!!!!!!!!!!!
c      Rec = xnum(j)
cori      Rec = xnumol(j1)
c      if (j.eq.ndt) nturb(j+1) = nturb(j)

      xnuvv_j = 0.d0
      dnuvvo_j = 0.d0
      dnuvvp_j = 0.d0
      dnuvvl_j = 0.d0
      xnuvv_j1 = 0.d0
      dnuvvo_j1 = 0.d0
      dnuvvp_j1 = 0.d0
      dnuvvl_j1 = 0.d0
      if (crz(j).gt.1.and.crz(j).ne.3.and.j.lt.ndt.and.
     &     j.lt.novlim(klenv,3)) then
c      if (crz(j).gt.0.and.j.lt.ndt) then
         if (iterh.gt.7) then

c          print *,'iterh',iterh,crz(j),abs(crz(j)),ndt,j,novlim(klenv,3)
c          print *,'j',j,'nturb(j+1)*nturb(j)',nturb(j+1)*nturb(j)
c     $         ,rkonvm(j),xnuvvm,Rec
            xnuvv_j = xnum(j)+nturb(j+1)*nturb(j)*xnuvv(j)
     & +Dondes(j)
            dnuvvo_j = nturb(j+1)*nturb(j)*dnuvvo(j)
            dnuvvp_j = nturb(j+1)*nturb(j)*dnuvvp(j)
            dnuvvl_j = nturb(j+1)*nturb(j)*dnuvvl(j)
         else if (Dv_prescr.ne.'Pr16') then
            if (xnuvvm.le.Rec.and.(rkonvm(j).lt.0.d0)) then
               nturb(j) = 0.d0
               xnuvv_j = xnum(j)+Dondes(j)
               dnuvvo_j = 0.d0
               dnuvvp_j = 0.d0
               dnuvvl_j = 0.d0
            else
               nturb(j) = 1.d0
               xnuvv_j = xnum(j)+(nturb(j+1)*xnuvv(j))+Dondes(j)
               dnuvvo_j = nturb(j+1)*dnuvvo(j)
               dnuvvp_j = nturb(j+1)*dnuvvp(j)
               dnuvvl_j = nturb(j+1)*dnuvvl(j)
            endif
         else
            nturb(j) = 1.d0
            xnuvv_j = xnum(j)+xnuvv(j)+Dondes(j)
            dnuvvo_j = dnuvvo(j)
            dnuvvp_j = dnuvvp(j)
            dnuvvl_j = dnuvvl(j)           
         endif
         dj = xnuvv_j
      endif
c      if (crz(j1).gt.1.and.crz(j).ne.3.and.j1.lt.
c     $     novlim(klenv,3)) then
      if (crz(j1).gt.0) then
         if (iterh.gt.7) then
            xnuvv_j1 = xnum(j1)+nturb(j1)*nturb(j)*xnuvv(j1)
     & +Dondes(j1)
            dnuvvo_j1 = nturb(j1)*nturb(j)*dnuvvo(j1)
            dnuvvp_j1 = nturb(j1)*nturb(j)*dnuvvp(j1)
            dnuvvl_j1 = nturb(j1)*nturb(j)*dnuvvl(j1)
         else if (Dv_prescr.ne.'Pr16') then
            if ((xnuvvm.le.Rec).and.(rkonvm(j).lt.0.d0)) then
               nturb(j1) = 0.d0
               xnuvv_j1 = xnum(j1)+Dondes(j1)
               dnuvvo_j1 = 0.d0
               dnuvvp_j1 = 0.d0
               dnuvvl_j1 = 0.d0
            else
               nturb(j) = 1.d0
               xnuvv_j1 = xnum(j1)+(nturb(j1)*xnuvv(j1))
     &              +Dondes(j1)
               dnuvvo_j1 = nturb(j1)*dnuvvo(j1)
               dnuvvp_j1 = nturb(j1)*dnuvvp(j1)
               dnuvvl_j1 = nturb(j1)*dnuvvl(j1)
            endif
         else
            nturb(j1) = 1.d0
            xnuvv_j1 = xnum(j1)+(nturb(j1)*xnuvv(j1))
     &           +Dondes(j1)
            dnuvvo_j1 = nturb(j1)*dnuvvo(j1)
            dnuvvp_j1 = nturb(j1)*dnuvvp(j1)
            dnuvvl_j1 = nturb(j1)*dnuvvl(j1)
         endif
         dj1 = xnuvv_j1
      endif
      
ccc   NL... additional viscosity ... nuadd
c      print *, 'Additional viscosity : Eggenberger et al. 2016, 2019'
      if (j.lt.ndt)then
       dj1=xnuvv_j1+3.5d4 
       dj=xnuvv_j+3.5d4
      endif 

      if (crz(j1).le.1.or.j1.ge.novlim(klenv,4)) then
c      if (crz(j1).lt.0) then
         xnuvv_j1 = 0.d0
         dnuvvo_j1 = 0.d0
         dnuvvp_j1 = 0.d0
         dnuvvl_j1 = 0.d0
c         if (crz(j1).le.0) then
            dj1 = Dconv(j1)
c         else if (crz(j1).eq.1.and.(lover.eq.33)) then                      ! modif AP TD 11/2019
c            dj1 = Dbar(j1)+3.d4
c         else if (crz(j1).eq.1.and.(lover.ge.34.and.lover.le.36)) then
c            dj1 = Dkyle(j1)+3.d4
c            dj1 = xnuvv_j1+3.d4
c         else if (crz(j1).eq.1.and.novopt.eq.1) then
c            dj1 = Dconv(j1)
c         endif                                                            ! modif AP TD 11/2019
      endif
      if (crz(j).le.1.or.j.ge.novlim(klenv,4)) then
c      if (crz(j).lt.0) then
         xnuvv_j = 0.d0
         dnuvvo_j = 0.d0
         dnuvvp_j = 0.d0
         dnuvvl_j = 0.d0
c         if (crz(j).le.0) then
            dj = Dconv(j)
c         else if (crz(j).eq.1.and.(lover.eq.33)) then ! modif AP TD 11/2019
c            dj = Dbar(j)+3.d4
c         else if (crz(j).eq.1.and.(lover.ge.34.and.lover.le.36)) then
c     dj = Dkyle(j)+3.d4
c            dj = xnuvv_j+3.d4
c         else if (crz(j).eq.1.and.novopt.eq.1) then
c            dj = Dconv(j)
c         endif                  ! modif AP TD 11/2019     
      endif
      
c$$$      write(552,'(2(1x,i4),4(1x,1pe11.4))'),j,crz(j),Dconv(j)
c$$$     $     ,Dkyle(j),dj,dj1

      if (j.eq.ndt) dj = Dconv(ndt)
c      dnuvv_j = 0.d0
      dnuvvo_j = 0.d0
      dnuvvp_j = 0.d0
      dnuvvl_j = 0.d0
      dnuvvo_j1 = 0.d0
      dnuvvp_j1 = 0.d0
      dnuvvl_j1 = 0.d0


      cte1 = rdensm(j1)*rray(j1)*rray(j1)*grav(j1)
      cte2 = rdensm(j)*rray(j)*rray(j)*grav(j)
      tamp0 = -1.5d0*pim4*dt*rtot**2/dm(j1)
      tamp1a = cte1*deltaKSm(j1)*xpsi(j1)/om(j1)*dj1
      tamp1b = cte2*deltaKSm(j)*xpsi(j)/om(j)*dj

      eqj(elin(4),neqj1) = eqj(elin(4),neqj1)+tamp0*(tamp1a-tamp1b)
     
      eqj(elin(4),pcol(1)) = eqj(elin(4),pcol(1))+tamp0*cte2*
     &     deltaKSm(j)*xpsi(j)/om(j)*(dj/om(j)-dnuvvo_j)
      eqj(elin(4),neqj+pcol(1)) = eqj(elin(4),neqj+pcol(1))-tamp0*
     &     cte1*deltaKSm(j1)*xpsi(j1)/om(j1)*(dj1/om(j1)-dnuvvo_j1)

      eqj(elin(4),pcol(4)) = eqj(elin(4),pcol(4))-tamp0*
     &     cte2*deltaKSm(j)*(dj+xpsi(j)*dnuvvp_j)/om(j)
      eqj(elin(4),neqj+pcol(4)) = eqj(elin(4),neqj+pcol(4))+tamp0*
     &     cte1*deltaKSm(j1)*(dj1+xpsi(j1)*dnuvvp_j1)/om(j1)

* Couple en surface
      if (j.gt.ndtenv.and.idiffvr.ne.0) then
c        fonction amortissement
c         tamp0 = 0.15d0*(rray(nmod)-rray(ndtenv))
c         tamp1 = -4.d0*((rray(j)-rray(nmod))/tamp0)**2
c         tamp2 = 2.d0/(sqrt(pi)*tamp0)*exp(tamp1)
c         tamp2 = tamp2/rray(j)**2

c        autre fonction
c         tamp2 = abs((rray(j)-rray(j1))/(rray(ndtenv)**3-rray(nmod)**3))
c         tamp2 = (rray(j)-rray(ndtenv))/(rray(nmod)-rray(ndtenv))
c         tamp2 = 2.d5*(rray(j)-rray(ndtenv))/
c     &        (rray(nmod)-rray(ndtenv))**2
c         tamp2 = 1.d0

      endif

C xnuvv : coefficient de diffusion pour le moment
C         le coefficient de diffusion pour la chimie
C         ne tient pas compte de la viscosite radiative qui transporte
C         le moment mais non les especes
C        Ici on utilise xnuvv et pas dj car pour la diffusion des especes
C        chimiques dans DIFFUSION, on somme Dhd et Dconv.

C Modif CC ondes (12/04/07) -->
C Dans le cas ou le calcul est effectue avec diffusion atomique,
C xnumol est deja pris en compte 
C et ne doit pas etre ajoute dans le calcul de Dhd

C      if (microdiffus) then 
C          if (j.ne.ndt) then
C             Dhd(j) = nturb(j+1)*nturb(j)*xnuvv(j)
C          else
C             Dhd(j) = nturb(j)*xnuvv(j)
C          endif
C          dift(j) = Dhd(j)
C          if (iterh.gt.1) then
C             Dhd(j1) = nturb(j1)*nturb(j)*xnuvv(j1)
C             dift(j1)= Dhd(j1)
C          endif
C      else
C          if (j.ne.ndt) then
C             Dhd(j) = xnumol(j)+nturb(j+1)*nturb(j)*xnuvv(j)
C          else
C             Dhd(j) = xnumol(j)+nturb(j)*xnuvv(j)
C          endif
C          dift(j) = Dhd(j)
C          if (iterh.gt.1) then
C             Dhd(j1) = xnumol(j1)+nturb(j1)*nturb(j)*xnuvv(j1)
C             dift(j1)= Dhd(j1)
C          endif
C      endif
      if (microdiffus) then 
         Dhd(j) = dshear(j)
         dift(j) = Dhd(j)
         if (iterh.gt.1) then
            Dhd(j1) = dshear(j1)
            dift(j1)= Dhd(j1)
         endif
      else
         Dhd(j) = xnumol(j)+dshear(j)
         dift(j) = Dhd(j)
         if (iterh.gt.1) then
            Dhd(j1) = xnumol(j1)+dshear(j1)
            dift(j1)= Dhd(j1)
         endif
      endif
C <--

c      Dhd(j) = xnuvv_j+Dsc(j)
c      dift(j) = xnuvv_j+Dsc(j)
c      Dhd(j1) = xnuvv_j1+Dsc(j1)
c      dift(j1) = xnuvv_j1+Dsc(j1)

*---------------------------------------------------------------------
***   Equation d'evolution de lambda ; eq. (5.28) these A. Palacios
***   Troisieme ligne du bloc matriciel si neqj = 5
***   N'existe pas si neqj < 5
*---------------------------------------------------------------------

      if (zgradmu) then

         xlambdam = 0.5d0*(xlambda(j)+xlambda(j1))
         xxnorm = xnorm*gs
         dlam = xlambda(j1) - xlambdaold(j1)
         l = max(j1,2)
         if (crz(l).le.1.or.l.ge.novlim(klenv,3)) then
            if (crz(l).eq.1.and.(lover.ge.33.and.lover.le.36)) then
               print *,'zgradmu jtransp malheureux !!'
               stop !modif AP TD 11/2019
c            diff = max(Dhold(l),Dkyle(l),Dbar(l))
c         else if (crz(l).lt.1) then
            diff = max(Dhold(l),Dconv(l))
         endif
         else
            diff = Dhold(l)
         endif
         
         tamp1 = diff/(rtot*ur_S)*6.d0/(rray(l)*rray(l))
 
         eqj(elin(5),neqj1) = dlam*xxnorm/ur_S + dt*(ur(j1)*rtot
     &        *abmurj(j1)+tamp1*xlambda(j1)*xxnorm)

         eqj(elin(5),neqj+pcol(5)) =  xxnorm/ur_S + tamp1*dt*xxnorm
         eqj(elin(5),neqj+pcol(3)) = dt*abmurj(j1)*rtot

c..    Ajout des derivees par rapport a lambda pour les equations concernees

* variable PSI
         tamp1 = -phiKS(j1)*xlambdam/omm
         eqj(elin(1),neqj1) = eqj(elin(1),neqj1)+tamp1

         eqj(elin(1),pcol(1)) = eqj(elin(1),pcol(1))-0.5d0*tamp1/omm
         eqj(elin(1),neqj+pcol(1)) = eqj(elin(1),neqj+pcol(1))-
     &        0.5d0*tamp1/omm
         eqj(elin(1),pcol(5)) = -0.5d0*phiKS(j1)/omm
         eqj(elin(1),neqj+pcol(5)) = eqj(elin(1),pcol(5))

* variable AUX
         tamp3 = phiKS(j1)+khi_mu(j1)

         eqj(elin(2),neqj1) = eqj(elin(2),neqj1) + tamp3*xlambdam

         eqj(elin(2),pcol(5)) = 0.5d0*tamp3
         eqj(elin(2),neqj+pcol(5)) = eqj(elin(2),pcol(5))


* variable Ur
         tamp5a = -pw23*phiKS(j1)/raprho
         tamp5b = rapeps*(fepsi*(epsimu(j1)+phiKS(j1))-phiKS(j1))
         tamp5 = (tamp5a+tamp5b)*xlambdam
         tamp6 = -pw23*phiKS(j1)/raprho*omega_S*omega_S

         eqj(elin(3),neqj1) = eqj(elin(3),neqj1) + aalpha*(tamp4*omm*
     &        tamp6*xlambdam+tamp5)

         eqj(elin(3),pcol(1)) = eqj(elin(3),pcol(1)) + aalpha*tamp4
     &        *tamp6*xlambdam
         eqj(elin(3),neqj+pcol(1)) = eqj(elin(3),neqj+pcol(1)) + aalpha*
     &        tamp4*tamp6*xlambdam
         eqj(elin(3),pcol(5)) =  0.5d0*aalpha*(tamp5a+tamp5b+tamp4*omm
     &        *tamp6)
         eqj(elin(3),neqj+pcol(5)) = eqj(elin(3),pcol(5))

* variable OMEGA
         tamp2a = -cte1*phiKSm(j1)*xlambda(j1)/om(j1)*dj1
         tamp2b = -cte2*phiKSm(j)*xlambda(j)/om(j)*dj

         eqj(elin(4),neqj1) = eqj(elin(4),neqj1)+tamp0*(tamp2a-tamp2b)

         eqj(elin(4),pcol(1)) = eqj(elin(4),pcol(1))+tamp0*cte2*
     &        phiKSm(j)*xlambda(j)/om(j)*(dnuvvo_j-dj/om(j))
         eqj(elin(4),neqj+pcol(1)) = eqj(elin(4),neqj+pcol(1))-tamp0*
     &        cte1*phiKSm(j1)*xlambda(j1)/om(j1)*(dnuvvo_j1-dj1/om(j1))
         eqj(elin(4),pcol(5)) = tamp0*cte2*phiKSm(j)*
     &        (dj+xlambda(j)*dnuvvl_j)/om(j)
         eqj(elin(4),neqj+pcol(5)) = -tamp0*cte1*phiKSm(j1)*
     &        (dj1+xlambda(j1)*dnuvvl_j1)/om(j1)
      endif


*------------------------------------------------------------------------------
*** Remplissage de la matrice "pseudo-jacobienne" pour inversion dans HEN_OMEGA
*------------------------------------------------------------------------------

      do k = 1,neqj
         do l = 1,nclc
            gj(k,l) = eqj(k,neqj-nclc+l)
         enddo
         gj(k,neqj+nclc+1) = -eqj(k,neqj1)
         do mm = 1,nclc
            do l = 1,ncls
               gj(k,mm) = gj(k,mm)+umat(j,l,mm)*eqj(k,l)
            enddo
         enddo
         do mm = nclc+1,neqj+nclc
            l = neqj-nclc+mm
            gj(k,mm) = eqj(k,l)
            if (mm.gt.neqj.and.j1.gt.ndb) then
               gj(k,mm) = -eqj(k,l)
            endif
         enddo
         do l = 1,ncls
            gj(k,neqj+nclc+1) = gj(k,neqj+nclc+1)-umat(j,l,nclc+1)*
     &           eqj(k,l)
         enddo
      enddo


      return
      end




** ---------------------------------------------------------------------
*
      SUBROUTINE SURFOMEGA (ndt,neqj,ncls,pcol,dt,dtom,times,bj,xjmod,
     & depotbomega)
*
** ---------------------------------------------------------------------
* Calcul des coefficients necessaires dans la methode de henyey
*   lors du calcul de l'evolution du moment cinetique ;
*     CONDITION LIMITE DE SURFACE
*
* Version generalisee incluant au moins 2 conditions limites :
* 1) CL sur xpsi
* 2) CL sur l'eq. d'evolution du moment cinetique
* 3) CL sur lambda (si ZGRADMU = T)
*
* SURFOMEGA est appelee par HEN_OMEGA
*
* Arguments:
*     ndt = borne superieure de la zone de melange
*     neqj = nombre d'eq. du systeme a resoudre
*     ncls = nombre de CL de surface associees au syst. a resoudre
*     dt = dtom/rtot = pas de temps normalise au rayon
*                     (pour conditionnement de la matrice)
*
*     dtom = pas de temps de calcul
*     times = age du modele (en annees)
*     bj(neqj+ncls,2*neqj+1) = matrice pour le bloc limite surface
*     xjmod(nsh,neqj) = matrice regroupant les vecteurs associes aux variables
*             independantes du syst.
*
*     NOTA BENE: Par defaut, on aura xjmod(i,1) = om(i)
*                                    xjmod(i,2) = aux(i)
*                                    xjmod(i,3) = ur(i)
*                                    xjmod(i,4) = xpsi(i)
*                TOUTE NOUVELLE VARIABLE INDEPENDANTE SERA RAJOUTEE A LA SUITE
*                exemple: xjmod(i,5) = xlambda(i) si zgradmu = T
*
* Auteur: A.Palacios
*
C Modifs CC ondes (22/11/06)
*
* ATTENTION: Il y a une asymetrie entre les CL en surface et au centre.
*       En effet, pour Ur>0 (vitesse dirigee vers le centre de l'etoile
*            a l'equateur), le moment cinetique du coeur de l'etoile
*            va augmenter, alors que le moment de la zone convective
*            exterieure va lui diminuer.  D'ou le '-' dans l'equation
*            de b1, qui est un '+' dans l'equation de z1.
*
* ----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.rot'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.transp'
      include 'evolcom.var'
      include 'evolcom.eng'
      include 'evolcom.surf'

      integer ndt,j,j1,k,i
      integer neqj,neqj1,ncls
      integer pcol,blin,ndtenv
      integer idiffvr_sav
      integer mlp,mlp0
      integer ndbcore

      integer klcore,klenv,klpulse

      character omegaconv_sav*1
      logical init

      double precision xmoment_tot,xom
      double precision bj,xjmod
      double precision dt,om,ur,xpsi,xlambda,dtom,times
      double precision vomm,omm,dom
      double precision xmom_new_ce,xmom_old_ce
      double precision vit
      double precision vom,vrray,vxpsi
      double precision abmurj,abmuj
      double precision xlambdaold
      double precision dlam,xxnorm
      double precision Dhold,disctime
      double precision xnuvv_ndt
      double precision nturb
      double precision KK1,KK2,aexp,bexp,mexp,Bsat,Bsun,omsol,Mdotsun
     $     ,Mdotsat,B,mdot

      double precision dnuvvo_j,dnuvvp_j,dnuvvo_j1,dnuvvp_j1
      double precision tamp0,tamp1,tamp1a,tamp1b,tamp2,tamp3
      double precision dj,dj1
      double precision omdom
      double precision etapar,dmlinc,dms
      double precision Bcrit,fstar,Bequi
      double precision expm,expr,ff,fourm,cte1,cte2
      double precision mtini,totm,mdredgeup
      double precision omega_sat
      double precision depotbomega
      double precision chi,pexp,twom
      double precision tauc,taucSUN,taucnorm
      double precision faccrit,btorq,torq0
      double precision tc_hp,rc_hp,tg

      
      common /envel/ ndtenv
      common /dturbulente/ nturb(0:nsh)
      common /moment/ xmoment_tot
      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)
      common /lambdaold/ xlambdaold(nsh)
      common /gradmuj/ abmurj(nsh),abmuj(nsh)
      common /calcDh/ Dhold(nsh)
      common /disclocking/ disctime
      common /massloss/ etapar,dmlinc,dms,mlp,mlp0
      common /init_rot/ init,omegaconv_sav
      common /mass/ mtini,totm,mdredgeup
      common /boreas_Mloss/ Bcrit,fstar,Bequi
      common /overshoot/ klcore,klenv,klpulse

      common /tauc_hp/ tauc

      dimension bj(neqj+ncls,2*neqj+1),xjmod(nsh,neqj),om(nsh),
     &     ur(nsh),xpsi(nsh),xlambda(nsh)
      dimension pcol(neqj),blin(ncls)
      dimension omdom(nsh)

*-----------------------
***   Initialisations
*-----------------------
      


      neqj1 = 2*neqj+1
      xmom_new_ce = 0.d0
      xmom_old_ce = 0.d0
      tamp1 = 0.2d0*rray(ndt)**4*rdensm(ndt)*xLloc
      if (idiffvr.eq.4) then
         om(ndt:nmod) = 1.d0
      else
         om(nmod) = xjmod(nmod,pcol(1))
         om(ndt) = xjmod(ndt,pcol(1))
      endif
      om(ndt-1) = xjmod(ndt-1,pcol(1))
      ur(ndt) = xjmod(ndt,pcol(3))
      xpsi(ndt) = xjmod(ndt,pcol(4))
      if (zgradmu) xlambda(ndt) = xjmod(ndt,pcol(5))
      
      if (omegaconv.eq.'t') then

c$$$         ndbcore=2
c$$$         if (novlim(1,3).eq.1) ndbcore=novlim(1,4)
c$$$         do j = ndtenv,nmod
c$$$            j1 = j-1
c$$$            vomm = 0.5d0*(vom(j)+vom(j1))
c$$$            omm = 0.5d0*(om(j)+om(j1))
c$$$            xmom_new_ce = xmom_new_ce+rraym2(j)*dm(j1)*om(ndt)
c$$$            xmom_old_ce = xmom_old_ce+vrraym2(j)*dm(j1)*vomm
c$$$         enddo

c         xmom_new_ce = xmom_new_ce*om(ndt)

         vomm = 0.5d0*(vom(ndt)+vom(ndt-1))
         xmom_new_ce = xmom_new_ce+rraym2(ndt)*dm(ndt-1)*om(ndt)
         xmom_old_ce = xmom_old_ce+vrraym2(ndt)*dm(ndt-1)*vomm
      else
         do j = ndt+1,nmod
            j1 = j-1
            vomm = 0.5d0*(vom(j)+vom(j1))
            xmom_new_ce = xmom_new_ce+rraym2(j)*dm(j1)
            xmom_old_ce = xmom_old_ce+vrraym2(j)*dm(j1)*vomm
         enddo
         xmom_new_ce = xmom_new_ce*om(ndt)
      endif
      do j = 1,ncls+neqj
         do k = 1,neqj1
            bj(j,k) = 0.d0
         enddo
      enddo

      if (zgradmu) then
         blin(1) = 1
         blin(2) = 3
         blin(3) = 2
      else
         blin(1) = 1
         blin(2) = 2
      endif

c..   Boundary condition for Psi

      bj(blin(2),neqj1) = xpsi(ndt)
      bj(blin(2),pcol(4)) = 1.d0


c..   no coupling and solid rotation in the envelope
      if (idiffvr.eq.0) then

         if (omegaconv.eq.'t')  then
            bj(blin(1),neqj1) = ur(ndt)*ur_S*om(ndt)
            bj(blin(1),pcol(1)) = ur(ndt)*ur_S
            bj(blin(1),pcol(3)) = om(ndt)*ur_S
         else
            bj(blin(1),neqj1) = om(ndt)*ur(ndt)*ur_S*dt+(xmom_new_ce
     &           -xmom_old_ce)/tamp1-depotbomega/tamp1
            bj(blin(1),pcol(1)) = ur(ndt)*ur_S*dt+xmom_new_ce/(tamp1
     &           *om(ndt))
            bj(blin(1),pcol(3)) = om(ndt)*ur_S*dt
         endif
      else
 
c..   coupling and solid rotation in the envelope         
         if (idiffvr.ne.6) then

            if (idiffvr.eq.5.or.idiffvr.ge.8) then
***   Kawaler law for magneitc braking
*     Couple proportionnel a omega_S**3
*     + saturation a omega=.15 10^-4 (environ 10 km/s pour le soleil)
c     if (om(nmod)*omega_S.lt.1.5d-5) then
               if (idiffvr.eq.5) then
c                  omega_sat = 10.d0*omsun !2.8d-5 ! 10 omegasun
                  omega_sat = om_sat*omsun
                  if (om(nmod)*omega_S.lt.omega_sat) then
                     dmoment = -diffcst*dtom*(om(ndt)*omega_S)**3
     &                    *sqrt(solr/totm)
***begin
!Reiners & Mohanty 2012 1/3
cc                     dmoment = -diffcst*dtn/(omega_sat**4)*(om(ndt)
cc     &                *omega_S)**5*(solr)**(1.6d1/3.d0)*
cc     &                totm**(-pw23)
***end
                     xmoment_tot = xmom_tots+dmoment
                  else
                     dmoment = -diffcst*dtom*om(nmod)*omega_S* 
     &                    omega_sat**2*sqrt(solr/totm)
***begin
!Reiners & Mohanty 2012 2/3
cc                     dmoment = -diffcst*dtn*(om(ndt)*
cc     &                    omega_S)*(solr)**(1.6d1/3.d0)*totm**(-pw23)
***end

                     xmoment_tot = xmom_tots+dmoment
                  endif
***   Matt et al. 2012  equation 9 using outputs of the boreas routine 
***   from Cranmer & Saar 2011
               else if (idiffvr.eq.8.and.(mlp.eq.23.or.mlp.eq.24)) then
                  KK1 = 1.d0*(m(nmod)/msun)**2-4.3d0*m(nmod)/msun
     &                 +10.0d0

c.. According to G&B2013 
cc                  mexp = 0.22
c.. GB13 for fast rotator
cc                  KK1 = 1.7
c.. GB13 for slower rotator
cc                  KK1 = 1.8
                  if (KK1.lt.0.5d0) KK1 = 0.5d0 
                  KK2 = 0.0506d0

                  mexp = 0.17d0
                  fourm = 4.d0*mexp
                  expm = 1.d0-2.d0*mexp
                  expr = 5.d0*mexp+2.d0
                  aexp = 1.3d0
                  bexp = 1.65d0
c                  ff = 0.1d0
                  ff = (om(ndt)*omega_S)**2.d0*
     &                 (rray(nmod)*rsun)**3.d0/(2.d0*g*(m(nmod)))
                  ff = dsqrt(ff)
                  cte1 =  KK1**2/(2.d0*g)**mexp
                  cte2 = 1.d0/(KK2**2+0.5d0*ff**2)**mexp
                  mdot = dms*msun/sec
c                  B = Bcrit!*fstar
                  B = 1.13d0*Bequi*fstar
                  dmoment = -dtom*cte1*B**fourm*mdot**expm*
     $                 (rray(nmod)*rtot)**expr/(m(nmod))**mexp*
     $                 om(ndt)*omega_S*cte2


***   Matt et al. 2015
               elseif (idiffvr.eq.9) then
                  chi = om_sat
                  mexp = 0.22d0
c                  mexp = 0.24d0 ! test TD 10/2019
                  pexp = 2.1d0
c                  pexp = 1.5d0 ! test TD 11/2019
                  twom = -2.d0*mexp
c$$$                  tauc = 3.1424d2*exp(-teff/1952.5d0)*exp(-(teff/
c$$$     $                 6250.d0)**18)+2.d-3

                  tg = 0.d0
                  do k = novlim(klenv,3),novlim(klenv,4)
                    if (sconv(k).gt.0.d0) tg = tg+(r(k+1)-r(k))/sconv(k)
                  enddo

                  k = novlim(klenv,3)
                  do while (r(k).lt.r(novlim(klenv,3))+
     &                 0.5d0*hp(novlim(klenv,3)).and.k.lt.nmod)
                     k = k+1
                  enddo

                  if (k.ge.novlim(klenv,4).or.sconv(k).eq.0.d0) then
                     tauc = tg
                  else
                     tc_hp = alphac*hp(k)/sconv(k)
                     rc_hp = r(k)
                  endif
                     
                  tauc = tc_hp/8.64d4

!     print *,'tauc=',tauc
                  
                  taucSUN = 12.90430912750d0
                  taucnorm = tauc/taucSUN
                  faccrit = om(ndt)*omega_S*dsqrt(r(nmod)/
     &                 gmr(nmod)) !! Need to include the deformation by rotation
                  omega_sat = chi*omsun/taucnorm
                  btorq = dsqrt(1.d0+(faccrit/0.072)**2.d0)
                  torq0 = diffcst*(rtot/rsun)**3.1d0*(m(nmod)/msun)**0.5
     &                 *btorq**twom

!!                  if (om(nmod)*omega_S.lt.omega_sat) then
                  if (om(ndt)*omega_S.lt.omega_sat) then                     
                     dmoment = -dtom*torq0*taucnorm**pexp
     &                    *(abs(om(ndt))*omega_S/omsun)**(pexp+1.d0)
                  else
                     dmoment = -dtom*torq0*(chi**pexp)*(abs(om(ndt))
     &                   *omega_S/omsun)
                  endif

c$$$                  print *,'torq0,dtom,taucnorm',torq0,dtom,taucnorm,
c$$$     &                 -dtom*torq0*taucnorm**pexp*(abs(om(ndt))*
c$$$     &                 omega_S/omsun)**(pexp+1.d0)

               endif
               dmoment1 = dmoment
               if (omegaconv.eq.'t') then
                  dmoment = dmoment/xLtot
                  bj(blin(1),neqj1) = -dmoment+(xmom_new_ce
     $                 -xmom_old_ce)/pim4
                  bj(blin(1),pcol(1)) = (xmom_new_ce/pim4-3.d0*dmoment)
     &                 /om(ndt)
                  xmoment_tot = xmom_tots+dmoment
               else
                  dmoment = dmoment/(tamp1*xLtot)                  
                  bj(blin(1),neqj1) = om(ndt)*ur_S*ur(ndt)*dt-dmoment
     &                 +(xmom_new_ce-xmom_old_ce)/tamp1
     &                 -depotbomega/(tamp1)
c     Uniquement pour Kawaler, à adapter pour Matt 2012 (?)
                  
                  if ((om(nmod)*omega_S.lt.omega_sat).and.(idiffvr
     &                 .eq.5)) then
                     bj(blin(1),pcol(1)) = ur(ndt)*ur_S*dt-3.d0*dmoment
     &                    /om(ndt)
***
!Reiners & Mohanty 2012 3/3
cc                     bj(blin(1),pcol(1)) = ur(ndt)*ur_S*dt-5.d0*dmoment
cc     &                    /om(ndt)
***
     &                    +xmom_new_ce/(om(ndt)*tamp1)
                  elseif ((om(nmod)*omega_S.lt.omega_sat).and.(idiffvr
     &                    .eq.9)) then
                     bj(blin(1),pcol(1)) = ur(ndt)*ur_S*dt-(pexp+1.d0+
     &                    2.0*mexp)*dmoment/om(ndt)
     &                    +xmom_new_ce/(om(ndt)*tamp1)
                  else 
                     bj(blin(1),pcol(1)) = ur(ndt)*ur_S*dt-dmoment
     &                    /om(ndt)
     &                    +xmom_new_ce/(om(ndt)*tamp1)
                  endif
                  bj(blin(1),pcol(3)) = om(ndt)*ur_S*dt
                  xmoment_tot = xmom_tots+dmoment

               endif
               

               

c$$$            else if(idiffvr.eq.9) then
c$$$c... Magnetic field Nadège (Mai 2011)
c$$$	    
c$$$c	    	print *, 'idiffvr 9'
c$$$               stop "Check braking"
c$$$C... calcul vitesse d'échappement
c$$$c		vescap=(2*g*totm/rtot)**0.5 
c$$$		
c$$$c		write(*,*) vescap
c$$$C... Calcul de la vitesse des vents à l'infini 
c$$$c		if (teff.geq.21000)then 
c$$$c			vwinf=2.65*vescap
c$$$c		elseif(teff.leq.10000)then 
c$$$c		        vwinf=vescap
c$$$c		else 
c$$$c		        vwinf=1.4*vescap
c$$$c		endif
c$$$		
c$$$C... formule de Meynet,Eggenberger & Maeder 2010  
c$$$		
c$$$c		etamagn = diffcst**2*rtot**2/(dms*vwinf)
c$$$c		dmoment = dtom*(0.29+(etamagn+0.25)**(0.25))**2
c$$$c		dmoment =(2/3)*dms*om(nmod)*omega_S*(rtot)**2*dmoment
c$$$		
c$$$c		dmoment = dmoment/tamp1
c$$$		
c$$$c		bj(blin(1),neqj1)=om(ndt)*ur_S*ur(ndt)*dt-dmoment
c$$$c     &        +(xmom_new_ce-xmom_old_ce)/tamp1
c$$$     
c$$$c     		if (omegaconv.eq.'t') then
c$$$c		bj(blin(1),pcol(1)) = ur(ndt)*ur_S*dt+(xmom_new_ce/
c$$$c     &                 tamp1-dmoment)/om(ndt)
c$$$c               else
c$$$c                  bj(blin(1),pcol(1)) = ur(ndt)*ur_S*dt-dmoment/om(ndt)
c$$$		
c$$$c     	       endif 
c$$$c     		bj(blin(1),pcol(3)) = om(ndt)*ur_S*dt
c$$$     		
c$$$c.... fin magnetic field 


             
            else
*     Calcul de la vitesse de rotation de la surface selon la loi choisie
               call calc_omega (times,xom,vit)
               if (idiffvr.eq.4) then
                  omega_S = xom
                  xjmod(ndt,pcol(1)) = 1.d0
                  if (omegaconv.eq.'t')  then
                     bj(blin(1),neqj1) = ur(ndt)*ur_S*xom
                     bj(blin(1),pcol(1)) = ur(ndt)*ur_S
                     bj(blin(1),pcol(3)) = xom*ur_S
                  else
                     bj(blin(1),neqj1) = xom*ur(ndt)*ur_S*dt
     $                    +(xmom_new_ce-xmom_old_ce)/tamp1
                     bj(blin(1),pcol(1)) = ur(ndt)*ur_S*dt+xmom_new_ce
     $                    /(tamp1*xom)
                     bj(blin(1),pcol(3)) = xom*ur_S*dt
                  endif
               else
                  xom = xom/omega_S
*     Couple < 0; dmoment < 0
*     couple applique par unite de temps
*     dmoment=couple*dt
                  dmoment = xmom_new_ce*xom-xmom_old_ce+xom*ur(ndt)*ur_S
     &                *dt*tamp1
                  xmoment_tot = xmom_tots+dmoment*xLtot
                  bj(blin(1),neqj1) = om(ndt)-xom
                  bj(blin(1),pcol(1)) = 1.d0
               endif                   
            endif
            
         else
c            if (omegaconv.eq.'f') then
            xmom_new_ce = 0.d0
            xmom_old_ce = 0.d0
            xnuvv_ndt = xnumm(ndt)+nturb(ndt+1)*nturb(ndt)*xnuvv(ndt)
            
c..   constant specific env. AM (no torque applied)            
            if (idiffvr.eq.6) then
***   CL sur xpsi
               bj(blin(2),neqj1) = -pw43*om(ndt)**2*rray(ndt)/(grav(ndt)
     $              *deltaKS(ndt)) + xpsi(ndt)
               bj(blin(2),pcol(1)) = -pw43*2.d0*om(ndt)*rray(ndt)/
     &              (grav(ndt)*deltaKS(ndt))
               bj(blin(2),pcol(4)) = 1.d0
               
***   CL sur l'equation d'evolution du moment cinetique
               
               do j = ndt+1,nmod
                  j1 = j-1
                  vomm = 0.5d0*(vom(j)+vom(j1))
                  xmom_old_ce = xmom_old_ce+vrraym2(j)*vomm*dm(j1)
                  xmom_new_ce = xmom_new_ce+(2.d0+rray(j1)/rray(j)*(1.d0
     &                 +rray(j1)/rray(j))+rray(j)/rray(j1)*(1.d0+rray(j)
     &                 /rray(j1)))*dm(j1)
               enddo
               
               bj(blin(1),neqj1) = (0.2d0*rdensm(ndt)*rray(ndt)**4*
     &              ur(ndt)*om(ndt)*ur_S-2.d0*om(ndt)*rray(ndt)**3*
     &              rdensm(ndt)*xnuvv_ndt/rtot)*dt+(pw16*xmom_new_ce*
     &              om(ndt)*rray(ndt)**2-xmom_old_ce)/xLloc 
               bj(blin(1),pcol(1)) = (0.2d0*rdensm(ndt)*rray(ndt)**4*
     &              ur_S*ur(ndt)-2.d0*rray(ndt)**3*rdensm(ndt)*xnuvv_ndt
     &              /rtot)*dt+pw16*xmom_new_ce*rray(ndt)**2/xLloc
               
               bj(blin(1),pcol(3)) = 0.2d0*rdensm(ndt)*rray(ndt)**4*ur_S
     &              *om(ndt)*dt
               
c..   rotation profile derived from 3D hydro simulations (Brun & Palacios 2009)  
               
            elseif (idiffvr.eq.7) then
               
c               call ash_rgb_prof(ndt,ndt,om,omdom,xpsi)
               call ash_rgb_prof(ndt,om,omdom,xpsi)
               
***   Boundary Condition for xpsi
               bj(blin(2),neqj1) = -pw23*omdom(ndt)*rray(ndt)**2/
     $              (grav(ndt)*deltaKS(ndt)) + xpsi(ndt)
               bj(blin(2),pcol(4)) = 1.d0

***   Boundary Condition for angular momentum equation
!               do j = ndt+1,nmod
!                  j1 = j-1
!                  call ash_rgb_prof(ndt,j,om,omdom,xpsi)
!                  vomm = 0.5d0*(vom(j)+vom(j1))
!                  omm = 0.5d0*(om(j)+om(j1))
!                  xmom_old_ce = xmom_old_ce+vrraym2(j)*vomm*dm(j1)
!                  xmom_new_ce = xmom_new_ce+rraym2(j)*omm*dm(j1)
!               enddo
               bj(blin(1),neqj1) = (0.2d0*rdensm(ndt)*rray(ndt)**4*
     $              ur(ndt)*om(ndt)*ur_S-omdom(ndt)*rray(ndt)**4
     $              *rdensm(ndt)*xnuvv_ndt/(rtot*om(ndt)))
     $              *dt +(xmom_new_ce-xmom_old_ce)/xLloc 

               bj(blin(1),pcol(1)) = (0.2d0*rdensm(ndt)*rray(ndt)**4*
     $              ur_S*ur(ndt)+omdom(ndt)*rray(ndt)**4
     $              *rdensm(ndt)*xnuvv_ndt/(rtot*om(ndt)**2))*dt
               
               bj(blin(1),pcol(3)) = 0.2d0*rdensm(ndt)*rray(ndt)**4*ur_S
     &              *om(ndt)*dt
 
            endif
         endif
      endif
*----------------------------------------------------------------
***   Prise en compte des variations de poids moleculaire moyen
*----------------------------------------------------------------

      if (zgradmu) then
         xxnorm = xnorm*gs
         dlam = xlambda(ndt) - xlambdaold(ndt)
         tamp3 = 6.d0*Dhold(ndt)/(rtot*ur_S*rray(ndt)**2)

         bj(blin(3),neqj1) =  dlam*xxnorm/ur_S + dt*(ur(ndt)*
     &        rtot*abmurj(ndt)+tamp3*xlambda(ndt)*xxnorm)
         bj(blin(3),pcol(5)) =  xxnorm/ur_S + tamp3*dt*xxnorm
         bj(blin(3),pcol(3)) = dt*rtot*abmurj(ndt)
      endif


*-----------------------------------------------------------------------
***   Conditionnement pour remplissage du bloc matriciel dans HEN_OMEGA
*-----------------------------------------------------------------------

      do i = 1,ncls
         bj(blin(i),neqj1) = - bj(blin(i),neqj1)
      enddo
      do i = ncls+1,ncls+neqj
         bj(i,neqj1) = - bj(i,neqj1)
      enddo
     
      if (idiffvr.ge.8) omegaconv = omegaconv_sav 
c      idiffvr = idiffvr_sav

      return
      end


*-----------------------------------------------------------------------------
*
      SUBROUTINE HEN_OMEGA (ndb,ndt,neqj,ncls,nclc,pcol,dtom,alph1,
     &     times,xjmod,error)
*
*-----------------------------------------------------------------------------
* Resolution des equations aux derivees partielles du premier ordre
* decrivant l'evolution du moment cinetique par la methode de relaxation
* de Newton-Raphson pour un schema matriciel de type Henyey
*
* Resolution de la matrice de la surface vers le centre (indices decroissants)
*
* Version generalisee
*
* HEN_OMEGA est appelee par DIF_OMEGA
* et appelle SURFOMEGA, INTOMEGA, CENTOMEGA et NWRMAT
*
*
* Arguments:
*
*     ndt = borne superieure de la zone de melange
*     neqj = nombre d'eq. du syst. a resoudre pour calculer l'evolution
*            du moment cinetique
*     ncls = nombre de CL en surface associees au syst.
*     nclc = nombre de CL au centre associees au syst.
*     dtom = pas de temps de calcul
*     alph1 = parametre de convergence pour la methode de Newton-Raphson
*     mr(nsh) = masse relative
*     times = age du modele (en annees)
*     xjmod(nsh,neqj) = matrice regroupant les vecteurs associes aux variables
*             independantes du syst.
*
*
*     L'ordonnancement de la matrice "pseudo-jacobienne" est fixe par un
*     pointeur, PCOL, qui sera a initialiser pour chaque configuration.
*     Ce systeme permet en particulier d'arranger la matrice de facon a ce 
*     qu'elle ne soit pas singuliere.
*     ACTUELLEMENT:
*        pcol(1) associe a omega
*        pcol(2) associe a aux
*        pcol(3) associe a ur
*        pcol(4) associe a xpsi
*        pcol(5) associe a lambda
*
*    Toute nouvelle variable sera rangee dans pcol de facon a apparaitre a la
*    suite de ces 5 variables deja definies
*      # Cas du systeme a 4 eq.: xjmod(i,pcol(1)) = xjmod(i,1) = om(j)
*                                xjmod(i,pcol(2)) = xjmod(i,2) = aux(i)
*                                xjmod(i,pcol(3)) = xjmod(i,3) = ur(i)
*                                xjmod(i,pcol(4)) = xjmod(i,4) = xpsi(i)
*
*      # Cas du systeme a 5 eq.: xjmod(i,pcol(1)) = xjmod(i,1) = om(j)
*                                xjmod(i,pcol(2)) = xjmod(i,2) = aux(i)
*                                xjmod(i,pcol(3)) = xjmod(i,4) = ur(i)
*                                xjmod(i,pcol(4)) = xjmod(i,5) = xpsi(i)
*                                xjmod(i,pcol(5)) = xjmod(i,3) = xlambda(i)
*
*
* Auteur : A.Palacios (15/03/2000)
*
*
c..   Avril 2004: rajout du traitement moment specifique constant
c..   ----------  dans les zones convectives (cas idiffvr = 6)

*-----------------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.mod'
      include 'evolcom.teq'
      include 'evolcom.transp'
      include 'evolcom.igw'
c      include 'evolcom.transpondesexcit'

      logical implicit_wave

      integer i,j,k,j1,l,ll,ii
      integer itmax,iterh,iiterondes
      integer ndb,ndt
      integer ndim,ndim2
      integer jgd,jgg
      integer neqj,neqj1,nclc,ncls,error
      integer pcol
      integer ndtenv
      integer klcore,klenv,klpulse

      double precision om,aux,ur,xpsi,xlambda,dtom,alph1,times,nturb
      double precision dift,dife
      double precision eqj,bj,gj,zj
      double precision umat,bmat,zmat,hmat
      double precision dt,fagd,agdj,dxj,xjmod
      double precision agg,agmax,xmultiplicateur
      double precision vgdj
      double precision gg
      double precision gdj
      double precision gdtamp
      double precision fredj,fred
      double precision epsa
      double precision Dhold,disctime
      double precision omprint(nmod)

c      double precision depotwaves
      logical init
      character omegaconv_sav*1

      common /envel/ ndtenv
      common /difcirc/ dift(nsh),dife(nsh)
      common /dturbulente/ nturb(0:nsh)
      common /tolang/ epsa
      common /calcDh/ Dhold(nsh)
      common /disclocking/ disctime
      common /init_rot/init,omegaconv_sav
      common /overshoot/ klcore,klenv,klpulse

      dimension jgd(neqj),jgg(neqj)
      dimension om(nsh),aux(nsh),ur(nsh),xpsi(nsh),xlambda(nsh)
     &     ,gj(neqj,neqj+nclc+1),bmat(neqj+ncls,ncls+1),
     &     umat(nsh,neqj,nclc+1),zmat(neqj+nclc),hmat(neqj,ncls+1)
      dimension bj(neqj+ncls,2*neqj+1),zj(neqj+nclc,neqj+nclc+1)
      dimension xjmod(nsh,neqj),eqj(neqj,2*neqj+1)
      dimension agdj(neqj),fagd(neqj),dxj(nsh,2*neqj)
      dimension vgdj(neqj),gg(neqj),gdj(neqj),gdtamp(neqj),fred(neqj)
      dimension pcol(neqj)

c      dimension depotwaves(nsh)

*--------------------
***   Initialisations
*--------------------

*     Normalisation du pas de temps pour l'evolution du moment cinetique
*     afin de conserver les dimensions des equations.

      dt = dtom/rtot
      neqj1 = 2*neqj+1

*     Tolerances pour la convergence de la methode

      do j = 1,neqj
         fagd(j) = 1.d0
         agdj(j) = epsa*fagd(j)
         vgdj(j) = 1.d0
         jgg(j) = 1
         jgd(j) = 1
      enddo

c      agg = epsa*1.d-4
c      agg = epsa*1.d-1
c      agg = 1.d-4
C Modif CC ondes (5/4/7) -->
C      agg = 1.d-6
      agg = 1.d-4
C      agmax = 1.d25
      agmax = 1.d40
C <--
      fredj = 1.d0
      xmultiplicateur = 1.4d0
      itmax = 900

      do i = 1,nmod
        Dhd(i) = 0.d0
        dift(i) = 0.d0
        om(i) = xjmod(i,pcol(1))
        aux(i) = xjmod(i,pcol(2))
        ur(i) = xjmod(i,pcol(3))
        xpsi(i) = xjmod(i,pcol(4))
        if (zgradmu) then
           xlambda(i) = xjmod(i,pcol(5))
        else
           xlambda(i) = 0.d0
        endif
      enddo
c      if (ndb.gt.1) xpsi(ndb) = xpsi(ndb-1)

c..   On s'assure que la CL xpsi = 0 aux bornes de la zone de melange
c..   soit satisfaite
      if (idiffvr.lt.6.and.omegaconv.ne.'t') then
         if (ndt.lt.nmod) then
            xpsi(ndt) = xpsi(ndt+1)
            xjmod(ndt,pcol(4)) = xjmod(ndt+1,pcol(4))
         endif
      endif

*-----------------------
***   Processus iteratif
*-----------------------
      do i = 1,nmod
         do j = 1,2*neqj
            dxj(i,j) = 0.d0
         enddo
      enddo

      iterh = 0
 77   iterh = iterh+1
c      alph1 = xmultiplicateur*alph1
c      alph1 = min(alph1,1.d0,fredj)

*     Calcul du coefficient vertical de diffusion turbulente sur base
*     des profils a l'iteration precedente
      do i = 1,nmod
         xpsi(i) = xjmod(i,pcol(4))
         om(i) = xjmod(i,pcol(1))
         ur(i) = xjmod(i,pcol(3))
      enddo
      call turbu_dhold (ndb,ndt,xpsi,xlambda,om,ur)

C -----------------------
c Update tranport by IGW
C -----------------------


c      implicit_wave = .true.
       implicit_wave = .false.
      iiterondes = iterh
      if (igwrot) then
         if (implicit_wave.or.iiterondes.eq.1) 
     &      call transport_ondes (om,ndb,ndtenv,dtom,
     &        iiterondes)
      else
          do i = 1,nmod
             depotwaves(i) = 0.d0
          enddo
      endif


*--------------------------------------------------------
***   Remplissage et inversion du bloc limite surface B
***   contient neqj+ncls lignes (= equations)
***   cf. eq. (5.46) these A. Palacios
*--------------------------------------------------------

      j = ndt
      j1 = ndt-1

      do k = ndb,ndt
         do i = 1,neqj
            do ii = 1,nclc+1
               umat(k,i,ii) = 0.d0
            enddo
         enddo
      enddo
      do i = 1,neqj
         do ii = 1,neqj+nclc+1
            gj(i,ii) = 0.d0
         enddo
         gg(i) = 0.d0
         gdj(i) = 0.d0
      enddo



      call surfomega (j,neqj,ncls,pcol,dt,dtom,times,bj,xjmod,
     & depotwaves(ndt+1))
      
      call intomega (j,j1,iterh,ndb,ndt,neqj,ncls,nclc,pcol,dt,umat,
     &     xjmod,eqj,gj,depotwaves(j))
      
      do i = 1,neqj
         gg(i) = eqj(i,neqj1)
      enddo

      do i = ncls+1,ncls+neqj
         l = i-ncls
         do k = 1,neqj+nclc+1
            bj(i,k+ncls) = gj(l,k)
         enddo
         do k = 1,ncls
            bj(i,k) = eqj(i-ncls,k)
         enddo
      enddo

      ndim = ncls+neqj
      ndim2 = nclc+1

      call nwrmat (bj,bmat,ndim,ndim2,error)


      if (error.gt.0) then    
         error = 18
         print *,"probleme en surface Bj at iteraction ",iterh
         print *,ndim,ndim2,neqj,ncls,ncls+neqj,2*neqj+1
         do i = 1,neqj+ncls
            write (nout,667) i,(bj(i,k),k=1,2*neqj+1)
         enddo
         print *,"Gj en j =",j
         print *,'derive en j :',1,neqj
         print *,'derive en j1:',neqj,2*neqj
         print *,'equation    :',2*neqj+1
         print *,ndim,ndim2,neqj,ncls,ncls+neqj,2*neqj+1
         do i = 1,neqj
            write (nout,667) i,(eqj(i,k),k=1,2*neqj+1)
         enddo
         print *,ndim,ndim2,neqj,ncls,ncls+neqj,2*neqj+1
         do i = 1,neqj+ncls
            write (nout,667)i,(bj(i,k),k=1,2*neqj+1)
 667        format (i2,1x,11(1pe11.4,1x))
         enddo
         return
      endif

c      if (idiffvr.eq.4.and.time/sec.lt.disctime) xjmod(ndt:nmod,pcol(1))
c     $     = breaktime

*     Initialisation de la matrice umat, qui sera utilisee
*     pour construire les bloc matriciels ventraux "pseudo-jacobiens"
*     a inverser

      do i = 1,neqj+ncls
         if (i.le.neqj) then
            do k = 1,nclc+1
               umat(j,i,k) = bmat(i,k)
            enddo
         else
            do k = 1,nclc+1
               umat(j1,i-neqj,k) = bmat(i,k)
            enddo
         endif
      enddo

*---------------------------------------------------------------------
***   Remplissage et inversion des blocs matriciels ventraux
***   contiennent neqj lignes (= equations)
***   Matrices "pseudo-jacobienne" (cf. eq. (5.53) these A. Palacios)
***   Couches ndt-1 a 3
*---------------------------------------------------------------------

 101  j = j-1
      j1 = j-1

      call intomega (j,j1,iterh,ndb,ndt,neqj,ncls,nclc,pcol,dt,umat,
     &     xjmod,eqj,gj,depotwaves(j))
      
      do i = 1,neqj
         if (abs(gg(i)).lt.abs(eqj(i,neqj1))) then
            gg(i) = eqj(i,neqj1)
            jgg(i) = j
         endif
      enddo

      if (j1.eq.ndb) goto 102

      ndim = neqj
      ndim2 = nclc+1

*     Calcul de (Ui,Vi,Wi) par inversion

      call nwrmat (gj,hmat,ndim,ndim2,error)
      if (error.gt.0) then
         error = 18
         print *,'derive en j :',1,neqj
         print *,'derive en j1:',neqj,2*neqj
         print *,'equation    :',2*neqj+1
         print *,ndim,ndim2,neqj,ncls,ncls+neqj,2*neqj+1
         do i = 1,neqj
            write (nout,667) i,(eqj(i,k),k=1,2*neqj+1)
         enddo
         print *,'stop in nwrmat, shell #',j
         return
      endif

      do i = 1,neqj
         if (i.le.nclc) then
            do k = 1,nclc+1
               umat(j,i+ncls,k) = hmat(i,k)
            enddo
         else
            do k = 1,nclc+1
               umat(j1,i-nclc,k) = hmat(i,k)
            enddo
         endif
      enddo

      goto 101


*--------------------------------------------------------
***   Remplissage et inversion du bloc limite central Z
***   Contient neqj+nclc lignes (= equations)
***   cf. eq. (5.48) these A. Palacios)
*--------------------------------------------------------

 102  do i = 1,neqj
         do k = 1,neqj+nclc+1
            zj(i,k) = gj(i,k)
         enddo
      enddo

      call centomega (ndb,neqj,nclc,pcol,dt,zj,xjmod)

      ndim = neqj+nclc
      ndim2 = 1

      call nwrmat (zj,zmat,ndim,ndim2,error)
      if (error.gt.0) then
         error = 18
         print *,neqj,nclc,neqj+nclc
         do i = 1,neqj+nclc
            write (nout,667)i,(zj(i,k),k=1,2*neqj+1)
         enddo
         print *,'pb in nwrmat - central bloc, shell #',j
         return
      endif


*---------------------------------------
***   Corrections (dxj) au point central
*---------------------------------------

      do i = 1,neqj+nclc
         if (i.le.nclc) dxj(j,i+ncls) = zmat(i)
         if (i.gt.nclc) dxj(j1,i-nclc) = zmat(i)
      enddo


*--------------------------------------------------------------
***   Determination des variations relatives obtenues au centre
*--------------------------------------------------------------
      do i = 1,neqj
         
         if (abs(xjmod(ndb,i)).gt.1.d-25) then
            gdj(i) = dxj(j1,i)/xjmod(ndb,i)
         else
            gdj(i) = 0.d0
         endif
         do l = neqj+1,neqj+nclc
*     En cas de condition limite nulle pour la ieme variable
            if (abs(xjmod(ndb,i)).eq.abs(zj(l,neqj+nclc+1))) then
               if (abs(xjmod(ndb+1,i)).gt.1.d-25) then
                  gdj(i) = dxj(j,i)/xjmod(ndb+1,i)
               else 
                  gdj(i) = 0.d0
               endif
            endif
         enddo
      enddo

*------------------------------------------------------------------
***   Nouvelles valeurs pour les variables independantes au point 1
*------------------------------------------------------------------

      do i = 1,neqj
         xjmod(ndb,i) = xjmod(ndb,i)+alph1*dxj(j1,i)
      enddo

      do i = 1,ncls
         dxj(j,i) = 0.d0
         do k = 1,nclc
            dxj(j,i) = dxj(j,i)+umat(j,i,k)*dxj(j,k+ncls)
         enddo
         dxj(j,i) = dxj(j,i)+umat(j,i,nclc+1)
      enddo

*-----------------------------
****  Remontee des corrections
*-----------------------------

      do j = ndb+1,ndt

         if (j.gt.ndb+1) then
            do i = 1,neqj
               dxj(j,i) = dxj(j-1,i)
            enddo
         endif

*     Recherche de l'erreur relative maximale
         do i = 1,neqj
            if (abs(xjmod(j,i)).gt.1.d-15) then
               gdtamp(i) = dxj(j,i)/xjmod(j,i)
               if (abs(gdtamp(i)).gt.abs(gdj(i))) then
                  gdj(i) = gdtamp(i)
                  jgd(i) = j
               endif
            endif
         enddo

*     Recherche de l'erreur maximale
c         do i = 1,neqj
c            gdtamp(i) = dxj(j,i)
c            if (abs(gdtamp(i)).gt.abs(gdj(i))) then
c               gdj(i) = gdtamp(i)
c            endif
c         enddo

*     Correction des valeurs

         do i = 1,neqj
            xjmod(j,i) = xjmod(j,i)+alph1*dxj(j,i)
         enddo
         if (crz(j).le.1.or.j.ge.novlim(klenv,3)) then
c         if (crz(j).lt.0) then
            xjmod(j,pcol(3))= 0.d0
            if (zgradmu) xjmod(j,pcol(5))= 0.d0
         endif

*     Corrections pour le point suivant

         j1 = j+1

        if (j1.le.ndt) then
           do i = ncls+1,neqj
              dxj(j1,i) = dxj(j,i)
           enddo
           do i = ncls+1,neqj
              dxj(j,i) = 0.d0
              do k = 1,nclc
                 dxj(j,i) = dxj(j,i)+umat(j1,i,k)*dxj(j1,k+ncls)
              enddo
              dxj(j,i) = dxj(j,i)+umat(j1,i,nclc+1)
           enddo
        endif

        if (j1.le.ndt) then
           do ll = 1,ncls
              dxj(j,ll) = 0.d0
              if (j1.lt.ndt) then
                 do k = 1,nclc
                    dxj(j,ll) = dxj(j,ll)+umat(j1,ll,k)*dxj(j,k+ncls)
                 enddo
              else if (j1.eq.ndt) then
                 do k = 1,nclc
                    dxj(j,ll) = dxj(j,ll)+umat(j1,ll,k)*dxj(j1,k+ncls)
                 enddo
              endif
              dxj(j,ll) = dxj(j,ll)+umat(j1,ll,nclc+1)
           enddo
        endif
      enddo

*     Ajustement de alph en fonction de l'alternance ou non des
*     corrections

      if (iterh.eq.1) then
         do i = 1,neqj
            vgdj(i) = 0.d0
         enddo
      endif

      fredj = 1.d0
      do i = 1,neqj
         fred(i) = 0.5d0*(1.d0+abs(gdj(i)+vgdj(i))/(abs(gdj(i))+
     &        abs(vgdj(i))+1.d-30))
         fredj = min(fredj,fred(i))
      enddo
      do i = 1,neqj
         vgdj(i) = gdj(i)
      enddo


*-----------------------------
***   Test des corrections
*-----------------------------

*     Si une des corrections est superieure a la correction maximale
*     toleree, il faut effectuer une nouvelle iteration.

       do i = 1,neqj
c          if (abs(gdj(i)).ge.agdj(i).or.abs(gg(i)).ge.agg) goto 114
          if (abs(gdj(i)).ge.agdj(i).or.abs(gg(i)).ge.1.d50) goto 114
c          if (abs(gdj(i)).ge.agdj(i)) goto 114
       enddo
       if (neqj.eq.5) then
          write (nout,2007) iterh,(jgd(i),gdj(i),jgg(i),gg(i),i= 1,neqj)
       else if (neqj.eq.4) then
          write (nout,2008) iterh,(jgd(i),gdj(i),jgg(i),gg(i),i= 1,neqj)
       endif

c      xjmod(1,pcol(1)) = xjmod(2,pcol(1))
      if (xjmod(1,pcol(1)).lt.0.d0) then
         print *,'omega < 0 at shell 1'
         error = 18
         return
      endif        
      do i = max(2,ndb),ndt
         if (xjmod(i,pcol(1)).lt.0.d0) then
            j = min(i,nmod1)
            j = max(2,j)
            print *,'omega < 0 at shell',i,xjmod(j,pcol(1)),
     &           xjmod(j-1,pcol(1)),xjmod(j+1,pcol(1))
            error = 18
            return
         endif
c..   omega must be within a factor of 2 between consecutive shells
         if (nphase.gt.3.and.
     &       abs(log(xjmod(i,pcol(1))/(xjmod(i-1,pcol(1))+1.0d-30))).gt.
     &        1.d0.and.j.lt.nmod1) then
c     &        0.3d0) then
            j = min(i,nmod1)
            j = max(2,j)
            print *,'Too large variation of omega at shell',i,
     &           xjmod(j,pcol(1)),
     &           xjmod(j-1,pcol(1)),xjmod(j+1,pcol(1))
            error = 18
            return
         endif
      enddo

      write (nout,'(5x,a23)') 'convergence successful'

      return

*     On n'a pas obtenu la convergence

 114  if (iterh.lt.itmax) then
         do i = 1,neqj
cl            if (abs(gdj(i)).ge.agmax.or.abs(gg(i)).ge.agmax) goto 115
            if (abs(gdj(i)).ge.agmax.or.abs(gg(i)).gt.1.d50) goto 115
c            if (abs(gdj(i)).ge.agmax) goto 115
         enddo
         if (iterh.eq.1) then
            if (zgradmu) then
               write (nout,2001)
               write (nout,2002)
               write (nout,2003) (agdj(i),agg,i = 1,neqj)
            else
               write (nout,2004)
               write (nout,2005)
               write (nout,2006) (agdj(i),agg,i = 1,neqj)
            endif
         endif
         if (imodpr.gt.10) then
            if (zgradmu) then
               write (nout,2007) iterh,(jgd(i),gdj(i),jgg(i),gg(i),i = 1
     &              ,neqj)
            else
               write (nout,2008) iterh,(jgd(i),gdj(i),jgg(i),gg(i),i = 1
     &              ,neqj)
            endif
         endif
         omegaconv = omegaconv_sav
         goto 77

 115     error = 18
         if (zgradmu) then
            write (nout,2001)
            write (nout,2002)
            write (nout,2007) iterh,(jgd(i),gdj(i),jgg(i),gg(i),
     &           i = 1,neqj)
         else
            write (nout,2004)
            write (nout,2005)
            write (nout,2008) iterh,(jgd(i),gdj(i),jgg(i),gg(i),
     &           i = 1,neqj)
         endif
         print *,'max correction overwhelmed'
         return
      else
         if (zgradmu) then
            write (nout,2001)
            write (nout,2002)
            write (nout,2007) iterh,(jgd(i),gdj(i),jgg(i),gg(i),
     &           i = 1,neqj)
         else
            write (nout,2004)
            write (nout,2005)
            write (nout,2008) iterh,(jgd(i),gdj(i),jgg(i),gg(i),
     &           i = 1,neqj)
         endif
         error = 18
         print *,'max number of iterations reached'
         return
      endif

 2007 format (i4,5('|',2(i4,1x,1pe9.2)),'|')
 2008 format (i4,4('|',2(i4,1x,1pe9.2)),'|')

 2001 format (3x,'#','|',6x,'Variable om(m)',8x,'|',5x,
     &     'Variable aux(m)',8x,'|',5x,'Variable lambda(m)',5x,'|',6x,
     &     'Variable ur(m)',8x,'|',3x,'Variable xpsi(m)',9x,'|')
 2002 format (4x,5('|',1x,'convergence',2x,'|',3x,'quality',3x),'|')
 2003 format (4x,5('|',2x,1pe9.2,3x,'|',2x,1pe9.2,2x),'|',/)
 2004 format (3x,'#','|',6x,'Variable om(m)',8x,'|',5x,
     &     'Variable aux(m)',8x,'|',6x,'Variable ur(m)',8x,'|',3x,
     &     'Variable xpsi(m)',9x,'|')
 2005 format (4x,4('|',1x,'convergence',2x,'|',3x,'quality',3x),'|')
 2006 format (4x,4('|',2x,1pe9.2,3x,'|',2x,1pe9.2,2x),'|',/)

      end


** ---------------------------------------------------------------------
*
      SUBROUTINE NORM(om,ur,aux,xpsi,xlambda,vom,vxpsi,vrray,ivar,ndb,
     &     ndt)
*
** ---------------------------------------------------------------------
* Normalise les variables
*
* A noter que la vitesse de rotation est normalisee par le
*   omega en surface au pas precedant, alors que le rayon et
*   la gravite sont normalises par rapport au modele actuel
*
* Auteur: S.Talon
*
* Adaptation: A.Palacios (14/03/2000)
*
* Derniere version: 03 mai 2000
*
* Cas ivar = 3 :  premier modele avec rotation si on part d'un profil
*                 differentiel de omega (idiffvr = 6)
*
* Cas ivar = 2 : cas general pour initialisation des profils et
*                normalisation de omega et U
* ----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.transp'

      integer j,irotbin
      integer model,modeli,model_old
      integer ndt,ivar,ndb
      integer nphase
      integer klcore,klenv,klpulse

      double precision ur_smax,ur_smin
      double precision om,ur,aux,xpsi,xlambda
      double precision vom,vrray,vxpsi
      double precision divis,divis3
      double precision rhmoy,epsmoy
      double precision disctime

      logical locconsAM
      logical dup3,agbphase,relaxpulse,thermalpulse, rgbphase,superagb
     $     ,flame,urca

      common /conserv/locconsAM

      common /outp/ irotbin,model,modeli,model_old
      common /moy/ rhmoy(nsh),epsmoy(nsh)
      common /disclocking/ disctime
      common /evopha/ nphase,dup3,agbphase,relaxpulse,thermalpulse,
     &     rgbphase,superagb,flame,urca
      common /overshoot/ klcore,klenv,klpulse

      dimension om(nsh),ur(nsh),aux(nsh),xpsi(nsh),xlambda(nsh)
      dimension vom(nsh),vrray(nsh),vxpsi(nsh)

      model_old = model-1

      if (ur_Ss.eq.0.d0) ur_Ss = 1.d0
      if (.not.(ivar.eq.2.or.ivar.eq.3)) then
         print *,'call to norm useless : stop norm'
         stop 'norm'
      endif

      divis = 1.d0/rtot
      omega_S = vom(nmod)
c      if (nphase.lt.2.and.time/sec.lt.disctime.and.
c    $     (idiffvr.eq.5.or.idiffvr.eq.8)) omega_S = breaktime
      divis3 = 1.d0/omega_S
      do j = 1,nmod
         hht(j) = hht(j)*divis
         hhp(j) = hhp(j)*divis
         vom(j) = vom(j)*divis3
      enddo

      if (ivar.eq.3) then

*     Cas de rotation differentielle de l'EC (idiffvr = 6)
*     dans le binaire initial de rotation
         do j = 1,nmod
            aux(j) = 1.d-50
            xlambda(j) = 0.d0
         enddo
         if (model.ne.model_old) then
            om(1) = vom(1)
            xpsi(1) = vxpsi(1)
            do j = 2,nmod
               om(j) = vom(j)*vrray(j)*vrray(j)/(rray(j)*rray(j))
               xpsi(j) = vxpsi(j)*gs/(rtot*omega_S**2)
            enddo
         else
            do j = 1,nmod
               ur(j) = urs(j)
            enddo
         endif

      else if (ivar.eq.2) then

c..    Initialisation des profils et normalisation de omega et U

         ur_smax = abs(urs(ndt))
         do j = ndt-1,1,-1
            ur_smax = max(ur_smax,abs(urs(j)))
         enddo
         ur_S = ur_smax*ur_Ss/2.d0
         ur_smin = 1.d-8
         ur_S = max(ur_smin,ur_S)

         if (model.ne.model_old) then
*Meilleure convergence dans le cas ou la circulation est presque bloquee
*par les gradients de mu
            om(1) = vom(1)
c..   On initialise omega (guess pour Newton-Raphson)
            do j = 2,nmod
               om(j) = vom(j)*vrray(j)*vrray(j)/(rray(j)*rray(j))
            enddo
            if (locconsAM) om(1) = om(2)
            do j = 1,nmod
               ur(j) = urs(j)
               aux(j) = auxs(j)
               xpsi(j) = xpsis(j)
               xlambda(j) = xlambdas(j)
            enddo
         else
            do j = 1,nmod
               ur(j) = urs(j)
               if (om(j).eq.0.d0) om(j) = om(j+1)
            enddo
         endif
      endif

      xnorm = rtot*omega_S*omega_S/(gs*gs)
      model_old = model
      coeff = rtot*rtot
      xLtot = coeff*omega_S
      coeff = 4.d0*pi*rtot*coeff
      xLloc = coeff

      call val_moy

cl
      do j = 1,nmod
         if (crz(j).le.1.or.j.ge.novlim(klenv,3)) urs(j) = 0.d0
c         if (crz(j).lt.0) urs(j) = 0.d0
         ur(j) = urs(j)
      enddo

      return
      end


** ---------------------------------------------------------------------
*
      SUBROUTINE VAL_MOY
*
** ---------------------------------------------------------------------
* Calcul des valeurs intermediaires des parametres de l'etoile
*
* Auteur: S.Talon
*
* Adaptation: A.Palacios (14/03/2000)
*
* Derniere version: 25 avril 2000
* ----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.transp'
      include 'evolcom.var'

      integer j,j1

      double precision rhmoy,epsmoy
      double precision vom,vrray,vxpsi
      double precision fac

      common /moy/ rhmoy(nsh),epsmoy(nsh)
      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)

      coeff = 4.d0*pi*rtot**3
      do j = 1,nmod
        rraym(j) = 0.d0
        rraym2(j) = 0.d0
        hnenn(j) = 0.d0
        vrraym(j) = 0.d0
        vrraym2(j) = 0.d0
c        xKt(j) = khi(j)/(cp(j)*rro(j))
      enddo

      do j = 2,nmod
         j1 = j-1
         rraym(j) = sqrt((rray(j)**2+rray(j)*rray(j1)+rray(j1)**2)/3.d0)
         vrraym(j) = sqrt((vrray(j)**2+vrray(j)*vrray(j1)+vrray(j1)**2)
     &        /3.d0)
         rraym2(j) = rraym(j)*rraym(j)
         vrraym2(j) = vrraym(j)*vrraym(j)
         hnenn(j) = -1.d0/(rray(j)-rray(j1))
c..  interface variables
         gravm(j) = 0.5d0*(grav(j)+grav(j1))
         hhtm(j) = 0.5d0*(hht(j)+hht(j1))
         xlumlm(j) = 0.5d0*(lum(j)+lum(j1))
         rmassm(j) = 0.5d0*(m(j)+m(j1))
         rkonvLm(j) = 0.5d0*(rkonvL(j)+rkonvL(j1))
         rkonvm(j) = 0.5d0*(rkonv(j)+rkonv(j1))
         rrog(j) = rraym2(j)/gravm(j)
         xnumm(j) = 0.5d0*(xnum(j)+xnum(j1))
c..  centered variables
         xKtm(j) = wj(j)*xKt(j)+wi(j)*xKt(j1)
         khitm(j) = wj(j)*khit(j)+wi(j)*khit(j1)
         khi_mum(j) = wj(j)*khi_mu(j)+wi(j)*khi_mu(j1)
         rdensm(j) = wj(j)*rro(j)+wi(j)*rro(j1)
         epsitm(j) = wj(j)*epsit(j)+wi(j)*epsit(j1)
         epsimum(j) = wj(j)*epsimu(j)+wi(j)*epsimu(j1)
      enddo
      gravm(2) = grav(1)
      hhtm(2) = hht(1)
      rrog(2) = rraym2(2)/gravm(2)

      xKtm(1) = xKt(1)
      rmassm(1) = 0.d0
      xlumlm(1) = 0.d0
      rdensm(1) = rro(1)

      rhmoy(1) = 0.d0
      epsmoy(1) = 0.d0
      fac = pi*rtot**3
      if (idiffty.eq.13) then
         do j = 2,nmod
c            epsmoy(j) = epsmoy(j-1) + dm(j-1)*(enucl(j-1)+egrav(j-1))
            rhmoy(j) = rmassm(j)*pw34/(fac*rraym(j)**3)
            epsmoy(j) = lum(j)/m(j)
         enddo
      else
         do j = 2,nmod
c            epsmoy(j) = epsmoy(j-1) + dm(j-1)*enucl(j-1)
            rhmoy(j) = rmassm(j)*pw34/(fac*rraym(j)**3)
            epsmoy(j) = lum(j)/m(j)
         enddo
c         do j = 2,nmod
c            epsmoy(j) = epsmoy(j)/m(j)
c         enddo
      endif

      return
      end


** ---------------------------------------------------------------------
*
      SUBROUTINE ASSYMPT (ndt,neqj,ncls,nclc,times,error)
*
** ---------------------------------------------------------------------
* Cette routine permet de calculer le profil de rotation
*   asymptotique (pour un modele qui n'evolue plus)
*   (option IDIFFTY = 9)
*
* ASYMPT est appelee par DIFFUSION,
*            et appelle NORM, LISS, LAMBDA, HEN_OMEGA, CONS_L, DESMULT61
*
* ATTENTION : - Pour appeler cette routine, il faut avoir fait
*               evoluer le modele avec IDIFFTY = 8 avant
*             - Le modele de l'etoile ne tient pas compte de la
*               loi de rotation finale dans son equilibre mecanique
*
* Auteur : S.Talon
*
* Adaptation: A.Palacios
*
* Derniere version : 6 octobre 1995
* ----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.data'
      include 'evolcom.diff'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.transp'

      double precision ur_smax
      double precision fluxshe,fluxcir
      double precision vom,vrray,vxpsi
      double precision om,ur,times,xpsi,aux,xlambda
      double precision xmoment_tot
      double precision alph,ddmax
      double precision dtnmax,dtas
      double precision xflux_shear,xflux_cir
      double precision xlambdaold
      double precision rhmoy,epsmoy
      double precision dift,dife
      double precision dif
      double precision omega,vomega,k2conv,k2rad,vsurf,
     &     angradr,angconvr
      double precision xjmod

      integer ndb,ndt,ideux,error
      integer i,ii,j,neqj,ncls,nclc,pcol

      logical convergence

      common /lambdaold/ xlambdaold(nsh)
      common /difcirc/ dift(nsh),dife(nsh)
      common /fluxomega/ fluxshe(nsh),fluxcir(nsh)
      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)
      common /moy/ rhmoy(nsh),epsmoy(nsh)
      common /moment/ xmoment_tot
      common /rotvar/ omega(nsh),vomega(nsh),k2conv,k2rad,vsurf,
     &     angradr,angconvr

      dimension om(nsh),aux(nsh),ur(nsh),xpsi(nsh),xlambda(nsh)
      dimension dif(nsh),xjmod(nsh,neqj)
      dimension pcol(neqj)


      ideux = 2
      dtnmax = 10.d0*dtn
      dtas = dtn

c..   Initialisation du pointeur de colonnes pour ordonnancement
c..   des matrices

      if (zgradmu) then
         pcol(1) = 1
         pcol(2) = 2
         pcol(3) = 4
         pcol(4) = 5
         pcol(5) = 3
      else
         pcol(1) = 1
         pcol(2) = 2
         pcol(3) = 3
         pcol(4) = 4
      endif


* Normalisation des variables
      call norm (om,ur,aux,xpsi,xlambda,vom,vxpsi,vrray,ideux,ndb,ndt)

      do i = 1,nmod
         xpsi(i) = xpsis(i)
      enddo
      xmoment_tot = xmom_tots

* Premiere iteration ; comme dans difom

      do i = 1,nmod
         xjmod(i,pcol(1)) = om(i)
         xjmod(i,pcol(2)) = aux(i)
         xjmod(i,pcol(3)) = ur(i)
         xjmod(i,pcol(4)) = xpsi(i)
         if (zgradmu) xjmod(i,pcol(5)) = xlambda(i)
      enddo

      alph = 0.15d0
c      if (zgradmu) call sgsmooth (xlambdaold,ndt,1,17,21)

      call hen_omega (ndb,ndt,neqj,ncls,nclc,pcol,dtas,alph,times,xjmod,
     &     error)
      if (error.gt.0) return

      do i = 1,ndt
         om(i) = xjmod(i,pcol(1))
         aux(i) = xjmod(i,pcol(2))
         ur(i) = xjmod(i,pcol(3))
         xpsi(i) = xjmod(i,pcol(4))
         if (zgradmu) xlambda(i) = xjmod(i,pcol(5))
      enddo

      do j = ndt+1,nmod
         om(j) = om(ndt)
         ur(j) = 0.d0
         V_circ(j) = 0.d0
      enddo
      call cons_l (om,xmoment_tot,ideux,ndb,ndt)
      do j = 1,ndt
         dif(j) = om(j)-vom(j)
      enddo
      do j = 1,nmod
         vom(j) = om(j)/om(nmod)
         vrray(j) = rray(j)
      enddo
      omega_S = omega_S*om(nmod)
      ndtold = ndt

      if (zgradmu) then
c         call sgsmooth (xlambdaold,ndt,1,17,21)
         do i = 1,ndt
            xlambdaold(i) = xlambda(i)
         enddo
      endif

* Boucle d'iterations, sans modification de la structure de l'etoile
* (en particulier, la structure de l'etoile ne tiendra pas compte de
* la force centrifuge finale obtenue)
      convergence = .false.
      ii = 0
      alph = 0.15d0
      do while (.not.convergence)
         ii = ii+1
c         print *,ii
         do i = 1,nmod
            xjmod(i,pcol(1)) = om(i)
            xjmod(i,pcol(2)) = aux(i)
            xjmod(i,pcol(3)) = ur(i)
            xjmod(i,pcol(4)) = xpsi(i)
            if (zgradmu) xjmod(i,pcol(5)) = xlambda(i)
         enddo

         call hen_omega (ndb,ndt,neqj,ncls,nclc,pcol,dtas,alph,times,
     &        xjmod,error)

         if (error.gt.0) then
            write (nout,*) 'HEN_OMEGA non converge'
            goto 10
         endif

         do i = 1,ndt
            om(i) = xjmod(i,pcol(1))
            aux(i) = xjmod(i,pcol(2))
            ur(i) = xjmod(i,pcol(3))
            xpsi(i) = xjmod(i,pcol(4))
            if (zgradmu) xlambda(i) = xjmod(i,pcol(5))
         enddo

         if (dtas.lt.dtnmax) dtas = 1.1d0*dtas

         ur_Ss = ur_S
         ur_smax = abs(urs(ndt))
         do j = 1,ndt-1
            ur_smax = max(ur_smax,abs(urs(j)))
         enddo
         ur_S = ur_smax*ur_Ss*0.5d0
         do i = 1,ndt
            ur(i) = ur(i)*ur_Ss/ur_S
         enddo
         alph = 0.2d0
         do j = ndt+1,nmod
            om(j) = om(ndt)
         enddo
         ddmax = 0.d0
         do j = 1,ndt
            dif(j) = om(j)-vom(j)
            if (ddmax.lt.abs(dif(j))) ddmax = abs(dif(j))
         enddo

         if (zgradmu) then
            do i = 1,ndt
               xlambdaold(i) = xlambda(i)
            enddo
         endif
         do j = 1,nmod
            vom(j) = om(j)/om(nmod)
         enddo
         omega_S = omega_S*om(nmod)
         if (ddmax.lt.1.d-8) convergence = .true.
c     if (ddmax.lt.1.d-12) convergence = .true.
      enddo

      do i = 1,nmod
         omega(i) = om(i)*omega_S
         auxs(i) = aux(i)
         urs(i) = ur(i)
         xpsis(i) = xpsi(i)
         oms(i) = om(i)
         if (zgradmu) xlambdas(i) = xlambda(i)
      enddo

      do j = 1,ndt
*     attention: shear non denormalise, cir oui
c         if (zgradmu) then
            xflux_shear = 1.5d0*grav(j)*rray(j)*rray(j)*dift(j)*
     &           ((phiKSm(j)*xlambda(j)-deltaKSm(j)*xpsi(j))/om(j))
     &           *rro(j)*omega_S*rray(nmod)/(grav(nmod)*1.d36)
c         else
c            xflux_shear = -1.5d0*grav(j)*rray(j)*rray(j)*dift(j)*
c     &           deltaKSm(j)*xpsi(j)/om(j)*rro(j)*omega_S/
c     &           (grav(nmod)*1.d36)
c         endif
         xflux_cir = rro(j)*rray(j)**4*om(j)*ur(j)*ur_S
     &        *omega_S/5.d36*rtot**4
         fluxshe(j) = xflux_shear
         fluxcir(j) = xflux_cir
      enddo

 10   continue

      return
      end


** ---------------------------------------------------------------------
*
      SUBROUTINE ASSYMPTDV (ndt,xmoment_tot)
*
** ---------------------------------------------------------------------
*     Routine dans laquelle l'evolution du profil de omega n'est due
*     qu'aux reajustements structurels, et suppose une conservation
*     du moment specifique dans chaque couche i.e.: vr^2*vom = r^2*om
*     Le profil rotation est initialise a un profil differentiel
*     comme suggere par Denissenkov & VandenBerg (2003)
*
*     Appelee dans DIFFUSION si idiffty = 14
*
*     Le coefficient de diffusion applique aux elements chimiques
*     est le Dv de Maeder & Meynet (1996):
*     Dv = 2*K/Nt*(0.2*om^2*n^2-Nu^2)
*
*     Fait appel a la routine NORM
*     Est appelee par DIF_OMEGA
*
*     Auteur: A Palacios
*
* ----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.data'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.transp'
      include 'evolcom.eng'
      include 'evolcom.mod'

      double precision r,vr
      double precision vom,vrray,vxpsi
      double precision gdelmu,xNt,geffom
      double precision abmurj,abmuj
      double precision omega,vomega,k2conv,k2rad,vsurf,angradr,
     &     angconvr,xmoment_tot
      double precision dift,dife
c      double precision mom

      integer ndb,ndt,ivar
      integer j,j1,k
      integer klcore,klenv,klpulse

      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)
      common /gradmuj/ abmurj(nsh),abmuj(nsh)
      common /varrvr/ r(nsh),vr(nsh)
      common /rotvar/ omega(nsh),vomega(nsh),k2conv,k2rad,vsurf,
     &     angradr,angconvr
      common /difcirc/ dift(nsh),dife(nsh)
      common /overshoot/ klcore,klenv,klpulse


      write (nout,*) 'Conservation of the specific angular momentum',
     &     '  throughout the star'

***   During MS, specific angular momentum independent of time

c      if (nphase.le.2) then
         omega_S = vomega(nmod)
         do j = 2,nmod
            j1 = j-1
            omega(j) = vom(j)*vr(j)*vr(j)/(r(j)*r(j))
            oms(j) = omega(j)/omega_s
            vomega(j) = omega(j)
            vom(j) = vomega(j)
            rraym(j) = sqrt((rray(j)**2+rray(j)*rray(j1)+rray(j1)**2)/
     &           3.d0)
            vrraym(j) = sqrt((vrray(j)**2+vrray(j)*vrray(j1)+
     &           vrray(j1)**2)/3.d0)
            rraym2(j) = rraym(j)*rraym(j)
            vrraym2(j) = vrraym(j)*vrraym(j)
         enddo
c      else if (nphase.gt.2) then
c         do j = 2,ndt
c            omega(j) = vom(j)*vr(j)*vr(j)/(r(j)*r(j))
c            xKt(j) = khi(j)/(cp(j)*rro(j))
c            vomega(j) = omega(j)
c         enddo
c         mom = 0.d0
c         do j = ndt,nmod
c            mom = vom(j)*vr(j)**2*dm(j)+mom
c         enddo
c         do j = ndt,nmod
c            omega(j) = mom/(m(ndt+1)*r(j)**2)
c         enddo
         omega(1) = omega(2)
         vomega(1) = vomega(2)
         oms(1) = oms(2)
         vom(1) = vom(2)
c      endif

* Definition du coeffcient de diffusion turbulente vertical

      do k = 1,ndt
         gdelmu = -gmr(k)*phiKSm(k)*abmurj(k)
         geffom = deltaKSm(k)*(gmr(k)-pw23*r(k)*omega(k)**2)
         if (crz(k).le.1.or.k.ge.novlim(klenv,3)) then
            Dhd(k) = 1.d0
         else
            xNt = -geffom/hhp(k)*rkonv(k)
            Dhd(k) = 2.d0*xKtm(k)/xNt*(0.2d0*omega(k)**2*4.d0-gdelmu)
            Dhd(k) = max(Dhd(k),1.d0)
         endif
         dift(k) = Dhd(k)
      enddo

      vsurf = omega(nmod)*1.d-5*rtot

      if (inits) then
         ivar = 1
         call cons_l (oms,xmoment_tot,ivar,ndb,ndt)
         xmom_tots = xmoment_tot
      endif

      return
      end


** ---------------------------------------------------------------------
*
      SUBROUTINE ASSYMPT3 (ndt,om,ur)
*
** ---------------------------------------------------------------------
*     Routine permettant de calculer les coefficients de diffusion
*     necessaires au calcul du transport des especes chimiques
*     lorsque l'on ne tient plus compte du transport de moment cinetique.
*
*     Appelee a la fin de la sequence principale.
*
*     Equivaut au traitement applique par Meynet et Maeder dans leurs
*     articles.
*
*     Fait appel aux routines CDIFF,CAL_DH,TURBU_DHOLD et NORM
*     Est appelee par DIF_OMEGA
*
*     Auteur: A Palacios
*
* ----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.data'
      include 'evolcom.diff'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.transp'

      double precision vom,vrray,vxpsi
      double precision om,ur,xpsi,aux,xlambda
      double precision omm,omm1
      double precision xmoment_tot
      double precision dift,dife,mult
      double precision Rec,xnuvvm
      double precision xnurad,xnumol,nturb
      double precision Deff,Dhold

      integer ndb,ndt,ideux
      integer i,i1,j,j1

      common /calcDh/ Dhold(nsh)
      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)
      common /moment/ xmoment_tot
      common /difcirc/ dift(nsh),dife(nsh)
      common /dturbulente/ nturb(0:nsh)
      common /diffvisc/ xnumol(nsh),xnurad(nsh)

      dimension om(nsh),aux(nsh),ur(nsh),xpsi(nsh),xlambda(nsh)
      dimension xnuvvm(nsh)
      dimension Deff(nsh)

      write (nout,*)'Stationary treatment of meridional circulation'

      ideux = 2

* Normalisation des variables et definition de xpsi
      call norm (om,ur,aux,xpsi,xlambda,vom,vxpsi,vrray,ideux,ndb,ndt)

      do i = 2,ndt
         i1 = i-1
         if (i1.eq.1) then
            omm1 = om(1)
         else
            omm1 = 0.5d0*(om(i1)+om(i1-1))
         endif
         omm = 0.5d0*(om(i)+om(i1))
         xpsi(i) = (phiKSm(i)*xlambda(i)-(hnenn(i)*pw23*rrog(i)*
     &        (omm-omm1))*om(i))/deltaKSm(i)
      enddo
      xpsi(1) = xpsi(2)
      xmoment_tot = xmom_tots

* Definition des coeffs de diffusion turbulente horizontal et vertical
      call cal_dh (ndb,ndt,ur,om,xpsi,xlambda,1)
      call turbu_dhold (ndb,ndt,xpsi,xlambda,om,ur)

* Redefinition du coefficient de diffusion turbulente vertical en fonction
* de la stabilite vis-a-vis du critere de Reynolds
      xnuvvm(1) = xnuvv(1)
      do j = ndt,2,-1
         j1 = j-1
         xnuvvm(j) = 0.5d0*(xnuvv(j)+xnuvv(j1))
         if (j1.ge.2) xnuvvm(j1) = 0.5d0*(xnuvv(j1)+xnuvv(j1-1))
         Rec = (3.d0+pw13)*xnum(j)
         if ((xnuvvm(j).le.Rec).and.(rkonvm(j).lt.0.d0)) then
            nturb(j) = 0.d0
            if (j1.gt.1.and.(xnuvvm(j1).le.Rec).and.
     &           (rkonvm(j1).lt.0.d0)) then
               nturb(j1) = 0.d0
            else
               nturb(j1) = 1.d0
            endif
         else
            nturb(j) = 1.d0
            if (j1.gt.1.and.(xnuvvm(j1).le.Rec).and.(rkonvm(j1).
     &           lt.0.d0)) then
               nturb(j1) = 0.d0
            else
               nturb(j1) = 1.d0
            endif
         endif
         mult = min(nturb(j),nturb(j+1))
         Dhd(j) = xnumol(j)+mult*xnuvv(j)
         dift(j) = Dhd(j)
         if (j.gt.1) dift(j1) = Dhd(j1)
      enddo
      

* Definition de la vitesse de circulation meridienne dans le cas
* stationaire (Eq. (5) Meynet & Maeder, 2000, A&A 361, pp101-120)
      do i = 1,ndt-1
         i1 = i+1
         ur(i1) = -5.d0*dift(i1)/om(i1)*hnenn(i1)*(om(i)-om(i1))/
     &        (rtot*ur_S)
c         ur(i1) = -5.d0*xnuvv(i1)/om(i1)*hnenn(i1)*(om(i)-om(i1))/
c     &        (rtot*ur_S)
      enddo


*------------------------------------------------------
***   Computation of Deff for diffusion of chemicals
*------------------------------------------------------

      do i = 1,nmod
         Deff(i) = 0.d0
         dife(i) = 0.d0
      enddo

      call cal_dh (ndb,ndt,ur,om,xpsi,xlambda,2)
         
      do i = 2,ndt
         i1 = i-1
         if (Dhold(i).gt.0.d0) then
            Deff(i) = (rray(i)*ur(i)*ur_S*rtot)**2/(Dhold(i)*30.d0)
         else
            Deff(i) = 0.d0
         endif
         Dhd(i) = Dhd(i)+Deff(i)
         dife(i) = Deff(i)
      enddo
      Dhd(1) = Dhd(1)+Deff(2)
      dife(1) = Deff(2)
      Dhd(ndt) = Dhd(ndt)+Deff(ndt)
      dife(ndt) = Deff(ndt)

      return
      end


** ---------------------------------------------------------------------
*
      SUBROUTINE ROT_SOL (ndt,om,aux,ur,xpsi,vom,times)
*
** ---------------------------------------------------------------------
* Calcul de la vitesse de circulation et de Deff pour une rotation
*   uniforme
*
* Auteur: S.Talon
*
* Adaptation: A.Palacios (14/03/2000)
*
* Derniere version: 25 avril 2000
* ----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.grad'
      include 'evolcom.mass'
      include 'evolcom.rot'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.transp'
      include 'evolcom.var'

      integer ndt,ndb
      integer i,iun,j,j1,k
      integer klenv,klcore,klpulse

      double precision om,aux,ur,xpsi,vom,times
      double precision xmoment,xmoment_tot,momenvel,xmoment_vieux
      double precision rhmoy,epsmoy,rhmoym
      double precision divis,divis1
      double precision xmom,dmom
      double precision xom,vit
      double precision aalpha,bbeta
      double precision raprho,rapeps,tamp4
      double precision ur_smax
      double precision xmoment1,xmoment2,omm
      double precision Deff,Dhold
      double precision dift,dife
      double precision omega_sat

      double precision Bcrit,fstar,Bequi
      double precision expm,expr,ff,fourm,cte1,cte2
      double precision KK1,KK2,aexp,bexp,mexp,Bsat,Bsun,omsol,Mdotsun
     $     ,Mdotsat,B,mdot
      double precision chi,pexp,twom,TeffSUN
      double precision tauc,taucSUN,taucnorm
      double precision faccrit,btorq,torq0
      double precision tc_hp,rc_hp,tg
      
      common /calcDh/ Dhold(nsh)
      common /moy/ rhmoy(nsh),epsmoy(nsh)
      common /moment/ xmoment_tot
      common /momreel/ xmoment,momenvel,xmoment_vieux
      common /difcirc/ dift(nsh),dife(nsh)

      common /boreas_Mloss/ Bcrit,fstar,Bequi
      common /tauc_hp/ tauc


      dimension Deff(nsh)
      dimension om(nsh),aux(nsh),ur(nsh),xpsi(nsh),vom(nsh)


      divis1 = 1.d0/rtot
* Normalisation des variables
      do j = 1,nmod
         ur(j) = 0.d0
         aux(j) = 1.d-40
         xpsi(j) = 1.d-40
         hht(j) = hht(j)*divis1
         hhp(j) = hhp(j)*divis1
         om(j) = 1.d0
* On initialise un omega deja normalise
         vom(j) = 1.d0
      enddo
      coeff = 4.d0*pi*rtot**3
      call val_moy

      if (xmom_tots.eq.0.d0.or.omega_S.eq.0.d0) then
         omega_S = vomega(nmod)
         xLtot = rtot*rtot*omega_S
         iun = 1
         call cons_l (om,xmoment_tot,iun,ndb,ndt)
         xmom_tots = xmoment_tot
      else if (idiffvr.eq.5.or.idiffvr.eq.0.or.idiffvr.eq.6.or.
     &        idiffvr.ge.8) then
           if (idiffvr.eq.5) then
c              omega_sat = 10.d0*omsun !2.8d-5
              omega_sat = om_sat*omsun
              if (omega_S.lt.omega_sat) then
                 dmoment = -diffcst*dtn*omega_S**3*sqrt(solr/totm) 
                 xmoment_tot = xmom_tots+dmoment
                 xmom_tots = xmoment_tot
              else
                 dmoment = -diffcst*dtn*omega_S*omega_sat**2
     &                    *sqrt(solr/totm) 
                 xmoment_tot = xmom_tots+dmoment
                 xmom_tots = xmoment_tot
              endif
              dmoment1 = dmoment
***   Matt et al. 2012  equation 9 using outputs of the boreas routine 
***   from Cranmer & Saar 2011
           else if (idiffvr.eq.8) then
                  KK1 = 1.8d0

                  if (KK1.lt.0.5d0) KK1 = 0.5d0 
                  KK2 = 0.0506d0
                  mexp = 0.22d0
                  fourm = 4.d0*mexp
                  expm = 1.d0-2.d0*mexp
                  expr = 5.d0*mexp+2.d0
                  aexp = 1.3d0
                  bexp = 1.65d0
                  ff = 0.1d0
c                  ff = (om(ndt)*omega_S)**2.d0*27.d0*
c     &                 (rray(nmod)*rsun)**3.d0/(8.d0*g*(m(nmod)))
c                  ff = dsqrt(ff)
                  cte1 =  KK1**2/(2.d0*g)**mexp
                  cte2 = 1.d0/(KK2**2+0.5d0*ff**2)**mexp
                  mdot = dms*msun/sec
c                  B = Bcrit!*fstar
                  B = 1.13d0*Bequi*fstar
                  dmoment = -dtn*cte1*B**fourm*mdot**expm*
     $                 (rray(nmod)*rtot)**expr/(m(nmod))**mexp*
     $                 om(ndt)*omega_S*cte2
                  dmoment1 = dmoment
                  xmoment_tot = xmom_tots+dmoment
                  xmom_tots = xmoment_tot

***   Matt et al. 2015
           elseif (idiffvr.eq.9) then
              chi = om_sat
              mexp = 0.22d0
              !!pexp = 1.7d0
              pexp = 2.1d0
              twom = -2.d0*mexp
              
              tg = 0.d0

              k = novlim(klenv,3)
              do while (r(k).lt.r(novlim(klenv,3))+
     &             0.5d0*hp(novlim(klenv,3)).and.k.lt.nmod)
                 k = k+1
              enddo
              
              if (k.ge.novlim(klenv,4).or.sconv(k).eq.0.d0) then
                 tauc = tg
              else
                 tc_hp = alphac*hp(k)/sconv(k)
                 rc_hp = r(k)
                 tauc = tc_hp
              endif
              tauc = 3.1424d2*exp(-teff/1952.5d0)*exp(-(teff/
     $             6250.d0)**18)+2.d-3

!!              tauc = tauc/8.64d4
              taucSUN = 12.90430912750d0
              taucnorm = tauc/taucSUN
              faccrit = om(ndt)*omega_S*dsqrt(r(nmod)/gmr(nmod)) !! Need to include the deformation by rotation
              omega_sat = chi*omsun/taucnorm
              btorq = dsqrt(1.d0+(faccrit/0.072)**2.d0)
              torq0 = diffcst*(rtot/rsun)**3.1d0*(m(nmod)/msun)**0.5*
     &             btorq**twom
              if (omega_S.lt.omega_sat) then                     
                 dmoment = -dtn*torq0*taucnorm**pexp
     &                *(omega_S/omsun)**(pexp+1.d0)
              else
                 dmoment = -dtn*torq0*(chi**pexp)*(omega_S/omsun)
              endif
              print *,'torq0,dtom,taucnorm',torq0,dtn,omega_S,taucnorm
              dmoment1 = dmoment
!!              print *,xmom_tots,dmoment
              xmoment_tot = xmom_tots+dmoment
              xmom_tots = xmoment_tot

              print *,"dmoment1=",dmoment1
              
           else if (ndt.lt.nmod.and.idiffvr.eq.6) then
              xmoment1 = 0.d0
              xmoment2 = 0.d0
              do j = 2,ndt
                 j1 = j-1
                 omm = 0.5d0*(om(j)+om(j1))
                 dmom = rraym2(j)*omm*dm(j1)
                 xmoment1 = xmoment1+dmom
              enddo
              do j = ndt+1,nmod
                 j1 = j-1
                 xmoment2 = xmoment2+dm(j1)
              enddo
              xmoment2 = xmoment2*om(ndt)*rraym2(ndt)**2
              xmoment_tot = (xmoment1+xmoment2)*omega_S*rtot*rtot
              xmom_tots = xmoment_tot
           else
              xmoment_tot = xmom_tots
           endif
           xmom = 0.d0
           do i = 2,nmod
              dmom = rraym2(i)*dm(i-1)
              xmom = xmom+dmom
           enddo
           xmom = xmom*rtot*rtot
           omega_S = xmoment_tot/xmom
           xmoment = xmoment_tot
      else
* times defini ds diffusion.f, il faut le transmettre a dif_omega
* et a o_rotsol
         call calc_omega (times,xom,vit)
         omega_S = xom
         xLtot = rtot*rtot*omega_S
         iun = 1
         call cons_l (om,xmoment_tot,iun,ndb,ndt)
         xmom_tots = xmoment_tot
      endif

* Calcul de la vitesse de circulation initiale
      xLtot = rtot*rtot*omega_S
      xLloc = coeff
      xnorm = rtot*omega_S*omega_S/(gs*gs)

      do j = 2,ndt
!! à vérifier pour les étoiles massives         
c         if (crz(j).le.1.or.j.ge.novlim(klenv,3)) then
         if (crz(j).lt.0) then
            aalpha = 0.d0
         else
            divis = -grav(j)*m(j)*cpm(j)*tm(j)*rkonvL(j)
            aalpha = lum(j)*pm(j)/(rdensm(j)*divis)*xnorm
         endif

* voir remarque dans ci (if inutile)
         rhmoym = 0.5d0*rhmoy(j)+0.5d0*rhmoy(j+1)
         raprho = rdensm(j)/rhmoym
         if (j.eq.nmod) then
            rapeps = 0.d0
         else
            if (enucl(j).gt.1.d-40) then
               rapeps = (wj(j)*enucl(j)+wi(j)*enucl(j-1))/epsmoy(j)
            else
               rapeps = 0.d0
            endif
         endif            

         bbeta = 2.d0*(pw43-raprho)
         tamp4 = 1.d0-omega_S*omega_S*2.384989d6/rdensm(j)-rapeps
         ur(j) = aalpha*bbeta*rray(j)/grav(j)*tamp4
      enddo

      ur(1) = 0.d0
      if (ndt.eq.nmod) ur(ndt) = 0.d0

      ur_smax = abs(ur(ndt))
      do j = 1,ndt-1
         ur_smax = max(ur_smax,abs(ur(j)))
      enddo
      ur_S = ur_smax*0.5d0
      do i = 1,ndt
         ur(i) = ur(i)/ur_S
      enddo

      do i = 1,nmod
         Dhd(i) = 0.d0
      enddo

      do i = 2,ndt
         if (Dhold(i).gt.0.d0) then
            Deff(i) = (rray(i)*ur(i)*ur_S*rtot)**2/(Dhold(i)*30.d0)
         else
            Deff(i) = 0.d0
         endif
         Dhd(i) = Dhd(i)+Deff(i)
         dife(i) = Deff(i)
      enddo
      Dhd(1) = Dhd(1)+Deff(2)
      dife(1) = Deff(2)
      Dhd(ndt) = Dhd(ndt)+Deff(ndt)
      dife(ndt) = Deff(ndt)


      return
      end



** ---------------------------------------------------------------------
*
      SUBROUTINE TURBU_DHOLD (ndb,ndt,xpsi,xlambda,om,ur)
*
** ---------------------------------------------------------------------
* Calcul de la viscosite turbulente
*                                2*(du/dz)^2
* On prend Ric = 1/4, et Dv = -------------------
*                             Nt^2/(K+Dh)+Nu^2/Dh
*
* ou (du/dz)^2 = (r*sin(teta)*d omega/dr)^2, et
* _________
* (du/dz)^2 = 4/5*(r*d omega/dr)^2
*
*
* Mai2004: Estimation du coefficient de diffusion pour l'instabilite GSF
*
*
* Adaptation: A.Palacios (16/03/2000)
*
* Derniere version: 25 avril 2000
* ----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.grad'
      include 'evolcom.mod'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.transp'
      include 'evolcom.var'

      integer iend
      integer k,kp1,i,i1
      integer ndb,ndt
      integer erreur

      double precision xpsi,om,ur,xlambda,thetaj
      double precision Dhold
      double precision xNt,xNu,xNr,geffom
      double precision gdelmu
      double precision exc,dexc_l,dexc_p,dexc_o
      double precision divis,lognundt1,lognundt2
      double precision abmurj,abmuj
      double precision fac,enuclm,tamp
      double precision dGSF,vGSF,xHj,vES,dES,vgmu
      double precision glim,gmm,a1,a2,a3,a4
      double precision pmu(nmod),pu(nmod),aa(4)
      double precision z(2,3),w(8),xx(3)

      double precision dexc_p2,dexc_l2,dexc_o2
      double precision tamp1,tamp2,tamp3
      double precision RiPr,poly,poly3,gacoef,alcoef,becoef
      double precision xnumol,xnurad
      double precision costhetac, thetac, Acoeff, Bcoeff,Ccoeff
      double precision coefomega, sol1, coefA1, coefB1, coefA2, coefC2
      double precision dshear,sol2
      double precision omm,omm1

      common /calcDh/ Dhold(nsh)
      common /nnnk/ xNt(nsh),xNu(nsh),xNr(nsh),geffom(nsh)
      common /gradmuj/ abmurj(nsh),abmuj(nsh)
      common /GSF/ dGSF(nsh),vES(nsh),dES(nsh)
      common /diffvisc/ xnumol(nsh),xnurad(nsh)
      common /chemdiffshear/ dshear(nsh)

      dimension xpsi(nsh),xlambda(nsh),om(nsh),vGSF(nsh),xHj(nsh),
     &     ur(nsh)
      dimension vgmu(nsh)

      iend = min(nmod-1,ndt)

      do k = 1,nmod
         if (inits) Dhold(k) = 0.d0
         xNt(k) = 0.d0
         xNu(k) = 0.d0
         xNr(k) = 0.d0
         vGSF(k) = 0.d0
         dGSF(k) = 0.d0
         vES(k) = 0.d0
         dES(k) = 0.d0
         vgmu(k) = 0.d0
         xnuvv(k) = 0.d0
         dnuvvo(k) = 0.d0
         dnuvvp(k) = 0.d0
         dnuvvl(k) = 0.d0
         geffom(k) = 0.d0
         dshear(k) = 0.d0
      enddo


c$$$c... For Ric = 1/6
c$$$      if (Dv_prescr.eq.'TZ97') then
c$$$         fac = 0.24d0
c$$$      else if (Dv_prescr.eq.'Za92') then
c$$$         fac = 0.05d0
c$$$      endif
c$$$c... For Ric = 1/4
c$$$c      if (Dv_prescr = 'TZ97') then
c$$$c         fac = 0.45d0
c$$$c      else if (Dv_prescr = 'Za92') then
c$$$c         fac = 0.05d0
c$$$c      endif
c... For Ric = 1/6
      if (Dv_prescr.eq.'TZ97') then
         fac = 2.d0/15.d0
      else if (Dv_prescr.eq.'Za92') then
         fac = 1.d0/45.d0
      endif
c... For Ric = 1/4
c      if (Dv_prescr = 'TZ97') then
c         fac = 0.2d0
c      else if (Dv_prescr = 'Za92') then
c         fac = 1.d0/30.d0
c      endif

c      Essais pour la breche du lithium dans le code de Grenoble
***   Dv according to Maeder & Meynet, 2001, A&A 373
c         fac = 1.75669d0
***   Dv according to Talon & Zahn, 1997
c      fac = 0.45d0 ! Ric = 1/4 
c      fac = 0.3d0 ! Ric = 1/6

c..NEW TEST
      do k = max(2,ndb),iend
         kp1 = k+1
         gdelmu = -gs*grav(k)*abmurj(k)*phiKSm(k)

         geffom(k) = gs*grav(k)-pw23*rray(k)*rtot*(om(k)*omega_S)**2
         if (rkonv(k).ge.-del_zc) then
            xNt(k) = geffom(k)/hhp(k)/rtot*del_zc*deltaKSm(k)
         else
            xNt(k) = -geffom(k)/hhp(k)/rtot*rkonv(k)*deltaKSm(k)
         endif
         xNu(k) = -gdelmu

         xNr(k) = -hnenn(kp1)*omega_S**2*((om(kp1)*rray(kp1)*
     &        rray(kp1))**2-(om(k)*rray(k)*rray(k))**2)/rraym(k)**3
c$$$         tamp = (phiKSm(k)*xlambda(k)-deltaKSm(k)*xpsi(k))/om(k)
c$$$         if (tamp.ne.0.d0.and.rray(k).ne.0.d0) then 
c            exc = tamp**2*(grav(k)*omega_S/rray(k))**2
c$$$            exc = tamp**2*(geffom(k)*omega_S/(gs*rray(k)))**2
c$$$            dexc_l = 2.0d0*exc/tamp*phiKSm(k)/om(k)
c$$$            dexc_p = -2.0d0*exc/tamp*deltaKSm(k)/om(k)
c$$$            dexc_o = -2.d0*exc/om(k)
         omm = 0.5d0*(om(k)+om(k-1)) 
         omm1 = 0.5d0*(om(k)+om(k+1)) 
         if (omm.ne.0d0.and.omm1.ne.0.d0.and.rray(k).ne.0.d0) then
            exc = (rray(k)*(omm1-omm)*omega_s/
     &           ((rraym(k+1)-rraym(k))))**2
            dexc_l = 0.d0
            dexc_p = 0.d0
            dexc_o = 0.d0
         else
            exc = 0.d0
            dexc_l = 0.d0
            dexc_p = 0.d0
            dexc_o = 0.d0
         endif
         
         if (Dv_prescr.eq.'TZ97') then
            if (Dhold(k).gt.0.d0) then
               divis = xNt(k)/(xKtm(k)+Dhold(k))+xNu(k)/Dhold(k)
            else
               divis = xNt(k)/xKtm(k)
            endif

            tamp = 2.d0*fac/divis
            if (divis.gt.0.d0) then
               xnuvv(k) = tamp*exc
               dnuvvp(k)= tamp*dexc_p
               dnuvvl(k)= tamp*dexc_l
               dnuvvo(k)= tamp*dexc_o
               dshear(k) = 0.6d0*exc/divis
            endif


         else if (Dv_prescr.eq.'Za92') then
c$$$            fac = 0.05d0
            fac = 1.d0/45.d0
            divis = xNt(k)/xKtm(k)
            tamp = 2.d0*fac/divis
            if (divis.gt.0.d0) then
               xnuvv(k) = tamp*exc
               dnuvvp(k)= tamp*dexc_p
               dnuvvl(k)= tamp*dexc_l
               dnuvvo(k)= tamp*dexc_o
               dshear(k) = 0.125d0*exc/divis
            endif

         else if (Dv_prescr.eq.'Pr16') then
            if (tamp.ne.0.d0.and.
     &       rray(k).ne.0.d0.and.
     &        xNt(k).ne.0.d0) then 

               alcoef = 3.34e-2
               becoef = 1.88d1
               gacoef = -2.86d3
               coefomega = xKtm(k)/xNt(k)*9.d0/4.d0*(grav(k)*omega_S
     $           /rray(k))**2*tamp**2
               Acoeff = alcoef*coefomega
               Bcoeff = becoef*xnumol(k)
               Ccoeff = gacoef*(xnumol(k))**2/coefomega
               
               sol1 = 0.5d0*(-Bcoeff + dsqrt(Bcoeff**2-4.d0*Acoeff
     $              *Ccoeff))/Acoeff

               sol1 = min(sol1,1.d0)
               sol1 = max(sol1,0.d0)

               thetac = dasin(dsqrt(sol1))
               if (thetac.ge.0.5d0*pi) then
                  costhetac = 0.d0
                  coefC2 = 0.d0
               else
                  costhetac = cos(thetac)
                  coefC2 = log(1.d0/tan(0.5*thetac))
               endif

               coefA1 = (costhetac-pw23*costhetac**3+0.2d0*costhetac**5)
               coefB1 = costhetac-pw13*costhetac**3
               coefA2 = coefB1
               
***   Expression de Dv
               xnuvv(k) =1.5d0*(Acoeff*coefA1+Bcoeff*coefB1+Ccoeff
     &              *costhetac)
***   Expression de la viscosité qui vaut 0.8 du Dv
               xnuvv(k) = 0.8d0*xnuvv(k)
               
               dnuvvl(k) = 0.8d0*3.0d0*(Acoeff*coefA1*phiKSm(k)/(om(k)*
     $              tamp)-Ccoeff*costhetac*phiKSm(k)/tamp*om(k)**2)
               dnuvvp(k) = 0.8d0*3.0d0*(-Acoeff*coefA1*deltaKSm(k)/ 
     $              (om(k)*tamp)+Ccoeff*costhetac*deltaKSm(k)/tamp*
     &              om(k)**2)
               dnuvvo(k) = 0.8d0*3.d0*(-Acoeff*coefA1/om(k)+Ccoeff*
     $              costhetac/om(k))
               
***   Expression pour la diffusion des especes chimiques
               
               dshear(k) = Acoeff*coefA2 + Bcoeff*costhetac + Ccoeff*
     $              coefC2
            else
               xnuvv(k) = 0.d0
               dnuvvp(k)= 0.d0
               dnuvvl(k)= 0.d0
               dnuvvo(k)= 0.d0
               dshear(k) = 0.d0
            endif

            if (dshear(k) /= dshear(k)) then
               print *,'TURBU_DHOLD : Problem at shell',k,dshear(k)
               dshear(k) = 0.d0
            endif
            if (xnuvv(k) /= xnuvv(k)) then 
               xnuvv(k) = 0.d0
               dnuvvp(k)= 0.d0
               dnuvvl(k)= 0.d0
               dnuvvo(k)= 0.d0
            endif
           
c            print *,'k,dshear=',k,dshear(k),'xnuvv=',xnuvv(k)
              
         else if (Dv_prescr.eq.'Ma97') then

c Dv from Maeder   
c with alpha = 1        
C*=====
C Dv from Maeder 1997
C
C Dv = Kt/(phi/delta*nablamu+nablaad-nabla) * alpha*Hp/(g*delta)*(9pi/32*om*dlnom/dlnr)**2
C
C with alpha = 1 for Ri_crit = 1
C
C*====             
c            divis = 2.d0*(xNt(k)+xNu(k))/(xKtm(k))
c         endif
*Ajout de *6 - identique a ce qui a ete fait pour l'etoile de 9 Mo
*Calibration de Dturb pour obtenir le bon Li du cote gauche de la breche
c        xnuvv(k) = 12.d0*exc/divis

c      Essais pour la breche du lithium dans le code de Grenoble
***   Dv according to Maeder & Meynet, 2001, A&A 373
c         fac = 1.75669d0

c..   Maeder 1997
c..   Resolution de l'equation cubique (6.47) 
            xx(1:3) = 0.d0
            if (rkonv(k).ge.-del_zc) then
               pmu(k) = phiKSM(k)*abmuj(k)/deltaKSm(k)/del_zc
            else
               pmu(k) = -phiKSM(k)*abmuj(k)/deltaKSm(k)/rkonv(k)
            endif            
            pu(k) = exc/xNt(k)
            aa(1) = 12.d0*pmu(k)
            aa(2) = 2.d0+2.d0*pmu(k)-6.d0*pu(k)
            aa(3) = 6.d0+2.d0*pmu(k)- pu(k)
            aa(4) = -pu(k)

c..   Appel routine NAG c02agf
            erreur = 0
c..   Sélection des valeurs de Gamma associées à (pu,pmu) >= (0,0)
c..   qui correspondent à l'instabilité de cisaillement en ZR
            if (aa(1).ne.0.d0) then
               call c02agf(aa,3,.true.,z,w,erreur)
               if (erreur.ne.0) then
                  write(*,*)'Erreur dans librairie NAG'
                  stop
               else
c..   Sélection des solutions réelles et positives (Gamma doit etre positif)
                  do i = 1,3
c                     if (z(2,i).eq.0.d0.and.z(1,i).gt.0.d0) 
                     if (z(2,i).eq.0.d0.and.z(1,i).gt.0.d0.and.
     $                    pmu(k).ge.0.d0.and. pu(k).ge.0.d0)
     &                    xx(i) = z(1,i)
                  enddo
c..   Si on est dans la zone de solutions multiple, prendre celles < 0.483
                  if (xx(1)*xx(2).ne.0.d0.or.xx(1)*xx(3).ne.0.d0.or.
     &                 xx(3)*xx(2).ne.0.d0) then
                     glim = 0.483d0
                     do i = 1,3
                        if (xx(i).lt.glim) gmm = xx(i)
                     enddo
c..   Sinon prendre la solution réelle unique trouvée
                  else 
                     do i = 1,3
                        if (xx(i).ne.0.d0) gmm = xx(i)
                     enddo
                  endif
                  xnuvv(k) = 2.d0*xKtm(k)*gmm
               endif
            endif
         endif
      enddo
c      print *,'xnuvv=',xnuvv(k)
      xnuvv(iend) = xnuvv(iend-1)
      dnuvvp(iend) = dnuvvp(iend-1)
      dnuvvl(iend) = dnuvvl(iend-1)
      dnuvvo(iend) = dnuvvo(iend-1)



* Protection aux limites des ZC
      if (xnuvv(ndt-1).ne.0.d0) then
         lognundt1 = log10(abs(xnuvv(ndt-1)))
      else
         lognundt1 = 0.d0
      endif
      if (xnuvv(ndt-2).ne.0.d0) then
         lognundt2 = log10(abs(xnuvv(ndt-2)))
      else
         lognundt2 = 0.d0
      endif
      if (lognundt1-lognundt2.gt.0.5d0) xnuvv(ndt-1) = xnuvv(ndt-2)

***   Instabilite GSF

      do i = 2,iend
         i1 = i-1
         if (-xNr(i).gt.xnum(i)*xNt(i)/xKtm(i)) then
            thetaj = (phiKSm(i)*xlambda(i)-deltaKSm(i)*xpsi(i))/om(i)
            xHj(i) = (2.d0*rray(i)*om(i)+1.5d0*grav(i)*thetaj)/
     &           (rtot*rray(i)**2*om(i))
            vGSF(i) = 2.d0*hhtm(i)*xHj(i)*1.5d0*grav(i)*thetaj/
     &           (rray(i)*om(i))*ur(i)*ur_S*rtot
            dGSF(i) = abs(vGSF(i)**2*rtot/(hnenn(i)*(vGSF(i)-
     &           vGSF(i1)+1.d-30)))
         endif

***   Eddington Sweet (d'apres Sofia & Endal)
         enuclm = wj(i)*enucl(i)+wi(i)*enucl(i-1)
         vES(i) = -abm(i)/(deltaKSm(i)*rkonv(i))*(om(i)*omega_S)**2*
     &        rray(i)**3*rtot**3*lum(i)/(g*m(i))**2*(2.d0*enuclm*
     &        rray(i)**2*rtot**2/lum(i)-2.d0*rray(i)**2*rtot**2/m(i)-
     &        pw34/(pi*rdensm(i)*rray(i)*rtot))

         vgmu(i) = hhp(i)*rtot**2*rray(i)*lum(i)/(g*m(i)**2)*
     &        phiKSm(i)*abmuj(i)/(deltaKSm(i)*rkonv(i))
         vES(i) = max(abs(vES(i))-abs(vgmu(i)),0.d0)
         if (vES(i)-vES(i-1).ne.0.d0) then
            dES(i) = abs(vES(i)**2*rtot/(hnenn(i)*(vES(i)-vES(i-1))))
         endif
      enddo

      return
      end



** ---------------------------------------------------------------------
*
**      SUBROUTINE ASH_RGB_PROF (ndt,i,om,omdom,xpsi) 
      SUBROUTINE ASH_RGB_PROF (ndt,om,omdom,xpsi) 
*
** ---------------------------------------------------------------------
* Angular velocity profile in the convective envelope of RGB stars
* according to 3D hydr simulations done with ASH.
* Brun & Palacios, 2009, ApJ
*
* Profile fit to case RGB2 (1/50th solar rotation) :
*
*               356 + 9.5 * r^(-0.044)
* Omega(r)  =  ------------------------ * contfac
*                 1.5 + 29.7 * r^(1.29)
*
*              
* Omega(r)  = contfac * (238.8 - 1471.43 * r + 3616.5 * r^2 - 3000.59 * r^ 3)
*
*
*
* Use a normalization factor so as the omega profile is continuous at
* the CE lower edge.
*
* Called  in SURFOMEGA when IDIFFVR = 7
*
*
* INPUT : ndt = last radiative shell before convective envelope
*
* OUTPUT: om = omega profile in the convective envelope according to ASH simulations
*         omdom = omega x domega
*
* Author : A.Palacios 
*
* Last version : 03/02/2011
* ----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'
      include 'evolcom.transp'
      include 'evolcom.teq'
      include 'evolcom.therm'

      double precision a1,a2,a3,a4,e1,e2,e1m1,e2m1,de1m1,de2m1,contfac
     $     ,omdom(nsh),om(nsh),xpsi(nsh)
      integer i,ndt
      
c$$$      a1 = 356.d0
c$$$      a2 = 9.5d0
c$$$      a3 = 1.5d0
c$$$      a4 = 29.7d0
c$$$      e1 = -0.044d0
c$$$      e2 = 1.29d0
c$$$      e1m1 = e1 - 1.d0
c$$$      e2m1 = e2 - 1.d0
c$$$      de1m1 = 2.d0*e1 - 1.d0
c$$$      de2m1 = 2.d0*e2 - 1.d0
c$$$      contfac = om(ndt)*(a3 + a4 * rray(ndt)**e2)/(a1+a2*rray(ndt)**e1)
c$$$      
c$$$      print *,'ndt in ash_rgb_prof',ndt
c$$$      omdom(1:nmod) = 0.d0
c$$$      do i = ndt,nmod
c$$$         om(i) = contfac * (a1 + a2 * rray(i)**e1)/(a3+a4*rray(i)**e2)
c$$$c         omdom(i) = contfac**2*((a3+a4*rray(i)**e2)**2*(2.d0*a1*a2
c$$$c     $        *e1*rray(i)**e1m1+a2**2*2.d0*rray(i)**de1m1)
c$$$c     $        -(a1+a2*rray(i)**e1)**2*(2.d0*a3*a4 *e2*rray(i)**e2m1+a4
c$$$c     $        **2*2.d0*rray(i)**de2m1))/(a3+a4*rray(i)**e2)**4
c$$$         omdom(i) = contfac**2*(a1+a2*rray(i)**e1)*((a3+a4*rray(i)**e2)
c$$$     $        *e1*a2*rray(i)**e1m1-(a1+a2*rray(i)**e1)*e2*a4*rray(i)
c$$$     $        **e2m1)/(a3+a4*rray(i)**e2)**3
c$$$      enddo
      a1 = -69.8425d0
      a2 = 72.4075d0
cs      e1 = -0.412881d0
      e1 = -0.112881d0
      e1m1 = e1 - 1.d0
c      print *,'om(ndt) beginning ash_prof',om(ndt)
      contfac = om(ndt)/(a1 + a2 * rray(ndt)**e1)
      
      omdom(1:nmod) = -1.5d0*grav(1:nmod)*deltaKS(1:nmod)*xpsi(1:nmod)
     $     /rray(1:nmod)**2
c      omdom(ndt) = contfac**2*(a1+a2*rray(i)**e1)*a2*e1*rray(i)**e1m1
      
      do i = ndt+1,nmod
         om(i) = contfac * (a1 + a2 * rray(i)**e1)
         omdom(i) = contfac**2*(a1+a2*rray(i)**e1)*a2*e1*rray(i)**e1m1
      enddo

c      do i=1,nmod
c         write(889,'(i4,1x,3(1pe14.7,1x))'),i,rray(i),om(i),omdom(i)
c      enddo

      
c$$$      a1 = 238.67862d0
c$$$      a2 = -1471.4282d0
c$$$      a3 = 3616.5004d0
c$$$      a4 = -3000.5877d0
c$$$      contfac = om(ndt) /(a1 + a2*rray(ndt) + a3*rray(ndt)**2 + a4
c$$$     $     *rray(ndt)**3)
c$$$      omdom(1:nmod) = 0.d0
c$$$      do i = ndt,nmod
c$$$         om(i) = (a1 + a2*rray(i) + a3*rray(i)**2 + a4 * rray(i)**3) *
c$$$     $        contfac
c$$$         omdom(i) = contfac * (a2 + 2.d0*a3*rray(i) + 3.d0*a4*rray(i)
c$$$     $        **2) * om(i)
c$$$
c$$$      enddo
      return
      end


**************************************
***                                  *
***   NUMERICAL ROUTINES             *
***                                  *
**************************************


*----------------------------------------------

      SUBROUTINE SGSMOOTH (func1,msht,mshb,nl,nrr)

*----------------------------------------------
* Permet de lisser une fonction f en utilisant
* des filtres de Savitzky-Golay
* Appelle les routines: SAVGOL,LUBKSB,LUDCMP
*----------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.teq'

      integer np,nl,nrr,mm,msht,mshb
      integer ndb0,nfin1
      integer i,j,k
      integer nl2,nr2

      double precision cc,func1,func2

      dimension cc(nl+nrr+1),func1(nsh),func2(nsh)

      np = nl+nrr+1
      mm = 2
      ndb0 = max(mshb,nl+1)
c     if (ndb0.lt.nl) ndb0 = nl+1
c     if (ndt0+nl.gt.nmod) ndt0 = nmod-nl

      do i = 1,nmod
         func2(i) = 0.d0
      enddo

      do j = mshb,nl
         nl2 = j-1
         nr2 = np-nl2-1
         call savgol (cc,np,nl2,nr2,0,mm)
         func2(j) = cc(1)*func1(j)
         do k = 1,nr2
            func2(j) = func2(j)+cc(k+1)*func1(j+k)
         enddo
         do k = 1,nl2
            func2(j) = func2(j)+cc(np-k+1)*func1(j-k)
         enddo
      enddo

      call savgol (cc,np,nl,nrr,0,mm)


      nfin1 = ndb0+int((msht-ndb0)/2)
      do j = ndb0,nfin1
         func2(j) = cc(1)*func1(j)
         do k = 1,nrr
            func2(j) = func2(j)+cc(k+1)*func1(j+k)
         enddo
         do k = 1,nl
            func2(j) = func2(j)+cc(np-k+1)*func1(j-k)
         enddo
      enddo
      do j = nfin1+1,msht
         func2(j) = cc(1)*func1(j)
         do k = 1,nl
           if ((j+k).le.nmod) then
              func2(j) = func2(j)+cc(np-k+1)*func1(j+k)
           else
              func2(j) = func2(j)+cc(np-k+1)*func1(nmod)
           endif
         enddo
         do k = 1,nrr
            func2(j) = func2(j)+cc(k+1)*func1(j-k)
         enddo
      enddo

      do i = mshb,msht
         func1(i) = func2(i)
      enddo

      return
      end


*----------------------------------------------------

      SUBROUTINE SGSMOOTH2 (func1,msht,mshb,nelem)

*----------------------------------------------------
* Permet de lisser une fonction f en utilisant
* des filtres de Savitzky-Golay
* Appelle les routines: SAVGOL,LUBKSB,LUDCMP
*----------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.teq'

      integer np,nl,nrr,mm,msht,mshb,nelem
      integer ndt0,ndb0,nfin2
      integer i,j,k,l

      double precision cc,func1,func2

      parameter (mm = 2, nl = 15, nrr = 20, np = nl+nrr+1)

      dimension cc(np),func1(nsh,nsp),func2(nsh,nsp)

      ndt0 = msht-1
      ndb0 = mshb+1
      if (ndb0.lt.nl) ndb0 = nl+1
      if (ndt0+nl.gt.nmod) ndt0 = nmod-nl

      call savgol (cc,np,nl,nrr,0,mm)

      do j = 1,nelem
         do i = 1,nmod
            func2(i,j) = 0.d0
         enddo

         do k = 1,ndb0+5
            func1(k,j) = 0.d0
         enddo

         do k = ndt0-5,nmod
            func1(k,j) = 0.d0
         enddo

         nfin2 = ndb0+int((ndt0-ndb0)/2)
         do l = ndb0,nfin2
            func2(l,j) = cc(1)*func1(l,j)
            do k = 1,nrr
               func2(l,j) = func2(l,j)+cc(k+1)*func1(l+k,j)
            enddo
            do k = 1,nl
               func2(l,j) = func2(l,j)+cc(np-k+1)*func1(l-k,j)
            enddo
         enddo
         do l = nfin2+1,ndt0
            func2(l,j) = cc(1)*func1(l,j)
            do k = 1,nl
               func2(l,j) = func2(l,j)+cc(np-k+1)*func1(l+k,j)
            enddo
            do k = 1,nrr
               func2(l,j) = func2(l,j)+cc(k+1)*func1(l-k,j)
            enddo
         enddo

         do i = 1,nmod
            func1(i,j) = func2(i,j)
         enddo
      enddo

      return
      end


*-------------------------------------------------------------------

      SUBROUTINE LISSAGE (yy,nmax,nmin,id)

*-------------------------------------------------------------------
* Lissage de la quantite y par moyenne arithmetique
* Auteur: S.Talon
* Derniere version: 20 septembre 1995
*-------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      integer nmin,nmax
      integer i,id

      double precision ttt,yy

      dimension yy(nsh)
      dimension ttt(nsh)

      if (id.eq.3) then
         do i = nmin,nmax
            ttt(i) = (3.d0*yy(i)+yy(i-1)+yy(i+1))/5.d0
         enddo
      elseif (id.eq.5) then
         do i = nmin,nmax
            if (i.ge.3) ttt(i)=(3.d0*yy(i)+2.d0*yy(i-1)+
     &           2.d0*yy(i+1)+yy(i+2)+yy(i-2))/9.d0
         enddo
         ttt(2) = (2.d0*yy(i)+yy(i-1)+yy(i+1))/4.d0
      endif
      do i = nmin,nmax
        yy(i) = ttt(i)
      enddo

      return
      end
