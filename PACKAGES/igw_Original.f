C ----------------------------------------------------------------------
      SUBROUTINE TRANSPORT_ONDES (om,ndb,ndt,dtom,iterondes)
C ----------------------------------------------------------------------
C Routine appelee par DIF_OMEGA dans jtranspc_ondesexcit.f
C
C C.Charbonnel & S.Talon (Fevrier 2003)
C Version incluant le calcul de l'excitation a chaque modele evolutif :
C C.Charbonnel (Decembre 2006)
C Version incluant les ondes produites par le coeur convectif
C S.Talon & C.Charbonnel (Octobre 2007)
!ORIGINAL VERSION
C ----------------------------------------------------------------------
      implicit none

      include 'evolpar.star'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.grad'
      include 'evolcom.teq'
      include 'evolcom.transp'
c      include 'evolcom.transpondesexcit'
      include 'evolcom.igw'
      include 'evolcom.therm'
      include 'evolcom.var'
      include 'evolcom.diff'

      double precision dzeit
      double precision om
      double precision xmoment_tot
      double precision dtom
C Modif CC (28/10/09) Passage de depotondestot dans evolcom.transpondesexcit
C      double precision depotwaves,depotondestot
c      double precision depotwaves
      double precision sumdepottot

      integer i, j, mtu, npasr, ncouche, iterondes
      integer ndt,ndb

      dimension om(nsh)!,depotwaves(nsh)
C Modif CC (28/10/09) Passage de depotondestot dans evolcom.transpondesexcit
C      dimension depotondestot(nsh)



      dzeit = dtom
      print *,'dzeit = ',dzeit,' dtom=',dtom,' dtn =',dtn
      print *,'entree dans TRANSPORT_ONDES, igwrot=',
     & igwrot,'igwsurfchim=',igwsurfchim

C --------------------------------------------
C Lecture des donnees et inversion des couches 
C --------------------------------------------

      ncouche = nmod
      mtu = nmod-ndt+1
      npasr = nmod-ndb+1

      do i = 1,nmod
         om_o(nmod-i+1) = om(i)
      enddo
          
      do i=1,nmod
         depotondestot(i)=0.d0
         depotwaves(i)=0.d0
      enddo

C ------------
C Excitation
C ------------
C      print *,'avant appel a WAVEFLUX'
      call waveflux(ncouche,mtu,npasr,dzeit)
 
C      print *,'sortie de WAVEFLUX'

C Calcul du depot du moment par les ondes
C----------------------------------------
C
      call ond_omega(iterondes,mtu,npasr,ncouche,dzeit,
     &               omega_S)
C
      print *,'sortie de OND_OMEGA'

C Modifs 3/10/07
      do i = mtu-1,npasr
         depotondestot(i) = 2.d0*depottot(i)/(rtot**2*omega_S)
C      write(*,*)'ncouche, i, depottot(i), depotondes(i) ',
C     & ncouche,i,depottot(i),depotondes(i)
      enddo 

C Test a revoir
C depotondes(mtu-1) --> ce qui est depose dans la ZC de surface
C depotondes(npasr) --> ce qui est depose dans la ZC centrale

      do i = 1,nmod-1
C Test CC 14 mai 2004 : depotwaves(i) = 1.d02*depotondes(nmod-i+1)
         depotwaves(i) = depotondestot(nmod-i+1)
C         depotwaves(i) = 1.d01*depotondestot(nmod-i+1)
      enddo 
C      print *,'1.d01*depotondestot'
C Fin test CC

      print *,'rtot=',rtot,' omega_S=',omega_S

      print *,'depottot(mtu)=',depottot(mtu),
     & 'depotondestot(mtu)=',depotondestot(mtu), 
     & ' depotwaves(ndt)=',depotwaves(ndt)

      print *,'depottot(mtu-1)=',depottot(mtu-1),
     & 'depotondestot(mtu-1)=',depotondestot(mtu-1), 
     & ' depotwaves(ndt+1)=',depotwaves(ndt+1)

      print *,'depottot(mtu+100)=',depottot(mtu+100),
     & 'depotondestot(mtu+100)=',depotondestot(mtu+100),
     & ' depotwaves(ndt-100)=',depotwaves(ndt-100)

      print *,'depotondestot(mtu+200)=',depotondestot(mtu+200),
     & ' depotwaves(ndt-200)=',depotwaves(ndt-200)

      print *,'depotondestot(npasr)=',depotondestot(npasr),
     & 'depotwaves(nmod-npasr+1)=',depotwaves(nmod-npasr+1)

      return
      end

C ---------------------------------------------------------------------
      SUBROUTINE WAVEFLUX(ncouche,mtu,npasr,dzeit)
C ---------------------------------------------------------------------
C Appelle la routine calculflux, qui appelle excitation pour m=1 et m=l
C  et interpole les autres valeurs de m
C Appellee par la routine transp_ondesexcit
C
C Derniere version: 9 septembre 2009
C ---------------------------------------------------------------------

      implicit none
      include 'evolpar.star'
      include 'evolcom.cons'
      include 'evolcom.transp'
c      include 'evolcom.transpondesexcit'
      include 'evolcom.igw'

      double precision freq_min,freq_max,dfreq,diff
C CC (9/9/9): r_ondes,r_ondes_core,iondes,iondes_core ne sont plus utilises
C             Voir iwave
C      double precision dzeit,delta_om,r_ondes,r_ondes_core,delta_norm
      double precision dzeit,delta_om,delta_norm
C      integer ll,ii,mm,ncouche,mtu,npasr,iondes,iondes_core
      integer ll,ii,mm,ncouche,mtu,npasr
      integer ll_xlum_max,mm_xlum_max

      write (*,*) 'Entree waveflux, mtu=',mtu

CC Test de calcul de delta_om au point ou Dondes devient faible
CC      et inferieur a xnurad
CC Ce test est a re-ecrire : Dondes est dans le sens STAREVOL, xnurad_o dans le sens Geneve
CC      if (igwsurfchim) then 
CC          iondes  = mtu
CC          do while (Dondes(iondes).gt.xnurad_o(iondes))
CC             iondes = iondes+1
CC          enddo
CC          r_ondes = rray_o(iondes)
CC      else
C          r_ondes = rray_o(mtu)-3.d-2
C	  r_ondes_core = rray_o(npasr)+3.d-2
CC      write (*,*) 'r_ondes, rray_o(mtu), mtu=',r_ondes,rray_o(mtu),mtu
CC      endif

C      iondes  = mtu
C      do while (rray_o(iondes).gt.r_ondes)
C        iondes = iondes+1
C      enddo
C      write (*,*) 'r_ondes, iondes, rray_o(mtu), mtu ',r_ondes,iondes,
C     &            rray_o(mtu),mtu

C      iondes_core  = npasr
C      do while (rray_o(iondes).lt.r_ondes_core)
C        iondes_core = iondes_core-1
C      enddo
C      write (*,*) 'r_ondes_core, iondes_core, ',r_ondes_core,iondes_core
     
C Calcul des flux pour les modes sectoriaux (l) et m=1 (1)
C --------------------------------------------------------
C calculflux remplace la lecture des tables de flux
C calcul pour l'enveloppe et le coeur

      call calculflux(ncouche,mtu,npasr,freq_min,freq_max)
 
      print *,'Dans Waveflux, lmax,nfreq,freq_min,freq_max',
     &         llmax,nfreq,freq_min,freq_max
 
      freq_min=freq_min*2.d0*pi*1.d-6
      freq_max=freq_max*2.d0*pi*1.d-6
      dfreq=(freq_max-freq_min)/(nfreq-1)
 
      do ii=1,nfreq
         freqonde(ii)=freq_min+(ii-1)*dfreq
      enddo

C Interpolation pour les autres modes
C -----------------------------------
      do ii=1,nfreq
      do ll=3,llmax
         diff = (xlum_ond(ii,ll,ll)-xlum_ond(ii,ll,1))
     &          /dfloat(ll-1)
         do mm=2,ll-1
            xlum_ond(ii,ll,mm)=xlum_ond(ii,ll,1) 
     &          + diff*dfloat(mm)
         enddo
         diff = (xlum_ond_core(ii,ll,ll)-xlum_ond_core(ii,ll,1))
     &          /dfloat(ll-1)
         do mm=2,ll-1
            xlum_ond_core(ii,ll,mm)=xlum_ond_core(ii,ll,1) 
     &          + diff*dfloat(mm)
         enddo
      enddo
      enddo
 

C Introduction de la valeur de coupure
C ------------------------------------
      do ii=1,nfreq
         do ll=1,llmax
            do mm=1,ll
               if (xlum_ond(ii,ll,mm) .gt. xlum_max) then
                  xlum_ond(ii,ll,mm)=xlum_ond(ii,ll,mm)*dzeit*dfloat(mm)
               else
                  xlum_ond(ii,ll,mm)=0.d0
               endif
               if (xlum_ond_core(ii,ll,mm) .gt. xlum_max) then
                  xlum_ond_core(ii,ll,mm)=xlum_ond_core(ii,ll,mm)*dzeit
     &                 *dfloat(mm)
               else
C            print *,'xlum_ond_core(ii,ll,mm),xlum_max',
C     &               xlum_ond_core(ii,ll,mm),xlum_max
                  xlum_ond_core(ii,ll,mm)=0.d0
               endif
            enddo
         enddo
      enddo
      print *,'bis,xlum_ond(1,1,1),xlum_ond_core(1,1,1)=',
     &         xlum_ond(1,1,1),xlum_ond_core(1,1,1)
      print *,'bis,xlum_ond(nfreq,lmax,lmax),
     &         xlum_ond_core(nfreq,lmax,lmax)',
     &         xlum_ond(nfreq,llmax,llmax),
     &         xlum_ond_core(nfreq,llmax,llmax)
 
      return
      end
      
C ---------------------------------------------------------------------
      subroutine ond_int(mtu,npasr,ncouche,ll,mm,omegaloc,d_om)
C ---------------------------------------------------------------------
C Calcul de l'epaisseur optique pour une onde avec une frequence
C  qui varie avec la profondeur
C
C d_om est deja la differentielle de la valeur moyenne
C i.e. 0.5*(om_o(i)+om_o(i+1))-om_o(mtu)
C
C OND_INT est appelee par
C
C Auteur : S.Talon
C
C Derniere version : 11 novembre 2000
C ---------------------------------------------------------------------

      implicit none
      include 'evolpar.star'
c      include 'evolcom.transpondesexcit'
      include 'evolcom.igw'

      integer mtu,npasr,ncouche,ll,mm
      double precision omegaloc,d_om(nsh)

      double precision xintret(nsh)
      double precision ampret(nsh)
      double precision t1,t2

      double precision sigmaret,dintret,xN2,ll1
      double precision amp_crit,amp_1,x_depot,abslogamp_crit
      integer n_allerretour
      integer i,k,k1,kk,un,allerretour
      logical retrograde

      allerretour=0
      n_allerretour=0
      amp_crit=1.d-15
      abslogamp_crit=abs(log10(amp_crit))

      ll1=(dfloat(ll+1)*dfloat(ll))**1.5d0
      do k=1,ncouche
        xintret(k)=0.d0
        depotonde(k)=0.d0
      enddo

C Calcul de l'integrale tau (cf. 19, Zahn Talon Matias 1997)
C ----------------------------------------------------------
      k1=mtu
      ampret(k1)=1.d0
      un=1
      retrograde = .true.

C Calcul de l'amplitude de l'onde retrograde
C ------------------------------------------
      do while (retrograde)
        k=k1
        k1=k+un
        sigmaret=omegaloc+dfloat(mm)*d_om(k)
        if (sigmaret.lt.1.d-16) sigmaret=1.d-16

C Modifs CC (10/05/04)
C        xN2=dsqrt(xN_o(k)/(xN_o(k)-sigmaret**2))
        xN2=dsqrt(abs(xN_o(k)/(xN_o(k)-sigmaret**2)))
        dintret=ll1*fint(k)*xN2/sigmaret**4
c dintret < 0
        xintret(k1)=xintret(k)+dintret
c        print *,k,k1,dintret,fint(k)
c        print *,'xintret=',xintret(k1),'dintret=',dintret
C A-t-on vraiment besoin de ce test?
C En fait, quand xint devient trop grand, on met simplement l'amplitude a 0
C Verifier en execution
C         if (xintret(k1).lt.1.d-27) then 
c         if (dabs(xintret(k1)).lt.1.d-20) then 
        if (xintret(k1).lt.1.d-20) then 
C Test CC (21/05/04)
C             xintret(k1)=1.d-27
           xintret(k1)=1.d-20
         endif
        
         if (xintret(k1) .gt. abslogamp_crit) then
c        if (abs(xintret(k1)) .gt. abslogamp_crit) then
C Test sigmaret (24/11/06)
C          ampret(k1)=1.d-16
CC (18/09/09)
            ampret(k1)=0.d0
            retrograde = .false.
         else
            ampret(k1)=dexp(-xintret(k1))
         endif
C	kk=min(k,k1)
C        depotondes(kk)=depotondes(kk)-ampret(k)+ampret(k1)
         depotonde(k1)=depotonde(k1)-ampret(k)+ampret(k1)
         sigmaret=omegaloc+dfloat(mm)*d_om(k+un+un)
         
         if ((un.eq.1 .and. sigmaret.gt.brunt_o(k+un+un)) 
     &       .or. k1.eq.npasr.or.rray_o(k).lt.1.d-3) then
            un=-1
            allerretour=allerretour+1
            
         else if (k1.eq.mtu) then
            un=+1
            allerretour=allerretour+1
            retrograde = .false.
         endif
      enddo

C      print *,'Dans ond_int, apres 1ere boucle'

C Si l'onde n'est pas dissipee apres 1 aller-retour, on calcule d'aller-retour le nombre necessaire,
C et le depot associe
C      print *,'allerretour',allerretour
C      if (allerretour.eq.2) then

      if (allerretour.eq.2.and.ampret(k1).lt.0.9999d0) then
        amp_1=ampret(k1)
        if (dlog10(amp_1).ne.0.d0) then 
            n_allerretour=dint(dlog10(amp_crit)/dlog10(amp_1))
c            n_allerretour=int(log10(amp_1)/log10(amp_crit))
        endif
C Test CC (18/09/09) -->
C         print *,'allerretour,n_allerretour=',allerretour,n_allerretour
C <--
        if (n_allerretour.ne.0) then
           x_depot=1.d0
           do i=1,n_allerretour-1
              x_depot=x_depot+amp_1**i
           enddo
           do i=mtu,npasr
              depotonde(i)=depotonde(i)*x_depot
           enddo
        endif
      endif

      return
      end
      
C ---------------------------------------------------------------------
      subroutine ond_int_core(mtu,npasr,ncouche,ll,mm,omegaloc,d_om)
C ---------------------------------------------------------------------
C Calcul de l'epaisseur optique pour une onde avec une frequence
C  qui varie avec la profondeur
C
C d_om est deja la differentielle de la valeur moyenne
C i.e. 0.5*(om_o(i)+om_o(i+1))-om_o(mtu)
C
C OND_INT est appelee par
C
C Auteur : S.Talon
C
C Derniere version : 10 octobre 2007
C ---------------------------------------------------------------------

      implicit none
      include 'evolpar.star'
c      include 'evolcom.transpondesexcit'
      include 'evolcom.igw'

      integer mtu,npasr,ncouche,ll,mm
      double precision omegaloc,d_om(nsh)

      double precision xintret(nsh)
      double precision ampret(nsh)

      double precision sigmaret,dintret,xN2,ll1
      double precision amp_crit,amp_1,x_depot,abslogamp_crit
      integer n_allerretour
      integer i,k,k1,kk,un,allerretour
      logical retrograde

      allerretour=0
      amp_crit=1.d-15
      abslogamp_crit=abs(log10(amp_crit))

      ll1=(dfloat(ll+1)*dfloat(ll))**1.5d0
      do k=1,ncouche
        xintret(k)=0.d0
        depotonde(k)=0.d0
      enddo

C Calcul de l'integrale tau (cf. 19, Zahn Talon Matias 1997)
C ----------------------------------------------------------
      k1=npasr
      ampret(k1)=1.d0
      un=-1
      retrograde = .true.

C Calcul de l'amplitude de l'onde retrograde
C ------------------------------------------
      do while (retrograde)
        k=k1
        k1=k+un

        sigmaret=omegaloc+dfloat(mm)*d_om(k)
        if (sigmaret.lt.1.d-16) sigmaret=1.d-16

C Modifs CC (10/05/04)
C        xN2=dsqrt(xN_o(k)/(xN_o(k)-sigmaret**2))
        xN2=dsqrt(abs(xN_o(k)/(xN_o(k)-sigmaret**2)))
        dintret=ll1*fint(k)*xN2/sigmaret**4

        xintret(k1)=xintret(k)+dintret
C As-t-on vraiment besoin de ce test?
C En fait, quand xint devient trop grand, on met simplement l'amplitude a 0
C Verifier en execution
C         if (xintret(k1).lt.1.d-27) then 
         if (xintret(k1).lt.1.d-20) then 
C Test CC (21/05/04)
C             xintret(k1)=1.d-27
             xintret(k1)=1.d-20
C             xintret(k1)=1.d-20
         endif
        if (xintret(k1) .gt. abslogamp_crit) then
C Test sigmaret (24/11/06)
C          ampret(k1)=1.d-16
          ampret(k1)=0.d0
          retrograde = .false.
        else
          ampret(k1)=exp(-xintret(k1))
        endif
C	kk=min(k,k1)
C        depotondes(kk)=depotondes(kk)-ampret(k)+ampret(k1)
        depotonde(k1)=depotonde(k1)-ampret(k)+ampret(k1)
        sigmaret=omegaloc+dfloat(mm)*d_om(k+un+un)

        if ((un.eq.-1 .and. sigmaret.gt.brunt_o(k+un+un)) 
     &       .or. k1.eq.mtu) then
          un=1
          allerretour=allerretour+1
        else if (k1.eq.npasr) then
          un=-1
          allerretour=allerretour+1
          retrograde = .false.
        endif
      enddo

C Si l'onde n'est pas dissipee apres 1 aller-retour, on calcule d'aller-retour le nombre necessaire,
C et le depot associe
C      if (allerretour.eq.2) then
      if (allerretour.eq.2.and.ampret(k1).lt.0.9999d0) then
        amp_1=ampret(k1)
        if (log10(amp_1).ne.0.d0) then 
C            n_allerretour=int(log10(amp_crit)/log10(amp_1))
            n_allerretour=int(log10(amp_1)/log10(amp_crit))
        endif
        x_depot=1.d0
        do i=1,n_allerretour-1
          x_depot=x_depot+amp_1**i
        enddo
        do i=mtu,npasr
          depotonde(i)=depotonde(i)*x_depot
        enddo
      endif

      return
      end

C ---------------------------------------------------------------------
      subroutine ond_omega(iterondes,mtu,npasr,ncouche,dt,
     &                     omega_S)
C ---------------------------------------------------------------------
C Calcul du depot de moment par les ondes pour un profil de rotation donne
C
C OND_OMEGA est appelee par TRANSPORT_ONDES
C
C Auteur : S.Talon
C
C Derniere version : 10 mai 2010
C ---------------------------------------------------------------------

      implicit none
      include 'evolpar.star'
      include 'evolcom.diff'
c      include 'evolcom.transpondesexcit'
      include 'evolcom.igw'

      integer mtu,npasr,ncouche,i,k,ii,ll,mm,iterondes
      double precision dt,d_om(nsh)
      double precision omega_S
C Test CC (25/11/04)
      double precision sumdepottot
C Modif CC (28/10/09) 
C Passage de depottot_surf et depottot_core dans evolcom.transpondesexcit
C      double precision depottot_surf(nsh),depottot_core(nsh)

      print *,'entree dans OND_OMEGA'

C      if (iterondes.eq.1) then 
          do k=1,ncouche
            depottot(k)=0.d0
            depottot_surf(k)=0.d0
            depottot_core(k)=0.d0
            d_om(k)=0.d0
          enddo
C      endif

      do i=mtu,npasr
        d_om(i)=0.5d0*(om_o(i)+om_o(i+1))-om_o(mtu)
        d_om(i)=d_om(i)*omega_S
      enddo

      print *,'nfreq=',nfreq,' lmax=',llmax,' On va entrer dans ond_int'
C      if(iterondes.eq.1) then 
      if (igwsurfrot) then 
          do ii=1,nfreq
             do ll=1,llmax
                do mm=1,ll
                   if (xlum_ond(ii,ll,mm).gt.0.d0) then
                      call ond_int(mtu,npasr,ncouche,ll,mm,freqonde(ii)
     &                     ,d_om)
C               ----
C This is the physical angular momentum luminosity transported by waves
                      do k=mtu,ncouche-1
                         depottot_surf(k)=depottot_surf(k)
     &                        +depotonde(k)*xlum_ond(ii,ll,mm)
c                    if (depotonde(k).ne.0.d0) print *,'depot=',
c     &                   depotonde(k)
C                 print *,"dans ond_omega, k,depottot_surf(k)",
C     &           k,depottot_surf(k)
c                    print *,depottot_surf(k)
                      enddo
                   endif
                enddo
             enddo
          enddo
C       do k=mtu,npasr-1
C                 print *,"dans ond_omega, k,depottot_surf(k)",
C     &           k,depottot_surf(k)
C       enddo
C      endif
      endif
C      print *,'Ecriture intermediaire 1 dans ond_omega'
C
      print *,"sortie d'ond_int"


C      if(iterondes.eq.1) then 
      if (igwcorerot) then
          do ii=1,nfreq
          do ll=1,llmax
          do mm=1,ll
             if (xlum_ond_core(ii,ll,mm).gt.0.d0) then
                call ond_int_core(mtu,npasr,ncouche,ll,mm,freqonde(ii),
     &                            d_om)
C               ----
C This is the physical angular momentum luminosity transported by waves
                 do k=mtu,ncouche-1
                    depottot_core(k)=depottot_core(k)
     &                         +depotonde(k)*xlum_ond_core(ii,ll,mm)
C         print *,"dans ond_omega, k,depottot_core(k)",
C     &            k,depottot_core(k)
                 enddo
             endif
          enddo
          enddo
          enddo
C       do k=mtu,npasr-1
C                 print *,"dans ond_omega, k,depottot_core(k)",
C     &           k,depottot_core(k)
C       enddo
C      endif
      endif
C      print *,'Ecriture intermediaire 2 dans ond_omega'
C
C depottot(mtu-1): contient le moment depose dans la ZC de surface
C depottot(npasr): contient le moment depose dans la ZC centrale

      depottot(mtu-1) = 0.d0
      do k=mtu,npasr-1
         depottot(mtu-1)=depottot(mtu-1)-depottot_surf(k)
      enddo

      depottot(npasr) = 0.d0
      do k=mtu,npasr-1
         depottot(npasr)=depottot(npasr)-depottot_core(k)
      enddo

      do k=mtu,npasr-1
         depottot(k) = depottot_surf(k)+depottot_core(k)
C         print *,"dans ond_omega, k,depottot_surf(k),depottot_core(k)",
C     &            k,depottot_surf(k),depottot_core(k)
      enddo
C      stop

C doit valoir 0
      sumdepottot=0.d0
      do k=mtu-1,npasr
         sumdepottot = sumdepottot+depottot(k)
      enddo
C
      print *,' Sortie de ond_omega'

      print *,'ond_omega'
      print *,'depot ZC surf = ',depottot(mtu-1)
      print *,'depot ZC core = ',depottot(npasr)
      print *,'  sumdepottot = ',sumdepottot
C
      return
      end
      
! ----------------------------------------------------------------------
      SUBROUTINE CALCULFLUX(ncouche,mtu,npasr,freq_min,freq_max)
! ----------------------------------------------------------------------
C Ancien program rotation de ST
C Appele par transp_ondesexcit.f ou transp_ondesexcit_20.f
C
C Adaptation pour STAREVOL : C.Charbonnel (15/12/06)
C Contient modifs F.Pantillon (19/06/07)
C Modifs ST & CC pour excitation depuis coeur convectif (10/10/07)
C Derniere version 25 janvier 2010
! ----------------------------------------------------------------------

C      include 'evoldif.inc'
C      include 'd_common'

      implicit none

      include 'evolpar.star'
      include 'evolcom.cons'
      include 'evolcom.teq'
      include 'evolcom.transp'
c      include 'evolcom.transpondesexcit'
      include 'evolcom.igw'
C Modif CC ondes (19/10/07)
      include 'evolcom.diff'

      double precision xmoment_tot
      double precision freq_min,freq_max
      double precision momentum(nsh)
      double precision tot,npi,tamp
c      double precision amplitudemoyenne,amplitudemoyenne_core

      integer i,j
      integer ncouche,mtu,npasr
C Test CC excitation (17/04/07)
      integer minlog
      data minlog/-99./

      data npi/6.2831853d-6/

c      common/amp/amplitudemoyenne(-2:2,nfreq,llmax),
c     &           amplitudemoyenne_core(-2:2,nfreq,llmax)
 
! ----------------------------------------------------------------------

! --------------------------------------------
! Il y a un decalage dans la version inversee des variables
! --------------------------------------------
      do i=1,ncouche 
        xdm_o(i-1) = 0.d0
      enddo

      do i=2,ncouche
c         print *,dm_o(i)
         xdm_o(i-1) = dm_o(i)
      enddo

! -------------------------------
! Calcul de quantites physiques
! -------------------------------

      if (igwrot) then
          xmoment_tot = 0.d0
          do i=1,ncouche-1
C Modif F.Pantillon
C        xmoment_tot = xmoment_tot+rraym2_o(i)*xdm_o(i)
             xmoment_tot = xmoment_tot+2.d0/3.d0*rraym2_o(i)*xdm_o(i)
          enddo
         print *,"CALCULFLUX,  Moment d'inertie ",xmoment_tot
         
         do j=1,ncouche
            momentum(j)=0.d0
         enddo

         do j=1,ncouche-1
C Modif F.Pantillon
C        momentum(j)=0.5d0*(om_o(j)+om_o(j+1))*omega_S*rraym2_o(j)
           momentum(j)=0.5d0*(om_o(j)+om_o(j+1))*omega_S*
     &                 rraym2_o(j)*2.d0/3.d0
         enddo

         tot=0.d0
         do i=1,ncouche-1
           tot=tot+momentum(i)*xdm_o(i)
         enddo
         print 1000,tot
 1000 format ('CALCULFLUX, MOMENT CINETIQUE INIT. : ',1x,e16.9)
      endif

! ------------
! Excitation
! ------------
C      print *,'CALCULFLUX, call Excitation'
      if (igwsurfrot) call excitation(mtu)
!     ----
C Modif CC (30/10/09)
      if (igwcorerot.and.(npasr.ne.ncouche)) 
     &     call excitation_core(npasr,ncouche)
!     ----

      if (igwsurfrot) call damping(mtu,npasr,ncouche,xmoment_tot)
!     ----
C Modif CC (30/10/09)
      if (igwcorerot.and.(npasr.ne.ncouche)) 
     &     call damping_core(mtu,npasr,ncouche,xmoment_tot)
!     ----


! Spectre 0.03 sous la ZC
! -----------------------

C Test CC Excitation (11/04/07)
c      open ( 2,file='ampm_ml', status='unknown',form='formatted')
c      open ( 3,file='ampm_m1', status='unknown',form='formatted')
c      open ( 4,file='ampm_1', status='unknown',form='formatted')
c      open ( 5,file='ampm_l', status='unknown',form='formatted')
c      open ( 6,file='amp_1', status='unknown',form='formatted')
c      open ( 7,file='amp_l', status='unknown',form='formatted')
c      open ( 8,file='flux_1', status='unknown',form='formatted')
c      open ( 9,file='flux_l', status='unknown',form='formatted')

      freq_min = freqrefinert(1)/npi
      freq_max = freqrefinert(nfreq)/npi

c      write (2,*) llmax,nfreq,freq_min,freq_max
c      write (3,*) llmax,nfreq,freq_min,freq_max
c      write (4,*) llmax,nfreq,freq_min,freq_max
c      write (5,*) llmax,nfreq,freq_min,freq_max
c      write (6,*) llmax,nfreq,freq_min,freq_max
c      write (7,*) llmax,nfreq,freq_min,freq_max

c      write (8,*) llmax,nfreq,freq_min,freq_max
c      write (9,*) llmax,nfreq,freq_min,freq_max

      do i=1,nfreq
        amplitudemoyenne(-1,i,1)=amplitudemoyenne(-2,i,1)
        amplitudemoyenne( 1,i,1)=amplitudemoyenne( 2,i,1)

      do j=1,llmax
c        write (2,*) amplitudemoyenne(-2,i,j)
c        write (3,*) amplitudemoyenne(-1,i,j)
c        write (4,*) amplitudemoyenne( 1,i,j)
c        write (5,*) amplitudemoyenne( 2,i,j)

        tamp=(amplitudemoyenne(-1,i,j)-amplitudemoyenne(1,i,j))

        if (tamp.gt.0.d0 .and. fluxtotalafiltrer(i,j).gt.0.d0) then

           xlum_ond(i,j,j) = tamp*fluxtotalafiltrer(i,j)
c          write (6,*) log10(tamp)
c          write (8,*) log10(xlum_ond(i,j,j))
        else
          xlum_ond(i,j,j) = tamp
c          write (6,*) minlog
c          if (xlum_ond(i,j,j).gt.0.d0) then 
c              write (8,*) log10(xlum_ond(i,j,j))
c          else 
c              write (8,*) minlog
c          endif
        endif

        tamp=(amplitudemoyenne(-2,i,j)-amplitudemoyenne(2,i,j))
c        if (fluxtotalafiltrer(i,j).eq.0.d0)!stop 'flux?'
c        tamp = dabs(tamp)

        if (tamp.gt.0.d0 .and. fluxtotalafiltrer(i,j).gt.0.d0) then
           xlum_ond(i,j,1) = tamp*fluxtotalafiltrer(i,j)
c           write (7,*) tamp
c           write (9,*) log10(xlum_ond(i,j,1))
        else
          xlum_ond(i,j,1) = tamp
c          write (7,*) minlog
c          if (xlum_ond(i,j,1).gt.0.d0) then 
c              write (9,*) log10(xlum_ond(i,j,1))
c          else 
c              write (9,*) minlog
c          endif
        endif
      enddo

      enddo

c      close(2)
c      close(3)
c      close(4)
c      close(5)
c      close(6)
c      close(7)
c      close(8)
c      close(9)
      
! Spectre 0.03 au-dessus de la ZC centrale
! ----------------------------------------

c      open ( 2,file='ampm_ml_core', status='unknown',form='formatted')
c      open ( 3,file='ampm_m1_core', status='unknown',form='formatted')
c      open ( 4,file='ampm_1_core', status='unknown',form='formatted')
c      open ( 5,file='ampm_l_core', status='unknown',form='formatted')
c      open ( 6,file='amp_1_core', status='unknown',form='formatted')
c      open ( 7,file='amp_l_core', status='unknown',form='formatted')
c      open ( 8,file='flux_1_core', status='unknown',form='formatted')
c      open ( 9,file='flux_l_core', status='unknown',form='formatted')

      freq_min = freqrefinert(1)/npi
      freq_max = freqrefinert(nfreq)/npi

c      write (2,*) llmax,nfreq,freq_min,freq_max
c      write (3,*) llmax,nfreq,freq_min,freq_max
c      write (4,*) llmax,nfreq,freq_min,freq_max
c      write (5,*) llmax,nfreq,freq_min,freq_max
c      write (6,*) llmax,nfreq,freq_min,freq_max
c      write (7,*) llmax,nfreq,freq_min,freq_max

c      write (8,*) llmax,nfreq,freq_min,freq_max
c      write (9,*) llmax,nfreq,freq_min,freq_max

      do i=1,nfreq
        amplitudemoyenne_core(-1,i,1)=amplitudemoyenne_core(-2,i,1)
        amplitudemoyenne_core( 1,i,1)=amplitudemoyenne_core( 2,i,1)

        do j=1,llmax
c        write (2,*) amplitudemoyenne_core(-2,i,j)
c        write (3,*) amplitudemoyenne_core(-1,i,j)
c        write (4,*) amplitudemoyenne_core( 1,i,j)
c        write (5,*) amplitudemoyenne_core( 2,i,j)

           tamp=(amplitudemoyenne_core(-1,i,j)
     &          -amplitudemoyenne_core(1,i,j))
           if (tamp.gt.0.d0 .and. fluxtotalafiltrer_core(i,j).gt.0.d0) 
     &          then
              xlum_ond_core(i,j,j) = tamp*fluxtotalafiltrer_core(i,j)
c          write (6,*) log10(tamp)
c          write (8,*) log10(xlum_ond_core(i,j,j))
           else
              xlum_ond_core(i,j,j) = tamp
c          write (6,*) minlog
c          if (xlum_ond_core(i,j,j).gt.0.d0) then
c             write (8,*) log10(xlum_ond_core(i,j,j))
c          else
c             write (8,*) minlog
c          endif
           endif

           tamp=(amplitudemoyenne_core(-2,i,j)
     &          -amplitudemoyenne_core(2,i,j))
           if (tamp.gt.0.d0 .and. fluxtotalafiltrer_core(i,j).gt.0.d0) 
     &          then
              xlum_ond_core(i,j,1) = tamp*fluxtotalafiltrer_core(i,j)
c          write (7,*) tamp
c          write (9,*) log10(xlum_ond_core(i,j,1))
           else
              xlum_ond_core(i,j,1) = tamp
c          write (7,*) minlog
c          if (xlum_ond_core(i,j,1).gt.0.d0) then
c              write (9,*) log10(xlum_ond_core(i,j,1))
c          else 
c              write (9,*) minlog
c          endif
           endif
        enddo

      enddo

c      close(2)
c      close(3)
c      close(4)
c      close(5)
c      close(6)
c      close(7)
c      close(8)
c      close(9)

      return
      end

c ---------------------------------------------------------------------
      subroutine excitation(mtu)
C ---------------------------------------------------------------------
C Calcul de l'excitation des ondes par la turbulence "de volume"
C Tient compte que d'un spectre en frequence et en l
C Spectre GMK
C
C Appelee par calculflux.f et calculflux_chimsurf.f
C
C Adapte par S.Talon du programme de P. Kumar
C Adapte pour STAREVOL par C.Charbonnel
C
C Derniere version : 3 avril 2007
C Contient modifs F.Pantillon (19/06/07)
C ---------------------------------------------------------------------
C      include 'evoldif.inc'
C      include 'd_common'

      implicit none
      include 'evolpar.star'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.transp'
c      include 'evolcom.transpondesexcit'
      include 'evolcom.igw'

      double precision kr,L_H,eta,zeta,w1,dw,Edot,termold,rL2,sum,
     &  old_kr,tem,tau_l,tau_h,h_max,term0,term,h_hor
      double precision e_freq,e_tot,nu_min,d_nu,nu_rot,wrot,ratio
      double precision rc,Hp_c,Nc,vc,rhoc,omegac,lc,flux_c
      double precision teff

      integer mtu
      integer ll,ii,i

      write(*,*) 'Entree dans excitation'

      nu_min = 0.3d0                !(microHz)
      d_nu = 0.1d0                  !(microHz)
      
      do ii=1,nfreq
      do ll=1,llmax
C        xlum_ond(ii,ll)=0.d0
        fluxtotalafiltrer(ii,ll)=0.d0
      enddo
      enddo

C Calcul des quantites liees a la ZC
C Ces calculs ont seulement valeur de diagnostic 
C et ne sont pas necessaires pour le calcul de l'excitation

      rc = rtot*rray_o(mtu)
      print *,"Excitation, rc,rtot,rray_o(mtu)",rc,rtot,rray_o(mtu)
      Hp_c = hp_o(mtu)
      ii = mtu
      do while ( (rc-0.2d0*Hp_c) .lt. rray_o(ii)*rtot )
        ii = ii+1
      enddo
      print *,'xNt_o(ii)=',xNt_o(ii)
      Nc = dsqrt(xNt_o(ii))

      print *,"Excitation, Nc = ",Nc
      ii=mtu
      do while ( (rc+0.5d0*Hp_c) .gt. rray_o(ii)*rtot ) 
        ii = ii-1
      enddo
      vc = vconv_o(ii)
      rhoc = rdens_o(ii)
      omegac = 2.d0*pi*vc/Hp_c
      lc     = 2.d0*pi*rc/Hp_c
      flux_c = 0.d0
      do ll = ii,mtu
        flux_c = flux_c + rdens_o(ll)*vconv_o(ll)**3
      enddo
      flux_c = flux_c/(mtu-ii) ! Calcul du flux moyen

C      print *," 0.5 Hp au dessus de la ZC"
C      print *,"     v_c = ",vc
C      print 1400,flux_c
C 1400 format ('  flux_c = ',1x,E16.9)
C      print *,"     l_c = ",lc
C      print *," omega_c = ",omegac
C      print *,"    nu_c = ",omegac/6.2831853d-6
C      print *,"    Hp_c = ",Hp_c
C      print *,"  Hp_c/R = ",Hp_c/rtot
C      print *,""

C      print *,"     N_c = ",Nc
C      print *," log N^2 = ",log10(Nc*Nc)
C      print *,"       m = ",rmass_o(ii)
C      print *,""

c      ii=mtu
c      do while ( (rc-Hp_c) .lt. rray_o(ii)*rtot ) 
c        ii=ii+1
c      enddo

C      print *," Hp sous la ZC "
C      print *,"      ii = ",ii 
C      print 1300,xKm_o(ii)
C 1300 format ('   K_T_c =',1x,E16.9)

C      print *,"    M_zc = ",rmass_o(mtu)
C      print *,"log(M_zc/M) = ",log10(1.d0-rmass_o(mtu))
CC      teff=dsqrt(dsqrt(xluml(1)/4./3.14159/5.67e-5/rtot**2))
CC      print *,"       Teff = ",teff 

C Fin des calculs diagnostics

      eta =pi/2.0d0
      zeta=1.8d0
C zeta: is the aspect ratio of eddies i.e. horizontal_size/radial_size
C free parameter chosen to fit the p-mode energy

      w1 = nu_min*6.2831853d-6
      dw = d_nu*6.2831853d-6
      print 10,nu_min,(nu_min+(nfreq-1)*d_nu)
 10   format(' nu va de ',f5.2,' a ',f5.2,' microHz')

      e_tot = 0.d0
 
      do ii=1,nfreq
        e_freq = 0.d0
	ratio=1.d0 ! Date de l'epoque du "mauvais" Coriolis

        do ll=1,llmax
        Edot = 0.d0
        termold = 0.d0
        rL2 = dfloat(ll*(ll+1))
        sum = 0.d0
        old_kr = 0.d0

c        do i=mtu,1,-1
c           write (994,994),i,xNt_o(i),rL2,rray_o(i),w1,cs_o(i),
c     &          vconv_o(i),alphac,hp_o(i),eta,rdens_o(i),h_max,xNt_o(i)
c 994       format (i4,30(1x,1pe17.10))
c        enddo
c        stop
        do 100 i=mtu,1,-1
          if (xNt_o(i) .gt. 0.d0) goto 100
c calculate the radial wavenumber and its integral 
          tem = rL2/(rray_o(i)*rtot)**2 - w1*w1/cs_o(i)**2
          if (tem .le. 0.d0) goto 100
	  if (vconv_o(i) .le. 0.d0) goto 100
          kr = dsqrt((abs(-xNt_o(i)+w1*w1)*tem))/w1
          sum = sum + 0.5d0*rtot*(rray_o(i)-rray_o(i+1))*(kr+old_kr)
c          write (992,*),i,rtot*(rray_o(i)-rray_o(i+1)),kr+old_kr
          old_kr = kr

c calculate resonant eddy size etc. 
c L_H: vertical size of energy bearing eddies.
          L_H = alphac*hp_o(i)
          tau_L = L_H/vconv_o(i)
          tau_h = tau_L

          if (w1*tau_L .gt. eta) tau_h = eta/w1
          h_max = L_H*(tau_h/tau_L)**1.5d0
c h_max is the size of an eddy that is resonant with the mode.
          term0 = rdens_o(i)*(h_max*h_max)**2*(h_max/L_H)*vconv_o(i)**3
          term0 = zeta*zeta*term0
	  if (sum.gt.5.d1) goto 101
          term = term0*dsqrt(w1*w1-xNt_o(i))*rL2**1.5d0*
     &         exp(-2.d0*sum)
          term = term/(rray_o(i)*rtot*(w1*rray_o(i)*rtot)**2)
c next correct for the high L cutoff of emission from an eddy i.e.
c multiply the above term by a factor of exp(-h_hor^2*L*(L+1)/2r^2)
c where h_hor is the horizontal size of the eddy (h_hor = h_max*zeta;
c where zeta is the aspect ratio).
          h_hor = h_max
          term = term*exp( -rL2*(h_hor/rtot/rray_o(i))**2/2.0d0)
          Edot = Edot+(term+termold)*0.5d0*rtot*(rray_o(i)-rray_o(i+1))
          ! energy luminosity input rate per unit frequency
          termold = term
 100    continue
 101    continue
c        stop "exict"
C Flux total excite, a filtrer dans DAMPING
C        xlum_ond(ii,ll)=Edot*dw/w1*ratio
C Modif F.Pantillon
C        fluxtotalafiltrer(ii,ll)=Edot*dw/w1*ratio
        fluxtotalafiltrer(ii,ll)=6.d0/(8.d0*pi)*Edot*dw/w1*ratio


C Frequence des ondes dans le referentiel inertiel
C (omega dans routines originales de ST)

        freqrefinert(ii)=w1
        e_freq=e_freq+Edot*dw*(2.d0*ll+1.d0)

        enddo ! Boucle sur les l

C        print 200,w1/6.2831853e-6,ratio,e_freq,xlum_ond(ii,lmax)
C        print 200,w1/6.2831853d-6,ratio,e_freq,fluxtotalafiltrer(ii,lmax)
C 200    format(' nu=',f5.2,' ratio=',f5.2,' energy=',e11.4,
C     &         ' xlum=',e11.4)
        e_tot=e_tot+e_freq
        w1=w1+dw
      enddo   ! Boucle sur les frequences

c      do ii=1,nfreq
c         write (996,997) ii,freqrefinert(ii),fluxtotalafiltrer(ii,llmax)
c 997     format (i4,30(1x,1pe17.10))
c      enddo

c      stop "excit"
      write(*,*) 'Sortie de excitation'
      return
      end

C ---------------------------------------------------------------------
      subroutine damping(mtu,npasr,ncouche,xmoment_tot)
C ---------------------------------------------------------------------
C Adapte pour STAREVOL par C.Charbonnel & S.Talon
C Appele par calculflux.f
C Derniere version : 31 octobre 2009
C ---------------------------------------------------------------------

      implicit none
 
      include 'evolpar.star'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.grad'
      include 'evolcom.teq'
      include 'evolcom.transp'
c      include 'evolcom.transpondesexcit'
      include 'evolcom.igw'
      include 'evolcom.therm'
      include 'evolcom.var'
      include 'evolcom.eng'
C Modifs CC ondes (19/10/07)
      include 'evolcom.diff'

      double precision e_total
c      double precision amplitudemoyenne,amplitudemoyenne_core
      
      double precision tau,kr,ll1,e_tot,e_freq,xN2,domega,xmoment_tot
      double precision omm(nsh)
      double precision r_os,npi
C Test CC (31/10/09)
      double precision r_os_core
      double precision rlim,deltar,amp,ampp,ampm
      double precision ampli(nsh)
      double precision xlumj,xint,dintpro
      data npi/6.2831853d-6/

      integer mtu,npasr,ncouche,i,ii,ll,mm,k,k1,j,jj
      integer iwave

c      common/amp/amplitudemoyenne(-2:2,nfreq,llmax),
c     &           amplitudemoyenne_core(-2:2,nfreq,llmax)
      
      print *,"Entree dans DAMPING"
C
C -------------------------------------------------------------------
C Recherche de la frequence de Brunt-Vaisala 0.5Hp au dessous de la ZC
C Servira de valeur minimale a N
      r_os = rray_o(mtu)-0.5d0*hp_o(mtu)/rtot
      print *,"mtu, rray_o(mtu),hp_o(mtu),rtot",mtu,rray_o(mtu),
     & hp_o(mtu),rtot
      print *,"r_OS=",r_os
      ii = mtu+1
      do while (rray_o(ii).gt.r_os)
        ii = ii+1
      enddo
C      i_os=ii
      print *,'ii=',ii
      print *,"rmtu=",rray_o(mtu)
      print *,"  hp=",hp_o(mtu)/rtot
      print *,"  xN=",xN_o(ii)
      print *,"nu max=",dsqrt(xN_o(ii))/npi," microHz"

C Pres de la limite d'un ZC, N --> 0 et formellement aucune
C onde ne peut se propager. Des que l'on a un peu d'ov ce 
C pb est surmonte car N augmente tres rapidement.
C Ces lignes ont ete rajoutees pour eviter une situation ou
C les ondes ne pourraient demarrer.
C Pas de pb car on ne s'interesse qu'aux basses frequences 
C et qu'en pratique a la couche suivante on a deja N > sigma.
C On laisse ces lignes dans la version evolutive pour eviter 
C les cas pathologiques potentiels ou l'EC va migrer.
C
C --> Test ecriture frequence BV avant et apres la condition sur r_os
C CC (18/09/09)
C
C      open (unit = 176,file = 'testFreqBV',status = 'unknown')
C      write (176, 1760)
C      do i=1,ncouche
C         write (176,1762) i,rray_o(i),hp_o(i),xN_o(i),brunt_o(i)
C      enddo
C
C <--
      do i=mtu,ii-1
        xN_o(i)=xN_o(ii)
        brunt_o(i)=brunt_o(ii)
      enddo

C --> Test ecriture frequence BV avant et apres la condition sur r_os
C CC (18/09/09)
C
C      write (176, 1761)
C      do i=1,ncouche
C         write (176,1762) i,rray_o(i),hp_o(i),xN_o(i),brunt_o(i)
C      enddo
C   
C 1760 format (//,5x,'Freq. BV avant le test')
C 1761 format (//,5x,'Freq. BV apres le test')
C 1762 format (2x,i5,4(2x,1pe17.10))
C      close (176)
C      stop
C
CC ----------------------------------------------------------------------------
CC Test CC (31/10/09) -->
CC ----------------------------------------------------------------------------
CCC Recherche de la frequence de Brunt-Vaisala 0.5Hp au-dessus du CC
CCC Servira de valeur minimale a N
CC      r_os_core = rray_o(npasr-1)+0.25d0*hp_o(npasr-1)/rtot
CC      r_os_core = rray_o(npasr-1)+hp_o(npasr-1)/rtot
C      r_os_core = rray_o(npasr-1)+0.4d0*hp_o(npasr-1)/rtot
C      ii = npasr-2
C      do while (rray_o(ii).lt.r_os_core)
C        ii = ii-1
C      enddo
CC
C      print *,"ii=",ii
C      print *,"r_OS_core=",r_os_core
C      print *,"rnpasr-1=",rray_o(npasr-1)
C      print *,"  hp=",hp_o(npasr-1)/rtot
C      print *,"  xN=",xN_o(ii)
C      print *,"nu max=",dsqrt(xN_o(ii))/npi," microHz"
CC
C      do i=ii+1,npasr
C        xN_o(i)=xN_o(ii)
C        brunt_o(i)=brunt_o(ii)
C      enddo
CC
CC Test CC (31/10/09) <-- 
CC
C -------------------------------------------------------------------
C Recherche de la couche a laquelle on va calculer le filtre
C Ici, j'ai mis 0.05R, verifier si c'est adequat
C ATTENTION : CECI POSE UN PB DANS LES PHASES AVANCEES QUAND
C LA FREQUENCE DE B.V. AUGMENTE BEAUCOUP VERS LE COEUR
C A VERIFIER DANS CES PHASES !!!!!
C Modif CC (5/4/7)
      deltar=0.05d0
C      deltar=0.03d0
      r_os=rray_o(mtu)-deltar
      ii=mtu+1
      do while (rray_o(ii).gt.r_os)
         ii=ii+1
      enddo
      iwave=ii

C -------------------------------------------------------------------
C Elimination des ondes non-propagatives
C A garder ??? ST - oct 07
c       do i = 1,nmod
c         print *,"xN_o",i,xN_o(i)
c      enddo
c      stop "xN_o"
      e_tot=0.d0
      do ii=1,nfreq

        e_freq=0.d0
        do ll=1,llmax
c           print *,"damping",ii,ll,xN_o(mtu+1),freqrefinert(ii)**2
          tau=0.d0
          ll1=(dfloat(ll)+1.d0)*ll
          i=mtu+1
c          print *,'comp xN_o',i,xN_o(i),freqrefinert(ii)**2
          if (xN_o(i) .le. freqrefinert(ii)**2) then
C           L'onde ne peut se propager
            fluxtotalafiltrer(ii,ll)=0.d0
C            print *,"A: elimination de nu=",freqrefinert(ii),"l=",ll
          else
            do while (tau.lt.0.5d0.and.xN_o(i).gt.freqrefinert(ii)**2
     &                 .and.i.lt.npasr.and.i.lt.ncouche-2)
              xN2=dsqrt(abs(xN_o(i)/(xN_o(i)-freqrefinert(ii)**2)))
              tau=tau+ll1**1.5d0*fint(i)*xN2/freqrefinert(ii)**4
              i=i+1
            enddo
            kr=dsqrt((xN_o(i-1)/freqrefinert(ii)**2-1.d0)*ll1)
     &         /(rray_o(i)*rtot)
            domega=rray_o(mtu)-rray_o(i)
            if (domega.lt.(1.d0/(kr*rtot)).or.i.eq.(mtu+2)) then
               fluxtotalafiltrer(ii,ll)=0.d0
c               print *,"C: elimination de nu=",freqrefinert(ii),"l=",ll
            endif
          endif
          
          e_freq=e_freq+fluxtotalafiltrer(ii,ll)*freqrefinert(ii)
     &           *(2.d0*ll+1.d0)

        enddo ! Boucle sur les l

        e_tot=e_tot+e_freq
      enddo   ! Boucle sur les frequences

      print 10,e_tot     ! Total energy in gravity waves
 10   format("Total energy flux in waves =",e14.7)
      print 20,e_tot/xmoment_tot
 20   format("Relative energy flux       =",e14.7)
      e_total=e_tot

CC Flux apres avoir elimine les ondes non-propagatives
C      open ( 2,file='flux3', status='unknown',form='formatted')
C      write (2,*) lmax,nfreq,freqrefinert(1)/npi,freqrefinert(nfreq)/npi
C      do ii=1,nfreq
C      do ll=1,lmax
C        if (fluxtotalafiltrer(ii,ll).eq.0.d0) then
C          write (2,*) 0.d0
C        else
C          write (2,*) log10(fluxtotalafiltrer(ii,ll)*ll*(ll+1.d0))
C        endif
C      enddo
C      enddo
C      close(2)

      if (e_tot.eq.0.d0) then
        print *,"Aucune onde ne se propage"
C Modif CC (5/12/07)
C        call exit(0)
      go to 200
      endif

C -------------------------------------------------------------------
C Ici, on veut calculer le flux filtre
C Pour l'instant on calcule m = 1 et m = l et on interpole

      if (iwave .gt. npasr) print *,"ATTENTION PB iwave.gt.npasr"

      do i = 1,ncouche
         omm(i) = 0.d0
      enddo

      do i = 1,ncouche
        omm(i) = 0.5d0*(om_o(i)+om_o(i+1))-om_o(mtu)
        omm(i) = omm(i)*omega_S
      enddo

C -------------------------------------------------------------------
C sortie_ondes <--> 174 <--> sortie_ondes <--> $FILEXw
c      open (unit = 174,file = 'sortie_ondes',status = 'unknown')
c      write (174,1740)
c      write(174,*) r(nmod),gmr(nmod),nmod,omega_S
c      do i=1,ncouche
c         write(174,1704) rray_o(i),xpress_o(i),rdens_o(i),
c     &      grav_o(i),cs_o(i),vconv_o(i),rtemp_o(i),xcp_o(i),
c     &      khi_o(i),xluml_o(i),rmass_o(i),dlnmu_o(i),hp_o(i),
c     &      dm_o(i),rkonv_o(i),om_o(i)
C     &      dm_o(i),xnum_o(i),rkonv_o(i),om_o(i)
c      enddo

c      close (174)
c 1740 format (//,5x,'VARIABLES NEEDED FOR CALCULATIONS OF THE
c     &  INTERNAL GRAVITY WAVES :',//)
c 1704 format (17(1x,1pe17.10))

C ---------------------------------------------------------------
C Modifs ST 18 juin 2007 -->
C Ici, on calcule le depot local de moment par la somme de toutes
C les ondes - on ne considere pas la rotation differentielle
C Ca va servir au calcul de Dondes_chim
      do i=1,ncouche
        lumwave(i)=0.d0
        ampli(i) = 0.d0 
      enddo
      do ii=1,nfreq
      do ll=1,llmax
      do mm=1,ll
        xlumj=fluxtotalafiltrer(ii,ll)*mm
	xint=0.d0
	ampli(mtu)=1.d0
	k=mtu
	k1=mtu+1
	do while (k.lt.npasr)
	  k=k1
          k1=k+1
C Modif CC (11/07/07)
C          xN2=dsqrt(abs(xN_o(i)/(xN_o(i)-freqrefinert(ii)**2)))
          xN2=dsqrt(abs(xN_o(k)/(xN_o(k)-freqrefinert(ii)**2)))
          dintpro=ll1*fint(k)*xN2/freqrefinert(ii)**4
          xint=xint+dintpro
	  ampli(k)=exp(-xint)
	  if (xint.gt.100) then
	    go to 100
	  endif
        enddo

 100	do j=mtu,k-1
          lumwave(j)=lumwave(j)+(ampli(j)-ampli(j+1))*xlumj
C        print *,"lumwave(j) ",lumwave(j),xlumj,ampli(j+1),ampli(j)
	enddo
      enddo
      enddo
      enddo

      if (.not.igwrot) return

C ----------------------------------------------------------------------------
C Ici, on va calculer le flux de moment en tenant compte du profil de rotation
      e_tot=0.d0
      do ii=1,nfreq

        e_freq=0.d0
        do ll=1,llmax
          mm=1
          call ondintexcit(mtu,iwave,ll,mm,freqrefinert(ii),omm,amp,
     &         ampp,ampm)
	  amplitudemoyenne( 1,ii,ll)=ampp
	  amplitudemoyenne(-1,ii,ll)=ampm
          mm=ll
	  call ondintexcit(mtu,iwave,ll,mm,freqrefinert(ii),omm,amp,
     &         ampp,ampm)
	  amplitudemoyenne( 2,ii,ll)=ampp
	  amplitudemoyenne(-2,ii,ll)=ampm

        enddo ! Boucle sur les l

      enddo   ! Boucle sur les frequences

 200  return
      end

C ---------------------------------------------------------------------
      subroutine damping_core(mtu,npasr,ncouche,xmoment_tot)
C ---------------------------------------------------------------------
C Adapte pour STAREVOL par C.Charbonnel & S.Talon
C Appele par calculflux.f
C Derniere version : 16 septembre 2009
C ---------------------------------------------------------------------

      implicit none
 
      include 'evolpar.star'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.grad'
      include 'evolcom.teq'
      include 'evolcom.transp'
c      include 'evolcom.transpondesexcit'
      include 'evolcom.igw'
      include 'evolcom.therm'
      include 'evolcom.var'
      include 'evolcom.eng'

      double precision e_total
c      double precision amplitudemoyenne,amplitudemoyenne_core
      
      double precision tau,kr,ll1,e_tot,e_freq,xN2,domega,xmoment_tot
      double precision omm(nsh)
      double precision r_os,npi
      double precision rlim,deltar,amp,ampp,ampm
      double precision ampli(nsh)
      double precision xlumj,xint,dintpro
      data npi/6.2831853d-6/

      integer mtu,npasr,ncouche,i,ii,ll,mm,k,k1,j
      integer iwave

c      common/amp/amplitudemoyenne(-2:2,nfreq,llmax),
c     &           amplitudemoyenne_core(-2:2,nfreq,llmax)
      
      print *,"Entree dans DAMPING_CORE"

C ----------------------------------------------------------------------------
CC Recherche de la frequence de Brunt-Vaisala 0.5Hp au-dessus du CC
CC Servira de valeur minimale a N
CC Expression sous l'EC:      r_os = rray_o(mtu)-0.5d0*hp_o(mtu)/rtot
      print *,"mtu,npasr,ncouche",mtu,npasr,ncouche
      print *,"rray_o(npasr-1),hp_o(npasr-1),rtot",rray_o(npasr-1),
     &        hp_o(npasr-1),rtot
      r_os = rray_o(npasr-1)+0.25d0*hp_o(npasr-1)/rtot
C      ii = npasr-1
      ii = npasr-2
      do while (rray_o(ii).lt.r_os)
        ii = ii-1
      enddo
C
      print *,"ii=",ii
      print *,"r_OS=",r_os
      print *,"rnpasr-1=",rray_o(npasr-1)
      print *,"  hp=",hp_o(npasr-1)/rtot
      print *,"  xN=",xN_o(ii)
      print *,"nu max=",dsqrt(xN_o(ii))/npi," microHz"
C
C Pres de la limite d'un ZC, N --> 0 et formellement aucune
C onde ne peut se propager. Des que l'on a un peu d'ov ce 
C pb est surmonte car N augmente tres rapidement.
C Ces lignes ont ete rajoutees pour eviter une situation ou
C les ondes ne pourraient demarrer.
C Pas de pb car on ne s'interesse qu'aux basses frequences 
C et qu'en pratique a la couche suivante on a deja N > sigma.
C On laisse ces lignes dans la version evolutive pour eviter 
C les cas pathologiques potentiels ou l'EC va migrer.
C
C --> Test ecriture frequence BV avant et apres la condition sur r_os
C CC (18/09/09)
C
C      open (unit = 175,file = 'testFreqBV',status = 'unknown')
C      write (175, 1750)
C      do i=1,ncouche
C         write (175,1752) i,rray_o(i),hp_o(i),xN_o(i),brunt_o(i)
C      enddo
C
C <--
   
      do i=ii+1,npasr
        xN_o(i)=xN_o(ii)
        brunt_o(i)=brunt_o(ii)
      enddo

C --> Test ecriture frequence BV avant et apres la condition sur r_os
C CC (18/09/09)
C
C      write (175, 1751)
C      do i=1,ncouche
C         write (175,1752) i,rray_o(i),hp_o(i),xN_o(i),brunt_o(i)
C      enddo
C   
C 1750 format (//,5x,'Freq. BV avant le test')
C 1751 format (//,5x,'Freq. BV apres le test')
C 1752 format (2x,i5,4(2x,1pe17.10))
C      close (175)
C      stop
C <--

C ----------------------------------------------------------------------------
C Recherche de la couche a laquelle on va calculer le filtre
C Ici, j'ai mis 0.03R, verifier si c'est adequat
C ATTENTION : CECI POSE UN PB DANS LES PHASES AVANCEES QUAND
C LA FREQUENCE DE B.V. AUGMENTE BEAUCOUP VERS LE COEUR
C A VERIFIER DANS CES PHASES !!!!!
C Modif CC (5/4/7)
      deltar=0.05d0
C Modif CC (12/2/8)
C      r_os=rray_o(npasr)-deltar
      r_os=rray_o(npasr)+deltar
      ii=npasr-1
      do while (rray_o(ii).lt.r_os)
        ii=ii-1
      enddo
      iwave=ii

      print *,"iwave=ii=",iwave

C      stop 

C ----------------------------------------------------------------------------
C Elimination des ondes non-propagatives
C A garder ??? ST - oct 07
      e_tot=0.d0
      do ii=1,nfreq

        e_freq=0.d0
        do ll=1,llmax
          tau=0.d0
          ll1=(dfloat(ll)+1.d0)*ll
          i=npasr-1
          if (xN_o(i) .le. freqrefinert(ii)**2) then
C           L'onde ne peut se propager
            fluxtotalafiltrer_core(ii,ll)=0.d0
C            print *,"A: elimination de nu=",freqrefinert(ii),"l=",ll
          else
            do while (tau.lt.0.5d0.and.xN_o(i).gt.freqrefinert(ii)**2
     &                 .and.i.gt.mtu)
              xN2=dsqrt(abs(xN_o(i)/(xN_o(i)-freqrefinert(ii)**2)))
              tau=tau+ll1**1.5d0*fint(i)*xN2/freqrefinert(ii)**4
              i=i-1
            enddo
            kr=dsqrt((xN_o(i-1)/freqrefinert(ii)**2-1.d0)*ll1)
     &         /(rray_o(i)*rtot)
            domega=rray_o(i)-rray_o(npasr)
            if (domega.lt.(1.d0/(kr*rtot)).or.i.eq.(npasr-2)) then
              fluxtotalafiltrer_core(ii,ll)=0.d0
C              print *,"C: elimination de nu=",freqrefinert(ii),"l=",ll
            endif
          endif
          
          e_freq=e_freq+fluxtotalafiltrer_core(ii,ll)*freqrefinert(ii)
     &           *(2.d0*ll+1.d0)

        enddo ! Boucle sur les l

        e_tot=e_tot+e_freq
      enddo   ! Boucle sur les frequences

      print 10,e_tot     ! Total energy in gravity waves
 10   format("Total energy flux in waves =",e14.7)
      print 20,e_tot/xmoment_tot
 20   format("Relative energy flux       =",e14.7)
      e_total=e_tot

CC Flux apres avoir elimine les ondes non-propagatives
C      open ( 2,file='flux3', status='unknown',form='formatted')
C      write (2,*) lmax,nfreq,freqrefinert(1)/npi,freqrefinert(nfreq)/npi
C      do ii=1,nfreq
C      do ll=1,lmax
C        if (fluxtotalafiltrer_core(ii,ll).eq.0.d0) then
C          write (2,*) 0.d0
C        else
C          write (2,*) log10(fluxtotalafiltrer_core(ii,ll)*ll*(ll+1.d0))
C        endif
C      enddo
C      enddo
C      close(2)

      if (e_tot.eq.0.d0) then
        print *,"Aucune onde ne se propage"
C Modif CC (5/12/07)
C        call exit(0)
        stop
      go to 200
      endif

C ----------------------------------------------------------------------------
C Ici, on veut calculer le flux filtre
C Pour l'instant on calcule m = 1 et m = l et on interpole

      if (iwave .lt. mtu) print *,"ATTENTION PB iwave.lt.mtu"

      do i = 1,ncouche
         omm(i) = 0.d0
      enddo

C Verifier le signe qu'il faut mettre pour omm
      do i = 1,ncouche
        omm(i) = om_o(npasr)-0.5d0*(om_o(i)+om_o(i+1))
        omm(i) = omm(i)*omega_S
      enddo

C ----------------------------------------------------------------------------
C Modifs ST 18 juin 2007 -->
C Ici, on calcule le depot local de moment par la somme de toutes
C les ondes - on ne considere pas la rotation differentielle
C Ca va servir au calcul de Dondes_chim
      do i=1,nsh
        lumwave_core(i) = 0.d0
        ampli(i) = 0.d0
      enddo
      do ii=1,nfreq
      do ll=1,llmax
      do mm=1,ll
        xlumj=fluxtotalafiltrer_core(ii,ll)*mm
	xint=0.d0
	ampli(npasr)=1.d0
	k=npasr
	k1=npasr-1
	do while (k.gt.mtu)
	  k=k1
          k1=k-1
C Modif CC (11/07/07)
C          xN2=dsqrt(abs(xN_o(i)/(xN_o(i)-freqrefinert(ii)**2)))
          xN2=dsqrt(abs(xN_o(k)/(xN_o(k)-freqrefinert(ii)**2)))
          dintpro=ll1*fint(k)*xN2/freqrefinert(ii)**4
          xint=xint+dintpro
	  ampli(k)=exp(-xint)
	  if (xint.gt.100) then
	    go to 100
	  endif
        enddo
C on a mis l'inverse d'en haut mais on n'est pas certaines de tout ca....
C verifier les limites et tout le bazar
 100	do j=k+1,npasr
          lumwave_core(j)=lumwave_core(j)+(ampli(j)-ampli(j+1))*xlumj
	enddo
      enddo
      enddo
      enddo

C ----------------------------------------------------------------------------
C Ici, on va calculer le flux de moment en tenant compte du profil de rotation
      e_tot=0.d0
      do ii=1,nfreq

        e_freq=0.d0
        do ll=1,llmax
          mm=1
	  call ondintexcit_core(npasr,ncouche,iwave,ll,mm,
     &                          freqrefinert(ii),omm,amp,ampp,ampm)
	  amplitudemoyenne_core( 1,ii,ll)=ampp
	  amplitudemoyenne_core(-1,ii,ll)=ampm
          mm=ll
	  call ondintexcit_core(npasr,ncouche,iwave,ll,mm,
     &                          freqrefinert(ii),omm,amp,ampp,ampm)
	  amplitudemoyenne_core( 2,ii,ll)=ampp
	  amplitudemoyenne_core(-2,ii,ll)=ampm

        enddo ! Boucle sur les l

      enddo   ! Boucle sur les frequences

 200  return
      end


C ---------------------------------------------------------------------
      subroutine ondintexcit(mtu,iwave,ll,mm,omegaloc,omm,amp,ampp,ampm)
C ---------------------------------------------------------------------
C Cette routine calcule le flux net de moment cinetique pour une seule
C  onde en tenant compte du profil de rotation.
C
C iwave: numero de la couche ou on calcule le flux filtre
C amp:   amplitude nette, a multiplier avec le flux total
C ampp:  amplitude de l'onde prograde
C ampm:  amplitude de l'onde retrograde
C om:    differentielle de la rotation moyenne i.e. 0.5*(om(i)+om(i+1))-om(mtu)
C
C Correspond a la routine ond_int de S.Talon
C Auteur : S.Talon
C Adaptation pour STAREVOL : C.Charbonnel
C Derniere version : 3 octobre 2007
C ---------------------------------------------------------------------

      implicit none

      include 'evolpar.star'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.grad'
      include 'evolcom.teq'
      include 'evolcom.transp'
c      include 'evolcom.transpondesexcit'
      include 'evolcom.igw'
      include 'evolcom.therm'
      include 'evolcom.var'

      double precision amp_crit,amp_1,n_depot,abslogamp_crit
      double precision omegaloc,omm(nsh),amp,ampp,ampm
      double precision xintpro(nsh),xintret(nsh)
      double precision amppro(nsh),ampret(nsh)
      double precision sigmapro,sigmaret,dintpro,dintret,xN2,ll1

      integer mtu,iwave,ncouche,ll,mm
      integer i,k,k1,un,allerretour
      integer n_allerretour


      amp_crit=1.d-15
      abslogamp_crit=-log(amp_crit)

      ll1=(dfloat(ll+1)*dfloat(ll))**1.5d0

      do k=1,mtu
        xintpro(k)=0.d0
        xintret(k)=0.d0
      enddo

C Calcul de l'integrale tau (cf. 19, Zahn Talon Matias 1997)
C et des amplitudes relatives
C ----------------------------------------------------------
      k=mtu
      k1=mtu
      amppro(k1)=1.d0
      ampret(k1)=1.d0
      un=1

C Calcul de l'amplitude de l'onde retrograde
C ------------------------------------------
      do while (k.le.iwave)
        k=k1
        k1=k+un

        sigmaret=omegaloc+dfloat(mm)*omm(k)
C      write(*,*)'sigmaret= ',sigmaret
        if (sigmaret.lt.0.d0) sigmaret=1.d-15

C Il faut faire le test suivant:
C Normalement, cette expression ne peut etre negative, en fait si elle l'est
C ca signifie que l'onde doit etre reflechie.
C MAIS il ne faut pas laisser xN_o prendre une valeur trop faible,
C i.e. une valeur inferieure a 5.e-4

C Modif CC ondes (3/4/7)
C        xN2=dsqrt(xN_o(k)/(xN_o(k)-sigmaret**2))
        xN2=dsqrt(abs(xN_o(k)/(xN_o(k)-sigmaret**2)))
C      write(*,*)' xN_o(k),  sigmaret**2,  xN2=   ',xN_o(k),
C     &     sigmaret**2,xN2
        
           dintret=ll1*fint(k)*xN2/sigmaret**4
           xintret(k1)=xintret(k)+dintret

        if (xintret(k1).gt.abslogamp_crit) then
          ampret(k1)=0.d0
        else
          ampret(k1)=exp(-xintret(k1))
C      write(*,*)'ampret(k1) =  ',ampret(k1)
        endif
      enddo

      k=mtu
      k1=mtu

C Calcul de l'amplitude de l'onde prograde
C ----------------------------------------
      do while (k.le.iwave)
        k=k1
        k1=k+un

        sigmapro=omegaloc-dfloat(mm)*omm(k)
        if (sigmapro.lt.0.d0) sigmapro=1.d-15

C Modif CC ondes (3/4/7)
C        xN2=dsqrt(xN_o(k)/(xN_o(k)-sigmapro**2))
        xN2=dsqrt(abs(xN_o(k)/(xN_o(k)-sigmapro**2)))
        dintpro=ll1*fint(k)*xN2/sigmapro**4

        xintpro(k1)=xintpro(k)+dintpro
        if (xintpro(k1).gt.abslogamp_crit) then
          amppro(k1)=0.d0
        else
          amppro(k1)=exp(-xintpro(k1))
C      write(*,*)'amppro(k1) = ',amppro(k1)
        endif
      enddo
      
C Ici, il faudra verifier que le signe est conforme a ton programme
C -----------------------------------------------------------------
      amp=ampret(iwave)-amppro(iwave)
      ampp=amppro(iwave)
      ampm=ampret(iwave)

c      print *,'exit',iwave,ll,mm,amp,ampp,ampm

C      write(*,*)'amp, ampp, ampm ',amp,ampp,ampm

      return
      end

C ---------------------------------------------------------------------
      subroutine ondintexcit_core(npasr,ncouche,iwave,ll,mm,omegaloc,
     &                            omm,amp,ampp,ampm)
C ---------------------------------------------------------------------
C Cette routine calcule le flux net de moment cinetique pour une seule
C  onde en tenant compte du profil de rotation.
C
C iwave: numero de la couche ou on calcule le flux filtre
C amp:   amplitude nette, a multiplier avec le flux total
C ampp:  amplitude de l'onde prograde
C ampm:  amplitude de l'onde retrograde
C om:    differentielle de la rotation moyenne i.e. 0.5*(om(i)+om(i+1))-om(mtu)
C
C Correspond a la routine ond_int de S.Talon
C Auteur : S.Talon
C Derniere version : 3 octobre 2007
C ---------------------------------------------------------------------

      implicit none

      include 'evolpar.star'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.grad'
      include 'evolcom.teq'
      include 'evolcom.transp'
c      include 'evolcom.transpondesexcit'
      include 'evolcom.igw'
      include 'evolcom.therm'
      include 'evolcom.var'

      double precision amp_crit,amp_1,n_depot,abslogamp_crit
      double precision omegaloc,omm(nsh),amp,ampp,ampm
      double precision xintpro(nsh),xintret(nsh)
      double precision amppro(nsh),ampret(nsh)
      double precision sigmapro,sigmaret,dintpro,dintret,xN2,ll1

      integer mtu,iwave,ncouche,ll,mm,npasr
      integer i,k,k1,un,allerretour
      integer n_allerretour

      amp_crit=1.d-15
      abslogamp_crit=-log(amp_crit)

      ll1=(dfloat(ll+1)*dfloat(ll))**1.5

      do k=1,mtu
        xintpro(k)=0.d0
        xintret(k)=0.d0
      enddo

C Calcul de l'integrale tau (cf. 19, Zahn Talon Matias 1997)
C et des amplitudes relatives
C ----------------------------------------------------------
      k=npasr
      k1=npasr
      amppro(k1)=1.d0
      ampret(k1)=1.d0
      un=-1

C Calcul de l'amplitude de l'onde retrograde
C ------------------------------------------
      do while (k.ge.iwave)
        k=k1
        k1=k+un

        sigmaret=omegaloc+dfloat(mm)*omm(k)
C      write(*,*)'sigmaret= ',sigmaret
        if (sigmaret.lt.0.d0) sigmaret=1.d-15

C Il faut faire le test suivant:
C Normalement, cette expression ne peut etre negative, en fait si elle l'est
C ca signifie que l'onde doit etre reflechie.
C MAIS il ne faut pas laisser xN_o prendre une valeur trop faible,
C i.e. une valeur inferieure a 5.e-4

C Modif CC ondes (3/4/7)
C        xN2=dsqrt(xN_o(k)/(xN_o(k)-sigmaret**2))
        xN2=dsqrt(abs(xN_o(k)/(xN_o(k)-sigmaret**2)))
C      write(*,*)' xN_o(k),  sigmaret**2,  xN2=   ',xN_o(k),
C     &     sigmaret**2,xN2
        dintret=ll1*fint(k)*xN2/sigmaret**4

        xintret(k1)=xintret(k)+dintret
        if (xintret(k1).gt.abslogamp_crit) then
          ampret(k1)=0.d0
        else
          ampret(k1)=exp(-xintret(k1))
C      write(*,*)'ampret(k1) =  ',ampret(k1)
        endif
      enddo

      k=npasr
      k1=npasr

C Calcul de l'amplitude de l'onde prograde
C ----------------------------------------
      do while (k.ge.iwave)
        k=k1
        k1=k+un

        sigmapro=omegaloc-dfloat(mm)*omm(k)
        if (sigmapro.lt.0.d0) sigmapro=1.d-15

C Modif CC ondes (3/4/7)
C        xN2=dsqrt(xN_o(k)/(xN_o(k)-sigmapro**2))
        xN2=dsqrt(abs(xN_o(k)/(xN_o(k)-sigmapro**2)))
        dintpro=ll1*fint(k)*xN2/sigmapro**4

        xintpro(k1)=xintpro(k)+dintpro
        if (xintpro(k1).gt.abslogamp_crit) then
          amppro(k1)=0.d0
        else
          amppro(k1)=exp(-xintpro(k1))
C      write(*,*)'amppro(k1) = ',amppro(k1)
        endif
      enddo
      
C Ici, il faudra verifier que le signe est conforme a ton programme
C -----------------------------------------------------------------
      amp=ampret(iwave)-amppro(iwave)
      ampp=amppro(iwave)
      ampm=ampret(iwave)

C      write(*,*)'amp, ampp, ampm ',amp,ampp,ampm

      return
      end


c ---------------------------------------------------------------------
      subroutine excitation_core(npasr,ncouche)
C ---------------------------------------------------------------------
C Calcul de l'excitation des ondes par la turbulence "de volume"
C Tient compte que d'un spectre en frequence et en l
C Spectre GMK
C
C Appelee par calculflux.f et calculflux_chimsurf.f
C
C Adapte par S.Talon du programme de P. Kumar
C Adapte pour STAREVOL par C.Charbonnel
C Adaptation pour un coeur convectif (octobre 07)
C
C Derniere version : 3 octobre 2007
C Contient modifs F.Pantillon (19/06/07)
C ---------------------------------------------------------------------
C      include 'evoldif.inc'
C      include 'd_common'

      implicit none
      include 'evolpar.star'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.transp'
c      include 'evolcom.transpondesexcit'
      include 'evolcom.igw'

      double precision kr,L_H,eta,zeta,w1,dw,Edot,termold,rL2,sum,
     &  old_kr,tem,tau_l,tau_h,h_max,term0,term,h_hor
      double precision e_freq,e_tot,nu_min,d_nu,nu_rot,wrot,ratio
      double precision rc,Hp_c,Nc,vc,rhoc,omegac,lc,flux_c
      double precision teff

      integer npasr,ncouche
      integer ll,ii,i

      write(*,*) 'Entree dans excitation_coeur'

      nu_min = 0.3d0                !(microHz)
      d_nu = 0.1d0                  !(microHz)
      
      do ii=1,nfreq
      do ll=1,llmax
        fluxtotalafiltrer_core(ii,ll)=0.d0
      enddo
      enddo

C Calcul des quantites liees a la ZC
C Ces calculs ont seulement valeur de diagnostic 
C et ne sont pas necessaires pour le calcul de l'excitation

      rc = rtot*rray_o(npasr)
      Hp_c = hp_o(npasr)
      ii = npasr
      do while ( (rc+0.2d0*Hp_c) .gt. rray_o(ii)*rtot )
        ii = ii-1
      enddo
      Nc = dsqrt(xNt_o(ii))

      ii=npasr
      do while ( (rc-0.5d0*Hp_c) .lt. rray_o(ii)*rtot ) 
        ii = ii+1
      enddo
      vc = vconv_o(ii)
      rhoc = rdens_o(ii)
      omegac = 2.d0*pi*vc/Hp_c
      lc     = 2.d0*pi*rc/Hp_c
      flux_c = 0.d0
      do ll = npasr,ii
        flux_c = flux_c + rdens_o(ll)*vconv_o(ll)**3
      enddo
      flux_c = flux_c/(ii-npasr) ! Calcul du flux moyen

C      print *," 0.5 Hp au dessus de la ZC"
C      print *,"     v_c = ",vc
C      print 1400,flux_c
C 1400 format ('  flux_c = ',1x,E16.9)
C      print *,"     l_c = ",lc
C      print *," omega_c = ",omegac
C      print *,"    nu_c = ",omegac/6.2831853d-6
C      print *,"    Hp_c = ",Hp_c
C      print *,"  Hp_c/R = ",Hp_c/rtot
C      print *,""

C      print *,"     N_c = ",Nc
C      print *," log N^2 = ",log10(Nc*Nc)
C      print *,"       m = ",rmass_o(ii)
C      print *,""

      ii=npasr
      do while ( (rc+Hp_c) .gt. rray_o(ii)*rtot ) 
        ii=ii-1
      enddo

C      print *," Hp au-dessus de la ZC "
C      print *,"      ii = ",ii 
C      print 1300,xKm_o(ii)
C 1300 format ('   K_T_c =',1x,E16.9)

C      print *,"    M_zc = ",rmass_o(npasr)
C      print *,"log(M_zc/M) = ",log10(rmass_o(npasr))
CC      teff=dsqrt(dsqrt(xluml(1)/4./3.14159/5.67e-5/rtot**2))
CC      print *,"       Teff = ",teff 

C Fin des calculs diagnostics

      eta =pi/2.0d0
      zeta=1.8d0
C zeta: is the aspect ratio of eddies i.e. horizontal_size/radial_size
C free parameter chosen to fit the p-mode energy

      w1 = nu_min*6.2831853d-6
      dw = d_nu*6.2831853d-6
      print 10,nu_min,(nu_min+(nfreq-1)*d_nu)
 10   format(' nu va de ',f5.2,' a ',f5.2,' microHz')

      e_tot = 0.d0
 
      do ii=1,nfreq
        e_freq = 0.d0
	ratio=1.d0 ! Date de l'epoque du "mauvais" Coriolis

        do ll=1,llmax
        Edot = 0.d0
        termold = 0.d0
        rL2 = dfloat(ll*(ll+1))
        sum = 0.d0
        old_kr = 0.d0
        do 100 i=npasr,ncouche
          if (xNt_o(i) .gt. 0.d0) goto 100
c calculate the radial wavenumber and its integral 
          tem = rL2/(rray_o(i)*rtot)**2 - w1*w1/cs_o(i)**2
          if (tem .le. 0.d0) goto 100
	  if (vconv_o(i) .le. 0.d0) goto 100
          kr = dsqrt((-xNt_o(i)+w1*w1)*tem)/w1
          sum = sum + 0.5d0*rtot*(rray_o(i)-rray_o(i+1))*(kr+old_kr)
          old_kr = kr

c calculate resonant eddy size etc. 
c L_H: vertical size of energy bearing eddies.
          L_H = alphac*hp_o(i)
          tau_L = L_H/vconv_o(i)
          tau_h = tau_L
          if (w1*tau_L .gt. eta) tau_h = eta/w1
          h_max = L_H*(tau_h/tau_L)**1.5d0
c h_max is the size of an eddy that is resonant with the mode.
          term0 = rdens_o(i)*(h_max*h_max)**2*(h_max/L_H)*vconv_o(i)**3
          term0 = zeta*zeta*term0
	  if (sum.gt.5.d1) goto 101
          term = term0*dsqrt(w1*w1-xNt_o(i))*rL2**1.5d0*exp(-2.d0*sum)
          term = term/(rray_o(i)*rtot*(w1*rray_o(i)*rtot)**2)
c next correct for the high L cutoff of emission from an eddy i.e.
c multiply the above term by a factor of exp(-h_hor^2*L*(L+1)/2r^2)
c where h_hor is the horizontal size of the eddy (h_hor = h_max*zeta;
c where zeta is the aspect ratio).
          h_hor = h_max
          term = term*exp( -rL2*(h_hor/rtot/rray_o(i))**2/2.0d0)
          Edot = Edot+(term+termold)*0.5d0*rtot*(rray_o(i)-rray_o(i+1))
          ! energy luminosity input rate per unit frequency
          termold = term
 100    continue
 101    continue

C Flux total excite, a filtrer dans DAMPING
C        xlum_ond(ii,ll)=Edot*dw/w1*ratio
C Modif F.Pantillon
C        fluxtotalafiltrer(ii,ll)=Edot*dw/w1*ratio
        fluxtotalafiltrer_core(ii,ll)=6.d0/(8.d0*pi)*Edot*dw/w1*ratio

C Frequence des ondes dans le referentiel inertiel
C (omega dans routines originales de ST)

        freqrefinert(ii)=w1
        e_freq=e_freq+Edot*dw*(2.d0*ll+1.d0)

        enddo ! Boucle sur les l

C        print 200,w1/6.2831853e-6,ratio,e_freq,xlum_ond(ii,lmax)
C        print 200,w1/6.2831853d-6,ratio,e_freq,fluxtotalafiltrer(ii,lmax)
C 200    format(' nu=',f5.2,' ratio=',f5.2,' energy=',e11.4,
C     &         ' xlum=',e11.4)
        e_tot=e_tot+e_freq
        w1=w1+dw
      enddo   ! Boucle sur les frequences

      write(*,*) 'Sortie de excitation_core'
      return
      end

