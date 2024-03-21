      SUBROUTINE mchange

************************************************************************
* Modify the mass distribution due to mass-loss and/or accretion   
* WARNING : never use t it is not defined yet, use lnt !
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.ion'

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
      
      integer error
      integer imenv,ii1,jj,jj1,nn1,nold,mm
      integer ibot,itop,iqmrbot,iqmrtop
      integer mlp1,nmod0,ibind,mlsh
      integer i,ii,ij,j,k,nn,j1
      integer itmax,ishockb,ishockt
      integer nconvsh,ienvsh
      integer icrazy
c      integer nenv,nend

      logical arndt,blocker,limit
      character prescrip*120,zdep*50

      double precision :: mass_loss_deJager
      double precision :: mass_loss_Crowther
      double precision :: mass_loss_vanLoon
      double precision :: mass_loss_Grafener_Hamann
      double precision :: mass_loss_Nugis_Lamers
      double precision :: mass_loss_Vink
      double precision :: mass_loss_Sander_Vink
      double precision :: mass_loss_Wood

      double precision xmm,momlim,dmsmec,rapom,omegacrit
      double precision meini,Leini,omsurfini,dLrad,omegamax
      double precision mgammaedd
      double precision dlogm,dma
      double precision mominnew,mominold
      double precision plog,pcomp,mdlog,vesc,dmsmax,dmsreim
      double precision qmr,vqmr,xtot
      double precision rco,tteff,lleff,totmv,totmd,
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
      double precision ledd
      double precision geffrot,Rp2
      double precision kap,kapm
      double precision disctime
      double precision alfa,logdms1,logdms2,dms1,dms2,dT
      double precision tau1,tau2,taus,gammaeff,vinf,gammarot
     $     ,vesceff,dmsthick
      double precision Dclump,etawind
      double precision acoeff, cbd, logmdotoff, gammaeb,logdmsvms
      
!     TODO (VOJE) : All common blocks should be in include files.
      common /disclocking/ disctime
      common /hydrodyn/ Lmax,tmax,mackmax,enucmax,itmax,ishockb,ishockt
      common /meshconv/ menvconv,dmenvconv,mcore,nconvsh,ienvsh
      common /newmesh/ qmr(nsh),vqmr(nsh),iqmrbot,iqmrtop
      common /metal/ FeH
      common /boreas_Mloss/ Bcrit,fstar,Bequi
      common /geffective/ geffrot(nsh),Rp2(nsh)
      common /opa/ kap(nsh),kapm(nsh)
c     common /eddington/ ledd


      if (abs(mlp0).gt.49) then
         limit = .true.
      else
         limit = .false.
      endif
c$$$      if (mlp0.gt.0) then
c$$$         zscaling = .true.
c$$$      else
c$$$         zscaling = .false.
c$$$      endif
c$$$      zscaling = .false.        ! Modif TD et CC Juillet 2018 (com supr.)

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
c      if (dmaccr.eq.0.d0.and.(nphase.lt.2.or.mlp.eq.0)) then
      if (dmaccr.eq.0.d0.and.((nphase.lt.2.and.mlp.ne.23)
     $    .or.mlp.eq.0)) then
!!         if (mlp.eq.23.and.time/sec.le.disctime) then
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
               
*_________________________________
***   mass-loss rate prescriptions
*---------------------------------

***   Reimers (1975) prescription (MS and RGB phases)
            if (mlp.eq.1.or.mlp.eq.2) then
               dms = min(etapar*3.98d-13*soll*solr/totm,1.d-4)
               print *,'Reimers mass loss',dms,'solar masses / year'
            endif

!     Multiple prescriptions
!     if log(Teff)>3.7
!        de Jager et al. (1988,A&AS,72,259)
!     if log(Teff)<3.7
!        Crowther (2001)
!     Massive stars with RGB phase
            if (mlp.eq.3.or.mlp.eq.4) then
               if (tteff.gt.3.7d0) then
                  dms = mass_loss_deJager(tteff,lleff)
               else
                  dms = mass_loss_Crowther(lleff)
               endif
               if (limit) dms = min(dms,1.d-4)
            endif

***   Vassiliadis & Wood (1993, ApJ 413, 641) prescription (AGB phase)
            if (mlp.eq.5.or.mlp.eq.6.or.mlp.eq.55.or.mlp.eq.56) then
               dms = mass_loss_Wood(solr,totm,leff,mlp)
c..   WARNING : limitation of mass loss rate to 1d-4 Mo/yr
               if (limit) dms = min(dms,-4.d0)
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

!     Multiple prescriptions
!     if log(Teff)>3.9
!        Vink et al. (2001,A&A,369,574)
!     if log(Teff)<3.9
!        de Jager et al. (1988,A&AS,72,259)
!     if H_surf<0.4 and log(Teff)>4.0
!        Nugis & Lamers (2000,A&A,360,227)
!     Massive stars, MS, post-MS, RSG and WR phase
            if (mlp.eq.15.or.mlp.eq.16) then
               if (tteff.le.4.0d0) then
                  !write (nout,200)
                  dms = mass_loss_deJager(tteff,lleff)
               else if (xsp(nmod1,ih1).le.0.4d0.and.tteff.gt.4.d0) then
                  !write (nout,500)
                  !dms = mass_loss_Nugis_Lamers(soll,zkint,xsp,nmod)
                  dms = mass_loss_Sander_Vink(soll,totm,zkint,
     &                 xsp,nmod)
               else
                  dms = mass_loss_Vink(soll,teff,totm,zkint,xsp,nmod1,
     &                 clumpfac)
               endif
            endif

!     Crowther (2001)
!     red-super giants
            if (mlp.eq.19.or.mlp.eq.20) then
               dms = mass_loss_Crowther(lleff)
            endif

!     van Loon et al (2005,A&A,438,273)
!     AGB & red supergiants
            if (mlp.eq.21.or.mlp.eq.22) then
               dms = mass_loss_vanLoon(soll,teff)
            endif

!     Constant mass loss rate (see starevol.par)
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
                  logmdotoff = 0.23d0*log10(zeff/zsol)-2.61
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
            
***   ---------------------------------------------------------------------------------------------------------------
            
            if (mlp1.eq.15.or.mlp1.eq.16) mlp = mlp1
            
***   ---------------------------------------------------------------------------------------------------------------

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
               print *,'mgamma',mgamma,'mgammaedd',mgammaedd
               mrom = m(nmod)*pw34/(pi*r(nmod)**3)
               if (tteff.ge.4.35d0) malpha = 0.52d0
               if (tteff.ge.4.30d0.and.tteff.lt.4.35d0)
     &              malpha = 0.24d0
               if (tteff.ge.4.d0.and.tteff.lt.4.3d0)
     &              malpha = 0.17d0
               if (tteff.ge.3.9d0.and.tteff.lt.4.d0)
     &              malpha = 0.15d0
               mexp = 1.d0/malpha-1.d0
c               if (vsurf.eq.0.d0) then
c                  vsurf = r(nmod)*vomega(nmod)
c               else
c                  vsurf = vsurf*1.d5
c               endif
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
         

***   Metallicity dependence
         zdep = ' '
c         if (zscaling.and.((mlp.ge.11.and.mlp.ne.17.and.mlp.ne.18)
c     &        .or.mlp.le.4.or.mlp.eq.7.or.mlp.eq.8)) then
         if (zscaling.and..not.((mlp.eq.15.or.mlp.eq.16)
     &        .and.tteff.gt.3.9)) then
C           dms = dms*dsqrt(zeff*5.d1)
C     Mokiem et al 2007
C            dms = dms*dsqrt(zeff/zsol)
            dms = dms*(zeff/zsol)**0.8d0
c            dms = dms*3
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
         
         print*,"Mass loss rate (Msun/yr) :",dms
*_____________________________________________________
***   change shells mass distribution due to mass-loss
*-----------------------------------------------------

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

***   if mlp odd : mass removed in shells above m = m(mlsh)
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

c         print *,'omegasurf',vomega(nmod),'omegacrit',dsqrt(g*m(nmod)
c     $        /(1.5d0*rtot)**3), 'rapport',vomega(nmod)/dsqrt(g*m(nmod)
c     $        /(1.5d0*rtot)**3)

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
 200  format (5x,'WARNING : Change mass loss --> de Jager')
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
      
      end SUBROUTINE mchange

****************************************************************************

      double precision function mass_loss_Vink(soll,teff,totm,zeff,xsp,
     &     nmod1,clumpfac)
      include 'evolpar.star'
      include 'evolcom.nuc'
      integer, intent(in) :: nmod1
      double precision, intent(in) :: soll, teff, totm, zeff,
     &     xsp(nsh,nsp), clumpfac
!      double precision :: HeH, qsig
      double precision :: sigmae, gammae, brackro
      double precision :: teffjump1, teffjump2, alfa, ratio
      double precision :: dms1, dms2, logdms1, logdms2
      print*,'Mass loss prescription : '//
     &'Vink et al. (2001,A&A,369,574)'
!     Lamers & Leitherer (1993), ApJ 412, p771, eq. 2
!      HeH = xsp(nmod1,ihe4)/(xsp(nmod1,ihe4)+xsp(nmod1,ih1))
!      if (teff.lt.3.d4) then
!         qsig = 0d0
!      elseif (teff.ge.3.d4.and.teff.lt.3.5d4) then
!         qsig = 0.5d0
!      else
!         qsig = 1.d0
!      endif
!      sigmae = 0.401d0*(1.d0+qsig*HeH)/(1.d0+3.d0*HeH)
      sigmae = 0.325d0
      gammae = 7.66d-5*sigmae*soll/totm
      brackro = -14.94d0+3.1857d0*gammae+0.85d0*
     &     log10(zeff/0.019d0)
!     WARNING : Unused teffjump1 !!!
      teffjump1 = 1.d3*(61.2d0+2.59d0*brackro)
      teffjump2 = 1.d3*(1.d2+6.d0*brackro)
      if (teffjump1.le.teffjump2) then
         write (nout,'("unrealistic stellar parameters")')
         stop 'mchange'
      endif
      if (teff.ge.27500.d0) then
         ratio = 2.6d0
         logdms1 = -6.697d0+2.194d0*log10(soll*1.d-5)
     &        -1.313d0*log10(totm/30.d0)
     &        -1.226d0*log10(ratio*0.5d0)
     &        +9.33d-1*log10(teff/4.d4)
     &        -10.92d0*(log10(teff/4.d4))**2
     &        +0.85d0*log10(zeff/0.019d0)
         mass_loss_Vink = 10.d0**logdms1
      else if (teff.le.22500.d0) then
         if (teff.ge.teffjump2) then
            ratio = 1.3d0
            logdms2 = -6.688d0+2.210d0*log10(soll*1.d-5)
     &           -1.339d0*log10(totm/30.d0)
     &           -1.601d0*log10(ratio*0.5d0)
     &           +1.07d0*log10(teff/2.d4)
     &           +0.85d0*log10(zeff/0.019d0)
            mass_loss_Vink = 10.d0**logdms2
         else
            ratio = 0.7d0
            logdms2 = -5.99d0+2.210d0*log10(soll*1.d-5)
     &           -1.339d0*log10(totm/30.d0)
     &           -1.601d0*log10(ratio*0.5d0)
     &           +1.07d0*log10(teff/2.d4)
     &           +0.85d0*log10(zeff/0.019d0)
            mass_loss_Vink = 10.d0**logdms2
         endif
c..   Smoothing the transitions following the eval_Vink_wind routine from MESA
      else if (teff.gt.22500.d0.and.teff.lt.27500.d0) then   
         alfa = (teff-22500.d0)/(5.d3)
         logdms1 = -6.697d0+2.194d0*log10(soll*1.d-5)
     &        -1.313d0*log10(totm/30.d0)
     &        -1.226d0*log10(2.6d0*0.5d0)
     &        +9.33d-1*log10(teff/4.d4)
     &        -10.92d0*(log10(teff/4.d4))**2
     &        +0.85d0*log10(zeff/0.019d0)
         logdms2 = -6.688d0+2.210d0*log10(soll*1.d-5)
     &        -1.339d0*log10(totm/30.d0)
     &        -1.601d0*log10(1.3d0*0.5d0)
     &        +1.07d0*log10(teff/2.d4)
     &        +0.85d0*log10(zeff/0.019d0)
         dms1 = 10.d0**logdms1
         dms2 = 10.d0**logdms2
         mass_loss_Vink = (1-alfa)*dms2 + alfa*dms1
      endif
c Mass loss rate from Vink et al 2000 modified to account for the
c effects of rotation (factor 0.85)
c     dms = 0.8d0*dms
!     Mass loss rate modified to account for the effect of clumping!
      mass_loss_Vink = mass_loss_Vink*clumpfac
      print*,'Clumping factor :',clumpfac
      print*,'teffjump1 :',teffjump1
      print*,'teffjump2 :',teffjump2
      end function mass_loss_Vink

****************************************************************************
      
      double precision function mass_loss_deJager(tteff,lleff)
      double precision, intent(in) :: tteff, lleff
      integer :: i, j, ii, ii1, jj, jj1, nn, nn1
      double precision :: xmlos, ymlos, logdms
      double precision :: amlos(6,5)
!     parameters for the de Jager mass loss prescription
      data ((amlos(i,j),i = 1,6),j = 1,5)/6.34916d0,3.41678d0,
     &   -1.08683d0,0.13095d0,0.22427d0,0.11968d0,-5.0424d0,0.15629d0,
     &   0.41952d0,-0.09825d0,0.46591d0,0.d0,-0.83426d0,2.96244d0,
     &   -1.37272d0,0.13025d0,0.d0,0.d0,-1.13925d0,0.33659d0,-1.07493d0,
     &   0.d0,0.d0,0.d0,-0.12202d0,0.57576d0,0.d0,0.d0,0.d0,0.d0 /
      print*,'Mass loss prescription : '//
     &'de Jager et al. (1988,A&AS,72,259)'
      logdms = 0.d0
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
               logdms = logdms
     &              - amlos(ii,jj)*cos(dble(ii1)*
     &              acos(xmlos))*cos(dble(jj1)*acos(ymlos))
            endif
         enddo
      enddo
      mass_loss_deJager = 10**logdms
      end function mass_loss_deJager

****************************************************************************

      double precision function mass_loss_Sander_Vink(soll,totm,zeff,
     &     xsp,nmod)
      include 'evolpar.star'
      include 'evolcom.nuc'
      integer, intent(in) :: nmod
      double precision, intent(in) :: soll, totm, zeff, xsp(nsh,nsp)
      double precision :: acoeff, cbd, logmdotoff, gammaeb, logdms
      double precision :: zsol
      print*,'Mass loss prescription : '//
     &'Sander & Vink (2020,MNRAS,499,873)'
!     ToDo : Change that.
      zsol = 0.014
      acoeff = 2.932d0
      gammae = 10.d0**(-4.8125)*(1+xsp(nmod,ih1))*soll/totm
      gammaeb = -0.324d0*log10(zeff/zsol)+0.244
      logmdotoff = 0.23d0*log10(zeff/zsol)-2.61
      cbd = -0.44d0*log10(zeff/zsol)+9.15
      logdms = acoeff*log10(-log10(1-gammae))-log10(2.d0)
     $     *(gammaeb/gammae)**cbd+logmdotoff
      mass_loss_Sander_Vink = 10**logdms
      end function mass_loss_Sander_Vink

****************************************************************************
      
      double precision function mass_loss_Nugis_Lamers(soll,zeff,xsp,
     &     nmod)
      include 'evolpar.star'
      include 'evolcom.nuc'
      integer, intent(in) :: nmod
      double precision, intent(in) :: soll, zeff, xsp(nsh,nsp)
      double precision :: logdms
      print*,'Mass loss prescription : '//
     &'Nugis & Lamers (2000,A&A,360,227)'
      logdms =-11.d0+1.29d0*log10(soll)+
     &     1.73d0*log10(xsp(nmod,ihe4))+0.47d0*log10(zeff)
      mass_loss_Nugis_Lamers = 10**logdms
      end function mass_loss_Nugis_Lamers

****************************************************************************

      double precision function mass_loss_Crowther(lleff)
!     Mass loss rate for RSG based on Sylvester+1998 and van Loon+99 (LMC)
      double precision, intent(in) :: lleff
      double precision :: logdms
      print*,'Mass loss prescription : '//
     &'Crowther (2001,ASSL,264,215)'
      logdms = 1.7d0*lleff-13.83d0
      mass_loss_Crowther = 10.d0**logdms
      end function mass_loss_Crowther

****************************************************************************

      double precision function mass_loss_Grafener_Hamann(lum,m,soll,
     &     teff,xsp,nmod)
!     Grafener & Hamann (2008)
      include 'evolcom.cons'
      include 'evolpar.star'
      include 'evolcom.nuc'
      integer, intent(in) :: nmod
      double precision, intent(in) :: lum(nsh), m(nsh), soll, teff
      double precision, intent(in) :: xsp(nsh,nsp)
      double precision :: msigmae, mgamma, logdms
      msigmae = 0.2d0*(1.d0+xsp(nmod,ih1)+xsp(nmod,ih2))
      mgamma = msigmae*lum(nmod)/(4*pi*g*c*m(nmod))
      logdms = 10.046d0+1.727d0*log10(mgamma-0.326d0)
     &     -3.5d0*log10(teff)+0.42d0*log10(soll)
     &     -0.45d0*xsp(nmod,ih1)
      mass_loss_Grafener_Hamann = 10**logdms
      end function mass_loss_Grafener_Hamann

****************************************************************************
      
      double precision function mass_loss_vanLoon(soll,teff)
      double precision, intent(in) :: soll, teff
      double precision :: logdms
      print*,'Mass loss prescription : '//
     &'van Loon et al (2005,A&A,438,273)'
      logdms = -5.65d0+1.05d0*log10(soll*1.d-4)-6.3d0*
     &     log10(teff/3500.d0)
      mass_loss_vanLoon = 10.d0**logdms
      end function mass_loss_vanLoon

****************************************************************************

      double precision function mass_loss_Wood(solr,totm,leff,mlp)
      integer, intent(in) :: mlp
      double precision, intent(in) :: solr, totm, leff
      double precision :: plog, pcomp, mdlog, vesc, dmsmax 
      print*,'Mass loss prescription : '//
     &'Vassiliadis & Wood (1993,ApJ,413,641)'
      plog = -2.07d0+1.94d0*log10(solr)-0.9d0*log10(totm)
      pcomp = 10.d0**plog
      if (totm.le.2.5d0.and.(mlp.eq.55.or.mlp.eq.56)) then
c..   modification for M > 2.5Mo to delay the onset of the superwind phase
         print*,'Original prescription'
         mdlog = -11.4d0+0.0125d0*(pcomp-1.d2*(totm-2.5d0))
      else
         print*,'No delay of the onset of the super-wind phase'
         mdlog = -11.4d0+0.0125d0*pcomp
      endif
      vesc = -13.5d0+0.056d0*pcomp
      vesc = vesc*1.d5
      vesc = max(3.d5,vesc)
      vesc = min(1.5d6,vesc)
      dmsmax = leff/(c*vesc)
      dmsmax = dmsmax*sec/msun
      mass_loss_Wood = min(10.d0**mdlog,dmsmax)
      end function mass_loss_Wood

****************************************************************************

      
