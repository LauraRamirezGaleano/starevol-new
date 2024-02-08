      subroutine select_table(tau,gstruc,taulim,nmod,t_eff,Tatm
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
      double precision t_eff,Tatm(nmod),Fcatm(nmod),rhoecatm(nmod)
      double precision gstruc,lgstruc,FeH,FeHatm
      double precision rtau,rT,rhotab,Fconvtab
      double precision Ttabeff,gtabeff,Ztabeff
      double precision nouveauT(59),nouveauFc(59),nouveaurho(59)
!      double precision nouveauT(17),nouveauFc(17),nouveaurho(17)
      double precision teffinterp(nmod)
      double precision bip
      double precision tableT(nmod,59),tableFc(nmod,59),
     &     tablerho(nmod,59)
!      double precision tableT(nmod,17),tableFc(nmod,17),
!     &     tablerho(nmod,17)
      double precision tableTZ,tableFcZ,tablerhoZ
      double precision Fc_atm,rho_atm,interpT
      double precision stau(127),sT(127),srhotab(127),sFconvtab(127)
!      double precision stau(120),sT(120),srhotab(120),sFconvtab(120)
      double precision srhotabZ,sFconvtabZ,sTZ
      double precision T_geff(12),rho_geff(12),Fconv_geff(12)
!      double precision T_geff(3),rho_geff(3),Fconv_geff(3)
      double precision T_Zeff(6),rho_Zeff(6),Fconv_Zeff(6)

      double precision taucoupleatm,rayenv

      dimension sum(59,10)
      dimension tau(nsh)

      common /atmospheres/ rtau(127,12,59,6),rT(127,12,59,6),
     &     rhotab(127,12,59,6),Fconvtab(127,12,59,6),Ttabeff(59),
     &     gtabeff(12),Ztabeff(6),filesize,taumin,taumax
      common /atmospheres2/ tableTZ(127,59,12),tablerhoZ(127,59,12),
     &     tableFcZ(127,59,12)
c$$$      common /atmospheres/ rtau(120,3,17,3),rT(120,3,17,3),
c$$$     &     rhotab(120,3,17,3),Fconvtab(120,3,17,3),Ttabeff(17),
c$$$     &     gtabeff(3),Ztabeff(3),filesize,taumin,taumax
c$$$      common /atmospheres2/ tableTZ(120,17,3),tablerhoZ(120,17,3),
c$$$     &     tableFcZ(120,17,3)
      
      common /metal/ FeH
      common /overshoot/ klcore,klenv,klpulse
  

      lgstruc = LOG10(gstruc)
      stau(:) = rtau(:,1,1,1)

      taumax = 100!/mtini
      taumin = 10!/mtini
!      taumin = 30

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
            call splintatm (gtabeff,T_geff,DDT,nb_geff,lgstruc,
     &           sT(k))

            rho_geff(:) = tablerhoZ(k,i,:)
            call splineatm (gtabeff,rho_geff,nb_geff,1.d50,1.d50,DDT,
     &           DT)
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
         call spline (stau,sT,filesize,1.d50,1.d50,DDT,DT) !1.d50 to obtain 'natural' limit conditions : second derivative = 0 at 0 and N.
         call spline (stau,sFconvtab,filesize,1.d50,1.d50,DDT1,DT)
         call spline (stau,srhotab,filesize,1.d50,1.d50,DDT2,DT)
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
      
!!!  3.  Selection of the closest temperature track for tau=[taumin;taumax] and mean it. => t_eff
      Teffinterp = 0.d0
      t_eff = 0.d0
      nb_elmt = 0
      do k=nmod,1,-1
         if (tau(k)>=taumin) exit
      enddo
      rayenv = r(novlim(klenv,4))-r(novlim(klenv,3))
!!      print *,'rayenv=',rayenv/r(nmod)
      if (rayenv.gt.r(nmod)*1.d-2) then
         do while (tau(k)<=taumax)
            nouveauT(:) = tableT(k,:)
            call splineatm (nouveauT,Ttabeff,nb_Teff,1.d50,1.d50,DDT,
     &           DT)
            call splintatm (nouveauT,Ttabeff,DDT,nb_Teff,t(k),
     &           Teffinterp(k))
            t_eff = t_eff + Teffinterp(k)
            k = k-1
            nb_elmt = nb_elmt + 1
         enddo
         t_eff = t_eff/(nb_elmt)
      else
         do k=nmod,1,-1
            if (tau(k)>=taucoupleatm) exit
         enddo
         nouveauT(:) = tableT(k,:)
         call splineatm (nouveauT,Ttabeff,nb_Teff,1.d50,1.d50,DDT,DT)
         call splintatm (nouveauT,Ttabeff,DDT,nb_Teff,t(k),
     &        Teffinterp(k))
         t_eff = Teffinterp(k)
      endif
      write(*,*) "Interpolating at Teff =", t_eff
      
!!!  4.  Creation of the effective T(tau) (and others) track(s) selected to determine q(tau)
      k = nmod
      do while (tau(k)<=taumax)
         nouveauT(:) = tableT(k,:)
         nouveauFc(:) = tableFc(k,:)
         nouveaurho(:) = tablerho(k,:)
         call splineatm (Ttabeff,nouveauT,nb_Teff,1.d50,1.d50,DDT,DT)
         call splintatm (Ttabeff,nouveauT,DDT,nb_Teff,t_eff,Tatm(k))
         call splineatm (Ttabeff,nouveauFc,nb_Teff,1.d50,1.d50,DDT,DT)
         call splintatm (Ttabeff,nouveauFc,DDT,nb_Teff,t_eff,
     &        Fcatm(k))
         call splineatm (Ttabeff,nouveaurho,nb_Teff,1.d50,1.d50,DDT,DT)
         call splintatm (Ttabeff,nouveaurho,DDT,nb_Teff,t_eff,
     &        rhoecatm(k))
         k = k-1
      enddo

      end subroutine select_table
