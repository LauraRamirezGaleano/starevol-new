
************************************************************************
* This part has been developed by:                                     *
*    Corinne Charbonnel                                                *
*    Version 1.0: May 2001                                             *
************************************************************************

      SUBROUTINE DALPHA (omst,omss,omtt,xkt,xMs,xMt,xm,xn,xs,xt,D1,D2,
     &     alpha)

*-----------------------------------------------------------------------
*     DALPHA calcule le coefficient de diffusion de Paquette a
*     partir des integrales de collision.
*     Derniere version : 30 novembre 1996
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
*-----------------------------------------------------------------------

      implicit none

      double precision omst,omss,omtt,xkt,xMs,xMt,xm,xn,xs,xt,D1,D2,
     &     alpha
      double precision AA,BB,CC,EE
      double precision Ps,Pt,Pst,Qs,Qt,Qst,Ss,St
      double precision delta

      double precision Z_st,Z_st1,Z_st2,K_st 

      common /paqthoul/ Z_st,Z_st1,Z_st2,K_st                                                 ! Ajout TD Juin 2018 pour coef. Paquette

      dimension omst(2,3),omss(2,3),omtt(2,3)

      AA = omst(2,2)/(5.d0*omst(1,1))
      BB = (5.d0*omst(1,2)-omst(1,3))/(5.d0*omst(1,1))
      CC = (2.d0*omst(1,2))/(5.d0*omst(1,1))-1.d0
      EE = xkt/(8.d0*xMs*xMt*omst(1,1))

      Ps = 8.d0*xMs*EE*omss(2,2)/(5.d0*xkt)
      Pt = 8.d0*xMt*EE*omtt(2,2)/(5.d0*xkt)
      Pst = 3.d0*(xMs-xMt)**2+4.d0*xMs*xMt*AA
      Qs = Ps*(6.d0*xMt*xMt+5.d0*xMs*xMs-4.d0*xMs*xMs*BB+8.d0*xMs*xMt*
     &     AA)
      Qt = Pt*(6.d0*xMs*xMs+5.d0*xMt*xMt-4.d0*xMt*xMt*BB+8.d0*xMt*xMs*
     &     AA)
      Qst = 3.d0*(xMs-xMt)**2*(5.d0-4.d0*BB)+4.d0*xMs*xMt*AA*
     &     (11.d0-4.d0*BB)+2.d0*Ps*Pt
      Ss = xMs*Ps-xMt*(3.d0*(xMt-xMs)+4.d0*xMs*AA)
      St = xMt*Pt-xMs*(3.d0*(xMs-xMt)+4.d0*xMt*AA)

      D1 = 3.d0*EE/(2.d0*xn*xm)
      delta = 5.d0*CC*CC*(xMs*xMs*Ps*xs*xs+xMt*xMt*Pt*xt*xt+Pst*xs*xt)/
     &     (xs*xs*Qs+xt*xt*Qt+xs*xt*Qst)
      D2 = D1/(1.d0-delta)
      alpha = 5.d0*CC*(xs*Ss-xt*St)/(xs*xs*Qs+xt*xt*Qt+xs*xt*Qst)

*** Coeff pour Thoul      
      Z_st = -CC
      Z_st1 = -2.d0*BB+2.5d0
      Z_st2 = 5.d0*AA
      K_st = xs*xt*xn*xkt/D1

*** Fin Coeff pour Thoul      
c       write(66,*) 'Debut'
c         write(66,'(9(1x,1pe11.4))') CC,xs,
c     &   Ss,xt,St,Qs,Qt,Qst,alpha
c$$$      print *, omst(1,2),omst(1,1)
c$$$      print *, (2.d0*omst(1,2))/(5.d0*omst(1,1))-1.d0
c$$$      print *, 'Zdiff', Z_st, 'Zdiff1', Z_st1, 'Zdiff2', Z_st2,
c$$$     $     'Kdiff', K_st
c$$$      stop
      return
      end



c     SUBROUTINE DIFFMIC (ielem,ndt,ndb,vdiff,diff12,ration,anuc,znuc)
      SUBROUTINE DIFFMIC (ielem,ndt,ndb,vdiff,diff12,anuc,znuc)

*-----------------------------------------------------------------------
*     Calcule la vitesse de diffusion microscopique :
*     Approximation d'un pur potentiel de Coulomb
*     (Chapman et Cowling,1970).
*     DIFFMIC est appelee par DIFFEL ou par ELEMENTS_FINIS
*     Derniere version : 8 aout 1996
*-----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.transp'
      include 'evolcom.var'

      integer ielem,ndt,ndb
      integer l

c     double precision vdiff,diff12,ration,anuc,znuc
      double precision vdiff,diff12,anuc,znuc
      double precision e2
      double precision xma,za,xmb,xmab,rho,temp,zb,xlambdaD,xa12,
     &     xdiff,xvdiff,dlogT,dray,a12
      double precision VD,DT,VDT

      dimension diff12(nsh),vdiff(nsh)
c     dimension ration(nsh,nis,18)
      dimension anuc(nsp),znuc(nsp)

*     Initialisations.

      e2 = 2.307121d-19

      xma = 1.d0
      za = 1.d0
      do l = 1,nmod
         diff12(l) = 0.d0
         vdiff(l) = 0.d0
      enddo
      xmb = anuc(ielem)
      xmab = (xma+xmb)/(xma*xmb)*avn

      do l = ndb+1,ndt
         rho = ro(l)
         temp = t(l)
         zb = znuc(ielem)
         xlambdaD = dsqrt((boltz*temp)/(4.d0*pi*(rho*avn)*e2))
         xa12 = (4.d0*boltz*temp)/(za*e2)*xlambdaD
         xdiff = 3.d0/(8.d0*rho*avn)*dsqrt(boltz*temp*xmab/(2.d0*pi))
         xvdiff = -(xma/avn)/(boltz*temp)*g*m(l)/(r(l)*r(l))
         dlogT = log(t(l+1))-log(t(l-1))
         dray = r(l+1)-r(l-1)
         xa12 = xa12/zb
         a12 = log(1.d0+xa12*xa12)
         diff12(l) = xdiff*(2.d0*boltz*temp/(za*zb*e2))**2/a12

*     Diffusion Thermique

         VD = xvdiff*diff12(l)*(xmb-zb/2.d0-0.5d0)
         DT = 2.65d0*zb*zb+0.805d0*zb*(zb-1.d0)
         VDT = DT*dlogT
         VDT = VDT/dray
         VDT = VDT*diff12(l)
         vdiff(l) = VD+VDT 
      enddo

      return
      end



c     SUBROUTINE DIFFPAQ (ielem,ndt,ndb,vdiffp,diff13,abond,ration,
c    &     ratneh,anuc,znuc)
      SUBROUTINE DIFFPAQ (ielem,ndt,ndb,vdiffp,diff13,abond)!,ratneh)

*-----------------------------------------------------------------------
*     Calcule la vitesse de diffusion microscopique :
*     Expression de Montmerle & Michaud (1976)
*     Approximation de Paquette.
*
*     Derniere version : 6 janvier 2003
*-----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.nuc'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.transp'
      include 'evolcom.var'
      include 'evolcom.ion'

      integer ielem,ndt,ndb
      integer nlimitediff
      integer l
      integer i                 ! ion part TD Oct.2018

      double precision vdiffp,diff13,abond,ration!,ratneh!,anuc,znuc   ! ion part TD Oct.2018
c     double precision vdiffp,diff13,abond,ratneh 
      double precision xvdiffbe,xvdbe,xvthermbe,xdmicbe
      double precision xvdiffhe,xvdhe,xvthermhe,xdmiche,xdthermhe
     $     ,xvdiffhenoc
      double precision xvdiffli,xvdli,xvthermli,xdmicli
      double precision xvdiffc12,xvdc12,xvthermc12,xdmicc12,xvdiffc12noc
      double precision xvdiffn14,xvdn14,xvthermn14,xdmicn14,xvdiffn14noc       ! add by TD Fev.2018
      double precision xvdiffo16,xvdo16,xvthermo16,xdmico16,xvdiffo16noc
      double precision xvdiffne20,xvdne20,xvthermne20,xdmicne20
     $     ,xvdiffne20noc   
      double precision xvdiffna23,xvdna23,xvthermna23,xdmicna23
     $     ,xvdiffna23noc   ! add by TD Fev.2018
c      double precision xdthermc12,xdthermo16
      double precision dlnTdr
c     double precision altse,diffpaqt,ddth,dmic
      double precision altse                                      ! ion part TD Oct.2018
      double precision altei,altpi,altep,diffpaqt,ddth,dmic
      double precision xma,zielem,xmb,as,zs,rho,temp,xvdiffp,zb
      double precision yb,vtot                                     ! ion part TD Oct.2018
      double precision VD,vtherm
      double precision gammat,xKPE,xlambi,alphaT

      common /coefdiffpaqbe/ xvdiffbe(nsh),xvdbe(nsh),
     &     xvthermbe(nsh),xdmicbe(nsh)
      common /coefdiffpaqhe/ xvdiffhe(nsh),xvdhe(nsh),
     &     xvthermhe(nsh),xdmiche(nsh),xdthermhe(nsh),xvdiffhenoc(nsh)
      common /coefdiffpaqli/ xvdiffli(nsh),xvdli(nsh),
     &     xvthermli(nsh),xdmicli(nsh)
      common /coefdiffpaqc12/ xvdiffc12(nsh),xvdc12(nsh),
     &     xvthermc12(nsh),xdmicc12(nsh),xvdiffc12noc(nsh)
      common /coefdiffpaqn14/ xvdiffn14(nsh),xvdn14(nsh),       ! add by TD Fev.2018
     &     xvthermn14(nsh),xdmicn14(nsh),xvdiffn14noc(nsh)
      common /coefdiffpaqo16/ xvdiffo16(nsh),xvdo16(nsh),       
     &     xvthermo16(nsh),xdmico16(nsh),xvdiffo16noc(nsh)
      common /coefdiffpaqne20/ xvdiffne20(nsh),xvdne20(nsh),    ! add by TD Fev.2018
     &     xvthermne20(nsh),xdmicne20(nsh),xvdiffne20noc(nsh)
      common /coefdiffpaqna23/ xvdiffna23(nsh),xvdna23(nsh),    ! add by TD Fev.2018
     &     xvthermna23(nsh),xdmicna23(nsh),xvdiffna23noc(nsh)
      common /dlnTdrray/ dlnTdr(nsh)
c      common /ioni/ zmean(nsh,nsp)                   ! ion part TD Oct.2018
      
c      dimension anuc(nsp),znuc(nsp)
      dimension diff13(nsh),vdiffp(nsh)
c      dimension ration(nsh,nis,18)!,ratneh(nsh),znuc(nsp)
      dimension abond(nsh,nis),ration(nsh,nis,18)!,ratneh(nsh)
c      dimension xdthermc12(nsp),xdthermo16(nsp)

*     Initialisations.
      
      xma = 1.d0
      zielem = znuc(ielem)
      do l = 1,nmod
         diff13(l) = 0.d0
         vdiffp(l) = 0.d0
      enddo
      xmb = anuc(ielem)         ! xmb = at
      as = 1.d0
      zs = 1.d0

      if (ndb.eq.1) nlimitediff = 2
      if (ndb.ne.1) nlimitediff = ndb

c$$$      print *, 'Par Athena !'
c$$$      do i=1,nsp
c$$$         do l = nlimitediff,ndt
c$$$            write(666,*) i,l,t(l)
c$$$     $           ,zmean(l,i)
c$$$         enddo
c$$$      enddo
c      stop

c Traitement de l'helium 4
      if (ielem.eq.ihe4) then 
         do l = nlimitediff,ndt
            rho = ro(l)
            temp = t(l)
      
c         xvdiffp = -xma/(avn*boltz*temp)*g*m(l)/(r(l)*r(l))

*     Calcul avec ionisation partielle
            
c$$$            do i = 1,int(znuc(ielem))
c$$$               if (ration(l,ielem,i).gt.1.d-5) then
c$$$                  zb = dble(i)
c$$$                  yb = abond(l,ielem)*ration(l,ielem,i)
c$$$                  call paq (as,zs,abond(l,2),xmb,zb,yb,rho,temp,dmic
c$$$     $                 ,ddth,altse,altep,diffpaqt,ratneh(l),muizero(l)) ! modif TD mars 2018 : zielem par xmb
c$$$                  VD = xvdiffp*dmic*(xmb-zb/2.d0-0.5d0)
c$$$                  vtherm = dmic*altse*dlnTdr(l)
c$$$                  vtot = vtherm+VD
c$$$                  vdiffp(l) = vdiffp(l)+ration(l,ielem,i)*vtot
c$$$                  diff13(l) = diff13(l)+ration(l,ielem,i)*dmic
c$$$               endif
c$$$  enddo
            !ionisation partielle TD (10/2018)
c$$$            zb = zmean(l,ielem)
c$$$            call paqHe (as,zs,abond(l,2),xmb,zb,abond(l,ielem),rho, ! modif TD mars 2018 : zielem par xmb
c$$$     &           temp,dmic,ddth,altei,altep,altpi,diffpaqt,ratneh(l),
c$$$     &           muizero(l))
c$$$            gammat = abond(l,ielem)
c$$$            xKPE = ((1.d0+0.5d0*(zb+1.d0)*gammat)*(1.d0+xmb*zb*gammat))
c$$$     &           /((1.d0+xmb*gammat)*(1.d0+0.5d0*(zb+1.d0)*zb*gammat))
c$$$            
c$$$            xlambi=1.d0+gammat*(zb-1.d0)*(1.d0-gammat)
c$$$            
c$$$            alphaT=(altpi+zs*altei-zb*altep)/(zs+1.d0)
c$$$
c$$$            xvdiffp = -rho/p(l)*g*m(l)/(r(l)*r(l))
c$$$            VD=dmic*(1.d0+gammat)*(((xmb-1.d0)+xKPE*(xmb-zb))
c$$$     &           /((1.d0+xmb*gammat)*xlambi))*xvdiffp
c$$$               
c$$$            vtherm = dmic*(1.d0+gammat)*alphaT*dlnTdr(l) ! Thermal part
c$$$            vdiffp(l) = vtherm+VD
c$$$            diff13(l) = dmic
c$$$c Pour l'ecriture
c$$$            xvdiffhe(l) = vdiffp(l)
c$$$            xvdhe(l) = VD
c$$$            xvthermhe(l) = vtherm
c$$$            xdmiche(l) = dmic
c$$$            xdthermhe(l)=alphaT               
c$$$         enddo

*     Calcul pour un milieu completement ionise
   
            zb = znuc(ielem)
            call paqHe (as,zs,abond(l,2),xmb,zb,abond(l,ielem),rho,         ! modif TD mars 2018 : zielem par xmb
     &           temp,dmic,ddth,altei,altep,altpi,diffpaqt,ratneh(l),
     &           muizero(l))
c            print *,'par Toutatis !', dlnTdr 
c Vitesse de diffusion de Montmerle & Michaud (1976 APJS 31,489)
c     Expression A.5
            gammat = abond(l,ielem)
            xKPE = ((1.d0+0.5d0*(zb+1.d0)*gammat)*(1.d0+xmb*zb*gammat))
     &           /((1.d0+xmb*gammat)*(1.d0+0.5d0*(zb+1.d0)*zb*gammat))

            xlambi=1.d0+gammat*(zb-1.d0)*(1.d0-gammat)

            alphaT=(altpi+zs*altei-zb*altep)/(zs+1.d0)

            xvdiffp = -rho/p(l)*g*m(l)/(r(l)*r(l))
            VD=dmic*(1.d0+gammat)*(((xmb-1.d0)+xKPE*(xmb-zb))
     &           /((1.d0+xmb*gammat)*xlambi))*xvdiffp

            vtherm = dmic*(1.d0+gammat)*alphaT*dlnTdr(l) ! Thermal part
            vdiffp(l) = vtherm+VD
            diff13(l) = dmic
c Pour l'ecriture
            xvdiffhe(l) = vdiffp(l)
            xvdhe(l) = VD
            xvthermhe(l) = vtherm
            xdmiche(l) = dmic
            xdthermhe(l)=alphaT
c            print *, xdmiche(l)
c            write(551,'(1x,i4,7(1x,1pe11.4))'),l,r(l),xdmiche(l)
         enddo
c         stop
      endif

      if (ielem.ne.ihe4) then 
         do l = nlimitediff,ndt
            rho = ro(l)
            temp = t(l)
c         xvdiffp = -xma/(avn*boltz*temp)*g*m(l)/(r(l)*r(l))

*     Calcul avec ionisation partielle

c$$$            do i = 1,int(znuc(ielem))
c$$$c               if (ration(l,ielem,i).gt.1.d-10) print *,'et pourtant...'
c$$$               if (ration(l,ielem,i).gt.1.d-5) then
c$$$                  zb = dble(i)
c$$$                  yb = abond(l,ielem)*ration(l,ielem,i)
c$$$c                  print *, 'pk tant de haine',yb,abond(l,ielem)
c$$$                  call paq (as,zs,abond(l,2),xmb,zb,yb,rho,temp,dmic
c$$$     $                 ,ddth,altse,altep,diffpaqt,ratneh(l),muizero(l)) ! modif TD mars 2018 : zielem par xmb
c$$$                  VD = xvdiffp*dmic*(xmb-zb/2.d0-0.5d0)
c$$$                  vtherm = dmic*altse*dlnTdr(l)
c$$$                  vtot = vtherm+VD
c$$$                  vdiffp(l) = vdiffp(l)+ration(l,ielem,i)*vtot
c$$$                  diff13(l) = diff13(l)+ration(l,ielem,i)*dmic
c$$$               endif
c$$$  enddo
               !ionisation partielle TD (10/2018)
c$$$               if (ielem.eq.ic12.or.ielem.eq.in14.or.ielem.eq.io16) then
c$$$                  zb = zmean(l,ielem)
c$$$               else
c$$$                  zb = znuc(ielem)
c$$$               endif
c$$$               call paq (as,zs,abond(l,2),xmb,zb,abond(l,ielem),             ! modif TD mars 2018 : zielem par xmb
c$$$     &              rho,temp,dmic,ddth,altpi,altep,diffpaqt,ratneh(l),
c$$$     &              muizero(l))
c$$$                xvdiffp = -rho/p(l)*g*m(l)/(r(l)*r(l))
c$$$                VD = xvdiffp*dmic*(2.d0*xmb-zb-1.d0)
c$$$                vtherm = dmic*(altpi-0.5d0*altep*zb)*dlnTdr(l)
c$$$                vdiffp(l) = vtherm+VD
c$$$                diff13(l) = dmic

*     Calcul pour un milieu completement ionise
            zb = znuc(ielem)
            
c            NE = rhosh/(amu*muesh) ! vient de Thoul, voir si possible adapter à MM
            call paq (as,zs,abond(l,2),xmb,zb,abond(l,ielem),             ! modif TD mars 2018 : zielem par xmb
     &           rho,temp,dmic,ddth,altpi,altep,diffpaqt,ratneh(l),
     &           muizero(l))
c Vitesse de diffusion de Montmerle & Michaud (1976 APJS 31,489)
c Expression A.5
            xvdiffp = -rho/p(l)*g*m(l)/(r(l)*r(l))
            VD = xvdiffp*dmic*(2.d0*xmb-zb-1.d0)
            vtherm = dmic*(altpi-0.5d0*altep*zb)*dlnTdr(l)
c            if (ielem.eq.io16) then
c               write(67,'(1x,i4,1x,9(1x,1pe11.4))') l,dmic,altpi,
c     &              altep, zb, dlnTdr(l),vtherm
c            endif
            vdiffp(l) = vtherm+VD
            diff13(l) = dmic

c Pour ecriture                
            if (ielem.eq.ic12) then
               xvdiffc12(l)=vdiffp(l)
               xvdc12(l)=VD
               xvthermc12(l)=vtherm
               xdmicc12(l)=dmic
c            xdthermc12(l)=(altpi - 0.5d0*altep*zb)
            endif
            if (ielem.eq.in14) then                       ! add by TD Fev.2018
               xvdiffn14(l)=vdiffp(l)
               xvdn14(l)=VD
               xvthermn14(l)=vtherm
               xdmicn14(l)=dmic
c            xdthermn14(l)=(altpi - 0.5d0*altep*zb)
            endif
             if (ielem.eq.io16) then
               xvdiffo16(l)=vdiffp(l)
               xvdo16(l)=VD
               xvthermo16(l)=vtherm
               xdmico16(l)=dmic
c            xdthermo16(l)=(altpi - 0.5d0*altep*zb)
            endif
            if (ielem.eq.ine20) then                        ! add by TD Fev.2018
               xvdiffne20(l)=vdiffp(l)
               xvdne20(l)=VD
               xvthermne20(l)=vtherm
               xdmicne20(l)=dmic
c            xdthermne20(l)=(altpi - 0.5d0*altep*zb)
            endif
            if (ielem.eq.31) then                        ! add by TD Fev.2018
               xvdiffna23(l)=vdiffp(l)
               xvdna23(l)=VD
               xvthermna23(l)=vtherm
               xdmicna23(l)=dmic
c            xdthermna23(l)=(altpi - 0.5d0*altep*zb)
            endif

c Calcul de la force radiative
c         if(fradia) then
c             if (ielem.eq.6.or.ielem.eq.7) then
c                vrad=dmic*Frad(l)/(boltz*temp)
c                vdiffp(l)=vdiffp(l)+vrad
c             endif
c          endif
            if (ielem.eq.ili7) then 
               xvdiffli(l) = vdiffp(l)
               xvdli(l) = VD
               xvthermli(l) = vtherm
               xdmicli(l) = dmic
            endif
            if (ielem.eq.ibe9) then 
               xvdiffbe(l) = vdiffp(l)
               xvdbe(l) = VD
               xvthermbe(l) = vtherm
               xdmicbe(l) = dmic
            endif
         enddo
      endif

*     Ajout CC (Calcul precedent fait uniquement entre ndt et ndb+1,
*     pour eviter pb dans cette version quand il n'y a plus
*     de coeur convectif)
*     ----> Valeurs en ndb = valeurs en ndb+1

      if (ndb.eq.1) then 
         vdiffp(ndb) = vdiffp(ndb+1)
         diff13(ndb) = diff13(ndb+1)
      endif

      return
      end

************************************************************************
      SUBROUTINE DIFFTHOUL (jshell,xsp,anuc,znuc,tsh,rhosh,muesh,abond,
     $     nmod)
*--------------------------------------------------------------------------------------
* This routine inverses the Burgers equations for microscopic diffusion
* Adapted from the diffusion.f routine by Anne Thoul
* See Thoul et al, ApJ 421, p 828 (1994) for details
*
* Original routine retrieved from :
* http://www.astro.ulg.ac.be/orientation/asterosis/article/EleDif/explication.html
*
* Version May 2007 (Starevol V2.92)
*     Ana Palacios
* Maj May-Juin 2018 (Thibaut Dumont)      
*
* From original routine header :
*===============================
*
C The parameter M is the number of species considered.
C
C Fluid 1 is the hydrogen
C Fluid 2 is the helium
C Fluids 3 to M-1 are the heavy elements
C Fluid M is the electrons
C
C The vectors A,Z and X contain the atomic mass numbers, 
C the charges (ionization), and the mass fractions, of the elements.
C NOTE: Since M is the electron fluid, its mass and charge must be
C      A(M)=m_e/m_u
C      Z(M)=-1.
C
C The array CL contains the values of the Coulomb Logarithms.
C The vector AP, AT, and array AX contains the results for the diffusion 
C coefficients.
C The vector C contains the concentrations
C CC is the total concentration: CC=sum(C_s)
C AC is proportional to the mass density: AC=sum(A_s C_s)
C The arrays XX,Y,YY and K are various parameters which appear in 
C Burgers equations.
C The vectors and arrays alpha, nu, gamma, delta, and ga represent
C the "right- and left-hand-sides" of Burgers equations, and later 
C the diffusion coefficients.
C
*-------------------------------------------------------------------------------------
      implicit none

      include 'evolcom.cons'
c      include 'eoscom.cons'    ! Ajout TD
      include 'evolpar.star'
      include 'evolcom.ion'     ! Ajout TD
      
      integer i,j,l,jshell,mspec,n,mmax,nmax,crz,lim,nmod,nmod_td
     $     ,nmod_td2
c$$$      parameter (mspec = 32)
c$$$      parameter (mmax = 32,nmax = 66)     ! avec isotopes jusqu'à Mg
      parameter (mspec = 36)
      parameter (mmax = 36,nmax = 74)      ! avec isotpoes après Mg
      ! H He
c$$$      parameter (mspec = 3)                   ! modif mai 2018
c$$$  parameter (mmax = 3, nmax = 8)
      ! H He C N O
c$$$      parameter (mspec = 6)                   ! modif mai 2018
c$$$      parameter (mmax = 6, nmax = 14)
      ! H He C N O Ne Na Mg Al27 Si28 P31 S32 Cl35
c$$$      parameter (mspec = 14)                   ! modif mai 2018 : suppr. Cl36
c$$$      parameter (mmax = 14, nmax = 30)
      
      integer indx(nmax)
      integer ineutron,ih1,ih2,ihe3,ihe4,ili6,ili7,ibe7,ibe9,ib8,ib11,
     &     ic12,ic13,ic14,in13,in14,in15,io15,io16,io17,io18,if19,ine20,
     &     ina23,img24,ial27,isi28,ip31,is32,icl35,idxspc,idxsps,idxspci
      integer ina24,ina25,ine23,ial26g,ial26m,if18,if20,ina22,img27,
     &     is35,icl36,img25,ine21,ine22,img26,ib10
      integer elemthoul, Heindex
      integer fcs,nphase,novlim                ! first convective shell
      
      double precision anuc,znuc,tsh,rhosh,abond,muesh    ! modif Fev 2019 muesh
      
      double precision xkt,xtemp,sum1
      double precision a,z,xthoul,ap,at,ax,cl,abondthoul
      double precision dc,cc,ac,xx,y,yy,k
      double precision Zs,Zt,ys,yt,ys0,yt0,rho,ratnehk,xmuizero
      double precision rray,rkonv,rkonvL,muizero
      double precision e2,amu,xmh,xn,xne,xns,
     &     xnt,sumnz2,xlambdad,xlambdai,xlambda,gradx
      double precision alpha,nu,gamma,delta,ga
      double precision k0,temp,d,xij,xi,k0thoul,xior
      double precision dlntdr,dlnpdr,dlnxdr,dlncdr
      double precision xsp,masselec,NE,NI,AO,CZ,LAMBDAD,LAMBDA,ZXA,NIZI2 ! Ajout TD + modif Fev.2019
      double precision XE, ZXB
      double precision xvdiffhe,xdmiche,xvdiffhenoc
      double precision xvdiffli
      double precision xvdiffbe
      double precision xvdiffc12,xvdiffc12noc,xvdiffc13,xvdiffc13noc
      double precision xvdiffn14,xvdiffn14noc,xvdiffn15,xvdiffn15noc     
      double precision xvdiffo16,xvdiffo16noc,xvdiffo17,xvdiffo17noc
     $     ,xvdiffo18,xvdiffo18noc
      double precision xvdiffne20,xvdiffne20noc   
      double precision xvdiffna23,xvdiffna23noc ! add by TD Mai.2018
      double precision vdiffthoul, dmicthoul
      double precision Zdiff,Zdiff1,Zdiff2,Kdiff,kappa_st,Z_st,Z_st1
     $     ,Z_st2,K_st                                                         ! pour coef Paq 
      double precision altei,altpi,altep,diffpaqt,ddth,dmic  ! pour coef Paq
      double precision ratio_ion,fr ! test TD Jan.2019
      double precision chronos  ! Ajout chronos pour mesure du tps de calcul
      double precision e_ap,e_at,e_ax,g_ap,g_at,g_ax,g_tot,e_tot,sum2
     $     ,sum3    ! Calcul champ electrique + gravité
      
      parameter (masselec = 9.1093897d-28)        ! Ajout TD
c      common /xspnuc/ xsp(nsh,nsp) ! Ajout TD
      common /paramthoul/ NE,NI,AO,CZ,LAMBDAD,LAMBDA, NIZI2, XE                            ! Ajout TD + modif Fev.2019
      common /dlntdrray/ dlntdr(nsh),dlnpdr(nsh),dlnxdr(nsh,nsp)
      common /dmicthoul/ vdiffthoul(nsh,mspec), dmicthoul(nsh,mspec)  ! modif
      common /specindex/ ineutron,ih1,ih2,ihe3,ihe4,ili6,ili7,ibe7,ibe9,
     &     ib8,ib10,ib11,ic12,ic13,ic14,in13,in14,in15,io15,io16,io17,
     &     io18,
     &     if19,ine20,ine21,ine22,ina23,img24,ial27,isi28,ip31,is32,
     &     icl35,ina24,
     &     ina25,ine23,img26,ial26g,ial26m,if18,if20,ina22,img27,is35,
     &     icl36,img25,idxspc(nprintc),idxsps(nprint),idxspci(9)
      common /thoulindex/ elemthoul(mspec-1)                        ! Ajout TD AP
      common /star/ rray(nsh),rkonv(nsh),rkonvL(nsh),muizero(nsh)
      common /coefdiffpaqhe/ xvdiffhe(nsh),xvdiffhenoc(nsh),xdmiche(nsh)
      common /coefdiffpaqli/ xvdiffli(nsh)
      common /coefdiffpaqbe/ xvdiffbe(nsh)
      common /coefdiffpaqc12/ xvdiffc12(nsh),xvdiffc12noc(nsh)
     $     ,xvdiffc13(nsh),xvdiffc13noc(nsh)
      common /coefdiffpaqn14/ xvdiffn14(nsh),xvdiffn14noc(nsh)
     $     ,xvdiffn15(nsh),xvdiffn15noc(nsh)
      common /coefdiffpaqo16/ xvdiffo16(nsh),xvdiffo16noc(nsh)
     $     ,xvdiffo17(nsh),xvdiffo17noc(nsh),xvdiffo18(nsh)
     $     ,xvdiffo18noc(nsh)
      common /coefdiffpaqne20/ xvdiffne20(nsh),xvdiffne20noc(nsh)
      common /coefdiffpaqna23/ xvdiffna23(nsh),xvdiffna23noc(nsh) ! Ajout TD Mai 2018
      common /paqthoul/ Z_st,Z_st1,Z_st2,K_st                                          ! Ajout TD Juin 2018 pour coef. Paquette       
c      common /ioni/ zmean(nsh,nsp)                   ! ion part TD Oct.2018
      common /shelltype/ crz(nsh)
      common /evopha/ nphase                ! ion part pms TD Fev.2019
      common /ovsh/ novlim(nmaxconv,ntypeconv)
      
c Ajout TD (18/01/2019)
      common /diffutest/ zmeanO16(nsh),dens_elec(nsh) ,coulomb(nsh)
     $     ,Zpaq(nsh),Zpaq1(nsh),Zpaq2(nsh),Kpaq(nsh)
      double precision zmeanO16,dens_elec,coulomb
     $     ,Zpaq,Zpaq1,Zpaq2,Kpaq
c Fin ajout 
      
      dimension ratio_ion(nsh) ! Test TD Jan.2019
      dimension Zdiff(mspec,mspec),Zdiff1(mspec,mspec)
     &    ,Zdiff2(mspec,mspec),Kdiff(mspec,mspec),kappa_st(mspec,mspec)
      dimension xi(nsh,mspec), xsp(nsh,nis), abond(nsh,nis), xior(nsh
     $     ,mspec)
      dimension dlncdr(mspec)
      dimension anuc(nsp),znuc(nsp),cl(mspec,mspec),xthoul(mspec),
     $     a(mspec),z(mspec),ap(mspec),at(mspec),ax(mspec,mspec)
     $     ,abondthoul(mspec)
      dimension dc(mmax),xx(mmax,mmax),y(mmax,mmax),yy(mmax,mmax),
     &     k(mmax,mmax),gradx(mmax)
      dimension alpha(nmax),nu(nmax),gamma(nmax,nmax),delta(nmax,nmax)
     &     ,ga(nmax)
      dimension e_ax(mspec),g_ax(mspec)

c     data (elemthoul(i),i=1,mspec-1) / ih1,ihe4 /
c      data (elemthoul(i),i=1,mspec-1) / 2,5 /
c      data (elemthoul(i),i=1,mspec-1) / 2,5,13,17,20 /
c      data (elemthoul(i),i=1,mspec-1) / 2,3,4,5,6,7,13,17,18,20,26,27,28
c     &     ,31,33,40,41,44,45,49 /
      
c.    ih1,ih2,ihe3,ihe4,ili6,ili7,ibe9,ib10,ib11,ic12,ic13,in14,in15,io16,
c.    io17,io18,if19,ine20,ine21,ine22,ina23,img24,img25,img26,ial27,isi28,
c.    ip31,is32,icl35
c$$$      data (elemthoul(i),i=1,mspec-1) / 2,3,4,5,6,7,10,11,12,13,14,17,18
c$$$     $     ,20,21,22,24,26,27,28,31,33,35,36,40,41,44,45,48,49,51 /

c     ih1,ih2,ihe3,ihe4,ili6,ili7,ibe9,ib10,ib11,ic12,ic13,in14,in15,io16,
c.    io17,io18,if19,ine20,ine21,ine22,ina23,img24,img25,img26,ial27,isi28,isi29,
c.    isi30,ip31,is32,is33,is34,icl35,is36,icl37
      data (elemthoul(i),i=1,mspec-1) / 2,3,4,5,6,7,10,11,12,13,14,17,18
     $     ,20,21,22,24,26,27,28,31,33,35,36,40,41,42,43,44,45,46,47,49
     $     ,50,52 /    !35 elem
***   Definition of variable specific to this routine
      amu = 1.6605402d-24

***   initialisation matrix Vdiffthoul
      if (nsh.eq.1) vdiffthoul(1:nsh,1:mspec) = 0.d0
      i = 0.d0
      j = 0.d0
      l = 0.d0

      do i=1,mspec-1
         xthoul(i) = xsp(jshell,elemthoul(i))
         abondthoul(i) = abond(jshell,elemthoul(i))
         gradx(i) = dlnxdr(jshell,elemthoul(i))         
         a(i) = anuc(elemthoul(i)) !Associated atomic masses
c$$$!     Associated charges (H1,H2,He3,He4,Li6,Li7,C12,C13,N14,N15,O16
c$$$!     ,O17,O18,F19,Ne20,Ne21,Ne22,Na23,Mg24,Mg25,Mg26,Al27,Si28,P31,S32,CL35)
      enddo
      if (nphase.eq.1) then
         do i=1,mspec-1     
            if (elemthoul(i).eq.5.or.elemthoul(i).eq.13
     $           .or.elemthoul(i).eq.17.or.elemthoul(i).eq.20.or.
     $           elemthoul(i).eq.26) then
               z(i) = zmean(jshell,elemthoul(i)) ! ionisation partielle
               if (z(i).lt.0.01d0*znuc(elemthoul(i))) z(i) = 0.01d0
     $              *znuc(elemthoul(i))
            else
               z(i) = znuc(elemthoul(i))
            endif
         enddo
      else if (nphase.gt.1) then ! refaire les tests sur les ionisations
         do i=1,mspec-1
            xthoul(i) = xsp(jshell,elemthoul(i))
            abondthoul(i) = abond(jshell,elemthoul(i))
            gradx(i) = dlnxdr(jshell,elemthoul(i))         
            a(i) = anuc(elemthoul(i))
            if (elemthoul(i).eq.5.or.elemthoul(i).eq.13
     $           .or.elemthoul(i).eq.17.or.elemthoul(i).eq.20
     $           .or.elemthoul(i).eq.3.or.elemthoul(i)
     $           .eq.4.or.elemthoul(i).eq.6.or.elemthoul(i).eq.7.or.
c     $           elemthoul(i).eq.11.or.elemthoul(i).eq.12.or.    ! Be9,Be10,B11
c     $           elemthoul(i).eq.10.or.
     $           elemthoul(i).eq.26.or.!elemthoul(i).eq.2.or.
     $           elemthoul(i).eq.31.or.elemthoul(i).eq.33.or.
     $           elemthoul(i).eq.40.or.elemthoul(i).eq.14.or.
     $           elemthoul(i).eq.18.or.elemthoul(i).eq.21.or.
     $           elemthoul(i).eq.22) then !.or.elemthoul(i).eq.2) then  H1
c     $           .or.elemthoul(i).eq.44  ! P
c     $           .or.elemthoul(i).eq.45.or.elemthoul(i).eq.49) then  ! S Cl
c            .eq.35.or.elemthoul(i).eq.36) then  ! Mg25, Mg26
               z(i) = zmean(jshell,elemthoul(i)) ! ionisation partielle
               if (z(i).lt.0.01d0*znuc(elemthoul(i))) z(i) = 0.01d0
     $              *znuc(elemthoul(i))
            else
               z(i) = znuc(elemthoul(i))
            endif
         enddo
      endif
      
      xthoul(mspec) = 0.d0
      abondthoul(mspec) = 0.d0
      gradx(mspec) = 0.d0
      a(mspec) = masselec/amu
      z(mspec) = -1.d0
c      print *, jshell,z(4),z(10)
c      if (jshell.eq.800) stop
      !!!!!!MODIF THOUL AVRIL 2018
C     Initialize parameters:
      
      k0 = 2.d0
      n = 2*mspec+2
      do i = 1,mspec
         dc(i) = 0.d0
      enddo
      cc = 0.d0
      ac = 0.d0          

c     calculate concentrations from mass fractions:
      ZXA = 0.d0
      do i = 1,mspec-1
         ZXA = ZXA+z(i)*xthoul(i)/a(i)
      enddo

      do i = 1,mspec-1
         dc(i) = xthoul(i)/(a(i)*ZXA)
      enddo
      dc(mspec) = 1.d0

C Calculate CC and AC:
         
c$$$      do i = 1,mspec
c$$$         cc = cc+dc(i)
c$$$         ac = ac+a(i)*dc(i)
c$$$      enddo
      cc = sum(dc(:))       ! quick
      ac = sum(a(:)*dc(:))
C Calculate the mass fraction of electrons:
      xthoul(mspec) = a(mspec)/ac

c calculate density of electrons (NE) from mass density (RHO):
      NE=rhosh/(amu*ac)
      xthoul(mspec) = NE*a(mspec)
      abondthoul(mspec) = xthoul(mspec)*muizero(jshell)/a(mspec)
c$$$      stop
c      NE = rhosh/(amu*muesh)   ! pour Thoul modif + !!
c calculate interionic distance (AO): 
      NI=0.d0
      NIZI2 = 0.d0   ! modif ++ Fev.2019
      DO I=1,mspec-1
c         NI=NI+dc(I)*NE
         NI=NI+(xthoul(i)/anuc(i)) ! modif + Fev.2019
         NIZI2=NIZI2+(xthoul(i)/anuc(i))*Z(I)**2 ! modif + Fev.2019
      ENDDO
      NI = NI*rhosh/amu         ! modif + Fev.2019
      NIZI2 = NIZI2*rhosh/amu+NE ! modif + Fev.2019

c      AO=(0.23873/NI)**(1./3.)
      AO = (3.d0/(4.d0*pi*NI))**(1.d0/3.d0)
c calculate Debye length (LAMBDAD):
      CZ=0.d0
c$$$      DO I=1,mspec
c$$$	 CZ=CZ+dc(I)*z(I)**2
c$$$      ENDDO
      CZ = sum(dc(:)*z(:)**2)   ! quick
c      LAMBDAD=SQRT(tsh*boltz/(4.d0*pi*ech**2*NE*CZ))
      LAMBDAD=DSQRT(tsh*boltz/(4.d0*pi*ech**2*NIZI2)) ! modif + Fev.2019
c calculate LAMBDA to use in Coulomb logarithm:
      LAMBDA=MAX(LAMBDAD,AO)
c      if (jshell.eq.2) print *, tsh,boltz,pi,ech,LAMBDAD,rhosh
c     if (jshell.eq.2) stop
c$$$      do i=1,nmod
c$$$c$$$  write(65,*) ratneh(i)
c$$$         print *, muizero(i)
c$$$      enddo
c$$$  stop
c$$$      if (jshell.eq.315) print *, abond(315,2),rhosh, tsh
c$$$     $     ,muizero(jshell),NE
c$$$      if (jshell.eq.315) stop
      do i=1,mspec!-1
         do j=1,mspec           !-1
            call paq (a(j),z(j),abondthoul(j),a(i),z(i)
     $           ,abondthoul(i),rhosh,tsh,dmic,ddth
     $           ,altpi,altep,diffpaqt,ratneh(jshell)
     $           ,muizero(jshell),LAMBDA,abondthoul(:),z(:)
     $           ,NE)

            Zdiff(i,j) = Z_st
            Zdiff1(i,j) = Z_st1
            Zdiff2(i,j) = Z_st2
            Kdiff(i,j) = K_st
         enddo
      enddo
c$$$      do i=1,mspec!-1
c$$$         do j=1,mspec
c$$$            write(49,'(1x,i4,4(1x,1pe11.4))') jshell,Zdiff(i,j)
c$$$     $           ,Zdiff1(i,j),Zdiff2(i,j),Kdiff(i,j)
c$$$         enddo
c$$$      enddo
c$$$      print *, 'jshell',jshell,Zdiff(4,4),Zdiff1(4,4),Zdiff2(4,4)
c$$$     $     ,Kdiff(4,4)
c$$$      stop
!!!!!!FIN MODIF THOUL AVRIL 2018
c$$$      do jshell=1,nmod
c$$$         print *,jshell,'ratneh',ratneh(jshell)
c$$$         write(39,*) jshell,ratneh(jshell)
c$$$      enddo
c$$$      stop
*     Coulomb logarithms (see equation 9 Thoul et al. 1994)
      do i = 1,mspec
	 do j = 1,mspec
	    xij = 2.3939d3*tsh*LAMBDA/abs(z(i)*z(j))     ! 4*boltz/ech**2 = 2.3939d3 (cgs)
	    cl(i,j) = 0.81245d0*log(1.d0+0.18769d0*xij**1.2d0)
	 enddo
      enddo

***   Normalisation Paquette's coeffecient for Thoul
      k0thoul = 1.41d-25*tsh**(-1.5)*NE*NE ! from MESA (cgs)
c      k0thoul = 1.41d-25*tsh**(-1.5)*(NE*1d6)**2 ! from MESA (SI)
c      k0thoul = 1.144d-40*(tsh*1d-7)**(-1.5)*(NE*1d-2)**2 ! from THOUL (thoul units)
c      k0thoul = 1.144d-40*tsh**(-1.5)*NE**2    ! from THOUL (cgs units)
c      kappa_st = Kdiff/K0
c      print *, 'K0',K0,'kappa_st',kappa_st
***   End normalisation Paquette's coefficents

      do i = 1,mspec
         do j = 1,mspec
            xx(i,j) = a(j)/(a(i)+a(j))
            y(i,j) = a(i)/(a(i)+a(j))
c            yy(i,j) = 3.0d0*y(i,j)+1.3d0*xx(i,j)*a(j)/a(i) ! sans paquette
            yy(i,j) = 3.0d0*y(i,j)+Zdiff1(i,j)*xx(i,j)*a(j)/a(i)
            k(i,j) = 1.d0*cl(i,j)*sqrt(a(i)*a(j)/(a(i)+a(j)))*dc(i)*
     &           dc(j)*z(i)**2*z(j)**2 ! eq. 37 Thoul et al. 94
         enddo
      enddo

c$$$      do i = 1,mspec-1
c$$$         do j = 1,mspec-1
c$$$            kappa_st(i,j) = Kdiff(i,j)/k0thoul
c$$$            k(i,j) = kappa_st(i,j)
c$$$c     if (k(i,j).ne.k(i,j)) print *,'i',i,'j',j, 'thoul',k(i,j)
c$$$         enddo
c$$$      enddo

      do j = 1,mspec-1          ! quick
         do i = 1,mspec-1
            kappa_st(i,j) = Kdiff(i,j)/k0thoul
            k(i,j) = kappa_st(i,j)
         enddo
      enddo
     
c..   Definition of the He4 index
      do i = 1,mspec
         if (elemthoul(i).eq.5) Heindex = i
      enddo
      
c..   Initialize terms in Eq.(30) to (32)
      alpha(1:n) = 0.d0
      nu(1:n) = 0.d0
      gamma(1:n,1:n) = 0.d0

      do i = 1,mspec
         alpha(i) = dc(i)/cc
         do j = 1,mspec
c.. sum only if elt different from electrons (index mspec) - modif following Hu et al. 2011             
            if (j.ne.mspec) then
               gamma(i,j) = -dc(j)/cc
               if (j.eq.i) then
                  gamma(i,j) = gamma(i,j)+1.d0
               endif
               gamma(i,j) = gamma(i,j)*dc(i)/cc
            endif
         enddo
      enddo
      
      do i = mspec+1,n-2
         nu(i) = 2.5d0*dc(i-mspec)/cc
      enddo

c...  Intialize delta
c$$$      do i = 1,n
c$$$         do j = 1,n
c$$$            delta(i,j) = 0.d0
c$$$         enddo
c$$$      enddo
      delta(:,:) = 0.d0 ! quick
      
c..   Compute coefficients on the rhs of equation 21 from Thoul et al. 1994
c...  Set of Eq. (33)      
      do i = 1,mspec
         do j = 1,mspec
            if (j.eq.i) then
               do l = 1,mspec
                  if(l.ne.i) then
                     delta(i,j) = delta(i,j)-k(i,l)
                  endif
               enddo
            else
               delta(i,j) = k(i,j)
            endif
         enddo

         do j = mspec+1,n-2
            if (j-mspec.eq.i) then
               do l = 1,mspec
                  if (l.ne.i) then
c                     delta(i,j) = delta(i,j)+0.6d0*xx(i,l)*k(i,l)
                     delta(i,j) = delta(i,j)+Zdiff(i,l)*xx(i,l)*k(i,l)
                  endif
               enddo
            else
c               delta(i,j) = -0.6d0*y(i,j-mspec)*k(i,j-mspec)
               delta(i,j) = -Zdiff(i,j-mspec)*y(i,j-mspec)*k(i,j-mspec)
            endif
         enddo
         
         delta(i,n-1) = dc(i)*z(i)
         
         delta(i,n) = -dc(i)*a(i)
      enddo

c...  Set of Eq. (34)      
      do i = mspec+1,n-2
         do j = 1,mspec
            if (j.eq.i-mspec) then
               do l = 1,mspec
                  if (l.ne.i-mspec) then
c                     delta(i,j) = delta(i,j)+1.5d0*xx(i-mspec,l)*
c     &                    k(i-mspec,l)
                     delta(i,j) = delta(i,j)+2.5d0*Zdiff(i-mspec,l)*xx(i
     &                    -mspec,l)*k(i-mspec,l)
                  endif
               enddo
            else
c               delta(i,j) = -1.5d0*xx(i-mspec,j)*k(i-mspec,j)
               delta(i,j) = -2.5d0*Zdiff(i-mspec,j)*xx(i-mspec,j)*k(i
     &              -mspec,j)
            endif
         enddo
         
         do j = mspec+1,n-2
            if (j-mspec.eq.i-mspec) then
               do l = 1,mspec
                  if (l.ne.i-mspec) then
c                     delta(i,j) = delta(i,j)-y(i-mspec,l)*k(i-mspec,l)*
c     &                    (1.6d0*xx(i-mspec,l)+yy(i-mspec,l))
                     delta(i,j) = delta(i,j)-y(i-mspec,l)*k(i-mspec,l)*
     &                    (0.8d0*Zdiff2(i-mspec,l)*xx(i-mspec,l)+yy(i
     &                    -mspec,l))
                  endif
               enddo
c               delta(i,j) = delta(i,j)-0.8d0*k(i-mspec,i-mspec)
               delta(i,j) = delta(i,j)-0.4d0*Zdiff2(i-mspec,i-mspec)*k(i
     &              -mspec,i-mspec)
            else
c               delta(i,j) = 2.7d0*k(i-mspec,j-mspec)*xx(i-mspec,j-mspec)
c     &              *y(i-mspec,j-mspec)
               delta(i,j) = k(i-mspec,j-mspec)*xx(i-mspec,j-mspec)*y(i
     &              -mspec,j-mspec)*(3d0+Zdiff1(i-mspec,j-mspec)-0.8d0
     &              *Zdiff2(i-mspec,j-mspec))
            endif
         enddo
         
         delta(i,n-1) = 0.d0
         
         delta(i,n) = 0.d0
      enddo

c...  Set of Eq. (35)      
      do j = 1,mspec
         delta(n-1,j) = dc(j)*z(j)
      enddo
      do j = mspec+1,n
         delta(n-1,j) = 0.d0
      enddo

c...  Set of Eq. (36)      
      do j = 1,mspec
         delta(n,j) = dc(j)*a(j)
      enddo
      do j = mspec+1,n
         delta(n,j) = 0.d0
      enddo
      
c      if (jshell.eq.2) print *, 'delta',delta,'indx',indx,n,nmax,d
C Inverse the system for each possible right-hand-side, i.e.,
C if alpha is the r.h.s., we obtain the coefficient A_p
C if nu    ---------------------------------------- A_T
C if gamma(i,j) ----------------------------------- A_Cj
C 
C If I=1, we obtain the hydrogen diffusion velocity
C If I=2, ------------- helium   ------------------
C If I=3,mspec-1, --------- heavy element -------------
C If I=mspec, ------------- electrons -----------------
C For I=mspec,2mspec, we get the heat fluxes
C For I=N-1, we get the electric field
C For I=N, we get the gravitational force g
c      CALL MKL_SET_NUM_THREADS(4)
      CALL dgetrf(n,n,delta,nmax,indx,d)
c$$$      print *, 'delta',delta,'indx',indx,n,nmax,d
c$$$      stop
      CALL dgetrs('N',n,1,delta,nmax,indx,alpha,nmax,d)
      CALL dgetrs('N',n,1,delta,nmax,indx,nu,nmax,d)
      do j = 1,n
         do i = 1,n
            ga(i) = gamma(i,j)
         enddo
         CALL dgetrs('N',n,1,delta,nmax,indx,ga,nmax,d)
         do i = 1,n
            gamma(i,j) = ga(i)
         enddo
      enddo
c      CALL MKL_SET_NUM_THREADS(1)
c      if (jshell.eq.2) print *, alpha(4),dc(4),cc
c      if (jshell.eq.4) stop
C The results for the coefficients must be multiplied by p/K_0:

      do i = 1,mspec
         alpha(i) = alpha(i)*k0*ac*cc
         nu(i) = nu(i)*k0*ac*cc
         do j = 1,mspec
            gamma(i,j) = gamma(i,j)*k0*ac*cc
         enddo
      enddo

      ! Test ajout E et g
      e_ap = alpha(n-1)*k0*ac*cc
      g_ap = alpha(n)*k0*ac*cc

      e_at = nu(n-1)*k0*ac*cc
      g_at = nu(n)*k0*ac*cc

      do i = 1,mspec
         e_ax(i) = gamma(n-1,i)*k0*ac*cc 
         g_ax(i) = gamma(n,i)*k0*ac*cc 
      enddo
c      print *, 'eAP',e_ap,'eAT',e_at,'eAX',e_ax(:)
c      print *, 'gAP',g_ap,'gAT',g_at,'gAX',g_ax(:)
c      stop
! Fin test ajout E et g

      do i = 1,mspec
         ap(i) = alpha(i)
         at(i) = nu(i)
         do j = 1,mspec
            ax(i,j) = gamma(i,j)
         enddo
      enddo
      
      temp = 0.d0
      do i = 1,mspec-1
         temp = temp+z(i)*xthoul(i)/a(i)*gradx(i)
      enddo
      do j = 1,mspec-1
         dlncdr(j) = temp
      enddo
      temp = 0.d0
      do i = 1,mspec-1
         temp = temp+z(i)*xthoul(i)/a(i)
      enddo
      do j = 1,mspec-1
         dlncdr(j) = gradx(j)-dlncdr(j)/temp
      enddo
      !     Test ajout E et G (suite)
         e_tot = e_ap*dlnpdr(jshell)+e_at*dlntdr(jshell)
         g_tot = g_ap*dlnpdr(jshell)+g_at*dlntdr(jshell)
         sum2 = 0.d0
         sum3 = 0.d0
         do j = 1,mspec
            sum2 = sum2 + e_ax(j)*dlncdr(j)
            sum3 = sum3 + g_ax(j)*dlncdr(j)
         enddo
         e_tot = e_tot + sum2
         g_tot = g_tot + sum3
c         print *,e_tot,g_tot
c         stop
!     Fin test ajout E et G (suite)
         
      do i  =  1,mspec
         sum1 = 0.d0
         do j  =  1,mspec-1
c            if (i.ne.j) then   
               sum1  =  sum1 + ax(i,j)*dlncdr(j)
c            endif
         enddo
        
         xior(jshell,i) = ap(i)*dlnpdr(jshell)+at(i)*dlntdr(jshell)  ! pour comparaison modele montreal ie sans grad. c
         xi(jshell,i) = ap(i)*dlnpdr(jshell)+at(i)*dlntdr(jshell)+sum1
       
c      xi(jshell,i)  =  xi(jshell,i)*tsh**2.5d0/rhosh*8.08998904d-16
c     xi(jshell,i) = xi(jshell,i)*tsh**2.5d0/rhosh*??? !(((1e-7)^2.5)/1e-2)* (Rsol(cm)/6e13/365.25/24/3600) = ???
         ! Vitesse sans dimension
         xi(jshell,i) = xi(jshell,i)*(tsh*1d-7)**(2.5d0)/(rhosh*1d-2)
         xior(jshell,i) = xior(jshell,i)*(tsh*1d-7)**(2.5d0)/(rhosh*1d
     $        -2)
         ! Vitesse en cgs (avec tau0 = 6E-13 annees, see part 3 of Thoul et al. 1994)
         xi(jshell,i) = (xi(jshell,i)*rsun)/(6d13*(365.25*24*3600))
         xior(jshell,i) = (xior(jshell,i)*rsun)/(6d13*(365.25*24*3600))
         vdiffthoul(jshell,i) = xi(jshell,i)
         dmicthoul(jshell,i) = abs(ax(i,i)
     &        *(tsh*1d-7)**(2.5d0)*rsun**2
     &        /(rhosh*1d-2*6d13*(365.25*24*3600)))
c$$$         if (jshell.ne.1) then
c$$$            if (crz(jshell-1).eq.-3.and.crz(jshell).eq.-2) then
c$$$               vdiffthoul(jshell,i) = vdiffthoul(jshell-1,i)
c$$$            endif
c$$$            if (abs(vdiffthoul(jshell,i)).gt.10) vdiffthoul(jshell,i) =
c$$$     $           vdiffthoul(jshell-1,i)
c$$$         endif
         
c$$$          if (i.eq.4) write(669,'(1x,i4,4(1x,1pe11.4))'),jshell
c$$$  $        ,vdiffthoul(jshell,4),xior(jshell,4),ap(4),at(4)

!     Attribution of Mg 24 velocity to Mg 25 and Mg 26 in aim to consider Mg partial ionisation
!     Modif TD (24/05/2019)
         if (i.eq.35.or.i.eq.36) then
            xi(jshell,i) = xi(jshell,33)
            xior(jshell,i) = xior(jshell,33)
         endif
!     Fin modif    
! Add by TD Mai.2018
         if (i.eq.4) then
            xvdiffhe(jshell)=vdiffthoul(jshell,i)
            xvdiffhenoc(jshell)=xior(jshell,i)
            xdmiche(jshell) = dmicthoul(jshell,i)
         endif
c$$$         if (i.eq.4) write(669,'(1x,i4,5(1x,1pe11.4))'),jshell
c$$$     $        ,vdiffthoul(jshell,4),xior(jshell,4),xvdiffhe(jshell)
c$$$     $        ,xvdiffhenoc(jshell),xdmiche(jshell)
         if (i.eq.10) then
            xvdiffc12(jshell)=vdiffthoul(jshell,i)
            xvdiffc12noc(jshell)=xior(jshell,i)
         endif
         if (i.eq.12) then 
            xvdiffn14(jshell)=vdiffthoul(jshell,i)
            xvdiffn14noc(jshell)=xior(jshell,i)
         endif
         if (i.eq.14) then
            xvdiffo16(jshell)=vdiffthoul(jshell,i)
            xvdiffo16noc(jshell)=xior(jshell,i)
         endif
         if (i.eq.18) then 
            xvdiffne20(jshell)=vdiffthoul(jshell,i)
            xvdiffne20noc(jshell)=xior(jshell,i)
         endif
         if (i.eq.21) then  
            xvdiffna23(jshell)=vdiffthoul(jshell,i)
            xvdiffna23noc(jshell)=xior(jshell,i)
         endif
         if (i.eq.11) then  
            xvdiffc13(jshell)=vdiffthoul(jshell,i)
            xvdiffc13noc(jshell)=xior(jshell,i)
         endif
         if (i.eq.13) then  
            xvdiffn15(jshell)=vdiffthoul(jshell,i)
            xvdiffn15noc(jshell)=xior(jshell,i)
         endif
         if (i.eq.15) then  
            xvdiffo17(jshell)=vdiffthoul(jshell,i)
            xvdiffo17noc(jshell)=xior(jshell,i)
         endif
         if (i.eq.16) then  
            xvdiffo18(jshell)=vdiffthoul(jshell,i)
            xvdiffo18noc(jshell)=xior(jshell,i)
         endif
      enddo

      return
      
      end                                   

      
*-----------------------------------------------------------------------

      SUBROUTINE COEFONDES_SURF (ndb,ndt)

*-----------------------------------------------------------------------
*     Calcule le coefficient de transport associe aux ondes
C     emises par l'enveloppe convective (zone convective de surface)
*     qui doit etre applique au moment et aux especes chimiques.
*     COEFONDES est appelee par DIFFUSION.
*     Entrees : ndt  : premiere couche radiative utile pour le calcul
*     ndb   : derniere couche radiative utile pour le calcul
*     (on peut avoir ndt = ndt-1, ndb = ndb+j)
C Modifs CC ondes (22/10/04)
*     Sorties (par common) :
*        Dondes : coefficient de transport associe aux ondes applique a omega
*        Dondeschim : coefficient de transport associe aux ondes applique aux
*                     elements chimiques
*
* T. Decressin (octobre 2010) : varibles definies du centre vers surf. *
*
*     Derniere version : 29 octobre 2009
*-----------------------------------------------------------------------
      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.grad'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.transp'
      include 'evolcom.var'
      include 'evolcom.igw'
      include 'evolcom.mod'

      integer ndb,ndt
      integer j,i_os,jj

      double precision amp,ssigma,rmtu,dr,xlgD,dmin
      double precision rdensmi_os,xcpmi_os,xKmi_os

      write (*,*) 'Entree dans COEFONDES_SURF'

*     Remise a zero du coefficient de transport associe aux ondes

      do j = 1,nmod
         Dondes(j) = 0.d0
         Dondeschim(j) = 0.d0
      enddo

*     Expression de Talon & Charbonnel (2004)
*     Premiere approximation
c      rtot = r(nmod)
c      rmtu = r(ndt)/rtot
C
C Modif CC (07/07) : coefficient calibre sur le coefficient SLO depuis
C la zams jusqu'au sommet de la RGB
C Modif ST (03/07/07)
C lumwave(j)-lumwave(j+1) represente la quantite de moment cinetique
C                         deposee dans la couche - il faut divisier
C                         par la masse de la couche pour obtenir une
C                         quantite physiquement significative
C CC : Attention, lumwave est rempli dans le sens du code de Geneve,
C cad de la surface vers le centre. 
C Il faut donc inverser le tableau pour ecrire Dondes qui lui est
C rempli dans le sens de STAREVOL, cad du centre vers la surface.
      do j = ndb,ndt
         Dondeschim(j) = (dlog10(xKt(j))+2*dlog10(dabs(lumwave(j))+
     &        1e-30)-
     &        68.d0)*dexp(-(r(ndt)-r(j))/(0.0003*(r(nmod)-r(ndt))))
         if (Dondeschim(j).lt.-1d-10) then
            Dondeschim(j) = 0.d0
         else 
            Dondeschim(j) = 10**(Dondeschim(j))
         endif
c        print *,j,Dondeschim(j),-(r(ndt)-r(j))/(0.0003*(r(nmod)-r(ndt)))
c     &        ,(dlog10(xKt(j))+2*dlog10(dabs(lumwave(j)))-
c     &        68.d0),xKt(j),lumwave(j)
c     &        (r(ndt)-r(j))/(0.0003*(1.d0-r(ndt)))

c         Dondeschim(j) = dabs(0.001d0*((lumwave(j)-lumwave(j+1))/
c     &        dm(j))**2)
c         print *,j,Dondeschim(j)
c       print *,"j,Dondeschim(j),lumwave(j)",j,Dondeschim(j),lumwave(j)
      enddo
c      stop
C Fin modif CC <--
C Modifs CC ondes (22/10/04)
C CC (8/09/09) On laisse Dondes(j) = 0.d0
C      do j = ndb,ndt
C         Dondes(j) = Dondeschim(j)
C      enddo
C Fin Modifs CC ondes (22/10/04)

C      do j = ndb,ndt
C         if (Dondeschim(j).ge.1.d-10) then 
C             write(*,*) 'j, Dondeschim(j) = ',j,Dondeschim(j)
C         endif
C      enddo
C      write(*,*) 'Dondes(ndb), Dondes(ndt) = ',Dondes(ndb),Dondes(ndt)
C      write(*,*) 'Dondeschim(ndb), Dondeschim(ndt),Dondeschim(ndt-10) = 
C     &           ',Dondeschim(ndb),Dondeschim(ndt),Dondeschim(ndt-10)

      return
      end


*-----------------------------------------------------------------------

      SUBROUTINE DTURBU (ndb,ndt,xomega,cd)

*-----------------------------------------------------------------------
*     Calcule le coefficient D de l'equation de diffusion.
*     Integration effectuee en masse.
*     DTURBU est appelee par DIFFUSION.
*     Derniere version : 7 juillet 1995
*     Entrees : ndt  : premiere couche radiative utile pour le calcul 
*     ndb   : derniere couche radiative utile pour le calcul
*     (on peut avoir ndt = ndt-1, ndb = ndb+j)
*     xomega  : vitesse angulaire.
*     Sorties : cd   : coefficient de l'equation de diffusion.
*     ou 1-omega*omega/2/pi/G/ro < 0.
*-----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.grad'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.transp'
      include 'evolcom.var'

      integer ndb,ndt
      integer ndbp1
      integer j

      double precision xomega,cd
      double precision ggama,xcd,xxcd

      dimension cd(nsh)

*     Coefficient de diffusion Zahn (1987)

*     Initialisations.

*     L'initialisation de ggama sera a reprendre a partir de Li 
*     dans nouvelle version du Soleil

c     ggama = 1.d0/7.d0
      ggama = 1.d0

*     Remise a zero du coefficient de diffusion turbulente.

      do j = 1,nmod
         cd(j) = 0.d0
      enddo

*     Calcul de D (D = cd), coefficient de diffusion turbulente.

*     Expression de Zahn (1987)

      if (ndb.eq.1) ndbp1 = 2
      if (ndb.ne.1) ndbp1 = ndb

      do j = ndbp1,ndt
         xcd = (xomega/g)**2*(r(j)**6)
         xcd = xcd/(abad(j)-abrad(j))
         xcd = xcd*lum(j)/(m(j)**3)
         xxcd = ggama
         cd(j) = xcd*xxcd
         cd(j) = abs(cd(j))
      enddo

      return
      end



      SUBROUTINE DTURBURGB (ndbmel,ndt,diffstcd,cd)

*-----------------------------------------------------------------------
*     Coefficient D de l'equation de diffusion :
*	Parametrisation pour tests melange sur RGB,
*	Coefficient de diffusion constant pour RGB
*     Integration effectuee en masse.
*     DTURBURGB est appelee par DIFFUSION.
*     Derniere version : 4 mai 1999
*     Entrees : ndt  : premiere couche radiative utile pour le calcul 
*     ndb   : derniere couche radiative utile pour le calcul
*     (on peut avoir ndt = ndt-1, ndb = ndb+j)
*     diffstcd : valeur de diffst lue dans starevol.par, utilisee dans
*		 ce cas pour cd.
*     Sorties : cd   : coefficient de l'equation de diffusion.
*-----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.teq'
      include 'evolcom.transp'

      integer ndbmel,ndt
      integer ii

      double precision diffstcd,cd

      dimension cd(nsh)

      do ii = 1,ndbmel-1
         cd(ii) = 0.d0
      enddo
      do ii = ndbmel,ndt
         cd(ii) = diffstcd
      enddo
      do ii = ndt+1,nmod
         cd(ii) = 0.d0
      enddo

      return
      end



      SUBROUTINE DTURBUCC95 (ndt,ndb,xomega,cd)

C=======================================================================
C
C Calcule le coefficient D de l'equation de diffusion.
C
C Integration effectuee en masse.
C
C=======================================================================
C
C DTURBU est appelee par DIFFUSION.
C ndb est en fait isup.
C
C Auteur : C.Charbonnel.       
C Modifs : Y.Gaige, O.Richard
C Derniere version : 13 aout 1998
C=======================================================================
C
C Entrees : ndt        : premiere couche radiative.
C           ndb      : derniere couche radiative.
C           xomega     : vitesse angulaire.
C
C Sorties : cd         : coefficient de l'equation de diffusion.
C
C=======================================================================

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.therm'
      include 'evolcom.transp'
      include 'evolcom.var'

      integer ndt,ndb
c     integer i
      integer j

      double precision xomega,cd
      double precision coeffh,djdt,cdt,dlnr,dlnrom

      dimension cd(nsh)

C Remise a zero du coefficient de diffusion turbulente.

      do j = 1,nsh
         cd(j) = 0.d0
      enddo

C Calcul de D (D=cd), coefficient de diffusion turbulente.

      coeffh = 1.d0
      coeff = 3.d0*coeffh/(80.d0*pi)
C Pour M = 1.5 Msol ; dj/dt = (j(bump)-j(bas RGB))/(t(bump)-t(bas RGB))
C     djdt = 2.381773462d34  
C Pour M = 0.83 Msol ; dj/dt = (j(tip)-j(bas RGB))/(t(tip)-t(bas RGB))
C     djdt = 6.375527207d34
      djdt = 2.385300702d34
      do j = ndt,ndb+1,-1
         cdt = coeff*abs(djdt)/(xomega*ro(j)*r(j)**3)
         if (r(j-1).gt.0.d0) then
            dlnr = log(r(j-1))-log(r(j+1))
            dlnrom = log(xomega*r(j-1)**2)-
     &           log(xomega*r(j+1)**2)
            cdt = cdt*2.d0/dlnrom*dlnr
            cd(j) = cdt
         endif
c        if (cd(j).gt.cd(j+1).and.j.lt.ndt) then
c           do i = j,ndb,-1
c              cd(i) = cd(j+1)
c           enddo
c           goto 300
c        endif
c        cd(j) = 1.d-1*cd(j)
      enddo

      return
      end



      SUBROUTINE DTURBUDENISSENKOV (ndt,ndb,cd)

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.therm'

      integer ndt,ndb
      integer ibottom,idown
      integer i,j

      double precision cd
      double precision lgcd,deltam 

      dimension cd(nsh),lgcd(nsh),deltam(nsh)

      common /coeff_diff/ ibottom,idown

C Remise a zero du coefficient de diffusion turbulente.

      do j = 1,nsh
         cd(j) = 0.d0
      enddo

C Calcul de D (D=cd), coefficient de diffusion turbulente.
 
*
* mr(j) est donnee en masse relative, i.e en fraction de masse stellaire
*     
      do i = ndt,ndb,-1
         deltam(i) = (mr(i)-mr(idown))/(mr(ndt)-mr(idown))
      enddo

      do j = ndt,ndb+1,-1
         if (deltam(j).ge.1.d-1) then
            lgcd(j) = 10.d0/3.d0*deltam(j)+6.5d0-1.d0/3.d0
         else
            lgcd(j) = 4.d0+25.d0*deltam(j)
         endif
         cd(j) = 10.d0**lgcd(j)
      enddo

      do j = ndt,ndb+1,-1
         cd(j) = 1.d0*cd(j)
      enddo 

      return
      end



      SUBROUTINE DTURBULITHIUM (ndt,ndb,cd)

*------------------------------------------------------------
* coefficient de diffusion pique en r=0.003 Rtot, avec deux 
* parties: lineaire, puis decroissance exponentielle au centre
*
* fait le 27 mars 2000
*------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.teq'
      include 'evolcom.transp'
      include 'evolcom.var'

      integer ndt,ndb
      integer ibottom,idown
      integer i,j

      double precision cd
      double precision multi,rlim

      dimension cd(nsh)

      common /coeff_diff/ ibottom,idown
 
      multi = r(nmod)
      rlim = 3.d-3*multi  
      idown = 0

      do j = 1,nsh
        cd(j) = 0.d0
      enddo

      do i = ndt,ndb,-1
         if (r(i).gt.rlim) then
           cd(i) = 1.d12/(rlim-r(ndt))*(5.d0*r(i)+(rlim-6.d0*
     &           r(ndt)))
         else
           if (r(i).le.rlim.and.r(i+1).gt.rlim) ibottom = i
c          cd(i) = 6.d12*exp(-r(i)/(5.d-4*multi)) 
c          cd(i) = 6.d12*exp(r(i)/(25.d-5*multi))  
           cd(i) = 1.d-10*exp(r(i)*29.422781d0/rlim)
         endif
      enddo
      do i = ndt,ndb,-1
        cd(i) = 5.d1*cd(i)
      enddo

      return
      end


      SUBROUTINE FIT_F (ap,psi,epsilon,yomega)

*-----------------------------------------------------------------------
*     FIT_F calcule les integrales de collision.
*-----------------------------------------------------------------------

      implicit none

      integer iaj,k,k1
      integer i,j

      double precision psi,epsilon,yomega
      double precision f,epsi
      double precision paq_c,paq_d
      double precision d0,d1

      logical ap

      common /fitcd/ paq_c(3,8,50),paq_d(8,50)

      dimension f(2,3),yomega(2,3)

      iaj = 0
      if (ap) iaj = 4

*     fit dans le cas ou: -7.0 < = psi < = 3.0

      if (psi.lt.3.d0) then
         k = int(5.d0*psi+36.d0)
         if (k.gt.50) then
             k=50
         else
            if (k.lt.1) then
            k=1
            endif
         endif 
         k1 = k+1
         d0 = psi-(-7.2d0+0.2d0*k)
         d1 = (-7.2d0+0.2d0*k1)-psi
         do j = 1,3
            f(1,j) = exp(paq_c(j,1+iaj,k)*d1**3+paq_c(j,2+iaj,k)*d0**3+
     &           paq_c(j,3+iaj,k)*d1+paq_c(j,4+iaj,k)*d0)
         enddo
         f(2,2) = exp(paq_d(1+iaj,k)*d1**3+paq_d(2+iaj,k)*d0**3+
     &        paq_d(3+iaj,k)*d1+paq_d(4+iaj,k)*d0)
      else
         epsi = exp(psi)

*     fit dans le cas: 3.0 < psi < 4.0 pour le potentiel repulsif
*     et pour psi > 4.0

         if (.not.ap.or.psi.gt.4.d0) then
            f(1,1) = 1.00141d0*epsi-3.18209d0
            f(1,2) = 0.99559d0*epsi-1.29553d0
            f(1,3) = 1.99814d0*epsi-0.64413d0
            f(2,2) = 1.99016d0*epsi-4.56958d0
         else

*     fit dans le cas: 3.0 < psi < 4.0 pour le potentiel attractif

            f(1,1) = 1.01101d0*epsi-3.19815d0
            f(1,2) = 1.04230d0*epsi-1.89637d0
            f(1,3) = 2.15672d0*epsi-2.81038d0
            f(2,2) = 2.08699d0*epsi-5.81444d0
         endif
      endif

*     calcul des integrales de collision

      do i = 1,3
         yomega(1,i) = f(1,i)*epsilon
      enddo
      yomega(2,2) = f(2,2)*epsilon

      return
      end


      SUBROUTINE PAQ (As,Zs,ys0,At,Zt,yt0,rho,xtemp,Dmic,ddth,altpi
     &     ,altep,diffpaqt,ratnehk,xmuizero,lambda,xpaq,zpaq,NE)

*-----------------------------------------------------------------------
*     PAQ permet de calculer le coefficient de diffusion microscopique 
*     de l'element t qui diffuse par rapport a l'element s (ou
*     inversement) par la methode de Paquette.
*     Paquette et al. (86) evite un certaine nombre d'approximation,
*     necessaire a Chapman et Cowling pour resoudre de facon analytique
*     les integrales de collision. Et font une resolution numerique de
*     ces integrales.
*     PAQ est appelee par DIFFPAQ
*     Derniere version: 20 fevrier 1997
*-----------------------------------------------------------------------

      implicit none

      include 'evolcom.cons'

      integer i,mspec,lmicro
c      parameter (mspec = 32)    ! ajout 28/02/19
      parameter (mspec = 36)
      double precision As,Zs,ys,ys0,At,Zt,yt,yt0,rho,xtemp,Dmic,ddth,
     &     diffpaqt,ratnehk,xmuizero
      double precision omegast,omegass,omegatt
      double precision altpi,altep
      double precision omegase,omegaee
      double precision e2,amu,xkt,xmh,xmms,xmmt,xmst,xMsts,xMstt
      double precision xmme,xmse,xMste
      double precision epsilse,epsilee
      double precision epsilon,epsilss,epsiltt,epsilst,xn,xne,xns,
     &     xnt,sumnz2,xlambdad,xlambdai,xlambda,gamaee,gamase,
     &     gamaet,gamast,gamass,gamatt,psist,psiss,psitt,alphast,
     &     psise,psiee,Dst2,Dst1,Dste1,Dste2,alphase
      double precision lambda,zpaq,xpaq,sumnz3,NE           ! modif nov.18 +  ! ajout 28/02/19

      logical ap

      dimension omegast(2,3),omegass(2,3),omegatt(2,3)
      dimension omegase(2,3),omegaee(2,3)
      dimension zpaq(mspec), xpaq(mspec)         ! ajout 28/02/19

      ys = ys0
      yt = yt0
      e2 = 2.307121d-19
      amu = 1.d0/avn

      xkt = boltz*xtemp
      xmh = 1.00787d0*amu
      xmms = As*amu
      xmmt = At*amu
      xmme = 5.485799d0 * 1.d-4 * amu
      xmst = xmms+xmmt
      xmse = xmms + xmme
      xMsts = xmms/xmst
      xMstt = xmmt/xmst
      xMste = xmme / xmse
c      print *, xkt,xmh,xmms,xmmt,xmme
c      stop
*     calcul du parametre liant les integrales de
*     collision au developpement (formule 66)

      epsilon = 0.25d0*dsqrt(pi)*e2*e2*xkt**(-1.5d0)
      epsilss = epsilon/dsqrt(xmms)*Zs**4
      epsiltt = epsilon/dsqrt(xmmt)*Zt**4
      epsilst = epsilss/dsqrt(2.d0*xMstt)*(Zt/Zs)**2
      epsilse = epsilon/dsqrt(2.d0*xmms*xMste)*Zs**2
      epsilee = epsilon/dsqrt(xmme)
      
      xn = rho/(xmuizero*xmh)   ! number density of ions
c      if (lmicro.eq.4) then
         xne = NE               ! THOUL change
c      endif
c      if (lmicro.eq.3) then
c         xne = xn*ratnehk          ! M&M change
c         print *, 'xn',xn,'ratnehk',ratnehk
c      endif
      
      xns = xn*ys              ! number density of atoms of element Zs in gas
      xnt = xn*yt

c$$$      if (lmicro.eq.3) then
c         sumnz2 = xns*Zs*Zs+xnt*Zt*Zt+xne ! Montmerle&Michaud
*     calcul de la longueur de debye: xlambdad
c         xlambdad = dsqrt(xkt/abs(4.d0*pi*e2*sumnz2)) ! M&M change
c$$$      else if (lmicro.eq.4) then
         sumnz3 = sum(xn*xpaq(:)*Zpaq(:)*Zpaq(:)) ! THOUL
         sumnz3 = sumnz3 + xne ! THOUL
*     calcul de la longueur de debye: xlambdad         
         xlambdad = dsqrt(xkt/abs(4.d0*pi*e2*sumnz3)) ! pour THOUL change
c$$$      endif                         ! modif Fev.2019
c      xlambdad = lambda                              ! modif nov.18
*     calcul de la distance inter-ionique moyenne

      xlambdai = (3.d0/(4.d0*pi*xn))**(1.d0/3.d0)

*     determination de la distance d'ecran

      xlambda = max(xlambdad,xlambdai)
c      print *, xlambdad-xlambdai
*     calcul de gamma (formule 67)

      gamaee = 4.d0*xkt*xlambda/e2
      gamase = gamaee/Zs
c      gamaet = gamaee/Zt
      gamast = gamaee/(Zs*Zt)
      gamass = gamaee/(Zs*Zs)
      gamatt = gamaee/(Zt*Zt)
      
*     calcul du paramettre utilise dans le
*     developpement des integrales de collision
*     (formule 68)

      psist = log(log(1.d0+gamast*gamast))
      psiss = log(log(1.d0+gamass*gamass))
      psitt = log(log(1.d0+gamatt*gamatt))
      psise = log(log(1.d0+gamase*gamase))
      psiee = log(log(1.d0+gamaee*gamaee))
     
*     calcul des integrales de collision
c      print *,'gamma',gamaee,gamase,gamast,gamass,gamatt
      ap = .false.
      call fit_f (ap,psist,epsilst,omegast)
      call fit_f (ap,psiss,epsilss,omegass)
      call fit_f (ap,psitt,epsiltt,omegatt)
      call fit_f (ap,psiee,epsilee,omegaee)
      ap=.true.
      call fit_f (ap,psise,epsilse,omegase)
      ap=.false.
      
      xn = xns + xnt               ! modif + 15/11/18
*     calcul des coefficients alpha_se pour la
*     diffusion thermique
      
c      print *, 'appel dalpha therm', pour Montmerle&Michaud 
      call dalpha (omegase,omegass,omegaee,xkt,xMsts,xMste,xmse,
     &     xn,ys,ratnehk,Dste1,Dste2,alphase)
c      print *, 'fin appel dalpha therm'      

      altep=-1.d0*alphase

*     calcul des coefficients de diffusion de Paquette
c      print *,omegast,omegass,omegatt,xkt,xMsts,xMstt,xmst,
c     &     xn,ys,yt,Dst1,Dst2,alphast
c     stop
      call dalpha (omegast,omegass,omegatt,xkt,xMsts,xMstt,xmst,
     &     xn,ys,yt,Dst1,Dst2,alphast)

      Dmic = Dst2
c      ddth = ys*yt*alphast
      altpi = alphast
c      diffpaqt = Dst1

      return
      end



      SUBROUTINE PAQHE (As,Zs,ys0,At,Zt,yt0,rho,xtemp,Dmic,ddth,altei
     &     ,altep,altpi,diffpaqt,ratnehk,xmuizero)

*-----------------------------------------------------------------------
*     PAQ permet de calculer le coefficient de diffusion microscopique 
*     de l'element t qui diffuse par rapport a l'element s (ou
*     inversement) par la methode de Paquette.
*     Paquette et al. (86) evite un certaine nombre d'approximation,
*     necessaire a Chapman et Cowling pour resoudre de facon analytique
*     les integrales de collision. Et font une resolution numerique de
*     ces integrales.
*     PAQHE est appelee par DIFFPAQ
*     Derniere version: 6 janvier 2003
*-----------------------------------------------------------------------

      implicit none

      include 'evolcom.cons'

      double precision As,Zs,ys,ys0,At,Zt,yt,yt0,rho,xtemp,Dmic,ddth,
     &     diffpaqt,ratnehk,xmuizero
      double precision omegast,omegass,omegatt
      double precision altei,altpi,altep
      double precision omegase,omegaee,omegate
      double precision e2,amu,xkt,xmh,xmms,xmmt,xmst,xMsts,xMstt
      double precision xMstes
      double precision xmme,xmse,xmte,xMste,xMset,xMstte
      double precision epsilse,epsilee,epsilte
      double precision epsilon,epsilss,epsiltt,epsilst,xn,xne,xns,
     &     xnt,sumnz2,xlambdad,xlambdai,xlambda,gamaee,gamase,
     &     gamaet,gamast,gamass,gamatt,psist,psiss,psitt,alphast,
     &     psise,psiee,Dst1,Dst2,Dste1,Dste2,Dstet1,Dstet2,psite,
     &     alphase,alphate,gamate

      logical ap

      dimension omegast(2,3),omegass(2,3),omegatt(2,3)
      dimension omegase(2,3),omegate(2,3),omegaee(2,3)

      ys = ys0
      yt = yt0
      e2 = 2.307121d-19
      amu = 1.d0/avn

      xkt = boltz*xtemp
      xmh = 1.00787d0*amu
      xmms = As*amu
      xmmt = At*amu
      xmme = 5.485799d0 * 1.d-4 * amu
      xmst = xmms+xmmt
      xmse = xmms + xmme
      xmte = xmmt + xmme
      xMsts = xmms/xmst
      xMstes = xmms/xmse
      xMstt = xmmt/xmst
      xMste = xmme / xmse
      xMset = xmme / xmte
      xMstte = xmmt / xmte

*     calcul du parametre liant les integrales de
*     collision au developpement (formule 66)

      epsilon = 0.25d0*dsqrt(pi)*e2*e2*xkt**(-1.5d0)
      epsilss = epsilon/dsqrt(xmms)*Zs**4
      epsiltt = epsilon/dsqrt(xmmt)*Zt**4
      epsilst = epsilss/dsqrt(2.d0*xMstt)*(Zt/Zs)**2
      epsilse = epsilon/dsqrt(2.d0*xmms*xMste)*Zs**2
      epsilee = epsilon/dsqrt(xmme)
      epsilte = epsilon/dsqrt(2.d0*xmmt*xMset)*Zt**2

      xn = rho/(xmuizero*xmh)
      xne = xn*ratnehk
      xns = xn*ys
      xnt = xn*yt
      sumnz2 = xns*Zs*Zs+xnt*Zt*Zt+xne

*     calcul de la longueur de debye: xlambdad

      xlambdad = dsqrt(xkt/abs(4.d0*pi*e2*sumnz2))

*     calcul de la distance inter-ionique moyenne

      xlambdai = (3.d0/(4.d0*pi*xn))**(1.d0/3.d0)

*     determination de la distance d'ecran

      xlambda = max(xlambdad,xlambdai)

*     calcul de gamma (formule 67)

      gamaee = 4.d0*xkt*xlambda/e2
      gamase = gamaee/Zs
      gamate = gamaee/Zt
      gamast = gamaee/(Zs*Zt)
      gamass = gamaee/(Zs*Zs)
      gamatt = gamaee/(Zt*Zt)

*     calcul du paramettre utilise dans le
*     developpement des integrales de collision
*     (formule 68)

      psist = log(log(1.d0+gamast*gamast))
      psiss = log(log(1.d0+gamass*gamass))
      psitt = log(log(1.d0+gamatt*gamatt))
      psise = log(log(1.d0+gamase*gamase))
      psiee = log(log(1.d0+gamaee*gamaee))
      psite = log(log(1.d0+gamate*gamate))

*     calcul des integrales de collision

      ap = .false.
      call fit_f (ap,psist,epsilst,omegast)
      call fit_f (ap,psiss,epsilss,omegass)
      call fit_f (ap,psitt,epsiltt,omegatt)
      call fit_f (ap,psiee,epsilee,omegaee)
      ap=.true.
      call fit_f (ap,psise,epsilse,omegase)
      call fit_f (ap,psite,epsilte,omegate)
      ap=.false.

*     calcul des coefficients alpha_se pour la
*     diffusion thermique
 
      call dalpha (omegase,omegass,omegaee,xkt,xMstes,xMste,xmse,
     &     xn,ys,ratnehk,Dste1,Dste2,alphase)

      altep=-1.d0*alphase

*     calcul des coefficients alpha_et pour la
*     diffusion thermique
 
      call dalpha (omegate,omegatt,omegaee,xkt,xMstte,xMset,xmte,
     &     xn,yt,ratnehk,Dstet1,Dstet2,alphate)

      altei=-1.d0*alphate

*     calcul des coefficients de diffusion de Paquette

      call dalpha (omegast,omegass,omegatt,xkt,xMsts,xMstt,xmst,
     &     xn,ys,yt,Dst1,Dst2,alphast)

      Dmic = Dst2
c      ddth = ys*yt*alphast
      altpi = alphast
c      diffpaqt = Dst1

      return
      end



      SUBROUTINE TACHOCLINE (ndt,ndb,coefDtacho)
C
C Cette routine calcule le coefficient de melange dans la tachocline
C selon Brun et al. (1999, ApJ 525, 1032)

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.therm'
      include 'evolcom.teq'
      include 'evolcom.transp'
      include 'evolcom.var'

      integer ndt,ndb
      integer i
      integer iliminftacho

      double precision coefDtacho,dtacho,xKtlocal,rliminftacho
c      double precision lgcd,deltam 
      double precision NBV,Q4,mu4,omegaBZC,coefNuH,nuHlocal
      double precision deltalocal,xsilocal,explocal
      double precision cos2local

      dimension coefDtacho(nsh),xKtlocal(nsh),nuHlocal(nsh)
      dimension xsilocal(nsh),cos2local(nsh),explocal(nsh)


* Parametres de la tachocline
*     
*     dtacho = 2*h, h = tachocline thickness
      dtacho = 0.1d0*r(nmod)
c      dtacho = 0.2d0*r(nmod)
*     Brunt-Vaisala frequency
      NBV = 25.d0   ! muHz
      Q4 = -1.707d-2
      mu4 = 4.933d0
      omegaBZC = 0.415d0
      coefNuH = (r(ndt)/dtacho)**4
      coefNuH = 4.d0*coefNuH*(2.d0*omegaBZC/NBV)**2       ! see eq. 11
      deltalocal = (1.d0/180.d0)*(1.d0/4.d0)*(8.d0/3.d0)**2
      deltalocal = deltalocal*mu4**6
      deltalocal = deltalocal*Q4*Q4       

*     ndt = novlim(klenv,3)-1 dans diffusion.f

* Recherche de la couche limite inferieure de la tachocline
      rliminftacho = r(ndt)-dtacho
      do i = ndt,ndb,-1
         iliminftacho = i
         if (r(i).lt.rliminftacho) goto 10
      enddo
 10   continue

      do i = ndt,iliminftacho,-1
          xKtlocal(i) = khi(i)/cp(i)/ro(i)
          nuHlocal(i) = coefNuH*xKtlocal(i)
          xsilocal(i) = mu4*(r(ndt)-r(i))/dtacho
          explocal(i) = exp(-2.d0*xsilocal(i))
          cos2local(i) = cos(xsilocal(i))*cos(xsilocal(i))
          coefDtacho(i) = deltalocal*nuHlocal(i)
          coefDtacho(i) = coefDtacho(i)*(dtacho/r(ndt))**2
          coefDtacho(i) = coefDtacho(i)*explocal(i)*cos2local(i)  ! coefDtacho = Dt eq. 15
          write (29,100) coefDtacho(i)
      enddo

 100  format (1x,1pe11.4)

      return
      end


      SUBROUTINE ZONEDELTAXRGB (ndt,ndb,ndbmel)

      implicit none

      include 'evolpar.star'

      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.transp'

      integer ndt,ndb,ndbmel
      integer j

      double precision xlim

      xlim = xsp(nmod,2)-grad_crit
      do j = ndt,ndb,-1
         if (xsp(j,2).le.xlim) then
            ndbmel = j
            goto 11
         endif
      enddo

 11   return
      end



      SUBROUTINE ZONEMUCARBON (ndb,ibottom)

      implicit none

      include 'evolpar.star'

      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.transp'

      integer ndb,ibottom
      integer i

      double precision douzetreize

      logical change

      dimension douzetreize(nsh)

      ibottom = 1
c     douzetreize_crit = 10.

      do i = 1,nmod
         if (xsp(i,15).gt.0.d0) douzetreize(i) = xsp(i,14)/xsp(i,15)
         if (xsp(i,15).eq.0.d0) douzetreize(i) = 0.d0
      enddo

      do while (i.gt.ndb)

         change = (douzetreize(i).lt.4.d0.and.douzetreize(i).gt.3.d0
     &        .and.abs(douzetreize(i)-douzetreize(i-1)).lt.
     &        5.d-2)
c        change = (douzetreize(i).ge.10.and.douzetreize(i-1).lt.10)

         if (change) then
            ibottom = i
            goto 10
         else 
            ibottom = i
         endif
         i = i-1
      enddo

 10   return
      end



      SUBROUTINE ZONEMELANGE (ndt,ndb)

      implicit none

      include 'evolpar.star'

      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'

      integer ndt,ndb
      integer ibottom,idown
      integer i,j,k

      double precision deltam

      common /coeff_diff/ ibottom,idown

      dimension deltam(nsh)

      ibottom = 1
      idown = 1

      do j = ndt,ndb,-1
         if (xsp(j,2).le.1.d-4) then
            idown = j 
            goto 10
         endif
      enddo

 10   do i = ndt,ndb,-1
         deltam(i) = (mr(i)-mr(idown))/(mr(ndt)-mr(idown))
      enddo

      do k = ndt,ndb,-1
         if (deltam(k).le.5.d-2) then
            ibottom = k
            goto 20
         endif
      enddo

 20   return
      end

**********************************************************************

      SUBROUTINE TRIDIAG (A,B,C,r,u,ldb)

**********************************************************************
c     inverse tridiagonal matrix : numerical recipes
c     A,B,C lower, central and upper bands, given r solves for u

      implicit none

      integer j,ldb

      double precision a,b,c,r,u,gam,bet

      dimension gam(ldb),a(*),b(*),c(*),r(*),u(*)

      if (b(1).eq.0.d0) stop 'tridiag : b(1) = 0.d0'
      bet=b(1)
      u(1)=r(1)/bet
      do j=2,ldb
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if (bet.eq.0.d0) then
           print *,j,gam(j),c(j-1),bet,b(j),a(j)
           stop 'tridiag : bet'
        else
           u(j)=(r(j)-a(j)*u(j-1))/bet
        endif
      enddo
      do j=ldb-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
      enddo

      return
      end

