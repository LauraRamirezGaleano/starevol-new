C-----------------------------------------------------------------------
      subroutine corrwind (teffcorr,reffcorr)
C-----------------------------------------------------------------------
C Sous-routine calculant la correction a la temperature effective due
C au vent stellaire tenant compte de la diffusion par electrons libres
C et l'influence des raies. (Theorie CAK amelioree)
C  (D. Schaerer, Fin juillet 1990)               
C                                                 
C  Parametre change: TEFFCORR = log10 de la valeur corrigee de la Teff
C                              = log10(teff) dans le code.
C
C les rayons sont exrpimes en centimetres (reff, r0)
C la perte de masse est exprimee en g/sec (mlr = dms*msun/sec)
C les vitesses sont exprimees en cm/s (vinf)
C
C Routine adaptee du code de Geneve-Toylouse pour STAREVOL
C
C Les variables PATMOS et TATMOS representent le log de la pression
C resp. temperature a tau=2/3 de l'atmosphere du programme d'evolution.
C On s'en sert pour calculer la concentration d'electrons. De cette
C maniere on peut tenir compte du changement pour les etoiles WR.
C
C Premiere approche: on suppose se=cst=0.22 (Langer)
C                    ceci nous evite de calculer les nouveaux mu
C
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
C-----------------------------------------------------------------------

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.surf'
      include 'evolcom.var'

      integer i

      double precision pi4
      double precision telec1,telec0,tline,rho,gradv1,gradv0
      double precision se,mlr,vinf,bet,r0,ck,del,alp,vth,ex
      double precision step
      double precision teffcorr
      double precision tintlines,tintelec,tautot,tmcak,xold,tauelecold,
     &     tauold
      double precision xmin,xmax,xlen
      double precision x,gr,qu,gg,xh,xtwo,cf,rh
      double precision reffcorr

      common /paramwind /se,mlr,vinf,bet,r0,ck,del,alp,vth,ex

      external telec1,telec0,tline,rho,gradv1,gradv0

C----  Constantes:  ----------------------------------------------------
C msun: facteur par lequel on multiplie dms pour avoir une
C       expression de la perte de masse en unites 'Masse solaire par an'
C rsun:  Rayon solaire en 'cm'
C pi4:   4 * pi
      parameter (pi4=12.566371d0)

C  Pour pouvoir changer les parametres du vent et des raies en cours de
C      OPEN(23,FILE='WINDCORRIN',STATUS='OLD')    
C      READ(23,*) VINF                            
C      READ(23,*) BET                             
C      READ(23,*) CK                              
C      READ(23,*) ALP                             
C      READ(23,*) DEL                             
C      READ(23,*) R0                              
C      READ(23,*) XMAX                            
C      READ(23,*) STEP                            
C      CLOSE(23)                                  

C  Si l'on prend des valeurs fixes, p.e.:         
      vinf = 2.d3
      bet =  2.d0
      ck =   0.124d0
      alp =  0.64d0
      del =  0.07d0
      r0 =   1.0d0
      xmax = 1.d3
      step = 6.d2

C-----Conversions en bonnes unites:----------------------------------

c      reff = surfr
      mlr = dms*msun/sec
      vinf = vinf*1.d5
      r0 = r0*reff


C----------Initialisations:-------------------------------------------

C  Langer posait SE=0.22 dans la correction appliquee jusqu' a present.
C  Ici cependant on recalcule SE:                 
c      call ionpart(patmos,tatmos)
c      se=0.398*(1./vmy-1./vmol)
      se = 0.22d0

C  On met TEFFCORR = 0 au debut. Si la routine donne ces valeurs
C  ceci signalerait que tau<2/3 jusqu'a X=1.0001 * REFF.
C  Dans ce cas les erreurs deviendrait trop grandes avec ce champ de
C  vitesse adopte.                                
      teffcorr = 0.d0

      tintlines = 0.d0
      tintelec = 0.d0
      tautot = 0.d0
      tmcak = 0.d0
      xold = 0.d0
      tauelecold = 0.d0
      tauold = 0.d0

      ex = (2.d0*alp*bet-alp-bet+1.d0)
      vth = sqrt(1.649959d8*teff)

C------------------Methode MCAK avec integration:-------------------

C  Avec le choix de XLEN=XMAX+4. on s'approchera au maximum jusqu'a
C  1.0001 * REFF.                                 
      xmin = -4.d0
      xmax = log10(xmax)
      xlen = xmax+abs(xmin)

      do i = int(abs(step)),0,-1
        x = reff*(1.d0 +1.d1**( i*xlen/abs(step) + xmin))

        if (bet.ge.1.d0) then
          gr = gradv1(x)
        else
          gr = gradv0(x)
        endif
        tmcak = se*vth*rho(x)/gr

C Partie dependante de delta:                     
        if (bet.ge.2.d0) then
          qu = 22.5d0**del-1.d0
        elseif (bet.ge.1.d0) then
          qu = 7.5d0**del-1.d0
        elseif (bet.ge.0.7d0) then
          qu = 4.d0**del-1.d0
        elseif (bet.ge.0.5d0) then
          qu = 2.5d0**del-1.d0
        else
          qu = 1.18d0**del-1.d0
        endif
        gg = ((mlr/(pi4*reff*reff*vinf))**del)*((1.02436d+13)**del)*
     1    2.d0**del*(qu*(reff/x)**2.d0+1.d0)
C Partie du Correction factor:                    
        xh = ((x/reff)-1.d0)/bet
        xtwo = (x/reff)*(x/reff)
        cf = (xtwo/((alp+1.d0)*(1.d0-xh)))*(1.d0-(1.d0+
     &       (xh-1.d0)/xtwo)**(alp+1.d0))

        tmcak = (ck/(tmcak**alp)) * gg * cf
        rh = mlr/(pi4*x*x*vinf*(1.d0-r0/x)**bet)

        if (xold.ne.0.d0) then
          tintlines = tintlines+(xold-x)*rh*se*tmcak
C Si l'on utilise un champ de vitesse different il faudra aussi integrer
C la prof.opt. des electrons numeriquement (prochaine ligne).
C Ici on utilise la version analytique...         
C          TINTELEC=TINTELEC+(XOLD-X)*RH*SE       

          if (bet.ge.1.d0) then
            tintelec = telec1(x)
          else
            tintelec = telec0(x)
          endif

          tautot = tintlines+tintelec
        endif

C Quand l'integration a depassee tau=2/3 on trouve le REFFCORR par
C interpolation lineaire.
        if (tautot.ge.pw23.and.teffcorr.eq.0.d0) then
          reffcorr = xold-(xold-x)*(pw23-tauold)/(tautot-tauold)
          teffcorr = log10(teff) + 0.5d0*log10(reff/reffcorr)
       endif

C Comme comparaison on fait la meme chose pour la profondeur optique des
C electrons (sans raies) et puis on peut TERMINER l'integration.
        if (tintelec.ge.pw23) then
          reffcorr = xold-(xold-x)*(pw23-tauelecold)/
     &          (tintelec-tauelecold)

          goto 99
        endif

        tauelecold = tintelec
        tauold = tautot
        xold = x

      end do


C----------------fin MCAK-------------------------------------------

   99 continue

      end

*--------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION telec1(r)
*--------------------------------------------------------------------
* Fonction definissant la profondeur optique due aux electrons
* dans le cas ou beta > 1
*--------------------------------------------------------------------

      implicit none

      common /paramwind /se,mlr,vinf,bet,r0,ck,del,alp,vth,ex

      double precision r,se,mlr,vinf,bet,r0,ck,del,alp,vth,ex
      double precision pi4

      parameter (pi4=12.566371d0)

      telec1 = se*mlr/(pi4*vinf*r0*(1.d0-bet))*
     &     (1.d0-1.d0/(1.d0-r0/r)**(bet-1.d0))

      return
      end

*--------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION telec0(r)
*--------------------------------------------------------------------
* Fonction definissant la profondeur optique due aux electrons
* dans le cas ou beta < 1
*--------------------------------------------------------------------

      implicit none

      common /paramwind /se,mlr,vinf,bet,r0,ck,del,alp,vth,ex

      double precision r,se,mlr,vinf,bet,r0,ck,del,alp,vth,ex
      double precision pi4

      parameter (pi4=12.566371d0)

      telec0 = se*mlr/(pi4*vinf*r0*(1.d0-bet))*
     &     (1.d0-(1.d0-r0/r)**(1.d0-bet))

      return
      end

*--------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION tline (r)
*--------------------------------------------------------------------
* Fonction definissant la profondeur optique due aux electronsaux raies
*--------------------------------------------------------------------

      implicit none

      common /paramwind /se,mlr,vinf,bet,r0,ck,del,alp,vth,ex

      double precision r,se,mlr,vinf,bet,r0,ck,del,alp,vth,ex
      double precision pi4

      parameter (pi4=12.566371d0)

      tline = (ck*(se**(1.d0-alp))*((mlr/pi4)**(1.d0-alp))*
     &     (vinf**(2.d0*alp-1.d0))*(bet**alp)*(1.d0-(1.d0-r0/r)**ex))
     &     /((vth**alp)*(r0**(1.d0-alp))*ex)

      return
      end



*-------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION rho(r)
*-------------------------------------------------------------------
* Densite dans le vent
*-------------------------------------------------------------------
      implicit none

      common /paramwind /se,mlr,vinf,bet,r0,ck,del,alp,vth,ex

      double precision r,se,mlr,vinf,bet,r0,ck,del,alp,vth,ex
      double precision pi4

      parameter (pi4=12.566371d0)

      rho = mlr/(pi4*r*r)/(vinf*(1.d0-r0/r)**bet)

      return
      end
*-------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION gradv1(r)
*-------------------------------------------------------------------
* Gradient dans le vent pour beta > 1
*-------------------------------------------------------------------
      implicit none

      common /paramwind /se,mlr,vinf,bet,r0,ck,del,alp,vth,ex

      double precision r,se,mlr,vinf,bet,r0,ck,del,alp,vth,ex

      gradv1 = vinf*bet*((1.d0-r0/r)**(bet-1.d0))*r0/(r*r)

      return
      end

*-------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION gradv0(r)
*-------------------------------------------------------------------
* Gradient dans le vent pour beta < 1
*-------------------------------------------------------------------
      implicit none

      common /paramwind /se,mlr,vinf,bet,r0,ck,del,alp,vth,ex

      double precision r,se,mlr,vinf,bet,r0,ck,del,alp,vth,ex

      gradv0 = vinf*bet*((1.d0-r0/r)**(bet-1.d0))*r0/(r*r)

      return
      end

