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
