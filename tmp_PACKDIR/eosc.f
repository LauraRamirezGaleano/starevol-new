
************************************************************************
*                    EQUATION OF STATE COMPUTATION                     *
*                           Emmanuel DUFOUR                            *
*                                                                      *
*                             References:                              *
*                        MNRAS 274,964 (1995)                          *
*                          A&A 23,325 (1973)                           *
************************************************************************
*                       Version 2.90: Jun 2006                         *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************
* All the subroutine parameters are written explicitly except the four *
*  "include" described below:                                          *
*  eoscom.math     :  number of calculated derivatives                 *
*  eoscom.ion      :  maximum characteristics of ionic component of the*
*                     plasma (limited by some datas: energy levels,...)*
*  eoscom.phys     :  physical characteristics (E and g0)              *
*  eoscom.cons     :  local physical and mathematical constants        *
*  eoscom.par      :  somes parameters of EOS                          *
************************************************************************
*                                                                      *
*                  Characteristics of this eos :                       *
*                                                                      *
*  NBZ is the maximum number of a set of integers {1...nbz}. These     *
*   numbers are the index of the nbz chemical species composing the    *
*   plasma: N1..Ni..Nnbz.                                              *
*  The first Z found in znuc at index nspmin is refered to N1 and must *
*  be Hydrogen. The last Z found in znuc at index nspmax is refered to *
*  Nnbz and must be heavy.                                             *
*  For each Ni, you choose to compute ionization or not: ionized(nbz)  *
*  If totioniz or not , the others are considered always totaly        *
*  ionized or always neutral.                                          *
*  The caracteristics of each species are Z(nbz), X(nbz), Y(nbz) plus  *
*  lnmH, lnmZ(nbz-2), lnmheavy and are deduced from xsp, ysp, anuc,    *
*  and znuc.                                                           *
*  To check these datas, assign true to checkcomposition in the        *
*  parameter card.                                                     *
*  For computed ionization species you need partition function and     *
*  energy levels. The partition functions are often approximated to be *
*  the fundamental energy level degeneracy, except for H2 whose        *
*  rotational and vibrationnal partition function are always taken into*
*  account. If boolean compZspecie is true eos will use a more complete*
*  partition function for specie too. Energy level degeneracy and      *
*  energy levels are stored in lng0H(4), lng0Z(nbzmax-2), lng0heavy,   *
*  enerlevelH(4) enerlevelZ(nbionZmax). The number of datas wich have  *
*  been stored allow you to compute ionization for elements from He to *
*  Ca. If not totioniz, for not computed ionization species, you need  *
*  lng0 too, so you can not choose not totioniz if Z>Ca exist.         *
*                                                                      *
*                 Furthermore, It's quite easy to add                  *
*                                                                      *
* - Partition functions and their derivatives for all elements in ZZ.  *
* - Others non perfect corrections in nonperfect (for the moment       *
*   Coulomb shielding and ionization pressure are taken into account). *
************************************************************************
*
      SUBROUTINE eos(flag,error)
*
************************************************************************
*
      IMPLICIT NONE
*
************************************************************************
*                           DECLARATIONS                               *
************************************************************************
*
*     ----------Inputs----------
*
*     number of isotopes (nsp) and chemical elements computed (nsp)
      include 'evolpar.star'
*     Z of the differents isotopes (znuc)
      include 'evolcom.nuc'
*     Independant variables (lnf and T)
      include 'evolcom.var'
*     flag: if flag=0 then eos initialized its constants
      INTEGER flag
      include 'eoscom.ion'
      include 'eoscom.par'
      include 'eoscom.phys'
      include 'eoscom.cons'
*
*     ----------Outputs----------
*
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.eos'
      include 'evolcom.grad'
*
*     ----------Intermediate variables----------
*
      INTEGER error
************************************************************************
*                               BEGIN                                  *
************************************************************************
*
      IF (flag.eq.0) THEN
         CALL initpar
         CALL initcons
         nonideal = (compcoulshiel.OR.compionpres)
         CALL initnbion
c         IF (.NOT.compavmass) THEN
c            CALL initlnm
c         END IF
*        lng0(56Fe)=ln(25)
         lng0heavy = 3.2188758d0
      ELSE
         CALL eos_sod(error)
         if (error.gt.0) return
      END IF
      END

************************************************************************
*                   EQUATION OF STATE COMPUTATION                      *
*           with ANALYTICAL Sedond Order Derivatives (SOD) computation *
*                                                                      *
* corlnf and corlnT are the corrections to apply to lnf and T in order *
* to get the derivatives                                               *
************************************************************************

      SUBROUTINE eos_sod(error)
*     ------------------
      IMPLICIT NONE
*     ----------Inputs----------
*
*     number of isotopes (nsp) and chemical elements computed (nsp)
      include 'evolpar.star'
*     maximum number of shell (nsh)
      include 'evolcom.teq'
*     physical constants
      include 'evolcom.cons'
*     Z of the differents isotopes (znuc)
      include 'evolcom.nuc'
*     molar abundances of the different isotopes (ysp)
      include 'evolcom.spec'
*     Independant variables (lnf and T)
      include 'evolcom.var'
*     Considering atomic diffusion or not       ! TD Nov. 2018
      include 'evolcom.diff'
      include 'evolcom.mod'
*
*     ----------Outputs----------
*
      include 'evolcom.therm'
      include 'evolcom.ion'
      include 'evolcom.therm2'
      include 'evolcom.eos'
      include 'evolcom.grad'
*
*     ---------- Mathematics and informatics ----------
*
      INTEGER i,jsh
      INTEGER j,zz                     ! Ionisation partielle TD Oct.2018
*     Chemical elements to be ionized
      include 'eoscom.ion'
*     Parameters
      include 'eoscom.par'
*     Physical parameters
      include 'eoscom.phys'
*     Local constants (CGS!)
      include 'eoscom.cons'
*     Number of derivatives
      include 'eoscom.math'
*     Truncated exp function
      DOUBLE PRECISION truncexpo,truncinv,truncdiff
*     Intermediate variables
      DOUBLE PRECISION term1,term2
*
      INTEGER error
*     ---------- Reformulation of input variables----------
*
*     degeneracy parameter
      DOUBLE PRECISION etash(deriv2),degen
*     Temperature in k/(mc2) units, kT, beta
      DOUBLE PRECISION Tnodim,RT,betaT,betaev,lnbetaev,Tinv,lnTsh,Tsh
*     square root of 1+f, f/(2*(f+1))
      DOUBLE PRECISION rootf,halff1,rootfinv,fadd1,fadd1inv,lnfsh
*     Total number of electrons (free or bounded)
      DOUBLE PRECISION addy,addYZ
*     Charge averages
      DOUBLE PRECISION zavY,zavYZ
*     Abundances for Z fixed
      DOUBLE PRECISION Y(nbz),X(nbz)
*     Current shell number
      INTEGER shell
*
*     ---------- For electronic component ----------
*
*     variables of Fermi-Dirac polynomial approximation
      DOUBLE PRECISION fpoly
*     adimensionnal thermodynamic charateristics of eletronic component
*     and first and second derivatives versus f,T
      DOUBLE PRECISION ine(deriv2),ise(deriv2),ipe(deriv2),iue(deriv2)
*     dimensionnal thermodynamic charateristics of eletronic component
      DOUBLE PRECISION pelec(deriv2),selec(deriv2),uelec(deriv2)
*
*     ---------- For non perfect gas corrections ---------
*
*     Variables for the non-perfect gas correction functions
      DOUBLE PRECISION betanodim,amune,amune0
*     Total non-perfect gas correction function and derivatives
      DOUBLE PRECISION cortot(deriv3)
*     Idem without pressure ionization gas correction
      DOUBLE PRECISION cornpi(deriv3)
*     The same for completely ionized gas (only for pressure ionization)
      DOUBLE PRECISION corpi0(deriv3)
*     Non perfect gas corrections for thermodynamic quantities p,s and u
      DOUBLE PRECISION dpe(deriv2),dse(deriv2),due(deriv2)
*     Effective degeneracy taking into account non perfect corrections
      DOUBLE PRECISION etaeff(deriv2),etanpi(deriv2)
*     Intermediate variables
      DOUBLE PRECISION dmut1,dmut2,dmut1x,dmut1y,dmut2x
      DOUBLE PRECISION dst1,dst2,dst1x,dst2y,dst2x
      INTEGER num
*
*     ---------- For ionization ----------
*
*     Partition function and derivatives
      DOUBLE PRECISION lnparfunH(nbionHmax,deriv3)
      DOUBLE PRECISION lnparfunz(nbionZmax,deriv3)
*     Abundances of ions and electrons and derivatives
      DOUBLE PRECISION abundH(nbionHmax,deriv2),abundZ(nbionZmax,deriv2)
      DOUBLE PRECISION abundelec(deriv2),abundelecinv
*     number of free eletrons and derivatives
      DOUBLE PRECISION heavelec(deriv2)
*     Abundances of electrons for completly ionized plasme
      DOUBLE PRECISION abundelec0(deriv2)
*     Which solution of the quadratic equation f or H has been choosen
      LOGICAL solutionH
*     Intermediate variables from H ionization in order to compute kimu
      DOUBLE PRECISION sumh,diffh,ah,bh,dheavelecdY(nbz),sdelh
*
*     ---------- Radiative and gazeous component ----------
*
*     Gaseous component for thermodynamic quantities p,s and u
      DOUBLE PRECISION pgaz(deriv2),sgaz(deriv2),ugaz(deriv2)
*     Sum of the volumic abundances
      DOUBLE PRECISION sumNi(deriv2),sumNinv
*     Raditive component for thermodynamic quantities p,s and u
      DOUBLE PRECISION pradi(deriv2),sradi(deriv2),uradi(deriv2)
*
*     ---------- Computation of thermodynamic quantities ----------
*
*     Pressure derivatives in the set of variables (rho,T)
*     and derivatives versus f and T.
      DOUBLE PRECISION zTsh,zrhosh,dzTdfsh,dzTdTsh,dzrodTsh,dzrodfsh
      DOUBLE PRECISION zTshinv,zrhoshinv
*     Pressure derivatives in the set of variables (f,T)
      DOUBLE PRECISION zaidf,zaidT,zaidff,zaidTT,zaidfT,zaidfinv
*     dlnmu/dlnP at a given shell
      DOUBLE PRECISION khimush,ksirho,ksiT,khirho,khiT
*     p,s,u at a given shell - srsh is the reduced entropy (s*mu/R)
      DOUBLE PRECISION ssh(deriv2),srsh,psh(deriv2),ush(deriv2)
      DOUBLE PRECISION pshinv,pfshinv,sfshinv
*     computation of derivatives versus rho with approximations
      DOUBLE PRECISION proro,pTT,proT,zrodro,zrodT,zTdro,zTdT
      DOUBLE PRECISION cpro,cpT,etro,ett,qro,qT,cpadd
*     ratio Pe/Pe(non degenerated) at a given shell
      DOUBLE PRECISION degpesh
*     ratio Pgaz/Ptotal at a given shell
      DOUBLE PRECISION betapsh,betanicorsh
*     mean molecular weight and derivatives at a given shell
      DOUBLE PRECISION mush,mushe,mushi,mushf,mushT
      DOUBLE PRECISION mushinv,mushiinv,musheinv
*     specific heat and derivatives at a given shell
      DOUBLE PRECISION cpsh,cpshT,cpshf,cvsh
*     adiabatic gradient and sound velocity at a given shell
      DOUBLE PRECISION adadsh,adadshf,adadshT,vssh
*     factor q and derivatives at a given shell, gamma1
      DOUBLE PRECISION Qsh,Qfsh,QTsh,gamma1sh
*     quantities
      DOUBLE PRECISION rhpsish,ratnehsh
      DOUBLE PRECISION rationsh(nsp,zionmax),ration(nsh,nis,zionmax) !ration(nsh,nis,18)     ! Ionisation partielle TD Oct.2018
      DOUBLE PRECISION zmeansh(nsp)
*     Intermediate variables
      DOUBLE PRECISION int1,int1f,int1T,int2,ratio1,ratio2,int1inv
*     rho at a given shell
      DOUBLE PRECISION rho(deriv2),rhofinv,rhoinv,drhodT,drhodlnf
*

      double precision etak,pek
      common /thermkap/ etak(nsh),pek(nsh)
      common /abundelec/ abundelec
c      common /ioni/ ration(nsh,nis,18)!, zmean(nsh,nsp)

      logical totioniz0
************************************************************************
*                         I - MAIN LOOP                                *
************************************************************************
*

c$$$      CALL zmoy(xsp)  ! test TD Nov. 2018
      totioniz0 = totioniz
c      print *, 'diffusion atomique',microdiffus,'totioniz2',totioniz
      DO shell=1,nmod
c         print *,'shell',shell
         ionized(2) = ionizHe
c..   if not atomic diffusion, not considering partial ionization of not abundant element  ! TD Nov.2018
         if (model.le.2) microdiffus=.false.
   
c..   if T > 5.d6 or rho > 1d.4, impose complete ionization for CNO,Ne
c     if (T(shell).gt.5.d6.or.lnf(shell).gt.1.d2) then
         if(T(shell).gt.5.d6.or.lnf(shell).gt.1.d2.or..not.microdiffus)            ! modif. TD Nov.2018
     $        then
c$$$            ionized(6) = .false.
c$$$            ionized(7) = .false.
c$$$            ionized(8) = .false.
c$$$            ionized(10) = .false.
            do i=3,nbz-1        
               ionized(i) = .false.
            enddo
            totioniz = .true.
         else
c ion part
            if (nphase.gt.1) then
               do i = 2,17
                  ionized(i) = .true.
               enddo
               ionized(14) = .false.
            endif
c ion part            
            ionized(6) = ionizC
            ionized(7) = ionizN
            ionized(8) = ionizO
            ionized(10) = ionizNe
            totioniz = totioniz0
         endif
c         stop

c..   whatever happens, H and He fully ionized for T > 1.d7 or rho > 1.d4
c          if (ionizHe) then
c             if (T(shell).lt.1.d7.and.lnf(shell).lt.1.d2) then
c                ionized(1) = .true.
c                ionized(2) = .true.
c             else
c                ionized(1) = .false.
c                ionized(2) = .false.
c             endif
c          else
c             if (T(shell).lt.tmaxioH.and.lnf(shell).lt.1.d2) then
c                ionized(1) = .true.
c             else
c                ionized(1) = .false.
c             endif
c             if (T(shell).lt.tmaxioHe.and.lnf(shell).lt.1.d2) then
c                ionized(2) = .true.
c             else
c                ionized(2) = .false.
c             endif
c          endif

c..   if X+Y < 0.01 or T > 4.d8 or rho > 1.d4
         if (ionizHe) then
            if (T(shell).gt.4.d8.or.xsp(shell,ih1)+xsp(shell,ihe4).lt.
     &           1.d-2.or.lnf(shell).gt.1.d2) then
               ionized(1) = .false.
               ionized(2) = .false.
            else
               ionized(1) = .true.
               ionized(2) = .true.
            endif
         else
            if (T(shell).lt.tmaxioH.and.lnf(shell).lt.1.d2) then
               ionized(1) = .true.
            else
               ionized(1) = .false.
            endif
            if (T(shell).lt.tmaxioHe.and.lnf(shell).lt.1.d2) then
               ionized(2) = .true.
            else
               ionized(2) = .false.
            endif
         endif

c..   if H and He fully ionized, do not compute pressure ionization
         if (.not.ionized(1).and..not.ionized(2)) then
            compionpres = .false.
         else
            compionpres = addprio
         endif

         lnfsh = lnf(shell) - lnf0
         Tsh = T(shell)
         lnTsh = DLOG(Tsh)
         fpoly = truncexpo(lnfsh)
         fadd1 = 1.0d0 + fpoly
         fadd1inv = 1.0d0 / fadd1
         rootf = DSQRT(fadd1)
         rootfinv = 1.0d0 / rootf
         halff1 = 0.5d0 * fpoly * fadd1inv
         Tinv = 1.0d0 / Tsh
         RT = rk * Tsh
         betaT = 1.0d0 / (boltz*Tsh)
         betaev = ergev * betaT
         lnbetaev = DLOG(betaev)
         betanodim = khiH * betaev
         Tnodim = mc2inv * boltz * Tsh
*
************************************************************************
*               II - ADIMENSIONNAL ELECTRONIC COMPONENT COMPUTATION    *
************************************************************************
*
         etash(1) = degen(lnfsh,rootf)
         etash(2) = 0.0d0
         etash(3) = rootf
         etash(4) = 0.0d0
         etash(5) = halff1
         etash(6) = 0.0d0
         CALL electronic_sod(fpoly,Tnodim,etash,rootfinv,
     &fadd1inv,ine,ipe,ise,iue)
*
************************************************************************
*        III - GAZEOUS COMPONENT COMPUTATION                           *
************************************************************************
*
*        ---------- Computation of <Z>, sum(Y)... ----------
*
         CALL chemY(shell,xsp,ysp,X,Y,addY,addYZ,zavY,zavYZ,zmeansh)
c         IF (compavmass) THEN
            CALL initlnmav(shell,X,Y)
c         END IF
*
*        ---------- computation of effective eta taking into -----------
*        ---------- account non perfect effects ------------------------
*
         amune = ine(1) * m8pl3
         IF (nonideal) THEN
            CALL nonperfect_sod(fpoly,rootfinv,amune,betanodim,zavY,
     & zavYZ,cortot,dmut1,dmut2,dmut1x,dmut1y,dmut2x,dst1,dst2,dst1x,
     & dst2y,dst2x,halff1,ine,fadd1inv,cornpi)
         ELSE
            DO i=1,deriv3
               cortot(i) = 0.0d0
               cornpi(i) = 0.0d0
               corpi0(i) = 0.0d0
            END DO
         END IF
         CALL neweta_sod(cortot,ine,dmut1,dmut2,dmut1x,dmut1y,dmut2x,
     &cornpi,etash,etaeff,etanpi)
*
*        ---------- ionization computation -----------------------------
*
         CALL ZH_sod(Tsh,lnTsh,betaev,lnparfunH)
         CALL ZZ_sod(Tsh,lnTsh,lnparfunZ)
         CALL sahaZ_sod(shell,betaev,Y,abundZ,lnparfunZ,etaeff,error)
         if (error.gt.0) then
            jsh = max(shell,3)
            write (nout,*) 'shell =',jsh,',T (shell-2 -->shell+2) =',
     &           T(jsh-2),T(jsh-1),T(jsh),T(jsh+1),T(jsh+2)
            RETURN
         endif
         CALL freelectronZ_sod(Y,abundZ,heavelec,zmeansh,dheavelecdY)
         i=shell
         CALL sahaH_sod(i,Y,betaev,lnbetaev,heavelec,abundelecinv,
     &abundH,abundelec,lnparfunH,etaeff,ine,solutionH,zmeansh,Tsh,
     &etanpi,sumh,diffh,ah,bh,sdelh,error)
         if (error.gt.0) RETURN
*
*        ---------- rho computation ------------------------------------
*
         rho(1) = amune * abundelecinv
CC         if (rho(1).gt.1.d12.or.rho(1).lt.1.d-20) then
         if (rho(1).gt.1.d15.or.rho(1).lt.1.d-30) then
            write (nout,"('eos #',i4,' rho =',1pe11.3,', T =',1pe11.3)")
     &           shell,rho(1),Tsh
            error = 10
            return
         endif
         rhoinv = truncinv(rho(1))
         ratio1 = abundelec(2) * abundelecinv
         ratio2 = abundelec(3) * abundelecinv
         rho(2) = ine(2) - ratio1
         rho(3) = ine(3) - ratio2
         rhofinv = truncinv(rho(3))
*        d/dT, d/dlnf
         drhodT = rho(2) * rho(1) * Tinv
         drhodlnf = rho(3) * rho(1)
*        d(dln)/(dln*dln)
         rho(4) = ine(4) + ratio1 * ratio1 -
     &abundelec(4) * abundelecinv
         rho(5) = ine(5) + ratio2 * ratio2 -
     &abundelec(5) * abundelecinv
         rho(6) = ine(6) + ratio1 * ratio2 -
     &abundelec(6) * abundelecinv
*
*
*        ---------- p,s,u computation ----------------------------------
*
         sumNi(1) = addY - abundH(4,1)
         sumNinv = 1.0d0 / sumNi(1)
         sumNi(2) = -abundH(4,2)
         sumNi(3) = -abundH(4,3)
         sumNi(4) = -abundH(4,4)
         sumNi(5) = -abundH(4,5)
         sumNi(6) = -abundH(4,6)
         CALL gaz_sod(Y,RT,betaev,lnbetaev,sumNi,lnparfunH,lnparfunZ,
     &abundH,abundZ,pgaz,sgaz,ugaz,rho,sumNinv)
*
************************************************************************
*          III - NON PERFECT COMPONENT COMPUTATION                     *
************************************************************************
*
         DO i=1,deriv2
            dpe(i) = 0.0d0
            dse(i) = 0.0d0
            due(i) = 0.0d0
         END DO
         IF (nonideal) THEN
            num = 1
            CALL dpecalc_sod(num,cortot,amune,ine(3),ine(2),ine(5),
     &ine(4),ine(6),RT,dmut2,dmut1x,dmut2x,dpe)
            CALL dsecalc_sod(num,cortot,ine(3),ine(2),ine(4),
     &ine(5),ine(6),dst1,dst2,dst1x,dst2y,dst2x,abundelec,dse)
            CALL duecalc_sod(num,cortot,RT,abundelec,ine(3),ine(2),due)
            IF (compionpres) THEN
               abundelec0(1) = addYZ
               abundelec0(2) = 0.0d0
               abundelec0(3) = 0.0d0
               abundelec0(4) = 0.0d0
               abundelec0(5) = 0.0d0
               abundelec0(6) = 0.0d0
               amune0 = rho(1) * abundelec0(1)
               num = -1
               CALL nonperfect0_sod(amune0,betanodim,corpi0,dmut1,
     &dmut2,dmut1x,dmut1y,dmut2x,dst1,dst2,dst1x,dst2y,dst2x)
               CALL dpecalc_sod(num,corpi0,amune0,rho(3),rho(2),rho(5),
     &rho(4),rho(6),RT,dmut2,dmut1x,dmut2x,dpe)
               CALL dsecalc_sod(num,corpi0,rho(3),rho(2),rho(4),
     &rho(5),rho(6),dst1,dst2,dst1x,dst2y,dst2x,abundelec0,dse)
               CALL duecalc_sod(num,corpi0,RT,abundelec0,rho(3),rho(2),
     &due)
            END IF
         END IF
*
************************************************************************
*           IV - RADIATIVE AND ELECTRONIC COMPONENT COMPUTATION        *
************************************************************************
*
         CALL radiativ_sod(Tsh,pradi,sradi,uradi,rho,rhoinv)
         CALL dimelec_sod(abundelec,RT,ise,ipe,iue,pelec,selec,uelec)
*
************************************************************************
*    V  - BASIC THERMODYNAMICS QUANTITIES                              *
*         (suffix "sh" means "at a fixed shell")                       *
************************************************************************
*
         DO i=1,deriv2
            psh(i) = pradi(i) + pgaz(i) + pelec(i) + dpe(i) + zero
            ssh(i) = sradi(i) + sgaz(i) + selec(i) + dse(i) + zero
            ush(i) = uradi(i) + ugaz(i) + uelec(i) + due(i) + zero
         END DO
         if (psh(1).lt.0.d0) then
            DO i=1,deriv2
               dpe (i) = 0.d0
               dse (i) = 0.d0
               due (i) = 0.d0
               psh(i) = pradi(i) + pgaz(i) + pelec(i) + zero
               ssh(i) = sradi(i) + sgaz(i) + selec(i) + zero
               ush(i) = uradi(i) + ugaz(i) + uelec(i) + zero
            END DO
         endif

         pshinv = 1.0d0 / psh(1)
         pfshinv = 1.0d0 / psh(3)
         sfshinv = 1.0d0 / ssh(3)
*
************************************************************************
*           V - OTHER THERMODYNAMICS QUANTITIES                        *
************************************************************************
*
*        ---------- zf,zT,zro,zt ------------------------------
*
         zaidT = psh(2) * pshinv
         zaidf = psh(3) * pshinv
         CALL transdftodrho(.true.,zaidf,zaidT,zrhosh,zTsh,rho(2),
     &rho(3),rhofinv)
         zrhoshinv = truncinv(zrhosh)
         zTshinv = truncinv(zTsh)
*
*        ---------- mean molecular weigth ------------------------------
*
         musheinv = abundelec(1)
         mushe = abundelecinv
         mushiinv = sumNi(1)
         mushi = sumNinv
         mushinv = abundelec(1) + sumNi(1)
         mush = truncinv(mushinv)
*        dmu/dln
         mushf = -mush * mush * (sumNi(3)+abundelec(3))
         mushT = -mush * mush * (sumNi(2)+abundelec(2))
         term1 = mushf*mushinv
         term2 = mushT*mushinv
         CALL transdftodrho(.true.,term1,term2,ksirho,ksiT,rho(2),
     &rho(3),rhofinv)
*
*        ---------- khiro,khit,khimu ------------------------------
*
c         khimush = -pgaz(1) * mushi * mushinv * pshinv
         CALL compkimu_sod(ah,bh,sdelh,sumh,diffh,heavelec,
     &dheavelecdY,Y,mushi,mushe,mush,abundelec,sumNi,zrhosh,ksirho,
     &pgaz,pshinv,khimush)
         khirho = zrhosh - khimush * ksirho
         khiT = zTsh - khimush * ksiT

*        ---------- delta, phi from Kippenhahn & Weigert -------------
*        ---------- Needed for angular momentum transport ------------

         deltaKS(shell) = zTsh/zrhosh
         phiKS(shell) = -khimush/zrhosh

*
*        ---------- specific heat ... ----------------------------------
*
*        intermediate terms (int1, int2,...)
         zaidfinv = truncinv(zaidf)
         term1 = zaidf * ssh(2)
         term2 = zaidT * ssh(3)
         int1 = truncdiff(term1,term2)
         int1inv = truncinv(int1)
         term1 = rho(3) * ssh(2)
         term2 = rho(2) * ssh(3)
         int2 = truncdiff(term1,term2)
*        thermodynamics quantities
         cpsh = int1 * zaidfinv
         cvsh = int2 * rhofinv
         gamma1sh = truncinv(int2) * int1
         Qsh = -Tsh * rho(1) * ssh(3) * pfshinv
         adadsh = -ssh(3) * int1inv
         vssh = DSQRT(psh(1)*gamma1sh*rhoinv)
cu         gammash = truncinv(cvsh) * cpsh
*
*        ---------- Other quantities ----------
*
         srsh = ssh(1) * mush / rk
         degpesh = pelec(1) * abundelecinv * rhoinv / RT
         betapsh = (psh(1)-pradi(1)) * pshinv
         betanicorsh = (psh(1)-dpe(1)) * pshinv
         ratnehsh = abundelec(1) / MAX(abundH(1,1),abundZ(1,1))
c         print *, abundelec(1),abundH(1,1),abundZ(1,1)
         CALL fordiffusion_sod(abundH,abundZ,rationsh,Y)
         rhpsish = 2.0d0 * rootfinv *
     &((iue(1)*Tnodim+1.0d0)*ine(1)*ine(3)-3.0d0*ipe(1)*ipe(3)+
     &iue(3)*ine(1)*Tnodim) /
     &(ine(1)-3.0d0*ipe(1)+iue(1)*ine(1)*Tnodim+zero)
*
************************************************************************
*            VI - DERIVATIVES FOR NEWTON-RAPHSON METHOD                *
************************************************************************
*
*        d(zaidx)/dlny
         zaidTT = psh(4) * pshinv - zaidT * zaidT
         zaidff = psh(5) * pshinv - zaidf * zaidf
         zaidfT = psh(6) * pshinv - zaidf * zaidT
         int1T = zaidfT * ssh(2) + zaidf * ssh(4) -
     &zaidTT * ssh(3) - zaidT * ssh(6)
         int1f = zaidff * ssh(2) + zaidf * ssh(6) -
     &zaidfT * ssh(3) - zaidT * ssh(5)
*        derivatives of grad_ad
         term1 = ssh(6) * sfshinv
         term2 = int1T * int1inv
         adadshT = adadsh * truncdiff(term1,term2)
         term1 = ssh(5) * sfshinv
         term2 = int1f * int1inv
         adadshf = adadsh * truncdiff(term1,term2)
         IF (approxsod) THEN
*        Analytical computation of Cp, Q, zT and zro derivatives
*        in the set of variables (rho,T) with some approximations
*        (a term already derived is not derived a second time) and
*        conversion in derivatives (f,T) derivatives of u and p
            etro = -(4.d0*uradi(1)*Tinv+1.5d0*rk*ksirho*mushinv) *
     &rhoinv
            ett = 3.0d0 * Tinv * (4.0d0*uradi(1)*Tinv-rk*ksiT*mushinv)
            proT = (1.d0-ksiT) * (1.d0-ksirho) * (pgaz(1)+pelec(1)) *
     &rhoinv * Tinv
            pTT = (ksiT*(ksiT-1.0d0)*(pgaz(1)+pelec(1))+12.d0*pradi(1))*
     &Tinv * Tinv
            proro = ksirho * (ksirho-1.0d0) * (pgaz(1)+pelec(1)) *
     &rhoinv * rhoinv
*           derivatives of zro
            zrodro = zrhosh * (1.d0-zrhosh) * rhoinv + rho(1) * proro *
     &pshinv
            zrodT = -zrhosh * zTsh * Tinv + rho(1) * proT * pshinv
            dzrodfsh = 0.d0
            dzrodTsh = 0.d0
            CALL transdftodrho(.false.,dzrodfsh,dzrodTsh,zrodro,zrodT,
     &drhodT,drhodlnf,rhofinv)
            dzrodTsh = dzrodTsh * Tsh
*           derivatives of zT
            zTdro = -zrhosh * zTsh * rhoinv + Tsh * proT * pshinv
            zTdT = zTsh * (1.d0-zTsh) * Tinv + Tsh * pTT * pshinv
            dzTdTsh = dzTdTsh * Tsh
            dzTdfsh = 0.d0
            dzTdTsh = 0.d0
            CALL transdftodrho(.false.,dzTdfsh,dzTdTsh,zTdro,zTdT,
     &drhodT,drhodlnf,rhofinv)
*           derivatives of Cp
            cpadd = zTsh * psh(1) * zrhoshinv * Tinv * rhoinv
            cpro = etro + cpadd * (zTdro*zTshinv-zrodro*zrhoshinv+
     &zrhosh*rhoinv-rhoinv)
            cpT = ett + cpadd * (zTdT*zTshinv-zrodT*zrhoshinv+
     &zTsh*Tinv-Tinv)
            cpshf = 0.d0
            cpshT = 0.d0
            CALL transdftodrho(.false.,cpshf,cpshT,cpro,cpT,drhodT,
     &drhodlnf,rhofinv)
            cpshT = cpshT * Tsh
*           derivatives of Q
            qro = zTdro * zrhoshinv - zTsh * zrodro * zrhoshinv *
     &zrhoshinv
            qT = zTdT * zrhoshinv - zTsh * zrodT * zrhoshinv * zrhoshinv
            Qfsh = 0.d0
            QTsh = 0.d0
            CALL transdftodrho(.false.,Qfsh,QTsh,qro,qT,drhodT,drhodlnf,
     &rhofinv)
            QTsh = QTsh * Tsh
         ELSE
*        Analytical computation in the set of variables (f,T)
*        without any approximation (same as grad_ad previously computed)
*           derivatives of Cp
            term1 = int1T * int1inv
            term2 = zaidfT * zaidfinv
            cpshT = cpsh * truncdiff(term1,term2)
            term1 = int1f * int1inv
            term2 = zaidff * zaidfinv
            cpshf = cpsh * truncdiff(term1,term2)
*           derivatives of zT
            term1 = zaidfT
            term2 = ((zaidff-zaidf*rho(5)*rhofinv)*rho(2)+zaidf*rho(6))*
     &rhofinv
            dzTdfsh = truncdiff(term1,term2)
            term1 = zaidTT
            term2 = ((zaidfT-zaidf*rho(6)*rhofinv)*rho(2)+zaidf*rho(4))*
     &rhofinv
            dzTdTsh = truncdiff(term1,term2)
*           derivatives of zro
            term1 = zaidff
            term2 = zaidf * rho(5) * rhofinv
            dzrodfsh = truncdiff(term1,term2) * rhofinv
            term1 = zaidfT
            term2 = zaidf * rho(6) * rhofinv
            dzrodTsh = truncdiff(term1,term2) * rhofinv
*           derivatives of Q
            term1 = rho(3) + ssh(5) * sfshinv
            term2 = psh(5) * pfshinv
            Qfsh = Qsh * truncdiff(term1,term2)
            term1 = rho(2) + ssh(6) * sfshinv + 1.0d0
            term2 = psh(6)*pfshinv
            QTsh = Qsh * truncdiff(term1,term2)
         END IF
*
************************************************************************
*      VII - SAVE QUANTITIES TO BE USED IN STAREVOL                    *
************************************************************************
*
*        ---------- All derivatives are d/dln -----
*

cu : unused

         etak(shell) = etash(1)
cu         pek(shell) = pelec(1)
         ro(shell) = rho(1)
         drodt(shell) = rho(2) * rho(1)
         drodf(shell) = drhodlnf
         mu(shell) = mush
         mui(shell) = mushi
         mue(shell) = mushe
         mueinv(shell) = musheinv
         muiinv(shell) = mushiinv
         muinv(shell) = mushinv
         dmudf(shell) = mushf
         dmudT(shell) = mushT
         s(shell) = ssh(1)
         dsdf(shell) = ssh(3)
         dsdT(shell) = ssh(2)
c         sr(shell) = srsh
         e(shell) = ush(1)
         dedf(shell) = ush(3)
         dedT(shell) = ush(2)
         p(shell) = psh(1)
         dpdf(shell) = psh(3)
         dpdT(shell) = psh(2)
cu         dpdfT(shell) = psh(6)
cu         dpdTT(shell) = psh(4)
         beta(shell) = betapsh
cu         betanicor(shell) = betanicorsh
         degpe(shell) = degpesh
cu         gamma(shell) = gammash
         gamma1(shell) = gamma1sh
         cs(shell) = vssh
         abad(shell) = adadsh
         dabadf(shell) = adadshf
         dabadT(shell) = adadshT
         rhpsi(shell) = rhpsish
c         kimu(shell) = khimush
c         zro(shell) = zrhosh
c         dzrodf(shell) = dzrodfsh
c         dzrodT(shell) = dzrodTsh
c         dzrodT(shell) = dzrodTsh
         zT(shell) = zTsh
         dzTdf(shell) = dzTdfsh
         cp(shell) = cpsh
         cv(shell) = cvsh
         dcpdT(shell) = cpshT
         dcpdf(shell) = cpshf
         eta(shell) = etash(1)
         prad(shell) = pradi(1)
c         q(shell) = Qsh
c         deltaKS(shell) = Qsh
         deltaKSdf(shell) = Qfsh
         deltaKSdt(shell) = QTsh
         ratneh(shell) = ratnehsh ! test Fev 2019 pour avoir NE
c         print *, 'ratneh',ratneh(2)
c         if(ratneh(2).ne.0) stop
c         ratneh(shell) = amune/amu ! à modifier pour prendre une nouvelle variable
c          ratneh(shell) = amune0/amu ! amune0 ionisation total
c         ratneh(shell) = ro(shell)*mue(shell)/amu
         DO i=1,nsp
            zmean(shell,i) = zmeansh(i)
            SELECT CASE(INT(znuc(i)))
               CASE(1)
                  j=1
               CASE(2)
                  j=2
               CASE(6)
                  j=3
               CASE(7)
                  j=4
               CASE(8)
                  j=5
               CASE(10)
                  j=6
               CASE DEFAULT
                  j=0
            END SELECT
            IF (j.NE.0) THEN
               DO zz=1,INT(znuc(i))
                  ration(shell,j,zz) = rationsh(i,zz)
               END DO
            END IF
         END DO
      END DO
c       print *,'NE STAREVOL', amune/amu
c$$$      do i=1,nsp
c$$$         do j = 1,nmod
c$$$            write(666,'(1x,i4,1x,i4,9(1x,1pe11.4))') i,j,T(j)
c$$$     $           ,zmean(j,i)
c$$$         enddo
c$$$      enddo
c      stop
c      print *,'abundelec',abundelec,'amu',amu
      END


************************************************************************
* Calculate some useful quantities related to abundances               *
* Inputs: ysp(nsh,nsp)=isotopic abundances                             *
*         znuc(nsp) Z for each isotope                                 *
*         nsp=number of isotopes computed in the nuclear network       *
*         nbz=number of chemical elements computed in the code         *
*         shell=number of the shell where the calculus is performed    *
* Outputs: Y(nsh) = sum of Y for isotopes which have same Z.           *
*          X(nsh) = sum of X for isotopes which have same Z.           *
*          zavY = sum(ZY)/sum(Y)                                       *
*          zavYZ = sum(YZ2)/sum(YZ)                                    *
*          addY = sum(Y)                                               *
*          addYZ = sum(YZ)                                             *
************************************************************************
      SUBROUTINE chemY(shell,xsp,ysp,X,Y,addY,addYZ,zavY,zavYZ,zmeansh)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'evolcom.nuc'
      include 'eoscom.ion'
      include 'eoscom.par'
      include 'eoscom.cons'
      INTEGER shell
*     Intermediate variables
      INTEGER iso,i
      DOUBLE PRECISION addYZZ,addX,addYZ
*     ouputs
      DOUBLE PRECISION Y(nbz),X(nbz),addY,zavY,zavYZ,zmeansh(nsp)
      DOUBLE PRECISION xsp(nsh,nsp),ysp(nsh,nsp)
*     ----------initialisation----------
      addYZ = 0.0d0
      addYZZ = 0.0d0
      addY = 0.0d0
      addX = 0.0d0
      zavY = 0.0d0
      zavYZ = 0.0d0
      DO i=1,nbz
         X(i) = 0.0d0
         Y(i) = 0.0d0
      END DO
*     ----------calculus----------
      DO iso=nspmin,nspmax
         X(INT(znuc(iso))) = X(INT(znuc(iso))) + xsp(shell,iso)
         Y(INT(znuc(iso))) = Y(INT(znuc(iso))) + ysp(shell,iso)
         addX = addX + xsp(shell,iso)
         zmeansh(i) = 0.0d0
      END DO
c     IF (ABS(addX-1.0d0).GT.physicalimit) THEN
c        WRITE (NOUT,*) "ERROR: EOS - chemY - sum(xsp)=",addX,"<>1 !!!"
c        WRITE (NOUT,*) "Shell=",shell
c        STOP
c     END IF
*
*     DO NOT TAKE "HEAVY" INTO ACCOUNT IN Z AVERAGE
*     (there is no good Z for it!)
*      ==> nbz-1
*
      DO i=1,nbz-1
         IF ((shell.eq.1).AND.checkcomposition) THEN
            WRITE (NOUT,"(A,I2,A,D14.8)") "X(Z=",Z(i),") =",X(i)
            WRITE (NOUT,"(A,I2,A,D14.8)") "Y(Z=",Z(i),") =",Y(i)
         END IF
         addYZ = addYZ + Y(i) * Z(i)
         addYZZ = addYZZ + Y(i) * Z(i) * Z(i)
         addY = addY + Y(i)
      END DO
      IF ((shell.eq.1).AND.checkcomposition) THEN
         WRITE (NOUT,"(A,I2,A,D14.8)") "X(Z=",Z(nbz),") =",X(nbz)
         WRITE (NOUT,"(A,I2,A,D14.8)") "Y(Z=",Z(nbz),") =",Y(nbz)
      END IF
      zavY = addYZ / addY
      zavYZ = addYZZ / addYZ
*
*     BUT "HEAVY" MUST BE TAKEN INTO ACCOUNT IN THE MASS
*      (strictly speaking, neutrons too, but they are negligeable)
*
      addY = addY + Y(nbz)

      END
************************************************************************
* Choose the right solution among H1 and H2 with the two constraints:  *
*         (1) 0<H<Yh          (2) 0<H<lim                              *
* lim is choosen among limem (which ensure the charge conservation) and*
*  limH2 (which ensures the proton conservation) depending on diff     *
* One indicates if the first solution has been choosen                 *
************************************************************************
      SUBROUTINE choice(shell,H1,H2,limH2,limem,diff,Yh,H,lim,one,a,b,
     &c,rat,dev,error)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.ion'
      include 'eoscom.cons'
      include 'eoscom.par'
      DOUBLE PRECISION H1,H2,limH2,limem,diff,Yh,a,b,c,rat,dev(4)
      INTEGER shell,error
      LOGICAL egal
*     Outputs
      DOUBLE PRECISION H,lim
      LOGICAL one
*     ----------If the proposed value exceed of a little ----------
*     ----------amount Yh, it is assumed to be Yh ----------
      IF (egal(H1,Yh).AND.(H1*Yh.GE.0.d0).AND.(H1.GT.Yh)) THEN
         H1 = Yh
      END IF
      IF (egal(H2,Yh).AND.(H2*Yh.GE.0.d0).AND.(H2.GT.Yh)) THEN
         H2 = Yh
      END IF
*     ----------If H2 is conputed----------
      IF (compH2) THEN
*        ----------Limits to be checked----------
         IF (diff.GT.0.d0) THEN
            lim = limH2
         ELSE
            IF (limem.LT.limH2) THEN
               lim = limem
            ELSE
               lim = limH2
            END IF
         END IF
*        ----------Limit values to be recognized----------
         IF (egal(H1,limH2).AND.(H1*limH2.GE.0.d0).AND.(H1.GT.limH2))
     &        THEN
            H1 = limH2
         END IF
         IF (egal(H1,limem).AND.(H1*limem.GE.0.d0).AND.(H1.GT.limem))
     &        THEN
            H1 = limem
         END IF
         IF (egal(H2,limH2).AND.(H2*limH2.GE.0.d0).AND.(H2.GT.limH2))
     &        THEN
             H2 = limH2
         END IF
         IF (egal(H2,limem).AND.(H2*limem.GE.0.d0).AND.(H2.GT.limem))
     &        THEN
            H2 = limem
         END IF
*        ----------Two egal values----------
         IF (egal(H1,H2).AND.(H1*H2.GE.0.d0)) THEN
            H = H1
            H2 = H1
         ELSE
*     ----------Two different values----------
            IF ((H1.LT.0.d0).OR.(H1.GT.lim)) THEN
               IF ((H2.LT.0.d0).OR.(H2.GT.lim)) THEN
                  WRITE (NOUT,*) "ERROR: EOS - no physical value!"
                  WRITE (NOUT,*) "Shell=",shell," H1=",H1," H2=",H2
                  WRITE (NOUT,*) "limem=",limem," limH2=",limH2,"lim=",
     &                 lim
                  WRITE (NOUT,*) "a=",a," b=",b," c=",c
                  WRITE (NOUT,*) "rath2h=",rat
                  ERROR = 11
                  return
               ELSE
                  H = H2
               END IF
            ELSE
               IF ((H2.LT.0.d0).OR.(H2.GT.lim)) THEN
                  H = H1
               ELSE
                  IF ((H1.EQ.0.d0).AND.(H2.GT.0.d0)) THEN
                     H = H2
                  ELSE
                     IF ((H2.EQ.0.d0).AND.(H1.GT.0.d0)) THEN
                        H = H1
                     ELSE
                        WRITE (NOUT,*) "EOS - two possible values!"
                        WRITE (NOUT,*) "Shell=",shell," H1=",H1," H2=",
     &                       H2
                        WRITE (NOUT,*) "limem=",limem," limH2=",limH2
                        ERROR = 11
                        return
                     END IF
                  END IF
               END IF
            END IF
         END IF
*     ----------IF H2 is not computed----------
      ELSE
         H = limH2
         lim = limH2

      END IF
      one = (H.EQ.H1)
*     ----------To check the selected value----------
      IF (one.AND.dev(1).GT.numericalimit) THEN
         WRITE (NOUT,*) "ERROR - choice - possible bad solution choosen"
         WRITE (NOUT,*) "Shell=",shell,"H=",H1,"dev=",dev(1)
         error = 10
         return
      ENDIF
      IF (.not.one.AND.dev(2).GT.numericalimit) THEN
         WRITE (NOUT,*) "ERROR - choice - possible bad solution choosen"
         WRITE (NOUT,*) "Shell=",shell,"H=",H2,"dev=",dev(2)
         error = 10
         return
      ENDIF
      IF (H.LT.0.d0) THEN
         WRITE (NOUT,*) "ERROR: EOS - choice - Nh < 0 !"
         WRITE (NOUT,*) "shell=",shell," diff=",diff
         WRITE (NOUT,*) "limH2=",limH2," limem=",limem," H=",H
         error = 10
         return
      END IF
      IF (H.GT.Yh) THEN
         WRITE (NOUT,*) "ERROR: EOS - choice - Nh > Yh !"
         WRITE (NOUT,*) "shell=",shell," H1=",H1," H2=",H2
         WRITE (NOUT,*) "H=",H," Yh=",Yh
         error = 10
         return
      END IF
      IF (H.GT.lim) THEN
         WRITE (NOUT,*) "ERROR: EOS - choice - Ne- or NH2 < 0 !"
         WRITE (NOUT,*) "Shell=",shell," diff=",diff
         WRITE (NOUT,*) "limH2=",limH2," limem=",limem," H=",H
         error = 10
         return
      END IF
      END

************************************************************************
* Determination of the CNO partition functions                         *
*  Reference:                                                          *
************************************************************************

      SUBROUTINE cno_sod(T,ZZ,nZ,lnparfunZ)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.ion'
      include 'eoscom.phys'
      include 'eoscom.cons'
      DOUBLE PRECISION T
*     Intermediate variables
      INTEGER i,j,debutE,qeff,ZZ,indx,rangE,debutA,rangA,nZ
      DOUBLE PRECISION pu(12,24),puu(12,24),truncln,truncinv
      DOUBLE PRECISION arge,fac,ratio,denominv,root,ratioT,ratioTT
      DOUBLE PRECISION nu,nus2,nus,nus2T,nus2TT,term,terminv
      DOUBLE PRECISION qu(5,24),quu(5,24),itst(20),sum,sumT,sumTT
c     DOUBLE PRECISION sumTTT
      DOUBLE PRECISION term1,term2,term3,term3T,term2T,term2TT,term3TT
      DOUBLE PRECISION a1,a2,a3,truncexpo,dln10,pe,theta,theta10
      DOUBLE PRECISION psi,psir,psir1,psir2,psirT,psiT,psiTT,psirTT

      PARAMETER (dln10=2.3025851d0)
      PARAMETER (a1=5040.0d0)
      PARAMETER (a2=165.2247d0)
      PARAMETER (a3=31.321d0)
      PARAMETER (pe=1.0d1)
*     outputs
      include 'eoscom.math'
      DOUBLE PRECISION lnparfunZ(nbionZmax,deriv3)
*     --------------------
      DATA ((pu(j,i),i = 1,24),j = 1,12)/8.0158d0,4.0003d0,
     &10.0281d0,6.0057d0,15.7995d0,2.1000d0,1.0000d0,14.0499d0,
     &8.0462d0,4.0003d0,10.3289d0,6.0044d0,1.1000d0,2.1000d0,1.0000d0,
     &4.0029d0,12.7843d0,8.0703d0,4.0003d0,10.5563d0,2.1000d0,1.1000d0,
     &2.1000d0,1.0000d0,5.8833d0,17.0841d0,15.7574d0,23.5757d0,
     &36.2005d0,0.0d0,0.0d0,30.8008d0,6.2669d0,19.3533d0,14.5021d0,
     &23.5612d0,0.0d0,0.0d0,0.0d0,5.3656d0,5.6828d0,5.7144d0,21.2937d0,
     &13.2950d0,2.0000d0,0.0d0,0.0d0,0.0d0,33.7521d0,82.9154d0,
     &186.2109d0,76.4185d0,0.0d0,0.0d0,0.0d0,883.1443d0,17.8696d0,
     &80.6462d0,187.1624d0,76.4344d0,0.0d0,0.0d0,0.0d0,36.2853d0,
     &98.0919d0,84.1156d0,78.7058d0,188.1390d0,4.0000d0,0.0d0,0.0d0,
     &0.0d0,595.3432d0,15.9808d0,15.4127d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &10.0000d0,282.8084d0,13.0998d0,108.1615d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &1044.3447d0,829.4396d0,529.0927d0,12.8293d0,14.6560d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,48.2044d0,55.9559d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &16.0000d0,7.3751d0,19.6425d0,191.8383d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &131.0217d0,50.9878d0,5.6609d0,16.2730d0,129.4922d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,435.8093d0,243.6311d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &64.0000d0,33.1390d0,94.3035d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &868.9779d0,199.0120d0,28.9355d0,123.6578d0,470.8512d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &215.4829d0,370.9539d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,14.8533d0,
     &2.0000d0,111.3620d0,327.2396d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,16.0000d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,93.1466d0,6.0000d0,494.0413d0,
     &48.7883d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,38.0000d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,10.0000d0,45.5249d0,102.2117d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,10.0000d0,
     &134.4751d0,20.0000d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,30.0000d0,0.0d0,161.9903d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,50.0000d0,0.0d0,
     &28.4184d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/
*
      DATA ((puu(j,i),i = 1,24),j = 1,12)/0.004d0,.008d0,6.691d0,
     &8.005d0,303.772d0,0.0d0,0.0d0,2.554d0,0.014d0,0.022d0,8.693d0,
     &9.999d0,0.0d0,0.0d0,0.0d0,0.022d0,3.472d0,0.032d0,0.048d0,
     &10.747d0,0.0d0,0.0d0,0.0d0,0.0d0,1.359d0,16.546d0,25.034d0,
     &40.804d0,354.208d0,0.0d0,0.0d0,9.169d0,2.131d0,31.259d0,37.650d0,
     &60.991d0,0.0d0,0.0d0,0.0d0,2.019d0,7.437d0,2.760d0,50.089d0,
     &52.323d0,11.949d0,0.0d0,0.0d0,0.0d0,6.454d0,21.614d0,40.975d0,
     &54.492d0,0.0d0,0.0d0,0.0d0,13.651d0,15.745d0,41.428d0,65.479d0,
     &82.262d0,0.0d0,0.0d0,0.0d0,9.812d0,22.579d0,35.328d0,66.604d0,
     &94.976d0,12.015d0,0.0d0,0.0d0,0.0d0,10.376d0,5.688d0,17.604d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,12.353d0,24.949d0,7.212d0,61.155d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,13.087d0,32.035d0,48.277d0,8.954d0,27.405d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,15.801d0,36.180d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,13.784d0,6.376d0,15.228d0,79.196d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &13.804d0,27.774d0,7.662d0,18.031d0,86.350d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,26.269d0,47.133d0,0.0d0,0.0d0,0.0d0,0.0d0,14.874d0,
     &14.246d0,34.387d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,16.061d0,33.678d0,
     &16.786d0,57.755d0,109.917d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,29.465d0,46.708d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,14.293d0,28.118d0,42.657d0,72.594d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,46.475d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,16.114d0,31.019d0,
     &54.522d0,68.388d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,49.468d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,34.204d0,50.204d0,82.397d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,30.892d0,56.044d0,
     &31.960d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,33.189d0,0.0d0,76.876d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,36.181d0,0.0d0,75.686d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0/
*
      DATA ((qu(j,i),i = 1,24),j = 1,5) /12.0d0,2.0d0,4.0d0,2.0d0,4.0d0,
     &0.0d0,0.0d0,18.0d0,12.0d0,2.0d0,4.0d0,2.0d0,0.0d0,0.0d0,0.0d0,
     &8.0d0,18.0d0,12.0d0,2.0d0,4.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &18.0d0,12.0d0,0.0d0,0.0d0,0.0d0,0.0d0,10.0d0,24.0d0,18.0d0,12.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,20.0d0,10.0d0,24.0d0,18.0d0,12.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,6.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,12.0d0,2.0d0,20.0d0,
     &6.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 10.0d0,0.0d0,18.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,10.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0/
*
      DATA ((quu(j,i),i = 1,24),j = 1,5) /6.0d0,6.0d0,6.0d0,6.0d0,4.0d0,
     &0.0d0,0.0d0,6.1d0,5.0d0,6.0d0,6.0d0,6.0d0,0.0d0,0.0d0,0.0d0,8.0d0,
     &6.0d0,6.0d0,5.0d0,6.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,5.0d0,5.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,4.0d0,3.0d0,5.0d0,6.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,6.0d0,5.0d0,4.9d0,5.0d0,6.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,4.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,3.0d0,3.0d0,4.0d0,4.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,3.0d0,0.0d0,4.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,4.0d0,
     &0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/
*
      DATA itst  /1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,1.0d0,1.0d0,
     &1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0/
*     --------------------
*     Variables des fits
      theta = a1 / T
      theta10 = -theta * dln10
      psiT = DSQRT(theta*DSQRT(pe))
      psi = truncinv(psiT)
      psi = psi * a2
*     dpsi/dlnT
      psiT = 0.5d0 * psi
      psiTT = 0.5d0 * psiT
      nu = a3 * theta
      debutE = indx(2,ZZ-1)
      debutA = indx(6,ZZ-1)
      DO qeff = 1,ZZ+1
         rangA = debutA + qeff
         rangE = debutE + qeff
         sum = 0.0d0
         sumT = 0.0d0
         sumTT = 0.0d0
         DO j=1,12
            sum = sum + pu(j,rangA) * truncexpo(puu(j,rangA)*theta10)
            sumT = sumT - puu(j,rangA) * theta10
            sumTT = sumTT + puu(j,rangA) * theta10
         END DO
         IF (rangA.LE.20) THEN
            arge = (enerlevelZ(rangE+1)-enerlevelZ(rangE)) * theta10
            DO j = 1,5
               root = DSQRT(DBLE(j))
               nus = nu * DBLE(j*j)
               psir = psi * root
               psir1 = psir + 1.0d0
               psir2 = psir + 0.5d0
*              dpsir/dlnT
               psirT = psiT * root
               psirTT = psiTT * root
               nus2 = nus * (1.0d0+0.5d0*nus)
*              dnus2/dlnT
               nus2T = -nus2 - 0.5d0 * nus * nus
               nus2TT = -nus2T + nus * nu
               term3 = inv3 * psir * (psir+1.0d0) * (psir+0.5d0)
*              dterm3/dlnT
               term3T = inv3 * psirT * (psir1*psir2+psir*psir2+
     &psir*psir1)
               term3TT = inv3 * psirTT * (psir1*psir2+psir*psir2+
     &psir*psir1) + inv3 * psirT * psirT * 2.0d0 * (psir2+psir1+psir)
               term1 = inv3 * quu(j,rangA) *
     &(quu(j,rangA)-1.0d0) * (quu(j,rangA)-0.5d0)
               term2 = psir * (quu(j,rangA)-1.0d0)
               denominv = truncinv(term2)
               ratio = (psir1-quu(j,rangA)) * denominv
               ratioT = psirT * denominv * (1.0d0-ratio)
               ratioTT = denominv * (psirTT*(1.0d0-ratio)-
     &2.0d0*psirT*ratioT)
               term2 = nus2 * ratio
*              dterm2/dlnT
               term2T = nus2T * ratio - nus2 * ratioT
               term2TT = nus2TT * ratio - nus2 * ratioTT
               term = term3-term1+term2
               terminv = truncinv(term)
               fac = itst(rangA) * qu(j,rangA)
               ratio = (term3T+term2T) * terminv
               sum = sum + fac * term * truncexpo(arge)
               sumT = sumT + ratio - arge
               sumTT = sumTT + terminv * ((term3TT+term2TT)-
     &ratio*(term3T+term2T)) + arge
            END DO
         END IF
         lnparfunZ(nz+qeff,1) = truncln(truncexpo(lnparfunZ(rangA,1))+
     &sum)
         lnparfunZ(nz+qeff,2) = sumT
         lnparfunZ(nz+qeff,4) = sumTT
c        lnparfunZ(nz+qeff,7) = sumTTT
         lnparfunZ(nz+qeff,7) = 0.d0
      END DO
      END

************************************************************************
* Computation of khimu = dlnP/dlnmu|rho,T                              *
*     with three changes of variables:                                 *
*  (lnf,T,Yi)==>(lnf,T,mu)==> (rho,T,mu) <==(rho,T)                    *
* Kimu is computed using a specified specie whose index is elem.       *
************************************************************************

      SUBROUTINE compkimu_sod(ah,bh,sdelh,sumh,diffh,heavelec,
     &dheavelecdY,Y,mushi,mushe,mush,abundelec,sumNi,zrho,ksirho,pgaz,
     &pshinv,khimu)
*     --------------------
      IMPLICIT NONE
      include 'evolpar.star'
      include 'eoscom.math'

*     inputs
      DOUBLE PRECISION ah,bh,sdelh,sumh,diffh,heavelec(deriv2)
      DOUBLE PRECISION dheavelecdY(nbz),Y(nbz),abundelec(deriv2)
      DOUBLE PRECISION sumNi(deriv2),zrho,ksirho,pgaz(deriv2),pshinv
      DOUBLE PRECISION mushe,mushi,mush
*     intermediate variables
      DOUBLE PRECISION c1,c2
*     standards derivatives dX/dY
      DOUBLE PRECISION dHdY,dH2dY
*     dlnX/dlnY or dlnX/lnmu
      DOUBLE PRECISION dmuidY,dmuedY,dmudY,drhodY,dphidY,dphidmu,drhodmu
*     Output
      DOUBLE PRECISION khimu
*     --------------------

      c1 = (bh-sdelh) / (2.0d0*ah*sdelh)
      c2 = 1.0d0 / sdelh
      dHdY = c1 * (sumh*dheavelecdY(1)-diffh) +
     &c2 * (heavelec(1)+Y(1)*dheavelecdY(1))
      dH2dY = 0.5d0 * (1.d0-sumh*dHdY)
      dmuidY = -Y(1) * mushi * (Y(1)-dH2dY)
      dmuedY = -Y(1) * mushe *  (diffh*dHdY+dheavelecdY(1))
      dmudY = mush * (abundelec(1)*dmuedY+sumNi(1)*dmuidY)
      drhodY = dmuedY
      dphidY = pgaz(1) * pshinv * (drhodY-dmuidY)
      dphidmu = dphidy / dmudy
      drhodmu = drhody / dmudy
      khimu = (dphidmu-zrho*drhodmu) / (1.0d0-drhodmu*ksirho)
      END

************************************************************************
* Coulomb shielding corrections computation                            *
*  These effects are fitted by means of function g(x,y) (here coco) /  *
*   DF = - Ne * k * T * g(x,y), with x=ne*mh and y=13.6/kT             *
* Inputs: x=ne*mh y=13.6/kT,f,rootfinv,zavY,zavYZ,halff1               *
* Outputs: coco and derivatives versus x and y                         *
* Reference: Slattery, Doolen, De Witt (Phys rev A 21, 2087, 1980)     *
************************************************************************

      SUBROUTINE coulomb_sod(f,x,y,zavY,zavYZ,coco,cocox,cocoy,
     &cocoxx,cocoyy,cocoxy,cocoxxx,cocoyyy,cocoyyx,cocoxxy,etapx,
     &term,sum,sumx,halff1,ine,rootfinv,fadd1inv)
*     ----------declarations----------
      IMPLICIT NONE
*     inputs
      include 'eoscom.math'
      include 'evolpar.star'
      include 'eoscom.ion'
      include 'eoscom.cons'
      include 'eoscom.par'
      DOUBLE PRECISION x,y,zavY,zavYZ,f,ine(deriv2),rootfinv
      DOUBLE PRECISION halff1,fadd1inv,etapx,term,sum,sumx
*     ouputs
      DOUBLE PRECISION coco,cocox,cocoy,cocoxx,cocoyy,cocoxy
      DOUBLE PRECISION cocoxxx,cocoyyy,cocoyyx,cocoxxy
*     Intermediate variables
      DOUBLE PRECISION theta,thetax,thetay,thetaxx,thetaxy,thetaxxx
      DOUBLE PRECISION thetayy,thetayyy,thetaxxy,thetayyx
      DOUBLE PRECISION deno,deno1,deno2,deno1m,deno2m,denomm,denommm
      DOUBLE PRECISION denom,memb,membf,membT,truncinv
      DOUBLE PRECISION ratio1,ratio2,ratio1x,ratio2y,ratio3,ratio4
      DOUBLE PRECISION produ1,produ2,produ3,inv,diff,termf,produ1f
      DOUBLE PRECISION gama,gamax,gamay,gamaxx,gamayy,gamaxy,gamax2
      DOUBLE PRECISION gamaxxx,gamaxxy,gamayyy,gamayyx,gamay2
      DOUBLE PRECISION g,gm,gmm,gmmm,suminv,deninv,ineTTTf,ineffTT
      DOUBLE PRECISION inefff,ineTTT,ineffT,ineTTf,inefinv,ineffff
      DOUBLE PRECISION inefffT,ratio5,combf2,combf3,produ1T
*     coefficient
      DOUBLE PRECISION a1,a2,a3,a2inv
*     ----------constantes----------
*     pressure ionization fit coefficients
      PARAMETER (a1=0.897520d0)
      PARAMETER (a2=0.768d0)
      PARAMETER (a3=0.208d0)
      PARAMETER (a2inv=1.3020833d0)
*     ----------theta=dlnne/dlneta computation----------
      IF (fapprocoul) THEN
*        WARNING:Based on f approximation==>false!!
         theta = truncinv(etapx)
         thetax = term * theta
         thetay = 1.5d0 * thetax
         thetaxx = (term*term-sum) * theta
         thetaxy = 1.5d0 * thetaxx
         thetayy = 1.5d0 * thetaxy
         thetaxxx = (term*(term*term-3.0d0*sum)-sumx) * theta
         thetaxxy = 1.5d0 * thetaxxx
         thetayyx = 1.5d0 * thetaxxy
         thetayyy = 1.5d0 * thetayyx
      ELSE
*        Based on f given by main routine==>true!!
*        WARNING: The next derivatives are computed with the assumption
*         that the derivatives of Ine whose order higher than 2 are nul.
*         It doesn't change the convergence with the Newton-Raphson
*         method.
*        To limit the difference between numerical and analytical
*        derivatives at low f, inefff is choosen to be 1d-40.
         inefff = 1.0d-40
         ineTTT = 0.0d0
         ineffT = 0.0d0
         ineTTf = 0.0d0
         ineffff = 0.0d0
         inefffT = 0.0d0
         ineffTT = 0.0d0
         ineTTTf = 0.0d0
         inefinv = truncinv(ine(3))
         ratio1 = ine(5) * inefinv
         ratio2 = ine(6) * inefinv
         ratio3 = inefff * inefinv
         ratio4 = ine(2) * inefinv
         ratio5 = ineffT * inefinv
         combf2 = halff1 * fadd1inv
         combf3 = halff1 * rootfinv
         diff = ratio3 - ratio1 * ratio1
         sum = ratio5 - ratio1 * ratio2
         term = diff - combf2
         termf = ineffff * inefinv - ratio1 * ratio3 -
     &2.0d0 * ratio1 * diff - combf2 * (1.0d0-f) * fadd1inv
         memb = ratio2 * ine(2) - ine(4)
         membf = sum * ine(2) + ratio2 * ine(6) - ineTTf
         membT = ine(4) * ratio2 + ine(2) *
     *(ineTTf*inefinv-ratio2*ratio2)-ineTTT
         produ1 = ratio5 * ine(2) - ineTTf
         produ1f = inefffT * ratio4 + ineffT * (ratio2-ratio4*ratio1) -
     &ineffTT
         produ1T = ineffTT * ratio4 + ine(6) * ratio5 * (1.0d0-ratio4) -
     &ineTTTf
*        Theta and derivatives
         theta = ine(3) * rootfinv
         thetax = (ratio1-halff1) * rootfinv
         thetay = thetax * ine(2) - rootfinv * ine(6)
         thetaxx = inefinv * (rootfinv*term-halff1*thetax)
         thetaxy = thetaxx * ine(2) - rootfinv * sum
         thetayy = thetaxy * ine(2) + thetax * memb - rootfinv *
     &(ratio5*ine(2)-ineTTf)
         thetaxxx = inefinv * ((rootfinv*termf-combf3*term-
     &combf2*thetax)*inefinv-(halff1+ratio1)*thetaxx)
         thetaxxy = thetaxxx * ine(2) + thetaxx * ratio2 + combf3 *
     &(ratio5*inefinv-ratio1*ratio2) - rootfinv * inefinv *
     &(inefffT*inefinv-ineffT*ratio1*(1.0d0+inefinv)-ine(6)*diff)
         thetayyx = thetaxxy * ine(2) + thetaxy * ratio2 +
     &thetaxx * memb+
     &(thetax*membf+combf3*produ1+rootfinv*produ1f) * inefinv
         thetayyy = thetayyx * ine(2) + 2.0d0 * thetaxy * memb +
     &thetax * (membf*ratio4-membT) - rootfinv *(produ1f*ratio4-produ1T)
      END IF
*     ----------affectations----------
      sum = theta + zavYZ
      suminv = truncinv(sum)
      ratio1 = thetax * suminv
      ratio2 = thetay * suminv
      produ1 = ratio1 * ratio1
      produ2 = ratio2 * ratio2
      ratio1x = thetaxx * suminv
      ratio2y = thetayy * suminv
*     ----------gamma function and derivatives---------
      gama = gamcoeff * x**inv3 * y * zavY**(2*inv3) * sum
c      print *,gama
      gamax = inv3 + ratio1
      gamay = 1.0d0 + ratio2
      gamaxx = ratio1x - produ1
      gamaxy = thetaxy * suminv - ratio1 * ratio2
      gamayy = ratio2y - produ2
      gamaxxx = thetaxxx * suminv - ratio1 * (ratio1x+2.0d0*gamaxx)
      gamayyy = thetayyy * suminv - ratio2 * (ratio2y+2.0d0*gamayy)
      gamaxxy = thetaxxy * suminv - ratio2 * ratio1x -
     &2.0d0 * ratio1 * gamaxy
      gamayyx = thetayyx * suminv - ratio1 * ratio2y -
     &2.0d0 * ratio2 * gamaxy
*     ----------affectations----------
      produ1 = gama + a3
      inv = truncinv(produ1)
      produ1 = gama * inv
      deno1 = produ1**a2inv
      deno1m = a3 * inv * a2inv
      deno2 = (a1*DSQRT(3.0d0/gama))**a2inv
      deno2m = -0.5d0 * a2inv
      deno = deno1 + deno2
      deninv = truncinv(deno)
      produ2 = deno1 * deno1m
      produ3 = deno2m * deno2m
      denom = produ2 + deno2 * deno2m
      denomm = produ2 * (deno1m-produ1) + deno2 * produ3
      denommm = produ2 * (deno1m*deno1m-3.0d0*deno1m*produ1+
     &produ1*inv*(gama-a3)) + deno2 * produ3 * deno2m
      ratio1 = denom * deninv
      ratio2 = a2 * deninv
      diff = denomm - ratio1 * denom
*     ----------gc and derivatives versus gamma----------
*     WARNING: gm=dlng/dlngama;gmm=d2g/dlngama2;gmmm=d3g/dlngama3
      g = a1 * gama * deninv**a2 / zavY
      gm = 1.0d0 - a2 * ratio1
      produ1 = gm * gm
      gmm = g * (produ1-diff*ratio2)
      gmmm = 3.0d0 * gmm * gm - 2.0d0 * g * produ1 * gm +
     &g * (diff*ratio1*ratio2-ratio2*denommm+
     &2.0d0*ratio2*ratio1*denomm-ratio2*ratio1*ratio1*denom)
*     gm is no more dlng/dlngama but dg/dlngama
      gm = gm * g
*     ----------gc and derivatives versus x and y----------
      coco = g
      cocox = gm * gamax
      cocoy = gm * gamay
      gamax2 = gamax * gamax
      gamay2 = gamay * gamay
      cocoxx = gmm * gamax2 + gm * gamaxx
      cocoyy = gmm * gamay2 + gm * gamayy
      cocoxy = gmm * gamax * gamay + gm * gamaxy
      cocoxxx = gmmm * gamax2 * gamax + 3.0d0 * gmm * gamax * gamaxx +
     &gm * gamaxxx
      cocoyyy = gmmm * gamay2 * gamay + 3.0d0 * gmm * gamay * gamayy +
     &gm * gamayyy
      cocoyyx = gmmm * gamay2 * gamax + gmm *
     &(2.0d0*gamay*gamaxy+gamax*gamayy) + gm * gamayyx
      cocoxxy = gmmm * gamax2 * gamay + gmm *
     &     (2.0d0*gamax*gamaxy+gamay*gamaxx) + gm * gamaxxy
c      print *, 'ne*mh',x,mh
c      stop
      END
************************************************************************
* Function for degeneracy parameter calculus                           *
* Inputs : lnf=LN(f) and rootf=SQRT(1+f) to avoid a second calculus    *
* Output : eta, degeneracy parameter                                   *
************************************************************************
      DOUBLE PRECISION FUNCTION degen(lnf,rootf)
*     -------------------
      IMPLICIT NONE
      DOUBLE PRECISION rootf,lnf,truncln
*     -------------------
      degen = lnf + 2.0d0 * (rootf-truncln(rootf+1.0d0))
      END

************************************************************************
* Dimension electronic thermodynamics quantities                       *
* Inputs : abundelec=Ne*amu,RT=R*T,adimensionnal quantitites and       *
*         derivatives (ine,ipe,ise,iue,inef,inet...)                   *
* Outputs: pelec, selec, uelec (deriv2), all derivatives (first and    *
*          second order are d/dln                                      *
************************************************************************

      SUBROUTINE dimelec_sod(abundelec,RT,ise,ipe,iue,pelec,selec,uelec)
*     ----------declarations----------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'evolcom.cons'
      include 'eoscom.cons'
      include 'eoscom.math'
      DOUBLE PRECISION abundelec(deriv2),RT
      DOUBLE PRECISION ipe(deriv2),ise(deriv2),iue(deriv2)
*     ouputs
      DOUBLE PRECISION pelec(deriv2),selec(deriv2),uelec(deriv2)
*     --------------------
      pelec(1) = mc28pl3 * ipe(1)
      pelec(2) = pelec(1) * ipe(2)
      pelec(3) = pelec(1) * ipe(3)
      pelec(4) = (ipe(4)+ipe(2)*ipe(2)) * pelec(1)
      pelec(5) = (ipe(5)+ipe(3)*ipe(3)) * pelec(1)
      pelec(6) = (ipe(6)+ipe(3)*ipe(2)) * pelec(1)
*     --------------------
      selec(1) = rk * abundelec(1) * ise(1)
      selec(2) = rk * abundelec(2) * ise(1) + rk * abundelec(1) * ise(2)
      selec(3) = rk * abundelec(3) * ise(1) + rk * abundelec(1) * ise(3)
      selec(4) = rk * abundelec(4) * ise(1) +
     &2.0d0 * rk * abundelec(2) * ise(2) + rk * abundelec(1) * ise(4)
      selec(5) = rk * abundelec(5) * ise(1) +
     &2.0d0 * rk * abundelec(3) * ise(3) + rk * abundelec(1) * ise(5)
      selec(6) = rk * abundelec(6) * ise(1) + rk * abundelec(2) *ise(3)+
     &rk * abundelec(3) * ise(2) + rk * abundelec(1) * ise(6)
*     --------------------
      uelec(1) = RT * abundelec(1) * iue(1)
      uelec(2) = RT * (iue(1)*abundelec(2)+abundelec(1)*iue(2)) +
     &uelec(1)
      uelec(3) = RT * (iue(1)*abundelec(3)+abundelec(1)*iue(3))
      uelec(4) = 0.d0
      uelec(5) = 0.d0
      uelec(6) = 0.d0
      END

************************************************************************
* Non-perfect correction for pressure                                  *
* Inputs: g and derivatives in fitxy, derivatives of x=amune: xt,xf... *
*  intermediates quantities already calculated (dmut1,...)             *
* Ouputs: Dpe(deriv2) and derivatives d()/dln                          *
*  WARNING: Dpe are not initialized but modified. If num=1 this is an  *
*           addition and if num=-1, this is an difference.             *
************************************************************************

      SUBROUTINE dpecalc_sod(num,fitxy,x,xf,xt,xff,xtt,xft,RT,dmut2,
     &     dmut1x,dmut2x,dpe)
*     ----------declarations----------
      IMPLICIT NONE
*     inputs
      include 'eoscom.math'
      DOUBLE PRECISION fitxy(deriv3),x,xf,xff,xfT,xt,xtt,RT
      INTEGER num
      DOUBLE PRECISION dmut2,dmut1x,dmut2x
*     Intermediate variables
      DOUBLE PRECISION dmuf,dmuft,dmuff,pecor,petcor,pefcor
*     ouputs
      DOUBLE PRECISION dpe(deriv2)
*     --------------------
      dmuf = -xf * dmut2
      dmuff = -xff * dmut2 - xf * xf * dmut2x
      dmuft = xf * dmut1x - xft * dmut2 - xf * xt * dmut2x
      pecor = -x * fitxy(2) * RT
      pefcor = dmuf * x * RT
      petcor = pecor - x * (dmut2*xt-fitxy(6)) * RT
      dpe(1) = dpe(1) + num * pecor
      dpe(3) = dpe(3) + num * pefcor
      dpe(2) = dpe(2) + num * petcor
      dpe(5) = dpe(5) + num * x * (dmuff+dmuf*xf) * RT
      dpe(4) = dpe(4) + num * (petcor*(2.0d0+xt)-pecor*(xt+1.0d0)-
     &(((dmut2x*xt-fitxy(9)-dmut1x)*xt+dmut2*xtt+fitxy(10))*x*RT))
      dpe(6) = dpe(6) + num * (pefcor+(dmuft+xt*dmuf)*x*RT)
      END

************************************************************************
* Non-perfect correction for entropy                                   *
* Inputs: g and derivatives in fitxy, derivatives of x=nemh: xt,xf...  *
*  intermediates quantities already calculated (dmut1,...)             *
* Ouputs: Dse(deriv2) and derivatives d()/dln                          *
*  WARNING: Dse are not initialized but modified. If num=1 this is an  *
*           addition and if num=-1, this is an difference.             *
************************************************************************

      SUBROUTINE dsecalc_sod(num,fitxy,xf,xt,xtt,xff,xft,dst1,dst2,
     &dst1x,dst2y,dst2x,abundelec,dse)
*     ----------declarations----------
      IMPLICIT NONE
*     inputs
      include 'eoscom.math'
      include 'evolcom.cons'
      DOUBLE PRECISION fitxy(deriv3),xf,xt,xtt,xff,xft
      DOUBLE PRECISION abundelec(deriv2)
      DOUBLE PRECISION dst1,dst2,dst1x,dst2y,dst2x
      INTEGER num
*     Intermediate variables
      DOUBLE PRECISION dseNk,dseNkf,dseNkT,dseNkff,dseNkTT,dseNkfT
*     outputs
      DOUBLE PRECISION dse(deriv2)
*     -------------------
      dseNk = fitxy(1) - fitxy(3)
      dseNkf = dst1 * xf
      dseNkt = dst1 * xt + dst2
      dseNkff = dst1x * xf * xf + dst1 * xff
      dseNktt = dst1x * xt * xt + dst1 * xtt + 2.0d0 * dst2x * xt -
     &dst2y
      dseNkft = dst1x * xt * xf + dst2x * xf + dst1 * xft
      dse(1) = dse(1) + num * rk * dseNk * abundelec(1)
      dse(3) = dse(3) + num * rk * (dseNkf*abundelec(1)+
     &dseNk*abundelec(3))
      dse(2) = dse(2) + num * rk * (dseNkT*abundelec(1)+
     &dseNk*abundelec(2))
      dse(5) = dse(5) + num * rk * (dseNkff*abundelec(1)+
     &2.0d0*dseNkf*abundelec(3)+dseNk*abundelec(5))
      dse(4) = dse(4) + num * rk * (dseNkTT*abundelec(1)+
     &2.0d0*dseNkT*abundelec(2)+dseNk*abundelec(4))
      dse(6) = dse(6) + num * rk * (dseNkfT*abundelec(1)+
     &dseNkT*abundelec(3)+dseNkf*abundelec(2)+dseNk*abundelec(6))
      END

************************************************************************
* Non-perfect correction for internal energy                           *
* Inputs: g and derivatives in fitxy, derivatives of x=nemh: xt,xf...  *
*  intermediates quantities already calculated (dmut1,...)             *
* Ouputs: Due(deriv2) and derivatives d()/dln                          *
*  WARNING: Due is  not initialized but modified. If num=1 this is an  *
*           addition and if num=-1, this is an difference.             *
************************************************************************

      SUBROUTINE duecalc_sod(num,fitxy,RT,abundelec,xf,xT,due)
*     ----------declarations----------
      IMPLICIT NONE
*     inputs
      include 'eoscom.math'
      DOUBLE PRECISION fitxy(deriv3),RT,abundelec(deriv2),xf,xT
      INTEGER num
*     ouputs
      DOUBLE PRECISION due(deriv2)
*     Intermediate variables
      DOUBLE PRECISION du,duf,duT,du1
*     --------------------
      du = -fitxy(3)
      duT = fitxy(5) - xT * fitxy(6)
      duf = -fitxy(6) * xf
      du1 = du * abundelec(1) * RT
      due(1) = due(1) + num * du1
      due(2) = due(2) + num * (RT*(du*abundelec(2)+abundelec(1)*duT)+
     &du1)
      due(3) = due(3) + num * (RT*(du*abundelec(3)+abundelec(1)*duf))
      END

************************************************************************
* Subroutine for electronic component calculus using polynomial        *
*  approximations of Fermi-Dirac integrals                             *
*                                                                      *
*  Input : fp,Tnodim and eta,rootf,halff1,rootfinv,                    *
*          fadd1inv (to avoid a second calculus)                       *
*  Output: n / n = 8pi * Ine / (Compton wavelength)**3                 *
*          p / p = 8pi * me * c * c * Ipe / (Compton wavelength)**3    *
*          s / s = se / (k*ne)                                         *
*          u / u = 8pi * me * c * c * Iue / (Compton wavelength)**3    *
*          nt,nf,pt,pf are first order derivatives of form dlnx/lny    *
*          st,sf are first order derivatives of form dx/dlny           *
*          All the second order derivatives are of the form:           *
*            d(first order)/dlny                                       *
* Rem: iue-->1.5d0 if eta-->-inf                                       *
************************************************************************

      SUBROUTINE electronic_sod(fp,Tnodim,eta,rootfinv,
     &fadd1inv,ine,ipe,ise,iue)
*     ---------declarations-----------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.math'
      include 'eoscom.ion'
      include 'eoscom.par'
      include 'eoscom.cons'
      DOUBLE PRECISION fp,eta(deriv2),Tnodim
      DOUBLE PRECISION fadd1inv,rootfinv
*     outputs
      DOUBLE PRECISION ine(deriv2),ipe(deriv2),ise(deriv2),iue(deriv2)
*     degree of the polynomial approximation
      INTEGER degree
      PARAMETER (degree=4)
*     compteurs
      INTEGER i,j,k,l
*     intermediate variables for approximating the FD
      DOUBLE PRECISION vf,vg,uf,ug,fac,p,gdpodg,gp,dlnpodlng
      DOUBLE PRECISION fdpodf,ggdpodgg,fgdpodgf,poinv,dlnpodlnf
      DOUBLE PRECISION d(3,6),ff(degree),gg(degree)
      DOUBLE PRECISION qett,qeft,qe,qef,qeff,qet,truncinv
      DOUBLE PRECISION gpog,ratio,ineinv,po,fpof,ffdpodff
*     Coefficients of the polynomials
      DOUBLE PRECISION fitfd(degree*degree*3)
      DATA FITFD/ 2.315472d0, 7.128660d0, 7.504998d0, 2.665350d0,
     & 7.837752d0, 23.507934d0, 23.311317d0, 7.987465d0, 9.215560d0,
     & 26.834068d0, 25.082745d0, 8.020509d0, 3.693280d0, 10.333176d0,
     & 9.168960d0, 2.668248d0,
     & 2.315472d0, 6.748104d0, 6.564912d0, 2.132280d0, 7.837752d0,
     & 21.439740d0, 19.080088d0, 5.478100d0, 9.215560d0, 23.551504d0,
     & 19.015888d0, 4.679944d0, 3.693280d0, 8.859868d0, 6.500712d0,
     & 1.334124d0,
     & 1.157736d0, 3.770676d0, 4.015224d0, 1.402284d0, 8.283420d0,
     & 26.184486d0, 28.211372d0, 10.310306d0, 14.755480d0,
     & 45.031658d0, 46.909420d0, 16.633242d0, 7.386560d0, 22.159680d0,
     & 22.438048d0, 7.664928d0/
*     --------------------

*        ---------affectations-----------
      gp = Tnodim * eta(3)
      vf = fadd1inv
      vg = 1.0d0 / (1.0d0+gp)
      uf = 2.0d0 * eta(5)
      ug = gp * vg
      fac = gp * (gp+1.0d0)
      fac = uf * fac * DSQRT(fac)
      ff(1) = 1.0d0
      gg(1) = 1.0d0
      DO i = 2, degree
         ff(I) = fp * ff(i-1)
         gg(I) = gp * gg(i-1)
         fac = fac * vf * vg
      END DO
      l = 1
*        ----------Ine calculus----------
      DO i = 1,3
         po = 0.0d0
         gdpodg = 0.0d0
         fdpodf = 0.0d0
         ffdpodff = 0.0d0
         ggdpodgg = 0.0d0
         fgdpodgf = 0.0d0
         DO j = 1,degree
            DO k = 1,degree
               p = fitfd(l) * gg(j) * ff(k)
               po = po + p
               gpog = (j-1) * p
               fpof = (k-1) * p
               gdpodg = gdpodg + gpog
               fdpodf = fdpodf + fpof
               ggdpodgg = ggdpodgg + gpog * (j-2)
               ffdpodff = ffdpodff + fpof * (k-2)
               fgdpodgf = fgdpodgf + gpog * (k-1)
               l = l + 1
            END DO
         END DO
         poinv = truncinv(po)
         dlnpodlng = gdpodg * poinv
         dlnpodlnf = fdpodf * poinv
         d(i,1) = fac * po
         d(i,2) = dlnpodlng + 1.5d0 + (2.5d0-degree) * ug
         d(i,3) = dlnpodlnf + 1.0d0 + (0.5d0*d(i,2)-degree) * uf
         po = truncinv(gdpodg)
         d(i,4) = dlnpodlng * (1.0d0+ggdpodgg*po-dlnpodlng) +
     &        (2.5d0-degree) * ug * vg
         d(i,6) = eta(5) * d(i,4) +
     &        (fgdpodgf*poinv-dlnpodlng*dlnpodlnf)
         d(i,5) = dlnpodlnf * (1.0d0-dlnpodlnf-eta(5)*dlnpodlng) +
     &        ffdpodff * poinv + fgdpodgf * poinv * eta(5) +
     &        uf * vf * (0.5d0*d(i,2)-degree) + eta(5) * d(i,6)
      END DO
*        ----------ne----------
      ine(1) = d(1,1)
      ine(2) = d(1,2)
      ine(3) = d(1,3)
      ine(4) = d(1,4)
      ine(5) = d(1,5)
      ine(6) = d(1,6)
*        ----------pe----------
      ipe(1) = gp * d(2,1)
      ipe(2) = 1.0d0 + d(2,2)
      ipe(3) = eta(5) + d(2,3)
      ipe(4) = d(2,4)
      ipe(5) = eta(5) * vf + d(2,5)
      ipe(6) = d(2,6)
*        ----------se----------
      ineinv = truncinv(ine(1))
      qe = d(3,1) * rootfinv * ineinv
      qet = d(3,2) - d(1,2)
      qef = d(3,3) - d(1,3) - eta(5)
      qett = d(3,4) - d(1,4)
      qeff = d(3,5) - d(1,5) - eta(5) * vf
      qeft = d(3,6) - d(1,6)
      ise(1) = qe + 2.0d0 * eta(3) - eta(1)
      ise(2) = qe * qet
      ise(3) = qe * qef - rootfinv
      ise(4) = qe * (qet*qet+qett)
      ise(5) = qe * (qef*qef+qeff) + eta(5) * rootfinv
      ise(6) = qe * (qef*qet+qeft)
*        ----------ue----------
      ratio = ipe(1) * ineinv / Tnodim
      iue(1) = ise(1) + eta(1) - ratio
      iue(2) = ise(2) + ratio * (ine(2)-ipe(2)+1.0d0)
      iue(3) = ise(3) + ratio * (ine(3)-ipe(3)) + eta(3)
      END

************************************************************************
*  INITIALIZATIONS                                                     *
************************************************************************
      BLOCK DATA eosinit
*     --------------------
      IMPLICIT NONE
      include 'evolpar.star'
      include 'eoscom.ion'
      include 'eoscom.phys'
*
* Energy levels for H, H+, H- and H2
*
      DATA enerlevelH /-13.598d0, .0d0, -14.354d0, -31.677/
*
* Energy levels for He,He+,He++,Li,Li+,...,Ca20+
*
      DATA enerlevelZ / -79.003d0, -54.416d0, .0d0,
     &-203.481d0, -198.089d0, -122.451d0, .0d0,
     &-399.139d0, -389.817d0, -371.606d0, -217.713d0, .0d0,
     &-670.969d0, -662.671d0, -637.516d0, -599.586d0, -340.22d0,
     &  0.0d0,
     &-1030.082d0, -1018.822d0, -994.439d0, -946.552d0, -882.06d0,
     &  -489.98d0, .0d0,
     &-1486.035d0, -1471.501d0, -1441.9d0, -1394.452d0, -1316.98d0,
     &  -1219.09d0, -667.03d0, .0d0,
     &-2043.812d0, -2030.194d0, -1995.077d0, -1940.143d0, -1862.73d0,
     &  -1748.83d0, -1610.71d0, -871.39d0, .0d0,
     &-2715.807d0, -2698.385d0, -2663.415d0, -2600.708d0, -2513.57d0,
     &  -2399.33d0, -2242.17d0, -2056.99d0, -1103.1d0, .0d0,
     &-3511.576d0, -3490.012d0, -3449.05d0, -3385.6d0, -3288.49d0,
     &  -3162.28d0, -3004.35d0, -2797.09d0, -2558.d0, -1362.2d0, .0d0,
     &-4419.895d0, -4414.756d0, -4367.47d0, -4295.83d0, -4196.92d0,
     &  -4058.52d0, -3886.37d0, -3677.89d0, -3413.7d0, -3113.8d0,
     &  -1648.7d0, .0d0,
     &-5451.084d0, -5443.438d0, -5428.403d0, -5348.26d0, -5238.95d0,
     &  -5097.68d0, -4911.17d0, -4686.22d0, -4420.3d0, -4092.3d0,
     &  -3724.8d0, -1963.d0, .0d0,
     &-6604.3d0, -6598.314d0, -6579.488d0, -6551.04d0, -6431.05d0,
     &  -6277.3d0, -6086.83d0, -5845.39d0, -5560.8d0, -5230.6d0,
     &  -4832.d0, -4390.d0, -2304.d0, .0d0,
     &-7887.229d0, -7879.078d0, -7862.733d0, -7829.241d0, -7784.1d0,
     &  -7617.33d0, -7412.25d0, -7165.76d0, -6862.6d0, -6511.5d0,
     &  -6110.1d0, -5634.d0, -5111.d0, -2673.d0, .0d0,
     &-9305.531d0, -9295.045d0, -9275.32d0, -9245.14d0, -9193.72d0,
     &  -9128.7d0, -8908.25d0, -8644.97d0, -8335.6d0, -7963.9d0,
     &  -7539.5d0, -7060.d0, -6499.d0, -5887.d0, -3070.d0, .0d0,
     &-10857.79d0, -10847.43d0, -10824.1d0, -10789.27d0, -10741.97d0,
     &  -10669.29d0, -10581.24d0, -10301.23d0, -9972.9d0, -9593.8d0,
     &  -9146.7d0, -8642.d0, -8077.d0, -7425.d0, -6718.d0, -3494.d0,
     &  .0d0,
     &-12554.437d0, -12541.47d0, -12517.66d0, -12478.05d0, -12424.59d0,
     &  -12356.89d0, -12259.86d0, -12145.67d0, -11797.3d0, -11396.9d0,
     &  -10941.3d0, -10412.d0, -9820.d0, -9163.d0, -8413.d0, -7604.d0,
     &  -3946.d0, .0d0,
     &-14398.338d0, -14382.579d0, -14354.95d0, -14314.21d0, -14254.4d0,
     &  -14179.36d0, -14088.35d0, -13963.95d0, -13820.5d0, -13397.9d0,
     &  -12919.d0, -12380.d0, -11762.d0, -11076.d0, -10320.d0, -9465.d0,
     &  -8547.d0, -4426.d0, .0d0,
     &-16380.651d0, -16376.31d0, -16344.68d0, -16298.96d0, -16238.04d0,
     &  -16155.38d0, -16055.48d0, -15737.78d0, -15782.8d0, -15607.d0,
     &  -15103.4d0, -14539.d0, -13910.d0, -13196.d0, -12409.d0,
     &  -11547.d0, -10579.d0, -9545.d0, -4934.d0, .0d0,
     &-18507.954d0, -18501.841d0, -18489.97d0, -18439.06d0, -18371.91d0,
     & -18287.48d0, -18178.7d0, -18051.d0, -17903.6d0, -17714.9d0,
     &  -17503.6d0, -16912.d0, -16255.d0, -15529.d0, -14712.d0,
     &  -13817.d0, -12843.d0, -11756.d0, -10599.d0, -5470.d0, .0d0/
*
* Degeneracy of the fundamental state for H, H+, H-, H2
*
      DATA lng0H /0.6931471810d0, 0.0d0, 0.0d0, 0.0d0/
*
* Degeneracy of the fundamental state for He,He+,He++,Li,Li+,...,Ca20+
*
      DATA lng0Z /0.0d0, 0.6931471810d0, 0.0d0,
     &0.6931471810d0, 0.0d0, 0.6931471810d0, 0.0d0,
     &0.0d0, 0.6931471810d0, 0.0d0, 0.6931471810d0, 0.0d0,
     &1.791759470d0, 0.0d0, 0.6931471810d0, 0.0d0, 0.6931471810d0,
     &0.0d0,
     &2.197224580d0,1.791759470d0, 0.0d0, 0.6931471810d0, 0.0d0,
     &0.6931471810d0, 0.0d0,
     &1.386294370d0, 2.197224580d0, 1.791759470d0, 0.0d0,
     &0.6931471810d0, 0.0d0, 0.6931471810d0, 0.0d0,
     &2.197224580d0, 1.386294370d0, 2.197224580d0, 1.791759470d0,
     &0.0d0, 0.6931471810d0, 0.0d0, 0.6931471810d0, 0.0d0,
     &1.791759470d0,2.197224580d0, 1.386294370d0, 2.197224580d0,
     &1.791759470d0, 0.0d0, 0.6931471810d0, 0.0d0, 0.6931471810d0,
     &0.0d0,
     &0.0d0, 1.791759470d0, 2.197224580d0, 1.386294370d0, 2.19722458d0,
     &1.791759470d0, 0.0d0, 0.6931471810d0, 0.0d0, 0.6931471810d0,
     &0.0d0,
     &0.6931471810d0, 0.0d0, 1.791759470d0, 2.197224580d0, 1.38629437d0,
     &2.197224580d0, 1.791759470d0, 0.0d0, 0.6931471810d0, 0.0d0,
     &0.6931471810d0, 0.0d0,
     &0.0d0, 0.6931471810d0, 0.0d0, 1.791759470d0, 2.197224580d0,
     &1.386294370d0, 2.197224580d0, 1.791759470d0, 0.0d0, 0.693147181d0,
     &0.0d0, 0.6931471810d0, 0.0d0,
     &1.791759470d0, 0.0d0, 0.6931471810d0, 0.0d0, 1.791759470d0,
     &2.197224580d0, 1.386294370d0, 2.197224580d0, 1.791759470d0, 0.0d0,
     &0.693147181d0, 0.0d0, 0.6931471810d0, 0.0d0,
     &2.197224580d0, 1.791759470d0, 0.0d0, 0.6931471810d0, 0.0d0,
     &1.791759470d0, 2.197224580d0, 1.386294370d0, 2.197224580d0,
     &1.791759470d0, 0.0d0, 0.693147181d0, 0.0d0, 0.6931471810d0, 0.0d0,
     &1.386294370d0, 2.197224580d0, 1.791759470d0, 0.0d0,
     &0.6931471810d0, 0.0d0, 1.791759470d0, 2.197224580d0, 1.38629437d0,
     &2.197224580d0, 1.791759470d0, 0.0d0, 0.693147181d0, 0.0d0,
     &0.6931471810d0, 0.0d0,
     &2.197224580d0, 1.386294370d0, 2.197224580d0, 1.791759470d0, 0.0d0,
     &0.6931471810d0, 0.0d0, 1.791759470d0, 2.197224580d0, 1.38629437d0,
     &2.197224580d0, 1.791759470d0, 0.0d0, 0.693147181d0, 0.0d0,
     &0.6931471810d0, 0.0d0,
     &1.791759470d0, 2.197224580d0, 1.386294370d0, 2.197224580d0,
     &1.791759470d0, 0.0d0, 0.6931471810d0, 0.0d0, 1.791759470d0,
     &2.197224580d0, 1.38629437d0, 2.197224580d0, 1.791759470d0, 0.0d0,
     &0.693147181d0, 0.0d0, 0.6931471810d0, 0.0d0,
     &0.0d0, 1.791759470d0, 2.197224580d0, 1.386294370d0, 2.197224580d0,
     &1.791759470d0, 0.0d0, 0.6931471810d0, 0.0d0, 1.791759470d0,
     &2.197224580d0, 1.38629437d0, 2.197224580d0, 1.791759470d0, 0.0d0,
     &0.693147181d0, 0.0d0, 0.6931471810d0, 0.0d0,
     &0.693147181d0, 0.0d0, 1.791759470d0, 2.197224580d0, 1.38629437d0,
     & 2.197224580d0, 1.791759470d0, 0.0d0, 0.6931471810d0, 0.0d0,
     & 1.791759470d0, 2.197224580d0, 1.38629437d0, 2.197224580d0,
     &1.791759470d0, 0.0d0, 0.693147181d0, 0.0d0, 0.693147181d0, 0.0d0,
     &0.0d0, 0.693147181d0, 0.0d0, 1.791759470d0, 2.197224580d0,
     &1.38629437d0, 2.197224580d0, 1.791759470d0, 0.0d0, 0.693147181d0,
     &0.0d0,1.791759470d0, 2.197224580d0, 1.38629437d0, 2.197224580d0,
     &1.791759470d0, 0.0d0, 0.693147181d0, 0.0d0, 0.693147181d0, 0.0d0/
      END

************************************************************************
* Say "yes" if x=y numerically                                         *
************************************************************************

      LOGICAL FUNCTION egal(x,y)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.cons'
      DOUBLE PRECISION x,y
*     --------------------
      IF ((x.EQ.0.d0).AND.(y.EQ.0.d0)) THEN
         egal = .TRUE.
      ELSE
         IF ((x.EQ.0.d0).OR.(y.EQ.0.d0)) THEN
            egal = ABS(x-y).LT.zero
         ELSE
            IF (x.GT.y) THEN
               egal = ABS((ABS(y/x)-1.0d0)).LT.numericalimit
            ELSE
               egal = ABS((ABS(x/y)-1.0d0)).LT.numericalimit
            END IF
         END IF
      END IF
      END

************************************************************************
* Compute effective eta (taking into account dktmu) and derivatives by *
*  mean of non-perfect correction for chemical potential               *
* INPUTS:f,rootf,eta                                                   *
* OUTPUTS: ratio,etaeff,etaeffT,etaefff...                             *
************************************************************************

      SUBROUTINE neweta_sod(fitxy,x,dmut1,dmut2,dmut1x,
     &dmut1y,dmut2x,fitnpi,eta,etaeff,etanpi)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.cons'
      include 'eoscom.math'
      DOUBLE PRECISION eta(deriv2),x(deriv2),fitxy(deriv3)
      DOUBLE PRECISION fitnpi(deriv3),dmut1,dmut2,dmut1x,dmut1y,dmut2x
*     Intermediate variables
      DOUBLE PRECISION dktmu,dktmuf,dktmut,dktmuff,dktmutt,dktmuft
      DOUBLE PRECISION bmunpi,bmunpif,bmunpit,bmunpiff,bmunpitt,bmunpift
      DOUBLE PRECISION dmu1,dmu2,dmu1x,dmu1y,dmu2x
*     outputs
      DOUBLE PRECISION etanpi(deriv2),etaeff(deriv2)
*     --------------------
      dktmu = -(fitxy(1)+fitxy(2))
      dktmuT = dmut1 - x(2) * dmut2
      dktmuf = -x(3) * dmut2
      dktmuff = -x(5) * dmut2 - x(3) * x(3) * dmut2x
      dktmuTT = 2.0d0 * dmut1x * x(2) - dmut1y - x(4) * dmut2 -
     &x(2) * x(2) * dmut2x
      dktmufT = x(3) * dmut1x - x(6) * dmut2 - x(3) * x(2) * dmut2x
*     --------------------
      dmu1 = fitnpi(3) + fitnpi(6)
      dmu2 = fitnpi(2) + fitnpi(4)
      dmu2x = fitnpi(4) + fitnpi(7)
      dmu1x = fitnpi(6) + fitnpi(9)
      dmu1y = fitnpi(5) + fitnpi(10)
      bmunpi = -(fitnpi(1)+fitnpi(2))
      bmunpiT = dmu1 - x(2) * dmu2
      bmunpif = -x(3) * dmu2
      bmunpiff = -x(5) * dmu2 - x(3) * x(3) * dmu2x
      bmunpiTT = 2.0d0 * dmu1x * x(2) - dmu1y - x(4) * dmu2 -
     &x(2) * x(2) * dmu2x
      bmunpifT = x(3) * dmu1x - x(6) * dmu2 - x(3) * x(2) * dmu2x
*     --------------------
      etaeff(1) = eta(1) + dktmu
      etaeff(3) = eta(3) + dktmuf
      etaeff(2) = dktmut
      etaeff(5) = eta(3) * eta(5) + dktmuff
      etaeff(4) = dktmutt
      etaeff(6) = dktmuft
*     --------------------
      etanpi(1) = eta(1) + bmunpi
      etanpi(3) = eta(3) + bmunpif
      etanpi(2) = bmunpit
      etanpi(5) = eta(3) * eta(5) + bmunpiff
      etanpi(4) = bmunpitt
      etanpi(6) = bmunpift
      END

************************************************************************
* Degeneracy and derivatives computation from f approximation for      *
*  non-perfect corrections                                             *
* Inputs : rootfap,fapx,term,sum,lnfap from SUBROUTINE fapprox         *
* Outputs: etap,etapx,etapxx,etapxxx                                   *
************************************************************************

      SUBROUTINE etapprox_sod(rootfap,fapx,term,sum,lnfap,etap,
     &etapx,etapxx,etapxxx)
*     ----------declarations----------
      IMPLICIT NONE
*     inputs
      DOUBLE PRECISION rootfap,fapx,term,sum,lnfap
*     ouputs
      DOUBLE PRECISION etap,etapx,etapxx,etapxxx
*     Intermediate variables
      DOUBLE PRECISION degen
*     --------------------
      etap = degen(lnfap,rootfap)
      etapx = fapx * rootfap
      etapxx = -etapx * term
      etapxxx = etapx * (term*term+sum)
*     --------------------
      END

************************************************************************
* f approximation for non-perfect effects - correct when f>>1          *
* Inputs : x=ne*mh y=13.6/kT                                           *
* Outputs: rootfap=SQRT(1+fap),lnfap=LN(fap),fapx=d(fap)/dlnx          *
*          term=fapx*(dlnsqrt(1+fap)/dlnf)-dlnfap/dlnx                 *
*          sum=
************************************************************************

      SUBROUTINE fapprox_sod(x,y,rootfap,fapx,term,sum,sumx,lnfap)
*     ----------declarations----------
      IMPLICIT NONE
*     inputs
      DOUBLE PRECISION x,y
      include 'evolpar.star'
      include 'eoscom.cons'
*     ouputs
      DOUBLE PRECISION rootfap,fapx,term,sum,lnfap,sumx
*     Intermediate variables
      DOUBLE PRECISION axy,bxy,add1,add2,bxy3,add1inv,produ
      DOUBLE PRECISION fap,fapxxx,fapxx,f1,f1inv,ratio,ratiox
      DOUBLE PRECISION a,b,truncln,bxyinv,ratioxx,fapxxxx,truncinv
*     f approximation fit coefficients
      PARAMETER (a=1.076540d0)
      PARAMETER (b=0.5695560d0)
*     --------------------
      axy = a * x * y**1.5d0
      bxy = b * axy
      bxy3 = bxy * inv3
      add1 = 1.0d0 + bxy
      add1inv = 1.0d0 / add1
      add2 = 1.0d0 + 4.0d0 * bxy3
*     ----------approximated f computation----------
      fap = axy * add1**inv3
      fapx = add2 * add1inv
      fapxx = bxy3 * add1inv / add2
      bxyinv = truncinv(bxy)
      fapxxx = fapxx * (3.0d0*bxyinv-4.0d0*bxy)
      fapxxxx = fapxxx - (4.0d0*bxy+3.0d0*bxyinv) * fapxx *
     &truncinv(fapxxx)
*     ----------for next eta calculus-----------------------
      f1 = fap + 1.0d0
      f1inv = 1.0d0 / f1
      rootfap = DSQRT(f1)
      lnfap = truncln(fap)
*     ratio=-dln(rootfap)/dlnf
      ratio = -0.5d0 * fap * f1inv
      ratiox = fapx * f1inv
      produ = ratiox + fapxx
      term = ratio * fapx - fapxx
      sum = fapxx * fapxxx - ratio * fapx * produ
      ratioxx = ratiox * (fapxx-fap*ratiox)
      sumx = fapxx * fapxxx * (fapxxx+fapxxxx) - fapx * ratio *
     &(produ*produ+ratioxx+fapxx*fapxxx)
      END


************************************************************************
* Pour les calculs avec diffusion                                      *
************************************************************************

      SUBROUTINE fordiffusion_sod(abundH,abundZ,rationsh,Y)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'evolcom.nuc'
      include 'eoscom.ion'
      include 'eoscom.par'
      include 'eoscom.math'
      DOUBLE PRECISION abundH(nbionHmax,deriv2),abundZ(nbionZmax,deriv2)
      DOUBLE PRECISION Y(nbz)
*     outputs
      DOUBLE PRECISION rationsh(nsp,zionmax)
*     Intermediate variables
      INTEGER iso,qeff,debut,rang,zz
*     --------------------
      debut = 1
      DO iso=2,3
         DO qeff=1,INT(znuc(iso))
            rang = debut + qeff
            rationsh(iso,qeff) = abundH(rang,1) / Y(izsp(iso))
         END DO
      END DO
      debut = 1
      DO iso=4,nsp
         zz = INT(znuc(iso))
         IF (ionized(zz)) THEN
            DO qeff=1,zz
               rang = debut + qeff
               rationsh(iso,qeff) = abundZ(rang,1) / Y(izsp(iso))
            END DO
            debut = debut + zz + 1
         ELSE
            DO qeff=1,zz-1
               rationsh(iso,qeff) = 0.0d0
            END DO
            IF (totioniz) THEN
               rationsh(iso,zz) = 1.0d0
            ELSE
               rationsh(iso,zz) = 0.0d0
            END IF
         END IF
      END DO
      END

************************************************************************
* Calculus of the number of free electrons due to elements heavier     *
*  than hydrogen                                                       *
*  INPUT: Y,nsh,nbz,shell,nbzion,nbion,zion,abundZ                     *
*  OUTPUT: heavelec(1) => mh*ne-, heavelec(2) => mh * dne/dlnT         *
*       heavelec(3)=df,heavelec(4)=dT2,heavelec(5)=df2,heavelec(6)=dfT *
************************************************************************

      SUBROUTINE freelectronZ_sod(Y,abundZ,heavelec,zmean,dheavelecdY)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.ion'
      include 'eoscom.par'
      include 'eoscom.math'
      include 'eoscom.cons'
      DOUBLE PRECISION Y(nbz),abundZ(nbionZmax,deriv2),zmean(nsp)
*     Intermediate variables
      INTEGER i,rang,debut,qeff,iso
      DOUBLE PRECISION nbelec(nbz),qeffreal
*     ouputs
      DOUBLE PRECISION heavelec(deriv2),dheavelecdY(nbz)
*     --------------------
      DO i=1,deriv2
         heavelec(i) = 0.0d0
      END DO
      rang = 0
      debut = 1
      dheavelecdY(1) = 0.0d0
*
*     loop on ionization-computed species (without H and heavy)
*     DO NOT TAKE "HEAVY" INTO ACCOUNT IN Z AVERAGE
*     (there is no good Z for it!)
*      ==> nbz-1
*
ccc...............
      nbelec(1) = 0.d0
      nbelec(nbz) = 0.d0
c      print *, 'nbz',nbz
      DO i=2,nbz-1      ! nbzmax ??
         dheavelecdY(i) = 0.0d0
         nbelec(i) = 0.0d0
         IF (ionized(i)) THEN
*           loop on charged ions p. 97 DUFOUR
            DO qeff=1,Z(i)
               qeffreal = DBLE(qeff)
               rang = debut + qeff
               nbelec(i) = nbelec(i) + qeffreal * abundZ(rang,1)
               heavelec(2) = heavelec(2) + qeffreal * abundZ(rang,2)
               heavelec(3) = heavelec(3) + qeffreal * abundZ(rang,3)
               heavelec(4) = heavelec(4) + qeffreal * abundZ(rang,4)
               heavelec(5) = heavelec(5) + qeffreal * abundZ(rang,5)
               heavelec(6) = heavelec(6) + qeffreal * abundZ(rang,6)
            END DO
****

            debut = debut + Z(i) + 1
         ELSE
            IF (totioniz) THEN
               nbelec(i) = DBLE(Z(i)) * Y(i)
            ELSE
               nbelec(i) = 0.0d0
            END IF
         END IF
c     print *,'nbz',i,'heavelec(1)',heavelec(1),'nbelec(i)',nbelec(i)
c         print *,i
c         print *, 'heavelec(1) avant',heavelec(1)
         heavelec(1) = heavelec(1) + nbelec(i)
c         print *, 'heavelec(1) apres',heavelec(1)
         dheavelecdY(i) = nbelec(i) / Y(i)
      END DO
c      print *,'NE STAREVOL',heavelec(1)/amu
c      stop
      zmean(1) = 0.d0
      zmean(nsp-1) = 0.d0
      zmean(nsp) = 0.d0
      DO iso=2,nsp-2
         zmean(iso) = nbelec(izsp(iso)) / Y(izsp(iso))
      END DO
      END

************************************************************************
* Computation of the thermodynamics quantities of the gazeous component*
************************************************************************

      SUBROUTINE gaz_sod(Y,RT,betaev,lnbetaev,sumNi,lnparfunH,
     &lnparfunZ,abundH,abundZ,pgaz,sgaz,ugaz,rho,sumNinv)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'evolcom.cons'
      include 'eoscom.ion'
      include 'eoscom.phys'
      include 'eoscom.par'
      include 'eoscom.math'
      include 'eoscom.cons'
      DOUBLE PRECISION rho(deriv2)
      DOUBLE PRECISION lnparfunH(nbionHmax,deriv3)
      DOUBLE PRECISION lnparfunZ(nbionZmax,deriv3)
      DOUBLE PRECISION abundH(nbionHmax,deriv2),abundZ(nbionZmax,deriv2)
      DOUBLE PRECISION sumNi(deriv2),RT,betaev,lnbetaev,Y(nbz),sumNinv
*     Intermediate variables
      DOUBLE PRECISION c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,ratio1,ratio2,lnm
      DOUBLE PRECISION lnnQ,mrhoT,prhoT,khifgaz,khiTgaz
      INTEGER i,finA,finE,rangE,rangA,qeff,indx
      DOUBLE PRECISION xlnx(deriv2),truncln,coeff
*     ouputs
      DOUBLE PRECISION pgaz(deriv2),sgaz(deriv2),ugaz(deriv2)
*     --------------------
      coeff = 2.5d0 - truncln(rho(1)) - lnavn
      prhot = rho(2) + 1.0d0
      mrhoT = 1.5d0 - rho(2)
      ratio1 = sumNi(2) * sumNinv
      ratio2 = sumNi(3) * sumNinv
      khiTgaz = prhoT + ratio1
      khifgaz = rho(3) + ratio2
*     ---------- P ----------
      pgaz(1) = rho(1) * RT * sumNi(1)
      pgaz(2) = pgaz(1) * khiTgaz
      pgaz(3) = pgaz(1) * khifgaz
      pgaz(4) = pgaz(1) * (khiTgaz*khiTgaz+rho(4)+sumNi(4)*sumNinv-
     &ratio1*ratio1)
      pgaz(5) = pgaz(1) * (khifgaz*khifgaz+rho(5)+sumNi(5)*sumNinv-
     &ratio2*ratio2)
      pgaz(6) = pgaz(1) * (khiTgaz*khifgaz+rho(6)+sumNi(6)*sumNinv-
     &ratio1*ratio2)
*     ---------- S & U ----------
      sgaz(1) = coeff * sumNi(1)
      sgaz(2) = coeff * sumNi(2) + mrhoT * sumNi(1)
      sgaz(3) = coeff * sumNi(3) - rho(3) * sumNi(1)
      sgaz(4) = coeff * sumNi(4) + 2.0d0 * sumNi(2) * mrhoT -
     &rho(4) * sumNi(1)
      sgaz(5) = coeff * sumNi(5) - 2.0d0 * sumNi(3) * rho(3) -
     &rho(5) * sumNi(1)
      sgaz(6) = coeff * sumNi(6) - sumNi(2) * rho(3) + sumNi(3) * mrhoT-
     &rho(6) * sumNi(1)
      ugaz(1) = 1.5d0 * sumNi(1)
      ugaz(2) = 1.5d0 * (sumNi(2)+sumNi(1))
      ugaz(3) = 1.5d0 * sumNi(3)
      ugaz(4) = 0.d0
      ugaz(5) = 0.d0
      ugaz(6) = 0.d0
      rangA = 0
      finA = 0
*     ---------- elements Z>1, without "heavy" ----------
      DO i=2,nbz-1
         lnm = lnmZ(i-1) + lnamu
         finE = indx(2,Z(i)-1)
         IF (ionized(i)) THEN
            DO qeff=1,Z(i)+1
               rangA = finA + qeff
               rangE = finE + qeff
               c0 = lnparfunZ(rangA,2)
               c9 = lnparfunZ(rangA,4)
               c1 = lnparfunZ(rangA,1) + c0 + lnnQ(lnbetaev,lnm)
               CALL loga_sod(abundZ,nbionZmax,deriv2,rangA,1,2,3,4,5,6,
     &xlnx)
               c3 = c0 + c9
               c6 = c9 + lnparfunZ(rangA,7)
               c4 = betaev * enerlevelZ(rangE)
               c2 = lnparfunZ(rangA,6)
               c7 = lnparfunZ(rangA,3) + c2
               c8 = lnparfunZ(rangA,9) + c2
               c5 = lnparfunZ(rangA,10) + c2
               sgaz(1) = sgaz(1) + abundZ(rangA,1) * c1 - xlnx(1)
               sgaz(2) = sgaz(2) + abundZ(rangA,2) * c1 - xlnx(2) +
     &abundZ(rangA,1) * c3
               sgaz(3) = sgaz(3) + abundZ(rangA,3) * c1 - xlnx(3) +
     &abundZ(rangA,1) * c7
               sgaz(4) = sgaz(4) + abundZ(rangA,4) * c1 - xlnx(4) +
     &abundZ(rangA,2) * 2.0d0 * c3 + abundZ(rangA,1) * c6
               sgaz(5) = sgaz(5) + abundZ(rangA,5) * c1 - xlnx(5) +
     &abundZ(rangA,1) * c8 + 2.0d0 * abundZ(rangA,3) * c7
               sgaz(6) = sgaz(6) + abundZ(rangA,6) * c1 - xlnx(6) +
     &abundZ(rangA,3) * c3 + abundZ(rangA,1) * c5 + abundZ(rangA,2) * c7
               ugaz(1) = ugaz(1) + abundZ(rangA,1) * (c4+c0)
               ugaz(2) = ugaz(2) + abundZ(rangA,2) * (c4+c0) +
     &abundZ(rangA,1) * (c9+c0)
               ugaz(3) = ugaz(3) + abundZ(rangA,3) * (c4+c0) +
     &abundZ(rangA,1) * c2
            END DO
*
*     Substract a constant term changes the zero point of energy
*     and insure that it is nul in nonionised plasma
*
            ugaz(1) = ugaz(1) - Y(i) * enerlevelZ(finE+1)*betaev
            finA = finA + Z(i) + 1
         ELSE
            IF (totioniz) THEN
*
*     ---- 0.0=ln(Z) of these completly ionised elements (like H+) -----
*     ---------- Note that there is no contribution to u because -------
*     -----------the last level has a nul energy -----------------------
*
               c1 = 0.0d0 + lnnQ(lnbetaev,lnm)
               sgaz(1) = sgaz(1) + Y(i) * c1 - Y(i) * truncln(Y(i))
            ELSE
*
*     ---------- Note the contribution to u is nul because of the ------
*     -----------the choice of zero level in neutral plasma :-----------
*     -----------i.d. u contribution is betaev * Y(Z) * enerlevelZ(Z)---
*     ----------- and we substract the same term after------------------
*
               c1 = lng0Z(finE+1) + lnnQ(lnbetaev,lnm)
               sgaz(1) = sgaz(1) + Y(i) * c1 - Y(i) * truncln(Y(i))
               ugaz(1) = ugaz(1) + 0.0d0
            END IF
         END IF
      END DO
*     ---------- Heavy ----------
      lnm = lnmheavy + lnamu
      c1 = lng0heavy + lnnQ(lnbetaev,lnm)
      sgaz(1) = sgaz(1) + Y(nbz) * c1 - Y(nbz) * truncln(Y(nbz))
*     ---------- Hydrogene ----------
      DO i=1,nbionHmax
         lnm = lnmH(i) + lnamu
         c0 = lnparfunH(i,2)
         c9 = lnparfunH(i,4)
         c1 = lnparfunH(i,1) + c0 + lnnQ(lnbetaev,lnm)
         CALL loga_sod(abundH,nbionHmax,deriv2,i,1,2,3,4,5,6,xlnx)
         c3 = c0 + c9
         c4 = betaev * enerlevelH(i)
         c6 = c9 + lnparfunH(i,7)
         c2 = lnparfunH(i,6)
         c7 = lnparfunH(i,3) + c2
         c8 = lnparfunH(i,9) + c2
         c5 = lnparfunH(i,10) + c2
         sgaz(1) = sgaz(1) + abundH(i,1) * c1 - xlnx(1)
         sgaz(2) = sgaz(2) + abundH(i,2) * c1 - xlnx(2) + abundH(i,1)*c3
         sgaz(3) = sgaz(3) + abundH(i,3) * c1 - xlnx(3) + abundH(i,1)*c7
         sgaz(4) = sgaz(4) + abundH(i,4) * c1 - xlnx(4) +
     &abundH(i,2) * 2.0d0 * c3 + abundH(i,1) * c6
         sgaz(5) = sgaz(5) + abundH(i,5) * c1 - xlnx(5) +
     &abundH(i,1) * c8 + 2.0d0 * abundH(i,3) * c7
         sgaz(6) = sgaz(6) + abundH(i,6) * c1 - xlnx(6) +
     &abundH(i,3) * c3 + abundH(i,1) * c5 + abundH(i,2) * c7
         ugaz(1) = ugaz(1) + abundH(i,1) * (c4+c0)
         ugaz(2) = ugaz(2) + abundH(i,2) * (c4+c0) +
     &abundH(i,1) * (c9+c0)
         ugaz(3) = ugaz(3) + abundH(i,3) * (c4+c0) +
     &abundH(i,1) * c2
      END DO
      ugaz(1) = ugaz(1) - 0.5d0 * Y(1) * enerlevelH(4) * betaev
      DO i=1,deriv2
         sgaz(i) = sgaz(i) * rk
         ugaz(i) = ugaz(i) * RT
      END DO
      END

************************************************************************
* Determination of the H partition function                            *
*  Reference: Irwin ApJSS 45,621 (1981)                                *
************************************************************************

      SUBROUTINE h_sod(lnTsh,lnparfunH)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.ion'
      include 'eoscom.phys'
      include 'eoscom.cons'
      DOUBLE PRECISION lnTsh
*     outputs
      include 'eoscom.math'
      DOUBLE PRECISION lnparfunH(nbionHmax,deriv3)
*     Intermediate variables
      INTEGER i,j,k,dimfit
      PARAMETER (dimfit=6)
      PARAMETER (k=1)
      DOUBLE PRECISION d0(dimfit)
      DOUBLE PRECISION powerlnT(dimfit)
      DATA d0 /-2.61655891d2, 1.63428326d2, -4.06133526d1,
     &5.03282928d0, -3.10998364d-1, 7.66654594d-3/
*     --------------------
      powerlnT(1) = 1.0d0
      DO i=2,dimfit
         powerlnT(i) = lnTsh * powerlnT(i-1)
      END DO
      DO j=1,deriv3
         lnparfunH(k,j) = 0.0d0
      END DO
      DO j=1,dimfit
         lnparfunH(k,1) = lnparfunH(k,1) + d0(j) * powerlnT(j)
      END DO
      DO j=1,dimfit-1
         lnparfunH(k,2) = lnparfunH(k,2) + d0(j+1) * j * powerlnT(j)
      END DO
      DO j=1,dimfit-2
         lnparfunH(k,4) = lnparfunH(k,4) +
     &d0(j+2) * j * (j+1) * powerlnT(j)
      END DO
      DO j=1,dimfit-3
         lnparfunH(k,7) = lnparfunH(k,7) +
     &d0(j+3) * j * (j+1) * (j+2) * powerlnT(j)
      END DO
      END

************************************************************************
* Determination of the H2 partition function                           *
*  Reference: Irwin AA 182, 348 (1987)                                 *
************************************************************************

      SUBROUTINE h2Irwin_sod(lnTsh,lnparfunH)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.ion'
      include 'eoscom.phys'
      include 'eoscom.cons'
      DOUBLE PRECISION lnTsh
*     outputs
      include 'eoscom.math'
      DOUBLE PRECISION lnparfunH(nbionHmax,deriv3)
*     Intermediate variables
      INTEGER i,j,k,dimfit
      PARAMETER (dimfit=8)
      PARAMETER (k=4)
      DOUBLE PRECISION d0(dimfit)
      DOUBLE PRECISION powerlnT(dimfit)
      DATA d0 /1.67298118410d4, -1.49945289142d4, 5.74838863349d3,
     &-1.22210505066d3, 1.55637569965d2, -1.18744926193d1,
     &5.02617615447d-1, -9.10563051348d-3/
*     --------------------
      powerlnT(1) = 1.0d0
      DO i=2,dimfit
         powerlnT(i) = lnTsh * powerlnT(i-1)
      END DO
      DO j=1,deriv3
         lnparfunH(k,j) = 0.0d0
      END DO
      DO j=1,dimfit
         lnparfunH(k,1) = lnparfunH(k,1) + d0(j) * powerlnT(j)
      END DO
      DO j=1,dimfit-1
         lnparfunH(k,2) = lnparfunH(k,2) + d0(j+1) * j * powerlnT(j)
      END DO
      DO j=1,dimfit-2
         lnparfunH(k,4) = lnparfunH(k,4) +
     &d0(j+2) * j * (j+1) * powerlnT(j)
      END DO
      DO j=1,dimfit-3
         lnparfunH(k,7) = lnparfunH(k,7) +
     &d0(j+3) * j * (j+1) * (j+2) * powerlnT(j)
      END DO
      END

************************************************************************
* H2 partition function calculus (rotation+vibration)                  *
* It is described by a analytical function based on the pressure       *
* equilibrium constant Kp(H2).    (Vardya ApJS,4,281,1960)             *
* Input: beta=1/kT (ev-1)       Output: parfun(4,i)                    *
************************************************************************

      SUBROUTINE h2Vardya_sod(beta,lnparfunH)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.ion'
      include 'eoscom.phys'
      include 'eoscom.math'
      include 'eoscom.cons'
      DOUBLE PRECISION beta
*     Intermediate variables
      DOUBLE PRECISION dbeta,expmdb,dbeta2,beta2,sum
      DOUBLE PRECISION arg1,arg2,arg3,truncexpo,truncln
      DOUBLE PRECISION zheta,zhetaT,zhetaTT,zhetaTTT,inv
      DOUBLE PRECISION d0,d1,d22,d33,truncinv
*     ouputs
      DOUBLE PRECISION lnparfunH(nbionHmax,deriv3)
*     ----------affectations----------
      PARAMETER (d0=8.7961574d0)
      PARAMETER (d1=0.448d0)
      PARAMETER (d22=2.439844d-2)
      PARAMETER (d33=6.162951d-4)	
*     --------------------
      arg1 = d1 * beta
      beta2 = beta * beta
      arg2 = d22 * beta2
      arg3 = d33 * beta * beta2
      dbeta = dissH2 * beta
      dbeta2 = dbeta * dbeta
      expmdb = truncexpo(-dbeta)
*     --------------------
      zheta = 1.0d0 - (1.0d0+dbeta) * expmdb
      inv = truncinv(zheta)
      zhetaT = -dbeta2 * expmdb * inv
      sum = dbeta - 2.0d0 - zhetaT
      zhetaTT = sum * zhetaT
      zhetaTTT = zhetaTT * sum - zhetaT * (dbeta+zhetaTT)
*     --------------------
      lnparfunH(4,1) = d0 + truncln(zheta) - 2.5d0 * truncln(dbeta) +
     &arg1 - arg2 + arg3
      lnparfunH(4,2) = zhetaT + 2.5d0 - arg1 + 2.0d0 * arg2 -
     &3.0d0 * arg3
      lnparfunH(4,4) = zhetaTT + arg1 - 4.0d0 * arg2 + 9.0d0 * arg3
      lnparfunH(4,7) = zhetaTTT - arg1 + 8.0d0 * arg2 - 27.0d0 * arg3
      END

************************************************************************
* Determination of the He partition function                           *
*  Reference: Irwin ApJSS 45,621 (1981)                                *
************************************************************************

      SUBROUTINE he_sod(lnTsh,nZ,lnparfunZ)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.ion'
      include 'eoscom.phys'
      include 'eoscom.cons'
      DOUBLE PRECISION lnTsh
      INTEGER nZ
*     outputs
      include 'eoscom.math'
      DOUBLE PRECISION lnparfunZ(nbionZmax,deriv3)
*     Intermediate variables
      INTEGER i,j,k1,k2,k,dimfit
      PARAMETER (dimfit=6)
      DOUBLE PRECISION d0(dimfit),d1(dimfit),d2(dimfit)
      DOUBLE PRECISION powerlnT(dimfit)
      DATA d0 /-3.76575219d-1, 2.33951687d-1, -5.79755525d-2,
     &7.16333160d-3, -4.41302573d-4, 1.08442997d-5/
      DATA d1 /6.93147179d-1, 9.29636701d-10, -2.30049742d-10,
     &2.83829746d-11, -1.74590774d-12, 4.28355287d-14/
      DATA d2 /0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/
*     --------------------
      powerlnT(1) = 1.0d0
      DO i=2,dimfit
         powerlnT(i) = lnTsh * powerlnT(i-1)
      END DO
      k = nz + 1
      k1 = nz + 2
      k2 = nz + 3
      DO j=1,deriv3
         lnparfunZ(k,j) = 0.0d0
         lnparfunZ(k1,j) = 0.0d0
         lnparfunZ(k2,j) = 0.0d0
      END DO
      DO j=1,dimfit
         lnparfunZ(k,1) = lnparfunZ(k,1) + d0(j) * powerlnT(j)
         lnparfunZ(k1,1) = lnparfunZ(k1,1) + d1(j) * powerlnT(j)
         lnparfunZ(k2,1) = lnparfunZ(k2,1) + d2(j) * powerlnT(j)
      END DO
      DO j=1,dimfit-1
         lnparfunZ(k,2) = lnparfunZ(k,2) + d0(j+1) * j * powerlnT(j)
         lnparfunZ(k1,2) = lnparfunZ(k1,2) + d1(j+1) * j * powerlnT(j)
         lnparfunZ(k2,2) = lnparfunZ(k2,2) + d2(j+1) * J * powerlnT(j)
      END DO
      DO j=1,dimfit-2
         lnparfunZ(k,4) = lnparfunZ(k,4) +
     &d0(j+2) * j * (j+1) * powerlnT(j)
         lnparfunZ(k1,4) = lnparfunZ(k1,4) +
     &d1(j+2) * j * (j+1) * powerlnT(j)
         lnparfunZ(k2,4) = lnparfunZ(k2,4) +
     &d2(j+2) * j * (j+1) * powerlnT(j)
      END DO
      DO j=1,dimfit-3
         lnparfunZ(k,7) = lnparfunZ(k,7) +
     &d0(j+3) * j * (j+1) * (j+2) * powerlnT(j)
         lnparfunZ(k1,7) = lnparfunZ(k1,7) +
     &d1(j+3) * j * (j+1) * (j+2) * powerlnT(j)
         lnparfunZ(k2,7) = lnparfunZ(k2,7) +
     &d2(j+3) * j * (j+1) * (j+2) * powerlnT(j)
      END DO
      END

************************************************************************
*   Improve the real solutions of a four order polynomial with         *
*    Newton-Raphson iterations                                         *
************************************************************************

      SUBROUTINE IMPROVEREAL(a,b,c,d,e,x0,x)
*     --------------------
      IMPLICIT NONE
*     inputs
      DOUBLE PRECISION a,b,c,d,e,x0
*     Intermediate variables
      INTEGER itermax,iter
      DOUBLE PRECISION xold,curx,dx,eps
      DOUBLE PRECISION f,df,p1,p2,p3,p4,p5,y,po,ppo,dpo
      PARAMETER (itermax=20,eps=1.0d-14)
*     outputs
      DOUBLE PRECISION x
*     --------------------
      po (p1,p2,p3,p4,p5,y) = p1*y**4+p2*y**3+p3*y*y+p4*y+p5
      ppo (p2,p3,p4,p5,y) = p2*y**2*y+p3*y*y+p4*y+p5
      dpo (p1,p2,p3,p4,p5,y) = 4.0d0*p1*y**3+3.0d0*p2*y*y+2.0d0*p3*y+p4
*     --------------------
      xold = x0
      curx = xold
      dx = 1.d90
      iter = 1
      curx = min(curx,1.d20)
      df = dpo(a,b,c,d,e,curx)
      IF (df.NE.0.0d0.AND.curx.NE.0.0d0) THEN
         DO WHILE(dx.GE.eps.AND.iter.LE.itermax)
            if (a.eq.0.d0) then
               f = ppo(b,c,d,e,curx)
            else
               f = po(a,b,c,d,e,curx)
            endif
c           f = po(a,b,c,d,e,curx)
            df = dpo(a,b,c,d,e,curx)
            xold = curx
            curx = xold-f/df
            dx = abs(curx-xold)/abs(xold)
            iter = iter+1
         END DO
      END IF
      x = curx
      end

************************************************************************
* Compute the sum: SUM(Z+1) from a to b                                *
************************************************************************

      INTEGER FUNCTION indx(a,b)
*     --------------------
      IMPLICIT NONE
      INTEGER a,b
*     --------------------
      indx = INT(0.5*(b+1)*b) - INT(0.5*a*(a-1)) + b - a + 1
      END

************************************************************************
* Initialization of various constants                                  *
*    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                       *
*    !!!!!!!!!!!!!!! C G S !!!!!!!!!!!!!!!!!!!!!                       *
*    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                       *
************************************************************************

      SUBROUTINE initcons
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'evolcom.cons'
      include 'eoscom.ion'
      include 'eoscom.phys'
      include 'eoscom.par'
*     Intermediate variables
      INTEGER i
      DOUBLE PRECISION qe2
      PARAMETER (qe2 = 2.307079556d-19)
*     outputs
      include 'eoscom.cons'
*     ---------- Mathematical constant ----------
      inv3 = 1.0d0 / 3.0d0
      ln2 = DLOG(2.0d0)
      ln4 = DLOG(4.0d0)
      ln10 = DLOG(10.0d0)
      ln10inv = 1.0d0 / ln10
*     Warning: corXp-corXm <>0 !!!!!
      corlnfp = 0.1d0
      corlnfm = -0.1d0
      corTm = -0.1d0
      corTp = 0.1d0
      zero = 1.0d-90
      limexp = 227.0d0
      liminf = 99.0d0
      numericalimit = 1.0d-13
      physicalimit = 1.0d-8
      lnf0 = 100.d0
*     ---------- Values of physical constants ----------
      ergev = 1.60217733d-12
      masselec = 9.1093897d-28
      lnmasselec = DLOG(masselec)
      amu = 1.0d0 / avn
      lnamu = DLOG(amu)
      lnavn = DLOG(avn)
*     Ionization potential of H
      khiH = enerlevelH(2) - enerlevelH(1)
*     Dissociation energy for H2
      dissH2 = 2.0d0 * enerlevelH(1) - enerlevelH(4)
*     Compton wavelength
      lc = h / (masselec*c)
*     8*pi/(Compton wavelength)**3
      pi8l3 = pim8 / lc**3
*     8*pi*amu/(Compton wavelength)**3
      m8pl3 = pi8l3 * amu
*     8*pi*me*c*c/(Compton wavelength)**3
      mc28pl3 = pi8l3 * masselec * c * c
*     1/(me*c*c)
      mc2inv = 1.0d0 / (masselec*c*c)
*     constant for coulomb effect fitting
      gamcoeff = (pim4*inv3/amu)**inv3 * qe2 / (ergev*khiH)
*     constant for nQ computation
      lncoeffnQ = DLOG(pim2*ergev) - 2.0d0 * DLOG(h)
*     ---------- Approximation limit ----------
c     Tmaxirwin = 1.6d4
c     Tmaxvardya = 2.5d4
      Tmaxvardya = 2.5d6
      Tmaxirwin = 2.5d6
      IF (compH2) THEN
         IF (vardyaZH2) THEN
            TmaxH2 = Tmaxvardya
         ELSE
            TmaxH2 = Tmaxirwin
         END IF
      ELSE
         TmaxH2 = 0.0d0
      END IF
      END

************************************************************************
* Compute nbionZ, zion(i), nbzion that is to say the number of ions,   *
*  except H which are presents in the plasma, the atomic number of     *
*  chemical species whose ionisation is computed, and the number of    *
*  atomic number whose ionization is computed.                         *
* Inputs: ionized(i)=true if the ionization of the Z=i element must be *
*                    computed.                                         *
*          nbz=number of chemical elements whose nucleosynthesis is    *
*             computed in the code                                     *
************************************************************************

      SUBROUTINE initnbion
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'evolcom.nuc'
      include 'eoscom.ion'
      include 'eoscom.par'
*     Intermediate variables
      INTEGER i,nz,ziso,izmax,iso
      LOGICAL exist
      CHARACTER*25 text
*     ouputs
*     include 'eoscom.ion'
*     --------------------
      DO i=1,nbz
         Z(i) = 0
      END DO
      nz = 1
      izsp(1) = 0
      Z(1) = INT(znuc(nspmin))
      izsp(2) = 1
      izmax = 1
      DO iso=nspmin+1,nspmax
         ziso = INT(znuc(iso))
         exist = .false.
         DO i=1,nz
            exist = (exist.OR.(Z(i).EQ.ziso))
         END DO
         IF (.NOT.exist) THEN
            nz = nz + 1
            Z(nz) = ziso
            IF (Z(nz).GT.Z(izmax)) THEN
               izmax = nz
            END IF
         END IF
      END DO
      DO iso=nspmin+1,nspmax
         ziso = INT(znuc(iso))
         DO i=1,nbz
            IF (Z(i).EQ.INT(znuc(iso))) THEN
               izsp(iso) = i
            END IF
         END DO
      END DO
      IF (nz.NE.nbz) THEN
         WRITE (NOUT,*) 'EOS - initnbion - ',nz,' detected <> nbz !!!'
         STOP
      END IF
      IF (.NOT.ionized(1)) THEN
         WRITE (NOUT,*) 'EOS - initnbion - hydrogen not ionized !!!'
         STOP
      END IF
      IF (ionized(nbz)) THEN
         WRITE (NOUT,*) 'EOS - initnbion - heavy ionized !!!'
         STOP
      END IF
      IF (Z(1).NE.1) THEN
         WRITE (NOUT,*) 'EOS - initnbion - Z(1)<>1 !!!'
         STOP
      END IF
      IF (.NOT.totioniz) THEN
         IF (Z(izmax).GT.zionmax) THEN
            WRITE (NOUT,*) 'EOS - initnbion - Zmax>Zionmax & not ',
     &           'totioniz'
            STOP
         END IF
      END IF
      DO i=1,nbz
         IF (ionized(i)) THEN
            IF (i.EQ.1) THEN
               text = 'HYDROGEN ==> ionized'
            ELSE
               text = 'ionized'
            END IF
            IF (z(i).GT.zionmax) THEN
               WRITE (NOUT,*) 'EOS - initnbion - ionized>zionmax)',i,
     &              z(i)
               STOP
            END IF
         ELSE
            IF (i.EQ.nbz) THEN
               text = 'HEAVY ==> not ionized'
            ELSE
               IF (totioniz) THEN
                  text = 'always totaly ionized'
               ELSE
                  text = 'always neutral'
               END IF
            END IF
         END IF
         IF (checkcomposition) THEN
            WRITE (NOUT,'(A,I3,A,I3,A,A)') 'Z(',i,')=',Z(i),' ',text
         END IF
      END DO
      END

************************************************************************
* Initialization of some parameters                                    *
************************************************************************

      SUBROUTINE initpar
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'evolcom.eos'
      include 'eoscom.ion'
      include 'eoscom.par'

*     Intermediate variables
      INTEGER i
      CHARACTER *50 c
*     --------------------
*     ----------PARAMETER CARD FOR EOS----------
*     --------FERMI-DIRAC INTEGRALS--------
c      approxFermiDirac = .true.
*     ----------HYDROGEN SPECIES----------
      compH2 = addH2
      compHm = addHm
*     ----------IONIZATION----------
      totioniz = Ztotioni
      do i = 3,nbz-1        ! nbzmax
         ionized(i) = .true.
c$$$         ionized(i) = .false.         !ionisation totale
      enddo   
      ionized(1) = .true.
      ionized(2) = ionizHe
c$$$      ionized(3) = .false.
c$$$      ionized(4) = .false.
c$$$      ionized(5) = .false.
      ionized(6) = ionizC
      ionized(7) = ionizN
      ionized(8) = ionizO
c$$$      ionized(9) = .false.
c$$$      ionized(10) = ionizNe
c$$$      ionized(11) = .false.
c$$$      ionized(12) = .false.
c$$$      ionized(13) = .false.
      ionized(14) = .false.         ! Si
c$$$      ionized(15) = .false.
c$$$      ionized(16) = .false.
c$$$      ionized(17) = .false.
      ionized(18) = .false.
      ionized(19) = .false.
      ionized(20) = .false.
      ionized(21) = .false.
*     ----------NON IDEAL EFFECTS----------
      compionpres = addprio
      compcoulshiel = .true.
      fapprocoul = .false.
      nonideal = (compcoulshiel.OR.compionpres)
*     ----------PARTITION FUNCTIONS----------
      vardyaZH2 = .false.
      irwinZH2 = .true.
      compZH0 = addfitZ
      compZHe = addfitZ
      compZCNO = .false.
*     ----------ATOMIC MASSES----------
      compavmass = .true.
*     ----------DERIVATIVES----------
      eosderiv = numderiv
      approxsod = .false.
*     ----------MONITOR----------
      checkcomposition = .false.
      checkparacard = .true.
*     --------- END OF PARAMETER CARD -----------
      IF (.NOT.compcoulshiel) THEN
         fapprocoul = .FALSE.
      END IF
      IF (eosderiv) THEN
         approxsod = .FALSE.
      END IF
      IF (checkparacard) THEN
         WRITE (90,'(/,A,/)') ' PARAMETER CARD OF THE EOS'
         WRITE (90,'(A)') 'o FERMI-DIRAC INTEGRALS'
c         IF (approxFermiDirac) THEN
            c = 'yes'
c         ELSE
c            c = 'no'
c         END IF
         WRITE (90,'(A,A)') '  Fermi-Dirac integral approximations: ',c
         WRITE (90,'(A)') 'o HYDROGEN SPECIES'
         IF (compH2) THEN
            c = 'yes'
         ELSE
            c = 'no'
         END IF
         WRITE (90,'(A,A)') '  H2 computation: ',c
         IF (compHm) THEN
            c = 'yes'
         ELSE
            c = 'no'
         END IF
         WRITE (90,'(A,A)') '  H- computation: ',c
         WRITE (90,'(A)') 'o IONIZATION'
         IF (totioniz) THEN
            c = '  always totaly ionized'
         ELSE
            c = '  always neutral'
         END IF
         WRITE (90,'(A,A)') '  Non ionized elements are: ',c
         DO i=1,nbzmax
            IF (ionized(i)) THEN
               c = 'yes'
            ELSE
               c = 'no'
            END IF
            WRITE (90,'(3x,A,I2,A,A)') 'Element ',i,
     &           ' compute ionization: ',c
         END DO
         WRITE (90,'(A)') 'o NON IDEAL EFFECTS'
         IF (compionpres) THEN
            c = 'yes'
         ELSE
            c = 'no'
         END IF
         WRITE (90,'(A,A)') '  Add pressure ionization effects: ',c
         IF (compcoulshiel) THEN
            c = 'yes'
         ELSE
            c = 'no'
         END IF
         WRITE (90,'(A,A)') '  Add Coulomb shielding effects: ',c
         IF (fapprocoul) THEN
            c = 'yes'
         ELSE
            c = 'no'
         END IF
         WRITE (90,'(A,A)') '  Approximate f for Coulomb shielding: ',c
         WRITE (90,'(A)') 'o PARTITION FUNCTIONS'
         IF (VardyaZH2) THEN
            c = 'yes'
         ELSE
            c = 'no'
         END IF
         WRITE (90,'(A,A)') '  Use Vardya partition function for H2: ',c
         IF (IrwinZH2) THEN
            c = 'yes'
         ELSE
            c = 'no'
         END IF
         WRITE (90,'(A,A)') '  Use Irwin partition function for H2: ',c
         IF (compZH0) THEN
            c = 'yes'
         ELSE
            c = 'no'
         END IF
         WRITE (90,'(A,A)') '  Add fit for H partition function: ',c
         IF (compZHe) THEN
            c = 'yes'
         ELSE
            c = 'no'
         END IF
         WRITE (90,'(A,A)')
     &        '  Add fit for He, He+, He++ partition functions: ',c
         IF (compZCNO) THEN
            c = 'yes'
         ELSE
            c = 'no'
         END IF
         WRITE (90,'(A,A)') '  Add fits for CNO partition functions: ',c
         WRITE (90,'(A)') 'o ATOMIC MASSES'
         IF (compavmass) THEN
            c = 'isotopic averages'
         ELSE
           c = 'main isotope masses'
         END IF
         WRITE (90,'(A,A)') '  Atomic masses are: ',c
         WRITE (90,'(A)') 'o NEWTON-RAPHSON DERIVATIVES'
         IF (eosderiv) THEN
            c = 'numerically'
         ELSE
            c = 'analytically'
         END IF
         WRITE (90,'(A,A)') '  Second order derivative computation: ',c
         IF (approxsod) THEN
            c = 'yes'
         ELSE
            c = 'no'
         END IF
         WRITE (90,'(A,A)') 
     &        '  Approximates second order derivatives : ',c
         WRITE (90,'(A)') 'o DEBUGGING & MONITORING'
         IF (checkcomposition) THEN
            c = 'yes'
         ELSE
            c = 'no'
         END IF
         WRITE (90,'(A,A)') '  Monitor the composition to check: ',c
         IF (checkparacard) THEN
            c = 'yes'
         ELSE
            c = 'no'
         END IF
         WRITE (90,'(A,A)') '  Monitor the parameter card to check: ',c
         WRITE (90,'(/,53("-"),/)')
      END IF
      IF (compH2.AND.(.NOT.VardyaZH2).AND.(.NOT.IrwinZH2)) THEN
         WRITE (NOUT,*) 'ERROR: EOS - initpar'
         WRITE (NOUT,*) 'You have to choose a H2 partition function !'
         STOP
      END IF
      IF (compH2.AND.VardyaZH2.AND.IrwinZH2) THEN
         WRITE (NOUT,*) 'ERROR: EOS - initpar '
         WRITE (NOUT,*) 'Must choose only ONE H2 partition function!'
         STOP
      END IF
      END

************************************************************************
* Initialization of the ion masses                                     *
************************************************************************

      SUBROUTINE initlnm
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'evolcom.nuc'
      include 'eoscom.ion'
      include 'eoscom.cons'
      include 'eoscom.par'
*     outputs
      include 'eoscom.phys'
*     Intermediate variables
*     nbzmax-2 because H and heavy are tabulated separatly
      INTEGER nmainisoZ(nbz-2),i

      nmainisoZ(1) = ihe4
      nmainisoZ(2) = ili7
      nmainisoZ(3) = ibe9
      nmainisoZ(4) = ib11
      nmainisoZ(5) = ic12
      nmainisoZ(6) = in14
      nmainisoZ(7) = io16
      nmainisoZ(8) = if19
      nmainisoZ(9) = ine20
      nmainisoZ(10) = ina23
      nmainisoZ(11) = img24
      nmainisoZ(12) = ial27
      nmainisoZ(13) = isi28
      nmainisoZ(14) = ip31
      nmainisoZ(15) = is32
      nmainisoZ(16) = icl35
      do i=1,nbz-2
         if (nmainisoZ(i).eq.0) then
            write (nout,*) 'nmainisoZ badly initialized : check ',
     &           'definition of species index in rininet'
            stop 'eosc initlnm'
         endif
      enddo

*     ---------H-----------
      lnmH(1) = DLOG(anuc(nspmin))
      lnmH(2) = lnmH(1)
      lnmH(3) = lnmH(1)
      lnmH(4) = lnmH(1) + ln2
      IF (checkcomposition) THEN
         WRITE (nout,1000) 'm(Z=',Z(1),') =',DEXP(lnmH(1))
      END IF
*     ---------Z>1-----------
      DO i=2,nbz-1
         lnmZ(i-1) = DLOG(anuc(nmainisoZ(i-1)))
         IF (checkcomposition) THEN
            WRITE (nout,1000) 'm(Z=',Z(i),') =',DEXP(lnmZ(i-1))
         END IF
      END DO
*     ----------Heavy----------
*        lng0(56Fe)=ln(25)
c      lng0heavy = 3.2188758d0
      lnmheavy = DLOG(anuc(nspmax))
      IF (checkcomposition) THEN
         WRITE (NOUT,1000) 'm(Z=',Z(nbz),') =',DEXP(lnmheavy)
      END IF
 1000 FORMAT (A,I2,A,D14.8)
      END

************************************************************************
* Initialization of the ion masses (isotopic average)                  *
************************************************************************

      SUBROUTINE initlnmav(shell,X,Y)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.ion'
      include 'eoscom.cons'
      include 'eoscom.par'
      DOUBLE PRECISION X(nbz),Y(nbz)
      INTEGER shell
*     Intermediate variables
      INTEGER i
      DOUBLE PRECISION ratio,truncln
*     outputs
      include 'eoscom.phys'
*     ---------H-----------
      ratio = truncln(X(1)/Y(1))
      IF (ratio.ne.ratio) THEN
         lnmH(1) = zero
      ELSE
         lnmH(1) = ratio
      END IF
      lnmH(2) = lnmH(1)
      lnmH(3) = lnmH(1)
      lnmH(4) = lnmH(1) + ln2
      IF ((shell.EQ.1).AND.checkcomposition) THEN
         WRITE (NOUT,1000) '<m>(Z=',1,') =',DEXP(lnmH(1))
      END IF
*     ---------Z>1-----------
      DO i=2,nbz-1
         ratio = truncln(X(i)/Y(i))
         IF (ratio.ne.ratio) THEN
            lnmZ(i-1) = zero
         ELSE
            lnmZ(i-1) = ratio
         END IF
      IF ((shell.EQ.1).AND.checkcomposition) THEN
         WRITE (NOUT,1000) '<m>(Z=',Z(i),') =',DEXP(lnmZ(i-1))
      END IF
      END DO
*     ----------Heavy----------
      ratio = truncln(X(nbz)/Y(nbz))
      IF (ratio.ne.ratio) THEN
         lnmheavy = zero
      ELSE
         lnmheavy = ratio
      END IF
      IF ((shell.EQ.1).AND.checkcomposition) THEN
         WRITE (NOUT,1000) '<m>(Z=',Z(nbz),') =',DEXP(lnmheavy)
      END IF
 1000 FORMAT (A,I2,A,D14.8)
      END

************************************************************************
* Computation of the quantic densitie for specie of mass m at          *
*  temperature kt=1/beta (in eV)                                       *
************************************************************************

      DOUBLE PRECISION FUNCTION lnnQ(lnbetaev,lnm)
*     --------------------
      IMPLICIT NONE
*     inputs
      DOUBLE PRECISION lnbetaev,lnm
      include 'evolpar.star'
      include 'eoscom.cons'
*     --------------------
      lnnQ = 1.5d0 * (lncoeffnQ+lnm-lnbetaev)
      END

************************************************************************
* Compute xlnx and derivatives                                         *
*        compute xlnx(1)=x*lnx=x(m,n1)*ln(m,n1)                        *
*                xlnx(2)=d(x*lnx)/dlnT and xlnx(3)=d(x*lnx)/dlnf       *
* Inuts : x(mmax,nmax), index m,n1,n2,n3                              *
* Output : xlnx(deriv2)                                                *
************************************************************************

      SUBROUTINE loga_sod(x,mmax,nmax,m,n1,n2,n3,n4,n5,n6,xlnx)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'eoscom.math'
      INTEGER nmax,mmax,m,n1,n2,n3,n4,n5,n6
      DOUBLE PRECISION x(mmax,nmax)
*     outputs
      DOUBLE PRECISION xlnx(deriv2)
*     Intermediate variable
      DOUBLE PRECISION lnx,xinv,lnx1,truncln,var,truncinv
*     --------------------
      var = x(m,n1)
      lnx = truncln(var)
      lnx1 = lnx + 1.0d0
      xinv = truncinv(var)
      xlnx(1) = var * lnx
      xlnx(2) = x(m,n2) * lnx1
      xlnx(3) = x(m,n3) * lnx1
      xlnx(4) = x(m,n4) * lnx1 + x(m,n2) * x(m,n2) * xinv
      xlnx(6) = x(m,n6) * lnx1 + x(m,n2) * x(m,n3) * xinv
      xlnx(5) = x(m,n5) * lnx1 + x(m,n3) * x(m,n3) * xinv
      END

************************************************************************
* Computation of function g related to free energy by F=-NkTg(x,y)     *
*     and derivatives. g will be used to get the thermodynamics        *
*     corrections due to non-perfect effects.                          *
* Inputs: x=ne*mh, y=13.6/kT, zavY=<z>y  zavYZ=<yz>y f and rootfinv    *
*         (the last for Coulomb corrections)                           *
* Outputs: fitxy(deriv3) contains g and its derivatives                *
*          quantities dmu... and ds... will be used later              *
************************************************************************

      SUBROUTINE nonperfect_sod(f,rootfinv,x,y,zavY,zavYZ,fitxy,dmut1,
     &dmut2,dmut1x,dmut1y,dmut2x,dst1,dst2,dst1x,dst2y,dst2x,halff1,
     &ine,fadd1inv,fitnpi)
*     ----------declarations----------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.math'
      include 'eoscom.ion'
      include 'eoscom.par'
      DOUBLE PRECISION x,y,zavy,zavyz,rootfinv,f,halff1,fadd1inv
      DOUBLE PRECISION ine(deriv2)
*     ouputs
      DOUBLE PRECISION fitxy(deriv3),fitnpi(deriv3)
      DOUBLE PRECISION dmut1,dmut2,dmut1x,dmut1y,dmut2x
      DOUBLE PRECISION dst1,dst2,dst1x,dst2y,dst2x
*     Intermediate variables
      DOUBLE PRECISION coco,cocox,cocoy,cocoxx,cocoyy,cocoxy
      DOUBLE PRECISION cocoxxx,cocoyyy,cocoyyx,cocoxxy
      DOUBLE PRECISION copi,copix,copiy,copixy,copixx,copiyy
      DOUBLE PRECISION copixxx,copiyyy,copiyyx,copixxy
      DOUBLE PRECISION rootfap,fapx,term,sum,lnfap
      DOUBLE PRECISION etap,etapx,etapxx,etapxxx,sumx
*     --------------------
      IF (fapprocoul.OR.compionpres) THEN
         CALL fapprox_sod(x,y,rootfap,fapx,term,sum,sumx,lnfap)
         CALL etapprox_sod(rootfap,fapx,term,sum,lnfap,etap,
     &etapx,etapxx,etapxxx)
      END IF
      IF (compcoulshiel)THEN
         CALL coulomb_sod(f,x,y,zavY,zavYZ,coco,cocox,cocoy,cocoxx,
     &cocoyy,cocoxy,cocoxxx,cocoyyy,cocoyyx,cocoxxy,etapx,term,sum,
     &sumx,halff1,ine,rootfinv,fadd1inv)
      ELSE
         coco = 0.0d0
         cocox = 0.0d0
         cocoy = 0.0d0
         cocoxx = 0.0d0
         cocoyy = 0.0d0
         cocoxy = 0.0d0
         cocoxxx = 0.0d0
         cocoyyy = 0.0d0
         cocoxxy = 0.0d0
         cocoyyx = 0.0d0
      END IF
      IF (compionpres) THEN
         CALL pressioni_sod(x,y,etap,etapx,etapxx,etapxxx,copi,copix,
     &copiy,copixy,copixx,copiyy,copixxx,copiyyy,copiyyx,copixxy)
      ELSE
         copi = 0.0d0
         copix = 0.0d0
         copiy = 0.0d0
         copixx = 0.0d0
         copiyy = 0.0d0
         copixy = 0.0d0
         copixxx = 0.0d0
         copiyyy = 0.0d0
         copixxy = 0.0d0
         copiyyx = 0.0d0
      END IF
      fitxy(1) = copi + coco
      fitxy(2) = copix + cocox
      fitxy(3) = copiy + cocoy
      fitxy(4) = copixx + cocoxx
      fitxy(5) = copiyy + cocoyy
      fitxy(6) = copixy + cocoxy
      fitxy(7) = copixxx + cocoxxx
      fitxy(8) = copiyyy + cocoyyy
      fitxy(9) = copixxy + cocoxxy
      fitxy(10) = copiyyx + cocoyyx
      fitnpi(1) = coco
      fitnpi(2) = cocox
      fitnpi(3) = cocoy
      fitnpi(4) = cocoxx
      fitnpi(5) = cocoyy
      fitnpi(6) = cocoxy
      fitnpi(7) = cocoxxx
      fitnpi(8) = cocoyyy
      fitnpi(9) = cocoxxy
      fitnpi(10) = cocoyyx
*     --------------------
      dmut1 = fitxy(3) + fitxy(6)
      dmut2 = fitxy(2) + fitxy(4)
      dmut1x = fitxy(6) + fitxy(9)
      dmut1y = fitxy(5) + fitxy(10)
      dmut2x = fitxy(4) + fitxy(7)
*     --------------------
      dst1 = fitxy(2) - fitxy(6)
      dst2 = fitxy(5) - fitxy(3)
      dst1x = fitxy(4) - fitxy(9)
      dst2x = fitxy(10) - fitxy(6)
      dst2y = fitxy(8) - fitxy(5)
      END

************************************************************************
* Computation of function g related to free energy by F=-NkTg(x,y)     *
*     and derivatives in case of totaly ionized matter. g will be used *
*     to get the thermodynamics corrections due to non-perfect effects.*
*     Only Ionization pressure is calculated in this case because we   *
*     have to correct the corrections such as DFpitot-->0 when Ne-->Ne0*
* Inputs: x=ne0*mh y=13.6/kT                                           *
* Outputs: fitxy0(deriv3) contains g and its derivatives               *
*          quantities dmu... and ds... will be used later              *
************************************************************************

      SUBROUTINE nonperfect0_sod(x,y,fitxy0,dmut1,dmut2,dmut1x,dmut1y,
     &dmut2x,dst1,dst2,dst1x,dst2y,dst2x)
*     ----------declarations----------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.math'
      include 'eoscom.ion'
      include 'eoscom.par'
      DOUBLE PRECISION x,y
*     ouputs
      DOUBLE PRECISION fitxy0(deriv3)
      DOUBLE PRECISION copi,copix,copiy,copixy,copixx,copiyy
      DOUBLE PRECISION copixxx,copiyyy,copiyyx,copixxy
      DOUBLE PRECISION dmut1,dmut2,dmut1x,dmut1y,dmut2x
      DOUBLE PRECISION dst1,dst2,dst1x,dst2y,dst2x
*     Intermediate variables
      DOUBLE PRECISION rootfap,fapx,term,sum,lnfap
      DOUBLE PRECISION etap,etapx,etapxx,etapxxx,sumx
*     --------------------
      CALL fapprox_sod(x,y,rootfap,fapx,term,sum,sumx,lnfap)
      CALL etapprox_sod(rootfap,fapx,term,sum,lnfap,etap,
     &etapx,etapxx,etapxxx)
      CALL pressioni_sod(x,y,etap,etapx,etapxx,etapxxx,copi,copix,
     &copiy,copixy,copixx,copiyy,copixxx,copiyyy,copiyyx,copixxy)
      fitxy0(1) = copi
      fitxy0(2) = copix
      fitxy0(3) = copiy
      fitxy0(4) = copixx
      fitxy0(5) = copiyy
      fitxy0(6) = copixy
      fitxy0(7) = copixxx
      fitxy0(8) = copiyyy
      fitxy0(9) = copixxy
      fitxy0(10) = copiyyx
*     --------------------
      dmut1 = fitxy0(3) + fitxy0(6)
      dmut2 = fitxy0(2) + fitxy0(4)
      dmut1x = fitxy0(6) + fitxy0(9)
      dmut1y = fitxy0(5) + fitxy0(10)
      dmut2x = fitxy0(4) + fitxy0(7)
*     --------------------
      dst1 = fitxy0(2) - fitxy0(6)
      dst2 = fitxy0(5) - fitxy0(3)
      dst1x = fitxy0(4) - fitxy0(9)
      dst2x = fitxy0(10) - fitxy0(6)
      dst2y = fitxy0(8) - fitxy0(5)
      END

************************************************************************
* Pressure ionization corrections computation                          *
*  These effects are fitted by means of function gpi(x,y) (here copi)  *
*   DF = - Ne * k * T * g(x,y), with x=ne*mh and y=13.6/kT             *
*  Here, gpi = fac * sum.                                              *
* Inputs: x=ne*mh,y=13.6/kT,etap and derivatives from routine etapprox *
* Outputs: copi and derivatives versus x and y                         *
************************************************************************

      SUBROUTINE pressioni_sod(x,y,etap,etapx,etapxx,etapxxx,copi,copix,
     &copiy,copixy,copixx,copiyy,copixxx,copiyyy,copiyyx,copixxy)
*     ----------declarations----------
      IMPLICIT NONE
*     inputs
      DOUBLE PRECISION x,y
      DOUBLE PRECISION etapx,etapxx,etap,etapxxx
*     ouputs
      DOUBLE PRECISION copi,copix,copiy,copixy,copixx,copiyy
      DOUBLE PRECISION copixxx,copiyyy,copiyyx,copixxy
*     intermediate variables
      DOUBLE PRECISION produ,xc,lnarg,lnargx
      DOUBLE PRECISION exparg,cexparg,expargx,fac,facx,facxx,facxxx
      DOUBLE PRECISION sumy,sumxx,sumxy,sumxxx,sumyyy,sumxxy
      DOUBLE PRECISION sumyyx,sumyy,sum,sumx,truncexpo,truncln,truncinv
*     constantes
      DOUBLE PRECISION c1,c2,c3,c4inv
*     ----------constantes----------
*     pressure ionization fit coefficients
      PARAMETER (c1=3.0d0)
      PARAMETER (c2=0.25d0)
      PARAMETER (c3=2.0d0)
      PARAMETER (c4inv=1.d2/3.d0)
*     ----------sum and derivatives calculus----------
      xc = x * c4inv
      lnarg = 1.0d0 / (1.0d0+xc)
      lnargx = -xc * lnarg
      produ = -c3 * lnarg * lnargx
      sum = y + etap - c3 * truncln(lnarg)
      sumx = etapx - c3 * lnargx
      sumy = y + 1.5d0 * etapx
      sumxx = etapxx + produ
      sumxy = 1.5d0 * etapxx
      sumyy = y + 1.5d0 * sumxy
      sumxxx = etapxxx + produ * (lnarg+lnargx)
      sumxxy = 1.5d0 * etapxxx
      sumyyx = 1.5d0 * sumxxy
      sumyyy = y + 1.5d0 * sumyyx
*     ----------fac and derivatives calculus----------
      exparg = (c1*truncinv(x))**c2
      cexparg = c2 * exparg
      expargx = -cexparg
      fac = truncexpo(-exparg)
      facx = cexparg * fac
      cexparg = cexparg - c2
      facxx = cexparg * facx
      facxxx = cexparg * facxx + c2 * facx * expargx
*     ----------calculus of g and derivatives----------
      copi = fac * sum
      copix = fac * sumx + facx * sum
      copiy = fac * sumy
      copixx = fac * sumxx + 2.0d0 * facx * sumx + sum * facxx
      copiyy = fac * sumyy
      copixy = facx * sumy + fac * sumxy
      copixxx = fac * sumxxx + 3.0d0 * (facx*sumxx+facxx*sumx) +
     &facxxx * sum
      copiyyy = fac * sumyyy
      copixxy = facxx * sumy + 2.0d0 * facx * sumxy + fac * sumxxy
      copiyyx = facx * sumyy + fac * sumyyx
      END

************************************************************************
* Computation of the thermodynamics quantities of the radiativ         *
*  component                                                           *
************************************************************************

      SUBROUTINE radiativ_sod(T,pradi,sradi,uradi,rho,rhoinv)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'evolcom.cons'
      include 'eoscom.math'
      include 'eoscom.cons'
      DOUBLE PRECISION rho(deriv2),rhoinv
      DOUBLE PRECISION T
*     Intermediate variables
      DOUBLE PRECISION T3
*     outputs
      DOUBLE PRECISION pradi(deriv2),sradi(deriv2),uradi(deriv2)
*     --------------------
      T3 = 4.d0*sig/c * T**3 * inv3
      pradi(1) = T3 * T
      pradi(2) = 4.0d0 * pradi(1)
      pradi(3) = 0.0d0
      pradi(4) = 16.0d0 * pradi(1)
      pradi(5) = 0.0d0
      pradi(6) = 0.0d0
      sradi(1) = 4.0d0 * T3 * rhoinv
      sradi(2) = (3.0d0-rho(2)) * sradi(1)
      sradi(3) = -rho(3) * sradi(1)
      sradi(4) = -rho(4) * sradi(1) + (3.0d0-rho(2)) * sradi(2)
      sradi(5) = -rho(5) * sradi(1) - rho(3) * sradi(3)
      sradi(6) = -rho(6) * sradi(1) - rho(3) * sradi(2)
      uradi(1) = 3.0d0 * pradi(1) * rhoinv
      uradi(2) = uradi(1) * (4.0d0-rho(2))
      uradi(3) = -uradi(1) * rho(3)
      uradi(4) = 0.d0
      uradi(5) = 0.d0
      uradi(6) = 0.d0
      END

************************************************************************
* Solve SAHA equations for H,H+,H-,H2 and e- (coupled because of H2)   *
* INPUTS: Y,nsh,nbz,betaev,heavelec,lnparfunH                           *
* OUTPUTS:abundH(4,deriv2),adundelec(deriv2)                           *
************************************************************************

      SUBROUTINE sahaH_sod(shell,Y,betaev,lnbetaev,heavelec,
     &abundelecinv,abundH,abundelec,lnparfunH,etaeff,ine,solutionH,
     &zmean,T,etanpi,sumh,diffh,ah,bh,sdelh,error)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.math'
      include 'eoscom.ion'
      include 'eoscom.cons'
      include 'eoscom.phys'
      include 'eoscom.par'
      DOUBLE PRECISION heavelec(deriv2),betaev,Y(nbz),lnbetaev
      DOUBLE PRECISION etaeff(deriv2),ine(deriv2),T,etanpi(deriv2)
      DOUBLE PRECISION lnparfunH(nbionHmax,deriv3)
*     Intermediate variables
      DOUBLE PRECISION lnHmH,lnHpH,lnH2H,ratHpH,ratHmH,ratH2H
      DOUBLE PRECISION ratHpHt,ratHmHt,ratH2Ht,ratHpHf,ratHmHf,ratH2Hf
      DOUBLE PRECISION ratH2Htt,ratH2Hff,ratH2HfT,ratHmHtt,ratHmHff
      DOUBLE PRECISION ratHmHfT,ratHpHTT,ratHpHff,ratHpHfT
      DOUBLE PRECISION H1,H2,H,limem,limH2,lim,lnne,truncln
      DOUBLE PRECISION truncexpo,lnnQ,dlnf,dlnT,dlnff,dlnTT,dlnfT
      DOUBLE PRECISION khiHHm,khiHHp,dissocH2,denom,denominv
      DOUBLE PRECISION ta1,ta2,ta3,ta4,ta5,tb1,tb2,tb3,tb4,tb5
      DOUBLE PRECISION tc1,tc2,tc3,tc4,tc5,sum,diff,ch
      DOUBLE PRECISION dev(4),reel(4),complexe(4,2)
      INTEGER i,j,shell,nbr,nbc
      LOGICAL egal,compH2eff
*     ouputs
      DOUBLE PRECISION abundH(nbionHmax,deriv2),abundelec(deriv2)
      DOUBLE PRECISION abundelecinv,zmean(nsp),sumh,diffh,ah,bh,sdelh
      LOGICAL solutionH
      INTEGER error
*     ---------affectations----------
      DO i=1,nbionHmax
         DO j=1,deriv2
            abundH(i,j) = 0.0d0
         END DO
      END DO
      lnne = truncln(ine(1)*pi8l3)
      dissocH2 = dissH2 * betaev
      khiHHp = khiH * betaev
      khiHHm = (enerlevelH(3)-enerlevelH(1)) * betaev
      compH2eff = compH2 .AND. (T.LE.TmaxH2)
*     ----------ratios and derivatives d/dln, dd/dln2----------
*     ---------- H+ ----------
      lnHpH = lnparfunH(2,1) - lnparfunH(1,1) - etaeff(1) - khiHHp
      dlnT = lnparfunH(2,2) - lnparfunH(1,2) - etaeff(2) + khiHHp
      dlnTT = lnparfunH(2,4) - lnparfunH(1,4) - etaeff(4) - khiHHp
      dlnf = lnparfunH(2,3) - lnparfunH(1,3) - etaeff(3)
      dlnff = lnparfunH(2,5) - lnparfunH(1,5) - etaeff(5)
      dlnfT = lnparfunH(2,6) - lnparfunH(1,6) - etaeff(6)
      ratHpH = truncexpo(lnHpH)
      ratHpHt = dlnT * ratHpH
      ratHpHf = dlnf * ratHpH
      ratHpHtt = ratHpH * (dlnTT+dlnT*dlnT)
      ratHpHff = ratHpH * (dlnff+dlnf*dlnf)
      ratHpHfT = ratHpH * (dlnfT+dlnf*dlnT)
*     ---------- H- ----------
      IF (compHm) THEN
          lnHmH = lnparfunH(3,1) - lnparfunH(1,1) + etanpi(1) - khiHHm
         dlnT = lnparfunH(3,2) - lnparfunH(1,2) + etanpi(2) + khiHHm
         dlnTT = lnparfunH(3,4) - lnparfunH(1,4) + etanpi(4) - khiHHm
         dlnf = lnparfunH(3,3) - lnparfunH(1,3) + etanpi(3)
         dlnff = lnparfunH(3,5) - lnparfunH(1,5) + etanpi(5)
         dlnfT = lnparfunH(3,6) - lnparfunH(1,6) + etanpi(6)
         ratHmH = truncexpo(lnHmH)
         ratHmHt = dlnT * ratHmH
         ratHmHf = dlnf * ratHmH
         ratHmHtt = ratHmH * (dlnTT+dlnT*dlnT)
         ratHmHff = ratHmH * (dlnff+dlnf*dlnf)
         ratHmHfT = ratHmH * (dlnfT+dlnf*dlnT)
      ELSE
         ratHmH = 0.0d0
         ratHmHt = 0.0d0
         ratHmHtt = 0.0d0
         ratHmHf = 0.0d0
         ratHmHff = 0.0d0
         ratHmHfT = 0.0d0
      END IF
*     ---------- H2 ----------
      IF (compH2eff) THEN
         lnH2H = lnparfunH(4,1) - lnnQ(lnbetaev,lnmh(4)+lnamu-ln4) +
     &lnne - 2.0d0 * lnparfunH(1,1) + dissocH2
         dlnT = Ine(2) + lnparfunH(4,2) - 2.0d0 * lnparfunH(1,2) -
     &1.5d0 - dissocH2
         dlnTT = Ine(4) + lnparfunH(4,4) - 2.0d0 * lnparfunH(1,4) +
     &dissocH2
         dlnf = Ine(3) + lnparfunH(4,3) - 2.0d0 * lnparfunH(1,3)
         dlnff = Ine(5) + lnparfunH(4,5) - 2.0d0 * lnparfunH(1,5)
         dlnfT = Ine(6) + lnparfunH(4,6) - 2.0d0 * lnparfunH(1,6)
         ratH2H = truncexpo(lnH2H)
         ratH2Ht = dlnT * ratH2H
         ratH2Hf = dlnf * ratH2H
         ratH2Htt = ratH2H * (dlnTT+dlnT*dlnT)
         ratH2Hff = ratH2H * (dlnff+dlnf*dlnf)
         ratH2HfT = ratH2H * (dlnfT+dlnf*dlnT)
      ELSE
         ratH2H = 0.0d0
         ratH2Ht = 0.0d0
         ratH2Htt = 0.0d0
         ratH2Hf = 0.0d0
         ratH2Hff = 0.0d0
         ratH2HfT = 0.0d0
      END IF
*     ---------Constraints to keep Ne- and Nh2 > 0-----------
      IF (egal(ratHpH,ratHmH)) THEN
         ratHpH = ratHmH + zero
      END IF
      diffh = ratHpH - ratHmH
      sum = heavelec(1) + (heavelec(1)+Y(1)) * ratHpH +
     &     (heavelec(1)-Y(1)) * ratHmH
c      print * ,'sum',sum, 'heavelec(1)',heavelec(1),'Y(1)',Y(1)
c      print *,'ratHpH',ratHpH,'ratHmH',ratHmH
c      stop
      sumh = 1.0d0 + ratHpH + ratHmH
      limH2 = Y(1) / sumh
      limem = -heavelec(1) / diffh
      IF (egal(sum,Y(1)*diffh)) THEN
         sum = sumh * diffh + 2.0d0 * ratH2H
         IF (sum.LE.0.0d0) THEN
            WRITE (NOUT,*)
            WRITE (NOUT,*) 'ERROR - EOS - sahah_sod - non ',
     &           'physical H+/H and H-/H values !!!'
            WRITE (NOUT,'(" shell=",i5,", H+/H=",1pe11.3,", H-/H=",
     &           1pe11.3)') shell,ratHpH,ratHmH
            WRITE (NOUT,*) 'heavelec(1)=',heavelec(1),' Y(1)=',Y(1)
            WRITE (NOUT,*) 'H2*Ne/(H*H)=',ratH2H
            error = 11
            RETURN
         END IF
      ELSE
         IF (sum.LE.0.0d0) THEN
            WRITE (NOUT,*)
            WRITE (NOUT,*) 'ERROR - EOS - sahah_sod - non ',
     &           ' physical H+/H and H-/H values !!!'
            WRITE (NOUT,'(" shell=",i5,", H+/H=",1pe11.3,", H-/H=",
     &           1pe11.3,", T=",1pe11.3)') shell,ratHpH,ratHmH,T
            WRITE (NOUT,*) 'heavelec(1)=',heavelec(1),' Y(1)=',Y(1)
            error = 11
            RETURN
         END IF
      END IF
*     ---------- Solve a second order polynomial (!) ----------
      ah = sumh * diffh + 2.0d0 * ratH2H
      bh = heavelec(1) * sumh - Y(1) * diffh
      ch = -heavelec(1) * Y(1)
      CALL secondegree(ah,bh,ch,nbr,nbc,reel,complexe,dev,.true.,sdelh)
      IF (nbr.EQ.0) THEN
         WRITE (NOUT,*) 'ERROR - sahaH_sod - no solutions in R'
         WRITE (NOUT,*) 'Shell=',shell,' real solutions=',nbr
         WRITE (NOUT,*) 'a=',ah,' b=',bh,' c=',ch
         error = 11
         return
      ELSE
         H1 = reel(1)
         H2 = reel(2)
      END IF
*     ---------- Choice of the 'good' root ----------
      IF (.NOT.compH2eff) THEN
         H1 = limH2
         H2 = limem
      END IF
      CALL choice(shell,H1,H2,limH2,limem,diffh,Y(1),H,lim,solutionH,
     &ah,bh,ch,ratH2H,dev,error)
      if (error.gt.0) return
*     ----------Abundances----------
      abundH(1,1) = H
      abundH(2,1) = ratHpH * abundH(1,1)
      abundH(3,1) = ratHmH * abundH(1,1)
      IF (egal(abundH(3,1),heavelec(1)+abundH(2,1))) THEN
         IF (totioniz) THEN
            abundelec(1) = numericalimit
         ELSE
            abundelec(1) = zero
         END IF
      ELSE
         abundelec(1) = heavelec(1) + abundH(2,1) - abundH(3,1)
      END IF
      IF (compH2eff) THEN
         IF (egal(Y(1),abundH(1,1)*sumh)) THEN
            abundH(4,1) = numericalimit
         ELSE
            abundH(4,1) = 0.5d0 * (Y(1)-abundH(1,1)-abundH(2,1)-
     &abundH(3,1))
         END IF
      ELSE
         abundH(4,1) = 0.0d0
      END IF
*     ---------- First order derivatives ----------
      sum = abundelec(1) + 2.0d0 * abundH(4,1)
      diff = abundelec(1) - 2.0d0 * abundH(4,1)
      denom = 4.0d0 * ratH2H * abundH(1,1) + abundelec(1) +
     &sum * ratHmH + diff * ratHpH
      denominv = 1.0d0 / denom
      abundH(1,2) = denominv * (2.0d0*abundH(4,1)*heavelec(2)-
     &abundH(1,1)*(sum*ratHmHt+diff*ratHpHt+2.0d0*abundH(1,1)*ratH2Ht))
      abundH(1,3) = denominv * (2.0d0*abundH(4,1)*heavelec(3)-
     &abundH(1,1)*(sum*ratHmHf+diff*ratHpHf+2.0d0*abundH(1,1)*ratH2Hf))
      abundH(2,2) = abundH(1,1) * ratHpHt + ratHpH * abundH(1,2)
      abundH(2,3) = abundH(1,1) * ratHpHf + ratHpH * abundH(1,3)
      abundH(3,2) = abundH(1,1) * ratHmHT + ratHmH * abundH(1,2)
      abundH(3,3) = abundH(1,1) * ratHmHf + ratHmH * abundH(1,3)
      IF ((abundelec(1).EQ.numericalimit).OR.
     &(abundelec(1).EQ.zero)) THEN
         abundelec(2) = zero
         abundelec(3) = zero
      ELSE
         abundelec(2) = heavelec(2) + abundH(2,2) - abundH(3,2)
         abundelec(3) = heavelec(3) + abundH(2,3) - abundH(3,3)
      END IF
      IF (compH2eff) THEN
         IF (abundH(4,1).EQ.numericalimit) THEN
            abundH(4,2) = zero
            abundH(4,3) = zero
         ELSE
            abundH(4,2) = -0.5d0 * (abundH(1,2)+abundH(2,2)+
     &abundH(3,2))
            abundH(4,3) = -0.5d0 * (abundH(1,3)+abundH(2,3)+
     &abundH(3,3))
         END IF
      ELSE
         abundH(4,2) = 0.0d0
         abundH(4,3) = 0.0d0
      END IF
*     ---------- second order derivatives ----------
      ta1 = 2.0d0 * abundelec(2) * abundH(4,2) -
     &4.0d0 * abundH(1,1) * abundH(1,2) * ratH2HT -
     &abundH(1,1) * abundH(1,1) * ratH2HTT -
     &2.0d0 * ratH2H * abundH(1,2) * abundH(1,2)
      ta2 = 2.0d0 * ratHpHT * abundH(1,2) + abundH(1,1) * ratHpHTT
      ta3 = 2.0d0 * ratHmHT * abundH(1,2) + abundH(1,1) * ratHmHTT
      ta4 = 0.0d0
      ta5 = heavelec(4)
      abundH(1,4) = (2.0d0*(ta1+abundH(4,1)*ta5)+abundelec(1)*ta4-
     &sum*ta3-diff*ta2) * denominv
      tb1 = 2.0d0 * abundelec(3) * abundH(4,3) -
     &4.0d0 * abundH(1,1) * abundH(1,3) * ratH2Hf -
     &abundH(1,1) * abundH(1,1) * ratH2Hff -
     &2.0d0 * ratH2H * abundH(1,3) * abundH(1,3)
      tb2 = 2.0d0 * ratHpHf * abundH(1,3) + abundH(1,1) * ratHpHff
      tb3 = 2.0d0 * ratHmHf * abundH(1,3) + abundH(1,1) * ratHmHff
      tb4 = 0.0d0
      tb5 = heavelec(5)
      abundH(1,5) = (2.0d0*(tb1+abundH(4,1)*tb5)+abundelec(1)*tb4-
     &sum*tb3-diff*tb2) * denominv
      tc1 = abundelec(2) * abundH(4,3) + abundelec(3) * abundH(4,2) -
     &2.0d0 * abundH(1,1) * (abundH(1,3)*ratH2HT+abundH(1,2)*ratH2Hf) -
     &abundH(1,1) * abundH(1,1) * ratH2HfT -
     &2.0d0 * ratH2H * abundH(1,3) * abundH(1,2)
      tc2 = ratHpHf * abundH(1,2) + ratHpHT * abundH(1,3) +
     &abundH(1,1) * ratHpHfT
      tc3 = ratHmHf * abundH(1,2) + ratHmHT * abundH(1,3) +
     &abundH(1,1) * ratHmHfT
      tc4 = 0.0d0
      tc5 = heavelec(6)
      abundH(1,6) = (2.0d0*(tc1+abundH(4,1)*tc5)+abundelec(1)*tc4-
     &sum*tc3-diff*tc2) * denominv
      abundH(2,4) = abundH(1,1) * ratHpHTT +
     &2.0d0 * abundH(1,2) * ratHpHT + ratHpH * abundH(1,4)
      abundH(2,5) = abundH(1,1) * ratHpHff +
     &2.0d0 * abundH(1,3) * ratHpHf + ratHpH * abundH(1,5)
      abundH(2,6) = abundH(1,1) * ratHpHfT + abundH(1,3) * ratHpHT +
     &abundH(1,2) * ratHpHf + ratHpH * abundH(1,6)
      abundH(3,4) = abundH(1,1) * ratHmHTT +
     &2.0d0 * abundH(1,2) * ratHmHT + ratHmH * abundH(1,4)
      abundH(3,5) = abundH(1,1) * ratHmHff +
     &2.0d0 * abundH(1,3) * ratHmHf + ratHmH * abundH(1,5)
      abundH(3,6) = abundH(1,1) * ratHmHfT + ratHmHf * abundH(1,2) +
     &abundH(1,3) * ratHmHT + ratHmH * abundH(1,6)
      IF ((abundelec(1).EQ.numericalimit).OR.
     &(abundelec(1).EQ.zero)) THEN
         abundelec(4) = zero
         abundelec(5) = zero
         abundelec(6) = zero
      ELSE
         abundelec(4) = heavelec(4) + abundH(2,4) - abundH(3,4)
         abundelec(5) = heavelec(5) + abundH(2,5) - abundH(3,5)
         abundelec(6) = heavelec(6) + abundH(2,6) - abundH(3,6)
      END IF
      IF (compH2eff) THEN
         IF (abundH(4,1).EQ.numericalimit) THEN
            abundH(4,4) = zero
            abundH(4,5) = zero
            abundH(4,6) = zero
         ELSE
            abundH(4,4) = -0.5d0 * (abundH(1,4)+abundH(2,4)+
     &abundH(3,4))
            abundH(4,5) = -0.5d0 * (abundH(1,5)+abundH(2,5)+
     &abundH(3,5))
            abundH(4,6) = -0.5d0 * (abundH(1,6)+abundH(2,6)+
     &abundH(3,6))
         END IF
      ELSE
         abundH(4,4) = 0.0d0
         abundH(4,5) = 0.0d0
         abundH(4,6) = 0.0d0
      END IF
*     abundelecinv must not be truncated because it is used in rho
*      computation
      abundelecinv = 1.0d0 / abundelec(1)
*     ----------zmean----------
      zmean(2) = (abundH(2,1)-abundH(3,1)) / Y(1)
      zmean(3) = zmean(2)
*     ---------- Check the positivity ----------
      DO i=2,nbionHmax
         IF (abundH(i,1).LT.0.0d0) THEN
            WRITE (NOUT,*) 'ERROR - EOS - sahah_sod - abundance<0 !!!!'
            WRITE (NOUT,*) 'shell=',shell,' abundH(',i,',1)=',
     &           abundH(i,1)
            error = 11
            return
         END IF
      END DO
      IF (abundelec(1).LT.0.0d0) THEN
         WRITE (NOUT,*) 'ERROR - EOS - sahah_sod - abundelec<0 !!!!'
         WRITE (NOUT,*) 'shell=',shell,' abundelec(1)=',abundelec(1)
         error = 11
         return
      END IF
      END

************************************************************************
* Solve SAHA equations for Z>1                                         *
************************************************************************

      SUBROUTINE sahaZ_sod(shell,betaev,Y,abundZ,lnparfunZ,etaeff,
     &error)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.ion'
      include 'eoscom.math'
      include 'eoscom.phys'
      include 'eoscom.par'
      include 'eoscom.cons'
*
      DOUBLE PRECISION betaev,Y(nbz),lnparfunZ(nbionZmax,deriv3)
      DOUBLE PRECISION etaeff(deriv2)
*     Intermediate variables
      INTEGER nbprodmax
      PARAMETER (nbprodmax=nbz*(nbz-1)/2)
      DOUBLE PRECISION dE,truncexpo,lnrat(nbprodmax),prodinv(nbprodmax)
      DOUBLE PRECISION eta,etaf,etat,etaff,etatt,etaft
      DOUBLE PRECISION ioni(nbz+1),ioniT(nbz+1),ionif(nbz+1)
      DOUBLE PRECISION ioniff(nbz+1),ioniTT(nbz+1),ionifT(nbz+1)
      DOUBLE PRECISION prod(nbprodmax),prodT(nbprodmax),prodf(nbprodmax)
      DOUBLE PRECISION prodff(nbprodmax),prodTT(nbprodmax)
      DOUBLE PRECISION prodfT(nbprodmax),testinf,sgn,aphi2
      DOUBLE PRECISION phi(nbz+1),phiT(nbz+1),phif(nbz+1)
      DOUBLE PRECISION phiff(nbz+1),phiTT(nbz+1),phifT(nbz+1),phinv
      INTEGER i,j,m,rangE,debutE,rangA,debutA,q,indx,shell,rangT
      INTEGER rangI,debutI,mul,first,ZZ,rang
*     ouputs
      DOUBLE PRECISION abundZ(nbionZmax,deriv2)
      INTEGER error
*     --------------------
      rang(mul,first,ZZ) = INT((mul-1)*ZZ-0.5d0*(mul-1)*(mul-2)+
     &     first)
*     --------------------
      DO i=1,nbionZmax
         DO j=1,deriv2
            abundZ(i,j) = 0.0d0
         END DO
      END DO
      debutA = 1
      DO i=2,nbz-1
         IF (ionized(i)) THEN
            eta = etaeff(1)
            etaf = etaeff(3)
            etaT = etaeff(2)
            etaff = etaeff(5)
            etafT = etaeff(6)
            etaTT = etaeff(4)
            debutE = indx(2,Z(i)-1) + 1
*
*           computation of ionization ratios and derivatives
*
            DO q=1,Z(i)
               rangE = debutE + q
               rangA = debutA + q
               dE = enerlevelZ(rangE) - enerlevelZ(rangE-1)
               lnrat(q) = lnparfunZ(rangA,1) - lnparfunZ(rangA-1,1) -
     &              eta - dE * betaev
               ioni(q) = truncexpo(lnrat(q))
               IF (lnrat(q).GT.limexp) THEN
                  WRITE (NOUT,*)
                  WRITE (NOUT,*) 'EOS - WARNING - sahaZ_sod - shell=',
     &                 shell,',Z=',Z(i)
                  error = 12
                  RETURN
               END IF
               ioniT(q) = lnparfunZ(rangA,2) - lnparfunZ(rangA-1,2)
     &              - etaT + dE * betaev
               ionif(q) = lnparfunZ(rangA,3) - lnparfunZ(rangA-1,3)
     &              - etaf
               ioniTT(q) = lnparfunZ(rangA,4) - lnparfunZ(rangA-1,4)
     &              - etaTT - dE * betaev
               ioniff(q) = lnparfunZ(rangA,5) - lnparfunZ(rangA-1,5)
     &              - etaff
               ionifT(q) = lnparfunZ(rangA,6) - lnparfunZ(rangA-1,6)
     &              - etafT
            END DO
*           Initializations
            DO j=1,nbprodmax
               prod(j) = 1.0d0
               prodf(j) = 0.0d0
               prodT(j) = 0.0d0
               prodff(j) = 0.0d0
               prodTT(j) = 0.0d0
               prodfT(j) = 0.0d0
            END DO
            DO j=1,Z(i)+1
               phi(j) = 1.0d0
               phif(j) = 0.0d0
               phiT(j) = 0.0d0
               phiff(j) = 0.0d0
               phiTT(j) = 0.0d0
               phifT(j) = 0.0d0
            END DO
            rangT = 1
*
*           computation of products of ionization ratios used for abundances
*           and their derivatives
*
            DO m=1,Z(i)
               DO debutI=1,Z(i)-m+1
                  rangI = debutI
                  testinf = 0.0d0
                  DO WHILE ((rangI-debutI).LT.m)
                     prod(rangT) = prod(rangT) * ioni(rangI)
                     testinf = testinf + lnrat(rangI) * ln10inv
                     prodf(rangT) = prodf(rangT) + ionif(rangI)
                     prodT(rangT) = prodT(rangT) + ioniT(rangI)
                     prodff(rangT) = prodff(rangT) + ioniff(rangI)
                     prodTT(rangT) = prodTT(rangT) + ioniTT(rangI)
                     prodfT(rangT) = prodfT(rangT) + ionifT(rangI)
                     rangI = rangI + 1
                  END DO
*                 truncation at liminf to avoid INF or -INF
                  IF (ABS(testinf).GT.liminf) THEN
                     sgn = testinf / ABS(testinf)
                     prod(rangT) = 10**(liminf*sgn)
                  END IF
                  prodinv(rangT) = 1.0d0 / prod(rangT)
                  rangT = rangT + 1
               END DO
            END DO
*
*           computation of 'phi' functions and derivatives
*
            DO q=0,Z(i)
               DO j=0,q-1
                  rangT = rang(q-j,j+1,Z(i))
                  phi(q+1) = phi(q+1) + prodinv(rangT)
                  phif(q+1) = phif(q+1) - prodinv(rangT) * prodf(rangT)
                  phiT(q+1) = phiT(q+1) - prodinv(rangT) * prodT(rangT)
                  phiff(q+1) = phiff(q+1) + prodinv(rangT) *
     &                 (prodf(rangT)*prodf(rangT)-prodff(rangT))
                  phiTT(q+1) = phiTT(q+1) + prodinv(rangT) *
     &                 (prodT(rangT)*prodT(rangT)-prodTT(rangT))
                  phifT(q+1) = phifT(q+1) + prodinv(rangT) *
     &                 (prodf(rangT)*prodT(rangT)-prodfT(rangT))
               END DO
               DO j=q+1,Z(i)
                  rangT = rang(j-q,q+1,Z(i))
                  phi(q+1) = phi(q+1) + prod(rangT)
                  phif(q+1) = phif(q+1) + prod(rangT) * prodf(rangT)
                  phiT(q+1) = phiT(q+1) + prod(rangT) * prodT(rangT)
                  phiff(q+1) = phiff(q+1) + prod(rangT) *
     &                 (prodf(rangT)*prodf(rangT)+prodff(rangT))
                  phiTT(q+1) = phiTT(q+1) + prod(rangT) *
     &                 (prodT(rangT)*prodT(rangT)+prodTT(rangT))
                  phifT(q+1) = phifT(q+1) + prod(rangT) *
     &                 (prodf(rangT)*prodT(rangT)+prodfT(rangT))
               END DO
            END DO
*
*           computation of abundances
*
            DO q=1,Z(i)+1
               rangA = debutA + q - 1
               phinv = 1.0d0 / phi(q)
               abundZ(rangA,1) = Y(i) * phinv
               aphi2 = abundZ(rangA,1) * phinv
               abundZ(rangA,2) = -aphi2 * phiT(q)
               abundZ(rangA,3) = -aphi2 * phif(q)
               abundZ(rangA,4) = aphi2 * (2.0d0*phinv*phiT(q)*phiT(q)-
     &              phiTT(q))
               abundZ(rangA,5) = aphi2 * (2.0d0*phinv*phif(q)*phif(q)-
     &              phiff(q))
               abundZ(rangA,6) = aphi2 * (2.0d0*phinv*phif(q)*phiT(q)-
     &              phifT(q))
            END DO
            debutA = debutA + z(i) + 1
         END IF
      END DO
      END

************************************************************************
*     Solve CARREFULLY a second order polynomial equation              *
*                    ------ TESTED ------                              *
************************************************************************

      SUBROUTINE secondegree(a,b,c,nbr,nbc,reel,complexe,ecart,nr,sdel)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.cons'
      DOUBLE PRECISION a,b,c
      LOGICAL nr
*     Intermediate variables
      DOUBLE PRECISION ainv,b2,l,sdel,truncdiff
      DOUBLE PRECISION x,x2,xr,xi,x2r,x2i,d,fac,coeff,pow,dev,ac4
      LOGICAL a0,b0,c0
      INTEGER i,order,nbsol
      PARAMETER (order=4)
*     outputs
      INTEGER nbr,nbc
      DOUBLE PRECISION reel(4),complexe(4,2),ecart(4)
*     --------------------
      a0 = ABS(a).LT.zero
      b0 = ABS(b).LT.zero
      c0 = ABS(c).LT.zero
      IF (.NOT.a0) THEN
         ainv = 1.0d0 / a
      END IF
      IF (a0.AND.b0.AND.c0)THEN
         nbr = 0
         nbc = 0
      END IF
      IF (a0.AND.b0.AND.(.NOT.c0)) THEN
         nbr = 0
         nbc = 0
      END IF
      IF (a0.AND.c0.AND.(.NOT.b0)) THEN
         nbr = 1
         nbc = 0
         reel(1) = 0.0d0
      END IF
      IF (b0.AND.c0.AND.(.NOT.a0)) THEN
         nbr = 2
         nbc = 0
         reel(1) = 0.0d0
         reel(2) = reel(1)
      END IF
      IF (a0.AND.(.NOT.c0).AND.(.NOT.b0)) THEN
         nbr = 1
         nbc = 0
         reel(1) = -c / b
      END IF
      IF (b0.AND.(.NOT.c0).AND.(.NOT.a0)) THEN
         l = -c * ainv
         IF (l.LT.0.0d0) THEN
            nbr = 0
            nbc = 2
            l = DSQRT(-l)
            complexe(1,1) = 0.0d0
            complexe(1,2) = l
            complexe(2,1) = 0.0d0
            complexe(2,2) = -l
         ELSE
            nbr = 2
            nbc = 0
            reel(1) = DSQRT(l)
            reel(2) = -reel(1)
         END IF
      END IF
      IF (c0.AND.(.NOT.b0).AND.(.NOT.a0)) THEN
         nbr = 2
         nbc = 0
         reel(1) = 0.0d0
         reel(2) = -b * ainv
      END IF
      IF ((.NOT.a0).AND.(.NOT.b0).AND.(.NOT.c0)) THEN
         b2 = b * b
         ac4 = 4.0d0 * a * c
         d = ac4 / b2
         l = 0.5d0 * ainv
         IF (d.GT.1.0d0) THEN
            nbr = 0
            nbc = 2
            sdel = DSQRT(truncdiff(ac4,b2))
            complexe(1,1) = -l * b
            complexe(1,2) = l * sdel
            complexe(2,1) = complexe(1,1)
            complexe(2,2) = -complexe(1,2)
         ELSE
            nbr = 2
            nbc = 0
            sdel = DSQRT(truncdiff(b2,ac4))
            IF (ABS(d).LT.numericalimit) THEN
               fac = 1.0d0
               coeff = 1.0d0
               pow = 0.5d0
               dev = 0.0d0
               DO i=1,order
                  fac = fac * i
                  coeff = coeff * (pow-i+1)
                  dev = dev + (-1)**i * coeff * d**i / fac
               END DO
               reel(1) = b * dev * l
               reel(2) = -b * ainv * (1.0d0+0.5d0*dev)
            ELSE
               reel(1) = l * (-b+sdel)
               reel(2) = l * (-b-sdel)
            END IF
         END IF
      END IF
      IF (nr) THEN
         DO i=1,nbr
            call improvereal(0.0d0,0.0d0,a,b,c,reel(i),x)
            reel(i) = x
         END DO
      END IF
      nbsol = 0
      DO i=1,nbr
         nbsol = nbsol + 1
         x = reel(i)
         x2 = x * x
         ecart(nbsol) = ABS(a*x2+b*x+c)
      END DO
      DO i=1,nbc
         nbsol = nbsol + 1
         xr = complexe(i,1)
         xi = complexe(i,2)
         x2r = xr * xr - xi * xi
         x2i = xr * xi + xr * xi
         ecart(nbsol) = ABS(a*(x2r+x2i)+b*(xr+xi)+c)
      END DO
      END

************************************************************************
*     Solve CARREFULLY a third order polynomial equation               *
************************************************************************

      SUBROUTINE thirdegree(a,b,c,d,nbr,nbc,reel,complexe,ecart,nr)
*     --------------------
      IMPLICIT NONE
*     inputs
      DOUBLE PRECISION a,b,c,d
      LOGICAL nr
*     Intermediate variables
      DOUBLE PRECISION q,r,zero,ainv,ainv2,inv3,sqr3,racub1(2,2),q3
      DOUBLE PRECISION racub,q3inv,pi,b3inv,l,m,dum
      DOUBLE PRECISION x,x2,xr,xi,x2r,x2i,x3r,x3i
      LOGICAL r0,q0,a0,d0
      PARAMETER (zero=1.d-90)
      INTEGER i,nbsol
*     outputs
      INTEGER nbr,nbc
      DOUBLE PRECISION reel(4),complexe(4,2),ecart(4),ecasec(4)
*     --------------------
      nbr = 0
      nbc = 0
      a0 = ABS(a).LT.zero
      d0 = ABS(d).LT.zero
      IF (a0) THEN
         CALL secondegree(b,c,d,nbr,nbc,reel,complexe,ecasec,nr,dum)
      END IF
      IF (d0) THEN
         CALL secondegree(a,b,c,nbr,nbc,reel,complexe,ecasec,nr,dum)
         nbr = nbr + 1
         reel(nbr) = 0.0d0
      END IF
      IF (.NOT.a0.AND.(.NOT.d0)) THEN
         inv3 = 1.0d0 / 3.0d0
         sqr3 = DSQRT(3.0d0)
         racub1(1,1) = -0.5d0
         racub1(1,2) = sqr3 * 0.5d0
         racub1(2,1) = -0.5d0
         racub1(2,2) = -sqr3 * 0.5d0
         ainv = 1.0d0 / a
         ainv2 = ainv * ainv
         b3inv = b * inv3
         q = c * ainv - b * b3inv * ainv2
         r = ainv * d + ainv2 * b3inv * (2.0d0*b3inv*b3inv*ainv-c)
         q0 = ABS(q).LT.zero
         r0 = ABS(r).LT.zero
         IF (r0) THEN
            IF (q0) THEN
               nbr = 3
               nbc = 0
               reel(1) = 0.0d0
               reel(2) = 0.0d0
               reel(3) = 0.0d0
            ELSE
               CALL secondegree(1.0d0,0.0d0,q,nbr,nbc,reel,complexe,
     &              ecasec,nr,dum)
               nbr = nbr + 1
               reel(nbr) = 0.0d0
            END IF
         ELSE
            IF (q0) THEN
               nbr = 1
               nbc = 2
               reel(1) = racub(-r)
               complexe(1,1) = reel(1) * racub1(1,1)
               complexe(1,2) = reel(1) * racub1(1,2)
               complexe(2,1) = reel(1) * racub1(2,1)
               complexe(2,2) = reel(1) * racub1(2,2)
            ELSE
               q3 = -(q**3)/27.0d0
               CALL secondegree(1.0d0,r,q3,nbr,nbc,reel,complexe,ecasec,
     &              nr,dum)
               IF (nbc.EQ.0) THEN
                  nbr = 1
                  nbc = 2
                  l = racub(reel(1))
                  m = racub(reel(2))
                  reel(1) = l + m
                  complexe(1,1) = l * racub1(1,1) + m * racub1(2,1)
                  complexe(1,2) = l * racub1(1,2) + m * racub1(2,2)
                  complexe(2,1) = l * racub1(2,1) + m * racub1(1,1)
                  complexe(2,2) = l * racub1(2,2) + m * racub1(1,2)
               ELSE
                  nbr = 3
                  nbc = 0
                  pi = 3.1415926535d0
                  q3inv = -q * inv3
                  l = DSQRT(4.0d0*q3inv)
                  m = inv3 * DACOS(-0.5d0*r*q3inv**(-1.5d0))
                  reel(1) = l * DCOS(m)
                  reel(2) = l * DCOS(m+2.0d0*pi*inv3)
                  reel(3) = l * DCOS(m+4.0d0*pi*inv3)
               END IF
            END IF
         END IF
         nbsol = 0
         DO i=1,nbr
            nbsol = nbsol + 1
            reel(i) = reel(i) - inv3 * b * ainv
            IF (nr) THEN
               call improvereal(0.0d0,a,b,c,d,reel(i),x)
               reel(i) = x
            END IF
            x = reel(i)
            x2 = x * x
            ecart(nbsol) = ABS(a*x*x2+b*x2+c*x+d)
         END DO
         DO i=1,nbc
            nbsol = nbsol + 1
            xr = complexe(i,1) - inv3 * b * ainv
            complexe(i,1) = xr
            xi = complexe(i,2)
            x2r = xr * xr - xi * xi
            x2i = xr * xi + xr * xi
            x3r = x2r * xr - x2i * xi
            x3i = x2r * xi + x2i * xr
            ecart(nbsol) = ABS(a*(x3r+x3i)+b*(x2r+x2i)+c*(xr+xi)+d)
         END DO
      END IF
      END

************************************************************************
*     Cubic function also for negative arguments                       *
************************************************************************

      DOUBLE PRECISION FUNCTION racub(x)
      IMPLICIT NONE
      DOUBLE PRECISION x,inv3
      inv3 = 1.0d0 / 3.0d0
      IF (x.GE.0.0d0) THEN
         racub = x**inv3
      ELSE
         racub = -((-x)**inv3)
      END IF
      END

***********************************************************************
* Transform derivatives versus f and T into derivatives versus         *
* rho and T, or contrary, depending on logical SENS.                   *
* IN/OUT: derivatives of rho versus f and T (eoscom.rho)               *
*         dxdf and dxdT: derivatives of one quantitie expressed in the *
*          set of variables (f,T)                                      *
*         dydrho and dydT: the same quantitie expressed in the set of  *
*          variables (rho,T)                                           *
************************************************************************

      SUBROUTINE transdftodrho(sens,dxdf,dxdT,dydrho,dydT,rhoT,rhof,
     &rhofinv)
*     --------------------
      IMPLICIT NONE
*     inputs
      LOGICAL sens
      DOUBLE PRECISION dxdf,dxdT
      DOUBLE PRECISION rhof,rhoT,rhofinv
*     Intermediate variables
      DOUBLE PRECISION truncdiff
*     outputs
      DOUBLE PRECISION dydrho,dydT
*     --------------------
      IF (sens) THEN
         dydrho = dxdf * rhofinv
         dydT = truncdiff(dxdT,dydrho*rhoT)
      ELSE
         dxdf = dydrho * rhof
         dxdT = dydT + dydrho * rhoT
      END IF
      END

************************************************************************
* Carefully difference                                                 *
************************************************************************

      DOUBLE PRECISION FUNCTION truncdiff(x,y)
*     -------------------
      IMPLICIT NONE
      include 'evolpar.star'
      include 'eoscom.cons'
      DOUBLE PRECISION x,y
      LOGICAL egal
*     -------------------
      IF (egal(x,y)) THEN
         truncdiff = numericalimit
      ELSE
         truncdiff = x - y
      END IF
      END

************************************************************************
* Truncated exponential between exp(-150) and exp(150)                 *
************************************************************************

      DOUBLE PRECISION FUNCTION truncexpo(x)
      IMPLICIT NONE
      include 'evolpar.star'
      include 'eoscom.cons'
      DOUBLE PRECISION x
         truncexpo = EXP(MAX(-limexp,MIN(limexp,x)))
      END

************************************************************************
* Truncated logarithm at zero                                          *
************************************************************************

      DOUBLE PRECISION FUNCTION truncln(x)
      IMPLICIT NONE
      include 'evolpar.star'
      include 'eoscom.cons'
      DOUBLE PRECISION x
         truncln = DLOG(DABS(zero+x))
      END

************************************************************************
* Truncated inverse function at zero                                   *
************************************************************************

      DOUBLE PRECISION FUNCTION truncinv(x)
      IMPLICIT NONE
      include 'evolpar.star'
      include 'eoscom.cons'
      DOUBLE PRECISION x
         x = x + zero
         truncinv = 1.0d0 / x
      END

************************************************************************
* Partition function calculus for H                                    *
* All the partition function are approximate with the statistical      *
* weight of the fundamental state except H2                            *
*  Inputs: beta=1/kT (ev-1)                                            *
*  Outputs:partition function parfunH(i,j)                             *
*   with i=1 ==> H, i=2 ==> H+, i=3 ==> H- and i=4 ==> H2.             *
*   and j=1 ==>partition function Z, j=2,3==>dlnZ/dln and              *
*   j=4,5,6==>ddlnZ/dln2                                               *
************************************************************************

      SUBROUTINE ZH_sod(Tsh,lnTsh,beta,lnparfunH)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.ion'
      include 'eoscom.phys'
      include 'eoscom.math'
      include 'eoscom.par'
      include 'eoscom.cons'
      DOUBLE PRECISION lnTsh,beta,Tsh
*     Intermediate variables
      INTEGER i,j
*     ouputs
      DOUBLE PRECISION lnparfunH(nbionHmax,deriv3)
*     --------------------
      DO i=1,nbionHmax
         DO j=2,deriv3
            lnparfunH(i,j) = 0.0d0
         END DO
      END DO
*     --------------------
      IF (compH2.AND.(Tsh.LE.TmaxH2)) THEN
         IF (VardyaZH2) THEN
            CALL h2Vardya_sod(beta,lnparfunH)
         END IF
         IF (IrwinZH2) THEN
            CALL h2Irwin_sod(lnTsh,lnparfunH)
         END IF
      ELSE
         lnparfunH(4,1) = lng0H(4)
      END IF
*     --------------------
      IF (compZH0.AND.(Tsh.LE.Tmaxirwin)) THEN
         CALL h_sod(lnTsh,lnparfunH)
      ELSE
         lnparfunH(1,1) = lng0H(1)
      END IF
*     --------------------
      lnparfunH(2,1) = lng0H(2)
      lnparfunH(3,1) = lng0H(3)
      END

************************************************************************
* Partition function calculus for Z>1                                  *
* All the partition function are approximate with the statistical      *
* weight of the fundamental state                                      *
*  Inputs: beta=1/kT (erg-1)                                           *
*  Outputs:partition function parfunZ(nbion,j)                         *
*   with nbion=1,2,3... ==> He,He+,He++,C,C+,...                       *
*   and j=1 ==>partition function Z, j=2,3==>dlnZ/dln and ddlnZ/dln2   *
************************************************************************

      SUBROUTINE ZZ_sod(Tsh,lnTsh,lnparfunz)
*     --------------------
      IMPLICIT NONE
*     inputs
      include 'evolpar.star'
      include 'eoscom.math'
      include 'eoscom.ion'
      include 'eoscom.phys'
      include 'eoscom.par'
      include 'eoscom.cons'
      DOUBLE PRECISION lnTsh,Tsh
*     Intermediate variables
      INTEGER i,qeff,rangG,finG,rangA,debut,indx,j
*     ouputs
      DOUBLE PRECISION lnparfunZ(nbionZmax,deriv3)
*     ---------Z=statistical weigth of fundamental state----------
      debut = 0
      DO i=2,nbz-1
         IF (ionized(i)) THEN
            finG = indx(2,z(i)-1)
            DO qeff=1,z(i)+1
               rangA = debut + qeff
               rangG = finG + qeff
               lnparfunZ(rangA,1) = lng0Z(rangG)
               DO j=2,deriv3
                  lnparfunZ(rangA,j) = 0.0d0
               END DO
            END DO
*           Calcul des fonction de partition de C, N et O...reference ?
*            ==> voir Manuel.
            IF (compZCNO) THEN
               IF ((z(i).EQ.6).OR.(z(i).EQ.7).OR.(z(i).EQ.8)) THEN
                  CALL cno_sod(Tsh,z(i),debut,lnparfunZ)
               END IF
            END IF
            IF (compZHe.AND.(z(i).EQ.2).AND.(Tsh.LE.Tmaxirwin)) THEN
               CALL he_sod(lnTsh,debut,lnparfunZ)
            END IF
            debut = debut + Z(i) + 1
         END IF
      END DO
      END



c$$$*********************************************************************
c$$$      SUBROUTINE zmoy(xsp)
c$$$* Test TD. Nov 2018
c$$$      IMPLICIT NONE
c$$$*     Number of isotopes (nsp) and chemical elements computed (nsp)
c$$$      include 'evolpar.star'
c$$$*     maximum number of shell (nsh)
c$$$      include 'evolcom.teq'
c$$$*     Independant variables (lnf and T)
c$$$      include 'evolcom.var'
c$$$*     Z of the differents isotopes (znuc)
c$$$      include 'evolcom.nuc'
c$$$*     Parameters
c$$$c      include 'eoscom.par'
c$$$
c$$$      INTEGER i,j
c$$$
c$$$*     Current Shell number
c$$$      Integer shell,amu,boltz
c$$$      PARAMETER (amu= 1.6605402d-24,boltz=1.3807d-16)
c$$$
c$$$      DOUBLE PRECISION NI(nsh),zpart(nsh),nbrelec(nsh),xsp 
c$$$      DIMENSION xsp(nsh,nis)
c$$$
c$$$********
c$$$      
c$$$      do i=1,nmod
c$$$         NI(i) = 0.5
c$$$         nbrelec(i) = (lnf(i)/(amu)) * (xsp(i,ihe4)
c$$$     $        /anuc(ihe4)) * 24.59 !*(T(i) / 11605*boltz)
c$$$         zpart(i) = nbrelec(i)/NI(i)
c$$$         if (i.eq.300) then
c$$$            print *, nbrelec(i),zpart(i), lnf(i),xsp(i,ihe4),T(i),NI(i)
c$$$            stop
c$$$         endif
c$$$      enddo
c$$$      return
c$$$      end
