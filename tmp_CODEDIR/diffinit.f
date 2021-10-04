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
c         muizero(i) = muizero(i)+xsp(i,2)/anuc(2)+xsp(i,5)/anuc(5) 
         do j = 1,nis   
c           if (j.ne.io16) muizero(i) = muizero(i)+xsp(i,j)/anuc(j)
            muizero(i) = muizero(i)+xsp(i,j)/anuc(j)  ! 1 / mui physique
         enddo        
         muizero(i) = 1.d0/muizero(i)    ! mui
         khi(i) = 3.02427122d-4*t(i)**3/(kap(i)*ro(i))
      enddo
c      write(286,*) model
c      DO i = 1,nmod
c         write(286,'(1x,i4,1x,1pe16.6,8(1x,1pe16.6))') i,r(i),
c     &        muizero(i)
c      ENDDO
c      call sgsmooth (muizero,nmod,1,17,21)
c     call sgsmooth (muizero,nmod,1,16,16)

c      id = 5
c      nmodmdix = nmod-10
c      do k = 1,10
c         call lissage (muizero,nmod-10,10,id)
c      enddo
c      write(287,*) model
c      DO i = 1,nmod
c         write(287,'(1x,i4,1x,1pe16.6,8(1x,1pe16.6))') i,r(i),
c     &        muizero(i)
c      ENDDO
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
