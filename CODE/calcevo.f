      SUBROUTINE calcevo (*)

****************************************************************************
*     Compute the next model and, in case of success, save the results     *
*     If convergence fails, the procedure is repeated by changing the      *
*     current time-step                                                    *
*     09/03: Simplification of the procedure of re-reading initial model   *
*     in case of crashes ( see also modification of rinimod)               *
*     Modifs CC ondes (2/04/07 --> 23/11/07)                               *
*                                                                          *
*     $LastChangedDate:: 2016-05-13 11:10:12 +0200 (Ven, 13 mai 2016)    $ *
*     $Author:: palacios                                                 $ *
*     $Rev:: 73                                                          $ *
*                                                                          *
****************************************************************************
#ifdef GYRE
c........................Rajoute par Alice
      use evolstell_gyre
#endif
      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.atm'
      include 'evolcom.cons'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.ion'
      include 'evolcom.lum'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.opa'
      include 'evolcom.rot'
      include 'evolcom.conv'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.var'

#ifdef GYRE     
c......................Rajoute par Alice
      include 'evolcom.grad'
      include 'omp_lib.h'
#endif      
      integer imerk,icrasht,iflag
      integer ntprof0,itermax0
      integer neff1, nconvsh, ienvsh
      integer error,iexp
      integer i,j,k,l,kl,ksh
      integer klenv,klpulse,klcore,iHecore
      integer imin,imax,itmaxf,itmax,ishockb,ishockt
      integer fcs               ! first convective shell

      double precision alpha01,alpha00,alphmin0
      double precision ledd,lnuc,vlnuc,dabslum,dloglum,dfnuc
      double precision sum1,taueffl,taune1l,taunel,rne1l,rnel,tne1l,
     &     tnel,rone1l,ronel,lne1l,lnel,deltnel,deltes,dtna,tmaxf,tmax,
     &     Lmax,Mackmax,enucmax,renucl,vi
      double precision zzmean,amean,xamean,yeq,y0,x,y,irrad
      double precision eng,engdt,engdro,engnup,engnupdt,engnupdro,engro
     &     ,engrocap
      double precision elemsave,mconv
      double precision menvconv,dmenvconv,mcore
      double precision omegacritic,Rpolinv
      double precision mHecore
      double precision geffrot,Rp2
      double precision times,chronos ! Ajout chronos pour mesure du tps de calcul

      logical partialmix,pass
      logical vlconv,adapt

      character month*36,GetDateTimeStr*15,fmt*26
      parameter (month = 'JanFebMarAprMayJunJulAugSepOctNovDec')
      parameter (fmt = '(I2.2,A1,I2.2,I3,A3,I4)')
      integer itime(8)

c................Rajoute par Alice
      logical :: res, resdet
      integer, dimension(nsh) :: normx
      logical, dimension(nsh) :: mask
      double precision brunt
      double precision bruntv(nsh), bruntv2(nsh)      
      common /brunt/ brunt(nsh),bruntv2,bruntv
      double precision xmHb, xmHt
      common /hburn/ xmHb, xmHt
      logical rgbphase_gyre
c.......................Fin rajout
      
      common /saveMS/ elemsave
      common /convergence/ imerk,icrasht
      common /grid/ alpha00,alpha01,alphmin0,ntprof0,itermax0
      common /overshoot/ klcore,klenv,klpulse
      common /resolution/ vlconv
      common /hydrodyn/ Lmax,tmax,Mackmax,enucmax,itmax,ishockb,ishockt
      common /meshconv/ menvconv,dmenvconv,mcore,nconvsh,ienvsh
      common /geffective/ geffrot(nsh),Rp2(nsh)
c     common /eddington/ ledd

      dimension renucl(nsh),vi(nreac),y0(nsp),x(nsp),y(nsp)


      
      error = 0
c       if (nphase.eq.4.and.novlim(1,3).eq.1) then
c        iHecore = novlim(1,4)
c        mHecore = m(iHecore)
c       else
c        iHecore = 1
c        mHecore = 0.d0
c       endif
      model = model+1
      maxsh = abs(maxsh0)
      if (model.le.0) time = 0.d0
      if (no.eq.0) then
         itermax0 = itermax
         alpha00 = alpha0
         alpha01 = alpha
         ntprof0 = ntprof
         alphmax0 = alphmax
         alphmin0 = alphmin
      endif
      if (model.le.3.and.no.eq.2) then
         itermax = itermax0
         alphmax = alphmax0
         alphmin = alphmin0
         alpha = alpha01
         alpha0 = alpha00
         ntprof = ntprof0
      endif
      
      if (imodpr.eq.11.and.ireset.eq.0) call mycompo

*------------------------------
***   some specific adaptations
*------------------------------

      if (icrasht.ne.0) then
         icrash = 9
         ireverse = 1
         write (nout,300)
         write (90,300)
      else
         icrash = abs(icrash0)
         if (icrash0.lt.0) ireverse = -1
      endif

c..   before AGB phase
      if (nphase.le.4.and..not.rgbphase.and..not.igwrot)
     &     maxmod = min(maxmod,300) ! modif Fev 2019

c..   AGB phase

c..   increase structural time step & constrain the number of saved models
c..   during the pulse
      if (ifail.eq.0) then
         if (thermalpulse.and.imodpr.ne.0) then
            fts = min(fts0,2.d-2)
c..   pulse starts, save immediately
            if (maxmod.eq.maxmod0) then
               if (model-modeli.gt.10) then
                  maxmod = no+1
               else
                  maxmod = min(11-mod(model-modeli,10)+no,maxmod0)
               endif
               write (nout,600) fts,maxmod
            else
               write (nout,650) fts
            endif
         else
            fts = fts0
         endif
c..   further constrain structural time step during 3DUP
         if ((icorr.eq.'t'.or.icorr.eq.'a').and.dup3) then
c     if (dup3) then
            fts = min(fts0,2.d-2)
            if (imodpr.eq.33.and.maxmod.eq.maxmod0) then
               maxmod = min(41-mod(model-modeli,40)+no,maxmod0)
               write (nout,700) fts,maxmod
            else
               write (nout,750) fts
            endif
         else
            fts = fts0
         endif
      endif
c..   thibaut
c     if (diffherwig.and.agbphase.and.totm.lt.10.1d0.and.dup3.and.
c     &     dtn*seci.gt.0.1d0) then
c     dtn = 0.1d0*sec
c     write (nout,108) dtn*seci
c     write (90,108) dtn*seci
c     endif

      ifail = 0
      no = no+1
      
*______________________________
***   abundance renormalization
*------------------------------

c..   index of maximum for each row
      normx(:) = maxloc(xsp, 2)
c..   loop over rows of xsp
      do k = 1,nmod
c..   compute sum the row      
         sum1 = sum(xsp(k, 1:nsp))
c..   for the element of the row which is the maximum,
c     remove 1+(sum of the row) so that sum of the row = 1 now         
         xsp(k, normx(k)) = maxval(xsp(k,1:nsp))+(1.d0-sum1)
      enddo

*__________________________________________________________________
***   extrapolation of the next model and variable initializations
*------------------------------------------------------------------

      call guess
      
*______________________________________________________________________
***   calculation of the factor ft and fp for the treatment of rotation
*----------------------------------------------------------------------

      if (hydrorot) then
         call potential
      else
         do k = 1,nmod
            ft(k) = 1.d0
            fp(k) = 1.d0
         enddo
      endif
      
*___________________________________________________
***   calculation of the stellar structure evolution
*---------------------------------------------------

      if (ireset.eq.0) alpha = alpha0
      
      call nwtraf (*20,error)
      
*________________________
***   test He-core growth
*------------------------
c     if (nphase.eq.4.and.xsp(novlim(1,4),ihe4).gt.
c     &     vxsp(novlim(1,4),ihe4)) then
c     write(nout,*)'Breathing pulse developing : central He4',
c     &        ' abundance increased during this timestep'
c     error = 35
c     goto 20
c     endif

c     if (nphase.eq.4.and.mHecore.gt.0.d0.and.xsp(1,ihe4).lt.0.5d0.and.
c     &     xsp(1,ihe4).gt.0.03d0) then
c     if ((mHecore-m(novlim(1,4))).gt.0.1d0*msun.and.
c     &        (novlim(2,3)-novlim(1,4)).gt.5) then
c     c         if (abs(m(novlim(1,4))-mHecore).gt.0.1d0*msun) then
c     write (nout,75) m(novlim(1,4))/msun,novlim(1,4),
c     &           mHecore/msun,iHecore
c     if (xsp(1,ihe4)/vxsp(1,ihe4).gt.0.001d0) then
c     write (nout,76) vxsp(1,ihe4),novlim(1,3),xsp(1,ihe4)
c     $           ,novlim(1,3)
c     error = 35
c     goto 20
c     endif
c     endif

*______________________________
***   test temperature increase
*------------------------------

!     Setting unused values of T to zero to avoid bugs (06/12/2023)
      do i=1,nsh
         if (i.gt.nmod) then
            t(i) = 0.d0
         endif
      enddo

      if (ftst.gt.0.d0.and.model.gt.10) then
         mask(:)=m.gt.mcut     
c..   Apply mask for t corresponding to m(i)<mcut except for t(1)
c     with merge function   
c..   Find tmax for elements if m(index_element)<mcut
         tmaxf=maxval([t(1),merge(t(2:), 0.d0, mask(2:))])
c..   Find index of tmax     
         itmaxf=maxloc([t(1),merge(t(2:), 0.d0, mask(2:))],1)
         if (abs(log10(tmax/tmaxf)).gt.log10(1.d0+ftst).and.max(tmax,
     &        tmaxf).gt.2.d7) then
            icrasht = icrasht+1
            if (icrasht.le.iresmax) then
               write (nout,80) int(ftst*100.d0),tmax,itmax,m(itmax)/msun
     &              ,tmaxf,itmaxf,m(itmaxf)/msun
               error = 40
               ireset = max(ireset-1,0)
               goto 20
            else
               stop 'Temperature increase too large : modify ftst'
            endif
         endif
      endif

      
      call thermo (error)


c     if (imodpr.eq.-1) then
c     call initondes
c     write (*,'(" mv sortie_ondes sortie_ondes_",i5)') model
c     stop
c     endif
      if (error.gt.0) goto 20



*____________________________________________________________________
***   
***   beyond this point the STRUCTURE has converged otherwise goto 20
***   
*--------------------------------------------------------------------


      icrasht = 0
      maxsh = abs(maxsh0)
      mixopt = mixopt0
c     print *,coefssv(:)
c     coefssv(:) = coefsv(:)
c     print *,coefssv(:)

*____________________________________________________________
***   calculation of the energy loss by plasma neutrinos
***   nuclear energy production rates and neutron irradiation
***   in case lnucl = .false.
*------------------------------------------------------------

      if (.not.lnucl) then
         do k = 1,nmod
            do l = 1,nis-1
               x(l) = max(ymin*anuc(l),xsp(k,l))
               y(l) = ysp(k,l)
            enddo
            x(nis) = 1.d-50
            x(nsp) = 1.d-50
            y(nis) = 1.d-50
            y(nsp) = 1.d-50
            tconv(k) = 1.d99
            if (t(k).gt.1.d7) then
               zzmean = 0.d0
               xamean = 0.d0
               do j = 1,nsp
                  zzmean = zzmean+znuc(j)*y(j)
                  xamean = xamean+y(j)
               enddo
               amean = 1.d0/xamean
               zzmean = zzmean*amean
               call nulosc (t(k),ro(k),amean,zzmean,engnup,engnupdt,
     &              engnupdro)
            else
               engnup = 0.d0
               engnupdt = 0.d0
               engnupdro = 0.d0
            endif
            enupla(k) = engnup

c..   calculation of the nuclear energy production rates
c..   neutrino energy losses and neutron irradiation

            irrad = 7.68771025912336d0*ysp(k,1)*ro(k)*dtn*dsqrt(t(k))
            if (irrad.le.1.72d0) then
               signt(k) = 43.8372d0*irrad+10.4d0
            else
               if (irrad.lt.2.92d0) then
                  signt(k) = -69.8417d0*irrad+205.93d0
               else
                  signt(k) = 1.99d0
               endif
            endif
            ksh = k
            call denucl (t(k),ro(k),mueinv(k),y,ksh,eng,engdt,engdro)
            enucl(k) = eng-enupla(k)

         enddo
      endif
      
*_______________________________________________________________
***   computation of vHconv and of the new convective boundaries
*---------------------------------------------------------------

      vhconv(1) = 0.d0
      do i = 2,nmod
         vhconv(i) = 0.d0
         tm(i) = t(i-1)**wi(i)*t(i)**wj(i)
      enddo
      tm(1) = t(1)
      rom(1) = ro(1)
      if (vlconv) call convzone
      do kl = 1,nsconv
         imin = novlim(kl,3)
         imax = novlim(kl,4)
         if (vlconv) call mlt (imin,imax,error)
         do i = imin,imax
            vhconv(i) = pim4*r(i)*fconv(i)/(tm(i)*rom(i))
         enddo
      enddo

***   compute parametric overshoot (alpha * Hp)
      if (novopt.ge.1) call oversh

*________________________________________________________________
***   calculation of the associated nucleosynthesis and diffusion
*----------------------------------------------------------------

c     if (ondes) call ondes
      if (idiffcc.and.agbphase.and..not.idiffcc0) then
         write (nout,100) idiffcc,idiffcc0
         idiffcc = idiffcc0
      endif

      if (nucreg.ne.3) then
         klcore = 0
         klpulse = 0
         klenv = 0

***   Compute mixing coefficients and nucleosynthesis
         pass = .true.
c     x         if (flame.and.diffzc) pass = .false.

c..   CASE 1 : diffusion (if nmixd > 0, diffusion treated in netdiff)
c..   
C     Modif CC ondes (19/10/07)
         if ((microdiffus.or.rotation.or.diffusconv.or.difftacho.or.
     &        thermohaline.or.igw).and.nmixd.eq.0.and.pass) then
C     &        thermohaline).and.nmixd.eq.0.and.pass) then
C     <--
!     Tests activation diffusion atomique ou non (modif Td Fev.2019)
            if (microdiffus.and.nphase.gt.5) then
               microdiffus = .false.
c     print *, 'phase after 5 microdiffus false'
            else if (nphase.eq.1.and.time*seci.lt.2e7) then ! Ajout test TD Juillet.2019 si rotation
c     else if (nphase.eq.1.and.time*seci.lt.2e6) then ! Ajout test TD AP Fev.2019
               microdiffus = .false.
c     print *,'no diffusion this young star'
            else if (nphase.eq.1.and.novlim(1,3).eq.1.and.nsconv.eq.1)
     $              then   
               microdiffus = .false.
c     print *,'microdiffus false - no radiative core'
            endif
!     Fin test activation diffusion atomique
            if (idiffnuc.eq.2) then
               if (nucreg.ne.0) call chedif (error)
               if (error.gt.0) goto 20
               call diffusion (1,error)
               if (error.gt.0) goto 20
            else
               call diffusion (1,error)
               if (error.gt.0) goto 20
               if (nucreg.ne.0) call chedif (error)
               if (error.gt.0) goto 20
            endif
            if (idiffnuc.eq.3) then
               call diffusion (1,error)
               if (error.gt.0) goto 20
            endif

         else

c..   CASE 2 : no diffusion or full coupling
c..   
            if (nucreg.ne.0) then
               if (nsconv.gt.0) then
                  if (novlim(1,3).eq.1) klcore = 1
                  do kl = 1,nsconv
                     if ((tau(novlim(kl,4)).lt.1.d2.or.t(novlim(kl,4))
     &                    .lt.1.d6).and.mr(novlim(kl,4)).gt.0.7d0.and.
     &                    klenv.eq.0) klenv = kl
                     if (t(novlim(kl,3)).gt.2.d8.and.novlim(kl,3).gt.1.
     &                    and.nphase.eq.5.and.klpulse.eq.0) klpulse = kl
                  enddo
                  if (nsconv.eq.1.and.klcore.eq.1) then
                     klenv = 0
                     klpulse = 0
                  endif
               endif
               call chedif (error)
               if (error.gt.0) goto 20
            else

c..   CASE 3 : no nucleosynthesis  (nucreg = 0)
c..   restore initial abundances and do mixing only
               do l = 1,nsp
                  do k = 1,nmod
                     xsp(k,l) = vxsp(k,l)
                     ysp(k,l) = xsp(k,l)/anuc(l)
                  enddo
               enddo
               do kl = 1,nsconv
                  imin = novlim(kl,3)
                  imax = novlim(kl,4)
                  call mix (dm,imin,imax,kl,0,partialmix)
               enddo
            endif
         endif
      else
c..   nucleosynthesis computed during convergence process
c..   simply update composition
         do l = 1,nsp
            do k = 1,nmod
               vxsp(k,l) = xsp(k,l)
            enddo
         enddo
      endif


***   check that diffusion and nucleosynthesis have a small energetic impact
***   and compute neutron equilibrium abundances

      lnuc = 0.d0
      vlnuc = 0.d0
      if (no.eq.1.and.nuclopt.eq.'n') write (90,*) '*** Compute final ',
     &     'neutron equilibrium abundance after nucleosynthesis'
      do k = 1,nmod
         do l = 1,nsp
            x(l) = xsp(k,l)
            y(l) = xsp(k,l)/anuc(l)
         enddo
         vlnuc = vlnuc+enucl(k)*dm(k)
         call vit (t(k),ro(k),mueinv(k),k,vi,0)
         if (nuclopt.eq.'n') then
            do l = 1,nsp
               y0(l) = x(l)/anuc(l)
            enddo
            call yneutron (y0,vi,yeq)
            xsp(k,1) = yeq
            vxsp(k,1) = yeq
            x(1) = yeq
         endif
         ksh = k
         call nuceng (-1,y,vi,ksh,eng,engro,engrocap)
         renucl(k) = eng-enupla(k)
         lnuc = lnuc+renucl(k)*dm(k)
      enddo
      if (nuclopt.eq.'m') then
         do kl = 1,nsconv
            imin = novlim(kl,3)
            imax = novlim(kl,4)
            yeq = 0.d0
            mconv = 0.d0
            do i = imin,imax
               mconv = mconv+dm(i)
               yeq = yeq+xsp(i,1)*dm(i)
            enddo
            yeq = yeq/mconv
            do i = imin,imax
               xsp(i,1) = yeq
            enddo
         enddo
      endif

      dloglum = abs(log(abs(vlnuc/lnuc)))
      dabslum = min(exp(dloglum)-1.d0,1.d6)
      dfnuc = log(1.d0+ftnuc)
      if (nphase.gt.1.and.model.gt.10.and.
     &     dloglum.gt.min(log(2.d0),dfnuc)) then
         if (dloglum.gt.dfnuc) then
            write (90,1810) model,min(9999,int(dabslum*1.d2)),
     &           int(ftnuc*1.d2)
            write (nout,1810) model,min(9999,int(dabslum*1.d2)),
     &           int(ftnuc*1.d2)
            error = 30
c     mixopt = .true.
         else
            error = -30
         endif
         if (dloglum.gt.1.d0.and.nmixd.eq.0.and.nphase.gt.2.and..not.
     &        thermalpulse.and.(nuclopt.eq.'c'.or.nuclopt.eq.'u')) then
            error = 30
            tolnuc = max(tolnuc0*1.d2,1.d-6)
c..   increase number of saved models in case nmixd = 5
            if (imodpr.eq.33) then
               maxmod = min(11-mod(model-modeli,10)+no,maxmod0)
            else
               maxmod = min(51-mod(model-modeli,50)+no,maxmod0)
            endif
            write (90,1850) model,nmixd0,maxmod0,maxmod
            write (nout,1850) model,nmixd0,maxmod0,maxmod
            nmixd = 5
         endif
         if (error.gt.0) goto 20
      endif
c..   save additional models if H is modified by more than 0.05 during MS
c..   at least 10 models computed
      if (xsp(1,ih1).lt.elemsave.and.nphase.eq.2) then
         maxmod = min(model-modeli,maxmod0)
         maxmod = max(maxmod,10)
         write (nout,2000)
         write (90,2000)
      endif
c..   idem for He during core He burning
      if (xsp(1,ihe4).lt.elemsave.and.nphase.eq.4) then
         maxmod = min(model-modeli,maxmod0)
         maxmod = max(maxmod,10)
         write (nout,2100)
         write (90,2100)
      endif

c..   save more models in the Hertzsprung gap
c      if (nphase.eq.3.and..not.rgbphase) maxmod = min(maxmod,25) ! couper TD 03/2020

*_______________________________________________________________________
***   determination of the photosphere location, temperature and density
*-----------------------------------------------------------------------

      neff = 0
      do k = nmod-2,2,-1
         if (tau(k).gt.taueff.and.tau(k+1).ge.taueff) exit
      enddo
      neff = k+2

      neff = min(neff,nmod)
      neff1 = neff-1
      taueffl = log(abs(taueff))
      taune1l = log(abs(tau(neff1)))
      taunel = log(abs(tau(neff)))
      rne1l = lnr(neff1)
      rnel = lnr(neff)
      tne1l = lnt(neff1)
      tnel = lnt(neff)
      rone1l = log(abs(ro(neff1)))
      ronel = log(abs(ro(neff)))
      lne1l = log(abs(lum(neff1)))
      lnel = log(abs(lum(neff)))
      deltnel = tne1l-tnel
      tphot = tne1l-(taune1l-taueffl)*deltnel/(taune1l-taunel)
      deltes = tne1l-tphot
      reff = rne1l-deltes*(rne1l-rnel)/deltnel
      roeff = rone1l-deltes*(rone1l-ronel)/deltnel
      roeff = min(roeff,1.d0)
      roeff = max(roeff,-20.d0)
      leff = lne1l-deltes*(lne1l-lnel)/deltnel
      leff = exp(leff)
      tphot = exp(tphot)
      reff = exp(reff)
      roeff = exp(roeff)

c..   TEST ROTATION DECEMBRE 2014 For a rotating star, the expression of
c..   the Eddington luminosity is modified to account for the change in
c..   the gravitational potential that affects the mass

      if (hydrorot) then
         ledd = pim4*c*g*m(neff)/kap(neff)*(1.d0-pw23*vomega(neff)**2
     $        *Rp2(neff)**3/(g*m(neff)))
      else
         ledd = pim4*c*g*m(neff)/kap(neff)
      endif
      if (leff.gt.ledd) then
         write (nout,90) model,leff/lsun,ledd/lsun
         error = 66
         goto 20
c      nmod = neff
      endif

*______________________________________
***   calculation of surface quantities
*--------------------------------------
      if (hydrorot) then
         if (ntprof.eq.2) then
            geff = log10(geffrot(neff-1))
         else
            geff = log10(geffrot(neff))
         endif

      else
         geff = log10(g*m(neff)/(reff*reff))
      endif
      soll = leff/lsun
      solr = reff/rsun
c     teffy = tsurf
      teffy = teff
      if (time.lt.1.d0) then
         write (nout,1100) model,min(soll,9999999.d0),min(9999.d0,solr),
     &        teffy,totm,dtn,time
      else
         write (nout,1200) model,min(soll,9999999.d0),min(9999.d0,solr),
     &        teffy,totm,dtn*seci,time*seci
      endif

*__________________________________________________________
***   mixing atmosphere abundances with convective envelope
*----------------------------------------------------------

      if (nsconv.gt.0.and.tau(novlim(nsconv,4)).lt.1.d2) then
         imin = novlim(nsconv,4)
         do i = 1,nsp
            do j = imin,nmod
               xsp(j,i) = xsp(imin-1,i)
               ysp(j,i) = xsp(j,i)/anuc(i)
            enddo
         enddo
      endif


*___________________________
***   storage of the results
*---------------------------

      call cpu_time (timeend)
      cputime = min(9999.d0,(timeend-timebeg))

      call prvar

*___________________________
***   Rajoute par Alice
*---------------------------

#ifdef GYRE

      ! For first iteration
      if (model .eq. 1) then
         gyreprec = .false.
         tpgy = 0.d0
         lpgy = 0.d0
         deltat = 0.d0
      endif

      dtgy = dtn*seci
      mgy = totm
      tgy = t(nmod)
      lgy = lum(nmod)/lsun

      xspgy = xsp(1,2)

      rgbphase_gyre = rgbphase.and.xmHb/xmHt.ge.0.8d0
      
      iter_gyre = imodpr  ! by default
      n_time_gyre = model

      call evo_iter("criteria.in", res, resdet, nphase, rgbphase_gyre)
   
      if (res .or. resdet) then ! if gyre must run in summary or detail mode     

         !msun = 1.989d33
         !rsun = 6.9599d10
         !lsun = 3.86d33
         nversion_gyre = 101
         n_gyre = nmod
         mstar_gyre = totm*msun
         rstar_gyre = r(nmod)
         lstar_gyre = lum(nmod)
         age_gyre = time*seci

         if (allocated(omega_gyre)) then
            deallocate(omega_gyre)
         endif
         if (allocated(kapkt_gyre)) then
            deallocate(kapkt_gyre)
         endif
         if (allocated(kapkro_gyre)) then
            deallocate(kapkro_gyre)
         endif
         if (allocated(enuclt_gyre)) then
            deallocate(enuclt_gyre)
         endif
         if (allocated(enuclro_gyre)) then
            deallocate(enuclro_gyre)
         endif
         if (allocated(r_gyre)) then
            deallocate(r_gyre)
         endif
         if (allocated(mr_gyre)) then
            deallocate(mr_gyre)
         endif
         if (allocated(lum_gyre)) then
            deallocate(lum_gyre)
         endif
         if (allocated(rho_gyre)) then
            deallocate(rho_gyre)
         endif
         if (allocated(t_gyre)) then
            deallocate(t_gyre)
         endif
         if (allocated(p_gyre)) then
            deallocate(p_gyre)
         endif
         if (allocated(abla_gyre)) then
            deallocate(abla_gyre)
         endif
         if (allocated(bruntv2_gyre)) then
            deallocate(bruntv2_gyre)
         endif
         if (allocated(gamma1_gyre)) then
            deallocate(gamma1_gyre)
         endif
         if (allocated(deltaks_gyre)) then
            deallocate(deltaks_gyre)
         endif
         if (allocated(kap_gyre)) then
            deallocate(kap_gyre)
         endif
         if (allocated(enucl_gyre)) then
            deallocate(enucl_gyre)
         endif
         if (allocated(abad_gyre)) then
            deallocate(abad_gyre)
         endif


         allocate(omega_gyre(n_gyre))
         allocate(kapkt_gyre(n_gyre))
         allocate(kapkro_gyre(n_gyre))
         allocate(enuclt_gyre(n_gyre))
         allocate(enuclro_gyre(n_gyre))
         allocate(r_gyre(n_gyre))
         allocate(mr_gyre(n_gyre))
         allocate(lum_gyre(n_gyre))
         allocate(rho_gyre(n_gyre))
         allocate(t_gyre(n_gyre))
         allocate(p_gyre(n_gyre))
         allocate(abla_gyre(n_gyre))
         allocate(bruntv2_gyre(n_gyre))
         allocate(gamma1_gyre(n_gyre))
         allocate(deltaks_gyre(n_gyre))
         allocate(kap_gyre(n_gyre))
         allocate(enucl_gyre(n_gyre))
         allocate(abad_gyre(n_gyre))

         
         do i = 1,nmod
            r_gyre(i) = r(i)
            mr_gyre(i) = mr(i)*msun
            lum_gyre(i) = lum(i)
            p_gyre(i) = p(i)
            t_gyre(i) = t(i)
            rho_gyre(i) = ro(i)
            abla_gyre(i) = abla(i)
            bruntv2_gyre(i) = bruntv2(i)
            gamma1_gyre(i) = gamma1(i)
            abad_gyre(i) = abad(i)
            deltaks_gyre(i) = deltaks(i)
            kap_gyre(i) = kap(i)
            kapkt_gyre(i) = t(i)*dkapdt(i)
            kapkro_gyre(i) = ro(i)*dkapdro(i)
            enucl_gyre(i) = enucl(i)
            enuclt_gyre(i) = t(i)*denucldt(i)
            enuclro_gyre(i) = ro(i)*denucldro(i)
            omega_gyre(i) = omega(i)
         enddo

!     open(unit=32, file="verif.log", position=
!     & "append", status ="unknown")
         
!     write(32, *)
!     write(32, *) "Time iteration :", n_time_gyre
!     write(32, 5100) "mr_gyre", "lum_gyre",
!     & "bruntv2_gyre", "enucl_gyre", "omega_gyre",
!     & "dkapdt_gyre", 
!     & "dkapdro_gyre", "denucldt_gyre", "denucldro_gyre", 
!     & "deltaks_gyre"
!     write_loop : do j = 1, size(r_gyre)
!     write(32, 6000) mr_gyre(j), lum_gyre(j),
!     & bruntv2_gyre(j), enucl_gyre(j), omega_gyre(j), 
!     & kapkt_gyre(j), 
!     &   kapkro_gyre(j), enuclt_gyre(j),         
!     &   enuclro_gyre(j), deltaks_gyre(j)
!     end do write_loop

!     close(unit=32) 

!     5100     format(10(A25))
!     6000     format(10(E25.16E3))

         call evolstell_build_data
         call evolstell_run('evolgyre','detail', resdet)

         tpgy = tgy
         lpgy = lgy

      endif
#endif      
c....................................Fin rajout
      
      call resulpr (error)
      print *, 'error',error
      if (error.eq.-9) then
         write (90,*) '   time-dependent convection activated'
         write (nout,*) '   time-dependent convection activated'
      endif
      
      iacc = iacc0
      hydro = hydro0
      ifail = 0
      if (error.eq.-30) then
         write (90,1800) model,min(9999,int(dabslum*1.d2))
         write (nout,1800) model,min(9999,int(dabslum*1.d2))
         if (icorr.eq.'t') ifail = -1
      endif
      write (90,1900)

      if (no.lt.maxmod) then
         call flush (90)
         return 1
      endif

      call date_and_time (values=itime)
      write (GetDateTimeStr,FMT) itime(5),':',itime(6),itime(3),
     &     month(3*itime(2)-2:3*itime(2)),itime(1)


      write (90,*)
      write (90,*)
      write (90,'("stop normal : ",a15)') GetDateTimeStr
      write (90,*)
      write (nout,*)
      write (nout,*) "stop 'normal'"


      return



*____________________________________
***   Managing crash situations
*------------------------------------

 20   if (error.eq.18) then
         write (90,1400)
         write (nout,1400)
      else
         if (error.eq.66) then
            write (90,1666) model
            write (nout,1666) model
         endif
         if (error.eq.50) then
            write (90,1555) model
            write (nout,1555) model
         endif
         if (error.eq.99) then
            write (90,1300) model
            write (nout,1300) model
         endif
         if (error.eq.1) then
            write (90,3001)
            write (nout,3001)
         endif
         if (error.eq.2) then
            write (90,3002)
            write (nout,3002)
         endif
         if (error.eq.5) then
            write (90,3005)
            write (nout,3005)
         endif
         if (error.eq.6) then
            write (90,3006)
            write (nout,3006)
         endif
c..   OPACITY problem
         if (error.eq.7) then
            write (90,3007)
            write (nout,3007)
         endif
         if (error.eq.80) then
            write (90,3080)
            write (nout,3080)
         endif
c     if (error.eq.81) then
c     write (90,3081)
c     write (nout,3081)
c     endif
         if (error.eq.82) then
            write (90,3082)
            write (nout,3082)
         endif
         if (error.eq.83) then
            write (90,3083)
            write (nout,3083)
         endif
         if (error.eq.84) then
            write (90,3084)
            write (nout,3084)
         endif
c..   EOS problem
         if (error.eq.10) then
            write (90,3010)
            write (nout,3010)
         endif
         if (error.eq.11) then
            write (90,3011)
            write (nout,3011)
         endif
         if (error.eq.12) then
            write (90,3012)
            write (nout,3012)
         endif
c..   DIFFUSION problem
         if (error.eq.20) then
            write (90,3020)
            write (nout,3020)
         endif
         if (error.eq.21) then
            write (90,3021)
            write (nout,3021)
         endif
         if (error.eq.22) then
            write (90,3022)
            write (nout,3022)
         endif
         if (error.eq.23) then
            write (90,3023)
            write (nout,3023)
         endif
         if (error.eq.15) then
            write (90,3015)
            write (nout,3015)
         endif
         if (error.eq.16) then
            write (90,3016)
            write (nout,3016)
         endif
         if (error.eq.9) then
            write (90,3009)
            write (nout,3009)
         endif
c..   NUCLEAR problem
         if (error.eq.13) then
            write (90,3013)
            write (nout,3013)
         endif
         if (error.eq.25) then
            write (90,3025)
            write (nout,3025)
         endif
         if (error.eq.26) then
            write (90,3026)
            write (nout,3026)
         endif
         if (error.eq.27) then
            write (90,3027)
            write (nout,3027)
         endif
         if (error.eq.14) then
            write (90,3014)
            write (nout,3014)
         endif
         if (error.eq.3) then
            write (90,3003)
            write (nout,3003)
            maxsh = 0
         endif
         if (error.eq.4) then
            write (90,3004)
            write (nout,3004)
         endif
c$$$         if (error.eq.101) then
c$$$            write (90,3101)
c$$$            write (nout,3101)
c$$$         endif
      endif

      ifail = 1
      iexp = 2
      if (error.eq.40.or.error.eq.30) then
         dtn = dtn/(facdt**4)
      else
         if (ireverse.ge.0) then
            if (ireset.lt.abs(icrash)) dtn = dtn/(facdt**iexp)
            if (ireset.eq.abs(icrash)) dtn = dtn0*facdt
            if (ireset.gt.abs(icrash).and.ireset.lt.iresmax) dtn = dtn*
     &           facdt
         else
            if (ireset.lt.abs(icrash)) dtn = dtn*facdt
            if (ireset.eq.abs(icrash)) dtn = dtn0/(facdt**iexp)
            if (ireset.gt.abs(icrash).and.ireset.lt.iresmax) dtn = dtn/
     &           (facdt**iexp)
         endif
      endif
      ireset = ireset+1

      if (ireset.eq.iresmax) goto 60
      if (dtn.lt.dtmin) goto 50


***   in case of failure, adaptative parameter card changes

      if (icorr.ne.'f') then
         if (icorr.eq.'m'.or.icorr.eq.'i'.or.icorr.eq.'h'.or.
     &        icorr.eq.'a') then
            adapt = .false.
         else
            adapt = .true.
         endif

c..   after 1 crash, mesh is frozen
         if (ireset.ge.1.and.(icorr.eq.'m'.or.adapt))
     &        then
            maxsh = 0
c     if (ireset.eq.1) dtn = dtn0
            write (nout,400)
            write (90,400)
            ireverse = 1
         else
            maxsh = abs(maxsh0)
            if (maxsh0.lt.0) ireverse = -1
         endif

c..   after 2 crashes, convergence acceleration disabled
         if (ireset.ge.2.and.iacc0.and.iacc.and.(icorr.eq.'h'.or.
     &        adapt)) then
            iacc = .false.
            write (nout,500)
            write (90,500)
         endif

c..   after 2 crashes if end AGB, change tau0
c     if (ireset.eq.3.and.nphase.eq.5.and.ntprof.ne.2.and.(adapt.or.
c     &        icorr.eq.'h').and.(dmenvconv.lt.0.15d0*msun.or.
c     &        dmenvconv/mtini/msun.lt.0.1d0)) then
c     if (tau0.lt.5.d-4) then
c     tau0 = 1.d-3
c     else
c     tau0 = 1.d-4
c     endif
c     dtn = dtn0
c     write (nout,520) tau0
c     write (90,520) tau0
c     endif

c..   after 2 crashes, tolerence on acceleration reduced
         if (hydro.and.(icorr.eq.'i'.or.adapt)) then
            if (ireset.ge.2) then
               eps(1) = eps(1)*1.d1
               write (nout,550) eps(1)
               write (90,550) eps(1)
            else
               eps(1) = eps0(1)
            endif
            if (eps(1).gt.1.d2) then
               eps(1) = eps0(1)
               hydro = .false.
            endif
         endif
      endif

      dtna = dtn
      if (dtn.gt.1.d2) then
         write (90,1500) model,dtn*seci,alpha
         write (nout,1500) model,dtn*seci,alpha
      else
         write (90,1510) model,dtn,alpha
         write (nout,1510) model,dtn,alpha
      endif

      error = 0
      

c..   restore saved model and current time step
      rewind (92)
      rewind (94)
      iflag = ifail
      call rinimod (92,94,92,94,iflag)

      dtn = dtna
      no = no-1
      
      call flush (90)
      return 1

 50   write (nout,1600)
      write (90,1600)
      goto 70
 60   write (nout,1700) model
      write (90,1700) model

 70   write (90,*)
      write (90,*) "failed model:"
      write (90,*)

      write (90,*)
      write (90,*) "stop 'failed'"
      write (90,*)
      write (nout,*)

      return


 75   format (/,' He-core mass changes too much : ',0pf6.3,' [#',i4,
     &     '] --> ',0pf6.3,' [#',i4,']')
 76   format (/,' He-core abundance changes too much : ',0pf6.3,' [#',i4
     &     ,'] --> ',0pf6.3,' [#',i4,']')
 80   format (/,' Temperature increase exceeds ftst = ',i3,'% : T = ',
     &     1pe9.3,' [#',i4,' Mr = ',0pf5.3,'] --> ',1pe9.3,' [#',i4,
     &     ' Mr = ',0pf5.3,']')
 90   format (1x,'WARNING : model ',i8,', Lsurf (',1pe9.3,') > ',
     &     'Eddington limit (Ledd = ',1pe9.3,')')
 100  format (1x,'WARNING : parameter card change: idiffcc = ',l1,
     &     'set to',l1)
 300  format (1x,'WARNING : parameter card change: icrash = 9')
 400  format (1x,'WARNING : parameter card change: maxsh = 0')
 500  format (1x,'WARNING : parameter card change: iacc = .false.')
 520  format (1x,'WARNING : parameter card change: tau0 = ',0pf7.5,
     &     ' and time step reset')
 550  format (1x,'WARNING : parameter card change: u = ',1pe8.2)
 600  format (1x,' pulse : fts changed to ',1pe8.2,' and maxmod ',
     &     'decreased to ',i3)
 650  format (1x,'WARNING :  pulse : fts changed to ',1pe8.2)
 700  format (1x,'WARNING :  3DUP  : fts changed to ',1pe8.2,
     &     ' and maxmod decreased to ',i3)
 750  format (1x,'WARNING :  3DUP  : fts changed to ',1pe8.2)
 1100 format (/,3x,'model',5x,'L',8x,'R',6x,'Teff',7x,'M',9x,'dt (s)',
     &     8x,'t (s)',/,75('-'),/,i8,1x,f10.2,1x,f7.2,1x,f7.0,1x,
     &     f11.8,1x,1pe10.4,1x,1pe16.10,/)
 1200 format (/,3x,'model',5x,'L',8x,'R',6x,'Teff',7x,'M',8x,'dt (yr)',
     &     7x,'t (yr)',/,75('-'),/,i8,1x,f10.2,1x,f7.2,1x,f7.0,1x,
     &     f8.3,1x,1pe10.4,1x,1pe16.10,/)
 1300 format (/,10x,'model ',i8,' failed : too many iterations',/)
 1400 format (/,10x,'new model computation with a reduced time-step,',
     &     ' due to a convergence problem for angular momentum',
     &     ' transport',/)
 1500 format (/,1x,'retry model ',i8,': dtn = ',1pe11.5,
     &     ' yr, with alpha = ',0pf6.4,/)
 1510 format (/,1x,'retry model ',i8,': dtn = ',1pe11.5,
     &     ' s, with alpha = ',0pf6.4,/)
 1555 format (/,1x,'retry model ',i8,':Too large variation of Vsurf',/)
 1600 format (//,5x,'model ',i8,' failed : time-step too low !',/)
 1666 format  (/,1x,'retry model ',i8,'because L > Ledd',/)
 1700 format (//,5x,'model ',i8,' failed : maximum number of crashes ',
     &     ' reached',/)
 1800 format (/,1x,'WARNING : model ',i8,', large variation of the ',
     &     'TOTAL nuclear luminosity (',i4,'%) : time step should be ',
     &     'reduced                  !')
 1810 format (/,1x,'WARNING : model ',i8,', too large variation of ',
     &     'the TOTAL nuclear luminosity (',i4,'% > ftnuc = ',i3,'%),')
c     &     ' mixopt changed : t')
 1850 format (/,1x,'WARNING : model ',i8,', too large variation of ',
     &     'the TOTAL nuclear luminosity (> 500%), nmixd changed : ',
     &     i1,' --> 5, nmaxmod changed : ',i3,' --> ',i3)
 1900 format (132('='))
 2000 format ('save model because of too large H depletion')
 2100 format ('save model because of too large He depletion')

c..   error labels

 3001 format (3x,'MLT : too many convective zones')
 3002 format (3x,'max iterations reached for MLT with hro')
 3005 format (3x,'Pb in NWRMAT')
 3006 format (3x,'Corrections too big')
 3007 format (3x,'KAPPA : Mass fractions exceed unity')
 3080 format (3x,'KAPPA : pb table 2')
c     3081 format (3x,'KAPPA : pb table 3 - O rich')
 3082 format (3x,'KAPPA : pb table 3 - N rich')
 3083 format (3x,'KAPPA : pb table 4')
 3084 format (3x,'KAPPA : pb molecular opacity')
 3010 format (3x,'EOS : rho outside reasonable limits')
 3011 format (3x,'EOS : pb with ionization of H')
 3012 format (3x,'EOS : pb with ionization of heavy elements')
 3020 format (3x,'DIFFSOLVE : division by zero at center !')
 3021 format (3x,'DIFFSOLVE : division by zero in the interior')
 3022 format (3x,'DIFFSOLVE : division by zero at surface !')
 3023 format (3x,'DIFFSOLVE : too many iterations')
 3015 format (3x,'ELEMENTS_FINIS : no solution found')
 3016 format (3x,'DIFFUSION : tstep too large')
 3009 format (3x,'NWRMAT : pb with convective diffusion')
 3013 format (3x,'NETDIFF : maximum iterations and/or time-step ',
     &     'too low !')
 3025 format (3x,'NETDIFF : negative abundance')
 3026 format (3x,'NETDIFF : renormalisation too large')
 3027 format (3x,'NETDIFF : maximum iterations and/or non-convergence')
 3014 format (3x,'NUCSOLVE : maximum iterations and/or non-convergence')
 3003 format (3x,'Pb in NETWORK, maxsh set to ZERO')
 3004 format (3x,'Pb in CHEDIF')
! 3101 format (3x,'NWTRAF : NaN in xmod')

      end
