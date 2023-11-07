      SUBROUTINE thermo (error)

************************************************************************
* Compute thermodynamical quantities and optical depth                 *
*                                                                      *
* $LastChangedDate:: 2016-05-11 17:20:46 +0200 (Mer, 11 mai 2016)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 62                                                          $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.eng'
      include 'evolcom.eos'
      include 'evolcom.ion'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.opa'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer error
      integer neff1
      integer j,k,l,iopa,iopaf

      integer i1,i2,i3,i4

      integer ktmax,ksh
      double precision tmax

      double precision kap1,kapt1,kapr1,kapx1,kapy1,engnup,engnupdt,
     &     engnupdro,eng,engdt,engdro,rst,rstdf,rstdt,taueffl,taune1l,
     &     taunel,rne1l,rnel,tne1l,tnel,lne1l,lnel,deltnel,deltes,drl
      double precision reffcorr,teffcorr,irrad
      double precision zzmean,amean,xamean
      double precision x(nsp),y(nsp)
      
*______________________________________________________________
***   calculation of the local metallicity and molar abundances
*--------------------------------------------------------------

      do l = 1,nsp
         do k = 1,nmod
            ysp(k,l) = xsp(k,l)/anuc(l)
         enddo
      enddo

*_____________________________________________________________
***   calculation of the equation of state, and of the related
***   thermodynamic quantities (+ derivatives)
*-------------------------------------------------------------
      
      call eos (1,error)
      if (error.gt.0) return
            
*_________________________________________
***   calculation of the screening factors
*-----------------------------------------

      if (iter.le.1) call fscreen (t,ro,xsp,rhpsi,anuc,znuc)

*_________________________________________________________________
***   calculation of H2 equilibrium abundance in case of accretion
*-----------------------------------------------------------------

      if (dmaccr.gt.0.d0.and.xspacc(ih2).gt.1.d-10.and.accphase.eq.0)
     &     call ydequil

*_______________________________
***   local variable definitions
*-------------------------------

      ktmax = 0
      tmax = -1.d0
      do k = 1,nmod
         if (t(k).gt.tmax) then
            ktmax = k
            tmax = t(ktmax)
         endif
      enddo

      iopa = 0
      iopaf = 0
      do k = 1,nmod
         ksh = k
         do l = 1,nis-1
            x(l) = max(ymin*anuc(l),xsp(k,l))
            y(l) = ysp(k,l)
         enddo
         x(nis) = 1.d-50
         x(nsp) = 1.d-50
         y(nis) = 1.d-50
         y(nsp) = 1.d-50
         tconv(k) = 1.d99

*_________________________________________________
***   calculation of the total opacity coefficient
*-------------------------------------------------

         call kappa (t(k),ro(k),muiinv(k),x,ksh,kap1,kapr1,kapt1,
     &        kapx1,kapy1,nretry,iter,iopa,iopaf,error)

         if (error.gt.0) return
         kap(k) = kap1
         dkapdro(k) = kapr1
         dkapdf(k) = kapr1*drodf(k)
         dkapdt(k) = kapt1*t(k)+kapr1*drodt(k)
         dkapdx(k) = kapx1
         dkapdy(k) = kapy1

*_______________________________________________________
***   calculation of the energy loss by plasma neutrinos
*-------------------------------------------------------

         if (lnucl) then
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

*______________________________________________________________
***   calculation of the rate of produced energy by the nuclear
***   reactions and of the losses due to associated neutrinos
*--------------------------------------------------------------

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

            call denucl (t(k),ro(k),mueinv(k),y,ksh,eng,engdt,engdro)

            if (iter.le.1.and.ktmax.gt.0.and.imodpr.eq.11.and.
     &           ireset.eq.0) call myengen(k,ktmax)

            denucldro(k) = engdro-engnupdro
            denucldf(k) = denucldro(k)*drodf(k)
            denucldt(k) = (engdt-engnupdt)*t(k)+denucldro(k)*drodt(k)
            enucl(k) = eng-enupla(k)
         endif
      enddo

*___________________________________________________
***   smooth opacity profile at molecular transition
*---------------------------------------------------

      if (iopa.gt.0) then
         do k = iopa-10,nmod
            kap(k) = log(kap(k))
         enddo
         call lissage (kap,min(nmod1-1,iopa+10),iopa-3,5)
         call lissage (dkapdf,min(nmod1-1,iopa+10),iopa-3,5)
         call lissage (dkapdt,min(nmod1-1,iopa+10),iopa-3,5)
         do k = iopa-10,nmod
            kap(k) = exp(kap(k))
         enddo
      endif

*______________________________________________________________
***   smooth opacity profile at OPAL <--> Fergusson transition
*--------------------------------------------------------------
* This smoothing leads to a step in temperature at this opacity
* transition
* Commented out on 04/02/2014

c$$$      if (iopaf.gt.0) then
c$$$         do k = iopaf-10,nmod
c$$$            kap(k) = log(kap(k))
c$$$         enddo
c$$$         call lissage (kap,min(nmod1-1,iopaf+10),iopaf-3,5)
c$$$         call lissage (dkapdf,min(nmod1-1,iopaf+10),iopaf-3,5)
c$$$         call lissage (dkapdt,min(nmod1-1,iopaf+10),iopaf-3,5)
c$$$         do k = iopaf-10,nmod
c$$$            kap(k) = exp(kap(k))
c$$$         enddo
c$$$      endif

*_____________________________________________
***   calculation of the optical depth profile
*---------------------------------------------
      
      tau(nmod) = tau0
      dtaudf(nmod) = 0.d0
      dtaudt(nmod) = 0.d0
      l = nmod
!     Comments : Printing tau profile for debugging (13/06/23)
c$$$      write(*,*) ""
c$$$      write(*,*) "Tau profile generation :"
c$$$      write(*,"(1X,A,6X,6(A,9X))") "sh","tau","rst","kappa","rho","T",
c$$$     &     "drl"
c$$$      write(*,*) "--------------------------------------------------"//
c$$$     &     "--------------------------"
c$$$      write(*,"(1X,I4,6(1X,E11.5))") l, tau(l), rst, kap(l), ro(l), 
c$$$     &        t(l), drl
      do l = nmod1,2,-1
         drl = abs(r(l+1)-r(l))
         rst = kap(l)*ro(l)*drl
         rstdf = (kap(l)*drodf(l)+ro(l)*dkapdf(l))*drl
         rstdt = (kap(l)*drodt(l)+ro(l)*dkapdt(l))*drl
         tau(l) = tau(l+1)+rst
         dtaudf(l) = dtaudf(l+1)+rstdf
         dtaudt(l) = dtaudt(l+1)+rstdt
c$$$         write(*,"(1X,I4,6(1X,E11.5))") l, tau(l), rst, kap(l), ro(l), 
c$$$     &        t(l), drl
      enddo
      tau(1) = tau(2)
      dtaudf(1) = dtaudf(2)
      dtaudt(1) = dtaudt(2)
      l = 1
c$$$      write(*,"(1X,I4,6(1X,E11.5))") l, tau(l), rst, kap(l), ro(l), 
c$$$     &        t(l), drl
c$$$      write(*,*) "--------------------------------------------------"//
c$$$     &     "--------------------------"

*____________________________________________
***   calculation of the effective quantities
*--------------------------------------------

      neff = 0
      do k = nmod-2,2,-1
         if (tau(k).gt.taueff.and.tau(k+1).ge.taueff) exit
      enddo
      neff = k+2
      neff1 = neff-1
      taueffl = log(abs(taueff))
      taune1l = log(abs(tau(neff1)))
      taunel = log(abs(tau(neff)))
      rne1l = lnr(neff1)
      rnel = lnr(neff)
      tne1l = lnt(neff1)
      tnel = lnt(neff)
      lne1l = log(abs(lum(neff1)))
      lnel = log(abs(lum(neff)))
      deltnel = tne1l-tnel
***   interpolate quatities at tau = 2/3
      tphot = tne1l-(taune1l-taueffl)*deltnel/(taune1l-taunel)
      deltes = tne1l-tphot
      tphot = exp(tphot)
      if (deltnel.ne.0.d0) then
         reff = rne1l-deltes*(rne1l-rnel)/deltnel
         leff = lne1l-deltes*(lne1l-lnel)/deltnel
         reff = exp(reff)
         leff = exp(leff)
      else
         leff = lum(neff)
         reff = r(neff)
      endif
      teff = (leff/(pim4*sig*reff*reff))**0.25d0
      
*_____________________________________________________________
***   corrected atmospheric temperature profile for optically
***   thick winds (WR phase)
*-------------------------------------------------------------

c      if (xsp(nmod1,ih1).le.0.2d0.and.log10(teff).gt.4.d0.and.
c     &     m(nmod)/msun.gt.10.d0) then
c         call corrwind (teffcorr,reffcorr)
c         teff = 10.d0**teffcorr
c         reff = reffcorr
c      endif
      
*__________________________________________
***   calculation of the internal structure
*------------------------------------------

      call structure (error)
      
      return
      end
