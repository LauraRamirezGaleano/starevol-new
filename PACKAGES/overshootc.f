******************************************************************************
*  overshootc.f gives the routine for penetrative convectivon / overshoot
*  Herwig
*  Baraffe et al. 2017
*  Augustson & Augustson 2019
*  Korre et al. 2019
*  Version 1.0 (July 2021)
*  Thibaut Dumont
*******************************************************************************

      SUBROUTINE diffherwig (istart,iend,Dmlt,Dherw,fover,pwi,icz)

************************************************************************
*   Compute the diffusion coefficient associated with diffusive        *
*   (Herwig) overshooting                                              *
*   icz =  1 : downward overshoot                                      *
*   icz = -1 : upward overshoot                                        *
*                                                                      *
*   Dmlt = conv. diffusion coef at the base/top of the convective zone *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.mass'
      include 'evolcom.nuc'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      double precision tnuc9,vrate,fscr,t9
      double precision hvcd,factor,fover,hv,redge,zover,fcross
      double precision Dmlt,cd0,cd1,cd01,Dherw,pwi
      double precision dt9l,weight,tdest,vdest,drdiff,tdiffus,rtau

      integer jvit
      integer istart,iend,imin,imax,idest,icz
      integer i,j

      common /funcscr/ fscr(nsh,nscr)
      common /nucvit/ tnuc9(nvit),vrate(nvit,nre)

      dimension Dherw(nsh)


      if (fover.le.0.d0.or.Dmlt.le.0.d0) then
         iend = istart
         return
      endif
      cd0 = log(Dmlt)
      cd1 = log(1.d-2)
      cd01 = cd0-cd1

***   gravity waves
      if (icz.eq.2) then
         hv = hp(istart)
         redge = r(istart)
         i = istart
         do  while (r(i).le.0.5d0*hv+redge)
            i = i+1
         enddo
         i = max(istart,i-1)
         hvcd = 1.6d-9*hp(i)*cd01*r(nmod)/sconv(i)
         factor = 0.d0
         i = istart-1
c         do while (factor.lt.1.d0.and.crz(i).gt.0)
         do while (factor.lt.1.d0)
            zover = (redge-r(i))**2
            factor = zover/(hvcd*hp(i))
            Dherw(i) = cd01*(1.d0-factor)+cd1
            Dherw(i) = exp(Dherw(i))
            crz(i) = 1
            i = i-1
            if (i.eq.1) then
               iend = 1
               return
            endif
         enddo
         iend = i+1
      endif 
         

***   exponential diffusion coefficient
***   pwi = 1 corresponds to Herwig's case
***   cd1 defines the cut in cd: for ln(cd).lt.cd1 cd(i) = 0.d0

      if (icz.eq.1) then
         hv = fover*hp(istart)
         redge = r(istart)
         hvcd = 2.d0/(hv*cd01)
         factor = 0.d0
         i = istart-1
c         do while (factor.lt.1.d0.and.crz(i).gt.0)
         do while (factor.lt.1.d0)
            zover = abs(redge-r(i))
            factor = zover*hvcd
            Dherw(i) = cd01*(1.d0-factor)**pwi+cd1
            Dherw(i) = exp(Dherw(i))
            sconv(i) = sconv(istart)
            tconv(i) = (r(i+1)-r(i))**2/(Dherw(i)+1.d-50)
c            tconv(i) = hp(i)**2/(Dherw(i)+1.d-50)
            crz(i) = 1
            i = i-1
            if (i.eq.1) then
c               print *,' Pb diffherwig : fover must be too large'
               iend = 1
               return
            endif
         enddo
         iend = i+1
      endif

c..   upward overshoot
      if (icz.eq.-1) then
         hv = fover*hp(istart)
         redge = r(istart)
         hvcd = 2.d0/(hv*cd01)
         factor = 0.d0
         i = istart+1
         do while (factor.lt.1.d0.and.crz(i).gt.0)
            zover = abs(redge-r(i))
            factor = zover*hvcd
            Dherw(i) = cd01*(1.d0-factor)**pwi+cd1
            Dherw(i) = exp(Dherw(i))
            sconv(i) = sconv(istart)
            crz(i) = 1
            tconv(i) = (r(i+1)-r(i))**2/(Dherw(i)+1.d-50)
c            tconv(i) = hp(i)**2/(Dherw(i)+1.d-50)
            i = i+1
            if (i.eq.nmod1) then
c               print *,' Pb diffherwig : fover must be too large'
               iend = nmod1
               return
            endif
         enddo
         iend = i+1
      endif


ccccccccccccccccccccccc
      return
ccccccccccccccccccccccc

***   computation of the proton destruction time-scale by 12C(p,g)13N

      imin = min(iend,istart)
      imax = max(iend,istart)
      drdiff = (r(imin+1)-r(imin))**2
      tdiffus = drdiff/Dherw(imin+1)
      t9 = t(imin)*1.d-9
      do j = 1,nvit
         if (t9.lt.tnuc9(j)) goto 5
      enddo
 5    jvit = j-1
      dt9l = log10(t9/tnuc9(jvit))
      weight = dt9l/(log10(tnuc9(jvit+1)/tnuc9(jvit)))
      vdest = vrate(jvit,icpg)+weight*(vrate(jvit+1,icpg)-
     &     vrate(jvit,icpg))
      vdest = 10.d0**vdest*ro(imin)*fscr(imin,6)+1.d-99
      tdest = 1.d0/(vdest*ysp(imin,ic12))

      rtau = 5.d0*mtini*mtini
      if (zkint.ge.5.d-5) rtau = rtau*dsqrt(5.d-5/zkint)
      if (tdest.lt.rtau*tdiffus) then
         if (icz.eq.-1) then
            write (nout,*) 'WARNING : top of convective zone too hot !'
            return
         endif
         idest = imin+1
         do i = imin+1,imax-1
            drdiff = (r(i+1)-r(i))**2
            t9 = t(i)*1.d-9
            do j = 1,nvit
               if (t9.lt.tnuc9(j)) goto 6
            enddo
 6          jvit = j-1
            dt9l = log10(t9/tnuc9(jvit))
            weight = dt9l/(log10(tnuc9(jvit+1)/tnuc9(jvit)))
            vdest = vrate(jvit,icpg)+weight*(vrate(jvit+1,icpg)-
     &           vrate(jvit,icpg))
            vdest = 10.d0**vdest*ro(i)*fscr(i,6)+1.d-99
            tdest = 1.d0/(vdest*ysp(i,ic12))
c        print *,i,tdest/3.1557807d7,drdiff/Dherw(i)/3.1557807d7
            if (tdest.ge.rtau*drdiff/Dherw(i)) goto 12
         enddo
 12      idest = i
         if ((imax-idest).gt.3) then
            fcross = fover*(r(imax)-r(idest))/(r(imax)-r(imin))
            write (nout,1000) fover,imin,fcross,idest
            write (90,1000) fover,imin,fcross,idest
            hv = fcross*hp(imax)
            hvcd = 2.d0/(hv*cd01)
            do i = imax,idest,-1
               zover = abs(redge-r(i))
               factor = zover*hvcd
               Dherw(i) = cd01*(1.d0-factor)**pwi+cd1
               Dherw(i) = exp(Dherw(i))
            enddo
         else
            idest = imax
            write (nout,1001) tdiffus,tdest
            write (90,1001) tdiffus,tdest
         endif
         do i = imin,idest-1
            Dherw(i) = 0.d0
         enddo
         imin = idest
      endif


 1000 format (1x,'Diffusion in H-burning region! fover reduced: ',
     &     1pe10.4,' [shell ',i4,'] ---> ',1pe10.4,' [shell ',i4,']')
 1001 format (1x,'Base of convective envelope too hot! No diffusion!',
     &     /,1x,' tdiffus = ',1pe8.2,' >> tdest = ',1pe8.2)

      return
      end

************************************************************************

      SUBROUTINE diffbaraffe (istart,iend,Dmlt,Dbar,fover)

************************************************************************
*   Compute the diffusion coefficient associated with penetrative      *
* convection as formalised in Eq(1) of Baraffe et al. 2017, ApJL 845, L6*
*   icz =  1 : downward overshoot                                      *
*   icz = -1 : upward overshoot                                        *
*                                                                      *
*   Dmlt = conv. diffusion coef at the base/top of the convective zone *
*                                                                      *
*   fover used to limit the diffusion below the CE                     *
* $LastChangedDate:: 2019-10-28 11:33:00 +0100 (Mon, 28 Oct 2019)    $ *
* $Author:: Palacios                                                 $ *
* Modif: T.D. 2020                                                     *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.mass'
      include 'evolcom.nuc'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      double precision depth,fover,hv,redge,zover
      double precision Dmlt,Dbar,muBa,lamBa

      integer jvit
      integer istart,iend
      integer i,j


      dimension Dbar(nsh)


      if (fover.le.0.d0.or.Dmlt.le.0.d0) then
         iend = istart
         return
      endif
      muBa = 5.d-3
      lamBa = 6.d-3
      Dbar(1:nmod) = 0.d0

      print *,'in diffbaraffe'

***   exponential diffusion coefficient
***   the depth is limited by the fover factor which is a fraction of a 
***   pressure scale height

      hv = fover*hp(istart)
      redge = r(istart)
      i = istart!-1     !modif TD 07/01/20
      do while (r(istart)-r(i).lt.hv)
         print *,'depth, hv',depth,hv
c      do while (i.gt.1) 
         zover = (redge-r(i))
         Dbar(i) = Dmlt*(1.d0-exp(-exp(-(zover/r(nmod)-muBa)/lamBa)))
         sconv(i) = sconv(istart)
         tconv(i) = (r(i+1)-r(i))**2/(Dbar(i)+1.d-50)
c         crz(i) = 1   ! Modif TD 09/12/2019
         i = i-1
         if (i.eq.1) then
            iend = 1
            return
         endif
      enddo
      iend = i+1
      
      return
      end


      
************************************************************************

      SUBROUTINE DiffKyle(istart,icz,iend,Dmlt,fover,Dkyle,itype)

************************************************************************
*   Weibel distribution from Augustson & Mathis 2019                   *
*   rotating convection model of Augustson & Mathis 2019,              *
*                                                                      *
*     Author: Kyle Augustson (CEA Saclay)                              *
*     Date of creation: 26 October 2019.                               *
*     Adaptation to Starevol : TD (04/11/2019)                         *
************************************************************************     
!     istart beginning of overshoot region
!     iend end of overshoot region
!     Dmlt mixing length diffusion coefficient at istart
!     Dkyle the coefficient
!     fover overshoot coefficient
!     itype = 1 -> Augustson & Mathis 2019
!     itype = 2 -> Korre et al. 2019 (adapted, see Dumont et al. 2021a)
************************************************************************

      
      Implicit None

      include 'evolpar.star'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.mass'
      include 'evolcom.nuc'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      Integer, Intent(In) :: istart,icz,itype
      Integer, Intent(InOut) :: iend
      Integer :: i
      Real*8, Intent(In) :: Dmlt, fover
      Real*8, Intent(InOut) :: Dkyle(nsh)
      Real*8 :: hv,redge,factor,zz,eps_kyle

      double precision omega
      common /rotvar/ omega(nsh)

      if ((fover.le.0.d0).or.(Dmlt.le.0.d0)) then
         iend = istart
         Dkyle = 0d0
         return
      endif
      
!     Weibel double exponential diffusion coefficient
      call Quintic(omega(istart),sconv(istart),alpha_mlt_hp*hp(istart),
     &     zz)
      factor = (5d0/2d0)**(0.25d0)
      redge = r(istart)

      if (itype.eq.1) then
         hv = factor*fover*hp(istart)/zz**0.75d0
!     Diffusion below a convection zone
         if (icz.eq.1) then
            i = istart          !-1   !modif TD 27/11/19
            factor = 1d-2
            eps_kyle = 1d-16
            zz = 1d32
            do while (zz.gt.factor)
               zz = 5d0/6d0-dabs(r(i)-redge)/hv
               zz = 1d0-dexp(-dexp(zz))
               zz = Dmlt*zz
               sconv(i) = sconv(istart)
               tconv(i) = (r(i+1)-r(i))**2/(zz+eps_kyle)
c            crz(i) = 1    ! modif TD 29/11/19
c            i = i-1
c            if (i.eq.1) then
c               iend = 1
c               return
c            end if
               Dkyle(i) = zz
c            print *, i,crz(i),'Dkyle',Dkyle(i)
               i = i-1          ! modif TD 27/11/19
               if (i.eq.1) then
                  iend = 1
                  return
               endif
            enddo
            iend = i+1
         endif

!     Diffusion above a convection zone
         if (icz.eq.-1) then
            i = istart+1
            factor = 1d-2
            eps_kyle = 1d-16
            zz = 1d32
            do while (zz.gt.factor)
               zz = 5d0/6d0-dabs(r(i)-redge)/hv
               zz = 1d0-dexp(-dexp(zz))
               zz = Dmlt*zz
               sconv(i) = sconv(istart)
               tconv(i) = (r(i+1)-r(i))**2/(zz+eps_kyle)
c            crz(i) = 1
c            i = i+1
c            if (i.eq.1) then
c               iend = 1
c               return
c            end if
               Dkyle(i) = zz
               i = i+1
               if (i.eq.1) then
                  iend = 1
                  return
               endif
            enddo
            iend = i-1
         endif
      else if (itype.eq.2) then
         hv = factor*fover*hp(istart)/dsqrt(zz)
         hv = 2d0*hv*hv
         !     Diffusion below a convection zone
         if (icz.eq.1) then
            i = istart          !-1   !modif TD 27/11/19
            factor = 1d-2
            eps_kyle = 1d-16
            zz = 1d32
            do while (zz.gt.factor)
               zz = -(r(i)-redge)**2/hv
               zz = dexp(zz)
               zz = Dmlt*zz
               sconv(i) = sconv(istart)
               tconv(i) = (r(i+1)-r(i))**2/(zz+eps_kyle)
c     crz(i) = 1    ! modif TD 29/11/19
c     i = i-1
c     if (i.eq.1) then
c     iend = 1
c     return
c     end if
               Dkyle(i) = zz
c     print *, i,crz(i),'Dkyle',Dkyle(i)
               i = i-1          ! modif TD 27/11/19
               if (i.eq.1) then
                  iend = 1
                  return
               endif
            enddo
            iend = i+1
         endif

!     Diffusion above a convection zone
         if (icz.eq.-1) then
            i = istart+1
            factor = 1d-2
            eps_kyle = 1d-16
            zz = 1d32
            do while (zz.gt.factor)
               zz = -(r(i)-redge)**2/hv
               zz = dexp(zz)
               zz = Dmlt*zz
               sconv(i) = sconv(istart)
               tconv(i) = (r(i+1)-r(i))**2/(zz+eps_kyle)
c     crz(i) = 1
c     i = i+1
c     if (i.eq.1) then
c     iend = 1
c     return
c     end if
               Dkyle(i) = zz
               i = i+1
               if (i.eq.1) then
                  iend = 1
                  return
               endif
            enddo
            iend = i-1
         endif
      endif
      return
      end

