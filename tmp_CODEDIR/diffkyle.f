
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
      call Quintic(omega(istart),sconv(istart),alphac*hp(istart),zz)
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
