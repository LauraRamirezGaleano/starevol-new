

************************************************************************

      SUBROUTINE oversh

************************************************************************
*   Calculate parametric overshooting : artificial extension of CZ     *
*   Also computes the diffusion coefficient associated with            *
*   convection in the overshooted layers (needed in case of rotation)  *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.conv'
      include 'evolcom.cons'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer imin,imax
      integer lmin,lmax
      integer kl,kl1,kdw,kup,nc
      integer i,j,k

      double precision dmconv,rext,rext1,rext2,dr,dxm
      double precision lambdac

      dimension lambdac(nsh)

!     Commented out VOJE -- allow novopt > 2 to trigger only core overshoot (11.03.2024)
!     if (novopt.eq.3) novopt = 1
!     if (novopt.eq.4) novopt = 2
      do kl = 1,nsconv
         imin = novlim(kl,1)
         imax = novlim(kl,2)
***   treat core overshoot only
         if (imin.gt.1.and.novopt.gt.2) goto 30

***   overshoot below CZ
         dmconv = mr(imax)-mr(imin)
         rext = adwn*hp(imin)
         if (t(imin).lt.1.d6) then
            i = imin-1
            goto 10
         endif
         do i = imin,1,-1
            if ((r(imin)-r(i)).gt.rext.or.(mr(imin)-mr(i)).gt.
     &           dmconv*adwn)  goto 10
         enddo
 10      if (i.eq.1) i = 0
         lmin = i+1

***   overshoot top CZ
         rext1 = aup*hp(imax)
         if (t(imax).lt.1.d6) then
            i = imax+1
            goto 20
         endif
         rext2 = aup*r(imax)
         do i = imax,nmod
            dr = r(i)-r(imax)
            dxm = mr(i)-mr(imax)
            if (novopt.eq.1.or.novopt.eq.3.and.
     &           (dr.gt.rext1.or.dxm.gt.dmconv*aup)) goto 20
            if (novopt.eq.2.or.novopt.eq.4.and.
     &           (dr.gt.min(rext1,rext2))) goto 20
         enddo
 20      lmax = i-1

         if (imin.gt.lmin) then
            do j = lmin,imin-1
               sconv(j) = sconv(imin)
               tconv(j) = tconv(imin)
               fconv(j) = lum(j)/(pim4*r(j)*r(j))
               lambdac(j) = hp(j)*alphac
               Dconv(j) = sconv(j)*lambdac(j)/3.d0
               crz(j) = 1
            enddo
         endif
         if (lmax.gt.imax) then
            do j = imax+1,lmax
               sconv(j) = sconv(imax)
               tconv(j) = tconv(imax)
               fconv(j) = lum(j)/(pim4*r(j)*r(j))
               lambdac(j) = hp(j)*alphac
               Dconv(j) = sconv(j)*lambdac(j)/3.d0
               crz(j) = 1
            enddo
         endif
         novlim(kl,3) = lmin
         novlim(kl,4) = lmax
         novlim(kl,7) = lmin
         novlim(kl,8) = lmax
      enddo


***   redefine novlim if zones overlap
 30   nc = nsconv
      do kl = 1,nsconv-1
         kl1 = kl+1
         if (novlim(kl,4).ge.novlim(kl1,3)) then
            kdw = kl
            if (novlim(kl1,3).lt.novlim(kl,3)) kdw = kl1
            kup = kl1
            if (novlim(kl,4).gt.novlim(kl1,4)) kup = kl
            novlim(kl,3) = novlim(kdw,3)
            novlim(kl,4) = novlim(kup,4)
            novlim(kl,7) = novlim(kdw,7)
            novlim(kl,8) = novlim(kup,8)
            do k = kl1+1,nc
               novlim(k-1,3) = novlim(k,3)
               novlim(k-1,4) = novlim(k,4)
               novlim(k-1,7) = novlim(k,3)
               novlim(k-1,8) = novlim(k,4)
            enddo
            novlim(nc,3) = 0
            novlim(nc,4) = 0
            novlim(nc,7) = 0
            novlim(nc,8) = 0
            nc = nc-1
            if (nc.eq.nsconv-1) goto 40
         endif
      enddo
 40   nsconv = nc

      return
      end
