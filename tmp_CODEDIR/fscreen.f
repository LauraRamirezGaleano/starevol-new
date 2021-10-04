

************************************************************************

      SUBROUTINE fscreen (t,ro,xsp,rhpsi,anuc,znuc)

************************************************************************
* Calculation of the screening factors                                 *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.ion'
      include 'evolcom.teq'

      integer i,ii,k,l

      double precision t(nsh),ro(nsh),xsp(nsh,nsp),rhpsi(nsh),anuc(nsp),
     &     znuc(nsp)
      double precision fscr
      double precision screen
      double precision dzi
      double precision y(nsp),zx,zy,zbar,ztilde,l0
      double precision tburn,tnucmin,fnetdt,fnetdt0,tolnuc,tolnuc0
      double precision nelectronk,zeta,h0,kfact,hfact,zstar,kscreen

      common /funcscr/ fscr(nsh,nscr)
      common /nuc_param/ tburn,tnucmin,fnetdt,fnetdt0,tolnuc,tolnuc0

      external screen

      kfact = pim4*ech**2/boltz
      hfact = -0.4d0*pi*kfact

      fscr = 1.d0


***   Charged reaction factors
      do k = 1,nmod
         nelectronk = ro(k)*avn*mueinv(k)
         if (t(k).gt.tnucmin) then
            do l = 1,nsp
               y(l) = max(1.d-50,xsp(k,l))/anuc(l)
            enddo
            zx = 0.d0
            zbar = 0.d0
            ztilde = 0.d0            
            do i = 1,nis
               zx = zx+y(i)
               zy = znuc(i)*y(i)
               zbar = zbar+zy
               ztilde = ztilde+znuc(i)*zy
            enddo
c..   note zx = muiinv(k)
            zstar = ztilde/zbar
            kscreen = dsqrt(kfact*nelectronk*(zstar+rhpsi(k))/t(k))
            zeta = 3.d0*kscreen**3/(pim4*nelectronk)
            h0 = hfact*nelectronk**2/(t(k)*kscreen**5)
            zbar = zbar/zx
            ztilde = dsqrt(ztilde/zx+rhpsi(k)*zbar)
            l0 = 1.88d8*dsqrt(ro(k)*zx/t(k)**3)
c..   nbz-1 because nbz = heavy
            do i = 1,nbz-1
               ii = nbz-1+i
               dzi = dble(i)
               fscr(k,i) = screen (y,zbar,ztilde,zeta,l0,h0,dzi,1.d0)
               fscr(k,ii) = screen (y,zbar,ztilde,zeta,l0,h0,dzi,2.d0)
            enddo
            fscr(k,2*nbz-1) = screen(y,zbar,ztilde,zeta,l0,h0,6.d0,6.d0)
            fscr(k,2*nbz) = screen (y,zbar,ztilde,zeta,l0,h0,8.d0,8.d0)
            fscr(k,2*nbz+1) = screen(y,zbar,ztilde,zeta,l0,h0,8.d0,6.d0)
         endif
      enddo

      return
      end
