

************************************************************************

      DOUBLE PRECISION FUNCTION screen (y,zbar,ztilde,zeta,l0,h0,z1,z2)

************************************************************************
* Calculate the screening factor for the thermonuclear reaction rates  *
* Input : nuclide charges z1 and z2                                    *
* Graboske etal 1973,ApJ,181,457 and Cox & Guili 17.15 for definitions *
* For intermediate screening :                                         *
* Mitler 1977, ApJ 212, 513 and Cox & Guili 17.15 for definitions      *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.nuc'

      double precision z1,z2,z12
      double precision zetafac,zeta1,zeta2
      double precision ztilde,zbar,h120,l0,l12,h120a
      double precision zeta,h0
      double precision y,zx,zmean,l086

      dimension y(nsp)

      h120 = 0.d0
      l12 = l0*z1*z2*ztilde

c... Graboske et al
c *___________________
c ***   weak screening
c *-------------------
c       if (l12.lt.0.1d0) then
c          h120 = l12

c *___________________________
c ***   intermediate screening
c *---------------------------
c       elseif (l12.ge.0.1d0.and.l12.le.5.d0) then
c          zmean = 0.d0
c          do i = 1,nis
c             zmean = zmean+(znuc(i)**1.58d0*y(i))
c          enddo
c          zmean = zmean/zx
c          l086 = 0.38d0*l0**0.86d0
c          h120 = zmean*ztilde**(-0.58d0)*zbar**(-0.28d0)
c          h120 = h120*((z1+z2)**1.86d0-z1**1.86d0-z2**1.86d0)*l086
c       endif


*____________________________________
***   weak and intermediate screening
***   Mitler (1977, ApJ, 212, 513) 
*------------------------------------
      if (l12.le.5.d0) then
         zeta1 = zeta*z1
         zeta2 = zeta*z2
         zetafac = (zeta1+zeta2+1.d0)**pw53-(zeta1+1.d0)**pw53-(zeta2+
     &        1.d0)**pw53+1.d0
         h120 = abs(h0*zetafac)
      endif

*_____________________
***   strong screening
*---------------------
      if (l12.ge.2.d0) then
         z12 = z1+z2
         h120a = (z12**pw53-z1**pw53-z2**pw53)
     &        +.316d0*zbar**pw13*(z12**pw43-z1**pw43-z2**pw43)
     &        +0.737d0/(zbar*l0**pw23)*(z12**pw23-z1**pw23-z2**pw23)
         h120a = h120a*0.624d0*zbar**pw13*l0**pw23
         if (l12.le.5.d0) h120 = min(h120,h120a)
         if (l12.gt.5.d0) h120 = h120a
         h120 = min(h120,h120a)
      endif
      screen = exp(min(220.d0,h120))

      return
      end
