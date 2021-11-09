

************************************************************************

      SUBROUTINE shearnwt (paru,parm,gammash0,gammash)

************************************************************************
* Search of the solution of the third-order polynomial equation        *
* for coupled shear mixing                                             *
* with thermal motions and composition gradient                        *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      integer iter

      double precision paru,parm,gammash0,gammash
      double precision eps
      double precision par1,par2,par3
      double precision gam,gamold,delgam
      double precision f,f1

      parameter (eps = 1.d-12)

      par1 = 12.d0*parm
      par2 = 2.d0*(1.d0+parm-3.d0*paru)
      par3 = 6.d0+2.d0*parm-paru

      gamold = gammash0
      gam = gamold

      iter = 1
 10   f = par1*gam**3+par2*gam*gam+par3*gam-paru
      f1 = 3.d0*par1*gam*gam+2.d0*par2*gam+par3
      gamold = gam
      gam = gamold-f/f1
      if (gam.lt.0.d0) gam = abs(gam)
      delgam = abs(gam-gamold)/abs(gamold)
      if (delgam.lt.eps) goto 20
      iter = iter+1
      goto 10

 20   gammash = gam

      return
      end
