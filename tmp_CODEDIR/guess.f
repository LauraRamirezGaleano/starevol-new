

************************************************************************

      SUBROUTINE guess

************************************************************************
*  estimate the new stellar model at the first iteration               *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.mod'
      include 'evolcom.rot'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'

      double precision xmod,xa,eps
      double precision t,vt,r,vr

      common /varia/ xmod(nsh,neq),xa(nsh,neq),eps(neq)
      common /vartvt/ t(nsh),vt(nsh)
      common /varrvr/ r(nsh),vr(nsh)

      xa(1:nmod,1:neq) = xmod(1:nmod,1:neq)

c..   guess new structure
      xmod(1:nmod,1:neq) = xmod(1:nmod,1:neq)+(xmod(1:nmod,1:neq)-
     &     xa(1:nmod,1:neq))*phi

      t(1:nmod) = exp(xmod(1:nmod,4))
      vt(1:nmod) = exp(xa(1:nmod,4))
      r(1:nmod) = exp(xmod(1:nmod,2))
      vr(1:nmod) = exp(xa(1:nmod,2))
      vro(1:nmod) = ro(1:nmod)
      vp(1:nmod) = p(1:nmod)
      ve(1:nmod) = e(1:nmod)
      vs(1:nmod) = s(1:nmod)
      vxsp(1:nmod,1:nsp) = xsp(1:nmod,1:nsp)
      omega(1:nmod) = vomega(1:nmod)
      if (dmaccr.gt.0.d0) then
         eacc(1:nmod) = lacc*facc(1:nmod)/(1.d0+facc(1:nmod))
         eaccdr(1:nmod) = lacc*dfaccdr(1:nmod)/(1.d0+facc(1:nmod))**2
      else
         eacc(1:nmod) = 0.d0
         eaccdr(1:nmod) = 0.d0
      endif
      nnucacc(1:nmaxconv) = 0
      iprint(1:nmaxconv) = 0
      rnucacc = 0.d0

      return
      end
