

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
*                                                                      *
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
