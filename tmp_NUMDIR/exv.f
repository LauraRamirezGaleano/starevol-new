

************************************************************************

      DOUBLE PRECISION FUNCTION exv (x)

************************************************************************
* Truncate the exponential function                                    *
* Input : exponential factor                                           *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      double precision x

      if (x.gt.150.d0) goto 10
      if (x.lt.-150.d0) goto 20
      if (abs(x).lt.1.d-9) goto 30
      goto 40
   10 exv = (x-149.d0)*exp(150.d0)
      return
   20 exv = -exp(-150.d0)/(x+149.d0)
      return
   30 exv = 1.d0-x
      return
   40 exv = exp(x)

      return
      end
