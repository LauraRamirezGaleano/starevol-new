

************************************************************************

      SUBROUTINE interlin (t,u,d1,d2,y1,y2,y3,y4,func,funcr,funct)

************************************************************************
* Make a bilinear interpolation (+ derivatives)                        *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      double precision t,u,d1,d2,y1,y2,y3,y4,func,funcr,funct
      double precision t1,u1

      t1 = 1.d0-t
      u1 = 1.d0-u
      func = t1*u1*y1+t*u1*y2+t*u*y3+t1*u*y4
      funcr = (-u1*y1+u1*y2+u*y3-u*y4)/d1
      funct = (-t1*y1-t*y2+t*y3+t1*y4)/d2

      return
      end
