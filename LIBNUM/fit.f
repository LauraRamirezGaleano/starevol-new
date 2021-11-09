

************************************************************************

      SUBROUTINE fit (x,y,ndata,a,b)

************************************************************************
* Fit a line (y = a+b*x) on a set of given points                      *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      integer ndata
      integer i

      double precision x,y,a,b
      double precision sx,sy,st2,sxoss,t

      dimension x(ndata),y(ndata)

      sx = 0.d0
      sy = 0.d0
      st2 = 0.d0
      b = 0.d0
      do i = 1,ndata
         sx = sx+x(i)
         sy = sy+y(i)
      enddo
      sxoss = sx/dble(ndata)
      do i = 1,ndata
         t = x(i)-sxoss
         st2 = st2+t*t
         b = b+t*y(i)
      enddo
      b = b/st2
      a = (sy-sx*b)/dble(ndata)

      return
      end
