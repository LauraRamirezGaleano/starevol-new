

************************************************************************

      SUBROUTINE yequil (v,y,ytol,dtn,taumix)

************************************************************************
*   Set to low abundant elements (y<ytol) to their equilibrium abundance
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.nuc'

      integer i,j
      integer ii,ij,ik,il,nn,mm

      double precision rdest,sprod,am,dtdest,yeqmin,dtn
      double precision v,y,yeq,ytol,taumix

      dimension v(nreac),y(nsp)
      dimension nn(nsp)


      yeqmin = min(1.d-14,ytol*1.d-2)
      do i = 1,nis-1
         if (y(i).lt.yeqmin) then
            sprod = 0.d0
            rdest = 0.d0
            do mm = 1,nreac
               ii = k1(mm)
               nn(ii) = k2(mm)
               ij = k3(mm)
               nn(ij) = k4(mm)
               ik = k5(mm)
               il = k7(mm)
               if (i.eq.ii.or.i.eq.ij) then
                  if (i.eq.ij) then
                     j = ii
                  else
                     j = ij
                  endif
                  if (nn(i).eq.1.and.nn(j).eq.1) then
                     rdest = rdest+v(mm)*y(j)
                  else
                     am = v(mm)*nn(i)/(fact(nn(i))*fact(nn(j)))
                     rdest = rdest+am*y(j)**nn(j)
                  endif
               elseif (i.eq.ik.or.i.eq.il) then
                  if (nn(ii).eq.1.and.nn(ij).eq.1) then
                     sprod = sprod+v(mm)*y(ii)*y(ij)
                  else
                     am = v(mm)*nn(ii)/(fact(nn(ii))*fact(nn(ij)))
                     sprod = sprod+am*y(ii)**nn(ii)*y(ij)**nn(ij)
                  endif
               endif
            enddo
            if (rdest.gt.0.d0) then
               dtdest = 1.d0/rdest
               yeq = sprod/rdest
c..   set y to yeq if dtn > 10 dtdest & tau_mix > dtdest & y < yeqmin
               if (dtn.ge.10.d0*dtdest.and.yeq.lt.yeqmin.and.
     &              taumix.gt.dtdest) then
                  y(i) = yeq
               endif
            endif
         endif
      enddo


      return
      end
