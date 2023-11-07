

************************************************************************

      SUBROUTINE nuceng (nflu,y,v,ksh,ee,eerho,eerhocap)

************************************************************************
* Calculate the nuclear energy production rate                         *
* nflu =< 0 : rhot,T
* nflu =  1 : rhot,T+dT
* nflu =  2 : rhot+drho,T
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
      include 'evolcom.spec'

      integer nflu,ksh,ite,iroe
      integer i1,i2,i3,i4
      integer i,j,mm

      double precision ee,eerho,eerhocap,efn,flux,qinucap
      double precision v(nreac),vi(nreac),y(nsp)
      double precision a1,a2,b1,b2
      double precision vcapte,enuelec,rhoye,tcapte

      common /electcap/ rhoye(nroe),tcapte(nte),vcapte(nte,nroe,nrce),
     &     enuelec(nte,nroe,nrce)
      common /electcap_par/ ite,iroe,a1,a2,b1,b2

      common /engen/ efn(nreac),flux(nreac)

      ee = 0.d0
      eerho = 0.d0
      eerhocap = 0.d0
      do i = 1,nreac
         vi(i) = v(i)
      enddo

***   T dependence
*  reaction : k2  X(k1) + k4 X(k3) --> k8 X(k7) + k6 X(k5)
*        ex : 1    C12  + 1   He   --> 0  gamma + 1   O16

      if (nflu.lt.2) then
         do mm = 1,nreac
            i1 = k1(mm)
            i2 = k2(mm)
            i3 = k3(mm)
            i4 = k4(mm)
            if (i2.eq.1) then
               if (i4.eq.1) then
                  flux(mm) = vi(mm)*y(i1)*y(i3)
               else
                  flux(mm) = vi(mm)*y(i1)
               endif
            else
               flux(mm) = vi(mm)*y(i1)**i2/fact(i2)
            endif
            efn(mm) = flux(mm)*(qi(mm)-qinu(mm))
            ee = ee+efn(mm)
            if (nflu.le.0) then
               if (mm.le.nre) then
                  eerho = eerho+dble(i2+i4-1)*efn(mm)
               else
                  eerhocap = eerhocap+efn(mm)
               endif
            endif
         enddo
***   rho dependence
      else
         do j = nre+1,nreac
            i = j-nre
            qinucap = a1*a2*enuelec(ite+1,iroe+1,i)+b1*b2*
     &           enuelec(ite,iroe,i)+a1*b2*enuelec(ite+1,iroe,i)+
     &           b1*a2*enuelec(ite,iroe+1,i)
            i1 = k1(j)
            eerhocap = eerhocap+vi(j)*y(i1)*(qi(j)-qinucap)
         enddo
      endif
      ee = ee*econv
      eerho = eerho*econv
      eerhocap = eerhocap*econv

      return
      end
