

************************************************************************

      SUBROUTINE neutri (tk,rok,mueinvk,x,ksh,q,dqt,dqro)

************************************************************************
* Calculate the neutrino loss rates due to the nuclear reactions       *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.eng'
      include 'evolcom.nuc'
      include 'evolcom.spec'

      integer i,j,i1,i2,i4
      integer ksh,iroe,ite

      double precision tk,rok,x,mueinvk
      double precision q,q2,qrho,qrho2,dqt,dqro
      double precision enu,enu2,v,tkp,dtk,rokp,drok
      double precision ee,ee2,eerho,eerho2,eecap,eecap2
      double precision a1,a2,b1,b2,y
      double precision vcapte,enuelec,rhoye,tcapte,qinuelec

      common /electcap/ rhoye(nroe),tcapte(nte),vcapte(nte,nroe,nrce),
     &     enuelec(nte,nroe,nrce)
      common /electcap_par/ ite,iroe,a1,a2,b1,b2

      dimension enu(nreac),enu2(nreac),qinuelec(nrce),v(nreac),x(nsp),
     &     y(nsp)

      ee = 0.d0
      ee2 = 0.d0
      eerho = 0.d0
      eerho2 = 0.d0
      eecap = 0.d0
      eecap2 = 0.d0
      do i = 1,nsp
         y(i) = x(i)/anuc(i)
      enddo
***   reactions independent of rho*Ye
      do j = 1,nre
         i2 = k2(j)
         i4 = k4(j)
         enu(j) = flu(ksh,j)*qinu(j)
         enu2(j) = flu2(j)*qinu(j)
         eerho = eerho+dble(i2+i4-1)*enu(j)
         ee = ee+enu(j)
         ee2 = ee2+enu2(j)
      enddo

***   rho*Ye-T dependent reactions (electron captures only)

***   rho,T+dT
      dtk = tk*1.d-4
      tkp = tk+dtk
      call vit (tkp,rok,mueinvk,ksh,v,2)
      do j = nre+1,nreac
         i = j-nre
         qinuelec(i) = a1*a2*enuelec(ite+1,iroe+1,i)+b1*b2*
     &        enuelec(ite,iroe,i)+a1*b2*enuelec(ite+1,iroe,i)+
     &        b1*a2*enuelec(ite,iroe+1,i)
      enddo
      do j = nre+1,nreac
         i1 = k1(j)
         i = j-nre
         flu2(j) = y(i1)*v(j)
         eecap2 = eecap2+flu2(j)*qinuelec(i)
      enddo

***   rho,T
      call vit (tk,rok,mueinvk,ksh,v,2)
      do j = nre+1,nreac
         i = j-nre
         qinuelec(i) = a1*a2*enuelec(ite+1,iroe+1,i)+b1*b2*
     &        enuelec(ite,iroe,i)+a1*b2*enuelec(ite+1,iroe,i)+
     &        b1*a2*enuelec(ite,iroe+1,i)
      enddo
      do j = nre+1,nreac
         i1 = k1(j)
         i = j-nre
         flu(ksh,j) = y(i1)*v(j)
         eecap = eecap+flu(ksh,j)*qinuelec(i)
      enddo

***   rho+drho,T
      drok = rok*1.d-4
      rokp = rok+drok
      call vit (tk,rokp,mueinvk,ksh,v,2)
      do j = nre+1,nreac
         i = j-nre
         qinuelec(i) = a1*a2*enuelec(ite+1,iroe+1,i)+b1*b2*
     &           enuelec(ite,iroe,i)+a1*b2*enuelec(ite+1,iroe,i)+
     &        b1*a2*enuelec(ite,iroe+1,i)
      enddo
      do j = nre+1,nreac
         i1 = k1(j)
         i = j-nre
         eerho2 = eerho2+y(i1)*v(j)*qinuelec(i)
      enddo
      ee = ee+eecap
      ee2 = ee2+eecap2
      q = ee*econv
      q2 = ee2*econv
      dqt = (q2-q)/dtk
      qrho = eecap*econv
      qrho2 = eerho2*econv
      dqro = eerho*econv/rok+(qrho2-qrho)/drok

      return
      end
