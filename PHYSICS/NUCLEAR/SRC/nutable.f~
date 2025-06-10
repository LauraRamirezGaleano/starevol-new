
      PROGRAM NUTABLE

************************************************************************
*     tables for interpolation in romue and t8                         *
*     of rates of neutrino energy losses                               *
************************************************************************

      implicit double precision (a-h,o-z)

      parameter(nromue = 105, ntnu = 44)

      common /nurates/ romuenu(nromue),tnu8(ntnu),rnu(nromue,ntnu)

      open (unit = 10,file = 'nusub.f')

      write (10,1000)
      write (10,1200)
      write (10,1000)
      write (10,1100)
      write (10,1300)
      write (10,1100)
      write (10,1000)
      write (10,1400)
      write (10,1500)

      t8l = -1.15d0
      dt8l = 0.05d0
      do j = 1,ntnu
         t8l = t8l+dt8l
         tnu8(j) = 10.d0**t8l
         if (j.eq.1) then
            romuel = -3.2d0
            dromuel = 0.1d0
         endif
         do i = 1,nromue
            if (j.eq.1) then
               romuel = romuel+dromuel
               romuenu(i) = 10.d0**romuel
            endif

            call nulos (romuenu(i),tnu8(j),qpl,qphot,qpair)

            rnu(i,j) = log(qpl+qphot+qpair)
            if (rnu(i,j).lt.-99.d0) rnu(i,j) = -99.9d0
         enddo
      enddo

      write (10,1501)
      do ii = 1,8
         iinf = 5*ii-4
         isup = iinf+4
         write (10,1502) (tnu8(j),j = iinf,isup)
      enddo
      write (10,1503) (tnu8(j),j = 41,ntnu)
      write (10,1504)
      do ii = 1,20
         iinf = 5*ii-4
         isup = iinf+4
         write (10,1502) (romuenu(i),i = iinf,isup)
      enddo
      write (10,1505) (romuenu(i),i = 101,nromue)

      do j = 1,ntnu
         write (10,1600) tnu8(j)
         write (10,1700) j,j
         write (10,1800) (rnu(i,j),i = 1,7)
         write (10,1800) (rnu(i,j),i = 8,14)
         write (10,1800) (rnu(i,j),i = 15,21)
         write (10,1800) (rnu(i,j),i = 22,28)
         write (10,1800) (rnu(i,j),i = 29,35)
         write (10,1800) (rnu(i,j),i = 36,42)
         write (10,1800) (rnu(i,j),i = 43,49)
         write (10,1800) (rnu(i,j),i = 50,56)
         write (10,1800) (rnu(i,j),i = 57,63)
         write (10,1800) (rnu(i,j),i = 64,70)
         write (10,1800) (rnu(i,j),i = 71,77)
         write (10,1800) (rnu(i,j),i = 78,84)
         write (10,1800) (rnu(i,j),i = 85,91)
         write (10,1800) (rnu(i,j),i = 92,98)
         write (10,1900) (rnu(i,j),i = 99,nromue)
      enddo

      write (10,1000)
      write (10,2000)
      write (10,1000)

      close (10)

 1000 format (1x)
 1100 format (72('*'))
 1200 format (6x,'BLOCK DATA NUSUB')
 1300 format ('*',14x,'tables of rates of neutrino energy losses',
     &     15x,'*')
 1400 format (6x,'implicit none',//,
     &     6x,'double precision romuenu,tnu8,rnu',/,
     &     6x,'integer i,j',/)
 1500 format (6x,'common /nurates/ romuenu(105),tnu8(44),rnu(105,44)')
 1501 format (1x,/,6x,'data (tnu8(j),j = 1,44) /')
 1502 format (5x,'& ',5(1pd12.6,','))
 1503 format (5x,'& ',3(1pd12.6,','),1pd12.6,'/')
 1504 format (1x,/,6x,'data (romuenu(i),i = 1,105) /')
 1505 format (5x,'& ',4(1pd12.6,','),1pd12.6,'/')
 1600 format (1x,/,'***** rates at T8 = ',0pf8.5,' K')
 1700 format (6x,'data ((rnu(i,j),i = 1,105),j = ',i2,',',i2,') /')
 1800 format (5x,'& ',7(f8.4,','))
 1900 format (5x,'& ',6(f8.4,','),f8.4,'/')
 2000 format (6x,'end',/)

      end



************************************************************************

      SUBROUTINE NULOS (romue,tt8,qpl,qphot,qpair)

************************************************************************
* Calculate the production of neutrinos inside the plasma              *
* vnu = number of neutrinos other than electron neutrinos              *
* q = energy loss rate (ergs/s/gr)                                     *
************************************************************************

      implicit none

      integer i

      double precision tt8
      double precision anu,bnu,cnu,dnu,cvnu,ca2nu,cvpnu,banu,bbnu
      double precision f
      double precision vnu,t59,cv2,cvp2,cap2,cva,cvaa,cpl,romue,romue3,
     &     alam,alam2,alam3,alam52,xi,xi2,xi3,gl,fx,fa,tbr,tbr15,tbr2,
     &     tbr3,tbr5,qpl,sq1,sq2,sq3,sqpho,arg,qph1,qph2,qphot,
     &     sqp1,sqp2,sqp,qp1,qp2,qpair

      common /nudata/ anu(3,3),bnu(3,3),cnu(3),dnu(5),cvnu,ca2nu,cvpnu,
     &     banu(6),bbnu(5)

      dimension f(3)

      romue3 = romue**3
      alam = tt8/5.9302d1
      alam2 = alam*alam
      alam3 = alam2*alam
      alam52 = alam**(-2.5d0)
      tbr = tt8
      tbr15 = tbr**1.5d0
      tbr2 = tbr*tbr
      tbr3 = tbr2*tbr
      tbr5 = tbr2*tbr3

      vnu = 2.d0
      cv2 = cvnu*cvnu
      cvp2 = cvpnu*cvpnu
      cap2 = ca2nu
      cva = cv2+ca2nu+vnu*(cvp2+cap2)
      cvaa = cv2-ca2nu+vnu*(cvp2-cap2)
      cpl = cv2+vnu*cvp2

      xi = (romue*1.d-9)**(1.d0/3.d0)/alam
      xi2 = xi*xi
      xi3 = xi2*xi
      gl = dnu(1)+alam2*(dnu(2)+alam2*(dnu(3)+alam2*(dnu(4)+
     &     alam2*dnu(5))))
      do i = 1,3
         fx = (anu(i,1)+anu(i,2)*xi+anu(i,3)*xi2)*exp(-cnu(i)*xi)
         fa = xi3+bnu(i,1)/alam+bnu(i,2)/alam2+bnu(i,3)/alam3
         f(i) = fx/fa
      enddo

      bbnu(3) = 7.75d5*tbr15+247.d0*tbr**3.85d0
      bbnu(4) = 4.07d0+0.0240d0*tbr**1.40d0
      bbnu(5) = 4.59d-5*tbr**(-0.110d0)

*_________________________________________
***   contribution of the plasma neutrinos
*-----------------------------------------

      qpl = cpl*romue3*f(1)

*________________________________________
***   contribution of the photo neutrinos
*----------------------------------------

      sq1 = 0.666d0*(1.d0+2.045d0*alam)**(-2.066d0)
      sq2 = alam*1.875d8+alam2*1.d8*(1.653d0+8.449d0*alam-1.604d0*alam2)
      sq3 = 1.d0+romue/sq2
      sqpho = sq1/sq3
      arg = (0.556d0*xi**(4.48d0))/(150.d0+xi**(3.30d0))
      qph1 = 0.5d0*cva*(1.d0-cvaa*sqpho/cva)
      qph2 = 0.893d9*alam**8*(1.d0+143.8d0*alam**(3.555d0))**(0.3516d0)*
     &     xi2*xi*exp(arg)
      qphot = qph1*qph2*f(2)

*_______________________________________
***   contribution of the pair neutrinos
*---------------------------------------

      sqp1 = 1.d0/(10.748d0*alam2+0.3967d0*sqrt(alam)+1.005d0)
      sqp2 = 1.d0+romue/(7.692d7*alam3+9.715d6*sqrt(alam))
      sqp2 = sqp2**(-0.30d0)
      sqp = sqp1*sqp2
      qp1 = 0.5d0*cva*gl*exp(-2.d0/alam)
      qp2 = 1.d0+cvaa*sqp/cva
      qpair = qp1*qp2*f(3)

*_________________________________________________
***   contribution of the bremstrhahlung neutrinos
*-------------------------------------------------

***   cannot be tabulated

      return
      end



************************************************************************

      BLOCK DATA EVODAT

************************************************************************

      implicit none

      integer i

      double precision anu,bnu,cnu,dnu,cvnu,ca2nu,cvpnu,banu,bbnu

      common /nudata/ anu(3,3),bnu(3,3),cnu(3),dnu(5),cvnu,ca2nu,cvpnu,
     &     banu(6),bbnu(5)

      data (anu(1,i),i = 1,3)/2.32d-7,8.449d-8,1.787d-8/
      data (anu(2,i),i = 1,3)/4.886d10,7.58d10,6.023d10/
      data (anu(3,i),i = 1,3)/6.002d19,2.084d20,1.872d21/
      data (bnu(1,i),i = 1,3)/2.581d-2,1.734d-2,6.99d-4/
      data (bnu(2,i),i = 1,3)/6.29d-3,7.483d-3,3.061d-4/
      data (bnu(3,i),i = 1,3)/9.383d-1,-4.141d-1,5.829d-2/
      data (cnu(i),i = 1,3)/0.56457d0,1.5654d0,5.5924d0/
      data (dnu(i),i = 1,5)/1.d0,-13.04d0,133.5d0,1534.0d0,918.6d0/
      data cvnu,ca2nu,cvpnu/0.9638d0,0.25d0,0.0362d0/
      data (banu(i),i = 1,6)/23.5,6.83d4,7.81d8,230.,6.70d5,7.66d9/
      data (bbnu(i),i = 1,5)/1.47,0.0329,1.0,1.0,1.0/

      end

