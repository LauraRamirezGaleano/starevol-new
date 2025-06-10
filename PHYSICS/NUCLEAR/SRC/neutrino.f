
      PROGRAM NEUTRINO

      implicit none

      integer nrog,ntg
      integer qdom
      integer i,j

      double precision tgl,dtgl,rogl,roglsi,drogl
      double precision rok,romue,tk,t8
      double precision x12c,x16o,zzmean,xamean,amean
      double precision qplasm,qphot,qpair,qbrem,qtot

      open (unit = 10,file = 'neutrino.out')

      x12c = 0.20d0
      x16o = 1.d0-x12c
      zzmean = (x12c+x16o)/2.d0
      xamean = x12c/12.d0+x16o/16.d0
      zzmean = zzmean/xamean
      amean = 1.d0/xamean

      tgl = 7.4d0
      dtgl = 0.1d0
      ntg = 21
      drogl = 0.1d0
      nrog = 71

      do i = 1,ntg
         tgl = tgl+dtgl
         tk = 10.d0**tgl
         rogl = 2.9d0
         do j = 1,nrog
            rogl = rogl+drogl
            roglsi = rogl+6.d0
            rok = 10.d0**rogl
            romue = rok/2.0d0
            t8 = tk*1.d-8

*________________________
***   global contribution
*------------------------

            call nulos (romue,rok,t8,zzmean,amean,
     &                  qplasm,qphot,qpair,qbrem)

            qtot = (qplasm+qphot+qpair+qbrem)/rok

            qdom = 0
            if (qplasm.gt.qphot.and.qplasm.gt.qpair.and.
     &           qplasm.gt.qbrem) qdom = 1
            if (qphot.gt.qplasm.and.qphot.gt.qpair.and.
     &           qphot.gt.qbrem) qdom = 2
            if (qpair.gt.qplasm.and.qpair.gt.qphot.and.
     &           qpair.gt.qbrem) qdom = 3
            if (qbrem.gt.qplasm.and.qbrem.gt.qphot.and.
     &           qbrem.gt.qpair) qdom = 4

            write (10,1000) roglsi,tgl,qdom

         enddo
      enddo

      close (10)

      stop 'normal'

 1000 format (1x,0pf5.2,1x,0pf4.2,2x,i1)

      end



************************************************************************

      SUBROUTINE NULOS (romue,rok,t8,zzmean,amean,
     &                  qplasm,qphot,qpair,qbrem)

************************************************************************
* Calculate the production of neutrinos inside the plasma              *
* vnu = number of neutrinos other than electron neutrinos              *
* q = energy loss rate (ergs/s/gr)                                     *
************************************************************************

      implicit none

      integer i

      double precision romue,rok,t8,zzmean,amean
      double precision qplasm,qphot,qpair,qbrem
      double precision anu,bnu,cnu,dnu,cvnu,ca2nu,cvpnu,banu,bbnu
      double precision romue3,alam,alam2,alam3,alam52,tbr,tbr15,tbr2,
     &     tbr3,tbr5,vnu,cv2,cvp2,cap2,cva,cvaa,cpl,xi,xi2,xi3,xibr,
     &     xibr2,gl,fx,fa,f
      double precision sq1,sq2,sq3,sqpho,arg,qph1,qph2,sqp1,sqp2,sqp,
     &     qp1,qp2,fbr,gbr,zzbr

      common /nudata/ anu(3,3),bnu(3,3),cnu(3),dnu(5),cvnu,ca2nu,cvpnu,
     &     banu(6),bbnu(5)

      dimension f(3)

      romue3 = romue**3
      alam = t8/5.9302d1
      alam2 = alam*alam
      alam3 = alam2*alam
      alam52 = alam**(-2.5d0)
      tbr = t8
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
      xibr = rok/(7.05d6*tbr15+5.12d4*tbr3)
      xibr2 = xibr*xibr
      gl = dnu(1)+alam2*(dnu(2)+alam2*(dnu(3)+alam2*(dnu(4)+
     &     alam2*dnu(5))))
      bbnu(3) = 7.75d5*tbr15+247.d0*tbr**3.85d0
      bbnu(4) = 4.07d0+0.0240d0*tbr**1.40d0
      bbnu(5) = 4.59d-5*tbr**(-0.110d0)
      do i = 1,3
         fx = (anu(i,1)+anu(i,2)*xi+anu(i,3)*xi2)*exp(-cnu(i)*xi)
         fa = xi3+bnu(i,1)/alam+bnu(i,2)/alam2+bnu(i,3)/alam3
         f(i) = fx/fa
      enddo

*_________________________________________
***   contribution of the plasma neutrinos
*-----------------------------------------

      qplasm = cpl*romue3*f(1)

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

      fbr = 1.d0/(banu(1)+banu(2)/tbr2+banu(3)/tbr5)+(1.26d0*
     &     (1.d0+1.d0/xibr))/(1.d0+bbnu(1)/xibr+bbnu(2)/xibr2)
      gbr = 1.d0/((1.d0+rok*1.d-9)*(banu(4)+banu(5)/tbr2+banu(6)/
     &     tbr5))+1.d0/(bbnu(3)/rok+bbnu(4)+bbnu(5)*rok**0.656d0)
      zzbr = 0.5738d0*rok*zzmean*zzmean/amean*tbr5*tbr
      qbrem = zzbr*0.5d0*(cva*fbr-cvaa*gbr)

      return
      end



************************************************************************

      BLOCK DATA NUDAT

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
