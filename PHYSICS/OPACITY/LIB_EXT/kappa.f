

************************************************************************

      SUBROUTINE kappa (kap1,kapr1,kapt1,kapx1,kapy1,error)

************************************************************************
* Determine the total opacity coefficient and its derivatives          *
*     kap1 = kappa
*     kapr1 = d kappa/dro
*     kapt1 = d kappa/dT
*     kapx1 = d kappa/dX  (X = H mass fraction)
*     kapy1 = d kappa/dY  (Y = He mass fraction)
*     zkint : metallicity of the initial model
*     (must be specified in the parameter card)
************************************************************************

      implicit none

      include 'evolpar.chem'
      include 'evolpar.conv'
      include 'evolpar.kap'
      include 'evolpar.str'

      include 'evolcom.diff'
      include 'evolcom.kap'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.opa'
      include 'evolcom.shk'

      integer mt,nt,kdel,kxdel,kk
      integer ir,it,j
c      integer ii,iim,iin
      integer error

      double precision kap1,kapr1,kapt1,kapx1,kapy1
      double precision dxy,kap1x,kap1y
      double precision ck,ckr,ckt
c      double precision yopa,y1opa,y2opa,y12opa
      double precision tik,rotik,til,rotil,xx,xy,xc,xn,xo,xno,xsolc,
     &     xsoln,xsolo,xz,topa,uopa,d1opa,d2opa,xsi1,xsi2,xsi11,xsi12,
     &     yyopa1,yyopa2,yyopa3,yyopa4,kap4,kapr4,kapt4

      double precision opact,dopact,dopacr,dopactd,ftredge,fedge,fzedge

      dimension ck(4),ckr(4),ckt(4),kk(4)
c      dimension yopa(4),y1opa(4),y2opa(4),y12opa(4)

      common /extopac/ opact,dopact,dopacr,dopactd,fedge,ftredge,fzedge

      mt = 0
      nt = 0
      dxy = 1.d-2
      tik = tk*1.d-3
      rotik = log10(rok/((tik*1.d-3)**3))
c      rotik = max(rotik,-9.d0)
      rotik = min(rotik,1.d1)
      til = tk*1.d-6
      rotil = rok/til**3
c      rotil = max(rotil,1.d-9)
      rotil = min(rotil,1.d10)
      xx = x(ih1)+x(ih2)
      xy = x(ihe3)+x(ihe4)
      xsolc = zkint*0.173285d0
      xsoln = zkint*0.053152d0
      xsolo = zkint*0.482272d0
      xc = x(ic12)+x(ic13)+x(ic14)-xsolc
      xn = x(in13)+x(in14)+x(in15)-xsoln
      xo = x(io16)+x(io17)+x(io18)-xsolo
      xc = max(0.d0,xc)
      xn = max(0.d0,xn)
      xo = max(0.d0,xo)
c     xz = zkint
c      if ((xx.ge.1.d-6.and.xc.lt.2.d0*xsolc.and.xn.lt.2.d0*xsoln.and.
c     &     xo.lt.2.d0*xsolo).or.nphase.le.5) then
c     &     xo.lt.2.d0*xsolo).or.nphase.le.3) then
c         xz = 1.d0-xx-xy
c      else
         xz = zkint
c      endif


*------------------------------------------------------------------------
***   calculation of the opacity coefficient for H-rich and low-T regions
***   low temperature opacities : Ferguson et al. (2005)
*------------------------------------------------------------------------

      if (tk.le.tlhkap) then
         do ir = 1,mlth
            if (rklth(ir).gt.rotik) goto 10
         enddo
 10      mt = ir-1
         do it = 1,nlth
            if (tklth(it).le.tk) goto 20
         enddo
 20      nt = it-1
         if (mt.lt.1) mt = 1
         if (nt.lt.1) nt = 1
         if (mt.gt.(mlth-1)) mt = mlth-1
         if (nt.gt.(nlth-1)) nt = nlth-1
         topa = (rotik-rklth(mt))/(rklth(mt+1)-rklth(mt))
         uopa = (tk-tklth(nt))/(tklth(nt+1)-tklth(nt))
         d1opa = rklth(mt+1)-rklth(mt)
         d2opa = tklth(nt+1)-tklth(nt)

         if (xz.le.1.d-5) then
            kdel = 14
            xsi1 = xz/1.d-5
         elseif (xz.le.3.d-5) then
            kdel = 13
            xsi1 = (xz-1.d-5)/3.d-5
         elseif (xz.le.1.d-4) then
            kdel = 12
            xsi1 = (xz-3.d-5)/7.d-5
         elseif (xz.le.3.d-4) then
            kdel = 11
            xsi1 = (xz-1.d-4)/3.d-4
         elseif (xz.le.1.d-3) then
            kdel = 10
            xsi1 = (xz-3.d-4)/7.d-4
         elseif (xz.le.2.d-3) then
            kdel = 9
            xsi1 = (xz-1.d-3)/1.d-3
         elseif (xz.le.4.d-3) then
            kdel = 8
            xsi1 = (xz-2.d-3)*5.d2
         elseif (xz.le.1.d-2) then
            kdel = 7
            xsi1 = (xz-4.d-3)/6.d-3
         elseif (xz.le.2.d-2) then
            kdel = 6
            xsi1 = (xz-1.d-2)*1.d2
         elseif (xz.le.3.d-2) then
            kdel = 5
            xsi1 = (xz-2.d-2)*1.d2
         elseif (xz.le.4.d-2) then
            kdel = 4
            xsi1 = (xz-3.d-2)*1.d2
         elseif (xz.le.6.d-2) then
            kdel = 3
            xsi1 = (xz-4.d-2)*5.d1
         elseif (xz.le.8.d-2) then
            kdel = 2
            xsi1 = (xz-6.d-2)*5.d1
         elseif (xz.le.1.d-1) then
            kdel = 1
            xsi1 = (xz-8.d-2)*5.d1
         elseif (xz.gt.1.d-1) then
            kdel = 1
            xsi1 = 1.d0
         endif
            
         if (xx.le.0.1d0) then
            kxdel = 0
            xsi2 = 1.d0
         elseif (xx.le.0.2d0) then
            kxdel = 15
            xsi2 = (xx-0.1d0)*1.d1
         elseif (xx.le.0.35d0) then
            kxdel = 30
            xsi2 = (xx-0.2d0)/1.5d-1
         elseif (xx.le.0.5d0) then
            kxdel = 45
            xsi2 = (xx-0.35d0)/1.5d-1
         elseif (xx.le.0.7d0) then
            kxdel = 60
            xsi2 = (xx-0.5d0)*5.d0
         elseif (xx.le.0.8d0) then
            kxdel = 75
            xsi2 = (xx-0.7d0)*1.d1
         elseif (xx.le.0.9d0) then
            kxdel = 90
            xsi2 = (xx-0.8d0)*1.d1
         elseif (xx.gt.0.9d0) then
            kxdel = 105
            xsi2 = 0.d0
         endif

         kk(1) = kxdel+kdel
         kk(2) = kk(1)+1
         kk(3) = kk(1)+15
         kk(4) = kk(3)+1

         do 70 j = 1,4
            ck(j) = 0.d0
            ckr(j) = 0.d0
            ckt(j) = 0.d0
            if (kk(j).eq.0) goto 70
            yyopa1 = opaclth(mt,nt,kk(j))
            yyopa2 = opaclth(mt+1,nt,kk(j))
            yyopa3 = opaclth(mt+1,nt+1,kk(j))
            yyopa4 = opaclth(mt,nt+1,kk(j))
            call interlin (topa,uopa,d1opa,d2opa,yyopa1,yyopa2,
     &           yyopa3,yyopa4,ck(j),ckr(j),ckt(j))
c            goto 70
***   cubic interpolation deactivated
c            do ii = 1,4
c               if (ii.eq.1) then
c                  iim = mt
c                  iin = nt
c               endif
c               if (ii.eq.2) iim = mt+1
c               if (ii.eq.3) iin = nt+1
c               if (ii.eq.4) iim = mt
c               yopa(ii) = opaclth(iim,iin,kk(j))
c               y1opa(ii) = (opaclth(iim+1,iin,kk(j))-opaclth(iim-1,iin,
c     &              kk(j)))/(rklth(iim+1)-rklth(iim-1))
c               y2opa(ii) = (opaclth(iim,iin+1,kk(j))-opaclth(iim,iin-1,
c     &              kk(j)))/(tklth(iin+1)-tklth(iin-1))
c               y12opa(ii) = (opaclth(iim+1,iin+1,kk(j))-opaclth(iim+1,
c     &              iin-1,kk(j))-opaclth(iim-1,iin+1,kk(j))+
c     &              opaclth(iim-1,iin-1,kk(j)))/((rklth(iim+1)-
c     &              rklth(iim-1))*(tklth(iin+1)-tklth(iin-1)))
c            enddo
c            call intercub (yopa,y1opa,y2opa,y12opa,d1opa,d2opa,topa,
c     &           uopa,ck(j),ckr(j),ckt(j))
 70      continue
         xsi12 = 1.d0-xsi2
         xsi11 = 1.d0-xsi1
         kap4 = xsi12*xsi11*ck(2)+xsi2*xsi11*ck(4)+xsi1*xsi2*ck(3)+
     &        xsi12*xsi1*ck(1)
         kapr4 = xsi12*xsi11*ckr(2)+xsi2*xsi11*ckr(4)+xsi1*xsi2*ckr(3)+
     &        xsi12*xsi1*ckr(1)
         kapt4 = xsi12*xsi11*ckt(2)+xsi2*xsi11*ckt(4)+xsi1*xsi2*ckt(3)+
     &        xsi12*xsi1*ckt(1)
         kap1 = 10.d0**kap4
         kapr1 = kap1*kapr4/rok
         kapt1 = kap1*(-3.d0*kapr4/tk+kapt4*2.302585093d0)
         kapx1 = 0.d0
         kapy1 = 0.d0

         return

*-------------------------------------------------
*** OPAL opacities for different Z and CO mixtures
*-------------------------------------------------

      else

         call opac (xz,xx,xc,xo,til,rotil,error)

         if (opact.gt.1.d15.or.error.gt.0) then
            write (*,1002) ksh,xz,xx,xc,xo,opact
            error = 81
            return
         endif

         kap1 = 10.d0**opact
         kapr1 = kap1*dopacr/rok
         kapt1 = kap1*(-3.d0*dopacr+dopact)/tk
         kapx1 = 0.d0
         kapy1 = 0.d0

      endif

      kapx1 = 0.d0
      kapy1 = 0.d0
      if (mixing) then
         xno = xo
         call opac (xz,xx+xx*dxy,xc,xno,til,rotil,error)
         if (error.gt.0) then
            write (*,1004) ksh,xz,xx+xx*dxy,xc,xno,opact
            error = 81
            return
         endif
         kap1x = 10.d0**opact
         kapx1 = (kap1x-kap1)/(kap1*dxy)
         call opac (xz,xx-xx*dxy,xc,xno,til,rotil,error)
         if (error.gt.0) then
            write (*,1005) ksh,xz,xx-xx*dxy,xc,xno,opact
            error = 81
            return
         endif
         kap1y = 10.d0**opact
         kapy1 = xy*(kap1y-kap1)/(kap1*xx*dxy)
      endif

 1002 format (' error itab 3 : shell ',i4,', xZ =',1pe10.3,', xH = ',
     &     1pe10.3,', xC = ',1pe10.3,', xO = ',1pe10.3,', kap = ',
     &     1pe10.3)
 1004 format (' error itab 3 diffusion : shell ',i4,', xZ =',1pe10.3,
     &     ', xH = ',1pe10.3,', xC = ',1pe10.3,', xO = ',1pe10.3,
     &     ', kap = ',1pe10.3)
 1005 format (' error itab 3 diffusion : shell ',i4,', xZ =',1pe10.3,
     &     ', xH = ',1pe10.3,', xC = ',1pe10.3,', xO = ',1pe10.3,
     &     ', kap = ',1pe10.3)

      return
      end
