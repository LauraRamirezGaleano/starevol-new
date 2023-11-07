

************************************************************************

      SUBROUTINE kappa (tk,rok,muiinvk,x,ksh,kap1,kapr1,kapt1,kapx1
     &     ,kapy1,opamol,iter,iopa,iopaf,error)

************************************************************************
* Determine the total opacity coefficient and its derivatives          *
*     kap1 = kappa
*     kapr1 = d kappa/dro
*     kapt1 = d kappa/dT
*     kapx1 = d kappa/dX  (X = H mass fraction)
*     kapy1 = d kappa/dY  (Y = He mass fraction)
*     zkint : metallicity of the initial model
*     (must be specified in the parameter card)
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.diff'
c      include 'evolcom.kap'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.opa'

      integer mlth,nlth,itlth
      
      double precision rklth,tklth,opaclth
C     Warning: dimensions of opacity tables change with solar comp.
c..   Young+18 (TD: 30/01/2020)
c     common /lthopac/ opaclth(19,81,155),tklth(81),rklth(19),
c    &     mlth,nlth,itlth
       common /lthopac/ opaclth(19,85,155),tklth(85),rklth(19),
     &     mlth,nlth,itlth
c#ifdef GRID
c      parameter (mlth = 19,nlth = 63,itlth = 104)
c      common /lthopac/ rklth(19),tklth(63),opaclth(19,63,104)
c#elif AS09
c      parameter (mlth = 19,nlth = 85,itlth = 155)
c      common /lthopac/ rklth(19),tklth(85),opaclth(19,85,155)      
c#else
c      parameter (mlth = 19,nlth = 85,itlth = 120)
c      common /lthopac/ rklth(19),tklth(85),opaclth(19,85,120)
c#endif

      integer mt,nt,kdel,kxdel,kk(4)
      integer ir,it,j,iter
      integer ii,iim,iin
      integer error,opamol,ksh,iopa,iopaf

      double precision tk,rok,muiinvk
      double precision x(nsp)
      double precision kap1,kapr1,kapt1,kapx1,kapy1
      double precision dxy,kap1x,kap1y
      double precision ck(4),ckr(4),ckt(4)
      double precision tik,rotik,til,rotil,xx,xy,xc,xn,xo,xno,xsolc,
     &     xsoln,xsolo,xz,topa,uopa,d1opa,d2opa,xsi1,xsi2,xsi11,xsi12,
     &     yyopa1,yyopa2,yyopa3,yyopa4,kap4,kapr4,kapt4
      double precision yopa(4),y1opa(4),y2opa(4),y12opa(4)

      double precision opact,dopact,dopacr,dopactd,ftredge,fedge,fzedge
      double precision xspsol,zsol

      character*6 refsolar

      common /extopac/ opact,dopact,dopacr,dopactd,fedge,ftredge,fzedge
      common /solar/ xspsol(nis+6),zsol,refsolar

      mt = 0
      nt = 0
      dxy = 1.d-2
      tik = tk*1.d-3
      rotik = log10(rok/((tik*1.d-3)**3))
      rotik = max(rotik,-9.d0)
      rotik = min(rotik,1.d1)
      til = tk*1.d-6
      rotil = rok/til**3
      rotil = max(rotil,1.d-9)
      rotil = min(rotil,1.d10)
      xx = x(ih1)+x(ih2)
      xy = x(ihe3)+x(ihe4)
      xsolc = zkint*(xspsol(ic12)+xspsol(ic13)+xspsol(ic14))/zsol
      xsoln = zkint*(xspsol(in13)+xspsol(in14)+xspsol(in15))/zsol
      xsolo = zkint*(xspsol(io16)+xspsol(io17)+xspsol(io18))/zsol
      xc = x(ic12)+x(ic13)+x(ic14)-xsolc
      xn = x(in13)+x(in14)+x(in15)-xsoln
      xo = x(io16)+x(io17)+x(io18)-xsolo
      xc = max(0.d0,xc)
      xn = max(0.d0,xn)
      xo = max(0.d0,xo)
***   If solar calibration with diffusion, then metallicity should correspond to the actual
***   metallicity of the model
c      if (totm.eq.1.d0) then
c         xz =  1.d0-xx-xy
c      else
         xz = zkint
c      endif

*------------------------------------------------------------------------
***   calculation of the opacity coefficient for low-T regions and
***   variable mixture of C and O (including molecules of H20, OH, H2, CO
***   C2 CN and N2)
***   Only if the C/O ratio (in number) is greater than 0.8 (i.e. C/O in
***   mass greater than 0.6)
***   WARNING : below 1700K, dust not taken into account
*------------------------------------------------------------------------

      if (tk.le.5.d3.and.opamol.eq.1.and.x(ic12).gt.0.6d0*x(io16)) then
         if (iopa.eq.0) then 
            iopa = ksh
            if (iter.eq.1) write (nout,100) iopa
         endif
         if (tk.lt.1.7d3) write (nout,200)
         call opa_co (tk,rok,muiinvk,x,ksh,kap1,kapr1,kapt1,error)
         kapx1 = 0.d0
         kapy1 = 0.d0

         return

      endif
      
*------------------------------------------------------------------------
***   calculation of the opacity coefficient for H-rich and low-T regions
***   low temperature opacities : Ferguson et al. (2005)
*------------------------------------------------------------------------
      
      if (tk.le.tlhkap) then
         if (iopaf.eq.0) then 
            iopaf = ksh
            if (iter.eq.1.and.model.eq.modeli) write (nout,250) iopaf
         endif
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
         
         if (refsolar.eq.'GN93') then
            
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
               xsi1 = (xz-1.d-3)*1.d3
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
         else if
     $           (refsolar.eq.'AGS05'.or.refsolar.eq.'AGSS09'.or.
     $           refsolar.eq.'AY18') then
            if (xz.le.1.d-5) then
               kdel = 15
               xsi1 = xz*1.d5
            elseif (xz.le.3.d-5) then
               kdel = 14
               xsi1 = (xz-1.d-5)/3.d-5
            elseif (xz.le.1.d-4) then
               kdel = 13
               xsi1 = (xz-3.d-5)/7.d-5
            elseif (xz.le.3.d-4) then
               kdel = 12
               xsi1 = (xz-1.d-4)/3.d-4
            elseif (xz.le.1.d-3) then
               kdel = 11
               xsi1 = (xz-3.d-4)/7.d-4
            elseif (xz.le.2.d-3) then
               kdel = 10
               xsi1 = (xz-1.d-3)*1.d3
            elseif (xz.le.4.d-3) then
               kdel = 9
               xsi1 = (xz-2.d-3)*5.d2
            elseif (xz.le.1.d-2) then
               kdel = 8
               xsi1 = (xz-4.d-3)/6.d-3
            elseif (xz.le.2.d-2) then
               kdel = 7
               xsi1 = (xz-1.d-2)*1.d2
            elseif (xz.le.3.d-2) then
               kdel = 6
               xsi1 = (xz-2.d-2)*1.d2
            elseif (xz.le.4.d-2) then
               kdel = 5
               xsi1 = (xz-3.d-2)*1.d2
            elseif (xz.le.5.d-2) then
               kdel = 4
               xsi1 = (xz-4.d-2)*1.d2
            elseif (xz.le.6.d-2) then
               kdel = 3
               xsi1 = (xz-5.d-2)*5.d1
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
               kxdel = 16
               xsi2 = (xx-0.1d0)*1.d1
            elseif (xx.le.0.35d0) then
               kxdel = 32
               xsi2 = (xx-0.2d0)/1.5d-1
            elseif (xx.le.0.5d0) then
               kxdel = 48
               xsi2 = (xx-0.35d0)/1.5d-1
            elseif (xx.le.0.7d0) then
               kxdel = 64
               xsi2 = (xx-0.5d0)*5.d0
            elseif (xx.le.0.8d0) then
               kxdel = 80
               xsi2 = (xx-0.7d0)*1.d1
            elseif (xx.le.0.9d0) then
               kxdel = 96
               xsi2 = (xx-0.8d0)*1.d1
            elseif (xx.gt.0.9d0) then
               kxdel = 112
               xsi2 = 0.d0
            endif

            kk(1) = kxdel+kdel
            kk(2) = kk(1)+1
            kk(3) = kk(1)+16
            kk(4) = kk(3)+1
            
         elseif (refsolar.eq.'GRID') then
C Warning: dimension of opacitiy tables changes for GRID version 
C Indices for Asplund/Cunha solar comp.

            if (xz.le.1.d-4) then
               kdel = 12
               xsi1 = xz*1.d4
            elseif (xz.le.3.d-4) then
               kdel = 11
               xsi1 = (xz-1.d-4)/2.d-4
            elseif (xz.le.1.d-3) then
               kdel = 10
               xsi1 = (xz-3.d-4)/7.d-4
            elseif (xz.le.2.d-3) then
               kdel = 9
               xsi1 = (xz-1.d-3)*1.d3
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
               kxdel = 13
               xsi2 = (xx-0.1d0)*1.d1
            elseif (xx.le.0.35d0) then
               kxdel = 26
               xsi2 = (xx-0.2d0)/1.5d-1
            elseif (xx.le.0.5d0) then
               kxdel = 39
               xsi2 = (xx-0.35d0)/1.5d-1
            elseif (xx.le.0.7d0) then
               kxdel = 52
               xsi2 = (xx-0.5d0)*5.d0
            elseif (xx.le.0.8d0) then
               kxdel = 65
               xsi2 = (xx-0.7d0)*1.d1
            elseif (xx.le.0.9d0) then
               kxdel = 78
               xsi2 = (xx-0.8d0)*1.d1
            elseif (xx.gt.0.9d0) then
               kxdel = 91         
               xsi2 = 0.d0
            endif
            
            
            kk(1) = kxdel+kdel
            kk(2) = kk(1)+1
            kk(3) = kk(1)+13
            kk(4) = kk(3)+1
         endif

         do j = 1,4
            ck(j) = 0.d0
            ckr(j) = 0.d0
            ckt(j) = 0.d0
            if (kk(j).ne.0) then
               yyopa1 = opaclth(mt,nt,kk(j))
               yyopa2 = opaclth(mt+1,nt,kk(j))
               yyopa3 = opaclth(mt+1,nt+1,kk(j))
               yyopa4 = opaclth(mt,nt+1,kk(j))
               call interlin (topa,uopa,d1opa,d2opa,yyopa1,yyopa2,
     &              yyopa3,yyopa4,ck(j),ckr(j),ckt(j))
            endif
         enddo
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
            write (nout,300) ksh,xz,xx,xc,xo,tk,rok,rotik
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
      if (rotation) then
         xno = xo
         call opac (xz,xx+xx*dxy,xc,xno,til,rotil,error)
         if (error.gt.0) then
            write (nout,400) ksh,xz,xx+xx*dxy,xc,xno,tk,rok,rotik
            error = 81
            return
         endif
         kap1x = 10.d0**opact
         kapx1 = (kap1x-kap1)/(kap1*dxy)
         call opac (xz,xx-xx*dxy,xc,xno,til,rotil,error)
         if (error.gt.0) then
            write (nout,400) ksh,xz,xx-xx*dxy,xc,xno,tk,rok,rotik
            error = 81
            return
         endif
         kap1y = 10.d0**opact
         kapy1 = xy*(kap1y-kap1)/(kap1*xx*dxy)
      endif

 100  format (' molecular opacity computed from shell #',i4)
 200  format (2x,'WARNING : T < 1700K ! dust not taken into account ',
     &        'in opacity')
 250  format (' Fergusson opacity computed from shell #',i4)
 300  format (3x,'KAPPA : pb table 3 - O rich : #',i4,', xZ =',1pe10.3,
     &     ', xH = ',1pe9.3,', xC = ',1pe9.3,', xO = ',1pe9.3,', T = ',
     &     1pe9.3,', rho = ',1pe9.3,', R = ',0pf8.3)
 400  format (3x,'KAPPA rotation : pb table 3 - O rich : #',i4,', xZ =',
     &     1pe9.3,', xH = ',1pe9.3,', xC = ',1pe9.3,', xO = ',1pe9.3,
     &     ', rho = ',1pe9.3,', R = ',0pf8.3)

      return
      end
