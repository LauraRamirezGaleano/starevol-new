

************************************************************************

      SUBROUTINE tstepx (mtini,totm,nmod1,dtnx,idtnx,iburning)

************************************************************************
* Constraint the time-step by the abundances variations                *
* in the last two models                                               *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.diff'
      include 'evolcom.eos'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.conv'
      include 'evolcom.spec'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer nmod1,idtnx,iburning,inuctop
      integer ite,iroe,nphase0
      integer itmax,ishockb,ishockt
      integer i,k

      double precision tk,rok,mueinvk
      double precision mtini,totm,dtnx
      double precision ftsh0,ftshe0
      double precision v(nreac)
      double precision rcno
      double precision a1,a2,b1,b2
      double precision t9,t913,t923,vd,vdmean,dtnd,yp,yhe4,yc12,yne20,
     &     yn14,yo16,vpp,dtnH,vn14,v3a,vc12,vo16,vcc,vcca,vccp,voca,
     &     vocp,vocn,vne,vneinv,vHe4,dtHe4,dtcc,dtne,dtoo,dtsi
      double precision TburnC,TburnNe
      double precision Lmax,tmax,Mackmax,enucmax


      common /burning/TburnC,TburnNe
      common /electcap_par/ ite,iroe,a1,a2,b1,b2
      common /hydrodyn/ Lmax,tmax,Mackmax,enucmax,itmax,ishockb,ishockt

      nphase0 = nphase
      iburning = 1
      idtnx = 1
      dtnx = 1.d90
      dtnH = 1.d99
      dtHe4 = 1.d99
      dtcc = 1.d99
      dtne = 1.d99
      dtoo = 1.d99
      dtsi = 1.d99
      if (nphase.eq.2) then
         ftsh0 = min(ftsh,2.d-2)
      else
         ftsh0 = ftsh
      endif
      if (nphase.eq.4.or.thermalpulse) then
c         ftshe0 = min(ftshe,1.d-2*abs(1.d0-4.d0*
c     &        log10(0.02d0/(zkint+1.d-8))))
         ftshe0 = min(ftshe,5.d-2)
      else
         ftshe0 = ftshe
      endif

*_____________________________________________
***   determination of nuclearly active shells
*---------------------------------------------

      inuctop = nmod1
      do i = 2,nmod1
         if (t(i).lt.tnucmin) exit
      enddo
      inuctop = i-1

*________________________________________________________
***   time-step constrained by D-burning (PMS phase only)
*--------------------------------------------------------

      if (nphase.eq.1.and.xsp(1,ih2).gt.1.d-8.and.crz(1).lt.-1.and.
     &     iaccr.eq.0.and.totm.lt.3.d0.and.inuctop.gt.1) then
         vdmean = 0.d0
         do k = 1,inuctop
            t9 = t(k)*1.d-9
            t913 = t9**pw13
            t923 = t913*t913
            vd = (2.24d3/t923*exp(max(-220.d0,-3.72d0/t913))*(1.d0+
     &           0.112d0*t913+3.38d0*t923+2.65d0*t9))*ro(k)
            vdmean = vdmean+vd*dm(k)
         enddo
         vdmean = vdmean/m(inuctop)
         dtnd = 2.d0/(vdmean*xsp(1,ih1))
C         dtnd = ftsh/(vdmean*xsp(1,ih1))
         dtnd = max(dtnd,1.d3*sec/mtini,1.d3)
C         dtnd = max(dtnd,5.d2*sec)
         if (dtnd.le.dtnx) then
            dtnx = dtnd
            iburning = 1
         endif
      endif

*_____________________________
***   variable initializations
*-----------------------------

      do 20 k = 1,inuctop
         tk = t(k)
         rok = ro(k)
         mueinvk = 2.d0/(1.d0+xsp(inuctop,ih1))
         signt(k) = 10.4d0
         call vit (tk,rok,mueinvk,k,v,0)
         yp = xsp(k,ih1)
         yhe4 = xsp(k,ihe4)/anuc(ihe4)
         yc12 = xsp(k,ic12)/anuc(ic12)
         yn14 = xsp(k,in14)/anuc(in14)
         yo16 = xsp(k,io16)/anuc(io16)
         yne20 = xsp(k,ine20)/anuc(ine20)

*_______________________________________
***   time-step constrained by H-burning
*---------------------------------------

         if (tk.lt.5.d7.and.yhe4.gt.2.5d-4) then
***   2 proto ( 0 gamma, 0 nutri)  1 deutr
            vpp = v(ippd)*yp*0.5d0
            dtnH = ftsh0/vpp
            if (dtnH.lt.dtnx) then
               dtnx = dtnH
               idtnx = k
               iburning = 2
            endif
            goto 20
         endif
         if (tk.lt.1.5d8.and.yp.gt.1.d-3) then
***   1 C  12 ( 1 proto, 0 gamma)  1 N  13
            vc12 = v(icpg)*yc12
            dtnH = ftsh0/vc12
            if (dtnH.lt.dtnx) then
               dtnx = dtnH
               idtnx = k
               iburning = 3
            endif
***   1 N  14 ( 1 proto, 0 gamma)  1 O  15
            vn14 = v(inpg)*yn14
            dtnH = ftsh0/vn14
            if (dtnH.lt.dtnx) then
               dtnx = dtnH
               idtnx = k
               iburning = 4
            endif
            goto 20
         endif

*________________________________________
***   time-step constrained by He-burning
*----------------------------------------

         if (tk.gt.9.d7.and.tk.lt.TburnC.and.yhe4.gt.2.5d-4.and.
     &        yp.lt.1.d-10) then
***   3 alpha ( 0 gamma, 0 gamma)  1 C  12
            v3a = v(i3a)
***   1 C  12 ( 1 alpha, 0 gamma)  1 O  16
            vc12 = v(icag)
***   1 O  16 ( 1 alpha, 0 gamma)  1 Ne 20
            vo16 = v(ioag)
            vHe4 = (v3a*yhe4**2/6.d0+vc12*yc12+vo16*yo16)
            dtHe4 = ftshe0/vHe4
            if (dtHe4.le.dtnx) then
               dtnx = dtHe4
               idtnx = k
               iburning = 5
            endif
            goto 20
         endif

*_______________________________________
***   time-step constrained by C-burning
*---------------------------------------

         if (tk.ge.TburnC.and.yp.lt.1.d-10.and.yhe4.lt.2.5d-4) then
***   2 C  12 ( 0 gamma, 1 alpha)  1 Ne 20
            vcca = v(iccga)
***   2 C  12 ( 0 gamma, 1 proto)  1 Na 23
            vccp = v(iccgp)
***   1 O  16 ( 1 C  12, 1 NEUT )  1 AL 27
            vocn = v(iocn)
***   1 O  16 ( 1 C  12, 1 PROT )  1 AL 27
            vocp = v(iocp)
***   1 O  16 ( 1 C  12, 1 HE  4)  1 MG 24
            voca = v(ioca)
            vcc = abs((vcca+vccp)*yc12*0.5d0+(voca+vocp+vocn)*yo16)
            dtcc = ftsc/vcc
            if (dtcc.le.dtnx.and.yc12.gt.1.d-4) then
               dtnx = dtcc
               idtnx = k
               iburning = 6
            endif
         endif

*________________________________________
***   time-step constrained by Ne-burning
*----------------------------------------

         if (tk.gt.TburnNe) then
***   1 Ne 20 ( 0 gamma, 1 alpha)  1 O  16
            vne = v(inega)
***   1 O  16 ( 1 alpha, 0 gamma)  1 Ne 20
***   inverse only accounted for if T > 1.15 TburnNe     
            if (tk.gt.1.15d0*TburnNe) then
               vneinv = v(ioag)*yhe4*yo16/yne20
            else
               vneinv = 1.d-99
            endif
            dtne = ftsh0/abs(vne-vneinv)
            if (dtne.le.dtnx.and.yne20.gt.1.d-4) then
               dtnx = dtne
               idtnx = k
               iburning = 7
            endif
         endif

*__________________________________________
***   determine smallest burning time scale
*------------------------------------------

         if (k.eq.itmax) then
c           if (dtsi.lt.dtoo.and.dtsi.lt.dtne.and.dtsi.lt.dtcc)
c    &           nphase = 9
c           if (dtoo.lt.dtsi.and.dtoo.lt.dtne.and.dtoo.lt.dtcc)
c    &           nphase = 8
            if (dtne.lt.dtsi.and.dtne.lt.dtoo.and.dtne.lt.dtcc.and.
     &           tmax.gt.TburnNe) nphase = 7
            if (dtcc.lt.dtsi.and.dtcc.lt.dtoo.and.dtcc.lt.dtne.and.
     &           tmax.gt.TburnC) nphase = 6
         endif

 20   continue

*_________________________________________________________
***   determination of the current stellar evolution phase
*---------------------------------------------------------

      if (tmax.lt.TburnC) then
         if (xsp(1,ihe4).lt.1.d-4) then
            nphase = 5
         else
            rcno = xsp(1,ic12)+xsp(1,in14)+xsp(1,io16)
            rcno = rcno/(xsp(nmod1,ic12)+xsp(nmod1,in14)+
     &           xsp(nmod1,io16))
            if (xsp(1,ih1).lt.1.d-4.and.tmax.gt.8.d7.and.(rcno.gt.
     &           1.05d0.or.totm.gt.8.d0)) then
               nphase = 4
            else
               if (xsp(1,ih1).lt.1.d-7) then
                  nphase = 3
               else
                 if (xsp(1,ih1)/xsp(nmod1,ih1).lt.0.998d0.or.
     &                 t(1).gt.3.d7) then
                     nphase = 2
                  else
                     nphase = 1
                  endif
               endif
            endif
         endif
      endif
      if (nphase.lt.nphase0) then
         if (nphase.eq.5) then
c            nphase = nphase0
            write (nout,*) 'star entering the SUPER-AGB phase'
         else
            write (nout,'(" WARNING : PHASE decreasing !! : nphase =",
     &           i2," -->",i2)') nphase0,nphase
         endif
      endif
      if (mtini.gt.12.d0) nphase = max(nphase0,nphase,1)

      if (nphase.le.2.and.lmicro.ge.2.and.lmicro.le.4.and.idiffcc)         ! Modif AP TD 06/2018 puis TD 11/18
     &     microdiffus = .true.

      return
      end
