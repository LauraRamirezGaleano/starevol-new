      SUBROUTINE prvar

************************************************************************
* Write results in evolhr, evolvar* and evolch* evolution files        *
* Modifs CC ondes (23/11/07)                                           *
* $LastChangedDate:: 2018-01-22 10:23:29 +0000 (Mon, 22 Jan 2018)     $ *
* $Author:: amard                                                    $ *
* $Rev:: 112                                                         $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.eng'
      include 'evolcom.eos'
      include 'evolcom.grad'
      include 'evolcom.lum'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.rot'
      include 'evolcom.conv'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.transp'
      include 'evolcom.var'
      include 'evolcom.ion'

      include 'evolcom.igw'

      character omegaconv*1,omegaconv0*1

C modif NL
      integer nzoneth,itth,ibth
C... fin modif NL
      integer icore,ienv,ib,it,jtot,nnmax,nt,inuctop,jth,nnmaxth
      integer i,i1,j,k,l,kp1,kl,ii,ip
      integer ibC,imC,itC,ibNe,imNe,itNe
      integer ibH,imH,itH,ibHe,imHe,itHe,imax,nzone,idiffex,normx
      integer klenv,kbce,ktce,bce,k_he,kconv,pass,klpulse,klcore

      integer iromax,irohalf,ndtenv

      double precision xmm(nmaxconv,3),xrr(nmaxconv,3),xmmth(nmaxconv,3)
     $     ,xrrth(nmaxconv,3)
      double precision shlimH,shlimHe,shlimC,shlimNe,shlimO
      double precision tmax,xmtmax,xmHb,xtHb,xroHb,xrHb,xmHm,xmHt,xtHt
     &     ,xroHt,xrHt,xmHeb,xtHeb,xroHeb,xrHeb,xmHem,xmHet,xtHet
     &     ,xroHet,xrHet,emaxHe,emaxH,xmb,xtb,xrobenv,xrb,xmt,xtt
     &     ,xrotenv,xrt,xmenvb,xtenvb,xroenvb,xrenvb,xmCt,xtCt,xroCt
     &     ,xrCt,xmCb,xtCb,xroCb,xrCb,emaxC,enumaxC,etamC,xmNeb,xrNeb
     &     ,xtNeb,xroNeb,xmNet,xrNet,xtNet,xroNet,xmNem,emaxNe,enumaxNe
     &     ,etamNe,enumaxHe,etamHe,xmCm,xmxs,xrxs,summ,xnuc
      double precision xmxsth,xrxsth,xmthb,xmtht
     $     ,xrthb,xrtht
      double precision lpph,rhomax,irradmax,lnu,leloc,leconv,lcno
      double precision xmoment,momenvel,xmoment_vieux,varmom
      double precision  epp(nsh),eCNO(nsh),e3a(nsh),eC12(nsh),
     &     eNe20(nsh),eO16(nsh),eppcno(nsh)
      double precision depot,dekin,deint,eenutri,eephot,lnuc,etot,edyn,
     &     totequil,vlnuc,eenuc,xmrad
      double precision eepot,veepot,eekin,veekin,eeint,veeint,eevisc
      double precision y(nsp),v(nreac)
      double precision diffst,diffcst
      double precision tautime,Ereac,taureac,rok
      double precision rossby,tconv_midCE
      double precision jrad,jconv,jtrue,jobs

      double precision romax,rohalf
      double precision omegamax,Nmax,brunt,BVmax,Nu2max,somegam
      double precision omcore,rmid(nsh),rmax,rmin,rcav

C Modif LA Coupling timescale (01/15)
      integer ndt,nshmax
      double precision Dshear,Fcirc,Fshear,tautot,taucirc,taushear,tauJ

C Modif LA convective core tc (03/16)
      double precision tg_core,tc_max_core,rc_max_core,tc_m_core,
     &     rc_m_core,rc_r_core,tc_r_core,tc_hp_core,rc_hp_core,tc_core,
     &     rc_core,omega_core,mcc,rcc

C Modif CC energie ondes (11/10/07)
      double precision londes
      double precision hpbce,rbce,hptce,factdc,tc,rc,tc_hp,rc_hp,rce
      double precision tc_r,rc_r,mce,tc_m,rc_m,tc_max,rc_max,tg,secd
      double precision nt2(nsh),nu2(nsh),bruntv2(nsh),bruntv(nsh),sum
      double precision dnu,dnuech,nu_max,errr,t_tot,sum_bce,t_bce,sum_he
      double precision r_he,t_he,sum_dp,dp
      double precision omegacore,omegajrad
      double precision Bcrit,fstar,Bequi
      double precision tauc

      character*2 cseq
      logical diffzc,thermal_equilibrium

      common / rossby_number / tconv_midCE
      common /omega_convection/ omegaconv,omegaconv0
      common /overshoot/ klcore,klenv,klpulse
      common /momreel/ xmoment,momenvel,xmoment_vieux
      common /turbulence/ diffst,diffcst,idiffex,diffzc,
     &     thermal_equilibrium
      common /nucleaire/ tautime(nsh,nreac),Ereac(nsh,nreac),taureac(nsh
     $     ,nreac)
      common /envel/ ndtenv
      common /brunt/ brunt(nsh),bruntv2,bruntv
      common /boreas_Mloss/ Bcrit,fstar,Bequi
      common /tauc_hp/ tauc


*____________________________________
***   print banner of evolution files
*------------------------------------

      if (no.eq.1.and.maxmod.eq.1.and.imodpr.le.0) then
         write (10,100) code_version
         write (11,110) code_version
         write (12,120) code_version
         write (13,130) code_version
         write (14,140) code_version
         write (15,150) code_version !H
         write (16,160) code_version !He
         write (17,162) code_version !C
         write (18,164) code_version !Ne
         write (27,170) code_version
         write (28,180) code_version

         write (29,260) code_version
         write (30,260) code_version

         write (21,190) code_version,'Central '
         write (22,200) code_version,'Central '
         write (123,210) code_version,'Central '
         write (124,220) code_version,'Central '
         write (31,190) code_version,'Surface '
         write (32,200) code_version,'Surface '
         write (33,210) code_version,'Surface '
         write (34,220) code_version,'Surface '
         write (41,230) code_version
         write (42,240) code_version
         write (43,250) code_version
      endif
      cseq = '  '
      if (model.eq.modeli+1) cseq = ' @'

*__________________________________________________________________
***   restore the accretion equilibrium abundance of H1, H2 and He3
*------------------------------------------------------------------

      if (dmaccr.gt.0.d0.and.accphase.le.1) then
         do k = iaccbot,iacctop
            xsp(k,ih2) = vxsp(k,ih2)
            xsp(k,ihe3) = vxsp(k,ihe3)
            normx = 1
            summ = 0.d0
            do l = 1,nsp
               xsp(k,l) = max(xsp(k,l),1.d-50)
               if (xsp(k,l).gt.xsp(k,normx)) normx = l
               summ = summ+xsp(k,l)
            enddo
            xsp(k,normx) = xsp(k,normx)+(1.d0-summ)
         enddo
      endif


c..   evolhr
      if (model.eq.modeli+1) then
         if (soll.le.99.9999d0.and.solr.le.9.9999d0.and.r(nmod)/rsun.le.
     &        9.9999d0) then
            write (10,1000) model,nphase,soll,solr,r(nmod)/rsun,teffy,
     &           roeff,geff,dms,totm,dtn*seci,time*seci,min(iter,99),
     &           ireset,nmod,cputime,cseq
         else
            if (r(nmod)/rsun.lt.9999.9d0) then
               write (10,1100) model,nphase,soll,solr,r(nmod)/rsun,
     &              teffy,roeff,geff,dms,totm,dtn*seci,time*seci,
     &              min(iter,99),ireset,nmod,cputime,cseq
            else
!               write (10,1150) model,nphase,min(soll,9999999.d0),
               write (10,1150) model,nphase,soll,
     &              min(9999.d0,solr),r(nmod)/rsun,teffy,roeff,geff,
     &              dms,totm,dtn*seci,time*seci,min(iter,99),ireset,
     &              nmod,cputime,cseq
            endif
         endif
      else
         if (soll.le.99.9999d0.and.solr.le.9.9999d0.and.r(nmod)/rsun.le.
     &        9.9999d0) then
            write (10,1200) model,nphase,soll,solr,r(nmod)/rsun,teffy,
     &           roeff,geff,dms,totm,dtn*seci,time*seci,min(iter,99),
     &           ireset,nmod,cputime
         else
            if (r(nmod)/rsun.lt.9999.9d0) then
!               write (10,1300) model,nphase,min(soll,9999999.d0),solr,
               write (10,1300) model,nphase,soll,solr,
     &              r(nmod)/rsun,teffy,roeff,geff,dms,totm,dtn*seci,
     &              time*seci,min(iter,99),ireset,nmod,cputime
            else
!               write (10,1350) model,nphase,min(soll,9999999.d0),
               write (10,1350) model,nphase,soll,
     &              min(9999.d0,solr),r(nmod)/rsun,teffy,roeff,geff,
     &              dms,totm,dtn*seci,time*seci,min(iter,99),ireset,
     &              nmod,cputime
            endif
         endif
      endif

*__________________________________________
***   determination of the max. temperature
*------------------------------------------

      imax = 1
      tmax = t(1)
      xmtmax = 0.d0
      do i = 1,nmod1
         if (t(i).gt.tmax) then
            tmax = t(i)
            imax = i
         endif
      enddo
      xmtmax = m(imax)/msun
      rhomax = ro(imax)


c. Modif LA 15/01/2013
*______________________________________________
***   determination of the half-max. density
*----------------------------------------------

      romax = 0.d0
      do i = 1,nmod1
         if (ro(i).gt.romax) then
            romax = ro(i)
            iromax = i
         endif
      enddo
      ii = 0
      rohalf = romax
      do while (rohalf.gt.(romax/2.d0))
         ii = ii+1
         rohalf = ro(ii)
         irohalf = ii
      enddo
      if (irohalf.ge.ndtenv) irohalf = ndtenv

c. Fin Modif LA

*____________________________________________________________________
***   calculation of global energetic characteristics of the new star
*--------------------------------------------------------------------

      eepot = 0.d0
      veepot = 0.d0
      eekin = 0.d0
      veekin = 0.d0
      eeint = 0.d0
      veeint = 0.d0
      eevisc = 0.d0
      lnuc = 0.d0
C Modif CC energie ondes (11/10/07)
C*      londes = 0.d0
      vlnuc = 0.d0
      eenutri = 0.d0
      do k = 2,nmod1
         eepot = eepot-dm(k)*(m(k+1)/r(k+1)+m(k)/r(k))
         veepot = veepot-dm(k)*(m(k+1)/vr(k+1)+m(k)/vr(k))
         if (hydro) then
            eekin = eekin+(u(k)+u(k+1))**2*dm(k)
            veekin = veekin+(vu(k)+vu(k+1))**2*dm(k)
         endif
         eeint = eeint+e(k)*dm(k)
         veeint = veeint+ve(k)*dm(k)
         eevisc = eevisc+evisc(k)*dm(k)
         lnuc = lnuc+enucl(k)*dm(k)
         vlnuc = vlnuc+venucl(k)*dm(k)
         eenutri = eenutri+enupla(k)*dm(k)
C Modif CC energie ondes (11/10/07)
C*         londes = londes+eondes(k)*dm(k)
      enddo
      eepot = eepot*g*0.5d0
      veepot = veepot*g*0.5d0
      eekin = eekin*0.125d0
      veekin = veekin*0.125d0
      eenuc = 0.5d0*(lnuc+vlnuc)*dtn
      eenutri = eenutri*dtn
      eephot = 0.5d0*(lum(nmod)+vlum(nmod))*dtn
      depot = eepot-veepot
      dekin = eekin-veekin
      deint = eeint-veeint
      edyn = depot+dekin+deint
      etot = edyn+eephot
      totequil = abs(1.d0-eenuc/etot)
      if (abs(totequil).gt.99.99d0) totequil = 99.9999d0
c..   dEtot = d(Ekin+Epot+Eint+Ephot) = dEnuc
      write (90,1450) depot,dekin,deint,eephot,eenutri,eenuc,totequil

      tkh = m(nmod)*m(nmod)*g/abs(r(nmod)*lum(nmod))

*______________________________________________________________________
***   calculation of the momentum of inertia, of the gravothermal and
***   nuclear luminosities associated to H, He, C, Ne, O and Si-burning
***   and neutron irradiation (for He-burning zones)
*----------------------------------------------------------------------

      vlpp = lpp
      vlh = lh
      vlhe = lhe
      vlc = lc
      vlne = vlne
      vlo = vlo
      vlsi = vlsi

c     Determination of the nuclear energy contributions

c... modif NL 6/11/2008
	do k=1,nmod
           rok = ro(k)
            call vit (t(k),ro(k),mueinv(k),k,v,0)
            call nucprod (v,y,t(k),epp(k),eCNO(k),e3a(k),eC12(k),
     &           eNe20(k),eO16(k))
***   2 PROT  ( 0 OOOOO, 0 OOOOO)  1 DEUT
            if (v(2).ne.0.d0) then
               tautime(k,ippg)=(2.d0)/(v(2)*xsp(k,ih1)*sec)
            else
               tautime(k,ippg) = 1.d99
            endif
            tautime(k,ippg) = min(tautime(k,ippg),1.d20)
            Ereac(k,ippg)=xsp(k,ih1)**2*v(2)*qi(ippg)*avn/(2.d0)
            taureac(k,ippg)=avn*rok*xsp(k,ih1)/(tautime(k,ippg))
***   1 DEUT  ( 1 PROT , 0 OOOOO)  1 HE  3
            if (v(3).ne.0.d0) then
               tautime(k,ipdg)=2.d0/(v(3)*xsp(k,ih1)*sec)
            else
               tautime(k,ipdg) = 1.d99
            endif
            tautime(k,ipdg) = min(tautime(k,ipdg),1.d20)
            Ereac(k,ipdg)=xsp(k,ih1)*xsp(k,ih2)*v(3)*qi(ipdg)*avn/(2.d0)
            taureac(k,ipdg)=avn*rok*xsp(k,ih2)/(2*tautime(k,ipdg))
***   2 HE  3 ( 0 OOOOO, 2 PROT )  1 HE  4
            if (v(8).ne.0.d0) then
               tautime(k,i2he3)=3*2.d0/(v(8)*xsp(k,ihe3)*sec)
            else
               tautime(k,i2he3) = 1.d99
            endif
            tautime(k,i2he3) = min(tautime(k,i2he3),1.d20)
            Ereac(k,i2he3)=xsp(k,ihe3)**2*v(8)*qi(i2he3)*avn/(9*2.d0)
            taureac(k,i2he3)=avn*rok*xsp(k,ihe3)/(3*tautime(k,i2he3))
***   1 HE  4 ( 1 HE  3, 0 OOOOO)  1 BE  7
            if (v(10).ne.0.d0) then
               tautime(k,ihe3ag)=4.d0/(v(10)*xsp(k,ihe4)*sec)
            else
               tautime(k,ihe3ag) = 1.d99
            endif
             tautime(k,ihe3ag) = min(tautime(k,ihe3ag),1.d20)
            Ereac(k,ihe3ag)=xsp(k,ihe3)*xsp(k,ihe4)*v(10)*qi(ihe3ag)
     &                    *avn/(12.d0)
            taureac(k,i2he3)=avn*rok*xsp(k,ihe3)/(3*tautime(k,ihe3ag))
***   1 BE  7 ( 0 betap, 0 nutri)  1 LI  7
            if (v(ibe7beta).ne.0.d0) then
               tautime(k,ibe7beta)=1.d0/(v(ibe7beta)*sec)
            else
               tautime(k,ibe7beta) = 1.d99
            endif
            tautime(k,ibe7beta) = min(tautime(k,ibe7beta),1.d20)
            Ereac(k,ibe7beta)=xsp(k,ibe7)*v(ibe7beta)*qi(ibe7beta)
     &      *avn/(7.d0)
            taureac(k,ibe7beta)=avn*rok*xsp(k,ibe7)/(7*tautime(k
     $           ,ibe7beta))
***   1 LI  7 ( 1 PROT , 0 OOOOO)  2 HE  4
            if (v(14).ne.0.d0) then
               tautime(k,ili7pa)=1.d0/(v(14)*xsp(k,ih1)*sec)
            else
               tautime(k,ili7pa) = 1.d99
            endif
            tautime(k,ili7pa) = min(tautime(k,ili7pa),1.d20)
            Ereac(k,ili7pa)=xsp(k,ili7)*xsp(k,ih1)*v(14)*qi(ili7pa)*avn
     &                     /(7.d0)
            taureac(k,ili7pa)=avn*rok*xsp(k,ili7)/(7*tautime(k,ili7pa))
***   1 BE  7 ( 1 PROT , 0 OOOOO)  1 B   8
            if (v(17).ne.0.d0) then
            tautime(k,ibe7pg)=1.d0/(v(17)*xsp(k,ih1)*sec)
            else
               tautime(k,ibe7pg) = 1.d99
            endif
            tautime(k,ibe7pg) = min(tautime(k,ibe7pg),1.d20)
            Ereac(k,ibe7pg)=xsp(k,ibe7)*xsp(k,ih1)*v(17)*qi(ibe7pg)*avn
     &                     /(7.d0)
            taureac(k,ibe7pg)=avn*rok*xsp(k,ibe7)/(7*tautime(k,ibe7pg))
***   1 B   8 ( 0 OOOOO, 0 OOOOO)  2 HE  4
            if (v(ib8beta).ne.0.d0) then
               tautime(k,ib8beta)=1.d0/(v(ib8beta)*sec)
            else
               tautime(k,ib8beta) = 1.d99
            endif
            tautime(k,ib8beta) = min(tautime(k,ib8beta),1.d20)
            Ereac(k,ib8beta)=xsp(k,ib8)*v(ib8beta)*qi(ib8beta)*avn
     $           /(8.d0)
            taureac(k,ib8beta)=avn*rok*xsp(k,ib8)/(8*tautime(k,ib8beta))
***   1 C  13 ( 1 PROT , 0 OOOOO)  1 N  14
            if (v(35).ne.0.d0) then
               tautime(k,ic13pg)=1.d0/(v(35)*xsp(k,ih1)*sec)
            else
               tautime(k,ic13pg) = 1.d99
            endif
            tautime(k,ic13pg) = min(tautime(k,ic13pg),1.d20)
            Ereac(k,ic13pg)=xsp(k,ic13)*xsp(k,ih1)*v(35)*qi(ic13pg)*avn
     &                    /(13.d0)
            taureac(k,ic13pg)=avn*rok*xsp(k,in14)/(13*tautime(k,ic13pg))
***   1 N  14 ( 1 PROT , 0 OOOOO)  1 O  15
            if (v(48).ne.0.d0) then
               tautime(k,in14pg)=1.d0/(v(48)*xsp(k,ih1)*sec)
            else
               tautime(k,in14pg) = 1.d99
            endif
            tautime(k,in14pg) = min(tautime(k,in14pg),1.d20)
            Ereac(k,in14pg)=xsp(k,in14)*xsp(k,ih1)*v(48)*qi(in14pg)*avn
     &                    /(14.d0)
            taureac(k,in14pg)=avn*rok*xsp(k,in14)/(14*tautime(k,in14pg))
***   1 C  12 ( 1 PROT , 0 OOOOO)  1 N  13
            if (v(30).ne.0.d0) then
               tautime(k,icpg)=1.d0/(v(30)*xsp(k,ih1)*sec)
            else
               tautime(k,icpg) = 1.d99
            endif
            tautime(k,icpg) = min(tautime(k,icpg),1.d20)
            Ereac(k,icpg)=xsp(k,ic12)*xsp(k,ih1)*v(30)*qi(icpg)*avn
     &                    /(12.d0)
            taureac(k,icpg)=avn*rok*xsp(k,ic12)/(12*tautime(k,icpg))
c...fin de modif
	enddo

      inuctop = nmod1
      shlimH = 0.d0
      shlimHe = 0.d0
      shlimC = 0.d0
      shlimNe = 0.d0
      shlimO = 0.d0
      do k = 1,nmod1
         if (t(k).ge.tnucmin) then
            do l = 1,nsp
               y(l) = ysp(k,l)
            enddo
            call vit (t(k),ro(k),mueinv(k),k,v,0)
            call nucprod (v,y,t(k),epp(k),eCNO(k),e3a(k),eC12(k),
     &           eNe20(k),eO16(k))
            eppcno(k) = epp(k) + eCNO(k)
            shlimH = max(shlimH,eppcno(k))
            shlimHe = max(shlimHe,e3a(k))
            shlimC = max(shlimC,eC12(k))
            shlimNe = max(shlimNe,abs(eNe20(k)))
            shlimO = max(shlimO,eO16(k))
          else
            inuctop = k-1
            goto 111
         endif
      enddo

 111  irradmax = -1.d0
      do k = 1,inuctop
         if (irradmax.lt.exposure(k)) then
            irradmax = exposure(k)
            xmrad = m(k)/msun
         endif
      enddo

c..   evolvar1
      if (eta(1).lt.9999.9d0) then
         write (11,1400) model,t(1),tmax,xmtmax,ro(1),rhomax,p(1),
     &        beta(1),eta(1),degpe(1),enupla(1),enucl(1),egrav(1),cseq
      else
         write (11,1410) model,t(1),tmax,xmtmax,ro(1),rhomax,p(1),
     &        beta(1),eta(1),degpe(1),enupla(1),enucl(1),egrav(1),cseq
      endif


*___________________________________________________________
***   determination of the nuclear burning region boundaries
***   and properties
*-----------------------------------------------------------

      shlimH = min(shlimH,shlim)
      shlimHe = min(shlimHe,shlim)
      if (nphase.lt.2) shlimHe = 1.d10
      if (nphase.ge.5) then
         shlimC = min(shlimC,shlim)
         if (nphase.ge.6) then
            shlimNe = min(abs(shlimNe),shlim)
            shlimO = min(shlimO,shlim)
         endif
      else
         shlimC = shlim
         shlimNe = shlim
         shlimO = shlim
      endif
      shlimNe = max(shlimNe,1.d-10)

c Burning regions for PP, CNO, 3a, C12 and Ne20

      ibH = 0
      itH = 0
      imH = 0
      xmHm = 0.d0
      emaxH = 0.d0
c      enumaxH = 0.d0
c      etamH = 0.d0
      xmHb = 0.d0
      xrHb = 0.d0
      xtHb = 0.d0
      xroHb = 0.d0
      xmHt = 0.d0
      xrHt = 0.d0
      xtHt = 0.d0
      xroHt = 0.d0

      ibHe = 0
      itHe = 0
      imHe = 0
      xmHem = 0.d0
      emaxHe = 0.d0
      enumaxHe = 0.d0
      etamHe = 0.d0
      xmHeb = 0.d0
      xrHeb = 0.d0
      xtHeb = 0.d0
      xroHeb = 0.d0
      xmHet = 0.d0
      xrHet = 0.d0
      xtHet = 0.d0
      xroHet = 0.d0

      ibC = 0
      itC = 0
      imC = 0
      xmCm = 0.d0
      emaxC = 0.d0
      enumaxC = 0.d0
      etamC = 0.d0
      xmCb = 0.d0
      xrCb = 0.d0
      xtCb = 0.d0
      xroCb = 0.d0
      xmCt = 0.d0
      xrCt = 0.d0
      xtCt = 0.d0
      xroCt = 0.d0

      ibNe = 0
      itNe = 0
      imNe = 0
      xmNem = 0.d0
      emaxNe = 0.d0
      enumaxNe = 0.d0
      etamNe = 0.d0
      xmNeb = 0.d0
      xrNeb = 0.d0
      xtNeb = 0.d0
      xroNeb = 0.d0
      xmNet = 0.d0
      xrNet = 0.d0
      xtNet = 0.d0
      xroNet = 0.d0
      xmNet = 0.d0
      xrNet = 0.d0
      xtNet = 0.d0
      xroNet = 0.d0

      xnuc = 1.5d0

      do k = 1,inuctop-1
c..   location of maximum energy production rate
         if (eppcno(k).ge.shlimH.and.eppcno(k).gt.emaxH) then
            imH = k
            emaxH = eppcno(imH)
         endif
         if (e3a(k).ge.shlimHe.and.e3a(k).gt.emaxHe) then
            imHe = k
            emaxHe = e3a(imHe)
         endif
         if (nphase.ge.5) then
            if (eC12(k).ge.shlimC.and.eC12(k).gt.emaxC.and.eC12(k).gt.
     &           xnuc*abs(eNe20(k))) then
               imC = k
               emaxC = eC12(imC)
            endif
            if (abs(eNe20(k)).ge.shlimNe.and.abs(eNe20(k)).gt.emaxNe
     &           .and.((nphase.le.6.and.abs(eNe20(k)).gt.xnuc*eC12(k))
     &           .or.(abs(eNe20(k)).gt.xnuc*eO16(k).and.nphase.gt.6)))
     &           then
               imNe = k
               emaxNe = abs(eNe20(imNe))
            endif
         endif
c..      top of burning shells
         kp1 = k+1
         if (eppcno(k).ge.shlimH.and.eppcno(kp1).lt.shlimH) itH = k
         if (e3a(k).ge.shlimHe.and.e3a(kp1).lt.shlimHe) itHe = k
         if (eC12(k).ge.shlimC.and.eC12(kp1).lt.shlimC) itC = k
         if (abs(eNe20(k)).ge.shlimNe.and.abs(eNe20(k)).gt.xnuc*eC12(k))
     &        itNe = k
      enddo

      itH = max(itH,imH)
      itHe = max(itHe,imHe)
      itC = max(itC,imC)
      itNe = max(itNe,imNe)

      do k = 1,inuctop-1
        kp1 = k+1
c..     base of burning shells
        if (imH.gt.0) then
           if (eppcno(k).ge.shlimH.and.k.eq.1) then
              ibH = k
           else if (eppcno(k).lt.shlimH.and.eppcno(kp1).ge.shlimH.and.
     &             ibH.eq.0) then
              ibH = kp1
           endif
        endif
        if (e3a(k).ge.shlimHe.and.k.eq.1) then
           ibHe = k
        else if (e3a(k).lt.shlimHe.and.e3a(kp1).ge.shlimHe.and.
     &          ibHe.eq.0) then
           ibHe = kp1
        endif
        if (itC.gt.0.and.ibC.eq.0.and.eC12(k).ge.shlimC.and.
     &       eC12(k).ge.xnuc*abs(eNe20(k)).and.k.gt.itNe) then
           ibC = k
        endif
        if (itNe.gt.0.and.ibNe.eq.0.and.abs(eNe20(k)).ge.shlimNe.and.
     &       ((nphase.le.6.and.abs(eNe20(k)).gt.xnuc*eC12(k)).or.
     &       (abs(eNe20(k)).ge.xnuc*eO16(k).and.nphase.gt.6))) then
           ibNe = k
        endif
      enddo

      if (itHe.ne.0.and.ibH.ne.0) itHe = min(max(ibH,imHe),itHe)
c      if (itC.ne.0) itC = min(ibHe,itC)
      if (itNe.ne.0) itNe = min(ibC,itNe)

      if (ibH.eq.0.and.itH.ne.0) stop 'prvar : error HBS boundaries'
      if (ibHe.eq.0.and.itHe.ne.0.or.imHe.gt.itHe) then
         print *,'HeBS : [',ibHe,':',imHe,':',itHe,']'
         stop 'prvar : error HeBS boundaries'
      endif
      if (ibC.eq.0.and.itC.ne.0.or.imC.gt.itC) then
         stop 'prvar : error CBS boundaries'
      endif
      if (ibNe.eq.0.and.itNe.ne.0.or.imNe.gt.itNe)
     &     stop 'prvar : error NeBS boundaries'

      r(1) = 0.d0
      lh = 0.d0
      lpp = 0.d0
      lcno = 0.d0
      lhe = 0.d0
      lc = 0.d0
      lne = 0.d0
      lo = 0.d0
      lsi = 0.d0
      lnu = 0.d0

c hydrogen burning shell
      if (ibH.gt.0) then
         xmHb = m(ibH)/msun
         xrHb = r(ibH)/r(nmod)
         xtHb = t(ibH)
         xroHb = ro(ibH)
         xmHt = m(itH)/msun
         xrHt = r(itH)/r(nmod)
         xtHt = t(itH)
         xroHt = ro(itH)
         xmHm = m(imH)/msun
         emaxH = enucl(imH)
c         enumaxH = enupla(imH)
c         etamH = eta(imH)
         do i = ibH,itH
            lh = lh+(enucl(i)+enupla(i))*dm(i)
            lcno = lcno+eCNO(i)*dm(i)
            lpp = lpp+epp(i)*dm(i)
         enddo
      endif

c helium burning shell
      if (ibHe.gt.0) then
         xmHeb = m(ibHe)/msun
         xrHeb = r(ibHe)/r(nmod)
         xtHeb = t(ibHe)
         xroHeb = ro(ibHe)
         xmHet = m(itHe)/msun
         xrHet = r(itHe)/r(nmod)
         xtHet = t(itHe)
         xroHet = ro(itHe)
         xmHem = m(imHe)/msun
         emaxHe = enucl(imHe)
         enumaxHe = enupla(imHe)
         etamHe = eta(imHe)
         do i = ibHe,itHe
            lhe = lhe+(enucl(i)+enupla(i))*dm(i)
         enddo
      endif

c carbon burning shell
      if (ibC.gt.0) then
         xmCb = m(ibC)/msun
         xrCb = r(ibC)/r(nmod)
         xtCb = t(ibC)
         xroCb = ro(ibC)
         xmCt = m(itC)/msun
         xrCt = r(itC)/r(nmod)
         xtCt = t(itC)
         xroCt = ro(itC)
         xmCm = m(imC)/msun
         emaxC = enucl(imC)
         enumaxC = enupla(imC)
         etamC = eta(imC)
         do i = ibC,itC
            lc = lc+(enucl(i)+enupla(i))*dm(i)
         enddo
      endif

c neon burning shell
      if (ibNe.gt.0) then
         xmNeb = m(ibNe)/msun
         xrNeb = r(ibNe)/r(nmod)
         xtNeb = t(ibNe)
         xroNeb = ro(ibNe)
         xmNet = m(itNe)/msun
         xrNet = r(itNe)/r(nmod)
         xtNet = t(itNe)
         xroNet = ro(itNe)
         xmNem = m(imNe)/msun
         emaxNe = enucl(imNe)
         enumaxNe = enupla(imNe)
         etamNe = eta(imNe)
         do i = ibNe,itNe
            lne = lne+(enucl(i)+enupla(i))*dm(i)
         enddo
      endif

      totgrav = 0.d0
      totnucl = 0.d0
c      lshr = 0.d0

      do i = 1,nmod
         totnucl = totnucl+enucl(i)*dm(i)
         totgrav = totgrav+egrav(i)*dm(i)
         lnu = lnu+enupla(i)*dm(i)
c         lshr = lshr+eshr(i)*dm(i)*thacc
      enddo

      leloc = 0.d0
      leconv = 0.d0
      if (urca) then
         do i = 1,inuctop
            leloc = leloc+eloc(i)*dm(i)
            leconv = leconv+egconv(i)*dm(i)
         enddo
         leloc = leloc/lsun
         leconv = leconv/lsun
         leloc = sign(max(abs(leloc),1.d-30),leloc)
         leconv = sign(max(abs(leconv),1.d-30),leconv)
      endif

      totnucl = totnucl/lsun
      totgrav = totgrav/lsun
      lh = lh/lsun
      lpp = lpp/lsun
      lhe = lhe/lsun
      lc = lc/lsun
      lne = lne/lsun
c      lo = lo/lsun
c      lsi = lsi/lsun
      lnu = lnu/lsun
      lh = sign(max(abs(lh),1.d-30),lh)
      lhe = sign(max(abs(lhe),1.d-30),lhe)
      lc = sign(max(abs(lc),1.d-30),lc)
      lne = sign(max(abs(lne),1.d-30),lne)
c      lo = sign(max(abs(lo),1.d-30),lo)
c      lsi = sign(max(abs(lsi),1.d-30),lsi)
c      lshr = lshr/lsun

c..   evolvar2
      write (12,1500) model,lh,lhe,lc,lne,lo,xmrad,lnu,totnucl,
     &     totgrav,irradmax,cseq
c      write (12,1500) model,lh,lhe,lc,lne,lo,leloc,lnu,totnucl,
c     &     totgrav,leconv,cseq

c..   evolvar5
      if (lcno.eq.0.d0.and.lpp.eq.0.d0) then
         lpph = 0.d0
      else
         lpph = max(log(lpp/lcno),-99.9998d0)
      endif
      write (15,1600) model,xmHb,xrHb,xtHb,xroHb,xmHt,xrHt,xtHt,xroHt,
     &     xmHm,emaxH,lpph,cseq
c..   evolvar6
      write (16,1700) model,xmHeb,xrHeb,xtHeb,xroHeb,xmHet,xrHet,xtHet,
     &     xroHet,xmHem,emaxHe,enumaxHe,etamHe,cseq
c..   evolvar7
      write (17,1700) model,xmCb,xrCb,xtCb,xroCb,xmCt,xrCt,xtCt,
     &     xroCt,xmCm,emaxC,enumaxC,etamC,cseq
c..   evolvar8
      write (18,1700) model,xmNeb,xrNeb,xtNeb,xroNeb,xmNet,xrNet,xtNet,
     &     xroNet,xmNem,emaxNe,enumaxNe,etamNe,cseq
c..   evolvar9
c      write (19,1700) model,xmOb,xrOb,xtOb,xroOb,xmOt,xrOt,xtOt,
c     &     xroOt,xmOm,emaxO,enumaxO,etamO,cseq
cc..   evolvar10
c      write (20,1700) model,xmSib,xrSib,xtSib,xroSib,xmSit,xrSit,xtSit,
c     &     xroSit,xmSim,emaxSi,enumaxSi,etamSi,cseq



*_______________________________________________________
***   determination of the convective regions properties
*-------------------------------------------------------

      r(1) = 0.d0
      xmb = 0.d0
      xtb = 0.d0
      xrobenv = 0.d0
      xrb = 0.d0
      xmt = 0.d0
      xtt = 0.d0
      xrotenv = 0.d0
      xrt = 0.d0
      xmenvb = 0.d0
      xtenvb = 0.d0
      xroenvb = 0.d0
      xrenvb = 0.d0
      icore = 0
      ienv = 0
      if (nsconv.eq.0) goto 80
      do kl = 1,nsconv
c         if (novlim(kl,3).gt.1.and.mr(novlim(kl,4)).gt.0.7d0.and.
         if (novlim(kl,3).gt.1.and.
     &        (tau(novlim(kl,4)).lt.1.d2.or.t(novlim(kl,4)).lt.1.d6))
     &        then
            ienv = kl
            goto 60
         endif
      enddo
 60   if (ienv.gt.0) then
         ib = novlim(ienv,3)
         xmenvb = m(ib)/msun
         xtenvb = log10(t(ib))
         xroenvb = log10(ro(ib))
         xrenvb = r(ib)/r(nmod)
      endif
      if (nsconv.gt.0.and.ienv.ne.1) then
         ib = novlim(1,3)
         it = novlim(1,4)
         icore = 1
         xmb = m(ib)/msun
         xtb = log10(t(ib))
         xrobenv = log10(ro(ib))
         xrb = r(ib)/r(nmod)
         xmt = m(it+1)/msun
         xtt = log10(t(it))
         xrotenv = log10(ro(it))
         xrt = r(it)/r(nmod)
      endif

c..   evolvar3
 80   write (13,1900) model,xmb,xrb,xtb,xrobenv,xmt,xrt,xtt,xrotenv,
     &     xmenvb,xrenvb,xtenvb,xroenvb,cseq

*_____________________________________________________
***   selection of the remaining convective zones
***   and sorting as a function of their relative mass
*-----------------------------------------------------

      jtot = 0
      do i = 1,nsconv
         if (i.ne.icore.and.i.ne.ienv.and.novlim(i,3).ne.0) then
            jtot = jtot+1
            xmm(jtot,1) = m(novlim(i,3))/msun
            xmm(jtot,2) = m(novlim(i,4)+1)/msun
            xmm(jtot,3) = xmm(jtot,2)-xmm(jtot,1)
            xrr(jtot,1) = r(novlim(i,3))/rsun
            xrr(jtot,2) = r(novlim(i,4)+1)/rsun
            xrr(jtot,3) = xrr(jtot,2)-xrr(jtot,1)
         endif
      enddo
      do i = 1,jtot
         do j = i+1,jtot
            if (xmm(j,3).gt.xmm(i,3)) then
               xmxs = xmm(j,1)
               xmm(j,1) = xmm(i,1)
               xmm(i,1) = xmxs
               xmxs = xmm(j,2)
               xmm(j,2) = xmm(i,2)
               xmm(i,2) = xmxs
               xmxs = xmm(j,3)
               xmm(j,3) = xmm(i,3)
               xmm(i,3) = xmxs
               xrxs = xrr(j,1)
               xrr(j,1) = xrr(i,1)
               xrr(i,1) = xrxs
               xrxs = xrr(j,2)
               xrr(j,2) = xrr(i,2)
               xrr(i,2) = xrxs
               xrxs = xrr(j,3)
               xrr(j,3) = xrr(i,3)
               xrr(i,3) = xrxs
            endif
         enddo
      enddo

*_____________________________________________________________
***   selection of the next four most massive convective zones
*-------------------------------------------------------------

      nnmax = 5
      if (jtot.lt.nnmax) nnmax = jtot
      do i = 1,nnmax
         do j = i,nnmax
            if (xmm(j,1).lt.xmm(i,1)) then
               xmxs = xmm(j,1)
               xmm(j,1) = xmm(i,1)
               xmm(i,1) = xmxs
               xmxs = xmm(j,2)
               xmm(j,2) = xmm(i,2)
               xmm(i,2) = xmxs
               xmxs = xmm(j,3)
               xmm(j,3) = xmm(i,3)
               xmm(i,3) = xmxs
               xrxs = xrr(j,1)
               xrr(j,1) = xrr(i,1)
               xrr(i,1) = xrxs
               xrxs = xrr(j,2)
               xrr(j,2) = xrr(i,2)
               xrr(i,2) = xrxs
               xrxs = xrr(j,3)
               xrr(j,3) = xrr(i,3)
               xrr(i,3) = xrxs
            endif
         enddo
      enddo
      do i = nnmax+1,nmaxconv
         xmm(i,1) = 0.d0
         xmm(i,2) = 0.d0
         xmm(i,3) = 0.d0
         xrr(i,1) = 0.d0
         xrr(i,2) = 0.d0
         xrr(i,3) = 0.d0
      enddo
      if (nsconv.le.0) then
         nzone = 1
      else
         nzone = nsconv
      endif

c... modif NL  19/09/08

         xmthb=0.d0
         xmtht=0.d0
	 xrthb=0.d0
	 xrtht=0.d0

         if (novlim(1,11).ge.1) then
            xmthb=m(novlim(1,11))/msun
            xmtht=m(novlim(1,12)+1)/msun
            xrthb=r(novlim(1,11))/rsun
            xrtht=r(novlim(1,12)+1)/rsun
         endif
	 write(*,*) xmthb,xmtht
C	 stop

c..   evolvar4
      write (14,2000) model,xmm(1,1),xmm(1,2),xmm(2,1),xmm(2,2),
     &     xmm(3,1),xmm(3,2),xmm(4,1),xmm(4,2),xmthb,xmtht,
     &     min(nzone,99),scr,cseq
c.. evolvar12

       write (28,2050) model,xrr(1,1),xrr(1,2),xrr(2,1),xrr(2,2),
     &     xrr(3,1),xrr(3,2),xrr(4,1),xrr(4,2),xrthb,xrtht,
     &     cseq


*____________________________________________________
***   storage of the accretion and rotation variables
*----------------------------------------------------

      k2rad = 0.d0
      k2conv = 0.d0
      jrad = 0.d0
      jconv = 0.d0
      angconvr = 0.d0
      angradr = 0.d0
      jobs = 0.d0
      if (dmaccr.eq.0.d0) then
         raccbot = 1.d0
         maccbot = 1.d0
         lacc = 0.d0
         facc(1:nmod) = 0.d0
      endif
      if (irotbin.eq.0) then
         Fenerg = 0.d0
         vsurf = 0.d0
         ur_Ss = 0.d0
      endif

***   Computing the radius of gyration of the radiative zone (k2rad)
***   and of the convective envelope (k2conv) according to Rucinski 1988 (AJ 95)

c$$$      if (ienv.gt.0.and.(omegaconv.ne.'s')) then
      if (ienv.gt.0) then
         if (novlim(ienv,3).gt.0) then
            print *, 'ienv',ienv,'novlim(ienv,3)',novlim(ienv,3)
     $           ,'novlim(1,4)',novlim(1,4) 
            do i = novlim(ienv,3),novlim(ienv,4)
               i1 = i-1
               k2conv = k2conv+dm(i1)*(r(i)**2+r(i1)**2+r(i)*r(i1))
     &              /3.d0
            enddo

            if (nsconv.eq.1) then
               do i = 2,novlim(ienv,3)-1
                  i1 = i-1
                  k2rad = k2rad+dm(i1)*(r(i)**2+r(i1)**2+r(i)*r(i1))
     &                 /3.d0
                  jrad = jrad+dm(i1)*(r(i)**2+r(i1)**2+r(i)*r(i1))
     &                 /3.d0*0.5d0*(omega(i1)+omega(i))
                  omegacore = omegacore+omega(i1)+omega(i)
               enddo
               omegajrad = jrad/k2rad
               omegacore = omegacore/(2*novlim(ienv,3)-1)
               print *, 'nsconv',nsconv,'omegajrad',omegajrad
     $              ,'omegacore',omegacore
            else if (nsconv.ge.2) then
               do i = novlim(1,4)+1,novlim(ienv,3)-1 ! Test TD pour prendre en compte la présence d'un coeur convectif 
                  i1 = i-1
                  k2rad = k2rad+dm(i1)*(r(i)**2+r(i1)**2+r(i)*r(i1))
     &                 /3.d0
                  jrad = jrad+dm(i1)*(r(i)**2+r(i1)**2+r(i)*r(i1))
     &                 /3.d0*0.5d0*(omega(i1)+omega(i))
                  omegacore = omegacore+omega(i1)+omega(i)
               enddo
               omegajrad = jrad/k2rad
               omegacore = omegacore/(2*(novlim(ienv,3)-novlim(1,4)+1)
     $              -1)
               print *, 'nsconv',nsconv,'omegajrad',omegajrad
     $              ,'omegacore',omegacore
            endif
         else
            do i = 2,nmod
               i1 = i-1
               k2conv = k2conv+dm(i1)*(r(i)**2+r(i1)**2+r(i)*r(i1))
     &              /3.d0
            enddo
            omegacore = omega_S
         endif
      else
         do i = 2,nmod
            i1 = i-1
            k2conv = k2conv+dm(i1)*(r(i)**2+r(i1)**2+r(i)*r(i1))
     &           /3.d0
         enddo
         omegacore = omega_S
         omegajrad = omega_S
      endif

      jrad = pw23*jrad
      jtrue = (pw23*k2conv*omega_S+jrad)/(totm*msun)
      jobs = (pw23*(k2conv+k2rad))*omega_S/(totm*msun)

      k2conv = k2conv*2.d0/(3.d0*totm*msun*(solr*rsun)**2)
      k2rad = k2rad*2.d0/(3.d0*totm*msun*(solr*rsun)**2)
      varmom = 0.d0


      if (irotbin.eq.1) then
         if (.not.diffzc.and.novlim(nsconv,3)-1.le.1.and.nphase.eq.1)
     &        then
            varmom = 1.d0
         else
            momenvel = momenvel
            xmoment = xmoment
            varmom = 1.d0-xmoment_vieux/xmoment
         endif
c     Rossby number for rotating models : Prot / tconv
c     Ro = 2*pi*R/(vsurf*tc)*8.055 with R in units of Rsun, vsurf in units of km/s and tc in units of days
         rossby = pim2/(omega_S*tconv_midCE)
         if (tconv_midCE.eq.0.d0) rossby = 0.d0
c..   evolvar11
      else
         Fenerg = 0.d0
         omega_S = 0.d0
         rossby = 0.d0
         vsurf = 0.d0
         xmoment = 0.d0
         dmoment1 = 0.d0
      endif
      write (27,2100) model,massrate,raccbot,maccbot,lacc/lsun,k2conv
     &     ,k2rad,rossby,omega_S,vsurf,xmoment,Fenerg,dmoment1,cseq



!!! Mean the rotation rate below 0.025 Rtot for asteroseismic comparison

      Nt2=0.d0
      Nu2=0.d0
      do i = 1,nmod
         bruntV2(i) = 0.d0
         bruntV(i) = 0.d0
      enddo

      do i=1,nmod-1
         grav(i) = gmr(i)/gmr(nmod)
         gs = gmr(nmod)

         Nt2(i) = gs*grav(i)*deltaKS(i)*(abad(i)-abla(i))/(hp(i))
         Nu2(i) = gs*grav(i)*phiKS(i)*abmu(i)/(hp(i))

         bruntV2(i) = Nt2(i)+Nu2(i)
          
      enddo
     
      nshmax = 1
      do i=1,nmod
         if (i.lt.novlim(klenv,3).and.(Nu2(i).gt.Nu2max)) then
            Nu2max = Nu2(i)
            nshmax = i
         endif
      enddo
      
      somegam = 0.d0
      do i=1,nshmax
        if (omega(i).gt.0.d0) somegam = somegam + omega(i)*dm(i)
      enddo
      if (nshmax.ne.1) then 
         omegamax = somegam/m(nshmax)
      else 
         omegamax = 0.d0
      endif

c$$$      rmid(1)=r(1)
c$$$      omcore = 0.d0
c$$$      do i=2,nmod
c$$$         i1 = i-1
c$$$         rmid(i) = dsqrt(pw13*(r(i)**2+r(i1)**2+r(i)*r(i1)))
c$$$      enddo
c$$$      BVmax = 0.d0
c$$$      do i=1,novlim(klenv,3)
c$$$         if (bruntV2(i).gt.BVmax) then
c$$$            BVmax = bruntV2(i)
c$$$         endif
c$$$      enddo
c$$$
c$$$      i=0
c$$$      rmax = 0.d0
c$$$      rmin = r(nmod)
c$$$      do while (r(i).lt.0.10*r(nmod))
c$$$         if (bruntV2(i).ge.0.1d0*BVmax) then
c$$$            if (r(i).gt.rmax) rmax = r(i)
c$$$            if (r(i).lt.rmax) rmin = r(i)
c$$$            ip = i+1
c$$$            omcore=omcore+(rmid(ip)-rmid(i))*omega(i)
c$$$         endif      
c$$$         i = i+1
c$$$      enddo
c$$$      rcav = rmax-rmin
c$$$      omegamax = omcore/rcav
c$$$         
c$$$!!!
c$$$c..   evolvar13
c$$$      Nmax = 1.d-50
c$$$      omegamax = 1.d-50
c$$$      do i=1,nmod
c$$$         if (i.lt.novlim(klenv,3).and.(brunt(i).gt.Nmax))
c$$$     &        Nmax = brunt(i)
c$$$         if (omega(i).gt.omegamax) omegamax = omega(i)
c$$$      enddo
c$$$      print *,'omegamax2=',omegamax
c$$$
c$$$
c$$$      Nt2=0.d0
c$$$      Nu2=0.d0
c$$$      do i = 1,nmod
c$$$         bruntV2(i) = 0.d0
c$$$         bruntV(i) = 0.d0
c$$$      enddo
c$$$
c$$$      do i=1,nmod-1
c$$$         grav(i) = gmr(i)/gmr(nmod)
c$$$         gs = gmr(nmod)
c$$$
c$$$         Nt2(i) = gs*grav(i)*deltaKS(i)*(abad(i)-abla(i))/(hp(i))
c$$$         Nu2(i) = gs*grav(i)*phiKS(i)*abmu(i)/(hp(i))
c$$$
c$$$         bruntV2(i) = Nt2(i)+Nu2(i)
c$$$      enddo
c$$$      somegam = 0
c$$$      do i=1,nmod
c$$$         if (i.lt.novlim(klenv,3).and.(Nu2(i).gt.Nu2max)) then
c$$$            Nu2max = Nu2(i)
c$$$            nshmax = i
c$$$         endif
c$$$      enddo
c$$$      do i=1,nshmax
c$$$        somegam = somegam + omega(i)*dm(i)
c$$$      enddo
c$$$      omegamax = somegam/m(nshmax)
c$$$      print *,'omegamax1=',omegamax


      

cc Coupling between envelope and radiative core Timescale

      if ((novlim(ienv,3).gt.1).and.(novlim(klenv,3).gt.1)) then
         ndt = novlim(klenv,3)-1
         print *,'ndt=',ndt
         Dshear = 0.d0
         if (omega(ndt)/=omega(ndt-1))
     &        Dshear = 0.5d0*(xnuvv(ndt)+xnuvv(ndt-1))
         jconv = k2conv*omega_S*2.d0/(3.d0*totm*msun*(solr*rsun)**2)
         fcirc = -0.2d0*ro(ndt)*r(ndt)**4*omega(ndt)*urs(ndt)*ur_Ss
         print *,omega(ndt),omega(ndt-1)
         fshear = -ro(ndt)*r(ndt)**4*Dshear*
     &        (omega(ndt)-omega(ndt-1))/(r(ndt)-r(ndt-1))
         tauJ = seci*(k2conv*jrad-k2rad*jconv)/(k2conv+k2rad)
         if (fcirc.ne.0.d0.or.fshear.ne.0.d0) tautot = tauJ/
     &        (fcirc+fshear)
         if (fcirc.ne.0.d0) then
            taucirc = tauJ/fcirc
            if (fshear.ne.0.d0) then
               taushear = tauJ/fshear
            else
               taushear = 0.d0
            endif
         else
            taucirc = 0.d0
            if (fshear.eq.0.d0) then
               tautot = 0.d0
               taushear = 0.d0
            else
               taushear = tauJ/fshear
            endif
         endif
      else
         ndt = nmod
         fcirc = 0.d0
         fshear = 0.d0
         tauJ = 0.d0
         tautot = 0.d0
         taushear = 0.d0
         taucirc = 0.d0
      endif


c      if (irotbin.eq.1) then
c         omegacore = omega(irohalf)
c      else
c         omegacore = 0.d0
c      endif

c.. evolvar13
      write (29,2300) model,omegajrad,jtrue,jobs,jrad,Bequi,omegamax,
     &     Nmax,cseq


c$$$c.. evolvar14
c$$$      write (30,2600) model,omegacore,jtrue,jobs,jrad,Bequi,omegamax,
c$$$     &     Nmax,tautot,taucirc,taushear,cseq

*______________________________________
***   storage of the abundance profiles
*--------------------------------------

c..   evolchc*
      write (21,2200) model,(xsp(1,idxspc(i)),i = 1,11),cseq
      write (22,2200) model,(xsp(1,idxspc(i)),i = 12,22),cseq
      write (123,2200) model,(xsp(1,idxspc(i)),i = 23,33),cseq
      write (124,2200) model,(xsp(1,idxspc(i)),i = 34,44),cseq
c      write (25,2200) model,(xsp(1,idxspci(i)),i = 1,9),cseq
c      write (26) model,time,(xsp(1,i),i = 1,nsp)

      if (nsconv.gt.0) then
         nt = novlim(nsconv,4)
      else
         nt = neff
      endif
      if (tau(nt).ge.1.d2) nt = neff
c..   evolchs*
      write (31,2200) model,(xsp(nt,idxsps(i)),i = 1,11),cseq
      write (32,2200) model,(xsp(nt,idxsps(i)),i = 12,22),cseq
      write (33,2200) model,(xsp(nt,idxsps(i)),i = 23,33),cseq
      write (34,2200) model,(xsp(nt,idxsps(i)),i = 34,nprint),cseq

*__________________________________________________________________________
***   Storage of L,R,X and Y with maximum precision for solar calibration
*--------------------------------------------------------------------------

      if (time*seci.eq.4.5700000000d+09.and.totm.eq.1.d0) then
         open(unit=115,file='input_calibration',status='unknown')
         write(115,2400)soll,solr,xsp(nt,idxsps(2))+xsp(nt,idxsps(3)),
     &   xsp(nt,idxsps(4))+xsp(nt,idxsps(5))
         close(115)
      endif

*________________________________________
***   storage of convective turnover-time
*----------------------------------------

!!! Initialisation
      tg = 0.d0
      tc_max = 0.d0
      rc_max= 0.d0
      tc_m = 0.d0
      rc_m = 0.d0
      rc_r = 0.d0
      tc_r = 0.d0
      tc_hp = 0.d0
      rc_hp = 0.d0
      tc = 0.d0
      rc = 0.d0

      tg_core = 0.d0
      tc_max_core = 0.d0
      rc_max_core= 0.d0
      tc_m_core = 0.d0
      rc_m_core = 0.d0
      rc_r_core = 0.d0
      tc_r_core = 0.d0
      tc_hp_core = 0.d0
      rc_hp_core = 0.d0
      tc_core = 0.d0
      rc_core = 0.d0
      omega_core = 0.d0
!!!
      if (icore.eq.1) then
         print *,'convective core !'

C     turnover time at Hp/2 below top of CC
            k = it
            do while (r(k).lt.r(it)-0.5d0*hp(it)
     &           .and.k.gt.1)
               k = k-1
            enddo
            if (k.eq.0) k = 1
            if (k.lt.ib.or.sconv(k).eq.0.d0) then
               tc = 0.d0
               rc = 0.d0
            else
               tc_core = alphac*hp(k)/sconv(k)
               rc_core = r(k)
            endif

c     turnover time at Hp below top of CC
            k = it
            do while (r(k).gt.r(it)-hp(it).and.k.lt.nmod)
               k = k-1
            enddo
            if (k.eq.0) k = 1
            if (k.lt.ib.or.sconv(k).eq.0.d0) then
               tc_hp_core = 0.d0
               rc_hp_core = 0.d0
            else
               tc_hp_core = alphac*hp(k)/sconv(k)
               rc_hp_core = r(k)
            endif

c     turnover time at 1/2R_cc
            rcc = r(it)
            k = it
            do while (r(k).gt.r(it)-0.5d0*rcc.and.
     &           k.gt.ib)
               k = k-1
            enddo
            if (k.eq.0) k = 1
!            if (k.eq.ib) k = ib + 1
            if (sconv(k).eq.0.d0) then
               tc_r_core = 0.d0
               rc_r_core = 0.d0
            else
               tc_r_core = alphac*hp(k)/sconv(k)
               rc_r_core = r(k)
            endif

c     turnover time at 1/2M_cc
            mcc =  m(it)
            k = it
            do while (m(k).gt.0.5*mcc.and.
     &           k.gt.ib)
               k = k-1
            enddo
            if (k.eq.0) k = 1
!            if (k.eq.ib) k = ib + 1
            if (sconv(k).eq.0.d0) then
               tc_m_core = 0.d0
               rc_m_core = 0.d0
            else
               tc_m_core = alphac*hp(k)/sconv(k)
               rc_m_core = r(k)
            endif

c     turnover time with maximal value (bottom and top layers are not taken into account)
            tc_max_core = -1.d0
            do k = ib+1,it-1
               if (sconv(k).gt.0.d0.and.
     &              tc_max_core.lt.alphac*hp(k)/sconv(k)) then
                  tc_max_core = alphac*hp(k)/sconv(k)
                  rc_max_core = r(k)
               endif
            enddo
            if (rc_max_core.lt.1.d-30) rc_max = 1.d-30
            if (tc_max_core.eq.-1.d0) then
               tc_max_core = 0.d0
               rc_max_core = 0.d0
            endif

c     integreted turmover time (bottom and top layers are not taken into account)
            tg_core = 0.d0
            do k = ib+1,it-1
               if (sconv(k).gt.0.d0) then
                  tg_core = tg_core+(r(k+1)-r(k))/sconv(k)
               endif
            enddo
            if (irotbin.eq.1) omega_core = omega(it-1)
         endif

!!!! Envelope -----------------------------------------------
      if (ienv.gt.0.or.novlim(klenv,3)-1.le.1) then
         print *,'convective envelope !'
!         if (novlim(ienv,3).ne.0) then
            kbce = novlim(ienv,3)
            hpbce = hp(kbce)

            if (0.05d0*hp(novlim(ienv,3)).lt.r(nmod)) then
               do while (r(kbce).lt.r(novlim(ienv,3))+
     &              0.05d0*hp(novlim(ienv,3)).and.kbce.lt.nmod)
                  kbce = kbce+1
               enddo
            endif
            rbce = r(kbce)

            ktce = novlim(ienv,4)
            hptce = hp(ktce)
            if (nphase.gt.2) then
               factdc = 0.05d0
            else
               factdc = 0.1d0
            endif
            do while (r(ktce).gt.r(novlim(ienv,4))-
     &           factdc*hp(novlim(ienv,4))
     &           .and.kbce.gt.1)
               ktce = ktce-1
            enddo

C     turnover time at Hp/2 over bottom of CE
            k = novlim(ienv,3)
            do while (r(k).lt.r(novlim(ienv,3))+0.5d0*hp(novlim(ienv,3))
     &           .and.k.lt.nmod)
               k = k+1
            enddo
            if (k.eq.0) k = 1
            if (k.ge.novlim(ienv,4).or.sconv(k).eq.0.d0) then
               tc = 0.d0
               rc = 0.d0
            else
               tc = alphac*hp(k)/sconv(k)
               rc = r(k)
            endif

c     turnover time at Hp over bottom of CE
            k = novlim(ienv,3)
            do while (r(k).lt.r(novlim(ienv,3))+hp(novlim(ienv,3))
     &           .and.k.lt.nmod)
               k = k+1
            enddo
           if (k.eq.0) then
               k = 1
            else if (k.ge.novlim(ienv,4).or.sconv(k).eq.0.d0) then
               tc_hp = 0.d0
               rc_Hp = 0.d0
            else
               tc_hp = alphac*hp(k)/sconv(k)
               rc_hp = r(k)
            endif

c     turnover time at 1/2R_ce
            rce = r(novlim(ienv,4))-r(novlim(ienv,3))
            k = novlim(ienv,3)
            do while (r(k).lt.r(novlim(ienv,3))+0.5d0*rce.and.
     &           k.lt.novlim(klenv,4))
               k = k+1
            enddo
            if (k.eq.0) k = 1
            if (k.eq.novlim(ienv,4)) k = novlim(klenv,4)-1
            if (sconv(k).eq.0.d0) then
               tc_r = 0.d0
               rc_r = 0.d0
            else
               tc_r = alphac*hp(k)/sconv(k)
               rc_r = r(k)
            endif

c     turnover time at 1/2M_ce
            mce =  m(novlim(ienv,4))-m(novlim(ienv,3))
            k = novlim(ienv,3)
            do while (m(k).lt.m(novlim(ienv,3))+0.5d0*mce.and.
     &           k.lt.novlim(klenv,4))
               k = k+1
            enddo
            if (k.eq.0) k = 1
            if (k.eq.novlim(ienv,4)) k = novlim(klenv,4)-1
            if (sconv(k).eq.0.d0) then
               tc_m = 0.d0
               rc_m = 0.d0
            else
               tc_m = alphac*hp(k)/sconv(k)
               rc_m = r(k)
            endif

c     turnover time with maximal value (bottom and top layers are not taken into account)
            tc_max = -1.d0
            do k = kbce,ktce
               if (sconv(k).gt.0.d0.and.
     &              tc_max.lt.alphac*hp(k)/sconv(k)) then
                  tc_max = alphac*hp(k)/sconv(k)
                  rc_max = r(k)
               endif
            enddo
            if (rc_max.lt.1.d-30) rc_max = 1.d-30
            if (tc_max.eq.-1.d0) then
               tc_max = 0.d0
               rc_max = 0.d0
            endif


c     integrated turmover time (bottom and top layers are not taken into account)
            tg = 0.d0
            do k = kbce,ktce
               if (sconv(k).gt.0.d0) tg = tg+(r(k+1)-r(k))/sconv(k)
            enddo

!            write (41,2500) model,tc/secd,rc/rsun,tc_hp/secd,
!     &           rc_hp/rsun,tc_r/secd,rc_r/rsun,tauc,rc_m/rsun,
!     &           tc_max/secd,rc_max/rsun,tg/secd,cseq
!         else
!            write (41,2500) model,0.d0,0.d0,0.d0,
!     &           0.d0,0.d0,0.d0,0.d0,0.d0,
!     &           0.d0,0.d0,0.d0,cseq
!         endif
         endif
c..   evoltc1 : convective turnover in the envelope
            write (41,2500) model,tc*seci,rc/rsun,tc_hp*seci,
     &           rc_hp/rsun,tc_r*seci,rc_r/rsun,tauc,rc_m/rsun,
     &           tc_max*seci,rc_max/rsun,tg*seci,cseq
c..   evoltc2 : convective turnover in the core
            write (42,2510) model,tc_core*seci,rc_core/rsun,
     &           tc_hp_core*seci,rc_hp_core/rsun,tc_r_core*seci,
     &           rc_r_core/rsun,tauc,rc_m_core/rsun,tc_max_core*seci,
     &           rc_max_core/rsun,tg_core*seci,omega_core,cseq
!

*______________________________________
***   storage of asteroseimic relations
*--------------------------------------

*--------
* fréquence de brunt Vaisalla
*--------

      Nt2=0.d0
      Nu2=0.d0
      do i = 1,nmod
         bruntV2(i) = 0.d0
         bruntV(i) = 0.d0
      enddo

      do i=1,nmod-1
         grav(i) = gmr(i)/gmr(nmod)
         gs = gmr(nmod)

         Nt2(i) = gs*grav(i)*deltaKS(i)*(abad(i)-abla(i))/(hp(i))
         Nu2(i) = gs*grav(i)*phiKS(i)*abmu(i)/(hp(i))

         bruntV2(i) = Nt2(i)+Nu2(i)
      enddo

      do i=2,nmod-1	
         if (bruntV2(i).gt.0.d0.and.crz(i).ne.-2) then
            bruntV(i) = dsqrt(bruntV2(i))

         elseif (bruntV2(i).le.0.d0.and.bruntV2(i-1).gt.0.d0
     &           .and.bruntV2(i+1).gt.0.d0
     &           .and.crz(i).ne.-2) then
            bruntV(i) = dsqrt(-bruntV2(i))
            bruntV2(i) = -bruntV2(i)
         else  		
            bruntV(i) = 0.d0
            bruntV2(i) = 0.d0
         endif
      enddo

c.... grande séparation avec relation asymptotique
        sum = 0.d0
      do k = 1,nmod-1
         sum = sum+(r(k+1)-r(k))/cs(k)
      enddo
      Dnu = 1.d0/(2.d0*sum)

c.... grande séparation avec les relations d'échelle
c     Dnu_sun = 134.9d-6 microHz
      Dnuech = 134.9d-6*(m(nmod)/msun)**0.5d0*(reff/rsun)**(-1.5d0)

c.... frenquency at which oscilation modes are maximal
c     nu_max_sun = 3150.d-6 microHz
      nu_max = 3150.d-6*(m(nmod)/msun)*(reff/rsun)**(-2.d0)*
     &  (teff/5780)**(-0.5d0)

c.... difference entre les relations asymptotiques et d'échelle
      errr = (Dnu-Dnuech)/Dnu

c....rayon acoustique totale
      T_tot=1/(2*Dnu)

c... rayon acoustique à la base de l'env.conv
      sum_BCE = 0.d0
      if (ienv.gt.0) then
         BCE = novlim(ienv,3)
         do k = 1,BCE
            sum_BCE=sum_BCE+(r(k+1)-r(k))/cs(k)
         enddo
      endif
      t_BCE = sum_BCE

c... Rayon acoustique a la base de la zone d'ionisation de l'He++
      sum_He = 0.d0
      do k= 1,nmod-1
	 if (gamma1(k).lt.1.55.and.t(k).le.1e5.and.
     &        gamma1(k).lt.gamma1(k+1)) then
            k_He=k
            goto 23
         endif
      enddo
 23   do k = 1,k_He
         sum_He=sum_He+(r(k+1)-r(k))/cs(k)
      enddo
	
      r_He=r(k_He)
      t_He=sum_He

c.... séparation en période entre les modes g pour l=1
      sum_Dp=0.d0
      Dp=0.d0
      kconv = 2
      do k = 2,nmod
         if (Dconv(k-1).eq.0.d0.and.Dconv(k).gt.0.d0) then
            kconv=k
            goto 25
         endif
      enddo
 25   if (ienv.gt.0) then
         if (kconv.eq.novlim(ienv,3)) then
            kconv=2
         endif
      endif
	
      pass = 0
	
 26   do k=kconv,nmod-1
         if ((2.d0*pi*nu_max)**(2.d0).lt.bruntV2(k).and.
     &        (2.d0*pi*nu_max)**(2.d0).lt.(2.d0*(cs(k)/r(k))**(2.d0)))
     &        sum_Dp=sum_Dp+((r(k+1)-r(k))*bruntV(k)*(r(k))**(-1.d0))
      enddo
      if (sum_Dp.lt.1.d-3.and.pass.eq.0.and.kconv.ne.2) then
         sum_Dp = 0.d0
         kconv = 2
         pass = 1
         goto 26
      endif	
	
      if (sum_Dp.gt.0.d0) then
         Dp=((2.0d0)**(0.5)*pi**(2))/(sum_Dp)
      else
         Dp=-1.0d0
      endif

c...evolas
      write (43,2600) model,Dnu,Dnuech,errr,T_tot,t_BCE,t_He,nu_max,Dp
     &     ,cseq


 100  format (' AGB network : version ',0pf4.2,/,1x,
     &     'General informations:',//,
     &     1x,'model',1x,'phase',4x,'L',6x,'Reff',5x,'R',6x,'Teff',3x,
     &     'rhoeff',6x,'geff',6x,'mlos',10x,'M',9x,'dt',11x,'t',10x,
     &     'it',1x,'crash',1x,'nshell',1x,'cputime',/,132('-'),/)
 110  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Values at: center',//,2x,'model',7x,'Tc',6x,'Tmax',6x,
     &     'Tmax_Mr',4x,'rhoc',5x,'rho_max',4x,'Pc',6x,'betac',4x,
     &     'etac',4x,'degpec',5x,'Enu_Pla',6x,'Enucl',6x,
     &     'Egrav',/,130('-'),/)
 120  format (' AGB network : version ',0pf4.2,/,1x,'Energetics :',
     &     1x,'integrated values',//,
     &     2x,'model',7x,'LH',9x,'LHe',9x,
     &     'LC',9x,'LNe',9x,'LO',10x,'LSi',7x,'Lnu_Pla',7x,
     &     'Lnucl',7x,'Lgrav',3x,'irradmax',/,125('-'),/)
 130  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Convective zones I:',//,
     &     2x,'model',5x,'conv1_Mb',2x,'conv1_Rb',1x,'conv1_Tb',1x,
     &     'conv1_rob',1x,'conv1_Mt',2x,'conv1_Rt',1x,'conv1_Tt',1x,
     &     'conv1_rot',4x,'env_Mb',4x,'env_Rb',4x,'env_Tb',4x,'env_rob',
     &     /,128('-'),/)
 140  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Convective zones II:',//,
     &     2x,'model',5x,'conv2_Mb',3x,'conv2_Mt',3x,'conv3_Mb',3x,
     &     'conv3_Mt',3x,'conv4_Mb',3x,'conv4_Mt',3x,'conv5_Mb',3x,
     &     'conv5_Mt',3x,'thermoh_Mb',3x,'thermoh_Mt',2x,'nconvt',/,
c     &     'conv5_Mt',3x,'conv6_Mb',3x,'conv6_Mt',2x,'nconvt',/,
     &     127('-'),/)
 150  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Burning zones : H-burning',//
     &     2x,'model',4x,'Hburn_Mb',2x,'Hburn_Rb',2x,'Hburn_Tb',1x,
     &     'Hburn_rob',2x,'Hburn_Mt',2x,'Hburn_Rt',2x,'Hburn_Tt',1x,
     &     'Hburn_rot',2x,'Hburn_Mm',3x,'Hburn_em',5x,'Lpp',/,
     &     121('-'),/)
 160  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Burning zones : He-burning',//
     &     2x,'model',3x,'Heburn_Mb',1x,'Heburn_Rb',1x,'Heburn_Tb',1x,
     &     'Heburn_rob',1x,'Heburn_Mt',1x,'Heburn_Rt',1x,'Heburn_Tt',1x,
     &     'Heburn_rot',1x,'Heburn_Mm',1x,'Heburn_em',1x,'Heburn_enum',
     &     1x,'Heburn_etam',/,135('-'),/)
 162  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Burning zones : C-burning',//
     &     2x,'model',4x,'Cburn_Mb',2x,'Cburn_Rb',2x,'Cburn_Tb',1x,
     &     'Cburn_rob',2x,'Cburn_Mt',2x,'Cburn_Rt',2x,'Cburn_Tt',1x,
     &     'Cburn_rot',2x,'Cburn_Mm',2x,'Cburn_em',2x,'Cburn_enum',
     &     1x,'Cburn_etam',/,132('-'),/)
 164  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Burning zones : Ne-burning',//,
     &     2x,'model',2x,'Neburn_Mb',1x,'Neburn_Rb',1x,'Neburn_Tb',1x,
     &     'Neburn_rob',1x,'Neburn_Mt',1x,'Neburn_Rt',1x,'Neburn_Tt',1x,
     &     'Neburn_rot',1x,'Neburn_Mm',1x,'Neburn_em',1x,'Neburn_enum',
     &     1x,'Neburn_etam',/,134('-'),/)
 170  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Rotation and Accretion:',
     &     //,2x,'model',5x,'macc',5x,'Racc',5x,'Macc',4x,
     &     'Lshear',4x,'k2conv',3x,'k2rad',5x,'Jconv',5x,'Jrad',
     &     5x,'vsurf',7x,'Usurf',7x,'Fenerg',5x,'dmoment'/,117('-'),/)
 180  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Convective zones II:',//,
     &     2x,'model',5x,'conv2_Rb',4x,'conv2_Rt',4x,'conv3_Rb',4x,
     &     'conv3_Rt',4x,'conv4_Rb',4x,'conv4_Rt',4x,'conv5_Rb',4x,
     &     'conv5_Rt',4x,'therm_Rb',4x,'therm_Rt',/,129('-'),/)
c     &     'conv5_Rt',4x,'conv6_Rb',4x,'conv6_Rt',/,129('-'),/)
 190  format (' AGB network : version ',0pf4.2,/,1x,a7,
     &     ' chemical abundances:',
     &     //,2x,'model',6x,'n ',7x,' H1 ',7x,' H2 ',7x,' He3',7x,
     &     ' He4',7x,' Li6', 7x,' Li7',7x,' Be7',7x,' Be9',7x,' B10',
     &     7x,' B11', /,127('-'),/)
 200  format (' AGB network : version ',0pf4.2,/,1x,a7,
     &     ' chemical abundances:',
     &     //,2x,'model',4x,' C12',7x,' C13',7x,' C14',7x,' N14',7x,
     &     ' N15',7x,' O15', 7x,' O16',7x,' O17',7x,' O18',7x,' F19',
     &      7x, 'Ne20',/,127('-'),/)
 210  format (' AGB network : version ',0pf4.2,/,1x,a7,
     &     ' chemical abundances:',
     &     //,2x,'model',4x,'Ne21',7x,'Ne22',7x,'Na23',7x,'Mg24',7x,
     &     'Mg25',7x,'Mg26', 7x,'Al26m',6x,'Al26g',6x,'Al27',7x,'Si28',
     &     7x,'Si29',/,127('-'),/)
 220  format (' AGB network : version ',0pf4.2,/,1x,a7,
     &     ' chemical abundances:',
     &     //,2x,'model',4x, 'Si30',7x,' P31',7x,' S32',7x,' S33',7x,
     &     ' S34',7x,' S35',7x,'CL35',7x,' S36',7x,'CL36',7x,'CL37',6x,
     &     'heavy',/,127('-'),/)
 230  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Convective envelope turnover times:',//,
     &     2x,'model',4x,'tc',9x,'rc',9x,'tc_hp',6x,'rc_hp',6x,'tc_r',
     &     7x,'rc_r',7x,'tc_m',6x,'rc_m',6x,'tc_max',5x,'rc_max',5x,
     &     'tg',/,117('-'),/)
 240  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Convective core turnover times:',//,
     &     2x,'model',3x,'tc_cc',5x,'rc_cc',5x,'tc_hp_cc',3x,
     &     'rc_hp_cc',3x,'tc_r_cc',3x,'rc_r_cc',4x,'tc_m_cc',5x,
     &     'rc_m_cc',3x,'tc_max_cc',2x,'rc_max_cc',3x,'tg_cc',5x,
     &     'omegac',/,117('-'),/)
 250  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Asteroseimic relations:',//,
     &     2x,'model',5x,'Dnu',8x,'Dnu_ech',3x,'Dnu_error',5x,'Ttot',
     &     8x,'tbce',8x,'tHe',8x,'numax',7x,'Dpg ',/,117('-'),/)

 260  format (' AGB network : version ',0pf4.2,/,1x,
     &     'Rotation and Accretion (bis):',
     &     //,2x,'model',5x,'omegacore',4x,'Jtrue',5x,'Jobs',7x,'Jrad',
     &     5x,'Bequi',3x,'omegamax',3x,'Nmax',4x,'tautot',3x,'taucirc',
     &     2x,'taushear',/,117('-'),/)

 270  format (' AGB network : version ',0pf4.2,/,1x,
     &     'F. Gallet data :',
     &     //,2x,'model',5x,'L/Lsun',4x,'Teff',5x,'time(yr)',7x,'R/Rsun',
     &     5x,'k2conv',3x,'k2rad',3x,'Mrad/Msun',4x,'Rrad/Rsun',/,
     &     117('-'),/)

 1000 format (i8,1x,i1,1x,0pf11.6,1x,0pf7.5,1x,0pf7.5,1x,0pf7.0,1x,
     &     1pe9.3,1x,1pe9.2,1x,1pe10.4,1x,0pf12.8,1x,1pe9.3,1x,
     &     1pe16.10,1x,i2,1x,i2,1x,i4,1x,0pe8.2,' ####',a2)
 1100 format (i8,1x,i1,1x,0pf11.2,1x,0pf7.2,1x,0pf7.2,1x,0pf7.0,1x,
     &     1pe9.3,1x,1pe9.2,1x,1pe10.4,1x,0pf12.8,1x,1pe9.3,1x,
     &     1pe16.10,1x,i2,1x,i2,1x,i4,1x,0pe8.2,' ####',a2)
 1150 format (i8,1x,i1,1x,0pf11.2,1x,0pf7.2,1x,0pf7.0,1x,0pf7.0,1x,
     &     1pe9.3,1x,1pe9.2,1x,1pe10.4,1x,0pf12.8,1x,1pe9.3,1x,
     &     1pe16.10,1x,i2,1x,i2,1x,i4,1x,0pe8.2,' ####',a2)
 1200 format (i8,1x,i1,1x,0pf11.6,1x,0pf7.5,1x,0pf7.5,1x,0pf7.0,1x,
     &     1pe9.3,1x,1pe9.2,1x,1pe10.4,1x,0pf12.8,1x,1pe9.3,1x,
     &     1pe16.10,1x,i2,1x,i2,1x,i4,1x,0pe8.2)
 1300 format (i8,1x,i1,1x,0pf11.2,1x,0pf7.2,1x,0pf7.2,1x,0pf7.0,1x,
     &     1pe9.3,1x,1pe9.2,1x,1pe10.4,1x,0pf12.8,1x,1pe9.3,1x,
     &     1pe16.10,1x,i2,1x,i2,1x,i4,1x,0pe8.2)
 1350 format (i8,1x,i1,1x,0pf11.2,1x,0pf7.2,1x,0pf7.0,1x,0pf7.0,1x,
     &     1pe9.3,1x,1pe9.2,1x,1pe10.4,1x,0pf12.8,1x,1pe9.3,1x,
     &     1pe16.10,1x,i2,1x,i2,1x,i4,1x,0pe8.2)
 1400 format (i8,2x,2(1x,1pe9.3),1x,0pf9.6,3(1x,1pe9.3),1x,
     &     0pf6.4,1x,0pf8.3,1x,0pf8.3,3x,1pe10.4,2(1x,1pe11.4),
     &     1x,a2)
 1410 format (i8,2x,2(1x,1pe9.3),1x,0pf9.6,3(1x,1pe9.3),1x,
     &     0pf6.4,1x,0pf8.1,1x,0pf8.3,3x,1pe10.4,2(1x,1pe11.4),
     &     1x,a2)
 1450 format (1x,'dEpot =',1pe10.3,', dEkin =',1pe10.3,', dEint =',
     &     1pe10.3,', dEphot =',1pe10.3,', dEnupla =',1pe10.3,
     &     ', dEnuc =',1pe10.3,', equil =',1pe9.3)
 1500 format (i8,2x,2(1x,1pe10.4),7(1x,1pe11.4),1x,1pe11.4,a2)
 1600 format (i8,2x,2(1x,0pf9.6,1x,1pe9.3,1x,1pe9.3,1x,1pe9.3),1x,
     &     0pf9.6,1x,1pe11.4,1x,0pf8.4,a2)
 1700 format (i8,2x,2(1x,0pf9.6,1x,1pe9.3,1x,1pe9.3,1x,1pe9.3),1x,
     &     0pf9.6,1x,1pe10.3,1x,1pe10.3,0pf8.3,a2)
 1900 format (i8,2x,1x,0pf10.7,1x,1pe9.3,1x,0pf7.4,1x,0pf8.4,1x,
     &     f10.6,1x,1pe9.3,1x,0pf7.4,1x,0pf8.4,2x,1x,0pf10.6,1x,
     &     1pe9.3,1x,0pf9.5,1x,f9.5,a2)
 2000 format (i8,2x,10(1x,0pf10.6),2x,1x,i2,a1,a2)
 2050 format (i8,2x,10(1x,1pe11.5),a2)
 2100 format (i8,2x,1x,1pe8.2,1x,1pe9.3,1x,0pf6.4,1x,1pe10.4,2(1x,
     &     1pe10.4),1x,1pe11.4,1x,1pe8.2,1x,1pe11.4,1x,1pe10.4,1x,
     &     1pe10.4,1x,1pe11.4,a2)
 2200 format (i8,11(1x,1pe10.4),a2)

 2300 format (i8,2x,1pe8.2,6(1x,1pe8.2),a2)

 2400 format(4(1x,0pf20.12))
 1211 format (i8,1x,22(e15.4E4,1x))
 2500 format(i8,11(1x,1pe10.4),a2)
 2510 format(i8,12(1x,1pe10.4),a2)
 2600 format (i8,1x,1pe10.4,1x,1pe11.4,1x,1pe11.4,1x,1pe11.4,
     &   1x,1pe11.4,1x,1pe11.4,1x,1pe11.4,1x,1pe11.4,a2)
      return
      end
