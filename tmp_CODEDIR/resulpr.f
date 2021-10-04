      SUBROUTINE resulpr (error)

************************************************************************
* Save computed model and write display file                           *
* 18/06 : adding some writings into evoldisp in case of rotation       *
* Modifs CC ondes (2/04/07 -> 23/11/07)                                *
* $LastChangedDate:: 2016-05-11 17:20:46 +0200 (Mer, 11 mai 2016)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 62                                                          $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.grad'
      include 'evolcom.ion'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.opa'
      include 'evolcom.rot'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.transp'
      include 'evolcom.var'
C Modifs CC ondes (2/04/07)
      include 'evolcom.igw'
c      include 'evolcom.transpondesexcit'

      integer error
      integer imin,imax
      integer j,k,kl,l,kl1
C Modifs CC ondes (2/04/07) -->
      integer kk
      double precision lgdmic,lgdturb
      double precision dturb,dmic,vdmic
C <--

      double precision nturb
      double precision Dhold
      double precision dift,dife
      double precision vom,vrray,vxpsi
      double precision abmurj,abmuj
      double precision FeH
      double precision tautime,Ereac,taureac
      double precision Ushear(nsh),ULj(nsh)
      double precision brunt,lumwaves(nsh),bruntV2,bruntV

      double precision depottotn,depottot_surfn,depottot_coren
      double precision depotondestotn,vrn

      double precision abundelec     ! Ajout pour creer donnees Montreal (Mars 2018)
      
      character*1 crzc(nsh)
      logical partialmix

      common /calcDh/ Dhold(nsh)
      common /difcirc/ dift(nsh),dife(nsh)
      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)
      common /gradmuj/ abmurj(nsh),abmuj(nsh)
      common /metal/ FeH

      common /brunt/ brunt(nsh),bruntV2(nsh),bruntV(nsh)
      common /abundelec/ abundelec(nsh)       ! Ajout pour creer donnees Montreal (Mars 2018)
      
C Modifs CC ondes (2/04/07) -->
      common /coefdiffpaqhe/ xvdiffhe(nsh),xvdiffhenoc(nsh)
      common /coefdiff/ dturb(nsh),dmic(nsh),vdmic(nsh),lgdturb(nsh)
     &     ,lgdmic(nsh)
      double precision xvdiffhe,xvdiffhenoc
C     <--
c Ajout TD (18/01/2019)
      common /diffutest/ zmeanO16(nsh),dens_elec(nsh) ,coulomb(nsh)
     $     ,Zpaq(nsh),Zpaq1(nsh),Zpaq2(nsh),Kpaq(nsh)
      double precision zmeanO16,dens_elec,coulomb
     $     ,Zpaq,Zpaq1,Zpaq2,Kpaq
c Fin ajout      
c Ajout par TD (22/02/2018) + modif 25/04/2019    
      common /coefdiffpaqli/ xvdiffli(nsh),xvdli(nsh),
     &     xvthermli(nsh),xdmicli(nsh)
      double precision xvdiffli,xvdli,xvthermli,xdmicli
      common /coefdiffpaqc12/ xvdiffc12(nsh),xvdiffc12noc(nsh)
     $     ,xvdiffc13(nsh),xvdiffc13noc(nsh)
      double precision xvdiffc12,xvdiffc12noc
     $     ,xvdiffc13,xvdiffc13noc
      common /coefdiffpaqn14/ xvdiffn14(nsh),xvdiffn14noc(nsh)
     $     ,xvdiffn15(nsh),xvdiffn15noc(nsh)
      double precision xvdiffn14,xvdiffn14noc
     $     ,xvdiffn15,xvdiffn15noc      
      common /coefdiffpaqo16/ xvdiffo16(nsh),xvdiffo16noc(nsh)
     $     ,xvdiffo17(nsh),xvdiffo17noc(nsh),xvdiffo18(nsh)
     $     ,xvdiffo18noc(nsh)
      double precision xvdiffo16,xvdiffo16noc
     $     ,xvdiffo17,xvdiffo17noc,xvdiffo18,xvdiffo18noc
      common /coefdiffpaqne20/ xvdiffne20(nsh),xvdiffne20noc(nsh)
      double precision xvdiffne20,xvdiffne20noc
      common /coefdiffpaqna23/ xvdiffna23(nsh),xvdiffna23noc(nsh)
      double precision xvdiffna23,xvdiffna23noc
c$$$      common /coefdiffpaqna23/ xvdiffna23(nsh),xvdna23(nsh),   
c$$$     &     xvthermna23(nsh),xdmicna23(nsh),xvdiffna23noc(nsh)
c$$$      double precision xvdiffna23,xvdna23,xvthermna23,xdmicna23
c$$$     $     ,xvdiffna23noc
c Fin ajout
      
      common /dturbulente/ nturb(0:nsh)
      common /nucleaire/ tautime(nsh,nreac),Ereac(nsh,nreac),taureac(nsh
     $     ,nreac)
      
      dimension depottotn(nsh),depottot_surfn(nsh),depottot_coren(nsh)
      dimension depotondestotn(nsh),vrn(nsh)

      if (imodpr.eq.11.or.imodpr.eq.12)
     &     call myflamespeed (.false.,imodpr)

c$$$      do k=1,nsp
c$$$         do j = 1,nmod
c$$$            write(666,'(1x,i4,1x,i4,9(1x,1pe11.4))') k,j,t(j)
c$$$     $           ,zmean(j,k)
c$$$         enddo
c$$$      enddo
      
*________________________________________________________
***   storage of the calculated models in the binary file
*--------------------------------------------------------
      print *,'In resulpr xsp 1',xsp(1420,2),xsp(1450,2)

      if (error.eq.-30) neff = -neff

      if (nsconv.gt.0.and.turnover(1).eq.0.d0) then
c..   determine turnover timescale
         do kl = 1,nsconv
            imin = novlim(kl,3)
            imax = novlim(kl,4)
            call mix (dm,imin,imax,kl,5,partialmix)
         enddo
      endif
      
      do k = 1,nmod
         if (crz(k).eq.-2) crzc(k) = 'c'
         if (crz(k).eq.-1) crzc(k) = 'a'
         if (crz(k).eq.-3) crzc(k) = 's'
         if (crz(k).eq.3) crzc(k) = 'S'
         if (crz(k).eq.2) crzc(k) = 't'
         if (crz(k).eq.1) crzc(k) = 'o'
         if (crz(k).eq.4) crzc(k) = 'r'
      enddo

      call system ('cp -f nextini.bin oldini.bin')
      rewind (92)
      write (92) model,nphase,totm,time,dtn,neff,code_version,
     &     nmod,nis-1,mdredgeup,FeH,turnenv,rayon0,
     &     (crzc(k),u(k),r(k),lnf(k),lnt(k),lum(k),
     &     vu(k),vr(k),vlnf(k),vlnt(k),vlum(k),facc(k),
     &     exposure(k),mueinv(k),dm(k),vhconv(k),omega(k),
     &     ro(k),vro(k),p(k),vp(k),e(k),ve(k),enucl(k),s(k),
     &     vs(k),sconv(k),pvisc(k),tau(k),
     &     (xsp(k,l),l = 1,nis-1),k = 1,nmod),
     &     (mtotlos(j),j = 1,nis-1)
      close (92)

      print *,'In resulpr xsp 2',xsp(1420,2),xsp(1450,2)

      do k = 1,nmod
         venucl(k) = enucl(k)
         vsconv(k) = sconv(k)
         vpvisc(k) = pvisc(k)
         vmueinv(k) = mueinv(k)
      enddo
      if (irotbin.eq.1) then
c         call system ('cp -f nextang.bin modang.bin')
         rewind (94)
c         rewind(83)
         write (94) ndtold,ndbold,ur_Ss,xmom_tots,inits,
     &        (auxs(k),urs(k),vxpsi(k),xpsis(k),V_circ(k),
     &        xalpha(k),xlambdas(k),xKt(k),dift(k),dife(k),
     &        Dhold(k),Dconv(k),abmurj(k),xnuvv(k),xnum(k),
     &        k = 1,nmod),
     &        ((vxsp(k,l),l = 1,nis-1),k = 1,nmod)
C Modifs CC ondes (2/04/07) -->
c         print *,'urs=',urs(1:nmod)
         do k=1,nmod
         kk=nmod-k+1
            if (idiffcc) then
                lgdturb(k) = 0.d0
                lgdmic(k) = 0.d0
                if (dturb(k).gt.0.d0) lgdturb(k) = log10(dturb(k))
c                if (xdmiche(k).gt.0.d0) lgdmic(k) = log10(xdmiche(k))  ! modif Mai 2019
            else
                lgdturb(k) = 0.d0
                lgdmic(k) = 0.d0
            endif
         enddo
C <--
      endif
      open (unit = 92,file = 'nextini.bin',form = 'unformatted',
     &     status = 'unknown')


c      if (microdiffus.or.thermohaline.or.igw) then
cc         write (81,8100)
c         rewind(85)
c         do k = 1,nmod
c            if (dabs(Ereac(k,ibe7pg)).lt.1d-99) Ereac(k,ibe7pg) = 0.d0
c            if (dabs(lumondestot(k)).lt.1d-99) lumondestot(k) = 0.d0
c            if (dabs(lumondes_surf(k)).lt.1d-99) then
c               lumondes_surf(k) = 0.d0
c            endif
c            if (dabs(lumondes_core(k)).lt.1d-99) then
c               lumondes_core(k) = 0.d0
c            endif
c         enddo
c         if (igw) then
c            do k = 1,nmod
c               depottotn(k) = depottot(k)/lsun
c               depottot_surfn(k) = depottot_surf(k)/lsun
c               depottot_coren(k) = depottot_core(k)/lsun
c               depotondestotn(k) = depotondestot(k)/lsun
c               vrn(k) = vr(k)/rsun
c            enddo
c            write (85) (lgdmic(k),lgdturb(k),vom(k),
c     &           depottotn(k)/lsun,depottot_surfn(k)/lsun,
c     &           depottot_coren(k)/lsun,
c     &           depotondestotn(k)/lsun,Dondes(k),Dondeschim(k),
c     &           lumwave(k),lumwave_core(k),
c     &           lumondes_surf(k),lumondes_core(k),lumondestot(k),
c     &           brunt_o(k),vrn(k)/rsun,vxpsi(k),Fenerg,
c     &           Fenerg_core,
c     &           Dmicro(k),vmicro(k),Dthc(k),phiKS(k),
c     &           deltaKS(k),tautime(k,ippg),Ereac(k,ippg),tautime(k
c     &           ,ipdg), Ereac(k,ipdg),tautime(k,i2he3),Ereac(k
c     &           ,i2he3),tautime(k,ihe3ag), Ereac(k,ihe3ag),tautime(k
c     &           ,ibe7beta),Ereac(k,ibe7beta), tautime(k,ili7pa)
c     &           ,Ereac(k,ili7pa),tautime(k,ibe7pg), Ereac(k,ibe7pg)
c     &           ,tautime(k,ib8beta),Ereac(k,ib8beta), tautime(k
c     &           ,ic13pg),Ereac(k,ic13pg),tautime(k,in14pg), Ereac(k
c     &           ,in14pg),tautime(k,icpg),Ereac(k,icpg), k = 1,nmod)
c         else
c            write (85) (Dmicro(k),vmicro(k),Dthc(k),phiKS(k),
c     &           deltaKS(k),tautime(k,ippg),Ereac(k,ippg),
c     &           tautime(k,ipdg), Ereac(k,ipdg),tautime(k,i2he3),
c     &           Ereac(k,i2he3),tautime(k,ihe3ag),Ereac(k,ihe3ag),
c     &           tautime(k,ibe7beta),Ereac(k,ibe7beta),
c     &           tautime(k,ili7pa),Ereac(k,ili7pa),tautime(k,ibe7pg),
c     &           Ereac(k,ibe7pg),tautime(k,ib8beta),Ereac(k,ib8beta),
c     &           tautime(k,ic13pg),Ereac(k,ic13pg),tautime(k,in14pg),
c     &           Ereac(k,in14pg),tautime(k,icpg),Ereac(k,icpg),
c     &           k = 1,nmod)
c         endif
c      endif


*_________________
***   print banner
*-----------------

      write (nout,*)
      if (nsconv.gt.0) then
         write (nout,1100)
      endif
      if (.not.ledouxmlt) then
         do j = 1,nsconv
            novlim(j,1) = novlim(j,3)
            novlim(j,2) = novlim(j,4)
         enddo
         nschwar = nsconv
      endif
      do kl = 1,nschwar
         if (novlim(kl,2)-novlim(kl,1).gt.5) then
            if ((diffover.or.novopt.gt.0).and..not.ledouxmlt.and.
     &           (novlim(kl,8).gt.novlim(kl,2).or.
     &           novlim(kl,1).gt.novlim(kl,7))) then
               write (90,1200) novlim(kl,1),novlim(kl,2),
     &              m(novlim(kl,1))/msun,m(novlim(kl,2))/msun,
     &              turnover(kl)*seci,novlim(kl,7),novlim(kl,8),
     &              (m(novlim(kl,1))-m(novlim(kl,7)))/msun,
     &              (m(novlim(kl,8))-m(novlim(kl,2)))/msun
               write (nout,1200) novlim(kl,1),novlim(kl,2),
     &              m(novlim(kl,1))/msun,m(novlim(kl,2))/msun,
     &              turnover(kl)*seci,novlim(kl,7),novlim(kl,8),
     &              (m(novlim(kl,1))-m(novlim(kl,7)))/msun,
     &              (m(novlim(kl,8))-m(novlim(kl,2)))/msun
            else
               write (90,1210) novlim(kl,1),novlim(kl,2),
     &              m(novlim(kl,1))/msun,m(novlim(kl,2))/msun,
     &              turnover(kl)*seci
               write (nout,1210) novlim(kl,1),novlim(kl,2),
     &              m(novlim(kl,1))/msun,m(novlim(kl,2))/msun,
     &              turnover(kl)*seci
            endif
         endif
      enddo
      if ((diffover.or.novopt.gt.0).and.ledouxmlt) then
         do kl = 1,nledoux
            kl1 = 0
            do l = 1,nschwar
               if (novlim(kl,8).gt.novlim(kl,6).or.
     &              novlim(kl,5).gt.novlim(kl,7)) then
                  kl1 = l
                  exit
               endif
            enddo
            if (kl1.ne.0) then
               write (90,1240) novlim(kl,5),novlim(kl,6),
     &              m(novlim(kl,5))/msun,m(novlim(kl,6))/msun,
     &              turnover(kl)*seci,novlim(kl1,7),novlim(kl1,8),
     &              (m(novlim(kl,5))-m(novlim(kl1,7)))/msun,
     &              (m(novlim(kl1,8))-m(novlim(kl,6)))/msun
               write (nout,1240) novlim(kl,5),novlim(kl,6),
     &              m(novlim(kl,5))/msun,m(novlim(kl,6))/msun,
     &              turnover(kl)*seci,novlim(kl1,7),novlim(kl1,8),
     &              (m(novlim(kl,5))-m(novlim(kl1,7)))/msun,
     &              (m(novlim(kl1,8))-m(novlim(kl,6)))/msun
            else
               write (90,1230) novlim(kl,5),novlim(kl,6),
     &              m(novlim(kl,5))/msun,m(novlim(kl,6))/msun,
     &              turnover(kl)*seci
               write (nout,1230) novlim(kl,5),novlim(kl,6),
     &              m(novlim(kl,5))/msun,m(novlim(kl,6))/msun,
     &              turnover(kl)*seci
            endif
         enddo
      endif
      do kl = 1,nsemiconv
         if (novlim(kl,10)-novlim(kl,9).gt.1.or.
     &        crz(novlim(kl,10)+1).lt.0) then
            write (90,1250) novlim(kl,9),novlim(kl,10),
     &           m(novlim(kl,9))/msun,m(novlim(kl,10))/msun
            write (nout,1250) novlim(kl,9),novlim(kl,10),
     &           m(novlim(kl,9))/msun,m(novlim(kl,10))/msun
         endif
      enddo
      do kl = 1,nthermha
         if (novlim(kl,12)-novlim(kl,11).gt.1.or.
     &        crz(novlim(kl,12)+1).lt.0) then
            write (90,1260) novlim(kl,11),novlim(kl,12),
     &           m(novlim(kl,11))/msun,m(novlim(kl,12))/msun
            write (nout,1260) novlim(kl,11),novlim(kl,12),
     &           m(novlim(kl,11))/msun,m(novlim(kl,12))/msun
         endif
      enddo

c..   Storage of critical Reynolds number, critical Richardson number
c..   and prescription used for horizontal turbulent diffusion coefficient

      if (model.eq.(modeli+1).and.rotation) then
         write (90,1270) Dh_prescr
c        do k = 1,nmod
c           write(667,*),'nturb',k,nturb(k),xnum(k),xnuvv(k)
c        enddo
      endif


      write (90,1000) totm,dms,time*seci,dtn*seci,teff,soll

c     if (no.eq.maxmod.or.time.eq.4.6000000d9*sec) then
      if (no.eq.maxmod) then
         if (microdiffus.or.thermohaline.or.igw) then
            if (igw) then

C Modifs CC ondes (27/10/09)
C depottot_surf(kk) est calcule dans ond_omega
C                   apres appel de ond_int qui calcule depotondes(kk)
C depottot_core(kk) est calcule dans ond_omega
C                   apres appel de ond_int_core qui calcule depotondes(kk)
C depottot(kk)=depottot_surf(kk)+depottot_core(kk) est calcule dans ond_omega
C depotondestot(kk)=depottot(kk)/(rtot**2*omega_S) dans transport_ondes
C depotwaves(k)=depotondestot(nmod-kk+1) dans transport_ondes

               do k = 1,nmod
                  if (dabs(lumondestot(k)).lt.1d-99)
     &                 lumondestot(k) = 0.d0
                  if (dabs(lumondes_surf(k)).lt.1d-99)
     &                 lumondes_surf(k) = 0.d0
                  if (dabs(lumondes_core(k)).lt.1d-99)
     &                 lumondes_core(k) = 0.d0
                  brunt(k)=brunt_o(nmod-k+1)
                  lumwaves(k)=lumwave(nmod-k+1)
               enddo
               rewind (85)

               write (85) (lgdmic(k),lgdturb(k),vom(k),
     &              depottot(k)/lsun,depottot_surf(k)/lsun,
     &              depottot_core(k)/lsun,
     &              depotwaves(k)/lsun,Dondes(k),Dondeschim(k),
     &              lumwaves(k),lumwave_core(k),
     &              lumondes_surf(k),lumondes_core(k),lumondestot(k),
     &              brunt(k),vr(k)/rsun,vxpsi(k),           
     &              Dmicro(k),vmicro(k),Dthc(k),phiKS(k),
     &              deltaKS(k),tautime(k,ippg),Ereac(k,ippg),tautime(k
     &              ,ipdg), Ereac(k,ipdg),tautime(k,i2he3),Ereac(k
     &              ,i2he3),tautime(k,ihe3ag), Ereac(k,ihe3ag),tautime(k
     &              ,ibe7beta),Ereac(k,ibe7beta), tautime(k,ili7pa)
     &              ,Ereac(k,ili7pa),tautime(k,ibe7pg), Ereac(k,ibe7pg)
     &              ,tautime(k,ib8beta),Ereac(k,ib8beta), tautime(k
     &              ,ic13pg),Ereac(k,ic13pg),tautime(k,in14pg), Ereac(k
     &              ,in14pg),tautime(k,icpg),Ereac(k,icpg), k = 1,nmod)

            else
               write (85) (Dmicro(k),vmicro(k),Dthc(k),phiKS(k),
     $              Dturbul(k),Dbar(k),Dkyle(k),coefDtacho(k),                                ! dturbul 10/2019 + dbar
     $              xvdiffhe(k),xvdiffhenoc(k),xvdiffc12(k)    ! Ajout TD Fev.2018
     $              ,xvdiffc12noc(k),xvdiffc13(k)
     $              ,xvdiffc13noc(k), xvdiffn14(k),xvdiffn14noc(k)
     $              ,xvdiffn15(k),xvdiffn15noc(k),
     $              xvdiffo16(k),xvdiffo16noc(k)
     $              ,xvdiffo17(k),xvdiffo17noc(k),xvdiffo18(k)
     $              ,xvdiffo18noc(k), xvdiffne20(k),xvdiffne20noc(k),
     $              xvdiffna23(k),xvdiffna23noc(k),
     $              deltaKS(k),tautime(k,ippg),Ereac(k        ! Fin ajout TD Fev.2018
     $              ,ippg), tautime(k,ipdg), Ereac(k,ipdg),tautime(k
     $              ,i2he3), Ereac(k,i2he3),tautime(k,ihe3ag),Ereac(k
     $              ,ihe3ag), tautime(k,ibe7beta),Ereac(k,ibe7beta),
     $              tautime(k,ili7pa),Ereac(k,ili7pa),tautime(k,ibe7pg),
     $              Ereac(k,ibe7pg),tautime(k,ib8beta),Ereac(k,ib8beta),
     $              tautime(k,ic13pg),Ereac(k,ic13pg),tautime(k,in14pg),
     $              Ereac(k,in14pg),tautime(k,icpg),Ereac(k,icpg), k = 1
     $              ,nmod)
               print *,'Vmic He resultpr',xvdiffhe(1000)
c               print *, xvdiffc12(k),xvdiffc12noc(k),xvdiffc13(k)
c     $              ,xvdiffc13noc(k)
            endif
         endif
         if (astero) then
            rewind(86)
            write(86) (dkapdt(k),dkapdro(k),denucldt(k),denucldro(k),
     &           k=1,nmod)
         endif
      endif
c$$$      if (time*seci/1e9.ge.4.565.and.time*seci/1e9.le.4.566) then
c$$$         write(444,'()') 
c$$$         write(444,'()')
c$$$         write(444,'(6x,1pe12.6,11x,0pf7.0,9x,1pe16.6)')
c$$$     &        time*seci/1e9, teff, r(nmod)/6.9551e10
c$$$         write(444,'(8x,f8.4)') totm
c$$$         write(444,'()') 
c$$$         write(444,'()')
c$$$         write(444,'()') 
c$$$
c$$$         write(444,44)
c$$$         DO k = nmod,1,-1
c$$$            write(444,'(1x,i4,2(1x,1pe16.6),30(1x,1pe10.4))')
c$$$     &           abs(k-nmod-1),1-m(k)/m(nmod),r(k),xsp(k,ih1),
c$$$     &           xsp(k,ihe3),xsp(k,ihe4),xsp(k,ic12),xsp(k,ic13),
c$$$     &           xsp(k,in14),xsp(k,io16),xsp(k,ine20)+xsp(k,27)+
c$$$     &           xsp(k,28),xsp(k,ina23),xsp(k,img24)+xsp(k,img25)+
c$$$     &           xsp(k,36),xsp(k,ial27),xsp(k,isi28)+xsp(k,42)+
c$$$     &           xsp(k,43),xsp(k,ip31),xsp(k,is32)+xsp(k,46)+xsp(k,47)
c$$$     &           +xsp(k,50),
c$$$     &           xsp(k,icl35)+xsp(k,icl36),0.d0,0.d0,0.d0,0.d0,0.d0,
c$$$     &           0.d0,xsp(k,53),0.d0,0.d0,xsp(k,ih2),xsp(k,ili6),
c$$$     &           xsp(k,ili7),xsp(k,ibe9),xsp(k,11),xsp(k,ib11)
c$$$         ENDDO   
c$$$
c$$$ 44      format (2x,'mesh',3x,'dm/m',6x,'r(cm)',8x,'X(H)',6x,'X(3He)',
c$$$     &        6x,'X(4He)',6x,'X(12C)',6x,
c$$$     &        'X(13C)',5x,'X(14N)',5x,'X(16O)',5x,
c$$$     &        'X(Ne)',5x,'X(Na)',5x,'X(Mg)',5x,'X(Al)',5x,'X(Si)',
c$$$     &        5x,'X(P)',5x,'X(S)',5x,'X(Cl)',5x,'X(Ar)',5x,'X(K)',
c$$$     &        5x,'X(Ca)',5x,'X(Ti)',5x,'X(Cr)',5x,'X(Mn)',5x,
c$$$     &        'X(Fe) = Xheavy',5x,'X(Ni)',5x,'X(Cc)',5x,'X(2D)',5x,
c$$$     &        'X(6Li)',5x,'X(7Li)',5x,'X(9Be)',5x,'X(10B)',5x,'X(11B)')
c$$$      
c$$$         write(555,'()') 
c$$$         write(555,'()')
c$$$         write(555,'(6x,1pe12.6,11x,0pf7.0,9x,1pe16.6)')
c$$$     &        time*seci/1e9, teff, r(nmod)/6.9551e10
c$$$         write(555,'(8x,f8.4)') totm
c$$$         write(555,'()') 
c$$$         write(555,'()')
c$$$         write(555,'()') 
c$$$
c$$$         write(555,55)
c$$$         DO k = nmod,1,-1
c$$$            write(555,'(1x,i4,7(1x,1pe16.6),1x,1pe16.6,1x,1pe16.6,
c$$$     &      1x,1pe16.6,1x,1pe16.6,1x,1pe16.6,1x,1pe16.6,
c$$$     &      3(1x,1pe16.6),3(1x,1pe16.6),1x,1pe16.6,1x,i4)')
c$$$     &      abs(k-nmod-1),1-m(k)/m(nmod),r(k),lum(k)/lsun,t(k),
c$$$     &      P(k),lnf(k),ro(k),0.d0,abrad(k)-abad(k),abad(k),
c$$$     &      kap(k),bruntV2(k),abla(k),Egrav(k),Enucl(k),mu(k),0.d0
c$$$     &      ,0.d0,0.d0,gamma1(k),crz(k)
c$$$         ENDDO
c$$$          
c$$$         print *,'k',k
c$$$          
c$$$ 55      format (2x,'mesh',6x,'dm/m',8x,'r(cm)',8x,'L',8x,'Teff',8x,
c$$$     &        'P',8x,'lnf',8x,'rho',8x,'Ne',8x,'gr-ga',8x,
c$$$     &        'dT_ad',8x,'kappa',8x,'N**2',8x,'dlnT/dlnP',8x,'e_grav',
c$$$     &        8x,'e_nuc',8x,'mu',8x,'fHI',8x,'fHeI',8x,'fHeII',
c$$$     &        8x,'Gamma1',8x,'crz')
c$$$  endif

c      print *, 'Vdiff He', xvdiffhe(0:100)
c      stop
     
c..   Storage of values for atomic diffusion (modif AP TD Jan.2019)
c      if(model.eq.(modeli+1).and.microdiffus) then
c         do k=1,nmod
c            write(654,'(1x,i4,9(1x,1pe11.4))'),k,T(k),zmeanO16(k),xsp(k
c     $           ,io16),dens_elec(k),coulomb(k),Zpaq(k),Zpaq1(k)
c     $           ,Zpaq2(k),Kpaq(k)
c         enddo
c      endif
      
      
 1000 format (' Final mass = ',f11.8,', dm = ',f11.9,', age = ',1pe12.6,
     &     ' yr, dt = ',1pe13.7,' yr, Teff = ',0pf7.0,', L = ',1pe10.4)
 1100 format (5x,'convective limits:')
 1200 format (2x,'Schwar.: ',i4,' -->',i4,' (',1pe9.3,',',1pe9.3,
     &     ')   turnover: ',1pe9.3,'yr   Overshoot (dm): ',i4,' -->',i4,
     &     ' (',1pe9.3,',',1pe9.3,')')
 1210 format (2x,'Schwar.: ',i4,' -->',i4,' (',1pe9.3,',',1pe9.3,
     &     ')   turnover: ',1pe9.3,'yr')
 1230 format (2x,'Ledoux : ',i4,' -->',i4,' (',1pe9.3,',',1pe9.3,
     &     ')   turnover: ',1pe9.3,'yr')
 1240 format (2x,'Ledoux : ',i4,' -->',i4,' (',1pe9.3,',',1pe9.3,
     &     ')   turnover: ',1pe9.3,'yr   Overshoot (dm): ',i4,' -->',i4,
     &     ' (',1pe9.3,',',1pe9.3,')')
 1250 format (2x,'semiconvection : ',i4,' -->',i4,' (',1pe9.3,',',
     &     1pe9.3,')')
 1260 format (2x,'Thermohaline : ',i4,' -->',i4,' (',1pe9.3,',',
     &     1pe9.3,')')
 1270 format (//,'o PRESCRIPTIONS FOR ROTATIONNAL MIXING',/,3x,
     &     'Critical Richardson number Ric = 0.25',/,3x,
     &     'Prescription for horizontal turbulence diffusion',
     &     ' coefficient Dh : ',1x,a8,/)
C Modifs CC ondes (2/04/07) -->
c 8101 format (2x,i4,1x,3(1x,1pe17.10),2(1x,1pe11.4),(1x,1pe17.10),
c     &      7(1x,1pe11.4),(1x,1pe17.10),3(1x,1pe11.4),(1x,1pe17.10),
c     & 2(1x,1pe11.4),8(1x,1pe11.4),(1x,1pe17.10),(1x,1pe11.4))
C <--
 8100 format('# nsh lgdmic lgdturb vomega depottot depottot_surf',
     &     'depottot_core depotwaves Dondes Dondeschim lumwave',
     &     'lumwave_core brunt_o lumondes vr vpsi',
     &     'Dmicro',6x,'Vmicro',5x,'Dthc',6x,'phiKS',6x
     $     ,'deltaKS',6x,'Tpp',6x,'Epp',6x,'Tpd',6x,'Epd',6x,'T2he3',6x
     $     ,'E2he3',6x,'The3he4',6x,'Ehe3he4',6x,'Tbe7b',6x,'Ebe7b',6x
     $     ,'Tli7p' ,6x,'Eli7p' ,6x,'Tbe7p',6x,'Ebe7p',6x,'Tb8b',6x
     $     ,'Eb8b',6x ,'Tc13p',6x ,'Ec13p',6x ,'Tn14p',6x,'En14p',6x
     $     ,'Tc12p',6x ,'Ec12p')
 8101 format(4x,'i',4x,'Dmicro',6x,'Vmicro',5x,'Dthc',6x,'phiKS',6x
     $     ,'deltaKS',6x,'Tpp',6x,'Epp',6x,'Tpd',6x,'Epd',6x,'T2he3',6x
     $     ,'E2he3',6x,'The3he4',6x,'Ehe3he4',6x,'Tbe7b',6x,'Ebe7b',6x
     $     ,'Tli7p' ,6x,'Eli7p' ,6x,'Tbe7p',6x,'Ebe7p',6x,'Tb8b',6x
     $     ,'Eb8b',6x ,'Tc13p',6x ,'Ec13p',6x ,'Tn14p',6x,'En14p',6x
     $     ,'Tc12p',6x ,'Ec12p')


 8102 format (2x,i4,5(1x,1pe11.4),26(1x,1pe11.4),5(1x,1pe11.4),
     &     22(1x,e15.4E2))
 8103 format (2x,i4,5(1x,1pe11.4),22(1x,1pe14.6))

      return
      end
