
      PROGRAM BINLIST_EVOL

************************************************************************
* Create the ascii output (the listing) from binary output evol file   *
* (structure and chemical profiles)                                    *
************************************************************************

      implicit none

      include 'evolpar.chem'
      include 'evolpar.conv'
      include 'evolpar.str'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.ent'
      include 'evolcom.eos'
      include 'evolcom.grad'
      include 'evolcom.ion'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.opa'
      include 'evolcom.rot'
      include 'evolcom.scov'
      include 'evolcom.shk'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.transp'
      include 'evolcom.turb'
      include 'evolcom.var'

      character*50 name
      character varstr*5,cmodel*6
      character*37 react2
      character*1 rcrzc

      integer imodtot,imod,err,error
      integer lenc,lenci,lenca,sm
      integer idiffnuc,ispec
      integer nratsh,nsequence
      integer ienvb,ienvt,nrat
      integer nvarstrm,nvarstr
      integer i,ii,ic,j,k,kl,l,mm,mn
      integer crzi,nspr
      integer klenv,klpulse,klcore,itop,ibase
      integer ntilde4

      logical idiffcc,idiffcc0

      double precision dturb,dmic,vdmic,dtinv,vvro,cd,cd0,Dherw
      double precision shlimplt 
      double precision v,enc,dumarray
      double precision totmcor,totmenv,ebind
      double precision t9,eng,eerho,eerhocap,t9k
      double precision version
      double precision Dhold
      double precision xNt,xNu,xNr,geffom
      double precision abmurj,abmuj,mul
      double precision xspr
      double precision epot,eint,ebindenv,epotenv,eintenv
      double precision comp,heat,work,etherm
      double precision gg,tt,renucl
      double precision Dtot,Dsc,fover
      double precision abadm,phiKSm,deltaKSm,khim,KKtt,alphasc
      double precision geffs,tnucg,tff

      parameter (nvarstrm = 44)
      double precision dift,dife,summe,rpvisc,rsconv,dlumdt,denucdt
c      double precision summe,rpvisc,rsconv

      common /difcirc/ dift(nsh),dife(nsh)
      common /coefdiff/ dturb(nsh),dmic(nsh),vdmic(nsh)
      common /diffcode/ idiffnuc,idiffcc,idiffcc0
      common /nnnk/ xNt(nsh),xNu(nsh),xNr(nsh),geffom(nsh)
      common /calcDh/ Dhold(nsh)
      common /gradmuj/ abmurj(nsh),abmuj(nsh),mul(nsh)
      common /overshoot/ klcore,klenv,klpulse

      dimension v(nreac),enc(nreac),nratsh(nsh),react2(nreac),
     &     varstr(nvarstrm),summe(nsh),rpvisc(nsh),rsconv(nsh),vvro(nsh)
c      dimension xKt(nsh),dift(nsh),dife(nsh),Dhold(nsh),vtheta(nsh)
      dimension tt(nsh),gg(nsh),Dtot(nsh),Dsc(nsh),Dherw(nsh)

      dimension crzi(nsh),rcrzc(nsh)
      dimension ispec(0:nbz,2*nbz)
      dimension xspr(nsh,nsp)
      dimension comp(nsh),heat(nsh),work(nsh),etherm(nsh),dlumdt(nsh),
     &     denucdt(nsh),cd(nsh),renucl(nsh)

      pi = 3.1415926535d0
      sec = 3.1557807d7
      seci = 1.d0/sec
      c = 2.99792458d10
      g = 6.67259d-8
      boltz = 1.380658d-16
      h = 6.6260755d-27
      sig = 2.d0*pi**5*boltz**4/(15.d0*h**3*c*c)
      sigc = 4.d0*sig/c
      avn = 6.0221367d23
      rk = avn*boltz
      ech = 1.60217733d-19
      econv = ech*avn*1.d13
      msun = 1.9891d33
      rsun = 6.9599d10
      lsun = 3.846d33
      cmodel = ' '
c      code_version = 2.1d0
      code_version = 9.15d0

      pim2 = 2.d0*pi
      pim4 = 4.d0*pi
      pim8 = 8.d0*pi
      pw14 = 1.d0/4.d0
      pw13 = 1.d0/3.d0
      pw23 = 2.d0/3.d0
      pw34 = 3.d0/4.d0
      pw43 = 4.d0/3.d0
      pw53 = 5.d0/3.d0
      vlog10 = log(10.d0)
      error = 0

c..   imod : model number in the binary file (usually = 1)
c..   sm = 0 : generate ascii file
c..   sm = 1 : generate smfile : filename.p[i]

      read (5,999) imod,sm,name
c..   ex : 1 1 shock2007_0200b

      write (*,*) 'Welcome to the binlist qui marche pas...'

      write (*,*)'filename : ',trim(name),', imod = ',imod
      i=lenci(name)
      cmodel = name(i+1:i+5)
      call make_number (cmodel,nsequence)

      write (*,*) 'toto'

      if (imod.eq.9) then
         write (*,*) 'coucou', imod 
         open (unit = 90,file = '/home/decressin/STAREVOL/' //
     &        'VISUSM/' // 
     &        trim(name),form = 'unformatted')
         open (unit = 99,
     &        file ='/home/decressin/STAREVOL/VISUSM/starevol.par',
     &        status = 'unknown')
         read (99,2) version
 2       format (/,43x,f4.2) 
         rewind(99)
      else
         open (unit = 90,file = '/home/decressin/STAREVOL/VISUSM/' 
     &        // trim(name) // 'b',form = 'unformatted')
         open (unit = 99,file ='/home/decressin/STAREVOL/VISUSM/' 
     &        // trim(name) // 'd',status = 'unknown')

         read (99,3) version
 3       format (/////,43x,0pf4.2)
         backspace(99)
         backspace(99)
      endif
      write (*,22) version
 22   format(' STAREVOL version : ',f4.2)
c      open (unit = 92,file = '/tmp/nextini.bin',form = 'unformatted',
c     &     status = 'unknown')

*_______________________________
***   reading of the binary file
*-------------------------------
      
      if (version.lt.2.05d0) then
         call rmodpar_2_03
      else
         call rmodpar
      endif

      close(99)
      mixopt = .false.
      lnucl = .true.
      diffusconv = .true.
      if (idiffty.ge.8.and.idiffty.le.15) then
         if (imod.eq.9) then
            close (93)
            close (95)
            open (unit = 93,file = '/home/decressin/STAREVOL/' //
     &           'VISUSM/' // 
     &           name(1:lenca(name)) // 'ang.bin',form = 'unformatted',
     &           status = 'unknown')
         else
            open (unit = 93,file = '/home/decressin/STAREVOL/' //
     &           'VISUSM/' // 
     &           trim(name) // 'a',form = 'unformatted',
     &           status = 'unknown')
         endif
      endif

      if (version.ge.2.1d0) then
         open (unit = 99,
     &        file ='/home/decressin/STAREVOL/BATCH/' //
     &        'starevolnuc_evol2.10.par',
     &        status = 'unknown')
      else
         open (unit = 99,
     &        file ='/home/decressin/STAREVOL/BATCH/' //
     &        'starevolnuc_evolsn.par',
     &        status = 'unknown')
      endif
      read (99,5)
 5    format(1x,//)

************************************************************************
*  setting the directory where the opacity tables are to be read       *
************************************************************************

      call set_opal_dir( '/home/decressin/STAREVOL_2.30/' //
     &     'EVOL/LIB/' )

      call rinimod (90,93,92,95,0)
      print *,'Total number of shells : ',nmod

      call eos (0,error)

      no = 0

      dtinv = 1.d0/dtn
      do i = 1,nmod
         do l = 1,nsp
            xspr(i,l) = xsp(i,l)
         enddo
         rcrzc(i) = crzc(i)
         renucl(i) = venucl(i)
         rpvisc(i) = vpvisc(i)
         rsconv(i) = vsconv(i)
         t(i) = exp(lnt(i))
         vt(i) = exp(vlnt(i))
         vvro(i) = vro(i)
         if (vfacc(i).gt.0) iaccr = 2
         if (idiffty.ge.8.and.idiffty.le.15) then
            Dhd(i) = dift(i)+dife(i)
            Dtot(i) = Dhd(i)+convdiff(i)
         endif
      enddo
*______________________________________________________________________
***   calculation of the factor ft and fp for the treatment of rotation
*----------------------------------------------------------------------

      if (hydrorot) then
         call potential
      else
         do k = 1,nmod
            ft(k) = 1.d0
            fp(k) = 1.d0
         enddo
      endif

      if (crzc(1).ne.'c') Dtot(1) = Dtot(2)
      mprot = 1.6726231d-24
      saha = (sqrt(pim2*boltz*9.1093897d-28)/h)**3 !2.4147031864d15
      call thermo (error)
      call diffinit

      alphasc = 4.d-2
      do j = 1,nmod1
         if (rcrzc(j).eq.'s') then
            khim = 0.5d0*(khi(j)+khi(j+1))
            Dsc(j) = alphasc*khim/(6.d0*rom(j)*cpm(j))*
     &           (abrad(j)-abm(j))/(abled(j)-abrad(j))
         endif
      enddo

c..   diffusive overshoot : Herwig
      klenv = 0
      klcore = 0
      klpulse = 0
      if (diffover.and.nsconv.gt.0) then
         if (novlim(1,3).eq.1) klcore = 1
         do kl = 1,nsconv
            if ((tau(novlim(kl,4)).lt.1.d2.or.t(novlim(kl,4)).lt.1.d6)
     &           .and.klenv.eq.0) klenv = kl
            if (t(novlim(kl,3)).gt.2.d8.and.novlim(kl,3).gt.1.and.
     &           nphase.eq.5.and.klpulse.eq.0) klpulse = kl
         enddo
         if (idiffty.eq.24.or.idiffty.eq.25) klenv = 0
         if (idiffty.ne.26.and.idiffty.ne.27) klcore = 0
         if (idiffty.eq.23.or.idiffty.eq.27) klpulse = 0

c..   diffusion below convective envelope
         if (klenv.ne.0) then
            itop = novlim(klenv,3)
            fover = etaturb
            cd0 = convdiff(itop)
            call diffherwig (itop,ibase,cd0,Dherw,fover,1.d0,1)
            novlim(klenv,7) = ibase
         endif
c..   diffusion above pulse
         if ((idiffty.eq.24.or.idiffty.eq.26).and.klpulse.ne.0) then
            ibase = novlim(klpulse,4)
            fover = etaturb
            cd0 = convdiff(ibase)
            call diffherwig (ibase,itop,cd0,Dherw,fover,1.d0,-1)
            novlim(klpulse,8) = itop
         endif
c..   diffusion below pulse
         if ((idiffty.eq.25.or.idiffty.eq.26).and.klpulse.ne.0) then
            itop = novlim(klpulse,3)
            fover = 1.d-4
            cd0 = convdiff(itop)
            call diffherwig (itop,ibase,cd0,Dherw,fover,1.d0,1)
            novlim(klpulse,7) = ibase
         endif
c..   diffusion above core
         if (klcore.ne.0) then
            ibase = novlim(klcore,4)
            fover = etaturb
            cd0 = convdiff(itop)
            call diffherwig (ibase,itop,cd0,Dherw,fover,1.d0,-1)
            novlim(klcore,8) = itop
         endif
      endif
      do j = 1,nmod
         dturb(j) = convdiff(j)+Dsc(j)+Dherw(j)
      enddo

*____________________________________________________
***   calculation of the gravitational binding energy
*----------------------------------------------------

      do i = 1,nmod
         comp(i) = (ro(i)-vro(i))*dtinv/ro(i) !d ln rho/dt
         heat(i) = (lnt(i)-vlnt(i))*dtinv     !d ln T/dt 
         work(i) = -p(i)*(1.d0/ro(i)-1.d0/vro(i))*dtinv         !-Pd(1/rho)/dt
         etherm(i) = (ve(i)-e(i))*dtinv       !-de/dt
         dlumdt(i) = (lum(i)-vlum(i))*dtinv/(lum(i)+1.d-30)
         denucdt(i) = (s(i)-vs(i))*dtinv/(s(i))
c..   egrav = -Pd(1/rho)/dt-de/dt)
      enddo
      totmenv = 0.d0
      totmcor = 0.d0
      do i = 1,nmod
         if (rcrzc(i).ne.'r'.and.tot(i).lt.1.d50) then
            totmcor = totmcor+tot(i)*dm(i)
         else
            goto 10
         endif
      enddo
 10   if (i.gt.1) then
         totmcor = totmcor/(sec*m(i+1))
      endif
      ienvt = 0
      do i = nmod1,1,-1
         if (rcrzc(i).ne.'r'.and.tot(i).lt.1.d50) then
            totmenv = totmenv+tot(i)*dm(i)
         else
            if (tau(i).lt.1.d2.and.totmenv.lt.1.d-10) then
               ienvt = i-1
            else
               goto 20
            endif
         endif
      enddo
 20   ienvb = i+1
      if (ienvb.lt.ienvt) then
         totmenv = totmenv/(sec*(m(ienvt+1)-m(ienvb)))
      endif

      ebind = 0.d0
      epot = 0.d0
      eint = 0.d0
      ebindenv = 0.d0
      epotenv = 0.d0
      eintenv = 0.d0
      do k = 1,nmod-1
         dm(k) = mr(k+1)-mr(k)
         eint = eint+e(k)*dm(k)
         epot = epot+g*m(k)/r(k)*dm(k)
         if (nsconv.gt.0.and.k.ge.novlim(nsconv,3)) then
            eintenv = eintenv+e(k)*dm(k)
            epotenv = epotenv+g*m(k)/r(k)*dm(k)
         endif
      enddo
      ebind = eint-epot
      ebindenv = eintenv-epotenv
      write (*,110) ebind,eint,epot
      write (*,111) ebindenv,eintenv,epotenv
 110  format(/,'STAR     : Ebind =',1pe11.4,' erg, Eint =',1pe11.4
     &     ,' erg, Epot =',1pe11.4,' erg')
 111  format('envelope : Ebind =',1pe11.4,' erg, Eint =',1pe11.4
     &     ,' erg, Epot =',1pe11.4,' erg')

*____________________________________________________
***   storage of the corresponding ascii listing file
*----------------------------------------------------
      
      if (sm.eq.0) then
         geffs = log10(g*m(nmod)/(r(nmod)*r(nmod)))
         tnucg = 1.1d10*m(nmod)*lsun/abs(lum(nmod)*msun)
         tff = sqrt(r(nmod)**3/(g*m(nmod)))
         call make_char (model,cmodel)
c         write(cmodel,'(i4)') model
         open (unit = 10,file = name(1:lenci(name)) // cmodel //
     &        '.lst', status = 'unknown')
         write (10,*)
         if (nphase.le.1) write (10,*) " Phase 1: Pre-Central H Burning"
         if (nphase.eq.2) write (10,*) " Phase 2: Central H Burning"
         if (nphase.eq.3) write (10,*)
     &        " Phase 3: Post-Central H Burning"
         if (nphase.eq.4) write (10,*) " Phase 4: Central He Burning"
         if (nphase.eq.5) write (10,*)
     &        " Phase 5: Post-Central He Burning"
         if (nphase.eq.6) write (10,*) " Phase 6: Central C Burning"
         if (nphase.ge.7) write (10,*)
     &        " Phase 7: Post-Central C Burning"
         write (10,1000)
         write (10,1000) 1
         write (10,1100) model
         write (10,1200) totm,dma,time,dtn
         if (iaccr.gt.0) then
            write (10,1250) dmacc1,ddmacc,dtnacc
            write (10,1260) taudnuc0,taudnuc,tauaccr,itauac
            write (10,1270) ldacc,lacc
         endif
         write (10,1300) tau0
         write (10,1400) lum(nmod)/lsun,r(nmod)/rsun,t(nmod),ro(nmod),
     &        geffs,totm
         write (10,1500) taueff
         write (10,1600) soll,solr,tsurf,roeff,geff,m(neff)/msun
         write (10,1800) neff,nmod,mlsh
         write (10,1900) nsconv,alphac
         do kl = 1,nsconv
            write (10,2000) novlim(kl,1),novlim(kl,2),novlim(kl,3),
     &           novlim(kl,4)
         enddo
         write (10,2100) tff,tkh,trdiff,tnucg
         write (10,2200) totmcor
         write (10,2300) totmenv
         write (10,2400) lh,lhe
         write (10,2500) tirrm1,ddm1,tirrm2
c         write (10,2600) depot,dekin,deint,edyn,eephot,eeneu,eenuc,etot

         write (10,*)
         write (10,2700)
         do i = 1,nmod
            write (10,2800) i,rcrzc(i),r(i),t(i),ro(i),p(i),beta(i),
     &           eta(i),lnf(i),sr(i),lum(i),u(i),mr(i)
         enddo
         write (10,2900)
         do i = 1,nmod
            write (10,3000) i,rcrzc(i),tau(i),kap(i),mu(i),mue(i),
     &           abad(i),abmu(i),abrad(i),abla(i),cs(i),cp(i),
     &           gamma(i),e(i),rpvisc(i)+1.d-20
         enddo
         write (10,3100)
         do i = 1,nmod
            write (10,3200) i,rcrzc(i),tot(i),rsconv(i),fconv(i),
     &           abs(frad(i)),renucl(i),renucl(i)+enupla(i),egrav(i)
     &           ,hydrat(i),dturb(i),dmic(i),vdmic(i),enupla(i)
         enddo
         if (irotbin.gt.0) then
            write (10,3330)
            do i = 1,nmod
               write (10,3360) i,rcrzc(i),omega(i),urs(i),dift(i),
     &              dife(i)
            enddo
         endif
         if (iaccr.gt.0) then
            if (irotbin.eq.0) then
               write (10,3231)
               do i = 1,nmod
                  write (10,3261) i,rcrzc(i),omega(i),vfacc(i),macc(i),
     &                 eacc(i),eshr(i),tacc(i)
               enddo
            else
               write (10,3232)
               do i = 1,nmod
                  write (10,3262) i,rcrzc(i),vfacc(i),macc(i),eacc(i),
     &                 eshr(i),tacc(i)
               enddo
            endif
         endif

         write (10,3500)
         write (10,3600) (elem(l),l = 1,10)
         do i = 1,nmod
            write (10,3700) i,rcrzc(i),(xspr(i,l),l = 1,10)
         enddo
         write (10,3800)
         write (10,3900) (mtotlos(l),l = 1,10)
         write (10,3600) (elem(l),l = 11,20)
         do i = 1,nmod
            write (10,3700) i,rcrzc(i),(xspr(i,l),l = 11,20)
         enddo
         write (10,3800)
         write (10,3900) (mtotlos(l),l = 11,20)
         write (10,3600) (elem(l),l = 21,30)
         do i = 1,nmod
            write (10,3700) i,rcrzc(i),(xspr(i,l),l = 21,30)
         enddo
         write (10,3800)
         write (10,3900) (mtotlos(l),l = 21,30)
         write (10,3600) (elem(l),l = 31,40)
         do i = 1,nmod
            write (10,3700) i,rcrzc(i),(xspr(i,l),l = 31,40)
         enddo
         write (10,3800)
         write (10,3900) (mtotlos(l),l = 31,40)
         write (10,3600) (elem(l),l = 41,50)
         do i = 1,nmod
            write (10,3700) i,rcrzc(i),(xspr(i,l),l = 41,50)
         enddo
         write (10,3800)
         write (10,3900) (mtotlos(l),l = 41,50)
         write (10,4000) (elem(l),l = 51,nsp-1)
         do i = 1,nmod
            write (10,4100) i,rcrzc(i),(xspr(i,l),l = 51,nsp-1)
         enddo
         write (10,3800)
         write (10,4200) (mtotlos(l),l = 51,nsp-1)

      else
*______________________________________________
***   generation of ascii files for Super Mongo
*----------------------------------------------

         call make_char (model,cmodel)
         open(unit=10,file='/home/decressin/STAREVOL/SMONGO/DATA/'// 
     &        name(1:lenci(name)) //cmodel// '.p0',status = 'unknown')
         open(unit=11,file='/home/decressin/STAREVOL/SMONGO/DATA/'// 
     &        name(1:lenci(name)) //cmodel// '.p1',status = 'unknown')
         open(unit=12,file='/home/decressin/STAREVOL/SMONGO/DATA/'// 
     &        name(1:lenci(name)) //cmodel// '.p2',status = 'unknown')
         open(unit=13,file='/home/decressin/STAREVOL/SMONGO/DATA/'// 
     &        name(1:lenci(name)) //cmodel// '.p3',status = 'unknown')
         open(unit=14,file='/home/decressin/STAREVOL/SMONGO/DATA/'// 
     &        name(1:lenci(name)) //cmodel// '.p4',status = 'unknown')
         open(unit=15,file='/home/decressin/STAREVOL/SMONGO/DATA/'// 
     &        name(1:lenci(name)) //cmodel// '.p5',status = 'unknown')
         open(unit=16,file='/home/decressin/STAREVOL/SMONGO/DATA/'// 
     &        name(1:lenci(name)) //cmodel// '.p6',status = 'unknown')
         open(unit=17,file='/home/decressin/STAREVOL/SMONGO/DATA/'// 
     &        name(1:lenci(name)) //cmodel// '.p7',status = 'unknown')
         open(unit=18,file='/home/decressin/STAREVOL/SMONGO/DATA/'// 
     &        name(1:lenci(name)) //cmodel// '.p8',status = 'unknown')
         open(unit=19,file='/home/decressin/STAREVOL/SMONGO/DATA/'// 
     &        name(1:lenci(name)) //cmodel// '.p9',status = 'unknown')

         write (10,1205)
         write (10,1210) nmod,totm,lum(nmod)/lsun,r(nmod)/rsun,time
     &        ,iaccr,irotbin,nphase,nmod,nsp-1,model,dtn,nsequence
     &        ,int(version*1.d2)
         write (11,1211)
         write (12,1212)
         write (13,1213)
c..   chemical profiles
         if (version.lt.2.1d0) then
            write (*,*) ' Old network (version < 2.10) detected'
            write (14,2214)
            write (15,2215)
            write (16,2216)
            write (17,2217)
            write (18,2218)
         else
            write (*,*) ' URCA network (version >= 2.10) detected'
            write (14,1214)
            write (15,1215)
            write (16,1216)
            write (17,1217)
            write (18,1218)
         endif
         write (19,1219)
         do i = 1,nmod
            if (rcrzc(i).eq.'c') then
               crzi(i) = 1
            else
               crzi(i) = 0
            endif
            write (11,2810) i,crzi(i),r(i)/rsun,t(i),ro(i),p(i),beta(i),
     &           eta(i),lnf(i),sr(i),lum(i)/lsun,u(i),mr(i),accel(i)
         enddo
         do i = 1,nmod
            if (abrad(i).gt.1.d35) abrad(i) = 1.d35
            write (12,3010) i,tau(i),kap(i),mu(i),mue(i),
     &           abad(i),abmu(i),abrad(i),abla(i),cs(i),cp(i),
     &           gamma(i),e(i),rpvisc(i)+1.d-20
         enddo
         do i = 1,nmod
            if (tot(i).gt.1.d37) tot(i) = 1.d37
            write (13,3210) i,tot(i),abs(rsconv(i)),fconv(i),
     &           frad(i),renucl(i),enupla(i),egrav(i),hydrat(i),
     &           dturb(i),dmic(i),vdmic(i),enunucl(i)
         enddo

         do i = 1,nmod
            do l = 1,nsp
               xspr(i,l) = max(xspr(i,l),1.d-35)
            enddo
         enddo
         do i = 1,nmod
            write (14,3710) i,(xspr(i,l),l = 1,10)
         enddo
         do i = 1,nmod
            write (15,3710) i,(xspr(i,l),l = 11,20)
         enddo
         do i = 1,nmod
            write (16,3710) i,(xspr(i,l),l = 21,30)
         enddo
         do i = 1,nmod
            write (17,3710) i,(xspr(i,l),l = 31,40)
         enddo
         do i = 1,nmod
            write (18,3710) i,(xspr(i,l),l = 41,50)
         enddo
         do i = 1,nmod
            write (19,4110) i,(xspr(i,l),l = 51,nsp-1)
         enddo
         
         write(*,*)irotbin,'irotbin'
         if (iaccr.gt.0) then
            open(unit=20,file='/home/decressin/STAREVOL/' //
     &           'SMONGO/DATA/' // 
     &         name(1:lenci(name)) //cmodel//'.p10',status = 'unknown')
            if (irotbin.eq.0) then
               write (20,1221)
               do i = 1,nmod
                  write (20,3461) i,vomega(i),vfacc(i),macc(i),
     &                 eacc(i),eshr(i),tacc(i)
               enddo
            else
               open(unit=21,file='/home/decressin/STAREVOL/' //
     &              'SMONGO/DATA/' //
     &         name(1:lenci(name)) //cmodel//'.p11',status = 'unknown')
               write (20,1222)
               write (21,1223)
               do i = 1,nmod
                  write (20,3462) i,vfacc(i),macc(i),eacc(i),
     &                 eshr(i),tacc(i)
                  write (21,3460) i,vomega(i),urs(i)*ur_Ss,V_circ(i)
     &                 *ur_Ss,dift(i),dife(i),Dtot(i),Dhold(i),thetas(i)
     &                 ,xlambdas(i),abmuj(i),xnuvv(i),xnum(i),xKt(i)
               enddo
            endif
         elseif (irotbin.gt.0) then
            open(unit=20,file='/home/decressin/STAREVOL/' //
     &           'SMONGO/DATA/' // 
     &         name(1:lenci(name)) //cmodel//'.p10',status = 'unknown')
            write (20,1223)
            do i = 1,nmod
               write (20,3460)i,vomega(i),urs(i)*ur_Ss,dift(i),dife(i)
     &              ,Dtot(i),Dhold(i),thetas(i),xlambdas(i),abmuj(i)
     &              ,xnuvv(i),xnum(i),xKt(i)
            enddo
         endif
         open(unit=22,file='/home/decressin/STAREVOL/SMONGO/DATA/'// 
     &         name(1:lenci(name)) //cmodel//'.p12',status = 'unknown')
         write (22,1224)
         do i = 1,nmod
            write (22,3465) i,comp(i),heat(i),work(i),etherm(i),
     &           dlumdt(i),denucdt(i)
         enddo
         
         close(10)
         close(11)
         close(12)
         close(13)
         close(14)
         close(15)
         close(16)
         close(17)
         close(18)
         close(19)
         close(20)

         write (*,*)
         write (*,*) 'Files saved as : /home/decressin/STAREVOL/' //
     &        'SMONGO/DATA/' // 
     &        name(1:lenci(name)) //cmodel
         write (*,*) 'readsi ',name(1:lenci(name)) //cmodel
         stop 'sm files generated !'
      endif

*_________________________________________________________
***   determination of the shells where the nuclear energy
***   production rate is maximum and significant
*---------------------------------------------------------

      do i = 1,nmod
         mueinv(i) = 1.d0/mue(i)
      enddo

      tnucmin = 5.50d5
      do i = 1,nmod
         if (t(i).gt.tnucmin) then
            ksh = i
            rok = ro(i)
            tk = t(i)
            t9 = t(i)/1.d9
            mueinvk = mueinv(i)
            rhpsik = 1.d0
            do ic = 1,nsp
               x(ic) = max(1.d-40,xsp(i,ic))
            enddo
            call vit (t9,v)
            call nuceng (0,v,eng,eerho,eerhocap)
         endif
      enddo
      shlimplt = 1.0d1
      nrat = 1
      nratsh(1) = 1
      do i = 2,nmod-1
         if (renucl(i-1).lt.renucl(i).and.renucl(i+1).lt.
     &        renucl(i).and.renucl(i).gt.shlimplt) then
            nrat = nrat+1
            nratsh(nrat) = i
         endif
      enddo
      if (nphase.eq.5.and.nsconv.ge.2) then
         do j = 1,nsconv-1
            k = novlim(j,3)+1
            if (xsp(k,2).lt.1.d-5) then
               nrat = nrat+1
               nratsh(nrat) = k
               nrat = nrat+1
               nratsh(nrat) = novlim(j,4)-1
            endif
         enddo
      endif

*___________________________________________________________
***   storage of the nuclear fluxes at the shells just found
*-----------------------------------------------------------

      write (10,4300) nrat,shlimplt
      do i = 1,nrat
         k = nratsh(i)
         write (10,4400) k,mr(k),renucl(k)
      enddo
      do i = 1,nrat
         k = nratsh(i)
         write (10,4500) k
         write (10,4600) (j,flu(k,j),j = 1,7)
         write (10,4600) (j,flu(k,j),j = 8,14)
         write (10,4600) (j,flu(k,j),j = 15,21)
         write (10,4600) (j,flu(k,j),j = 22,28)
         write (10,4600) (j,flu(k,j),j = 29,35)
         write (10,4600) (j,flu(k,j),j = 36,42)
         write (10,4600) (j,flu(k,j),j = 43,49)
         write (10,4600) (j,flu(k,j),j = 50,56)
         write (10,4600) (j,flu(k,j),j = 57,63)
         write (10,4600) (j,flu(k,j),j = 64,70)
         write (10,4600) (j,flu(k,j),j = 71,77)
         write (10,4600) (j,flu(k,j),j = 78,84)
         write (10,4600) (j,flu(k,j),j = 85,91)
         write (10,4600) (j,flu(k,j),j = 92,98)
         write (10,4600) (j,flu(k,j),j = 99,105)
         write (10,4600) (j,flu(k,j),j = 106,112)
         write (10,4600) (j,flu(k,j),j = 113,119)
         write (10,4600) (j,flu(k,j),j = 120,126)
         write (10,4600) (j,flu(k,j),j = 127,133)
         write (10,4600) (j,flu(k,j),j = 134,140)
         write (10,4600) (j,flu(k,j),j = 141,147)
         write (10,4600) (j,flu(k,j),j = 148,154)
         write (10,4600) (j,flu(k,j),j = 155,161)
         write (10,4600) (j,flu(k,j),j = 162,168)
         write (10,4600) (j,flu(k,j),j = 169,175)
         write (10,4700) (j,flu(k,j),j = 176,nreac)
      enddo
      write (10,4800) totgrav,totnucl

*____________________________________________________________________
***   storage of the nuclear reactions that are the most important
***   contributors to the nuclear energy prodution at the same shells
*--------------------------------------------------------------------

      write (10,4900)
c      totnuclnu = 0.d0
c      do i = 1,nmod-1
c         dm(i) = (mr(i+1)-mr(i))*totm*msun
c         totnuclnu = totnuclnu+renucl(i)*dm(i)
c      enddo
      dm(nmod) = 0.d0
      do j = 1,nreac
         enc(j) = 0.d0
         do i = 1,nmod-1
            enc(j) = enc(j)+flu(i,j)*qi(j)*econv*dm(i)
         enddo
c         enc(j) = enc(j)/totnuclnu*1.d2
      enddo
      call sort (enc,react,react2)
      do j = nreac,1,-1
         if (enc(j).gt.1.d-2) write (10,5000) react2(j),enc(j)
      enddo
      do i = 1,nrat
         k = nratsh(i)
         write (10,5100) k,renucl(k)
         do j = 1,nreac
            enc(j) = flu(k,j)*qi(j)*econv/renucl(k)*1.d2
         enddo
         call sort (enc,react,react2)
         do j = nreac,1,-1
            if (enc(j).gt.1.d-2) write (10,5000) react2(j),enc(j)
         enddo
      enddo
      write (10,5200)

      close (10)

      stop 'binlist'

 877  format(/,'Arg : #1 model number (in file)',/,6x,
     &     '#2 output : 0 ascii (.lst)',/,19x,'1  smongo (.p0-9)',/,19x,
     &     '2  IDL (.l)',/,6x,'#3 filename',/,3x,/,
     &     '  --> ex: 5 1 m1.0z0_0005b')
 997  format (' Not so many models saved in that file!',/)
 998  format (' Not so many listing models!',/,
     &     ' The last model has been selected.')
 999  format (i1,1x,i1,1x,a50)

c..   ascii formats
 1000 format ('@',i2)
 1100 format (5x,'model # ',i6,/)
 1200 format (5x,'mass = ',f14.10,', dM = ',1pe12.6,'; age = ',
     &     1pe21.15,' yr, dt = ',1pe12.6,' yr',/)
 1250 format (5x,'accreted matter: dM(true) = ',1pe12.6,', ddM = ',
     &     1pe13.6,', dtnacc = ',1pe12.6,' yr')
 1260 format (22x,'tauD0 = ',1pe12.6,' yr, tauD = ',1pe12.6,' yr, ',
     &     'tauaccr = ',1pe12.6,' yr (min at # ',i4,')')
 1270 format (22x,'LDaccr = ',1pe12.6,' Laccr = ',1pe12.6)
 1300 format (5x,'at tau = ',f7.5,' (numerical surface):')
 1400 format (7x,'L* = ',f10.2,'  R* = ',f7.2,'  T(N) = ',f7.0,
     &     '  ro(N) = ',1pe11.5,'  log[g(N)] = ',1pe11.4,'  M* = ',
     &     0pf14.10)
 1500 format (5x,'at tau = ',f7.5,' (photosphere):')
 1600 format (7x,'L = ',f10.2,'  Reff = ',f7.2,'  Teff = ',f7.0,
     &     '  roeff = ',1pe11.5,'  log(geff) = ',1pe11.4,'  M = ',
     &     0pf14.10,/)
 1800 format (5x,'surface: neff = ',i4,'  and  n = ',i4,/,
     &     14x,'mass loss from shell',2x,i4,/)
 1900 format (5x,'# of convective zones: ',i2,' [alphac = ',f7.4,
     &     ']')
 2000 format (5x,'Schwar.: ',i4,'-->',i4,',  Oversh.: ',i4,'-->',i4,
     &     ' (',f5.2,',',f5.2,')')
 2100 format (/,5x,'free-fall time-scale',9x,' = ',1pe11.5,' yr',/,5x,
     &     'Kelvin-Helmholtz time-scale',2x,' = ',1pe11.5,' yr',/,5x,
     &     'thermal adjustment time-scale = ',1pe11.5,' yr',/,5x,
     &     'nuclear time-scale',11x,' = ',1pe11.5,' yr',/)
 2200 format (5x,'mean turnover time of the conv. core = ',1pe11.5,
     &     ' yr',/)
 2300 format (5x,'mean turnover time of the conv. env. = ',1pe11.5,
     &     ' yr',/)
 2400 format (5x,'H-burning Luminosity = ',f14.2,/,5x,
     &     'He-burning Luminosity = ',f14.2)
 2500 format (5x,'He-burning shell mean irradiation: total = ',
     &     1pe9.3,' (in ',1pe10.4,' sm)  -  pulse = ',1pe9.3)
 2600 format (/,5x,'global energetics:',/,5x,'epot = ',1pe11.4,
     &     ', ekin = ',1pe11.4,', eint = ',1pe11.4,
     &     ', edyn = epot+ekin+eint = ',1pe11.4,/,5x,'ephot = ',
     &     1pe11.4,', eneu = ',1pe11.4,', enuc = ',1pe11.4,
     &     ', etot = edyn+ephot+eneu = ',1pe11.4,' =? enuc')
 2700 format (//,3x,'#',7x,'r/R*',8x,'T',9x,'rho',9x,'P',7x,'beta',4x,
     &     'eta',5x,'lnf',9x,'s',11x,'l',10x,'u',12x,'m/M*',/)
 2800 format (1x,i4,a1,1x,0pf11.9,3(1x,1pe10.4),1x,0pf6.4,1x,0pf7.3,1x,
     &     0pf8.3,1x,1pe11.4,1x,1pe11.4,1x,1pe10.3,1x,0pf15.13)
 2900 format (//,3x,'#',6x,'tau',6x,'kappa',7x,'mu',5x,'mue',5x,
     &     'adgrad',3x,'mugrad',5x,'radgrad',6x,'grad',7x,'soundv',
     &     7x,'cp',6x,'gamma',5x,'eint',6x,'pvisc',/)
 3000 format (1x,i4,a1,1x,0pf9.5,1x,1pe10.4,1x,0pf7.4,1x,1pe9.3,1x,
     &     0pf7.5,1x,1pe10.3,1x,1pe11.4,1x,1pe11.4,1x,1pe10.4,1x,
     &     1pe10.4,1x,0pf7.4,1x,1pe10.4,1x,1pe10.3)
 3100 format (//,3x,'#',4x,'tauconv',5x,'sconv',6x,'lconv/l',3x,
     &     'lrad/l',7x,'enuc',4x,'enuc-enu',5x,'egrav',5x,'hydrat',
     &     6x,'Dconv',7x,'dmicro',5x,'vdmicro',4x,'enupla',/)
 3200 format (1x,i4,a1,1x,1pe10.4,1x,1pe10.4,1x,1pe10.3,1x,
     &     1pe9.3,1x,1pe10.4,1x,1pe11.4,1x,1pe11.4,1x,1pe9.3,
     &     4(1x,1pe11.4))
 3231 format (//,3x,'#',4x,'omega',7x,'facc',7x,'macc',8x,'eacc',8x,
     &     'eshr',7x,'tacc',/)
 3261 format (1x,i4,a1,1x,1pe11.5,1x,0pf9.7,1x,1pe12.5,1x,1pe11.4,1x,
     &     1pe11.4,1x,1pe11.5)
 3232 format (//,3x,'#',4x,'facc',7x,'macc',8x,'eacc',8x,'eshr',8x,
     &     'tacc',/)
 3262 format (1x,i4,a1,1x,0pf9.7,1x,1pe12.5,1x,1pe11.4,1x,
     &     1pe11.4,1x,1pe11.5)
 3300 format (//,1x,'kinetic to convective fluxes ratio:')
 3330 format (//,3x,'#',6x,'omega',8x,'urs',8x,'dift',8x,'dife',8x,
     &     'Dhd',8x,'theta',8x,'lambda',8x,'abmuj',8x,'Dshear',8x,'Nu'/)
 3360 format (1x,i4,a1,1x,1pe11.5,1x,0pf11.4,1x,1pe11.4,1x,1pe10.4)
 3500 format (//,43x,'chemical abundance profiles:')
 3600 format (//,3x,'#',2x,10(3x,a5,3x),/)
 3700 format (1x,i4,a1,10(1x,1pe10.4))
 3800 format (/,43x,'mass lost (in sm) for each nuclide:',/)
 3900 format (6x,10(1x,1pe10.4))
 4000 format (//,3x,'#',2x,3(3x,a5,3x),4x,'Sum X',/)
 4100 format (1x,i4,a1,3(1x,1pe10.4),1x,1pe11.4)
 4200 format (6x,3(1x,1pe10.4))


c..   sm files formats

 1205 format ('#  nsh',8x,'MM',13x,'LL',10x,'RR',11x,'times',13x,
     &     'iacc irot npha  shell  nxsp   modl   tstep',7x,
     &     'nseq version ')
 1210 format (1x,i5,2x,f14.10,2x,1pe11.4,2x,1pe10.3,2x,1pe21.15,3(4x,i1)
     &     ,2x,i5,2x,i4,2x,i6,2x,1pe9.3,1x,i6,3x,i4,5x,'AGB')

 1211 format ('# nsh yzi',9x,'r',11x,'T',9x,'rho',11x,'P',7x,'beta',4x,
     &     'eta',6x,'lnf',8x,'s',11x,'Lr',9x,'u',11x,'xmr',14x,'dm')
 2810 format (1x,i4,2x,i1,1x,0pf14.9,3(1x,1pe11.4),1x,0pf6.4,1x,0pf7.3,
     &     1x,0pf8.3,1x,1pe11.4,1x,1pe11.4,1x,1pe10.3,2(1x,0pf15.13))

 1212 format ('# nsh',5x,'tau',7x,'kap',8x,'mu',7x,'mue',5x,'abadd',5x,
     &     'abmu',7x,'abrad',8x,'abla',8x,'cs',9x,'cp',6x,'gamma1',4x,
     &     'eint',7x,'pvisc')
 3010 format (1x,i4,1x,1pe10.3,1x,1pe10.4,1x,0pf7.4,1x,1pe10.4,1x,
     &     0pf7.5,1x,1pe10.3,1x,1pe11.4,1x,1pe11.4,1x,1pe10.4,1x,
     &     1pe10.4,1x,0pf7.4,1x,1pe10.4,1x,1pe10.3)

 1213 format ('# nsh',3x,'tconv',6x,'Vconv',6x,'rfconv',5x,'rfrad',6x,
     &     'enucl',7x,'enupla',7x,'egrav',5x,'hydrat',5x,'Dconv',7x,
     &     'Dmicro',6x,'Vmicro',6x,'enunucl')
 3210 format (1x,i4,1x,1pe10.4,1x,1pe10.4,1x,1pe10.3,1x,
     &     1pe9.3,1x,1pe11.4,1x,1pe11.4,1x,1pe11.4,1x,1pe9.3,
     &     4(1x,1pe11.4))

 1221 format ('# nsh',4x,'omega',7x,'facc',8x,'macc',8x,'eacc',8x,
     &     'eshr',7x,'tacc')
 3461 format (1x,i4,1x,1pe11.5,1x,0pf9.7,1x,1pe12.5,1x,1pe11.4,1x,
     &     1pe11.4,1x,1pe11.5)

 1222 format ('# nsh',4x,'facc',7x,'macc',7x,'eacc',7x,
     &     'eshr',7x,'tacc')
 3462 format (1x,i4,1x,0pf9.7,1x,1pe12.5,1x,1pe11.4,1x,
     &     1pe11.4,1x,1pe11.5)

 1223 format ('# nsh',3x,'omega',9x,'Ucirc',8x,'Dshear',7x,'Dcirc',
     &     8x,'Dtot',10x,'Dh',9 x,'theta',8x,'xlambda',7x,
     &     'abmuj',7x,'Nuturb',8x,'Numol',9x, 'Kt')
 3460 format (1x,i4,1x,1pe11.5,1x,1pe12.5,1x,1pe12.5,1x,
     &     1pe12.5,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5,
     &     1x,1pe12.5,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5)

 1224 format ('# nsh',5x,'comp',9x,'heat',9x,'work',9x,'etherm',6x,
     &     'dlumdt',7x,'denucdt')
 3465 format (1x,i4,6(1x,1pe12.5))

c..   chemical profiles
c..   new version >= 2.10
 1214 format ('# nsh',7x,'n',11x,'H1',11x,'H2',10x,'He3',10x,'He4',10x,
     &     'Li6',10x,'Li7',10x,'Be7',10x,'B8',11x,'Be9')
 1215 format ('# nsh',5x,'B10',10x,'B11',10x,'C12',10x,'C13',10x,'N13',
     &     10x,'C14',10x,'N14',10x,'N15',10x,'O15',10x,'O16')
 1216 format ('# nsh',6x,'O17',10x,'O18',10x,'F18',10x,'F19',10x,'F20',
     &     9x,'Ne20',9x,'Ne21',9x,'Ne22',9x,'Na22',9x,'Ne23')
 1217 format ('# nsh',5x,'Na23',9x,'Na24',9x,'Mg24',9x,'Na25',9x,
     &     'Mg25',9x,'Mg26',8x,'Alm26',8x,'Alg26',9x,'Mg27',9x,
     &     'Al27')
 1218 format ('# nsh',5x,'Si28',9x,'Si29',9x,'Si30',10x,'P31',10x,
     &     'S32',10x,'S33',10x,'S34',10x,'S35',9x,'Cl35',10x,'S36')
 1219 format ('# nsh',4x,'Cl36',7x,'Cl37',6x,'heavy',9x,'sumX')
 3710 format (1x,i4,10(1x,1pe12.6))
 4110 format (1x,i4,3(1x,1pe10.4),1x,1pe11.4)

c..   old version < 2.10
 2214 format ('# nsh',5x,'n',10x,'H1',9x,'H2',8x,'He3',8x,'He4',8x,
     &     'Li6',8x,'Li7',8x,'Be7',9x,'B8',8x,'Be9')
 2215 format ('# nsh',4x,'B10',8x,'B11',8x,'C11',8x,'C12',8x,'C13',8x,
     &     'N13',8x,'C14',8x,'N14',8x,'N15',8x,'O15')
 2216 format ('# nsh',4x,'O16',8x,'O17',8x,'O18',8x,'F18',8x,'F19',7x,
     &     'Ne20',7x,'Ne21',7x,'Ne22',7x,'Na22',7x,'Na23')
 2217 format ('# nsh',4x,'Mg24',7x,'Mg25',7x,'Mg26',6x,'Alm26',6x,
     &     'Alg26',7x,'Al27',7x,'Si28',7x,'Si29',7x,'Si30',7x,'Si31')
 2218 format ('# nsh',5x,'P31',7x,'Si32',8x,'P32',8x,'S32',8x,'P33',
     &     8x,'S33',8x,'S34',8x,'S35',7x,'Cl35',8x,'S36')


c..   not used

 4300 format (//,1x,i2,' shells with a maximum in nuclear energy ',
     &     'production greater than ',1pe7.1,' erg/g/sec:',
     &     //,2x,'#',7x,'m/M*',9x,'enuc',/)
 4400 format (1x,i4,1x,f12.10,1x,1pe12.6)
 4500 format (//,1x,'nuclear fluxes at the burning shell #',i4,':',/)
 4600 format (7('#',i4,': ',1pe11.5))
 4700 format (6('#',i4,': ',1pe11.5))
 4800 format (//,1x,'energy contributions to the total luminosity:',/,
     &     1x,'gravitation: ',f13.5,'  -  nuclear : ',f13.5)
 4900 format (//,1x,'global nuclear reaction contributions to the ',
     &     'total nuclear energy production:',/)
 5000 format (1x,a37,': ',f9.6,' %')
 5100 format (//,1x,'nuclear reaction contributions to the nuclear ',
     &     'energy production by the burning shell # ',i4,
     &     ' [enuc = ',1pe12.6,']:',/)
 5200 format (//,'*-*-*-*-*-*-*-*-  End of Listing  -*-*-*-*-*-*-*-*',
     &     //)

      end


************************************************************************

      SUBROUTINE SORT (xxor,react,react2)

************************************************************************
* Sort the array enc into ascending numerical order, with              *
* corresponding re-ordering of the string react2 (initially = react)   *
************************************************************************

      implicit none

      include 'evolpar.chem'

      character react*37,react2*37,rreact*37

      integer i,j,l

      double precision xxor
      double precision xxori

      dimension xxor(nreac),react(nreac),react2(nreac)

      do l = 1,nreac
         react2(l) = react(l)
      enddo

      do j = 2,nreac
         xxori = xxor(j)
         rreact = react2(j)
         do i = j-1,1,-1
            if (xxor(i).le.xxori) goto 10
            xxor(i+1) = xxor(i)
            react2(i+1) = react2(i)
         enddo
         i = 0
 10      xxor(i+1) = xxori
         react2(i+1) = rreact
      enddo

      return
      end


************************************************************************

      INTEGER FUNCTION LENCI (chain)

************************************************************************
*   Last non-blank character position in string a                      *
************************************************************************

      character*(*) chain

      integer i,n

      n = len(chain)
      do i = n,1,-1
         if (chain(i:i).eq.'_') then
            lenci = i
            return
         endif
      enddo

      lenci = 1

      return
      end


************************************************************************

      INTEGER FUNCTION LENCA (chain)

************************************************************************
*   Last non-blank character position in string a                      *
************************************************************************

      character*(*) chain

 
      integer i,n

      n = len(chain)
      do i = 1,n
         if (chain(i:i).eq.'i') then
            lenca = i-1
            return
         endif
      enddo

      lenca = 0

      return
      end


************************************************************************

      SUBROUTINE MAKE_CHAR (number,chain)

************************************************************************

      character*(*) chain
      character*1 car
      integer i,a,number
      double precision fact,nn

      external car

      nn = dble(number)
      do i = 1,6
         fact = 10.**(i-1)/1.d5
         a = int(nn*fact)
         chain(i:i) = car(a)
         nn = nn - dble(a)/fact
      enddo

      return
      end


************************************************************************

      SUBROUTINE MAKE_NUMBER (chain,number)

************************************************************************

      character*(*) chain
      character*1 car
      integer i,a,num,number
      double precision fact,nn

      external num

      number = 0
      do i = 1,4
         fact = 10**(4-i)
         nn = num(chain(i:i))
         number = number + fact*nn
      enddo

      return
      end


************************************************************************

         CHARACTER FUNCTION CAR (a)

************************************************************************

         integer a

         if (a.eq.0) car = '0'
         if (a.eq.1) car = '1'
         if (a.eq.2) car = '2'
         if (a.eq.3) car = '3'
         if (a.eq.4) car = '4'
         if (a.eq.5) car = '5'
         if (a.eq.6) car = '6'
         if (a.eq.7) car = '7'
         if (a.eq.8) car = '8'
         if (a.eq.9) car = '9'

         return
         end

************************************************************************

         INTEGER FUNCTION NUM (a)

************************************************************************

         character*1 a

         num = 0
         if (a.eq.'0') num = 0
         if (a.eq.'1') num = 1
         if (a.eq.'2') num = 2
         if (a.eq.'3') num = 3
         if (a.eq.'4') num = 4
         if (a.eq.'5') num = 5
         if (a.eq.'6') num = 6
         if (a.eq.'7') num = 7
         if (a.eq.'8') num = 8
         if (a.eq.'9') num = 9

         return
         end
