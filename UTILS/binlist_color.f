
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
      include 'evolcom.eos'
      include 'evolcom.grad'
      include 'evolcom.ion'
      include 'evolcom.lum'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.opa'
      include 'evolcom.rot'
      include 'evolcom.shk'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.transp'
      include 'evolcom.var'

      character*50 name
      character*100 dirmongo,diresults,dirbatch,opalib
      character cmodel*6
      character*37 react2
      character*1 rcrzc

      integer icourant,itmax,ishockt,ishockb,imack
      integer imod,error
      integer lenci,lenca,sm
      integer nratsh,nsequence
      integer ienvb,ienvt,nrat
      integer i,ic,j,k,kl,l,ip
      integer crzi
      integer klenv,klpulse,klcore,itop,ibase
      integer inshock,outshock

      double precision dturb,dmic,vdmic,dtinv,vvro,cd0,Dherw
      double precision shlimplt,trdiff
      double precision v,enc,dumarray,sr
      double precision totmcor,totmenv,ebind
      double precision t9,eng,eerho,eerhocap
      double precision version
      double precision Dhold
      double precision xNt,xNu,xNr,geffom
      double precision abmurj,abmuj,mul
      double precision xspr
      double precision epot,eint,ebindenv,epotenv,eintenv
      double precision comp,heat,work,etherm,vve
      double precision renucl
      double precision Dtot,fover
      double precision khim,alphasc
      double precision geffs,tnucg,tff
      double precision mack,Mackmax,Mackmin,tshock
      double precision dtcour,lmax,tmax,enucmax

      double precision dift,dife,rpvisc,rsconv,dlumdt,denucdt
c      double precision summe,rpvisc,rsconv

      logical shockdetect

      common /difcirc/ dift(nsh),dife(nsh)
      common /coefdiff/ dturb(nsh),dmic(nsh),vdmic(nsh)
      common /nnnk/ xNt(nsh),xNu(nsh),xNr(nsh),geffom(nsh)
      common /calcDh/ Dhold(nsh)
      common /gradmuj/ abmurj(nsh),abmuj(nsh),mul(nsh)
      common /overshoot/ klcore,klenv,klpulse
      common /hydrodyn/ lmax,tmax,Mackmax,enucmax,itmax,ishockb,ishockt

      dimension v(nreac),enc(nreac),nratsh(nsh),react2(nreac),
     &     rpvisc(nsh),rsconv(nsh),vvro(nsh)
c      dimension xKt(nsh),dift(nsh),dife(nsh),Dhold(nsh),vtheta(nsh)
      dimension Dtot(nsh),Dherw(nsh)
      dimension mack(nsh),sr(nsh)

      dimension crzi(nsh),rcrzc(nsh),vve(nsh)
      dimension xspr(nsh,nsp)
      dimension comp(nsh),heat(nsh),work(nsh),etherm(nsh),dlumdt(nsh),
     &     denucdt(nsh),renucl(nsh)

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
      mprot = 1.6726231d-24
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

      write (*,*)'filename : ',trim(name),', imod = ',imod
      i=lenci(name)
      cmodel = name(i+1:i+5)
      call make_number (cmodel,nsequence)
      call getenv('DIR_SMONGO',dirmongo)
      call getenv('DIR_RESULTS',diresults)
      call getenv('DIR_BATCH',dirbatch)

      if (imod.eq.9) then
         open (unit = 90,file = trim(diresults) // trim(name),
     &        form = 'unformatted')
         open (unit = 99,file = trim(diresults) // '/starevol.par',
     &        status = 'unknown')
         read (99,2) version
 2       format (/,43x,f4.2) 
         rewind(99)
      else
         open (unit = 90,file = trim(diresults) // trim(name) // 'b',
     &        form = 'unformatted')
         open (unit = 99,file = trim(diresults) // trim(name) // 'd',
     &        status = 'unknown')

         read (99,3) version
 3       format (/////,43x,0pf4.2)
         backspace(99)
         backspace(99)
      endif
      write (*,22) version
 22   format(' STAREVOL version : ',f4.2)
      open (unit = 92,file = '/tmp/nextini.bin',form = 'unformatted',
     &     status = 'unknown')

*_______________________________
***   reading of the binary file
*-------------------------------
      
c      if (version.lt.2.05d0) then
c         call rmodpar_2_03
c      else
         call rmodpar
c      endif

      close(99)
      mixopt = .false.
      lnucl = .true.
      nuclopt = 'u'
      diffusconv = .true.
      if (idiffty.ge.8.and.idiffty.le.15) then
         if (imod.eq.9) then
            open (unit = 93,file = trim(diresults) // 
     &           name(1:lenca(name)) // 'ang.bin',form = 'unformatted',
     &           status = 'unknown')
         else
            open (unit = 93,file = trim(diresults) // trim(name) // 'a',
     &           form = 'unformatted',status = 'unknown')
         endif
      endif

      if (version.ge.2.1d0) then
         open (unit = 99, file = trim(dirbatch) // 
     &        '/starevolnuc_evol.par',status = 'unknown')
      else
         open (unit = 99, file = trim(dirbatch) // 
     &        '/starevolnuc_old_evol.par',status = 'unknown')
      endif
      read (99,5)
 5    format(1x,//)


      call rinimod (90,93,92,95,0)
      print *,'Total number of shells : ',nmod

      call eos (0,error)

      no = 0

*______________________________________________
***   generation of ascii files for Super Mongo
*----------------------------------------------

      call make_char (model,cmodel)
      open(unit=10,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p0',status = 'unknown')
      open(unit=14,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p4',status = 'unknown')
      open(unit=15,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p5',status = 'unknown')
      open(unit=16,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p6',status = 'unknown')
      open(unit=17,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p7',status = 'unknown')
      open(unit=18,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p8',status = 'unknown')
      open(unit=19,file=trim(dirmongo) // 'DATA/' // 
     &     name(1:lenci(name)) //cmodel// '.p9',status = 'unknown')

      write (10,1205)
      write (10,1210) nmod,totm,lum(nmod)/lsun,r(nmod)/rsun,time
     &     ,iaccr,irotbin,nphase,nmod,nsp-1,model,dtn,nsequence
     &     ,int(version*1.d2)
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
      
      close(10)
      close(14)
      close(15)
      close(16)
      close(17)
      close(18)
      close(19)

      write (*,*)
      write (*,*) 'Files saved as : ',trim(dirmongo) // 'DATA/' //
     &     name(1:lenci(name)) //cmodel
      write (*,*) 'readsi ',name(1:lenci(name)) //cmodel
      stop 'sm files generated !'

      stop 'binlist'

 877  format(/,'Arg : #1 model number (in file)',/,6x,
     &     '#2 output : 0 ascii (.lst)',/,19x,'1  smongo (.p0-9)',/,19x,
     &     '2  IDL (.l)',/,6x,'#3 filename',/,3x,/,
     &     '  --> ex: 5 1 m1.0z0_0005b')
 997  format (' Not so many models saved in that file!',/)
 998  format (' Not so many listing models!',/,
     &     ' The last model has been selected.')
 999  format (i1,1x,i1,1x,a50)

c..   sm files formats

 1205 format ('#  nsh',8x,'MM',13x,'LL',10x,'RR',11x,'times',13x,
     &     'iacc irot npha  shell  nxsp   modl   tstep',7x,
     &     'nseq version ')
 1210 format (1x,i5,2x,f14.10,2x,1pe11.4,2x,1pe10.3,2x,1pe21.15,3(4x,i1)
     &     ,2x,i5,2x,i4,2x,i6,2x,1pe9.3,1x,i6,3x,i4,5x,'AGB')

 1211 format ('# nsh yzi',9x,'r',11x,'T',9x,'rho',11x,'P',7x,'beta',4x,
     &     'eta',6x,'lnf',8x,'s',11x,'Lr',9x,'u',11x,'xmr',14x,'accel')
 2810 format (1x,i4,2x,i1,1x,0pf14.9,3(1x,1pe11.4),1x,0pf6.4,1x,0pf8.3,
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
 3461 format (1x,i4,1x,1pe12.5,1x,0pf9.7,1x,1pe12.5,1x,1pe11.4,1x,
     &     1pe11.4,1x,1pe11.5)

 1222 format ('# nsh',4x,'facc',7x,'macc',7x,'eacc',7x,
     &     'eshr',7x,'tacc')
 3462 format (1x,i4,1x,0pf9.7,1x,1pe12.5,1x,1pe11.4,1x,
     &     1pe11.4,1x,1pe11.5)

 1223 format ('# nsh',3x,'omega',9x,'Ucirc',8x,'Dshear',7x,'Dcirc',
     &     8x,'Dtot',10x,'Dh',9 x,'theta',8x,'xlambda',7x,
     &     'abmuj',7x,'Nuturb',8x,'Numol',9x, 'Kt')
 3460 format (1x,i4,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5,1x,
     &     1pe12.5,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5,
     &     1x,1pe12.5,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5)

 1224 format ('# nsh',5x,'comp',9x,'heat',9x,'work',9x,'etherm',6x,
     &     'dlumdt',7x,'denucdt',7x,'eloc',8x,'egconv')
 3465 format (1x,i4,8(1x,1pe12.5))

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
      integer i,num,number
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
