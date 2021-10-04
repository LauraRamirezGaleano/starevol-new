      SUBROUTINE rininet

************************************************************************
* Read nuclear network and isotopes properties                         *
* 09/03: All nuclei with mass fractions greater than 1d-15 diffuse     *
*                                                                      *
* $LastChangedDate:: 2016-05-17 18:06:34 +0200 (Mar, 17 mai 2016)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 83                                                          $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.data'
      include 'evolcom.nuc'
      include 'evolcom.spec'

      integer normx,na,nb,nc,nd0,ia,ib,ic,id,iznuc,ianuc
      integer idum,kbeta,kbetaec,kneut
      integer ispec(0:nbz,2*nbz+3)
      integer i,j,k,l

      double precision sum1

      logical diffelem

      common /logicdiff/ diffelem(nsp)

      integer nz,lirn,ina,jna
      parameter (nz = 370,lirn=nsp)

c      double precision v(nreac),y(nsp),eqd(nsp,3*nsp+1)
c      common /nuc_matrix/ jna(nz),ina(lirn)

*___________________________________________
***   determination of numerical constraints
*-------------------------------------------

      ymin = 1.d-50
      ytmin = 1.d-12*tolnuc
      dnuc = 1.d-10*tolnuc

*_________________________________________________
***   determination of the nuclear matrix features
*-------------------------------------------------

      fact(0) = 1.d0
      fact(1) = 1.d0
      fact(2) = 2.d0
      fact(3) = 6.d0
      fact(4) = 24.d0
      fact(5) = 120.d0
      fact(6) = 720.d0

      do l = 1,2*nbz+3
         do j = 0,nbz
            ispec(j,l) = -1
         enddo
      enddo

*__________________
***   read nuclides
*------------------

      ial26g = -1
      ial26m = -1
      ndecay = 0
      do i = 1,nsp
         read (99,1000) idum,elem(i),anuc(i),znuc(i),xspacc(i),istate(i)
         if (istate(i).eq.'F'.and.i.lt.nis) then
            ndecay = ndecay+1
            kdecay(ndecay) = i
         endif
         if (i.lt.nsp) ispec(int(znuc(i)),int(anuc(i))) = i
         if (elem(i).eq.'AL26M'.or.elem(i).eq.'AL26m') ial26m = i
         if (elem(i).eq.'AL26G'.or.elem(i).eq.'AL26g') ial26g = i
      enddo
      if (ial26m.eq.-1) then
          write (nout,*) 'ial26m not defined'
         stop 'STOP rininet'
      endif
      if (ial26g.eq.-1) then
          write (nout,*) 'ial26g not defined'
         stop 'STOP rininet'
      endif
      if (ndecay.ne.10) then
         write (nout,*) ndecay
         stop 'rininet : wrong starevol.par, check istate'
      endif

*_______________________________________________________
***   initialize molecule index for opacity calculations
*-------------------------------------------------------

      idxmol(1,1) = iih2
      idxmol(1,2) = iih
      idxmol(1,3) = iih
      idxmol(2,1) = iih2o
      idxmol(2,2) = iih
      idxmol(2,3) = iio
      idxmol(3,1) = iioh
      idxmol(3,2) = iih
      idxmol(3,3) = iio
      idxmol(4,1) = iico
      idxmol(4,2) = iic
      idxmol(4,3) = iio
      idxmol(5,1) = iicn
      idxmol(5,2) = iic
      idxmol(5,3) = iin
      idxmol(6,1) = iic2
      idxmol(6,2) = iic
      idxmol(6,3) = iic
      idxmol(7,1) = iin2
      idxmol(7,2) = iin
      idxmol(7,3) = iin

*_____________________________
***   initialize species index
*-----------------------------

      ineutron = ispec(0,1)
      ih1 = ispec(1,1)
      ih2 = ispec(1,2)
      ihe3 = ispec(2,3)
      ihe4 = ispec(2,4)
      ili6 = ispec(3,6)
      ili7 = ispec(3,7)
      ibe7 = ispec(4,7)
      ibe9 = ispec(4,9)
      ib8 = ispec(5,8)
      ib11 = ispec(5,11)
      ic12 = ispec(6,12)
      ic13 = ispec(6,13)
      ic14 = ispec(6,14)
      in13 = ispec(7,13)
      in14 = ispec(7,14)
      in15 = ispec(7,15)
      io15 = ispec(8,15)
      io16 = ispec(8,16)
      io17 = ispec(8,17)
      io18 = ispec(8,18)
      if18 = ispec(9,18)
      if19 = ispec(9,19)
      if20 = ispec(9,20)
      ine20 = ispec(10,20)
      ine23 = ispec(10,23)
      ina22 = ispec(11,22)
      ina23 = ispec(11,23)
      ina24 = ispec(11,24)
      ina25 = ispec(11,25)
      img24 = ispec(12,24)
      img25 = ispec(12,25)
      img27 = ispec(12,27)
      ial27 = ispec(13,27)
      isi28 = ispec(14,28)
      ip31 = ispec(15,31)
      is32 = ispec(16,32)
      is35 = ispec(16,35)
      icl35 = ispec(17,35)
      icl36 = ispec(17,36)
c     Ajout pour diffusion Thoul                                 ! modif Thoul
      ine21 = ispec(10,21)
      ine22 = ispec(10,22)
      img26 = ispec(12,26)
      ib10 = ispec(5,10)
      
      zsol = 1.d0-xspsol(ih1)-xspsol(ih2)-xspsol(ihe3)-xspsol(ihe4)


*______________________________________________________
***   initialize desintegration half-lifes (in seconds)
*------------------------------------------------------

      do i = 1,nsp
         tdecay(i) = 1.d99
      enddo
      tdecay(ibe7)   = 4650048.d0
      tdecay(ib8)    = 0.770d0
      tdecay(in13)   = 598.2d0
      tdecay(ic14)   = 1.803478d11
      tdecay(io15)   = 122.2d0
      tdecay(if18)   = 6588.d0
      tdecay(if20)   = 11.d0
      tdecay(ina22)  = 82205792.d0
      tdecay(ine23)  = 37.2d0
      tdecay(ina24)  = 5.385d4
      tdecay(ina25)  = 59.3d0
      tdecay(ial26m) = 6.3452d0
      tdecay(ial26g) = 2.272098d13
      tdecay(img27)  = 567.d0
      tdecay(is35)   = 7534080.d0
      tdecay(icl36)  = 9.498634d14


*___________________________________________________
***   initialize index idxspc and idxsps for outputs
*---------------------------------------------------

      i = 0
      do k = 1,nsp-2
         iznuc = int(znuc(k))
         ianuc = int(anuc(k))
         if (.not.((iznuc.eq.5.and.ianuc.eq.8).or.k.eq.in13.or.
     &        (iznuc.eq.9.and.(ianuc.eq.18.or.ianuc.eq.20)).or.
     &        (iznuc.eq.10.and.ianuc.eq.23).or.(iznuc.eq.11.and.
     &        (ianuc.eq.22.or.ianuc.eq.24.or.ianuc.eq.25)).or.
     &        (iznuc.eq.12.and.ianuc.eq.27))) then
            i = i+1
            if (i.gt.nprint) then
               write (nout,*) 'nprint (',nprint,') must be set to',i,
     &              ' in evolpar.star'
               stop ' rininet : idxsp[c,s] dimension problem'
            endif
            idxspc(i) = k
            idxsps(i) = k
         endif
      enddo
c..   index for .c5 file
      idxspci(1) = in13
      idxspci(2) = if18
      idxspci(3) = if20
      idxspci(4) = ine23
      idxspci(5) = ina22
      idxspci(6) = ina24
      idxspci(7) = ina25
      idxspci(8) = img25
      idxspci(9) = img27

      do i = 1,nsp
         diffelem(i) = .true.
      enddo

*_______________________________________
***   accreted abundance renormalization
*---------------------------------------

      if (iaccr.ge.1) then
         normx = 1
         sum1 = 0.d0
         do l = 1,nsp
            if (xspacc(l).gt.xspacc(normx)) normx = l
            sum1 = sum1+xspacc(l)
         enddo
         if (abs(sum1-1.d0).gt.1.d-6) then
            write (nout,*) 'composition of accreted matter not ',
     &           'normalized : sumxsp-1 = ',sum1-1.d0
            stop 'rininet'
         endif
         xspacc(normx) = xspacc(normx)+(1.d0-sum1)

         muiacc = 0.d0
         do l = 1,nis
            muiacc = muiacc+xspacc(l)/anuc(l)
         enddo
      endif

      read (99,1200)

*________________________________________________________________
***   read nuclear reaction network
*  reaction : kna X(kia) + knb X(kib) --> knd X(kid) + knc X(kic)
***      ex :  1    C12  +  1    He   -->  0  gamma  +  1   O16
*----------------------------------------------------------------

      iddn = -1
      ippd = -1
      inpg = -1
      i3a = -1
      icag = -1
      ioag = -1
      iccga = -1
      iccgp = -1
      ioca = -1
      iocp = -1
      iocn = -1
      iooga = -1
      ioogp = -1
      inega = -1
      ippg = -1
      ipdg = -1
      i2he3 = -1
      ihe3ag = -1
      ili6pa = -1
      ili7pa = -1
      ili7d = -1
      ibe7beta = -1
      ibe7pg = -1
      ib8beta = -1
      ib11pa = -1
      icpg = -1
      ic13pg = -1
      in13beta = -1
      in14pg = -1
      in15pa = -1
      in15pg = -1
      io15beta = -1
      io16pg = -1
      io17pa = -1
      ine20ag = -1
      img24ag = -1
      isi28ag = -1

      kbeta = 0
      kbetaec = 0
      kneut = 0

      do l = 1,nreac
         read (99,1100) idum,na,nb,nd0,nc,qi(l),ia,ib,id,ic,qinu(l)
***   index of beta decay reactions
         if (ib.eq.0.and.id.eq.0) then
            kbeta = kbeta+1
            if (kbeta.gt.nbeta) goto 10
            jbeta(kbeta) = l
            if (l.gt.nre) kbetaec = kbetaec+1
         endif
***   index of neutron capture reactions
         if (ib.eq.1) then
            kneut = kneut+1
            if (kneut.gt.nneut) goto 10
            jneut(kneut) = l
         endif

         if (ib.eq.0) ib = nsp
         if (id.eq.0) id = nsp
         if (ia.eq.ih2.and.na.eq.2.and.ic.eq.ihe3) iddn = l
         if (ia.eq.ili7.and.ib.eq.ih2.and.ic.eq.ihe4) ili7d = l
         if (ia.eq.ih1.and.na.eq.2.and.ic.eq.ih2) ippd = l
         if (ia.eq.in14.and.ib.eq.ih1.and.ic.eq.io15) inpg = l
         if (ia.eq.ic12.and.ib.eq.ih1.and.ic.eq.in13) icpg = l
         if (ia.eq.ihe4.and.ic.eq.ic12) i3a = l
         if (ia.eq.ic12.and.ib.eq.ihe4.and.ic.eq.io16) icag = l
         if (ia.eq.io16.and.ib.eq.ihe4.and.ic.eq.ine20) ioag = l
         if (ia.eq.ic12.and.id.eq.ihe4.and.ic.eq.ine20) iccga = l
         if (ia.eq.ic12.and.id.eq.ih1.and.ic.eq.ina23) iccgp = l
         if (ia.eq.io16.and.id.eq.ihe4.and.ic.eq.isi28) iooga = l
         if (ia.eq.io16.and.id.eq.ih1.and.ic.eq.ip31) ioogp = l
         if (ia.eq.io16.and.id.eq.ihe4.and.ic.eq.img24) ioca = l
         if (ia.eq.io16.and.id.eq.ih1.and.ic.eq.ial27) iocp = l
         if (ia.eq.io16.and.id.eq.1.and.ic.eq.ial27) iocn = l
         if (ia.eq.ine20.and.id.eq.ihe4.and.ic.eq.io16) inega = l
         if (ia.eq.ih1.and.na.eq.2.and.ic.eq.ih2) ippg = l
         if (ia.eq.ih2.and.ib.eq.ih1.and.ic.eq.ihe3) ipdg = l
         if (ia.eq.ihe3.and.na.eq.2.and.ic.eq.ihe4) i2he3 = l
         if (ia.eq.ihe4.and.ib.eq.ihe3.and.ic.eq.ibe7) ihe3ag = l
         if (ia.eq.ili6.and.ib.eq.ih1.and.ic.eq.ihe4) ili6pa = l
         if (ia.eq.ili7.and.ib.eq.ih1.and.ic.eq.ihe4) ili7pa = l
         if (ia.eq.ibe7.and.ib.eq.nsp.and.ic.eq.ili7) ibe7beta = l
         if (ia.eq.ibe7.and.ib.eq.ih1.and.ic.eq.ib8) ibe7pg = l
         if (ia.eq.ib8.and.ib.eq.nsp.and.ic.eq.ihe4) ib8beta = l
         if (ia.eq.ib11.and.ib.eq.ih1.and.ic.eq.ihe4) ib11pa = l
         if (ia.eq.ic13.and.ib.eq.ih1.and.ic.eq.in14) ic13pg = l
         if (ia.eq.in13.and.ib.eq.nsp.and.ic.eq.ic13) in13beta = l
         if (ia.eq.in14.and.ib.eq.ih1.and.ic.eq.io15) in14pg = l
         if (ia.eq.in15.and.ib.eq.ih1.and.ic.eq.ic12) in15pa = l
         if (ia.eq.in15.and.ib.eq.ih1.and.ic.eq.io16) in15pg = l
         if (ia.eq.io15.and.ib.eq.nsp.and.ic.eq.in15) io15beta = l
         if (ia.eq.io16.and.ib.eq.ih1.and.ic.eq.io17) io16pg = l
         if (ia.eq.io17.and.ib.eq.ih1.and.id.eq.ihe4) io17pa = l
         if (ia.eq.ine20.and.ib.eq.ihe4.and.ic.eq.img24) ine20ag = l
         if (ia.eq.img24.and.ib.eq.ihe4.and.ic.eq.isi28) img24ag = l
         if (ia.eq.isi28.and.ib.eq.ihe4.and.ic.eq.is32) isi28ag = l

         k1(l) = ia
         k2(l) = na
         if (na.eq.0) then
            write (nout,*) ' Network error for reaction : ',idum
            write (nout,*) ' Wrong stoechiometric coefficients'
            stop 'rininet'
         endif
         k3(l) = ib
         k4(l) = nb
         k5(l) = ic
         k6(l) = nc
         k7(l) = id
         k8(l) = nd0
      enddo
      qdeut = qi(2)*econv
      if (kneut+1.ne.idneut) then
         write (nout,*) ' idneut = ',idneut,' not equal to kneut =',
     &        kneut+1
         write (nout,*) ' set idneut = ',kneut+1,' in evolpar.star'
         stop 'STOP rininet'
      endif
      if (ippd.eq.-1) then
         write (nout,*) 'ippd not defined'
         stop 'STOP rininet'
      endif
      if (iddn.eq.-1) then
         write (nout,*) 'iddn not defined'
         stop 'STOP rininet'
      endif
      if (inpg.eq.-1) then
         write (nout,*) 'inpg not defined'
         stop 'STOP rininet'
      endif
      if (i3a.eq.-1) then
         write (nout,*) 'i3a not defined'
         stop 'STOP rininet'
      endif
      if (icag.eq.-1) then
         write (nout,*) 'icag not defined'
         stop 'STOP rininet'
      endif
      if (ioag.eq.-1) then
         write (nout,*) 'ioag not defined'
         stop 'STOP rininet'
      endif
      if (iccga.eq.-1) then
         write (nout,*) 'iccga not defined'
         stop 'STOP rininet'
      endif
      if (iccgp.eq.-1) then
         write (nout,*) 'iccgp not defined'
         stop 'STOP rininet'
      endif
      if (iocp.eq.-1) then
         write (nout,*) 'iocp not defined'
         stop 'STOP rininet'
      endif
      if (ioca.eq.-1) then
         write (nout,*) 'ioca not defined'
         stop 'STOP rininet'
      endif
      if (iocn.eq.-1) then
         write (nout,*) 'iocn not defined'
         stop 'STOP rininet'
      endif
      if (inega.eq.-1) then
         write (nout,*) 'inega not defined'
         stop 'STOP rininet'
      endif
      if (ippg.eq.-1) then
         write (nout,*) 'ippg not defined'
         stop 'STOP rininet'
      endif
      if (ipdg.eq.-1) then
         write (nout,*) 'ipdg not defined'
         stop 'STOP rininet'
      endif
      if (i2he3.eq.-1) then
         write (nout,*) 'i2he3 not defined'
         stop 'STOP rininet'
      endif
      if (ihe3ag.eq.-1) then
         write (nout,*) 'ihe3ag not defined'
         stop 'STOP rininet'
      endif
      if (ili7d.eq.-1) then
         write (nout,*) 'ili7d not defined'
         stop 'STOP rininet'
      endif
      if (ili6pa.eq.-1) then
         write (nout,*) 'ili6pa not defined'
         stop 'STOP rininet'
      endif
      if (ili7pa.eq.-1) then
         write (nout,*) 'ili7pa not defined'
         stop 'STOP rininet'
      endif
      if (ibe7beta.eq.-1) then
         write (nout,*) 'ibe7beta not defined'
         stop 'STOP rininet'
      endif
      if (ibe7pg.eq.-1) then
         write (nout,*) 'ibe7pg not defined'
         stop 'STOP rininet'
      endif
      if (ib8beta.eq.-1) then
         write (nout,*) 'ib8beta not defined'
         stop 'STOP rininet'
      endif
      if (ib11pa.eq.-1) then
         write (nout,*) 'ib11pa not defined'
         stop 'STOP rininet'
      endif
      if (icpg.eq.-1) then
         write (nout,*) 'icpg not defined'
         stop 'STOP rininet'
      endif
      if (ic13pg.eq.-1) then
         write (nout,*) 'ic13pg not defined'
         stop 'STOP rininet'
      endif
      if (in13beta.eq.-1) then
         write (nout,*) 'in13beta not defined'
         stop 'STOP rininet'
      endif
      if (in14pg.eq.-1) then
         write (nout,*) 'in14pg not defined'
         stop 'STOP rininet'
      endif
      if (in15pa.eq.-1) then
         write (nout,*) 'in15pa not defined'
         stop 'STOP rininet'
      endif
      if (in15pg.eq.-1) then
         write (nout,*) 'in15pg not defined'
         stop 'STOP rininet'
      endif
      if (io15beta.eq.-1) then
         write (nout,*) 'io15beta not defined'
         stop 'STOP rininet'
      endif
      if (io16pg.eq.-1) then
         write (nout,*) 'ippd not defined'
         stop 'STOP rininet'
      endif
      if (io17pa.eq.-1) then
         write (nout,*) 'io17pa not defined'
         stop 'STOP rininet'
      endif
      if (ine20ag.eq.-1) then
         write (nout,*) 'ine20ag not defined'
         stop 'STOP rininet'
      endif
      if (img24ag.eq.-1) then
         write (nout,*) 'img24ag not defined'
         stop 'STOP rininet'
      endif
      if (isi28ag.eq.-1) then
         write (nout,*) 'isi28ag not defined'
         stop 'STOP rininet'
      endif

*________________________________
***   read nuclear reactions name
*--------------------------------

      rewind (99)
      do i = 1,33+nsp+1
         read (99,1200)
      enddo
      do j = 1,nreac
         read (99,1300) react(j)
      enddo

      close (99)  !  starevol.par

*______________________________________
***   initialize nuclear reaction types
*--------------------------------------

      do l = 1,nreac
***   index of neutron emission reactions
         if (k7(l).eq.1) then
            kneut = kneut+1
            jneut(kneut) = l
         endif
      enddo

***   LiBeB reactions
      jlibeb(1) = ili7pa
      jlibeb(2) = ib11pa
      jlibeb(3) = ili6pa

***   check dimension compatibility
 10   if (kneut.ne.nneut.or.kbeta.ne.nbeta.or.kbetaec.ne.nbetaec) then
         write (nout,*) 'Wrong dimension for vectors jneut/jbeta'
         if (kneut.gt.nneut) write (nout,*) 'The number of detected ',
     &        'reactions involving neutrons is kneut =',kneut,', it',
     &        ' must equal nneut =',nneut,' as defined in evolpar.star'
         if (kbeta.gt.nbeta) write (nout,*) 'The number of detected ',
     &        'beta decay reactions is kbeta =',kbeta,', it must equal',
     &        ' nbeta =',nbeta,' as defined in evolpar.star'
         if (kbetaec.gt.nbetaec) write (nout,*) 'The number of ',
     &        'detected electron capture reactions is kbetaec =',
     &        kbetaec,', it must equal nbetaec =',nbetaec,' as defined',
     &        ' in evolpar.star'
         write (nout,*) 'Check starevol.par (version >= 2.60 - 2.80)'
         stop 'STOP in rininet'
      endif


*_________________________________________________
***   initialize pointer for sparse nuclear matrix
*-------------------------------------------------

C$$$      v(1:nreac) = 1.d2
C$$$      y(1:nsp) = 1.d-2
C$$$      call freac (v,eqd,y,1.d-2,.true.)
C$$$      na = 0
C$$$      ia = 1
C$$$      ina(ia) = 1
C$$$      do  i = 1, nis
C$$$        do  j = 1, nis
C$$$           if (eqd(i,j) .ne. 0.d0) then
C$$$              na = na+1
C$$$              jna(na) = j
C$$$c              print *,na,j
C$$$           endif
C$$$        enddo
C$$$        ia = ia+1
C$$$        ina(ia) = na+1
C$$$      enddo
C$$$c      print *,ina

 1000 format (1x,i3,1x,a5,1x,f9.5,1x,f4.0,3x,1pe11.5,2x,a1)
 1100 format (i4,3x,i2,8x,i2,7x,i2,8x,i2,13x,f7.3,7x,4(i4),8x,f7.3)
 1200 format (1x)
 1300 format (7x,a37)

      return
      end
