

************************************************************************

      SUBROUTINE nucprod(v,y,tk,epp,eCNO,e3a,eC12,eNe20,eO16)

************************************************************************
* Compute energy generation rate for specific reactions                *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.nuc'

      double precision epp,eCNO,e3a,eC12,eNe20,eO16
      double precision ylim,tk
      double precision v(nreac),y(nsp)

      ylim = 1.d-8
      epp = 0.d0
      eCNO = 0.d0
      e3a = 0.d0
      eC12 = 0.d0
      eNe20 = 0.d0
      eO16 = 0.d0

c PP chains
c.. exclude beta decay reactions

c      if (y(ih1).ge.ylim.or.nmixd.gt.0) then
      if (y(ih1).ge.ylim) then
         epp = y(ih1)**2*v(ippg)*qi(ippg)*0.5d0+
     &   y(ih1)*y(ih2)*v(ipdg)*qi(ipdg)+
     &   y(ihe3)**2*v(i2he3)*qi(i2he3)*0.5d0+
     &   y(ihe4)*y(ihe3)*v(ihe3ag)*qi(ihe3ag)+
c     &   y(ibe7)*v(ibe7beta)*qi(ibe7beta)+
     &   y(ili7)*y(ih1)*v(ili7pa)*qi(ili7pa)+
     &   y(ibe7)*y(ih1)*v(ibe7pg)*qi(ibe7pg)
c     &   y(ib8)*v(ib8beta)*qi(ib8beta)
         epp = epp*econv

c CNO bicycle

         eCNO = y(ic12)*y(ih1)*v(icpg)*qi(icpg)+
c     &   y(in13)*v(in13beta)*qi(in13beta)+
     &   y(ic13)*y(ih1)*v(ic13pg)*qi(ic13pg)+
     &   y(in14)*y(ih1)*v(in14pg)*qi(in14pg)+
c     &   y(io15)*v(io15beta)*qi(io15beta)+
     &   y(in15)*y(ih1)*v(in15pa)*qi(in15pa)+
     &   y(in15)*y(ih1)*v(in15pg)*qi(in15pg)+
     &   y(io16)*y(ih1)*v(io16pg)*qi(io16pg)+
     &   y(io17)*y(ih1)*v(io17pa)*qi(io17pa)
         eCNO = eCNO*econv
      endif

c 3a and C12(a,g)
      if (y(ihe4).gt.ylim) then
         e3a = y(ihe4)**3*v(i3a)*qi(i3a)/6.d0+
     &   y(ihe4)*y(ic12)*v(icag)*qi(icag)
         e3a = e3a*econv
      endif

c Carbon burning 12C(12C,g)24Mg is not taken into account in the network
      if (y(ic12).gt.ylim.and.tk.gt.3.7d8) then
         eC12 = y(ic12)**2*v(iccga)*qi(iccga)*0.5d0+
     &        y(ic12)**2*v(iccgp)*qi(iccgp)*0.5d0
         eC12 = eC12*econv
      endif

c Neon photodisintegration/Burning
      if (y(ine20).gt.ylim.and.tk.gt.7.d8) then
         eNe20 = y(ine20)*v(inega)*qi(inega)
c         if (eNe20.gt.shlim.or.nphase.eq.7) then
c            eNe20 = eNe20+y(ine20)*y(ihe4)*v(ine20ag)*qi(ine20ag)+
c     &      y(img24)*y(ihe4)*v(img24ag)*qi(img24ag)+
c     &      y(isi28)*y(ihe4)*v(isi28ag)*qi(isi28ag)
c         endif
         eNe20 = eNe20*econv
      endif

c Oxygen Burning
      if (y(io16).gt.ylim.and.tk.gt.1.d9) then
         eO16 = 0.5d0*y(io16)**2*(v(iooga)*qi(iooga)+v(ioogp)*qi(ioogp))
         eO16 = eO16*econv
      endif

      return
      end
