
************************************************************************

      SUBROUTINE vit (tk,rok,mueinvk,ksh,v,iflag)

************************************************************************
* Calculate the nuclear reaction rates                                 *
* (rho*Navogadro*(sigma*vit)Maxwell, in 1/sec)                         *
*  iflag = 0 : defaut, compute all rates                               *
*  iflag = 1 : compute only rates associated to LiBeB (for HBB)        *
*  iflag = 2 : compute only electron capture rate (Ye dependence)      *
*                                                                      *
* $LastChangedDate:: 2014-02-04 15:45:03 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 11                                                          $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.mod'
      include 'evolcom.nuc'


      integer ivit,ite,iroe,iflag,idv0,iiv0
      integer i,j,jj,ksh

      double precision tburn,tnucmin,fnetdt,fnetdt0,tolnuc,tolnuc0
      double precision tk,rok,mueinvk,v
      double precision t9
      double precision tnuc9,vrate
      double precision fscr
      double precision dt9l,weight,t912,t913,t923,vsec,vmax,dh
      double precision a1,a2,b1,b2,vcapte,enuelec,rhoye,tcapte,romue

      common /funcscr/ fscr(nsh,nscr)
      common /nucvit/ tnuc9(nvit),vrate(nvit,nre)
      common /electcap/ rhoye(nroe),tcapte(nte),vcapte(nte,nroe,nrce),
     &     enuelec(nte,nroe,nrce)
      common /electcap_par/ ite,iroe,a1,a2,b1,b2
      common /nuc_param/ tburn,tnucmin,fnetdt,fnetdt0,tolnuc,tolnuc0

      dimension v(nreac),idv0(nreac)

      t9 = tk*1.d-9

*____________________________
***   Treatment of LiBeB only
*----------------------------

      if (iflag.eq.1) then
         do i = 1,nvit
            if (t9.lt.tnuc9(i)) goto 5
         enddo
 5       ivit = i-1
         if (ivit.eq.0) then
            v(ili7pa) = 1.d-99
            v(ib11pa) = 1.d-99
            v(ili6pa) = 1.d-99
            return
         else
            if (ivit.le.(nvit-1)) then
               dt9l = log10(t9/tnuc9(ivit))
               weight = dt9l/(log10(tnuc9(ivit+1)/tnuc9(ivit)))
               do i = 1,nlibeb
                  j = jlibeb(i)
                  v(j) = vrate(ivit,j)+weight*(vrate(ivit+1,j)-
     &                 vrate(ivit,j))
                  v(j) = 10.d0**v(j)
               enddo
            else
               do i = 1,nlibeb
                  j = jlibeb(i)
                  v(j) = 10.d0**vrate(nvit,j)
               enddo
            endif
         endif
c..   1 LI  6 ( 1 PROT , 1 HE  3)  1 HE  4
         v(ili6pa) = v(ili6pa)*rok*fscr(ksh,3)
c..   1 LI  7 ( 1 PROT , 0 OOOOO)  2 HE  4
         v(ili7pa) = v(ili7pa)*rok*fscr(ksh,3)
c..   1 B  11 ( 1 PROT , 0 OOOOO)  3 HE  4
         v(ib11pa) = v(ib11pa)*rok*fscr(ksh,5)
         return
      endif


*________________________________________
***   Treatment of electron captures only
*----------------------------------------
      if (iflag.eq.2) goto 25


*________________________________________________
***   interpolation of the nuclear reaction rates
*------------------------------------------------

      do i = 1,nvit
         if (t9.lt.tnuc9(i)) goto 10
      enddo
 10   ivit = i-1
      if (ivit.eq.0) then
         do j = 1,nre
            v(j) = 0.d0
         enddo
         do j = 1,nbeta-nbetaec
            jj = jbeta(j)
            v(jj) = vrate(1,jj)
            v(jj) = 10.d0**v(jj)
         enddo
***   1 BE  7 ( 0 betap, 0 nutri)  1 LI  7
         t912 = dsqrt(t9)
         t913 = t9**pw13
         t923 = t913*t913
         vsec = (1.34d-10/t912*(1.d0-0.537d0*t913+3.86d0*t923+2.7d-3/t9*
     &        exp(min(220.d0,2.515d-3/t9))))*rok*mueinvk
         if (t9.le.1.d-3) then
            vmax = 1.5032d-7
            vsec = min(vsec,vmax)
         endif
         v(ibe7beta) = vsec
         if (nphase.gt.2.and.tk.lt.tnucmin) goto 25
         do j = 1,idneut
            jj = jneut(j)
            v(jj) = vrate(1,jj)
            v(jj) = 10.d0**v(jj)
         enddo
         goto 20
      endif
      if (ivit.gt.(nvit-1)) then
         do j = 1,nre
            v(j) = vrate(nvit,j)
            v(j) = 10.d0**v(j)
         enddo
         goto 20
      endif
      dt9l = log10(t9/tnuc9(ivit))
      weight = dt9l/(log10(tnuc9(ivit+1)/tnuc9(ivit)))
      do j = 1,nre
         v(j) = vrate(ivit,j)+weight*(vrate(ivit+1,j)-vrate(ivit,j))
         v(j) = 10.d0**v(j)
      enddo

 20   continue

      iiv0 = 0
      do j = 1,nre
         if (v(j).lt.1.d-97) then
            iiv0 = iiv0+1
            idv0(iiv0) = j
         endif
      enddo

***   1 PROT  ( 1 NEUT , 0 OOOOO)  1 DEUT 
      v(1) = v(1)*rok
***   2 PROT  ( 0 OOOOO, 0 OOOOO)  1 DEUT 
      v(2) = v(2)*rok*fscr(ksh,1)
***   1 DEUT  ( 1 PROT , 0 OOOOO)  1 HE  3
      v(3) = v(3)*rok*fscr(ksh,1)
***   2 DEUT  ( 0 OOOOO, 0 OOOOO)  1 HE  4
      v(4) = v(4)*rok*fscr(ksh,1)
***   2 DEUT  ( 0 OOOOO, 1 NEUT )  1 HE  3
      v(5) = v(5)*rok*fscr(ksh,1)
***   1 HE  3 ( 1 NEUT , 0 OOOOO)  1 HE  4
      v(6) = v(6)*rok
***   1 HE  3 ( 1 DEUT , 1 PROT )  1 HE  4
      v(7) = v(7)*rok*fscr(ksh,2)
***   2 HE  3 ( 0 OOOOO, 2 PROT )  1 HE  4
      v(8) = v(8)*rok*fscr(ksh,19)
***   1 HE  4 ( 1 DEUT , 0 OOOOO)  1 LI  6
      v(9) = v(9)*rok*fscr(ksh,2)
***   1 HE  4 ( 1 HE  3, 0 OOOOO)  1 BE  7
      v(10) = v(10)*rok*fscr(ksh,19)
***   3 HE  4 ( 0 OOOOO, 0 OOOOO)  1 C  12
      v(11) = v(11)*rok*rok*fscr(ksh,19)*fscr(ksh,21)
***   1 LI  6 ( 1 PROT , 0 OOOOO)  1 BE  7
      v(12) = v(12)*rok*fscr(ksh,3)
***   1 LI  6 ( 1 PROT , 1 HE  3)  1 HE  4
      v(13) = v(13)*rok*fscr(ksh,3)
***   1 LI  7 ( 1 PROT , 0 OOOOO)  2 HE  4
      v(14) = v(14)*rok*fscr(ksh,3)
***   1 LI  7 ( 1 DEUT , 1 NEUT )  2 HE  4
      v(15) = v(15)*rok*fscr(ksh,3)
***   1 BE  7 ( 0 betap, 0 nutri)  1 LI  7
      if (ivit.ne.0) then
         vsec = v(ibe7beta)*rok*mueinvk
         if (t9.le.1.d-3) then
            vmax = 1.5032d-7
            vsec = min(vsec,vmax)
         endif
         v(ibe7beta) = vsec
      endif
***   1 BE  7 ( 1 PROT , 0 OOOOO)  1 B   8
      v(17) = v(17)*rok*fscr(ksh,4)
***   1 BE  7 ( 1 DEUT , 1 PROT )  2 HE  4
      v(18) = v(18)*rok*fscr(ksh,4)
***   1 B   8 ( 0 OOOOO, 0 OOOOO)  2 HE  4

***   1 B   8 ( 0 OOOOO, 1 PROT )  1 BE  7

***   1 BE  9 ( 1 PROT , 0 OOOOO)  1 B  10
      v(21) = v(21)*rok*fscr(ksh,4)
***   1 BE  9 ( 1 PROT , 1 DEUT )  2 HE  4
      v(22) = v(22)*rok*fscr(ksh,4)
***   1 BE  9 ( 1 PROT , 1 HE  4)  1 LI  6
      v(23) = v(23)*rok*fscr(ksh,4)
***   1 B  10 ( 1 PROT , 1 HE  4)  1 BE  7
      v(24) = v(24)*rok*fscr(ksh,5)
***   1 B  10 ( 1 PROT , 0 OOOOO)  1 B  11
      v(25) = v(25)*rok*fscr(ksh,5)
***   1 B  11 ( 1 PROT , 0 OOOOO)  1 C  12
      v(26) = v(26)*rok*fscr(ksh,5)
***   1 B  11 ( 1 PROT , 0 OOOOO)  3 HE  4
      v(27) = v(27)*rok*fscr(ksh,5)
***   1 B  11 ( 1 HE  4, 1 PROT )  1 C  14
      v(28) = v(28)*rok*fscr(ksh,22)
***   1 C  12 ( 1 NEUT , 0 OOOOO)  1 C  13
      v(29) = v(29)*rok
***   1 C  12 ( 1 PROT , 0 OOOOO)  1 N  13
      v(30) = v(30)*rok*fscr(ksh,6)
***   1 C  12 ( 1 HE  4, 0 OOOOO)  1 O  16
      v(31) = v(31)*rok*fscr(ksh,23)
***   2 C  12 ( 0 OOOOO, 1 HE  4)  1 NE 20
      v(32) = v(32)*rok*fscr(ksh,35)
***   2 C  12 ( 0 OOOOO, 1 PROT )  1 NA 23
      v(33) = v(33)*rok*fscr(ksh,35)
***   1 C  13 ( 1 NEUT , 0 OOOOO)  1 C  14
      v(34) = v(34)*rok
***   1 C  13 ( 1 PROT , 0 OOOOO)  1 N  14
      v(35) = v(35)*rok*fscr(ksh,6)
***   1 C  13 ( 1 HE  4, 1 NEUT )  1 O  16
      v(36) = v(36)*rok*fscr(ksh,23)
***   1 N  13 ( 0 OOOOO, 0 OOOOO)  1 C  13

***   1 N  13 ( 1 NEUT , 1 PROT )  1 C  13
      v(38) = v(38)*rok
***   1 N  13 ( 1 PROT , 0 OOOOO)  1 N  14
      v(39) = v(39)*rok*fscr(ksh,7)
***   1 C  14 ( 1 NEUT , 0 OOOOO)  1 N  15
      v(40) = v(40)*rok
***   1 C  14 ( 1 PROT , 0 OOOOO)  1 N  15
      v(41) = v(41)*rok*fscr(ksh,6)
***   1 C  14 ( 1 PROT , 1 NEUT )  1 N  14
      v(42) = v(42)*rok*fscr(ksh,6)
***   1 C  14 ( 1 PROT , 1 HE  4)  1 B  11
      v(43) = v(43)*rok*fscr(ksh,6)
***   1 C  14 ( 1 HE  4, 1 NEUT )  1 O  17
      v(44) = v(44)*rok*fscr(ksh,23)
***   1 C  14 ( 1 HE  4, 0 OOOOO)  1 O  18
      v(45) = v(45)*rok*fscr(ksh,23)
***   1 N  14 ( 1 NEUT , 0 OOOOO)  1 N  15
      v(46) = v(46)*rok
***   1 N  14 ( 1 NEUT , 1 PROT )  1 C  14
      v(47) = v(47)*rok
***   1 N  14 ( 1 PROT , 0 OOOOO)  1 O  15
      v(48) = v(48)*rok*fscr(ksh,7)
***   1 N  14 ( 1 HE  4, 0 OOOOO)  1 F  18
      v(49) = v(49)*rok*fscr(ksh,24)
***   1 N  15 ( 1 NEUT , 0 OOOOO)  1 O  16
      v(50) = v(50)*rok
***   1 N  15 ( 1 PROT , 1 HE  4)  1 C  12
      v(51) = v(51)*rok*fscr(ksh,7)
***   1 N  15 ( 1 PROT , 0 OOOOO)  1 O  16
      v(52) = v(52)*rok*fscr(ksh,7)
***   1 N  15 ( 1 HE  4, 0 OOOOO)  1 F  19
      v(53) = v(53)*rok*fscr(ksh,24)
***   1 O  15 ( 0 OOOOO, 0 OOOOO)  1 N  15

***   1 O  15 ( 1 NEUT , 1 HE  4)  1 C  12
      v(55) = v(55)*rok
***   1 O  15 ( 1 NEUT , 1 PROT )  1 N  15
      v(56) = v(56)*rok
***   1 O  16 ( 1 NEUT , 0 OOOOO)  1 O  17
      v(57) = v(57)*rok
***   1 O  16 ( 1 PROT , 0 OOOOO)  1 O  17
      v(58) = v(58)*rok*fscr(ksh,8)
***   1 O  16 ( 1 HE  4, 0 OOOOO)  1 NE 20
      v(59) = v(59)*rok*fscr(ksh,25)
***   1 O  16 ( 1 C  12, 1 NEUT )  1 AL 27
      v(60) = v(60)*rok*fscr(ksh,37)
***   1 O  16 ( 1 C  12, 1 PROT )  1 AL 27
      v(61) = v(61)*rok*fscr(ksh,37)
***   1 O  16 ( 1 C  12, 1 HE  4)  1 MG 24
      v(62) = v(62)*rok*fscr(ksh,37)
***   2 O  16 ( 0 OOOOO, 1 PROT )  1 P  31
      v(63) = v(63)*rok*fscr(ksh,36)
***   2 O  16 ( 0 OOOOO, 1 HE  4)  1 SI 28
      v(64) = v(64)*rok*fscr(ksh,36)
***   1 O  17 ( 1 NEUT , 1 HE  4)  1 C  14
      v(65) = v(65)*rok
***   1 O  17 ( 1 NEUT , 0 OOOOO)  1 O  18
      v(66) = v(66)*rok
***   1 O  17 ( 1 PROT , 0 OOOOO)  1 F  18
      v(67) = v(67)*rok*fscr(ksh,8)
***   1 O  17 ( 1 PROT , 1 HE  4)  1 N  14
      v(68) = v(68)*rok*fscr(ksh,8)
***   1 O  17 ( 1 HE  4, 1 NEUT )  1 NE 20
      v(69) = v(69)*rok*fscr(ksh,25)
***   1 O  17 ( 1 HE  4, 0 OOOOO)  1 NE 21
      v(70) = v(70)*rok*fscr(ksh,25)
***   1 O  18 ( 1 NEUT , 0 OOOOO)  1 F  19
      v(71) = v(71)*rok
***   1 O  18 ( 1 PROT , 0 OOOOO)  1 F  19
      v(72) = v(72)*rok*fscr(ksh,8)
***   1 O  18 ( 1 PROT , 1 HE  4)  1 N  15
      v(73) = v(73)*rok*fscr(ksh,8)
***   1 O  18 ( 1 HE  4, 0 OOOOO)  1 NE 22
      v(74) = v(74)*rok*fscr(ksh,25)
***   1 O  18 ( 1 HE  4, 1 NEUT )  1 NE 21
      v(75) = v(75)*rok*fscr(ksh,25)
***   1 F  18 ( 0 OOOOO, 0 OOOOO)  1 O  18

***   1 F  18 ( 1 NEUT , 1 PROT )  1 O  18
      v(77) = v(77)*rok
***   1 F  18 ( 1 NEUT , 1 HE  4)  1 N  15
      v(78) = v(78)*rok
***   1 F  18 ( 1 HE  4, 1 PROT )  1 NE 21
      v(79) = v(79)*rok*fscr(ksh,26)
***   1 F  19 ( 1 NEUT , 0 OOOOO)  1 NE 20
      v(80) = v(80)*rok
***   1 F  19 ( 1 PROT , 0 OOOOO)  1 NE 20
      v(81) = v(81)*rok*fscr(ksh,9)
***   1 F  19 ( 1 PROT , 1 HE  4)  1 O  16
      v(82) = v(82)*rok*fscr(ksh,9)
***   1 F  19 ( 1 HE  4, 1 PROT )  1 NE 22
      v(83) = v(83)*rok*fscr(ksh,26)
***   1 NE 20 ( 1 NEUT , 0 OOOOO)  1 NE 21
      v(84) = v(84)*rok
***   1 NE 20 ( 1 PROT , 0 OOOOO)  1 NE 21
      v(85) = v(85)*rok*fscr(ksh,10)
***   1 NE 20 ( 0 OOOOO, 1 HE  4)  1 O  16

***   1 NE 20 ( 1 HE  4, 0 OOOOO)  1 MG 24
      v(87) = v(87)*rok*fscr(ksh,27)
***   1 NE 20 ( 1 HE  4, 1 PROT )  1 NA 23
      v(88) = v(88)*rok*fscr(ksh,27)
***   1 NE 21 ( 1 NEUT , 0 OOOOO)  1 NE 22
      v(89) = v(89)*rok
***   1 NE 21 ( 1 NEUT , 1 HE  4)  1 O  18
      v(90) = v(90)*rok
***   1 NE 21 ( 1 PROT , 0 OOOOO)  1 NA 22
      v(91) = v(91)*rok*fscr(ksh,10)
***   1 NE 21 ( 1 HE  4, 0 OOOOO)  1 MG 25
      v(92) = v(92)*rok*fscr(ksh,27)
***   1 NE 21 ( 1 HE  4, 1 NEUT )  1 MG 24
      v(93) = v(93)*rok*fscr(ksh,27)
***   1 NE 22 ( 1 NEUT , 0 OOOOO)  1 NA 23
      v(94) = v(94)*rok
***   1 NE 22 ( 1 PROT , 0 OOOOO)  1 NA 23
      v(95) = v(95)*rok*fscr(ksh,10)
***   1 NE 22 ( 1 HE  4, 1 NEUT )  1 MG 25
      v(96) = v(96)*rok*fscr(ksh,27)
***   1 NE 22 ( 1 HE  4, 0 OOOOO)  1 MG 26
      v(97) = v(97)*rok*fscr(ksh,27)
***   1 NA 22 ( 0 OOOOO, 0 OOOOO)  1 NE 22

***   1 NA 22 ( 1 NEUT , 0 OOOOO)  1 NA 23
      v(99) = v(99)*rok
***   1 NA 22 ( 1 NEUT , 1 PROT )  1 NE 22
      v(100) = v(100)*rok
***   1 NA 22 ( 1 PROT , 0 OOOOO)  1 NA 23
      v(101) = v(101)*rok*fscr(ksh,11)
***   1 NA 23 ( 1 NEUT , 0 OOOOO)  1 MG 24
      v(102) = v(102)*rok
***   1 NA 23 ( 1 PROT , 1 HE  4)  1 NE 20
      v(103) = v(103)*rok*fscr(ksh,11)
***   1 NA 23 ( 1 PROT , 0 OOOOO)  1 MG 24
      v(104) = v(104)*rok*fscr(ksh,11)
***   1 NA 23 ( 1 HE  4, 1 PROT )  1 MG 26
      v(105) = v(105)*rok*fscr(ksh,28)
***   1 MG 24 ( 1 NEUT , 0 OOOOO)  1 MG 25
      v(106) = v(106)*rok
***   1 MG 24 ( 1 PROT , 0 OOOOO)  1 MG 25
      v(107) = v(107)*rok*fscr(ksh,12)
***   1 MG 24 ( 1 HE  4, 1 PROT )  1 AL 27
      v(108) = v(108)*rok*fscr(ksh,29)
***   1 MG 24 ( 1 HE  4, 0 OOOOO)  1 SI 28
      v(109) = v(109)*rok*fscr(ksh,29)
***   1 MG 25 ( 1 NEUT , 0 OOOOO)  1 MG 26
      v(110) = v(110)*rok
***   1 MG 25 ( 1 PROT , 0 OOOOO)  1 AL26m
      v(111) = v(111)*rok*fscr(ksh,12)
***   1 MG 25 ( 1 PROT , 0 OOOOO)  1 AL26g
      v(112) = v(112)*rok*fscr(ksh,12)
***   1 MG 25 ( 1 HE  4, 1 PROT )  1 SI 28
      v(113) = v(113)*rok*fscr(ksh,29)
***   1 MG 25 ( 1 HE  4, 1 NEUT )  1 SI 28
      v(114) = v(114)*rok*fscr(ksh,29)
***   1 MG 25 ( 1 HE  4, 0 OOOOO)  1 SI 29
      v(115) = v(115)*rok*fscr(ksh,29)
***   1 MG 26 ( 1 NEUT , 0 OOOOO)  1 AL 27
      v(116) = v(116)*rok
***   1 MG 26 ( 1 PROT , 0 OOOOO)  1 AL 27
      v(117) = v(117)*rok*fscr(ksh,12)
***   1 MG 26 ( 1 HE  4, 0 OOOOO)  1 SI 30
      v(118) = v(118)*rok*fscr(ksh,29)
***   1 MG 26 ( 1 HE  4, 1 PROT )  1 SI 29
      v(119) = v(119)*rok*fscr(ksh,29)
***   1 MG 26 ( 1 HE  4, 1 NEUT )  1 SI 29
      v(120) = v(120)*rok*fscr(ksh,29)
***   1 AL26m ( 0 OOOOO, 0 OOOOO)  1 MG 26

***   1 AL26m ( 1 NEUT , 1 PROT )  1 MG 26
      v(122) = v(122)*rok
***   1 AL26m ( 1 NEUT , 1 HE  4)  1 NA 23
      v(123) = v(123)*rok
***   1 AL26m ( 1 NEUT , 0 OOOOO)  1 AL 27
      v(124) = v(124)*rok
***   1 AL26m ( 1 PROT , 0 OOOOO)  1 AL 27
      v(125) = v(125)*rok*fscr(ksh,13)
***   1 AL26g ( 0 OOOOO, 0 OOOOO)  1 MG 26

***   1 AL26g ( 1 NEUT , 0 OOOOO)  1 AL 27
      v(127) = v(127)*rok
***   1 AL26g ( 1 NEUT , 1 PROT )  1 MG 26
      v(128) = v(128)*rok
***   1 AL26g ( 1 NEUT , 1 HE  4)  1 NA 23
      v(129) = v(129)*rok
***   1 AL26g ( 1 PROT , 0 OOOOO)  1 AL 27
      v(130) = v(130)*rok*fscr(ksh,13)
***   1 AL26g ( 0 OOOOO, 0 OOOOO)  1 AL26m

***   1 AL26m ( 0 OOOOO, 0 OOOOO)  1 AL26g

***   1 AL 27 ( 1 PROT , 0 OOOOO)  1 SI 28
      v(133) = v(133)*rok*fscr(ksh,13)
***   1 AL 27 ( 1 PROT , 1 HE  4)  1 MG 24
      v(134) = v(134)*rok*fscr(ksh,13)
***   1 AL 27 ( 1 HE  4, 1 NEUT )  1 SI 30
      v(135) = v(135)*rok*fscr(ksh,30)
***   1 AL 27 ( 1 HE  4, 1 PROT )  1 SI 30
      v(136) = v(136)*rok*fscr(ksh,30)
***   1 SI 28 ( 1 NEUT , 0 OOOOO)  1 SI 29
      v(137) = v(137)*rok
***   1 SI 28 ( 1 PROT , 0 OOOOO)  1 SI 29
      v(138) = v(138)*rok*fscr(ksh,14)
***   1 SI 28 ( 1 HE  4, 0 OOOOO)  1 S  32
      v(139) = v(139)*rok*fscr(ksh,31)
***   1 SI 29 ( 1 NEUT , 0 OOOOO)  1 SI 30
      v(140) = v(140)*rok
***   1 SI 29 ( 1 PROT , 0 OOOOO)  1 SI 30
      v(141) = v(141)*rok*fscr(ksh,14)
***   1 SI 30 ( 1 NEUT , 0 OOOOO)  1 P  31
      v(142) = v(142)*rok
***   1 SI 30 ( 1 PROT , 0 OOOOO)  1 P  31
      v(143) = v(143)*rok*fscr(ksh,14)
***   1 P  31 ( 1 NEUT , 0 OOOOO)  1 S  32
      v(144) = v(144)*rok
***   1 P  31 ( 1 PROT , 1 HE  4)  1 SI 28
      v(145) = v(145)*rok*fscr(ksh,15)
***   1 P  31 ( 1 PROT , 0 OOOOO)  1 S  32
      v(146) = v(146)*rok*fscr(ksh,15)
***   1 S  32 ( 1 NEUT , 0 OOOOO)  1 S  33
      v(147) = v(147)*rok
***   1 S  32 ( 1 NEUT , 1 HE  4)  1 SI 29
      v(148) = v(148)*rok
***   1 S  32 ( 1 PROT , 0 OOOOO)  1 S  33
      v(149) = v(149)*rok*fscr(ksh,16)
***   1 S  33 ( 1 NEUT , 0 OOOOO)  1 S  34
      v(150) = v(150)*rok
***   1 S  33 ( 1 NEUT , 1 HE  4)  1 SI 30
      v(151) = v(151)*rok
***   1 S  33 ( 1 PROT , 0 OOOOO)  1 S  34
      v(152) = v(152)*rok*fscr(ksh,16)
***   1 S  34 ( 1 NEUT , 0 OOOOO)  1 S  35
      v(153) = v(153)*rok
***   1 S  34 ( 1 PROT , 0 OOOOO)  1 CL 35
      v(154) = v(154)*rok*fscr(ksh,16)
***   1 S  35 ( 0 OOOOO, 0 OOOOO)  1 CL 35

***   1 S  35 ( 1 NEUT , 0 OOOOO)  1 S  36
      v(156) = v(156)*rok
***   1 S  35 ( 1 PROT , 0 OOOOO)  1 CL 36
      v(157) = v(157)*rok*fscr(ksh,16)
***   1 CL 35 ( 1 NEUT , 0 OOOOO)  1 CL 36
      v(158) = v(158)*rok
***   1 CL 35 ( 1 PROT , 0 OOOOO)  1 HEAVY
      v(159) = v(159)*rok*fscr(ksh,17)
***   1 S  36 ( 1 NEUT , 0 OOOOO)  1 CL 37
      v(160) = v(160)*rok
***   1 S  36 ( 1 PROT , 0 OOOOO)  1 CL 37
      v(161) = v(161)*rok*fscr(ksh,16)
***   1 CL 36 ( 0 OOOOO, 0 OOOOO)  1 HEAVY

***   1 CL 36 ( 1 NEUT , 0 OOOOO)  1 CL 37
      v(163) = v(163)*rok
***   1 CL 36 ( 1 PROT , 0 OOOOO)  1 CL 37
      v(164) = v(164)*rok*fscr(ksh,17)
***   1 CL 37 ( 1 NEUT , 0 OOOOO)  1 HEAVY
      v(165) = v(165)*rok
***   1 CL 37 ( 1 PROT , 0 OOOOO)  1 HEAVY
      v(166) = v(166)*rok*fscr(ksh,17)
***   1 HEAVY ( 1 NEUT , 0 OOOOO)  1 captn
      v(167) = (5.11404d0+signt(ksh))*1.44d5*rok

      do j = 1,iiv0
         v(idv0(j)) = 1.d-99
      enddo


*______________________________________________
***   electron capture rates (rho dependence)
* 166    1 N  14 ( 0 capte, 0 nutri)  1 C  14
* 167    1 C  14 ( 0 betam, 0 nutri)  1 N  14
* 168    1 NE 20 ( 0 capte, 0 nutri)  1 F  20
* 169    1 F  20 ( 0 betam, 0 nutri)  1 NE 20
* 170    1 NA 23 ( 0 capte, 0 nutri)  1 NE 23
* 171    1 NE 23 ( 0 betam, 0 nutri)  1 NA 23
* 172    1 MG 24 ( 0 capte, 0 nutri)  1 NA 24
* 173    1 NA 24 ( 0 betam, 0 nutri)  1 MG 24
* 174    1 MG 25 ( 0 capte, 0 nutri)  1 NA 25
* 175    1 NA 25 ( 0 betam, 0 nutri)  1 MG 25
* 176    1 AL 27 ( 0 capte, 0 nutri)  1 MG 27
* 177    1 MG 27 ( 0 betam, 0 nutri)  1 AL 27
*----------------------------------------------

 25   romue = log10(rok*mueinvk)
      jj = 0

***   extended table (high T and rho)
***   determine index "ite" for temperature
      if (t9.le.tcapte(1)) then
         ite = 1
         jj = 1
         a1 = 0.d0
         b1 = 1.d0
      endif
      if (t9.ge.tcapte(nte)) then
         ite = nte-1
         jj = 2
         a1 = 1.d0
         b1 = 0.d0
      endif
      if (jj.eq.0) then
         do i = 2,nte
            if (t9.lt.tcapte(i)) goto 50
         enddo
 50      ite = i-1
         dh = tcapte(ite+1)-tcapte(ite)
         a1 = (t9-tcapte(ite))/dh
         b1 = 1.d0-a1
      endif
***   determine index "iroe" for rho*Ye
      if (romue.le.rhoye(1)) then
         iroe = 1
         jj = 3
         a2 = 0.d0
         b2 = 1.d0
      endif
      if (romue.ge.rhoye(nroe)) then
         iroe = nroe-1
         jj = 4
         a2 = 1.d0
         b2 = 0.d0
      endif
      if (jj.lt.3) then
         do i = 2,nroe
            if (romue.lt.rhoye(i)) exit
         enddo
         iroe = i-1
         dh = rhoye(iroe+1)-rhoye(iroe)
         a2 = (romue-rhoye(iroe))/dh
         b2 = 1.d0-a2
      endif

***   determine reaction rates (bilinear interpolation)
      do j = nre+1,nreac
         i = j-nre
         v(j) = a1*a2*vcapte(ite+1,iroe+1,i)+b1*b2*
     &        vcapte(ite,iroe,i)+a1*b2*vcapte(ite+1,iroe,i)+
     &        b1*a2*vcapte(ite,iroe+1,i)
         v(j) = 10.d0**v(j)
      enddo

      return
      end
