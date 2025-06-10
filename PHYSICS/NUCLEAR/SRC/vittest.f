
      PROGRAM VITTEST

************************************************************************
* Test of the tabulated nuclear reaction rates at various temperatures *
************************************************************************

      implicit double precision (a-h,o-z)
      parameter (nre = 180)
      parameter (ntest = 31)
      dimension v(nre)

      open (unit = 10,file = 'vit.out')

      rok = 1.d0

      write (10,1000)

      tlogd = 6.d0
      tlogu = 9.d0
      dtlog = (tlogu-tlogd)/(ntest-1)
      tlog = tlogd-dtlog
      do i = 1,ntest
        tlog = tlog+dtlog
        t9 = 10.d0**tlog*1.d-9
        call vit (rok,t9,v)

        write (10,1100) tlog
        write (10,1200) (m,v(m),m = 1,10)
        write (10,1200) (m,v(m),m = 11,20)
        write (10,1200) (m,v(m),m = 21,30)
        write (10,1200) (m,v(m),m = 31,40)
        write (10,1200) (m,v(m),m = 41,50)
        write (10,1200) (m,v(m),m = 51,60)
        write (10,1200) (m,v(m),m = 61,70)
        write (10,1200) (m,v(m),m = 71,80)
        write (10,1200) (m,v(m),m = 81,90)
        write (10,1200) (m,v(m),m = 91,100)
        write (10,1200) (m,v(m),m = 101,110)
        write (10,1200) (m,v(m),m = 111,120)
        write (10,1200) (m,v(m),m = 121,130)
        write (10,1200) (m,v(m),m = 131,140)
        write (10,1200) (m,v(m),m = 141,150)
        write (10,1200) (m,v(m),m = 151,160)
        write (10,1200) (m,v(m),m = 161,170)
        write (10,1200) (m,v(m),m = 171,180)

      enddo

 1000 format (/,1x,'Nuclear reaction rates from subr. VIT, at ',
     &       'various temperatures',/)
 1100 format (/,1x,'*** at T = ',1pe10.4,' K:',/)
 1200 format (10(1x,'#',i3,'=',1pe9.3))

      close (10)

      end



************************************************************************

      SUBROUTINE VIT (rho,t9,v)

************************************************************************
* Calculate the nuclear reaction rates                                 *
* (rho*Navogadro*(sigma*vit)Maxwell, in 1/sec)                         *
************************************************************************

      implicit double precision (a-h,o-z)
      double precision mueinvk
      parameter (nre = 180)
      parameter (nvit = 250)
      common /nucvit/ vrate(nvit,nre)
      dimension v(nre),fscr(nre)
 
      signtk = 0.d0
      mueinvk = 1.d0

      reft9l = -3.3142341d0
      delt9l = 1.3204119d-2
      do i = 1,nvit
         reft9l = reft9l+delt9l
         reft9 = 10.d0**reft9l
         if (t9.lt.reft9) goto 10
      enddo
 10   ivit = i-1
      if (ivit.eq.0) then
         do j = 1,nre
            v(j) = vrate(1,j)
            v(j) = 10.d0**v(j)
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
      dt9l = log10(t9)-reft9l+delt9l
      weight = dt9l/delt9l
      do j = 1,nre
         v(j) = vrate(ivit,j)+weight*(vrate(ivit+1,j)-vrate(ivit,j))
         v(j) = 10.d0**v(j)
      enddo

*_________________________________________
***   calculation of the screening factors
*-----------------------------------------

 20   nscr = 0
c     fft1 = screen (1.d0,1.d0)
c     fft2 = screen (6.d0,1.d0)
c     fft3 = screen (2.d0,2.d0)
c     fft4 = screen (6.d0,2.d0)
c     if (fft1.gt.1.10d0.or.fft2.gt.1.10d0) nscr = 1
c     if (fft3.gt.1.10d0.or.fft4.gt.1.10d0) nscr = 2
c     if (nscr.ge.1) then
c        do k = 1,17
c           kk = 17+k
c           dk = dble(k)
c           fscr(k) = screen (dk,1.d0)
c           fscr(kk) = 1.d0
c           if (nscr.eq.2) fscr(kk) = screen (dk,2.d0)
c        enddo
c     else
         do k = 1,34
            fscr(k) = 1.d0
         enddo
c     endif
      fscr(35) = 1.d0
c     if (nscr.eq.2) fscr(35) = screen (6.d0,6.d0)

***   2 proto ( 0 gamma, 0 nutri)  1 deutr
      v(1) = v(1)*rho*fscr(1)
***   1 deutr ( 1 proto, 0 gamma)  1 He  3
      v(2) = v(2)*rho*fscr(1)
***   2 deutr ( 0 gamma, 0 gamma)  1 alpha
      v(3) = v(3)*rho*fscr(1)
***   2 deutr ( 0 gamma, 1 neutr)  1 He  3
      v(4) = v(4)*rho*fscr(1)
***   1 He  3 ( 1 deutr, 1 proto)  1 alpha
      v(5) = v(5)*rho*fscr(2)
***   2 He  3 ( 0 gamma, 2 proto)  1 alpha
      v(6) = v(6)*rho*fscr(19)
***   1 alpha ( 1 deutr, 0 gamma)  1 Li  6
      v(7) = v(7)*rho*fscr(2)
***   1 alpha ( 1 He  3, 0 gamma)  1 Be  7
      v(8) = v(8)*rho*fscr(19)
***   3 alpha ( 0 gamma, 0 gamma)  1 C  12
      v(9) = v(9)*rho*rho*fscr(19)*fscr(21)
***   1 Li  6 ( 1 proto, 0 gamma)  1 Be  7
      v(10) = v(10)*rho*fscr(3)
***   1 Li  6 ( 1 proto, 1 He  3)  1 alpha
      v(11) = v(11)*rho*fscr(3)
***   1 Li  7 ( 1 proto, 1 alpha)  1 alpha
      v(12) = v(12)*rho*fscr(3)
***   1 Li  7 ( 1 deutr, 1 neutr)  2 alpha
      v(13) = v(13)*rho*fscr(3)
***   1 Be  7 ( 0 betap, 0 nutri)  1 Li  7
      vsec = v(14)*rho*mueinvk
      if (t9.le.1.d-3) then
         vmax = 1.5032d-7
         vsec = min(vsec,vmax)
      endif
      v(14) = vsec
***   1 Be  7 ( 1 proto, 0 gamma)  1 B   8
      v(15) = v(15)*rho*fscr(4)
***   1 Be  7 ( 1 deutr, 1 proto)  2 alpha
      v(16) = v(16)*rho*fscr(4)
***   1 B   8 ( 0 betap, 0 nutri)  2 alpha

***   1 Be  9 ( 1 proto, 0 gamma)  1 B  10
      v(18) = v(18)*rho*fscr(4)
***   1 Be  9 ( 1 proto, 1 deutr)  2 alpha
      v(19) = v(19)*rho*fscr(4)
***   1 Be  9 ( 1 proto, 1 alpha)  1 Li  6
      v(20) = v(20)*rho*fscr(4)
***   1 B  10 ( 1 proto, 0 gamma)  1 C  11
      v(21) = v(21)*rho*fscr(5)
***   1 B  10 ( 1 proto, 1 alpha)  1 Be  7
      v(22) = v(22)*rho*fscr(5)
***   1 B  11 ( 1 proto, 0 gamma)  1 C  12
      v(23) = v(23)*rho*fscr(5)
***   1 B  11 ( 1 proto, 1 alpha)  2 alpha
      v(24) = v(24)*rho*fscr(5)
***   1 B  11 ( 1 alpha, 1 proto)  1 C  14
      v(25) = v(25)*rho*fscr(22)
***   1 C  11 ( 0 betap, 0 nutri)  1 B  11

***   1 C  11 ( 1 proto, 0 gamma)  1 C  12
      v(27) = v(27)*rho*fscr(6)
***   1 C  12 ( 1 neutr, 0 gamma)  1 C  13
      v(28) = v(28)*rho
***   1 C  12 ( 1 proto, 0 gamma)  1 N  13
      v(29) = v(29)*rho*fscr(6)
***   1 C  12 ( 1 alpha, 0 gamma)  1 O  16
      v(30) = v(30)*rho*fscr(23)
***   2 C  12 ( 0 gamma, 1 alpha)  1 Ne 20
      v(31) = v(31)*rho*fscr(35)
***   2 C  12 ( 0 gamma, 1 proto)  1 Na 23
      v(32) = v(32)*rho*fscr(35)
***   1 C  13 ( 1 neutr, 0 gamma)  1 C  14
      v(33) = v(33)*rho
***   1 C  13 ( 1 proto, 0 gamma)  1 N  14
      v(34) = v(34)*rho*fscr(6)
***   1 C  13 ( 1 alpha, 1 neutr)  1 O  16
      v(35) = v(35)*rho*fscr(23)
***   1 N  13 ( 0 betap, 0 nutri)  1 C  13

***   1 N  13 ( 1 neutr, 1 proto)  1 C  13
      v(37) = v(37)*rho
***   1 N  13 ( 1 proto, 0 gamma)  1 N  14
      v(38) = v(38)*rho*fscr(7)
***   1 C  14 ( 0 betam, 0 nutri)  1 N  14

***   1 C  14 ( 1 neutr, 0 gamma)  1 N  15
      v(40) = v(40)*rho
***   1 C  14 ( 1 proto, 0 gamma)  1 N  15
      v(41) = v(41)*rho*fscr(6)
***   1 C  14 ( 1 proto, 1 neutr)  1 N  14
      v(42) = v(42)*rho*fscr(6)
***   1 C  14 ( 1 proto, 1 alpha)  1 B  11
      v(43) = v(43)*rho*fscr(6)
***   1 C  14 ( 1 alpha, 1 neutr)  1 O  17
      v(44) = v(44)*rho*fscr(23)
***   1 C  14 ( 1 alpha, 0 gamma)  1 O  18
      v(45) = v(45)*rho*fscr(23)
***   1 N  14 ( 1 neutr, 0 gamma)  1 N  15
      v(46) = v(46)*rho
***   1 N  14 ( 1 neutr, 1 proto)  1 C  14
      v(47) = v(47)*rho
***   1 N  14 ( 1 proto, 0 gamma)  1 O  15
      v(48) = v(48)*rho*fscr(7)
***   1 N  14 ( 1 alpha, 0 gamma)  1 F  18
      v(49) = v(49)*rho*fscr(24)
***   1 N  15 ( 1 neutr, 0 gamma)  1 O  16
      v(50) = v(50)*rho
***   1 N  15 ( 1 proto, 1 alpha)  1 C  12
      v(51) = v(51)*rho*fscr(7)
***   1 N  15 ( 1 proto, 0 gamma)  1 O  16
      v(52) = v(52)*rho*fscr(7)
***   1 N  15 ( 1 alpha, 0 gamma)  1 F  19
      v(53) = v(53)*rho*fscr(24)
***   1 O  15 ( 0 betap, 0 nutri)  1 N  15

***   1 O  15 ( 1 neutr, 1 alpha)  1 C  12
      v(55) = v(55)*rho
***   1 O  15 ( 1 neutr, 1 proto)  1 N  15
      v(56) = v(56)*rho
***   1 O  16 ( 1 neutr, 0 gamma)  1 O  17
      v(57) = v(57)*rho
***   1 O  16 ( 1 proto, 0 gamma)  1 O  17
      v(58) = v(58)*rho*fscr(8)
***   1 O  16 ( 1 alpha, 0 gamma)  1 Ne 20
      v(59) = v(59)*rho*fscr(25)
***   1 O  17 ( 1 neutr, 1 alpha)  1 C  14
      v(60) = v(60)*rho
***   1 O  17 ( 1 neutr, 0 gamma)  1 O  18
      v(61) = v(61)*rho
***   1 O  17 ( 1 proto, 0 gamma)  1 F  18
      v(62) = v(62)*rho*fscr(8)
***   1 O  17 ( 1 proto, 1 alpha)  1 N  14
      v(63) = v(63)*rho*fscr(8)
***   1 O  17 ( 1 alpha, 1 neutr)  1 Ne 20
      v(64) = v(64)*rho*fscr(25)
***   1 O  17 ( 1 alpha, 0 gamma)  1 Ne 21
      v(65) = v(65)*rho*fscr(25)
***   1 O  18 ( 1 neutr, 0 gamma)  1 F  19
      v(66) = v(66)*rho
***   1 O  18 ( 1 proto, 0 gamma)  1 F  19
      v(67) = v(67)*rho*fscr(8)
***   1 O  18 ( 1 proto, 1 alpha)  1 N  15
      v(68) = v(68)*rho*fscr(8)
***   1 O  18 ( 1 alpha, 0 gamma)  1 Ne 22
      v(69) = v(69)*rho*fscr(25)
***   1 O  18 ( 1 alpha, 1 neutr)  1 Ne 21
      v(70) = v(70)*rho*fscr(25)
***   1 F  18 ( 0 betap, 0 nutri)  1 O  18

***   1 F  18 ( 1 neutr, 1 proto)  1 O  18
      v(72) = v(72)*rho
***   1 F  18 ( 1 neutr, 1 alpha)  1 N  15
      v(73) = v(73)*rho
***   1 F  18 ( 1 alpha, 1 proto)  1 Ne 21
      v(74) = v(74)*rho*fscr(26)
***   1 F  19 ( 1 neutr, 0 gamma)  1 Ne 20
      v(75) = v(75)*rho
***   1 F  19 ( 1 proto, 0 gamma)  1 Ne 20
      v(76) = v(76)*rho*fscr(9)
***   1 F  19 ( 1 proto, 1 alpha)  1 O  16
      v(77) = v(77)*rho*fscr(9)
***   1 F  19 ( 1 alpha, 1 proto)  1 Ne 22
      v(78) = v(78)*rho*fscr(26)
***   1 Ne 20 ( 1 neutr, 0 gamma)  1 Ne 21
      v(79) = v(79)*rho
***   1 Ne 20 ( 1 proto, 0 gamma)  1 Ne 21
      v(80) = v(80)*rho*fscr(10)
***   1 Ne 20 ( 1 alpha, 0 gamma)  1 Mg 24
      v(81) = v(81)*rho*fscr(27)
***   1 Ne 20 ( 1 alpha, 1 proto)  1 Na 23
      v(82) = v(82)*rho*fscr(27)
***   1 Ne 21 ( 1 neutr, 0 gamma)  1 Ne 22
      v(83) = v(83)*rho
***   1 Ne 21 ( 1 neutr, 1 alpha)  1 O  18
      v(84) = v(84)*rho
***   1 Ne 21 ( 1 proto, 0 gamma)  1 Na 22
      v(85) = v(85)*rho*fscr(10)
***   1 Ne 21 ( 1 alpha, 0 gamma)  1 Mg 25
      v(86) = v(86)*rho*fscr(27)
***   1 Ne 21 ( 1 alpha, 1 neutr)  1 Mg 24
      v(87) = v(87)*rho*fscr(27)
***   1 Ne 22 ( 1 neutr, 0 gamma)  1 Na 23
      v(88) = v(88)*rho
***   1 Ne 22 ( 1 proto, 0 gamma)  1 Na 23
      v(89) = v(89)*rho*fscr(10)
***   1 Ne 22 ( 1 alpha, 1 neutr)  1 Mg 25
      v(90) = v(90)*rho*fscr(27)
***   1 Ne 22 ( 1 alpha, 0 gamma)  1 Mg 26
      v(91) = v(91)*rho*fscr(27)
***   1 Na 22 ( 0 betap, 0 nutri)  1 Ne 22

***   1 Na 22 ( 1 neutr, 0 gamma)  1 Na 23
      v(93) = v(93)*rho
***   1 Na 22 ( 1 neutr, 1 proto)  1 Ne 22
      v(94) = v(94)*rho
***   1 Na 22 ( 1 proto, 0 gamma)  1 Na 23
      v(95) = v(95)*rho*fscr(11)
***   1 Na 23 ( 1 neutr, 0 gamma)  1 Mg 24
      v(96) = v(96)*rho
***   1 Na 23 ( 1 proto, 1 alpha)  1 Ne 20
      v(97) = v(97)*rho*fscr(11)
***   1 Na 23 ( 1 proto, 0 gamma)  1 Mg 24
      v(98) = v(98)*rho*fscr(11)
***   1 Na 23 ( 1 alpha, 1 proto)  1 Mg 26
      v(99) = v(99)*rho*fscr(28)
***   1 Mg 24 ( 1 neutr, 0 gamma)  1 Mg 25
      v(100) = v(100)*rho
***   1 Mg 24 ( 1 proto, 0 gamma)  1 Mg 25
      v(101) = v(101)*rho*fscr(12)
***   1 Mg 24 ( 1 alpha, 1 proto)  1 Al 27
      v(102) = v(102)*rho*fscr(29)
***   1 Mg 24 ( 1 alpha, 0 gamma)  1 Si 28
      v(103) = v(103)*rho*fscr(29)
***   1 Mg 25 ( 1 neutr, 0 gamma)  1 Mg 26
      v(104) = v(104)*rho
***   1 Mg 25 ( 1 proto, 0 gamma)  1 Al26m
      v(105) = v(105)*rho*fscr(12)
***   1 Mg 25 ( 1 proto, 0 gamma)  1 Al26g
      v(106) = v(106)*rho*fscr(12)
***   1 Mg 25 ( 1 alpha, 1 proto)  1 Si 28
      v(107) = v(107)*rho*fscr(29)
***   1 Mg 25 ( 1 alpha, 1 neutr)  1 Si 28
      v(108) = v(108)*rho*fscr(29)
***   1 Mg 25 ( 1 alpha, 0 gamma)  1 Si 29
      v(109) = v(109)*rho*fscr(29)
***   1 Mg 26 ( 1 neutr, 0 gamma)  1 Al 27
      v(110) = v(110)*rho
***   1 Mg 26 ( 1 proto, 0 gamma)  1 Al 27
      v(111) = v(111)*rho*fscr(12)
***   1 Mg 26 ( 1 alpha, 0 gamma)  1 Si 30
      v(112) = v(112)*rho*fscr(29)
***   1 Mg 26 ( 1 alpha, 1 proto)  1 Si 29
      v(113) = v(113)*rho*fscr(29)
***   1 Mg 26 ( 1 alpha, 1 neutr)  1 Si 29
      v(114) = v(114)*rho*fscr(29)
***   1 Al26m ( 0 betap, 0 nutri)  1 Mg 26

***   1 Al26m ( 1 neutr, 1 proto)  1 Mg 26
      v(116) = v(116)*rho
***   1 Al26m ( 1 neutr, 1 alpha)  1 Na 23
      v(117) = v(117)*rho
***   1 Al26m ( 1 neutr, 0 gamma)  1 Al 27
      v(118) = v(118)*rho
***   1 Al26m ( 1 proto, 0 gamma)  1 Al 27
      v(119) = v(119)*rho*fscr(13)
***   1 Al26g ( 0 betap, 0 nutri)  1 Mg 26

***   1 Al26g ( 1 neutr, 0 gamma)  1 Al 27
      v(121) = v(121)*rho
***   1 Al26g ( 1 neutr, 1 proto)  1 Mg 26
      v(122) = v(122)*rho
***   1 Al26g ( 1 neutr, 1 alpha)  1 Na 23
      v(123) = v(123)*rho
***   1 Al26g ( 1 proto, 0 gamma)  1 Al 27
      v(124) = v(124)*rho*fscr(13)
***   1 Al26g ( 0 gamma, 0 gamma)  1 Al26m

***   1 Al 27 ( 1 proto, 0 gamma)  1 Si 28
      v(126) = v(126)*rho*fscr(13)
***   1 Al 27 ( 1 proto, 1 alpha)  1 Mg 24
      v(127) = v(127)*rho*fscr(13)
***   1 Al 27 ( 1 alpha, 1 neutr)  1 Si 30
      v(128) = v(128)*rho*fscr(30)
***   1 Si 28 ( 1 neutr, 0 gamma)  1 Si 29
      v(129) = v(129)*rho
***   1 Si 28 ( 1 proto, 0 gamma)  1 Si 29
      v(130) = v(130)*rho*fscr(14)
***   1 Si 29 ( 1 neutr, 0 gamma)  1 Si 30
      v(131) = v(131)*rho
***   1 Si 29 ( 1 proto, 0 gamma)  1 Si 30
      v(132) = v(132)*rho*fscr(14)
***   1 Si 30 ( 1 neutr, 0 gamma)  1 Si 31
      v(133) = v(133)*rho
***   1 Si 30 ( 1 proto, 0 gamma)  1 P  31
      v(134) = v(134)*rho*fscr(14)
***   1 Si 31 ( 0 betam, 0 nutri)  1 P  31

***   1 Si 31 ( 1 neutr, 0 gamma)  1 Si 32
      v(136) = v(136)*rho
***   1 Si 31 ( 1 proto, 1 neutr)  1 P  31
      v(137) = v(137)*rho*fscr(14)
***   1 Si 31 ( 1 proto, 0 gamma)  1 P  32
      v(138) = v(138)*rho*fscr(14)
***   1 P  31 ( 1 neutr, 1 proto)  1 Si 31
      v(139) = v(139)*rho
***   1 P  31 ( 1 neutr, 0 gamma)  1 P  32
      v(140) = v(140)*rho
***   1 P  31 ( 1 proto, 1 alpha)  1 Si 28
      v(141) = v(141)*rho*fscr(15)
***   1 P  31 ( 1 proto, 0 gamma)  1 S  32
      v(142) = v(142)*rho*fscr(15)
***   1 Si 32 ( 0 betam, 0 nutri)  1 P  32

***   1 Si 32 ( 1 neutr, 0 gamma)  1 P  33
      v(144) = v(144)*rho
***   1 Si 32 ( 1 proto, 1 neutr)  1 P  32
      v(145) = v(145)*rho*fscr(14)
***   1 Si 32 ( 1 proto, 0 gamma)  1 P  33
      v(146) = v(146)*rho*fscr(14)
***   1 P  32 ( 0 betam, 0 nutri)  1 S  32

***   1 P  32 ( 1 neutr, 0 gamma)  1 P  33
      v(148) = v(148)*rho
***   1 P  32 ( 1 neutr, 1 proto)  1 Si 32
      v(149) = v(149)*rho
***   1 P  32 ( 1 proto, 0 gamma)  1 S  33
      v(150) = v(150)*rho*fscr(15)
***   1 P  32 ( 1 proto, 1 alpha)  1 Si 29
      v(151) = v(151)*rho*fscr(15)
***   1 P  32 ( 1 proto, 1 neutr)  1 S  32
      v(152) = v(152)*rho*fscr(15)
***   1 S  32 ( 1 neutr, 1 proto)  1 P  32
      v(153) = v(153)*rho
***   1 S  32 ( 1 neutr, 0 gamma)  1 S  33
      v(154) = v(154)*rho
***   1 S  32 ( 1 neutr, 1 alpha)  1 Si 29
      v(155) = v(155)*rho
***   1 S  32 ( 1 proto, 0 gamma)  1 S  33
      v(156) = v(156)*rho*fscr(16)
***   1 P  33 ( 0 betam, 0 nutri)  1 S  33

***   1 P  33 ( 1 neutr, 0 gamma)  1 S  34
      v(158) = v(158)*rho
***   1 P  33 ( 1 proto, 0 gamma)  1 S  34
      v(159) = v(159)*rho*fscr(15)
***   1 P  33 ( 1 proto, 1 alpha)  1 Si 30
      v(160) = v(160)*rho*fscr(15)
***   1 P  33 ( 1 proto, 1 neutr)  1 S  33
      v(161) = v(161)*rho*fscr(15)
***   1 S  33 ( 1 neutr, 0 gamma)  1 S  34
      v(162) = v(162)*rho
***   1 S  33 ( 1 neutr, 1 alpha)  1 Si 30
      v(163) = v(163)*rho
***   1 S  33 ( 1 neutr, 1 proto)  1 P  33
      v(164) = v(164)*rho
***   1 S  33 ( 1 proto, 0 gamma)  1 S  34
      v(165) = v(165)*rho*fscr(16)
***   1 S  34 ( 1 neutr, 0 gamma)  1 S  35
      v(166) = v(166)*rho
***   1 S  34 ( 1 proto, 0 gamma)  1 Cl 35
      v(167) = v(167)*rho*fscr(16)
***   1 S  35 ( 0 betam, 0 nutri)  1 Cl 35

***   1 S  35 ( 1 neutr, 0 gamma)  1 S  36
      v(169) = v(169)*rho
***   1 S  35 ( 1 proto, 0 gamma)  1 Cl 36
      v(170) = v(170)*rho*fscr(16)
***   1 Cl 35 ( 1 neutr, 0 gamma)  1 Cl 36
      v(171) = v(171)*rho
***   1 Cl 35 ( 1 proto, 0 gamma)  1 heavy
      v(172) = v(172)*rho*fscr(17)
***   1 S  36 ( 1 neutr, 0 gamma)  1 Cl 37
      v(173) = v(173)*rho
***   1 S  36 ( 1 proto, 0 gamma)  1 Cl 37
      v(174) = v(174)*rho*fscr(16)
***   1 Cl 36 ( 0 betam, 0 nutri)  1 heavy

***   1 Cl 36 ( 1 neutr, 0 gamma)  1 Cl 37
      v(176) = v(176)*rho
***   1 Cl 36 ( 1 proto, 0 gamma)  1 Cl 37
      v(177) = v(177)*rho*fscr(17)
***   1 Cl 37 ( 1 neutr, 0 gamma)  1 heavy
      v(178) = v(178)*rho
***   1 Cl 37 ( 1 proto, 0 gamma)  1 heavy
      v(179) = v(179)*rho*fscr(17)
***   1 heavy ( 1 neutr, 0 gamma)  1 captn
      v(180) = (5.11404d0+signtk)*1.44d5*rho

      return
      end

