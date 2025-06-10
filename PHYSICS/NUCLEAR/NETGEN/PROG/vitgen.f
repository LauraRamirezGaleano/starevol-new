
      PROGRAM VITGEN

************************************************************************
*       tables for interpolation in T of nuclear reaction rates        *
************************************************************************

      implicit none

      integer i,j,k,nre,ngrid,ngrid1,nvit
      integer iinf,isup,ii
      double precision t9min,t9max,dt9,t9,tgrid,flag,ln10
      double precision v,vrate,tnuc9,lnS,t,S
      double precision dt8,t8,t8grid
      double precision t1,t2,t3,t4,dt12,dt23,dt34
      integer n12,n23,n34,nmax,nprint

      character cc*40
      double precision t9l,del

      parameter(nre = 185, nvit = 90)
c      parameter(nre = 180, nvit = 250)

      dimension v(nre),flag(nre)
      dimension tnuc9(nvit),vrate(nvit,nre)
      dimension t8grid(nvit)


      ln10 = log(10.d0)

      open (unit = 10,file = 'vitsub.f')
      open (unit = 20,file = 'vit_global.dat')

      ngrid = nvit
      read (20,100)
      do k = 1,nre-1
         read (20,200) cc
         read (20,300) flag(k)
         write(*,*) k,flag(k),cc
         if (flag(k).gt.0) stop 'flag'
         read (20,*) (vrate(j,k), j = 1,ngrid)
      enddo

c 100  format (44(1x,/))
c 100  format (36(1x,/))
 100  format (17(1x,/))
 200  format (/,8x,a40)
 300  format (1x,f13.3)
      
      t8grid(1) = 1.d-2
      t1 = 1.d-2
      t2 = 1.d0
      t3 = 1.d1
      t4 = 5.d1
      n12 = 20
      n23 = 50
      n34 = 20
      nmax = n12+n23+n34
      nprint = int(nmax/5)
      dt12 = (t2-t1)/dble(n12)
      dt23 = (t3-t2)/dble(n23)
      dt34 = (t4-t3)/dble(n34-1)
      t8grid(nmax) = 5.d1
c      dt8 = (log10(t8grid(nvit))-log10(t8grid(1)))/(nvit-1)
c      t8 = log10(t8grid(1))
      do k = 2,n12+1
         t8grid(k) = t8grid(k-1)+dt12
      enddo
      do k = n12+2,n12+n23+1
         t8grid(k) = t8grid(k-1)+dt23
      enddo
      do k = n12+n23+2,nmax
         t8grid(k) = t8grid(k-1)+dt34
      enddo

      do k=1,ngrid
         tnuc9(k) = t8grid(k)*1.d-1
      enddo

      write (10,1000)
      write (10,1200)
      write (10,1000)
      write (10,1100)
      write (10,1300)
      write (10,1100)
      write (10,1000)
      write (10,1400) nvit,nre
      write (10,1500)

      write (10,1501)
      do ii = 1,17
         iinf = 5*ii-4
         isup = iinf+4
         write (10,1502) (tnuc9(i),i = iinf,isup)
      enddo      
      iinf = 5*ii-4
      isup = iinf+4
      write (10,1503) (tnuc9(i),i = iinf,isup)
      write (10,1000)

c S is the conversion factor: sig(mbarn)->Na<sig v>
c      S=7.7629d+04*dsqrt(t8grid(k))

      do i = 1,ngrid
         t = tnuc9(i)*1.d9
         S = 7.7629d+04*dsqrt(10.d0*tnuc9(i))
         do k = 1,nre-1
            if (flag(k).eq.0..or.flag(k).eq.-12.) v(k) = log(vrate(i,k))
            if (flag(k).eq.-1.) v(k) = log(S*vrate(i,k))
            if (flag(k).eq.-10..or.flag(k).eq.-100.or.flag(k).eq.-14.
     &           .or.flag(k).eq.-11.) v(k) = vrate(i,k)*ln10
            if (v(k).lt.-99.d0) v(k) = -99.d0
         enddo         
c neutron sink
         v(nre) = 13.60745264d0
         write (10,1600) t
         write (10,1700) i,i
         write (10,1800) (v(j),j = 1,7)
         write (10,1800) (v(j),j = 8,14)
         write (10,1800) (v(j),j = 15,21)
         write (10,1800) (v(j),j = 22,28)
         write (10,1800) (v(j),j = 29,35)
         write (10,1800) (v(j),j = 36,42)
         write (10,1800) (v(j),j = 43,49)
         write (10,1800) (v(j),j = 50,56)
         write (10,1800) (v(j),j = 57,63)
         write (10,1800) (v(j),j = 64,70)
         write (10,1800) (v(j),j = 71,77)
         write (10,1800) (v(j),j = 78,84)
         write (10,1800) (v(j),j = 85,91)
         write (10,1800) (v(j),j = 92,98)
         write (10,1800) (v(j),j = 99,105)
         write (10,1800) (v(j),j = 106,112)
         write (10,1800) (v(j),j = 113,119)
         write (10,1800) (v(j),j = 120,126)
         write (10,1800) (v(j),j = 127,133)
         write (10,1800) (v(j),j = 134,140)
         write (10,1800) (v(j),j = 141,147)
         write (10,1800) (v(j),j = 148,154)
         write (10,1800) (v(j),j = 155,161)
         write (10,1800) (v(j),j = 162,168)
         write (10,1800) (v(j),j = 169,175)
c         write (10,1900) (v(j),j = 176,180)
         write (10,1800) (v(j),j = 176,182)
         write (10,1900) (v(j),j = 183,185)
      enddo
      write (10,1000)
      write (10,2000)
      write (10,1000)

      close (10)
      close (20)

 1000 format (1x)
 1100 format (72('*'))
 1200 format (6x,'BLOCK DATA VITSUB')
 1300 format ('*',11x,'tables of nuclear reaction rates in ',
     &     'function of T',10x,'*',/,'*',17x,'(see starevol.par for ',
     &     'the references)',16x,'*')
 1400 format (6x,'implicit none',//,
     &     6x,'integer i,j',/,6x,'integer nvit,nre',/,
     &     6x,'double precision tnuc9,vrate',//,
     &     6x,'parameter (nvit = ',i3,', nre = ',i3,')')
 1500 format (6x,'common /nucvit/ tnuc9(nvit),vrate(nvit,nre)')
 1501 format (1x,/,6x,'data (tnuc9(i),i = 1,nvit) /')
 1502 format (5x,'& ',5(1pd12.6,','))
 1503 format (5x,'& ',4(1pd12.6,','),1pd12.6,'/')
 1600 format (1x,/,'***** nav * <sig*v> / rho  at T = ',d10.4,' K')
 1700 format (6x,'data ((vrate(i,j),i = ',i3,',',i3,'),j = 1,nre) /')
 1800 format (5x,'& ',7(f8.4,','))
 1900 format (5x,'& ',2(f8.4,','),f8.4,' /')
c 1900 format (5x,'& ',4(f8.4,','),f8.4,' /')
 2000 format (6x,'end',/)

      end

