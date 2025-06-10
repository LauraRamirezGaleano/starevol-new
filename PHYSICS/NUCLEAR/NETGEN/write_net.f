      program write_net

      implicit none

      integer i,j,nsp,icode,ireac,nreac,nsave,coef,fs

      double precision anuc,znuc,zcharge,mass,itype,qinu

      character species*5,creac*5

      dimension species(100),icode(100),anuc(100),znuc(100),creac(4),
     &     ireac(4),coef(4),itype(4)

      open (unit=99,file='starevol.par')
      open (unit=88,file='vit_temp.f')

      read(99,10)
 10   format(32/)

      nsp = 0

 100  nsp = nsp+1
      read(99,111) icode(nsp),species(nsp),anuc(nsp),znuc(nsp)
 111  format(i4,1x,a5,1x,f9.5,1x,f4.0)
      if (icode(nsp).lt.icode(nsp-1)) goto 300
      goto 100

 300  species(nsp) = 'nutri'
      icode(nsp) = 0
      anuc(nsp) = 0.d0
      znuc(nsp) = 0.d0
      nsp = nsp+1
      species(nsp) = 'betap'
      icode(nsp) = 0
      anuc(nsp) = 0.d0
      znuc(nsp) = 1.d0
      nsp = nsp+1
      species(nsp) = 'betam'
      icode(nsp) = 0
      anuc(nsp) = 0.d0
      znuc(nsp) = -1.d0
      nsave = 0

 400  read(99,450,end=1000) nreac,coef(1),creac(1),coef(2),creac(2)
     &     ,coef(3),creac(3),coef(4),creac(4),qinu
 450  format(i4,4x,i1,1x,a5,3x,i1,1x,a5,2x,i1,1x,a5,3x,i1,1x,a5,46x,
     &     f7.3)

      do 350 i=1,4
         do j=1,nsp
            if (creac(i).eq.species(j)) then
               ireac(i) = icode(j)
               itype(i) = j
               goto 350
            endif
         enddo
 350  continue
      fs = int(coef(1)*znuc(itype(1))+17.d0*coef(2)*(znuc(itype(2))-1))
      if (coef(1).eq.2.and.coef(2).eq.0) then         
         fs = int(znuc(itype(1))+17*int(znuc(itype(1))-1.d0))
      endif
      if (znuc(itype(1)).eq.6.d0.and.coef(1).eq.2) fs = 35
      if (znuc(itype(1)).eq.8.d0.and.coef(1).eq.2) fs = 36  
      if (coef(2).eq.0.and.coef(1).eq.1) then
         write (88,480) coef(1),creac(1),coef(2),creac(2),coef(3),
     &        creac(3),coef(4),creac(4)
      else
         if (itype(2).eq.1) then
            write (88,486) coef(1),creac(1),coef(2),creac(2),coef(3),
     &           creac(3),coef(4),creac(4),nreac,nreac
         else
            if (coef(1).eq.3) fs = -99
            write (88,485) coef(1),creac(1),coef(2),creac(2),coef(3),
     &           creac(3),coef(4),creac(4),nreac,nreac,fs
         endif
      endif
 480  format('***',3x,i1,1x,a5,' ( ',i1,1x,a5,', ',i1,1x,a5,')  ',i1,1x,
     &     a5,/)
 485  format('***',3x,i1,1x,a5,' ( ',i1,1x,a5,', ',i1,1x,a5,')  ',i1,1x,
     &     a5,/,6x,'v(',i3,') = v(',i3,')*rok*fscr(ksh,',i3,')')
 486  format('***',3x,i1,1x,a5,' ( ',i1,1x,a5,', ',i1,1x,a5,')  ',i1,1x,
     &     a5,/,6x,'v(',i3,') = v(',i3,')*rok')
      if (nreac.gt.nsave) then
         zcharge = coef(1)*znuc(itype(1))+coef(2)*znuc(itype(2))-coef(3)
     &        *znuc(itype(3))-coef(4)*znuc(itype(4))
         mass = coef(1)*anuc(itype(1))+coef(2)*anuc(itype(2))-coef(3)
     &        *anuc(itype(3))-coef(4)*anuc(itype(4))
c..   test mass, charge conservation, qinu
c         if ((zcharge.ne.0.d0.or.mass.ne.0.d0).and.qinu.eq.0.d0) 
cc         if ((zcharge.ne.0.d0.or.mass.ne.0.d0)) 
c     &        write(*,700) nreac,ireac(1),ireac(2),ireac(3),ireac(4)
c     &        ,int(mass),int(zcharge),species(itype(2)),species(itype(3)
c     &        ),qinu
         write(*,500) ireac(1),ireac(2),ireac(3),ireac(4)
 500     format(4(i4))
 700     format(i4,6(1x,i3),3x,a5,',',a5,f7.3)
         nsave = nreac
         goto 400
      endif

 1000 continue
      end
