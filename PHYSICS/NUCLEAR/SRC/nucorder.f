
      PROGRAM NUCORDER

************************************************************************
* Reorder the nuclide list, the reaction network and the subroutine VIT*
************************************************************************

      implicit double precision (a-h,o-z)
      parameter(nisolin = 500,nlreact = 500,nvit = 3000)
      character vitesse(nvit)*80,reactio(nvit)*37,vitlc*1
      character*5 react1(nlreact),react3(nlreact),react5(nlreact),
     &            react7(nlreact),ctemp,isotope(nisolin)
      character*7 refs(nlreact),cref
      common /alpha/ q(nlreact),l1(nlreact),l2(nlreact),l3(nlreact),
     &               l4(nlreact),l5(nlreact),l6(nlreact),l7(nlreact),
     &               l8(nlreact)
      common /beta/ iverif(nisolin),na(nisolin),nz(nisolin),
     &              isotope
      common /delta/ reactio
      common /nucvit/ vitesse
      common /mem/ vitlc(nvit,80),nnc(nvit)
      common /netw/ react1,react3,react5,react7

      open (unit = 10,file = 'nucl.dat')
      open (unit = 20,file = 'nucl.out')
      open (unit = 30,file = 'nucl.err')

***   reading of the nuclide list
      nisot = 0
   10 nisot = nisot+1
      read (10,1000) i,isotope(nisot),na(nisot),nz(nisot)
      if (isotope(nisot).ne.'* end') goto 10
      nisot = nisot-1
***   ordering of the nuclide list
      do i = 1,nisot-1
        do 20 j = 1,nisot-i
          if (na(j).lt.na(j+1)) goto 20
          if ((na(j).eq.na(j+1)).and.(nz(j).lt.nz(j+1))) goto 20
          ctemp = isotope(j)
          isotope(j) = isotope(j+1)
          isotope(j+1) = ctemp
          itemp = na(j)
          na(j) = na(j+1)
          na(j+1) = itemp
          itemp = nz(j)
          nz(j) = nz(j+1)
          nz(j+1) = itemp
   20   continue
      enddo
***   storage of the nuclide list
      do i = 1,nisot
        write (20,2100) i,isotope(i),na(i),nz(i)
      enddo
      write (20,*) '    * end'

***   reading of the nuclear reactions
      nreact = 0
   30 nreact = nreact+1
      read (10,1100) l2(nreact),react1(nreact),l4(nreact),
     &               react3(nreact),l6(nreact),react5(nreact),
     &               l8(nreact),react7(nreact),q(nreact),refs(nreact)
      if (react1(nreact).ne.'**end') goto 30
      nreact = nreact-1
***   putting codes to the nuclear reactions
      do 40 i = 1,nreact
        l1(i) = 0
        l3(i) = 0
        l5(i) = 0
        l7(i) = 0
        do j = 1,nisot
          if (react1(i).eq.isotope(j)) l1(i) = j
          if (react3(i).eq.isotope(j)) l3(i) = j
          if (react5(i).eq.isotope(j)) l5(i) = j
          if (react7(i).eq.isotope(j)) l7(i) = j
        enddo
        if ((l1(i)*l3(i)*l5(i)*l7(i)).ne.0) goto 40
        write (30,*) 'the following reaction uses an undefined isotope:'
        write (30,2200) i,l2(i),react1(i),l4(i),react3(i),l6(i),
     &                  react5(i),l8(i),react7(i),q(i)
   40 continue
***   verification of the nuclide list
      do i = 1,nisot
        iverif(i) = 0
      enddo
      do i = 1,nreact
        iverif(l1(i)) = 1
        iverif(l3(i)) = 1
        iverif(l5(i)) = 1
        iverif(l7(i)) = 1
      enddo
      do 50 i = 1,nisot
        if (iverif(i).ne.0) goto 50
        write (30,1300) isotope(i)
   50 continue
***   searching for of a created isotope with zero initial abundance
      do 60 i = 1,nisot-1
        do j = 1,nreact
          if (l5(j).eq.i) goto 60
          if (l7(j).eq.i) goto 60
        enddo
        write (30,1400) isotope(i)
   60 continue
***   searching for produced isotopes never destroyed
      do i = 1,nisot
        iverif(i) = 0
      enddo
      do i = 1,nreact
        iverif(l5(i)) = 1
        iverif(l7(i)) = 1
      enddo
      do i = 1,nreact
        iverif(l1(i)) = 0
        iverif(l3(i)) = 0
      enddo
      do 70 i = 1,nisot-1
        if (iverif(i).eq.0) goto 70
        write (30,1500) isotope(i)
   70 continue
***   verification of the baryonic number conservation
      na(nisot) = 0
      do 80 i = 1,nreact
        iresidu = l2(i)*na(l1(i))+l4(i)*na(l3(i))-l6(i)*na(l5(i))
     &            -l8(i)*na(l7(i))
        if (iresidu.eq.0) goto 80
        write (30,*) 'the following reaction does not conserve the'
     &               ,' barionic number'
        write (30,2200) i,l2(i),react1(i),l4(i),react3(i),l6(i),
     &                  react5(i),l8(i),react7(i),q(i)
   80 continue
***   ordering of the nuclear reactions
      do i = 1,nreact-1
        do 90 j = 1,nreact-i
          if (l1(j).lt.l1(j+1)) goto 90
          if ((l1(j).eq.l1(j+1)).and.(l3(j).lt.l3(j+1))) goto 90
          ctemp = react1(j)
          react1(j) = react1(j+1)
          react1(j+1) = ctemp
          ctemp = react3(j)
          react3(j) = react3(j+1)
          react3(j+1) = ctemp
          ctemp = react5(j)
          react5(j) = react5(j+1)
          react5(j+1) = ctemp
          ctemp = react7(j)
          react7(j) = react7(j+1)
          react7(j+1) = ctemp
          itemp = l1(j)
          l1(j) = l1(j+1)
          l1(j+1) = itemp
          itemp = l2(j)
          l2(j) = l2(j+1)
          l2(j+1) = itemp
          itemp = l3(j)
          l3(j) = l3(j+1)
          l3(j+1) = itemp
          itemp = l4(j)
          l4(j) = l4(j+1)
          l4(j+1) = itemp
          itemp = l5(j)
          l5(j) = l5(j+1)
          l5(j+1) = itemp
          itemp = l6(j)
          l6(j) = l6(j+1)
          l6(j+1) = itemp
          itemp = l7(j)
          l7(j) = l7(j+1)
          l7(j+1) = itemp
          itemp = l8(j)
          l8(j) = l8(j+1)
          l8(j+1) = itemp
          rtemp = q(j)
          q(j) = q(j+1)
          q(j+1) = rtemp
          cref = refs(j)
          refs(j) = refs(j+1)
          refs(j+1) = cref
   90   continue
      enddo
***   storage of the nuclear reactions
      do i = 1,nreact
        write (20,2000) i,l2(i),react1(i),l4(i),react3(i),l6(i),
     &                  react5(i),l8(i),react7(i),q(i),l1(i),l3(i),
     &                  l5(i),l7(i),refs(i)
      enddo
      write (20,*) '         **end'

***   reading of the subroutine vit
      nligne = 0
  100 nligne = nligne+1
      read (10,1700) vitesse(nligne)
      if (vitesse(nligne)(:27).ne.'***   end of reaction rates')
     &   goto 100
      rewind (20)
      do i = 1,nisot+1
        read (20,1600) reactio(i)
      enddo
      do i = 1,nreact+1
        read (20,1600) reactio(i)
      enddo
***   storage of the beginning of the subroutine vit
      do m = 1,nligne
        read (vitesse(m),2300) (vitlc(m,kk),kk = 1,80)
      enddo
      do m = 1,nligne
        do kk = 80,1,-1
          if (vitlc(m,kk).ne.' ') goto 110
        enddo
  110   nnc (m) = kk
        if (kk.le.0) nnc(m) = 1
      enddo
      i = 1
  120 write (20,2300) (vitlc(i,kk),kk = 1,nnc(i))
      i = i+1
      if (vitesse(i)(1:3).ne.'** ') goto 120
***   ordering and storage of the nuclear reaction rates
      do 160 i = 1,nreact
        j = 0
  130   j = j+1
        if (j.ne.nligne) goto 140
        write (30,1800) reactio(i)
        goto 160
  140   if (reactio(i).ne.vitesse(j)(6:42)) goto 130
  150   if ((vitesse(j)(:8).eq.'      v(').and.(vitesse(j)(13:16)
     &     .eq.') = ')) then
          write (20,1900) vitesse(j)(:8),i,(vitlc(j,kk),kk = 13,nnc(j))
        else
          write (20,2300) (vitlc(j,kk),kk = 1,nnc(j))
        endif
        j = j+1
        if (vitesse(j)(1:2).ne.'**') goto 150
  160 continue
      write (20,1200)

 1000 format (i4,1x,a5,1x,i3,8x,i2)
 1100 format (5x,i4,1x,a5,2x,i2,1x,a5,1x,i2,1x,a5,2x,i2,1x,a5,7x,f7.3,
     &       25x,a7)
 1200 format ('***   end of reaction rates',//,'      return',
     &       /,'      end',/)
 1300 format (1x,a5,' is not used in the reaction network')
 1400 format (1x,a5,' is destroyed but is not present initially ',
     &       'and is not created',/)
 1500 format (1x,a5,' is created but never destroyed')
 1600 format (7x,a37)
 1700 format (a80)
 1800 format (1x,a37,' has not been found',/)
 1900 format (a8,i4,68(a1))
 2000 format (i4,1x,i4,1x,a5,1x,'(',i2,1x,a5,',',i2,1x,a5,')',1x,i2,
     &       1x,a5,' : q = ',f7.3,' mev',3x,4(i4),2x,a10)
 2100 format (i4,1x,a5,1x,i3,'.00000',2x,i2,'.')
 2200 format (i4,1x,i4,1x,a5,1x,'(',i2,1x,a5,',',i2,1x,a5,')',1x,i2,
     &       1x,a5,' : q = ',f7.3,' mev')
 2300 format(80(a1))

      close (10)
      close (20)
      close (30)

      stop
      end

