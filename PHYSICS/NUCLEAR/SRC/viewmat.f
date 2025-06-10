
      PROGRAM VIEWMAT

************************************************************************
* View of the nuclear network matrix non-zero elements                 *
************************************************************************

      parameter (nsp = 55,nre = 180)
      integer*2 l1,l3,l5,l7,jp
      character*80 ifile,ofile
      character*1 a
      common /a/ l1(nre),l3(nre),l5(nre),l7(nre),ibmin(10),ibmax(10)
      common /b/ a(nsp,nsp)
      common /c/ jp(nsp)

      open (unit = 10,file = 'nuc.netw')
      open (unit = 11,file = 'matrix.nuc')

      do i = 1,nsp
        kkk = kkk+1
        if (kkk.eq.10) kkk = 0
        jp(i) = kkk
      enddo
      do i = 1,nsp
        do j = 1,nsp
          a(j,i)='.'
        enddo
      enddo

      do i = 1,nre
        read (10,1000) l1(i),l3(i),l5(i),l7(i)
      enddo

      do m = 1,nre
        ii=l1(m)
        ij=l3(m)
        ik=l5(m)
        il=l7(m)
        a(ii,ii)='*'
        a(ij,ii)='*'
        a(ij,ij)='*'
        a(ii,ij)='*'
        a(ii,ik)='*'
        a(ij,ik)='*'
        a(ii,il)='*'
        a(ij,il)='*'
      enddo
      do k = 1,nsp
        a(k,k) = '@'
      enddo

      write (11,1100)
      nbl = int(nsp/128+1)
      ibmin(nbl) = nsp-128+1
      ibmin(1) = 1
      ibmax(1) = 128
      ibmax(nbl) = nsp
      if (nbl.le.2) goto 10
      do j = 2,nbl-1
        ibmin(j) = ibmin(j-1)+128
        ibmax(j) = ibmax(j-1)+128
      enddo
   10 do j = 1,nbl
        write (11,1200)(jp(i),i = ibmin(j),ibmax(j))
        do i = ibmin(j),ibmax(j)
          write (11,1300) i,(a(k,i),k = ibmin(j),ibmax(j))
        enddo
      enddo

 1000 format (65x,4(i4))
 1100 format (/,1x,'matrix morphology of the nuclear reaction ',
     &       'network',//)
 1200 format (/,5x,128(i1),/)
 1300 format (1x,i3,1x,128(a1))

      close (10)
      close (11)

      stop

      end

