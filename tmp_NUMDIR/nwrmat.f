
************************************************************************

      SUBROUTINE nwrmat (a,b,nm,mm,error)

************************************************************************
* Calculate the Henyey matrix for the stellar structure                *
* a(n,m) <-> b(n*m)  ==> b((j-1)*n+i) = a(i,j)
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      integer nm,mm
      integer error
      integer npm,nj,jj,j1,jm,ij,i1,i2,ik,jk
      integer i,j,k

      double precision a,b
      double precision amax,zwi,factor

      dimension a(nm*(nm+mm)),b(nm*mm)

      npm = nm+mm
      do j = 1,nm
        nj = (j-1)*nm
        jj = nj+j
        j1 = j+1
        amax = abs(a(jj))
        jm = j
        if (j1.le.nm) then
          do i = j1,nm
            ij = nj + i
            if (abs(a(ij)).gt.amax) then
              amax = abs (a(ij))
              jm = i
            endif
          enddo
          if (jm.lt.j) goto 40
          if (jm.gt.j) then
            i1 = jm + nj
            i2 = jj
            do i = j,npm
              zwi = a(i1)
              a(i1) = a(i2)
              a(i2) = zwi
              i1 = i1+nm
              i2 = i2+nm
            enddo
          endif
        endif
        if (a(jj).eq.0.d0) goto 40
        do i = 1,nm
          if (i.ne.j) then
            ij = nj+i
            ik = nj+i
            jk = jj
            factor = -a(ij)/a(jj)
            do k = j1,npm
              jk = jk+nm
              ik = ik+nm
              a(ik) = a(ik)+factor*a(jk)
            enddo
          endif
        enddo
        jk = jj
        factor = 1.d0/a(jj)
        do k = j1,npm
          jk = jk+nm
          a(jk) = a(jk)*factor
        enddo
      enddo

      i1 = nm*mm
      i2 = nm*nm
      do i = 1,i1
        i2 = i2+1
        b(i) = a(i2)
      enddo
      return

   40 error = 5
      i = jj-nm*(int(jj/nm))
      j = int(jj/nm)+1
      if (i.eq.0) then
         i = nm
         j = j-1
      endif
      print *,'zero in matrix a[',nm,'x (',nm,'+',mm,')] at ',
     &     'position : ',jj,'which corresponds to position (',
     &     i,',',j,')'


      return
      end
