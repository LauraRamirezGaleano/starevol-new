      PROGRAM TGRID

      implicit none
      
      double precision dt8,t8,t8grid
      double precision t1,t2,t3,t4,dt12,dt23,dt34
      integer nvit,i,k,iinf,isup
      integer n12,n23,n34,nmax,nprint

      parameter (nvit = 200)
      dimension t8grid(nvit)

      t8grid(1) = 1.d-2
      t1 = 1.d-2
      t2 = 1.d0
      t3 = 1.d1
      t4 = 1.d2
      n12 = 20
      n23 = 35
      n34 = 35
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

      open (unit = 12,file = 'Tgrid.in')
      write(12,100) nmax
      do k = 1,nprint
         iinf = 5*k-4
         isup = iinf+4
         write(12,200) (t8grid(i),i = iinf,isup) 
      enddo
  
 100  format (i3,/,' 1. 1.')
 200  format (5(1x,1pd12.6))

      end
