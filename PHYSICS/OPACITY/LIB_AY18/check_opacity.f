      program checkopa

      implicit none

      integer i,j,k,km,l,nn,al
      
      double precision tk,rk,kappa(76,39,10),rklth(39),tklth(76)
      double precision x(10),y(10),z(10)
      character*255 opalib
      character*300 dum,ttt

      call getenv('DIR_OPAAY18',opalib)

      open(unit=110,file=trim(opalib) // 'AY18hz',status = 'old'
     $     ,form = 'formatted')

      
      data (rklth(km),km = 1,39)/ -9.000d0, -8.500d0, -8.000d0,-7.500d0,
     $     -7.000d0,-6.500d0,-6.000d0,-5.500d0,-5.000d0,-4.500d0,
     $     -4.000d0,-3.500d0,-3.000d0,-2.500d0,-2.000d0,-1.500d0,
     $     -1.000d0,-0.500d0, 0.000d0, 0.500d0,1.000d0,1.500d0,2.000d0
     $     ,2.500d0,3.000d0,3.500d0,4.000d0,4.500d0, 5.000d0,5.500d0
     $     ,6.000d0,6.500d0,7.000d0,7.500d0,8.000d0,8.500d0,9.000d0,
     $     9.500d0, 10.000d0 /

      data (tklth(j),j=1,76)/3.7500d0,  3.800d0 ,  3.8500d0,  3.90000d0
     $     ,  3.9500d0,  4. 0000d0 ,  4.0500d0,  4.10000d0 ,  4.1500d0,
     $     4.20000d0 ,  4.2500d0,  4.30000d0 ,  4.3500d0,  4.40000d0 ,
     $     4.4500d0,  4.50000d0 ,  4.5500d0,  4.60000d0, 4.6500d0,
     $     4.70000d0 ,  4.7500d0,  4.80000d0 ,  4.8500d0,  4.90000d0 ,
     $     4.9500d0,  5.0000d0  ,  5.0500d0, 5.10000d0 ,  5.1500d0,
     $     5.20000d0 ,  5.2500d0,  5.3000d0 ,  5.3500d0,  5.40000d0 ,
     $     5.4500d0,  5.50000d0 , 5.5500d0,  5.60000d0 ,  5.6500d0,
     $     5.70000d0 ,  5.7500d0,  5.80000d0 ,  5.8500d0,  5.90000d0 ,
     $     5.9500d0, 6.0000d0 ,  6.10000d0,  6.20000d0,  6.30000d0,
     $     6.40000d0,  6.50000d0,  6.60000d0,  6.70000d0,  6.80000d0,
     $     6.90000d0,  7.0000d0 ,7.10000d0,  7.20000d0,  7.30000d0,
     $     7.40000d0,  7.50000d0,  7.60000d0,  7.70000d0,  7.80000d0,
     $     7.90000d0,  8.0000d0 ,  8.10000d0, 8.30000d0,  8.50000d0,
     $     8.70000d0,  8.90000d0,  9.10000d0,  9.30000d0,  9.50000d0,
     $     9.70000d0, 9.9000d0 /


      l = 1
      read (110,'(738/)')
      read(110,'(36x,3(f6.4,3x))') x(l),y(l),z(l)
      read(110,'(4/)')
      do i = 1,76
         read(110,'(5x,39f7.3)')(kappa(i,k,l),k=1,39)
      enddo
      l = l+1
      rewind(110)
      read (110,'(1817/)')
      read(110,'(36x,3(f6.4,3x))') x(l),y(l),z(l)
      read(110,'(4/)')
      do i = 1,76
         read(110,'(5x,39f7.3)')(kappa(i,k,l),k=1,39)
      enddo
      l =l+1
      rewind(110)
      read (110,'(2896/)')
      read(110,'(36x,3(f6.4,3x))') x(l),y(l),z(l)
      read(110,'(4/)')
      do i = 1,76
         read(110,'(5x,39f7.3)')(kappa(i,k,l),k=1,39)
      enddo
      l =l+1
      rewind(110)
      read (110,'(3975/)')
      read(110,'(36x,3(f6.4,3x))')x(l),y(l),z(l) 
      read(110,'(4/)')
      do i = 1,76
         read(110,'(5x,39f7.3)')(kappa(i,k,l),k=1,39)
      enddo
      l =l+1
      rewind(110)
      read (110,'(5054/)')
      read(110,'(36x,3(f6.4,3x))') x(l),y(l),z(l)
      read(110,'(4/)')
      do i = 1,76
         read(110,'(5x,39f7.3)')(kappa(i,k,l),k=1,39)
      enddo
      l =l+1
      rewind(110)
      read (110,'(6133/)')
      read(110,'(36x,3(f6.4,3x))') x(l),y(l),z(l)
      read(110,'(4/)')
      do i = 1,76
         read(110,'(5x,39f7.3)')(kappa(i,k,l),k=1,39)
      enddo
      l =l+1
      rewind(110)
      read (110,'(7212/)')
      read(110,'(36x,3(f6.4,3x))') x(l),y(l),z(l)
      read(110,'(4/)')
      do i = 1,76
         read(110,'(5x,39f7.3)')(kappa(i,k,l),k=1,39)
      enddo
      l =l+1
      rewind(110)
      read (110,'(8291/)')
      read(110,'(36x,3(f6.4,3x))')x(l),y(l),z(l)
      read(110,'(4/)')
      do i = 1,76
         read(110,'(5x,39f7.3)')(kappa(i,k,l),k=1,39)
      enddo
      l =l+1
      rewind(110)
      read (110,'(9370/)')
      read(110,'(36x,3(f6.4,3x))') x(l),y(l),z(l)
      read(110,'(4/)')
      do i = 1,76
         read(110,'(5x,39f7.3)')(kappa(i,k,l),k=1,39)
      enddo
      l =l+1
      rewind(110)
      read (110,'(10200/)')
      read(110,'(36x,3(f6.4,3x))') x(l),y(l),z(l)
      read(110,'(4/)')
      do i = 1,76
         read(110,'(5x,39f7.3)')(kappa(i,k,l),k=1,39)
      enddo
      rewind(110)

      al = 700
      do l = 1,10
         write(al,'(3(f6.4,3x))') x(l),y(l),z(l)
         do i = 1,76
            write(al,'(f5.3,39f7.3)')tklth(i),(kappa(i,k,l),k=1,39)
         enddo
         al = al+1 
      enddo
      end
           
