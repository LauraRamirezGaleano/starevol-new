	program ferg04bext
        implicit double precision (a-h,o-z)
c
c	purpose:
c		read fergusson 2004 opacity tables (Asplund et al. 2005 solar composition) 
c               and convert values from ascii to binary format .
c
c
c       History: 
c
c       18/08/97:  creation from s/r alex95bext.f
c
c       last modification: 21/08/97
c
c	Variables:
c	lun ....... logical unit number of opacity-table-file
c	ntab ...... nr. of different x-values
c	nval ...... max. nr. of op.values pro line (= # of rlg-values)
c	nlin ...... max. nr. of op.table lines (= # of log(T)=values)
c       nzva ...... nr. of different z-values
c                   (z = mass fraction of heavy elements)
c
c	rlg ....... array[1..nval,1..ntab], 
c		    decade log of r=density(gm/cm**3)/t6**3
c	tlg ....... array[1..nlin,1..ntab].
c                   decade log of temperature
c	opa ....... array[1..nvar,1..nlin,1..ntab],
c                   opacity values
c
c       setup array-dimensions
c
        parameter(ntab=23,ntab1=8,nzva=16)
        parameter(nval=19,nlin=85)
        parameter(lun=31)
c
        dimension       ir (     ntab,nzva)
	dimension	rlg(nval,ntab,nzva)
	dimension	tlg(nlin,ntab,nzva)
	dimension	opa(nval,nlin,ntab,nzva)
        dimension       xtab(ntab),ztab(nzva)
        dimension       xi(2),di(2)

	character*132 	line
	character*21	ifname(ntab,nzva),ofname
	character*5     Xvalues(ntab1),Zvalues(nzva),tt
c      initialize values of X
	data Xvalues /'0','1','2','35','5','7','8','9'/
	
c      initialize values of Z
	data Zvalues /'0','00001','00003','0001','0003','001','002','004',
     &    '01','02','03','04','05','06','08','1'/


c
c-----initialize output filename and table
      data ofname
     .     /'ferg04.bin'/

	do i = 1,nzva
	   do j = 1,ntab
	      do k = 1,nlin
		 tlg(k,j,i) = -99.d0
		 do l = 1,nval
		    opa(l,k,j,i) = -99.d0
		 enddo
	      enddo
	      do l = 1,nval
		 rlg(l,j,i) = -99.d0
	      enddo
	   enddo
	enddo
c
c     open outputfile
      lun1=lun+1
      open(lun1,file=ofname,status='unknown',
     .          form='unformatted',err=9011)
c
c
c
      do k=1,ntab1
c         open inputfiles
c---------read tables
	 do l=1,nzva
	    ifname(k,l) = 'ags04.' // trim(Xvalues(k)) // '.' // 
     &     trim(Zvalues(l)) // '.tron'
	     open(lun,file=ifname(k,l),status='old',err=9001)
	     print *,lun,ifname(k,l),'k',k
c-----------read 1st line 
c	     read(lun,'(a)',err=9002) line
	     read(lun,'(43x,f8.6,8x,f8.6)')xtab(k),ztab(l)
c       read blank line
	      read(lun,'(40x,a5)')tt
c       read blank line
	      read(lun,'(a)') line

c-----------read data tlg,rlg and opacity-values
c	    read/write string 'logT' + rlg-values
	     read(lun,   '(a6,19f7.3)') ct6,(rlg(i,k,l),i=1,nval)
	     write(lun+1              )     (rlg(i,k,l),i=1,nval)
c-----------read/write logT & opacity-values
	     do j=1,nlin
		read(lun,'(f5.3,2x,f6.3,18(1x,f6.3))',err=9002)
     &           tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
		write(lun+1,             err=9020)
     &           tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	     enddo
	  enddo
          close(lun)
	enddo

c       Table X = 0.92, Z = 0.08
c	k = k+1
	print *,'k',k,xtab(k)
	l = nzva-1
	open(lun,file='ags04.92.08.tron',status='old',err=9001)
	read(lun,'(a)',err=9002) line
	read(line,'(43x,f8.6,8x,f8.6)')xtab(k),ztab(l)
	     write(*,'(a4,1x,43x,f8.6,8x,f8.6)')'tutu',xtab(k),ztab(l)
c       read blank line
	      read(lun,'(a)') line
c       read blank line
	      read(lun,'(a)') line

c-----------read data tlg,rlg and opacity-values
c	    read/write string 'logT' + rlg-values
	read(lun,   '(a6,19f7.3)') ct6,(rlg(i,k,l),i=1,nval)
	write(lun1              )     (rlg(i,k,l),i=1,nval)
c-----------read/write logT & opacity-values
	do j=1,nlin
	read(lun,'(f5.3,2x,f6.3,18(1x,f6.3))',err=9002)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	   write(lun1,             err=9020)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	enddo
	close(lun)

c	print *,'k X = 1.0',k
c       Table X = 0.94, Z = 0.06

	k = k+1
	l = nzva-2
	open(lun,file='ags04.94.06.tron',status='old',err=9001)
	read(lun,'(a)',err=9002) line
	read(line,'(43x,f8.6,8x,f8.6)')xtab(k),ztab(l)
	     write(*,'(a4,1x,43x,f8.6,8x,f8.6)')'tutu',xtab(k),ztab(l)
c       read blank line
	      read(lun,'(a)') line
c       read blank line
	      read(lun,'(a)') line

	print *,'k',k,xtab(k)
c-----------read data tlg,rlg and opacity-values
c	    read/write string 'logT' + rlg-values
	read(lun,   '(a6,19f7.3)') ct6,(rlg(i,k,l),i=1,nval)
	write(lun1              )     (rlg(i,k,l),i=1,nval)
c-----------read/write logT & opacity-values
	do j=1,nlin
	read(lun,'(f5.3,2x,f6.3,18(1x,f6.3))',err=9002)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	   write(lun1,             err=9020)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	enddo
	close(lun)
	
c	print *,'k X = 1.0',k
c       Tables @ X = 0.95

      k = k+1
c         open inputfiles
c---------read tables
	 do l=1,nzva-3
	    ifname(k,l) = 'ags04.95.' // trim(Zvalues(l)) // '.tron'
	     print *,ifname(k,l)
	     open(lun,file=ifname(k,l),status='old',err=9001)
c-----------read 1st line 
	     read(lun,'(a)',err=9002) line
	     read(line,'(43x,f8.6,8x,f8.6)')xtab(k),ztab(l)
c       read blank line
	      read(lun,'(a)') line
c       read blank line
	      read(lun,'(a)') line

	print *,'k',k,xtab(k)
c-----------read data tlg,rlg and opacity-values
c	    read/write string 'logT' + rlg-values
	     read(lun,   '(a6,19f7.3)') ct6,(rlg(i,k,l),i=1,nval)
	     write(lun+1              )     (rlg(i,k,l),i=1,nval)
c-----------read/write logT & opacity-values
	     do j=1,nlin
		read(lun,'(f5.3,2x,f6.3,18(1x,f6.3))',err=9002)
     &           tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
		write(lun+1,             err=9020)
     &           tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	     enddo
	     close(lun)
	  enddo

c	print *,'k X = 1.0',k
c       Table X = 0.96, Z = 0.04
	k = k+1
	l = nzva-3
	open(lun,file='ags04.96.04.tron',status='old',err=9001)
	read(lun,'(a)',err=9002) line
	read(line,'(43x,f8.6,8x,f8.6)')xtab(k),ztab(l)
c       read blank line
	      read(lun,'(a)') line
c       read blank line
	      read(lun,'(a)') line

	print *,'k',k,xtab(k)
c-----------read data tlg,rlg and opacity-values
c	    read/write string 'logT' + rlg-values
	read(lun,   '(a6,19f7.3)') ct6,(rlg(i,k,l),i=1,nval)
	write(lun1              )     (rlg(i,k,l),i=1,nval)
c-----------read/write logT & opacity-values
	do j=1,nlin
	   read(lun,'(f5.3,2x,f6.3,18(1x,f6.3))',err=9002)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	   write(lun1,             err=9020)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	enddo
	close(lun)

c	print *,'k X = 1.0',k
c       Table X = 0.97, Z = 0.03

	k = k+1
	l = l-1
	open(lun,file='ags04.97.03.tron',status='old',err=9001)
	read(lun,'(a)',err=9002) line
	read(line,'(43x,f8.6,8x,f8.6)')xtab(k),ztab(l)
c       read blank line
	      read(lun,'(a)') line
c       read blank line
	      read(lun,'(a)') line

	print *,'k',k,xtab(k)
c-----------read data tlg,rlg and opacity-values
c	    read/write string 'logT' + rlg-values
	read(lun,   '(a6,19f7.3)') ct6,(rlg(i,k,l),i=1,nval)
	write(lun1              )     (rlg(i,k,l),i=1,nval)
c-----------read/write logT & opacity-values
	do j=1,nlin
	   read(lun,'(f5.3,2x,f6.3,18(1x,f6.3))',err=9002)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	   write(lun1,             err=9020)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	enddo
	close(lun)

c	print *,'k X = 1.0',k
c       Table X = 0.98, Z = 0.02

	k = k+1
	l = l-1
	open(lun,file='ags04.98.02.tron',status='old',err=9001)
	read(lun,'(a)',err=9002) line
	read(line,'(43x,f8.6,8x,f8.6)')xtab(k),ztab(l)
c       read blank line
	      read(lun,'(a)') line
c       read blank line
	      read(lun,'(a)') line

	print *,'k',k,xtab(k)
c-----------read data tlg,rlg and opacity-values
c	    read/write string 'logT' + rlg-values
	read(lun,   '(a6,19f7.3)') ct6,(rlg(i,k,l),i=1,nval)
	write(lun1              )     (rlg(i,k,l),i=1,nval)
c-----------read/write logT & opacity-values
	do j=1,nlin
	   read(lun,'(f5.3,2x,f6.3,18(1x,f6.3))',err=9002)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	   write(lun1,             err=9020)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	enddo
	close(lun)

c	print *,'k X = 1.0',k
c       Table X = 0.99, Z = 0.01

	k = k+1
	l = l-1
	open(lun,file='ags04.99.01.tron',status='old',err=9001)
	read(lun,'(a)',err=9002) line
	read(line,'(43x,f8.6,8x,f8.6)')xtab(k),ztab(l)
c       read blank line
	      read(lun,'(a)') line
c       read blank line
	      read(lun,'(a)') line

	print *,'k',k,xtab(k)
c-----------read data tlg,rlg and opacity-values
c	    read/write string 'logT' + rlg-values
	read(lun,   '(a6,19f7.3)') ct6,(rlg(i,k,l),i=1,nval)
	write(lun1              )     (rlg(i,k,l),i=1,nval)
c-----------read/write logT & opacity-values
	do j=1,nlin
	   read(lun,'(f5.3,2x,f6.3,18(1x,f6.3))',err=9002)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	   write(lun1,             err=9020)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	enddo
	close(lun)

c	print *,'k X = 1.0',k
c       Table X = 0.996, Z = 0.004

	k = k+1
	l = l-1
	open(lun,file='ags04.996.004.tron',status='old',err=9001)
	read(lun,'(a)',err=9002) line
	read(line,'(43x,f8.6,8x,f8.6)')xtab(k),ztab(l)
c       read blank line
	      read(lun,'(a)') line
c       read blank line
	      read(lun,'(a)') line

	print *,'k',k,xtab(k)
c-----------read data tlg,rlg and opacity-values
c	    read/write string 'logT' + rlg-values
	read(lun,   '(a6,19f7.3)') ct6,(rlg(i,k,l),i=1,nval)
	write(lun1              )     (rlg(i,k,l),i=1,nval)
c-----------read/write logT & opacity-values
	do j=1,nlin
	   read(lun,'(f5.3,2x,f6.3,18(1x,f6.3))',err=9002)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	   write(lun1,             err=9020)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	enddo
	close(lun)

c	print *,'k X = 1.0',k
c       Table X = 0.998, Z = 0.002

	k = k+1
	l = l-1
	open(lun,file='ags04.998.002.tron',status='old',err=9001)
	read(lun,'(a)',err=9002) line
	read(line,'(43x,f8.6,8x,f8.6)')xtab(k),ztab(l)
c       read blank line
	      read(lun,'(a)') line
c       read blank line
	      read(lun,'(a)') line

	print *,'k',k,xtab(k)
c-----------read data tlg,rlg and opacity-values
c	    read/write string 'logT' + rlg-values
	read(lun,   '(a6,19f7.3)') ct6,(rlg(i,k,l),i=1,nval)
	write(lun1              )     (rlg(i,k,l),i=1,nval)
c-----------read/write logT & opacity-values
	do j=1,nlin
	   read(lun,'(f5.3,2x,f6.3,18(1x,f6.3))',err=9002)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	   write(lun1,             err=9020)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	enddo
	close(lun)

c	print *,'k X = 1.0',k
c       Table X = 0.999, Z = 0.001

	k = k+1
	l = l-1
	open(lun,file='ags04.999.001.tron',status='old',err=9001)
	read(lun,'(a)',err=9002) line
	read(line,'(43x,f8.6,8x,f8.6)')xtab(k),ztab(l)
c       read blank line
	      read(lun,'(a)') line
c       read blank line
	      read(lun,'(a)') line

	print *,'k',k,xtab(k)
c-----------read data tlg,rlg and opacity-values
c	    read/write string 'logT' + rlg-values
	read(lun,   '(a6,19f7.3)') ct6,(rlg(i,k,l),i=1,nval)
	write(lun1              )     (rlg(i,k,l),i=1,nval)
c-----------read/write logT & opacity-values
	do j=1,nlin
	   read(lun,'(f5.3,2x,f6.3,18(1x,f6.3))',err=9002)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	   write(lun1,             err=9020)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	enddo
	close(lun)

c	print *,'k X = 1.0',k
c       Table X = 0.9997, Z = 0.0003

	k = k+1
	l = l-1
	open(lun,file='ags04.9997.0003.tron',status='old',err=9001)
	read(lun,'(a)',err=9002) line
	read(line,'(43x,f8.6,8x,f8.6)')xtab(k),ztab(l)
c       read blank line
	      read(lun,'(a)') line
c       read blank line
	      read(lun,'(a)') line

	print *,'k',k,xtab(k)
c-----------read data tlg,rlg and opacity-values
c	    read/write string 'logT' + rlg-values
	read(lun,   '(a6,19f7.3)') ct6,(rlg(i,k,l),i=1,nval)
	write(lun1              )     (rlg(i,k,l),i=1,nval)
c-----------read/write logT & opacity-values
	do j=1,nlin
	   read(lun,'(f5.3,2x,f6.3,18(1x,f6.3))',err=9002)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	   write(lun1,             err=9020)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	enddo
	close(lun)


c	print *,'k X = 1.0',k
c       Table X = 0.9999, Z = 0.0001

	k = k+1
	l = l-1
	open(lun,file='ags04.9999.0001.tron',status='old',err=9001)
	read(lun,'(a)',err=9002) line
	read(line,'(43x,f8.6,8x,f8.6)')xtab(k),ztab(l)
c       read blank line
	      read(lun,'(a)') line
c       read blank line
	      read(lun,'(a)') line

	print *,'k',k,xtab(k)
c-----------read data tlg,rlg and opacity-values
c	    read/write string 'logT' + rlg-values
	read(lun,   '(a6,19f7.3)') ct6,(rlg(i,k,l),i=1,nval)
	write(lun1              )     (rlg(i,k,l),i=1,nval)
c-----------read/write logT & opacity-values
	do j=1,nlin
	   read(lun,'(f5.3,2x,f6.3,18(1x,f6.3))',err=9002)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	   write(lun1,             err=9020)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	enddo
	close(lun)

c	print *,'k X = 1.0',k
c       Table X = 0.99997, Z = 0.00003

	k = k+1
	l = l-1
	open(lun,file='ags04.99997.00003.tron',status='old',err=9001)
	read(lun,'(a)',err=9002) line
	read(line,'(43x,f8.6,8x,f8.6)')xtab(k),ztab(l)
c       read blank line
	      read(lun,'(a)') line
c       read blank line
	      read(lun,'(a)') line

	print *,'k',k,xtab(k)
c-----------read data tlg,rlg and opacity-values
c	    read/write string 'logT' + rlg-values
	read(lun,   '(a6,19f7.3)') ct6,(rlg(i,k,l),i=1,nval)
	write(lun1              )     (rlg(i,k,l),i=1,nval)
c-----------read/write logT & opacity-values
	do j=1,nlin
	   read(lun,'(f5.3,2x,f6.3,18(1x,f6.3))',err=9002)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	   write(lun1,             err=9020)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	enddo
	close(lun)

c	print *,'k X = 1.0',k
c       Table X = 0.99999, Z = 0.00001

	k = k+1
	l = l-1
	open(lun,file='ags04.99999.00001.tron',status='old',err=9001)
	read(lun,'(a)',err=9002) line
	read(line,'(43x,f8.6,8x,f8.6)')xtab(k),ztab(l)
c       read blank line
	      read(lun,'(a)') line
c       read blank line
	      read(lun,'(a)') line

	print *,'k',k,xtab(k)
c-----------read data tlg,rlg and opacity-values
c	    read/write string 'logT' + rlg-values
	read(lun,   '(a6,19f7.3)') ct6,(rlg(i,k,l),i=1,nval)
	write(lun1              )     (rlg(i,k,l),i=1,nval)
c-----------read/write logT & opacity-values
	do j=1,nlin
	   read(lun,'(f5.3,2x,f6.3,18(1x,f6.3))',err=9002)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	   write(lun1,             err=9020)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	enddo
	close(lun)

c       Table X = 1.0, Z = 0.

	k = k+1
	l = l-1
	open(lun,file='ags04.one.tron',status='old',err=9001)
	read(lun,'(a)',err=9002) line
	read(line,'(43x,f8.6,8x,f8.6)')xtab(k),ztab(l)
c       read blank line
	      read(lun,'(a)') line
c       read blank line
	      read(lun,'(a)') line

	print *,'k',k,xtab(k)
c-----------read data tlg,rlg and opacity-values
c	    read/write string 'logT' + rlg-values
	read(lun,   '(a6,19f7.3)') ct6,(rlg(i,k,l),i=1,nval)
	write(lun1              )     (rlg(i,k,l),i=1,nval)
c-----------read/write logT & opacity-values
	do j=1,nlin
	   read(lun,'(f5.3,2x,f6.3,18(1x,f6.3))',err=9002)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	   write(lun1,             err=9020)
     &     tlg(j,k,l)           ,(opa(i,j,k,l),i=1,nval)
	enddo
	close(lun)


C c
C c
C c       output data in binary format
C c       --------------------------------
C c
C         write(lun1) nzva,(ztab(l),l=nzva,1,-1)
C         write(lun1) ntab,(xtab(k),k=1,ntab)
C         do 4001 l=nzva,1,-1
C            do 5001 k=1,ntab
C               print '(a,f4.2,a,f6.4)', 
C      +              ' write X= ',xtab(k),' Z= ',ztab(l)
C               write(lun1)ir(k,l),(rlg(i,k,l),i=1,ir(k,l))
C               write(lun1)nlin-4 ,(tlg(j,k,l),j=nlin,5,-1)!'-4' ignore 4 largest tlg
C               write(lun1)((opa(i,j,k,l),i=1,ir(k,l)),j=nlin,5,-1)
C  5001      continue
C  4001   continue
C c
        close(lun1)
	stop
c
 9001   print *,'ferg04bext: error in opening ',ifname(k,l),' lun= ',lun
	stop
 9011   print *,'ferg04bext: error in opening ',ofname,' lun= ',lun
	stop
 9002	print *,'ferg04bext: error in reading: line: ',j
	stop
 9020   print *,'ferg04bext: error in writing ',ofname,'lun= ',lun+1 
        stop
c
	end
