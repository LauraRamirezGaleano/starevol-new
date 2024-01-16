
c***************************************************************************

	program calibration



C       Program to calibrate the solar model at 4.57 Gyrs.
C       Based on a routine developed by  P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
C	for the CESAM code version 4
	
C       The quantities that vary are :
C       al = alpha MLT parameter
C       x = initial hydrogen mass fraction
C       y = initial helium mass fraction
C       zsx = initial metals mass fraction over hydrogen mass fraction

C       Input parameters
C       =================
C
C       
C       Calibration matrix containing the values of alpha, X and Y for 4 models:
C       1st model : reference model (al(1), x(1), y(1))
C       2d model : change alpha (al(2), x(2) = x(1), y(2) = y(1))
C       3d model : change X but maintain the ration Z/X of the reference model (al(3) = al(1), x(3), y(3))
C       4th model : change Y without maintaining the Z/X of the reference model (al(4) = al(1), x(4) = x(1), y(4))
C
C       solar age, luminosity and radius
C       solar Z/X ratio considered (GN93 or AGS06)
C       wanted calibration error
C
C       File containing the ouputs for the models of the calibration matrix @ the considered solar age
C       final luminosity, final radius, final hydrogen mass fraction, final Z/X ratio 
C
C       File containing the characteristics of the currently computed model
C       final radius, final luminosity, initial X, initial Y, alpha, final X, final Y
C
C       Ouputs
C       ======
C       New values for X,Y,Z Z/X and alpha to be adopted for the next calibration step
C
C       Version November 2017 for automated script
C
C ===============================================================================================================

	implicit none
	
	integer ppara
	parameter(ppara=3)	!Total number of parameters to be set

	
	integer un,i,j,long,ncouche,iconst,ivar,ivers

	double precision l,r,m,alfa,cal,zsx_sol,rsol,lsol,alpha,
     &    x_zc,y_zc,z_zc,zsx_zc,bid,zsx0,x0_c,y0_c,z0_c,zsx0_c,x0,z0,y0,
     &   a(ppara,ppara),b(ppara),al(0:ppara),d(ppara),
     &   x(0:ppara),y(0:ppara),beta(0:ppara),rf(0:ppara),lf(0:ppara),
     &   betaf(0:ppara),xf(0:ppara),yf(0:ppara),dl,dr,dbeta

	character*80 fname33,no,fname23
	character*2 sc
	
	double precision rsrsol_obj,lslsol_obj,zsx_obj
	double precision xlumsol,raysol,zsurxsol
        data rsrsol_obj/1.0000d0/
        data lslsol_obj/1.0000d0/
C	Solar data 
	parameter (cal = 1.d-5,xlumsol = 3.828d33, raysol = 6.9599d10 )



C	Checking whether the model is already calibrated 
	read(5,*)dl,dr,dbeta
	write(*,*)dl,dr,dbeta
	if ((abs(dl).le.1.d-5).and.(abs(dr).le.1.d-5)
     &      .and.(abs(dbeta).le.1.d-4)) then
	   print *,' Your model is calibrated : dr = ',dr,' dl =',dl,
     &	   ' d(Z/X) =',dbeta
	   stop
	endif
	
C       Checking the reference for the solar chemical composition
	read(5,*)sc
	if (trim(sc).eq.'1') then 
	   zsurxsol =  0.0245d0
	else if (trim(sc).eq.'2') then 
	   zsurxsol =  0.0165d0
	else if (trim(sc).eq.'3') then 
	   zsurxsol =  0.018207d0
	endif
	
	print*,zsurxsol
	
	
C       Reading the output of 4 solar models to be used to build the calibration matrix
	
	fname33 = 'calib_donnees_diff'
	open (33,file=fname33,status='old')
	
        read(33,*) xf(0),yf(0),lf(0),rf(0),x(0),y(0),al(0)
        read(33,*) xf(1),yf(1),lf(1),rf(1),x(1),y(1),al(1)
        read(33,*) xf(2),yf(2),lf(2),rf(2),x(2),y(2),al(2)
        read(33,*) xf(3),yf(3),lf(3),rf(3),x(3),y(3),al(3)

	print *, 'valeur entree modele 0 calib_diff'
	print *, xf(0),yf(0),lf(0),rf(0),x(0),y(0),al(0)
	do i=0,3
	   beta(i)=(1.d0-x(i)-y(i))/x(i) !beta =Z/X
	   betaf(i) = (1.d0-xf(i)-yf(i))/xf(i)
	enddo
	

C       Matrix to be reverted
	
	do i=1,ppara
	   do j=1,ppara
	      if(i.ne.j)then
		 a(i,j)=0.d0
	      else
		 a(i,i)=1.d0
	      endif
	   enddo
	enddo
	
	a(1,1)=(rf(1)-rf(0))/(al(1)-al(0)) !dR / d alpha
	a(1,2)=(rf(2)-rf(0))/(x(2)-x(0)) !dR / d X
C	a(1,3)=(rf(3)-rf(0))/(beta(3)-beta(0))	!dR / d beta
	a(1,3) = 0.d0
	
	a(2,1)=(lf(1)-lf(0))/(al(1)-al(0)) !dL / d alpha
	a(2,2)=(lf(2)-lf(0))/(x(2)-x(0)) !dL / d X
C	a(2,3)=(lf(3)-lf(0))/(beta(3)-beta(0))	!dL / d beta
	a(2,3) = 0.d0
	
	a(3,1)=(betaf(1)-betaf(0))/(al(1)-al(0)) !d beta* / d alpha 
	a(3,2)=(betaf(2)-betaf(0))/(x(2)-x(0)) !d beta* / d X
C	a(3,3)=(betaf(3)-betaf(0))/(beta(3)-beta(0))	!d beta* / d beta
	a(3,3) = 1.d0
	
	
***     Iteration 1
	r = rf(0)
	l = lf(0)
	x0 = x(0)
	z0 = 1.d0-x(0)-y(0)
	y0 = y(0)
	alfa = al(0)
	x_zc = xf(0)
	y_zc = yf(0)


c      	z0 = 1.d0-x0-y0  
	z_zc = 1.d0-x_zc-y_zc
	zsx_zc = z_zc/x_zc
	zsx0 = z0/x0
	
	write(6,10)x0,1.d0-x0-z0,z0,zsx0,alfa
10	format(1x,'modele initial:',/,
     1	t2,'Xini=',1pd10.3,' Yini=',1pd10.3,' Zini=',1pd10.3,
     2	' Zini/Xini=',1pd10.3,' altini=',1pd10.3,/,1x,/)
	
	write(6,11)r,l,zsx_zc,x_zc,y_zc,z_zc
11	format(1x,'modele evolue:',/,
     1	t2,'R=',1pd12.5,' L=',1pd12.5,' Z/X=',1pd12.5,
     3	t2,'Xzc=',1pd12.5,' Yzc=',1pd12.5,' Zzc=',1pd12.5,//)
	
	write(6,12)raysol,xlumsol,zsurxsol
12	format(t2,'Rsol, Lsol, Z/Xsol a obtenir pour calibrer',
     1	/,t2,'R sol=',1pd12.5,' L sol=',1pd12.5,' Z/X sol=',1pd12.5,/)

C       Right hand side matrix
	
	b(1) = r-rsrsol_obj
	b(2) = l-lslsol_obj
c	b(3) = zsx_zc-zsx_sol
	b(3) = zsx_zc-zsurxsol
	do i=1,3
	   d(i)=abs(b(i))
	enddo
c	d(3)=d(3)/zsx_sol
	d(3)=d(3)/zsurxsol
	
	print*
	print*,'calibration'
	write(6,2)(b(i),i=1,3)
 2	format(1x,'d_r=',1pd9.2,' d_l=',1pd9.2,1x,'d_Zsx=',1pd9.2,/)
	
C	Inversion fo tghe matrix to retrieve delta alpha, delta X and delta beta
	bid = 0
	do i = 1,3
	   bid = max(d(i),bid)
	enddo
	if(bid.lt.cal)then
	   write(6,*)'modele calibre'
	   un = 0
	else
	   call simq(a,d,3)
	endif
	
c$$$	alfa = alfa-b(1)
c$$$	x0_c = x0-b(2)
c$$$	zsx0 = zsx0-b(3)
c$$$	z0 = x0*zsx0
c$$$	y0_c = 1.d0-x0_c-z0

c$$$	alfa = alfa-abs(b(1))
c$$$	x0_c = x0-abs(b(2))
c$$$	zsx0 = zsx0-abs(b(3))
c$$$	zsx0_c = zsx0-abs(b(3))
c$$$	z0_c = x0_c*zsx0_c
c$$$	y0_c = 1.d0-x0_c-z0_c

	alfa = alfa-d(1)
	x0_c = x0-d(2)
	zsx0_c = zsx0-d(3)
	z0_c = x0_c*zsx0_c
	y0_c = 1.d0-x0_c-z0_c

	print *, 'alfa',alfa,'x0_c',x0_c,'y0_c',y0_c,'z0_c',z0
	print *, 'd',d(1), d(2), d(3)
	print *, 'b', b(1), b(2), b(3)
	
C       Print new values in the output_calibration_i file
	 
	 print*
c	 write(6,103)x0_c,y0_c,z0,zsx0,alfa
c 3	 format(t2,'nouvelles valeurs pour calibration',/,
c     1	 t2,'X0='0pf15.13,' Y0=',0pf15.13,' Z0=',0pf15.13,
c     2	' Z0/X0=',1pd16.9,' alpha=',0pf15.13,/,t2,/)

         write (6,103)alfa,x0_c,y0_c,z0_c
 103	 format('alpha =',0pf15.13,1x,'X0 =',0pf15.13,1x,'Y0 =',0pf15.13,
     &   1x,'Z0 =',0pf15.13)

         close(33)
	 
C       If matrix inversion resulted in a non-physical value of the paralmeters, exit
	 if(alfa .lt. 0.d0)then
	  print*,'Erreur l/Hp negatif'
	  stop
	 elseif(x0_c .lt. 0.d0)then
	  print*,'Erreur X < 0'
	  stop
	 elseif(x0_c .gt. 1.d0)then
	  print*,'Erreur X > 1'
	  stop
	 elseif(y0_c .lt. 0.d0)then
	  print*,'Erreur Y < 0'
	  stop
	 elseif(y0_c .gt. 1.d0)then
	  print*,'Erreur Y > 1'
	  stop
	 elseif(z0_c .lt. 0.d0)then
	  print*,'Erreur Z < 0'
	  stop
	 elseif(z0_c .gt. 1.d0)then
	  print*,'Erreur Z > 1'
	  stop
	 endif

	 alpha=alfa
	 x0=x0_c
	 y0=y0_c
	stop

	end

c*************************************************************

	subroutine simq(a,b,n)

c	resolution d'un systeme lineaire -ssp 360

c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91

	implicit none

	integer*4 n
	integer*4 jj,j,jy,it,i,ij,imax,i1,i2,iqs,ix,ixj,jx,ixjx,
     1jjx,ny,ia,ib,ic,k

	double precision a(1),b(1)
	double precision tol,biga,save

c        forward solution
	tol=-1
	jj=-n
	do 65 j=1,n
      jy=j+1
      jj=jj+n+1
      biga=0.
      it=jj-j
      do 30 i=j,n
c        search for maximum coefficient in column
      ij=it+i
      if(abs(biga)-abs(a(ij)))20,30,30
   20 biga=a(ij)
      imax=i
   30 continue
c        test for pivot less than tolerance (singular matrix)
      if(abs(biga)-tol) 35,35,40
 35   print 36,biga,tol
 36   format(2x,'dans simq biga=',1pe10.3,' < tol=',1pe10.3)
      stop
c        interchange rows if necessary
   40 i1=j+n*(j-2)
      it=imax-j
      do 50 k=j,n
      i1=i1+n
      i2=i1+it
      save=a(i1)
      a(i1)=a(i2)
      a(i2)=save
c        divide equation by leading coefficient
   50 a(i1)=a(i1)/biga
      save=b(imax)
      b(imax)=b(j)
      b(j)=save/biga
c        eliminate next variable
      if(j-n) 55,70,55
   55 iqs=n*(j-1)
      do 65 ix=jy,n
      ixj=iqs+ix
      it=j-ix
      do 60 jx=jy,n
      ixjx=n*(jx-1)+ix
      jjx=ixjx+it
   60 a(ixjx)=a(ixjx)-(a(ixj)*a(jjx))
   65 b(ix)=b(ix)-(b(j)*a(ixj))
c        back solution
   70 ny=n-1
      it=n*n
      do 80 j=1,ny
      ia=it-j
      ib=n-j
      ic=n
      do 80 k=1,j
      b(ib)=b(ib)-a(ia)*b(ic)
      ia=ia-n
   80 ic=ic-1
      return
      end
