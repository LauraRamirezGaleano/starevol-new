
c***************************************************************************

	program calib_zsx255

c	calibration zsx pour le soleil avec une eventuelle variation de la
c	composition en X et Y dans la ZC externe au cours de l'evolution
c	en tenant compte de la perte de masse et de la rotation

c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice

c	CESAM4

c	le Z/X est fixe dans ce programme: zsx=0.0245 par exemple

c	cf.  A. Noel & N. Grevesse,   "Proceedings of the first SOHO
c	workshop, Anapolis, Maryland, USA,
c	25-28 August 1992 (ESA SP-348, november 1992)".

c	la calibration donne X et Y d'ou Z=1-X-Y qu'il faut ajuster a zsx

	implicit none
	
	integer ppara
	parameter(ppara=3)	!nombre total de parameteres a fixer

	
	integer un,i,j,long,ncouche,iconst,ivar,ivers

	double precision l,r,m,alfa,cal,zsx_sol,rsol,lsol,alpha,
     &    x_zc,y_zc,z_zc,zsx_zc,bid,zsx0,x0_c,y0_c,x0,z0_c,z0,y0,
     &   a(ppara,ppara),b(ppara),al(0:ppara),d(ppara),
     &   x(0:ppara),y(0:ppara),beta(0:ppara),rf(0:ppara),lf(0:ppara),
     &   betaf(0:ppara),xf(0:ppara),yf(0:ppara)
	double precision dl,dr,dbeta

	character*80 fname33,no,fname23
	character*2 sc
	double precision rsrsol_obj,lslsol_obj,zsx_obj
        data rsrsol_obj/1.0000d0/
        data lslsol_obj/1.0000d0/

C	les donnees

	cal=1.d-5		!niveau de calibration
	rsol = 6.9599d10
	lsol = 3.846d33


c.. Check that model is not calibrated already
	read(5,*)dl,dr,dbeta
	write(*,*)dl,dr,dbeta
	if ((abs(dl).le.1.d-5).and.(abs(dr).le.1.d-5)
     &     .and.(abs(dbeta).le.1.d-4)) then
	   print *,' Your model is calibrated : dr = ',dr,' dl =',dl,
	1	' d(Z/X) =',dbeta
	   stop
	endif

c.. Adapt reference Z/X value according to reference solar abundances considered
	read(5,*)sc
	if (trim(sc).eq.'1') then 
	   zsx_sol =  0.0245d0
	else if (trim(sc).eq.'2') then 
	   zsx_sol =  0.0165d0
	else if (trim(sc).eq.'3') then 
	   zsx_sol =  0.0180d0
	endif

	print*,zsx_sol
	
C	valeurs avant evolution pour le calcul des coeff. de calibration	
	
	fname33 = 'calib_donnees_diff'
	open (33,file=fname33,status='old')
 
        read(33,*) xf(0),yf(0),rf(0),lf(0),x(0),y(0),al(0)
        read(33,*) xf(1),yf(1),rf(1),lf(1),x(1),y(1),al(1)
        read(33,*) xf(2),yf(2),rf(2),lf(2),x(2),y(2),al(2)
        read(33,*) xf(3),yf(3),rf(3),lf(3),x(3),y(3),al(3)

		
	do i=0,3
	 beta(i)=(1.d0-x(i)-y(i))/x(i)		!beta =Z/X
	 betaf(i) = (1.d0-xf(i)-yf(i))/xf(i)
	enddo
	
	do i=1,ppara
	 do j=1,ppara
	  if(i .ne. j)then
	   a(i,j)=0.d0
	  else
	   a(i,i)=1.d0
	  endif
	 enddo
	enddo
	
	a(1,1)=(rf(1)-rf(0))/(al(1)-al(0))	!dR / d alpha
	a(1,2)=(rf(2)-rf(0))/(x(2)-x(0))	!dR / d X
	a(1,3)=(rf(3)-rf(0))/(beta(3)-beta(0))	!dR / d beta
c	a(1,3) = 0.d0
	
	a(2,1)=(lf(1)-lf(0))/(al(1)-al(0))	!dL / d alpha
	a(2,2)=(lf(2)-lf(0))/(x(2)-x(0))	!dL / d X
	a(2,3)=(lf(3)-lf(0))/(beta(3)-beta(0))	!dL / d beta
c	a(2,3) = 0.d0
	
	a(3,1)=(betaf(1)-betaf(0))/(al(1)-al(0))	!d beta* / d alpha 
	a(3,2)=(betaf(2)-betaf(0))/(x(2)-x(0))		!d beta* / d X
	a(3,3)=(betaf(3)-betaf(0))/(beta(3)-beta(0))	!d beta* / d beta
c	a(3,3) = 1.d0
	
	
***     Iteration 1
	r = rf(0)
	l = lf(0)
	x0 = x(0)
	z0 = 1.d0-x(0)-y(0)
	y0 = y(0)
	alfa = al(0)
	x_zc = xf(0)
	y_zc = yf(0)


      	z0 = 1.d0-x0-y0  
	z_zc=1.d0-x_zc-y_zc
	zsx_zc=z_zc/x_zc
	zsx0=z0/x0
	
	
	b(1)=r-rsrsol_obj
	b(2)=l-lslsol_obj
	b(3)=zsx_zc-zsx_sol	
	do i=1,3
	 d(i)=abs(b(i))
	enddo
	d(3)=d(3)/zsx_sol
	
c	print*
c	print*,'calibration'
c	write(6,2)(b(i),i=1,3)
c2	format(1x,'d_r=',1pd9.2,' d_l=',1pd9.2,1x,'d_Zsx=',1pd9.2,/)
	
	bid=0
	do i=1,3
	 bid=max(d(i),bid)
	enddo
	if(bid .lt. cal)then
	   print *,' Your model is calibrated : dr = ',dr,' dl =',dl,
     &  	' d(Z/X) =',dbeta
	   un=0
	else
	   call simq(a,b,3)
	endif
	 	 
	 alfa=alfa-b(1)
	 x0_c=x0-b(2)
	 zsx0=zsx0-b(3)
	 z0=x0_c*zsx0
	 y0_c=1.d0-x0_c-z0

	write(*,103) alfa,x0_c,y0_c,z0


 103    format('alpha =',0pf15.13,1x,'X0 =',0pf15.13,1x,'Y0 =',0pf15.13,
     &   1x,'Z0 =',0pf15.13)

	
c	 if(alfa .lt. 0.d0)then
c	  print*,'Erreur l/Hp negatif'
c	  stop
c	 elseif(x0_c .lt. 0.d0)then
c	  print*,'Erreur X < 0'
c	  stop
c	 elseif(x0_c .gt. 1.d0)then
c	  print*,'Erreur X > 1'
c	  stop
c	 elseif(y0_c .lt. 0.d0)then
c	  print*,'Erreur Y < 0'
c	  stop
c	 elseif(y0_c .gt. 1.d0)then
c	  print*,'Erreur Y > 1'
c	  stop
c	 elseif(z0 .lt. 0.d0)then
c	  print*,'Erreur Z < 0'
c	  stop
c	 elseif(z0 .gt. 1.d0)then
c	  print*,'Erreur Z > 1'
c	  stop
c	 endif

	 alpha=alfa
	 x0=x0_c
	 y0=y0_c

	close(33)
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
      tol=0.
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
