c==============================================================
c
	program calibration
c
c==============================================================
c
c Calcul Alpha et Yinit du modele le plus proche du soleil. Pour
c cela il faut avoir entre les valeurs de L,R,Y,alpha des trois 
c modeles deja calcules dans le fichier calib_donnees (L en erg/s,
c R en cm).
c La calibration utilise les valeurs du rayon et de la luminosite
c du soleil de Guenther (92) (L=0.38515*10^34 erg.s^-1 et
c R=0.69598*10^11 cm).
c Alpha, Yinit, sont alors ecrits dans le fichier modele_calibre.
c
c
c       Le 19 juin 2006: calcul avec Z = 1.739765022D-02
c	Le 11 novembre 2017 : calibration avec AGSS09
c================================================================
c
        implicit none


	integer i,j
	integer indx,index


	double precision xlum,ray,xinit,yinit,alpha,xsurf,ysurf,beta,
	1    betaf,dxlum,dray,d_Zsx,A,B,absB,dl,dr,dbeta,d_l,d_r,bid,
	2    cal,xinit_cal,yinit_cal,beta_cal,alpha_cal,dzsx,
	3    x0,z0,y0,r0,l0,zsx0,alfa0,x_zc,y_zc,z_zc,zsx_zc
	double precision xlumsol,raysol,zsurxsol
	double precision d

        dimension xlum(4),ray(4),xinit(4),yinit(4),alpha(4)
	dimension xsurf(4),ysurf(4),beta(4),betaf(4)
        dimension dxlum(4),dray(4),dzsx(4)
	dimension A(3,3),B(3),absB(3),indx(3),index(3)

C	Solar data 
c	parameter (cal = 1.d-5,xlumsol = 3.846d33, raysol = 6.9599d10 )
	parameter (cal = 1.d-5,xlumsol = 3.828d33, raysol = 6.9599d10 )

        character*40 fname33
	character*2 sc
c
c----------------------------------------------------------------

        fname33 = 'calib_donnees_diff'
        open (33,file=fname33,status='old')
 
        read(33,*) xsurf(1),ysurf(1),xlum(1),ray(1),xinit(1),yinit(1)
	1    ,alpha(1)
        read(33,*) xsurf(2),ysurf(2),xlum(2),ray(2),xinit(2),yinit(2)
	2    ,alpha(2)
        read(33,*) xsurf(3),ysurf(3),xlum(3),ray(3),xinit(3),yinit(3)
	3    ,alpha(3)
        read(33,*) xsurf(4),ysurf(4),xlum(4),ray(4),xinit(4),yinit(4)
	4    ,alpha(4)
c
	do i = 1,4
	   beta(i)=(1.d0-xinit(i)-yinit(i))/xinit(i) !beta =Z/X
	   betaf(i)=(1.d0-xsurf(i)-ysurf(i))/xsurf(i) !beta =Z/X
	enddo

	read(5,*)dl,dr,dbeta
	write(*,*)dl,dr,dbeta
	if ((abs(dl).le.1.d-5).and.(abs(dr).le.1.d-5)
     &     .and.(abs(dbeta).le.1.d-4)) then
	   print *,' Your model is calibrated : dr = ',dr,' dl =',dl,
     &	   ' d(Z/X) =',dbeta
	   stop
	endif

	read(5,*)sc
	if (trim(sc).eq.'1') then 
	   zsurxsol =  0.0245d0
	else if (trim(sc).eq.'2') then 
	   zsurxsol =  0.0165d0
	else if (trim(sc).eq.'3') then 
	   zsurxsol =  0.018207d0
	endif

	print*,zsurxsol

C       Matrix to be reverted

	do i = 1,3
	   do j = 1,3
	      if (i.ne.j) then
		 A(i,j)=0.d0
	      else
		 A(i,i)=1.d0
	      endif
	   enddo
	enddo

	A(1,1) = (ray(2) - ray(1))/(alpha(2)-alpha(1))  !dR / d alpha
        A(1,2) = (ray(3) - ray(1))/(xinit(3)-xinit(1))  !dR / d X
c        A(1,3) = (ray(4) - ray(1))/(beta(4)-beta(1))   !dR / d beta
        A(1,3) = 0.d0
	
        A(2,1) = (xlum(2) - xlum(1))/(alpha(2)-alpha(1)) !dL / d alpha
	A(2,2) = (xlum(3) - xlum(1))/(xinit(3)-xinit(1)) !dL / d X
c        A(2,3) = (xlum(4) - xlum(1))/(beta(4)-beta(1))  !dL / d beta
        A(2,3) = 0.d0


        A(3,1) = (betaf(2) - betaf(1))/(alpha(2)-alpha(1)) !d beta* / d alpha 
        A(3,2) = (betaf(3) - betaf(1))/(xinit(3)-xinit(1)) !d beta* / d X
c        A(3,3) = (betaf(4) - betaf(1))/(beta(4)-beta(1))  !d beta* / d beta
        A(3,3) = 1.d0

***     Iteration 1
	r0 = ray(1)
	l0 = xlum(1)
	x0 = xinit(1)
	z0 = 1.d0-xinit(1)-yinit(1)
	y0 = yinit(1)
	alfa0 = alpha(1)
	x_zc = xsurf(1)
	y_zc = ysurf(1)


      	z0 = 1.d0-x0-y0  
	z_zc = 1.d0-x_zc-y_zc
	zsx_zc = z_zc/x_zc
	zsx0 = z0/x0

        write(6,10)x0,1.d0-x0-z0,z0,zsx0,alfa0
10	format(1x,'modele initial:',/,
     1	t2,'X0=',1pd10.3,' Y0=',1pd10.3,' Z0=',1pd10.3,
     2	' Z0/X0=',1pd10.3,' alpha=',1pd10.3,/,1x,/)
	
	write(6,11)r0,l0,zsx_zc,x_zc,y_zc,z_zc
11	format(1x,'modele evolue:',/,
     1	t2,'R=',1pd12.5,' L=',1pd12.5,' Z/X=',1pd12.5,
     3	t2,'Xzc=',1pd12.5,' Yzc=',1pd12.5,' Zzc=',1pd12.5,//)
	
	write(6,12)raysol,xlumsol,zsurxsol
12	format(t2,'Rsol, Lsol, Z/Xsol a obtenir pour calibrer',
     1	/,t2,'R sol=',1pd12.5,' L sol=',1pd12.5,' Z/X sol=',1pd12.5,/)

	

C       Result matrix
	
        B(1) = ray(1) - 1.d0
        B(2) = xlum(1) - 1.d0
        B(3) = betaf(1) - zsurxsol


        do i = 1,3
           absB(i) = abs(B(i))
        enddo
        absB(3) = absB(3)/zsurxsol

c$$$c..   Resolution du systeme pour trouver les corrections a appliquer
c$$$c..   en utilisant la methode LU (routines Numerical Recipes)
c$$$
c$$$        call ludcmp(A,3,3,indx,d)
c$$$    call lubksb(A,3,3,index,absB)

c       MODIF TD - Modele AP
	
C	Inversion fo tghe matrix to retrieve delta alpha, delta X and delta beta
	bid = 0
	do i = 1,3
	   bid = max(absB(i),bid)
	enddo
	if(bid .lt. cal)then
	   write(6,*)'modele calibre'
c	   un = 0
	else
	   call simq(A,B,3)
	endif
	
c	alfa = alfa-b(1)
c	x0_c = x0-b(2)
c	zsx0 = zsx0-b(3)
c	z0 = x0*zsx0
c	y0_c = 1.d0-x0_c-z0


c FIN MODIF TD	
	

c..   New values of alpha, X and Y to be used.

        alpha_cal = alpha(1)-absB(1)
        xinit_cal = xinit(1)-absB(2)
        beta_cal = beta(1)-absB(3)
        
        yinit_cal = 1.d0-xinit_cal*(1.d0+beta_cal)
	
c..     Error on the solution
	d_l = B(1)
	d_r = B(2)
	d_Zsx = B(3)
	
c	write(*,102)d_r,d_l,d_Zsx,alpha_cal,xinit_cal,yinit_cal,
c     &   beta_cal*xinit_cal
	write(*,103)alpha_cal,xinit_cal,yinit_cal,beta_cal*xinit_cal


c 102	format('d_r = ',1pd8.1,1x,'d_l = ',1pd8.1,1x,'d_Zsx = ',1pd8.1,
c     &  1x,'alpha = ',1pf15.13,1x,'X0 =',0pf15.13,1x,'Y0 =',0pf15.13,1x,
c     &  'Z0 =',0pf15.13)

 103    format('alpha =',0pf15.13,1x,'X0 =',0pf15.13,1x,'Y0 =',0pf15.13,
     &   1x,'Z0 =',0pf15.13)

	close(33)

C       If matrix inversion resulted in a non-physical value of the paralmeters, exit
	 if(alpha_cal .lt. 0.d0)then
	  print*,'Erreur l/Hp negatif'
	  stop
	 elseif(xinit_cal .lt. 0.d0)then
	  print*,'Erreur X < 0'
	  stop
	 elseif(xinit_cal .gt. 1.d0)then
	  print*,'Erreur X > 1'
	  stop
	 elseif(yinit_cal .lt. 0.d0)then
	  print*,'Erreur Y < 0'
	  stop
	 elseif(yinit_cal .gt. 1.d0)then
	  print*,'Erreur Y > 1'
	  stop
	 elseif(beta_cal*xinit_cal .lt. 0.d0)then
	  print*,'Erreur Z < 0'
	  stop
	 elseif(beta_cal*xinit_cal .gt. 1.d0)then
	  print*,'Erreur Z > 1'
	  stop
	 endif

	 alpha=alfa0
	 x0=xinit_cal
	 y0=yinit_cal
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

*----------------------------------------------------------

      SUBROUTINE LUBKSB (A,N,NP,INDX,B)

*----------------------------------------------------------

      implicit none

      INTEGER N,NP,INDX
      INTEGER I,II,J,LL

      DOUBLE PRECISION A,B,SUM

      DIMENSION A(NP,NP),INDX(N),B(N)

      II = 0
      DO I = 1,N
         LL = INDX(I)
         SUM = B(LL)
         B(LL) = B(I)
         IF (II.NE.0) THEN
            DO J = II,I-1
               SUM = SUM-A(I,J)*B(J)
            ENDDO
         ELSE IF (SUM.NE.0.D0) THEN
            II = I
         ENDIF
         B(I) = SUM
      ENDDO
      DO I = N,1,-1
         SUM = B(I)
         IF (I.LT.N) THEN
            DO J = I+1,N
               SUM = SUM-A(I,J)*B(J)
            ENDDO
         ENDIF
         B(I) = SUM/A(I,I)
      ENDDO

      RETURN
      END

*----------------------------------------------------------

      SUBROUTINE LUDCMP (A,N,NP,INDX,D)

*----------------------------------------------------------

      implicit none

      INTEGER N,NP,INDX
      INTEGER NMAX,IMAX
      INTEGER I,J,K

      DOUBLE PRECISION A,D,VV,AAMAX,SUM,DUM,TINY

      PARAMETER (NMAX = 100,TINY = 1.D-20)

      DIMENSION A(NP,NP),INDX(N),VV(NMAX)

      D = 1.D0
      DO I = 1,N
         AAMAX = 0.D0
         DO J = 1,N
            IF (ABS(A(I,J)).GT.AAMAX) AAMAX = ABS(A(I,J))
         ENDDO
         IF (AAMAX.EQ.0.D0) STOP 'Singular matrix.'
         VV(I) = 1.D0/AAMAX
      ENDDO
      DO J = 1,N
         IF (J.GT.1) THEN
            DO I = 1,J-1
               SUM = A(I,J)
               IF (I.GT.1) THEN
                  DO K = 1,I-1
                     SUM = SUM-A(I,K)*A(K,J)
                  ENDDO
                  A(I,J) = SUM
               ENDIF
            ENDDO
         ENDIF
         AAMAX = 0.D0
         DO I = J,N
            SUM = A(I,J)
            IF (J.GT.1) THEN
               DO K = 1,J-1
                  SUM = SUM-A(I,K)*A(K,J)
               ENDDO
               A(I,J) = SUM
            ENDIF
            DUM = VV(I)*ABS(SUM)
            IF (DUM.GE.AAMAX) THEN
               IMAX = I
               AAMAX = DUM
            ENDIF
         ENDDO
         IF (J.NE.IMAX) THEN
            DO K = 1,N
               DUM = A(IMAX,K)
               A(IMAX,K) = A(J,K)
               A(J,K) = DUM
            ENDDO
            D = -D
            VV(IMAX) = VV(J)
         ENDIF
         INDX(J) = IMAX
         IF (J.NE.N) THEN
            IF (A(J,J).EQ.0.D0) A(J,J) = TINY
            DUM = 1.D0/A(J,J)
            DO I = J+1,N
               A(I,J) = A(I,J)*DUM
            ENDDO
         ENDIF
      ENDDO
      IF (A(N,N).EQ.0.D0) A(N,N) = TINY

      RETURN
      END


