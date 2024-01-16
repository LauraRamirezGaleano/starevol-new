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
c Le 19 juin 2006: calcul avec Z = 1.739765022D-02
c================================================================
c
        implicit double precision(a-h,o-z)

        dimension xlum(3),ray(3),yinit(3),alpha(3)
        dimension dxlum(3),dray(3)

        character*40 fname33
c
c----------------------------------------------------------------
c
c        fname33 = 'calib_donnees_std_gag291'
c        fname33 = 'calib_donnees_stdCE_shell_gag291'
c        fname33 = 'calib_donnees_stdCE_H2HmHe_gag291'
        fname33 = 'calib_donnees'
        open (33,file=fname33,status='old')
 
        read(33,*) xlum(1),ray(1),yinit(1),alpha(1)
        read(33,*) xlum(2),ray(2),yinit(2),alpha(2)
        read(33,*) xlum(3),ray(3),yinit(3),alpha(3)
c

	read(5,*)dl,dr
	write(*,*)dl,dr
	if (abs(dl).le.1.d-5.and.abs(dr).le.1.d-5) then
	   print *,' Your model is calibrated : dr = ',dr,' dl =',dl
	   stop
	endif



c	raysol = 6.9599d10
c	xlumsol = 3.846d33

c         dxlum(1)=xlum(1)/xlumsol
c         dxlum(2)=xlum(2)/xlumsol
c         dxlum(3)=xlum(3)/xlumsol
c
c         dray(1)=ray(1)/raysol
c         dray(2)=ray(2)/raysol
c         dray(3)=ray(3)/raysol

         dxlum(1)=xlum(1)
         dxlum(2)=xlum(2)
         dxlum(3)=xlum(3)
c
         dray(1)=ray(1)
         dray(2)=ray(2)
         dray(3)=ray(3)
c
c
         dy2=yinit(1)-yinit(2)
         dy3=yinit(1)-yinit(3)
         dal2=alpha(1)-alpha(2)
         dal3=alpha(1)-alpha(3)
c
c calcul des coefficients a et b du systeme d'equations
c------------------------------------------------------
c
        if (yinit(1).ne.yinit(2)) then
         a2=dxlum(1)-dxlum(2)
         a2=a2*dy3/dy2
         a2=dxlum(1)-dxlum(3)-a2
         denom=dal2*dy3/dy2
         denom=dal3-denom
         a2=a2/denom
c
         a1=(dxlum(1)-dxlum(2)-a2*dal2)/dy2
c
         a0=dxlum(1)-a1*yinit(1)-a2*alpha(1)
c
         b2=dray(1)-dray(3)-((dray(1)-dray(2))*dy3)/dy2
         b2=b2/(dal3-(dal2*dy3)/dy2)
c
         b1=(dray(1)-dray(2)-b2*dal2)/dy2
c
         b0=dray(1)-b1*yinit(1)-b2*alpha(1)
        else
         a2=dxlum(1)-dxlum(3)
         a2=a2*dy2/dy3
         a2=dxlum(1)-dxlum(2)-a2
         denom=dal3*dy2/dy3
         denom=dal2-denom
         a2=a2/denom
c
         a1=(dxlum(1)-dxlum(3)-a2*dal3)/dy3
c
         a0=dxlum(1)-a1*yinit(1)-a2*alpha(1)
c
         b2=dray(1)-dray(2)-((dray(1)-dray(3))*dy2)/dy3
         b2=b2/(dal2-(dal3*dy2)/dy3)
c
         b1=(dray(1)-dray(3)-b2*dal3)/dy3
c
         b0=dray(1)-b1*yinit(1)-b2*alpha(1)
c
        endif
c
c calcul de alpha et Yinit pour le modele calibre
c------------------------------------------------
c
         xalpha=1.-b0-(1.-a0)*b1/a1
         xalpha=xalpha/(b2-b1*a2/a1)
c
         iialpha=int(xalpha*1e+6)
         xx=xalpha*1e+6-iialpha
         xxx=dble(iialpha)
         if (xx.lt.0.5) then
           xxalpha=1e-6*xxx 
         else
           xxalpha=1e-6*(xxx+1.)
         endif
c
         xyinit=(1.-a0-a2*xalpha)/a1
c
         iiyinit=int(xyinit*1e+6)
         xx=xyinit*1e+6-iiyinit
         xxx=dble(iiyinit)
         if (xx.lt.0.5) then
           xxyinit=1e-6*xxx
         else
           xxyinit=1e-6*(xxx+1.)
         endif
         write(*,*)'New Y0 =', xxyinit,'New alpha =',xxalpha
      close(33)
       end

