c
c  get the rates from vit.dat
c
       read(nr,*) (tgrid(i),i=1,ngrid)
c     locate the position of current T within the grid
c    and compute step, a, b, c and d parameters for cubic spline interpolation
       klo=1
       khi=ngrid
   20  if(khi-klo.gt.1) then
	k=(khi+klo)/2
	if(tgrid(k).gt.t8) then
	  khi=k
        else
	  klo=k
        endif
        goto 20
       endif
c  test the grid
       if(t8.lt.tgrid(1)) write(*,1005) t8,tgrid(1)
       if(t8.gt.tgrid(ngrid)) write(*,1006) t8,tgrid(ngrid)
c
       h=tgrid(khi)-tgrid(klo)
       a=(tgrid(khi) - t8)/h
       b=(t8 - tgrid(klo))/h
       c=a**3 - a
       d=b**3 - b
c
c  locate the position of current rho within the density grid (1.,3.,10.,30.)
c
       if (ane.le.3.) then
         anegridhigh=3.
         deltagrid=2.
	 jlo=1
         jhi=2
       else if (ane.le.10.) then
         anegridhigh=10.
         deltagrid=7.
	 jlo=2
	 jhi=3
       else
	 anegridhigh=30.
         deltagrid=20.
	 jlo=3
	 jhi=4
       endif
       do 29 j=1,ngrid
       vsecond(j)=0.
       do 29 k=1,4
   29  vsecondbeta(j,k)=0.
c
    1  continue
       read(nr,100,end=4) sym1,ka1,typem,rtp,sym2,ka2,typem2
  100  format(11x,a2,i3,a1,2x,a16,5x,a2,i3,a1)
       do 30 iz=0,103
       if (sym(iz).eq.sym1) goto 40
   30  continue
       if (iz.ge.103) then 
         write(*,*) 'Problem in identification of reaction'
         write(*,100) sym1,ka1,typem,rtp,sym2,ka2,typem2
	 stop
       endif
   40  continue
       kz=iz
       i=ka1+iminz(kz)-iaminz(kz)
       i=i+i/isomz(kz)
       if (typem.eq.'m'.and.i+1.eq.isomz(kz)) i=i+1
c
       read(nr,105) test
  105  format(a20)
       if (test.eq.'                    ') goto 1
       read(nr,110) flagi
  110  format(1x,f4.0)
c           read the reaction kind:
c                0: beta from exp or GT2    -> keep as it is   and lin interpol.
c               -1: (n,g) from beer         -> * RHO * S       and lin interpol.
c              -12: (n,g) from SMOKER       -> * RHO           and lin interpol.
c        1,3,10,30: (beta rho sensitive)    -> interpolate rho and log interpol.
c              -10: logarithmic interpol.   -> * RHO           and log interpol.
c 
       if(flagi.le.0.) then
            read(nr,*) (vdum(j), j = 1,ngrid)
c           read(nr,*) (vsecond(j), j = 1,ngrid)
       else 
	    do k=1,4
	    read(nr,*) (vdumbeta(j,k), j=1,ngrid)
c	    read(nr,*) (vsecondbeta(j,k), j=1,ngrid)
            if(k.ne.4) read(nr,*) dummy 
            enddo
       endif
c
c if residual nucleus is an isomer: get the index m of the isomer 
       if (typem.eq.'m'.and.i.ne.isomz(kz)) goto 1
       if (typem.eq.'g'.and.i+1.ne.isomz(kz)) goto 1
       if (typem.eq.' '.and.i.eq.isomz(kz)) goto 1
       if (typem.eq.' '.and.i+1.eq.isomz(kz)) goto 1
       iresidual=0
       do 70 m=1,niso
        if (sym(kzm(m)).eq.sym2.and.kam(m).eq.ka2) then
	   iresidual=1
	   goto 80
	endif
   70  continue
   80  continue
       if ((typem2.eq.'m'.or.typem2.eq.'g').and.iresidual.eq.0) goto 1
       if (typem2.eq.' '.and.iresidual.eq.1) goto 1
c  Isomer not found as residual nucleus or residual nucleus is isomer -> skip!'
c
       if (kz.lt.kzmin.or.kz.gt.kzmax+1) goto 1
       if (i.lt.iminz(kz).or.i.gt.iminz(kz+1)-1) goto 1
c  compute cubic spline interpolation
c
       if(flagi.le.0.) rate= a*vdum(klo)+b*vdum(khi)+
     &                   (c*vsecond(klo)+d*vsecond(khi))*(h*h)/6. 
c
       if (flagi.eq.-1.) rate=rate*S
       if (flagi.eq.-10.) rate=10.**rate*rho
       if (flagi.eq.-11.) rate=10.**rate
       if (flagi.eq.-12.) rate=rate*rho
       if (typem.eq.'m'.or.typem.eq.'g') rcorr(i)=1.
       if (rtp.eq.symng) then
          if (flagi.gt.-11.5) iexpng(i)=1
      	  if (typem2.ne.'m') then
	    signg(i)=rate/rho/avog*rcorr(i)
          else
	    signgm(m)=rate/rho/avog*rcorr(i)
          endif
       endif
       if (rtp.eq.symna) signa(i)=rate/rho/avog
       if (rtp.eq.symga) sigga(i)=rate
       if (rtp.eq.symb) then
            if (flagi.le.0.) goto 50
c  reaction has to be interpolated in density
c   spline interpolation in temperature for two density grid points
c   (note: log of beta decay rate is handled)    
            vbetalow=a*vdumbeta(klo,jlo)    
     &                  +b*vdumbeta(khi,jlo) 
     &                  +(c*vsecondbeta(klo,jlo) 
     &                  +d*vsecondbeta(khi,jlo))*(h*h)/6. 
            vbetahigh=a*vdumbeta(klo,jhi)    
     &                  +b*vdumbeta(khi,jhi) 
     &                  +(c*vsecondbeta(klo,jhi) 
     &                  + d*vsecondbeta(khi,jhi))*(h*h)/6. 
            vbetalow=10.**vbetalow
	    vbetahigh=10.**vbetahigh
c  linear interpolation in density 
            rate=(vbetahigh-vbetalow)/deltagrid*(ane-anegridhigh)
     &           +vbetahigh
c  extrapolation often leads to negative rates 
c  check if so, and then put it to the closest density grid point
c  (supplementary table values would be needed for those cases)
            if(rate.lt.0.) then
               if(ane-anegridhigh.lt.0.) then
c  extrapolation at lower end of table
                  rate=vbetalow
               else
c  extrapolation at upper end of table
                  rate=vbetahigh
               endif    
            endif
  50        continue
            if (sym2.eq.sym(kz+1)) then
		if (typem2.ne.'m') then
c!
c     	  if (kz.eq.55.and.ka1.eq.134) rate=rate/3.
		  blamda(i)=rate
		else 
		  blamdam(m)=rate
		endif
	    endif
            if (sym2.eq.sym(kz-1)) then
		if (typem2.ne.'m') then
		  blamdap(i)=rate
		else 
		  blamdapm(m)=rate
		endif
	    endif
            if (sym2.ne.sym(kz+1).and.sym2.ne.sym(kz-1)) then
             write(*,*) 'Problem in identification of beta-decay'
             write(*,100) sym1,ka1,rtp,sym2,ka2
	     stop
            endif
       endif
