
      PROGRAM STAREVOLMAN

************************************************************************
* Model file manager, for starevolc.f (Version 2.80 Dec 2005       (LS)*
* Module allowing to transform solar abundances in log(X/H)+12         *
* scale into solar mass fractions ( solx, Dec 2005)                (AP)*
*                                                                      *
* $LastChangedDate:: 2012-06-26 14:05:25 +0200 (Tue, 26 Jun 2012)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 144                                                         $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.data'
      include 'evolcom.mod'

      integer izchange,iheavy,refsolchem
      integer idiffvr

      character*1 opt
      character*1 rotationtype

      double precision diffvr,system
      double precision xheavysol,xFesol,fracFe,FeHsol
      double precision znew,mini,rini,pindex,FeH!,sigc
      double precision xini,yini,zini,xsol,ysol,zsolar,zcalc
      double precision breaktime
      double precision xhe4

      common /param/ znew,mini,rini,zini,pindex,FeH
      common /bigbang/xini,yini,xsol,ysol,zsolar,xFesol,
     &     FeHsol,fracFe,iheavy
      common /rotlaw/ diffvr,breaktime,idiffvr

      call getenv('SOLCOMP',refsolar)
      refsolar = trim(refsolar)
      if (refsolar.eq.'GN93') then
         refsolchem = 1
         xFesol = 1.147417671020920d-3 ! = xFe/xHeavy
      else if (refsolar.eq.'AGS05') then
         refsolchem = 2
         xFesol = 1.0666D-03
      else if (refsolar.eq.'AGSS09') then
         refsolchem = 3
         xFesol = 1.1953D-03
      else if (refsolar.eq.'GRID') then
         refsolchem = 4
         xFesol = 1.147417671020920d-3
      else
         stop 'Wrong solar composition, check $SOLCOMP variable'
      endif


      iheavy = nis-1 !53
***   Solar composition : WARNING xspsol (1-->iheavy)

      xsol = xspref(refsolchem,2)
      ysol = xspref(refsolchem,5)
      zsolar = 1.d0-xspref(refsolchem,2)-xspref(refsolchem,3)
     &     -xspref(refsolchem,4)-xspref(refsolchem,5)
      xheavysol = xspref(refsolchem,iheavy)
c  From Grevesse, 1998: solar Fe mass fraction = xFesol = 1.2589082d-3
c  According to Grevesse, Noels & Sauval 1996 (ASP Conf. Ser. 99, 117) XFesol = 1.2588215d-3
c  According to Grevesse & Noels 1993 (Origin and Evolution of Elements, ed. Prantzos, 
c  Vangioni-Flam & Casse) XFesol = 1.147417671020920d-3
c      xFesol = 1.2588215d-3  ! = xFe/xHeavy

      fracFe = xFesol/xheavysol
c According to Grevesse, Noels & Sauval 1996 (ASP Conf. Ser. 99, 117)
c log[N(Fe)/N(H)]sun = -4.5
c According to Grevesse & Noels 1993 (Origin and Evolution of Elements, ed. Prantzos, 
c  Vangioni-Flam & Casse, 15)
c log[N(Fe)/N(H)]sun = -4.5
      FeHsol = log10(xFesol/(56.d0*xsol)) !-4.5d0


***   primordial abundances

c   WMAP : Coc et al.(2004, apj, 600, 544) 
      xini = 0.7521d0
      yini = 0.2479d0
      zini = 0.d0
c   Coc et al.(2004, apj, 600, 544) : Y_0 = 0.2479(+-0.0004)
c   Cyburt et al.(2003 astro-ph/0302431) : Y_0 = 0.2484(+0.0004-0.0005)
c   Cuoco et al. (2003 astro-ph/0307213) : Y_0 = 0.2474(+0.0008-0.0005)
c   Scott Burles and David Tytler, 1998, ApJ et Izotov, 1999 : xini = 0.765d0
c   Bonifacio & Molaro, 1997 : yini = 0.235d0

***   constants
      pi = 3.1415926535d0
      pim2 = 2.d0*pi
      pim4 = 4.d0*pi
      pim8 = 8.d0*pi
      sec = 3.1557807d7
      c = 2.99792458d10
      g = 6.67259d-8
      boltz = 1.380658d-16
      h = 6.6260755d-27
      sig = 2.d0*pi**5*boltz**4/(15.d0*h**3*c*c)
      sigc = 4.d0*sig/c
      avn = 6.0221367d23
      rk = avn*boltz
      ech = 1.60217733d-19
      econv = ech*avn*1.d13

      pw23 = 2.d0/3.d0

      msun = 1.9891d33
      rsun = 6.959d10
      lsun = 3.851d33

      code_version = 2.80d0

      opt = '9'
      open (unit = 20, file = 'modini.bin',
     &     form = 'unformatted', status = 'old')
      
      read (5,*) znew,xhe4,diffvr,rotationtype
      print *,znew,xhe4,diffvr,rotationtype
      open (unit = 30, file = 'modini1.bin', form = 'unformatted',
     &     status = 'unknown') 
      call zchange (xhe4,20,30,opt)
      print *,'  input : modini.bin  --->  output : modini1.bin,',
     &     ' zinit.out'
      print *
      if (diffvr.gt.0.d0) then
         call binrot(rotationtype)
      endif

      end



************************************************************************

      SUBROUTINE EXTRACT

************************************************************************
* Extract binary data                                                  *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.rot'
      include 'evolcom.var'

      integer k,l,nshell
      double precision znew,mini,rini,zini,pindex,FeH

      character crzc*1

      common /param/ znew,mini,rini,zini,pindex,FeH


      write (*,100)
 100  format (/,2x,'enter shell number')
      read (5,*) nshell

      open (unit=20,file = 'modini.bin',form = 'unformatted',
     &     status = 'unknown')

      rewind (20)
      call rmod(20)
      close (20)

      write (*,*)
      write (*,300) model,totm,time/sec
 300  format ('  model #',i7,': mass = ',f8.5,', age = ',1pe10.4)
      write (*,*)

      open (unit=30,file = 'modini.shell',status = 'unknown')

      k = nshell
      if (crz(k).eq.-2) crzc = 'c'
      if (crz(k).eq.-1) crzc = 'a'
      if (crz(k).eq.-3) crzc = 's'
      if (crz(k).eq.3) crzc = 'S'
      if (crz(k).eq.2) crzc = 't'
      if (crz(k).eq.1) crzc = 'o'
      if (crz(k).eq.4) crzc = 'r'
      write (30,*) nshell
      write (30,*) model,nphase,totm,time,dtn
      write (30,*) neff,bin_version,nmod,nsp,mdredgeup,FeH,turnenv,
     &     rayon0
      write (30,*) k,crzc,u(k),r(k),lnf(k),lnt(k),lum(k)
      write (30,*) vu(k),vr(k),vlnf(k),vlnt(k),vlum(k),vfacc(k)
      write (30,*) m(k),dm(k),vhconv(k),vomega(k)
      write (30,*) ro(k),vro(k),p(k),vp(k),e(k),ve(k),venucl(k)
      write (30,*) s(k),vs(k),vsconv(k),vpvisc(k),tau(k)
      do l = 1,nis-1
         write (30,*) l,xsp(k,l)
      enddo

      write (*,*)
      print *,' input : modini.bin  --->  output : modini.shell'
      write (*,*)

      close (30)

      return
      end 



************************************************************************

      SUBROUTINE BINASC

************************************************************************
* Create an ascii file from a binary one                               *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.ion'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.rot'
      include 'evolcom.var'

      integer j,k,l
      double precision znew,mini,rini,zini,pindex,FeH

      character*1 crzc(nsh)

      common /param/ znew,mini,rini,zini,pindex,FeH


      open (unit=20,file = 'modini.bin',form = 'unformatted',
     &     status = 'unknown')

      rewind (20)
      call rmod(20)
      close (20)

      write (*,*)
      write (*,300) model,totm,time/sec
 300  format ('  model #',i7,': mass = ',f8.5,', age = ',1pe10.4)
      write (*,*)

      do k = 1,nmod
         if (crz(k).eq.-2) crzc(k) = 'c'
         if (crz(k).eq.-1) crzc(k) = 'a'
         if (crz(k).eq.-3) crzc(k) = 's'
         if (crz(k).eq.3) crzc(k) = 'S'
         if (crz(k).eq.2) crzc(k) = 't'
         if (crz(k).eq.1) crzc(k) = 'o'
         if (crz(k).eq.4) crzc(k) = 'r'
      enddo

      open (unit=30,file = 'modini1.asc',status = 'unknown')

      write (30,*)
      write (30,*) model,nphase,totm,time,dtn
      write (30,*) neff,bin_version,nmod,nsp,mdredgeup,FeH,turnenv,
     &     rayon0
      do k = 1,nmod
         write (30,*) k,crzc(k),u(k),r(k),lnf(k),lnt(k),lum(k)
         write (30,*) vu(k),vr(k),vlnf(k),vlnt(k),vlum(k),vfacc(k)
         write (30,*) m(k),vmueinv(k),dm(k),vhconv(k),vomega(k)
         write (30,*) ro(k),vro(k),p(k),vp(k),e(k),ve(k),venucl(k)
         write (30,*) s(k),vs(k),vsconv(k),vpvisc(k),tau(k)
         do l = 1,nis
            write (30,*) xsp(k,l)
         enddo
      enddo
      do j = 1,nis
         write (30,*) mtotlos(j)
      enddo
      write (*,*)
      print *,' input : modini.bin  --->  output : modini1.asc'
      write (*,*)

      close (30)

      return
      end 



************************************************************************

      SUBROUTINE BINBIN (flag)

************************************************************************
* Transform a binary file to another one                               *
* flag = 1 : upgrade binary file to version >= 2.80                    *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'

      integer k,flag

      double precision ynewsol,yoldsol,rathe,dely
      double precision zreal
      double precision znew,mini,rini,zini,pindex,FeH

      common /param/ znew,mini,rini,zini,pindex,FeH


      open (unit=20,file = 'modini.bin',form = 'unformatted',
     &     status = 'unknown')

      rewind (20)
      call rmod(20)
      close (20)

      write (*,*)
      write (*,300) model,totm,time/sec
 300  format ('  model #',i7,': mass = ',f8.5,', age = ',1pe10.4)
      write (*,*)

***   desired modifications follow if flag <> 1

      if (flag.ne.1) then
         write (6,*) 'Give the new total Y value:'
         read (5,*) ynewsol
         if (ynewsol.lt.1.d-5) goto 100
         yoldsol = xsp(1,4)+xsp(1,5)
         rathe = xsp(1,5)/yoldsol
         dely = ynewsol-yoldsol
         do k = 1,nmod
            xsp(k,2) = xsp(k,2)-dely
            xsp(k,4) = ynewsol*(1.d0-rathe)
            xsp(k,5) = ynewsol*rathe
         enddo
      else
         m(1) = 0.d0
         do k = 2,nmod
            if (dm(k-1).lt.0.d0.or.m(k).lt.m(k-1)) write (6,200) k,dm(k)
            dm(k-1) = abs(dm(k-1))
            m(k) = dm(k-1)+m(k-1)
         enddo
         dm(nmod) = 0.d0
      endif
***   
      zreal = 1.d0-xsp(nmod,2)-xsp(nmod,3)-xsp(nmod,4)-xsp(nmod,5)
      call zFeH (1,zreal,FeH,xsp(nmod,nis-1))

 100  open (unit=30,file = 'modini1.bin',form = 'unformatted',
     &     status = 'unknown')
      call wmod(30)
      close (30)

      write (6,150) zreal,FeH
      write (*,*)

 150  format (/,' input : modini.bin  --->  output : modini1.bin',/,
     &     ' Z = ',1pe10.4,', [Fe/H] = ',1pe11.4)
 200  format (' WARNING : shell [',i4,'] TOO SMALL !!!!  dm = ',1pe13.6) 
      return
      end 


************************************************************************

c      SUBROUTINE BINROT (diffvr,rotationtype)
      SUBROUTINE BINROT (rotationtype)

************************************************************************
* Generate inital rotation profile
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.transp'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.rot'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer iheavy
      integer i1,i,k,l,count,iBCE,type

 
      character*1 rotationtype

      double precision aa,zreal,xheavy
c      double precision auxs(nsh),urs(nsh),vtheta(nsh),thetas(nsh)
c     &     ,V_circ(nsh),xalpha(nsh),xlambdas(nsh),xKt(nsh),dift(nsh)
c     &     ,dife(nsh),Dhold(nsh),convdiff(nsh),abmuj(nsh),xnuvv(nsh)
c     &     ,xnum(nsh),ur_Ss,xmom_tots
      double precision convdiff(nsh)
      double precision omm,dmom
      double precision gmr(nsh)
      double precision znew,mini,rini,zini,pindex,FeH
      double precision xini,yini,xsol,ysol,zsolar,xFesol,
     &     FeHsol,fracFe
      double precision dift,dife,Dhold
      double precision abmurj,abmuj
      double precision vxpsi(nsh)
      double precision a1,a2,e1,e1m1

      common /difcirc/ dift(nsh),dife(nsh)
      common /calcDh/ Dhold(nsh)
      common /gradmuj/ abmurj(nsh),abmuj(nsh)
      common /param/ znew,mini,rini,zini,pindex,FeH
      common /bigbang/xini,yini,xsol,ysol,zsolar,xFesol,
     &     FeHsol,fracFe,iheavy

      
      
      open (unit = 20, file = 'modini1.bin',form = 'unformatted',
     &     status = 'unknown')

      call rmod(20)
      close (20)
      do k = 1,nmod
         if (crz(k).eq.-2.and.crz(k-1).ne.-2) iBCE = k
      enddo
      zreal = 1.d0-xsp(nmod,2)-xsp(nmod,3)-xsp(nmod,4)-xsp(nmod,5)
      xheavy = xsp(nmod,iheavy)
      call zFeH (1,zreal,FeH,xheavy)

      write (*,*)
      write (*,300) model,totm,time/sec,FeH
 300  format ('  model #',i7,': mass = ',f8.5,', age = ',1pe10.4,
     &     ', [Fe/H] = ',0pf7.4)
      write (*,*)
      
      omega_S = diffvr*1.d5/r(nmod)
c      print*,omega_S,diffvr,r(nmod)

      if (rotationtype.eq.'s') then
c..   Solid boby rotation
        print *,'Solid body rotation, omega =',omega_S
         do k = 1,nmod
            vomega(k) = omega_S
c            thetas(k) = 0.d0
            xpsis(k) = 0.d0
         enddo
         inits = .true.
      else if (rotationtype.eq.'d') then
         open (unit = 21, file = 'modang1.bin',form = 'unformatted',
     &     status = 'unknown')
         call rmodrot(21)
         close (21)

c..   Differential rotation
c..   According to Denissenkov & Vandenberg, 2003, ApJ, 598, 1254
cc         aa = r(nmod)*(diffvr*1.d5)**2/(3.d0*g*m(nmod))
c         aa = 0.24d0*1.d-3
c         open (unit = 26,file = 'prof_j_ZAMS_try',status = 'unknown')
c         do k = 2,nmod
c            gmr(k) = g*m(k)/(r(k)*r(k))
c            vomega(k) = sqrt(aa*3.d0*g*m(k)/r(k)**3)
c            jzams(k) = 2.d0*vomega(k)*r(k)**2/3.d0
c            print *,'k,omega',k,vomega(k),r(k)**2,jzams(k)
c            write(26,1003)k,jzams(k)
cc            vomega(k) = sqrt(aa*g*m(k)/r(k)**3)
c         enddo
c 1003    format (1x,i4,1x,1pe14.7)
c         close (26)
         do k = 2,nmod
            gmr(k) = g*m(k)/(r(k)*r(k))
         enddo
         print *,'Specific angular momentum conservation [1]',
     &        ' or hydro sim-based profile [2]?'
         read(5,'(i1)') type
         if (type.eq.2) then
            a1 = -69.8425d0
            a2 = 72.4075d0
            e1 = -0.412881d0
            e1m1 = e1 - 1.d0
         endif

         do k = 2,nmod
            grav(k) = gmr(k)/gmr(nmod)
c..   If starting at the end of 1st DUP (M = 0.85 Msun, Z = 0.0005)
c            vomega(k) = 8.d-6*(0.955869913d0*6.96d10/r(k))**2
c..   If starting at the bump (M = 0.85 Msun, Z = 0.0005)
c            vomega(k) = 8.d-6*(1.067377916d0*6.96d10/r(k))**2
c...  General case - Modified 21 january 2011
            if (k.ge.iBCE) then
               if (type.eq.1) then ! Specific angular momentum conservation
                  if (k.eq.iBCE) omega_S = (r(iBCE)**2*vomega(iBCE))
     $                 /r(nmod)**2
                  vomega(k) = omega_S *(r(nmod)/r(k))**2
                  xpsis(k) = -2.d0/3.d0*r(k)**2/grav(k)*vomega(k)*
     $                 (vomega(k)-vomega(k-1))/((r(k)-r(k-1))*r(nmod)
     $                 *omega_S**2)
               else if (type.eq.2) then ! hydro-sim-based profile
                  vomega(k) = vomega(iBCE)*(a1 + a2 * (r(k)/r(nmod))
     $                 **e1)/(a1 + a2 * (r(iBCE)/r(nmod))**e1)

               endif
            endif
         enddo
c         vomega(1) = vomega(2)
c         do k = 2,nmod
c            grav(k) = gmr(k)/gmr(nmod)
c            thetas(k) = 2.d0/3.d0*r(k)**2/grav(k)*(vomega(k)-
c     &           vomega(k-1))/((r(k)-r(k-1))*r(nmod))
c             write(662,*)k,vomega(k),thetas(k)
c
c         enddo
c         thetas(1) = thetas(2)
         inits = .false.
c         inits = .true.
      endif

      open (unit = 30, file = 'modini2.bin',form = 'unformatted',
     &     status = 'unknown')
      call wmod(30)
      close(30)

      if (crz(1).gt.0) then
         ndbold = 1
      else
         k = 2
         do while (crz(k).lt.-1) 
            k = k+1
         enddo
         ndbold = k
      endif

      count = 0
      do k = nmod,2,-1
         count = count+1
         if (crz(k).lt.-1.and.crz(k-1).gt.0.and.count.gt.1) 
     &        ndtold = k
      enddo

      xmom_tots = 0.d0

      do i = 2,nmod
         i1 = i-1
         omm = 0.5d0*(vomega(i)+vomega(i1))
         dmom = (r(i)**2+r(i1)**2+r(i)*r(i1))/3.d0*omm*dm(i1)
         xmom_tots = xmom_tots+dmom
         ur_Ss = 0.d0
         auxs(i) = 0.d0
         urs(i) = 0.d0
c         vtheta(i) = thetas(i)
c         vtheta(i) = 0.d0
c         thetas(i) = 0.d0
         V_circ(i) = 0.d0
         xalpha(i) = 0.d0
         xlambdas(i) = 0.d0
         xKt(i) = 0.d0
         dift(i) = 0.d0
         dife(i) = 0.d0
         Dhold(i) = 0.d0
         convdiff(i)= 0.d0
         abmuj(i) = 0.d0
         xnuvv(i) = 0.d0
         xnum(i) = 0.d0
      enddo
      auxs(1) = 0.d0
      urs(1) = 0.d0
c      vtheta(1) = 0.d0
c      thetas(1) = 0.d0
c      vtheta(1) = vtheta(2)
      V_circ(1) = 0.d0
      xalpha(1) = 0.d0
      xlambdas(1) = 0.d0
      xKt(1) = 0.d0
      dift(1) = 0.d0
      dife(1) = 0.d0
      Dhold(1) = 0.d0
      convdiff(1)= 0.d0
      abmuj(1) = 0.d0
      xnuvv(1) = 0.d0
      xnum(1) = 0.d0

      open (unit = 40, file = 'modang2.bin',form = 'unformatted',
     &     status = 'unknown')

      write (40) ndtold,ndbold,ur_Ss,xmom_tots,inits,
c     &     (auxs(k),urs(k),vtheta(k),thetas(k),V_circ(k),
     &     (auxs(k),urs(k),vxpsi(k),xpsis(k),V_circ(k),
     &     xalpha(k),xlambdas(k),xKt(k),dift(k),dife(k),
     &     Dhold(k),convdiff(k),abmuj(k),xnuvv(k),xnum(k),
     &     k = 1,nmod),
     &     ((xsp(k,l),l=1,nis-1),k=1,nmod)
      close(40)

      write (*,*)
      print *,'  input modini1.bin,modang1.bin -->', 
     &     'output modini2.bin, modang2.bin'
      write (*,*)

      end


************************************************************************

      SUBROUTINE HOMOL

************************************************************************
* Calculate an initial model by homology                               *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.data'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.rot'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'


      integer last,nfile,iheavy
      integer i,k

      character modbini*60,cgzip*160,dirinit*100
      character name*100,nameout*100,cmass*5,cname*100

      logical lowmetal

      double precision t0,vt0,ro0,vro0,zreal
      double precision cM0,cM,cR0,cR,rM,rR,rmue,aa,bb,x1,x2
      double precision rmass(100)
      double precision znew,mini,rini,zini,pindex,FeH,xheavy
      double precision xini,yini,xsol,ysol,zsolar,xFesol,
     &     FeHsol,fracFe

      common /param/ znew,mini,rini,zini,pindex,FeH
      common /bigbang/xini,yini,xsol,ysol,zsolar,xFesol,
     &     FeHsol,fracFe,iheavy

      dimension t0(nsh),vt0(nsh),ro0(nsh),vro0(nsh)
      dimension name(100),nameout(100)

      lowmetal = zini.gt.0.d0.and.zini.le.1.d-5

      call getenv('DIR_INIT',dirinit)
      if (lowmetal) then
         dirinit = trim(dirinit) // 'Z0/'
         call system ("cd $DIR_INIT/Z0 ; ls sini* > /tmp/sini.dir")
      else
         dirinit = trim(dirinit)
         call system (
     &        "cd $DIR_INIT ; ls sini* | grep -v Z0 > /tmp/sini.dir")
      endif
      open (unit = 50, file = '/tmp/sini.dir', status = 'old')

      nfile = 0
      do 10 i = 1,100
         read (50,*,end=20) cname
         
         if (cname(7:7).eq.'.') then
            cmass = cname(6:9)
         else
            cmass = cname(6:11)
         endif
         nfile = nfile+1
         read (cmass,'(0pf5.2)') rmass(nfile)
         name(nfile) = cname
 10   continue

      
 20   call sort1 (rmass,name,nameout,nfile)
      i = 1
      if (mini.le.rmass(1)) goto 30
      do i = 1,nfile-1
         if (mini.ge.rmass(i).and.mini.lt.rmass(i+1)) goto 30
      enddo

 30   modbini = nameout(i)
      write (*,*)
      write (*,'("  homology from model : ",a25)') modbini

      last = 60
      do while (last.gt.1.and.modbini(last:last).eq.' ' )
         last = last - 1
      enddo
      if (modbini(last:last) .eq. ' ' ) stop 'blank file name'
      if (modbini(last-1:last).eq.'gz') then
         last = last-3
         cgzip = 'gunzip -f ' // trim(dirinit) // trim(modbini)
         call system (cgzip)
      endif
      open (unit = 40, file = trim(dirinit) // modbini(1:last),
     &     form = 'unformatted',status = 'unknown')
      call rmod(40)
      close (40)
      cgzip = 'gzip ' // trim(dirinit) // modbini(1:last)
      call system (cgzip)
      call system ("rm -f /tmp/sini.dir")

      model = 0
      nphase = 1
c      times0 = 0.d0
      mdredgeup = totm*msun
      rayon0 = 0.d0
      turnenv = 0.d0
      taccend = 0.d0
      cM0 = totm
      cM = mini
      totm = cM
      time = 0.d0
      dtn = 1.d3*sec
      cR0 = r(nmod)
      cR = cR0*sqrt(cM/cM0)
      if (cM0.gt.15.d0) cR = cR0*(cM/cM0)**(2.d0/5.d0)
      rM = cM/cM0
      rR = cR/cR0
      rmue = 1.d0
      aa = 1.d0
      bb = 3.5d0
      nmod1 = nmod-1
      time = 0.d0
      do i = 1,nmod
         m(i) = m(i)*rM
         u(i) = u(i)*rR
         vu(i) = vu(i)*rR
         r(i) = r(i)*rR
         vr(i) = vr(i)*rR
         t0(i) = exp(lnt(i))
         t(i) = t0(i)*rmue*rM/rR
         lnt(i) = log(t(i))
         vt0(i) = exp(vlnt(i))
         vt(i) = vt0(i)*rmue*rM/rR
         vlnt(i) = log(vt(i))
         if (t0(i).lt.5.d4) then
            aa = 0.74d0
            bb = 10.7d0
         endif
         x1 = 3.d0+bb-aa
         x2 = 3.d0*aa-bb
         lum(i) = lum(i)*rM**x1*rR**x2
         vlum(i) = vlum(i)*rM**x1*rR**x2
         ro0(i) = ro(i)
         ro(i) = ro0(i)*rM/(rR*rR*rR)
         vro0(i) = vro(i)
         vro(i) = vro0(i)*rM/(rR*rR*rR)
         p(i) = p(i)*rM*rM/(rR*rR)**2
         vp(i) = vp(i)*rM*rM/(rR*rR)**2
         lnf(i) = lnf(i)+log(ro(i)/ro0(i))+1.5d0*log(t0(i)/t(i))
         vlnf(i) = vlnf(i)+log(vro(i)/vro0(i))+1.5d0*log(vt0(i)/vt(i))
c         mue(i) = mue0(i)*rmue
c         mui = 1.282401013417d0
c         mue(i) = 1.d0/mu(i)-1.d0/mui
c         mue(i) = 1.d0/mue(i)
c         mue(i) = min(1.d10,abs(mue(i)))
c         lnf(i) = log(ro(i))-1.5d0*log(t(i))+18.03d0-log(mue(i))
c         vlnf(i) = log(vro(i))-1.5d0*log(vt(i))+18.03d0-log(mue(i))
c         lnf(i) = lnf(i)+100.d0
c         vlnf(i) = vlnf(i)+100.d0
         vomega(i) = 0.d0
         pvisc(i) = 0.d0
      enddo
      do i = 1,nmod1
         dm(i) = m(i+1)-m(i)
      enddo
      dm(nmod) = 0.d0
      zreal = 1.d0-xsp(nmod,2)-xsp(nmod,3)-xsp(nmod,4)-xsp(nmod,5)
      xheavy = xsp(nmod,iheavy)      
      call zFeH (1,zreal,FeH,xheavy)

      write (*,1000) totm,time/sec,dtn/sec,FeH

      open (unit = 60, file = 'modini1.bin',
     &     form = 'unformatted',status = 'unknown')
c     &     form = 'unformatted',convert = 'big_endian',

      bin_version = 2.36d0
      call wmod(60)
      znew = zini
      if (znew.lt.0) then
         write (*,401)
 401     format (/,5x,'Choose your option:',//
     &        3x,'o Anders & Grevesse 1989            [AG89]',/,
     &        3x,'o Grevesse & Noels 1993             [GN93]',/,
     &        3x,'o Grevesse, Noesl & Sauval 1996     [GNS96]',/,
     &        3x,'o Grevesse & Sauval 1998            [GS98]',/,
     &        3x,'o Asplund, Grevesse & Sauval 2005   [AGS05]',//,
     &        '    ---> Enter option',5x,a5)
         read (5,'(a6)') refsolar
         call solx(refsolar,znew)
      endif
      call zchange (60,60,'1')
      close(60)

c      do k = 1,nmod
c         crz(k) = 4
c         if (crz(k).lt.-1) crz(k) = -2
c      enddo

      print *,'   --->  output : modini1.bin'


 1000 format (/,2x,'new model: Mini = ',f8.5,'; age = ',1pe10.4,
     &     ' yr, dt = ',1pe10.4,' yr, [Fe/H] = ',0pf7.4)

      end 


************************************************************************

      SUBROUTINE ZCHANGE (xhe4,iu1,iu2,opt)

************************************************************************
* Change of global metallicity for a given binary file                 *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.data'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.spec'
      include 'evolcom.teq'


      integer iheavy,iu1,iu2
      integer i,k,l
      integer testLi
      integer refsolchem

      character*1 deut,calpha
      character*1 opt, opt2
      double precision xini,yini,zslope,x,y
      double precision sumx,zreal
      double precision xFesol,FeHsol,fracFe,xFenew
      double precision coefalpha, coefnoalpha, zalpha, znoalpha
      double precision OFechoice,NaFechoice,CFechoice,numalpha
      double precision znew,mini,rini,zini,pindex,FeH
      double precision xsol,ysol,zsolar,h2h1,he3he4,hxsp
      double precision xhe4,xo16,xc12,xna23,xli7,nLi,xh1
      double precision sum

      common /param/ znew,mini,rini,zini,pindex,FeH
      common /bigbang/xini,yini,xsol,ysol,zsolar,xFesol, 
     &     FeHsol,fracFe,iheavy

      dimension hxsp(nsp)

      ih1 = 2
      ih2 = 3
      ihe3 = 4
      ihe4 = 5

      numalpha = 1.d0

      open (unit = 31,file = 'zini.out',status = 'unknown')
 
      if (refsolar.eq.'GN93') refsolchem = 1
      if (refsolar.eq.'AGS05') refsolchem = 2
      if (refsolar.eq.'AGSS09') refsolchem = 3
      if (refsolar.eq.'GRID') refsolchem = 4
     
      write (31,99) refsolar
 99   format ('  o solar reference composition : ',a6)
      calpha = 'n'
      deut = 'n'

      rewind (iu1)
      call rmod(IU1)

      write (*,*)
      write (*,300) model,totm,time/sec
 300  format ('  model #',i7,': mass = ',f8.5,', age = ',1pe10.4)
      write (*,*)

      if (znew.le.1.d-5) then
         h2h1 = 2.6d-5
         he3he4 = 1.04d-5*xini/yini
c            Coc et al.(2003 astro-ph/0309480) : 
c                D/H = 2.6d-5
c                3He/H = 1.04d-5
c            Cyburt et al.(2003 astro-ph/0302431) : 
c                D/H = 2.74(+0.26-0.16)d-5
c                3He/H = (9.30(+1.00-0.67)d-6
c            Cuoco et al. (2003 astro-ph/0307213) : 
c                D/H = 2.56(+0.35-0.24)d-5
c                3He/H = 0(0.99+0.08-0.07)d-5
c            Scott Burles and David Tytler, 1998, ApJ 
c                h2h1 = 6.79d-5
c            Scott Burles and David Tytler, 1998, ApJ 
c                he3he4 = 5.1d-5*xini/yini
      else
         h2h1 = 6.8002312d-05
         he3he4 = 1.0644011d-04
      endif

      zslope = (xsol-xini)/zsolar

      x = zslope*znew+xini
      y = 1.d0-x-znew
      hxsp(1) = 1.d-50
      hxsp(2) = x/(h2h1+1.d0)
      hxsp(3) = x*h2h1/(h2h1+1.d0)
      hxsp(4) = y*he3he4/(he3he4+1.d0)
      hxsp(5) = y/(he3he4+1.d0)
      if (deut.eq.'y') then
         hxsp(2) = hxsp(2)-hxsp(3)
         hxsp(4) = hxsp(4)+hxsp(3)
         hxsp(3) = 1.d-30
      endif

      ili7 = 7
      ic12 = 13
      io16 = 20
      ine20 = 26
      img24 = 33
      isi28 = 41
      ina23 = 31
      
      testLi = 0
      if ((znew/x).lt.7.81245d-03) then
         testLi = 1
         nLi = 4.14954d-10
         xli7 = 7.0d0*x*nLi
         hxsp(7) = xli7
         write (*,111)
         write (31,111)
 111     format (' o Impose Li abundance to the Spite"s Plateau ',
     &        'value = 2.618 for [Fe/H]<-0.5')
      else
         hxsp(7) = xspref(refsolchem,7)*znew/zsolar
      endif

      hxsp(6) = xspref(refsolchem,6)*znew/zsolar

      zreal = hxsp(6)+hxsp(7)
      do i = 8,iheavy
         hxsp(i) = xspref(refsolchem,i)*znew/zsolar
         if (xspref(refsolchem,i).lt.1.d-30.or.hxsp(i).lt.1.d-30) 
     &        hxsp(i) = 1.d-50
         zreal = zreal+hxsp(i)
      enddo
      hxsp(nis) = 1.d-50
      hxsp(nsp) = 1.d-50
      hxsp(5) =  1.d0-zreal-hxsp(2)-hxsp(3)-hxsp(4)
c..   zreal is different from znew if Li has been set to its Plateau value


*________________________
***   Alpha enhancement
*------------------------
      
c     In the solar-scaled mixture
      xFenew = fracFe*hxsp(iheavy)
      print *,'%%%%%%% xFenew %%%%%%%%%',xFenew
      FeH = log10(xFenew/(hxsp(2)*56.d0))-FeHsol

c.. Alpha enhancement
       
      if (calpha.eq.'y') then
         zalpha = 0.d0
         znoalpha = 0.d0
         do i = ihe3,iheavy-1
            if (i.eq.io16.or.i.eq.ine20.or.i.eq.img24.or.i.eq.isi28) 
     &           then
               zalpha = zalpha+hxsp(i)
            else
               znoalpha = znoalpha+hxsp(i)
            endif
         enddo
         if (testLi.eq.1) then
            xli7 = hxsp(ili7)
            znoalpha = znoalpha-xli7
         endif

c  mass fractions of non-alpha elements decreased by a factor coefnoalpha 
c  mass fractions of alpha elements increased by a factor coefalpha 
c  mass fraction of H and Heavy(= Fe) are NOT modified so that 
c  [alpha/Fe] and [Fe/H] are kept constant but Z changes
         coefalpha = numalpha
         coefnoalpha = (znoalpha+zalpha*(1.d0-numalpha))/znoalpha
         if (coefnoalpha.le.0.d0) then
            print *,'too large enrichment, mass fractions < 0 !!'
            return
         endif
         sumx = 0.d0
         do i = ihe3,iheavy-1
            if (i.eq.io16.or.i.eq.ine20.or.i.eq.img24.or.i.eq.isi28)
     &           then
               hxsp(i) = coefalpha*hxsp(i)
            else
               hxsp(i) = coefnoalpha*hxsp(i)
            endif
            sumx = sumx+hxsp(i)
         enddo
         if (testLi.eq.1) then
            sumx = sumx-hxsp(ili7)+xli7
            hxsp(ili7) = xli7
         endif
c..  zreal has changed because He mass fraction has been modified
         zreal = sumx+hxsp(iheavy)-hxsp(ihe3)-hxsp(ihe4)
      endif


*__________________________________________________________
***   Modification of chosen mass fractions for comparison 
***   with observations
*----------------------------------------------------------

c      xhe4 = hxsp(ihe4)
      xc12 = hxsp(ic12)
      xo16 = hxsp(io16)
      xna23 = hxsp(ina23)
c$$$      if (opt.eq.'9') then
c$$$ 400     write (*,401)
c$$$ 401     format (/,5x,'Choose your option:',//
c$$$     &        3x,'o Change 4He fraction of mass          [1]',/,
c$$$     &        3x,'o Change [12C/Fe] abundance            [2]',/,
c$$$     &        3x,'o Change [16O/Fe] abundance            [3]',/,
c$$$     &        3x,'o Change [23Na/Fe] abundance           [4]',/,
c$$$     &        3x,'o Exit                                 [5]',//,
c$$$     &        '    ---> Enter option',5x,a1)
c$$$         read (5,200) opt2
c$$$         
c$$$         if (opt2.ne.'1'.and.opt2.ne.'2'.and.opt2.ne.'3'.and.
c$$$     &        opt2.ne.'4'.and.opt2.ne.'5') goto 400
c$$$         if (opt2.eq.'1') then
c$$$ 409        write (*,410)
c$$$ 410        format (5x,'Choose 4He fraction of mass:',5x,a1)
c$$$            read (5,*) xhe4
c$$$            if (xhe4.gt.1.d0-zreal) then
c$$$               write (*,435)
c$$$               goto 409
c$$$            endif
c$$$            goto 400
c$$$         endif
c$$$         if (opt2.eq.'2') then
c$$$ 419        write (*,420)
c$$$ 420        format (5x,'Choose [12C/Fe] abundance:',5x,a1)
c$$$            read (5,*) CFechoice
c$$$            xc12 = 10.d0**CFechoice*xspref(refsolchem,ic12)*xFenew
c$$$     &           /xFesol
c$$$            if (xc12.gt.hxsp(5)) then
c$$$               write (*,440)
c$$$               goto 419
c$$$            endif
c$$$            goto 400
c$$$         endif
c$$$         if (opt2.eq.'3') then
c$$$ 429        write (*,430)
c$$$ 430        format (5x,'Choose [16O/Fe] abundance:',5x,a1)
c$$$            read (5,*) OFechoice
c$$$            xo16 = 10.d0**OFechoice*xspref(refsolchem,io16)*xFenew
c$$$     &           /xFesol
c$$$            if (xo16.gt.hxsp(5)) then
c$$$               write (*,440)
c$$$               goto 429
c$$$            endif
c$$$            goto 400
c$$$         endif
c$$$         if (opt2.eq.'4') then
c$$$ 431        write (*,432)
c$$$ 432        format (5x,'Choose [23Na/Fe] abundance:',5x,a1)
c$$$            read (5,*) NaFechoice
c$$$            xna23 = 10.d0**NaFechoice*xspref(refsolchem,ina23)*
c$$$     &           xFenew/xFesol
c$$$            if (xna23.gt.hxsp(5)) then
c$$$               write (*,440)
c$$$               goto 431
c$$$            endif
c$$$            goto 400
c$$$         endif
c$$$         
c$$$ 435     format (5x,'Invalid choice for 4He abundance!')
c$$$ 440     format (5x,'Abundance is larger than chosen Z, choose new 
c$$$     &        abundance!')

c..   If option 1 not selected : adjust Y according to the chosen modifications 
c..     so all [X/Fe] and [X/H] are conserved. 
c..   else H is modified. In all cases the metallicity Z is modified

c$$$         if (opt2.eq.'5') then
      print *
      if (xc12.ne.hxsp(ic12)) then
         write (*,1201) CFechoice
         write (31,1201) CFechoice
      endif
      if (xo16.ne.hxsp(io16)) then
         write (*,1202) OFechoice
         write (31,1202) OFechoice
      endif
      if (xna23.ne.hxsp(ina23)) then
         write (*,1203) NaFechoice
         write (31,1203) NaFechoice
      endif
      if (xhe4.eq.hxsp(ihe4)) then
         hxsp(ihe4) = hxsp(ihe4)+hxsp(io16)-xo16+hxsp(ic12)-xc12+
     &        hxsp(ina23)-xna23
         hxsp(io16) = xo16
         hxsp(ic12) = xc12
         hxsp(ina23) = xna23
      else
         xh1 = hxsp(ih1)
         hxsp(ih1) = hxsp(ih1)+hxsp(io16)-xo16+hxsp(ic12)-xc12+
     &        hxsp(ina23)-xna23+hxsp(ihe4)-xhe4
         hxsp(ihe4) = xhe4 
         hxsp(io16) = xo16
         hxsp(ic12) = xc12
         hxsp(ina23) = xna23
         write (*,1204) xh1,hxsp(ih1)
         write (31,1204) xh1,hxsp(ih1)
         write (*,1205) xhe4,hxsp(ihe4)
         write (31,1205) xhe4,hxsp(ihe4)
      endif

 1201 format (' o [12C/Fe]  imposed to ',0pf6.3)
 1202 format (' o [16C/Fe]  imposed to ',0pf6.3)
 1203 format (' o [23Na/Fe] imposed to ',0pf6.3)
 1204 format (' o X(1H)  : ',1pe11.4,', --> ',1pe11.4)
 1205 format (' o X(4He) : ',1pe11.4,', --> ',1pe11.4)

      zreal = 0.d0
      do k = 1,nis
         zreal = zreal+hxsp(k)
         mtotlos(k) = 1.d-50
      enddo
      if (abs(1.d0-zreal).gt.1.d-15) then
         print *,'Normalisation problem in Zchange : sumx = ',zreal
         stop
      endif
      zreal = zreal-hxsp(2)-hxsp(3)-hxsp(4)-hxsp(5)

      FeH = log10(xFenew/(hxsp(2)*56.d0))-FeHsol
      write (*,1000) zreal,FeH
      write (*,*)

      call fracmass (hxsp,xFenew)

      write (*,8999) xsol,ysol,zsolar,log10(hxsp(7)/hxsp(2)/7.d0)+12.d0

      write (31,8999) xsol,ysol,zsolar,log10(hxsp(7)/hxsp(2)/7.d0)+12.d0
      write (31,9000) bin_version,totm,time/sec,zreal,FeH
      write (31,*)'                                X            Xsol'
      write (31,9001) hxsp(1),xspref(refsolchem,1)
      write (31,9002) hxsp(2),xspref(refsolchem,2)
      write (31,9003) hxsp(3),xspref(refsolchem,3)
      write (31,9004) hxsp(4),xspref(refsolchem,4)
      write (31,9005) hxsp(5),xspref(refsolchem,5)
      write (31,9006) hxsp(6),xspref(refsolchem,6)
      write (31,9007) hxsp(7),xspref(refsolchem,7)
      write (31,9008) hxsp(8),xspref(refsolchem,8)
      write (31,9009) hxsp(9),xspref(refsolchem,9)
      write (31,9010) hxsp(10),xspref(refsolchem,10)
      write (31,9011) hxsp(11),xspref(refsolchem,11)
      write (31,9012) hxsp(12),xspref(refsolchem,12)
      write (31,9013) hxsp(13),xspref(refsolchem,13)
      write (31,9014) hxsp(14),xspref(refsolchem,14)
      write (31,9015) hxsp(15),xspref(refsolchem,15)
      write (31,9016) hxsp(16),xspref(refsolchem,16)
      write (31,9017) hxsp(17),xspref(refsolchem,17)
      write (31,9018) hxsp(18),xspref(refsolchem,18)
      write (31,9019) hxsp(19),xspref(refsolchem,19)
      write (31,9020) hxsp(20),xspref(refsolchem,20)
      write (31,9021) hxsp(21),xspref(refsolchem,21)
      write (31,9022) hxsp(22),xspref(refsolchem,22)
      write (31,9023) hxsp(23),xspref(refsolchem,23)
      write (31,9024) hxsp(24),xspref(refsolchem,24)
      write (31,9025) hxsp(25),xspref(refsolchem,25)
      write (31,9026) hxsp(26),xspref(refsolchem,26)
      write (31,9027) hxsp(27),xspref(refsolchem,27)
      write (31,9028) hxsp(28),xspref(refsolchem,28)
      write (31,9029) hxsp(29),xspref(refsolchem,29)
      write (31,9030) hxsp(30),xspref(refsolchem,30)
      write (31,9031) hxsp(31),xspref(refsolchem,31)
      write (31,9032) hxsp(32),xspref(refsolchem,32)
      write (31,9033) hxsp(33),xspref(refsolchem,33)
      write (31,9034) hxsp(34),xspref(refsolchem,34)
      write (31,9035) hxsp(35),xspref(refsolchem,35)
      write (31,9036) hxsp(36),xspref(refsolchem,36)
      write (31,9037) hxsp(37),xspref(refsolchem,37)
      write (31,9038) hxsp(38),xspref(refsolchem,38)
      write (31,9039) hxsp(39),xspref(refsolchem,39)
      write (31,9040) hxsp(40),xspref(refsolchem,40)
      write (31,9041) hxsp(41),xspref(refsolchem,41)
      write (31,9042) hxsp(42),xspref(refsolchem,42)
      write (31,9043) hxsp(43),xspref(refsolchem,43)
      write (31,9044) hxsp(44),xspref(refsolchem,44)
      write (31,9045) hxsp(45),xspref(refsolchem,45)
      write (31,9046) hxsp(46),xspref(refsolchem,46)
      write (31,9047) hxsp(47),xspref(refsolchem,47)
      write (31,9048) hxsp(48),xspref(refsolchem,48)
      write (31,9049) hxsp(49),xspref(refsolchem,49)
      write (31,9050) hxsp(50),xspref(refsolchem,50)
      write (31,9051) hxsp(51),xspref(refsolchem,51)
      write (31,9052) hxsp(52),xspref(refsolchem,52)
      write (31,9053) hxsp(53),xspref(refsolchem,53)
      write (31,9054) hxsp(54)
      write (31,9055) hxsp(55)
      close (31)

      sum = 0.d0
      do l = 1,nis
         sum = sum+hxsp(l)
         do k = 1,nmod
            xsp(k,l) = hxsp(l)
         enddo
      enddo
      print *,'sum elets',sum

      rewind (iu2)
      call wmod(iu2)

      close (iu2)
      close (iu1)


 1000 format (3x,'computed Z = ',1pe13.7,', [Fe/H] = ',0pf7.4)
 8999 format (3x,'Xsol = ',f8.6,', Ysol = ',f8.6,', Zsolar = ',f8.6,
     &     ', N(Li) = ',f8.4)
 9000 format (/,2x,'INITIAL CHEMICAL SPECIES - version ',f5.3,/,2x,
     &     'mass : ',f8.4,', age = ',1pe11.5,', Z = ',1pe12.6,
     &     ', [Fe/H] = ',0pf7.4/)
 9001 format ('   1 NEUT    1.00000   0.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9002 format ('   2 PROT    1.00000   1.',3x,1pe11.5,3x,1pe11.5)
 9003 format ('   3 DEUT    2.00000   1.',3x,1pe11.5,3x,1pe11.5)
 9004 format ('   4 HE  3   3.00000   2.',3x,1pe11.5,3x,1pe11.5)
 9005 format ('   5 HE  4   4.00000   2.',3x,1pe11.5,3x,1pe11.5)
 9006 format ('   6 LI  6   6.00000   3.',3x,1pe11.5,3x,1pe11.5)
 9007 format ('   7 LI  7   7.00000   3.',3x,1pe11.5,3x,1pe11.5)
 9008 format ('   8 BE  7   7.00000   4.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9009 format ('   9 B   8   8.00000   5.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9010 format ('  10 BE  9   9.00000   4.',3x,1pe11.5,3x,1pe11.5)
 9011 format ('  11 B  10  10.00000   5.',3x,1pe11.5,3x,1pe11.5)
 9012 format ('  12 B  11  11.00000   5.',3x,1pe11.5,3x,1pe11.5)
 9013 format ('  13 C  12  12.00000   6.',3x,1pe11.5,3x,1pe11.5)
 9014 format ('  14 C  13  13.00000   6.',3x,1pe11.5,3x,1pe11.5)
 9015 format ('  15 N  13  13.00000   7.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9016 format ('  16 C  14  14.00000   6.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9017 format ('  17 N  14  14.00000   7.',3x,1pe11.5,3x,1pe11.5)
 9018 format ('  18 N  15  15.00000   7.',3x,1pe11.5,3x,1pe11.5)
 9019 format ('  19 O  15  15.00000   8.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9020 format ('  20 O  16  16.00000   8.',3x,1pe11.5,3x,1pe11.5)
 9021 format ('  21 O  17  17.00000   8.',3x,1pe11.5,3x,1pe11.5)
 9022 format ('  22 O  18  18.00000   8.',3x,1pe11.5,3x,1pe11.5)
 9023 format ('  23 F  18  18.00000   9.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9024 format ('  24 F  19  19.00000   9.',3x,1pe11.5,3x,1pe11.5)
 9025 format ('  25 F  20  20.00000   9.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9026 format ('  26 NE 20  20.00000  10.',3x,1pe11.5,3x,1pe11.5)
 9027 format ('  27 NE 21  21.00000  10.',3x,1pe11.5,3x,1pe11.5)
 9028 format ('  28 NE 22  22.00000  10.',3x,1pe11.5,3x,1pe11.5)
 9029 format ('  29 NA 22  22.00000  11.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9030 format ('  30 NE 23  23.00000  10.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9031 format ('  31 NA 23  23.00000  11.',3x,1pe11.5,3x,1pe11.5)
 9032 format ('  32 NA 24  24.00000  11.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9033 format ('  33 MG 24  24.00000  12.',3x,1pe11.5,3x,1pe11.5)
 9034 format ('  34 NA 25  25.00000  11.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9035 format ('  35 MG 25  25.00000  12.',3x,1pe11.5,3x,1pe11.5)
 9036 format ('  36 MG 26  26.00000  12.',3x,1pe11.5,3x,1pe11.5)
 9037 format ('  37 AL26M  26.00000  13.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9038 format ('  38 AL26G  26.00000  13.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9039 format ('  39 MG 27  27.00000  12.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9040 format ('  40 AL 27  27.00000  13.',3x,1pe11.5,3x,1pe11.5)
 9041 format ('  41 SI 28  28.00000  14.',3x,1pe11.5,3x,1pe11.5)
 9042 format ('  42 SI 29  29.00000  14.',3x,1pe11.5,3x,1pe11.5)
 9043 format ('  43 SI 30  30.00000  14.',3x,1pe11.5,3x,1pe11.5)
 9044 format ('  44 P  31  31.00000  15.',3x,1pe11.5,3x,1pe11.5)
 9045 format ('  45 S  32  32.00000  16.',3x,1pe11.5,3x,1pe11.5)
 9046 format ('  46 S  33  33.00000  16.',3x,1pe11.5,3x,1pe11.5)
 9047 format ('  47 S  34  34.00000  16.',3x,1pe11.5,3x,1pe11.5)
 9048 format ('  48 S  35  35.00000  16.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9049 format ('  49 CL 35  35.00000  17.',3x,1pe11.5,3x,1pe11.5)
 9050 format ('  50 S  36  36.00000  16.',3x,1pe11.5,3x,1pe11.5)
 9051 format ('  51 CL 36  36.00000  17.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9052 format ('  52 CL 37  37.00000  17.',3x,1pe11.5,3x,1pe11.5)
 9053 format ('  53 HEAVY  38.00000  18.',3x,1pe11.5,3x,1pe11.5)
 9054 format ('  54 CAPTN  39.00000  18.',3x,1pe11.5,3x,
     &     '1.00000d-50  *')
 9055 format ('  55 OOOOO 999.99999  99.',3x,1pe11.5,3x,
     &     '1.00000d-50  *')

      return
      end


************************************************************************

      SUBROUTINE ZFEH (flag,z,FeH,xheavy)

************************************************************************
* This routine  computes the metallicity [Fe/H] corresponding to       *
* a given Z and xHeavy (flag = 1) or the reverse (flag = 2)            *
*                                                                      *
* The assumed formulae to perform this computation are                 *
*                                                                      *
* o if the mass fraction of heavy is unknown one assumes Z/Fe = cte    *
*      Z = XINI * 10**[Fe/H] / (alpha * 10**[Fe/H] + Xsol/Zsolar)      *
* o if the mass fraction of heavy is known one assumes xHeavy/Fe = cte *
*   Z = XINI * 10**[Fe/H] / (alpha * 10**[Fe/H] + Xsol/Zsolar) , where *
*    alpha = (XINI - Xsol) / Zsolar , and                              *
*    XINI = Hydrogen mass fraction at t = 0 (when Z = 0),              *
*    Xsol and Zsolar being the Hydrogen mass fraction and metallicity Z*
*       at the present Sun (in case of no diffusion), respectively.    *
*                                                                      *
* This program needs, as input: iz and [Fe/H] or Z/xHeavy              *
* The only output is: Z or [Fe/H], depending on the choice for iz      *
************************************************************************

      implicit none

      integer iz,iiz,flag,iheavy

      double precision FeH,z,xheavy,xFe
      double precision xini,yini,xsol,ysol,zsolar,xFesol,fracFe,FeHsol
      double precision FeH10,div,x,alpha

      common /bigbang/xini,yini,xsol,ysol,zsolar,xFesol,FeHsol,fracFe,
     &     iheavy

      alpha = (xini-xsol)/zsolar

*** Input values:

      iiz = 0
      iz = flag
      if (flag.eq.1.or.flag.eq.2) goto 10
      write (*,1000)
      read (*,1100) iiz

      if (iiz.eq.1) then
         write (*,1200)
         read (*,1300) FeH
      endif
      if (iiz.eq.2) then
         write (*,1400)
         read (*,1500) z
      endif
      if (iiz.ne.1.and.iiz.ne.2) then
         write (*,999)
         return
      endif
      iz = abs(iiz-3)

*** Computations:
 10   if (iz.eq.1.or.iz.eq.2) then
         if (xheavy.lt.0.d0) then
            if (iz.eq.1) then
               div = zsolar*(xini-alpha*z)
               FeH10 = xsol*z/div
               FeH = log10(FeH10)
            elseif (iz.eq.2) then
               FeH10 = 10.d0**FeH
               div = alpha*FeH10+xsol/zsolar
               z = xini*FeH10/div
            endif
         else            
            if (iz.eq.1) then
               xFe = xHeavy*fracFe 
               x = (xsol-xini)*xFe/xFesol+xini
               FeH = log10(xFe/(56.d0*x))-FeHsol
            elseif (iz.eq.2) then
               FeH10 = 56.d0*10.d0**(FeH+FeHsol)
               div = xFesol-(xsol-xini)*FeH10
               z = zsolar*xini*FeH10/div
            endif
         endif
      endif
      if (iiz.eq.0) return

*** Write output value:

      if (iiz.eq.1) then
         write (*,2100) z
      endif
      if (iiz.eq.2) then
         write (*,2200) FeH
      endif

 999  format (1x,'*** Unacceptable choice for iz !')

 1000 format (/,5x,'Calculate Z from a given [Fe/H]  [1]',
     &     /,5x,'Calculate [Fe/H] from a given Z  [2]')
 1100 format (i1)
 1200 format (1x,'Input [Fe/H] value = ?')
 1300 format (d13.6)
 1400 format (1x,'Input Z value = ?')
 1500 format (d12.6)

c 2100 format (5x,'Z = ',1pe12.6)
 2100 format (5x,'Z = ',f8.6)
 2200 format (5x,'[Fe/H] = ',1pe13.6)

      return
      end


************************************************************************

      SUBROUTINE PSITOTHETA

************************************************************************
*   convert angular variable PSI --> THETA
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.rot'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.transp'
      include 'evolcom.therm'

      integer i,error,k,l,iret

      double precision system
      double precision rklth,tklth,opaclth
      double precision phim,deltam
      double precision vom,vrray,vxpsi
      double precision abmurj,abmuj
      double precision dift,dife,Dhold
      double precision r,vr

      common /lthopac/ rklth(19),tklth(85),opaclth(19,85,120)
      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)
      common /gradmuj/ abmurj(nsh),abmuj(nsh)
      common /difcirc/ dift(nsh),dife(nsh)
      common /calcDh/ Dhold(nsh)
      common /varrvr/ r(nsh),vr(nsh)


      error = 0

      print *,' Input : starevol.par and modini.bin (V >= 2.80)'
      open (unit = 99,err = 20, file = 'starevol.par',status = 'old')
      open (unit = 90,file = '_tmpdisp',status = 'unknown')
      open (unit = 91,file = 'modini1.bin',form = 'unformatted',
     &      status = 'unknown')
      open (unit = 92,file = '_tmp.bin',form = 'unformatted',
     &      status = 'unknown')

************************************************************************
*  read and rearrangements of all the data tables + read of parameters *
************************************************************************

      call rmodpar
         
************************************************************************
*     read of the initial model                                        *
************************************************************************

      call rinimod (91,93,92,94,0)
      iret=system("rm -f _tmp.bin" )
      iret=system("rm -f _tmpdisp.bin" )

      if (bin_version.lt.2.80d0) then
         write (*,*) ' Binary file already in THETA - no change made'
         close (91)
         close (93)
         close (92)
         close (95)
         return
      endif

      do l = 1,nsp
         do k = 1,nmod
            ysp(k,l) = xsp(k,l)/anuc(l)
         enddo
      enddo
      call eos (1,error)

      turnover(1) = 0.d0

      if (numeric.eq.2.or.numeric.eq.3) then
c.. first order accuracy in the spatial derivatives
         do i = 1,nmod
            wi(i) = 0.5d0
            wj(i) = 0.5d0
         enddo
      else
c.. second order accuracy in the spatial deriatives
         wi(1) = 0.5d0
         wj(1) = 0.5d0
         do i = 2,nmod1
            wi(i) = dm(i)/(dm(i)+dm(i-1))
            wj(i) = 1.d0-wi(i)
         enddo
         wi(nmod) = 0.5d0
         wj(nmod) = 0.5d0
      endif
      do k = 2,nmod
         deltam = deltaKS(k-1)*wi(k)+deltaKS(k)*wj(k)
         phim = phiKS(k-1)*wi(k)+phiKS(k)*wj(k)
         xpsis(k) = (phim*xlambdas(k)-deltam*xpsis(k))/vomega(k)*
     &        vomega(nmod)
         vxpsi(k) = (phim*xlambdas(k)-deltam*vxpsi(k))/vomega(k)*
     &        vomega(nmod)
         omega(k) = vomega(k)
      enddo
      omega(1) = vomega(1)
      xpsis(1) = xpsis(2)
      vxpsi(1) = vxpsi(2)

      write (95) ndtold,ndbold,ur_Ss,xmom_tots,inits,
     &     (auxs(k),urs(k),vxpsi(k),xpsis(k),V_circ(k),
     &     xalpha(k),xlambdas(k),xKt(k),dift(k),dife(k),
     &     Dhold(k),Dconv(k),abmurj(k),xnuvv(k),xnum(k),k = 1,
     &     nmod),((vxsp(k,l),l = 1,nis-1),k = 1,nmod)

      close (90)  !  evoldisp
      close (91)  !  modini.bin
      close (92)  !  nextini.bin
      if (irotbin.eq.1) then
         close (93)  !  modang.bin
         close (95)  !  nextang.bin
      endif

      print *,'  input : modini.bin, modang.bin --->  output : ',
     &     'nextang.bin,'
      return
 20   print *,'error in opening starevol.par, may be not found'

      return
      end



************************************************************************

      SUBROUTINE THERMO (error)

************************************************************************

      implicit none

      integer error

      return
      end

************************************************************************

      SUBROUTINE FSCREEN (t,ro,xsp,rhpsi,anuc,znuc)

************************************************************************

      implicit none

      include 'evolpar.star'

      double precision y,t,ro,xsp,rhpsi,anuc,znuc

      dimension y(nsp),t(nsh),ro(nsh),xsp(nsh,nsp),rhpsi(nsh),anuc(nsp),
     &     znuc(nsp)

      return
      end


************************************************************************

      SUBROUTINE DTMS

************************************************************************
* This code  computes the typical maximum time-step suited to          *
* compute the main sequence phase of a star of initial mass M and      *
* of initial metallicity Z in roughly 200 evolution models;            *
* this corresponds to a maximum time-step of 4.61d7 for a M = 1.00 and *
* Z = 0.020 star (dtref)                                               *
*                                                                      *
* The assumed formula to perform this computation is                   *
* dt(M,Z) = dt(1,0.02) * (M/1)**-2.9 * (Z/0.02)**0.1 , where           *
*    M is given in solar units                                         *
*                                                                      *
* This program needs, as input: M and Z                                *
* The only output is: dt(M,Z) to put as dtmax, in starevol.par         *
************************************************************************

      implicit none

      double precision m,z
      double precision alpha,beta,dtref
      double precision dtmax

*** Assumed values:

      alpha = -2.0d0
      beta = 0.10d0
      dtref = 4.61d7

*** Input values:

      write (*,1000)
      read (*,*) m,z

*** Computations:

      dtmax = dtref*m**alpha*(z/0.020d0)**beta

*** Writing output value:

      write (6,2100) dtmax

 1000 format (/,5x,'Enter initial M (in Mo) and metallicity Z')

 2100 format (5x,'dtmax = ',1pe12.6,/)

      end


************************************************************************

      subroutine rmodrot(ureadrot)

************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.ion'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.rot'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.therm'
      include 'evolcom.teq'
      include 'evolcom.var'
      include 'evolcom.nuc'
      include 'evolcom.transp'

      integer j,k,l,ureadrot,nspr

      character*1 crzc(nsh)

      double precision y(nsp),xc11,xsi31,xsi32,xp32
      double precision znew,mini,rini,zini,pindex,FeH
      double precision dift,dife,Dhold
      double precision abmurj,abmuj
      double precision vxpsi(nsh),vxspr(nmod,nsp)

      common /param/ znew,mini,rini,zini,pindex,FeH
      common /difcirc/ dift(nsh),dife(nsh)
      common /calcDh/ Dhold(nsh)
      common /gradmuj/ abmurj(nsh),abmuj(nsh)

      read (ureadrot) ndtold,ndbold,ur_Ss,xmom_tots,inits,
     &        (auxs(k),urs(k),vxpsi(k),xpsis(k),V_circ(k),
     &        xalpha(k),xlambdas(k),xKt(k),dift(k),dife(k),
     &        Dhold(k),Dconv(k),abmurj(k),xnuvv(k),xnum(k),k = 1,
     &        nmod),((vxspr(k,l),l = 1,nspr),k = 1,nmod)

      end

************************************************************************

      subroutine rmod(uread)

************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.ion'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.rot'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.therm'
      include 'evolcom.teq'
      include 'evolcom.var'
      include 'evolcom.nuc'

      integer j,k,l,uread,nspr

      character*1 crzc(nsh)

      double precision y(nsp),xc11,xsi31,xsi32,xp32
      double precision znew,mini,rini,zini,pindex,FeH

      common /param/ znew,mini,rini,zini,pindex,FeH

      read (uread) model,nphase,totm,time,dtn,neff,bin_version,nmod,
     &     nspr,mdredgeup,FeH,turnenv,rayon0,
     &     (crzc(k),u(k),r(k),lnf(k),lnt(k),lum(k),
     &     vu(k),vr(k),vlnf(k),vlnt(k),vlum(k),vfacc(k),
     &     exposure(k),vmueinv(k),dm(k),vhconv(k),vomega(k),
     &     ro(k),vro(k),p(k),vp(k),e(k),ve(k),venucl(k),s(k),vs(k),
     &     vsconv(k),vpvisc(k),tau(k),
     &     (xsp(k,l),l = 1,nspr),k = 1,nmod),
     &     (mtotlos(j),j = 1,nspr)
      do k = 1,nmod
         xsp(k,nis) = 1.d-50
         xsp(k,nsp) = 1.d-50
         if (crzc(k).eq.'c') crz(k) = -2 
         if (crzc(k).eq.'a') crz(k) = -1 
         if (crzc(k).eq.'s') crz(k) = -3 
         if (crzc(k).eq.'S') crz(k) = 3 
         if (crzc(k).eq.'t') crz(k) = 2 
         if (crzc(k).eq.'o') crz(k) = 1 
         if (crzc(k).eq.'r') crz(k) = 4 
      enddo
      m(1) = 0.d0
      do k = 2,nmod
         m(k) = m(k-1)+dm(k-1)
      enddo

      mtotlos(nis) = 1.d-50
      mtotlos(nsp) = 1.d-50

***   make correspondance between old and new binary files
      if (bin_version.gt.5.d2.or.bin_version.eq.0.d0) then
         write (*,*) 'species conversion in rinimod'
         do k=1,nsp
            y(k) = xsp(1,k)
         enddo
         do k = 1,nmod
            xc11 = xsp(k,13)
            xsi31 = xsp(k,40)
            xsi32 = xsp(k,42)
            xp32 = xsp(k,43)
            xsp(k,45) = xsp(k,44)+xp32+xsi32
            xsp(k,44) = xsp(k,41)+xsi31
            xsp(k,43) = xsp(k,39)
            xsp(k,42) = xsp(k,38)
            xsp(k,41) = xsp(k,37)
            xsp(k,40) = xsp(k,36)
            xsp(k,39) = 1.d-50   ! Mg 27
            xsp(k,38) = xsp(k,35)
            xsp(k,37) = xsp(k,34)
            xsp(k,36) = xsp(k,33)
            xsp(k,35) = xsp(k,32)
            xsp(k,34) = 1.d-50   ! Na 25
            xsp(k,33) = xsp(k,31)
            xsp(k,32) = 1.d-50   ! Na 24
            xsp(k,31) = xsp(k,30)
            xsp(k,30) = 1.d-50   ! Ne 23
            do j = 13,24
               xsp(k,j) = xsp(k,j+1)
            enddo           
            xsp(k,25) = 1.d-50   ! F  20
            xsp(k,12) = xsp(k,12)+xc11
         enddo
         do k=1,nsp
            write (*,20) k,k,y(k),k,xsp(1,k)
 20         format(' new species index :',i2,', old xsp(1,',i2,') = ',
     &           1pe9.3,' --> new xsp(1,',i2,') = ',1pe9.3)
         enddo
         bin_version = 2.61d0
      endif

      end


************************************************************************

      subroutine wmod(uread)

************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.ion'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.rot'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.therm'
      include 'evolcom.teq'
      include 'evolcom.var'

      integer j,k,l,uread
      double precision znew,mini,rini,zini,pindex,FeH

      character*1 crzc(nsh)

      common /param/ znew,mini,rini,zini,pindex,FeH

      do k = 1,nmod
         if (crz(k).eq.-2) crzc(k) = 'c'
         if (crz(k).eq.-1) crzc(k) = 'a'
         if (crz(k).eq.-3) crzc(k) = 's'
         if (crz(k).eq.3) crzc(k) = 'S'
         if (crz(k).eq.2) crzc(k) = 't'
         if (crz(k).eq.1) crzc(k) = 'o'
         if (crz(k).eq.4) crzc(k) = 'r'
      enddo

      write (uread) model,nphase,totm,time,dtn,neff,bin_version,nmod,
     &     nis-1,mdredgeup,FeH,turnenv,rayon0,
     &     (crzc(k),u(k),r(k),lnf(k),lnt(k),lum(k),
     &     vu(k),vr(k),vlnf(k),vlnt(k),vlum(k),vfacc(k),
     &     m(k),vmueinv(k),dm(k),vhconv(k),vomega(k),
     &     ro(k),vro(k),p(k),vp(k),e(k),ve(k),venucl(k),s(k),vs(k),
     &     vsconv(k),vpvisc(k),tau(k),
     &     (xsp(k,l),l = 1,nis-1),k = 1,nmod),
     &     (mtotlos(j),j = 1,nis-1)

      end



************************************************************************

      SUBROUTINE mypause()

************************************************************************

      implicit none

      character a*1

      print *
      print *, " Press <c> to continue"
      read (*,*) a
      if (a.eq.'x'.or.a.eq.'q') stop

      END


************************************************************************

      SUBROUTINE SORT1 (xxor,filename,filename2,nfile)

************************************************************************
* Sort array filename into ascending numerical order (output filename2)*
* following xxor ranking                                               *
************************************************************************

      implicit none

      character filename*100,filename2*100,rfilename*100

      integer i,j,l,nfile

      double precision xxor
      double precision xxori

      dimension xxor(nfile),filename(nfile),filename2(nfile)

      do l = 1,nfile
         filename2(l) = filename(l)
      enddo

      do j = 2,nfile
         xxori = xxor(j)
         rfilename = filename2(j)
         do i = j-1,1,-1
            if (xxor(i).le.xxori) goto 10
            xxor(i+1) = xxor(i)
            filename2(i+1) = filename2(i)
         enddo
         i = 0
 10      xxor(i+1) = xxori
         filename2(i+1) = rfilename
      enddo

      return
      end


************************************************************************

      SUBROUTINE fracmass(hxsp,xFenew)
 
************************************************************************       
* On utilise pour le calcul la formule suivante:                       *
* [Z/Fe]=log10(hxsp(Z)*mmol(Fe)/hxsp(Fe)*mmol(Z))-log10(hxsp(Z)*       *
* mmol(Fe)/hxsp(Fe)*mmol(Z))[Soleil]                                   *
************************************************************************

      
      implicit none
      
      include 'evolpar.star'

      include 'evolcom.data'
      
      integer refsolchem
      double precision xFenew,fac,hxsp(nsp)
      double precision OFe, MgFe, AlFe, NaFe      
      

c..   Iron mass fraction according to adopted reference solar composition
      if (refsolar.eq.'GN93') then 
         fac = 1.147417671020920d-3
         refsolchem = 1
      endif
      if (refsolar.eq.'AGS05') then 
         fac = 1.0666D-03
         refsolchem = 2
      endif
      if (refsolar.eq.'AGSS09') then 
         fac = 1.1953D-03
         refsolchem = 3
      endif
      if (refsolar.eq.'GRID') then 
         fac = 1.147417671020920d-3
         refsolchem = 4
      endif

      OFe = log10(hxsp(20)*56.d0/(xFenew*16.d0))-
     &      log10(xspref(refsolchem,20)*56.d0/(fac*16.d0))
      
      
      MgFe = log10(hxsp(33)*56.d0/(xFenew*24.d0))-
     &       log10(xspref(refsolchem,33)*56.d0/(fac*24.d0))
      
      
      AlFe = log10(hxsp(40)*56.d0/(xFenew*27.d0))-
     &       log10(xspref(refsolchem,40)*56.d0/(fac*27.d0))
      
      
      NaFe = log10(hxsp(31)*56.d0/(xFenew*23.d0))-
     &       log10(xspref(refsolchem,31)*56.d0/(fac*23.d0))
      
      
      write (*,500) OFe, NaFe, MgFe, AlFe
  
 500  format (3x,'[O/Fe] = ',1pe11.4,', [Na/Fe] = ',1pe11.4,
     &  ', [Mg/Fe] = ',1pe11.4,', [Al/Fe] = ',1pe11.4)   
      
      return
      end

************************************************************************

      subroutine SOLX (initsolar,zinisol)

************************************************************************
* Computation of solar isotopic mass fractions from abundances in      *
* the log(N_H) = 12.00 scale.                                          *
* Based on meteoritic solar chemical composition from                  *
* Grevesse & Sauval, 1998, Space Sci. Rev. 85, 161                     *
*                                                                      *
* Conversion formulae from Arnett "Supernovae and Nucleosynthesis" book*
*                                                                      *
* epsilon = log(N_X/N_H)+12 withh N_H = 1d12 so that N_X = 10**epsilon *
* where N_X is the number density of element X, all isotopes included  *
* The number density for each isotope j of element X is then           *
* N_j = N_X*fraciso(j)                                                 *
* The mass fraction of isotope Xj is equal to (N_j*A_j)/(sum_k N_k*A_k)*
* where k ranges over all the stable isotopes considered, and          *
* A_j is the atomic weight of isotope j.                               *
*                                                                      *
* One can choose the reference for the chemical solar composition:     *
*      initsolar = 'AG89' => Anders & Grevesse 1989                    *
*      initsolar = 'GNS96' => Grevesse, Noels & Sauval 1996            *
*      initsolar = 'GS98' => Grevesse & Sauval 1998                    *
*      initsolar = 'AGS05' => Asplund, Grevesse & Sauval, 2005         * 
*      initsolar = 'AGSS09' => Asplund, Grevesse, Sauval & Scott 2009  * 
************************************************************************

      implicit none 
      integer nele,niso,nam,nnetwa,nsp,idum,nis,ni,nelt,nisotop,
     &     neltaccount,iheavy
      integer i,j,l,ij
      integer refsolchem

      parameter (nele = 83, niso = 286)
      parameter (nam = 240)
      parameter (nnetwa = 35, nsp = 55,nis = 54)

      character*2 name
      character*3 adum
      character*6 initsolar,refsolar
      character*100 utils
      logical noaccounteos

      double precision element,am,z,aiso,fracx,abund,abiso,xspref,zsol
      double precision dum,relah1,rela,sum,sumtot,fracxz,sumtotal,
     &     zinisol,xspsol
      double precision xini,yini,zini,xsol,ysol,zsolar,xheavysol
      double precision xFesol,FeHsol,fracFe
      double precision sumeos,sumnoeos,YNi,YFe,YAr,YTi,YCa,YCr,YMn

      dimension element(nele),fracx(niso)
      common /abund/ abund(nele)
      common /atomic/ am(niso),z(nele),ni(nele),name(nele)
      common /isorat/ abiso(niso),relah1(niso),rela(niso)

      common /isoabund/ aiso(niso)
      common /composition/ xspref(4,nis+6)
      common /solar/ xspsol(nis+6),zsol,refsolar
      common /bigbang/xini,yini,xsol,ysol,zsolar,xFesol,
     &     FeHsol,fracFe,iheavy

C       write(*,*) ' Set of chemical solar abundances?'
C       write(*,*) '(AG89, GN93,GNS96, GS98 or AGS05)'
C       read(5,'(a5)')initsolar
      print *,initsolar,zinisol

      call getenv('DIR_PROG',utils)

      if (trim(initsolar).eq.'GN93') then 
         nelt = 28
         nisotop = 70
         neltaccount = nnetwa
      else
         nelt = nele
         nisotop = niso
         neltaccount = nnetwa
      endif


      if (trim(initsolar).eq.'GN93') refsolchem = 1
      if (trim(initsolar).eq.'AGS05') refsolchem = 2
      if (trim(initsolar).eq.'AGSS09') refsolchem = 3
      if (trim(initsolar).eq.'GRID') refsolchem = 4

      print *,trim(initsolar)
      open (unit = 10,file = trim(utils) // trim(initsolar) //'.dat')

      open (unit = 21,file = 'solar_data') 
      

      if (trim(initsolar).eq.'AG89'.or.trim(initsolar).eq.'AGS05'.or.
     &     trim(initsolar).eq.'AGSS09') read (10,1005) 
      do l = 1,nelt
         if (trim(initsolar).ne.'GNS96') then
            read (10,1000) idum,adum,abund(l),dum
         else
            read (10,1100) idum,adum,dum,abund(l)
         endif
      enddo

      do i = 1,niso
         aiso(i) = 0.d0
         fracx(i) = 0.d0
      enddo
      do i = 1,nis+6
         xspref(refsolchem,i) = 0.d0
      enddo

      do i = 1,nelt
         element(i) = 10.d0**(abund(i)-12.d0)*1.d12
      enddo

      sumtot = 0.d0
      ij = 0
      do i = 1,nelt
         do j = 1,ni(i)
            ij = ij+1
            aiso(ij) = element(i)*abiso(ij)
            sumtot = sumtot+aiso(ij)*am(ij)
         enddo
      enddo

      sum = 0.d0
      do i = 1,nisotop
         fracx(i) = aiso(i)*am(i)/sumtot
         sum = sum+fracx(i)
      enddo
      write (*,111) trim(initsolar),fracx(62)
 111  format (/,'  XFesol for ',a5,' solar chemical composition : ',
     &     1pd10.4)
      xFesol = fracx(62)


* Computation of the mass fraction associated with the 'heavy' element
      fracxz = 0.d0
      do j = neltaccount+1,nisotop
         fracxz = fracxz+fracx(j)
      enddo
c$$$* Determine the contribution of Ar, Ca, Ti, Cr, Mn, Fe and Ni to heavy
c$$$         noaccounteos = ((j.ge.39.and.j.le.41).or.j.eq.48.or.j.eq.54.or.
c$$$     &     j.eq.55.or.j.eq.65)
c$$$         if (noaccounteos) sumnoeos = sumnoeos+fracx(j)
c$$$      enddo
c$$$      YAr = (fracx(36)+fracx(37)+fracx(38))/fracxz
c$$$      YCa = (fracx(42)+fracx(43)+fracx(44)+fracx(45)+fracx(46)+fracx(47)
c$$$     &     )/fracxz
c$$$      YTi = (fracx(49)+fracx(50)+fracx(51)+fracx(52)+fracx(53))/fracxz
c$$$      YCr = (fracx(56)+fracx(57)+fracx(58)+fracx(59))/fracxz
c$$$      YMn = fracx(60)/fracxz
c$$$      YFe = (fracx(61)+fracx(62)+fracx(63)+fracx(64))/fracxz
c$$$      YNi = (fracx(66)+fracx(67)+fracx(68)+fracx(69)+fracx(70))/fracxz
c$$$      print *,'sum heavy not included in eos',YAr,YCa,YTi,YCr,YMn,YFe,
c$$$     &     YNi,sumnoeos/fracxz
c$$$      
c$$$      YAr = YAr/(1.d0-sumnoeos/fracxz)
c$$$      YCa = YCa/(1.d0-sumnoeos/fracxz)
c$$$      YTi = YTi/(1.d0-sumnoeos/fracxz)
c$$$      YCr = YCr/(1.d0-sumnoeos/fracxz)
c$$$      YMn = YMn/(1.d0-sumnoeos/fracxz)
c$$$      YFe = YFe/(1.d0-sumnoeos/fracxz)
c$$$      YNi = YNi/(1.d0-sumnoeos/fracxz)
c$$$
c$$$      YFe = YFe + (1.d0-(YAr+YCa+YCr+YMn+YFe+YNi+YTi))
c$$$      print *,'sum heavy not included in eos',YAr,YCa,YTi,YCr,YMn,YFe
c$$$     &   ,YNi,YAr+YCa+YCr+YMn+YFe+YNi+YTi , 1.d0-sumnoeos/fracxz 
c$$$
c$$$      stop


      sum = 0.d0
      do j = 1,neltaccount
         sum = sum+fracx(j)
      enddo

      sumtotal = sum+fracxz
      print *,sumtotal
      
      if (abs(sumtotal-1.d0).gt.1.d-10) then
         print *,'renomalization : ',1.d0-fracxz-sum
* Renormalization of the mass fraction sum on hydrogen
         fracx(1) = fracx(1)+(1.d0-fracxz-sum)
      endif

* Checking that the sum of solar mass fractions is equal to unity
      sum = 0.d0
      do i = 1,neltaccount
         sum = sum + fracx(i)
      enddo
      sum = sum+fracxz
      if (abs(sum-1.d0).gt.1.d-10) then
         print *,' SOLAR ABUNDANCE NOT NORMALIZED !! - stop solx'
         stop
      endif

      xspref(refsolchem,1) = 1.d-50         !   1 NEUT  
      xspref(refsolchem,2) = fracx(1)       !   2 PROT  
      xspref(refsolchem,3) = fracx(2)       !   3 DEUT  
      xspref(refsolchem,4) = fracx(3)       !   4 HE  3 
      xspref(refsolchem,5) = fracx(4)       !   5 HE  4 
      xspref(refsolchem,6) = fracx(5)       !   6 LI  6 
      xspref(refsolchem,7) = fracx(6)       !   7 LI  7 
      xspref(refsolchem,8) = 1.d-50         !   8 BE  7 *
      xspref(refsolchem,9) = 1.d-50         !   9 B   8 *
      xspref(refsolchem,10) = fracx(7)      !  10 BE  9  
      xspref(refsolchem,11) = fracx(8)      !  11 B  10  
      xspref(refsolchem,12) = fracx(9)      !  12 B  11  
      xspref(refsolchem,13) = fracx(10)     !  13 C  12  
      xspref(refsolchem,14) = fracx(11)     !  14 C  13 
      xspref(refsolchem,15) = 1.d-50        !  15 N  13 *
      xspref(refsolchem,16) = 1.d-50        !  16 C  14 *
      xspref(refsolchem,17) = fracx(12)     !  17 N  14 
      xspref(refsolchem,18) = fracx(13)     !  18 N  15 
      xspref(refsolchem,19) = 1.d-50        !  19 O  15 *
      xspref(refsolchem,20) = fracx(14)     !  20 O  16 
      xspref(refsolchem,21) = fracx(15)     !  21 O  17 
      xspref(refsolchem,22) = fracx(16)     !  22 O  18 
      xspref(refsolchem,23) = 1.d-50        !  23 F  18 *
      xspref(refsolchem,24) = fracx(17)     !  24 F  19 
      xspref(refsolchem,25) = 1.d-50        !  25 F  20 *
      xspref(refsolchem,26) = fracx(18)     !  26 NE 20 
      xspref(refsolchem,27) = fracx(19)     !  27 NE 21 
      xspref(refsolchem,28) = fracx(20)     !  28 NE 22 
      xspref(refsolchem,29) = 1.d-50        !  29 NA 22 *
      xspref(refsolchem,30) = 1.d-50        !  30 NE 23 *
      xspref(refsolchem,31) = fracx(21)     !  31 NA 23 
      xspref(refsolchem,32) = 1.d-50        !  32 NA 24 *
      xspref(refsolchem,33) = fracx(22)     !  33 MG 24 
      xspref(refsolchem,34) = 1.d-50        !  34 NA 25 *
      xspref(refsolchem,35) = fracx(23)     !  35 MG 25 
      xspref(refsolchem,36) = fracx(24)     !  36 MG 26 
      xspref(refsolchem,37) = 1.d-50        !  37 AL26m *
      xspref(refsolchem,38) = 1.d-50        !  38 AL26g *
      xspref(refsolchem,39) = 1.d-50        !  39 MG 27 *
      xspref(refsolchem,40) = fracx(25)     !  40 AL 27 
      xspref(refsolchem,41) = fracx(26)     !  41 SI 28 
      xspref(refsolchem,42) = fracx(27)     !  42 SI 29 
      xspref(refsolchem,43) = fracx(28)     !  43 SI 30 
      xspref(refsolchem,44) = fracx(29)     !  44 P  31 
      xspref(refsolchem,45) = fracx(30)     !  45 S  32 
      xspref(refsolchem,46) = fracx(31)     !  46 S  33 
      xspref(refsolchem,47) = fracx(32)     !  47 S  34 
      xspref(refsolchem,48) = 1.d-50        !  48 S  35 *
      xspref(refsolchem,49) = fracx(33)     !  49 CL 35 
      xspref(refsolchem,50) = fracx(34)     !  50 S  36 
      xspref(refsolchem,51) = 1.d-50        !  51 CL 36 *
      xspref(refsolchem,52) = fracx(35)     !  52 CL 37 
      xspref(refsolchem,53) = fracxz        !  53 HEAVY 
      if (refsolchem.eq.1) then
c.. Grevesse & Noels 93
c..    heavy = 83.103% Fe + 5.652% Ar + 4.911% Ni +
c..   4.319% Ca + 1.144% Cr + 0.635% Mn + 0.236% Ti
         xspref(refsolchem,54) = 5.652d-2*fracxz  ! Ar
         xspref(refsolchem,55) = 4.319d-2*fracxz    ! Ca
         xspref(refsolchem,56) = 0.236d-2*fracxz    ! Ti
         xspref(refsolchem,57) = 1.144d-2*fracxz    ! Cr 
         xspref(refsolchem,58) = 0.635d-2*fracxz    ! Mn
         xspref(refsolchem,59) = 8.3103d-1*fracxz   ! Fe 
         xspref(refsolchem,60) = 4.911d-2*fracxz   ! Ni
      else if (refsolchem.eq.2) then
c..   Asplund,Grevesse & Sauval 2005
c..   heavy = 85.579% Fe + 2.976% Ar + 4.927% Ni +
c..   4.235% Ca + 1.202% Cr + 0.879% Mn + 0.202% Ti 
         xspref(refsolchem,54) = 2.976d-2*fracxz  ! Ar
         xspref(refsolchem,55) = 4.235d-2*fracxz    ! Ca
         xspref(refsolchem,56) = 0.202d-2*fracxz    ! Ti
         xspref(refsolchem,57) = 1.202d-2*fracxz    ! Cr 
         xspref(refsolchem,58) = 0.879d-2*fracxz    ! Mn
         xspref(refsolchem,59) = 8.5579d-1*fracxz    ! Fe 
         xspref(refsolchem,60) = 4.927d-2*fracxz   ! Ni
      else if (refsolchem.eq.3) then
c..   Asplund, Grevesse, Sauval & Scott 2009
c..   heavy = 84.143% Fe + 5.652% Ar + 4.664% Ni +
c..   4.197% Ca + 1.087% Cr + 0.708% Mn + 0.204% Ti 
         xspref(refsolchem,54) = 5.652d-2*fracxz  ! Ar
         xspref(refsolchem,55) = 4.197d-2*fracxz    ! Ca
         xspref(refsolchem,56) = 0.204d-2*fracxz    ! Ti
         xspref(refsolchem,57) = 1.087d-2*fracxz    ! Cr 
         xspref(refsolchem,58) = 0.708d-2*fracxz    ! Mn
         xspref(refsolchem,59) = 8.4143d-1*fracxz    ! Fe 
         xspref(refsolchem,60) = 4.664d-2*fracxz   ! Ni
      endif
     
      do l = 1,nis
         write(667,*) xspref(refsolchem,l)
      enddo

      write (21,3000) (xspref(refsolchem,l),l = 1,2)
      write (21,3100) (xspref(refsolchem,l),l = 3,7)
      write (21,3100) (xspref(refsolchem,l),l = 8,12)
      write (21,3100) (xspref(refsolchem,l),l = 13,17)
      write (21,3100) (xspref(refsolchem,l),l = 18,22)
      write (21,3100) (xspref(refsolchem,l),l = 23,27)
      write (21,3100) (xspref(refsolchem,l),l = 28,32)
      write (21,3100) (xspref(refsolchem,l),l = 33,37)
      write (21,3100) (xspref(refsolchem,l),l = 38,42)
      write (21,3100) (xspref(refsolchem,l),l = 43,47)
      write (21,3100) (xspref(refsolchem,l),l = 48,52)
      write (21,3100) (xspref(refsolchem,l),l = 53,57)
      write (21,3200) (xspref(refsolchem,l),l = 58,nis+6)      
      zinisol = 1.d0-xspref(refsolchem,2)-xspref(refsolchem,3)
     & 	-xspref(refsolchem,4)-xspref(refsolchem,5)
      zsolar = zinisol
      xsol = xspref(refsolchem,2)
      ysol = xspref(refsolchem,5)
      xheavysol = fracxz
      fracFe = xFesol/xheavysol
      FeHsol = log10(xFesol/(56.d0*xsol))

      write(*,*) ' Output is ', trim(utils) // 'solar_data'

 1000 format (i2,1x,a2,2x,f5.2,2x,f11.7)
 1005 format (4(/))
 1100 format (i2,1x,a3,1x,f11.6,1x,f5.2)

 3000 format (/,6x,'data (xspref(i),i = 1,nis+6) /',1x,1pe6.0,',',
     &     0pe21.15,',')
 3100 format (5x,'&',1x,5(1pd11.5,','))
 3200 format (5x,'&',1x,2(1pd11.5,','),1pd11.5,' /')

      return
      end


************************************************************************
      subroutine CHANGEXSP(iu1,iu2)
************************************************************************
* Change composition to Asplund-Cunha for GVA GRID
************************************************************************
      implicit none

      
      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.data'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.spec'
      include 'evolcom.teq'

      integer iu1,iu2
      integer j,k,refsolchem
      
      double precision sum,heavy
*_____________________________
***   initialize species index
*-----------------------------

      ih1 =  2
      ih2 =  3
      ihe3 = 4
      ihe4 = 5
      ili6 = 6
      ili7 = 7
      ibe7 = 8
      ibe9 = 10
      ib8 =  9
      ib11 = 12
      ic12 = 13
      ic13 = 14
      ic14 = 16
      in13 = 14
      in14 = 17
      in15 = 18
      io15 = 19
      io16 = 20
      io17 = 21
      io18 = 22
      if18 = 23
      if19 = 24
      if20 = 25
      ine20 = 26
      ine23 = 30
      ina22 = 29
      ina23 = 31
      ina24 = 32
      ina25 = 34
      img24 = 33
      img25 = 35
      img27 = 39
      ial27 = 40
      isi28 = 41
      ip31 = 44
      is32 = 45
      is35 = 48
      icl35 = 49
      icl36 = 51


      rewind (iu1)
      call rmod(iu1)

      do k = 1,nmod
         do j = 1,nsp
            xsp(k,j) = 1.d-50
         enddo
      enddo

      do k = 1,nmod
         xsp(k,ih1)      = 7.199510383353603D-01
         xsp(k,ih2)      = 4.896166464000003D-05
         xsp(k,ihe3)     = 4.309361250122357D-05
         xsp(k,ihe4)     = 2.660000000000000D-01
         xsp(k,ili6)     = 6.348070062382616D-10
         xsp(k,ili7)     = 8.996284335074931D-09
         xsp(k,ibe9)     = 1.687361752301596D-10
         xsp(k,11)       = 8.090584217540142D-10 !B10 
         xsp(k,ib11)     = 3.941942292861476D-09
         xsp(k,ic12)     = 2.265305429944870D-03
         xsp(k,ic13)     = 3.631168314692027D-05
         xsp(k,in14)     = 6.562993361397868D-04
         xsp(k,in15)     = 2.341829587990781D-06
         xsp(k,io16)     = 5.690639357507690D-03
         xsp(k,io17)     = 3.820235567086054D-06
         xsp(k,io18)     = 1.284118370119140D-05
         xsp(k,if19)     = 5.383193113495001D-07
         xsp(k,ine20)    = 1.786658475578834D-03
         xsp(k,27)       = 5.698737469568906D-06 !Ne21
         xsp(k,28)       = 2.396650381755236D-04 !Ne22
         xsp(k,ina23)    = 2.653734474502985D-05
         xsp(k,img24)    = 4.991580001943164D-04
         xsp(k,35)       = 6.693082946261805D-05 !Mg25
         xsp(k,36)       = 7.674648424542478D-05 !Mg26
         xsp(k,ial27)    = 4.936163560320215D-05
         xsp(k,isi28)    = 5.973604724455358D-04
         xsp(k,42)       = 6.651698767600542D-05 !Si29
         xsp(k,43)       = 4.814167611933788D-05 !Si30
         xsp(k,ip31)     = 5.537539284595751D-06
         xsp(k,is32)     = 3.240687095754634D-04
         xsp(k,46)       = 3.480864395059802D-06 !S33
         xsp(k,47)       = 1.804006047947540D-05 !S34
         xsp(k,icl35)    = 6.540577506206110D-06
         xsp(k,52)       = 2.148960358406137D-06 !Cl37
      enddo
      sum = 0.d0
      do j = 1,nsp
         sum = sum+xsp(1,j)
      enddo
      print *,sum,1-sum
      heavy = 1.d0-sum
      do k = 1,nmod 
         xsp(k,53) = heavy
      enddo

      sum = 0.d0
      do j = 6,nsp
         sum = sum+xsp(1,j)
      enddo
      print *,'New Z',sum,1-xsp(1,2)-xsp(1,3)-xsp(1,4)-xsp(1,5),
     &     sum-1.4d-2,xsp(1,2)+xsp(1,3)

         
      if (refsolar.eq.'GN93') refsolchem = 1
      if (refsolar.eq.'AGS05') refsolchem = 2
      if (refsolar.eq.'AGSS09') refsolchem = 3
      if (refsolar.eq.'GRID') refsolchem = 4

      rewind (iu2)
      call wmod(iu2)

      close (iu2)
      close (iu1)

      open (unit = 31,file = 'solar.out',status = 'unknown')
      write (31,9000) bin_version,totm,time/sec,sum
      write (31,*)'                           Xsolnew       Xsolold'
      write (31,9001) xsp(1,1),xspref(refsolchem,1) 
      write (31,9002) xsp(1,2),xspref(refsolchem,2)
      write (31,9003) xsp(1,3),xspref(refsolchem,3)
      write (31,9004) xsp(1,4),xspref(refsolchem,4)
      write (31,9005) xsp(1,5),xspref(refsolchem,5)
      write (31,9006) xsp(1,6),xspref(refsolchem,6)
      write (31,9007) xsp(1,7),xspref(refsolchem,7)
      write (31,9008) xsp(1,8),xspref(refsolchem,8)
      write (31,9009) xsp(1,9),xspref(refsolchem,9)
      write (31,9010) xsp(1,10),xspref(refsolchem,10)
      write (31,9011) xsp(1,11),xspref(refsolchem,11)
      write (31,9012) xsp(1,12),xspref(refsolchem,12)
      write (31,9013) xsp(1,13),xspref(refsolchem,13)
      write (31,9014) xsp(1,14),xspref(refsolchem,14)
      write (31,9015) xsp(1,15),xspref(refsolchem,15)
      write (31,9016) xsp(1,16),xspref(refsolchem,16)
      write (31,9017) xsp(1,17),xspref(refsolchem,17)
      write (31,9018) xsp(1,18),xspref(refsolchem,18)
      write (31,9019) xsp(1,19),xspref(refsolchem,19)
      write (31,9020) xsp(1,20),xspref(refsolchem,20)
      write (31,9021) xsp(1,21),xspref(refsolchem,21)
      write (31,9022) xsp(1,22),xspref(refsolchem,22)
      write (31,9023) xsp(1,23),xspref(refsolchem,23)
      write (31,9024) xsp(1,24),xspref(refsolchem,24)
      write (31,9025) xsp(1,25),xspref(refsolchem,25)
      write (31,9026) xsp(1,26),xspref(refsolchem,26)
      write (31,9027) xsp(1,27),xspref(refsolchem,27)
      write (31,9028) xsp(1,28),xspref(refsolchem,28)
      write (31,9029) xsp(1,29),xspref(refsolchem,29)
      write (31,9030) xsp(1,30),xspref(refsolchem,30)
      write (31,9031) xsp(1,31),xspref(refsolchem,31)
      write (31,9032) xsp(1,32),xspref(refsolchem,32)
      write (31,9033) xsp(1,33),xspref(refsolchem,33)
      write (31,9034) xsp(1,34),xspref(refsolchem,34)
      write (31,9035) xsp(1,35),xspref(refsolchem,35)
      write (31,9036) xsp(1,36),xspref(refsolchem,36)
      write (31,9037) xsp(1,37),xspref(refsolchem,37)
      write (31,9038) xsp(1,38),xspref(refsolchem,38)
      write (31,9039) xsp(1,39),xspref(refsolchem,39)
      write (31,9040) xsp(1,40),xspref(refsolchem,40)
      write (31,9041) xsp(1,41),xspref(refsolchem,41)
      write (31,9042) xsp(1,42),xspref(refsolchem,42)
      write (31,9043) xsp(1,43),xspref(refsolchem,43)
      write (31,9044) xsp(1,44),xspref(refsolchem,44)
      write (31,9045) xsp(1,45),xspref(refsolchem,45)
      write (31,9046) xsp(1,46),xspref(refsolchem,46)
      write (31,9047) xsp(1,47),xspref(refsolchem,47)
      write (31,9048) xsp(1,48),xspref(refsolchem,48)
      write (31,9049) xsp(1,49),xspref(refsolchem,49)
      write (31,9050) xsp(1,50),xspref(refsolchem,50)
      write (31,9051) xsp(1,51),xspref(refsolchem,51)
      write (31,9052) xsp(1,52),xspref(refsolchem,52)
      write (31,9053) xsp(1,53),xspref(refsolchem,53)
      write (31,9054) xsp(1,54)
      write (31,9055) xsp(1,55)
      close (31)


 1000 format (3x,'computed Z = ',1pe13.7)
 9000 format (/,2x,'INITIAL CHEMICAL SPECIES - version ',f5.3,/,2x,
     &     'mass : ',f8.4,', age = ',1pe11.5,', Z = ',1pe12.6)
 9001 format ('   1 NEUT    1.00000   0.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9002 format ('   2 PROT    1.00000   1.',3x,1pe11.5,3x,1pe11.5)
 9003 format ('   3 DEUT    2.00000   1.',3x,1pe11.5,3x,1pe11.5)
 9004 format ('   4 HE  3   3.00000   2.',3x,1pe11.5,3x,1pe11.5)
 9005 format ('   5 HE  4   4.00000   2.',3x,1pe11.5,3x,1pe11.5)
 9006 format ('   6 LI  6   6.00000   3.',3x,1pe11.5,3x,1pe11.5)
 9007 format ('   7 LI  7   7.00000   3.',3x,1pe11.5,3x,1pe11.5)
 9008 format ('   8 BE  7   7.00000   4.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9009 format ('   9 B   8   8.00000   5.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9010 format ('  10 BE  9   9.00000   4.',3x,1pe11.5,3x,1pe11.5)
 9011 format ('  11 B  10  10.00000   5.',3x,1pe11.5,3x,1pe11.5)
 9012 format ('  12 B  11  11.00000   5.',3x,1pe11.5,3x,1pe11.5)
 9013 format ('  13 C  12  12.00000   6.',3x,1pe11.5,3x,1pe11.5)
 9014 format ('  14 C  13  13.00000   6.',3x,1pe11.5,3x,1pe11.5)
 9015 format ('  15 N  13  13.00000   7.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9016 format ('  16 C  14  14.00000   6.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9017 format ('  17 N  14  14.00000   7.',3x,1pe11.5,3x,1pe11.5)
 9018 format ('  18 N  15  15.00000   7.',3x,1pe11.5,3x,1pe11.5)
 9019 format ('  19 O  15  15.00000   8.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9020 format ('  20 O  16  16.00000   8.',3x,1pe11.5,3x,1pe11.5)
 9021 format ('  21 O  17  17.00000   8.',3x,1pe11.5,3x,1pe11.5)
 9022 format ('  22 O  18  18.00000   8.',3x,1pe11.5,3x,1pe11.5)
 9023 format ('  23 F  18  18.00000   9.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9024 format ('  24 F  19  19.00000   9.',3x,1pe11.5,3x,1pe11.5)
 9025 format ('  25 F  20  20.00000   9.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9026 format ('  26 NE 20  20.00000  10.',3x,1pe11.5,3x,1pe11.5)
 9027 format ('  27 NE 21  21.00000  10.',3x,1pe11.5,3x,1pe11.5)
 9028 format ('  28 NE 22  22.00000  10.',3x,1pe11.5,3x,1pe11.5)
 9029 format ('  29 NA 22  22.00000  11.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9030 format ('  30 NE 23  23.00000  10.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9031 format ('  31 NA 23  23.00000  11.',3x,1pe11.5,3x,1pe11.5)
 9032 format ('  32 NA 24  24.00000  11.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9033 format ('  33 MG 24  24.00000  12.',3x,1pe11.5,3x,1pe11.5)
 9034 format ('  34 NA 25  25.00000  11.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9035 format ('  35 MG 25  25.00000  12.',3x,1pe11.5,3x,1pe11.5)
 9036 format ('  36 MG 26  26.00000  12.',3x,1pe11.5,3x,1pe11.5)
 9037 format ('  37 AL26M  26.00000  13.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9038 format ('  38 AL26G  26.00000  13.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9039 format ('  39 MG 27  27.00000  12.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9040 format ('  40 AL 27  27.00000  13.',3x,1pe11.5,3x,1pe11.5)
 9041 format ('  41 SI 28  28.00000  14.',3x,1pe11.5,3x,1pe11.5)
 9042 format ('  42 SI 29  29.00000  14.',3x,1pe11.5,3x,1pe11.5)
 9043 format ('  43 SI 30  30.00000  14.',3x,1pe11.5,3x,1pe11.5)
 9044 format ('  44 P  31  31.00000  15.',3x,1pe11.5,3x,1pe11.5)
 9045 format ('  45 S  32  32.00000  16.',3x,1pe11.5,3x,1pe11.5)
 9046 format ('  46 S  33  33.00000  16.',3x,1pe11.5,3x,1pe11.5)
 9047 format ('  47 S  34  34.00000  16.',3x,1pe11.5,3x,1pe11.5)
 9048 format ('  48 S  35  35.00000  16.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9049 format ('  49 CL 35  35.00000  17.',3x,1pe11.5,3x,1pe11.5)
 9050 format ('  50 S  36  36.00000  16.',3x,1pe11.5,3x,1pe11.5)
 9051 format ('  51 CL 36  36.00000  17.',3x,1pe11.5,3x,1pe11.5,2x,'*')
 9052 format ('  52 CL 37  37.00000  17.',3x,1pe11.5,3x,1pe11.5)
 9053 format ('  53 HEAVY  38.00000  18.',3x,1pe11.5,3x,1pe11.5)
 9054 format ('  54 CAPTN  39.00000  18.',3x,1pe11.5,3x,
     &     '1.00000d-50  *')
 9055 format ('  55 OOOOO 999.99999  99.',3x,1pe11.5,3x,
     &     '1.00000d-50  *')

      



      return
      end

************************************************************************

      BLOCK DATA

************************************************************************

      implicit double precision (a-h , o-z)

      character*2 name

      common /atomic/ am(286),z(83),ni(83),name(83)
      common /isorat/ abiso(286),relah1(286),rela(286)

      data (name(k),k = 1,83) /
     & 'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al',
     & 'Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe',
     & 'Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ',
     & 'Zr','Nb','Mo','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ',
     & 'Xe','Cs','Ba','La','Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho',
     & 'Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     & 'Tl','Pb','Bi','Th','U ' /

      data (ni(k),k = 1,83) /
     & 2,2,2,1,2,2,2,3,1,3,1,3,1,3,1,4,2,3,3,6,1,5,2,4,1,4,1,5,2,5,2,5,
     & 1,6,2,6,2,4,1,5,1,7,7,1,6,2,8,2,10,2,8,1,9,1,7,2,4,1,7,7,2,7,1,
     & 7,1,6,1,7,2,6,2,5,2,7,2,6,1,7,2,4,1,1,2 /

      data (z(k),k = 1,83) /
     & 1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,7.d0,8.d0,9.d0,10.d0,11.d0,12.d0,
     & 13.d0,14.d0,15.d0,16.d0,17.d0,18.d0,19.d0,20.d0,21.d0,22.d0,
     & 23.d0,24.d0,25.d0,26.d0,27.d0,28.d0,29.d0,30.d0,31.d0,32.d0,
     & 33.d0,34.d0,35.d0,36.d0,37.d0,38.d0,39.d0,40.d0,41.d0,42.d0,
     & 44.d0,45.d0,46.d0,47.d0,48.d0,49.d0,50.d0,51.d0,52.d0,53.d0,
     & 54.d0,55.d0,56.d0,57.d0,58.d0,59.d0,60.d0,62.d0,63.d0,64.d0,
     & 65.d0,66.d0,67.d0,68.d0,69.d0,70.d0,71.d0,72.d0,73.d0,74.d0,
     & 75.d0,76.d0,77.d0,78.d0,79.d0,80.d0,81.d0,82.d0,83.d0,90.d0,
     & 92.d0 /

      data (am(l),l = 1,286) /
     & 1.d0,2.d0,3.d0,4.d0,6.d0,7.d0,9.d0,10.d0,11.d0,12.d0,13.d0,14.d0,
     & 15.d0,16.d0,17.d0,18.d0,19.d0,20.d0,21.d0,22.d0,23.d0,24.d0,
     & 25.d0,26.d0,27.d0,28.d0,29.d0,30.d0,31.d0,32.d0,33.d0,34.d0,
     & 36.d0,35.d0,37.d0,36.d0,38.d0,40.d0,39.d0,40.d0,41.d0,40.d0,
     & 42.d0,43.d0,44.d0,46.d0,48.d0,45.d0,46.d0,47.d0,48.d0,49.d0,
     & 50.d0,50.d0,51.d0,50.d0,52.d0,53.d0,54.d0,55.d0,54.d0,56.d0,
     & 57.d0,58.d0,59.d0,58.d0,60.d0,61.d0,62.d0,64.d0,63.d0,65.d0,
     & 64.d0,66.d0,67.d0,68.d0,70.d0,69.d0,71.d0,70.d0,72.d0,73.d0,
     & 74.d0,76.d0,75.d0,74.d0,76.d0,77.d0,78.d0,80.d0,82.d0,79.d0,
     & 81.d0,78.d0,80.d0,82.d0,83.d0,84.d0,86.d0,85.d0,87.d0,84.d0,
     & 86.d0,87.d0,88.d0,89.d0,90.d0,91.d0,92.d0,94.d0,96.d0,93.d0,
     & 92.d0,94.d0,95.d0,96.d0,97.d0,98.d0,100.d0,96.d0,98.d0,99.d0,
     & 100.d0,101.d0,102.d0,104.d0,103.d0,102.d0,104.d0,105.d0,106.d0,
     & 108.d0,110.d0,107.d0,109.d0,106.d0,108.d0,110.d0,111.d0,112.d0,
     & 113.d0,114.d0,116.d0,113.d0,115.d0,112.d0,114.d0,115.d0,116.d0,
     & 117.d0,118.d0,119.d0,120.d0,122.d0,124.d0,121.d0,123.d0,120.d0,
     & 122.d0,123.d0,124.d0,125.d0,126.d0,128.d0,130.d0,127.d0,124.d0,
     & 126.d0,128.d0,129.d0,130.d0,131.d0,132.d0,134.d0,136.d0,133.d0,
     & 130.d0,132.d0,134.d0,135.d0,136.d0,137.d0,138.d0,138.d0,139.d0,
     & 136.d0,138.d0,140.d0,142.d0,141.d0,142.d0,143.d0,144.d0,145.d0,
     & 146.d0,148.d0,150.d0,144.d0,147.d0,148.d0,149.d0,150.d0,152.d0,
     & 154.d0,151.d0,153.d0,152.d0,154.d0,155.d0,156.d0,157.d0,158.d0,
     & 160.d0,159.d0,156.d0,158.d0,160.d0,161.d0,162.d0,163.d0,164.d0,
     & 165.d0,162.d0,164.d0,166.d0,167.d0,168.d0,170.d0,169.d0,168.d0,
     & 170.d0,171.d0,172.d0,173.d0,174.d0,176.d0,175.d0,176.d0,174.d0,
     & 176.d0,177.d0,178.d0,179.d0,180.d0,180.d0,181.d0,180.d0,182.d0,
     & 183.d0,184.d0,186.d0,185.d0,187.d0,184.d0,186.d0,187.d0,188.d0,
     & 189.d0,190.d0,192.d0,191.d0,193.d0,190.d0,192.d0,194.d0,195.d0,
     & 196.d0,198.d0,197.d0,196.d0,198.d0,199.d0,200.d0,201.d0,202.d0,
     & 204.d0,203.d0,205.d0,204.d0,206.d0,207.d0,208.d0,209.d0,232.d0,
     & 235.d0,238.d0 /

c     data (aref(l),l = 1,286) /
c    & 0.279d11,0.949d6,0.386d6,0.272d10,4.28d0,52.8d0,0.730d0,4.22d0,
c    & 17.d0,0.999d7,0.111d6,0.312d7,0.115d5,0.237d8,0.904d4,0.476d5,
c    & 843.d0,0.320d7,0.777d4,0.234d6,0.574d5,0.848d6,0.107d6,0.118d6,
c    & 0.849d5,0.922d6,0.467d5,0.310d5,0.104d5,0.489d6,0.386d4,0.217d5,
c    & 103.d0,0.286d4,913.d0,0.850d5,0.160d5,26.d0,0.352d4,5.48d0,
c    & 254.d0,0.592d5,395.d0,82.5d0,0.127d4,2.4d0,114.d0,34.2d0,192.d0,
c    & 175.d0,0.177d4,132.d0,130.d0,0.732d0,292.d0,587.d0,0.113d5,
c    & 0.128d4,319.d0,0.955d4,0.522d5,0.825d6,0.198d5,0.252d4,0.225d4,
c    & 0.337d5,0.129d5,557.d0,0.177d4,449.d0,361.d0,161.d0,613.d0,
c    & 352.d0,51.7d0,236.d0,7.8d0,22.7d0,15.1d0,24.4d0,32.6d0,9.28d0,
c    & 43.4d0,9.28d0,6.56d0,0.55d0,5.6d0,4.7d0,14.7d0,30.9d0,5.7d0,
c    & 5.98d0,5.82d0,0.153d0,0.999d0,5.15d0,5.16d0,25.7d0,7.84d0,
c    & 5.12d0,2.11d0,0.132d0,2.32d0,1.51d0,19.4d0,4.64d0,5.87d0,1.28d0,
c    & 1.96d0,1.98d0,0.32d0,0.698d0,0.378d0,0.236d0,0.406d0,0.425d0,
c    & 0.244d0,0.615d0,0.246d0,0.103d0,0.35d-1,0.236d0,0.234d0,0.316d0,
c    & 0.588d0,0.348d0,0.344d0,0.142d-1,0.155d0,0.31d0,0.38d0,0.368d0,
c    & 0.163d0,0.252d0,0.234d0,0.201d-1,0.143d-1,0.201d0,0.206d0,
c    & 0.388d0,0.197d0,0.463d0,0.121d0,0.79d-2,0.176d0,0.372d-1,
c    & 0.252d-1,0.129d-1,0.555d0,0.293d0,0.925d0,0.328d0,1.25d0,0.177d0,
c    & 0.221d0,0.177d0,0.132d0,0.43d-2,0.124d0,0.428d-1,0.229d0,
c    & 0.342d0,0.909d0,1.53d0,1.63d0,0.9d0,0.571d-2,0.509d-2,0.103d0,
c    & 1.28d0,0.205d0,1.02d0,1.24d0,0.459d0,0.373d0,0.372d0,0.476d-2,
c    & 0.453d-2,0.109d0,0.296d0,0.353d0,0.504d0,3.22d0,0.4d-3,0.446d0,
c    & 0.216d-2,0.283d-2,1.d0,0.126d0,0.167d0,0.225d0,0.1d0,0.197d0,
c    & 0.687d-1,0.142d0,0.477d-1,0.467d-1,0.8d-2,0.399d-1,0.292d-1,
c    & 0.356d-1,0.191d-1,0.689d-1,0.586d-1,0.465d-1,0.508d-1,0.66d-3,
c    & 0.719d-2,0.488d-1,0.676d-1,0.516d-1,0.82d-1,0.721d-1,0.603d-1,
c    & 0.22d-3,0.37d-3,0.922d-2,0.745d-1,0.101d0,0.982d-1,0.111d0,
c    & 0.889d-1,0.35d-3,0.404d-2,0.843d-1,0.576d-1,0.672d-1,0.374d-1,
c    & 0.378d-1,0.32d-3,0.756d-2,0.354d-1,0.543d-1,0.4d-1,0.788d-1,
c    & 0.315d-1,0.357d-1,0.103d-2,0.24d-3,0.793d-2,0.287d-1,0.42d-1,
c    & 0.21d-1,0.541d-1,0.248d-5,0.207d-1,0.17d-3,0.35d-1,0.19d-1,
c    & 0.408d-1,0.38d-1,0.193d-1,0.351d-1,0.12d-3,0.107d-1,0.807d-2,
c    & 0.898d-1,0.109d0,0.178d0,0.277d0,0.247d0,0.414d0,0.17d-3,
c    & 0.105d-1,0.441d0,0.453d0,0.338d0,0.963d-1,0.187d0,0.52d-3,
c    & 0.339d-1,0.574d-1,0.785d-1,0.448d-1,0.101d0,0.233d-1,0.543d-1,
c    & 0.13d0,0.611d-1,0.593d0,0.644d0,1.83d0,0.144d0,0.42d-1,0.573d-2,
c    & 0.181d-1 /



c* Isotopic ratios from Grevesse & Noels 1993 (see solar_abund file in
c* this directory)

c      data (abiso(l),l = 1,286) /
c     &  .9999660, .0000340, .0001419, .9998581, .0749825, .9250175,
c     & 1.0000000, .1988690, .8011310, .9890110, .0109890, .9963276,
c     &  .0036724, .9976158, .0003805, .0020037,1.0000000, .9297542,
c     &  .0022576, .0679883,1.0000000, .7903075, .0997204, .1099720,
c     & 1.0000000, .9222767, .0467140, .0310093,1.0000000, .9501363,
c     &  .0075001, .0421635, .0002001, .7580175, .2419825, .8413676,
c     &  .1583751, .0002574, .9313451, .0014499, .0672050, .9694762,
c     &  .0064686, .0013510, .0207979, .0000393, .0018669,1.0000000,
c     &  .0800333, .0729471, .7378074, .0550229, .0541892, .0025006,
c     &  .9974994, .0435266, .8379060, .0949132, .0236542,1.0000000,
c     &  .0580309, .9171558, .0220117, .0028015,1.0000000, .6825178,
c     &  .2612605, .0112808, .0358474, .0090935, .6915709, .3084291,
c     &  .4863150, .2792543, .0410155, .1872273, .0061880, .6005291,
c     &  .3994709, .2051110, .2740417, .0780094, .3648285, .0780094,
c     & 1.0000000, .0088496, .0901046, .0756235, .2365245, .4971842,
c     &  .0917136, .5067797, .4932203, .0033998, .0221990, .1144394,
c     &  .1146616, .5710857, .1742145, .7081604, .2918396, .0056502,
c     &  .0993066, .0646349, .8304084,1.0000000, .5144610, .1121823,
c     &  .1717791, .1735320, .0280456,1.0000000, .1482353, .0925490,
c     &  .1592157, .1666667, .0956863, .2411765, .0964706, .0553763,
c     &  .0188172, .1268817, .1258065, .1698925, .3161290, .1870968,
c     & 1.0000000, .0102144, .1114947, .2229895, .2733420, .2647101,
c     &  .1172493, .5185185, .4814815, .0124814, .0088798, .1248137,
c     &  .1279185, .2409339, .1223299, .2875062, .0751366, .0429581,
c     &  .9570419, .0097273, .0065894, .0033732, .1451246, .0766153,
c     &  .2418743, .0857673, .3268572, .0462830, .0577884, .5728155,
c     &  .4271845, .0008938, .0257737, .0088961, .0475983, .0710856,
c     &  .1889381, .3180146, .3387999,1.0000000, .0012173, .0010851,
c     &  .0219579, .2728746, .0437026, .2174469, .2643472, .0978511,
c     &  .0795174,1.0000000, .0010598, .0010086, .0242692, .0659053,
c     &  .0785966, .1122172, .7169432, .0008961, .9991039, .0019098,
c     &  .0025022, .8841811, .1114068,1.0000000, .2720348, .1209044,
c     &  .2381816, .0830613, .1716842, .0576714, .0564623, .0308523,
c     &  .1538758, .1126109, .1372927, .0736599, .2657154, .2259931,
c     &  .4779034, .5220966, .0020003, .0217912, .1479012, .2048795,
c     &  .1563873, .2485225, .2185180,1.0000000, .0005577, .0009379,
c     &  .0233708, .1888419, .2560138, .2489164, .2813617,1.0000000,
c     &  .0013950, .0161027, .3360038, .2295827, .2678465, .1490693,
c     & 1.0000000, .0012909, .0304986, .1428110, .2190576, .1613684,
c     &  .3178958, .1270776, .9719575, .0280425, .0015587, .0515035,
c     &  .1863999, .2727804, .1363902, .3513671, .0001198, .9998802,
c     &  .0012785, .2632173, .1428894, .3068361, .2857787, .3547794,
c     &  .6452206, .0001784, .0159063, .0119966, .1334939, .1620360,
c     &  .2646093, .4117796, .3736762, .6263238, .0001270, .0078418,
c     &  .3293576, .3383198, .2524328, .0719210,1.0000000, .0015320,
c     &  .0998763, .1691120, .2312769, .1319899, .2975664, .0686465,
c     &  .2946283, .7053717, .0195326, .1895719, .2058758, .5850197,
c     & 1.0000000,1.0000000, .2404532, .7595468 /


* Isotopic ratios from Asplund et al 2009 (see solar_abund file in
* this directory)

      data (abiso(l),l = 1,286) /
     &  .9999800, .0000200, .0001660, .9998340, .0790000, .9241000,
     & 1.0000000, .1990000, .8010000, .9889380, .0110620, .9977100,
     &  .0022900, .9976210, .0003790, .0020000,1.0000000, .9294310,
     &  .0022280, .0684100,1.0000000, .7899000, .1000000, .1101000,
     & 1.0000000, .9222970, .0468320, .0308720,1.0000000, .9493000,
     &  .0076000, .0429000, .0002000, .7578000, .2422000, .8459460,
     &  .1538080, .0002460, .9313200, .0014700, .0672100, .9694100,
     &  .0064700, .0013500, .0208600, .0000400, .0018700,1.0000000,
     &  .0825000, .0744000, .7372000, .0541000, .0518000, .0025000,
     &  .9975000, .0434500, .8378900, .0950100, .0236500,1.0000000,
     &  .0584500, .9175400, .0211900, .0028200,1.0000000, .6807690,
     &  .2622310, .0113990, .0363450, .0092560, .6917000, .3083000,
     &  .4863000, .2790000, .0410000, .1875000, .0062000, .6010800,
     &  .3989200, .2084000, .2754000, .0773000, .3628000, .0761000,
     & 1.0000000, .0089000, .0937000, .0763000, .2877000, .4961000,
     &  .0873000, .5069000, .4931000, .0036200, .0232600, .1165500,
     &  .1154600, .5690300, .1720800, .7084400, .2915600, .0055800,
     &  .0986780, .0689610, .8267810,1.0000000, .5145000, .1122000,
     &  .1715000, .1738000, .0280000,1.0000000, .1452500, .0915100,
     &  .1583800, .1667200, .0959900, .2439100, .0982400, .0554000,
     &  .0187000, .1276000, .1260000, .1706000, .3155000, .1862000,
     & 1.0000000, .0102000, .1114000, .2233000, .2733000, .2646000,
     &  .1172000, .5183900, .4816100, .0125000, .0089000, .1249000,
     &  .1280000, .2413000, .1222000, .2873000, .0749000, .0429000,
     &  .9571000, .0097000, .0066000, .0034000, .1454000, .0768000,
     &  .2422000, .0859000, .3258000, .0463000, .0579000, .5721000,
     &  .4279000, .0009000, .0255000, .0089000, .0474000, .0707000,
     &  .1884000, .3174000, .3408000,1.0000000, .0012200, .0010800,
     &  .0218800, .2725500, .0437600, .2169300, .2651400, .0979000,
     &  .0795400,1.0000000, .0010600, .0010100, .0241700, .0659200,
     &  .0785400, .1123200, .7169800, .0009100, .9999090, .0018500,
     &  .0025100, .8845000, .1111400,1.0000000, .2704400, .1202300,
     &  .2372900, .0876300, .1713000, .0571600, .0559600, .0307000,
     &  .1499000, .1124000, .1382000, .0738000, .2675000, .2275000,
     &  .4781000, .5219000, .0020000, .0218000, .1480000, .2047000,
     &  .1565000, .2484000, .2186000,1.0000000, .0005600, .0009500,
     &  .0232900, .1888900, .2547500, .2489600, .2826000,1.0000000,
     &  .0013900, .0160100, .3350300, .2286900, .2697800, .1491000,
     & 1.0000000, .0012000, .0298000, .1409000, .2169000, .1610000,
     &  .3203000, .1300000, .9717950, .0282050, .0016200, .0520600,
     &  .1860600, .2729700, .1362900, .3510000, .0001200, .9998800,
     &  .0012000, .2650000, .1431000, .3064000, .2843000, .3466200,
     &  .6433800, .0002000, .0159800, .0127100, .1333700, .1626100,
     &  .2644400, .4107000, .3730000, .6270000, .0001400, .0078200,
     &  .3296700, .3383200, .2524200, .0716300,1.0000000, .0015000,
     &  .0997000, .1687000, .2310000, .1318000, .2986000, .0687000,
     &  .2952400, .7047600, .0199700, .1858200, .2056300, .5885800,
     & 1.0000000,1.0000000, .2428600, .7571200 /

      end
