

************************************************************************

      SUBROUTINE convzone

************************************************************************
*                     Delineate the convective boundaries
* novlim(kl,1-2)   -->  Schwarzschild criterion
* novlim(kl,3-4)   -->  CURRENT LIMITS and parametric overshoot   ! note td : 3 - base ZC et 4 - haut ZC (nsconv : nbre ZC)
* novlim(kl,5-6)   -->  Ledoux criterion
* novlim(kl,7-8)   -->  diffusive overshoot and mesh
* novlim(kl,9-10)  -->  semiconvective zones
* novlim(kl,11-12) -->  thermohaline zones
*
* crzc = c  -2  : convective shell      abled < abrad or abadd < abrad
* crzc = a  -1  : atmospheric shell     tau < 10
* crzc = s  -3  : semi-convective shell abadd < abrad < abled  (Ledoux)
*                 shell mixed : semiconv+diffzc
*        S   3  : if Ledoux criterion used, shell *not* mixed (semiconv)
* crzc = t   2  : thermohaline shell    abled < abrad < abadd
* crzc = o   1  : overshoot layers
* crzc = r   4  : radiative shell       abrad < abled or abrad < abadd
*
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
C... modif NL 
      include 'evolcom.diff'
      include 'evolcom.therm2'
      include 'evolcom.transp'
      include 'evolcom.nuc'
      include 'evolcom.spec'
C...fin modif 
      include 'evolcom.grad'
      include 'evolcom.mod'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer kl,nconv1,nconv2,nrad,nshdel
      integer delnconv,nc,ncz
      integer nsconv0,novlim0
      integer i,j,k,jl,it
      integer ishift,kc,kc0,kc1,imin
      integer klenv,klpulse,klcore
C... modif NL 
      integer nmodmdix,id,ip,nshell,ibth,itth
C.. fin modif 
      

      double precision extabrad
      double precision delabrad
      double precision delmconv
C... modif NL
      double precision abmuliss,abmax,abmin,gradledoux
C...fin modif 

      logical isconv,isconvS,isconvT,isconvTnuc,isconvL,
     &       tbot,ttop,neutrality

      dimension isconv(nsh),isconvL(nsh),isconvS(nsh),isconvT(nsh),
     &     novlim0(nmaxconv,ntypeconv),isconvTnuc(nsh)
C... modif NL
      dimension abmuliss(nsh),gradledoux(nsh)
C...fin modif      
      
      common /overshoot/ klcore,klenv,klpulse

***   define gradient differences
      do i = 1,nmod1
         gradledoux(i) = abled(i)-abrad(i)
         isconv(i) = abrad(i).gt.abm(i)   !Schwarzchild
         isconvL(i) = abrad(i).gt.abled(i).and.isconv(i)  !Ledoux
         isconvS(i) = isconv(i).and.abmu(i).gt.0.d0  !semiconvection
	 isconvT(i) = .not.isconv(i).and.abmu(i).lt.0.d0.
     &                  and.gradledoux(i).gt.0.d0  !thermohaline
	 isconvTnuc(i) = .not.isconv(i).and.abmu(i).lt.0.d0.
     &                   and.gradledoux(i).gt.0.d0.
     &                   and.xsp(i,ihe3).ge.1.d-08 !thermohaline
     
     
      enddo
      isconv(nmod) = .false.
      isconvS(nmod) = .false.
      isconvT(nmod) = .false.
      isconvTnuc(nmod) = .false.
      isconvL(nmod) = .false. 

      scr = 'r'
      do i = 1,nmod
         crz(i) = 4
      enddo

***   define novlim
      nsconv0 = nsconv
      do kl = 1,ntypeconv
         do j = 1,nmaxconv
            novlim0(j,kl) = novlim(j,kl)
            novlim(j,kl) = 0
         enddo
      enddo
      nsconv = 0
      nschwar = 0
      nledoux = 0
      nsemiconv = 0
      nthermha = 0

***   define Schwarzschild and Ledoux regions
      do i = 1,nmod
         if (isconv(i)) crz(i) = -2
      enddo
      it = -3
      if (ledouxmlt) it = 3
      do j = 1,2
         kl = 1
         jl = 1
         if (j.eq.2) then
            do i = 1,nmod
               if (isconvS(i)) crz(i) = it
            enddo
            jl = 5
         endif
         if (isconv(1)) novlim(kl,jl) = 1
         do k = 2,nmod1
            tbot = crz(k).ne.-2.and.crz(k+1).eq.-2
            ttop = crz(k).eq.-2.and.crz(k+1).ne.-2
            if (tbot) novlim(kl,jl) = k+1
            if (ttop.and.k-novlim(kl,jl).gt.nczm) then
               novlim(kl,jl+1) = k
               kl = kl+1
               if (kl.gt.nmaxconv) then
                  write (nout,*) 'too many convective zones !',j
                  goto 10
               endif
            endif
         enddo
 10      if (j.eq.1) nschwar = kl-1
         if (j.eq.2) nledoux = kl-1
      enddo

***   define semiconvective and thermohaline zones
      do i = 1,nmod
         if (isconvT(i)) crz(i) = 2
      enddo
      do j = 1,2
         kl = 1
         jl = 2*j+7
         if (j.eq.2) it = 2
         if (crz(1).eq.it) novlim(kl,jl) = 1
         do k = 2,nmod1
            tbot = crz(k).ne.it.and.crz(k+1).eq.it
            ttop = crz(k).eq.it.and.crz(k+1).ne.it
            if (tbot) novlim(kl,jl) = k+1
            if (ttop.and.k-novlim(kl,jl).ge.nczm) then
               novlim(kl,jl+1) = k
               kl = kl+1
               if (kl.gt.nmaxconv) then
                  write (nout,*) kl,'too many unstable zones !',j
                  goto 20
               endif
            endif
         enddo
 20      if (j.eq.1) nsemiconv = kl-1
         if (j.eq.2) nthermha = kl-1
      enddo

      if (ledouxmlt) then
         nsconv = nledoux
         do kl = 1,nsconv
            novlim(kl,3) = novlim(kl,5)
            novlim(kl,4) = novlim(kl,6)
         enddo
      else
         nsconv = nschwar
         do kl = 1,nsconv
            novlim(kl,3) = novlim(kl,1)
            novlim(kl,4) = novlim(kl,2)
         enddo
      endif

C... modif NL 
      do i= 2,nmod 
         if (.not.isconvTnuc(i-1).and.isconvTnuc(i)) then 
	     novlim(1,11)= i
         goto 22 
	 endif 
      enddo
 22   do i= 1,nmod-1 
	 if (.not.isconvTnuc(i+1).and.isconvTnuc(i)) then 
	     novlim(1,12)= i
         goto 25 
	 endif
      enddo       
 

*___________________________________________________________________
***   suppress convective zones whose width is less than nczm shells
***   additional suppression during agbphase if shell mass too small
*-------------------------------------------------------------------

c      if (nphase.ge.5.and.nsconv.ge.2.and.nczm.gt.0) then
 25   ncz = abs(nczm)
      if (nsconv.ge.2.and.nczm.ne.0) then
         k = 0
 50      k = k+1
         delmconv = 1.d99
         if (agbphase) delmconv = (m(novlim(k,4))-m(novlim(k,3)))/msun
         if (nczm.gt.0.or.t(novlim(k,4)).lt.1.d6) then
            delnconv = novlim(k,4)-novlim(k,3)
         else
            delnconv = ncz+1
         endif
         if (delmconv.lt.delmcz.or.delnconv.le.ncz) then
c..   the suppressed convective shells are assumed to be radiative
            do i = novlim(k,3),novlim(k,4)
c               crz(i) = 4
               if (crz(i).eq.-2) crz(i) = 4
            enddo
            do kl = k,nsconv-1
               novlim(kl,3) = novlim(kl+1,3)
               novlim(kl,4) = novlim(kl+1,4)
            enddo
            novlim(nsconv,4) = 0
            novlim(nsconv,3) = 0
            k = k-1
            nsconv = nsconv-1
         endif
         if (k.lt.nsconv) goto 50
      endif

*____________________________________________________
***   merging of adjacent convective zones,
***   if separated by less than nczm radiative shells
*----------------------------------------------------

      if (nphase.gt.1.and.nsconv.gt.1.and.nczm.ne.0) then
         nc = nsconv
         do kl = 1,nsconv-1
            nconv1 = novlim(kl,4)-novlim(kl,3)+1
            nconv2 = novlim(kl+1,4)-novlim(kl+1,3)+1
            if (nczm.gt.0.or.t(novlim(kl,4)).lt.1.d6) then
               nrad = novlim(kl+1,3)-novlim(kl,4)-1
            else
               nrad = ncz+1
            endif
            if (nrad.lt.nconv1.and.nrad.le.ncz.and.nrad.le.nconv2)
     &           then
               delabrad = abrad(novlim(kl+1,3))-abrad(novlim(kl,4))
               nshdel = novlim(kl+1,3)-novlim(kl,4)
               do k = novlim(kl,4)+1,novlim(kl+1,3)-1
                  abrad(k) = abrad(novlim(kl,4))+delabrad/dble(nshdel)*
     &                 dble(k-novlim(kl,4))
c                  abrad(k) = max(abrad(k),abm(k)*1.001d0)
                  crz(k) = 4
               enddo
               novlim(kl,4) = novlim(kl+1,4)
               do j = kl+1,nc-1
                  novlim(j,3) = novlim(j+1,3)
                  novlim(j,4) = novlim(j+1,4)
               enddo
               novlim(nc,4) = 0
               novlim(nc,3) = 0
               nc = nc-1
               if (nc.eq.nsconv-1) goto 60
            endif
         enddo
 60      nsconv = nc
      endif

      do kl = 1,nsconv
         novlim(kl,7) = novlim(kl,3)
         novlim(kl,8) = novlim(kl,4)
      enddo
      if (nsconv.ge.1) then
         if (ledouxmlt) then
            scr = 'l'
         else
            scr = 'c'
         endif
      endif

***   define atmospheric layers
      if (nsconv.gt.0) then
         if (t(novlim(nsconv,4)).lt.1.d6) then
            imin = novlim(nsconv,4)
             do i = imin,nmod
               crz(i) = -1
               Dconv(i) = Dconv(imin-1)
            enddo
         endif
      endif

*____________________________________________________________________
***   search for gradient neutrality (abrad = abad) at all convective 
***   boundaries during the first itermix iterations
*--------------------------------------------------------------------

      neutrality = idup.eq.1.or.idup.eq.3  
c      if (nphase.eq.2) neutrality = .true. ! modif Ana Palacios 09/2018
      if (nsconv.gt.0.and.neutrality.and.iter.le.abs(itermix)) then
         do kl = 3,4
            if (kl.eq.3) then
               ishift = 1
            else
               ishift = -1
            endif
            do k = 1,nsconv
               kc = novlim(k,kl)
c..  do not treat upper part of conv. envelope
               if (k.eq.klenv.and.kl.eq.4) kc = 0
               if (kc.gt.1) then
                  kc1 = kc+ishift
                  kc0 = kc-ishift
                  if (ishift.eq.1) then
                     delabrad = dm(kc0)/dm(kc)
                  else
                     delabrad = dm(kc)/dm(kc1)
                  endif
                  extabrad = abrad(kc)-delabrad*(abrad(kc1)-abrad(kc))
c                  if (extabrad.gt.abm(kc0).and.extabrad.lt.abrad(kc))
c     &                 then
	          if (extabrad.gt.abm(kc0)) then
                     novlim(k,kl) = kc0
                     novlim(k,kl-2) = novlim(k,kl)
                     !abrad(kc0) = extabrad
                     crz(kc0) = -2
                  endif
               endif
            enddo
         enddo
      endif

c..   check if convective boundaries have changed
      if (chkmix) then
         if (nsconv.ne.nsconv0) then
            lmix = .true.
         else
            do j = 1,nsconv
               do kl = 3,4
                  if (novlim0(j,kl).ne.novlim(j,kl)) then
                     lmix = .true.
                     return
                  endif
               enddo
            enddo
         endif
      endif

 2000 format (i8,2x,10(1x,0pf10.7),2x,1x,i2,a1,a2)

      return
      end
