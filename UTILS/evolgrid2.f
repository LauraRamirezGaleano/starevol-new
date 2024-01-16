**********************************************************************

      program evolgrid2

***********************************************************************
*     Shorten evolution files (all phases)
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *

      implicit none

      integer k,l,i,j,idum,model0
      integer nbr,nevolfile,entrunit,exitunit
      integer nmod,kmod,cmod
      integer ind,ind2
      integer Discret

      parameter (kmod=1000000)
      parameter (Discret=300)

      integer modsave(kmod)

      character*140 string
      character*80 name,namin,namout,nameout
      character*5 charintab(19),charoutab(19)

      logical Cburning

      double precision Lmin,Lmax,Teffmin,TEffmax,rhocmin,rhocmax,Tcmin,
     &     Tcmax,H1min,H1max,He4min,He4max,tmin,tmax
      double precision dL,dTeff,drhoc,dTc,dH1,dHe4,dt
      double precision Lold,Teffold,rhocold,Tcold,H1old,He4old,told
      double precision dum,LHemin,LHemax
      double precision Ls(kmod),Teff(kmod),t(kmod),rhoc(kmod),tc(kmod),
     &     H1(kmod),He4(kmod)



      data (charintab(i),i=1,19) /'.v1','.v2','.v3','.v4',
     &     '.v5','.v6','.v7','.v8','.v11','.v12','.c1','.c2',
     &     '.c3','.c4','.s1','.s2','.s3','.s4','.hr'/
      data (charoutab(i),i=1,19) /'.v1','.v2','.v3','.v4',
     &     '.v5','.v6','.v7','.v8','.v11','.v12','.c1',
     &     '.c2','.c3','.c4','.s1','.s2','.s3',
     &     '.s4','.hr'/

Ccc Frequence number of models to save during the interpulse
************************************************************************
      nbr = 70

Ccc Star filename
      write (*,*) 'Enter star filename :'
      read (5,*) name
************************************************************************
      
Ccc read L, Teff and t
************************************************************************
      namin = trim(name)//'.hr'

      open (unit=1,file=namin, form='formatted',
     &     status='old', action='read') 

      string = '1234'
      cmod = 0
      do while(string(1:4) .ne. '----')
         read (1,40) string
         cmod = cmod+1
      enddo
 40   format(A4)

      l=1
      do k=cmod,kmod
         read (1,*,err=50,end=50) idum,dum,Ls(l),dum,dum,Teff(l),dum,dum
     &        ,dum,dum,dum,t(l)
         if (l.eq.1) model0 = idum
         l=l+1
      enddo
 50   close (1)
      nmod = l-1
      write (*,60) nmod,name
 60   format (/,i6,' lines read in file ',a,/)
      if (nmod.gt.kmod) then
         write (*,'("Too many lines : nmod > kmod")')
         stop
      endif

Ccc read rhoc and Tc
***********************************************************************
      namin = trim(name)//'.v1'

      open (unit=1,file=namin,status='unknown',form='formatted') 

      string = '1234'
      cmod = 0
      do while(string(1:4) .ne. '----')
         read (1,40) string
         cmod = cmod+1
      enddo

      do k=1,nmod
         read (1,*,err=65,end=65) idum,Tc(k),dum,dum,rhoc(k)
      enddo
 65   close (1)


Ccc read Xc and Yc
***********************************************************************
      namin = trim(name)//'.c1'

      open (unit=1,file=namin,status='unknown',form='formatted') 

      string = '1234'
      cmod = 0
      do while(string(1:4) .ne. '----')
         read (1,40) string
         cmod = cmod+1
      enddo

      do k=1,nmod
         read (1,*,err=70,end=70) idum,dum,H1(k),dum,dum,He4(k)
      enddo
 70   close (1)

Ccc Select models to be saved
C   ************************************************   
      l = 1

      Lmin = 1.d99
      Lmax = -1.d0
      Teffmin = 1.d99
      Teffmax = -1.d0
      rhocmin = 1.d99
      rhocmax = -1.d0
      Tcmin = 1.d99
      Tcmax = -1.d0
      tmin = 1.d99
      tmax = -1.d0
      H1min = 1.d99
      H1max = -1.d0
      He4min = 1.d99
      He4max = -1.d0

      do k = 1,nmod
         if (Ls(k).le.Lmin) Lmin = Ls(k)
         if (Ls(k).ge.Lmax) Lmax = Ls(k)
         if (Teff(k).le.Teffmin) Teffmin = Teff(k)
         if (Teff(k).ge.Teffmax) Teffmax = Teff(k)
         if (rhoc(k).le.rhocmin) rhocmin = rhoc(k)
         if (rhoc(k).ge.rhocmax) rhocmax = rhoc(k)
         if (Tc(k).le.Tcmin) Tcmin = Tc(k)
         if (Tc(k).ge.Tcmax) Tcmax = Tc(k)
         if (H1(k).le.H1min) H1min = H1(k)
         if (H1(k).ge.H1max) H1max = H1(k)
         if (He4(k).le.He4min) He4min = He4(k)
         if (He4(k).ge.He4max) He4max = He4(k)
         if (t(k).le.tmin) tmin = t(k)
         if (t(k).ge.tmax) tmax = t(k)
      enddo

      dL = dabs(dlog10(Lmax)-dlog10(Lmin))/dble(Discret)
      dTeff = dabs(dlog10(Teffmax)-dlog10(Teffmin))/dble(Discret)
      drhoc = dabs(dlog10(rhocmax)-dlog10(rhocmin))/dble(Discret)
      dTc = dabs(dlog10(Tcmax)-dlog10(Tcmin))/dble(Discret)
      dH1 = dabs(H1max-H1min)/dble(Discret)
      dHe4 = dabs(He4max-He4min)/dble(Discret)
      dt = dabs(tmax-tmin)/dble(Discret)
      
      Lold = Ls(1)
      Teffold = Teff(1)
      rhocold = rhoc(1)
      Tcold = Tc(1)
      H1old = H1(1)
      He4old = He4(1)
      told = t(1)

      ind = 1
      modsave(ind) = 1
      do k = 2,nmod
         if (dabs(dlog10(Ls(k))-dlog10(Lold))>=dL .or.
     &        dabs(dlog10(Teff(k))-dlog10(Teffold))>=dTeff .or.
     &        dabs(dlog10(rhoc(k))-dlog10(rhocold))>=drhoc .or.
     &        dabs(dlog10(Tc(k))-dlog10(Tcold))>=dTc .or.
     &        dabs(H1(k)-H1old)>=dH1 .or.
     &        dabs(He4(k)-He4old)>=dHe4 .or.
     &        dabs(t(k)-told)>=dt) then
            ind = ind+1
            modsave(ind) = k

            Lold = Ls(k)
            Teffold = Teff(k)
            rhocold = rhoc(k)
            Tcold = Tc(k)
            H1old = H1(k)
            He4old = He4(k)
            told = t(k)
         endif
      enddo
 
Ccc Read/Write files
C   ************************************************

      nevolfile = 19
      nameout = trim(name)//'_gridhrd'

      do k = 1,nevolfile,1

         namin = trim(name)//charintab(k)
         namout = trim(nameout)//charoutab(k)

         entrunit = entrunit+15
         exitunit = entrunit+50

         write (*,'("processing file : ",A)') namin
         open (unit = entrunit,file = trim(namin),
     &        status='old', form= 'formatted', action='read')

         open  (unit = exitunit,file = trim(namout), 
     &        status = 'unknown', form= 'formatted')

         string = '1234'
         cmod = 0
         do while(string(1:4) .ne. '----')
            read (entrunit,110,err=100,end=100) string
            write (exitunit,110) string
            cmod = cmod+1
         enddo

         read (entrunit,110,err=100) string
         write (exitunit,110) string
         read (entrunit,110,err=100) string
         write (exitunit,110) string

         l = 1
         j = 1
         do i = 1,nmod
            read (entrunit,111,err=100,end=100) idum,string
            if (i.ge.modsave(l)) then
               write (exitunit,111) l,string
               l = l+1
            endif
         enddo

 100     close(entrunit)
         close(exitunit)
         if (k.eq.1) write (*,120) l+j

      enddo

      write (*,150) l+j

 110  format (A140)
 111  format (i8,1x,a139)
 120  format (i6,' lines recorded')
 150  format (/,5x,'number of lines',i7)

      end

************************************************************************

      INTEGER FUNCTION LENCI (chain)

************************************************************************
*   Last non-blank character position in string a                      *
************************************************************************

      character*(*) chain

      integer i,n

      n = len(chain)
      do i = 1,n
         if (chain(i:i).eq.' ') then
            lenci = i-1
            return
         end if
      enddo

      lenci = 0

      return
      end
