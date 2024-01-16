**********************************************************************

      program evolcut

***********************************************************************
*     Shorten evolution files during AGB phase
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *

      implicit none

      integer k,l,i,j,idum,earlyAGB,firstpulse,model0,iCend
      integer nbr,nevolfile,entrunit,exitunit
      integer nmod,kmod,cmod
      integer ind,ind2

      parameter (kmod=1000000)

      integer ipmod(kmod),smod(kmod),nphase(kmod),nconv(kmod)

      character*140 string
      character*80 name,namin,namout,nameout
      character*5 charintab(19),charoutab(19)

      logical Cburning

      double precision dum,LHemin,LHemax
      double precision LHe(kmod),LC(kmod),xC(kmod)



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
      
Ccc read LHe and LC
************************************************************************
      namin = trim(name)//'.v2'

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
         read (1,*,err=50,end=50) idum,dum,LHe(l),LC(l)
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

Ccc read XC12
***********************************************************************
      namin = trim(name)//'.c2'

      open (unit=1,file=namin,status='unknown',form='formatted') 

      string = '1234'
      cmod = 0
      do while(string(1:4) .ne. '----')
         read (1,40) string
         cmod = cmod+1
      enddo

      Cburning = .false.
      earlyAGB = 0
      iCend = 0
      do k=1,nmod
         read (1,*,err=65,end=65) idum,xc(k)
      enddo
 65   close (1)


Ccc read nphase
***********************************************************************
      namin = trim(name)//'.hr'

      open (unit=1,file=namin,status='unknown',form='formatted') 

      string = '1234'
      cmod = 0
      do while(string(1:4) .ne. '----')
         read (1,40) string
         cmod = cmod+1
      enddo

      Cburning = .false.
      earlyAGB = 0
      iCend = 0
      do k=1,nmod
         read (1,*,err=70,end=70) idum,nphase(k)
         if (nphase(k).eq.6.and..not.Cburning) Cburning = .true.
         if (earlyAGB.eq.0.and.Cburning.and.xc(k).lt.2.d-2.and.
     &        LC(k).lt.1.d2) then
            if ((nphase(k).eq.6.and.LC(k).lt.1.d2*LHe(k)).or.
     &           nphase(k).eq.5) earlyAGB = k
         endif
c        if (earlyAGB.eq.0.and.Cburning.and.nphase(k).eq.5) earlyAGB = k
      enddo
 70   close (1)


Ccc read nconv 
************************************************************************
      namin = trim(name)//'.v4'

      open (unit=1,file=namin, form='formatted',
     &     status='old', action='read') 

      string = '1234'
      cmod = 0
      do while(string(1:4) .ne. '----')
         read (1,40) string
         cmod = cmod+1
      enddo

      firstpulse = 0
      LHemin = 1.d10
      do k=1,nmod
         read (1,'(122x,i2)',err=80,end=80) nconv(k)
         if (firstpulse.eq.0.and.k.gt.earlyAGB.and.nconv(k).eq.2)
     &        firstpulse = max(1,k-100)
         if (firstpulse.gt.0) then
            LHemin = min(LHe(k),LHemin)
         endif
      enddo
 80   close (1)

      write (*,90) earlyAGB,firstpulse,LHemin
 90   format (' o early AGB phase starts at model ',i7,/
     &     ' o first pulse occurs at model',i8,' ( 0 means no pulse)',/
     &     ' o minimum He luminosity ',1pe9.3' Lsun',/)


Ccc Index of models before the AGB phase and during the pulse (smod) 
Ccc during the inter-pulse phase (ipmod)
************************************************************************
      ind = 1
      ind2 = 1

c       do k = 1,nmod
c          if (LHe(k).lt.1.d3.and.LHe(k).gt.0.d0.and.nphase(k).ge.5) then
c             ipmod(ind) = k
c             ind = ind+1
c          elseif (LHe(k).ge.1.d3.or.nphase(k).lt.5) then
c             smod(ind2) = k
c             ind2 = ind2+1
c          endif
c       enddo
      LHemax = max(1.d3,LHemin*2.d0)
      do k = 1,nmod
         if (nphase(k).ge.5) then
            if (LHe(k).lt.LHemax.or.(nconv(k).eq.1.and.LHe(k).lt.LHemax)
     &           .or.(k.gt.earlyAGB.and.nconv(k).eq.1)) then
               ipmod(ind) = k
               ind = ind+1
            elseif (LHe(k).ge.LHemax) then
               smod(ind2) = k
               ind2 = ind2+1
            endif
         else
            smod(ind2) = k
            ind2 = ind2+1
         endif
      enddo

c..   select models to be saved
      l = 1
      do k = 1,ind,nbr
         ipmod(l) = ipmod(k)
         l = l+1
      enddo

Ccc Read/Write files
C   ************************************************

      nevolfile = 19
      nameout = trim(name)//'_cut'

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
            read (entrunit,110,err=100,end=100) string
            if (i.ge.ipmod(l)) then
               write (exitunit,110) string
               l = l+1
            elseif (i.ge.smod(j)) then
               write (exitunit,110) string
               j = j+1
            endif
         enddo

 100     close(entrunit)
         close(exitunit)
         if (k.eq.1) write (*,120) l+j

      enddo

      write (*,150) l+j

 110  format (A140)
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
