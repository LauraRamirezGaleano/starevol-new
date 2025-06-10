************************************************************************

      PROGRAM evoldeline

************************************************************************
c     Output files truncation between user defined sequences
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *

      implicit none

      integer idiffnuc,modi,modf,model,entrunit,exitunit,
     &     lenci,i,ncount,ftype,ipass

      integer shell1,iold,k,k0,k1,jdx,imod,isave
      double precision version

      character*120 name,namin,namout
      character*150 string
      character*7 cmodel
      character*5 extsn(21),extsnx(21)
      character*5 extevol(23),extevolx(23)
      character*5 extmass(20),extmassx(20)
      character*5 extalpha(12),extalphax(12)
      character*5 charintab(21),charoutab(21)

      logical list

      data (extsn(i),i=1,21) /'.v1','.v2','.v3','.v4',
     &     '.v5','.v6','.v7','.v8','.v9','.v10','.v11','.v12','.c1',
     &     '.c2','.c3','.c4','.s1','.s2','.s3','.s4','.hr'/
      data (extsnx(i),i=1,21) /'.v1s','.v2s','.v3s','.v4s',
     &     '.v5s','.v6s','.v7s','.v8s','.v9s','.v10s','.v11s','.v12s',
     &     '.c1s','.c2s','.c3s','.c4s','.s1s','.s2s','.s3s','.s4s',
     &     '.hrs'/
      data (extevol(i),i=1,23) /'.v1','.v2','.v3','.v4',
     &     '.v5','.v6','.v7','.v8','.v11','.v12','.v13','.c1','.c2',
     &     '.c3','.c4','.s1','.s2','.s3','.s4','.hr','.as','.tc1',
     &     '.tc2'/
      data (extevolx(i),i=1,23) /'.v1s','.v2s','.v3s','.v4s',
     &     '.v5s','.v6s','.v7s','.v8s','.v11s','.v12s','.v13s','.c1s',
     &     '.c2s','.c3s','.c4s','.s1s','.s2s','.s3s',
     &     '.s4s','.hrs','.ass','.tc1s','.tc2s'/
      data (extmass(i),i=1,20) /'.v1','.v2','.v3','.v4',
     &     '.v5','.v6','.v7','.v8','.v11','.v12','.c1','.c2',
     &     '.c3','.c4','.c5','.s1','.s2','.s3','.s4','.hr'/
      data (extmassx(i),i=1,20) /'.v1s','.v2s','.v3s','.v4s',
     &     '.v5s','.v6s','.v7s','.v8s','.v11s','.v12s','.c1s',
     &     '.c2s','.c3s','.c4s','.c5s','.s1s','.s2s','.s3s',
     &     '.s4s','.hrs'/
c       data (extalpha(i),i=1,13) /'.v1','.v2','.v3','.v4',
c      &     '.v5','.v6','.v7','.v8','.c1','.c2','.s1',
c      &     '.s2','.hr'/
c       data (extalphax(i),i=1,13) /'.v1s','.v2s','.v3s','.v4s',
c      &     '.v5s','.v6s','.v7s','.v8s','.c1s','.c2s','.s1s',
c      &     '.s2s','.hrs'/
      data (extalpha(i),i=1,12) /'.v1','.v2','.v3','.v4',
     &     '.v5','.v6','.v7','.c1','.c2','.s1','.s2','.hr'/
      data (extalphax(i),i=1,12) /'.v1s','.v2s','.v3s','.v4s',
     &     '.v5s','.v6s','.v7s','.c1s','.c2s','.s1s','.s2s','.hrs'/


      ftype = 0
      version = 0.d0
c     modf=10000
c      write (*,*) 'Enter star directory/filename (within quote)' //
c     &     ' starting and ending model :'
c      write (*,*) 'Beware input files will be changed!!!'
      read (5,*) modi,modf,name
c     name='test/m15z02'
c     modi=10
c     modf=20

      namin= trim(name)//'.hr'
      print *, namin
      entrunit = 10
      open (unit = entrunit, file = trim(namin), 
     &     form= 'formatted',status='old',action='read')

      read (entrunit,'(A150)') string
      write (*,*) string
      if (string(1:4).eq.' sn '.or.string(1:4).eq.' SN ') then
         ftype = 1
         ncount = 21
         backspace(entrunit)
         read (entrunit,100) version
 100     format (22x,0pf4.2)
      endif
      if (string(1:4).eq.' evo'.or.string(1:4).eq.' AGB') then
         ftype = 2
         backspace(entrunit)
         read (entrunit,101) version
         ncount = 23
 101     format (23x,0pf4.2)
      endif
      if (string(1:4).eq.' alp'.or.string(1:4).eq.' ALP') then
         ftype = 3
         ncount = 12
         backspace(entrunit)
         read (entrunit,102) version
 102     format (24x,0pf4.2)
      endif
      if (string(1:4).eq.' MAS') then
         ftype = 4
         ncount = 20
         backspace(entrunit)
         read (entrunit,103) version
 103     format (27x,0pf4.2)
      endif

      if (ftype.eq.0) then
         write (*,*) 'Banner not recognized'
         stop ' evoldeline'
      endif

      write (*,'("  Recognize version : ",0pf4.2,
     &     ", number output files :",i2)') version,ncount
      if (ftype.eq.1) then
         do k = 1,ncount
            charintab(k) = extsn(k)
            charoutab(k) = extsnx(k)
         enddo
      endif
      if (ftype.eq.2) then
         do k = 1,ncount
            charintab(k) = extevol(k)
            charoutab(k) = extevolx(k)
         enddo
      endif
      if (ftype.eq.3) then
         do k = 1,ncount
            charintab(k) = extalpha(k)
            charoutab(k) = extalphax(k)
         enddo
      endif
      if (ftype.eq.4) then
         do k = 1,ncount
            charintab(k) = extmass(k)
            charoutab(k) = extmassx(k)
         enddo
      endif
      close(entrunit)

      do k = 1,ncount,1

         entrunit=entrunit+1
         exitunit=entrunit+50

c         namin= name(1:lenci(name))//charintab(k)
c         namout= name(1:lenci(name))//charoutab(k)

         namin= trim(name)//charintab(k)
         namout= trim(name)//charoutab(k)

         write (*,'("processing file : ",A)') namin
         open (unit = entrunit,file = trim(namin), 
     &        form= 'formatted',status='old',action='read')

         open (unit = exitunit,file = trim(namout), 
     &     form= 'formatted',status='unknown')

         list=.true.

         do while(string(1:4) .ne. '----')
            read (entrunit,'(A150)',err=150,end=200) string
            write (exitunit,'(A150)') string
         enddo
         read (entrunit,'(A150)',err=150) string
         write (exitunit,'(A150)') string
         read (entrunit,'(A150)',err=150) string
         write (exitunit,'(A150)') string

         ipass = 0
 50      read (entrunit,'(A7)',err=120,end=200) cmodel
         read (cmodel,'(i7)') model
         backspace(entrunit)
         read (entrunit,'(A150)',err=150,end=200) string
         if (model.eq.0) model = modi
         list = .true.
         if (model.ge.modi.and.model.le.modf.and.ipass.eq.0) 
     &        ipass = 1
         if (model.ge.modi.and.model.le.modf.and.ipass.eq.1) 
     &        list=.false.
         if (model.ge.modf.and.ipass.eq.1) ipass = 2
         if (list) write (exitunit,'(A150)') string
         goto 50

 120     write (*,'("read error - model number :",A7)') cmodel
         stop
 150     write (*,'("read error - line :",A150)') string
         stop

 200     close(entrunit)
         close(exitunit)
         if (ipass.eq.0) write (*,*) '   no lines removed'
      enddo

      if (ftype.eq.1) call delc6(modi,modf,name)

      end


************************************************************************

      SUBROUTINE delc6(modi,modf,name)

************************************************************************
c     truncate the file name.c6 after the shell nb : shell1

      implicit none

      integer imod,nmod,entrunit,entrunitsav,nsp
      integer i,k,dummy,ipass,lenci
      integer model,modi,modf

      double precision xsp,time,system

      character*120 name,namein,namesave,cmd

      logical list

      parameter (nsp = 1000, nmod = 500000)

      dimension xsp(nsp)

      !external system


*___________________
***   mode evolution
*-------------------

      entrunit = 26
      entrunitsav = 27
      namein = trim(name) // '.c6'
      namesave= trim(name) // '.sav6'
      cmd = 'cp -f ' // trim(namein) // ' ' // trim(namesave)
      write (*,*) cmd 
      !dummy=system(cmd)
      open (unit = entrunit,file = trim(namein),form = 'unformatted')
      open (unit = entrunitsav,file = trim(namesave),
     &     form = 'unformatted')

      ipass = 0
      imod = 0
      list=.true.
      do 55 k=1,nmod
         read (entrunitsav,err=55,end=10) model,time,(xsp(i),i = 1,nsp)
         list = .true.
         if (model.ge.modi.and.model.le.modf.and.ipass.eq.0) 
     &        ipass = 1
         if (model.ge.modi.and.model.le.modf.and.ipass.eq.1) 
     &        list=.false.
         if (model.gt.modf.and.ipass.eq.1) ipass = 2
         if (list) then
            write (entrunit) model,time,(xsp(i),i = 1,nsp)
            imod = imod+1
         endif
 55   continue

c      iold = 0
c      k0 = 0
c      jdx=-1
c      do k=1,nmod
c         k1 = k-k0
c         isave = imod
c         read (26,end=10,err=10) imod,time(k1),(xsp(k1,i),i = 1,nsp)
c         write (*,*) imod,k,k1
c         if (imod.gt.shell1) goto 10
c         write (27) imod,time(k1),(xsp(k1,i),i = 1,nsp)
c         if (k.gt.1.and.imod.le.iold) then
c            k0 = iold-imod+1
c            time(k1-k0) = time(k1)
c            xmod(k1-k0) = dble(imod)
c            iold = imod
c            do i=1,nsp
c               xsp(k1-k0,i) = xsp(k1,i)
c            enddo
c         else
c            iold = imod
c            xmod(k1) = dble(imod)
c         endif
c      enddo

 10   close(entrunit)
      close(entrunitsav)
      
      write (*,*) 'last recorded/read model :',imod
      write (*,*) 'original file : ', name(1:lenci(name)) // '.sav6'
      write (*,*) 'new file      : ', name(1:lenci(name)) // '.c6'
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












