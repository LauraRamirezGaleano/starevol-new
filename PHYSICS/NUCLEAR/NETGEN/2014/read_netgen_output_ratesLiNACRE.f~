      program main
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     version history
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c 03/07/2001
c     Bug fixed in statements 
c
c               if(aflag(j)(7:8).eq.'--') then
c                   flag(ireac) = -10.
c               else if(aflag(j)(6:8).eq.'---') then
c                   flag(ireac) = -100.
c               else if(aflag(j)(5:8).eq.'----') then
c                   flag(ireac) = -200.
c               else if(aflag(j)(8:8).eq.'+') then
c                   flag(ireac) = -11.
c               else if(aflag(j)(7:8).eq.'++') then
c                   flag(ireac) = -13.
c
c               Their order must be reversed!!
c
c 06/02/2003    Routine rewritten to handle density-dependent rates
c 22/04/2014    Routine adapted to comply with version 3.30 of STAREVOL (AP)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



      implicit double precision (a-h,o-z)

c.. Actual number of reactions = 2 (see file Li7palphaBe7pgammaNACREIImax.network)
c      parameter(ngrid=250,nre=10000) 
      parameter(ngrid=250,nre=2,nreac=167) 

c.. Lines added to compute reaction rates for the desired temperature grid
      parameter (nvit =  90)
      double precision tnuc9
      integer echem
      integer sreacnum(nre)

      common/cvit/v(NRE),tgrid(ngrid),aNegrid(nre,11)
      common/cvgrid/vgrid(ngrid,nre,11),flag(nre),ndata(nre)
      common/network/n1(nre),n2(nre),n3(nre),n4(nre),
     &               z1(nre),z2(nre),z3(nre),z4(nre),
     &               a1(nre),a2(nre),a3(nre),a4(nre),
     &               reaction(nre),ireac,kgrid
      common/Q/Qrad(nre),Qnu(nre)

      double precision anuc,znuc
      dimension tnuc9(nvit)
      dimension species(100),icode(100),anuc(100),znuc(100),echem(30,60)
      dimension ireacode(nre,4)
      dimension ireac(4)
      character z1*6,z2*6,z3*6,z4*6
      character species*5
      character reaction*37
      character reference*172
      dimension reference(nre)
      dimension f(0:4),vv(nvit,nre)
      double precision vr(167)

c.. Definition of the desired temperature grid using a unit of 10^9 K for T
      data (tnuc9(i),i = 1,nvit) /
     & 5.000000D-04,8.800000D-04,1.270000D-03,1.670000D-03,2.080000D-03,
     & 2.500000D-03,2.930000D-03,3.370000D-03,3.820000D-03,4.280000D-03,
     & 4.750000D-03,5.230000D-03,5.720000D-03,6.220000D-03,6.730000D-03,
     & 7.250000D-03,7.780000D-03,8.320000D-03,8.870000D-03,9.430000D-03,
     & 1.000000D-02,1.058000D-02,1.152000D-02,1.282000D-02,1.448000D-02,
     & 1.650000D-02,1.888000D-02,2.162000D-02,2.472000D-02,2.818000D-02,
     & 3.200000D-02,3.618000D-02,4.072000D-02,4.562000D-02,5.088000D-02,
     & 5.650000D-02,6.248000D-02,6.882000D-02,7.552000D-02,8.258000D-02,
     & 9.000000D-02,9.778000D-02,1.068000D-01,1.170600D-01,1.285590D-01,
     & 1.412980D-01,1.552780D-01,1.704970D-01,1.869550D-01,2.046540D-01,
     & 2.235920D-01,2.437710D-01,2.651890D-01,2.878470D-01,3.117450D-01,
     & 3.368820D-01,3.632600D-01,3.908770D-01,4.197340D-01,4.498310D-01,
     & 4.811680D-01,5.137450D-01,5.475610D-01,5.826180D-01,6.189140D-01,
     & 6.564500D-01,6.952250D-01,7.352410D-01,7.764970D-01,8.189920D-01,
     & 8.627270D-01,9.077020D-01,9.539170D-01,1.001370D+00,1.050070D+00,
     & 1.100000D+00,1.151170D+00,1.232180D+00,1.343020D+00,1.483690D+00,
     & 1.654190D+00,1.854530D+00,2.084700D+00,2.344700D+00,2.634530D+00,
     & 2.954190D+00,3.303690D+00,3.683020D+00,4.092180D+00,4.531170D+00/

c..   Index of the associated reaction within the STAREVOL v3.30 network
      data (sreacnum(i),i=1,nre) /14, 17/

c===========================================
c     factorials

      f(0) = 1.
      f(1) = 1.
      f(2) = 2.
      f(3) = 6.
      f(4) = 24.

c...  First call to vit to read the netgen output files using 
c...  an ad-hoc value for T and the density (iread = 0)

      T = 9.0D7
      rho = 100.D0
      iread = 0
      pme = 0.3D0

      call vit(T,rho,pme,iread)

c..   Retrieve the information on the reactions to plug it back
c..   into starevol.par

c..   Read in starevol.par the species in the network and their number code
      open (unit = 99,file = 'starevol_3.31.par')
      read(99,40)
 40   format(32/)
      nsp = 0

 41   nsp = nsp+1
      read(99,110) icode(nsp),species(nsp),anuc(nsp),znuc(nsp)
      echem(int(znuc(nsp)),int(anuc(nsp))) =  icode(nsp)
 110  format(i4,1x,a5,1x,f9.5,1x,f4.0)
      if (icode(nsp).lt.icode(nsp-1)) goto 43
      goto 41

 43   species(nsp) = 'NUTRI'
      icode(nsp) = 0
      anuc(nsp) = 0.d0
      znuc(nsp) = 0.d0
      nsp = nsp+1
      species(nsp) = 'BETAP'
      icode(nsp) = 0
      anuc(nsp) = 0.d0
      znuc(nsp) = 1.d0
      nsp = nsp+1
      species(nsp) = 'BETAM'
      icode(nsp) = 0
      anuc(nsp) = 0.d0
      znuc(nsp) = -1.d0
      nsp = nsp+1
      species(nsp) = 'CAPTE'
      icode(nsp) = 55
      anuc(nsp) = 999.99999d0
      znuc(nsp) = 99.d0
      close(99)

      do i = 1,nre
         do j = 1,nsp
            if (z1(i).eq.species(j)) ireacode(i,1) = icode(j)
            if (z2(i).eq.species(j)) ireacode(i,2) = icode(j)
            if (z3(i).eq.species(j)) ireacode(i,3) = icode(j)
            if (z4(i).eq.species(j)) ireacode(i,4) = icode(j)
         enddo
c         write(*,444) i,reaction(i) ,Qrad(i)+Qnu(i),ireacode(i,1)
c     $        ,ireacode(i,2),ireacode(i,3),ireacode(i,4),Qnu(i)
      enddo
            

c..   Part added to compute the reaction rates for the entire temperature grid.
c..   Density is maintained constant and will be substracted from the final data
c..   For each rate, return to the Nav * <sig*v> / rho value needed for starevol
c..   according to the type of reaction considered (use flag to disentangle)

      iread = 1

      do i = 1,nvit
         T = tnuc9(i)*1.d9
         T8 = tnuc9(i)*10.d0
         call vit(T,rho,pme,iread)

         do j = 1,nre
            if (flag(j).eq.-14.) then
c     electron capture on Be7
               v(j) = v(j) * pme / rho
            else if (flag(j).eq.-13.) then
c     electron capture
               v(j) = v(j) * pme / RHO 
            else if (flag(j).eq.-11.) then
c     photodisintegration or beta-decay
               v(j) = v(j)
            else if (flag(j).eq.-10.) then
c     two-particle reaction
c     !if identical particles: factorials! 
               v(j) = v(j) * f(n1(j))* f(n2(j)) / RHO
            else if (flag(j).eq.-100.) then
c     three-particle reactions
c     !if identical particles: factorials! 
               v(j) = v(j) * f(n1(j))* f(n2(j)) / (RHO * RHO)
            else if (flag(j).eq.-200.) then
c     four-particle reactions
c     !if identical particles: factorials! 
               v(j) = v(j) * f(n1(j))* f(n2(j)) / (RHO * RHO * RHO)
            endif
            if (v(j).gt.1.d-99) then 
               vv(i,j) = log10(v(j))
            else
               vv(i,j) = -99.0d0
            endif
         enddo
c         write(6,*) 
c         write(6,*) tnuc9(i),(vv(i,k),k=1,nre)
      enddo

c..  Create the new list of reactions and associated bibliographic 
c..  references to modify STAREVOL.PAR

c..   Read bibliographic references for the selected network
c..   The "modified" file is created using the reference file given
c..   by NETGEN, in which we have extracted the lines containing
c..   the string " 1  ", which corresponds to the beginning of the line 
c...  containing the reference for the adopted rate

      open(unit=16,file='updated_Li7palphaBe7pgammaNACREIImax.refs')
      do i = 1,nre
         read(16,'(a172)')reference(i)
      enddo


      open(unit=98,file='nuclear_network_LiNACREIImax')
      do i = 1,nre
         write(98,444)sreacnum(i),reaction(i),Qrad(i)+Qnu(i),ireacode(i
     $        ,1),ireacode(i,2),ireacode(i,3),ireacode(i,4),Qnu(i)
      enddo
      write(98,446)
      do i = 1,nre
         write(98,445)sreacnum(i),reference(i)
      enddo
      
c..   Create a new vitsub that reproduces the data in the original
c..   vistub_167c.f for the rates that have not been modified, and  
c..   introduces the new rates for the updated reactions

      open (unit = 10,file = 'vitsub_167e.f')
      open (unit = 20,file = 'vitsub_167d.f',status='old')

      write (10,1000)
      write (10,1200)
      write (10,1000)
      write (10,1100)
      write (10,1300)
      write (10,1100)
      write (10,1000)
      write (10,1400) nvit,nreac
      write (10,1500)

      write (10,1501)
      do ii = 1,17
         iinf = 5*ii-4
         isup = iinf+4
         write (10,1502) (tnuc9(i),i = iinf,isup)
      enddo      
      iinf = 5*ii-4
      isup = iinf+4
      write (10,1503) (tnuc9(i),i = iinf,isup)
      write (10,1000)
      
      read (20,500)
      do i = 1,nvit
         t = tnuc9(i)*1.d9
         read (20,1800) (vr(j),j = 1,7)
         read (20,1800) (vr(j),j = 8,14)
         read (20,1800) (vr(j),j = 15,21)
         read (20,1800) (vr(j),j = 22,28)
         read (20,1800) (vr(j),j = 29,35)
         read (20,1800) (vr(j),j = 36,42)
         read (20,1800) (vr(j),j = 43,49)
         read (20,1800) (vr(j),j = 50,56)
         read (20,1800) (vr(j),j = 57,63)
         read (20,1800) (vr(j),j = 64,70)
         read (20,1800) (vr(j),j = 71,77)
         read (20,1800) (vr(j),j = 78,84)
         read (20,1800) (vr(j),j = 85,91)
         read (20,1800) (vr(j),j = 92,98)
         read (20,1800) (vr(j),j = 99,105)
         read (20,1800) (vr(j),j = 106,112)
         read (20,1800) (vr(j),j = 113,119)
         read (20,1800) (vr(j),j = 120,126)
         read (20,1800) (vr(j),j = 127,133)
         read (20,1800) (vr(j),j = 134,140)
         read (20,1800) (vr(j),j = 141,147)
         read (20,1800) (vr(j),j = 148,154)
         read (20,1800) (vr(j),j = 155,161)
         read (20,1801) (vr(j),j = 162,167)
         do k = 1,nreac
            do j = 1,nre
               if (sreacnum(j).eq.k) then
                  if (i.eq.1) print *,'k',k,vv(i,j),vr(k)
                  vr(k) = vv(i,j)
               endif
            enddo
         enddo
         write (10,1600) t
         write (10,1700) i,i
         write (10,1800) (vr(j),j = 1,7)
         write (10,1800) (vr(j),j = 8,14)
         write (10,1800) (vr(j),j = 15,21)
         write (10,1800) (vr(j),j = 22,28)
         write (10,1800) (vr(j),j = 29,35)
         write (10,1800) (vr(j),j = 36,42)
         write (10,1800) (vr(j),j = 43,49)
         write (10,1800) (vr(j),j = 50,56)
         write (10,1800) (vr(j),j = 57,63)
         write (10,1800) (vr(j),j = 64,70)
         write (10,1800) (vr(j),j = 71,77)
         write (10,1800) (vr(j),j = 78,84)
         write (10,1800) (vr(j),j = 85,91)
         write (10,1800) (vr(j),j = 92,98)
         write (10,1800) (vr(j),j = 99,105)
         write (10,1800) (vr(j),j = 106,112)
         write (10,1800) (vr(j),j = 113,119)
         write (10,1800) (vr(j),j = 120,126)
         write (10,1800) (vr(j),j = 127,133)
         write (10,1800) (vr(j),j = 134,140)
         write (10,1800) (vr(j),j = 141,147)
         write (10,1800) (vr(j),j = 148,154)
         write (10,1800) (vr(j),j = 155,161)
         write (10,1802) (vr(j),j = 162,167)         
      enddo
      write (10,1000)
      write (10,2000)
      write (10,1000)
      
      close(10)
      close(16)
      close(98)
      close(20)



 6003 format(10(1x,e9.4))
 444     format(i4,4x,a37,': Q =',f8.3,' MeV',i7,i4,i4,i4,
     &        '  QINU =',f7.3,' MeV')
 445  format(i5,1x,a172)
 446  format(/,'nuclear network references:',/)
 500  format(///////////////////////////////////////)
 1000 format (1x)
 1100 format (72('*'))
 1200 format (6x,'BLOCK DATA VITSUB')
 1300 format ('*',11x,'tables of nuclear reaction rates in ',
     &     'function of T',10x,'*',/,'*',17x,'(see starevol.par for ',
     &     'the references)',16x,'*')
 1400 format (6x,'implicit none',//,
     &     6x,'integer i,j',/,6x,'integer nvit,nre',/,
     &     6x,'double precision tnuc9,vrate',//,
     &     6x,'parameter (nvit = ',i3,', nre = ',i3,')')
 1500 format (6x,'common /nucvit/ tnuc9(nvit),vrate(nvit,nre)')
 1501 format (1x,/,6x,'data (tnuc9(i),i = 1,nvit) /')
 1502 format (5x,'& ',5(1pd12.6,','))
 1503 format (5x,'& ',4(1pd12.6,','),1pd12.6,'/')
 1600 format (1x,/,'***** nav * <sig*v> / rho  at T = ',d10.4,' K')
 1700 format (6x,'data ((vrate(i,j),i = ',i3,',',i3,'),j = 1,nre) /')
 1800 format (5x,'& ',7(f8.4,','))
 1801 format (5x,'& ',6(f8.4,','),///)
 1802 format (5x,'& ',5(f8.4,','),f8.4,' /')
 1900 format (5x,'& ',2(f8.4,','),f8.4,' /')
c 1900 format (5x,'& ',4(f8.4,','),f8.4,' /')
 2000 format (6x,'end',/)


      end

c============================================================================
      subroutine vit(T,rho,pme,iread)

c     This routine provides the reaction rates v(i) (i=1,nreac) 
c     for temperature T and density RHO 
c     from a linear interpolation of the log of the grid-point rates read
c     from the file 'vit.dat' containing the output from Netgen.
c
c     Note that the factorials accounting for identical particles have
c     already been included in v(i).
c     To obtain the evolution dYj/dt of species j, simply multiply v(i) by
c     the molar fractions of the reacting species j1 + j2: Yj1^n1 Yj2^n2,
c     where n1, n2 are the stoechiometric factors (stored in the vectors
c     with the same name).
c
c     The last lines of subroutine vit print all variables to standard output,
c     to check that it worked properly
c
c     The first call to vit should be done with 
c     iread = 0: to read the vit.dat file, and then store them in array vgrid
c
c     Subsequent calls may be done with
c     iread = 1: to compute rates at other temperatures, after vgrid
c                has been initialized by the first call with iread = 0     
 
c     The parameters 
c          NGRID = maximum number of temperature grid points 
c          NRE   = maximum number of reactions
c
c     may be lowered to better match your network.
c
c     The actual number of grid points and reaction rates are computed
c     automatically, and stored in variables IREAC and KGRID
c
c     PME = mean molecular weight per electron
c     T = temperature in K
c     RHO = density in g/cm3
c
c
c     The rates on grid points are stored in the array 
c         vgrid(ngrid,nre,11)
c     where the last index refers to the density grid for the weak rates
c     dependent upon electron number density:
c     vgrid(ngrid,nre,i) with i=1,11 
c             corresponds to the rate at different electron densities Ne 
c             as given by the variable aNegrid(nre,i)
c     For the other reactions (with flags < 0), 
c     only vgrid(ngrid,nre,1) is relevant 
c
c          
c
c     Coding of reactions:
c
c       n1 to n4: stoechiometric factors
c       z1 to z4: chemical symbol
c       a1 to a4: atomic mass

c     The following table allows to convert proton number into
c     chemical symbol, if required:

c      CHARACTER*2 SYM(0:105)*2
c      DATA SYM/'n ',
c     & 'H ','HE','LI','BE','B ','C ','N ','O ','F ','NE',
c     1 'NA','MG','AL','SI','P ','S ','CL','AR','K ','CA',
c     2 'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',
c     3 'GA','GE','AS','SE','BR','KR','RB','SR','Y ','ZR', 
c     4 'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN',
c     5 'SB','TE','I ','XE','CS','BA','LA','CE','PR','ND',
c     6 'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB',
c     7 'LU','HF','TA','W ','RE','OS','IR','PT','AU','HG',
c     8 'TL','PB','BI','PO','AT','RN','FR','RA','AC','TH',
c     9 'PA','U ','NP','PU','AM','CM','BK','CF','ES','FM',
c     & 'MD','NO','LR','RF','HA'/
c
c     The lines between ======= should be duplicated to the main program 

c     ====================================================================
c     Duplicate the following lines in the main program
c     
      implicit double precision (a-h,o-z)

c.. Actual number of reactions = 57 (see file network_22avril2014_NACRE_NACREII_Illiadis2010_ppchains_CNONeNaMgAl)
c      parameter(ngrid=250,nre=10000)       
      parameter(ngrid=250,nre=2)

      common/cvit/v(NRE),tgrid(ngrid),aNegrid(nre,11)
      common/cvgrid/vgrid(ngrid,nre,11),flag(nre),ndata(nre)
      common/network/n1(nre),n2(nre),n3(nre),n4(nre),
     &               z1(nre),z2(nre),z3(nre),z4(nre),
     &               a1(nre),a2(nre),a3(nre),a4(nre),
     &               reaction(nre),ireac,kgrid
      common/Q/Qrad(nre),Qnu(nre) 

      character z1*6,z2*6,z3*6,z4*6,deut*6,prot*6
      character reaction*37

c     ======================================================================

      character longline(7)*132,dummy*1
      character aflag(11)*8
      dimension nn1(11),nn2(11),nn3(11),nn4(11)
      character zz1(11)*6,zz2(11)*6,zz3(11)*6,zz4(11)*6
      dimension QQrad(11),QQnu(11)
      dimension f(0:4)

c     factorials

      f(0) = 1.
      f(1) = 1.
      f(2) = 2.
      f(3) = 6.
      f(4) = 24.

       T8 = T * 1.E-08

c     ane is the electron density
c     PME is the mean molecular weight per electron 
c         =  [ SUM (XZ/A) ]**(-1)

      ane   = rho * 6.02E+23 / PME

c     if iread = 0, start by reading the data file, then compute the rates

      if (iread.eq.0) then

         open(unit=10,file='Li7palphaBe7pgammaNACREIImax.rates')
       
c     read the data table 
         
         ireac = 0
         do ll = 1,50000
            
c     initialize the string
            do i = 1,7
               do j = 1,132
                  longline(i)(j:j)=' '
               enddo
            enddo
            
c     read the header
            read(10,1100,end=9999) (longline(j),j=1,7)
            ireac = ireac + 1
            
c     l = number of data records on line
            leng = index(longline(5),'         ') 
            if (leng.eq.0) leng = 132
            l = (leng - 11)/11
            ndata(ireac) = l
            
c     write(6,*) longline(5),l,leng
            read(longline(7),1200) (aflag(j),j=1,l)
            aNegrid(ireac,1) = 1.
c     l > 1: flag corresponds to electron densities
c     read electron densities
            if(l.gt.1) then
               read(longline(7),1210) (aNegrid(ireac,j),j=1,l)
               flag(ireac) = aNegrid(ireac,1)
            endif
            
            read(longline(1),1201,end=1) (nn1(j),j=1,l)
 1          read(longline(1),1202,end=2) (zz1(j),j=1,l)
 2          read(longline(2),1201,end=3) (nn2(j),j=1,l)
 3          read(longline(2),1202,end=4) (zz2(j),j=1,l)
 4          read(longline(3),1201,end=5) (nn3(j),j=1,l)
 5          read(longline(3),1202,end=6) (zz3(j),j=1,l)
 6          read(longline(4),1201,end=7) (nn4(j),j=1,l)
 7          read(longline(4),1202,end=8) (zz4(j),j=1,l)
            
 8          read(longline(5),1205,end=11) (QQrad(j),j=1,l)
 11         read(longline(6),1205,err=20,end=12) (QQnu(j), j=1,l)
 12         continue
            
c     read the data
            read(10,*)
            read(10,*)
            
c     check the grid size
            if (ireac.eq.1) then 
               do i = 1,ngrid
                  read(10,1206,err=21) 
     &                 dummy,tgrid(i),(vgrid(i,ireac,j),j=1,l)
                  if (dummy(1:1) .eq. '#') goto 13
                  kgrid = i
               enddo
            else 
               do i = 1,kgrid
                  read(10,1206,end=13,err=21) dummy,
     &                 tgrid(i),(vgrid(i,ireac,j),j=1,l)
               enddo
            endif
            read(10,*)
 13         continue

            n1(ireac) = nn1(1)
            n2(ireac) = nn2(1)
            n3(ireac) = nn3(1)
            n4(ireac) = nn4(1)
            z1(ireac) = zz1(1)
            z2(ireac) = zz2(1)
            z3(ireac) = zz3(1)
            z4(ireac) = zz4(1)
            Qrad(ireac)=QQrad(1)
            Qnu(ireac) =QQnu(1)
            
            if (aflag(1)(5:8).eq.'----') then
               flag(ireac) = -200.
            else if (aflag(1)(6:8).eq.'---') then
               flag(ireac) = -100.
            else if (aflag(1)(7:8).eq.'--') then
               flag(ireac) = -10.
            else if (aflag(1)(6:8).eq.'+++') then
               flag(ireac) = -14.
            else if (aflag(1)(7:8).eq.'++') then
               flag(ireac) = -13.
            else if (aflag(1)(8:8).eq.'+') then
               flag(ireac) = -11.
            endif        
            if (zz1(1)(1:4).eq.'NEUT'.or.zz1(1)(1:4).eq.'PROT') then 
               a1(ireac) = 1.
            else if (zz1(1)(1:5).eq.'OOOOO') then
               a1(ireac) = 0.
            else if (zz1(1)(1:4).eq.'DEUT') then
               a1(ireac) = 2.
            else if (zz1(1)(1:4).eq.'TRIT') then
               a1(ireac) = 3.
            else 
               read (zz1(1)(3:5),'(i3)') i1
               a1(ireac) = float(i1)
            endif
            if (zz2(1)(1:4).eq.'NEUT'.or.zz2(1)(1:4).eq.'PROT') then 
               a2(ireac) = 1.
            else if (zz2(1)(1:5).eq.'OOOOO') then
               a2(ireac) = 0.
            else if (zz2(1)(1:4).eq.'DEUT') then
               a2(ireac) = 2.
            else if (zz2(1)(1:4).eq.'TRIT') then
               a2(ireac) = 3.
            else if (zz2(1)(1:5).eq.'HE  4') then
               a2(ireac) = 4.
            else if (zz2(1)(1:5).eq.'HE  3') then
               a2(ireac) = 3.
            endif
            if (zz3(1)(1:4).eq.'NEUT'.or.zz3(1)(1:4).eq.'PROT') then 
               a3(ireac) = 1.
            else if (zz3(1)(1:5).eq.'OOOOO') then
               a3(ireac) = 0.
            else if (zz3(1)(1:4).eq.'DEUT') then
               a3(ireac) = 2.
            else if (zz3(1)(1:4).eq.'TRIT') then
               a3(ireac) = 3.
            else if (zz3(1)(1:5).eq.'HE  4') then
               a3(ireac) = 4.
            else if (zz3(1)(1:5).eq.'HE  3') then
               a3(ireac) = 3.
            endif
            if (zz4(1)(1:4).eq.'NEUT'.or.zz4(1)(1:4).eq.'PROT') then 
               a4(ireac) = 1.
            else if (zz4(1)(1:5).eq.'OOOOO') then
               a4(ireac) = 0.
            else if (zz4(1)(1:4).eq.'DEUT') then
               a4(ireac) = 2.
            else if (zz4(1)(1:4).eq.'TRIT') then
               a4(ireac) = 3.
            else 
               read(zz4(1)(3:5),'(i3)') i4
               a4(ireac) = float(i4)
            endif
            
c            write(reaction(ireac),1204) 
c     &           nn1(1),zz1(1)(1:6),nn2(1),zz2(1)(1:5),
c     &           nn3(1),zz3(1)(1:5),nn4(1),zz4(1)(1:6)
         enddo

 1100    format(4(a132,//),3(a132,/))
 1101    format (i3,7(a132,/))
 1200    format(14x,11(a8,3x))
 1201    format(14x,11(i1,10x))
 1202    format(14x,11(2x,a6,3x))
 1203    format(14x,11(4x,i3,4x))
c 1204    format(i1,1x,a6,'( ',i1,1x,a5,',',1x,i1,1x,a5,')',2x,i1,1x,a6)
 1205    format(14x,11(f8.3,3x))
 1206    format(a1,2x,f8.4,11(1x,e10.4))
 1210    format(14x,11(e8.1,3x))
         
      endif
      
 9999 continue
c     end of iread=0 loop
      close(10)


c..   Read the network that has been input to the netgen server
      if (iread.eq.0) then
         open (unit=15,file=
     &        'Li7palphaBe7pgammaNACREIImax.network')
         deut = 'H   2'
         prot = 'H   1'
         do i = 1,nre
            read(15,1211) z1(i),z2(i),z3(i),z4(i)
            if (z1(i).eq.deut) z1(i) = 'DEUT '
            if (z1(i).eq.prot) z1(i) = 'PROT '
            if (z2(i).eq.deut) z2(i) = 'DEUT '
            if (z2(i).eq.deut) z2(i) = 'PROT '
            if (z3(i).eq.deut) z3(i) = 'DEUT '
            if (z3(i).eq.deut) z3(i) = 'PROT '
            if (z4(i).eq.deut) z4(i) = 'DEUT '
            if (z4(i).eq.deut) z4(i) = 'PROT '

            write(reaction(i),1204) 
     &           n1(i),z1(i)(1:6),n2(i),z2(i)(1:5),
     &           n3(i),z3(i)(1:5),n4(i),z4(i)(1:6)
         enddo
      endif

 1204 format(i1,1x,a6,'( ',i1,1x,a5,',',1x,i1,1x,a5,')',2x,i1,1x,a6)
 1211 format(11x,a6,4x,a5,4x,a5,5x,a6)

c        interpolate reaction rate

c     locate the position of current T within the grid
c     and compute step, a, b, c and d parameters for cubic spline interpolation


c..   In order to prevent false extrapolation at low temperatures, 
c..   it is necessary to define the vgrid table OUTSIDE its bonds, 
c..   for index O in temperature
      vgrid(0,1:nre,1) = 0.d0

c...  Definition of interpolation bounds
      
      k = 1
      do while (tgrid(k).le.T8.and.k.lt.ngrid)
         k = k+1
      enddo
      klo = k-1
      khi = k
      
      h=tgrid(khi)-tgrid(klo)
      b=(T8 - tgrid(klo))/h

      if (T8.eq.0.0127d0) write(6,*) T8,klo,khi
c      write(17,1003) (tgrid(ll),ll=1,ngrid)
c      write(17,*) h,a,b,c,d

c     perform linear interpolation of log(v)

      do i=1,nre

         if (vgrid(klo,i,1).GT.0.D0 .and. vgrid(khi,i,1).gt.0.D0) then 
            
            if (flag(i).le.0.) then 
c     reaction is not electron-density-dependent beta-decay
               
               v(i)= log10(vgrid(klo,i,1)) + b * (log10(vgrid(khi,i,1)) 
     &              - log10(vgrid(klo,i,1))) 
               
c     test the reaction kind

               if (flag(i).eq.-14.) then
c     electron capture on Be7
                  v(i) = 10.**v(i) * RHO / PME
                  if (v(i).gt.1.51e-7.and.T8.lt.0.01) then
                     v(i) = 1.51e-7
                  endif
               else if (flag(i).eq.-13.) then
c     electron capture
                  v(i) = 10.**v(i) * RHO / PME
               else if (flag(i).eq.-11.) then
c     photodisintegration or beta-decay
                  v(i) = 10.**v(i)
               else if (flag(i).eq.-10.) then
c     two-particle reaction
c     !if identical particles: factorials! 
                  v(i) = 10.**v(i) * RHO /f(n1(i))/f(n2(i))
               else if (flag(i).eq.-100.) then
c     three-particle reactions
c     !if identical particles: factorials! 
                  v(i) = 10.**v(i) * RHO * RHO  /f(n1(i))/f(n2(i))
               else if (flag(i).eq.-200.) then
c     four-particle reactions
c     !if identical particles: factorials! 
                  v(i) = 10.**v(i) * RHO * RHO * RHO / f(n1(i))/f(n2(i))
               endif
               
            else
               
c     beta-decay rate has to be interpolated in density
               
c     locate the position of current rho within the density grid aNegrid
               
               
               jlo = 1
               jhi = ndata(i)
 101           if(jhi-jlo.gt.1) then
                  j = (jhi+jlo)/2
                  if(aNegrid(i,j).gt.ane) then
                     jhi=j
                  else
                     jlo=j
                  endif
                  goto 101
               endif
               deltagrid=aNegrid(i,jhi)-aNegrid(i,jlo)
               
c     linear interpolation in temperature for two density grid points
c     (note: log of beta decay rate is handled)    
               
               vbetalow    = log10(vgrid(klo,i,jlo))    
     &              + b * 
     &              (log10(vgrid(khi,i,jlo))  
     &              - log10(vgrid(klo,i,jlo))) 
               vbetahigh   = log10(vgrid(klo,i,jhi))    
     &              + b * 
     &              (log10(vgrid(khi,i,jhi)) 
     &              - log10(vgrid(klo,i,jhi))) 
               
               vbetalow    = 10.**vbetalow
               vbetahigh   = 10.**vbetahigh
               
c     write(6,*) vbetalow,vbetahigh,jlo,jhi,ndata(i)
               
c     linear interpolation in density 
               
               v(i) = (vbetahigh-vbetalow)/deltagrid
     &              *(ane-aNegrid(i,jhi))
     &              +  vbetahigh
               
c     extrapolation often leads to negative rates 
c     check if so, and then put it to the closest density grid point
c     or to zero
c     (supplementary table values would be needed for those cases)
               
               if(v(i).lt.0.) then
                  if(ane-aNegrid(i,jhi).lt.0.) then
c     extrapolation at lower end of table
                     v(i) = 0.
c     v(i) = vbetalow
                  else
c     extrapolation at upper end of table
                     v(i) = 0.
c     v(i) = vbetahigh
                  endif    
               endif
               
            endif
            
         else
            v(i) = 0.d0
         endif
         
      enddo
      
c     write reaction rates to standard output
      if (iread.eq.0) then
         
         do i = 1,ireac
            write(6,6000) reaction(i),flag(i)
            write(6,6002) n1(i),a1(i),z1(i),
     &                 n2(i),a2(i),z2(i),
     &                 n3(i),a3(i),z3(i),
     &                 n4(i),a4(i),z4(i)
            write(6,6001) Qrad(i),Qnu(i)
            do l=1,ndata(i)
               write(6,6003) (vgrid(k,i,l),k=1,kgrid)
            enddo
            write(6,6004) T8, rho,ane,v(i)
         enddo
      endif


 6000    format(1x,a37,1x,'flag= ',f5.0)
 6001    format(1x,'Qrad (MeV) = ',f7.3,1x,'Qnu (MeV) = ',f7.3)
 6002    format(1x,'i1 = ',i1,1x,'a1 = ',f4.0,1x,'z1 = ',a6,/,
     &          1x,'i2 = ',i1,1x,'a2 = ',f4.0,1x,'z2 = ',a6,/,
     &          1x,'i3 = ',i1,1x,'a3 = ',f4.0,1x,'z3 = ',a6,/,
     &          1x,'i4 = ',i1,1x,'a4 = ',f4.0,1x,'z4 = ',a6)

 6004    format(' T8 = ',f5.2,1x,'rho =',e8.3,1x,
     &          'electron density =',e7.2,'cm-3',/,
     &          ' rate = ',e9.4,/)

 6003    format(11(1x,e9.4))


         return

c        branch here if NaN was detected in Qnu line
 20      continue

         write(6,7000) longline(1)(15:21),longline(2)(15:21),
     &                 longline(3)(15:21),longline(4)(15:21)

 7000    format('** WARNING: NaN have been detected among Qnu values',
     &          ' for reaction ',a7,' + ',a7,' = ',a7,' + ',a7,/,
     &          ' This may be because Qnu is density-dependent.',
     &          ' The Netgen log file provides a complete table',
     &          ' of Qnu(T,rho).',/,' Edit the vit.dat file to remove',
     &          ' NaNs and try again')  

         stop

c        branch here if NaN was detected among V
 21      continue

         write(6,7001) longline(1)(15:21),longline(2)(15:21),
     &                 longline(3)(15:21),longline(4)(15:21)

 7001    format('** WARNING: NaN has been detected in reaction rate',
     &          ' for reaction ',a7,' + ',a7,' = ',a7,' + ',a7,/,
     &          'This may be because your grid extends beyond the'
     &          ' valid data range for this reaction',
     &          ' (see Netgen log file)',/,
     &          ' Either change your temperature or density grid,'
     &          ' or remove NaN in vit.dat by extrapolating the',
     &          ' existing data (risky!)')  

         stop

         end

