
      PROGRAM EVOLTPAGB

************************************************************************
* Automatic detection of thermal pulses along the AGB phase, and       *
* strage of all the pertinent features about them                      *
*           version 2.70                                               *
* improve merge small convective zones, fix bug with pulse definition  *
* new definition for lambda, include dilution factor                   *
************************************************************************

      implicit none

      character*100 seqname
      character*5 chemi

      integer nmod,nspi,npul,no,modseq,nseq1,nseq2

      parameter (nmod = 900000, nspi = 44)
      parameter (npul = 500)

      double precision leff(nmod),reff(nmod),rtot(nmod),teff(nmod),
     &    roeff(nmod),geff(nmod),mloss(nmod),mtot(nmod),dtn(nmod),
     &    time(nmod),tcpu(nmod),ttot
      double precision lh(nmod),lhe(nmod),lc(nmod),lnucl(nmod),
     &     lgrav(nmod)
      double precision emaxHe(nmod),xrHem(nmod),xtHem(nmod),
     &     xroHem(nmod),emaxH(nmod),xrHm(nmod),xtHm(nmod),xroHm(nmod)
      double precision xmHeb(nmod),xtHeb(nmod),xmHem(nmod),xmHet(nmod),
     &     xtHet(nmod),xmHb(nmod),xtHb(nmod),xmHm(nmod),xmHt(nmod),
     &     xtHt(nmod)
      double precision xrHeb(nmod),xroHeb(nmod),xrHet(nmod),
     &     xroHet(nmod),xrHb(nmod),xroHb(nmod),xrHt(nmod),xroHt(nmod)
      double precision xmb(nmod),xrb(nmod),xtb(nmod),xrob(nmod),
     &     xmt(nmod),xrt(nmod),xtt(nmod),xrot(nmod),xmenvb(nmod),
     &     xrenvb(nmod),xtenvb(nmod),xroenvb(nmod)
      double precision xmb2(nmod),xmt2(nmod),xmb3(nmod),xmt3(nmod),
     &     xmb4(nmod),xmt4(nmod),xmb5(nmod),xmt5(nmod),xmb6(nmod),
     &     xmt6(nmod)
      double precision xsps(nmod,nspi),xsp0(nspi),yields(nspi)

      double precision t0,dt0,dum,xmenv1,xmenv2
      double precision lhemax(npul),mtpul(npul),mbpul(npul),
     &     mHtn,dtp,lhmax(npul)
      double precision tbeg(npul),dtpul(npul),dtintpul(npul),
     &     dup_Mb(npul),dup_dm(npul),dmHpul(npul),lambda(npul),
     &     mtotbeg(npul),lbeg(npul),rbeg(npul),dup_Tb(npul),
     &     mlbeg(npul),tbmax(npul),ttmax(npul),robmax(npul),
     &     rotmax(npul),Mcore_pul(npul),dmpulmax(npul),foverlap(npul),
     &     fdilu(npul)

      integer dupH(npul)
      logical icz(nmod),itp(nmod)

      integer lenc,idum
      integer nmods,npuls,iendC,iAGB,iHflash
      integer model(nmod),nphase(nmod)
      integer mcar,nseq(nmod)
      integer ipw
      integer istartAGB,ibegtp,imod,mbeg(npul),mtop(npul),mend(npul),
     &     mdup(npul),mdupcn(npul),mbot(npul)
      integer dnmod
      integer mHa(npul),modbeg(npul),modtop(npul),modend(npul),
     &     moddup(npul)
      integer i,i1,i2,j,m,m1,m2,jj

      dimension chemi(nspi),modseq(nmod)

      data (chemi(j),j = 1,nspi) /
     & 'neutr','proto','deutr','He  3','alpha','Li  6','Li  7','Be  7',
     & 'Be  9','B  10','B  11','C  12','C  13','C  14','N  14','N  15',
     & 'O  15','O  16','O  17','O  18','F  19','Ne 20',
     & 'Ne 21','Ne 22','Na 23','Mg 24','Mg 25','Mg 26','Al26m','Al26g',
     & 'Al 27','Si 28','Si 29','Si 30','P  31','S 32','S  33','S  34',
     & 'S  35','Cl 35','S  36','Cl 36','Cl 37','heavy' /

      write (*,*) 'filename ? '
      read (5,'(a100)') seqname

      open (unit = 10, file = 
     &     seqname(1:lenc(seqname)) // '.hr', status = 'old')
      open (unit = 11, file = 
     &     seqname(1:lenc(seqname)) // '.v1', status = 'old')
      open (unit = 12, file = 
     &     seqname(1:lenc(seqname)) // '.v2', status = 'old')
      open (unit = 13, file = 
     &     seqname(1:lenc(seqname)) // '.v3', status = 'old')
      open (unit = 14, file = 
     &     seqname(1:lenc(seqname)) // '.v4', status = 'old')
      open (unit = 15, file = 
     &     seqname(1:lenc(seqname)) // '.v5', status = 'old')
      open (unit = 16, file = 
     &     seqname(1:lenc(seqname)) // '.v6', status = 'old')
      open (unit = 17, file = 
     &     seqname(1:lenc(seqname)) // '.v11', status = 'old')
      open (unit = 18, file = 
     &     seqname(1:lenc(seqname)) // '.v12', status = 'old')
      open (unit = 21, file = 
     &     seqname(1:lenc(seqname)) // '.s1', status = 'old')
      open (unit = 22, file = 
     &     seqname(1:lenc(seqname)) // '.s2', status = 'old')
      open (unit = 23, file = 
     &     seqname(1:lenc(seqname)) // '.s3', status = 'old')
      open (unit = 24, file = 
     &     seqname(1:lenc(seqname)) // '.s4', status = 'old')

      open (unit = 30, file = 
     &     seqname(1:lenc(seqname)) // '.tpagb0', status = 'unknown')
      open (unit = 31, file = 
     &     seqname(1:lenc(seqname)) // '.tpagb1', status = 'unknown')
      open (unit = 32, file = 
     &     seqname(1:lenc(seqname)) // '.tpagb2', status = 'unknown')
      open (unit = 33, file = 
     &     seqname(1:lenc(seqname)) // '.tpagb3', status = 'unknown')
c      open (unit = 34, file = 
c     &     seqname(1:lenc(seqname)) // '.yields', status = 'unknown')
      open (unit = 35, file = 
     &     seqname(1:lenc(seqname)) // '.out', status = 'unknown')

************************************************************************
***   Reading of the required evolutionary files
************************************************************************

      read (10,990)
      read (11,990)
      read (12,990)
      read (13,990)
      read (14,990)
      read (15,990)
      read (16,990)
      read (17,990)
      read (18,990)
      read (21,990)
      read (22,990)
      read (23,990)
      read (24,990)

      ttot = 0.d0
      do m = 1,nmod
         read (10,*,end = 999) model(m),nphase(m),leff(m),reff(m),
     &        rtot(m),teff(m),roeff(m),geff(m),mloss(m),mtot(m),dtn(m),
     &        time(m),idum,idum,idum,tcpu(m)
         ttot = ttot+tcpu(m)
      enddo
 999  nmods = m-1
      ttot = ttot/3600.d0
      rewind (10)
      read (10,990)
      do m = 1,nmods
         read (10,1003) no,mcar
c         if (m.eq.1.and.mcar.eq.0) stop '1st sequence number not found'
         if (mcar.ne.0) then
            nseq(m) = mcar
            modseq(mcar) = m
         else
            nseq(m) = nseq(m-1)
         endif
      enddo

      print *,'Number of lines : ',nmods
      do m = 1,nmods
         read (12,*) idum,lh(m),lhe(m),lc(m),dum,dum,dum,dum,lnucl(m)
     &        ,lgrav(m)
         read (13,*) idum,xmb(m),xrb(m),xtb(m),xrob(m),xmt(m),xrt(m),
     &        xtt(m),xrot(m),xmenvb(m),xrenvb(m),xtenvb(m),xroenvb(m)
         read (14,*) idum,xmb2(m),xmt2(m),xmb3(m),xmt3(m),xmb4(m),
     &        xmt4(m),xmb5(m),xmt5(m),xmb6(m),xmt6(m)
         read (15,*) idum,xmHb(m),xrHb(m),xtHb(m),xroHb(m),xmHt(m),
     &        xrHt(m),xtHt(m),xroHt(m),xmHm(m)
         read (16,*) idum,xmHeb(m),xrHeb(m),xtHeb(m),xroHeb(m),xmHet(m),
     &        xrHet(m),xtHet(m),xroHet(m),xmHem(m)
         read (21,*) idum,(xsps(m,j),j = 1,11)
         read (22,*) idum,(xsps(m,j),j = 12,22)
         read (23,*) idum,(xsps(m,j),j = 23,33)
         read (24,*) idum,(xsps(m,j),j = 34,44)
      enddo

************************************************************************
***   Yields computation                                            
************************************************************************

c       do j = 1,nspi
c          xsp0(j) = xsps(1,j)
c          yields(j) = 0.d0
c       enddo
c       do m = 1,nmods
c          do j = 1,nspi
c             yields(j) = yields(j)+mloss(m)*(xsps(m,j)-xsp0(j))*dtn(m)
c          enddo
c       enddo

************************************************************************
***   Thermal pulse detection and computation of related quantities
************************************************************************

      iendC = 0
      iAGB = 1
      iHflash = 0
      do m = 1,nmods
c..  merge small convective zones in the pulse region
         if (lgrav(m).le.0.d0) then
            if (xmt4(m).gt.1.d-11.and.xmt4(m).lt.xmHb(m).and.(xmb4(m)-
     &           xmt3(m)).lt.1.d-2) then
               xmt3(m) = xmt4(m)
c               xmb(m) = xmb4(m)
            endif
            if (xmt3(m).gt.1.d-11.and.xmt3(m).lt.xmHb(m).and.(xmb3(m)-
     &           xmt2(m)).lt.1.d-2) then
               xmt2(m) = xmt3(m)
c               xmb(m) = xmb3(m)
            endif
            if (xmt2(m).gt.1.d-11.and.xmt2(m).lt.xmHb(m).and.(xmb2(m)-
     &           xmt(m)).lt.1.d-2) then
               xmt(m) = xmt2(m)
c               xmb(m) = xmb2(m)
            endif
         endif

c..  define presence of pulse
c         itp(m) = lhe(m).gt.1.d3.and.(xtb(m).gt.8.176d0.or.xmb(m)
c     &        .gt.1.d-11).and.xmb(m).lt.xmHet(m)
         itp(m) = xmb(m).lt.xmHet(m).and.xmb(m).gt.0.d0
         icz(m) = .false.
c         if (xmb2(m).gt.1.d-11.and.xmt2(m).lt.xmHb(m).and.lh(m).lt.1.d8)
c     &        icz(m)=.true.
         if (xmb2(m).gt.1.d-11.and.xmt2(m).lt.xmHb(m))  icz(m)=.true.
         if (lc(m).gt.lhe(m).and.lc(m).gt.1.d1) iendC = m         
c         if (lh(m).gt.1.d7.and.nphase(m).gt.4) iHflash = m
      enddo

c..  if superAGB phase
      iAGB = max(iAGB,iendC)
      if (iendC.gt.0) then
         do m = 1,iendC
            itp(m) = .false.
            icz(m)= .false.
         enddo
      endif

      npuls = 0

      istartAGB = nmods+1
      ibegtp = nmods+1
      m = iAGB
c..   find beginning of AGB phase
      if (nphase(m).gt.4) goto 101
      do m = iAGB,nmods-1
         m1 = m-1
         m2 = m+1
         if (nphase(m).gt.4.and.nphase(m1).le.4.and.nphase(m2).gt.4)
     &        goto 101
      enddo
 101  istartAGB = m
      if (istartAGB.ge.(nmods-1)) then
         write (6,*) 'No AGB phase !'
         stop ' 1'
      endif
      imod = istartAGB
 102  do m = imod+1,nmods-1
         m1 = m-1
         m2 = m+1
         if (xmb(m).ge.1.d-11.and.xmb(m1).lt.1.d-11.and.xmb(m2).ge.
     &        1.d-11.and.xmb(m).lt.xmHet(m).and.lh(m).gt.1.d2.and.
     &        lh(m).lt.1.d8) goto 103
c         if (xmb(m).ge.1.d-11.and.xmb(m1).lt.1.d-11.and.xmb(m2).ge.
c     &        1.d-11.and.xmb(m).lt.xmHet(m).and.lh(m).gt.1.d2.and.
c     &        icz(m)) goto 103
      enddo

 103  npuls = npuls+1
      if (npuls.le.1.and.m.ge.(nmods-1)) then
         write (*,*) 'No TP-AGB phase !'
         stop
      endif
      if (npuls.eq.1) then
         ibegtp = m
         t0 = time(ibegtp-10)
      endif
      mbeg(npuls) = m
      do m = mbeg(npuls)+1,nmods-1
         m1 = m-1
         m2 = m+1
c         if (xmb(m).lt.1.d-11.and.xmb(m1).ge.1.d-11.and.xmb(m2).lt.
c     &        1.d-11) goto 104
c         dtp = time(m)-time(mbeg(npuls))
         if (itp(m).and..not.itp(m2)) goto 104
      enddo

 104  mend(npuls) = m
      if (m.ge.(nmods-1).or.m-mbeg(npuls).lt.4) then
         npuls = npuls-1
         goto 105
      endif
      mtpul(npuls) = -1.d99
      mtop(npuls) = -1
      do m = mbeg(npuls)+1,mend(npuls)-1
         if (xmt(m).gt.mtpul(npuls).and..not.icz(m)) then
            mtpul(npuls) = xmt(m)
            mtop(npuls) = m
         endif
      enddo
      mbpul(npuls) = 1.d99
      mbot(npuls) = -1
      do m = mbeg(npuls)+1,mend(npuls)-1
         if (xmb(m).lt.mbpul(npuls).and.xmb(m).gt.1.d-11.and.
     &        .not.icz(m)) then
            mbpul(npuls) = xmb(m)
            mbot(npuls) = m
         endif
      enddo
      lhemax(npuls) = -1.d99
      do m = mbeg(npuls)+1,mend(npuls)-1
         if (lhe(m).gt.lhemax(npuls)) then
            lhemax(npuls) = lhe(m)
            lhmax(npuls) = lh(m)
         endif
      enddo
      if (mend(npuls).lt.(nmods-100)) then
         imod = mend(npuls)+1
         goto 102
      else
         goto 105
      endif

 105  mbeg(npuls+1) = nmods
      do i = 1,npuls
         i1 = i+1
         i2 = i-1

*------------------------
*** pulse characteristics
*------------------------

         modbeg(i) = model(mbeg(i))
         modtop(i) = model(mtop(i))
         modend(i) = model(mend(i))
         tbeg(i) = time(mbeg(i))
c..  stellar parameters at pulse start
         mtotbeg(i) = mtot(mbeg(i))
         lbeg(i) = leff(mbeg(i))
         rbeg(i) = reff(mbeg(i))
         mlbeg(i) = mloss(mbeg(i))
c..  pulse duration
         dtpul(i) = time(mend(i))-time(mbeg(i))
         dnmod = int((mbeg(i1)-mend(i))*0.002)
         if (dnmod.lt.2) dnmod = 2
c..  rho/T at max pulse extension
         tbmax(i) = 10.d0**(xtb(mtop(i)))
         ttmax(i) = 10.d0**(xtt(mtop(i)))
         robmax(i) = 10.d0**(xrob(mtop(i)))
         rotmax(i) = 10.d0**(xrot(mtop(i)))
c.. core mass at the beginning of the pulse
         Mcore_pul(i) = xmHeb(mbeg(i))
         if (Mcore_pul(i).lt.1.d-2) Mcore_pul(i) = xmHb(mbeg(i))
c.. maximum extension of pulse (Mpulse_base,min - M_pulse_top,max)
         dmpulmax(i) = mtpul(i)-mbpul(i)

c.. interpulse duration
         if (i.eq.npuls) then
            dtintpul(i) = 1.d0
         else
            dtintpul(i) = time(mbeg(i1))-time(mend(i))
         endif

*---------------------------
*** envelope characteristics
*---------------------------

c.. determine DUP characteristics
         dup_Mb(i) = 1.d99
         do j = mend(i)+1,mbeg(i1)-1
            if (xmenvb(j).lt.dup_Mb(i)) then
               mdup(i) = j
               moddup(i) = model(mdup(i))
               dup_Mb(i) = xmenvb(j)
               dup_Tb(i) = 10.d0**xtenvb(j)
            endif
         enddo
         do j = mdup(i)+1,mbeg(i1)-1
            if (lh(j).gt.10.d0*lhe(j)) goto 106
         enddo
 106     mdupcn(i) = j
c..  is 3DUP activated ?
         dupH(i) = 0
         if (mtpul(i).gt.dup_Mb(i)) dupH(i) = 1
c..  mass if the 3DUP
         if (dupH(i).eq.1) then
            dup_dm(i) = mtpul(i)-dup_Mb(i)
c            dup_dm(i) = xmenvb((mbeg(i)-dnmod))-dup_Mb(i)
c            dup_dm(i) = xmHm((mbeg(i)-dnmod))-dup_Mb(i)
         else
            dup_dm(i) = 0.d0
         endif

c..  compute lambda
c.. true definition : lambda = delta M_H/delta M_DUP
c c..  reactivation of HBS
c          if (i.gt.1.and.mHa(i).gt.0) then
c            if (i.gt.1) then
c cl              do j = mend(i2)+1,mbeg(i)-1
c               do j = mdup(i2)+1,mbeg(i)-1
c                  if (lh(j).gt.1.d1.and.lh(j-1).lt.1.d1) goto 107
c               enddo
c    107        mHa(i) = j
c            else
c               mHa(i) = -1
c            endif
c            dmHpul(i) = xmHm((mbeg(i)-dnmod))-xmHm(mHa(i))
c            if (dmHpul(i).lt.1.d-10) dmHpul(i) = xmHm((mbeg(i)-dnmod))-
c     &           dup_Mb(i2)
c            lambda(i) = dup_dm(i)/dmHpul(i)
         if (i.gt.1) then
            dmHpul(i) = mtpul(i)-dup_Mb(i2)
         else
            dmHpul(i) = 0.d0
         endif
         if (dupH(i).eq.1) then
            lambda(i) = dup_dm(i)/dmHpul(i)
         else
            lambda(i) = 0.d0
         endif
c..  dilution factor
         if (dupH(i).eq.1) then
            fdilu(i) = dup_dm(i)/dmpulmax(i)
         else
            fdilu(i) = 0.d0
         endif
c..  pulse overlap factor
         if (i.gt.1) then
            foverlap(i) = (min(mtpul(i2),dup_Mb(i2))-mbpul(i))/
     &           dmpulmax(i)
         else
            foverlap(i) = 0.d0
         endif
      enddo
      if(npuls.gt.1) then
         ipw = int(log10(dtintpul(npuls-1)))
      else
         ipw = 0
      endif
      dt0 = 10.d0**(ipw)


************************************************************************
***   Storage of the results
************************************************************************

      write (30,3000) t0,dt0,ttot
      write (31,3100)
      write (32,3200)
      write (33,3300)
      do i = 1,npuls
         nseq1 = nseq(mdup(i))
         xmenv1 = xmenvb(modseq(nseq1))
         nseq2 = nseq1+1
         if (nseq2.gt.mcar) then
            xmenv2 = xmenvb(nmods)
         else
            xmenv2 = xmenvb(modseq(nseq2))
         endif
         if (xmenv2.lt.xmenv1) nseq1 = nseq2
         write (31,3101) i,modbeg(i),nseq(mbeg(i)),modtop(i),
     &        nseq(mtop(i)),modend(i),nseq(mend(i)),moddup(i),
     &        nseq1,tbeg(i),dtpul(i),dtintpul(i)
         write (32,3201) i,mtotbeg(i),lbeg(i),rbeg(i),mlbeg(i),
     &        tbmax(i),ttmax(i),robmax(i),rotmax(i)
         write (33,3301) i,Mcore_pul(i),mbpul(i),lhemax(i),dmpulmax(i),
     &        fdilu(i),foverlap(i),dupH(i),dup_dm(i),lambda(i),
     &        dmHpul(i),dup_Tb(i)
      enddo
c      write (34,3400)
c      do j = 1,nspi
c         write (34,3401) chemi(j),yields(j)
c      enddo
      do i = 1,npuls-1
         mHtn = xmHt(mbeg(i+1)-1)
         if (mHtn.lt.1.d-2) mHtn = xmenvb(mbeg(i+1)-1)
         write (35,3501) i,nseq(mbeg(i))-1,nseq(mend(i))+1,
     &        nseq(mdup(i))+1,nseq(mdupcn(i))+1,nseq(mbeg(i+1))-1,
     &        mbpul(i),mtpul(i),dup_Mb(i),mbpul(i+1),mHtn
      enddo

      close (10)
      close (11)
      close (12)
      close (13)
      close (14)
      close (15)
      close (16)
      close (17)
      close (21)
      close (22)
      close (23)
      close (24)
      close (25)
      close (30)
      close (31)
      close (32)
      close (33)
c      close (34)
      close (35)

 990  format (6(1x,/))

 1000 format (i6,1x,i1,1x,f10.6,1x,f7.5,1x,f7.5,1x,f7.0,1x,
     &     e9.3,1x,e9.2,1x,e10.4,1x,f12.8,1x,e9.3,1x,e16.10,12x,f5.0)
 1001 format (132x,i5)
 1002 format (1x)
c 1003 format (132x,a1)
 1003 format (i7,125x,i5)
 1200 format (9x,e10.4,1x,e10.4)
 1300 format (9x,f10.7,1x,e9.3,1x,f7.4,1x,f8.4,1x,f10.7,1x,e9.3,1x,
     &     f7.4,1x,f8.4,3x,f10.7,1x,e9.3,1x,f9.5,1x,f9.5)
 1400 format (8x,10(1x,f10.7))
 1500 format (8x,2(1x,f9.6,1x,e9.3,1x,e9.3,1x,e9.3),1x,f9.6)
 1600 format (8x,2(1x,f9.6,1x,e9.3,1x,e9.3,1x,e9.3),1x,f9.6)
 2100 format (6x,11(1x,e10.4))
 2500 format (6x,9(1x,e10.4))

 3000 format (2x,1pe15.9,3x,1pe7.1,/,
     &     '#Total integrated CPU time used : ',0pf8.1,' h')
 3100 format (' npul',2x,'mpulse seqpulse',2x,'mtop seqtop',
     &     2x,'mend seqend',3x,'mdup seqdup',5x,'tp_start',
     &     6x,'dtpulse',3x,'dtinter')
 3101 format (1x,i3,3x,4(i7,' ',i4,' ',1x),1pe15.9,1x,0pf9.4,1x,1pe11.5)

 3200 format (' npul',6x,'Mtp',10x,'Ltp',11x,'Rtp',11x,'Mlosstp',9x,
     &     'pulse_Tb',3x,'pulse_Tt',4x,'pulse_rob',4x,'pulse_rot')
 3201 format (1x,i3,2x,0pf11.8,4x,0pf9.2,6x,0pf7.2,8x,1pe10.4,
     &     6x,1pe10.4,1x,1pe10.4,2x,1pe10.4,3x,1pe10.4)

 3300 format (' npul',1x,'pulse_Mcore',3x,'pulse_Mbm',4x,
     &     'pulse_LHem',1x,'pulse_dmmax',3x,'fdilu',2x,
     &     'fovlap',1x,'dupH',3x,'dup_dm',3x,'lambda',3x,
     &     'dmHpul',4x,'dup_Tb')
 3301 format (1x,i3,2x,0pf10.7,3x,0pf10.7,4x,1pe10.4,3x,0pf8.6,4x,
     &     0pf6.4,1x,0pf6.4,3x,i1,3x,0pf8.6,1x,0pf7.4,2x,0pf8.6,2x,
     &     1pe10.4)

c  3100 format (1x,'Total integrated CPU time used : ',0pf8.1,' h',/,/,
c      &     1x,'npul',1x,'mod(beg)[seq]',1x,'mod(top)[seq]',
c      &     1x,'mod(end)[seq]',1x,'mod(dup)[seq]',3x,'time(beg; yr)',
c      &     2x,'delt(pul; yr)',1x,'delt(intpul; yr)')
c  3101 format (1x,i3,3x,4(i7,'[',i4,']',1x),1pe15.9,4x,0pf7.2,6x,1pe11.5)
c  3200 format (1x,'npul',1x,'mtot(beg; sm)',1x,'ltot(beg; sl)',
c      &     1x,'rtot(beg; sr)',1x,'mloss(beg; sm/yr)',2x,'Tbmax (K)',
c      &     2x,'Ttmax (K)',2x,'robmax (cgs)',1x,'rotmax (cgs)')
c  3201 format (1x,i3,2x,0pf11.8,4x,0pf9.2,6x,0pf7.2,8x,1pe10.4,
c      &     6x,1pe10.4,1x,1pe10.4,2x,1pe10.4,3x,1pe10.4)
c  3300 format (1x,'npul',1x,'mc(beg; sm)',1x,'mpul(min; sm)',1x,
c      &     'LHe(max; sl)',1x,'dmpulmax(sm)',1x,'fdilu',2x,
c      &     'fovlap',1x,'dupH',1x,'dup_dm(sm)',1x,'lambda',2x,
c      &     'dmHpul(sm)',1x,'dup_Tb(K)')
c  3301 format (1x,i3,2x,0pf10.7,3x,0pf10.7,4x,1pe10.4,3x,0pf8.6,4x,
c      &     0pf6.4,1x,0pf6.4,3x,l1,3x,0pf8.6,1x,0pf7.4,2x,0pf8.6,2x,
c      &     1pe10.4)

c 3400 format (1x,'Yields (in sm) :')
c 3401 format (1x,a5,' = ',1pe12.5)
 3501 format (1x,i3,1x,5(1x,i4),1x,5(1x,0pf12.10))

      end


************************************************************************

      INTEGER FUNCTION LENC (a)

************************************************************************
* Give last non-blank character position in string a                   *
************************************************************************

      implicit none

      character*(*) a

      integer i,n

      n = len(a)
      do 10 i = n,1,-1
         if (a(i:i).ne.' ') then
            lenc = i
            return
         end if
 10   continue

      lenc = 0

      return
      end
