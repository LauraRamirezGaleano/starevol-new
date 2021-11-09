      
*******************************************************************
      SUBROUTINE myflamespeed (beg,iflag)
*******************************************************************

*     Following Garcia-Berro, Ritossa & Iben ApJ 485:765-784 (1997) we
*     use the words "flame front" to mean the base of the convective
*     shell, we divide the burning zone or "flame" into two parts: a
*     "precursor flame" (ahead of the front) and the main body of the
*     flame (behind the front), we use the distance between the base of
*     the convective shell and the minimum in the luminosity profile in
*     the radiative region ahead of the front as a measure of the width
*     of the precursor flame.

*     beg = .true. at the begining of the time step
*     beg = .false. at the end of the time step

*     At the begining of the time step we compute the position of the
*     convective zones (if any) [stored in variables prefixed with
*     "v"] and the flame front theoretical speed (following Timmes,
*     Woosley & Taam ApJ 420:348-363 (1994))

*     At the end of the time step (once the model has converged) we
*     compute the new position of the convective zones (if still any)
*     and the effective flame front speed during this time step. We
*     therefore make use of the "v" prefixed variables computed at
*     the begining of the time step.
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *


      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.eng'
      include 'evolcom.ion'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.opa'
      include 'evolcom.conv'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      include 'evolcom.flame'

      integer i,j,jmax,jmin,ir,it,ix,ir1,it1,ix1,iflag
      logical beg,expectbcz,expecttcz,lincrease,cbfinished
      double precision rCO,rho,t9,diff

      integer itopCOcore,dnreal,ilntmax,ilum0,ipiclum,nbilum0,iminlum,
     &   ilummincz,iminff,ilimpiclum
      double precision drmi,drma
      double precision pfmwidth,pfrwidth,rff,mff,tempff,roff,xc12ff,
     &   xo16ff,dr,drmin,drmax,drbcz,drtcz,drav,rctime,lummaxcz,
     &   lummax,lntmax,drreal,dmreal,lminlumprof,mminlumprof,mlummaxcz,
     &   lmin,lmax,lummincz,mlummincz,mminff,absb,abst,prod,xc12cbcz,
     &   lnucbcz,lnuccbcz,lgravcbcz,lnurz,lnucrz,lgravrz,taulum
      dimension pfmwidth(maxnbcz),pfrwidth(maxnbcz),rff(maxnbcz),
     &   mff(maxnbcz),tempff(maxnbcz),roff(maxnbcz),xc12ff(maxnbcz),
     &   xo16ff(maxnbcz),drmin(maxnbcz),lminlumprof(maxnbcz),
     &   drmax(maxnbcz),drbcz(maxnbcz),drtcz(maxnbcz),lummaxcz(maxnbcz),
     &   drav(maxnbcz),dnreal(maxnbcz),diff(maxnbcz),drreal(maxnbcz),
     &   dmreal(maxnbcz),mminlumprof(maxnbcz),mlummaxcz(maxnbcz),
     &   mminff(maxnbcz),lummincz(maxnbcz),ilummincz(maxnbcz),
     &   mlummincz(maxnbcz),iminff(maxnbcz),absb(maxnbcz),abst(maxnbcz),
     &   prod(maxnbcz),xc12cbcz(maxnbcz),lnucbcz(maxnbcz),
     &   lnuccbcz(maxnbcz),lgravcbcz(maxnbcz),lnurz(maxnbcz),
     &   lnucrz(maxnbcz),lgravrz(maxnbcz),taulum(maxnbcz)

      integer nrhoflame,nt9flame,nxCOflame
      parameter (nrhoflame=12, nt9flame=6, nxCOflame=3)
      double precision rhoflame,t9flame,xCOflame,flamespeed
      dimension rhoflame(nrhoflame),t9flame(nt9flame),
     &   xCOflame(nxCOflame),flamespeed(nxCOflame,nrhoflame,nt9flame)
      common /twtflames/ rhoflame,t9flame,xCOflame,flamespeed

      logical lprint

      double precision ddm,interp3dlog

      if (.not.(iflag.eq.11.or.iflag.eq.12)) return
      if (iflag.eq.12) then
         lprint = .false.
      else
         lprint = .true.
      endif
      nbcz = 0
      lpreflame = .false.

      if (nphase.lt.6) return

***   Find the upper limit of the CO core (if any)

      do i = nmod,1,-1
         if (xsp(i,ic12).ge.xsp(i,ihe4)) goto 10
      enddo
 10   itopCOcore = i


***   Find if carbon burning has ended or not

      cbfinished = .false.
      if (itopCOcore.gt.0.and.xsp(1,ic12).lt.xsp(1,ine20)) then
         cbfinished = .true.
      endif
      if (cbfinished) return

***   Find the limits of the convective zones in the CO core (if any)

      expectbcz = .true.
      expecttcz = .false.
      do i = 1,itopCOcore
         if (expectbcz.and.crz(i).lt.-1.and.crz(i+1).lt.-1) then
            if (nbcz.ge.maxnbcz) then
               if (lprint) print *,'WARNING: more than ',maxnbcz,
     &              ' convective zones in CO core !'
               goto 20
            endif
            nbcz = nbcz + 1
            ibCOcz(nbcz) = i
            expectbcz = .false.
            expecttcz = .true.
         else if (expecttcz.and.crz(i).gt.0) then
            if (i-1-ibCOcz(nbcz).lt.4) then
               ibCOcz(nbcz) = 0
               nbcz = nbcz - 1
            else
               itCOcz(nbcz) = i-1
            endif
            expectbcz = .true.
            expecttcz = .false.
         endif
      enddo
 20   if (expecttcz) then
         if (lprint) print *, 'Did not find top of convective zone #',
     &        nbcz,'-> considering first ',nbcz-1,' zones instead'
         nbcz = nbcz-1 
      endif

      if (itopCOcore.gt.0) then
         if (nbcz.gt.0.and.lprint) then
            print *,nbcz,' convective zone(s) found in CO core:'
            do i=1,nbcz
               write (nout,111) i,ibCOcz(i),m(ibCOcz(i))/msun,
     &              itCOcz(i),m(itCOcz(i))/msun,itCOcz(i)-ibCOcz(i),
     &              (m(itCOcz(i))-m(ibCOcz(i)))/msun
            enddo
         endif
      endif
      if (lprint) then
         do i=1,nbcz
            if (ibCOcz(i).eq.1) write (nout,100)
         enddo
      endif

***   Look for preflame evidences

      if (itopCOcore.gt.0.and.nbcz.eq.0) then

*     Look for shell where temperature is at maximum

         ilntmax = 0
         lntmax = 0.d0
         do i=itopCOcore,1,-1
            if (lnt(i).gt.lntmax) then
               lntmax = lnt(i)
               ilntmax = i
            endif
         enddo
         if (ilntmax.eq.0.and.lprint) then
            print *,'Did not find max temperature between shell #1 ',
     &         'and shell #',itopCOcore
         endif

*     Look for lum=0

         nbilum0 = 0
         ilum0 = 0
         do i=itopCOcore-1,1,-1
            if (lum(i+1).ge.0.d0.and.lum(i).lt.0.d0) then
               if (lprint) write (nout,001) i,i+1,lum(i),lum(i+1),
     &              m(i)/msun,r(i)/1.d5
               ilum0 = i
               nbilum0 = nbilum0+1
            endif
         enddo
         if (lprint) then
            if (nbilum0.gt.1) then
               print *,'More than one luminosity sign change'
            else if  (nbilum0.eq.1.and.ilum0.lt.ilntmax) then
               print *,'WARNING: lum=0 not above max temp point ',
     &              '(shell#',ilntmax,')'
            endif
         endif

         if (nbilum0.eq.1) then

*     Look for a luminosity pic above max temperature shell
*     La technique utilisee ci-dessous pour trouver la position du gros
*     pic de luminosite au dessus de ilum0 ne marche pas si il existe du
*     bruit numerique qui fait apparaitre entre ilum0 et le gros pic
*     (genre trois points...) un petit pic tel que lum(ilimpiclum_{petit
*     pic}) < lum(ilimpiclum_{gros pic})

            ipiclum = 0
            ilimpiclum = 0
            height = 0.d0
            do i=ilum0+1,itopCOcore-1
               if (lum(i+1).lt.lum(i)) goto 25
            enddo
  25        if (i.eq.ilum0+1) then
               if (lprint) print *,'Luminosity not increasing above ',
     &              'lum=0 point'
            else if (i.eq.itopCOcore) then
               if (lprint) write (nout,003)
            else
               ipiclum = i
c               do i=ipiclum,itopCOcore-1
c                  if (lum(i+1).gt.lum(i)) goto 28
c               enddo               
c  28           if (i.eq.itopCOcore) then
c                  write (nout,003)
c               else
c                  ilimpiclum = i
c                  height = lum(ipiclum)/lum(ilimpiclum)
c                     write (nout,002) ilum0+1,ipiclum,ilimpiclum,height,
c     &                  (m(ilimpiclum)-m(ilum0+1))/msun,
c     &                  (r(ilimpiclum)-r(ilum0+1))/1.d5
c                endif
               ilimpiclum = ipiclum
               do i=ipiclum+1,itopCOcore
                  if (lum(i).lt.lum(ilimpiclum)) ilimpiclum = i
               enddo  
               do i=ilimpiclum-1,ilum0,-1
                  if (lum(i).gt.lum(ipiclum)) ipiclum = i
               enddo
               height = lum(ipiclum)/lum(ilimpiclum)
               if (lprint) then
                  if (ilimpiclum.gt.(ipiclum+1).and.height.gt.1.1d0)
     &                 then
                     write (nout,002) ilum0+1,ipiclum,ilimpiclum,height,
     &                    (m(ilimpiclum)-m(ilum0+1))/msun,
     &                    (r(ilimpiclum)-r(ilum0+1))/1.d5
                  else
                     write (nout,005) ilum0+1,ipiclum,ilimpiclum,height
                  endif
               endif
            endif

*     Look for a minimum in the luminosity profile under lum=0

c            iminlum = 0
c            do i=ilum0,2,-1
c               if (lum(i-1).gt.lum(i)) goto 29
c            enddo
c  29        if (i.eq.ilum0) then
c               print *,'Luminosity not decreasing under lum=0 point'
c            else if (i.eq.1) then
c                  print *,'No minimum in the luminosity profile under ',
c     &            'lum=0 point'
c            else
c               iminlum = i
c               lummax = -1.d99
c               drpreflame = r(ilum0)-r(iminlum)
c               do j=1,nmod
c                  if (lum(j).gt.lummax) then
c                     lummax = lum(j)
c                  endif
c               enddo
c                  write (nout,006) iminlum,lum(iminlum),
c     &            lum(iminlum)/lummax,(m(ilum0)-m(iminlum))/msun,
c     &            drpreflame/1.d5,ilum0-iminlum,lummax
c            endif

            iminlum = ilum0
            do i=ilum0-1,1,-1
               if (lum(i).lt.lum(iminlum)) iminlum = i
            enddo
            if (iminlum.lt.ilum0) then
               lummax = -1.d99
               drpreflame = r(ilum0)-r(iminlum)
               do j=1,nmod
                  if (lum(j).gt.lummax) then
                     lummax = lum(j)
                  endif
               enddo
               if (lprint) write (nout,006) iminlum,lum(iminlum),
     &            lum(iminlum)/lummax,(m(ilum0)-m(iminlum))/msun,
     &            drpreflame/1.d5,ilum0-iminlum,lummax
            elseif (lprint) then
               print *,'No minimum in the luminosity profile under ',
     &               'lum = 0 point'
            endif

*     Discretization diagnostic for future precursor flame (if any)

            if (iminlum.gt.0.and.iminlum.lt.ilum0) then
               drmi = 1.d99
               drma = 0.d0
               do i=iminlum,ilum0-1
                  if ((r(i+1)-r(i)).gt.drma) then
                     drma = r(i+1)-r(i)
                  else if ((r(i+1)-r(i)).lt.drmi) then
                     drmi = r(i+1)-r(i)
                  endif
               enddo
               if (lprint) write (nout,007) drmi/1.d5,drma/1.d5,
     &               (r(ilum0)-r(ilum0-1))/1.d5
            endif         
         endif

*     Warn user for pre-flame detection, if any

c         if (nbilum0.eq.1.and.ipiclum.gt.0.and.ilimpiclum.gt.0.and.
c     &      iminlum.gt.0) then
         if (nbilum0.eq.1
     &      .and.ipiclum.gt.0.and.ilimpiclum.gt.(ipiclum+1)
     &      .and.height.gt.1.1d0
     &      .and.iminlum.gt.0.and.iminlum.lt.ilum0) then
c            lpreflame = .true.
            if (lprint) then
               if (lpreflame) then
                  write (nout,004)
               else
                  print *,'***** WARNING: preflame detected but ',
     &                 'lpreflame set to .false. *****'
               endif
            endif
         endif
      endif

***   Get the flames properties

      do i=1,nbcz

*     Flame front position, temperature, density and 
*     chemical composition
         rff(i) = r(ibCOcz(i))
         mff(i) = m(ibCOcz(i))
c        WARNING: at the end of mesh, t(i) is not updated: we have to use exp(lnt(i)) !!
         if (ibCOcz(i).gt.1) then
            tempff(i) = exp(lnt(ibCOcz(i)-1))
            roff(i) = ro(ibCOcz(i)-1)
            xc12ff(i) = xsp(ibCOcz(i)-1,ic12)
            xo16ff(i) = xsp(ibCOcz(i)-1,io16)
         else
            tempff(i) = exp(lnt(ibCOcz(i)))
            roff(i) = ro(ibCOcz(i))
            xc12ff(i) = xsp(ibCOcz(i),ic12)
            xo16ff(i) = xsp(ibCOcz(i),io16)
         endif
         if (lprint) then
            print *,'Flame front #',i,':'
            print *,'   r           =',rff(i)/1.d5,' km'
            print *,'   m           =',mff(i)/msun,' Msun'
            print *,'   T8          =',tempff(i)/1.d8
            print *,'   rho6        =',roff(i)/1.d6
            print *,'   xC12        =',xc12ff(i)
            print *,'   xO16        =',xo16ff(i)
            print *,'   xC12 + xO16 =',xc12ff(i)+xo16ff(i)
         endif

*     Minimum in the luminosity profile under flame front
         iminlumprof(i) = 0
         lminlumprof(i) = 0.d0
         mminlumprof(i) = 0.d0
         if (ibCOcz(i).gt.1) then
            jmax = ibCOcz(i)
            if (i.eq.1) then
               jmin = 1
            else
               jmin = itCOcz(i-1)
            endif
            do j=jmax-1,jmin,-1
               if (lum(j).gt.lum(j+1)) goto 30
            enddo
   30       if (j.eq.(jmin-1)) then
               if (lprint) print *,'WARNING: Did not find the minimum ',
     &              'of the luminosity profile under flame front #',i
            else if (j.eq.(jmax-1)) then
               if (lprint) print *,'WARNING: Luminosity not decreasing',
     &              ' under flame front #',i
            else
               if (lum(j).gt.0.d0) then
                  if (lprint) print *, 'WARNING: Minimum of the ',
     &                 'luminosity profile under flame front #',i,
     &                 'is not <0 ! min(lum) = ',lum(j)
               else
                  iminlumprof(i) = j
                  lminlumprof(i) = lum(j)
                  mminlumprof(i) = m(j)
               endif
            endif
         endif

         iminff(i) = 0
         if (iminlumprof(i).gt.0) then
            jmax = iminlumprof(i)
            if (i.eq.1) then
               jmin = 1
            else
               jmin = itCOcz(i-1)
            endif
            do j=jmax,jmin,-1
               if (abs(lum(j)).lt.abs(lum(iminlumprof(i))/1.d3)) then
                  iminff(i) = j
                  mminff(i) = m(j)
                  if (lprint) print *,' lum=lum(iminlumprof)/1.d3 at',
     &                 ' shell # ',iminff(i),' , m = ',mminff(i)/msun
                  goto 35
               endif
            enddo
   35       if (iminff(i).eq.0.and.lprint) then
               print *,' Did not find lum=lum(iminlumprof)/1.d3 under ',
     &         'iminlumprof'
            endif
         endif

*     Precursor flame width
         if (iminlumprof(i).gt.0) then
            pfmwidth(i) = m(ibCOcz(i)) - m(iminlumprof(i))
            pfrwidth(i) = r(ibCOcz(i)) - r(iminlumprof(i))
            if (lprint)  then
               print *,'Width of precursor flame #',i,':'
               print *,'   di = ',ibCOcz(i) - iminlumprof(i),' shells'
               print *,'   dm = ',pfmwidth(i)/msun,' Msun'
               print *,'   dr = ',pfrwidth(i)/1.d5,' km'
               print *,'Discretization diagnostic for precursor flame #'
     &              ,i,':'
            endif
            drmin(i) = 1.d99
            drmax(i) = 0.d0
            do j= iminlumprof(i),ibCOcz(i)-1
               dr = r(j+1)-r(j)
               if (dr.gt.drmax(i)) drmax(i) = dr
               if (dr.lt.drmin(i)) drmin(i) = dr
            enddo
            drmax(i) = drmax(i)/1.d5
            drmin(i) = drmin(i)/1.d5
            drav(i) = (pfrwidth(i)/(ibCOcz(i)-iminlumprof(i)))/1.d5
            drbcz(i) = (r(ibCOcz(i))-r(ibCOcz(i)-1))/1.d5
            if (lprint)  then
               print *,'   dr_bff = ', drbcz(i),' km'
               print *,'   dr_min = ', drmin(i),' km'
               print *,'   dr_max = ', drmax(i),' km'
               print *,'   <dr>   = ', drav(i) ,' km'
            endif
         else
            pfmwidth(i) = 0
            pfrwidth(i) = 0
            drmin(i) = 0.d0
            drmax(i) = 0.d0
            drav(i) = 0.d0
            if (ibCOcz(i).gt.1) then
               drbcz(i) = (r(ibCOcz(i))-r(ibCOcz(i)-1))/1.d5
               if (lprint) then
                  print *,'   dr_bff = ', drbcz(i),' km'
                  print *,'Could not determine precursor flame width'
               endif
            else
               drbcz(i) = 0.d0
            endif
         endif

*     Characteristics of top of convective zone

         drtcz(i) = (r(itCOcz(i))-r(itCOcz(i)-1))/1.d5
         if (lprint) print *,'Convective zone #',i,':'
         if (lprint) print *,'   dr_tcz = ', drtcz(i),' km'
         if (lum(ibCOcz(i)+1).gt.lum(ibCOcz(i))) then
            lincrease = .true.
            if (lprint) print *,'   Luminosity increases above bcz ',
     &           '(shell #',ibCOcz(i),')'
         else
            lincrease = .false.
            if (lprint) print *,'   Luminosity decreases above bcz ',
     &           '(shell #',ibCOcz(i),')'
         endif
         do j=ibCOcz(i)+1,itCOcz(i)-1
            if (lum(j+1).lt.lum(j).and.lincrease) then
               if (lprint) print *
     &              ,'   Luminosity decreases above shell #',j
               lincrease = .false.
            else if (lum(j+1).gt.lum(j).and..not.lincrease) then
               if (lprint) print *
     &              ,'   Luminosity increases above shell #',j
               lincrease = .true.
            endif
         enddo
         if (lprint) print *,'   Top of convective zone at shell #'
     &        ,itCOcz(i)
         ilummaxcz(i) = 0
         lmax = -1.d99
         do j=ibCOcz(i),itCOcz(i)
            if (lum(j).gt.lmax) then
               lmax = lum(j)
               lummaxcz(i) = lmax
               ilummaxcz(i) = j
               mlummaxcz(i) = m(j)
            endif
         enddo
         lummax = -1.d99
         do j=1,nmod
            if (lum(j).gt.lummax) then
               lummax = lum(j)
            endif
         enddo
         if (lprint) then
            print *,'   Max luminosity in this star: lum(max*) = ',
     &           lummax
            print *,'   Luminosity at top of this cz: shell #',
     &           itCOcz(i),' lum = ',lum(itCOcz(i))
            print *,'   lum/lum(max*) = ',lum(itCOcz(i))/lummax
            print *,'   Max luminosity in this cz: shell #',
     &           ilummaxcz(i),' lummaxcz = ',lummaxcz(i)
            print *,'   |lminlumprof/lummaxcz| = ',
     &           abs(lminlumprof(i)/lummaxcz(i)),
     &           ' [lminlumprof = ',lminlumprof(i),' ]'
         endif
         ilummincz(i) = 0
         if (ilummaxcz(i).gt.0) then
            lmin = 1.d99
            jmin = ilummaxcz(i)
            if (i.lt.nbcz) then
               jmax = ibCOcz(i+1)
            else
               jmax = itopCOcore
            endif
            do j=jmin,jmax
               if (lum(j).lt.lmin) then
                  lmin = lum(j)
                  lummincz(i) = lmin
                  ilummincz(i) = j
                  mlummincz(i) = m(j)
               endif
            enddo
         endif

         absb(i) = 0.d0
         abst(i) = 0.d0
         prod(i) = 0.d0
         if (iminff(i).gt.0.and.ilummincz(i).gt.0) then
            absb(i) = abs(lminlumprof(i))/(mminlumprof(i)-mminff(i))
            abst(i) = (lummaxcz(i)-lummincz(i))/
     &         (mlummincz(i)-mlummaxcz(i))
            prod(i) = (lummaxcz(i)-lminlumprof(i))/
     &         (mlummaxcz(i)-mminlumprof(i))
            if (lprint) print *,' absb = ',absb(i),' abst =
     &           ',abst(i) ,' prod = ',prod(i)
         endif

*     Computes XC12m, Lnu, Lgrav and carbon luminosity 
*     inside carbon burning convective zone

         xc12cbcz(i) = 0.d0 
         lnucbcz(i) = 0.d0
         lnuccbcz(i) = 0.d0
         lgravcbcz(i) = 0.d0
         ddm = 1.d0
         if (modeli.gt.0.and.model.gt.modeli) then
            do j=ibCOcz(i),itCOcz(i)
               ddm = ddm+dm(j)
               xc12cbcz(i) = xc12cbcz(i)+xsp(j,ic12)*dm(j)
               lnucbcz(i) = lnucbcz(i)+enupla(j)*dm(j)
               lnuccbcz(i) = lnuccbcz(i)+(enucl(j)+enupla(j))*dm(j)
               lgravcbcz(i) = lgravcbcz(i)+(egrav(j)*dm(j))
            enddo
         endif
         xc12cbcz(i) = xc12cbcz(i)/ddm
         if (lprint) then            
            print *,'Carbon burning convective zone #',i,':'
            print *,'  CZ XC12    = ',xc12cbcz(i)
            print *,'  CZ Lnu_pla = ',lnucbcz(i)/lsun,' Lsun'
            print *,'  CZ LC      = ',lnuccbcz(i)/lsun,' Lsun'
            print *,'  CZ Lgrav   = ',lgravcbcz(i)/lsun,' Lsun'
         endif

*     Lnu, Lgrav and carbon luminosity under carbon burning
*     convective zone

         lnurz(i) = 0.d0
         lnucrz(i) = 0.d0
         lgravrz(i) = 0.d0
         jmax = ibCOcz(i)
         if (i.eq.1) then
            jmin = 1
         else
            jmin = itCOcz(i-1)
         endif
         if (modeli.gt.0.and.model.gt.modeli) then
            do j=jmin,jmax
               lnurz(i) = lnurz(i)+enupla(j)*dm(j)
               lnucrz(i) = lnucrz(i)+(enucl(j)+enupla(j))*dm(j)
               lgravrz(i) = lgravrz(i)+(egrav(j)*dm(j))
            enddo
         endif
         if (lprint) then
            print *,'Under carbon burning convective zone #',i,':'
            print *,'  RZ Lnu_pla = ',lnurz(i)/lsun,' Lsun'
            print *,'  RZ LC      = ',lnucrz(i)/lsun,' Lsun'
            print *,'  RZ Lgrav   = ',lgravrz(i)/lsun,' Lsun'
         endif
      enddo

***   Different behavior depending on wether it is the begining or 
*     the end of the time step


*********************************
*     Beginning of the time step
      if (beg) then
*********************************

*     Compute theoretical flame front speed

         do i=1,nbcz

            urff(i) = 1.d-99
            thurff(i) = 0.d0
            myerrx(i) = ' '
            myerrt(i) = ' '
            myerrr(i) = ' '
            if (ibCOcz(i).gt.1) then
               if (xO16ff(i).gt.0.d0) then
                  rCO = xc12ff(i)/xO16ff(i)
               else
                  if (lprint) then
                     print *,'ERROR: unexpected value for ',
     &                    'X(O16) at flame front #',i,': X(O16)=',
     &                    xO16ff(i)
                     print *,'Can t compute theoretical flame speed'
                  endif
                  goto 99
               endif
               do j=1,nxCOflame
                  if (rCO.lt.xCOflame(j)) goto 40
               enddo
   40          if (j.eq.1) then
                  if (lprint) write (nout,550) i,rCO,xCOflame(1)
                  ix = 1
                  ix1 = 2
                  myerrx(i) = 'x'
               else if (j.eq.(nxCOflame+1)) then
                  if (lprint) write (nout,551) i,rCO,xCOflame(nxCOflame)
                  ix = nxCOflame-1
                  ix1 = nxCOflame
                  myerrx(i) = 'X'
               else
                  ix = j-1
                  ix1 = j
                  myerrx(i) = ' '
               endif

               rho = roff(i)
               do j=1,nrhoflame
                  if (rho.lt.rhoflame(j)) goto 50
               enddo
   50          if (j.eq.1) then
                  if (lprint) write (nout,552) i,rho,rhoflame(1)
                  ir = 1
                  ir1 = 2
                  myerrr(i) = 'r'
               else if (j.eq.(nrhoflame+1)) then
                  if (lprint) write (nout,553) i,rho,rhoflame(nrhoflame)
                  ir = nrhoflame-1
                  ir1 = nrhoflame
                  myerrr(i) = 'R'
               else
                  ir = j-1
                  ir1 = j
                  myerrr(i) = ' '
               endif

               t9 = tempff(i)/1.d9
               do j=1,nt9flame
                  if (t9.lt.t9flame(j)) goto 60
               enddo
   60          if (j.eq.1) then
                  if (lprint) write (nout,554) i,t9,t9flame(1)
                  it = 1
                  it1 = 2
                  myerrt(i) = 't'
               else if (j.eq.(nt9flame+1)) then
                  if (lprint) write (nout,554) i,t9,t9flame(nt9flame)
                  it = nt9flame-1
                  it1 = nt9flame
                  myerrt(i) = 'T'
               else
                  it = j-1
                  it1 = j
                  myerrt(i) = ' '
               endif


c log interpolation enables to avoid thurff(i)<0...
c [furthermore, the plot log[v_theo(rho,T9)] looks much smoother
c than the plot v_theo(rho,T9).]
               thurff(i) = interp3dlog (rCO,rho,t9,
     &               xCOflame(ix),xCOflame(ix1),
     &               rhoflame(ir),rhoflame(ir1),
     &               t9flame(it),t9flame(it1),
     &               flamespeed(ix,ir,it),flamespeed(ix,ir1,it),
     &               flamespeed(ix,ir,it1),flamespeed(ix,ir1,it1),
     &               flamespeed(ix1,ir,it),flamespeed(ix1,ir1,it),
     &               flamespeed(ix1,ir,it1),flamespeed(ix1,ir1,it1))

               if (lprint) write (nout,666) i,-thurff(i)
            endif
         enddo

*     Compute theoretical flame front advance

         do i=1,nbcz
            dntheo(i) = 0
            dmtheo(i) = 0.d0
            if (ibCOcz(i).gt.1) then
               drtheo(i) = (thurff(i)*dtn)/1.d5
               if (drtheo(i).gt.0.d0) then
                  if (drtheo(i).lt.drbcz(i)) then
                     if (lprint) write (nout,333) i,drtheo(i),drbcz(i)
                  else
                     jmax = ibCOcz(i)
                     if (i.eq.1) then
                        jmin = 1
                     else
                        jmin = itCOcz(i-1)
                     endif
                     do j=jmax-1,jmin,-1
                        if ((r(jmax)-r(j))/1.d5.gt.drtheo(i)) goto 90
                     enddo
   90                dntheo(i) = (j+1)-jmax
                     if (lprint) write (nout,334) i,drtheo(i),drbcz(i)
     &                    ,dntheo(i)
                     dmtheo(i) = m(ibCOcz(i))-m(j+1)
                  endif
               else
                  if (lprint) write (nout,*) ' -> Flame front #',i,
     &                 ' speed > 0, not supposed to go to the center'
               endif
            endif
         enddo

*     Save the values

   99    vnbcz = nbcz
         do i=1,nbcz
            vtempff(i) = tempff(i)
            vrff(i) = rff(i)
            vmff(i) = mff(i)
            vibCOcz(i) = ibCOcz(i)
            vpfrwidth(i) = pfrwidth(i)
         enddo

*********************************
*     End of the time step
      else
*********************************

         if (model.eq.modeli+1) write (40,881)

         taulum(1) = 0.d0
         dnreal(1) = 0
         drreal(1) = 0.d0
         urff(1) = 0.d0
         umff(1) = 0.d0
         diff(1) = 0.d0

         if (vnbcz.ne.nbcz) then

            if (nbcz.gt.0.and.lprint) then
               print *,'WARNING: Not as many convective zones as ',
     &           'during previous time step -> ',
     &           'No computation of flames speeds'
            endif

c            do i=1,nbcz
c               do j=1,vnbcz
c                  if (rff(i).gt.(vrff(j)-thurff(j)*dtn).and.
c     &                  rff(i).lt.(vrff(j)+thurff(j)*dtn)) then
c                     print *,'New ff #',i,' seems to correspond to ',
c     &                     'previous ff #',j
c                  endif
c               enddo
c            enddo

         else

            do i=1,nbcz
               taulum(i) = dtn*abs(lum(ibCOCz(i)-1))/
     &               abs(lum(ibCOCz(i)-1)-vlum(vibCOCz(i)-1))
               if (lprint) print *,'Flame front #',i,' taulum = '
     &              ,taulum(i)*seci,' yr'
            enddo

            do i=1,nbcz
               if (vibCOcz(i).gt.1) then
                  dnreal(i) = ibCOcz(i)-vibCOcz(i)
                  if (lprint) write (nout,444) i,dnreal(i),dntheo(i)
*     Compare flame advance with pf width
                  drreal(i) = rff(i)-vrff(i)
                  if (lprint) write (nout,445) i,drreal(i)/1.d5,
     &                 -drtheo(i)
                  dmreal(i) = mff(i)-vmff(i)
                  if (lprint) write (nout,446) i,dmreal(i)/msun,
     &                 -dmtheo(i)
     &                 /msun
                  if (drreal(i).lt.0) then
                     if (vpfrwidth(i).gt.0.d0.and.lprint) then
                        print *,'ff r advance/vpfrwidth = ',
     &                     drreal(i)/vpfrwidth(i)
                     endif
                     if (pfrwidth(i).gt.0.d0.and.lprint) then
                        print *,'ff r advance/pfrwidth = ',
     &                     drreal(i)/pfrwidth(i)
                     endif
                  endif
               endif
            enddo

            do i=1,nbcz
               if (vibCOcz(i).gt.1) then
*     Compute flame front speed
                  urff(i) = (rff(i)-vrff(i))/dtn
                  umff(i) = (mff(i)-vmff(i))/dtn
*     Compare it with theoretical flame front speed
                  if (urff(i).lt.0.d0.and.thurff(i).gt.0.d0) then
                     diff(i) = thurff(i)/(-urff(i))
                  endif
                  if (lprint) write (nout,222) i,umff(i)/(msun*seci),
     &                  urff(i),-thurff(i),diff(i)
                  if (((mff(i)-vmff(i))*msun).gt.5d-2) then
                     if (lprint) print *,'   -> Seems strange...'
                  endif
*     Compute expected time to reach the center
                  if (urff(i).lt.0.d0) then
                     rctime = abs(rff(i)/urff(i))
                     if (lprint) write (nout,777) rctime*seci,int(rctime
     &                    /dtn)
                  endif
               endif
            enddo

         endif

*     Store results

         if (nbcz.gt.0) then
            write (40,888) model,time*seci,dtn*seci,nbcz,
     &         m(ibCOcz(1))/msun,m(itCOcz(1))/msun,
     &         (m(itCOcz(1))-m(ibCOcz(1)))/msun,
     &         iminlumprof(1),
     &         ibCOcz(1),itCOcz(1),itCOcz(1)-ibCOcz(1),
     &         rff(1)/1.d5,mff(1)/msun,tempff(1)/1.d8,roff(1)/1.d6,
     &         xc12ff(1),xo16ff(1),xc12ff(1)+xo16ff(1),
     &         pfmwidth(1)/msun,pfrwidth(1)/1.d5,
     &         drmin(1),drmax(1),drav(1),drbcz(1),drtheo(1),
     &         -thurff(1),urff(1),diff(1),lminlumprof(1)/lsun,
     &         lummaxcz(1)/lsun,
     &         dntheo(1),dnreal(1),xc12cbcz(1),lnucbcz(1)/lsun,
     &         lnuccbcz(1)/lsun,lgravcbcz(1)/lsun,lnurz(1)/lsun,
     &         lnucrz(1)/lsun,lgravrz(1)/lsun,taulum(1)*seci
         endif

      endif

      if (lprint) write (nout,*)

  001 format(1x,'Luminosity sign changes between shell #',i4,
     &   ' and shell #',i4,': lum = ',1pd11.4,' -> lum = ',1pd11.4,
     &   ' [ m = ',1pd10.4,' Msun, r = ',1pd10.4,' km ]') 
  002 format(1x,'Luminosity pic: shells #',i4,' - ',i4,' - ',i4,
     &   ' [ height = ',1pd8.2,', dm = ',1pd10.4,' Msun, dr = ',1pd10.4,
     &   ' km ]')
  003 format(1x,'No luminosity pic above lum=0 point')
  004 format(
     &   1x,'**************************************',/,
     &   1x,'!!! WARNING !!! Pre-flame detected !!!',/,
     &   1x,'**************************************')
  005 format(1x,'Shells #',i4,' - ',i4,' - ',i4,', Height = ',
     &   1pd8.2,' -> Does not look like a true luminosity pic')
  006 format(1x,'Minimum in the luminosity profile under ',
     &   'lum=0 point is at shell #',i4,' lum = ',1pd11.4,
     &   ' lum/lummax = ',1pd11.4,' dmpfw = ',1pd10.4,
     &   ' Msun  drpfw = ',1pd10.4,' Km  dipfw = ',i4,
     &   ' [ Max luminosity in this star: lummax = ',1pd11.4,' ]')
  007 format(1x,'Disc diag for future pf: ',
     &   'drmin = ',1pd11.4,' km - drmax = ',1pd11.4,
     &   ' km - dr_lum0 = ',1pd11.4,' km')
  100 format(
     &   1x,'****************************',/,
     &   1x,'!!! Convection at center !!!',/,
     &   1x,'****************************')
  111 format(3x,i4,':',i4,' (',1pd12.6,' Msun) - ',i4,
     &   ' (',1pd12.6,' Msun) -- width: ',i4,' shells (',1pd12.6,
     &   ' Msun)')         
  222 format(3x,'Speed of flame front #',i4,' : v_m = ',1pd11.4,
     &   ' Msun/yr   v_r = ',1pd9.2,
     &   ' cm/s',' (v_theo = ',1pd9.2,' cm/s, |v_theo/v_r| = ',
     &   1pd11.4,')')
  333 format(1x,'Theoretical advance of flame front #',i4,': ',1pd11.4,
     &   ' km... smaller than width of shell ',
     &   'under flame front (',1pd10.4,' km) -> Don t expect it ',
     &   'to move toward the center')
  334 format(1x,'Theoretical advance of flame front #',i4,': ',1pd11.4,
     &   ' km... greater than width of shell ',
     &   'under flame front (',1pd10.4,' km) -> Expect it to ',
     &   'move ',i4,' shell(s)')
  444 format(1x,'Flame front #',i2,' moved ',i4,' shells (expected ',
     &   i4,')')
  445 format(1x,'Flame front #',i2,' moved ',1pd11.4,' km (expected ',
     &   1pd11.4,' km)')
  446 format(1x,'Flame front #',i2,' moved ',1pd11.4,' msun (expected ',
     &   1pd11.4,' msun)')
  550 format(1x,'WARNING: Chemical composition at flame ',
     &  'front #',i4,' (C/O=',1pd10.4,') lower than first ',
     &  'tabulated composition (C/O=',1pd10.4,')')
  551 format(1x,'WARNING: Chemical composition at flame ',
     &  'front #',i4,' (C/O=',1pd10.4,') larger than last ',
     &  'tabulated composition (C/O=',1pd10.4,')')
  552 format(1x,'WARNING: Density at flame front #',i2,
     &   ' (rho=',1pd10.4,') lower than first ',
     &   'tabulated density (rho=',1pd10.4,')')
  553 format(1x,'WARNING: Density at flame front #',i2,
     &   ' (rho=',1pd10.4,') larger than last ',
     &   'tabulated density (rho=',1pd10.4,')')
  554 format(1x,'WARNING: Temperature at flame front #',i2,
     &   ' (T9=',1pd10.4,') lower than first ',
     &   'tabulated temperature (T9=',1pd10.4,')')
  666 format(1x,'Theoretical flame front #',i2,' speed: ',1pd9.2,
     &   ' cm/s')
  777 format(1x,'With current speed, expected time for the flame ',
     &   'front to reach the center: ',1pd10.4,' yr. With current ',
     &   'time step, this means in ',i8,' models.')
  881 format(
     &   '#  1             2                3      4     5          ',
     &   '6           7       8    9   10    11     12         ',
     &   '13        14         15        16         17         ',
     &   '18         19         20          21         22         ',
     &   '23         24         25         26          27          ',
     &   '28           29         30         31    32      ',
     &   '33         34          35         36           ',
     &   '37          38          39        40',/,
     &   '#model         time              dtn   nbcz   m_b1       ',
     &   'm_t1     m_t1-m_b1  ',
     &   'il1  ib1  it1 it1-ib1  rff1       mff1      t8ff1     ',
     &   'ro6ff1    xc12ff1    xo16ff1   xc12+xo16   ',
     &   'dmw(pf1)   drw(pf1)  drmin(pf1) drmax(pf1) <dr(pf1)>   ',
     &   'dr_bff1    dr_theo1  u_theo(ff1) u_real(ff1) utheo/u_re     ',
     &   'lmin     lummaxcz   ditheo direal XC12cbcz Lnu_pla_cbcz  ',
     &   'LC_cbcz    Lgrav_cbcz  Lnu_pla_rz    LC_rz      Lgrav_rz   ',
     &   'taulum'/,
     &   '#',274('-'))
  888 format(i5,1x,1pd22.16,1x,1pd10.4,1x,i2,1x,
     &   1pd10.4,1x,1pd10.4,1x,
     &   1pd10.4,1x,
     &   i4,1x,
     &   i4,1x,i4,1x,i4,1x,
     &   1pd10.4,1x,1pd10.4,1x,1pd10.4,1x,1pd10.4,1x,
     &   1pd10.4,1x,1pd10.4,1x,1pd10.4,1x,
     &   1pd10.4,1x,1pd10.4,1x,
     &   1pd10.4,1x,1pd10.4,1x,1pd10.4,1x,1pd10.4,1x,1pd11.4,1x,
     &   1pd11.4,1x,1pd11.4,1x,1pd10.4,1x,1pd11.4,1x,
     &   1pd11.4,1x,
     &   i5,1x,i5,1x,1pd10.4,1x,1pd11.4,1x,
     &   1pd11.4,1x,1pd11.4,1x,1pd11.4,1x,
     &   1pd11.4,1x,1pd11.4,1x,1pd10.4,1x)
      end


*******************************************************************
      BLOCK DATA FLAMES
*******************************************************************

      implicit none

      integer i

      integer nrhoflame,nt9flame,nxCOflame
      parameter (nrhoflame=12, nt9flame=6, nxCOflame=3)

      double precision rhoflame,t9flame,xCOflame,flamespeed
      dimension rhoflame(nrhoflame),t9flame(nt9flame),
     &   xCOflame(nxCOflame),flamespeed(nxCOflame,nrhoflame,nt9flame)
      common /twtflames/ rhoflame,t9flame,xCOflame,flamespeed

      data (xCOflame(i),i=1,nxCOflame) /
     &   0.25d0, 0.42857143d0, 0.66666667d0/

      data (rhoflame(i),i=1,nrhoflame) /
     &   1.0d6, 1.2d6, 1.5d6, 2.0d6, 3.0d6, 4.0d6, 5.0d6, 6.0d6,
     &   7.0d6, 8.0d6, 9.0d6, 1.0d7/

      data (t9flame(i),i=1,nt9flame) /
     &   0.6d0, 0.7d0, 0.8d0, 0.9d0, 1.0d0, 1.1d0/

***  X(C12)=0.2 X(O16)=0.8, decreasing density

c First line corrected by Gwen because second value probably false,
c 1.38e-4 -> 1.38e-3 makes the log(v(rho,T9)) plot look smoother:
      data (flamespeed(1,12,i),i=1,nt9flame) /
c     &   1.65d-4, 1.38d-4, 7.09d-3, 3.13d-2, 1.23d-1, 3.96d-1/
     &   1.65d-4, 1.38d-3, 7.09d-3, 3.13d-2, 1.23d-1, 3.96d-1/
      data (flamespeed(1,11,i),i=1,nt9flame) /
     &   1.60d-4, 1.32d-3, 6.94d-3, 3.02d-2, 1.20d-1, 3.82d-1/
      data (flamespeed(1,10,i),i=1,nt9flame) /
     &   1.54d-4, 1.25d-3, 6.72d-3, 2.95d-2, 1.15d-1, 3.68d-1/
      data (flamespeed(1,9,i),i=1,nt9flame) /
     &   1.22d-4, 1.19d-3, 6.46d-3, 2.93d-2, 1.11d-1, 3.52d-1/
      data (flamespeed(1,8,i),i=1,nt9flame) /
     &   1.15d-4, 1.13d-3, 6.21d-3, 2.86d-2, 1.07d-1, 3.31d-1/
      data (flamespeed(1,7,i),i=1,nt9flame) /
     &   1.00d-4, 1.06d-3, 5.99d-3, 2.75d-2, 1.00d-1, 3.09d-1/
      data (flamespeed(1,6,i),i=1,nt9flame) /
     &   9.89d-5, 9.59d-4, 5.66d-3, 2.56d-2, 9.40d-2, 2.83d-1/
      data (flamespeed(1,5,i),i=1,nt9flame) /
     &   9.00d-5, 9.12d-4, 5.40d-3, 2.39d-2, 8.48d-2, 2.48d-1/
      data (flamespeed(1,4,i),i=1,nt9flame) /
     &   9.00d-5, 8.04d-4, 4.82d-3, 2.10d-2, 7.46d-2, 2.05d-1/
      data (flamespeed(1,3,i),i=1,nt9flame) /
     &   8.55d-5, 7.59d-4, 4.65d-3, 2.00d-2, 6.60d-2, 1.63d-1/
      data (flamespeed(1,2,i),i=1,nt9flame) /
     &   7.89d-5, 6.64d-4, 4.40d-3, 1.85d-2, 5.87d-2, 1.03d-1/
c Last line corrected by Gwen because first and second value look strange
c on a log(v(rho,T9)) plot + pb in linear extrapolation (v<0)
c Replacement values come from the extrapolation of the 
c (rho=1e7,T9=0.5) values given in the grids for the 2 other 
c chemical compositions:
      data (flamespeed(1,1,i),i=1,nt9flame) /
c     &   1.00d-7, 7.60d-6, 4.18d-3, 1.73d-2, 5.08d-2, 6.08d-2/
     &   8.46d-5, 7.43d-4, 4.18d-3, 1.73d-2, 5.08d-2, 6.08d-2/

***  X(C12)=0.3 X(O16)=0.7, decreasing density

      data (flamespeed(2,12,i),i=1,nt9flame) /
     &   2.00d-4, 2.04d-3, 1.03d-2, 4.60d-2, 1.88d-1, 5.89d-1/
      data (flamespeed(2,11,i),i=1,nt9flame) /
     &   1.85d-4, 1.95d-3, 1.00d-2, 4.53d-2, 1.85d-1, 5.77d-1/
      data (flamespeed(2,10,i),i=1,nt9flame) /
     &   1.71d-4, 1.88d-3, 9.87d-3, 4.49d-2, 1.83d-1, 5.73d-1/
      data (flamespeed(2,9,i),i=1,nt9flame) /
     &   1.56d-4, 1.78d-3, 9.62d-3, 4.43d-2, 1.76d-1, 5.77d-1/
      data (flamespeed(2,8,i),i=1,nt9flame) /
     &   1.42d-4, 1.73d-3, 9.36d-3, 4.36d-2, 1.70d-1, 5.44d-1/
      data (flamespeed(2,7,i),i=1,nt9flame) /
     &   1.27d-4, 1.64d-3, 9.20d-3, 4.25d-2, 1.60d-1, 5.06d-1/
      data (flamespeed(2,6,i),i=1,nt9flame) /
     &   1.11d-4, 1.52d-3, 8.87d-3, 4.05d-2, 1.50d-1, 4.70d-1/
      data (flamespeed(2,5,i),i=1,nt9flame) /
     &   9.55d-5, 1.41d-3, 8.25d-3, 3.77d-2, 1.39d-1, 4.30d-1/
      data (flamespeed(2,4,i),i=1,nt9flame) /
     &   9.35d-5, 1.26d-3, 7.22d-3, 3.29d-2, 1.25d-1, 3.83d-1/
      data (flamespeed(2,3,i),i=1,nt9flame) /
     &   9.32d-5, 1.23d-3, 7.20d-3, 3.23d-2, 1.16d-1, 3.42d-1/
      data (flamespeed(2,2,i),i=1,nt9flame) /
     &   9.22d-5, 1.14d-3, 6.89d-3, 3.03d-2, 1.06d-1, 3.11d-1/
      data (flamespeed(2,1,i),i=1,nt9flame) /
     &   9.15d-5, 1.09d-3, 6.78d-3, 2.95d-2, 1.02d-1, 2.77d-1/

***  X(C12)=0.4 X(O16)=0.6, decreasing density

c      data (flamespeed(3,12,i),i=1,nt9flame) /
c     &   2.00d-4, 2.04d-3, 1.03d-2, 4.60d-2, 1.88d-1, 5.89d-1/ 
      data (flamespeed(3,12,i),i=1,nt9flame) /
     &   2.25d-4, 2.88d-3, 1.50d-2, 6.71d-2, 2.61d-1, 8.63d-1/
      data (flamespeed(3,11,i),i=1,nt9flame) /
     &   2.00d-4, 2.75d-3, 1.45d-2, 6.46d-2, 2.54d-1, 8.40d-1/
      data (flamespeed(3,10,i),i=1,nt9flame) /
     &   1.91d-4, 2.61d-3, 1.40d-2, 6.30d-2, 2.46d-1, 8.10d-1/
      data (flamespeed(3,9,i),i=1,nt9flame) /
     &   1.75d-4, 2.49d-3, 1.34d-2, 6.18d-2, 2.39d-1, 7.77d-1/
      data (flamespeed(3,8,i),i=1,nt9flame) /
     &   1.65d-4, 2.35d-3, 1.29d-2, 6.01d-2, 2.29d-1, 7.36d-1/
      data (flamespeed(3,7,i),i=1,nt9flame) /
     &   1.40d-4, 2.23d-3, 1.25d-2, 5.81d-2, 2.18d-1, 6.96d-1/
      data (flamespeed(3,6,i),i=1,nt9flame) /
     &   1.50d-4, 2.02d-3, 1.18d-2, 5.49d-2, 2.06d-1, 6.48d-1/
      data (flamespeed(3,5,i),i=1,nt9flame) /
     &   1.25d-4, 1.94d-3, 1.14d-2, 5.15d-2, 1.89d-1, 5.95d-1/
      data (flamespeed(3,4,i),i=1,nt9flame) /
     &   1.25d-4, 1.76d-3, 1.03d-2, 4.66d-2, 1.74d-1, 5.42d-1/
      data (flamespeed(3,3,i),i=1,nt9flame) /
     &   1.10d-4, 1.73d-3, 1.01d-2, 4.51d-2, 1.64d-1, 5.00d-1/
      data (flamespeed(3,2,i),i=1,nt9flame) /
     &   1.00d-4, 1.65d-3, 9.69d-3, 4.30d-2, 1.56d-1, 4.68d-1/
      data (flamespeed(3,1,i),i=1,nt9flame) /
     &   9.90d-5, 1.60d-3, 9.42d-3, 4.18d-2, 1.47d-1, 4.33d-1/

      end

*******************************************************************
      SUBROUTINE mylocalcompo (nc,lprint)
*******************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.eos'
      include 'evolcom.nuc'
      include 'evolcom.spec'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer i,j,nc,permute,kcxsp,ktmp
      double precision cxsp,tmp,cxsppart
      logical lprint

      dimension cxsp(nsp),kcxsp(nsp)

      if (lprint) write (nout,444) nc,m(nc)/msun,exp(lnt(nc)),ro(nc),
     &     eta(nc)

      do i = 1,nsp
         cxsp(i) = xsp(nc,i)
         kcxsp(i) = i
      enddo

c     cxsp(i) est l'element classé i pour l'abondance
c     kcxsp(i) est le numero de l'element classé i pour l'abondance

      do j = 1,nsp
         permute = 0
         do i = 1,nsp-1
            if (cxsp(i+1).gt.cxsp(i)) then
               tmp = cxsp(i+1)
               cxsp(i+1) = cxsp(i)
               cxsp(i) = tmp
               ktmp = kcxsp(i+1)
               kcxsp(i+1) = kcxsp(i)
               kcxsp(i) = ktmp
               permute = 1
            endif
         enddo
      enddo

      j = 0
      cxsppart = 0.d0
      do i = 1,nsp
         cxsppart = cxsppart + cxsp(i)
         if (cxsppart.gt.0.99d0) then
            j = i
            if (lprint) write (nout,888) nc,j,cxsppart*1.d2
            goto 30
         endif
      enddo

 30   if (j.eq.0.and.lprint) write (nout,*) 'Did not find the number ',
     &     'of most abundant nuclei'
      if (lprint) write (nout,22) (elem(kcxsp(i)),cxsp(i),i=1,j)
 22   format (4(2x,a5,' (',1pe12.6,') '))

 444  format (3x,'At shell #',i4,': m/msun = ',1pd10.4,', T = ',1pd8.2,
     &     ', rho = ',1pd8.2,', eta = ',1pd9.2)

 888  format (3x,'At shell #',i4,': ',i3,
     &      ' element(s) = ',f6.2,'% of the shell abundance:')
      end


*********************************************************************
      DOUBLE PRECISION FUNCTION interp3dlog (x,y,z,x1,x2,y1,y2,z1,z2,
     &      fx1y1z1,fx1y2z1,fx1y1z2,fx1y2z2,
     &      fx2y1z1,fx2y2z1,fx2y1z2,fx2y2z2)
*********************************************************************

      implicit none

      include 'evolpar.star'

      double precision x,y,z,x1,x2,y1,y2,z1,z2,
     &      px1,px2,py1,py2,pz1,pz2,
     &      fx1y1z1,fx1y2z1,fx1y1z2,fx1y2z2,
     &      fx2y1z1,fx2y2z1,fx2y1z2,fx2y2z2
      

      if (x1.eq.x2) then
         write (nout,*) 'interp3dlog ERROR: x1 = x2 !'
         interp3dlog = 0.d0
         return
      endif

      if (y1.eq.y2) then
         write (nout,*) 'interp3dlog ERROR: y1 = y2 !'
         interp3dlog = 0.d0
         return
      endif

      if (z1.eq.z2) then
         write (nout,*) 'interp3dlog ERROR: z1 = z2 !'
         interp3dlog = 0.d0
         return
      endif

      if (fx1y1z1.le.0) then
         write (nout,*) 'interp3dlog ERROR: fx1y1z1 = ',fx1y1z1,' <= 0'
         interp3dlog = 0.d0
         return
      endif

      if (fx1y1z2.le.0) then
         write (nout,*) 'interp3dlog ERROR: fx1y1z2 = ',fx1y1z2,' <= 0'
         interp3dlog = 0.d0
         return
      endif

      if (fx1y2z1.le.0) then
         write (nout,*) 'interp3dlog ERROR: fx1y2z1 = ',fx1y2z1,' <= 0'
         interp3dlog = 0.d0
         return
      endif

      if (fx1y2z2.le.0) then
         write (nout,*) 'interp3dlog ERROR: fx1y2z2 = ',fx1y2z2,' <= 0'
         interp3dlog = 0.d0
         return
      endif

      if (fx2y1z1.le.0) then
         write (nout,*) 'interp3dlog ERROR: fx2y1z1 = ',fx2y1z1,' <= 0'
         interp3dlog = 0.d0
         return
      endif

      if (fx2y1z2.le.0) then
         write (nout,*) 'interp3dlog ERROR: fx2y1z2 = ',fx2y1z2,' <= 0'
         interp3dlog = 0.d0
         return
      endif

      if (fx2y2z1.le.0) then
         write (nout,*) 'interp3dlog ERROR: fx2y2z1 = ',fx2y2z1,' <= 0'
         interp3dlog = 0.d0
         return
      endif

      if (fx1y1z1.le.0) then
         write (nout,*) 'interp3dlog ERROR: fx1y1z1 = ',fx1y1z1,' <= 0'
         interp3dlog = 0.d0
         return
      endif

      if (fx2y2z2.le.0) then
         write (nout,*) 'interp3dlog ERROR: fx2y2z2 = ',fx2y2z2,' <= 0'
         interp3dlog = 0.d0
         return
      endif

      px2 = (x-x1)/(x2-x1)
      py2 = (y-y1)/(y2-y1)
      pz2 = (z-z1)/(z2-z1)

      px1 = 1.d0-px2
      py1 = 1.d0-py2
      pz1 = 1.d0-pz2

      interp3dlog = exp(
     &              px1*py1*pz1*log(fx1y1z1) +
     &              px1*py1*pz2*log(fx1y1z2) +
     &              px1*py2*pz1*log(fx1y2z1) +
     &              px1*py2*pz2*log(fx1y2z2) +
     &              px2*py1*pz1*log(fx2y1z1) +
     &              px2*py1*pz2*log(fx2y1z2) +
     &              px2*py2*pz1*log(fx2y2z1) +
     &              px2*py2*pz2*log(fx2y2z2)
     &              )

      end

*******************************************************************
      DOUBLE PRECISION FUNCTION dtnflame (currentdt)
*******************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'  
      include 'evolcom.therm'
      include 'evolcom.var'
      include 'evolcom.flame'  
      include 'evolcom.mod'

      integer i
      double precision facdtnflame,currentdt,typicalspeed,
     &      dtnvt,dtnvr,dtnvm

      integer iterp
      double precision vdtn
      common /viter/ vdtn,iterp

      dtnflame = 1.d99

* Choose:
*     ftacc = 1/6 at the begining
*     ftacc = 1/4 if convection first disappears
*     ftacc = 0 for rff <~ 50 km
      facdtnflame = 1.d0/ftacc
      typicalspeed = 1.d-2

      if (ftacc.gt.1.d30) then
         print *,'Time step not constrained by flame speed'
         return
      endif

      if (model.eq.modeli.and.nbcz.gt.0) then
         dtnflame = vdtn
c         dtnflame = vdtn*1.3d0
         write (nout,300) dtnflame*seci
         return
      endif

      if (lpreflame) then
         dtnflame = (drpreflame/typicalspeed)/(facdtnflame*height)
         print *,'dtn_preflame = ',dtnflame*seci,' yr'
         if (currentdt.gt.dtnflame) then
            write (nout,223) currentdt*seci,dtnflame*seci
         endif
      endif

      do i=1,nbcz

         if (iminlumprof(i).gt.0.and.ibCOcz(i).gt.iminlumprof(i)) then

            if (thurff(i).gt.0.d0) then
               dtnvt = (r(ibCOcz(i))-r(iminlumprof(i)))/
     &            (facdtnflame*thurff(i))
               if (dtnflame.gt.dtnvt) then
                  dtnflame = dtnvt
                  write (nout,110) i,dtnflame*seci
               endif
            endif 

            if (urff(i).lt.0.d0) then
               dtnvr = (r(ibCOcz(i))-r(iminlumprof(i)))/
     &            (facdtnflame*abs(urff(i)))
               if (dtnflame.gt.dtnvr) then
                  dtnflame = dtnvr
                  write (nout,111) i,dtnflame*seci
               endif
            endif

            if (umff(i).lt.0.d0) then
               dtnvm = (m(ibCOcz(i))-m(iminlumprof(i)))/
     &            (facdtnflame*abs(umff(i)))
               if (dtnflame.gt.dtnvm) then
                  dtnflame = dtnvm
                  write (nout,112) i,dtnflame*seci
               endif
            endif

            if (currentdt.gt.dtnflame) then
               write (nout,222) i,currentdt*seci,dtnflame*seci
            endif

         endif
      enddo

  110 format('Flame front #',i2,' dtnflame = ',1pd10.4,' yr ',
     &   '[ based on vrtheo ]')
  111 format('Flame front #',i2,' dtnflame = ',1pd10.4,' yr ',
     &   '[ based on vr ]')
  112 format('Flame front #',i2,' dtnflame = ',1pd10.4,' yr ',
     &   '[ based on vm ]')
  222 format('WARNING: Time step must be constrained by front flame #',
     &   i2,': dtn = ',1pd10.4,' yr -> ',1pd10.4,' yr')
  223 format('WARNING: Time step must be constrained by pre-flame',
     &   ': dtn = ',1pd10.4,' yr -> ',1pd10.4,' yr')
  300 format(' model == modeli && nbcz > 0 -> dtnflame = vdtn = ',
     &   1pd10.4,' yr')
      end
