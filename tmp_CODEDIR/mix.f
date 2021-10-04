

***********************************************************************

      SUBROUTINE mix (dm,imin,imax,kl,imix,partialmix)

************************************************************************
* Instantaneous mixing of species in convective zones                  *
* imix = 0 : default, do the mixing
* imix = 1 : special treatment of deuterium (accretion)
* imix = 2 : before nucleosynthesis, do the mixing
* imix = 3 : after nucleosynthesis, do the mixing
* imix = 5 : compute turnover timescales only
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.nuc'
      include 'evolcom.spec'
      include 'evolcom.teq'

c..   input
      double precision dm(nsh)

      integer imin,imax,imix
      integer ielabond
      integer i,ii,j,k,l,kl
      integer klenv,klpulse,klcore

      logical partialmix,irecede

      double precision mconv,drconv,tmix,dtni,tconvmix,dmeff
      double precision x,xk,xj
      double precision x0,x1
      double precision t,vt
      double precision r,vr

      dimension tmix(nsh)
      dimension ielabond(7)

      common /vartvt/ t(nsh),vt(nsh)
      common /varrvr/ r(nsh),vr(nsh)
      common /overshoot/ klcore,klenv,klpulse

      if (imin.ge.imax.or.dtn.lt.1.d-12) return

      tconvmix = 0.d0
      partialmix = .false.

*____________________________________________
***   compute convective turnover time scales
*--------------------------------------------

      if (imix.ne.3) then
         do j = imin,imax-1
            tconvmix = tconvmix+(r(j+1)-r(j))/sconv(j)
         enddo
         turnover(kl) = tconvmix
         if (nsconv.eq.1) then
            turnenv = turnover(kl)
         else
            if (t(imin).gt.1.d6.and.r(imax)/r(nmod).gt.0.9d0) then
               turnenv = turnover(kl)
            endif
         endif
      else
         tconvmix = turnover(kl)
      endif
      if (nuclopt.eq.'p'.and.tconvmix.gt.dtn) partialmix = .true.

      if (imix.eq.5) return

c..   do not mix before nucleosynthesis (imix=2) if dtnmix = t
c..   or tconvmix >> dtn (partialmix = .true.)
      if (imix.eq.2.and.(dtnmix.eq.'t'.or.partialmix)) return


*__________________________________
***   Instantaneous mixing
***   if dtnmix = f or tconv >> dtn
*----------------------------------

      if (dtnmix.eq.'f'.or..not.partialmix) then
         mconv = 0.d0
         do j = imin,imax
            mconv = mconv+dm(j)
         enddo
!$OMP PARALLEL REDUCTION(+:x)
!$OMP DO
         do l = 2,nis
            x = 0.d0
            do j = imin,imax
               x = x+vxsp(j,l)*dm(j)
            enddo
            xj = x/mconv
            do j = imin,imax
               xsp(j,l) = xj
               ysp(j,l) = xj/anuc(l)
            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL

*___________________________________________________
***   special treatment of LiBeB during HBB
***   restore old (equilibrium) abundances for LiBeB
*---------------------------------------------------

         if (imix.eq.2.or.imix.eq.3) write (nout,200) imin,imax
         if (hbb.and.kl.eq.klenv.and.partialmix) then
            do k = imin,imax
               xsp(k,ili7) = vxsp(k,ili7)
               xsp(k,ili6) = vxsp(k,ili6)
               xsp(k,ib11) = vxsp(k,ib11)
               ysp(k,ili7) = xsp(k,ili7)/anuc(ili7)
               ysp(k,ili6) = xsp(k,ili6)/anuc(ili6)
               ysp(k,ib11) = xsp(k,ib11)/anuc(ib11)
            enddo
            if (imix.eq.2.or.imix.eq.3) write (nout,250) imin,imax
         endif

      else

*_______________________________________________
***   time-dependent mixing when tauconv < dtn
***   if dtnmix = t and partialmix = .true
***
*** Sparks Endal, 1980, ApJ 237, 130
*** Chieffi et al. 1998, ApJ 502,737
*-----------------------------------------------

c..   check that the convective zone has moved into regions
c..   of different chemical composition
         if (iter.eq.1) then
            irecede = .false.
         elseif (mixopt) then
            irecede = .true.
            ielabond(1) = ih1
            ielabond(2) = ili7
            ielabond(3) = ic12
            ielabond(4) = io16
            ielabond(5) = ine20
            ielabond(6) = img24
            ielabond(7) = isi28
            do i = 1,7
               ii = ielabond(i)
               if (log10(abs(xsp(imin,ii)/vxsp(imin,ii))).gt.1.d-2.or.
     &              log10(abs(xsp(imax,ii)/vxsp(imax,ii))).gt.1.d-2) 
     &              then
                  irecede = .false.
                  goto 10
               endif
            enddo
         else
            irecede = .false.
         endif
 10      if (irecede.and.imix.ne.3) return

         write (nout,100) imin,imax,tconvmix*seci,dtn*seci
         dtni = 1.d0/dtn
         do k = imin,imax
            tmix(k) = 0.d0
            do j = k-1,imin,-1
               drconv = r(j+1)-r(j)
               tmix(j) = tmix(j+1)+drconv/sconv(j)
            enddo
            do j = k+1,imax
               drconv = r(j)-r(j-1)
               tmix(j) = tmix(j-1)+drconv/sconv(j-1)
            enddo
            do i = 2,nis
               xk = 0.d0
               mconv = 0.d0
               do l = imin,imax
                  dmeff = dm(l)*exp(-tmix(l)*dtni)
                  xk = xk+(vxsp(l,i)-vxsp(k,i))*dmeff
                  mconv = mconv+dmeff
               enddo
               xsp(k,i) = vxsp(k,i)+xk/mconv
            enddo
         enddo
         return
      endif

*_______________________________________________________________
***   specific mixing of light elements, in case of call ydequil
*---------------------------------------------------------------

      if (imix.eq.1) then
         imin = max(imin,idcur)
         x0 = 0.d0
         x1 = 0.d0
         mconv = 0.d0
         do j = imin,imax
            x0 = x0+vyd0(j)*dm(j)
            x1 = x1+yd0(j)*dm(j)
            mconv = mconv+dm(j)
         enddo
         x0 = x0/mconv
         x1 = x1/mconv
         do j = imin,imax
            vyd1(j) = x0
            yd1(j) = x1
         enddo
         do l = 2,4
            x = 0.d0
            do j = imin,imax
               x = x+vxsp(j,l)*dm(j)
            enddo
            x = x/mconv
            do j = imin,imax
               vxsp(j,l) = x
            enddo
         enddo
      endif

 100  format (' partial mixing in shells : [',i4,',',i4,'] , ',
     &     'tau_conv = ',1pe9.3,' > dtn = ',1pe9.3)
 200  format ('  instantaneous mixing in convective zone [',i4,',',
     &     i4,']')
 250  format ('   HBB : special mixing for Li6, Li7 and B11 in the ',
     &     'envelope [',i4,',',i4,']')

      return
      end
