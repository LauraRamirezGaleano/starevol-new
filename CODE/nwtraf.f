

************************************************************************

      SUBROUTINE nwtraf (*,error)

************************************************************************
* Apply the Newton-Raphson method for the five equations               *
* of the stellar structure evolution                                   *
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
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.eos'
      include 'evolcom.ion'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'

      include 'evolcom.grad'

      integer iiter
      integer ncm(neq)
      integer neq2
      integer iitmax,i1,imerk,icrasht
      integer klenv,klpulse,klcore
      integer i,j,k,kl,l,mm
      integer error
      integer itmax,ishockt,ishockb

      parameter (neq2 = 2*neq)

      double precision xmod,xa,eps
      double precision t,vt,r,vr
      double precision dx(neq),dv(neq),w(nsh,neq),a(neq,neq2),
     &     um(neq,neq),v(nsh,neq,neq),corm(neq),alpha1(neq)
      double precision xishk,xx
      double precision Lmax,tmax,mackmax,enucmax
      double precision dp,ag

      logical corrmax

      common /calcul/ iiter
      common /convergence/ imerk,icrasht
      common /overshoot/ klcore,klenv,klpulse
      common /varia/ xmod(nsh,neq),xa(nsh,neq),eps(neq)
      common /vartvt/ t(nsh),vt(nsh)
      common /varrvr/ r(nsh),vr(nsh)

      common /hydrodyn/ Lmax,tmax,mackmax,enucmax,itmax,ishockb,ishockt
      
      iitmax = itermax

      corrmax = .true.
      imerk = -1
      iter = 0
      error = 0

 10   iter = iter+1
      iiter = iter
      if (iter.le.1) then
         do i = 1,neq
            alpha1(i) = alpha*0.5d0
         enddo
      endif
      if (iter.eq.2) then
         do i = 1,neq
            alpha1(i) = alpha
         enddo
      endif
c      print *, 'dtn nwtraf',dtn*seci
*________________________________________________
***   calculation of all the constitutive physics
*------------------------------------------------
      
      call thermo (error)
      if (error.gt.0) return 1
      
*___________________________________________________
***   application of the central boundary conditions
*---------------------------------------------------

      ish = 1

      call centeq

      a(1:neq,1:neq) = eq(1:neq,1:neq)
      a(1:neq,neq+1:neq2) = 0.d0
      do j = 1,neq
c         do mm = 1,neq
c            a(mm,j) = eq(mm,j)
c            a(mm,neq+j) = 0.d0
c         enddo
         a(j,neq+j) = 1.d0
      enddo
      call nwrmat (a,um,neq,neq,error)
      if (error.gt.0) then
         write (nout,'("error nwtraf, central shell 1",i4)') crz(ish)
         do l = 1,neq
            do i = 1,neq
               write (nout,'(i2,i2,3(1x,1pe10.3))') i,l,eq(i,l),
     &              eq(i,neq+l),eq(i,2*neq+l)
            enddo
         enddo
         return 1
      endif
      do mm = 1,neq
         w(ish,mm) = 0.d0
         do k = 1,neq
            w(ish,mm) = w(ish,mm)+um(mm,k)*eq(k,neq1)
            v(ish,mm,k) = 0.d0
            do j = 1,neq
               v(ish,mm,k) = v(ish,mm,k)+um(mm,j)*eq(j,neq+k)
            enddo
         enddo
      enddo
      
*__________________________________________________
***   determination of the internal shell structure
*--------------------------------------------------

      do ish = 2,nmod1
         i1 = ish-1

         call inteq

         do j = 1,neq
            do l = 1,neq
               a(j,l) = eq(j,neq+l)
               do k = 1,neq
                  a(j,l) = a(j,l)-eq(j,k)*v(i1,k,l)
               enddo
               a(j,neq+l) = 0.d0
            enddo
            a(j,neq+j) = 1.d0
         enddo
         call nwrmat (a,um,neq,neq,error)
         if (error.gt.0) then
            write (nout,'("error nwtraf, shell #",i4,i3)') ish,crz(ish)
            do l = 1,neq
               do i = 1,neq
                  write (nout,'(i2,i2,3(1x,1pe10.3))') i,l,eq(i,l),
     &                 eq(i,neq+l),eq(i,2*neq+l)
               enddo
            enddo
            return 1
         endif
         do l = 1,neq
            dv(l) = -eq(l,neq1)
            do k = 1,neq
               dv(l) = dv(l)+eq(l,k)*w(i1,k)
            enddo
         enddo
         do k = 1,neq
            w(ish,k) = 0.d0
            do l = 1,neq
               w(ish,k) = w(ish,k)-um(k,l)*dv(l)
               v(ish,k,l) = 0.d0
               do j = 1,neq
                  v(ish,k,l) = v(ish,k,l)+um(k,j)*eq(j,neq2+l)
               enddo
            enddo
         enddo
      enddo

*___________________________________________________
***   application of the surface boundary conditions
*---------------------------------------------------

      ish = nmod
      i1 = ish-1
      
      call surfeq
      
      do j = 1,neq
         do l = 1,neq
            a(j,l) = eq(j,neq+l)
            do k = 1,neq
               a(j,l) = a(j,l)-eq(j,k)*v(i1,k,l)
            enddo
            a(j,neq+l) = 0.d0
         enddo
         a(j,neq+j) = 1.d0
      enddo
      call nwrmat (a,um,neq,neq,error)
      if (error.gt.0) then
         write (nout,'("error nwtraf, surface #",i4,a1)') ish,crz(ish)
         do l = 1,neq
            do i = 1,neq
               write (nout,'(i2,i2,3(1x,1pe10.3))') i,l,eq(i,l),
     &              eq(i,neq+l),eq(i,2*neq+l)
            enddo
         enddo
         return 1
      endif
      do j = 1,neq
         dv(j) = -eq(j,neq1)
         do k = 1,neq
            dv(j) = dv(j)+eq(j,k)*w(i1,k)
         enddo
      enddo
*_______________________________________________________________
***   determination of the evolved structure (iterative process)
*---------------------------------------------------------------

      do mm = 1,neq
         dx(mm) = 0.d0
         do j = 1,neq
            dx(mm) = dx(mm)+um(mm,j)*dv(j)
         enddo
      enddo
      do i = 1,neq
         corm(i) = 0.d0

         if (xmod(ish,i).ne.0.d0.and.((numeric.eq.1.and.t(ish).gt
     &        .tmaxioHe).or.numeric.ne.1)) corm(i) = dx(i)/xmod(ish,i)
c... norbert
c         if (i.ne.5.and.xmod(ish,i).ne.0.d0) corm(i) = dx(i)/xmod(ish,i)
c         if (i.eq.5) corm(i) = dx(i)
c         epsmax(5) = 1.d50
         ncm(i) = ish

cccc---------------------------------------------
ccc TEST : small corrections in the upper layers
c          if (t(ish).lt.1.d5) then
c             do l = 1,neq
c                alpha1(l) = 0.1d0
c             enddo
c          endif
cccc---------------------------------------------
         xmod(ish,i) = xmod(ish,i)+alpha1(i)*dx(i)
      enddo
      do ish = nmod1,1,-1
cccc---------------------------------------------
ccc TEST : small corrections in the upper layers
C          if (t(ish).lt.1.d5) then
C             do l = 1,neq
C                alpha1(l) = 0.1d0
C             enddo
C          else
C             do l = 1,neq
C                alpha1(l) = alpha2(l)
C             enddo
C          endif
cccc---------------------------------------------
         do k = 1,neq
            dv(k) = -w(ish,k)
            do l = 1,neq
               dv(k) = dv(k)-v(ish,k,l)*dx(l)
            enddo
         enddo
         
         do k = 1,neq
            dx(k) = dv(k)
            xishk = xmod(ish,k)
            if (ish.eq.1.or.xishk.eq.0.d0) goto 20
            xx = dx(k)/xmod(ish,k)
            if (k.eq.5.and.abs(xmod(ish,5)).lt.1.d-15*Lmax) xx = 0.d0
            if (abs(xx).le.abs(corm(k)).or.(numeric.eq.1.and.t(ish).lt.
     &           tmaxioHe)) goto 20
            corm(k) = xx
            ncm(k) = ish
 20         xmod(ish,k) = xmod(ish,k)+alpha1(k)*dx(k)
         enddo
         
      enddo
      
      xmod(1,2) = max(xmod(1,2),log(1.d-99))
      t(1) = exp(min(xmod(1,4),23.d0))
      r(1) = exp(xmod(1,2))
      do ish = 2,nmod
         t(ish) = exp(min(xmod(ish,4),23.d0))
         r(ish) = exp(min(xmod(ish,2),50.d0))
      enddo
      if (icorr.eq.'t') then
         do ish = neff,2,-1
c..  check if pressure gradient is a decreasing function of mass 
c..  except in the outer layers
            if (abs(corm(2)).lt.eps(2).or.(numeric.eq.1.and.t(ish).gt.
     &           tmaxioHe)) then
               dp = p(ish)-p(ish-1)
c            accel(ish) =(xmod(ish,1)-xa(ish,1)-psi(ish)*(xmod(ish,1)-
c     &           xmod(ish-1,1)))/dtn
               ag = accel(ish)+g*m(ish)/(r(ish)*r(ish))
               if (dp*ag.ge.0.d0) then
                  corm(2) = 9.99d0
                  ncm(2) = ish
                  exit
               endif
            endif
         enddo
      endif
      
      if (iter.eq.1) write (nout,1000)
      if (nsconv.le.1) write (nout,1001) iter,(alpha1(i),ncm(i),corm(i),
     &     i = 1,neq),(novlim(1,i),i = 3,4),nsconv,scr
      if (nsconv.eq.2) write (nout,1002) iter,(alpha1(i),ncm(i),corm(i),
     &     i = 1,neq),((novlim(kl,i),i = 3,4),kl = 1,2),max(nsconv,1),
     &     scr
      if (nsconv.eq.3) write (nout,1003) iter,(alpha1(i),ncm(i),corm(i),
     &     i = 1,neq),((novlim(kl,i),i = 3,4),kl = 1,3),max(nsconv,1),
     &     scr
      if (nsconv.ge.4) write (nout,1004) iter,(alpha1(i),ncm(i),corm(i),
     &     i = 1,neq),((novlim(kl,i),i = 3,4),kl = 1,4),max(nsconv,1),
     &     scr
      if (no.eq.1) then
         if (iter.eq.1) write (90,1000)
         if (nsconv.le.1) write (90,1001) iter,(alpha1(i),ncm(i),
     &        corm(i),i = 1,neq),(novlim(1,i),i = 3,4),nsconv,scr
         if (nsconv.eq.2) write (90,1002) iter,(alpha1(i),ncm(i),
     &        corm(i),i = 1,neq),((novlim(kl,i),i = 3,4),kl = 1,2),
     &        max(nsconv,1),scr
         if (nsconv.eq.3) write (90,1003) iter,(alpha1(i),ncm(i),
     &        corm(i),i = 1,neq),((novlim(kl,i),i = 3,4),kl = 1,3),
     &        max(nsconv,1),scr
         if (nsconv.ge.4) write (90,1004) iter,(alpha1(i),ncm(i),
     &        corm(i),i = 1,neq),((novlim(kl,i),i = 3,4),kl = 1,4),
     &        max(nsconv,1),scr
      endif
      
c.... NEW
      do i = 1,5
         corm(i) = abs(corm(i))
      enddo

      if (corm(1).gt.epsmax(1).or.corm(2).gt.epsmax(2).or.corm(3).gt.
     &     epsmax(3).or.corm(4).gt.epsmax(4).or.corm(5).gt.epsmax(5))
     &     then
         error = 6
         return 1
      endif
***   determine correction threshold in case nucreg = 3
      corrmax = abs(corm(2)).gt.eps(2)*1.d2.or.abs(corm(3)).gt.
     &     eps(3)*10.d0.or.abs(corm(4)).gt.eps(4)*10.d0.or.
     &     abs(corm(5)).gt.eps(5)*1.d5.or.abs(corm(1)).gt.eps(1)
*___________________________________________________________________
***   calculation of the associated nucleosynthesis and diffusion,
***   if required at each Newton-Raphson iteration for the structure
***   case : nucreg = 3
*-------------------------------------------------------------------

      if (nucreg.eq.3.and.iter.ge.(itermin-1).and..not.corrmax) then
         klcore = 0
         klpulse = 0
         klenv = 0
         if ((microdiffus.or.rotation.or.diffusconv.or.difftacho.or.
     &        thermohaline).and.nmixd.eq.0) then
c..   in case nmixd > 0, diffusion treated independently in netdiff
            if (idiffnuc.eq.2) then
               if (nucreg.ne.0) call chedif (error)
               if (error.gt.0) return 1
               call diffusion (1,error)
               if (error.gt.0) return 1
            else
               call diffusion (1,error)
               if (error.gt.0) return 1
               if (nucreg.ne.0) call chedif (error)
               if (error.gt.0) return 1
            endif
            if (idiffnuc.eq.3) then
               write (nout,*) 'Processing diffusive mixing'
               if (no.eq.1) write (90,*) 'Processing diffusive mixing'
               call diffusion (1,error)
               if (error.gt.0) return 1
            endif
         else
            if (nmixd.gt.0.and.nsconv.gt.0) then
               if (novlim(1,3).eq.1) klcore = 1
               do kl = 1,nsconv
                  if ((tau(novlim(kl,4)).lt.1.d2.or.t(novlim(kl,4))
     &                 .lt.1.d6).and.mr(novlim(kl,4)).gt.0.7d0.and.
     &                 klenv.eq.0) klenv = kl
                  if (t(novlim(kl,3)).gt.2.d8.and.novlim(kl,3).gt.1
     &                 .and.nphase.eq.5.and.klpulse.eq.0) klpulse = kl
               enddo
               if (nsconv.eq.1.and.klcore.eq.1) then
                  klenv = 0
                  klpulse = 0
               endif
               if (nmixd.eq.1) klenv = 0
               if (nmixd.eq.2) klcore = 0
               if (nmixd.le.3) klpulse = 0
            endif
            call chedif (error)
            if (error.gt.0) return 1
         endif
      endif
      if (iter.gt.iitmax) then
         error = 99
         return 1
      endif
      if (iacc) call nwracc (corm,iter,alpha1)
c..   if lmix (convective boundaries have moved and mixing has been performed) 
c..   do another iteration
      if (iter.lt.itermin.or.lmix) goto 10

      
      imerk = 0
      do i = 1,neq
         if (abs(corm(i)).gt.eps(i)) imerk = imerk+1
      enddo
      if (imerk.gt.0) goto 10

 1000 format (/,1x,'#','|',3x,'Variable u(m)',4x,'|',3x,'Variable r(m)',
     &     4x,'|',2x,'Variable lnf(m)',3x,'|',2x,'Variable lnT(m)',3x,
     &     '|',3x,'Variable l(m)',4x,'|',1x,'first conv. zones',/)
 1001 format (i2,5('|',0pf5.3,1x,i4,1x,1pe9.2),'|',(1x,i4,'-',i4),1x,
     &     i2,a1)
 1002 format (i2,5('|',0pf5.3,1x,i4,1x,1pe9.2),'|',2(1x,i4,'-',i4),1x,
     &     i2,a1)
 1003 format (i2,5('|',0pf5.3,1x,i4,1x,1pe9.2),'|',3(1x,i4,'-',i4),1x,
     &     i2,a1)
 1004 format (i2,5('|',0pf5.3,1x,i4,1x,1pe9.2),'|',4(1x,i4,'-',i4),1x,
     &     i2,a1)

      return
      end
