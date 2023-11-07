

************************************************************************

      SUBROUTINE rinimod (uread,ureadrot,uwrite,uwriterot,flag)

************************************************************************
* Read the initial stellar model                                       *
* 09/03: mass fraction of photons & CAPTN not stored in binary anymore *
* 23/09: arguments added to the subroutine name to manage reading and  *
*        writing of binary files, including crash cases                *
*  flag  = 0 : first reading                                           *
*  flag <> 0 : failed computation, reload previous model               *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.ion'
      include 'evolcom.lum'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.opa'
      include 'evolcom.rot'
      include 'evolcom.conv'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.transp'
      include 'evolcom.var'

      integer nspr
      integer error,last_model
      integer i,j,k,l,kl,jl,it
      integer uread,ureadrot,uwrite,uwriterot,flag
      integer ishockb,ishockt,itmax
      integer klenv,klpulse,klcore

      character*1 crzc(nsh)
      logical tbot,ttop,ierr,test

      double precision vom,vrray,vxpsi
      double precision dift,dife,Dhold
      double precision xc11,xsi31,xsi32,xp32,xp33,xar36,FeH
      double precision abmurj,abmuj
      double precision zstart,y,xspr,vxspr,mtotlosr
      double precision Lmax,tmax,Mackmax,enucmax,elemsave
      double precision phim,deltam
      double precision vvar(4,nsh)
      double precision eta,degpe,rhpsi
      double precision fscr

      common /saveMS/ elemsave
      common /hydrodyn/ Lmax,tmax,Mackmax,enucmax,itmax,ishockb,ishockt
      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)
      common /difcirc/ dift(nsh),dife(nsh)
      common /calcDh/ Dhold(nsh)
      common /gradmuj/ abmurj(nsh),abmuj(nsh)
      common /metal/ FeH
      common /edegen/ eta(nsh),degpe(nsh),rhpsi(nsh)
      common /funcscr/ fscr(nsh,nscr)
      common /overshoot/ klcore,klenv,klpulse

      dimension y(nsp),xspr(nsh,100),vxspr(nsh,100),mtotlosr(100)

***   read binary models
      
      rewind (uread)

      read (uread) model,nphase,totm,time,dtn,neff,bin_version,nmod,
     &     nspr,mdredgeup,FeH,turnenv,rayon0,
     &     (crzc(k),u(k),r(k),lnf(k),lnt(k),lum(k),
     &     vu(k),vr(k),vlnf(k),vlnt(k),vlum(k),vfacc(k),
     &     exposure(k),vmueinv(k),dm(k),vhconv(k),vomega(k),
     &     ro(k),vro(k),p(k),vp(k),e(k),ve(k),venucl(k),s(k),vs(k),
     &     vsconv(k),vpvisc(k),tau(k),
     &     (xspr(k,l),l = 1,nspr),k = 1,nmod),
     &     (mtotlosr(j),j = 1,nspr)
      if (irotbin.eq.1) then
         read (ureadrot) ndtold,ndbold,ur_Ss,xmom_tots,inits,
     &        (auxs(k),urs(k),vxpsi(k),xpsis(k),V_circ(k),
     &        xalpha(k),xlambdas(k),xKt(k),dift(k),dife(k),
     &        Dhold(k),Dconv(k),abmurj(k),xnuvv(k),xnum(k),k = 1,
     &        nmod),((vxspr(k,l),l = 1,nspr),k = 1,nmod)
      endif
      ierr = .false.
      if (neff.lt.0) then
         ierr = .true.
         neff = abs(neff)
      endif
      
c..   compatibility check
      if (nphase.gt.1.and.mdredgeup.eq.0.d0) mdredgeup = totm*msun
      if (bin_version.lt.2.95d0) then
         write (*,*) '*** Initialize neutron exposure ***'
         forall (i = 1:nmod) exposure(i) = 0.d0
         mdredgeup = totm*msun
      endif
      do k = 1,nmod
         crz(k) = 0
         if (crzc(k).eq.'c') crz(k) = -2 
         if (crzc(k).eq.'a') crz(k) = -1 
         if (crzc(k).eq.'s') crz(k) = -3 
         if (crzc(k).eq.'l'.or.crzc(k).eq.'S') crz(k) = 3 
         if (crzc(k).eq.'t') crz(k) = 2 
         if (crzc(k).eq.'o') crz(k) = 1 
         if (crzc(k).eq.'r') crz(k) = 4 
         if (crz(k).eq.0) then
            write (*,'("Pb shell ",i4," crzc = ",a1," not recognized")')
     &           k,crzc(k)
            stop
         endif
      enddo
      
****************************************************
* check chemical composition of first sequence model
****************************************************

      if (flag.eq.0) then
         inquire (file = 'last_model', exist = test)
         if (test) then
            open (unit = 9, file = 'last_model')
            read (9,*) last_model
            close (9)
            if (model.ne.last_model) then
               write (*,*) 'The previous and new model numbers are ',
     &              'not consecutive : check binary/sequence name'
               write (90,*) 'The previous and new model numbers are ',
     &              'not consecutive : check binary/sequence name'
               stop 'rinimod'
            endif            
         endif
         write (nout,100) zkint,FeH
         write (90,100) zkint,FeH
 100     format (1x,' INITIAL COMPOSITION : Z = ',0pf8.6,' , [Fe/H] = ',
     &        0pf7.4)
      endif
      if (nspr+2.gt.nsp) then
         write (nout,*) 'binary conversion : MASSIVE --> EVOL format'
         do k = 1,nmod
            do l = 1,51
               xsp(k,l) = xspr(k,l)
            enddo
            xar36 = xspr(k,52)
            xsp(k,52) = xspr(k,53)
            xsp(k,53) = xar36
            do l = 54,nspr
               xsp(k,53) = xspr(k,53)+xspr(k,l)
            enddo
            if (irotbin.eq.1) then
               do l = 1,51
                  vxsp(k,l) = vxspr(k,l)
               enddo
               xar36 = vxspr(k,52)
               vxsp(k,52) = vxspr(k,53)
               vxsp(k,53) = xar36
               do l = 54,nspr
                  vxsp(k,53) = vxspr(k,53)+vxspr(k,l)
               enddo
            endif
         enddo
         do l = 1,51
            mtotlos(l) = mtotlosr(l)
         enddo
         xar36 = mtotlosr(52)
         mtotlosr(52) = mtotlosr(53)
         mtotlosr(53) = xar36
         do l = 1,nspr
            mtotlos(53) = mtotlosr(53)+mtotlosr(l)
         enddo
      else
         do l = 1,nspr
            mtotlos(l) = mtotlosr(l)
         enddo
         do l = 1,nis-1
            do k = 1,nmod
               xsp(k,l) = xspr(k,l)
            enddo
         enddo
         if (irotbin.eq.1) then
            do l = 1,nis-1
               do k = 1,nmod
                  vxsp(k,l) = vxspr(k,l)
               enddo
            enddo
         endif
         do l = 1,nis-1
            mtotlos(l) = mtotlosr(l)
         enddo
      endif

      mtotlos(nis) = 1.d-50
      mtotlos(nsp) = 1.d-50
      do k = 1,nmod
         xsp(k,nis) = 1.d-50
         xsp(k,nsp) = 1.d-50
      enddo


***   make correspondance between old and new binary files
      if (bin_version.gt.5.d2.or.bin_version.eq.0.d0) then
         write (nout,*) 'species conversion in binary : NETWORK change'
         do k = 1,nsp
            y(k) = xsp(1,k)
         enddo
         do k = 1,nmod
            xc11 = xsp(k,13)
            xsi31 = xsp(k,40)
            xsi32 = xsp(k,42)
            xp32 = xsp(k,43)
            xp33 = xsp(k,45)
            xsp(k,46) = xsp(k,46)+xp33
            xsp(k,45) = xsp(k,44)+xp32+xsi32
            xsp(k,44) = xsp(k,41)+xsi31
            xsp(k,43) = xsp(k,39) ! SI 30
            xsp(k,42) = xsp(k,38) ! SI 29
            xsp(k,41) = xsp(k,37) ! SI 28
            xsp(k,40) = xsp(k,36) ! AL 27
            xsp(k,39) = 1.d-50    ! Mg 27
            xsp(k,38) = xsp(k,35)
            xsp(k,37) = xsp(k,34)
            xsp(k,36) = xsp(k,33)
            xsp(k,35) = xsp(k,32)
            xsp(k,34) = 1.d-50    ! Na 25
            xsp(k,33) = xsp(k,31)
            xsp(k,32) = 1.d-50    ! Na 24
            xsp(k,31) = xsp(k,30)
            xsp(k,30) = 1.d-50    ! Ne 23
            do j = 13,24
               xsp(k,j) = xsp(k,j+1)
            enddo
            xsp(k,25) = 1.d-50    ! F  20
            xsp(k,12) = xsp(k,12)+xc11
         enddo
         do k = 1,nsp
            write (nout,200) k,k,y(k),k,xsp(1,k)
 200        format(' # :',i2,', old xsp(1,',i2,') = ',
     &           1pe9.3,' --> new xsp(1,',i2,') = ',1pe9.3)
         enddo
      endif
      
c..   initializations and variables definition

      r(1) = 1.d-99
      vr(1) = 1.d-99
      lnr(1) = -227.d0
      vlnr(1) = -227.d0
      enucl(1) = venucl(1)
      m(1) = 0.d0
      do k = 2,nmod
         m(k) = m(k-1)+dm(k-1)
         lnr(k) = log(r(k))
         vlnr(k) = log(vr(k))
         enucl(k) = venucl(k)
      enddo

      write(899,*)'IN RINIMOD'
      do i = 1,nmod
         write(899,*) i,vomega(i)
      enddo

      if (flag.eq.0) then
         rewind (uwrite)
         write (uwrite) model,nphase,totm,time,dtn,neff,bin_version,
     &        nmod,nis-1,mdredgeup,FeH,turnenv,rayon0,
     &        (crzc(k),u(k),r(k),lnf(k),lnt(k),lum(k),
     &        vu(k),vr(k),vlnf(k),vlnt(k),vlum(k),vfacc(k),
     &        exposure(k),vmueinv(k),dm(k),vhconv(k),vomega(k),
     &        ro(k),vro(k),p(k),vp(k),e(k),ve(k),venucl(k),s(k),vs(k),
     &        vsconv(k),vpvisc(k),tau(k),
     &        (xsp(k,l),l = 1,nis-1),k = 1,nmod),
     &        (mtotlos(j),j = 1,nis-1)
         if (irotbin.eq.1) then
            rewind (uwriterot)
            write (uwriterot) ndtold,ndbold,ur_Ss,xmom_tots,inits,
     &           (auxs(k),urs(k),vxpsi(k),xpsis(k),V_circ(k),
     &           xalpha(k),xlambdas(k),xKt(k),dift(k),dife(k),
     &           Dhold(k),Dconv(k),abmurj(k),xnuvv(k),xnum(k),
     &           k = 1,nmod),((vxsp(k,l),l = 1,nis-1),k = 1,nmod)
         endif
                  
c         if (time.gt.4.76d17) stop 'rinimod : age exceeds 15 Gyr !'
      endif

*____________________________
***   various initializations
*----------------------------

c..   stop microscopic diffusion if nphase > 4
      if (nphase.gt.4.and.lmicro.ne.0) then
         lmicro = 0
         microdiffus = .false.
         print *, 'microdiffus become false nphase > 4'
      endif

      nmod1 = nmod-1
      dmaccr = 0.d0
      error = 0
      mr(1) = 0.d0
      forall (k = 1:nmod)
         mr(k) = m(k)/m(nmod)
         t(k) = exp(lnt(k))
         vt(k) = exp(vlnt(k))
      end forall
      

***   initialize network parameters + thermodynamical variables
      if (flag.eq.0) then
         call rininet
         call eos (0,error)

c.. define initial H and He abundances for additional save
         if (nphase.le.2) then
            elemsave = max(xsp(1,ih1)-0.05d0,0.d0)
         elseif (nphase.le.4) then
            elemsave = max(xsp(1,ihe4)-0.05d0,0.d0)
         else
            elemsave = 0.d0
         endif

***   initialize model parameters
         no = 0      ! current model between [1,nmaxmod]
         if (ierr.and.icorr.eq.'t') then
            ifail = -1
         else
            ifail = 0
         endif

         iter = 1
c         if (nphase.eq.1) ishtest = .false.
         if (nphase.eq.2.or.nphase.eq.4.or.nphase.eq.6) ishtest = 't'
         if (nphase.lt.2) then
            if (abs(log(mtini/totm)).gt.log(1.02d0).and.iaccr.eq.0) then
               write (nout,300) totm
 300           format ('change mtini in starevol.par, does not ',
     &              'correspond to stellar mass model = ',0pf8.5)
               stop ' rinimod'
            endif
         endif
c..   check if metallicity table is correct
         Zstart = 1.d0-xsp(nmod,ih1)-xsp(nmod,ih2)-xsp(nmod,ihe3)-
     &        xsp(nmod,ihe4)
         
         if (((abs(log(zstart/zkint)).gt.log(1.15d0).and.Zstart.gt.
     &        1.d-6).or.(Zstart.lt.1.d-6.and.zkint.gt.1.d-6)).and.
     &        nphase.lt.3.and.(.not.microdiffus.or.crzc(1).eq.'c')) then
            write (*,400) zkint,zstart
 400        format (/,2x,'WARNING : Zkint in starevol.par (',f8.6,')'
     &           ' does not correspond to stellar composition ',f8.6,/)
c            if (nphase.le.2) stop ' rinimod'
         endif

***   initialize nuclear luminosities (needed for init)
         totgrav = 0.d0
         totnucl = 0.d0
         vlpp = 0.d0
         vlh = 0.d0
         vlhe = 0.d0
         vlc = 0.d0
         vlne = 0.d0
         vlo = 0.d0
         vlsi = 0.d0
         lpp = 0.d0
         lh = 0.d0
         lhe = 0.d0
         lc = 0.d0
         lne = 0.d0
         lo = 0.d0
         lsi = 0.d0
          
***   initialize thermo (needed for screening factors, rotation)
         forall (l = 1:nsp, k = 1:nmod) ysp(k,l) = xsp(k,l)/anuc(l)
         forall (k = 1:nmod )
            vvar(1,k) = e(k)
            vvar(2,k) = p(k)
            vvar(3,k) = s(k)
            vvar(4,k) = ro(k)
         end forall
         call eos (1,error)
         forall (k = 1:nmod )
            e(k) = vvar(1,k)
            p(k) = vvar(2,k)
            s(k) = vvar(3,k)
            ro(k) = vvar(4,k)
         end forall
c..   initialize screening factors
         call fscreen (t,ro,xsp,rhpsi,anuc,znuc)

      endif

***   initialize vmueinv for URCA process if old binary is used
      if (vmueinv(1).eq.0.d0) then
         write (nout,*) ' WARNING : binary change, compute mueinv'
         write (90,*) ' WARNING : binary change, compute mueinv'
         forall (k = 1:nmod) vmueinv(k) = mueinv(k)
      endif

***   initialize variables for ROTATION
      if (rotation.and.bin_version.lt.2.8d0) then
         write (nout,*) ' WARNING : binary modification, compute PSI'
         write (90,*) ' WARNING : binary modification, compute PSI'
         if (numeric.eq.2.or.numeric.eq.3) then
c.. first order accuracy in the spatial derivatives
            forall (i = 1:nmod)        
               wi(i) = 0.5d0
               wj(i) = 0.5d0
            end forall
         else
c.. second order accuracy in the spatial deriatives
            wi(1) = 0.5d0
            wj(1) = 0.5d0
            forall (i = 2:nmod1)
               wi(i) = dm(i)/(dm(i)+dm(i-1))
               wj(i) = 1.d0-wi(i)
            end forall
            wi(nmod) = 0.5d0
            wj(nmod) = 0.5d0
         endif
         do k = 2,nmod
            deltam = deltaKS(k-1)*wi(k)+deltaKS(k)*wj(k)
            phim = phiKS(k-1)*wi(k)+phiKS(k)*wj(k)
            xpsis(k) = (phim*xlambdas(k)-xpsis(k)*vomega(k)
     &           /vomega(nmod))/deltam
         enddo
         xpsis(1) = xpsis(2)
      endif

***   define novlim (needed for mesh)
      novlim(1:nmaxconv,1:ntypeconv) = 0
      nsconv = 0
      nschwar = 0
      nledoux = 0
      nsemiconv = 0
      nthermha = 0

***   define Schwarzschild and Ledoux boundaries
      kl = 1
      if (crz(1).lt.-1) novlim(kl,1) = 1
      do k = 2,nmod1
         tbot = crz(k).gt.0.and.crz(k+1).lt.-1
         ttop = crz(k).lt.-1.and.crz(k+1).gt.-2
         if (tbot) novlim(kl,1) = k+1
         if (ttop) then
            if (novlim(kl,1).lt.k) then
               novlim(kl,2) = k
               kl = kl+1
               if (kl.gt.nmaxconv) goto 20
            endif
         endif
      enddo
 20   nschwar = kl-1
      if (nschwar.gt.0) then
         if (novlim(nschwar,2).eq.0) then
            write (*,*)'Problem in determining Schwarzschild boundaries'
            stop 
         endif
      endif
      kl = 1
      it = -3
      if (ledouxmlt) it = 3
      if (crz(1).lt.-1) novlim(kl,5) = 1
      do k = 2,nmod1
         tbot = crz(k).ne.-2.and.crz(k+1).eq.-2
         ttop = crz(k).eq.-2.and.crz(k+1).ne.-2
         if (tbot) novlim(kl,5) = k+1
         if (ttop) then
            if (novlim(kl,5).lt.k) then
               novlim(kl,6) = k
               kl = kl+1
               if (kl.gt.nmaxconv) goto 30
            endif
         endif
      enddo
 30   nledoux = kl-1
      if (nledoux.gt.0) then
         if (novlim(nledoux,6).eq.0) then
            write (*,*) 'Problem in determining Ledoux boundaries'
            stop 
         endif
      endif
***   define semiconvective and thermohaline zones
      do j = 1,2
         kl = 1
         jl = 2*j+7
         if (j.eq.2) it = 2
         if (crz(1).eq.it) novlim(kl,jl) = 1
         do k = 2,nmod1
            tbot = crz(k).ne.it.and.crz(k+1).eq.it
            ttop = crz(k).eq.it.and.crz(k+1).ne.it
            if (tbot) novlim(kl,jl) = k+1
            if (ttop) then
               if (novlim(kl,jl).lt.k) then
                  novlim(kl,jl+1) = k
                  kl = kl+1
                  if (kl.gt.nmaxconv) goto 40
               endif
            endif
         enddo
 40      if (j.eq.1) nsemiconv = kl-1
         if (j.eq.2) nthermha = kl-1
      enddo

      if (ledouxmlt) then
         nsconv = nledoux
         novlim(1:nsconv,3) = novlim(1:nsconv,5)
         novlim(1:nsconv,4) = novlim(1:nsconv,6)
      else
         nsconv = nschwar
         novlim(1:nsconv,3) = novlim(1:nsconv,1)
         novlim(1:nsconv,4) = novlim(1:nsconv,2)
      endif
      scr = 'r'
      if (nsconv.gt.0) scr = 'c'

*** define phase
      teff = exp(lnt(neff))
      rgbphase = nphase.eq.3.and.teff.lt.teffagb.and.totm.le.2.5d0
      agbphase = nphase.eq.5.and.teff.lt.teffagb.and.totm.le.10.d0
      superagb = nphase.eq.5.and.teff.lt.teffagb.and.totm.ge.7.d0.and.
     &     totm.le.12.d0.and.xsp(1,ic12).lt.2d-2

***   initialize variables for accretion
      dmaccr = 0.d0
      if (iaccr.ge.1) then
         forall (l = 1:nsp, k = 1:nmod) vxsp(k,l) = xsp(k,l)
         if (numeric.eq.2.or.numeric.eq.3) then
c.. first order accuracy in the spatial derivatives
            forall (i = 1:nmod)        
               wi(i) = 0.5d0
               wj(i) = 0.5d0
            end forall
         else
c.. second order accuracy in the spatial deriatives
            forall (i = 2:nmod1)
               wi(i) = dm(i)/(dm(i)+dm(i-1))
               wj(i) = 1.d0-wi(i)
            end forall
            wi(1) = 0.5d0
            wj(1) = 0.5d0
            wi(nmod) = 0.5d0
            wj(nmod) = 0.5d0
         endif
         ishockb = nmod1
         ishockt = 1
         sigma0 = xiaccr*dsqrt(g*r(nmod)*m(nmod))
         iacctop = 1
         iaccbot = nmod
         do k = 1,nmod
            if (facc(k).gt.0.d0.and.iaccbot.eq.nmod) iaccbot = k
            if (facc(k).gt.0.d0) iacctop = max(iacctop,k)
         enddo
         forall (i = 1:nmod)
            facc(i) = vfacc(i)
            omi(i) = 0.d0
            psi(i) = 0.d0
            ft(i) = 1.d0
            fp(i) = 1.d0
         end forall
         hbb = .false.
         klenv = nsconv
         if (accphase.eq.0) call thermo (error)
      endif
      
      return
      end
