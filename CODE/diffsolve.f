*************************************************************************

      SUBROUTINE diffsolve (abond,ratneh,anuc,znuc,dmix,f1,f2,f3,dmdr,
     &     ielem,nshell,icall,error)

*************************************************************************
*     solve diffusion equation : dX/dt = d/dm(DdX/dt-vdiff*V)
*                                                                      *
* $LastChangedDate:: 2014-02-04 15:45:03 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 11                                                          $ *
*                                                                      *
*************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.diff'
      include 'evolcom.mass'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer, intent (in) :: ielem
      integer ij,ndt,ndb,error
      integer iiter,iitermin,iitermax,ip,im
      integer isolve,icall
      integer i,ijk,nshell,icor,icormin,jshell
      integer mspec             ! Ajout pour form. thoul
      integer elemthoul
      parameter (mspec = 32)    ! Ajout pour form. thoul

      integer ineutron,ih1,ih2,ihe3,ihe4,ili6,ili7,ibe7,ibe9,ib8,ib10, ! Ajout pour form. thoul
     &     ib11,ic12,ic13,ic14,in13,in14,in15,io15,io16,io17,io18,if19,
     &     ine20,ine21,ine22,ina23,img24,img26,ial27,isi28,ip31,is32,
     &	   icl35,idxspc,idxsps,idxspci
      integer ina24,ina25,ine23,ial26g,ial26m,if18,if20,ina22,img27,
     &     is35,icl36,img25

      double precision abond,ratneh,anuc,znuc,tsh,rhosh !Ajout des deux derniers pour form. thoul
      double precision vdiffp,diff13,ration  ! M&M ioni part 09/18
      double precision rray,rkonv,rkonvL,muizero
      double precision, save :: cd(nsh),cdm(nsh)
      double precision xa,d12,dmix,vdiff,vv  !,vdiffthoul !Ajout des deux derniers pour form. thoul
      double precision dmdr,x1,x2,x3,x
      double precision v,w,f1,f2,f3
      double precision tol,alpha,corm,cormin,dx,a1,eqd,dv,cdi,cdip
      double precision vdiffthoul, dmicthoul
      
      logical indiffus

      common /star/ rray(nsh),rkonv(nsh),rkonvL(nsh),muizero(nsh)
      common /thoulindex/ elemthoul(mspec-1)                        ! Ajout TD AP
      common /specindex/ ineutron,ih1,ih2,ihe3,ihe4,ili6,ili7,ibe7,ibe9,                    ! Ajout pour form. thoul
     &     ib8,ib10,ib11,ic12,ic13,ic14,in13,in14,in15,io15,io16,io17
     &     ,io18,if19,ine20,ine21,ine22,ina23,img24,ial27,isi28,ip31,
     &     is32,icl35,ina24,ina25,ine23,img26,ial26g,ial26m,if18,if20,
     &     ina22,img27,is35,icl36,img25,idxspc(nprintc),idxsps(nprint)
     &	   ,idxspci(9)
      common /dmicthoul/ vdiffthoul(nsh,mspec), dmicthoul(nsh,mspec) ! modif

      dimension abond(nsh,nis),xa(nshell),x(nsh),ratneh(nsh),
     &     anuc(nsp),znuc(nsp)      
      dimension d12(nsh),dmix(nsh),eqd(4)
      dimension vdiff(nsh),vv(nshell)
      dimension v(nshell),w(nshell),f1(nsh),f2(nsh),f3(nsh),dmdr(nsh)

      isolve = 1
      ij = ielem
      tol = 1.d-14
      alpha = 0.75d0
      iitermax = 200

c..   initialisations
      if (ielem.eq.2.or.(ielem.eq.3.and.microdiffus)) then
         dmix(nmod) = dmix(nmod1)
         do i = 1,nmod1
            ip = i+1
            cd(i) = dmix(i)
            cdm(ip) = dmix(i)**wi(ip)*dmix(ip)**wj(ip)
         enddo
      endif
         
 100  iiter = 0
      iitermin = 3
      ndb = 1
      ndt = nmod-1
      indiffus = .false.
      icormin = 1
      cormin = 1.d99
      do i = 1,nmod1
         vv(i) = 0.d0
         xa(i) = abond(i,ielem)
         x(i) = xa(i)
      enddo

*_____________________________
***   Newton Raphson procedure
*-----------------------------

  200 continue
      icor = 1
      corm = 0.d0
      iiter=iiter+1

*_____________________________________________
***   Compute microscopic diffusion velocities
*---------------------------------------------

      if (microdiffus.and.icall.eq.1) then
         indiffus = .true.
c..   Diffusion "Chapman & Cowling"
         if (lmicro.eq.2) then
            call diffmic (ij,ndt,ndb,vdiff,d12,anuc,znuc)
c..   Diffusion "Paquette et al. (1986)"
         else if (lmicro.eq.3.or.lmicro.eq.5) then
            call diffpaq (ij,ndt,ndb,vdiff,d12,abond,ratneh)
c..   Diffusion "Thoul et al. (1994)"
         else if (lmicro.eq.4.or.lmicro.eq.6) then
            do i = 2,nmod1            
               if (crz(i).gt.0) then
                  cd(i) = dmix(i)
               endif
               do ijk = 1,mspec-1
                  if (ij.eq.elemthoul(ijk)) then
                     vv(i)=vdiffthoul(i,ijk)*dmdr(i)
                  endif
               enddo
               cdm(i) = cd(i-1)**wi(i)*cd(i)**wj(i)
            enddo
         endif
         Dmicro(1:nsh) = d12(1:nsh)
         vmicro(1:nsh) = vdiff(1:nsh)
         if (lmicro.eq.2.or.lmicro.eq.3.or.lmicro.eq.5) then
            do i = 2,nmod1
               if (crz(i).gt.0) then
                  cd(i) = dmix(i)+d12(i)*dmdr(i)*dmdr(i)
                  vv(i) = vdiff(i)*dmdr(i)
               endif
               cdm(i) = cd(i-1)**wi(i)*cd(i)**wj(i)
            enddo
         endif
      endif

c..   center
      cdip = cdm(2)
      x1 = f1(1)*vv(2)
      x2 = f2(1)*cdip
      eqd(4) = x(1)-xa(1)-x2*(x(2)-x(1))+x1*x(2)
      eqd(1) = 0.d0
      eqd(2) = 1.d0+x2
      eqd(3) = -x2+x1
      if (eqd(2).eq.0.d0) then
         error = 20
         return
      endif
      w(1) = eqd(4)/eqd(2)
      v(1) = eqd(3)/eqd(2)

c..   interior
      do i = 2,nmod1-1
         ip = i+1
         im = i-1
         cdi = cdm(i)
         cdip = cdm(ip)
         x3 = f3(i)*cdi
         x2 = f2(i)*cdip
         eqd(4) = x(i)-xa(i)-x2*(x(ip)-x(i))+x3*(x(i)-x(im))+
     &        f1(i)*(vv(ip)*x(ip)-vv(i)*x(i))
         eqd(1) = -x3
         eqd(2) = 1.d0+x3+x2-f1(i)*vv(i)
         eqd(3) = -x2+f1(i)*vv(ip)
         if (isolve.eq.0) then
            a1 = eqd(2)-eqd(1)*v(im)
            if (a1.eq.0.d0) a1 = -eqd(3)-eqd(1)*(1.d0+v(im))+1.d0
         elseif (isolve.eq.1) then
            a1 = -eqd(3)-eqd(1)*(1.d0+v(im))+1.d0
            if (a1.eq.0.d0) a1 = eqd(2)-eqd(1)*v(im)
         elseif (isolve.eq.2) then
            a1 = eqd(2)-eqd(1)*v(im)
            a1 = 0.5d0*(a1-eqd(3)-eqd(1)*(1.d0+v(im)))+0.5d0
         endif
         if (a1.eq.0.d0) then
            iiter = iitermax
            goto 300
         endif
         dv = eqd(1)*w(im)-eqd(4)
         w(i) = -dv/a1
         v(i) = eqd(3)/a1
      enddo

c..   surface
      i = nmod1
      cdi = cdm(i)
      x1 = f1(i)*vv(i)
      x3 = f3(i)*cdi
      eqd(4) = x(i)-xa(i)+x3*(x(i)-x(i-1))-x1*x(i)
      eqd(1) = -x3
      eqd(2) = 1.d0+x3-x1
      eqd(3) = 0.d0
      a1 = eqd(2)-eqd(1)*v(i-1)
      if (a1.eq.0.d0) a1 = -eqd(3)-eqd(1)*(1.d0+v(i-1))+1.d0
      if (a1.eq.0.d0) then
         write (*,'("a1 = 0 in diffsolve at surface #",i4," for species"
     &        ,i3,", nshell =",i4," ["i2"]")') i,ielem,NSHELL,crz(i)
         write (*,'("shell",5x,"X",8x,"Xold",8x,"VV",9x,"cd",9x,"f1",
     &        9x,"f2",9x,"f3")')
         error = 22
         stop
      endif
      dx = (eqd(1)*w(i-1)-eqd(4))/a1

      x(i) = x(i)+alpha*dx
      if (x(i).lt.0.d0) x(i) = 0.5d0*(x(i)-alpha*dx)

      if (x(i).gt.10.d0*tol) then
         corm = dx/x(i)
      else
         corm = dx
      endif

      do i = nmod1-1,1,-1
         dx = -w(i)-v(i)*dx
         x(i) = x(i)+alpha*dx
         if (x(i).lt.0.d0) x(i) = 0.5d0*(x(i)-alpha*dx)

         if (x(i).gt.10.d0*tol) then
            if (abs(dx/x(i)).gt.abs(corm)) then
               icor = i
               corm = dx/x(i)
            endif
         else
            if (abs(dx).gt.abs(corm)) corm = dx
         endif
      enddo
      if (abs(corm).le.cormin) then
         cormin = abs(corm)
         icormin = icor
      endif

      if (indiffus) then
         x(nmod) = x(nmod1)
         do i = 1,nmod1
            abond(i,ielem) = x(i)
         enddo
      endif

 300  if (iiter.ge.iitermax) then
c.. change convergence acceleration
         if (isolve.eq.1) then
            alpha = 0.2d0
            iitermax = 2000
            isolve = 0
            goto 100
         endif
c.. change resolution of equation ? don't remember the origin
         if (isolve.eq.0) then
            isolve = 2
            goto 100
         endif
c.. no convergence
         write (nout,'(3x,"Convergence problem species #",i2)') ielem
         error = 23
         return
      endif

      if (crz(icor).ne.4.and.(abs(corm).gt.tol.or.iiter.lt.iitermin))
     &     goto 200

      if (.not.indiffus) then
         x(nmod) = x(nmod1)
         do i = 1,nmod1
            abond(i,ielem) = x(i)
         enddo
      endif


      return
      end
