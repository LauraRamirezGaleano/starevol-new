

************************************************************************

      SUBROUTINE netdiff (imin,imax,cd,dtn,error)

************************************************************************
* Compute the abundance evolution resulting from the coupling of the   *
* nucleosynthesis and diffusion equations                              *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eos'
      include 'evolcom.ion'
      include 'evolcom.therm2'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.spec'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer i,j,k,l,imin,imax,ielem,error
      integer iterd,neqd,neqd1,neqd2,mm,id,i1
      integer ncmmax,idmax,ntry,istep,itermax

      logical nequil,xequil

      double precision vv(nsh,nreac),vi(nreac)
      double precision f1(nsh),f2(nsh),f3(nsh),dd(nsh),cd(nsh),dmdr(nsh)
      double precision alphad,xx,xidk,xtol
      double precision eqd,xmod,xa,y,yeq
      double precision dx,dxmax,dv,w,a,um,v,alpha1,cormmax
      double precision xm_sum,dt
      double precision tk,rok,dtn,dtnuc,dtold,fac
      
      dimension eqd(nsp,3*nsp+1),xmod(nsh,nsp),xa(nsh,nsp),y(nsp)
      dimension dx(nsp),dv(nsp),w(nsh,nsp),a(nsp,2*nsp),um(nsp,nsp),
     &     v(nsh,nsp,nsp),alpha1(nsp)

      xequil = .false.
      nequil = .true.

      itermax = 1000
      icrash = 0
      ntry = 0
      istep = 0
      dt = 0.d0
      dtnuc = dtn*1.d-1

c..   computation nuclear reaction rates in all shells
      do k = imin,imax
         tk = 0.5d0*(vt(k)+t(k))
         rok = 0.5d0*(vro(k)+ro(k))
         call vit (tk,rok,mueinv(k),k,vi,0)
         do j = 1,nreac
            vv(k,j) = vi(j)
         enddo
      enddo

c..   initialisations
      alphad = 1.d0
      xtol = 1.d-8*tolnuc
      neqd = nsp
      neqd2 = 2*neqd
      neqd1 = 3*neqd+1

      do j = 1,nsp
         do i = imin,imax
            xmod(i,j) = xsp(i,j)/anuc(j)
            xa(i,j) = xmod(i,j)
         enddo
      enddo
      do i = 1,neqd
         alpha1(i) = alphad
      enddo

      do k = imin,imax
         if (k.eq.1) then
            dmdr(1) = 0.d0
            f1(1) = dtnuc/dm(1)
            f2(1) = 2.d0*f1(1)/(dm(1)+dm(2))
            f3(1) = 0.d0
         else
            dmdr(k) = pim4*r(k)*r(k)*rom(k)
            f1(k) = dtnuc/dm(k)
            f2(k) = 2.d0*f1(k)/(dm(k)+dm(k+1))
            f3(k) = 2.d0*f1(k)/(dm(k)+dm(k-1))
         endif
         dd(k) = cd(k)*dmdr(k)*dmdr(k)
      enddo
      dd(1) = dd(2)
      do k = imax-1,imin,-1
         dd(k+1) = dd(k)**wi(k+1)*dd(k+1)**wj(k+1)
      enddo
      dd(1) = dd(2)

*_______________________________________________________________________
***   calculation of the coupled nucleosynthesis and diffusion equations
***   Newton Raphson procedure
*-----------------------------------------------------------------------

      dtold = dtnuc
      iterd = 0
 10   iterd = iterd+1
      if (iterd.ge.itermax.or.dtnuc.lt.1.d-10) goto 20


c..   central boundary conditions
*--------------------------------

      id = imin
      do l = 1,neqd1
         do k = 1,neqd
            eqd(k,l) = 0.d0
         enddo
      enddo
      do j = 1,neqd
         y(j) = xmod(id,j)
      enddo
      do j = 1,nreac
         vi(j) = vv(id,j)
      enddo
c..   set neutrons to equilibrium abundance
      if (nequil) then
         call yneutron (y,vi,yeq)
         y(1) = yeq
      endif
c..   set  equilibrium abundance for elements with y < yeqmin
      if (xequil) call yequil (vi,y,ytminc,dtnuc,tconv(id))

      call freac (vi,eqd,y,dtnuc,.true.)
      if (kcentd.ne.3) call fdiff (eqd,dd,xmod,f1,f2,f3,id,imin,imax)

      do j = 1,neqd
         if (kcentd.ne.3) then
            eqd(j,neqd1) = eqd(j,neqd1)-xmod(id,j)+xa(id,j)
            eqd(j,j) = eqd(j,j)-1.d0
         else
            eqd(j,neqd1) = -xmod(id,j)+xmod(id-1,j)
            eqd(j,j) = -1.d0
         endif
         do l = 1,neqd
            a(l,j) = eqd(l,j)
            a(l,neqd+j) = 0.d0
         enddo
         a(j,neqd+j) = 1.d0
      enddo
      call nwrmat (a,um,neqd,neqd,error)
      if (error.gt.0) then
         error = 9
         return
      endif
      do mm = 1,neqd
         w(id,mm) = 0.d0
         do k = 1,neqd
            w(id,mm) = w(id,mm)+um(mm,k)*eqd(k,neqd1)
            v(id,mm,k) = 0.d0
            do j = 1,neqd
               v(id,mm,k) = v(id,mm,k)+um(mm,j)*eqd(j,neqd+k)
            enddo
         enddo
      enddo

      do l = 1,neqd1
         do k = 1,neqd
            eqd(k,l) = 0.d0
         enddo
      enddo

c..   internal shell diffusion equations
*---------------------------------------

      do id = imin+1,imax-1
         i1 = id-1

         do j = 1,neqd
            y(j) = xmod(id,j)
         enddo
         do j = 1,nreac
            vi(j) = vv(id,j)
         enddo
c..   set neutrons to equilibrium abundance
         if (nequil) then
            call yneutron (y,vi,yeq)
            y(1) = yeq
         endif
c..   set  equilibrium abundance for elements with y < yeqmin
         if (xequil) call yequil (vi,y,ytminc,dtnuc,tconv(id))

         call freac (vi,eqd,y,dtnuc,.false.)
         call fdiff (eqd,dd,xmod,f1,f2,f3,id,imin,imax)
         
         do j = 1,neqd
            eqd(j,neqd1) = eqd(j,neqd1)-xmod(id,j)+xa(id,j)
            eqd(j,neqd+j) = eqd(j,neqd+j)-1.d0
            do l = 1,neqd
               a(j,l) = eqd(j,neqd+l)
               do k = 1,neqd
                  a(j,l) = a(j,l)-eqd(j,k)*v(i1,k,l)
               enddo
               a(j,neqd+l) = 0.d0
            enddo
            a(j,neqd+j) = 1.d0
         enddo
         call nwrmat (a,um,neqd,neqd,error)
         if (error.gt.0) then
            error = 9
            return
         endif
         do l = 1,neqd
            dv(l) = -eqd(l,neqd1)
            do k = 1,neqd
               dv(l) = dv(l)+eqd(l,k)*w(i1,k)
            enddo
         enddo
         do k = 1,neqd
            w(id,k) = 0.d0
            do l = 1,neqd
               w(id,k) = w(id,k)-um(k,l)*dv(l)
               v(id,k,l) = 0.d0
               do j = 1,neqd
                  v(id,k,l) = v(id,k,l)+um(k,j)*eqd(j,neqd2+l)
               enddo
            enddo
         enddo
         do k = 1,neqd
            do l = neqd+1,neqd2
               eqd(k,l) = 0.d0
            enddo
            eqd(k,neqd1) = 0.d0
            eqd(k,k) = 0.d0
            eqd(k,k+neqd2) = 0.d0
         enddo
      enddo

c..   surface boundary conditions
*--------------------------------

      id = imax
      i1 = id-1

      do j = 1,neqd
         y(j) = xmod(id,j)
      enddo
      do j = 1,nreac
         vi(j) = vv(id,j)
      enddo

      call freac (vi,eqd,y,dtnuc,.false.)
      if (ksurfd.ne.3) call fdiff (eqd,dd,xmod,f1,f2,f3,id,imin,imax)

      do j = 1,neqd
         if (ksurfd.ne.3) then
            eqd(j,neqd1) = eqd(j,neqd1)-xmod(id,j)+xa(id,j)
            eqd(j,neqd+j) = eqd(j,neqd+j)-1.d0
         else
            eqd(j,neqd1) = -xmod(id,j)+xmod(id+1,j)
            eqd(j,neqd+j) = -1.d0
         endif
         do l = 1,neqd
            a(j,l) = eqd(j,neqd+l)
            do k = 1,neqd
               a(j,l) = a(j,l)-eqd(j,k)*v(i1,k,l)
            enddo
            a(j,neqd+l) = 0.d0
         enddo
         a(j,neqd+j) = 1.d0
      enddo
      call nwrmat (a,um,neqd,neqd,error)
      if (error.gt.0) then
         error = 9
         return
      endif
      do j = 1,neqd
         dv(j) = -eqd(j,neqd1)
         do k = 1,neqd
            dv(j) = dv(j)+eqd(j,k)*w(i1,k)
         enddo
      enddo

*______________________________________________________
***   determination of the abundance profile variations
***   (iterative process)
*------------------------------------------------------

      do mm = 1,neqd
         dx(mm) = 0.d0
         do j = 1,neqd
            dx(mm) = dx(mm)+um(mm,j)*dv(j)
         enddo
         xmod(id,mm) = xmod(id,mm)+alpha1(mm)*dx(mm)
      enddo

      cormmax = 0.d0
      ncmmax = 1
      idmax = 2

      do id = imax-1,imin,-1
         do k = 1,neqd
            dv(k) = -w(id,k)
            do l = 1,neqd
               dv(k) = dv(k)-v(id,k,l)*dx(l)
            enddo
         enddo
         do k = 1,neqd
            dx(k) = dv(k)
            xmod(id,k) = xmod(id,k)+alpha1(k)*dx(k)
            xidk = max(abs(xa(id,k)),abs(xmod(id,k)))
            xx = abs(dx(k))/xidk
            if (abs(xmod(id,k)).gt.xtol.and.abs(xa(id,k)).gt.xtol) then
               if (abs(xx).gt.abs(cormmax)) then
                  cormmax = xx
                  dxmax = dx(k)
                  ncmmax = id
                  idmax = k
c                  write (*,'(i4,i3,1x,a5,5(1x,1pe10.3))') id,idmax,
c     &               elem(k),xmod(id,k),xa(id,k),dx(k),cormmax,xx
               endif
            endif
         enddo
      enddo
cc..  check convergence
      if (abs(cormmax).gt.rprecd) then
         icrash = icrash+1
         if (mod(iterd,2).eq.0) then
            dtold = dtnuc
            dtnuc = dtnuc*0.33d0
         endif
         write(nout,'("shell",i4,", dt=",1pe10.4,1x,a5,", dY = ",
     &        1pe13.6,", Y =",1pe13.6,", dY/Y =",1pe13.6)') 
     &        ncmmax,dtnuc,elem(idmax),dxmax,xa(ncmmax,idmax),
     &        cormmax
         fac = dtnuc/dtold
         if (mod(iterd,2).eq.0) then
            do i = imin,imax
c                     do j = 1,nsp
c                        xmod(i,j) = xa(i,j)
c                     enddo
               f1(i) = f1(i)*fac
               f2(i) = f2(i)*fac
               f3(i) = f3(i)*fac
            enddo
         endif
         goto 10
      endif

c..   solution accepted, update abundances and increase time step 
      dt = dt+dtnuc
      if (dt.lt.dtn) then
         istep = istep+1
         dtold = dtnuc
         dtnuc = min(dtnuc*1.78d0,dtn-dt)
c         dtnuc = min(dtnuc*max(1.58d0,cormmax),dtn-dt)
         fac = dtnuc/dtold
         write (nout,100) min(istep,999),cormmax,elem(idmax),ncmmax,
     &        xa(ncmmax,idmax)*anuc(idmax),xmod(ncmmax,idmax)*
     &        anuc(idmax),icrash,dtnuc,1.d2*dt/dtn
         do i = imin,imax
            xmod(i,nis-1) = xmod(i,nis-1)+xmod(i,nis)*anuc(nis)/
     &           anuc(nis-1)
            xmod(i,nis) = 1.d-55
            xmod(i,nsp) = 1.d-55
            do j = 1,nsp
               xa(i,j) = xmod(i,j)
            enddo
            f1(i) = f1(i)*fac
            f2(i) = f2(i)*fac
            f3(i) = f3(i)*fac
         enddo
         iterd = 0
         icrash = 0
         goto 10
      else
c..   time step completed, return               
         do j = imin,imax
            xm_sum = 0.d0
            ielem = 2
            xmod(j,nis-1) = xmod(j,nis-1)+xmod(j,nis)*anuc(nis)/
     &           anuc(nis-1)
            xmod(j,nis) = 1.d-55
            xmod(j,nsp) = 1.d-55
            do k = 1,nsp
               xsp(j,k) = max(xmod(j,k)*anuc(k),1.d-50)
               vxsp(j,k) = xsp(j,k)
               xm_sum = xm_sum+xsp(j,k)
               if (xsp(j,k).gt.xsp(j,ielem)) ielem = k
            enddo
            if (xsp(j,ielem).lt.0.d0) then
               error = 25
               return
            endif            
            if (abs(1.d0-xm_sum).gt.xtol) then
               error = 26
               return
            endif            
            xsp(j,ielem) = xsp(j,ielem)+(1.d0-xm_sum)
            if (.not.nequil) then
               do l = 1,neqd
                  y(l) = xmod(j,l)
               enddo
               do i = 1,nreac
                  vi(i) = vv(j,i)
               enddo
               call yneutron (y,vi,yeq)
               xsp(j,1) = yeq
               vxsp(j,1) = yeq
            endif
            vxsp(j,ielem) = xsp(j,ielem)
         enddo
         return
      endif


***   Manage crash situations

 20   if (ntry.lt.2.and.error.eq.0) then
         ntry = ntry+1
c..   1st try : set equilibrium abundance
         xequil = .true.
c..   2nd try : set free neutron abundance
         if (ntry.eq.2) nequil = .not.nequil
         if (nequil) then
            error = 27
            return
         endif
         write (nout,*) iterd,dtnuc,xequil,nequil
         iterd = 0
         dtold = dtnuc
         dt = 0.d0
         dtnuc = dtn*0.125d0
         fac = dtnuc/dtold
         do k = imin,imax
            f1(k) = f1(k)*fac
            f2(k) = f2(k)*fac
            f3(k) = f3(k)*fac
            do j = 1,nsp
               xa(i,j) = xmod(i,j)
            enddo
         enddo
         goto 10
      else
         error = 27
         return
      endif


 100  format(i3,"| dY/Y=",1pe10.3,", for ",a5,", shell =",i4,", X0 =",
     &     1pe10.3,", X =",1pe10.3,", ifail=",i3,", dt=",1pe9.3,
     &     ", time=",0pf10.7,' %')

      return
      end
