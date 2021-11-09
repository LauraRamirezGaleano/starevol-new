

************************************************************************

      SUBROUTINE nucsolve (y,v,dtn,irrad0,nsink,ksh,error)

************************************************************************
*     Solve the nucleosynthesis equations
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.nuc'

      integer i,j,l
      integer istart,error,normx,negflag
      integer neqd,neqd1,ksh,ntry
      integer iterd,itermax
      integer indx(nis),info

      double precision d,suma,sumd,sumx
      double precision yeq,dtn,dtnuc,dt,irrad0
      double precision v(nreac),eqd(nsp,3*nsp+1)
      double precision y(nsp),x(nsp),y0(nsp),ysav(nsp)
      double precision p(nis),fjac(nis,nis)

      integer nz,lirn,ina,jna,licn
      parameter (nz = 370,lirn=nsp,licn=nz+1)
c      integer iw(8*nis),indx_s(5*nis)
c      double precision fjac_s(nz),w(nis)

      logical nequil,xequil,nsink,normalize

      common /nuc_matrix/ jna(nz),ina(lirn)

      ntry = 0
      dt = 0.d0
      if (irrad0.lt.1.d-15) then
         dtnuc = dtn
      else
         dtnuc = dtn*5.d-2
      endif
      irrad0 = 0.d0
      itermax = 10000
      neqd1 = 3*nsp+1
      neqd = nis

      normalize = .false.
      xequil = .false.
      nequil = .true.
      yeq = 1.d-50

c..   for extremely metal poor stars (Z = 0), no neutrons sink
c      if (.not.nsink) nequil = .false.

      if (error.eq.-1) then
         nequil = .false.
         xequil = .true.
      endif
      ysav(1:nsp) = y(1:nsp)
      y0(1:nsp) = y(1:nsp)
      eqd(1:nsp,1:neqd1) = 0.d0
      if (nequil) then
         istart = 2
      else
         istart = 1
      endif

c..   set neutrons to equilibrium abundance
      if (nequil) then
         call yneutron (y,v,yeq)
         y(1) = yeq
      endif

      iterd = 0
 10   iterd = iterd+1
      if (iterd.ge.itermax.or.dtnuc.lt.1.d-10) goto 20

c..   set  equilibrium abundance for elements with y < yeqmin
      if (xequil) call yequil (v,y,ytminc,dtnuc,1.d99)
      call freac (v,eqd,y,dtnuc,.true.)
      do j = 1,neqd
c.withp         p(j) = y(j)-y0(j)-eqd(j,neqd1)
         p(j) = -eqd(j,neqd1)
         do l = 1,neqd
            fjac(l,j) = eqd(l,j)
         enddo
         fjac(j,j) = fjac(j,j)-1.d0
      enddo
c.. leqs
c      call leqs(fjac,p,neqd,neqd)
c.. numerical recipee
c      call ludcmp(fjac,neqd,neqd,indx,d)
c      call lubksb(fjac,neqd,neqd,indx,p)
c.. lapack
cx!$OMP PARALLEL PRIVATE(info)
      call dgetrf (neqd,neqd, fjac, neqd, indx, info)
      call dgetrs ('N',neqd,1,fjac,neqd,indx,p,neqd,info)
cx!$OMP END PARALLEL
c..  sparse matrix
c      p = 1.d0
c      j = 1
c      do i = 1,nz
c         if (i.eq.ina(j+1)) j = j+1
c         print *,i,j,jna(i)
c         fjac_s(i) = fjac(j,jna(i))
c      enddo
c      call ma28ad(neqd,nz,fjac_s,nz,ina,lirn,jna,d,indx_s,iw,w,info)
c      call ma28cd(neqd,fjac_s,nz,jna,indx_s,p,w,1)
c.. gift
c      do j = 1,neqd
c         fjac(j,neqd+1) = p(j)
c      enddo
c      call gift_test(fjac,neqd,neqd+1)
c      do j = 1,neqd
c         p(j) = fjac(j,neqd+1)
c      enddo
      sumd = 0.d0
      suma = 0.d0
      do i = istart,neqd
         sumd = sumd + p(i)*anuc(i) 
         suma = suma + abs(p(i))*anuc(i)
c         write (*,'(i4,4(1x,1pe11.3))') i,suma,sumd,p(i)*anuc(i),y(i)
      enddo
c      stop
      negflag = 0
      do i = istart,neqd
         y(i) = y(i)+p(i)
         if (y(i).lt.-2.d16) negflag = i
c..  check linearization approximation is valid : dY/dt * dtnuc << Y
         if ((abs(p(i)/y0(i)).gt.0.1d0.and.y(i).gt.ytminc.and.
     &        y0(i).gt.ytminc).or.negflag.gt.0.or.suma.ge.5.d-2.or.
     &        abs(sumd).ge.1.d-10) then
            dtnuc = dtnuc*0.23d0
            do l = 1,nsp
               y(l) = y0(l)
               do j = 1,nsp
                  eqd(j,l) = 0.d0
               enddo
               eqd(l,neqd1) = 0.d0
            enddo
            goto 10
         endif
      enddo

c..   solution accepted, update abundances and increase time step 
      dt = dt+dtnuc
      if (dt.lt.dtn) then
         do i = 1,nsp
            y0(i) = max(y(i),1.d-50)
            do j = 1,nsp
               eqd(j,i) = 0.d0
            enddo
            eqd(i,neqd1) = 0.d0
         enddo
         if (normalize) then
            sumx = 0.d0
            normx = 2
            do i = 1,nsp
               x(i) = max(anuc(i)*y(i),1.d-50)
               if (x(i).gt.x(normx)) normx = i
               sumx = sumx+x(i)
               y(normx) = (x(normx)+1.d0-sumx)/anuc(normx)
               y0(normx) = y(normx)
            enddo
         endif
c..   set neutrons to equilibrium abundance
         if (nequil) then
            call yneutron (y,v,yeq)
            y(1) = yeq
         endif
         irrad0 = irrad0+dtnuc*yeq
         dtnuc = min(dtnuc*1.58d0,dtn-dt)
         goto 10
      else
c..   time step completed, return
         if (nsink) then
            call yneutron (y,v,yeq)
            y(1) = yeq
         endif
         irrad0 = irrad0+dtnuc*yeq
         error = 0
         return
      endif


***   Manage crash situations

 20   if (ntry.lt.2.and.error.eq.0) then
         ntry = ntry+1
c..   1st try : set equilibrium abundance
         xequil = .true.
c..   2nd try : set free neutron abundance
         if (ntry.eq.2) nequil = .not.nequil
         if (.not.nsink.and.nequil) then
            error = 14
            return
         endif
         if (nequil) then
            istart = 2
         else
            istart = 1
         endif
         write (nout,100) ksh,iterd,dtnuc,xequil,nequil
         iterd = 0
         dt = 0.d0
         dtnuc = dtn*0.125d0
         y(1:nsp) = ysav(1:nsp)
         eqd(1:nsp,1:nsp) = 0.d0
         eqd(1:nsp,neqd1) = 0.d0
         goto 10
      else
         error = 14
         return
      endif

 100  format (' NUCSOLVE problem at shell #',i4,' : performed ',i5,
     &     ' iterations, final time step = ',1pe9.3,
     &     '(s), new try with Xequil = ',l1,', Nequil = ',l1)

      return
      end
