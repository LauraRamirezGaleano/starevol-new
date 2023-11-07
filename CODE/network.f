

************************************************************************

      SUBROUTINE network (cshell,xi,mconv,dt,tnucmin,imin,imax,error)

************************************************************************
* Calculate the chemical evolution during the current model time-step  *
* cshell = c   : convective zone
* cshell = r   : radiative zone
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.ion'
      include 'evolcom.mod'
      include 'evolcom.nuc'
c      include 'evolcom.conv'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer error,imin,imax,ish
      integer normx,k
      integer i,j,l

      double precision tk,rok,mconv
      double precision xi,dt,y0,ysav,v,vi,irrad0,
     &     yeqmean,dmam,zlocal
      double precision dnucc,yeq,sumxsp,sen,tmean,tnucmin

      logical ipass,nsink

      character cshell*1

      dimension xi(nsp),y0(nsp),ysav(nsp),v(nreac),vi(nreac)

      ipass = .false.
      nsink = .true.

***   special treatment of neutrons for extremely metal poor stars
      if (zkint.lt.5.d-6) then
         zlocal = 0.d0
         do i = 2,ib11
            zlocal = zlocal+xi(i)
         enddo
         zlocal = 1.d0-zlocal
         if (zlocal.lt.1.d-10) then
            k7(iddn) = 2        !  2D(g,n)He3   --> 2D(g,p)He3
            k7(ili7d) = 2       !  Li7(D,n)2He4 --> Li7(D,p)2He4
            nsink = .false.
         else
            k7(iddn) = 1        !  2D(g,n)He3
            k7(ili7d) = 1       !  Li7(D,n)2He4
         endif
      endif

***   set up nuclear network tolerence in advanced evolutionary stages
      if (cshell.eq.'r') then
         ish = imin
         tmean = t(ish)
      else
         tmean = dsqrt(t(imin)*t(imax))
      endif
      if (tmean.gt.4.d8) then
         if (tmean.gt.1.2d9) then
            dnucc = min(1.d-5,dnuc*1.d4)
            ytminc = min(1.d-7,ytmin*1.d4)
         else
            dnucc = min(1.d-7,dnuc*1.d2)
            ytminc = min(1.d-9,ytmin*1.d2)
         endif
      else
         dnucc = dnuc
         ytminc = ytmin
      endif

***   initialisations

      irrad0 = 0.d0
      do l = 1,nsp
c         x(l) = xi(l)
         y0(l) = xi(l)/anuc(l)
         ysav(l) = y0(l)
      enddo

*_____________________________________________________________________
***   determine neutron equilibrium abundance, temperature and density
*---------------------------------------------------------------------

***   radiative shell
      if (cshell.eq.'r') then
         tk = (t(ish)+vt(ish))*0.5d0
         rok = (ro(ish)+vro(ish))*0.5d0
         call vit (tk,rok,mueinv(ish),ish,v,0)
         if (nsink) then 
            call yneutron (y0,v,yeq)
            y0(1) = yeq
            ysav(1) = yeq
         endif
         irrad0 = 7.68771025912336d0*yeq*rok*dsqrt(tk)
      endif

***   convective shell
      if (cshell.eq.'c') then
         do i = 1,nreac
            v(i) = 0.d0
         enddo
         yeqmean = 1.d-50
         do k = imin,imax
            tk = (vt(k)+t(k))*0.5d0
            rok = (vro(k)+ro(k))*0.5d0
            dmam = dm(k)*mconv
            call vit (tk,rok,mueinv(k),k,vi,0)
            if (nphase.le.2.or.tk.gt.tnucmin) then
               if (nsink) then 
                  call yneutron (y0,vi,yeq)
                  y0(1) = yeq
               endif
               irrad0 = irrad0+7.68771025912336d0*yeq*rok*dsqrt(tk)*dmam
               do j = 1,nreac
                  if (k1(j).eq.1.or.k3(j).eq.1) then
                     v(j) = v(j)+vi(j)*dmam*y0(1)
                  else
                     v(j) = v(j)+vi(j)*dmam
                  endif
               enddo
            else
               do i = 1,nbeta
                  j = jbeta(i)
                  v(j) = v(j)+vi(j)*dmam
               enddo
            endif
            yeqmean = yeqmean+y0(1)*dmam
         enddo
         y0(1) = yeqmean
         ysav(1) = y0(1)
         do j = 1,nreac
            if (k1(j).eq.1.or.k3(j).eq.1) v(j) = v(j)/yeqmean
         enddo
      endif
      irrad0 = irrad0*dt

*_______________________________________
***   calculation of the nucleosynthesis
*---------------------------------------

 10   call nucsolve (y0,v,dt,irrad0,nsink,imin,error)
      if (error.gt.0) return

      y0(nis-1) = y0(nis-1)+y0(nis)*anuc(nis)/anuc(nis-1)
      y0(nis) = 1.d-55
      y0(nsp) = 1.d-55
      sumxsp = 0.d0
      normx = 1
      do i = 1,nsp
         xi(i) = max(anuc(i)*y0(i),1.d-50)
         if (xi(i).gt.xi(normx)) normx = i
         sumxsp = sumxsp+xi(i)
      enddo
      sen = 1.d0-sumxsp
      if (abs(sen).gt.dnucc) then
         if (.not.ipass) then
            do l = 1,nsp
               y0(l) = ysav(l)
            enddo
            ytminc = ytminc*1.d2
            dnucc = dnucc*1.d1
            error = -1
            ipass = .true.
            goto 10
         else
            write (nout,1100) imin,sen,dnucc,ytminc,tk*1.d-9
            error = 3
         endif
      else
c..   convergence and normalization
         xi(normx) = xi(normx)+sen
c..   compute neutron irradiation
c..   irrad = y0(1)*rok*avn*dt*dsqrt(2.d0*boltz*tk*avn*0.98d0)*1.d-27
         if (cshell.eq.'c') then
            do k = imin,imax
               exposure(k) = exposure(k)+7.68771025912336d0*ro(k)*
     &              dsqrt(t(k))*irrad0
            enddo
         else
            exposure(imin) = exposure(ish)+7.68771025912336d0*ro(ish)*
     &           dsqrt(t(ish))*irrad0
         endif
      endif


 1000 format (/,5x,'shell #',i4,', 1-sumxsp = ',1pd13.6,' > dnucc = ',
     &     1pd12.6,' with ytminc = ',1pd10.3,', t9 = ',1pd9.3,/,5x,
     &     'Try again with smaller tolerance : ytminc = ',1pd10.3,/)
 1100 format (/,5x,'shell #',i4,', 1-sumxsp = ',1pd13.6,' > dnucc = ',
     &     1pd12.6,' with ytminc = ',1pd10.3,', t9 = ',1pd9.3)


      return
      end
