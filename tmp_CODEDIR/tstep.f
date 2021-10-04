

************************************************************************

      SUBROUTINE tstep (ng,dtnt,dtnf,dtnro,dtnr,dtnu,dtnl,idtnt,
     &     idtnf,idtnro,idtnr,idtnu,idtnl)

************************************************************************
* Find the smallest time-step according to the changes of              *
* the temperature and density profiles in the last two models          *
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
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.transp'
      include 'evolcom.var'

      integer idtnt,idtnf,idtnro,idtnr,idtnu,idtnl
      integer k,ng,ni,nf

      double precision vom,vrray,vxpsi
      double precision dtnt,dtnf,dtnro,dtnr,dtnu,dtnl
      double precision dtt,dtf,dtro,dtr,dtu,dtl
      double precision ft

      common /oldstar/ vom(nsh),vrray(nsh),vxpsi(nsh)

      dtnacc = 1.d99
      idtnu = 1
      dtnu = 1.d99
      idtnr = 1
      dtnr = 1.d99
      idtnf = 1
      dtnf = 1.d99
      idtnt = 1
      dtnt = 1.d99
      idtnro = 1
      dtnro = 1.d99
      idtnl = 1
      dtnl = 1.d99

      if (dtmin.eq.dtmax0) return

      ft = fts*dtn

      ni = 1
      nf = ng

***   Mass cut : only above Mcut is the time step constrained (hydro)
      if (mcut.lt.m(ng).and.hydro.and.ivisc.gt.0) then
         do k = ng,2,-1
            if (m(k).lt.mcut) goto 12
         enddo
 12      ni = k
      endif

*________________________
***   accretion time-step
*------------------------

      if (iaccr.gt.0.and.massrate.gt.0.d0) then
         if (accphase.eq.6) then
            dtnacc = ftacc*(m(nmod)-menv)*sec/(msun*massrate)
         elseif (accphase.eq.7) then
            dtnacc = ftacc*menv*sec/(msun*massrate)
         elseif (accphase.eq.8) then
            dtnacc = ftacc*m(nmod)*(1.d0-menv)*sec/(msun*massrate)
         else
            dtnacc = ftacc*(m(iacctop)-m(iaccbot))*sec/(msun*massrate)
         endif
         if (dtnacc.lt.0.d0) stop 'problem in mass accretion regim'
      endif

*________________________________________________________
***   timestep constrained in the shell region : ni to nf
*--------------------------------------------------------

      do k = ni+1,nf

         if (hydro.and.u(k).ne.vu(k)) then
            dtu = ft*abs(u(k)/(u(k)-vu(k)))
         else
            dtu = 1.d99
         endif
         dtr = ft*abs(r(k)/(r(k)-vr(k)))
         dtf = ft*abs(lnf(k)/(lnf(k)-vlnf(k)))
         dtt = ft*t(k)/abs(t(k)-vt(k))
         dtro = ft*ro(k)/abs(ro(k)-vro(k))
         dtl = ft*abs(lum(k)/(lum(k)-vlum(k)))

         if (dtu.lt.dtnu) then
            idtnu = k
            dtnu = dtu
         endif
         if (dtr.lt.dtnr) then
            idtnr = k
            dtnr = dtr
         endif
         if (dtf.lt.dtnf) then
            idtnf = k
            dtnf = dtf
         endif
         if (dtt.lt.dtnt) then
            idtnt = k
            dtnt = dtt
         endif
         if (dtro.lt.dtnro) then
            idtnro = k
            dtnro = dtro
         endif
         if (dtl.lt.dtnl) then
            idtnl = k
            dtnl = dtl
         endif
      enddo

      return
      end
