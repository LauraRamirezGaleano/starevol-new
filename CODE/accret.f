

************************************************************************

      SUBROUTINE accret

************************************************************************
* calculate the mass deposition and angular velocity in each shell due *
* to accretion                                                         *
*                                                                      *
*               PARAMETRIC ACCRETION : accphase values                 *
*                                                                      *
*  0  : accretion model, including D burning (suited for PMS phase)    *
* 1,4 : accretion inside the star (planet accretion)                   *
*  5  : uniform accretion from the surface with facc=ric (menv unknown)*
*  6  : uniform accretion from the surface M* to menv (facc unknown)   *
*  7  : uniform accretion from the surface M* to M*-menv (facc unknown)*
*  8  : uniform accretion from the surface M* to menv*M  (facc unknown)*
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
      include 'evolcom.eng'
      include 'evolcom.grad'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer itnr,itime
      integer iexit,icrasha,icenter,idet
      integer i,j,k
      integer iaccmax,itmax,ishockb,ishockt

      integer ibotC,itopC,ibotH,itopH,itopHe,ibotHe
      integer lb,lt,ideb

      double precision pts,shlim0,Lmax,tmax,Mackmax,enucmax
      double precision schwar(nsh)
      double precision dfdr
      double precision totmv,sinaccr,fdtacc0,dmaccrg,fmax,fmin0,
     &     maccr1,maccr,hx,y,fmin,rmid,dtacc0,dlim,fmean,flinea,fmx,
     &     faccmax,dminf

c     parameter (itnr = 10,pts = 1.d-2)
      parameter (itnr = 4,pts = 1.d-2)

      common /hydrodyn/ Lmax,tmax,Mackmax,enucmax,itmax,ishockb,ishockt

      external dfdr

      if (iaccr.eq.0) goto 30

      faccmax = 1.d0
      totmv = m(nmod)
      dlim = 1.d0
      lacc = 0.d0
      itime = 0
      iaccbot = 2
      iaccmax = 50
      iexit = 0
      idet = 0
      icrasha = 0
      fdtacc0 = fdtacc
      dmaccr = 0.d0
      sigma0 = xiaccr*dsqrt(g*r(nmod)*totmv)

      if (dtn*seci.lt.1.d-2.or.totmv.gt.massend.or.(tmax.gt.1.d8.and.
     &     xsp(itmax,ihe4).gt.0.5d0)) then
         iaccr = 0
         write (nout,*) 'accretion stopped'
         goto 30
      endif

***   initialization of accretion variables

      if (model.lt.0) goto 30

      dtacc0 = dtn
      dmaccr = massrate*dtn*seci
      dmaccrg = dmaccr*msun

*---------------------
***   Planet accretion
*---------------------

      if (accphase.gt.0.and.accphase.le.4) then

***   determination of the nuclear burning region boundaries

         sinthac = 1.d0
         shlim0 = shlim-1.d0
 11      shlim0 = shlim0+1.d0
         itopH = 0
         itopHe = 0
         ibotH = 0
         ibotHe = 0
         itopC = 0
         ibotC = 0
         nnulim(1:nmaxenu,1:3) = 0
         k = 1
         if (enucl(1).ge.shlim0) then
            nnulim(1,1) = 1
            nnulim(1,2) = 1
            do i = 2,nmod
               if (enucl(i).lt.shlim0) goto 101
            enddo
 101        nnulim(1,3) = i
            k = 2
         endif
         do 201 i = 2,nmod
            if (enucl(i-1).lt.shlim0.and.enucl(i).ge.shlim0) then
               nnulim(k,1) = i
               nnulim(k,2) = i
            endif
            if (nnulim(k,2).eq.0) goto 201
            if (enucl(i).gt.enucl(i-1).and.enucl(i).gt.
     &           enucl(i+1).and.enucl(i).gt.enucl(nnulim(k,2)))
     &           nnulim(k,2) = i
            if (enucl(i).lt.shlim0.and.enucl(i-1).ge.shlim0) then
               nnulim(k,3) = i
               k = k+1
               if (k.gt.nmaxenu) goto 301
            endif
 201     continue
 301     ideb = 1
         if (nnulim(1,1).eq.1) ideb = 2
         do 401 i = ideb,nmaxenu
            if (nnulim(i,1).eq.0) goto 501
            if (enucl(nnulim(i,2)).lt.shlim0) goto 401
            lb = nnulim(i,1)
            lt = nnulim(i,3)
            if (xsp(lb,ih1).lt.xsp(lt,ih1).and.xsp(lb,ihe4).gt.
     &           xsp(lt,ihe4)) then
               itopH = nnulim(i,3)
               ibotH = nnulim(i,1)
            endif
            if (xsp(lb,ihe4).lt.xsp(lt,ihe4).and.itopHe.eq.0) then
               ibotHe = nnulim(i,1)
               itopHe = nnulim(i,3)
            endif
 401     continue
 501     if (agbphase.and.k.eq.2.and.itopH.eq.0) goto 11
         do i = nsconv,1,-1
            if (novlim(i,3).ne.0) goto 601
         enddo
 601     if (i.lt.1) i = 1
         ibotC = novlim(i,3)
         itopC = novlim(i,4)
         if (nsconv.gt.1.and.mr(ibotC).gt.0.95d0) then
            ibotC = novlim(i-1,3)
            iacctop = ibotC+5
         endif
         iacctop = ibotC+5
         if (itopH+1.gt.iacctop) iacctop = itopH
         if (accphase.eq.1) then
            if (itopH+1.lt.iacctop) then
               iaccbot = itopH+1
            else
               iaccbot = ibotH
c              iaccbot = ibotC
            endif
            if (itopH.eq.0) then
               iaccbot = ibotC-1
               iacctop = iaccbot+10
            endif
         endif
         if (accphase.eq.2) iaccbot = ibotH
         if (accphase.eq.3) then
            if (itopH.ne.0) then
               iaccbot = itopHe+1
            else
               iaccbot = ibotHe
            endif
            if (iaccbot.gt.iacctop) iaccbot = ibotHe
         endif
         facc(1:nmod) = 0.d0
         macc(1:nmod) = 0.d0
c        if (itacc.eq.3) then
c           fmean = dmaccrg/(sinthac*(m(iacctop+1)-m(iaccbot)))
c           do k = iaccbot,iacctop
c              facc(k) = fmean
c            enddo
c        endif
 701     fmean = dmaccrg/(sinthac*(m(iacctop+1)-m(iaccbot)))
         flinea = 2.d0*fmean/(m(iacctop+1)-m(iaccbot))
         fmx = flinea*((m(iacctop+1)+m(iacctop))*0.5d0-m(iaccbot))
         if (fmx.gt.faccmax.and.mr(iacctop).lt.0.9d0) then
            iacctop = iacctop+10
            goto 701
         endif
         do k = iaccbot,iacctop
            facc(k) = flinea*((m(k+1)+m(k))*0.5d0-m(iaccbot))
            macc(k+1) = macc(k)+dm(k)*facc(k)*sinthac
         enddo
         icenter = 0
         iaccbot = iaccbot-1
         dmaccr = dmaccrg/msun
         raccbot = r(iaccbot)/r(nmod)
         maccbot = m(iaccbot)/m(nmod)
         if (facc(iacctop).gt.1.d0/sinthac) then
            if (iaccr.eq.1) iaccr = 2
            if (iaccr.eq.3) iaccr = 4
         endif
         goto 25
      else
         iacctop = nmod-2
      endif

      if (totmv.gt.6.d0) faccmax = 3.d0

*--------------------------------------------------
***   Novae explosion : pill-up mass on top of star
*--------------------------------------------------

      if (accphase.ge.5) then
         sinthac = 1.d0
***   accrete with f = cste
         if (accphase.eq.5) then
            fmean = abs(ric)
            dminf = m(nmod)-dmaccrg/(sinthac*fmean)
         endif
***   accrete from a fixed mass coordinate (Menv)
         if (accphase.eq.6) then
            if (menv.gt.m(nmod)) then
               write (nout,100)
               stop 'accret : accphase = 6'
            endif
            dminf = menv
            fmean = (m(nmod)-menv)/(sinthac*dmaccrg)
         endif
***   accretion with a constant mass depth
         if (accphase.eq.7) then
            if (menv.gt.m(nmod)) then
               write (nout,100)
               stop 'accret : accphase = 7'
            endif
            dminf = m(nmod)-menv
            fmean = menv/(sinthac*dmaccrg)
         endif
***   accretion with a constant relative mass depth
         if (accphase.eq.8) then
            if (menv.gt.1.d0) then
               write (nout,200)
               stop 'accret : accphase = 8'
            endif
            dminf = m(nmod)*menv
            fmean = m(nmod)*(1.d0-menv)/(sinthac*dmaccrg)
         endif
         do k = nmod,2,-1
            if (m(k).le.dminf) then
               iaccbot = k
               iacctop = nmod
               goto 1
            endif
         enddo
 1       if (accphase.eq.5) iaccbot = min(iaccbot,nmod-10)
         do k = 1,iaccbot-1
            facc(k) = 0.d0
            macc(k+1) = 0.d0
         enddo
         fmean = dmaccrg/(sinthac*(m(nmod)-m(iaccbot)))
c        flinea = 2.d0*fmean/(m(nmod)-m(iaccbot))
c        fmx = flinea*((m(nmod)+m(nmod1))*0.5d0-m(iaccbot))
         do k = iaccbot,nmod1
c           facc(k) = flinea*((m(k+1)+m(k))*0.5d0-m(iaccbot))
            facc(k) = fmean
            macc(k+1) = macc(k)+dm(k)*facc(k)*sinthac
         enddo
         macc(1) = 0.d0
         facc(nmod) = facc(nmod1)
         icenter = 0
         iaccbot = iaccbot-1
         raccbot = r(iaccbot)/r(nmod)
         maccbot = m(iaccbot)/m(nmod)
         dmaccr = dmaccrg/msun
         maccr = macc(nmod)
         time = time-dtn
         dtn = maccr*sec/(massrate*msun)
         time = time+dtn
         dmaccr = dmaccrg/msun
         if (iaccr.eq.1) iaccr = 2
         if (iaccr.eq.3) iaccr = 4
         goto 25
      endif


*-----------------------------
***   Accretion onto PMS stars
*-----------------------------

      sinaccr = dsqrt(1.d0-accrw**rir)
      thacc = 1.d0-sinaccr*sinaccr*pw13
      sinthac = sinaccr*thacc

***   Initializations

      do j = 2,nmod1
         schwar(j) = abad(j)-abla(j)
      enddo
      schwar(1) = schwar(2)
      schwar(nmod) = schwar(nmod1)
      do j = 1,itnr
         do i = 2,nmod1
            if (abs((schwar(i+1)-schwar(i))/schwar(i)).gt.pts.or.
     &           abs((schwar(i)-schwar(i-1))/schwar(i)).gt.pts)
     &           schwar(i) = abs(schwar(i+1)+schwar(i-1))*0.5d0
         enddo
      enddo

 5    if (icrasha.gt.iaccmax) then
         if (iexit.gt.10) stop 'accret'
         write (nout,300) iaccbot,dtn*seci,maccr
         write (90,300) iaccbot,dtn*seci,maccr
         iexit = iexit+1
         time = time-dtn
         if (iexit.eq.6) then
            dtn = dtacc0
            fdtacc0 = 1.d0/fdtacc
         endif
         dtn = dtn/fdtacc0
         time = time+dtn
         dmaccr = massrate*dtn*seci
         dmaccrg = dmaccr*msun
      endif

      iaccbot = 1
      fmax = dmaccrg/(totmv*sinthac)
      fmin0 = fmax/fdtacc
      facc(iaccbot) = 0.d0
      macc(iaccbot) = 0.d0
      icenter = 0
      maccr1 = 0.d0

 10   iaccbot = iaccbot+1
      facc(iaccbot) = 0.d0
      macc(iaccbot) = 0.d0
      macc(iaccbot+1) = 0.d0
      icrasha = 0

 20   facc(1) = facc(iaccbot)
      if (icenter.eq.1) then
         if (iaccr.eq.2.or.iaccr.eq.4) maccr = m(2)*sinthac*facc(1)
         if (iaccr.eq.1.or.iaccr.eq.3) maccr = m(2)*sinthac*facc(1)/
     &        (1.d0-sinthac*facc(1))
      else
         maccr = 0.d0
      endif
      macc(2) = maccr
      idet = 1
      icrasha = icrasha+1
      if (icrasha.gt.iaccmax) goto 5
      do i = iaccbot+1,iacctop
         hx = r(i)-r(i-1)
         y = hx*dfdr (i-1,0.d0,facc(i-1),schwar(i-1),idet)
         facc(i) = facc(i-1)+y/6.d0
         y = hx*dfdr (i-1,hx*0.5d0,facc(i-1)+y*0.5d0,schwar(i-1),idet)
         facc(i) = facc(i)+y/3.d0
         y = hx*dfdr (i-1,hx*0.5d0,facc(i-1)+y*0.5d0,schwar(i-1),idet)
         facc(i) = facc(i)+y/3.d0
         y = hx*dfdr (i-1,hx,facc(i-1)+y,schwar(i-1),idet)
         facc(i) = facc(i)+y/6.d0

         if (((iaccr.eq.1.or.iaccr.eq.3).and.facc(i).gt.1.d0).or.
     &        iaccbot.ge.(nmod-10)) then
            icrasha = 101
            goto 5
         endif
         if (facc(i).lt.0.d0.or.idet.eq.0) then
            goto 10
         endif
         if ((iaccr.eq.2.or.iaccr.eq.4).and.facc(i).gt.faccmax) then
            facc(i) = faccmax
         endif
         if (iaccr.eq.1.or.iaccr.eq.3) maccr = maccr+dm(i)*sinthac*
     &        facc(i)/(1.d0-sinthac*facc(i))
         if (iaccr.eq.2.or.iaccr.eq.4) maccr = maccr+dm(i)*sinthac*
     &        facc(i)
         macc(i+1) = maccr
         if (itime.eq.0.and.maccr.gt.(1.05d0*dmaccrg)) then
            if (icenter.eq.0) then
               maccr1 = maccr
               goto 10
            else
               if (icrasha.eq.1.and.facc(iaccbot).eq.fmin0) then
                  icrasha = 0
                  facc(iaccbot) = facc(iaccbot)/fdtacc
                  fmin0 = facc(iaccbot)
               else
                  fmax = facc(iaccbot)
                  facc(iaccbot) = facc(iaccbot)-(fmax-fmin)*0.5d0
               endif
               goto 20
            endif
         endif
      enddo

      if (maccr1.eq.0.d0.and.icenter.eq.0) then
         icenter = 1
         iaccbot = 2
         icrasha = 0
         fmin = facc(iaccbot)
         facc(iaccbot) = fmin0
         goto 20
      endif

      if (maccr.lt.(0.95d0*dmaccrg).and.icenter.eq.1) then
         fmin = facc(iaccbot)
         facc(iaccbot) = facc(iaccbot)+(fmax-fmin)*0.5d0
         goto 20
      endif

      time = time-dtn
      dtn = maccr*sec/(massrate*msun)
      dmaccr = maccr/msun
      time = time+dtn
      raccbot = r(iaccbot)/r(nmod)
      maccbot = m(iaccbot)/m(nmod)
      lacc = alphaccr*g*totmv*thacc/(r(nmod)*dtn)
      dlim = dtacc0/dtn
      macc(nmod) = macc(nmod1)

      if (itime.eq.0.and.icrasha.eq.1.and.icenter.eq.0.and.dlim.gt.
     &     facdt) then
         itime = 1
         iaccbot = iaccbot-2
         goto 10
      endif

***   end accretion : define accretion variables

 25   do i = nmod1,1,-1
         if (i.gt.iaccbot.and.i.lt.iacctop) then
            rmid = (r(i)-r(i-1))*0.5d0
            dfaccdr(i) = dfdr (i-1,rmid,facc(i),schwar(i),idet)
            if (i.eq.iaccbot+1.and.idet.eq.0) dfaccdr(i) = 0.d0
            tacc(i) = (macc(nmod)-macc(i))/dmaccr
            if (accphase.eq.0.and.tacc(i).gt.0.d0) then
               vacc(i) = alpha_mlt_hp*hp(i)/tacc(i)
            else
               vacc(i) = 0.d0
            endif
         else
            if (icenter.eq.0) facc(i) = 0.d0
            dfaccdr(i) = 0.d0
            tacc(i) = macc(nmod)/dmaccr
            vacc(i) = 1.d-37
         endif
      enddo
      tacc(nmod1) = tacc(nmod1-1)
      tacc(nmod) = tacc(nmod1)
      vacc(nmod1) = vacc(nmod1-1)
      vacc(nmod) = vacc(nmod1)
      dfaccdr(nmod) = dfaccdr(nmod1)
      if (accphase.eq.0.or.accphase.ge.5) then
         write (nout,400) iaccbot+1,iacctop,facc(max(iacctop-1,1))*
     &        sinthac,t(max(iaccbot-1,1)),m(max(iaccbot-1,1))/msun
      else
         write (nout,500) iaccbot+1,iacctop,facc(iacctop)*sinthac,
     &        ibotHe,itopHe,ibotH,itopH,ibotC,itopC,shlim0
      endif

      return

***   no accretion, initializations

 30   massrate = 0.d0
      dmaccr = 0.d0
      maccbot = 1.d0
      raccbot = 1.d0
      lacc = 0.d0
      facc(1:nmod) = 0.d0
      vfacc(1:nmod) = 0.d0
      macc(1:nmod) = 0.d0
      dfaccdr(1:nmod) = 0.d0
      tacc(1:nmod) = 1.d37
      vacc(1:nmod) = 1.d-37


 100  format ('** CHANGE PARAMETER CARD : menv > m(nmod)')
 200  format ('** CHANGE PARAMETER CARD : menv > 1.d0')
 300  format (5x,'non-convergence of f [',i3,']  -->  dtn = ',1pe11.5,
     &     ' yr, Macc =',1pe11.5,/)
 400  format (' accretion : ',i4,' - ',i4,' ; fmax = ',1pe8.2,
     &     ', Tbot = ',1pe8.2,', Mbot = ',0pf7.4,/)
 500  format (' accretion : ',i4,' - ',i4,' ; fmax = ',1pe8.2,
     &     ' | He shell : ',i4,' - ',i4,' | H shell : ',i4,
     &     ' - ',i4,' | envelope ',i4,' - ',i4,' | limit : ',0pf5.0,/)

      return
      end
