

***********************************************************************

      SUBROUTINE chedif (error)

************************************************************************
* Calculate the chemical evolution due to nuclear reactions, mixing    *
* and slow transport processes                                         *
* 17/09/03 In case of rotation, nucleosynthesis is computed in a       *
* radiative way in the ENTIRE star. Homogeneisation of CZ is done      *
* afterwards in the DIFFUSION routine                                  *
*  Coupling nucleosynthesis-diffusion equations                        *
*  nmixd = 1 : core only                                               *
*  nmixd = 2 : enveloppe only                                          *
*  nmixd = 3 : core+enveloppe                                          *
*  nmixd = 4 : all convective zones                                    *
*  nmixd = 5 : entire star                                             *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.conv'
      include 'evolcom.diff'
      include 'evolcom.eos'
      include 'evolcom.ion'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer imin,imax,mmean
      integer kl,i,l,k,ij
      integer error
      integer klenv,klcore,klpulse

      double precision xi(nsp),xx(nsh,nsp)
      double precision dturb,dmic,vdmic,lgdturb,lgdmic

      logical nucl(nsh),ipmix(nmaxconv),pass,partialmix

      common /coefdiff/ dturb(nsh),dmic(nsh),vdmic(nsh),lgdturb(nsh)
     &     ,lgdmic(nsh)
      common /overshoot/ klcore,klenv,klpulse

      print*,""
      if (t(1).eq.t(2)) then
         error = 4
         return
      endif

*___________________________________________________________________
***   various abundance initializations (in case of call to ydequil)
*-------------------------------------------------------------------

      if (dmaccr.gt.0.d0.and.xspacc(ih2).gt.1.d-10.and.accphase.eq.0)
     &     then
         do k = iaccbot,iacctop
            vxsp(k,ih1) = xsp(k,ih1)
            vxsp(k,ih2) = xsp(k,ih2)
            vxsp(k,ihe3) = xsp(k,ihe3)
         enddo
      endif
      do k = 1,nmod
         nucl(k) = .false.
      enddo

*-----------------------------------------------------------------------*
***  coupling nucleosynthesis and time-dependent mixing in the all star *
*-----------------------------------------------------------------------*

c..  compute diffusion coefficients
      if (nmixd.gt.0) then
         call diffusion (1,error)
         if (error.gt.0) return
      endif
      if (nmixd.eq.5) then
         kcentd = 1
         ksurfd = 1
         write (nout,*) 'Coupling diffusion/nucleosynthesis equations ',
     &        'throughout the entire structure'
         call netdiff (1,nmod1,dturb,dtn,error)
         return
      endif


*_________________________________________
***   nucleosynthesis in convective shells
*-----------------------------------------


***   Instantaneous mixing in convective zones if NO DIFFUSION
      pass = .not.diffzc
cx      if (flame.and.diffzc) pass = .true.
      if (nsconv.gt.0.and.nmixd.lt.4.and.pass) then
c..   restore initial abundances
         do l = 1,nsp
            do k = 1,nmod
               xsp(k,l) = vxsp(k,l)
               ysp(k,l) = xsp(k,l)/anuc(l)
            enddo
         enddo
         do 10 kl = 1,nsconv
            if ((nmixd.eq.1.and.kl.eq.klcore).or.
     &           (nmixd.eq.2.and.kl.eq.klenv).or.
     &           (nmixd.eq.3.and.(kl.eq.klcore.or.kl.eq.klenv))) goto 10
            imin = novlim(kl,3)
            imax = novlim(kl,4)
            am(kl) = 0.d0
            do i = imin,imax
               am(kl) = am(kl)+dm(i)
            enddo
            am(kl) = 1.d0/am(kl)
            if (nucreg.eq.2.and..not.flame) then
               call mix (dm,imin,imax,kl,5,partialmix)
               ipmix(kl) = .true.
            else
               call mix (dm,imin,imax,kl,2,partialmix)
               ipmix(kl) = partialmix
            endif
            if (.not.ipmix(kl)) then
               if (no.eq.1) write (90,100) imin,imax
               write (nout,100) imin,imax
               mmean = int((imin+imax)/2)
               do l = 1,nsp
                  xi(l) = xsp(mmean,l)
               enddo
               call network ('c',xi,am(kl),dtn,tnucmin,imin,imax,error)
               if (error.gt.0) return
               do k = imin,imax
                  nucl(k) = .true.
                  do l = 1,nsp
                     xsp(k,l) = xi(l)
                  enddo
c..  update composition
                  if (nucreg.ne.3) then
                     do l = 1,nsp
                        vxsp(k,l) = xi(l)
                     enddo
                  endif
               enddo
            else
***   in case of partial mixing, convective zone treated as radiative
               if (no.eq.1) write (90,200) imin,imax
               write (nout,200) imin,imax
               do k = imin,imax
                  ij = k
                  do l = 1,nsp
                     xi(l) = xsp(k,l)
                  enddo
                  call network ('r',xi,0.d0,dtn,tnucmin,ij,ij,error)
                  if (error.gt.0) return
                  nucl(k) = .true.
                  if (nucreg.eq.3) then
                     do l = 1,nsp
                        xx(k,l) = vxsp(k,l)
                     enddo
                  endif
c..  update composition
                  do l = 1,nsp
                     xsp(k,l) = xi(l)
                     vxsp(k,l) = xi(l)
                  enddo
               enddo
               call mix (dm,imin,imax,kl,3,partialmix)
               if (nucreg.eq.3) then
                  do l = 1,nsp
                     do k = imin,imax
                        vxsp(k,l) = xx(k,l)
                     enddo
                  enddo
               endif
             endif
 10      continue
      endif


*-----------------------------------------------------------------------*
***  nucleosynthesis calculation with time-dependent convective mixing  *
*-----------------------------------------------------------------------*


      if (nmixd.gt.0.and.nsconv.gt.0) then
         kcentd = 1
         ksurfd = 1
         if ((nmixd.eq.1.or.nmixd.eq.3).and.klcore.eq.1) then
            write (nout,400) novlim(1,8)
            imax = novlim(1,8)
            call netdiff (1,imax,dturb,dtn,error)
            do k = 1,imax
               nucl(k) = .true.
            enddo
         endif
         if ((nmixd.eq.2.or.nmixd.eq.3).and.klenv.gt.0) then
            write (nout,500) novlim(klenv,7),nmod1
            imin = novlim(klenv,7)
            call netdiff (imin,nmod1,dturb,dtn,error)
            do k = imin,nmod
               nucl(k) = .true.
            enddo
         endif
         if (nmixd.eq.4) then
            do kl = 1,nsconv
               imin = max(1,novlim(kl,7))
               imax = min(nmod1,novlim(kl,8)+1)
               if (kl.eq.klenv) imax = nmod1
               write (nout,600) imin,imax
               call netdiff (imin,imax,dturb,dtn,error)
               do k = imin,imax
                  nucl(k) = .true.
               enddo
            enddo
         endif
      endif

*_________________________________________________________
***   nucleosynthesis in radiative and thermohaline shells
*---------------------------------------------------------

      if (no.eq.1) write (90,*) 'Processing radiative nucleosynthesis'
      write (nout,*) 'Processing radiative nucleosynthesis'
      do k = 1,nmod
         if (.not.nucl(k)) then
            ij = k
            do l = 1,nsp
               xi(l) = xsp(k,l)
            enddo
            call network ('r',xi,0.d0,dtn,tnucmin,ij,ij,error)
            if (error.gt.0) return
            nucl(k) = .true.
            do l = 1,nsp
               xsp(k,l) = xi(l)
            enddo
c..  update composition
            if (nucreg.ne.3) then
               do l = 1,nsp
                  vxsp(k,l) = xi(l)
               enddo
            endif
         endif
      enddo


 100  format (' Processing one-zone convective nucleosynthesis in [',
     &     i4,',',i4,']')
 200  format (' Processing convective zone : [',i4,',',i4,
     &     '] as a radiative one')
 400  format (' Coupling diffusion/nucleosynthesis equations in the ',
     &     'core + overshoot regions : [ 1,',i4,']')
 500  format (' Coupling diffusion/nucleosynthesis equations in the ',
     &     'envelope + overshoot regions : [',i4,',',i4,']')
 600  format (' Coupling diffusion/nucleosynthesis equations in ',
     &     'convective zone + overshoot regions : [',i4,',',i4,']')

      return
      end
