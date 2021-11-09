      SUBROUTINE myengen(nc,ncmax)

*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.diff'
      include 'evolcom.eng'
      include 'evolcom.eos'
      include 'evolcom.grad'
      include 'evolcom.ion'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.opa'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      double precision efn,msol,cefn,cefntot,rap,
     &         tmp,cefnpart,cxsp,cxsppart,flux,cflux,cfluxtot,
     &         cfluxpart,rapen,rapflux
      integer i,j,k,l,nc,ncmax,ktmp,kcefn,kcflux,permute,kcxsp

      common /engen/ efn(nreac),flux(nreac)
      
      parameter (msol=1.9892d33)

      dimension cefn(nreac),kcefn(nreac),cxsp(nsp),kcxsp(nsp),
     &      cflux(nreac),kcflux(nreac)
      
      
      if (nc.ne.1.and.(nc.ne.ncmax.or.ncmax.eq.1)
c     &      .and.nc.ne.75.and.nc.ne.220.and.nc.ne.35
     &   ) return 

      write (*,*)

c    determination de la composition chimique locale
c----------------------------------------------------

      write (*,444) nc,m(nc)/msol,t(nc),ro(nc),eta(nc)

      do i = 1,nsp
         cxsp(i) = xsp(nc,i)
         kcxsp(i) = i
      enddo
      
c     cxsp(i) est l'element classé i pour l'abondance
c     kcxsp(i) est le numero de l'element classé i pour l'abondance
      
      do j = 1,nsp
         permute = 0
         do i = 1,nsp-1
            if (cxsp(i+1).gt.cxsp(i)) then
               tmp = cxsp(i+1)
               cxsp(i+1) = cxsp(i)
               cxsp(i) = tmp
               ktmp = kcxsp(i+1)
               kcxsp(i+1) = kcxsp(i)
               kcxsp(i) = ktmp
               permute = 1
            endif
         enddo
      enddo

      j = 0
      cxsppart = 0.d0
      do i = 1,nsp
         cxsppart = cxsppart + cxsp(i)
         if (cxsppart.gt.0.99d0) then
            j = i
            write (*,888) nc,j,cxsppart*1.d2
            goto 30
         endif
      enddo
   30 if (j.eq.0) 
     &   write (*,*) 'Did not find the number of most abundant nuclei'
      write (*,22) (elem(kcxsp(i)),cxsp(i),i=1,j)
 22   format (4(2x,a5,' (',1pe12.6,') '))

c    determination des reactions les plus energetiques
c-----------------------------------------------------

      do i = 1,nreac
         cefn(i) = abs(efn(i))
         kcefn(i) = i
      enddo

c     cefn(i) est la reac classee i pour l'energetique
c     kcefn(i) est le numero de la reac classee i pour l'energetique
      
      do j = 1,nreac
         permute = 0
         do i = 1,nreac-1
            if (cefn(i+1).gt.cefn(i)) then
               tmp = cefn(i+1)
               cefn(i+1) = cefn(i)
               cefn(i) = tmp
               ktmp = kcefn(i+1)
               kcefn(i+1) = kcefn(i)
               kcefn(i) = ktmp
               permute = 1
            endif
         enddo
      enddo
            
      cefntot = 0.d0
      do i = 1,nreac
         cefntot = cefntot + cefn(i)
      enddo
      if (cefntot.eq.0.d0) then
c        If the following message appears, it may indicate
c        that common /engen/ has not been introduced in file nuceng.f
         write (*,*) 'Myengen: No nuclear reactions ???'
         goto 999
      endif

      cefnpart = 0.d0
      j = 0
      if (nphase.lt.5) then
         rapen = 0.99d0
      else
         rapen = 0.98d0
      endif
      do i = 1,nreac
         cefnpart = cefnpart + cefn(i)
         rap = cefnpart/cefntot
         if (rap.ge.rapen) then
            j = i
            write (*,222) nc,j,rap*1.d2
            goto 60
         endif
      enddo
      
   60 do i=1,j
         write (*,333) react(kcefn(i)),cefn(i),cefn(1)/cefn(i)
      enddo
      
      if (nphase.ge.5) then
c 2 C12 reactions 
         do k = iccga,iccgp
            l = 0
            do i=1,nreac
               if (kcefn(i).eq.k) then
                  l = i
                  goto 10
               endif
            enddo
   10       if (l.eq.0) then
               write (*,*) 'Did not find reaction # ',k,' in the list'
            else
               if (l.gt.j) then
                  write (*,555) react(k),l,cefn(l),cefn(1)/cefn(l)
               endif
            endif
         enddo

c reactions on Ne20
         do k = inega,inega
            l = 0
            do i=1,nreac
               if (kcefn(i).eq.k) then
                  l = i
                  goto 20
               endif
            enddo
   20       if (l.eq.0) then
               write (*,*) 'Did not find reaction # ',k,' in the list'
            else
               if (l.gt.j) then
                  write (*,555) react(k),l,cefn(l),cefn(1)/cefn(l)
               endif
            endif
         enddo
      endif
      

c      print *,'protons:',xsp(nc,ih1)/flux(29),xsp(nc,ih1),xsp(nc,ic13),
c     &   xsp(nc,ineutron)*avn*ro(nc)
      
c skip max fluxes:
      goto 999
      
c    determination des fluxs les plus importants
c------------------------------------------------

      do i = 1,nreac
         cflux(i) = abs(flux(i))
         kcflux(i) = i
      enddo

c     cflux(i) est la reac classee i pour l'energetique
c     kcflux(i) est le numero de la reac classee i pour l'energetique
      
      do j = 1,nreac
         permute = 0
         do i = 1,nreac-1
            if (cflux(i+1).gt.cflux(i)) then
               tmp = cflux(i+1)
               cflux(i+1) = cflux(i)
               cflux(i) = tmp
               ktmp = kcflux(i+1)
               kcflux(i+1) = kcflux(i)
               kcflux(i) = ktmp
               permute = 1
            endif
         enddo
      enddo
            
      cfluxtot = 0.d0
      do i = 1,nreac
c         if (cflux(i).gt.0.d0) cfluxtot = cfluxtot + cflux(i)
c          if (cflux(i).gt.0.d0) print *, 'flux<0'
          cfluxtot = cfluxtot + cflux(i)
      enddo

      cfluxpart = 0.d0
      j = 0
      if (nphase.lt.5) then
         rapflux = 0.99d0
      else
         rapflux = 0.98d0
      endif
      do i = 1,nreac
         cfluxpart = cfluxpart + cflux(i)
         rap = cfluxpart/cfluxtot
         if (rap.ge.rapflux) then
            j = i
            write (*,666) nc,j,rap*1.d2
            goto 90
         endif
      enddo
      
   90 if (j.eq.0) then
         write (*,*) 'myengen: No flux ???'
      else
         do i=1,j
            write (*,333) react(kcflux(i)),cflux(i),cflux(1)/cflux(i)
         enddo
      endif
      
      if (nphase.ge.5) then
c 2 C12 reactions 
         do k = 31,32
            l = 0
            do i=1,nreac
               if (kcflux(i).eq.k) then
                  l = i
                  goto 50
               endif
            enddo
   50       if (l.eq.0) then
               write (*,*) 'Did not find reaction #',k,' in the list'
            else
               if (l.gt.j) then
                  write (*,555) react(k),l,cflux(l),cflux(1)/cflux(l)
               endif
            endif
         enddo

c reactions on Ne20
c         do k = 81,85
         do k = 83,83
            l = 0
            do i=1,nreac
               if (kcflux(i).eq.k) then
                  l = i
                  goto 40
               endif
            enddo
   40       if (l.eq.0) then
               write (*,*) 'Did not find reaction #',k,' in the list'
            else
               if (l.gt.j) then
                  write (*,555) react(k),l,cflux(l),cflux(1)/cflux(l)
               endif
            endif
         enddo
      endif

  999 write (*,*)



  444 format (1x,'At shell #',i4,': m/msol = ',1pd10.4,', T = ',1pd8.2,
     &      ', rho = ',1pd8.2,', eta = ',1pd9.2)
  555 format (2x,'Reac ',a,' is at rank ',i3,' [',
     &      1pd9.2,'] R=',d9.2)
  888 format (1x,'At shell #',i4,': ',i3,
     &      ' element(s) = ',f6.2,'% of the shell abundance:')
  222 format (1x,'At shell #',i4,': ',i3,' reaction(s) provide(s) ',
     &      f6.2,'% of the nuclear energy')
  666 format (1x,'At shell #',i4,': ',
     &      i3,' reaction(s) contribute to ',f6.2,
     &      '% of the total shell flux')
  333 format (1x,a37,' [',1pd8.2,'] R=',d8.2)
      return
      end
