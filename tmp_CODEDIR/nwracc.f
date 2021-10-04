

************************************************************************

      SUBROUTINE nwracc (corm,iter,alpha1)

************************************************************************
* Accelerate the convergence of the Henyey's method                    *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.mod'

      integer iter
      integer niter,ncor,nav,nsame,nrecover,nsamemax
      integer k,l,ii

      double precision corm,alpha1
      double precision co,lnco,al
      double precision newalpha1
      double precision fac,f1,f2,facinv,ff1,f1f,ff2,f2f
      double precision f3,f4
      double precision eps,med
      double precision av,last,diflast,difmed,alp1,alp2,sum,a,b,c,d,e
      double precision factor1,factor2
      double precision step
      double precision s1,s2,s3,s4,s5
      double precision pic
      double precision p1,p2,p3,p4

      logical tochange,interpth
      logical nowait

      parameter (niter = 4, ncor = niter-1, nsamemax = 4)
c     parameter (f3 = 0.45d0, f4 = 1/dsqrt(f3))
      parameter (f3 = 0.45d0, f4 = 1.490712d0)
      parameter (fac = 3.8d0, f1 = 2.d0, f2 = 4.d0)
      parameter (facinv = 1.d0/fac, f1f = f1/fac, ff1 = fac/f1)
      parameter (f2f = f2/fac, ff2 = fac/f2)
      parameter (eps = 0.5d0)
      parameter (nowait = .true.)

      common /accelnwr/ co(neq,niter),lnco(neq,niter),
     &     al(neq,niter),nsame(neq),nrecover(neq)

      dimension corm(neq),alpha1(neq),newalpha1(neq)
      dimension tochange(neq)

      save tochange

      step (s1,s2,s3,s4,s5) = s1/(exp(s2*s5**3)+s3)+s4
      pic (p1,p2,p3,p4) = p1*exp(-p4*p4/p2)+p3

      if (iter.le.1) then
         do k = 1,neq
            nsame(k) = 1
            nrecover(k) = 0
            do l = 1,niter
               co(k,l) = 0.d0
               lnco(k,l) = 0.d0
               al(k,l) = 0.d0
            enddo
            tochange(k) = .false.
         enddo
      endif

***   Test for an interpulse period (AGB phase)
      interpth = (agbphase.and.(corm(5).lt.1.d0))

      ii = 2
      if (hydro) ii = 1
      do l = niter,2,-1
         do k = ii,neq
            co(k,l) = co(k,l-1)
            lnco(k,l) = lnco(k,l-1)
            al(k,l) = al(k,l-1)
         enddo
      enddo
      do k = ii,neq
         co(k,1) = corm(k)
         lnco(k,1) = log(abs(corm(k)))
         al(k,1) = alpha1(k)
      enddo

      if (iter.ge.ncor) then
         if (iter.ge.niter) then
            nav = niter
         else
            nav = ncor
         endif
         do k = ii,neq
            tochange(k) = .not.tochange(k)
            sum = 0.d0
            do l = 1,niter
               sum = sum+lnco(k,l)
            enddo
            med = lnco(k,1)-lnco(k,ncor)
            av = sum/nav
            last = lnco(k,1)
            diflast = av-last
            difmed = av-med
            if (difmed.le.(-abs(eps))) then
               c = (ff1-fac)/(facinv-ff1)
               a = (facinv-fac)*c
               d = fac
               b = 1.d0
            endif
            if (abs(difmed).lt.abs(eps)) then
               c = (1.d0-ff2)/(f2f-1.d0)
               a = (f2f-ff2)*c
               b = 1.d0
               d = ff2
            endif
            if (difmed.ge.abs(eps)) then
               c = (f1f-ff1)/(facinv-f1f)
               a = (facinv-ff1)*c
               b = 1.d0
               d = ff1
            endif
            e = min(diflast,3.d0)
            factor1 = step (a,b,c,d,e)
            if (diflast.le.(-abs(eps))) then
               if (difmed.gt.0.d0) then
                  a = f2f-facinv
                  b = 8.d0
                  c = facinv
                  d = 0.d0
                  factor2 = pic (a,b,c,difmed)
               endif
               if (difmed.le.0.d0) then
                  a = f2f-f1f
                  b = 8.d0
                  c = f1f
                  d = 0.d0
                  factor2 = pic (a,b,c,difmed)
               endif
            endif
            if (abs(diflast).lt.abs(eps)) then
               c = (1.d0-f1f)/(ff1-1.d0)
               a = (ff1-f1f)*c
               b = 1.d0
               d = f1f
               e = min(difmed,3.d0)
               factor2 = step (a,b,c,d,e)
            endif
            if (diflast.ge.abs(eps)) then
               if (difmed.gt.0.d0) then
                  a = ff2-ff1
                  b = 8.d0
                  c = ff1
                  d = 0.d0
                  factor2 = pic (a,b,c,difmed)
               endif
               if (difmed.le.0.d0) then
                  a = ff2-fac
                  b = 8.d0
                  c = fac
                  d = 0.d0
                  factor2 = pic (a,b,c,difmed)
               endif
            endif
            alp1 = alpha1(k)*factor1
            alp2 = alpha1(k)*factor2
            if (interpth) then
               newalpha1(k) = max(alp1,alp2)
            else
               newalpha1(k) = 0.5d0*(alp1+alp2)
            endif
            newalpha1(k) = min(newalpha1(k),alphmax)
            newalpha1(k) = max(newalpha1(k),alphmin)
** 3 consecutive decrease in the corrections
            if (lnco(k,1).lt.lnco(k,2).and.lnco(k,2).lt.lnco(k,3).and.
     &           newalpha1(k).le.alpha1(k))
     &           newalpha1(k) = min(1.2d0*alpha1(k),alphmax)
** 3 consecutive increase in the corrections
            if (lnco(k,1).gt.lnco(k,2).and.lnco(k,2).gt.lnco(k,3).and.
     &           newalpha1(k).ge.alpha1(k))
     &           newalpha1(k) = max(f3*alpha1(k),alphmin)
** Try to avoid a cycling in the luminosity
c           if (k.eq.5) then
               if (newalpha1(k).eq.alpha1(k)) then
                  nsame(k) = nsame(k)+1
               else
                  nsame(k) = 0
               endif
               if (nsame(k).eq.nsamemax) then
                  newalpha1(k) = max(f3*alpha1(k),alphmin)
                  nrecover(k) = 2
               endif
               if (nrecover(k).gt.0.and.nsame(k).ne.nsamemax) then
                  newalpha1(k) = min(f4*alpha1(k),alphmax)
                  nrecover(k) = nrecover(k)-1
               endif
c           endif
            if (nowait) then
               alpha1(k) = newalpha1(k)
            else
               if (tochange(k)) then
                  alpha1(k) = newalpha1(k)
               else
                  alpha1(k) = alpha1(k)
               endif
            endif
         enddo
      endif

      return
      end
