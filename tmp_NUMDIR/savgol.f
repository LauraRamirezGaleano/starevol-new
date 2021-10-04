*-----------------------------------------------------------------

      SUBROUTINE SAVGOL (cc,np,nl,nrr,ld,mn)

*-----------------------------------------------------------------
*     USES lubksb,ludcmp
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
*-----------------------------------------------------------------

      implicit none

      INTEGER ld,mn,nl,np,nrr,MMAX

      DOUBLE PRECISION cc(np)

      PARAMETER (MMAX = 6)

      INTEGER imj,ipj,j,k,kk,mm,indx(MMAX+1)

      DOUBLE PRECISION d,fac,sum,a(MMAX+1,MMAX+1),b(MMAX+1)

      if (np.lt.nl+nrr+1.or.nl.lt.0.or.nrr.lt.0.or.ld.gt.mn.or.mn.gt.
     &     MMAX.or.nl+nrr.lt.mn) stop 'bad args in savgol'
      do ipj = 0,2*mn
         sum = 0.d0
         if (ipj.eq.0) sum = 1.d0
         do  k = 1,nrr
            sum = sum+dble(k)**ipj
         enddo
         do k = 1,nl
            sum = sum+dble(-k)**ipj
         enddo
         mm = min(ipj,2*mn-ipj)
         do imj = -mm,mm,2
            a(1+(ipj+imj)/2,1+(ipj-imj)/2) = sum
         enddo
      enddo
      call ludcmp (a,mn+1,MMAX+1,indx,d)
      do j = 1,mn+1
         b(j) = 0.d0
      enddo
      b(ld+1) = 1.d0
      call lubksb (a,mn+1,MMAX+1,indx,b)
      do kk = 1,np
         cc(kk) = 0.d0
      enddo
      do k = -nl,nrr
         sum = b(1)
         fac = 1.d0
         do mm = 1,mn
            fac = fac*dble(k)
            sum = sum+b(mm+1)*fac
         enddo
         kk = mod(np-k,np)+1
         cc(kk) = sum
      enddo

      return
      end

