!=======================================================================
      subroutine a02aaf(xxr,xxi,yr,yi)
!-----------------------------------------------------------------------
!     mark 2a release.  nag copyright 1973
!     mark 4.5 revised
!     mark 5c revised
!     mark 11c revised. ier-467 (mar 1985)
!     mark 11.5(f77) revised. (sept 1985.)
!     computes the square root of a complex number
!-----------------------------------------------------------------------
      implicit none
      include 'evolcom.num'

!     .. scalar arguments ..
      real(8),intent(in):: xxi,xxr
      real(8),intent(out):: yi,yr
!     .. local scalars ..
      real(8):: h,xi,xr
      real(8):: a02abf
!-----------------------------------------------------------------------
      h=0.d0
      xr = abs(xxr)
      xi = xxi

      if (xr > one) h = sqrt(xr*half+a02abf(xr*half,xi*half))
      if (xr <= one) h = sqrt(xr+a02abf(xr,xi))*sqrt(half)
      if (xi /= zero) xi = xi/(h+h)

      if (xxr >= zero) then
         yr = h
         yi = xi
         return
      else
         if (xi >= zero) then
            yr = xi
            yi = h
            return
         else
            yr = -xi
            yi = -h
            return
         endif
      endif

      return

      end subroutine a02aaf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      real(8) function a02abf(xxr,xxi)
!-----------------------------------------------------------------------
!     nag copyright 1975
!     mark 4.5 revised
!     mark 5c revised
!     mark 11.5(f77) revised. (sept 1985.)

!     returns the absolute value of a complex number via routine
!     name
!-----------------------------------------------------------------------
      implicit none
      include 'evolcom.num'

!     .. scalar arguments ..
      real(8),intent(in):: xxi,xxr

      real(8):: h,xi,xr
!-----------------------------------------------------------------------
      xr = abs(xxr)
      xi = abs(xxi)

      if (xi > xr) then
         h = xr
         xr = xi
         xi = h
      endif
      if (xi == zero) then
         a02abf = xr
         return
      else
         h = xr*sqrt(one+(xi/xr)**2)
         a02abf = h
         return
      endif

      return

      end function a02abf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine a02acf(xxr,xxi,yyr,yyi,zr,zi)
!-----------------------------------------------------------------------
!     mark 2a release.  nag copyright 1973
!     mark 4.5 revised
!     mark 5c revised
!     mark 11.5(f77) revised. (sept 1985.)

!     divides one complex number by a second
!-----------------------------------------------------------------------
      implicit none
      include 'evolcom.num'

!     .. scalar arguments ..
      real(8),intent(in):: xxi,xxr,yyi,yyr
      real(8),intent(out):: zi,zr
      real(8):: a,h
!-----------------------------------------------------------------------
      if (abs(yyr) > abs(yyi)) then
         h = yyi/yyr
         a = one/(h*yyi+yyr)
         zr = (xxr+h*xxi)*a
         zi = (xxi-h*xxr)*a
         return
      else
         h = yyr/yyi
         a = one/(h*yyr+yyi)
         zr = (h*xxr+xxi)*a
         zi = (h*xxi-xxr)*a
         return
      endif

      return

      end subroutine a02acf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine c02agf(a,n,scale,z,work,ifail)
!-----------------------------------------------------------------------
!     mark 13 release. nag copyright 1988.
!     mark 15b revised. ier-945 (nov 1991).

!     c02agf attempts to find all the roots of the nth order real
!     polynomial equation

!     a(0)*z**n + a(1)*z**(n-1) + ... + a(n-1)*z + a(n) = 0.

!     the zeros of polynomials of degree 1 and 2 are calculated using
!     the "standard" closed formulas
!     z = -a(1)/a(0) and
!     z = (-a(1) +/- sqrt(disc))/(2*a(0)) respectively, where
!     disc = a(1)**2 - 4*a(0)*a(2).
!     for n >= 3, the roots are located iteratively using a variant of
!     laguerre's method, which is cubically convergent for isolated
!     zeros (real or complex) and linearly convergent for multiple
!     zeros.

!     c02agf itself is essentially a dummy routine whose function is to
!     partition the work array work for use by c02agz.
!     work is partitioned into 2 arrays each of size (n + 1).
!-----------------------------------------------------------------------
      implicit none
      include 'evolcom.num'

!     .. scalar arguments ..
      integer,intent(in):: n
      integer,intent(out):: ifail
      logical,intent(in):: scale
!     .. array arguments ..
      real(8),dimension(0:n),intent(in):: a
      real(8),dimension(2*(n+1)),intent(out):: work
      real(8),dimension(2,n),intent(out):: z(2,n)
!     .. local scalars ..
      real(8):: big
      integer:: i, ier, ndeg, nrec
      logical:: sc
!-----------------------------------------------------------------------
      ier = ifail
      sc = scale
      nrec = 0
      ndeg = n
      if (n < 1 .or. a(0) == zero) then
         rewind(222)
         write(222,*) 'crash in nag c02agf'
         stop
      end if
!     initialize z to be -infinity.
      big = 1.0d0/(sqrt(2.0d0)*safemin)
      do i = 1, n
         z(1,i) = -big
         z(2,i) = -big
      enddo
      call c02agz(a,ndeg,sc,z,work(1),work(n+2),ier)
      if (ier == 0) then
         ifail = 0
         return
      end if
      if (ier == 2 .or. ier == 3) then
         rewind(222)
         write(222,*) 'crash in nag c02agz'
         stop
      endif

      return
!     
      end subroutine c02agf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine c02agt(a0,b0,c0,zsm,zlg)
!-----------------------------------------------------------------------
!     mark 13 release. nag copyright 1988.
!     based on the routine  qdrtc, written by brian t. smith

!     this subroutine determines the roots of the quadratic equation

!     a0*x**2 + b0*x + c0

!     where a0, b0, and c0 are real coefficients, and zsm and zlg
!     the smallest and largest root in magnitude, respectively.

!     the roots are computed to within a relative error of a few
!     units in the last place (depending on the accuracy of the
!     basic arithmetic operations) except when underflow or overflow
!     occurs in which case the true roots are within a few units in
!     the last place of the underflow or overflow threshold.

!     if the leading coefficient is zero, the larger root is
!     set to the largest machine representable number and the
!     overflow flag ovflow is set true.  if all three coefficients are
!     zero, the overflow flag is set, both roots are set to the largest
!     representable number, but no divide check is created.

!     this program is designed to take advantage of systems that report
!     overflow and underflow conditions in an efficient way.  that is,
!     if, whenever an overflow or underflow occurs, certain flags are
!     set (that is, the logical variables ovflow and unflow in the
!     common block ac02ag), c02agt can use these indicators to indicate
!     that the roots overflow or underflow and cannot be represented.

!     however, as implemented in the nag library, the routine simply
!     assumes that the machine terminates on overflow and ignores
!     underflow.
!     
!     c02agx -- determine the exponent of a number in terms of the
!     model.
!     c02agr -- form a number with a given mantissa and exponent
!     precision.
!-----------------------------------------------------------------------
      implicit none
      include 'evolcom.num'

!     .. scalar arguments ..
      real(8),intent(in):: a0,b0,c0
!     .. array arguments ..
      real(8),dimension(2),intent(out):: zlg,zsm
!     .. local scalars ..
      real(8):: a, b, c, d, sc, sqrtd
      integer:: expbsq, sclexp, c02agx
      real(8):: c02agy
!-----------------------------------------------------------------------
!     initialize local variables with the input coefficients.
      a = a0
      b = -b0
      c = c0

!     check for  a = zero  or  c = zero.
      if (a /= zero) then
         if (c /= zero) then    ! at this point, a and c are non-zero.

!     scale the coefficients so that the product a * c is near
!     1.0d0 in magnitude.  this avoids spurious underflow/overflow
!     conditions when the true results are within range.

!     the scale factor is a power of the base near to
!     sqrt(abs(a*c)).  this choice avoids unnecessary rounding
!     errors but is expensive to compute when floating point
!     manipulative functions are not available in machine code.

            sclexp = (c02agx(a)+c02agx(c))/2

!     the scale factor is  base ** sclexp.  if a and c are scaled
!     using this scale factor as a dividend, then the
!     the scaled product a'*c' is between base**(-2) and
!     base in magnitude, where base is the base for model numbers
!     of the type of a.

!     but before performing the scaling, check to see if it is
!     necessary -- that is, if b is so large in magnitude that
!     b**2 exceeds abs(4*a*c) by more than the relative machine
!     precision for the real(8)          data type,
!     the discriminant is in effect b and no scaling is required.
!     however, if b is so small in magnitude that abs(4*a*c)
!     exceeds b**2 in magnitude by more than this same relative
!     machine precision, b is in effect zero, but a and c are
!     still scaled to avoid spurious underflows/overflows.

!     compute the exponent of the square of the scaled b.

            if (abs(b) /= zero) then
               expbsq = 2*(c02agx(b)-sclexp)
            else
               expbsq = -2*expdep
            end if

!     check if b**2 is too big.

            if (expbsq <= expdep) then ! b**2 is not too big.  scaling will be performed.

!     a and c should be scaled using the usual scale
!     manipulation function but for efficiency, the
!     scaling is performed by division.

               sclexp = min(sclexp+1,emaxm1)
               sclexp = max(sclexp,eminm1)
               sc = c02agy(one,sclexp-c02agx(one))

!     check if it is too small.

               if (expbsq < -expdep) then ! b is too small.  set it to zero.

                  b = zero
               else             ! b is neither too large nor too small.  scale it.

                  b = (b/sc)*half
               end if
               a = a/sc
               c = c/sc
               d = b*b - a*c
               sqrtd = sqrt(abs(d))
               if (d <= zero) then ! the roots are complex.

                  zlg(1) = b/a
                  zlg(2) = abs(sqrtd/a)
                  zsm(1) = zlg(1)
                  zsm(2) = -zlg(2)
               else             ! the roots are real and sqrtd is not zero.

                  b = sign(sqrtd,b) + b
                  zsm(1) = c/b
                  zsm(2) = zero
                  zlg(1) = b/a
                  zlg(2) = zero

!     because of rounding errors in the square root and
!     divisions above (particularly on machines that
!     truncate and only when b is small), the real roots may
!     be improperly ordered -- set them so that the smaller
!     one is opposite in sign to the larger one.

                  if (abs(zlg(1)) < abs(zsm(1))) then
                     zsm(1) = -zlg(1)
                     zsm(2) = -zlg(2)
                  end if
               end if
            else                ! expbsq > expdep

!     at this point, b is very large; in this case, the
!     coefficients need not be scaled as the discriminant
!     is essentially b.

               zsm(1) = c/b
               zsm(2) = zero
               zlg(1) = b/a
               zlg(2) = zero

!     because of rounding errors in the square root and
!     divisions above (particularly on machines that truncate
!     and only when b is small), the real roots may be
!     improperly ordered -- set them so that the smaller one
!     is opposite in sign to the larger one.

               if (abs(zlg(1)) < abs(zsm(1))) then
                  zsm(1) = -zlg(1)
                  zsm(2) = -zlg(2)
               end if
            end if              ! expbsq
         else                   ! c is zero, but a is not.

            zsm(1) = zero
            zsm(2) = zero
            zlg(1) = b/a
            zlg(2) = zero
         end if                 ! c
      else                      ! a is zero.  indicate that at least one root has overflowed.

         ovflow = .true.
         zlg(1) = finity
         zlg(2) = zero
         if (b == zero .and. c /= zero) then

!     a and b are zero, but c is not.  set the roots to infinity
!     but of opposite sign to indicate this.

            zsm(1) = -zlg(1)
            zsm(2) = -zlg(2)
         else
            if (b == zero) then ! all coefficients are zero.  set both roots to + infinity.

               zsm(1) = zlg(1)
               zsm(2) = zlg(2)
            else

!     a is zero, but b is not.  compute the smaller root.

               zsm(1) = c/b
               zsm(2) = zero
            end if
         end if
      end if

      return

      end subroutine c02agt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine c02agu(a0,b0,c0,zsm,zlg)
!-----------------------------------------------------------------------
!     mark 13 release. nag copyright 1988.
!     based on the routine  cqdrtc, written by brian t. smith

!     this subroutine determines the roots of the quadratic equation

!     a0*x**2 + b0*x + c0

!     where a0, b0, and c0 are complex coefficients, and zsm and zlg
!     the smallest and largest root in magnitude, respectively.

!     the roots are computed to within a relative error of a few
!     units in the last place (depending on the accuracy of the
!     basic arithmetic operations) except when underflow or overflow
!     occurs in which case the true roots are within a few units in
!     the last place of the underflow or overflow threshold.

!     if the leading coefficient is zero, the larger root is
!     set to the largest machine representable number and the
!     overflow flag ovflow is set true.  if all three coefficients are
!     zero, the overflow flag is set, both roots are set to the largest
!     representable number, but no divide check is created.

!     this program is designed to take advantage of systems that report
!     overflow and underflow conditions in an efficient way.  that is,
!     if, whenever an overflow or underflow occurs, certain flags are
!     set (that is, the logical variables ovflow and unflow in the
!     common block ac02ag), c02agu can use these indicators to indicate
!     that the roots overflow or underflow and cannot be represented.

!     however, as implemented in the nag library, the routine simply
!     assumes that the machine terminates on overflow and ignores
!     underflow.

!     c02agx -- determine the exponent of a number in terms of the
!     model.
!     c02agr -- form a number with a given mantissa and exponent
!     precision.
!-----------------------------------------------------------------------
      implicit none
      include 'evolcom.num'

!     .. array arguments ..
      real(8),dimension(2),intent(in):: a0,b0,c0
      real(8),dimension(2),intent(out):: zlg,zsm
!     .. local scalars ..
      real(8):: sc,xc1,xc2
      integer:: expbsq, sclexp, c02agx
!     .. local arrays ..
      real(8),dimension(2):: a,b,c,ct,d
!     .. statement functions ..
      real(8):: appabs,c02agy
!     .. statement function definitions ..
      appabs(xc1,xc2) = max(abs(xc1),abs(xc2))
!-----------------------------------------------------------------------
!     initialize local variables with the input coefficients.
      a(1) = a0(1)
      a(2) = a0(2)
      b(1) = -b0(1)
      b(2) = -b0(2)
      c(1) = c0(1)
      c(2) = c0(2)

!     check for  a = cmplx(zero, zero)  or  c = cmplx(zero, zero).

      if (appabs(a(1),a(2)) /= zero) then
         if (appabs(c(1),c(2)) /= zero) then

!     at this point, a and c are non-zero.

!     scale the coefficients so that the product a * c is near
!     1.0d0 in magnitude.  this avoids spurious underflow/overflow
!     conditions when the true results are within range.

!     the scale factor is a power of the base near to
!     sqrt(abs(a*c)).  this choice avoids unnecessary rounding
!     errors but is expensive to compute when floating point
!     manipulative functions are not available in machine code.

            sclexp = (c02agx(appabs(a(1),a(2)))+c02agx(appabs(c(1)
     $           ,c(2))))/2

!     the scale factor is  base ** sclexp.  if a and c are scaled
!     using this scale factor as a dividend, then the
!     the scaled product a'*c' is between base**(-2) and
!     base in magnitude, where base is the base for model numbers
!     of the type of a.

!     but before performing the scaling, check to see if it is
!     necessary -- that is, if b is so large in magnitude that
!     b**2 exceeds abs(4*a*c) by more than the relative machine
!     precision for the real(8)          data type,
!     the discriminant is in effect b and no scaling is required.
!     however, if b is so small in magnitude that abs(4*a*c)
!     exceeds b**2 in magnitude by more than this same relative
!     machine precision, b is in effect zero, but a and c are
!     still scaled to avoid spurious underflows/overflows.

!     compute the exponent of the square of the scaled b.

            if (appabs(b(1),b(2)) /= zero) then
               expbsq = 2*(c02agx(appabs(b(1),b(2)))-sclexp)
            else
               expbsq = -2*expdep
            end if

!     check if b**2 is too big.

            if (expbsq <= expdep) then

!     b**2 is not too big.  scaling will be performed.

!     a and c should be scaled using the usual scale
!     manipulation function but for efficiency, the
!     scaling is performed by division.

               sclexp = min(sclexp+1,emaxm1)
               sclexp = max(sclexp,eminm1)
               sc = c02agy(one,sclexp-c02agx(one))

!     check if it is too small.

               if (expbsq < -expdep) then

!     b is too small.  set it to zero.

                  b(1) = zero
                  b(2) = zero
               else

!     b is neither too large nor too small.  scale it.

                  b(1) = (b(1)/sc)*half
                  b(2) = (b(2)/sc)*half
               end if
               a(1) = a(1)/sc
               a(2) = a(2)/sc
               c(1) = c(1)/sc
               c(2) = c(2)/sc

!     the magnitude of the discriminant will not underflow
!     or overflow -- however, a component of it may underflow.

               ct(1) = b(1)*b(1) - b(2)*b(2) - a(1)*c(1) + a(2)*c(2)
               ct(2) = two*b(2)*b(1) - a(2)*c(1) - a(1)*c(2)
               call a02aaf(ct(1),ct(2),d(1),d(2))

!     in order to ensure that the larger root is assigned to
!     zlg, select the sign of d so that b+d is larger in
!     magnitude than b-d.  (this condition reduces to the
!     condition that real(b)*real(d)+aimag(b)*aimag(d)>0.)

               if (d(1)*b(1)+d(2)*b(2) <= zero) then
                  d(1) = -d(1)
                  d(2) = -d(2)
               end if
               b(1) = b(1) + d(1)
               b(2) = b(2) + d(2)
            end if              ! expbsq

!     at this point, b is either very large or moderate; in case
!     it is moderate, the coefficients have been scaled and b
!     represents the sum of the scaled input coefficent b0 and the
!     discriminant; in case b is very large, the coefficients need
!     not be scaled (the discriminant is essentially b), and the
!     roots can be computed with the same quotients as in the
!     former case.

            call a02acf(b(1),b(2),a(1),a(2),zlg(1),zlg(2))
            call a02acf(c(1),c(2),b(1),b(2),zsm(1),zsm(2))
         else                   ! c is zero, but a is not.

            zsm(1) = zero
            zsm(2) = zero
            call a02acf(b(1),b(2),a(1),a(2),zlg(1),zlg(2))
         end if                 ! c
      else                      ! a is zero.  indicate that at least one root has overflowed.
!     
         ovflow = .true.
         zlg(1) = finity
         zlg(2) = finity
         if (appabs(b(1),b(2)) == zero .and. appabs(c(1),c(2)) /= zero)
     $        then
!     
!     a and b are zero, but c is not.  set the roots to infinity
!     but of opposite sign to indicate this.
!     
            zsm(1) = -zlg(1)
            zsm(2) = -zlg(2)
         else
            if (appabs(b(1),b(2)) == zero) then
!     
!     all coefficients are zero.  set both roots to + infinity.
!     
               zsm(1) = zlg(1)
               zsm(2) = zlg(2)
            else
!     
!     a is zero, but b is not.  compute the smaller root.
!     
               call a02acf(c(1),c(2),b(1),b(2),zsm(1),zsm(2))
            end if
         end if
      end if                    ! a

      return

      end subroutine c02agu

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine c02agv(dx,ndeg,a,p,pprime,pdprim,error,deflat)
!-----------------------------------------------------------------------
!     mark 13 release. nag copyright 1988.
!     based on the routine  rpolyr, written by brian t. smith

!     this subroutine evaluates a polynomial of degree ndeg with real
!     coefficients at a real point, its first derivative and half its
!     second derivative at that point and an estimate of the error in
!     the polynomial value.  the estimate is a guaranteed upper bound,
!     assuming the coefficients are exact, and the floating point
!     arithmetic is 'well-behaved'.
!     that is,

!     fl(c+d) = (c+d)(1+e),  abs(e) <= deps * (small integer)
!     fl(c*d) = (c*d)(1+f),  abs(f) <= deps * (small integer)

!     where fl(.) is the machine value for the enclosed operation
!     and  deps  is a small constant, usually equal to the
!     relative machine precision.

!     the polynomial is assumed to be of the form
!     p(x) = sum(a(ndeg-i)*x**i, i=0,1,...,ndeg )
!     where  a(i) and  x  are real(8)         .
!     note: a(0) is the coefficient of x**ndeg and  a(ndeg) is the
!     constant coefficient.

!     usage: this program can be used to evaluate a real polynomial
!     at a real point and to check whether the polynomial value
!     is essentially zero.  that is, if  abs(p) <= error, the
!     point x is indistinguishable from a zero of p, because of
!     rounding error introduced by the floating point arithmetic.

!     as a second usage, the circle centered at x and of radius
!     r = ndeg*(abs(p)+error)/abs(pprime)  is known to contain a
!     zero of p.  if it is known that there is no complex zero
!     in this region, then there is a real zero in the closed
!     interval  [x-r,x+r].

!     history: this code was taken from the program  zerpol, a zero-
!     finding algorithm for polynomials, univ. of toronto,
!     computer science department, 1965.  the code has been
!     modified by b. t. smith to be portable to fortran 77
!     systems.

!     this program is designed to take advantage of systems that report
!     overflow and underflow conditions in an efficient way.  that is,
!     if, whenever an overflow or underflow occurs, certain flags are
!     set (that is, the logical variables ovflow and unflow in the
!     common block ac02ag), c02agv can use these indicators to indicate
!     that the coefficients of the polynomial need to be scaled.  such
!     scaling permits the determination of the roots of the polynomial
!     without intermediate underflow/overflow contaminating the
!     computed roots.

!     however, as implemented in the nag library, the routine simply
!     assumes that the machine terminates on overflow and ignores
!     underflow.
!-----------------------------------------------------------------------
      implicit none
      include 'evolcom.num'

!     .. scalar arguments ..
      real(8),intent(in):: dx
      real(8),intent(out):: error,p,pdprim,pprime
      integer,intent(in):: ndeg
!     .. array arguments ..
      real(8),dimension(0:ndeg),intent(in):: a
      real(8),dimension(0:ndeg),intent(out):: deflat
!     .. local scalars ..
      real(8):: absp,absx,dv,sx,w
      integer:: i
!-----------------------------------------------------------------------
      if (ndeg >= 1) then

!     assume the polynomial is scaled so that no over/underflow
!     occurs.  if the largest coefficient is larger than the
!     sqrt(finity), underflow can probably be ignored.  however,
!     overflow cannot be ignored in these calculations in general.

         sx = dx
         absx = abs(sx)
         w = zero
         dv = a(0)
         deflat(0) = dv
         deflat(1) = a(1) + dx*deflat(0)
         do i = 2,ndeg
            w= dv + sx*w
            dv = deflat(i-1) + dx*dv
            deflat(i) = a(i) + dx*deflat(i-1)
         enddo
         p = deflat(ndeg)

!     if ovflow, the coefficients must be scaled down.
!     this rescaling is to be performed by the calling program.
!     the need to rescale is indicated by the flag  ovflow.

!     for a machine that uses ieee arithmetic, an effective but
!     incomplete test for overflow is to check whether the result
!     p is equal to infinity.  this may work on other machines if
!     an overflow does not terminate the computation and is replaced
!     with a value that persists through subsequent computations.

         absp = abs(p)
         ovflow = absp >= finity
         if (.not. ovflow) then
            error = to3rds*abs(a(0))
            do i = 1,ndeg-1
               error = abs(deflat(i)) + absx*error
            enddo
            error = sxteen*deps*(abs(deflat(ndeg))+three*absx*error)

!     check to see if the error bound computation has overflowed.
!     if so, multiply it by ndeg to insure it is large and ignore
!     the overflow condition.

!     for a machine that uses ieee arithmetic, an effective but
!     incomplete test for overflow is to check whether the result
!     error is equal to infinity.  this may work on other machines
!     if an overflow does not terminate the computation and is
!     replaced with a value that persists through subsequent
!     computations.

            ovflow = error >= finity
            if (ovflow) then
               p = finity
               pprime = finity
               pdprim = finity
               error = tiny
               ovflow = .false.
            else

!     compute the first and second derivatives.

               pprime = dv
               pdprim = w
            end if
         end if
      else if (ndeg == 0) then
         p = a(0)
         error = zero
         pprime = zero
         pdprim = zero
      end if

      return

      end subroutine c02agv

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine c02agw(xr,xi,ndeg,a,p,pprime,pdprim,error,deflat)
!-----------------------------------------------------------------------
!     mark 13 release. nag copyright 1988.
!     based on the routine  cpolyr, written by brian t. smith

!     this subroutine evaluates a polynomial of degree ndeg with real
!     coefficients at a complex point (xr, xi), its first derivative
!     and half its second derivative at that point and an estimate of
!     the error in the polynomial value.  the estimate is a guaranteed
!     upper bound, assuming the coefficients are exact, and the floating
!     point arithmetic is 'well-behaved'.
!     that is,

!     fl(c+d) = (c+d)(1+e),  abs(e) <= deps * (small integer)
!     fl(c*d) = (c*d)(1+f),  abs(f) <= deps * (small integer)

!     where fl(.) is the machine value for the enclosed operation
!     and  deps  is a small constant, usually equal to the
!     relative machine precision.

!     the polynomial is assumed to be of the form
!     p(x) = sum(a(ndeg-i)*x**i, i=0,1,...,ndeg )
!     where  a(i)  and  x  (the indeterminate) are of type complex.
!     note: a(0) is the coefficient of x**ndeg and  a(ndeg) is the
!     constant coefficient.

!     usage: this program can be used to evaluate a real polynomial
!     at a complex point and to check whether the polynomial
!     value is essentially zero.  that is, if  abs(p) <= error,
!     the point x is indistinguishable from a zero of p, because
!     of rounding error introduced by the floating point
!     arithmetic.

!     as a second usage, the circle centered at x and of radius
!     r = ndeg*(abs(p)+error)/abs(pprime)  is known to contain a
!     zero of p.  if it is known that there is no complex zero
!     in this region, then there is a real zero in the closed
!     interval  [x-r,x+r].

!     history: this code was taken from the program  zerpol, a zero-
!     finding algorithm for polynomials, univ. of toronto,
!     computer science department, 1965.  the code has been
!     modified by b. t. smith to be portable to fortran 77
!     systems.

!     this program is designed to take advantage of systems that report
!     overflow and underflow conditions in an efficient way.  that is,
!     if, whenever an overflow or underflow occurs, certain flags are
!     set (that is, the logical variables ovflow and unflow in the
!     common block ac02ag), c02agw can use these indicators to indicate
!     that the coefficients of the polynomial need to be scaled. such
!     scaling permits the determination of the roots of the polynomial
!     without intermediate underflow/overflow contaminating the
!     computed roots.

!     however, as implemented in the nag library, the routine simply
!     assumes that the machine terminates on overflow and ignores
!     underflow.
!-----------------------------------------------------------------------
      implicit none
      include 'evolcom.num'

!     .. scalar arguments ..
      integer,intent(in):: ndeg
      real(8),intent(in):: xi,xr
      real(8),intent(out):: error
!     .. array arguments ..
      real(8),dimension(0:ndeg),intent(in):: a
      real(8),dimension(0:ndeg),intent(out):: deflat
      real(8),dimension(2),intent(out):: p,pdprim,pprime
!     .. local scalars ..
      real(8):: absp,absx,dr,dsc,dt,dt1,dv,dx,dx2,dy,s,s1,sc,v,xc1,xc2
      integer:: i
!     .. local arrays ..
      real(8),dimension(2):: cx
!     .. statement functions ..
      real(8):: apxabs, a02abf
!     .. statement function definitions ..
      apxabs(xc1,xc2) = abs(xc1) + abs(xc2)
!-----------------------------------------------------------------------
      if (ndeg >= 1) then

!     assume the polynomial is scaled so that no over/underflow
!     occurs.  if the largest coefficient is larger than the
!     sqrt(finity), underflow can probably be ignored.  however,
!     overflow cannot be ignored in these calculations in general.

         s = zero
         s1 = zero
         dt1 = zero
         dt = a(0)
         deflat(0) = a(0)
         dx = xr
         dy = xi

!     sc  is computed to check if scaling is needed to avoid
!     spurious overflow or underflow.  the scale factor needs to be
!     an upper bound for the modulus of the root to avoid such range
!     exceptions.  however, the exact modulus of the root is needed
!     to compute an estimate of the error in the polynomial value.

         cx(1) = dx
         cx(2) = dy
         absx = a02abf(cx(1),cx(2))
         sc = apxabs(cx(1),cx(2))

!     if  abs(cmplx(x,y)).le.sqrt(smallest no.), scale up   x and y.
!     if  abs(cmplx(x,y)).ge.sqrt(largest  no.), scale down x and y.

         if (sc < sqrtfy .and. sc >= sqrtty) then

!     scaling of  dx2  and  dr  is unnecessary.

            dx2 = dx + dx
            dr = dx*dx + dy*dy
            if (ndeg >= 3) then
               deflat(1) = a(1) + dx2*a(0)
               deflat(2) = a(2) + (dx2*deflat(1)-dr*a(0))
               do i = 3, ndeg - 1
                  v = s1*dr
                  s1 = s
                  s = dt1 + (dx2*s-v)
                  dv = dt1*dr
                  dt1 = dt
                  dt = (dx2*dt-dv) + deflat(i-2)
                  deflat(i) = a(i) + (dx2*deflat(i-1)-dr*deflat(i-2))
               enddo
               v = s1*dr
               s1 = s
               s = dt1 + (dx*s-v)
               dv = dt1*dr
               dt1 = dt
               dt = (dx*dt-dv) + deflat(ndeg-2)
               deflat(ndeg) = a(ndeg) + (dx*deflat(ndeg-1)-dr
     $              *deflat(ndeg-2))
            else if (ndeg == 2) then
               deflat(1) = a(1) + dx2*a(0)
               deflat(2) = a(2) + (dx*deflat(1)-dr*a(0))
            else

!     ndeg is 1.

               dt = zero
               deflat(1) = a(1) + dx*a(0)
            end if
         else                   ! sc >= sqrtfy .or. sc < sqrtty

!     scale  dx  and  dy  lest  dr  overflows or underflows.

            dsc = sc
            dx = dx/dsc
            dy = dy/dsc

!     dr  cannot overflow, fortunately.

            dr = (dx*dx+dy*dy)*dsc
            dx2 = dx + dx
            if (ndeg >= 3) then
               deflat(1) = a(1) + (dx2*a(0))*dsc
               deflat(2) = a(2) + (dx2*deflat(1)-dr*a(0))*dsc
               do i = 3, ndeg - 1
                  v = s1*dr
                  s1 = s
                  s = dt1 + (dx2*s-v)*dsc
                  dv = dt1*dr
                  dt1 = dt
                  dv = dx2*dt - dv
                  dt = deflat(i-2) + dv*dsc
                  deflat(i) = a(i) + (dx2*deflat(i-1)-dr*deflat(i-2))
     $                 *dsc
               enddo
               v = s1*dr
               s1 = s
               s = dt1 + (dx*s-v)*dsc
               dv = dt1*dr
               dt1 = dt
               dv = dx*dt - dv
               dt = deflat(ndeg-2) + dv*dsc
               deflat(ndeg) = a(ndeg) + (dx*deflat(ndeg-1)-dr
     $              *deflat(ndeg-2))*dsc
            else if (ndeg == 2) then
               deflat(1) = a(1) + (dx2*a(0))*dsc
               deflat(2) = a(2) + (dx*deflat(1)-dr*a(0))*dsc
            else

!     ndeg is 1.

               dt = zero
               deflat(1) = a(1) + (dx*a(0))*dsc
            end if
            dy = dy*dsc
         end if

!     compute the polynomial value.

         p(1) = deflat(ndeg)
         p(2) = dy*deflat(ndeg-1)

!     if ovflow, the coefficients must be scaled down.
!     this rescaling is to be performed by the calling program.
!     the need to rescale is indicated by the flag  ovflow.

!     for a machine that uses ieee arithmetic, an effective but
!     incomplete test for overflow is to check whether the result
!     p is equal to infinity.  this may work on other machines if
!     an overflow does not terminate the computation and is replaced
!     with a value that persists through subsequent computations.

         absp = a02abf(p(1),p(2))
         ovflow = absp >= finity
         if (.not. ovflow) then
            error = abs(a(0))
            do i = 2, ndeg - 1
               error = abs(deflat(i-1)) + absx*error
            enddo
            error = sxteen*deps*((nine*error*absx+three*abs(deflat(ndeg
     $           -1)))*absx+abs(deflat(ndeg)))

!     check to see if the error bound computation has overflowed.
!     if so, multiply it by ndeg to insure it is large and ignore
!     the overflow condition.

!     for a machine that uses ieee arithmetic, an effective but
!     incomplete test for overflow is to check whether the result
!     error is equal to infinity.  this may work on other machines
!     if an overflow does not terminate the computation and is
!     replaced with a value that persists through subsequent
!     computations.

            ovflow = error >= finity
            if (ovflow) then
               error = ndeg*error
               ovflow = .false.
            end if

!     compute the first and second derivatives.  if overflow
!     occurs in this computation, the flag ovflow will be set.

            dv = two*dy
            pprime(1) = deflat(ndeg-1) - dv*dt1*dy
            pprime(2) = dv*dt
            pdprim(1) = dt - dv*(dv*s)
            pdprim(2) = dy*(three*dt1-dv*(dv*s1))
         end if
      else if (ndeg == 0) then
         error = zero
         p(1) = a(0)
         p(2) = zero
         pprime(1) = zero
         pprime(2) = zero
         pdprim(1) = zero
         pdprim(2) = zero
      end if

      return

      end subroutine c02agw

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      integer function c02agx(dx)
!-----------------------------------------------------------------------
!     mark 13 release. nag copyright 1988.
!     mark 15 revised. ier-892 (apr 1991).
!     mark 16 revised. ier-969 (jun 1993).
!     mark 18 revised (thread safety). (sep 1996).

!     based on the routine  dexpnt, written by brian t. smith

!     this function computes the exponent e where x is represented
!     as  0.0 or s * f * b**e  where  s  is a sign (plus or minus one),
!     f  is a fraction, either zero or satisfies  1/dbase <= f < 1,
!     and  e  satisfies  -1021 <= e <= 1024.
!-----------------------------------------------------------------------
      implicit none
      include 'evolcom.num'

!     .. scalar arguments ..
      real(8),intent(in):: dx
!     .. local scalars ..
      real(8):: a,absx,depsc,dpnewu,temp, c02agy
      integer:: e,mnexp,mxexp
!-----------------------------------------------------------------------
      dbase = 2
      depsc = realprecision
      mnexp = -1021
      mxexp = 1024
      dpnewl = safemin
      dpnewu = safemax
      newu = mxexp - 1
      newl = mnexp - 1
      temp = dble(dbase)*(one-depsc)
      fact = dpnewu/temp
      if (dx /= zero) then      ! dx is in the range

!     dbase**(mnexp-1) <= abs(dx) < dbase**mxexp.

         absx = abs(dx)
         e = int(log(absx)/log(dble(dbase)))
         if (e >= mxexp) then
            e = mxexp - 1
         else if (e < mnexp) then
            e = mnexp
         end if
         a = absx/c02agy(one,e)
         do
            if (a >= one) then
               e = e + 1
               a = a/dbase
            else if (a < one/dbase) then
               e = e - 1
               a = a*dbase
            else
               exit
            end if
         enddo
      else
         e = 0
      end if
      c02agx = e

      return

      end function c02agx

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      real(8) function c02agy(dx,exp)
!-----------------------------------------------------------------------
!     mark 13 release. nag copyright 1988.
!     mark 15 revised. ier-893 (apr 1991).
!     based on the routine  dscale, written by brian t. smith

!     this function computes the scaled product  dx * b ** exp
!     where  b  is the base of entities of type  real(8)         .
!-----------------------------------------------------------------------
      implicit none
      include 'evolcom.num'

!     .. scalar arguments ..
      real(8),intent(in):: dx
      integer,intent(in):: exp
!     .. local scalars ..
      real(8):: dpe,dsc,power
      integer:: e
!-----------------------------------------------------------------------
      e = exp
      dsc = dx

!     if the exponent scaling is out of range for the precomputed
!     powers of the base, scale repetitively by the largest
!     precomputed power until e is within range.
!     check for e too large.
      do while (e > newu)
         dsc = dsc*fact
         e = e - newu
      enddo

!     check for e too small.
      do while (e < newl)
         dsc = dsc*dpnewl
         e = e - newl
      enddo

!     scale by the remaining scaling factor.
!     set dpe = dbase**e.

      if (e == 0) then
         dpe = one
      else
         if (e < 0) then
            e = -e
            power = one/dbase
         else
            power = dbase
         end if
         dpe = one
         do
            if (mod(e,2) == 1) then
               dpe = dpe*power
            endif
            e= e/2
            if (e <= 0) then
               exit
            end if
            power = power*power
         enddo
      end if
      c02agy = dsc*dpe

      return

      end function c02agy

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine c02agz(a,ndeg,scale,z,du,deflat,ier)
!-----------------------------------------------------------------------
!     mark 13 release. nag copyright 1988.
!     mark 14 revised. ier-709 (dec 1989).
!     mark 16 revised. ier-970 (jun 1993).
!     mark 17 revised. ier-1630 (jun 1995).
!     based on the routine  zerpol, written by brian t. smith

!     this subroutine computes the n zeros of the real polynomial
!     a(0)*z**n + a(1)*z**(n-1) + ... + a(n-1)*z + a(n) = 0
!     given by cmplx(z(1,j),z(2,j)), where j = 1,2,...,n.

!     gama, theta and phi are arbitrary parameters which for this
!     implementation have been set to 1.0, 2.0 and 0.2 respectively.
!     there is no inherent limitation on the degree other than as the
!     degree of a polynomial increases, its roots become ill-conditioned
!     as a function of the coefficients.

!     this program is designed to take advantage of systems that report
!     overflow and underflow conditions in an efficient way.  that is,
!     if, whenever an overflow or underflow occurs, certain flags are
!     set (that is, the logical variables ovflow and unflow in the
!     common block ac02ag), c02agz can use these indicators to optimally
!     scale the coefficients of the polynomial.  the optimal scaling
!     permits the determination of the roots of the polynomial without
!     intermediate underflow/overflow contaminating the computed roots.

!     however, as implemented in the nag library, the routine simply
!     assumes that the machine terminates on overflow and ignores
!     underflow.

!     53 -- number of digits in the mantissa of model numbers of
!     type real(8)         .
!     x02ajf -- relative machine precision for entities of type
!     real(8)         .
!     1024 -- maximum exponent of entities of type real(8)         .
!     c02agy -- scale the first argument by a value with an exponent
!     equal to the second argument.
!     c02agx -- determine the exponent of a number in terms of the
!     model.
!     c02agw -- compute the polynomial value, first derivative, and
!     second derivative at a complex point.
!     c02agu -- determine the roots of a quadratic equation with complex
!     coefficients.
!     c02agt -- determine the roots of a quadratic equation with real
!     coefficients.
!     c02agv -- compute the polynomial value, first derivative, and
!     second derivative at a real point.
!     c02ags -- returns true if the argument is unnormalized or zero.
!-----------------------------------------------------------------------
      implicit none
      include 'evolcom.num'

!     .. parameters ..
      real(8),parameter:: gama=1.0d0,theta=2.0d0,phi=0.2d0,small=1.0d-3
     $     ,bigone=1.0001d0,smlone=0.99999d0,rconst=1.445d0,onepqt
     $     =1.25d0
!     .. scalar arguments ..
      integer,intent(in):: ndeg
      integer,intent(inout):: ier
      logical,intent(out):: scale
!     .. array arguments ..
      real(8),dimension(0:ndeg),intent(in):: a
      real(8),dimension(0:ndeg),intent(out):: deflat,du
      real(8),dimension(2,ndeg),intent(out):: z
!     .. local scalars ..
      real(8):: abdir,abdiro,abscl,dx,dz0i,dz0r,dzni,dznr,e,f,f0,fejer
     $     ,fn,g,lowerb,mxcoef,r,ratio,rtn,s,t,upperb,v,w, 
     &     x2n,x2n1,xn,xn1,xn2,xn2n
      integer:: i,iers,ihalf,ispir,iter,k,mxcfex,n,scbyex, c02agx
      logical:: cauchy,contin,ovf,savo,savu,spiral,startd,unf
!     .. local arrays ..
      real(8),dimension(2):: c,cdir,cdiro,cf,cf1,cf2,cl,cspir,ctemp
      real(8):: c02agy, a02abf, f06blf
!-----------------------------------------------------------------------
      tiny = safemin
      sqrtty = sqrt(tiny)
      finity = safemax
      sqrtfy = sqrt(finity)
      expdep = 54
      eminm1 = -1022
      emaxm1 = 1023
      lrgexp = emaxm1 + 1 - expdep
      deps = realprecision
      iers = ier
      ier = 0
      iter = 0
      ihalf = 0
      ispir = 0
      n = ndeg
      abdiro = 0.d0
      abscl = 0.d0
      fejer = 0.d0
      lowerb = 0.d0
      upperb = 0.d0
      rtn = 0.d0
      x2n = 0.d0
      x2n1 = 0.d0
      xn = 0.d0
      xn1 = 0.d0
      xn2n = 0.d0
      cspir(:)=0.d0

!     save overflow/underflow indicators of the calling program.
      savo = ovflow
      savu = unflow
      ovf = .false.
      unf = .false.

!     move the coefficients a(i) to du(i) and determine the largest
!     coefficient in magnitude.
      mxcoef = zero
      do i = 0,n
         du(i) = a(i)
         mxcoef = max(mxcoef,abs(a(i)))
      enddo
      if (mxcoef == zero) then
         do i = 1,n
            z(1,i) = finity
            z(2,i) = zero
         enddo
         n = 0
         ovf = .true.
      else

!     determine a scaling for the coefficients so that the largest
!     coefficient in magnitude is large in magnitude -- that is, near
!     deps * base ** 1024, unless the largest coefficient in
!     magnitude is larger than this quantity.  in this case, set
!     scale to .false. and do not scale the coefficients.

         mxcfex = c02agx(mxcoef)
         if (mxcfex > lrgexp) then
            scbyex = 0
            scale = .false.
         else
            scbyex = lrgexp - mxcfex
         end if
      end if

!     indicate that the cauchy region containing the smallest zeros
!     of the current polynomial has not been computed.

      cauchy = .false.

!     do loop while n>2.

      l0: do while (n > 2)
!     if scale = .true., scale the coefficients so that the largest
!     coefficient in magnitude is large.

      if (scale) then
         if (scbyex /= 0) then
            do i = 0, n
               du(i) = c02agy(du(i),scbyex)
            enddo
            scbyex = 0
         end if
      end if
      unf = unflow .or. unf
      
!     find the number i of consecutive leading coefficients equal to
!     zero.

      do i = 0,n-1
         if ((du(i)+zero) == zero) then

!     each vanished leading coefficient yields an infinite
!     zero.

            z(1,n-i) = finity
            z(2,n-i) = zero
         else

!     exit the loop on i -- the first non-zero coefficient is
!     the i-th coefficient.

            exit
         end if
      enddo
      if (i /= 0) then

!     slide back the coefficients and declare overflow.

         do k = i,n
            du(k-i) = du(k)
         enddo
         n = n - i

!     give an error message if the vanishing leading coefficients
!     have occurred because the polynomial has been scaled down.
!     this can only happen if overflow was detected during
!     the computation with the original input coefficients.

         if (scbyex == -expdep) then
            ier = 3
            if (iers /= 1) then
               write (*,'("** c02agf p(z) evaluation overflow.")')
               write (*,'(" ** if this message occurs contact nag.")')
            end if
            return
         end if
         
!     signal overflow and cycle the do loop on n.
         
         ovf = .true.
         cycle l0
      end if
      
!     find the number i of consecutive trailing coefficients equal to
!     zero.

      do i = n,1,-1
         if ((du(i)+zero) == zero) then

!     extract roots (if any) from the origin = (0., 0.)

            z(1,i) = zero
            z(2,i) = zero
         else

!     exit the loop on i -- the first non-zero coefficient is
!     the i-th coefficient.

            exit
         end if
      enddo
      if (i /= n) then

!     reduce the degree by the number of zero roots and
!     then cycle the do loop on n.

         n = i
         cycle l0
      end if

!     initialize logical underflow/overflow condition status
!     variables.

      ovflow = .false.
      unflow = .false.

!     henceforth  n .gt. 2,  du0 .ne. 0. , and  du(n) .ne. 0.
!     check to see whether the cauchy bounds need to be computed.

      if (.not. cauchy) then

!     initialize some useful constants.

         xn = dble(n)
         xn1 = dble(n-1)
         xn2 = dble(n-2)
         x2n = two/xn
         x2n1 = x2n/xn1
         xn2n = xn2/xn
         rtn = sqrt(xn)

!     calculate  g, an upper bound for the smallest zero.
!     start with  g = abs( geometric mean of the zeros).

         g = exp((log(abs(du(n)))-log(abs(du(0))))/xn+small)

!     calculate laguerre-step  cdir  and  fejer-bound, which is
!     an upper bound for the smallest zero.
!     calculation of the laguerre step involves the square of
!     reciprocal of newton's step.  since it can easily overflow,
!     the fejer bound is calculated with no such overflows and the
!     laguerre step is calculated from it.

         ovflow = .false.
         r = f06blf(du(n-1),du(n),ovflow)

!     if ovflow, a root of polynomial is within
!     n * base ** (-1022) of  zero.

         if (ovflow) then

!     thus, assume a root is zero, by assuming du(n) is zero.

            z(1,n) = zero
            z(2,n) = zero
            n = n - 1
            
!     cycle the do loop on n.
            
            cycle l0
         end if

!     the laguerre step and fejer bounds are computed from the
!     smaller root of a quadratic polynomial.

         call c02agt(x2n1*du(n-2),x2n*du(n-1),du(n),c,cf1)
         r= xn2n*r
         ctemp(1) = c(1)*r + xn1
         ctemp(2) = c(2)*r
         call a02acf(c(1),c(2),ctemp(1),ctemp(2),cdiro(1),cdiro(2))
         abdiro = a02abf(cdiro(1),cdiro(2))
         g= min(g,bigone*min(a02abf(c(1),c(2)),rtn*abdiro))
         
!     calculate the cauchy-lower bound  r  for the smallest zero
!     by solving for the root  r  of the polynomial equation
!     abs(du(n)) = sum( abs(du(i))*r**(n-i), i = 0, n-1 )
!     using  newton's method.

         r= g
         s= bigone*g
         unflow = .false.

!     newton iteration loop for the cauchy lower bound r.

         do while (r < s)
            t = abs(du(0))
            s = zero
            ovflow = .false.
            do i = 1, n - 1
               s = r*s + t
               t = r*t + abs(du(i))
            enddo
            s = r*s + t
            
!     it can be proved that s cannot underflow.
            
            t = (r*t-abs(du(n)))/s
            s = r
            r = r - t
         enddo

         if (ovflow) then

!     the coefficients are too large;  scale them down and
!     then cycle the do loop on n.

            scbyex = -expdep
            cycle l0
         end if

!     abs( smallest root ) < r/(2**(1/n) - 1 ) <  1.445*n*r.
!     thus, 1.445*n*r is another upper bound and the cauchy bound
!     has been computed, so set

         cauchy = .true.
         upperb = min(rconst*xn*r,g)
         lowerb = smlone*s
         unf = unflow .or. unf
      end if
      
!     now   lowerb < abs( smallest zero ) < upperb
!     initialize the iteration to begin at the origin.
!     (in the code below, f0 is initialized but its value never
!     usefully referenced -- it avoids reference to an uninitialized
!     variable in the test to accept the next iterate when the
!     iteration is not started.)

      fejer = upperb
      g = upperb
      cdir(1) = cdiro(1)
      cdir(2) = cdiro(2)
      abdir = abdiro
      ratio = abdir/g
      dznr = zero
      dzni = zero
      fn = abs(du(n))
      f0 = fn
      spiral = .false.
      startd = .false.
      contin = .true.

!     do while (contin) loop, searching for a real root,
!     or pair of complex roots.  the next iterate is
!     zn=cmplx(dznr , dzni).

      l1:do while (contin)
      iter = iter + 1

!     re-entry point to accept, modify, or reject the
!     laguerre step.

!     reject  cdir  if  abs(cdir) > theta*g .

      if (ratio > theta) then

!     current laguerre step is not acceptable.
!     if startd, reduce previous laguerre step by half.

         if (startd) then
            ihalf = ihalf + 1
            abscl = half*abscl
            cl(1) = half*cl(1)
            cl(2) = half*cl(2)

!     has the step become negligible ?

            dx = abs(dznr) + abs(dzni)
            if (dx+abscl /= dx) then
               dznr = dz0r + cl(1)
               dzni = dz0i + cl(2)
            else

!     otherwise, c02agf has hung-up.

               if (fn >= e*xn**2) then
                  ier = 2
                  if (iers /= 1) then
                     write (*,'(" ** the method has failed.")')
                   write (*,'("this error is very unlikely to occur.")')
                  end if
                  return        !   <== exit
               end if

!     exit the iteration loop  do while(contin).

               contin = .false.
               cycle l1
            end if
         else

!     if .not. startd, has zn been on the inner cauchy
!     radius ?

            ispir = ispir + 1
            if (spiral) then
               c(1) = cspir(1)*dznr - cspir(2)*dzni
               c(2) = cspir(2)*dznr + cspir(1)*dzni
            else

!     set spiral to  .true..  put  zn  on the inner
!     circle of the annulus containing the smallest
!     zero in the direction of the laguerre step.

               spiral = .true.
               cspir(1) = -onepqt/xn
               cspir(2) = one
               abscl = lowerb/xn**2
               ctemp(1) = cdir(1)/abdir
               ctemp(2) = cdir(2)/abdir
               c(1) = ctemp(1)*lowerb
               c(2) = ctemp(2)*lowerb
            end if

!     set  zn  to the next point on the spiral.

            dznr = c(1)
            dzni = c(2)
         end if
      else

!     cdir  at the origin is in the direction of decreasing
!     function value, so

         startd = .true.

!     accept  cdir  if  abs(cdir) <= gama*g.

         if (ratio > gama .and. (startd .or. spiral .or. lowerb <= gama
     $        *g)) then
            ratio = gama/ratio
            cdir(1) = cdir(1)*ratio
            cdir(2) = cdir(2)*ratio
            abdir = abdir*ratio
         end if

!     accept the previous iterate.
!     save the data associated with the current iterate.

         g = fejer
         cl(1) = cdir(1)
         cl(2) = cdir(2)
         abscl = abdir
         f0 = fn
         dz0r = dznr
         dz0i = dzni
         dznr = dz0r + cl(1)
         dzni = dz0i + cl(2)
      end if

!     is  zn  close to the real axis relative to step size ?

      if (abs(dzni) > phi*abscl) then

!     zn  is complex.  thus, attempt to divide the polynomial
!     by the quadratic factor  (z**2-x2*z+r); that is,

!     sum(du(i)*z**(n-i)) =
!     (z**2-x2*z+r) * sum( d(i)*z**(n-i-2), i=0,n-2) +
!     d(n-1)*(z-x) + d(n)   for all z,

!     where (x,y) and (x,-y) are zeros of  z**2-x2*z+r.
!     in the code below,
!     e  is error bound for the value of polynomial,
!     cf is the  value of the polynomial at (x,y),
!     cf1 is the first  derivitive of polynomial at (x,y),
!     2*cf2 is the second derivative of polynomial at (x,y),
!     and d(i) are the coefficients of quotient polynomial
!     deflat(i).

!     be sure that the overflow indicator is turned off.

         unflow = .false.
         ovflow = .false.
         call c02agw(dznr,dzni,n,du,cf,cf1,cf2,e,deflat)
         fn = a02abf(cf(1),cf(2))

!     check for overflow.

         if (ovflow) then

!     indicate that the polynomial needs to be scaled down
!     and cycle the do loop on n.  note: cauchy is not reset
!     as the cauchy bounds need not be recomputed.

            scbyex = -expdep
            cycle l0
         end if

!     check to see if  zn  is a zero or if underflow has
!     occurred.

         if (fn <= e .or. unflow) then
            if (unflow) then
               ier = 3
               if (iers /= 1) then
                  write (*,'(" ** c02agf p(z) evaluation underflow.")')
                  write (*,'("**if this message occurs contact nag.")')
               end if
               return           !   <== exit
            end if

!     a root has been found -- exit the iteration
!     loop  do while(contin).

            contin = .false.
            cycle l1
         end if
         unf = unflow .or. unf

!     has the function value decreased ?

         if (fn >= f0 .and. startd) then

!     no, it has not.  indicate that the laguerre step is
!     unacceptable.  (a ratio larger than theta indicates
!     that the laguerre step should be shortened.)

            ratio = bigone*theta

!     cycle iteration loop  do while(contin).

            cycle l1
         end if

!     find the laguerre step at  zn.

         ovflow = .false.
         call a02acf(cf1(1),cf1(2),cf(1),cf(2),c(1),c(2))

!     if ovflow, a root of polynomial is within  n*2**(-1022)
!     of  zero.

         if (ovflow) then
            unf = .true.

!     a root has been found -- exit the iteration loop.

            contin = .false.
            cycle l1
         end if

!     compute the laguerre step  cdir  and the bound  fejer
!     at  zn.  the laguerre step and fejer bounds are computed
!     from the smaller root of a quadratic polynomial.

         cf2(1) = cf2(1)*x2n1
         cf2(2) = cf2(2)*x2n1
         ctemp(1) = cf1(1)*x2n
         ctemp(2) = cf1(2)*x2n
         call c02agu(cf2,ctemp,cf,cdir,cf1)
         fejer = a02abf(cdir(1),cdir(2))
         ctemp(1) = c(1)*xn2n
         ctemp(2) = c(2)*xn2n
         c(1) = ctemp(1)*cdir(1) - ctemp(2)*cdir(2)
         c(2) = ctemp(2)*cdir(1) + ctemp(1)*cdir(2)
         c(1) = c(1) + xn1
         ctemp(1) = cdir(1)
         ctemp(2) = cdir(2)
         call a02acf(ctemp(1),ctemp(2),c(1),c(2),cdir(1),cdir(2))
         abdir = a02abf(cdir(1),cdir(2))
         ratio = abdir/g
         fejer = min(rtn*abdir,fejer)

!     is the step size negligible ?  (this test may be provably
!     redundant on some well-behaved arithmetics.)

         dx = abs(dznr) + abs(dzni)
         if (dx+abdir == dx) then

!     the step is negligible.  assume  zn=(dznr,dzni) is a
!     root. exit the iteration loop  do while(contin).

            contin = .false.
            cycle l1
         end if

!     now determine whether  cdir  is acceptable.

      else

!     zn  is real.  thus, attempt to divide the polynomial
!     by the linear factor (z-x); that is,

!     sum(du(i)*z**(n-i)) =
!     (z-x) * sum(d(i)*z**(n-i-1), i=0,n-1) + d(n)
!     for all z,

!     where x is a zero of  z-x.  in the code below,
!     e    is error bound for the value of polynomial at x,
!     f    is the value of the polynomial at x,
!     v    is the first derivitive of polynomial at x,
!     2*w  is the second derivative of polynomial at x, and
!     d(i) are the coefficients of quotient polynomial
!     deflat(i).

!     be sure that the overflow indicator is turned off.

         ovflow = .false.
         unflow = .false.

!     the iterate is taken as real. set the imaginary part to
!     zero.

         dzni = zero
         call c02agv(dznr,n,du,f,v,w,e,deflat)
         fn = abs(f)

!     check for overflow.

         if (ovflow) then

!     indicate that the polynomial needs to be scaled down
!     and cycle the do loop on n.  note: cauchy is not reset
!     as the cauchy bounds need not be recomputed.

            scbyex = -expdep
            cycle l0
         end if

!     check to see if  zn  is a zero or if underflow has
!     occurred.

         if (fn <= e .or. unflow) then
            if (unflow) then
               ier = 3
               if (iers /= 1) then
                  write (*,'(" ** c02agf p(z) evaluation underflow.")')
                 write (*,'(" ** if this message occurs contact nag.")')
               end if
               return           !   <== exit
            end if

!     a root has been found -- exit the iteration
!     loop  do while(contin).

            contin = .false.
            cycle l1
         end if
         unf = unflow .or. unf

!     has the function value decreased ?

         if (fn >= f0 .and. startd) then

!     no, it has not.  indicate that the laguerre step is
!     unacceptable.  (a ratio larger than theta indicates
!     that the laguerre step should be shortened.)

            ratio = bigone*theta

!     cycle iteration loop  do while(contin).

            cycle l1
         end if
         ovflow = .false.

!     find the laguerre step at  dznr.

         r = f06blf(v,f,ovflow)

!     if ovflow,  a root of polynomial is within
!     4 * n * base ** (-1022)  of  zn.

         if (ovflow) then
            unf = .true.

!     a root has been found -- exit the iteration
!     loop  do while(contin).

            contin = .false.
            cycle l1
         end if

!     compute the laguerre step  cdir  and the bound  fejer
!     at  zn .  the laguerre step and fejer bounds are computed
!     from the smaller root of a quadratic polynomial.

         call c02agt(x2n1*w,x2n*v,f,c,cf1)
         fejer = a02abf(c(1),c(2))
         r = xn2n*r
         ctemp(1) = c(1)*r + xn1
         ctemp(2) = c(2)*r
         call a02acf(c(1),c(2),ctemp(1),ctemp(2),cdir(1),cdir(2))
         abdir = a02abf(cdir(1),cdir(2))
         ratio = abdir/g
         fejer = min(rtn*abdir,fejer)

!     is the step size negligible ?

         dx = abs(dznr)
         if (dx+abdir == dx) then

!     the step is negligible.  assume  zn=(dznr,dzni) is
!     a root. exit the iteration loop  do while(contin).

            contin = .false.
            cycle l1
         end if

!     now determine whether  cdir  is acceptable.

      end if

!     repeat the iteration loop do while(contin).

      enddo l1

!     a root has been computed.  deflate the polynomial.

      if (dzni /= zero) then

!     accept zn as a complex root and deflate for a complex root.
!     put coefficients of the quotient polynomial in the du array.
!     du(0) is unchanged for the deflated polynomial.

         do i=1,n-2
            du(i) = deflat(i)
         enddo
         z(1,n) = dznr
         z(2,n) = dzni
         z(1,n-1) = z(1,n)
         z(2,n-1) = -z(2,n)
         n= n - 2
      else
         
!     accept zn as a real root and deflate for a real root.
!     put coefficients of the quotient polynomial in the du array.
!     du(0) is unchanged for the deflated polynomial.

         do i=1,n-1
            du(i) = deflat(i)
         enddo
         z(1,n) = dznr
         z(2,n) = zero
         n= n - 1
      end if

!     indicate that the cauchy region containing the smallest zeros
!     of the current polynomial has not been computed.

      cauchy = .false.

!     repeat the loop while n>2 for decreasing n.

      enddo l0

!     the polynomial is now of degree 2 or less.  determine the
!     remaining roots directly rather than iteratively.

      ovflow = .false.
      unflow = .false.
      if (n == 2) then
         call c02agt(du(0),du(1),du(2),ctemp,c)
         z(1,1) = c(1)
         z(2,1) = c(2)
         z(1,2) = ctemp(1)
         z(2,2) = ctemp(2)
      else if (n == 1) then
         z(1,1) = -du(1)/du(0)
         z(2,1) = zero
      else
         ovf = ovf .or. ovflow
         unf = unf .or. unflow

!     restore overflow and underflow indicators and enable message.

         ovflow = savo
         unflow = savu

!     provide only the relevant over/underflow messages.

         if (ovf) r = finity*finity
         if (unf) r = tiny*tiny
      end if

      return

      end subroutine c02agz

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e02acf(x,y,n,aa,m1,ref)
!-----------------------------------------------------------------------
!     mark 1 release.  nag copyright 1971
!     mark 4.5 revised
!     mark 5c revised
!     mark 9b revised. ier-361 (jan 1982)
!     mark 11.5(f77) revised. (sept 1985.)
!     mark 14c revised. ier-877 (nov 1990).
!     calculates a minimax polynomial fit to a set of data points
!     as a
!     series of chebyshev polynomials.
!-----------------------------------------------------------------------
      implicit none
      include 'evolcom.num'

!     .. scalar arguments ..
      real(8):: ref
      integer,intent(in):: m1,n
!     .. array arguments ..
      real(8),dimension(m1):: aa
      real(8),dimension(n):: x,y
!     .. local scalars ..
      real(8):: abshi,ai,ai1,d,denom,h,hi,himax,hmax,prevh,rhi,rhi1,xi
     $     ,xj,xnexth
      integer:: i,i1,i2,ij1,imax,iri,irj,j,j1,k,m,m2
      logical:: loopexit
!     .. local arrays ..
      real(8),dimension(100):: a,rh,rx
      integer,dimension(100):: ir
!-----------------------------------------------------------------------
!     enforced hard fail for out-of-bounds parameters
!if (m1 >= n .or. m1 >= 100) ifail = p01abf(0,1,srname,0,p01rec)
      if (m1 >= n .or. m1 >= 100) then
         rewind(222)
         write(222,*) 'crash in nag e02acf'
         stop
      endif
      do i = 2, n
!if (x(i) <= x(i-1)) ifail = p01abf(0,2,srname,0,p01rec)
         if (x(i) <= x(i-1)) then
            rewind(222)
            write(222,*) 'crash in nag e02acf'
            stop
         endif
      enddo
      m2 = m1 + 1
      m = m1 - 1
      prevh = -one
      ir(1) = 1
      ir(m2) = n
      d = dble(n-1)/dble(m1)
      h = d
      if (m /= 0) then
         do i = 2, m1
            ir(i) = int(h+half) + 1
            h = h + d
         enddo
      endif
      l0: do
      h = -one
      do i = 1, m2
         iri = ir(i)
         rx(i) = x(iri)
         a(i) = y(iri)
         rh(i) = -h
         h = -h
      enddo
      l1:do j = 1, m1
      i1 = m2
      ai1 = a(i1)
      rhi1 = rh(i1)
      i = m2
      do
         i= i - 1
         ij1 = i - j + 1
         denom = rx(i1) - rx(ij1)
         ai = a(i)
         rhi = rh(i)
         a(i1) = (ai1-ai)/denom
         rh(i1) = (rhi1-rhi)/denom
         i1 = i
         ai1 = ai
         rhi1 = rhi
         if (i-j <= 0) then
            cycle l1
         else
            cycle
         endif
      enddo
      enddo l1
      h = -a(m2)/rh(m2)
      do i = 1, m2
         a(i) = a(i) + rh(i)*h
      enddo
      if (m /= 0) then
         j= m1
         do
            j = j - 1
            xj = rx(j)
            i = j
            ai = a(i)
            j = j + 1
            do i1 = j, m1
               ai1 = a(i1)
               a(i) = ai - xj*ai1
               ai = ai1
               i = i1
            enddo
            j = j - 1
            if (j-1 <= 0) then
               exit
            endif
         enddo
      endif
      hmax = abs(h)
      if (hmax <= prevh) then
         a(m2) = -hmax
         do i = 1, m1
            aa(i) = a(i)
         enddo
         ref = a(m2)
         return                 !   <== exit
      endif
      a(m2) = hmax
      prevh = hmax
      imax = ir(1)
      himax = h
      j = 1
      irj = ir(j)
      do i = 1, n
         if (i /= irj) then
            xi = x(i)
            hi = zero
            k = m2
            do
               k = k - 1
               hi = hi*xi + a(k)
               if (k-1 <= 0) then
                  exit
               endif
            enddo
            hi=hi-y(i)
            abshi = abs(hi)
            if (abshi <= hmax) then
               cycle
            endif
            hmax = abshi
            himax = hi
            imax = i
            cycle
         endif
         if (j >= m2) then
            cycle
         endif
         j = j + 1
         irj = ir(j)
      enddo
      if (imax == ir(1)) then
         do i = 1, m1
            aa(i) = a(i)
         enddo
         ref = a(m2)
         return                 !   <== exit
      endif
      loopexit=.false.
      do i = 1, m2
         if (imax < ir(i)) then
            loopexit=.true.
            exit
         endif
      enddo
      if (.not. loopexit) then
         i = m2
      endif
      i2 = int(dble(i)*half)
      i2 = i - 2*i2
      xnexth = h
      if (i2 == 0) xnexth = -h
      if (himax*xnexth >= 0.0d0) then
         ir(i) = imax
         cycle l0
      endif
      if (imax < ir(1)) then
         j1 = m2
         j= m2
         do
            j = j - 1
            ir(j1) = ir(j)
            j1 = j
            if (j-1 <= 0) then
               exit
            endif
         enddo
         ir(1) = imax
         cycle l0
      endif
      if (imax > ir(m2)) then
         j= 1
         do j1 = 2, m2
            ir(j) = ir(j1)
            j = j1
         enddo
         ir(m2) = imax
         cycle l0
      endif
      ir(i-1) = imax
      cycle l0
      enddo l0

      return

      end subroutine e02acf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      real(8) function f06blf(a,b,fail)
!-----------------------------------------------------------------------
!     mark 12 release. nag copyright 1986.
!     mark 18 revised (thread safety). (sep 1996).
!-----------------------------------------------------------------------
      implicit none
      include 'evolcom.num'

!     .. scalar arguments ..
      real(8),intent(in):: a,b
      logical,intent(out):: fail
!     ..

!     f06blf returns the value div given by

!     div = ( a/b                 if a/b does not overflow,
!     (
!     ( 0.0                 if a .eq. 0.0,
!     (
!     ( sign( a/b )*flmax   if a .ne. 0.0  and a/b would overflow,

!     where  flmax  is a large value, via the function name. in addition if
!     a/b would overflow then  fail is returned as true, otherwise  fail is
!     returned as false.

!     note that when  a and b  are both zero, fail is returned as true, but
!     div  is returned as  0.0. in all other cases of overflow  div is such
!     that  abs( div ) = flmax.

!     when  b = 0  then  sign( a/b )  is taken as  sign( a ).

!     nag fortran 77 o( 1 ) basic linear algebra routine.

!     -- written on 26-october-1982.
!     sven hammarling, nag central office.

!     .. parameters ..
      real(8),parameter:: flmax =4.49423283715580d+307,flmin
     $     =2.22507385850721d-308
!     .. local scalars ..
      real(8):: absb,div
!-----------------------------------------------------------------------
      if (a == zero) then
         div = zero
         if (b == zero) then
            fail = .true.
         else
            fail = .false.
         end if
      else

         if (b == zero) then
            div  =  sign(flmax,a)
            fail = .true.
         else
            absb = abs(b)
            if (absb >= one) then
               fail = .false.
               if (abs(a) >= absb*flmin) then
                  div = a/b
               else
                  div = zero
               end if
            else
               if (abs(a) <= absb*flmax) then
                  fail = .false.
                  div  =  a/b
               else
                  fail = .true.
                  div  = flmax
                  if (((a < zero) .and. (b > zero)) .or. ((a > zero)
     $                 .and. (b < zero))) div = -div
               end if
            end if
         end if
      end if

      f06blf = div
      return

      end function f06blf
