
************************************************************************

      SUBROUTINE Quintic(omkyle,vmlt,lmlt,val)

************************************************************************
*   Newton-Raphson solver for the solution to the quintic nondiffusive *
*   rotating convection model of Augustson & Mathis 2019,              *
*                                                                      *
*     Author: Kyle Augustson (CEA Saclay)                              *
*     Date of creation: 26 October 2019.                               *
*     Adaptation to Starevol : TD (04/11/2019)                         *
************************************************************************ 
!     omkyle is the local angular velocity
!     vmlt is the local mixing length velocity
!     lmlt is the local mixing length
!     val is the return value of the 3/2 root of the total wave vector
!
!     Note TD: routine to compute parameter z (see eq. 69, 70)      
************************************************************************
      
      Implicit None
      
      Real*8, Intent(In) :: omkyle, vmlt, lmlt
      Real*8, Intent(InOut) :: val
      Real*8 :: c, err, tol, x, guess, tmp, fun, funp, pi, Roc
      Integer :: n, nmax

      pi = 3.141592653589793238462643383279d0
      
!     Convective Rossby number (eq. 22)      
      Roc = vmlt/(2d0*omkyle*lmlt)
! c is the third term of equation (46)      
      c = 18d0/(5d0*pi*Roc)**2
      
      if (c.gt.20d0) then
         guess = 0.65d0*(c**(0.2d0))
      else
         guess = 1.1d0
      end if
      nmax = 100
      x = guess
!     Newton-Raphson algorithm
      tol = 1d-15
      err = 1d0
      n=0
      do while ((err>tol) .and. (n<nmax))
         fun = c+5d0*x**2-2d0*x**5
         funp = 10d0*x*(1d0-x**3)
         tmp = x - fun/funp
         err = abs(tmp-x)
         x = tmp
         n = n + 1
      end do
      val = tmp

      return
      end
