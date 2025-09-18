      SUBROUTINE Quintic(omkyle,vmlt,lmlt,val)
C***********************************************************************
C  Safeguarded solver for z from:  c + 5 x^2 - 2 x^5 = 0
C  with  c = 18 / ( (5*pi*Roc)^2 ),  Roc = vmlt / (2*omkyle*lmlt)
C  Returns positive real root x = z.  Robust to small/large inputs.
C***********************************************************************
      IMPLICIT NONE
      REAL*8, INTENT(IN)    :: omkyle, vmlt, lmlt
      REAL*8, INTENT(INOUT) :: val

      REAL*8 :: pi, Roc, c, denom, x, fun, funp
      REAL*8 :: a, b, fa, fb, xnew, err
      REAL*8 :: tol, z0, tiny, rocmin, bmax
      INTEGER :: n, nmax

      pi    = 3.141592653589793238462643383279D0
      z0    = (5.0D0/2.0D0)**(1.0D0/3.0D0)   ! ~1.357, Ω -> 0 limit
      tiny  = 1.0D-300
      rocmin= 1.0D-12
      bmax  = 1.0D6
      tol   = 1.0D-12
      nmax  = 100

C----- Degenerate inputs: return z0 (non-rotating limit)
      IF (lmlt .LE. 0.0D0) THEN
         val = z0
         RETURN
      ENDIF

C----- Build Roc safely
      denom = 2.0D0*omkyle*lmlt
      IF (DABS(denom) .LT. tiny) THEN
C        effectively Ω -> 0  (or vanishing mixing length denominator)
         Roc = 1.0D300               ! treat as infinite Roc => c ~ 0
      ELSE
         Roc = DABS(vmlt/denom)      ! Roc >= 0, sign irrelevant for c
      ENDIF
      IF (Roc .LT. rocmin) Roc = rocmin

C----- c = 18 / ( (5*pi*Roc)^2 ), safe against overflow/underflow
      c = 18.0D0 / ( (5.0D0*pi*Roc)**2 )

C----- Bracket the positive root: f(0)=c>0, find b with f(b)<0
      a  = 0.0D0
      fa = c

C----- Initial guess scaling ~ c^(1/5) behaviour for large c
      x  = 1.1D0
      IF (c .GT. 2.0D1) x = 0.65D0 * c**0.2D0
      b  = DMAX1(x, 1.0D0)

      fb = c + 5.0D0*b*b - 2.0D0*b**5
      DO WHILE (fb .GT. 0.0D0 .AND. b .LT. bmax)
         b  = 2.0D0*b
         fb = c + 5.0D0*b*b - 2.0D0*b**5
      END DO

C----- If we failed to bracket (pathological), fall back to z0
      IF (fb .GT. 0.0D0) THEN
         val = z0
         RETURN
      ENDIF

C----- Ensure starting x is inside [a,b]
      IF (x .LE. a .OR. x .GE. b) x = 0.5D0*(a+b)

C----- Safeguarded Newton with bisection fallback
      DO n = 1, nmax
         fun  = c + 5.0D0*x*x - 2.0D0*x**5
         funp = 10.0D0*x*(1.0D0 - x**3)

C-------- Newton step
         IF (funp .NE. 0.0D0) THEN
            xnew = x - fun/funp
         ELSE
            xnew = 0.5D0*(a+b)
         ENDIF

C-------- If step leaves bracket or is non-finite/negative, bisect
         IF (xnew .LE. a .OR. xnew .GE. b .OR. xnew .NE. xnew
     &       .OR. xnew .LE. 0.0D0) THEN
            xnew = 0.5D0*(a+b)
         ENDIF

C-------- Update bracket
         fun = c + 5.0D0*xnew*xnew - 2.0D0*xnew**5
         IF (fun .GT. 0.0D0) THEN
            a  = xnew
            fa = fun
         ELSE
            b  = xnew
            fb = fun
         ENDIF

         err = DABS(b - a)
         x   = xnew
         IF (err .LE. tol*(1.0D0 + x)) EXIT
      END DO

      val = x

C----- Final safety: ensure positive, finite z
      IF (val .NE. val .OR. val .LE. 0.0D0 .OR. val .LT. z0) THEN
         val = z0
      ENDIF

      RETURN
      END
