      SUBROUTINE splineatmb(x,y,n,yp1,ypn,y2,y1)
***********************************************************************************************
      implicit none
      INTEGER n,NMAX
      double precision yp1,ypn,x(n),y(n),y2(n)
      double precision a1,a2,a3,y1(n)
      PARAMETER (NMAX=5000)
!            Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with x1 < x2 < ... < xN, and given values yp1 and ypn for the first derivative of the inter- polating function at points 1 and n, respectively, this routine ret
!            Parameter: NMAX is the largest anticipated value of n.
      INTEGER i,k
      double precision p,qn,sig,un,u(n)
      if (yp1.gt..99e30) then
         y2(1)=0.
         u(1)=0.
      else
         y2(1)=-0.5
!     The lower boundary condition is set either to be “natural”
!     or else to have a specified first derivative.
!     j
         u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      y1(1) = 0.d0
      y1(n) = (y(n)-y(n-1))/(x(n)-x(n-1))

      do i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.
         y2(i)=(sig-1.)/p
!la         if (x(i)==x(i-1)) then
         if (x(i)==x(i-1)) then
            a1 = (x(i)-x(i+1))/(x(i-1)-x(i+1))
            a2 = (2.*x(i)-x(i-1)-x(i+1))/(x(i)-x(i+1))
            a3 = 0.d0
            u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i)))/
     &           (x(i+1)-x(i-1))-sig*u(i-1))/p
         else
            a1 = (x(i)-x(i+1))/((x(i-1)-x(i))*(x(i-1)-x(i+1)))
            a2 = (2.*x(i)-x(i-1)-x(i+1))/((x(i)-x(i-1))*(x(i)-x(i+1)))
            a3 = (x(i)-x(i-1))/((x(i+1)-x(i))*(x(i+1)-x(i-1)))
            u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &           /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
         endif
         y1(i) = a1*y(i-1)+a2*y(i)+a3*y(i+1)

      enddo
      if (ypn.gt..99e30) then
         qn=0.
         un=0.
      else
         qn=0.5
!     The upper boundary condition is set either to be “natural”
!     or else to have a specified first derivative.
         un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      return
      END
