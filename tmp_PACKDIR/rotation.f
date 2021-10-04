
***********************************************************************

      SUBROUTINE potential

************************************************************************
*     Subroutines computing the correction factors fp and ft to be     *
*     applied to the stellar structure equations in case of rotation   *
*                                                                      *
*     Source : Endal & Sofia, 1976, ApJ 210, 184 Appendix A            *
*     -------                                                          *
*                                                                      *
* $LastChangedDate:: 2016-05-13 16:35:05 +0200 (Ven, 13 mai 2016)    $ *
* $Author:: palacios                                                 $ *
* $Rev:: 80                                                          $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.rot'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'
      include 'evolcom.eng'

      integer i,ip,k

      double precision rom
      double precision gp,gpi,ri,rip
      double precision r1,r2,epsmin,h1,hmin
      double precision r0,rtsafe,psiint,gpsi,spgp,spgm
      double precision rae,rov,eta2,wr,rpsi,mpsi,aeta
      double precision intfac,rrs,rrs1
      double precision sfequi,intsin3th
      double precision geffrot,Rp2,angleL

      common /radi/ rae(intmax),rov,eta2,wr,rpsi,mpsi,aeta
      common /spgp/ spgp(nsh),sfequi(nsh)
c..   test rotaton 27/11/2014
      common /geffective/ geffrot(nsh),Rp2(nsh)

      dimension r0(nsh),gp(nsh),gpi(nsh),spgm(nsh)

      save

      forall (i = 1:nsh) geffrot(i) = 0.d0
      forall (i = 1:nsh) Rp2(i) = 0.d0

      angleL = sqrt(pw23)

c..   main loop
      eta2 = -1.00d-20
      intfac = 0.d0
      epsmin = 1.d-4
      r0(1) = 0.d0

      do i = 1,nmod1,1
         ip = i+1
         rom = m(ip)*0.75d0/(pi*r(ip)**3)
         rov = ro(i)/rom
         ri = lnr(i)
         rip = lnr(ip)
         h1 = (rip-ri)/8.d0
         hmin = (rip-ri)/300.d0

c..   solve for the radau equation : compute eta2
         call radau (eta2,ri,rip,rov,epsmin,h1,hmin)

c..   compute r0 (find zero of function apot)
         mpsi = m(ip)
         rpsi = r(ip)
         wr = omega(ip)
         r1 = 0.01d0*rpsi
         r2 = 2.00d0*rpsi
         r0(ip) = rtsafe(r1,r2,1.0d-8,i)

c..   compute the integral factor needed for the potential
         intfac = psiint(ro(ip),ro(i),mpsi,wr,eta2,r0(i),r0(ip))+intfac

c..   compute g and 1/g for the different values of the angle theta
         do k = 1,intmax
            gp(k) = gpsi(r0(ip),intfac,k)
            gpi(k) = 1.d0/gp(k)
         enddo

c..   TEST ROTATION 27/11/2014
         geffrot(i) = gp(intmax)
         Rp2(i) =  r0(i)*(1.d0-aeta*plegendre(62))

c  compute Spsi*<g> and Spsi*<1/g>
         rrs = rae(1)**2*sint(1)

         rrs1 = rae(intmax)**2*sint(intmax)
         sfequi(ip) = (rrs+rrs1)*dtheta2
         spgp(ip) = (gp(1)*rrs+gp(intmax)*rrs1)*dtheta2
         spgm(ip) = (gpi(1)*rrs+gpi(intmax)*rrs1)*dtheta2
c         gmr(ip) = (gp(1)+gp(intmax))*dtheta2
c _LA_ à revoir         gmr(ip) = (gp(intmax)*sint(intmax)**3)*dtheta2
         intsin3th = dtheta2
         do k = 2,intmax-1
            rrs = rae(k)*rae(k)*sint(k)*dtheta
            spgp(ip) = spgp(ip)+gp(k)*rrs
            spgm(ip) = spgm(ip)+gpi(k)*rrs
            sfequi(ip) = sfequi(ip)+rrs
c _LA_ à revoir            gmr(ip) = gmr(ip)+gp(k)*sint(k)**3*dtheta
         enddo
c _LA_ à revoir         gmr(ip) = gmr(ip)*1.5d0
c         gmr(ip) = gmr(ip)/(0.5d0*pi)
         sfequi(ip) = pim4*sfequi(ip)
         

c  compute fp and ft


         fp(ip) = rpsi**4/(g*mpsi*spgm(ip))
         ft(ip) = rpsi**4/(spgp(ip)*spgm(ip))
c         ft(ip) = min(ft(ip),1.d0)
c         fp(ip) = min(fp(ip),1.d0)

      enddo

      ft(2) = ft(3)
      fp(2) = fp(3)
      ft(1) = ft(2)
      fp(1) = fp(2)
      ft(nmod-2) = ft(nmod-3)
      ft(nmod1) = ft(nmod-2)
      ft(nmod) = ft(nmod1)
      fp(nmod-2) = fp(nmod-3)
      fp(nmod1) = fp(nmod-2)
      fp(nmod) = fp(nmod1)
c      forall (i = 1:nmod) fp(i) = 0.98d0
c      forall (i = 1:nmod) ft(i) = 0.99d0
c      print *,'fprotation = ',fp(1:nmod)
c      print *,'ftrotation = ',ft(1:nmod)


      return
      end


***********************************************************************

      SUBROUTINE radau (ystart,x1,x2,rov,epsmin,h1,hmin)

***********************************************************************
*     radau equation solved by runge-kutta method

      implicit none

      integer i,nok,nbad,maxstp

      double precision ystart,x1,x2,rov,epsmin,h1,hmin
      double precision x,h,hnext,hdid,y,yscal,dydx,tiny

      parameter (maxstp = 10000,tiny = 1.d-30)

      x = x1
      h = dsign(h1,x2-x1)
      nok = 0
      nbad = 0
      y = ystart
      do i = 1,maxstp
         dydx = 6.d0*(1.d0-rov*(y+1.d0))-y*(y-1.d0)
         yscal = dabs(y)+dabs(h*dydx)+tiny
         if ((x+h-x2)*(x+h-x1).gt.0.d0) h = x2-x
         call rkqc (y,dydx,x,h,epsmin,yscal,hdid,hnext,rov)
         if (hdid.eq.h)then
            nok = nok+1
         else
            nbad = nbad+1
         endif
         if ((x-x2)*(x2-x1).ge.0.d0) then
            ystart = y
            return
         endif
         if (dabs(hnext).lt.hmin) hnext = hmin
         h = hnext
      enddo

      return
      end


***********************************************************************

      SUBROUTINE rkqc (y,dydx,x,htry,epsmin,yscal,hdid,hnext,rov)

***********************************************************************
*     Runge-Kutta driver with quality control

*     Compares results from one big step and 2 small steps
*     Evaluates the next time step (or the current one if too big)
*     on how different the results are.

*     Source: http://optics.unige.ch/vdm/marel_files/numerical_recipies/rkqc.for

      implicit none

      double precision pgrow,pshrnk,fcor,h,hh,ytemp,ysav,dysav,safety
      double precision xsav,errmax,errcon
      double precision y,dydx,x,htry,epsmin,yscal,hdid,hnext,rov

      parameter (safety = 0.9d0,errcon = 6.d-4)

      pgrow = -0.20d0
      pshrnk = -0.25d0
      fcor = 2.d-1/3.d0
      xsav = x

      ysav = y
      dysav = dydx
      h = htry
 
c..   Solution with 2 consecutive small steps
 10   hh = 0.5d0*h
      call rk4 (ysav,dysav,xsav,hh,ytemp,rov)
      x = xsav+hh
      dydx = 6.d0*(1.d0-rov*(ytemp+1.d0))-ytemp*(ytemp-1.d0)
      call rk4 (ytemp,dydx,x,hh,y,rov)
      x = xsav+h
c..   Solution with one big step (h = 2*hh) 
      call rk4 (ysav,dysav,xsav,h,ytemp,rov)
      errmax = 0.d0
c..   Modifying current or next step according to
c..   the relative difference of the results y and ytemp
      ytemp = y-ytemp
      errmax = max(errmax,dabs(ytemp/yscal))
      errmax = errmax/epsmin
      if (errmax.gt.1.d0) then
         h = safety*h*(errmax**pshrnk)
         goto 10
      else
         hdid = h
         if (errmax.gt.errcon) then
            hnext = safety*h*(errmax**pgrow)
         else
            hnext = 4.d0*h
         endif
      endif
      y = y+ytemp*fcor

      return
      end


***********************************************************************

      SUBROUTINE rk4(y,dydx,x,h,yout,rov)

***********************************************************************
*     runge-kutta solver


      implicit none

      double precision y,yt,dydx,x,h,yout,hh,h6,xh,dyt,dym,rov

      hh = h*0.5d0
      h6 = h/6.d0
      xh = x+hh
      yt = y+hh*dydx
      dyt = 6.d0*(1.d0-rov*(yt+1.d0))-yt*(yt-1.d0)
      yt = y+hh*dyt
      dym = 6.d0*(1.d0-rov*(yt+1.d0))-yt*(yt-1.d0)
      yt = y+h*dym
      dym = dyt+dym
      dyt = 6.d0*(1.d0-rov*(yt+1.d0))-yt*(yt-1.d0)
      yout = y+h6*(dydx+dyt+2.d0*dym)

      return
      end


***********************************************************************

      SUBROUTINE apot (r0,dr0,drpsidr0)

***********************************************************************
*     compute variable A (Eq. A10) and its derivatives

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'

      double precision r0,dr0,drpsidr0,rfk,coeff,a1
      double precision rae,rov,eta2,wr,rpsi,mpsi,aeta

      common /radi/ rae(intmax),rov,eta2,wr,rpsi,mpsi,aeta

      if (r0.gt.1.d30) r0 = 1.d30
      a1 = wr*wr*r0/(0.6d0*g*mpsi*(2.d0+eta2))
      aeta = r0*a1
      coeff = max(1.0d-20,1.d0+(0.6d0-pw235*aeta)*aeta**2)**pw13
      rfk = r0*coeff
      dr0 = rfk-rpsi
      drpsidr0 = coeff+r0*(1.2d0-pw635*aeta)*pw23*a1*aeta/coeff**2
      if (drpsidr0.gt.1.d30) drpsidr0 = 1.d30


      return
      end


***************************************************************************

      FUNCTION rtsafe (x1,x2,racc,ishell)

***************************************************************************
*   search for r0, solution of rpsi = r0*[1+3/5*A(ro)^2-2/35*A(r0)^3]^1/3
*   Using a combination of Newton-Raphson and bisection, RTSAFE finds the
*   root of a function bracketed between x1 and x2
*   Source: Numerical Recipes pp 359  

      implicit none

      integer j,maxit,ishell

      double precision x1,x2,f,df,fl,fh,dx,xacc,xh,xl,swap,rtsafe,racc,
     &     temp,dxold

      parameter (maxit = 300)

      call apot (x1,fl,df)
      call apot (x2,fh,df)
      xacc = racc*(dabs(x1)+dabs(x2))
      if (fl*fh.ge.0.d0) then
         print *,'pb shell',ishell
         stop 'rtsafe : root must be bracketed'
      endif
      if (fl.lt.0.d0) then
         xl = x1
         xh = x2
      else
         xh = x1
         xl = x2
         swap = fl
         fl = fh
         fh = swap
      endif
      rtsafe = 0.5d0*(x1+x2)
      dxold = dabs(x2-x1)
      dx = dxold
      call apot (rtsafe,f,df)
      do j = 1,maxit
         if (((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0.d0
     &        .or.dabs(2.d0*f).gt.dabs(dxold*df)) then
            dxold = dx
            dx = 0.5d0*(xh-xl)
            rtsafe = xl+dx
            if (xl.eq.rtsafe) return
         else
            dxold = dx
            dx = f/df
            temp = rtsafe
            rtsafe = rtsafe-dx
            if (temp.eq.rtsafe) return
         endif
         if (dabs(dx).lt.xacc) return
         call apot (rtsafe,f,df)
         if (f.lt.0.d0) then
            xl = rtsafe
            fl = f
         else
            xh = rtsafe
            fh = f
         endif
      enddo


      return
      end


***********************************************************************

      FUNCTION psiint (roip,roim,xmip,omegaip,etaip,rnim,rnip)

***********************************************************************

      implicit none

      include 'evolcom.cons'

      double precision roip,roim,xmip,omegaip,etaip,rnim,rnip,psiint
      double precision rnipv,rnimv

      rnipv = rnip*1.d-5
      rnimv = rnim*1.d-5
c..   ori
c      psiint = 0.5d0*(roip+roim)*(rnipv**8-rnimv**8)/xmip
c      psiint = psiint*1.d40
c      psiint = psiint*(5.d0+etaip)/(2.d0+etaip)*0.125d0*omegaip*omegaip
      psiint = 0.5d0*(roip+roim)*(rnipv**7-rnimv**7)/xmip
      psiint = psiint*1.d35
      psiint = psiint*(5.d0+etaip)/(2.d0+etaip)*pw17*omegaip*omegaip

      return
      end


**********************************************************************

      FUNCTION gpsi(r0,intfac,itheta)

***********************************************************************
*     compute local effective gravity (as a function of theta)

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.rot'

      integer itheta

      double precision dpsidr,dpsidtheta,r,ri,ri2,wr2,r0,intfac,gpsi
      double precision rae,rov,eta2,wr,rpsi,mpsi,aeta

      common /radi/ rae(intmax),rov,eta2,wr,rpsi,mpsi,aeta

      r = r0*(1.d0-aeta*plegendre(itheta))
      ri = 1.0d0/r
      ri2 = ri*ri
      wr2 = wr*wr
            
      rae(itheta) = r
      dpsidr = -g*mpsi*ri2+pim4*ri2*ri2*plegendre(itheta)*intfac+wr2*r*
     &     sint2(itheta)
      dpsidtheta = (pim4*ri2*ri*intfac+wr2*r*r)*cost(itheta)*
     &     sint(itheta)
      gpsi = dsqrt(dpsidr*dpsidr+dpsidtheta*dpsidtheta*ri2)

c..   Meynet and Maeder (1997, A&A, 321, 465)
c     dwrdpsi = wr*r*r*sint2(itheta)
c     gpsi = (1.d0-r*r*sint2(itheta)*wr*dwrdpsi)*gpsi

      end
