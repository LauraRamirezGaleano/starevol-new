      program factor

c cd /home/siess/EVOL/CODE/GAG
c f90 -I./INCL factor.f rmodpar.o rininet.o -o factor.e

      implicit none

      include 'evolpar.chem'
      include 'evolpar.conv'
      include 'evolpar.str'

      include 'evolcom.cons'
      include 'evolcom.eos'
      include 'evolcom.ion'
      include 'evolcom.nuc'
      include 'evolcom.shk'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer i,ii,k,l,nro

      double precision fscr,fscr_ec
      double precision rs
      double precision screen,screen_ec
      double precision dzi
      double precision y,zx,zy,zbar,ztilde,l0,l086

      common /screening/ y(nsp),zbar,ztilde,zx,l0,l086
      common /funcscr/ fscr(nsh,nscr)

      dimension fscr_ec(nsh,nbetaec)

      external screen
c      external screen_ec

      pw17 = 1.d0/7.d0
      pw14 = 1.d0/4.d0
      pw13 = 1.d0/3.d0
      pw16 = 1.d0/6.d0
      pw23 = 2.d0/3.d0
      pw34 = 3.d0/4.d0
      pw43 = 4.d0/3.d0
      pw53 = 5.d0/3.d0
      nmod = 50
      nro = 11

      print *,nmod,nro
      open (unit = 99,file = 'starevol_gag261.par',status = 'unknown')
      call rmodpar
      call rininet

      do k = 1,nmod
         do l = 1,nsp
            xsp(k,l) = 1.D-50
         enddo
         xsp(k,ic12) = 1.d0
      enddo

      do k = 1,nmod
         rhpsi(k) = 1.d0
         t(k) = 10.d0**(dble(k)*4.d-2+7)
      enddo

      do i = 1,nscr
         do k = 1,nmod
            fscr(k,i) = 1.d0
         enddo
      enddo

***   Charged reaction factors      
      do ii = 1,nro
         ro(ii) = 10.d0**(9d0-dble(ii-1)*0.5d0)         
         do k = 1,nmod
            if (t(k).gt.tnucmin) then
               ksh = k
               rok = ro(ii)
               tk = t(k)
               rhpsik = rhpsi(k)
               muiinvk = muinv(k)
               do l = 1,nsp
                  x(l) = max(1.d-50,xsp(k,l))
                  y(l) = x(l)/anuc(l)
               enddo
               zx = 0.d0
               zbar = 0.d0
               ztilde = 0.d0            
               do i = 1,nis
                  zx = zx+y(i)
                  zy = znuc(i)*y(i)
                  zbar = zbar+zy
                  ztilde = ztilde+znuc(i)*zy
               enddo
               zbar = zbar/zx
               ztilde = sqrt(ztilde/zx+rhpsik*zbar)
               l0 = 1.88d8*sqrt(rok*zx/tk**3)
               l086 = 0.38d0*l0**0.86d0
c..   nbz-1 because nbz = heavy
c            do i = 1,nbz-1
c               ii = nbz-1+i
c               dzi = dble(i)
c               fscr(k,i) = screen (dzi,1.d0)
c               fscr(k,ii) = screen (dzi,2.d0)
c            enddo
c            fscr(k,2*nbz-1) = screen (6.d0,6.d0)
c            fscr(k,2*nbz) = screen (8.d0,8.d0)
               fscr(k,1) = screen (6.d0,6.d0)
            endif
         enddo
      enddo

***   Electron capture reactions
c       do k = 1,nmod
c          if (t(k).gt.tnucmin) then
c             do i = 1,nbetaec
c                ii = k1(jbeta(nbeta+i))
c                rs = 1.388d-2*(anuc(ii)*ro(k)/znuc(ii))**pw13
c                fscr_ec(k,i) = screen_ec (rs)
c             enddo
c          endif
c       enddo

      end


************************************************************************

      DOUBLE PRECISION FUNCTION screen (z1,z2)

************************************************************************
* Calculate the screening factor for the thermonuclear reaction rates  *
* Input : nuclide charges z1 and z2                                    *
* Graboske etal 1973,ApJ,181,457 and Cox & Guili 17.15 for definitions *
************************************************************************

      implicit none

      include 'evolpar.chem'
      include 'evolpar.conv'
      include 'evolpar.str'

      include 'evolcom.cons'
      include 'evolcom.nuc'
      include 'evolcom.shk'

      integer i

      double precision z1,z2,z12,zx,tfermi,tion
      double precision a1,a2,tij,gij,f,beta,aij,y1,y2
      double precision y,zbar,ztilde,zmean,h120,l0,l086,l12,h120a

      character p*1

      common /screening/ y(nsp),zbar,ztilde,zx,l0,l086


      h120 = 0.d0
      l12 = l0*z1*z2*ztilde

      p = 'x'
*___________________
***   weak screening
*-------------------
      if (l12.lt.0.1d0) then
         h120 = l12
         p = 'w'
*___________________________
***   intermediate screening
*---------------------------
      elseif (l12.ge.0.1d0.and.l12.le.5.d0) then
         zmean = 0.d0
         do i = 1,nis
            zmean = zmean+(znuc(i)**1.58d0*y(i))
         enddo
         zmean = zmean/zx
         h120 = zmean*ztilde**(-0.58d0)*zbar**(-0.28d0)
         h120 = h120*((z1+z2)**1.86d0-z1**1.86d0-z2**1.86d0)*l086
         p = 'i'
      endif
            
*_____________________
***   strong screening
*---------------------
      if (l12.ge.2.d0) then
         z12 = z1+z2
         h120a = (z12**pw53-z1**pw53-z2**pw53)
     &        +.316d0*zbar**pw13*(z12**pw43-z1**pw43-z2**pw43)
     &        +0.737d0/(zbar*l0**pw23)*(z12**pw23-z1**pw23-z2**pw23)
         h120a = h120a*0.624d0*zbar**pw13*l0**pw23
         if (l12.le.5.d0.and.h120.lt.h120a) then
            p = 'i'
         else
            p = 's'
         endif
         if (l12.le.5.d0) h120 = min(h120,h120a)
         if (l12.gt.5.d0) h120 = h120a
      endif
      screen = exp(min(220.d0,h120))

c     return

*____________________________
***   extremely dense plasmas
*----------------------------

c     condition : Tion < Tk << Tfermi
      Tfermi = 5.9302d9*(sqrt(1.0d0+1.018d4*(rok/zbar)**pw23)-1.0d0)
      Tion = 2.173d2*rok**pw23/(max(z1,z2)*2.d0)**pw53
      gij = 6.93361274d-1*zbar**2*l0**pw23
c..   approximate Ai = 2*Zi
      if (z1.eq.1.d0) then
         a1 = 1.d0
      else
         a1 = 2.d0*z1
      endif
      if (z2.eq.1.d0) then
         a2 = 1.d0
      else
         a2 = 2.d0*z2
      endif
      aij = a1*a2/(a1+a2)
      tij = 4.2487d3*(aij*(z1*z2)**2/tk)**pw13
      beta = 3*gij/tij
      f = (0.0455d0*beta+0.348d0*beta**3+9.49d0*beta**6-0.123d0*beta**12
     &     +0.101d0*beta**13)/(1.d0+1.d2*beta**4+0.267d0*beta**12)
      write (*,100) rok,tk,log10(gij),log10(tij),log10(beta)
     &     ,log10(screen),log10(exp(1.25d0*gij-tij*f)),log10(exp(1.25d0
     &     *gij-tij*f)/screen)

c     screen = exp(1.25d0*gij-tij*f)

 100  format (8(1x,1pe11.4))

      return
      end

