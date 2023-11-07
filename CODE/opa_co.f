

************************************************************************

      SUBROUTINE opa_co (tk,rok,muiinvk,x,ksh,kap1,kapr1,kapt1,error)

************************************************************************
*     compute opacity and derivatives in cool surface layers           *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none
      
      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.data'
      include 'evolcom.ion'
      include 'evolcom.nuc'

      logical logicNR,debug

      integer i,j,k,ij,ksh
      integer order(neqchimmax),indx(neqchimmax+natomsmax)
      integer error,neqo,n1,info
      integer imol,im,ip,jcur,jmol
      integer ir1,nr1,ir2,nr2,ip1
      integer iterelec

C for test
c      double precision rapCO
C 
C for print
c      double precision nnh,nno2,nno,nch,nn2,no,nn,nc
C
      double precision x(nsp),tk,rok,muiinvk
      double precision xxh,xxhe,xxz
      double precision theta,rhompinv
      double precision xmol(neqchimmax+natomsmax),
     &     xmol0(neqchimmax+natomsmax)
      double precision logKp(neqchimmax),logK(neqchimmax)
      double precision kap1,kapr1,kapt1
      double precision logKt,tempnbr
      double precision kdiss,temp,tempbis,nel0im
      double precision errf,tolf,errx,tolx
      double precision fjac(neqchimmax+natomsmax,neqchimmax+natomsmax),
     &     fvec(neqchimmax+natomsmax),p(neqchimmax+natomsmax)
      double precision kdisso(neqchimmax)
      double precision nelec,zscale,lntk,logt,betaev,delta,nelecb
      double precision y,a,b,pel
      double precision nh,nh2,nh2o,noh,nco,ncn,nc2
      double precision Kion(nioniz),ntot(nioniz),lnQ(2),xnsp(nioniz)
      double precision expq
      double precision nhtot
      double precision tt,tm1,t2,tm2,t3,t4,tm4,t5,tm5,t6,tm6,t7,t10
      double precision tm7,t9,t05
      double precision rm1,r05,r025,rm075
      double precision temp1,temp2,temp3,temp4,temp5,temp6
      double precision temp1r,temp3r
      double precision temp1t,temp2t,temp3t,temp4t,temp5t,temp6t
      double precision temp3a
      double precision temp7,temp7r,temp7t
      double precision temp8,temp8r,temp9r,temp9t
      double precision temp10,temp10t,temp11
      double precision temp13r
      double precision temp11t,temp12,temp12r,temp12a,temp12t,temp13a
      double precision temp13b,temp13,temp9
      double precision num2,den2inv,dden2t
      double precision num3,num3inv,den3inv,dnum3r,dden3r,dnum3t,dden3t
      double precision den4inv,dden4t
      double precision num5,den5inv,dden5t
      double precision num6,den6inv,dden6t
      double precision num7,den7inv,dnum7t,dden7t
      double precision num9,den9inv,dden9t
      double precision num10,den10inv,dden10t
      double precision num11,den11inv,dden11t
      double precision expdeninv,expfrac
      double precision tkm1,thetatm1

c..  initial isotopic number abondances (number densities)
      rhompinv = rok/mprot
      do j = 1,nspecies
         xmol(j) = 0.d0
         xmol0(j) = 0.d0
      enddo  
      xmol(iih) = (x(ih1)/anuc(ih1)+x(ih2)/anuc(ih2))*rhompinv
      xmol(iic) = (x(ic12)/anuc(ic12)+x(ic13)/anuc(ic13)+
     &     x(ic14)/anuc(ic14))*rhompinv
      xmol(iin) = (x(in13)/anuc(in13)+x(in14)/anuc(in14)+
     &     x(in15)/anuc(in15))*rhompinv
      xmol(iio) = (x(io15)/anuc(io15)+x(io16)/anuc(io16)+
     &     x(io17)/anuc(io17)+x(io18)/anuc(io18))*rhompinv

C --> For test
c      RapCO = 0.99d0
c      xmol(iic) = RapCO*xmol(iio)
c      xmol(iin) = xmol(iic)*0.4d0
C <--

      xmol0(iih) = xmol(iih)
      xmol0(iic) = xmol(iic)
      xmol0(iin) = xmol(iin)
      xmol0(iio) = xmol(iio)

C need for opacity calculation
      xxh = x(ih1)+x(ih2)
      xxhe = x(ihe3)+x(ihe4)
      xxz = 1-xxh-xxhe
           
*-------------------------------------------------------------
***   calculation of the dissociation constants K'(T) and K(T)
*-------------------------------------------------------------

      theta = 5040.d0/tk
      logkT = log10(boltz*tk)

      do i = 1,natoms
         logKp(i) = 0.d0
      enddo
      do i = 1,neqchim
         logKp(i+natoms) = 0.d0
         do j = 1,5
            logKp(i+natoms) = logKp(i+natoms)+akap(i,j)*theta**(j-1)
         enddo
         if (ireac(i,1).eq.1) then
            tempnbr = ireac(i,3)-1.d0
         elseif (ireac(i,1).eq.2) then 
            tempnbr = ireac(i,3)+ireac(i,5)-1.d0
         endif
         logK(i+natoms) = logKp(i+natoms)-tempnbr*logkT
      enddo

c..   sort logK by increasing order of magnitude
      do i = 1,neqchim
         ij = natoms+1
         do j = natoms+1,natoms+neqchim
            if (logK(j).lt.logK(i+natoms)) ij = ij+1
         enddo
         order(ij-natoms) = i+natoms
         kdisso(i+natoms) = 10.d0**logK(i+natoms)
      enddo     

*----------------------------------------------------------------
***   compute molecular and atomic number densities : xmol (cm-3)
*----------------------------------------------------------------

C find a first solution
      do j = 1,neqchim
         jcur = order(j)
         kdiss = kdisso(jcur)
         jmol = jcur-natoms

C Molecules type X2
         if (ireac(jmol,1).eq.1) then 
            if (ireac(jmol,3).ne.2) stop 'Pb in opa_co 1'
            im = ireac(jmol,2)
            imol = ireac(jmol,4)
            temp = dsqrt(1.d0+8.d0*xmol(im)/kdiss)-1.d0
            xmol(im) = 0.25d0*kdiss*temp
            xmol(imol) = 0.25d0*xmol(im)*temp
C Molecules type XY
         elseif (ireac(jmol,1).eq.2) then 
c     find the most abundant species 
            if (xmol(ireac(jmol,2)).gt.xmol(ireac(jmol,4))) then
               ip = ireac(jmol,2)
               im = ireac(jmol,4)
            else
               ip = ireac(jmol,4)
               im = ireac(jmol,2)
            endif
            imol = ireac(jmol,6)
            if (ireac(jmol,3)+ireac(jmol,5).eq.2) then !molecule type XY
               temp = kdiss+xmol(ip)-xmol(im)
               nel0im = xmol(im)
               xmol(im) = 0.5d0*temp*(dsqrt(1.d0+4.d0*kdiss*xmol(im)/
     &              temp**2)-1.d0)
               xmol(imol) = abs(nel0im-xmol(im))
               xmol(ip) = xmol(ip)-xmol(imol)
            else                !molecule type H20
               if (ireac(jmol,3).eq.2) then
                  ip = ireac(jmol,2)
                  im = ireac(jmol,4)
               else
                  ip = ireac(jmol,4)
                  im = ireac(jmol,2)
               endif
                  imol = ireac(jmol,6)
               if (xmol(im).lt.xmol(ip)) then ! (nH>>n0)
                  xmol(im) = kdiss/(kdiss+xmol(ip)**2)*xmol(im)
                  xmol(imol) = xmol(ip)**2*xmol(im)/kdiss
                  xmol(ip) = xmol(ip)-2.d0*xmol(imol)
               else             ! (nH<<nO)
                  xmol(imol) = xmol(ip)
                  xmol(ip) = 0.d0
                  xmol(im) = xmol(im)-xmol(imol)
               endif
            endif
         else
            print *,jmol,j,nreac
            print *,ireac(jmol,1),ireac(jmol,2),ireac(jmol,3),
     &           ireac(jmol,4),ireac(jmol,5),ireac(jmol,6)
            stop 'Pb in opa_co 2'
         endif
      enddo

      logicNR = .true.
      debug = .false.
      neqo = neqchimmax+natomsmax
      n1 = natoms+neqchim

c.. Start newtom-raphson iterations

      if (logicNR) then
         tolf = 1.d-12
         tolx = 1.d-3
         do k = 1,100

C initialisation
            fvec(1:neqchimmax) = 0.d0
            fjac(1:neqchimmax,1:neqchimmax) = 0.d0


C --> Compute vectors and jacobien matrix
c     atoms
            do i = 1,natoms
               fvec(i) = xmol(i)-xmol0(i)
               fjac(i,i) = 1.d0
            enddo
            do i = 1,neqchim
               ir1 = ireac(i,2)
               nr1 = ireac(i,3)
c     monoatomic molecules
               if (ireac(i,1).eq.1) then
                  ip1 = ireac(i,4)

                  fvec(ir1) = fvec(ir1)+nr1*xmol(ip1)
                  fvec(ip1) = xmol(ir1)**nr1-kdisso(ip1)*xmol(ip1)

                  fjac(ir1,ip1) = fjac(ir1,ip1)+nr1
                  if (fjac(ip1,ip1).ne.0.d0) stop 'PB in usrfun 1'
                  fjac(ip1,ip1) = -kdisso(ip1)
                  fjac(ip1,ir1) = nr1*xmol(ir1)**(nr1-1)
c     diatomic molecules
               elseif (ireac(i,1).eq.2) then
                  ir1 = ireac(i,2)
                  nr1 = ireac(i,3)
                  ir2 = ireac(i,4)
                  nr2 = ireac(i,5)
                  ip1 = ireac(i,6)

                  fvec(ir1) = fvec(ir1)+nr1*xmol(ip1)
                  fvec(ir2) = fvec(ir2)+nr2*xmol(ip1)
                  fvec(ip1) = xmol(ir1)**nr1*xmol(ir2)**nr2-kdisso(ip1)*
     &                 xmol(ip1)

                  fjac(ir1,ip1) = fjac(ir1,ip1)+nr1
                  fjac(ir2,ip1) = fjac(ir2,ip1)+nr2
                  if (fjac(ip1,ip1).ne.0.d0) then
                     print *,ip1,fjac(ip1,ip1)
                     stop 'PB in usrfun 2'
                  endif
                  fjac(ip1,ip1) = -kdisso(ip1)
                  if (nr1.eq.1) then 
                     fjac(ip1,ir1) = xmol(ir2)**nr2
                  else
                     fjac(ip1,ir1) = nr1*xmol(ir1)**(nr1-1)*
     &                    xmol(ir2)**nr2
                  endif
                  if (nr2.eq.1) then 
                     fjac(ip1,ir2) = xmol(ir1)**nr1
                  else
                     fjac(ip1,ir2) = nr2*xmol(ir2)**(nr2-1)*
     &                    xmol(ir1)**nr1
                  endif
c     no triatomic molecules
               else
                  stop 'PB in usrfun 3'
               endif
            enddo

            if (debug) then
               write (nout,*) '/****  VECT  ****/'
               write (nout,100) (fvec(j),j=1,neqchim+natoms)
               write (nout,*) '/****  JACO  ****/'
               do i = 1,neqchim+natoms
                  write (nout,100) (fjac(i,j),j=1,neqchim+natoms)
 100              format (30(1x,1pe8.1))
               enddo
            endif

            errf = 0.d0
            do i = 1,nspecies     !Check function convergence.
               errf = errf+abs(fvec(i))
            enddo
            if(errf.le.tolf) then
               print *,'1',k
               goto 10
            endif
            do i = 1,natoms+neqchim !Right-hand side of linear equations
               p(i) = -fvec(i) 
            enddo

c --> Solve linear equations
c.. leqs
c            call leqs(fjac,p,n1,neqo)
c.. lapack
            call dgetrf (neqo,natoms+neqchim,fjac,neqo,indx,info)
            call dgetrs ('N',n1,1,fjac,neqo,indx,p,n1,info)

            errx = 0.d0           !Check root convergence.
            do i = 1,nspecies     !Update solution.
               if (xmol(i).ne.0.d0) errx = errx+abs(p(i)/xmol(i))
               if (.not.(p(i).le.0.d0).and..not.(p(i).ge.0.d0)) then
                  print *,i,p(i),ksh,tk,rok
                  stop 'PB !!!!!!!'
               endif
               xmol(i) = xmol(i)+p(i)
            enddo
            if(errx.le.tolx) then
               goto 10
            endif
         enddo

      endif
C verification de la solution converg�e :
 10   if (k.ge.100) then
         print *,'Too many iteration in NR procedure: T =',tk,
     &        ' rho =',rok
         error = 84
         return
c         stop "N-R failed"
      endif
      do i = 1,nmole+4
         if (xmol(i).lt.0.d0) then
            if (xmol(i).gt.-1.d-5) then
               xmol(i) = 0.d0
            else
               write (nout,200) imol,ksh
 200           format ('error in abundances determination for ',i2,
     &              ', Abund < 0 at shell :',i4)
               write (nout,*) 'abundance :',xmol(i)
               error = 84
               return
c               stop "Abundance < 0 in opa_co"
            endif
         endif
      enddo

*-------------------------------
***  determine electron pressure
*-------------------------------

      nelec = 0.d0
      xnsp(1:nioniz) = 1.d-50
c...  compute abundances according to molecules formation
      do i = 2,nbz-1
         if (i.lt.6.or.i.gt.8) then
            do j = 1,nis
               if (znuc(j).eq.i) xnsp(i) = xnsp(i)+x(j)/anuc(j)
            enddo
            xnsp(i) = xnsp(i)*rhompinv
         endif
      enddo
      xnsp(1) = xmol(1) !H
      xnsp(6) = xmol(2) !C
      xnsp(7) = xmol(3) !N
      xnsp(8) = xmol(4) !O
c..   scale from Grevesse 95 for elements > Cl
      zscale = zkint*rhompinv/zsol
      do k = 1,7
         xnsp(k+17) = xspsol(nis-1+k)/aheavy(k)*zscale
      enddo

c...  compute partition functions
      lntk = log(tk)
      logT = log10(5040.d0/tk)
      betaev = -1.d0/(tk*8.61173d-5)

      ij = 0
      do j = 1,nioniz
c..   Irwin table only considers F,Ca,Ti,Cr,Mn,Fe,Ni (i.e. j=9,19..24)
         if (j.eq.9.or.j.ge.19) then
            ij = ij+1
            lnQ(1) = aIrwin(ij,1)
            lnQ(2) = aIrwin(ij,7)
            do i = 2,6
               lnQ(1) = lnQ(1)+aIrwin(ij,i)*lntk**(i-1)
               lnQ(2) = lnQ(2)+aIrwin(ij,i+6)*lntk**(i-1)
            enddo
            expQ = exp(lnQ(1)-lnQ(2))

         else
            lnQ(1) = aSauval(j,1)
            lnQ(2) = aSauval(j,6)
            do i = 2,5
               lnQ(1) = lnQ(1)+aSauval(j,i)*logT**(i-1)
               lnQ(2) = lnQ(2)+aSauval(j,i+5)*logT**(i-1)
            enddo
            expQ = 10.d0**(lnQ(1)-lnQ(2))

         endif
         ntot(j) = xnsp(j)
         Kion(j) = 2.d0*expQ*saha*tk**1.5d0*exp(Eion(j)*betaev)
      enddo

      delta = 1.d-3
      nelecb = rhompinv*muiinvk*10**(-1.d1+tk*1.d-3)
      iterelec = 0
      y = 2.d0*delta
      do while (abs(y).gt.delta)
         iterelec = iterelec+1
         nelec = nelecb
         y = nelec
         a = 1
         do j = 1,nioniz
            temp = 1.d0/(nelec+Kion(j))
            tempbis = Kion(j)*ntot(j)*temp
            y = y-tempbis
            a = a+tempbis*temp
         enddo
         b = y-a*nelec
         nelecb = -b/a
c         if (iterelec.ge.1d2) print *,iterelec,y,delta,abs(y).gt.delta
         if (iterelec.eq.1d2) goto 20
      enddo

 20   pel = nelec*((1.d0-nelec/(muiinvk*rhompinv)))*boltz*tk

      nh = xmol(iih)
c      nc = xmol(iic)
c      nn = xmol(iin)
c      no = xmol(iio)
      nh2 = xmol(iih2)
      nh2o = xmol(iih2o)
      noh = xmol(iioh)
      nco = xmol(iico)
      ncn = xmol(iicn)
      nc2 = xmol(iic2)
c      nn2 = xmol(iin2)
c      nch = xmol(iich)
c      nno = xmol(iino)
c      no2 = xmol(iio2)
c      nnh = xmol(iinh)

Ccc... Verifications
c      write (111,11) ksh,tk,logK(natoms+1),logK(natoms+2),logK(natoms+3)
c     &     ,logK(natoms+4),logK(natoms+5),logK(natoms+6),logK(natoms+7),
c     &     nh,nc,nn,no,nh2,nh2o,noh,nco,ncn,nc2,rok,nn2,nch,nno,no2,nnh,
c     &     xmolhplus,nelec

c 11   format (1x,i4,40(1x,2pe11.4))


*------------------
***   compute kappa
*------------------

      nhtot = nh+2*nh2+2*nh2o+noh

      tt = tk*1.d-4
      tm1 = 1.d0/tt
      t2 = tt**2
      tm2 = tm1**2
      t3 = tt*t2
      t4 = t2**2
      tm4 = tm2**2
      t5 = t2*t3
      tm5 = tm1*tm4
      t6 = t3**2
      tm6 = tm1*tm5
      t7 = t3*t4
      tm7 = tm1*tm6
      t9 = t7*t2
      t10 = tt*t9
      t05 = dsqrt(tt)

      rm1 = 1.d0/rok
      r05 = dsqrt(rok)
      r025 = dsqrt(r05)
      rm075 = rm1**0.75d0

      kap1 = 0.d0
      kapr1 = 0.d0
      kapt1 = 0.d0


c Keeley 70
c electrons
      temp1 = 5.4d-13*rm1*tm1
      temp1t = -temp1*tm1
      temp1r = -temp1*rm1

      num2 = t05
      den2inv = 1.d0/(2.d6*tm4+2.1d0*t6)
      dden2t = -8.d6*tm5+12.6d0*t5
      temp2 = num2*den2inv
      temp2t= temp2*(0.5d0*tm1-dden2t*den2inv)

      temp3a = 4.d-3*r025+2.d-4*t4
      num3 = temp3a
      num3inv = 1.d0/num3
      den3inv = 1.d0/(4.5d0*t6*temp3a+r025*t3)
      dnum3r = 1.d-3*rm075
      dden3r = rm075*t3*(0.25d0+4.5d-3*t3)
      dnum3t = 8.d-4*t3
      dden3t = t2*(r025*(3.d0+0.108d0*t3)+9.d-3*t7)
      temp3 = (1.d0-2.d0*nh2/nhtot)*num3*den3inv
      temp3r = temp3*(dnum3r*num3inv-dden3r*den3inv)
      temp3t = temp3*(dnum3t*num3inv-dden3t*den3inv)

      den4inv = 1.d0/(1.4d3*tt+t6)
      dden4t = 1.4d3+6.d0*t5
      temp4 = den4inv
      temp4t = -temp4*dden4t*den4inv

      num5 = 1.5d0
      den5inv = 1.d0/(1.d6+0.1d0*t6)
      dden5t = 0.6d0*t5
      temp5 = num5*den5inv
      temp5t = -temp5*dden5t*den5inv

      num6 = t05
      den6inv = 1.d0/(20.d0*tt+5.d0*t4+t5)
      dden6t = 20.d0+20.d0*t3+5.d0*t4
      temp6 = num6*den6inv
      temp6t = temp6*(0.5d0*tm1-dden6t*den6inv)

      kap1 = pel*(temp1+xxh*(temp2+temp3)+xxhe*(temp4+temp5)+xxz*temp6)
      kapr1 = pel*(temp1r+xxh*temp3r)
      kapt1 = pel*(temp1t+xxh*(temp2t+temp3t)+xxhe*(temp4t+temp5t)+
     &     xxz*temp6t)

c H et H2
      num7 = 5.55d-27*t4
      den7inv = 1.d0/(1.d0+10.d0*t6+3.42d-5*tm6)
      dnum7t = 2.22d-26*t3
      dden7t = 60.d0*t5-2.052d-4*tm7
      temp7 = (nh+nh2)*rm1*num7*den7inv
      temp7r = -temp7*rm1
      temp7t = temp7*(dnum7t/num7-dden7t*den7inv)

c CO
      temp8 = 2.75d-26*nco*rm1
      temp8r = -temp8*rm1

c OH
      num9 = 1.4d-21*t6
      den9inv = 1.d0/(0.1d0+t6)
      dden9t = 6.d0*t5
      temp9 = noh*rm1*num9*den9inv
      temp9r = -temp9*rm1
      temp9t = temp9*(6.d0*tm1-dden9t*den9inv)

c H20
c Attention: changement du facteur temp11 suivant Marigo (2002)
      num10 = 2.6d-27
      den10inv = 1.d0/(4.23d-4+t4)
      dden10t = 4.d0*t3
      temp10 = num10*den10inv
      temp10t = -temp10*dden10t*den10inv

      expdeninv = 1.d0/(tt+0.37d0)
      expfrac = 3.2553d0*expdeninv
      num11 = 9.72d-21*exp(-expfrac)
      den11inv = 1.d0/(1.d0+3.78d3*t10)
      dden11t = 3.78d4*t9
      temp11 = num11*den11inv
      temp11t = temp11*(expfrac*expdeninv-dden11t*den11inv)

      temp12a = nh2o*rm1
      temp12 = temp12a*(temp10+temp11)
      temp12r = -temp12*rm1
      temp12t = temp12a*(temp10t+temp11t)

      kap1 = kap1+temp7+temp8+temp9+temp12
      kapr1 = kapr1+temp7r+temp8r+temp9r+temp12r
      kapt1 = 1.d-4*(kapt1+temp7t+temp9t+temp12t)

c Scalo & Ulrich 75 (CN) and  Querci et al 71 (C2)
      tkm1 = 1.d0/tk
      theta = 5040.d0*tkm1
      thetatm1 = theta*tkm1
      temp13a = -19.212d0+2.2479d0*theta-2.8069d0*theta**2+0.76d0*
     &     theta**3-0.078384d0*theta**4
      temp13b = -2.2479d0+2.d0*2.8069d0*theta-3.d0*
     &     0.76d0*theta**2+4.d0*0.078384d0*theta**3

      temp13 = 10.d0**temp13a*(ncn+nc2)*rm1
      temp13r = -temp13*rm1

      kap1 = kap1+temp13
      kapr1 = kapr1+temp13r
      kapt1 = kapt1+log(10.d0)*temp13*temp13b*thetatm1

      return
      end
