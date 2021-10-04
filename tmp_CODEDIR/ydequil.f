
************************************************************************

      SUBROUTINE ydequil

************************************************************************
* Calculate the accretion equilibrium abundance of H1, H2 and He3      *
* * 1  :    2 proto ( 0 gamma, 0 nutri)  1 deutr   -->      used   +   *
* * 2  :    1 deutr ( 1 proto, 0 gamma)  1 he  3   -->      used   -   *
* * 3  :    2 deutr ( 0 gamma, 0 gamma)  1 alpha   --> not  used   -   *
* * 4  :    2 deutr ( 0 gamma, 1 neutr)  1 he  3   --> not  used   -   *
* * 5  :    1 he  3 ( 1 deutr, 1 proto)  1 alpha   --> not  used   -   *
* * 7  :    1 alpha ( 1 deutr, 0 gamma)  1 li  6   --> not  used   -   *
* * 13 :    1 li  7 ( 1 deutr, 1 neutr)  2 alpha   --> not  used   -   *
* * 16 :    1 be  7 ( 1 deutr, 1 proto)  2 alpha   --> not  used   -   *
* * 19 :    1 be  9 ( 1 proto, 1 deutr)  2 alpha   -->      used   +   *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.nuc'
      include 'evolcom.conv'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      integer mmean,ishell,iend,istart,nconv1,kl0
      integer i,j,k,l,kl,ideut,ivit,ndeut
      integer normx,imin,imax

      double precision sum1
      double precision exv
      double precision tdest,vdest,vdmean,prod,sprod
c     double precision fprod
      double precision xprot,t9,fscrd,fscrbe,xbe,taumean,dtdtaud,xdd,
     &     xdeq,xheq,xpeq,sdest,tauac,taud,xdconv,xdeut,xdburn,dxm,v
c    &     aa0,bb0,edx,edf,slope,v
      double precision dt9l,weight,tnuc9,vrate

      logical partialmix

      parameter (ndeut = 3)

      common /nucvit/ tnuc9(nvit),vrate(nvit,nre)

      dimension tdest(nsh),prod(nsh),taumean(nmaxconv),xdeq(nmaxconv),
     &     xheq(nmaxconv),xpeq(nmaxconv),xdd(nmaxconv),v(ndeut),
     &     ideut(ndeut)

      external exv

      data (ideut(k),k = 1,ndeut)/1,2,19/

      if (nsconv.eq.0) return

*______________________
***   convective region
*----------------------

      iprint(1) = 0
      tauac = 1.d99
      if (nsconv.ge.1) then
c..   restore initial abundances
         do l = 1,nsp
            do k = 1,nmod
               xsp(k,l) = vxsp(k,l)
               ysp(k,l) = xsp(k,l)/anuc(l)
            enddo
         enddo
         do kl = 1,nsconv
            mmin(kl) = novlim(kl,3)
            mmax(kl) = novlim(kl,4)
            tauaccr = 1.d99
            taud = 1.d99
            taumean(kl) = 1.d99
            itauac = -1
            mmean = int((mmin(kl)+mmax(kl))*0.5)
            vdmean = 0.d0
c           fprod = 0.d0
            sprod = 0.d0
            if (itacc.eq.2) then
               idcur = nnucacc(kl)
            else
               idcur = mmin(kl)
            endif
            imin = mmin(kl)
            imax = mmax(kl)
            call mix (dm,imin,imax,kl,1,partialmix)
            nnucacc(kl) = mmin(kl)-1
            if (t(imin).ge.tnucmin) then
               xprot = ysp(mmean,ih1)
               xbe = ysp(mmean,10)
               do 30 k = mmin(kl),mmax(kl)
                  if (t(k).gt.tnucmin) then

*________________________________________________
***   interpolation of the nuclear reaction rates
*------------------------------------------------

                     t9 = t(k)*1.d-9
                     do i = 1,nvit
                        if (t9.lt.tnuc9(i)) goto 10
                     enddo
 10                  ivit = i-1
                     if (ivit.eq.0) then
                        do j = 1,ndeut
                           v(j) = vrate(1,ideut(j))
                           v(j) = 10.d0**v(j)
                        enddo
                        goto 20
                     endif
                     if (ivit.gt.(nvit-1)) then
                        do j = 1,ndeut
                           v(j) = vrate(nvit,ideut(j))
                           v(j) = 10.d0**v(j)
                        enddo
                        goto 20
                     endif
                     dt9l = log10(t9/tnuc9(ivit))
                     weight = dt9l/(log10(tnuc9(ivit+1)/tnuc9(ivit)))
                     do j = 1,ndeut
                        v(j) = vrate(ivit,ideut(j))+(vrate(ivit+1,
     &                       ideut(j))-vrate(ivit,ideut(j)))*weight
                        v(j) = 10.d0**v(j)
                     enddo

 20                  fscrd = 1.d0
                     vdest = v(2)*ro(k)*fscrd
                     tdest(k) = 1.d0/(vdest*xprot)
                     fscrbe = 1.d0
                     prod(k) = (v(3)*xbe*fscrbe+v(1)*xprot*fscrd)*ro(k)
                     if (tdest(k).le.5.d0*tconv(k).and.itacc.lt.3)
     &                    then
                        nnucacc(kl) = k
                     else
                        sprod = sprod+prod(k)*dm(k)
                        vdmean = vdmean+vdest*dm(k)
                     endif
                  endif
 30            continue

               if (itacc.ne.1.or.itacc.ne.2.or.nnucacc(kl).lt.mmin(kl))
     &              then
                  am(kl) = 1.d0/(m(mmax(kl)+1)-m(mmin(kl)))
                  amacc(kl) = macc(mmax(kl)+1)-macc(mmin(kl))
                  if (amacc(kl).gt.0.d0) tauaccr = dtn/(am(kl)*
     &                 amacc(kl))
                  taumean(kl) = tauaccr
               else
                  am(kl) = 1.d0/(m(mmax(kl)+1)-m(nnucacc(kl)+1))
                  amacc(kl) = macc(mmax(kl)+1)-macc(nnucacc(kl)+1)
                  if (amacc(kl).gt.0.d0) tauaccr = dtn/(am(kl)*
     &                 amacc(kl))
                  if (macc(nnucacc(kl)+1).gt.0.d0) taumean(kl) = dtn*
     &                 (m(nnucacc(kl)+1)-m(mmin(kl)))/(macc(nnucacc(kl)+
     &                 1)-macc(mmin(kl)))
               endif
               if (nnucacc(kl).lt.mmax(kl)) then
                  vdmean = vdmean/(m(mmax(kl)+1)-m(nnucacc(kl)+1))
                  if (vdmean.gt.0.d0) taud = 1.d0/(vdmean*xprot)
                  taudnuc = min(taud,taudnuc)
                  dtdtaud = dtn/taud
c                 fprod = sprod/(m(mmax(kl)+1)-m(mmin(kl)))
                  xdd(kl) = vyd1(mmean)*(1.d0-am(kl)*amacc(kl))
                  if (dtdtaud.gt.1.d-9) then
                     xdeq(kl) = xdd(kl)*exv(-dtdtaud)+xspacc(ih2)*
     &                    (1.d0-exv(-dtdtaud))*taud/tauaccr
                  else
                     xdeq(kl) = xdd(kl)+xspacc(ih2)*dtn/tauaccr
                  endif
c                 xdeq(kl) = xprot*taud*fprod*anuc(ih2)+xdd(kl)*
c    &                 exv(-dtdtaud)+xspacc(ih2)*(1.d0-exv(-dtdtaud))*
c    &                 taud/tauaccr
                  xdburn = (xdd(kl)+xspacc(ih2)*taud/tauaccr)*(1.d0-
     &                 exv(-dtdtaud))
                  xpeq(kl) = vxsp(mmean,ih1)-xdburn*anuc(ih1)/
     &                 anuc(ih2)
                  xheq(kl) = vxsp(mmean,ihe3)+xdburn*anuc(ihe3)/
     &                 anuc(ih2)
                  if (itacc.eq.1) xdd(kl) = vyd1(mmean)*(1.d0-dtn/
     &                 taumean(kl))
                  if (nnucacc(kl).le.mmin(kl).or.itacc.ge.3) then
                     do k = mmin(kl),mmax(kl)
                        xsp(k,ih1) = xpeq(kl)
                        xsp(k,ih2) = xdeq(kl)
                        xsp(k,ihe3) = xheq(kl)
                     enddo
                  endif
               endif
            endif
         enddo
      endif

      ldacc = massrate*msun*qdeut*xspacc(ih2)/(anuc(ih2)*sec)

*_____________________
*     radiative region
*---------------------

      ishell = nsconv
      if (crz(1).gt.0.or.nnucacc(1).gt.0) ishell = nsconv+1
      do kl = 1,ishell
         iend = nmod
         if (crz(1).gt.0.or.nnucacc(1).gt.0) then
            if (kl.eq.1) then
               istart = 1
               kl0 = 1
               iend = nnucacc(kl0)
            else
               istart = mmax(kl-1)+1
               kl0 = kl
               if (kl.lt.ishell) iend = nnucacc(kl0)
            endif
         else
            istart = mmax(kl)+1
            kl0 = kl+1
            if (kl.lt.nsconv) iend = nnucacc(kl0)
         endif
         do k = istart,iend
            sdest = 0.d0
            sprod = 0.d0
            taud = 1.d99
            tauac = 1.d99
            if (t(k).gt.tnucmin) then
               if (crz(k).ge.-1.or.itacc.eq.2) then
                  if (facc(k).gt.0.d0) then
                     if (iaccr.eq.2.or.iaccr.eq.4) tauac = dtn*(facc(k)*
     &                    sinthac+1.d0)/(facc(k)*sinthac)
                     if (iaccr.eq.1.or.iaccr.eq.3) tauac = dtn/(facc(k)*
     &                    sinthac)
                  endif
                  if (tauac.lt.tauaccr) then
                     itauac = k
                     tauaccr = tauac
                  endif
               endif
               if (crz(k).ge.-1) then
                  xdeut = yd1(k)

*________________________________________________
***   interpolation of the nuclear reaction rates
*------------------------------------------------

                  t9 = t(k)*1.d-9
                  do i = 1,nvit
                     if (t9.lt.tnuc9(i)) goto 40
                  enddo
 40               ivit = i-1
                  if (ivit.eq.0) then
                     do j = 1,ndeut
                        v(j) = vrate(1,ideut(j))
                        v(j) = 10.d0**v(j)
                     enddo
                     goto 50
                  endif
                  if (ivit.gt.(nvit-1)) then
                     do j = 1,ndeut
                        v(j) = vrate(nvit,ideut(j))
                        v(j) = 10.d0**v(j)
                     enddo
                     goto 50
                  endif
                  dt9l = log10(t9/tnuc9(ivit))
                  weight = dt9l/(log10(tnuc9(ivit+1)/tnuc9(ivit)))
                  do j = 1,ndeut
                     v(j) = vrate(ivit,ideut(j))+(vrate(ivit+1,
     &                    ideut(j))-vrate(ivit,ideut(j)))*weight
                     v(j) = 10.d0**v(j)
                  enddo

 50               fscrd = 1.d0
                  sdest = v(2)*ro(k)*ysp(k,ih1)*fscrd
                  tdest(k) = 1.d0/sdest
                  taudnuc = min(taud,taudnuc)
                  fscrbe = 1.d0
                  xbe = ysp(k,10)
                  sprod = (v(3)*xbe*fscrbe+v(1)*ysp(k,ih1)*fscrd)*ro(k)
               else
                  if (itacc.le.1) then
                     tauac = taumean(kl0)
                     xdeut = xdd(kl0)
                  endif
                  if (itacc.eq.2) xdeut = yd1(k)
                  sprod = prod(k)
               endif
               taud = tdest(k)
               dtdtaud = dtn/taud
               xsp(k,ih2) = xspacc(ih2)*(1.d0-exv(-dtdtaud))*taud/tauac+
     &              xdeut*exv(-dtdtaud)
c    &              xdeut*exv(-dtdtaud)+xsp(k,ih1)*taud*sprod*anuc(ih2)/
c    &              anuc(ih1)
               xdburn = (xdeut+xspacc(ih2)*taud/tauac)*(1.d0-
     &              exv(-dtdtaud))
               xsp(k,ihe3) = vxsp(mmean,ihe3)+xdburn*anuc(ihe3)/
     &              anuc(ih2)
               xsp(k,ih1) = vxsp(mmean,ih1)-xdburn*anuc(ih1)/
     &              anuc(ih2)
            else
               tdest(k) = 1.d50
            endif
         enddo
      enddo

*__________________________
***   linear interpolation
*__________________________

      do 70 kl = 1,nsconv
         if (nnucacc(kl).gt.novlim(kl,3)) rnucacc = r(nnucacc(kl))/
     &        r(nmod)
         if (nnucacc(kl).eq.0.or.itacc.eq.3.or.(nnucacc(kl).lt.
     &        mmin(kl).and.itacc.lt.3).or.t(nnucacc(kl)).lt.
     &        tnucmin) goto 70
         if (facc(mmax(kl)).gt.0.d0.and.tdest(nnucacc(kl)+1).lt.
     &        1.d2*dtn) then
            nconv1 = min(nnucacc(kl)+1,mmax(kl))
            if ((nconv1-nnucacc(kl)).lt.2) goto 70
            dxm = m(nconv1)-m(nnucacc(kl))
            xdconv = (xdeq(kl)*(m(mmax(kl))-m(nnucacc(kl)))-0.5d0*
     &           dxm*xsp(nnucacc(kl),ih2))/(m(mmax(kl))-m(nconv1)+
     &           0.5d0*dxm)
            do k = nnucacc(kl),nconv1
               xsp(k,ih2) = (xdconv-xsp(nnucacc(kl),ih2))*(m(k)-
     &              m(nnucacc(kl)))/dxm+xsp(nnucacc(kl),ih2)
               xsp(k,ihe3) = vxsp(k,ihe3)+anuc(ihe3)/anuc(ih2)*
     &              (vxsp(k,ih2)-xsp(k,ih2))
               xsp(k,ih1) = vxsp(k,ih1)-anuc(ih1)/anuc(ih2)*
     &              (vxsp(k,ih2)-xsp(k,ih2))
            enddo
            if (mmax(kl).gt.nconv1) then
               do k = nconv1+1,mmax(kl)
                  xsp(k,ih2) = xsp(nconv1,ih2)
                  xsp(k,ihe3) = xsp(nconv1,ihe3)
                  xsp(k,ih1) = xsp(nconv1,ih1)
               enddo
            endif
         endif
 70   continue

      if (iter.eq.1) taudnuc0 = taudnuc

*______________________________
***   abundance renormalization
*------------------------------

      do k = 1,nmod
         normx = 1
         sum1 = 0.d0
         do j = 1,nsp
            if (xsp(k,j).gt.xsp(k,normx)) normx = j
            sum1 = sum1+xsp(k,j)
         enddo
         xsp(k,normx) = xsp(k,normx)+(1.d0-sum1)
      enddo

      return
      end
