
************************************************************************

      SUBROUTINE fdiff (eqd,cd,x,f1,f2,f3,id,imin,imax)

************************************************************************
*     compute the diffusion matrix elements and derivatives
*     eq[i,neqd1] = dYi = dt*Fi(Yj)  (eq of nucleosynthesis at shell id)
*                 = dt * d (D dYi/dm)/dm
*     eq[i,j=1,neqd]        = [dFi/dYj]_{id-1} (jacobian)
*     eq[i,j=neqd+1,2neqd]  = [dFi/dYj]_{id}   (jacobian)
*     eq[i,j=2neqd+1,3neqd] = [dFi/dYj]_{id+1} (jacobian)
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.diff'
      include 'evolcom.conv'
      include 'evolcom.teq'

      integer j,id,neqd,neqd1,neqd2,imin,imax

      double precision x
      double precision cd,f1,f2,f3,eqd
      double precision xbcz,x2,x3

      dimension cd(nsh),f1(nsh),f2(nsh),f3(nsh),eqd(nsp,3*nsp+1),
     &     x(nsh,nsp)

      neqd = nsp
      neqd2 = 2*neqd
      neqd1 = 3*neqd+1
      kcentd = 1
      ksurfd = 1

c..   start filling matrix elements at j=2 because neutrons don't diffuse
***   interior
      if (id.gt.imin.and.id.lt.imax) then
         x3 = f3(id)*cd(id)
         x2 = f2(id)*cd(id+1)
         do j = 2,neqd
            eqd(j,neqd1) = eqd(j,neqd1)+x2*(x(id+1,j)-x(id,j))-
     &           x3*(x(id,j)-x(id-1,j))
            eqd(j,j) = eqd(j,j)+x3
            eqd(j,neqd+j) = -x2-x3+eqd(j,neqd+j)
            eqd(j,neqd2+j) = eqd(j,neqd2+j)+x2
         enddo
         return
      endif


***   inner boundary condition
      if (id.eq.imin) then
c..   center or continuity equation
         if (kcentd.eq.1.or.id.eq.1) then
            x2 = f2(id)*cd(id+1)
            do j = 2,neqd
               eqd(j,neqd1) = eqd(j,neqd1)+x2*(x(id+1,j)-x(id,j))
               eqd(j,j) = eqd(j,j)-x2
               eqd(j,neqd+j) = eqd(j,neqd+j)+x2
            enddo
            return
         endif
c..   reservoir
         if (kcentd.eq.2) then
            do j = 1,neqd
               xbcz = cd(id+1)*yconvlc(j)*f1(id)/delmc
               eqd(j,neqd1) = eqd(j,neqd1)-xbcz*(x(id+1,j)/x(id,j)-1.d0)
               eqd(j,j) = eqd(j,j)+xbcz*x(id+1,j)/(x(id,j)*x(id,j))
               eqd(j,neqd+j) = eqd(j,neqd+j)-xbcz/x(id,j)
            enddo
         endif
c..   Y = cste
c         if (kcentd.eq.3) then
c            do j = 2,neqd
c               eqd(j,neqd1) = -x(id,j)+x(id-1,j)
c               eqd(j,j) = -1.d0
c            enddo
c         endif
      endif

***   outer boundary condition
      if (id.eq.imax) then
c..   surface or continuity equation
         if (ksurfd.eq.1.or.id.eq.nmod1) then
            x3 = f3(id)*cd(id)
            do j = 2,neqd
               eqd(j,neqd1) = eqd(j,neqd1)-x3*(x(id,j)-x(id-1,j))
               eqd(j,j) = eqd(j,j)+x3
               eqd(j,neqd+j) = eqd(j,neqd+j)-x3
            enddo
            return
         endif
c..   reservoir
         if (ksurfd.eq.2) then
            do j = 2,neqd
               xbcz = cd(id)*yconvls(j)*f1(id)/delms
               eqd(j,neqd1) = eqd(j,neqd1)-xbcz*(1.d0-x(id-1,j)/x(id,j))
               eqd(j,j) = eqd(j,j)+xbcz/x(id,j)
               eqd(j,neqd+j) = eqd(j,neqd+j)-xbcz*x(id-1,j)/
     &              (x(id,j)*x(id,j))
            enddo
         endif
c..   Y = cste
c          if (ksurfd.eq.3) then
c             do j = 2,neqd
c                eqd(j,neqd1) = -x(id,j)+x(id+1,j)
c                eqd(j,neqd+j) = -1.d0
c             enddo
c          endif
      endif


      return
      end
