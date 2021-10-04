

************************************************************************

      SUBROUTINE yneutron (y,v,yeq)

************************************************************************
*     Calculate the neutron equilibrium abundance yeq                  *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.nuc'

      integer mm
      integer ii,ij,j

      double precision y,v,yeq,zlocal
      double precision sprod,sdest

      dimension y(*),v(nreac)

***   no neutron captures for extremely metal poor stars (Z = 0)
      if (zkint.lt.5.d-6) then
         zlocal = 0.d0
         do j = 2,ib11
            zlocal = zlocal+y(j)*anuc(j)
         enddo
         zlocal = 1.d0-zlocal
         if (zlocal.lt.1.d-10) then
            yeq = 1.d-50
            return
         endif
      endif

      sdest = 0.d0
!$OMP PARALLEL
!$OMP DO PRIVATE(j) REDUCTION(+:sdest)
      do j = 1,idneut-1
         mm = jneut(j)
         ii = k1(mm)
         sdest = sdest+v(mm)*y(ii)
      enddo
!$OMP END DO

      sprod = 0.5d0*v(iddn)*y(ih2)*y(ih2)
!$OMP DO PRIVATE(j) REDUCTION(+:sprod)
      do j = idneut+1,nneut
         mm = jneut(j)
         ii = k1(mm)
         ij = k3(mm)
         sprod = sprod+v(mm)*y(ii)*y(ij)
      enddo
!$OMP END DO
!$OMP END PARALLEL

      if (sdest.gt.0.d0) then
         yeq = sprod/sdest
      else
         yeq = y(1)
      endif
      yeq = max(abs(yeq),1.d-50)
      return
      end
