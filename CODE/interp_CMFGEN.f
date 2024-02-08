
      SUBROUTINE interp_CMFGEN(Mdot,gs,interp_T,interp_rho,
     &     interp_r,interp_ne,interp_v)
      
************************************************************************
*                                                                      *
* Interpolate within the CMFGEN atmosphere tables to produce the       *
* atmospheric profiles used in the computation of the internal         *
* structure for massive stars.                                         *
*                                                                      *
*                                                                      *
* References                                                           *
* ----------                                                           *
*                                                                      *
* [Shepard, 1968] Shepard, D. (1968). A two-dimensional interpolation  *
*     function for irregularly-spaced data. In Proceedings of the 1968 *
*     23rd ACM national conference, page 517–524, New York, NY, USA.   *
*     Association for Computing Machinery. DOI : 10.1145/800186.810616 *
*                                                                      *
************************************************************************
*                                                                      *
* $LastChangedDate:: 2023-06-06 16:39:00 +0000 (Tue, 06 Jun 2023)    $ *
* $Author:: VOJE                                                     $ *
* $Rev:: 110                                                         $ *
*                                                                      *
************************************************************************

      implicit none

************************************************************************
*                              INCLUDES                                *
************************************************************************

      include 'evolpar.star'
      include 'evolcom.atm'
      include 'evolcom.surf'
      include 'evolcom.teq'

************************************************************************
*                             DEFINITIONS                              *
************************************************************************

***   Inputs/outputs.
************************************************************************
      
***   Mdot : Value of log(Mdot) (mass loss) at the current time-step.
***   gs   : Value of log(gs) (effective gravity) at the current
***          time-step.
      double precision, intent(in) :: Mdot, gs

***   interp_T   : Interpolated T profile.
***   interp_rho : Interpolated rho profile.
***   interp_r   : Interpolated r profile.
***   interp_ne  : Interpolated ne profile.
***   interp_v   : Interpolated v profile.
      double precision, intent(out) :: interp_T(nsh), interp_rho(nsh),
     &     interp_r(nsh), interp_ne(nsh), interp_v(nsh)

***   Temporary variables.
************************************************************************
      
      integer :: i, j, k, l, n, o, nlines, id, stat
      double precision :: i_int, j_int, k_int
      double precision :: wn, sum_wn, wn_arr(100), used_CMFGEN(8,4),
     &     pp
      double precision :: tmp_tau_C(100), tmp_T_C(100), tmp_rho_C(100),
     &     tmp_r_C(100), tmp_ne_C(100), tmp_v_C(100)
      double precision :: DDT(100), DT(100), DDrho(100), Drho(100),
     &     DDr(100), Dr(100), DDne(100), Dne(100), DDv(100), Dv(100)
      double precision, allocatable :: tmp_T(:,:,:,:),
     &     tmp_rho(:,:,:,:), tmp_r(:,:,:,:), tmp_ne(:,:,:,:),
     &     tmp_v(:,:,:,:)

************************************************************************
*                       CMFGEN INTERPOLATION                           *
************************************************************************
      
***   1. Allocate arrays.
************************************************************************
      
      allocate(tmp_T(nb_Teff_CMFGEN,nb_Mdot_CMFGEN,nb_gs_CMFGEN,nsh))
      allocate(tmp_rho(nb_Teff_CMFGEN,nb_Mdot_CMFGEN,nb_gs_CMFGEN,nsh))
      allocate(tmp_r(nb_Teff_CMFGEN,nb_Mdot_CMFGEN,nb_gs_CMFGEN,nsh))
      allocate(tmp_ne(nb_Teff_CMFGEN,nb_Mdot_CMFGEN,nb_gs_CMFGEN,nsh))
      allocate(tmp_v(nb_Teff_CMFGEN,nb_Mdot_CMFGEN,nb_gs_CMFGEN,nsh))
      
***   2. Put everything on the same grid.
************************************************************************
      
      do id=1,max_id
         i = avail_CMFGEN(id,1)
         j = avail_CMFGEN(id,2)
         k = avail_CMFGEN(id,3)
         nlines = avail_CMFGEN(id,4)
         tmp_tau_C(:) = tau_CMFGEN(i,j,k,:)
         tmp_T_C(:) = T_CMFGEN(i,j,k,:)
         tmp_rho_C(:) = rho_CMFGEN(i,j,k,:)
         tmp_r_C(:) = r_CMFGEN(i,j,k,:)
         tmp_ne_C(:) = ne_CMFGEN(i,j,k,:)
         tmp_v_C(:) = v_CMFGEN(i,j,k,:)
         call splineatm (tmp_tau_C,tmp_T_C,nlines,1.d50,
     &        1.d50,DDT,DT)
         call splineatm (tmp_tau_C,tmp_rho_C,nlines,1.d50,
     &        1.d50,DDrho,Drho)
         call splineatm (tmp_tau_C,tmp_r_C,nlines,1.d50,
     &        1.d50,DDr,Dr)
         call splineatm (tmp_tau_C,tmp_ne_C,nlines,1.d50,
     &        1.d50,DDne,Dne)
         call splineatm (tmp_tau_C,tmp_v_C,nlines,1.d50,
     &        1.d50,DDv,Dv)
         l = nmod
         do while(tau(l).lt.100)
            n = nmod-l+1
            call splintatm (tmp_tau_C,
     &           tmp_T_C,DDT,nlines,
     &           tau(l),tmp_T(i,j,k,n))
            call splintatm (tmp_tau_C,
     &           tmp_rho_C,DDrho,nlines,
     &           tau(l),tmp_rho(i,j,k,n))
            call splintatm (tmp_tau_C,
     &           tmp_r_C,DDr,nlines,
     &           tau(l),tmp_r(i,j,k,n))
            call splintatm (tmp_tau_C,
     &           tmp_ne_C,DDne,nlines,
     &           tau(l),tmp_ne(i,j,k,n))
            call splintatm (tmp_tau_C,
     &           tmp_v_C,DDv,nlines,
     &           tau(l),tmp_v(i,j,k,n))
            l = l-1
         end do
      end do

***   3. Interpolate in (Teff, Mdot, gs) using an IDW interpolation
***   .. method [Shepard, 1968].
************************************************************************

      pp = 2.d0
      i_int = (teff-2.0d4)/1.d3
      j_int = -(Mdot+4.0)/1.d-1
      k_int = (gs-2.4d0)/5.d-2

***   3.1 Select the 8 nearest neighbours.
      wn_arr = 0.d0
      do id=1,max_id
         i = avail_CMFGEN(id,1)
         j = avail_CMFGEN(id,2)
         k = avail_CMFGEN(id,3)
         wn = ((i-i_int)*(i-i_int) + (j-j_int)*(j-j_int) +
     &        (k-k_int)*(k-k_int))**-pp/2.d0
         wn_arr(id) = wn
      end do
      call sort(100,wn_arr)
      do id=1,max_id
         i = avail_CMFGEN(id,1)
         j = avail_CMFGEN(id,2)
         k = avail_CMFGEN(id,3)
         wn = ((i-i_int)*(i-i_int) + (j-j_int)*(j-j_int) +
     &        (k-k_int)*(k-k_int))**-pp/2.d0
         do o=1,8
            if (ABS(wn-wn_arr(92+o))/wn.lt.1d-3) then
               used_CMFGEN(o,1) = i
               used_CMFGEN(o,2) = j
               used_CMFGEN(o,3) = k
               used_CMFGEN(o,4) = wn
            end if
         end do
      end do

***   3.2 Apply the IDW algorithm considering only the 8 nearest
***   ... neighbours.
      l = nmod
      do while (tau(l).lt.100)
         n = nmod-l+1
         sum_wn = 0.d0
         interp_T(l) = 0.d0
         interp_rho(l) = 0.d0
         interp_r(l) = 0.d0
         interp_ne(l) = 0.d0
         interp_v(l) = 0.d0
         do o=1,8
            i = used_CMFGEN(o,1)
            j = used_CMFGEN(o,2) 
            k = used_CMFGEN(o,3)
            wn = used_CMFGEN(o,4) 
            sum_wn = sum_wn + wn
            interp_T(l) = interp_T(l) + wn*tmp_T(i,j,k,n)
            interp_rho(l) = interp_rho(l) + wn*tmp_rho(i,j,k,n)
            interp_r(l) = interp_r(l) + wn*tmp_r(i,j,k,n)
            interp_ne(l) = interp_ne(l) + wn*tmp_ne(i,j,k,n)
            interp_v(l) = interp_v(l) + wn*tmp_v(i,j,k,n)
         end do
         interp_T(l) = interp_T(l)/sum_wn
         interp_rho(l) = interp_rho(l)/sum_wn
         interp_r(l) = interp_r(l)/sum_wn
         interp_ne(l) = interp_ne(l)/sum_wn
         interp_v(l) = interp_v(l)/sum_wn
         l = l-1
      end do

************************************************************************
*                               END                                    *
************************************************************************

      END SUBROUTINE interp_CMFGEN
