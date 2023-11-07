      SUBROUTINE interp_TLUSTY(gs,interp_T,interp_rho,t_eff)
      
************************************************************************
*                                                                      *
* Interpolate T and rho profiles with the TLUSTY atmosphere tables to  *
* produce the atmosphere profiles used as boundary conditions in the   *
* Newton-Raphson computation.                                          *
*                                                                      *
************************************************************************
*                                                                      *
* Inputs/outputs                                                       *
* --------------                                                       *
*                                                                      *
* gs (DOUBLE PRECISION)                                                *
*     value of gs (gravity field at the surface) at the current        *
*     time-step                                                        *
*                                                                      *   
* interp_T (DOUBLE PRECISION, (nsh))                                   *
*    interpolated T profile                                            *
*                                                                      *
* interp_rho (DOUBLE PRECISION, (nsh))                                 *
*    interpolated rho profile                                          *
*                                                                      *
* t_eff (DOUBLE PRECISION)                                             *
*    value of <Teff>, mean effective temperature for tau inside the    *
*    interval [30:100]                                                 *
*                                                                      *
************************************************************************
*                                                                      *
* $LastChangedDate:: 2023-04-05 10:41:00 +0000 (Wen, 05 Apr 2023)    $ *
* $Author:: VOJE                                                     $ *
* $Rev:: 110                                                         $ *
*                                                                      *
************************************************************************

      IMPLICIT NONE

************************************************************************
*                              Includes                                *
************************************************************************

      INCLUDE 'evolpar.star'
      INCLUDE 'evolcom.atm'
      INCLUDE 'evolcom.mod'
      INCLUDE 'evolcom.surf'
      INCLUDE 'evolcom.teq'
      INCLUDE 'evolcom.var'

************************************************************************
*                             Definitions                              *
************************************************************************

!     Inputs/outputs
      DOUBLE PRECISION, INTENT(IN) :: gs
      DOUBLE PRECISION, INTENT(OUT) :: t_eff
      DOUBLE PRECISION, INTENT(OUT) :: interp_T(nsh), interp_rho(nsh)

!     Temporary variables
      INTEGER :: i, j, k, l, m, stat, nb_elem, nb_g_tmp
      DOUBLE PRECISION :: Teffinterp
      DOUBLE PRECISION :: tmp_tau_C(50), tmp_T_C(50), tmp_rho_C(50)
      DOUBLE PRECISION :: DDT(50), DT(50), DDrho(50), Drho(50)
      DOUBLE PRECISION, ALLOCATABLE :: DDT2(:), DT2(:), DDrho2(:),
     &     Drho2(:)
      DOUBLE PRECISION, ALLOCATABLE :: DDT3(:), DT3(:),
     &     DDrho3(:), Drho3(:)
      DOUBLE PRECISION, ALLOCATABLE :: new_T(:,:), tmp_T_2D(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: new_tmp_T(:), tmp_T(:)
      DOUBLE PRECISION, ALLOCATABLE :: new_rho(:,:), tmp_rho_2D(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: new_tmp_rho(:), tmp_rho(:)
      DOUBLE PRECISION, ALLOCATABLE :: g_list_2(:), tmp_T_2(:),
     &     tmp_rho_2(:)

!     Teff and g lists corresponding to the TLUSTY atmosphere tables
      DOUBLE PRECISION, ALLOCATABLE :: Teff_list(:), g_list(:)

!     COMMONS : TLUSTY atmosphere tables and associated numbers
      INTEGER :: nb_Teff_TLUSTY, nb_g_TLUSTY
      DOUBLE PRECISION :: tau_TLUSTY,
     &     T_TLUSTY, rho_TLUSTY
      COMMON /TLUSTY_mass_atm/ tau_TLUSTY(12,9,50),
     &     T_TLUSTY(12,9,50), rho_TLUSTY(12,9,50)
      COMMON /TLUSTY_mass_atm_nb/ nb_Teff_TLUSTY, nb_g_TLUSTY

************************************************************************
*                   TLUSTY interpolation function                      *
************************************************************************
      
!     1. Create the Teff and g lists
      nb_Teff_TLUSTY = 12
      nb_g_TLUSTY = 9
      ALLOCATE(Teff_list(nb_Teff_TLUSTY))
      ALLOCATE(g_list(nb_g_TLUSTY))
      DO i=1,nb_Teff_TLUSTY
         Teff_list(i) = 25000.d0 + i*2500.d0
      END DO
      DO i=1,nb_g_TLUSTY
         g_list(i) = 2.75d0 + 0.25d0*i
      END DO
      
!     2. Allocate arrays
      ALLOCATE(tmp_T_2D(nb_g_TLUSTY,nsh))
      ALLOCATE(tmp_T(nb_g_TLUSTY))
      ALLOCATE(new_T(nb_Teff_TLUSTY,nsh))
      ALLOCATE(new_tmp_T(nb_Teff_TLUSTY))
      ALLOCATE(DT2(nb_g_TLUSTY))
      ALLOCATE(DDT2(nb_g_TLUSTY))
      ALLOCATE(DT3(nb_Teff_TLUSTY))
      ALLOCATE(DDT3(nb_Teff_TLUSTY))
      ALLOCATE(tmp_rho_2D(nb_g_TLUSTY,nsh))
      ALLOCATE(tmp_rho(nb_g_TLUSTY))
      ALLOCATE(new_rho(nb_Teff_TLUSTY,nsh))
      ALLOCATE(new_tmp_rho(nb_Teff_TLUSTY))
      ALLOCATE(Drho2(nb_g_TLUSTY))
      ALLOCATE(DDrho2(nb_g_TLUSTY))
      ALLOCATE(Drho3(nb_Teff_TLUSTY))
      ALLOCATE(DDrho3(nb_Teff_TLUSTY))

!     3. Put everything on the same grid and interpolate in g      
!     .. 3.1. Put everything on the same grid
      DO i=1,nb_Teff_TLUSTY
         DO j=1,nb_g_TLUSTY
            tmp_tau_C(:) = tau_TLUSTY(i,j,:)
            tmp_T_C(:) = T_TLUSTY(i,j,:)
            tmp_rho_C(:) = rho_TLUSTY(i,j,:)
            CALL splineatm (tmp_tau_C,tmp_T_C,50,1.d50,
     &           1.d50,DDT,DT)
            CALL splineatm (tmp_tau_C,tmp_rho_C,50,1.d50,
     &           1.d50,DDrho,Drho)
            k = nmod
            DO WHILE(tau(k).LT.100)
               m = nmod-k+1
               CALL splintatm (tmp_tau_C,
     &              tmp_T_C,DDT,50,
     &              tau(k),tmp_T_2D(j,m))
               CALL splintatm (tmp_tau_C,
     &              tmp_rho_C,DDrho,50,
     &              tau(k),tmp_rho_2D(j,m))
               k = k-1
            END DO
         END DO
         
!     .. 3.2. Interpolate in g
         k = nmod
         new_T(i,:) = 0.d0
         new_rho(i,:) = 0.d0
         DO WHILE(tau(k).LT.100)
            m = nmod-k+1
            tmp_T(:) = tmp_T_2D(:,m)
            tmp_rho(:) = tmp_rho_2D(:,m)
            nb_g_tmp = 0
            DO j=1,nb_g_TLUSTY
               IF (tmp_T(j).GT.1d-20) THEN
                  nb_g_tmp = nb_g_tmp+1
               END IF
            END DO
            ALLOCATE(g_list_2(nb_g_tmp))
            ALLOCATE(tmp_T_2(nb_g_tmp))
            ALLOCATE(tmp_rho_2(nb_g_tmp))
            l = 1
            DO j=1,nb_g_TLUSTY
               IF (tmp_T(j).GT.1d-20) THEN
                  g_list_2(l) = g_list(j)
                  tmp_T_2(l) = tmp_T(j)
                  tmp_rho_2(l) = tmp_rho(j)
                  l = l+1
               END IF
            END DO
            DEALLOCATE(DT2)
            DEALLOCATE(DDT2)
            DEALLOCATE(Drho2)
            DEALLOCATE(DDrho2)
            ALLOCATE(DT2(nb_g_tmp))
            ALLOCATE(DDT2(nb_g_tmp))
            ALLOCATE(Drho2(nb_g_tmp))
            ALLOCATE(DDrho2(nb_g_tmp))
            CALL splineatm (g_list_2,tmp_T_2,nb_g_tmp,1.d50,
     &           1.d50,DDT2,DT2)
            CALL splintatm (g_list_2,tmp_T_2,DDT2,nb_g_tmp,
     &           gs,new_T(i,m))
            CALL splineatm (g_list_2,tmp_rho_2,nb_g_tmp,1.d50,
     &           1.d50,DDrho2,Drho2)
            CALL splintatm (g_list_2,tmp_rho_2,DDrho2,nb_g_tmp,
     &           gs,new_rho(i,m))
            DEALLOCATE(g_list_2)
            DEALLOCATE(tmp_T_2)
            DEALLOCATE(tmp_rho_2)
            k = k-1
         END DO
      END DO

!     4. Compute <Teff>, mean Teff for tau in [30:100]
      t_eff = 0.d0
      nb_elem = 0
      k = nmod
      DO WHILE (tau(k).LT.1)
         k = k-1
      END DO
      DO WHILE (tau(k).LT.100)
         m = nmod-k+1
         new_tmp_T(:) = new_T(:,m)
         CALL splineatm (new_tmp_T,Teff_list,nb_Teff_TLUSTY,1.d50,
     &        1.d50,DDT3,DT3)
         CALL splintatm (new_tmp_T,Teff_list,DDT3,nb_Teff_TLUSTY,
     &        T(k),Teffinterp)
         t_eff = t_eff + Teffinterp
         nb_elem = nb_elem + 1
         k = k-1
      END DO
      t_eff = t_eff/nb_elem
      IF (t_eff/t_eff.NE.1.d0) THEN
         WRITE(*,*) ""
         WRITE(*,*) "ERROR in interp_TLUSTY : "
     &        //"Cannot interpolate in Teff = NaN."
         WRITE(*,*) "nb_elem = ", nb_elem
         t_eff = 1d6
      END IF
      WRITE(*,*) ""
      WRITE(*,FMT="(1X,A,F8.2,A)") "Interpolating in TLUSTY tables at"//
     &     " Teff = ", t_eff, " K"
      WRITE(*,*) ""

!     5. Compute atmospheric T(tau) by interpolating in <Teff>
      k = nmod
      DO WHILE(tau(k).LT.100)
         m = nmod-k+1
         new_tmp_T(:) = new_T(:,m)
         new_tmp_rho(:) = new_rho(:,m)
         CALL splineatm (Teff_list,new_tmp_T,nb_Teff_TLUSTY,1.d50,
     &        1.d50,DDT3,DT3)
         CALL splintatm (Teff_list,new_tmp_T,DDT3,nb_Teff_TLUSTY,
     &        t_eff,interp_T(k))
         CALL splineatm (Teff_list,new_tmp_rho,nb_Teff_TLUSTY,1.d50,
     &        1.d50,DDrho3,Drho3)
         CALL splintatm (Teff_list,new_tmp_rho,DDrho3,nb_Teff_TLUSTY,
     &        t_eff,interp_rho(k))
         k = k-1
      END DO
      
************************************************************************
*                               End                                    *
************************************************************************

      END SUBROUTINE interp_TLUSTY
