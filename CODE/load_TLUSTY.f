      SUBROUTINE load_TLUSTY(folder_path)
      
************************************************************************
*                                                                      *
* Load the TLUSTY atmosphere tables and store them in the COMMON       *
* /TLUSTY_mass_atm/ arrays.                                            *
*                                                                      *
************************************************************************
*                                                                      *
* Inputs/outputs                                                       *
* --------------                                                       *
*                                                                      *
* folder_path (CHARACTER)                                              *
*    path of the folder where TLUSTY files are stored                  *
*                                                                      *
************************************************************************
*                                                                      *
* $LastChangedDate:: 2023-04-04 16:59:00 +0000 (Tue, 04 Apr 2023)    $ *
* $Author:: VOJE                                                     $ *
* $Rev:: 110                                                         $ *
*                                                                      *
************************************************************************

      IMPLICIT NONE
      
************************************************************************
*                           Definitions                                *
************************************************************************

!     Inputs/outputs
      CHARACTER(LEN=51), INTENT(IN) :: folder_path

!     Temporary variables
      CHARACTER(LEN=16) :: atm_file_name
      CHARACTER(LEN=5) :: char_Teff
      CHARACTER(LEN=3) :: char_g
      INTEGER :: i, j, k
      INTEGER :: stat

!     Teff and Mdot lists corresponding to the TLUSTY atmosphere tables
      INTEGER, ALLOCATABLE :: Teff_list(:), g_list(:)

!     COMMONS : TLUSTY atmosphere tables and associated numbers
      INTEGER :: nb_Teff_TLUSTY, nb_g_TLUSTY
      DOUBLE PRECISION :: tau_TLUSTY,
     &     T_TLUSTY, rho_TLUSTY
      COMMON /TLUSTY_mass_atm/ tau_TLUSTY(12,9,50),
     &     T_TLUSTY(12,9,50), rho_TLUSTY(12,9,50)
      COMMON /TLUSTY_mass_atm_nb/ nb_Teff_TLUSTY, nb_g_TLUSTY

************************************************************************
*                   Load TLUSTY tables function                        *
************************************************************************

!     1. Create the Teff and Mdot lists
      nb_Teff_TLUSTY = 12
      nb_g_TLUSTY = 9
      ALLOCATE(Teff_list(nb_Teff_TLUSTY))
      ALLOCATE(g_list(nb_g_TLUSTY))
      DO i=1,nb_Teff_TLUSTY
         Teff_list(i) = 25000 + i*2500
      END DO
      DO i=1,nb_g_TLUSTY
         g_list(i) = 275 + 25*i
      END DO

!     2. Display a console message
      WRITE(*,*) "Reading TLUSTY atmosphere files"
      WRITE(*,*) ""

!     3. Fill the TLUSTY tables
      DO i=1,nb_Teff_TLUSTY
         DO j=1,nb_g_TLUSTY
            WRITE(char_Teff,FMT="(I5)") Teff_list(i)
            WRITE(char_g,FMT="(I3)") g_list(j)
            atm_file_name = "G"//char_Teff//"g"//char_g//
     &           "v10.11"
            stat = 1
            OPEN (UNIT=12,FILE=folder_path//atm_file_name,
     &           IOSTAT=stat,STATUS="OLD")

!     .. 3.1. If the file doesn't exist, fill the table with incorrect
!     ....... profiles
            IF (.NOT.stat.EQ.0) THEN
               CLOSE(12)
               DO k=1,50
                  tau_TLUSTY(i,j,k) = 8.498d-6*exp(k*1.628d-1)
                  T_TLUSTY(i,j,k) = 0.d0
                  rho_TLUSTY(i,j,k) = 0.d0
               END DO
               
!     .. 3.2. Else, read the file and fill the tables with the
!     ....... corresponding TLUSTY profiles 
            ELSE
               DO k=1,50
                  READ(12,FMT="(20X,E10.1,14X,E10.1,14X,E10.1)") 
     &                 tau_TLUSTY(i,j,k), T_TLUSTY(i,j,k),
     &                 rho_TLUSTY(i,j,k)
               END DO
               CLOSE(12)
               
            END IF
         END DO
      END DO
      
************************************************************************
*                               End                                    *
************************************************************************
      
      END SUBROUTINE load_TLUSTY
