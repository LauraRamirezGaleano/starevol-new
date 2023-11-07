
      SUBROUTINE load_CMFGEN(folder_path)
      
************************************************************************
*                                                                      *
* Load the CMFGEN atmosphere tables and store them in the COMMON       *
* /CMFGEN_mass_atm/ arrays.                                            *
*                                                                      *
************************************************************************
*                                                                      *
* Inputs/outputs                                                       *
* --------------                                                       *
*                                                                      *
* folder_path (CHARACTER)                                              *
*    path of the folder where CMFGEN files are stored                  *
*                                                                      *
************************************************************************
*                                                                      *
* $LastChangedDate:: 2023-06-06 16:39:00 +0000 (Tue, 06 Jun 2023)    $ *
* $Author:: VOJE                                                     $ *
* $Rev:: 110                                                         $ *
*                                                                      *
************************************************************************

      IMPLICIT NONE

************************************************************************
*                              INCLUDES                                *
************************************************************************
      
      INCLUDE 'evolpar.star'
      INCLUDE 'evolcom.atm'

************************************************************************
*                           DEFINITIONS                                *
************************************************************************

***   Inputs/outputs.
************************************************************************
      
      CHARACTER(LEN=58), INTENT(IN) :: folder_path

***   Temporary variables.
************************************************************************

      CHARACTER(LEN=21) :: atm_file_name
      CHARACTER(LEN=5) :: char_Teff
      CHARACTER(LEN=3) :: char_Mdot, char_gs
      INTEGER :: i, j, k, l, n, nlines, id
      INTEGER :: stat
      LOGICAL :: flag

***   Teff, Mdot and gs lists corresponding to the CMFGEN atmosphere
***   tables.
************************************************************************
      
      INTEGER, ALLOCATABLE :: Teff_list(:)
      DOUBLE PRECISION, ALLOCATABLE :: Mdot_list(:), gs_list(:)

************************************************************************
*                   LOAD CMFGEN TABLES FUNCTION                        *
************************************************************************

***   1. Create the Teff, Mdot and gs lists.
************************************************************************
      
      nb_Teff_CMFGEN = 30
      nb_Mdot_CMFGEN = 24
      nb_gs_CMFGEN = 38
      ALLOCATE(Teff_list(nb_Teff_CMFGEN))
      ALLOCATE(Mdot_list(nb_Mdot_CMFGEN))
      ALLOCATE(gs_list(nb_gs_CMFGEN))
      DO i=1,nb_Teff_CMFGEN
         Teff_list(i) = 20000 + i*1000
      END DO
      DO i=1,nb_Mdot_CMFGEN
         Mdot_list(i) = 4.0 + 0.1*i
      END DO
      DO i=1,nb_gs_CMFGEN
         gs_list(i) = 2.4 + 0.05*i
      END DO

***   2. Display a console message.
************************************************************************
      
      WRITE(*,*) "Reading CMFGEN atmosphere files"
      WRITE(*,*) ""

***   3. Fill the CMFGEN tables.
************************************************************************
      
      id = 1
      DO i=1,nb_Teff_CMFGEN
         DO j=1,nb_Mdot_CMFGEN
            DO k=1,nb_gs_CMFGEN
               flag = .TRUE.
               WRITE(char_Teff,FMT="(I5)") Teff_list(i)
               WRITE(char_Mdot,FMT="(F3.1)") Mdot_list(j)
               WRITE(char_gs,FMT="(F3.1)") gs_list(k)
               atm_file_name = "CMFGEN_"//char_Teff//"_-"//char_Mdot//
     &              "_"//char_gs
               DO l=1,50
                  IF (k.EQ.2*l+1) THEN
                     flag = .FALSE.
                  END IF
               END DO
               stat = 1
               OPEN (UNIT=121,FILE=folder_path//atm_file_name,
     &              IOSTAT=stat,STATUS="OLD")
               IF (stat.EQ.0.AND.flag) THEN
                  CALL system("wc -l "//folder_path//
     &                 atm_file_name//" | awk '{print $1}' > nlines")
                  OPEN(UNIT=131,FILE="nlines",STATUS="OLD")
                  READ(131,*) nlines
                  nlines = nlines-2
                  CLOSE(131)
                  avail_CMFGEN(id,:) = (/i,j,k,nlines/)
                  id = id+1
                  DO l=1,2
                     READ(121,*)
                  END DO
                  DO l=1,nlines
                     READ(121,FMT="(E11.1, 5(1X,E11.1))") 
     &                    tau_CMFGEN(i,j,k,l), r_CMFGEN(i,j,k,l), 
     &                    T_CMFGEN(i,j,k,l), rho_CMFGEN(i,j,k,l),
     &                    ne_CMFGEN(i,j,k,l), v_CMFGEN(i,j,k,l)
                     r_CMFGEN(i,j,k,l) = r_CMFGEN(i,j,k,l)*6.957d10
                  END DO
               END IF
               CLOSE(121)
            END DO
         END DO
      END DO
      max_id = id-1
      CALL system("rm nlines")
      
************************************************************************
*                               END                                    *
************************************************************************
      
      END SUBROUTINE load_CMFGEN
