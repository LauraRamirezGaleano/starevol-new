
      SUBROUTINE load_CMFGEN(folder_path)
      
************************************************************************
*                                                                      *
* Load the CMFGEN atmosphere tables and store them in the COMMON       *
* /CMFGEN_mass_atm/ arrays.                                            *
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

************************************************************************
*                           DEFINITIONS                                *
************************************************************************

***   Inputs/outputs.
************************************************************************
      
***   folder_path : path of the folder where CMFGEN files are stored
      character(len=58), intent(in) :: folder_path

***   Temporary variables.
************************************************************************

      character(len=21) :: atm_file_name
      character(len=5) :: char_Teff
      character(len=3) :: char_Mdot, char_gs
      integer :: i, j, k, l, n, nlines, id
      integer :: stat
      logical :: flag

***   Teff, Mdot and gs lists corresponding to the CMFGEN atmosphere
***   tables.
************************************************************************
      
      integer, allocatable :: Teff_list(:)
      double precision, allocatable :: Mdot_list(:), gs_list(:)

************************************************************************
*                   LOAD CMFGEN TABLES FUNCTION                        *
************************************************************************

***   1. Create the Teff, Mdot and gs lists.
************************************************************************
      
      nb_Teff_CMFGEN = 30
      nb_Mdot_CMFGEN = 24
      nb_gs_CMFGEN = 38
      allocate(Teff_list(nb_Teff_CMFGEN))
      allocate(Mdot_list(nb_Mdot_CMFGEN))
      allocate(gs_list(nb_gs_CMFGEN))
      do i=1,nb_Teff_CMFGEN
         Teff_list(i) = 20000 + i*1000
      end do
      do i=1,nb_Mdot_CMFGEN
         Mdot_list(i) = 4.0 + 0.1*i
      end do
      do i=1,nb_gs_CMFGEN
         gs_list(i) = 2.4 + 0.05*i
      end do

***   2. Display a console message.
************************************************************************
      
      write(*,*) "Reading CMFGEN atmosphere files"
      write(*,*) ""

***   3. Fill the CMFGEN tables.
************************************************************************
      
      id = 1
      do i=1,nb_Teff_CMFGEN
         do j=1,nb_Mdot_CMFGEN
            do k=1,nb_gs_CMFGEN
               flag = .true.
               write(char_Teff,FMT="(I5)") Teff_list(i)
               write(char_Mdot,FMT="(F3.1)") Mdot_list(j)
               write(char_gs,FMT="(F3.1)") gs_list(k)
               atm_file_name = "CMFGEN_"//char_Teff//"_-"//char_Mdot//
     &              "_"//char_gs
               do l=1,50
                  if (k.eq.2*l+1) then
                     flag = .false.
                  end if
               end do
               stat = 1
               OPEN (UNIT=121,FILE=folder_path//atm_file_name,
     &              IOSTAT=stat,STATUS="OLD")
               if (stat.eq.0.and.flag) then
                  call system("wc -l "//folder_path//
     &                 atm_file_name//" | awk '{print $1}' > nlines")
                  open(unit=131,file="nlines",status="OLD")
                  read(131,*) nlines
                  nlines = nlines-2
                  close(131)
                  avail_CMFGEN(id,:) = (/i,j,k,nlines/)
                  id = id+1
                  do l=1,2
                     read(121,*)
                  end do
                  do l=1,nlines
                     read(121,FMT="(E11.1, 5(1X,E11.1))") 
     &                    tau_CMFGEN(i,j,k,l), r_CMFGEN(i,j,k,l), 
     &                    T_CMFGEN(i,j,k,l), rho_CMFGEN(i,j,k,l),
     &                    ne_CMFGEN(i,j,k,l), v_CMFGEN(i,j,k,l)
                     r_CMFGEN(i,j,k,l) = r_CMFGEN(i,j,k,l)*6.957d10
                  end do
               end if
               close(121)
            end do
         end do
      end do
      max_id = id-1
      call system("rm nlines")
      
************************************************************************
*                               END                                    *
************************************************************************
      
      END SUBROUTINE load_CMFGEN
