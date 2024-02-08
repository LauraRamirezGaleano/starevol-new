      Subroutine interpZ
*******************************************************************************************

      implicit none

      include 'evolpar.star'
      include 'evolcom.atm'

      integer i,j,k
      integer filesize
      double precision rtau,rT,rhotab,Fconvtab
      double precision T_Zeff(6),rho_Zeff(6),Fconv_Zeff(6)
      double precision Ttabeff,gtabeff,Ztabeff
      double precision FeH,DT,DDT,sTZ,srhotabZ,sFconvtabZ

      double precision tableTZ,tablerhoZ,tableFcZ

      dimension DT(512),DDT(512)

      common /atmospheres/ rtau(127,12,59,6),rT(127,12,59,6),
     &     rhotab(127,12,59,6),Fconvtab(127,12,59,6),Ttabeff(59),
     &     gtabeff(12),Ztabeff(6),filesize
      common /atmospheres2/ tableTZ(127,59,10),tablerhoZ(127,59,10),
     &     tableFcZ(127,59,10)
      common /metal/ FeH

      print *,FeH

      sTZ = 0.d0
      srhotabZ = 0.d0
      sFconvtabZ = 0.d0


      do i = 1,nb_Teff
         do j = 1,nb_geff
            do k = 1,filesize

               T_Zeff(:) = rT(k,j,i,:)
               call splineatm (Ztabeff,T_Zeff,nb_Zeff,1.d50,1.d50,
     &              DDT,DT)
               call splintatm (Ztabeff,T_Zeff,DDT,nb_Zeff,FeH,sTZ)
               tableTZ(k,i,j) = sTZ

               rho_Zeff(:) = rhotab(k,j,i,:)
               call splineatm (Ztabeff,rho_Zeff,nb_Zeff,1.d50,1.d50,
     &              DDT,DT)
               call splintatm (Ztabeff,rho_Zeff,DDT,nb_Zeff,FeH,
     &              srhotabZ)
               tablerhoZ(k,i,j) = srhotabZ

               Fconv_Zeff(:) = Fconvtab(k,j,i,:)
               call splineatm (Ztabeff,Fconv_Zeff,nb_Zeff,1.d50,1.d50,
     &              DDT,DT)
               call splintatm (Ztabeff,Fconv_Zeff,DDT,nb_Zeff,FeH,
     &              sFconvtabZ)
               tableFcZ(k,i,j) = sFconvtabZ

            enddo
         enddo
      enddo

      return
      end
