C ----------------------------------------------------------------------
      SUBROUTINE DIFFINITONDESEXCIT
C ----------------------------------------------------------------------
C Routine appelee par DIFFUSION
C Permet l'initialisation des quantites physiques necessaires 
C au calcul du transport par les ondes
C
C C.Charbonnel (17 mars 2003)
C C.Charbonnel (1 decembre 2006) : Modifications pour calcul de l'excitation
C                                  et des flux 
C
C S.Talon, 25 novembre 2004
C Attention: dans le code d'evolution, certaines variables sont definies
C            aux limites des couches et d'autres, au centre. 
C CENTRE:  p,ro,cs,t,cp,khi
C LIMITES: rray,grav,lum,m,hp,ht,rkonvL,rkonv

C CC A VERIFIER :
C Centre ou limites : xluml 
C
C exemple pour nmod=7
C haut: numerotation evolution (i)  bas: numerotation ondes (j)
C                               nmod
C       1   2   3   4   5   6   7
C       | 1 | 2 | 3 | 4 | 5 | 6 |
C       | 6 | 5 | 4 | 3 | 2 | 1 |
C       7   6   5   4   3   2   1
C ncouche
C aux points de grille on a: i=nmod-j+1
C            au centre on a: i=nmod-j
C
C ATTENTION: Dans le code de circulation meridienne la notation au centre 
C            des couches est decalee
C       | 6 | 5 | 4 | 3 | 2 | 1 |  --> ondes       (i)
C       | 2 | 3 | 4 | 5 | 6 | 7 |  --> circulation (k)
C aux points de grille on a: i=k
C            au centre on a: i=nmod-k+1
C ----------------------------------------------------------------------
      implicit none

      include 'evolpar.star'
      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.grad'
      include 'evolcom.teq'
      include 'evolcom.transp'
c      include 'evolcom.transpondesexcit'
      include 'evolcom.igw'

      include 'evolcom.therm'
      include 'evolcom.therm2'
      include 'evolcom.var'
C Modif CC ondes (19/10/07)
      include 'evolcom.eng'

      integer i, j, ncouche

      double precision xnumol,ddmax
      double precision phiKSm_o(nsh),deltaKSm_o(nsh)

C --------------------------------------------
C Lecture des donnees et inversion des couches 
C --------------------------------------------

      print *,"In DIFFINITONDESEXCIT"

      ncouche = nmod
      rtot = r(nmod)
      gs = gmr(nmod)

C Remise a zero de toutes les variables utiles au calcul des ondes
C ----------------------------------------------------------------
      do i = 1,nmod
         rray(i) = r(i)/rtot
C         vrray(i) = vr(i)/rtot
         grav(i) = gmr(i)/gmr(nmod)
      enddo

C      do i = nmod,1,-1
      do i = 1,nmod
         rray_o(i) = 0.d0
         grav_o(i) = 0.d0
         rmass_o(i) = 0.d0
         dlnmu_o(i) = 0.d0
         hp_o(i) = 0.d0
         xnuvv_o(i) = 0.d0
         vconv_o(i) = 0.d0
         xluml_o(i) = 0.d0
C         xnum_o(i) = 0.d0
         rkonv_o(i) = 0.d0
         rraym_o(i) = 0.d0
         rraym2_o(i) = 0.d0
         hnenn_o(i) = 0.d0
	 dm_o(i) = 0.d0
         xNt_o(i) = 0.d0
         xNmu_o(i) = 0.d0
         xN_o(i) = 0.d0
         brunt_o(i)= 0.d0
         rdens_o(i) = 0.d0
         rtemp_o(i) = 0.d0
         xpress_o(i) = 0.d0 
         cs_o(i) = 0.d0 
         khi_o(i) = 0.d0
         xnurad_o(i) = 0.d0
         xcp_o(i) = 0.d0
         xKm_o(i) = 0.d0
         fint(i) = 0.d0
      enddo

C Quantites definies aux couches dans le code d'evolution stellaire
C -----------------------------------------------------------------
      do i = 1,nmod
         j=nmod-i+1
         rray_o(j) = rray(i)
         grav_o(j) = grav(i)
         rmass_o(j) = m(i)
         dlnmu_o(j) = abmu(i)

         deltaKSm_o(j) = deltaKSm(i)
         phiKSm_o(j) = phiKSm(i)

         hp_o(j) = hp(i)
         xnuvv_o(j) = xnuvv(i)
         vconv_o(j) = sconv(i)
         xluml_o(j) = lum(i)
C         xnum_o(j) = xnum(i)
      enddo

      rkonv_o(1)=abla(nmod)-abad(nmod)
      do i = 2,nmod
         j=nmod-i+1
         rkonv_o(j) = abla(i)-dsqrt(abad(i)*abad(i-1))
C      print *,"rkonv_o",rkonv_o(j),abla(i),abad(i),abad(i-1)
      enddo
      
C Calcul de quantites physiques
C -----------------------------
      do i=1,ncouche-1
        rraym_o(i) = 0.5d0*(rray_o(i)+rray_o(i+1))
        rraym2_o(i) = rraym_o(i)*rraym_o(i)*rtot*rtot
        hnenn_o(i) = 1.d0/(rray_o(i+1)-rray_o(i))
	dm_o(i) = rmass_o(i)-rmass_o(i+1)
        xNt_o(i) = -grav_o(i)*gs/hp_o(i)*rkonv_o(i)*deltaKSm_o(i)
        xNmu_o(i) = grav_o(i)*gs/hp_o(i)*dlnmu_o(i)*phiKSm_o(i)
        xN_o(i) = xNt_o(i)+xNmu_o(i)
c        write (998,998) ncouche-i+1,xNt_o(i),grav_o(i),gs,hp_o(i),
c     &       rkonv_o(i),abla(ncouche-i+1),abrad(ncouche-i+1),
c     &       abm(ncouche-i+1)
c 998    format (i4,33(1x,1pe17.10))
        if (xN_o(i).gt.0.d0) then 
            brunt_o(i) = dsqrt(xN_o(i))
        else
            brunt_o(i) = 0.d0
        endif
      enddo 
      i=ncouche    
        xNt_o(i) = -grav_o(i)*gs/hp_o(i)*rkonv_o(i)*deltaKSm_o(i)
        xNmu_o(i) = grav_o(i)*gs/hp_o(i)*dlnmu_o(i)*phiKSm_o(i)
        xN_o(i) = xNt_o(i)+xNmu_o(i)
        if (xN_o(i).gt.0.d0) then
            brunt_o(i) = dsqrt(xN_o(i))
        else
            brunt_o(i) = 0.d0
        endif
        
c        do i = ncouche,1,-1
c           print *,i,ncouche-i,xNt_o(i),grav_o(i),gs,hp_o(i),rkonv_o(i)
c           print *,i,ncouche-i,xNmu_o(i),grav_o(i),gs/hp_o(i),dlnmu_o(i)
c        enddo
c        stop "init"

C Quantites definies aux inter-couches dans le code d'evolution stellaire
C Ces quantites doivent etre reconstruites aux points de grille
C -----------------------------------------------------------------------
      do i = 2,nmod-1
         dm_m1(i)=dm_o(i)/(dm_o(i-1)+dm_o(i))
         dm_p1(i)=dm_o(i-1)/(dm_o(i-1)+dm_o(i))
      enddo
      
      do i = 2,nmod-1
         j=nmod-i+1
         rdens_o(j) = ro(i-1)*dm_m1(j)+ ro(i)*dm_p1(j)
         rtemp_o(j) =  t(i-1)*dm_m1(j)+  t(i)*dm_p1(j)
         xpress_o(j) = p(i-1)*dm_m1(j)+  p(i)*dm_p1(j)
         cs_o(j) = cs(i-1)*dm_m1(j)+ cs(i)*dm_p1(j)
         khi_o(j)   = khi(i-1)*dm_m1(j)+khi(i)*dm_p1(j)
         xnurad_o(j) = khi_o(j)*rtemp_o(j)/rdens_o(j)/4.4937758d21
C  val num : 5*c**2
         ddmax = 1.5964182d-8*rtemp_o(j)**1.5d0/dsqrt(rdens_o(j))
C  val num : 1.5/e**3*(mp*k**3/pi)**0.5
         xnumol = 2.1688529d-15*rtemp_o(j)**2.5d0/rdens_o(j)/dlog(ddmax)
C         xnum_o(j) = xnumol+xnurad_o(j)
         xcp_o(j) = cp(i-1)*dm_m1(j)+ cp(i)*dm_p1(j)
      enddo
         
C -------------------------------
C Calcul de quantites physiques
C -------------------------------
      do i=1,ncouche-1
        j=nmod-i
        xKm_o(i)=khi(j)/cp(j)/ro(j)          ! Diff. thermique
        fint(i)= -(xKm_o(i)+xnuvv_o(i))*brunt_o(i)*xNt_o(i)/hnenn_o(i)/
     &          rraym_o(i)**3/rtot**2
        
       write (990,990)i,fint(i),xKm_o(i),xnuvv_o(i),brunt_o(i),xNt_o(i),
     &           hnenn_o(i),rraym_o(i),rtot
 990    format (i4,30(1x,1pe17.10))
      enddo
c      stop
CC sortie_ondes <--> 174 <--> sortie_ondes <--> $FILEXw
C      open (unit = 174,file = 'sortie_ondes',status = 'unknown')
C      write (174,1740)
C      write(174,*) r(nmod),gmr(nmod),nmod,omega_S
C      do i=1,ncouche
C         write(174,1704) rray_o(i),xpress_o(i),rdens_o(i),
C     &      grav_o(i),cs_o(i),vconv_o(i),rtemp_o(i),xcp_o(i),
C     &      khi_o(i),xluml_o(i),rmass_o(i),dlnmu_o(i),hp_o(i),
C     &      dm_o(i),xnum_o(i),rkonv_o(i)
CC     &      dm_o(i),xnum_o(i),rkonv_o(i),om_o(i)
C      enddo

C      close (174)
C 1740 format (//,5x,'VARIABLES NEEDED FOR CALCULATIONS OF THE
C     &  INTERNAL GRAVITY WAVES :',//)
C 1704 format (16(1x,1pe17.10))
CC 1704 format (17(1x,1pe17.10))

C      do i=1,ncouche
C        print *,"rray_o(i),grav_o(i)",rray_o(i),grav_o(i)
C     enddo
C
C      write (*,*) 'Dans diffinitondesexcit.f'
C      do i = 1,ncouche
C       write (*,*) 'omega_S,rraym2_o(i),om_o(i),dm_o(i)',
C     &   omega_S,rraym2_o(i),om_o(i),dm_o(i)
C      enddo


      print *,"Sortie de DIFFINITONDESEXCIT"
      return
      end
