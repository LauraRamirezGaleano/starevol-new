

************************************************************************

      SUBROUTINE diffturbul (lover,rho_bcz,DHE,Dturbul)

************************************************************************
*     Compute the diffusion coefficient associated with turbulence     *
*     (Richer et al. 2000; Richard et al.2005)                         *
*                                                                      *
* $LastChangedDate:: 2019-09-11    $ *
* $Author:: TD                                                 $ *
* $Rev:: 0                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.cons'
      include 'evolcom.conv'
      include 'evolcom.mass'
      include 'evolcom.nuc'
      include 'evolcom.spec'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.var'

      double precision DHE,Dturbul1,Dturbul2,DHE1,Dturbul3,Dturbul

c      integer istart,iend,imin,imax,idest,icz
       integer i,j,jshell
       integer om_turbul,n
       integer fix_layer
       integer lover
       parameter (om_turbul = 400, n = 3)

       double precision rho0, DHE0, Tfix,rho_bcz,DHE10,T0
       
       dimension Dturbul(nsh)
       dimension DHE(nsh),Dturbul1(nsh),Dturbul2(nsh)
     $      ,DHE1(nsh),Dturbul3(nsh)

      fix_layer = 0
*     Diffusion coefficient of turbulente put to zero
      do j = 1,nmod
         Dturbul(j) = 0.d0
* Analytical approximation (Richer et al. 2000) footnote 2
c         DHE1(j) = 3.3e-15*t(j)**(2.5)/(4*ro(j)*log(1
c     $        +1.125e-16*t(j)**3/ro(j)))
      enddo
      
      Tfix = 10**(6.42)         ! test 30/01/2020 - > solar twins
c      Tfix = 10**(6.50)
c       Tfix = 10**(6.60)
c     Tfix = 10**(6.425)                  ! fixation point of the turbulence function of temperature -> Soleil
      do jshell = 1,nmod
         if (t(jshell).le.Tfix) then
            if (fix_layer.ne.0) exit
            fix_layer = jshell
c            print *, 'fixation layer',fix_layer
            rho0 = ro(jshell)   ! density corresponding to the fixation point
            T0 = t(jshell)
            DHE0 = DHE(jshell)  ! Coefficient of diffusion of helium at the point of fixation
c            DHE10 = DHE1(jshell) ! Analytical Coefficient of diffusion of helium at the point of fixation
         endif
* Analytical approximation (Richer et al. 2000) footnote 2         
         DHE10 = 3.3e-15*T0**(2.5)/(4*rho0*log(1
     $        +1.125e-16*T0**3/rho0))
c      print *, 'DHE01',DHE0   
      enddo

      print *, 'rho0',rho0,'T0',T0, fix_layer
      print *, 'tsh',t(fix_layer),'rhosh',ro(fix_layer),'rhobcz',
     $     rho_bcz,ro(341)
c      stop
*     Calcul of the turbulence coefficient (eq. 1 and 2 Richer et al.2000)
      do jshell = 1,nmod 
         Dturbul(jshell) = om_turbul*DHE0*(rho0/ro(jshell))**n    ! eq. 2 Richard et al. 2005
c     Dturbul1(jshell) = om_turbul*DHE1(jshell)*(rho0/ro(jshell))**n
         Dturbul1(jshell) = om_turbul*DHE10*(rho0/ro(jshell))**n ! formule à partir DHE analytique
c         Dturbul2(jshell) = 12500*(rho_bcz/ro(jshell))**n ! eq.3 Richard et al. 2005
         Dturbul2(jshell) = 5000*(rho_bcz/ro(jshell))**n ! test 04/02/2020
         Dturbul3(jshell) = om_turbul*DHE(jshell)*(rho0/ro(jshell))**n ! formule à partir profil DHE Thoul
c     print *, 'jshell',jshell,'Dturbul',log(abs(Dturbul(jshell)))
c$$$         write(550,'(1x,i4,11(1x,1pe11.4))'),jshell,r(jshell),t(jshell)
c$$$     $        ,log10(abs(Dturbul(jshell))),log10(abs(Dturbul1(jshell)))
c$$$     $        ,log10(abs(Dturbul2(jshell))) ,DHE(jshell),DHE0
c$$$     $        ,DHE1(jshell) ,log10(abs(Dturbul3(jshell))),DHE10
c$$$  $        ,ro(jshell)
         if (lover.eq.60.or.lover.eq.70.or.lover.eq.72) then
            Dturbul(jshell) = Dturbul1(jshell)
         else if (lover.eq.61.or.lover.eq.71.or.lover.eq.73) then
            Dturbul(jshell) = Dturbul2(jshell)
         endif
c$$$         write(550,'(1x,i4,11(1x,1pe11.4))'),jshell,r(jshell),t(jshell)
c$$$     $        ,ro(jshell),log10(abs(Dturbul(jshell)))
      enddo
c      stop
      return
      end
