


************************************************************************

      SUBROUTINE denucl (tk,rok,mueinvk,x,ksh,eng,engdt,engdro)

************************************************************************
* Compute the nuclear energy production rates and its derivatives      *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      integer nflu,ksh

      double precision tk,rok,mueinvk,x
      double precision tkp,rokp,eng,engro,engdt,engdro,engrocap
      double precision v
      double precision dtk,drok,engro1,eng1,engrocap1,engro2,eng2,
     &     engrocap2

      dimension x(nsp),v(nreac)

      if (tk.ge.1.d5) then
***   rho,T+dT
         dtk = tk*1.d-4
         tkp = tk+dtk
         nflu = 1
         call vit (tkp,rok,mueinvk,ksh,v,0)
         call nuceng (nflu,x,v,ksh,eng1,engro1,engrocap1)
***   rho+drho,T
         nflu = 2
         drok = rok*1.d-4
         rokp = rok+drok
         call vit (tk,rokp,mueinvk,ksh,v,2)
         call nuceng (nflu,x,v,ksh,eng2,engro2,engrocap2)
***   rho,T
         nflu = 0
         call vit (tk,rok,mueinvk,ksh,v,0)
         call nuceng (nflu,x,v,ksh,eng,engro,engrocap)
         engdt = (eng1-eng)/dtk
         engdro = engro/rok+(engrocap2-engrocap)/drok
      else
***   rho,T
         nflu = 0
         call vit (tk,rok,mueinvk,ksh,v,0)
         call nuceng (nflu,x,v,ksh,eng,engro,engrocap)
         engdt = 0.d0
         engdro = 0.d0
      endif

      return
      end
