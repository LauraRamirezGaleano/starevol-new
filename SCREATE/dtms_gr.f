
      PROGRAM DTMS

************************************************************************
* This code is computing the typical maximum time-step suited to       *
* compute the main sequence phase of a star of initial mass M and      *
* of initial metallicity Z in roughly 200 evolution models;            *
* this corresponds to a maximum time-step of 4.61d7 for a M = 1.00 and *
* Z = 0.020 star (dtref)                                               *
*                                                                      *
* The assumed formula to perform this computation is                   *
* dt(M,Z) = dt(1,0.02) * (M/1)**-2.9 * (Z/0.02)**0.1 , where           *
*    M is given in solar units                                         *
*                                                                      *
* This program needs, as input: M and Z                                *
* The only output is: dt(M,Z) to put as dtmax, in starevol.par         *
************************************************************************

      implicit none

      double precision m,z
      double precision alpha,beta,dtref
      double precision dtmax

*** Assumed values:

      alpha = -2.90d0
      beta = 0.10d0
      dtref = 1.00d7

*** Input values:

      read (*,*) m,z
c      write(*,*) 'm,z',m,z

*** Computations:

      dtmax = dtref*m**alpha*(z/0.020d0)**beta
      
*** Writing output value:

      write (*,2100) dtmax

 2100 format (1pe9.3)

      end
