* -*-fortran-*-
*                                                                      *
* $LastChangedDate:: 2014-02-04 19:38:17 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 18                                                          $ *
*                                                                      *
***   chemicals, nuclear
      integer nbz,nis,nsp,nre,nscr,nreac,nprint,nprintc
      integer nbeta,nbetaec,nvit,nlibeb
***   electron captures
      integer nte,nroe,nrce,nte1,nroe1,nrce1
c..   idneut: neutron capture reactions, nneut : reactions involving n
      integer idneut,nneut

c..   index nbz corresponds to the charge of heavy	
      parameter (nbz = 18,nis = 54,nsp = nis+1)
c      parameter (nbz = 21,nis = 54,nsp = nis+1) ! modif TD Fev.2019
      parameter (nbeta = 16,nbetaec = 6,idneut = 50,nneut = 62)
      parameter (nlibeb = 3)
      parameter (nrce1 = 0,nte1 = 10,nroe1 = 10) 
      parameter (nrce = 12,nte = 11,nroe = 10) 
      parameter (nscr = 2*nbz+1,nvit = 90,nre = 167,nreac=nre+nrce)
      parameter (nprint = 44,nprintc = nprint)

***   metallicity
      double precision zkint
      common /metallicity/ zkint

***   convection
      integer nmaxconv,nmaxenu,ntypeconv
      parameter (nmaxconv = 500,ntypeconv = 12,nmaxenu = 25)

***   opacity tables
c      integer mlth,nlth,itlth
      integer nth,nrh,itabh
c..   molecules included in opacity calculation
      integer nmole
      integer neqchimmax,natomsmax
      parameter (neqchimmax = 20, natomsmax = 5)

c..   atomes included in ionization calculation for opacity 
      integer nioniz

c      parameter (mlth = 19,nlth = 85,itlth = 120)
c      parameter (mlth = 19,nlth = 85,itlth = 155)
C Warning: dimension of opacity tables change with solar comp.
c#ifdef GRID 
c      parameter (mlth = 19,nlth = 63,itlth = 104)
c#else
c      parameter (mlth = 19,nlth = 85,itlth = 120)
c#endif
      parameter (nrh = 17,nth = 11,itabh= 77)
      parameter (nmole = 7, nioniz = 24)

***   shells
c Note: nsh is also defined in evolcom.igw. Update both files if this value changes.
      integer nsh,neq,neq1,intmax

      parameter (nsh = 4500, neq = 5, neq1 = 3*neq+1)
c      parameter (nsh = 8000, neq = 5, neq1 = 3*neq+1)
      parameter (intmax = 101)

***   input logical unit
      integer nout
      common /logicunit/ nout

 
