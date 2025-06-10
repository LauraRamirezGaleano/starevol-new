      program c12c12
      implicit double precision (a-h,o-z)
      parameter (ngrid=60)
      dimension t8(ngrid)
      dimension rn(ngrid),rp(ngrid),ra(ngrid)
      dimension rinvn(ngrid),rinvp(ngrid),rinva(ngrid)
      rmin=1.e-99
      open(unit=10,file='tgridc')
      read(10,*)
      read(10,*)
      read(10,*) (t8(i),i=1,ngrid)
      close(10)
      qn=-2.598
      qp=2.240
      qa=4.617
      g12=1.
      gn=2.
      gp=2.
      ga=1.
      gmg23=4.
      gna23=4.
      do 50 i=1,ngrid
      t9=t8(i)/10.
      t9a=t9/(1.+0.0396*t9)
      r=4.27e26*t9a**(5./6.)/t9**1.5*
     &	 exp(-84.165/t9a**(1./3.)-2.12e-3*t9**3)
      yn=0
      yp=0.44
      ya=0.56
      if (t9.lt.3.3.and.t9.ge.1.75) then
	yn=0.05
	yp=0.45
	ya=0.50
      endif
      if (t9.ge.3.3) then
	yn=0.07
	yp=0.40
	ya=0.53
      endif
      rn(i)=r*yn
      rp(i)=r*yp
      ra(i)=r*ya
      rinvn(i)=rn(i)/gn/gmg23*(12.*12./23.)**1.5/2.*texp(-11.605*qn/t9)
      rinvp(i)=rp(i)/gp/gna23*(12.*12./23.)**1.5/2.*texp(-11.605*qp/t9)
      rinva(i)=ra(i)/ga*(12.*12./4./20.)**1.5/2.*texp(-11.605*qa/t9)
      if (rn(i).le.rmin) rn(i)=rmin
      if (rp(i).le.rmin) rp(i)=rmin
      if (ra(i).le.rmin) ra(i)=rmin
      if (rinvn(i).le.rmin) rinvn(i)=rmin
      if (rinvp(i).le.rmin) rinvp(i)=rmin
      if (rinva(i).le.rmin) rinva(i)=rmin
   50 continue
      open(unit=11,file='c12c12.out')
      write(11,999) (t8(i),i=1,ngrid)
      write(11,*) '2 C12 --> MG23 + neut: Q=-2.598 MeV'
      write(11,999) (log10(rn(i)),i=1,ngrid)
      write(11,*) 'MG23 + neut --> 2 C12: Q= 2.598 MeV'
      write(11,999) (log10(rinvn(i)),i=1,ngrid)
      write(11,*) '2 C12 --> NA23 + prot: Q= 2.240 MeV'
      write(11,999) (log10(rp(i)),i=1,ngrid)
      write(11,*) 'NA23 + prot --> 2 C12: Q=-2.240 MeV'
      write(11,999) (log10(rinvp(i)),i=1,ngrid)
      write(11,*) '2 C12 --> NE20 + alph: Q= 4.617 MeV'
      write(11,999) (log10(ra(i)),i=1,ngrid)
      write(11,*) 'NE20 + alph --> 2 C12: Q=-4.617 MeV'
      write(11,999) (log10(rinva(i)),i=1,ngrid)
  999 format(6e13.6)
      close(11)
      end
c
      function texp(x)
      implicit double precision (a-h,o-z)
      texp=0.
      if (x.gt.200..or.texp.lt.-200.) return
      texp=dexp(x)
      return
      end
