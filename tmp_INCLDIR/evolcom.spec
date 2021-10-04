
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:36:36 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 7                                                           $ *
*                                                                      *

      character nuclopt*1,nuclopt0*1

      double precision vxsp
      double precision enunucl,enupla,flu2
      double precision mtotlos
      double precision flu,xsp,ysp,tburn,tnucmin,fnetdt,
     &     fnetdt0,tolnuc,tolnuc0

      common /chinit/ vxsp(nsh,nsp)
      common /nue/ enunucl(nsh),enupla(nsh),flu2(nreac)
      common /xsplos/ mtotlos(nsp)
      common /xspnuc/ flu(nsh,nreac),xsp(nsh,nsp),ysp(nsh,nsp)
      common /nuc_param/ tburn,tnucmin,fnetdt,fnetdt0,tolnuc,tolnuc0
      common /nuc_c_param/ nuclopt,nuclopt0
