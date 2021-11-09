
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:36:36 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 7                                                           $ *
*                                                                      *
      double precision tlhkap
      double precision kap,kapm
      double precision dkapdt,dkapdro,dkapdf,dkapdx,dkapdy

      common /kaptab/ tlhkap
      common /opa/ kap(nsh),kapm(nsh)
      common /opax/ dkapdt(nsh),dkapdro(nsh),dkapdf(nsh),
     &     dkapdx(nsh),dkapdy(nsh)

