
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:36:36 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 7                                                           $ *
*                                                                      *
      character cphase*4,phase*6

      integer nphase
      integer irotbin,model,modeli,model_old
      integer imodbin,maxmod,maxmod0,imodpr,imodpr0
      integer no
      integer ivisc,nq0
      integer icrash,icrash0,iresmax,ireset,ifail,ireverse,numeric,
     &     nretry

      double precision alphmin,alphmax,alphmax0
      double precision epsmax,code_version,bin_version
      double precision timebeg,timeend,cputime
      double precision phi,phi0,alpha,alpha0,sigma,dynfac
      double precision dtmin,dtmax,dtmax0,q0,mcut,menv
      double precision facdt,fkhdt,shlim,fts,ftsh,ftshe,ftsc,ftacc,ftst,
     &     fts0,facdt0,fkhdt0,ftnuc

      logical hydro,hydrorot,iacc,dup3,agbphase,relaxpulse,thermalpulse,
     &     rgbphase,superagb,flame,urca,iacc0,hydro0,iatm
      logical astero
      character*1 icorr,icorr0

      common /accel/ alphmin,alphmax,alphmax0,sigma,dynfac,iacc,hydro,
     &     hydrorot,iacc0,hydro0
      common /asterosismo/ astero
      common /converg/ epsmax(neq)
      common /viscosity/ q0,mcut,menv,ivisc,nq0
      common /cpu/ timebeg,timeend,cputime,code_version,bin_version
      common /evopha/ nphase,dup3,agbphase,relaxpulse,thermalpulse,
     &     rgbphase,superagb,flame,urca
      common /evolphase/ phase(9),cphase(9)
      common /outp/ irotbin,model,modeli,model_old
      common /parmod/ phi,phi0,alpha,alpha0,imodbin,maxmod,imodpr,
     &     imodpr0,maxmod0
      common /print/ no
      common /resetopt/ icrash,icrash0,iresmax,ireset,ifail,ireverse,
     &	   nretry,numeric,iatm
      common /autocorrect/ icorr,icorr0
      common /timed/ dtmin,dtmax,facdt,fkhdt,shlim,fts,ftsh,ftshe,
     &     ftsc,ftacc,ftst,ftnuc,fts0,facdt0,fkhdt0,dtmax0
