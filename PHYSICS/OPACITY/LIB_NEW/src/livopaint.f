***---------------------------------------------------------------------------
*
      subroutine opac(z,xh,xxci,xxoi,t6,r,error)
*     ==========================================
*
*..... The purpose of this subroutine is to interpolate log kappa
*          (and obtain smooth derivatives)
*      in hydrogen (if X>0) and in C/O abundance and T6, R, i.e. (X,Xc,Xo,T6,R)
*
*  It merely takes the logs of T6 and R, then calls OPAL to do the actual work.
*
*      z = Z = metallicity (this should always be the same)
*      xh = X = hydrogen mass fraction
*      xxci = Xc = carbon mass fraction (excess over what is in Z)
*      xxoi = Xo = oxygen mass fraction (excess over what is in Z)
*      t6 = T6 = temperature in millions of degrees kelvin
*      r = R = density(g/cm**3)/T6**3
*===
* Version : may 2003 (from April 2001 updated OPAL routine z14xcotrin21.f)
*===
***---------------------------------------------------------------------------

      implicit double precision (a-h, o-z)

      integer error


      if ( t6 .le. 0.d0 .or. r .le. 0.d0 ) then
         write(6,8437) t6,r
 8437    format(' '/' STOP -- OPAC: non-positive value of T6=',
     $        1p,e11.3,' or R=',e11.3)
         error = 8
         return
      endif
c
      slt = log10(t6)
      slr = log10(r)
c
      
      call opal(z,xh,xxci,xxoi,slt,slr,error)
c
      return
      end


*
***-----------------------------------------------------------------------------

      subroutine opal(z,xh,xxci,xxoi,slt,slr,error)
*     =======================================
*
*..... The purpose of this subroutine is to interpolate log kappa
*          (and obtain smooth derivatives)
*      in hydrogen (if X>0) and in C/O abundance and T6, R, i.e. (X,Xc,Xo,T6,R)
*
*      z = Z = metallicity (this should always be the same)
*      xh = X = hydrogen mass fraction
*      xxci = Xc = carbon mass fraction (excess over what is in Z)
*      xxoi = Xo = oxygen mass fraction (excess over what is in Z)
*      slt = logT6 = Log10{temperature in millions of degrees kelvin}
*      slr = logR = Log10{density(g/cm**3)/T6**3}
*
*..... to use OPAC or OPAL insert common/extopac/ in the calling routine.
*      This common contains interpolated values for log kappa and its 
*      first derivatives, and "out-of-table" indicators fedge,ftredge,fzedge.
*
***-----------------------------------------------------------------------------
      implicit double precision (a-h, o-z)

      integer error

c PARAMETERS to specify opacity storage matrices: see OPALINIT
c

c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,

c!      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=28, nrb=1, nre=28,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
c PARAMETERS nrdel and ntdel give matrix position differences: see OPALINIT
c
      parameter ( nrdel=nrb-1, ntdel=ntb-1 )
c
c PARAMETERS: offsets for Z, for X, min low-T X offset: see OPALINIT
c
      parameter ( zdel=0.001d0, xdel=0.03d0, xdelmin=0.001d0 )
c
c PARAMETERS: for use in call to READZEXCO, if it is called from here in OPAL:
c  iu_low = 50 = lowest FORTRAN I/O unit number to be used in reading opacities
c                 (units 50 through 56 may be used)
c  k_hz = 1 = khighz value for READZEXCO ("use 'GN93hz' in Z-interpolation")
c  ofe_brack = 0.0 = [O/Fe] value for READZEXCO
c
      parameter ( iu_low=50, k_hz=1, ofe_brack=0.0d0 )
c
c PARAMETER badlogkval = 1.e+35 is stored to indicate missing Log(kappa) values
c
      parameter ( badlogkval=1.d+35 )
c
c PARAMETERS used during tests for high-T boundary of opacity storage matrices
c
      parameter ( ntm_m3=ntm-3, nt_m2=nt-2 )
c
c PARAMETERS used to specify logT6 tabulation values
c
c!!      parameter ( k81=nt-3, k80=k81-1, k60=k81-21, ks59=k60-1+ntdel )
      parameter ( k81=nt-9, k80=k81-1, k60=k81-21, ks59=k60-1+ntdel )
      parameter ( flt81m6=8.1d0-6.d0, flt60m6=6.0d0-6.d0, 
     &     flt370m6=3.70d0-6.d0 )
c
c PARAMETERS defining the storage for the additional X-values from 'GN93hz':
c
      parameter ( mx_hi=2*mx, mo_m1=mo-1, mx_hi_nz=mx_hi*nz )
c
c COMMON /xhi_opal_z/ : auxiliary matrices for additional 'GN93hz' X-values:
c
      common /xhi_opal_z/ xhi_in(mx_hi), xhi_use(mx_hi,nz),
     $     xxx_hi(mx_hi), nx_hi(nz), ireq_hi(mx_hi), khighx(nz),
     $     kavail_xhi, kuse_xhi, kdo_xhi
      save /xhi_opal_z/
c
c COMMON /a_opal_z/ : matrices for opacity storage: see OPALINIT
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime(mx),cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
c COMMON /b_opal_z/ : high and low temperature limits, Z-values: see OPALINIT
c
      common/b_opal_z/ nta(0:nrm_p1),zz(mx,nz),ntax0(0:nrm),
     $     ntax03(0:nrm), sltlo, slthi, dltlo_inv, dlthi_inv,
     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
      save /b_opal_z/
c
c COMMON /bb_opal_z/ : some indices & abundances for T6,R and C,O interpolation
c
      common/bb_opal_z/ l1,l2,l3,l4,k1,k2,k3,k4,ip,iq(4),xodp,xcdp,
     $     xxco,cxx,oxx,kzf,kzg,kzh,kzf2
      save /bb_opal_z/
c
c COMMON /recoin_opal_z/ : see OPALINIT
c
      common/recoin_opal_z/ itimeco,mxzero,mx03,kope,igznotgx
      save /recoin_opal_z/
c
c COMMON /c_opal_ctrl_smooth/ : flags to control the smoothing:
c
      common/c_opal_ctrl_smooth/ init_smo, low_CO_smo, interp_CO_smo
      save /c_opal_ctrl_smooth/
c
c COMMON /c_level_err_opal_z/  error-checking level, set by SET_ERR_CHECK
c
c /c_level_err_opal_z/: --> data{level_err}
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c
c COMMON /EXTOPAC/ : Return variables : see also instructions above
c
c.... OPACT - opacity obtained from a quadratic interpolation at fixed
c      log T6 at three values of log R; followed by quadratic interpolation
c      along log T6.  Results smoothed by mixing overlapping quadratics.
c.... DOPACT - is Dlog(k)/Dlog(T6) at constant R smoothed by mixing quadratics.
c.... DOPACR - is  Dlog(k)/Dlog(R) at constant T smoothed by mixing quadratics.
c.... DOPACTD - is Dlog(k)/Dlog(T6) at constant rho.
c.... FEDGE = 1.0 inside T6,R,Z boundaries, goes to zero when too far outside
c      for extrapolation (in which case opacity is not calculated).
c.... FTREDGE = 1.0 inside T,R boundaries, goes to zero when too far outside.
c.... FZEDGE = 1.0 inside Z limits, goes to zero when too far outside.
c
      common/extopac/ opvals(4),fedge,ftredge,fzedge
      save /extopac/
      equivalence (opvals(1),opact), (opvals(2),dopact),
     $     (opvals(3),dopacr), (opvals(4),dopactd)
c
c-debug[
c-debug;      common/outdeb/ ioudeb,oudebl,koudeb
c-debug]
c-test-xdel[              ! test various xdel values in opacity X-interpolation
c-test-xdel;
c-test-xdel;      parameter ( n_xdel_test=9, n_xdel_test_m1=n_xdel_test-1 )
c-test-xdel;      common/otherxdel/ opdxi(n_xdel_test),xdel_test(n_xdel_test)
c-test-xdel;      equivalence (opdxi(1),opdx001), (opdxi(2),opdx002),
c-test-xdel;     $      (opdxi(3),opdx0035), (opdxi(4),opdx005),
c-test-xdel;     $      (opdxi(5),opdx0075), (opdxi(6),opdx01),
c-test-xdel;     $      (opdxi(7),opdx02), (opdxi(n_xdel_test_m1),opdx03),
c-test-xdel;     $      (opdxi(n_xdel_test),opdxuse)
c-test-xdel;      data xdel_test/.001,.002,.0035,.005,.0075,.01,.02,.03,.03/
c-test-xdel]
c___
      dimension iqx(4),ntaxat(0:nrm)
c===					! we have not yet gotten a good opacity
      opact = badlogkval
      fedge = 0.d0
      ftredge = 0.d0
      fzedge = 0.d0
c
c..... set-up C/O axis points: making xxci & xxoi arguments of OPAL is simpler
c      than resetting xxc & xxo after changing them, and also avoids an error
c      if constants are used for these values in the calling program.
c
      xxc = xxci
      xxo = xxoi
      xxco = xxc + xxo
      if ( z+xh+xxco-1.d0-6 .gt. 1.0d0 .or. z .lt. -1.d-8 .or.
     $     min(xh,xxc+z,xxo+z,xxco+z) .le. -1.d-6 ) then
         write(6,4397) z,xh,xxc,xxo
 4397    format(' '/' STOP -- OPAL: bad value(s) Z',f10.7,' X',f10.7,
     $        ' C',f10.7,' O',f10.7)
         error = 8
         return
      endif

c
c..... set X indices: if xh = a table X-value, then only use that X-table:
c
      if ( abs( xh - xa(1) ) .lt. 9.d-7 ) then
         mf = 1
         mg = 1
         mh = 1
         mf2 = 1
      else if ( mx .eq. 1 ) then
         write(6,4396) xh,xa(1)
 4396    format(' '/' STOP -- OPAL: mx=1, but X=',f10.7,
     $        ' differs from table value',f10.7)
         stop
c		    ! else: find X-indices: results in ihi=5 (for xh > 0.35) or
c		    !        ihi=4 (for xh > 0.1) or ihi=3 (for xh < or = 0.1):
      else
         ilo = 2
         ihi = mx
         do while( ihi-ilo .gt. 1 )
            imd = (ihi+ilo)/2
            if ( xh .le. xa(imd) ) then
               ihi = imd
            else
               ilo = imd
            endif
         enddo
c						    ! if xh = a table X-value
         if ( abs( xh - xa(ilo) ) .lt. 9.d-7 ) then
            mf = ilo
            mg = ilo
            mh = ilo
            mf2 = ilo
         else if ( abs( xh - xa(ihi) ) .lt. 9.d-7 ) then
            mf = ihi
            mg = ihi
            mh = ihi
            mf2 = ihi
c				! else: will need to interpolate
         else
            mf = ilo - 1
            mg = ilo
            mh = ihi
            if ( xh .le. xa(2) ) then
               mf2 = mh
            else
               mf2 = min( ihi + 1 , mx )
            endif
         endif
      endif
c
c  If X > 0.03 and C+O is not too large, then the more-numerous X-values from
c  'GN93hz' can help, if use-flag is set: set f_xhi for later use.  Only used
c  if C+O < 0.3 and ( kdo_xhi = 1 and X > .7 ) or ( kdo_xhi = 2 and X > .03 );
c  FULL opacity shifts are only used if C+O < 0.2, partial for 0.2 < C+O < 0.3:
c
      if ( kdo_xhi .le. 0 .or. mf2 .lt. 4 .or. mf .gt. 3 .or.
     $     ( kdo_xhi .eq. 1 .and. xh .le. xa(mx) ) ) then
         f_xhi = -1.d0
      else
         f_xhi = 1.d0 - 100.d0 * max( xxci + xxoi - 0.2d0 , 0.0d0 )**2
      endif
c
c..... If necessary, read data files and do initializations, using READZEXCO
c
      ircall = 0
      do m = mf,mf2
         if ( itime(m) .ne. 12345678 ) ircall = 1
      enddo

      if ( ircall .ne. 0 )
     $     call readzexco(nz,-1.d0,z,-1.d0,k_hz,iu_low,ofe_brack)
c
c..... set Z indices
c							! is Z out of range?
      if ( z .le. zlo_ex .or. z .ge. zhi_ex ) then
         if ( level_err .ge. 2 ) then
            write(6,10) z, zlo_ex, zhi_ex
 10         format(' '/' OPAL: Z=',f11.8,
     $           ' is outside max extrap range (',f11.8,',',f10.8,')')
            stop ' STOP -- OPAL: bad Z value. '
         endif
         return
      endif
c							! check Z-extrapolation
      if ( z .lt. zlow ) then
         fzedge = max( 0.0d0 , ( z - zlo_ex ) / ( zlow - zlo_ex ) )
      else if ( z .gt. zhigh ) then
         fzedge = max( 0.0d0 , ( zhi_ex - z ) / ( zhi_ex - zhigh ) )
      else
         fzedge = 1.d0
      endif

c				   ! this shouldn't happen, but just in case...
      if ( fzedge .le. 0.d0 ) then
         if ( level_err .ge. 2 ) then
            write(6,10) z, zlo_ex, zhi_ex
            stop ' STOP -- OPAL: bad Z value. '
         endif
         return
      endif
c						! check for Z-table value:
      if ( numz .eq. 1 ) then
         ihi = 1
      else
         ihi = 0
         if ( numz .le. 3 ) then
            do i = 1, numz
               if ( abs( zsto(i) - z ) .le. zacc(i) ) ihi = i
            enddo
         else if ( abs( zsto(1) - z ) .le. zacc(1) ) then
            ihi = 1
         endif
      endif
c						! get Z-table indices:
      if ( ihi .gt. 0 ) then
         kzf = ihi
         kzg = ihi
         kzh = ihi
         kzf2 = ihi
      else if ( numz .le. 3 ) then
         kzf = 1
         kzg = 2
         kzh = numz
         kzf2 = numz
      else
         ilo = 2
         ihi = numz
         do while( ihi-ilo .gt. 1 )
            imd = (ihi+ilo)/2
            if ( z .le. zsto(imd) ) then
               ihi = imd
            else
               ilo = imd
            endif
         enddo
         if ( abs( zsto(ihi) - z ) .le. zacc(ihi) ) then
            kzf = ihi
            kzg = ihi
            kzh = ihi
            kzf2 = ihi
         else if ( abs( zsto(ilo) - z ) .le. zacc(ilo) ) then
            kzf = ilo
            kzg = ilo
            kzh = ilo
            kzf2 = ilo
         else
            kzf = ilo - 1
            kzg = ilo
            kzh = ihi
            if ( z .le. zsto(2) ) then
               kzf2 = kzh
            else
               kzf2 = min( ihi + 1 , numz )
            endif
         endif
      endif
c		! note that xxh is not used (except perhaps in calling routine)
      xxh = xh
c							     ! check T-R edges:
      ftredge = max( 0.d0 , min(slr-slrlo,0.d0)*dlrlo_inv + 1.d0 )
     $     * max( 0.d0 , min(slrhi-slr,0.d0)*dlrhi_inv + 1.d0 )
     $     * max( 0.d0 , min(slt-sltlo,0.d0)*dltlo_inv + 1.d0 )
     $     * max( 0.d0 , min(slthi-slt,0.d0)*dlthi_inv + 1.d0 )
      fedge = max( 0.d0 , ftredge * fzedge )
c						   ! if too far outside, return
      if ( fedge .le. 0.d0 ) then
         if ( level_err .ge. 2 ) then
            write(6,20) slt+6.d0, slr
 20         format(' '/' OPAL: logT=',f9.6,', logRHO=',f11.6,
     $           ' lies outside the opacity matrix')
            stop ' STOP -- OPAL: bad T or RHO value. '
         endif
         return
      endif
c
c..... convert xh to logarithmic shifted by xdel (C and O are converted later)
c
      xxx = log10(xdel+xh)
c
c..... Determine log R and log T6 grid points to use in the interpolation.
c
      if ( slt .ge. flt81m6 ) then
         k2sat = min(int((slt-flt81m6)*dfs(nt))+k81,nt)+ntdel
      else if ( slt .ge. flt60m6 ) then
         k2sat = min(int((slt-flt60m6)*dfs(k81))+k60,k80)+ntdel
      else if ( k60 .le. 0 ) then
         k2sat = k60+ntdel
      else
         k2sat = min(max(int((slt-flt370m6)*dfs(k60)),ntdel),ks59)
      endif
c
      l2sat = max(min(int((slr-alr(1))*dfsr(nr)+1.d0),nr),0)+nrdel
c
      k1x = -99
      l1x = -99
      k3sat = k2sat+1
      k2 = max(k2sat-ntdel,1)
      k3 = min(k3sat-ntdel,nt)
      l3sat = l2sat+1
      l4sat = min(l3sat,nrm)+1
      l2 = max(l2sat-nrdel,1)
      l3 = min(l3sat-nrdel,nr)
      if ( min(k3,l3) .le. 0 .or. k2 .gt. nt .or. l2 .gt. nr ) then
         ftredge = 0.d0
         fedge = 0.d0
         if ( level_err .ge. 2 ) then
            write(6,20) slt+6.d0, slr
            stop ' STOP -- OPAL: bad T or RHO value. '
         endif
         return
      endif
c		! initial assumption: 3x3 grid in T,R
      ipx = 2
      do i = 1,4
         iqx(i) = 2
      enddo
c					! Check upper-right ragged cut-out:
      if ( l2sat .gt. 0 ) then
c							. . .
c						1.	. .      too far out to
c							. .  *   extrapolate
         if ( k2sat .gt. nta(l2sat) ) then
            ftredge = 0.d0
            fedge = 0.d0
c							. . .    extrapolate
c						2.	. . .    T and R out
c							. .  *   from a corner
         else if ( k2sat .eq. nta(l2sat) .and.
     $           k2sat .gt. nta(l3sat) ) then
            ft = max(min(alt(k2)-slt,0.d0)*dfs(k3)+1.d0,0.d0)
            fr = max(min(alr(l2)-slr,0.d0)*dfsr(l3)+1.d0,0.d0)
            if ( l2 .lt. nr ) ftredge = ftredge*fr
            if ( k2 .lt. nt ) ftredge = ftredge*ft
            k1x = k2-2
            l1x = l2-2
c						. . .	. . .   extrapolate (in
c					3.	. .*	. . .   either T or R)
c						. .	. .*    in a corner
         else if ( k2sat .lt. nta(l2sat) .and.
     $           k2sat .eq. nta(l3sat) ) then
            ft = max(min(alt(k2)-slt,0.d0)*dfs(k3)+1.d0,0.d0)
            fr = max(min(alr(l2)-slr,0.d0)*dfsr(l3)+1.d0,0.d0)
            if ( ft .gt. fr ) then
               ftredge = ftredge*ft
               k1x = k2-2
               l1x = l2-1
               do i = 1,3
                  if ( k2sat-3+i .le. nta(l4sat) ) iqx(i) = 3
               enddo
            else
               ftredge = ftredge*fr
               k1x = k2-1
               l1x = l2-2
               if ( k3sat .lt. nta(l2sat) ) ipx = 3
            endif
c						. . .	. . .     extrapolate R
c					4.	. .	. . .*    out from a
c						. .*	. . .*    high-R edge
         else if ( k2sat .lt. nta(l2sat) .and.
     $           k2sat .gt. nta(l3sat) ) then
            if ( l2 .lt. nr ) ftredge = 
     $           ftredge*max(min(alr(l2)-slr,0.d0)*dfsr(l3)+1.d0,0.d0)
            if ( k3sat .lt. nta(l2sat) .and. k2sat .gt. ntb ) ipx = 3
            k1x = max(k2-1,1)
            l1x = l2-2
c						. . .	. . . .  interpolate,
c					5.	. .*.	. .*.*.  1-space inside
c						.*.*.	.*.*.    high-R,T edges
         else if ( k3sat .le. nta(l3sat) .and.
     $           k3sat .ge. nta(l4sat) ) then
            if ( k3sat .lt. nta(l3sat) .and. k2sat .gt. ntb ) ipx = 3
            k1x = max(k2-1,1)
            l1x = max(l2-1,1)
            if ( l2sat .gt. nrb ) then
               do i = 1,3
                  if ( k1x-1+ntdel+i .le. nta(l4sat) ) iqx(i) = 3
               enddo
            endif
         endif
      endif
c
      fedge = max( ftredge * fzedge , 0.d0 )
c
c					! if too far out to extrapolate, return
      if ( fedge .le. 0.d0 ) then
         if ( level_err .ge. 2 ) then
            write(6,20) slt+6.d0, slr
            stop ' STOP -- OPAL: bad T or RHO value. '
         endif
         return
      endif
c							 . . . .   Outside or
c							 . . . .   1-space
c						6.	*. . .     inside max
c							* * *      high-T edge
      if ( k1x .lt. 0 .and. k3sat .ge. ntm ) then
         k1x = nt_m2
         l1x = max(l2-1,1)
         if ( l2sat .gt. nrb ) then
            do i = 1,3
               if ( ntm_m3+i .le. nta(l4sat) ) iqx(i) = 3
            enddo
         endif
c				   7.   Anywhere except high-T or high-R edges:
      else if ( k1x .lt. 0 ) then
         k1x = max(k2-1,1)
         l1x = max(l2-1,1)
         if ( k2sat .gt. ntb ) ipx = 3
         if ( l2sat .gt. nrb ) then
            do i = 1,4
               iqx(i) = 3
            enddo
         endif
      endif
c		       ! check low-T,low-R corner for X=0; avoid it if possible
      ichgr = 1
      if ( mf .eq. mxzero .and. k1x+ntdel .lt. ntax0(l1x+nrdel) ) then
c									! avoid
         if ( mf2 .eq. mf+3 ) then
            mf = mf2-2
            mg = mf2-1
            mh = mf2
c						   ! like region 1. too far out
         else if ( k3sat .lt. ntax0(l3sat) ) then
            ftredge = 0.d0
            fedge = 0.d0
            if ( level_err .ge. 2 ) then
               write(6,20) slt+6.d0, slr
               stop ' STOP -- OPAL: bad T or RHO value. '
            endif
            return
c		     ! else, will need to revise T,R indices for first m (= mf)
         else
            ichgr = 0
         endif
      endif
c						    ! check similarly for X=.03
      if ( ( mf .eq. mx03 .or. mg .eq. mx03 ) .and.
     $     k1x+ntdel .lt. ntax03(l1x+nrdel) ) then
c									! avoid
         if ( mf2 .eq. mf+3 .and. mf .eq. mx03 ) then
            mf = mf2-2
            mg = mf2-1
            mh = mf2
c						   ! like region 1. too far out
         else if ( k3sat .lt. ntax03(l3sat) ) then
            ftredge = 0.d0
            fedge = 0.d0
            if ( level_err .ge. 2 ) then
               write(6,20) slt+6.d0, slr
               stop ' STOP -- OPAL: bad T or RHO value. '
            endif
            return
c					! if need to revise T,R indices for mf:
         else if ( mf .eq. mx03 ) then
            ichgr = 0
         endif
      endif
c	     ! xhemxi: subtract Z to prevent out-of-range C+O values at small X
c
      xhemxi = max( 1.d0 - xh - z , 0.d0 )
      ftrbeg = ftredge
      if ( kzf2 .gt. kzf ) then
         zlogd = log10( z + zdel )
      else
         zlogd = 0.0d0
      endif
c				! ------------------------- loop over X-mixes m
      do m = mf,mf2
c				  ! set (or restore) grid-indices, if necessary
         if ( ichgr .ne. 0 ) then
            do i = 1,4
               iq(i) = iqx(i)
            enddo
            ip = ipx
            k1 = k1x
            l1 = l1x
            k2 = k1+1
            k3 = k1+2
            k4 = k1+3
            l2 = l1+1
            l3 = l1+2
            l4 = l1+3
         endif
c				! check for low-T,low-R cutout at X=0 or X=.03:
         ichgr = 0
         if ( m .eq. mxzero .and.
     $        k1x+ntdel .lt. ntax0(l1x+nrdel) ) then
            ichgr = 1
            do i = max(l2sat-2,0),min(l4sat+2,nrm)
               ntaxat(i) = ntax0(i)
            enddo
         else if ( m .eq. mx03 .and.
     $           k1x+ntdel .lt. ntax03(l1x+nrdel) ) then
            ichgr = 1
            do i = max(l2sat-2,0),min(l4sat+2,nrm)
               ntaxat(i) = ntax03(i)
            enddo
         endif
c				  ! change grid indices, if in low-{R,T} cutout
         if ( ichgr .ne. 0 ) then
            k3 = min(k3sat-ntdel,nt)
            l3 = min(l3sat-nrdel,nr)
            l1sat = max(l2sat-1,0)
            ip = 2
            do i = 1,4
               iq(i) = 2
            enddo
            ftrprev = ftredge
            ftredge = ftrbeg
c					       ! 2. extrapolate T,R from corner
            if ( k3sat .eq. ntaxat(l3sat) .and.
     $           k3sat .lt. ntaxat(l2sat) ) then
               if ( l4sat .ge. nre ) then
                  ftredge = 0.d0
               else
                  ft = max(min(slt-alt(k3),0.d0)*dfs(k3)+1.d0,0.d0)
                  fr = max(min(slr-alr(l3),0.d0)*dfsr(l3)+1.d0,0.d0)
                  if ( l3 .gt. 1 ) ftredge = ftredge*fr
                  if ( k3 .gt. 1 ) ftredge = ftredge*ft
                  k1 = k3
                  l1 = l3
               endif
c						   ! 3. extrap T or R in corner
            else if ( k3sat .gt. ntaxat(l3sat) .and.
     $              k3sat .eq. ntaxat(l2sat) ) then
               ft = max(min(slt-alt(k3),0.d0)*dfs(k3)+1.d0,0.d0)
               fr = max(min(slr-alr(l3),0.d0)*dfsr(l3)+1.d0,0.d0)
               if ( ft .gt. fr .or. l4sat .ge. nre ) then
                  ftredge = ftredge*ft
                  k1 = k3
                  l1 = l3-2
                  if ( k3sat .lt. ntaxat(l1+nrdel) ) then
                     l1 = l3-1
                     if ( l3 .ge. nr ) ftredge = 0.d0
                  else if ( l3 .lt. nr ) then
                     do i = 1,3
                        iq(i) = 3
                     enddo
                  endif
               else
                  ftredge = ftredge*fr
                  k1 = max(ntaxat(l3sat)-ntdel,k3-2)
                  l1 = l3
                  if ( k1 .eq. k3-2 ) ip = 3
               endif
c						     ! 4. extrapolate R
            else if ( k3sat .gt. ntaxat(l3sat) .and.
     $              k3sat .lt. ntaxat(l2sat) ) then
               if ( l4sat .ge. nre ) then
                  ftredge = 0.d0
               else
                  k1 = max(ntaxat(l3sat)-ntdel,k3-2)
                  l1 = l3
                  if ( k1 .eq. k3-2 ) ip = 3
                  if ( l3 .gt. 1 ) ftredge = 
     $                 ftredge*max(min(slr-alr(l3),0.d0)*dfsr(l3)+1.d0,
     $                 0.d0)
               endif
c						     ! 8. extrapolate T
            else if ( k3sat .eq. ntaxat(l3sat) .and.
     $              k3sat .eq. ntaxat(l2sat) ) then
               k1 = k3
               l1 = l3-2
               if ( k3sat .lt. ntaxat(l1+nrdel) ) then
                  l1 = l3-1
                  if ( l3 .ge. nr ) ftredge = 0.d0
               else if ( l3 .lt. nr ) then
                  do i = 1,3
                     iq(i) = 3
                  enddo
               endif
               if ( k3 .gt. 1 ) ftredge = ftredge
     $              *max(min(slt-alt(k3),0.d0)*dfs(k3)+1.d0,0.d0)
c							    ! 5. inside an edge
            else if ( k2sat .ge. ntaxat(l2sat) .and.
     $              k2sat .le. ntaxat(l1sat) ) then
               if ( k2sat .eq. ntaxat(l2sat) .or. l3 .eq. nr ) then
                  l1 = l3-2
                  k1 = k3-1
                  if ( k2sat .lt. ntaxat(l1+nrdel) ) then
                     if ( l3 .lt. nr ) then
                        l1 = l3-1
                     else if ( k3sat .lt. ntaxat(l1+nrdel) ) then
                        ftredge = 0.d0
                     else
                        k1 = k3
                        ftredge = ftredge
     $                       *max(min(slt-alt(k3),0.d0)*dfs(k3)+1.d0,
     $                       0.d0)
                     endif
                  else if ( l3 .lt. nr ) then
                     do i = 1,3
                        iq(i) = 3
                     enddo
                  endif
               else
                  k1 = max(k3-2,1)
                  l1 = l3-1
                  if ( k1 .eq. k3-2 ) ip = 3
               endif
            endif
c					! if too far out to extrapolate, return
            if ( ftredge .le. 0.d0 ) then
               fedge = 0.d0
               opact = badlogkval
               if ( level_err .ge. 2 ) then
                  write(6,20) slt+6.d0, slr
                  stop ' STOP -- OPAL: bad T or RHO value. '
               endif
               return
            endif
c					   ! use smaller of X=0 or X=.03 values
            ftredge = min(ftredge,ftrprev)
            fedge = max( ftredge * fzedge , 0.d0 )
c
c					   ! get rest of revised grid indices
            k2 = k1+1
            k3 = k1+2
            k4 = k1+3
            l2 = l1+1
            l3 = l1+2
            l4 = l1+3
c			! end of revised grid in low-T,R cutout
         endif
c			     ! ------------------------- loop over Z-values kz:
         do kz = kzf, kzf2
c				! xhemx: subtract Z to prevent out-of-range C+O
c						! values at small X
            xhemx = 1.d0 - xa(m) - zsto(kz)
c							! If no X or Z interp
            if ( kzf2 .eq. kzf .and. mf2 .eq. mf ) then
c
               xxc = xxci
               xxo = xxoi
c				! Else, if we will be interpolating in X or Z:
            else
c
c............. C and O  fractions determined by the ray through the origin that
c              also passes through the point (Xc,Xo). Specific interpolation 
c              values are determined by tabulated X values; i.e., xa(m).
c              Interpolation along the ray gives log(kappa(Xc,Xo)).
c              (Advantage of method: keeps indices within table boundaries.)
c
               if ( xhemxi .gt. 1.d-6 ) then
                  cmod = xhemx / xhemxi
               else
                  cmod = 0.d0
               endif
               if ( xxci .gt. 0.d0 ) then
                  xxc = cmod * xxci
               else if ( xxci .ge. -1.d-8 .or. z .lt. 1.d-8 ) then
                  xxc = 0.d0
               else
                  xxc = max( xxci / z , -1.d0 ) * zsto(kz)
               endif
               if ( xxoi .gt. 0.d0 ) then
                  xxo = cmod * xxoi
               else if ( xxoi .ge. -1.d-8 .or. z .lt. 1.d-8 ) then
                  xxo = 0.d0
               else
                  xxo = max( xxoi / z , -1.d0 ) * zsto(kz)
               endif
c
            endif
c
            xxco = xxc + xxo
c
c..... convert xxc and xxo to logarithmic shifted by Z+zdel
c
            cxx = log10(zzz(kz)+xxc)
            oxx = log10(zzz(kz)+xxo)
c				    ! set up table C,O abundances for this m,kz
            nc = n(m,1,kz)
            no = nc
            do i = 1,nc-1
               xc(i) = xcs(i)
               xo(i) = xos(i)
            enddo
            xc(nc) = xhemx
            xo(nc) = xhemx
c
            do i = 1,nc
               ox(i) = oxf(m,i,kz)
               cx(i) = cxf(m,i,kz)
               xcd(i) = xcdf(m,i,kz)
               xod(i) = xodf(m,i,kz)
               cxd(i) = cxdf(m,i,kz)
               oxd(i) = oxdf(m,i,kz)
            enddo
c
            xodp = max(-xxc+xc(nc),0.d0)
            xcdp = max(-xxo+xo(no),0.d0)
c-debug[
c-debug;            if ( m .eq. mf .and. kz .eq. kzf .and. ioudeb .gt. 0 )
c-debug;     $           write(6,9409) z,xh,xxci,xxoi,10.**slt,slt,slr,
c-debug;     $           xxc,xxo,mf,mf2,kzf,kzf2,k1,ip,l1,(iq(i),i=1,ip+1)
c-debug; 9409       format(' '/' OPAL: Z',f10.7,' X',f10.7,' C',f10.7,
c-debug;     $           ' O',f10.7,' T6',f12.7,' logT6',f11.7,' logR',f11.7,
c-debug;     $           ' xxc',f10.7,' xxo',f10.7,' m',2i2,' kz',2i3,
c-debug;     $           ' k',i3,'+',i1,' l',i3,'+',4i1)
c-debug]
c
c  Interpolate in C and O: COINTSMO is better, and was more thoroughly tested:
c
            if ( interp_CO_smo .gt. 0 ) then
               call cointsmo(xxc,xxo,kz)
            else
               call cointerp(xxc,xxo,kz)
            endif
c				! ---------------- end of loop over Z-values kz
         enddo
c
c....... Interpolate in Z, if necessary, mixing overlapping quadratics.
c
         call qzlog4int( zlogd )
c
c....... completed C,O,Z interpolation. Now interpolate T6 and log R, usually
c        on a 4x4 grid. (log(T6(i)),i=i1,i1+3),log(R(j)),j=j1,j1+3)).  Grid may
c        differ between X=0, X=.03, and X>.03 mixes, under some conditions.
c        Procedure mixes overlapping quadratics to obtain smoothed derivatives.
c
         call t6rinterp(slr,slt)
c				  ! ---------------- end of loop over X-mixes m
      enddo
c
c  Completed C,O,Z,T6,R interpolation; interpolate logKappa & derivatives in X
c
c			! for low T with 0.0 < X < 0.1, may need to reduce xdel
      xdelat = xdel
      if ( mf .eq. mxzero .and. mg .eq. mx03 .and. mh .eq. mf+2 ) then
         delhi = opk(mh,1)-opk(mg,1)
         dello = opk(mg,1)-opk(mf,1)
         if ( delhi .gt. 0.02d0 .and. delhi .lt. dello ) then
            xdelat = max( xdel*(delhi/dello)**2 , xdelmin )
            if ( delhi .lt. 0.1d0 )
     $           xdelat = xdelat + (xdel-xdelat)*((0.1d0-delhi)*
     $           12.5d0)**2
            if ( xdelat .lt. xdel ) then
               is = 0
c			 ! get (mf,mg,mh) interpolated values with revised xdel
               do i = 1,4
                  opvals(i) = qzinter(is,1,xh,2,opk(mf,i),opk(mg,i),
     $                 opk(mh,i),0.0d0,xa(mf),xa(mg),xa(mh),0.0d0,
     $                 xdelat)
                  is = 1
               enddo
            endif
         endif
      endif
c
      is = 0
c			           ! if use only one X-table
      if ( mf .eq. mh ) then
         do i = 1,4
            opvals(i) = opk(mf,i)
         enddo
c-test-xdel[
c-test-xdel;         do i = 1,n_xdel_test
c-test-xdel;            opdxi(i) = opact
c-test-xdel;         enddo
c-test-xdel]
c			          ! 2 tables: interpolate linearly in X
      else if ( mg .eq. mh ) then
         dixr = (xx(mg)-xxx)*dfsx(mg)
         do i = 1,4
            opvals(i) = opk(mf,i)*dixr + opk(mg,i)*(1.d0-dixr)
         enddo
c-test-xdel[
c-test-xdel;         do i = 1,n_xdel_test_m1
c-test-xdel;            opdxi(i) = ( opk(mf,1)
c-test-xdel;     $           * log10((xa(mg)+xdel_test(i))/(xh+xdel_test(i)))
c-test-xdel;     $           + opk(mg,1)
c-test-xdel;     $           * log10((xh+xdel_test(i))/(xa(mf)+xdel_test(i))) )
c-test-xdel;     $           / log10((xa(mg)+xdel_test(i))
c-test-xdel;     $           /(xa(mf)+xdel_test(i)))
c-test-xdel;         enddo
c-test-xdel;         opdxuse = opact
c-test-xdel]
c			           ! 3 tables: interpolate in X using quadratic
      else if ( mh .eq. mf2 ) then
c				      ! if revised xdel was NOT used (usually!)
         if ( xdelat .ge. xdel ) then
            do i = 1,4
               opvals(i) = quad(is,1,xxx,opk(mf,i),opk(mg,i),opk(mh,i),
     $              xx(mf),xx(mg),xx(mh))
               is = 1
            enddo
         endif
c-test-xdel[
c-test-xdel;         do i = 1,n_xdel_test_m1
c-test-xdel;            opdxi(i) = quad(0,1,log10(xh+xdel_test(i)),
c-test-xdel;     $           opk(mf,1),opk(mg,1),opk(mh,1),
c-test-xdel;     $           log10(xa(mf)+xdel_test(i)),
c-test-xdel;     $           log10(xa(mg)+xdel_test(i)),
c-test-xdel;     $           log10(xa(mh)+xdel_test(i)))
c-test-xdel;         enddo
c-test-xdel;         opdxuse = opact
c-test-xdel]
c		   ! 4 tables: interpolate X between two overlapping quadratics
      else
         dixr = (xx(mh)-xxx)*dfsx(mh)
c				      ! if revised xdel was NOT used (usually!)
         if ( xdelat .ge. xdel ) then
            do i = 1,4
               opvals(i) = quad(is,1,xxx,opk(mf,i),opk(mg,i),opk(mh,i),
     $              xx(mf),xx(mg),xx(mh))*dixr
     $              + quad(is,2,xxx,opk(mg,i),opk(mh,i),opk(mf2,i),
     $              xx(mg),xx(mh),xx(mf2))*(1.d0-dixr)
               is = 1
            enddo
c		  ! else, if revised xdel was used, combine it with (mg,mh,mf2)
         else
            do i = 1,4
               opvals(i) = opvals(i)*dixr
     $              + quad(is,2,xxx,opk(mg,i),opk(mh,i),opk(mf2,i),
     $              xx(mg),xx(mh),xx(mf2))*(1.d0-dixr)
               is = 1
            enddo
         endif
c-test-xdel[
c-test-xdel;         do i = 1,n_xdel_test_m1
c-test-xdel;            opdxi(i) = ( quad(0,1,log10(xh+xdel_test(i)),opk(mf,1),
c-test-xdel;     $           opk(mg,1),opk(mh,1),log10(xa(mf)+xdel_test(i)),
c-test-xdel;     $           log10(xa(mg)+xdel_test(i)),
c-test-xdel;     $           log10(xa(mh)+xdel_test(i)))
c-test-xdel;     $           * log10((xa(mh)+xdel_test(i))/(xh+xdel_test(i)))
c-test-xdel;     $           + quad(0,2,log10(xh+xdel_test(i)),opk(mg,1),
c-test-xdel;     $           opk(mh,1),opk(mf2,1),log10(xa(mg)+xdel_test(i)),
c-test-xdel;     $           log10(xa(mh)+xdel_test(i)),
c-test-xdel;     $           log10(xa(mf2)+xdel_test(i)))
c-test-xdel;     $           * log10((xh+xdel_test(i))/(xa(mg)+xdel_test(i))) )
c-test-xdel;     $           / log10((xa(mh)+xdel_test(i))
c-test-xdel;     $           /(xa(mg)+xdel_test(i)))
c-test-xdel;         enddo
c-test-xdel;         opdxuse = opact
c-test-xdel]
      endif
c
c  If the 'GN93hz' X-values are not available, just check for X > 0.76:
c
      if ( f_xhi .le. 0.0d0 ) then
c
         if ( xh .gt. 0.76d0 .and. kdo_xhi .le. 0 ) then
            fedge = max( 0.0d0 , 1.d0 - ( xh - 0.76d0 ) / 0.04d0 )
            if ( fedge .le. 0.0d0 .and. level_err .ge. 2 ) then
               write(6,30) xh
 30            format(' '/' X=',f10.6,
     $              ' > 0.8, but GN93hz X-values unavailable')
               stop ' STOP -- OPAL: Error: X too large. '
            endif
         endif
c			
c  If X > 0.03 and C+O is not too large, then the more-numerous X-values from
c  'GN93hz' can help; note that f_xhi was set above, according to this.
c
      else
c					! 0.03 < X < 0.1: set new mf,mg,mh,mf2
c					! so that only upper quadratic shifted:
         if ( mf .eq. mxzero ) then
c					! set new mf,mg,mh,mf2 so that only
c				! upper X-interpolation quadratic is shifted
            mf = 1
            mg = 1
            mh = 2
            mf2 = 3
c								! X = 1-Z-C-O:
         else if ( xh .gt. 0.999999d0 - z - xxci - xxoi ) then
c								! high-X edge
            mf = mx_hi
            mg = mx_hi
            mh = mx_hi
            mf2 = mx_hi
c					    ! X > 0.9: possibly fancy Z-interp:
         else if ( xh .ge. 0.9000009d0 ) then
c						! first: set new mf,mg,mh,mf2
            mf2 = nx_hi(kzf)
            x_4 = 1.d0 - z
            mg = min( mx_hi - 2 , mf2 - 1 )
            mf = mg - 1
            mh = mg + 1
c							! X = 0.95 case
            if ( abs( xhi_in(mh) - xh ) .lt. 9.d-7 .and.
     $           xhi_use(mf2,kzf2) .ge. xhi_in(mh) - 1.d-6 ) then
               mh = mg + 1
               mf2 = mh
               mg = mh
               mf = mh
               x_3 = 0.0d0
c						! 3-X-pt: Z > Zsto(kzf) > 0.05
            else if ( mf2 .eq. mg + 1  ) then
               mh = mf2
               x_3 = x_4
c						  ! 4-X-pt: Zsto(kzf2) < 0.05
            else if ( nx_hi(kzf2) .eq. mf2 ) then
               x_3 = xhi_use(mh,kzf2)
c							! 3-X-pt: X > 0.95
               if ( xh .ge. x_3 - 9.d-7 ) mf = mg
c							      ! Otherwise
            else
               if ( z .lt. zsto(kzh) - zacc(kzh) ) then
                  k_k = kzg
               else
                  k_k = kzh
               endif
               if ( nx_hi(k_k) .lt. mf2 ) then
c						! kzf2 *     * *     Z |
c						!      |     |  \      |
c						! kzh  *     *   *     +---> X
c						!      |     | x  \
c						! kzg  *     *     *    (use 3
c						!      |     |     |\    X-pts)
c						! kzf  *     *     * *
c						!      mf   mg    mh  mf2
                  mh = mf2
                  x_3 = x_4
               else
c						! kzf2 *     * *     Z |
c						!      |     |  \      |
c						! kzh  *     *   *     +---> X
c						!      |     |    \
c						!      |     | x  |\  (usually
c						! kzg  *     *    * *   use 4
c						!      |     |    |  \   X-pts)
c						! kzf  *     *    *   *
c						!      mf   mg   mh  mf2
c
                  if ( xhi_use(mh,kzh) .ge. xhi_in(mh) - 1.d-6 ) then
                     x_3 = xhi_in(mh)
                  else
c						  ! mh: curved Z-interpolation
                     x_3 = qzinter(0,1,z,kzf2-kzf,
     $                    xhi_use(mh,kzf),xhi_use(mh,kzg),
     $                    xhi_use(mh,kzh),xhi_use(mh,kzf2),zsto(kzf),
     $                    zsto(kzg),zsto(kzh),zsto(kzf2),zdel)
                  endif
c							! 3-X-pt: X > 0.95
                  if ( xh .ge. x_3 - 9.d-7 ) mf = mg
               endif
            endif
c					! ELSE: general case: 0.1 < X < 0.9:
         else
            mg = 2
            mh = mx_hi - 2
            do while ( mh - mg .gt. 1 )
               imd = ( mg + mh ) / 2
               if ( xh .le. xhi_in(imd) ) then
                  mh = imd
               else
                  mg = imd
               endif
            enddo
c							       ! exact X-value:
            if ( abs( xh - xhi_in(mh) ) .lt. 9.d-7 ) then
               mf = mh
               mg = mh
               mf2 = mh
            else if ( abs( xh - xhi_in(mg) ) .lt. 9.d-7 ) then
               mf = mg
               mh = mg
               mf2 = mg
c					     ! or general 4-pt X-interpolation:
            else
               mf = mg - 1
               mf2 = min( mh + 1 , nx_hi(kzf) )
            endif
c				      ! get the X-values for 3rd and 4th X-pts:
            x_3 = xhi_in(mh)
c					! 3-X-pt case
            if ( mf2 .eq. mh ) then
               x_4 = x_3
c						! x4 = 0.9 or smaller
            else if ( mf2 .le. mx_hi - 2 ) then
               x_4 = xhi_in(mf2)
c							! x4 = 0.95 or 1-Z
            else if ( nx_hi(kzf2) .eq. nx_hi(kzf) ) then
               x_4 = max( x_3 , min( xhi_use(mf2,kzf) , 1.d0 - z ) )
c								   ! otherwise
            else
               if ( z .lt. zsto(kzh) - zacc(kzh) ) then
                  k_k = kzg
               else
                  k_k = kzh
               endif
               if ( nx_hi(k_k) .lt. mf2 ) then
c						! kzf2 *     * *     Z |
c						!      |     |  \      |
c						! kzh  *     *   *     +---> X
c						!      |  x  |    \
c						! kzg  *     *     *   x4 = 1-Z
c						!      |     |     |\
c						! kzf  *     *     * *
c						!      mg   mh        mf2
                  mf2 = mx_hi
                  x_4 = 1.d0 - z
               else
c						! kzf2 *     * *     Z |
c						!      |     |  \      |
c						! kzh  *     *   *     +---> X
c						!      |     |    \
c						!      |  x  |    |\   x4 from
c						! kzg  *     *    * *  the line
c						!      |     |    |  \   "mf2"
c						! kzf  *     *    *   *
c						!      mg   mh   mf2
c
                  if ( xhi_use(mf2,kzh) .ge. xhi_in(mf2) - 1.d-6 ) then
                     x_4 = xhi_in(mf2)
                  else
c						  ! mf2: curved Z-interpolation
                     x_4 = qzinter(0,1,z,kzf2-kzf,
     $                    xhi_use(mf2,kzf),xhi_use(mf2,kzg),
     $                    xhi_use(mf2,kzh),xhi_use(mf2,kzf2),zsto(kzf),
     $                    zsto(kzg),zsto(kzh),zsto(kzf2),zdel)
                  endif
               endif
            endif
         endif
c		! m is the temporary-opacity-storage X-index
         m = 0
c			  ! loop ix over 'GN93hz'-opacity-shift X-indices
         do ix = mf, mf2
c								    ! ix valid?
            if ( ix .le. mg .or. ix .eq. mh .or. ix .eq. mf2 ) then
c
               m = m + 1
c					     ! if X available from 'Gz???.x??':
               if ( ireq_hi(ix) .eq. 0 ) then
c						! no opacity shifts at this X
                  do i = 1, 4
                     opk(m,i) = 0.0d0
                  enddo
c				   ! ELSE: if X not available from 'Gz???.x??':
               else
c					! get opacity shifts: loop over Z, T, R
                  do kz = kzf, kzf2
                     ixm = min( ix , nx_hi(kz) )
                     if ( ixm .le. 5 ) then
                        io = mo
                     else
                        ixm = ixm - 5
                        io = mo_m1
                     endif
                     do it = k1,k1+ip
                        do ir = l1,l1+iq(it-k1+1)
                           opl(it,ir,kz) = co(ixm,mc,io,it,ir,kz)
                        enddo
                     enddo
                  enddo
c						! interpolate over Z
                  call qzlog4int( zlogd )
c						! interpolate over T and R
                  call t6rinterp(slr,slt)
c
               endif
c
            endif
c
         enddo
c
c  Now add the just-computed 'GN93hz' added-X-value opacity shifts and their
c  derivatives to the original opacity and derivative values:
c
         is = 0
c				    ! if 0.03 < X < 0.1 (1st quadratic absent):
c
         if ( mg .eq. 1 .and. mf2 .eq. 3 ) then
            f_xhi = f_xhi * (1.d0-dixr)
            do i = 1, 4
               opvals(i) = opvals(i) + f_xhi * quad(is,1,xxx,opk(1,i),
     $              opk(2,i),opk(3,i),xxx_hi(1),xxx_hi(2),xxx_hi(3))
               is = 1
            enddo
c					! if use only one X-table
         else if ( mf .eq. mh ) then
            do i = 1, 4
               opvals(i) = opvals(i) + f_xhi * opk(1,i)
            enddo
c			           ! 3 tables: interpolate in X using quadratic
         else if ( mh .eq. mf2 ) then
            if ( x_3 .eq. xhi_in(mh) ) then
               x_3 = xxx_hi(mh)
            else
               x_3 = log10( x_3 + xdel )
            endif
            do i = 1, 4
               opvals(i) = opvals(i) + f_xhi * qchk(is,1,xxx,opk(1,i),
     $              opk(2,i),opk(3,i),xxx_hi(mf),xxx_hi(mg),x_3)
               is = 1
            enddo
c					! 3 tables at high-X end of matrix
         else if ( mf .eq. mg ) then
            if ( x_3 .eq. xhi_in(mh) ) then
               x_3 = xxx_hi(mh)
            else
               x_3 = log10( x_3 + xdel )
            endif
            if ( x_4 .eq. xhi_in(mf2) ) then
               x_4 = xxx_hi(mf2)
            else
               x_4 = log10( x_4 + xdel )
            endif
            do i = 1, 4
               opvals(i) = opvals(i) + f_xhi * qchk(is,1,xxx,opk(1,i),
     $              opk(2,i),opk(3,i),xxx_hi(mg),x_3,x_4)
               is = 1
            enddo
c					! 3 tables: x4 = x3 (should not happen)
         else if ( x_3 .ge. x_4 ) then
            if ( x_4 .eq. xhi_in(mf2) ) then
               x_4 = xxx_hi(mf2)
            else
               x_4 = log10( x_4 + xdel )
            endif
            do i = 1, 4
               opvals(i) = opvals(i) + f_xhi * qchk(is,1,xxx,opk(1,i),
     $              opk(2,i),opk(4,i),xxx_hi(mf),xxx_hi(mg),x_4)
               is = 1
            enddo
c		   ! 4 tables: interpolate X between two overlapping quadratics
         else
            if ( x_3 .eq. xhi_in(mh) ) then
               x_3 = xxx_hi(mh)
            else
               x_3 = log10( x_3 + xdel )
            endif
            if ( x_4 .eq. xhi_in(mf2) ) then
               x_4 = xxx_hi(mf2)
            else
               x_4 = log10( x_4 + xdel )
            endif
            dixr = ( x_3 - xxx ) / ( x_3 - xxx_hi(mg) )
            do i = 1,4
               opvals(i) = opvals(i) + f_xhi * ( qchk(is,1,xxx,opk(1,i),
     $              opk(2,i),opk(3,i),xxx_hi(mf),xxx_hi(mg),x_3) * dixr
     $              + qchk(is,2,xxx,opk(2,i),opk(3,i),opk(4,i),
     $              xxx_hi(mg),x_3,x_4) * (1.d0-dixr) )
               is = 1
            enddo
         endif
c
      endif
c
c-debug[
c-debug;      ichk = 0
c-debug;      do m = mf,mf2
c-debug;         if ( .not. abs(opk(m,1)) .le. oudebl ) ichk = 1
c-debug;      enddo
c-debug;      if ( ichk .gt. 0 .or. .not.
c-debug;     $     abs(opact) .le. oudebl .or. ioudeb .gt. 1 ) then
c-debug;         koudeb = koudeb+1
c-debug;         write(6,8415) mf,mf2,kzf,kzf2,,k1,ip,l1,iq(1),iq(2),
c-debug;     $        iq(3),iq(ip+1),z,xh,xxci,xxoi,slt,slr
c-debug; 8415    format(' '/' opk(X): m',2i2,' kz',2i3,' k1',i3,'+',i1,
c-debug;     $        ' l1',i3,'+',4i1,' Z',f10.7,' X',f10.7,' C',f10.7,
c-debug;     $        ' O',f10.7,' logT6',f12.7,' logR',f12.7)
c-debug;         do m = mf,mf2
c-debug;            write(6,8473) '    ',m,xa(m),(opk(m,i),i=1,4)
c-debug; 8473       format(a4,' (x',i1,') X=',f10.7,' logK=',g15.7,
c-debug;     $           ' DT=',g15.7,' DR=',g15.7,' DTro=',g15.7,a4)
c-debug;         enddo
c-debug;         write(6,8473) ' ==>',0,xh,(opvals(i),i=1,4),' <=='
c-debug;c-test-xdel[
c-debug;c-test-xdel;         if ( mf .ne. mh ) write(6,9387) xh,
c-debug;c-test-xdel;     $        (opdxi(i),i=1,n_xdel_test_m1),
c-debug;c-test-xdel;     $        (xdel_test(i),i=1,n_xdel_test_m1)
c-debug;c-test-xdel; 9387    format('          X=',f10.7,' logK=',8g15.7/
c-debug;c-test-xdel;     $        '                   for delX=',f9.4,7f15.4)
c-debug;c-test-xdel]
c-debug;      endif
c-debug]
c
      return
      end



*-------------------------------------------------------------------------------
      block data opacitydat
*
*     Initially at the beginning of the OPALINIT subroutine
*-------------------------------------------------------------------------------
      
      implicit double precision (a-h,o-z)


c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=28, nrb=1, nre=28,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )

      parameter ( nrdel=nrb-1, ntdel=ntb-1 )

      parameter ( zdel=0.001d0, xdel=0.03d0, xdelmin=0.001d0 )

      parameter ( badlogkval=1.d+35, badlogklim=20.d0 )
c

      parameter ( mx_hi=2*mx, mo_m1=mo-1, mx_hi_nz=mx_hi*nz )
c

c
c /xhi_opal_z/: --> data{ALL}
      common /xhi_opal_z/ xhi_in(mx_hi), xhi_use(mx_hi,nz),
     $     xxx_hi(mx_hi), nx_hi(nz), ireq_hi(mx_hi), khighx(nz),
     $     kavail_xhi, kuse_xhi, kdo_xhi
      save /xhi_opal_z/
c


c /a_opal_z/: --> data{indx,xcs,xos,xa,itime}
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime(mx),cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/

c
c /b_opal_z/: --> data{ALL BUT zz}
      common/b_opal_z/ nta(0:nrm_p1),zz(mx,nz),ntax0(0:nrm),
     $     ntax03(0:nrm), sltlo, slthi, dltlo_inv, dlthi_inv,
     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
      save /b_opal_z/
c
c COMMON /extopac/  : variables returning opacity values, as described above

      common /extopac/ opact,dopact,dopacr,dopactd,fedge,ftredge,fzedge
      save /extopac/
c

c
c /recoin_opal_z/: --> data{ALL}
      common/recoin_opal_z/ itimeco,mxzero,mx03,kope,igznotgx
      save /recoin_opal_z/
c

c
      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
c
c COMMON /opalmixes/ : makeup of Z in opacity mixtures (described above):
c
c /opalmixes/: --> data{ALL BUT xofe_opalmixes}
c-implicit;      real*4 xiz_mix,fninz_mix,bracketife_mix,bracketofe_opalmixes,
c-implicit;     $     xofe_opalmixes,xiz_opalmixes,fninz_opalmixes
      character*2 cel_opalmixes(nel_zmix)
      character*8 cfile_opalmixes(n_zmixes)
      common/opalmixes/ xiz_mix(nel_zmix),fninz_mix(nel_zmix),
     $     bracketife_mix(nel_zmix),bracketofe_opalmixes(n_zmixes),
     $     xofe_opalmixes(n_zmixes),xiz_opalmixes(nel_zmix,n_zmixes),
     $     fninz_opalmixes(nel_zmix,n_zmixes),
     $     cel_opalmixes,cfile_opalmixes
      save /opalmixes/
c
c PARAMETERS defining the matrices used for the Z-interpolation:
c  mz = 8 = number of metallicities Z available in the 'Gz???.x??' files
c  mzhi = 11 = number of Z-indices for the mix-number matrix nofz (see below)
c  mzal = 13 = number of metallicities Z available in the 'GN93hz' file
c  nzm = 14 = combined total number of metallicities Z available in files
c  nadd_zavail = 6 = number of metallicities Z besides those in 'Gz???.x??'
c
      parameter ( mz=8, mz_m1=mz-1, mz_m2=mz-2, mzhi=11, mzal=13,
     $     nzm=mzal+1, nadd_zavail=nzm-mz )
c
c COMMON /zinter_opal_z/ : values used in Z-interpolation:
c  zvalhi(mzhi) = Z-range limits for the mix-number matrix nofz (see below)
c  nofz(mzhi,5,mo) = the number of different C-tables at each O-tabulation and
c                     X-tabulation value, for each relevant range of Z
c  mnofz(mx) = X-table m-index in nofz corresponding to the m-index in xa(mx):
c               if the X-table abundances in xa are unchanged,  mnofz(i) = i
c  zval(mz) = Z-tabulation values available in the 'Gz???.x??' files
c  zalval(mzal) = Z-tabulation values available in the 'GN93hz' file
c  zavail(nzm) = combined Z-tabulation values available in the files
c  iadd_zavail(nadd_zavail) = best order in which to reduce the intervals in
c                              zval() by adding metallicities from zavail()
c
c /zinter_opal_z/: --> data{ALL}
      common/zinter_opal_z/ zvalhi(mzhi),nofz(mzhi,5,mo),mnofz(mx),
     $     zval(mz),zalval(mzal),zavail(nzm),iadd_zavail(nadd_zavail)
      save /zinter_opal_z/
c
c COMMON /czinte_opal_z/ : the X- and Z- parts of the 'Gx??z*' file names;
c                          also used to specify the 'Gz???.x??' file names:
c
c /czinte_opal_z/: --> data{ALL}
      character*4 cxfil(5),czfil(mz)
      common/czinte_opal_z/ cxfil,czfil
      save /czinte_opal_z/
c
c COMMON /c_opal_ctrl_smooth/ : flags to control the opacity smoothing:
c  init_smo = 0 : do not smooth the opacities on input
c           = 1 (default): on input, subroutine OPALTAB smooths the opacities
c  low_CO_smo = 0 : do not perform this CO-direction smoothing
c             = 1 (default): on input, a few opacities in the 3 mixes adjacent
c                    to the C=O=0.0 mix (i.e., the 3 mixes with C or O = 0.01,
c                    C+O no more than 0.02) are smoothed in the C-O direction,
c                    if opacity changes between mixes with C,O = 0.0, 0.03, 0.1
c                    are monotonic but the opacity at C,O = 0.01 does not fit
c                    the trend; the resulting adjustments are small, and only
c                    occur at a small minority of the (T6,R) points
c  interp_CO_smo = 0 : use the old subroutine COINTERP for interpolating among
c                       CO-rich opacities when OPAC or OPAL is called
c                = 1 (default): use the new subroutine COINTSMO instead, for
c                       smoother interpolation among CO-rich opacities
c
c /c_opal_ctrl_smooth/: --> data{ALL}
      common/c_opal_ctrl_smooth/ init_smo, low_CO_smo, interp_CO_smo
      save /c_opal_ctrl_smooth/
c
c COMMON /opdir/ : copdir = the name of the directory holding the opacity files
c  Note that in other routines copdir is defined to be CHARACTER*80 ; here it
c  is merely being initialized to be blank, meaning "use the current directory"
c
c /opdir/: --> data{copdir}
      character*80 copdir
      common/opdir/ copdir
      save /opdir/
c
      character*1 chpdir(80)
      equivalence (chpdir(1),copdir)
c
c COMMON /chkpoc/ : character(s) allowed to terminate a directory name
c
      character*1 cb(6)
      common/chkpoc/cb
      save /chkpoc/
c
c COMMON /alink_opal_z/ : contains data needed for smoothing routine OPALTAB
c
      common/alink_opal_z/ NTEMP,NSM,nrlow,nrhigh,RLE,t6arr(100),
     $     coff(100,nrm)
      save /alink_opal_z/
c
c COMMON /d_opal_z/ :
c  dkap = derivative of value returned by quadratic interpolation function QDER
c
c /d_opal_z/: --> data{dkap}
      common/d_opal_z/ dkap
      save /d_opal_z/
c
c COMMON /c_level_err_opal_z/ :
c  level_err = error-checking level, set by SET_ERR_CHECK; allowed values are:
c                0 : Do not check input Nzin, Zlo, Z, Zhi in READZEXCO.
c                1 : (Default): Do check these; fatal error if checks fail.
c                2 : In this case, it is also a fatal error if FEDGE = 0.
c
c /c_level_err_opal_z/: --> data{level_err}
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/

c /xhi_opal_z/ data:
c
      data xhi_in / 0., 0.1, 0.2, 0.35, 0.5, 0.7, 0.8, 0.9, 0.95, 1. /
      data xhi_use / mx_hi_nz * -1.0 /, xxx_hi / mx_hi * -9.0 /
      data nx_hi / nz * 0 /
      data ireq_hi / 0, 0, 1, 0, 1, 0, 1, 1, 1, 1 /
      data khighx / nz * 0 /
      data kavail_xhi / 0 /, kuse_xhi / 2 /, kdo_xhi / 0 /
c
c /a_opal_z/ data:
c								! indx(1:101)
      data indx/1,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4,
     $     4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,
     $     6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
     $     7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7/
c								! X(1:mx) mx=5
      data (xa(i),i=1,5)/ 0.0, 0.03, 0.1, 0.35, 0.7 /
c								! C,O(1:mc)
      data xcs/0.0,0.01,0.03,0.1,0.2,0.4,0.6,1.0/
      data xos/0.0,0.01,0.03,0.1,0.2,0.4,0.6,1.0/
c								! init-flags
      data itime/mx*0/
c
c /b_opal_z/ data:
c								! nta(0:nrm_p1)
c!!      data nta/57, 70,70,70,70,70, 70,70,70,70,70,
c!!     $     70,70,70,70,69, 64,60,58,57, -99/
cc								! ntax0(0:nrm)
c!!      data ntax0/999, 6,5,5,5,4, 4,4,3,1,1, 1,1,1,1,1, 1,1,1,1/
cc								! ntax03(0:nrm)
c!!      data ntax03/999, 5,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1, 1,1,1,1/

c                                                               ! nta(0:nrm_p1)
      data nta/57, 76,76,76,76,76, 76,76,76,76,76, 76,76,76,76,
     $     76,76,76,76,76, 76,76,76,76,76, 76,76,76,76, -99/
c								! ntax0(0:nrm)
      data ntax0/999, 1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,
     $     1,1,1,1,1, 1,1,1,1/
c								! ntax03(0:nrm)
      data ntax03/999, 5,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,
     $     1,1,1,1,1, 1,1,1,1/
c
      data sltlo/-99.d0/, slthi/-99.d0/, dltlo_inv/-99.d0/, 
     $     dlthi_inv/-99.d0/

      data slrlo/-99.d0/, slrhi/-99.d0/, dlrlo_inv/-99.d0/, 
     $     dlrhi_inv/-99.d0/

      data init_trvals/0/
c
c /extopac/ data:
c
      data opact/1.d35/, dopact/0.d0/, dopacr/0.d0/, dopactd/0.d0/,
     $     fedge/0.d0/, ftredge/0.d0/, fzedge/0.d0/
c
c /recoin_opal_z/ data:
c
      data itimeco/0/, mxzero/1/, mx03/2/, kope/0/, igznotgx/0/
c
c /opalmixes/ data:
c
      data cfile_opalmixes/'GN93hz  ','Alrd96a2','C95hz   ','W95hz   ',
     $     '        '/
      data cel_opalmixes/'C ','N ','O ','Ne','Na','Mg','Al','Si',
     $     'P ','S ','Cl','Ar','K ','Ca','Ti','Cr','Mn','Fe','Ni'/
      data xiz_opalmixes/
     $     0.173285d0,0.053152d0,0.482273d0,0.098668d0,0.001999d0,
     $     0.037573d0,0.003238d0,0.040520d0,0.000355d0,0.021142d0,
     $     0.000456d0,0.005379d0,0.000210d0,0.003734d0,0.000211d0,
     $     0.001005d0,0.000548d0,0.071794d0,0.004459d0,
     $     0.102693d0,0.031499d0,0.570253d0,0.116656d0,0.002363d0,
     $     0.044428d0,0.000962d0,0.047912d0,0.000420d0,0.024999d0,
     $     0.000539d0,0.006360d0,0.000124d0,0.004415d0,0.000245d0,
     $     0.000595d0,0.000230d0,0.042538d0,0.002769d0,
     $     0.091924d0,0.028196d0,0.642620d0,0.052341d0,0.001060d0,
     $     0.050066d0,0.001718d0,0.053992d0,0.000188d0,0.028172d0,
     $     0.000242d0,0.002853d0,0.000279d0,0.004975d0,0.000275d0,
     $     0.000533d0,0.000116d0,0.038085d0,0.002365d0,
     $     0.076451d0,0.023450d0,0.672836d0,0.084869d0,0.000882d0,
     $     0.041639d0,0.001428d0,0.035669d0,0.000157d0,0.019942d0,
     $     0.000201d0,0.002373d0,0.000092d0,0.005209d0,0.000387d0,
     $     0.000443d0,0.000242d0,0.031675d0,0.002056d0,
     $     0.173285d0,0.053152d0,0.482273d0,0.098668d0,0.001999d0,
     $     0.037573d0,0.003238d0,0.040520d0,0.000355d0,0.021142d0,
     $     0.000456d0,0.005379d0,0.000210d0,0.003734d0,0.000211d0,
     $     0.001005d0,0.000548d0,0.071794d0,0.004459d0/
      data xiz_mix/
     $     0.173285d0,0.053152d0,0.482273d0,0.098668d0,0.001999d0,
     $     0.037573d0,0.003238d0,0.040520d0,0.000355d0,0.021142d0,
     $     0.000456d0,0.005379d0,0.000210d0,0.003734d0,0.000211d0,
     $     0.001005d0,0.000548d0,0.071794d0,0.004459d0/
      data fninz_opalmixes/
     $     0.245518d0,0.064578d0,0.512966d0,0.083210d0,0.001479d0,
     $     0.026308d0,0.002042d0,0.024552d0,0.000195d0,0.011222d0,
     $     0.000219d0,0.002291d0,0.000091d0,0.001586d0,0.000075d0,
     $     0.000329d0,0.000170d0,0.021877d0,0.001293d0,
     $     0.147909d0,0.038904d0,0.616594d0,0.100010d0,0.001778d0,
     $     0.031622d0,0.000617d0,0.029512d0,0.000234d0,0.013490d0,
     $     0.000263d0,0.002754d0,0.000055d0,0.001906d0,0.000089d0,
     $     0.000198d0,0.000072d0,0.013177d0,0.000816d0,
     $     0.131157d0,0.034498d0,0.688325d0,0.044451d0,0.000790d0,
     $     0.035301d0,0.001091d0,0.032945d0,0.000104d0,0.015059d0,
     $     0.000117d0,0.001224d0,0.000122d0,0.002127d0,0.000099d0,
     $     0.000176d0,0.000036d0,0.011687d0,0.000691d0,
     $     0.108211d0,0.028462d0,0.714945d0,0.071502d0,0.000652d0,
     $     0.029125d0,0.000900d0,0.021591d0,0.000086d0,0.010575d0,
     $     0.000096d0,0.001010d0,0.000040d0,0.002210d0,0.000137d0,
     $     0.000145d0,0.000075d0,0.009642d0,0.000595d0,
     $     0.245518d0,0.064578d0,0.512966d0,0.083210d0,0.001479d0,
     $     0.026308d0,0.002042d0,0.024552d0,0.000195d0,0.011222d0,
     $     0.000219d0,0.002291d0,0.000091d0,0.001586d0,0.000075d0,
     $     0.000329d0,0.000170d0,0.021877d0,0.001293d0/
      data fninz_mix/
     $     0.245518d0,0.064578d0,0.512966d0,0.083210d0,0.001479d0,
     $     0.026308d0,0.002042d0,0.024552d0,0.000195d0,0.011222d0,
     $     0.000219d0,0.002291d0,0.000091d0,0.001586d0,0.000075d0,
     $     0.000329d0,0.000170d0,0.021877d0,0.001293d0/
      data bracketife_mix/nel_zmix*0.d0/
      data bracketofe_opalmixes/0.0d0,0.3d0,0.4d0,0.5d0,0.d0/
c
c /zinter_opal_z/ data: (note: zavail, iadd_zavail are computed in get_zavail)
c
      data zavail/nzm*0.d0/, iadd_zavail/nadd_zavail*0/
c
      data zvalhi/0.,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1/
c
      data zval/0.,0.001,0.004,0.01,0.02,0.03,0.05,0.1/
c
      data zalval/0.,0.0001,0.0003,0.001,0.002,0.004,0.01,0.02,
     $     0.03,0.04,0.06,0.08,0.1/
c
      data (mnofz(i),i=1,mx)/1,2,3,4,5/
c
      data nofz/
     $     8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,8,8,8,8,8,8,
     $     8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,7,7,7,7,7,7,
     $     6,6,6,6,6,6,6,6,6,6,5,
     $     8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,8,8,8,8,8,8,
     $     8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,7,7,7,7,7,7,7,
     $     6,6,6,6,6,6,6,6,6,5,5,
     $     8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,8,8,8,8,8,8,
     $     8,8,8,8,8,8,8,8,8,8,8, 8,8,7,7,7,7,7,7,7,7,7,
     $     6,6,6,6,6,6,6,5,5,5,5,
     $     8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,8,8,8,8,8,8,
     $     8,8,8,8,8,8,8,8,8,8,8, 7,7,7,7,7,7,7,7,7,7,7,
     $     5,5,5,5,5,5,5,5,5,5,4,
     $     8,8,8,8,8,8,8,8,8,8,8, 8,8,8,8,8,8,8,8,8,8,8,
     $     8,8,8,8,8,8,8,8,8,8,7, 7,7,7,7,7,6,6,6,6,6,6,
     $     4,4,4,4,4,4,4,3,3,2,1,
     $     7,7,7,7,7,7,7,7,7,7,7, 7,7,7,7,7,7,7,7,7,7,7,
     $     7,7,7,7,7,7,7,7,7,7,6, 6,6,6,6,6,5,5,5,5,5,5,
     $     0,0,0,0,0,0,0,0,0,0,0,
     $     6,6,6,6,6,6,6,6,6,6,6, 6,6,6,6,6,6,6,6,6,6,6,
     $     6,6,6,6,6,6,6,6,6,6,5, 4,4,3,3,2,1,0,0,0,0,0,
     $     0,0,0,0,0,0,0,0,0,0,0,
     $     0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,
     $     0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,
     $     0,0,0,0,0,0,0,0,0,0,0/
c
c /czinte_opal_z/ data:
c
      data cxfil/'Gx00','Gx03','Gx10','Gx35','Gx70'/
      data czfil/'z0  ','z001','z004','z01 ',
     $     'z02 ','z03 ','z05 ','z10 '/
c
c /c_opal_ctrl_smooth/ data:
c
      data init_smo/1/, low_CO_smo/1/, interp_CO_smo/1/
c
c /opdir/ data:
c
      data chpdir/80*' '/
c
c /d_opal_z/ data:
c
      data dkap/0.0/
c
c /c_level_err_opal_z/ data:
c
      data level_err / 1 /

      data cb/ '/', '/', '_', '~', '+', '-' /


      end

*-------------------------------------------------------------------------------
      subroutine opalinit(khighz,iulow,ofebrack,z,kz)
*     ===============================================
*
*  INITIALIZATIONS AND OPACITY FILE SETUP:
*
*  This subroutine performs some initializations that would otherwise be done
*  at the beginning of subroutine READZEXCO.  These do some grid set-up, look
*  for the user-supplied non-zero [O/Fe] file if khighz = 5, calculate the
*  [O/Fe] values for each of the possible mixes, and find the OPAL-opacity-file
*  directory name.
*
*
*-------------------------------------------------------------------------------
      
      implicit double precision (a-h,o-z)

c PARAMETERS defining the matrices used to hold the opacity values:
c  nz = 14 = maximum number of Z-tabulation values (see arrays zavail, zsto)
c  mx = 5 = number of X-tabulation values (see array xa)
c  mc = mo = 8 = number of C- and O-tabulation values (see arrays xcs and xos)
c  nrm = 19 = maximum number of R-tabulation values (see array alrf)
c  nrb = 1 = index of lowest R-value to be stored         \  default: store
c  nre = 19 = index of highest R-value to be stored        }  all 19 R-values
c  nr = nre-nrb+1 = 19 = number of R-values to be stored  /
c  ntm = 70 = maximum number of T-tabulation values (see array flogtin)
c  ntb = 1 = index of lowest T-value to be stored         \  default: store
c  nt = ntm-ntb+1 = 70 = number of T-values to be stored  /   all 70 T-values
c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=28, nrb=1, nre=28,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
c PARAMETERS nrdel and ntdel give matrix position differences that result from
c  any reduction of nr or nt below nrm or ntm, respectively
c
      parameter ( nrdel=nrb-1, ntdel=ntb-1 )
c
c PARAMETERS:
c  zdel = 0.001 = offset for Z, Z+C, and Z+O, to make log interpolation behave
c                  reasonably at small Z values (zdel is only used in READEXCO)
c  xdel = 0.03 = usual (high-T) offset for X, to make log interpolation behave
c                 reasonably at small X; note that 0.03 works better than 0.005
c  xdelmin = 0.001 = lowest value of X offset ever used (at low temperatures)
c
      parameter ( zdel=0.001d0, xdel=0.03d0, xdelmin=0.001d0 )
c
c PARAMETERS: value used for "bad" (missing) logKappa values:
c
      parameter ( badlogkval=1.d+35, badlogklim=20.d0 )
c
c PARAMETERS defining the storage for the additional X-values from 'GN93hz':
c  mx_hi = 10 = number of X-values contained in 'GN93hz'.
c  mo_m1 = mo - 1 = 7 : used for the position in the matrix  co()  below where
c                         some of the 'GN93hz' opacity-shifts will be stored
c                         (see COMMON /a_opal_z/ below, for this matrix)
c  mx_hi_nz = mx_hi * nz : the number of initialization values requred for the
c                            matrix  xhi_use()  in COMMON /xhi_opal_z/ below.
c
      parameter ( mx_hi=2*mx, mo_m1=mo-1, mx_hi_nz=mx_hi*nz )
c
c COMMON /xhi_opal_z/ : auxiliary matrices for additional 'GN93hz' X-values:
c  xhi_in(mx_hi) = { 0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0 },
c                    the X-values contained in 'GN93hz' (for the case Z=0.0)
c  xhi_use(mx_hi,nz) = the 'GN93hz' X-values available for each stored Z-value
c                        (indexed by kz = 1, ..., numz),  with the exception
c                        that xhi_use(1,kz) = 0.03 (from 'Gz???.x03' files);
c                        note that in each case the highest value is 1 - Z
c  xxx_hi(mx_hi) = log10( xhi_in() + xdel ) = logarithmic 'GN93hz' X-values;
c                    note exception that xxx_hi(1) = log10( 0.03 + xdel )
c  nx_hi(nz) = number of 'GN93hz' X-values in xhi_use() at each stored Z-value
c  ireq_hi(mx_hi) = flags to tell whether the corresponding 'GN93hz' X-values
c                     are unavailable from the 'Gz???.x??' files
c  khighx(nz) = flag to tell whether the 'GN93hz' opacities were read in, for
c                 each of the stored Z-values
c  kavail_xhi = flag to tell whether the 'GN93hz' opacities are available
c  kuse_xhi = flag to tell whether the 'GN93hz' X-values shoud be used for
c               X-interpolation (see description of subroutine SET_KXHI above)
c  kdo_xhi = internal flag controlling use of the 'GN93hz' X-values
c
c /xhi_opal_z/: --> data{ALL}
      common /xhi_opal_z/ xhi_in(mx_hi), xhi_use(mx_hi,nz),
     $     xxx_hi(mx_hi), nx_hi(nz), ireq_hi(mx_hi), khighx(nz),
     $     kavail_xhi, kuse_xhi, kdo_xhi
      save /xhi_opal_z/
c
c COMMON /a_opal_z/ : matrices for opacity storage and interpolation:
c  indx(101) is used to get tabulation-index from excess C or O values: INDX
c              (see data) refers to the index i of abundance grid points xc(i)
c              or xo(i): e.g., i = indx( int( 100. * max( C , 0. ) + 1. ) )
c  t6list(nt) = stored-grid T6 values (obtained from logT values)
c  alr(nr) = stored-grid log(R) values (at log(R)-intervals of 0.5)
c  n(mx,mo,nz) = the number of different C-tables present at each X-tabulation,
c                 O-tabulation, and Z-tabulation value (see also array  nofz
c                 in common/zinter_opal_z/ )
c  alt(nt) = stored-grid log(T6) values; note log(T6) = log(T) - 6.
c  dfs(nt) = stored-grid inverse log(T6)-spacing: dfs(i) = 1./[alt(i)-alt(i-1)]
c             (note that dfs(1) = dfs(2), for convenience)
c  dfsr(nr) = stored-grid inverse log(R)-spacing (unlike dfs, the dfsr values
c              are all equal): dfsr(i) = 1./[alr(i)-alr(i-1)] = 2.
c  b(3) is used to hold intermediate values during C-O interpolation
c  m = the index of the current X-table
c  mf = the lowest (first) index of the X-table(s) needed for X-interpolation
c  xa(8) = the tabulated X-values: actually there are only five: see "data"
c  alrf(nrm) = opacity-file-grid log(R) values (note that this grid, with nrm
c               log(R) values, may be larger than the stored-grid, with nr)
c  flogtin(ntm) = opacity-file-grid log(T) values (again, ntm may be > nt)
c  dfsx(mx) = inverse logarithmic-X-grid spacings: dfsx(i) = 1./[xx(i)-xx(i-1)]
c  oxf(mx,mc,nz) = logarithmic-O-grid tabulation values for a given C value:
c                   for each X-table index m and Z-table index k,  oxf(m,i,k) =
c                   log10[ min{ xos(i) , 1-xa(m)-zsto(k) } + zsto(k) + 0.001 ]
c  cxf(mx,mo,nz) = logarithmic-C-grid tabulation values similarly
c  xcdf(mx,mo,nz) = maximum possible C value for a given O value:
c                    for each X-table index m and Z-table index kz,
c                    xcdf(m,i,kz) = 1 - xa(m) - zsto(kz) - xo(i)
c  xodf(mx,mc,nz) = maximum possible O value for a given C value, similarly
c  itime(mx) = "opacities available" flags for each X-table index m: each
c                initialization-flag itime(m) is initially zero, and is set to
c                12345678 after opacities have been read in and initializations
c                performed for that X-index: see "data"
c  cxdf(mx,mo,nz) = logarithmic maximum C value for a given O value:
c                    for each X-table index m and Z-table index kz,
c                    cxdf(m,i,kz) = log10[ xcdf(m,i,kz) + zsto(kz) + 0.001 ]
c  oxdf(mx,mc,nz) = logarithmic maximum O value for a given C value, similarly
c  q(4) = temporary: opacity-derivative at each T, in T-R interpolation
c  h(4) = temporary: opacities log10(Kappa) at each T, in T-R interpolation
c  xcd(mc) = maximum possible C at present m and kz:  xcd(i) = xcdf(m,i,kz)
c  xod(mc) = maximum possible O at present m and kz:  xod(i) = xodf(m,i,kz)
c  xc(mc) = C-tabulation at present m and kz: xc(i) = min{ xcs(i) , 1-xa(m)-Z }
c  xo(mo) = O-tabulation at present m and kz: xo(i) = min{ xos(i) , 1-xa(m)-Z }
c  xcs(mc) = C-tabulation values: see "data"
c  xos(mo) = O-tabulation values: see "data"
c  cxd(mc) = logarithmic maximum C value for a given O value at present X-table
c             index m and Z-table index kz: cxd(i) = cxdf(m,i,kz)
c  oxd(mo) = logarithmic maximum O value for a given C value, similarly
c  cx(mc) = logarithmic-C-grid tabulation values at present X-table index m and
c            Z-table index kz: cx(i) = cxf(m,i,kz)
c  ox(mo) = logarithmic-O-grid tabulation values at present m, similarly
c  zzz(nz) = shifted Z-value (for logarithmic interpolation purposes):
c             for each Z-table index kz,  zzz(kz) = zsto(kz) + 0.001
c  xxh = xh = stored value of desired X value (xxh is never actually used)
c  xx(mx) = logarithmic-X-grid tabulation values: xx(i) = log10[ xa(i) + 0.03 ]
c            (note that previous program versions added 0.005 rather than 0.03;
c            the latter works better for X near zero at log(T) > 5.)
c  nc = n(m,1,kz) = number of C-grid values available at X,Z-table indices m,kz
c  no = nc = number of O-grid values, similarly
c  zsto(nz) = stored Z-values available for Z-interpolation
c  zvint(nz) = logartihmic stored Z-values: zvint(i) = log10[ zsto(i) + 0.001 ]
c  dfsz(nz) = inverse log-Z-grid spacings: dfsz(i) = 1./[zvint(i)-zvint(i-1)]
c  zacc(nz) = accuracy with which Z must match the corresponding stored Z-value
c              in order to be considered equal to it
c  zlow,zmiddle,zhigh = lowest, "typical", and highest stored Z-values
c  zlo_ex,zhi_ex = extreme Z-limits for Z-extrapolation
c  numz = number of stored Z-values available for Z-interpolation
c  co(mx,mc,mo,nt,nr,nz) = stored opacities log10(Kappa): co(m,i,j,k,l,kz) =
c                           log10[ Kappa(X_m,C_i,O_j,T_k,R_l,Z_kz) ] , where
c                           X_m = xa(m), C_i = min{xcs(i),1-xa(m)-Z-xos(j)},
c                           O_j = xos(j), T_k = alt(k), R_l = alr(l), and
c                           Z_kz = zsto(kz); except that, for j = mo, the
c                           "diagonal" tables are stored, with C_i = xcs(i) and
c                           O_j = 1-xa(m)-zsto(kz)-xcs(i).  Note that not quite
c                           all (m,i,j) locations are used; unused locations
c                           (m,mc,mo,...) and (m,mc,mo-1,...) are used for
c                           temporary storage for opacities from 'GN93hz' and
c                           the file with non-zero [O/Fe], if these are used.
c******************************************************************************
c  Note that the old arrays "diag" & "diago" are not needed; here:
c  co(m,n(m,i,kz),i,it,ir,kz) = diag(m,i,it,ir,kz)  and
c  co(m,i,mo,it,ir,kz) = diago(m,n(m,1,kz)-i,it,ir,kz) = diago(m,no-i,it,ir,kz)
c******************************************************************************
c  opk(mx,4) = temporary: for each X-table index m used in the X-interpolation,
c               opk(m,1:4) holds the log10(Kappa) value and the derivatives
c               with respect to T (at constant R), to R (at constant T), and to
c               T (at constant density), for that m value (already interpolated
c               in C, O, Z, T, and R)
c  opl(nt,nr,nz) = temporary: for each stored-grid (T_k,R_l,Z_kz) value used in
c                   the T-R-Z interpolation,  opl(k,l,kz)  holds the opacity
c                   log10(Kappa) at that T_k, R_l, and Z_kz (already
c                   interpolated in C and O); the Z-interpolated values at each
c                   (T_k,R_l) are subsequently stored in  opl(k,l,1)
c  cof(nt,nr) = temporary: in the subroutine READEXCO, cof is used temporarily
c                to hold the opacities for non-zero [O/Fe] where they will not
c                be overwritten when reading in the 'GN93hz' opacities; in the
c                subroutine COINTSMO, cof is used temporarily to hold some
c                logarithmic C and O grid-values
c
c /a_opal_z/: --> data{indx,xcs,xos,xa,itime}
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime(mx),cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
c COMMON /b_opal_z/ : high and low logT6 and logR limits, and mix Z-values:
c  nta(0:nrm+1) = maximum T-index where opacity values exist for each R-index;
c                  if nre < nrm, then nta(nre+1:nrm) will be set to -99, to
c                  indicate that no opacities are present there (values at
c                  positions < nrb are NOT reset, however); see "data"
c  zz(mx,nz) = Z values used for each X-table index m and each stored Z-value;
c               those at differing m but the same Z-index must all be the same!
c  ntax0(0:nrm) = minimum T-index where opacity values exist in the X=0.0 (m=1)
c                  tables; if nrb > 1, then ntax0(0:nrb-1) will be set to 999,
c                  to indicate that no opacities are present there (values at
c                  positions > nre are NOT reset, however).  The worst-case
c                  values given below (see "data") are reset to the actual
c                  values when the opacities are read in by READEXCO
c  ntax03(0:nrm) = minimum T-index in the X=0.03 (m=2) tables, similarly
c  sltlo, slthi = low and high logT6 limits: a logT6 value outside this range
c                  is considered to require extrapolation (by default these
c                  these limits lie at the boundaries of the matrix; they may
c                  be reset to lie inside it, but not outside it)
c  dltlo_inv, dlthi_inv = (inverse of) allowed extrapolation delta_logT6 beyond
c                          the above limits sltlo, slthi: by default, one grid
c                          spacing beyond the edge of the matrix (they can be
c                          reset, but in no case will extrapolation be allowed
c                          more than one grid spacing beyond the matrix edge)
c  slrlo, slrhi = low and high logR limits (handled similarly to logT6 limits)
c  dlrlo_inv, dlrhi_inv = (inverse of) allowed extrapolation delta_logR
c
c /b_opal_z/: --> data{ALL BUT zz}
      common/b_opal_z/ nta(0:nrm_p1),zz(mx,nz),ntax0(0:nrm),
     $     ntax03(0:nrm), sltlo, slthi, dltlo_inv, dlthi_inv,
     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
      save /b_opal_z/
c
c COMMON /extopac/  : variables returning opacity values, as described above

      common /extopac/ opact,dopact,dopacr,dopactd,fedge,ftredge,fzedge
      save /extopac/
c
c COMMON /recoin_opal_z/ :
c  itimeco = 12345678 after initializations have been carried out (init-flag)
c  mxzero = 1 is set to the X-index of the mix with X=0
c  mx03 = 2 is set to the X-index of the mix with X=0.03
c  kope is the length of the directory-name (where the opacity files are)
c  igznotgx is a flag to tell whether the OPAL opacity file names are in the
c            new form  Gz???.x??  rather than the old form  Gx??z*
c            (initially,  igznotgx = 0  meaning "look for either")
c
c /recoin_opal_z/: --> data{ALL}
      common/recoin_opal_z/ itimeco,mxzero,mx03,kope,igznotgx
      save /recoin_opal_z/
c
c PARAMETERS defining the matrices used to hold the mix compositions:
c  nel_zmix = 19 = number of heavy elements in the opacity mixture (C thru Ni)
c  n_zmixes = 5 = number of "hz" opacity files available (the fifth one is for
c                  a user-supplied "hz" opacity file, with non-zero [O/Fe])
c  kel_o = 3 = position of oxygen (O) in the list of mix-elements
c  kel_fe = nel_zmix-1 = 18 = position of iron (Fe) in the list of mix-elements
c
      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
c
c COMMON /opalmixes/ : makeup of Z in opacity mixtures (described above):
c
c /opalmixes/: --> data{ALL BUT xofe_opalmixes}
c-implicit;      real*4 xiz_mix,fninz_mix,bracketife_mix,bracketofe_opalmixes,
c-implicit;     $     xofe_opalmixes,xiz_opalmixes,fninz_opalmixes
      character*2 cel_opalmixes(nel_zmix)
      character*8 cfile_opalmixes(n_zmixes)
      common/opalmixes/ xiz_mix(nel_zmix),fninz_mix(nel_zmix),
     $     bracketife_mix(nel_zmix),bracketofe_opalmixes(n_zmixes),
     $     xofe_opalmixes(n_zmixes),xiz_opalmixes(nel_zmix,n_zmixes),
     $     fninz_opalmixes(nel_zmix,n_zmixes),
     $     cel_opalmixes,cfile_opalmixes
      save /opalmixes/
c
c PARAMETERS defining the matrices used for the Z-interpolation:
c  mz = 8 = number of metallicities Z available in the 'Gz???.x??' files
c  mzhi = 11 = number of Z-indices for the mix-number matrix nofz (see below)
c  mzal = 13 = number of metallicities Z available in the 'GN93hz' file
c  nzm = 14 = combined total number of metallicities Z available in files
c  nadd_zavail = 6 = number of metallicities Z besides those in 'Gz???.x??'
c
      parameter ( mz=8, mz_m1=mz-1, mz_m2=mz-2, mzhi=11, mzal=13,
     $     nzm=mzal+1, nadd_zavail=nzm-mz )
c
c COMMON /zinter_opal_z/ : values used in Z-interpolation:
c  zvalhi(mzhi) = Z-range limits for the mix-number matrix nofz (see below)
c  nofz(mzhi,5,mo) = the number of different C-tables at each O-tabulation and
c                     X-tabulation value, for each relevant range of Z
c  mnofz(mx) = X-table m-index in nofz corresponding to the m-index in xa(mx):
c               if the X-table abundances in xa are unchanged,  mnofz(i) = i
c  zval(mz) = Z-tabulation values available in the 'Gz???.x??' files
c  zalval(mzal) = Z-tabulation values available in the 'GN93hz' file
c  zavail(nzm) = combined Z-tabulation values available in the files
c  iadd_zavail(nadd_zavail) = best order in which to reduce the intervals in
c                              zval() by adding metallicities from zavail()
c
c /zinter_opal_z/: --> data{ALL}
      common/zinter_opal_z/ zvalhi(mzhi),nofz(mzhi,5,mo),mnofz(mx),
     $     zval(mz),zalval(mzal),zavail(nzm),iadd_zavail(nadd_zavail)
      save /zinter_opal_z/
c
c COMMON /czinte_opal_z/ : the X- and Z- parts of the 'Gx??z*' file names;
c                          also used to specify the 'Gz???.x??' file names:
c
c /czinte_opal_z/: --> data{ALL}
      character*4 cxfil(5),czfil(mz)
      common/czinte_opal_z/ cxfil,czfil
      save /czinte_opal_z/
c
c COMMON /c_opal_ctrl_smooth/ : flags to control the opacity smoothing:
c  init_smo = 0 : do not smooth the opacities on input
c           = 1 (default): on input, subroutine OPALTAB smooths the opacities
c  low_CO_smo = 0 : do not perform this CO-direction smoothing
c             = 1 (default): on input, a few opacities in the 3 mixes adjacent
c                    to the C=O=0.0 mix (i.e., the 3 mixes with C or O = 0.01,
c                    C+O no more than 0.02) are smoothed in the C-O direction,
c                    if opacity changes between mixes with C,O = 0.0, 0.03, 0.1
c                    are monotonic but the opacity at C,O = 0.01 does not fit
c                    the trend; the resulting adjustments are small, and only
c                    occur at a small minority of the (T6,R) points
c  interp_CO_smo = 0 : use the old subroutine COINTERP for interpolating among
c                       CO-rich opacities when OPAC or OPAL is called
c                = 1 (default): use the new subroutine COINTSMO instead, for
c                       smoother interpolation among CO-rich opacities
c
c /c_opal_ctrl_smooth/: --> data{ALL}
      common/c_opal_ctrl_smooth/ init_smo, low_CO_smo, interp_CO_smo
      save /c_opal_ctrl_smooth/
c
c COMMON /opdir/ : copdir = the name of the directory holding the opacity files
c  Note that in other routines copdir is defined to be CHARACTER*80 ; here it
c  is merely being initialized to be blank, meaning "use the current directory"
c
c /opdir/: --> data{copdir}
      character*80 copdir
      common/opdir/ copdir
      save /opdir/
c
      character*1 chpdir(80)
      equivalence (chpdir(1),copdir)
c
c COMMON /chkpoc/ : character(s) allowed to terminate a directory name
c
      character*1 cb(6)
      common/chkpoc/cb
      save /chkpoc/
c
c COMMON /alink_opal_z/ : contains data needed for smoothing routine OPALTAB
c
      common/alink_opal_z/ NTEMP,NSM,nrlow,nrhigh,RLE,t6arr(100),
     $     coff(100,nrm)
      save /alink_opal_z/
c
c COMMON /d_opal_z/ :
c  dkap = derivative of value returned by quadratic interpolation function QDER
c
c /d_opal_z/: --> data{dkap}
      common/d_opal_z/ dkap
      save /d_opal_z/
c
c COMMON /c_level_err_opal_z/ :
c  level_err = error-checking level, set by SET_ERR_CHECK; allowed values are:
c                0 : Do not check input Nzin, Zlo, Z, Zhi in READZEXCO.
c                1 : (Default): Do check these; fatal error if checks fail.
c                2 : In this case, it is also a fatal error if FEDGE = 0.
c
c /c_level_err_opal_z/: --> data{level_err}
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c___
      logical lxst
c
c COMMON /outdeb/ : controls the extent of error and debugging output, if the
c                    debug statements are un-commented by removing "c-debug;"
c                    at the beginning of the relevant lines 
c   if a log10(opacity) value is > oudebl, then the relevant debug output is
c                    performed (provided that it is not commented out)
c   else if ioudeb > 5 , then all debug output controlled by ioudeb is always
c                    performed (provided that it has not been commented out),
c                    i.e., initial-call, CO-, Z-, R-, T-, and X-interp output
c   else if ioudeb > 4 , then initial-call, Z-, R-, T-, and X-interp output
c   else if ioudeb > 3 , then initial-call, R-, T-, and X-interp output
c   else if ioudeb > 2 , then initial-call, T-, and X-interp output
c   else if ioudeb > 1 , then initial-call and X-interp output
c   else if ioudeb > 0 , then initial-call output is performed
c   koudeb counts how many debug outputs have been performed (e.g., so that
c           you can set ioudeb = -1 or increase oudebl after a given number
c           of debug outputs have been performed)
c-debug[
c-debug;      common/outdeb/ ioudeb,oudebl,koudeb
c-debug;c
c-debug;      data ioudeb/-1/,oudebl/7./,koudeb/0/
c-debug]
c
c NOTE that any lines commented out with "c-debug-chk;" correspond to cases
c  where output is NOT controlled by the variables in common/outdeb/ ; if
c  these statements are un-commented by removing "c-debug-chk;" then the
c  relevant output will occur whenever the subroutine READZEXCO is called.
c
c NOTE that lines commented out with "c-test-xdel;" correspond to code which
c  tests the effect of changing the value of xdel, the offset for logarithmic
c  X-interpolation.  There are some such cases in the subroutine OPAL, for
c  tests of interpolation among the mixes with X-values xa(m), and some in
c  the subroutine READZEXCO, for tests of interpolation among the 'GN93hz'
c  mixes with X = 0.0, 0.1, and 0.2 in order to get the mix with X = 0.03
c  (the latter produce output similar to that from "c-debug-chk;" lines).
c
c---

c					! These initializations only done once:
      if ( itimeco .ne. 12345678 ) then
c							! check some parameters
c!!         if ( nrm .ne. 19 .or. ntm .ne. 70 ) stop
         if ( nrm .ne. 28 .or. ntm .ne. 76 ) stop
c!!     $        ' STOP -- OPAL: NRM .ne. 19 or NTM .ne. 70 . '
     $        ' STOP -- OPAL: NRM .ne. 28 or NTM .ne. 76 . '
         if ( nr .lt. 6 ) stop
     $        ' STOP -- OPAL: Too few R values; NRE-NRB < 5 . '
         if ( nrb .le. 0 .or. nrb .gt. 12 .or. nre .gt. nrm ) stop
     $        ' STOP -- OPAL: NRB < 1 or NRB > 12 or NRE > NTM . '
         if ( mc .ne. mo .or. mc .ne. 8 ) stop
     $        ' STOP -- OPAL: MC .ne. MO or MC .ne. 8 . '
         if ( nt .lt. 8 .or. ntb .le. 0 .or. ntb .gt. nta(nre)-3 ) stop
     $        ' STOP -- OPAL: NT < 8 or NTB < 1 or NTB > NTA(NRE)-3. '
         if ( mx .le. 0 .or. mx .gt. 5 ) stop
     $        ' STOP -- OPAL: MX < 1 or MX > 5 . '
         if ( nz .le. 0 .or. nz .gt. nzm ) stop
     $        ' STOP -- OPAL: NZ < 1 or NZ > 14 . '
c						   ! initialize T,R values
         if ( init_trvals .le. 0 ) call get_trvals
c
c			  ! get combined Z-tabulation values available in files
         call get_zavail
c					  ! get log10( X_GH93hz + xdel ) values
         xxx_hi(1) = log10( 0.03d0 + xdel )
         do ix = 2, mx_hi
            xxx_hi(ix) = log10( xhi_in(ix) + xdel )
         enddo
c				! have now done once-only initializations
         itimeco = 12345678
c			     ! end of initializations that are done only once
      endif
c		   ! obtain the directory specification for the Gz???.x?? files
      kope = 80
      do while( kope .gt. 1 .and. copdir(kope:kope) .eq. ' ' )
         kope = kope - 1
      enddo
      if ( kope .eq. 1 .and. copdir(1:1) .eq. ' ' ) kope = 0
c								! possibly new?
      if ( kope .gt. 0 ) then
c						 ! check for [O/Fe]-file name
         if ( copdir(kope:kope) .ne. cb(1) .and.
     $        copdir(kope:kope) .ne. cb(2) ) then
            do while ( kope .gt. 1 .and.
     $           copdir(kope:kope) .ne. cb(1) .and.
     $           copdir(kope:kope) .ne. cb(2) )
               kope = kope - 1
            enddo
            if ( kope .eq. 1 .and.
     $           copdir(kope:kope) .ne. cb(1) .and.
     $           copdir(kope:kope) .ne. cb(2) ) kope = 0
            cfile_opalmixes(n_zmixes) = copdir(kope+1:)
            copdir(kope+1:) = ' '
         endif
c
      endif
c
      if ( kope .gt. 71 ) then
         write(6,10) (copdir(i:i),i=1,kope)
 10      format(
     $        ' STOP -- READCO: OPAL directory name > 71 characters:'/
     $        ' ',80a1)
         stop
      endif
c
c  NOTE that some systems return FALSE for the existence of a directory, so
c  one cannot check for the directory's existence.
c
c-dir;      if ( kope .gt. 0 ) then
c-dir;         call inqfil( copdir, lxst )
c-dir;         if ( .not. lxst ) then
c-dir;            write(6,20) (copdir(i:i),i=1,kope)
c-dir; 20         format(' STOP -- READCO: OPAL directory does not exist:'/
c-dir;     $           ' ',80a1)
c-dir;            stop
c-dir;         endif
c-dir;      endif
c		    ! just in case mx = 1 (i.e., if there is only one X-value)
      dfsx(2) = 1.d0
c
      igznotgx = 0
      mxzero = 0
      mx03 = 0
c		     ! indices of X=0 and X=.03 mixes, and X-part of file names
      do i = 1, mx
c					! loop over X-index (i, not m, here!)
         if ( xa(i) .eq. 0.0d0 ) then
            mxzero = i
            cxfil(i) = 'Gx00'
            mnofz(i) = 1
         else if ( abs(xa(i)-0.03d0) .lt. 1.d-6 ) then
            xa(i) = 0.03d0
            mx03 = i
            cxfil(i) = 'Gx03'
            mnofz(i) = 2
         else if ( abs(xa(i)-0.1d0) .lt. 1.d-6 ) then
            xa(i) = 0.1d0
            cxfil(i) = 'Gx10'
            mnofz(i) = 3
         else if ( abs(xa(i)-0.35d0) .lt. 1.d-6 ) then
            xa(i) = 0.35d0
            cxfil(i) = 'Gx35'
            mnofz(i) = 4
         else if ( abs(xa(i)-0.7d0) .lt. 1.d-6 ) then
            xa(i) = 0.7d0
            cxfil(i) = 'Gx70'
            mnofz(i) = 5
         else
            stop ' STOP -- OPAL: bad X value in array  xa(mx) . '
         endif
c					! initialize xx, for X-interpolations
         xx(i) = log10(xdel+xa(i))
         if ( i .ge. 2 ) then
            dfsx(i) = 1.d0/(xx(i)-xx(i-1))
            if ( dfsx(i) .le. 0.d0 ) stop
     $           ' STOP -- OPAL: bad X order in array  xa(mx) . '
         endif
c					! set all mixes as "not read in yet"
         itime(i) = 0
c			! have not yet read any opacities for this Z-value:
c
         if ( kz .gt. 0 .and. kz .le. nz ) then
            do mq = 1, nr
               do il = 1, nt
                  do k = 1, mo
                     do j = 1, mc
                        co(i,j,k,il,mq,kz) = badlogkval
                     enddo
                  enddo
               enddo
            enddo
         endif
c		! end of X-loop
      enddo
c
      dfsx(1) = dfsx(2)
c						! set khizat, as in READEXCO
      kzbelo = mz
      do while( kzbelo .gt. 1 .and. z .le. zval(kzbelo)-1.d-6 )
         kzbelo = kzbelo-1
      enddo
c
      khizat = min( max( khighz , 0 ) , n_zmixes )
      if ( ofebrack .eq. 0. .or. z .lt. 1.d-8 ) khizat = min(khizat,1)
      if ( khizat .eq. 1 .and. ( z .lt. 1.d-8 .or.
     $     ( z .ge. 0.01d0 .and. z .le. 0.02d0 ) .or.
     $     ( abs(zval(kzbelo)-z) .le. 1.d-6 .and. z .ge. 1.d-5 ) ) )
     $     khizat = 0
c					! should use the input [O/Fe] filename?
      if ( khizat .ge. n_zmixes ) then
c								! it exists?
         if ( cfile_opalmixes(n_zmixes) .eq. '        ' ) stop
     $        ' STOP -- READCO: no user-specified [O/Fe]-file. '
c
c			 ! obtain mix specifications for the input [O/Fe] file
         igetzxi = 1
         i_rewind = 0
         itab_dum = 0
c						     ! use copdir temporarily
         copdir(kope+1:) = cfile_opalmixes(n_zmixes)
         call open_chk_zip( iulow, copdir, i_gzip,
     $        'STOP -- Error: user-specified [O/Fe]-file not found.' )
         ifound = mixfind(iulow,n_zmixes,igetzxi,i_rewind,itab_dum,
     $        0.0d0,0.0d0,0.0d0,0.0d0)
         if ( ifound .eq. 0 ) stop
     $        ' STOP -- READCO: bad user-specified [O/Fe]-file. '
         call close_chk_zip( iulow, copdir, i_gzip )
c
c					  ! remove filename from directory name
         copdir(kope+1:) = '        '
c
      endif
c				! get mix Z-composition specifications (these
c				! will be recomputed for any mix read in later)
      do i = 1, n_zmixes
         xofe_opalmixes(i) = xiz_opalmixes(kel_o,i)
     $        /max(xiz_opalmixes(kel_fe,i),1.d-36)
         bracketofe_opalmixes(i) = log10(xofe_opalmixes(i)
     $        /xofe_opalmixes(1))
      enddo
c				! Reset current-mix data.  If no [O/Fe] shift:
      if ( khizat .le. 1 ) then
         do i = 1,nel_zmix
            xiz_mix(i) = xiz_opalmixes(i,1)
            fninz_mix(i) = fninz_opalmixes(i,1)
            bracketife_mix(i) = 0.d0
         enddo
c		! Else, if there is the [O/Fe] shift (also done in READEXCO):
      else
c		! get interpolation factors fofe (for GN93hz) and omfofe=1-fofe
c
         xofe = 10.d0**ofebrack * xofe_opalmixes(1)
         fofe = ( xiz_opalmixes(kel_o,khizat)
     $        - xofe * xiz_opalmixes(kel_fe,khizat) )
     $        / ( ( xiz_opalmixes(kel_fe,1)
     $        - xiz_opalmixes(kel_fe,khizat) ) * xofe
     $        + xiz_opalmixes(kel_o,khizat) - xiz_opalmixes(kel_o,1) )
         omfofe = 1.d0-fofe
c					! get Z-composition of interpolated mix
         do i = 1,nel_zmix
            xiz_mix(i) = fofe*xiz_opalmixes(i,1)
     $           + omfofe*xiz_opalmixes(i,khizat)
            fninz_mix(i) = fofe*fninz_opalmixes(i,1)
     $           + omfofe*fninz_opalmixes(i,khizat)
         enddo
         do i = 1,nel_zmix
            bracketife_mix(i) = log10(
     $           ( max(xiz_mix(i),1.d-36) * xiz_opalmixes(kel_fe,1) )
     $           / ( max(xiz_mix(kel_fe),1.d-36)
     $           * xiz_opalmixes(i,1) ) )
         enddo
      endif
c					! end of initializations
      return
      end
c

*------------------------------------------------------------------------------
*
      subroutine get_zavail
*     =====================
*
*  Obtain combined Z-tabulation values available in the files, and best order
*  in which to enhance the 'Gz???.x??' metallicity table.
*
*-------------------------------------------------------------------------------

      implicit double precision (a-h , o-z)


c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=28, nrb=1, nre=28,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      parameter ( zdel=0.001d0, xdel=0.03d0, xdelmin=0.001d0 )
c
      parameter ( mz=8, mz_m1=mz-1, mz_m2=mz-2, mzhi=11, mzal=13,
     $     nzm=mzal+1, nadd_zavail=nzm-mz )
c
      common/zinter_opal_z/ zvalhi(mzhi),nofz(mzhi,5,mo),mnofz(mx),
     $     zval(mz),zalval(mzal),zavail(nzm),iadd_zavail(nadd_zavail)
      save /zinter_opal_z/
c
      common/recoin_opal_z/ itimeco,mxzero,mx03,kope,igznotgx
      save /recoin_opal_z/
c___
      dimension z_l(nzm), idz(nzm)
c===					   ! (if it has not already been done):
      if ( iadd_zavail(1) .eq. 0 ) then
c
c  Get the combined Z-array zavail() by combining zval() and zalval(); this
c  should be the same as just adding the value Z=0.05 to the array zalval():
c
         k_l = 1
         k_h = 1
         k_o = 0
         k_a = 0
c
         do while ( k_l .le. mz .or. k_h .le. mzal )
c
            k_a = k_a + 1
            if ( k_a .gt. nzm ) stop
     $           ' STOP -- READCO: combined files have > 14 Z-values. '
            if ( k_l .gt. mz ) then
               zavail(k_a) = zalval(k_h)
            else if ( k_h .gt. mzal ) then
               zavail(k_a) = zval(k_l)
            else
               zavail(k_a) = min( zval(k_l) , zalval(k_h) )
            endif
            z_l(k_a) = log10( zavail(k_a) + zdel )
            idz(k_a) = 0
            if ( k_l .le. mz ) then
               if ( zval(k_l) .lt. zavail(k_a) + 1.d-6 ) then
                  idz(k_a) = k_a - k_o
                  k_o = k_a
                  k_l = k_l + 1
               endif
            endif
            if ( k_h .le. mzal ) then
               if ( zalval(k_h) .lt. zavail(k_a) + 1.d-6 )
     $              k_h = k_h + 1
            endif
c
         enddo
c
         if ( k_a .lt. nzm ) stop
     $        ' STOP -- READCO: combined files have < 14 Z-values. '
c
c  Get the best order to add values from zavail() to those of zval(), in order
c  to minimize the size of the largest interval at each step; this should
c  result in array values in iadd_zavail() of 5, 3, 13, 10, 12, 2:
c
         k_a = 0
c
         do while ( k_a .lt. nadd_zavail )
c								! next step:
            k_a = k_a + 1
c				! handle special cases where Z-range endpoints
c					! differ (this should never occur!!!):
            if ( idz(1) .eq. 0 ) then
c					! extend range to low Z, if necessary
               iadd_zavail(k_a) = 1
               k_h = 2
               do while ( k_h .lt. nzm .and. idz(k_h) .eq. 0 )
                  k_h = k_h + 1
               enddo
               if ( idz(k_h) .eq. 0 )
     $              stop ' STOP -- READCO: mz = 0 cannot happen. '
               idz(k_h) = k_h - 1
               idz(1) = 1
c
            else if ( idz(nzm) .eq. 0 ) then
c					     ! or extend to high Z if necessary
               iadd_zavail(k_a) = nzm
               k_h = nzm - 1
               do while ( k_h .gt. 1 .and. idz(k_h) .eq. 0 )
                  k_h = k_h - 1
               enddo
               if ( idz(k_h) .eq. 0 )
     $              stop ' STOP -- READCO: this REALLY cannot happen. '
               idz(nzm) = nzm - k_h
c
            else
c		      ! GENERALLY: find largest remaining subdividable interval
               k_h = 0
               dz_max = 0.d0
               do i = 2, nzm
                  if ( idz(i) .gt. 1 ) then
                     d_z = z_l(i) - z_l(i-idz(i))
                     if ( d_z .gt. dz_max ) then
                        dz_max = d_z
                        k_h = i
                     endif
                  endif
               enddo
               if ( k_h .eq. 0 )
     $              stop ' STOP -- READCO: k_h = 0 cannot happen. '
c
c					   ! find best subdivision of interval
               if ( idz(k_h) .eq. 2 ) then
                  k_l = k_h - 2
                  k_o = k_h - 1
               else
                  k_l = k_h - idz(k_h)
                  dz_max = 0.d0
                  k_o = 0
                  do i = k_l + 1, k_h - 1
                     d_z = ( z_l(k_h) - z_l(i) )
     $                    / ( z_l(i) - z_l(k_l) )
                     if ( d_z .gt. 1.d0 ) d_z = 1.d0 / d_z
                     if ( d_z .gt. dz_max ) then
                        dz_max = d_z
                        k_o = i
                     endif
                  enddo
                  if ( k_o .eq. 0 )
     $                 stop ' STOP -- READCO: k_o = 0 cannot happen. '
               endif
c						! store this subdivision
               iadd_zavail(k_a) = k_o
               idz(k_o) = k_o - k_l
               idz(k_h) = k_h - k_o
c
            endif
c
         enddo
c
      endif
c
      return
      end
c


*------------------------------------------------------------------------------
*
      subroutine get_trvals
*     =====================
*------------------------------------------------------------------------------

      implicit double precision (a-h , o-z)


c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=28, nrb=1, nre=28,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      parameter ( nrdel=nrb-1, ntdel=ntb-1 )
c
c PARAMETERS: positions of DlogT-change in table, low logT and logR values:
c
c!!      parameter ( ks81=ntm-3, ks83=ks81+1, ks60=ks81-21, ks61=ks60+1,
c!!     $     alrlo=-8.0d0, flogtlo=3.75d0, flogt60=6.0d0, flogt81=8.1d0 )
      parameter ( ks81=ntm-9, ks83=ks81+1, ks60=ks81-21, ks61=ks60+1,
     $     alrlo=-9.0d0, flogtlo=3.75d0, flogt60=6.0d0, flogt81=8.1d0 )
c
c COMMON /a_opal_z/ : matrices for opacity storage and interpolation:
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime(mx),cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
c COMMON /b_opal_z/ : high and low logT6 and logR limits, and mix Z-values:
c
      common/b_opal_z/ nta(0:nrm_p1),zz(mx,nz),ntax0(0:nrm),
     $     ntax03(0:nrm), sltlo, slthi, dltlo_inv, dlthi_inv,
     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
      save /b_opal_z/
c
c COMMON /alink_opal_z/ : contains data needed for smoothing routine OPALTAB
c
      common/alink_opal_z/ NTEMP,NSM,nrlow,nrhigh,RLE,t6arr(100),
     $     coff(100,nrm)
      save /alink_opal_z/
c===					! only initialize once
      if ( init_trvals .gt. 0 ) return
c
      init_trvals = 1
c				! intialize logR and 1/delta(logR) values
      do i = 1,nrm
         alrf(i) = (i-1)*0.5d0+alrlo
      enddo
      do i = 1,nr
         alr(i) = alrf(i+nrdel)
         dfsr(i) = 2.d0
      enddo
c				! intialize logT, logT6, T6, and 1/delta(logT6)
      flogtin(1) = flogtlo
      do i = 2,ks60
         flogtin(i) = (i-1)*0.05d0+flogtlo
         if ( i .ge. ntb ) dfs(i-ntdel) = 20.d0
      enddo
      if ( abs(flogtin(ks60)-flogt60) .gt. 5.d-6 ) stop
     $     ' STOP -- READCO: initialization error. '
      flogtin(ks60) = flogt60
      do i = ks61,ks81
         flogtin(i) = (i-ks60)*0.1d0+flogt60
         if ( i .ge. ntb ) dfs(i-ntdel) = 10.d0
      enddo
      if ( abs(flogtin(ks81)-flogt81) .gt. 5.d-6 ) stop
     $     ' STOP -- READCO: initialization error. '
      flogtin(ks81) = flogt81
      do i = ks83,ntm
         flogtin(i) = (i-ks81)*0.2d0+flogt81
         if ( i .ge. ntb ) dfs(i-ntdel) = 5.d0
      enddo
      do i=1,ntm
         t6arr(i) = 10.d0**(flogtin(i)-6.d0)
      enddo
      do i=1,nt
         alt(i) = flogtin(i+ntdel)-6.d0
         t6list(i) = t6arr(i+ntdel)
      enddo
c
c-done-above;      do i = 2,nt
c-done-above;         dfs(i) = 1./(alt(i)-alt(i-1))
c-done-above;      enddo
c-done-above;      do i = 2,nr
c-done-above;         dfsr(i) = 1./(alr(i)-alr(i-1))
c-done-above;      enddo
c					  ! For extrapolation at low R and T6
      dfsr(1) = dfsr(2)
      dfs(1) = dfs(2)
c						      ! R-extrapolation limits
      slrlo = alr(1)
      slrhi = alr(nr)
      dlrlo_inv = dfsr(1)
      dlrhi_inv = dfsr(nr)
c						       ! T-extrapolation limits
      sltlo = alt(1)
      slthi = alt(nt)
      dltlo_inv = dfs(1)
      dlthi_inv = dfs(nt)
c
      return
      end


*
***-----------------------------------------------------------------------------
*
      subroutine set_opal_dir( cdirin )
*     =================================
*
*     The input character variable  cdirin  can be used to specify the 
*     directory where the OPAL opacity files will subsequently be looked for 
*     (default is the local directory, which can also be specified by supplying
*     a blank argument to SET_OPAL_DIR). Note that the total length of the 
*     name MUST NOT exceed 71 characters.
*
***-----------------------------------------------------------------------------
      character*(*) cdirin

      common/recoin_opal_z/ itimeco,mxzero,mx03,kope,igznotgx
      save /recoin_opal_z/

      character*80 copdir
      common/opdir/ copdir
      save /opdir/

      character*1 cb(6)
      common/chkpoc/cb
      save /chkpoc/
c---
      character*255 cnamenew
      character*1 chname(255)
      equivalence (chname(1),cnamenew)
c
      logical lxst
c===
      cnamenew = cdirin

      kope = 255
      do while ( kope .gt. 1 .and. chname(kope) .eq. ' ' )
         kope = kope - 1
      enddo

      if ( kope .gt. 1 .and. chname(1) .eq. ' ' ) then
         idel = 1
         do while ( idel .lt. kope .and. chname(idel) .eq. ' ' )
            idel = idel + 1
         enddo
         idel = idel - 1
         do i = idel + 1, kope
            chname(i-idel) = chname(i)
            chname(i) = ' '
         enddo
         kope = kope - idel
      endif

      if ( kope .eq. 1 .and. chname(kope) .eq. ' ' ) then
         kope = 0
      else
         do i = 1, kope
            if ( chname(i) .eq. ' ' ) stop
     $           ' STOP -- SET_OPAL_DIR: blanks in directory name. '
         enddo
         if ( kope .lt. 80 .and. chname(kope) .ne. cb(1) .and.
     $        chname(kope) .ne. cb(2) ) then
            kope = kope + 1
            chname(kope) = cb(1)
         endif
      endif

      if ( kope .gt. 71 ) then
         write(6,10) (chname(i),i=1,kope)
 10      format(' STOP -- SET_OPAL_DIR:',
     $        ' OPAL directory name > 71 characters:'/' ',80a1)
         stop
      endif

      copdir = cnamenew(1:80)
c
c  NOTE that some systems return FALSE for the existence of a directory, so
c  one cannot check for the directory's existence.
c
c-dir;      if ( kope .gt. 0 ) then
c-dir;         call inqfil( copdir, lxst )
c-dir;         if ( .not. lxst ) then
c-dir;            write(6,20) (chname(i),i=1,kope)
c-dir; 20         format(' STOP -- SET_OPAL_DIR:',
c-dir;     $           ' OPAL directory does not exist:'/' ',80a1)
c-dir;            stop
c-dir;         endif
c-dir;      endif
c
      return
      end
cc
cc******************************************************************************
cc
c      subroutine set_ofe_file( cfileofe )
cc     ===================================
cc
c      character*(*) cfileofe
cc
c      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
cc
c      character*2 cel_opalmixes(nel_zmix)
c      character*8 cfile_opalmixes(n_zmixes)
c      common/opalmixes/ xiz_mix(nel_zmix),fninz_mix(nel_zmix),
c     $     bracketife_mix(nel_zmix),bracketofe_opalmixes(n_zmixes),
c     $     xofe_opalmixes(n_zmixes),xiz_opalmixes(nel_zmix,n_zmixes),
c     $     fninz_opalmixes(nel_zmix,n_zmixes),
c     $     cel_opalmixes,cfile_opalmixes
c      save /opalmixes/
cc
c      character*1 cb(6)
c      common/chkpoc/cb
c      save /chkpoc/
cc---
c      character*255 cnamenew
c      character*1 chname(255)
c      equivalence (chname(1),cnamenew)
cc===
c      cnamenew = cfileofe
cc
c      last = 255
c      do while ( last .gt. 1 .and. chname(last) .eq. ' ' )
c         last = last - 1
c      enddo
cc
c      if ( last .gt. 1 .and. chname(1) .eq. ' ' ) then
c         idel = 1
c         do while ( idel .lt. last .and. chname(idel) .eq. ' ' )
c            idel = idel + 1
c         enddo
c         idel = idel - 1
c         do i = idel + 1, last
c            chname(i-idel) = chname(i)
c            chname(i) = ' '
c         enddo
c         last = last - idel
c      endif
cc
c      cfile_opalmixes(n_zmixes) = cnamenew(1:8)
cc
c      return
c      end
c
cc******************************************************************************
cc
c      subroutine set_xhi( kxhi )
cc     ==========================
cc
c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
cc
c      parameter ( mx_hi=2*mx, mo_m1=mo-1, mx_hi_nz=mx_hi*nz )
cc
c      common /xhi_opal_z/ xhi_in(mx_hi), xhi_use(mx_hi,nz),
c     $     xxx_hi(mx_hi), nx_hi(nz), ireq_hi(mx_hi), khighx(nz),
c     $     kavail_xhi, kuse_xhi, kdo_xhi
c      save /xhi_opal_z/
cc						! set high-X flag value
c      kuse_xhi = max( 0 , min( 2 , kxhi ) )
c      if ( kavail_xhi .le. 0 ) then
c         kdo_xhi = 0
c      else
c         kdo_xhi = kuse_xhi
c      endif
cc
c      return
c      end
cc
cc******************************************************************************
cc
c      subroutine ask_xhi( kxhi, kavail )
cc     ==================================
cc
c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
cc
c      parameter ( mx_hi=2*mx, mo_m1=mo-1, mx_hi_nz=mx_hi*nz )
cc
c      common /xhi_opal_z/ xhi_in(mx_hi), xhi_use(mx_hi,nz),
c     $     xxx_hi(mx_hi), nx_hi(nz), ireq_hi(mx_hi), khighx(nz),
c     $     kavail_xhi, kuse_xhi, kdo_xhi
c      save /xhi_opal_z/
cc						! return high-X flag value
c      kxhi = kuse_xhi
cc						! return availability flag
c      kavail = max( 0 , min( 1 , kavail_xhi ) )
cc
c      return
c      end
cc
cc******************************************************************************
cc
c      subroutine set_err_check( level )
cc     =================================
cc						! set error-checking level
c      common /c_level_err_opal_z/ level_err
c      save /c_level_err_opal_z/
cc
c      level_err = max( 0 , min( 2 , level ) )
cc
c      return
c      end
cc
cc******************************************************************************
cc
c      subroutine ask_err_check( level )
cc     =================================
cc						! return error-checking level
c      common /c_level_err_opal_z/ level_err
c      save /c_level_err_opal_z/
cc
c      level = level_err
cc
c      return
c      end
cc
cc******************************************************************************
cc
c      subroutine set_logt6_limits( vlo, dvlo, vhi, dvhi )
cc     ===================================================
cc
c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
cc
c      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
c     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
c     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
c     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime(mx),cxdf(mx,mo,nz),
c     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
c     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
c     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
c     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
c     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
c      save /a_opal_z/
cc
c      common/b_opal_z/ nta(0:nrm_p1),zz(mx,nz),ntax0(0:nrm),
c     $     ntax03(0:nrm), sltlo, slthi, dltlo_inv, dlthi_inv,
c     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
c      save /b_opal_z/
cc===						! set T,R-values, if necessary
c      if ( init_trvals .le. 0 ) call get_trvals
cc
cc  Set the logT6 limits, according to the input values:
cc
c      if ( vlo .gt. -90. ) sltlo = max( alt(1) , min( alt(nt) , vlo ) )
cc
c      if ( dvlo .lt. -90. ) then
c         dltlo_inv = max( 1. / ( sltlo - alt(1) + 1. / dfs(1) ) ,
c     $        dltlo_inv )
c      else if ( dvlo .lt. 0. ) then
c         dltlo_inv = dfs(1)
c      else
c         dltlo_inv = 1. / min( sltlo - alt(1) + 1. / dfs(1) ,
c     $        max( dvlo , 1.e-6 ) )
c      endif
cc
c      if ( vhi .gt. -90. ) slthi = max( alt(1) , min( alt(nt) , vhi ) )
cc
c      if ( dvhi .lt. -90. ) then
c         dlthi_inv = max( 1. / ( alt(nt) - slthi + 1. / dfs(nt) ) ,
c     $        dlthi_inv )
c      else if ( dvhi .lt. 0. ) then
c         dlthi_inv = dfs(nt)
c      else
c         dlthi_inv = 1. / min( alt(nt) - slthi + 1. / dfs(nt) ,
c     $        max( dvhi , 1.e-6 ) )
c      endif
cc
c      return
c      end
cc
cc******************************************************************************
cc
c      subroutine set_logr_limits( vlo, dvlo, vhi, dvhi )
cc     ==================================================
cc
c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
cc
c      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
c     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
c     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
c     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime(mx),cxdf(mx,mo,nz),
c     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
c     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
c     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
c     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
c     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
c      save /a_opal_z/
cc
c      common/b_opal_z/ nta(0:nrm_p1),zz(mx,nz),ntax0(0:nrm),
c     $     ntax03(0:nrm), sltlo, slthi, dltlo_inv, dlthi_inv,
c     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
c      save /b_opal_z/
cc===						! set T,R-values, if necessary
c      if ( init_trvals .le. 0 ) call get_trvals
cc
cc  Set the logR limits, according to the input values:
cc
c      if ( vlo .gt. -90. ) slrlo = max( alr(1) , min( alr(nr) , vlo ) )
cc
c      if ( dvlo .lt. -90. ) then
c         dlrlo_inv = max( 1. / ( slrlo - alr(1) + 1. / dfsr(1) ) ,
c     $        dlrlo_inv )
c      else if ( dvlo .lt. 0. ) then
c         dlrlo_inv = dfsr(1)
c      else
c         dlrlo_inv = 1. / min( slrlo - alr(1) + 1. / dfsr(1) ,
c     $        max( dvlo , 1.e-6 ) )
c      endif
cc
c      if ( vhi .gt. -90. ) slrhi = max( alr(1) , min( alr(nr) , vhi ) )
cc
c      if ( dvhi .lt. -90. ) then
c         dlrhi_inv = max( 1. / ( alr(nr) - slrhi + 1. / dfsr(nr) ) ,
c     $        dlrhi_inv )
c      else if ( dvhi .lt. 0. ) then
c         dlrhi_inv = dfsr(nr)
c      else
c         dlrhi_inv = 1. / min( alr(nr) - slrhi + 1. / dfsr(nr) ,
c     $        max( dvhi , 1.e-6 ) )
c      endif
cc
c      return
c      end
cc
cc******************************************************************************
cc
c      subroutine reset_z_limits( vlo, dvlo, vhi, dvhi )
cc     =================================================
cc
c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
cc
c      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
c     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
c     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
c     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime(mx),cxdf(mx,mo,nz),
c     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
c     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
c     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
c     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
c     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
c      save /a_opal_z/
cc				! check whether opacities have been read in
c      igot = 0
c      do i = 1, mx
c         if ( itime(i) .eq. 12345678 ) igot = 1
c      enddo
cc				! if no opacities were read in, cannot reset!
c      if ( igot .eq. 0 ) return
cc
c      if ( min( vlo , vhi ) .gt. -1.e-6 ) then
c         zlow = max( zsto(1) , min( vlo , vhi , zsto(numz) ) )
c         zhigh = min( zsto(numz) , max( vhi , vlo , zsto(1) ) )
c      else if ( vlo .gt. -1.e-6 ) then
c         zlow = max( zsto(1) , min( vlo , zhigh , zsto(numz) ) )
c      else if ( vhi .gt. -1.e-6 ) then
c         zhigh = min( zsto(numz) , max( vhi , zlow , zsto(1) ) )
c      endif
cc
c      if ( dvlo .gt. -1.e-6 ) then
c         if ( zlow .le. zsto(1) + zacc(1) ) then
c            zlo_ex = zlow - max( dvlo , zacc(1) )
c         else
c            zlo_ex = zlow - max( dvlo , zacc(numz) )
c         endif
c      else if ( zlow .le. zsto(1) + zacc(1) ) then
c         zlo_ex = min( zlo_ex , zlow - zacc(1) )
c      else
c         zlo_ex = min( zlo_ex , zlow - zacc(numz) )
c      endif
cc
c      if ( dvhi .gt. -1.e-6 ) then
c         zhi_ex = zhigh + max( dvhi , zacc(numz) )
c      else
c         zhi_ex = max( zhi_ex , zhigh + zacc(numz) )
c      endif
cc
c      return
c      end
cc
cc******************************************************************************
cc
c      subroutine ask_logt6_limits( vlo, dvlo, vhi, dvhi )
cc     ===================================================
cc
c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
cc
c      common/b_opal_z/ nta(0:nrm_p1),zz(mx,nz),ntax0(0:nrm),
c     $     ntax03(0:nrm), sltlo, slthi, dltlo_inv, dlthi_inv,
c     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
c      save /b_opal_z/
cc===						! set T,R-values, if necessary
c      if ( init_trvals .le. 0 ) call get_trvals
cc
c      vlo = sltlo
c      dvlo = 1. / dltlo_inv
c      vhi = slthi
c      dvhi = 1. / dlthi_inv
cc
c      return
c      end
cc
cc******************************************************************************
cc
c      subroutine ask_logr_limits( vlo, dvlo, vhi, dvhi )
cc     ==================================================
cc
c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
cc
c      common/b_opal_z/ nta(0:nrm_p1),zz(mx,nz),ntax0(0:nrm),
c     $     ntax03(0:nrm), sltlo, slthi, dltlo_inv, dlthi_inv,
c     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
c      save /b_opal_z/
cc===						! set T,R-values, if necessary
c      if ( init_trvals .le. 0 ) call get_trvals
cc
c      vlo = slrlo
c      dvlo = 1. / dlrlo_inv
c      vhi = slrhi
c      dvhi = 1. / dlrhi_inv
cc
c      return
c      end
cc
cc******************************************************************************
cc
c      subroutine ask_z_limits( nzmax, zmin, zmax )
cc     ============================================
cc
cc  Returns NZ (maximum allowed stored Z-values) and Z-limits (0.0 and 0.1)
cc
c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
cc===
c      nzmax = nz
c      zmin = 0.0
c      zmax = 0.1
cc
c      return
c      end
cc
cc******************************************************************************
cc
c      subroutine ask_z_use( nzuse, zlo, zmid, zhi, zloex, zhiex )
cc     ===========================================================
cc
cc  Returns the current values of numz, zlow, zmiddle, zhigh, zlo_ex, zhi_ex
cc
c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
cc
c      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
c     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
c     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
c     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime(mx),cxdf(mx,mo,nz),
c     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
c     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
c     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
c     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
c     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
c      save /a_opal_z/
cc===
c      igot = 0
c      do i = 1, mx
c         if ( itime(i) .eq. 12345678 ) igot = 1
c      enddo
cc
c      if ( igot .eq. 0 ) then
c         nzuse = 0
c         zlo = 0.
c         zmid = 0.02
c         zhi = 0.1
c         zloex = -1.e-8
c         zhiex = 0.12
c      else
c         nzuse = numz
c         zlo = zlow
c         zmid = zmiddle
c         zhi = zhigh
c         zloex = zlo_ex
c         zhiex = zhi_ex
c      endif
cc
c      return
c      end
cc
cc******************************************************************************
cc
c      subroutine ask_z_array( kzstart, karraystart, zarray, narray )
cc     ==============================================================
cc
cc  Returns Z-values from zsto(), starting with element kzstart, in the
cc  array zarray(narray), starting with element karraystart; any excess
cc  elements in zarray() after nums is reached are filled with values of -1.
cc
c      dimension zarray(narray)
cc
c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
cc
c      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
c     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
c     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
c     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime(mx),cxdf(mx,mo,nz),
c     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
c     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
c     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
c     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
c     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
c      save /a_opal_z/
cc===
c      igot = 0
c      do i = 1, mx
c         if ( itime(i) .eq. 12345678 ) igot = 1
c      enddo
cc
c      if ( igot .eq. 0 ) then
c         nzuse = 0
c      else
c         nzuse = numz
c      endif
cc
c      k_z = max( kzstart , 1 )
c      k_a = karraystart
c      do while ( k_z .le. nzuse .and. k_a .le. narray )
c         zarray(k_a) = zsto(k_z)
c         k_a = k_a + 1
c         k_z = k_z + 1
c      enddo
c      do while ( k_a .le. narray )
c         zarray(k_a) = -1.
c         k_a = k_a + 1
c      enddo
cc
c      return
c      end
cc
cc******************************************************************************
cc
c      subroutine set_smooth( initsmooth, lowCOsmooth, interpCOsmooth )
cc     ================================================================
cc
c      common/c_opal_ctrl_smooth/ init_smo, low_CO_smo, interp_CO_smo
c      save /c_opal_ctrl_smooth/
cc
c      if ( initsmooth .ge. 0 ) init_smo = min( initsmooth , 1 )
cc
c      if ( lowCOsmooth .ge. 0 ) low_CO_smo = min( lowCOsmooth , 1 )
cc
c      if ( interpCOsmooth .ge. 0 )
c     $     interp_CO_smo = min( interpCOsmooth , 1 )
cc
c      return
c      end
cc
cc******************************************************************************
cc
c      subroutine ask_smooth( initsmooth, lowCOsmooth, interpCOsmooth )
cc     ================================================================
cc
c      common/c_opal_ctrl_smooth/ init_smo, low_CO_smo, interp_CO_smo
c      save /c_opal_ctrl_smooth/
cc
c      initsmooth = init_smo
cc
c      lowCOsmooth = low_CO_smo
cc
c      interpCOsmooth = interp_CO_smo
cc
c      return
c      end
cc

*
***----------------------------------------------------------------------------
c
      subroutine readco(z,kallrd,khighz,iulow)
c     ========================================
c
c..... The purpose of this subroutine is to read the data tables; actually,
c      it just calls READEXCO to do the work, setting [O/Fe] = 0.0
c
c Z is the metallicity; opacities will be interpolated (quadratically) in Z if
c   necessary, with values of Z from 0.0 to 0.1 being allowed.
c if kallrd = 0 , then only opacities for m-th X value are read in;
c   else, opacities for all mx X values are read in.
c if khighz = 0 , then the file GN93hz is not used;
c   else, it may be used for the case C=O=0.0 to get improved opacities.
c iulow is the unit number from which the lowest Z-value files are read;
c   units iulow thru iulow+3 may be needed.
c
***----------------------------------------------------------------------------

      implicit double precision (a-h, o-z)

      call readexco(z,kallrd,min(khighz,1),iulow,0.0d0)
c
      return
      end


*
***----------------------------------------------------------------------------
*
      subroutine readexco(z,kallrd,khighz,iulow,ofebrack)
c     ===================================================
c
c..... The purpose of this subroutine is to read the data tables
c
c Z is the metallicity; opacities will be interpolated (quadratically) in Z if
c   necessary, with values of Z from 0.0 to 0.1 being allowed.
c if kallrd = 0 , then only opacities for m-th X value are read in;
c   else, opacities for all mx X-values are read in.
c if khighz < 0 , behave the same as for a value of 4 (or 1, if ofebrack = 0.0)
c if khighz = 0 , then the file GN93hz is not used;
c           = 1 or more, then GN93hz is used if it will improve Z-interpolation
c           = 2 , then Alrd96a2 and GN93hz are used to interpolate in [O/Fe]
c           = 3 , then C95hz and GN93hz are used to interpolate in [O/Fe]
c           = 4 , then W95hz and GN93hz are used to interpolate in [O/Fe]
c           = 5 , then GN93hz and a user-specified file (appearing after the
c                  directory specification in the character variable copdir)
c                  are used to interpolate in [O/Fe]
c iulow is the unit number from which the lowest Z-value files are read;
c   units iulow thru iulow+3 may be needed.
c ofebrack is [O/Fe], logarithmic oxygen (or alpha-element) enhancement factor,
c   relative to the Sun: ofebrack = log10{ (O_mix/Fe_mix) / (O_sun/Fe_sun) } .
c   If khighz is 0 or 1, then ofebrack is ignored; otherwise, interpolate
c   logKappa linearly between mix 1 and mix khighz, the interpolation factors
c   being such as to yield the desired [O/Fe] from combining these mixes.
c   Note that GN93hz has [O/Fe] = 0.0, Alrd96a2 has [O/Fe] = 0.3, C95hz has
c   [O/Fe] = 0.4, and W95hz has [O/Fe] = 0.5, but they have different patterns
c   of enhancements for elements other than oxygen.
c
***-----------------------------------------------------------------------------

      implicit double precision (a-h, o-z)


c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=28, nrb=1, nre=28,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime(mx),cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c===
      numz = 1
      if ( z .lt. -1.d-6 ) then
         zat = 0.02d0
      else if ( z .lt. 1.d-8 ) then
         zat = 0.d0
      else
         zat = z
      endif
c
      call read_kz(1,zat,kallrd,khighz,iulow,ofebrack)
c
      zlow = zat - zacc(1)
      zmiddle = zat
      zhigh = zat + zacc(1)
      zlo_ex = zlow - zacc(1)
      zhi_ex = zhigh + zacc(1)
      dfsz(1) = 1.d0
c
      return
      end


*
***------------------------------------------------------------------------
*
      subroutine readzexco(nzin,zlo,z,zhi,khighz,iulow,ofebrack)
*     ==========================================================
*
*     This subroutine is used to readin the OPAL opacity files, allowing the
*     user to control whether and how opacities will subsequently be 
*     interpolated in Z. 
*  The controlling input variables are:
*
*    Nzin      INTEGER: the number of metallicity values Z_i to be stored,
*                for subsequent Z-interpolation when OPAL or OPAC is called;
*                this is discussed just below.
*    Zlo       DOUBLE-PRECISION REAL: the lowest metallicity value that will
*                be required; this is discussed just below.
*    Z         DOUBLE-PRECISION REAL: the "typical" or central metallicity
*                value; this is discussed just below.
*    Zhi       DOUBLE-PRECISION REAL: the highest metallicity value that will
*                be required; this is discussed just below.
*    khighz    INTEGER: controls the use of the C=O=0.0 opacity file 'GN93hz',
*                and of the similar files having [O/Fe] > 0.0:
*                 khighz < 0: treat the same as  khighz = 4  below (or, for the
*                              case when ofebrack = 0.0, treat as khighz = 1).
*                 khighz = 0: use of the file 'GN93hz' is prevented; only for
*                              X < 0.75 is accurate X-interpolation available.
*                 khighz = 1: the file 'GN93hz' is used to obtain opacities for
*                              the C=O=0.0 mixes (i.e., opacities with better
*                              Z-interpolation), including the added X-values
*                              X={0.2,0.5,0.8,0.9,0.95,1-Z} (i.e., allowing
*                              accurate X-interpolation up to X = 1-Z); for
*                              the mixes with C+O > 0.0, corresponding opacity
*                              shifts are applied, for consistency.
*                 khighz = 2: file 'Alrd96a2' with [O/Fe] = 0.3 \  is used in
*                 khighz = 3: file 'C95hz' with [O/Fe] = 0.4     } addition to
*                 khighz = 4: file 'W95hz' with [O/Fe] = 0.5    /  'GN93hz',
*                              if READZEXCO was called with a non-zero value of
*                              ofebrack, in order to interpolate in the excess
*                              oxygen/alpha-element enrichment [O/Fe].
*                 khighz = 5: the name of a file with non-zero [O/Fe] must have
*                              been set already, by calling the subroutine
*                              SET_OFE_FILE described below; it will be used
*                              instead of 'Alrd96a2', 'C95hz', or 'W95hz' when
*                              interpolating in [O/Fe] (its [O/Fe] value will
*                              be computed when it is read in; if it actually
*                              has [O/Fe] = 0.0, the resulting behavior is not
*                              defined and will surely be erroneous).
*
*     default value used for khighz: 1
*
*    iulow     INTEGER: the beginning Fortran unit number for reading opacity
*                files; Fortran units  iulow  through  iulow + 3  may be used.
*                (Note: unless the use explicitly calls READZEXCO, READEXCO, or
*                READCO, the default-setup call to READZEXCO in OPAL will be
*                executed; it has iulow = 23).
*    ofebrack  DOUBLE-PRECISION REAL: the value of [O/Fe], the logarithmic
*                oxygen (or alpha-element) enhancement factor, relative to the
*                Sun:  ofebrack = log10{ (O_mix/Fe_mix) / (O_sun/Fe_sun) } ,
*                where O_mix, Fe_mix, O_sun, and Fe_sun are number densities.
*                If  khighz  is 0 or 1, then ofebrack is ignored; otherwise,
*                READZEXCO interpolates (or extrapolates) log(Kappa) linearly
*                between mix 1 and mix  khighz , the interpolation factors
*                being such as to yield the desired [O/Fe] from combining these
*                mixes.  Note that 'GN93hz' has [O/Fe] = 0.0 by definition,
*                'Alrd96a2' has [O/Fe] = 0.3, 'C95hz' has [O/Fe] = 0.4, and
*                'W95hz' has [O/Fe] = 0.5, but they have different patterns of
*                enhancements for elements other than oxygen; their elemental
*                abundances and the corresponding opacity shifts are discussed
*                further below.
*
*
***-------------------------------------------------------------------------

      implicit double precision (a-h, o-z) 
      
c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=28, nrb=1, nre=28,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime(mx),cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
      parameter ( zdel=0.001d0, xdel=0.03d0, xdelmin=0.001d0 )
c
      parameter ( mz=8, mz_m1=mz-1, mz_m2=mz-2, mzhi=11, mzal=13,
     $     nzm=mzal+1, nadd_zavail=nzm-mz )
c
      common/zinter_opal_z/ zvalhi(mzhi),nofz(mzhi,5,mo),mnofz(mx),
     $     zval(mz),zalval(mzal),zavail(nzm),iadd_zavail(nadd_zavail)
      save /zinter_opal_z/
c
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c___
      dimension z_use(nzm)
c===
c
c  Get the number of Z-values and the "typical" value of Z
c
      if ( level_err .gt. 0 .and.
     $     ( nzin .le. 0 .or. nzin .gt. nz ) ) then
         write(6,10) nzin, nz
 10      format(' '/' STOP -- READZEXCO: bad Nzin =',i12,
     $        ' (should lie in range 1 to',i3,').')
         stop ' STOP -- READZEXCO: Error: bad Nzin value. '
      endif
c
      numz = max( 1 , min( nzin , nz , nzm ) )
c
      if ( z .ge. -1.d-6 ) then
         zat = z
      else if ( zlo .gt. -1.d-6 .and. zhi .ge. zlo ) then
         if ( nzin .eq. 2 ) then
            zat = ( max( zlo , 0.0d0 ) + max( zhi , 0.0d0 ) ) * 0.5d0
         else
            zat = exp( ( log( max( zlo , 0.0d0 ) + zdel )
     $           + log( max( zhi , 0.0d0 ) + zdel ) ) * 0.5d0 ) - zdel
         endif
         if ( zat .lt. 1.d-8 ) zat = 0.d0
      else if ( max( zlo , zhi ) .gt. -1.d-6 ) then
         zat = max( zlo , zhi , 0.0d0 )
      else
         zat = 0.02d0
      endif
      if ( zat .lt. 1.d-8 ) zat = 0.d0
c
c  If there is only one Z-value to store, call readexco instead, and return:
c

      if ( numz .eq. 1 ) then
c
         call readexco(zat,1,khighz,iulow,ofebrack)
c
         if ( nzin .eq. 1 ) then
            if ( zlo .ge. -1.d-6 ) then
               zlow = min( zlow , zlo )
               zlo_ex = min( zlo_ex , max( 2.d0 * zlow - zmiddle , 
     $              0.d0 ) )
            endif
            zhigh = max( zhigh , zhi )
            zlo_ex = min( zlo_ex , zlow - zacc(1) )
            zhi_ex = max( zhi_ex , 2.d0 * zhigh - zmiddle ,
     $           zhigh + zacc(1) )
         endif
c
         return
c
      endif
c
c  If there is more than one Z-value to store:
c					      ! make sure of available Z-values
      call get_zavail
c								! check Z > 0.1
      if ( max( zlo , zat , zhi ) .gt. zavail(nzm) + 1.d-6 ) then
         write(6,20) nzin, zlo, z, zhi, 'Z > 0.1 is not allowed!!!'
 20      format(' '/' STOP -- READZEXCO(Nzin=',i2,
     $        '): bad Z values',3f12.8/'   ',a)
         stop ' STOP -- READZEXCO: Error: Z > 0.1 is not allowed!!! '
      endif
c					! check for Zlo, Z, Zhi out of order
      zlo_t = max( zlo , 0.0d0 )
      zhi_t = max( zhi , 0.0d0 )
      if ( zlo .lt. -1.d-6 ) zlo_t = zat
      if ( zhi .lt. -1.d-6 ) zhi_t = zat
      if ( level_err .gt. 0 .and. ( zat .gt. zhi_t + 1.d-6 .or.
     $     zlo_t .gt. min( zat , zhi_t ) + 1.d-6 ) ) then
         write(6,20) nzin, zlo, z, zhi,
     $        'Zlo, Z, Zhi should be in increasing order'
         stop ' STOP -- READZEXCO: Error: bad Z values. '
      endif
c
      zmiddle = min( zat , zavail(nzm) )
c
c  If there should be two Z-values, then handle this case and return:

      if ( numz .eq. 2 ) then
c							! error checking
         if ( level_err .gt. 0 .and.
     $        max( zlo , zhi ) .gt. -1.d-6 ) then
            if ( min( zlo , zhi ) .gt. -1.d-6 ) then
               if ( zhi_t - zlo_t .lt. 1.d-5 .or.
     $              zlo_t * 1.01d0 .gt. zhi_t + 1.d-6 ) then
                  write(6,20) nzin, zlo, z, zhi,
     $                 'Zlo and Zhi too close together for Nzin=2'
                  stop ' STOP -- READZEXCO: Error: bad Z values. '
               endif
            endif
            if ( zhi_t .lt. 1.d-8 ) then
               rat_t = 0.d0
            else
               rat_t = zlo_t / zhi_t
            endif
            if ( rat_t .lt. 0.6d0 .and.
     $           zhi_t - zlo_t .gt. 0.00020001d0 ) then
               write(6,20) nzin, zlo, z, zhi,
     $              'Zlo and Zhi too far apart for Nzin=2'
               stop ' STOP -- READZEXCO: Error: bad Z values. '
            endif
         endif
c
c  (will be linear interpolation in log[Z+0.001]); default range is plus/minus
c  10 percent in Z or 2.e-5, minimum range is at least 1 percent or 1.e-5
c
         if ( nzin .gt. 2 .or. zlo .lt. -1.d-6 .or. zlo .gt. zat ) then
            if ( nzin .gt. 2 .or. zhi .lt. zat ) then
               zlow = max( 0.0d0 , min( 0.9d0 * zmiddle ,
     $              zmiddle - 1.d-5 , 0.8182d0 * zavail(nzm) ) )
               zhigh = min( zavail(nzm) , max( 1.1d0 * zmiddle ,
     $              zmiddle + 1.d-5 , 1.1d0 * zlow / 0.9d0 ) )
            else
               zhigh = max( zmiddle , 1.d-5 ,
     $              min( zhi , zavail(nzm) ) )
               zlow = max( 0.0d0 , min( zmiddle ,
     $              0.9d0 * zhigh / 1.1d0 , zhigh - 2.d-5 ) )
            endif
         else if ( zhi .lt. zat ) then
            zlow = max( 0.0d0 ,
     $           min( zlo , zmiddle , zavail(nzm) - 0.01d0 ) )
            zhigh = min( zavail(nzm) , max( zmiddle ,
     $           1.1d0 * zlow / 0.9d0 , zlow + 2.d-5 ) )
         else if (zhi - zlo .lt. 1.d-5 .or. zlo * 1.01d0 .gt. zhi) then
            zlow = max( 0.0d0 , min( zlo , zmiddle - 5.d-6 ,
     $           zmiddle / 1.005d0 , zavail(nzm) - 0.01d0 ) )
            zhigh = min( zavail(nzm) , max( zhi , zmiddle + 5.d-6 ,
     $           zmiddle * 1.005d0 , 1.d-5 ) )
         else
            zlow = max( 0.0d0 , min( zlo , zavail(nzm) - 0.01d0 ) )
            zhigh = min( zavail(nzm) , zhi )
         endif
c
         zat = zlow
         call read_kz(1,zat,1,khighz,iulow,ofebrack)
c
         zat = zhigh
         call read_kz(2,zat,1,khighz,iulow,ofebrack)
c
         zlo_ex = zlow - 0.5d0 * ( zhigh - zlow )
         zhi_ex = zhigh + 0.5d0 * ( zhigh - zlow )
         dfsz(2) = 1.d0 / ( zvint(2) - zvint(1) )
         dfsz(1) = dfsz(2)
c
         return
c
c  Else if all available Z-values should be used, handle this case and return:
c
      else if ( numz .eq. nzm ) then

         zlow = 0.0d0
         zmiddle = zat
         zhigh = zavail(nzm)
         zlo_ex = -1.d-6
         zhi_ex = 2.d0 * zavail(nzm) - zavail(nzm-1)
c
         do kz = 1, nzm
            zat = zavail(kz)
            call read_kz(kz,zat,1,khighz,iulow,ofebrack)
            if ( kz .gt. 1 )
     $           dfsz(kz) = 1.d0 / ( zvint(kz) - zvint(kz-1) )
         enddo
c
         dfsz(1) = dfsz(2)
c
         return
c
      endif
c
c  If there should be at least three Z-values:
c						     ! check if Z-range too big
      if ( level_err .gt. 0 ) then
         j = mz
         do while ( j .gt. 2 .and. zhi_t .gt. zval(j-1) - 1.d-6 )
            j = j - 1
         enddo
         i = 1
         do while ( i .lt. j - 1 .and. zlo_t .lt. zval(i+1) - 1.d-6 )
            i = i + 1
         enddo
         if ( nzin .lt. j - i ) then
            write(6,20) nzin, zlo, z, zhi,
     $           'Zlo and Zhi too far apart for given Nzin value'
            stop ' STOP -- READZEXCO: Error: bad Z values. '
         endif
      endif
c						! check for input Z-range
      ilodel = 1
      if ( zlo .lt. -1.d-6 .or.
     $     ( zlo .gt. zhi .and. zhi .ge. -1.d-6 ) .or.
     $     ( zlo .gt. z .and. z .ge. -1.d-6 ) ) then
         if ( z .ge. -1.d-6 ) then
            zlow = z
         else
            zlow = zat
         endif
      else
         zlow = zlo
         ilodel = 0
      endif
      if ( zlow .lt. 1.d-8 ) then
         zlow = 0.d0
         ilo1 = 1
         ilo2 = 1
      else
         zlow = min( zlow , zavail(nzm-2) )
         ilo2 = nzm - 2
         do while ( ilo2 .gt. 2 .and.
     $        zlow .lt. zavail(ilo2-1) + 1.d-6 )
            ilo2 = ilo2 - 1
         enddo
         ilo1 = 1
         do while ( ilo1 .lt. ilo2 .and.
     $        zlow .gt. zavail(ilo1+1) - 1.d-6 )
            ilo1 = ilo1 + 1
         enddo
      endif
c
      ihidel = 1
      if ( zhi .ge. max( z , zlo , -1.d-6 ) ) then
         zhigh = zhi
         ihidel = 0
      else if ( z .gt. -1.d-6 ) then
         zhigh = z
      else
         zhigh = zat
      endif
      zhigh = max( zavail(3) , min( zavail(nzm) , zhigh ) )
      if ( ilodel .eq. 0 .and. ilo1 .ge. nzm + 1 - numz ) then
         ihi2 = nzm
         ilo1 = nzm + 1 - numz
      else
         ihi2 = nzm
         do while ( ihi2 .gt. 3 .and.
     $        zhigh .lt. zavail(ihi2-1) + 1.d-6 )
            ihi2 = ihi2 - 1
         enddo
         ihi1 = 3
         do while ( ihi1 .lt. ihi2 .and.
     $        zhigh .gt. zavail(ihi1+1) - 1.d-6 )
            ihi1 = ihi1 + 1
         enddo
         if ( ihidel .eq. 0 .and. ihi2 .le. numz ) then
            ilo1 = 1
            ihi2 = numz
         endif
      endif
c
c  If the number of Z-values to be used "numz"  is sufficent to encompass the
c  input Z-range, then handle this case and return:
c
      if ( ihi2 - ilo1 .lt. numz ) then
c					! get range -> numz, alternately adding
c							! low and high Z-values
         if ( numz .eq. 3 .and. ihi2 - ilo1 .eq. 1 .and.
     $        zavail(ihi2) - z .lt. z - zavail(ilo1) .and.
     $        ( ihidel .gt. 0 .or. ilodel .eq. 0 ) )
     $        ihi2 = min( ihi2 + 1 , nzm )
         do while ( ihi2 - ilo1 + 1 .lt. numz .and.
     $        ( ( ilo1 .gt. 1 .and. ilodel .gt. 0 ) .or.
     $        ( ihi2 .lt. nzm .and. ihidel .gt. 0 ) ) )
            ilo1 = max( ilo1 - ilodel , 1 )
            if ( ihi2 - ilo1 + 1 .lt. numz )
     $           ihi2 = min( ihi2 + ihidel , nzm )
         enddo
         do while ( ihi2 - ilo1 + 1 .lt. numz )
            ilo1 = max( ilo1 - 1 , 1 )
            if ( ihi2 - ilo1 + 1 .lt. numz )
     $           ihi2 = min( ihi2 + 1 , nzm )
         enddo
c
         zlow = zavail(ilo1)
         if ( numz .eq. 3 ) then
            zmiddle = zavail(ilo1+1)
         else if ( z .lt. -1.d-6 ) then
            zmiddle = ( zhi + zlo ) * 0.5d0
         endif
         zhigh = zavail(ihi2)
         if ( ilo1 .eq. 1 ) then
            zlo_ex = -1.d-6
         else
            zlo_ex = zavail(ilo1-1)
         endif
         zhi_ex = 2.d0 * zavail(ihi2) - zavail(ihi2-1)
c
         do kz = 1, numz
            zat = zavail(ilo1+kz-1)
            call read_kz(kz,zat,1,khighz,iulow,ofebrack)
            if ( kz .gt. 1 )
     $           dfsz(kz) = 1.d0 / ( zvint(kz) - zvint(kz-1) )
         enddo
c
         dfsz(1) = dfsz(2)
c
         return
c
      endif
c
c  The input Z-range does not fit between a set of numz values from zavail();
c  if the input value of numz was 3, then use the input Z-range even if it is
c  VERY wide (check only for unbalanced intervals in log[Z+0.001]), and return:
c
      if ( numz .eq. 3 .and. nzin .eq. 3 ) then
c
         zl_1 = log10( zlow + zdel )
         zl_2 = log10( zmiddle + zdel )
         zl_3 = log10( zhigh + zdel )
         if ( 2.d0 * ( zl_3 - zl_2 ) .lt. zl_2 - zl_1 .or.
     $        2.6d0 * ( zl_2 - zl_1 ) .lt. zl_3 - zl_2 ) then
            zl_2 = ( zl_1 + zl_3 ) * 0.5d0
            zmiddle = 10.d0**zl_2 - zdel
         endif
c
         zat = zlow
         call read_kz(1,zat,1,khighz,iulow,ofebrack)
         zat = zmiddle
         call read_kz(2,zat,1,khighz,iulow,ofebrack)
         zat = zhigh
         call read_kz(3,zat,1,khighz,iulow,ofebrack)
c
         dfsz(3) = 1.d0 / ( zvint(3) - zvint(2) )
         dfsz(2) = 1.d0 / ( zvint(2) - zvint(1) )
         dfsz(1) = dfsz(2)
c
         zlo_ex = min( zlow - zacc(1) ,
     $        10.d0**( zl_1 - min( zl_2 - zl_1 ,
     $        log10( ( zavail(ilo1+1) + zdel )
     $        / ( zavail(ilo1) + zdel ) ) ) ) - zdel )
         if ( zlo_ex .lt. 1.d-8 .and. zlo_ex .gt. 0. ) zlo_ex = 0.d0
         zhi_ex = min( 10.d0**( 2.d0 * zl_3 - zl_2 ) - zdel ,
     $        zhigh + zavail(ihi2) - zavail(ihi2-1) )
c
         return
c
      endif
c
c  The input Z-range does not fit between a set of numz Z-values from zavail(),
c  and the input value of numz was greater than 3; find the corresponding
c  positions in the array coarser array zval() of Z-values available from the
c  'Gz???.x??' files, and check whether this suffices:
c
      if ( zlow .lt. 1.d-8 ) then
         jlo1 = 1
         jlo2 = 1
      else
         jlo2 = mz - 2
         do while ( jlo2 .gt. 2 .and. zlow .lt. zval(jlo2-1) + 1.d-6 )
            jlo2 = jlo2 - 1
         enddo
         jlo1 = 1
         do while ( jlo1 .lt. jlo2 .and.
     $        zlow .gt. zval(jlo1+1) - 1.d-6 )
            jlo1 = jlo1 + 1
         enddo
      endif
      jhi2 = mz
      do while ( jhi2 .gt. 3 .and. zhigh .lt. zval(jhi2-1) + 1.d-6 )
         jhi2 = jhi2 - 1
      enddo
      jhi1 = 3
      do while ( jhi1 .lt. jhi2 .and. zhigh .gt. zval(jhi1+1) - 1.d-6 )
         jhi1 = jhi1 + 1
      enddo
c
      nuse = jhi2 + 1 - jlo1
c
c  If this coarser spacing works, then use it, after shifting the endpoints
c  in to the closest encompassing Z-values from zavail(), and adding as many
c  Z-values from zavail() as is allowed by the value of numz:
c
      if ( nuse .le. numz ) then
c
         do i = jlo1, jhi2
            z_use(i+1-jlo1) = zval(i)
         enddo
         z_use(1) = zavail(ilo1)
         z_use(nuse) = zavail(ihi2)
         if ( ihi2 .eq. nzm .and. ilo1 .eq. nzm - 4 .and. numz .eq. 3 )
     $        z_use(2) = zavail(nzm-2)
c
         k_a = 0
         do while ( nuse .lt. numz .and. k_a .lt. nadd_zavail )
            k_a = k_a + 1
            zat = zavail( iadd_zavail(k_a) )
            if ( zat .gt. z_use(1) + 1.d-6 .and.
     $           zat .lt. z_use(nuse) - 1.d-6 ) then
               i = 2
               do while ( i .lt. nuse .and. z_use(i) + 1.d-6 .lt. zat )
                  i = i + 1
               enddo
               do j = nuse, i, -1
                  z_use(j+1) = z_use(j)
               enddo
               z_use(i) = zat
               nuse = nuse + 1
            endif
         enddo
         if ( nuse .ne. numz ) stop
     $        ' STOP -- READEXCO Error: nuse < numz cannot happen. '
c
c  Else, if the coarser spacing still does not suffice, use an equal-interval
c  spacing in log(Z+0.001):
c
      else
c
c-dont;         z_use(1) = zavail(ilo1)
c-dont;         z_use(numz) = zavail(ihi2)
         z_use(1) = zlow
         z_use(numz) = zhigh
         z_1 = log( z_use(1) + zdel )
         z_2 = log( z_use(numz) + zdel )
         z_3 = ( z_2 - z_1 ) / ( numz - 1 )
         do i = 2, numz - 1
            z_use(i) = exp( ( i - 1 ) * z_3 + z_1 ) - zdel
         enddo
c
      endif
c
c  In either of the above cases, read in the corresponding opacities:
c
      do kz = 1, numz
         call read_kz(kz,z_use(kz),1,khighz,iulow,ofebrack)
         if ( kz .gt. 1 ) dfsz(kz) = 1.d0 / ( zvint(kz) - zvint(kz-1) )
      enddo
c
      dfsz(1) = dfsz(2)
c
      zlow = z_use(1)
      zhigh = z_use(numz)
c
      zlo_ex = min( zlow - zacc(1) ,
     $     10.d0**( zvint(1) - log10( ( zavail(ilo1+1) + zdel )
     $     / ( zavail(ilo1) + zdel ) ) ) - zdel )
      zhi_ex = zhigh + zavail(ihi2) - zavail(ihi2-1)
c							! we are done: return
      return
      end




*
***--------------------------------------------------------------------------

      subroutine read_kz(kz,z,kallrd,khighz,iulow,ofebrack)
*     =====================================================
*
* PARAMETERS to control the offset from zero for the logarithmic interpolation:
*  zdel = 0.001 (must be the same value as in OPAL)
*  xdel = 0.03 (must be the same value as in OPAL)
*  xdelgn93 = 0.005 = xdel value for use with X-interpolation in 'GN93hz' file
*                      among X-tables 0.0, 0.1, 0.2 (to get X = 0.03 mix); note
*                      that 0.005 works slightly better for this than 0.03
*
*
***---------------------------------------------------------------------------

      implicit double precision (a-h, o-z)


      parameter ( zdel=0.001d0, xdel=0.03d0, xdelgn93=0.005d0 )
c
c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=28, nrb=1, nre=28,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      parameter ( mx_hi=2*mx, mo_m1=mo-1 )
c
      common /xhi_opal_z/ xhi_in(mx_hi), xhi_use(mx_hi,nz),
     $     xxx_hi(mx_hi), nx_hi(nz), ireq_hi(mx_hi), khighx(nz),
     $     kavail_xhi, kuse_xhi, kdo_xhi
      save /xhi_opal_z/
c
      parameter ( nrdel=nrb-1, ntdel=ntb-1 )
      parameter ( nrm_m2=nrm-2, nt_m1=nt-1, nre_p1=nre+1, nre_m1=nre-1 )
      parameter ( badlogkval=1.d+35, badlogklim=20.d0 )
c!      parameter ( ks81=ntm-3, ks83=ks81+1, ks60=ks81-21, ks61=ks60+1,
c!     $     alrlo=-8.0d0, flogtlo=3.75d0, flogt60=6.0d0, flogt81=8.1d0 )
      parameter ( ks81=ntm-9, ks83=ks81+1, ks60=ks81-21, ks61=ks60+1,
     $     alrlo=-9.0d0, flogtlo=3.75d0, flogt60=6.0d0, flogt81=8.1d0 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime(mx),cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
c
      character*2 cel_opalmixes(nel_zmix)
      character*8 cfile_opalmixes(n_zmixes)
      common/opalmixes/ xiz_mix(nel_zmix),fninz_mix(nel_zmix),
     $     bracketife_mix(nel_zmix),bracketofe_opalmixes(n_zmixes),
     $     xofe_opalmixes(n_zmixes),xiz_opalmixes(nel_zmix,n_zmixes),
     $     fninz_opalmixes(nel_zmix,n_zmixes),
     $     cel_opalmixes,cfile_opalmixes
      save /opalmixes/
c
      common/b_opal_z/ nta(0:nrm_p1),zz(mx,nz),ntax0(0:nrm),
     $     ntax03(0:nrm), sltlo, slthi, dltlo_inv, dlthi_inv,
     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
      save /b_opal_z/
c
c COMMON /c_opal_ctrl_smooth/ : flags to control the opacity smoothing:
c
      common/c_opal_ctrl_smooth/ init_smo, low_CO_smo, interp_CO_smo
      save /c_opal_ctrl_smooth/
c
      common/alink_opal_z/ NTEMP,NSM,nrlow,nrhigh,RLE,t6arr(100),
     $     coff(100,nrm)
      save /alink_opal_z/
c
      COMMON/CST_OPAL_Z/ NRL,RLS,nset,tmax
      save /CST_OPAL_Z/
c
      common/extopac/ opact,dopact,dopacr,dopactd,fedge,ftredge,fzedge
      save /extopac/
c
      common/recoin_opal_z/ itimeco,mxzero,mx03,kope,igznotgx
      save /recoin_opal_z/
c
      parameter ( mz=8, mz_m1=mz-1, mz_m2=mz-2, mzhi=11, mzal=13,
     $     nzm=mzal+1, nadd_zavail=nzm-mz )
c
      common/zinter_opal_z/ zvalhi(mzhi),nofz(mzhi,5,mo),mnofz(mx),
     $     zval(mz),zalval(mzal),zavail(nzm),iadd_zavail(nadd_zavail)
      save /zinter_opal_z/
c
      character*4 cxfil(5),czfil(mz)
      common/czinte_opal_z/ cxfil,czfil
      save /czinte_opal_z/
c
c  The following common block /opdir/ contains a character file giving the
c  directory where the Gz???.x?? files are to be found.
c
      character*80 copdir
      common/opdir/ copdir
      save /opdir/
c___
c
c  copfil = the full opacity filename (including directory), changed as needed
c  copalt = temporary (needed only for error message: new form of opacity file)
c
      character*80 copfil,copalt
c
c  cin / ch  holds a line read in from the opacity file (with opacity values)
c
c!!      character*137 cin
c!!      character*1 ch(137)
      character*200 cin
      character*1 ch(200)
      equivalence (ch(1),cin)
c
c  lxst is a logical variable used to check whether the new form Gz???.x?? of
c        opacity files exists, rather than the old form Gx??z*
c
      logical lxst
c
c LOCAL ARRAYS:
c  cofzat(4,nrm) = temporary storage for input opacity values as a function of
c                   R, at up to four Z-table values (these will be interpolated
c                   in Z if necessary, prior to being stored)
c  xzalat(3) = X-table value(s) for relevant mix(es) in 'GN93hz'.  For m = mx03
c               (X=0.03), xzalat = { 0.0 , 0.1 , 0.2 }; else, xzalat(1) = xa(m)
c  kzvalnof(4) = Z-indexes in number-of-mixes table nofz for each of the (up
c                 to 4) Z-tables being interpolated among (in Gz???.x?? files)
c  nofmjz(4) = temporary number-of-C-mixes for each of the (up to 4) Z-tables
c               being interpolated among: for a given X-table index m and
c               O-table index j,  nofmjz(i) = nofz(kzvalnof(i),mnofz(m),j)
c  iz_get(4) = flag, for each of the (up to 4) Z-tables being interpolated
c               among, to tell the function mixfind whether to read in or check
c               the Z-composition; 1 = read it in, -1 = check it, 0 = neither
c               (i.e., not presently positioned at the beginning of the file).
c               The function mixfind resets this flag to zero.
c  iz_rew(4) = flag, for each of the (up to 4) Z-tables being interpolated
c               among, to tell the function mixfind whether to rewind the file
c               being read in before looking for the next (Z,X,C,O) mix;
c               1 = rewind, 0 = do not rewind.  The function mixfind resets
c               this flag to zero, unless the mix being looked for is not found
c               or does not follow consecutively after the previous mix, in
c               which case it is reset to unity ("rewind next time")
c  iz_tab(4) = mix-table index in Gz???.x?? opacity files, for each of the (up
c               to 4) Z-tables being interpolated among; it is set by mixfind
c  ico_got(mc,mo) = matrix of flags telling whether the corresponding (C,O)
c                    opacity table was read in/interpolated successfully (for
c                    the current X-table index m); 1 = succeeded, 0 = failed.
c                    In a few cases, near the line C+O = 1-X-Z, a mix may not
c                    be present at enough Z-table values to allow it to be
c                    interpolated in Z, in which case this flag will be set to
c                    zero; rather than extrapolating such a mix in Z, it is
c                    subsequently interpolated in C or O (note that it will
c                    lie near a mix on the line C+O = 1-X-Z, so there will not
c                    be much of an interpolation).
c  nthx0(0:nrm) = temporary version of ntax0(0:nrm), used for "hz"-files
c  nthx03(0:nrm) = temporary version of ntax03(0:nrm), used for "hz"-files
c
      dimension cofzat(4,nrm),xzalat(3),kzvalnof(4),nofmjz(4),
     $     iz_get(4),iz_rew(4),iz_tab(4),ico_got(mc,mo),
     $     nthx0(0:nrm),nthx03(0:nrm)
c
c LOCAL ARRAYS: cofzhi(ntm,nrm,6) is used for temporary storage for opacities
c  from 'GN93hz' or from file with non-zero [O/Fe]; xhi_look(6) and zhi_look(6)
c  are used to hold the 'GN93hz' X-values, including those not found in the
c  'Gz???.x??' files, and the corresponding Z-values to be looked for.
c
      dimension cofzhi(ntm,nrm,6), xhi_look(6), zhi_look(6)
c
c Storage for the compression-flags and opacity file names that are opened:
c
      dimension igzip_sto(0:3)
      character*80 cop_sto(0:3)
c
c-debug-chk[
c-debug;      dimension chk_max(20),chk_min(20),chk_sum(20),chk_ssq(20),
c-debug;     $     n_chk(20)
c-debug-chk]
c-test-xdel[                              ! test xdel values in GN93hz X-interp
c-test-xdel;      parameter ( n_xdtst=4 )
c-test-xdel;      dimension cof_tst(ntm,nrm,n_xdtst),dif_tst(4,n_xdtst),
c-test-xdel;     $     xdel_tst(n_xdtst)
c-test-xdel;      data xdel_tst/0.03,0.01,0.001,0.0001/
c-test-xdel]
c
c===

      if ( kz .le. 0 .or. kz .gt. nz ) stop
     $     ' STOP -- READCO: Error: Z-index out of range. '
c
c If initializations have not yet been performed, or if a new Z-value has been
c specified, then check some parameter values, and do some initializations.
c
      if ( khighz .ge. 0 ) then
         khizat = min( khighz , n_zmixes )
      else if ( ofebrack .eq. 0.0d0 ) then
         khizat = 1
      else
         khizat = 4
      endif
      if ( itimeco .ne. 12345678 .or. kallrd .ne. 0 ) 
     $     call opalinit(khizat,iulow,ofebrack,z,kz)
c							     ! Z out of range?
      if ( z .ge. zval(mz) + 1.d-6 .or. z .le. -1.d-8 ) then
         write(6,8703) z
 8703    format(' '/' STOP -- READCO: Z=',f10.7,
     $        ': Z > 0.1 (or < 0) not allowed!!!')
         stop
      endif
c						      ! accuracy to match table
c!!      zacc(kz) = min( 1.d-6 , max( 1.d-8 , 0.01d0 * z ) )
c!! Taking 1.d-6 for the accuracy 
      zacc(kz) = min( 1.d-3 , max( 1.d-8 , 0.01d0 * z ) )

c
c  Find the correct Z value(s) to be read in, in the range 0.0 to 0.1
c
      kzbelo = mz
      do while( kzbelo .gt. 1 .and. z .le. zval(kzbelo) - 1.d-6 )
         kzbelo = kzbelo-1
      enddo

c
c  If Z is very close to a tabulated Z-value, don't need to interpolate,
c  unless Z --> 0 ; note that nmorez ("number of extra Z-values") indicates
c  the presence and type of interpolation:
c    nmorez = 0 if the tabulated Z-value can be used (with no Z-interpolation),
c    nmorez = 2 if quadratic interpolation among 3 Z-values will be performed,
c    nmorez = 3 if overlapping quadratices will be used among 4 Z-values.
c
      if ( abs( zval(kzbelo) - z ) .le. zacc(kz) ) then
         zat = zval(kzbelo)
         kzlow = kzbelo
         k1zfac = kzbelo
         nmorez = 0
c		     ! else: closest 3 or 4 table Z-values to interpolate among
      else
         zat = z
         kzlow = min( max( 1 , kzbelo - 1 ) , mz_m2 )
         k1zfac = min( kzbelo , mz_m1 )
         nmorez = min( kzbelo + 2 , mz ) - kzlow
      endif
c
      kzlow_m1 = kzlow - 1
      k2zfac = min( k1zfac + 1 , kzlow + nmorez )
c							! find position in nofz
      do i = kzlow, kzlow + nmorez
         kzvalnof(i-kzlow_m1) = int( 100.d0  * zval(i) + 1.01d0 )
      enddo
      kznof = min( int( 100.d0 * z + 1.0001d0 ) , mzhi )
      lznof = max( int( 100.d0 * z + 0.9999d0  ) , 1 )
c
c  Check if need C+O=0.0 "hz" tables: if [O/Fe] = 0 or Z = 0, there will be no
c  need to interpolate in [O/Fe] (set khizat = 1), and if there is no need to
c  interpolate in [O/Fe], then for Z equal to a Z-table value or .01 < Z < .02
c  the "hz" tables will yield no improvement and are not used (set khizat = 0).
c
      if ( ofebrack .eq. 0. .or. z .eq. 0. ) khizat = min( khizat , 1 )
      if ( khizat .eq. 1 .and. ( ( zval(k1zfac) .ge. 0.01d0 .and.
     $     zval(k2zfac) .le. 0.02d0 ) .or. nmorez .eq. 0 ) ) khizat = 0
c
c  If needed, get position in C+O=0.0 "hz"-tables, which have more Z values.
c  Note that nzalmo (for GN93hz) is analogous to nmorez (for Gz???.x?? files)
c
      if ( khizat .gt. 0 ) then
         kzalbe = mzal
         do while( kzalbe .gt. 1 .and. z .le. zalval(kzalbe)-1.d-6 )
            kzalbe = kzalbe-1
         enddo
         if ( abs( zalval(kzalbe) - z ) .le. zacc(kz) ) then
            zat = zalval(kzalbe)
            kzalow = kzalbe
            nzalmo = 0
         else
            kzalow = max( 1 , kzalbe - 1 )
            nzalmo = min( kzalbe + 2 , mzal ) - kzalow
         endif
      endif
c			       ! set the directory-part of the opacity filename
      if ( kope .eq. 0 ) then
         copfil = ' '
      else
         copfil = copdir(1:kope)
      endif
c				! store present m-value; get m-range to read in
      mstore = m
c
      khighx(kz) = 0
      if ( kallrd .le. 0 ) then
         if ( m .le. 0 .or. m .gt. mx ) stop
     $        ' STOP -- READCO: bad index m for X-table to read in. '
         mstart = m
         mend = m
      else
         mstart = 1
         mend = mx
         if ( khighz .gt. 0 .and. mx .eq. 5 .and.
     $        max( abs( xhi_in(1) - xa(1) ) , abs( 0.03d0 - xa(2) ) ,
     $        abs( xhi_in(2) - xa(3) ) , abs( xhi_in(4) - xa(4) ) ,
     $        abs( xhi_in(6) - xa(5) ) ) .le. 1.d-6 ) khighx(kz) = 1
      endif
c			   ! get shifted z-value Z+zdel, for log interpolation
      zzz(kz) = zat + zdel
      zvint(kz) = log10( zzz(kz) )
      zsto(kz) = zat
c			 ! should read in Z-composition from first opacity file
      igetzxi = 1
c
c {--------------------------------------- Begin loop over m values (X-tables):
c
      do mvalue = mstart, mend
c					     ! set X-index m
         m = mvalue
c					     ! later Z-compositions: just check
         if ( m .ne. mstart ) igetzxi = -1
c					     ! get C,O compositions for this m
         xhemx = 1.d0 - xa(m) - zat
         do i = 1, mc
            nc = i
            no = i
            xc(i) = xcs(i)
            xo(i) = xos(i)
c						  ! allow some round-off error:
            if ( xcs(i) .ge. xhemx - 1.d-6 ) then
               xc(i) = xhemx
               xo(i) = xhemx
               goto 3
            endif
         enddo
 3       continue
c				! check that number of C-mixes is correct
c
         if ( nc .ne. nofz(kznof,mnofz(m),1) .and.
     $        nc .ne. nofz(lznof,mnofz(m),1) ) then
            write(6,8704) m,lznof,kznof,nofz(lznof,mnofz(m),1),
     $           nofz(kznof,mnofz(m),1),nc
 8704       format(' '/' STOP -- READCO(m=',i1,
     $           ') bad nc value: nofz({',i2,' or',i3,'},m,1)={',
     $           i1,' or',i2,'} .ne. nc=',i1)
            stop
         endif
c			! check that Z is the same for all m already read in
         zz(m,kz) = zat
         if ( kallrd .eq. 0 ) then
            do i = 1, mx
               if ( i .ne. m .and. itime(i) .eq. 12345678 ) then
                  if ( zz(i,kz) .ne. zat ) then
c-old;                     write(6,8707) m,zz(m,kz),i,zz(i,kz)
c-old; 8707                format(' '/' STOP -- READCO(m=',i1,'): Z(m)=',
c-old;     $                    f10.7,' not equal to Z(',i1,')=',f10.7)
c-old;                     stop
                     itime(i) = 0
                  endif
               endif
            enddo
         endif
c
c........ initialization: for itime, oxf...oxdf, n(m,*,kz)  (xx was done above)
c
         itime(m) = 12345678
c
c........ this is the first time through this m and kz.  Calculate the decadic
c         log of the perimeter points shifted by Z+zdel (to avoid divergence 
c         at origin: zdel=0.001); m refers to xa(m), the hydrogen table value.
c
         do i = 1, nc
c						    ! O,C values for each X,C,O
            oxf(m,i,kz) = log10( zzz(kz) + xo(i) )
            cxf(m,i,kz) = log10( zzz(kz) + xc(i) )
            xcdf(m,i,kz) = xo(no) - xo(i)
            xodf(m,i,kz) = xc(nc) - xc(i)
            cxdf(m,i,kz) = log10( zzz(kz) + xcdf(m,i,kz) )
            oxdf(m,i,kz) = log10( zzz(kz) + xodf(m,i,kz) )
c							   ! present C,O values
            ox(i) = oxf(m,i,kz)
            cx(i) = cxf(m,i,kz)
            xcd(i) = xcdf(m,i,kz)
            xod(i) = xodf(m,i,kz)
            cxd(i) = cxdf(m,i,kz)
            oxd(i) = oxdf(m,i,kz)
         enddo
c				! set and check number-of-mixes table n(m,j,kz)
         do j = 1, nc - 1
            do i = 1, nc
               if ( xcd(j) .ge. xc(i) ) then
                  n(m,j,kz) = i + 1
                  if ( xcd(j) .lt. xc(i) + 1.d-6 ) n(m,j,kz) = i
               endif
            enddo
            if ( n(m,1,kz) .ne. nc .or.
     $           ( n(m,j,kz) .ne. nofz(kznof,mnofz(m),j)
     $           .and. n(m,j,kz) .ne. nofz(lznof,mnofz(m),j) ) ) then
               write(6,8705) m,nc,j,n(m,j,kz),lznof,kznof,j,
     $              nofz(lznof,mnofz(m),j),nofz(kznof,mnofz(m),j)
 8705          format(' '/' STOP -- READCO(m=',i1,',nc=',i1,
     $              ') bad value of n(m,',i1,')=',i1,
     $              ' .ne. nofz({',i2,' or',i3,'},m,',i1,')={',
     $              i1,' or',i2,'}')
               stop
            endif
         enddo
c		      ! nc-th elements sometimes needed, though this may not be
         do j = nc, mc
            n(m,j,kz) = 0
         enddo
c		   		   ! initialize boundaries at low-X,low-R,low-T
         if ( kz .eq. 1 ) then
            if ( m .eq. mxzero ) then
               do i = nrb, nre
                  ntax0(i) = ntb
               enddo
            else if ( m .eq. mx03 ) then
               do i = nrb, nre
                  ntax03(i) = ntb
               enddo
            endif
         endif
c			! end of initializations for this m value
c
c If it will increase accuracy, first read in C+O=0.0 "hz" table(s) from GN93hz
c  (and possibly from a file with [O/Fe] > 0), for the present m value.
c
c Note that 'Gz???.x??' files contain Z=0.05, while 'GN93hz' contains Z=0.04
c  and 0.06; thus, for Z near 0.05, the 'Gz???.x??' opacities are more accurate
c  than the 'GN93hz' opacities.  For .04 < Z < .05 or for .05 < Z < .06 , the
c  'GN93hz' opacities are read in, but their effects are reduced (according to
c  how close Z is to 0.05) by setting the factor  facxhz  to a value less than
c  unity.
c
c Note that X=0.03 (m=2) requires X-interpolation in the 'GN93hz' tables: this
c  is only done if interpolating in [O/Fe] as well, and the effect of the m=2
c  'GN93hz' opacities is nullified by setting  facxhz  to zero.  (The 'GN93hz'
c  opacity shifts for m=2 may be obtained later by interpolating opacity shifts
c  among m=1,3,4 if this is possible.)
c
c					! read 'GN93hz' only if necessary:
         if ( khizat .gt. 0 .and.
     $        ( m .ne. mx03 .or. khizat .gt. 1 ) ) then
c
c				   ! for m = mx03 (X=.03), X-interpolation in
c				   ! GN93hz is less accurate than the Gz???.x??
c				   ! opacities: the GN93hz opacities are needed
c				    ! only for [O/Fe] shifts, so set facxhz=0.
            if ( m .eq. mx03 ) then
               facxhz = 0.d0
               nxdo = 3
               xzalat(1) = 0.d0
               xzalat(2) = 0.1d0
               xzalat(3) = 0.2d0
               do i = nrb, nre
                  nthx03(i) = ntb
               enddo
c				! Else (for m not mx03): do need GN93hz shifts:
            else
               facxhz = 1.d0
c			       ! but near Z = 0.05, the Gz???.x?? opacities are
c						        ! better: reduce facxhz
               if ( zat .gt. 0.04d0 .and. zat .lt. 0.06d0 )
     $              facxhz = 1.d0 - 100.d0 * min( zat - 0.04d0 , 
     $              0.06d0 - zat )
               nxdo = 1
               xzalat(1) = xa(m)
               if ( m .eq. mxzero ) then
                  do i = nrb, nre
                     nthx0(i) = ntb
                  enddo
               endif
            endif
c			      ! number of "excess" X-tables (for interpolation)
            nmorex = nxdo - 1
c			      ! indices for GN93hz Z-tabulation array zalval
            j1 = kzalow
            j2 = j1 + min( 1 , nzalmo )
            j3 = j1 + min( 2 , nzalmo )
            j4 = j1 + min( 3 , nzalmo )
c
            is = 0
            isx = 0
            iu = iulow
c		       ! read C+O=0.0 "hz"-table(s): GN93hz, and if khizat > 1,
c						       ! a file with [O/Fe] > 0
            do iofe = 1, khizat, max( khizat - 1 , 1 )
c								! get filename
               copfil(kope+1:80) = cfile_opalmixes(iofe)
c								! open file
               call open_chk_zip( iu, copfil, igzip,
     $              'STOP -- Error: hz-file (C+O=0.0) not found.' )
c
c				     ! dummy table-number; initial cofzhi index
               itab_dum = 0
               kzsto = 0
c					! get Z-composition(s) from "hz"-files,
               if ( m .eq. mstart ) then
                  igetzxi = 1
c			       ! or just check them if they were already gotten
               else
                  igetzxi = -1
               endif
c					! loop over X values, if more than one
               do kx = 1, nxdo
c					! increment cofzhi mix-store-position
                  kzsto = kzsto + 1
c						  ! loop over file Z values
                  do iz = kzalow, kzalow + nzalmo
c						  ! kat is Z-index in cofzhi
                     kat = kzsto + iz - kzalow
c					       ! find mix; stop if not found
                     i_rewind = 0
                     ifound = mixfind(iu,iofe,igetzxi,i_rewind,
     $                    itab_dum,zalval(iz),xzalat(kx),0.0d0,0.0d0)
                     if ( ifound .eq. 0 ) then
                        write(6,1791) zalval(iz),xzalat(kx),0.0d0,0.0d0,
     $                       cfile_opalmixes(iofe)
 1791                   format(' '/' READCO: Error finding mix Z=',
     $                       f6.4,' X=',f6.4,' C=',f6.4,' O=',f6.4,
     $                       ' file=',a8/' ')
                        stop ' STOP -- READCO: error reading hz-mix. '
                     endif
c								! check [O/Fe]
                     if ( kx .eq. 1 .and. iz .eq. kzalow ) then
                        if ( iofe .eq. 1 ) then
c						! (this cannot happen)
                           if ( bracketofe_opalmixes(1) .ne. 0.d0 ) stop
     $                          ' STOP -- READCO: non-0 [O/Fe]_GN93hz '
c						! is file [O/Fe] too small?
                        else if ( abs(bracketofe_opalmixes(iofe)) .lt.
     $                          min(0.1d0*abs(ofebrack),0.001d0) ) then
                           write(6,2631) bracketofe_opalmixes(iofe),
     $                          cfile_opalmixes(iofe),ofebrack
 2631                      format(' '/' STOP -- READCO: [O/Fe] =',
     $                          f10.6,' in ',a8,
     $                          ' too small relative to',f10.6)
                           stop
                        endif
                     endif
c				  ! loop over logT values, to read in opacities
                     do k = 1, ntm
c					   ! read logT, & logKappa(R) for all R
                        read(iu,7300) cin
c!! 7300                   format(a137)
 7300                   format(a200)                  
                       read(cin,7140) flt, (cofzhi(k,il,kat),il=1,nrm)
c!! 7140                   format(f4.2,19f7.3)
 7140                   format(f4.2,28f7.3)
c								   ! bad logT ?
                        if ( abs(flogtin(k)-flt) .gt. 1.d-5 ) stop
     $                       ' STOP -- READCO: bad logT value. '
c							      ! logKappa(R) is:
                        do il = nrm, 1, -1
c								       ! absent
                           if ( cin(7*il-2:7*il+4) .eq.
     $                          '       ' ) then
                              if ( k .le. max(nta(il),nta(0)) ) then 
                                 stop' STOP -- READCO: bad upper edge. '
                              endif
c
c						  ! get value, for smoothing
c
                              cofzhi(k,il,kat) = 2.d0*cofzhi(k-1,il,kat)
     $                             - cofzhi(k-2,il,kat)
c							     ! should be absent
                           else if ( k .gt. nta(il) .and.
     $                             il .ge. nrb .and. il .le. nre ) then
                              stop ' STOP -- READCO: bad upper edge. '
c									! 9.999
                           else if ( cofzhi(k,il,kat) .gt. 9.d0 ) then
                              if ( ( m .ne. mxzero .and. m .ne. mx03 )
     $                             .or. il .ge. nrm_m2 .or.
     $                             flt .ge. 4.2d0 ) then
                                 stop ' STOP -- READCO: bad low edge. '
                              else if ( m .eq. mxzero ) then
                                 nthx0(il) = max( nthx0(il) , k + 1 )
                              else if ( m .eq. mx03 ) then
                                 nthx03(il) = max( nthx03(il) , k + 1 )
                              endif
c						  ! get value, for smoothing
c
                              cofzhi(k,il,kat) = 2.d0*cofzhi(k,il+1,kat)
     $                             - cofzhi(k,il+2,kat)
                           endif
c						! end of check-logKappa(R) loop
                        enddo
c						! end of T-loop
                     enddo
c						! end of Z-loop
                  enddo
c						! interpolate in Z, if needed
                  if ( nzalmo .gt. 0 ) then
                     kdelzhi = kzalow - kzsto
                     do k = 1, ntm
                        do il = 1, nrm
                           cofzhi(k,il,kzsto) = qzinter(is,1,zat,
     $                          nzalmo,cofzhi(k,il,j1-kdelzhi),
     $                          cofzhi(k,il,j2-kdelzhi),
     $                          cofzhi(k,il,j3-kdelzhi),
     $                          cofzhi(k,il,j4-kdelzhi),
     $                          zalval(j1),zalval(j2),zalval(j3),
     $                          zalval(j4),zdel)
                           is = 1
                        enddo
                     enddo
                  endif
c						! end of X-loop
               enddo
c			 ! close hz-file (have all necessary opacities from it)
c
               call close_chk_zip( iu, copfil, igzip )
c
c						! interpolate in X if necessary
               if ( nxdo .eq. 3 ) then
                  do k = 1,ntm
                     do il = 1,nrm
c-test-xdel[
c-test-xdel;                        do ij = 1,n_xdtst
c-test-xdel;                           cof_tst(k,il,ij) = qzinter(isx,ij+2,
c-test-xdel;     $                          xa(m),nmorex,cofzhi(k,il,1),
c-test-xdel;     $                          cofzhi(k,il,2),cofzhi(k,il,3),
c-test-xdel;     $                          0.0,xzalat(1),xzalat(2),xzalat(3),
c-test-xdel;     $                          0.0,xdel_tst(ij))
c-test-xdel;                        enddo
c-test-xdel]
                        cofzhi(k,il,1) = qzinter(isx,2,xa(m),nmorex,
     $                       cofzhi(k,il,1),cofzhi(k,il,2),
     $                       cofzhi(k,il,3),0.0d0,xzalat(1),xzalat(2),
     $                       xzalat(3),0.0d0,xdelgn93)
                        isx = 1
                     enddo
                  enddo
               endif
c-test-xdel[
c-test-xdel;               do ij = 1,n_xdtst
c-test-xdel;                  do il = 1,nrm
c-test-xdel;                     do k = 1,ntm
c-test-xdel;                        coff(k,il) = cof_tst(k,il,ij)
c-test-xdel;                     enddo
c-test-xdel;                  enddo
c-test-xdel;                  if ( init_smo .gt. 0 ) then
c-test-xdel;                     tmax = 10.
c-test-xdel;                     nset = ks81
c-test-xdel;                     NSM = 1
c-test-xdel;                     RLS = alrf(1)
c-test-xdel;                     RLE = alrf(nrm)
c-test-xdel;                     nrhigh = int(dfsr(nr)*(RLE-RLS)+1.00001)
c-test-xdel;                     nrlow = 1
c-test-xdel;                     call opaltab
c-test-xdel;                  endif
c-test-xdel;                  do il = 1,nrm
c-test-xdel;                     do k = 1,ntm
c-test-xdel;                        cof_tst(k,il,ij) = coff(k,il)
c-test-xdel;                     enddo
c-test-xdel;                  enddo
c-test-xdel;               enddo
c-test-xdel]
c			    ! transfer opacities from Z,X-interpolation storage
               do il = 1,nrm
                  do k = 1,ntm
                     coff(k,il) = cofzhi(k,il,1)
                  enddo
               enddo
c					 ! smooth hz-opacities, if init_smo > 0
               if ( init_smo .gt. 0 ) then
                  tmax = 10.d0
                  nset = ks81
                  NSM = 1
c					 ! note: MUST have all dfsr(i) = 1./0.5
                  RLS = alrf(1)
                  RLE = alrf(nrm)
                  nrhigh = int( dfsr(nr) * ( RLE - RLS ) + 1.00001d0 )
                  nrlow = 1
c				  ! fit and smooth OPAL kappas, up to T6 = tmax
                  call opaltab
c				! end of hz-opacity smoothing
               endif
c					    ! set any missing values to 1.0E+35
               do il = nre, nrb, -1
                  if ( nta(il) .lt. ntm ) then
                     do k = nta(il) + 1, ntm
                        coff(k,il) = badlogkval
                     enddo
                  endif
                  if ( m .eq. mxzero ) then
                     if ( il .lt. nre )
     $                    nthx0(il) = max( nthx0(il) , nthx0(il+1) )
                     if ( nthx0(il) .gt. ntb ) then
                        do k = ntb, nthx0(il) - 1
                           coff(k,il) = badlogkval
                        enddo
                     endif
                  else if ( m .eq. mx03 ) then
                     if ( il .lt. nre )
     $                    nthx03(il) = max( nthx03(il) , nthx03(il+1) )
                     if ( nthx03(il) .gt. ntb ) then
                        do k = ntb, nthx03(il) - 1
                           coff(k,il) = badlogkval
                        enddo
                     endif
                  endif
               enddo
c					     ! IF have all needed hz-opacities,
c					     !  then store them in matrix CO
               if ( iofe .eq. khizat ) then
                  moat = mo
c					     ! if have two hz-opacity sets,
c					     !  store the results of first set,
                  if ( iofe .gt. 1 ) then
                     do il = 1,nr
                        do k = 1,nt
                           co(m,mc,mo,k,il,kz) = cof(k,il)
                        enddo
                     enddo
c				  ! and update storage-position for second set
                     moat = mo_m1
                  endif
c					     ! store present hz-opacity set
                  do il = 1,nr
                     jl = il+nrdel
                     do k = 1,nt
                        co(m,mc,moat,k,il,kz) = coff(k+ntdel,jl)
                     enddo
                  enddo
c				! else, if have only first of 2 hz-files:
               else
c				  ! store it where it will not be overwritten
                  do il = 1,nr
                     jl = il+nrdel
                     do k = 1,nt
                        cof(k,il) = coff(k+ntdel,jl)
                     enddo
                  enddo
               endif
c				! end of loop reading in C+O=0.0 table(s)
            enddo
c			    ! only check (don't store) Z-composition, in future
            igetzxi = -1
c				! end of obtaining C+O=0.0 'GN93hz'-table(s)
         endif
c
c Read in opacities from 'Gz???.x??' files, interpolating in Z if necessary,
c  for files with the present m value (i.e., with the present X value xa(m)).
c UPDATE: 25 May 1999: look for newer format  Gz???.x??  first; if this is not
c  found, then look for older format  Gx??z*  instead.
c
c			! turn off all "mix-acquired" flags
         do i = 1, mc
            do j = 1, mo
               ico_got(i,j) = 0
            enddo
         enddo
c					! get filename(s) and open file(s):
         do iu = iulow,iulow+nmorez
            if ( igznotgx .ge. 0 ) then
               copfil(kope+1:kope+1) = 'G'
               copfil(kope+2:kope+5) = czfil(kzlow+iu-iulow)
               if ( copfil(kope+4:kope+4) .eq. ' ' )
     $              copfil(kope+4:kope+4) = '0'
               if ( copfil(kope+5:kope+5) .eq. ' ' )
     $              copfil(kope+5:kope+5) = '0'
               copfil(kope+6:80) = cxfil(m)
               copfil(kope+6:kope+6) = '.'
            else
               copfil(kope+1:kope+4) = cxfil(m)
               copfil(kope+5:80) = czfil(kzlow+iu-iulow)
            endif
            if ( igznotgx .eq. 0 ) then
               call inqfil(copfil,lxst)
               if ( .not. lxst )
     $               call inqfil( copfil(1:kope+9) // '.gz' , lxst )
               if ( .not. lxst )
     $              call inqfil( copfil(1:kope+9) // '.Z' , lxst )
               if ( lxst ) then
                  igznotgx = 1
               else
                  copalt = copfil
                  copfil(kope+1:kope+4) = cxfil(m)
                  copfil(kope+5:80) = czfil(kzlow+iu-iulow)
                  call inqfil(copfil,lxst)
                  if ( .not. lxst ) then
                     k_e = kope + 8
                     if ( copfil(k_e:k_e) .eq. ' ' ) k_e = k_e - 1
                     if ( copfil(k_e:k_e) .eq. ' ' ) k_e = k_e - 1
                     call inqfil( copfil(1:k_e) // '.gz' , lxst )
                     if ( .not. lxst )
     $                    call inqfil( copfil(1:k_e) // '.Z' , lxst )
                  endif
                  if ( lxst ) then
                     igznotgx = -1
                  else
                     write(6,7399) copalt,copfil
 7399                format(' '/' STOP -- READCO: neither Gz???.x??',
     $                    ' nor Gx??z* OPAL opacity files found:'/
     $                    ' ',a80/' ',a80)
                     stop
                  endif
               endif
            endif
            cop_sto(iu-iulow) = copfil
            call open_chk_zip( iu, copfil, igzip_sto(iu-iulow),
     $           'STOP -- Error: opacity file not found.' )
         enddo
c				! read in Z-composition only for 1st file
         do i = 1, nmorez + 1
            iz_get(i) = -1
            iz_rew(i) = 0
            iz_tab(i) = 0
         enddo
         iz_get(1) = igetzxi
c				! Z-position indices j1 to j4, for array zval
         j1 = kzlow
         j2 = j1 + min( 1 , nmorez )
         j3 = j1 + min( 2 , nmorez )
         j4 = j1 + min( 3 , nmorez )
         is = 0
c				! loop over dXo (excess oxygen) index j
         do j = 1, no - 1
c				! number of dXc values at Z for this dXo
            ihi = n(m,j,kz)
            do jz = 1, nmorez + 1
               nofmjz(jz) = nofz(kzvalnof(jz),mnofz(m),j)
            enddo
c				! loop over dXc (excess carbon) index i
            do i = 1, ihi
c				! number of "extra" Z-values, for interpolation
               nmorat = nmorez
c				! note: the case i=1,j=1 will ALWAYS have is1=1
               is1 = 1
c					     ! loop over Z values: find tables
               do iz = kzlow, kzlow + nmorez
c					     ! other Z-index jz starts at 1
                  jz = iz - kzlow_m1
                  iu = iulow + iz - kzlow
c					  ! if a mix (with higher Z-value) is
c					  ! missing a needed C-O table, then
c					  ! rewinding may work (if needed table
c					  ! duplicates an earlier C-O table)
c
                  if ( i .gt. nofmjz(jz) ) iz_rew(jz) = 1
                  cget = xcs(i)
                  if ( i .eq. ihi )
     $                 cget = min( cget , 1.d0-xa(m)-xos(j)-zval(iz) )
                  ifound = mixfind(iu,1,iz_get(jz),iz_rew(jz),
     $                 iz_tab(jz),zval(iz),xa(m),cget,xos(j))
c							      ! if table is not
c					      ! present in this file, it cannot
c					      ! be used in the Z-interpolation,
                  if ( ifound .eq. 0 ) then
c							     ! so reduce nmorat
                     nmorat = min( iz - kzlow - 1 , nmorat )
c							     ! (cannot happen):
                     if ( nmorat .lt. 0 ) then
                        write(6,1791) zval(iz),xa(m),cget,xos(j)
                        stop ' STOP -- READCO: error reading mix. '
c-debug-chk[
c-debug-chk;                     else
c-debug-chk;                        write(6,1878) m,z,i,j,jz,nofmjz(jz),
c-debug-chk;     $                       iz_rew(jz),nmorat
c-debug-chk; 1878                   format(' m=',i1,' Z=',f9.7,
c-debug-chk;     $                       ' cannot find i=',i1,' j=',i1,
c-debug-chk;     $                       ' Z(jz=',i1,'): Nofmjz=',i1,
c-debug-chk;     $                       ' irew=',i2,' nmorat=',i2)
c-debug-chk]
                     endif
                     is1 = 2
                     is = 0
                  endif
               enddo
c		       ! if needed table exists at enough Z-values, read it in:
c
               if ( nmorat .eq. nmorez .or. nmorat .eq. 2 ) then
c								 ! loop over T:
                  do k = 1, ntm
c					        ! loop over Z values: read line
                     do iz = kzlow, kzlow + nmorat
                        jz = iz - kzlow_m1
                        iu = iulow + iz - kzlow
c					   ! read logT, & logKappa(R) for all R
                        read(iu,7300) cin
                        read(cin,7140) flt,(cofzat(jz,il),il=1,nrm)
c								   ! bad logT ?
                        if ( abs(flogtin(k)-flt) .gt. 1.d-5 ) then 
                           stop ' STOP -- READCO: bad logT value. '
                        endif
c								  ! store logT
                        if ( k .ge. ntb ) alt(k-ntdel) = flt - 6.d0
                        flogtin(k) = flt
c							      ! logKappa(R) is:
                        do il = nrm, 1, -1
c								       ! absent
                           if ( cin(7*il-2:7*il+4) .eq.
     $                          '       ' ) then
                              if ( k .le. max(nta(il),nta(0)) ) stop
     $                             ' STOP -- READCO: bad upper edge. '
                              cofzat(1,il) = badlogkval
c								 ! or should be
                           else if ( k .gt. nta(il) .and.
     $                             il .ge. nrb .and. il .le. nre ) then
                              stop ' STOP -- READCO: bad upper edge. '
c								     ! or 9.999
                           else if ( cofzat(jz,il) .gt. 9.d0 ) then
                              if ( ( m .ne. mxzero .and. m .ne. mx03 )
     $                             .or. il .ge. nrm_m2 .or.
     $                             flt .ge. 4.2d0 ) then
                                 stop ' STOP -- READCO: bad low edge. '
c
c							    ! set lower bounds:
                              else if ( m .eq. mxzero ) then
                                 ntax0(il) = max( ntax0(il) , k + 1 )
                              else if ( m .eq. mx03 ) then
                                 ntax03(il) = max( ntax03(il) , k + 1 )
                              endif
c							       ! for smoothing:
                              cofzat(jz,il) = 2.d0*cofzat(jz,il+1)
     $                             -cofzat(jz,il+2)
                           endif
c						! end of check-logKappa(R) loop
                        enddo
c					! end of Z-loop
                     enddo
c				  ! interpolate logKappa(R) in Z; store in COFF
                     do il = 1, nrm
c				    ! if opacity missing, extrapolate value for
c								    ! smoothing
                        if ( abs(cofzat(1,il)) .gt. badlogklim ) then
                           coff(k,il) = 2.d0*coff(k-1,il) - coff(k-2,il)
c
c						      ! else if table-Z is O.K.
                        else if ( nmorez .eq. 0 ) then
                           coff(k,il) = cofzat(1,il)
c							! else, Z-interpolation
                        else
                           coff(k,il) = qzinter(is,is1,zat,nmorat,
     $                          cofzat(j1-kzlow_m1,il),
     $                          cofzat(j2-kzlow_m1,il),
     $                          cofzat(j3-kzlow_m1,il),
     $                          cofzat(j4-kzlow_m1,il),zval(j1),
     $                          zval(j2),zval(j3),zval(j4),zdel)
                           is = 1
                        endif
                     enddo
c						! end of T-loop
                  enddo
c					    ! smooth opacities, if init_smo > 0
                  if ( init_smo .gt. 0 ) then
                     tmax = 10.d0
                     nset = ks81
                     NSM = 1
c					 ! note: MUST have all dfsr(i) = 1./0.5
                     RLS = alrf(1)
                     RLE = alrf(nrm)
                     nrhigh = int( dfsr(nr) * ( RLE - RLS ) + 1.00001d0)
                     nrlow = 1
c				  ! fit and smooth OPAL kappas, up to T6 = tmax
                     call opaltab
c				  ! end of opacity smoothing
                  endif
c				  ! store opacities
                  do il = 1, nr
                     jl = il + nrdel
                     do k = 1, nt
                        co(m,i,j,k,il,kz) = coff(k+ntdel,jl)
                     enddo
                  enddo
c				   ! set flag indicating this table was read in
                  ico_got(i,j) = 1
c				   ! end of reading in table
               endif
c			! end of loop over dXc (excess carbon) index
            enddo
c			! end of loop over dXo (excess oxygen) index
         enddo
c
c........ Read remaining diagonal tables (along line Y=0 in dXc,dXo plane)
c
         do jz = 1, nmorez + 1
            nofmjz(jz) = nofz(kzvalnof(jz),mnofz(m),1) - 1
            if ( nofmjz(jz) .lt. no - 1 ) iz_rew(jz) = 1
         enddo
c			! loop over dXc (excess carbon) inverted-index j; note
c			  ! that table being read in will be stored at i=(no-j)
         do j = 1, no - 1
            nmorat = nmorez
            is1 = 1
            nomj = no - j
c					     ! loop over Z values: find tables
            do iz = kzlow, kzlow + nmorat
               iu = iulow + iz - kzlow
               jz = iz - kzlow_m1
               oget = 1.d0 - xa(m) - xcs(nomj) - zval(iz)
               ifound = mixfind(iu,1,iz_get(jz),iz_rew(jz),
     $              iz_tab(jz),zval(iz),xa(m),xcs(nomj),oget)
               if ( ifound .eq. 0 ) then
                  nmorat = min( iz - kzlow - 1 , nmorat )
                  if ( nmorat .lt. 0 ) then
                     write(6,1791) zval(iz),xa(m),xcs(nomj),oget
                     stop ' STOP -- READCO: error reading mix. '
c-debug-chk[
c-debug-chk;                  else
c-debug-chk;                     write(6,2878) m,z,j,mo,jz,nofmjz(jz),
c-debug-chk;     $                    iz_rew(jz),nmorat
c-debug-chk; 2878                format(' m=',i1,' Z=',f9.7,
c-debug-chk;     $                    ' cannot find i=',i1,' j=',i1,
c-debug-chk;     $                    ' Z(jz=',i1,'): Nofmjz=',i1,
c-debug-chk;     $                    ' irew=',i2,' nmorat=',i2)
c-debug-chk]
                  endif
                  is1 = 2
                  is = 0
               endif
            enddo
c		       ! if needed table exists at enough z-values, read it in:
c
            if ( nmorat .eq. nmorez .or. nmorat .eq. 2 ) then
c								 ! loop over T:
               do k = 1, ntm
c					        ! loop over Z values: read line
                  do iz = kzlow, kzlow + nmorat
                     jz = iz - kzlow_m1
                     iu = iulow + iz - kzlow
c					   ! read logT, & logKappa(R) for all R
                     read(iu,7300) cin
                     read(cin,7140) flt,(cofzat(jz,il),il=1,nrm)
c								   ! bad logT ?
                     if ( abs(flogtin(k)-flt) .gt. 1.d-5 ) stop
     $                    ' STOP -- READCO: bad logT value. '
c							      ! logKappa(R) is:
                     do il = nrm, 1, -1
c								       ! absent
                        if ( cin(7*il-2:7*il+4) .eq. '       ' ) then
                           if ( k .le. max(nta(il),nta(0)) ) stop
     $                          ' STOP -- READCO: bad upper edge. '
                           cofzat(1,il) = badlogkval
c								 ! or should be
                        else if ( k .gt. nta(il) .and.
     $                          il .ge. nrb .and. il .le. nre ) then
                           stop ' STOP -- READCO: bad upper edge. '
c								     ! or 9.999
                        else if ( cofzat(jz,il) .gt. 9.d0 ) then
                           if ( ( m .ne. mxzero .and. m .ne. mx03 )
     $                          .or. il .ge. nrm_m2 .or.
     $                          flt .ge. 4.2d0 ) then
                              stop ' STOP -- READCO: bad low edge. '
                           else if ( m .eq. mxzero ) then
                              ntax0(il) = max( ntax0(il) , k + 1 )
                           else if ( m .eq. mx03 ) then
                              ntax03(il) = max( ntax03(il) , k + 1 )
                           endif
c								! for smoothing
                           cofzat(jz,il) = 2.d0*cofzat(jz,il+1)
     $                          - cofzat(jz,il+2)
                        endif
c						! end of check-logKappa(R) loop
                     enddo
c					! end of Z-loop
                  enddo
c				   ! interpolate in Z; store in COFF
                  do il = 1, nrm
                     if ( abs(cofzat(1,il)) .gt. badlogklim ) then
                        coff(k,il) = 2.d0*coff(k-1,il) - coff(k-2,il)
                     else if ( nmorez .eq. 0 ) then
                        coff(k,il) = cofzat(1,il)
                     else
                        coff(k,il) = qzinter(is,is1,zat,nmorat,
     $                       cofzat(j1-kzlow_m1,il),
     $                       cofzat(j2-kzlow_m1,il),
     $                       cofzat(j3-kzlow_m1,il),
     $                       cofzat(j4-kzlow_m1,il),zval(j1),
     $                       zval(j2),zval(j3),zval(j4),zdel)
                        is = 1
                     endif
                  enddo
c						! end of T-loop
               enddo
c					    ! smooth opacities, if init_smo > 0
               if ( init_smo .gt. 0 ) then
                  tmax = 10.d0
                  nset = ks81
                  NSM = 1
c					 ! note: MUST have all dfsr(i) = 1./0.5
                  RLS = alrf(1)
                  RLE = alrf(nrm)
                  nrhigh = int( dfsr(nr) * ( RLE - RLS ) + 1.00001d0 )
                  nrlow = 1
c				  ! fit and smooth OPAL kappas, up to T6 = tmax
                  call opaltab
c				  ! end of opacity smoothing
               endif
c				  ! store opacities
               do il = 1,nr
                  jl = il+nrdel
                  do k = 1,nt
                     co(m,nomj,mo,k,il,kz) = coff(k+ntdel,jl)
                  enddo
               enddo
c				   ! set flag indicating this table was read in
               ico_got(nomj,mo) = 1
c				    ! end of reading in table
            endif
c			! end of loop over dXc (excess carbon) inverted-index
         enddo
c					! close 'Gz???.x??' files
         do iu = iulow, iulow + nmorez
            call close_chk_zip( iu, cop_sto(iu-iulow),
     $           igzip_sto(iu-iulow) )
         enddo
c				    ! for X=0 or .03, ensure low-R,low-T corner
c				    !  has no steps in the wrong direction
         if ( m .eq. mxzero ) then
            do il = nre_m1, nrb, -1
               ntax0(il) = max( ntax0(il) , ntax0(il+1) )
            enddo
         else if ( m .eq. mx03 ) then
            do il = nre_m1,nrb,-1
               ntax03(il) = max( ntax03(il) , ntax03(il+1) )
            enddo
         endif
c			! Set any missing opacity values (high-T,R or low-T,R,X
c			!  corners) to badlogkval = 1.0E+35
         do il = 1, nr
            jl = il + nrdel
            if ( m .eq. mxzero ) then
               khi = ntax0(jl) - ntb
            else if ( m .eq. mx03 ) then
               khi = ntax03(jl) - ntb
            else
               khi = 0
            endif
            if ( khi .gt. 0 ) then
               do j = 1,mo
                  if ( j .lt. no ) then
                     ihi = n(m,j,kz)
                  else if ( j .eq. mo ) then
                     ihi = no - 1
                  else
                     ihi = 0
                  endif
                  if ( ihi .gt. 0 ) then
                     do k = 1, khi
                        do i = 1, ihi
                           co(m,i,j,k,il,kz) = badlogkval
                        enddo
                     enddo
                  endif
               enddo
            endif
            if ( nta(jl) .lt. ntm ) then
               do j = 1, mo
                  if ( j .lt. no ) then
                     ihi = n(m,j,kz)
                  else if ( j .eq. mo ) then
                     ihi = no - 1
                  else
                     ihi = 0
                  endif
                  if ( ihi .gt. 0 ) then
                     do k = nta(jl) + 1 - ntdel, nt
                        do i = 1, ihi
                           co(m,i,j,k,il,kz) = badlogkval
                        enddo
                     enddo
                  endif
               enddo
            endif
         enddo
c
c........ Interpolate any missing opacity tables in dXc or dXo; these will be
c          at or near the line Y=0, arising from the fact that files for higher
c          Z-values may not have had as many dXc or dXo tables in them as are
c          needed for the input (interpolated) Z-value.
c							  ! first, main tables:
         do j = 1, no - 1
            ihi = n(m,j,kz)
            do i = 1, ihi
c					      ! if flag indicates missing table
               if ( ico_got(i,j) .eq. 0 ) then
                  oat = xos(j)
                  cat = min( xcs(i) , 1.d0 - xa(m) - zat - oat )
c-debug-chk[
c-debug-chk;                  write(6,1973) zat,m,i,j,ihi,cat,oat
c-debug-chk; 1973             format(' '/'     Z=',f9.7,' --- interpolate',
c-debug-chk;     $                 ' mix (m=',i1,',i=',i1,',j=',i1,
c-debug-chk;     $                 ') with ihi=',i1,' C=',f10.7,' O=',f10.7)
c-debug-chk;                  difmax = -9.999999
c-debug-chk;                  difmin = 9.999999
c-debug-chk;                  sumdif = 0.
c-debug-chk;                  sumsq = 0.
c-debug-chk;                  numsq = 0
c-debug-chk;                  diflmax = -9.999999
c-debug-chk;                  diflmin = 9.999999
c-debug-chk;                  sumldif = 0.
c-debug-chk;                  sumlsq = 0.
c-debug-chk]
c					   ! if C > or = O in missing table,
                  if ( cat .ge. oat ) then
c					   ! then the only C-value that can be
c					   ! missing is the next-to-highest one
c						       ! at ihi-1 = n(m,j,kz)-1
c
                     if ( ico_got(ihi,j) .le. 0 .and. i .ne. ihi )
     $                    write(6,1873) zat,m,ihi,j,ihi,
     $                    min(xcs(ihi),1.d0-xa(m)-zat-oat),oat
                     if ( i .ne. ihi-1 .or. i .lt. 3 )
     $                    write(6,1873) zat,m,i,j,ihi,
     $                    min(xcs(i),1.d0-xa(m)-zat-oat),oat
 1873                format(' '/'     Z=',f9.7,' ??? CANNOT miss',
     $                    ' mix (m=',i1,',i=',i1,',j=',i1,
     $                    ') with ihi=',i1,' C=',f10.7,' O=',f10.7)
                     if ( ico_got(ihi,j) .le. 0 .or.
     $                    i .lt. 3 .or. i .ne. ihi - 1 ) stop
     $                    ' STOP -- READCO: mix CANNOT be missing. '
c
                     im2 = i - 2
                     im1 = i - 1
                     cxhi = log10( zzz(kz) + min( xcs(ihi) ,
     $                    1.d0 - xa(m) - zat - oat ) )
c-debug-chk[
c-debug-chk;                     write(6,1974) i,j,'C',xcs(i),'C',cx(i),
c-debug-chk;     $                    im2,j,im1,j,ihi,j,'C',xcs(im2),xcs(im1),
c-debug-chk;     $                    min(xcs(ihi),1.-xa(m)-zat-oat),
c-debug-chk;     $                    'C',cx(im2),cx(im1),cxhi
c-debug-chk; 1974                format('      --- interpolate (',i1,',',i1,
c-debug-chk;     $                    ') in ',a1,'=',f10.6,' log',a1,'=',f10.6,
c-debug-chk;     $                    ' among (',i1,',',i1,') (',i1,',',i1,
c-debug-chk;     $                    ') (',i1,',',i1,'): ',a1,'=',3f10.6,
c-debug-chk;     $                    ' log',a1,'=',3f10.6)
c-debug-chk]
c				! interpolate in C to get missing table:
                     is = 0
                     do il = 1, nr
                        do k = 1, nt
                           if ( abs( co(m,im1,j,k,il,kz) ) .lt.
     $                          badlogklim ) then
                              co(m,i,j,k,il,kz) = quad(is,1,cx(i),
     $                             co(m,im2,j,k,il,kz),
     $                             co(m,im1,j,k,il,kz),
     $                             co(m,ihi,j,k,il,kz),
     $                             cx(im2),cx(im1),cxhi)
                              is = 1
c-debug-chk[
c-debug-chk;                              if ( t6list(k) .gt. 0.09999 ) then
c-debug-chk;                                 dif = co(m,i,j,k,il,kz)
c-debug-chk;     $                                -co(m,ihi,j,k,il,kz)
c-debug-chk;                                 difl = co(m,i,j,k,il,kz)
c-debug-chk;     $                                -co(m,im1,j,k,il,kz)
c-debug-chk;                                 difmax = max(dif,difmax)
c-debug-chk;                                 difmin = min(dif,difmin)
c-debug-chk;                                 sumdif = sumdif+dif
c-debug-chk;                                 sumsq = sumsq+dif**2
c-debug-chk;                                 numsq = numsq+1
c-debug-chk;                                 diflmax = max(difl,diflmax)
c-debug-chk;                                 diflmin = min(difl,diflmin)
c-debug-chk;                                 sumldif = sumldif+difl
c-debug-chk;                                 sumlsq = sumlsq+difl**2
c-debug-chk;                              endif
c-debug-chk]
                           endif
                        enddo
                     enddo
                     ico_got(i,j) = -1
c					! else, if C < O in missing table, but
c					       ! it is not on the diagonal Y=0:
                  else if ( i .lt. ihi ) then
c					      ! then the only O-value that can
c					      ! be missing is next-to-highest
c						      ! one, at j = n(m,i,kz)-1
c
                     if ( ico_got(i,mo) .le. 0 )
     $                    write(6,2873) z,m,i,mo,n(m,1,kz)-1,
     $                    xcs(i),1.d0-xa(m)-zat-xcs(i)
                     if ( j .lt. 3 .or. j .ne. n(m,i,kz)-1 )
     $                    write(6,2873) z,m,i,j,ihi,
     $                    min(xcs(i),1.d0-xa(m)-zat-oat),oat,
     $                    ' n(m,i)=',n(m,i,kz)
 2873                format(' '/'     Z=',f9.7,' ??? CANNOT miss',
     $                    ' mix (m=',i1,',i=',i1,',j=',i1,
     $                    ') with ihi=',i1,' C=',f10.7,' O=',f10.7,
     $                    a8,i1)
                     if ( ico_got(i,mo) .le. 0 .or.
     $                    j .lt. 3 .or. j .ne. n(m,i,kz)-1 ) stop
     $                    ' STOP -- READCO: mix CANNOT be missing. '
c
                     jm2 = j - 2
                     jm1 = j - 1
                     oxhi = log10( 1.d0 - xa(m) - zat - xcs(i)
     $                    + zzz(kz) )
c-debug-chk[
c-debug-chk;                     write(6,2974) i,j,'O',xos(j),'O',ox(j),
c-debug-chk;     $                    i,jm2,i,jm1,i,mo,'O',xos(jm2),xos(jm1),
c-debug-chk;     $                    1.-xa(m)-zat-xcs(i),
c-debug-chk;     $                    'O',ox(jm2),ox(jm1),oxhi
c-debug-chk; 2974                format('      --- interpolate (',i1,',',i1,
c-debug-chk;     $                    ') in ',a1,'=',f10.6,' log',a1,'=',f10.6,
c-debug-chk;     $                    ' among (',i1,',',i1,') (',i1,',',i1,
c-debug-chk;     $                    ') (',i1,',',i1,'): ',a1,'=',3f10.6,
c-debug-chk;     $                    ' log',a1,'=',3f10.6)
c-debug-chk]
c				! interpolate in O to get missing table:
                     is = 0
                     do il = 1, nr
                        do k = 1, nt
                           if ( abs( co(m,i,jm1,k,il,kz) ) .lt.
     $                          badlogklim ) then
                              co(m,i,j,k,il,kz) = quad(is,1,ox(j),
     $                             co(m,i,jm2,k,il,kz),
     $                             co(m,i,jm1,k,il,kz),
     $                             co(m,i,mo,k,il,kz),
     $                             ox(jm2),ox(jm1),oxhi)
                              is = 1
c-debug-chk[
c-debug-chk;                              if ( t6list(k) .gt. 0.09999 ) then
c-debug-chk;                                 dif = co(m,i,j,k,il,kz)
c-debug-chk;     $                                -co(m,i,mo,k,il,kz)
c-debug-chk;                                 difl = co(m,i,j,k,il,kz)
c-debug-chk;     $                                -co(m,i,jm1,k,il,kz)
c-debug-chk;                                 difmax = max(dif,difmax)
c-debug-chk;                                 difmin = min(dif,difmin)
c-debug-chk;                                 sumdif = sumdif+dif
c-debug-chk;                                 sumsq = sumsq+dif**2
c-debug-chk;                                 numsq = numsq+1
c-debug-chk;                                 diflmax = max(difl,diflmax)
c-debug-chk;                                 diflmin = min(difl,diflmin)
c-debug-chk;                                 sumldif = sumldif+difl
c-debug-chk;                                 sumlsq = sumlsq+difl**2
c-debug-chk;                              endif
c-debug-chk]
                           endif
                        enddo
                     enddo
                     ico_got(i,j) = -1
c					! else, if C < O in missing table, and
c					!  it is on the diagonal Y=0 (note that
c					!  this never actually happens):
                  else
                     nmorat = 3
                     j3 = 3
                     do while( j3 .lt. no - 1 .and. cat .gt. xcs(j3)
     $                    .and. ico_got(j3+1,mo) .gt. 0 )
                        j3 = j3 + 1
                     enddo
                     j4 = j3 + 1
                     if ( j4 .ge. no .or. ico_got(j4,mo) .le. 0 ) then
                        j4 = j3
                        nmorat = 2
                     endif
                     j2 = j3 - 1
                     j1 = j2 - 1
c
                     if ( ico_got(j1,mo) .le. 0 ) write(6,4873)
     $                    zat,m,j1,mo,no-1,xcs(j1),1.d0-xa(m)-zat-
     $                    xcs(j1)
                     if ( ico_got(j2,mo) .le. 0 ) write(6,4873)
     $                    zat,m,j2,mo,no-1,xcs(j2),1.d0-xa(m)-zat-
     $                xcs(j2)
                     if ( ico_got(j3,mo) .le. 0 ) write(6,4873)
     $                    zat,m,j3,mo,no-1,xcs(j3),1.d0-xa(m)-zat-
     $                xcs(j3)
 4873                format(' '/'     Z=',f9.7,' ??? CANNOT miss',
     $                    ' mix (m=',i1,',i=',i1,',j=',i1,
     $                    ') with ihi=',i1,' C=',f10.7,' O=',f10.7)
                     if ( ico_got(j1,mo) .le. 0 .or.
     $                    ico_got(j2,mo) .le. 0 .or.
     $                    ico_got(j3,mo) .le. 0 ) stop
     $                    ' STOP -- READCO: mix CANNOT be missing. '
c
c-debug-chk[
c-debug-chk;                     write(6,1975) i,j,'C',cat,'C',cx(i),
c-debug-chk;     $                    nmorat+1,j1,mo,j2,mo,j3,mo,j4,mo,
c-debug-chk;     $                    'C',xcs(j1),xcs(j2),xcs(j3),xcs(j4),
c-debug-chk;     $                    'C',cx(j1),cx(j2),cx(j3),cx(j4)
c-debug-chk; 1975                format(' (',i1,',',i1,') ',a1,'=',f10.6,
c-debug-chk;     $                    ' log',a1,'=',f10.6,' among',i2,
c-debug-chk;     $                    ' of (',i1,',',i1,') (',i1,',',i1,
c-debug-chk;     $                    ') (',i1,',',i1,') (',i1,',',i1,'): ',
c-debug-chk;     $                    a1,'=',4f10.6,' log',a1,'=',4f10.6)
c-debug-chk]
                     is = 0
                     do il = 1, nr
                        do k = 1, nt
                           if ( abs( co(m,j1,mo,k,il,kz) ) .lt.
     $                          badlogklim ) then
                              co(m,i,j,k,il,kz) = qzinter(is,1,cat,
     $                             nmorat,co(m,j1,mo,k,il,kz),
     $                             co(m,j2,mo,k,il,kz),
     $                             co(m,j3,mo,k,il,kz),
     $                             co(m,j4,mo,k,il,kz),xcs(j1),
     $                             xcs(j2),xcs(j3),xcs(j4),zzz(kz))
                              is = 1
c-debug-chk[
c-debug-chk;                              if ( t6list(k) .gt. 0.09999 ) then
c-debug-chk;                                 dif = co(m,i,j,k,il,kz)
c-debug-chk;     $                                -co(m,j3,mo,k,il,kz)
c-debug-chk;                                 difl = co(m,i,j,k,il,kz)
c-debug-chk;     $                                -co(m,j2,mo,k,il,kz)
c-debug-chk;                                 difmax = max(dif,difmax)
c-debug-chk;                                 difmin = min(dif,difmin)
c-debug-chk;                                 sumdif = sumdif+dif
c-debug-chk;                                 sumsq = sumsq+dif**2
c-debug-chk;                                 numsq = numsq+1
c-debug-chk;                                 diflmax = max(difl,diflmax)
c-debug-chk;                                 diflmin = min(difl,diflmin)
c-debug-chk;                                 sumldif = sumldif+difl
c-debug-chk;                                 sumlsq = sumlsq+difl**2
c-debug-chk;                              endif
c-debug-chk]
                           endif
                        enddo
                     enddo
                     ico_got(i,j) = -1
c-debug-chk[
c-debug-chk;                     write(6,8712) numsq,diflmin,diflmax,
c-debug-chk;     $                    sumldif/max(numsq,1),
c-debug-chk;     $                    sqrt(sumlsq/max(numsq,1)),
c-debug-chk;     $                    numsq,difmin,difmax,sumdif/max(numsq,1),
c-debug-chk;     $                    sqrt(sumsq/max(numsq,1))
c-debug-chk; 8712                format(' '/'      ',
c-debug-chk;     $                    ' --- result: relative DIFmid:',
c-debug-chk;     $                    i5,'[',f10.6,' ,',f10.6,' ]ave',f10.6,
c-debug-chk;     $                    ' rms',f10.6,'  DIFhi:',i5,'[',f10.6,
c-debug-chk;     $                    ' ,',f10.6,' ]ave',f10.6,' rms',f10.6)
c-debug-chk]
                  endif
               endif
            enddo
         enddo
c			   ! next, remaining diagonal tables (at Y=0):
c			   !  (only the table at no-1 can possibly be missing)
         do j = 1, no - 2
            if ( ico_got(j,mo) .eq. 0 ) write(6,3873)
     $           zat,m,j,mo,no-1,xcs(j),1.d0-xa(m)-zat-xcs(j)
 3873       format(' '/'     Z=',f9.7,' ??? CANNOT miss',
     $           ' mix (m=',i1,',i=',i1,',j=',i1,
     $           ') with ihi=',i1,' C=',f10.7,' O=',f10.7)
            if ( ico_got(j,mo) .eq. 0 ) stop
     $           ' STOP -- READCO: mix CANNOT be missing. '
         enddo
c
         j = no - 1
c						! if table at no-1 is missing:
         if ( ico_got(j,mo) .eq. 0 ) then
            oat = 1.d0 - xa(m) - zat - xcs(j)
c-debug-chk[
c-debug-chk;            write(6,4973) m,zat,j,mo,no-1,xcs(j),oat
c-debug-chk; 4973       format(' '/'     Z=',f9.7,' --- interpolate',
c-debug-chk;     $           ' mix (m=',i1,',i=',i1,',j=',i1,
c-debug-chk;     $           ') with ihi=',i1,' C=',f10.7,' O=',f10.7)
c-debug-chk;            difmax = -9.999999
c-debug-chk;            difmin = 9.999999
c-debug-chk;            sumdif = 0.
c-debug-chk;            sumsq = 0.
c-debug-chk;            numsq = 0
c-debug-chk;            diflmax = -9.999999
c-debug-chk;            diflmin = 9.999999
c-debug-chk;            sumldif = 0.
c-debug-chk;            sumlsq = 0.
c-debug-chk]
c			! may use quadratic, or two overlapping quadratics
            nmorat = 3
            j3 = 3
            do while( j3 .lt. no - 1 .and. oat .gt. xos(j3) .and.
     $           ico_got(max(n(m,j3+1,kz),1),j3+1) .gt. 0 )
               j3 = j3+1
            enddo
            j4 = j3+1
            if ( j4 .ge. no .or.
     $           ico_got(j4,max(n(m,j4,kz),1)) .le. 0 ) then
               j4 = j3
               nmorat = 2
            endif
            j2 = j3-1
            j1 = j2-1
c-noneed[                                                     ! (checked above)
c-noneed;            if ( ico_got(n(m,j1,kz),j1) .le. 0 .or.
c-noneed;     $           ico_got(n(m,j2,kz),j2) .le. 0 .or.
c-noneed;     $           ico_got(n(m,j3,kz),j3) .le. 0 ) stop
c-noneed;     $           ' STOP -- READCO: mix CANNOT be missing. '
c-noneed]
c-debug-chk[
c-debug-chk;            write(6,2975) j,mo,'O',oat,'O',log10(oat+zzz(kz)),
c-debug-chk;     $           nmorat+1,n(m,j1,kz),j1,n(m,j2,kz),j2,
c-debug-chk;     $           n(m,j3,kz),j3,n(m,j4,kz),j4,
c-debug-chk;     $           'O',xos(j1),xos(j2),xos(j3),xos(j4),
c-debug-chk;     $           'O',ox(j1),ox(j2),ox(j3),ox(j4)
c-debug-chk; 2975       format(' (',i1,',',i1,') ',a1,'=',f10.6,
c-debug-chk;     $           ' log',a1,'=',f10.6,' among',i2,
c-debug-chk;     $           ' of (',i1,',',i1,') (',i1,',',i1,
c-debug-chk;     $           ') (',i1,',',i1,') (',i1,',',i1,'): ',
c-debug-chk;     $           a1,'=',4f10.6,' log',a1,'=',4f10.6)
c-debug-chk]
c			! interpolate along the diagonal (using O-abundance)
            is = 0
            do il = 1, nr
               do k = 1, nt
                  if ( abs( co(m,n(m,j1,kz),j1,k,il,kz) ) .lt.
     $                 badlogklim ) then
                     co(m,j,mo,k,il,kz) = qzinter(is,1,oat,nmorat,
     $                    co(m,n(m,j1,kz),j1,k,il,kz),
     $                    co(m,n(m,j2,kz),j2,k,il,kz),
     $                    co(m,n(m,j3,kz),j3,k,il,kz),
     $                    co(m,n(m,j4,kz),j4,k,il,kz),
     $                    xos(j1),xos(j2),xos(j3),xos(j4),zzz(kz))
                     is = 1
c-debug-chk[
c-debug-chk;                     if ( t6list(k) .gt. 0.09999 ) then
c-debug-chk;                        dif = co(m,j,mo,k,il,kz)
c-debug-chk;     $                       - co(m,n(m,j3,kz),j3,k,il,kz)
c-debug-chk;                        difl = co(m,j,mo,k,il,kz)
c-debug-chk;     $                       - co(m,n(m,j2,kz),j2,k,il,kz)
c-debug-chk;                        difmax = max(dif,difmax)
c-debug-chk;                        difmin = min(dif,difmin)
c-debug-chk;                        sumdif = sumdif+dif
c-debug-chk;                        sumsq = sumsq+dif**2
c-debug-chk;                        numsq = numsq+1
c-debug-chk;                        diflmax = max(difl,diflmax)
c-debug-chk;                        diflmin = min(difl,diflmin)
c-debug-chk;                        sumldif = sumldif+difl
c-debug-chk;                        sumlsq = sumlsq+difl**2
c-debug-chk;                     endif
c-debug-chk]
                  endif
               enddo
            enddo
            ico_got(j,mo) = -1
c-debug-chk[
c-debug-chk;            write(6,7712) numsq,diflmin,diflmax,
c-debug-chk;     $           sumldif/max(numsq,1),
c-debug-chk;     $           sqrt(sumlsq/max(numsq,1)),
c-debug-chk;     $           numsq,difmin,difmax,sumdif/max(numsq,1),
c-debug-chk;     $           sqrt(sumsq/max(numsq,1))
c-debug-chk; 7712       format(' '/'      ',
c-debug-chk;     $           ' --- result: relative DIFmid:',
c-debug-chk;     $           i5,'[',f10.6,' ,',f10.6,' ]ave',f10.6,
c-debug-chk;     $           ' rms',f10.6,'  DIFhi:',i5,'[',f10.6,
c-debug-chk;     $           ' ,',f10.6,' ]ave',f10.6,' rms',f10.6)
c-debug-chk]
         endif
c                                                      0.10-- @ +   +     @
c  If possible, make mixes next to the C=O=0.0       C
c   mix smooth for CO-interpolation (but only if
c   low_CO_smo > 0 in common/c_opal_ctrl_smooth/).     0.03-- @ +   @     +
c   The diagram at right, of the lower-left corner
c   of the C-O plane, shows the mixes that may be      0.01-- * *   +     +
c   interpolated as "*" , the mixes interpolated       0.00-- @ *   @     @
c   among as "@" , and unused mixes as "+".                   | |   |     |
c                                                            0. |  .03   .10
c                                                              .01           O
         if ( low_CO_smo .gt. 0 ) then
            a1 = oxf(m,1,kz)
            a2 = oxf(m,2,kz)
            a3 = oxf(m,3,kz)
            a4 = oxf(m,4,kz)
            is = 0
c-debug-chk[
c-debug-chk;            difmax = -9.999999
c-debug-chk;            difmin = 9.999999
c-debug-chk;            sumdif = 0.
c-debug-chk;            sumsq = 0.
c-debug-chk;            numsq = 0
c-debug-chk;            diflmax = -9.999999
c-debug-chk;            diflmin = 9.999999
c-debug-chk;            sumldif = 0.
c-debug-chk;            sumlsq = 0.
c-debug-chk;            numlsq = 0
c-debug-chk]
c				! loop over the three mixes next to C=O=0.0 mix
            do imix = 1, 3
               ifac = ( 4 - imix ) / 2
               jfac = imix / 2
               do il = 1, nr
                  do k = 1, nt
                     v1 = co(m,1,1,k,il,kz)
                     if ( abs(v1) .lt. badlogklim ) then
                        v2 = co(m,1+ifac,1+jfac,k,il,kz)
                        v3 = co(m,1+2*ifac,1+2*jfac,k,il,kz)
                        v4 = co(m,1+3*ifac,1+3*jfac,k,il,kz)
                        cofmin = min(0.8d0*v1+0.2d0*v3 , 0.2d0*v1+
     $                       0.8d0*v3)
                        cofmax = max(0.8d0*v1+0.2d0*v3 , 0.2d0*v1+
     $                       0.8d0*v3)
                        if ( (v4-v3)*(v3-v1) .gt. 0.d0 .and. 
     $                       ( v2 .lt. cofmin .or.
     $                       v2 .gt. cofmax ) ) then
                           co(m,1+ifac,1+jfac,k,il,kz) = max( min(
     $                          quad(is,1,a2,v1,v3,v4,a1,a3,a4) ,
     $                          cofmax ) , cofmin )
                           is = 1
c-debug-chk[
c-debug-chk;                           dif = co(m,1+ifac,1+jfac,k,il,kz)-v2
c-debug-chk;                           if ( t6list(k) .gt. 0.09999 ) then
c-debug-chk;                              difmax = max(dif,difmax)
c-debug-chk;                              difmin = min(dif,difmin)
c-debug-chk;                              sumdif = sumdif+dif
c-debug-chk;                              sumsq = sumsq+dif**2
c-debug-chk;                              numsq = numsq+1
c-debug-chk;                           else
c-debug-chk;                              diflmax = max(dif,diflmax)
c-debug-chk;                              diflmin = min(dif,diflmin)
c-debug-chk;                              sumldif = sumldif+dif
c-debug-chk;                              sumlsq = sumlsq+dif**2
c-debug-chk;                              numlsq = numlsq+1
c-debug-chk;                           endif
c-debug-chk]
                        endif
                     endif
                  enddo
               enddo
            enddo
c-debug-chk[
c-debug-chk;            write(6,8733) m,zat,numlsq,diflmin,diflmax,
c-debug-chk;     $           sumldif/max(numlsq,1),
c-debug-chk;     $           sqrt(sumlsq/max(numlsq,1)),
c-debug-chk;     $           numsq,difmin,difmax,sumdif/max(numsq,1),
c-debug-chk;     $           sqrt(sumsq/max(numsq,1))
c-debug-chk; 8733       format(' '/' m=',i1,' Z=',f9.7,
c-debug-chk;     $           ' fix C,O=[1,2] T6<0.1:',
c-debug-chk;     $           i5,'[',f10.6,' ,',f10.6,' ]ave',f10.6,
c-debug-chk;     $           ' rms',f10.6,' T6>0.1:',i5,'[',f10.6,' ,',
c-debug-chk;     $           f10.6,' ]ave',f10.6,' rms',f10.6)
c-debug-chk]
         endif
c
c  Peform any specified opacity shifts, from GN93hz and [O/Fe]
c
c						     ! If m=2 shift & no [O/Fe]
         if ( khizat .eq. 1 .and. m .eq. mx03 ) then
c						     !  set all shifts to zero
c						     !   (m=2 GN93hz shift may
c						     !   be interpolated later)
            do il = 1, nr
               do k = 1, nt
                  co(m,mc,mo_m1,k,il,kz) = 0.d0
                  co(m,mc,mo,k,il,kz) = 0.d0
               enddo
            enddo
c					       ! Else, if there are any shifts:
         else if ( khizat .gt. 0 ) then
c					! If there is no [O/Fe] shift:
            if ( khizat .eq. 1 ) then
c					!  then set all [O/Fe] shifts to zero
               do il = 1, nr
                  do k = 1, nt
                     co(m,mc,mo_m1,k,il,kz) = 0.d0
                  enddo
               enddo
c					! Else, if there is the [O/Fe] shift:
            else
c		! get interpolation factors fofe (for GN93hz) and omfofe=1-fofe
c
               xofe = 10.d0**ofebrack * xofe_opalmixes(1)
               fofe = ( xiz_opalmixes(kel_o,khizat)
     $              - xofe * xiz_opalmixes(kel_fe,khizat) )
     $              / ( ( xiz_opalmixes(kel_fe,1)
     $              - xiz_opalmixes(kel_fe,khizat) ) * xofe
     $              + xiz_opalmixes(kel_o,khizat)
     $              - xiz_opalmixes(kel_o,1) )
               omfofe = 1.d0 - fofe
c					! get Z-composition of interpolated mix
               if ( m .eq. mstart ) then
                  do i = 1, nel_zmix
                     xiz_mix(i) = fofe*xiz_opalmixes(i,1)
     $                    + omfofe*xiz_opalmixes(i,khizat)
                     fninz_mix(i) = fofe*fninz_opalmixes(i,1)
     $                    + omfofe*fninz_opalmixes(i,khizat)
                  enddo
                  do i = 1, nel_zmix
                     bracketife_mix(i) = log10(
     $                    ( max(xiz_mix(i),1.d-36)
     $                    * xiz_opalmixes(kel_fe,1) )
     $                    / ( max(xiz_mix(kel_fe),1.d-36)
     $                    * xiz_opalmixes(i,1) ) )
                  enddo
c-debug-chk[
c-debug-chk;                  write(6,2377) ofebrack,khizat,
c-debug-chk;     $                 bracketofe_opalmixes(khizat),
c-debug-chk;     $                 fofe,khizat,omfofe,khizat,khizat
c-debug-chk; 2377             format(' '/' To get mix with [O/Fe] =',f11.7,
c-debug-chk;     $                 ' from mix',i2,' with [O/Fe] =',f11.7,
c-debug-chk;     $                 ': f_(1) =',f11.7,' , f_(',i1,') =',f11.7,
c-debug-chk;     $                 ':'/' '/'  i   Xi/Z_(1)   Ni/Nz_(1)',
c-debug-chk;     $                 '   Xi/Z_(',i1,')   Ni/Nz_(',i1,')',
c-debug-chk;     $                 '   Xi/Z_mix   Ni/Nz_mix   [i/Fe]'/
c-debug-chk;     $                 ' ==  ========== ==========  ==========',
c-debug-chk;     $                 ' ==========  ========== ==========',
c-debug-chk;     $                 ' ==========')
c-debug-chk;                  do i = 1,nel_zmix
c-debug-chk;                     write(6,2376) cel_opalmixes(i),
c-debug-chk;     $                    xiz_opalmixes(i,1),
c-debug-chk;     $                    fninz_opalmixes(i,1),
c-debug-chk;     $                    xiz_opalmixes(i,khizat),
c-debug-chk;     $                    fninz_opalmixes(i,khizat),xiz_mix(i),
c-debug-chk;     $                    fninz_mix(i),bracketife_mix(i)
c-debug-chk; 2376                format(' ',a2,3(f12.7,f11.7),f11.7)
c-debug-chk;                  enddo
c-debug-chk]
               endif
c			   ! compute [O/Fe] shifts relative to GN93hz opacities
c-debug-chk[
c-debug-chk;               sumsq = 0.
c-debug-chk;               sumdif = 0.
c-debug-chk;               difmax = -9.999999
c-debug-chk;               difmin = 9.999999
c-debug-chk;               numsq = 0
c-debug-chk]
               do il = 1, nr
                  do k = nt, 1, -1
                     if ( abs( co(m,mc,mo_m1,k,il,kz) ) .lt. badlogklim
     $                    .and. abs( co(m,mc,mo,k,il,kz) ) .lt.
     $                    badlogklim ) then
                        dif = ( co(m,mc,mo_m1,k,il,kz)
     $                       - co(m,mc,mo,k,il,kz) ) * omfofe
                        co(m,mc,mo_m1,k,il,kz) = dif
c-debug-chk[
c-debug-chk;                        if ( t6list(k) .gt. 0.009999 ) then
c-debug-chk;                           sumdif = sumdif+dif
c-debug-chk;                           sumsq = sumsq+dif**2
c-debug-chk;                           difmax = max(dif,difmax)
c-debug-chk;                           difmin = min(dif,difmin)
c-debug-chk;                           numsq = numsq+1
c-debug-chk;                        endif
c-debug-chk]
                     else if ( k .lt. nt ) then
                        co(m,mc,mo_m1,k,il,kz) =
     $                       co(m,mc,mo_m1,k+1,il,kz)
                     else
                        co(m,mc,mo_m1,k,il,kz) = 0.d0
                     endif
                  enddo
               enddo
c-debug-chk[
c-debug-chk;               write(6,2379) m,zat,numsq,difmin,difmax,
c-debug-chk;     $              sumdif/max(numsq,1),sqrt(sumsq/max(numsq,1)),
c-debug-chk;     $              sqrt(max(sumsq-sumdif**2/max(numsq,1),0.)
c-debug-chk;     $              /max(numsq-1,1))
c-debug-chk; 2379          format(' '/' m=',i1,' Z=',f9.7,
c-debug-chk;     $              ' [O/Fe] deltas for T6>0.01: N=',
c-debug-chk;     $              i4,' DEL[',f10.6,' ,',f10.6,' ] DELave=',f10.6,
c-debug-chk;     $              ' DELrms=',f10.6,' sig',f10.6)
c-debug-chk]
            endif
c			! compute GN93hz shifts relative to Gz???.x?? opacities
c-debug-chk[
c-debug-chk;            difmax = -9.999999
c-debug-chk;            difmin = 9.999999
c-debug-chk;            sumdif = 0.
c-debug-chk;            sumsq = 0.
c-debug-chk;            numsq = 0
c-debug-chk]
c-test-xdel[
c-test-xdel;            do ij = 1,n_xdtst
c-test-xdel;               dif_tst(1,ij) = -9.999999
c-test-xdel;               dif_tst(2,ij) = 9.999999
c-test-xdel;               dif_tst(3,ij) = 0.
c-test-xdel;               dif_tst(4,ij) = 0.
c-test-xdel;            enddo
c-test-xdel]
            do il = 1, nr
               do k = nt, 1, -1
                  if ( abs( co(m,mc,mo,k,il,kz) ) .lt. badlogklim .and.
     $                 abs( co(m,1,1,k,il,kz) ) .lt. badlogklim ) then
                     dif = co(m,mc,mo,k,il,kz) - co(m,1,1,k,il,kz)
c
c							! (reduce GN93hz shifts
c							! by the factor facxhz)
                     co(m,mc,mo,k,il,kz) = dif * facxhz
c-debug-chk[
c-debug-chk;                     if ( t6list(k) .gt. 0.009999 ) then
c-debug-chk;                        difmax = max(dif,difmax)
c-debug-chk;                        difmin = min(dif,difmin)
c-debug-chk;                        sumdif = sumdif+dif
c-debug-chk;                        sumsq = sumsq+dif**2
c-debug-chk;                        numsq = numsq+1
c-debug-chk;                     endif
c-debug-chk]
c-test-xdel[
c-test-xdel;                     if ( t6list(k) .gt. 0.009999 ) then
c-test-xdel;                        if ( nxdo .eq. 3 ) then
c-test-xdel;                           do ij = 1,n_xdtst
c-test-xdel;                              dtst = cof_tst(k+ntdel,il+nrdel,ij)
c-test-xdel;     $                             -co(m,1,1,k,il,kz)
c-test-xdel;                              dif_tst(1,ij) = max(dif_tst(1,ij),
c-test-xdel;     $                             dtst)
c-test-xdel;                              dif_tst(2,ij) = min(dif_tst(2,ij),
c-test-xdel;     $                             dtst)
c-test-xdel;                              dif_tst(3,ij) = dif_tst(3,ij)+dtst
c-test-xdel;                              dif_tst(4,ij) = dif_tst(4,ij)+dtst**2
c-test-xdel;                           enddo
c-test-xdel;                        endif
c-test-xdel;                     endif
c-test-xdel]
                  else if ( k .lt. nt ) then
                     co(m,mc,mo,k,il,kz) = co(m,mc,mo,k+1,il,kz)
                  else
                     co(m,mc,mo,k,il,kz) = 0.d0
                  endif
               enddo
            enddo
c-debug-chk[
c-debug-chk;            write(6,2378) m,zat,numsq,difmin,difmax,
c-debug-chk;     $           sumdif/max(numsq,1),sqrt(sumsq/max(numsq,1)),
c-debug-chk;     $           sqrt(max(sumsq-sumdif**2/max(numsq,1),0.)
c-debug-chk;     $           /max(numsq-1,1)),facxhz
c-debug-chk; 2378       format(' '/' m=',i1,' Z=',f9.7,
c-debug-chk;     $           ' GN93hz deltas for T6>0.01: N=',
c-debug-chk;     $           i4,' DEL[',f10.6,' ,',f10.6,' ] DELave=',f10.6,
c-debug-chk;     $           ' DELrms=',f10.6,' sig',f10.6,
c-debug-chk;     $           ' reduced by facxhz=',f10.7)
c-debug-chk]
c-test-xdel[
c-test-xdel;            if ( nxdo .eq. 3 ) then
c-test-xdel;               do ij = 1,n_xdtst
c-test-xdel;                  write(6,5817) numsq,dif_tst(1,ij),dif_tst(2,ij),
c-test-xdel;     $                 dif_tst(3,ij)/max(numsq,1),
c-test-xdel;     $                 sqrt(dif_tst(4,ij)/max(numsq,1)),
c-test-xdel;     $                 sqrt(max(dif_tst(4,ij)-dif_tst(3,ij)**2
c-test-xdel;     $                 /max(numsq,1),0.)/max(numsq-1,1)),
c-test-xdel;     $                 xdel_tst(ij)
c-test-xdel; 5817             format('                ',
c-test-xdel;     $                 ' GN93hz deltas for T6>0.01: N=',i4,
c-test-xdel;     $                 ' DEL[',f10.6,' ,',f10.6,' ] DELave=',f10.6,
c-test-xdel;     $                 ' DELrms=',f10.6,' sig',f10.6,
c-test-xdel;     $                 ' for Xdel=',f6.4)
c-test-xdel;               enddo
c-test-xdel;            endif
c-test-xdel]
         endif
c-debug-chk[
c-debug-chk;         do i = 1,no-1
c-debug-chk;            oat = 1.-xa(m)-zat-xcs(i)
c-debug-chk;            io = -1
c-debug-chk;            do j = 1,no-1
c-debug-chk;               ihi = n(m,j,kz)
c-debug-chk;               cat = min(xcs(ihi),1.-xa(m)-zat-xos(j))
c-debug-chk;               if ( max( abs(xcs(i)-cat) , abs(oat-xos(j)) )
c-debug-chk;     $              .lt. 0.0011 ) then
c-debug-chk;                  io = ihi
c-debug-chk;                  jo = j
c-debug-chk;               endif
c-debug-chk;            enddo
c-debug-chk;            if ( io .gt. 0 ) then
c-debug-chk;               difmax = -9.999999
c-debug-chk;               difmin = 9.999999
c-debug-chk;               sumdif = 0.
c-debug-chk;               sumsq = 0.
c-debug-chk;               numsq = 0
c-debug-chk;               do il = 1,nr
c-debug-chk;                  do k = 6,nta(il+nrdel)-ntdel
c-debug-chk;                     dif = co(m,io,jo,k,il,kz)-co(m,i,mo,k,il,kz)
c-debug-chk;                     difmax = max(dif,difmax)
c-debug-chk;                     difmin = min(dif,difmin)
c-debug-chk;                     sumdif = sumdif+dif
c-debug-chk;                     sumsq = sumsq+dif**2
c-debug-chk;                     numsq = numsq+1
c-debug-chk;                  enddo
c-debug-chk;               enddo
c-debug-chk;               write(6,1598) m,zat,io,jo,i,mo,numsq,difmin,difmax,
c-debug-chk;     $              sumdif/max(numsq,1),sqrt(sumsq/max(numsq,1)),
c-debug-chk;     $              min(xcs(io),1.-xa(m)-z-xos(jo)),xos(jo),
c-debug-chk;     $              xcs(i),1.-xa(m)-z-xcs(i)
c-debug-chk; 1598          format(' '/' m=',i1,' Z=',f9.7,' d:(',i1,',',i1,
c-debug-chk;     $              ')-(',i1,',',i1,') for T6>0.01: N=',i4,' DIF[',
c-debug-chk;     $              f10.6,' ,',f10.6,' ] DIFave=',f10.6,' DIFrms=',
c-debug-chk;     $              f10.6,' CO',2f10.7,' &',2f10.7)
c-debug-chk;            endif
c-debug-chk;         enddo
c-debug-chk]
c			! End of loop over m values
      enddo
c
c }--------------------------------------- End of loop over m values (X-tables)
c
c-debug-chk[
c-debug-chk;      write(6,8418) (i,(n(i,j,kz),j=1,mo),' ',i=mstart,mend)
c-debug-chk; 8418 format(' '/' -- n(m,j): ',5(' (m=',i1,')',8i2,a1))
c-debug-chk]
c		! interpolate GN93hz opacity shifts for m=2, if possible; note
c		! that other shifts being interpolated among already contain
c		! the factor of facxhz. No need to revise any m=2 [O/Fe] shift.
c
      if ( kallrd .ne. 0 .and. khizat .ge. 1 .and.
     $     mx .ge. 4 .and. mxzero .eq. 1 .and. mx03 .eq. 2 ) then
c
c-debug-chk[
c-debug-chk;         sumsq = 0.
c-debug-chk;         sumdif = 0.
c-debug-chk;         difmax = -9.999999
c-debug-chk;         difmin = 9.999999
c-debug-chk;         numsq = 0
c-debug-chk]
         is = 0
c			! for all densities and temperatures
         do il = 1, nr
            do k = nt, 1, -1
c					! if it is possible to interpolate
c
               if ( abs( co(1,1,1,k,il,kz) ) .lt. badlogklim .and.
     $              abs( co(2,1,1,k,il,kz) ) .lt. badlogklim ) then
c
c							     ! new GN93hz shift
                  dif = quad(is,1,xx(2),co(1,mc,mo,k,il,kz),
     $                 co(3,mc,mo,k,il,kz),co(4,mc,mo,k,il,kz),
     $                 xx(1),xx(3),xx(4))
                  is = 1
c-debug-chk[
c-debug-chk;                  if ( t6list(k) .gt. 0.009999 ) then
c-debug-chk;                     sumdif = sumdif+dif
c-debug-chk;                     sumsq = sumsq+dif**2
c-debug-chk;                     difmax = max(dif,difmax)
c-debug-chk;                     difmin = min(dif,difmin)
c-debug-chk;                     numsq = numsq+1
c-debug-chk;                  endif
c-debug-chk]
                  co(2,mc,mo,k,il,kz) = dif
               else if ( k .lt. nt ) then
                  co(2,mc,mo,k,il,kz) = co(2,mc,mo,k+1,il,kz)
               else
                  co(2,mc,mo,k,il,kz) = 0.d0
               endif
c
            enddo
         enddo
c-debug-chk[
c-debug-chk;         if ( facxhz .gt. 0. .and. facxhz .lt. 1.  ) then
c-debug-chk;            difmin = difmin/facxhz
c-debug-chk;            difmax = difmax/facxhz
c-debug-chk;            sumdif = sumdif/facxhz
c-debug-chk;            sumsq = sumsq/facxhz**2
c-debug-chk;         else
c-debug-chk;            facxhz = 1.
c-debug-chk;         endif
c-debug-chk;         write(6,2371) z,numsq,difmin,difmax,
c-debug-chk;     $        sumdif/max(numsq,1),sqrt(sumsq/max(numsq,1)),facxhz
c-debug-chk; 2371    format(' '/' m=2 Z=',f9.7,
c-debug-chk;     $        ' GN93hz alt-deltas T6>0.01: N=',
c-debug-chk;     $        i4,' DIF[',f10.6,' ,',f10.6,' ] DIFave=',f10.6,
c-debug-chk;     $        ' DIFrms=',f10.6,' reduced by facxhz=',f10.7)
c-debug-chk]
c		 ! end of interpolation of GN93hz opacity shifts for m=2
      endif
c				! apply all opacity shifts calculated above
      if ( khizat .gt. 0 ) then
c					! Begin loop over m values:
         do mvalue = mstart, mend
c
            m = mvalue
c-debug-chk[
c-debug-chk;            difmax = -9.999999
c-debug-chk;            difmin = 9.999999
c-debug-chk;            sumdif = 0.
c-debug-chk;            sumsq = 0.
c-debug-chk;            numsq = 0
c-debug-chk;            difcmax = -9.999999
c-debug-chk;            difcmin = 9.999999
c-debug-chk;            sumcdif = 0.
c-debug-chk;            sumcsq = 0.
c-debug-chk;            numcsq = 0
c-debug-chk;            diflmax = -9.999999
c-debug-chk;            diflmin = 9.999999
c-debug-chk;            sumldif = 0.
c-debug-chk;            sumlsq = 0.
c-debug-chk;            numlsq = n(m,1,kz)-2
c-debug-chk;            do j = 1,n(m,1,kz)-1
c-debug-chk;               numlsq = numlsq+n(m,j,kz)
c-debug-chk;            enddo
c-debug-chk;            difomax = -9.999999
c-debug-chk;            difomin = 9.999999
c-debug-chk;            difdmax = 0.
c-debug-chk;            numosq = 0
c-debug-chk;            num1 = 0
c-debug-chk]
c				! perform the opacity shifts computed above
            do il = 1, nr
               do k = 1, nta(il+nrdel) - ntdel
c
                  dif = co(m,mc,mo,k,il,kz) + co(m,mc,mo_m1,k,il,kz)
c-debug-chk[
c-debug-chk;                  if ( t6list(k) .gt. 0.009999 ) then
c-debug-chk;                     difcmax = max(dif,difcmax)
c-debug-chk;                     difcmin = min(dif,difcmin)
c-debug-chk;                     sumcdif = sumcdif+dif
c-debug-chk;                     sumcsq = sumcsq+dif**2
c-debug-chk;                     numcsq = numcsq+1
c-debug-chk;                     if ( dif .eq. 0. ) then
c-debug-chk;                        numsq = numsq+numlsq
c-debug-chk;                        num1 = num1+numlsq
c-debug-chk;                        diflmax = max(0.,diflmax)
c-debug-chk;                        diflmin = min(0.,diflmin)
c-debug-chk;                        difmax = max(0.,difmax)
c-debug-chk;                        difmin = min(0.,difmin)
c-debug-chk;                     endif
c-debug-chk;                  endif
c-debug-chk]
                  if ( dif .ne. 0.d0 ) then
c
                     if ( abs(dif) .gt. 0.01d0 ) then
                        diffac = 10.d0**dif - 1.d0
                     else
c						! (more accurate for small dif)
                        difl = 2.302585d0 * dif
                        diffac = ((difl*.16666667d0+0.5d0)*difl+1.d0)*
     $                       difl
                     endif
c				! first, do shifts for all except C=O=0.0 mix:
                     ilo = 2
                     do j = 1, mo
                        if ( j .lt. n(m,1,kz) ) then
                           ihi = n(m,j,kz)
                        else if ( j .eq. mo ) then
                           ihi = n(m,1,kz) - 1
                        else
                           ihi = 0
                        endif
                        if ( ihi .gt. 0 ) then
                           do i = ilo, ihi
c					   ! If Kappa(COrich) does not exceed
c					   ! Kappa(C=O=0): then apply shift as
c					   ! delta{Log(Kappa)}, the fractional
c					   ! shift in Kappa, which may be a
c					   ! smaller shift delta'{Kappa} in the
c							 ! opacity Kappa itself
                              if ( co(m,i,j,k,il,kz) .le.
     $                             co(m,1,1,k,il,kz) ) then
                                 difl = dif
c					   ! Else Kappa(COrich) > Kappa(C=O=0):
c					   ! convert delta{Log(Kappa)} to an
c					   ! absolute shift delta{Kappa} and
c					   ! apply this (taking the log again);
c					   ! this will be a smaller fractional
c					   ! shift delta'{log(Kappa)}
                              else
                                 difl = log10(10.d0**( co(m,1,1,k,il,kz)
     $                                - co(m,i,j,k,il,kz) ) * diffac
     $                                + 1.d0 )
c-debug-chk[
c-debug-chk;                                 if ( abs(difl) .gt.
c-debug-chk;     $                                abs(dif) ) then
c-debug-chk;                                    difdmax = max(
c-debug-chk;     $                                   abs(difl)-abs(dif),
c-debug-chk;     $                                   difdmax)
c-debug-chk;                                    difomax = max(difomax,dif)
c-debug-chk;                                    difomin = min(difomin,dif)
c-debug-chk;                                    numosq = numosq+1
c-debug-chk;                                 endif
c-debug-chk]
c					! this can only happen to the extent of
c							       ! roundoff error
                                 if ( abs(difl) .gt. abs(dif) )
     $                                difl = dif
                              endif
c								 ! apply shift
                              co(m,i,j,k,il,kz) =
     $                             co(m,i,j,k,il,kz) + difl
c-debug-chk[
c-debug-chk;                              if ( t6list(k) .gt. 0.009999 ) then
c-debug-chk;                                 difmax = max(difl,difmax)
c-debug-chk;                                 difmin = min(difl,difmin)
c-debug-chk;                                 sumdif = sumdif+difl
c-debug-chk;                                 sumsq = sumsq+difl**2
c-debug-chk;                                 numsq = numsq+1
c-debug-chk;                                 if ( difl .eq. dif )
c-debug-chk;     $                                num1 = num1+1
c-debug-chk;                                 difl = (difl-dif)/dif
c-debug-chk;                                 diflmax = max(difl,diflmax)
c-debug-chk;                                 diflmin = min(difl,diflmin)
c-debug-chk;                                 sumldif = sumldif+difl
c-debug-chk;                                 sumlsq = sumlsq+difl**2
c-debug-chk;                              endif
c-debug-chk]
                           enddo
                        endif
                        ilo = 1
                     enddo
c						       ! now shift C=O=0.0 mix:
c
                     co(m,1,1,k,il,kz) = co(m,1,1,k,il,kz) + dif
c
                  endif
c
               enddo
            enddo
c-debug-chk[
c-debug-chk;            write(6,9782) m,numcsq,difcmin,difcmax,sumcdif
c-debug-chk;     $           /max(numcsq,1),sqrt(sumcsq/max(numcsq,1)),
c-debug-chk;     $           sqrt(max(sumcsq-sumcdif**2/max(numcsq,1),0.)
c-debug-chk;     $           /max(numcsq-1,1)),z
c-debug-chk; 9782       format(' '/' m=',i1,' total deltas C+O=0, T6>0.01:',
c-debug-chk;     $           i6,' [',f10.6,' ,',f10.6,' ]ave',f10.6,
c-debug-chk;     $           ' rms',f10.6,' sig',f10.6,'  for Z=',f10.7)
c-debug-chk;            write(6,8782) m,numsq,difmin,difmax,
c-debug-chk;     $           sumdif/max(numsq,1),sqrt(sumsq/max(numsq,1)),
c-debug-chk;     $           diflmin+1.,num1,diflmax+1.,
c-debug-chk;     $           sumldif/max(numsq,1)+1.,sqrt(max(sumlsq
c-debug-chk;     $           -sumldif**2/max(numsq,1),0.)/max(numsq-1,1))
c-debug-chk; 8782       format(' '/' m=',i1,' total deltas C+O>0, T6>0.01:',
c-debug-chk;     $           i6,' [',f10.6,' ,',f10.6,' ]ave',f10.6,
c-debug-chk;     $           ' rms',f10.6,'  freduce[',f10.6,' ,',i6,':',
c-debug-chk;     $           f10.6,' ]ave',f10.6,' sig',f10.6)
c-debug-chk;            if ( numosq .gt. 0 ) write(6,8783) numosq,difdmax,
c-debug-chk;     $           difomin,difomax
c-debug-chk; 8783       format(' '/i23,
c-debug-chk;     $           ' Kco > K0 cases where log(linear delta)',
c-debug-chk;     $           ' > old log delta, by up to',f13.9,
c-debug-chk;     $           ' for deltas as large as [',f13.9,' ,',f13.9,' ]')
c-debug-chk]
c					! end of loop over m-values
         enddo
c					! end of opacity shifts
      endif
c
c  If  khighx > 0  then one should set the corrections at the 'GN93hz' X-values
c  that are not present in the 'Gz???.x??' files.
c
      if ( khighx(kz) .le. 0 ) then
c					! set flags showing high-X unavailable
         kavail_xhi = 0
         kdo_xhi = 0
c
      else
c				! set flags showing whether high-X is available
         kavail_xhi = 1
         do i = 1, numz
            kavail_xhi = min( khighx(i) , kavail_xhi )
         enddo
         if ( kavail_xhi .le. 0 ) then
            kdo_xhi = 0
         else
            kdo_xhi = kuse_xhi
         endif
c-debug-chk[
c-debug-chk;         do i = 1, 20
c-debug-chk;            chk_max(i) = -9.
c-debug-chk;            chk_min(i) = 9.
c-debug-chk;            chk_sum(i) = 0.
c-debug-chk;            chk_ssq(i) = 0.
c-debug-chk;            n_chk(i) = 0
c-debug-chk;         enddo
c-debug-chk]
c		 ! get the 'GN93hz' Z-indices (may not have been done above)
         zat = z
         kzalbe = mzal
         do while( kzalbe .gt. 1 .and. z .le. zalval(kzalbe)-1.d-6 )
            kzalbe = kzalbe-1
         enddo
         if ( abs( zalval(kzalbe) - z ) .le. zacc(kz) ) then
            zat = zalval(kzalbe)
            kzalow = kzalbe
            nzalmo = 0
         else
            kzalow = max( 1 , kzalbe - 1 )
            nzalmo = min( kzalbe + 2 , mzal ) - kzalow
         endif
         int_hi_z = 0
c			       ! set the directory-part of the opacity filename
         if ( kope .eq. 0 ) then
            copfil = ' '
         else
            copfil = copdir(1:kope)
         endif
c								! get filename
         copfil(kope+1:80) = cfile_opalmixes(1)
         iu = iulow
c								! open file
         call open_chk_zip( iu, copfil, igzip,
     $        'STOP -- Error: hz-file (C+O=0.0,[O/Fe]=0) not found.' )
c
         itab_dum = 0
         khighx(kz) = 1
c								! [O/Fe] > 0 ?
         if ( khighz .gt. 1 .and. khighz .le. 5 .and.
     $        kzalow + nzalmo .gt. 1 ) then
            khighx(kz) = 2
            copfil(kope+1:80) = cfile_opalmixes(khighz)
            iu_ofe = iu + 1
            call open_chk_zip( iu_ofe, copfil, igzip_ofe,
     $           'STOP -- Error: hz-file (C+O=0.0,[O/Fe]>0) not found.'
     $           )
            itab_dum_ofe = 0
         endif
c			! how many X-values to store for the present value of Z
         mx_use = mx_hi
         do while ( xhi_in(mx_use-1) .gt. 0.999999d0 - zat )
            mx_use = mx_use - 1
         enddo
         nx_hi(kz) = mx_use
c
         ix = 0
         m = 0
         io = mo
         iz_hi = nzalmo
c						! loop over 'GN93hx' X-values:
c
         do while ( ix .lt. mx_use .and. iz_hi .lt. 5 )
c							! get position in co()
            ix = ix + 1
            m = m + 1
            if ( m .gt. mx ) then
               io = io - 1
               m = 1
            endif
c				! get Z and X values to look for in 'GN93hx'
            iz_hi = nzalmo
            if ( ix .eq. mx_use ) then
               do iz = 0, nzalmo
                  zhi_look(iz+1) = zalval(kzalow+iz)
                  xhi_look(iz+1) = 1.d0 - zhi_look(iz+1)
               enddo
            else
               do iz = 0, nzalmo
                  zhi_look(iz+1) = zalval(kzalow+iz)
                  xhi_look(iz+1) = min( xhi_in(ix) ,
     $                 1.d0 - zhi_look(iz+1) )
               enddo
c			! check for X-column bifurcation at Z = 0.05, X = 0.95
c
               if ( ix .eq. mx_hi - 1 .and. nzalmo .eq. 3 ) then
                  if ( zat .gt. 0.03d0 .and. zat .lt. 0.04d0 ) then
                     iz_hi = 2
                     int_hi_z = 0
                  else if (zat .gt. 0.04d0 .and. zat .lt. 0.05d0) then
                     iz_hi = 5
                     int_hi_z = 0
                     zhi_look(5) = zhi_look(3)
                     xhi_look(5) = xhi_look(3)
                     zhi_look(6) = zhi_look(4)
                     xhi_look(6) = xhi_look(4)
                     zhi_look(3) = zalval(kzalow)
                     xhi_look(3) = 1.d0 - zhi_look(3)
                     zhi_look(4) = zalval(kzalow+1)
                     xhi_look(4) = 1.d0 - zhi_look(4)
                  endif
               endif
            endif
c				! loop over the required Z-values for this X
            do iz = 0, iz_hi
c
               kat = iz + 1
c					       ! find mix; stop if not found
               i_rewind = 0
               igetzxi = 0
               ifound = mixfind(iu,1,igetzxi,i_rewind,
     $              itab_dum,zhi_look(kat),xhi_look(kat),0.0d0,0.0d0)
               if ( ifound .eq. 0 ) then
                  i_rewind = 1
                  igetzxi = 0
                  ifound = mixfind(iu,1,igetzxi,i_rewind,
     $                 itab_dum,zhi_look(kat),xhi_look(kat),0.0d0,0.0d0)
                  if ( ifound .eq. 0 ) then
                     write(6,1791) zhi_look(kat),xhi_look(kat),
     $                    0.0d0,0.0d0,cfile_opalmixes(1)
                     stop ' STOP -- READCO: error reading hz-mix. '
                  endif
               endif
               if ( khighx(kz) .gt. 1 ) then
                  i_rewind = 0
                  igetzxi = 0
                  ifound = mixfind(iu_ofe,1,igetzxi,i_rewind,
     $                 itab_dum_ofe,zhi_look(kat),xhi_look(kat),
     $                 0.0d0,0.0d0)
                  if ( ifound .eq. 0 ) then
                     i_rewind = 1
                     igetzxi = 0
                     ifound = mixfind(iu_ofe,1,igetzxi,i_rewind,
     $                    itab_dum_ofe,zhi_look(kat),xhi_look(kat),
     $                    0.0d0,0.0d0)
                     if ( ifound .eq. 0 ) then
                        write(6,1791) zhi_look(kat),xhi_look(kat),
     $                       0.0d0,0.0d0,cfile_opalmixes(khighz)
                        stop ' STOP -- READCO: error reading hz-mix. '
                     endif
                  endif
               endif
c				  ! loop over logT values, to read in opacities
               do k = 1, ntm
c					   ! read logT, & logKappa(R) for all R
                  read(iu,8300) cin
c!! 8300             format(a137)
 8300             format(a200)
                  read(cin,8140) flt, (cofzhi(k,il,kat),il=1,nrm)
c!! 8140             format(f4.2,19f7.3)
 8140             format(f4.2,28f7.3)
c								   ! bad logT ?
                  if ( abs(flogtin(k)-flt) .gt. 1.d-5 ) stop
     $                 ' STOP -- READCO: bad logT value. '
c
                  il_lo = 1
                  il_hi = nrm
c							      ! logKappa(R) is:
                  do il = nrm, 1, -1
c								       ! absent
                     if ( cin(7*il-2:7*il+4) .eq. '       ' ) then
                        if ( k .le. max(nta(il),nta(0)) ) stop
     $                       ' STOP -- READCO: bad upper edge. '
                        il_hi = min( il_hi , il - 1 )
c							     ! should be absent
                     else if ( k .gt. nta(il) .and.
     $                       il .ge. nrb .and. il .le. nre ) then
                        stop ' STOP -- READCO: bad upper edge. '
c									! 9.999
                     else if ( cofzhi(k,il,kat) .gt. 9.d0 ) then
                        if ( ix .ne. 1 ) stop
     $                       ' STOP -- READCO: bad low edge [O/Fe]=0. '
                        il_lo = max( il_lo , il + 1 )
                     endif
                  enddo
c					      ! also read [O/Fe] > 0, if needed
                  if ( khighx(kz) .gt. 1 ) then
                     read(iu_ofe,8300) cin
                     read(cin,8140) flt, (coff(k,il),il=1,nrm)
                     if ( abs(flogtin(k)-flt) .gt. 1.d-5 ) stop
     $                    ' STOP -- READCO: bad logT value. '
                     do il = nrm, 1, -1
                        if ( cin(7*il-2:7*il+4) .eq. '       ' ) then
                           if ( k .le. max(nta(il),nta(0)) ) stop
     $                          ' STOP -- READCO: bad upper edge. '
                           il_hi = min( il_hi , il - 1 )
                        else if ( k .gt. nta(il) .and.
     $                          il .ge. nrb .and. il .le. nre ) then
                           stop ' STOP -- READCO: bad upper edge. '
                        else if ( coff(k,il) .gt. 9.d0 ) then
                           if ( ix .ne. 1 ) stop
     $                          ' STOP -- READCO: bad low edge. '
                           il_lo = max( il_lo , il + 1 )
                        endif
                     enddo
                     do il = 1, nrm
                        cofzhi(k,il,kat) = fofe * cofzhi(k,il,kat)
     $                       + omfofe * coff(k,il)
                     enddo
                  endif
c								! for smoothing
                  if ( il_lo .gt. 1 .or. il_hi .lt. nrm ) then
                     do il = nrm, 1, -1
                        if ( il .gt. il_hi ) then
                           cofzhi(k,il,kat) = 2.d0 * cofzhi(k-1,il,kat)
     $                          - cofzhi(k-2,il,kat)
                        else if ( il .lt. il_lo ) then
                           cofzhi(k,il,kat) = 2.d0 * cofzhi(k,il+1,kat)
     $                          - cofzhi(k,il+2,kat)
                        endif
                     enddo
                  endif
c			! (end of loop to read opacities at all T-values):
               enddo
c			! (end of loop over required Z-values for this X):
            enddo
c							 ! actual X at Zsto(kz)
            xhi_use(ix,kz) = min( xhi_in(ix) , 1.d0 - zat )
c							     ! Z-interpolation:
            if ( iz_hi .le. 3 ) then
c					! standard case: for all T,R:
               do k = 1, ntm
                  do il = 1, nrm
c							       ! logK at Zsto,X
                     coff(k,il) = qzinter(int_hi_z,1,zat,iz_hi,
     $                    cofzhi(k,il,1),cofzhi(k,il,2),
     $                    cofzhi(k,il,3),cofzhi(k,il,4),
     $                    zhi_look(1),zhi_look(2),zhi_look(3),
     $                    zhi_look(4),zdel)
                     int_hi_z = 1
                  enddo
               enddo
c					    ! ELSE: bifurcation:
            else
c					    ! do both X = 1-Z and X = 0.95
               xhi_use(mx_hi,kz) = 1.d0 - zat
c						! for all T,R:
               do k = 1, ntm
                  do il = 1, nrm
c							   ! logK at Zsto,X=1-Z
                     coff(k,il) = qzinter(int_hi_z,1,zat,3,
     $                    cofzhi(k,il,3),cofzhi(k,il,4),
     $                    cofzhi(k,il,5),cofzhi(k,il,6),
     $                    zhi_look(3),zhi_look(4),zhi_look(5),
     $                    zhi_look(6),zdel)
c							  ! temp: Z=0.05,X=0.95
                     cof_tmp = qzinter(int_hi_z,2,0.05d0,3,
     $                    cofzhi(k,il,3),cofzhi(k,il,4),
     $                    cofzhi(k,il,5),cofzhi(k,il,6),
     $                    zhi_look(3),zhi_look(4),zhi_look(5),
     $                    zhi_look(6),zdel)
c							      ! Zsto(kz),X=0.95
                     cofzhi(k,il,1) = qzinter(int_hi_z,3,zat,2,
     $                    cofzhi(k,il,1),cofzhi(k,il,2),
     $                    cof_tmp,0.0d0,
     $                    zhi_look(1),zhi_look(2),0.05d0,0.0d0,zdel)
                     int_hi_z = 1
                  enddo
               enddo
c					 ! smooth hz-opacities, if init_smo > 0
               if ( init_smo .gt. 0 ) then
                  tmax = 10.d0
                  nset = ks81
                  RLS = alrf(1)
                  RLE = alrf(nrm)
                  nrhigh = int( dfsr(nr) * ( RLE - RLS ) + 1.00001d0 )
                  nrlow = 1
                  call opaltab
               endif
c					     ! store X=1-Z hz-opacity set
               do il = 1, nr
                  jl = il + nrdel
                  do k = 1, nt
                     co(mx,mc,mo_m1,k,il,kz) = coff(k+ntdel,jl)
                  enddo
               enddo
c			     ! prepare to smooth present X=0.95 hz-opacity set
               do k = 1, ntm
                  do il = 1, nrm
                     coff(k,il) = cofzhi(k,il,1)
                  enddo
               enddo
c
            endif
c					 ! smooth hz-opacities, if init_smo > 0
            if ( init_smo .gt. 0 ) then
               tmax = 10.d0
               nset = ks81
               RLS = alrf(1)
               RLE = alrf(nrm)
               nrhigh = int( dfsr(nr) * ( RLE - RLS ) + 1.00001d0 )
               nrlow = 1
               call opaltab
            endif
c					     ! store present hz-opacity set
            do il = 1, nr
               jl = il + nrdel
               do k = 1, nt
                  co(m,mc,io,k,il,kz) = coff(k+ntdel,jl)
               enddo
            enddo
c			! (end of loop over 'GN93hx' X-values):
         enddo
c						! close 'GN93hz' file(s)
         if ( khighx(kz) .gt. 1 ) then
            call close_chk_zip( iu_ofe, copfil, igzip_ofe )
            copfil(kope+1:80) = cfile_opalmixes(1)
         endif
         call close_chk_zip( iu, copfil, igzip )
c							! some needed values
         xxx_max = log10( 1.d0 - zat + xdel )
         xxx_03 = log10( 0.03d0 + xdel )
         f_3 = ( xx(4) - xxx_hi(3) ) * dfsx(4)
         omf_3 = 1.d0 - f_3
c					! for convenience, define ALL xhi_use
         if ( mx_use .lt. mx_hi ) then
            do ix = mx_use + 1, mx_hi
               xhi_use(ix,kz) = xhi_use(mx_use,kz)
            enddo
         endif
c			! get the dlogKappa values for 'GN93hx' X-values:
         int_hi_z = 0
         int_hi_1 = 0
         int_hi_2 = 0
         int_hi_3 = 0
c		       ! loop over all densities and temperatures
         do il = 1, nr
            jl = il + nrdel
            do k = 1, nt
c					     ! if in high-T,R cutout: no shift:
               if ( k .gt. nta(jl) ) then
                  do ix = 1, mx
                     co(ix,mc,mo,k,il,kz) = 0.0d0
                     co(ix,mc,mo_m1,k,il,kz) = 0.0d0
                  enddo
c				       ! else: compute shifts for all X-values:
               else
c								   ! bad logK ?
                  if ( max( co(3,1,1,k,il,kz) , co(4,1,1,k,il,kz) ,
     $                 co(5,1,1,k,il,kz) ) .gt. 9.d0 ) stop
     $                 ' STOP -- Error: bad co(3:5,1,1,*,*) cannot. '
c-debug-chk[
c-debug-chk;                  if ( k+ntdel .gt. 5 ) then
c-debug-chk;                     cof_del = co(2,mc,mo,k,il,kz)
c-debug-chk;     $                       - co(3,1,1,k,il,kz)
c-debug-chk;                     n_chk(12) = n_chk(12) + 1
c-debug-chk;                     chk_max(12) = max( chk_max(12) , cof_del )
c-debug-chk;                     chk_min(12) = min( chk_min(12) , cof_del )
c-debug-chk;                     chk_sum(12) = chk_sum(12) + cof_del
c-debug-chk;                     chk_ssq(12) = chk_ssq(12) + cof_del**2
c-debug-chk;                     cof_del = co(4,mc,mo,k,il,kz)
c-debug-chk;     $                       - co(4,1,1,k,il,kz)
c-debug-chk;                     n_chk(14) = n_chk(14) + 1
c-debug-chk;                     chk_max(14) = max( chk_max(14) , cof_del )
c-debug-chk;                     chk_min(14) = min( chk_min(14) , cof_del )
c-debug-chk;                     chk_sum(14) = chk_sum(14) + cof_del
c-debug-chk;                     chk_ssq(14) = chk_ssq(14) + cof_del**2
c-debug-chk;                     cof_del = co(1,mc,mo_m1,k,il,kz)
c-debug-chk;     $                    - co(5,1,1,k,il,kz)
c-debug-chk;                     n_chk(16) = n_chk(16) + 1
c-debug-chk;                     chk_max(16) = max( chk_max(16) , cof_del )
c-debug-chk;                     chk_min(16) = min( chk_min(16) , cof_del )
c-debug-chk;                     chk_sum(16) = chk_sum(16) + cof_del
c-debug-chk;                     chk_ssq(16) = chk_ssq(16) + cof_del**2
c-debug-chk;                  endif
c-debug-chk]
c							! if no logK at X=0.03:
                  if ( co(2,1,1,k,il,kz) .gt. 9.d0 .or.
     $                 k+ntdel .lt. ntax03(nr+nrdel) ) then
c							    ! 3-X-pt at X=0.2
                     cof_tmp = quad(int_hi_1,1,xxx_hi(3),
     $                    co(2,mc,mo,k,il,kz),co(4,mc,mo,k,il,kz),
     $                    co(1,mc,mo_m1,k,il,kz),
     $                    xxx_hi(2),xxx_hi(4),xxx_hi(6))
                     co(3,mc,mo,k,il,kz) = co(3,mc,mo,k,il,kz)
     $                    - cof_tmp
                     int_hi_1 = 1
c-debug-chk[
c-debug-chk;                     cof_o = quad(1,1,xxx_hi(3),co(3,1,1,k,il,kz),
c-debug-chk;     $                    co(4,1,1,k,il,kz),co(5,1,1,k,il,kz),
c-debug-chk;     $                    xxx_hi(2),xxx_hi(4),xxx_hi(6))
c-debug-chk;                     cof_del = cof_tmp - cof_o
c-debug-chk]
c					   ! else: if logK available at X=0.03:
                  else
c-debug-chk[
c-debug-chk;                     if ( k+ntdel .gt. 5 .and.
c-debug-chk;     $                    k+ntdel .lt. ntax0(nr+nrdel) .and.
c-debug-chk;     $                    co(1,1,1,k,il,kz) .lt. 9. ) then
c-debug-chk;                        cof_del = co(1,mc,mo,k,il,kz)
c-debug-chk;     $                       - co(1,1,1,k,il,kz)
c-debug-chk;                        n_chk(11) = n_chk(11) + 1
c-debug-chk;                        chk_max(11) = max( chk_max(11) , cof_del )
c-debug-chk;                        chk_min(11) = min( chk_min(11) , cof_del )
c-debug-chk;                        chk_sum(11) = chk_sum(11) + cof_del
c-debug-chk;                        chk_ssq(11) = chk_ssq(11) + cof_del**2
c-debug-chk;                     endif
c-debug-chk]
c							      ! 4-X-pt at X=0.2
                     cof_tmp = f_3 * quad(int_hi_3,2,xxx_hi(3),
     $                    co(2,1,1,k,il,kz),co(2,mc,mo,k,il,kz),
     $                    co(4,mc,mo,k,il,kz),
     $                    xx(2),xxx_hi(2),xxx_hi(4))
     $                    + omf_3 * quad(int_hi_3,3,xxx_hi(3),
     $                    co(2,mc,mo,k,il,kz),co(4,mc,mo,k,il,kz),
     $                    co(1,mc,mo_m1,k,il,kz),
     $                    xxx_hi(2),xxx_hi(4),xxx_hi(6))
                     co(3,mc,mo,k,il,kz) = co(3,mc,mo,k,il,kz)
     $                    - cof_tmp
                     int_hi_3 = 1
c-debug-chk[
c-debug-chk;                     cof_o = f_3 * quad(1,2,xxx_hi(3),
c-debug-chk;     $                    co(2,1,1,k,il,kz),co(3,1,1,k,il,kz),
c-debug-chk;     $                    co(4,1,1,k,il,kz),
c-debug-chk;     $                    xx(2),xxx_hi(2),xxx_hi(4))
c-debug-chk;     $                    + omf_3 * quad(1,3,xxx_hi(3),
c-debug-chk;     $                    co(3,1,1,k,il,kz),co(4,1,1,k,il,kz),
c-debug-chk;     $                    co(5,1,1,k,il,kz),
c-debug-chk;     $                    xxx_hi(2),xxx_hi(4),xxx_hi(6))
c-debug-chk;                     cof_del = cof_tmp - cof_o
c-debug-chk]
                  endif
c-debug-chk[
c-debug-chk;                  if ( k+ntdel .gt. 5 ) then
c-debug-chk;                     n_chk(13) = n_chk(13) + 1
c-debug-chk;                     chk_max(13) = max( chk_max(13) , cof_del )
c-debug-chk;                     chk_min(13) = min( chk_min(13) , cof_del )
c-debug-chk;                     chk_sum(13) = chk_sum(13) + cof_del
c-debug-chk;                     chk_ssq(13) = chk_ssq(13) + cof_del**2
c-debug-chk;                  endif
c-debug-chk]
c							! 3-X-pt dlogK at X=0.5
                  cof_tmp = quad(int_hi_z,4,xxx_hi(5),
     $                 co(2,mc,mo,k,il,kz),co(4,mc,mo,k,il,kz),
     $                 co(1,mc,mo_m1,k,il,kz),
     $                 xxx_hi(2),xxx_hi(4),xxx_hi(6))
                  co(5,mc,mo,k,il,kz) = co(5,mc,mo,k,il,kz) - cof_tmp
c-debug-chk[
c-debug-chk;                  if ( k+ntdel .gt. 5 ) then
c-debug-chk;                     cof_o = quad(1,4,xxx_hi(5),
c-debug-chk;     $                    co(3,1,1,k,il,kz),co(4,1,1,k,il,kz),
c-debug-chk;     $                    co(5,1,1,k,il,kz),
c-debug-chk;     $                    xxx_hi(2),xxx_hi(4),xxx_hi(6))
c-debug-chk;                     cof_del = cof_tmp - cof_o
c-debug-chk;                     n_chk(15) = n_chk(15) + 1
c-debug-chk;                     chk_max(15) = max( chk_max(15) , cof_del )
c-debug-chk;                     chk_min(15) = min( chk_min(15) , cof_del )
c-debug-chk;                     chk_sum(15) = chk_sum(15) + cof_del
c-debug-chk;                     chk_ssq(15) = chk_ssq(15) + cof_del**2
c-debug-chk;                  endif
c-debug-chk]
c				     ! 3-X-pt dlogK at X = 0.8, 0.9, 0.95, 1-Z:
                  do ix = 7, mx_use
                     if ( xhi_in(ix) .lt. 1.000001d0 - zat ) then
                        xxx_at = xxx_hi(ix)
                     else
                        xxx_at = log10( 1.d0 - zat + xdel )
                     endif
                     cof_tmp = quad(int_hi_z,ix,xxx_at,
     $                    co(2,mc,mo,k,il,kz),co(4,mc,mo,k,il,kz),
     $                    co(1,mc,mo_m1,k,il,kz),
     $                    xxx_hi(2),xxx_hi(4),xxx_hi(6))
                     co(ix-5,mc,mo_m1,k,il,kz) =
     $                    co(ix-5,mc,mo_m1,k,il,kz) - cof_tmp
c-debug-chk[
c-debug-chk;                     if ( k+ntdel .gt. 5 ) then
c-debug-chk;                        cof_o = quad(1,ix,xxx_at,
c-debug-chk;     $                       co(3,1,1,k,il,kz),co(4,1,1,k,il,kz),
c-debug-chk;     $                       co(5,1,1,k,il,kz),
c-debug-chk;     $                       xxx_hi(2),xxx_hi(4),xxx_hi(6))
c-debug-chk;                        cof_del = cof_tmp - cof_o
c-debug-chk;                        n_chk(ix+10) = n_chk(ix+10) + 1
c-debug-chk;                        chk_max(ix+10) = max( chk_max(ix+10) ,
c-debug-chk;     $                       cof_del )
c-debug-chk;                        chk_min(ix+10) = min( chk_min(ix+10) ,
c-debug-chk;     $                       cof_del )
c-debug-chk;                        chk_sum(ix+10) = chk_sum(ix+10) + cof_del
c-debug-chk;                        chk_ssq(ix+10) = chk_ssq(ix+10)
c-debug-chk;     $                       + cof_del**2
c-debug-chk;                     endif
c-debug-chk]
                  enddo
                  int_hi_z = 1
c				 ! dlogK=0.0 at X = 0.03, 0.1, 0.35, 0.7: these
c					   ! are available in 'Gz???.x??' files
                  co(1,mc,mo,k,il,kz) = 0.0d0
                  co(2,mc,mo,k,il,kz) = 0.0d0
                  co(4,mc,mo,k,il,kz) = 0.0d0
                  co(1,mc,mo_m1,k,il,kz) = 0.0d0
c-debug-chk[
c-debug-chk;                  if ( k+ntdel .gt. 5 ) then
c-debug-chk;                     m = 0
c-debug-chk;                     io = mo
c-debug-chk;                     do ix = 1, mx_use
c-debug-chk;                        m = m + 1
c-debug-chk;                        if ( m .gt. 5 ) then
c-debug-chk;                           m = 1
c-debug-chk;                           io = io - 1
c-debug-chk;                        endif
c-debug-chk;                        cof_del = co(m,mc,io,k,il,kz)
c-debug-chk;                        n_chk(ix) = n_chk(ix) + 1
c-debug-chk;                        chk_max(ix) = max( chk_max(ix) , cof_del )
c-debug-chk;                        chk_min(ix) = min( chk_min(ix) , cof_del )
c-debug-chk;                        chk_sum(ix) = chk_sum(ix) + cof_del
c-debug-chk;                        chk_ssq(ix) = chk_ssq(ix) + cof_del**2
c-debug-chk;                     enddo
c-debug-chk;                  endif
c-debug-chk]
               endif
            enddo
         enddo
c-debug-chk[
c-debug-chk;         write(6,6273) kz, zat, ofebrack, mx_use, iz_hi, kzalow,
c-debug-chk;     $        (n_chk(ix),ix=1,20),
c-debug-chk;     $        (chk_min(ix),ix=1,20), (chk_max(ix),ix=1,20),
c-debug-chk;     $        (chk_sum(ix)/max(n_chk(ix),1),ix=1,20),
c-debug-chk;     $        (sqrt(chk_ssq(ix)/max(n_chk(ix),1)),ix=1,20)
c-debug-chk; 6273    format(' '/' kz =',i3,'  Z =',f10.6,'  [O/Fe] =',f6.3,
c-debug-chk;     $        '  mx_use =',i3,'  iz_hi =',i2,'  kzalow =',i3,
c-debug-chk;     $        ' : X_hi deltas:'/' '/' N',20i10/' min',20f10.6/
c-debug-chk;     $        ' max',20f10.6/' ave',20f10.6/' rms',20f10.6/' ')
c-debug-chk]
c			      ! remember to set X_1 = 0.03 for 'GN93hz' shifts
         xhi_use(1,kz) = 0.03d0
c
      endif
c			! restore old m value, stored at beginning of READEXCO
      m = mstore
c			! and return.
      return
      end
*
******************************************************************************
*
      subroutine cointsmo(xxc,xxo,kz)
*     ===============================
*
*  The purpose of COINTSMO is to interpolate smoothly in C and O abundances.
*
*  This subroutine yields smoother opacities than alternate COINTERP below.
*
*  Note that the quadratic-interpolation function quad has been replaced with
*  the function qchk here; the latter function checks whether two of the
*  interpolation points are nearly coincident (which would magnify the effect
*  of any uncertainties in the tabulated opacities), and uses something more
*  nearly linear if so.  This is sometimes necessary to prevent wildly wrong
*  opacity values for certain Z-values above Z=0.03, and also in some cases to
*  allow linear interpolation along lines where only two opacity values exist.
*  For the special case where C or O is slightly negative (slight depletion in
*  C or O), the function qchk does a linear extrapolation using a combination
*  of the lowest three C or O gridpoints.
*
***-----------------------------------------------------------------------------
      implicit double precision (a-h, o-z)

c
c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=28, nrb=1, nre=28,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime(mx),cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
      common/bb_opal_z/ l1,l2,l3,l4,k1,k2,k3,k4,ip,iq(4),xodp,xcdp,
     $     xxco,cxx,oxx,kzf,kzg,kzh,kzf2
      save /bb_opal_z/
c-debug[
c-debug;      common/outdeb/ ioudeb,oudebl,koudeb
c-debug]
c___
      dimension j2m(mc),j3m(mc),j4m(mc)
c===
      is = 0
c							 ! IF C+O = 0: trivial:
      if ( max( abs(xxc) , abs(xxo) ) .lt. 1.d-6 ) then
         do it = k1,k1+ip
            do ir = l1,l1+iq(it-k1+1)
               opl(it,ir,kz) = co(m,1,1,it,ir,kz)
            enddo
         enddo
         return
      endif
c							 ! ELSE: if C,O not 0:
      i3 = min( indx(int(100.d0*max(xxc,0.d0))+1) + 1 , nc )
      i4 = min( i3+1 , nc )
      i1 = max( i3-2 , 1 )
      i2 = i3-1
      j3 = min( indx(int(100.d0*max(xxo,0.d0))+1) + 1 , no )
      j4 = min( j3+1 , no )
      j1 = max( j3-2 , 1 )
      j2 = j3-1
c
      n2 = i1+1
      n3 = i1+2
      m2 = j1+1
      m3 = j1+2
c			      ! if C > or = O: then j3 < no unless m=5, C=O=0.1
      if ( xxc .ge. xxo ) then
         if ( i4 .gt. n3 ) cfac = (cxx-cx(i2))/(cx(i3)-cx(i2))
c							       ! if O = 0.0:
         if ( abs(xxo) .lt. 1.d-6 ) then
            if ( i4 .le. n3 ) then
               do it = k1,k1+ip
                  do ir = l1,l1+iq(it-k1+1)
                     opl(it,ir,kz) = qchk(is,1,cxx,co(m,i1,1,it,ir,kz),
     $                    co(m,n2,1,it,ir,kz),co(m,n3,1,it,ir,kz),
     $                    cx(i1),cx(n2),min(cx(n3),cxd(1)))
                     is = 1
c-debug[
c-debug;                     if ( ioudeb .gt. 5 .or. .not.
c-debug;     $                    abs(opl(it,ir,kz)) .le. oudebl ) then
c-debug;                        if ( ioudeb .le. 5 .or. ( it .eq. k1 .and.
c-debug;     $                       ir .eq. l1 ) ) write(6,*) ' '
c-debug;                        write(6,9414) m,it,ir,kz,'O=0','C',cxx,'K',
c-debug;     $                       opl(it,ir,kz),0.,(j,1,min(cx(j),cxd(1)),
c-debug;     $                       co(m,j,1,it,ir,kz),j=i1,n3)
c-debug; 9414                   format(' COINTSMO(x',i1,',t',i2.2,',r',i2.2,
c-debug;     $                       ',z',i2.2,')',a3,' ',a1,f11.7,' : ',a1,
c-debug;     $                       ' =',g15.7,' <--(f=',f10.7,
c-debug;     $                       ') kc,ko,CorO,K:',4(i3,i2,2f11.7))
c-debug;                     endif
c-debug]
                  enddo
               enddo
            else
               do it = k1,k1+ip
                  do ir = l1,l1+iq(it-k1+1)
                     opl(it,ir,kz) = (1.d0-cfac) * qchk(is,1,cxx,
     $                    co(m,i1,1,it,ir,kz),co(m,i2,1,it,ir,kz),
     $                    co(m,i3,1,it,ir,kz),cx(i1),cx(i2),cx(i3))
     $                    + cfac * qchk(is,2,cxx,co(m,i2,1,it,ir,kz),
     $                    co(m,i3,1,it,ir,kz),co(m,i4,1,it,ir,kz),
     $                    cx(i2),cx(i3),min(cx(i4),cxd(1)))
                     is = 1
c-debug[
c-debug;                     if ( ioudeb .gt. 5 .or. .not.
c-debug;     $                    abs(opl(it,ir,kz)) .le. oudebl ) then
c-debug;                        if ( ioudeb .le. 5 .or. ( it .eq. k1 .and.
c-debug;     $                       ir .eq. l1 ) ) write(6,*) ' '
c-debug;                        write(6,9414) m,it,ir,kz,'O=0','C',cxx,'K',
c-debug;     $                       opl(it,ir,kz),cfac,(j,1,min(cx(j),cxd(1)),
c-debug;     $                       co(m,j,1,it,ir,kz),j=i1,i4)
c-debug;                     endif
c-debug]
                  enddo
               enddo
            endif
            return
         endif
c						! else if C > O, and O not 0.0
         icomax = -1
         if ( xxco .gt. xcd(1) - 1.d-6 ) icomax = 1
         if ( j1 .lt. j2 ) ofac = (oxx-ox(j2))/(ox(j3)-ox(j2))
         if ( icomax .lt. 0 ) then
            i4 = min( i4 , n(m,j2,kz) )
            cof(i1,1) = cx(i1)
            cof(i2,1) = cx(i2)
            cof(i3,1) = cx(i3)
            cof(i4,1) = cx(i4)
            if ( i3 .ge. n(m,j3,kz) .and. xxo .ge. xod(i3) ) then
               i4 = i3
               icomax = 0
               cof(i4,1) = log10(zzz(kz)+xcdp)
            else if ( i4 .ge. n(m,j3,kz) ) then
               icomax = 0
               cof(i4,1) = log10(zzz(kz)+xcdp)
            endif
            do i = i1,i4
               j2m(i) = m2
               if ( m2 .ge. n(m,i,kz) ) j2m(i) = mo
               j3m(i) = m3
               if ( m3 .ge. n(m,i,kz) ) j3m(i) = mo
               j4m(i) = j4
               if ( j4 .ge. n(m,i,kz) ) j4m(i) = mo
            enddo
         endif
         ihi = i4
         if ( icomax .ge. 0 ) then
            ihi = i4-1
            if ( j4 .lt. no ) then
               j2m(i4) = n(m,j4,kz)
               j4m(i4) = j4
               cof(4,4) = ox(j4)
            else
               j2m(i4) = max( n(m,j4-1,kz) - 2 , 1 )
               j4m(i4) = mo
               cof(4,4) = oxd(j2m(i4))
            endif
            j4n = 0
            if ( xxo .gt. xod(nc-1) + 1.d-6 ) then
               if ( icomax .gt. 0 ) cof(i4,1) = log10(zzz(kz)+xcdp)
               bfac = min( 0.5d0 , ( xxo - xod(nc-1) )
     $              / max( xod(1)-2.d0*xod(nc-1) , 1.d-6 ) )
               ib3 = min( indx(int(100.d0*max(xcdp,0.d0))+1) + 1 , nc )
               ib4 = min( ib3 + 1 , nc )
               ib1 = max( ib3 - 2 , 1 )
               ib2 = ib3-1
               nb2 = ib1+1
               nb3 = ib1+2
               if ( ib4 .gt. nb3 )
     $              afac = (cof(i4,1)-cx(ib2))/(cx(ib3)-cx(ib2))
               if ( ib4 .lt. nc ) then
                  j2n = ib4
                  j4n = mo
                  cof(5,5) = cx(ib4)
               else
                  j4n = max( n(m,ib4-1,kz) - 2 , 1 )
                  j2n = n(m,j4n,kz)
                  cof(5,5) = cxd(j4n)
               endif
            endif
         endif
c
         do it = k1,k1+ip
            do ir = l1,l1+iq(it-k1+1)
               if ( icomax .ge. 0 ) then
                  if ( j4 .le. m3 ) then
                     cof(i4,2) = qchk(is,4,oxx,
     $                    co(m,n(m,j1,kz),j1,it,ir,kz),
     $                    co(m,n(m,m2,kz),m2,it,ir,kz),
     $                    co(m,j2m(i4),j4m(i4),it,ir,kz),
     $                    ox(j1),ox(m2),cof(4,4))
                  else
                     cof(i4,2) = (1.d0-ofac) * qchk(is,4,oxx,
     $                    co(m,n(m,j1,kz),j1,it,ir,kz),
     $                    co(m,n(m,j2,kz),m2,it,ir,kz),
     $                    co(m,n(m,j3,kz),m3,it,ir,kz),
     $                    ox(j1),ox(j2),ox(j3))
     $                    + ofac * qchk(is,8,oxx,
     $                    co(m,n(m,j2,kz),j2,it,ir,kz),
     $                    co(m,n(m,j3,kz),j3,it,ir,kz),
     $                    co(m,j2m(i4),j4m(i4),it,ir,kz),
     $                    ox(j2),ox(j3),cof(4,4))
                  endif
                  if ( j4n .gt. 0 ) then
                     if ( ib4 .le. nb3 ) then
                        cof(i4,2) = (1.d0-bfac) * cof(i4,2)
     $                       + bfac * qchk(is,11,cof(i4,1),
     $                       co(m,ib1,mo,it,ir,kz),
     $                       co(m,nb2,mo,it,ir,kz),
     $                       co(m,j2n,j4n,it,ir,kz),
     $                       cx(ib1),cx(nb2),cof(5,5))
                     else
                        cof(i4,2) = (1.d0-bfac) * cof(i4,2) + bfac
     $                       * ( (1.d0-afac) * qchk(is,11,cof(i4,1),
     $                       co(m,ib1,mo,it,ir,kz),
     $                       co(m,ib2,mo,it,ir,kz),
     $                       co(m,ib3,mo,it,ir,kz),cx(ib1),cx(ib2),
     $                       cx(ib3)) + afac * qchk(is,12,cof(i4,1),
     $                       co(m,ib2,mo,it,ir,kz),
     $                       co(m,ib3,mo,it,ir,kz),
     $                       co(m,j2n,j4n,it,ir,kz),
     $                       cx(ib2),cx(ib3),cof(5,5)) )
                     endif
                  endif
               endif
               if ( icomax .gt. 0 ) then
                  opl(it,ir,kz) = cof(i4,2)
               else
                  iw = 0
                  do i = i1,ihi
                     iw = iw+1
                     if ( j4m(i) .le. j3m(i) ) then
                        cof(i,2) = qchk(is,iw,oxx,co(m,i,j1,it,ir,kz),
     $                       co(m,i,j2m(i),it,ir,kz),
     $                       co(m,i,j3m(i),it,ir,kz),ox(j1),
     $                       min(ox(m2),oxd(i)),min(ox(m3),oxd(i)))
                     else
                        cof(i,2) = (1.d0-ofac) * qchk(is,iw,oxx,
     $                       co(m,i,j1,it,ir,kz),co(m,i,j2,it,ir,kz),
     $                       co(m,i,j3,it,ir,kz),ox(j1),ox(j2),ox(j3))
     $                       + ofac * qchk(is,iw+4,oxx,
     $                       co(m,i,j2,it,ir,kz),
     $                       co(m,i,j3,it,ir,kz),
     $                       co(m,i,j4m(i),it,ir,kz),
     $                       ox(j2),ox(j3),min(ox(j4),oxd(i)))
                     endif
                  enddo
                  if ( i4 .le. n3 ) then
                     opl(it,ir,kz) = qchk(is,9,cxx,cof(i1,2),cof(n2,2),
     $                    cof(n3,2),cof(i1,1),cof(n2,1),cof(n3,1))
                  else
                     opl(it,ir,kz) = (1.d0-cfac) * qchk(is,9,cxx,
     $                    cof(i1,2),cof(i2,2),cof(i3,2),
     $                    cof(i1,1),cof(i2,1),cof(i3,1))
     $                    + cfac * qchk(is,10,cxx,cof(i2,2),cof(i3,2),
     $                    cof(i4,2),cof(i2,1),cof(i3,1),cof(i4,1))
                  endif
               endif
               is = 1
c-debug[
c-debug;               if ( ioudeb .gt. 5 .or.
c-debug;     $              .not. abs(opl(it,ir,kz)) .le. oudebl .or.
c-debug;     $              .not. abs(cof(i4,2)) .le. oudebl .or.
c-debug;     $              ( icomax .le. 0 .and.
c-debug;     $              ( .not. abs(cof(i1,2)) .le. oudebl .or.
c-debug;     $              .not. abs(cof(i2,2)) .le. oudebl .or.
c-debug;     $              .not. abs(cof(i3,2)) .le. oudebl ) ) ) then
c-debug;                  write(6,*) ' '
c-debug;                  if ( icomax .le. 0 ) then
c-debug;                     do i = i1,ihi
c-debug;                        if ( j4m(i) .le. j3m(i) ) then
c-debug;                           write(6,9414) m,it,ir,kz,'C>O','O',oxx,
c-debug;     $                          'K',cof(i,2),0.,i,j1,ox(j1),
c-debug;     $                          co(m,i,j1,it,ir,kz),
c-debug;     $                          i,j2m(i),min(ox(m2),oxd(i)),
c-debug;     $                          co(m,i,j2m(i),it,ir,kz),
c-debug;     $                          i,j3m(i),min(ox(m3),oxd(i)),
c-debug;     $                          co(m,i,j3m(i),it,ir,kz)
c-debug;                        else
c-debug;                           write(6,9414) m,it,ir,kz,'C>O','O',oxx,
c-debug;     $                          'K',cof(i,2),ofac,(i,j,ox(j),
c-debug;     $                          co(m,i,j,it,ir,kz),j=j1,j3),
c-debug;     $                          i,j4m(i),min(ox(j4),oxd(i)),
c-debug;     $                          co(m,i,j4m(i),it,ir,kz)
c-debug;                        endif
c-debug;                     enddo
c-debug;                  endif
c-debug;                  if ( j4 .le. m3 ) ofac = 0.
c-debug;                  if ( icomax .ge. 0 ) write(6,9414) m,it,ir,kz,'C>O',
c-debug;     $                 'o',oxx,'K',cof(i4,2),ofac,(n(m,j,kz),j,ox(j),
c-debug;     $                 co(m,n(m,j,kz),j,it,ir,kz),j=j1,j4-1),j2m(i4),
c-debug;     $                 j4m(i4),cof(4,4),co(m,j2m(i4),j4m(i4),it,ir,kz)
c-debug;                  if ( ib4 .le. nb3 ) afac = 0.
c-debug;                  if ( icomax .ge. 0 .and. j4n .gt. 0 ) write(6,9414)
c-debug;     $                 m,it,ir,kz,'C>O','c',cof(i4,1),'b',bfac,afac,
c-debug;     $                 (j,mo,cx(j),co(m,j,mo,it,ir,kz),j=ib1,ib4-1),
c-debug;     $                 j2n,j4n,cof(5,5),co(m,j2n,j4n,it,ir,kz)
c-debug;                  if ( i4 .le. n3 ) cfac = 0.
c-debug;                  if ( icomax .le. 0 ) write(6,9415) m,it,ir,kz,
c-debug;     $                 'C>O','C',cxx,opl(it,ir,kz),cfac,
c-debug;     $                 (' ',j,cof(j,1),cof(j,2),j=i1,i4)
c-debug; 9415             format(' COINTSMO(x',i1,','t,i2.2,',r',i2.2,
c-debug;     $                 ',z',i2.2,')',a3,' ',a1,f11.7,
c-debug;     $                 ' : K =',g15.7,' <--(f=',f10.7,
c-debug;     $                 ') kc,ko,CorO,K:',4(a1,' (',i1,')',2f11.7))
c-debug;               endif
c-debug]
            enddo
         enddo
c					! else if C < O: then i3 < nc
      else
         if ( j4 .gt. m3 ) ofac = (oxx-ox(j2))/(ox(j3)-ox(j2))
c							       ! if C = 0.0:
         if ( abs(xxc) .lt. 1.d-6 ) then
            if ( j4 .le. m3 ) then
               j3m(1) = m3
               if ( m3 .ge. no ) j3m(1) = mo
               do it = k1,k1+ip
                  do ir = l1,l1+iq(it-k1+1)
                     opl(it,ir,kz) = qchk(is,1,oxx,co(m,1,j1,it,ir,kz),
     $                    co(m,1,m2,it,ir,kz),co(m,1,j3m(1),it,ir,kz),
     $                    ox(j1),ox(m2),min(ox(m3),oxd(1)))
                     is = 1
c-debug[
c-debug;                     if ( ioudeb .gt. 5 .or. .not.
c-debug;     $                    abs(opl(it,ir,kz)) .le. oudebl ) then
c-debug;                        if ( ioudeb .le. 5 .or. ( it .eq. k1 .and.
c-debug;     $                       ir .eq. l1 ) ) write(6,*) ' '
c-debug;                        write(6,9414) m,it,ir,kz,'C=0','O',oxx,'K',
c-debug;     $                       opl(it,ir,kz),0.,(1,j,ox(j),
c-debug;     $                       co(m,1,j,it,ir,kz),j=j1,m2),1,j3m(1),
c-debug;     $                       co(m,1,j3m(1),it,ir,kz),min(ox(j3),oxd(1))
c-debug;                     endif
c-debug]
                  enddo
               enddo
            else
               j4m(1) = j4
               if ( j4 .ge. no ) j4m(1) = mo
               do it = k1,k1+ip
                  do ir = l1,l1+iq(it-k1+1)
                     opl(it,ir,kz) = (1.d0-ofac) * qchk(is,1,oxx,
     $                    co(m,1,j1,it,ir,kz),co(m,1,j2,it,ir,kz),
     $                    co(m,1,j3,it,ir,kz),ox(j1),ox(j2),ox(j3))
     $                    + ofac * qchk(is,2,oxx,co(m,1,j2,it,ir,kz),
     $                    co(m,1,j3,it,ir,kz),co(m,1,j4m(1),it,ir,kz),
     $                    ox(j2),ox(j3),min(ox(j4),oxd(1)))
                     is = 1
c-debug[
c-debug;                     if ( ioudeb .gt. 5 .or. .not.
c-debug;     $                    abs(opl(it,ir,kz)) .le. oudebl ) then
c-debug;                        if ( ioudeb .le. 5 .or. ( it .eq. k1 .and.
c-debug;     $                       ir .eq. l1 ) ) write(6,*) ' '
c-debug;                        write(6,9414) m,it,ir,kz,'C=0','O',oxx,'K',
c-debug;     $                       opl(it,ir,kz),ofac,(1,j,ox(j),
c-debug;     $                       co(m,1,j,it,ir,kz),j=j1,j3),
c-debug;     $                       1,j4m(1),co(m,1,j4m(1),it,ir,kz),
c-debug;     $                       min(ox(j4),oxd(1))
c-debug;                     endif
c-debug]
                  enddo
               enddo
            endif
            return
         endif
c						! else if O > C, and C not 0.0:
         icomax = -1
         if ( xxco .gt. xcd(1)-1.d-6 ) icomax = 1
         if ( i1 .lt. i2 ) cfac = (cxx-cx(i2))/(cx(i3)-cx(i2))
         if ( icomax .lt. 0 ) then
            j4 = min( j4 , n(m,i2,kz) )
            cof(j1,1) = ox(j1)
            cof(j2,1) = ox(j2)
            cof(j3,1) = ox(j3)
            cof(j4,1) = ox(j4)
            if ( j3 .ge. n(m,i3,kz) .and. xxc .ge. xcd(j3) ) then
               j4 = j3
               icomax = 0
               cof(j4,1) = log10(zzz(kz)+xodp)
            else if ( j4 .ge. n(m,i3,kz) ) then
               icomax = 0
               cof(j4,1) = log10(zzz(kz)+xodp)
            endif
         endif
         ihi = j4
         if ( icomax .ge. 0 ) then
            ihi = j4-1
            if ( i4 .lt. nc ) then
               j2m(4) = i4
               j4m(4) = mo
               cof(4,4) = cx(i4)
            else
               j4m(4) = max( n(m,i4-1,kz) - 2 , 1 )
               j2m(4) = n(m,j4m(4),kz)
               cof(4,4) = cxd(j4m(4))
            endif
            j4n = 0
            if ( xxc .gt. xcd(no-1) + 1.d-6 ) then
               if ( icomax .gt. 0 ) cof(j4,1) = log10(zzz(kz)+xodp)
               bfac = min((xxc-xcd(no-1))/max(xcd(1)-2.d0*xcd(no-1),
     $              1.d-6),0.5d0)
               jb3 = min( indx(int(100.d0*max(xodp,0.d0))+1) + 1 , no )
               jb4 = min( jb3+1 , no )
               jb1 = max( jb3-2 , 1 )
               jb2 = jb3-1
               mb2 = jb1+1
               mb3 = jb1+2
               if ( jb4 .gt. mb3 )
     $              afac = (cof(j4,1)-ox(jb2))/(ox(jb3)-ox(jb2))
               if ( jb4 .lt. no ) then
                  j2n = n(m,jb4,kz)
                  j4n = jb4
                  cof(5,5) = ox(jb4)
               else
                  j2n = max( n(m,jb4-1,kz) - 2 , 1 )
                  j4n = mo
                  cof(5,5) = oxd(j2n)
               endif
            endif
         endif
c
         do it = k1,k1+ip
            do ir = l1,l1+iq(it-k1+1)
               if ( icomax .ge. 0 ) then
                  if ( i4 .le. n3 ) then
                     cof(j4,2) = qchk(is,4,cxx,
     $                    co(m,i1,mo,it,ir,kz),co(m,n2,mo,it,ir,kz),
     $                    co(m,j2m(4),j4m(4),it,ir,kz),
     $                    cx(i1),cx(n2),cof(4,4))
                  else
                     cof(j4,2) = (1.d0-cfac) * qchk(is,4,cxx,
     $                    co(m,i1,mo,it,ir,kz),co(m,i2,mo,it,ir,kz),
     $                    co(m,i3,mo,it,ir,kz),cx(i1),cx(i2),cx(i3))
     $                    + cfac * qchk(is,8,cxx,
     $                    co(m,i2,mo,it,ir,kz),co(m,i3,mo,it,ir,kz),
     $                    co(m,j2m(4),j4m(4),it,ir,kz),
     $                    cx(i2),cx(i3),cof(4,4))
                  endif
                  if ( j4n .gt. 0 ) then
                     if ( jb4 .le. mb3 ) then
                        cof(j4,2) = (1.d0-bfac) * cof(j4,2)
     $                       + bfac * qchk(is,11,cof(j4,1),
     $                       co(m,n(m,jb1,kz),jb1,it,ir,kz),
     $                       co(m,n(m,mb2,kz),mb2,it,ir,kz),
     $                       co(m,j2n,j4n,it,ir,kz),
     $                       ox(jb1),ox(mb2),cof(5,5))
                     else
                        cof(j4,2) = (1.d0-bfac) * cof(j4,2) + bfac
     $                       * ( (1.d0-afac) * qchk(is,11,cof(j4,1),
     $                       co(m,n(m,jb1,kz),jb1,it,ir,kz),
     $                       co(m,n(m,jb2,kz),mb2,it,ir,kz),
     $                       co(m,n(m,jb3,kz),mb3,it,ir,kz),
     $                       ox(jb1),ox(jb2),ox(jb3))
     $                       + afac * qchk(is,12,cof(j4,1),
     $                       co(m,n(m,jb2,kz),jb2,it,ir,kz),
     $                       co(m,n(m,jb3,kz),jb3,it,ir,kz),
     $                       co(m,j2n,j4n,it,ir,kz),
     $                       ox(jb2),ox(jb3),cof(5,5)) )
                     endif
                  endif
               endif
               if ( icomax .gt. 0 ) then
                  opl(it,ir,kz) = cof(j4,2)
               else
                  iw = 0
                  do i = j1,ihi
                     iw = iw+1
                     if ( i4 .le. n3 .or. i4 .gt. n(m,i,kz) ) then
                        cof(i,2) = qchk(is,iw,cxx,co(m,i1,i,it,ir,kz),
     $                       co(m,n2,i,it,ir,kz),
     $                       co(m,min(n3,n(m,i,kz)),i,it,ir,kz),
     $                       cx(i1),min(cx(n2),cxd(i)),
     $                       min(cx(n3),cxd(i)))
                     else
                        cof(i,2) = (1.d0-cfac) * qchk(is,iw,cxx,
     $                       co(m,i1,i,it,ir,kz),co(m,i2,i,it,ir,kz),
     $                       co(m,i3,i,it,ir,kz),cx(i1),cx(i2),cx(i3))
     $                       + cfac * qchk(is,iw+4,cxx,
     $                       co(m,i2,i,it,ir,kz),co(m,i3,i,it,ir,kz),
     $                       co(m,i4,i,it,ir,kz),
     $                       cx(i2),cx(i3),min(cx(i4),cxd(i)))
                     endif
                  enddo
                  if ( j4 .le. m3 ) then
                     opl(it,ir,kz) = qchk(is,9,oxx,cof(j1,2),cof(m2,2),
     $                    cof(m3,2),cof(j1,1),cof(m2,1),cof(m3,1))
                  else
                     opl(it,ir,kz) = (1.d0-ofac) * qchk(is,9,oxx,
     $                    cof(j1,2),cof(j2,2),cof(j3,2),
     $                    cof(j1,1),cof(j2,1),cof(j3,1))
     $                    + ofac * qchk(is,10,oxx,cof(j2,2),cof(j3,2),
     $                    cof(j4,2),cof(j2,1),cof(j3,1),cof(j4,1))
                  endif
               endif
               is = 1
c-debug[
c-debug;               if ( ioudeb .gt. 5 .or.
c-debug;     $              .not. abs(opl(it,ir,kz)) .le. oudebl .or.
c-debug;     $              .not. abs(cof(j4,2)) .le. oudebl .or.
c-debug;     $              ( icomax .le. 0 .and.
c-debug;     $              ( .not. abs(cof(j1,2)) .le. oudebl .or.
c-debug;     $              .not. abs(cof(j2,2)) .le. oudebl .or.
c-debug;     $              .not. abs(cof(j3,2)) .le. oudebl ) ) ) then
c-debug;                  write(6,*) ' '
c-debug;                  if ( icomax .le. 0 ) then
c-debug;                     do i = j1,ihi
c-debug;                        if ( i4 .le. n3 .or. i4 .gt. n(m,i,kz) ) then
c-debug;                           write(6,9414) m,it,ir,kz,'C<O','C',cxx,'K',
c-debug;     $                          cof(i,2),0.,i1,i,cx(i1),
c-debug;     $                          co(m,i1,i,it,ir,kz),n2,i,
c-debug;     $                          min(cx(n2),cxd(i)),co(m,n2,i,it,ir,kz),
c-debug;     $                          min(n3,n(m,i,kz)),i,min(cx(n3),cxd(i)),
c-debug;     $                          co(m,min(n3,n(m,i,kz)),i,it,ir,kz)
c-debug;                        else
c-debug;                           write(6,9414) m,it,ir,kz,'C<O','C',cxx,'K',
c-debug;     $                          cof(i,2),cfac,(j,i,min(cx(j),cxd(i)),
c-debug;     $                          co(m,j,i,it,ir,kz),j=i1,i4)
c-debug;                        endif
c-debug;                     enddo
c-debug;                  endif
c-debug;                  if ( i4 .le. n3 ) cfac = 0.
c-debug;                  if ( icomax .ge. 0 ) write(6,9414) m,it,ir,kz,
c-debug;     $                 'C<O','c',cxx,'K',cof(j4,2),cfac,(j,mo,cx(j),
c-debug;     $                 co(m,j,mo,it,ir,kz),j=i1,i4-1),j2m(4),
c-debug;     $                 j4m(4),cof(4,4),co(m,j2m(4),j4m(4),it,ir,kz)
c-debug;                  if ( jb4 .le. mb3 ) afac = 0.
c-debug;                  if ( icomax .ge. 0 .and. j4n .gt. 0 ) write(6,9414)
c-debug;     $                 m,it,ir,kz,'C<O','o',cof(j4,1),'b',bfac,afac,
c-debug;     $                 (n(m,j,kz),j,ox(j),co(m,n(m,j,kz),j,it,ir,kz),
c-debug;     $                 j=jb1,jb4-1),j2n,j4n,cof(5,5),
c-debug;     $                 co(m,j2n,j4n,it,ir,kz)
c-debug;                  if ( j4 .le. m3 ) ofac = 0.
c-debug;                  if ( icomax .le. 0 ) write(6,9415) m,it,ir,kz,
c-debug;     $                 'C<O','O',oxx,opl(it,ir,kz),ofac,
c-debug;     $                 (' ',j,cof(j,1),cof(j,2),j=j1,j4)
c-debug;               endif
c-debug]
            enddo
         enddo
      endif
c
      return
      end


*
***-------------------------------------------------------------------------

      subroutine cointerp(xxc,xxo,kz)
*     ===============================
*
*     Interpolate in C and O abundances.
*
*  This subroutine yields less-smooth opacities than alternate
*  subroutine COINTSMO above.
*
*  Note that the quadratic-interpolation function quad has been
*  replaced with the function qchk here; the latter function checks
*  whether two of the interpolation points are nearly coincident
*  (which would magnify the effect of any uncertainties in the
*  tabulated opacities), and uses something more nearly linear
*  if so.  This is sometimes necessary to prevent wildly wrong
*  opacity values for certain Z-values above Z=0.03.
*
*
***-------------------------------------------------------------------------

      implicit double precision (a-h, o-z)


c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=28, nrb=1, nre=28,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime(mx),cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
      common/bb_opal_z/ l1,l2,l3,l4,k1,k2,k3,k4,ip,iq(4),xodp,xcdp,
     $     xxco,cxx,oxx,kzf,kzg,kzh,kzf2
      save /bb_opal_z/

c-debug[
c-debug;      common/outdeb/ ioudeb,oudebl,koudeb
c-debug]
c
c===
c
c-debug[
c-debug; 9413 format(' '/' COINTERP: m=',i2,' kz=',i3,' xxc',f10.7,
c-debug;     $     ' xxo',f10.7,' xxc+xxo',f10.7,' xxco',f10.7,
c-debug;     $     ' COmax',f10.7,' nc',i2,' n(m,1:nc)',8i2)
c-debug]
c
      is = 0
c							! if C+O = 0:
      if ( max( abs(xxc) , abs(xxo) ) .lt. 1.d-6 ) then
         do it = k1,k1+ip
            do ir = l1,l1+iq(it-k1+1)
               opl(it,ir,kz) = co(m,1,1,it,ir,kz)
            enddo
         enddo
c
         return
c         ! include boundaries (later, fix any possible division by 0)
c
      else if ( xxc .gt. xcd(3)-1.d-6 ) then
c__________
         oxdp = log10(zzz(kz)+xodp)
c					  ! interpolation in region c1:
c
c     ! include boundaries (qchk & fac fix any possible division by 0)
c
         if ( xxc .gt. xcd(2)-1.d-6 ) then
c				     ! handle possibility that xodp = 0
            fac = max(min((oxx-ox(1))/max(oxdp-ox(1),1.d-6),1.d0),0.d0)
c
            do it = k1,k1+ip
               do ir = l1,l1+iq(it-k1+1)
                  b(1) = qchk(is,1,cxx,co(m,nc-2,1,it,ir,kz),
     $                 co(m,nc-1,1,it,ir,kz),co(m,nc,1,it,ir,kz),
     $                 cx(nc-2),cx(nc-1),cx(nc))
                  b(2) = qchk(is,2,cxx,co(m,nc,1,it,ir,kz),
     $                 co(m,n(m,2,kz),2,it,ir,kz),
     $                 co(m,n(m,3,kz),3,it,ir,kz),
     $                 cxd(1),cxd(2),cxd(3))
c				     ! handle possibility that xodp = 0
                  opl(it,ir,kz) = b(1)+(b(2)-b(1))*fac
                  is = 1
c-debug[
c-debug;                  if ( .not. abs(opl(it,ir,kz)) .le. oudebl .or.
c-debug;     $                 .not. abs(b(1)) .le. oudebl .or.
c-debug;     $                 .not. abs(b(2)) .le. oudebl .or.
c-debug;     $                 ioudeb .gt. 5 ) then
c-debug;                     write(6,9413) m,kz,xxc,xxo,xxc+xxo,xxco,
c-debug;     $                    xc(nc),nc,(n(m,j,kz),j=1,nc)
c-debug;                     koudeb = koudeb+1
c-debug;                     write(6,9414) 1,'co(m,',nc-2,1,it,ir,
c-debug;     $                    co(m,nc-2,1,it,ir,kz),'co(m,',nc-1,
c-debug;     $                    1,it,ir,co(m,nc-1,1,it,ir,kz),
c-debug;     $                    'diag(',m,1,it,ir,
c-debug;     $                    co(m,nc,1,it,ir,kz),'cxx',cxx,
c-debug;     $                    ' cx(',nc-2,cx(nc-2),' cx(',nc-1,
c-debug;     $                    cx(nc-1),' cx(',nc,cx(nc)
c-debug; 9414                format(' b',i1,3(' ',a5,i1,',',i1,',',
c-debug;     $                    i2,',',i2,')',f12.7),' ',a3,f12.7,
c-debug;     $                    3(' ',a4,i1,')',f12.7))
c-debug;                     write(6,9414) 2,'diag(',m,1,it,ir,
c-debug;     $                    co(m,nc,1,it,ir,kz),'diag(',m,2,it,
c-debug;     $                    ir,co(m,n(m,2,kz),2,it,ir,kz),
c-debug;     $                    'diag(',m,3,it,ir,
c-debug;     $                    co(m,n(m,3,kz),3,it,ir,kz),'cxx',
c-debug;     $                    cxx,'cxd(',1,cxd(1),'cxd(',2,
c-debug;     $                    cxd(2),'cxd(',3,cxd(3)
c-debug;                     write(6,9415) b(1),b(2),0.,'oxx',oxx,
c-debug;     $                    ' ox(',1,ox(1),'oxdp',0,oxdp,'----',
c-debug;     $                    0,0.,it,ir,opl(it,ir,kz),'c1  '
c-debug; 9415                format(' --> b:',3f12.7,' ',a3,f12.7,
c-debug;     $                    3(' ',a4,i1,')',f12.7),' --> opl(',
c-debug;     $                    i2,',',i2,') =',f12.7,'  region ',
c-debug;     $                    a4)
c-debug;                  endif
c-debug]
               enddo
            enddo
c					  ! interpolation in region c2:
         else
            do it = k1,k1+ip
               do ir = l1,l1+iq(it-k1+1)
                  b(1) = qchk(is,1,cxx,co(m,nc-2,1,it,ir,kz),
     $                 co(m,nc-1,1,it,ir,kz),co(m,nc,1,it,ir,kz),
     $                 cx(nc-2),cx(nc-1),cx(nc))
                  b(2) = qchk(is,2,cxx,co(m,n(m,2,kz)-2,2,it,ir,kz),
     $                 co(m,n(m,2,kz)-1,2,it,ir,kz),
     $                 co(m,n(m,2,kz),2,it,ir,kz),
     $                 cx(n(m,2,kz)-2),cx(n(m,2,kz)-1),cxd(2))
                  b(3) = qchk(is,3,cxx,co(m,nc,1,it,ir,kz),
     $                 co(m,n(m,2,kz),2,it,ir,kz),
     $                 co(m,n(m,3,kz),3,it,ir,kz),
     $                 cxd(1),cxd(2),cxd(3))
                  opl(it,ir,kz) = qchk(is,4,oxx,b(1),b(2),b(3),
     $                 ox(1),ox(2),oxdp)
                  is = 1
c-debug[
c-debug;                  if ( .not. abs(opl(it,ir,kz)) .le. oudebl .or.
c-debug;     $                 .not. abs(b(1)) .le. oudebl .or.
c-debug;     $                 .not. abs(b(2)) .le. oudebl .or.
c-debug;     $                 .not. abs(b(3)) .le. oudebl .or.
c-debug;     $                 ioudeb .gt. 5 ) then
c-debug;                     write(6,9413) m,kz,xxc,xxo,xxc+xxo,xxco,
c-debug;     $                    xc(nc),nc,(n(m,j,kz),j=1,nc)
c-debug;                     koudeb = koudeb+1
c-debug;                     write(6,9414) 1,'co(m,',nc-2,1,it,ir,
c-debug;     $                    co(m,nc-2,1,it,ir,kz),'co(m,',nc-1,
c-debug;     $                    1,it,ir,co(m,nc-1,1,it,ir,kz),
c-debug;     $                    'diag(',m,1,it,ir,
c-debug;     $                    co(m,nc,1,it,ir,kz),'cxx',cxx,
c-debug;     $                    ' cx(',nc-2,cx(nc-2),' cx(',nc-1,
c-debug;     $                    cx(nc-1),' cx(',nc,cx(nc)
c-debug;                     write(6,9414) 2,'co(m,',n(m,2,kz)-2,2,it,
c-debug;     $                    ir,co(m,n(m,2,kz)-2,2,it,ir,kz),
c-debug;     $                    'co(m,',n(m,2,kz)-1,2,it,ir,
c-debug;     $                    co(m,n(m,2,kz)-1,2,it,ir,kz),
c-debug;     $                    'diag(',m,2,it,ir,
c-debug;     $                    co(m,n(m,2,kz),2,it,ir,kz),
c-debug;     $                    'cxx',cxx,
c-debug;     $                    ' cx(',n(m,2,kz)-2,cx(n(m,2,kz)-2),
c-debug;     $                    ' cx(',n(m,2,kz)-1,cx(n(m,2,kz)-1),
c-debug;     $                    'cxd(',2,cxd(2)
c-debug;                     write(6,9414) 3,'diag(',m,1,it,ir,
c-debug;     $                    co(m,nc,1,it,ir,kz),'diag(',m,2,it,
c-debug;     $                    ir,co(m,n(m,2,kz),2,it,ir,kz),
c-debug;     $                    'diag(',m,3,it,ir,
c-debug;     $                    co(m,n(m,3,kz),3,it,ir,kz),'cxx',
c-debug;     $                    cxx,'cxd(',1,cxd(1),'cxd(',2,
c-debug;     $                    cxd(2),'cxd(',3,cxd(3)
c-debug;                     write(6,9415) b(1),b(2),b(3),'oxx',oxx,
c-debug;     $                    ' ox(',1,ox(1),' ox(',2,ox(2),
c-debug;     $                    'oxdp',0,oxdp,it,ir,opl(it,ir,kz),
c-debug;     $                    'c2  '
c-debug;                  endif
c-debug]
               enddo
            enddo
         endif
c
         return
c
      else if ( nc .ge. 5 ) then
c__________			   ! interpolation in regions c3 to c6:
         do i = 4,nc-1
c
c	   ! do not go beyond middle (where c3-c6 overlaps o3-o6), and
c	   ! include boundaries (qchk fixes any possible division by 0)
c
            if ( xxc .gt. xcd(i)-1.d-6 .and. xxo .gt. xo(i-1)-1.d-6
     $           .and. xcd(i-1) .gt. xc(i-1) ) then
c
               oxdp = log10(zzz(kz)+xodp)
               m1 = i-1
               m2 = i-2
c
               do it = k1,k1+ip
                  do ir = l1,l1+iq(it-k1+1)
                     b(1) = qchk(is,1,cxx,
     $                    co(m,n(m,m2,kz)-2,m2,it,ir,kz),
     $                    co(m,n(m,m2,kz)-1,m2,it,ir,kz),
     $                    co(m,n(m,m2,kz),m2,it,ir,kz),
     $                    cx(n(m,m2,kz)-2),cx(n(m,m2,kz)-1),cxd(m2))
                     b(2) = qchk(is,2,cxx,
     $                    co(m,n(m,m1,kz)-2,m1,it,ir,kz),
     $                    co(m,n(m,m1,kz)-1,m1,it,ir,kz),
     $                    co(m,n(m,m1,kz),m1,it,ir,kz),
     $                    cx(n(m,m1,kz)-2),cx(n(m,m1,kz)-1),cxd(m1))
                     b(3) = qchk(is,3,cxx,
     $                    co(m,n(m,m2,kz),m2,it,ir,kz),
     $                    co(m,n(m,m1,kz),m1,it,ir,kz),
     $                    co(m,n(m,i,kz),i,it,ir,kz),
     $                    cxd(m2),cxd(m1),cxd(i))
                     opl(it,ir,kz) = qchk(is,4,oxx,b(1),b(2),b(3),
     $                    ox(m2),ox(m1),oxdp)
                     is = 1
c-debug[
c-debug;                     if ( .not. abs(opl(it,ir,kz)) .le. oudebl .or.
c-debug;     $                    .not. abs(b(1)) .le. oudebl .or.
c-debug;     $                    .not. abs(b(2)) .le. oudebl .or.
c-debug;     $                    .not. abs(b(3)) .le. oudebl .or.
c-debug;     $                    ioudeb .gt. 5 ) then
c-debug;                        write(6,9413) m,kz,xxc,xxo,xxc+xxo,xxco,
c-debug;     $                       xc(nc),nc,(n(m,j,kz),j=1,nc)
c-debug;                        koudeb = koudeb+1
c-debug;                        write(6,9414) 1,'co(m,',n(m,m2,kz)-2,
c-debug;     $                       m2,it,ir,
c-debug;     $                       co(m,n(m,m2,kz)-2,m2,it,ir,kz),
c-debug;     $                       'co(m,',n(m,m2,kz)-1,m2,it,ir,
c-debug;     $                       co(m,n(m,m2,kz)-1,m2,it,ir,kz),
c-debug;     $                       'diag(',m,m2,it,ir,
c-debug;     $                       co(m,n(m,m2,kz),m2,it,ir,kz),
c-debug;     $                       'cxx',cxx,' cx(',n(m,m2,kz)-2,
c-debug;     $                       cx(n(m,m2,kz)-2),' cx(',
c-debug;     $                       n(m,m2,kz)-1,cx(n(m,m2,kz)-1),
c-debug;     $                       'cxd(',m2,cxd(m2)
c-debug;                        write(6,9414) 2,'co(m,',n(m,m1,kz)-2,
c-debug;     $                       m1,it,ir,
c-debug;     $                       co(m,n(m,m1,kz)-2,m1,it,ir,kz),
c-debug;     $                       'co(m,',n(m,m1,kz)-1,m1,it,ir,
c-debug;     $                       co(m,n(m,m1,kz)-1,m1,it,ir,kz),
c-debug;     $                       'diag(',m,m1,it,ir,
c-debug;     $                       co(m,n(m,m1,kz),m1,it,ir,kz),
c-debug;     $                       'cxx',cxx,' cx(',n(m,m1,kz)-2,
c-debug;     $                       cx(n(m,m1,kz)-2),' cx(',
c-debug;     $                       n(m,m1,kz)-1,cx(n(m,m1,kz)-1),
c-debug;     $                       'cxd(',m1,cxd(m1)
c-debug;                        write(6,9414) 3,'diag(',m,m2,it,ir,
c-debug;     $                       co(m,n(m,m2,kz),m2,it,ir,kz),
c-debug;     $                       'diag(',m,m1,it,ir,
c-debug;     $                       co(m,n(m,m1,kz),m1,it,ir,kz),
c-debug;     $                       'diag(',m,i,it,ir,
c-debug;     $                       co(m,n(m,i,kz),i,it,ir,kz),
c-debug;     $                       'cxx',cxx,'cxd(',m2,cxd(m2),
c-debug;     $                       'cxd(',m1,cxd(m1),'cxd(',i,cxd(i)
c-debug;                        write(6,9415) b(1),b(2),b(3),'oxx',
c-debug;     $                       oxx,' ox(',i-2,ox(i-2),' ox(',
c-debug;     $                       i-1,ox(i-1),'oxdp',0,oxdp,it,ir,
c-debug;     $                       opl(it,ir,kz),'c3-6'
c-debug;                     endif
c-debug]
                  enddo
               enddo
c
               return
c
            endif
c
         enddo
      endif
c
c	   ! include boundaries (later, fix any possible division by 0)
c
      if ( xxo .gt. xod(3)-1.d-6 ) then
c__________
         cxdp = log10(zzz(kz)+xcdp)
c					  ! interpolation in region o1:
c
c     ! include boundaries (qchk & fac fix any possible division by 0)
c
         if ( xxo .gt. xod(2)-1.d-6 ) then
c				     ! handle possibility that xcdp = 0
            fac = max(min((cxx-cx(1))/max(cxdp-cx(1),1.d-6),1.d0),0.d0)
c
            do it = k1,k1+ip
               do ir = l1,l1+iq(it-k1+1)
                  b(1) = qchk(is,1,oxx,co(m,1,no-2,it,ir,kz),
     $                 co(m,1,no-1,it,ir,kz),co(m,1,mo,it,ir,kz),
     $                 ox(no-2),ox(no-1),ox(no))
                  b(2) = qchk(is,2,oxx,co(m,1,mo,it,ir,kz),
     $                 co(m,2,mo,it,ir,kz),co(m,3,mo,it,ir,kz),
     $                 oxd(1),oxd(2),oxd(3))
c				     ! handle possibility that xcdp = 0
                  opl(it,ir,kz) = b(1)+(b(2)-b(1))*fac
                  is = 1
c-debug[
c-debug;                  if ( .not. abs(opl(it,ir,kz)) .le. oudebl .or.
c-debug;     $                 .not. abs(b(1)) .le. oudebl .or.
c-debug;     $                 .not. abs(b(2)) .le. oudebl .or.
c-debug;     $                 ioudeb .gt. 5 ) then
c-debug;                     write(6,9413) m,kz,xxc,xxo,xxc+xxo,xxco,
c-debug;     $                    xc(nc),nc,(n(m,j,kz),j=1,nc)
c-debug;                     koudeb = koudeb+1
c-debug;                     write(6,9414) 1,'co(m,',1,no-2,it,ir,
c-debug;     $                    co(m,1,no-2,it,ir,kz),'co(m,',1,
c-debug;     $                    no-1,it,ir,co(m,1,no-1,it,ir,kz),
c-debug;     $                    'digo(',m,no-1,it,ir,
c-debug;     $                    co(m,1,mo,it,ir,kz),'oxx',
c-debug;     $                    oxx,' ox(',no-2,ox(no-2),' ox(',
c-debug;     $                    no-1,ox(no-1),' ox(',no,ox(no)
c-debug;                     write(6,9414) 2,'digo(',m,no-1,it,ir,
c-debug;     $                    co(m,1,mo,it,ir,kz),'digo(',m,no-2,
c-debug;     $                    it,ir,co(m,2,mo,it,ir,kz),'digo(',m,
c-debug;     $                    no-3,it,ir,co(m,3,mo,it,ir,kz),
c-debug;     $                    'oxx',oxx,'oxd(',1,oxd(1),'oxd(',2,
c-debug;     $                    oxd(2),'oxd(',3,oxd(3)
c-debug;                     write(6,9415) b(1),b(2),0.,'cxx',cxx,
c-debug;     $                    ' cx(',1,cx(1),'cxdp',0,cxdp,'----',
c-debug;     $                    0,0.,it,ir,opl(it,ir,kz),'o1  '
c-debug;                  endif
c-debug]
               enddo
            enddo
c					  ! interpolation in region o2:
         else
            do it = k1,k1+ip
               do ir = l1,l1+iq(it-k1+1)
                  b(1) = qchk(is,1,oxx,co(m,1,no-2,it,ir,kz),
     $                 co(m,1,no-1,it,ir,kz),co(m,1,mo,it,ir,kz),
     $                 ox(no-2),ox(no-1),ox(no))
                  b(2) = qchk(is,2,oxx,co(m,2,n(m,2,kz)-2,it,ir,kz),
     $                 co(m,2,n(m,2,kz)-1,it,ir,kz),
     $                 co(m,2,mo,it,ir,kz),
     $                 ox(n(m,2,kz)-2),ox(n(m,2,kz)-1),oxd(2))
                  b(3) = qchk(is,3,oxx,co(m,1,mo,it,ir,kz),
     $                 co(m,2,mo,it,ir,kz),co(m,3,mo,it,ir,kz),
     $                 oxd(1),oxd(2),oxd(3))
                  opl(it,ir,kz) = qchk(is,4,cxx,b(1),b(2),b(3),
     $                 cx(1),cx(2),cxdp)
                  is = 1
c-debug[
c-debug;                  if ( .not. abs(opl(it,ir,kz)) .le. oudebl .or.
c-debug;     $                 .not. abs(b(1)) .le. oudebl .or.
c-debug;     $                 .not. abs(b(2)) .le. oudebl .or.
c-debug;     $                 .not. abs(b(3)) .le. oudebl .or.
c-debug;     $                 ioudeb .gt. 5 ) then
c-debug;                     write(6,9413) m,kz,xxc,xxo,xxc+xxo,xxco,
c-debug;     $                    xc(nc),nc,(n(m,j,kz),j=1,nc)
c-debug;                     koudeb = koudeb+1
c-debug;                     write(6,9414) 1,'co(m,',1,no-2,it,ir,
c-debug;     $                    co(m,1,no-2,it,ir,kz),'co(m,',1,
c-debug;     $                    no-1,it,ir,co(m,1,no-1,it,ir,kz),
c-debug;     $                    'digo(',m,no-1,it,ir,
c-debug;     $                    co(m,1,mo,it,ir,kz),'oxx',
c-debug;     $                    oxx,' ox(',no-2,ox(no-2),' ox(',
c-debug;     $                    no-1,ox(no-1),' ox(',no,ox(no)
c-debug;                     write(6,9414) 2,'co(m,',2,n(m,2,kz)-2,it,
c-debug;     $                    ir,co(m,2,n(m,2,kz)-2,it,ir,kz),
c-debug;     $                    'co(m,',2,n(m,2,kz)-1,it,ir,
c-debug;     $                    co(m,2,n(m,2,kz)-1,it,ir,kz),
c-debug;     $                    'digo(',m,no-2,it,ir,
c-debug;     $                    co(m,2,mo,it,ir,kz),'oxx',oxx,
c-debug;     $                    ' ox(',n(m,2,kz)-2,ox(n(m,2,kz)-2),
c-debug;     $                    ' ox(',n(m,2,kz)-1,ox(n(m,2,kz)-1),
c-debug;     $                    'oxd(',2,oxd(2)
c-debug;                     write(6,9414) 3,'digo(',m,no-1,it,ir,
c-debug;     $                    co(m,1,mo,it,ir,kz),'digo(',m,no-2,
c-debug;     $                    it,ir,co(m,2,mo,it,ir,kz),'digo(',m,
c-debug;     $                    nc-3,it,ir,co(m,3,mo,it,ir,kz),
c-debug;     $                    'oxx',oxx,'oxd(',1,oxd(1),'oxd(',2,
c-debug;     $                    oxd(2),'oxd(',3,oxd(3)
c-debug;                     write(6,9415) b(1),b(2),b(3),'cxx',cxx,
c-debug;     $                    ' cx(',1,cx(1),' cx(',2,cx(2),
c-debug;     $                    'cxdp',0,cxdp,it,ir,opl(it,ir,kz),
c-debug;     $                    'o2  '
c-debug;                  endif
c-debug]
               enddo
            enddo
         endif
c
         return
c
      else if ( no .ge. 5 ) then
c__________			   ! interpolation in regions o3 to o6:
         do i = 4,no-1
c
c	   ! do not go beyond middle (where o3-o6 overlaps c3-c6), and
c	   ! include boundaries (qchk fixes any possible division by 0)
c
            if ( xxo .gt. xod(i)-1.d-6 .and. xxc .ge. xc(i-1)-1.d-6
     $           .and. xod(i-1) .gt. xo(i-1)-1.d-6 ) then
c
               cxdp = log10(zzz(kz)+xcdp)
               m2 = i-2
               m1 = i-1
c
               do it = k1,k1+ip
                  do ir = l1,l1+iq(it-k1+1)
                     b(1) = qchk(is,1,oxx,
     $                    co(m,m2,n(m,m2,kz)-2,it,ir,kz),
     $                    co(m,m2,n(m,m2,kz)-1,it,ir,kz),
     $                    co(m,m2,mo,it,ir,kz),
     $                    ox(n(m,m2,kz)-2),ox(n(m,m2,kz)-1),oxd(m2))
                     b(2) = qchk(is,2,oxx,
     $                    co(m,m1,n(m,m1,kz)-2,it,ir,kz),
     $                    co(m,m1,n(m,m1,kz)-1,it,ir,kz),
     $                    co(m,m1,mo,it,ir,kz),
     $                    ox(n(m,m1,kz)-2),ox(n(m,m1,kz)-1),oxd(m1))
                     b(3) = qchk(is,3,oxx,co(m,m2,mo,it,ir,kz),
     $                    co(m,m1,mo,it,ir,kz),co(m,i,mo,it,ir,kz),
     $                    oxd(m2),oxd(m1),oxd(i))
                     opl(it,ir,kz) = qchk(is,4,cxx,b(1),b(2),b(3),
     $                    cx(m2),cx(m1),cxdp)
                     is = 1
c-debug[
c-debug;                     if ( .not. abs(opl(it,ir,kz)) .le. oudebl .or.
c-debug;     $                    .not. abs(b(1)) .le. oudebl .or.
c-debug;     $                    .not. abs(b(2)) .le. oudebl .or.
c-debug;     $                    .not. abs(b(3)) .le. oudebl .or.
c-debug;     $                    ioudeb .gt. 5 ) then
c-debug;                        write(6,9413) m,kz,xxc,xxo,xxc+xxo,xxco,
c-debug;     $                       xc(nc),nc,(n(m,j,kz),j=1,nc)
c-debug;                        koudeb = koudeb+1
c-debug;                        write(6,9414) 1,'co(m,',m2,
c-debug;     $                       n(m,m2,kz)-2,it,ir,
c-debug;     $                       co(m,m2,n(m,m2,kz)-2,it,ir,kz),
c-debug;     $                       'co(m,',m2,n(m,m2,kz)-1,it,ir,
c-debug;     $                       co(m,m2,n(m,m2,kz)-1,it,ir,kz),
c-debug;     $                       'digo(',m,no-m2,it,ir,
c-debug;     $                       co(m,m2,mo,it,ir,kz),'oxx',oxx,
c-debug;     $                       ' ox(',n(m,m2,kz)-2,
c-debug;     $                       ox(n(m,m2,kz)-2),' ox(',
c-debug;     $                       n(m,m2,kz)-1,ox(n(m,m2,kz)-1),
c-debug;     $                       'oxd(',m2,oxd(m2)
c-debug;                        write(6,9414) 2,'co(m,',m1,
c-debug;     $                       n(m,m1,kz)-2,it,ir,
c-debug;     $                       co(m,m1,n(m,m1,kz)-2,it,ir,kz),
c-debug;     $                       'co(m,',m1,n(m,m1,kz)-1,it,ir,
c-debug;     $                       co(m,m1,n(m,m1,kz)-1,it,ir,kz),
c-debug;     $                       'digo(',m,no-m1,it,ir,
c-debug;     $                       co(m,m1,mo,it,ir,kz),'oxx',oxx,
c-debug;     $                       ' ox(',n(m,m1,kz)-2,
c-debug;     $                       ox(n(m,m1,kz)-2),' ox(',
c-debug;     $                       n(m,m1,kz)-1,ox(n(m,m1,kz)-1),
c-debug;     $                       'oxd(',m1,oxd(m1)
c-debug;                        write(6,9414) 3,'digo(',m,no-m2,it,ir,
c-debug;     $                       co(m,m2,mo,it,ir,kz),'digo(',
c-debug;     $                       m,no-m1,it,ir,
c-debug;     $                       co(m,m1,mo,it,ir,kz),
c-debug;     $                       'digo(',m,no-i,it,ir,
c-debug;     $                       co(m,i,mo,it,ir,kz),'oxx',oxx,
c-debug;     $                       'oxd(',m2,oxd(m2),'oxd(',m1,
c-debug;     $                       oxd(m1),'oxd(',i,oxd(i)
c-debug;                        write(6,9415) b(1),b(2),b(3),'cxx',
c-debug;     $                       cxx,' cx(',m2,cx(m2),' cx(',m1,
c-debug;     $                       cx(m1),'cxdp',0,cxdp,it,ir,
c-debug;     $                       opl(it,ir,kz),'o3-6'
c-debug;                     endif
c-debug]
                  enddo
               enddo
c
               return
c
            endif
c
         enddo
c
      endif
c__________		! else, interpolation in lower left of C-O grid
c
c.....find index of C grid
c                 (must also allow index = nc, to avoid extrapolation)
c
      ie = 100 * max( int(xxc) , int(0.d0) ) + 1
      i3 = max( min( indx(ie) + 1 , nc ) , 3 )
      i1 = i3-2
      i2 = i3-1
c
c.....find index of O grid:
c                   must also allow index = no, to avoid extrapolation
c
      ie = 100 * max( int(xxo) , int(0.d0) ) + 1
      j3 = max( min( indx(ie) + 1 , no ) , 3 )
      j1 = j3-2
      j2 = j3-1
c			! lower-O part of grid: interpolate C before O:
c
      if ( j3 .lt. no .and. i3 .le. n(m,j3,kz) .and. 
     $     ( xxc .lt. xcd(j3)+1.d-6 .or. xxc .ge. xxo ) ) then
c
         do it = k1,k1+ip
            do ir = l1,l1+iq(it-k1+1)
               iw = 0
               do jx = j1,j1+2
                  iw = iw+1
c		 ! if i3 = n(m,jx,kz), must replace cx(i3) with cxd(jx)
                  b(iw) = qchk(is,iw,cxx,co(m,i1,jx,it,ir,kz),
     $                 co(m,i2,jx,it,ir,kz),co(m,i3,jx,it,ir,kz),
     $                 cx(i1),cx(i2),min(cx(i3),cxd(jx)))
               enddo
               iw = iw+1
               opl(it,ir,kz) = qchk(is,iw,oxx,b(1),b(2),b(3),
     $              ox(j1),ox(j2),ox(j3))
               is = 1
c-debug[
c-debug;               if ( .not. abs(opl(it,ir,kz)) .le. oudebl .or.
c-debug;     $              .not. abs(b(1)) .le. oudebl .or.
c-debug;     $              .not. abs(b(2)) .le. oudebl .or.
c-debug;     $              .not. abs(b(3)) .le. oudebl .or.
c-debug;     $              ioudeb .gt. 5 ) then
c-debug;                  write(6,9413) m,kz,xxc,xxo,xxc+xxo,xxco,
c-debug;     $                 xc(nc),nc,(n(m,j,kz),j=1,nc)
c-debug;                  koudeb = koudeb+1
c-debug;                  iw = 0
c-debug;                  do jx = j1,j1+2
c-debug;                     iw = iw+1
c-debug;                     if ( cx(i3) .le. cxd(jx) ) then
c-debug;                        write(6,9414) iw,'co(m,',i1,jx,it,ir,
c-debug;     $                       co(m,i1,jx,it,ir,kz),'co(m,',i2,
c-debug;     $                       jx,it,ir,co(m,i2,jx,it,ir,kz),
c-debug;     $                       'co(m,',i3,jx,it,ir,
c-debug;     $                       co(m,i3,jx,it,ir,kz),
c-debug;     $                       'cxx',cxx,' cx(',i1,cx(i1),
c-debug;     $                       ' cx(',i2,cx(i2),' cx(',i3,cx(i3)
c-debug;                     else
c-debug;                        write(6,9414) iw,'co(m,',i1,jx,it,ir,
c-debug;     $                       co(m,i1,jx,it,ir,kz),'co(m,',i2,
c-debug;     $                       jx,it,ir,co(m,i2,jx,it,ir,kz),
c-debug;     $                       'co(m,',i3,jx,it,ir,
c-debug;     $                       co(m,i3,jx,it,ir,kz),'cxx',cxx,
c-debug;     $                       ' cx(',i1,cx(i1),' cx(',i2,
c-debug;     $                       cx(i2),'cxd(',jx,cxd(jx)
c-debug;                     endif
c-debug;                  enddo
c-debug;                  write(6,9415) b(1),b(2),b(3),'oxx',oxx,
c-debug;     $                 ' ox(',j1,ox(j1),' ox(',j2,ox(j2),
c-debug;     $                 ' ox(',j3,ox(j3),it,ir,opl(it,ir,kz),
c-debug;     $                 'CloO'
c-debug;               endif
c-debug]
            enddo
         enddo
c	      ! else: high-O part of grid: must interpolate O before C:
      else
         do it = k1,k1+ip
            do ir = l1,l1+iq(it-k1+1)
               iw = 0
               do ix = i1,i1+2
                  iw = iw+1
                  if ( j3 .lt. n(m,ix,kz) ) then
                     b(iw) = qchk(is,iw,oxx,co(m,ix,j1,it,ir,kz),
     $                    co(m,ix,j2,it,ir,kz),co(m,ix,j3,it,ir,kz),
     $                    ox(j1),ox(j2),ox(j3))
                  else
                     b(iw) = qchk(is,iw,oxx,co(m,ix,j1,it,ir,kz),
     $                    co(m,ix,j2,it,ir,kz),co(m,ix,mo,it,ir,kz),
     $                    ox(j1),ox(j2),oxd(ix))
                  endif
               enddo
               iw = iw+1
               opl(it,ir,kz) = qchk(is,iw,cxx,b(1),b(2),b(3),
     $              cx(i1),cx(i2),cx(i3))
               is = 1
c-debug[
c-debug;               if ( .not. abs(opl(it,ir,kz)) .le. oudebl .or.
c-debug;     $              .not. abs(b(1)) .le. oudebl .or.
c-debug;     $              .not. abs(b(2)) .le. oudebl .or.
c-debug;     $              .not. abs(b(3)) .le. oudebl .or.
c-debug;     $              ioudeb .gt. 5 ) then
c-debug;                  write(6,9413) m,kz,xxc,xxo,xxc+xxo,
c-debug;     $                 xxco,xc(nc),nc,(n(m,j,kz),j=1,nc)
c-debug;                  koudeb = koudeb+1
c-debug;                  iw = 0
c-debug;                  do ix = i1,i1+2
c-debug;                     iw = iw+1
c-debug;                     if ( j3 .lt. n(m,ix,kz) ) then
c-debug;                        write(6,9414) iw,'co(m,',ix,j1,it,ir,
c-debug;     $                       co(m,ix,j1,it,ir,kz),'co(m,',ix,
c-debug;     $                       j2,it,ir,co(m,ix,j2,it,ir,kz),
c-debug;     $                       'co(m,',ix,j3,it,ir,
c-debug;     $                       co(m,ix,j3,it,ir,kz),
c-debug;     $                       'oxx',oxx,' ox(',j1,ox(j1),
c-debug;     $                       ' ox(',j2,ox(j2),' ox(',j3,ox(j3)
c-debug;                     else
c-debug;                        write(6,9414) iw,'co(m,',ix,j1,it,ir,
c-debug;     $                       co(m,ix,j1,it,ir,kz),'co(m,',ix,
c-debug;     $                       j2,it,ir,co(m,ix,j2,it,ir,kz),
c-debug;     $                       'digo(',m,no-ix,it,ir,
c-debug;     $                       co(m,ix,mo,it,ir,kz),'oxx',oxx,
c-debug;     $                       ' ox(',j1,ox(j1),' ox(',j2,
c-debug;     $                       ox(j2),'oxd(',ix,oxd(ix)
c-debug;                     endif
c-debug;                  enddo
c-debug;                  write(6,9415) b(1),b(2),b(3),'cxx',cxx,
c-debug;     $                 ' cx(',i1,cx(i1),' cx(',i2,cx(i2),
c-debug;     $                 ' cx(',i3,cx(i3),it,ir,opl(it,ir,kz),
c-debug;     $                 'hi-O'
c-debug;               endif
c-debug]
            enddo
         enddo
      endif
c
      return
      end

*
***----------------------------------------------------------------------------
*
      subroutine t6rinterp(slr,slt)
*     =============================
*
*     The purpose of this subroutine is to interpolate in logT6 and logR
*     NOTE THAT for 2-dimensional quadratic interpolation, IT DOES NOT MATTER 
*     which direction is interpolated first, horizontal or vertical: the
*     result is the same for interpolated value and derivatives.
*
***----------------------------------------------------------------------------

      implicit double precision (a-h, o-z)

c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=28, nrb=1, nre=28,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime(mx),cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
c COMMON /d_opal_z/ :
c  dkap = derivative of value returned by quadratic interpolation function QDER
c
      common/d_opal_z/ dkap
      save /d_opal_z/
c
      common/bb_opal_z/ l1,l2,l3,l4,k1,k2,k3,k4,ip,iq(4),xodp,xcdp,
     $     xxco,cxx,oxx,kzf,kzg,kzh,kzf2
      save /bb_opal_z/
c-debug[
c-debug;      common/outdeb/ ioudeb,oudebl,koudeb
c-debug]
c
c===
      dix2 = max(min((alr(l3)-slr)*dfsr(l3),1.d0),0.d0)
      is = 0
      iu = 0
c			   ! for each of the 3 or 4 T values, interpolate in R:
      do kx = k1,k1+ip
c			! interpolate quadratically among (first) 3 R-values
         iu = iu+1
         h(iu) = qder(is,1,slr,opl(kx,l1,1),opl(kx,l2,1),opl(kx,l3,1),
     $        alr(l1),alr(l2),alr(l3))
         q(iu) = dkap
c				  ! if are 4 R-values at this T: mix quadratics
         if ( iq(iu) .eq. 3 ) then
            h(iu) = dix2*h(iu)+(1.d0-dix2)*qder(is,2,slr,opl(kx,l2,1),
     $           opl(kx,l3,1),opl(kx,l4,1),alr(l2),alr(l3),alr(l4))
            q(iu) = dix2*q(iu)+(1.d0-dix2)*dkap
         endif
         is = 1
c-debug[
c-debug;         if ( ioudeb .gt. 3 .or. .not. abs(h(iu)) .le. oudebl .or.
c-debug;     $        .not. abs(q(iu)) .le. oudebl ) then
c-debug;            if ( ioudeb .le. 3 .or. iu .eq. 1 ) write(6,*) ' '
c-debug;            write(6,8912) m,kx,slr,h(iu),q(iu),
c-debug;     $           (kx,j,alr(j),opl(kx,j,1),j=l1,l1+iq(iu))
c-debug; 8912       format(' T6RINTERP(x',i1,',t',i2.2') R',f10.6,
c-debug;     $           ' : K,dKdR =',2g15.7,' <-- it,ir,R,K:',
c-debug;     $           4(i4,i3,f10.6,f11.7))
c-debug;            koudeb = koudeb+1
c-debug;         endif
c-debug]
      enddo
c
c  interpolate in (first) 3 T-values: get opacity, T-derivative, R-derivative:
c  note that QDER and QUAD share the same storage, so calling one sets up for
c  the other as well, and quad(1,...) may correctly follow qder(0,...).
c
      opk(m,1) = qder(0,3,slt,h(1),h(2),h(3),alt(k1),alt(k2),alt(k3))
      opk(m,2) = dkap
      opk(m,3) = quad(1,3,slt,q(1),q(2),q(3),alt(k1),alt(k2),alt(k3))
c
c				    ! if there are 4 T-values: mix quadratics
      if ( ip .eq. 3 ) then
         dix = max(min((alt(k3)-slt)*dfs(k3),1.d0),0.d0)
         opk(m,1) = opk(m,1)*dix+(1.d0-dix)
     $        *qder(0,4,slt,h(2),h(3),h(4),alt(k2),alt(k3),alt(k4))
         opk(m,2) = opk(m,2)*dix+dkap*(1.d0-dix)
         opk(m,3) = opk(m,3)*dix+(1.d0-dix)
     $        *quad(1,4,slt,q(2),q(3),q(4),alt(k2),alt(k3),alt(k4))
      endif
      opk(m,4) = opk(m,2)-3.d0*opk(m,3)
c-debug[
c-debug;      if ( .not. abs(opk(m,1)) .le. oudebl .or. ioudeb .gt. 2 ) then
c-debug;         write(6,8913) m, (1.-dix)*(ip-2), (1.-dix2)*(iq(1)-2),
c-debug;     $        slt,k1,ip,slr,l1,iq(1),iq(2),iq(3),iq(ip+1),
c-debug;     $        opk(m,1),opk(m,2),opk(m,3),opk(m,4),
c-debug;     $        (j,alt(j),h(j-k1+1),q(j-k1+1),j=k1,k1+ip)
c-debug; 8913    format(' '/' T6RINTERP(x',i1,') fhiT',f9.6,' fhiR',f9.6,
c-debug;     $        ' logT6',f10.6,' (',i2,'+',i1,') logR',f10.6,' (',i2,'+',
c-debug;     $        4i1,') logK',g15.7,' DT',f12.7,' DR',f12.7,' DTro',f12.7/
c-debug;     $        ' T6RINTERP(x1) fhiT',f9.6,' fhiR',f9.6,
c-debug;     $        '             <-- it,T,K,dKdR',4(i4,f10.6,2f11.7))
c-debug;         koudeb = koudeb+1
c-debug;      endif
c-debug]
      if ( opk(m,1) .gt. 1.d+5 ) then
         write(6,10) m,10.d0**slt,k1,ip,slr,l1,iq(1),iq(2),iq(3),
     $        iq(ip+1),opk(m,1)
 10      format(' '/' T6RINTERP(m=',i1,') T6=',f15.9,':',i2,'+',i1,
     $        ' logR',f12.7,':',i2,'+',4i1,' logK',e10.3/
     $        '    STOP -- interpolation indices out of range:',
     $        ' PLEASE REPORT CONDITIONS')
         stop
      endif
c
      return
      end



*
***----------------------------------------------------------------------------

      subroutine qzlog4int( zlogd )
*     =============================
*
*..... this subroutine performs bi-quadratic interpolation of logKappa in the
*      log10(Z_i+zdel) values stored in the array zvint(nz), for each of the
*      relevant positions in the C,O-interpolated opacity matrix opl(nt,nr,nz)
*      (given the input values Z and zdel).  Note that this subroutine uses
*      the quadratic-interpolation function quad.  Depending on the number of
*      Z-values to interpolate among, single-quadratic or linear interpolation
*      may be used instead.  Note that zlogd = log10(Z+zdel).
*
*      NOTE that since errors in the opacities may be large compared to the
*      opacity differences between opacities at adjacent Z-values, quadratic
*      interpolation is forced to be monotonic by using values at adjacent
*      tabulated Z-values as upper and lower limits.  No such restriction is
*      placed on extrapolated values.
***----------------------------------------------------------------------------

      implicit double precision (a-h, o-z)

c PARAMETERS to specify opacity storage matrices (see OPALINIT):
c
c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=28, nrb=1, nre=28,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
c COMMON /a_opal_z/ : matrices for opacity storage (see OPALINIT):
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime(mx),cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
c COMMON /bb_opal_z/ : some indices & abundances for T6,R and C,O interpolation
c
      common/bb_opal_z/ l1,l2,l3,l4,k1,k2,k3,k4,ip,iq(4),xodp,xcdp,
     $     xxco,cxx,oxx,kzf,kzg,kzh,kzf2
      save /bb_opal_z/
c
c-debug[
c-debug;      common/outdeb/ ioudeb,oudebl,koudeb
c-debug]
c===
c						! bi-quadratic interpolation:
      if ( kzf2 .gt. kzh ) then
c
         v1 = zvint(kzf)
         v2 = zvint(kzg)
         v3 = zvint(kzh)
         v4 = zvint(kzf2)
         f = ( v3 - zlogd ) / ( v3 - v2 )
         omf = 1.d0 - f
c
         is = 0
         do it = k1, k1 + ip
            do ir = l1, l1 + iq(it-k1+1)
               vkl = opl(it,ir,kzg)
               vkh = opl(it,ir,kzh)
               v0 = max( min(vkl,vkh) , min( max(vkl,vkh) ,
     $              f * quad(is,15,zlogd,opl(it,ir,kzf),vkl,vkh,
     $              v1,v2,v3) + omf * quad(is,30,zlogd,
     $              vkl,vkh,opl(it,ir,kzf2),v2,v3,v4) ) )
               is = 1
c-debug[
c-debug;               if ( ioudeb .gt. 4 .or. .not. abs(v0) .le. oudebl ) then
c-debug;                  if ( ioudeb .le. 4 .or. ( it .eq. k1 .and.
c-debug;     $                 ir .eq. l1 ) ) write(6,*) ' '
c-debug;                  write(6,8912) m,it,ir,zlogd,v0,
c-debug;     $                 (j,zvint(j),opl(it,ir,j),j=kzf,kzf2)
c-debug; 8912             format(' QZLOG4INT(x',i1,',t',i2.2,',r',i2.2,') Z',
c-debug;     $                 f10.6,' : K =',g15.7,' <-- iz,Z,K:',
c-debug;     $                 4(i4,f10.6,f11.7))
c-debug;                  koudeb = koudeb+1
c-debug;               endif
c-debug]
               opl(it,ir,1) = v0
            enddo
         enddo
c						! quadratic interpolation:
      else if ( kzh .gt. kzg ) then
c
         v1 = zvint(kzf)
         v2 = zvint(kzg)
         v3 = zvint(kzh)
         is = 0
         if ( zlogd .le. v2 .and. zlogd .ge. v1 ) then
            do it = k1, k1 + ip
               do ir = l1, l1 + iq(it-k1+1)
                  vkl = opl(it,ir,kzf)
                  vkh = opl(it,ir,kzg)
                  v0 = max( min(vkl,vkh) , min( max(vkl,vkh) ,
     $                 quad(is,30,zlogd,
     $                 vkl,vkh,opl(it,ir,kzh),v1,v2,v3) ) )
                  is = 1
c-debug[
c-debug;                  if ( ioudeb .gt. 4 .or.
c-debug;     $                 .not. abs(v0) .le. oudebl ) then
c-debug;                     if ( ioudeb .le. 4 .or. ( it .eq. k1 .and.
c-debug;     $                    ir .eq. l1 ) ) write(6,*) ' '
c-debug;                     write(6,8912) m,it,ir,zlogd,v0,
c-debug;     $                    (j,zvint(j),opl(it,ir,j),j=kzf,kzh)
c-debug;                     koudeb = koudeb+1
c-debug;                  endif
c-debug]
                  opl(it,ir,1) = v0
               enddo
            enddo
         else if ( zlogd .ge. v2 .and. zlogd .le. v3 ) then
            do it = k1, k1 + ip
               do ir = l1, l1 + iq(it-k1+1)
                  vkl = opl(it,ir,kzg)
                  vkh = opl(it,ir,kzh)
                  v0 = max( min(vkl,vkh) , min( max(vkl,vkh) ,
     $                 quad(is,30,zlogd,
     $                 opl(it,ir,kzf),vkl,vkh,v1,v2,v3) ) )
                  is = 1
c-debug[
c-debug;                  if ( ioudeb .gt. 4 .or.
c-debug;     $                 .not. abs(v0) .le. oudebl ) then
c-debug;                     if ( ioudeb .le. 4 .or. ( it .eq. k1 .and.
c-debug;     $                    ir .eq. l1 ) ) write(6,*) ' '
c-debug;                     write(6,8912) m,it,ir,zlogd,v0,
c-debug;     $                    (j,zvint(j),opl(it,ir,j),j=kzf,kzh)
c-debug;                     koudeb = koudeb+1
c-debug;                  endif
c-debug]
                  opl(it,ir,1) = v0
               enddo
            enddo
         else
            do it = k1, k1 + ip
               do ir = l1, l1 + iq(it-k1+1)
                  v0 = quad(is,30,zlogd,opl(it,ir,kzf),
     $                 opl(it,ir,kzg),opl(it,ir,kzh),v1,v2,v3)
c-debug[
c-debug;                  if ( ioudeb .gt. 4 .or.
c-debug;     $                 .not. abs(v0) .le. oudebl ) then
c-debug;                     if ( ioudeb .le. 4 .or. ( it .eq. k1 .and.
c-debug;     $                    ir .eq. l1 ) ) write(6,*) ' '
c-debug;                     write(6,8912) m,it,ir,zlogd,v0,
c-debug;     $                    (j,zvint(j),opl(it,ir,j),j=kzf,kzh)
c-debug;                     koudeb = koudeb+1
c-debug;                  endif
c-debug]
                  opl(it,ir,1) = v0
                  is = 1
               enddo
            enddo
         endif
c						! linear interpolation:
      else if ( kzg .gt. kzf ) then
c
         f = ( zvint(kzg) - zlogd )
     $        / ( zvint(kzg) - zvint(kzf) )
         omf = 1.d0 - f
         do it = k1, k1 + ip
            do ir = l1, l1 + iq(it-k1+1)
               v0 = f * opl(it,ir,kzf) + omf * opl(it,ir,kzg)
c-debug[
c-debug;               if ( ioudeb .gt. 4 .or. .not. abs(v0) .le. oudebl ) then
c-debug;                  if ( ioudeb .le. 4 .or. ( it .eq. k1 .and.
c-debug;     $                 ir .eq. l1 ) ) write(6,*) ' '
c-debug;                  write(6,8912) m,it,ir,zlogd,v0,
c-debug;     $                 (j,zvint(j),opl(it,ir,j),j=kzf,kzg)
c-debug;                  koudeb = koudeb+1
c-debug;               endif
c-debug]
                  opl(it,ir,1) = v0
            enddo
         enddo
c						! or no interpolation:
      else if ( kzf .ne. 1 ) then
c
         do it = k1, k1 + ip
            do ir = l1, l1 + iq(it-k1+1)
               opl(it,ir,1) = opl(it,ir,kzf)
            enddo
         enddo
c
      endif
c
      return
      end
c


*
***----------------------------------------------------------------------------
*
      double precision function quad(ic,i,x,y1,y2,y3,x1,x2,x3)
*     ========================================================
*
*..... this function performs a quadratic interpolation.
*
*  Storage for dx_i values that need not be computed on each call (see "ic");
*  NOTE that this same storage is used by all of QUAD, QDER, and QCHK.
*

***----------------------------------------------------------------------------

      implicit double precision (a-h, o-z)

      common/coquad_opal_z/ xx12(30),xx13(30),xx23(30),xx1pxx2(30)
      save /coquad_opal_z/
c___
      dimension xx(3),yy(3)
c===
      xx(1) = x1
      yy(1) = y1
      yy(2) = y2
      yy(3) = y3
c			    ! quad may be called many times with same x1,x2,x3;
c			    ! compute & store X-deltas only if flag ic says so:
      if ( ic .eq. 0 ) then
         xx(2) = x2
         xx(3) = x3
         xx12(i) = 1.d0/(xx(1)-xx(2))
         xx13(i) = 1.d0/(xx(1)-xx(3))
         xx23(i) = 1.d0/(xx(2)-xx(3))
         xx1pxx2(i) = xx(1)+xx(2)
      endif
c
      c3 = ( (yy(1)-yy(2))*xx12(i) - (yy(2)-yy(3))*xx23(i) ) * xx13(i)
      c2 = (yy(1)-yy(2))*xx12(i) - xx1pxx2(i) * c3
      c1 = yy(1) - xx(1) * ( c2 + xx(1) * c3 )
      quad = c1+x*(c2+x*c3)
      return
      end



*
***----------------------------------------------------------------------------
*
      double precision function qder(ic,i,x,y1,y2,y3,x1,x2,x3)
*     ========================================================
*
*..... this function performs a quadratic interpolation; it is identical to the
*      function quad, except that it also computes the derivative dkap of the
*      quadratic at the given position x (see  common /d_opal_z/  below).
*
*  COMMON /d_opal_z/ : dkap returns the derivative (in interpolation-direction)
*
***-----------------------------------------------------------------------------

      implicit double precision (a-h, o-z)


      common/d_opal_z/ dkap
      save /d_opal_z/
c
c  Storage for dx_i values that need not be computed on each call (see "ic");
c  NOTE that this same storage is used by all of QUAD, QDER, and QCHK.
c
      common/coquad_opal_z/ xx12(30),xx13(30),xx23(30),xx1pxx2(30)
      save /coquad_opal_z/
c___
      dimension xx(3),yy(3)
c===
      xx(1) = x1
      yy(1) = y1
      yy(2) = y2
      yy(3) = y3
c			    ! qder may be called many times with same x1,x2,x3;
c			    ! compute & store X-deltas only if flag ic says so:
      if ( ic .eq. 0 ) then
         xx(2) = x2
         xx(3) = x3
         xx12(i) = 1.d0/(xx(1)-xx(2))
         xx13(i) = 1.d0/(xx(1)-xx(3))
         xx23(i) = 1.d0/(xx(2)-xx(3))
         xx1pxx2(i) = xx(1)+xx(2)
      endif
c
      c3 = ( (yy(1)-yy(2))*xx12(i) - (yy(2)-yy(3))*xx23(i) ) * xx13(i)
      c2 = (yy(1)-yy(2))*xx12(i) - xx1pxx2(i) * c3
      c1 = yy(1) - xx(1) * ( c2 + xx(1) * c3 )
      dkap = c2+(x+x)*c3
      qder = c1+x*(c2+x*c3)
      return
      end



*
***----------------------------------------------------------------------------
c
      double precision function qchk(ic,i,x,y1,y2,y3,x1,x2,x3)
c     ========================================================
c
c..... this function calls quad(ic,i,x,y1,y2,y3,x1,x2,x3) to perform a
c      quadratic interpolation, but first checks whether any pair of x-values
c      is too close together to make a quadratic interpolation reasonable;
c      if this is the case, something more nearly linear is used instead.
c      Also, for C or O < 0.0 (i.e., x < x1 < x3, as can occur for CNO-depleted
c      matter), linear extrapolation is used.  Note that opacity derivatives
c      are not needed for C-O interpolation, and are not computed.
c
c QCHK is really neaded only for Z slightly less than 0.02, 0.05, or 0.07, or
c   for 0.03 < Z < 0.05 or 0.08 < Z < 0.1, where the C+O=1-X-Z line for one
c   or more of the X-values passes very close above one of the usual C-O grid
c   points; this can result in quadratic interpolation errors in the opacities
c   of more than an order of magnitude.  The solution is to avoid using a
c   quadratic fit if two of the three x-values are too close together.
c
c QCHK is used ONLY for interpolating in C and O; it is not needed elsewhere.
c
c      NOTE that if a quadratic is fitted through 3 points with a large
c      interval R=(x2-x1) with values differing by D=(y2-y1), next to a much
c      smaller interval r=(x3-x2) with values differing by d=(y3-y2), and
c      the close-together points y2 and y3 have a relative error E, then
c      at the middle of the large interval R this error is magnified by
c      a factor of (1/4)(R/r).  At the middle of the interval R, the
c      difference between a linear and a quadratic is (1/4)[D-(R/r)d];
c      if this is less than the magnified error (1/4)(R/r)E, i.e.,
c      if E > | (r/R)D - d | , then the linear fit is better.  For Z < 0.04,
c      the opacity errors should be a few percent, and the RELATIVE error
c      bewteen adjacent nearly-identical compositions may be much smaller:
c      for example, in G91x35z03, compare the following tables:
c TABLE #  7  Grvss'91 (12/92) X=0.3500 Y=0.0200 Z=0.0300 dXc=0.6000 dXo=0.0000
c TABLE # 15  Grvss'91 (12/92) X=0.3500 Y=0.0100 Z=0.0300 dXc=0.6000 dXo=0.0100
c TABLE # 16  Grvss'91 (12/92) X=0.3500 Y=0.0000 Z=0.0300 dXc=0.6100 dXo=0.0100
c      Systematic opacity differences between these tables appear to be of
c      order 0.01 or less in general, and RANDOM error in these differences
c      appears to be at most of order 0.001 (i.e., 0.2 percent).  In QCHK, a
c      quadratic-error magnification of nearly 3 is allowed (R/r=11.5) before
c      beginning to switch over to linear interpolation; the switch-over is
c      complete at R/r=24.  The ratios used in the code are actually r/(r+R)
***----------------------------------------------------------------------------
      
      implicit double precision (a-h, o-z)


      parameter ( ratbeg=0.08d0, ratful=0.04d0, 
     $     ratdel=1.d0/(ratbeg-ratful) )
c
c  Storage for factors that need not be computed on each call (see "ic");
c  NOTE that QCHK calls QUAD: the QUAD/QDER storage is controlled by "ic" too.
c
      common/coqchk_opal_z/ facq(30),iokq(30),iloq(30),i1q(30),i2q(30),
     $     dxinvq(30),xloq(30),flin2(30),flin3(30),lin(30)
      save /coqchk_opal_z/
c___
      dimension xx(3),yy(3),r(3)
c===
c			    ! qchk may be called many times with same x1,x2,x3;
c			    ! compute & store X-deltas only if flag ic says so:
      if ( ic .eq. 0 ) then
         xx(1) = x1
         xx(2) = x2
         xx(3) = x3
         r(1) = abs(xx(3)-xx(2))
         r(2) = abs(xx(3)-xx(1))
         r(3) = abs(xx(2)-xx(1))
         if ( xx(3) - xx(1) .gt. 1.d-6 ) then
            lin(i) = 1
            omf = max( 0.d0 , min( 1.d0 ,
     $           ( xx(3) - xx(1) ) / max( xx(2) - xx(1) , 1.d-6 ) ) )
            flin3(i) = omf / ( xx(3) - xx(1) )
            flin2(i) = ( 1.d0 - omf ) / max( xx(2) - xx(1) , 1.d-6 )
         else
            lin(i) = 0
         endif
         dxrat = min(r(1),r(2),r(3))/max(r(1),r(2),r(3))
         if ( dxrat .ge. ratbeg ) then
            iokq(i) = 1
         else
            if ( r(3) .lt. min(r(1),r(2)) ) then
               iloq(i) = 3
            else if ( r(2) .lt. r(1) ) then
               iloq(i) = 2
            else
               iloq(i) = 1
            endif
            i1q(i) = mod(iloq(i),3)+1
            i2q(i) = 6-i1q(i)-iloq(i)
            xloq(i) = xx(iloq(i))
            dxinvq(i) = 1.d0/((xx(i1q(i))+xx(i2q(i)))*0.5d0-xloq(i))
            if ( dxrat .gt. ratful ) then
               iokq(i) = 0
               facq(i) = (dxrat-ratful)*ratdel
            else
               iokq(i) = -1
               facq(i) = 0.d0
            endif
         endif
      endif
c
      if ( x .lt. x1 ) then
         if ( lin(i) .gt. 0 ) then
            qchk = ( ( y2 - y1 ) * flin2(i) + ( y3 - y1 ) * flin3(i) )
     $           * ( x - x1 ) + y1
            return
         endif
      endif
c
      if ( iokq(i) .gt. 0 ) then
         qchk = quad(ic,i,x,y1,y2,y3,x1,x2,x3)
      else
         yy(1) = y1
         yy(2) = y2
         yy(3) = y3
         if ( iokq(i) .lt. 0 ) then
            qchk = ((yy(i1q(i))+yy(i2q(i)))*0.5d0-yy(iloq(i)))
     $           *(x-xloq(i))*dxinvq(i)+yy(iloq(i))
         else
            qchk = (((yy(i1q(i))+yy(i2q(i)))*0.5d0-yy(iloq(i)))
     $           *(x-xloq(i))*dxinvq(i)+yy(iloq(i)))*(1.d0-facq(i))
     $           +facq(i)*quad(ic,i,x,y1,y2,y3,x1,x2,x3)
         endif
      endif
      return
      end


*
***----------------------------------------------------------------------------
c
      double precision function qzinter(ic,i,z,nmorez,f1,f2,f3,f4,z1,z2,
     $     z3,z4,zdel)
c     ==================================================================
c
c..... this function performs linear, quadratic, or bi-quadratic interpolation,
c      of logKappa in log(Z+zdel), for nmorez = 1, 2, or 3, respectively;
c      inputs are Z, nmorez = one less than the number of Z-values to
c      interpolate among, logKappa values f1 thru f4, Z-values z1 thru z4, and
c      zdel = 0.001 to make things work correctly near Z = 0.  Note that this
c      function is also sometimes used to interpolate in X or C or O.  It makes
c      use of the quadratic-interpolation function quad.
c
c  Storage for values that need not be computed on each call:
***------------------------------------------------------------------------------

      implicit double precision (a-h, o-z)

      common/qzint_opal_z/ v(15,5),f(15),omf(15)
      save /qzint_opal_z/
c===
      if ( ic .eq. 0 ) then
         if ( nmorez .gt. 0 ) then
            if ( zdel .lt. 1.d-5 .or. zdel .gt. 0.1011d0) stop
     $           ' STOP -- QZINTER: bad Zdel value. '
            v(i,1) = log10(z1+zdel)
            v(i,2) = log10(z2+zdel)
            v(i,5) = log10(z+zdel)
            if ( nmorez .eq. 1 ) then
               f(i) = (v(i,2)-v(i,5))/(v(i,2)-v(i,1))
               omf(i) = 1.d0-f(i)
            else
               v(i,3) = log10(z3+zdel)
               if ( nmorez .ge. 3 ) then
                  v(i,4) = log10(z4+zdel)
                  f(i) = (v(i,3)-v(i,5))/(v(i,3)-v(i,2))
                  omf(i) = 1.d0-f(i)
               endif
            endif
         endif
c-debug[
c-debug;      else if ( nmorez .gt. 0 ) then
c-debug;         if ( max( abs( v(i,1) - log10(z1+zdel) ) ,
c-debug;     $        abs( v(i,2) - log10(z2+zdel) ) ,
c-debug;     $        abs( v(i,5) - log10(z+zdel) ) ) .gt. 1.d-5 ) stop
c-debug;     $        ' STOP -- QZINTER: Error: expected same X-values. '
c-debug;         if ( nmorez .eq. 1 ) then
c-debug;            if ( abs( f(i) - (v(i,2)-v(i,5))/(v(i,2)-v(i,1)) ) .gt.
c-debug;     $           1.d-5 ) stop
c-debug;     $           ' STOP -- QZINTER: Error: expected same X-values. '
c-debug;         else
c-debug;            if ( abs( v(i,3) - log10(z3+zdel) ) .gt. 1.d-5 ) stop
c-debug;     $           ' STOP -- QZINTER: Error: expected same X-values. '
c-debug;            if ( nmorez .ge. 3 ) then
c-debug;               if ( max( abs( v(i,4) - log10(z4+zdel) ) ,
c-debug;     $              abs( f(i) - (v(i,3)-v(i,5))/(v(i,3)-v(i,2)) ) )
c-debug;     $              .gt. 1.d-5 ) stop
c-debug;     $              ' STOP -- QZINTER: Error: expected same X-values. '
c-debug;            endif
c-debug;         endif
c-debug]
      endif
c
      if ( nmorez .le. 0 ) then
         qzinter = f1
      else if ( nmorez .eq. 1 ) then
         qzinter = max( min( f(i)*f1 + omf(i)*f2 , max(f1,f2) ) ,
     $        min(f1,f2) )
      else if ( nmorez .eq. 2 ) then
         if ( v(i,5) .lt. v(i,2) ) then
            vlo = min(f1,f2)
            vhi = max(f1,f2)
         else
            vlo = min(f2,f3)
            vhi = max(f2,f3)
         endif
         qzinter = max( min ( quad(ic,i,v(i,5),f1,f2,f3,
     $        v(i,1),v(i,2),v(i,3)) , vhi ) , vlo )
      else
         qzinter = max( min ( f(i)*quad(ic,i,v(i,5),f1,f2,f3,
     $        v(i,1),v(i,2),v(i,3)) + omf(i)*quad(ic,i+15,
     $        v(i,5),f2,f3,f4,v(i,2),v(i,3),v(i,4)) ,
     $        max(f2,f3) ) , min(f2,f3) )
      endif
c
      return
      end

*
***-----------------------------------------------------------------------
*
      function mixfind(iu,iofe,igetzxi,irew,itab,zget,xget,cget,oget)
*     ===============================================================
*
***-----------------------------------------------------------------------

      implicit double precision (a-h, o-z)

c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=28, nrb=1, nre=28,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
c!      parameter ( ks81=ntm-3, ks83=ks81+1, ks60=ks81-21, ks61=ks60+1,
c!     $     alrlo=-8.0d0, flogtlo=3.75d0, flogt60=6.0d0, flogt81=8.1d0 )
      parameter ( ks81=ntm-9, ks83=ks81+1, ks60=ks81-21, ks61=ks60+1,
     $     alrlo=-9.0d0, flogtlo=3.75d0, flogt60=6.0d0, flogt81=8.1d0 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime(mx),cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
c
      character*2 cel_opalmixes(nel_zmix)
      character*8 cfile_opalmixes(n_zmixes)
      common/opalmixes/ xiz_mix(nel_zmix),fninz_mix(nel_zmix),
     $     bracketife_mix(nel_zmix),bracketofe_opalmixes(n_zmixes),
     $     xofe_opalmixes(n_zmixes),xiz_opalmixes(nel_zmix,n_zmixes),
     $     fninz_opalmixes(nel_zmix,n_zmixes),
     $     cel_opalmixes,cfile_opalmixes
      save /opalmixes/
c___
      character*90 cin
      character*1 ch(90)
      equivalence (ch(1),cin)
c===
      ifound = 0
      cin = '                '
c					! if must rewind to beginning of file:
      if ( irew .ne. 0 ) then
         rewind(iu)
c					! else, if must get Xi/Z values:
      else if ( igetzxi .ne. 0 ) then
c
         do while( cin(1:31) .ne. ' Element   Abundance - relative' )
            read(iu,10,end=900) cin
 10         format(a90)
         enddo
         kel = 1
c					! begin loop to get Xi/Z values
         do while( kel .le. nel_zmix )
            read(iu,10,end=900) cin
            if ( cin(1:16) .eq. ' Table Summaries' ) goto 50
            ke = 90
            do while( ke .gt. 1 .and. ch(ke) .eq. ' ' )
               ke = ke-1
            enddo
            kb = 1
            do while( kb .le. ke .and. ch(kb) .eq. ' ' )
               kb = kb+1
            enddo
            if ( cin(kb:kb+1) .eq. cel_opalmixes(kel) ) then
               if ( ke .lt. kb+20 ) goto 50
               read(cin(ke-7:ke),30,err=50) vx
 30            format(f8.6)
               ke = ke-8
               do while( ke .gt. 1 .and. ch(ke) .eq. ' ' )
                  ke = ke-1
               enddo
               if ( ke .lt. kb+11 ) goto 50
               read(cin(ke-7:ke),30,err=50) vn
               if ( igetzxi .gt. 0 ) then
                  xiz_opalmixes(kel,iofe) = vx
                  fninz_opalmixes(kel,iofe) = vn
               else if ( abs(fninz_opalmixes(kel,iofe)-vn) .gt. 5.d-7
     $                 .or. abs(xiz_opalmixes(kel,iofe)-vx) .gt. 1.5d-6
     $                 ) then
                  write(6,40) cel_opalmixes(kel),cfile_opalmixes(iofe),
     $                 vx,xiz_opalmixes(kel,iofe),
     $                 vn,fninz_opalmixes(kel,iofe)
 40               format(' '/' READCO: Warning: new Xi/Z for ',a2,
     $                 ' in ',a8,' mix differs from stored value:'/
     $                 '          new Xi/Z',f9.6,' vs.',f9.6,
     $                 ' ,  new Ni/Nz',f9.6,' vs.',f9.6/' ')
c-dont;                  goto 60
               endif
               kel = kel+1
            endif
c					! end of loop to get Xi/Z values
         enddo
c							! get xO/xFe and [O/Fe]
         if ( igetzxi .gt. 0 ) then
            xofe_opalmixes(iofe) = xiz_opalmixes(kel_o,iofe)
     $           /max(xiz_opalmixes(kel_fe,iofe),1.d-36)
            bracketofe_opalmixes(iofe) = log10(xofe_opalmixes(iofe)
     $           /xofe_opalmixes(1))
         endif
c					! no read error: jump to continuation
         goto 60
c				! if error reading Xi/Z values, say so
 50      write(6,20) kel
 20      format(' '/' READCO: Warning: error reading',
     $        ' Z-abundance fractions at element',i3/' ')
c								! continuation
 60      continue
c				! end of reading Xi/Z values in file header
      endif
c							! find start of tables
      if ( irew .ne. 0 .or. igetzxi .ne. 0 ) then
         do while( cin(1:30) .ne. '******************************' )
            read(iu,10,end=900) cin
         enddo
         igetzxi = 0
      endif
c				  ! look for mix with required composition:
      do while ( ifound .eq. 0 )
         read(iu,10,end=900) cin
         if ( cin(1:7) .eq. 'TABLE #' ) then
            ke = 90
            do while( ke .gt. 1 .and. ch(ke) .eq. ' ' )
               ke = ke-1
            enddo
            if ( ke .lt. 60 ) goto 900
            read(cin(ke-48:ke),100) xat,yat,zat,cat,oat
 100        format(3(3x,f6.4),2(5x,f6.4))
            if ( max(abs(zat-zget),abs(xat-xget),
     $           abs(cat-cget),abs(oat-oget)) .lt. 1.d-6 ) ifound = 1
         endif
      enddo
c				  ! found required mix: read its table number
      read(cin(8:10),105) itabat
 105  format(i3)
c		 ! if it does not consecutively follow previous table, may need
c		 !  to rewind back to beginning of file for next composition
      irew = 0
      if ( itabat .ne. itab+1 ) irew = 1
      itab = itabat
c					! check log R values in table head
      do i = 1,3
         read(iu,10,end=900) cin
      enddo
c      read(iu,*,err=900,end=900) cin(1:4),(alrf(i),i=1,nrm)
c      write(*,*)cin(1:4),(alrf(i),i=1,nrm)
      read(iu,110,err=900,end=900) cin(1:4),(alrf(i),i=1,nrm)
c!! 110  format(a4,f6.1,18f7.1)
 110  format(a4,f6.1,27f7.1)
      if ( cin(1:4) .ne. 'logT' .or.
     $     abs(alrf(1)-alrlo) .gt. 1.d-5 ) goto 900
      do k = 2,nrm
         if ( abs(alrf(k)-alrf(k-1)-0.5d0) .gt. 1.d-5 ) stop
     $        ' STOP -- READCO: bad  log R  value in table read in. '
      enddo
c				! read blank line before first table line
      read(iu,10,end=900) cin
c					! return
 900  mixfind = ifound
c-debug-chk[
c-debug-chk;      if ( ifound .eq. 0 )
c-debug-chk;     $     write(6,1739) iu,itab,irew,zget,xget,cget,oget
c-debug-chk; 1739 format(' '/' MIXFIND: unit',i3,' after TABLE',i3,
c-debug-chk;     $     ', irew=',i2,': could not find mix Z=',f10.7,
c-debug-chk;     $     ' X=',f10.7,' C=',f10.7,' O=',f10.7)
c-debug-chk]
      if ( ifound .eq. 0 ) irew = 1
      return
      end
c



***----------------------------------------------------------------
c
      SUBROUTINE SPLINE(X,Y,N,Y2)
c     ===========================
c
***----------------------------------------------------------------
      implicit double precision (a-h , o-z)

      PARAMETER ( NMAX=100 )
C
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
C
C     FIRST DERIVATIVES AT END POINTS USING CUBIC FIT
C===
      YP1 = ((Y(3)-Y(1))*(X(2)-X(1))**2
     $     -(Y(2)-Y(1))*(X(3)-X(1))**2)/
     $     ((X(3)-X(1))*(X(2)-X(1))*(X(2)-X(3)))
      YPN = ((Y(N-2)-Y(N))*(X(N-1)-X(N))**2
     $     -(Y(N-1)-Y(N))*(X(N-2)-X(N))**2)/
     $     ((X(N-2)-X(N))*(X(N-1)-X(N))*(X(N-1)-X(N-2)))
C
      Y2(1) = -0.5d0
      U(1) = (3.d0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      DO I = 2,N-1
         SIG = (X(I)-X(I-1))/(X(I+1)-X(I-1))
         P = SIG*Y2(I-1)+2.d0
         Y2(I) = (SIG-1.d0)/P
         U(I) = (6.d0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     $        /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      ENDDO
      QN = 0.5d0
      UN = (3.d0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      Y2(N) = (UN-QN*U(N-1))/(QN*Y2(N-1)+1.d0)
      DO K = N-1,1,-1
         Y2(K) = Y2(K)*Y2(K+1)+U(K)
      ENDDO
      RETURN
      END


*
***----------------------------------------------------------------------

      SUBROUTINE SPLINT(XA,YA,N,Y2A,X,Y,YP)
*     =====================================
***-----------------------------------------------------------------------

      implicit double precision (a-h, o-z)

      DIMENSION XA(N),YA(N),Y2A(N)
C===
      KLO = 1
      KHI = N
      do while( KHI-KLO .GT. 1 )
         K = (KHI+KLO)/2
         IF ( XA(K) .GT. X ) THEN
            KHI = K
         ELSE
            KLO = K
         ENDIF
      enddo
      H = XA(KHI)-XA(KLO)
      IF ( H .EQ. 0.d0 ) STOP ' STOP -- SPLINT: Bad XA input. '
      A = (XA(KHI)-X)/H
      B = (X-XA(KLO))/H
      Y = A*YA(KLO)+B*YA(KHI)+
     $     ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.d0
      YP = 0.05d0*  (  (-YA(KLO)+YA(KHI))/H
     $     +      ( -(3*A**2-1)*Y2A(KLO)
     $     +(3*B**2-1)*Y2A(KHI) )*H/6.d0 )
      RETURN
      END

*
***----------------------------------------------------------------------
      SUBROUTINE FITY
c     ===============
C
C  THIS ROUTINE MAKES SPLINE FITS FOR F AND FX, AND OBTAINS
C  FY AND FXY
C						! modified:
*
***----------------------------------------------------------------------


      implicit double precision (a-h, o-z)

      COMMON/CST_OPAL_Z/ NRL,RLS,nset,tmax
      save /CST_OPAL_Z/
C
      PARAMETER ( IPR=28 )
C
      COMMON/CF_OPAL_Z/ F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)
      save /CF_OPAL_Z/
c___
      DIMENSION A(IPR),B(IPR),AD(IPR),BD(IPR)
C===								    ! modified:
      DO I = 1,nset
         DO J = 1,NRL
            A(J) = F(I,J)
            B(J) = FX(I,J)
         ENDDO
C
         CALL GETD(A,NRL,AD,AP1,APN)
         CALL GETD(B,NRL,BD,BP1,BPN)
C
         FY(I,1) = AP1
         FY(I,NRL) = APN
         FXY(I,1) = BP1
         FXY(I,NRL) = BPN
         DO J = 2,NRL-1
            FY(I,J) =  -A(J)+A(J+1)-2.d0*AD(J)-AD(J+1)
            FXY(I,J) = -B(J)+B(J+1)-2.d0*BD(J)-BD(J+1)
         ENDDO
      ENDDO
C
      RETURN
      END

*
***----------------------------------------------------------------------
      SUBROUTINE FITX
c     ===============
C
C  THIS ROUTINE IS USED ONLY AFTER SMOOTHING.
C  ITS FUNCTION IS TO RECOMPUTE FX USING SMOOTHED F.
C
*
***----------------------------------------------------------------------
      implicit double precision (a-h, o-z)

      PARAMETER ( IPR=28 )
C					! modified:
      COMMON/CST_OPAL_Z/ NRL,RLS,nset,tmax
      save /CST_OPAL_Z/
c
      COMMON/CF_OPAL_Z/ F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)
      save /CF_OPAL_Z/
C___
      DIMENSION A(85),D(85)
C===
      DO J = 1,NRL
C					! modified:
         DO I = 1,nset
            A(I) = F(I,J)
         ENDDO
C					! modified:
         CALL GETD(A,nset,D,AP1,APN)
         FX(1,J) = AP1
C					! modified:
         FX(nset,J) = APN
C					! modified:
         DO I = 2,nset-1
            FX(I,J) = -A(I)+A(I+1)-2.d0*D(I)-D(I+1)
         ENDDO
      ENDDO
C
      RETURN
      END

*
***----------------------------------------------------------------------
      SUBROUTINE GETD(F,N,D,FP1,FPN)
c     ==============================
C
C  SIMPLIFIED CODE FOR SPLINE COEFFICIENTS, FOR CASE OF INTERVALS OF UNITY.
*
***----------------------------------------------------------------------


      implicit double precision (a-h, o-z)

      DIMENSION F(N),D(N),T(85)
C===
      FP1 = (-11.d0*F(1)+18.d0*F(2)-9.d0*F(3)+2.d0*F(4))/6.d0
      FPN = (11.d0*F(N)-18.d0*F(N-1)+9.d0*F(N-2)-2.d0*F(N-3))/6.d0
C
      D(1) = -0.5d0
      T(1) = 0.5d0*(-F(1)+F(2)-FP1)
C
      DO J = 2,N-1
         D(J) = -1.d0/(4.d0+D(J-1))
         T(J) = -D(J)*(F(J-1)-2.d0*F(J)+F(J+1)-T(J-1))
      ENDDO
C
      D(N) = (FPN+F(N-1)-F(N)-T(N-1))/(2.d0+D(N-1))
C
      DO J = N-1,1,-1
         D(J) = D(J)*D(J+1)+T(J)
      ENDDO
C
      RETURN
      END


*
***------------------------------------------------------------------------
c
      SUBROUTINE INTERPOPA(FLT,FLRHO,G,DGDT,DGDRHO,IERR)
c     ===============================================
C
C  GIVEN F,FX,FY AND FXY ON THE GRID POINTS, THIS ROUTINE
C  DOES BI-CUBIC INTERPOLATIONS USING METHODS DESCRIBED IN
C  Numerical Recipes, PP. 118 TO 120
C
*
***------------------------------------------------------------------------

      implicit double precision (a-h, o-z)

      PARAMETER ( IPR=28 )
C						! modified:
      COMMON/CST_OPAL_Z/ NRL,RLS,nset,tmax
      save /CST_OPAL_Z/
C
      COMMON/CF_OPAL_Z/ F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)
      save /CF_OPAL_Z/
c___
      DIMENSION B(16)
      LOGICAL IERR
C
C  EXTREME LIMITS ALLOWED ARE:-
C     (3.800-0.0125) TO (8.000+0.0125) FOR LOG10(T)
C     (RLS-0.125) TO (RLE+0.1254) FOR LOG10(R)
C     (ALLOWING FOR SMALL EXTRAPOLATIONS BEYOND TABULAR VALUES)
C
C     FUNCTION DEFINITIONS FOR CUBIC EXPANSION
C===
      FF(S,T) =   B( 1)+T*(B( 2)+T*(B( 3)+T*B( 4)))
     $     +S*(   B( 5)+T*(B( 6)+T*(B( 7)+T*B( 8)))
     $     +S*(   B( 9)+T*(B(10)+T*(B(11)+T*B(12)))
     $     +S*(   B(13)+T*(B(14)+T*(B(15)+T*B(16))) )))
C
      FFX(S,T) =   B( 5)+T*(B( 6)+T*(B( 7)+T*B( 8)))
     $     +S*( 2*(B( 9)+T*(B(10)+T*(B(11)+T*B(12))))
     $     +S*( 3*(B(13)+T*(B(14)+T*(B(15)+T*B(16)))) ))
C
      FFY(S,T) =   B( 2)+S*(B( 6)+S*(B(10)+S*B(14)))
     $     +T*( 2*(B( 3)+S*(B( 7)+S*(B(11)+S*B(15))))
     $     +T*( 3*(B( 4)+S*(B( 8)+S*(B(12)+S*B(16)))) ))
C
C  Note that statement function FFXY is never used, and thus can be omitted.
C
c-noneed[						 ! FFXY is never used!
c-noneed;      FFXY(S,T) =  B( 6)+T*(2*B( 7)+3*T*B( 8))
c-noneed;     $     +S*(  2*B(10)+T*(4*B(11)+6*T*B(12))
c-noneed;     $     +S*(  3*B(14)+T*(6*B(15)+9*T*B(16)) ))
c-noneed]
C
C     BEGINNING OF EXECUTABLE STATEMENTS
C===
      IERR = .FALSE.
C
      X = 20.d0*(FLT-3.800d0)+1.d0
      FLR = FLRHO+18.d0-3.d0*FLT
      Y = 2.d0*( FLR - RLS )+1.d0
C
      IF ( X .LT. 2.d0 ) THEN
         IF ( X .LT. 0.75d0 ) THEN
            IERR = .TRUE.
         ELSE
            I = 1
         ENDIF
      ELSE IF ( X .GT. 84.d0 ) THEN
         IF ( X .GT. 85.25d0 ) THEN
            IERR = .TRUE.
         ELSE
            I = 84
         ENDIF
      ELSE
         I = INT(X)
      ENDIF
      U = X-I
C
      IF ( Y .LT. 2.d0 ) THEN
         IF ( Y .LT. 0.75d0 ) THEN
            IERR = .TRUE.
         ELSE
            J = 1
         ENDIF
      ELSE IF ( Y .GT. NRL-1 ) THEN
         IF ( Y .GT. NRL+.25 ) THEN
            IERR = .TRUE.
         ELSE
            J = NRL-1
         ENDIF
      ELSE
         J = INT(Y)
      ENDIF
      V = Y-J
C
      IF ( IERR ) THEN
         G = 9.999d0
         DGDT = 9.999d0
         DGDRHO = 9.999d0
         RETURN
      ENDIF
C
C  GIVEN FUNCTIONS AND DERIVATIVES AT GRID POINTS, COMPUTE COEFFICIENTS.
c
      B(1) = F(I,J)
      B(2) = FY(I,J)
      B(3) = 3*(-F(I,J)+F(I,J+1))-2*FY(I,J)-FY(I,J+1)
      B(4) = 2*(F(I,J)-F(I,J+1))+FY(I,J)+FY(I,J+1)
C
      B(5) = FX(I,J)
      B(6) = FXY(I,J)
      B(7) = 3*(-FX(I,J)+FX(I,J+1))-2*FXY(I,J)-FXY(I,J+1)
      B(8) = 2*(FX(I,J)-FX(I,J+1))+FXY(I,J)+FXY(I,J+1)
C
      B(9) = 3*(-F(I,J)+F(I+1,J))-2*FX(I,J)-FX(I+1,J)
      B(10) = 3*(-FY(I,J)+FY(I+1,J))-2*FXY(I,J)-FXY(I+1,J)
      B(11) = 9*(F(I,J)-F(I+1,J)+F(I+1,J+1)-F(I,J+1))
     $     +6*(FX(I,J)-FX(I,J+1)+FY(I,J)-FY(I+1,J))
     $     +4*FXY(I,J)
     $     +3*(FX(I+1,J)-FX(I+1,J+1)-FY(I+1,J+1)+FY(I,J+1))
     $     +2*(FXY(I,J+1)+FXY(I+1,J))
     $     +FXY(I+1,J+1)
      B(12) = 6*(-F(I,J)+F(I+1,J)-F(I+1,J+1)+F(I,J+1))
     $     +4*(-FX(I,J)+FX(I,J+1))
     $     +3*(-FY(I,J)+FY(I+1,J)+FY(I+1,J+1)-FY(I,J+1))
     $     +2*(-FX(I+1,J)+FX(I+1,J+1)-FXY(I,J)-FXY(I,J+1))
     $     -FXY(I+1,J)-FXY(I+1,J+1)
C
      B(13) = 2*(F(I,J)-F(I+1,J))+FX(I,J)+FX(I+1,J)
      B(14) = 2*(FY(I,J)-FY(I+1,J))+FXY(I,J)+FXY(I+1,J)
      B(15) = 6*(-F(I,J)+F(I+1,J)-F(I+1,J+1)+F(I,J+1))
     $     +4*(-FY(I,J)+FY(I+1,J))
     $     +3*(-FX(I,J)-FX(I+1,J)+FX(I+1,J+1)+FX(I,J+1))
     $     +2*(FY(I+1,J+1)-FY(I,J+1)-FXY(I,J)-FXY(I+1,J))
     $     -FXY(I+1,J+1)-FXY(I,J+1)
      B(16) = 4*(F(I,J)-F(I+1,J)+F(I+1,J+1)-F(I,J+1))
     $     +2*(FX(I,J)+FX(I+1,J)-FX(I+1,J+1)-FX(I,J+1)
     $     +FY(I,J)-FY(I+1,J)-FY(I+1,J+1)+FY(I,J+1))
     $     +FXY(I,J)+FXY(I+1,J)+FXY(I+1,J+1)+FXY(I,J+1)
C
C  GET G=LOG10(ROSS), DGDT=d LOG10(ROSS)/d LOG10(T),
C      DGDRHO=d LOG10(ROSS)/d LOG10(RHO)
c
      G = FF(U,V)
      DGDT = 20.d0*FFX(U,V)-6.d0*FFY(U,V)
      DGDRHO = 2.d0*FFY(U,V)
C
      RETURN
      END

*
***-----------------------------------------------------------------------
c
      SUBROUTINE SMOOTH
c     =================
C
C  THIS SUBROUTINE USES A 2-DIMENSIONAL GENERALISATION OF THE SMOOTHING
C  TECHNIQUES DESCRIBED ON PP. 644 TO 649 OF Numerical Recipes.
C
C  CONSIDER THE 25 POINTS DEFINED BY
C       I+n, n=-2,-1,0,1,2 AND J+m, m=-2,-1,0,1,2.
C  THE FUNCTION TO BE SMOOTHED IS FITTED TO A BI-CUBIC, INVOLVING
C  16 COEFFICIENTS, USING TECHNIQUES OF LEAST-SQUARES. THE SMOOTHED
C  FUNCTION (TEMPORARILY STORED IN FXY) IS GIVEN BY THE FITTED VALUE
C  AT THE POINT I AND J.
C
C  THE FITTING IS SHIFTED FOR POINTS CLOSE TO BOUNDARIES.
C
*
***-----------------------------------------------------------------------

      implicit double precision (a-h, o-z)

      PARAMETER ( IPR=28 )
C						! modified
      COMMON/CST_OPAL_Z/ NRL,RLS,nset,tmax
      save /CST_OPAL_Z/
c
      COMMON/CF_OPAL_Z/ F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)
      save /CF_OPAL_Z/
C___
      DIMENSION GAM(6),BET(11),ALP(11)
c---

      DATA GAM/+0.0073469388d0,-0.0293877551d0,-0.0416326531d0,
     $     +0.1175510204d0,+0.1665306122d0,+0.2359183673d0/
      DATA BET/
     $     -0.0048979592d0,-0.0661224490d0,-0.0293877551d0,
     $     0.0195918367d0,0.2644897959d0,0.1175510204d0,-0.0783673469d0,
     $     0.0277551020d0,0.3746938776d0,0.1665306122d0,-0.1110204082d0/
      DATA ALP/
     $     -0.0844897959d0,-0.0048979592d0,0.0073469388d0,
     $     0.0012244898d0,0.3379591837d0,0.0195918367d0,-0.0293877551d0,
     $     0.4787755102d0,0.0277551020d0,-0.0416326531d0,-0.0069387755d0
     $      /

      DO I = 3,nset-2 
C
         J = 1
         FXY(I,J) = 
     $        ALP(1)*( F(I-2,J  )+F(I+2,J  ) )
     $        +ALP(2)*( F(I-2,J+1)+F(I+2,J+1)+F(I-2,J+3)+F(I+2,J+3)
     $        +F(I-1,J+4)+F(I+1,J+4) )
     $        +ALP(3)*( F(I-2,J+2)+F(I+2,J+2) )
     $        +ALP(4)*( F(I-2,J+4)+F(I+2,J+4) )
     $        +ALP(5)*( F(I-1,J  )+F(I+1,J  ) )
     $        +ALP(6)*( F(I-1,J+1)+F(I+1,J+1)+F(I-1,J+3)+F(I+1,J+3) )
     $        +ALP(7)*( F(I-1,J+2)+F(I+1,J+2) )
     $        +ALP(8)*  F(I  ,J  )
     $        +ALP(9)*( F(I  ,J+1)+F(I  ,J+3) )
     $        +ALP(10)* F(I  ,J+2) +ALP(11)*F(I  ,J+4)
C
         J = 2
         FXY(I,J) = 
     $        BET(1)*( F(I-2,J-1)+F(I+2,J-1)+F(I-2,J+3)+F(I+2,J+3) )
     $        +BET(2)*( F(I-2,J  )+F(I+2,J  ) )
     $        +BET(3)*( F(I-2,J+1)+F(I+2,J+1) )
     $        +BET(4)*( F(I-2,J+2)+F(I+2,J+2)+F(I-1,J-1)+F(I+1,J-1)
     $        +F(I-1,J+3)+F(I+1,J+3) )
     $        +BET(5)*( F(I-1,J  )+F(I+1,J  ) )
     $        +BET(6)*( F(I-1,J+1)+F(I+1,J+1) )
     $        +BET(7)*( F(I-1,J+2)+F(I+1,J+2) )
     $        +BET(8)*( F(I  ,J-1)+F(I  ,J+3) )
     $        +BET(9)*F(I  ,J  ) +BET(10)*F(I  ,J+1) +BET(11)*F(I  ,J+2)
C
         DO J = 3,NRL-2
            FXY(I,J) = 
     $           GAM(1)*( F(I-2,J-2)+F(I-2,J+2)+F(I+2,J-2)+F(I+2,J+2) )
     $           +GAM(2)*( F(I-2,J+1)+F(I-2,J-1)+F(I-1,J-2)+F(I-1,J+2)
     $           +F(I+1,J-2)+F(I+1,J+2)+F(I+2,J-1)+F(I+2,J+1) )
     $           +GAM(3)*( F(I-2,J  )+F(I+2,J  )+F(I  ,J-2)+F(I  ,J+2) )
     $           +GAM(4)*( F(I-1,J-1)+F(I-1,J+1)+F(I+1,J-1)+F(I+1,J+1) )
     $           +GAM(5)*( F(I-1,J  )+F(I  ,J-1)+F(I  ,J+1)+F(I+1,J  ) )
     $           +GAM(6)*  F(I  ,J  )
         ENDDO
C
         J = NRL-1
         if (J.le.0) then            
c 		print *, I,J,NRL,FXY(I,J),rle,rls
		stop 'livopaint: array bounds exceded'
         endif
         FXY(I,J) = 
     $        BET(1)*( F(I-2,J+1)+F(I+2,J+1)+F(I-2,J-3)+F(I+2,J-3) )
     $        +BET(2)*( F(I-2,J  )+F(I+2,J  ) )
     $        +BET(3)*( F(I-2,J-1)+F(I+2,J-1) )
     $        +BET(4)*( F(I-2,J-2)+F(I+2,J-2)+F(I-1,J+1)+F(I+1,J+1)
     $        +F(I-1,J-3)+F(I+1,J-3) )
     $        +BET(5)*( F(I-1,J  )+F(I+1,J  ) )
     $        +BET(6)*( F(I-1,J-1)+F(I+1,J-1) )
     $        +BET(7)*( F(I-1,J-2)+F(I+1,J-2) )
     $        +BET(8)*( F(I  ,J+1)+F(I  ,J-3) )
     $        +BET(9)*F(I  ,J  ) +BET(10)*F(I  ,J-1) +BET(11)*F(I  ,J-2)
C
         J = NRL
         FXY(I,J) = 
     $        ALP(1)*( F(I-2,J  )+F(I+2,J  ) )
     $        +ALP(2)*( F(I-2,J-1)+F(I+2,J-1)+F(I-2,J-3)+F(I+2,J-3)
     $        +F(I-1,J-4)+F(I+1,J-4) )
     $        +ALP(3)*( F(I-2,J-2)+F(I+2,J-2) )
     $        +ALP(4)*( F(I-2,J-4)+F(I+2,J-4) )
     $        +ALP(5)*( F(I-1,J  )+F(I+1,J  ) )
     $        +ALP(6)*( F(I-1,J-1)+F(I+1,J-1)+F(I-1,J-3)+F(I+1,J-3) )
     $        +ALP(7)*( F(I-1,J-2)+F(I+1,J-2) )
     $        +ALP(8)*  F(I  ,J  )
     $        +ALP(9)*( F(I  ,J-1)+F(I  ,J-3) )
     $        +ALP(10)* F(I  ,J-2) +ALP(11)*F(I  ,J-4)
C
      ENDDO
C			! modified
      DO I = 3,nset-2
         DO J = 1,NRL
            F(I,J) = FXY(I,J)
         ENDDO
      ENDDO
C
      RETURN
      END


C
c******************************************************************************
c
      subroutine opaltab
c     ==================
C
C  CODE FOR FITTING AND SMOOTHING OPAL DATA. ADAPTED FROM A CODE
C     WRITTEN BY MIKE SEATON (obtained june 1993)
C
C     OPAL DATA.
C     ASSUMES FIRST T6 = .0056341325 = 10.**(3.75-6.) ; LAST T6 = tmax = 10.
C     USES RECTANGULAR ARRAY FOR VARIABLES T6 AND LOG10(R)
C
C     (1) NSM=NUMBER OF PASSES THROUGH SMOOTHING FILTER.
C     USE OF NSM=1 OR 2 IS RECOMMENDED.
C     NO SMOOTHING WITH NSM=0
C     (2) RANGE FOR LOG10(R),
C     RLS=FIRST VALUE, RLE=LAST VALE
C     (RLS MUST BE FIRST VALUYE IN TABLE)
C
C  SUBROUTINE INTERPOPA
C     AFTER PROCESSING, DATA ARE IN A FORM FOR USE OF
C               SUBROUTINE INTERPOPA
C     WHICH GIVES LOG(ROSS) AND TWO FIRST DERIVATIVES FOR ANY
C     VALUES OF LOG(T) AND LOG(RHO). SEE BELOW FOR FURTHER
C     EXPLANATION.
C
C  OUTPUT FOR THE CASE OF NSM .GT. 0.
C     INTERP IS USED TO OBTAIN SMOOTHED DATA INTERPOLATED
C     BACK TO THE ORIGINAL OPAL MESH.
C
C  THE SUBROUTINES SPLINE AND SPLINT ARE ADAPTED FROM THOSE GIVE BY
C  W.H. Press, S.A. Teulolsky, W.T. Vettering and B.P. Flannery,
C  "Numerical Recipes in FORTRAN", 2nd edn., 1992, C.U.P.
C  OTHER REFERENCES ARE MADE TO METHODS DESCRIBED IN THAT BOOK.
***----------------------------------------------------------------

      implicit double precision (a-h, o-z)
C
      PARAMETER ( IP=100, IPR=28 )
c
c      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
c!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=1, mx=5, mc=8, mo=8, nrm=28, nrb=1, nre=28,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/alink_opal_z/ NTEMP,NSM,nrlow,nrhigh,RLE,t6arr(100),
     $     coff(100,nrm)
      save /alink_opal_z/
c						! modified:
      COMMON/CST_OPAL_Z/ NRL,RLS,nset,tmax
      save /CST_OPAL_Z/
c
      COMMON/CF_OPAL_Z/ F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)
      save /CF_OPAL_Z/
c___
      DIMENSION U(IP),ROSSL(IP,IPR),V(IP),V2(IP)
c
c-noneed;      CHARACTER*1 HEAD(100)
c
      LOGICAL IERR
C===
      NRL = int(2.d0*(RLE-RLS)+1.00001d0)
      if (NRL.gt.IPR) then
         write(*,*)  ' change IPR =',ipr,' to IPR =',nrl
         stop
      endif
C
C     STORE LOG10(T) IN U AND LOG10(ROSS) IN ROSSL
C     CHECK FIRST VALUE OF T6
c
      do j = 1,NRL
         ROSSL(1,j) = coff(1,j)
      enddo
c
      T6 = t6arr(1)
      U(1) = 6.d0+LOG10(T6)
c
C     SET ROSSL UP TO T6=t6arr(nset)
c
      NTEMP = 1
      do while( T6 .LT. Tmax )
         NTEMP = NTEMP+1
         do j = 1,NRL
            ROSSL(NTEMP,j) = coff(NTEMP,j)
         enddo
         T6 = t6arr(NTEMP)
         U(NTEMP) = 6.d0+LOG10(T6)
      enddo
c
      IF ( NTEMP .GT. IP ) THEN
         PRINT*,' OPALTAB: REQUIRE PARAMETER IP OF AT LEAST ',NTEMP
         STOP ' STOP -- OPALTAB: NTEMP > IP . '
      ENDIF
C
C
C     DEFINE VARIABLES
C         X=20.0*(LOG10(T)-3.80)+1
C         Y=2.0*(LOG10(R)-RLS)+1
C     USE INDICES I=1 TO nset AND J=1 TO NRL
C     X AND Y ARE SUCH THAT, ON MESH-POINT (I,J), X=I AND Y=J
C     OBTAIN:-
C         F(I,J)=LOG10(ROSS)
C         FX(I,J)=dF/dX
C         FY(I,J)=dF/dY
C         FXY(I,J)=ddF/dXdY
C
C
C     FIRST GET F AND FX, INTERPOLATING FROM OPAL T6 TO
C     INTERVAL OF 0.05 IN LOG10(T).
c
      DO J = 1,NRL
c
C        FOR EACH LOG10(R), STORE LOG10(ROSS) IN V(I)
c
         DO I = 1,NTEMP
            V(I) = ROSSL(I,J)
         ENDDO
C
C        GET FIRST DERIVATIVES AT END POINTS: done in SPLINE, using cubic fit
C
C        GET SECOND DERIVATIVES FOR SPLINE FIT: done by SPLINE
c
         CALL SPLINE(U,V,NTEMP,V2)
C
C        INTERPOLATE TO LOG10(T)=FLT, FLT=3.8(0.05)8.0
c							! modified:
         DO I = 1,nset
            FLT = 3.75d0+0.05d0*I
            CALL SPLINT(U,V,NTEMP,V2,FLT,F(I,J),FX(I,J))
         ENDDO
C
      ENDDO
C
C  OPTION FOR SMOOTHING
C
      IF ( NSM .GT. 0 ) THEN
         DO NS = 1,NSM
            CALL SMOOTH
         ENDDO
         CALL FITX
      ENDIF
C
C  GET FY AND FXY
C
      CALL FITY
C
C  THE ARRAYS F, FX, FY AND FXY ARE NOW STORED
C
C  CAN NOW DO INTERPOLATIONS USING
C       CALL INTERPOPA(FLT,FLRHO,G,DGDT,DGDRHO,IERR)
C       INPUT IS FLT=LOG10(T), FLRHO=LOG10(RHO)
C       OUTPUT IS G=LOG10(ROSS)
C              DGDT=dG/d(LOG10(T))
C            DGDRHO=dG/d(LOG10(RHO))
C              IERR=.TRUE. IF INPUT FLT, FLRHO ARE OUT-OF-RANGE,
C                          ELSE IERR=.FALSE.
C
C INTERPOLATE BACK TO OPAL POINTS
C
      IF ( NSM .GT. 0 ) THEN
c
         do il = 1,NRL
            coff(1,il) = ROSSL(1,il)
         enddo
c
         DO K = 2,NTEMP
            FLT = U(K)
            DO IL = nrlow,nrhigh
               FLR = RLS+0.5d0*(IL-1)
               FLRHO = FLR-18.d0+3.d0*FLT
               CALL INTERPOPA(FLT,FLRHO,G,DGDT,DGDRHO,IERR)
               IF ( IERR ) THEN
                  stop ' STOP -- OPALTAB: INTERP T/rho range error. '
               ENDIF
               V(IL) = G
            ENDDO
            do il = nrlow,nrhigh
               coff(K,il) = V(il)
            enddo
         ENDDO
c
      ENDIF
C
      return
      END
c


c***********************************************************************
c  NOTE that the subroutines below have alternate parts, all but one
c  commented out, for various flavors of UNIX and for VMS (the Linux
c  versions should actually work for any flavor of UNIX).
c***********************************************************************
*
***-----------------------------------------------------------------------
c
      subroutine open_chk_zip( iu, fname, igzip, cmsg )
c     -------------------------------------------------

***-----------------------------------------------------------------------
c
      character*(*) fname
      character*(*) cmsg
c
      character*160 fnalt
      character*180 cgzip
      logical lxst
c
      call inqfil( fname, lxst )
c
      if ( lxst ) then
         igzip = 0
      else
c-vms[
c-vms;         fnalt = fname
c-vms]
c-sun-iris-linux[
         if ( fname(1:2) .eq. '~/' ) then
            call getenv( 'HOME' , fnalt )
            last = 1
            do while ( last .lt. 160 .and. fnalt(last:last) .ne. ' ' )
               last = last + 1
            enddo
            fnalt(last:) = fname(2:)
         else
            fnalt = fname
         endif
c-sun-iris-linux]
         last = 160
         do while ( last .gt. 1 .and. fnalt(last:last) .eq. ' ' )
            last = last - 1
         enddo
         if ( fnalt(last:last) .eq. ' ' ) stop
     $        ' STOP -- Error: blank file name. '
         cgzip = 'gunzip ' // fnalt(1:last) // '.gz'
         call inqfil( cgzip(8:), lxst )
         if ( lxst ) then
            igzip = 1
            call system( cgzip )
         else
            cgzip = 'uncompress ' // fnalt(1:last) // '.Z'
            call inqfil( cgzip(12:), lxst )
            if ( lxst ) then
               igzip = -1
               call system( cgzip )
            else
               write(6,'(a)') cmsg
               stop
            endif
         endif
      endif
c
      call opoldr( iu, fname )
c
      return
      end


*
***----------------------------------------------------
c
      subroutine close_chk_zip( iu, fname, igzip )
c     --------------------------------------------

***-----------------------------------------------------
c
      character*(*) fname
c
      character*160 fnalt
c
      close(iu)
c
      if ( igzip .ne. 0 ) then
c-vms[
c-vms;         fnalt = fname
c-vms]
c-sun-iris-linux[
         if ( fname(1:2) .eq. '~/' ) then
            call getenv( 'HOME' , fnalt )
            last = 1
            do while ( last .lt. 160 .and. fnalt(last:last) .ne. ' ' )
               last = last + 1
            enddo
            fnalt(last:) = fname(2:)
         else
            fnalt = fname
         endif
c-sun-iris-linux]
c
         if ( igzip .gt. 0 ) then
c
            call system( 'gzip ' // fnalt )
c
         else if ( igzip .lt. 0 ) then
c
            call system( 'compress ' // fnalt )
c
         endif
c
      endif
c
      return
      end


*
***---------------------------------------------------

      subroutine opoldr(iu,fname)
*     ---------------------------
*
* Open an old formatted file:
*
***----------------------------------------------------

      character*80 fname

c
c-linux[
      character*160 fnalt
c-linux]
c
c-vms[                                                !  For VMS:
c-vms;      data cb/ ':', ']', ';', ';', ';', ';' /
c-vms]
c-sun-iris-linux[                                     ! For UNIX:
c-sun-iris-linux]
c
c  For Linux: get home directory name if necessary, and open the file
c  with the err= keyword to prevent coredump on file open error
c  (actually, this should work on any Unix system, provided that the
c  environment variable HOME is correctly defined as the home directory):
c
c-linux[
      if ( fname(1:2) .eq. '~/' ) then
         call getenv( 'HOME' , fnalt )
         i = 1
         do while ( i .lt. 160 .and. fnalt(i:i) .ne. ' ' )
            i = i + 1
         enddo
         fnalt(i:) = fname(2:)
      else
         fnalt = fname
      endif
      open(iu,file=fnalt,form='FORMATTED',status='OLD',
     $     iostat=ioperr,err=900)

      i = 1
      do while ( i .lt. 160 .and. fnalt(i:i) .ne. ' ' )
         i = i + 1
      enddo
      write(*,*) 'reading opacity file: ',fnalt(1:i)
c      write(*,*) 'reading opacity file: ',fnalt(1:30)
      return
c
 900  write(6,910) ioperr,iu,fnalt
 910  format(' '/' Error',i12,' opening unit',i3,' with old file:'/
     $     ' ',a160)
      stop ' STOP -- OPOLDR: Error opening old file. '
c-linux]
c
c  For Sun UNIX: open the file:
c
c-sun[
c-sun;      open(iu,file=fname,form='FORMATTED',status='OLD')
c-sun;      return
c-sun]
c
c  For VMS, or for SGI Iris UNIX: open the file as read-only:
c
c-vms-iris[
c-vms-iris;      open(iu,file=fname,form='FORMATTED',status='OLD',
c-vms-iris;     $     readonly)
c-vms-iris;      return
c-vms-iris]
c
      end
*
***----------------------------------------------------------------------
*
      subroutine inqfil(fname,lxst)
*     -----------------------------
***----------------------------------------------------------------------

      character*(*) fname
      logical lxst
c
c  For Linux: get home directory name, if necessary
c  (actually, this should work on any Unix system, provided that the
c  environment variable HOME is correctly defined as the home directory):
c
c-linux[
      character*250 fnalt
c
      fnalt = fname
c
      if ( fnalt(1:2) .eq. '~/' ) then
         fnalt = ' '
         call getenv( 'HOME' , fnalt )
         i = 1
         do while ( i .lt. 250 .and. fnalt(i:i) .ne. ' ' )
            i = i + 1
         enddo
         fnalt(i:) = fname(2:)
      endif
c

      inquire( file = fnalt, exist = lxst )


c-linux]
c
c  Anything except Linux: just look for filename as is:
c
c-sun-vms-iris[
c-sun-vms-iris;      inquire( file = fname , exist = lxst )
c-sun-vms-iris]
c
      return
      end
c
c************************************************************************
c
