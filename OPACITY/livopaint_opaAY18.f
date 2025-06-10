c******************************************************************************
c
      block data opal_opac_data
c     =========================
***   Modified for STAREVOL
c
c  COMMON BLOCK DATA INITIALIZATIONS:
c  ----------------------------------
c
c PARAMETERS defining the matrices used to hold the opacity values:
c  nz = 14 = maximum number of Z-tabulation values (see arrays zavail, zsto)
c  mx = 5 = number of X-tabulation values (see array xa)
c  mc = mo = 8 = number of C- and O-tabulation values (see arrays xcs and xos)
c  nrm = 39 = maximum number of R-tabulation values (see array alrf)
c  nrb = 1 = index of lowest R-value to be stored         \  default: store
c  nre = 39 = index of highest R-value to be stored        }  all 39 R-values
c  nr = nre-nrb+1 = 39 = number of R-values to be stored  /
c  ntm = 76 = maximum number of T-tabulation values (see array flogtin)
c  ntb = 1 = index of lowest T-value to be stored         \  default: store
c  nt = ntm-ntb+1 = 76 = number of T-values to be stored  /   all 76 T-values
c

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!    $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
c PARAMETERS nrdel and ntdel give matrix position differences that result from
c  any reduction of nr or nt due to increased nrb or ntb values, respectively
c
      parameter ( nrdel=nrb-1, ntdel=ntb-1 )
c
c PARAMETERS:
c  zdel = 0.001 = offset for Z, Z+C, and Z+O, to make log interpolation behave
c                  reasonably at small Z values: Z-interpolation is performed
c                  using log(Z+zdel), while the CO-interpolation is performed
c                  using log(C+Z+zdel) and log(O+Z+zdel)
c  xdel = 0.03 = usual (high-T) offset for X, to make log interpolation behave
c                 reasonably at small X; note that 0.03 works better than 0.005
c  xdelmin = 0.001 = lowest value of X offset ever used (at low temperatures)
c
      parameter ( zdel=0.001d0, xdel=0.03d0, xdelmin=0.001d0 )
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
c  xcno_use(mx_hi,nz) = the 'GN93hz' X-values available for each stored Z-value
c                         (indexed by kz = 1, ..., numz); note that for each
c                         value of kz the highest X-value is 1 - Z(kz)
c  xhi_use(mx_hi,nz) = the 'GN93hz' X-values same as xcno_use(mx_hi,nz), except
c                        that xhi_use(1,kz) = 0.03 (from 'Gz???.x03' files)
c  xxx_cno(mx_hi) = log10( xhi_in() + xdel ) = logarithmic 'GN93hz' X-values
c  xxx_hi(mx_hi) = logarithmic 'GN93hz' X-values same as xxx_cno(mx_hi), except
c                    that xxx_hi(1) = log10( 0.03 + xdel )
c  nx_hi(nz) = number of 'GN93hz' X-values in xhi_use() at each stored Z-value
c  ireq_hi(mx_hi) = flags to tell whether the corresponding 'GN93hz' X-values
c                     are unavailable from the 'Gz???.x??' files
c  khighx(nz) = flag to tell whether the 'GN93hz' opacities were read in, for
c                 each of the stored Z-values
c  kavail_xhi = flag to tell whether the 'GN93hz' opacities are available
c  kuse_xhi = flag to tell whether the 'GN93hz' X-values should be used for
c               X-interpolation (see description of subroutine SET_KXHI above)
c  kdo_xhi = internal flag controlling use of the 'GN93hz' X-values
c  kavail_cno = flag to tell whether the CNO-interpolation deltas are available
c  kuse_cno = flag to tell whether CNO-interpolation should be performed
c  kdo_cno = internal flag controlling use of CNO-interpolation
c  kavail_user = flag to tell whether user-interpolation deltas are available
c  kuse_user = flag to tell whether user-interpolation should be performed
c  kdo_user = internal flag controlling use of user-interpolation
c
c /xhi_opal_z/: --> data{ALL}
      common /xhi_opal_z/ xhi_in(mx_hi), xcno_use(mx_hi,nz),
     $     xhi_use(mx_hi,nz), xxx_cno(mx_hi), xxx_hi(mx_hi),
     $     nx_hi(nz), ireq_hi(mx_hi), khighx(nz), kavail_xhi, kuse_xhi,
     $     kdo_xhi, kavail_cno, kuse_cno, kdo_cno, kavail_user,
     $     kuse_user, kdo_user
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
c  itime = "opacities available" flag (initially 0): itime is set to 12345678
c             when all opacities have been read in.
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
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
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
      common/b_opal_z/ nta(0:nrm_p1),ntax0(0:nrm),
     $     ntax03(0:nrm), sltlo, slthi, dltlo_inv, dlthi_inv,
     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
      save /b_opal_z/
c
c COMMON /e_opal_z/ : variables returning opacity values, as described above
c
c /e_opal_z/: --> data{ALL}
      common/e_opal_z/ opact,dopact,dopacr,dopactd,fedge,ftredge,fzedge
      save /e_opal_z/
c
c COMMON /x_opal_z/ : variables containing stored OPAL_F_CNOU input values
c  
c /x_opal_z/: --> data{ALL}
      common /x_opal_z/ z_opal, x_opal, xc_opal, xo_opal, slt_opal,
     $     slr_opal, fcn_opal, fcon_opal, fcnone_opal, fu_opal
      save /x_opal_z/
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
c COMMON /alt_change_opal_z/ :
c  main_alt_change = 0 unless a new file has been set to replace 'GS98hz'
c  iulow = 23 (by default) = the beginning Fortran unit number for reading
c            opacity files; Fortran units iulow through iulow + 3 may be used.
c  khighz_in = khighz value used when reading opacities.
c  ofebrack_in = [O/Fe] value used when reading opacities.
c
c /alt_change_opal_z/: --> data{ALL}
      common /alt_change_opal_z/ main_alt_change, iulow, khighz_in,
     $     ofebrack_in
      save /alt_change_opal_z/
c
c PARAMETERS defining the matrices used to hold the mix compositions:
c  nel_zmix = 19 = number of heavy elements in the opacity mixture (C thru Ni)
c  n_zmixes = 8 = number of "hz" opacity files available (the fifth one is for
c                  a user-supplied "hz" opacity file, with non-zero [O/Fe])
c  kel_o = 3 = position of oxygen (O) in the list of mix-elements
c  kel_fe = nel_zmix-1 = 18 = position of iron (Fe) in the list of mix-elements
c
!!      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
!!      parameter ( nel_zmix=19, n_zmixes=7, kel_o=3, kel_fe=nel_zmix-1 )
      parameter ( nel_zmix=19, n_zmixes=8, kel_o=3, kel_fe=nel_zmix-1 )
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
c  n_totmix = 10 = n_zmixes + 5 = number of "GS98hz" opacity files available,
c                                   plus 5 more files for CNO-interpolation and
c                                   user-specified interpolation.
c  n_cnobeg = 6 = n_zmixes + 1 = position of first CNO-interpolation filename
c
      parameter ( n_totmix = n_zmixes + 5, n_cnobeg = n_zmixes + 1 )
c
c COMMON /opalGS98mixes/ : makeup of Z in GS98 opacity mixtures (see above):
c
c /opalGS98mixes/: --> data{ALL BUT xofe_opalGS98}
c-implicit;      real*4 bracketofe_opalGS98, xofe_opalGS98, xiz_opalGS98,
c-implicit;     $     fninz_opalGS98, atwt_opalGS98
      character*255 cfile_opalGS98(n_totmix)
      common /opalGS98mixes/ bracketofe_opalGS98(n_totmix),
     $     xofe_opalGS98(n_totmix),xiz_opalGS98(nel_zmix,n_totmix),
     $     fninz_opalGS98(nel_zmix,n_totmix),atwt_opalGS98(nel_zmix),
     $     cfile_opalGS98
      save /opalGS98mixes/
c
c COMMON /ext_CNO_opal_z/ : default extensions for CNO-interpolation files:
c
c /ext_CNO_opal_z/: --> data{ALL}
      character*8 cdef_CNO_ext(n_cnobeg:n_totmix)
      common /ext_CNO_opal_z/ len_def_CNO_ext(n_cnobeg:n_totmix),
     $     cdef_CNO_ext
      save /ext_CNO_opal_z/
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
c           = 1 : on input, subroutine OPALTAB smooths the opacities, but do
c                  NOT smooth the CNO-interpolation opacity deltas
c           = 2 (default): on input, subroutine OPALTAB smooths the opacities,
c                  including those used to get the CNO-interpolation deltas
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
c COMMON /opdir/ : copdir = the name of the directory holding opacitythe files;
c  here it is initialized to be blank, meaning "use the current directory".
c
c /opdir/: --> data{copdir}
      character*255 copdir
      common/opdir/ copdir
      save /opdir/
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
c                3 : In this case, it is also a fatal error if CNO-interp files
c                      have been read in and you call OPAL_X_CNO_U with a metal
c                      composition array of size other than 19.
c
c /c_level_err_opal_z/: --> data{level_err}
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c
c COMMON /c_ctab/
c  ctab = character*1 variable containing a tab character: needed only in order
c           to allow tabs in input files (and in character string argumentss)
c           to be interpreted as whitespace
c
c /c_ctab/: --> data{ctab}
      character*1 ctab
      common /c_ctab/ ctab
      save /c_ctab/
c
c COMMON /chkpoc/ :
c  cb = character(s) allowed to terminate a directory name [cb(1) and cb(2)],
c         and non-alphanumeric characters allowed in a filename [cb(3) thru 6].
c         DEFAULT (sun/iris/linux): cb(1:2) = '/', cb(3:6) = '_', '~', '+', '-'
c         VMS (must edit data statement): cb(1:2) = ':', ']', cb(3:6) = ';'
c
c /chkpoc/: --> data{cb}
      character*1 cb(6)
      common/chkpoc/cb
      save /chkpoc/
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
c
c /xhi_opal_z/ data:
c
      data xhi_in / 0.d0, 0.1d0, 0.2d0, 0.35d0, 0.5d0, 0.7d0, 0.8d0,
     &     0.9d0, 0.95d0, 1.d0 /
      data xcno_use / mx_hi_nz * -1.0d0 /, xhi_use / mx_hi_nz * -1.0d0 /
      data xxx_cno / mx_hi * -9.0d0 /, xxx_hi / mx_hi * -9.0d0 /
      data nx_hi / nz * 0 /
      data ireq_hi / 0, 0, 1, 0, 1, 0, 1, 1, 1, 1 /
      data khighx / nz * 0 /
      data kavail_xhi / 0 /, kuse_xhi / 2 /, kdo_xhi / 0 /
      data kavail_cno / 0 /, kuse_cno / 1 /, kdo_cno / 0 /
      data kavail_user / 0 /, kuse_user / 1 /, kdo_user / 0 /
c
c /a_opal_z/ data:
c								! indx(1:101)
      data indx/1,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4,
     $     4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,
     $     6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
     $     7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7/
c								! X(1:mx) mx=5
      data xa/ 0.0d0, 0.03d0, 0.1d0, 0.35d0, 0.7d0, 3 * 0.0d0 /
c								! C,O(1:mc)
      data xcs/0.0d0,0.01d0,0.03d0,0.1d0,0.2d0,0.4d0,0.6d0,1.0d0/
      data xos/0.0d0,0.01d0,0.03d0,0.1d0,0.2d0,0.4d0,0.6d0,1.0d0/
c								! init-flag
      data itime / 0 /
c
c /b_opal_z/ data:
c								! nta(0:nrm_p1)
!!      data nta/57, 70,70,70,70,70, 70,70,70,70,70,
!!     $     70,70,70,70,69, 64,60,58,57, -99/
c								! ntax0(0:nrm)

      data nta/57,76,76,76,76,76,76,76,76,76,76,76,76,76,76,     ! nta(0:nrm_p1)
     $     76,76,76,76,76,76,76,76,76,76,76,76,76,76, 
     &     76,76,76,76,76,76,76,76,76,76,76, -99/
c
c								! ntax0(0:nrm)
!!      data ntax0/999, 6,5,5,5,4, 4,4,3,1,1, 1,1,1,1,1, 1,1,1,1/
      data ntax0/999, 1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,
     $     1,1,1,1,1, 1,1,1,1 ,1,1,1,1,1 ,1,1,1,1,1,1/  
c								! ntax03(0:nrm)
!!      data ntax03/999, 5,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1, 1,1,1,1/
      data ntax03/999, 5,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,
     $     1,1,1,1,1, 1,1,1,1 ,1,1,1,1,1 ,1,1,1,1,1,1/
c
      data sltlo/-99.d0/, slthi/-99.d0/, dltlo_inv/-99.d0/, 
     &     dlthi_inv/-99.d0/
      data slrlo/-99.d0/, slrhi/-99.d0/, dlrlo_inv/-99.d0/, 
     &     dlrhi_inv/-99.d0/
      data init_trvals/0/
c
c /e_opal_z/ data:
c
      data opact/1.d35/, dopact/0.d0/, dopacr/0.d0/, dopactd/0.d0/,
     $     fedge/0.d0/, ftredge/0.d0/, fzedge/0.d0/
c
c /x_opal_z/ data:
c
      data z_opal / -1.0d0 /, x_opal / -1.0d0 /, xc_opal / -1.0d0 /
      data xo_opal / -1.0d0 /, slt_opal / -99.0d0 /, slr_opal / -99.0d0/
      data fcn_opal / 0.0d0 /, fcon_opal / 0.0d0 /, fcnone_opal / 0.0d0/
      data fu_opal / 0.0d0 /
c
c /recoin_opal_z/ data:
c
      data itimeco/0/, mxzero/1/, mx03/2/, kope/0/, igznotgx/0/
c
c /alt_change_opal_z/ data:
c
      data main_alt_change / 0 /, iulow / 23 /, khighz_in / 0 /,
     $     ofebrack_in / 0.0d0 /
c
c /opalmixes/ data:
c
      data cfile_opalmixes/'AY18hz','AGSS09hz','AGS05hz ','GN93hz  '
     $     ,'Alrd96a2','C95hz   ','W95hz   ','        '/
      data cel_opalmixes/'C ','N ','O ','Ne','Na','Mg','Al','Si',
     $     'P ','S ','Cl','Ar','K ','Ca','Ti','Cr','Mn','Fe','Ni'/
      data xiz_opalmixes/
     $     0.170552d0,0.049872d0,0.413275d0,0.127717d0,0.002100d0,
     $     0.050988d0,0.004004d0,0.047897d0,0.000420d0,0.022259d0,
     $     0.000589d0,0.005273d0,0.000220d0,0.004634d0,0.000226d0,
     $     0.001195d0,0.000779d0,0.092849d0,0.005151d0, 
     $     0.176857d0,0.051822d0,0.428745d0,0.093972d0,0.002186d0,
     $     0.052950d0,0.004160d0,0.049731d0,0.000436d0,0.023130d0,
     $     0.000613d0,0.005490d0,0.000229d0,0.004797d0,0.000234d0,
     $     0.001242d0,0.000809d0,0.096628d0,0.005331d0,
     $     0.177050d0,0.050682d0,0.439156d0,0.083833d0,0.002041d0,
     $     0.049455d0,0.003798d0,0.054575d0,0.000426d0,0.026575d0,
     $     0.000674d0,0.003632d0,0.000281d0,0.004914d0,0.000230d0,
     $     0.001363d0,0.000809d0,0.094518d0,0.005987d0,
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
c     $     0.173285d0,0.053152d0,0.482273d0,0.098668d0,0.001999d0,
c     $     0.037573d0,0.003238d0,0.040520d0,0.000355d0,0.021142d0,
c     $     0.000456d0,0.005379d0,0.000210d0,0.003734d0,0.000211d0,
c     $     0.001005d0,0.000548d0,0.071794d0,0.004459d0/
c AGS04
c     $     0.177050d0,0.050682d0,0.439156d0,0.083833d0,0.002041d0,
c     $     0.049455d0,0.003798d0,0.054575d0,0.000426d0,0.026575d0,
c     $     0.000674d0,0.003632d0,0.000281d0,0.004914d0,0.000230d0,
c     $     0.001363d0,0.000809d0,0.094518d0,0.005987d0/
c AGSS09
c     $     0.176857d0,0.051822d0,0.428745d0,0.093972d0,0.002186d0,
c     $     0.052950d0,0.004160d0,0.049731d0,0.000436d0,0.023130d0,
c     $     0.000613d0,0.005490d0,0.000229d0,0.004797d0,0.000234d0,
c     $     0.001242d0,0.000809d0,0.096628d0,0.005331d0/
c AY18      
     $     0.170552d0,0.049872d0,0.413275d0,0.127717d0,0.002100d0,
     $     0.050988d0,0.004004d0,0.047897d0,0.000420d0,0.022259d0,
     $     0.000589d0,0.005273d0,0.000220d0,0.004634d0,0.000226d0,
     $     0.001195d0,0.000779d0,0.092849d0,0.005151d0/
      data fninz_opalmixes/
     $     0.250290d0,0.062761d0,0.455302d0,0.111561d0,0.001610d0,
     $     0.036977d0,0.002616d0,0.030060d0,0.000239d0,0.012238d0,
     $     0.000293d0,0.002327d0,0.000099d0,0.002038d0,0.000083d0,
     $     0.000405d0,0.000250d0,0.029305d0,0.001547d0, 
     $     0.258495d0,0.064931d0,0.470383d0,0.081743d0,0.001669d0,
     $     0.038234d0,0.002707d0,0.031078d0,0.000247d0,0.012661d0,
     $     0.000304d0,0.002412d0,0.000103d0,0.002101d0,0.000086d0,
     $     0.000419d0,0.000258d0,0.030371d0,0.001594d0,
     $     0.258191d0,0.063378d0,0.480774d0,0.072768d0,0.001555d0,
     $     0.035640d0,0.002466d0,0.034062d0,0.000241d0,0.014519d0,
     $     0.000333d0,0.001592d0,0.000126d0,0.002147d0,0.000084d0,
     $     0.000459d0,0.000258d0,0.029644d0,0.001786d0,
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
c     $     0.245518d0,0.064578d0,0.512966d0,0.083210d0,0.001479d0,
c     $     0.026308d0,0.002042d0,0.024552d0,0.000195d0,0.011222d0,
c     $     0.000219d0,0.002291d0,0.000091d0,0.001586d0,0.000075d0,
c     $     0.000329d0,0.000170d0,0.021877d0,0.001293d0/
c AGS04
c     $     0.258191d0,0.063378d0,0.480774d0,0.072768d0,0.001555d0,
c     $     0.035640d0,0.002466d0,0.034062d0,0.000241d0,0.014519d0,
c     $     0.000333d0,0.001592d0,0.000126d0,0.002147d0,0.000084d0,
c     $     0.000459d0,0.000258d0,0.029644d0,0.001786d0/
c AGSS09
c     $     0.258495d0,0.064931d0,0.470383d0,0.081743d0,0.001669d0,
c     $     0.038234d0,0.002707d0,0.031078d0,0.000247d0,0.012661d0,
c     $     0.000304d0,0.002412d0,0.000103d0,0.002101d0,0.000086d0,
c     $     0.000419d0,0.000258d0,0.030371d0,0.001594d0/
c AY18      
     $     0.250290d0,0.062761d0,0.455302d0,0.111561d0,0.001610d0,
     $     0.036977d0,0.002616d0,0.030060d0,0.000239d0,0.012238d0,
     $     0.000293d0,0.002327d0,0.000099d0,0.002038d0,0.000083d0,
     $     0.000405d0,0.000250d0,0.029305d0,0.001547d0/ 
      data bracketife_mix/nel_zmix*0.0d0/
c      data bracketofe_opalmixes/0.0d0,0.3d0,0.4d0,0.5d0,0.0d0,0.0d0,
c     &     0.0D0/
      data bracketofe_opalmixes/0.0d0,0.3d0,0.4d0,0.5d0,0.0d0,0.0d0,
     &     0.0D0,0.0D0/
c
c /opalGS98mixes/ data:
c
c      data bracketofe_opalGS98 / 0.0d0, 0.3d0, 0.4d0, 0.5d0, 
c     $     8 * 0.0d0 /
      data bracketofe_opalGS98 / 0.0d0, 0.3d0, 0.4d0, 0.5d0, 
     $     9 * 0.0d0 /
      data xiz_opalGS98 /
     $     0.171836d0, 0.050335d0, 0.467356d0, 0.104831d0, 0.002090d0,
     $     0.039924d0, 0.003603d0, 0.044057d0, 0.000423d0, 0.023513d0,
     $     0.000292d0, 0.004335d0, 0.000228d0, 0.003896d0, 0.000195d0,
     $     0.001117d0, 0.000779d0, 0.076433d0, 0.004757d0,
     $     0.101930d0, 0.029858d0, 0.553139d0, 0.124072d0, 0.002473d0,
     $     0.047252d0, 0.001071d0, 0.052144d0, 0.000501d0, 0.027829d0,
     $     0.000346d0, 0.005131d0, 0.000135d0, 0.004611d0, 0.000226d0,
     $     0.000663d0, 0.000327d0, 0.045339d0, 0.002955d0,
     $     0.091638d0, 0.026843d0, 0.626052d0, 0.055905d0, 0.001114d0,
     $     0.053480d0, 0.001921d0, 0.059017d0, 0.000226d0, 0.031497d0,
     $     0.000156d0, 0.002312d0, 0.000305d0, 0.005218d0, 0.000256d0,
     $     0.000596d0, 0.000165d0, 0.040761d0, 0.002537d0,
     $     0.076359d0, 0.022368d0, 0.656744d0, 0.090832d0, 0.000929d0,
     $     0.044563d0, 0.001601d0, 0.039063d0, 0.000188d0, 0.022339d0,
     $     0.000130d0, 0.001927d0, 0.000101d0, 0.005474d0, 0.000362d0,
     $     0.000496d0, 0.000346d0, 0.033965d0, 0.002214d0,
     $     0.171836d0, 0.050335d0, 0.467356d0, 0.104831d0, 0.002090d0,
     $     0.039924d0, 0.003603d0, 0.044057d0, 0.000423d0, 0.023513d0,
     $     0.000292d0, 0.004335d0, 0.000228d0, 0.003896d0, 0.000195d0,
     $     0.001117d0, 0.000779d0, 0.076433d0, 0.004757d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,  
     $     nel_zmix * 0.0d0, nel_zmix * 0.0d0, nel_zmix * 0.0d0,
     $     nel_zmix * 0.0d0, nel_zmix * 0.0d0 /
      data fninz_opalGS98 /
     $     0.245825d0, 0.061748d0, 0.501922d0, 0.089265d0, 0.001562d0,
     $     0.028224d0, 0.002294d0, 0.026954d0, 0.000235d0, 0.012602d0,
     $     0.000141d0, 0.001865d0, 0.000100d0, 0.001670d0, 0.000070d0,
     $     0.000369d0, 0.000244d0, 0.023517d0, 0.001393d0,
     $     0.148069d0, 0.037193d0, 0.603216d0, 0.107280d0, 0.001877d0,
     $     0.033921d0, 0.000693d0, 0.032394d0, 0.000282d0, 0.015145d0,
     $     0.000170d0, 0.002241d0, 0.000060d0, 0.002007d0, 0.000082d0,
     $     0.000222d0, 0.000104d0, 0.014165d0, 0.000878d0,
     $     0.131883d0, 0.033128d0, 0.676395d0, 0.047890d0, 0.000838d0,
     $     0.038036d0, 0.001231d0, 0.036324d0, 0.000126d0, 0.016982d0,
     $     0.000076d0, 0.001000d0, 0.000135d0, 0.002251d0, 0.000092d0,
     $     0.000198d0, 0.000052d0, 0.012616d0, 0.000747d0,
     $     0.108877d0, 0.027349d0, 0.702986d0, 0.077089d0, 0.000692d0,
     $     0.031401d0, 0.001016d0, 0.023820d0, 0.000104d0, 0.011933d0,
     $     0.000063d0, 0.000826d0, 0.000044d0, 0.002339d0, 0.000129d0,
     $     0.000163d0, 0.000108d0, 0.010416d0, 0.000646d0,
     $     0.245825d0, 0.061748d0, 0.501922d0, 0.089265d0, 0.001562d0,
     $     0.028224d0, 0.002294d0, 0.026954d0, 0.000235d0, 0.012602d0,
     $     0.000141d0, 0.001865d0, 0.000100d0, 0.001670d0, 0.000070d0,
     $     0.000369d0, 0.000244d0, 0.023517d0, 0.001393d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0,
     &     0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     $     nel_zmix * 0.0d0, nel_zmix * 0.0d0, nel_zmix * 0.0d0,
     $     nel_zmix * 0.0d0, nel_zmix * 0.0d0 /
      data atwt_opalGS98 /
     $     12.01100d0, 14.00670d0, 15.99940d0, 20.17900d0, 22.98977d0,
     $     24.30500d0, 26.98154d0, 28.08550d0, 30.97376d0, 32.06000d0,
     $     35.45300d0, 39.94800d0, 39.09830d0, 40.08000d0, 47.90000d0,
     $     51.99600d0, 54.93800d0, 55.84700d0, 58.70000d0 /
      data cfile_opalGS98 / 'GS98hz', 'GS98hz_OFe.3_Alrd96a2',
     $     'GS98hz_OFe.4_C95', 'GS98hz_OFe.5_W95', 9 * ' ' /
c
c /ext_CNO_opal_z/ data:
c
      data len_def_CNO_ext / 0, 5, 6, 8, 5 /
      data cdef_CNO_ext / '        ',
     $     '.CtoN   ', '.COtoN  ', '.CNOtoNe', '.user   ' /
c
c /zinter_opal_z/ data: (note: zavail, iadd_zavail are computed in get_zavail)
c
      data zavail / nzm * 0.0d0 /, iadd_zavail / nadd_zavail * 0 /
c
      data zvalhi /0.d0,0.01d0,0.02d0,0.03d0,0.04d0,0.05d0,0.06d0,0.07d0
     $    ,0.08d0,0.09d0,0.1d0/
c
      data zval / 0.d0,0.001d0,0.004d0,0.01d0,0.02d0,0.03d0,0.05d0,
     $     0.1d0 /
c
      data zalval / 0.d0,0.0001d0,0.0003d0,0.001d0,0.002d0,0.004d0,
     $     0.01d0,0.02d0,0.03d0,0.04d0,0.06d0,0.08d0,0.1d0 /
c
      data (mnofz(i),i=1,mx) / 1, 2, 3, 4, 5 /
c
      data nofz /
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
     $     0,0,0,0,0,0,0,0,0,0,0 /
c
c /czinte_opal_z/ data:
c
!!      data cxfil / 'Gx00', 'Gx03', 'Gx10', 'Gx35', 'Gx70' /
      data cxfil / 'Ax00', 'Ax03', 'Ax10', 'Ax35', 'Ax70' /
      data czfil / 'z0  ', 'z001', 'z004', 'z01 ',
     $     'z02 ', 'z03 ', 'z05 ', 'z10 '/
c
c /c_opal_ctrl_smooth/ data:
c
      data init_smo / 2 /, low_CO_smo / 1 /, interp_CO_smo / 1 /
c
c /opdir/ data:
c
      data copdir / ' ' /
c
c /d_opal_z/ data:
c
      data dkap / 0.0d0 /
c
c /c_level_err_opal_z/ data:
c
*****      data level_err / 1 /
      data level_err / 0 /
c
c /c_ctab/ data:
c
      data ctab / '	' /
c
c /chkpoc/ data:
c
c-vms[                                                !  For VMS:
c-vms;      data cb / ':', ']', ';', ';', ';', ';' /
c-vms]
c-sun-iris-linux[                                     ! For UNIX:
      data cb / '/', '/', '_', '~', '+', '-' /
c-sun-iris-linux]
c
      end
c
c******************************************************************************
c
      subroutine opalinit( khighz, ofebrack, z, kz )
c     ==============================================
c
c  INITIALIZATIONS AND OPACITY FILE SETUP:
c
c  This subroutine performs some initializations that would otherwise be done
c  at the beginning of subroutine READZEXCO.  These do some grid set-up, look
c  for the user-supplied non-zero [O/Fe] file if khighz = 5 or -5, calculate
c  [O/Fe] values for each of the possible mixes, and find the OPAL-opacity-file
c  directory name.
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      parameter ( nrdel=nrb-1, ntdel=ntb-1 )
c
c PARAMETERS:
c  zdel = 0.001 = offset for Z, Z+C, and Z+O, to make log interpolation behave
c                  reasonably at small Z values: Z-interpolation is performed
c                  using log(Z+zdel), while the CO-interpolation is performed
c                  using log(C+Z+zdel) and log(O+Z+zdel)
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
      parameter ( mx_hi=2*mx, mo_m1=mo-1, mx_hi_nz=mx_hi*nz )
c
      common /xhi_opal_z/ xhi_in(mx_hi), xcno_use(mx_hi,nz),
     $     xhi_use(mx_hi,nz), xxx_cno(mx_hi), xxx_hi(mx_hi),
     $     nx_hi(nz), ireq_hi(mx_hi), khighx(nz), kavail_xhi, kuse_xhi,
     $     kdo_xhi, kavail_cno, kuse_cno, kdo_cno, kavail_user,
     $     kuse_user, kdo_user
      save /xhi_opal_z/
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
      common/b_opal_z/ nta(0:nrm_p1),ntax0(0:nrm),
     $     ntax03(0:nrm), sltlo, slthi, dltlo_inv, dlthi_inv,
     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
      save /b_opal_z/
c
      common/e_opal_z/ opact,dopact,dopacr,dopactd,fedge,ftredge,fzedge
      save /e_opal_z/
c
      common /x_opal_z/ z_opal, x_opal, xc_opal, xo_opal, slt_opal,
     $     slr_opal, fcn_opal, fcon_opal, fcnone_opal, fu_opal
      save /x_opal_z/
c
      common/recoin_opal_z/ itimeco,mxzero,mx03,kope,igznotgx
      save /recoin_opal_z/
c
      common /alt_change_opal_z/ main_alt_change, iulow, khighz_in,
     $     ofebrack_in
      save /alt_change_opal_z/
c
!!      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
!!      parameter ( nel_zmix=19, n_zmixes=7, kel_o=3, kel_fe=nel_zmix-1 )
      parameter ( nel_zmix=19, n_zmixes=8, kel_o=3, kel_fe=nel_zmix-1 )
c
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
      parameter ( n_totmix = n_zmixes + 5, n_cnobeg = n_zmixes + 1 )
c
c-implicit;      real*4 bracketofe_opalGS98, xofe_opalGS98, xiz_opalGS98,
c-implicit;     $     fninz_opalGS98, atwt_opalGS98
      character*255 cfile_opalGS98(n_totmix)
      common /opalGS98mixes/ bracketofe_opalGS98(n_totmix),
     $     xofe_opalGS98(n_totmix),xiz_opalGS98(nel_zmix,n_totmix),
     $     fninz_opalGS98(nel_zmix,n_totmix),atwt_opalGS98(nel_zmix),
     $     cfile_opalGS98
      save /opalGS98mixes/
c
      character*8 cdef_CNO_ext(n_cnobeg:n_totmix)
      common /ext_CNO_opal_z/ len_def_CNO_ext(n_cnobeg:n_totmix),
     $     cdef_CNO_ext
      save /ext_CNO_opal_z/
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
      common/c_opal_ctrl_smooth/ init_smo, low_CO_smo, interp_CO_smo
      save /c_opal_ctrl_smooth/
c
      character*255 copdir
      common/opdir/ copdir
      save /opdir/
c
      character*1 cb(6)
      common/chkpoc/cb
      save /chkpoc/
c
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c___
c-dir;      logical lxst
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
c-debug]
c
c===
c					! These initializations only done once:
      if ( itimeco .ne. 12345678 ) then
c							! check some parameters
!!         if ( nrm .ne. 19 .or. ntm .ne. 70 ) stop
!!     $        ' STOP -- OPAL: NRM .ne. 19 or NTM .ne. 70 . '
         if ( nrm .ne. 39 .or. ntm .ne. 76 ) stop
     $        ' STOP -- OPAL: NRM .ne. 39 or NTM .ne. 76 . '
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
         do ix = 1, mx_hi
            xxx_cno(ix) = log10( xhi_in(ix) + xdel )
            xxx_hi(ix) = xxx_cno(ix)
         enddo
         xxx_hi(1) = log10( 0.03d0 + xdel )
c
c				! have now done once-only initializations
         itimeco = 12345678
c			     ! end of initializations that are done only once
      endif
c		   ! obtain the directory specification for the Gz???.x?? files
c
      kope = lnblnk(copdir)
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
      if ( kope .gt. 246 ) then
         write(6,10) copdir(:kope)
 10      format(
     $        ' STOP -- READCO: OPAL directory name > 246 characters:'/
     $        ' ',a)
         stop
      endif
c
c  NOTE that some systems return FALSE for the existence of a directory, so
c  one cannot check for the directory's existence.
c
c-dir;      if ( kope .gt. 0 ) then
c-dir;         call inqfil( copdir, lxst )
c-dir;         if ( .not. lxst ) then
c-dir;            write(6,20) copdir(:kope)
c-dir; 20         format(' STOP -- READCO: OPAL directory does not exist:'/
c-dir;     $           ' ',a)
c-dir;            stop
c-dir;         endif
c-dir;      endif
c		    ! just in case mx = 1 (i.e., if there is only one X-value)
      dfsx(2) = 1.d0
c
      itime = 0
      igznotgx = 0
      mxzero = 0
      mx03 = 0
c		     ! indices of X=0 and X=.03 mixes, and X-part of file names
      do i = 1, mx
c					! loop over X-index (i, not m, here!)
         if ( xa(i) .eq. 0.0d0 ) then
            mxzero = i
!!            cxfil(i) = 'Gx00'
            cxfil(i) = 'Ax00'
            mnofz(i) = 1
         else if ( abs(xa(i)-0.03d0) .lt. 1.d-6 ) then
            xa(i) = 0.03d0
            mx03 = i
!!            cxfil(i) = 'Gx03'
            cxfil(i) = 'Ax03'
            mnofz(i) = 2
         else if ( abs(xa(i)-0.1d0) .lt. 1.d-6 ) then
            xa(i) = 0.1d0
!!            cxfil(i) = 'Gx10'
            cxfil(i) = 'Ax10'
            mnofz(i) = 3
         else if ( abs(xa(i)-0.35d0) .lt. 1.d-6 ) then
            xa(i) = 0.35d0
!!            cxfil(i) = 'Gx35'
            cxfil(i) = 'Ax35'
            mnofz(i) = 4
         else if ( abs(xa(i)-0.7d0) .lt. 1.d-6 ) then
            xa(i) = 0.7d0
!!            cxfil(i) = 'Gx70'
            cxfil(i) = 'Ax70'
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
      if ( z .eq. 0.d0 ) then
         khizat = 0
         klozat = 0
      else if ( khighz .lt. 0 ) then
         khizat = 1
         if ( ofebrack .eq. 0.d0 ) then
            klozat = 1
         else
            klozat = min( mod( abs(khighz), 10 ) , n_zmixes )
         endif
      else
         klozat = 0
         khizat = min( mod( khighz, 10 ) , n_zmixes )
         if ( ofebrack .eq. 0.d0 ) khizat = min(khizat,1)
         if ( khizat .eq. 1 .and.
     $        ( ( z .ge. 0.01d0 .and. z .le. 0.02d0 ) .or.
     $        ( abs(zval(kzbelo)-z) .le. 1.d-6 .and. z .ge. 1.d-5 ) ) )
     $        khizat = 0
      endif
c				! check length of GS98 filenames to be used
      if ( klozat .gt. 0 ) then
         do k = 1, klozat, max( klozat - 1 , 1 )
            if ( kope + lnblnk( cfile_opalGS98(k) ) .gt. 255 ) stop
     $           ' STOP -- READCO: GS98 file name too long. '
         enddo
      endif
c					! should use the input [O/Fe] filename?
      if ( khizat .ge. n_zmixes ) then
c								! it exists?
         if ( cfile_opalmixes(n_zmixes) .eq. '        ' ) stop
     $        ' STOP -- READCO: no user-specified [O/Fe]-file. '
c
c			 ! obtain mix specifications for the input [O/Fe] file
         igetzxi = 9
         i_rewind = 0
         itab_dum = 0
         line = 0
c						     ! use copdir temporarily
         copdir(kope+1:) = cfile_opalmixes(n_zmixes)
         call open_chk_zip( iulow, copdir, i_gzip,
     $        'STOP -- Error: user-specified [O/Fe]-file not found.' )
         ifound = mixfind(iulow,n_zmixes,igetzxi,i_rewind,itab_dum,
     $        line,0.0d0,0.0d0,0.0d0,0.0d0)
         if ( ifound .eq. 0 ) stop
     $        ' STOP -- READCO: bad user-specified [O/Fe]-file. '
         call close_chk_zip( iulow, copdir, i_gzip )
c
c					  ! remove filename from directory name
         copdir(kope+1:) = ' '
c					   ! or use GS98 input [O/Fe] filename?
      else if ( klozat .ge. n_zmixes ) then
c								! it exists?
         if ( cfile_opalGS98(n_zmixes) .eq. ' ' ) stop
     $        ' STOP -- READCO: no user-specified GS98 [O/Fe]-file. '
c								      ! length?
         if ( kope + lnblnk( cfile_opalGS98(n_zmixes) ) .gt. 255 ) stop
     $        ' STOP -- READCO: user GS98 [O/Fe]-file name too long. '
c
c			 ! obtain mix specifications for the input [O/Fe] file
         igetzxi = 9
         i_rewind = 0
         itab_dum = 0
         line = 0
c						       ! use copdir temporarily
         copdir(kope+1:) = cfile_opalGS98(n_zmixes)
         call open_chk_zip( iulow, copdir, i_gzip,
     $        'STOP -- Error: user GS98 [O/Fe]-file not found.' )
         ifound = mixfind(iulow,-n_zmixes,igetzxi,i_rewind,itab_dum,
     $        line,0.0d0,0.0d0,0.0d0,0.0d0)
         if ( ifound .eq. 0 ) stop
     $        ' STOP -- READCO: bad user-specified GS98 [O/Fe]-file. '
         call close_chk_zip( iulow, copdir, i_gzip )
c
c					  ! remove filename from directory name
         copdir(kope+1:) = ' '
c
      endif
c							     ! changed 'GS98hz'
      if ( khighz .lt. 0 .and. main_alt_change .gt. 0 ) then
c								! it exists?
         if ( cfile_opalGS98(1) .eq. ' ' ) stop
     $        ' STOP -- READCO: no main alternate [O/Fe]=0.0 file. '
c								      ! length?
         if ( kope + lnblnk( cfile_opalGS98(1) ) .gt. 255 ) stop
     $        ' STOP -- READCO: alternate [O/Fe]=0.0 name too long. '
c
c			 ! obtain mix specifications for input [O/Fe]=0 file
         igetzxi = 9
         i_rewind = 0
         itab_dum = 0
         line = 0
c						       ! use copdir temporarily
         copdir(kope+1:) = cfile_opalGS98(1)
         call open_chk_zip( iulow, copdir, i_gzip,
     $        'STOP -- Error: alternate [O/Fe]=0.0 file not found.' )
         ifound = mixfind(iulow,-1,igetzxi,i_rewind,itab_dum,
     $        line,0.0d0,0.0d0,0.0d0,0.0d0)
         if ( ifound .eq. 0 ) stop
     $        ' STOP -- READCO: bad alternate [O/Fe]=0.0 file. '
         call close_chk_zip( iulow, copdir, i_gzip )
c							! have read it now
         main_alt_change = main_alt_change - 2
c
c					  ! remove filename from directory name
         copdir(kope+1:) = ' '
c
      endif
c				! get mix Z-composition specifications (these
c				! will be recomputed for any mix read in later)
      do i = 1, n_zmixes
         xofe_opalGS98(i) = fninz_opalGS98(kel_o,i)
     $        / max( fninz_opalGS98(kel_fe,i) , 1.d-36 )
         bracketofe_opalGS98(i) = log10( xofe_opalGS98(i)
     $        / xofe_opalGS98(1) )
         xofe_opalmixes(i) = fninz_opalmixes(kel_o,i)
     $        / max( fninz_opalmixes(kel_fe,i) , 1.d-36 )
         bracketofe_opalmixes(i) = log10( xofe_opalmixes(i)
     $        / xofe_opalmixes(1) )
      enddo
c			       ! Reset current-mix data.  If GS98 [O/Fe] shift:
      if ( klozat .gt. 1 ) then
c				! get interpolation factors fofe (for GS98hz) &
c						  ! omfofe=1-fofe (other file)
         xofe = 10.d0**ofebrack * xofe_opalGS98(1)
         fofe = ( fninz_opalGS98(kel_o,klozat)
     $        - xofe * fninz_opalGS98(kel_fe,klozat) )
     $        / ( ( fninz_opalGS98(kel_fe,1)
     $        - fninz_opalGS98(kel_fe,klozat) ) * xofe
     $        + fninz_opalGS98(kel_o,klozat)
     $        - fninz_opalGS98(kel_o,1) )
         omfofe = 1.d0 - fofe
c					! get Z-composition of interpolated mix
         sum_niai = 0.0d0
         do i = 1, nel_zmix
            fninz_mix(i) = fofe * fninz_opalGS98(i,1)
     $           + omfofe * fninz_opalGS98(i,klozat)
            xiz_mix(i) = fninz_mix(i) * atwt_opalGS98(i)
            sum_niai = sum_niai + xiz_mix(i)
         enddo
         do i = 1, nel_zmix
            xiz_mix(i) = xiz_mix(i) / sum_niai
            bracketife_mix(i) = log10( ( max( fninz_mix(i) , 1.d-36 )
     $           * fninz_opalGS98(kel_fe,1) )
     $           / ( max( fninz_mix(kel_fe) , 1.d-36 )
     $           * fninz_opalGS98(i,1) ) )
         enddo
c				      ! Else, if use GS98 but no [O/Fe] shift:
      else if ( khighz .lt. 0 ) then
c
         do i = 1,nel_zmix
            xiz_mix(i) = xiz_opalGS98(i,1)
            fninz_mix(i) = fninz_opalGS98(i,1)
            bracketife_mix(i) = 0.d0
         enddo
c				      ! Else, if there is no GN93 [O/Fe] shift:
      else if ( khizat .le. 1 ) then
c
         do i = 1,nel_zmix
            xiz_mix(i) = xiz_opalmixes(i,1)
            fninz_mix(i) = fninz_opalmixes(i,1)
            bracketife_mix(i) = 0.d0
         enddo
c		  ! Else, if there is the [O/Fe] shift (also done in READEXCO):
      else
c		! get interpolation factors fofe (for GN93hz) and omfofe=1-fofe
c
         xofe = 10.d0**ofebrack * xofe_opalmixes(1)
         fofe = ( fninz_opalmixes(kel_o,khizat)
     $        - xofe * fninz_opalmixes(kel_fe,khizat) )
     $        / ( ( fninz_opalmixes(kel_fe,1)
     $        - fninz_opalmixes(kel_fe,khizat) ) * xofe
     $        + fninz_opalmixes(kel_o,khizat)
     $        - fninz_opalmixes(kel_o,1) )
         omfofe = 1.d0 - fofe
c					! get Z-composition of interpolated mix
         sum_niai = 0.0d0
         do i = 1, nel_zmix
            fninz_mix(i) = fofe * fninz_opalmixes(i,1)
     $           + omfofe * fninz_opalmixes(i,khizat)
            xiz_mix(i) = fninz_mix(i) * atwt_opalGS98(i)
            sum_niai = sum_niai + xiz_mix(i)
         enddo
         do i = 1, nel_zmix
            xiz_mix(i) = xiz_mix(i) / sum_niai
            bracketife_mix(i) = log10( ( max( fninz_mix(i) , 1.d-36 )
     $           * fninz_opalmixes(kel_fe,1) )
     $           / ( max( fninz_mix(kel_fe) , 1.d-36 )
     $           * fninz_opalmixes(i,1) ) )
         enddo
      endif
c					! end of initializations
      return
      end
c
c******************************************************************************
c
      subroutine get_zavail
c     =====================
c
c  Obtain combined Z-tabulation values available in the files, and best order
c  in which to enhance the 'Gz???.x??' metallicity table.
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
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
               if ( zval(k_l) .lt. zavail(k_a) + 1.e-6 ) then
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
                     if ( d_z .gt. 1. ) d_z = 1. / d_z
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
c******************************************************************************
c
      subroutine get_trvals
c     =====================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      parameter ( nrdel=nrb-1, ntdel=ntb-1 )
c
c PARAMETERS: positions of DlogT-change in table, low logT and logR values:
c
      parameter ( ks81=ntm-3, ks83=ks81+1, ks60=ks81-21, ks61=ks60+1,
     $     alrlo=-9.0d0, flogtlo=3.75d0, flogt60=6.6d0, flogt81=9.3d0 )
!!     $     alrlo=-8.0, flogtlo=3.75, flogt60=6.0, flogt81=8.1 )
c
c COMMON /a_opal_z/ : matrices for opacity storage and interpolation:
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
c COMMON /b_opal_z/ : high and low logT6 and logR limits, and mix Z-values:
c
      common/b_opal_z/ nta(0:nrm_p1),ntax0(0:nrm),
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
         if ( i .le. 46 ) then
            flogtin(i) = (i-1)*0.05d0+flogtlo
         else if ( i .gt. 46 ) then
            flogtin(i) = flogtin(i-1)+0.1d0
         endif            
         if ( i .ge. ntb ) dfs(i-ntdel) = 20.d0
      enddo
      if ( abs(flogtin(ks60)-flogt60) .gt. 5.d-6 )
     &   stop ' STOP -- READCO: initialization error. '
      flogtin(ks60) = flogt60
      do i = ks61,ks81
         if ( i .le. 67) then
            flogtin(i) = (i-ks60)*0.1d0+flogt60
         else
            flogtin(i) = flogtin(i-1)+0.2d0
         endif
         if ( i .ge. ntb ) dfs(i-ntdel) = 10.d0
      enddo
      if ( abs(flogtin(ks81)-flogt81) .gt. 5.d-6 ) then 
         print *,'abs(flogtin(ks81)-flogt81)',flogtin(ks81),flogt81
         stop ' STOP -- READCO: initialization error. '
      endif
      flogtin(ks81) = flogt81
      do i = ks83,ntm
         flogtin(i) = (i-ks81)*0.2d0+flogt81
         if ( i .ge. ntb ) dfs(i-ntdel) = 5.d0
      enddo
      do i=1,ntm
         t6arr(i) = 10.d0**(flogtin(i)-6.0d0)
      enddo
      do i=1,nt
         alt(i) = flogtin(i+ntdel)-6.0d0
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
c						       ! R-extrapolation limits
      slrlo = alr(1)
      slrhi = alr(nr)
c				     ! make 1-grid-pt extrap just within limits
      dlrlo_inv = dfsr(1) * 0.999999d0
      dlrhi_inv = dfsr(nr) * 0.999999d0
c						       ! T-extrapolation limits
      sltlo = alt(1)
      slthi = alt(nt)
c				     ! make 1-grid-pt extrap just within limits
      dltlo_inv = dfs(1) * 0.999999d0
      dlthi_inv = dfs(nt) * 0.999999d0
c
      return
      end
c
c******************************************************************************
c
      subroutine z_fcno(xh,xmet,nmet,fu,z,xci,xoi,fcn,fcon,fcnone,fuse)
c     =================================================================
c
c.....Given  XH, XMET(NMET), and FU  : calculates and returns the values of
c     Z, XCI, XOI, FCN, FCON, FCNONE, FUSE  (used by OPAL_X_CNO_FU below).
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      dimension xmet(nmet)
c
!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      parameter ( mx_hi=2*mx, mo_m1=mo-1, mo_m2=mo-2 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
      common /xhi_opal_z/ xhi_in(mx_hi), xcno_use(mx_hi,nz),
     $     xhi_use(mx_hi,nz), xxx_cno(mx_hi), xxx_hi(mx_hi),
     $     nx_hi(nz), ireq_hi(mx_hi), khighx(nz), kavail_xhi, kuse_xhi,
     $     kdo_xhi, kavail_cno, kuse_cno, kdo_cno, kavail_user,
     $     kuse_user, kdo_user
      save /xhi_opal_z/
c
!!      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
!!      parameter ( nel_zmix=19, n_zmixes=7, kel_o=3, kel_fe=nel_zmix-1 )
      parameter ( nel_zmix=19, n_zmixes=8, kel_o=3, kel_fe=nel_zmix-1 )
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
      parameter ( n_totmix = n_zmixes + 5, n_cnobeg = n_zmixes + 1 )
c
      character*255 cfile_opalGS98(n_totmix)
      common /opalGS98mixes/ bracketofe_opalGS98(n_totmix),
     $     xofe_opalGS98(n_totmix),xiz_opalGS98(nel_zmix,n_totmix),
     $     fninz_opalGS98(nel_zmix,n_totmix),atwt_opalGS98(nel_zmix),
     $     cfile_opalGS98
      save /opalGS98mixes/
c
      common /cno_delta_opal_z/ fcno_mul(4), fninz_cno(nel_zmix,5),
     $     xiz_cno(nel_zmix,5), d_fninz_user(nel_zmix),
     $     fcno_fac(0:3,4), fninz_heavy, xiz_heavy, d_fninz_u_heavy,
     $     s_ninzai_mix, ds_ninzai_u, fn_o_over_cno, fninz_co_mix
      save /cno_delta_opal_z/
c
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c___
      dimension fninz_u(nel_zmix), xiz_u(nel_zmix), fninz(nel_zmix)
c
c PARAMETERS: for use in call to READZEXCO, if it is called from here:
c  k_hz = 1 = khighz value for READZEXCO ("use 'GN93hz' in Z-interpolation,
c               but no CNO-interpolation")
c  ofe_brack = 0.0 = [O/Fe] value for READZEXCO
c
      parameter ( k_hz=1, ofe_brack=0.0d0 )
c===
c
c  Perform checks on the inputs, and compute total mass Ztot of metals.
c
      if ( nmet .lt. 5 ) then
         write(6,10) nmet
 10      format(' '/' STOP -- OPAL_X_CNO_FU: metals array',
     $        ' Xmet(Nmet): size too small: Nmet =',i2,' < 5')
         stop ' STOP -- OPAL: bad composition (Xmet array size). '
      else if ( level_err .gt. 2 .and.
     $        nmet .ne. nel_zmix .and. kdo_cno .gt. 0 ) then
         write(6,15) nmet, nel_zmix
 15      format(' '/' STOP -- OPAL_X_CNO_FU: metals array',
     $        ' Xmet(Nmet): Nmet =',i5,', but nel_zmix =',i3)
         stop ' STOP -- OPAL: bad composition (Xmet array size). '
      endif
c
      ximin = min( xh , xmet(1) , xmet(2) , xmet(3) , xmet(4) )
      Zheavy = 0.0d0
      do i = 5, nmet
         Zheavy = Zheavy + xmet(i)
         ximin = min( ximin , xmet(i) )
      enddo
      Ztot = Zheavy + xmet(1) + xmet(2) + xmet(3) + xmet(4)
c
      if ( xh + Ztot - 1.d-6 .gt. 1.0d0 .or. Zheavy .lt. -1.d-8 .or.
     $     min( Ztot , ximin ) .lt. -1.d-6 ) then
         write(6,20) xh, Ztot
 20      format(' '/' STOP -- OPAL_X_CNO_FU: bad X',
     $        f10.7,' Ztot',f10.7,' or Zheavy',f11.8)
         stop ' STOP -- OPAL: bad composition. '
      endif
c
c..... If necessary, read data files and do initializations, using READZEXCO
c      (this performs initializations necessary for calculating quantities
c      further below); note that xiz_mix is initialized by a data statement
c      to 'GN93hz' relative-metal-abundance values
c
      if ( itime .ne. 12345678 ) then
c					 ! first get approximate Z using the
c					 ! 'GN93hz' relative metal abundances
         if ( Ztot .le. 1.d-6 ) then
            z = max( Ztot , 0.0d0 )
         else
            xiz_heavy_u = 0.0d0
            do i = 5, nel_zmix
               xiz_heavy_u = xiz_heavy_u + xiz_mix(i)
            enddo
            z = max( Zheavy / xiz_heavy_u , 1.d-6 )
         endif
c		       ! read opacities: use 'GN93hz', but no CNO-interpolation
c							    ! will be possible:
         call readzexco(nz,-1.,z,-1.d0,k_hz,-1,ofe_brack)
c
      endif
c
c  If there are essentially no metals (except perhaps C,N,O,Ne), there is no
c  point in trying to do anything fancy with interpolation in the make-up of Z.
c
      if ( Ztot .le. 1.d-6 .or. Zheavy .le. 1.d-8 ) then
c
         fuse = 0.0d0
         fcn = 0.0d0
         fcon = 0.0d0
         fcnone = 0.0d0
c				      ! for very small Z, this is good enough:
         if ( Ztot .le. 1.e-6 ) then
            z = max( Ztot , 0.0d0 )
         else
            z = max( Zheavy / xiz_heavy , 1.d-6 )
         endif
c					    ! any excess C,N,O,Ne --> XCI, XOI:
         xci = xmet(1) - z * xiz_mix(1)
     $        + 0.5d0 * ( xmet(2) - z * xiz_mix(2) )
         xoi = xmet(3) + xmet(4)
     $        - z * ( xiz_mix(3) + xiz_mix(4) )
     $        + 0.5d0 * ( xmet(2) - z * xiz_mix(2) )
c
c  OTHERWISE: if the total amount of metals is not insignificant:
c
      else
c
c  Obtain the OPAL reference mix (relative to Z) in both number and mass
c  fractions, and the fraction of these that correspond to elements heavier
c  than Ne; note that this reference mix can be affected by the composition
c  of the user-specified opacity-shift file, if this is being used.
c
         if ( kdo_user .le. 0 .or. abs(fu) .lt. 1.d-6 ) then
c								! no fu effect:
            fuse = 0.0d0
            do i = 1, nel_zmix
               fninz_u(i) = fninz_mix(i)
               xiz_u(i) = xiz_mix(i)
            enddo
            fninz_heavy_u = fninz_heavy
            xiz_heavy_u = xiz_heavy
c				   ! else: include fu effect, if it is non-zero
         else
c			! don't allow too much heavy-element reduction from fu
            fuse = fu
            fninz_heavy_u = fninz_heavy + fuse * d_fninz_u_heavy
            if ( fninz_heavy_u .lt. 0.5d0 * fninz_heavy ) then
               fuse = fu * 0.5d0 * fninz_heavy
     $              / ( fninz_heavy - fninz_heavy_u )
               fninz_heavy_u = 0.5d0 * fninz_heavy
            endif
c				! get shifted composition using fu-factor fuse
            sum_niai = 0.0d0
            do i = 1, nel_zmix
               fninz_u(i) = fninz_mix(i) + fuse * d_fninz_user(i)
               xiz_u(i) = fninz_u(i) * atwt_opalGS98(i)
               sum_niai = sum_niai + xiz_u(i)
            enddo
c
            xiz_u(1) = xiz_u(1) / sum_niai
            xiz_u(2) = xiz_u(2) / sum_niai
            xiz_u(3) = xiz_u(3) / sum_niai
            xiz_u(4) = xiz_u(4) / sum_niai
c
            xiz_heavy_u = 0.0d0
            do i = 5, nel_zmix
               xiz_u(i) = xiz_u(i) / sum_niai
               xiz_heavy_u = xiz_heavy_u + xiz_u(i)
            enddo
c
         endif
c
c  If CNO-interpolation is not available, then all we have to do now is to
c  compute the metallicity Z and the excess C and O amounts:
c
         if ( kdo_cno .le. 0 ) then
c					! no CNO-interpolation: factors are 0.0
            fcn = 0.0d0
            fcon = 0.0d0
            fcnone = 0.0d0
c						  ! heavies give metallicity Z
            z = max( Zheavy / xiz_heavy_u , 0.0d0 )
c						 ! excess C,N,O,Ne --> XCI, XOI
            xci = xmet(1) - z * xiz_mix(1)
     $           + 0.5d0 * ( xmet(2) - z * xiz_mix(2) )
            xoi = xmet(3) + xmet(4)
     $           - z * ( xiz_mix(3) + xiz_mix(4) )
     $           + 0.5d0 * ( xmet(2) - z * xiz_mix(2) )
c
c  ELSE, if CNO-interpolation is available: determine the CNO-interpolation
c  factors as well as the metallicity Z and the excess C and O amounts:
c
         else
c
c.....Get the factor f_nz that converts fractions of Z to the actual mass
c     fractions, and use it to convert the input mass fractions to adjusted
c     number fractions of Z (the CNO-mix) --- if the user-fraction  fu  is
c     non-zero, subtract off the effects of the difference in composition
c     between the user-opacity-file and the OPAL mix (multiplied by fu).
c
c  NOTE: for OPAL-mix(Z): SUM_z{ Xmix_i } = Z, SUM_z{ fninz_mix(i) } = 1.0,
c
c                 SUM_h{ Xmix_i } = Zheavy, SUM_h{ fninz_mix(i) } = fninz_heavy
c
c             where  fninz_mix(i) = ( Xmix_i / A_i ) / SUM_z{ Xmix_j / A_j } so
c
c                fninz_heavy = SUM_h{ Xmix_i / A_i } / SUM_z{ Xmix_j / A_j } .
c
c					    ! If user supplies OPAL-element Xi:
            if ( nmet .eq. nel_zmix ) then
c
c        For input-mix with same number of elements as OPAL-mix, assume that
c                  SUM_h{ fninz(i) } = fninz_heavy  is what determines Z, i.e.,
c                    fninz_heavy = SUM_h{ X_i / A_i } / SUM_z{ Xmix_j / A_j }
c                  i.e., assume  SUM_h{ X_i / A_i } = SUM_h{ Xmix_i / A_i } .
c             But  SUM_z{ fninz(i) }  is NOT necessarily equal to 1.0, since
c             there may be excess C and O.  We therefore set
c                f_nz = SUM_h{ X_i / A_i } / fninz_heavy
c                     = SUM_h{ X_i / A_i }
c                          / [ SUM_h{ Xmix_i / A_i } / SUM_z{ Xmix_j / A_j } ]
c                     = SUM_z{ Xmix_j / A_j }
c             due to our assumption  SUM_h{ X_i / A_i } = SUM_h{ Xmix_i / A_i }
c             and thus we set
c                fninz(i) = X_i / ( f_nz * A_i )
c                         = ( X_i / A_i ) / SUM_z{ Xmix_i / A_i }
c             which is what we want.
c						! get metallicity factor f_nz
               f_nz = 0.0d0
               do i = 5, nel_zmix
                  f_nz = f_nz + xmet(i) / atwt_opalGS98(i)
               enddo
               f_nz = f_nz / fninz_heavy_u
c						! and use it to get CNO-mix
               do i = 1, nel_zmix
                  fninz(i) = xmet(i) / ( f_nz * atwt_opalGS98(i) )
     $                 - fuse * d_fninz_user(i)
               enddo
c					     ! Else: non-OPAL-Xi, but fu = 0.0:
            else if ( fuse .eq. 0.0d0 ) then
c
c       For input-mix with different number of elements from OPAL-mix, assume
c                 SUM_h{ X_i } = SUM_h{ Xmix_i } = Zheavy .
c             But  xiz_heavy = SUM_h{ Xmix_i / Z } = Zheavy / Z , and
c                  s_ninzai_mix = SUM_z{ fninz_mix(i) * A_i }
c                               = SUM_z{ Xmix_i } / SUM_z{ Xmix_j / A_j }
c                               = Z / SUM_z{ Xmix_j / A_j }
c             We therefore set
c                f_nz = Zheavy / ( xiz_heavy * s_ninzai_mix )
c                     = Zheavy
c                          / ( [ Zheavy / Z ] * [ Z / SUM_z{ Xmix_j / A_j } ] )
c                     = SUM_z{ Xmix_j / A_j }
c             and thus for the non-heavies we set
c                fninz(i) = X_i / ( f_nz * A_i )
c                         = ( X_i / A_i ) / SUM_z{ Xmix_i / A_i }
c             which is what we want (for the heavies, we assume that they are
c             distributed as in the OPAL-mix).
c								! get f_nz
               f_nz = Zheavy / ( xiz_heavy_u * s_ninzai_mix )
c								! and CNO-mix
               do i = 1, 4
                  fninz(i) = xmet(i) / ( f_nz * atwt_opalGS98(i) )
               enddo
               do i = 5, nel_zmix
                  fninz(i) = fninz_mix(i)
               enddo
c					   ! Else: non-OPAL-Xi and non-zero fu:
            else
c								! get f_nz
               f_nz = Zheavy / ( xiz_heavy_u
     $              * ( s_ninzai_mix + fuse * ds_ninzai_u ) )
c								! and CNO-mix
               do i = 1, 4
                  fninz(i) = xmet(i) / ( f_nz * atwt_opalGS98(i) )
     $                 - fuse * d_fninz_user(i)
               enddo
               do i = 5, nel_zmix
                  fninz(i) = fninz_mix(i)
               enddo
c
            endif
c
c.....Get the number fraction differences between the adjusted CNO-mix and the
c     original OPAL-mix.  Any excess in the sum of N + Ne abundances must come
c     from C and/or O depletion, so use the former to get these latter amounts.
c     Note that any Ne increase is expected to come from both C and O, so set
c     a preliminary O depletion to be the amount of the Ne increase multiplied
c     by the initial ratio of O to the sum of C,N,O.  If this is larger than
c     the initial O abundance, reduce the excess O depletion (i.e., the amount
c     of O-depletion beyond the initial O abundance) by a factor of 10 (this
c     is equivalent to assuming that the third dredge-up yields an order of
c     magnitude larger C increase than the O increase).  Most or all of the
c     difference between the total N + Ne increase and this (prelimiary) O
c     depletion should be assigned to C depletion (but the O depletion may need
c     a slight readjustment).  Adjust the excess C and O amounts according to
c     these C and O depletions; but if the total excess CO is positive, do not
c     allow either of the C excess or the O excess to be negative --- note that
c     this condition may cause a readjustment of the C and O depletions.
c
c					   ! get CNO-mix - OPAL-mix differences
            dn_c = fninz(1) - fninz_mix(1)
            dn_n = fninz(2) - fninz_mix(2)
            dn_o = fninz(3) - fninz_mix(3)
            dn_ne = fninz(4) - fninz_mix(4)
c					    ! total increase in N + Ne
            deln_nne = dn_ne + dn_n
c						! total excess CO (by number)
            fn_co_ex = dn_c + dn_o + deln_nne
c						! If have negligible excess CO:
            if ( fn_co_ex .lt. 1.e-4 ) then
c						! just divide it equally in C,O
               fn_c_ex = 0.5d0 * fn_co_ex
               fn_o_ex = fn_c_ex
               deln_c = fn_c_ex - dn_c
               deln_o = fn_o_ex - dn_o
c						! Else: have some excess CO:
            else
c						 ! If negligible increase N+Ne:
               if ( deln_nne .lt. 1.d-5 ) then
c						! just divide it equally in C,O
                  deln_c = 0.5d0 * deln_nne
                  deln_o = deln_c
c				       ! Else: both C,O --> N,Ne and excess CO:
               else
c					     ! any Ne increase is partly from O
                  if ( dn_ne .gt. 0.0d0 ) then
                     deln_o = fn_o_over_cno * dn_ne
                     if ( deln_o .gt. fninz_mix(3) )
     $                    deln_o = fninz_mix(3)
     $                    + 0.1d0 * ( deln_o - fninz_mix(3) )
                     deln_rem = deln_nne - deln_o
                  else
                     deln_o = 0.0d0
                     deln_rem = deln_nne
                  endif
c						     ! If not much (more) N+Ne:
c
                  if ( deln_rem .le. fninz_mix(1) ) then
c							 ! assign it to C
                     deln_c = deln_rem
c				      ! Else if N+Ne increase is < initial C+O:
c
                  else if ( deln_nne .le. fninz_co_mix ) then
c								    ! divide it
                     deln_o = max( deln_nne - fninz_mix(1) , deln_o )
                     deln_c = deln_nne - deln_o
c						 ! Else if large N+Ne increase:
                  else
c							   ! most is from C
                     deln_o = max( fninz_mix(3) , deln_o )
                     deln_c = deln_nne - deln_o
c
                  endif
c
               endif
c				! adjust the excess C and O amounts according
c					! to the above C and O depletions
               fn_c_ex = dn_c + deln_c
               fn_o_ex = dn_o + deln_o
c					! but total excess CO is > 0, so do not
c						! allow negative excess C or O
               if ( fn_c_ex .lt. 0.0d0 ) then
                  deln_c = deln_c - fn_c_ex
                  deln_o = deln_o + fn_c_ex
                  fn_o_ex = fn_o_ex + fn_c_ex
                  fn_c_ex = 0.0d0
               else if ( fn_o_ex .lt. 0.0d0 ) then
                  deln_c = deln_c + fn_o_ex
                  deln_o = deln_o - fn_o_ex
                  fn_c_ex = fn_c_ex + fn_o_ex
                  fn_o_ex = 0.0d0
               endif
c
            endif
c
c.....The excess C and O are not part of Z: subtract them off, then use the
c     previously-computed  f_nz  factor to calculate the metallicity Z.  Obtain
c     the excess C and O mass fractions XCI and XOI from the excess C and O
c     number fractions.  Use pre-calculated interpolation factors to get the
c     CNO-interpolation factors FCN, FCON, FCNONE from the C, N, and O number
c     fractions (relative to Z) that have just been computed (if all three of
c     these CNO-interpolation factors are very small, set them to zero).
c
c					  ! subtract off excess C, O by number
            fninz(1) = fninz(1) - fn_c_ex
            fninz(3) = fninz(3) - fn_o_ex
c					  ! compute metallicity mass fraction Z
            z = 0.0d0
            do i = 1, nel_zmix
               z = z + ( fninz(i) + fuse * d_fninz_user(i) )
     $              * atwt_opalGS98(i)
            enddo
            z = max( z * f_nz , 0.0d0 )
c						    ! excess C,O mass fractions
            xci = fn_c_ex * atwt_opalGS98(1) * f_nz
            xoi = fn_o_ex * atwt_opalGS98(3) * f_nz
c							     ! for CNO-interp:
            fcn = ( fcno_fac(0,2) + fcno_fac(1,2) * fninz(1)
     $           + fcno_fac(2,2) * fninz(2)
     $           + fcno_fac(3,2) * fninz(3) ) * fcno_mul(2)
            fcon = ( fcno_fac(0,3) + fcno_fac(1,3) * fninz(1)
     $           + fcno_fac(2,3) * fninz(2)
     $           + fcno_fac(3,3) * fninz(3) ) * fcno_mul(3)
            fcnone = ( fcno_fac(0,4) + fcno_fac(1,4) * fninz(1)
     $           + fcno_fac(2,4) * fninz(2)
     $           + fcno_fac(3,4) * fninz(3) ) * fcno_mul(4)
c								! very small?
            if ( max( abs(fcn) , abs(fcon) , abs(fcnone) ) .lt.
     $           1.d-5 .or. z .lt. 1.d-6 ) then
               fcn = 0.0d0
               fcon = 0.0d0
               fcnone = 0.0d0
            endif
c
         endif
c
      endif
c
c  Check for xci, xoi too large or too negative.
c
      del_sum = ( xh + z + xci + xoi - 1.0d0 ) * 0.5d0
      if ( del_sum .gt. 0.0d0 ) then
         xci = xci - del_sum
         xoi = xoi - del_sum
      endif
c
      if ( xci .lt. -0.5d0 * z ) then
         xoi = max( -0.5d0 * z , xoi + ( xci + 0.5d0 * z ) )
         xci = -0.5d0 * z
      else if ( xoi .lt. -0.5d0 * z ) then
         xci = max( -0.5d0 * z , xci + ( xoi + 0.5d0 * z ) )
         xoi = -0.5d0 * z
      endif
c
      return
      end
c
c******************************************************************************
c
      subroutine opac(z,xh,xci,xoi,t6,r,error)
c     ==================================
c
c.....This is just an alternate interface to OPAL_F_CNOU below, which it calls
c     after taking the log of T6 and R and setting "NO CNO/user interpolation";
c        temperature-input  T6 = temperature in millions of degrees kelvin
c        density-parameter-input  R = density(g/cm**3) / T6**3
c===
      implicit none

      double precision z,xh,xci,xoi,t6,r,slt,slr

      integer error

      error = 0

      if ( t6 .le. 0.d0 .or. r .le. 0.d0 ) then
         write(6,8437) t6,r
 8437    format(' '/' STOP -- OPAC: non-positive value of T6=',
     $        1p,e11.3,' or R=',e11.3)
         stop
      endif
c
      slt = log10(t6)
      slr = log10(r)
c
c      print *, 'z',z,'x',xh,'xci',xci,'xoi',xoi
      call opal_f_cnou(z,xh,xci,xoi,slt,slr,0.0d0,0.0d0,0.0d0,0.0d0)
c      print *, 'z ap',z,'x ap',xh,'xci ap',xci,'xoi ap',xoi
c
      return
      end
c
c******************************************************************************
c
      subroutine opal(z,xh,xci,xoi,slt,slr,error)
c     =====================================
c
c.....This is just an alternate interface to OPAL_F_CNOU below, which it calls
c     after setting "NO CNO/user interpolation":
c===
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      integer error
      common/extopac/ opact,dopact,dopacr,dopactd,fedge,ftredge,fzedge

      call opal_f_cnou(z,xh,xci,xoi,slt,slr,0.0d0,0.0d0,0.0d0,0.0d0)
c
      return
      end
c
c******************************************************************************
c
      subroutine opal_x_cno_fu(xh,slt,slr,xmet,nmet,fu)
c     =================================================
c
c.....This is an alternate interface to OPAL_F_CNOU below
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      dimension xmet(nmet)
c===
c
c  Get the CNO-interpolation factors, the metallicity Z, and the excess C,O:
c
      call z_fcno(xh,xmet,nmet,fu,z,xci,xoi,fcn,fcon,fcnone,fuse)
c
c  Get the opacities by calling OPAL_F_CNOU:
c
      call opal_f_cnou(z,xh,xci,xoi,slt,slr,fcn,fcon,fcnone,fuse)
c
      return
      end
c
c******************************************************************************
c
      subroutine opal_f_cnou(z,xh,xci,xoi,slt,slr,fcn,fcon,fcnone,fu)
c     ===============================================================
c
c..... The purpose of this subroutine is to interpolate log kappa
c          (and obtain smooth derivatives)
c      in hydrogen (if X>0) and in C/O abundance and T6, R, i.e. (X,Xc,Xo,T6,R)
c      Interpolation in CNO abundances is allowed by user-specified fractions.
c
c      z = Z = metallicity (this should always be the same)
c      xh = X = hydrogen mass fraction
c      xci = Xc = carbon mass fraction (excess over what is in Z)
c      xoi = Xo = oxygen mass fraction (excess over what is in Z)
c      slt = logT6 = Log10{temperature in millions of degrees kelvin}
c      slr = logR = Log10{density(g/cm**3)/T6**3}
c      fcn = fraction to be applied of opacity shift from standard composition
c              to a composition with most or all of C converted to N
c      fcon = fraction to be applied of opacity shift from standard composition
c               to a composition with most or all of C and O converted to N
c      fcnone = fraction to be applied of the opacity shift from the standard
c                 composition to a composition with all C,N,O converted to Ne
c      fu = fraction to be applied of opacity shift from standard composition
c             to the composition of the user-specified opacity file
c
c..... to use OPAL_F_CNOU, insert common/e_opal_z/ in the calling routine.
c      This common contains interpolated values for log kappa and its 
c      first derivatives, and "out-of-table" indicators fedge,ftredge,fzedge.
c
c PARAMETERS to specify opacity storage matrices: see OPALINIT
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
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
c  k_hz = 1 = khighz value for READZEXCO ("use 'GN93hz' in Z-interpolation,
c               but no CNO-interpolation")
c  ofe_brack = 0.0 = [O/Fe] value for READZEXCO
c
      parameter ( k_hz=1, ofe_brack=0.0d0 )
c
c PARAMETER badlogkval = 1.e+35 is stored to indicate missing Log(kappa) values
c
      parameter ( badlogkval=1.d+35 )
c
c PARAMETERS used during tests for high-T boundary of opacity storage matrices
c
      parameter ( ntm_m3=ntm-3, nt_m2=nt-2 )
c
c PARAMETERS used to specify logT6 and logR tabulation values
c
CANA
      parameter ( k81=nt-3, k80=k81-1, k60=k81-21, ks59=k60-1+ntdel )
!!      parameter ( flt81m6=8.1d0-6.d0, flt60m6=6.0d0-6.d0, 
!!    &     flt370m6=3.70d0-6.d0 )
      parameter ( flt81m6=9.3d0-6.d0, flt60m6=6.6d0-6.d0, 
     &     flt370m6=3.70d0-6.d0 )
      parameter ( flrmid = -3.0d0 )
c
c PARAMETERS defining the storage for the additional X-values from 'GN93hz':
c
      parameter ( mx_hi=2*mx, mo_m1=mo-1, mx_hi_nz=mx_hi*nz )
c
c COMMON /xhi_opal_z/ : auxiliary matrices for additional 'GN93hz' X-values:
c
      common /xhi_opal_z/ xhi_in(mx_hi), xcno_use(mx_hi,nz),
     $     xhi_use(mx_hi,nz), xxx_cno(mx_hi), xxx_hi(mx_hi),
     $     nx_hi(nz), ireq_hi(mx_hi), khighx(nz), kavail_xhi, kuse_xhi,
     $     kdo_xhi, kavail_cno, kuse_cno, kdo_cno, kavail_user,
     $     kuse_user, kdo_user
      save /xhi_opal_z/
c
c COMMON /a_opal_z/ : matrices for opacity storage: see OPALINIT
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
c COMMON /b_opal_z/ : high and low temperature limits, Z-values: see OPALINIT
c
      common/b_opal_z/ nta(0:nrm_p1),ntax0(0:nrm),
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
c COMMON /E_OPAL_Z/ : Return variables : see also instructions above
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
      common /extopac/ opvals(4),fedge,ftredge,fzedge
      save /extopac/
      equivalence (opvals(1),opact), (opvals(2),dopact),
     $     (opvals(3),dopacr), (opvals(4),dopactd)
c
c  Stored values of input parameters for this routine, in case it was called
c  from another interface
c
      common /x_opal_z/ z_opal, x_opal, xc_opal, xo_opal, slt_opal,
     $     slr_opal, fcn_opal, fcon_opal, fcnone_opal, fu_opal
      save /x_opal_z/
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
      z_opal = z
      x_opal = xh
      xc_opal = xci
      xo_opal = xoi
      slt_opal = slt
      slr_opal = slr
c
c..... set-up C/O axis points: making xci & xoi arguments of OPAL is simpler
c      than resetting xxc & xxo after changing them, and also avoids an error
c      if constants are used for these values in the calling program.
c
      xxc = xci
      xxo = xoi
      xxco = xxc + xxo
      if ( z+xh+xxco-1.d-6 .gt. 1.0d0 .or. z .lt. -1.d-8 .or.
     $     min(xh,xxc+z,xxo+z,xxco+z) .le. -1.d-6 ) then
         write(6,4397) z,xh,xxc,xxo
 4397    format(' '/' STOP -- OPAL: bad value(s) Z',f11.8,' X',f10.7,
     $        ' C',f10.7,' O',f10.7)
         stop ' STOP -- OPAL: bad composition. '
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
c..... If necessary, read data files and do initializations, using READZEXCO
c
      if ( itime .ne. 12345678 )
     $     call readzexco(nz,-1.,z,-1.,k_hz,-1,ofe_brack)
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
         f_xhi = 1.d0 - 100.d0 * max( xci + xoi - 0.2d0 , 0.0d0 )**2.d0
      endif
c
c  Check whether CNO- and/or user-interpolation will be needed
c
      if ( z .lt. 1.d-8 ) then
c
         need_cno = 0
         need_user = 0
c
      else
c
         if ( kdo_cno .gt. 0 .and. max( abs(fcn) ,
     $        abs(fcon) , abs(fcnone) ) .gt. 1.d-8 ) then
            need_cno = 1
         else
            need_cno = 0
         endif
c
         if ( kdo_user .gt. 0 .and. abs(fu) .gt. 1.d-8 ) then
            need_user = 1
         else
            need_user = 0
         endif
c
      endif
c
      fcn_opal = fcn * need_cno
      fcon_opal = fcon * need_cno
      fcnone_opal = fcnone * need_cno
      fu_opal = fu * need_user
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
            write(6,20) slt+6., slr
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
c      Try to avoid overestimating the extent of extrapolation into any cutout.
c
      if ( slt .gt. flt81m6 ) then
         k2sat = ntdel + min( nt , k81 + max( 0 ,
     $        int( ( slt - flt81m6 ) * dfs(nt) - 1.d-6 ) ) )
      else if ( slt .gt. flt60m6 ) then
         k2sat = ntdel + min( k80 , k60 + max( 0 ,
     $        int( ( slt - flt60m6 ) * dfs(k81) - 1.d-6 ) ) )
      else if ( k60 .le. 0 ) then
         k2sat = ntdel + k60
      else
         k2sat = min( ks59 , max( ntdel ,
     $        int( ( slt - flt370m6 ) * dfs(k60) + 1.d-6 ) ) )
      endif
c
      if ( slr .gt. flrmid ) then
         l2sat = nrdel + max( 0 , min( nr ,
     $        int( ( slr - alr(1) ) * dfsr(nr) + 0.999999d0 ) ) )
      else
         l2sat = nrdel + max( 0 , min( nr ,
     $        int( ( slr - alr(1) ) * dfsr(nr) + 1.000001d0 ) ) )
      endif
c
      k1x = -99
      l1x = -99
      k3sat = k2sat+1
      k2 = max(k2sat-ntdel,1)
      k3 = min(k3sat-ntdel,nt)
      l3sat = l2sat+1
      l4sat = min(l3sat,nre)+1
      l2 = max(l2sat-nrdel,1)
      l3 = min(l3sat-nrdel,nr)
      if ( min(k3,l3) .le. 0 .or. k2 .gt. nt .or. l2 .gt. nr ) then
         ftredge = 0.d0
         fedge = 0.d0
         if ( level_err .ge. 2 ) then
            write(6,20) slt+6., slr
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
            ft = max(min(alt(k2)-slt,0.d0)*dfs(k3)+1.000001d0,0.d0)
            fr = max(min(alr(l2)-slr,0.d0)*dfsr(l3)+1.000001d0,0.d0)
            if ( l2 .lt. nr ) ftredge = ftredge*fr
            if ( k2 .lt. nt ) ftredge = ftredge*ft
            k1x = k2-2
            l1x = l2-2
c						. . .	. . .   extrapolate (in
c					3.	. .*	. . .   either T or R)
c						. .	. .*    in a corner
         else if ( k2sat .lt. nta(l2sat) .and.
     $           k2sat .eq. nta(l3sat) ) then
            ft = max(min(alt(k2)-slt,0.d0)*dfs(k3)+1.000001d0,0.d0)
            fr = max(min(alr(l2)-slr,0.d0)*dfsr(l3)+1.000001d0,0.d0)
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
            if ( l2 .lt. nr ) ftredge = ftredge
     $           * max(min(alr(l2)-slr,0.d0)*dfsr(l3)+1.000001d0,0.d0)
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
                  ft = max(min(slt-alt(k3),0.d0)*dfs(k3)+1.000001d0,
     &                 0.d0)
                  fr = max(min(slr-alr(l3),0.d0)*dfsr(l3)+1.000001d0,
     &                 0.d0)
                  if ( l3 .gt. 1 ) ftredge = ftredge*fr
                  if ( k3 .gt. 1 ) ftredge = ftredge*ft
                  k1 = k3
                  l1 = l3
               endif
c						   ! 3. extrap T or R in corner
            else if ( k3sat .gt. ntaxat(l3sat) .and.
     $              k3sat .eq. ntaxat(l2sat) ) then
               ft = max(min(slt-alt(k3),0.d0)*dfs(k3)+1.000001d0,0.d0)
               fr = max(min(slr-alr(l3),0.d0)*dfsr(l3)+1.000001d0,0.d0)
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
                  if ( l3 .gt. 1 ) ftredge = ftredge
     $             * max(min(slr-alr(l3),0.d0)*dfsr(l3)+1.000001d0,0.d0)
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
               if ( k3 .gt. 1 ) ftredge = ftredge * max( 0.d0 ,
     $              min(slt-alt(k3),0.d0)*dfs(k3)+1.000001d0 )
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
                        ftredge = ftredge * max( 0.d0 ,
     $                       min(slt-alt(k3),0.d0)*dfs(k3)+1.000001d0 )
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
               xxc = xci
               xxo = xoi
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
               if ( xci .gt. 0.d0 ) then
                  xxc = cmod * xci
               else if ( xci .ge. -1.d-8 .or. z .lt. 1.d-8 ) then
                  xxc = 0.d0
               else
                  xxc = max( xci / z , -1. ) * zsto(kz)
               endif
               if ( xoi .gt. 0.d0 ) then
                  xxo = cmod * xoi
               else if ( xoi .ge. -1.d-8 .or. z .lt. 1.d-8 ) then
                  xxo = 0.d0
               else
                  xxo = max( xoi / z , -1.d0 ) * zsto(kz)
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
c-debug;     $           write(6,9409) z,xh,xci,xoi,10.**slt,slt,slr,
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
            xdelat = max( xdel*(delhi/dello)**2.d0, xdelmin )
            if ( delhi .lt. 0.1d0 )
     $           xdelat = xdelat + (xdel-xdelat)*((0.1d0-delhi)*
     &           12.5d0)**2.d0
            if ( xdelat .lt. xdel ) then
               is = 0
c			 ! get (mf,mg,mh) interpolated values with revised xdel
               do i = 1,4
                  opvals(i) = qzinter(is,1,xh,2,opk(mf,i),opk(mg,i),
     $                opk(mh,i),0.0d0,xa(mf),xa(mg),xa(mh),0.0d0,xdelat)
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
c  If the 'GN93hz' X-indices will be needed, obtain them:
c
      if ( f_xhi .gt. 0.d0 .or. max( need_cno , need_user ) .gt. 0 )then
c
c								! X < 0.1:
         if ( mf .eq. mxzero ) then
c					! set new mf,mg,mh,mf2 so that only the
c					! upper X-interp quadratic is shifted
            mf = 1
            mg = 1
            if ( mf2 .gt. 1 ) then
               mh = 2
               mf2 = 3
            else
               mh = 1
               mf2 = 1
            endif
c								! X = 1-Z:
         else if ( xh .gt. 0.999999d0 - z ) then
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
c
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
               xx_3 = xxx_hi(mh)
            else
               xx_3 = log10( x_3 + xdel )
            endif
            do i = 1, 4
               opvals(i) = opvals(i) + f_xhi * qchk(is,1,xxx,opk(1,i),
     $              opk(2,i),opk(3,i),xxx_hi(mf),xxx_hi(mg),xx_3)
               is = 1
            enddo
c					! 3 tables at high-X end of matrix
         else if ( mf .eq. mg ) then
            if ( x_3 .eq. xhi_in(mh) ) then
               xx_3 = xxx_hi(mh)
            else
               xx_3 = log10( x_3 + xdel )
            endif
            if ( x_4 .eq. xhi_in(mf2) ) then
               xx_4 = xxx_hi(mf2)
            else
               xx_4 = log10( x_4 + xdel )
            endif
            do i = 1, 4
               opvals(i) = opvals(i) + f_xhi * qchk(is,1,xxx,opk(1,i),
     $              opk(2,i),opk(3,i),xxx_hi(mg),xx_3,xx_4)
               is = 1
            enddo
c					! 3 tables: x4 = x3 (should not happen)
         else if ( x_3 .ge. x_4 ) then
            if ( x_4 .eq. xhi_in(mf2) ) then
               xx_4 = xxx_hi(mf2)
            else
               xx_4 = log10( x_4 + xdel )
            endif
            do i = 1, 4
               opvals(i) = opvals(i) + f_xhi * qchk(is,1,xxx,opk(1,i),
     $              opk(2,i),opk(4,i),xxx_hi(mf),xxx_hi(mg),xx_4)
               is = 1
            enddo
c		   ! 4 tables: interpolate X between two overlapping quadratics
         else
            if ( x_3 .eq. xhi_in(mh) ) then
               xx_3 = xxx_hi(mh)
            else
               xx_3 = log10( x_3 + xdel )
            endif
            if ( x_4 .eq. xhi_in(mf2) ) then
               xx_4 = xxx_hi(mf2)
            else
               xx_4 = log10( x_4 + xdel )
            endif
            dixr = ( xx_3 - xxx ) / ( xx_3 - xxx_hi(mg) )
            do i = 1,4
               opvals(i) = opvals(i) + f_xhi * ( qchk(is,1,xxx,opk(1,i),
     $              opk(2,i),opk(3,i),xxx_hi(mf),xxx_hi(mg),xx_3) * dixr
     $              + qchk(is,2,xxx,opk(2,i),opk(3,i),opk(4,i),
     $              xxx_hi(mg),xx_3,xx_4) * (1.d0-dixr) )
               is = 1
            enddo
         endif
c
      endif
c
c  If the CNO-interpolation is needed, perform it.
c
      if ( max( need_cno , need_user ) .gt. 0 ) then
c
         f_2 = fcn * need_cno
         f_3 = fcon * need_cno
         f_4 = fcnone * need_cno
         f_5 = fu * need_user
c				 ! m is the temporary-opacity-storage X-index
         m = 0
c			  ! loop ix over 'GN93hz'-opacity-shift X-indices
         do ix = mf, mf2
c								    ! ix valid?
            if ( ix .le. mg .or. ix .eq. mh .or. ix .eq. mf2 ) then
c
               m = m + 1
c					! get opacity shifts: loop over Z, T, R
               do kz = kzf, kzf2
c
                  ixm = min( ix , nx_hi(kz) )
c
                  call index_co_deltas( 2, ixm, ix2, ic2, io2 )
                  call index_co_deltas( 3, ixm, ix3, ic3, io3 )
                  call index_co_deltas( 4, ixm, ix4, ic4, io4 )
                  call index_co_deltas( 5, ixm, ix5, ic5, io5 )
c
                  do it = k1, k1 + ip
                     do ir = l1, l1 + iq(it-k1+1)
                        opl(it,ir,kz) = f_2 * co(ix2,ic2,io2,it,ir,kz)
     $                       + f_3 * co(ix3,ic3,io3,it,ir,kz)
     $                       + f_4 * co(ix4,ic4,io4,it,ir,kz)
     $                       + f_5 * co(ix5,ic5,io5,it,ir,kz)
                     enddo
                  enddo
c
               enddo
c						! interpolate over Z
               call qzlog4int( zlogd )
c						! interpolate over T and R
               call t6rinterp(slr,slt)
c
            endif
c
         enddo
c
c  Now add the just-computed 'GN93hz' CNO-interpolation opacity shifts and
c  their derivatives to the original opacity and derivative values:
c
         is = 0
c				    ! if 0.0 < X < 0.1 (1st quadratic absent):
c
         if ( mg .eq. 1 .and. mf2 .eq. 3 ) then
            do i = 1, 4
               opvals(i) = opvals(i) + quad(is,1,xxx,opk(1,i),opk(2,i),
     $              opk(3,i),xxx_cno(1),xxx_cno(2),xxx_cno(3))
               is = 1
            enddo
c					! if use only one X-table
         else if ( mf .eq. mh ) then
            do i = 1, 4
               opvals(i) = opvals(i) + opk(1,i)
            enddo
c			           ! 3 tables: interpolate in X using quadratic
         else if ( mh .eq. mf2 ) then
            if ( x_3 .eq. xhi_in(mh) ) then
               xx_3 = xxx_cno(mh)
            else
               xx_3 = log10( x_3 + xdel )
            endif
            do i = 1, 4
               opvals(i) = opvals(i) + qchk(is,1,xxx,opk(1,i),
     $              opk(2,i),opk(3,i),xxx_cno(mf),xxx_cno(mg),xx_3)
               is = 1
            enddo
c					! 3 tables at high-X end of matrix
         else if ( mf .eq. mg ) then
            if ( x_3 .eq. xhi_in(mh) ) then
               xx_3 = xxx_cno(mh)
            else
               xx_3 = log10( x_3 + xdel )
            endif
            if ( x_4 .eq. xhi_in(mf2) ) then
               xx_4 = xxx_cno(mf2)
            else
               xx_4 = log10( x_4 + xdel )
            endif
            do i = 1, 4
               opvals(i) = opvals(i) + qchk(is,1,xxx,opk(1,i),
     $              opk(2,i),opk(3,i),xxx_cno(mg),xx_3,xx_4)
               is = 1
            enddo
c					! 3 tables: x4 = x3 (should not happen)
         else if ( x_3 .ge. x_4 ) then
            if ( x_4 .eq. xhi_in(mf2) ) then
               xx_4 = xxx_cno(mf2)
            else
               xx_4 = log10( x_4 + xdel )
            endif
            do i = 1, 4
               opvals(i) = opvals(i) + qchk(is,1,xxx,opk(1,i),
     $              opk(2,i),opk(4,i),xxx_cno(mf),xxx_cno(mg),xx_4)
               is = 1
            enddo
c		   ! 4 tables: interpolate X between two overlapping quadratics
         else
            if ( x_3 .eq. xhi_in(mh) ) then
               xx_3 = xxx_cno(mh)
            else
               xx_3 = log10( x_3 + xdel )
            endif
            if ( x_4 .eq. xhi_in(mf2) ) then
               xx_4 = xxx_cno(mf2)
            else
               xx_4 = log10( x_4 + xdel )
            endif
            dixr = ( xx_3 - xxx ) / ( xx_3 - xxx_cno(mg) )
            do i = 1,4
               opvals(i) = opvals(i) + qchk(is,1,xxx,opk(1,i),opk(2,i),
     $              opk(3,i),xxx_cno(mf),xxx_cno(mg),xx_3) * dixr
     $              + qchk(is,2,xxx,opk(2,i),opk(3,i),opk(4,i),
     $              xxx_cno(mg),xx_3,xx_4) * ( 1.d0 - dixr )
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
c-debug;         write(6,8415) mf,mf2,kzf,kzf2,k1,ip,l1,iq(1),iq(2),
c-debug;     $        iq(3),iq(ip+1),z,xh,xci,xoi,slt,slr
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
c
c******************************************************************************
c
      subroutine ask_last_opac( op, dopt, dopr, doptd, fe, ftre, fze )
c     ================================================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      common/extopac/ opact,dopact,dopacr,dopactd,fedge,ftredge,fzedge
c
      op = opact
      dopt = dopact
      dopr = dopacr
      doptd = dopactd
      fe = fedge
      ftre = ftredge
      fze = fzedge
c
      return
      end
c
c******************************************************************************
c
      subroutine ask_last_xcnou(z,x,xc,xo,slt,slr,fcn,fcon,fcnone,fu)
c     ===============================================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      common /x_opal_z/ z_opal, x_opal, xc_opal, xo_opal, slt_opal,
     $     slr_opal, fcn_opal, fcon_opal, fcnone_opal, fu_opal
c
      z = z_opal
      x = x_opal
      xc = xc_opal
      xo = xo_opal
      slt = slt_opal
      slr = slr_opal
      fcn = fcn_opal
      fcon = fcon_opal
      fcnone = fcnone_opal
      fu = fu_opal
c
      return
      end
c
c******************************************************************************
c
      subroutine set_opal_dir( cdirin )
c     =================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      character*(*) cdirin
c
      common/recoin_opal_z/ itimeco,mxzero,mx03,kope,igznotgx
      save /recoin_opal_z/
c
      character*255 copdir
      common/opdir/ copdir
      save /opdir/
c
      character*1 cb(6)
      common/chkpoc/cb
      save /chkpoc/
c
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c---
c-dir;      logical lxst
c===
      last = lnblnk( cdirin )
      ibeg = max( 1 , non_blank_begin(cdirin) )
c
      kope = last - ibeg + 1
c
      if ( kope .eq. 0  ) then
c
         copdir = ' '
c
      else
c
         if ( level_err .gt. 0 ) then
            iblank = num_blanks_contained( cdirin )
            if ( iblank .gt. 0 ) then
               write(6,10) iblank, cdirin(ibeg:last)
 10            format(' WARNING:',i5,' blanks contained in OPAL',
     $              ' directory name:'/' ',a)
               if ( level_err .ge. 2 ) stop
     $              ' STOP -- SET_OPAL_DIR Error: blanks in name. '
            endif
         endif
c
         if ( cdirin(last:last) .ne. cb(1) .and.
     $        cdirin(last:last) .ne. cb(2) ) kope = kope + 1
c
         if ( kope .gt. 246 ) then
            write(6,20) kope, cdirin(ibeg:last)
 20         format(' Error: length',i5,
     $           ' exceeds 246 for OPAL directory name:'/' ',a)
            stop ' STOP -- SET_OPAL_DIR Error: name too long. '
         endif
c
         copdir = cdirin(ibeg:)
c
         if ( kope .gt. last - ibeg + 1 ) copdir(kope:kope) = cb(1)
c
c  NOTE that some systems return FALSE for the existence of a directory, so
c  one cannot check for the directory's existence.
c
c-dir;         call inqfil( copdir, lxst )
c-dir;         if ( .not. lxst ) then
c-dir;            write(6,30) copdir(:kope)
c-dir; 30         format(' STOP -- SET_OPAL_DIR:',
c-dir;     $           ' OPAL directory does not exist:'/' ',a)
c-dir;            stop
c-dir;         endif
c
      endif
c
      return
      end
c
c******************************************************************************
c
      subroutine set_ofe_file( cfileofe )
c     ===================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      character*(*) cfileofe
c
!!      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
!!      parameter ( nel_zmix=19, n_zmixes=7, kel_o=3, kel_fe=nel_zmix-1 )
      parameter ( nel_zmix=19, n_zmixes=8, kel_o=3, kel_fe=nel_zmix-1 )
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
      character*1 cb(6)
      common/chkpoc/cb
      save /chkpoc/
c
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c===
      last = lnblnk( cfileofe )
      ibeg = max( 1 , non_blank_begin(cfileofe) )
c
      if ( level_err .gt. 0 ) then
c
         if ( last - ibeg .ge. 8 ) then
            write(6,10) last - ibeg + 1, cfileofe(ibeg:last)
 10         format(' WARNING: length',i5,
     $           ' exceeds 8 for filename for khighz = 5:'/' ',a)
            if ( level_err .ge. 2 ) stop
     $           ' STOP -- SET_OFE_FILE Error: name too long. '
         endif
c
         iblank = num_blanks_contained( cfileofe )
c
         if ( iblank .gt. 0 ) then
            write(6,20) iblank, cfileofe(ibeg:last)
 20         format(' WARNING:',i5,' blanks contained in filename',
     $           ' for khighz = -5:'/' ',a)
            if ( level_err .ge. 2 ) stop
     $           ' STOP -- SET_OFE_FILE Error: blanks in name. '
         endif
c
      endif
c
      cfile_opalmixes(n_zmixes) = cfileofe(ibeg:)
c
      return
      end
c
c******************************************************************************
c
      subroutine set_altmix_ofe_file( cfileofe )
c     ==========================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      character*(*) cfileofe
c
!!      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
!!      parameter ( nel_zmix=19, n_zmixes=7, kel_o=3, kel_fe=nel_zmix-1 )
      parameter ( nel_zmix=19, n_zmixes=8, kel_o=3, kel_fe=nel_zmix-1 )
c
      parameter ( n_totmix = n_zmixes + 5, n_cnobeg = n_zmixes + 1 )
c
      character*255 cfile_opalGS98(n_totmix)
      common /opalGS98mixes/ bracketofe_opalGS98(n_totmix),
     $     xofe_opalGS98(n_totmix),xiz_opalGS98(nel_zmix,n_totmix),
     $     fninz_opalGS98(nel_zmix,n_totmix),atwt_opalGS98(nel_zmix),
     $     cfile_opalGS98
      save /opalGS98mixes/
c
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c===
      last = lnblnk( cfileofe )
c
      ibeg = max( 1 , non_blank_begin(cfileofe) )
c
      if ( level_err .gt. 0 ) then
c
         if ( last - ibeg .ge. 255 ) then
            write(6,10) last - ibeg + 1, cfileofe(ibeg:last)
 10         format(' WARNING: length',i5,
     $           ' exceeds 255 for filename for khighz = -5:'/' ',a)
            if ( level_err .ge. 2 ) stop
     $           ' STOP -- SET_ALTMIX_OFE_FILE Error: name too long. '
         endif
c
         iblank = num_blanks_contained( cfileofe )
c
         if ( iblank .gt. 0 ) then
            write(6,20) iblank, cfileofe(ibeg:last)
 20         format(' WARNING:',i5,' blanks contained in filename',
     $           ' for khighz = -5:'/' ',a)
            if ( level_err .ge. 2 ) stop
     $           ' STOP -- SET_ALTMIX_OFE_FILE Error: blanks in name. '
         endif
c
      endif
c
      cfile_opalGS98(n_zmixes) = cfileofe(ibeg:)
c
      return
      end
c
c******************************************************************************
c
      subroutine set_altmix_main_file( cfile_hz )
c     ===========================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      character*(*) cfile_hz
c
!!      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
!!      parameter ( nel_zmix=19, n_zmixes=7, kel_o=3, kel_fe=nel_zmix-1 )
      parameter ( nel_zmix=19, n_zmixes=8, kel_o=3, kel_fe=nel_zmix-1 )
c
      parameter ( n_totmix = n_zmixes + 5, n_cnobeg = n_zmixes + 1 )
c
      character*255 cfile_opalGS98(n_totmix)
      common /opalGS98mixes/ bracketofe_opalGS98(n_totmix),
     $     xofe_opalGS98(n_totmix),xiz_opalGS98(nel_zmix,n_totmix),
     $     fninz_opalGS98(nel_zmix,n_totmix),atwt_opalGS98(nel_zmix),
     $     cfile_opalGS98
      save /opalGS98mixes/
c
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c
      common /alt_change_opal_z/ main_alt_change, iulow, khighz_in,
     $     ofebrack_in
      save /alt_change_opal_z/
c===
      last = lnblnk( cfile_hz )
c
      ibeg = max( 1 , non_blank_begin(cfile_hz) )
c
      if ( level_err .gt. 0 ) then
c
         if ( last - ibeg .ge. 255 ) then
            write(6,10) last - ibeg + 1, cfile_hz(ibeg:last)
 10         format(' WARNING: length',i5,
     $           ' exceeds 255 for filename for khighz = -1:'/' ',a)
            if ( level_err .ge. 2 ) stop
     $           ' STOP -- SET_ALTMIX_MAIN_FILE Error: name too long. '
         endif
c
         iblank = num_blanks_contained( cfile_hz )
c
         if ( iblank .gt. 0 ) then
            write(6,20) iblank, cfile_hz(ibeg:last)
 20         format(' WARNING:',i5,' blanks contained in filename',
     $           ' for khighz = -5:'/' ',a)
            if ( level_err .ge. 2 ) stop
     $           ' STOP -- SET_ALTMIX_MAIN_FILE Error: blanks in name '
         endif
c
      endif
c
      if ( last .eq. 0 ) then
         cfile_opalGS98(1) = 'GS98hz'
         main_alt_change = 2
      else
         cfile_opalGS98(1) = cfile_hz(ibeg:)
         main_alt_change = 1
      endif
c
      return
      end
c
c******************************************************************************
c
      subroutine set_cno_files( cf_hz, cf_c, cf_o, cf_n, cf_user )
c     ============================================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      character*(*) cf_hz, cf_c, cf_o, cf_n, cf_user
c
!!      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
!!      parameter ( nel_zmix=19, n_zmixes=7, kel_o=3, kel_fe=nel_zmix-1 )
      parameter ( nel_zmix=19, n_zmixes=8, kel_o=3, kel_fe=nel_zmix-1 )
c
      parameter ( n_totmix = n_zmixes + 5, n_cnobeg = n_zmixes + 1 )
c
      character*255 cfile_opalGS98(n_totmix)
      common /opalGS98mixes/ bracketofe_opalGS98(n_totmix),
     $     xofe_opalGS98(n_totmix),xiz_opalGS98(nel_zmix,n_totmix),
     $     fninz_opalGS98(nel_zmix,n_totmix),atwt_opalGS98(nel_zmix),
     $     cfile_opalGS98
      save /opalGS98mixes/
c
      character*8 cdef_CNO_ext(n_cnobeg:n_totmix)
      common /ext_CNO_opal_z/ len_def_CNO_ext(n_cnobeg:n_totmix),
     $     cdef_CNO_ext
      save /ext_CNO_opal_z/
c
      parameter ( nbegp1 = n_cnobeg + 1, nbegp2 = n_cnobeg + 2 )
c
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c===
      l_max = 0
      iwarn = max( 2 - level_err , 0 )
c
      do k = n_cnobeg, n_totmix
c
         if ( k .eq. n_cnobeg ) then
            last = lnblnk( cf_hz )
            ibeg = max( 1 , non_blank_begin(cf_hz) )
            cfile_opalGS98(k) = cf_hz(ibeg:)
         else if ( k .eq. nbegp1 ) then
            last = lnblnk( cf_c )
            ibeg = max( 1 , non_blank_begin(cf_c) )
            cfile_opalGS98(k) = cf_c(ibeg:)
         else if ( k .eq. nbegp2 ) then
            last = lnblnk( cf_o )
            ibeg = max( 1 , non_blank_begin(cf_o) )
            cfile_opalGS98(k) = cf_o(ibeg:)
         else if ( k .eq. n_totmix ) then
            last = lnblnk( cf_user )
            ibeg = max( 1 , non_blank_begin(cf_user) )
            cfile_opalGS98(k) = cf_user(ibeg:)
         else
            last = lnblnk( cf_n )
            ibeg = max( 1 , non_blank_begin(cf_n) )
            cfile_opalGS98(k) = cf_n(ibeg:)
         endif
c
         if ( k .eq. n_cnobeg ) then
            l_hz = last - ibeg + 1
         else if ( last .eq. 0 .and. l_hz .gt. 0 ) then
            cfile_opalGS98(k) =
     $           cfile_opalGS98(n_cnobeg)(:min(l_hz,255)) //
     $           cdef_CNO_ext(k)
            last = l_hz + len_def_CNO_ext(k)
         else if ( last .gt. 0 .and. l_hz .eq. 0 .and.
     $           iwarn .eq. 0 ) then
            write(6,10)
 10         format(' WARNING: SET_CNO_FILES: blank  cfile_hz',
     $           '  but other input filename(s) non-blank')
            iwarn = iwarn + 1
         endif
c
         l_max = max( l_max , last - ibeg + 1 )
c
      enddo
c
      if ( l_max .gt. 255 .and. level_err .gt. 0 ) then
         write(6,20) l_max
 20      format(' WARNING: SET_CNO_FILES: largest filename length',i5,
     $        ' exceeds 255')
         if ( level_err .ge. 2 ) stop
     $        ' STOP -- SET_CNO_FILES Error: filename too long. '
      endif
c
      return
      end
c
c******************************************************************************
c
      subroutine ask_khighz_ofe( khighz_used, ofebrack_used )
c     =======================================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      common /alt_change_opal_z/ main_alt_change, iulow, khighz_in,
     $     ofebrack_in
      save /alt_change_opal_z/
c
      khighz_used = khighz_in
      ofebrack_used = ofebrack_in
c
      return
      end
c
c******************************************************************************
c
      subroutine ask_main_opal_file( cf_main_used )
c     =============================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      character*(*) cf_main_used
c
      common /alt_change_opal_z/ main_alt_change, iulow, khighz_in,
     $     ofebrack_in
      save /alt_change_opal_z/
c
!!      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
!!      parameter ( nel_zmix=19, n_zmixes=7, kel_o=3, kel_fe=nel_zmix-1 )
      parameter ( nel_zmix=19, n_zmixes=8, kel_o=3, kel_fe=nel_zmix-1 )
c
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
      parameter ( n_totmix = n_zmixes + 5, n_cnobeg = n_zmixes + 1 )
c
c-implicit;      real*4 bracketofe_opalGS98, xofe_opalGS98, xiz_opalGS98,
c-implicit;     $     fninz_opalGS98, atwt_opalGS98
      character*255 cfile_opalGS98(n_totmix)
      common /opalGS98mixes/ bracketofe_opalGS98(n_totmix),
     $     xofe_opalGS98(n_totmix),xiz_opalGS98(nel_zmix,n_totmix),
     $     fninz_opalGS98(nel_zmix,n_totmix),atwt_opalGS98(nel_zmix),
     $     cfile_opalGS98
      save /opalGS98mixes/
c
      if ( khighz_in .ge. 0 ) then
         cf_main_used = cfile_opalmixes(1)
      else
         cf_main_used = cfile_opalGS98(1)
      endif
c
      return
      end
c
c******************************************************************************
c
      subroutine set_xhi( kxhi )
c     ==========================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      parameter ( mx_hi=2*mx, mo_m1=mo-1, mx_hi_nz=mx_hi*nz )
c
      common /xhi_opal_z/ xhi_in(mx_hi), xcno_use(mx_hi,nz),
     $     xhi_use(mx_hi,nz), xxx_cno(mx_hi), xxx_hi(mx_hi),
     $     nx_hi(nz), ireq_hi(mx_hi), khighx(nz), kavail_xhi, kuse_xhi,
     $     kdo_xhi, kavail_cno, kuse_cno, kdo_cno, kavail_user,
     $     kuse_user, kdo_user
      save /xhi_opal_z/
c						! set high-X flag value
      kuse_xhi = max( 0 , min( 2 , kxhi ) )
      if ( kavail_xhi .le. 0 ) then
         kdo_xhi = 0
      else
         kdo_xhi = kuse_xhi
      endif
c
      return
      end
c
c******************************************************************************
c
      subroutine ask_xhi( kxhi, kavail )
c     ==================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      parameter ( mx_hi=2*mx, mo_m1=mo-1, mx_hi_nz=mx_hi*nz )
c
      common /xhi_opal_z/ xhi_in(mx_hi), xcno_use(mx_hi,nz),
     $     xhi_use(mx_hi,nz), xxx_cno(mx_hi), xxx_hi(mx_hi),
     $     nx_hi(nz), ireq_hi(mx_hi), khighx(nz), kavail_xhi, kuse_xhi,
     $     kdo_xhi, kavail_cno, kuse_cno, kdo_cno, kavail_user,
     $     kuse_user, kdo_user
      save /xhi_opal_z/
c						! return high-X flag value
      kxhi = kuse_xhi
c						! return availability flag
      kavail = max( 0 , min( 1 , kavail_xhi ) )
c
      return
      end
c
c******************************************************************************
c
      subroutine set_cno_interp( kcno, kuser )
c     ========================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      parameter ( mx_hi=2*mx, mo_m1=mo-1, mx_hi_nz=mx_hi*nz )
c
      common /xhi_opal_z/ xhi_in(mx_hi), xcno_use(mx_hi,nz),
     $     xhi_use(mx_hi,nz), xxx_cno(mx_hi), xxx_hi(mx_hi),
     $     nx_hi(nz), ireq_hi(mx_hi), khighx(nz), kavail_xhi, kuse_xhi,
     $     kdo_xhi, kavail_cno, kuse_cno, kdo_cno, kavail_user,
     $     kuse_user, kdo_user
      save /xhi_opal_z/
c						! set CNO-interpolation flags
      kuse_cno = max( min( kcno , 1 ) , 0 )
      kuse_user = max( min( kuser , 1 ) , 0 )
      kdo_cno = kuse_cno * kavail_cno
      kdo_user = kuse_user * kavail_user
c
      return
      end
c
c******************************************************************************
c
      subroutine ask_cno_interp( kcno, kuser, kcno_avail, kuser_avail )
c     =================================================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      parameter ( mx_hi=2*mx, mo_m1=mo-1, mx_hi_nz=mx_hi*nz )
c
      common /xhi_opal_z/ xhi_in(mx_hi), xcno_use(mx_hi,nz),
     $     xhi_use(mx_hi,nz), xxx_cno(mx_hi), xxx_hi(mx_hi),
     $     nx_hi(nz), ireq_hi(mx_hi), khighx(nz), kavail_xhi, kuse_xhi,
     $     kdo_xhi, kavail_cno, kuse_cno, kdo_cno, kavail_user,
     $     kuse_user, kdo_user
      save /xhi_opal_z/
c					! return CNO-interpolation flags
      kcno = kuse_cno
      kuser = kuse_user
      kcno_avail = kavail_cno
      kuser_avail = kavail_user
c
      return
      end
c
c******************************************************************************
c
      subroutine set_err_check( level )
c     =================================
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

c						! set error-checking level
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c
      level_err = max( 0 , min( 3 , level ) )
c
      return
      end
c
c******************************************************************************
c
      subroutine ask_err_check( level )
c     =================================
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

c						! return error-checking level
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c
      level = level_err
c
      return
      end
c
c******************************************************************************
c
      subroutine set_logt6_limits( vlo, dvlo, vhi, dvhi )
c     ===================================================
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

c
!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
      common/b_opal_z/ nta(0:nrm_p1),ntax0(0:nrm),
     $     ntax03(0:nrm), sltlo, slthi, dltlo_inv, dlthi_inv,
     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
      save /b_opal_z/
c===						! set T,R-values, if necessary
      if ( init_trvals .le. 0 ) call get_trvals
c
c  Set the logT6 limits, according to the input values; by default, make
c  1-grid-point extrapolation beyond the matrix edge just within the limits:
c
      if ( vlo .gt. -90.d0 ) sltlo = max( alt(1) , min( alt(nt) , vlo ))
c
      if ( dvlo .lt. -90.d0 ) then
         dltlo_inv = max( dltlo_inv ,
     $        1.d0 / ( sltlo - alt(1) + 0.999999d0 / dfs(1) ) )
      else if ( dvlo .lt. 0.d0 ) then
         dltlo_inv = dfs(1) * 0.999999d0
      else
         dltlo_inv = 1.d0 / min( max( dvlo , 1.d-6 ) ,
     $        sltlo - alt(1) + 0.999999d0 / dfs(1) )
      endif
c
      if ( vhi .gt. -90.d0 ) slthi = max( alt(1) , min( alt(nt) , vhi ))
c
      if ( dvhi .lt. -90.d0) then
         dlthi_inv = max( dlthi_inv ,
     $        1.d0 / ( alt(nt) - slthi + 0.999999d0 / dfs(nt) ) )
      else if ( dvhi .lt. 0.d0 ) then
         dlthi_inv = dfs(nt) * 0.999999d0
      else
         dlthi_inv = 1.d0 / min( max( dvhi , 1.d-6 ) ,
     $        alt(nt) - slthi + 0.999999d0 / dfs(nt) )
      endif
c
      return
      end
c
c******************************************************************************
c
      subroutine set_logr_limits( vlo, dvlo, vhi, dvhi )
c     ==================================================
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

c
!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
      common/b_opal_z/ nta(0:nrm_p1),ntax0(0:nrm),
     $     ntax03(0:nrm), sltlo, slthi, dltlo_inv, dlthi_inv,
     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
      save /b_opal_z/
c===						! set T,R-values, if necessary
      if ( init_trvals .le. 0 ) call get_trvals
c
c  Set the logR limits, according to the input values; by default, make
c  1-grid-point extrapolation beyond the matrix edge just within the limits:
c
      if ( vlo .gt. -90.d0 ) slrlo = max( alr(1) , min( alr(nr) , vlo ))
c
      if ( dvlo .lt. -90.d0 ) then
         dlrlo_inv = max( dlrlo_inv ,
     $        1.d0 / ( slrlo - alr(1) + 0.999999d0 / dfsr(1) ) )
      else if ( dvlo .lt. 0.d0 ) then
         dlrlo_inv = dfsr(1) * 0.999999d0
      else
         dlrlo_inv = 1.d0 / min( max( dvlo , 1.d-6 ) ,
     $        slrlo - alr(1) + 0.999999d0 / dfsr(1) )
      endif
c
      if ( vhi .gt. -90.d0 ) slrhi = max( alr(1) , min( alr(nr) , vhi ))
c
      if ( dvhi .lt. -90.d0 ) then
         dlrhi_inv = max( dlrhi_inv ,
     $        1.d0 / ( alr(nr) - slrhi + 0.999999d0 / dfsr(nr) ) )
      else if ( dvhi .lt. 0.d0 ) then
         dlrhi_inv = dfsr(nr) * 0.999999d0
      else
         dlrhi_inv = 1.d0 / min( max( dvhi , 1.d-6 ) ,
     $        alr(nr) - slrhi + 0.999999d0 / dfsr(nr) )
      endif
c
      return
      end
c
c******************************************************************************
c
      subroutine reset_z_limits( vlo, dvlo, vhi, dvhi )
c     =================================================
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

c
!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c				! if no opacities were read in, cannot reset!
c
      if ( itime .ne. 12345678 ) return
c
      if ( min( vlo , vhi ) .gt. -1.d-6 ) then
         zlow = max( zsto(1) , min( vlo , vhi , zsto(numz) ) )
         zhigh = min( zsto(numz) , max( vhi , vlo , zsto(1) ) )
      else if ( vlo .gt. -1.d-6 ) then
         zlow = max( zsto(1) , min( vlo , zhigh , zsto(numz) ) )
      else if ( vhi .gt. -1.d-6 ) then
         zhigh = min( zsto(numz) , max( vhi , zlow , zsto(1) ) )
      endif
c
      if ( dvlo .gt. -1.d-6 ) then
         if ( zlow .le. zsto(1) + zacc(1) ) then
            zlo_ex = zlow - max( dvlo , zacc(1) )
         else
            zlo_ex = zlow - max( dvlo , zacc(numz) )
         endif
      else if ( zlow .le. zsto(1) + zacc(1) ) then
         zlo_ex = min( zlo_ex , zlow - zacc(1) )
      else
         zlo_ex = min( zlo_ex , zlow - zacc(numz) )
      endif
c
      if ( dvhi .gt. -1.d-6 ) then
         zhi_ex = zhigh + max( dvhi , zacc(numz) )
      else
         zhi_ex = max( zhi_ex , zhigh + zacc(numz) )
      endif
c
      return
      end
c
c******************************************************************************
c
      subroutine ask_logt6_limits( vlo, dvlo, vhi, dvhi )
c     ===================================================
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

c
!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/b_opal_z/ nta(0:nrm_p1),ntax0(0:nrm),
     $     ntax03(0:nrm), sltlo, slthi, dltlo_inv, dlthi_inv,
     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
      save /b_opal_z/
c===						! set T,R-values, if necessary
      if ( init_trvals .le. 0 ) call get_trvals
c
      vlo = sltlo
      dvlo = 1.d0 / dltlo_inv
      vhi = slthi
      dvhi = 1.d0 / dlthi_inv
c
      return
      end
c
c******************************************************************************
c
      subroutine ask_logr_limits( vlo, dvlo, vhi, dvhi )
c     ==================================================
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

c
!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/b_opal_z/ nta(0:nrm_p1),ntax0(0:nrm),
     $     ntax03(0:nrm), sltlo, slthi, dltlo_inv, dlthi_inv,
     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
      save /b_opal_z/
c===						! set T,R-values, if necessary
      if ( init_trvals .le. 0 ) call get_trvals
c
      vlo = slrlo
      dvlo = 1.d0 / dlrlo_inv
      vhi = slrhi
      dvhi = 1.d0 / dlrhi_inv
c
      return
      end
c
c******************************************************************************
c
      subroutine ask_z_limits( nzmax, zmin, zmax )
c     ============================================
c
c  Returns NZ (maximum allowed stored Z-values) and Z-limits (0.0 and 0.1)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c===
      nzmax = nz
      zmin = 0.0d0
      zmax = 0.1d0
c
      return
      end
c
c******************************************************************************
c
      subroutine ask_z_use( nzuse, zlo, zmid, zhi, zloex, zhiex )
c     ===========================================================
c
c  Returns the current values of numz, zlow, zmiddle, zhigh, zlo_ex, zhi_ex
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!    $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c===
      if ( itime .ne. 12345678 ) then
         nzuse = 0
         zlo = 0.d0
         zmid = 0.02d0
         zhi = 0.1d0
         zloex = -1.d-8
         zhiex = 0.12d0
      else
         nzuse = numz
         zlo = zlow
         zmid = zmiddle
         zhi = zhigh
         zloex = zlo_ex
         zhiex = zhi_ex
      endif
c
      return
      end
c
c******************************************************************************
c
      subroutine ask_z_array( kzstart, karraystart, zarray, narray )
c     ==============================================================
c
c  Returns Z-values from zsto(), starting with element kzstart, in the
c  array zarray(narray), starting with element karraystart; any excess
c  elements in zarray() after nums is reached are filled with values of -1.
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      dimension zarray(narray)
c
!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c===
      if ( itime .ne. 12345678 ) then
         nzuse = 0
      else
         nzuse = numz
      endif
c
      k_z = max( kzstart , 1 )
      k_a = karraystart
      do while ( k_z .le. nzuse .and. k_a .le. narray )
         zarray(k_a) = zsto(k_z)
         k_a = k_a + 1
         k_z = k_z + 1
      enddo
      do while ( k_a .le. narray )
         zarray(k_a) = -1.d0
         k_a = k_a + 1
      enddo
c
      return
      end
c
c******************************************************************************
c
      subroutine set_smooth( initsmooth, lowCOsmooth, interpCOsmooth )
c     ================================================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      common/c_opal_ctrl_smooth/ init_smo, low_CO_smo, interp_CO_smo
      save /c_opal_ctrl_smooth/
c
      if ( initsmooth .ge. 0 ) init_smo = min( initsmooth , 2 )
c
      if ( lowCOsmooth .ge. 0 ) low_CO_smo = min( lowCOsmooth , 1 )
c
      if ( interpCOsmooth .ge. 0 )
     $     interp_CO_smo = min( interpCOsmooth , 1 )
c
      return
      end
c
c******************************************************************************
c
      subroutine ask_smooth( initsmooth, lowCOsmooth, interpCOsmooth )
c     ================================================================
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

c
      common/c_opal_ctrl_smooth/ init_smo, low_CO_smo, interp_CO_smo
      save /c_opal_ctrl_smooth/
c
      initsmooth = init_smo
c
      lowCOsmooth = low_CO_smo
c
      interpCOsmooth = interp_CO_smo
c
      return
      end
c
c******************************************************************************
c
      subroutine readco(z,kallrd,khighz,iu_lo)
c     ========================================
c
c..... The purpose of this subroutine is to read the data tables; actually,
c      it just calls READEXCO to do the work, setting [O/Fe] = 0.0
c
c Z is the metallicity; opacities will be interpolated (quadratically) in Z if
c   necessary, with values of Z from 0.0 to 0.1 being allowed.
c kallrd is ignored (it is present only for backward compatibility).
c if khighz = 0 , then the file GN93hz is not used; else, GN93hz (or GS98hz,
c   for khighz < 0) may be used for the case C=O=0.0 to get improved opacities.
c iu_lo is the unit number from which the lowest Z-value files are read;
c   units iu_lo thru iu_lo+3 may be needed.
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

c===
      call readexco(z,kallrd,max(min(khighz,1),-1),iu_lo,0.0d0)
c
      return
      end
c
c******************************************************************************
c
      subroutine readexco(z,kallrd,khighz,iu_lo,ofebrack)
c     ===================================================
c
c..... The purpose of this subroutine is to read the data tables
c
c Z is the metallicity; opacities will be interpolated (quadratically) in Z if
c   necessary, with values of Z from 0.0 to 0.1 being allowed.
c kallrd is ignored (it is present only for backward compatibility).
c if khighz = 0 , then the file GN93hz is not used;
c           = 1 or more, then GN93hz is used if it will improve Z-interpolation
c           = 2 , then Alrd96a2 and GN93hz are used to interpolate in [O/Fe]
c           = 3 , then C95hz and GN93hz are used to interpolate in [O/Fe]
c           = 4 , then W95hz and GN93hz are used to interpolate in [O/Fe]
c           = 5 , then GN93hz and a user-specified file are used to interpolate
c                  in [O/Fe]
c           = -1 thru -5 : use GS98hz (and its variants with [O/Fe] > 0.0, for
c                  khighz < -1) to get alternate-mix opacities (note that the
c                  file GN93hz is still read in, for reference purposes)
c           = -11 thru -15, 11 thru 15: same as khighz = -1 thru -5 and 1 thru
c                  5, except that CNO-interpolation opacity files are read in
c           = -21 thru -25, 21 thru 25: same as khighz = -1 thru -5 and 1 thru
c                  5, except that a user-specified opacity interpolation file
c                  is read in
c           = -31 thru -35, 31 thru 35: same as khighz = -1 thru -5 and 1 thru
c                  5, except that BOTH the CNO- and user-specified opacity
c                  interpolation files are read in
c iu_lo is the unit number from which the lowest Z-value files are read;
c   units iu_lo thru iu_lo+3 may be needed.
c ofebrack is [O/Fe], logarithmic oxygen (or alpha-element) enhancement factor,
c   relative to the Sun: ofebrack = log10{ (O_mix/Fe_mix) / (O_sun/Fe_sun) } .
c   If khighz is 0 or 1, then ofebrack is ignored; otherwise, interpolate
c   logKappa linearly between mix 1 and mix khighz, the interpolation factors
c   being such as to yield the desired [O/Fe] from combining these mixes.
c   Note that GN93hz has [O/Fe] = 0.0, Alrd96a2 has [O/Fe] = 0.3, C95hz has
c   [O/Fe] = 0.4, and W95hz has [O/Fe] = 0.5, but they have different patterns
c   of enhancements for elements other than oxygen.  For khighz < 0, it is the
c   new file GS98hz that is defined to have [O/Fe] = 0.0, and its variants (by
c   default, GS98hz_OFe.3_Alrd96a2, GS98hz_OFe.4_C95, and GS98hz_OFe.5_W95) are
c   the ones used to interpolate in [O/Fe].
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
      parameter ( mx_hi=2*mx )
c
      common /xhi_opal_z/ xhi_in(mx_hi), xcno_use(mx_hi,nz),
     $     xhi_use(mx_hi,nz), xxx_cno(mx_hi), xxx_hi(mx_hi),
     $     nx_hi(nz), ireq_hi(mx_hi), khighx(nz), kavail_xhi, kuse_xhi,
     $     kdo_xhi, kavail_cno, kuse_cno, kdo_cno, kavail_user,
     $     kuse_user, kdo_user
      save /xhi_opal_z/
c===
      numz = 1
      if ( z .lt. -1.d-6 ) then
         zat = 0.02d0
      else if ( z .lt. 1.d-8 ) then
         zat = 0.d0
      else
         zat = z
      endif
      kavail_cno = 1
      kavail_user = 1
c
      call read_kz(1,zat,kallrd,khighz,iu_lo,ofebrack)
c
      zlow = zat - zacc(1)
      zmiddle = zat
      zhigh = zat + zacc(1)
      zlo_ex = zlow - zacc(1)
      zhi_ex = zhigh + zacc(1)
      dfsz(1) = 1.d0
c
      call finish_cno
c
      return
      end
c
c******************************************************************************
c
      subroutine readzexco(nzin,zlo,z,zhi,khighz,iu_lo,ofebrack)
c     ==========================================================
c
c  Similar to READEXCO, except that a range of  nzin  Z-values from  zlo  to
c  zhi  are used: see comments near the beginning of this file.
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
      parameter ( mx_hi=2*mx )
c
      common /xhi_opal_z/ xhi_in(mx_hi), xcno_use(mx_hi,nz),
     $     xhi_use(mx_hi,nz), xxx_cno(mx_hi), xxx_hi(mx_hi),
     $     nx_hi(nz), ireq_hi(mx_hi), khighx(nz), kavail_xhi, kuse_xhi,
     $     kdo_xhi, kavail_cno, kuse_cno, kdo_cno, kavail_user,
     $     kuse_user, kdo_user
      save /xhi_opal_z/
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
      common /alt_change_opal_z/ main_alt_change, iulow, khighz_in,
     $     ofebrack_in
      save /alt_change_opal_z/
c
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c___
      dimension z_use(nzm)
c===
      kavail_cno = 1
      kavail_user = 1
c
      if ( level_err .gt. 0 .and. main_alt_change .ne. 0 .and.
     $     khighz .lt. -1 .and. mod( abs(khighz), 10 ) .gt. 1 .and.
     $     mod( abs(khighz), 10 ) .lt. 5 ) then
         write(6,5) khighz
 5       format(' WARNING: khighz=',i4,' used after',
     $        ' ''GS98hz'' file replaced: may yield bad [O/Fe].')
         if ( level_err .gt. 1 ) stop
     $        ' STOP -- READZEXCO: Error: non-GS98hz, khighz=-2,-3,-4 '
      endif
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
c  If there is only one Z-value to store, call readexco instead, and return
c  (note that readexco will call read_kz and finish_cno):
c
      if ( numz .eq. 1 ) then
c
         call readexco(zat,1,khighz,iu_lo,ofebrack)
c
         if ( nzin .eq. 1 ) then
            if ( zlo .ge. -1.d-6 ) then
               zlow = min( zlow , zlo )
               zlo_ex = min(zlo_ex , max( 2.d0 * zlow - zmiddle , 0.d0))
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
         write(6,20) nzin, zlo, zat, zhi, 'Z > 0.1 is not allowed!!!'
 20      format(' '/' STOP -- READZEXCO(Nzin=',i2,
     $        '): bad Z values',3f12.8/'   ',a,a)
         stop ' STOP -- READZEXCO: Error: Z > 0.1 is not allowed!!! '
      endif
c					! check for Zlo, Z, Zhi out of order
      zlo_t = max( zlo , 0.0d0 )
      zhi_t = max( zhi , 0.0d0 )
      if ( zlo .lt. -1.e-6 ) zlo_t = zat
      if ( zhi .lt. -1.e-6 ) zhi_t = zat
      if ( level_err .gt. 0 .and. ( zat .gt. zhi_t + 1.d-6 .or.
     $     zlo_t .gt. min( zat , zhi_t ) + 1.e-6 ) ) then
         write(6,20) nzin, zlo_t, zat, zhi_t,
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
               if ( zhi_t - zlo_t .lt. 1.d-5 ) then
                  write(6,20) nzin, zlo_t, z, zhi_t,
     $                 'Zhi - Zlo < 1.E-5: Zlo and Zhi',
     $                 ' too close together for Nzin = 2'
                  stop ' STOP -- READZEXCO: Error: bad Z values. '
               else if ( zlo_t * 1.01d0 .gt. zhi_t + 1.d-6 ) then
                  write(6,20) nzin, zlo_t, z, zhi_t,
     $                 'Zhi / Zlo < 1.01: Zlo and Zhi',
     $                 ' too close together for Nzin = 2'
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
               write(6,20) nzin, zlo_t, z, zhi_t,
     $              'Zlo / Zhi < 0.6 (and Zhi-Zlo>0.0002):',
     $              'Zlo and Zhi too far apart for Nzin = 2'
               stop ' STOP -- READZEXCO: Error: bad Z values. '
            endif
         endif
c
c  (will be linear interpolation in log[Z+0.001]); default range is plus/minus
c  10 percent in Z or 2.e-5, minimum range is at least 1 percent or 1.e-5
c
         if ( nzin .gt. 2 .or. zlo .lt. -1.e-6 .or. zlo .gt. zat ) then
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
         else if ( zhi - zlo .lt. 1.d-5 .or. zlo * 1.01d0 .gt. zhi) then
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
         call read_kz(1,zat,1,khighz,iu_lo,ofebrack)
c
         zat = zhigh
         call read_kz(2,zat,1,khighz,iu_lo,ofebrack)
c
         zlo_ex = zlow - 0.5d0 * ( zhigh - zlow )
         zhi_ex = zhigh + 0.5d0 * ( zhigh - zlow )
         dfsz(2) = 1.d0 / ( zvint(2) - zvint(1) )
         dfsz(1) = dfsz(2)
c
         call finish_cno
c
         return
c
c  Else if all available Z-values should be used, handle this case and return:
c
      else if ( numz .eq. nzm ) then
c
         zlow = 0.0d0
         zmiddle = zat
         zhigh = zavail(nzm)
         zlo_ex = -1.d-6
         zhi_ex = 2.d0 * zavail(nzm) - zavail(nzm-1)
c
         do kz = 1, nzm
            zat = zavail(kz)
            call read_kz(kz,zat,1,khighz,iu_lo,ofebrack)
            if ( kz .gt. 1 )
     $           dfsz(kz) = 1.d0 / ( zvint(kz) - zvint(kz-1) )
         enddo
c
         dfsz(1) = dfsz(2)
c
         call finish_cno
c
         return
c
      endif
c
c  If there should be at least three Z-values:
c						     ! check if Z-range too big
      if ( level_err .gt. 0 ) then
         j = mz
         do while ( j .gt. 2 .and. zhi_t .lt. zval(j-1) + 1.d-6 )
            j = j - 1
         enddo
         i = 1
         do while ( i .lt. j - 1 .and. zlo_t .gt. zval(i+1) - 1.d-6 )
            i = i + 1
         enddo
         if ( nzin .lt. j - i ) then
            write(6,20) nzin, zlo_t, z, zhi_t,
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
      if ( zlow .lt. 1.e-8 ) then
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
            call read_kz(kz,zat,1,khighz,iu_lo,ofebrack)
            if ( kz .gt. 1 )
     $           dfsz(kz) = 1.d0 / ( zvint(kz) - zvint(kz-1) )
         enddo
c
         dfsz(1) = dfsz(2)
c
         call finish_cno
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
         call read_kz(1,zat,1,khighz,iu_lo,ofebrack)
         zat = zmiddle
         call read_kz(2,zat,1,khighz,iu_lo,ofebrack)
         zat = zhigh
         call read_kz(3,zat,1,khighz,iu_lo,ofebrack)
c
         dfsz(3) = 1.d0 / ( zvint(3) - zvint(2) )
         dfsz(2) = 1.d0 / ( zvint(2) - zvint(1) )
         dfsz(1) = dfsz(2)
c
         zlo_ex = min( zlow - zacc(1) ,
     $        10.d0**( zl_1 - min( zl_2 - zl_1 ,
     $        log10( ( zavail(ilo1+1) + zdel )
     $        / ( zavail(ilo1) + zdel ) ) ) ) - zdel )
         if ( zlo_ex .lt. 1.d-8 .and. zlo_ex .gt. 0.d0 ) zlo_ex = 0.d0
         zhi_ex = min( 10.d0**( 2.d0 * zl_3 - zl_2 ) - zdel ,
     $        zhigh + zavail(ihi2) - zavail(ihi2-1) )
c
         call finish_cno
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
     $        zlow .gt. zval(jlo1+1) - 1.e-6 )
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
         call read_kz(kz,z_use(kz),1,khighz,iu_lo,ofebrack)
         if ( kz .gt. 1 ) dfsz(kz) = 1. / ( zvint(kz) - zvint(kz-1) )
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
c
      call finish_cno
c							! we are done: return
      return
      end
c
c******************************************************************************
c
      subroutine dump_opal_opac( iu_out, cf_d )
c     =========================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      character*(*) cf_d
c
!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/c_opal_ctrl_smooth/ init_smo, low_CO_smo, interp_CO_smo
      save /c_opal_ctrl_smooth/
c
      common /alt_change_opal_z/ main_alt_change, iulow, khighz_in,
     $     ofebrack_in
      save /alt_change_opal_z/
c
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
      parameter ( mx_hi=2*mx, mo_m1=mo-1, mo_m2=mo-2 )
c
      common /xhi_opal_z/ xhi_in(mx_hi), xcno_use(mx_hi,nz),
     $     xhi_use(mx_hi,nz), xxx_cno(mx_hi), xxx_hi(mx_hi),
     $     nx_hi(nz), ireq_hi(mx_hi), khighx(nz), kavail_xhi, kuse_xhi,
     $     kdo_xhi, kavail_cno, kuse_cno, kdo_cno, kavail_user,
     $     kuse_user, kdo_user
      save /xhi_opal_z/
c
!!      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
!!      parameter ( nel_zmix=19, n_zmixes=7, kel_o=3, kel_fe=nel_zmix-1 )
      parameter ( nel_zmix=19, n_zmixes=8, kel_o=3, kel_fe=nel_zmix-1 )
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
      parameter ( n_totmix = n_zmixes + 5, n_cnobeg = n_zmixes + 1 )
c
      character*255 cfile_opalGS98(n_totmix)
      common /opalGS98mixes/ bracketofe_opalGS98(n_totmix),
     $     xofe_opalGS98(n_totmix),xiz_opalGS98(nel_zmix,n_totmix),
     $     fninz_opalGS98(nel_zmix,n_totmix),atwt_opalGS98(nel_zmix),
     $     cfile_opalGS98
      save /opalGS98mixes/
c
      character*8 cdef_CNO_ext(n_cnobeg:n_totmix)
      common /ext_CNO_opal_z/ len_def_CNO_ext(n_cnobeg:n_totmix),
     $     cdef_CNO_ext
      save /ext_CNO_opal_z/
c
      common /cno_delta_opal_z/ fcno_mul(4), fninz_cno(nel_zmix,5),
     $     xiz_cno(nel_zmix,5), d_fninz_user(nel_zmix),
     $     fcno_fac(0:3,4), fninz_heavy, xiz_heavy, d_fninz_u_heavy,
     $     s_ninzai_mix, ds_ninzai_u, fn_o_over_cno, fninz_co_mix
      save /cno_delta_opal_z/
c
      common/b_opal_z/ nta(0:nrm_p1),ntax0(0:nrm),
     $     ntax03(0:nrm), sltlo, slthi, dltlo_inv, dlthi_inv,
     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
      save /b_opal_z/
c
      parameter ( zdel=0.001, xdel=0.03, xdelmin=0.001 )
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
c
      common/alink_opal_z/ NTEMP,NSM,nrlow,nrhigh,RLE,t6arr(100),
     $     coff(100,nrm)
      save /alink_opal_z/
c
      character*4 cxfil(5),czfil(mz)
      common/czinte_opal_z/ cxfil,czfil
      save /czinte_opal_z/
c
      character*255 copdir
      common/opdir/ copdir
      save /opdir/
c===
      if ( itime .ne. 12345678 ) stop
     $     ' STOP -- DUMP_OPAL_OPAC: Error: no opacities to dump. '
c
      last = lnblnk( cf_d )
c
      if ( last .le. 0 ) stop
     $     ' STOP -- DUMP_OPAL_OPAC: Error: blank dumpfile name. '
c
      if ( iu_out .gt. 0 .and. iu_out .le. 99 .and.
     $     iu_out .ne. 5 .and. iu_out .ne. 6 ) then
         iu = iu_out
      else
         iu = iulow
      endif
c
      call opneuf( iu, cf_d )
c
      write(iu) init_smo, low_CO_smo, interp_CO_smo,
     $     main_alt_change, iulow, khighz_in,
     $     ofebrack_in, level_err, nz, mx
c
      write(iu) indx,t6list,alr,n
      write(iu) alt,dfs,dfsr,b,m,mf,xa,alrf
      write(iu) flogtin,dfsx,oxf,cxf
      write(iu) xcdf,xodf,itime,cxdf
      write(iu) oxdf,q,h,xcd,xod,xc,xo
      write(iu) xcs,xos,cxd,oxd,cx,ox,zzz,xxh
      write(iu) xx,nc,no,zsto,zvint,dfsz,zacc
      write(iu) zlow,zmiddle,zhigh,zlo_ex,zhi_ex, numz
c
      do kz = 1, numz
         do kr = 1, nr
            do kt = 1, nt
               write(iu) ( ( ( co(ix,ic,io,kt,kr,kz), ix = 1, mx ),
     $              ic = 1, mc ), io = 1, mo )
            enddo
         enddo
      enddo
c
      write(iu) xhi_in, xcno_use,
     $     xhi_use, xxx_cno, xxx_hi,
     $     nx_hi, ireq_hi, khighx, kavail_xhi, kuse_xhi,
     $     kdo_xhi, kavail_cno, kuse_cno, kdo_cno, kavail_user,
     $     kuse_user, kdo_user
c
      write(iu) xiz_mix,fninz_mix,
     $     bracketife_mix,bracketofe_opalmixes,
     $     xofe_opalmixes,xiz_opalmixes,
     $     fninz_opalmixes,
     $     cel_opalmixes,cfile_opalmixes
c
      write(iu) bracketofe_opalGS98,
     $     xofe_opalGS98,xiz_opalGS98,
     $     fninz_opalGS98,atwt_opalGS98,
     $     cfile_opalGS98
c
      write(iu) len_def_CNO_ext,
     $     cdef_CNO_ext
c
      write(iu) fcno_mul, fninz_cno,
     $     xiz_cno, d_fninz_user,
     $     fcno_fac, fninz_heavy, xiz_heavy, d_fninz_u_heavy,
     $     s_ninzai_mix, ds_ninzai_u, fn_o_over_cno, fninz_co_mix
c
      write(iu) nta, ntax0,
     $     ntax03, sltlo, slthi, dltlo_inv, dlthi_inv,
     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
c
      write(iu) zvalhi,nofz,mnofz,
     $     zval,zalval,zavail,iadd_zavail
c
      write(iu) itimeco,mxzero,mx03,kope,igznotgx
c
      write(iu) NTEMP,NSM,nrlow,nrhigh,RLE,t6arr
c
      write(iu) cxfil,czfil, copdir
c
      close(iu)
c
      return
      end
c
c******************************************************************************
c
      subroutine read_opal_dump( iu_in, cf_d )
c     ========================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      character*(*) cf_d
c
!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/c_opal_ctrl_smooth/ init_smo, low_CO_smo, interp_CO_smo
      save /c_opal_ctrl_smooth/
c
      common /alt_change_opal_z/ main_alt_change, iulow, khighz_in,
     $     ofebrack_in
      save /alt_change_opal_z/
c
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
      parameter ( mx_hi=2*mx, mo_m1=mo-1, mo_m2=mo-2 )
c
      common /xhi_opal_z/ xhi_in(mx_hi), xcno_use(mx_hi,nz),
     $     xhi_use(mx_hi,nz), xxx_cno(mx_hi), xxx_hi(mx_hi),
     $     nx_hi(nz), ireq_hi(mx_hi), khighx(nz), kavail_xhi, kuse_xhi,
     $     kdo_xhi, kavail_cno, kuse_cno, kdo_cno, kavail_user,
     $     kuse_user, kdo_user
      save /xhi_opal_z/
c
!!      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
!!      parameter ( nel_zmix=19, n_zmixes=7, kel_o=3, kel_fe=nel_zmix-1 )
      parameter ( nel_zmix=19, n_zmixes=8, kel_o=3, kel_fe=nel_zmix-1 )
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
      parameter ( n_totmix = n_zmixes + 5, n_cnobeg = n_zmixes + 1 )
c
      character*255 cfile_opalGS98(n_totmix)
      common /opalGS98mixes/ bracketofe_opalGS98(n_totmix),
     $     xofe_opalGS98(n_totmix),xiz_opalGS98(nel_zmix,n_totmix),
     $     fninz_opalGS98(nel_zmix,n_totmix),atwt_opalGS98(nel_zmix),
     $     cfile_opalGS98
      save /opalGS98mixes/
c
      character*8 cdef_CNO_ext(n_cnobeg:n_totmix)
      common /ext_CNO_opal_z/ len_def_CNO_ext(n_cnobeg:n_totmix),
     $     cdef_CNO_ext
      save /ext_CNO_opal_z/
c
      common /cno_delta_opal_z/ fcno_mul(4), fninz_cno(nel_zmix,5),
     $     xiz_cno(nel_zmix,5), d_fninz_user(nel_zmix),
     $     fcno_fac(0:3,4), fninz_heavy, xiz_heavy, d_fninz_u_heavy,
     $     s_ninzai_mix, ds_ninzai_u, fn_o_over_cno, fninz_co_mix
      save /cno_delta_opal_z/
c
      common/b_opal_z/ nta(0:nrm_p1),ntax0(0:nrm),
     $     ntax03(0:nrm), sltlo, slthi, dltlo_inv, dlthi_inv,
     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
      save /b_opal_z/
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
c
      common/alink_opal_z/ NTEMP,NSM,nrlow,nrhigh,RLE,t6arr(100),
     $     coff(100,nrm)
      save /alink_opal_z/
c
      character*4 cxfil(5),czfil(mz)
      common/czinte_opal_z/ cxfil,czfil
      save /czinte_opal_z/
c
      character*255 copdir
      common/opdir/ copdir
      save /opdir/
c___
      logical lxst
c===
      last = lnblnk( cf_d )
c
      if ( last .le. 0 ) stop
     $     ' STOP -- READ_OPAL_DUMP: Error: blank dumpfile name. '
c
      call inqfil( cf_d, lxst )
c
      if ( .not. lxst ) then
         write(6,'(" READ_OPAL_DUMP: dumpfile not found:"/" ",a)')
     $        cf_d(1:last)
         stop ' STOP -- READ_OPAL_DUMP: Error: dumpfile not found. '
      endif
c
      if ( iu_in .gt. 0 .and. iu_in .le. 99 .and.
     $     iu_in .ne. 5 .and. iu_in .ne. 6 ) then
         iu = iu_in
      else
         iu = iulow
      endif
c
      call opoluf( iu, cf_d )
c
      read(iu) init_smo, low_CO_smo, interp_CO_smo,
     $     main_alt_change, iulow, khighz_in,
     $     ofebrack_in, level_err, nz_sto, mx_sto
c
      if ( nz .ne. nz_sto ) then
         close(iu)
         write(6,
     $        '(" READ_OPAL_DUMP: nz=",i3," .ne. nz=",i3," in"/" ",a)')
     $        nz, nz_sto, cf_d(1:last)
         stop ' STOP -- READ_OPAL_DUMP: Error: bad parameter nz. '
      else if ( mx .ne. mx_sto ) then
         close(iu)
         write(6,
     $        '(" READ_OPAL_DUMP: mx=",i3," .ne. mx=",i3," in"/" ",a)')
     $        mx, mx_sto, cf_d(1:last)
         stop ' STOP -- READ_OPAL_DUMP: Error: bad parameter mx. '
      endif
c
      read(iu) indx,t6list,alr,n
      read(iu) alt,dfs,dfsr,b,m,mf,xa,alrf
      read(iu) flogtin,dfsx,oxf,cxf
      read(iu) xcdf,xodf,itime,cxdf
      read(iu) oxdf,q,h,xcd,xod,xc,xo
      read(iu) xcs,xos,cxd,oxd,cx,ox,zzz,xxh
      read(iu) xx,nc,no,zsto,zvint,dfsz,zacc
      read(iu) zlow,zmiddle,zhigh,zlo_ex,zhi_ex, numz
c
      do kz = 1, numz
         do kr = 1, nr
            do kt = 1, nt
               read(iu) ( ( ( co(ix,ic,io,kt,kr,kz), ix = 1, mx ),
     $              ic = 1, mc ), io = 1, mo )
            enddo
         enddo
      enddo
c
      read(iu) xhi_in, xcno_use,
     $     xhi_use, xxx_cno, xxx_hi,
     $     nx_hi, ireq_hi, khighx, kavail_xhi, kuse_xhi,
     $     kdo_xhi, kavail_cno, kuse_cno, kdo_cno, kavail_user,
     $     kuse_user, kdo_user
c
      read(iu) xiz_mix,fninz_mix,
     $     bracketife_mix,bracketofe_opalmixes,
     $     xofe_opalmixes,xiz_opalmixes,
     $     fninz_opalmixes,
     $     cel_opalmixes,cfile_opalmixes
c
      read(iu) bracketofe_opalGS98,
     $     xofe_opalGS98,xiz_opalGS98,
     $     fninz_opalGS98,atwt_opalGS98,
     $     cfile_opalGS98
c
      read(iu) len_def_CNO_ext,
     $     cdef_CNO_ext
c
      read(iu) fcno_mul, fninz_cno,
     $     xiz_cno, d_fninz_user,
     $     fcno_fac, fninz_heavy, xiz_heavy, d_fninz_u_heavy,
     $     s_ninzai_mix, ds_ninzai_u, fn_o_over_cno, fninz_co_mix
c
      read(iu) nta, ntax0,
     $     ntax03, sltlo, slthi, dltlo_inv, dlthi_inv,
     $     slrlo, slrhi, dlrlo_inv, dlrhi_inv, init_trvals
c
      read(iu) zvalhi,nofz,mnofz,
     $     zval,zalval,zavail,iadd_zavail
c
      read(iu) itimeco,mxzero,mx03,kope,igznotgx
c
      read(iu) NTEMP,NSM,nrlow,nrhigh,RLE,t6arr
c
      read(iu) cxfil,czfil, copdir
c
      close(iu)
c
      return
      end
c
c******************************************************************************
c
      subroutine read_kz(kz,z,kallrd,khighz,iu_lo,ofebrack)
c     =====================================================
c
c NOTE: kallrd is ignored (it is present only for backward compatibility).
c
c PARAMETERS to control the offset from zero for the logarithmic interpolation:
c  zdel = 0.001 (must be the same value as in OPAL)
c  xdel = 0.03 (must be the same value as in OPAL)
c  xdelgn93 = 0.005 = xdel value for use with X-interpolation in 'GN93hz' file
c                      among X-tables 0.0, 0.1, 0.2 (to get X = 0.03 mix); note
c                      that 0.005 works slightly better for this than 0.03
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
C       implicit none
      
C       integer i,j,nxdo,nmorex,j1,j2,j3,j4,is,isx,iu,moat,iofe,itab_dum,
C      &     kzsto,kx,kat,i_rewind,ifound,mixfind,k,il,kdelzhi,jl,is1 ,iz
C      &     ,jz,nomj,khi,ihi,im2,im1,jm2,jm1,nmorat,imix,ifac,jfac, il_lo
C      &     ,il_hi,lnblnk,itab_dum_ofe, iu_ofe,jo,ko,kc,jc,jx,jjx,jjc,jjo
C      &     ,kko,kkx,kkc,iset,k_e,igzip
C       integer kz,kallrd,khighz,iu_lo,nz,mx,mc,mo,nrm,ntm,nrb,nre,ntb,nr,
C      &     nt,nrm_p1,mx_hi,mo_m1,mo_m2,nx_hi,ireq_hi,khighx,kavail_xhi,
C      &     kuse_xhi, kdo_xhi, kavail_cno, kuse_cno, kdo_cno, kavail_user
C      &     ,kuse_user, kdo_user,ks81,ks83,ks60,ks61,indx,n,m,mf,int_hi_z
C       integer nrdel, ntdel,nrm_m2,nt_m1,nre_p1,nre_m1,numz,last
C       integer nel_zmix, n_zmixes, kel_o, kel_fe,mstore,ix,kfil
C       integer  n_totmix,n_cnobeg,nc,no,mofe,ilo,mx_use,iz_hi,io
C       integer nta,ntax0,ntax03,init_trvals,init_smo,interp_CO_smo
C      &     ,low_CO_smo,len_def_CNO_ext,kdel
C       integer NTEMP,NSM,nrlow,nrhigh,NRL,RLS,nset,itimeco,mxzero,mx03
C      &     ,kope,igznotgx
C       integer mz, mz_m1, mz_m2, mzhi, mzal,nzm, nadd_zavail
C       integer nofz,mnofz,iadd_zavail,igetzxi,int_hi_1,int_hi_2
C       integer main_alt_change, iulow, khighz_in,level_err,int_hi_3
C       integer kzvalnof,nofmjz,iz_get,iz_rew,iz_tab,ico_got,nthx0,nthx03
C       integer igzip_sto,line,iulow,kzbelo,kzlow,k1zfac,nmorez,kzlow_m1,
C      &     k2zfac,kznof,lznof,khighz_index,khighz_cno,klozat,khizat
C      &     ,kzalbe,kzalow, nzalmo
C       double precision zdel,xdel,xdelgn93,z,ofebrack,xhi_in,xcno_use
C      &     ,xxx_cno,xxx_hi,badlogkval,badlogklim,alrlo,flogtlo,flogt60,
C      &     flogt81,t6list,alr,alt,dfs,dfsx,dfsr,b,alrf,flogtin,oxf,cxf,
C      &     oxdf,q,h,xcd,xod,xc,xo,xcs,xos,cxd,oxd,cx,ox,zzz,xxh,zlow
C      &     ,zmiddle,zhigh,zlo_ex,zhi_ex,co,opk,opl,cof,qzinter,oat,cat,
C      &     cxhi,oxhi,a1,a2,a3,a4,v1,v2,v3,v4,cofmin,cofmax,sum_niai,dif,
C      &     difl,diffac,cget
C       double precision xiz_mix,fninz_mix,bracketife_mix
C      &     ,bracketofe_opalmixes,xofe_opalmixes,xiz_opalmixes
C      &     ,fninz_opalmixes,xhemx,xhi_use,xx,xxx_at
C       double precision bracketofe_opalGS98,xofe_opalGS98,oget,
C      &     xiz_opalGS98,fninz_opalGS98,atwt_opalGS98
C       double precision sltlo, slthi, dltlo_inv, dlthi_inv,slrlo, slrhi,
C      &     dlrlo_inv, dlrhi_inv,rle,t6arr,coff,tmax,opact,dopact,dopacr
C      &     ,dopactd,fedge,ftredge,fzedge
C       double precision zvalhi,zval,zalval,zavail, ofebrack_in,cofzat,
C      &     xzalat,cofzhi,xhi_look,zhi_look,zacc,zat,facxhz,xofe,fofe
C      &     ,omfofe,cof_tmp,xxx_max,xxx_03,f_3,omf_3

c
      parameter ( zdel=0.001d0, xdel=0.03d0, xdelgn93=0.005d0 )
c
!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      parameter ( mx_hi=2*mx, mo_m1=mo-1, mo_m2=mo-2 )
c
      common /xhi_opal_z/ xhi_in(mx_hi), xcno_use(mx_hi,nz),
     $     xhi_use(mx_hi,nz), xxx_cno(mx_hi), xxx_hi(mx_hi),
     $     nx_hi(nz), ireq_hi(mx_hi), khighx(nz), kavail_xhi, kuse_xhi,
     $     kdo_xhi, kavail_cno, kuse_cno, kdo_cno, kavail_user,
     $     kuse_user, kdo_user
      save /xhi_opal_z/
c
      parameter ( nrdel=nrb-1, ntdel=ntb-1 )
      parameter ( nrm_m2=nrm-2, nt_m1=nt-1, nre_p1=nre+1, nre_m1=nre-1 )
      parameter ( badlogkval=1.d+35, badlogklim=20.d0 )
      parameter ( ks81=ntm-3, ks83=ks81+1, ks60=ks81-21, ks61=ks60+1,
     $     alrlo=-9.0,d0 flogtlo=3.75d0, flogt60=6.6,d0 flogt81=9.3d0 )
!!     $     alrlo=-8.0, flogtlo=3.75, flogt60=6.0, flogt81=8.1 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
!!      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
!!      parameter ( nel_zmix=19, n_zmixes=7, kel_o=3, kel_fe=nel_zmix-1 )
      parameter ( nel_zmix=19, n_zmixes=8, kel_o=3, kel_fe=nel_zmix-1 )
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
      parameter ( n_totmix = n_zmixes + 5, n_cnobeg = n_zmixes + 1 )
c
      character*255 cfile_opalGS98(n_totmix)
      common /opalGS98mixes/ bracketofe_opalGS98(n_totmix),
     $     xofe_opalGS98(n_totmix),xiz_opalGS98(nel_zmix,n_totmix),
     $     fninz_opalGS98(nel_zmix,n_totmix),atwt_opalGS98(nel_zmix),
     $     cfile_opalGS98
      save /opalGS98mixes/
c
      character*8 cdef_CNO_ext(n_cnobeg:n_totmix)
      common /ext_CNO_opal_z/ len_def_CNO_ext(n_cnobeg:n_totmix),
     $     cdef_CNO_ext
      save /ext_CNO_opal_z/
c
      common/b_opal_z/ nta(0:nrm_p1),ntax0(0:nrm),
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
      common /alt_change_opal_z/ main_alt_change, iulow, khighz_in,
     $     ofebrack_in
      save /alt_change_opal_z/
c
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c
      character*4 cxfil(5),czfil(mz)
      common/czinte_opal_z/ cxfil,czfil
      save /czinte_opal_z/
c
c  The following common block /opdir/ contains a character file giving the
c  directory where the Gz???.x?? files are to be found.
c
      character*255 copdir
      common/opdir/ copdir
      save /opdir/
c___
c
c  copfil = the full opacity filename (including directory), changed as needed
c  copalt = temporary (needed only for error message: new form of opacity file)
c
      character*255 copfil,copalt
c
c  cin / ch  holds a line read in from the opacity file (with opacity values)
c
c      character*137 cin
c      character*1 ch(137)
      character*277 cin
      character*1 ch(277)
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
c Storage for the compression-flags and opacity file names that are opened,
c and the line numbers in these opacity files:
c
      dimension igzip_sto(0:3)
      character*255 cop_sto(0:3)
      dimension line(0:3)
c
c-debug-chk[
c-debug-chk;      dimension chk_max(20),chk_min(20),chk_sum(20),chk_ssq(20),
c-debug-chk;     $     n_chk(20)
c-debug-chk;c
c-debug-chk;      common /readkz_opal_debug_chk/ iout_debug_chk_ofe
c-debug-chk]
c-test-xdel[                              ! test xdel values in GN93hz X-interp
c-test-xdel;      parameter ( n_xdtst=4 )
c-test-xdel;      dimension cof_tst(ntm,nrm,n_xdtst),dif_tst(4,n_xdtst),
c-test-xdel;     $     xdel_tst(n_xdtst)
c-test-xdel;      data xdel_tst/0.03,0.01,0.001,0.0001/
c-test-xdel]
c-debug-chk[
c-debug-chk;      data iout_debug_chk_ofe / 99999 /
c-debug-chk]
c
c===


      if ( kz .le. 0 .or. kz .gt. nz ) stop
     $     ' STOP -- READCO: Error: Z-index out of range: cannot be! '
c
c Check input unit number iu_lo, and use it to set iulow.
c
      if ( iu_lo .ge. 7 .and. iu_lo .le. 96 ) then
         iulow = iu_lo
      else if ( iu_lo .ge. 0 .and. level_err .gt. 0 ) then
         stop ' STOP -- READCO: Error: bad iulow value. '
      endif
c
      if ( level_err .gt. 99 .and. kallrd .le. 0 ) write(6,
     $     '(" ***WARNING: READCO: kallrd < 1 ignored.")')
c
c Check some parameter values, and do some initializations.
c
      call opalinit( khighz, ofebrack, z, kz )
c
      khighz_in = khighz
      ofebrack_in = ofebrack
c							     ! Z out of range?
      if ( z .ge. zval(mz) + 1.d-6 .or. z .le. -1.d-8 ) then
         write(6,8703) z
 8703    format(' '/' STOP -- READCO: Z=',f10.7,
     $        ': Z > 0.1 (or < 0) not allowed!!!')
         stop
      endif
c						      ! accuracy to match table
      zacc(kz) = min( 1.d-6 , max( 1.d-8 , 0.01d0 * z ) )
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
         kzvalnof(i-kzlow_m1) = int( 100.d0 * zval(i) + 1.01d0 )
      enddo
      kznof = min( int( 100.d0 * z + 1.0001d0 ) , mzhi )
      lznof = max( int( 100.d0 * z + 0.9999d0 ) , 1 )
c
c  NOTE that if Z > 0 and the GS98 mix is being used (khighz < 0), then both
c  the GN93hz and GS98hz files will always be needed; if [O/Fe] is non-zero
c  and khighz < -1, then a GS98 [O/Fe]-file is also needed.  OTHERWISE:
c  Check if need C+O=0.0 "hz" tables: if [O/Fe] = 0 or Z = 0, there will be no
c  need to interpolate in [O/Fe] (set khizat = 1), and if there is no need to
c  interpolate in [O/Fe], then for Z equal to a Z-table value or .01 < Z < .02
c  the "hz" tables will yield no improvement and are not used (set khizat = 0).
c
      khighz_index = min( mod( abs(khighz), 10 ) , n_zmixes )
      if ( mx .eq. 5 ) then
         khighz_cno = min( abs(khighz) / 10 , 3 )
         if ( khighz_cno .ge. 2 .and. kavail_user .le. 0 )
     $        khighz_cno = khighz_cno - 2
         if ( mod( khighz_cno, 2 ) .eq. 1 .and. kavail_cno .le. 0 )
     $        khighz_cno = khighz_cno - 1
      else
         khighz_cno = 0
      endif
c
      if ( z .lt. zacc(kz) ) then
         klozat = 0
         khizat = 0
      else if ( khighz .lt. 0 ) then
         khizat = 1
         if ( ofebrack .eq. 0.d0 ) then
            klozat = 1
         else
            klozat = khighz_index
         endif
      else
         klozat = 0
         khizat = khighz_index
         if ( ofebrack .eq. 0.d0 ) khizat = min(khizat,1)
         if ( khizat .eq. 1 .and. ( ( zval(k1zfac) .ge. 0.01d0 .and.
     $        zval(k2zfac) .le. 0.02d0 ) .or. nmorez .eq. 0 ) )
     $        khizat = 0
      endif
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
         copfil = copdir(:kope)
      endif
c				! store present m-value (should be unnecessary)
      mstore = m
c
      khighx(kz) = 0
      if ( khighz_index .ne. 0 .and. mx .eq. 5 .and.
     $     max( abs( xhi_in(1) - xa(1) ) , abs( 0.03d0 - xa(2) ) ,
     $     abs( xhi_in(2) - xa(3) ) , abs( xhi_in(4) - xa(4) ) ,
     $     abs( xhi_in(6) - xa(5) ) ) .le. 1.d-6 ) khighx(kz) = 1
c
c			   ! get shifted z-value Z+zdel, for log interpolation
      zzz(kz) = zat + zdel
      zvint(kz) = log10( zzz(kz) )
      zsto(kz) = zat
c			 ! should read in Z-composition from first opacity file
      igetzxi = 1
c
c {--------------------------------------- Begin loop over m values (X-tables):
c
      do m = 1, mx
c					     ! later Z-compositions: just check
         if ( m .ne. 1 ) igetzxi = -1
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
c
c........ initialization: for itime, oxf...oxdf, n(m,*,kz)  (xx was done above)
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
c					     ! read 'GN93hz' only if necessary:
c
         if ( khizat .gt. 0 .and. ( m .ne. mx03 .or.
     $        khizat .gt. 1 .or. klozat .gt. 0 ) ) then
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
     $              facxhz = 1.d0- 100.d0*min( zat - 0.04d0,0.06d0-zat)
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
            moat = mo
            iofe = 1
c						! read C+O=0.0 "hz"-table(s):
            do while ( iofe .ne. 0 )
c							! get filename
               if ( iofe .gt. 0 ) then
                  copfil(kope+1:) = cfile_opalmixes(iofe)
               else
                  copfil(kope+1:) = cfile_opalGS98(-iofe)
               endif
c								! open file
               call open_chk_zip( iu, copfil, igzip,
     $              'STOP -- Error: hz-file (C+O=0.0) not found.' )
c
               line(1) = 0
c				     ! dummy table-number; initial cofzhi index
               itab_dum = 0
               kzsto = 0
c					! get Z-composition(s) from "hz"-files,
               if ( m .eq. 1 ) then
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
     $                    itab_dum,line(1),
     $                    zalval(iz),xzalat(kx),0.0d0,0.0d0)
                     if ( ifound .eq. 0 ) then
                        write(6,1791) zalval(iz),xzalat(kx),0.0d0,0.0d0,
     $                          copfil(:lnblnk(copfil))
 1791                   format(' '/' READCO: Error finding mix Z=',
     $                       f6.4,' X=',f6.4,' C=',f6.4,' O=',f6.4,
     $                       '  from file:'/' ',a/' ')
                        stop ' STOP -- READCO: error reading hz-mix. '
                     endif
c								! check [O/Fe]
                     if ( kx .eq. 1 .and. iz .eq. kzalow ) then
                        if ( iofe .eq. 1 ) then
c						! (this cannot happen)
                           if ( abs( bracketofe_opalmixes(1) ) .gt.
     $                          3.e-7 ) write(6,1783) 'GN93hz',
     $                          bracketofe_opalmixes(1)
 1783                      format(/' READCO: Error: ',a,
     $                          ' file has [O/Fe] =',f12.7/20x,
     $                          '(its [O/Fe] should be 0.0)',5x,
     $                          'THIS SHOULD NOT HAPPEN!'/)
                           if ( abs( bracketofe_opalmixes(1) ) .gt.
     $                          3.e-7 ) stop
     $                          ' STOP -- READCO: non-0 [O/Fe]_GN93hz '
c						! is file [O/Fe] too small?
                        else if ( iofe .gt. 0 ) then
                           if ( abs(bracketofe_opalmixes(iofe)) .lt.
     $                          max(0.1d0*abs(ofebrack),0.001d0) ) then
                              write(6,2631) ofebrack,
     $                             bracketofe_opalmixes(iofe),
     $                             iofe, copfil(:lnblnk(copfil))
 2631                         format(' '/' READCO: [O/Fe] =',f10.6,
     $                             ': cannot get from [O/Fe] =',
     $                             f10.6,' in file',i3,':'/' ',a/' ')
                              stop ' STOP -- READCO: bad [O/Fe] file. '
                           endif
                        else if ( iofe .eq. -1 ) then
c						! (this cannot happen)
                           if ( abs( bracketofe_opalGS98(1) ) .gt.
     $                          3.d-7 ) write(6,1783)
     $                          'alternate-solar (e.g., GS98hz)',
     $                          bracketofe_opalGS98(1)
                           if ( abs( bracketofe_opalGS98(1) ) .gt.
     $                          3.d-7 ) stop
     $                          ' STOP -- READCO: non-0 [O/Fe]_GS98hz '
c						! is file [O/Fe] too small?
                        else 
                           if ( abs(bracketofe_opalGS98(-iofe)) .lt.
     $                          max(0.1d0*abs(ofebrack),0.001d0) ) then
                              write(6,2631) ofebrack,
     $                             bracketofe_opalGS98(-iofe),
     $                             iofe, copfil(:lnblnk(copfil))
                              stop ' STOP -- READCO: bad [O/Fe] file. '
                           endif
                        endif
                     endif
c				  ! loop over logT values, to read in opacities
                     do k = 1, ntm
c					      ! read logT,{logKappa(R) @ all R}
                        line(1) = line(1) + 1
                        read(iu,7300) cin
c 7300                   format(a137)
 7300                   format(a277)
                        read(cin,7140) flt, (cofzhi(k,il,kat),il=1,nrm)
!! 7140                   format(f4.2,19f7.3)
 7140                   format(f4.2,39f7.3)
c								   ! bad logT ?
                        if ( abs(flogtin(k)-flt) .gt. 1.d-5 ) then
                           write(6,1734) flt, flogtin(k),
     $                          copfil(:lnblnk(copfil)), line(1),
     $                          cin(:max(1,lnblnk(cin))),
     $                          zalval(iz),xzalat(kx),0.0d0,0.0d0
 1734                      format(/' Error reading logT value =',f10.6,
     $                          ' should be',f10.6,
     $                          ' from opacity file:'/' ',a/
     $                          ' at line',i8,', which contained:'/
     $                          ' "',a,'"'/
     $                          ' while reading mix [Z=',f6.4,
     $                          ' X=',f6.4,' C=',f6.4,' O=',f6.4,']'/
     $                          ' *****THIS SHOULD NOT HAPPEN.'/)
                           stop ' STOP -- READCO: bad logT value. '
                        endif
c							      ! logKappa(R) is:
                        do il = nrm, 1, -1
c								       ! absent
                           if ( cin(7*il-2:7*il+4) .eq.
     $                          '       ' ) then
                              if ( k .le. max(nta(il),nta(0)) ) stop
     $                             ' STOP -- READCO: bad upper edge. '
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
c				! store present hz-opacity set in free space
               do il = 1,nr
                  jl = il+nrdel
                  do k = 1,nt
                     co(m,mc,moat,k,il,kz) = coff(k+ntdel,jl)
                  enddo
               enddo
c				! position to store next hz-opacity set
               moat = moat - 1
c						! get next iofe value, if any
               if ( iofe .eq. -klozat ) then
                  iofe = 0
               else if ( iofe .eq. -1 ) then
                  iofe = -klozat
               else if ( iofe .eq. khizat ) then
                  iofe = max( -klozat , -1 )
               else
                  iofe = khizat
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
c               copfil(kope+1:kope+1) = 'G'
               copfil(kope+1:kope+1) = 'A'
               copfil(kope+2:kope+5) = czfil(kzlow+iu-iulow)
               if ( copfil(kope+4:kope+4) .eq. ' ' )
     $              copfil(kope+4:kope+4) = '0'
               if ( copfil(kope+5:kope+5) .eq. ' ' )
     $              copfil(kope+5:kope+5) = '0'
               copfil(kope+6:) = cxfil(m)
               copfil(kope+6:kope+6) = '.'
            else
               copfil(kope+1:kope+4) = cxfil(m)
               copfil(kope+5:) = czfil(kzlow+iu-iulow)
            endif
            if ( igznotgx .eq. 0 ) then
               call inqfil(copfil,lxst)
               if ( .not. lxst )
     $              call inqfil( copfil(1:kope+9) // '.gz' , lxst )
               if ( .not. lxst )
     $              call inqfil( copfil(1:kope+9) // '.Z' , lxst )
               if ( lxst ) then
                  igznotgx = 1
               else
                  copalt = copfil
                  copfil(kope+1:kope+4) = cxfil(m)
                  copfil(kope+5:) = czfil(kzlow+iu-iulow)
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
                     write(6,7399) copalt(:lnblnk(copalt)),
     $                    copfil(:lnblnk(copfil)),copal,copfil
!! 7399                format(' '/' STOP -- READCO: neither Gz???.x??',
!!     $                    ' nor Gx??z* OPAL opacity files found:'/
!! 7399                format(' '/' STOP -- READCO: neither AGS0z???.x??',
!!     $                    ' nor AGSx??z* OPAL opacity files found:'/
 7399                format(' '/' STOP -- READCO: neither AGSz???.x??',
     $                    ' nor AGSS09x??z* OPAL opacity files found:'/
     $                    ' ',a/' ',a/' ',a/' ',a)
!!                     stop ' STOP -- READCO: Gz???.x?? file not found. '
!!                     stop ' STOP -- READCO: AGSz???.x?? file not found.'
                     stop ' STOP -- READCO: A09z???.x?? file not found.'
                  endif
               endif
            endif
            cop_sto(iu-iulow) = copfil
            call open_chk_zip( iu, copfil, igzip_sto(iu-iulow),
     $        'STOP -- Error: AGSS09z???.x?? opacity file not found.' )
!!     $           'STOP -- Error: AGSz???.x?? opacity file not found.' )
!!     $           'STOP -- Error: Gz???.x?? opacity file not found.' )
            line(iu-iulow) = 0
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
     $                 iz_tab(jz),line(iu-iulow),
     $                 zval(iz),xa(m),cget,xos(j))
c							      ! if table is not
c					      ! present in this file, it cannot
c					      ! be used in the Z-interpolation,
                  if ( ifound .eq. 0 ) then
c							     ! so reduce nmorat
                     nmorat = min( iz - kzlow - 1 , nmorat )
c							     ! (cannot happen):
                     if ( nmorat .lt. 0 ) then
                        write(6,1791) zval(iz),xa(m),cget,xos(j),
     $                       cop_sto(iu
     $                       -iulow)(:lnblnk(cop_sto(iu-iulow)))
                        stop ' STOP -- READCO: error reading mix. '

c-debug-chk[                     else
c-debug-chk;                        write(6,1878) m,z,i,j,jz,nofmjz(jz),
c-debug-chk;     $                       iz_rew(jz),nmorat
c-debug-chk; 1878                   format(' m=',i1,' Z=',f9.7,
c-debug-chk;     $                       ' cannot find i=',i1,' j=',i1,
c-debug-chk;     $                       ' Z(jz=',i1,'): Nofmjz=',i1,
c-debug-chk]     $                       ' irew=',i2,' nmorat=',i2)

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
                        line(iu-iulow) = line(iu-iulow) + 1
                        read(iu,7300) cin
                        read(cin,7140) flt,(cofzat(jz,il),il=1,nrm)
c								   ! bad logT ?
                        if ( abs(flogtin(k)-flt) .gt. 1.d-5 ) then
                           write(6,1734) flt, flogtin(k), cop_sto(iu
     $                          -iulow)(:lnblnk(cop_sto(iu-iulow))),
     $                          line(iu-iulow),
     $                          cin(:max(1,lnblnk(cin))),
     $                          zval(iz),xa(m),cget,xos(j)
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
c
            nmorat = nmorez
            is1 = 1
            nomj = no - j
c					     ! loop over Z values: find tables
            do iz = kzlow, kzlow + nmorat
               iu = iulow + iz - kzlow
               jz = iz - kzlow_m1
               oget = 1.d0 - xa(m) - xcs(nomj) - zval(iz)
               ifound = mixfind(iu,1,iz_get(jz),iz_rew(jz),
     $              iz_tab(jz),line(iu-iulow),
     $              zval(iz),xa(m),xcs(nomj),oget)
               if ( ifound .eq. 0 ) then
                  nmorat = min( iz - kzlow - 1 , nmorat )
                  if ( nmorat .lt. 0 ) then
                     write(6,1791) zval(iz),xa(m),xcs(nomj),oget,
     $                    cop_sto(iu-iulow)(:lnblnk(cop_sto(iu-iulow)))
                     stop ' STOP -- READCO: error reading mix. '

c-debug-chk[                  else
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
                     line(iu-iulow) = line(iu-iulow) + 1
                     read(iu,7300) cin
                     read(cin,7140) flt,(cofzat(jz,il),il=1,nrm)
c								   ! bad logT ?
                     if ( abs(flogtin(k)-flt) .gt. 1.d-5 ) then
                        write(6,1734) flt, flogtin(k), cop_sto(iu
     $                       -iulow)(:lnblnk(cop_sto(iu-iulow))),
     $                       line(iu-iulow), cin(:max(1,lnblnk(cin))),
     $                       zval(iz),xa(m),xcs(nomj),oget
                        stop ' STOP -- READCO: bad logT value. '
                     endif
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
     $                          flt .ge. 4.2d0) then
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
     $                    ' STOP -- READCO: missing mix: CANNOT be. '
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
     $                    ' STOP -- READCO: missing mix: CANNOT be. '
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
     $                   zat,m,j1,mo,no-1,xcs(j1),1.d0-xa(m)-zat-xcs(j1)
                     if ( ico_got(j2,mo) .le. 0 ) write(6,4873)
     $                   zat,m,j2,mo,no-1,xcs(j2),1.d0-xa(m)-zat-xcs(j2)
                     if ( ico_got(j3,mo) .le. 0 ) write(6,4873)
     $                   zat,m,j3,mo,no-1,xcs(j3),1.d0-xa(m)-zat-xcs(j3)
 4873                format(' '/'     Z=',f9.7,' ??? CANNOT miss',
     $                    ' mix (m=',i1,',i=',i1,',j=',i1,
     $                    ') with ihi=',i1,' C=',f10.7,' O=',f10.7)
                     if ( ico_got(j1,mo) .le. 0 .or.
     $                    ico_got(j2,mo) .le. 0 .or.
     $                    ico_got(j3,mo) .le. 0 ) stop
     $                    ' STOP -- READCO: missing mix: CANNOT be. '
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
     $           ' STOP -- READCO: missing mix: CANNOT be. '
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
                        cofmin = min( 0.8d0*v1+0.2d0*v3,0.2d0*v1+0.8d0*
     &                       v3 )
                        cofmax = max( 0.8d0*v1+0.2d0*v3,0.2d0*v1+0.8d0*
     &                       v3 )
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
         if ( m .eq. 1 ) then
c				     ! if needed, get Z-composition of GS98 mix
            if ( klozat .eq. 1 ) then
c
               do i = 1, nel_zmix
                  xiz_mix(i) = xiz_opalGS98(i,1)
                  fninz_mix(i) = fninz_opalGS98(i,1)
                  bracketife_mix(i) = 0.d0
               enddo
c							       ! or of GN93 mix
            else if ( khizat .eq. 1 .and. klozat .eq. 0 ) then
c
               do i = 1, nel_zmix
                  xiz_mix(i) = xiz_opalmixes(i,1)
                  fninz_mix(i) = fninz_opalmixes(i,1)
                  bracketife_mix(i) = 0.d0
               enddo
c
            endif
c
         endif
c						    ! If m=2 shift & no [O/Fe]:
         if ( khizat .eq. 1 .and. m .eq. mx03 .and.
     $        klozat .eq. 0 ) then
c						! set all shifts to zero
c			  ! (m=2 GN93hz shift may  be interpolated later)
            do il = 1, nr
               do k = 1, nt
                  co(m,mc,mo_m1,k,il,kz) = 0.d0
                  co(m,mc,mo,k,il,kz) = 0.d0
               enddo
            enddo
c					       ! Else, if there are any shifts:
         else if ( khizat .gt. 0 ) then
c					       ! If there is no [O/Fe] shift:
c
            if ( khizat .eq. 1 .and. klozat .eq. 0 ) then
c							  ! then set it to zero
               do il = 1, nr
                  do k = 1, nt
                     co(m,mc,mo_m1,k,il,kz) = 0.d0
                  enddo
               enddo
c			! Else, if there is the [O/Fe] or the GS98-GN93 shift:
            else
c		! get interpolation factors fofe (for GN93hz) and omfofe=1-fofe
c					 ! & Z-composition of interpolated mix
               if ( klozat .gt. 1 ) then
c							  ! GS98 + [O/Fe] shift
                  xofe = 10.d0**ofebrack * xofe_opalGS98(1)
                  fofe = ( fninz_opalGS98(kel_o,klozat)
     $                 - xofe * fninz_opalGS98(kel_fe,klozat) )
     $                 / ( ( fninz_opalGS98(kel_fe,1)
     $                 - fninz_opalGS98(kel_fe,klozat) ) * xofe
     $                 + fninz_opalGS98(kel_o,klozat)
     $                 - fninz_opalGS98(kel_o,1) )
                  omfofe = 1.d0 - fofe
                  mofe = mo_m2
                  moat = mo_m1
c					! get Z-composition of interpolated mix
                  if ( m .eq. 1 ) then
                     sum_niai = 0.0d0
                     do i = 1, nel_zmix
                        fninz_mix(i) = fofe * fninz_opalGS98(i,1)
     $                       + omfofe * fninz_opalGS98(i,klozat)
                        xiz_mix(i) = fninz_mix(i) * atwt_opalGS98(i)
                        sum_niai = sum_niai + xiz_mix(i)
                     enddo
                     do i = 1, nel_zmix
                        xiz_mix(i) = xiz_mix(i) / sum_niai
                        bracketife_mix(i) = log10(
     $                       ( max( fninz_mix(i) , 1.d-36 )
     $                       * fninz_opalGS98(kel_fe,1) )
     $                       / ( max( fninz_mix(kel_fe) , 1.d-36 )
     $                       * fninz_opalGS98(i,1) ) )
                     enddo
c-debug-chk[
c-debug-chk;                     if ( iout_debug_chk_ofe .gt. 0 ) then
c-debug-chk;                        write(6,2377) ofebrack,-klozat,
c-debug-chk;     $                       bracketofe_opalGS98(klozat),
c-debug-chk;     $                       fofe,klozat,omfofe,klozat,klozat
c-debug-chk; 2377                   format(' '/' To get mix with [O/Fe] =',
c-debug-chk;     $                       f11.7,' from mix',i2,' with [O/Fe] =',
c-debug-chk;     $                       f11.7,': f_(1) =',f11.7,' , f_(',i1,
c-debug-chk;     $                       ') =',f11.7,':'/' '/
c-debug-chk;     $                       '  i   Xi/Z_(1)   Ni/Nz_(1)',
c-debug-chk;     $                       '   Xi/Z_(',i1,')   Ni/Nz_(',i1,')',
c-debug-chk;     $                       '   Xi/Z_mix   Ni/Nz_mix   [i/Fe]'/
c-debug-chk;     $                       ' ==  ========== ==========',
c-debug-chk;     $                       '  ========== ==========  ==========',
c-debug-chk;     $                       ' ========== ==========')
c-debug-chk;                        do i = 1,nel_zmix
c-debug-chk;                           write(6,2376) cel_opalmixes(i),
c-debug-chk;     $                          xiz_opalGS98(i,1),
c-debug-chk;     $                          fninz_opalGS98(i,1),
c-debug-chk;     $                          xiz_opalGS98(i,klozat),
c-debug-chk;     $                          fninz_opalGS98(i,klozat),
c-debug-chk;     $                          xiz_mix(i),fninz_mix(i),
c-debug-chk;     $                          bracketife_mix(i)
c-debug-chk; 2376                      format(' ',a2,3(f12.7,f11.7),f11.7)
c-debug-chk;                        enddo
c-debug-chk;                        write(6,'(" ")')
c-debug-chk;                        iout_debug_chk_ofe = iout_debug_chk_ofe - 1
c-debug-chk;                     endif
c-debug-chk]
                  endif
c
               else if ( khizat .gt. 1 ) then
c							   ! [O/Fe] shift only
                  xofe = 10.d0**ofebrack * xofe_opalmixes(1)
                  fofe = ( fninz_opalmixes(kel_o,khizat)
     $                 - xofe * fninz_opalmixes(kel_fe,khizat) )
     $                 / ( ( fninz_opalmixes(kel_fe,1)
     $                 - fninz_opalmixes(kel_fe,khizat) ) * xofe
     $                 + fninz_opalmixes(kel_o,khizat)
     $                 - fninz_opalmixes(kel_o,1) )
                  omfofe = 1.d0 - fofe
                  mofe = mo_m1
                  moat = mo
c					! get Z-composition of interpolated mix
                  if ( m .eq. 1 ) then
                     sum_niai = 0.0d0
                     do i = 1, nel_zmix
                        fninz_mix(i) = fofe * fninz_opalmixes(i,1)
     $                       + omfofe * fninz_opalmixes(i,khizat)
                        xiz_mix(i) = fninz_mix(i) * atwt_opalGS98(i)
                        sum_niai = sum_niai + xiz_mix(i)
                     enddo
                     do i = 1, nel_zmix
                        xiz_mix(i) = xiz_mix(i) / sum_niai
                        bracketife_mix(i) = log10(
     $                       ( max( fninz_mix(i) , 1.d-36 )
     $                       * fninz_opalmixes(kel_fe,1) )
     $                       / ( max( fninz_mix(kel_fe) , 1.d-36 )
     $                       * fninz_opalmixes(i,1) ) )
                     enddo
c-debug-chk[
c-debug-chk;                     if ( iout_debug_chk_ofe .gt. 0 ) then
c-debug-chk;                        write(6,2377) ofebrack,khizat,
c-debug-chk;     $                       bracketofe_opalmixes(khizat),
c-debug-chk;     $                       fofe,khizat,omfofe,khizat,khizat
c-debug-chk;                        do i = 1,nel_zmix
c-debug-chk;                           write(6,2376) cel_opalmixes(i),
c-debug-chk;     $                          xiz_opalmixes(i,1),
c-debug-chk;     $                          fninz_opalmixes(i,1),
c-debug-chk;     $                          xiz_opalmixes(i,khizat),
c-debug-chk;     $                          fninz_opalmixes(i,khizat),
c-debug-chk;     $                          xiz_mix(i),fninz_mix(i),
c-debug-chk;     $                          bracketife_mix(i)
c-debug-chk;                        enddo
c-debug-chk;                        write(6,'(" ")')
c-debug-chk;                        iout_debug_chk_ofe = iout_debug_chk_ofe - 1
c-debug-chk;                     endif
c-debug-chk]
                  endif
c
               else
c				! GS98 shift only
                  mofe = mo_m1
                  moat = mo_m1
                  fofe = 0.d0
                  omfofe = 1.d0
c
               endif
c		     ! compute [O/Fe],GS98 shifts relative to GN93hz opacities
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
     $                    badlogklim .and. abs( co(m,mc,mofe,k,il,kz) )
     $                    .lt. badlogklim ) then
                        dif = ( co(m,mc,mofe,k,il,kz) * omfofe
     $                       + co(m,mc,moat,k,il,kz) * fofe )
     $                       - co(m,mc,mo,k,il,kz)
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
                     else if ( il .gt. 1 ) then
                        co(m,mc,mo_m1,k,il,kz) = 
     $                       co(m,mc,mo_m1,k,il-1,kz)
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
c				! note: facxhz=0.0 for m=2, <1.0 for .04<Z<.06
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
                  else if ( il .gt. 1 ) then
                     co(m,mc,mo,k,il,kz) = co(m,mc,mo,k,il-1,kz)
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
c-debug-chk;      write(6,8418) (i,(n(i,j,kz),j=1,mo),' ',i=1,mx)
c-debug-chk; 8418 format(' '/' -- n(m,j): ',5(' (m=',i1,')',8i2,a1))
c-debug-chk]
c		! interpolate GN93hz opacity shifts for m=2, if possible; note
c		! that other shifts being interpolated among already contain
c		! the factor of facxhz. No need to revise any m=2 [O/Fe] shift.
c
      if ( khizat .gt. 0 .and. mx .ge. 4 .and.
     $     mxzero .eq. 1 .and. mx03 .eq. 2 ) then
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
               else if ( il .gt. 1 ) then
                  co(2,mc,mo,k,il,kz) = co(2,mc,mo,k,il-1,kz)
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
         do m = 1, mx
c
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
                        diffac = ((difl*0.16666667d0+0.5d0)*difl+1.d0)*
     &                       difl
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
                                 difl = log10( 10.d0**(co(m,1,1,k,il,kz)
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
c			! how many GN93hz X-values for the present value of Z
      mx_use = mx_hi
      do while ( xhi_in(mx_use-1) .gt. 0.999999d0 - zat )
         mx_use = mx_use - 1
      enddo
      nx_hi(kz) = mx_use
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
         do i = 1, kz
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
            cop_sto(1) = ' '
            cop_sto(2) = ' '
         else
            cop_sto(1) = copdir(:kope)
            cop_sto(2) = copdir(:kope)
         endif
c								! get filename
         if ( klozat .eq. 0 ) then
            cop_sto(1)(kope+1:) = cfile_opalmixes(1)
         else
            cop_sto(1)(kope+1:) = cfile_opalGS98(1)
         endif
         iu = iulow
c								! open file
         call open_chk_zip( iu, cop_sto(1), igzip,
     $        'STOP -- Error: hz-file (C+O=0.0,[O/Fe]=0) not found.' )
c
         line(1) = 0
         line(2) = 0
         itab_dum = 0
         khighx(kz) = 1
c							! Z > 0 & [O/Fe] > 0 ?
c
         if ( khighz_index .gt. 1 .and. kzalow + nzalmo .gt. 1 ) then
            khighx(kz) = 2
            if ( khighz .gt. 0 ) then
               cop_sto(2)(kope+1:) = cfile_opalmixes(khighz_index)
            else
               cop_sto(2)(kope+1:) = cfile_opalGS98(khighz_index)
            endif
            iu_ofe = iu + 1
            call open_chk_zip( iu_ofe, cop_sto(2), igzip_ofe,
     $           'STOP -- Error: hz-file (C+O=0.0,[O/Fe]>0) not found.'
     $           )
            itab_dum_ofe = 0
         endif
c
         ix = 0
         m = 0
         io = mo
         iz_hi = nzalmo
c						! loop over 'GN93hz' X-values:
c
         do while ( ix .lt. mx_use .and. iz_hi .lt. 5 )
c							! get position in co()
            ix = ix + 1
            m = m + 1
            if ( m .gt. mx ) then
               io = io - 1
               m = 1
            endif
c				! get Z and X values to look for in 'GN93hz'
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
                  else if ( zat .gt. 0.04d0 .and. zat .lt. 0.05d0 ) then
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
               ifound = mixfind(iu,1,igetzxi,i_rewind,itab_dum,
     $              line(1),zhi_look(kat),xhi_look(kat),0.0d0,0.0d0)
               if ( ifound .eq. 0 ) then
                  i_rewind = 1
                  igetzxi = 0
                  ifound = mixfind(iu,1,igetzxi,i_rewind,itab_dum,
     $                 line(1),zhi_look(kat),xhi_look(kat),0.0d0,0.0d0)
                  if ( ifound .eq. 0 ) then
                     write(6,1791) zhi_look(kat),xhi_look(kat),
     $                    0.0d0,0.0d0,cop_sto(1)(:lnblnk(cop_sto(1)))
                     stop ' STOP -- READCO: error reading hz-mix. '
                  endif
               endif
               if ( khighx(kz) .gt. 1 ) then
                  i_rewind = 0
                  igetzxi = 0
                  ifound = mixfind(iu_ofe,1,igetzxi,i_rewind,
     $                 itab_dum_ofe,line(2),
     $                 zhi_look(kat),xhi_look(kat),0.0d0,0.0d0)
                  if ( ifound .eq. 0 ) then
                     i_rewind = 1
                     igetzxi = 0
                     ifound = mixfind(iu_ofe,1,igetzxi,i_rewind,
     $                    itab_dum_ofe,line(2),
     $                    zhi_look(kat),xhi_look(kat),0.0d0,0.0d0)
                     if ( ifound .eq. 0 ) then
                        write(6,1791) zhi_look(kat),xhi_look(kat),
     $                       0.0d0,0.0d0,cop_sto(2)(:lnblnk(cop_sto(2)))
                        stop ' STOP -- READCO: error reading hz-mix. '
                     endif
                  endif
               endif
c				  ! loop over logT values, to read in opacities
               do k = 1, ntm
c					   ! read logT, & logKappa(R) for all R
                  line(1) = line(1) + 1
                  read(iu,8300) cin
c 8300             format(a137)
 8300             format(a277)
                  read(cin,8140) flt, (cofzhi(k,il,kat),il=1,nrm)
!! 8140             format(f4.2,19f7.3)
 8140             format(f4.2,39f7.3)
c								   ! bad logT ?
                  if ( abs(flogtin(k)-flt) .gt. 1.d-5 ) then
                     write(6,1734) flt, flogtin(k),
     $                    cop_sto(1)(:lnblnk(cop_sto(1))),
     $                    line(1), cin(:max(1,lnblnk(cin))),
     $                    zhi_look(kat),xhi_look(kat),0.0d0,0.0d0
                     stop ' STOP -- READCO: bad logT value. '
                  endif
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
                     line(2) = line(2) + 1
                     read(iu_ofe,8300) cin
                     read(cin,8140) flt, (coff(k,il),il=1,nrm)
                     if ( abs(flogtin(k)-flt) .gt. 1.d-5 ) then
                        write(6,1734) flt, flogtin(k),
     $                       cop_sto(2)(:lnblnk(cop_sto(2))),
     $                       line(2), cin(:max(1,lnblnk(cin))),
     $                       zhi_look(kat),xhi_look(kat),0.0d0,0.0d0
                        stop ' STOP -- READCO: bad logT value. '
                     endif
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
            xhi_use(ix,kz) = min( xhi_in(ix) , 1. - zat )
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
                     cof_tmp = qzinter(int_hi_z,2,0.05,3,
     $                    cofzhi(k,il,3),cofzhi(k,il,4),
     $                    cofzhi(k,il,5),cofzhi(k,il,6),
     $                    zhi_look(3),zhi_look(4),zhi_look(5),
     $                    zhi_look(6),zdel)
c							      ! Zsto(kz),X=0.95
                     cofzhi(k,il,1) = qzinter(int_hi_z,3,zat,2,
     $                    cofzhi(k,il,1),cofzhi(k,il,2),
     $                    cof_tmp,0.0d0,
     $                    zhi_look(1),zhi_look(2),0.05,0.0d0,zdel)
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
c			! (end of loop over 'GN93hz' X-values):
         enddo
c						! close 'GN93hz' file(s)
         if ( khighx(kz) .gt. 1 ) then
            call close_chk_zip( iu_ofe, cop_sto(2), igzip_ofe )
         endif
         call close_chk_zip( iu, cop_sto(1), igzip )
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
c			! get the dlogKappa values for 'GN93hz' X-values:
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
     $                 ' STOP -- Error: bad co(3:5,1,1,*,*) cannot be '
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
c
c  If required, read in CNO- and/or user-interpolation opacity tables
c
      if ( khighz_cno .eq. 0 ) then
c
         kavail_cno = 0
         kavail_user = 0
c
      else
c
         kdel = 1
         if ( khighz_cno .ge. 2 ) then
            khi = n_totmix
            if ( khighz_cno .eq. 2 ) then
               kdel = khi - n_cnobeg
               kavail_cno = 0
            endif
         else
            khi = n_totmix - 1
            kavail_user = 0
         endif
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
            cop_sto(1) = ' '
            cop_sto(2) = ' '
         else
            cop_sto(1) = copdir(:kope)
            cop_sto(2) = copdir(:kope)
         endif
c								! get filename
         if ( cfile_opalGS98(n_cnobeg) .ne. ' ' ) then
            cop_sto(1)(kope+1:) = cfile_opalGS98(n_cnobeg)
         else if ( khighz .gt. 0 ) then
            cop_sto(1)(kope+1:) = cfile_opalmixes(1)
         else
            cop_sto(1)(kope+1:) = cfile_opalGS98(1)
         endif
         last = lnblnk( cop_sto(1) )
         iu = iulow
c		    ! get indices where X=1-Z 'GN93hz' opacities will be stored
c
         call index_co_deltas( 5, mx_hi, kkx, kkc, kko )
c							 ! loop over CNO-files
         do kfil = n_cnobeg, khi, kdel
c
            iset = kfil - n_zmixes
            if ( iset .eq. 1 ) iset = 5
c
            if ( cfile_opalGS98(kfil) .ne. ' ' ) then
               cop_sto(2)(kope+1:) = cfile_opalGS98(kfil)
            else
               cop_sto(2) = cop_sto(1)(1:last) // cdef_CNO_ext(kfil)
            endif
c								! open file
            call open_chk_zip( iu, cop_sto(2), igzip,
     $           'STOP -- Error: hz-file (C+O=0.0,defCNO) not found.' )
c
            line(2) = 0
            itab_dum = 0
            ix = 0
            iz_hi = nzalmo
            igetzxi = 1
c								! X=1-Z indices
            call index_co_deltas( iset, mx_hi, jjx, jjc, jjo )
c
c						! loop over 'GN93hz' X-values:
c
            do while ( ix .lt. mx_use .and. iz_hi .lt. 5 )
c							  ! get indices in co()
               ix = ix + 1
               call index_co_deltas( iset, ix, jx, jc, jo )
c							    ! & stored 'GN93hz'
               call index_co_deltas( 5, ix, kx, kc, ko )
c
c				! get Z and X values to look for in CNO-file
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
     $                    1.d0 - zhi_look(iz+1) )
                  enddo
c			! check for X-column bifurcation at Z = 0.05, X = 0.95
c
                  if ( ix .eq. mx_hi - 1 .and. nzalmo .eq. 3 ) then
                     if ( zat .gt. 0.03d0 .and. zat .lt. 0.04d0 ) then
                        iz_hi = 2
                        int_hi_z = 0
                     else if (zat .gt. 0.04d0 .and. zat .lt. 0.05d0)then
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
                  ifound = mixfind(iu,-kfil,igetzxi,i_rewind,
     $                 itab_dum,line(2),
     $                 zhi_look(kat),xhi_look(kat),0.0d0,0.0d0)
                  igetzxi = 0
                  if ( ifound .eq. 0 ) then
                     i_rewind = 1
                     ifound = mixfind(iu,-kfil,igetzxi,i_rewind,
     $                    itab_dum,line(2),
     $                    zhi_look(kat),xhi_look(kat),0.0d0,0.0d0)
                     igetzxi = 0
                     if ( ifound .eq. 0 ) then
                        write(6,1791) zhi_look(kat),xhi_look(kat),
     $                       0.0d0,0.0d0,cop_sto(2)(:last)
                        stop ' STOP -- READCO: error reading CNO-mix. '
                     endif
                  endif
c				  ! loop over logT values, to read in opacities
                  do k = 1, ntm
c					   ! read logT, & logKappa(R) for all R
                     line(2) = line(2) + 1
c                     read(iu,'(a137)') cin
                     read(iu,'(a277)') cin
!!                     read(cin,'(f4.2,19f7.3)') flt,
                     read(cin,'(f4.2,39f7.3)') flt,
     $                    (cofzhi(k,il,kat),il=1,nrm)
c								   ! bad logT ?
                     if ( abs(flogtin(k)-flt) .gt. 1.d-5 ) then
                        write(6,1734) flt, flogtin(k),
     $                       cop_sto(2)(:lnblnk(cop_sto(2))),
     $                       line(2), cin(:max(1,lnblnk(cin))),
     $                       zhi_look(kat),xhi_look(kat),0.0d0,0.0d0
                        stop ' STOP -- READCO: bad logT value. '
                     endif
c
                     il_lo = 1
                     il_hi = nrm
c							      ! logKappa(R) is:
                     do il = nrm, 1, -1
c								       ! absent
                        if ( cin(7*il-2:7*il+4) .eq. '       ' ) then
                           il_hi = min( il_hi , il - 1 )
c									! 9.999
                        else if ( cofzhi(k,il,kat) .gt. 9.d0 ) then
                           il_lo = max( il_lo , il + 1 )
                        endif
                     enddo
c								  ! for deltas
                     if ( il_lo .gt. 1 .or. il_hi .lt. nrm ) then
                        do il = nrm, 1, -1
                           if ( il .gt. il_hi ) then
                              cofzhi(k,il,kat) = 2.d0 * cofzhi(k-1,il,
     &                             kat)
     $                             - cofzhi(k-2,il,kat)
                           else if ( il .lt. il_lo ) then
                              cofzhi(k,il,kat) = 2.d0 * cofzhi(k,il+1,
     &                             kat)
     $                             - cofzhi(k,il+2,kat)
                           endif
                        enddo
                     endif
c			   ! (end of loop to read opacities at all T-values):
                  enddo
c			! (end of loop over required Z-values for this X):
               enddo
c							      ! actual X at Zat
               xcno_use(ix,kz) = min( xhi_in(ix) , 1.d0 - zat )
c							     ! Z-interpolation:
               if ( iz_hi .le. 3 ) then
c					! standard case: for all T,R:
                  do k = 1, ntm
                     do il = 1, nrm
c							       ! logK at Zsto,X
                        coff(k,il) = qzinter(int_hi_z,1,zat,iz_hi,
     $                       cofzhi(k,il,1),cofzhi(k,il,2),
     $                       cofzhi(k,il,3),cofzhi(k,il,4),
     $                       zhi_look(1),zhi_look(2),zhi_look(3),
     $                       zhi_look(4),zdel)
                        int_hi_z = 1
                     enddo
                  enddo
c					    ! ELSE: bifurcation:
               else
c					        ! do both X = 1-Z and X = 0.95
                  xcno_use(mx_hi,kz) = 1.d0 - zat
c						! for all T,R:
                  do k = 1, ntm
                     do il = 1, nrm
c							   ! logK at Zsto,X=1-Z
                        coff(k,il) = qzinter(int_hi_z,1,zat,3,
     $                       cofzhi(k,il,3),cofzhi(k,il,4),
     $                       cofzhi(k,il,5),cofzhi(k,il,6),
     $                       zhi_look(3),zhi_look(4),zhi_look(5),
     $                       zhi_look(6),zdel)
c							  ! temp: Z=0.05,X=0.95
                        cof_tmp = qzinter(int_hi_z,2,0.05d0,3,
     $                       cofzhi(k,il,3),cofzhi(k,il,4),
     $                       cofzhi(k,il,5),cofzhi(k,il,6),
     $                       zhi_look(3),zhi_look(4),zhi_look(5),
     $                       zhi_look(6),zdel)
c							      ! Zsto(kz),X=0.95
                        cofzhi(k,il,1) = qzinter(int_hi_z,3,zat,2,
     $                       cofzhi(k,il,1),cofzhi(k,il,2),
     $                       cof_tmp,0.0d0,
     $                       zhi_look(1),zhi_look(2),0.05d0,0.0d0,zdel)
                        int_hi_z = 1
                     enddo
                  enddo
c					     ! smooth CNO-opac, if init_smo > 1
                  if ( init_smo .ge. 2 ) then
                     tmax = 10.d0
                     nset = ks81
                     RLS = alrf(1)
                     RLE = alrf(nrm)
                     nrhigh = int( dfsr(nr) * ( RLE - RLS ) + 1.00001d0)
                     nrlow = 1
                     call opaltab
                  endif
c						 ! store X=1-Z hz-opacity set
                  if ( kfil .eq. n_cnobeg ) then
                     do il = 1, nr
                        jl = il + nrdel
                        do k = 1, nt
                           co(jjx,jjc,jjo,k,il,kz) = coff(k+ntdel,jl)
                        enddo
                     enddo
                  else
                     do il = 1, nr
                        jl = il + nrdel
                        do k = 1, nt
                           co(jjx,jjc,jjo,k,il,kz) = coff(k+ntdel,jl)
     $                          - co(kkx,kkc,kko,k,il,kz)
                        enddo
                     enddo
                  endif
c				! get present X=0.95 hz-opacity set
                  do k = 1, ntm
                     do il = 1, nrm
                        coff(k,il) = cofzhi(k,il,1)
                     enddo
                  enddo
c
               endif
c					   ! smooth CNO-opac, if init_smo > 1
               if ( init_smo .ge. 2 ) then
                  tmax = 10.d0
                  nset = ks81
                  RLS = alrf(1)
                  RLE = alrf(nrm)
                  nrhigh = int( dfsr(nr) * ( RLE - RLS ) + 1.00001d0 )
                  nrlow = 1
                  call opaltab
               endif
c					      ! store present hz-opacity set
               if ( kfil .eq. n_cnobeg ) then
                  do il = 1, nr
                     jl = il + nrdel
                     do k = 1, nt
c-debug[
c-debug;                        if ( co(jx,jc,jo,k,il,kz) .lt. badlogklim )
c-debug;     $                       stop ' STOP -- Error: CNO overwrite. '
c-debug]
                        co(jx,jc,jo,k,il,kz) = coff(k+ntdel,jl)
                     enddo
                  enddo
               else
                  do il = 1, nr
                     jl = il + nrdel
                     do k = 1, nt
c-debug[
c-debug;                        if ( kfil .lt. n_totmix .and.
c-debug;     $                       co(jx,jc,jo,k,il,kz) .lt. badlogklim )
c-debug;     $                       stop ' STOP -- Error: CNO overwrite. '
c-debug]
                        co(jx,jc,jo,k,il,kz) = coff(k+ntdel,jl)
     $                       - co(kx,kc,ko,k,il,kz)
                     enddo
                  enddo
               endif
c			! (end of loop over 'GN93hz' X-values):
            enddo
c							! close CNO-file
            call close_chk_zip( iu, cop_sto(2), igzip )
c							! (end CNO-file loop):
         enddo
c
      endif
c			! restore old m value (this should not be necessary)
      m = mstore
c			! and return.
      return
      end
c
c******************************************************************************
c
      subroutine cointsmo(xxc,xxo,kz)
c     ===============================
c
c  The purpose of COINTSMO is to interpolate smoothly in C and O abundances.
c
c  This subroutine yields smoother opacities than alternate COINTERP below.
c
c  Note that the quadratic-interpolation function quad has been replaced with
c  the function qchk here; the latter function checks whether two of the
c  interpolation points are nearly coincident (which would magnify the effect
c  of any uncertainties in the tabulated opacities), and uses something more
c  nearly linear if so.  This is sometimes necessary to prevent wildly wrong
c  opacity values for certain Z-values above Z=0.03, and also in some cases to
c  allow linear interpolation along lines where only two opacity values exist.
c  For the special case where C or O is slightly negative (slight depletion in
c  C or O), the function qchk does a linear extrapolation using a combination
c  of the lowest three C or O gridpoints.
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
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
                     opl(it,ir,kz) = (1.-cfac) * qchk(is,1,cxx,
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
     $              / max( xod(1)-2.*xod(nc-1) , 1.e-6 ) )
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
c-debug; 9415             format(' COINTSMO(x',i1,',t',i2.2,',r',i2.2,
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
                     opl(it,ir,kz) = (1.-ofac) * qchk(is,1,oxx,
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
c
c******************************************************************************
c
      subroutine cointerp(xxc,xxo,kz)
c     ===============================
c
c     Interpolate in C and O abundances.
c
c  This subroutine yields less-smooth opacities than alternate
c  subroutine COINTSMO above.
c
c  Note that the quadratic-interpolation function quad has been
c  replaced with the function qchk here; the latter function checks
c  whether two of the interpolation points are nearly coincident
c  (which would magnify the effect of any uncertainties in the
c  tabulated opacities), and uses something more nearly linear
c  if so.  This is sometimes necessary to prevent wildly wrong
c  opacity values for certain Z-values above Z=0.03.
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
      common/bb_opal_z/ l1,l2,l3,l4,k1,k2,k3,k4,ip,iq(4),xodp,xcdp,
     $     xxco,cxx,oxx,kzf,kzg,kzh,kzf2
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
      ie = 100 * max( int(xxc) , 0 ) + 1
      i3 = max( min( indx(ie) + 1 , nc ) , 3 )
      i1 = i3-2
      i2 = i3-1
c
c.....find index of O grid:
c                   must also allow index = no, to avoid extrapolation
c
      ie = 100 * max( int(xxo) , 0 ) + 1
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
c
c******************************************************************************
c
      subroutine t6rinterp(slr,slt)
c     =============================
c
c     The purpose of this subroutine is to interpolate in logT6 and logR
c     NOTE THAT for 2-dimensional quadratic interpolation, IT DOES NOT MATTER 
c     which direction is interpolated first, horizontal or vertical: the
c     result is the same for interpolated value and derivatives.
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
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
      dix2 = max(min((alr(l3)-slr)*dfsr(l3),1.),0.)
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
c
c******************************************************************************
c
      subroutine qzlog4int( zlogd )
c     =============================
c
c..... this subroutine performs bi-quadratic interpolation of logKappa in the
c      log10(Z_i+zdel) values stored in the array zvint(nz), for each of the
c      relevant positions in the C,O-interpolated opacity matrix opl(nt,nr,nz)
c      (given the input values Z and zdel).  Note that this subroutine uses
c      the quadratic-interpolation function quad.  Depending on the number of
c      Z-values to interpolate among, single-quadratic or linear interpolation
c      may be used instead.  Note that zlogd = log10(Z+zdel).
c
c      NOTE that since errors in the opacities may be large compared to the
c      opacity differences between opacities at adjacent Z-values, quadratic
c      interpolation is forced to be monotonic by using values at adjacent
c      tabulated Z-values as upper and lower limits.  No such restriction is
c      placed on extrapolated values.
c
c PARAMETERS to specify opacity storage matrices (see OPALINIT):
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
c COMMON /a_opal_z/ : matrices for opacity storage (see OPALINIT):
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
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
c******************************************************************************
c
      function quad(ic,i,x,y1,y2,y3,x1,x2,x3)
c     =======================================
c
c..... this function performs a quadratic interpolation.
c
c  Storage for dx_i values that need not be computed on each call (see "ic");
c  NOTE that this same storage is used by all of QUAD, QDER, and QCHK.
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

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
c
c******************************************************************************
c
      double precision function qder(ic,i,x,y1,y2,y3,x1,x2,x3)
c     =======================================
c
c..... this function performs a quadratic interpolation; it is identical to the
c      function quad, except that it also computes the derivative dkap of the
c      quadratic at the given position x (see  common /d_opal_z/  below).
c
c  COMMON /d_opal_z/ : dkap returns the derivative (in interpolation-direction)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

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
c
c******************************************************************************
c
      double precision function qchk(ic,i,x,y1,y2,y3,x1,x2,x3)
c     =======================================
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
c      complete at R/r=24.  The ratios used in the code are actually r/(r+R).
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

c
      parameter (ratbeg=0.08d0,ratful=0.04d0,
     &     ratdel=1.d0/(ratbeg-ratful))
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
c
c******************************************************************************
c
      double precision function qzinter(ic,i,z,nmorez,f1,f2,f3,f4,z1,z2
     &     ,z3,z4,zdel)
c     ============================================================
c
c..... this function performs linear, quadratic, or bi-quadratic interpolation,
c      of logKappa in log(Z+zdel), for nmorez = 1, 2, or 3, respectively;
c      inputs are Z, nmorez = one less than the number of Z-values to
c      interpolate among, logKappa values f1 thru f4, Z-values z1 thru z4, and
c      zdel = 0.001 to make things work correctly near Z = 0.  Note that this
c      function is also sometimes used to interpolate in X or C or O.  It makes
c      use of the quadratic-interpolation function quad.
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

c  Storage for values that need not be computed on each call:
c
      common/qzint_opal_z/ v(15,5),f(15),omf(15)
      save /qzint_opal_z/
c===
      if ( ic .eq. 0 ) then
         if ( nmorez .gt. 0 ) then
            if ( zdel .lt. 1.d-5 .or. zdel .gt. 0.1011d0 ) stop
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
c-debug;     $        abs( v(i,5) - log10(z+zdel) ) ) .gt. 1.e-5 ) stop
c-debug;     $        ' STOP -- QZINTER: Error: expected same X-values. '
c-debug;         if ( nmorez .eq. 1 ) then
c-debug;            if ( abs( f(i) - (v(i,2)-v(i,5))/(v(i,2)-v(i,1)) ) .gt.
c-debug;     $           1.e-5 ) stop
c-debug;     $           ' STOP -- QZINTER: Error: expected same X-values. '
c-debug;         else
c-debug;            if ( abs( v(i,3) - log10(z3+zdel) ) .gt. 1.e-5 ) stop
c-debug;     $           ' STOP -- QZINTER: Error: expected same X-values. '
c-debug;            if ( nmorez .ge. 3 ) then
c-debug;               if ( max( abs( v(i,4) - log10(z4+zdel) ) ,
c-debug;     $              abs( f(i) - (v(i,3)-v(i,5))/(v(i,3)-v(i,2)) ) )
c-debug;     $              .gt. 1.e-5 ) stop
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
c
c******************************************************************************
c
      function mixfind(iu,iofe,igetzxi,irew,itab,l,zget,xget,cget,oget)
c     =================================================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      parameter ( ks81=ntm-3, ks83=ks81+1, ks60=ks81-21, ks61=ks60+1,
     $     alrlo=-9.0d0, flogtlo=3.75d0, flogt60=6.6d0, flogt81=9.3d0 )
!!     $     alrlo=-8.0, flogtlo=3.75, flogt60=6.0, flogt81=8.1 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
!!      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
!!      parameter ( nel_zmix=19, n_zmixes=7, kel_o=3, kel_fe=nel_zmix-1 )
      parameter ( nel_zmix=19, n_zmixes=8, kel_o=3, kel_fe=nel_zmix-1 )
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
      parameter ( n_totmix = n_zmixes + 5, n_cnobeg = n_zmixes + 1 )
c
      character*255 cfile_opalGS98(n_totmix)
      common /opalGS98mixes/ bracketofe_opalGS98(n_totmix),
     $     xofe_opalGS98(n_totmix),xiz_opalGS98(nel_zmix,n_totmix),
     $     fninz_opalGS98(nel_zmix,n_totmix),atwt_opalGS98(nel_zmix),
     $     cfile_opalGS98
      save /opalGS98mixes/
c___
      character*255 cin
      character*1 ch(255)
      equivalence (ch(1),cin)
c-debug-chk[
c-debug-chk;      common /mixfind_opal_debug_chk/ iout_debug_chk(-10:10)
c-debug-chk;c---
c-debug-chk;      data iout_debug_chk / 21 * 99999 /
c-debug-chk]
c===
      ifound = 0
      cin = '                '
c					! if must rewind to beginning of file:
      if ( irew .ne. 0 ) then
c
         rewind(iu)
         l = 0
c					! else, if must get Xi/Z values:
      else if ( igetzxi .ne. 0 ) then
c
         do while( cin(1:31) .ne. ' Element   Abundance - relative' )
            l = l + 1
            read(iu,'(a255)',end=900) cin
         enddo
         sum_X = 0.d0
         sum_XoverA = 0.d0
         kel = 1
c					! begin loop to get Xi/Z values
         do while( kel .le. nel_zmix )
            l = l + 1
            read(iu,'(a255)',end=900) cin
            if ( cin(1:16) .eq. ' Table Summaries' ) then 
               print *,'Table summaries',kel,nel_zmix
               goto 50
            endif
            ke = 255
            do while( ke .gt. 1 .and. ch(ke) .eq. ' ' )
               ke = ke-1
            enddo
            kb = 1
            do while( kb .le. ke .and. ch(kb) .eq. ' ' )
               kb = kb+1
            enddo
            if ( cin(kb:kb+1) .eq. cel_opalmixes(kel) ) then
c								! get Ni
               if ( ke .lt. kb+20 ) then 
                  print *,'ke < kb+20',ke,kb+20,cin(kb:kb+1)
     $                 ,cel_opalmixes(kel)
                  goto 50
               endif
               read(cin(ke-8:ke),'(f9.6)',err=50) vx
c							! last may be atomic wt
               if ( vx .gt. 1.000001d0 ) then
                  ke = ke-9
                  do while( ke .gt. 1 .and. ch(ke) .eq. ' ' )
                     ke = ke-1
                  enddo
                  if ( ke .lt. kb+20 ) goto 50
                  read(cin(ke-8:ke),'(f9.6)',err=50) vx
               endif
               sum_X = sum_X + vx
               sum_XoverA = sum_XoverA + vx / atwt_opalGS98(kel)
c								   ! get Xi
               ke = ke-9
               do while( ke .gt. 1 .and. ch(ke) .eq. ' ' )
                  ke = ke-1
               enddo
               if ( ke .lt. kb+11 ) goto 50
               read(cin(ke-7:ke),'(f8.6)',err=50) vn
c
               if ( iofe .lt. 0 ) then
                  if ( igetzxi .gt. 0 ) then
c-debug-chk[
c-debug-chk;                     if ( iout_debug_chk(iofe) .gt. 0 ) then
c-debug-chk;                        if ( kel .eq. 1 ) write(6,3813) iofe,
c-debug-chk;     $                       cfile_opalGS98(-iofe)
c-debug-chk; 3813                   format(' '/' Mix in GS98hz and in file',i3,
c-debug-chk;     $                       ' = ',a40/
c-debug-chk;     $                       '  i  el      GS98:Ni   GS98:Xi',
c-debug-chk;     $                       ' stored:Ni stored:Xi  read:Ni ',
c-debug-chk;     $                       '  read:Xi sto[i/Fe]'/
c-debug-chk;     $                       ' --  --    --------- ---------',
c-debug-chk;     $                       ' --------- --------- ---------',
c-debug-chk;     $                       ' --------- --------')
c-debug-chk;                        write(6,3814) kel, cel_opalmixes(kel),
c-debug-chk;     $                       fninz_opalGS98(kel,1),
c-debug-chk;     $                       xiz_opalGS98(kel,1),
c-debug-chk;     $                       fninz_opalGS98(kel,-iofe),
c-debug-chk;     $                       xiz_opalGS98(kel,-iofe), vn, vx,
c-debug-chk;     $                       log10( (
c-debug-chk;     $                       max(fninz_opalGS98(kel,-iofe),1.e-36)
c-debug-chk;     $                       * fninz_opalGS98(kel_fe,1) )
c-debug-chk;     $                       / ( max(fninz_opalGS98(kel_fe,-iofe),
c-debug-chk;     $                       1.e-36) * fninz_opalGS98(kel,1) ) )
c-debug-chk; 3814                   format(i3,': ',a2,' --',6f10.7,f9.5)
c-debug-chk;                     endif
c-debug-chk]
                     xiz_opalGS98(kel,-iofe) = vx
                     fninz_opalGS98(kel,-iofe) = vn
                  else if ( max( abs( fninz_opalGS98(kel,-iofe) - vn )
     $                    , abs( xiz_opalGS98(kel,-iofe) - vx ) )
     $                    .gt. 1.5d-5 ) then
                     write(6,35) cel_opalmixes(kel),
     $                    cfile_opalGS98(-iofe)(1:
     $                    lnblnk(cfile_opalGS98(-iofe))),
     $                    vx,xiz_opalGS98(kel,-iofe),
     $                    vn,fninz_opalGS98(kel,-iofe)
 35                  format(' '/' READCO: Warning: stored value',
     $                    ' differs from new Xi/Z for ',a2,
     $                    ' in GS98 file:'/' ',a/
     $                    '          new Xi/Z',f9.6,' vs.',f9.6,
     $                    ' ,  new Ni/Nz',f9.6,' vs.',f9.6/' ')
                  endif
               else if ( igetzxi .gt. 0 ) then
                  xiz_opalmixes(kel,iofe) = vx
                  fninz_opalmixes(kel,iofe) = vn
c-debug-chk[
c-debug-chk;                  if ( iout_debug_chk(iofe) .gt. 0 ) then
c-debug-chk;                     if ( kel .eq. 1 ) write(6,3813) iofe,
c-debug-chk;     $                    cfile_opalmixes(iofe)
c-debug-chk;                     write(6,3814) kel, cel_opalmixes(kel),
c-debug-chk;     $                    fninz_opalmixes(kel,1),
c-debug-chk;     $                    xiz_opalmixes(kel,1),
c-debug-chk;     $                    fninz_opalmixes(kel,iofe),
c-debug-chk;     $                    xiz_opalmixes(kel,iofe), vn, vx,
c-debug-chk;     $                    log10(
c-debug-chk;     $                    ( max(fninz_opalmixes(kel,iofe),1.e-36)
c-debug-chk;     $                    * fninz_opalmixes(kel_fe,1) )
c-debug-chk;     $                    / ( max(fninz_opalmixes(kel_fe,iofe),
c-debug-chk;     $                    1.e-36) * fninz_opalmixes(kel,1) ) )
c-debug-chk;                  endif
c-debug-chk]
               else if ( max( abs( fninz_opalmixes(kel,iofe) - vn ) ,
     $                 abs( xiz_opalmixes(kel,iofe) - vx ) )
     $                 .gt. 1.5d-5 ) then
                  write(6,40) cel_opalmixes(kel),cfile_opalmixes(iofe),
     $                 vx,xiz_opalmixes(kel,iofe),
     $                 vn,fninz_opalmixes(kel,iofe)
 40               format(' '/' READCO: Warning: new Xi/Z for ',a2,
     $                 ' in ',a8,' mix differs from stored value:'/
     $                 '          new Xi/Z',f9.6,' vs.',f9.6,
     $                 ' ,  new Ni/Nz',f9.6,' vs.',f9.6/' ')
c
c-dont;                  goto 60
               endif
               kel = kel+1
            endif
c					! end of loop to get Xi/Z values
         enddo
c					! check Xi vs Ni; get xO/xFe and [O/Fe]
         if ( igetzxi .gt. 0 ) then
            kel_err = 0
            if ( iofe .lt. 0 ) then
               do kel = 1, nel_zmix
                  if ( abs( xiz_opalGS98(kel,-iofe)
     $                 / ( atwt_opalGS98(kel) * sum_XoverA )
     $                 - fninz_opalGS98(kel,-iofe)) .gt. 0.00001d0) then
                     kel_err = kel_err + 1
                     write(6,4613) kel, cel_opalmixes(kel), iofe,
     $                    fninz_opalGS98(kel,-iofe),
     $                    xiz_opalGS98(kel,-iofe)
     $                    / ( atwt_opalGS98(kel) * sum_XoverA )
 4613                format(' READCO: Error in element',i3,
     $                    ' = "',a2,'" when reading in mix',i3,':'/
     $                    '   number fraction',f10.6,' does not match',
     $                    f10.6,' = (Xi/Ai) / [Sum{Xi/Ai}]')
                  endif
                  xiz_opalGS98(kel,-iofe) =
     $                 xiz_opalGS98(kel,-iofe) / sum_X
                  fninz_opalGS98(kel,-iofe) = xiz_opalGS98(kel,-iofe)
     $                 / ( atwt_opalGS98(kel) * sum_XoverA )
               enddo
c-debug-chk[
c-debug-chk;               if ( iout_debug_chk(iofe) .gt. 0 ) then
c-debug-chk;                  write(6,3815) xofe_opalGS98(1),
c-debug-chk;     $                 bracketofe_opalGS98(1),
c-debug-chk;     $                 xofe_opalGS98(-iofe),
c-debug-chk;     $                 bracketofe_opalGS98(-iofe),
c-debug-chk;     $                 fninz_opalGS98(kel_o,-iofe)
c-debug-chk;     $                 / max(fninz_opalGS98(kel_fe,-iofe),1.e-36),
c-debug-chk;     $                 log10( ( fninz_opalGS98(kel_o,-iofe)
c-debug-chk;     $                 / max(fninz_opalGS98(kel_fe,-iofe),1.e-36) )
c-debug-chk;     $                 / xofe_opalGS98(1) )
c-debug-chk; 3815             format('O/Fe,[O/Fe]',6f10.6/' ')
c-debug-chk;                  iout_debug_chk(iofe) = iout_debug_chk(iofe) - 1
c-debug-chk;               endif
c-debug-chk]
               xofe_opalGS98(-iofe) = fninz_opalGS98(kel_o,-iofe)
     $              / max( fninz_opalGS98(kel_fe,-iofe) , 1.d-36 )
               bracketofe_opalGS98(-iofe) = log10( max( 1.d-36 ,
     $              xofe_opalGS98(-iofe) / xofe_opalGS98(1) ) )
            else
               do kel = 1, nel_zmix
                  if ( abs( xiz_opalmixes(kel,iofe)
     $                 / ( atwt_opalGS98(kel) * sum_XoverA )
     $                 - fninz_opalmixes(kel,iofe)) .gt. 0.00001d0) then
                     kel_err = kel_err + 1
                     write(6,4613) kel, cel_opalmixes(kel), iofe,
     $                    fninz_opalmixes(kel,iofe),
     $                    xiz_opalmixes(kel,iofe)
     $                    / ( atwt_opalGS98(kel) * sum_XoverA )
                  endif
                  xiz_opalmixes(kel,iofe) =
     $                 xiz_opalmixes(kel,iofe) / sum_X
                  fninz_opalmixes(kel,iofe) = xiz_opalmixes(kel,iofe)
     $                 / ( atwt_opalGS98(kel) * sum_XoverA )
               enddo
c-debug-chk[
c-debug-chk;               if ( iout_debug_chk(iofe) .gt. 0 ) then
c-debug-chk;                  write(6,3815) xofe_opalmixes(1),
c-debug-chk;     $                 bracketofe_opalmixes(1),
c-debug-chk;     $                 xofe_opalmixes(iofe),
c-debug-chk;     $                 bracketofe_opalmixes(iofe),
c-debug-chk;     $                 fninz_opalmixes(kel_o,iofe)
c-debug-chk;     $                 / max(fninz_opalmixes(kel_fe,iofe),1.e-36),
c-debug-chk;     $                 log10( ( fninz_opalmixes(kel_o,iofe)
c-debug-chk;     $                 / max(fninz_opalmixes(kel_fe,iofe),1.e-36) )
c-debug-chk;     $                 / xofe_opalmixes(1) )
c-debug-chk;                  iout_debug_chk(iofe) = iout_debug_chk(iofe) - 1
c-debug-chk;               endif
c-debug-chk]
               xofe_opalmixes(iofe) = fninz_opalmixes(kel_o,iofe)
     $              / max( fninz_opalmixes(kel_fe,iofe) , 1.d-36 )
               bracketofe_opalmixes(iofe) = log10( max( 1.d-36 ,
     $              xofe_opalmixes(iofe) / xofe_opalmixes(1) ) )
            endif
            if ( kel_err .ne. 0 ) stop
     $           ' STOP -- READCO: Incompatible Ni vs. Xi read in. '
         endif
c					! no read error: jump to continuation
         goto 60
c				! if error reading Xi/Z values, say so
 50      write(6,20) kel
 20      format(' '/' READCO: Warning: error reading',
     $        ' Z-abundance fractions at element',i3/' ')
         if ( igetzxi .ge. 9 ) stop
     $        ' STOP -- READCO: Cannot get mix from user [O/Fe]-file. '
c
c								! continuation
 60      continue
c				! end of reading Xi/Z values in file header
      endif
c							! find start of tables
      if ( irew .ne. 0 .or. igetzxi .ne. 0 ) then
         do while( cin(1:30) .ne. '******************************' )
            l = l + 1
            read(iu,'(a255)',end=900) cin
         enddo
         igetzxi = 0
      endif
c				  ! look for mix with required composition:
      do while ( ifound .eq. 0 )
         l = l + 1
         read(iu,'(a255)',end=900) cin
         if ( cin(1:7) .eq. 'TABLE #' ) then
            ke = 90
            do while( ke .gt. 1 .and. ch(ke) .eq. ' ' )
               ke = ke-1
            enddo
            if ( ke .lt. 60 ) goto 900
            read(cin(ke-48:ke),100) xat,yat,zat,cat,oat
 100        format(3(3x,f6.4),2(5x,f6.4))
            if ( max(abs(zat-zget),abs(xat-xget),
     $           abs(cat-cget),abs(oat-oget)) .lt. 1.d-6 ) ifound = -1
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
         l = l + 1
         read(iu,'(a255)',end=900) cin
      enddo
      l = l + 1
      read(iu,110,err=900,end=900) cin(1:4),(alrf(i),i=1,nrm)
!! 110  format(a4,f6.1,18f7.1)
 110  format(a4,f6.1,38f7.1)
c					! this may or may not be useful/correct
      if ( cin(1:4) .ne. 'logT' .or.
     $     abs(alrf(1)-alrlo) .gt. 1.d-5 ) goto 900
c
      do k = 2,nrm
         if ( abs(alrf(k)-alrf(k-1)-0.5) .gt. 1.e-5 ) stop
     $        ' STOP -- READCO: bad  log R  value in table read in. '
      enddo
c				     ! read blank line before first table line
      l = l + 1
      read(iu,'(a255)',end=900) cin
c				     ! table header lines appear correct:
      ifound = iabs( ifound )
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
c******************************************************************************
c
      subroutine index_co_deltas( iset, kxhz, jx, jc, jo )
c     ====================================================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      parameter ( mx_hi=2*mx, mx_m1=mx-1, mo_m1=mo-1, mo_m2=mo-2,
     $     mo_m10=mo-mx_hi, mo_p1=mo+1, mc_m1=mc-1, mc_m2=mc-2,
     $     mc_p7=mc+mx+2 )
c
      if ( iset .le. 0 .or. iset .gt. 5 .or. kxhz .le. 0 .or.
     $     kxhz .gt. mx_hi .or. mx .ne. 5 ) stop
     $     ' INDEX_DELTAS: Error: bad inputs: cannot happen. '
c
      if ( iset .eq. 2 ) then
         if ( kxhz .le. 5 ) then
            jx = kxhz
            jc = mc
            jo = mo_m2
         else
            jx = kxhz - 5
            jc = mc_m1
            jo = mo_m1
         endif
      else if ( iset .eq. 3 ) then
         jx = mx
         if ( kxhz .le. 5 ) then
            jc = mc
            jo = kxhz
         else
            jc = mc_m1
            jo = kxhz - 5
         endif
      else if ( iset .eq. 4 ) then
         jx = mx
         if ( kxhz .le. 5 ) then
            jc = kxhz
            jo = mo_m1
         else
            jc = kxhz - 5
            jo = mo_m2
         endif
      else if ( iset .eq. 5 ) then
         if ( kxhz .le. 5 ) then
            jx = mx
            jc = mc_m2
            jo = mo_p1 - kxhz
         else
            jx = mx_m1
            jc = min( mc_p7 - kxhz , mc )
            jo = min( mo_m10 + kxhz , mo_m1 )
         endif
      else
         jc = mc
         if ( kxhz .le. 5 ) then
            jx = kxhz
            jo = mo
         else
            jx = kxhz - 5
            jo = mo_m1
         endif
      endif
c
      return
      end
c
c******************************************************************************
c
      subroutine finish_cno
c     =====================
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
     $     nr=nre+1-nrb, ntm=76, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
c
      parameter ( mx_hi=2*mx, mo_m1=mo-1, mo_m2=mo-2 )
c
      common/a_opal_z/ indx(101),t6list(nt),alr(nr),n(mx,mo,nz),
     $     alt(nt),dfs(nt),dfsr(nr),b(3),m,mf,xa(8),alrf(nrm),
     $     flogtin(ntm),dfsx(mx),oxf(mx,mc,nz),cxf(mx,mo,nz),
     $     xcdf(mx,mo,nz),xodf(mx,mc,nz),itime,cxdf(mx,mo,nz),
     $     oxdf(mx,mc,nz),q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),
     $     xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz(nz),xxh,
     $     xx(mx),nc,no,zsto(nz),zvint(nz),dfsz(nz),zacc(nz),
     $     zlow,zmiddle,zhigh,zlo_ex,zhi_ex,numz,
     $     co(mx,mc,mo,nt,nr,nz),opk(mx,4),opl(nt,nr,nz),cof(nt,nr)
      save /a_opal_z/
c
      common /xhi_opal_z/ xhi_in(mx_hi), xcno_use(mx_hi,nz),
     $     xhi_use(mx_hi,nz), xxx_cno(mx_hi), xxx_hi(mx_hi),
     $     nx_hi(nz), ireq_hi(mx_hi), khighx(nz), kavail_xhi, kuse_xhi,
     $     kdo_xhi, kavail_cno, kuse_cno, kdo_cno, kavail_user,
     $     kuse_user, kdo_user
      save /xhi_opal_z/
c
!!      parameter ( nel_zmix=19, n_zmixes=5, kel_o=3, kel_fe=nel_zmix-1 )
!!      parameter ( nel_zmix=19, n_zmixes=7, kel_o=3, kel_fe=nel_zmix-1 )
      parameter ( nel_zmix=19, n_zmixes=8, kel_o=3, kel_fe=nel_zmix-1 )
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
      parameter ( n_totmix = n_zmixes + 5, n_cnobeg = n_zmixes + 1 )
c
      character*255 cfile_opalGS98(n_totmix)
      common /opalGS98mixes/ bracketofe_opalGS98(n_totmix),
     $     xofe_opalGS98(n_totmix),xiz_opalGS98(nel_zmix,n_totmix),
     $     fninz_opalGS98(nel_zmix,n_totmix),atwt_opalGS98(nel_zmix),
     $     cfile_opalGS98
      save /opalGS98mixes/
c
      common /cno_delta_opal_z/ fcno_mul(4), fninz_cno(nel_zmix,5),
     $     xiz_cno(nel_zmix,5), d_fninz_user(nel_zmix),
     $     fcno_fac(0:3,4), fninz_heavy, xiz_heavy, d_fninz_u_heavy,
     $     s_ninzai_mix, ds_ninzai_u, fn_o_over_cno, fninz_co_mix
      save /cno_delta_opal_z/
c
      common /c_level_err_opal_z/ level_err
      save /c_level_err_opal_z/
c___
      parameter ( ncno2 = n_cnobeg + 1, ncno3 = ncno2 + 1,
     $     ncno4 = ncno3 + 1 )
c
      dimension cno_del_chk(2:4,4)
c===
c
c  If user-specified opacity shifts were read in and the error-check level = 2,
c  then check the composition difference relative to the standard opacity file:
c
      if ( kavail_user .gt. 0 ) then
c
         del_max = 0.0d0
         do i = 1, nel_zmix
            fninz_cno(i,5) = fninz_opalGS98(i,n_totmix)
            xiz_cno(i,5) = xiz_opalGS98(i,n_totmix)
            d_fninz_user(i) = fninz_opalGS98(i,n_totmix)
     $           - fninz_opalGS98(i,n_cnobeg)
            del_max = max( del_max , abs( d_fninz_user(i) ) )
         enddo
c
         if ( del_max .lt. 1.5d-5 .and. level_err .gt. 0 ) then
            write(6,10)
 10         format(' '/' WARNING: user-specified OPAL opacity',
     $           ' interpolation file has a'/
     $           '          composition identical to that',
     $           ' of the standard opacity file.'/' ')
         else if ( del_max .lt. 1.d-3 .and. level_err .ge. 2 ) then
            write(6,20) del_max
 20         format(' '/' WARNING: user-specified OPAL opacity',
     $           ' interpolation file has a composition'/
     $           '          similar to that of the standard',
     $           ' opacity file: max delta =',1p,e9.2/' ')
         endif
c
      endif
c
c  If CNO-interpolation opacity shifts were read in, check them and compute
c  some useful multiplicative factors:
c
      if ( kavail_cno .gt. 0 ) then
c				     ! check for similar C,N,O,Ne compositions:
         del_12 = 0.0d0
         del_13 = 0.0d0
         del_14 = 0.0d0
         del_23 = 0.0d0
         del_34 = 0.0d0
         del_24 = 0.0d0
c
         do i = 1, 4
c							! for linear-dep check
            cno_del_chk(2,i) = fninz_opalGS98(i,ncno2)
     $           - fninz_opalGS98(i,n_cnobeg)
            cno_del_chk(3,i) = fninz_opalGS98(i,ncno3)
     $           - fninz_opalGS98(i,n_cnobeg)
            cno_del_chk(4,i) = fninz_opalGS98(i,ncno4)
     $           - fninz_opalGS98(i,n_cnobeg)
c
c				! check for very similar pairs of compositions
c
            del_12 = max( del_12 , abs( cno_del_chk(2,i) ) )
            del_13 = max( del_13 , abs( cno_del_chk(3,i) ) )
            del_14 = max( del_14 , abs( cno_del_chk(4,i) ) )
            del_23 = max( del_23 , abs( fninz_opalGS98(i,ncno3)
     $           - fninz_opalGS98(i,ncno2) ) )
            del_34 = max( del_34 , abs( fninz_opalGS98(i,ncno4)
     $           - fninz_opalGS98(i,ncno3) ) )
            del_24 = max( del_24 , abs( fninz_opalGS98(i,ncno2)
     $           - fninz_opalGS98(i,ncno4) ) )
c						! are all of C,N,O,Ne varied?
            if ( max( abs( cno_del_chk(2,i) ) ,
     $           abs( cno_del_chk(3,i) ) ,
     $           abs( cno_del_chk(4,i) ) ) .lt. 0.05d0 ) kavail_cno = 0
c
         enddo
c
         if ( min( del_12 , del_13 , del_14 ,
     $        del_23 , del_34 , del_24 ) .lt. 0.05d0 ) kavail_cno = 0
c
c				! check for linear dependence among CNO-mixes
         lindep = 0
c
         if ( kavail_cno .gt. 0 ) then
c
            do k = 3, 4
               i = 1
               do while ( i .lt. kel_o .and.
     $              abs( cno_del_chk(2,i) ) .lt. 0.009d0 )
                  i = i + 1
               enddo
               f = cno_del_chk(k,i) / cno_del_chk(2,i)
               if ( max(
     $              abs( f * cno_del_chk(2,1) - cno_del_chk(k,1) ) ,
     $              abs( f * cno_del_chk(2,2) - cno_del_chk(k,2) ) ,
     $              abs( f * cno_del_chk(2,kel_o)
     $              - cno_del_chk(k,kel_o) ) ) .lt. 0.005d0 ) lindep = 1
            enddo
c
            if ( lindep .eq. 0 ) then
c
               g = cno_del_chk(2,1) * cno_del_chk(3,2)
     $              - cno_del_chk(3,1) * cno_del_chk(2,2)
               if ( abs( g ) .lt. 2.5d-5 ) then
                  lindep = 1
               else
                  f = ( cno_del_chk(3,2) * cno_del_chk(4,1)
     $                 - cno_del_chk(4,2) * cno_del_chk(3,1) ) / g
                  g = ( cno_del_chk(2,1) * cno_del_chk(4,2)
     $                 - cno_del_chk(4,1) * cno_del_chk(2,2) ) / g
                  if ( abs( f * cno_del_chk(2,3) + g * cno_del_chk(3,3)
     $                 - cno_del_chk(4,3) ) .lt. 0.005d0 ) lindep = 1
               endif
c
            endif
c
         endif
c								    ! Bad CNO?
         if ( kavail_cno .eq. 0 ) then
c
            if ( level_err .gt. 0 ) write(6,30)
 30         format(' WARNING: CNO-interpolation in OPAL',
     $           ' opacities is NOT POSSIBLE: the C,N,O,Ne'/
     $           ' abundances are too similar in the',
     $           ' specified CNO-interpolation opacity files.')
            if ( level_err .ge. 2 ) stop
     $           ' STOP -- READCO Error: bad CNO-interpolation files. '
c
c							     ! or low main CNO?
         else if ( fninz_opalGS98(1,n_cnobeg) .lt. 0.05d0 .or.
     $           fninz_opalGS98(2,n_cnobeg) .lt. 0.01d0 .or.
     $           fninz_opalGS98(kel_o,n_cnobeg) .lt. 0.2d0 ) then
c
            kavail_cno = 0
            if ( level_err .gt. 0 ) write(6,40)
 40         format(' WARNING: CNO-interpolation in OPAL',
     $           ' opacities is NOT POSSIBLE: the "standard"'/
     $           ' CNO-opacity-file has low abundance(s) of C,N,O =',
     $           1p,3e9.2)
            if ( level_err .ge. 2 ) stop
     $           ' STOP -- READCO Error: bad "standard"-CNO-file. '
c
c							      ! enhanced C,O or
c								 ! depleted Ne?
         else if ( fninz_opalGS98(1,n_cnobeg) .lt. 0.999d0 * max(
     $           fninz_opalGS98(1,ncno2) , fninz_opalGS98(1,ncno3) ,
     $           fninz_opalGS98(1,ncno4) ) .or.
     $           fninz_opalGS98(kel_o,n_cnobeg) .lt. 0.999d0 * max(
     $           fninz_opalGS98(kel_o,ncno2) ,
     $           fninz_opalGS98(kel_o,ncno3) ,
     $           fninz_opalGS98(kel_o,ncno4) ) .or.
     $           fninz_opalGS98(4,n_cnobeg) .gt. 1.001d0 * min(
     $           fninz_opalGS98(4,ncno2) , fninz_opalGS98(4,ncno3) ,
     $           fninz_opalGS98(4,ncno4) ) ) then
c
            kavail_cno = 0
            if ( level_err .gt. 0 ) write(6,50)
 50         format(' WARNING: CNO-interpolation in OPAL',
     $           ' opacities is NOT POSSIBLE: these CNO-'/
     $           ' interpolation files should NOT have',
     $           ' C or O enhancements or Ne depletions.')
            if ( level_err .ge. 2 ) stop
     $           ' STOP -- READCO Error: bad CNO-interpolation files. '
c
         else if ( lindep .gt. 0 ) then
c
            kavail_cno = 0
            if ( level_err .gt. 0 ) write(6,60)
 60         format(' WARNING: CNO-interpolation in OPAL',
     $           ' opacities is NOT POSSIBLE: compositions of'/
     $           '            the CNO-interpolation files',
     $           ' are linearly dependent in {C,N,O}-space.')
            if ( level_err .ge. 2 ) stop
     $           ' STOP -- READCO Error: bad CNO-interpolation files. '
c
         else
c			   ! Else: O.K. so far; check elements heavier than Ne:
            del_12 = 0.0d0
            del_13 = 0.0d0
            del_14 = 0.0d0
            del_23 = 0.0d0
            del_34 = 0.0d0
            del_24 = 0.0d0
            del_max = 0.0d0
c
            do i = 5, nel_zmix
c								! sums of diffs
               del_12 = del_12 + fninz_opalGS98(i,n_cnobeg)
     $              - fninz_opalGS98(i,ncno2)
               del_13 = del_13 + fninz_opalGS98(i,n_cnobeg)
     $              - fninz_opalGS98(i,ncno3)
               del_14 = del_14 + fninz_opalGS98(i,n_cnobeg)
     $              - fninz_opalGS98(i,ncno4)
               del_23 = del_23 + fninz_opalGS98(i,ncno3)
     $              - fninz_opalGS98(i,ncno2)
               del_34 = del_34 + fninz_opalGS98(i,ncno4)
     $              - fninz_opalGS98(i,ncno3)
               del_24 = del_24 + fninz_opalGS98(i,ncno2)
     $              - fninz_opalGS98(i,ncno4)
c								! max diff
               del_max = max( abs( fninz_opalGS98(i,n_cnobeg)
     $              - fninz_opalGS98(i,ncno2) ) ,
     $              abs( fninz_opalGS98(i,n_cnobeg)
     $              - fninz_opalGS98(i,ncno3) ) ,
     $              abs( fninz_opalGS98(i,n_cnobeg)
     $              - fninz_opalGS98(i,ncno4) ) , del_max )
c
            enddo
c
            del_sum = max( abs(del_12) , abs(del_13) , abs(del_14) ,
     $           abs(del_23) , abs(del_34) , abs(del_24) )
c								 ! Bad heavies?
            if ( del_sum .gt. 0.001d0 .or. del_max .gt. 0.05d0 ) then
c
               kavail_cno = 0
               if ( level_err .gt. 0 ) write(6,70) del_sum, del_max
 70            format(' WARNING: CNO-interpolation in OPAL',
     $              ' opacities is NOT POSSIBLE: C+N+O+Ne sums'/
     $              ' differ by',1p,e9.2,
     $              ' > 0.001, OR max heavy-element-delta of',e9.2,
     $              ' > 0.05')
               if ( level_err .ge. 2 ) stop
     $              ' STOP -- READCO Error: bad CNO-interp-files. '
c								   ! Else: O.K.
            else
c			  ! composition of main-CNO file same as OPAL mix used?
               if ( max(
     $              abs( fninz_opalGS98(1,n_cnobeg) - fninz_mix(1) ) ,
     $              abs( fninz_opalGS98(2,n_cnobeg) - fninz_mix(2) ) ,
     $              abs( fninz_opalGS98(kel_o,n_cnobeg)
     $              - fninz_mix(kel_o) ) ,
     $              abs( fninz_opalGS98(4,n_cnobeg) - fninz_mix(4) ) )
     $              .lt. 0.001d0 ) then
c					! no CNO-modifications necessary
                  do k = 1, 4
c
                     fcno_mul(k) = 1.0d0
c
                     j = k + n_zmixes
c
                     do i = 1, nel_zmix
                        fninz_cno(i,k) = fninz_opalGS98(i,j)
                        xiz_cno(i,k) = xiz_opalGS98(i,j)
                     enddo
c
                  enddo
c			  ! Else: if composition differs from OPAL mix used:
               else
c					! get modified CNO-mixes
                  do i = 1, nel_zmix
                     fninz_cno(i,1) = fninz_mix(i)
                     fninz_cno(i,2) = fninz_mix(i)
                     fninz_cno(i,3) = fninz_mix(i)
                     fninz_cno(i,4) = fninz_mix(i)
                  enddo
c					! and CNO-modification factors
                  do k = 2, 4
c
                     j = k + n_zmixes
c
                     f_c = min( 1.0d0 , fninz_opalGS98(1,j)
     $                    / fninz_opalGS98(1,n_cnobeg) )
                     f_n = min( 1.0d0 , fninz_opalGS98(2,j)
     $                    / fninz_opalGS98(2,n_cnobeg) )
                     f_o = min( 1.0d0 , fninz_opalGS98(kel_o,j)
     $                    / fninz_opalGS98(kel_o,n_cnobeg) )
c
                     fninz_cno(1,k) = f_c * fninz_mix(1)
                     fninz_cno(kel_o,k) = f_o * fninz_mix(kel_o)
c
                     if ( f_n .lt. 1.0 ) then
c
                        fninz_cno(2,k) = f_n * fninz_mix(2)
                        fninz_cno(4,k) = fninz_mix(4)
     $                       + ( 1.d0 - f_c ) * fninz_mix(1)
     $                       + ( 1.d0 - f_n ) * fninz_mix(2)
     $                       + ( 1.d0 - f_o ) * fninz_mix(kel_o)
c
                        fcno_mul(k) = ( fninz_cno(4,k) - fninz_mix(4) )
     $                       / ( ( 1.d0 - f_c )
     $                       * fninz_opalGS98(1,n_cnobeg)
     $                       + ( 1.d0 - f_n )
     $                       * fninz_opalGS98(2,n_cnobeg)
     $                       + ( 1.d0 - f_o )
     $                       * fninz_opalGS98(kel_o,n_cnobeg) )
c
                     else
c
                        del_co_orig =
     $                       ( 1.d0 - f_c ) * fninz_opalGS98(1,n_cnobeg)
     $                       + ( 1.d0 - f_o )
     $                       * fninz_opalGS98(kel_o,n_cnobeg)
                        fad_ne = max( 0.0d0 , min( 1.0d0 ,
     $                       ( fninz_opalGS98(4,j)
     $                       - fninz_opalGS98(4,n_cnobeg) )
     $                       / del_co_orig ) )
                        del_co = ( 1.d0 - f_c ) * fninz_mix(1)
     $                       + ( 1.d0 - f_o ) * fninz_mix(kel_o)
                        fninz_cno(2,k) = fninz_mix(2)
     $                       + ( 1.d0 - fad_ne ) * del_co
                        fninz_cno(4,k) = fninz_mix(4) + fad_ne * del_co
c
                        fcno_mul(k) = del_co / del_co_orig
c
                     endif
c
                     sum_aini = 0.0d0
                     do i = 1, nel_zmix
                        xiz_cno(i,k) = fninz_cno(i,k) * atwt_opalGS98(i)
                        sum_aini = sum_aini + xiz_cno(i,k)
                     enddo
                     do i = 1, nel_zmix
                        xiz_cno(i,k) = xiz_cno(i,k) / sum_aini
                     enddo
c
                  enddo
c
               endif
c
               x1 = fninz_cno(1,1)
               x2 = fninz_cno(1,2)
               x3 = fninz_cno(1,3)
               x4 = fninz_cno(1,4)
c
               y1 = fninz_cno(2,1)
               y2 = fninz_cno(2,2)
               y3 = fninz_cno(2,3)
               y4 = fninz_cno(2,4)
c
               z1 = fninz_cno(3,1)
               z2 = fninz_cno(3,2)
               z3 = fninz_cno(3,3)
               z4 = fninz_cno(3,4)
c
               d = ( x2 - x1 ) * ( ( y3 - y1 ) * ( z4 - z1 )
     $                             - ( y4 - y1 ) * ( z3 - z1 ) )
     $              + ( x3 - x1 ) * ( ( y4 - y1 ) * ( z2 - z1 )
     $                             - ( y2 - y1 ) * ( z4 - z1 ) )
     $              + ( x4 - x1 ) * ( ( y2 - y1 ) * ( z3 - z1 )
     $                             - ( y3 - y1 ) * ( z2 - z1 ) )
c
               if ( d .eq. 0.0d0 ) stop
     $              ' STOP -- READCO Error: CNO-interp: D = 0. '
c
               fcno_fac(0,1) = ( x2 * ( y3 * z4 - y4 * z3 )
     $              + x3 * ( y4 * z2 - y2 * z4 )
     $              + x4 * ( y2 * z3 - y3 * z2 ) ) / d
               fcno_fac(1,1) = ( y2 * ( z4 - z3 )
     $              + y3 * ( z2 - z4 ) + y4 * ( z3 - z2 ) ) / d
               fcno_fac(2,1) = ( x2 * ( z3 - z4 )
     $              + x3 * ( z4 - z2 ) + x4 * ( z2 - z3 ) ) / d
               fcno_fac(3,1) = ( x2 * ( y4 - y3 )
     $              + x3 * ( y2 - y4 ) + x4 * ( y3 - y2 ) ) / d
c
               fcno_fac(0,2) = ( x1 * ( y4 * z3 - y3 * z4 )
     $              + x3 * ( y1 * z4 - y4 * z1 )
     $              + x4 * ( y3 * z1 - y1 * z3 ) ) / d
               fcno_fac(1,2) = ( y1 * ( z3 - z4 )
     $              + y3 * ( z4 - z1 ) + y4 * ( z1 - z3 ) ) / d
               fcno_fac(2,2) = ( x1 * ( z4 - z3 )
     $              + x3 * ( z1 - z4 ) + x4 * ( z3 - z1 ) ) / d
               fcno_fac(3,2) = ( x1 * ( y3 - y4 )
     $              + x3 * ( y4 - y1 ) + x4 * ( y1 - y3 ) ) / d
c
               fcno_fac(0,3) = ( x1 * ( y2 * z4 - y4 * z2 )
     $              + x2 * ( y4 * z1 - y1 * z4 )
     $              + x4 * ( y1 * z2 - y2 * z1 ) ) / d
               fcno_fac(1,3) = ( y1 * ( z4 - z2 )
     $              + y2 * ( z1 - z4 ) + y4 * ( z2 - z1 ) ) / d
               fcno_fac(2,3) = ( x1 * ( z2 - z4 )
     $              + x2 * ( z4 - z1 ) + x4 * ( z1 - z2 ) ) / d
               fcno_fac(3,3) = ( x1 * ( y4 - y2 )
     $              + x2 * ( y1 - y4 ) + x4 * ( y2 - y1 ) ) / d
c
               fcno_fac(0,4) = ( x1 * ( y3 * z2 - y2 * z3 )
     $              + x2 * ( y1 * z3 - y3 * z1 )
     $              + x3 * ( y2 * z1 - y1 * z2 ) ) / d
               fcno_fac(1,4) = ( y1 * ( z2 - z3 )
     $              + y2 * ( z3 - z1 ) + y3 * ( z1 - z2 ) ) / d
               fcno_fac(2,4) = ( x1 * ( z3 - z2 )
     $              + x2 * ( z1 - z3 ) + x3 * ( z2 - z1 ) ) / d
               fcno_fac(3,4) = ( x1 * ( y2 - y3 )
     $              + x2 * ( y3 - y1 ) + x3 * ( y1 - y2 ) ) / d
c
            endif
c
         endif
c
      endif
c			! get number and mass fraction heavier than Ne (in Z)
      fninz_heavy = 0.0d0
      xiz_heavy = 0.0d0
      d_fninz_u_heavy = 0.0d0
      s_ninzai_mix = fninz_mix(1) * atwt_opalGS98(1)
     $     + fninz_mix(2) * atwt_opalGS98(2)
     $     + fninz_mix(3) * atwt_opalGS98(3)
     $     + fninz_mix(4) * atwt_opalGS98(4)
      ds_ninzai_u = d_fninz_user(1) * atwt_opalGS98(1)
     $     + d_fninz_user(2) * atwt_opalGS98(2)
     $     + d_fninz_user(3) * atwt_opalGS98(3)
     $     + d_fninz_user(4) * atwt_opalGS98(4)
      do i = 5, nel_zmix
         fninz_heavy = fninz_heavy + fninz_mix(i)
         xiz_heavy = xiz_heavy + xiz_mix(i)
         d_fninz_u_heavy = d_fninz_u_heavy + d_fninz_user(i)
         s_ninzai_mix = s_ninzai_mix + fninz_mix(i) * atwt_opalGS98(i)
         ds_ninzai_u = ds_ninzai_u + d_fninz_user(i) * atwt_opalGS98(i)
      enddo
      fn_o_over_cno = fninz_mix(3)
     $     / ( fninz_mix(1) + fninz_mix(2) + fninz_mix(3) )
      fninz_co_mix = fninz_mix(1) + fninz_mix(3)
c					        ! set internal CNO-interp flags
      kdo_cno = kavail_cno * kuse_cno
      kdo_user = kavail_user * kuse_user
c					 ! set flag "finshed reading opacities"
      itime = 12345678
c
      return
      end
c
c******************************************************************************
c
      SUBROUTINE SPLINE(X,Y,N,Y2)
c     ===========================
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c___
      PARAMETER ( NMAX=100 )
C
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
C
C     FIRST DERIVATIVES AT END POINTS USING CUBIC FIT
C===
      YP1 = ((Y(3)-Y(1))*(X(2)-X(1))**2.d0
     $     -(Y(2)-Y(1))*(X(3)-X(1))**2.d0)/
     $     ((X(3)-X(1))*(X(2)-X(1))*(X(2)-X(3)))
      YPN = ((Y(N-2)-Y(N))*(X(N-1)-X(N))**2.d0
     $     -(Y(N-1)-Y(N))*(X(N-2)-X(N))**2.d0)/
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
c
c******************************************************************************
c
      SUBROUTINE SPLINT(XA,YA,N,Y2A,X,Y,YP)
c     =====================================
c___
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

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
     $     ((A**3.d0-A)*Y2A(KLO)+(B**3.d0-B)*Y2A(KHI))*(H**2.d0)/6.d0
      YP = 0.05d0*((-YA(KLO)+YA(KHI))/H
     $     +      ( -(3.d0*A**2.d0-1.d0)*Y2A(KLO)
     $     +(3.d0*B**2.d0-1.d0)*Y2A(KHI) )*H/6.d0 )
      RETURN
      END
c
c******************************************************************************
c
      SUBROUTINE FITY
c     ===============
C
C  THIS ROUTINE MAKES SPLINE FITS FOR F AND FX, AND OBTAINS
C  FY AND FXY
C     ! modified:

      integer IPR,i,j,nrl,nset
      double precision A,B,AD,BD,APN,AP1,BP1,BPN,F,FX,FY,FXY
      double precision tmax,rls
      COMMON/CST_OPAL_Z/ NRL,RLS,nset,tmax
      save /CST_OPAL_Z/
C
!!      PARAMETER ( IPR=20 )
      PARAMETER ( IPR=39 )
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
c
c******************************************************************************
c
      SUBROUTINE FITX
c     ===============
C
C  THIS ROUTINE IS USED ONLY AFTER SMOOTHING.
C  ITS FUNCTION IS TO RECOMPUTE FX USING SMOOTHED F.
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
C
!!      PARAMETER ( IPR=20 )
      PARAMETER ( IPR=39 )
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
C
c******************************************************************************
c
      SUBROUTINE GETD(F,N,D,FP1,FPN)
c     ==============================
C
C  SIMPLIFIED CODE FOR SPLINE COEFFICIENTS, FOR CASE OF INTERVALS OF UNITY.
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
C___
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
C
c******************************************************************************
c
      SUBROUTINE INTERP(FLT,FLRHO,G,DGDT,DGDRHO,IERR)
c     ===============================================
C
C  GIVEN F,FX,FY AND FXY ON THE GRID POINTS, THIS ROUTINE
C  DOES BI-CUBIC INTERPOLATIONS USING METHODS DESCRIBED IN
C  Numerical Recipes, PP. 118 TO 120
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
C
!!      PARAMETER ( IPR=20 )
      PARAMETER ( IPR=39 )
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
         I = X
      ENDIF
      U = X-I
C
      IF ( Y .LT. 2. ) THEN
         IF ( Y .LT. 0.75d0 ) THEN
            IERR = .TRUE.
         ELSE
            J = 1
         ENDIF
      ELSE IF ( Y .GT. NRL-1 ) THEN
         IF ( Y .GT. NRL+0.25d0 ) THEN
            IERR = .TRUE.
         ELSE
            J = NRL-1
         ENDIF
      ELSE
         J = Y
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
      B(3) = 3.d0*(-F(I,J)+F(I,J+1))-2.d0*FY(I,J)-FY(I,J+1)
      B(4) = 2.d0*(F(I,J)-F(I,J+1))+FY(I,J)+FY(I,J+1)
C
      B(5) = FX(I,J)
      B(6) = FXY(I,J)
      B(7) = 3.d0*(-FX(I,J)+FX(I,J+1))-2.d0*FXY(I,J)-FXY(I,J+1)
      B(8) = 2.d0*(FX(I,J)-FX(I,J+1))+FXY(I,J)+FXY(I,J+1)
C
      B(9) = 3.d0*(-F(I,J)+F(I+1,J))-2.d0*FX(I,J)-FX(I+1,J)
      B(10) = 3.d0*(-FY(I,J)+FY(I+1,J))-2.d0*FXY(I,J)-FXY(I+1,J)
      B(11) = 9.d0*(F(I,J)-F(I+1,J)+F(I+1,J+1)-F(I,J+1))
     $     +6.d0*(FX(I,J)-FX(I,J+1)+FY(I,J)-FY(I+1,J))
     $     +4.d0*FXY(I,J)
     $     +3.d0*(FX(I+1,J)-FX(I+1,J+1)-FY(I+1,J+1)+FY(I,J+1))
     $     +2.d0*(FXY(I,J+1)+FXY(I+1,J))
     $     +FXY(I+1,J+1)
      B(12) = 6.d0*(-F(I,J)+F(I+1,J)-F(I+1,J+1)+F(I,J+1))
     $     +4.d0*(-FX(I,J)+FX(I,J+1))
     $     +3.d0*(-FY(I,J)+FY(I+1,J)+FY(I+1,J+1)-FY(I,J+1))
     $     +2.d0*(-FX(I+1,J)+FX(I+1,J+1)-FXY(I,J)-FXY(I,J+1))
     $     -FXY(I+1,J)-FXY(I+1,J+1)
C
      B(13) = 2.d0*(F(I,J)-F(I+1,J))+FX(I,J)+FX(I+1,J)
      B(14) = 2.d0*(FY(I,J)-FY(I+1,J))+FXY(I,J)+FXY(I+1,J)
      B(15) = 6.d0*(-F(I,J)+F(I+1,J)-F(I+1,J+1)+F(I,J+1))
     $     +4.d0*(-FY(I,J)+FY(I+1,J))
     $     +3.d0*(-FX(I,J)-FX(I+1,J)+FX(I+1,J+1)+FX(I,J+1))
     $     +2.d0*(FY(I+1,J+1)-FY(I,J+1)-FXY(I,J)-FXY(I+1,J))
     $     -FXY(I+1,J+1)-FXY(I,J+1)
      B(16) = 4.d0*(F(I,J)-F(I+1,J)+F(I+1,J+1)-F(I,J+1))
     $     +2.d0*(FX(I,J)+FX(I+1,J)-FX(I+1,J+1)-FX(I,J+1)
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
c
c******************************************************************************
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
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
C
C
!!      PARAMETER ( IPR=20 )
      PARAMETER ( IPR=39 )
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
     $     +0.0195918367d0,0.2644897959d0,+0.1175510204d0,
     $     -0.0783673469d0,+0.0277551020d0,0.3746938776d0,
     $     +0.1665306122d0,-0.1110204082d0/
      DATA ALP/
     $     -0.0844897959d0,-0.0048979592d0,+0.0073469388d0,
     $     +0.0012244898d0,0.3379591837d0,+0.0195918367d0,
     $     -0.0293877551d0,+0.4787755102d0,0.0277551020d0,
     $     -0.0416326531d0,-0.0069387755d0/
C===
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
C  SUBROUTINE INTERP
C     AFTER PROCESSING, DATA ARE IN A FORM FOR USE OF
C               SUBROUTINE INTERP
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
      implicit none
c      implicit double precision (a-h,o-z)
c      implicit integer (i-n)
C

      integer ip,ipr,i,j,ns,il,k,nz,mx,mc,mo,nrm,nrl,nrb,nre,nr,ntm,ntb
     &     ,nt,nrm_p1,nset,ntemp,nsm,nrlow,nrhigh
      double precision rle,t6,t6arr,coff,tmax,f,fx,fy,fxy,u,rossl,v,v2
     &     ,rls,flt,flrho,dgdt,g,dgdrho,flr

!!      PARAMETER ( IP=100, IPR=20 )
      PARAMETER ( IP=100, IPR=39 )
c
!!      parameter ( nz=14, mx=5, mc=8, mo=8, nrm=19, nrb=1, nre=19,
!!     $     nr=nre+1-nrb, ntm=70, ntb=1, nt=ntm+1-ntb, nrm_p1=nrm+1 )
      parameter ( nz=5, mx=5, mc=8, mo=8, nrm=39, nrb=1, nre=39,
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
            FLT = 3.75d0+0.05d0*dble(I)
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
C       CALL INTERP(FLT,FLRHO,G,DGDT,DGDRHO,IERR)
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
               FLR = RLS+0.5d0*dble(IL-1)
               FLRHO = FLR-18.d0+3.d0*FLT
               CALL INTERP(FLT,FLRHO,G,DGDT,DGDRHO,IERR)
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
c
      subroutine open_chk_zip( iu, fname, igzip, cmsg )
c     -------------------------------------------------
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      character*(*) fname
      character*(*) cmsg
c
      character*511 ctmp
c
      logical lxst
c===
      last = lnblnk( fname )
c
      if ( last .le. 0 ) stop
     $     ' STOP -- OPEN_CHK_ZIP Error: blank file name. '
c
      call inqfil( fname, lxst )
c
      if ( lxst ) then
         igzip = 0
      else if ( last .gt. 508 ) then
         stop ' STOP -- OPEN_CHK_ZIP Error: file name too long. '
      else
         ctmp = fname(1:last) // '.gz'
         call inqfil( ctmp, lxst )
         if ( lxst ) then
            igzip = 1
            ctmp = 'gunzip ' // fname(1:last) // '.gz'
            call system( ctmp )
         else
            ctmp = fname(1:last) // '.Z'
            call inqfil( ctmp, lxst )
            if ( lxst ) then
               igzip = -1
               ctmp = 'uncompress ' // fname(1:last) // '.Z'
               call system( ctmp )
            else
               write(6,'(" ",a)') cmsg
c               print *,'tututut'
               write(6,'(" ",a)') fname
               stop ' STOP -- READCO Error: file not found. '
            endif
         endif
      endif
c
      call opoldr( iu, fname )
c
      return
      end
c
c************************************************************************
c
      subroutine close_chk_zip( iu, fname, igzip )
c     --------------------------------------------
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      character*(*) fname
c
      character*511 ctmp
c
      close(iu)
c
      if ( igzip .gt. 0 ) then
c
         ctmp = 'gzip ' // fname
         call system( ctmp )
c
      else if ( igzip .lt. 0 ) then
c
         ctmp = 'compress ' // fname
         call system( ctmp )
c
      endif
c
      return
      end
c
c************************************************************************
c
      integer function num_blanks_contained( fname )
c     --------------------------------------
c
c  Count the number of blanks between first and last non-blank characters
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      character*(*) fname
c
      character*1 ctab
      common /c_ctab/ ctab
      save /c_ctab/
c
      iblank = 0
c
      if ( lnblnk(fname) .gt. 0 ) then
c
         do i = non_blank_begin(fname), lnblnk(fname)
            if ( fname(i:i) .eq. ' ' .or. fname(i:i) .eq. ctab )
     $           iblank = iblank + 1
         enddo
c
      endif
c
      num_blanks_contained = iblank
c
      return
      end
c
c************************************************************************
c
      integer function non_blank_begin( fname )
c     ---------------------------------
c
c  Find the first non-blank character in the input character variable
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      character*(*) fname
c
      character*1 ctab
      common /c_ctab/ ctab
      save /c_ctab/
c
      last = lnblnk( fname )
c
      if ( last .le. 1 ) then
         i = last
      else
         i = 1
         do while ( i .lt. last .and.
     $        ( fname(i:i) .eq. ' ' .or. fname(i:i) .eq. ctab ) )
            i = i + 1
         enddo
      endif
c
      non_blank_begin = i
c
      return
      end
c
c************************************************************************
c  NOTE that the subroutines below have alternate parts, all but one
c  commented out, for various flavors of UNIX and for VMS (the Linux
c  versions should actually work for any flavor of UNIX).
c***********************************************************************
c
      subroutine opoldr(iu,fname)
c     ---------------------------
c
c Open an old formatted file:
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      character*(*) fname
c
      character*1 cb(6)
      common/chkpoc/cb
      save /chkpoc/
c
c-linux[
      character*255 fnalt
c-linux]
c
c  For Linux: get home directory name if necessary, and open the file
c  with the err= keyword to prevent coredump
c  (actually, this should work on any Unix system, provided that the
c  environment variable HOME is correctly defined as the home directory):
c-linux[
      call linux_get_home_dir(fname,fnalt,ialt)
c
      if ( ialt .gt. 0 ) then
         open(iu,file=fnalt,form='FORMATTED',status='OLD',
     $        iostat=ioperr,err=900)
      else
         open(iu,file=fname,form='FORMATTED',status='OLD',
     $        iostat=ioperr,err=900)
      endif
c
      write(*,*) 'reading opacity file: ',trim(fname)
      return
c
 900  write(6,910) ioperr,iu,fname(:lnblnk(fname))
 910  format(' '/' Error',i12,' opening unit',i3,' with old file:'/
     $     ' ',a)
      stop ' STOP -- OPOLDR: Error opening old file. '
c-linux]
c
c  For Sun UNIX: open the file:
c-sun[
c-sun;      open(iu,file=fname,form='FORMATTED',status='OLD')
c-sun;      return
c-sun]
c
c  For VMS, or for Iris UNIX: open the file as read-only:
c-vms-iris[
c-vms-iris;      open(iu,file=fname,form='FORMATTED',status='OLD',
c-vms-iris;     $     readonly)
c-vms-iris;      return
c-vms-iris]
c
      end
c
c************************************************************************
c
      subroutine opoluf(iu,fname)
c     ---------------------------
c
c Open an old unformatted file:
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      character*(*) fname
c
c  For Linux: open the file, with the err= keyword to prevent coredump:
c-linux[
      character*255 fnalt
c
      call linux_get_home_dir(fname,fnalt,ialt)
c
      if ( ialt .gt. 0 ) then
         open(iu,file=fnalt,form='UNFORMATTED',status='OLD',
     $        iostat=ioperr,err=900)
      else
         open(iu,file=fname,form='UNFORMATTED',status='OLD',
     $        iostat=ioperr,err=900)
      endif
c
      return
c
 900  write(6,910) ioperr,iu,fname(:lnblnk(fname))
 910  format(' '/' Error',i12,' opening unit',i3,
     $     ' with old unformatted file:'/' ',a)
      stop ' STOP -- OPOLUF: Error opening old unformatted file. '
c-linux]
c
c  For Sun UNIX: open the file:
c-sun[
c-sun;      open(iu,file=fname,form='UNFORMATTED',status='OLD')
c-sun;      return
c-sun]
c
c  For VMS or Iris UNIX: open the file as read-only:
c-vms-iris[
c-vms-iris;      open(iu,file=fname,form='UNFORMATTED',status='OLD',
c-vms-iris;     $     readonly)
c-vms-iris;      return
c-vms-iris]
c
      end
c
c************************************************************************
c
      subroutine opneuf(iu,fname)
c     ---------------------------
c
c Open a new unformatted file:
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      character*(*) fname
c
c  For Linux: open the file, with the err= keyword to prevent coredump:
c-linux[
      character*255 fnalt
c
      call linux_get_home_dir(fname,fnalt,ialt)
c
      if ( ialt .gt. 0 ) then
         open(iu,file=fnalt,form='UNFORMATTED',status='UNKNOWN',
     $        iostat=ioperr,err=900)
      else
         open(iu,file=fname,form='UNFORMATTED',status='UNKNOWN',
     $        iostat=ioperr,err=900)
      endif
c
      return
c
 900  write(6,910) ioperr,iu,fname(:lnblnk(fname))
 910  format(' '/' Error',i12,' opening unit',i3,
     $     ' with new unformatted file:'/' ',a)
      stop ' STOP -- OPNEUF: Error opening new unformatted file. '
c-linux]
c
c  For UNIX: open the file status UNKNOWN so not an error if file exists:
c-sun-iris[
c-sun-iris;      open(iu,file=fname,form='UNFORMATTED',status='UNKNOWN')
c-sun-iris;      return
c-sun-iris]
c
c  For VMS: open the file status NEW:
c-vms[
c-vms;      open(iu,file=fname,form='UNFORMATTED',status='NEW')
c-vms;      return
c-vms]
c
      end
c
c************************************************************************
c
      subroutine inqfil(fname,lxst)
c     -----------------------------
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      character*(*) fname
      logical lxst
c
c  For Linux: get home directory name, if necessary
c  (actually, this should work on any Unix system, provided that the
c  environment variable HOME is correctly defined as the home directory):
c-linux[
      character*255 fnalt
c
      call linux_get_home_dir(fname,fnalt,ialt)
c
      if ( ialt .gt. 0 ) then
         inquire( file = fnalt, exist = lxst )
      else
         inquire( file = fname, exist = lxst )
      endif
c-linux]
c
c  Anything else: just look for filename as is:
c-sun-vms-iris[
c-sun-vms-iris;      inquire( file = fname , exist = lxst )
c-sun-vms-iris]
c
      return
      end
c
c************************************************************************
c
c-linux[
      subroutine linux_get_home_dir(fname,fnalt,ialt)
c     -----------------------------------------------
c
c  For Linux, at least with fort77, the prefix '~' in a filename is not
c  recognized as "home directory": get it from HOME environment variable.
c  (Actually, this should work on any Unix system, provided that the
c  environment variable HOME is correctly defined as the home directory.)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      character*(*) fname
      character*(*) fnalt
c
      if ( len(fname) .ge. 2 .and. fname(1:2) .eq. '~/' ) then
c
         ialt = 1
         call getenv( 'HOME' , fnalt )
         i = lnblnk(fnalt)
         if ( lnblnk(fname) + i - 1 .gt. len(fnalt) ) then
            write(6,900) lnblnk(fname) + i - 1, len(fnalt),
     $           fnalt(1:i), fname(2:lnblnk(fname))
 900        format(' '/' Error: filename has',i6,
     $           ' characters >',i6,' --- too long:'/' ',a,a)
            stop ' STOP -- LINUX_GET_HOME_DIR: filename too long. '
         endif
         fnalt(i+1:) = fname(2:)
c
      else
c
         ialt = 0
c
      endif
c
      return
      end
c-linux]
c
c************************************************************************
c
c-linux[
      integer function lnblnk(fname)
c     ----------------------
c
c  Needed for fort77 under Linux, since fort77 linker can't find  lnblnk_
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      character*(*) fname
c
      character*1 ctab
      common /c_ctab/ ctab
      save /c_ctab/
c
      i = len(fname)
c
      do while( i .gt. 1 .and.
     $     ( fname(i:i) .eq. ' ' .or. fname(i:i) .eq. ctab ) )
         i = i - 1
      enddo
      if ( fname(i:i) .eq. ' ' .or. fname(i:i) .eq. ctab ) i = i - 1
c
      lnblnk = i
c
      return
      end
c-linux]
c
c************************************************************************
c
