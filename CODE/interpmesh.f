

************************************************************************

      SUBROUTINE interpmesh (il,ir,ls)

************************************************************************
* Interpolate variables in case of shell addition/removal              *
*                                                                      *
* delete shell                                                         *
* ------------                                                         *
* before    |      |     |    |        |                               *
*          il-1    il  il+1   ir     ir+1                              *
* after     |      |          |        |                               *
*          il-1    il         ir     ir+1                              *
*                                                                      *
* add shell                                                            *
* ---------                                                            *
* before    |      |          |        |                               *
*          il-1    il         ir     ir+1                              *
* after     |      |     |    |        |                               *
*          il-1    il  il+1   ir     ir+1                              *
*                                                                      *
* ls = 0 : shell shifting                                              *
* ls = 1 : shell suppression                                           *
* ls = 2 : shell addition, interpolate il+1 from il and ir (default)   *
* ls = 3 : shell addition, interpolate il+1 from il-1 and ir           *
* ls = 4 : shell addition, interpolate il+1 from il-1 and il           *
*   treatment of mass loss/accretion                                   *
* ls = 5 : rescaling the pseudo-lagrangian coordinate                  *
*   treatment of shock fronts                                          *
* ls = 6 : shell addition, interpolate il from ir and ir+1             *
* ls = 7 : shell addition, special treatment of central shell          *
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.acc'
      include 'evolcom.cons'
      include 'evolcom.eng'
      include 'evolcom.ion'
      include 'evolcom.mass'
      include 'evolcom.mod'
      include 'evolcom.nuc'
      include 'evolcom.rot'
      include 'evolcom.conv'
      include 'evolcom.spec'
      include 'evolcom.surf'
      include 'evolcom.teq'
      include 'evolcom.therm'
      include 'evolcom.transp'
      include 'evolcom.var'
      include 'evolcom.igw' !LR 20250617

      integer ls,il,il0,il1,ir,ir1,ir2
      integer ia,ib,ic,id,iab,icd,ka,kb
      integer iqmrbot,iqmrtop
      integer k,l
      integer ja,jb,ixa,ixb,ixc,ixd

      double precision ai,bi,ci,di,del,ak,bk
      double precision midmr,vmidmr,qmr,vqmr
      double precision r2dm0,r2dm1,r2dm2,r2dm3,r2dm4
      double precision zj1,zj2,zj3,rm1,rm2,rm3
      double precision xai,xbi,xci,xdi
      double precision domj,vconv

      common /newmesh/ qmr(nsh),vqmr(nsh),iqmrbot,iqmrtop


      il0 = il-1
      il1 = il+1
      ir1 = ir+1
      ir2 = ir+2


*__________________________
*
***      SHELL SHIFTING
*
*--------------------------


      if (ls.eq.0) then
***   interface variables
         vu(il) = vu(ir)
         vlum(il) = vlum(ir)
         vr(il) = vr(ir)

         u(il) = u(ir)
         r(il) = r(ir)
         lum(il) = lum(ir)

         vsconv(il) = vsconv(ir)
         vhconv(il) = vhconv(ir)
         mr(il) = mr(ir)
         m(il) = m(ir)

***   centered variables
         lnf(il) = lnf(ir)
         lnt(il) = lnt(ir)
         ro(il) = ro(ir)
         p(il) = p(ir)
         e(il) = e(ir)
         s(il) = s(ir)

         vlnf(il) = vlnf(ir)
         vlnt(il) = vlnt(ir)
         vro(il) = vro(ir)
         vp(il) = vp(ir)
         ve(il) = ve(ir)
         vs(il) = vs(ir)

         crz(il) = crz(ir)
         dm(il) = dm(ir)
         tau(il) = tau(ir)
         vmueinv(il) = vmueinv(ir)
         venucl(il) = venucl(ir)
         vpvisc(il) = vpvisc(ir)
         exposure(il) = exposure(ir)

         do l = 1,nsp
            xsp(il,l) = xsp(ir,l)
         enddo
         vomega(il) = vomega(ir)
         if (irotbin.eq.1) then
            auxs(il) = auxs(ir)
            urs(il) = urs(ir)
            xpsis(il) = xpsis(ir)
            xlambdas(il) = xlambdas(ir)
         endif
         if (iaccr.gt.0) then
            facc(il) = facc(ir)
            macc(il) = macc(ir)
            vfacc(il) = vfacc(ir)
            dfaccdr(il) = dfaccdr(ir)
            vacc(il) = vacc(ir)
            tacc(il) = tacc(ir)
         endif
      endif



*___________________________________________________________
***
***                       SUPPRESSION
***
***   adjustment of shells adjacent to a suppressed one
***   only applies to centered variables (.i.e. not lum,r,u)
***   new version
*-----------------------------------------------------------


      if (ls.eq.1) then
         ai = dm(il)/(dm(il)+dm(ir))
         bi = 1.d0-ai

         lnf(il) = ai*lnf(il)+bi*lnf(ir)
         lnt(il) = ai*lnt(il)+bi*lnt(ir)
         p(il) = exp(ai*log(p(il))+bi*log(p(ir)))
         e(il) = exp(ai*log(e(il))+bi*log(e(ir)))
         if (s(il).gt.0.d0.and.s(ir).gt.0.d0) then
            s(il) = exp(ai*log(s(il))+bi*log(s(ir)))
         else
            s(il) = ai*s(il)+bi*s(ir)
         endif
         ro(il) = exp(ai*log(ro(il))+bi*log(ro(ir)))

         vlnf(il) = ai*vlnf(il)+bi*vlnf(ir)
         vlnt(il) = ai*vlnt(il)+bi*vlnt(ir)
         vp(il) = exp(ai*log(vp(il))+bi*log(vp(ir)))
         ve(il) = exp(ai*log(ve(il))+bi*log(ve(ir)))
         if (vs(il).gt.0.d0.and.vs(ir).gt.0.d0) then
            vs(il) = exp(ai*log(vs(il))+bi*log(vs(ir)))
         else
            vs(il) = ai*vs(il)+bi*vs(ir)
         endif
         vro(il) = exp(ai*log(vro(il))+bi*log(vro(ir)))

         tau(il) = ai*tau(il)+bi*tau(ir)
         vmueinv(il) = ai*vmueinv(il)+bi*vmueinv(ir)
         venucl(il) = ai*venucl(il)+bi*venucl(ir)
         vpvisc(il) = ai*vpvisc(il)+bi*vpvisc(ir)
c         exposure(il) = ai*exposure(il)+bi*exposure(ir)
         exposure(il) = max(exposure(il),exposure(ir))

         if (xsp(il,2).ne.xsp(ir,2)) then
            do l = 1,nsp
               xsp(il,l) = ai*xsp(il,l)+bi*xsp(ir,l)
            enddo
         endif

***   Interpolation ensures angular momentum conservation

         if (irotbin.eq.1) then
c           vomega(il) = ai*vomega(il)+bi*vomega(ir)
            if (il.gt.1) then 
               r2dm0 = dm(il0)*(r(il0)**2+r(il)*r(il0)+r(il)**2)/3.d0
            else
               r2dm0 = 0.d0
            endif
            r2dm3 = (dm(il)+dm(ir))*(r(il)**2+r(il)*r(ir1)+r(ir1)**2)
     &           /3.d0
            r2dm1 = dm(il)*(r(il)**2+r(il)*r(ir)+r(ir)**2)/3.d0
            r2dm2 = dm(ir)*(r(ir)**2+r(ir)*r(ir1)+r(ir1)**2)/3.d0
            r2dm4 = dm(ir1)*(r(ir2)**2+r(ir2)*r(ir1)+r(ir1)**2)/3.d0

c..   Increment domj is applied to omega at shells il and ir1 in order
c..   to allow angular momentum conservation when shell ir is suppressed
c            domj = (r2dm1*(vomega(il)+vomega(ir))+r2dm2*(vomega(ir)+
c     &           vomega(ir1))-r2dm3*(vomega(il)+vomega(ir1)))/
c     &           (r2dm0+r2dm4+2.d0*r2dm3)
c            vomega(il) = vomega(il)+domj
c            vomega(ir1) = vomega(ir1)+domj

c..   New formulation for the redistribution of omega
            domj = (vomega(ir)*(r2dm2+r2dm1)+vomega(il)*(r2dm1-r2dm3)+
     &           vomega(ir1)*(r2dm2-r2dm3))/(vomega(ir1)*(r2dm3+r2dm4)+
     &           vomega(il)*(r2dm3+r2dm0))

            vomega(il) = vomega(il)*(1.d0+domj)
            vomega(ir1) = vomega(ir1)*(1.d0+domj)

            auxs(il) = ai*auxs(il)+bi*auxs(ir)
            urs(il) = ai*urs(il)+bi*urs(ir)
            xpsis(il) = ai*xpsis(il)+bi*xpsis(ir)
            xlambdas(il) = ai*xlambdas(il)+bi*xlambdas(ir)
         endif
      endif



*______________________________________________________
***
***                       ADDITION
***
***   interpolation of the variables in the added shell
***   define interpolation coefficients
*------------------------------------------------------



      if (ls.eq.2.or.ls.eq.3.or.ls.eq.5.or.ls.eq.6.or.ls.eq.7) then

***   interpolation at il+1 with variables defined at il and ir
*     default
         if (ls.eq.2) then
            ka = il
            kb = ir
            ak = 0.5d0
            bk = 0.5d0
c..   centered variables
            iab = il1
            bi = dm(il1)/(dm(ir)+2.d0*dm(il1))
            ai = 1.d0-bi
            ia = il
            ib = ir
            icd = il
            di = dm(il)/(dm(il0)+2.d0*dm(il))
            ci = 1.d0-di
            ic = il
            id = il0
c..   chemical composition
            xai = ai
            xbi = bi
            xci = ci
            xdi = di
            ixa = ia
            ixb = ib
            ixc = ic
            ixd = id
         endif

***   interpolation at il+1 with variables defined at il-1 and ir
         if (ls.eq.3) then
            ka = il
            kb = ir
            ak = 0.5d0
            bk = 0.5d0
c..   centered variables
            del = dm(il0)+4.d0*dm(il)+dm(ir)
            del = 1.d0/del
            iab = il1
            ai = (3.d0*dm(il)+dm(ir))*del
            bi = 1.d0-ai
            ia = il0
            ib = ir
            icd = il
            ci = (dm(il)+dm(ir))*del
            di = 1.d0-ci
            ic = il0
            id = ir
c..   chemical composition
            xai = ai
            xbi = bi
            xci = ci
            xdi = di
            ixa = ia
            ixb = ib
            ixc = ic
            ixd = id
         endif

***   special treatment of inner shock-front
***   extrapolation at ir from variables defined at il-1 and il
         if (ls.eq.5) then
            ka = il0
            kb = il
            ak = -dm(il)/dm(il0)
            bk = 1.d0-ak
c..   centered variables
            iab = il
            bi = dm(il)/(dm(il0)+2.d0*dm(il))
            ai = 1.d0-bi
            ia = il
            ib = il0
            icd = il1
            ci = -2.d0*dm(il)/(dm(il)+dm(il0))
            di = 1.d0-ci
            ic = il0
            id = il
c..   chemical composition
            xai = ai
            xbi = bi
            xci = ci
            xdi = di
            ixa = ia
            ixb = ib
            ixc = ic
            ixd = id
         endif

***   special treatment of upper shock-front
***   extrapolation at il from variables defined at ir and ir+1
         if (ls.eq.6) then
            ka = ir
            kb = ir+1
            bk = -dm(il)/dm(ir)
            ak = 1.d0-bk
c..   centered variables
            iab = il1
            ai = dm(il1)/(dm(ir)+2.d0*dm(il1))
            bi = 1.d0-ai
            ia = ir
            ib = il
            icd = il
            ci = -dm(il1)/(2.d0*dm(il1)+dm(ir))
            di = 1.d0-ci
            ic = ir
            id = il
c..   chemical composition
            xai = ai
            xbi = bi
            xci = ci
            xdi = di
            ixa = ia
            ixb = ib
            ixc = ic
            ixd = id
         endif

***   special treatment of central quantities
         if (ls.eq.7) then
            ka = il
            kb = ir
            ak = 0.5d0
            bk = 0.5d0
c..   centered variables
            iab = il1
            ia = ir
            ib = il
            ai = dm(iab)/(dm(ia)+2.d0*dm(iab))
            bi = 1.d0-ai
            icd = il
            ic = il
            id = il1
            ci = 2.d0
            di = -1.d0
c..   chemical composition
            xai = ai
            xbi = bi
            xci = ci
            xdi = di
            ixa = ia
            ixb = ib
            ixc = ic
            ixd = id
         endif

*------------------------------------
*
***    INTERPOLATION (shell addition)
*
*------------------------------------

c..   interface variables
         vr(il1) = (0.5d0*(vr(il)**3+vr(ir)**3))**pw13
         vu(il1) = ak*vu(ka)+bk*vu(kb)
         vlum(il1) = ak*vlum(ka)+bk*vlum(kb)

         r(il1) = (0.5d0*(r(il)**3+r(ir)**3))**pw13
         u(il1) = ak*u(ka)+bk*u(kb)
         lum(il1) = ak*lum(ka)+bk*lum(kb)

         vconv = ak*vsconv(ka)+bk*vsconv(kb)
         if (vconv.gt.0.d0.or.ls.eq.2) then
            vsconv(il1) = vconv
            vhconv(il1) = ak*vhconv(ka)+bk*vhconv(kb)
         else
            vsconv(il1) = 0.5d0*(vsconv(il)+vsconv(ir))
            vhconv(il1) = 0.5d0*(vhconv(il)+vhconv(ir))
         endif
         crz(il1) = crz(il)

***   coefficients for omega interpolation with
***   angular momentum conservation
         if (irotbin.eq.1) then
            rm1 = (r(il)**2+r(il)*r(ir)+r(ir)**2)/3.d0
            rm2 = (r(il)**2+r(il)*r(il1)+r(il1)**2)/3.d0
            rm3 = (r(ir)**2+r(ir)*r(il1)+r(il1)**2)/3.d0
            zj1 = rm1*(dm(il)+dm(il1))
            zj2 = rm2*dm(il)
            zj3 = rm3*dm(il1)
            ja = il
            jb = ir
            vomega(il1) = (vomega(ja)*(zj1-zj2)+vomega(jb)*(zj1-zj3))/
     &           (zj2+zj3)
            auxs(il1) = ak*auxs(ka)+bk*auxs(kb)
            urs(il1) = ak*urs(ka)+bk*urs(kb)
            xpsis(il1) = ak*xpsis(ka)+bk*xpsis(kb)
            xlambdas(il1) = ak*xlambdas(ka)+bk*xlambdas(kb)
         else
            vomega(il1) = ak*vomega(ka)+bk*vomega(kb)
         endif

***   centered variables
***   WARNING : order matters, first interpolate at iab and then at icd
         lnf(iab) = ai*lnf(ia)+bi*lnf(ib)
         lnt(iab) = ai*lnt(ia)+bi*lnt(ib)
         p(iab) = exp(ai*log(p(ia))+bi*log(p(ib)))
         e(iab) = exp(ai*log(e(ia))+bi*log(e(ib)))
         s(iab) = exp(ai*log(s(ia))+bi*log(s(ib)))
         ro(iab) = exp(ai*log(ro(ia))+bi*log(ro(ib)))

         vlnf(iab) = ai*vlnf(ia)+bi*vlnf(ib)
         vlnt(iab) = ai*vlnt(ia)+bi*vlnt(ib)
         vp(iab) = exp(ai*log(vp(ia))+bi*log(vp(ib)))
         ve(iab) = exp(ai*log(ve(ia))+bi*log(ve(ib)))
         vs(iab) = exp(ai*log(vs(ia))+bi*log(vs(ib)))
         vro(iab) = exp(ai*log(vro(ia))+bi*log(vro(ib)))

         tau(iab) = ai*tau(ia)+bi*tau(ib)
         vmueinv(iab) = ai*vmueinv(ia)+bi*vmueinv(ib)
         venucl(iab) = ai*venucl(ia)+bi*venucl(ib)
         if (vpvisc(ia).gt.0.d0.and.vpvisc(ib).gt.0.d0) then
            vpvisc(iab) = exp(ai*log(vpvisc(ia))+bi*log(vpvisc(ib)))
         else
            vpvisc(iab) = ai*vpvisc(ia)+bi*vpvisc(ib)
         endif
         exposure(iab) = exposure(ia)

         lnf(icd) = ci*lnf(ic)+di*lnf(id)
         lnt(icd) = ci*lnt(ic)+di*lnt(id)
         p(icd) = exp(ci*log(p(ic))+di*log(p(id)))
         e(icd) = exp(ci*log(e(ic))+di*log(e(id)))
         if (s(ic).gt.0.d0.and.s(id).gt.0.d0) then
            s(icd) = exp(ci*log(s(ic))+di*log(s(id)))
         else
            s(icd) = ci*s(ic)+di*s(id)
         endif
         ro(icd) = exp(ci*log(ro(ic))+di*log(ro(id)))

         vlnf(icd) = ci*vlnf(ic)+di*vlnf(id)
         vlnt(icd) = ci*vlnt(ic)+di*vlnt(id)
         vp(icd) = exp(ci*log(vp(ic))+di*log(vp(id)))
         ve(icd) = exp(ci*log(ve(ic))+di*log(ve(id)))
         if (vs(ic).gt.0.d0.and.vs(id).gt.0.d0) then
            vs(icd) = exp(ci*log(vs(ic))+di*log(vs(id)))
         else
            s(icd) = ci*vs(ic)+di*vs(id)
         endif
         vro(icd) = exp(ci*log(vro(ic))+di*log(vro(id)))

         tau(icd) = ci*tau(ic)+di*tau(id)
         vmueinv(icd) = ci*vmueinv(ic)+di*vmueinv(id)
         venucl(icd) = ci*venucl(ic)+di*venucl(id)
         if (vpvisc(ic).gt.0.d0.and.vpvisc(id).gt.0.d0) then
            vpvisc(icd) = exp(ci*log(vpvisc(ic))+di*log(vpvisc(id)))
         else
            vpvisc(icd) = ci*vpvisc(ic)+di*vpvisc(id)
         endif
         exposure(icd) = exposure(iab)

***   conservative laws (mass fraction conserved)
         if (xsp(ixb,io16).ne.xsp(ixd,io16)) then
            do l = 1,nsp
c..   default
               xsp(iab,l) = xai*xsp(ixa,l)+xbi*xsp(ixb,l)
               xsp(icd,l) = xci*xsp(ixc,l)+xdi*xsp(ixd,l)
c..   mesh2
c               xsp(iab,l) = xsp(ixa,l)
c               xsp(icd,l) = xsp(iab,l)
c..   mesh4
c               xsp(iab,l) = xsp(ixb,l)
c               xsp(icd,l) = xsp(ixd,l)
c..   mesh3
c               if (crz(ia).eq.crz(ixb)) then
c                  xsp(iab,l) = xai*xsp(ixa,l)+xbi*xsp(ixb,l)
c               else
c                  xsp(iab,l) = xsp(ixa,l)
c               endif
c               if (crz(ic).eq.crz(ixd)) then
c                  xsp(icd,l) = xci*xsp(ixc,l)+xdi*xsp(ixd,l)
c               else
c                  xsp(icd,l) = xsp(ixa,l)
c               endif
c..   mesh1
c               if (crz(il0).eq.crz(ir)) then
c                  xsp(iab,l) = xai*xsp(ixa,l)+xbi*xsp(ixb,l)
c                  xsp(icd,l) = xci*xsp(ixc,l)+xdi*xsp(ixd,l)
c               else
c                  xsp(iab,l) = xsp(ixa,l)
c                  xsp(icd,l) = xsp(ixa,l)
c               endif
c..   mesh5
c              if (crz(il0).eq.crz(ir)) then
c                 xsp(iab,l) = xai*xsp(ixa,l)+xbi*xsp(ixb,l)
c                 xsp(icd,l) = xci*xsp(ixc,l)+xdi*xsp(ixd,l)
c              else
c                 xsp(iab,l) = xsp(ixb,l)
c                 xsp(icd,l) = xsp(ixd,l)
c              endif
            enddo
         endif

         if (iaccr.gt.0) then
            facc(iab) = ai*facc(ia)+bi*facc(ib)
            macc(iab) = ai*macc(ia)+bi*macc(ib)
            vfacc(iab) = ai*vfacc(ia)+bi*vfacc(ib)
            dfaccdr(iab) = ai*dfaccdr(ia)+bi*dfaccdr(ib)
            vacc(iab) = ai*vacc(ia)+bi*vacc(ib)
            tacc(iab) = ai*tacc(ia)+bi*tacc(ib)

            facc(icd) = ci*facc(ic)+di*facc(id)
            macc(icd) = ci*macc(ic)+di*macc(id)
            vfacc(icd) = ci*vfacc(ic)+di*vfacc(id)
            dfaccdr(icd) = ci*dfaccdr(ic)+di*dfaccdr(id)
            vacc(icd) = ci*vacc(ic)+di*vacc(id)
            tacc(icd) = ci*tacc(ic)+di*tacc(id)
         endif
      endif



*__________________________________________________________________
***
***               TREATMENT OF MASS LOSS/ACCRETION
***
***   change of independent variable in case of accretion/mass-loss
***   and interpolation of old variables at the new mesh point qmr
*------------------------------------------------------------------


      if (ls.eq.4) then
         ia = ir-1
         ib = ir
         ic = ir-1
         id = ir
         bi = (qmr(il)-vqmr(ir-1))/(vqmr(ir)-vqmr(ir-1))
         ai = 1.d0-bi
***   special treatment for centered variables
         midmr = 0.5d0*(qmr(il)+qmr(il+1))
         vmidmr = 0.5d0*(vqmr(ia)+vqmr(ir))
         if (midmr.lt.vmidmr) then
            ic = ia-1
            id = ir-1
            goto 20
         endif
         if (ir.lt.nmod) then
            vmidmr = 0.5d0*(vqmr(ir)+vqmr(ir+1))
            if (midmr.gt.vmidmr) then
               do k = ir+1,nmod1
                  vmidmr = 0.5d0*(vqmr(k)+vqmr(k+1))
                  if (vmidmr.gt.midmr) then
                     id = k
                     ic = id-1
                     goto 20
                  endif
               enddo
            endif
         else
            id = nmod1
            ic = nmod1
         endif
 20      di = (qmr(il)+qmr(il+1)-vqmr(ic)-vqmr(ic+1))/(vqmr(id+1)-
     &        vqmr(ic))
         ci = 1.d0-di

***   centered variable
         lnf(il) = (ci*lnf(ic)+di*lnf(id))
         lnt(il) = (ci*lnt(ic)+di*lnt(id))
         p(il) = exp(ci*log(p(ic))+di*log(p(id)))
         e(il) = exp(ci*log(e(ic))+di*log(e(id)))
         s(il) = exp(ci*log(s(ic))+di*log(s(id)))
         ro(il) = exp(ci*log(ro(ic))+di*log(ro(id)))

         vlnf(il) = (ci*vlnf(ic)+di*vlnf(id))
         vlnt(il) = (ci*vlnt(ic)+di*vlnt(id))
         vp(il) = exp(ci*log(vp(ic))+di*log(vp(id)))
         ve(il) = exp(ci*log(ve(ic))+di*log(ve(id)))
         vs(il) = exp(ci*log(vs(ic))+di*log(vs(id)))
         vro(il) = exp(ci*log(vro(ic))+di*log(vro(id)))

         tau(il) = (ci*tau(ic)+di*tau(id))
         vmueinv(il) = (ci*vmueinv(ic)+di*vmueinv(id))
         venucl(il) = (ci*venucl(ic)+di*venucl(id))
         if (vpvisc(ic).gt.0.d0.and.vpvisc(id).gt.0.d0) then
            vpvisc(il) = exp(ci*log(vpvisc(ic))+di*log(vpvisc(id)))
         else
            vpvisc(il) = ci*vpvisc(ic)+di*vpvisc(id)
         endif
         exposure(il) = (ci*exposure(ic)+di*exposure(id))

***   variables defined at interface
         u(il) = ai*u(ia)+bi*u(ib)
         r(il) = (ai*r(ia)**3+bi*r(ib)**3)**pw13
         lum(il) = ai*lum(ia)+bi*lum(ib)

         vu(il) = ai*vu(ia)+bi*vu(ib)
         vr(il) = (ai*vr(ia)**3+bi*vr(ib)**3)**pw13
         vlum(il) = ai*vlum(ia)+bi*vlum(ib)

         vsconv(il) = ai*vsconv(ia)+bi*vsconv(ib)
         vhconv(il) = ai*vhconv(ia)+bi*vhconv(ib)

         vqmr(ia) = qmr(il)

      endif


      return
      end
