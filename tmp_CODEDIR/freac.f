

************************************************************************

      SUBROUTINE freac (vv,eqd,y,dt,lflag)

************************************************************************
*  compute the nuclear matrix elements and derivatives
*  eq[i,neqd1] = dt*dYi/dt = dt*Fi(Yj)     (nucleosynthesis equation)
*  eq[i,j=neqd+1,2neqd] = dt*d(dYi/dt)/dYj      (jacobian)
*  lflag = .true.  : nucleosynthesis alone
*  lflag = .false. : coupling nucleosynthesis-mixing
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
************************************************************************

      implicit none

      include 'evolpar.star'

      include 'evolcom.nuc'

      integer neqd,neqd1
      integer kna,kia,knb,kib,knd,kid,knc,kic
      integer i,nn

      double precision dt,z,za,zb
      double precision vv(nreac),eqd(nsp,3*nsp+1),y(nsp)

      logical lflag

      neqd = nsp
      neqd1 = 3*neqd+1


c..  reaction : kna X(kia) + knb X(kib) --> knd X(kid) + knc X(kic)
c..  reaction : k2  X(k1)  + k4  X(k3)  --> k8  X(k7)  + k6  X(k5)
c          ex :  1    C12  +  1    He   -->  0  gamma  +  1   O16

      if (lflag) then
         nn = 0
      else
         nn = neqd
      endif
      do i = 1,nreac
         kna = k2(i)
         kia = k1(i)
         knb = k4(i)
         kib = k3(i)
         knd = k6(i)
         kid = k5(i)
         knc = k8(i)
         kic = k7(i)
         if (kna.eq.1.and.knb.eq.1) then
            z = vv(i)*y(kia)*y(kib)
            za = vv(i)*y(kib)
            zb = vv(i)*y(kia)
         else
            if (knb.eq.0) then
               if (kna.eq.1) then
                  z = vv(i)*y(kia)
                  za = vv(i)
               else
                  z = vv(i)*y(kia)**kna/fact(kna)
                  za = vv(i)*y(kia)**(kna-1)/fact(kna-1)
               endif
               zb = 0.d0
            else
               z = vv(i)*y(kia)**kna*y(kib)**knb/(fact(kna)*fact(knb))
               za = vv(i)*y(kia)**(kna-1)*y(kib)**knb/
     &              (fact(kna-1)*fact(knb))
               zb = vv(i)*y(kia)**kna*y(kib)**(knb-1)/
     &              (fact(kna)*fact(knb-1))
            endif
         endif
         z = z*dt
         za = za*dt
         zb = zb*dt

         eqd(kia,neqd1) = eqd(kia,neqd1)-kna*z
         eqd(kib,neqd1) = eqd(kib,neqd1)-knb*z
         eqd(kic,neqd1) = eqd(kic,neqd1)+knc*z
         eqd(kid,neqd1) = eqd(kid,neqd1)+knd*z

c..  jacobian
         eqd(kia,nn+kia) = eqd(kia,nn+kia)-kna*za
         eqd(kia,nn+kib) = eqd(kia,nn+kib)-kna*zb
         eqd(kib,nn+kia) = eqd(kib,nn+kia)-knb*za
         eqd(kib,nn+kib) = eqd(kib,nn+kib)-knb*zb
         eqd(kic,nn+kia) = eqd(kic,nn+kia)+knc*za
         eqd(kic,nn+kib) = eqd(kic,nn+kib)+knc*zb
         eqd(kid,nn+kia) = eqd(kid,nn+kia)+knd*za
         eqd(kid,nn+kib) = eqd(kid,nn+kib)+knd*zb

      enddo

      return
      end
