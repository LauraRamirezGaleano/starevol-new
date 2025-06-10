
      PROGRAM BINDING

***********************************************************************
* Computation of nuclear binding energies                             *
* from analytic formula of Myers and Swiatecki, 1966                  *
***********************************************************************

      implicit none

      double precision c1,c2,c3,c4,a02,x,e,f
      double precision a1,a2,a3,a4
      double precision bvol,bsurf,bsym,rc,isym,e2
      double precision aa,nn,zz
      double precision eb,eba

      integer nnucl
      integer a,n,z
      integer i

      parameter (nnucl = 286)

      dimension a(nnucl),n(nnucl),z(nnucl)
      dimension eb(nnucl),eba(nnucl)

      open (unit = 10,file = 'binding.dat')
      open (unit = 20,file = 'binding.out')

      do i = 1,nnucl
         read (10,1000) a(i),z(i)
         n(i) = a(i)-z(i)
      enddo

      do i = 1,nnucl

         aa = dble(a(i))
         nn = dble(n(i))
         zz = dble(z(i))

         a1 = 16.9177d0
         a2 = 19.120d0
         a3 = 0.76278d0
         a4 = 101.777d0

         c1 = 15.677d0*(1.d0-1.79d0*((nn-zz)/aa)**2)
         c2 = 18.56d0*(1.d0-1.79d0*((nn-zz)/aa)**2)
         c3 = 0.717d0
         c4 = 1.21129d0
         a02 = 0.3645d0*aa**(-2.d0/3.d0)
         x = c3*zz*zz/(2.d0*c2*aa)
         e = 2.d0/5.d0*c2*aa**(2.d0/3.d0)*(1.d0-x)*a02
         f = 4.d0/105.d0*c2*aa**(2.d0/3.d0)*(1.d0+2.d0*x)*
     &        a02**(3.d0/2.d0)

         bvol = 15.56d0
         bsurf = 17.23d0
         bsym = 46.57d0
         rc = 1.24d0*aa**(1.d0/3.d0)
         isym = nn-zz
         e2 = 1.44d0

c        eb(i) = c1*aa-c2*aa**(2.d0/3.d0)-c3*(zz*zz*aa**(-1.d0/3.d0))+
c    &        c4*zz*zz/aa-(4.d0*e**3/(9.d0*f*f))+(8.d0*e**3/(27.d0*f*f))
c        eb(i) = a1*aa-a2*aa**(2.d0/3.d0)-a3*(zz*zz*aa**(-1.d0/3.d0))-
c    &        a4/(4.d0*aa)*(aa-2.d0*zz)**2
         eb(i) = bvol*aa-bsurf*aa**(2.d0/3.d0)-0.5d0*bsym*isym*isym
     &        /aa-3.d0/5.d0*zz*zz*e2/rc
         eba(i) = eb(i)/aa

      enddo

      do i = 1,nnucl
         write (20,2000) i,a(i),z(i),n(i),eba(i)
      enddo

 1000 format (9x,i3,1x,i2)
 2000 format (1x,i3,1x,3(1x,i3),2x,f8.5)

      end
