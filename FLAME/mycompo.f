      SUBROUTINE mycompo

*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *
      implicit none

      include 'evolpar.star'

      include 'evolcom.nuc'
      include 'evolcom.spec'
      include 'evolcom.teq'

      integer tab,itab,ntab,i,j,i1,elems,nelem,ielem
      double precision xspmin

      parameter (ntab=10,xspmin=1.d-1,nelem=5)
      dimension tab(ntab),elems(nelem)

      elems(1) = ih1
      elems(2) = ihe4
      elems(3) = ic12
      elems(4) = io16
      elems(5) = ine20
      
      do j = 1,nelem
         
         ielem = elems(j)

         do i = 1,ntab
            tab(i) = 0
         enddo

         itab = 0
         if (xsp(1,ielem).gt.xspmin) then
            itab = itab + 1
            tab(itab) = 1
         endif
         do i = 2,nmod-1
            i1 = i-1
            if (xsp(i,ielem).ge.xspmin.and.xsp(i1,ielem).lt.xspmin) then
               itab = itab + 1
               tab(itab) = i
            else if (xsp(i,ielem).lt.xspmin.and.xsp(i1,ielem).ge.xspmin) 
     &         then
               itab = itab + 1
               tab(itab) = i1
            endif         
            if (itab.ge.(ntab-1)) goto 10
         enddo
         if (xsp(nmod,ielem).gt.xspmin) then
            itab = itab + 1
            tab(itab) = nmod
         endif

   10    if (itab.gt.0) then
            write (*,100) elem(ielem),(tab(i),i=1,itab)
         endif
   
      enddo
   
      if (itab.gt.0) write (*,*)

  100 format (1x,a,': shells',5(' ',i4,'-',i4))
     
      return
      end
