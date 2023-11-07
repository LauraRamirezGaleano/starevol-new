
      SUBROUTINE sort(n,arr)

************************************************************************
*                                                                      *
* Sort the array arr of length n using the Quicksort algorithm.        *
* From Press, W.H. et al. (1992). Numerical recipes in FORTRAN. The    *
* art of scientific computing. ISBN : 0-521-43064-X                    *
*                                                                      *
************************************************************************
*                                                                      *
* Note : If this algorithm is too slow you might want to consider      *
*        using the Heapsort algorithm described in the same source.    *
*                                                                      *
************************************************************************

      IMPLICIT NONE

************************************************************************
*                             DEFINITIONS                              *
************************************************************************
      
      INTEGER n,M,NSTACK
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      DOUBLE PRECISION arr(n)
      DOUBLE PRECISION a,temp

************************************************************************
*                         QUICKSORT ALGORITHM                          *
************************************************************************

      jstack=0
      l=1
      ir=n
 1    if(ir-l.lt.M)then 
         do j=l+1,ir
            a=arr(j)
            do i=j-1,l,-1
               if(arr(i).le.a)goto 2
               arr(i+1)=arr(i)
            enddo
            i=l-1
 2          arr(i+1)=a
         enddo
         if(jstack.eq.0)return
         ir=istack(jstack) 
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2 
         temp=arr(k)
         arr(k)=arr(l+1)
         arr(l+1)=temp
         if(arr(l).gt.arr(ir))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
         endif
         if(arr(l+1).gt.arr(ir))then
            temp=arr(l+1)
            arr(l+1)=arr(ir)
            arr(ir)=temp
         endif
         if(arr(l).gt.arr(l+1))then
            temp=arr(l)
            arr(l)=arr(l+1)
            arr(l+1)=temp
         endif
         i=l+1 
         j=ir
         a=arr(l+1) 
 3       continue 
         i=i+1
         if(arr(i).lt.a)goto 3
 4       continue
         j=j-1 
         if(arr(j).gt.a)goto 4
         if(j.lt.i)goto 5 
         temp=arr(i) 
         arr(i)=arr(j)
         arr(j)=temp
         goto 3 
 5       arr(l+1)=arr(j) 
         arr(j)=a
         jstack=jstack+2
         if(jstack.gt.NSTACK) pause 'NSTACK too small in sort'
         if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
      goto 1
      
************************************************************************
*                               END                                    *
************************************************************************

      END SUBROUTINE sort
