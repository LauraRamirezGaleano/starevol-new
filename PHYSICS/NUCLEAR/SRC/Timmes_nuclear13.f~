      program timenet
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
      include 'burn_common.dek'
      include 'network.dek'


c..this program exercises the aprox13 network

c..declare
      character*40     hfile
      logical          lfile
      integer          i,j,nok,nbad,iind(13)
      double precision tstep,conserv,
     1                 tin,din,ein,xin(13),tout,dout,eout,xout(13),
     2                 ydum(13),abar,zbar



c..initialize the network 
       call init_aprox13


c...write the time evolution and set the output file name
       lfile = .true.
       hfile = 'aprox13_a'


c..set the burning mode
      one_step    = .false.
      hydrostatic = .true.
      expansion   = .false.
      self_heat   = .false.
      bbang       = .false.



c..set the initial temperature and density
      tin   = 6.0e8
      din   = 1.0e3



c..set the initial composition
      do i=1,ionmax
       xin(i) = 1.0d-30
      enddo

      xin(ihe4) = 1.0d0



c..set the ending time
      tstep = 1.0e6





c..set the peak expansion temperature if needed
c..psi =  1 is an adiabatic expansion
c..psi = -1 is an adiabatic implosion
      if (expansion) then
       psi       = 1.0d0
       den0      = din
       temp0     = tin
       temp_stop = 1.0d7
       if ( (psi .eq. 1.0  .and. temp_stop .ge. tin)  .or.
     1      (psi .eq. -1.0 .and. temp_stop .le. tin)) 
     2    stop 'bad adiabatic temp_stop in routine burner'

      else
       psi       = 0.0d0
       temp_stop = 1.0d30
      end if




c..screening and table options
       screen_on   = 1
       use_tables  = 0



c..halt the integration when the mass fraction 
c..of name_stop drops below xmass_stop 

      name_stop = 'he4 '
      xmass_stop = -1.0d30


c..be sure the isotope is in the network
      do i=1,ionmax
       if (ionam(i) .eq. name_stop) then
        id_stop = i
        goto 11
       end if
      enddo
      stop 'name_stop not in network' 
 11   continue




c..a message
        write(6,*) 
        write(6,*) 'starting integration'
        write(6,*) 




c..get the internal energy corresponding to the
c..initial temperature,density, and composition

      call azbar(xin,aion,zion,ionmax,
     1           ydum,abar,zbar)

c..call an eos 
       temp_row(1) = tin
       den_row(1)  = din
       abar_row(1) = abar
       zbar_row(1) = zbar
       jlo_eos = 1
       jhi_eos = 1

       call helmeos

       ein = etot_row(1)   



c..burn it
        call burner(tstep,tin,din,ein,xin,tout,dout,eout,xout,
     1              conserv,nok,nbad,
     2              lfile,hfile)




c..output 
       write(6,*) ' '
       write(6,04) netname,
     1             ' tin =',tin,' din =',din,' ein =',ein,
     2             ' tout=',tout,' dout=',dout,' eout=',eout,
     3             ' sum =',conserv,' nbad=',nbad,' nok=',nok
04     format(1x,a,':',/,
     1        1x,3(a,1pe20.12),/,
     2        1x,3(a,1pe20.12),/,
     3        1x,a,1pe11.3,2(a,i5))
       write(6,*) ' '



c..write out the biggest mass fractions
       call indexx(ionmax,xout,iind)
       j = min(20,ionmax)
       write(6,09) 'top ',j,' mass fractions:'
09     format(1x,a,i2,a)
       j = max(ionmax-19,1)
       write(6,02) (ionam(iind(i)),xout(iind(i)), i=ionmax,j,-1)
02    format(1x,a4,'=',1pe10.3,' ',a4,'=',1pe10.3,' ',
     1          a4,'=',1pe10.3,' ',a4,'=',1pe10.3,' ',
     2          a4,'=',1pe10.3)


      end      










      subroutine burner(tstep,tin,din,ein,xin,tout,dout,eout,xout,
     1                  conserv,nok,nbad,
     2                  lfile,hfile)
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
      include 'burn_common.dek'
      include 'network.dek'

c..given time step tstep, temperature tin, density din, thermal
c..energy ein, and the composition xin, this routine returns the 
c..burned composition xout, final temperature tout, final density dout, 
c..and the final thermal energy eout.

c..declare the pass
      character*40     hfile
      logical          lfile
      integer          nok,nbad
      double precision tstep,tin,din,ein,xin(1),
     1                 tout,dout,eout,xout(1),conserv


c..local variables
      integer          i
      double precision xcons,yex,ydum(neqs),
     1                 dydt_dum(neqs),xdum(neqs),abar,zbar

      double precision conv
      parameter        (conv = ev2erg*1.0d6*avo)

      external         aprox13,daprox13
      external         stifbs_leqs




c..tdim sets the size of the storage arrays
c..keep tdim small if evolutions are not being written out
c..and make tdim larger if evolution are being written out
      integer          tdim,kount,iprint,kstore
      parameter        (tdim=1000, iprint=0, kstore=tdim)



c..for the integration driver
      double precision beg,stptry,stpmin,tend,dxsav,
     1                 ys2(neqs),ttime(tdim),elem(neqs,tdim)
      parameter        (dxsav   = 0.0d0)


c..for easy tolerance changes 
      double precision odescal,tol
      parameter        (tol     = 1.0d-6, 
     1                  odescal = 1.0d-6)


c..for writing an evolution history 
      character*80     string
      integer          j,lop,rem,ilop,jrem,kb,ke,k
      double precision sum




c..popular format statements
01    format(1x,'*',t9,a,t20,a,t31,a,t42,a,t53,a,t64,a,   
     1              t75,a,t86,a,t97,a,t108,a,t119,a)
02    format(1x,a4,' =',1pe11.3,'   ',a4,' =',1pe11.3,'   ',
     1          a4,' =',1pe11.3,'   ',a4,' =',1pe11.3)
03    format(a30,i8,a4)
04    format(1x,1p15e11.3)
05    format(1x,i4,1p11e11.3)
06    format(1x,'* ',a,': ',2(a,i5))
07    format(1x,'* ',a,4(a,1pe11.3))



 

c..load the mass fractions
      do i=1,ionmax
       xmass(i) = xin(i)
      enddo


c..get abar, zbar and a few other composition variables
      call azbar(xmass,aion,zion,ionmax,
     1           ymass,abar,zbar)


c..stuff the initial conditions into ys2
      do i=1,ionmax
       ys2(i) = ymass(i)
      enddo


c..energy, temperature, density
      ys2(iener) = ein
      ys2(itemp) = tin
      ys2(iden)  = din






c..single step (tend=tstep), hydrostatic, or expansion ending times
c..the variable tstep has two meanings here. tstep in single step mode
c..is the size of the time step to try. tstep in hydrostatic or expansion
c..mode is the ending integration time. the integration driver really
c..gets some exercise if tstep is large in single step mode.

      beg = 0.0d0
      tend = tstep
      if (one_step) then
       stptry = tstep 
       stpmin = tstep * 1.0d-20
      else
       stptry = 1.0d-16
       stpmin = stptry * 1.0d-12
      end if


       


c..zero the output array
      do j=1,tdim
       do i=1,neqs
        elem(i,j) = 0.0d0
       enddo
      enddo 


c..integrate the aprox13 network
       call netint(beg,stptry,stpmin,tend,ys2,   
     1             tol,dxsav,kstore,   
     2             ttime,elem,tdim,neqs,tdim,neqs,
     3             nok,nbad,kount,odescal,iprint,
     4             aprox13,daprox13,stifbs_leqs)



c..set the output composition
      do i=1,ionmax
       xout(i) = ys2(i) * aion(i)
      enddo



c..output temperature, density, and thermal energy 
      tout = elem(itemp,kount)
      dout = elem(iden,kount)
      eout = elem(iener,kount)



c..set the mass non-conservation
      conserv = 0.0d0
      do i=1,ionmax
       conserv = conserv + xout(i)
      enddo
      conserv = 1.0d0 - conserv





c..if the evolution history is to be written out
      if (lfile) then


c..write the above quantities
       write(string,03) hfile,0,'.dat'
       call sqeeze(string)
       open (unit=2, file=string, status='unknown')

       if (one_step) then
        write(2,07) 'one_step:','  btemp=',btemp,' bden=',bden

       else if (hydrostatic) then
        write(2,07) 'hydrostatic:','  btemp=',btemp,' bden=',bden

       else if (expansion) then
        write(2,07) 'expansion:','  temp0=',temp0,' den0=',den0,
     1              ' temp_stop=',temp_stop

       else if (self_heat) then
        write(2,07) 'self_heat:','   temp0=',elem(itemp,1),
     1                           '   den0=',elem(iden,1)

       end if


       write(2,06) netname,'  nbad=',nbad,'  nok=',nok

       write(2,01) 'time','temp','den','ener','sdot','sneut',
     1             's-snu','ye','1-sum'



c..loop through the output points
       do i=1,kount


c..form the mass fractions
        do j=1,ionmax
         xdum(j) = min(1.0d0,max(elem(j,i) * aion(j),1.0d-30))
        enddo


c..mass conservation
        sum = 0.0d0
        do j=1,ionmax
         sum = sum + xdum(j)
        enddo
        sum = 1.0d0 - sum
        xcons = sum


c..get ye using normalized mass fractions
        sum = 0.0d0
        do j=1,ionmax
         sum = sum + xdum(j)
        enddo
        sum = 1.0d0/sum
        do j=1,ionmax
         xdum(j) = min(1.0d0,max(sum * xdum(j),1.0d-30))
        enddo

c..get abar, zbar and a few other composition variables
        call azbar(xdum,aion,zion,ionmax,
     1             ydum,abar,zbar)

        yex = zbar/abar


c..finish filling the dummy array
        do j=1,ionmax
         ydum(j) = elem(j,i)
        enddo
        ydum(iener) = elem(iener,i)
        ydum(itemp) = elem(itemp,i)  
        ydum(iden)  = elem(iden,i)


c..get the right hand sides, exact energy generation rate and so on 
        call aprox13(ttime(i),ydum,dydt_dum)



c..and write what we found

        write(2,05) i,ttime(i),elem(itemp,i),elem(iden,i),
     1              elem(iener,i),sdot,sneut,dydt_dum(iener),
     2              yex,xcons


c..end of kount loop 
        enddo
        close(unit=2)





c..write out the isotopic mass fractions in blocks of 8
c..lop is how many groups of 8 exist; jrem is the remainder
       lop  = ionmax/8
       jrem  = ionmax - 8*lop
       do ilop = 1,lop+1
        kb = 1 + 8*(ilop-1)
        ke = 8 + 8*(ilop-1)
        if (ilop .eq. lop+1  .and. jrem .eq. 0) goto 60
        if (ilop .eq. lop+1) ke = ionmax

c..open the output file
        write(string,03) hfile,ilop,'.dat'
        call sqeeze(string)
        open (unit=2, file=string, status='unknown')

        if (one_step) then
         write(2,07) 'one_step:','  btemp=',btemp,' bden=',bden
 
        else if (hydrostatic) then
         write(2,07) 'hydrostatic:','  btemp=',btemp,' bden=',bden
 
        else if (expansion) then
         write(2,07) 'expansion:','  temp0=',temp0,' den0=',den0,
     1                            ' temp_stop=',temp_stop

        else if (self_heat) then
         write(2,07) 'self_heat:','   temp0=',elem(itemp,1),
     1                            '   den0=',elem(iden,1)

        end if


        write(2,06) netname,'  nbad=',nbad,'  nok=',nok
        write(2,01) 'time',(ionam(k), k=kb,ke)
        do i=1,kount
          write(2,04) ttime(i),(elem(k,i)*aion(k), k=kb,ke)
        enddo

        close(unit=2)
60      continue
       enddo

c..end of lfile if
      end if


      return
      end







      subroutine netint(start,stptry,stpmin,stopp,bc, 
     1                  eps,dxsav,kkmax,  
     2                  xrk,yrk,xphys,yphys,xlogi,ylogi, 
     3                  nok,nbad,kount,odescal,iprint,
     4                  derivs,jakob,steper)   
      include 'implno.dek'
      include 'network.dek'


c..ode integrator for stiff odes
c..tuned for nnuclear reacton networks 

c..input:
c..start   = beginning integration point
c..stptry  = suggested first step size
c..stpmin  = minimum allowable step size
c..stopp   = ending integration point
c..bc      = initial conditions, array of physical dimension yphys
c..eps     = desired fraction error during the integration
c..dxsav   = incremental vale of indep variable at which to store the solution
c..          if zero, solution is stored at every intermediate point
c..          if not zero, solution is done and saved at every dxsav point
c..kkmax   = maximum number of solution points to store, kkmax < xphys
c..odescal = error scaling factor 
c..iprint  = integer to determines if the solution is printed as it evolves
c..derivs  = name of the routine that contains the odes
c..jakob   = name of the routine that contains the jacobian of the odes
c..steper  = name of the routine that will take a single step

c..output:
c..nok    = number of succesful steps taken
c..nbad   = number of bad steps taken, bad but retried and then succesful
c..kount  = total number of steps stored in arrays xrk and yrk
c..xrk    = the independent variable solution 
c..         array of physical dimension xphys, logical dimension xlogi
c..yrk    = the dependent variable solution 
c..         array of physical dimension (yphys,xphys) with 
c..                  logical  dimension (ylogi,xlogi)



c..declare the pass
      external         derivs,jakob,steper 
      integer          xphys,yphys,xlogi,ylogi,kkmax,nok,nbad,kount,
     1                 iprint
      double precision start,stptry,stpmin,stopp,bc(yphys),eps,
     1                 dxsav,xrk(xphys),yrk(yphys,xphys),
     2                 odescal


c..local variables
      integer          nmax,stpmax,i,j,nstp 
      parameter        (nmax = 500, stpmax=200000)   
      double precision yscal(nmax),y(nmax),dydx(nmax),   
     1                 x,xsav,h,hdid,hnext,sum,xdum(nmax),tiny
      parameter        (tiny=1.0d-15)



c..for smooth plot timesteps
      double precision ratio,xfloor,ychangemax,ynew,yold,yy,dtx





c..here are the format statements for printouts as we integrate 
100   format(1x,i6,1p4e10.2) 
101   format(1x,1p12e10.2) 


c..initialize    
      if (ylogi  .gt. yphys) stop 'ylogi > yphys in routine netint' 
      if (yphys  .gt. nmax)  stop 'yphys > nmax in routine netint' 
      if (kkmax  .gt. xphys) stop 'kkmax > xphys in routine netint'
      x     = start    
      h     = sign(stptry,stopp-start)  
      nok   = 0  
      nbad  = 0 
      kount = 0    


c..store the first step  
      do i=1,ylogi 
       y(i) = bc(i)   
      enddo
      xsav = x - 2.0d0 * dxsav 

c..take at most stpmax steps 
      do nstp=1,stpmax 


c..positive definite abundance fractions
       do i=1,ionmax
        y(i) = min(1.0d0, max(y(i),1.0d-30))
       enddo



c..get the right hand sides
       call derivs(x,y,dydx) 


c..scaling vector used to monitor accuracy 
       do i=1,ylogi 
        yscal(i) = max(odescal,abs(y(i))) 
       enddo


c..store intermediate results    
       if (kkmax .gt. 0) then 
        if ( (abs(dxsav) - abs(x-xsav)) .le. tiny) then  
         if ( kount .lt. (kkmax-1) ) then   
          kount         = kount+1   
          xrk(kount)    = x    
          do i=1,ylogi  
           yrk(i,kount) = y(i) 
          enddo
          if (iprint .eq. 1) then 
c           write(6,100) kount,xrk(kount),yrk(1,kount) 
          write(6,100) kount,xrk(kount),yrk(itemp,kount),
     1                 yrk(iden,kount),yrk(iener,kount) 
c           write(6,100) kount,xrk(kount),yrk(iener,kount) 
c           write(6,101) (yrk(j,kount), j=1,ylogi) 
          end if 
          call flush(6) 
          xsav=x  
         end if 
        end if   
       end if 


c..if the step can overshoot the stop point or the dxsav increment then cut it 
       if ((x+h-stopp)*(x+h-start) .gt. 0.0d0) h = stopp - x   
       if (dxsav .ne. 0.0d0 .and. h.gt.(xsav-x+dxsav)) h=xsav+dxsav-x 


c..do an integration step 
       call steper(y,dydx,ylogi,x,h,eps,yscal,hdid,hnext,
     1             derivs,jakob,nstp)

       if (hdid.eq.h) then 
        nok = nok+1    
       else  
        nbad = nbad+1  
       end if 



c..this is the normal exit point, save the final step 
       if ( (nstp .eq. stpmax)                         .or. 
     1      (x-stopp)*(stopp-start).ge. 0.0d0          .or.
     2      (psi .eq.  1.0 .and. y(itemp) .lt. temp_stop)  .or.
     3      (psi .eq. -1.0 .and. y(itemp) .gt. temp_stop)  .or.
     4      (y(id_stop)*aion(id_stop) .lt. xmass_stop)   ) then 

        do i=1,ylogi   
         bc(i) = y(i)  
        enddo
        if (kkmax.ne.0) then    
         kount         = kount+1   
         xrk(kount)    = x    
         do i=1,ylogi  
          yrk(i,kount) = y(i)  
         enddo
         if (iprint .eq. 1) then 
c          write(6,100) kount,xrk(kount),yrk(1,kount) 
           write(6,100) kount,xrk(kount),yrk(itemp,kount),
     1                  yrk(iden,kount),yrk(iener,kount) 
c          write(6,101) (yrk(j,kount), j=1,ylogi) 
         end if 
         call flush(6) 
        end if 
        return   
       end if 


c..set the step size for the next iteration; stay above stpmin 

       h = hnext 

       if (abs(hnext).lt.stpmin) stop 'hnext < stpmin in netint' 



c..back for another iteration or death 
      enddo
      stop 'more than stpmax steps required in netint'  
      end 








      subroutine stifbs_leqs(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,
     1                        derivs,jakob,nstp)
      include 'implno.dek'
      include 'dense_matrix.dek'


c..for dense analytic jacobians, lu decomposition linear algebra 

c..semi-implicit extrapolation step for integrating stiff ode's with monitoring 
c..of local truncation error to adjust stepsize. inputs are the dependent  
c..variable vector y(1:nv) and its derivative dydx(1:nv) at the starting of the 
c..independent variable x. also input are the stepsize to be attempted htry, 
c..the required accuracy eps, and the vector yscal against which the error is 
c..scaled. on output, y and x are replaced by their new values, hdid is the  
c..stepsize actually accomplished, and hnext is the estimated next stepsize. 
c..dervs is a user supplied function that computes the right hand side of 
c..the equations.
c..
c..declare  
      external         derivs,jakob
      logical          first,reduct
      integer          nv,nmax,kmaxx,imax
      parameter        (nmax  = iodemax, 
     1                  kmaxx = 7, 
     2                  imax  = kmaxx+1)
      integer          i,iq,j,k,kk,km,kmax,kopt,nvold,nseq(imax),ifirst,
     1                 iat
      double precision y(nv),dydx(nv),x,htry,eps,yscal(nv),hdid,hnext,
     1                 eps1,epsold,errmax,fact,h,red,scale,work,wrkmin,
     2                 xest,xnew,a(imax),alf(kmaxx,kmaxx),err(kmaxx),
     3                 yerr(nmax),ysav(nmax),yseq(nmax),safe1,safe2,
     4                 redmax,redmin,tiny,scalmx
      parameter        (safe1 = 0.25d0, safe2 = 0.7d0, redmax=1.0d-5,
     1                  redmin = 0.7d0, tiny = 1.0d-30, 
     2                  scalmx = 0.1d0)

c..for jacobian pictures
      character*20     string
      integer          nstp

c..assume that the independent variable is not explicit in the odes
      data             first/.true./, epsold/-1.0d0/, nvold/-1/
      data             nseq /2, 6, 10, 14, 22, 34, 50, 70/
      data             ifirst/0/   



c..a new tolerance or a new number, so reinitialize
      if (eps .ne. epsold  .or.  nv .ne. nvold) then
       hnext = -1.0e29
       xnew  = -1.0e29
       eps1  = safe1 * eps


c..compute the work coefficients a_k
       a(1)  = nseq(1) + 1
       do k=1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
       enddo

c..compute alf(k,q)
       do iq=2,kmaxx
        do k=1,iq-1
         alf(k,iq) = eps1**((a(k+1) - a(iq+1)) /
     1               ((a(iq+1) - a(1) + 1.0d0) * (2*k + 1)))
        enddo
       enddo
       epsold = eps
       nvold  = nv

c..add cost of jacobians to work coefficients
       a(1)   = nv + a(1)
       do k=1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
       enddo

c..determine optimal row number for convergence
       do kopt=2,kmaxx-1
        if (a(kopt+1) .gt. a(kopt)*alf(kopt-1,kopt)) go to 01
       enddo
01     kmax = kopt
      end if

c..save the starting values 
      h    = htry
      do i=1,nv  
       ysav(i)  = y(i)
      enddo

c..get the dense jacobian in dens_dfdy
      call jakob(x,y,dens_dfdy,nv,iodemax) 



c..a new stepsize or a new integration, re-establish the order window
      if (h .ne. hnext  .or.  x .ne. xnew) then
       first = .true.
       kopt = kmax
      end if
      reduct = .false.

c..evaluate the sequence of semi implicit midpoint rules
02    do 18 k=1,kmax

       xnew = x + h

       if (xnew .eq. x) stop 'stepsize too small in routine stiffbs'

       call simpr_leqs(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
       xest = (h/nseq(k))**2 
       call net_pzextr(k,xest,yseq,y,yerr,nv) 


c..compute normalized error estimate
       if (k .ne. 1) then
        errmax = tiny
        do i=1,nv

         errmax = max(errmax,abs(yerr(i)/yscal(i)))

        enddo
        errmax   = errmax/eps   


        km = k - 1
        err(km) = (errmax/safe1)**(1.0d0/(2*km+1))
       end if

c..in order window
       if (k .ne. 1  .and. (k .ge. kopt-1  .or. first)) then


c..converged
        if (errmax .lt. 1.0) go to 04


c..possible step size reductions
        if (k .eq. kmax  .or.  k .eq. kopt + 1) then
         red = safe2/err(km)
         go to 03
        else if (k .eq. kopt) then
         if (alf(kopt-1,kopt) .lt. err(km)) then
          red = 1.0d0/err(km)
          go to 03
         end if
        else if (kopt .eq. kmax) then
         if (alf(km,kmax-1) .lt. err(km)) then
          red = alf(km,kmax-1) * safe2/err(km)
          go to 03
         end if
        else if (alf(km,kopt) .lt. err(km)) then
         red = alf(km,kopt-1)/err(km)
         go to 03
        end if
       end if
18    continue


c..reduce stepsize by at least redmin and at most redmax
03    red    = min(red,redmin)
      red    = max(red,redmax)
      h      = h * red
      reduct = .true.
      go to 2


c..successful step; get optimal row for convergence and corresponding stepsize
04    x = xnew
      hdid = h
      first = .false.
      wrkmin = 1.0e35
      do kk=1,km
       fact = max(err(kk),scalmx)
       work = fact * a(kk+1)
       if (work .lt. wrkmin) then
        scale  = fact
        wrkmin = work
        kopt   = kk + 1
       end if
      enddo
c..
c..check for possible order increase, but not if stepsize was just reduced
      hnext = h/scale
      if (kopt .ge. k  .and.  kopt .ne. kmax  .and.  .not.reduct) then
       fact = max(scale/alf(kopt-1,kopt),scalmx)
       if (a(kopt+1)*fact .le. wrkmin) then
        hnext = h/fact
        kopt = kopt + 1
       end if
      end if
      return
      end







      subroutine simpr_leqs(y,dydx,n,xs,htot,nstep,yout,derivs) 
      include 'implno.dek'
      include 'dense_matrix.dek'

c..an implicit midpoint stepper, for lu decomp dense linear algebra. 

c..declare 
      external         derivs 
      integer          nmax,n,nstep,nmaxx
      parameter        (nmaxx=iodemax) 
      integer          i,j,nn
      double precision y(n),dydx(n),xs,htot, 
     1                 yout(n),h,x,del(nmaxx),ytemp(nmaxx)

c..for the leqs linear algebra
      double precision dsav(iodemax,iodemax)


c..stepsize this trip, and make the a matrix 
      h = htot/nstep 
      do j=1,n
       do i=1,n
        dmat(i,j) = -h * dens_dfdy(i,j) 
       enddo
      enddo
      do i=1,n
       dmat(i,i) = 1.0d0 + dmat(i,i) 
      end do

      do j=1,n
       do i=1,n
        dsav(i,j) = dmat(i,j)
       enddo
      enddo


c..use yout as temporary storage; the first step 
      do i=1,n 
       yout(i) = h * dydx(i)
      enddo

      do j=1,n
       do i=1,n
        dmat(i,j) = dsav(i,j)
       enddo
      enddo
      call leqs(dmat,yout,n,iodemax) 

      do 14 i=1,n 
       del(i)   = yout(i) 
       ytemp(i) = y(i) + del(i) 
14    continue 
      x = xs + h 
      call derivs(x,ytemp,yout) 
c.. 
c..use yout as temporary storage; general step 
      do 17 nn=2,nstep 
       do 15 i=1,n 
        yout(i) = h*yout(i) - del(i) 
15     continue

       do j=1,n
        do i=1,n
         dmat(i,j) = dsav(i,j)
        enddo
       enddo
       call leqs(dmat,yout,n,iodemax) 

       do 16 i=1,n 
        del(i)   = del(i) + 2.0d0 * yout(i) 
        ytemp(i) = ytemp(i) + del(i) 
16     continue 
       x = x + h 
       call derivs(x,ytemp,yout) 
17    continue 
c.. 
c..take the last step 
      do 18 i=1,n 
       yout(i) = h * yout(i) - del(i)  
18    continue

      do j=1,n
       do i=1,n
        dmat(i,j) = dsav(i,j)
       enddo
      enddo
      call leqs(dmat,yout,n,iodemax) 

      do 19 i=1,n 
       yout(i) = ytemp(i) + yout(i) 
19    continue       
      return 
      end 









      subroutine net_pzextr(iest,xest,yest,yz,dy,nv) 
      include 'implno.dek'

c..use polynomial extrapolation to evaluate nv functions at x=0 by fitting 
c..a polynomial to a sequence of estimates with progressively smaller values 
c..x=xest, and corresponding function vectors yest(1:nv). the call is number  
c..iest in the sequence of calls. extrapolated function values are output as  
c..yz(1:nv), and their estimated error is output as dy(1:nv) 


c..declare 
      integer          iest,nv,j,k1,nmax,imax 
      parameter        (nmax=500, imax=13) 
      double precision xest,dy(nv),yest(nv),yz(nv),delta,f1,f2,q, 
     1                 d(nmax),qcol(nmax,imax),x(imax) 

c..sanity checks
      if (iest .gt. imax) stop 'iest > imax in net_pzextr'
      if (nv .gt. nmax) stop 'nv > nmax in net_pzextr'


c..save current independent variables 
      x(iest) = xest 
      do j=1,nv 
       dy(j) = yest(j) 
       yz(j) = yest(j) 
      enddo

c..store first estimate in first column 
      if (iest .eq. 1) then 
       do j=1,nv 
        qcol(j,1) = yest(j) 
       enddo
      else 
       do j=1,nv 
        d(j) = yest(j) 
       enddo
       do k1=1,iest-1 
        delta = 1.0d0/(x(iest-k1) - xest) 
        f1    = xest * delta 
        f2    = x(iest-k1) * delta 

c..propagate tableu 1 diagonal more 
        do j=1,nv 
         q          = qcol(j,k1) 
         qcol(j,k1) = dy(j) 
         delta      = d(j) - q 
         dy(j)      = f1*delta 
         d(j)       = f2*delta 
         yz(j)      = yz(j) + dy(j) 
        enddo
       enddo
       do j=1,nv 
        qcol(j,iest) = dy(j) 
       enddo
      end if 
      return 
      end 






      subroutine leqs(a,b,n,np) 
      include 'implno.dek' 

c..solves a linear system of equations a x = b via gauss jordan elimination 
c..with no pivoting. a is destroyed on exit. on input b contains the 
c..right hand side, on exit it is the solution. 

c..declare 
      integer          n,np,n1,i,j,k,l,imax,jj 
      double precision a(np,np),b(np),r,c 


c..for each column 
      n1 = n-1 
      do i=1,n 

c..find maximum element in each row 
       r = abs(a(i,1)) 
       do j=2,n 
        c = abs(a(i,j)) 
        if (r.lt.c) r=c 
       enddo

c..divide that row and the right hand side by the maximum element 
       do j=1,n 
        a(i,j) = a(i,j)/r 
       enddo
       b(i) = b(i)/r 
      enddo


c..for each column, do the elimination; bail if the element is zero 
      do j=1,n1 
       l = j + 1 
       do i=l,n 
        r = -a(i,j) 
        if (r.eq.0.0) go to 50 
        r = r/a(j,j) 
        do k=l,n 
         a(i,k) = a(i,k) + r*a(j,k) 
        enddo
        b(i) = b(i) + r*b(j) 
50      continue 
       enddo
      enddo


c..the matrix is now in triangular form; back sub 
      b(n) = b(n)/a(n,n) 
      do l=1,n1 
       i    = n-l 
       r    = 0.0d0 
       imax = i + 1 
       do j=imax,n 
        jj = i + n + 1 - j 
        r  = r + a(i,jj)*b(jj) 
       enddo
       b(i) = (b(i) - r) / a(i,i) 
      enddo
      return 
      end 





cxxx




c..this file contains aprox13 network 

c..routine aprox13 sets up the odes 
c..routine faprox13 gets the flows from aprox13
c..routine daprox13 sets up the dense aprox13 jacobian
c..routine aprox13rat generates the reaction rates for routine aprox13
c..routine aprox13tab generates the raw rates using table interpolation
c..routine screen_aprox13 applies screening corrections to the raw rates
c..routine init_aprox13 initializes the aprox13 network




      subroutine aprox13(tt,y,dydt)   
      include 'implno.dek'
      include 'const.dek'
      include 'burn_common.dek' 
      include 'network.dek'
      include 'vector_eos.dek'

c..this routine sets up the system of ode's for the aprox13 nuclear reactions.
c..this is an alpha chain + heavy ion network with (a,p)(p,g) links
c..
c..isotopes: he4,  c12,  o16,  ne20, mg24, si28, s32,
c..          ar36, ca40, ti44, cr48, fe52, ni56



c..declare the pass
      double precision tt,y(1),dydt(1)


c..local variables
      integer          i
      double precision r1,s1,t1,u1,v1,w1,x1,y1,denom,
     1                 enuc,taud,taut,z,sum1,sum2,
     2                 zbarxx,ytot1,abar,zbar,
     3                 snudt,snudd,snuda,snudz

      double precision conv,sixth
      parameter        (conv = ev2erg*1.0d6*avo,
     1                  sixth = 1.0d0/6.0d0)



c..positive definite mass fractions
      do i=1,ionmax
       y(i) = min(1.0d0,max(y(i),1.0d-30))
      enddo


c..positive definite temperatures and densities
      y(itemp) = min(1.0d11,max(y(itemp),1.0d4))
      y(iden)  = min(1.0d11,max(y(iden),1.0d-10))



c..set the common block temperature and density
      btemp = y(itemp)
      bden  = y(iden)


c..generate abar and zbar for this composition
      zbarxx  = 0.0d0
      ytot1   = 0.0d0
      do i=1,ionmax
       ytot1    = ytot1 + y(i)
       zbarxx   = zbarxx + zion(i) * y(i)
      enddo
      abar  = 1.0d0/ytot1
      zbar  = zbarxx * abar



c..get the reaction rates
      if (use_tables .eq. 1) then
       call aprox13tab
      else
       call aprox13rat
      end if



c..do the screening here because the corrections depend on the composition

      call screen_aprox13(y)



c..branching ratios for various dummy proton links
      r1     = 0.0d0
      denom  = ratdum(iralpa) + ratdum(iralpg)
      if (denom .ne. 0.0) r1 = ratdum(iralpa)/denom

      s1     = 0.0d0
      denom  = ratdum(irppa) + ratdum(irppg)
      if (denom .ne. 0.0) s1 = ratdum(irppa)/denom

      t1     = 0.0d0
      denom  = ratdum(irclpa) + ratdum(irclpg)
      if (denom .ne. 0.0) t1 = ratdum(irclpa)/denom

      u1     = 0.0d0
      denom  = ratdum(irkpa) + ratdum(irkpg)
      if (denom .ne. 0.0) u1 = ratdum(irkpa)/denom

      v1     = 0.0d0
      denom  = ratdum(irscpa) + ratdum(irscpg)
      if (denom .ne. 0.0) v1 = ratdum(irscpa)/denom

      w1     = 0.0d0
      denom  = ratdum(irvpa) + ratdum(irvpg)
      if (denom .ne. 0.0) w1 = ratdum(irvpa)/denom

      x1     = 0.0d0
      denom  = ratdum(irmnpa) + ratdum(irmnpg)
      if (denom .ne. 0.0) x1 = ratdum(irmnpa)/denom

      y1     = 0.0d0
      denom  = ratdum(ircopa) + ratdum(ircopg)
      if (denom .ne. 0.0) y1 = ratdum(ircopa)/denom





c..set up the system of odes: 
c..he4 reactions
c..heavy ion reactions
      dydt(ihe4) =  0.5d0 * y(ic12) * y(ic12) * ratdum(ir1212)      
     1            + 0.5d0 * y(ic12) * y(io16) * ratdum(ir1216)
     2            + 0.56d0 * 0.5d0 * y(io16) * y(io16) * ratdum(ir1616) 
     3            + 0.34d0 * s1 * 0.5d0 * y(io16)*y(io16)*ratdum(ir1616)

c..(a,g) and (g,a) reactions
      dydt(ihe4) =  dydt(ihe4)
     1            - 0.5d0 * y(ihe4) * y(ihe4) * y(ihe4) * ratdum(ir3a) 
     2            + 3.0d0 * y(ic12) * ratdum(irg3a)       
     3            - y(ihe4)  * y(ic12) * ratdum(ircag)       
     4            + y(io16)  * ratdum(iroga)
     5            - y(ihe4)  * y(io16) * ratdum(iroag)       
     6            + y(ine20) * ratdum(irnega)            
     7            - y(ihe4)  * y(ine20) * ratdum(irneag) 
     8            + y(img24) * ratdum(irmgga)
     9            - y(ihe4)  * y(img24)* ratdum(irmgag)
     &            + y(isi28) * ratdum(irsiga) 
     1            - y(ihe4)  * y(isi28)*ratdum(irsiag)
     2            + y(is32)  * ratdum(irsga) 

      dydt(ihe4) =  dydt(ihe4)
     1            - y(ihe4)  * y(is32) * ratdum(irsag)
     2            + y(iar36) * ratdum(irarga) 
     3            - y(ihe4)  * y(iar36)*ratdum(irarag)
     4            + y(ica40) * ratdum(ircaga) 
     5            - y(ihe4)  * y(ica40)*ratdum(ircaag)
     6            + y(iti44) * ratdum(irtiga) 
     7            - y(ihe4)  * y(iti44)*ratdum(irtiag)
     8            + y(icr48) * ratdum(ircrga) 
     9            - y(ihe4)  * y(icr48)*ratdum(ircrag)
     &            + y(ife52) * ratdum(irfega) 
     1            - y(ihe4)  * y(ife52) * ratdum(irfeag)
     2            + y(ini56) * ratdum(irniga)


c..(a,p)(p,g) and (g,p)(p,a) reactions
      dydt(ihe4) =  dydt(ihe4)
     1            - y(ihe4)  * y(img24) * ratdum(irmgap) * (1.0d0-r1)
     2            + y(isi28) * ratdum(irsigp) * r1
     3            - y(ihe4)  * y(isi28) * ratdum(irsiap) * (1.0d0-s1)
     4            + y(is32)  * ratdum(irsgp) * s1
     5            - y(ihe4)  * y(is32) * ratdum(irsap) * (1.0d0-t1)   
     6            + y(iar36) * ratdum(irargp) * t1 
     7            - y(ihe4)  * y(iar36) * ratdum(irarap) * (1.0d0-u1)
     8            + y(ica40) * ratdum(ircagp) * u1
     9            - y(ihe4)  * y(ica40) * ratdum(ircaap) * (1.0d0-v1)   
     &            + y(iti44) * ratdum(irtigp) * v1
     1            - y(ihe4)  * y(iti44) * ratdum(irtiap) * (1.0d0-w1)
     2            + y(icr48) * ratdum(ircrgp) * w1
     3            - y(ihe4)  * y(icr48) * ratdum(ircrap) * (1.0d0-x1)   
     4            + y(ife52) * ratdum(irfegp) * x1 
     5            - y(ihe4)  * y(ife52) * ratdum(irfeap) * (1.0d0-y1)
     6            + y(ini56) * ratdum(irnigp) * y1


c..c12 reactions   
      dydt(ic12) = -y(ic12) * y(ic12) * ratdum(ir1212) 
     1            - y(ic12) * y(io16) * ratdum(ir1216)       
     2            + sixth * y(ihe4) * y(ihe4) * y(ihe4) * ratdum(ir3a)
     3            - y(ic12) * ratdum(irg3a)             
     4            - y(ic12) * y(ihe4) * ratdum(ircag)
     5            + y(io16) * ratdum(iroga)

c..o16 reactions 
      dydt(io16) = -y(ic12) * y(io16) * ratdum(ir1216)       
     1            - y(io16) * y(io16) * ratdum(ir1616) 
     2            + y(ic12) * y(ihe4) * ratdum(ircag)
     3            - y(io16) * y(ihe4) * ratdum(iroag)
     4            - y(io16) * ratdum(iroga)             
     5            + y(ine20) * ratdum(irnega)

c..ne20 reactions 
      dydt(ine20) =  0.5d0 * y(ic12) * y(ic12) * ratdum(ir1212)      
     1             + y(io16) * y(ihe4) * ratdum(iroag)       
     2             - y(ine20) * y(ihe4) * ratdum(irneag)
     3             - y(ine20) * ratdum(irnega)
     4             + y(img24) * ratdum(irmgga)            

c..mg24 reactions
      dydt(img24)  = 0.5d0 * y(ic12) * y(io16) * ratdum(ir1216)
     1              + y(ine20) * y(ihe4) * ratdum(irneag)
     2              - y(img24) * y(ihe4) * ratdum(irmgag) 
     4              - y(img24) * ratdum(irmgga) 
     5              + y(isi28) * ratdum(irsiga)
     6              - y(img24) *  y(ihe4) *  ratdum(irmgap) * (1.0d0-r1)
     7              + y(isi28) * r1 * ratdum(irsigp)

c..si28 reactions
      dydt(isi28)  =  0.5d0 * y(ic12) * y(io16) * ratdum(ir1216) 
     1              + 0.56d0 * 0.5d0 * y(io16)*y(io16) * ratdum(ir1616)
     2              + 0.34d0 * 0.5d0 * y(io16)*y(io16)*s1*ratdum(ir1616)
     3              + y(img24) * y(ihe4) * ratdum(irmgag) 
     4              - y(isi28) * y(ihe4) * ratdum(irsiag) 
     5              - y(isi28) * ratdum(irsiga)
     6              + y(is32)  * ratdum(irsga) 
     7              + y(img24) * y(ihe4) * ratdum(irmgap) * (1.0d0-r1)
     8              - y(isi28) * r1 * ratdum(irsigp) 
     9              - y(isi28) * y(ihe4) * ratdum(irsiap) * (1.0d0-s1)
     &              + y(is32)  * s1 * ratdum(irsgp)

c..s32 reactions
      dydt(is32) =   0.34d0*0.5d0*y(io16)**2 * ratdum(ir1616)*(1.0d0-s1)
     1             + 0.1d0 * 0.5d0 * y(io16) * y(io16) * ratdum(ir1616) 
     2             +  y(isi28) * y(ihe4) * ratdum(irsiag) 
     3             - y(is32) * y(ihe4) * ratdum(irsag) 
     4             - y(is32) * ratdum(irsga) 
     5             + y(iar36) * ratdum(irarga)
     6             + y(isi28) * y(ihe4) * ratdum(irsiap) * (1.0d0-s1)
     7             - y(is32) * s1 * ratdum(irsgp)
     8             - y(is32) * y(ihe4) * ratdum(irsap) * (1.0d0-t1)
     9             + y(iar36) * t1 * ratdum(irargp)


c..ar36 reactions
      dydt(iar36) =  y(is32)  * y(ihe4) * ratdum(irsag)
     1             - y(iar36) * y(ihe4) * ratdum(irarag)
     2             - y(iar36) * ratdum(irarga) 
     3             + y(ica40) * ratdum(ircaga)
     4             + y(is32)  * y(ihe4) * ratdum(irsap) * (1.0d0-t1) 
     5             - y(iar36) * t1 * ratdum(irargp)
     6             - y(iar36) * y(ihe4) * ratdum(irarap) * (1.0d0-u1)
     7             + y(ica40) * ratdum(ircagp) * u1


c..ca40 reactions
      dydt(ica40) =  y(iar36) * y(ihe4) * ratdum(irarag)
     1             - y(ica40) * y(ihe4) * ratdum(ircaag)
     2             - y(ica40) * ratdum(ircaga) 
     3             + y(iti44) * ratdum(irtiga) 
     4             + y(iar36) * y(ihe4) * ratdum(irarap) * (1.0d0-u1)
     5             - y(ica40) * ratdum(ircagp) * u1
     6             - y(ica40) * y(ihe4) * ratdum(ircaap) * (1.0d0-v1)
     7             + y(iti44) * ratdum(irtigp) * v1

c..ti44 reactions
      dydt(iti44) =  y(ica40) * y(ihe4) * ratdum(ircaag)
     1             - y(iti44) * y(ihe4) * ratdum(irtiag)
     2             - y(iti44) * ratdum(irtiga) 
     3             + y(icr48) * ratdum(ircrga)
     4             + y(ica40) * y(ihe4) * ratdum(ircaap)*(1.0d0-v1)
     5             - y(iti44) * v1 * ratdum(irtigp)
     6             - y(iti44) * y(ihe4) * ratdum(irtiap) * (1.0d0-w1)
     7             + y(icr48) * w1 * ratdum(ircrgp)

c..cr48 reactions
      dydt(icr48) =  y(iti44) * y(ihe4) * ratdum(irtiag)
     1             - y(icr48) * y(ihe4) * ratdum(ircrag)
     2             - y(icr48) * ratdum(ircrga) 
     3             + y(ife52) * ratdum(irfega) 
     4             + y(iti44) * y(ihe4) * ratdum(irtiap)*(1.0d0-w1)
     5             - y(icr48) * w1 * ratdum(ircrgp)
     6             - y(icr48) * y(ihe4) * ratdum(ircrap) * (1.0d0-x1)
     7             + y(ife52) * x1 * ratdum(irfegp)

c..fe52 reactions
      dydt(ife52) =  y(icr48) * y(ihe4) * ratdum(ircrag)
     1             - y(ife52) * y(ihe4) * ratdum(irfeag)
     2             - y(ife52) * ratdum(irfega)
     3             + y(ini56) * ratdum(irniga)
     4             + y(icr48) * y(ihe4) * ratdum(ircrap) * (1.0d0-x1) 
     5             - y(ife52) * x1 * ratdum(irfegp)
     6             - y(ife52) * y(ihe4) * ratdum(irfeap) * (1.0d0-y1)
     7             + y(ini56) * y1 * ratdum(irnigp)

c..ni56 reactions
      dydt(ini56) =  y(ife52) * y(ihe4) * ratdum(irfeag)
     1             - y(ini56) * ratdum(irniga)
     2             + y(ife52) * y(ihe4) * ratdum(irfeap) * (1.0d0-y1)  
     3             - y(ini56) * y1 * ratdum(irnigp)




c..instantaneous energy generation rate
      enuc = 0.0d0
      do i=1,ionmax
       enuc = enuc + dydt(i) * bion(i)
      enddo
      enuc = enuc * conv


c..get the neutrino losses
      call sneut5(btemp,bden,abar,zbar,
     1            sneut,snudt,snudd,snuda,snudz)


c..append an energy equation
      sdot        = enuc
      dydt(iener) = enuc - sneut




c..the type of temperature and density ode's depend
c..on the burn mode:


c..hydrostatic or single step cases
      if (hydrostatic  .or.  one_step) then
       dydt(itemp) = 0.0d0
       dydt(iden)  = 0.0d0



c..adiabatic expansion or contraction
      else if (expansion) then

       taud = 446.0d0/sqrt(den0) 
       taut = 3.0d0 * taud  

       dydt(itemp) = -psi * y(itemp)/taut
       dydt(iden)  = -psi * y(iden)/taud




c..self heating
      else if (self_heat) then


c..call an eos 
       temp_row(1) = btemp
       den_row(1)  = bden
       abar_row(1) = abar
       zbar_row(1) = zbar
       jlo_eos = 1
       jhi_eos = 1

       call helmeos


c..density equation
       dydt(iden) = 0.0d0


c..sum1 is d(ener)/d(yi)*d(yi)/dt = d(ener)/d(abar)*d(abar)/d(yi)*d(yi)/dt
       sum1 = 0.0d0
       do i=1,ionmax
        sum1 = sum1 - dydt(i)
       enddo
       sum1 = sum1 * dea_row(1)*abar*abar


c..sum2 is d(ener)/d(yi)*d(yi)/dt = d(ener)/d(zbar)*d(zbar)/d(yi)*d(yi)/dt
       sum2 = 0.0d0
       do i=1,ionmax
        sum2 = sum2 + (zion(i) - zbar)*dydt(i)
       enddo
       sum2 = sum2 * dez_row(1)*abar


c..temperature equation that is self-consistent with an eos 
       z           = 1.0d0/cv_row(1)
       dydt(itemp) = z*(dydt(iener) - ded_row(1)*dydt(iden) 
     1                   - sum1 - sum2)

      end if 


      return
      end   






      subroutine daprox13(tt,y,dfdy,nlog,nphys)   
      include 'implno.dek'
      include 'const.dek'
      include 'burn_common.dek'
      include 'network.dek'
      include 'vector_eos.dek'

c..this routine sets up the dense aprox13 jacobian

c..declare the pass
      integer          nlog,nphys
      double precision tt,y(1),dfdy(nphys,nphys)

c..local variables
      integer          i,j
      double precision denom,denomdt,denomdd,sum1,sum2,
     1                 r1,r1dt,r1dd,s1,s1dt,s1dd,t1,t1dt,t1dd,
     2                 u1,u1dt,u1dd,v1,v1dt,v1dd,w1,w1dt,w1dd,
     3                 x1,x1dt,x1dd,y1,y1dt,y1dd,zz,xx

      double precision zbarxx,ytot1,abar,zbar,taud,taut,
     1                 snudt,snudd,snuda,snudz

      double precision conv,sixth
      parameter        (conv = ev2erg*1.0d6*avo,
     1                  sixth = 1.0d0/6.0d0)



c..zero the jacobian
      do j=1,nlog
       do i=1,nlog
        dfdy(i,j) = 0.0d0
       enddo
      enddo


c..positive definite mass fractions
      do i=1,ionmax
       y(i) = min(1.0d0,max(y(i),1.0d-30))
      enddo


c..positive definite temperatures and densities
      y(itemp) = min(1.0d11,max(y(itemp),1.0d4))
      y(iden)  = min(1.0d11,max(y(iden),1.0d-10))



c..set the common block temperature and density
      btemp = y(itemp)
      bden  = y(iden)


c..generate abar and zbar for this composition
      zbarxx  = 0.0d0
      ytot1   = 0.0d0
      do i=1,ionmax
       ytot1    = ytot1 + y(i)
       zbarxx   = zbarxx + zion(i) * y(i)
      enddo
      abar  = 1.0d0/ytot1
      zbar  = zbarxx * abar



c..get the reaction rates
      if (use_tables .eq. 1) then
       call aprox13tab
      else
       call aprox13rat
      end if



c..do the screening here because the corrections depend on the composition

      call screen_aprox13(y)



c..branching ratios for various dummy proton links
      r1       = 0.0d0
      r1dt     = 0.0d0
      r1dd     = 0.0d0
      denom    = ratdum(iralpa) + ratdum(iralpg)
      denomdt  = dratdumdt(iralpa) + dratdumdt(iralpg)
      denomdd  = dratdumdd(iralpa) + dratdumdd(iralpg)
      if (denom .ne. 0.0) then
       zz   = 1.0d0/denom
       r1   = ratdum(iralpa)*zz
       r1dt = (dratdumdt(iralpa) - r1*denomdt)*zz
       r1dd = (dratdumdd(iralpa) - r1*denomdd)*zz
      end if

      s1       = 0.0d0
      s1dt     = 0.0d0
      s1dd     = 0.0d0
      denom    = ratdum(irppa) + ratdum(irppg)
      denomdt  = dratdumdt(irppa) + dratdumdt(irppg)
      denomdd  = dratdumdd(irppa) + dratdumdd(irppg)
      if (denom .ne. 0.0) then
       zz   = 1.0d0/denom
       s1   = ratdum(irppa)*zz
       s1dt = (dratdumdt(irppa) - s1*denomdt)*zz
       s1dd = (dratdumdd(irppa) - s1*denomdd)*zz
      end if

      t1       = 0.0d0
      t1dt     = 0.0d0
      t1dd     = 0.0d0
      denom    = ratdum(irclpa) + ratdum(irclpg)
      denomdt  = dratdumdt(irclpa) + dratdumdt(irclpg)
      denomdd  = dratdumdd(irclpa) + dratdumdd(irclpg)
      if (denom .ne. 0.0) then
       zz   = 1.0d0/denom
       t1   = ratdum(irclpa)*zz
       t1dt = (dratdumdt(irclpa) - t1*denomdt)*zz
       t1dd = (dratdumdd(irclpa) - t1*denomdd)*zz
      end if

      u1       = 0.0d0
      u1dt     = 0.0d0
      u1dd     = 0.0d0
      denom    = ratdum(irkpa) + ratdum(irkpg)
      denomdt  = dratdumdt(irkpa) + dratdumdt(irkpg)
      denomdd  = dratdumdd(irkpa) + dratdumdd(irkpg)
      if (denom .ne. 0.0) then
       zz   = 1.0d0/denom
       u1   = ratdum(irkpa)*zz
       u1dt = (dratdumdt(irkpa) - u1*denomdt)*zz
       u1dd = (dratdumdd(irkpa) - u1*denomdd)*zz
      end if 

      v1       = 0.0d0
      v1dt     = 0.0d0
      v1dd     = 0.0d0
      denom    = ratdum(irscpa) + ratdum(irscpg)
      denomdt  = dratdumdt(irscpa) + dratdumdt(irscpg)
      denomdd  = dratdumdd(irscpa) + dratdumdd(irscpg)
      if (denom .ne. 0.0) then
       zz   = 1.0d0/denom
       v1   = ratdum(irscpa)*zz
       v1dt = (dratdumdt(irscpa) - v1*denomdt)*zz
       v1dd = (dratdumdd(irscpa) - v1*denomdd)*zz
      end if 

      w1       = 0.0d0
      w1dt     = 0.0d0
      w1dd     = 0.0d0
      denom    = ratdum(irvpa) + ratdum(irvpg)
      denomdt  = dratdumdt(irvpa) + dratdumdt(irvpg)
      denomdd  = dratdumdd(irvpa) + dratdumdd(irvpg)
      if (denom .ne. 0.0) then
       zz   = 1.0d0/denom
       w1   = ratdum(irvpa)*zz
       w1dt = (dratdumdt(irvpa) - w1*denomdt)*zz
       w1dd = (dratdumdd(irvpa) - w1*denomdd)*zz
      end if

      x1       = 0.0d0
      x1dt     = 0.0d0
      x1dd     = 0.0d0
      denom    = ratdum(irmnpa) + ratdum(irmnpg)
      denomdt  = dratdumdt(irmnpa) + dratdumdt(irmnpg)
      denomdd  = dratdumdd(irmnpa) + dratdumdd(irmnpg)
      if (denom .ne. 0.0) then
       zz   = 1.0d0/denom
       x1   = ratdum(irmnpa)*zz
       x1dt = (dratdumdt(irmnpa) - x1*denomdt)*zz
       x1dd = (dratdumdd(irmnpa) - x1*denomdd)*zz
      endif

      y1       = 0.0d0
      y1dt     = 0.0d0
      y1dd     = 0.0d0
      denom    = ratdum(ircopa) + ratdum(ircopg)
      denomdt  = dratdumdt(ircopa) + dratdumdt(ircopg)
      denomdd  = dratdumdd(ircopa) + dratdumdd(ircopg)
      if (denom .ne. 0.0) then
       zz   = 1.0d0/denom
       y1   = ratdum(ircopa)*zz
       y1dt = (dratdumdt(ircopa) - y1*denomdt)*zz
       y1dd = (dratdumdd(ircopa) - y1*denomdd)*zz
      end if



c..he4 jacobian elements
      dfdy(ihe4,ihe4)  = -1.5d0 * y(ihe4) * y(ihe4) * ratdum(ir3a)   
     1                   - y(ic12)  * ratdum(ircag)      
     2                   - y(io16)  * ratdum(iroag)      
     3                   - y(ine20) * ratdum(irneag) 
     4                   - y(img24) * ratdum(irmgag)
     5                   - y(isi28) * ratdum(irsiag) 
     6                   - y(is32)  * ratdum(irsag) 
     7                   - y(iar36) * ratdum(irarag)
     8                   - y(ica40) * ratdum(ircaag)
     9                   - y(iti44) * ratdum(irtiag)
     &                   - y(icr48) * ratdum(ircrag)
     1                   - y(ife52) * ratdum(irfeag) 

      dfdy(ihe4,ihe4)  = dfdy(ihe4,ihe4)
     1                   - y(img24) * ratdum(irmgap) * (1.0d0-r1)
     2                   - y(isi28) * ratdum(irsiap) * (1.0d0-s1)
     3                   - y(is32) * ratdum(irsap)   * (1.0d0-t1)
     4                   - y(iar36) * ratdum(irarap) * (1.0d0-u1)
     5                   - y(ica40) * ratdum(ircaap) * (1.0d0-v1) 
     6                   - y(iti44) * ratdum(irtiap) * (1.0d0-w1)
     7                   - y(icr48) * ratdum(ircrap) * (1.0d0-x1) 
     8                   - y(ife52) * ratdum(irfeap) * (1.0d0-y1)


      dfdy(ihe4,ic12)  = y(ic12) * ratdum(ir1212)   
     1                   + 0.5d0 * y(io16) * ratdum(ir1216)   
     2                   + 3.0d0 * ratdum(irg3a)  
     3                   - y(ihe4) * ratdum(ircag) 

      dfdy(ihe4,io16)  = 0.5d0 * y(ic12) * ratdum(ir1216)
     1                   + 1.12d0 * 0.5d0 * y(io16) * ratdum(ir1616) 
     2                   + 0.68d0 * s1 * 0.5d0*y(io16) * ratdum(ir1616) 
     3                   + ratdum(iroga)
     4                   - y(ihe4) * ratdum(iroag)          

      dfdy(ihe4,ine20) =  ratdum(irnega)
     1                  - y(ihe4) * ratdum(irneag)

      dfdy(ihe4,img24) =   ratdum(irmgga)
     1                   - y(ihe4) * ratdum(irmgag) 
     2                   - y(ihe4) * ratdum(irmgap) * (1.0d0-r1)

      dfdy(ihe4,isi28) =   ratdum(irsiga)
     1                   - y(ihe4) * ratdum(irsiag)
     2                   - y(ihe4) * ratdum(irsiap) * (1.0d0-s1)
     3                   + r1 * ratdum(irsigp) 


      dfdy(ihe4,is32)  =   ratdum(irsga)
     1                   - y(ihe4) * ratdum(irsag)
     2                   - y(ihe4) * ratdum(irsap) * (1.0d0-t1)
     3                   + s1 * ratdum(irsgp)

      dfdy(ihe4,iar36) =   ratdum(irarga)
     1                   - y(ihe4) * ratdum(irarag)
     2                   - y(ihe4) * ratdum(irarap) * (1.0d0-u1)
     3                   + t1 * ratdum(irargp)

      dfdy(ihe4,ica40) =   ratdum(ircaga)
     1                   - y(ihe4) * ratdum(ircaag)
     2                   - y(ihe4) * ratdum(ircaap) * (1.0d0-v1)
     3                   + u1 * ratdum(ircagp)

      dfdy(ihe4,iti44) =   ratdum(irtiga)
     1                   - y(ihe4) * ratdum(irtiag)
     2                   - y(ihe4) * ratdum(irtiap) * (1.0d0-w1)
     3                   + v1 * ratdum(irtigp)

      dfdy(ihe4,icr48) =   ratdum(ircrga)
     1                   - y(ihe4) * ratdum(ircrag)
     2                   - y(ihe4) * ratdum(ircrap) * (1.0d0-x1)
     3                   + w1 * ratdum(ircrgp)

      dfdy(ihe4,ife52) =   ratdum(irfega)
     1                   - y(ihe4) * ratdum(irfeag) 
     2                   - y(ihe4) * ratdum(irfeap) * (1.0d0-y1)
     3                   + x1 * ratdum(irfegp) 

      dfdy(ihe4,ini56) =   ratdum(irniga)
     1                   + y1 * ratdum(irnigp) 


c..c12 jacobian elements
      dfdy(ic12,ihe4) = 0.5d0 * y(ihe4) * y(ihe4) * ratdum(ir3a)
     1                 - y(ic12) * ratdum(ircag)   

      dfdy(ic12,ic12) = -2.0d0 * y(ic12) * ratdum(ir1212) 
     1                  - y(io16) * ratdum(ir1216)  
     2                  - ratdum(irg3a)
     3                  - y(ihe4) * ratdum(ircag)

      dfdy(ic12,io16) = -y(ic12) * ratdum(ir1216)   
     1                 + ratdum(iroga)



c..o16 jacobian elements
      dfdy(io16,ihe4) = y(ic12)*ratdum(ircag)
     1                - y(io16)*ratdum(iroag) 

      dfdy(io16,ic12) = -y(io16)*ratdum(ir1216) 
     1                 + y(ihe4)*ratdum(ircag)

      dfdy(io16,io16) = - y(ic12) * ratdum(ir1216)  
     1                  - 2.0d0 * y(io16) * ratdum(ir1616) 
     2                  - y(ihe4) * ratdum(iroag)
     3                  - ratdum(iroga) 

      dfdy(io16,ine20) = ratdum(irnega)



c..ne20 jacobian elements
      dfdy(ine20,ihe4)  = y(io16) * ratdum(iroag)  
     1                  - y(ine20) * ratdum(irneag)       

      dfdy(ine20,ic12)  = y(ic12) * ratdum(ir1212) 

      dfdy(ine20,io16)  = y(ihe4) * ratdum(iroag)

      dfdy(ine20,ine20) = -y(ihe4) * ratdum(irneag) 
     1                   - ratdum(irnega)

      dfdy(ine20,img24) = ratdum(irmgga)



c..mg24 jacobian elements
      dfdy(img24,ihe4)  = y(ine20) * ratdum(irneag)
     1                   -y(img24) * ratdum(irmgag) 
     2                   -y(img24) * ratdum(irmgap) * (1.0d0-r1) 

      dfdy(img24,ic12)  = 0.5d0 * y(io16) * ratdum(ir1216)

      dfdy(img24,io16)  = 0.5d0 * y(ic12) * ratdum(ir1216)

      dfdy(img24,ine20) = y(ihe4) * ratdum(irneag)

      dfdy(img24,img24) = -y(ihe4) * ratdum(irmgag) 
     1                   - ratdum(irmgga) 
     2                   - y(ihe4) * ratdum(irmgap) * (1.0d0-r1)  

      dfdy(img24,isi28) = ratdum(irsiga) 
     1                  + r1 * ratdum(irsigp)



c..si28 jacobian elements
      dfdy(isi28,ihe4)  = y(img24) * ratdum(irmgag)
     1                  - y(isi28) * ratdum(irsiag) 
     2                  + y(img24) * ratdum(irmgap) * (1.0d0-r1)
     3                  - y(isi28) * ratdum(irsiap) * (1.0d0-s1)

      dfdy(isi28,ic12)  = 0.5d0 * y(io16) * ratdum(ir1216) 

      dfdy(isi28,io16)  =   0.5d0 * y(ic12) * ratdum(ir1216) 
     1                    + 1.12d0 * 0.5d0*y(io16) * ratdum(ir1616)  
     2                    + 0.68d0 * 0.5d0*y(io16) * s1 * ratdum(ir1616)     

      dfdy(isi28,img24) = y(ihe4) * ratdum(irmgag) 
     1                  + y(ihe4) * ratdum(irmgap) * (1.0d0-r1)

      dfdy(isi28,isi28) = -y(ihe4) * ratdum(irsiag) 
     1                   - ratdum(irsiga)
     2                   - r1 * ratdum(irsigp) 
     3                   - y(ihe4) * ratdum(irsiap) * (1.0d0-s1)

      dfdy(isi28,is32)  = ratdum(irsga) 
     1                  + s1 * ratdum(irsgp)



c..s32 jacobian elements
      dfdy(is32,ihe4)  = y(isi28) * ratdum(irsiag)
     1                 - y(is32) * ratdum(irsag)
     2                 + y(isi28) * ratdum(irsiap) * (1.0d0-s1)
     3                 - y(is32) * ratdum(irsap) * (1.0d0-t1)

      dfdy(is32,io16)  = 0.68d0*0.5d0*y(io16)*ratdum(ir1616)*(1.0d0-s1)
     1                   + 0.2d0 * 0.5d0*y(io16) * ratdum(ir1616) 

      dfdy(is32,isi28) = y(ihe4) * ratdum(irsiag) 
     1                  + y(ihe4) * ratdum(irsiap) * (1.0d0-s1)

      dfdy(is32,is32)  = -y(ihe4) * ratdum(irsag)
     1                  - ratdum(irsga) 
     2                  - s1 * ratdum(irsgp)
     3                  - y(ihe4) * ratdum(irsap) * (1.0d0-t1)

      dfdy(is32,iar36) = ratdum(irarga) 
     1                 + t1 * ratdum(irargp)



c..ar36 jacobian elements
      dfdy(iar36,ihe4)  = y(is32)  * ratdum(irsag)
     1                  - y(iar36) * ratdum(irarag)
     2                  + y(is32)  * ratdum(irsap) * (1.0d0-t1) 
     3                  - y(iar36) * ratdum(irarap) * (1.0d0-u1)

      dfdy(iar36,is32)  = y(ihe4) * ratdum(irsag)
     1                   + y(ihe4) * ratdum(irsap) * (1.0d0-t1)

      dfdy(iar36,iar36) = -y(ihe4) * ratdum(irarag) 
     1                   - ratdum(irarga) 
     2                   - t1 * ratdum(irargp)
     3                   - y(ihe4) * ratdum(irarap) * (1.0d0-u1)

      dfdy(iar36,ica40) = ratdum(ircaga) 
     1                  + ratdum(ircagp) * u1



c..ca40 jacobian elements
      dfdy(ica40,ihe4)   = y(iar36) * ratdum(irarag)
     1                   - y(ica40) * ratdum(ircaag) 
     2                   + y(iar36) * ratdum(irarap)*(1.0d0-u1)
     3                   - y(ica40) * ratdum(ircaap)*(1.0d0-v1)

      dfdy(ica40,iar36)  = y(ihe4) * ratdum(irarag)
     1                   + y(ihe4) * ratdum(irarap)*(1.0d0-u1)

      dfdy(ica40,ica40)  = -y(ihe4) * ratdum(ircaag)
     1                    - ratdum(ircaga) 
     2                    - ratdum(ircagp) * u1
     3                    - y(ihe4) * ratdum(ircaap)*(1.0d0-v1)

      dfdy(ica40,iti44)  = ratdum(irtiga) 
     1                    + ratdum(irtigp) * v1



c..ti44 jacobian elements
      dfdy(iti44,ihe4)   = y(ica40) * ratdum(ircaag)
     1                   - y(iti44) * ratdum(irtiag)
     2                   + y(ica40) * ratdum(ircaap)*(1.0d0-v1)
     3                   - y(iti44) * ratdum(irtiap)*(1.0d0-w1)

      dfdy(iti44,ica40)  = y(ihe4) * ratdum(ircaag)
     1                   + y(ihe4) * ratdum(ircaap)*(1.0d0-v1)

      dfdy(iti44,iti44)  = -y(ihe4) * ratdum(irtiag)
     1                    - ratdum(irtiga) 
     2                    - v1 * ratdum(irtigp)
     3                    - y(ihe4) * ratdum(irtiap)*(1.0d0-w1)

      dfdy(iti44,icr48)  = ratdum(ircrga) 
     1                   + w1 * ratdum(ircrgp)



c..cr48 jacobian elements
      dfdy(icr48,ihe4)  = y(iti44) * ratdum(irtiag)
     1                  - y(icr48) * ratdum(ircrag)
     2                  + y(iti44) * ratdum(irtiap)*(1.0d0-w1)
     3                  - y(icr48) * ratdum(ircrap)*(1.0d0-x1)

      dfdy(icr48,iti44) = y(ihe4) * ratdum(irtiag)
     1                  + y(ihe4) * ratdum(irtiap)*(1.0d0-w1)

      dfdy(icr48,icr48) = -y(ihe4) * ratdum(ircrag) 
     1                   - ratdum(ircrga) 
     2                   - w1 * ratdum(ircrgp)
     3                   - y(ihe4) * ratdum(ircrap)*(1.0d0-x1)

      dfdy(icr48,ife52) = ratdum(irfega) 
     1                  + x1 * ratdum(irfegp)



c..fe52 jacobian elements
      dfdy(ife52,ihe4)  = y(icr48) * ratdum(ircrag)
     1                  - y(ife52) * ratdum(irfeag)
     2                  + y(icr48) * ratdum(ircrap) * (1.0d0-x1) 
     3                  - y(ife52) * ratdum(irfeap) * (1.0d0-y1)

      dfdy(ife52,icr48) = y(ihe4) * ratdum(ircrag)
     1                  + y(ihe4) * ratdum(ircrap) * (1.0d0-x1) 

      dfdy(ife52,ife52) = - y(ihe4) * ratdum(irfeag)
     1                    - ratdum(irfega)
     2                    - x1 * ratdum(irfegp)
     3                    - y(ihe4) * ratdum(irfeap) * (1.0d0-y1)

      dfdy(ife52,ini56) = ratdum(irniga)
     1                  + y1 * ratdum(irnigp)



c..ni56 jacobian elements
      dfdy(ini56,ihe4)  = y(ife52) * ratdum(irfeag)
     1                  + y(ife52) * ratdum(irfeap) * (1.0d0-y1)

      dfdy(ini56,ife52) = y(ihe4) * ratdum(irfeag)
     1                  + y(ihe4) * ratdum(irfeap) * (1.0d0-y1)

      dfdy(ini56,ini56) = -ratdum(irniga)
     1                   - y1 * ratdum(irnigp)



c..append the temperature derivatives of the rate equations

c..he4 reactions and heavy ion reactions
      dfdy(ihe4,itemp) =  0.5d0 * y(ic12) * y(ic12) * dratdumdt(ir1212)      
     1            + 0.5d0 * y(ic12) * y(io16) * dratdumdt(ir1216)
     2            + 0.56d0*0.5d0*y(io16) * y(io16) * dratdumdt(ir1616) 
     3            + 0.34d0 *s1*0.5d0*y(io16)*y(io16)*dratdumdt(ir1616)
     3            + 0.34d0 *s1dt*0.5d0*y(io16)*y(io16)*ratdum(ir1616)

c..(a,g) and (g,a) reactions
      dfdy(ihe4,itemp) =  dfdy(ihe4,itemp)
     1            - 0.5d0 * y(ihe4)*y(ihe4)*y(ihe4) * dratdumdt(ir3a) 
     2            + 3.0d0 * y(ic12) * dratdumdt(irg3a)       
     3            - y(ihe4)  * y(ic12) * dratdumdt(ircag)       
     4            + y(io16)  * dratdumdt(iroga)
     5            - y(ihe4)  * y(io16) * dratdumdt(iroag)       
     6            + y(ine20) * dratdumdt(irnega)            
     7            - y(ihe4)  * y(ine20) * dratdumdt(irneag) 
     8            + y(img24) * dratdumdt(irmgga)
     9            - y(ihe4)  * y(img24)* dratdumdt(irmgag)
     &            + y(isi28) * dratdumdt(irsiga) 
     1            - y(ihe4)  * y(isi28)*dratdumdt(irsiag)
     2            + y(is32)  * dratdumdt(irsga) 

      dfdy(ihe4,itemp) =  dfdy(ihe4,itemp)
     1            - y(ihe4)  * y(is32) * dratdumdt(irsag)
     2            + y(iar36) * dratdumdt(irarga) 
     3            - y(ihe4)  * y(iar36)*dratdumdt(irarag)
     4            + y(ica40) * dratdumdt(ircaga) 
     5            - y(ihe4)  * y(ica40)*dratdumdt(ircaag)
     6            + y(iti44) * dratdumdt(irtiga) 
     7            - y(ihe4)  * y(iti44)*dratdumdt(irtiag)
     8            + y(icr48) * dratdumdt(ircrga) 
     9            - y(ihe4)  * y(icr48)*dratdumdt(ircrag)
     &            + y(ife52) * dratdumdt(irfega) 
     1            - y(ihe4)  * y(ife52) * dratdumdt(irfeag)
     2            + y(ini56) * dratdumdt(irniga)


c..(a,p)(p,g) and (g,p)(p,a) reactions
      dfdy(ihe4,itemp) =  dfdy(ihe4,itemp)
     1            - y(ihe4)  * y(img24) * dratdumdt(irmgap) * (1.0d0-r1)
     2            + y(isi28) * dratdumdt(irsigp) * r1
     3            - y(ihe4)  * y(isi28) * dratdumdt(irsiap) * (1.0d0-s1)
     4            + y(is32)  * dratdumdt(irsgp) * s1
     5            - y(ihe4)  * y(is32) * dratdumdt(irsap) * (1.0d0-t1)   
     6            + y(iar36) * dratdumdt(irargp) * t1 
     7            - y(ihe4)  * y(iar36) * dratdumdt(irarap) * (1.0d0-u1)
     8            + y(ica40) * dratdumdt(ircagp) * u1
     9            - y(ihe4)  * y(ica40) * dratdumdt(ircaap) * (1.0d0-v1)   
     &            + y(iti44) * dratdumdt(irtigp) * v1
     1            - y(ihe4)  * y(iti44) * dratdumdt(irtiap) * (1.0d0-w1)
     2            + y(icr48) * dratdumdt(ircrgp) * w1
     3            - y(ihe4)  * y(icr48) * dratdumdt(ircrap) * (1.0d0-x1)   
     4            + y(ife52) * dratdumdt(irfegp) * x1 
     5            - y(ihe4)  * y(ife52) * dratdumdt(irfeap) * (1.0d0-y1)
     6            + y(ini56) * dratdumdt(irnigp) * y1


      dfdy(ihe4,itemp) =  dfdy(ihe4,itemp)
     1            + y(ihe4)  * y(img24) * ratdum(irmgap)* r1dt
     2            + y(isi28) * ratdum(irsigp) * r1dt
     3            + y(ihe4)  * y(isi28) * ratdum(irsiap) * s1dt
     4            + y(is32)  * ratdum(irsgp) * s1dt
     5            + y(ihe4)  * y(is32) * ratdum(irsap) * t1dt   
     6            + y(iar36) * ratdum(irargp) * t1dt 
     7            + y(ihe4)  * y(iar36) * ratdum(irarap) * u1dt
     8            + y(ica40) * ratdum(ircagp) * u1dt
     9            + y(ihe4)  * y(ica40) * ratdum(ircaap) * v1dt   
     &            + y(iti44) * ratdum(irtigp) * v1dt
     1            + y(ihe4)  * y(iti44) * ratdum(irtiap) * w1dt
     2            + y(icr48) * ratdum(ircrgp) * w1dt
     3            + y(ihe4)  * y(icr48) * ratdum(ircrap) * x1dt   
     4            + y(ife52) * ratdum(irfegp) * x1dt 
     5            + y(ihe4)  * y(ife52) * ratdum(irfeap) * y1dt
     6            + y(ini56) * ratdum(irnigp) * y1dt


c..c12 reactions   
      dfdy(ic12,itemp) = 
     1             -y(ic12) * y(ic12) * dratdumdt(ir1212) 
     1            - y(ic12) * y(io16) * dratdumdt(ir1216)       
     2            + sixth * y(ihe4)*y(ihe4)*y(ihe4) * dratdumdt(ir3a)
     3            - y(ic12) * dratdumdt(irg3a)             
     4            - y(ic12) * y(ihe4) * dratdumdt(ircag)
     5            + y(io16) * dratdumdt(iroga)

c..o16 reactions 
      dfdy(io16,itemp) = 
     1             -y(ic12) * y(io16) * dratdumdt(ir1216)       
     1            - y(io16) * y(io16) * dratdumdt(ir1616) 
     2            + y(ic12) * y(ihe4) * dratdumdt(ircag)
     3            - y(io16) * y(ihe4) * dratdumdt(iroag)
     4            - y(io16) * dratdumdt(iroga)             
     5            + y(ine20) * dratdumdt(irnega)

c..ne20 reactions 
      dfdy(ine20,itemp) =  
     1               0.5d0 * y(ic12) * y(ic12) * dratdumdt(ir1212)      
     1             + y(io16) * y(ihe4) * dratdumdt(iroag)       
     2             - y(ine20) * y(ihe4) * dratdumdt(irneag)
     3             - y(ine20) * dratdumdt(irnega)
     4             + y(img24) * dratdumdt(irmgga)            

c..mg24 reactions
      dfdy(img24,itemp)  = 
     1                0.5d0 * y(ic12) * y(io16) * dratdumdt(ir1216)
     1              + y(ine20) * y(ihe4) * dratdumdt(irneag)
     2              - y(img24) * y(ihe4) * dratdumdt(irmgag) 
     4              - y(img24) * dratdumdt(irmgga) 
     5              + y(isi28) * dratdumdt(irsiga)
     6              - y(img24) *  y(ihe4) *dratdumdt(irmgap)*(1.0d0-r1)
     7              + y(isi28) * r1 * dratdumdt(irsigp)
     6              + y(img24) *  y(ihe4) *ratdum(irmgap)*r1dt
     7              + y(isi28) * r1dt * ratdum(irsigp)

c..si28 reactions
      dfdy(isi28,itemp)  =  
     1                0.5d0 * y(ic12) * y(io16) * dratdumdt(ir1216) 
     1              + 0.56d0*0.5d0*y(io16)*y(io16) * dratdumdt(ir1616)
     2              + 0.34d0*0.5d0*y(io16)*y(io16)*s1*dratdumdt(ir1616)
     2              + 0.34d0*0.5d0*y(io16)*y(io16)*s1dt*ratdum(ir1616)
     3              + y(img24) * y(ihe4) * dratdumdt(irmgag) 
     4              - y(isi28) * y(ihe4) * dratdumdt(irsiag) 
     5              - y(isi28) * dratdumdt(irsiga)
     6              + y(is32)  * dratdumdt(irsga) 
     7              + y(img24) * y(ihe4) * dratdumdt(irmgap)*(1.0d0-r1)
     8              - y(isi28) * r1 * dratdumdt(irsigp) 
     7              - y(img24) * y(ihe4) * ratdum(irmgap)*r1dt
     8              - y(isi28) * r1dt * ratdum(irsigp) 
     9              - y(isi28) * y(ihe4) * dratdumdt(irsiap)*(1.0d0-s1)
     &              + y(is32)  * s1 * dratdumdt(irsgp)
     9              + y(isi28) * y(ihe4) * ratdum(irsiap)*s1dt
     &              + y(is32)  * s1dt * ratdum(irsgp)

c..s32 reactions
      dfdy(is32,itemp) =  
     1              0.34d0*0.5d0*y(io16)**2*dratdumdt(ir1616)*(1.0d0-s1)
     1             - 0.34d0*0.5d0*y(io16)**2 * ratdum(ir1616)*s1dt
     1             + 0.1d0*0.5d0 * y(io16) * y(io16) * dratdumdt(ir1616) 
     2             +  y(isi28) * y(ihe4) * dratdumdt(irsiag) 
     3             - y(is32) * y(ihe4) * dratdumdt(irsag) 
     4             - y(is32) * dratdumdt(irsga) 
     5             + y(iar36) * dratdumdt(irarga)
     6             + y(isi28) * y(ihe4) * dratdumdt(irsiap) * (1.0d0-s1)
     7             - y(is32) * s1 * dratdumdt(irsgp)
     6             - y(isi28) * y(ihe4) * ratdum(irsiap) * s1dt
     7             - y(is32) * s1dt * ratdum(irsgp)
     8             - y(is32) * y(ihe4) * dratdumdt(irsap) * (1.0d0-t1)
     9             + y(iar36) * t1 * dratdumdt(irargp)
     8             - y(is32) * y(ihe4) * ratdum(irsap) * t1dt
     9             + y(iar36) * t1dt * ratdum(irargp)


c..ar36 reactions
      dfdy(iar36,itemp) =  
     1               y(is32)  * y(ihe4) * dratdumdt(irsag)
     1             - y(iar36) * y(ihe4) * dratdumdt(irarag)
     2             - y(iar36) * dratdumdt(irarga) 
     3             + y(ica40) * dratdumdt(ircaga)
     4             + y(is32)  * y(ihe4) * dratdumdt(irsap) * (1.0d0-t1) 
     5             - y(iar36) * t1 * dratdumdt(irargp)
     4             - y(is32)  * y(ihe4) * ratdum(irsap) * t1dt 
     5             - y(iar36) * t1dt * ratdum(irargp)
     6             - y(iar36) * y(ihe4) * dratdumdt(irarap) * (1.0d0-u1)
     7             + y(ica40) * dratdumdt(ircagp) * u1
     6             + y(iar36) * y(ihe4) * ratdum(irarap) * u1dt
     7             + y(ica40) * ratdum(ircagp) * u1dt


c..ca40 reactions
      dfdy(ica40,itemp) =  
     1               y(iar36) * y(ihe4) * dratdumdt(irarag)
     1             - y(ica40) * y(ihe4) * dratdumdt(ircaag)
     2             - y(ica40) * dratdumdt(ircaga) 
     3             + y(iti44) * dratdumdt(irtiga) 
     4             + y(iar36) * y(ihe4) * dratdumdt(irarap) * (1.0d0-u1)
     5             - y(ica40) * dratdumdt(ircagp) * u1
     4             - y(iar36) * y(ihe4) * ratdum(irarap) * u1dt
     5             - y(ica40) * ratdum(ircagp) * u1dt
     6             - y(ica40) * y(ihe4) * dratdumdt(ircaap) * (1.0d0-v1)
     7             + y(iti44) * dratdumdt(irtigp) * v1
     6             + y(ica40) * y(ihe4) * ratdum(ircaap) * v1dt
     7             + y(iti44) * ratdum(irtigp) * v1dt


c..ti44 reactions
      dfdy(iti44,itemp) =  
     1               y(ica40) * y(ihe4) * dratdumdt(ircaag)
     1             - y(iti44) * y(ihe4) * dratdumdt(irtiag)
     2             - y(iti44) * dratdumdt(irtiga) 
     3             + y(icr48) * dratdumdt(ircrga)
     4             + y(ica40) * y(ihe4) * dratdumdt(ircaap)*(1.0d0-v1)
     5             - y(iti44) * v1 * dratdumdt(irtigp)
     4             - y(ica40) * y(ihe4) * ratdum(ircaap)*v1dt
     5             - y(iti44) * v1dt * ratdum(irtigp)
     6             - y(iti44) * y(ihe4) * dratdumdt(irtiap) * (1.0d0-w1)
     7             + y(icr48) * w1 * dratdumdt(ircrgp)
     6             + y(iti44) * y(ihe4) * ratdum(irtiap) * w1dt
     7             + y(icr48) * w1dt * ratdum(ircrgp)


c..cr48 reactions
      dfdy(icr48,itemp) =  
     1               y(iti44) * y(ihe4) * dratdumdt(irtiag)
     1             - y(icr48) * y(ihe4) * dratdumdt(ircrag)
     2             - y(icr48) * dratdumdt(ircrga) 
     3             + y(ife52) * dratdumdt(irfega) 
     4             + y(iti44) * y(ihe4) * dratdumdt(irtiap)*(1.0d0-w1)
     5             - y(icr48) * w1 * dratdumdt(ircrgp)
     4             - y(iti44) * y(ihe4) * ratdum(irtiap)*w1dt
     5             - y(icr48) * w1dt * ratdum(ircrgp)
     6             - y(icr48) * y(ihe4) * dratdumdt(ircrap) * (1.0d0-x1)
     7             + y(ife52) * x1 * dratdumdt(irfegp)
     6             + y(icr48) * y(ihe4) * ratdum(ircrap) * x1dt
     7             + y(ife52) * x1dt * ratdum(irfegp)

c..fe52 reactions
      dfdy(ife52,itemp) =  
     1               y(icr48) * y(ihe4) * dratdumdt(ircrag)
     1             - y(ife52) * y(ihe4) * dratdumdt(irfeag)
     2             - y(ife52) * dratdumdt(irfega)
     3             + y(ini56) * dratdumdt(irniga)
     4             + y(icr48) * y(ihe4) * dratdumdt(ircrap) * (1.0d0-x1) 
     5             - y(ife52) * x1 * dratdumdt(irfegp)
     4             - y(icr48) * y(ihe4) * ratdum(ircrap) * x1dt 
     5             - y(ife52) * x1dt * ratdum(irfegp)
     6             - y(ife52) * y(ihe4) * dratdumdt(irfeap) * (1.0d0-y1)
     7             + y(ini56) * y1 * dratdumdt(irnigp)
     6             + y(ife52) * y(ihe4) * ratdum(irfeap) * y1dt
     7             + y(ini56) * y1dt * ratdum(irnigp)

c..ni56 reactions
      dfdy(ini56,itemp) =  
     1               y(ife52) * y(ihe4) * dratdumdt(irfeag)
     1             - y(ini56) * dratdumdt(irniga)
     2             + y(ife52) * y(ihe4) * dratdumdt(irfeap) * (1.0d0-y1)  
     3             - y(ini56) * y1 * dratdumdt(irnigp)
     2             - y(ife52) * y(ihe4) * ratdum(irfeap) * y1dt  
     3             - y(ini56) * y1dt * ratdum(irnigp)




c..append the density derivatives of the rate equations

c..he4 reactions and heavy ion reactions
      dfdy(ihe4,iden) =  0.5d0 * y(ic12) * y(ic12) * dratdumdd(ir1212)      
     1            + 0.5d0 * y(ic12) * y(io16) * dratdumdd(ir1216)
     2            + 0.56d0*0.5d0*y(io16) * y(io16) * dratdumdd(ir1616) 
     3            + 0.34d0 *s1*0.5d0*y(io16)*y(io16)*dratdumdd(ir1616)
     3            + 0.34d0 *s1dd*0.5d0*y(io16)*y(io16)*ratdum(ir1616)

c..(a,g) and (g,a) reactions
      dfdy(ihe4,iden) =  dfdy(ihe4,iden)
     1            - 0.5d0 * y(ihe4)*y(ihe4)*y(ihe4) * dratdumdd(ir3a) 
     2            + 3.0d0 * y(ic12) * dratdumdd(irg3a)       
     3            - y(ihe4)  * y(ic12) * dratdumdd(ircag)       
     4            + y(io16)  * dratdumdd(iroga)
     5            - y(ihe4)  * y(io16) * dratdumdd(iroag)       
     6            + y(ine20) * dratdumdd(irnega)            
     7            - y(ihe4)  * y(ine20) * dratdumdd(irneag) 
     8            + y(img24) * dratdumdd(irmgga)
     9            - y(ihe4)  * y(img24)* dratdumdd(irmgag)
     &            + y(isi28) * dratdumdd(irsiga) 
     1            - y(ihe4)  * y(isi28)*dratdumdd(irsiag)
     2            + y(is32)  * dratdumdd(irsga) 

      dfdy(ihe4,iden) =  dfdy(ihe4,iden)
     1            - y(ihe4)  * y(is32) * dratdumdd(irsag)
     2            + y(iar36) * dratdumdd(irarga) 
     3            - y(ihe4)  * y(iar36)*dratdumdd(irarag)
     4            + y(ica40) * dratdumdd(ircaga) 
     5            - y(ihe4)  * y(ica40)*dratdumdd(ircaag)
     6            + y(iti44) * dratdumdd(irtiga) 
     7            - y(ihe4)  * y(iti44)*dratdumdd(irtiag)
     8            + y(icr48) * dratdumdd(ircrga) 
     9            - y(ihe4)  * y(icr48)*dratdumdd(ircrag)
     &            + y(ife52) * dratdumdd(irfega) 
     1            - y(ihe4)  * y(ife52) * dratdumdd(irfeag)
     2            + y(ini56) * dratdumdd(irniga)


c..(a,p)(p,g) and (g,p)(p,a) reactions
      dfdy(ihe4,iden) =  dfdy(ihe4,iden)
     1            - y(ihe4)  * y(img24) * dratdumdd(irmgap) * (1.0d0-r1)
     2            + y(isi28) * dratdumdd(irsigp) * r1
     3            - y(ihe4)  * y(isi28) * dratdumdd(irsiap) * (1.0d0-s1)
     4            + y(is32)  * dratdumdd(irsgp) * s1
     5            - y(ihe4)  * y(is32) * dratdumdd(irsap) * (1.0d0-t1)   
     6            + y(iar36) * dratdumdd(irargp) * t1 
     7            - y(ihe4)  * y(iar36) * dratdumdd(irarap) * (1.0d0-u1)
     8            + y(ica40) * dratdumdd(ircagp) * u1
     9            - y(ihe4)  * y(ica40) * dratdumdd(ircaap) * (1.0d0-v1)   
     &            + y(iti44) * dratdumdd(irtigp) * v1
     1            - y(ihe4)  * y(iti44) * dratdumdd(irtiap) * (1.0d0-w1)
     2            + y(icr48) * dratdumdd(ircrgp) * w1
     3            - y(ihe4)  * y(icr48) * dratdumdd(ircrap) * (1.0d0-x1)   
     4            + y(ife52) * dratdumdd(irfegp) * x1 
     5            - y(ihe4)  * y(ife52) * dratdumdd(irfeap) * (1.0d0-y1)
     6            + y(ini56) * dratdumdd(irnigp) * y1


      dfdy(ihe4,iden) =  dfdy(ihe4,iden)
     1            + y(ihe4)  * y(img24) * ratdum(irmgap)* r1dd
     2            + y(isi28) * ratdum(irsigp) * r1dd
     3            + y(ihe4)  * y(isi28) * ratdum(irsiap) * s1dd
     4            + y(is32)  * ratdum(irsgp) * s1dd
     5            + y(ihe4)  * y(is32) * ratdum(irsap) * t1dd   
     6            + y(iar36) * ratdum(irargp) * t1dd 
     7            + y(ihe4)  * y(iar36) * ratdum(irarap) * u1dd
     8            + y(ica40) * ratdum(ircagp) * u1dd
     9            + y(ihe4)  * y(ica40) * ratdum(ircaap) * v1dd   
     &            + y(iti44) * ratdum(irtigp) * v1dd
     1            + y(ihe4)  * y(iti44) * ratdum(irtiap) * w1dd
     2            + y(icr48) * ratdum(ircrgp) * w1dd
     3            + y(ihe4)  * y(icr48) * ratdum(ircrap) * x1dd   
     4            + y(ife52) * ratdum(irfegp) * x1dd 
     5            + y(ihe4)  * y(ife52) * ratdum(irfeap) * y1dd
     6            + y(ini56) * ratdum(irnigp) * y1dd


c..c12 reactions   
      dfdy(ic12,iden) = 
     1             -y(ic12) * y(ic12) * dratdumdd(ir1212) 
     1            - y(ic12) * y(io16) * dratdumdd(ir1216)       
     2            + sixth * y(ihe4)*y(ihe4)*y(ihe4) * dratdumdd(ir3a)
     3            - y(ic12) * dratdumdd(irg3a)             
     4            - y(ic12) * y(ihe4) * dratdumdd(ircag)
     5            + y(io16) * dratdumdd(iroga)

c..o16 reactions 
      dfdy(io16,iden) = 
     1             -y(ic12) * y(io16) * dratdumdd(ir1216)       
     1            - y(io16) * y(io16) * dratdumdd(ir1616) 
     2            + y(ic12) * y(ihe4) * dratdumdd(ircag)
     3            - y(io16) * y(ihe4) * dratdumdd(iroag)
     4            - y(io16) * dratdumdd(iroga)             
     5            + y(ine20) * dratdumdd(irnega)

c..ne20 reactions 
      dfdy(ine20,iden) =  
     1               0.5d0 * y(ic12) * y(ic12) * dratdumdd(ir1212)      
     1             + y(io16) * y(ihe4) * dratdumdd(iroag)       
     2             - y(ine20) * y(ihe4) * dratdumdd(irneag)
     3             - y(ine20) * dratdumdd(irnega)
     4             + y(img24) * dratdumdd(irmgga)            

c..mg24 reactions
      dfdy(img24,iden)  = 
     1                0.5d0 * y(ic12) * y(io16) * dratdumdd(ir1216)
     1              + y(ine20) * y(ihe4) * dratdumdd(irneag)
     2              - y(img24) * y(ihe4) * dratdumdd(irmgag) 
     4              - y(img24) * dratdumdd(irmgga) 
     5              + y(isi28) * dratdumdd(irsiga)
     6              - y(img24) *  y(ihe4) *dratdumdd(irmgap)*(1.0d0-r1)
     7              + y(isi28) * r1 * dratdumdd(irsigp)
     6              + y(img24) *  y(ihe4) *ratdum(irmgap)*r1dd
     7              + y(isi28) * r1dd * ratdum(irsigp)

c..si28 reactions
      dfdy(isi28,iden)  =  
     1                0.5d0 * y(ic12) * y(io16) * dratdumdd(ir1216) 
     1              + 0.56d0*0.5d0*y(io16)*y(io16) * dratdumdd(ir1616)
     2              + 0.34d0*0.5d0*y(io16)*y(io16)*s1*dratdumdd(ir1616)
     2              + 0.34d0*0.5d0*y(io16)*y(io16)*s1dd*ratdum(ir1616)
     3              + y(img24) * y(ihe4) * dratdumdd(irmgag) 
     4              - y(isi28) * y(ihe4) * dratdumdd(irsiag) 
     5              - y(isi28) * dratdumdd(irsiga)
     6              + y(is32)  * dratdumdd(irsga) 
     7              + y(img24) * y(ihe4) * dratdumdd(irmgap)*(1.0d0-r1)
     8              - y(isi28) * r1 * dratdumdd(irsigp) 
     7              - y(img24) * y(ihe4) * ratdum(irmgap)*r1dd
     8              - y(isi28) * r1dd * ratdum(irsigp) 
     9              - y(isi28) * y(ihe4) * dratdumdd(irsiap)*(1.0d0-s1)
     &              + y(is32)  * s1 * dratdumdd(irsgp)
     9              + y(isi28) * y(ihe4) * ratdum(irsiap)*s1dd
     &              + y(is32)  * s1dd * ratdum(irsgp)

c..s32 reactions
      dfdy(is32,iden) =  
     1              0.34d0*0.5d0*y(io16)**2*dratdumdd(ir1616)*(1.0d0-s1)
     1             - 0.34d0*0.5d0*y(io16)**2 * ratdum(ir1616)*s1dd
     1             + 0.1d0*0.5d0 * y(io16) * y(io16) * dratdumdd(ir1616) 
     2             +  y(isi28) * y(ihe4) * dratdumdd(irsiag) 
     3             - y(is32) * y(ihe4) * dratdumdd(irsag) 
     4             - y(is32) * dratdumdd(irsga) 
     5             + y(iar36) * dratdumdd(irarga)
     6             + y(isi28) * y(ihe4) * dratdumdd(irsiap) * (1.0d0-s1)
     7             - y(is32) * s1 * dratdumdd(irsgp)
     6             - y(isi28) * y(ihe4) * ratdum(irsiap) * s1dd
     7             - y(is32) * s1dd * ratdum(irsgp)
     8             - y(is32) * y(ihe4) * dratdumdd(irsap) * (1.0d0-t1)
     9             + y(iar36) * t1 * dratdumdd(irargp)
     8             - y(is32) * y(ihe4) * ratdum(irsap) * t1dd
     9             + y(iar36) * t1dd * ratdum(irargp)


c..ar36 reactions
      dfdy(iar36,iden) =  
     1               y(is32)  * y(ihe4) * dratdumdd(irsag)
     1             - y(iar36) * y(ihe4) * dratdumdd(irarag)
     2             - y(iar36) * dratdumdd(irarga) 
     3             + y(ica40) * dratdumdd(ircaga)
     4             + y(is32)  * y(ihe4) * dratdumdd(irsap) * (1.0d0-t1) 
     5             - y(iar36) * t1 * dratdumdd(irargp)
     4             - y(is32)  * y(ihe4) * ratdum(irsap) * t1dd 
     5             - y(iar36) * t1dd * ratdum(irargp)
     6             - y(iar36) * y(ihe4) * dratdumdd(irarap) * (1.0d0-u1)
     7             + y(ica40) * dratdumdd(ircagp) * u1
     6             + y(iar36) * y(ihe4) * ratdum(irarap) * u1dd
     7             + y(ica40) * ratdum(ircagp) * u1dd


c..ca40 reactions
      dfdy(ica40,iden) =  
     1               y(iar36) * y(ihe4) * dratdumdd(irarag)
     1             - y(ica40) * y(ihe4) * dratdumdd(ircaag)
     2             - y(ica40) * dratdumdd(ircaga) 
     3             + y(iti44) * dratdumdd(irtiga) 
     4             + y(iar36) * y(ihe4) * dratdumdd(irarap) * (1.0d0-u1)
     5             - y(ica40) * dratdumdd(ircagp) * u1
     4             - y(iar36) * y(ihe4) * ratdum(irarap) * u1dd
     5             - y(ica40) * ratdum(ircagp) * u1dd
     6             - y(ica40) * y(ihe4) * dratdumdd(ircaap) * (1.0d0-v1)
     7             + y(iti44) * dratdumdd(irtigp) * v1
     6             + y(ica40) * y(ihe4) * ratdum(ircaap) * v1dd
     7             + y(iti44) * ratdum(irtigp) * v1dd

c..ti44 reactions
      dfdy(iti44,iden) =  
     1               y(ica40) * y(ihe4) * dratdumdd(ircaag)
     1             - y(iti44) * y(ihe4) * dratdumdd(irtiag)
     2             - y(iti44) * dratdumdd(irtiga) 
     3             + y(icr48) * dratdumdd(ircrga)
     4             + y(ica40) * y(ihe4) * dratdumdd(ircaap)*(1.0d0-v1)
     5             - y(iti44) * v1 * dratdumdd(irtigp)
     4             - y(ica40) * y(ihe4) * ratdum(ircaap)*v1dd
     5             - y(iti44) * v1dd * ratdum(irtigp)
     6             - y(iti44) * y(ihe4) * dratdumdd(irtiap) * (1.0d0-w1)
     7             + y(icr48) * w1 * dratdumdd(ircrgp)
     6             + y(iti44) * y(ihe4) * ratdum(irtiap) * w1dd
     7             + y(icr48) * w1dd * ratdum(ircrgp)

c..cr48 reactions
      dfdy(icr48,iden) =  
     1               y(iti44) * y(ihe4) * dratdumdd(irtiag)
     1             - y(icr48) * y(ihe4) * dratdumdd(ircrag)
     2             - y(icr48) * dratdumdd(ircrga) 
     3             + y(ife52) * dratdumdd(irfega) 
     4             + y(iti44) * y(ihe4) * dratdumdd(irtiap)*(1.0d0-w1)
     5             - y(icr48) * w1 * dratdumdd(ircrgp)
     4             - y(iti44) * y(ihe4) * ratdum(irtiap)*w1dd
     5             - y(icr48) * w1dd * ratdum(ircrgp)
     6             - y(icr48) * y(ihe4) * dratdumdd(ircrap) * (1.0d0-x1)
     7             + y(ife52) * x1 * dratdumdd(irfegp)
     6             + y(icr48) * y(ihe4) * ratdum(ircrap) * x1dd
     7             + y(ife52) * x1dd * ratdum(irfegp)

c..fe52 reactions
      dfdy(ife52,iden) =  
     1               y(icr48) * y(ihe4) * dratdumdd(ircrag)
     1             - y(ife52) * y(ihe4) * dratdumdd(irfeag)
     2             - y(ife52) * dratdumdd(irfega)
     3             + y(ini56) * dratdumdd(irniga)
     4             + y(icr48) * y(ihe4) * dratdumdd(ircrap) * (1.0d0-x1) 
     5             - y(ife52) * x1 * dratdumdd(irfegp)
     4             - y(icr48) * y(ihe4) * ratdum(ircrap) * x1dd 
     5             - y(ife52) * x1dd * ratdum(irfegp)
     6             - y(ife52) * y(ihe4) * dratdumdd(irfeap) * (1.0d0-y1)
     7             + y(ini56) * y1 * dratdumdd(irnigp)
     6             + y(ife52) * y(ihe4) * ratdum(irfeap) * y1dd
     7             + y(ini56) * y1dd * ratdum(irnigp)

c..ni56 reactions
      dfdy(ini56,iden) =  
     1               y(ife52) * y(ihe4) * dratdumdd(irfeag)
     1             - y(ini56) * dratdumdd(irniga)
     2             + y(ife52) * y(ihe4) * dratdumdd(irfeap) * (1.0d0-y1)  
     3             - y(ini56) * y1 * dratdumdd(irnigp)
     2             - y(ife52) * y(ihe4) * ratdum(irfeap) * y1dd  
     3             - y(ini56) * y1dd * ratdum(irnigp)





c..append the energy generation rate jacobian elements
      do j=1,ionmax
       do i=1,ionmax
        dfdy(iener,j) = dfdy(iener,j) + dfdy(i,j)*bion(i)
       enddo
       dfdy(iener,j) = dfdy(iener,j)*conv 
       dfdy(iener,itemp) = dfdy(iener,itemp) + dfdy(j,itemp)*bion(j)
       dfdy(iener,iden)  = dfdy(iener,iden) + dfdy(j,iden)*bion(j)
      enddo
      dfdy(iener,itemp) = dfdy(iener,itemp) * conv
      dfdy(iener,iden)  = dfdy(iener,iden) * conv




c..account for the neutrino losses
      call sneut5(btemp,bden,abar,zbar,
     1            sneut,snudt,snudd,snuda,snudz)

      do j=1,ionmax
       dfdy(iener,j) = dfdy(iener,j)  
     1               - (-abar*abar*snuda + (zion(j) - zbar)*abar*snudz)
      enddo
      dfdy(iener,itemp) = dfdy(iener,itemp) - snudt
      dfdy(iener,iden)  = dfdy(iener,iden)  - snudd



c..for hydrostatic or one step burns all the temperature and density 
c..jacobian elements are zero, so there is nothing to do.


c..adiabatic expansion
      if (expansion) then
       taud = 446.0d0/sqrt(den0) 
       taut = 3.0d0 * taud  
       dfdy(itemp,itemp) = -psi/taut
       dfdy(iden,iden)   = -psi/taud


c..for self-heating, we need the specific heat at constant volume
      else if (self_heat) then


c..call an eos 
       temp_row(1) = btemp
       den_row(1)  = bden
       abar_row(1) = abar
       zbar_row(1) = zbar
       jlo_eos = 1
       jhi_eos = 1

       call helmeos

       zz = 1.0d0/cv_row(1)


c..d(itemp)/d(yi)
      do j=1,ionmax
       dfdy(itemp,j) = zz*dfdy(iener,j)
      enddo

      xx = dea_row(1)*abar*abar*zz
      do j=1,ionmax
       do i=1,ionmax
        dfdy(itemp,j) = dfdy(itemp,j) - dfdy(i,j)*xx
       enddo
      enddo

      xx = dez_row(1)*abar*zz
      do j=1,ionmax
       do i=1,ionmax
        dfdy(itemp,j) = dfdy(itemp,j) - dfdy(i,j)*(zion(i)-zbar)*xx
       enddo
      enddo


c..d(itemp)/d(temp)
       sum1 = 0.0d0
       do i=1,ionmax
        sum1 = sum1 - dfdy(i,itemp)
       enddo
       sum1 = sum1 * dea_row(1)*abar*abar

       sum2 = 0.0d0
       do i=1,ionmax
        sum2 = sum2 + (zion(i) - zbar)*dfdy(i,itemp)
       enddo
       sum2 = sum2 * dez_row(1)*abar

       dfdy(itemp,itemp) = zz*(dfdy(iener,itemp) - sum1 - sum2)


c..d(itemp)/d(den)
       sum1 = 0.0d0
       do i=1,ionmax
        sum1 = sum1 - dfdy(i,iden)
       enddo
       sum1 = sum1 * dea_row(1)*abar*abar

       sum2 = 0.0d0
       do i=1,ionmax
        sum2 = sum2 + (zion(i) - zbar)*dfdy(i,iden)
       enddo
       sum2 = sum2 * dez_row(1)*abar

       dfdy(itemp,iden) = zz*(dfdy(iener,iden) - sum1 - sum2)

      end if



c..shut down the temperature and density derivatives
c      do i=1,ionmax
c       dfdy(i,itemp) = 0.0d0
c       dfdy(i,iden) = 0.0d0
c      enddo 
c      dfdy(iener,itemp) = 0.0d0
c      dfdy(iener,iden) = 0.0d0


      return
      end   





      subroutine aprox13rat
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

c..this routine generates nuclear reaction rates for the aprox13 network.

c..declare  
      integer          i
      double precision rrate,drratedt,drratedd


c..zero the rates
      do i=1,nrat
       ratraw(i) = 0.0d0
      enddo
      do i=1,nrat
       dratrawdt(i) = 0.0d0
      enddo
      do i=1,nrat
       dratrawdd(i) = 0.0d0
      enddo

      if (btemp .lt. 1.0e6) return


c..get the temperature factors
      call tfactors(btemp)     



c..c12(a,g)o16
      call rate_c12ag(btemp,bden,
     1     ratraw(ircag),dratrawdt(ircag),dratrawdd(ircag),
     2     ratraw(iroga),dratrawdt(iroga),dratrawdd(iroga))

c..triple alpha to c12 
      call rate_tripalf(btemp,bden,
     1     ratraw(ir3a),dratrawdt(ir3a),dratrawdd(ir3a),
     2     ratraw(irg3a),dratrawdt(irg3a),dratrawdd(irg3a))

c..c12 + c12 
      call rate_c12c12(btemp,bden,
     1     ratraw(ir1212),dratrawdt(ir1212),dratrawdd(ir1212),
     2     rrate,drratedt,drratedd)

c..c12 + o16 
      call rate_c12o16(btemp,bden,
     1     ratraw(ir1216),dratrawdt(ir1216),dratrawdd(ir1216),
     2     rrate,drratedt,drratedd)

c..o16 + o16 
      call rate_o16o16(btemp,bden,
     1     ratraw(ir1616),dratrawdt(ir1616),dratrawdd(ir1616),
     2     rrate,drratedt,drratedd)

c..o16(a,g)ne20
      call rate_o16ag(btemp,bden,
     1     ratraw(iroag),dratrawdt(iroag),dratrawdd(iroag),
     2     ratraw(irnega),dratrawdt(irnega),dratrawdd(irnega))

c..ne20(a,g)mg24 
      call rate_ne20ag(btemp,bden,
     1     ratraw(irneag),dratrawdt(irneag),dratrawdd(irneag),
     2     ratraw(irmgga),dratrawdt(irmgga),dratrawdd(irmgga))

c..mg24(a,g)si28 
      call rate_mg24ag(btemp,bden,
     1     ratraw(irmgag),dratrawdt(irmgag),dratrawdd(irmgag),
     2     ratraw(irsiga),dratrawdt(irsiga),dratrawdd(irsiga))

c..mg24(a,p)al27 
      call rate_mg24ap(btemp,bden,
     1     ratraw(irmgap),dratrawdt(irmgap),dratrawdd(irmgap),
     2     ratraw(iralpa),dratrawdt(iralpa),dratrawdd(iralpa))

c..al27(p,g)si28
      call rate_al27pg(btemp,bden,
     1     ratraw(iralpg),dratrawdt(iralpg),dratrawdd(iralpg),
     2     ratraw(irsigp),dratrawdt(irsigp),dratrawdd(irsigp))

c..si28(a,g)s32 
      call rate_si28ag(btemp,bden,
     1     ratraw(irsiag),dratrawdt(irsiag),dratrawdd(irsiag),
     2     ratraw(irsga),dratrawdt(irsga),dratrawdd(irsga))

c..si28(a,p)p31 
      call rate_si28ap(btemp,bden,
     1     ratraw(irsiap),dratrawdt(irsiap),dratrawdd(irsiap),
     2     ratraw(irppa),dratrawdt(irppa),dratrawdd(irppa))

c..p31(p,g)s32 
      call rate_p31pg(btemp,bden,
     1     ratraw(irppg),dratrawdt(irppg),dratrawdd(irppg),
     2     ratraw(irsgp),dratrawdt(irsgp),dratrawdd(irsgp))

c..s32(a,g)ar36 
      call rate_s32ag(btemp,bden,
     1     ratraw(irsag),dratrawdt(irsag),dratrawdd(irsag),
     2     ratraw(irarga),dratrawdt(irarga),dratrawdd(irarga))

c..s32(a,p)cl35 
      call rate_s32ap(btemp,bden,
     1     ratraw(irsap),dratrawdt(irsap),dratrawdd(irsap),
     2     ratraw(irclpa),dratrawdt(irclpa),dratrawdd(irclpa))

c..cl35(p,g)ar36
      call rate_cl35pg(btemp,bden,
     1     ratraw(irclpg),dratrawdt(irclpg),dratrawdd(irclpg),
     2     ratraw(irargp),dratrawdt(irargp),dratrawdd(irargp))

c..ar36(a,g)ca40 
      call rate_ar36ag(btemp,bden,
     1     ratraw(irarag),dratrawdt(irarag),dratrawdd(irarag),
     2     ratraw(ircaga),dratrawdt(ircaga),dratrawdd(ircaga))

c..ar36(a,p)k39
      call rate_ar36ap(btemp,bden,
     1     ratraw(irarap),dratrawdt(irarap),dratrawdd(irarap),
     2     ratraw(irkpa),dratrawdt(irkpa),dratrawdd(irkpa))

c..k39(p,g)ca40 
      call rate_k39pg(btemp,bden,
     1     ratraw(irkpg),dratrawdt(irkpg),dratrawdd(irkpg),
     2     ratraw(ircagp),dratrawdt(ircagp),dratrawdd(ircagp))

c..ca40(a,g)ti44 
      call rate_ca40ag(btemp,bden,
     1     ratraw(ircaag),dratrawdt(ircaag),dratrawdd(ircaag),
     2     ratraw(irtiga),dratrawdt(irtiga),dratrawdd(irtiga))

c..ca40(a,p)sc43 
      call rate_ca40ap(btemp,bden,
     1     ratraw(ircaap),dratrawdt(ircaap),dratrawdd(ircaap),
     2     ratraw(irscpa),dratrawdt(irscpa),dratrawdd(irscpa))

c..sc43(p,g)ti44 
      call rate_sc43pg(btemp,bden,
     1     ratraw(irscpg),dratrawdt(irscpg),dratrawdd(irscpg),
     2     ratraw(irtigp),dratrawdt(irtigp),dratrawdd(irtigp))

c..ti44(a,g)cr48
      call rate_ti44ag(btemp,bden,
     1     ratraw(irtiag),dratrawdt(irtiag),dratrawdd(irtiag),
     2     ratraw(ircrga),dratrawdt(ircrga),dratrawdd(ircrga))

c..ti44(a,p)v47 
      call rate_ti44ap(btemp,bden,
     1     ratraw(irtiap),dratrawdt(irtiap),dratrawdd(irtiap),
     2     ratraw(irvpa),dratrawdt(irvpa),dratrawdd(irvpa))

c..v47(p,g)cr48 
      call rate_v47pg(btemp,bden,
     1     ratraw(irvpg),dratrawdt(irvpg),dratrawdd(irvpg),
     2     ratraw(ircrgp),dratrawdt(ircrgp),dratrawdd(ircrgp))

c..cr48(a,g)fe52
      call rate_cr48ag(btemp,bden,
     1     ratraw(ircrag),dratrawdt(ircrag),dratrawdd(ircrag),
     2     ratraw(irfega),dratrawdt(irfega),dratrawdd(irfega))

c..cr48(a,p)mn51 
      call rate_cr48ap(btemp,bden,
     1     ratraw(ircrap),dratrawdt(ircrap),dratrawdd(ircrap),
     2     ratraw(irmnpa),dratrawdt(irmnpa),dratrawdd(irmnpa))

c..mn51(p,g)fe52 
      call rate_mn51pg(btemp,bden,
     1     ratraw(irmnpg),dratrawdt(irmnpg),dratrawdd(irmnpg),
     2     ratraw(irfegp),dratrawdt(irfegp),dratrawdd(irfegp))

c..fe52(a,g)ni56
      call rate_fe52ag(btemp,bden,
     1     ratraw(irfeag),dratrawdt(irfeag),dratrawdd(irfeag),
     2     ratraw(irniga),dratrawdt(irniga),dratrawdd(irniga))

c..fe52(a,p)co55 
      call rate_fe52ap(btemp,bden,
     1     ratraw(irfeap),dratrawdt(irfeap),dratrawdd(irfeap),
     2     ratraw(ircopa),dratrawdt(ircopa),dratrawdd(ircopa))

c..co55(p,g)ni56 
      call rate_co55pg(btemp,bden,
     1     ratraw(ircopg),dratrawdt(ircopg),dratrawdd(ircopg),
     2     ratraw(irnigp),dratrawdt(irnigp),dratrawdd(irnigp))


c..write out the rates
c      write(6,133) btemp,bden
c 133  format(1x,1pe12.4)
c      do i=1,nrat
c       write(6,134) ratnam(i),ratraw(i)
c 134   format(1x,a,'  ',1pe14.6,1pe11.3,1pe14.6)
c      enddo
c      read(5,*)



c..for a strict alpha chain with only (a,g) and (g,a) reactions
c..shut down the (a,p) (p,g) and (g,p) (p,a) rates 
c      ratraw(irmgap) = 1.0d-40
c      ratraw(iralpa) = 0.0d0
c      ratraw(iralpg) = 1.0d-40
c      ratraw(irsigp) = 1.0d-40

c      ratraw(irsiap) = 1.0d-40
c      ratraw(irppa)  = 0.0d0
c      ratraw(irppg)  = 1.0d-40
c      ratraw(irsgp)  = 1.0d-40

c      ratraw(irsap)  = 1.0d-40
c      ratraw(irclpa) = 0.0d0
c      ratraw(irclpg) = 1.0d-40
c      ratraw(irargp) = 1.0d-40

c      ratraw(irarap) = 1.0d-40
c      ratraw(irkpa)  = 0.0d0
c      ratraw(irkpg)  = 1.0d-40
c      ratraw(ircagp) = 1.0d-40

c      ratraw(ircaap) = 1.0d-40
c      ratraw(irscpa) = 0.0d0
c      ratraw(irscpg) = 1.0d-40
c      ratraw(irtigp) = 1.0d-40

c      ratraw(irtiap) = 1.0d-40
c      ratraw(irvpa)  = 0.0d0
c      ratraw(irvpg)  = 1.0d-40
c      ratraw(ircrgp) = 1.0d-40

c      ratraw(ircrap) = 1.0d-40
c      ratraw(irmnpa) = 0.0d0
c      ratraw(irmnpg) = 1.0d-40
c      ratraw(irfegp) = 1.0d-40

c..shutting down this (a,p)(p,g) sequence 
c..will make aprox13 exactly like iso13
c      ratraw(irfeap) = 1.0d-40
c      ratraw(ircopa) = 0.0d0
c      ratraw(ircopg) = 1.0d-40
c      ratraw(irnigp) = 1.0d-40

      return
      end   





      subroutine aprox13tab
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

c..uses tables instead of analytical expressions to evaluate the 
c..raw reaction rates. a cubic polynomial is hardwired for speed.

      integer          i,j,imax,iat,mp,ifirst
      parameter        (mp = 4)
      double precision tlo,thi,tstp,bden_sav,btemp_sav,
     1                 x,x1,x2,x3,x4,a,b,c,d,e,f,g,h,p,q,
     2                 alfa,beta,gama,delt
      data             ifirst/0/


c..make the table
      if (ifirst .eq. 0) then
       ifirst = 1


c..set the log temperature loop limits
       thi  = 10.0d0
       tlo  = 6.0d0
       imax = int(thi-tlo)*120 + 1
       if (imax .gt. nrattab) stop 'imax too small in aprox13tab'
       tstp = (thi - tlo)/float(imax-1)


c..save the input
       btemp_sav = btemp
       bden_sav  = bden


c..form the table
       bden = 1.0d0
       do i=1,imax
        btemp = tlo + float(i-1)*tstp
        btemp = 10.0d0**(btemp)
        call aprox13rat
        ttab(i) = btemp
        do j=1,nrat
         rattab(j,i)    = ratraw(j)
         drattabdt(j,i) = dratrawdt(j)
         drattabdd(j,i) = dratrawdd(j)
        enddo
       enddo

c..restore the input
       bden  = bden_sav
       btemp = btemp_sav
      end if
 

c..normal execution starts here
c..set the density dependence vector
      dtab(ircag)  = bden 
      dtab(iroga)  = 1.0d0
      dtab(ir3a)   = bden*bden
      dtab(irg3a)  = 1.0d0
      dtab(ir1212) = bden 
      dtab(ir1216) = bden 
      dtab(ir1616) = bden 
      dtab(iroag)  = bden
      dtab(irnega) = 1.0d0
      dtab(irneag) = bden 
      dtab(irmgga) = 1.0d0
      dtab(irmgag) = bden 
      dtab(irsiga) = 1.0d0
      dtab(irmgap) = bden
      dtab(iralpa) = bden 
      dtab(iralpg) = bden 
      dtab(irsigp) = 1.0d0
      dtab(irsiag) = bden 
      dtab(irsga)  = 1.0d0
      dtab(irppa)  = bden
      dtab(irsiap) = bden
      dtab(irppg)  = bden
      dtab(irsgp)  = 1.0d0
      dtab(irsag)  = bden
      dtab(irarga) = 1.0d0
      dtab(irsap)  = bden
      dtab(irclpa) = bden
      dtab(irclpg) = bden
      dtab(irargp) = 1.0d0
      dtab(irarag) = bden
      dtab(ircaga) = 1.0d0
      dtab(irarap) = bden
      dtab(irkpa)  = bden
      dtab(irkpg)  = bden
      dtab(ircagp) = 1.0d0
      dtab(ircaag) = bden
      dtab(irtiga) = 1.0d0
      dtab(ircaap) = bden
      dtab(irscpa) = bden
      dtab(irscpg) = bden
      dtab(irtigp) = 1.0d0
      dtab(irtiag) = bden
      dtab(ircrga) = 1.0d0
      dtab(irtiap) = bden
      dtab(irvpa)  = bden
      dtab(irvpg)  = bden
      dtab(ircrgp) = 1.0d0
      dtab(ircrag) = bden
      dtab(irfega) = 1.0d0
      dtab(ircrap) = bden
      dtab(irmnpa) = bden
      dtab(irmnpg) = bden
      dtab(irfegp) = 1.0d0
      dtab(irfeag) = bden
      dtab(irniga) = 1.0d0
      dtab(irfeap) = bden
      dtab(ircopa) = bden
      dtab(ircopg) = bden
      dtab(irnigp) = 1.0d0


c..hash locate the temperature
      iat = int((log10(btemp) - tlo)/tstp) + 1
      iat = max(1,min(iat - mp/2 + 1,imax - mp + 1))

c..setup the lagrange interpolation coefficients for a cubic
      x  = btemp
      x1 = ttab(iat)
      x2 = ttab(iat+1)
      x3 = ttab(iat+2)
      x4 = ttab(iat+3)
      a  = x - x1
      b  = x - x2
      c  = x - x3
      d  = x - x4
      e  = x1 - x2
      f  = x1 - x3
      g  = x1 - x4
      h  = x2 - x3
      p  = x2 - x4
      q  = x3 - x4
      alfa =  b*c*d/(e*f*g)
      beta = -a*c*d/(e*h*p)
      gama =  a*b*d/(f*h*q)
      delt = -a*b*c/(g*p*q)

c..crank off the raw reaction rates
      do j=1,nrat
       ratraw(j) = (alfa*rattab(j,iat)
     1            + beta*rattab(j,iat+1)
     2            + gama*rattab(j,iat+2) 
     3            + delt*rattab(j,iat+3) 
     4              ) * dtab(j)

       dratrawdt(j) = (alfa*drattabdt(j,iat)
     1               + beta*drattabdt(j,iat+1)
     2               + gama*drattabdt(j,iat+2) 
     3               + delt*drattabdt(j,iat+3) 
     4                 ) * dtab(j)

       dratrawdd(j) =  alfa*drattabdd(j,iat)
     1               + beta*drattabdd(j,iat+1)
     2               + gama*drattabdd(j,iat+2) 
     3               + delt*drattabdd(j,iat+3) 

      enddo

c..hand finish the three body reactions
      dratrawdd(ir3a) = bden * dratrawdd(ir3a)

      return
      end






      subroutine screen_aprox13(y)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

c..this routine computes the screening factors
c..and applies them to the raw reaction rates,
c..producing the final reaction rates used by the
c..right hand sides and jacobian matrix elements

c..this routine assumes screen_on = 1 or = 0 has been set at a higher 
c..level presumably in the top level driver


c..declare
      integer          i,j,k,jscr,init
      double precision y(*),sc1a,sc1adt,sc1add,sc2a,sc2adt,sc2add,
     1                 sc3a,sc3adt,sc3add,
     2                 abar,zbar,z2bar,ytot1,zbarxx,z2barxx
      data             init/1/


c..roll all of them 
      do i=1,nrat
       ratdum(i)    = ratraw(i)
       dratdumdt(i) = dratrawdt(i)
       dratdumdd(i) = dratrawdd(i)
       scfac(i)     = 1.0d0
       dscfacdt(i)  = 0.0d0
       dscfacdt(i)  = 0.0d0
      end do

c..if screening is off
      if (screen_on .eq. 0) return


c..screening is on
c..with the passed composition, compute abar,zbar and other variables
      zbarxx  = 0.0d0
      z2barxx = 0.0d0
      ytot1   = 0.0d0
      do i=1,ionmax
       ytot1    = ytot1 + y(i)
       zbarxx   = zbarxx + zion(i) * y(i)
       z2barxx  = z2barxx + zion(i) * zion(i) * y(i)
      enddo
      abar   = 1.0d0/ytot1
      zbar   = zbarxx * abar
      z2bar  = z2barxx * abar



c..first the always fun triple alpha and its inverse
      jscr = 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             zion(ihe4),aion(ihe4),zion(ihe4),aion(ihe4),
     2             jscr,init,sc1a,sc1adt,sc1add)

      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             zion(ihe4),aion(ihe4),4.0d0,8.0d0,
     2             jscr,init,sc2a,sc2adt,sc2add)

      sc3a   = sc1a * sc2a                     
      sc3adt = sc1adt*sc2a + sc1a*sc2adt       
      sc3add = sc1add*sc2a + sc1a*sc2add       

      ratdum(ir3a)    = ratraw(ir3a) * sc3a
      dratdumdt(ir3a) = dratrawdt(ir3a)*sc3a + ratraw(ir3a)*sc3adt
      dratdumdd(ir3a) = dratrawdd(ir3a)*sc3a + ratraw(ir3a)*sc3add

      scfac(ir3a)     = sc3a
      dscfacdt(ir3a)  = sc3adt
      dscfacdd(ir3a)  = sc3add


c..c12 to o16 
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             zion(ic12),aion(ic12),zion(ihe4),aion(ihe4),
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ircag)     = ratraw(ircag) * sc1a
      dratdumdt(ircag)  = dratrawdt(ircag)*sc1a + ratraw(ircag)*sc1adt
      dratdumdd(ircag)  = dratrawdd(ircag)*sc1a + ratraw(ircag)*sc1add

      scfac(ircag)      = sc1a
      dscfacdt(ircag)   = sc1adt
      dscfacdt(ircag)   = sc1add


c..c12 + c12
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             zion(ic12),aion(ic12),zion(ic12),aion(ic12),
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ir1212)    = ratraw(ir1212) * sc1a
      dratdumdt(ir1212) = dratrawdt(ir1212)*sc1a + ratraw(ir1212)*sc1adt
      dratdumdd(ir1212) = dratrawdd(ir1212)*sc1a + ratraw(ir1212)*sc1add

      scfac(ir1212)     = sc1a
      dscfacdt(ir1212)  = sc1adt
      dscfacdd(ir1212)  = sc1add



c..c12 + o16
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             zion(ic12),aion(ic12),zion(io16),aion(io16),
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ir1216)    = ratraw(ir1216) * sc1a
      dratdumdt(ir1216) = dratrawdt(ir1216)*sc1a + ratraw(ir1216)*sc1adt
      dratdumdd(ir1216) = dratrawdd(ir1216)*sc1a + ratraw(ir1216)*sc1add

      scfac(ir1216)     = sc1a
      dscfacdt(ir1216)  = sc1adt
      dscfacdd(ir1216)  = sc1add



c..o16 + o16
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             zion(io16),aion(io16),zion(io16),aion(io16),
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ir1616)    = ratraw(ir1616) * sc1a
      dratdumdt(ir1616) = dratrawdt(ir1616)*sc1a + ratraw(ir1616)*sc1adt
      dratdumdd(ir1616) = dratrawdd(ir1616)*sc1a + ratraw(ir1616)*sc1add

      scfac(ir1616)     = sc1a
      dscfacdt(ir1616)  = sc1adt
      dscfacdd(ir1616)  = sc1add



c..o16 to ne20
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             zion(io16),aion(io16),zion(ihe4),aion(ihe4),
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(iroag)    = ratraw(iroag) * sc1a 
      dratdumdt(iroag) = dratrawdt(iroag)*sc1a + ratraw(iroag)*sc1adt 
      dratdumdd(iroag) = dratrawdd(iroag)*sc1a + ratraw(iroag)*sc1add 

      scfac(iroag)     = sc1a
      dscfacdt(iroag)  = sc1adt
      dscfacdd(iroag)  = sc1add



c..ne20 to mg24
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             zion(ine20),aion(ine20),zion(ihe4),aion(ihe4),
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irneag)    = ratraw(irneag) * sc1a
      dratdumdt(irneag) = dratrawdt(irneag)*sc1a + ratraw(irneag)*sc1adt
      dratdumdd(irneag) = dratrawdd(irneag)*sc1a + ratraw(irneag)*sc1add

      scfac(irneag)     = sc1a
      dscfacdt(irneag)  = sc1adt
      dscfacdd(irneag)  = sc1add


c..mg24 to si28
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             zion(img24),aion(img24),zion(ihe4),aion(ihe4),
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irmgag)    = ratraw(irmgag) * sc1a
      dratdumdt(irmgag) = dratrawdt(irmgag)*sc1a + ratraw(irmgag)*sc1adt
      dratdumdd(irmgag) = dratrawdd(irmgag)*sc1a + ratraw(irmgag)*sc1add

      scfac(irmgag)     = sc1a
      dscfacdt(irmgag)  = sc1adt
      dscfacdd(irmgag)  = sc1add

      ratdum(irmgap)    = ratraw(irmgap) * sc1a
      dratdumdt(irmgap) = dratrawdt(irmgap)*sc1a + ratraw(irmgap)*sc1adt
      dratdumdd(irmgap) = dratrawdd(irmgap)*sc1a + ratraw(irmgap)*sc1add

      scfac(irmgap)     = sc1a
      dscfacdt(irmgap)  = sc1adt
      dscfacdd(irmgap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             13.0d0,27.0d0,1.0d0,1.0d0,
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(iralpa)    = ratraw(iralpa) * sc1a
      dratdumdt(iralpa) = dratrawdt(iralpa)*sc1a + ratraw(iralpa)*sc1adt
      dratdumdd(iralpa) = dratrawdd(iralpa)*sc1a + ratraw(iralpa)*sc1add

      scfac(iralpa)     = sc1a
      dscfacdt(iralpa)  = sc1adt
      dscfacdd(iralpa)  = sc1add

      ratdum(iralpg)    = ratraw(iralpg) * sc1a
      dratdumdt(iralpg) = dratrawdt(iralpg)*sc1a + ratraw(iralpg)*sc1adt
      dratdumdd(iralpg) = dratrawdd(iralpg)*sc1a + ratraw(iralpg)*sc1add

      scfac(iralpg)     = sc1a
      dscfacdt(iralpg)  = sc1adt
      dscfacdd(iralpg)  = sc1add



c..si28 to s32
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             zion(isi28),aion(isi28),zion(ihe4),aion(ihe4),
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irsiag)    = ratraw(irsiag) * sc1a
      dratdumdt(irsiag) = dratrawdt(irsiag)*sc1a + ratraw(irsiag)*sc1adt
      dratdumdd(irsiag) = dratrawdd(irsiag)*sc1a + ratraw(irsiag)*sc1add

      scfac(irsiag)     = sc1a
      dscfacdt(irsiag)  = sc1adt
      dscfacdd(irsiag)  = sc1add


      ratdum(irsiap)    = ratraw(irsiap) * sc1a
      dratdumdt(irsiap) = dratrawdt(irsiap)*sc1a + ratraw(irsiap)*sc1adt
      dratdumdd(irsiap) = dratrawdd(irsiap)*sc1a + ratraw(irsiap)*sc1add

      scfac(irsiap)     = sc1a
      dscfacdt(irsiap)  = sc1adt
      dscfacdd(irsiap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             15.0d0,31.0d0,1.0d0,1.0d0,
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irppa)     = ratraw(irppa) * sc1a 
      dratdumdt(irppa)  = dratrawdt(irppa)*sc1a  + ratraw(irppa)*sc1adt 
      dratdumdd(irppa)  = dratrawdd(irppa)*sc1a  + ratraw(irppa)*sc1add 

      scfac(irppa)      = sc1a
      dscfacdt(irppa)   = sc1adt
      dscfacdd(irppa)   = sc1add

      ratdum(irppg)     = ratraw(irppg) * sc1a
      dratdumdt(irppg)  = dratrawdt(irppg)*sc1a + ratraw(irppg)*sc1adt
      dratdumdd(irppg)  = dratrawdd(irppg)*sc1a + ratraw(irppg)*sc1add

      scfac(irppg)      = sc1a
      dscfacdt(irppg)   = sc1adt
      dscfacdd(irppg)   = sc1add



c..s32 to ar36
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             zion(is32),aion(is32),zion(ihe4),aion(ihe4),
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irsag)     = ratraw(irsag) * sc1a 
      dratdumdt(irsag)  = dratrawdt(irsag)*sc1a + ratraw(irsag)*sc1adt 
      dratdumdd(irsag)  = dratrawdd(irsag)*sc1a + ratraw(irsag)*sc1add 

      scfac(irsag)      = sc1a
      dscfacdt(irsag)   = sc1adt
      dscfacdd(irsag)   = sc1add

      ratdum(irsap)     = ratraw(irsap) * sc1a
      dratdumdt(irsap)  = dratrawdt(irsap)*sc1a + ratraw(irsap)*sc1adt
      dratdumdd(irsap)  = dratrawdd(irsap)*sc1a + ratraw(irsap)*sc1add

      scfac(irsap)      = sc1a
      dscfacdt(irsap)   = sc1adt
      dscfacdd(irsap)   = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             17.0d0,35.0d0,1.0d0,1.0d0,
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irclpa)    = ratraw(irclpa) * sc1a
      dratdumdt(irclpa) = dratrawdt(irclpa)*sc1a + ratraw(irclpa)*sc1adt
      dratdumdd(irclpa) = dratrawdd(irclpa)*sc1a + ratraw(irclpa)*sc1add

      scfac(irclpa)     = sc1a
      dscfacdt(irclpa)  = sc1adt
      dscfacdt(irclpa)  = sc1add

      ratdum(irclpg)    = ratraw(irclpg) * sc1a
      dratdumdt(irclpg) = dratrawdt(irclpg)*sc1a + ratraw(irclpg)*sc1adt
      dratdumdd(irclpg) = dratrawdd(irclpg)*sc1a + ratraw(irclpg)*sc1add

      scfac(irclpg)     = sc1a
      dscfacdt(irclpg)  = sc1adt
      dscfacdd(irclpg)  = sc1add



c..ar36 to ca40
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             zion(iar36),aion(iar36),zion(ihe4),aion(ihe4),
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irarag)    = ratraw(irarag) * sc1a
      dratdumdt(irarag) = dratrawdt(irarag)*sc1a + ratraw(irarag)*sc1adt
      dratdumdd(irarag) = dratrawdd(irarag)*sc1a + ratraw(irarag)*sc1add

      scfac(irarag)     = sc1a
      dscfacdt(irarag)  = sc1adt
      dscfacdd(irarag)  = sc1add

      ratdum(irarap)    = ratraw(irarap) * sc1a
      dratdumdt(irarap) = dratrawdt(irarap)*sc1a + ratraw(irarap)*sc1adt
      dratdumdd(irarap) = dratrawdd(irarap)*sc1a + ratraw(irarap)*sc1add

      scfac(irarap)     = sc1a
      dscfacdt(irarap)  = sc1adt
      dscfacdd(irarap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             19.0d0,40.0d0,1.0d0,1.0d0,
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irkpa)     = ratraw(irkpa) * sc1a 
      dratdumdt(irkpa)  = dratrawdt(irkpa)*sc1a  + ratraw(irkpa)*sc1adt
      dratdumdd(irkpa)  = dratrawdd(irkpa)*sc1a  + ratraw(irkpa)*sc1add

      scfac(irkpa)      = sc1a
      dscfacdt(irkpa)   = sc1adt
      dscfacdd(irkpa)   = sc1add

      ratdum(irkpg)     = ratraw(irkpg) * sc1a 
      dratdumdt(irkpg)  = dratrawdt(irkpg)*sc1a  + ratraw(irkpg)*sc1adt
      dratdumdd(irkpg)  = dratrawdd(irkpg)*sc1a  + ratraw(irkpg)*sc1add

      scfac(irkpg)      = sc1a
      dscfacdt(irkpg)   = sc1adt
      dscfacdd(irkpg)   = sc1add



c..ca40 to ti44
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             zion(ica40),aion(ica40),zion(ihe4),aion(ihe4),
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ircaag)    = ratraw(ircaag) * sc1a
      dratdumdt(ircaag) = dratrawdt(ircaag)*sc1a + ratraw(ircaag)*sc1adt
      dratdumdd(ircaag) = dratrawdd(ircaag)*sc1a + ratraw(ircaag)*sc1add

      scfac(ircaag)     = sc1a
      dscfacdt(ircaag)  = sc1adt
      dscfacdd(ircaag)  = sc1add

      ratdum(ircaap)    = ratraw(ircaap) * sc1a
      dratdumdt(ircaap) = dratrawdt(ircaap)*sc1a + ratraw(ircaap)*sc1adt
      dratdumdd(ircaap) = dratrawdd(ircaap)*sc1a + ratraw(ircaap)*sc1add

      scfac(ircaap)     = sc1a
      dscfacdt(ircaap)  = sc1adt
      dscfacdd(ircaap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             21.0d0,45.0d0,1.0d0,1.0d0,
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irscpa)    = ratraw(irscpa) * sc1a
      dratdumdt(irscpa) = dratrawdt(irscpa)*sc1a + ratraw(irscpa)*sc1adt
      dratdumdd(irscpa) = dratrawdd(irscpa)*sc1a + ratraw(irscpa)*sc1add

      scfac(irscpa)     = sc1a
      dscfacdt(irscpa)  = sc1adt
      dscfacdd(irscpa)  = sc1add

      ratdum(irscpg)    = ratraw(irscpg) * sc1a
      dratdumdt(irscpg) = dratrawdt(irscpg)*sc1a + ratraw(irscpg)*sc1adt
      dratdumdd(irscpg) = dratrawdd(irscpg)*sc1a + ratraw(irscpg)*sc1add

      scfac(irscpg)     = sc1a
      dscfacdt(irscpg)  = sc1adt
      dscfacdd(irscpg)  = sc1add



c..ti44 to cr48
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             zion(iti44),aion(iti44),zion(ihe4),aion(ihe4),
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irtiag)    = ratraw(irtiag) * sc1a
      dratdumdt(irtiag) = dratrawdt(irtiag)*sc1a + ratraw(irtiag)*sc1adt
      dratdumdd(irtiag) = dratrawdd(irtiag)*sc1a + ratraw(irtiag)*sc1add

      scfac(irtiag)     = sc1a
      dscfacdt(irtiag)  = sc1adt
      dscfacdd(irtiag)  = sc1add

      ratdum(irtiap)    = ratraw(irtiap) * sc1a
      dratdumdt(irtiap) = dratrawdt(irtiap)*sc1a + ratraw(irtiap)*sc1adt
      dratdumdd(irtiap) = dratrawdd(irtiap)*sc1a + ratraw(irtiap)*sc1add

      scfac(irtiap)  = sc1a
      dscfacdt(irtiap)  = sc1adt
      dscfacdd(irtiap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             23.0d0,47.0d0,1.0d0,1.0d0,
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irvpa)     = ratraw(irvpa) * sc1a 
      dratdumdt(irvpa)  = dratrawdt(irvpa)*sc1a  + ratraw(irvpa)*sc1adt
      dratdumdd(irvpa)  = dratrawdd(irvpa)*sc1a  + ratraw(irvpa)*sc1add

      scfac(irvpa)      = sc1a
      dscfacdt(irvpa)   = sc1adt
      dscfacdd(irvpa)   = sc1add

      ratdum(irvpg)     = ratraw(irvpg) * sc1a 
      dratdumdt(irvpg)  = dratrawdt(irvpg)*sc1a  + ratraw(irvpg)*sc1adt
      dratdumdd(irvpg)  = dratrawdd(irvpg)*sc1a  + ratraw(irvpg)*sc1add

      scfac(irvpg)      = sc1a
      dscfacdt(irvpg)   = sc1adt
      dscfacdd(irvpg)   = sc1add



c..cr48 to fe52
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             zion(icr48),aion(icr48),zion(ihe4),aion(ihe4),
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ircrag)    = ratraw(ircrag) * sc1a
      dratdumdt(ircrag) = dratrawdt(ircrag)*sc1a + ratraw(ircrag)*sc1adt
      dratdumdd(ircrag) = dratrawdd(ircrag)*sc1a + ratraw(ircrag)*sc1add

      scfac(ircrag)     = sc1a
      dscfacdt(ircrag)  = sc1adt
      dscfacdd(ircrag)  = sc1add

      ratdum(ircrap)    = ratraw(ircrap) * sc1a
      dratdumdt(ircrap) = dratrawdt(ircrap)*sc1a + ratraw(ircrap)*sc1adt
      dratdumdd(ircrap) = dratrawdd(ircrap)*sc1a + ratraw(ircrap)*sc1add

      scfac(ircrap)     = sc1a
      dscfacdt(ircrap)  = sc1adt
      dscfacdd(ircrap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             25.0d0,51.0d0,1.0d0,1.0d0,
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irmnpa)    = ratraw(irmnpa) * sc1a
      dratdumdt(irmnpa) = dratrawdt(irmnpa)*sc1a + ratraw(irmnpa)*sc1adt
      dratdumdd(irmnpa) = dratrawdd(irmnpa)*sc1a + ratraw(irmnpa)*sc1add

      scfac(irmnpa)     = sc1a
      dscfacdt(irmnpa)  = sc1adt
      dscfacdd(irmnpa)  = sc1add

      ratdum(irmnpg)    = ratraw(irmnpg) * sc1a
      dratdumdt(irmnpg) = dratrawdt(irmnpg)*sc1a + ratraw(irmnpg)*sc1adt
      dratdumdd(irmnpg) = dratrawdd(irmnpg)*sc1a + ratraw(irmnpg)*sc1add

      scfac(irmnpg)     = sc1a
      dscfacdt(irmnpg)  = sc1adt
      dscfacdd(irmnpg)  = sc1add


c..fe52 to ni56
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             zion(ife52),aion(ife52),zion(ihe4),aion(ihe4),
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irfeag)    = ratraw(irfeag) * sc1a
      dratdumdt(irfeag) = dratrawdt(irfeag)*sc1a + ratraw(irfeag)*sc1adt
      dratdumdd(irfeag) = dratrawdd(irfeag)*sc1a + ratraw(irfeag)*sc1add

      scfac(irfeag)     = sc1a
      dscfacdt(irfeag)  = sc1adt
      dscfacdd(irfeag)  = sc1add

      ratdum(irfeap) = ratraw(irfeap) * sc1a
      dratdumdt(irfeap) = dratrawdt(irfeap)*sc1a + ratraw(irfeap)*sc1adt
      dratdumdd(irfeap) = dratrawdd(irfeap)*sc1a + ratraw(irfeap)*sc1add

      scfac(irfeap)     = sc1a
      dscfacdt(irfeap)  = sc1adt
      dscfacdd(irfeap)  = sc1add

      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar,
     1             27.0d0,55.0d0,1.0d0,1.0d0,
     2             jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ircopa)    = ratraw(ircopa) * sc1a
      dratdumdt(ircopa) = dratrawdt(ircopa)*sc1a + ratraw(ircopa)*sc1adt
      dratdumdd(ircopa) = dratrawdd(ircopa)*sc1a + ratraw(ircopa)*sc1add

      scfac(ircopa)     = sc1a
      dscfacdt(ircopa)  = sc1adt
      dscfacdd(ircopa)  = sc1add

      ratdum(ircopg)    = ratraw(ircopg) * sc1a
      dratdumdt(ircopg) = dratrawdt(ircopg)*sc1a + ratraw(ircopg)*sc1adt
      dratdumdd(ircopg) = dratrawdd(ircopg)*sc1a + ratraw(ircopg)*sc1add

      scfac(ircopg)     = sc1a
      dscfacdt(ircopg)  = sc1adt
      dscfacdd(ircopg)  = sc1add


c..reset the screen initialization flag
      init = 0


c..debugs
c      do i=1,nrat
c       if (ratdum(i) .lt. 0.0) then
c        write(6,110) i,ratnam(i),ratraw(i),scfac(i),ratdum(i)
c 110    format(1x,i4,' ',a,' ',1p3e12.4)
c        stop 'negative rate'
c       end if
c      enddo 

c      do i=1,nrat
c       write(6,111) i,ratnam(i),ratraw(i),scfac(i),ratdum(i)
c 111   format(1x,i4,' ',a,' ',1p3e12.4)
c      enddo 
c      read(5,*)

      return
      end





      subroutine faprox13(tt,xin)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'
c..
c..this routine computes the flows of the torch network
c..input is xin, the mass fraction vector, and the time tt.
c..output through the common blocks in network.dek 
c..are the names flonam and magnitude of each reaction term flowx.
c..
c..declare
      character*80     string
      integer          i,j,k,inu,in,ip,ia
      double precision tt,xin(1),y(ionmax),denom,
     1                 r1,s1,t1,u1,v1,w1,x1,y1,a1,a2,abar,zbar


c..zero the flows and names
      nflowx = 0
      do i=1,nrat
       flonam(i) = ' '
       flofor(i) = 0.0d0
       florev(i) = 0.0d0
       flonet(i) = 0.0d0
      enddo


c..positive definite mass fractions
      do i=1,ionmax
       xin(i) = min(1.0d0,max(xin(i),1.0d-30))
      enddo

c..copy the mass fractions
      do i=1,ionmax
       xmass(i) = xin(i)
      enddo

c..get abar, zbar and a few other composition variables
      call azbar(xmass,aion,zion,ionmax,
     1           ymass,abar,zbar)


c..copy the molar abundances
      do i=1,ionmax
       y(i) = ymass(i)
      enddo


c..get the reaction rates
      if (use_tables .eq. 1) then
       call aprox13tab
      else
       call aprox13rat
      end if



c..do the screening here because the corrections depend on the composition

      call screen_aprox13(y)




c..branching ratios for various dummy proton links
      r1     = 0.0d0
      denom  = ratdum(iralpa) + ratdum(iralpg)
      if (denom .ne. 0.0) r1 = ratdum(iralpa)/denom

      s1     = 0.0d0
      denom  = ratdum(irppa) + ratdum(irppg)
      if (denom .ne. 0.0) s1 = ratdum(irppa)/denom

      t1     = 0.0d0
      denom  = ratdum(irclpa) + ratdum(irclpg)
      if (denom .ne. 0.0) t1 = ratdum(irclpa)/denom

      u1     = 0.0d0
      denom  = ratdum(irkpa) + ratdum(irkpg)
      if (denom .ne. 0.0) u1 = ratdum(irkpa)/denom

      v1     = 0.0d0
      denom  = ratdum(irscpa) + ratdum(irscpg)
      if (denom .ne. 0.0) v1 = ratdum(irscpa)/denom

      w1     = 0.0d0
      denom  = ratdum(irvpa) + ratdum(irvpg)
      if (denom .ne. 0.0) w1 = ratdum(irvpa)/denom

      x1     = 0.0d0
      denom  = ratdum(irmnpa) + ratdum(irmnpg)
      if (denom .ne. 0.0) x1 = ratdum(irmnpa)/denom

      y1     = 0.0d0
      denom  = ratdum(ircopa) + ratdum(ircopg)
      if (denom .ne. 0.0) y1 = ratdum(ircopa)/denom





c..set up the system of odes: 
c..he4 reactions
c..heavy ion reactions
      a1 = y(ic12) * y(ic12) * ratdum(ir1212)      
      nflowx = nflowx + 1
      string = '2'//ionam(ic12)//'(g,a)'//ionam(ine20)
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = 0.0d0
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ic12) * y(io16) * ratdum(ir1216)
      nflowx = nflowx + 1
      string = 'c12+o16(g,a)'//ionam(img24)
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = 0.0d0
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(io16) * y(io16) * ratdum(ir1616) 
      nflowx = nflowx + 1
      string = '2'//ionam(io16)//'(g,a)'//ionam(isi28)
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = 0.0d0
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = s1 * y(io16) * y(io16) * ratdum(ir1616) 
      nflowx = nflowx + 1
      string = 'p31(p,a)'//ionam(isi28)
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = 0.0d0
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)


c..(a,g) and (g,a) reactions
      a1 = y(ihe4) * y(ihe4) * y(ihe4) * ratdum(ir3a) 
      a2 = y(ic12) * ratdum(irg3a)       
      nflowx = nflowx + 1
      string = 'aa(a,g)c12'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(ic12) * ratdum(ircag)       
      a2 = y(io16)  * ratdum(iroga)
      nflowx = nflowx + 1
      string = 'c12(a,g)o16'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(io16) * ratdum(iroag)       
      a2 = y(ine20) * ratdum(irnega)            
      nflowx = nflowx + 1
      string = 'o16(a,g)ne20'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(ine20) * ratdum(irneag) 
      a2 = y(img24) * ratdum(irmgga)
      nflowx = nflowx + 1
      string = 'ne20(a,g)mg24'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(img24)* ratdum(irmgag)
      a2 = y(isi28) * ratdum(irsiga) 
      nflowx = nflowx + 1
      string = 'mg24(a,g)si28'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(isi28)*ratdum(irsiag)
      a2 = y(is32)  * ratdum(irsga) 
      nflowx = nflowx + 1
      string = 'si28(a,g)s32'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(is32) * ratdum(irsag)
      a2 = y(iar36) * ratdum(irarga) 
      nflowx = nflowx + 1
      string = 's32(a,g)ar36'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(iar36)*ratdum(irarag)
      a2 = y(ica40) * ratdum(ircaga) 
      nflowx = nflowx + 1
      string = 'ar36(a,g)ca40'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(ica40)*ratdum(ircaag)
      a2 = y(iti44) * ratdum(irtiga) 
      nflowx = nflowx + 1
      string = 'ca40(a,g)ti44'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(iti44)*ratdum(irtiag)
      a2 = y(icr48) * ratdum(ircrga) 
      nflowx = nflowx + 1
      string = 'ti44(a,g)cr48'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(icr48)*ratdum(ircrag)
      a2 = y(ife52) * ratdum(irfega) 
      nflowx = nflowx + 1
      string = 'cr46(a,g)fe52'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(ife52) * ratdum(irfeag)
      a2 = y(ini56) * ratdum(irniga)
      nflowx = nflowx + 1
      string = 'fe52(a,g)ni56'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)


c..(a,p)(p,g) and (g,p)(p,a) reactions
      a1 = y(ihe4)  * y(img24) * ratdum(irmgap) * (1.0d0-r1)
      a2 = y(isi28) * ratdum(irsigp) * r1
      nflowx = nflowx + 1
      string = 'mg24(ap,pg)si28'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(isi28) * ratdum(irsiap) * (1.0d0-s1)
      a2 = y(is32)  * ratdum(irsgp) * s1
      nflowx = nflowx + 1
      string = 'si28(ap,pg)s32'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(is32) * ratdum(irsap) * (1.0d0-t1)   
      a2 = y(iar36) * ratdum(irargp) * t1 
      nflowx = nflowx + 1
      string = 's32(ap,pg)ar36'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(iar36) * ratdum(irarap) * (1.0d0-u1)
      a2 = y(ica40) * ratdum(ircagp) * u1
      nflowx = nflowx + 1
      string = 'ar36(ap,pg)ca40'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(ica40) * ratdum(ircaap) * (1.0d0-v1)   
      a2 = y(iti44) * ratdum(irtigp) * v1
      nflowx = nflowx + 1
      string = 'ca40(ap,pg)ti44'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(iti44) * ratdum(irtiap) * (1.0d0-w1)
      a2 = y(icr48) * ratdum(ircrgp) * w1
      nflowx = nflowx + 1
      string = 'ti44(ap,pg)cr48'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(icr48) * ratdum(ircrap) * (1.0d0-x1)   
      a2 = y(ife52) * ratdum(irfegp) * x1 
      nflowx = nflowx + 1
      string = 'cr48(ap,pg)fe52'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      a1 = y(ihe4)  * y(ife52) * ratdum(irfeap) * (1.0d0-y1)
      a2 = y(ini56) * ratdum(irnigp) * y1
      nflowx = nflowx + 1
      string = 'fe52(ap,pg)ni56'
      call sqeeze(string)
      flonam(nflowx) = string
      flofor(nflowx) = a1
      florev(nflowx) = a2
      flonet(nflowx) = flofor(nflowx) - florev(nflowx)

      return
      end






      subroutine init_aprox13
      include 'implno.dek'
      include 'network.dek'
c..
c..this routine initializes stuff for the aprox13 network
c..
c..declare
      integer          i


c..for easy zeroing of the isotope pointers
      integer          isotp(nisotp) 
      equivalence      (isotp(1),ih1)


c..zero all the isotope pointers
      do i=1,nisotp
       isotp(i)   = 0
      enddo

c..set the size of the network and the number of rates
      idnet   = idaprox13
      ionmax  = 13
      iener   = ionmax + 1
      itemp   = ionmax + 2
      iden    = ionmax + 3
      neqs    = iden
      nrat    = 59
      netname = 'aprox13'


c..set the id numbers of the elements

      ihe4  = 1
      ic12  = 2
      io16  = 3
      ine20 = 4
      img24 = 5
      isi28 = 6
      is32  = 7
      iar36 = 8
      ica40 = 9
      iti44 = 10
      icr48 = 11
      ife52 = 12
      ini56 = 13 

      ionbeg = 1
      ionend = 13


c..set the names of the elements
      ionam(ihe4)  = 'he4 '
      ionam(ic12)  = 'c12 '
      ionam(io16)  = 'o16 '
      ionam(ine20) = 'ne20'
      ionam(img24) = 'mg24'
      ionam(isi28) = 'si28'
      ionam(is32)  = 's32 '
      ionam(iar36) = 'ar36'
      ionam(ica40) = 'ca40'
      ionam(iti44) = 'ti44'
      ionam(icr48) = 'cr48'
      ionam(ife52) = 'fe52'
      ionam(ini56) = 'ni56'

c..set the number of nucleons in the element
      aion(ihe4)  = 4.0d0
      aion(ic12)  = 12.0d0   
      aion(io16)  = 16.0d0   
      aion(ine20) = 20.0d0   
      aion(img24) = 24.0d0  
      aion(isi28) = 28.0d0  
      aion(is32)  = 32.0d0  
      aion(iar36) = 36.0d0  
      aion(ica40) = 40.0d0  
      aion(iti44) = 44.0d0  
      aion(icr48) = 48.0d0  
      aion(ife52) = 52.0d0  
      aion(ini56) = 56.0d0  

c..set the number of protons in the element
      zion(ihe4)  = 2.0d0
      zion(ic12)  = 6.0d0
      zion(io16)  = 8.0d0
      zion(ine20) = 10.0d0   
      zion(img24) = 12.0d0  
      zion(isi28) = 14.0d0  
      zion(is32)  = 16.0d0  
      zion(iar36) = 18.0d0  
      zion(ica40) = 20.0d0  
      zion(iti44) = 22.0d0  
      zion(icr48) = 24.0d0  
      zion(ife52) = 26.0d0  
      zion(ini56) = 28.0d0  

c..set the number of neutrons
       do i=1,ionmax
        nion(i) = aion(i) - zion(i)
       enddo

c..set the binding energy of the element
      bion(ihe4)  =  28.29603d0 
      bion(ic12)  =  92.16294d0
      bion(io16)  = 127.62093d0 
      bion(ine20) = 160.64788d0 
      bion(img24) = 198.25790d0 
      bion(isi28) = 236.53790d0 
      bion(is32)  = 271.78250d0 
      bion(iar36) = 306.72020d0 
      bion(ica40) = 342.05680d0  
      bion(iti44) = 375.47720d0 
      bion(icr48) = 411.46900d0 
      bion(ife52) = 447.70800d0
      bion(ini56) = 484.00300d0 

c..set the partition functions - statistical weights, ground-state only here
      do i=1,ionmax
       wpart(i) = 1.0d0
      enddo

c..set the id numbers of the reaction rates
      ir3a   = 1
      irg3a  = 2
      ircag  = 3
      ir1212 = 4
      ir1216 = 5
      ir1616 = 6
      iroga  = 7
      iroag  = 8
      irnega = 9
      irneag = 10
      irmgga = 11
      irmgag = 12
      irsiga = 13
      irmgap = 14
      iralpa = 15
      iralpg = 16
      irsigp = 17
      irsiag = 18
      irsga  = 19
      irsiap = 20
      irppa  = 21
      irppg  = 22
      irsgp  = 23
      irsag  = 24
      irarga = 25
      irsap  = 26
      irclpa = 27
      irclpg = 28
      irargp = 29
      irarag = 30
      ircaga = 31
      irarap = 32
      irkpa  = 33
      irkpg  = 34
      ircagp = 35
      ircaag = 36
      irtiga = 37
      ircaap = 38
      irscpa = 39
      irscpg = 40
      irtigp = 41
      irtiag = 42
      ircrga = 43
      irtiap = 44
      irvpa  = 45
      irvpg  = 46
      ircrgp = 47
      ircrag = 48
      irfega = 49
      ircrap = 50
      irmnpa = 51
      irmnpg = 52
      irfegp = 53
      irfeag = 54
      irniga = 55
      irfeap = 56
      ircopa = 57
      ircopg = 58
      irnigp = 59

c..set the names of the reaction rates
      ratnam(ir3a)   = 'r3a  '
      ratnam(irg3a)  = 'rg3a '
      ratnam(ircag)  = 'rcag '
      ratnam(ir1212) = 'r1212'
      ratnam(ir1216) = 'r1216'
      ratnam(ir1616) = 'r1616'
      ratnam(iroga)  = 'roga '
      ratnam(iroag)  = 'roag '
      ratnam(irnega) = 'rnega'
      ratnam(irneag) = 'rneag'
      ratnam(irmgga) = 'rmgga'
      ratnam(irmgag) = 'rmgag'
      ratnam(irsiga) = 'rsiga'
      ratnam(irmgap) = 'rmgap'
      ratnam(iralpa) = 'ralpa'
      ratnam(iralpg) = 'ralpg'
      ratnam(irsigp) = 'rsigp'
      ratnam(irsiag) = 'rsiag'
      ratnam(irsga)  = 'rsga '
      ratnam(irsiap) = 'rsiap'
      ratnam(irppa)  = 'rppa '
      ratnam(irppg)  = 'rppg '
      ratnam(irsgp)  = 'rsgp '
      ratnam(irsag)  = 'rsag '
      ratnam(irarga) = 'rarga'
      ratnam(irsap)  = 'rsap '
      ratnam(irclpa) = 'rclpa'
      ratnam(irclpg) = 'rclpg'
      ratnam(irargp) = 'rargp'
      ratnam(irarag) = 'rarag'
      ratnam(ircaga) = 'rcaga'
      ratnam(irarap) = 'rarap'
      ratnam(irkpa)  = 'rkpa '
      ratnam(irkpg)  = 'rkpg '
      ratnam(ircagp) = 'rcagp'
      ratnam(ircaag) = 'rcaag'
      ratnam(irtiga) = 'rtiga'
      ratnam(ircaap) = 'rcaap'
      ratnam(irscpa) = 'rscpa'
      ratnam(irscpg) = 'rscpg'
      ratnam(irtigp) = 'rtigp'
      ratnam(irtiag) = 'rtiag'
      ratnam(ircrga) = 'rcrga'
      ratnam(irtiap) = 'rtiap'
      ratnam(irvpa)  = 'rvpa '
      ratnam(irvpg)  = 'rvpg '
      ratnam(ircrgp) = 'rcrgp'
      ratnam(ircrag) = 'rcrag'
      ratnam(irfega) = 'rfega'
      ratnam(ircrap) = 'rcrap'
      ratnam(irmnpa) = 'rmnpa'
      ratnam(irmnpg) = 'rmnpg'
      ratnam(irfegp) = 'rfegp'
      ratnam(irfeag) = 'rfeag'
      ratnam(irniga) = 'rniga'
      ratnam(irfeap) = 'rfeap'
      ratnam(ircopa) = 'rcopa'
      ratnam(ircopg) = 'rcopg'
      ratnam(irnigp) = 'rnigp'

      return
      end



cxxx


c..reaction rate library

c..torch rates
c..li7(t,n)   a(an,g)    be9(p,d)    be9(p,n)    b10(a,n)   b11(a,n)
c..n14(p,a)   c11(p,g)   c12(a,n)    c13(a,n)    c13(p,n)   c14(a,g)
c..c14(p,n)   c14(p,g)   o16(p,a)    n14(p,n)    n14(a,n)   n15(p,n)
c..n15(a,n)   n15(a,g)   o14(a,g)    o17(a,g)    o17(a,n)   o18(a,g)
c..o18(a,n)   ne20(p,a)  f18(p,g)    f19(p,g)    f19(p,n)   f19(a,p)   
c..na22(n,a)  ne20(p,g)  na23(p,a)   ne20(n,g)   ne21(p,g)  ne21(a,g)  
c..ne22(p,g)  ne22(a,g)  na22(n,p)   ne22(a,n)   na21(p,g)  mg24(p,a)
c..ne21(a,n)  na22(p,g)  na23(p,g)   na23(p,n)   mg24(p,g)  al27(p,a)
c..mg25(p,g)  mg25(a,p)  mg25(a,g)   mg25(a,n)   mg26(p,g)  mg26(a,g)
c..mg26(a,n)  al25(p,g)  al26(p,g)   al27(a,n)   si27(p,g)  si28(p,g)
c..si29(p,g)  si30(p,g)

c..bigbang rates:
c..n(e-nu)p   p(e-,nu)n  d(p,n)      d(n,g)      d(d,p)     d(d,n)      
c..t(p,n)     d(d,g)     t(p,g)      t(d,n)      t(t,2n)    he3(d,p)   
c..he3(t,d)   he3(t,np)  he4(np,g)   he4(d,g)    he4(t,n)   li6(p,he3) 
c..li6(n,g)   li7(d,n)   lit(t,2n)   li7(he3,np) li6(p,g)   li7(p,n)
c..be7(d,p)   be7(t,np)  be7(3he,2p) li6(a,g)    li7(a,n)   be9(p,g)   
c..b10(p,a)   li7(a,g)   b11(p,a)    be7(a,g)    b11(p,n)   b8(a,p)
c..b10(p,g)   c11(n,a)   be9(a,n)    b11(p,g)    b11(a,p)   

c..pp123 rates:
c..p(p,e+nu)  p(n,g)     d(p,g)      he3(n,g)    he3+he3    he3(a,g)    
c..be7(e-,nu) be7(p,g)   li7(p,g)    li7(p,a)    b8(e+,nu)

c..cno rates:
c..c12(p,g)   n13(e-nu)  c13(p,g)    n14(p,g)    o15(e-nu)  n14(a,g)
c..n15(p,g)   n15(p,a)   o16(p,g)    o17(p,a)    o17(p,g)   o18(p,a)   
c..o18(p,g)   f17(e-nu)  f18(e-nu)   f19(p,a)

c..hot cno rates
c..n13(p,g)   o14(e-nu)  o14(a,p)    o15(a,g)    f17(p,g)   ne18(e-nu)
c..f18(p,a)   ne18(a,p)  ne19(p,g)   si26(a,p)

c..alfa chain rates:
c..a(aa,g)    c12(a,g)   c12+c12     c12+o16     o16+o16    o16(a,g)    
c..ne20(a,g)  ne20(a,g)  mg24(a,g)   mg24(a,p)   al27(p,g)  si28(a,g)  
c..si28(a,p)  p31(p,g)   s32(a,g)    s32(a,p)    cl35(p,g)  ar36(a,g)
c..ar36(a,p)  k39(p,g)   ca40(a,g)   ca40(a,p)   sc43(p,g)  ti44(a,g)
c..ti44(a,p)  v47(p,g)   cr48(a,g)   cr(a,p)     mn51(p,g)  fe52(a,g)  
c..fe52(a,p)  co55(p,g)

c..photodisintegraation rates:
c..fe52(n,g) fe53(n,g)  fe54(p,g)






      subroutine tfactors(temp)
      include 'implno.dek'
      include 'tfactors.dek'

c..sets various populat temperature factors into common block
c..this routine must be called before any of the rates are called

c..declare the pass
      double precision temp

c..all these are in common block

      t9    = temp * 1.0d-9
      t92   = t9*t9
      t93   = t9*t92
      t94   = t9*t93
      t95   = t9*t94
      t96   = t9*t95

      t912  = sqrt(t9)
      t932  = t9*t912
      t952  = t9*t932
      t972  = t9*t952

      t913  = t9**oneth
      t923  = t913*t913
      t943  = t9*t913
      t953  = t9*t923
      t973  = t953*t923
      t9113 = t973*t943

      t914  = t9**(0.25d0)
      t934  = t914*t914*t914
      t954  = t9*t914
      t974  = t9*t934

      t915  = t9**onefif
      t935  = t915*t915*t915
      t945  = t915 * t935
      t965  = t9 * t915

      t917  = t9**onesev
      t927  = t917*t917
      t947  = t927*t927

      t918  = sqrt(t914)
      t938  = t918*t918*t918
      t958  = t938*t918*t918

      t9i   = 1.0d0/t9
      t9i2  = t9i*t9i
      t9i3  = t9i2*t9i

      t9i12 = 1.0d0/t912
      t9i32 = t9i*t9i12
      t9i52 = t9i*t9i32
      t9i72 = t9i*t9i52

      t9i13 = 1.0d0/t913
      t9i23 = t9i13*t9i13
      t9i43 = t9i*t9i13
      t9i53 = t9i*t9i23

      t9i14 = 1.0d0/t914
      t9i34 = t9i14*t9i14*t9i14
      t9i54 = t9i*t9i14

      t9i15 = 1.0d0/t915
      t9i35 = t9i15*t9i15*t9i15
      t9i45 = t9i15 * t9i35
      t9i65 = t9i*t9i15

      t9i17 = 1.0d0/t917
      t9i27 = t9i17*t9i17 
      t9i47 = t9i27*t9i27

      t9i18 = 1.0d0/t918
      t9i38 = t9i18*t9i18*t9i18
      t9i58 = t9i38*t9i18*t9i18

      return
      end





      subroutine rate_aan(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd,ddd

c..he4(an,g)be9
      aa  = 1.0d0 + 0.344*t9
      bb  = t92 * aa
      dbb = 2.0d0 * t9 * aa + t92*0.344

      cc  = 1.0d0/bb
      dcc = -cc*cc*dbb

      dd  = 2.59e-6 * exp(-1.062*t9i)
      ddd = dd*1.062*t9i2

      term    = cc * dd
      dtermdt = dcc*dd + cc*ddd

c..rates
      fr    = den * den * term 
      dfrdt = den * den * dtermdt * 1.0d-9
      dfrdd = 2.0d0 * den * term

      rev      = 5.84e19 * t93 * exp(-18.260*t9i)
      drevdt   = rev*(3.0d0*t9i + 18.260*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_be9pd(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd


c..be9(p,d)be8 =>2a
      aa  = 2.11e+11 * t9i23 * exp(-10.359*t9i13 - t92/0.2704)
      daa = aa*(-twoth*t9i + oneth*10.359*t9i43 - 2.0d0*t9/0.2704)

      bb  = 1.0d0  + 0.04*t913 + 1.09*t923 + 0.307*t9
     1      + 3.21*t943 + 2.30*t953
      dbb = oneth*0.04*t9i23 + twoth*1.09*t9i13 + 0.307
     1      + fourth*3.21*t913 + fiveth*2.30*t923 

      cc  = 5.79e+08 * t9i * exp(-3.046*t9i)
      dcc = cc*(-t9i + 3.046*t9i2)

      dd  = 8.50e+08 * t9i34 * exp(-5.800*t9i)
      ddd = dd*(-0.75d0*t9i + 5.800*t9i2)

      term    = aa*bb + cc + dd
      dtermdt = daa*bb + aa*dbb + dcc + ddd

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 8.07e-11 * t9i32 *exp(-7.555*t9i)
      drevdt   = rev*(-1.5d0*t9i + 7.555*t9i2)

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_be9pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,zz,dzz


c..be9(p,n)b9
      aa  = 5.58e7*(1.0d0 + 0.042*t912 + 0.985*t9)
      daa = 5.58e7*(0.5d0*0.042*t9i12 + 0.985)
    
      zz  = exp(-21.473*t9i)
      dzz = zz*21.473*t9i2

      bb  = aa * zz
      dbb = daa*zz + aa*dzz

      cc  = 1.02e+09 * t9i32 * exp(-26.725*t9i)
      dcc = cc*(-1.5d0*t9i + 26.725*t9i2)

      term    = bb + cc
      dtermdt = dbb + dcc

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term


      bb  = 0.998 * aa 
      dbb = 0.998 * daa  

      cc  = 0.998 * 1.02e+09 * t9i32 * exp(-5.252*t9i)
      dcc = cc*(-1.5d0*t9i + 5.252*t9i2)

      term    = bb + cc
      dtermdt = dbb + dcc

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_b10an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt


c..b10(a,n)n13
      term    = 1.20e+13 * t9i23 * exp(-27.989*t9i13 - t92/91.948921)
      dtermdt = -twoth*term*t9i 
     1          + term*(oneth*27.989*t9i43 - 2.0d0*t9/91.948921)

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.34 * exp(-12.287*t9i)
      drevdt   = rev*12.287*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev*term

      return
      end





      subroutine rate_b11an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff


c..b11(a,n)n14
      aa  = 6.97e+12 * t9i23 * exp(-28.234*t9i13 - t92/0.0196)
      daa = aa*(-twoth*t9i + oneth*28.234*t9i43 - 2.0d0*t9/0.0196)

      bb  = 1.0d0 + 0.015*t913 + 8.115*t923 + 0.838*t9  
     1      + 39.804*t943 + 10.456*t953
      dbb = oneth*0.015*t9i23 + twoth*8.115*t9i13 + 0.838
     1      + fourth*39.804*t913 + fiveth*10.456*t923 

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.79 * t9i32 * exp(-2.827*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.827*t9i2)

      ee   = 1.71e+03 * t9i32 * exp(-5.178*t9i)
      dee  = ee*(-1.5d0*t9i + 5.178*t9i2)

      ff   = 4.49e+06 * t935 * exp(-8.596*t9i)
      dff  = ff*(0.6d0*t9i + 8.596*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.67 * exp(-1.835*t9i)
      drevdt   = rev*1.835*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev*term

      return
      end





      subroutine rate_n14pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd,ddd,
     1                 t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz


c..n14(p,a)b11
       aa     = 1.0d0 + 0.0478*t9
       bb     = aa**twoth
       dbb    = twoth*bb/aa*0.0478

       zz     = 1.0d0/bb
       cc     = aa + 7.56e-03*t953*zz
       dcc    = 0.0478 + (fiveth*7.56e-3*t923 - 7.56e-3*t953*zz*dbb)*zz

       zz     = 1.0d0/cc
       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*dcc)*zz

       zz      = 1.0d0/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz*dt9a

       t9a56   = t9a**fivsix
       dt9a56  = fivsix * t9a56*zz*dt9a

       dd      = 2.63e+16 * t9a56 * t9i32 * exp(-31.883/t9a13)
       ddd     = dd*(dt9a56/t9a56 - 1.5d0*t9i
     1           + 31.883/t9a13**2 * dt9a13)

       term    = dd * exp(-33.915*t9i)
       dtermdt = term*(ddd/dd + 33.915*t9i2)

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 0.272 * dd
      drevdt   = 0.272 * ddd

      rr    = den * rev 
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev 

      return
      end





      subroutine rate_c11pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd


c..c11(p,g)n12
      aa  = 4.24e+04 * t9i23 * exp(-13.658*t9i13 - t92/2.647129)
      daa = aa*(-twoth*t9i + oneth*13.658*t9i43 - 2.0d0*t9/2.647129)

      bb  = 1.0d0  + 0.031*t913 + 3.11*t923 + 0.665*t9
     1      + 4.61*t943 + 2.50*t953
      dbb = oneth*0.031*t9i23 + twoth*3.11*t9i13 + 0.665
     1      + fourth*4.61*t913 + fiveth*2.50*t923 

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 8.84e+03 * t9i32 * exp(-7.021*t9i)
      ddd  = dd*(-1.5d0*t9i + 7.021*t9i2)

      term    = cc + dd 
      dtermdt = dcc + ddd 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.33e+10 * t932 * exp(-6.975*t9i)
      drevdt   = rev*(1.5d0*t9i + 6.975*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_c12an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,bb,dbb


c..c12(a,n)o15
      aa  = 2.48e7 * (1.0d0 + 0.188*t912 + 0.015*t9)
      daa = 2.48e7 * (0.5d0*0.188*t9i12 + 0.015)

      bb  = exp(-98.661*t9i)
      dbb = bb*98.661*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb
      

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = den * 1.41 * aa
      drrdt = den * 1.41 * daa  * 1.0d-9
      drrdd = 1.41 * aa

      return
      end





      subroutine rate_c13an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff,gg,dgg


c..c13(a,n)o16
      aa  = 6.77e+15 * t9i23 * exp(-32.329*t9i13 - t92/1.648656)
      daa = aa*(-twoth*t9i + oneth*32.329*t9i43 - 2.0d0*t9/1.648656)

      bb  = 1.0d0 + 0.013*t913 + 2.04*t923 + 0.184*t9
      dbb = oneth*0.013*t9i23 + twoth*2.04*t9i13 + 0.184

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 3.82e+05 * t9i32 * exp(-9.373*t9i)
      ddd  = dd*(-1.5d0*t9i + 9.373*t9i2)

      ee   = 1.41e+06 * t9i32 * exp(-11.873*t9i)
      dee  = ee*(-1.5d0*t9i + 11.873*t9i2)

      ff   = 2.0e+09 * t9i32 * exp(-20.409*t9i)
      dff  = ff*(-1.5d0*t9i + 20.409*t9i2)

      gg   = 2.92e+09 * t9i32 * exp(-29.283*t9i)
      dgg  = gg*(-1.5d0*t9i + 29.283*t9i2)

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 5.79e+00 * exp(-25.711*t9i)
      drevdt = rev*25.711*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_c13pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


c..c13(p,n)n13
      aa  = 1.88e+08*(1.0d0 - 0.167*t912 + 0.037*t9)
      daa = 1.88e+08*(0.037 - 0.5d0*0.167*t9i12)

      bb  = exp(-34.846*t9i)
      dbb = bb*34.846*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = den * 0.998 * aa
      drrdt = den * 0.998 * daa * 1.0d-9
      drrdd = 0.998 * aa

      return
      end






      subroutine rate_c14ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff


c..c14(a,g)o18
      aa  = 1.528e+09 * t9i23 * exp(-32.513*t9i13 - t92/7.086244)
      daa = aa*(-twoth*t9i + oneth*32.513*t9i43 - 2.0d0*t9/7.086244)

      bb  = 1.0d0 + 0.0128*t913 - 0.869*t923 - 0.0779*t9
     1      + 0.321*t943 + 0.0732*t953
      dbb = oneth*0.0128*t9i23 - twoth*0.869*t9i13 - 0.0779
     1      + fourth*0.321*t913 + fiveth*0.0732*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 3.375e+08 * t9i2 * exp(-32.513*t9i13)
      ddd  = dd*(-2.0d0*t9i + oneth*32.513*t9i43)

      ee   = 9.29e-08 * t9i32 * exp(-2.048*t9i)
      dee  = ee*(-1.5d0*t9i + 2.048*t9i2)

      ff   = 2.77e+03 * t9i45 * exp(-9.876*t9i)
      dff  = ff*(-0.8d0*t9i + 9.876*t9i2)

      term    = cc + dd + ee + ff 
      dtermdt = dcc + ddd + dee + dff 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.42e+10 * t932 * exp(-72.262*t9i)
      drevdt   = rev*(1.5d0*t9i + 72.262*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_c14pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,zz,dzz


c..c14(p,n)n14
      aa  = 7.19e+05*(1.0d0 + 0.361*t912 + 0.502*t9)
      daa = 7.19e+05*(0.5d0*0.361*t9i12 + 0.502)

      zz  = exp(-7.263*t9i)
      dzz = zz*7.263*t9i2

      bb  = aa * zz
      dbb = daa*zz + aa*dzz

      cc  = 3.34e+08 * t9i12 * exp(-12.246*t9i)
      dcc = cc*(-0.5d0*t9i + 12.246*t9i2)

      term    = bb + cc
      dtermdt = dbb + dcc  

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      cc  = 3.34e+08 * t9i12 * exp(-4.983*t9i)
      dcc = cc*(-0.5d0*t9i + 4.983*t9i2)

      rr    = den * 0.333 * (aa + cc)
      drrdt = den * 0.333 * (daa + dcc) * 1.0d-9
      drrdd = 0.333 * (aa + cc)

      return
      end






      subroutine rate_c14pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee


c..c14(p,g)n14
      aa  = 6.80e+06 * t9i23 * exp(-13.741*t9i13 - t92/32.729841)
      daa = aa*(-twoth*t9i + oneth*13.741*t9i43 - 2.0d0*t9/32.729841)

      bb  = 1.0d0 + 0.03*t913 + 0.503*t923 + 0.107*t9
     1      + 0.213*t943 + 0.115*t953
      dbb = oneth*0.03*t9i23 + twoth*0.503*t9i13 + 0.107
     1      + fourth*0.213*t913 + fiveth*0.115*t923 

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 5.36e+03 * t9i32 * exp(-3.811*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.811*t9i2)

      ee   = 9.82e+04 * t9i13 * exp(-4.739*t9i)
      dee  = ee*(-oneth*t9i + 4.739*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.00e+09 * t932 * exp(-118.452*t9i)
      drevdt   = rev*(1.5d0*t9i + 118.452*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_o16pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd,ddd,
     1                 t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz


c..o16(p,a)n13
       aa     = 1.0d0 + 0.0776*t9
       bb     = aa**twoth
       dbb    = twoth*bb/aa*0.0776

       zz     = 1.0d0/bb
       cc     = aa + 0.0264*t953*zz
       dcc    = 0.0776 + (fiveth*0.0264*t923 - 0.0264*t953*zz*dbb)*zz

       zz     = 1.0d0/cc
       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*dcc)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix * t9a56*zz

       dd      = 1.88e+18 * t9a56 * t9i32 * exp(-35.829/t9a13) 
       ddd     = dd*(dt9a56/t9a56 - 1.5d0*t9i
     1           + 35.829/t9a13**2 * dt9a13)

       term    = dd * exp(-60.561*t9i)
       dtermdt = term*(ddd/dd + 60.561*t9i2)

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 0.172 * dd
      drevdt   = 0.172 * ddd

      rr    = den * rev 
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev 

      return
      end





      subroutine rate_n14pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,bb,dbb


c..n14(p,n)o14
      aa  = 6.74e+07 * (1.0d0 + 0.658*t912 + 0.379*t9)
      daa = 6.74e+07 * (0.5d0*0.658*t9i12 + 0.379)

      bb  = exp(-68.762*t9i)
      dbb = bb*68.762*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = den * 2.99 * aa
      drrdt = den * 2.99 * daa * 1.0d-9
      drrdd = 2.99 * aa 

      return
      end





      subroutine rate_n14an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,zz,dzz


c..n14(a,n)f17
      aa  = 5.24e9*(1.0d0 - 1.15*t912 + 0.365*t9) 
      daa = 5.24e9*(0.365 - 0.5d0*1.15*t9i12) 

      zz  = exp(-t92/7.828804)
      dzz = -zz*2.0d0*t9/7.828804

      bb  = aa * zz
      dbb = daa*zz + aa*dzz

      cc   = 3.28e10 * t9i32 * exp(-1.5766e1*t9i)
      dcc  = cc*(-1.5d0*t9i + 1.5766e1*t9i2)

      term     = (bb + cc) * exp(-54.942*t9i)
      dtermdt  = term*((dbb+dcc)/(bb+cc) + 54.942*t9i2)

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term     = 1.48 * (bb + cc)
      dtermdt  = 1.48 * (dbb + dcc)

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_n15pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,t9a,aa,daa,bb,dbb


c..n15(p,n)o15
      t9a = min(t9,10.0d0)
      aa  = 3.51e+08 * (1.0d0 + 0.452*t912 - 0.191*t9a)
      if (t9a .eq. 10.0) then
       daa = 3.51e+08 * 0.5d0*0.452*t9i12 
      else 
       daa = 3.51e+08 * (0.5d0*0.452*t9i12 - 0.191)
      end if

      bb  = exp(-41.032*t9i)
      dbb = bb*41.032*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term    = 0.998 * aa
      dtermdt = 0.998 * daa 

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_n15an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,bb,dbb


c..n15(a,n)f18
      aa  = 3.14e8 * (1.0d0 - 0.641*t912 + 0.108*t9)
      daa = 3.14e8 * (0.108 - 0.5d0*0.641*t9i12)

      bb  = exp(-74.479*t9i)
      dbb = bb*74.479*t9i2

      term     = aa * bb
      dtermdt  = daa*bb + aa*dbb

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term     = 2.0d0 * aa
      dtermdt  = 2.0d0 * daa

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_n15ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee


c..n15(a,g)f19
      aa  = 2.54e+10 * t9i23 * exp(-36.211*t9i13 - t92/0.379456)
      daa = aa*(-twoth*t9i + oneth*36.211*t9i43 - 2.0d0*t9/0.379456)

      bb  = 1.0d0 + 0.012*t913 + 1.69*t923 + 0.136*t9
     1      + 1.91*t943 + 0.391*t953
      dbb = oneth*0.012*t9i23 + twoth*1.69*t9i13 + 0.136
     1      + fourth*1.91*t913 + fiveth*0.391*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 9.83e-03 * t9i32 * exp(-4.232*t9i)
      ddd  = dd*(-1.5d0*t9i + 4.232*t9i2)

      ee   = 1.52e+03 * t9 * exp(-9.747*t9i)
      dee  = ee*(t9i + 9.747*t9i2)

      term    = cc + dd + ee 
      dtermdt = dcc + ddd + dee 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 5.54e+10 * t932 * exp(-46.578*t9i)
      drevdt = rev*(1.5d0*t9i + 46.578*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_o14ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff


c..o14(a,g)ne18
      aa  = 9.47e+08 * t9i23 * exp(-39.388*t9i13 - t92/0.514089)
      daa = aa*(-twoth*t9i + oneth*39.388*t9i43 - 2.0d0*t9/0.514089)

      bb  = 1.0d0 + 0.011*t913 + 1.974*t923 + 0.146*t9
     1      + 3.036*t943 + 0.572*t953
      dbb = oneth*0.011*t9i23 + twoth*1.974*t9i13 + 0.146
     1      + fourth*3.036*t913 + fiveth*0.572*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.16e-01 * t9i32 * exp(-11.733*t9i)
      ddd  = dd*(-1.5d0*t9i + 11.733*t9i2)

      ee   = 3.39e+01 * t9i32 * exp(-22.609*t9i)
      dee  = ee*(-1.5d0*t9i + 22.609*t9i2)

      ff   = 9.10e-03 * t95 * exp(-12.159*t9i)
      dff  = ff*(5.0d0*t9i + 12.159*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 5.42e+10 * t932 * exp(-59.328*t9i)
      drevdt = rev*(1.5d0*t9i + 59.328*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_o17ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,
     2                 ft9a,dft9a,fpt9a,dfpt9a,gt9x,dgt9x,zz 


c..o17(a,g)ne21
       aa    = 1.0d0 + 0.1646*t9
       zz    = 1.0d0/aa 
       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.1646)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 0.786/t9a
       daa    = -aa*zz
       bb     = aa**3.51
       dbb    = 3.51*bb/aa * daa
       ft9a   = exp(-bb)
       dft9a  = -ft9a*dbb

       aa     = t9a/1.084
       bb     = aa**1.69
       dbb    = 1.69*bb/aa * dt9a/1.084 
       fpt9a  = exp(-bb)
       dfpt9a = -fpt9a*dbb

       aa     = oneth*exp(-10.106*t9i)
       daa    = aa*10.106*t9i2
       gt9x   = 1.0d0 + aa
       dgt9x  = daa

       zz     = 1.0d0/gt9x
       aa     = 1.73e17 * fpt9a*zz
       daa    = (1.73e17*dfpt9a - aa*dgt9x)*zz

       bb     = 3.50e15 * ft9a*zz
       dbb    = (3.50e15*dft9a - bb*dgt9x)*zz

       term    = (aa+bb) * t9a56 * t9i32 * exp(-39.914/t9a13)
       dtermdt = term*((daa+dbb)/(aa+bb) + dt9a56/t9a56 
     1           - 1.5d0*t9i + 39.914/t9a13**2 * dt9a13)

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 8.63e+10 * t932 * exp(-85.305*t9i)
      drevdt   = rev*(1.5d0*t9i + 85.305*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_o17an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,dd,
     1                 t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,gt9x,dgt9x,zz 


c..o17(a,n)ne20
       aa     = 1.0d0 + 0.0268*t9
       bb     = aa**twoth
       dbb    = twoth*bb/aa*0.0268

       zz     = 1.0d0/bb
       cc     = aa + 0.0232*t953*zz
       dcc    = 0.0268 + (fiveth*0.0232*t923 - 0.0232*t953*zz*dbb)*zz

       zz     = 1.0d0/cc
       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*dcc)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix * t9a56*zz

       dd     = oneth*exp(-10.106*t9i)
       gt9x   = 1.0d0 + dd
       dgt9x  = dd*10.106*t9i2

       term      = 1.03e+18/gt9x * t9a56 * t9i32 * exp(-39.914/t9a13)
       dtermdt   = term*(-dgt9x/gt9x + dt9a56/t9a56 
     1             - 1.5d0*t9i + 39.914/t9a13**2 * dt9a13)

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.86e+01 * exp(-6.852*t9i)
      drevdt   = rev*6.852*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_o18ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,theta
      parameter        (theta = 0.1d0)


c..o18(a,g)ne22
c..giessen et al 1994 nuc phys a 567, 146 for t9 less than 0.3
c..cf88 otherwise

      if (t9.lt.0.3) then
       aa   = 1.066d-41 * t9i32 * exp(-5.507d-01*t9i)
       daa  = aa*(-1.5d0*t9i + 5.507d-1*t9i2)

       bb   = 1.852d-13 * t9i32 * exp(-2.070*t9i)
       dbb  = bb*(-1.5d0*t9i + 2.070*t9i2)

       cc   = 1.431d-02 * t9i32 * exp(-4.462*t9i)
       dcc  = cc*(-1.5d0*t9i + 4.462*t9i2)

       dd   = 2.055d-04 * t9i32 * exp(-5.374*t9i)
       ddd  = dd*(-1.5d0*t9i + 5.374*t9i2)

       ee   = 5.332d+00 * t9i32 * exp(-6.285*t9i)
       dee  = ee*(-1.5d0*t9i + 6.285*t9i2)

       ff   = 1.457d+00 * t9i32 * exp(-7.121*t9i)
       dff  = ff*(-1.5d0*t9i + 7.121*t9i2)

       gg   = 3.121d-02 * t9i32 * exp(-7.292*t9i)
       dgg  = gg*(-1.5d0*t9i + 7.292*t9i2)

       hh   = 6.23d+03 * t9 * exp(-16.987*t9i)
       dhh  = hh*(t9i + 16.987*t9i2)

       term    = aa + bb + cc + dd + ee + ff + gg + hh
       dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg + dhh


      else 
       aa  = 1.82d+12 * t9i23 * exp(-40.057*t9i13 - t92/0.117649)
       daa = aa*(-twoth*t9i + oneth*40.057*t9i43 - 2.0d0*t9/0.117649)

       bb  = 1.0d0 + 0.01*t913 + 0.988*t923 + 0.072*t9 
     1       + 3.17*t943 + 0.586*t953
       dbb = oneth*0.01*t9i23 + twoth*0.988*t9i13 + 0.072
     1      + fourth*3.17*t913 + fiveth*0.586*t923

       cc   = aa * bb
       dcc  = daa*bb + aa*dbb

       dd   = 7.54 * t9i32 * exp(-6.228*t9i)
       ddd  = dd*(-1.5d0*t9i + 6.228*t9i2)

       ee   = 34.8 * t9i32 * exp(-7.301*t9i)
       dee  = ee*(-1.5d0*t9i + 7.301*t9i2)

       ff   = 6.23d+03 * t9 * exp(-16.987*t9i)
       dff  = ff*(t9i + 16.987*t9i2)

       gg   = theta * 1.0d-11 * t9i32 * exp(-1.994*t9i)
       dgg  = gg*(-1.5d0*t9i + 1.994*t9i2)

       term    = cc + dd + ee + ff + gg
       dtermdt = dcc + ddd + dee + dff + dgg
      end if

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.85d+10 * t932 * exp(-112.208*t9i)
      drevdt   = rev*(1.5d0*t9i + 112.208*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_o18an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,dd,
     1                 ee,dee,ff,dff,gg,dgg,hh,dhh,ft9a,dft9a,gt9,dgt9, 
     2                 t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,gt9i,zz


c..o18(a,n)ne21
      aa     = 1.0d0 + 0.0483*t9
      bb     = aa**twoth
      dbb    = twoth*bb/aa*0.0483

      zz     = 1.0d0/bb
      cc     = aa + 0.00569*t953*zz
      dcc    = 0.0483 + (fiveth*0.00569*t923 - 0.00569*t953*zz*dbb)*zz

      zz     = 1.0d0/cc
      t9a    = t9*zz
      dt9a   = (1.0d0 - t9a*dcc)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56   = t9a**fivsix
      dt9a56  = fivsix * t9a56*zz

      dd     = 5.0d0 * exp(-23.002*t9i)
      gt9    = 1.0d0 + dd
      gt9i   = 1.0d0/gt9
      dgt9   = dd*23.002*t9i2

      ee     = 0.431/t9a
      dee    = -ee*zz
      ff     = ee**3.89
      dff    = 3.89*ff/ee*dee
      ft9a   = exp(-ff)
      dft9a  = -ft9a*dff

      gg     = 7.22e+17 * ft9a*gt9i * t9a56 * t9i32 * exp(-40.056/t9a13)
      dgg    = gg*(-dff - gt9i*dgt9 + dt9a56/t9a56 - 1.5d0*t9i
     1         + 40.056/t9a13**2 *dt9a13) 

      hh     = 150.31 / gt9 * exp(-8.045*t9i)
      dhh    = hh*(-gt9i*dgt9 + 8.045*t9i2)

      term     = gg + hh
      dtermdt  = dgg + dhh

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

c..must protect the 8.045*t9i from overflow, so write it this way

      gg     = 7.22e+17*gt9i * t9a56 * t9i32 
     1         * exp(-ff - 40.056/t9a13 + 8.045*t9i)
      dgg    = gg*(-gt9i*dgt9 + dt9a56/t9a56 - 1.5d0*t9i
     1          - dff + 40.056/t9a13**2*dt9a13 - 8.045*t9i2)

      hh     = 150.31 * gt9i 
      dhh    = -hh*gt9i*dgt9 

      term    = 0.784 * (gg + hh)
      dtermdt = 0.784 * (dgg + dhh)

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_ne20pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd,ddd,
     1                 ee,dee,ff,dff,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,
     2                 zz,t9b


c..ne20(p,a)f17
       aa     = 1.0d0 + 0.0612*t9
       bb     = aa**twoth
       dbb    = twoth*bb/aa*0.0612

       zz     = 1.0d0/bb
       cc     = aa + 0.013*t953*zz
       dcc    = 0.0612 + (fiveth*0.013*t923 - 0.013*t953*zz*dbb)*zz

       zz     = 1.0d0/cc
       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*dcc)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix * t9a56*zz

       t9b     = min(t9,10.0d0)
       dd      = 5.31 + 0.544*t9b - 0.0523*t9b*t9b
       ddd     = 0.544 - 2.0d0*0.0523*t9b 
       if (t9b .eq. 10.0) ddd     = 0.0d0
       
       ee      = 3.25e19 * dd * t9a56 * t9i32 * exp(-43.176/t9a13)
       dee     = ee*(ddd/dd + dt9a56/t9a56 - 1.5d0*t9i
     1           +   43.176/t9a13**2 * dt9a13)

       ff      = exp(-47.969*t9i)
       dff     = ff*47.969*t9i2

       term    = ee * ff
       dtermdt = dee*ff + ee*dff


c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 0.0537 * ee
      drevdt   = 0.0537 * dee

      rr    = den * rev 
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev 

      return
      end





      subroutine rate_f18pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff,gg,dgg


c..f18(p,g)ne19
c..wiescher and kettner, apj 263, 891 1982
      aa  = 1.658e+7 * t9i23 * exp(-18.06*t9i13)
      daa = aa*(-twoth*t9i + oneth*18.06*t9i43)

      bb  = 4.604 + 0.106*t913 + 0.053*t923 + 0.009*t9
     1      - 0.036*t943 - 0.015*t953
      dbb = oneth*0.106*t9i23 + twoth*0.053*t9i13 + 0.009
     1      - fourth*0.036*t913 - fiveth*0.015*t923 

c..for temps greater than about t9 = 20, bb goes negative
      if (bb .le. 0.0) then
       bb = 0.0d0
       dbb = 0.0d0
      end if

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 4.55e-14 * t9i32* exp(-0.302*t9i) 
      ddd  = dd*(-1.5d0*t9i + 0.302*t9i2)

      ee   = 327.0 * t9i32 * exp(-3.84*t9i)
      dee  = ee*(-1.5d0*t9i + 3.84*t9i2)

      ff   = 1.32e+04 * t9i32 * exp(-5.22*t9i)
      dff  = ff*(-1.5d0*t9i + 5.22*t9i2)

      gg   = 93.0 * t9i32 * exp(-4.29*t9i)
      dgg  = gg*(-1.5d0*t9i + 4.29*t9i2)

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.73e+10 * t932 * exp(-74.396*t9i)
      drevdt   = rev*(1.5d0*t9i + 74.396*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_f19pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,zz


c..f19(p,g)ne20
      aa  = 6.04e+07 * t9i23 * exp(-18.113*t9i13 - t92/0.173056)
      daa = aa*(-twoth*t9i + oneth*18.113*t9i43 - 2.0d0*t9/0.173056)

      bb  = 1.0d0 + 0.023*t913 + 2.06*t923 + 0.332*t9  
     1      + 3.16*t943 + 1.30*t953
      dbb = oneth*0.023*t9i23 + twoth*2.06*t9i13 + 0.332
     1      + fourth*3.16*t913 + fiveth*1.30*t923 

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 6.32e+02 * t9i32 * exp(-3.752*t9i)  
      ddd  = dd*(-1.5d0*t9i + 3.752*t9i2)

      ee   = 7.56e+04 * t9i27 * exp(-5.722*t9i)
      dee  = ee*(-twosev*t9i + 5.722*t9i2)

      ff   = 7.0*exp(-16.44*t9i)
      dff  = ff*16.44*t9i2

      gg   = 4.0 * exp(-2.09*t9i) 
      dgg  = gg*2.09*t9i2

      hh   = 1.0d0 + ff + gg
      dhh  = dff + dgg

      zz      = 1.0d0/hh
      term    = (cc + dd + ee)*zz
      dtermdt = (dcc + ddd + dee - term*dhh)*zz

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.7e+10 * t932 * exp(-149.093*t9i)
      drevdt   = rev*(1.5d0*t9i + 149.093*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_f19pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,t9a,aa,daa,bb,dbb


c..f19(p,n)ne19
      aa  = 1.27e+08 * (1.0d0 - 0.147*t912 + 0.069*t9)
      daa = 1.27e+08 * (0.069 - 0.5d0*0.147*t9i12)

      bb  = exp(-46.659*t9i)
      dbb = bb*46.659*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term    = 0.998 * aa
      dtermdt = 0.998 * daa 

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_f19ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


c..f19(a,p)ne22
      aa  = 4.50e+18 * t9i23 * exp(-43.467*t9i13  - t92/0.405769)
      daa = -twoth*aa*t9i + aa*(oneth*43.467*t9i43 - 2.0d0*t9/0.405769)

      bb   = 7.98e+04 * t932 * exp(-12.760*t9i)
      dbb  = 1.5d0*bb*t9i + bb*12.760*t9i2

      term    = aa + bb
      dtermdt = daa + dbb 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.36 * exp(-19.439*t9i)
      drevdt   = rev*19.439*t9i2

      rr    = den * rev * term 
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_na22na(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 t9b,t9b2,t9b3


c..na22(n,a)f19
      t9b  = min(t9,10.0d0)
      t9b2 = t9b*t9b
      t9b3 = t9b2*t9b
      aa  = 1.0d0 + 0.8955*t9b - 0.05645*t9b2 + 7.302e-04*t9b3
      daa = 0.8955 - 2.0d0*0.05645*t9b + 3.0d0*7.302e-4*t9b2
      if (t9b .eq. 10.0) daa = 0.0d0

      term     = 1.21e6 * exp(aa)
      dtermdt  = term*daa

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.10 * exp(-22.620*t9i)
      drevdt   = rev*22.620*t9i2

      rr    = den * rev * term 
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_ne20pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff,gg,dgg,zz


c..ne20(p,g)na21
      aa  = 9.55e+06 * exp(-19.447*t9i13) 
      daa = aa*oneth*19.447*t9i43

      bb  = 1.0d0 + 0.0127*t9i23
      dbb = -twoth*0.0127*t9i53
   
      cc  = t92 * bb * bb
      dcc = 2.0d0*cc*t9i + 2.0d0*t92*bb*dbb

      zz  = 1.0d0/cc 
      dd  = aa*zz
      ddd = (daa - dd*dcc)*zz

      aa  = 2.05e+08 * t9i23 * exp(-19.447*t9i13)  
      daa = aa*(-twoth*t9i + oneth*19.447*t9i43)

      bb  = sqrt (t9/0.21)
      dbb = 0.5d0/(bb * 0.21)

      cc  = 2.67 * exp(-bb)
      dcc = -cc*dbb

      ff  = 1.0d0 + cc
      
      gg  = aa*ff
      dgg = daa*ff + aa*dcc 


      aa  = 18.0 * t9i32 * exp(-4.242*t9i)  
      daa = aa*(-1.5d0*t9i + 4.242*t9i2)

      bb  = 10.2 * t9i32 * exp(-4.607*t9i)
      dbb = bb*(-1.5d0*t9i + 4.607*t9i2)

      cc  = 3.6e+04 * t9i14 * exp(-11.249*t9i)
      dcc = cc*(-0.25d0*t9i + 11.249*t9i2)

      term    = dd + gg + aa + bb + cc
      dtermdt = ddd + dgg + daa + dbb + dcc


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.63e+09 * t932 * exp(-28.216*t9i)
      drevdt   = rev*(1.5d0*t9i + 28.216*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_na23pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,theta
      parameter        (theta = 0.1d0)


c..na23(p,a)ne20
c..el eid & champagne 1995 
      if (t9.le.2.0) then
       aa  = 1.26d+10 * t9i23 * exp(-20.758*t9i13 - t92/0.0169)
       daa = -twoth*aa*t9i + aa*(oneth*20.758*t9i43 - 2.0d0*t9/0.0169)

       bb  = 1.0d0  + 0.02*t913 - 13.8*t923 - 1.93*t9 
     1       + 234.0*t943 + 83.6*t953
       dbb = oneth*0.02*t9i23 - twoth*13.8*t9i13 - 1.93
     1       + fourth*234.0*t913 + fiveth*83.6*t923

       cc   = aa * bb
       dcc  = daa*bb + aa*dbb

       dd   = 4.38 * t9i32 * exp(-1.979*t9i)
       ddd  = -1.5d0*dd*t9i + dd*1.979*t9i2

       ee   = 6.50d+06 * (t9**(-1.366)) * exp(-6.490*t9i)
       dee  = -1.366d0*ee*t9i + ee*6.490*t9i2

       ff   = 1.19d+08 * (t9**1.055) * exp(-11.411*t9i)
       dff  = 1.055*ff*t9i + ff*11.411*t9i2

       gg   = theta * 9.91d-14 * t9i32 * exp(-0.418*t9i)
       dgg  = -1.5d0*gg*t9i + gg*0.418*t9i2      

       term    = cc + dd + ee + ff + gg 
       dtermdt = dcc + ddd + dee + dff + dgg 



c..cf88 + one term from gorres, wiesher & rolfs 1989, apj 343, 365
      else 
       aa  = 8.56d+09 * t9i23 * exp(-20.766*t9i13 - t92/0.017161)
       daa = -twoth*aa*t9i + aa*(oneth*20.766*t9i43 - 2.0d0*t9/0.017161)

       bb  = 1.0d0  + 0.02*t913 + 8.21*t923 + 1.15*t9 
     1       + 44.36*t943 + 15.84*t953
       dbb = oneth*0.02*t9i23 + twoth*8.21*t9i13 + 1.15
     1       + fourth*44.36*t913 + fiveth*15.84*t923

       cc   = aa * bb
       dcc  = daa*bb + aa*dbb

       dd   = 4.02 * t9i32 * exp(-1.99*t9i)
       ddd  = -1.5d0*dd*t9i + dd*1.99*t9i2

       ee   = 1.18d+04 * t9i54 * exp(-3.148*t9i)
       dee  = -1.25d0*ee*t9i + ee*3.148*t9i2

       ff   = 8.59d+05 * t943 * exp(-4.375*t9i)
       dff  = fourth*ff*t9i + ff*4.375*t9i2

       gg   = theta * 3.06d-12 * t9i32 * exp(-0.447*t9i)
       dgg  = -1.5d0*gg*t9i + gg*0.447*t9i2

       hh   = theta * 0.820 * t9i32 * exp(-1.601*t9i)
       dhh  = -1.5d0*hh*t9i + hh*1.601*t9i2

       term    = cc + dd + ee + ff + gg + hh
       dtermdt = dcc + ddd + dee + dff + dgg + dhh
      end if


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.25 * exp(-27.606*t9i)
      drevdt   = rev*27.606*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_ne20ng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc


c..ne20(n,g)ne21
c..wm88 Apj 239, 943; fit over range of experimental data, constant otherwise

      if (t9 .lt. 5.8025d-2) then
       term    = 5.449d+03
       dtermdt = 0.0d0

      else if (t9 .gt. 1.1605) then
       term    = 6.977d+04
       dtermdt = 0.0d0
  
      else if (t9 .ge. 5.8025d-2 .and. t9 .le. 2.9012d-1) then
       term    = 4.7219d+3 + 2.5248d+4*t9 - 2.7448d+5*t92
     1           + 9.2848d+5*t93
       dtermdt = 2.5248d+4 - 2.0d0*2.7448d+5*t9
     1           + 3.0d0*9.2848d+5*t92
       
      else

       aa  = 1.802d+04 * (t9/0.348)**4.43
       daa = 4.43 * aa * t9i

       bb  = -5.931 * (t9-0.348) + 1.268 * (t9-0.348)**2
       dbb = -5.931 + 2.0d0*1.268*(t9 - 0.348)

       cc  = exp(bb)
       dcc = cc*dbb

       term    = aa * cc
       dtermdt = daa*cc + aa*dcc

      end if

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      =  4.650d+09 * t932 * exp(-78.46*t9i)
      drevdt   = rev*(1.5d0*t9i + 78.46*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ne21pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,xx,dxx,zz,
     2                 theta
      parameter        (theta = 0.1d0)


c..ne21(p,g)na22

c..el eid & champagne 1995 

      if (t9.le.2.0) then
       aa  = 3.4d+08 * t9i23 * exp(-19.41*t9i13)
       daa = aa*(-twoth*t9i + oneth*19.41*t9i43)

       bb  = (16.7*t9 - 1.0)**2
       dbb = 2.0d0*(16.7*t9 - 1.0)*16.7

       cc  = 0.56 * exp(-bb)
       dcc = -cc*dbb

       dd  = 1.0d0 + cc
       ddd = dcc

       ee   = aa * dd
       dee  = daa*dd + aa*ddd

       ff   = 6.12 * t9i32 * exp(-1.403*t9i)
       dff  = ff*(-1.5d0*t9i + 1.403*t9i2)

       gg   = 1.35d+04 * t9i32 * exp(-3.008*t9i)
       dgg  = gg*(-1.5d0*t9i + 3.008*t9i2)      

       aa   = t9**0.67
       daa  = 0.67*aa*t9i
       zz   = 1.0d0/aa

       hh   = 3.12d+06 * t9**(-0.72) * exp(-8.268*zz)
       dhh  = hh*(-0.72d0*t9i + 8.268*zz*zz*daa)

       xx   = theta * 1.1d-03 * t9i32 * exp(-1.114*t9i)
       dxx  = xx*(-1.5d0*t9i + 1.114*t9i2)      

       term    = ee + ff + gg + hh + xx
       dtermdt = dee + dff + dgg + dhh + dxx


c..cf88 
      else 

       aa  = theta * 2.95d+08 * t9i23 * exp(-19.462*t9i13 -t92/0.003364)
       daa = aa*(-twoth*t9i + oneth*19.462*t9i43 - 2.0d0*t9/0.003364)

       bb  = 1.0d0 + 0.021*t913 + 13.29*t923 + 1.99*t9 
     1       + 124.1*t943 + 47.29*t953
       dbb = oneth*0.021*t9i23 + twoth*13.29*t9i13 + 1.99
     1       + fourth*124.1*t913 + fiveth*47.29*t923

       cc   = aa * bb
       dcc  = daa*bb + aa*dbb

       dd   = theta * 7.80d-01 * t9i32 * exp(-1.085*t9i)
       ddd  = dd*(-1.5d0*t9i + 1.085*t9i2)

       ee   = 4.37d+08 * t9i23 * exp(-19.462*t9i13)
       dee  = ee*(-twoth*t9i + oneth*19.462*t9i43)

       ff   = 5.85 * t9i32 * exp(-1.399*t9i)
       dff  = ff*(-1.5d0*t9i + 1.399*t9i2)

       gg   = 1.29d+04 * t9i32 * exp(-3.009*t9i)
       dgg  = gg*(-1.5d0*t9i + 3.009*t9i2)

       hh   = 3.15d+05 * t9i35 * exp(-5.763*t9i)
       dhh  = hh*(-0.6d0*t9i + 5.763*t9i2)

       term    = cc + dd + ee + ff + gg + hh
       dtermdt = dcc + ddd + dee + dff + dgg + dhh
      end if


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.06d+10 * t932 * exp(-78.194*t9i)
      drevdt   = rev*(1.5d0*t9i + 78.194*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_ne21ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,gg,dgg,hh,dhh,zz,
     2                 t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56


c..ne21(a,g)mg25
       aa    = 1.0d0 + 0.0537*t9
       zz    = 1.0d0/aa

       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.0537)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 8.72e-03*t9 - 6.87e-04*t92 + 2.15e-05*t93
       daa    = 8.72e-3 - 2.0d0*6.87e-4*t9 + 3.0d0*2.15e-5*t92

       bb     = 1.52e-04 * exp(-46.90*t9i13*aa)
       dbb    = bb*46.90*(oneth*t9i43*aa - t9i13*daa) 

       cc     = 1.5*exp(-4.068*t9i)
       dcc    =  cc*4.068*t9i2

       gg     = 2.0 * exp(-20.258*t9i)
       dgg    = gg*20.258*t9i2

       hh     = 1.0d0 + cc + gg
       dhh    = dcc + dgg

       zz     = 1.0d0/hh
       dd     = bb*zz
       ddd    = (dbb - dd*dhh)*zz

       aa     = 4.94e+19 * t9a56 * t9i32 * exp(-46.89/t9a13)
       daa    = aa*(dt9a56/t9a56 - 1.5d0*t9i  
     1              + 46.89/t9a13**2 * dt9a13)

       bb     =  2.66e+07 * t9i32 * exp(-22.049*t9i)
       dbb    = bb*(-1.5d0*t9i + 22.049*t9i2)

       cc     = aa + bb
       dcc    = daa + dbb

       term    = dd * cc
       dtermdt = ddd*cc + dd*dcc 


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 4.06e+10 * t932 * exp(-114.676*t9i)
      drevdt = rev*(1.5d0*t9i + 114.676*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_ne21an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,zz,
     2                 t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56


c..ne21(a,n)mg24
       aa    = 1.0d0 + 0.0537*t9
       zz    = 1.0d0/aa

       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.0537)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 4.94e+19 * t9a56 * t9i32 * exp(-46.89/t9a13)
       daa    = aa*(dt9a56/t9a56 - 1.5d0*t9i  
     1              + 46.89/t9a13**2 * dt9a13)

       bb     =  2.66e+07 * t9i32 * exp(-22.049*t9i)
       dbb    = bb*(-1.5d0*t9i + 22.049*t9i2)

       cc     = 2.0d0*exp(-20.258*t9i)
       dcc    = cc*20.258*t9i2
   
       dd     = 1.5*exp(-4.068*t9i)
       ddd    = dd*4.068*t9i2

       ee     = 1.0d0 + cc + dd
       dee    = dcc + ddd

       zz      = 1.0d0/ee
       term    = (aa + bb)*zz
       dtermdt = (daa + dbb - term*dee)*zz


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 12.9 * exp(-29.606*t9i)
      drevdt   = rev*29.606*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_ne22pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff,gg,dgg,theta
      parameter        (theta = 0.1d0)


c..ne22(p,g)na23

c..el eid & champagne 1995 

      if (t9.le.2.0) then
       aa  = 1.05d+09 * t9i23 * exp(-19.431*t9i13)
       daa = aa*(-twoth*t9i + oneth*19.431*t9i43)

       bb  = 1.24d-09 * t9i32 * exp(-0.414*t9i)
       dbb = bb*(-1.5d0*t9i + 0.414*t9i2)

       cc  = 2.90d-02 * t9i32 * exp(-1.752*t9i)
       dcc = cc*(-1.5d0*t9i + 1.752*t9i2)

       dd  = 9.30d+04 * t9**(-1.174) * exp(-5.100*t9i)
       ddd = dd*(-1.174*t9i + 5.100*t9i2)

       ee   = 5.71d+05 * t9**(0.249) * exp(-7.117*t9i)
       dee  = ee*(0.249*t9i + 7.117*t9i2)

       ff   = theta * 3.25d-04 * t9i32 * exp(-0.789*t9i)
       dff  = ff*(-1.5d0*t9i + 0.789*t9i2)

       gg   = theta * 0.10 * t9i32 * exp(-1.161*t9i)
       dgg  = gg*(-1.5d0*t9i + 1.161*t9i2)      

       term    = aa + bb + cc + dd + ee + ff + gg
       dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg


c..cf88 
      else 

       aa  = 1.15d+09 * t9i23 * exp(-19.475*t9i13)
       daa = aa*(-twoth*t9i + oneth*19.475*t9i43)

       bb  = 9.77d-12 * t9i32 * exp(-0.348*t9i)
       dbb = bb*(-1.5d0*t9i + 0.348*t9i2)

       cc   = 8.96d+03 * t9i32 * exp(-4.84*t9i)
       dcc  = cc*(-1.5d0*t9i + 4.84*t9i2)

       dd   = 6.52d+04 * t9i32 * exp(-5.319*t9i)
       ddd  = dd*(-1.5d0*t9i + 5.319*t9i2)

       ee   = 7.97d+05 * t9i12 * exp(-7.418*t9i)
       dee  = ee*(-0.5d0*t9i + 7.418*t9i2)

       ff   = theta * 1.63d-01 * t9i32 * exp(-1.775*t9i)
       dff  = ff*(-1.5d0*t9i + 1.775*t9i2)

       term    = aa + bb + cc + dd + ee + ff 
       dtermdt = daa + dbb + dcc + ddd + dee + dff 

      end if


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.67d+09 * t932 * exp(-102.048*t9i)
      drevdt   = rev*(1.5d0*t9i + 102.048*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ne22ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,res1,dres1,
     2                 ft9a,dft9a,fpt9a,dfpt9a,gt9x,dgt9x,
     3                 t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,
     4                 rdmass,res2,zz
      parameter        (rdmass = 22.0d0*4.0d0/26.0d0,
     1                  res2   = -11.604d0 * 22.0d0/26.0d0)


c..ne22(a,g)mg26
c..kappeler 1994 apj 437, 396 

      if (t9 .lt. 1.25) then

       res1 = 1.54d-01*(t9*rdmass)**(-1.5)
       dres1 = -1.5d0 * res1 * t9i

       aa    = 1.7d-36 * res1 * exp(res2*t9i*0.097)
       daa   = aa/res1*dres1 - aa*res2*0.097*t9i2

       bb    = 1.5d-7 * res1 * exp(res2*t9i*0.400)
       dbb   = bb/res1*dres1 - bb * res2 * 0.400 * t9i2

       cc    = 0.5 * res1 * 3.7d-2 * exp(res2*t9i*0.633)  
       dcc   = cc/res1*dres1 - cc*res2*0.633*t9i2

       dd    = res1 * 3.6d+1 * exp(res2*t9i*0.828) 
       ddd   = dd/res1*dres1 - dd*res2*0.828*t9i2

       term    = aa + bb + cc + dd
       dtermdt = daa + dbb + dcc + ddd


c..cf88
      else
       aa    = 1.0d0 + 0.0548*t9
       zz    = 1.0d0/aa
    
       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.0548)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 0.197/t9a
       daa    = -aa*zz
       bb     = aa**4.82
       dbb    = 4.82*bb/aa * daa
       ft9a   = exp(-bb)
       dft9a  = -ft9a*dbb

       aa     = t9a/0.249
       bb     = aa**2.31
       dbb    = 2.31*bb/aa * dt9a/0.249
       fpt9a  = exp(-bb)
       dfpt9a = -fpt9a*dbb

       aa     = 5.0d0*exp(-14.791*t9i)
       daa    = aa*14.791*t9i2
       gt9x   = 1.0d0 + aa
       dgt9x  = daa

       zz     = 1.0d0/gt9x
       aa     = 4.16e19 * fpt9a*zz
       daa    = (4.16e19*dfpt9a - aa*dgt9x)*zz

       bb     = 2.08e16 * ft9a*zz
       dbb    = (2.08e16*dft9a - bb*dgt9x)*zz

       term    = (aa+bb) * t9a56 * t9i32 * exp(-47.004/t9a13)
       dtermdt = term*((daa+dbb)/(aa+bb) + dt9a56/t9a56
     1                 - 1.5d0*t9i + 47.004/t9a13**2 * dt9a13)
      end if


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 6.15d+10 * t932 * exp(-123.151*t9i)
      drevdt = rev*(1.5d0*t9i + 123.151*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_na22np(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa


c..na22(n,p)ne22
      aa  = 1.0d0 - 3.037e-02*t9 + 8.380e-03*t92 - 7.101e-04*t93
      daa =  -3.037e-02 + 2.0d0*8.380e-03*t9 - 3.0d0*7.101e-04*t92

      term    = 1.24e+08 * exp(aa)
      dtermdt = term*daa

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev     = 7.01*exp(-42.059*t9i)
      drevdt  = rev*42.059*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_ne22an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,gg,dgg,ft9a,dft9a,gt9x,dgt9x,
     2                 t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,res1,res2,
     3                 zz 
      parameter        (res1 = 2.4731857150793075d-2, 
     1                  res2 = -9.8187694549560547d0)

c..note: res1=1.54d-1*(88./26.)**(-1.5)   res2=-11.604*(22./26.)

c..ne22(a,n)mg25
c..kappeler 1994 apj 437, 396 ; wiescher suggest only 828 kev, ignore 633 kev

      if (t9 .lt. 0.6) then
       term = res1*1.64d+02 * t9i32 * exp(t9i*0.828*res2)
       dtermdt = -1.5d0*term*t9i - term*res2*0.828*t9i2 

c..cf88
      else 
       aa    = 1.0d0 + 0.0548*t9
       zz    = 1.0d0/aa

       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.0548)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 0.197/t9a
       daa    = -aa*zz
       bb     = aa**4.82
       dbb    = 4.82*bb/aa * daa
       ft9a   = exp(-bb)
       dft9a  = -ft9a*dbb

       gg     = bb
       dgg    = dbb

       aa     = 5.0d0*exp(-14.791*t9i)
       daa    = aa*14.791*t9i2
       gt9x   = 1.0d0 + aa
       dgt9x  = daa

       zz     = 1.0d0/gt9x
       aa     = ft9a*zz
       daa    = (dft9a - aa*dgt9x)*zz

       bb     = 4.16e+19 * t9a56 * t9i32 * exp(-47.004/t9a13)
       dbb    = bb*(dt9a56/t9a56 - 1.5d0*t9i
     1              + 47.004/t9a13**2 * dt9a13)

       cc     = aa*bb
       dcc    = daa*bb + aa*dbb

       dd      = 1.44e-04*zz * exp(-5.577*t9i) 
       ddd     = -dd*zz*dgt9x + dd*5.577*t9i2

       term    = cc + dd
       dtermdt = dcc + ddd
      end if


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 7.833d-5
      drevdt = 0.0d0
      if (t9 .gt. 0.008) then
       rev    = 0.544 * exp(5.577*t9i)
       drevdt = -rev*5.577*t9i2
      end if

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_na21pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd


c..na21(p,g)mg22
      aa  = 1.41e+05 * t9i23 * exp(-20.739*t9i13 -  t92/0.133956)
      daa = aa*(-twoth*t9i + oneth*20.739*t9i43 - 2.0d0*t9/0.133956)

      bb  = 1.0d0 + 0.020*t913 + 4.741*t923 + 0.667*t9
     1      + 16.380*t943 + 5.858*t953
      dbb = oneth*0.020*t9i23 + twoth*4.741*t9i13 + 0.667
     1      + fourth*16.380*t913 + fiveth*5.858*t923 

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 6.72e+02 * t9i34 * exp(-2.436*t9i)
      ddd  = dd*(-0.75d0*t9i + 2.436*t9i2)

      term    = cc + dd 
      dtermdt = dcc + ddd 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.44e+10 * t932 * exp(-63.790*t9i)
      drevdt   = rev*(1.5d0*t9i + 63.790*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg24pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,gg,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz


c..mg24(p,a)na21
       aa     = 1.0d0 + 0.127*t9
       zz     = 1.0d0/aa

       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*0.127)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix * t9a56*zz

       gg      = min(t9,12.0d0)
       aa      = 4.43 + 3.31*gg - 0.229*gg*gg
       daa     = 3.31 - 2.0d0*0.229*gg
       if (gg .eq. 12.0) daa = 0.0d0

       bb      = 1.81e21 * t9a56 * t9i32 * exp(-49.967/t9a13)
       dbb     = bb*(dt9a56/t9a56 - 1.5d0*t9i
     1               + 49.967/t9a13**2 * dt9a13)

       cc      = aa*bb
       dcc     = daa*bb + aa*dbb

       dd      = exp(-79.843*t9i)
       ddd     = dd*79.843*t9i2

       term    = cc * dd
       dtermdt = dcc*dd + cc*ddd

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 0.0771 * cc
      drevdt   = 0.0771 * dcc

      rr    = den * rev 
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev 

      return
      end






      subroutine rate_na22pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


c..na22(p,g)mg23
      aa  = 9.63e-05 * t932 * exp(-0.517*t9i)
      daa = aa*(1.5d0*t9i + 0.517*t9i2)

      bb  = 2.51e+04 * t9 * exp(-2.013*t9i)
      dbb = bb*(t9i + 2.013*t9i2)

      term    = aa + bb 
      dtermdt = daa + dbb 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.27e+10 * t932 * exp(-87.933*t9i)
      drevdt   = rev*(1.5d0*t9i + 87.933*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_na23pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff,gg,dgg,hh,hhi,xx,dxx,theta
      parameter        (theta = 0.1d0)


c..na23(p,g)mg24

c..el eid & champagne 1995 
      if (t9.le.2.0) then
       aa  = 2.47d+09 * t9i23 * exp(-20.758*t9i13)
       daa = aa*(-twoth*t9i + oneth*20.758*t9i43)

       bb  = 9.19d+01 * t9i32 * exp(-2.789*t9i)
       dbb = bb*(-1.5d0*t9i + 2.789*t9i2)

       cc  = 1.72d+04 * t9i32 * exp(-3.433*t9i)
       dcc = cc*(-1.5d0*t9i + 3.433*t9i2)

       dd  = 3.44d+04 * t9**0.323 * exp(-5.219*t9i)
       ddd = dd*(0.323*t9i + 5.219*t9i2)

       ee   = theta * 2.34d-04 * t9i32 * exp(-1.590*t9i)
       dee  = ee*(-1.5d0*t9i + 1.590*t9i2)

       term    = aa + bb + cc + dd + ee 
       dtermdt = daa + dbb + dcc + ddd + dee 


c..cf88 + gorres, wiesher & rolfs 1989, apj 343, 365
      else 
 
       aa  = 2.93d+08 * t9i23 * exp(-20.766*t9i13 - t92/0.088209)
       daa = aa*(-twoth*t9i + oneth*20.766*t9i43 - 2.0d0*t9/0.088209)

       bb  = 1.0d0 + 0.02*t913 + 1.61*t923 + 0.226*t9 
     1       + 4.94*t943 + 1.76*t953
       dbb = oneth*0.02*t9i23 + twoth*1.61*t9i13 + 0.226
     1      + fourth*4.94*t913 + fiveth*1.76*t923 

       xx  = aa * bb
       dxx = daa*bb + aa*dbb

       cc   = 9.34d+01 * t9i32 * exp(-2.789*t9i)
       dcc  = cc*(-1.5d0*t9i + 2.789*t9i2)

       dd   = 1.89d+04 * t9i32 * exp(-3.434*t9i)
       ddd  = dd*(-1.5d0*t9i + 3.434*t9i2)

       ee   = 5.1d+04 * t915 * exp(-5.51*t9i)
       dee  = ee*(0.2d0*t9i + 5.51*t9i2)

       ff   = theta * 0.820 * t9i32 * exp(-1.601*t9i)
       dff  = ff*(-1.5d0*t9i + 1.601*t9i2)

       gg   = 1.5 * exp(-5.105*t9i)
       dgg  = gg*5.105*t9i2

       hh   = 1.0d0 + gg
       hhi  = 1.0d0/hh

       term    = (xx + cc + dd + ee + ff) * hhi
       dtermdt = (dxx + dcc + ddd + dee + dff - term*dgg)*hhi

      end if


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.49d+10 * t932 * exp(-135.665*t9i)
      drevdt   = rev*(1.5d0*t9i + 135.665*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_na23pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 t9a,dt9a,t9a32,dt9a32,zz


c..na23(p,n)mg24
      aa  = 1.0d0 + 0.141*t9
      zz  = 1.0d0/aa

      t9a = t9*zz
      dt9a = (1.0d0 - t9a*0.141)*zz

      aa    = sqrt(t9a)
      t9a32 = t9a * aa
      dt9a32 = 1.5d0 * aa * dt9a

      bb   = 9.29d8 * (1.0d0 - 0.881d0 * t9a32 * t9i32)
      dbb  = -9.29d8 * 0.881d0 * t9i32*(dt9a32 - 1.5d0*t9a32*t9i)

      cc   = exp(-56.173*t9i)
      dcc  = cc*56.173*t9i2

      term    = bb * cc
      dtermdt = dbb*cc + bb*dcc

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term    = 0.998 * bb
      dtermdt = 0.998 * dbb

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end






      subroutine rate_mg24pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff,gg,ggi


c..mg24(p,g)al25
      aa  = 5.60e+08 * t9i23 * exp(-22.019*t9i13) 
      daa = aa*(-twoth*t9i + oneth*22.019*t9i43)

      bb  = 1.0d0 + 0.019*t913 - 0.173*t923 - 0.023*t9
      dbb = oneth*0.019*t9i23 - twoth*0.173*t9i13 - 0.023

c..stop negative rates above t9 = 10
      if (bb .le. 0.0) then
       bb  = 0.0d0
       dbb = 0.0d0
      end if 

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd   = 1.48e+03 * t9i32 * exp(-2.484*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.484*t9i2)

      ee   = 4.00e+03 * exp(-4.180*t9i)
      dee  = ee*4.180*t9i2

      ff   = 5.0 * exp(-15.882*t9i)
      dff  = ff*15.882*t9i2

      gg   = 1.0d0 + ff
      ggi  = 1.0d0/gg

      term    = (cc + dd + ee) * ggi
      dtermdt = (dcc + ddd + dee - term*dff)*ggi


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.13e+09 * t932 * exp(-26.358*t9i)
      drevdt   = rev*(1.5d0*t9i + 26.358*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_al27pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff,gg,dgg,theta
      parameter        (theta = 0.1d0) 


c..al27(p,a)mg24
c..champagne 1996
      aa  = 4.71d+05 * t9i23 * exp(-23.25*t9i13 - 3.57*t92)
      daa = -twoth*aa*t9i + aa*(oneth*23.25*t9i43 - 2.0d0*3.57*t9)

      bb  = 1.0d0 + 0.018*t913 - 7.29*t923 - 0.914*t9 
     1      + 77.2*t943 + 24.6*t953
      dbb = oneth*0.018*t9i23 - twoth*7.29*t9i13 - 0.914
     1      + fourth*77.2*t913 + fiveth*24.6*t923 

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 2.23d+04 * (t9**3.989) * exp(-2.148 * t9**(-1.293))
      ddd = 3.989*dd*t9i + 1.293*dd*2.148*t9**(-2.293)

      ee   = 0.17 * 1.29d-09 * t9i32 * exp(-0.836*t9i)
      dee  = -1.5d0*ee*t9i + ee*0.836*t9i2

      ff   = theta * 2.73d-03 * t9i32 * exp(-2.269*t9i)
      dff  = -1.5d0*ff*t9i + ff*2.269*t9i2

      gg   = theta * 2.60d-02 * t9i32 * exp(-2.492*t9i)
      dgg  = -1.5d0*gg*t9i + gg*2.492*t9i2

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 1.81*exp(-18.572*t9i)
      drevdt   = rev*18.572*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_mg25pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff,gg,dgg


c..mg25(p,g)al26
      aa  = 3.57e+09 * t9i23 * exp(-22.031*t9i13 - t92/0.0036)
      daa = aa*(-twoth*t9i + oneth*22.031*t9i43 - 2.0d0*t9/0.0036)

      bb  = 1.0d0 + 0.019*t913 + 7.669*t923 + 1.015*t9 
     1      + 167.4*t943 + 56.35*t953
      dbb = oneth*0.019*t9i23 + twoth*7.669*t9i13 + 1.015
     1      + fourth*167.4*t913 + fiveth*56.35*t923 

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 3.07e-13 * t9i32 * exp(-0.435*t9i)  
      ddd = dd*(-1.5d0*t9i + 0.435*t9i2)

      ee   = 1.94e-07 * t9i32 * exp(-0.673*t9i)
      dee  = ee*(-1.5d0*t9i + 0.673*t9i2)

      ff   = 3.15e-05 * t9**(-3.40)* exp(-1.342*t9i - t92/169.0)
      dff  = ff*(-3.40d0*t9i + 1.342*t9i2 - 2.0d0*t9/169.0)

      gg   = 1.77e+04 * t958 * exp(-3.049*t9i - t92/169.0)
      dgg  = gg*(0.625*t9i + 3.049*t9i2 - 2.0d0*t9/169.0)

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 1.03e+10 * t932 * exp(-73.183*t9i)
      drevdt   = rev*(1.5d0*t9i + 73.183*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_mg25ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc


c..mg25(a,p)si28
      aa  = -23.271*t9i13 + 6.46*t9 - 2.39*t92 + 0.506*t93 
     1      - 6.04e-2*t94 + 3.75e-3*t95 - 9.38e-5*t96 
 
      daa = oneth*23.271*t9i43 + 6.46 - 2.0d0*2.39*t9 + 3.0d0*0.506*t92 
     1      - 4.0d0*6.04e-2*t93 + 5.0d0*3.75e-3*t94 - 6.0d0*9.38e-5*t95 

      bb  = 3.23e8 * t9i23 * exp(aa)
      dbb  = -twoth*bb*t9i + bb*daa

c..dbb/bb
      cc   = -twoth*t9i + daa

      term    = bb * exp(-13.995*t9i)
      dtermdt = term*cc + term*13.995*t9i2

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 2.86 * bb
      drevdt   = 2.86 * dbb

      rr    = den * rev 
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end





      subroutine rate_mg25ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,gt9x,dgt9x,t9a,dt9a,t9a13,dt9a13,
     2                 t9a56,dt9a56,zz


c..mg25(a,g)si29
      aa    = 1.0d0 + 0.0630*t9
      zz    = 1.0d0/aa

      t9a   = t9*zz
      dt9a  = (1.0d0 - t9a*0.0630)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix * t9a56*zz

      aa     = oneth*10.0d0*exp(-13.180*t9i)
      daa    = aa*13.180*t9i2
      gt9x   = 1.0d0 + aa
      dgt9x  = daa

      bb     = 1.0d0/gt9x
      dbb    = -bb*bb*dgt9x

      cc     = 3.59e+20 * bb * t9a56 * t9i32 * exp(-53.41/t9a13)
      dcc    = cc*(dbb*gt9x + dt9a56/t9a56 - 1.5d0*t9i
     1             + 53.41/t9a13**2 * dt9a13) 

      dd     = 0.0156*t9 - 1.79e-03*t92 + 9.08e-05*t93
      ddd    = 0.0156 - 2.0d0*1.79e-03*t9 + 3.0d0*9.08e-05*t92

      ee     = 5.87e-04*exp(-53.42*t9i13*dd)
      dee    = ee*53.42*(oneth*t9i43*dd - t9i13*ddd)

      term    = cc * ee
      dtermdt = dcc*ee + cc*dee

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 1.90e+11 * t932 * exp(-129.128*t9i)
      drevdt = rev*(1.5d0*t9i + 129.128*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg25an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,zz,
     1                 gt9x,dgt9x,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56


c..mg25(a,n)si28
      aa    = 1.0d0 + 0.0630*t9
      zz    = 1.0d0/aa

      t9a   = t9*zz
      dt9a  = (1.0d0 - t9a*0.0630)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix * t9a56*zz

      aa     = oneth*10.0d0*exp(-13.180*t9i)
      daa    = aa*13.180*t9i2
      gt9x   = 1.0d0 + aa
      dgt9x  = daa

      bb     = 1.0d0/gt9x
      dbb    = -bb*bb*dgt9x

      term    = 3.59e+20 * bb * t9a56 * t9i32 * exp(-53.41/t9a13)
      dtermdt = term*(dbb*gt9x + dt9a56/t9a56 - 1.5d0*t9i
     1          + 53.41/t9a13**2 * dt9a13) 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 20.0*exp(-30.792*t9i)
      drevdt   = rev*30.792*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_mg26pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,theta
      parameter        (theta = 0.1d0)


c..mg26(p,g)al27
c..champagne 1996 
      aa  = 8.54d-12 * t9i32 * exp(-0.605*t9i)
      daa = aa*(-1.5d0*t9i + 0.605*t9i2)

      bb  = 2.75d-06 * t9i32 * exp(-1.219*t9i)
      dbb = bb*(-1.5d0*t9i + 1.219*t9i2)

      cc  = 1.30d-02 * t9i32 * exp(-1.728*t9i)
      dcc = cc*(-1.5d0*t9i + 1.728*t9i2)

      dd  = 8.06d+00 * t9i32 * exp(-2.537*t9i) 
      ddd = dd*(-1.5d0*t9i + 2.537*t9i2)

      ee   = 1.45d+03 * t9i32 * exp(-3.266*t9i) 
      dee  = ee*(-1.5d0*t9i + 3.266*t9i2)

      ff   = 4.03d+04 * t9i32 * exp(-3.784*t9i)
      dff  = ff*(-1.5d0*t9i + 3.784*t9i2)

      gg   = 8.82d+04 * t9**(-0.21) * exp(-4.194*t9i)
      dgg  = gg*(-0.21*t9i + 4.194*t9i2)

      hh   = theta * 1.93d-05 * t9i32 * exp(-1.044*t9i)
      dhh  = hh*(-1.5d0*t9i + 1.044*t9i2)

      term    = aa + bb + cc + dd + ee + ff + gg + hh
      dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg + dhh


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.14d+09 * t932 * exp(-95.99*t9i)
      drevdt   = rev*(1.5d0*t9i + 95.99*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg26ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,gt9x,dgt9x,t9a,dt9a,t9a13,dt9a13,
     2                 t9a56,dt9a56,zz


c..mg26(a,g)si30
      aa    = 1.0d0 + 0.0628*t9
      zz    = 1.0d0/aa

      t9a   = t9*zz
      dt9a  = (1.0d0 - t9a*0.0628)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix * t9a56*zz

      aa     = 5.0d0*exp(-20.990*t9i)
      daa    = aa*20.990*t9i2
      gt9x   = 1.0d0 + aa
      dgt9x  = daa

      bb     = 1.0d0/gt9x
      dbb    = -bb*bb*dgt9x

      cc     = 2.93e+20 * bb * t9a56 * t9i32 * exp(-53.505/t9a13)
      dcc    = cc*(dbb*gt9x + dt9a56/t9a56 - 1.5d0*t9i
     1             + 53.505/t9a13**2 * dt9a13) 

      dd     = 0.0751*t9 - 0.0105*t92 + 5.57e-04*t93
      ddd    = 0.0751 - 2.0d0*0.0105*t9 + 3.0d0*5.57e-04*t92

      ee     = 4.55e-2 * exp(-53.51*t9i13*dd)
      dee    = ee*53.51*(oneth*t9i43*dd - t9i13*ddd)

      term    = cc * ee
      dtermdt = dcc*ee + cc*dee

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 6.38e+10 * t932 * exp(-123.52*t9i)
      drevdt = rev*(1.5d0*t9i + 123.52*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg26an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,zz,
     1                 gt9x,dgt9x,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56


c..mg26(a,n)si29
      aa    = 1.0d0 + 0.0628*t9
      zz    = 1.0d0/aa

      t9a   = t9*zz
      dt9a  = (1.0d0 - t9a*0.0628)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix * t9a56*zz

      aa     = 5.0d0*exp(-20.990*t9i)
      daa    = aa*20.990*t9i2
      gt9x   = 1.0d0 + aa
      dgt9x  = daa

      bb     = 1.0d0/gt9x
      dbb    = -bb*bb*dgt9x

      term   = 2.93e+20 * bb * t9a56 * t9i32 * exp(-53.505/t9a13)
      dtermdt= term*(dbb*gt9x + dt9a56/t9a56 - 1.5d0*t9i
     1         + 53.505/t9a13**2 * dt9a13) 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.68*exp(-0.401*t9i)
      drevdt   = rev*0.401*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_al25pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee

c..al25(p,g)si26
c..coc et al 1995 a&a 299, 479 , case b

      aa  = 8.98d+1 * t9i32 * exp(-4.874*t9i)
      daa = aa*(-1.5d0*t9i + 4.874*t9i2)

      bb  = 1.568d+3 * t9i32 * exp(-9.632*t9i)
      dbb = bb*(-1.5d0*t9i + 9.632*t9i2)

      cc  = 2.42d+8 * t9i23 * exp(-23.18*t9i13)
      dcc = cc*(-twoth*t9i + oneth*23.18*t9i43)

      dd  = 4.10d-02 * t9i32 * exp(-1.741*t9i)
      ddd = dd*(-1.5d0*t9i + 1.741*t9i2)

      ee  = 2.193d+3 * t9i32 * exp(-4.642*t9i)
      dee = ee*(-1.5d0*t9i + 4.642*t9i2)

      term    = aa + bb + cc + dd + ee 
      dtermdt = daa + dbb + dcc + ddd + dee 


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.117d+11 * t932 * exp(-64.048*t9i)
      drevdt   = rev*(1.5d0*t9i + 64.048*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_al26pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff,theta
      parameter        (theta = 0.1d0)


c..al26(p,g)si27
c..coc et al 1995 a&a 299, 479 

      aa  = 1.53d+9 * t9**(-1.75) * exp(-23.19*t9i13)
      daa = aa*(-1.75*t9i + oneth*23.19*t9i43)

      bb  = theta*8.7d-7 * t9i32 * exp(-0.7845*t9i)
      dbb = bb*(-1.5d0*t9i + 0.7845*t9i2)

      cc  = theta*1.00d-3 * t9i32 * exp(-1.075*t9i)
      dcc = cc*(-1.5d0*t9i + 1.075*t9i2)

      dd  = 9.00d+00 * t9i32 * exp(-2.186*t9i)
      ddd = dd*(-1.5d0*t9i + 2.186*t9i2)

      ee  = 5.05d+02 * t9i32 * exp(-3.209*t9i)
      dee = ee*(-1.5d0*t9i + 3.209*t9i2)

      ff  = 9.45d+03 * t9i * exp(-4.008*t9i)
      dff = ff*(-t9i + 4.008*t9i2)

      term    = aa + bb + cc + dd + ee + ff
      dtermdt = daa + dbb + dcc + ddd + dee + dff 


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.46d+10 * t932 * exp(-86.621*t9i)
      drevdt   = rev*(1.5d0*t9i + 86.621*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_al27an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


c..al27(a,n)p30

      aa    = 8.2e+04*exp(-30.588*t9i)
      daa   = aa*30.588*t9i2

      bb    = 5.21e+05 * t974 * exp(-33.554*t9i)
      dbb   = 1.75d0*bb*t9i + bb*33.554*t9i2

      term    = aa + bb
      dtermdt = daa + dbb

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term


      aa  = 5.21e+05 * t974 * exp(-2.966*t9i)
      daa = aa*(1.75d0*t9i + 2.966*t9i2) 

      rev      = 6.75d0 * (8.20e4 + aa)
      drevdt   = 6.75d0 * daa

      rr    = den * rev 
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev 

      return
      end





      subroutine rate_si27pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd


c..si27(p,g)p28

      aa  = 1.64e+09 * t9i23 * exp(-24.439*t9i13)
      daa = aa*(-twoth*t9i + oneth*24.439*t9i43)

      bb  = 2.00e-08 * t9i32 * exp(-0.928*t9i) 
      dbb = bb*(-1.5d0*t9i + 0.928*t9i2)

      cc  = 1.95e-02 * t9i32 * exp(-1.857*t9i)
      dcc = cc*(-1.5d0*t9i + 1.857*t9i2)

      dd  = 3.70e+02 * t9i47 * exp(-3.817*t9i)
      ddd = dd*(-foursev*t9i + 3.817*t9i2)

      term    = aa + bb + cc + dd 
      dtermdt = daa + dbb + dcc + ddd 


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.62e+10 * t932 * exp(-23.960*t9i)
      drevdt   = rev*(1.5d0*t9i + 23.960*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si28pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,xx,dxx



c..si28(p,g)p29

c..champagne et al 96  

      if (t9.le.5.0) then

       aa  = 8.44d+08 * t9i23 * exp(-24.389*t9i13 - t92/8.4681)
       daa = aa*(-twoth*t9i + oneth*24.389*t9i43 - 2.0d0*t9/8.4681)

       bb  = 1.0d0 + 0.17*t913 + 0.113*t923 + 0.0135*t9
     1       + 0.194*t943 + 0.0591*t953
       dbb = oneth*0.17*t9i23 + twoth*0.113*t9i13 + 0.0135
     1      + fourth*0.194*t913 + fiveth*0.0591*t923 

       xx  = aa * bb
       dxx = daa*bb + aa*dbb

       cc   = 2.92d+02 * t9i32 * exp(-4.157*t9i)
       dcc  = cc*(-1.5d0*t9i + 4.157*t9i2)

       dd   = 4.30d+05 * t9i32 * exp(-18.51*t9i)
       ddd  = dd*(-1.5d0*t9i + 18.51*t9i2)

       ee   = 6.05d+03 * t9i32 * exp(-18.17*t9i)
       dee  = ee*(-1.5d0*t9i + 18.17*t9i2)

       term    = xx + cc + dd + ee 
       dtermdt = dxx + dcc + ddd + dee 


c..cf88 
      else 

       aa  = 1.64d+08 * t9i23 * exp(-24.449*t9i13 - t92/8.4681)
       daa = aa*(-twoth*t9i + oneth*24.449*t9i43 - 2.0d0*t9/8.4681)

       bb  = 1.0d0 + 0.017*t913 - 4.11*t923 - 0.491*t9
     1       + 5.22*t943 + 1.58*t953
       dbb = oneth*0.017*t9i23 - twoth*4.11*t9i13 - 0.491
     1      + fourth*5.22*t913 + fiveth*1.58*t923 

       xx  = aa * bb
       dxx = daa*bb + aa*dbb

       cc   = 3.52d+02 * t9i32 * exp(-4.152*t9i)
       dcc  = cc*(-1.5d0*t9i + 4.152*t9i2)

       dd   = 6.3d+05 * t9i32 * exp(-18.505*t9i)
       ddd  = dd*(-1.5d0*t9i + 18.505*t9i2)

       ee   = 1.69d+03 * exp(-14.518*t9i)
       dee  = ee*14.518*t9i2

       term    = xx + cc + dd + ee 
       dtermdt = dxx + dcc + ddd + dee 

      end if


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.46d+09 * t932 * exp(-31.879*t9i)
      drevdt   = rev*(1.5d0*t9i + 31.879*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si29pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,xx,dxx



c..si29(p,g)p30

      aa  = 3.26e+09 * t9i23 * exp(-24.459*t9i13 - t92/0.065536)
      daa = aa*(-twoth*t9i + oneth*24.459*t9i43 - 2.0d0*t9/0.065536)

      bb  = 1.0d0 + 0.017*t913 + 4.27*t923 + 0.509*t9
     1      + 15.40*t943 + 4.67*t953
      dbb = oneth*0.017*t9i23 + twoth*4.27*t9i13 + 0.509
     1     + fourth*15.40*t913 + fiveth*4.67*t923 

      xx  = aa * bb
      dxx = daa*bb + aa*dbb

      cc   = 2.98e+03 * t9i32 * exp(-3.667*t9i) 
      dcc  = cc*(-1.5d0*t9i + 3.667*t9i2)

      dd   = 3.94e+04 * t9i32 * exp(-4.665*t9i)
      ddd  = dd*(-1.5d0*t9i + 4.665*t9i2)

      ee   = 2.08e+04 * t912 * exp(-8.657*t9i)
      dee  = ee*(0.5d0*t9i + 8.657*t9i2)

      term    = xx + cc + dd + ee 
      dtermdt = dxx + dcc + ddd + dee 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.26e+10 * t932 * exp(-65.002*t9i)
      drevdt   = rev*(1.5d0*t9i + 65.002*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si30pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,xx,dxx



c..si30(p,g)p31

      aa  = 4.25e8 * t9i23 * exp(-24.468*t9i13 - t92/0.4489)  
      daa = aa*(-twoth*t9i + oneth*24.468*t9i43 - 2.0d0*t9/0.4489)  

      bb  = 1.0d0 + 0.017*t913 + 0.150*t923 + 0.018*t9 
     1      + 5.53*t943 + 1.68*t953  
      dbb = oneth*0.017*t9i23 + twoth*0.150*t9i13 + 0.018
     1     + fourth*5.53*t913 + fiveth*1.68*t923 

      xx  = aa * bb
      dxx = daa*bb + aa*dbb

      cc   = 1.86e4 * t9i32 * exp(-5.601*t9i)
      dcc  = cc*(-1.5d0*t9i + 5.601*t9i2)

      dd   = 3.15e5 * t9i32 * exp(-6.961*t9i)        
      ddd  = dd*(-1.5d0*t9i + 6.961*t9i2)

      ee   = 2.75e5 * t9i12 * exp(-10.062*t9i)        
      dee  = ee*(-0.5d0*t9i + 10.062*t9i2)

      term    = xx + cc + dd + ee 
      dtermdt = dxx + dcc + ddd + dee 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.50e9 * t932 * exp(-84.673*t9i)
      drevdt   = rev*(1.5d0*t9i + 84.673*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end








      subroutine rate_weaknp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,
     1                 zm1,zm2,zm3,zm4,zm5,
     2                 c1,c2
      parameter        (c1 = 1.0d0/5.93d0,
     1                  c2 = 0.98d0/886.7d0)


c..n(e-nu)p and p(e-,nu)n
c..from schramm and wagoner 1977
c..currently accepted best value for the neutron lifetime, 
c..886.7 (+/- 1.9) seconds. P.R. Huffman et al., Nature, 6 January 2000.

      zm1   = t9 * c1
      zm2   = zm1*zm1 
      zm3   = zm1*zm2 
      zm4   = zm1*zm3 
      zm5   = zm1*zm4 

      aa   = 27.512*zm5 + 36.492*zm4 + 11.108*zm3 
     1       - 6.382*zm2 + 0.565*zm1 + 1.0d0
      daa  = (5.0d0*27.512*zm4 + 4.0d0*36.492*zm3 + 3.0d0*11.108*zm2
     1       - 2.0d0*6.382*zm1 + 0.565)*c1

c..n=>p 
      fr    = c2 * aa
      dfrdt = c2 * daa * 1.0d-9
      dfrdd = 0.0d0


      aa  = 27.617*zm5 + 34.181*zm4 + 18.059*zm3 
     1      - 16.229*zm2 + 5.252*zm1
      daa = (5.0d0*27.617*zm4 + 4.0d0*34.181*zm3 + 3.0d0*18.059*zm2
     1      - 2.0d0*16.229*zm1 + 5.252)*c1

      bb = exp(-2.531d0/zm1)
      dbb = bb*2.531d0/zm2*c1

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

c..p=>n
      rr    = c2 * cc
      drrdt = c2 * dcc * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_dpn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc


c..d(p,n)2p
       aa   = 3.35e7 * exp(-3.720*t9i13)
       daa  = aa*oneth*3.720d0*t9i43 

       bb   = 1.0d0 + 0.784*t913 + 0.346*t923 + 0.690*t9 
       dbb  = oneth*0.784*t9i23 + twoth*0.346*t9i13 + 0.690

       term    = aa * bb
       dtermdt = daa * bb + aa * dbb

c..rate 
      cc = exp(-25.815*t9i)
      dcc = cc*25.815*t9i2

      fr    = den * cc * term 
      dfrdt = den * (dcc*term + cc*dtermdt) * 1.0d-9
      dfrdd = cc * term 

      rev      =  4.24e-10 * t9i32
      drevdt   = -1.5d0*rev*t9i

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_dng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,c1
      parameter        (c1 = 66.2d0*18.9d0)


c..d(n,g)t
      term    = 66.2 * (1.0d0  + 18.9*t9)
      dtermdt = c1

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      =  1.63e+10 * t9i32 * exp(-72.62*t9i)
      drevdt   = rev*(-1.5d0*t9i + 72.62*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_ddp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb



c..d(d,p)t
      aa  = 4.13e8 * t9i23 * exp(-4.258*t9i13)
      daa = -twoth*aa*t9i + oneth*aa*4.258*t9i43

      bb  = 1.0d0 + 0.098*t913 + 4.39e-2*t923 + 3.01e-2*t9
     1      + 0.543*t943 + 0.946*t953
      dbb = oneth*0.098*t9i23 + twoth*4.39e-2*t9i13 + 3.01e-2
     1      + fourth*0.543*t913 + fiveth*0.946*t923 

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 1.73 * exp(-46.798*t9i) 
      drevdt   = rev*46.798*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_ddn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb



c..d(d,n)he3
      aa  = 3.88e8 * t9i23 * exp(-4.258*t9i13)  
      daa = -twoth*aa*t9i + oneth*aa*4.258*t9i43

      bb  = 1.0d0 + 0.098*t913 + 0.418*t923 + 0.287*t9  
     1      + 0.638*t943 + 1.112*t953 
      dbb = oneth*0.098*t9i23 + twoth*0.418*t9i13 + 0.287
     1      + fourth*0.638*t913 + fiveth*1.112*t923 

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 1.730 * exp(-37.935*t9i) 
      drevdt   = rev*37.935*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_tpn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa


c..t(p,n)he3
      term    = 7.07e8 * (1.0d0 - 0.15*t912 + 0.098*t9)
      dtermdt = 7.07e8 * (-0.5d0*0.15*t9i12 + 0.098)

      aa  = exp(-8.863*t9i)
      daa = aa*8.863*t9i2

c..rate 
      fr    = den * aa * term 
      dfrdt = den * (daa*term + aa*dtermdt) * 1.0d-9
      dfrdd = aa * term 

      rev      = 0.998
      drevdt   = 0.0d0

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_ddg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


c..d(d,g)he4
      aa  = 4.84e+01 * t9i23 * exp(-4.258*t9i13)  
      daa = aa*(-twoth*t9i + oneth*4.258*t9i43)

      bb  = 1.0d0 + 0.098*t913 - 0.203*t923 - 0.139*t9  
     1      + 0.106*t943 + 0.185*t953 
      dbb = oneth*0.098*t9i23 - twoth*0.203*t9i13 - 0.139
     1      + fourth*0.106*t913 + fiveth*0.185*t923 

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 4.53e+10 * t932 * exp(-276.729*t9i)
      drevdt   = rev*(1.5d0*t9i + 276.729*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_tpg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


c..t(p,g)he4
      aa  = 2.20e+04 * t9i23 * exp(-3.869*t9i13)
      daa = aa*(-twoth*t9i + oneth*3.869*t9i43)

      bb  = 1. + 0.108*t913 + 1.68*t923 + 1.26*t9
     1      + 0.551*t943 + 1.06*t953 
      dbb = oneth*0.108*t9i23 + twoth*1.68*t9i13 + 1.26
     1      + fourth*0.551*t913 + fiveth*1.06*t923 

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 2.61e+10 * t932 * exp(-229.932*t9i)
      drevdt   = rev*(1.5d0*t9i + 229.932*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_tdn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc


c..t(d,n)he4 ; the "dt" reaction
      aa  = 8.09e+10 * t9i23 * exp(-4.524*t9i13 - t92/0.0144)
      daa = -twoth*aa*t9i + aa*(oneth*4.524*t9i43 - 2.0d0*t9/0.0144)

      bb  = 1.0d0 + 0.092*t913 + 1.80*t923 + 1.16*t9
     1      + 10.52*t943 + 17.24*t953
      dbb = oneth*0.092*t9i23 + twoth*1.80*t9i13 + 1.16
     1      + fourth*10.52*t913 + fiveth*17.24*t923 

      cc  = 8.73e+08 * t9i23 * exp(-0.523*t9i)
      dcc = -twoth*cc*t9i + cc*0.523*t9i2
      
      term    = aa * bb + cc
      dtermdt = daa*bb + aa*dbb + dcc

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 5.54*exp(-204.117*t9i)   
      drevdt   = rev*204.117*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_tt2n(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


c..t(t,2n)he4 
      aa  = 1.67e+09 * t9i23 * exp(-4.872*t9i13)  
      daa = aa*(-twoth*t9i + oneth*4.872*t9i43)

      bb  = 1.0d0 + 0.086*t913 - 0.455*t923 - 0.272*t9  
     1      + 0.148*t943 + 0.225*t953
      dbb = oneth*0.086*t9i23 - twoth*0.455*t9i13 - 0.272
     1      + fourth*0.148*t913 + fiveth*0.225*t923 

      term    = aa * bb 
      dtermdt = daa*bb + aa*dbb 

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 3.38e-10 * t9i32 * exp(-131.504*t9i)
      drevdt   = rev*(-1.5d0*t9i + 131.504*t9i2)

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_he3dp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc


c..he3(d,p)he4 
      aa  = 5.86e+10 * t9i23 * exp(-7.181*t9i13 - t92/0.099225)
      daa = -twoth*aa*t9i + aa*(oneth*7.181*t9i43 - 2.0d0*t9/0.099225)

      bb  = 1.0d0 + 0.058*t913 + 0.142*t923 + 0.0578*t9
     1      + 2.25*t943 + 2.32*t953
      dbb = oneth*0.058*t9i23 + twoth*0.142*t9i13 + 0.0578
     1      + fourth*2.25*t913 + fiveth*2.32*t923 

      cc  = 4.36e+08 * t9i12 * exp(-1.72*t9i)
      dcc = -0.5d0*cc*t9i + cc*1.72*t9i2
      
      term    = aa * bb + cc
      dtermdt = daa*bb + aa*dbb + dcc

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 5.55*exp(-212.980*t9i)   
      drevdt   = rev*212.980*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_he3td(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,t9a,dt9a,
     1                 t9a13,dt9a13,t9a56,dt9a56,zz


c..he3(t,d)he4 
      aa       = 1.0d0 + 0.128*t9
      zz       = 1.0d0/aa

      t9a      = t9*zz
      dt9a     = (1.0d0 - t9a*0.128)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56    = t9a**fivsix
      dt9a56   = fivsix*t9a56*zz

      term     = 5.46e+09 * t9a56 * t9i32 * exp(-7.733/t9a13)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i
     1                + 7.733/t9a13**2 * dt9a13)

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 1.60*exp(-166.182*t9i)
      drevdt   = rev*166.182*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_he3tnp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,t9a,dt9a,
     1                 t9a13,dt9a13,t9a56,dt9a56,zz


c..he3(t,np)he4 
      aa       = 1.0d0 + 0.115*t9
      zz       = 1.0d0/aa

      t9a      = t9*zz
      dt9a     = (1.0d0 - t9a*0.115)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56    = t9a**fivsix
      dt9a56   = fivsix*t9a56*zz

      term     = 7.71e+09 * t9a56 * t9i32 * exp(-7.733/t9a13)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i
     1                + 7.733/t9a13**2 * dt9a13)

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 3.39e-10*t9i32 * exp(-140.367*t9i)
      drevdt   = rev*(-1.5d0*t9i + 140.367*t9i2)

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_he4npg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


c..he4(np,g)li6
      aa  = 4.62e-6 * t9i2 * exp(-19.353*t9i)
      daa = aa*(-2.0d0*t9i + 19.353*t9i2)

      bb  = 1.0d0 + 0.075*t9
      dbb = 0.075

      term    = aa * bb 
      dtermdt = daa*bb + aa*dbb 

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 7.22e19 * t93 * exp(-42.933*t9i)
      drevdt   = rev*(3.0d0*t9i + 42.933*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_he4dg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc


c..he4(d,g)li6
      aa  = 3.01e1 * t9i23 * exp(-7.423*t9i13) 
      daa = aa*(-twoth*t9i + oneth*7.423*t9i43)

      bb  = 1.0d0 + 0.056*t913 - 4.85*t923 + 8.85*t9 
     1      - 0.585*t943 - 0.584*t953 
      dbb = oneth*0.056*t9i23 - twoth*4.85*t9i13 + 8.850
     1      - fourth*0.585*t913 - fiveth*0.584*t923

c..rate goes negative for t9 greater than about 15, so try this 
      if (bb .le. 0.0) then
       bb = 0.0d0
       dbb = 0.0d0
      end if

      cc =  8.55e1 * t9i32 * exp(-8.228*t9i)
      dcc = cc*(-1.5d0*t9i + 8.228*t9i2)

      term    = aa * bb + cc
      dtermdt = daa*bb + aa*dbb + dcc


c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 1.53e10 * t932 * exp(-17.1180*t9i)  
      drevdt   = rev*(1.5d0*t9i + 17.1180*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_he4tn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,t9a,dt9a,t9a32,dt9a32,zz


c..he4(t,n)li6
      aa     = 1.0d0 + 49.180*t9   
      zz     = 1.0d0/aa
     
      t9a   = t9*zz
      dt9a   = (1.0d0 - t9a*49.180)*zz

      t9a32  = t9a * sqrt(t9a) 
      dt9a32 = 1.5d0*t9a32/t9a * dt9a 

      aa     = 1.80e8 * exp(-55.4940*t9i) 
      daa    = aa*55.4940*t9i2

      bb     = 1.0d0 - 0.2610 * t9a32 * t9i32  
      dbb    = -0.2610*(-1.5d0*t9a32*t9i52 + dt9a32*t9i32) 

      cc     = aa*bb
      dcc    = daa*bb + aa*dbb 

      dd     = 2.72e9 * t9i32 * exp(-57.8840*t9i)   
      ddd    = dd*(-1.5d0*t9i + 57.8840*t9i2)

      term    = cc + dd
      dtermdt = dcc + ddd

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 2.72e9 * t9i32 * exp(-2.39*t9i) 
      drevdt   = rev*(-1.5d0*t9i + 2.39*t9i2)

      term = 0.935*(1.80e8*bb + rev)
      dtermdt = 0.935*(1.80e8*dbb + drevdt)

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_li6phe3(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd 


c..li6(p,he3)he4 
      aa  = 3.73e10 * t9i23 * exp(-8.413*t9i13 - t92/30.25) 
      daa = -twoth*aa*t9i + aa*(oneth*8.413*t9i43 - 2.0d0*t9/30.25)

      bb  = 1.0d0 + 0.050*t913 - 0.061*t923 - 0.0210*t9  
     1      + 0.0060*t943 + 0.0050*t953  
      dbb = oneth*0.050*t9i23 - twoth*0.061*t9i13 - 0.0210
     1      + fourth*0.0060*t913 + fiveth*0.0050*t923 

      cc  = 1.33e10 * t9i32 * exp(-17.7630*t9i) 
      dcc = -1.5d0*cc*t9i + cc*17.7630*t9i2

      dd  =  1.29e9 * t9i * exp(-21.82*t9i)
      ddd = -dd*t9i + dd*21.82*t9i2
      
      term    = aa * bb + cc + dd
      dtermdt = daa*bb + aa*dbb + dcc + ddd

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 1.070 * exp(-46.6310*t9i)
      drevdt   = rev*46.6310*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_li6ng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd 


c..li6(n,g)li7
c..malaney-fowler 1989

      term    = 5.10e3
      dtermdt = 0.0d0

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 1.19e10 * t932 * exp(-84.17*t9i) 
      drevdt   = rev*(1.5d0*t9i + 84.17*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_he4tg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


c..he4(t,g)li7
      aa  = 8.67e5 * t9i23 * exp(-8.08*t9i13) 
      daa = aa*(-twoth*t9i + oneth*8.08*t9i43)

      bb  = 1.0d0 + 0.052*t913 - 0.448*t923 - 0.165*t9   
     1      + 0.144*t943 + 0.134*t953  
      dbb = oneth*0.052*t9i23 - twoth*0.448*t9i13 - 0.165
     1      + fourth*0.144*t913 + fiveth*0.134*t923 

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 1.11e10 * t932 * exp(-28.64*t9i)
      drevdt   = rev*(1.5d0*t9i + 28.64*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_li7dn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt


c..li7(d,n)2a
      term    = 2.92e11 * t9i23 * exp(-10.259*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*10.259*t9i43)

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 9.95e-10 * t9i32 * exp(-175.476*t9i)
      drevdt   = rev*(-1.5d0*t9i + 175.476*t9i2)

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_li7tn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt


c..li7(t,n)be9
c..malaney and fowler (apjl, 345, l5, 1989)
      term    = 1.46d+11 * t9i23 * exp(-11.333*t9i13)
      dtermdt = -twoth*term*t9i + oneth*term*11.333*t9i43

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 0.0d0
      drevdt   = 0.0d0

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_li7t2n(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt


c..li7(t,2n)2a
      term    = 8.81e11 * t9i23 * exp(-11.333*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*11.333*t9i43)

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 1.22e-19 * t9i3 * exp(-102.864*t9i)
      drevdt   = rev*(-3.0d0*t9i + 102.864*t9i2)

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_li7he3np(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt


c..li7(he3,np)2a
      term    = 1.11e13 * t9i23 * exp(-17.989*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*17.989*t9i43)

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 6.09e-20 * t9i3 * exp(-111.727*t9i)
      drevdt   = rev*(-3.0d0*t9i + 111.727*t9i2)

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_li6pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,
     1                 t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz


c..li6(p,g)be7
      if (t9 .gt. 10.0) then
       t9a  = 1.0d0
       dt9a = 0.0d0
      else  
       aa   = 1.0d0 - 0.0969*t9

       bb   = aa**(-twoth)
       dbb  = twoth*bb/aa*0.0969

       cc   = aa + 0.0284*t953*bb
       dcc  = -0.0969 + 0.0284*(fiveth*t923*bb + t953*dbb)

       zz   = 1.0d0/cc
       t9a  = t9*zz
       dt9a = (1.0d0 - t9a*dcc)*zz
      end if

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix*t9a56*zz

      term    = 6.69e+05 * t9a56 * t9i32 * exp(-8.413/t9a13)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i
     1                + 8.413/t9a13**2 * dt9a13)

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 1.19e+10 * t932 * exp(-65.054*t9i)
      drevdt   = rev*(1.5d0*t9i + 65.054*t9i2)

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term


      return
      end




      subroutine rate_li7pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


c..li7(p,n)be7
      aa  = 5.15e+09 * exp(-1.167*t913 - 19.081*t9i) 
      daa = aa*(-oneth*1.167*t9i23 + 19.081*t9i2)

      bb  = 7.84e+09 * t9i32 * exp(-22.832*t9i)
      dbb = -1.5d0*bb*t9i + bb*22.832*t9i2

      term    = aa + bb
      dtermdt = daa + dbb

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      aa  = 5.15e+09 * exp(-1.167*t913) 
      daa = -aa*oneth*1.167*t9i23

      bb  = 0.998 * 7.84e+09 * t9i32 * exp(-3.751*t9i)
      dbb = -1.5d0*bb*t9i + bb*3.751*t9i2

      term    = aa + bb
      dtermdt = daa + dbb
      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_li7ng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa


c..li7(n,g)li8
c..apj 372, 1
      aa      = 4.26d+03 * t9i32 * exp(-2.576*t9i)
      daa     = aa*(-1.5d0*t9i + 2.576*t9i2)

      term    = 3.144d+03 + aa
      dtermdt = daa

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev    = 1.2923d+10 * t932 * exp(-2.359d+01*t9i)
      drevdt = rev*(1.5d0*t9i + 2.359d+01*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_be7dp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt


c..be7(d,p)2a
      term    = 1.07e12 * t9i23 * exp(-12.428*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*12.428*t9i43)

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 9.97e-10 * t9i32 * exp(-194.557*t9i)
      drevdt   = rev*(-1.5d0*t9i + 194.557*t9i2)

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_be7tnp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt


c..be7(t,np)2a
      term    = 2.91e12 * t9i23 * exp(-13.729*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*13.729*t9i43)

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 6.09e-20 * t9i3 * exp(-121.944*t9i)
      drevdt   = rev*(-3.0d0*t9i + 121.944*t9i2)

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end



      subroutine rate_be7he32p(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt


c..be7(he3,2p)2a
      term    = 6.11e13 * t9i23 * exp(-21.793*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*21.793*t9i43)

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 1.22e-19 * t9i3 * exp(-130.807*t9i)
      drevdt   = rev*(-3.0d0*t9i + 130.807*t9i2)

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_be9pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee


c..be9(p,a)li6
      aa  = 2.11e11 * t9i23 * exp(-10.359*t9i13 - t92/0.2704)  
      daa = -twoth*aa*t9i + aa*(oneth*10.359*t9i43 - 2.0d0*t9/0.2704)  

      bb   = 1.0d0 + 0.04*t913 + 1.09*t923 + 0.307*t9 
     1      + 3.21*t943 + 2.30*t953
      dbb  = oneth*0.04*t9i23 + twoth*1.09*t9i13 + 0.307
     1       + fourth*3.21*t913 + fiveth*2.30*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 4.51e8 * t9i * exp(-3.046*t9i)  
      ddd  = -dd*t9i + dd*3.046*t9i2

      ee   = 6.70e8 * t9i34 * exp(-5.160*t9i)
      dee  = -0.75d0*ee*t9i + ee*5.160*t9i2

      term    = cc + dd + ee 
      dtermdt = dcc + ddd + dee 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.18e-1 * exp(-24.674*t9i) 
      drevdt   = rev*24.674*t9i2

      rr    = den * rev * term 
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_li6ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee


c..li6(a,g)b10
      aa  = 4.06e6 * t9i23 * exp(-18.790*t9i13 - t92/1.758276)
      daa = aa*(-twoth*t9i + oneth*18.790*t9i43 - 2.0d0*t9/1.758276)

      bb   = 1.0d0 + 0.022*t913 + 1.54*t923 + 0.239*t9  
     1       +  2.20*t943 + 0.869*t953 
      dbb  = oneth*0.022*t9i23 + twoth*1.54*t9i13 + 0.239
     1       + fourth*2.20*t913 + fiveth*0.869*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1910.0 * t9i32 * exp(-3.484*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.484*t9i2)

      ee   = 1.01e4 * t9i * exp(-7.269*t9i)
      dee  = ee*(-t9i + 7.269*t9i2)

      term    = cc + dd + ee 
      dtermdt = dcc + ddd + dee 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.58e10 * t932 * exp(-51.753*t9i)
      drevdt   = rev*(1.5d0*t9i + 51.753*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_li7an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt

c..li7(a,n)b10
      term    = 3.84e8 * exp(-32.382*t9i)  
      dtermdt = term*32.382*t9i2

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.84e8*1.32
      drevdt   = 0.0d0

      rr    = den * rev 
      drrdt = 0.0d0
      drrdd = rev

      return
      end





      subroutine rate_be9pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee


c..be9(p,g)b10
      aa  = 1.33e7 * t9i23 * exp(-10.359*t9i13 - t92/0.715716)
      daa = aa*(-twoth*t9i + oneth*10.359*t9i43 - 2.0d0*t9/0.715716)

      bb   = 1.0d0 + 0.04*t913 + 1.52*t923 + 0.428*t9  
     1       + 2.15*t943 + 1.54*t953
      dbb  = oneth*0.04*t9i23 + twoth*1.52*t9i13 + 0.428
     1       + fourth*2.15*t913 + fiveth*1.54*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 9.64e4 * t9i32 * exp(-3.445*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.445*t9i2)

      ee   = 2.72e6 * t9i32 * exp(-10.62*t9i)
      dee  = ee*(-1.5d0*t9i + 10.62*t9i2)

      term    = cc + dd + ee 
      dtermdt = dcc + ddd + dee 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.73e9 * t932 * exp(-76.427*t9i)
      drevdt   = rev*(1.5d0*t9i + 76.427*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_b10pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd

c..b10(p,a)li7
      aa  = 1.26e11 * t9i23 * exp(-12.062*t9i13 - t92/19.377604)
      daa = -twoth*aa*t9i + aa*(oneth*12.062*t9i43 - 2.0d0*t9/19.377604)

      bb   = 1.0d0 + 0.035*t913 - 0.498*t923 - 0.121*t9  
     1       + 0.3*t943 + 0.184*t953
      dbb  = oneth*0.035*t9i23 - twoth*0.498*t9i13 - 0.121
     1       + fourth*0.3*t913 + fiveth*0.184*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 2.59e9 * t9i * exp(-12.260*t9i)
      ddd  = -dd*t9i + dd*12.260*t9i2

      term    = cc + dd 
      dtermdt = dcc + ddd 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.54e-01 * exp(-13.301*t9i)
      drevdt   = rev*13.301*t9i2

      rr    = den * rev * term 
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_li7ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee


c..li7(a,g)b11
      aa  = 3.55e7 * t9i23 * exp(-19.161*t9i13 -t92/17.598025)
      daa = aa*(-twoth*t9i + oneth*19.161*t9i43 - 2.0d0*t9/17.598025)

      bb   = 1.0d0 + 0.022*t913 + 0.775*t923 + 0.118*t9
     1       + 0.884*t943 + 0.342*t953
      dbb  = oneth*0.022*t9i23 + twoth*0.775*t9i13 + 0.118
     1       + fourth*0.884*t913 + fiveth*0.342*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 3.33e2 * t9i32 * exp(-2.977*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.977*t9i2)

      ee   = 4.10e4 * t9i * exp(-6.227*t9i)
      dee  = ee*(-t9i + 6.227*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.02e10 * t932 * exp(-100.538*t9i)
      drevdt   = rev*(1.5d0*t9i + 100.538*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_b11pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee,ff,dff

c..b11(p,a)be8=>2a
      aa  = 2.20e12 * t9i23 * exp(-12.095*t9i13 - t92/2.702736)
      daa = aa*(-twoth*t9i + oneth*12.095*t9i43 - 2.0d0*t9/2.702736)

      bb   = 1.0d0  + 0.034*t913 + 0.14*t923 + 0.034*t9
     1       + 0.19*t943 + 0.116*t953
      dbb  = oneth*0.034*t9i23 + twoth*0.14*t9i13 + 0.034
     1       + fourth*0.19*t913 + fiveth*0.116*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 4.03e6 * t9i32 * exp(-1.734*t9i)
      ddd  = dd*(-1.5d0*t9i + 1.734*t9i2)

      ee   = 6.73e9 * t9i32 * exp(-6.262*t9i)
      dee  = ee*(-1.5d0*t9i + 6.262*t9i2)

      ff   = 3.88e9*t9i * exp(-14.154*t9i)
      dff  = ff*(-t9i + 14.154*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd +dee + dff

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.5e-10* t9i32 *exp(-100.753*t9i)
      drevdt   = rev*(-1.5d0*t9i + 100.753*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_be7ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee

c..be7(a,g)c11
      aa  = 8.45e+07 * t9i23 * exp(-23.212*t9i13 - t92/22.743361)
      daa = aa*(-twoth*t9i + oneth*23.212*t9i43 - 2.0d0*t9/22.743361)

      bb   = 1.0d0 + 0.018*t913 + 0.488*t923 + 0.061*t9
     1       + 0.296*t943 + 0.095*t953
      dbb  = oneth*0.018*t9i23 + twoth*0.488*t9i13 + 0.061
     1       + fourth*0.296*t913 + fiveth*0.095*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.25e+04 * t9i32 * exp(-6.510*t9i)
      ddd  = dd*(-1.5d0*t9i + 6.510*t9i2)

      ee   = 1.29e+05 * t9i54 * exp(-10.039*t9i)
      dee  = ee*(-1.25d0*t9i + 10.039*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.02e+10 * t932 * exp(-87.539*t9i)
      drevdt   = rev*(1.5d0*t9i + 87.539*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_b11pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


c..b11(p,n)c11
      aa  = 1.69e8*(1.0d0 - 0.048*t912 + 0.010*t9)  
      daa = 1.69e8*(-0.5d0*0.048*t9i12 + 0.010)

      bb  = exp(-32.080*t9i)
      dbb = bb*32.080*t9i2

      term    = aa*bb
      dtermdt = daa*bb + aa*dbb

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      term     = 0.998*aa
      dtermdt  = 0.998*daa

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_b8ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,
     2                 pp,dpp,qq,dqq,rr,drr 


c..b8(a,p)c11
      aa  = 1.67e-09 * t9i32 * exp(-1.079*t9i) 
      daa = -1.5d0*aa*t9i + aa*1.079*t9i2

      bb  = 9.55e-08 * t9i32 * exp(-1.311*t9i)
      dbb = -1.5d0*bb*t9i + bb*1.311*t9i2

      cc  = 1.98e-01 * t9i32 * exp(-2.704*t9i) 
      dcc = -1.5d0*cc*t9i + cc*2.704*t9i2

      dd  = 1.34e+00 * t9i32 * exp(-4.282*t9i)
      ddd = -1.5d0*dd*t9i + dd*4.282*t9i2

      ee  = 3.22e+04 * t9i32 * exp(-6.650*t9i) 
      dee = -1.5d0*ee*t9i + ee*6.650*t9i2

      ff  = 2.33e+05 * t9i32 * exp(-8.123*t9i)
      dff = -1.5d0*ff*t9i + ff*8.123*t9i2

      gg  = 2.55e+06 * t9i32 * exp(-11.99*t9i) 
      dgg = -1.5d0*gg*t9i + gg*11.99*t9i2

      hh  = 9.90e+06 * t9i32 * exp(-13.50*t9i)
      dhh = -1.5d0*hh*t9i + hh*13.50*t9i2

      pp  = 1.41e+06 * t9i32 * exp(-16.51*t9i) 
      dpp = -1.5d0*pp*t9i + pp*16.51*t9i2

      qq  = 1.99e+07 * t9i32 * exp(-18.31*t9i)
      dqq = -1.5d0*qq*t9i + qq*18.31*t9i2

      rr  = 6.01e+07 * t9i32 * exp(-20.63*t9i) 
      drr = -1.5d0*rr*t9i + rr*20.63*t9i2

      term    = aa + bb + cc + dd + ee + ff + gg + hh + pp + qq + rr
      dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg + dhh 
     1          + dpp + dqq + drr

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 3.101 * exp(-85.95*t9i)
      drevdt   = rev*85.95*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev*term

      return
      end





      subroutine rate_b10pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee


c..b10(p,g)c11
      aa  = 4.61e+05 * t9i23 * exp(-12.062*t9i13 - t92/19.377604)
      daa = aa*(-twoth*t9i + oneth*12.062*t9i43 - 2.0d0*t9/19.377604)

      bb   = 1.0d0 + 0.035*t913 + 0.426*t923 + 0.103*t9
     1       + 0.281*t943 + 0.173*t953
      dbb  = oneth*0.035*t9i23 + twoth*0.426*t9i13 + 0.103
     1       + fourth*0.281*t913 + fiveth*0.173*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.93e+05 * t9i32 * exp(-12.041*t9i)
      ddd  = dd*(-1.5d0*t9i + 12.041*t9i2)

      ee   = 1.14e+04 * t9i32 * exp(-16.164*t9i)
      dee  = ee*(-1.5d0*t9i + 16.164*t9i2)

      term    = cc + dd + ee 
      dtermdt = dcc + ddd + dee 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.03e+10 * t932 * exp(-100.840*t9i)
      drevdt   = rev*(1.5d0*t9i + 100.840*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_c11na(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd


c..c11(n,a)be8=>2a
c..hauser feshbach calculation by woosley on aug 26, 1988.

      fr    = den * 7.0e4
      dfrdt = 0.0d0
      dfrdd = 7.0e4

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_be9an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh


c..be9(a,n)c12
c..Wrean 94 Phys Rev C (1994) vol 49, #2, 1205
      aa  = 6.476d+13 * t9i23 * exp(-23.8702*t9i13)
      daa = -twoth*aa*t9i + oneth*aa*23.8702*t9i43

      bb  = (1.0d0 - 0.3270*t913)
      dbb = -oneth*0.3270*t9i23

c..rate goes negative for t9 greater than about 15, so try this 
      if (bb .le. 0.0) then
       bb  = 0.0d0
       dbb = 0.0d0
      end if

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 6.044d-3*t9i32*exp(-1.401*t9i) 
      ddd = -1.5d0*dd*t9i + dd*1.401*t9i2

      ee  = 7.268*t9i32*exp(-2.063*t9i) 
      dee = -1.5d0*ee*t9i + ee*2.063*t9i2

      ff  = 3.256d+4*t9i32*exp(-3.873*t9i) 
      dff = -1.5d0*ff*t9i + ff*3.873*t9i2

      gg  = 1.946d+5*t9i32*exp(-4.966*t9i) 
      dgg = -1.5d0*gg*t9i + gg*4.966*t9i2

      hh  = 1.838e9*t9i32*exp(-15.39*t9i) 
      dhh = -1.5d0*hh*t9i + hh*15.39*t9i2

      term    = cc + dd + ee + ff + gg + hh
      dtermdt = dcc + ddd + dee + dff + dgg + dhh


c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term 

      rev      = 10.3 * exp(-66.160*t9i)
      drevdt   = rev*66.160*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term


      return
      end




      subroutine rate_b11pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee


c..b11(p,g)c12
      aa  = 4.62e+07 * t9i23 * exp(-12.095*t9i13 - t92/0.57121)
      daa = aa*(-twoth*t9i + oneth*12.095*t9i43 - 2.0d0*t9/0.57121)

      bb   = 1.0d0 + 0.035*t913 + 3.0*t923 + 0.723*t9
     1       + 9.91*t943 + 6.07*t953
      dbb  = oneth*0.035*t9i23 + twoth*3.0*t9i13 + 0.723
     1       + fourth*9.91*t913 + fiveth*6.07*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 7.89e+03 * t9i32 * exp(-1.733*t9i)
      ddd  = dd*(-1.5d0*t9i + 1.733*t9i2)

      ee   = 9.68e+04 * t9i15 * exp(-5.617*t9i)
      dee  = ee*(-0.2d0*t9i + 5.617*t9i2)

      term    = cc + dd + ee 
      dtermdt = dcc + ddd + dee 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.01e+10 * t932 * exp(-185.173*t9i)
      drevdt   = rev*(1.5d0*t9i + 185.173*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_b11ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee,ff,dff


c..b11(a,p)c14
      aa  = 5.37e+11 * t9i23 * exp(-28.234*t9i13 - t92/0.120409)
      daa = aa*(-twoth*t9i + oneth*28.234*t9i43 - 2.0d0*t9/0.120409)

      bb   = 1.0d0 + 0.015*t913 + 5.575*t923 + 0.576*t9
     1       + 15.888*t943 + 4.174*t953
      dbb  = oneth*0.015*t9i23 + twoth*5.575*t9i13 + 0.576
     1       + fourth*15.888*t913 + fiveth*4.174*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 5.44e-03 * t9i32 * exp(-2.827*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.827*t9i2)

      ee   = 3.36e+02 * t9i32 * exp(-5.178*t9i)
      dee  = ee*(-1.5d0*t9i + 5.178*t9i2)

      ff   = 5.32e+06 * t9i38 * exp(-11.617*t9i)
      dff  = ff*(-0.375*t9i + 11.617*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.10e+01 * exp(-9.098*t9i)
      drevdt   = rev*9.098*t9i2

      rr    = den * rev * term 
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_pp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,bb,dbb


c..p(p,e+nu)d
      if (t9 .le. 3.0) then
       aa   = 4.01d-15 * t9i23 * exp(-3.380d0*t9i13) 
       daa  = aa*(-twoth*t9i + oneth*3.380d0*t9i43) 

       bb   = 1.0d0 + 0.123d0*t913 + 1.09d0*t923 + 0.938d0*t9 
       dbb  = oneth*0.123d0*t9i23 + twoth*1.09d0*t9i13 + 0.938d0

       term    = aa * bb
       dtermdt = daa * bb + aa * dbb

      else
       term    = 1.1581136d-15
       dtermdt = 0.0d0
      end if

c..rate 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_png(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa


c..p(n,g)d
c..smith,kawano,malany 1992

       aa      = 1.0d0 - 0.8504*t912 + 0.4895*t9
     1           - 0.09623*t932 + 8.471e-3*t92 
     2           - 2.80e-4*t952

       daa     =  -0.5d0*0.8504*t9i12 + 0.4895
     1           - 1.5d0*0.09623*t912 + 2.0d0*8.471e-3*t9 
     2           - 2.5d0*2.80e-4*t932
  
       term    = 4.742e4 * aa
       dtermdt = 4.742e4 * daa 


c..wagoner,schramm 1977
c      aa      = 1.0d0 - 0.86*t912 + 0.429*t9
c      daa     =  -0.5d0*0.86*t9i12 + 0.429

c      term    = 4.4d4 * aa
c      dtermdt = 4.4d4 * daa



c..rates 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.71d+09 * t932 * exp(-25.82*t9i)
      drevdt   = rev*(1.5d0*t9i + 25.82*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_dpg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


c..d(p,g)he3
      aa      = 2.24d+03 * t9i23 * exp(-3.720*t9i13) 
      daa     = aa*(-twoth*t9i + oneth*3.720*t9i43)

      bb      = 1.0d0 + 0.112*t913 + 3.38*t923 + 2.65*t9
      dbb     = oneth*0.112*t9i23 + twoth*3.38*t9i13 + 2.65

      term    = aa * bb 
      dtermdt = daa * bb + aa * dbb  


c..rates 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.63d+10 * t932 * exp(-63.750*t9i)
      drevdt   = rev*(1.5d0*t9i + 63.750*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_he3ng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt


c..he3(n,g)he4
      term    = 6.62 * (1.0d0 + 905.0*t9)
      dtermdt = 5.9911d3

c..rates 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.61d+10 * t932 * exp(-238.81*t9i)
      drevdt   = rev*(1.5d0*t9i + 238.81*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_he3he3(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


c..he3(he3,2p)he4
      aa   = 6.04d+10 * t9i23 * exp(-12.276*t9i13) 
      daa  = aa*(-twoth*t9i + oneth*12.276*t9i43)
      
      bb   = 1.0d0 + 0.034*t913 - 0.522*t923 - 0.124*t9 
     1       + 0.353*t943 + 0.213*t953
      dbb  = oneth*0.034*t9i23 - twoth*0.522*t9i13 - 0.124
     1       + fourth*0.353*t913 + fiveth*0.213*t923

      term    = aa * bb 
      dtermdt = daa*bb + aa*dbb

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.39e-10 * t9i32 * exp(-149.230*t9i)
      drevdt   = rev*(-1.5d0*t9i + 149.230*t9i2)

      rr    = den * rev * term 
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term 

      return
      end





      subroutine rate_he3he4(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,t9a,dt9a,
     1                 t9a13,dt9a13,t9a56,dt9a56,zz


c..he3(he4,g)be7
      aa      = 1.0d0 + 0.0495*t9
      daa     = 0.0495    

      zz      = 1.0d0/aa
      t9a     = t9*zz
      dt9a    = (1.0d0 - t9a*daa)*zz

      zz      = dt9a/t9a 
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56   = t9a**fivsix
      dt9a56  = fivsix*t9a56*zz

      term    = 5.61d+6 * t9a56 * t9i32 * exp(-12.826/t9a13)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i
     1                + 12.826/t9a13**2 * dt9a13)

c..rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.11e+10 * t932 * exp(-18.423*t9i)
      drevdt   = rev*(1.5d0*t9i + 18.423*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_be7em(temp,den,ye,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,ye,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,bb,dbb

c..be7(e-,nu+g)li7
      if (t9 .le. 3.0) then
       aa  = 0.0027 * t9i * exp(2.515e-3*t9i) 
       daa = -aa*t9i - aa*2.515e-3*t9i2

       bb  = 1.0d0 - 0.537*t913 + 3.86*t923 + aa
       dbb = -oneth*0.537*t9i23 + twoth*3.86*t9i13 + daa

       term    = 1.34e-10 * t9i12 * bb
       dtermdt = -0.5d0*term*t9i + 1.34e-10*t9i12*dbb

      else
       term    = 0.0d0
       dtermdt = 0.0d0  
      endif

c..rates
      fr    = ye * den * term
      dfrdt = ye * den * dtermdt * 1.0d-9
      dfrdd = ye * term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_be7pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


c..be7(p,g)b8
      aa      = 3.11e+05 * t9i23 * exp(-10.262*t9i13)  
      daa     = aa*(-twoth*t9i + oneth*10.262*t9i43)

      bb      = 2.53e+03 * t9i32 * exp(-7.306*t9i)  
      dbb     = bb*(-1.5d0*t9i + 7.306*t9i2)

      term    = aa + bb 
      dtermdt = daa + dbb 


c..rates 
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.30e+10 * t932 * exp(-1.595*t9i)
      drevdt   = rev*(1.5d0*t9i + 1.595*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_li7pag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,
     1                 t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz,
     2                 term1,dterm1,term2,dterm2,rev,drevdt



c..7li(p,g)8be=>2a  
      aa   = 1.56e+05 * t9i23 * exp(-8.472*t9i13 - t92/2.876416)
      daa  = aa*(-twoth*t9i + oneth*8.472*t9i43 - 2.0d0*t9/2.876416)

      bb   = 1.0d0 + 0.049*t913 + 2.498*t923 + 0.86*t9
     1       + 3.518*t943 + 3.08*t953
      dbb  = oneth*0.049*t9i23 + twoth*2.498*t9i13 + 0.86
     1       + fourth*3.518*t913 + fiveth*3.08*t923

      cc   = aa*bb
      dcc  = daa*bb + aa*dbb

      dd   =  1.55e+06 * t9i32 * exp(-4.478*t9i)
      ddd  = dd*(-1.5d0*t9i + 4.478*t9i2)

      term1  = cc + dd
      dterm1 = dcc + ddd

      rev    = 6.55e+10 * t932 * exp(-200.225*t9i)
      drevdt = rev*(1.5d0*t9i + 200.225*t9i2) 


c..7li(p,a)a 
      aa     = 1.0d0 + 0.759*t9

      zz     = 1.0d0/aa
      t9a    = t9*zz
      dt9a   = (1.0d0 - t9a*0.759)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix*t9a56*zz

      aa     = 1.096e+09 * t9i23 * exp(-8.472*t9i13)
      daa    = aa*(-twoth*t9i + oneth*8.472*t9i43)

      bb     = -4.830e+08 * t9a56 * t9i32 * exp(-8.472/t9a13)
      dbb    = bb*(dt9a56/t9a56 - 1.5d0*t9i + 8.472/t9a13**2*dt9a13) 
    
      cc     = 1.06e+10 * t9i32 * exp(-30.442*t9i) 
      dcc    = cc*(-1.5d0*t9i + 30.442*t9i2) 

      term2   = aa + bb + cc
      dterm2  = daa + dbb + dcc 

      rev    = 4.69 * exp(-201.291*t9i)
      drevdt = aa*201.291*t9i2 


c..sum of these two rates
      term = term1 + term2
      dtermdt = dterm1 + dterm2 


c..rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_b8ep(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0,
     1                  halflife = 0.77d0,
     2                  con      = lntwo/halflife)


c..b8(e+,nu)be8 => 2a

      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_c12pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee


c..c12(p,g)13n
      aa   = 2.04e+07 * t9i23 * exp(-13.69*t9i13 - t92/2.25)
      daa  = aa*(-twoth*t9i + oneth*13.69*t9i43 - 2.0d0*t9/2.25)

      bb   = 1.0d0 + 0.03*t913 + 1.19*t923 + 0.254*t9 
     1       + 2.06*t943 + 1.12*t953
      dbb  = oneth*0.03*t9i23 + twoth*1.19*t9i13 + 0.254
     1       + fourth*2.06*t913 + fiveth*1.12*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.08e+05 * t9i32 * exp(-4.925*t9i)
      ddd  = dd*(-1.5d0*t9i + 4.925*t9i2)

      ee   = 2.15e+05 * t9i32 * exp(-18.179*t9i)
      dee  = ee*(-1.5d0*t9i + 18.179*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 8.84e+09 * t932 * exp(-22.553*t9i)
      drevdt   = rev*(1.5d0*t9i + 22.553*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_n13em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0,
     1                  halflife = 597.9d0,
     2                  con      = lntwo/halflife)

c..n13(e-nu)c13
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_c13pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd


c..c13(p,g)13n
      aa   = 8.01e+07 * t9i23 * exp(-13.717*t9i13 - t92/4.0)
      daa  = aa*(-twoth*t9i + oneth*13.717*t9i43 - 0.5d0*t9)

      bb   = 1.0d0 + 0.030*t913 + 0.958*t923 + 0.204*t9
     1       + 1.39*t943 + 0.753*t953
      dbb  = oneth*0.030*t9i23 + twoth*0.958*t9i13 + 0.204
     1       + fourth*1.39*t913 + fiveth*0.753*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.21e+06 * t9i65 * exp(-5.701*t9i)
      ddd  = dd*(-sixfif*t9i + 5.701*t9i2)

      term    = cc + dd 
      dtermdt = dcc + ddd 


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.19e+10 * t932 * exp(-87.621*t9i)
      drevdt   = rev*(1.5d0*t9i + 87.621*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_n14pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee


c..n14(p,g)o15
      aa  = 4.90e+07 * t9i23 * exp(-15.228*t9i13 - t92/10.850436)
      daa = aa*(-twoth*t9i + oneth*15.228*t9i43 - 2.0d0*t9/10.850436)

      bb   = 1.0d0 + 0.027*t913 - 0.778*t923 - 0.149*t9 
     1       + 0.261*t943 + 0.127*t953
      dbb  = oneth*0.027*t9i23 - twoth*0.778*t9i13 - 0.149
     1       + fourth*0.261*t913 + fiveth*0.127*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 2.37e+03 * t9i32 * exp(-3.011*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.011*t9i2)

      ee   = 2.19e+04 * exp(-12.530*t9i)
      dee  = ee*12.530*t9i2

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 2.70e+10 * t932 * exp(-84.678*t9i)
      drevdt = rev*(1.5d0*t9i + 84.678*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_o15em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0,
     1                  halflife = 122.24d0,
     2                  con      = lntwo/halflife)

c..o15(e-nu)n15
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_n14ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee,ff,dff


c..n14(a,g)f18
      aa  = 7.78d+09 * t9i23 * exp(-36.031*t9i13- t92/0.776161) 
      daa = aa*(-twoth*t9i + oneth*36.031*t9i43 - 2.0d0*t9/0.776161)

      bb   = 1.0d0 + 0.012*t913 + 1.45*t923 + 0.117*t9 
     1       + 1.97*t943 + 0.406*t953
      dbb  = oneth*0.012*t9i23 + twoth*1.45*t9i13 + 0.117
     1       + fourth*1.97*t913 + fiveth*0.406*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 2.36d-10 * t9i32 * exp(-2.798*t9i)  
      ddd  = dd*(-1.5d0*t9i + 2.798*t9i2)

      ee   = 2.03 * t9i32 * exp(-5.054*t9i) 
      dee  = ee*(-1.5d0*t9i + 5.054*t9i2)

      ff   = 1.15d+04 * t9i23 * exp(-12.310*t9i)
      dff  = ff*(-twoth*t9i + 12.310*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.42e+10 * t932 * exp(-51.236*t9i)
      drevdt   = rev*(1.5d0*t9i + 51.236*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_n15pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee,ff,dff


c..n15(p,g)o16
      aa  = 9.78e+08 * t9i23 * exp(-15.251*t9i13 - t92/0.2025)
      daa = aa*(-twoth*t9i + oneth*15.251*t9i43 - 2.0d0*t9/0.2025)

      bb   = 1.0d0  + 0.027*t913 + 0.219*t923 + 0.042*t9 
     1       + 6.83*t943 + 3.32*t953
      dbb  = oneth*0.027*t9i23 + twoth*0.219*t9i13 + 0.042
     1       + fourth*6.83*t913 + fiveth*3.32*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.11e+04*t9i32*exp(-3.328*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.328*t9i2)

      ee   = 1.49e+04*t9i32*exp(-4.665*t9i)
      dee  = ee*(-1.5d0*t9i + 4.665*t9i2)

      ff   = 3.8e+06*t9i32*exp(-11.048*t9i)
      dff  = ff*(-1.5d0*t9i + 11.048*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.62e+10 * t932 * exp(-140.734*t9i)
      drevdt   = rev*(1.5d0*t9i + 140.734*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_n15pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,
     2                 theta
      parameter        (theta = 0.1d0)


c..n15(p,a)c12
      aa  = 1.08d+12*t9i23*exp(-15.251*t9i13 - t92/0.272484)
      daa = aa*(-twoth*t9i + oneth*15.251*t9i43 - 2.0d0*t9/0.272484)

      bb   = 1.0d0 + 0.027*t913 + 2.62*t923 + 0.501*t9
     1       + 5.36*t943 + 2.60*t953
      dbb  = oneth*0.027*t9i23 + twoth*2.62*t9i13 + 0.501
     1       + fourth*5.36*t913 + fiveth*2.60*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.19d+08 * t9i32 * exp(-3.676*t9i) 
      ddd  = dd*(-1.5d0*t9i + 3.676*t9i2)

      ee   = 5.41d+08 * t9i12 * exp(-8.926*t9i)
      dee  = ee*(-0.5d0*t9i + 8.926*t9i2)

      ff   = theta * 4.72d+08 * t9i32 * exp(-7.721*t9i) 
      dff  = ff*(-1.5d0*t9i + 7.721*t9i2)

      gg   = theta * 2.20d+09 * t9i32 * exp(-11.418*t9i)
      dgg  = gg*(-1.5d0*t9i + 11.418*t9i2)

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.06d-01*exp(-57.625*t9i)
      drevdt   = rev*57.625*t9i2

      rr    = den * rev * term 
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_o16pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee,zz


c..o16(p,g)f17
      aa  = exp(-0.728*t923)
      daa = -twoth*aa*0.728*t9i13

      bb  = 1.0d0 + 2.13 * (1.0d0 - aa)
      dbb = -2.13*daa

      cc  = t923 * bb
      dcc = twoth*cc*t9i + t923*dbb

      dd   = exp(-16.692*t9i13)
      ddd  = oneth*dd*16.692*t9i43

      zz   = 1.0d0/cc
      ee   = dd*zz
      dee  = (ddd - ee*dcc)*zz

      term    = 1.50d+08 * ee
      dtermdt = 1.50d+08 * dee


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.03e+09*t932*exp(-6.968*t9i)
      drevdt   = rev*(1.5d0*t9i + 6.968*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_o17pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,res1,dres1,res2,dres2,res3,dres3,
     2                 res4,dres4,res5,dres5,res6,dres6,zz,
     3                 theta
      parameter        (theta = 0.1d0)



c..o17(p,a)n14
c..rate from jeff blackmons thesis, includes terms from fowler 75,
c..landre 1990 (a&a 240, 85), and new results
c..use rev factor from cf88 rate

      aa  = 1.53d+07 * t9i23 * exp(-16.712*t9i13 - t92/0.319225)
      daa = aa*(-twoth*t9i + oneth*16.712*t9i43 - 2.0d0*t9/0.319225)

      bb   = 1.0d0 + 0.025*t913 + 5.39*t923 + 0.940*t9
     1       + 13.5*t943 + 5.98*t953
      dbb  = oneth*0.025*t9i23 + twoth*5.39*t9i13 + 0.940
     1       + fourth*13.5*t913 + fiveth*5.98*t923

      res1  = aa * bb
      dres1 = daa*bb + aa*dbb

      res2  = 2.92d+06 * t9 * exp(-4.247*t9i)
      dres2 = res2*(t9i + 4.247*t9i2)


      aa    = 0.479 * t923 + 0.00312
      daa   = twoth*0.479*t9i13

      bb    = aa*aa
      dbb   = 2.0d0 * aa * daa

      cc    =  1.78d+05 * t9i23 * exp(-16.669*t9i13)
      dcc   = cc*(-twoth*t9i + oneth*16.669*t9i43)

      zz    = 1.0d0/bb 
      res3  = cc*zz
      dres3 = (dcc - res3*dbb)*zz

      res4  = 8.68d+10 * t9 * exp(-16.667*t9i13 - t92/0.0016)
      dres4 = res4*(t9i + oneth*16.667*t9i43 - 2.0d0*t9/0.0016)

      res5  = 9.22d-04 * t9i32 * exp(-0.767*t9i)
      dres5 = res5*(-1.5d0*t9i + 0.767*t9i2)

      res6  = theta * 98.0 * t9i32 * exp(-2.077*t9i)
      dres6 = res6*(-1.5d0*t9i + 2.077*t9i2)

      term    = res1 + res2 + res3 + res4 + res5 + res6
      dtermdt = dres1 + dres2 + dres3 + dres4 + dres5 + dres6

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.676 * exp(-13.825*t9i)
      drevdt   = rev*13.825*t9i2

      rr    = den * rev * term 
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_o17pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,
     2                 t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,
     3                 zz,theta
      parameter        (theta = 0.1d0)


c..o17(p,g)18f
c..from landre et al 1990 a&a 240, 85
      aa     = 1.0d0 + 2.69*t9
      zz     = 1.0d0/aa
 
      t9a    = t9*zz
      dt9a   = (1.0d0 - t9a*2.69)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix*t9a56*zz

      aa  = 7.97d+07 * t9a56 * t9i32 * exp(-16.712/t9a13)
      daa = aa*(dt9a56/t9a56 - 1.5d0*t9i + 16.712/t9a13**2*dt9a13)

      bb  = 1.0d0  + 0.025*t913 - 0.051*t923 - 8.82d-3*t9
      dbb = oneth*0.025*t9i23 - twoth*0.051*t9i13 - 8.82d-3
      if (bb .le. 0.0) then
       bb  = 0.0d0
       dbb = 0.0d0
      end if

      cc  = 1.51d+08 * t9i23 * exp(-16.712*t9i13)
      dcc = cc*(-twoth*t9i + oneth*16.712*t9i43)

      dd  = bb*cc
      ddd = dbb*cc + bb*dcc  

      ee  = 1.56d+5 * t9i * exp(-6.272*t9i)
      dee = ee*(-t9i + 6.272*t9i2)

      ff  = 2.0d0 * theta * 3.16d-05 * t9i32 * exp(-0.767*t9i)
      dff = ff*(-1.5d0*t9i + 0.767*t9i2)

      gg  = theta * 98.0 * t9i32 * exp(-2.077*t9i)
      dgg = gg*(-1.5d0*t9i + 2.077*t9i2)

      term    = aa + dd + ee + ff + gg
      dtermdt = daa + ddd + dee + dff + dgg

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.66d+10 * t932 * exp(-65.061*t9i)
      drevdt   = rev*(1.5d0*t9i + 65.061*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_o18pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg


c..o18(p,a)n15
      aa  = 3.63e+11 * t9i23 * exp(-16.729*t9i13 - t92/1.852321)
      daa = -twoth*aa*t9i + aa*(oneth*16.729*t9i43 - 2.0d0*t9/1.852321)

      bb  = 1.0d0 + 0.025*t913 + 1.88*t923 + 0.327*t9
     1      + 4.66*t943 + 2.06*t953
      dbb = oneth*0.025*t9i23 + twoth*1.88*t9i13 + 0.327
     1      + fourth*4.66*t913 + fiveth*2.06*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 9.90e-14 * t9i32 * exp(-0.231*t9i)
      ddd = -1.5d0*dd*t9i + dd*0.231*t9i2

      ee  = 2.66e+04 * t9i32 * exp(-1.670*t9i)
      dee = -1.5d0*ee*t9i + ee*1.670*t9i2

      ff  = 2.41e+09 * t9i32 * exp(-7.638*t9i)
      dff = -1.5d0*ff*t9i + ff*7.638*t9i2

      gg  = 1.46e+09 * t9i * exp(-8.310*t9i)
      dgg = -gg*t9i + gg*8.310*t9i2

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.66e-01 * exp(-46.191*t9i)
      drevdt   = rev*46.191*t9i2

      rr    = den * rev * term 
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_o18pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee,ff,dff


c..o18(p,g)19f
      aa  = 3.45e+08 * t9i23 * exp(-16.729*t9i13 - t92/0.019321)
      daa = aa*(-twoth*t9i + oneth*16.729*t9i43 - 2.0d0*t9/0.019321)

      bb  = 1.0d0 + 0.025*t913 + 2.26*t923 + 0.394*t9
     1      + 30.56*t943 + 13.55*t953
      dbb = oneth*0.025*t9i23 + twoth*2.26*t9i13 + 0.394
     1      + fourth*30.56*t913 + fiveth*13.55*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb  

      dd  = 1.25e-15 * t9i32 * exp(-0.231*t9i) 
      ddd = dd*(-1.5d0*t9i + 0.231*t9i2)

      ee  = 1.64e+02 * t9i32 * exp(-1.670*t9i)
      dee = ee*(-1.5d0*t9i + 1.670*t9i2)

      ff  = 1.28e+04 * t912 * exp(-5.098*t9i)
      dff = ff*(0.5d0*t9i + 5.098*t9i2)

      term    = cc + dd + ee + ff 
      dtermdt = dcc + ddd + dee + dff 


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.20e+09 * t932 * exp(-92.769*t9i)
      drevdt   = rev*(1.5d0*t9i + 92.769*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_f17em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0,
     1                  halflife = 64.49d0,
     2                  con      = lntwo/halflife)

c..f17(e-nu)o17
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_f18em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0,
     1                  halflife = 6586.2d0,
     2                  con      = lntwo/halflife)

c..f18(e-nu)o18
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_f19pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh


c..f19(p,a)o16
      aa  = 3.55d+11 * t9i23 * exp(-18.113*t9i13 - t92/0.714025)
      daa = -twoth*aa*t9i + aa*(oneth*18.113*t9i43 - 2.0d0*t9/0.714025)

      bb  = 1.0d0 + 0.023*t913 + 1.96*t923 + 0.316*t9
     1      + 2.86*t943 + 1.17*t953
      dbb = oneth*0.023*t9i23 + twoth*1.96*t9i13 + 0.316
     1      + fourth*2.86*t913 + fiveth*1.17*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 3.67d+06 * t9i32 * exp(-3.752*t9i) 
      ddd = -1.5d0*dd*t9i + dd*3.752*t9i2

      ee  = 3.07d+08 * exp(-6.019*t9i)
      dee = ee*6.019*t9i2

      ff  = 4.0*exp(-2.090*t9i) 
      dff = ff*2.090*t9i2

      gg  = 7.0*exp(-16.440*t9i)
      dgg = gg*16.440*t9i2

      hh  = 1.0d0 + ff + gg
      dhh = dff + dgg

      term    = (cc + dd + ee)/hh
      dtermdt = ((dcc + ddd + dee) - term*dhh)/hh


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.54e-01 * exp(-94.159*t9i)
      drevdt   = rev*94.159*t9i2

      rr    = den * rev * term 
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_n13pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd


c..n13(p,g)o14
c..Keiner et al 1993 Nucl Phys A552, 66
      aa  = -1.727d+7 * t9i23 * exp(-15.168*t9i13 - t92/0.69288976)
      daa = aa*(-twoth*t9i + oneth*15.168*t9i43 -2.0d0*t9/0.69288976)

      bb  = 1.0d0 + 0.027*t913 - 17.54*t923 - 3.373*t9
     1      + 0.0176*t943 + 0.766d-2*t953
      dbb = oneth*0.027*t9i23 - twoth*17.54*t9i13 - 3.373
     1      + fourth*0.0176*t913 + fiveth*0.766d-2*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb  

      dd  = 3.1d+05 * t9i32 * exp(-6.348*t9i)
      ddd = dd*(-1.5d0*t9i + 6.348*t9i2)

      term    = cc + dd
      dtermdt = dcc + ddd

c..goes negative below about t7=1.5
c..note cf88 rate stays positive
      if (term .lt. 0.0) then
       term    = 0.0d0
       dtermdt = 0.0d0
      end if 


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.57d+10*t932*exp(-53.706*t9i)
      drevdt   = rev*(1.5d0*t9i + 53.706*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_o14em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0,
     1                  halflife = 70.606d0,
     2                  con      = lntwo/halflife)

c..f18(e-nu)o18
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_o14ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee,ff,dff


c..o14(a,p)f17
      aa  = 1.68e+13 * t9i23 * exp(-39.388*t9i13- t92/0.514089)
      daa = -twoth*aa*t9i + aa*(oneth*39.388*t9i43 - 2.0d0*t9/0.514089)

      bb  = 1.0d0 + 0.011*t913 + 13.117*t923 + 0.971*t9 
     1      + 85.295*t943 + 16.061*t953
      dbb = oneth*0.011*t9i23 + twoth*13.117*t9i13 + 0.971
     1      + fourth*85.295*t913 + fiveth*16.061*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 3.31e+04 * t9i32 * exp(-11.733*t9i)
      ddd = -1.5d0*dd*t9i + dd*11.733*t9i2

      ee  = 1.79e+07 * t9i32 * exp(-22.609*t9i) 
      dee = -1.5d0*ee*t9i + ee*22.609*t9i2

      ff  = 9.00e+03 * t9113 * exp(-12.517*t9i)
      dff = elvnth*ff*t9i + ff*12.517*t9i2

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff     


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.93e-01*exp(-13.820*t9i)
      drevdt   = rev*13.820*t9i2

      rr    = den * rev * term 
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_o15ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh


c..o15(a,g)ne19

      aa  = 3.57d+11 * t9i23 * exp(-39.584d+0*t9i13 - t92/9.0d0)
      daa = aa*(-twoth*t9i + oneth*39.584d0*t9i43 - 2.0d0*t9/9.0d0)

      bb  = 1.0d0 + 0.011*t913 - 0.273*t923 - 0.020*t9
      dbb = oneth*0.011*t9i23 - twoth*0.273*t9i13 - 0.020

      cc  = aa*bb
      dcc = daa*bb + aa*dbb  

      dd  = 5.10d+10 * t9i23 * exp(-39.584d+0*t9i13 - t92/3.751969)
      ddd = dd*(-twoth*t9i + oneth*39.584*t9i43 - 2.0d0*t9/3.751969)

      ee  = 1.0d0 + 0.011*t913 + 1.59*t923 + 0.117*t9
     1      + 1.81*t943 + 0.338*t953
      dee = oneth*0.011*t9i23 + twoth*1.59*t9i13 + 0.117
     1      + fourth*1.81*t913 + fiveth*0.338*t923

      ff  = dd*ee
      dff = ddd*ee + dd*dee  

      gg  = 3.95d-1 * t9i32 * exp(-5.849*t9i)
      dgg = gg*(-1.5d0*t9i + 5.849*t9i2)

      hh  = 1.90d+1 * t9**2.85 * exp(-7.356*t9i - t92/64.0)
      dhh = hh*(2.85*t9i + 7.356*t9i2 - t9/32.0)


      term    = cc + ff + gg + hh
      dtermdt = dcc + dff + dgg + dhh


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.54e+10 * t932 * exp(-40.957*t9i)
      drevdt   = rev*(1.5d0*t9i + 40.957*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0
 
      return
      end






      subroutine rate_f17pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee


c..f17(p,g)ne18
c..wiescher and kettner, ap. j., 263, 891 (1982)

      aa  = 1.66e+07 * t9i23 * exp(-18.03*t9i13)
      daa = aa*(-twoth*t9i + oneth*18.03*t9i43)

      bb  = 2.194 + 0.050*t913 - 0.376*t923 - 0.061*t9
     1      + 0.026*t943 + 0.011*t953
      dbb = oneth*0.050*t9i23 - twoth*0.376*t9i13 - 0.061
     1      + fourth*0.026*t913 + fiveth*0.011*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb  

      dd  = 839.0 * t9i32 * exp(-6.93*t9i)  
      ddd = dd*(-1.5d0*t9i + 6.93*t9i2)

      ee  = 33.56 * t9i32 * exp(-7.75*t9i)
      dee = ee*(-1.5d0*t9i + 7.75*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.087e+11 * t932 * exp(-45.501*t9i)
      drevdt   = rev*(1.5d0*t9i + 45.501*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_ne18em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0,
     1                  halflife = 1.672d0,
     2                  con      = lntwo/halflife)

c..ne18(e-nu)f18
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_f18pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee,ff,dff


c..f18(p,a)o15
c..wiescher and kettner, ap. j., 263, 891 (1982)

      aa  = 1.66e-10 * t9i32 * exp(-0.302*t9i)
      daa = aa*(-1.5d0*t9i + 0.302*t9i2)

      bb  = 1.56e+05 * t9i32 * exp(-3.84*t9i)
      dbb = bb*(-1.5d0*t9i + 3.84*t9i2)

      cc  = 1.36e+06 * t9i32 * exp(-5.22*t9i)
      dcc = cc*(-1.5d0*t9i + 5.22*t9i2)

      dd  = 8.1e-05 * t9i32 * exp(-1.05*t9i)
      ddd = dd*(-1.5d0*t9i + 1.05*t9i2)

      ee  = 8.9e-04 * t9i32 * exp(-1.51*t9i)
      dee = ee*(-1.5d0*t9i + 1.51*t9i2)

      ff  = 3.0e+05 * t9i32 * exp(-4.29*t9i)
      dff = ff*(-1.5d0*t9i + 4.29*t9i2)

      term    = aa + bb + cc + dd + ee + ff
      dtermdt = daa + dbb + dcc + ddd + dee + dff


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.93e-01 * exp(-33.433*t9i)
      drevdt   = rev*33.433*t9i2

      rr    = den * rev * term 
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_ne18ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee,zz

      double precision z1,a1,ztot,ared,r,c1,c2,c3,c4
      parameter        (z1   = 10.0d0,
     1                  a1   = 18.0d0,
     2                  ztot = 2.0d0 * z1,
     3                  ared = 4.0d0*a1/(4.0d0 + a1),
     4                  r    = 5.1566081196876965d0,
     5                  c1   = 4.9080044545315392d10,
     6                  c2   = 4.9592784569936502d-2,
     7                  c3   = 1.9288564401521285d1,
     8                  c4   = 4.6477847042196437d1)

c..note:
c      r    = 1.09 * a1**oneth + 2.3
c      c1   = 7.833e9 * 0.31 * ztot**fourth/(ared**fivsix)
c      c2   = 0.08617 * 0.1215 * sqrt(ared*r**3/ztot)
c      c3   = 2.0d0 * 0.52495 * sqrt(ared*r*ztot)
c      c4   = 4.2487 * (ztot**2*ared)**oneth


c..ne18ap(a,p)na21
c..was a call to aprate

      aa  = 1.0d0 + c2*t9
      zz  = c2/aa

      bb  = aa**fivsix
      dbb = fivsix*bb*zz

      cc  = t923 * bb
      dcc = twoth*cc*t9i + t923 * dbb 

      dd = aa**oneth
      ddd = oneth*dd*zz

      ee  = t9i13 * dd
      dee = -oneth*ee*t9i + t9i13 * ddd

      zz      = 1.0d0/cc 
      term    = c1*zz * exp(c3 - c4*ee)
      dtermdt = -term*(zz*dcc + c4*dee) 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 0.0d0
      drevdt = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_ne19pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg


c..ne19(p,g)na20

      aa  = 1.71d+6 * t9i23 * exp(-19.431d0*t9i13)
      daa = aa*(-twoth*t9i + oneth*19.431*t9i43)

      bb  = 1.0d0 + 0.021*t913 + 0.130*t923 + 1.95d-2*t9
     1      + 3.86d-2*t943 + 1.47d-02*t953 
      dbb = oneth*0.021*t9i23 + twoth*0.130*t9i13 + 1.95d-2
     1      + fourth*3.86d-2*t913 + fiveth*1.47d-2*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb  


      dd  = 1.89d+5 * t9i23 * exp(-19.431d0*t9i13 - t92/1.304164)
      ddd = dd*(-twoth*t9i + oneth*19.431*t9i43 - 2.0d0*t9/1.304164)

      ee  = 1.0d0 + 0.021*t913 + 2.13*t923 + 0.320*t9 
     1      + 2.80*t943 + 1.07*t953
      dee = oneth*0.021*t9i23 + twoth*2.13*t9i13 + 0.320
     1      + fourth*2.80*t913 + fiveth*1.07*t923

      ff  = dd*ee
      dff = ddd*ee + dd*dee  

      gg  = 8.45d+3 * t9i54 * exp(-7.64d0*t9i)
      dgg = gg*(-fivfour*t9i + 7.64d0*t9i2)


      term    = cc + ff + gg
      dtermdt = dcc + dff + dgg


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.39e+09 * t932 * exp(-25.519*t9i)
      drevdt   = rev*(1.5d0*t9i + 25.519*t9i2)

      rr    = rev * term 
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_si26ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,
     1                 cc,dcc,dd,ddd,ee,dee,zz

      double precision z1,a1,ztot,ared,r,c1,c2,c3,c4
      parameter        (z1   = 14.0d0,
     1                  a1   = 26.0d0,
     2                  ztot = 2.0d0 * z1,
     3                  ared = 4.0d0*a1/(4.0d0 + a1),
     4                  r    = 5.5291207145640335d0, 
     5                  c1   = 7.3266779970543091d10,
     6                  c2   = 4.7895369289991982d-02,
     7                  c3   = 2.4322657793918662d1,
     8                  c4   = 5.9292366232997814d1) 

c..note:
c      r    = 1.09 * a1**oneth + 2.3
c      c1   = 7.833e9 * 0.31 * ztot**fourth/(ared**fivsix)
c      c2   = 0.08617 * 0.1215 * sqrt(ared*r**3/ztot)
c      c3   = 2.0d0 * 0.52495 * sqrt(ared*r*ztot)
c      c4   = 4.2487 * (ztot**2*ared)**oneth



c..si26ap(a,p)p29
c..was a call to aprate

      aa  = 1.0d0 + c2*t9
      zz  = c2/aa

      bb  = aa**fivsix
      dbb = fivsix*bb*zz

      cc  = t923 * bb
      dcc = twoth*cc*t9i + t923 * dbb 

      dd = aa**oneth
      ddd = oneth*dd*zz

      ee  = t9i13 * dd
      dee = -oneth*ee*t9i + t9i13 * ddd

      zz      = 1.0d0/cc
      term    = c1*zz * exp(c3 - c4*ee)
      dtermdt = -term*(zz*dcc + c4*dee) 

c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.0d0
      drevdt   = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_c12ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,f1,df1,f2,df2,
     2                 zz


c..c12(a,g)o16
      aa   = 1.0d0 + 0.0489d0*t9i23
      daa  = -twoth*0.0489d0*t9i53

      bb   = t92*aa*aa
      dbb  = 2.0d0*(bb*t9i + t92*aa*daa)

      cc   = exp(-32.120d0*t9i13 - t92/12.222016d0)
      dcc  = cc * (oneth*32.120d0*t9i43 - 2.0d0*t9/12.222016d0)

      dd   = 1.0d0 + 0.2654d0*t9i23
      ddd  = -twoth*0.2654d0*t9i53

      ee   = t92*dd*dd
      dee  = 2.0d0*(ee*t9i + t92*dd*ddd)

      ff   = exp(-32.120d0*t9i13)
      dff  = ff * oneth*32.120d0*t9i43

      gg   = 1.25d3 * t9i32 * exp(-27.499*t9i)
      dgg  = gg*(-1.5d0*t9i + 27.499*t9i2)

      hh   = 1.43d-2 * t95 * exp(-15.541*t9i)
      dhh  = hh*(5.0d0*t9i + 15.541*t9i2)

      zz   = 1.0d0/bb
      f1   = cc*zz
      df1  = (dcc - f1*dbb)*zz

      zz   = 1.0d0/ee
      f2   = ff*zz
      df2  = (dff - f2*dee)*zz

      term    = 1.04d8*f1  + 1.76d8*f2 + gg + hh
      dtermdt = 1.04d8*df1 + 1.76d8*df2 + dgg + dhh


c..1.7 times cf88 value
      term     = 1.7d0 * term
      dtermdt  = 1.7d0 * dtermdt

      fr    = term * den
      dfrdt = dtermdt * den * 1.0d-9
      dfrdd = term

      rev    = 5.13d10 * t932 * exp(-83.111*t9i)
      drevdt = rev*(1.5d0*t9i + 83.111*t9i2)

      rr     = rev * term
      drrdt  = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd  = 0.0d0
 
      return
      end






      subroutine rate_tripalf(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,r2abe,dr2abedt,rbeac,
     1                 drbeacdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee,
     2                 ff,dff,xx,dxx,yy,dyy,zz,dzz,uu,vv,f1,df1,rc28
      parameter        (rc28   = 0.1d0)



c..triple alfa to c12
c..this is a(a,g)be8
      aa    = 7.40d+05 * t9i32 * exp(-1.0663*t9i) 
      daa   = aa*(-1.5d0*t9i  + 1.0663*t9i2)

      bb    = 4.164d+09 * t9i23 * exp(-13.49*t9i13 - t92/0.009604)
      dbb   = bb*(-twoth*t9i + oneth*13.49*t9i43 - 2.0d0*t9/0.009604)

      cc    = 1.0d0 + 0.031*t913 + 8.009*t923 + 1.732*t9  
     1        + 49.883*t943 + 27.426*t953
      dcc   = oneth*0.031*t9i23 + twoth*8.009*t9i13 + 1.732
     1        + fourth*49.883*t913 + fiveth*27.426*t923

      r2abe    = aa + bb * cc
      dr2abedt = daa + dbb*cc + bb*dcc 


c..this is be8(a,g)c12
      dd    = 130.0d0 * t9i32 * exp(-3.3364*t9i)  
      ddd   = dd*(-1.5d0*t9i + 3.3364*t9i2)

      ee    = 2.510d+07 * t9i23 * exp(-23.57*t9i13 - t92/0.055225) 
      dee   = ee*(-twoth*t9i + oneth*23.57*t9i43 - 2.0d0*t9/0.055225)

      ff    = 1.0d0 + 0.018*t913 + 5.249*t923 + 0.650*t9 + 
     1        19.176*t943 + 6.034*t953
      dff   = oneth*0.018*t9i23 + twoth*5.249*t9i13 + 0.650
     1        + fourth*19.176*t913 + fiveth*6.034*t923

      rbeac    = dd + ee * ff
      drbeacdt = ddd + dee * ff + ee * dff


c..a factor
      xx    = rc28 * 1.35d-07 * t9i32 * exp(-24.811*t9i)
      dxx   = xx*(-1.5d0*t9i + 24.811*t9i2)


c..high temperature rate
      if (t9.gt.0.08) then
       term    = 2.90d-16 * r2abe * rbeac + xx
       dtermdt =   2.90d-16 * dr2abedt * rbeac  
     1           + 2.90d-16 * r2abe * drbeacdt  
     2           + dxx

c..low temperature rate
      else
       uu   = 0.8d0*exp(-(0.025*t9i)**3.263) 
       yy   = 0.01 + 0.2d0 + uu
       dyy  = uu * 3.263*(0.025*t9i)**2.263 * (0.025*t9i2)
       vv   = 4.0d0*exp(-(t9/0.025)**9.227) 
       zz   = 1.0d0 + vv
       dzz  = vv * 9.227*(t9/0.025)**8.227 * 40.0d0
       aa   = 1.0d0/zz
       f1   = yy * aa
       df1  = (dyy - f1*dzz)*aa              
       term = 2.90d-16 * r2abe * rbeac * f1 +  xx 
       dtermdt =   2.90d-16 * dr2abedt * rbeac * f1 
     1           + 2.90d-16 * r2abe * drbeacdt * f1 
     2           + 2.90d-16 * r2abe * rbeac * df1 
     3           + dxx
      end if


c..rates
      fr    = term * den * den
      dfrdt = dtermdt * den * den * 1.0d-9
      dfrdd = 2.0d0 * term * den

      rev    = 2.00d+20*t93*exp(-84.424*t9i)
      drevdt = rev*(3.0d0*t9i + 84.424*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_c12c12(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,
     1                 aa,zz


c..c12 + c12 reaction
      aa      = 1.0d0 + 0.0396*t9
      zz      = 1.0d0/aa

      t9a     = t9*zz
      dt9a    = (1.0d0 -  t9a*0.0396)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56   = t9a**fivsix
      dt9a56  = fivsix*t9a56*zz

      term    = 4.27d+26 * t9a56 * t9i32 * 
     1          exp(-84.165/t9a13 - 2.12d-03*t93)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i
     1                + 84.165/t9a13**2*dt9a13 - 6.36d-3*t92)

c..rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0
 
      return
      end






      subroutine rate_c12c12npa(temp,den,
     1                fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd,
     2                fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd,
     3                fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,
     1                fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd,
     2                fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd,
     3                fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd

c..locals
      double precision term,dtermdt,rev,drevdt,t9a,dt9a,t9a13,dt9a13,
     1                 t9a56,dt9a56,aa,daa,bb,dbb,cc,dcc,dd,ddd,zz,
     2                 b24n,db24n,b24p,db24p,b24a,db24a


c..c12(c12,n)mg23 
c..c12(c12,p)na23 
c..c12(c12,a)ne20 


      aa      = 1.0d0 + 0.0396*t9
      zz      = 1.0d0/aa

      t9a     = t9 * zz
      dt9a    = (1.0d0 -  t9a*0.0396) * zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56   = t9a**fivsix
      dt9a56  = fivsix*t9a56*zz

      aa = 4.27d+26 * t9a56 * t9i32 * 
     1     exp(-84.165/t9a13 - 2.12d-03*t93)

      daa = aa * (dt9a56/t9a56 - 1.5d0*t9i
     1          + 84.165/t9a13**2*dt9a13 - 6.36d-3*t92)


c..neutron branching from dayras switkowski and woosley 1976
      if (t9 .ge. 1.5) then

       bb    =  0.055 * exp(0.976 - 0.789*t9)
       dbb   = -bb*0.789

       b24n  = 0.055  - bb
       db24n = -dbb

      else 

       bb    = 1.0d0 + 0.0789*t9 + 7.74*t92
       dbb   = 0.0789 + 2.0d0*7.74*t9

       cc    = 0.766*t9i3
       dcc   = -3.0d0*cc*t9i

       dd    = bb * cc
       ddd   = dbb*cc + bb*dcc 

       b24n  = 0.859*exp(-dd)
       db24n = -b24n*ddd
      end if


c..proton branching ratio
      if (t9.gt.3.) then
        b24p  = oneth*(1.0d0 - b24n)
        db24p = -oneth*db24n

        b24a  = 2.0d0 * b24p
        db24a = 2.0d0 * db24p

       else
        b24p  = 0.5d0*(1.0d0 - b24n)
        db24p = -0.5d0*db24n

        b24a  = b24p
        db24a = db24p

       end if


c..rates 

c..c12(c12,n)mg23
      term    = aa * b24n
      dtermdt = daa*b24n + aa*db24n
      fr1     = den * term
      dfr1dt  = den * dtermdt * 1.0d-9
      dfr1dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 3.93 * exp(30.16100515d0*t9i)
       drevdt = -rev*30.16100515d0*t9i2
      end if
      rr1    = den * rev * term
      drr1dt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drr1dd = rev*term


c..c12(c12,p)na23
      term    = aa * b24p
      dtermdt = daa*b24p + aa*db24p
      fr2     = den * term
      dfr2dt  = den * dtermdt * 1.0d-9
      dfr2dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 3.93 * exp(-25.98325915d0*t9i)
       drevdt = rev*25.98325915d0*t9i2
      end if
      rr2    = den * rev * term
      drr2dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr2dd = rev*term


c..c12(c12,a)ne20
      term    = aa * b24a
      dtermdt = daa*b24a + aa*db24a
      fr3     = den * term
      dfr3dt  = den * dtermdt * 1.0d-9
      dfr3dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 2.42 * exp(-53.576110995d0*t9i)
       drevdt = rev*53.576110995d0*t9i2
      end if
      rr3    = den * rev * term
      drr3dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr3dd = rev*term

      return
      end






      subroutine rate_c12o16(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,t9a,dt9a,t9a13,dt9a13,t9a23,dt9a23,
     1                 t9a56,dt9a56,aa,daa,bb,dbb,cc,dcc,zz


c..c12 + o16 reaction; see cf88 references 47-4

      if (t9.ge.0.5) then
       aa     = 1.0d0 + 0.055*t9
       zz     = 1.0d0/aa

       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*0.055)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a23  = t9a13*t9a13
       dt9a23 = 2.0d0 * t9a13 * dt9a13

       t9a56  = t9a**fivsix
       dt9a56 = fivsix*t9a56*zz

       aa      = exp(-0.18*t9a*t9a) 
       daa     = -aa * 0.36 * t9a * dt9a

       bb      = 1.06d-03*exp(2.562*t9a23)
       dbb     = bb * 2.562 * dt9a23

       cc      = aa + bb
       dcc     = daa + dbb

       zz      = 1.0d0/cc
       term    = 1.72d+31 * t9a56 * t9i32 * exp(-106.594/t9a13) * zz
       dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i
     1                 + 106.594/t9a23*dt9a13 - zz*dcc)

      else
       term    = 2.6288035d-29
       dtermdt = 0.0d0
      endif


c..the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0      
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_c12o16npa(temp,den,
     1                fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd,
     2                fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd,
     3                fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,
     1                fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd,
     2                fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd,
     3                fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd

c..locals
      double precision term,dtermdt,rev,drevdt,t9a,dt9a,t9a13,dt9a13,
     1                 t9a23,dt9a23,t9a56,dt9a56,aa,daa,bb,dbb,cc,dcc,
     2                 dd,ddd,b27n,b27p,b24a,zz


c..c12(o16,n)si27
c..c12(o16,p)al27
c..c12(o16,a)mg24


      if (t9.ge.0.5) then
       aa     = 1.0d0 + 0.055*t9
       zz     = 1.0d0/aa

       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*0.055)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a23  = t9a13*t9a13
       dt9a23 = 2.0d0 * t9a13 * dt9a13

       t9a56  = t9a**fivsix
       dt9a56 = fivsix*t9a56*zz

       aa     = exp(-0.18*t9a*t9a) 
       daa    = -aa * 0.36 * t9a * dt9a

       bb     = 1.06d-03*exp(2.562*t9a23)
       dbb    = bb * 2.562 * dt9a23

       cc     = aa + bb
       dcc    = daa + dbb

       zz     = 1.0d0/cc
       dd     = 1.72d+31 * t9a56 * t9i32 * exp(-106.594/t9a13) *zz
       ddd    = dd*(dt9a56/t9a56 - 1.5d0*t9i
     1           + 106.594/t9a23 * dt9a13 - dcc*zz)

      else
c       dd     = 2.6288035d-29
       dd     = 0.0d0
       ddd    = 0.0d0
      endif


c..branching ratios from pwnsz data
        b27n = 0.1d0
        b27p = 0.5d0
        b24a = 0.4d0


c..rates 

c..c12(o16,n)si27
      term    = dd * b27n
      dtermdt = ddd * b27n 
      fr1     = den * term
      dfr1dt  = den * dtermdt * 1.0d-9
      dfr1dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 1.58d0 * exp(4.8972467d0*t9i)
       drevdt = -rev*4.8972467d0*t9i2
      end if
      rr1    = den * rev * term
      drr1dt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drr1dd = rev*term


c..c12(o16,p)al27
      term    = dd * b27p
      dtermdt = ddd * b27p 
      fr2     = den * term
      dfr2dt  = den * dtermdt * 1.0d-9
      dfr2dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 1.58d0 * exp(-59.9970745d0*t9i)
       drevdt = rev*59.9970745d0*t9i2
      end if
      rr2    = den * rev * term
      drr2dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr2dd = rev*term


c..c12(o16,a)mg24
      term    = dd * b24a
      dtermdt = ddd * b24a 
      fr3     = den * term
      dfr3dt  = den * dtermdt * 1.0d-9
      dfr3dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 2.83d0 * exp(-78.5648345d0*t9i)
       drevdt = rev*78.5648345d0*t9i2
      end if
      rr3    = den * rev * term
      drr3dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr3dd = rev*term


      return
      end








      subroutine rate_o16o16(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt


c..o16 + o16
      term  = 7.10d36 * t9i23 *
     1        exp(-135.93 * t9i13 - 0.629*t923 
     2             - 0.445*t943 + 0.0103*t9*t9)

      dtermdt = -twoth*term*t9i
     1          + term * (oneth*135.93*t9i43 - twoth*0.629*t9i13
     2                    - fourth*0.445*t913 + 0.0206*t9)


c..rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0      
      drrdt = 0.0d0
      drrdd = 0.0d0
 
      return
      end






      subroutine rate_o16o16npad(temp,den,
     1                fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd,
     2                fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd,
     3                fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd,
     4                fr4,dfr4dt,dfr4dd,rr4,drr4dt,drr4dd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,
     1                fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd,
     2                fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd,
     3                fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd,
     4                fr4,dfr4dt,dfr4dd,rr4,drr4dt,drr4dd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,
     1                 b32n,db32n,b32p,db32p,b32a,db32a,b32d,db32d,
     2                 ezro,dezro,dlt,ddlt,xxt,dxxt,thrs,dthrs


c..o16(o16,n)s31
c..o16(o16,p)p31
c..o16(o16,a)si28
c..o16(o16,d)p30

      aa  = 7.10d36 * t9i23 *
     1        exp(-135.93 * t9i13 - 0.629*t923 
     2             - 0.445*t943 + 0.0103*t9*t9)

      daa = -twoth*aa*t9i
     1       + aa * (oneth*135.93*t9i43 - twoth*0.629*t9i13
     2                    - fourth*0.445*t913 + 0.0206*t9)


c..branching ratios highly uncertain;  guessed using fcz 1975
c..deuteron channel is endoergic. apply error function cut-off.
       ezro = 3.9*t923
       dezro = twoth*ezro*t9i

       dlt  = 1.34*t9**fivsix
       ddlt = fivsix*dlt*t9i

       xxt  = 2.0d0*(2.406 - ezro)/dlt
       dxxt = -(2.0d0*dezro + xxt*ddlt)/dlt

       call fowthrsh(xxt,thrs,dthrs)

       b32d  = 0.05d0*thrs
       db32d = 0.05d0*dthrs*dxxt

       b32n  = 0.1d0
       db32n = 0.0d0

       b32a  = 0.25d0
       db32a = 0.0d0

       b32p  = 1.0d0 - b32d - b32a - b32n
       db32p = -db32d



c..rates 

c..o16(o16,n)s31
      term    = aa * b32n
      dtermdt = daa*b32n + aa*db32n
      fr1     = den * term
      dfr1dt  = den * dtermdt * 1.0d-9
      dfr1dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 5.92 * exp(-16.8038228d0*t9i)
       drevdt = rev*16.8038228d0*t9i2
      end if
      rr1    = den * rev * term
      drr1dt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drr1dd = rev*term

c..o16(o16,p)p31
      term    = aa * b32p
      dtermdt = daa*b32p + aa*db32p
      fr2     = den * term
      dfr2dt  = den * dtermdt * 1.0d-9
      dfr2dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 5.92*exp(-89.0788286d0*t9i)
       drevdt = rev*89.0788286d0*t9i2
      end if
      rr2    = den * rev * term
      drr2dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr2dd = rev*term


c..o16(o16,a)si28
      term    = aa * b32a
      dtermdt = daa*b32a + aa*db32a
      fr3     = den * term
      dfr3dt  = den * dtermdt * 1.0d-9
      dfr3dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 3.46*exp(-111.3137212d0*t9i)
       drevdt = rev*111.3137212d0*t9i2
      end if
      rr3    = den * rev * term
      drr3dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr3dd = rev*term


c..o16(o16,d)p30
      term    = aa * b32d
      dtermdt = daa*b32d + aa*db32d
      fr4     = den * term
      dfr4dt  = den * dtermdt * 1.0d-9
      dfr4dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 0.984*exp(27.9908982d0*t9i)
       drevdt = -rev*27.9908982d0*t9i2
      end if
      rr4    = den * rev * term
      drr4dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr4dd = rev*term

      return
      end






      subroutine fowthrsh(x,thrs,dthrs)
      include 'implno.dek'

c..fowler threshold fudge function. 
c..err func rational (abramowitz p.299)7.1.25 and its derivative

c..declare
      double precision x,thrs,dthrs,
     1                 ag,z,z2,t,t2,t3,tt,er,aa,daa,dt,dtt,der

      ag   = sign(1.0d0,x)
      z    = abs(x)
      z2   = z*z
      aa   = 1.0d0 + 0.47047d0*z

      t    = 1.0d0/aa
      dt   = -t*t*0.47047*ag

      t2   = t*t
      t3   = t2*t

      tt   = 0.3480242d0*t - 0.0958798d0*t2 + 0.7478556d0*t3
      dtt  = dt * (0.3480242d0 - 2.0d0*0.0958798d0*t 
     1             + 3.0d0*0.7478556d0*t2)

      thrs  = 0.5d0
      dthrs = -5.6452433937900004d-1

      if (z .ne. 0) then
       aa   = exp(-z2)
       daa  = -2.0d0*aa*z*ag

       er   = 1.0d0 - tt * aa
       der  = -dtt*aa - tt*daa

       thrs  = 0.5d0 * (1.0d0 - ag*er)
       dthrs = -0.5d0*ag*der
      end if
      return
      end





      subroutine rate_o16ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,term1,dterm1,aa,daa,bb,dbb,
     1                 cc,dcc,term2,dterm2,rev,drevdt



c..o16(a,g)ne20
      term1   = 9.37d9 * t9i23 * exp(-39.757*t9i13 - t92/2.515396)
      dterm1  = term1*(-twoth*t9i 
     1          + oneth*39.757*t9i43 - 2.0d0*t9/2.515396)

      aa      = 62.1 * t9i32 * exp(-10.297*t9i)  
      daa     = aa*(-1.5d0*t9i + 10.297*t9i2)

      bb      = 538.0d0 * t9i32 * exp(-12.226*t9i)  
      dbb     = bb*(-1.5d0*t9i + 12.226*t9i2)

      cc      = 13.0d0 * t92 * exp(-20.093*t9i)
      dcc     = cc*(2.0d0*t9i + 20.093*t9i2)

      term2   = aa + bb + cc 
      dterm2  = daa + dbb + dcc

      term    = term1 + term2
      dtermdt = dterm1 + dterm2  


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.65d+10*t932*exp(-54.937*t9i)
      drevdt   = rev*(1.5d0*t9i + 54.937*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0
 
      return
      end




      subroutine rate_ne20ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,term1,dterm1,aa,daa,bb,dbb,
     1                 term2,dterm2,term3,dterm3,rev,drevdt,zz,rc102
      parameter        (rc102 = 0.1d0)



c..ne20(a,g)mg24

      aa   = 4.11d+11 * t9i23 * exp(-46.766*t9i13 - t92/4.923961) 
      daa  = aa*(-twoth*t9i + oneth*46.766*t9i43 - 2.0d0*t9/4.923961)

      bb   = 1.0d0 + 0.009*t913 + 0.882*t923 + 0.055*t9  
     1       + 0.749*t943 + 0.119*t953
      dbb  = oneth*0.009*t9i23 + twoth*0.882*t9i13 + 0.055
     1       + fourth*0.749*t913 + fiveth*0.119*t923

      term1  = aa * bb
      dterm1 = daa * bb + aa * dbb


      aa   = 5.27d+03 * t9i32 * exp(-15.869*t9i)  
      daa  = aa*(-1.5d0*t9i + 15.869*t9i2)

      bb   = 6.51d+03 * t912 * exp(-16.223*t9i)
      dbb  = bb*(0.5d0*t9i + 16.223*t9i2)
 
      term2  = aa + bb
      dterm2 = daa + dbb   


      aa   = 42.1 * t9i32 * exp(-9.115*t9i) 
      daa  = aa*(-1.5d0*t9i + 9.115*t9i2)

      bb   =  32.0 * t9i23 * exp(-9.383*t9i)
      dbb  = bb*(-twoth*t9i + 9.383*t9i2)

      term3  = rc102 * (aa + bb)
      dterm3 = rc102 * (daa + dbb)   


      aa  = 5.0d0*exp(-18.960*t9i)
      daa = aa*18.960*t9i2

      bb  = 1.0d0 + aa
      dbb = daa 

      zz      = 1.0d0/bb
      term    = (term1 + term2 + term3)*zz
      dtermdt = ((dterm1 + dterm2 + dterm3) - term*dbb)*zz


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.01d+10 * t932 * exp(-108.059*t9i)
      drevdt   = rev*(1.5d0*t9i + 108.059*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg24ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee,
     1                 ff,dff,gg,dgg,hh,hhi,rev,drevdt,rc121
      parameter        (rc121 = 0.1d0)



c..24mg(a,g)28si 

      aa    = 4.78d+01 * t9i32 * exp(-13.506*t9i) 
      daa   = aa*(-1.5d0*t9i + 13.506*t9i2)

      bb    =  2.38d+03 * t9i32 * exp(-15.218*t9i)
      dbb   = bb*(-1.5d0*t9i + 15.218*t9i2)  

      cc    = 2.47d+02 * t932 * exp(-15.147*t9i) 
      dcc   = cc*(1.5d0*t9i + 15.147*t9i2) 

      dd    = rc121 * 1.72d-09 * t9i32 * exp(-5.028*t9i)
      ddd   = dd*(-1.5d0*t9i + 5.028*t9i2)

      ee    = rc121* 1.25d-03 * t9i32 * exp(-7.929*t9i)
      dee   = ee*(-1.5d0*t9i + 7.929*t9i2)

      ff    = rc121 * 2.43d+01 * t9i * exp(-11.523*t9i)
      dff   = ff*(-t9i + 11.523*t9i2)

      gg    = 5.0d0*exp(-15.882*t9i)
      dgg   = gg*15.882*t9i2
   
      hh    = 1.0d0 + gg
      hhi   = 1.0d0/hh

      term    = (aa + bb + cc + dd + ee + ff) * hhi
      dtermdt = (daa + dbb + dcc + ddd + dee + dff - term*dgg) * hhi


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.27d+10 * t932 * exp(-115.862*t9i)
      drevdt   = rev*(1.5d0*t9i + 115.862*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_mg24ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee,
     1                 ff,dff,gg,dgg,term1,dterm1,term2,dterm2,
     2                 rev,drevdt,rc148
      parameter        (rc148 = 0.1d0)



c..24mg(a,p)al27
      aa     = 1.10d+08 * t9i23 * exp(-23.261*t9i13 - t92/0.024649) 
      daa    = -twoth*aa*t9i + aa*(23.261*t9i43 - 2.0d0*t9/0.024649)

      bb     =  1.0d0 + 0.018*t913 + 12.85*t923 + 1.61*t9  
     1         + 89.87*t943 + 28.66*t953
      dbb    = oneth*0.018*t9i23 + twoth*12.85*t9i13 + 1.61
     1          + fourth*89.87*t913 + fiveth*28.66*t923

      term1  = aa * bb
      dterm1 = daa * bb + aa * dbb  

      aa     = 129.0d0 * t9i32 * exp(-2.517*t9i) 
      daa    = -1.5d0*aa*t9i + aa*2.517*t9i2

      bb     = 5660.0d0 * t972 * exp(-3.421*t9i) 
      dbb    = 3.5d0*bb*t9i +  bb*3.421*t9i2

      cc     = rc148 * 3.89d-08 * t9i32 * exp(-0.853*t9i)  
      dcc    = -1.5d0*cc*t9i + cc*0.853*t9i2

      dd     = rc148 * 8.18d-09 * t9i32 * exp(-1.001*t9i)
      ddd    = -1.5d0*dd*t9i + dd*1.001*t9i2

      term2  = aa + bb + cc + dd
      dterm2 = daa + dbb + dcc + ddd  

      ee     = oneth*exp(-9.792*t9i)
      dee    = ee*9.792*t9i2

      ff     =  twoth * exp(-11.773*t9i)
      dff    = ff*11.773*t9i2

      gg     = 1.0d0 + ee + ff
      dgg    = dee + dff 

      term    = (term1 + term2)/gg
      dtermdt = ((dterm1 + dterm2) - term*dgg)/gg


c..the rates 
      rev      = 1.81 * exp(-18.572*t9i)
      drevdt   = rev*18.572*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end






      subroutine rate_al27pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,
     1                 dd,ddd,ee,dee,ff,dff,gg,dgg


c..al27(p,g)si28
c..champagne 1996 

      aa  = 1.32d+09 * t9i23 * exp(-23.26*t9i13)
      daa = aa*(-twoth*t9i + oneth*23.26*t9i43)

      bb  = 3.22d-10 * t9i32 * exp(-0.836*t9i)*0.17
      dbb = bb*(-1.5d0*t9i + 0.836*t9i2)

      cc  = 1.74d+00 * t9i32 * exp(-2.269*t9i)
      dcc = cc*(-1.5d0*t9i + 2.269*t9i2)

      dd  = 9.92d+00 * t9i32 * exp(-2.492*t9i)
      ddd = dd*(-1.5d0*t9i + 2.492*t9i2)

      ee  = 4.29d+01 * t9i32 * exp(-3.273*t9i)
      dee = ee*(-1.5d0*t9i + 3.273*t9i2)

      ff  = 1.34d+02 * t9i32 * exp(-3.654*t9i)
      dff = ff*(-1.5d0*t9i + 3.654*t9i2)

      gg  = 1.77d+04 * (t9**0.53) * exp(-4.588*t9i)
      dgg = gg*(0.53*t9i + 4.588*t9i2)

      term    = aa + bb + cc + dd + ee + ff + gg
      dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg


c..rates
      fr    = den * term 
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.13d+11 * t932 * exp(-134.434*t9i)
      drevdt   = rev*(1.5d0*t9i + 134.434*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_al27pg_old(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee,
     1                 ff,dff,gg,dgg,hh,dhh,xx,dxx,yy,dyy,zz,dzz,pp,
     2                 rev,drevdt,rc147
      parameter        (rc147 = 0.1d0)



c..27al(p,g)si28  cf88

      aa  = 1.67d+08 * t9i23 * exp(-23.261*t9i13 - t92/0.024025)  
      daa = aa*(-twoth*t9i + oneth*23.261*t9i43 - 2.0d0*t9/0.024025) 

      bb  = 1.0d0 + 0.018*t913 + 5.81*t923 + 0.728*t9 
     1      + 27.31*t943 + 8.71*t953  
      dbb = oneth*0.018*t9i23 + twoth*5.81*t9i13 + 0.728
     1      + fourth*27.31*t913 + fiveth*8.71*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb 

      dd  = 2.20d+00 * t9i32 * exp(-2.269*t9i)  
      ddd = dd*(-1.5d0*t9i + 2.269*t9i2)

      ee  = 1.22d+01 * t9i32 * exp(-2.491*t9i) 
      dee = ee*(-1.5d0*t9i + 2.491*t9i2)

      ff  =  1.50d+04 * t9 * exp(-4.112*t9i)  
      dff = ff*(t9i + 4.112*t9i2)

      gg  = rc147 * 6.50d-10 * t9i32 * exp(-0.853*t9i) 
      dgg = gg*(-1.5d0*t9i + 0.853*t9i2)

      hh  = rc147 * 1.63d-10 * t9i32 * exp(-1.001*t9i)
      dhh = hh*(-1.5d0*t9i + 1.001*t9i2)

      xx     = oneth*exp(-9.792*t9i)
      dxx    = xx*9.792*t9i2

      yy     =  twoth * exp(-11.773*t9i)
      dyy    = yy*11.773*t9i2

      zz     = 1.0d0 + xx + yy
      dzz    = dxx + dyy 

      pp      = 1.0d0/zz
      term    = (cc + dd + ee + ff + gg + hh)*pp
      dtermdt = ((dcc + ddd + dee + dff + dgg + dhh) - term*dzz)*pp


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.13d+11*t932*exp(-134.434*t9i)
      drevdt   = rev*(1.5d0*t9i + 134.434*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_si28ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..si28(a,g)s32
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 6.340d-2*z + 2.541d-3*z2 - 2.900d-4*z3
      if (z .eq. 10.0) then
       daa = 0
      else
       daa   = 6.340d-2 + 2.0d0*2.541d-3*t9 - 3.0d0*2.900d-4*t92
      end if

      term    = 4.82d+22 * t9i23 * exp(-61.015 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 61.015*t9i13*(oneth*t9i*aa - daa)) 

c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.461d+10 * t932 * exp(-80.643*t9i)
      drevdt   = rev*(1.5d0*t9i + 80.643*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si28ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..si28(a,p)p31
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 2.798d-3*z + 2.763d-3*z2 - 2.341d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 2.798d-3 + 2.0d0*2.763d-3*t9 - 3.0d0*2.341d-4*t92
      end if

      term    = 4.16d+13 * t9i23 * exp(-25.631 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*25.631*t9i13*(oneth*t9i*aa - daa) 


c..the rates 
      rev      = 0.5825d0 * exp(-22.224*t9i)
      drevdt   = rev*22.224*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_p31pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..p31(p,g)s32
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.928d-1*z - 1.540d-2*z2 + 6.444d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.928d-1 - 2.0d0*1.540d-2*t9 + 3.0d0*6.444d-4*t92
      end if

      term    = 1.08d+16 * t9i23 * exp(-27.042 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 27.042*t9i13*(oneth*t9i*aa - daa)) 


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.764d+10 * t932 * exp(-102.865*t9i)
      drevdt   = rev*(1.5d0*t9i + 102.865*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_s32ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..s32(a,g)ar36
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 4.913d-2*z + 4.637d-3*z2 - 4.067d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 4.913d-2 + 2.0d0*4.637d-3*t9 - 3.0d0*4.067d-4*t92
      end if

      term    = 1.16d+24 * t9i23 * exp(-66.690 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 66.690*t9i13*(oneth*t9i*aa - daa)) 


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.616d+10 * t932 * exp(-77.080*t9i)
      drevdt   = rev*(1.5d0*t9i + 77.080*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_s32ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..s32(a,p)cl35
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.041d-1*z - 1.368d-2*z2 + 6.969d-4*z3
      if (z .eq. 10) then
       daa = 0.0d0
      else
       daa   = 1.041d-1 - 2.0d0*1.368d-2*t9 + 3.0d0*6.969d-4*t92
      end if

      term    = 1.27d+16 * t9i23 * exp(-31.044 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*31.044*t9i13*(oneth*t9i*aa - daa) 


c..the rates 
      rev      = 1.144 * exp(-21.643*t9i)
      drevdt   = rev*21.643*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_cl35pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt


c..cl35(p,g)ar36
      aa    = 1.0d0 + 1.761d-1*t9 - 1.322d-2*t92 + 5.245d-4*t93
      daa   = 1.761d-1 - 2.0d0*1.322d-2*t9 + 3.0d0*5.245d-4*t92


      term    =  4.48d+16 * t9i23 * exp(-29.483 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 29.483*t9i13*(oneth*t9i*aa - daa)) 


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.568d+10*t932*exp(-98.722*t9i)
      drevdt   = rev*(1.5d0*t9i + 98.722*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_ar36ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..ar36(a,g)ca40
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.458d-1*z - 1.069d-2*z2 + 3.790d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.458d-1 - 2.0d0*1.069d-2*t9 + 3.0d0*3.790d-4*t92
      end if

      term    = 2.81d+30 * t9i23 * exp(-78.271 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 78.271*t9i13*(oneth*t9i*aa - daa)) 


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.740d+10 * t932 * exp(-81.711*t9i)
      drevdt   = rev*(1.5d0*t9i + 81.711*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ar36ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..ar36(a,p)k39
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 4.826d-3*z - 5.534d-3*z2 + 4.021d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 4.826d-3 - 2.0d0*5.534d-3*t9 + 3.0d0*4.021d-4*t92
      end if

      term    = 2.76d+13 * t9i23 * exp(-34.922 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*34.922*t9i13*(oneth*t9i*aa - daa) 


c..the rates 
      rev      = 1.128*exp(-14.959*t9i)
      drevdt   = rev*14.959*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_k39pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..k39(p,g)ca40
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.622d-1*z - 1.119d-2*z2 + 3.910d-4*z3
      if (z .eq. 10) then
       daa = 0.0d0
      else 
       daa   = 1.622d-1 - 2.0d0*1.119d-2*t9 + 3.0d0*3.910d-4*t92
      end if

      term    = 4.09d+16 * t9i23 * exp(-31.727 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 31.727*t9i13*(oneth*t9i*aa - daa)) 


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.600d+10 * t932 * exp(-96.657*t9i)
      drevdt   = rev*(1.5d0*t9i + 96.657*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ca40ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..ca40(a,g)ti44
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.650d-2*z + 5.973d-3*z2 - 3.889d-04*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.650d-2 + 2.0d0*5.973d-3*t9 - 3.0d0*3.889d-4*t92
      end if

      term    = 4.66d+24 * t9i23 * exp(-76.435 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 76.435*t9i13*(oneth*t9i*aa - daa)) 


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.843d+10 * t932 * exp(-59.510*t9i)
      drevdt   = rev*(1.5d0*t9i + 59.510*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ca40ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..ca40(a,p)sc43
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 - 1.206d-2*z + 7.753d-3*z2 - 5.071d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = -1.206d-2 + 2.0d0*7.753d-3*t9 - 3.0d0*5.071d-4*t92
      end if 

      term    = 4.54d+14 * t9i23 * exp(-32.177 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*32.177*t9i13*(oneth*t9i*aa - daa) 


c..the rates 
      rev      = 2.229 * exp(-40.966*t9i)
      drevdt   = rev*40.966*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_sc43pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..sc43(p,g)ca40
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.023d-1*z - 2.242d-3*z2 - 5.463d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.023d-1 - 2.0d0*2.242d-3*t9 - 3.0d0*5.463d-5*t92
      end if

      term    = 3.85d+16 * t9i23 * exp(-33.234 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 33.234*t9i13*(oneth*t9i*aa - daa)) 


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.525d+11 * t932 * exp(-100.475*t9i)
      drevdt   = rev*(1.5d0*t9i + 100.475*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_ti44ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..ti44(a,g)cr48
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.066d-1*z - 1.102d-2*z2 + 5.324d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.066d-1 - 2.0d0*1.102d-2*t9 + 3.0d0*5.324d-4*t92
      end if

      term    = 1.37d+26 * t9i23 * exp(-81.227 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 81.227*t9i13*(oneth*t9i*aa - daa)) 


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.928d+10*t932*exp(-89.289*t9i)
      drevdt   = rev*(1.5d0*t9i + 89.289*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ti44ap(temp,den,
     1                      fr,dfrdt,dfrdd,
     2                      rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,
     1                 fr,dfrdt,dfrdd,
     2                 rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..ti44(a,p)v47
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 2.655d-2*z - 3.947d-3*z2 + 2.522d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else 
       daa   = 2.655d-2 - 2.0d0*3.947d-3*t9 + 3.0d0*2.522d-4*t92
      end if

      term    = 6.54d+20 * t9i23 * exp(-66.678 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*66.678*t9i13*(oneth*t9i*aa - daa) 


c..the rates 
      rev      = 1.104 * exp(-4.723*t9i)
      drevdt   = rev*4.723*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_v47pg(temp,den,
     1                      fr,dfrdt,dfrdd,
     2                      rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,
     1                 fr,dfrdt,dfrdd,
     2                 rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..v47(p,g)cr48
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 9.979d-2*z - 2.269d-3*z2 - 6.662d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 9.979d-2 - 2.0d0*2.269d-3*t9 - 3.0d0*6.662d-5*t92
      end if

      term    = 2.05d+17 * t9i23 * exp(-35.568 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 35.568*t9i13*(oneth*t9i*aa - daa)) 


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.649d+10*t932*exp(-93.999*t9i)
      drevdt   = rev*(1.5d0*t9i + 93.999*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_cr48ag(temp,den,
     1                      fr,dfrdt,dfrdd,
     2                      rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,
     1                 fr,dfrdt,dfrdd,
     2                 rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..cr48(a,g)fe52
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 6.325d-2*z - 5.671d-3*z2 + 2.848d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 6.325d-2 - 2.0d0*5.671d-3*t9 + 3.0d0*2.848d-4*t92
      end if

      term    = 1.04d+23 * t9i23 * exp(-81.420 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 81.420*t9i13*(oneth*t9i*aa - daa)) 


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.001d+10 * t932 * exp(-92.177*t9i)
      drevdt   = rev*(1.5d0*t9i + 92.177*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_cr48ap(temp,den,
     1                      fr,dfrdt,dfrdd,
     2                      rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,
     1                 fr,dfrdt,dfrdd,
     2                 rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..cr48(a,p)mn51
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.384d-2*z + 1.081d-3*z2 - 5.933d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else 
       daa   = 1.384d-2 + 2.0d0*1.081d-3*t9 - 3.0d0*5.933d-5*t92
      end if

      term    = 1.83d+26 * t9i23 * exp(-86.741 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*86.741*t9i13*(oneth*t9i*aa - daa) 


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.6087*exp(-6.510*t9i)
      drevdt   = rev*6.510*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_mn51pg(temp,den,
     1                      fr,dfrdt,dfrdd,
     2                      rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,
     1                 fr,dfrdt,dfrdd,
     2                 rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..mn51(p,g)fe52
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 8.922d-2*z - 1.256d-3*z2 - 9.453d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 8.922d-2 - 2.0d0*1.256d-3*t9 - 3.0d0*9.453d-5*t92
      end if

      term    = 3.77d+17 * t9i23 * exp(-37.516 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 37.516*t9i13*(oneth*t9i*aa - daa)) 


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.150d+11*t932*exp(-85.667*t9i)
      drevdt   = rev*(1.5d0*t9i + 85.667*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_fe52ag(temp,den,
     1                      fr,dfrdt,dfrdd,
     2                      rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,
     1                 fr,dfrdt,dfrdd,
     2                 rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..fe52(a,g)ni56
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 7.846d-2*z - 7.430d-3*z2 + 3.723d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 7.846d-2 - 2.0d0*7.430d-3*t9 + 3.0d0*3.723d-4*t92
      end if

      term    = 1.05d+27 * t9i23 * exp(-91.674 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 91.674*t9i13*(oneth*t9i*aa - daa)) 


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.064d+10*t932*exp(-92.850*t9i)
      drevdt   = rev*(1.5d0*t9i + 92.850*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_fe52ap(temp,den,
     1                      fr,dfrdt,dfrdd,
     2                      rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,
     1                 fr,dfrdt,dfrdd,
     2                 rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..fe52(a,p)co55
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.367d-2*z + 7.428d-4*z2 - 3.050d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else 
       daa   = 1.367d-2 + 2.0d0*7.428d-4*t9 - 3.0d0*3.050d-5*t92
      end if

      term    = 1.30d+27 * t9i23 * exp(-91.674 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*91.674*t9i13*(oneth*t9i*aa - daa) 


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.4597*exp(-9.470*t9i)
      drevdt   = rev*9.470*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_co55pg(temp,den,
     1                      fr,dfrdt,dfrdd,
     2                      rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,
     1                 fr,dfrdt,dfrdd,
     2                 rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


c..co55(p,g)ni56
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 9.894d-2*z - 3.131d-3*z2 - 2.160d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 9.894d-2 - 2.0d0*3.131d-3*t9 - 3.0d0*2.160d-5*t92
      end if

      term    = 1.21d+18 * t9i23 * exp(-39.604 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 39.604*t9i13*(oneth*t9i*aa - daa)) 


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.537d+11*t932*exp(-83.382*t9i)
      drevdt   = rev*(1.5d0*t9i + 83.382*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_fe52ng(temp,den,
     1                      fr,dfrdt,dfrdd,
     2                      rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,
     1                 fr,dfrdt,dfrdd,
     2                 rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,tq2


c..fe52(n,g)fe53
      tq2     = t9 - 0.348d0
      term    = 9.604d+05 * exp(-0.0626*tq2)
      dtermdt = -term*0.0626

c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.43d+09 * t932 * exp(-123.951*t9i) 
      drevdt   = rev*(1.5d0*t9i + 123.951*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_fe53ng(temp,den,
     1                      fr,dfrdt,dfrdd,
     2                      rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,
     1                 fr,dfrdt,dfrdd,
     2                 rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,tq1,tq10,dtq10,tq2


c..fe53(n,g)fe54
      tq1   = t9/0.348
      tq10  = tq1**0.10
      dtq10 = 0.1d0*tq10/(0.348*tq1)
      tq2   = t9 - 0.348d0

      term    = 1.817d+06 * tq10 * exp(-0.06319*tq2)
      dtermdt = term/tq10*dtq10 - term*0.06319

c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.56d+11 * t932 * exp(-155.284*t9i)
      drevdt   = rev*(1.5d0*t9i + 155.284*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_fe54pg(temp,den,
     1                      fr,dfrdt,dfrdd,
     2                      rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


c..declare the pass
      double precision temp,den,
     1                 fr,dfrdt,dfrdd,
     2                 rr,drrdt,drrdd

c..locals
      double precision term,dtermdt,rev,drevdt,aa,daa,z,z2,z3


c..fe54(p,g)co55
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 9.593d-2*z - 3.445d-3*z2 + 8.594d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 9.593d-2 - 2.0d0*3.445d-3*t9 + 3.0d0*8.594d-5*t92
      end if

      term    = 4.51d+17 * t9i23 * exp(-38.483 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 38.483*t9i13*(oneth*t9i*aa - daa)) 


c..the rates 
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.400d+09 * t932 * exp(-58.605*t9i)
      drevdt   = rev*(1.5d0*t9i + 58.605*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end









c..function wien1
c..function dwien1dx
c..function wien2
c..function dwien2dx
c..function func1
c..function dfunc1dx
c..function func2
c..function dfunc2dx
c..function bb_qromb
c..function bb_trapzd
c..routine bb_polint does polynomial interpolation





      double precision function wien1(x)
      include 'implno.dek'
      include 'const.dek'
      
c..this is the function given in 
c..weinberg's "gravitation and cosmology" page 537, equation 15.6.40

c..declare the pass
      double precision x

c..communicate xcom via common block
      double precision xcom
      common /tes1/    xcom

c..local variables
      external         func1
      double precision func1,f1,con


c..the integration limits ylo and yhi, along with the integration
c..tolerance tol, will give at least 9 significant figures of precision

      double precision ylo,yhi,tol
      parameter        (ylo = 1.0d-6,
     1                  yhi = 50.0d0,
     2                  tol = 1.0d-10,
     3                  con = 45.0d0/(2.0d0*pi*pi*pi*pi))


c..for quadrature
      integer          nquad,ifirst
      parameter        (nquad = 100)
      double precision xquad(nquad),wquad(nquad)      
      data             ifirst/0/



c..initialization of the quadrature abcissas and weights
      if (ifirst .eq. 0) then
       ifirst = 1
       call bb_gauleg(ylo,yhi,xquad,wquad,nquad)
      end if


c..don't do any integration if x is large enough
      if (x .gt. 50.0) then
       wien1 = 1.0d0

c..do the integration
      else
       xcom = x
c       call bb_qromb(func1,ylo,yhi,tol,f1)  
       call bb_qgaus(func1,xquad,wquad,nquad,f1)
       wien1 = 1.0d0 + con * f1
      end if

      return
      end





      double precision function dwien1dx(x)
      include 'implno.dek'
      include 'const.dek'
      
c..this is the derivative with respect to x of the function given in 
c..weinberg's "gravitation and cosmology" page 537, equation 15.6.40

c..declare the pass
      double precision x

c..communicate xcom via common block
      double precision xcom
      common /tes1/    xcom

c..local variables
      external         dfunc1dx
      double precision dfunc1dx,df1,con


c..the integration limits ylo and yhi, along with the integration
c..tolerance tol, will give at least 9 significant figures of precision

      double precision ylo,yhi,tol
      parameter        (ylo = 1.0d-6,
     1                  yhi = 50.0d0,
     2                  tol = 1.0d-10,
     3                  con = 45.0d0/(2.0d0*pi*pi*pi*pi))


c..for quadrature
      integer          nquad,ifirst
      parameter        (nquad = 100)
      double precision xquad(nquad),wquad(nquad)      
      data             ifirst/0/


c..initialization of the quadrature abcissas and weights
      if (ifirst .eq. 0) then
       ifirst = 1
       call bb_gauleg(ylo,yhi,xquad,wquad,nquad)
      end if


c..don't do any integration if x is large enough
      if (x .gt. 50.0) then
       dwien1dx = 0.0d0

c..do the integration
      else
       xcom = x
c       call bb_qromb(dfunc1dx,ylo,yhi,tol,df1)  
       call bb_qgaus(dfunc1dx,xquad,wquad,nquad,df1)
       dwien1dx = con * df1
      end if

      return
      end







      double precision function wien2(x)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'
      
c..this is the function given in 
c..weinberg's "gravitation and cosmology" page 537, equation 15.6.40

c..declare the pass
      double precision x


c..communicate xcom via common block
      double precision xcom
      common /tes1/    xcom


c..communicate the number of neutrino families
c..using 2 families of neutrinos duplicates the time-temperature
c..table in wienberg's "gravitation and cosmology", page 540, table 15.4

c..brought in through network.dek
c      double precision xnnu
c      common /nufam/   xnnu


c..local variables
      external         func2
      double precision func2,f2,wien1


c..the integration limits ylo and yhi, along with the integration
c..tolerance tol, will give at least 9 significant figures of precision

      double precision ylo,yhi,tol,con1,con2,con3,fthirds
      parameter        (ylo     = 1.0d-6,
     1                  yhi     = 50.0d0,
     2                  tol     = 1.0d-10,
     3                  con2    = 4.0d0/11.0d0,
     4                  con3    = 30.0d0/(pi*pi*pi*pi),
     5                  fthirds = 4.0d0/3.0d0)


c..for quadrature
      integer          nquad,ifirst
      parameter        (nquad = 100)
      double precision xquad(nquad),wquad(nquad)      
      data             ifirst/0/


c..initialization of the quadrature abcissas and weights
      if (ifirst .eq. 0) then
       ifirst = 1
       call bb_gauleg(ylo,yhi,xquad,wquad,nquad)

c..a constant that depends on the number of neutrino families
       con1 = xnnu * 7.0d0/8.0d0
      end if



c..don't do any integration if x is large enough
      if (x .gt. 50.0) then
       wien2 = 1.0d0 + con1*con2**fthirds

c..do the integration
      else 
       xcom = x
c       call bb_qromb(func2,ylo,yhi,tol,f2)  
       call bb_qgaus(func2,xquad,wquad,nquad,f2)
       wien2 = 1.0d0 + con1 * (con2 * wien1(x))**fthirds + con3 * f2
      end if

      return
      end






      double precision function dwien2dx(x)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'
      
c..this is the derivative with respect to x of the function given in 
c..weinberg's "gravitation and cosmology" page 537, equation 15.6.40

c..declare the pass
      double precision x



c..communicate xcom via common block
      double precision xcom
      common /tes1/    xcom


c..communicate the number of neutrino families
c..using 2 families of neutrinos duplicates the time-temperature
c..table in wienberg's "gravitation and cosmology", page 540, table 15.4

c..brought in through network.dek
c      double precision xnnu
c      common /nufam/   xnnu


c..local variables
      external         dfunc2dx
      double precision dfunc2dx,df2,wien1,w1,dwien1dx,dw1


c..the integration limits ylo and yhi, along with the integration
c..tolerance tol, will give at least 9 significant figures of precision

      double precision ylo,yhi,tol,con1,con2,con3,fthirds,third
      parameter        (ylo     = 1.0d-6,
     1                  yhi     = 50.0d0,
     2                  tol     = 1.0d-10,
     3                  con2    = 4.0d0/11.0d0,
     4                  con3    = 30.0d0/(pi*pi*pi*pi),
     5                  fthirds = 4.0d0/3.0d0,
     6                  third   = 1.0d0/3.0d0)


c..for quadrature
      integer          nquad,ifirst
      parameter        (nquad = 100)
      double precision xquad(nquad),wquad(nquad)      
      data             ifirst/0/


c..initialization of the quadrature abcissas and weights
      if (ifirst .eq. 0) then
       ifirst = 1
       call bb_gauleg(ylo,yhi,xquad,wquad,nquad)

c..a constant that depends on the number of neutrino families
       con1 = xnnu * 7.0d0/8.0d0
      end if



c..don't do any integration if x is large enough
      if (x .gt. 50.0) then
       dwien2dx = 0.0d0

c..do the integration
      else 
       xcom = x
c       call bb_qromb(dfunc2dx,ylo,yhi,tol,df2)  
       call bb_qgaus(dfunc2dx,xquad,wquad,nquad,df2)
       w1   = wien1(x)
       dw1  = dwien1dx(x)
c       w2   = 1.0d0 + con1*(con2*w1)**fthirds + con3*f2
       dwien2dx = fthirds*con1*(con2*w1)**third * con2*dw1 + con3*df2
      end if

      return
      end






      double precision function func1(y)  
      include 'implno.dek'
      include 'const.dek'

c..this is the integrand of the function given in 
c..weinberg's "gravitation and cosmology" page 537, equation 15.6.40

c..declare the pass
      double precision y

c..local variables 
      double precision y2,x2,aa,aalim
      parameter        (aalim = 200.0d0)


c..communicate xcom via common block
      double precision xcom
      common /tes1/    xcom

      func1 = 0.0d0
      if (aa .le. aalim) then
       y2    = y * y  
       x2    = xcom * xcom
       aa    = sqrt(y2 + x2)
       func1 = (aa + y2/(3.0d0*aa)) * y2 / (exp(aa) + 1.0d0)
      end if
      return
      end   






      double precision function dfunc1dx(y)  
      include 'implno.dek'
      include 'const.dek'

c..this is the derivative with respect to x of the integrand in the function 
c..given by weinberg's "gravitation and cosmology" page 537, equation 15.6.40

c..declare the pass
      double precision y

c..local variables 
      double precision y2,x2,aa,daa,zz,denom,ddenom,f1,aalim
      parameter        (aalim = 200.0d0)


c..communicate xcom via common block
      double precision xcom
      common /tes1/    xcom

      dfunc1dx = 0.0d0
      if (aa .le. aalim) then
       y2       = y * y  
       x2       = xcom * xcom
       aa       = sqrt(y2 + x2)
       daa      = xcom/aa
       zz       = exp(aa)
       denom    = zz + 1.0d0
       ddenom   = zz * daa
       f1       = (aa + y2/(3.0d0*aa)) * y2 / denom
       dfunc1dx = (1.0d0 - y2/(3.0d0*aa**2)) * daa * y2 / denom
     1            - f1/denom * ddenom
      end if
      return
      end   





      double precision function func2(y)  
      include 'implno.dek'
      include 'const.dek'

c..this is the integrand of the function given in 
c..weinberg's "gravitation and cosmology" page 539, equation 15.6.48

c..declare the pass
      double precision y

c..local variables 
      double precision y2,x2,aa,aalim
      parameter        (aalim = 200.0d0)


c..communicate xcom via common block
      double precision xcom
      common /tes1/    xcom

      func2 = 0.0d0
      if (aa .le. aalim) then
       y2    = y * y  
       x2    = xcom * xcom
       aa    = sqrt(y2 + x2)
       func2 = aa * y2 / (exp(aa) + 1.0d0)
      end if
      return
      end   

      




      double precision function dfunc2dx(y)  
      include 'implno.dek'
      include 'const.dek'

c..this is the derivative with respect to x of the integrand 
c..of the function given in weinberg's "gravitation and cosmology" 
c..page 539, equation 15.6.48

c..declare the pass
      double precision y

c..local variables 
      double precision y2,x2,aa,daa,zz,denom,ddenom,aalim,f2
      parameter        (aalim = 200.0d0)


c..communicate xcom via common block
      double precision xcom
      common /tes1/    xcom

      dfunc2dx = 0.0d0
      if (aa .le. aalim) then
       y2       = y * y  
       x2       = xcom * xcom
       aa       = sqrt(y2 + x2)
       daa      = xcom/aa
       zz       = exp(aa)
       denom    = zz + 1.0d0
       ddenom   = zz * daa
       f2       = aa * y2 / denom
       dfunc2dx = (daa*y2 - f2*ddenom)/denom
      end if
      return
      end   












      subroutine bb_qromb(func,a,b,eps,ss)  
      include 'implno.dek' 

c..returns as ss the integral of the function func from a to b with fractional 
c..accuracy eps. integration by romberg's method of order 2k where e.g k=2 is  
c..simpson's rule.  
c..  
c..jmax limits the total number of steps; k is the 
c..the number of points used in the extrapolation; arrays s and h store the  
c..trapazoidal approximations and their relative step sizes. 

c..declare 
      external          func 
      integer           j,jmax,jmaxp,k,km    
      parameter         (jmax=20, jmaxp=jmax+1, k=5, km=k-1) 
      double precision  a,b,ss,s(jmaxp),h(jmaxp),eps,dss,func 

      h(1) = 1.0d0 
      do j=1,jmax 
       call bb_trapzd(func,a,b,s(j),j)  
       if (j .ge. k) then    
        call bb_polint(h(j-km),s(j-km),k,0.0d0,ss,dss)    
        if (abs(dss) .le. eps*abs(ss)) return    
       end if    
       s(j+1) = s(j) 
       h(j+1) = 0.25d0 * h(j)  
      enddo

c      write(6,*) ' after ',jmax,' iterations ' 
c      write(6,*) ' of trying to integrate between ',a,' and ',b 
c      write(6,*) ' and fractional accuracy ',eps 
c      write(6,*) ' the integral is ',ss 
c      write(6,*) ' and error estimate ',dss 
c      write(6,*) ' so that abs(dss) ',abs(dss), 
c     1           ' > eps*abs(ss)',eps*abs(ss) 
      ss = 0.0d0 
c      stop       'too many steps in qromb'  
      return 
      end    





      subroutine bb_trapzd(func,a,b,s,n)    
      include 'implno.dek' 

c..this routine computes the n'th stage of refinement of an extended 
c..trapazoidal rule. func is input as the name of a function to be   
c..integrated between limits a and b. when n=1 the routine returns as s  
c..the crudest estimate of the integral of func(x)dx from a to b.    
c..subsequent calls with n=2,3... will improve the accuracy of s by adding   
c..2**(n-2) additional interior points. s should not be modified between 
c..sequential calls. 
c.. 
c..this routine is the workhorse of all the following closed formula 
c..integration routines. 
c..  
c..local it  is the number of points to be added on the next call    
c..local del is the step size.   
c.. 
c..declare  
      external          func 
      integer           n,it,j   
      double precision  func,a,b,s,del,x,sum,tnm 

c..go 
      if (n.eq.1) then   
       s  = 0.5d0 * (b-a) * ( func(a) + func(b) )   
      else   
       it  = 2**(n-2) 
       tnm = it  
       del = (b-a)/tnm   
       x   = a + (0.5d0 *del)    
       sum = 0.0d0 
       do j=1,it  
        sum = sum + func(x)  
        x   = x + del  
       enddo
       s  = 0.5d0 * (s + (b-a)*sum/tnm) 
      end if 
      return     
      end    







      subroutine bb_polint(xa,ya,n,x,y,dy)
      include 'implno.dek'

c..given arrays xa and ya of length n and a value x, this routine returns a 
c..value y and an error estimate dy. if p(x) is the polynomial of degree n-1
c..such that ya = p(xa) then the returned value is y = p(x) 

c..declare
      integer          n,nmax,ns,i,m
      parameter        (nmax=10)
      double precision xa(n),ya(n),x,y,dy,c(nmax),d(nmax),dif,dift,
     1                 ho,hp,w,den


c..find the index ns of the closest table entry; initialize the c and d tables
      ns  = 1
      dif = abs(x - xa(1))
      do i=1,n
       dift = abs(x - xa(i))
       if (dift .lt. dif) then
        ns  = i
        dif = dift
       end if
       c(i)  = ya(i)
       d(i)  = ya(i)
      enddo

c..first guess for y
      y = ya(ns)

c..for each column of the table, loop over the c's and d's and update them
      ns = ns - 1
      do m=1,n-1
       do i=1,n-m
        ho   = xa(i) - x
        hp   = xa(i+m) - x
        w    = c(i+1) - d(i)
        den  = ho - hp
        if (den .eq. 0.0) stop ' 2 xa entries are the same in polint'
        den  = w/den
        d(i) = hp * den
        c(i) = ho * den
       enddo

c..after each column is completed, decide which correction c or d, to add
c..to the accumulating value of y, that is, which path to take in the table
c..by forking up or down. ns is updated as we go to keep track of where we
c..are. the last dy added is the error indicator.
       if (2*ns .lt. n-m) then
        dy = c(ns+1)
       else
        dy = d(ns)
        ns = ns - 1
       end if
       y = y + dy
      enddo
      return
      end





      subroutine bb_qgaus(func,x,w,n,ss)
      include 'implno.dek'

c..returns as ss the quadrature summation of the function func as determined 
c..by the abcissas x and weights w.

c..declare
      external          func
      integer           i,n
      double precision  func,ss,x(n),w(n)

      ss = 0.0d0
      do i=1,n
       ss = ss + w(i)*func(x(i))
      enddo
      return
      end   





      subroutine bb_gauleg(x1,x2,x,w,n)
      include 'implno.dek'
      include 'const.dek'

c..given the lower and upper limits of integration x1 and x2, and given n, 
c..this routine returns arrays x and w of length n, containing the 
c..abscissas and weights of the gauss-legendre n-point quadrature formula

c..declare
      integer            i,m,j,n
      double precision   x1,x2,x(n),w(n),eps,xm,xl,p1,p2,p3,pp,z,z1 
      parameter          (eps=1.0e-14)   


c..roots are symmetric in the interval so we only have to find half of them
      m = (n+1)/2   
      xm = 0.5d0 * (x2 + x1)
      xl = 0.5d0 * (x2 - x1)

c..loop over the desired roots and make a slick guess at each one   
      do i=1,m   
       z = cos(3.141592653589d0 * (i-0.25d0)/(n + 0.5d0))

c..newton do while loop  
1      continue 
       p1 = 1.0d0   
       p2 = 0.0d0   

c..loop the recurrence relation to get the legendre polynomial at z 
       do j=1,n  
        p3 = p2 
        p2 = p1 
        p1 = ((2.0d0 * j - 1.0d0) * z * p2 - (j - 1.0d0)*p3)/j
       enddo

c..p1 is now the desired legendre polynomial. pp is the derivative.
       pp = n * (z*p1 - p2)/(z*z - 1.0d0) 
       z1 = z   
       z  = z1 - p1/pp   
       if (abs(z-z1) .gt. eps) goto  1   

c..scale to the users interval  
       x(i)     = xm - xl*z 
       x(n+1-i) = xm + xl * z   
       w(i)     = 2.0d0 * xl/((1.0d0 - z*z)*pp*pp)  
       w(n+1-i) = w(i)  
      enddo
      return
      end










c..this file contains auxillary network routine

c..routine screen5 computes screening factors
c..routine sneut5 computes neutrino loss rates 
c..routine ifermi12 does an inverse fermi integral of order 1/2
c..routine zfermim12 does an inverse fermi integral of order -1/2

c..routine ecapnuc computes electron capture rates  








      subroutine screen5(temp,den,zbar,abar,z2bar,
     1                   z1,a1,z2,a2,jscreen,init,
     2                   scor,scordt,scordd)
      include 'implno.dek'

c..this subroutine calculates screening factors and their derivatives
c..for nuclear reaction rates in the weak, intermediate and strong regimes. 
c..based on graboske, dewit, grossman and cooper apj 181 457 1973 for 
c..weak screening. based on alastuey and jancovici apj 226 1034 1978, 
c..with plasma parameters from itoh et al apj 234 1079 1979, for strong 
c..screening. 

c..input:
c..temp    = temperature
c..den     = density
c..zbar    = mean charge per nucleus
c..abar    = mean number of nucleons per nucleus 
c..z2bar   = mean square charge per nucleus
c..z1 a1   = charge and number in the entrance channel
c..z2 a2   = charge and number in the exit channel
c..jscreen = counter of which reaction is being calculated 
c..init    = flag to compute the more expensive functions just once

c..output:
c..scor    = screening correction
c..scordt  = derivative of screening correction with temperature
c..scordd  = derivative of screening correction with density


c..declare the pass
      integer          jscreen,init
      double precision temp,den,zbar,abar,z2bar,z1,a1,z2,a2,
     1                 scor,scordt,scordd


c..local variables
      integer          i
      double precision aa,daadt,daadd,bb,cc,dccdt,dccdd,
     1                 pp,dppdt,dppdd,qq,dqqdt,dqqdd,rr,drrdt,drrdd,
     2                 ss,dssdt,dssdd,tt,dttdt,dttdd,uu,duudt,duudd,
     3                 vv,dvvdt,dvvdd,a3,da3,tempi,dtempi,deni,
     2                 qlam0z,qlam0zdt,qlam0zdd,
     3                 h12w,dh12wdt,dh12wdd,h12,dh12dt,dh12dd,
     4                 taufac,taufacdt,gamp,gampdt,gampdd,gampi,
     5                 gamef,gamefdt,gamefdd,
     6                 tau12,tau12dt,alph12,alph12dt,alph12dd,
     7                 xlgfac,dxlgfacdt,dxlgfacdd,
     8                 gamp14,gamp14dt,gamp14dd,
     9                 xni,dxnidd,ytot,
     &                 temp_old,den_old,zbar_old,abar_old


c..screening variables
c..zs13    = (z1+z2)**(1./3.)
c..zhat    = combination of z1 and z2 raised to the 5/3 power
c..zhat2   = combination of z1 and z2 raised to the 5/12 power
c..lzav    = log of effective charge
c..aznut   = combination of a1,z1,a2,z2 raised to 1/3 power

      integer          abigrat
      parameter        (abigrat   = 6000)

      double precision zs13(abigrat),zhat(abigrat),
     1                 zhat2(abigrat),lzav(abigrat),
     2                 aznut(abigrat),zs13inv(abigrat)


c..parameter fact is the cube root of 2 
      double precision  x13,x14,x53,x532,x512,fact,co2
      parameter        (x13   = 1.0d0/3.0d0,  
     1                  x14   = 1.0d0/4.0d0, 
     3                  x53   = 5.0d0/3.0d0, 
     4                  x532  = 5.0d0/32.0d0,
     5                  x512  = 5.0d0/12.0d0,
     6                  fact  = 1.25992104989487d0,
     7                  co2   = x13 * 4.248719d3)


      data     temp_old/-1.0d0/, den_old/-1.0d0/,
     1         zbar_old/-1.0d0/, abar_old/-1.0d0/ 




c..compute and store the more expensive screening factors
      if (init .eq. 1) then
       if (jscreen .gt. abigrat) stop 'jscreen > abigrat in screen5'
       zs13(jscreen)    = (z1 + z2)**x13
       zs13inv(jscreen) = 1.0d0/zs13(jscreen)
       zhat(jscreen)    = (z1 + z2)**x53  - z1**x53 - z2**x53
       zhat2(jscreen)   = (z1 + z2)**x512 - z1**x512 -z2**x512
       lzav(jscreen)    = x53 * log(z1*z2/(z1 + z2))
       aznut(jscreen)   = (z1**2 * z2**2 * a1*a2 / (a1 + a2))**x13
      endif


c..calculate average plasma, if need be
      if (temp_old .ne. temp .or.
     1    den_old  .ne. den  .or.  
     2    zbar_old  .ne. zbar  .or.  
     3    abar_old  .ne. abar ) then

       temp_old = temp
       den_old  = den
       zbar_old  = zbar
       abar_old  = abar

       ytot     = 1.0d0/abar
       rr       = den * ytot
       tempi   = 1.0d0/temp
       dtempi  = -tempi*tempi
       deni    = 1.0d0/den

       pp       = sqrt(rr*tempi*(z2bar + zbar)) 
       qq       = 0.5d0/pp *(z2bar + zbar) 
       dppdt    = qq*rr*dtempi
       dppdd    = qq*ytot*tempi

       qlam0z   = 1.88d8 * tempi * pp
       qlam0zdt = 1.88d8 * (dtempi*pp + tempi*dppdt)
       qlam0zdd = 1.88d8 * tempi * dppdd

       taufac   = co2 * tempi**x13 
       taufacdt = -x13*taufac*tempi

       qq      = rr*zbar
       xni     = qq**x13
       dxnidd  = x13 * xni * deni

       aa     = 2.27493d5 * tempi * xni
       daadt  = 2.27493d5 * dtempi * xni
       daadd  = 2.27493d5 * tempi * dxnidd
      end if


c..calculate individual screening factors 
      bb       = z1 * z2
      gamp     = aa
      gampdt   = daadt
      gampdd   = daadd

      qq       = fact * bb * zs13inv(jscreen)       
      gamef    = qq * gamp 
      gamefdt  = qq * gampdt
      gamefdd  = qq * gampdd

      tau12    = taufac * aznut(jscreen) 
      tau12dt  = taufacdt * aznut(jscreen)

      qq       = 1.0d0/tau12 
      alph12   = gamef * qq
      alph12dt = (gamefdt - alph12*tau12dt) * qq
      alph12dd = gamefdd * qq
      


c..limit alph12 to 1.6 to prevent unphysical behavior.  
c..this should really be replaced by a pycnonuclear reaction rate formula 
      if (alph12 .gt. 1.6) then 
       alph12   = 1.6d0 
       alph12dt = 0.0d0
       alph12dd = 0.0d0

       gamef    = 1.6d0 * tau12 
       gamefdt  = 1.6d0 * tau12dt
       gamefdd  = 0.0d0 

       qq       = zs13(jscreen)/(fact * bb) 
       gamp     = gamef * aa
       gampdt   = gamefdt * aa
       gampdd   = 0.0d0
      end if 



c..weak screening regime 
      h12w    = bb * qlam0z 
      dh12wdt = bb * qlam0zdt
      dh12wdd = bb * qlam0zdd

      h12     = h12w 
      dh12dt  = dh12wdt
      dh12dd  = dh12wdd



c..intermediate and strong sceening regime
      if (gamef .gt. 0.3) then 

       gamp14   = gamp**x14
       rr       = 1.0d0/gamp
       qq       = 0.25d0*gamp14*rr 
       gamp14dt = qq * gampdt
       gamp14dd = qq * gampdd

       cc       =   0.896434d0 * gamp * zhat(jscreen) 
     1            - 3.44740d0  * gamp14 * zhat2(jscreen)  
     2            - 0.5551d0   * (log(gamp) + lzav(jscreen)) 
     3            - 2.996d0 

       dccdt    =   0.896434d0 * gampdt * zhat(jscreen) 
     1            - 3.44740d0  * gamp14dt * zhat2(jscreen)  
     2            - 0.5551d0*rr*gampdt

       dccdd    =   0.896434d0 * gampdd * zhat(jscreen) 
     1            - 3.44740d0  * gamp14dd * zhat2(jscreen)  
     2            - 0.5551d0*rr*gampdd

       a3     = alph12 * alph12 * alph12 
       da3    = 3.0d0 * alph12 * alph12

       qq     = 0.014d0 + 0.0128d0*alph12
       dqqdt  = 0.0128d0*alph12dt
       dqqdd  = 0.0128d0*alph12dd

       rr     = x532 - alph12*qq
       drrdt  = -(alph12dt*qq + alph12*dqqdt) 
       drrdd  = -(alph12dd*qq + alph12*dqqdd) 

       ss     = tau12*rr
       dssdt  = tau12dt*rr + tau12*drrdt
       dssdd  = tau12*drrdd

       tt     =  -0.0098d0 + 0.0048d0*alph12
       dttdt  = 0.0048d0*alph12dt
       dttdd  = 0.0048d0*alph12dd
       
       uu     =  0.0055d0 + alph12*tt
       duudt  = alph12dt*tt + alph12*dttdt
       duudd  = alph12dd*tt + alph12*dttdd

       vv   = gamef * alph12 * uu  
       dvvdt= gamefdt*alph12*uu + gamef*alph12dt*uu + gamef*alph12*duudt  
       dvvdd= gamefdd*alph12*uu + gamef*alph12dd*uu + gamef*alph12*duudd  

       h12     = cc - a3 * (ss + vv)
       rr      = da3 * (ss + vv)
       dh12dt  = dccdt - rr*alph12dt - a3*(dssdt + dvvdt)
       dh12dd  = dccdd - rr*alph12dd - a3*(dssdd + dvvdd)

       rr     =  1.0d0 - 0.0562d0*a3
       ss     =  -0.0562d0*da3
       drrdt  = ss*alph12dt
       drrdd  = ss*alph12dd

       if (rr .ge. 0.77d0) then
        xlgfac    = rr
        dxlgfacdt = drrdt 
        dxlgfacdd = drrdd 
       else
        xlgfac    = 0.77d0
        dxlgfacdt = 0.0d0
        dxlgfacdd = 0.0d0
       end if 


       h12    = log(xlgfac) + h12 
       rr     = 1.0d0/xlgfac
       dh12dt = rr*dxlgfacdt + dh12dt 
       dh12dd = rr*dxlgfacdd + dh12dd 


       if (gamef .le. 0.8) then 
        rr     = 2.0d0*(0.8d0-gamef) 
        drrdt  = -2.0d0*gamefdt
        drrdd  = -2.0d0*gamefdd

        ss     = 2.0d0*(gamef-0.3d0) 
        dssdt  = 2.0d0*gamefdt  
        dssdd  = 2.0d0*gamefdd  

        vv     = h12 
        h12    = h12w*rr + vv*ss
        dh12dt = dh12wdt*rr + h12w*drrdt + dh12dt*ss + vv*dssdt
        dh12dd = dh12wdd*rr + h12w*drrdd + dh12dd*ss + vv*dssdd
       end if 

c..end of intermediate and strong screening if
      end if 


c..machine limit the output
      h12    = max(min(h12,300.0d0),0.0d0) 
      scor   = exp(h12) 
      if (h12 .eq. 300.0d0) then
       scordt = 0.0d0
       scordd = 0.0d0
      else 
       scordt = scor * dh12dt
       scordd = scor * dh12dd
      end if

c      write(6,111) 'weak =',h12w,' total =',h12,
c     1             ' 1-ratio =',1.0d0-h12w/h12,' correction',scor
c 111  format(1x,4(a,1pe13.6))
c      read(5,*)

      return 
      end 









      double precision function snupp(yp,ratepp,ybe7,ratebeec,
     1                                yb8,rateb8epnu)
      include 'implno.dek'
      include 'const.dek'

c..computes approximate neutrino losses from pp chain reactions
c..see page 142 of astro 289j notes for these loss formulas

c..input:
c..yp         = proton molar abbundance
c..ratepp     = pp reaction rate
c..ybe7       = be7 molar abundance
c..ratebeec   = be7 electron capture reaction rate
c..yb8        = b8 molar abundance
c..rateb8epnu = b8 decay reaction rate


c..declare the pass
      double precision yp,ratepp,ybe7,ratebeec,yb8,rateb8epnu


c..local variables
      double precision pp1nu,pp2nu,pp3nu,conv
      parameter        (conv = ev2erg*1.0d6*avo)


c..nu losses from p(p,e-nu)h2
      pp1nu  = yp*yp*ratepp * 0.5d0 * 0.263d0


c..nu losses from be7(n=>p)li7  
      pp2nu  = ybe7 * ratebeec * 0.81d0


c..nu losses from b8(p=>n)be8=>2a
      pp3nu  = yb8 * rateb8epnu * 7.73d0

c..sum the pp-chain neutrino losses and convert to erg/g/s
      snupp  = (pp1nu + pp2nu + pp3nu) * conv

      return
      end




      double precision function snucno(yn13,bc13,bn13,yo14,bn14,bo14,
     1                                 yo15,bn15,bo15,yf17,bo17,bf17,
     2                                 yf18,bo18,bf18)
      include 'implno.dek'
      include 'const.dek'

c..computes approximate neutrino losses from cno cycle  reactions
c..see page 142 of astro 289j notes for these loss formulas

c..input:
c..yn13 = n13 molar abundance
c..bc13 = c13 binding energy in mev
c..bn13 = n13 binding energy in mev
c..yo14 = o14 molar abundance
c..bn14 = n14 binding energy in mev
c..bo14 = o14 binding energy in mev
c..yo15 = o15 molar abundance
c..bn15 = n15 binding energy in mev
c..bo15 = o15 binding energy in mev
c..yf17 = f17 molar abundance
c..bo17 = o17 binding energy in mev
c..bf17 = f17 binding energy in mev
c..yf18 = f18 molar abundance
c..bo18 = o18 binding energy in mev
c..bf18 = f18 binding energy in mev


c..declare the pass
      double precision yn13,bc13,bn13,yo14,bn14,bo14,
     1                 yo15,bn15,bo15,yf17,bo17,bf17,
     2                 yf18,bo18,bf18

c..local variables
      double precision sum,sum2,enu13n,enu14o,enu15o,enu17f,enu18f,
     1                 conv,lntwo,tm1,tm2,tm3,tm4,tm5
      parameter        (conv  = ev2erg*1.0d6*avo,
     1                  lntwo = 0.693147181d0,
     2                  tm1   = lntwo/597.9d0, 
     3                  tm2   = lntwo/70.606d0,
     4                  tm3   = lntwo/124.0,
     5                  tm4   = lntwo/64.49,
     6                  tm5   = lntwo/6586.2)


c..13n(e+nu)13c 
      sum    = bc13 - bn13 - 0.782d0 - 1.022d0
      sum    = 1.0d0 + sum/0.511d0
      sum2   = sum*sum
      enu13n = 0.5d0 * sum * 0.511d0 * (1.0d0 - 1.0d0/sum2) 
     1         * (1.0d0 - 1.0d0/(4.0d0*sum) - 1.0d0/(9.0d0*sum2))
      enu13n = yn13 * enu13n * tm1


c..hot cno cycle 14o(e+nu)14n 
      sum    = bn14 - bo14 - 0.782d0 - 1.022d0
      sum    = 1.0d0 + sum/0.511d0
      sum2   = sum*sum
      enu14o = 0.5d0 * sum * 0.511d0 * (1.0d0 - 1.0d0/sum2) 
     1         * (1.0d0 - 1.0d0/(4.0d0*sum) - 1.0d0/(9.0d0*sum2))
      enu14o = yo14 * enu14o * tm2


c..15o(e+nu)15n 
      sum    = bn15 - bo15 - 0.782d0 - 1.022d0
      sum    = 1.0d0 + sum/0.511d0
      sum2   = sum*sum
      enu15o = 0.5d0 * sum * 0.511d0 * (1.0d0 - 1.0d0/sum2) 
     1         * (1.0d0 - 1.0d0/(4.0d0*sum) - 1.0d0/(9.0d0*sum2))
      enu15o = yo15 * enu15o * tm3


c..17f(e+nu)17o 
      sum    = bo17 - bf17 - 0.782d0 - 1.022d0
      sum    = 1.0d0 + sum/0.511d0
      sum2   = sum*sum
      enu17f = 0.5d0 * sum * 0.511d0 * (1.0d0 - 1.0d0/sum2) 
     1         * (1.0d0 - 1.0d0/(4.0d0*sum) - 1.0d0/(9.0d0*sum2))
      enu17f  = yf17 * enu17f * tm4


c..18f(e+nu)18o 
      sum    = bo18 - bf18 - 0.782d0 - 1.022d0
      sum    = 1.0d0 + sum/0.511d0
      sum2   = sum*sum
      enu18f = 0.5d0 * sum * 0.511d0 * (1.0d0 - 1.0d0/sum2) 
     1         * (1.0d0 - 1.0d0/(4.0d0*sum) - 1.0d0/(9.0d0*sum2))
      enu18f = yf18 * enu18f * tm5


c..sum the cno cycle losses and convert to erg/g/s
      snucno = (enu13n + enu14o + enu15o + enu17f + enu18f) * conv

      return 
      end








      subroutine sneut5(temp,den,abar,zbar,
     1                  snu,dsnudt,dsnudd,dsnuda,dsnudz)
      include 'implno.dek'
      include 'const.dek'

c..this routine computes neutrino losses from the analytic fits of
c..itoh et al. apjs 102, 411, 1996, and also returns their derivatives. 

c..input:
c..temp = temperature 
c..den  = density
c..abar = mean atomic weight
c..zbar = mean charge

c..output:
c..snu    = total neutrino loss rate in erg/g/sec
c..dsnudt = derivative of snu with temperature
c..dsnudd = derivative of snu with density
c..dsnuda = derivative of snu with abar
c..dsnudz = derivative of snu with zbar


c..declare the pass
      double precision temp,den,abar,zbar,
     1                 snu,dsnudt,dsnudd,dsnuda,dsnudz

c..local variables
      integer          i
      double precision spair,spairdt,spairdd,spairda,spairdz,
     1                 splas,splasdt,splasdd,splasda,splasdz,
     2                 sphot,sphotdt,sphotdd,sphotda,sphotdz,
     3                 sbrem,sbremdt,sbremdd,sbremda,sbremdz,
     4                 sreco,srecodt,srecodd,srecoda,srecodz

      double precision t9,xl,xldt,xlp5,xl2,xl3,xl4,xl5,xl6,xl7,xl8,xl9,
     1                 xlmp5,xlm1,xlm2,xlm3,xlm4,xlnt,cc,den6,tfermi,
     2                 a0,a1,a2,a3,b1,b2,b3,c00,c01,c02,c03,c04,c05,c06,
     3                 c10,c11,c12,c13,c14,c15,c16,c20,c21,c22,c23,c24,
     4                 c25,c26,dd00,dd01,dd02,dd03,dd04,dd05,dd11,dd12,
     5                 dd13,dd14,dd15,dd21,dd22,dd23,dd24,dd25,b,c,d,f0,
     6                 f1,deni,tempi,abari,zbari,f2,f3,z,deni,xmue,ye,
     7                 rp1,rn1,dum,dumdt,dumdd,dumda,dumdz,
     8                 gum,gumdt,gumdd,gumda,gumdz


c..pair production
      double precision rm,rmdd,rmda,rmdz,rmi,gl,gldt,
     1                 zeta,zetadt,zetadd,zetada,zetadz,zeta2,zeta3,
     2                 xnum,xnumdt,xnumdd,xnumda,xnumdz,
     3                 xden,xdendt,xdendd,xdenda,xdendz,
     4                 fpair,fpairdt,fpairdd,fpairda,fpairdz,
     5                 qpair,qpairdt,qpairdd,qpairda,qpairdz

c..plasma 
      double precision gl2,gl2dt,gl2dd,gl2da,gl2dz,gl12,gl32,gl72,gl6,
     1                 ft,ftdt,ftdd,ftda,ftdz,fl,fldt,fldd,flda,fldz,
     2                 fxy,fxydt,fxydd,fxyda,fxydz

c..photo
      double precision tau,taudt,cos1,cos2,cos3,cos4,cos5,sin1,sin2,
     1                 sin3,sin4,sin5,last,xast,
     2                 fphot,fphotdt,fphotdd,fphotda,fphotdz,
     3                 qphot,qphotdt,qphotdd,qphotda,qphotdz

c..brem
      double precision t8,t812,t832,t82,t83,t85,t86,t8m1,t8m2,t8m3,t8m5,
     1                 t8m6,
     2                 eta,etadt,etadd,etada,etadz,etam1,etam2,etam3,
     3                 fbrem,fbremdt,fbremdd,fbremda,fbremdz,
     4                 gbrem,gbremdt,gbremdd,gbremda,gbremdz,
     5                 u,gm1,gm2,gm13,gm23,gm43,gm53,v,w,fb,gt,gb,
     6                 fliq,fliqdt,fliqdd,fliqda,fliqdz, 
     7                 gliq,gliqdt,gliqdd,gliqda,gliqdz 

c..recomb
      double precision ifermi12,zfermim12,nu,nudt,nudd,nuda,nudz,
     1                 nu2,nu3,bigj,bigjdt,bigjdd,bigjda,bigjdz



c..numerical constants
      double precision fac1,fac2,fac3,oneth,twoth,con1,sixth,iln10
      parameter        (fac1   = 5.0d0 * pi / 3.0d0,
     2                  fac2   = 10.0d0 * pi,
     3                  fac3   = pi / 5.0d0,
     4                  oneth  = 1.0d0/3.0d0,
     5                  twoth  = 2.0d0/3.0d0,
     6                  con1   = 1.0d0/5.9302d0,
     7                  sixth  = 1.0d0/6.0d0,
     8                  iln10  = 4.342944819032518d-1)


c..theta is sin**2(theta_weinberg) = 0.2319 plus/minus 0.00005 (1996)
c..xnufam is the number of neutrino flavors = 3.02 plus/minus 0.005 (1998)
c..change theta and xnufam if need be, and the changes will automatically
c..propagate through the routine. cv and ca are the vector and axial currents.

      double precision theta,xnufam,cv,ca,cvp,cap,tfac1,tfac2,tfac3,
     1                 tfac4,tfac5,tfac6
      parameter        (theta  = 0.2319d0,
     1                  xnufam = 3.0d0,
     2                  cv     = 0.5d0 + 2.0d0 * theta,
     3                  cvp    = 1.0d0 - cv,
     4                  ca     = 0.5d0,
     5                  cap    = 1.0d0 - ca,
     6                  tfac1  = cv*cv + ca*ca + 
     7                           (xnufam-1.0d0) * (cvp*cvp+cap*cap),
     8                  tfac2  = cv*cv - ca*ca + 
     9                           (xnufam-1.0d0) * (cvp*cvp - cap-cap),
     &                  tfac3  = tfac2/tfac1,
     1                  tfac4  = 0.5d0 * tfac1,
     2                  tfac5  = 0.5d0 * tfac2,
     3                  tfac6  = cv*cv + 1.5d0*ca*ca + (xnufam - 1.0d0)*
     4                           (cvp*cvp + 1.5d0*cap*cap))



c..initialize 
      spair   = 0.0d0
      spairdt = 0.0d0
      spairdd = 0.0d0
      spairda = 0.0d0
      spairdz = 0.0d0

      splas   = 0.0d0
      splasdt = 0.0d0
      splasdd = 0.0d0
      splasda = 0.0d0
      splasdz = 0.0d0

      sphot   = 0.0d0
      sphotdt = 0.0d0
      sphotdd = 0.0d0
      sphotda = 0.0d0
      sphotdz = 0.0d0

      sbrem   = 0.0d0
      sbremdt = 0.0d0
      sbremdd = 0.0d0
      sbremda = 0.0d0
      sbremdz = 0.0d0

      sreco   = 0.0d0
      srecodt = 0.0d0
      srecodd = 0.0d0
      srecoda = 0.0d0
      srecodz = 0.0d0

      snu     = 0.0d0
      dsnudt  = 0.0d0
      dsnudd  = 0.0d0
      dsnuda  = 0.0d0
      dsnudz  = 0.0d0

      if (temp .lt. 1.0e7) return


c..to avoid lots of divisions
      deni  = 1.0d0/den
      tempi = 1.0d0/temp
      abari = 1.0d0/abar
      zbari = 1.0d0/zbar


c..some composition variables
      ye    = zbar*abari
      xmue  = abar*zbari




c..some frequent factors
      t9     = temp * 1.0d-9
      xl     = t9 * con1
      xldt   = 1.0d-9 * con1
      xlp5   = sqrt(xl)
      xl2    = xl*xl
      xl3    = xl2*xl
      xl4    = xl3*xl
      xl5    = xl4*xl
      xl6    = xl5*xl
      xl7    = xl6*xl
      xl8    = xl7*xl
      xl9    = xl8*xl
      xlmp5  = 1.0d0/xlp5
      xlm1   = 1.0d0/xl
      xlm2   = xlm1*xlm1
      xlm3   = xlm1*xlm2
      xlm4   = xlm1*xlm3

      rm     = den*ye
      rmdd   = ye
      rmda   = -rm*abari
      rmdz   = den*abari
      rmi    = 1.0d0/rm

      a0     = rm * 1.0d-9
      a1     = a0**oneth 
      zeta   = a1 * xlm1
      zetadt = -a1 * xlm2 * xldt
      a2     = oneth * a1*rmi * xlm1
      zetadd = a2 * rmdd 
      zetada = a2 * rmda
      zetadz = a2 * rmdz
      
      zeta2 = zeta * zeta
      zeta3 = zeta2 * zeta




c..pair neutrino section
c..for reactions like e+ + e- => nu_e + nubar_e 

c..equation 2.8 
      gl   = 1.0d0 - 13.04d0*xl2 +133.5d0*xl4 +1534.0d0*xl6 +918.6d0*xl8
      gldt = xldt*(-26.08d0*xl +534.0d0*xl3 +9204.0d0*xl5 +7348.8d0*xl7)

c..equation 2.7

      a1     = 6.002d19 + 2.084d20*zeta + 1.872d21*zeta2
      a2     = 2.084d20 + 2.0d0*1.872d21*zeta

      if (t9 .lt. 10.0) then
       b1     = exp(-5.5924d0*zeta)
       b2     = -b1*5.5924d0
      else
       b1     = exp(-4.9924d0*zeta)
       b2     = -b1*4.9924d0
      end if
      
      xnum   = a1 * b1
      c      = a2*b1 + a1*b2
      xnumdt = c*zetadt
      xnumdd = c*zetadd
      xnumda = c*zetada
      xnumdz = c*zetadz

      if (t9 .lt. 10.0) then
       a1   = 9.383d-1*xlm1 - 4.141d-1*xlm2 + 5.829d-2*xlm3
       a2   = -9.383d-1*xlm2 + 2.0d0*4.141d-1*xlm3 - 3.0d0*5.829d-2*xlm4
      else
       a1   = 1.2383d0*xlm1 - 8.141d-1*xlm2 
       a2   = -1.2383d0*xlm2 + 2.0d0*8.141d-1*xlm3 
      end if

      b1   = 3.0d0*zeta2

      xden   = zeta3 + a1
      xdendt = b1*zetadt + a2*xldt
      xdendd = b1*zetadd
      xdenda = b1*zetada
      xdendz = b1*zetadz

      a1      = 1.0d0/xden
      fpair   = xnum*a1
      fpairdt = (xnumdt - fpair*xdendt)*a1
      fpairdd = (xnumdd - fpair*xdendd)*a1
      fpairda = (xnumda - fpair*xdenda)*a1
      fpairdz = (xnumdz - fpair*xdendz)*a1


c..equation 2.6
      a1     = 10.7480d0*xl2 + 0.3967d0*xlp5 + 1.005d0
      a2     = xldt*(2.0d0*10.7480d0*xl + 0.5d0*0.3967d0*xlmp5) 
      xnum   = 1.0d0/a1
      xnumdt = -xnum*xnum*a2

      a1     = 7.692d7*xl3 + 9.715d6*xlp5
      a2     = xldt*(3.0d0*7.692d7*xl2 + 0.5d0*9.715d6*xlmp5)

      c      = 1.0d0/a1
      b1     = 1.0d0 + rm*c

      xden   = b1**(-0.3d0)

      d      = -0.3d0*xden/b1
      xdendt = -d*rm*c*c*a2
      xdendd = d*rmdd*c 
      xdenda = d*rmda*c 
      xdendz = d*rmdz*c 

      qpair   = xnum*xden
      qpairdt = xnumdt*xden + xnum*xdendt
      qpairdd = xnum*xdendd
      qpairda = xnum*xdenda
      qpairdz = xnum*xdendz



c..equation 2.5
      a1    = exp(-2.0d0*xlm1)
      a2    = a1*2.0d0*xlm2*xldt

      spair   = a1*fpair
      spairdt = a2*fpair + a1*fpairdt
      spairdd = a1*fpairdd
      spairda = a1*fpairda
      spairdz = a1*fpairdz

      a1      = spair
      spair   = gl*a1
      spairdt = gl*spairdt + gldt*a1
      spairdd = gl*spairdd
      spairda = gl*spairda
      spairdz = gl*spairdz

      a1      = tfac4*(1.0d0 + tfac3 * qpair)
      a2      = tfac4*tfac3

      a3      = spair
      spair   = a1*a3
      spairdt = a1*spairdt + a2*qpairdt*a3
      spairdd = a1*spairdd + a2*qpairdd*a3
      spairda = a1*spairda + a2*qpairda*a3
      spairdz = a1*spairdz + a2*qpairdz*a3




c..plasma neutrino section 
c..for collective reactions like gamma_plasmon => nu_e + nubar_e
c..equation 4.6

      a1   = 1.019d-6*rm
      a2   = a1**twoth
      a3   = twoth*a2/a1

      b1   =  sqrt(1.0d0 + a2)
      b2   = 1.0d0/b1
  
      c00  = 1.0d0/(temp*temp*b1)

      gl2   = 1.1095d11 * rm * c00

      gl2dt = -2.0d0*gl2*tempi
      d     = rm*c00*b2*0.5d0*b2*a3*1.019d-6
      gl2dd = 1.1095d11 * (rmdd*c00  - d*rmdd)
      gl2da = 1.1095d11 * (rmda*c00  - d*rmda)
      gl2dz = 1.1095d11 * (rmdz*c00  - d*rmdz)
      

      gl    = sqrt(gl2)
      gl12  = sqrt(gl)
      gl32  = gl * gl12
      gl72  = gl2 * gl32
      gl6   = gl2 * gl2 * gl2


c..equation 4.7
      ft   = 2.4d0 + 0.6d0*gl12 + 0.51d0*gl + 1.25d0*gl32
      gum  = 1.0d0/gl2
      a1   =(0.25d0*0.6d0*gl12 +0.5d0*0.51d0*gl +0.75d0*1.25d0*gl32)*gum
      ftdt = a1*gl2dt
      ftdd = a1*gl2dd
      ftda = a1*gl2da
      ftdz = a1*gl2dz


c..equation 4.8
      a1   = 8.6d0*gl2 + 1.35d0*gl72
      a2   = 8.6d0 + 1.75d0*1.35d0*gl72*gum

      b1   = 225.0d0 - 17.0d0*gl + gl2
      b2   = -0.5d0*17.0d0*gl*gum + 1.0d0

      c    = 1.0d0/b1
      fl   = a1*c

      d    = (a2 - fl*b2)*c       
      fldt = d*gl2dt
      fldd = d*gl2dd
      flda = d*gl2da
      fldz = d*gl2dz
     

c..equation 4.9 and 4.10
      cc   = log10(2.0d0*rm)
      xlnt = log10(temp)

      xnum   = sixth * (17.5d0 + cc - 3.0d0*xlnt)
      xnumdt = -iln10*0.5d0*tempi
      a2     = iln10*sixth*rmi
      xnumdd = a2*rmdd
      xnumda = a2*rmda 
      xnumdz = a2*rmdz 

      xden   = sixth * (-24.5d0 + cc + 3.0d0*xlnt)
      xdendt = iln10*0.5d0*tempi
      xdendd = a2*rmdd
      xdenda = a2*rmda 
      xdendz = a2*rmdz 


c..equation 4.11
      if (abs(xnum) .gt. 0.7d0  .or.  xden .lt. 0.0d0) then
       fxy   = 1.0d0
       fxydt = 0.0d0
       fxydd = 0.0d0
       fxydz = 0.0d0
       fxyda = 0.0d0

      else 

       a1  = 0.39d0 - 1.25d0*xnum - 0.35d0*sin(4.5d0*xnum)
       a2  = -1.25d0 - 4.5d0*0.35d0*cos(4.5d0*xnum)

       b1  = 0.3d0 * exp(-1.0d0*(4.5d0*xnum + 0.9d0)**2)
       b2  = -b1*2.0d0*(4.5d0*xnum + 0.9d0)*4.5d0

       c   = min(0.0d0, xden - 1.6d0 + 1.25d0*xnum)
       if (c .eq. 0.0) then
        dumdt = 0.0d0
        dumdd = 0.0d0
        dumda = 0.0d0
        dumdz = 0.0d0
       else
        dumdt = xdendt + 1.25d0*xnumdt
        dumdd = xdendd + 1.25d0*xnumdd
        dumda = xdenda + 1.25d0*xnumda
        dumdz = xdendz + 1.25d0*xnumdz
       end if

       d   = 0.57d0 - 0.25d0*xnum
       a3  = c/d
       c00 = exp(-1.0d0*a3**2)

       f1  = -c00*2.0d0*a3/d
       c01 = f1*(dumdt + a3*0.25d0*xnumdt)
       c02 = f1*(dumdd + a3*0.25d0*xnumdd)
       c03 = f1*(dumda + a3*0.25d0*xnumda)
       c04 = f1*(dumdz + a3*0.25d0*xnumdz)

       fxy   = 1.05d0 + (a1 - b1)*c00
       fxydt = (a2*xnumdt -  b2*xnumdt)*c00 + (a1-b1)*c01
       fxydd = (a2*xnumdd -  b2*xnumdd)*c00 + (a1-b1)*c02
       fxyda = (a2*xnumda -  b2*xnumda)*c00 + (a1-b1)*c03
       fxydz = (a2*xnumdz -  b2*xnumdz)*c00 + (a1-b1)*c04

      end if



c..equation 4.1 and 4.5
      splas   = (ft + fl) * fxy
      splasdt = (ftdt + fldt)*fxy + (ft+fl)*fxydt
      splasdd = (ftdd + fldd)*fxy + (ft+fl)*fxydd
      splasda = (ftda + flda)*fxy + (ft+fl)*fxyda
      splasdz = (ftdz + fldz)*fxy + (ft+fl)*fxydz

      a2      = exp(-gl)
      a3      = -0.5d0*a2*gl*gum

      a1      = splas
      splas   = a2*a1
      splasdt = a2*splasdt + a3*gl2dt*a1
      splasdd = a2*splasdd + a3*gl2dd*a1
      splasda = a2*splasda + a3*gl2da*a1
      splasdz = a2*splasdz + a3*gl2dz*a1

      a2      = gl6
      a3      = 3.0d0*gl6*gum

      a1      = splas
      splas   = a2*a1
      splasdt = a2*splasdt + a3*gl2dt*a1
      splasdd = a2*splasdd + a3*gl2dd*a1
      splasda = a2*splasda + a3*gl2da*a1
      splasdz = a2*splasdz + a3*gl2dz*a1


      a2      = 0.93153d0 * 3.0d21 * xl9
      a3      = 0.93153d0 * 3.0d21 * 9.0d0*xl8*xldt

      a1      = splas
      splas   = a2*a1
      splasdt = a2*splasdt + a3*a1
      splasdd = a2*splasdd 
      splasda = a2*splasda 
      splasdz = a2*splasdz 




c..photoneutrino process section  
c..for reactions like e- + gamma => e- + nu_e + nubar_e
c..                   e+ + gamma => e+ + nu_e + nubar_e
c..equation 3.8 for tau, equation 3.6 for cc,
c..and table 2 written out for speed
      if (temp .ge. 1.0d7  .and. temp .lt. 1.0d8) then
       tau  =  log10(temp * 1.0d-7)
       cc   =  0.5654d0 + tau
       c00  =  1.008d11
       c01  =  0.0d0
       c02  =  0.0d0
       c03  =  0.0d0
       c04  =  0.0d0
       c05  =  0.0d0
       c06  =  0.0d0
       c10  =  8.156d10
       c11  =  9.728d8
       c12  = -3.806d9
       c13  = -4.384d9
       c14  = -5.774d9
       c15  = -5.249d9
       c16  = -5.153d9
       c20  =  1.067d11
       c21  = -9.782d9 
       c22  = -7.193d9
       c23  = -6.936d9
       c24  = -6.893d9
       c25  = -7.041d9
       c26  = -7.193d9
       dd01 =  0.0d0
       dd02 =  0.0d0
       dd03 =  0.0d0
       dd04 =  0.0d0
       dd05 =  0.0d0
       dd11 = -1.879d10
       dd12 = -9.667d9
       dd13 = -5.602d9
       dd14 = -3.370d9
       dd15 = -1.825d9
       dd21 = -2.919d10
       dd22 = -1.185d10
       dd23 = -7.270d9
       dd24 = -4.222d9
       dd25 = -1.560d9

      else if (temp .ge. 1.0d8  .and. temp .lt. 1.0d9) then
       tau   =  log10(temp * 1.0d-8)
       cc   =  1.5654d0
       c00  =  9.889d10 
       c01  = -4.524d8
       c02  = -6.088d6 
       c03  =  4.269d7 
       c04  =  5.172d7 
       c05  =  4.910d7 
       c06  =  4.388d7
       c10  =  1.813d11
       c11  = -7.556d9 
       c12  = -3.304d9  
       c13  = -1.031d9
       c14  = -1.764d9  
       c15  = -1.851d9
       c16  = -1.928d9
       c20  =  9.750d10
       c21  =  3.484d10
       c22  =  5.199d9  
       c23  = -1.695d9  
       c24  = -2.865d9  
       c25  = -3.395d9  
       c26  = -3.418d9
       dd01 = -1.135d8   
       dd02 =  1.256d8   
       dd03 =  5.149d7   
       dd04 =  3.436d7   
       dd05 =  1.005d7
       dd11 =  1.652d9  
       dd12 = -3.119d9  
       dd13 = -1.839d9  
       dd14 = -1.458d9  
       dd15 = -8.956d8
       dd21 = -1.549d10  
       dd22 = -9.338d9  
       dd23 = -5.899d9  
       dd24 = -3.035d9  
       dd25 = -1.598d9

      else if (temp .ge. 1.0d9) then
       tau  =  log10(t9)
       cc   =  1.5654d0
       c00  =  9.581d10
       c01  =  4.107d8
       c02  =  2.305d8   
       c03  =  2.236d8   
       c04  =  1.580d8   
       c05  =  2.165d8   
       c06  =  1.721d8
       c10  =  1.459d12
       c11  =  1.314d11
       c12  = -1.169d11  
       c13  = -1.765d11  
       c14  = -1.867d11  
       c15  = -1.983d11  
       c16  = -1.896d11
       c20  =  2.424d11
       c21  = -3.669d9
       c22  = -8.691d9  
       c23  = -7.967d9  
       c24  = -7.932d9  
       c25  = -7.987d9  
       c26  = -8.333d9
       dd01 =  4.724d8
       dd02 =  2.976d8   
       dd03 =  2.242d8   
       dd04 =  7.937d7   
       dd05 =  4.859d7
       dd11 = -7.094d11
       dd12 = -3.697d11
       dd13 = -2.189d11  
       dd14 = -1.273d11  
       dd15 = -5.705d10
       dd21 = -2.254d10
       dd22 = -1.551d10
       dd23 = -7.793d9
       dd24 = -4.489d9
       dd25 = -2.185d9
      end if

      taudt = iln10*tempi


c..equation 3.7, compute the expensive trig functions only one time
      cos1 = cos(fac1*tau)
      cos2 = cos(fac1*2.0d0*tau)
      cos3 = cos(fac1*3.0d0*tau)
      cos4 = cos(fac1*4.0d0*tau)
      cos5 = cos(fac1*5.0d0*tau)
      last = cos(fac2*tau)

      sin1 = sin(fac1*tau)
      sin2 = sin(fac1*2.0d0*tau)
      sin3 = sin(fac1*3.0d0*tau)
      sin4 = sin(fac1*4.0d0*tau)
      sin5 = sin(fac1*5.0d0*tau)
      xast = sin(fac2*tau)

      a0 = 0.5d0*c00 
     1     + c01*cos1 + dd01*sin1 + c02*cos2 + dd02*sin2
     2     + c03*cos3 + dd03*sin3 + c04*cos4 + dd04*sin4
     3     + c05*cos5 + dd05*sin5 + 0.5d0*c06*last

      f0 =  taudt*fac1*(-c01*sin1 + dd01*cos1 - c02*sin2*2.0d0 
     1     + dd02*cos2*2.0d0 - c03*sin3*3.0d0 + dd03*cos3*3.0d0 
     2     - c04*sin4*4.0d0 + dd04*cos4*4.0d0
     3     - c05*sin5*5.0d0 + dd05*cos5*5.0d0) 
     4     - 0.5d0*c06*xast*fac2*taudt

      a1 = 0.5d0*c10 
     1     + c11*cos1 + dd11*sin1 + c12*cos2 + dd12*sin2
     2     + c13*cos3 + dd13*sin3 + c14*cos4 + dd14*sin4
     3     + c15*cos5 + dd15*sin5 + 0.5d0*c16*last

      f1 = taudt*fac1*(-c11*sin1 + dd11*cos1 - c12*sin2*2.0d0 
     1     + dd12*cos2*2.0d0 - c13*sin3*3.0d0 + dd13*cos3*3.0d0 
     2     - c14*sin4*4.0d0 + dd14*cos4*4.0d0 - c15*sin5*5.0d0 
     3     + dd15*cos5*5.0d0) - 0.5d0*c16*xast*fac2*taudt

      a2 = 0.5d0*c20 
     1     + c21*cos1 + dd21*sin1 + c22*cos2 + dd22*sin2
     2     + c23*cos3 + dd23*sin3 + c24*cos4 + dd24*sin4
     3     + c25*cos5 + dd25*sin5 + 0.5d0*c26*last

      f2 = taudt*fac1*(-c21*sin1 + dd21*cos1 - c22*sin2*2.0d0 
     1     + dd22*cos2*2.0d0 - c23*sin3*3.0d0 + dd23*cos3*3.0d0 
     2     - c24*sin4*4.0d0 + dd24*cos4*4.0d0 - c25*sin5*5.0d0 
     3     + dd25*cos5*5.0d0) - 0.5d0*c26*xast*fac2*taudt

c..equation 3.4
      dum   = a0 + a1*zeta + a2*zeta2
      dumdt = f0 + f1*zeta + a1*zetadt + f2*zeta2 + 2.0d0*a2*zeta*zetadt
      dumdd = a1*zetadd + 2.0d0*a2*zeta*zetadd
      dumda = a1*zetada + 2.0d0*a2*zeta*zetada
      dumdz = a1*zetadz + 2.0d0*a2*zeta*zetadz

      z      = exp(-cc*zeta)

      xnum   = dum*z
      xnumdt = dumdt*z - dum*z*cc*zetadt
      xnumdd = dumdd*z - dum*z*cc*zetadd
      xnumda = dumda*z - dum*z*cc*zetada
      xnumdz = dumdz*z - dum*z*cc*zetadz

      xden   = zeta3 + 6.290d-3*xlm1 + 7.483d-3*xlm2 + 3.061d-4*xlm3

      dum    = 3.0d0*zeta2
      xdendt = dum*zetadt - xldt*(6.290d-3*xlm2 
     1         + 2.0d0*7.483d-3*xlm3 + 3.0d0*3.061d-4*xlm4)
      xdendd = dum*zetadd
      xdenda = dum*zetada
      xdendz = dum*zetadz

      dum      = 1.0d0/xden
      fphot   = xnum*dum
      fphotdt = (xnumdt - fphot*xdendt)*dum
      fphotdd = (xnumdd - fphot*xdendd)*dum
      fphotda = (xnumda - fphot*xdenda)*dum
      fphotdz = (xnumdz - fphot*xdendz)*dum
  

c..equation 3.3
      a0     = 1.0d0 + 2.045d0 * xl
      xnum   = 0.666d0*a0**(-2.066d0)
      xnumdt = -2.066d0*xnum/a0 * 2.045d0*xldt

      dum    = 1.875d8*xl + 1.653d8*xl2 + 8.449d8*xl3 - 1.604d8*xl4
      dumdt  = xldt*(1.875d8 + 2.0d0*1.653d8*xl + 3.0d0*8.449d8*xl2 
     1         - 4.0d0*1.604d8*xl3)

      z      = 1.0d0/dum
      xden   = 1.0d0 + rm*z
      xdendt =  -rm*z*z*dumdt
      xdendd =  rmdd*z
      xdenda =  rmda*z
      xdendz =  rmdz*z

      z      = 1.0d0/xden
      qphot = xnum*z
      qphotdt = (xnumdt - qphot*xdendt)*z
      dum      = -qphot*z
      qphotdd = dum*xdendd
      qphotda = dum*xdenda
      qphotdz = dum*xdendz

c..equation 3.2
      sphot   = xl5 * fphot
      sphotdt = 5.0d0*xl4*xldt*fphot + xl5*fphotdt
      sphotdd = xl5*fphotdd
      sphotda = xl5*fphotda
      sphotdz = xl5*fphotdz

      a1      = sphot
      sphot   = rm*a1
      sphotdt = rm*sphotdt  
      sphotdd = rm*sphotdd + rmdd*a1  
      sphotda = rm*sphotda + rmda*a1  
      sphotdz = rm*sphotdz + rmdz*a1  

      a1      = tfac4*(1.0d0 - tfac3 * qphot)
      a2      = -tfac4*tfac3

      a3      = sphot
      sphot   = a1*a3
      sphotdt = a1*sphotdt + a2*qphotdt*a3
      sphotdd = a1*sphotdd + a2*qphotdd*a3
      sphotda = a1*sphotda + a2*qphotda*a3
      sphotdz = a1*sphotdz + a2*qphotdz*a3

      if (sphot .le. 0.0) then
       sphot   = 0.0d0
       sphotdt = 0.0d0
       sphotdd = 0.0d0
       sphotda = 0.0d0
       sphotdz = 0.0d0
      end if





c..bremsstrahlung neutrino section 
c..for reactions like e- + (z,a) => e- + (z,a) + nu + nubar
c..                   n  + n     => n + n + nu + nubar
c..                   n  + p     => n + p + nu + nubar
c..equation 4.3

      den6   = den * 1.0d-6
      t8     = temp * 1.0d-8
      t812   = sqrt(t8)
      t832   = t8 * t812
      t82    = t8*t8
      t83    = t82*t8 
      t85    = t82*t83
      t86    = t85*t8
      t8m1   = 1.0d0/t8
      t8m2   = t8m1*t8m1
      t8m3   = t8m2*t8m1
      t8m5   = t8m3*t8m2
      t8m6   = t8m5*t8m1


      tfermi = 5.9302d9*(sqrt(1.0d0+1.018d0*(den6*ye)**twoth)-1.0d0)

c.."weak" degenerate electrons only
      if (temp .gt. 0.3d0 * tfermi) then

c..equation 5.3
       dum   = 7.05d6 * t832 + 5.12d4 * t83
       dumdt = (1.5d0*7.05d6*t812 + 3.0d0*5.12d4*t82)*1.0d-8

       z     = 1.0d0/dum
       eta   = rm*z
       etadt = -rm*z*z*dumdt
       etadd = rmdd*z
       etada = rmda*z
       etadz = rmdz*z

       etam1 = 1.0d0/eta
       etam2 = etam1 * etam1
       etam3 = etam2 * etam1


c..equation 5.2
       a0    = 23.5d0 + 6.83d4*t8m2 + 7.81d8*t8m5
       f0    = (-2.0d0*6.83d4*t8m3 - 5.0d0*7.81d8*t8m6)*1.0d-8
       xnum  = 1.0d0/a0

       dum   = 1.0d0 + 1.47d0*etam1 + 3.29d-2*etam2
       z     = -1.47d0*etam2 - 2.0d0*3.29d-2*etam3
       dumdt = z*etadt
       dumdd = z*etadd
       dumda = z*etada
       dumdz = z*etadz

       c00   = 1.26d0 * (1.0d0+etam1)
       z     = -1.26d0*etam2
       c01   = z*etadt
       c02   = z*etadd
       c03   = z*etada
       c04   = z*etadz
       
       z      = 1.0d0/dum
       xden   = c00*z
       xdendt = (c01 - xden*dumdt)*z
       xdendd = (c02 - xden*dumdd)*z
       xdenda = (c03 - xden*dumda)*z
       xdendz = (c04 - xden*dumdz)*z

       fbrem   = xnum + xden
       fbremdt = -xnum*xnum*f0 + xdendt
       fbremdd = xdendd
       fbremda = xdenda
       fbremdz = xdendz


c..equation 5.9
       a0    = 230.0d0 + 6.7d5*t8m2 + 7.66d9*t8m5
       f0    = (-2.0d0*6.7d5*t8m3 - 5.0d0*7.66d9*t8m6)*1.0d-8

       z     = 1.0d0 + rm*1.0d-9 
       dum   = a0*z
       dumdt = f0*z
       z     = a0*1.0d-9
       dumdd = z*rmdd
       dumda = z*rmda
       dumdz = z*rmdz

       xnum   = 1.0d0/dum
       z      = -xnum*xnum
       xnumdt = z*dumdt
       xnumdd = z*dumdd
       xnumda = z*dumda
       xnumdz = z*dumdz

       c00   = 7.75d5*t832 + 247.0d0*t8**(3.85d0)
       dd00  = (1.5d0*7.75d5*t812 + 3.85d0*247.0d0*t8**(2.85d0))*1.0d-8

       c01   = 4.07d0 + 0.0240d0 * t8**(1.4d0)
       dd01  = 1.4d0*0.0240d0*t8**(0.4d0)*1.0d-8

       c02   = 4.59d-5 * t8**(-0.110d0)
       dd02  = -0.11d0*4.59d-5 * t8**(-1.11d0)*1.0d-8

       z     = den**(0.656d0)
       dum   = c00*rmi  + c01  + c02*z 
       dumdt = dd00*rmi + dd01 + dd02*z
       z     = -c00*rmi*rmi
       dumdd = z*rmdd + 0.656d0*c02*den**(-0.454d0)
       dumda = z*rmda 
       dumdz = z*rmdz 

       xden  = 1.0d0/dum
       z      = -xden*xden
       xdendt = z*dumdt
       xdendd = z*dumdd
       xdenda = z*dumda
       xdendz = z*dumdz

       gbrem   = xnum + xden
       gbremdt = xnumdt + xdendt
       gbremdd = xnumdd + xdendd
       gbremda = xnumda + xdenda
       gbremdz = xnumdz + xdendz


c..equation 5.1
       dum    = 0.5738d0*zbar*ye*t86*den
       dumdt  = 0.5738d0*zbar*ye*6.0d0*t85*den*1.0d-8
       dumdd  = 0.5738d0*zbar*ye*t86
       dumda  = -dum*abari
       dumdz  = 0.5738d0*2.0d0*ye*t86*den

       z       = tfac4*fbrem - tfac5*gbrem
       sbrem   = dum * z
       sbremdt = dumdt*z + dum*(tfac4*fbremdt - tfac5*gbremdt) 
       sbremdd = dumdd*z + dum*(tfac4*fbremdd - tfac5*gbremdd) 
       sbremda = dumda*z + dum*(tfac4*fbremda - tfac5*gbremda) 
       sbremdz = dumdz*z + dum*(tfac4*fbremdz - tfac5*gbremdz) 




c..liquid metal with c12 parameters (not too different for other elements)
c..equation 5.18 and 5.16

      else
       u     = fac3 * (log10(den) - 3.0d0)
       a0    = iln10*fac3*deni

c..compute the expensive trig functions of equation 5.21 only once
       cos1 = cos(u)
       cos2 = cos(2.0d0*u)
       cos3 = cos(3.0d0*u)
       cos4 = cos(4.0d0*u)
       cos5 = cos(5.0d0*u)

       sin1 = sin(u)
       sin2 = sin(2.0d0*u)
       sin3 = sin(3.0d0*u)
       sin4 = sin(4.0d0*u)
       sin5 = sin(5.0d0*u)

c..equation 5.21
       fb =  0.5d0 * 0.17946d0  + 0.00945d0*u + 0.34529d0   
     1       - 0.05821d0*cos1 - 0.04969d0*sin1
     2       - 0.01089d0*cos2 - 0.01584d0*sin2
     3       - 0.01147d0*cos3 - 0.00504d0*sin3
     4       - 0.00656d0*cos4 - 0.00281d0*sin4
     5       - 0.00519d0*cos5 

       c00 =  a0*(0.00945d0 
     1       + 0.05821d0*sin1       - 0.04969d0*cos1
     2       + 0.01089d0*sin2*2.0d0 - 0.01584d0*cos2*2.0d0
     3       + 0.01147d0*sin3*3.0d0 - 0.00504d0*cos3*3.0d0
     4       + 0.00656d0*sin4*4.0d0 - 0.00281d0*cos4*4.0d0
     5       + 0.00519d0*sin5*5.0d0) 

      
c..equation 5.22
       ft =  0.5d0 * 0.06781d0 - 0.02342d0*u + 0.24819d0
     1       - 0.00944d0*cos1 - 0.02213d0*sin1
     2       - 0.01289d0*cos2 - 0.01136d0*sin2
     3       - 0.00589d0*cos3 - 0.00467d0*sin3
     4       - 0.00404d0*cos4 - 0.00131d0*sin4
     5       - 0.00330d0*cos5 

       c01 = a0*(-0.02342d0  
     1       + 0.00944d0*sin1       - 0.02213d0*cos1
     2       + 0.01289d0*sin2*2.0d0 - 0.01136d0*cos2*2.0d0
     3       + 0.00589d0*sin3*3.0d0 - 0.00467d0*cos3*3.0d0
     4       + 0.00404d0*sin4*4.0d0 - 0.00131d0*cos4*4.0d0
     5       + 0.00330d0*sin5*5.0d0) 


c..equation 5.23
       gb =  0.5d0 * 0.00766d0 - 0.01259d0*u + 0.07917d0
     1       - 0.00710d0*cos1 + 0.02300d0*sin1
     2       - 0.00028d0*cos2 - 0.01078d0*sin2
     3       + 0.00232d0*cos3 + 0.00118d0*sin3
     4       + 0.00044d0*cos4 - 0.00089d0*sin4
     5       + 0.00158d0*cos5

       c02 = a0*(-0.01259d0
     1       + 0.00710d0*sin1       + 0.02300d0*cos1
     2       + 0.00028d0*sin2*2.0d0 - 0.01078d0*cos2*2.0d0
     3       - 0.00232d0*sin3*3.0d0 + 0.00118d0*cos3*3.0d0
     4       - 0.00044d0*sin4*4.0d0 - 0.00089d0*cos4*4.0d0
     5       - 0.00158d0*sin5*5.0d0)


c..equation 5.24
       gt =  -0.5d0 * 0.00769d0  - 0.00829d0*u + 0.05211d0
     1       + 0.00356d0*cos1 + 0.01052d0*sin1
     2       - 0.00184d0*cos2 - 0.00354d0*sin2
     3       + 0.00146d0*cos3 - 0.00014d0*sin3
     4       + 0.00031d0*cos4 - 0.00018d0*sin4
     5       + 0.00069d0*cos5 

       c03 = a0*(-0.00829d0
     1       - 0.00356d0*sin1       + 0.01052d0*cos1
     2       + 0.00184d0*sin2*2.0d0 - 0.00354d0*cos2*2.0d0
     3       - 0.00146d0*sin3*3.0d0 - 0.00014d0*cos3*3.0d0
     4       - 0.00031d0*sin4*4.0d0 - 0.00018d0*cos4*4.0d0
     5       - 0.00069d0*sin5*5.0d0) 


       dum   = 2.275d-1 * zbar * zbar*t8m1 * (den6*abari)**oneth
       dumdt = -dum*tempi
       dumdd = oneth*dum*deni
       dumda = -oneth*dum*abari
       dumdz = 2.0d0*dum*zbari
     
       gm1   = 1.0d0/dum
       gm2   = gm1*gm1
       gm13  = gm1**oneth
       gm23  = gm13 * gm13
       gm43  = gm13*gm1
       gm53  = gm23*gm1


c..equation 5.25 and 5.26
       v  = -0.05483d0 - 0.01946d0*gm13 + 1.86310d0*gm23 - 0.78873d0*gm1
       a0 = oneth*0.01946d0*gm43 - twoth*1.86310d0*gm53 + 0.78873d0*gm2

       w  = -0.06711d0 + 0.06859d0*gm13 + 1.74360d0*gm23 - 0.74498d0*gm1
       a1 = -oneth*0.06859d0*gm43 - twoth*1.74360d0*gm53 + 0.74498d0*gm2


c..equation 5.19 and 5.20
       fliq   = v*fb + (1.0d0 - v)*ft
       fliqdt = a0*dumdt*(fb - ft)
       fliqdd = a0*dumdd*(fb - ft) + v*c00 + (1.0d0 - v)*c01 
       fliqda = a0*dumda*(fb - ft)
       fliqdz = a0*dumdz*(fb - ft)

       gliq   = w*gb + (1.0d0 - w)*gt
       gliqdt = a1*dumdt*(gb - gt)
       gliqdd = a1*dumdd*(gb - gt) + w*c02 + (1.0d0 - w)*c03
       gliqda = a1*dumda*(gb - gt)
       gliqdz = a1*dumdz*(gb - gt)


c..equation 5.17
       dum    = 0.5738d0*zbar*ye*t86*den
       dumdt  = 0.5738d0*zbar*ye*6.0d0*t85*den*1.0d-8
       dumdd  = 0.5738d0*zbar*ye*t86
       dumda  = -dum*abari
       dumdz  = 0.5738d0*2.0d0*ye*t86*den

       z       = tfac4*fliq - tfac5*gliq
       sbrem   = dum * z
       sbremdt = dumdt*z + dum*(tfac4*fliqdt - tfac5*gliqdt) 
       sbremdd = dumdd*z + dum*(tfac4*fliqdd - tfac5*gliqdd) 
       sbremda = dumda*z + dum*(tfac4*fliqda - tfac5*gliqda) 
       sbremdz = dumdz*z + dum*(tfac4*fliqdz - tfac5*gliqdz) 

      end if




c..recombination neutrino section
c..for reactions like e- (continuum) => e- (bound) + nu_e + nubar_e
c..equation 6.11 solved for nu
      xnum   = 1.10520d8 * den * ye /(temp*sqrt(temp))
      xnumdt = -1.50d0*xnum*tempi
      xnumdd = xnum*deni
      xnumda = -xnum*abari
      xnumdz = xnum*zbari

c..the chemical potential
      nu   = ifermi12(xnum)

c..a0 is d(nu)/d(xnum)
      a0 = 1.0d0/(0.5d0*zfermim12(nu))
      nudt = a0*xnumdt
      nudd = a0*xnumdd
      nuda = a0*xnumda
      nudz = a0*xnumdz

      nu2  = nu * nu
      nu3  = nu2 * nu

c..table 12
      if (nu .ge. -20.0  .and. nu .lt. 0.0) then
       a1 = 1.51d-2
       a2 = 2.42d-1
       a3 = 1.21d0
       b  = 3.71d-2
       c  = 9.06e-1
       d  = 9.28d-1
       f1 = 0.0d0
       f2 = 0.0d0
       f3 = 0.0d0
      else if (nu .ge. 0.0  .and. nu .le. 10.0) then
       a1 = 1.23d-2
       a2 = 2.66d-1
       a3 = 1.30d0
       b  = 1.17d-1
       c  = 8.97e-1
       d  = 1.77d-1
       f1 = -1.20d-2
       f2 = 2.29d-2
       f3 = -1.04d-3
      end if


c..equation 6.7, 6.13 and 6.14
      if (nu .ge. -20.0  .and.  nu .le. 10.0) then

       zeta   = 1.579d5*zbar*zbar*tempi
       zetadt = -zeta*tempi
       zetadd = 0.0d0
       zetada = 0.0d0
       zetadz = 2.0d0*zeta*zbari

       c00    = 1.0d0/(1.0d0 + f1*nu + f2*nu2 + f3*nu3)  
       c01    = f1 + f2*2.0d0*nu + f3*3.0d0*nu2
       dum    = zeta*c00
       dumdt  = zetadt*c00 + zeta*c01*nudt
       dumdd  = zeta*c01*nudd
       dumda  = zeta*c01*nuda
       dumdz  = zetadz*c00 + zeta*c01*nudz

     
       z      = 1.0d0/dum
       dd00   = dum**(-2.25) 
       dd01   = dum**(-4.55)
       c00    = a1*z + a2*dd00 + a3*dd01
       c01    = -(a1*z + 2.25*a2*dd00 + 4.55*a3*dd01)*z
    

       z      = exp(c*nu)  
       dd00   = b*z*(1.0d0 + d*dum)        
       gum    = 1.0d0 + dd00
       gumdt  = dd00*c*nudt + b*z*d*dumdt  
       gumdd  = dd00*c*nudd + b*z*d*dumdd  
       gumda  = dd00*c*nuda + b*z*d*dumda  
       gumdz  = dd00*c*nudz + b*z*d*dumdz  


       z   = exp(nu)  
       a1  = 1.0d0/gum

       bigj   = c00 * z * a1
       bigjdt = c01*dumdt*z*a1 + c00*z*nudt*a1 - c00*z*a1*a1 * gumdt
       bigjdd = c01*dumdd*z*a1 + c00*z*nudd*a1 - c00*z*a1*a1 * gumdd
       bigjda = c01*dumda*z*a1 + c00*z*nuda*a1 - c00*z*a1*a1 * gumda
       bigjdz = c01*dumdz*z*a1 + c00*z*nudz*a1 - c00*z*a1*a1 * gumdz


c..equation 6.5
       z     = exp(zeta + nu)
       dum   = 1.0d0 + z
       a1    = 1.0d0/dum
       a2    = 1.0d0/bigj

       sreco   = tfac6 * 2.649d-18 * ye * zbar**13 * den * bigj*a1
       srecodt = sreco*(bigjdt*a2 - z*(zetadt + nudt)*a1)
       srecodd = sreco*(1.0d0*deni + bigjdd*a2 - z*(zetadd + nudd)*a1)
       srecoda = sreco*(-1.0d0*abari + bigjda*a2 - z*(zetada+nuda)*a1)
       srecodz = sreco*(14.0d0*zbari + bigjdz*a2 - z*(zetadz+nudz)*a1)

      end if 


c..convert from erg/cm^3/s to erg/g/s 
c..comment these out to duplicate the itoh et al plots

      spair   = spair*deni
      spairdt = spairdt*deni
      spairdd = spairdd*deni - spair*deni
      spairda = spairda*deni
      spairdz = spairdz*deni  

      splas   = splas*deni
      splasdt = splasdt*deni
      splasdd = splasdd*deni - splas*deni
      splasda = splasda*deni
      splasdz = splasdz*deni  

      sphot   = sphot*deni
      sphotdt = sphotdt*deni
      sphotdd = sphotdd*deni - sphot*deni
      sphotda = sphotda*deni
      sphotdz = sphotdz*deni  

      sbrem   = sbrem*deni
      sbremdt = sbremdt*deni
      sbremdd = sbremdd*deni - sbrem*deni
      sbremda = sbremda*deni
      sbremdz = sbremdz*deni  

      sreco   = sreco*deni
      srecodt = srecodt*deni
      srecodd = srecodd*deni - sreco*deni
      srecoda = srecoda*deni
      srecodz = srecodz*deni  


c..the total neutrino loss rate
      snu    =  splas + spair + sphot + sbrem + sreco
      dsnudt =  splasdt + spairdt + sphotdt + sbremdt + srecodt  
      dsnudd =  splasdd + spairdd + sphotdd + sbremdd + srecodd 
      dsnuda =  splasda + spairda + sphotda + sbremda + srecoda 
      dsnudz =  splasdz + spairdz + sphotdz + sbremdz + srecodz 

      return
      end






      double precision function ifermi12(f)
      include 'implno.dek'

c..this routine applies a rational function expansion to get the inverse
c..fermi-dirac integral of order 1/2 when it is equal to f.
c..maximum error is 4.19d-9.   reference: antia apjs 84,101 1993

c..declare
      integer          i,m1,k1,m2,k2
      double precision f,an,a1(12),b1(12),a2(12),b2(12),rn,den,ff,
     1                 z,drn


c..load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /0.5d0, 4, 3, 6, 5/
      data  (a1(i),i=1,5)/ 1.999266880833d4,   5.702479099336d3,
     1     6.610132843877d2,   3.818838129486d1,
     2     1.0d0/
      data  (b1(i),i=1,4)/ 1.771804140488d4,  -2.014785161019d3,
     1     9.130355392717d1,  -1.670718177489d0/
      data  (a2(i),i=1,7)/-1.277060388085d-2,  7.187946804945d-2, 
     1                    -4.262314235106d-1,  4.997559426872d-1,
     2                    -1.285579118012d0,  -3.930805454272d-1,
     3     1.0d0/
      data  (b2(i),i=1,6)/-9.745794806288d-3,  5.485432756838d-2,
     1                    -3.299466243260d-1,  4.077841975923d-1,
     2                    -1.145531476975d0,  -6.067091689181d-2/


      if (f .lt. 4.0d0) then
       rn  = f + a1(m1)
       do i=m1-1,1,-1
        rn  = rn*f + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*f + b1(i)
       enddo
       ifermi12 = log(f * rn/den)

      else
       ff = 1.0d0/f**(1.0d0/(1.0d0 + an))
       rn = ff + a2(m2)
       do i=m2-1,1,-1
        rn = rn*ff + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*ff + b2(i)
       enddo
       ifermi12 = rn/(den*ff)
      end if
      return
      end






      double precision function zfermim12(x)
      include 'implno.dek'

c..this routine applies a rational function expansion to get the fermi-dirac
c..integral of order -1/2 evaluated at x. maximum error is 1.23d-12.
c..reference: antia apjs 84,101 1993

c..declare
      integer          i,m1,k1,m2,k2
      double precision x,an,a1(12),b1(12),a2(12),b2(12),rn,den,xx

c..load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /-0.5d0, 7, 7, 11, 11/
      data  (a1(i),i=1,8)/ 1.71446374704454d7,    3.88148302324068d7,
     1                     3.16743385304962d7,    1.14587609192151d7,
     2                     1.83696370756153d6,    1.14980998186874d5,
     3                     1.98276889924768d3,    1.0d0/
      data  (b1(i),i=1,8)/ 9.67282587452899d6,    2.87386436731785d7,
     1                     3.26070130734158d7,    1.77657027846367d7,
     2                     4.81648022267831d6,    6.13709569333207d5,
     3                     3.13595854332114d4,    4.35061725080755d2/
      data (a2(i),i=1,12)/-4.46620341924942d-15, -1.58654991146236d-12,
     1                    -4.44467627042232d-10, -6.84738791621745d-8,
     2                    -6.64932238528105d-6,  -3.69976170193942d-4,
     3                    -1.12295393687006d-2,  -1.60926102124442d-1,
     4                    -8.52408612877447d-1,  -7.45519953763928d-1,
     5                     2.98435207466372d0,    1.0d0/
      data (b2(i),i=1,12)/-2.23310170962369d-15, -7.94193282071464d-13,
     1                    -2.22564376956228d-10, -3.43299431079845d-8,
     2                    -3.33919612678907d-6,  -1.86432212187088d-4,
     3                    -5.69764436880529d-3,  -8.34904593067194d-2,
     4                    -4.78770844009440d-1,  -4.99759250374148d-1,
     5                     1.86795964993052d0,    4.16485970495288d-1/


      if (x .lt. 2.0d0) then
       xx = exp(x)
       rn = xx + a1(m1)
       do i=m1-1,1,-1
        rn = rn*xx + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*xx + b1(i)
       enddo
       zfermim12 = xx * rn/den
c..
      else
       xx = 1.0d0/(x*x)
       rn = xx + a2(m2)
       do i=m2-1,1,-1
        rn = rn*xx + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*xx + b2(i)
       enddo
       zfermim12 = sqrt(x)*rn/den
      end if
      return
      end






      subroutine mazurek(btemp,bden,y56,ye,rn56ec,sn56ec) 
      include 'implno.dek'

c..this routine evaluates mazurel's 1973 fits for the ni56 electron 
c..capture rate rn56ec and neutrino loss rate sn56ec 

c..input: 
c..y56 = nickel56 molar abundance
c..ye  = electron to baryon number, zbar/abar

c..output:
c..rn56ec = ni56 electron capture rate
c..sn56ec = ni56 neutrino loss rate

c..declare 
      integer          ifirst,jp,kp,jr,jd,ii,ik,ij,j,k 
      double precision btemp,bden,y56,ye,rn56ec,sn56ec,
     1                 rnt(2),rne(2,7),datn(2,6,7), 
     2                 tv(7),rv(6),rfdm(4),rfd0(4),rfd1(4),rfd2(4), 
     3                 tfdm(5),tfd0(5),tfd1(5),tfd2(5), 
     4                 t9,r,rfm,rf0,rf1,rf2,dfacm,dfac0,dfac1,dfac2, 
     5                 tfm,tf0,tf1,tf2,tfacm,tfac0,tfac1,tfac2

c..initialize 
      data  rv /6.0, 7.0, 8.0, 9.0, 10.0, 11.0/ 
      data  tv /2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0/ 
      data (((datn(ii,ik,ij),ik=1,6),ij=1,7),ii=1,1) / 
     1    -3.98, -2.84, -1.41,  0.20,  1.89,  3.63, 
     2    -3.45, -2.62, -1.32,  0.22,  1.89,  3.63, 
     3    -2.68, -2.30, -1.19,  0.27,  1.91,  3.62, 
     4    -2.04, -1.87, -1.01,  0.34,  1.94,  3.62, 
     5    -1.50, -1.41, -0.80,  0.45,  1.99,  3.60, 
     6    -1.00, -0.95, -0.54,  0.60,  2.06,  3.58, 
     7    -0.52, -0.49, -0.21,  0.79,  2.15,  3.55 / 
      data (((datn(ii,ik,ij),ik=1,6),ij=1,7),ii=2,2) / 
     1    -3.68, -2.45, -0.80,  1.12,  3.13,  5.19, 
     2    -2.91, -2.05, -0.64,  1.16,  3.14,  5.18, 
     3    -1.95, -1.57, -0.40,  1.24,  3.16,  5.18, 
     4    -1.16, -0.99, -0.11,  1.37,  3.20,  5.18, 
     5    -0.48, -0.40,  0.22,  1.54,  3.28,  5.16, 
     6     0.14,  0.19,  0.61,  1.78,  3.38,  5.14, 
     7     0.75,  0.78,  1.06,  2.07,  3.51,  5.11 / 
      data  ifirst /0/ 

c..first time; calculate the cubic interp parameters for ni56 electron capture 
      if (ifirst .eq. 0) then 
       ifirst = 1 
       do k=2,4 
        rfdm(k)=1./((rv(k-1)-rv(k))*(rv(k-1)-rv(k+1))*(rv(k-1)-rv(k+2))) 
        rfd0(k)=1./((rv(k)-rv(k-1))*(rv(k)-rv(k+1))*(rv(k)-rv(k+2))) 
        rfd1(k)=1./((rv(k+1)-rv(k-1))*(rv(k+1)-rv(k))*(rv(k+1)-rv(k+2))) 
        rfd2(k)=1./((rv(k+2)-rv(k-1))*(rv(k+2)-rv(k))*(rv(k+2)-rv(k+1))) 
       enddo
       do j=2,5 
        tfdm(j)=1./((tv(j-1)-tv(j))*(tv(j-1)-tv(j+1))*(tv(j-1)-tv(j+2))) 
        tfd0(j)=1./((tv(j)-tv(j-1))*(tv(j)-tv(j+1))*(tv(j)-tv(j+2))) 
        tfd1(j)=1./((tv(j+1)-tv(j-1))*(tv(j+1)-tv(j))*(tv(j+1)-tv(j+2))) 
        tfd2(j)=1./((tv(j+2)-tv(j-1))*(tv(j+2)-tv(j))*(tv(j+2)-tv(j+1))) 
       enddo
      end if 

c..calculate ni56 electron capture and neutrino loss rates 
      rn56ec = 0.0 
      sn56ec = 0.0 
      if ( (btemp .lt. 2.0e9) .or. (bden*ye .lt. 1.0e6)) return 
      t9    = max(btemp,1.4d10) * 1.0d-9 
      r     = max(6.0d0,min(11.0d0,log10(bden*ye))) 
      jp    = min(max(2,int(0.5d0*t9)),5) 
      kp    = min(max(2,int(r)-5),4) 
      rfm   = r - rv(kp-1) 
      rf0   = r - rv(kp) 
      rf1   = r - rv(kp+1) 
      rf2   = r - rv(kp+2) 
      dfacm = rf0*rf1*rf2*rfdm(kp) 
      dfac0 = rfm*rf1*rf2*rfd0(kp) 
      dfac1 = rfm*rf0*rf2*rfd1(kp) 
      dfac2 = rfm*rf0*rf1*rfd2(kp) 
      tfm   = t9 - tv(jp-1) 
      tf0   = t9 - tv(jp) 
      tf1   = t9 - tv(jp+1) 
      tf2   = t9 - tv(jp+2) 
      tfacm = tf0*tf1*tf2*tfdm(jp) 
      tfac0 = tfm*tf1*tf2*tfd0(jp) 
      tfac1 = tfm*tf0*tf2*tfd1(jp) 
      tfac2 = tfm*tf0*tf1*tfd2(jp) 

c..evaluate the spline fits
      do jr = 1,2 
       do jd = jp-1,jp+2 
        rne(jr,jd) =   dfacm*datn(jr,kp-1,jd) + dfac0*datn(jr,kp,jd) 
     1               + dfac1*datn(jr,kp+1,jd) + dfac2*datn(jr,kp+2,jd) 
       enddo
       rnt(jr) =  tfacm*rne(jr,jp-1) + tfac0*rne(jr,jp) 
     1          + tfac1*rne(jr,jp+1) + tfac2*rne(jr,jp+2) 
      enddo

c..set the output
      rn56ec = 10.0d0**rnt(1) 
      sn56ec = 6.022548d+23 * 8.18683d-7 * y56 * 10.0d0**rnt(2) 
      return 
      end 








      subroutine ecapnuc(etakep,temp,rpen,rnep,spen,snep)
      include 'implno.dek'

c..given the electron degeneracy parameter etakep (chemical potential
c..without the electron's rest mass divided by kt) and the temperature temp,
c..this routine calculates rates for 
c..electron capture on protons rpen (captures/sec/proton),
c..positron capture on neutrons rnep (captures/sec/neutron), 
c..and their associated neutrino energy loss rates 
c..spen (ergs/sec/proton) and snep (ergs/sec/neutron)

c..declare
      double precision etakep,temp,rpen,rnep,spen,snep

      integer          iflag
      double precision t9,t5,qn,etaef,etael,zetan,eta,etael2,
     1                 etael3,etael4,f1l,f2l,f3l,f4l,f5l,f1g,
     2                 f2g,f3g,f4g,f5g,exmeta,eta2,eta3,eta4,
     3                 fac0,fac1,fac2,fac3,rie1,rie2,facv0,facv1,
     4                 facv2,facv3,facv4,rjv1,rjv2,spenc,snepc,
     5                 pi2,exeta,zetan2,f0,etael5,
     6                 qn1,ft,twoln,cmk5,cmk6,bk,pi,qn2,c2me,
     7                 xmp,xmn,qndeca,tmean
      parameter        (qn1    = -2.0716446d-06,
     1                  ft     = 1083.9269d0,
     2                  twoln  = 0.6931472d0,
     3                  cmk5   = 1.3635675d-49,
     4                  cmk6   = 2.2993864d-59,
     5                  bk     = 1.38062e-16,
     6                  pi     = 3.1415927d0,
     7                  pi2    = pi * pi,
     8                  qn2    = 2.0716446d-06,
     9                  c2me   = 8.1872665d-07,
     &                  xmp    = 1.6726485d-24,
     1                  xmn    = 1.6749543d-24,
     2                  qndeca = 1.2533036d-06,
     3                  tmean  = 886.7d0)
c     3                  tmean  = 935.14d0)
      


c..tmean and qndeca are the mean lifetime and decay energy of the neutron
c..xmp,xnp are masses of the p and n in grams.
c..c2me is the constant used to convert the neutrino energy
c..loss rate from mec2/s (as in the paper) to ergs/particle/sec.

c..initialize
      rpen  = 0.0d0
      rnep  = 0.0d0
      spen  = 0.0d0
      snep  = 0.0d0
      t9    = temp * 1.0d-9
      iflag = 0
      qn    = qn1


c..chemical potential including the electron rest mass
      etaef = etakep + c2me/bk/temp


c..iflag=1 is for electrons,  iflag=2 is for positrons
502   iflag = iflag + 1
      if (iflag.eq.1) etael = qn2/bk/temp
      if (iflag.eq.2) etael = c2me/bk/temp
      if (iflag.eq.2) etaef = -etaef

      t5    = temp*temp*temp*temp*temp
      zetan = qn/bk/temp
      eta   = etaef - etael

c..protect from overflowing with large eta values
      if (eta .le. 6.8e+02) then
       exeta = exp(eta)
      else 
       exeta = 0.0d0
      end if
      etael2 = etael*etael
      etael3 = etael2*etael
      etael4 = etael3*etael
      etael5 = etael4*etael
      zetan2 = zetan*zetan
      if (eta .le. 6.8e+02) then
       f0 = log(1.0d0 + exeta)
      else
       f0 = eta
      end if

c..if eta le. 0., the following fermi integrals apply
      f1l = exeta
      f2l = 2.0d0   * f1l
      f3l = 6.0d0   * f1l
      f4l = 24.0d0  * f1l
      f5l = 120.0d0 * f1l

c..if eta gt. 0., the following fermi integrals apply:
      f1g = 0.0d0
      f2g = 0.0d0
      f3g = 0.0d0
      f4g = 0.0d0
      f5g = 0.0d0
      if (eta .gt. 0.0) then
       exmeta = dexp(-eta)
       eta2   = eta*eta
       eta3   = eta2*eta
       eta4   = eta3*eta
       f1g = 0.5d0*eta2 + 2.0d0 - exmeta
       f2g = eta3/3.0d0 + 4.0d0*eta + 2.0d0*exmeta
       f3g = 0.25d0*eta4 + 0.5d0*pi2*eta2 + 12.0d0 - 6.0d0*exmeta
       f4g = 0.2d0*eta4*eta + 2.0d0*pi2/3.0d0*eta3 + 48.0d0*eta
     1       + 24.0d0*exmeta
       f5g = eta4*eta2/6.0d0 + 5.0d0/6.0d0*pi2*eta4 
     2       + 7.0d0/6.0d0*pi2*eta2  + 240.0d0 -120.d0*exmeta
       end if

c..factors which are multiplied by the fermi integrals
      fac3 = 2.0d0*zetan + 4.0d0*etael
      fac2 = 6.0d0*etael2 + 6.0d0*etael*zetan + zetan2
      fac1 = 4.0d0*etael3 + 6.0d0*etael2*zetan + 2.0d0*etael*zetan2
      fac0 = etael4 + 2.0d0*zetan*etael3 + etael2*zetan2

c..electron capture rates onto protons with no blocking
      rie1 = f4l + fac3*f3l + fac2*f2l + fac1*f1l + fac0*f0
      rie2 = f4g + fac3*f3g + fac2*f2g + fac1*f1g + fac0*f0

c..neutrino emission rate for electron capture:
      facv4 = 5.0d0*etael + 3.0d0*zetan
      facv3 = 10.0d0*etael2 + 12.0d0*etael*zetan + 3.0d0*zetan2
      facv2 = 10.0d0*etael3 + 18.0d0*etael2*zetan
     1        + 9.0d0*etael*zetan2 + zetan2*zetan
      facv1 = 5.0d0*etael4 + 12.0d0*etael3*zetan 
     1        + 9.0d0*etael2*zetan2 + 2.0d0*etael*zetan2*zetan
      facv0 = etael5 + 3.0d0*etael4*zetan
     1        + 3.0d0*etael3*zetan2 + etael2*zetan2*zetan
      rjv1  = f5l + facv4*f4l + facv3*f3l
     1        + facv2*f2l + facv1*f1l + facv0*f0
      rjv2  = f5g + facv4*f4g + facv3*f3g
     1        + facv2*f2g + facv1*f1g + facv0*f0

c..for electrons capture onto protons
      if (iflag.eq.2) go to 503
      if (eta.gt.0.) go to 505
      rpen  = twoln*cmk5*t5*rie1/ft
      spen  = twoln*cmk6*t5*temp*rjv1/ft
      spenc = twoln*cmk6*t5*temp*rjv1/ft*c2me
      go to 504
505   rpen = twoln*cmk5*t5*rie2/ft
      spen = twoln*cmk6*t5*temp*rjv2/ft
      spenc = twoln*cmk6*t5*temp*rjv2/ft*c2me
504   continue
      qn = qn2
      go to 502

c..for positrons capture onto neutrons
503   if (eta.gt.0.) go to 507
      rnep  = twoln*cmk5*t5*rie1/ft
      snep  = twoln*cmk6*t5*temp*rjv1/ft
      snepc = twoln*cmk6*t5*temp*rjv1/ft*c2me
c      if (rho.lt.1.0e+06) snep=snep+qndeca*xn(9)/xmn/tmean
      go to 506
507   rnep  = twoln*cmk5*t5*rie2/ft
      snep  = twoln*cmk6*t5*temp*rjv2/ft
      snepc = twoln*cmk6*t5*temp*rjv2/ft*c2me
c      if (rho.lt.1.0e+06) snep=snep+qndeca*xn(9)/xmn/tmean
506   continue
      return
      end







c..this file contains routines that sort, search and select parts of arrays: 
c.. 
c..index and rank makers: 
c..routine indexx constructs a sort index for a real array





      subroutine indexx(n,arr,indx) 
      include 'implno.dek' 
c.. 
c..indexes an array arr(1:n). that is it outputs the array indx(1:n) such 
c..that arr(indx(j)) is in ascending order for j=1...n. the input quantities 
c..are not changed. 
c.. 
c..declare 
      integer          n,indx(n),m,nstack 
      parameter        (m=7, nstack = 50) 
      integer          i,indxt,ir,itemp,j,jstack,k,l,istack(nstack) 
      double precision arr(n),a 
c.. 
c..initialize 
      do 11 j=1,n 
       indx(j) = j 
11    continue 
      jstack = 0 
      l      = 1 
      ir     = n 
c.. 
c..insertion sort when subbarray small enough 
1     if (ir - l .lt. m) then 
       do 13 j=l+1,ir 
        indxt = indx(j) 
        a     = arr(indxt) 
        do 12 i=j-1,l,-1 
         if (arr(indx(i)) .le. a) go to 2 
         indx(i+1) = indx(i) 
12      continue 
        i = l - 1 
2       indx(i+1) = indxt 
13     continue 
c.. 
c..pop stack and begin a new round of partitioning 
       if (jstack .eq. 0) return 
       ir     = istack(jstack) 
       l      = istack(jstack-1) 
       jstack = jstack - 2 
c.. 
c..choose median of left, center and right elements as partitioning element 
c..also rearrange so that a(l+1) < a(l) < a(ir) 
      else 
       k         = (l + ir)/2 
       itemp     = indx(k) 
       indx(k)   = indx(l+1) 
       indx(l+1) = itemp 
 
       if (arr(indx(l)) .gt. arr(indx(ir))) then 
        itemp    = indx(l) 
        indx(l)  = indx(ir) 
        indx(ir) = itemp 
       end if 
 
 
       if(arr(indx(l+1)).gt.arr(indx(ir)))then 
        itemp=indx(l+1) 
        indx(l+1)=indx(ir) 
        indx(ir)=itemp 
       endif 
       if(arr(indx(l)).gt.arr(indx(l+1)))then 
        itemp=indx(l) 
        indx(l)=indx(l+1) 
        indx(l+1)=itemp 
       endif 
 
c.. 
c..initialize pointers for partitioning 
       i     = l + 1 
       j     = ir 
       indxt = indx(l+1) 
       a     = arr(indxt) 
3      continue 
       i = i + 1 
       if (arr(indx(i)) .lt. a) go to 3 
4      continue 
       j = j - 1 
       if (arr(indx(j)) .gt. a) go to 4 
       if (j .lt. i) go to 5 
       itemp   = indx(i) 
       indx(i) = indx(j) 
       indx(j) = itemp 
       go to 3 
c.. 
5      indx(l+1) = indx(j) 
       indx(j)   = indxt 
       jstack    = jstack + 2 
c.. 
c..push pointers to larger subarray on stack 
       if (jstack .gt. nstack) stop 'jstack > nstack in routine indexx' 
       if (ir - i + 1  .ge.  j - l) then 
        istack(jstack)   = ir 
        istack(jstack-1) = i 
        ir               = j - 1  
       else 
        istack(jstack)   = j-1 
        istack(jstack-1) = l 
        l                = i 
       end if 
      end if 
      go to 1 
      end 







c..some system and glue utility routines


      integer function lenstr(string,istrln)
      include 'implno.dek'
c..
c..lenstr returns the non blank length length of the string.
c..
c..declare
      integer       istrln,i
      character*(*) string


      lenstr=0
      do i=istrln,1,-1
       if (string(i:i).ne. ' ') then
        if (ichar(string(i:i)).ne. 0 )then
         lenstr=i
         goto  20
        end if
       end if
      enddo
20    return
      end




      subroutine sqeeze(line)
      include 'implno.dek'
c..
c..this routine takes line and removes all blanks, such as
c..those from writing to string with fortran format statements
c..
c..declare
      character*(*)  line
      character*1    achar
      integer        l,n,k,lend,lsiz,lenstr


c..find the end of the line
      lsiz = len(line)
      lend = lenstr(line,lsiz)
      n    = 0
      l    = 0

c..do the compression in place
10    continue
      l = l + 1
      achar = line(l:l)
      if (achar .eq. ' ') goto 10
      n = n + 1
      line(n:n) = achar
      if (l .lt. lend) goto 10

c..blank the rest of the line
      do k=n+1,lsiz
       line(k:k) = ' '
      enddo
      return
      end






      subroutine azbar(xmass,aion,zion,ionmax,
     1                 ymass,abar,zbar)
      include 'implno.dek'

c..this routine calculates composition variables for an eos routine

c..input:
c..mass fractions     = xmass(1:ionmax)
c..number of nucleons = aion(1:ionmax)
c..charge of nucleus  = zion(1:ionmax)
c..number of isotopes = ionmax
c..
c..output:
c..molar abundances        = ymass(1:ionmax), 
c..mean number of nucleons = abar
c..mean nucleon charge     = zbar

c..declare
      integer          i,ionmax
      double precision xmass(ionmax),aion(ionmax),zion(ionmax),
     1                 ymass(ionmax),abar,zbar,zbarxx,ytot1

      zbarxx  = 0.0d0
      ytot1   = 0.0d0
      do i=1,ionmax
       ymass(i) = xmass(i)/aion(i)
       ytot1    = ytot1 + ymass(i)
       zbarxx   = zbarxx + zion(i) * ymass(i)
      enddo
      abar   = 1.0d0/ytot1
      zbar   = zbarxx * abar
      return
      end












      subroutine helmeos
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'


c..given a temperature temp [K], density den [g/cm**3], and a composition 
c..characterized by abar and zbar, this routine returns most of the other 
c..thermodynamic quantities. of prime interest is the pressure [erg/cm**3], 
c..specific thermal energy [erg/gr], the entropy [erg/g/K], along with 
c..their derivatives with respect to temperature, density, abar, and zbar.
c..other quantites such the normalized chemical potential eta (plus its
c..derivatives), number density of electrons and positron pair (along 
c..with their derivatives), adiabatic indices, specific heats, and 
c..relativistically correct sound speed are also returned.
c..
c..this routine assumes planckian photons, an ideal gas of ions,
c..and an electron-positron gas with an arbitrary degree of relativity
c..and degeneracy. interpolation in a table of the helmholtz free energy
c..is used to return the electron-positron thermodynamic quantities.
c..all other derivatives are analytic.
c..
c..references: cox & giuli chapter 24 ; timmes & swesty apj 1999


c..declare
      integer          i,j,k
      double precision temp,den,abar,zbar,ytot1,ye,
     1                 x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida,
     2                 dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt,
     3                 dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt,
     4                 deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt,
     5                 dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion,
     6                 sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd,
     7                 dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp,
     8                 gam1,gam2,gam3,chit,chid,nabad,sound,etaele,
     9                 detadt,detadd,xnefer,dxnedt,dxnedd,s,
     &                 sioncon,forth,forpi,kergavo,ikavo,asoli3,light2

      parameter        (sioncon = (2.0d0 * pi * amu * kerg)/(h*h),
     1                  forth   = 4.0d0/3.0d0,
     2                  forpi   = 4.0d0 * pi,
     3                  kergavo = kerg * avo, 
     4                  ikavo   = 1.0d0/kergavo,
     5                  asoli3  = asol/3.0d0,
     6                  light2  = clight * clight)

c..for the abar derivatives
      double precision dpradda,deradda,dsradda,
     1                 dpionda,deionda,dsionda,
     2                 dpepda,deepda,dsepda,
     3                 dpresda,denerda,dentrda,
     4                 detada,dxneda


c..for the zbar derivatives
      double precision dpraddz,deraddz,dsraddz,
     1                 dpiondz,deiondz,dsiondz,
     2                 dpepdz,deepdz,dsepdz,
     3                 dpresdz,denerdz,dentrdz,
     4                 detadz,dxnedz

c..for the multipliers
      double precision radmult,ionmult,elemult,coulmult


c..for the tables, in general
      character*132    string,word
      logical          ibhere
      integer          imax,jmax

c..normal table
      parameter        (imax = 211, jmax = 71)

c..bigger table 
c      parameter        (imax = 251, jmax = 91)

      double precision d(imax),t(jmax)



c..for the helmholtz free energy tables
      double precision f(imax,jmax),fd(imax,jmax),
     1                 ft(imax,jmax),fdd(imax,jmax),ftt(imax,jmax),
     2                 fdt(imax,jmax),fddt(imax,jmax),fdtt(imax,jmax),
     3                 fddtt(imax,jmax)


c..for the pressure derivative with density ables
      double precision dpdf(imax,jmax),dpdfd(imax,jmax),
     1                 dpdft(imax,jmax),dpdfdd(imax,jmax),
     2                 dpdftt(imax,jmax),dpdfdt(imax,jmax)


c..for chemical potential tables
      double precision ef(imax,jmax),efd(imax,jmax),
     1                 eft(imax,jmax),efdd(imax,jmax),eftt(imax,jmax),
     2                 efdt(imax,jmax)


c..for the number density tables
      double precision xf(imax,jmax),xfd(imax,jmax),
     1                 xft(imax,jmax),xfdd(imax,jmax),xftt(imax,jmax),
     2                 xfdt(imax,jmax)


c..for the interpolations
      integer          iat,jat
      double precision tlo,thi,tstp,tstpi,dlo,dhi,dstp,dstpi,
     1                 tsav,dsav,free,df_d,df_t,df_dd,df_tt,df_dt
      double precision dth,dt2,dti,dt2i,dd,dd2,ddi,dd2i,xt,xd,mxt,mxd,
     1                 si0t,si1t,si2t,si0mt,si1mt,si2mt,
     2                 si0d,si1d,si2d,si0md,si1md,si2md,
     3                 dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt,
     4                 dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md,
     5                 ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt,
     6                 ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md,
     7                 z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2,
     8                 dpsi2,ddpsi2,din,h5,fi(36),
     9                 xpsi0,xdpsi0,xpsi1,xdpsi1,h3,
     1                 w0t,w1t,w2t,w0mt,w1mt,w2mt,
     2                 w0d,w1d,w2d,w0md,w1md,w2md

c..for storing the differences
      double precision dt_sav(jmax),dt2_sav(jmax),
     1                 dti_sav(jmax),dt2i_sav(jmax),
     2                 dd_sav(imax),dd2_sav(imax),
     3                 ddi_sav(imax),dd2i_sav(imax)

c..for partial ionization
      double precision denion,dydz,dxdt,r,drdz,p,dpdt,dpdz,
     1                 fx,dfxdd,dfxdt,dfxdz,w,dsdt,dsdz,
     1                 zeff,dzeffdd,dyedd,pcon1,pcon2
      parameter        (pcon1 = hion * ev2erg,
     1                  pcon2 = (h*h)/(2.0d0 * pi * me) )


c..for the uniform background coulomb correction
      double precision dsdd,dsda,lami,inv_lami,lamida,lamidd,
     1                 plasg,plasgdd,plasgdt,plasgda,plasgdz,
     3                 ecoul,decouldd,decouldt,decoulda,decouldz,
     4                 pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz,
     5                 scoul,dscouldd,dscouldt,dscoulda,dscouldz,
     6                 a1,b1,c1,d1,e1,a2,b2,c2,third,esqu
      parameter        (a1    = -0.898004d0, 
     1                  b1    =  0.96786d0, 
     2                  c1    =  0.220703d0, 
     3                  d1    = -0.86097d0,
     4                  e1    =  2.5269d0, 
     5                  a2    =  0.29561d0, 
     6                  b2    =  1.9885d0,    
     7                  c2    =  0.288675d0,
     8                  third =  1.0d0/3.0d0,
     9                  esqu  =  qe * qe)

c..for initialization
      integer          ifirst
      data             ifirst/0/ 


c..quintic hermite polynomial statement functions
c..psi0 and its derivatives
      psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
      dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
      ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)


c..psi1 and its derivatives
      psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
      dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
      ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)


c..psi2  and its derivatives
      psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
      dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
      ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)


c..biquintic hermite polynomial statement function
      h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)=
     1       fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t
     2     + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt
     3     + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t
     4     + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt
     5     + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t
     6     + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt
     7     + fi(13) *w1d*w0t   + fi(14) *w1md*w0t
     8     + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt
     9     + fi(17) *w2d*w0t   + fi(18) *w2md*w0t
     &     + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt
     1     + fi(21) *w1d*w1t   + fi(22) *w1md*w1t
     2     + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt
     3     + fi(25) *w2d*w1t   + fi(26) *w2md*w1t
     4     + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt
     5     + fi(29) *w1d*w2t   + fi(30) *w1md*w2t
     6     + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt
     7     + fi(33) *w2d*w2t   + fi(34) *w2md*w2t
     8     + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt



c..cubic hermite polynomial statement functions
c..psi0 & derivatives
      xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
      xdpsi0(z) = z * (6.0d0*z - 6.0d0)


c..psi1 & derivatives
      xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
      xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0


c..bicubic hermite polynomial statement function
      h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = 
     1       fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t 
     2     + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt
     3     + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t 
     4     + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt
     5     + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t 
     6     + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt
     7     + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t 
     8     + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt



c..popular format statements
01    format(1x,5(a,1pe11.3))
02    format(1x,a,1p4e16.8)
03    format(1x,4(a,1pe11.3))
04    format(1x,4(a,i4))


c..do this stuff once
      if (ifirst .eq. 0) then
       ifirst = 1
       open(unit=19,file='helm_table.dat',status='old')
c       open(unit=19,file='helm_table_big.dat',status='old')
c       open(unit=19,file='helm_table_test.dat',status='old')
c       open(unit=19,file='helm_table_ionize.dat',status='old')


c..read the normal helmholtz free energy table
       tlo   = 4.0d0
       thi   = 11.0d0
       tstp  = (thi - tlo)/float(jmax-1)
       tstpi = 1.0d0/tstp
       dlo   = -10.0d0
       dhi   = 11.0d0
       dstp  = (dhi - dlo)/float(imax-1)
       dstpi = 1.0d0/dstp

c..for the bigger  table
c       tlo   = 2.0d0
c       thi   = 11.0d0
c       tstp  = (thi - tlo)/float(jmax-1)
c       tstpi = 1.0d0/tstp
c       dlo   = -10.0d0
c       dhi   = 15.0d0
c       dstp  = (dhi - dlo)/float(imax-1)
c       dstpi = 1.0d0/dstp


       do j=1,jmax
        tsav = tlo + (j-1)*tstp
        t(j) = 10.0d0**(tsav)
        do i=1,imax
         dsav = dlo + (i-1)*dstp
         d(i) = 10.0d0**(dsav)
         read(19,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j),
     1            fddt(i,j),fdtt(i,j),fddtt(i,j)
        enddo
       enddo


c..read the pressure derivative with density table
       do j=1,jmax
        do i=1,imax
         read(19,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
        enddo
       enddo

c..read the electron chemical potential table
       do j=1,jmax
        do i=1,imax
         read(19,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
        enddo
       enddo

c..read the number density table
       do j=1,jmax
        do i=1,imax
         read(19,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
        enddo
       enddo

c..construct the temperature and density deltas and their inverses 
       do j=1,jmax-1
        dth          = t(j+1) - t(j)
        dt2         = dth * dth
        dti         = 1.0d0/dth
        dt2i        = 1.0d0/dt2
        dt_sav(j)   = dth
        dt2_sav(j)  = dt2
        dti_sav(j)  = dti
        dt2i_sav(j) = dt2i
       end do
       do i=1,imax-1
        dd          = d(i+1) - d(i)
        dd2         = dd * dd
        ddi         = 1.0d0/dd
        dd2i        = 1.0d0/dd2
        dd_sav(i)   = dd
        dd2_sav(i)  = dd2
        ddi_sav(i)  = ddi
        dd2i_sav(i) = dd2i
       enddo

       close(unit=19)
       write(6,*)
       write(6,*) 'finished reading eos table'
       write(6,04) 'imax=',imax,' jmax=',jmax
       write(6,03) 'temp(1)   =',t(1),' temp(jmax)   =',t(jmax)
       write(6,03) 'ye*den(1) =',d(1),' ye*den(imax) =',d(imax)
       write(6,*)


c..set the multipliers
      radmult  = 1.0d0
      ionmult  = 1.0d0
      elemult  = 1.0d0
      coulmult = 1.0d0


      end if



c..start of pipeline loop, normal executaion starts here
      eosfail = .false.
      do j=jlo_eos,jhi_eos

       if (temp_row(j) .le. 0.0) stop 'temp less than 0 in helmeos'
       if (den_row(j)  .le. 0.0) stop 'den less than 0 in helmeos'

       temp  = temp_row(j)
       den   = den_row(j)
       abar  = abar_row(j)
       zbar  = zbar_row(j)
       ytot1 = 1.0d0/abar
       ye    = ytot1 * zbar



c..initialize
       deni    = 1.0d0/den
       tempi   = 1.0d0/temp 
       kt      = kerg * temp
       ktinv   = 1.0d0/kt


c..radiation section:
       if (radmult .ne. 0) then
        prad    = asoli3 * temp * temp * temp * temp
        dpraddd = 0.0d0
        dpraddt = 4.0d0 * prad*tempi
        dpradda = 0.0d0
        dpraddz = 0.0d0

        erad    = 3.0d0 * prad*deni
        deraddd = -erad*deni
        deraddt = 3.0d0 * dpraddt*deni
        deradda = 0.0d0
        deraddz = 0.0d0

        srad    = (prad*deni + erad)*tempi
        dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
        dsraddt = (dpraddt*deni + deraddt - srad)*tempi
        dsradda = 0.0d0
        dsraddz = 0.0d0
       end if


c..ion section:
        xni     = avo * ytot1 * den
        dxnidd  = avo * ytot1
        dxnida  = -xni * ytot1

       if (ionmult .ne. 0) then
        pion    = xni * kt
        dpiondd = dxnidd * kt
        dpiondt = xni * kerg
        dpionda = dxnida * kt 
        dpiondz = 0.0d0

        eion    = 1.5d0 * pion*deni
        deiondd = (1.5d0 * dpiondd - eion)*deni
        deiondt = 1.5d0 * dpiondt*deni
        deionda = 1.5d0 * dpionda*deni
        deiondz = 0.0d0
    

c..sackur-tetrode equation for the ion entropy of 
c..a single ideal gas characterized by abar
        x       = abar*abar*sqrt(abar) * deni/avo
        s       = sioncon * temp
        z       = x * s * sqrt(s)
        y       = log(z)
        sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
        dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi
     1             - kergavo * deni * ytot1
        dsiondt = (dpiondt*deni + deiondt)*tempi - 
     1            (pion*deni + eion) * tempi*tempi 
     2            + 1.5d0 * kergavo * tempi*ytot1
        x       = avo*kerg/abar
        dsionda = (dpionda*deni + deionda)*tempi 
     1            + kergavo*ytot1*ytot1* (2.5d0 - y)
        dsiondz = 0.0d0
       end if



c..electron-positron section:
       if (elemult .ne. 0) then


c..assume complete ionization 
        xnem    = xni * zbar


c..enter the table with ye*den
        din = ye*den 


c..bomb proof the input
        if (temp .gt. t(jmax)) then
         write(6,01) 'temp=',temp,' t(jmax)=',t(jmax)
         write(6,*) 'temp too hot, off grid'       
         write(6,*) 'setting eosfail to true and returning'
         call flush(6)
         eosfail = .true.
         return
        end if
        if (temp .lt. t(1)) then
         write(6,01) 'temp=',temp,' t(1)=',t(1)
         write(6,*) 'temp too cold, off grid'
         write(6,*) 'setting eosfail to true and returning'
         call flush(6)
         eosfail = .true.
         return
        end if
        if (din  .gt. d(imax)) then
         write(6,01) 'den*ye=',din,' d(imax)=',d(imax)
         write(6,*) 'ye*den too big, off grid'
         write(6,*) 'setting eosfail to true and returning'
         call flush(6)
         eosfail = .true.
         return
        end if
        if (din  .lt. d(1)) then
         write(6,01) 'ye*den=',din,' d(1)=',d(1)
         write(6,*) 'ye*den too small, off grid'
         write(6,*) 'setting eosfail to true and returning'
         call flush(6)
         eosfail = .true.
         return
        end if

c..hash locate this temperature and density
        jat = int((log10(temp) - tlo)*tstpi) + 1
        jat = max(1,min(jat,jmax-1))
        iat = int((log10(din) - dlo)*dstpi) + 1
        iat = max(1,min(iat,imax-1))


c..access the table locations only once
        fi(1)  = f(iat,jat)
        fi(2)  = f(iat+1,jat)
        fi(3)  = f(iat,jat+1)
        fi(4)  = f(iat+1,jat+1)
        fi(5)  = ft(iat,jat)
        fi(6)  = ft(iat+1,jat)
        fi(7)  = ft(iat,jat+1)
        fi(8)  = ft(iat+1,jat+1)
        fi(9)  = ftt(iat,jat)
        fi(10) = ftt(iat+1,jat)
        fi(11) = ftt(iat,jat+1)
        fi(12) = ftt(iat+1,jat+1)
        fi(13) = fd(iat,jat)
        fi(14) = fd(iat+1,jat)
        fi(15) = fd(iat,jat+1)
        fi(16) = fd(iat+1,jat+1)
        fi(17) = fdd(iat,jat)
        fi(18) = fdd(iat+1,jat)
        fi(19) = fdd(iat,jat+1)
        fi(20) = fdd(iat+1,jat+1)
        fi(21) = fdt(iat,jat)
        fi(22) = fdt(iat+1,jat)
        fi(23) = fdt(iat,jat+1)
        fi(24) = fdt(iat+1,jat+1)
        fi(25) = fddt(iat,jat)
        fi(26) = fddt(iat+1,jat)
        fi(27) = fddt(iat,jat+1)
        fi(28) = fddt(iat+1,jat+1)
        fi(29) = fdtt(iat,jat)
        fi(30) = fdtt(iat+1,jat)
        fi(31) = fdtt(iat,jat+1)
        fi(32) = fdtt(iat+1,jat+1)
        fi(33) = fddtt(iat,jat)
        fi(34) = fddtt(iat+1,jat)
        fi(35) = fddtt(iat,jat+1)
        fi(36) = fddtt(iat+1,jat+1)
 

c..various differences
        xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
        xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
        mxt = 1.0d0 - xt
        mxd = 1.0d0 - xd

c..the six density and six temperature basis functions
        si0t =   psi0(xt)
        si1t =   psi1(xt)*dt_sav(jat)
        si2t =   psi2(xt)*dt2_sav(jat)

        si0mt =  psi0(mxt)
        si1mt = -psi1(mxt)*dt_sav(jat)
        si2mt =  psi2(mxt)*dt2_sav(jat)

        si0d =   psi0(xd)
        si1d =   psi1(xd)*dd_sav(iat)
        si2d =   psi2(xd)*dd2_sav(iat)

        si0md =  psi0(mxd)
        si1md = -psi1(mxd)*dd_sav(iat)
        si2md =  psi2(mxd)*dd2_sav(iat)

c..derivatives of the weight functions
        dsi0t =   dpsi0(xt)*dti_sav(jat)
        dsi1t =   dpsi1(xt)
        dsi2t =   dpsi2(xt)*dt_sav(jat)

        dsi0mt = -dpsi0(mxt)*dti_sav(jat)
        dsi1mt =  dpsi1(mxt)
        dsi2mt = -dpsi2(mxt)*dt_sav(jat)

        dsi0d =   dpsi0(xd)*ddi_sav(iat)
        dsi1d =   dpsi1(xd)
        dsi2d =   dpsi2(xd)*dd_sav(iat)

        dsi0md = -dpsi0(mxd)*ddi_sav(iat)
        dsi1md =  dpsi1(mxd)
        dsi2md = -dpsi2(mxd)*dd_sav(iat)

c..second derivatives of the weight functions
        ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
        ddsi1t =   ddpsi1(xt)*dti_sav(jat)
        ddsi2t =   ddpsi2(xt)
 
        ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
        ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
        ddsi2mt =  ddpsi2(mxt)

c        ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
c        ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
c        ddsi2d =   ddpsi2(xd)

c        ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
c        ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
c        ddsi2md =  ddpsi2(mxd)


c..the free energy
        free  = h5(iat,jat,
     1          si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
     2          si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to density
        df_d  = h5(iat,jat,
     1          si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
     2          dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

c..derivative with respect to temperature
        df_t = h5(iat,jat,
     1          dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,
     2          si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to density**2
c        df_dd = h5(iat,jat,
c     1          si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
c     2          ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

c..derivative with respect to temperature**2
        df_tt = h5(iat,jat,
     1        ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt,
     2          si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to temperature and density
        df_dt = h5(iat,jat,
     1          dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,
     2          dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)



c..now get the pressure derivative with density, chemical potential, and 
c..electron positron number densities
c..get the interpolation weight functions
        si0t   =  xpsi0(xt)
        si1t   =  xpsi1(xt)*dt_sav(jat)

        si0mt  =  xpsi0(mxt)
        si1mt  =  -xpsi1(mxt)*dt_sav(jat)

        si0d   =  xpsi0(xd)
        si1d   =  xpsi1(xd)*dd_sav(iat)

        si0md  =  xpsi0(mxd)
        si1md  =  -xpsi1(mxd)*dd_sav(iat)


c..derivatives of weight functions
        dsi0t  = xdpsi0(xt)*dti_sav(jat)
        dsi1t  = xdpsi1(xt)

        dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
        dsi1mt = xdpsi1(mxt)

        dsi0d  = xdpsi0(xd)*ddi_sav(iat)
        dsi1d  = xdpsi1(xd)

        dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
        dsi1md = xdpsi1(mxd)


c..look in the pressure derivative only once
        fi(1)  = dpdf(iat,jat)
        fi(2)  = dpdf(iat+1,jat)
        fi(3)  = dpdf(iat,jat+1)
        fi(4)  = dpdf(iat+1,jat+1)
        fi(5)  = dpdft(iat,jat)
        fi(6)  = dpdft(iat+1,jat)
        fi(7)  = dpdft(iat,jat+1)
        fi(8)  = dpdft(iat+1,jat+1)
        fi(9)  = dpdfd(iat,jat)
        fi(10) = dpdfd(iat+1,jat)
        fi(11) = dpdfd(iat,jat+1)
        fi(12) = dpdfd(iat+1,jat+1)
        fi(13) = dpdfdt(iat,jat)
        fi(14) = dpdfdt(iat+1,jat)
        fi(15) = dpdfdt(iat,jat+1)
        fi(16) = dpdfdt(iat+1,jat+1)

c..pressure derivative with density
        dpepdd  = h3(iat,jat,
     1                 si0t,   si1t,   si0mt,   si1mt,
     2                 si0d,   si1d,   si0md,   si1md)
        dpepdd  = max(ye * dpepdd,1.0d-30)



c..look in the electron chemical potential table only once
        fi(1)  = ef(iat,jat)
        fi(2)  = ef(iat+1,jat)
        fi(3)  = ef(iat,jat+1)
        fi(4)  = ef(iat+1,jat+1)
        fi(5)  = eft(iat,jat)
        fi(6)  = eft(iat+1,jat)
        fi(7)  = eft(iat,jat+1)
        fi(8)  = eft(iat+1,jat+1)
        fi(9)  = efd(iat,jat)
        fi(10) = efd(iat+1,jat)
        fi(11) = efd(iat,jat+1)
        fi(12) = efd(iat+1,jat+1)
        fi(13) = efdt(iat,jat)
        fi(14) = efdt(iat+1,jat)
        fi(15) = efdt(iat,jat+1)
        fi(16) = efdt(iat+1,jat+1)


c..electron chemical potential etaele
        etaele  = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2               si0d,   si1d,   si0md,   si1md)


c..derivative with respect to density
        x       = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2              dsi0d,  dsi1d,  dsi0md,  dsi1md)
        detadd  = ye * x  

c..derivative with respect to temperature
        detadt  = h3(iat,jat,
     1              dsi0t,  dsi1t,  dsi0mt,  dsi1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to abar and zbar
       detada = -x * din * ytot1
       detadz =  x * den * ytot1



c..look in the number density table only once
        fi(1)  = xf(iat,jat)
        fi(2)  = xf(iat+1,jat)
        fi(3)  = xf(iat,jat+1)
        fi(4)  = xf(iat+1,jat+1)
        fi(5)  = xft(iat,jat)
        fi(6)  = xft(iat+1,jat)
        fi(7)  = xft(iat,jat+1)
        fi(8)  = xft(iat+1,jat+1)
        fi(9)  = xfd(iat,jat)
        fi(10) = xfd(iat+1,jat)
        fi(11) = xfd(iat,jat+1)
        fi(12) = xfd(iat+1,jat+1)
        fi(13) = xfdt(iat,jat)
        fi(14) = xfdt(iat+1,jat)
        fi(15) = xfdt(iat,jat+1)
        fi(16) = xfdt(iat+1,jat+1)

c..electron + positron number densities
       xnefer   = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to density
       x        = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2              dsi0d,  dsi1d,  dsi0md,  dsi1md)
       x = max(x,1.0d-30)
       dxnedd   = ye * x

c..derivative with respect to temperature
       dxnedt   = h3(iat,jat,
     1              dsi0t,  dsi1t,  dsi0mt,  dsi1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to abar and zbar
       dxneda = -x * din * ytot1
       dxnedz =  x  * den * ytot1



c..the desired electron-positron thermodynamic quantities

c..dpepdd at high temperatures and low densities is below the
c..floating point limit of the subtraction of two large terms.
c..since dpresdd doesn't enter the maxwell relations at all, use the
c..bicubic interpolation done above instead of the formally correct expression
        x       = din * din
        pele    = x * df_d
        dpepdt  = x * df_dt
c        dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)
        s       = dpepdd/ye - 2.0d0 * din * df_d
        dpepda  = -ytot1 * (2.0d0 * pele + s * din)
        dpepdz  = den*ytot1*(2.0d0 * din * df_d  +  s)


        x       = ye * ye
        sele    = -df_t * ye
        dsepdt  = -df_tt * ye
        dsepdd  = -df_dt * x
        dsepda  = ytot1 * (ye * df_dt * din - sele)
        dsepdz  = -ytot1 * (ye * df_dt * den  + df_t)


        eele    = ye*free + temp * sele
        deepdt  = temp * dsepdt
        deepdd  = x * df_d + temp * dsepdd
        deepda  = -ye * ytot1 * (free +  df_d * din) + temp * dsepda
        deepdz  = ytot1* (free + ye * df_d * den) + temp * dsepdz

c..end of elemult
       end if 




c..coulomb section:
c..initialize
       if (coulmult .ne. 0) then
        pcoul    = 0.0d0
        dpcouldd = 0.0d0
        dpcouldt = 0.0d0
        dpcoulda = 0.0d0
        dpcouldz = 0.0d0
        ecoul    = 0.0d0
        decouldd = 0.0d0
        decouldt = 0.0d0
        decoulda = 0.0d0
        decouldz = 0.0d0
        scoul    = 0.0d0
        dscouldd = 0.0d0
        dscouldt = 0.0d0
        dscoulda = 0.0d0
        dscouldz = 0.0d0


c..uniform background corrections only 
c..from yakovlev & shalybkov 1989 
c..lami is the average ion seperation
c..plasg is the plasma coupling parameter
        z        = forth * pi
        s        = z * xni
        dsdd     = z * dxnidd
        dsda     = z * dxnida

        lami     = 1.0d0/s**third
        inv_lami = 1.0d0/lami
        z        = -third * lami
        lamidd   = z * dsdd/s
        lamida   = z * dsda/s

        plasg    = zbar*zbar*esqu*ktinv*inv_lami
        z        = -plasg * inv_lami 
        plasgdd  = z * lamidd
        plasgda  = z * lamida
        plasgdt  = -plasg*ktinv * kerg
        plasgdz  = 2.0d0 * plasg/zbar


c..yakovlev & shalybkov 1989 equations 82, 85, 86, 87
         if (plasg .ge. 1.0) then
          x        = plasg**(0.25d0) 
          y        = avo * ytot1 * kerg 
          ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
          pcoul    = third * den * ecoul
          scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x
     1              + d1 * (log(plasg) - 1.0d0) - e1)

          y        = avo*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
          decouldd = y * plasgdd 
          decouldt = y * plasgdt + ecoul/temp
          decoulda = y * plasgda - ecoul/abar
          decouldz = y * plasgdz

          y        = third * den
          dpcouldd = third * ecoul + y*decouldd
          dpcouldt = y * decouldt
          dpcoulda = y * decoulda
          dpcouldz = y * decouldz


          y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x+1.25d0*c1/x+d1)
          dscouldd = y * plasgdd
          dscouldt = y * plasgdt
          dscoulda = y * plasgda - scoul/abar
          dscouldz = y * plasgdz


c..yakovlev & shalybkov 1989 equations 102, 103, 104
         else if (plasg .lt. 1.0) then
          x        = plasg*sqrt(plasg)
          y        = plasg**b2
          z        = c2 * x - third * a2 * y
          pcoul    = -pion * z
          ecoul    = 3.0d0 * pcoul/den
          scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)

          s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
          dpcouldd = -dpiondd*z - pion*s*plasgdd
          dpcouldt = -dpiondt*z - pion*s*plasgdt
          dpcoulda = -dpionda*z - pion*s*plasgda
          dpcouldz = -dpiondz*z - pion*s*plasgdz

          s        = 3.0d0/den
          decouldd = s * dpcouldd - ecoul/den
          decouldt = s * dpcouldt
          decoulda = s * dpcoulda
          decouldz = s * dpcouldz

          s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x-a2*(b2-1.0d0)*y)
          dscouldd = s * plasgdd
          dscouldt = s * plasgdt
          dscoulda = s * plasgda - scoul/abar
          dscouldz = s * plasgdz
         end if



c..bomb proof
        x   = prad + pion + pele + pcoul
        if (x .le. 0.0) then

c         write(6,*) 
c         write(6,*) 'coulomb corrections are causing a negative pressure'
c         write(6,*) 'setting all coulomb corrections to zero'
c         write(6,*) 

         pcoul    = 0.0d0
         dpcouldd = 0.0d0
         dpcouldt = 0.0d0
         dpcoulda = 0.0d0
         dpcouldz = 0.0d0
         ecoul    = 0.0d0
         decouldd = 0.0d0
         decouldt = 0.0d0
         decoulda = 0.0d0
         decouldz = 0.0d0
         scoul    = 0.0d0
         dscouldd = 0.0d0
         dscouldt = 0.0d0
         dscoulda = 0.0d0
         dscouldz = 0.0d0
        end if
       end if




c..sum all the components
       pres    = prad + pion + pele + pcoul
       ener    = erad + eion + eele + ecoul
       entr    = srad + sion + sele + scoul

       dpresdd = dpraddd + dpiondd + dpepdd + dpcouldd 
       dpresdt = dpraddt + dpiondt + dpepdt + dpcouldt
       dpresda = dpradda + dpionda + dpepda + dpcoulda
       dpresdz = dpraddz + dpiondz + dpepdz + dpcouldz

       denerdd = deraddd + deiondd + deepdd + decouldd
       denerdt = deraddt + deiondt + deepdt + decouldt
       denerda = deradda + deionda + deepda + decoulda
       denerdz = deraddz + deiondz + deepdz + decouldz

       dentrdd = dsraddd + dsiondd + dsepdd + dscouldd
       dentrdt = dsraddt + dsiondt + dsepdt + dscouldt
       dentrda = dsradda + dsionda + dsepda + dscoulda
       dentrdz = dsraddz + dsiondz + dsepdz + dscouldz


c..the temperature and density exponents (c&g 9.81 9.82) 
c..the specific heat at constant volume (c&g 9.92)
c..the third adiabatic exponent (c&g 9.93)
c..the first adiabatic exponent (c&g 9.97) 
c..the second adiabatic exponent (c&g 9.105)
c..the specific heat at constant pressure (c&g 9.98) 
c..and relativistic formula for the sound speed (c&g 14.29)
       zz    = pres*deni
       zzi   = den/pres
       chit  = temp/pres * dpresdt
       chid  = dpresdd*zzi
       cv    = denerdt
       x     = zz * chit/(temp * cv)
       gam3  = x + 1.0d0
       gam1  = chit*x + chid
       nabad = x/gam1
       gam2  = 1.0d0/(1.0d0 - nabad)
       cp    = cv * gam1/chid
       z     = 1.0d0 + (ener + light2)*zzi
       sound = clight * sqrt(gam1/z)


c..maxwell relations; each is zero if the consistency is perfect
       x   = den * den

       dse = temp*dentrdt/denerdt - 1.0d0

       dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0

       dsp = -dentrdd*x/dpresdt - 1.0d0


c..store this row
        ptot_row(j)   = pres
        dpt_row(j)    = dpresdt
        dpd_row(j)    = dpresdd
        dpa_row(j)    = dpresda   
        dpz_row(j)    = dpresdz

        etot_row(j)   = ener
        det_row(j)    = denerdt
        ded_row(j)    = denerdd
        dea_row(j)    = denerda   
        dez_row(j)    = denerdz

        stot_row(j)   = entr 
        dst_row(j)    = dentrdt
        dsd_row(j)    = dentrdd
        dsa_row(j)    = dentrda        
        dsz_row(j)    = dentrdz

        prad_row(j)   = prad
        erad_row(j)   = erad
        srad_row(j)   = srad 

        pion_row(j)   = pion
        eion_row(j)   = eion
        sion_row(j)   = sion 
        xni_row(j)    = xni

        pele_row(j)   = pele
        ppos_row(j)   = 0.0d0
        dpept_row(j)  = dpepdt
        dpepd_row(j)  = dpepdd
        dpepa_row(j)  = dpepda  
        dpepz_row(j)  = dpepdz

        eele_row(j)   = eele
        epos_row(j)   = 0.0d0
        deept_row(j)  = deepdt
        deepd_row(j)  = deepdd
        deepa_row(j)  = deepda   
        deepz_row(j)  = deepdz

        sele_row(j)   = sele 
        spos_row(j)   = 0.0d0
        dsept_row(j)  = dsepdt 
        dsepd_row(j)  = dsepdd 
        dsepa_row(j)  = dsepda        
        dsepz_row(j)  = dsepdz

        xnem_row(j)   = xnem
        xne_row(j)    = xnefer
        dxnet_row(j)  = dxnedt
        dxned_row(j)  = dxnedd
        dxnea_row(j)  = dxneda
        dxnez_row(j)  = dxnedz
        xnp_row(j)    = 0.0d0
        zeff_row(j)   = zbar

        etaele_row(j) = etaele
        detat_row(j)  = detadt
        detad_row(j)  = detadd
        detaa_row(j)  = detada
        detaz_row(j)  = detadz
        etapos_row(j) = 0.0d0

        pcou_row(j)   = pcoul
        ecou_row(j)   = ecoul
        scou_row(j)   = scoul 
        plasg_row(j)  = plasg

        dse_row(j)    = dse
        dpe_row(j)    = dpe
        dsp_row(j)    = dsp

        cv_row(j)     = cv
        cp_row(j)     = cp
        gam1_row(j)   = gam1
        gam2_row(j)   = gam2
        gam3_row(j)   = gam3
        cs_row(j)     = sound

c..end of pipeline loop
      enddo
      return
      end



