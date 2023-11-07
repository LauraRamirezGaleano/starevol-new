      PROGRAM runevol

#ifdef GYRE
      use evolstell_gyre
#endif
      implicit none 

      integer i, i_mod, ibeg, iend, filesize
      character (len=5) :: istr, iendstr, ibegstr, bintype, phase
      character(len=200) :: filename, phdir, subinis
      character(len=24) :: starevolog, stopage
      character(len=25) :: filex, hostname, root_name, v_file
      character(len=100) :: stardate
      character(len=10) :: date
      integer error, size_error
      character(len=30) :: seqa, seqaz, seqb, seqbz, seqd, seqdz
      integer sizetot, sizea, sizeaz, sizeb, sizebz, sized, sizedz

      ! Initialize GYRE
#ifdef GYRE
      call evolstell_startup('gyre.in','evolgyre')
#endif
      ! If runevol.e launched by itself
      if ( COMMAND_ARGUMENT_COUNT() == 0 ) then

	  call starevol_init

      ! If it is launched by batchesun_3.40
      else if ( COMMAND_ARGUMENT_COUNT() == 10 ) then
      
         ! Collect arguments
         call getarg(1,filename)
         starevolog=trim(filename)//'.rec'
         call getarg(2,phdir)   ! character string
         call getarg(3,subinis) ! character string         
         call getarg(4,ibegstr) ! convert to integer
         read(ibegstr,'(I5)') ibeg 
         call getarg(5,iendstr) ! convert to integer
         read(iendstr,'(I5)') iend
         call getarg(6,phase)   ! character string
         call getarg(7, stopage) ! character string         
         call getarg(8, bintype) ! get it in char
         call getarg(9, root_name)
         call getarg(10, v_file)
      
         i=ibeg
         do while ( i <= iend )
         
            i_mod=i + 10000

            write(istr,'(I5)') i_mod 
            istr=istr(2:5)
            
            ! Run bash_before         
            call system('$DIR_BATCH/bash_before '//trim(filename)
     &        //' '//trim(ibegstr)//' '//trim(iendstr)//' '//trim(istr))

            ! Retrieve the error from bash_before if there was one
            inquire(file='error.log', size=size_error)
            if ( size_error /=0 ) then
               open(unit=1, file='error.log', action='read') 
               read(1,*) error
               close(1)
               if ( error /= 0) then 
                 call exit(error)
               endif
            endif
            
            ! If the computation is run for the first time, tables must be read
            if ( i == ibeg ) then
               call starevol_init
               
            ! If it is not, tables are already in memory 
            else
               call starevol
      
            end if

            ! Read execution date to use it in bash_after_sun
            open(unit=1, file='date.log', action='read')
            read(1, *) stardate
            close(1)
            
            ! Run bash_after_sun      
            call system('$DIR_BATCH/bash_after_sun '//trim(filename)//
     &        ' '//trim(phdir)//' '//trim(subinis)//' '//trim(ibegstr)//
     &        ' '//istr//' '//trim(iendstr)//' '//trim(stardate)//' '
     &        //trim(phase)//' '//trim(stopage))

            ! Retrieve the error from bash_after_sun if there was one
            inquire(file='error.log', size=size_error)
            if ( size_error /=0 ) then
               open(unit=1, file='error.log', action='read') 
               read(1,*) error
               close(1)
               if ( error /= 0) then
                  call exit(error)
               endif
            endif
            
            ! Run bash_evolfile           
            call system('$DIR_BATCH/bash_evolfile '//
     &        trim(filename)//' '//trim(phdir)//' ' //trim(subinis)//
     &        ' '//trim(root_name)//' '//trim(v_file))
            
            ! Retrieve the error from bash_after_sun if there was one
            inquire(file='error.log', size=size_error)
            if ( size_error /=0 ) then
               open(unit=1, file='error.log', action='read') 
               read(1,*) error
               close(1)
               if ( error /= 0) then
                  call exit(error)
               endif
            endif
            
            ! Go to next iteration
            i=i+1
         end do
         
      ! If it is launched by batche_3.40
      else 
      
         ! Collect arguments
         call getarg(1,filename)
         starevolog=trim(filename)//'.rec'
         call getarg(2,phdir)   ! character string
         call getarg(3,subinis) ! character string
         call getarg(4,ibegstr) ! convert to integer
         read(ibegstr,'(I5)') ibeg
         call getarg(5,iendstr) ! convert to integer
         read(iendstr,'(I5)') iend
         call getarg(6,phase)   ! character string
         call getarg(7, bintype) ! get it in char
         call getarg(8, root_name)
         call getarg(9, v_file)
      
         i=ibeg
         do while ( i <= iend )
         
            i_mod=i + 10000

            write(istr,'(I5)') i_mod
            istr=istr(2:5)
         
            ! Run bash_before
            call system('$DIR_BATCH/bash_before '//trim(filename)
     &        //' '//trim(ibegstr)//' '//trim(iendstr)//' '//trim(istr))

            ! Retrieve the error from bash_before if there was one
            inquire(file='error.log', size=size_error)
            if ( size_error /=0 ) then
               open(unit=1, file='error.log', action='read') 
               read(1,*) error
               close(1)
               if ( error /= 0) then
                  call exit(error)
               endif
            endif
            
            ! If the computation is run for the first time, tables must be read
            if ( i == ibeg ) then
               call starevol_init
            
            ! If it is not, tables are already in memory 
            else
               call starevol
         
            end if
            
            ! Read execution date to use it in bash_after_sun
            open(unit=1, file='date.log', action='read')
            read(1, *) stardate
            close(1)
            
            ! Run bash_after      
            call system('$DIR_BATCH/bash_after '//trim(filename)//
     &        ' '//trim(phdir)//' '//trim(subinis)//' '//trim(ibegstr)//
     &        ' '//istr//' '//trim(iendstr)//' '//trim(stardate)//' '
     &        //trim(phase))
     
            ! Retrieve the error from bash_after if there was one
            inquire(file='error.log', size=size_error)
            if ( size_error /=0 ) then
               open(unit=1, file='error.log', action='read') 
               read(1,*) error
               close(1)
               if ( error /= 0) then
                  call exit(error)
               endif
            endif
       
            ! Run bash_evolfile        
            call system('$DIR_BATCH/bash_evolfile '//
     &        trim(filename)//' '//trim(phdir)//' ' //trim(subinis)//
     &        ' '//trim(root_name)//' '//trim(v_file))
     
            ! Retrieve the error from bash_evolfile if there was one
            inquire(file='error.log', size=size_error)
            if ( size_error /=0 ) then
               open(unit=1, file='error.log', action='read') 
               read(1,*) error
               close(1)
               if ( error /= 0) then
                  call exit(error)
               endif
            endif
            
            ! Go to next iteration
            i=i+1
         end do
      
      end if
      
      ! Finalize GYRE
#ifdef GYRE      
      call evolstell_final
#endif
      end program runevol
      
      
