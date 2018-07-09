  !! ----------------------------------------------------------------------- !!
  !>
  !! Set station location.
  !<
  block

    integer :: i, j
    
    open(21,file = "example/gages.dat", status = 'unknown') 
    
    read(21,*)nst
    
    allocate( ist(nst), jst(nst), xst(nst), yst(nst) )
    allocate( wav(0:int(nt / dt), nst+1) )
   
    do i=1,nst
       read(21,'(i3, i5, i5)')j, ist(i), jst(i)
    enddo
    
    wav(:,:) = 0.0
    
    close(21)
       
  end block
