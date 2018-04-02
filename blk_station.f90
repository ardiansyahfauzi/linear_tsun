  !! ----------------------------------------------------------------------- !!
  !>
  !! Set station location.
  !! nst, xst(:), yst(:), ist(:), jst(:) should be set here with memory alloc.
  !<
  block

    integer :: i, j
    
    open(21,file = 'gages.dat', status = 'unknown') 
    
    read(21,*)nst
    
    allocate( ist(nst), jst(nst), xst(nst), yst(nst) )
    allocate( wav(nst,0:nt) )
   
    do i=1,nst
       read(21,'(i3, i5, i5)')j, ist(i), jst(i)
    enddo
    
    wav(:,:) = 0.0
    
    close(21)
       
  end block
