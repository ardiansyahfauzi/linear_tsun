  !! ----------------------------------------------------------------------- !!
  !>
  !! Bathymetry depth hh(:,:) ( in m-unit ) must be filled here. 
  !<
  block

    integer :: i, j
    !! --
   
    ! modified by Fauzi
    open(21,file='bath_in_new.dat', status='unknown')
    
    do i=1,nx
      read(21,*)(hh(i,j), j=1,ny)
    enddo
    close(21)
    !-------------------

  end block
