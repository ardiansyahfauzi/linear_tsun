  !! ----------------------------------------------------------------------- !!
  !>
  !! Tsunami initial height.
  !! Tsunami height eta(:,:) ( in m-unit ) must be defined here.
  !<
  !! --
  block
   
    integer :: i, j
    
    open(22,file="example/eta.dat", status='old')
    
    do i=1,nx
      read(22,*)(eta(i,j), j=1,ny)
    enddo
    
    close(22)

  end block
