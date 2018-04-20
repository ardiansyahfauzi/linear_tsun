!! ------------------------------------------------------------------------- !!
!>
!!
!! Linear Long Wave (LLW) for tsunami propagation in spherical coordinate system.
!! This code is modification of Maeda (2015) code which is written in
!! cartesian coordinate system.
!! 
!! The governing equation of LLW is adopted from Nakamura's model.
!! This code is partially modified by A. Fauzi
!! Last modified on 2018.2.1
!!
!<
!! ------------------------------------------------------------------------- !!
program llw

  implicit none

  !! fixed parameters
  integer, parameter :: STDERR  = 0            !< Standard Error
  integer, parameter :: STDOUT  = 6            !< Standard Output
  real,    parameter :: pi      = atan(1.0)*4
  real,    parameter :: g0      = 9.80665      !< gravity
  integer, parameter :: na      = 20           !< absorber thickness
  real,    parameter :: r       = 6371012.0    !< earth radius (m)
  real,    parameter :: omega   = 2.0 * acos(-1.0) / 86400 !< earth rotational velocity (rad/s)
  character(len=20)  :: fld
  
  !! control parameters
  integer            :: nx,  ny                !< model size
  integer            :: nt, tint               !< time step size and gauge interval
  real               :: dx,  dy                !< grid width
  real               :: dt                     !< time step
  integer            :: mov                    !< snapshots

  !! derived global variables
  real               :: dxi, dyi               !< inverse grid width
  real               :: xbeg, ybeg             !< model edge
  integer            :: nst                    !< number of stations
  
  !! time
  integer            :: start, finish
  integer            :: clock_rate, clock_max

  !! arrays
  real,    allocatable :: eta(:,:)             !< tsunami height
  real,    allocatable :: mm(:,:)              !< tsunami vel
  real,    allocatable :: nn(:,:)              !< tsunami vel
  real,    allocatable :: hh(:,:)              !< bathymetry
  real,    allocatable :: hm(:,:)              !< x-averaged bathymetry
  real,    allocatable :: hn(:,:)              !< y-averaged bathymetry
  real,    allocatable :: fh(:,:)              !< land filter
  real,    allocatable :: fm(:,:)              !< land filter
  real,    allocatable :: fn(:,:)              !< land filter
  real,    allocatable :: g1x(:,:), g2x(:,:)   !< PML absorber
  real,    allocatable :: g1y(:,:), g2y(:,:)   !< PML absorber
  real,    allocatable :: eta_x(:,:)           !< PML eta(x)
  real,    allocatable :: wav(:,:)             !< waveform at stations
  integer, allocatable :: ist(:), jst(:)       !< station location (i,j)
  real,    allocatable :: xst(:), yst(:)       !< station location (x,y)
  real,    allocatable :: theta(:), thetal(:)  !< latitude coordinate

  !! ----------------------------------------------------------------------- !!

  !! parameter input
  include "blk_param.f90"

  !! ----------------------------------------------------------------------- !!
  !>
  !! memory allocation
  !<
  !! --
  block

    integer :: j
    real    :: degrad, secrad

    allocate( eta(1:nx,1:ny) )
    allocate( mm (1:nx,1:ny) )
    allocate( nn (1:nx,1:ny) )
    allocate( hh (1:nx,1:ny) )
    allocate( hm (1:nx,1:ny) )
    allocate( hn (1:nx,1:ny) )
    allocate( fh (1:nx,1:ny) )
    allocate( fm (1:nx,1:ny) )
    allocate( fn (1:nx,1:ny) )
    allocate( g1x(2,nx), g2x(2,nx), g1y(2,ny), g2y(2,ny) )
    allocate( eta_x(nx,ny) )
    allocate( theta(ny), thetal(ny) )

    call system_clock(start,clock_rate,clock_max)
    eta(:,:) = 0.0
    mm (:,:) = 0.0
    nn (:,:) = 0.0
    
    degrad = pi / 180.0
    secrad = pi / (180.0 * 3600)
    
    ybeg = ybeg * degrad
    
    dy   = dy * secrad
    dx   = dx * secrad
!$OMP PARALLEL PRIVATE(j)
!$OMP DO
    do j=1,ny
       theta (j) = ybeg + (j-1) * dy
       thetal(j) = ybeg + (j-1) * dy - (0.5 * dy)
    enddo
!$OMP END DO
!$OMP END PARALLEL

    dxi = 1.0 / dx
    dyi = 1.0 / dy 

  end block

  include 'blk_initheight.f90'

  include 'blk_bathymetry.f90'

  include 'blk_station.f90'
  
  !! Create /out folder
  !! --
  if(mov.eq.0)then
     write(fld,'(a)')'out'
     call system("mkdir " // trim(adjustl(fld)))
  endif
  
  !! ----------------------------------------------------------------------- !!
  !>
  !! PML absorber settings
  !<
  !! --
  block

    real :: d1, d2, b1, b2
    integer :: i, j
    real :: d0, b0, R0
    real :: xxb, xxc
    real :: hx
    integer, parameter :: pd = 2
    integer, parameter :: pb = 2
    real,    parameter :: c0 = 500.0 !! assumed phase speed
    !! --


    !! damping profile constants
    hx = na * dx * r
    R0 = 10**( - ( log10( real(na) ) - 1 ) / log10( 2.0 )  - 3.0 )
    d0 = - ( 1.0 / (2.0*hx) ) * ( pd +1 ) * c0 * log( R0 )
    b0 = 2.0

    !! initialize eta: assume initial amplitude wihtin PML layer is zero
    eta_x(1:nx,1:ny) = 0.0

    !! initialize
!$OMP PARALLEL SHARED (g1x,g2x,g1y,g2y) PRIVATE(i,j)
!$OMP DO
    do i=1, 2
       do j=1, nx
          g1x(i,j) = 1.0
          g2x(i,j) = 1.0
       end do
       do j=1, ny
          g1y(i,j) = 1.0
          g2y(i,j) = 1.0
       end do
    end do
!$OMP END DO
!$OMP END PARALLEL    

    !! set absorber

    do i=1, na

      !! distance from PML boundary
       xxb = max(hx  - (i    ) * dx * r,0.)
       xxc = max(hx  - (i+0.5) * dx * r,0.)

       d1 = d0 * ( xxb / hx )**pd
       b1 = 1.0 + ( b0 - 1.0 ) * ( xxb / hx )**pb

       d2 = d0 * ( xxc / hx )**pd
       b2 = 1.0 + ( b0 - 1.0 ) * ( xxc / hx )**pb

       g1x(1,i) = (b1 - d1/2.0*dt) / (b1 + d1/2*dt)
       g1x(2,i) = (b1 - d2/2.0*dt) / (b1 + d2/2*dt)
       g2x(1,i) =  1.0             / (b1 + d1/2*dt)
       g2x(2,i) =  1.0             / (b1 + d2/2*dt)

       g1y(:,i) = g1x(:,i)
       g2y(:,i) = g2x(:,i)

       ! Opposite Side; Exchange Staggered Grid
       g1x(1,nx-i+1) = g1x(2,i)
       g1x(2,nx-i+1) = g1x(1,i)
       g2x(1,nx-i+1) = g2x(2,i)
       g2x(2,nx-i+1) = g2x(1,i)

       g1y(1,ny-i+1) = g1x(2,i)
       g1y(2,ny-i+1) = g1x(1,i)
       g2y(1,ny-i+1) = g2x(2,i)
       g2y(2,ny-i+1) = g2x(1,i)

    end do

  end block


  !! ----------------------------------------------------------------------- !!
  !>
  !! Use 10-m cut off for linear long wave
  !<
  block

    integer :: i, j
    real, parameter :: CUTOFF_DEPTH = 10.0
!$OMP PARALLEL SHARED(eta,hh) PRIVATE(i,j)
!$OMP DO
    do j=1, ny
       do i=1, nx
          if( hh(i,j) < 0.0 ) then
             eta(i,j) = 0.0
          else if ( hh(i,j) < CUTOFF_DEPTH ) then
             hh(i,j) = 10.0
          end if
       end do
    end do
!$OMP END DO
!$OMP END PARALLEL 
  end block


  !! ----------------------------------------------------------------------- !!
  !>
  !! bathymetry averaging for staggered grid calculation
  !<
  block

    integer :: i, j

    ! x-averaged bathymetry
!$OMP PARALLEL SHARED(hm,hn) PRIVATE(i,j)
!$OMP DO
    do j=1, ny
       hm(1,j) = max( hh(1,j), 0.0 )
       do i=2, nx
          hm(i,j) = ( hh(i,j)+hh(i-1,j) ) / 2.0
          if( hh(i,j) <= 0.0  .or.  hh(i-1,j) <= 0.0 ) hm(i,j) = 0.0
       end do
    end do

    ! y-averaged bathymetry
!$OMP DO
    do i=1, nx
       do j=2, ny
          hn(i,j) = ( hh(i,j) + hh(i,j-1) ) / 2.0
          if( hh(i,j) <= 0.0  .or.  hh(i,j-1) <= 0.0 ) hn(i,j) = 0.0
       end do
       hn(i,1) = max( hh(i,1), 0.0 )
    end do
!$OMP END DO
!$OMP END PARALLEL 
  end block


  !! ----------------------------------------------------------------------- !!
  !>
  !! define land-filter matrix for reflection boundaries
  !<
  block

    integer :: i, j

    fm(:,:) = 1.0
    fn(:,:) = 1.0
    fh(:,:) = 1.0
!$OMP PARALLEL SHARED(fm,fn,fh) PRIVATE(i,j)
!$OMP DO
    do j=1, ny
       do i=1, nx
          if( hm(i,j) < 0.0 ) fm(i,j) = 0.0
          if( hn(i,j) < 0.0 ) fn(i,j) = 0.0
          if( hh(i,j) < 0.0 ) fh(i,j) = 0.0
       end do
    end do
!$OMP END DO
!$OMP END PARALLEL 
  end block


  !! ----------------------------------------------------------------------- !!
  !>
  !! CFL stability condition check
  !<
  block

    if( min(dx * r,dy * r)/dt < sqrt( 2*g0*maxval(hh(:,:)) ) ) then
       write(STDERR,*) "ERROR: stability condition is violated: ",  &
            min(dx * r,dy * r)/dt,  sqrt( 2*g0*maxval(hh(:,:)) )
       stop
    end if

  end block



  !! ----------------------------------------------------------------------- !!
  !>
  !! time stepping
  !<
  block

    integer :: it
    real :: dxeta(nx,ny), dyeta(nx,ny)
    real :: dxm(nx,ny), dyn(nx,ny), ttt, elapsed, remaining, prog
    !! --

    do it= 1, nt

       if( mod(it*dt, real(60)) == 0 ) then
          ttt  = it * dt / 60.0
          prog = real(it) / real(nt)
          call system_clock(finish, clock_rate, clock_max)
          elapsed   = real(finish-start)/real(clock_rate)/60.0
          remaining = ((1 - prog) / prog) * elapsed !* clock_rate/60.0
          write(*,'(a, f7.2, a, a, f7.2, a, f7.2, a)')'Elapsed time = ',elapsed,&
                ' min; ', 'Remaining time = ', remaining, ' min; ', prog * 100, '/100%'
       end if

       !! ------------------------------------------------------------------ !!
       !>
       !! difference of eta with respect to x, y
       !<
       !! --
       block

         integer :: i, j
         !! --

!$OMP PARALLEL SHARED(dxeta,dyeta) PRIVATE(i,j)
!$OMP DO
         do j=1, ny
            do i=2, nx
               dxeta(i,j) = ( eta(i,j) - eta(i-1,j) ) * dxi
            end do

            dxeta(1,j) = ( eta(1,j) -        0.0 ) * dxi
         end do
!$OMP DO
         do i=1, nx
            do j=2, ny
               dyeta(i,j) = ( eta(i,j) - eta(i,j-1) ) * dyi
            end do

            dyeta(i,1) = ( eta(i,1) -        0.0 ) * dyi
         end do
!$OMP END DO
!$OMP END PARALLEL 
!write(*,*)dxeta(1,100), dxeta(780,382), dyeta(100,1), dyeta(780,382)
!pause
       end block
       
       !! ------------------------------------------------------------------ !!
       !>
       !! update velocity in the interior domain
       !<
       !!--
       update_vel: block

         integer :: i, j
         
         !! --

!$OMP PARALLEL SHARED(mm,nn) PRIVATE(i,j)
!$OMP DO
         do j=na+1, ny-na
            do i=na+1, nx-na

               mm(i,j) = mm(i,j) - (g0*hm(i,j)*dxeta(i,j)*dt*(1/(cos(theta(j)) * r)) &
                       + 0.0) * fm(i,j)
               nn(i,j) = nn(i,j) - (g0*hn(i,j)*dyeta(i,j)*dt*(1/r) &
                       - 0.0) * fn(i,j)
                       
            end do   
         end do
!$OMP END DO
!$OMP END PARALLEL 
       end block update_vel



       !! ------------------------------------------------------------------ !!
       !>
       !! update velocity in the PML layer
       !<
       !! --
       block

         integer :: i, j
         !! --

         !! South
!$OMP PARALLEL SHARED(mm,nn) PRIVATE(i,j)
!$OMP DO
         do j=1, na
            do i=1, nx
               mm(i,j) = g1x(1,i) * mm(i,j) &
                       - g2x(1,i) * g0 * hm(i,j) * dxeta(i,j) * dt * (1/(cos(theta(j)) * r)) * fm(i,j)
               nn(i,j) = g1y(1,j) * nn(i,j) &
                       - g2y(1,j) * g0 * hn(i,j) * dyeta(i,j) * dt * (1/r) * fn(i,j)
            end do
         end do

         !! North
!$OMP DO
         do j=ny-na+1,ny
            do i=1, nx
               mm(i,j) = g1x(1,i) * mm(i,j) &
                       - g2x(1,i) * g0 * hm(i,j) * dxeta(i,j) * dt * (1/(cos(theta(j)) * r)) * fm(i,j)
               nn(i,j) = g1y(1,j) * nn(i,j) &
                       - g2y(1,j) * g0 * hn(i,j) * dyeta(i,j) * dt * (1/r) * fn(i,j)
            end do
         end do

         !! West
!$OMP DO
         do j=na+1, ny-na
            do i=1, na
               mm(i,j) = g1x(1,i) * mm(i,j) &
                       - g2x(1,i) * g0 * hm(i,j) * dxeta(i,j) * dt * (1/(cos(theta(j)) * r)) * fm(i,j)
               nn(i,j) = g1y(1,j) * nn(i,j) &
                       - g2y(1,j) * g0 * hn(i,j) * dyeta(i,j) * dt * (1/r) * fn(i,j)
            end do
         end do

         !! East
!$OMP DO
         do j=na+1, ny-na
            do i=nx-na+1,nx
               mm(i,j) = g1x(1,i) * mm(i,j) &
                       - g2x(1,i) * g0 * hm(i,j) * dxeta(i,j) * dt * (1/(cos(theta(j)) * r)) * fm(i,j)
               nn(i,j) = g1y(1,j) * nn(i,j) &
                       - g2y(1,j) * g0 * hn(i,j) * dyeta(i,j) * dt * (1/r) * fn(i,j)
            end do
         end do         
!$OMP END DO
!$OMP END PARALLEL
       end block



       !! ------------------------------------------------------------------ !!
       !>
       !! difference of velocities with respect to x, y
       !<
       !! --
       block

         integer :: i, j
         !! --
!$OMP PARALLEL SHARED(dxm,dyn) PRIVATE(i,j)
!$OMP DO
         do j=1, ny
            dxm(nx,j) = ( 0.0 - mm(nx,j)  ) * dxi
            do i=1, nx-1
               dxm(i,j) = ( mm(i+1,j) - mm(i,j) ) * dxi
            end do
         end do
!$OMP DO
         do i=1, nx
            dyn(i,ny) = ( 0.0 - nn(i,ny)*cos(thetal(ny)) ) * dyi
            do j=1, ny-1
               dyn(i,j) = ( nn(i,j+1)*cos(thetal(j+1)) - nn(i,j)*cos(thetal(j)) ) * dyi
            end do
         end do
!$OMP END DO
!$OMP END PARALLEL
       end block



       !! ------------------------------------------------------------------ !!
       block

         integer i, j
         real :: dtxi, dtyi
         !! --

         dtxi = dt * dxi
         dtyi = dt * dyi
!$OMP PARALLEL SHARED(eta) PRIVATE(i,j)
!$OMP DO
         do j=na+1, ny-na
            do i=na+1, nx-na
               eta(i,j) =  eta(i,j) - ( dxm(i,j) +  dyn(i,j) )*dt* (1/(cos(theta(j)) * r)) * fh(i,j)
            end do
         end do
!$OMP END DO
!$OMP END PARALLEL
       end block



       !! ------------------------------------------------------------------ !!
       block

         integer :: i, j
         real :: eta_y
         !! --
!$OMP PARALLEL SHARED(eta,eta_x) PRIVATE(i,j,eta_y)
!$OMP DO
         do j=1, na
            do i=1, nx
               eta_y = eta(i,j) - eta_x(i,j)

               eta_x(i,j) = g1x(2,i) * eta_x(i,j) - g2x(2,i) * dxm(i,j) * (1/(cos(theta(j)) * r))* dt
               eta_y      = g1y(2,j) * eta_y      - g2y(2,j) * dyn(i,j) * (1/(cos(thetal(j)) * r))* dt

               eta(i,j) = ( eta_x(i,j) + eta_y ) * fh(i,j)

            end do
         end do
!$OMP DO
         do j=ny-na+1,ny
            do i=1, nx
               eta_y = eta(i,j) - eta_x(i,j)

               eta_x(i,j) = g1x(2,i) * eta_x(i,j) - g2x(2,i) * dxm(i,j) * (1/(cos(theta(j)) * r))* dt
               eta_y      = g1y(2,j) * eta_y      - g2y(2,j) * dyn(i,j) * (1/(cos(thetal(j)) * r))* dt

               eta(i,j) = ( eta_x(i,j) + eta_y ) * fh(i,j)

            end do
         end do
!$OMP DO
         do j=na+1, ny-na
            do i=1, na
               eta_y = eta(i,j) - eta_x(i,j)

               eta_x(i,j) = g1x(2,i) * eta_x(i,j) - g2x(2,i) * dxm(i,j) * (1/(cos(theta(j)) * r))* dt
               eta_y      = g1y(2,j) * eta_y      - g2y(2,j) * dyn(i,j) * (1/(cos(thetal(j)) * r))* dt

               eta(i,j) = ( eta_x(i,j) + eta_y ) * fh(i,j)

            end do
         end do
!$OMP DO
         do j=na+1, ny-na
            do i=nx-na+1,nx
               eta_y = eta(i,j) - eta_x(i,j)

               eta_x(i,j) = g1x(2,i) * eta_x(i,j) - g2x(2,i) * dxm(i,j) * (1/(cos(theta(j)) * r))* dt
               eta_y      = g1y(2,j) * eta_y      - g2y(2,j) * dyn(i,j) * (1/(cos(thetal(j)) * r))* dt

               eta(i,j) = ( eta_x(i,j) + eta_y ) * fh(i,j)

            end do
         end do
!$OMP END DO
!$OMP END PARALLEL
       end block
       
       !! ------------------------------------------------------------------ !!
       block

         integer :: i
         
         if (mod(int(dt*it), tint)==0.0) then
!$OMP PARALLEL SHARED(wav) PRIVATE(i)
!$OMP DO
            do i=1, nst
               wav(i,it) = eta( ist(i), jst(i) )
            end do
!$OMP END DO
!$OMP END PARALLEL
         endif

       end block
       
       !! ----------------------------
       ! Modified by Fauzi
       ! To create output file of eta
       block
      
          integer :: i,j
          character(len=20):: fileout, fld
        
          if (mov.eq.0)then      
          if(mod(int(it*dt), 60 )==0) then
             write(fileout,'(a, i5.5, a)')"out/eta_", int(it*dt), ".dat"
             open(21,file=trim(fileout), status='unknown')     

             do i=1,nx
                write(21,'(10000f9.2)')(eta(i,j), j=1,ny)
             enddo

             close(21)
         endif
         endif
        
       end block
      !! ----------------------------

    end do

  end block

  !! ----------------------------------------------------------------------- !!
  block

    integer :: i, it
    character(256) :: fn
    integer, parameter :: io = 100
    !! --

    do i=1, nst
       write(fn,'(A,I3.3,A)') 'tsunami_', i, '.dat'
       open(io, file=fn, action='write', status='unknown' )
       do it=0, nt
          if (mod(int(it*dt), tint) == 0.0) then
             write(io,*) it*dt, wav(i,it)
          endif
       end do
       close(io)
    end do

  end block

  !! ----------------------------------------------------------------------- !!
  block

    deallocate( eta )
    deallocate( mm  )
    deallocate( nn  )
    deallocate( hh  )
    deallocate( hm  )
    deallocate( hn  )
    deallocate( fh  )
    deallocate( fm  )
    deallocate( fn  )
    deallocate( ist, jst )
    deallocate( wav )
    deallocate( g1x, g2x, g1y, g2y )
    deallocate( eta_x )
    deallocate( theta, thetal )

  end block

end program llw
!! ------------------------------------------------------------------------- !!
