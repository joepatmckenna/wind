program main
  use wind
  implicit none

  integer i
  integer, parameter :: N_LAT = 628, N_LON = 1440
  real, parameter :: PI = 2*acos(0.)
  real, parameter :: MIN_LAT = (-78.375+90)*PI/180, MAX_LAT = (78.375+90)*PI/180
  real, parameter :: MIN_LON = (-179.875+180)*PI/180, MAX_LON = (179.85+180)*PI/180

  integer,parameter :: mx=628, my=1440
  real :: x(mx), y(my), r(mx*my), z(mx*my)
  integer nx, ny
  real, allocatable :: tx(:), ty(:), c(:)
  integer ier

  x(:) = (/(MIN_LAT+i*(MAX_LAT-MIN_LAT)/(mx-1),i=0,mx-1)/)
  y(:) = (/(MIN_LON+i*(MAX_LON-MIN_LON)/(my-1),i=0,my-1)/)

  open (10, file='./u.txt', action='read')
  do i=1,N_LAT
     read(10,*) r(1+(i-1)*N_LON:i*N_LON)
  end do
  close(10)

  call rep(mx,x,my,y,r,nx,tx,ny,ty,c,ier)

  call ev(tx,nx,ty,ny,c,x,mx,y,my,z,ier)

  ! print*, 'err', sum(abs(z-r))

  ! from example
  ! call spgrid(iopt,ider,mu,u,mv,v,r,r0,r1,s,nuest,nvest, &
  !      nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
  ! call bispev(tu,nu,tv,nv,c,3,3,u,mu,v,mv,f,wk,100,iw,25,ier)
end program main

! reset
! rm *.o *.mod;
! gfortran -c header.f90
! gfortran spline_test.f90 header.o dierckx/*.f -o spline_test
! ./spline_test


! f2py -c -m wind header.f90 dierckx/*.f

