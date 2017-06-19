program main
  implicit none

  integer, parameter :: N_LAT = 628, N_LON = 1440
  real, parameter :: PI = 2*acos(0.)
  real, parameter :: MIN_LAT = -78.375, MAX_LAT = 78.375
  real, parameter :: MIN_LON = -179.875, MAX_LON = 179.85

  ! construct spline parameters
  integer, parameter :: iopt(3) = (/0,0,0/), ider(4)=(/-1,0,-1,0/)
  integer, parameter :: mu=N_LAT, mv=N_LON
  integer, parameter :: nuest=mu+6, nvest=mv+7
  integer, parameter :: lwrk=12+nuest*(mv+nvest+3)+nvest*24+4*mu+8*mv+mv+nvest
  integer, parameter :: kwrk=5+mu+mv+nuest+nvest

  real :: u(mu), v(mv), r(mu*mv), r0, r1, s=0
  integer iwrk(kwrk)
  real wrk(lwrk)
  integer nu, nv
  real tu(nuest), tv(nvest), c((nuest-4)*(nvest-4)), fp
  integer ier

  ! evaluate spline parameters
  integer :: mx, my
  real, allocatable :: x(:),y(:),z(:)
  real, allocatable :: wk(:)
  integer, allocatable :: iwk(:)
  integer lwk, kwk

  ! dummy
  integer i

  u=PI*((/(MIN_LAT+i*(MAX_LAT-MIN_LAT)/(N_LAT-1),i=0,N_LAT-1)/)+90)/180
  v=PI*((/(MIN_LON+i*(MAX_LON-MIN_LON)/(N_LON-1),i=0,N_LON-1)/)+180)/180

  open (10, file='./u.txt', action='read')
  do i=1,N_LAT
     read(10,*) r(1+(i-1)*N_LON:i*N_LON)
  end do

  ! construct spline
  call spgrid(iopt,ider,mu,u,mv,v,r,r0,r1,s,nuest,nvest,nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
  print*, ier
  print*, 'pi', PI
  print*, 'mu test', mu >= 4-min(1,ider(1)+1)-min(1,ider(3)+1)-ider(2)-ider(4)
  print*, 'u between 0 and PI', (u(1)>0).and.(u(mu)<PI), minval(u), maxval(u)
  print*, 'mv test', mv > 3
  print*, 'v test', (v(1)>=-PI.and.v(1)<PI), (v(mv)<v(1)+2*PI), minval(v), maxval(v)
  print*, 'large enough nuest/nvest',mu+6+iopt(2)+iopt(3), mv+7
  print*, 'num u/v knots', nu, nv
  print*, 'fp', fp
  print*, 'lwrk test', lwrk, lwrk >= 12+nuest*(mv+nvest+3)+nvest*24+4*mu+8*mv+max((mv+nvest),nuest)
  print*, 'kwkrk test', kwrk, kwrk >= 5+mu+mv+nuest+nvest
  print*, -1<=iopt(1), iopt(1)<=1, 0<=iopt(2), iopt(2)<=1, 0<=iopt(3), iopt(3)<=1
  print*, -1<=ider(1), ider(1)<=1, 0<=ider(2), ider(2)<=1, ider(2)==0
  print*, -1<=ider(3), ider(3)<=1, 0<=ider(4), ider(4)<=1, ider(4)==0
  print*, mv >= 4, nuest >=8, nvest >= 8
  print*, kwrk>=5+mu+mv+nuest+nvest
  print*, lwrk >= 12+nuest*(mv+nvest+3)+nvest*24+4*mu+8*mv+max(nuest,mv+nvest)
  print*, s>=0
  print*, nuest>=mu+6+iopt(2)+iopt(3), nuest, mu+6+iopt(2)+iopt(3)
  print*, nvest>=mv+7, nvest, mv+7


  mx=N_LAT; my=N_LON
  lwk=mx*4+my*4; kwk=mx+my

  allocate(x(mx),y(my),z(mx*my),wk(lwk),iwk(kwk))
  x=PI*((/(MIN_LAT+i*(MAX_LAT-MIN_LAT)/(mx-1),i=0,mx-1)/)+90)/180
  y=PI*((/(MIN_LON+i*(MAX_LON-MIN_LON)/(my-1),i=0,my-1)/)+180)/180

  do i=2,10!mx,my
     print*, tu(4) <= x(i-1), x(i-1) <= x(i), x(i) <= tu(nu-3), x(i), tu(nu-3)
     print*, tv(4) <= y(i-1), y(i-1) <= y(i), y(i) <= tv(nv-3), y(i), tv(nv-3)
  end do

  ! evaluate spline
  call bispev(tu,nu,tv,nv,c,3,3,x,mx,y,my,z,wk,lwk,iwk,kwk,ier)
  print*, ier

  print*, 'err:', sum(abs(z-r))

  deallocate(x,y,z,wk,iwk)

end program main

! gfortran spline_test.bak.f90 dierckx/spgrid.f dierckx/fpspgr.f dierckx/fpchec.f dierckx/fpchep.f dierckx/fpknot.f dierckx/fpopsp.f dierckx/fprati.f dierckx/fpgrsp.f dierckx/fpsysy.f dierckx/fpback.f dierckx/fpbacp.f dierckx/fpbspl.f dierckx/fpcyt1.f dierckx/fpcyt2.f dierckx/fpdisc.f dierckx/fpgivs.f dierckx/fprota.f dierckx/bispev.f dierckx/fpbisp.f

