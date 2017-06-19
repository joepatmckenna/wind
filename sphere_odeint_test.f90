program sphere_odeint_test
  implicit none

  integer,parameter :: mx=628, my=1440
  real*8 :: x(mx), y(my), z(mx*my), z0, z1
  real*8 :: xx(1), yy(1)
  real*8 :: s=0.
  integer,parameter :: nxest=mx+6, nyest=my+7
  integer unx1,unx2,vnx1,vnx2
  integer uny1,uny2,vny1,vny2
  real*8 :: utx1(nxest),uty1(nyest),utx2(nxest),uty2(nyest)
  real*8 :: vtx1(nxest),vty1(nyest),vtx2(nxest),vty2(nyest)
  real*8 :: uc1((nxest-4)*(nyest-4)),uc2((nxest-4)*(nyest-4))
  real*8 :: vc1((nxest-4)*(nyest-4)),vc2((nxest-4)*(nyest-4))
  real*8 fp
  integer,parameter :: lwrk=12+nxest*(my+nyest+3)+nyest*24+4*mx+8*my+max(my+nyest,nxest)
  integer,parameter :: kwrk=5+mx+my+nxest+nyest
  real*8 :: wrk(lwrk)
  integer :: iwrk(kwrk)
  integer :: ier

  integer i
  integer, parameter :: N_LAT = 628, N_LON = 1440
  real, parameter :: PI = 2*acos(0.)
  real, parameter :: MIN_LAT = (-78.375+90)*PI/180, MAX_LAT = (78.375+90)*PI/180
  real, parameter :: MIN_LON = (-179.875+180)*PI/180, MAX_LON = (179.85+180)*PI/180

  integer, parameter :: t0=0, t1=0, tn=10, dt=1
  real*8 :: lat0=0,lon0=0
  real*8 :: traj((tn-t1)/dt,2)

  x = (/(MIN_LAT+i*(MAX_LAT-MIN_LAT)/(mx-1),i=0,mx-1)/)
  y = (/(MIN_LON+i*(MAX_LON-MIN_LON)/(my-1),i=0,my-1)/)
  open (10, file='./u.txt', action='read')
  do i=1,N_LAT
     read(10,*) z(1+(i-1)*N_LON:i*N_LON)
  end do
  close(10)

  call spgrid((/0,0,0/),(/-1,0,-1,0/),mx,x,my,y,z,z0,z1,s,nxest,nyest, &
       unx1,utx1,uny1,uty1,uc1,fp,wrk,lwrk,iwrk,kwrk,ier)
  print*, ier
  call spgrid((/0,0,0/),(/-1,0,-1,0/),mx,x,my,y,z,z0,z1,s,nxest,nyest, &
       unx2,utx2,uny2,uty2,uc2,fp,wrk,lwrk,iwrk,kwrk,ier)
  print*, ier
  call spgrid((/0,0,0/),(/-1,0,-1,0/),mx,x,my,y,z,z0,z1,s,nxest,nyest, &
       vnx1,vtx1,vny1,vty1,vc1,fp,wrk,lwrk,iwrk,kwrk,ier)
  print*, ier
  call spgrid((/0,0,0/),(/-1,0,-1,0/),mx,x,my,y,z,z0,z1,s,nxest,nyest, &
       vnx2,vtx2,vny2,vty2,vc2,fp,wrk,lwrk,iwrk,kwrk,ier)
  print*, ier

  call sphere_forward_euler(t0,t1,tn,dt,lat0,lon0, &
       utx1,unx1,uty1,uny1,uc1, &
       utx2,unx2,uty2,uny2,uc2, &
       vtx1,vnx1,vty1,vny1,vc1, &
       vtx2,vnx2,vty2,vny2,vc2,3,3,traj)

end program sphere_odeint_test
