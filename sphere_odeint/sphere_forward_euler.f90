subroutine sphere_forward_euler(t0,t1,tn,dt,lat,lon, &
     utx1,unx1,uty1,uny1,uc1, &
     utx2,unx2,uty2,uny2,uc2, &
     vtx1,vnx1,vty1,vny1,vc1, &
     vtx2,vnx2,vty2,vny2,vc2,kx,ky)
  implicit none

  integer t0, t1, tn, dt
  real*8 :: lat,lon
  integer unx1, uny1
  real*8 utx1(unx1), uty1(uny1)
  real*8 :: uc1((unx1-kx-1)*(uny1-ky-1))
  integer unx2, uny2
  real*8 utx2(unx2), uty2(uny2)
  real*8 :: uc2((unx2-kx-1)*(uny2-ky-1))
  integer vnx1, vny1
  real*8 vtx1(vnx1), vty1(vny1)
  real*8 :: vc1((vnx1-kx-1)*(vny1-ky-1))
  integer vnx2, vny2
  real*8 vtx2(vnx2), vty2(vny2)
  real*8 :: vc2((vnx2-kx-1)*(vny2-ky-1))
  integer,parameter :: lwrk=8, kwrk=2
  ! integer,parameter :: lwrk=372, kwrk=22
  real*8 :: wrk(lwrk)
  integer :: iwrk(kwrk)
  real*8 u1(1), u2(1), v1(1), v2(1)
  integer :: ier
  integer :: kx,ky
  ! real*8 traj(1+(tn-t1)/dt,2)

  real*8, parameter :: SEC_PER_HR=3600
  real*8, parameter :: PI = 2*acos(0.)
  real*8, parameter :: TWO_PI = 2*PI
  real*8, parameter :: HALF_PI = .5*PI
  real*8, parameter :: EARTH_RADIUS=6371.008
  real*8, parameter :: MIN_LAT = -78.375/180*PI
  real*8, parameter :: MAX_LAT = 78.375/180*PI
  integer :: t
  real*8 t_dur, dt_hr
  real*8 sin_lat,cos_lat
  real*8 sin_lon,cos_lon
  real*8 p(3),dp(3),lat_tan(3),lon_tan(3),w(2)
  ! integer :: usher

  dt_hr = dt/SEC_PER_HR
  t_dur=tn-t0

  t = t1-dt
  ! usher = 1

  do while (t < tn)

     ! precompute trig vals
     sin_lat = sin(lat)
     cos_lat = cos(lat)
     sin_lon = sin(lon)
     cos_lon = cos(lon)

     ! project up to surface of earth
     p(1) = cos_lat*cos_lon
     p(2) = cos_lat*sin_lon
     p(3) = sin_lat
     p = EARTH_RADIUS*p

     ! vectors tangent to lat and lon
     lat_tan(1) = -sin_lon
     lat_tan(2) = cos_lon
     lat_tan(3) = 0
     lon_tan(1) = -sin_lat*cos_lon
     lon_tan(2) = -sin_lat*sin_lon
     lon_tan(3) = cos_lat

     ! flow forward
     lat = lat + HALF_PI
     lon = lon + PI

     call bispev(utx1,unx1,uty1,uny1,uc1,3,3, &
          (/lat/),1,(/lon/),1,u1,wrk,lwrk,iwrk,kwrk,ier)
     call bispev(utx2,unx2,uty2,uny2,uc2,3,3, &
          (/lat/),1,(/lon/),1,u2,wrk,lwrk,iwrk,kwrk,ier)
     call bispev(vtx1,vnx1,vty1,vny1,vc1,3,3, &
          (/lat/),1,(/lon/),1,v1,wrk,lwrk,iwrk,kwrk,ier)
     call bispev(vtx2,vnx2,vty2,vny2,vc2,3,3, &
          (/lat/),1,(/lon/),1,v2,wrk,lwrk,iwrk,kwrk,ier)

     w(1) = tn-t
     w(2) = t-t0
     w = w / t_dur
     dp = (w(1)*u1(1) + w(2)*u2(1))*lat_tan &
          + (w(1)*v1(1) + w(2)*v2(1))*lon_tan
     p = p + dt_hr*dp

     ! nondimensionalize
     p = p / norm2(p)
     lat = asin(p(3))
     lon = atan2(p(2),p(1))

     ! truncate lats out of bounds
     lat = max(lat,MIN_LAT)
     lat = min(lat,MAX_LAT)

     ! traj(usher,1) = lat
     ! traj(usher,2) = lon

     t = t + dt
     ! usher = usher + 1

  end do

  ! print*, idx
  ! print*, traj(1,2)

end subroutine sphere_forward_euler
