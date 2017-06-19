subroutine spherical_forward_euler(t0,tn,dt,x0,x)
  implicit none

  real(8) dt
  real(8), dimension(2) :: w
  real(8), dimension(n_lat,n_lon,2) :: x0, u, v
  real(8), dimension(n_lat,n_lon,2) :: x

  !f2py intent(in) :: n_lat,n_lon,dt,w,x0,u,v
  !f2py intent(out) :: x

  real(8), parameter :: PI = 2*asin(1.)
  real(8), parameter :: TWO_PI = 2*PI
  real(8), parameter :: EARTH_RADIUS = 6371.008
  real(8), parameter :: MIN_LAT_RAD = -78.375/180*PI
  real(8), parameter :: MAX_LAT_RAD = 78.375/180*PI

  integer lat_idx,lon_idx
  real(8) lat,lon,sin_lat,cos_lat,sin_lon,cos_lon
  real(8), dimension(3) :: p,dp,lat_tan,lon_tan

  print*, n_lat, n_lon
  do lat_idx=1,n_lat
     do lon_idx=1,n_lon

        ! precompute trig vals
        lat = x0(lat_idx,lon_idx,1)
        lon = x0(lat_idx,lon_idx,2)
        sin_lat = sin(lat)
        cos_lat = cos(lat)
        sin_lon = sin(lon)
        cos_lon = cos(lon)

        ! project up to surface of earth
        p = EARTH_RADIUS*(/cos_lat*cos_lon,cos_lat*sin_lon,sin_lat/)

        ! vectors tangent to lat and lon
        lat_tan = (/-sin_lon,cos_lon,real(0.,8)/)
        lon_tan = (/-sin_lat*cos_lon,-sin_lat*sin_lon,cos_lat/)

        ! flow forward
        dp = (w(1)*u(lat_idx,lon_idx,1) + w(2)*u(lat_idx,lon_idx,2))*lat_tan
        dp = dp + (w(1)*v(lat_idx,lon_idx,1) + w(2)*v(lat_idx,lon_idx,2))*lon_tan
        dp = dp*dt
        p = p + dp

        ! nondimensionalize
        p = p/norm2(p)
        lon = asin(p(3))
        lat = atan2(p(2),p(1))

        ! truncate lats out of bounds and wrap lons around
        lat = max(lat,MIN_LAT_RAD)
        lat = min(lat,MAX_LAT_RAD)
        if (lon>PI) then
           lon = lon-TWO_PI
        else if (lon<-PI) then
           lon = lon+TWO_PI
        end if

        x(lat_idx,lon_idx,:)=(/lat,lon/)
     end do
  end do

end subroutine spherical_forward_euler
