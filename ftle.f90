function dlat(lat1, lat2)
  implicit none
  real*8 lat1,lat2,dlat
  dlat = acos(cos(lat1)*cos(lat2)+sin(lat1)*sin(lat2))
end function dlat

function dlon(lon1, lon2)
  implicit none
  real*8 lon1,lon2,dlon
  dlon = acos(cos(lon1-lon2))
end function dlon

subroutine ftle(n_lat,n_lon,x,y,z,T,e)
  implicit none

  integer n_lat,n_lon
  real*8 x(n_lat,n_lon),y(n_lat,n_lon),z(n_lat,n_lon,2)
  real*8 T
  real*8 e(n_lat-2,n_lon)

  real*8, external :: dlat,dlon

  integer lat_idx,lon_idx
  integer lat_idx_pos,lon_idx_pos
  integer lat_idx_neg,lon_idx_neg
  real*8 J(2,2),a,b,max_ev

  do lat_idx=2,n_lat-1
     do lon_idx=1,n_lon

        lat_idx_pos = lat_idx + 1
        lat_idx_neg = lat_idx - 1
        lon_idx_pos = lon_idx + 1
        if (lon_idx_pos > n_lon) then
           lon_idx_pos = 1
        end if
        lon_idx_neg = lon_idx - 1
        if (lon_idx_neg < 1) then
           lon_idx_neg = n_lon
        end if

        J(1,1) = dlat(z(lat_idx_pos,lon_idx,1), z(lat_idx_neg,lon_idx,1))
        J(1,2) = dlat(z(lat_idx,lon_idx_pos,1), z(lat_idx,lon_idx_neg,1))

        J(2,1) = dlon(z(lat_idx_pos,lon_idx,2), z(lat_idx_neg,lon_idx,2))
        J(2,2) = dlon(z(lat_idx,lon_idx_pos,2), z(lat_idx,lon_idx_neg,2))

        J(:,1) = J(:,1) / dlat(x(lat_idx_pos,lon_idx), x(lat_idx_neg,lon_idx))
        J(:,2) = J(:,2) / dlon(y(lat_idx,lon_idx_pos), y(lat_idx,lon_idx_neg))

        J = matmul(transpose(J),J)

        a = J(1,1)+J(2,2)
        b = J(1,1)*J(2,2)-J(1,2)*J(2,1)
        a = .5*a
        b = sqrt(a**2-b)
        max_ev = max(a+b,a-b)

        print*,a+b,a-b

        e(lat_idx-1,lon_idx) = log(sqrt(max_ev))/T

     end do
  end do

end subroutine ftle
