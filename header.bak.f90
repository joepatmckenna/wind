module wind
  implicit none

contains

  subroutine rep(mx,x,my,y,r,nx,tx,ny,ty,c,ier)
    implicit none

    integer :: iopt(3)=(/0,0,0/), ider(4)=(/-1,0,-1,0/)
    integer mx, my
    real :: x(mx), y(my), r(mx*my)
    real r0,r1
    real :: s=0
    integer :: nxest, nyest
    integer nx, ny
    real, allocatable :: tx(:), ty(:), c(:)
    real fp
    integer lwrk, kwrk
    real, allocatable :: wrk(:)
    integer, allocatable :: iwrk(:)
    integer ier

    nxest=mx+6
    nyest=my+7
    lwrk=12+nxest*(my+nyest+3)+nyest*24+4*mx+8*my+max(my+nyest,nxest)
    kwrk=5+mx+my+nxest+nyest
    allocate(tx(nxest), ty(nyest), c((nxest-4)*(nyest-4)), wrk(lwrk), iwrk(kwrk))

    ! print*, 'mx test', mx >= 4-min(1,ider(1)+1)-min(1,ider(3)+1)-ider(2)-ider(4)
    ! print*, 'my test', my >= 4
    ! print*, 'nx/ny test', nxest >=8, nyest >= 8
    ! print*, 'kwrk test', kwrk>=5+mx+my+nxest+nyest
    ! print*, 'lwrk test', lwrk>=12+nxest*(my+nxest+3)+nyest*24+4*mx+8*my+max(nxest,my+nyest)
    ! if (s==0) then
    !    print*, 'nxest/nyest test', nxest>=mx+6+iopt(2)+iopt(3), nyest>=my+7
    ! end if

    call spgrid(iopt,ider,mx,x,my,y,r,r0,r1,s,nxest,nyest, &
         nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)

    deallocate(wrk, iwrk)

    ! print*, 'rep'
    ! print*, 'ier', ier
    ! print*, 'fp', fp
    ! print*, 'nx/ny', nx, ny

  end subroutine rep

  subroutine ev(tx,nx,ty,ny,c,x,mx,y,my,z,ier)
    implicit none

    integer nx, ny
    real tx(nx), ty(ny)
    real c((nx-4)*(ny-4))
    integer mx, my
    real x(mx), y(my), z(mx*my)
    integer lwrk,kwrk
    real, allocatable :: wrk(:)
    integer, allocatable :: iwrk(:)
    integer ier

    lwrk=4*mx+4*my
    kwrk=mx+my

    allocate(wrk(lwrk), iwrk(kwrk))

    call bispev(tx,nx,ty,ny,c,3,3,x,mx,y,my,z,wrk,lwrk, &
         iwrk,kwrk,ier)

    ! print*, 'rep'
    ! print*, 'ier', ier
    ! print*, 'fp', fp
    ! print*, 'nx/ny', nx, ny

    deallocate(wrk,iwrk)

  end subroutine ev

end module wind
