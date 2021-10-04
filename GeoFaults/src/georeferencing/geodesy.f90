module geodesy

    implicit none

    ! declaration of Earth ellipsoid parameters
    ! from Longley et al., 2001, table 4.2, WGS84 ellipsoid case
    real (kind=r8b), parameter :: Earth_major_radius = 6378.137 !km
    ! from http://home.online.no/~sigurdhu/WGS84_Eng.html, and also eq. 6.29, p. 225 in Leick, 1995
    real (kind=r8b), parameter :: eccentricity_sq = 0.00669437999014



contains


end module geodesy