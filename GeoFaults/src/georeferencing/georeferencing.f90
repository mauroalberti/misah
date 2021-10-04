module georeferencing

    implicit none

    type :: spherical_loc
        sequence
        real(kind=r8b) :: lat,lon, depth ! lat and lon in decimal degrees, depth in km
    end type spherical_loc



    contains


end module georeferencing